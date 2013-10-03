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

/*******************************************************************************
   We now include this header file if it was not included before.
*******************************************************************************/
#ifndef sse2_macro_gcc_h_included
#define sse2_macro_gcc_h_included

/*
SSE2-ified version of PAIR_SQUARE_4. Data enter in [__tAr, ~tAr], [__tAi, ~tAi], pairs, where the imaginary part
of each input pair is assumed offset +0x10 in memory from the real part, i.e. needs no explicit pointer reference.

For the sincos twiddles: using the notation of the scalar PAIR_SQUARE_4() macro,"__c" means [c0,s1], "__s" means [s0,c1].
For these, due to the buterfly indexing pattern, we cannot assume that __s = __c + 0x10, so feed both pointers explicitly.

We use shufpd xmm, xmm, 1 to swap lo and hi doubles of an xmm register for the various operations with one swapped input.
*/
	/* Complex multiply of 2 roots of unity - use e.g. for "multiply up" of sincos twiddles. */
#ifdef USE_AVX

	#define SSE2_CMUL_EXPO(XcA,XcB,XcAmB,XcApB)\
	{\
	__asm__ volatile (\
		"movq	%[__cA]		,%%rax\n\t"\
		"movq	%[__cB]		,%%rbx\n\t"\
		"movq	%[__cAmB]	,%%rcx\n\t"\
		"movq	%[__cApB]	,%%rdx\n\t"\
		"\n\t"\
		"vmovaps	    (%%rax),%%ymm0\n\t"\
		"vmovaps	0x20(%%rax),%%ymm2\n\t"\
		"vmovaps	    (%%rbx),%%ymm4\n\t"\
		"vmovaps	0x20(%%rbx),%%ymm5\n\t"\
		"vmovaps	%%ymm0,%%ymm1\n\t"\
		"vmovaps	%%ymm2,%%ymm3\n\t"\
		"\n\t"\
		"vmulpd	%%ymm4,%%ymm0,%%ymm0\n\t"\
		"vmulpd	%%ymm5,%%ymm1,%%ymm1\n\t"\
		"vmulpd	%%ymm4,%%ymm2,%%ymm2\n\t"\
		"vmulpd	%%ymm5,%%ymm3,%%ymm3\n\t"\
		"vmovaps	%%ymm0,%%ymm4\n\t"\
		"vmovaps	%%ymm1,%%ymm5\n\t"\
		"vaddpd	%%ymm3,%%ymm0,%%ymm0\n\t"\
		"vsubpd	%%ymm2,%%ymm1,%%ymm1\n\t"\
		"vsubpd	%%ymm3,%%ymm4,%%ymm4\n\t"\
		"vaddpd	%%ymm2,%%ymm5,%%ymm5\n\t"\
		"vmovaps	%%ymm0,    (%%rcx)\n\t"\
		"vmovaps	%%ymm1,0x20(%%rcx)\n\t"\
		"vmovaps	%%ymm4,    (%%rdx)\n\t"\
		"vmovaps	%%ymm5,0x20(%%rdx)\n\t"\
		:					/* outputs: none */\
		: [__cA]  "m" (XcA)	/* All inputs from memory addresses here */\
		 ,[__cB]  "m" (XcB)\
		 ,[__cAmB] "m" (XcAmB)\
		 ,[__cApB] "m" (XcApB)\
		: "cc","memory","rax","rbx","rcx","rdx","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5"		/* Clobbered registers */\
	);\
	}

	// AVX version has shufpd immediate = 5 = 0101_2, which is the doubled analog of the SSE2 imm8 = 1 = 01_2:
	#define PAIR_SQUARE_4_SSE2(XtAr, XtBr, XtCr, XtDr, Xc, Xs, Xforth)\
	{\
	__asm__ volatile (\
		"movq	%[__tDr]	,%%rdx		\n\t"\
		"movq	%[__tAr]	,%%rax		\n\t"\
		"\n\t"\
		"vmovaps	    (%%rdx)	,%%ymm6		\n\t"\
		"vmovaps	0x20(%%rdx)	,%%ymm7		\n\t"\
		"vmovaps	    (%%rax)	,%%ymm0		\n\t"\
		"vmovaps	0x20(%%rax)	,%%ymm3		\n\t"\
		"vshufpd	$5,%%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vshufpd	$5,%%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vmovaps	    (%%rax)	,%%ymm2		\n\t"\
		"vmovaps	0x20(%%rax)	,%%ymm1		\n\t"\
		"\n\t"\
		"vmulpd	%%ymm6,%%ymm0,%%ymm0		\n\t"\
		"vmulpd	%%ymm7,%%ymm3,%%ymm3		\n\t"\
		"vmulpd	%%ymm6,%%ymm1,%%ymm1		\n\t"\
		"vmulpd	%%ymm7,%%ymm2,%%ymm2		\n\t"\
		"vaddpd	%%ymm3,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	%%ymm2,%%ymm1,%%ymm1		\n\t"\
		"vaddpd	%%ymm0,%%ymm0,%%ymm0		\n\t"\
		"vaddpd	%%ymm1,%%ymm1,%%ymm1		\n\t"\
		"\n\t"\
		"movq	%[__tCr]	,%%rcx		\n\t"\
		"movq	%[__tBr]	,%%rbx		\n\t"\
		"\n\t"\
		"vmovaps	    (%%rcx)	,%%ymm6		\n\t"\
		"vmovaps	0x20(%%rcx)	,%%ymm7		\n\t"\
		"vmovaps	    (%%rbx)	,%%ymm2		\n\t"\
		"vmovaps	0x20(%%rbx)	,%%ymm5		\n\t"\
		"vshufpd	$5,%%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vshufpd	$5,%%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vmovaps	    (%%rbx)	,%%ymm4		\n\t"\
		"vmovaps	0x20(%%rbx)	,%%ymm3		\n\t"\
		"\n\t"\
		"vmulpd	%%ymm6,%%ymm2,%%ymm2		\n\t"\
		"vmulpd	%%ymm7,%%ymm5,%%ymm5		\n\t"\
		"vmulpd	%%ymm6,%%ymm3,%%ymm3		\n\t"\
		"vmulpd	%%ymm7,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	%%ymm5,%%ymm2,%%ymm2		\n\t"\
		"vsubpd	%%ymm4,%%ymm3,%%ymm3		\n\t"\
		"vaddpd	%%ymm2,%%ymm2,%%ymm2		\n\t"\
		"vaddpd	%%ymm3,%%ymm3,%%ymm3		\n\t"\
		"\n\t"\
		"vmovaps	    (%%rax)	,%%ymm4		\n\t"\
		"vmovaps	0x20(%%rax)	,%%ymm5		\n\t"\
		"vsubpd	%%ymm5,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	%%ymm4,%%ymm5,%%ymm5		\n\t"\
		"vmulpd	%%ymm5,%%ymm4,%%ymm4		\n\t"\
		"\n\t"\
		"vmovaps	    (%%rax)	,%%ymm5		\n\t"\
		"vmulpd	0x20(%%rax)	,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	%%ymm5		,%%ymm5,%%ymm5		\n\t"\
		"vmovaps	%%ymm4	,    (%%rax)	\n\t"\
		"vmovaps	%%ymm5	,0x20(%%rax)	\n\t"\
		"\n\t"\
		"vsubpd	%%ymm4,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	%%ymm5,%%ymm1,%%ymm1		\n\t"\
		"\n\t"\
		"vmovaps	    (%%rbx)	,%%ymm6		\n\t"\
		"vmovaps	0x20(%%rbx)	,%%ymm7		\n\t"\
		"vsubpd	%%ymm7,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	%%ymm6,%%ymm7,%%ymm7		\n\t"\
		"vmulpd	%%ymm7,%%ymm6,%%ymm6		\n\t"\
		"\n\t"\
		"vmovaps	    (%%rbx)	,%%ymm7		\n\t"\
		"vmulpd	0x20(%%rbx)	,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	%%ymm7		,%%ymm7,%%ymm7		\n\t"\
		"vmovaps	%%ymm6	,    (%%rbx)	\n\t"\
		"vmovaps	%%ymm7	,0x20(%%rbx)	\n\t"\
		"\n\t"\
		"vsubpd	%%ymm6,%%ymm2,%%ymm2		\n\t"\
		"vsubpd	%%ymm7,%%ymm3,%%ymm3		\n\t"\
		"\n\t"\
		"vmovaps	    (%%rdx)	,%%ymm4		\n\t"\
		"vmovaps	0x20(%%rdx)	,%%ymm5		\n\t"\
		"vsubpd	%%ymm5,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	%%ymm4,%%ymm5,%%ymm5		\n\t"\
		"vmulpd	%%ymm5,%%ymm4,%%ymm4		\n\t"\
		"\n\t"\
		"vmovaps	    (%%rdx)	,%%ymm5		\n\t"\
		"vmulpd	0x20(%%rdx)	,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	%%ymm5		,%%ymm5,%%ymm5		\n\t"\
		"vmovaps	%%ymm4	,    (%%rdx)	\n\t"\
		"vmovaps	%%ymm5	,0x20(%%rdx)	\n\t"\
		"vshufpd	$5,%%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vshufpd	$5,%%ymm5,%%ymm5,%%ymm5		\n\t"\
		"\n\t"\
		"vsubpd	%%ymm4,%%ymm0,%%ymm0		\n\t"\
		"vaddpd	%%ymm5,%%ymm1,%%ymm1		\n\t"\
		"\n\t"\
		"vmovaps	    (%%rcx)	,%%ymm6		\n\t"\
		"vmovaps	0x20(%%rcx)	,%%ymm7		\n\t"\
		"vsubpd	%%ymm7,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	%%ymm6,%%ymm7,%%ymm7		\n\t"\
		"vmulpd	%%ymm7,%%ymm6,%%ymm6		\n\t"\
		"\n\t"\
		"vmovaps	    (%%rcx)	,%%ymm7		\n\t"\
		"vmulpd	0x20(%%rcx)	,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	%%ymm7		,%%ymm7,%%ymm7		\n\t"\
		"vmovaps	%%ymm6	,    (%%rcx)	\n\t"\
		"vmovaps	%%ymm7	,0x20(%%rcx)	\n\t"\
		"vshufpd	$5,%%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vshufpd	$5,%%ymm7,%%ymm7,%%ymm7		\n\t"\
		"\n\t"\
		"vsubpd	%%ymm6,%%ymm2,%%ymm2		\n\t"\
		"vaddpd	%%ymm7,%%ymm3,%%ymm3		\n\t"\
		"\n\t"\
		"/*** Can do in || with above segment ***/\n\t"\
		"movq	%[__c]		,%%rax		\n\t"\
		"movq	%[__s]		,%%rbx		\n\t"\
		"movq	%[__forth]	,%%rdx		\n\t"\
		"vmovaps	%%ymm0		,%%ymm4		\n\t"\
		"vmovaps	%%ymm1		,%%ymm5		\n\t"\
		"vmulpd	(%%rax)	,%%ymm0,%%ymm0		\n\t"\
		"vmulpd	(%%rax)	,%%ymm1,%%ymm1		\n\t"\
		"vaddpd	%%ymm4	,%%ymm0,%%ymm0		\n\t"\
		"vaddpd	%%ymm5	,%%ymm1,%%ymm1		\n\t"\
		"vmulpd	(%%rbx)	,%%ymm4,%%ymm4		\n\t"\
		"vmulpd	(%%rbx)	,%%ymm5,%%ymm5		\n\t"\
		"vsubpd	%%ymm5	,%%ymm0,%%ymm0		\n\t"\
		"vaddpd	%%ymm4	,%%ymm1,%%ymm1		\n\t"\
		"vmulpd	(%%rdx)	,%%ymm0,%%ymm0		\n\t"\
		"vmulpd	(%%rdx)	,%%ymm1,%%ymm1		\n\t"\
		"\n\t"\
		"/*** Can do in || wjth above segment ***/\n\t"\
		"vmovaps	%%ymm2	,%%ymm6		\n\t"\
		"vmovaps	%%ymm3	,%%ymm7		\n\t"\
		"vmulpd	(%%rbx)	,%%ymm2,%%ymm2		\n\t"\
		"vmulpd	(%%rbx)	,%%ymm3,%%ymm3		\n\t"\
		"vsubpd	%%ymm6	,%%ymm2,%%ymm2		\n\t"\
		"vsubpd	%%ymm7	,%%ymm3,%%ymm3		\n\t"\
		"vmulpd	(%%rax)	,%%ymm6,%%ymm6		\n\t"\
		"vmulpd	(%%rax)	,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	%%ymm7	,%%ymm2,%%ymm2		\n\t"\
		"vsubpd	%%ymm6	,%%ymm3,%%ymm3		\n\t"\
		"vmulpd	(%%rdx)	,%%ymm2,%%ymm2		\n\t"\
		"vmulpd	(%%rdx)	,%%ymm3,%%ymm3		\n\t"\
		"\n\t"\
		"movq	%[__tAr]	,%%rax		\n\t"\
		"movq	%[__tBr]	,%%rbx		\n\t"\
		"\n\t"\
		"vmovaps	    (%%rax)	,%%ymm4		\n\t"\
		"vmovaps	0x20(%%rax)	,%%ymm5		\n\t"\
		"vmovaps	    (%%rbx)	,%%ymm6		\n\t"\
		"vmovaps	0x20(%%rbx)	,%%ymm7		\n\t"\
		"vaddpd	%%ymm0,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	%%ymm1,%%ymm5,%%ymm5		\n\t"\
		"vsubpd	%%ymm2,%%ymm6,%%ymm6		\n\t"\
		"vsubpd	%%ymm3,%%ymm7,%%ymm7		\n\t"\
		"vmovaps	%%ymm4	,    (%%rax)	\n\t"\
		"vmovaps	%%ymm5	,0x20(%%rax)	\n\t"\
		"vmovaps	%%ymm6	,    (%%rbx)	\n\t"\
		"vmovaps	%%ymm7	,0x20(%%rbx)	\n\t"\
		"\n\t"\
		"movq	%[__tCr]	,%%rcx		\n\t"\
		"movq	%[__tDr]	,%%rdx		\n\t"\
		"\n\t"\
		"vshufpd	$5,%%ymm0,%%ymm0,%%ymm0		\n\t"\
		"vshufpd	$5,%%ymm1,%%ymm1,%%ymm1		\n\t"\
		"vshufpd	$5,%%ymm2,%%ymm2,%%ymm2		\n\t"\
		"vshufpd	$5,%%ymm3,%%ymm3,%%ymm3		\n\t"\
		"\n\t"\
		"vmovaps	    (%%rdx)	,%%ymm4		\n\t"\
		"vmovaps	0x20(%%rdx)	,%%ymm5		\n\t"\
		"vmovaps	    (%%rcx)	,%%ymm6		\n\t"\
		"vmovaps	0x20(%%rcx)	,%%ymm7		\n\t"\
		"vaddpd	%%ymm0,%%ymm4,%%ymm4		\n\t"\
		"vsubpd	%%ymm1,%%ymm5,%%ymm5		\n\t"\
		"vsubpd	%%ymm2,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	%%ymm3,%%ymm7,%%ymm7		\n\t"\
		"vmovaps	%%ymm4	,    (%%rdx)	\n\t"\
		"vmovaps	%%ymm5	,0x20(%%rdx)	\n\t"\
		"vmovaps	%%ymm6	,    (%%rcx)	\n\t"\
		"vmovaps	%%ymm7	,0x20(%%rcx)	\n\t"\
		:					/* outputs: none */\
		: [__tAr] "m" (XtAr)	/* All inputs from memory addresses here */\
		 ,[__tBr] "m" (XtBr)\
		 ,[__tCr] "m" (XtCr)\
		 ,[__tDr] "m" (XtDr)\
		 ,[__c] "m" (Xc)\
		 ,[__s] "m" (Xs)\
		 ,[__forth] "m" (Xforth)\
		: "cc","memory","rax","rbx","rcx","rdx","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
	);\
	}

	#define SSE2_RADIX_03_DFT_X2(Xcc0, Xi0,Xi1,Xi2, Xo0,Xo1,Xo2, Xj0,Xj1,Xj2, Xu0,Xu1,Xu2)\
	{\
	__asm__ volatile (\
		"movq	%[__i0],%%rax				\n\t	movq	%[__j0],%%r10		\n\t"\
		"movq	%[__i1],%%rbx				\n\t	movq	%[__j1],%%r11		\n\t"\
		"movq	%[__i2],%%rcx				\n\t	movq	%[__j2],%%r12		\n\t"\
		"movq	%[__cc0],%%rdx				\n\t"\
		"vmovaps	    (%%rbx),%%ymm2		\n\t	vmovaps	    (%%r11),%%ymm10	\n\t"\
		"vmovaps	0x20(%%rbx),%%ymm3		\n\t	vmovaps	0x20(%%r11),%%ymm11	\n\t"\
		"vmovaps	    (%%rax),%%ymm0		\n\t	vmovaps	    (%%r10),%%ymm8 	\n\t"\
		"vmovaps	0x20(%%rax),%%ymm1		\n\t	vmovaps	0x20(%%r10),%%ymm9 	\n\t"\
		"vmovaps	    (%%rcx),%%ymm6		\n\t	vmovaps	    (%%r12),%%ymm14	\n\t"\
		"vmovaps	0x20(%%rcx),%%ymm7		\n\t	vmovaps	0x20(%%r12),%%ymm15	\n\t"\
		"vmovaps	%%ymm2,%%ymm4			\n\t	vmovaps	%%ymm10,%%ymm12		\n\t"\
		"vmovaps	%%ymm3,%%ymm5			\n\t	vmovaps	%%ymm11,%%ymm13		\n\t"\
		"movq	%[__o0],%%rax				\n\t	movq	%[__u0],%%r10		\n\t"\
		"movq	%[__o1],%%rbx				\n\t	movq	%[__u1],%%r11		\n\t"\
		"movq	%[__o2],%%rcx				\n\t	movq	%[__u2],%%r12		\n\t"\
		"vaddpd	%%ymm6,%%ymm2,%%ymm2		\n\t	vaddpd	%%ymm14,%%ymm10,%%ymm10		\n\t"\
		"vaddpd	%%ymm7,%%ymm3,%%ymm3		\n\t	vaddpd	%%ymm15,%%ymm11,%%ymm11		\n\t"\
		"vsubpd	%%ymm6,%%ymm4,%%ymm4		\n\t	vsubpd	%%ymm14,%%ymm12,%%ymm12		\n\t"\
		"vsubpd	%%ymm7,%%ymm5,%%ymm5		\n\t	vsubpd	%%ymm15,%%ymm13,%%ymm13		\n\t"\
		"vaddpd	%%ymm2,%%ymm0,%%ymm0		\n\t	vaddpd	%%ymm10,%%ymm8 ,%%ymm8 		\n\t"\
		"vaddpd	%%ymm3,%%ymm1,%%ymm1		\n\t	vaddpd	%%ymm11,%%ymm9 ,%%ymm9 		\n\t"\
		"vmovaps	    (%%rdx),%%ymm6		\n\t"\
		"vmovaps	0x20(%%rdx),%%ymm7		\n\t"\
		"vmovaps	%%ymm0,    (%%rax)		\n\t	vmovaps	%%ymm8 ,    (%%r10)	\n\t"\
		"vmovaps	%%ymm1,0x20(%%rax)		\n\t	vmovaps	%%ymm9 ,0x20(%%r10)	\n\t"\
		"vmulpd	%%ymm6,%%ymm2,%%ymm2		\n\t	vmulpd	%%ymm6 ,%%ymm10,%%ymm10		\n\t"\
		"vmulpd	%%ymm6,%%ymm3,%%ymm3		\n\t	vmulpd	%%ymm6 ,%%ymm11,%%ymm11		\n\t"\
		"vmulpd	%%ymm7,%%ymm4,%%ymm4		\n\t	vmulpd	%%ymm7 ,%%ymm12,%%ymm12		\n\t"\
		"vmulpd	%%ymm7,%%ymm5,%%ymm5		\n\t	vmulpd	%%ymm7 ,%%ymm13,%%ymm13		\n\t"\
		"vaddpd	%%ymm0,%%ymm2,%%ymm2		\n\t	vaddpd	%%ymm8 ,%%ymm10,%%ymm10		\n\t"\
		"vaddpd	%%ymm1,%%ymm3,%%ymm3		\n\t	vaddpd	%%ymm9 ,%%ymm11,%%ymm11		\n\t"\
		"vmovaps	%%ymm2,%%ymm0			\n\t	vmovaps	%%ymm10,%%ymm8 		\n\t"\
		"vmovaps	%%ymm3,%%ymm1			\n\t	vmovaps	%%ymm11,%%ymm9 		\n\t"\
		"vsubpd	%%ymm5,%%ymm2,%%ymm2		\n\t	vsubpd	%%ymm13,%%ymm10,%%ymm10		\n\t"\
		"vaddpd	%%ymm4,%%ymm3,%%ymm3		\n\t	vaddpd	%%ymm12,%%ymm11,%%ymm11		\n\t"\
		"vaddpd	%%ymm5,%%ymm0,%%ymm0		\n\t	vaddpd	%%ymm13,%%ymm8 ,%%ymm8 		\n\t"\
		"vsubpd	%%ymm4,%%ymm1,%%ymm1		\n\t	vsubpd	%%ymm12,%%ymm9 ,%%ymm9 		\n\t"\
		"vmovaps	%%ymm2,    (%%rbx)		\n\t	vmovaps	%%ymm10,    (%%r11)	\n\t"\
		"vmovaps	%%ymm3,0x20(%%rbx)		\n\t	vmovaps	%%ymm11,0x20(%%r11)	\n\t"\
		"vmovaps	%%ymm0,    (%%rcx)		\n\t	vmovaps	%%ymm8 ,    (%%r12)	\n\t"\
		"vmovaps	%%ymm1,0x20(%%rcx)		\n\t	vmovaps	%%ymm9 ,0x20(%%r12)	\n\t"\
		:					/* outputs: none */\
		: [__cc0] "m" (Xcc0)	/* All inputs from memory addresses here */\
		 ,[__i0] "m" (Xi0)\
		 ,[__i1] "m" (Xi1)\
		 ,[__i2] "m" (Xi2)\
		 ,[__o0] "m" (Xo0)\
		 ,[__o1] "m" (Xo1)\
		 ,[__o2] "m" (Xo2)\
		 ,[__j0] "m" (Xj0)\
		 ,[__j1] "m" (Xj1)\
		 ,[__j2] "m" (Xj2)\
		 ,[__u0] "m" (Xu0)\
		 ,[__u1] "m" (Xu1)\
		 ,[__u2] "m" (Xu2)\
		: "cc","memory","rax","rbx","rcx","rdx","r10","r11","r12","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"		/* Clobbered registers */\
	);\
	}

	#define SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(Xadd0, Xadd1, Xadd2, Xadd3, Xtmp, Xstride)\
	{\
	__asm__ volatile (\
		"movq	%[__tmp]   ,%%rax	\n\t"\
		"movq	%[__stride],%%rsi	\n\t"\
		"movq	%%rax,%%rbx			\n\t"\
		"addq	%%rsi,%%rbx			/* add_in1  */\n\t"\
		"shlq	$1,%%rsi			/* stride*2 */\n\t"\
		"vmovaps	    (%%rax),%%ymm0	\n\t"\
		"vmovaps	    (%%rbx),%%ymm2	\n\t"\
		"vmovaps	0x20(%%rax),%%ymm1	\n\t"\
		"vmovaps	0x20(%%rbx),%%ymm3	\n\t"\
		"vmovaps	    (%%rax),%%ymm4	\n\t"\
		"vmovaps	    (%%rbx),%%ymm6	\n\t"\
		"vmovaps	0x20(%%rax),%%ymm5	\n\t"\
		"vmovaps	0x20(%%rbx),%%ymm7	\n\t"\
		"addq	%%rsi,%%rax			/* add_in2  */\n\t"\
		"addq	%%rsi,%%rbx			/* add_in3  */\n\t"\
		"vaddpd	    (%%rax),%%ymm0,%%ymm0	\n\t"\
		"vaddpd	    (%%rbx),%%ymm2,%%ymm2	\n\t"\
		"vaddpd	0x20(%%rax),%%ymm1,%%ymm1	\n\t"\
		"vaddpd	0x20(%%rbx),%%ymm3,%%ymm3	\n\t"\
		"vsubpd	    (%%rax),%%ymm4,%%ymm4	\n\t"\
		"vsubpd	    (%%rbx),%%ymm6,%%ymm6	\n\t"\
		"vsubpd	0x20(%%rax),%%ymm5,%%ymm5	\n\t"\
		"vsubpd	0x20(%%rbx),%%ymm7,%%ymm7	\n\t"\
		"/* Finish radix-4 butterfly and store results into main-array slots: */\n\t"\
		"movq	%[__add0],%%rax		\n\t"\
		"movq	%[__add1],%%rbx		\n\t"\
		"movq	%[__add2],%%rcx		\n\t"\
		"movq	%[__add3],%%rdx		\n\t"\
		"vsubpd	%%ymm2,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	%%ymm7,%%ymm4,%%ymm4		\n\t"\
		"vsubpd	%%ymm3,%%ymm1,%%ymm1		\n\t"\
		"vsubpd	%%ymm6,%%ymm5,%%ymm5		\n\t"\
		"vmovaps	%%ymm0,    (%%rbx)	\n\t"\
		"vmovaps	%%ymm4,    (%%rcx)	\n\t"\
		"vmovaps	%%ymm1,0x20(%%rbx)	\n\t"\
		"vmovaps	%%ymm5,0x20(%%rdx)	\n\t"\
		"vaddpd	%%ymm2,%%ymm2,%%ymm2		\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	%%ymm3,%%ymm3,%%ymm3		\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	%%ymm0,%%ymm2,%%ymm2		\n\t"\
		"vaddpd	%%ymm4,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	%%ymm1,%%ymm3,%%ymm3		\n\t"\
		"vaddpd	%%ymm5,%%ymm6,%%ymm6		\n\t"\
		"vmovaps	%%ymm2,    (%%rax)	\n\t"\
		"vmovaps	%%ymm7,    (%%rdx)	\n\t"\
		"vmovaps	%%ymm3,0x20(%%rax)	\n\t"\
		"vmovaps	%%ymm6,0x20(%%rcx)	\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__add1] "m" (Xadd1)\
		 ,[__add2] "m" (Xadd2)\
		 ,[__add3] "m" (Xadd3)\
		 ,[__tmp] "m" (Xtmp)\
		 ,[__stride] "e" (Xstride)\
		: "cc","memory","rax","rbx","rcx","rdx","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
	);\
	}

	/* DIF radix-4 subconvolution, sans twiddles, inputs in __i0-3, outputs in __o0-3, possibly coincident with inputs: */
	#define SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_E(Xi0,Xi1,Xi2,Xi3, Xo0,Xo1,Xo2,Xo3)\
	{\
	__asm__ volatile (\
		"movq	%[__i0],%%rax		\n\t"\
		"movq	%[__i1],%%rbx		\n\t"\
		"movq	%[__i2],%%rcx		\n\t"\
		"movq	%[__i3],%%rdx		\n\t"\
		"vmovaps	    (%%rax),%%ymm0		\n\t"\
		"vmovaps	    (%%rbx),%%ymm4		\n\t"\
		"vmovaps	0x20(%%rax),%%ymm1		\n\t"\
		"vmovaps	0x20(%%rbx),%%ymm5		\n\t"\
		"vmovaps	%%ymm0,%%ymm2			\n\t"\
		"vmovaps	%%ymm4,%%ymm6			\n\t"\
		"vmovaps	%%ymm1,%%ymm3			\n\t"\
		"vmovaps	%%ymm5,%%ymm7			\n\t"\
		"vaddpd	    (%%rcx),%%ymm0,%%ymm0	\n\t"\
		"vaddpd	    (%%rdx),%%ymm4,%%ymm4	\n\t"\
		"vaddpd	0x20(%%rcx),%%ymm1,%%ymm1	\n\t"\
		"vaddpd	0x20(%%rdx),%%ymm5,%%ymm5	\n\t"\
		"vsubpd	    (%%rcx),%%ymm2,%%ymm2	\n\t"\
		"vsubpd	    (%%rdx),%%ymm6,%%ymm6	\n\t"\
		"vsubpd	0x20(%%rcx),%%ymm3,%%ymm3	\n\t"\
		"vsubpd	0x20(%%rdx),%%ymm7,%%ymm7	\n\t"\
		"/* Finish radix-4 butterfly and store results into main-array slots: */\n\t"\
		"movq	%[__o0],%%rax		\n\t"\
		"movq	%[__o1],%%rbx		\n\t"\
		"movq	%[__o2],%%rcx		\n\t"\
		"movq	%[__o3],%%rdx		\n\t"\
		"vsubpd	%%ymm4,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	%%ymm7,%%ymm2,%%ymm2		\n\t"\
		"vsubpd	%%ymm5,%%ymm1,%%ymm1		\n\t"\
		"vsubpd	%%ymm6,%%ymm3,%%ymm3		\n\t"\
		"vmovaps	%%ymm0,    (%%rbx)		\n\t"\
		"vmovaps	%%ymm2,    (%%rcx)		\n\t"\
		"vmovaps	%%ymm1,0x20(%%rbx)		\n\t"\
		"vmovaps	%%ymm3,0x20(%%rdx)		\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	%%ymm0,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	%%ymm2,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	%%ymm1,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	%%ymm3,%%ymm6,%%ymm6		\n\t"\
		"vmovaps	%%ymm4,    (%%rax)		\n\t"\
		"vmovaps	%%ymm7,    (%%rdx)		\n\t"\
		"vmovaps	%%ymm5,0x20(%%rax)		\n\t"\
		"vmovaps	%%ymm6,0x20(%%rcx)		\n\t"\
		:					/* outputs: none */\
		: [__i0] "m" (Xi0)	/* All inputs from memory addresses here */\
		 ,[__i1] "m" (Xi1)\
		 ,[__i2] "m" (Xi2)\
		 ,[__i3] "m" (Xi3)\
		 ,[__o0] "m" (Xo0)\
		 ,[__o1] "m" (Xo1)\
		 ,[__o2] "m" (Xo2)\
		 ,[__o3] "m" (Xo3)\
		: "cc","memory","rax","rbx","rcx","rdx","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
	);\
	}

	#define SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(Xadd0, Xadd1, Xadd2, Xadd3, Xtmp, Xstride)\
	{\
	__asm__ volatile (\
		"movq	%[__add0],%%rax		\n\t"\
		"movq	%[__add1],%%rbx		\n\t"\
		"movq	%[__add2],%%rcx		\n\t"\
		"movq	%[__add3],%%rdx		\n\t"\
		"vmovaps	    (%%rax),%%ymm0	\n\t"\
		"vmovaps	    (%%rcx),%%ymm4	\n\t"\
		"vmovaps	0x20(%%rax),%%ymm1	\n\t"\
		"vmovaps	0x20(%%rcx),%%ymm5	\n\t"\
		"vmovaps	%%ymm0,%%ymm2			\n\t"\
		"vmovaps	%%ymm4,%%ymm6			\n\t"\
		"vmovaps	%%ymm1,%%ymm3			\n\t"\
		"vmovaps	%%ymm5,%%ymm7			\n\t"\
		"movq	%[__tmp]   ,%%rax	\n\t"\
		"movq	%[__stride],%%rcx	\n\t"\
		"vaddpd	    (%%rbx),%%ymm0,%%ymm0	\n\t"\
		"vaddpd	    (%%rdx),%%ymm4,%%ymm4	\n\t"\
		"vaddpd	0x20(%%rbx),%%ymm1,%%ymm1	\n\t"\
		"vaddpd	0x20(%%rdx),%%ymm5,%%ymm5	\n\t"\
		"vsubpd	    (%%rbx),%%ymm2,%%ymm2	\n\t"\
		"vsubpd	    (%%rdx),%%ymm6,%%ymm6	\n\t"\
		"vsubpd	0x20(%%rbx),%%ymm3,%%ymm3	\n\t"\
		"vsubpd	0x20(%%rdx),%%ymm7,%%ymm7	\n\t"\
		"movq	%%rax,%%rbx			\n\t"\
		"addq	%%rcx,%%rbx			\n\t"\
		"movq	%%rbx,%%rdx			\n\t"\
		"addq	%%rcx,%%rcx			\n\t"\
		"addq	%%rcx,%%rdx			\n\t"\
		"addq	%%rax,%%rcx			\n\t"\
		"/* Finish radix-4 butterfly and store results into temp-array slots: */\n\t"\
		"vsubpd	%%ymm4,%%ymm0,%%ymm0			\n\t"\
		"vsubpd	%%ymm7,%%ymm2,%%ymm2			\n\t"\
		"vsubpd	%%ymm5,%%ymm1,%%ymm1			\n\t"\
		"vsubpd	%%ymm6,%%ymm3,%%ymm3			\n\t"\
		"vmovaps	%%ymm0,    (%%rcx)	\n\t"\
		"vmovaps	%%ymm2,    (%%rdx)	\n\t"\
		"vmovaps	%%ymm1,0x20(%%rcx)	\n\t"\
		"vmovaps	%%ymm3,0x20(%%rbx)	\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4			\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7			\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5			\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6			\n\t"\
		"vaddpd	%%ymm0,%%ymm4,%%ymm4			\n\t"\
		"vaddpd	%%ymm2,%%ymm7,%%ymm7			\n\t"\
		"vaddpd	%%ymm1,%%ymm5,%%ymm5			\n\t"\
		"vaddpd	%%ymm3,%%ymm6,%%ymm6			\n\t"\
		"vmovaps	%%ymm4,    (%%rax)	\n\t"\
		"vmovaps	%%ymm7,    (%%rbx)	\n\t"\
		"vmovaps	%%ymm5,0x20(%%rax)	\n\t"\
		"vmovaps	%%ymm6,0x20(%%rdx)	\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__add1] "m" (Xadd1)\
		 ,[__add2] "m" (Xadd2)\
		 ,[__add3] "m" (Xadd3)\
		 ,[__tmp] "m" (Xtmp)\
		 ,[__stride] "e" (Xstride)\
		: "cc","memory","rax","rbx","rcx","rdx","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
	);\
	}

	/* DIT radix-4 subconvolution, sans twiddles, inputs in __i0-3, outputs in __o0-3, possibly coincident with inputs: */
	#define SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_E(Xi0,Xi1,Xi2,Xi3, Xo0,Xo1,Xo2,Xo3)\
	{\
	__asm__ volatile (\
		"movq	%[__i0],%%rax		\n\t"\
		"movq	%[__i1],%%rbx		\n\t"\
		"movq	%[__i2],%%rcx		\n\t"\
		"movq	%[__i3],%%rdx		\n\t"\
		"vmovaps	    (%%rax),%%ymm0		\n\t"\
		"vmovaps	    (%%rcx),%%ymm4		\n\t"\
		"vmovaps	0x20(%%rax),%%ymm1		\n\t"\
		"vmovaps	0x20(%%rcx),%%ymm5		\n\t"\
		"vmovaps	%%ymm0,%%ymm2			\n\t"\
		"vmovaps	%%ymm4,%%ymm6			\n\t"\
		"vmovaps	%%ymm1,%%ymm3			\n\t"\
		"vmovaps	%%ymm5,%%ymm7			\n\t"\
		"vaddpd	    (%%rbx),%%ymm0,%%ymm0	\n\t"\
		"vaddpd	    (%%rdx),%%ymm4,%%ymm4	\n\t"\
		"vaddpd	0x20(%%rbx),%%ymm1,%%ymm1	\n\t"\
		"vaddpd	0x20(%%rdx),%%ymm5,%%ymm5	\n\t"\
		"vsubpd	    (%%rbx),%%ymm2,%%ymm2	\n\t"\
		"vsubpd	    (%%rdx),%%ymm6,%%ymm6	\n\t"\
		"vsubpd	0x20(%%rbx),%%ymm3,%%ymm3	\n\t"\
		"vsubpd	0x20(%%rdx),%%ymm7,%%ymm7	\n\t"\
		"/* Finish radix-4 butterfly and store results into output-array slots: */\n\t"\
		"movq	%[__o0],%%rax		\n\t"\
		"movq	%[__o1],%%rbx		\n\t"\
		"movq	%[__o2],%%rcx		\n\t"\
		"movq	%[__o3],%%rdx		\n\t"\
		"vsubpd	%%ymm4,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	%%ymm7,%%ymm2,%%ymm2		\n\t"\
		"vsubpd	%%ymm5,%%ymm1,%%ymm1		\n\t"\
		"vsubpd	%%ymm6,%%ymm3,%%ymm3		\n\t"\
		"vmovaps	%%ymm0,    (%%rcx)		\n\t"\
		"vmovaps	%%ymm2,    (%%rdx)		\n\t"\
		"vmovaps	%%ymm1,0x20(%%rcx)		\n\t"\
		"vmovaps	%%ymm3,0x20(%%rbx)		\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	%%ymm0,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	%%ymm2,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	%%ymm1,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	%%ymm3,%%ymm6,%%ymm6		\n\t"\
		"vmovaps	%%ymm4,    (%%rax)		\n\t"\
		"vmovaps	%%ymm7,    (%%rbx)		\n\t"\
		"vmovaps	%%ymm5,0x20(%%rax)		\n\t"\
		"vmovaps	%%ymm6,0x20(%%rdx)		\n\t"\
		:					/* outputs: none */\
		: [__i0] "m" (Xi0)	/* All inputs from memory addresses here */\
		 ,[__i1] "m" (Xi1)\
		 ,[__i2] "m" (Xi2)\
		 ,[__i3] "m" (Xi3)\
		 ,[__o0] "m" (Xo0)\
		 ,[__o1] "m" (Xo1)\
		 ,[__o2] "m" (Xo2)\
		 ,[__o3] "m" (Xo3)\
		: "cc","memory","rax","rbx","rcx","rdx","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
	);\
	}

	#define SSE2_RADIX_05_DFT_0TWID_X2(Xcc1, Xi0,Xi1,Xi2,Xi3,Xi4, Xo0,Xo1,Xo2,Xo3,Xo4, Yi0,Yi1,Yi2,Yi3,Yi4, Yo0,Yo1,Yo2,Yo3,Yo4)\
	{\
	__asm__ volatile (\
		"movq	%[__i0],%%rsi				\n\t"\
		"movq	%[__i1],%%rax				\n\t"\
		"movq	%[__i2],%%rbx				\n\t"\
		"movq	%[__i3],%%rcx				\n\t"\
		"movq	%[__i4],%%rdx				\n\t"\
		"movq	%[__o0],%%rdi				\n\t"\
		"vmovaps	    (%%rax),%%ymm0		\n\t		movq	%[__i5],%%r10				\n\t"\
		"vmovaps	0x20(%%rax),%%ymm1		\n\t		movq	%[__i6],%%r11				\n\t"\
		"vmovaps	    (%%rbx),%%ymm2		\n\t		movq	%[__i7],%%r12				\n\t"\
		"vmovaps	0x20(%%rbx),%%ymm3		\n\t		movq	%[__i8],%%r13				\n\t"\
		"vmovaps	    (%%rcx),%%ymm4		\n\t		movq	%[__i9],%%r14				\n\t"\
		"vmovaps	0x20(%%rcx),%%ymm5		\n\t		movq	%[__o5],%%r15				\n\t"\
		"vmovaps	    (%%rdx),%%ymm6		\n\t		vmovaps	    (%%r11),%%ymm8 			\n\t"\
		"vmovaps	0x20(%%rdx),%%ymm7		\n\t		vmovaps	0x20(%%r11),%%ymm9 			\n\t"\
		"vsubpd	%%ymm6,%%ymm0,%%ymm0		\n\t		vmovaps	    (%%r12),%%ymm10			\n\t"\
		"vsubpd	%%ymm7,%%ymm1,%%ymm1		\n\t		vmovaps	0x20(%%r12),%%ymm11			\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6		\n\t		vmovaps	    (%%r13),%%ymm12			\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7		\n\t		vmovaps	0x20(%%r13),%%ymm13			\n\t"\
		"vaddpd	%%ymm0,%%ymm6,%%ymm6		\n\t		vmovaps	    (%%r14),%%ymm14			\n\t"\
		"vaddpd	%%ymm1,%%ymm7,%%ymm7		\n\t		vmovaps	0x20(%%r14),%%ymm15			\n\t"\
		"vsubpd	%%ymm4,%%ymm2,%%ymm2		\n\t		vsubpd	%%ymm14,%%ymm8 ,%%ymm8 		\n\t"\
		"vsubpd	%%ymm5,%%ymm3,%%ymm3		\n\t		vsubpd	%%ymm15,%%ymm9 ,%%ymm9 		\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4		\n\t		vaddpd	%%ymm14,%%ymm14,%%ymm14		\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5		\n\t		vaddpd	%%ymm15,%%ymm15,%%ymm15		\n\t"\
		"vaddpd	%%ymm2,%%ymm4,%%ymm4		\n\t		vaddpd	%%ymm8 ,%%ymm14,%%ymm14		\n\t"\
		"vaddpd	%%ymm3,%%ymm5,%%ymm5		\n\t		vaddpd	%%ymm9 ,%%ymm15,%%ymm15		\n\t"\
		"movq	%[__cc1],%%rax		/* Sincos data shared by both streams */	\n\t"\
		"vsubpd	%%ymm4,%%ymm6,%%ymm6		\n\t		vsubpd	%%ymm12,%%ymm10,%%ymm10		\n\t"\
		"vsubpd	%%ymm5,%%ymm7,%%ymm7		\n\t		vsubpd	%%ymm13,%%ymm11,%%ymm11		\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4		\n\t		vaddpd	%%ymm12,%%ymm12,%%ymm12		\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5		\n\t		vaddpd	%%ymm13,%%ymm13,%%ymm13		\n\t"\
		"vaddpd	%%ymm6,%%ymm4,%%ymm4		\n\t		vaddpd	%%ymm10,%%ymm12,%%ymm12		\n\t"\
		"vaddpd	%%ymm7,%%ymm5,%%ymm5		\n\t		vaddpd	%%ymm11,%%ymm13,%%ymm13		\n\t"\
		"vaddpd	    (%%rsi),%%ymm4,%%ymm4	\n\t		vsubpd	%%ymm12,%%ymm14,%%ymm14		\n\t"\
		"vaddpd	0x20(%%rsi),%%ymm5,%%ymm5	\n\t		vsubpd	%%ymm13,%%ymm15,%%ymm15		\n\t"\
		"vmovaps	%%ymm4,    (%%rdi)		\n\t		vaddpd	%%ymm12,%%ymm12,%%ymm12		\n\t"\
		"vmovaps	%%ymm5,0x20(%%rdi)		\n\t		vaddpd	%%ymm13,%%ymm13,%%ymm13		\n\t"\
		"vmulpd	0x20(%%rax),%%ymm6,%%ymm6	\n\t		vaddpd	%%ymm14,%%ymm12,%%ymm12		\n\t"\
		"vmulpd	0x20(%%rax),%%ymm7,%%ymm7	\n\t		vaddpd	%%ymm15,%%ymm13,%%ymm13		\n\t"\
		"vsubpd	    (%%rsi),%%ymm4,%%ymm4	\n\t		vaddpd	    (%%r10),%%ymm12,%%ymm12	\n\t"\
		"vsubpd	0x20(%%rsi),%%ymm5,%%ymm5	\n\t		vaddpd	0x20(%%r10),%%ymm13,%%ymm13	\n\t"\
		"vmulpd	    (%%rax),%%ymm4,%%ymm4	\n\t		vmovaps	%%ymm12,    (%%r15)			\n\t"\
		"vmulpd	    (%%rax),%%ymm5,%%ymm5	\n\t		vmovaps	%%ymm13,0x20(%%r15)			\n\t"\
		"vaddpd	    (%%rdi),%%ymm4,%%ymm4	\n\t		vmulpd	0x20(%%rax),%%ymm14,%%ymm14	\n\t"\
		"vaddpd	0x20(%%rdi),%%ymm5,%%ymm5	\n\t		vmulpd	0x20(%%rax),%%ymm15,%%ymm15	\n\t"\
		"vsubpd	%%ymm6,%%ymm4,%%ymm4		\n\t		vsubpd	    (%%r10),%%ymm12,%%ymm12	\n\t"\
		"vsubpd	%%ymm7,%%ymm5,%%ymm5		\n\t		vsubpd	0x20(%%r10),%%ymm13,%%ymm13	\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6		\n\t		vmulpd	    (%%rax),%%ymm12,%%ymm12	\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7		\n\t		vmulpd	    (%%rax),%%ymm13,%%ymm13	\n\t"\
		"vaddpd	%%ymm4,%%ymm6,%%ymm6		\n\t		vaddpd	    (%%r15),%%ymm12,%%ymm12	\n\t"\
		"vaddpd	%%ymm5,%%ymm7,%%ymm7		\n\t		vaddpd	0x20(%%r15),%%ymm13,%%ymm13	\n\t"\
		"vmovaps	%%ymm4,    (%%rsi)		\n\t		vsubpd	%%ymm14,%%ymm12,%%ymm12		\n\t"\
		"vmovaps	%%ymm5,0x20(%%rsi)		\n\t		vsubpd	%%ymm15,%%ymm13,%%ymm13		\n\t"\
		"vmovaps	%%ymm0,%%ymm4			\n\t		vaddpd	%%ymm14,%%ymm14,%%ymm14		\n\t"\
		"vmovaps	%%ymm1,%%ymm5			\n\t		vaddpd	%%ymm15,%%ymm15,%%ymm15		\n\t"\
		"vsubpd	%%ymm2     ,%%ymm0,%%ymm0	\n\t		vaddpd	%%ymm12,%%ymm14,%%ymm14		\n\t"\
		"vsubpd	%%ymm3     ,%%ymm1,%%ymm1	\n\t		vaddpd	%%ymm13,%%ymm15,%%ymm15		\n\t"\
		"vmulpd	0x40(%%rax),%%ymm0,%%ymm0	\n\t		vmovaps	%%ymm12,    (%%r10)			\n\t"\
		"vmulpd	0x40(%%rax),%%ymm1,%%ymm1	\n\t		vmovaps	%%ymm13,0x20(%%r10)			\n\t"\
		"vmulpd	0x60(%%rax),%%ymm2,%%ymm2	\n\t		vmovaps	%%ymm8 ,%%ymm12				\n\t"\
		"vmulpd	0x60(%%rax),%%ymm3,%%ymm3	\n\t		vmovaps	%%ymm9 ,%%ymm13				\n\t"\
		"vmulpd	0x80(%%rax),%%ymm4,%%ymm4	\n\t		vsubpd	%%ymm10,%%ymm8 ,%%ymm8 		\n\t"\
		"vmulpd	0x80(%%rax),%%ymm5,%%ymm5	\n\t		vsubpd	%%ymm11,%%ymm9 ,%%ymm9 		\n\t"\
		"vaddpd	%%ymm0     ,%%ymm2,%%ymm2	\n\t		vmulpd	0x40(%%rax),%%ymm8 ,%%ymm8 	\n\t"\
		"vaddpd	%%ymm1     ,%%ymm3,%%ymm3	\n\t		vmulpd	0x40(%%rax),%%ymm9 ,%%ymm9 	\n\t"\
		"vsubpd	%%ymm4     ,%%ymm0,%%ymm0	\n\t		vmulpd	0x60(%%rax),%%ymm10,%%ymm10	\n\t"\
		"vsubpd	%%ymm5     ,%%ymm1,%%ymm1	\n\t		vmulpd	0x60(%%rax),%%ymm11,%%ymm11	\n\t"\
		"vmovaps	    (%%rsi),%%ymm4		\n\t		vmulpd	0x80(%%rax),%%ymm12,%%ymm12	\n\t"\
		"vmovaps	0x20(%%rsi),%%ymm5		\n\t		vmulpd	0x80(%%rax),%%ymm13,%%ymm13	\n\t"\
		"/*------ End of shared-trig-data accessed via rax register ------*/	\n\t"\
		"movq	%[__o1],%%rax				\n\t		vaddpd	%%ymm8 ,%%ymm10,%%ymm10		\n\t"\
		"movq	%[__o4],%%rdx				\n\t		vaddpd	%%ymm9 ,%%ymm11,%%ymm11		\n\t"\
		"vsubpd	%%ymm3,%%ymm6,%%ymm6		\n\t		vsubpd	%%ymm12,%%ymm8 ,%%ymm8 		\n\t"\
		"vsubpd	%%ymm2,%%ymm7,%%ymm7		\n\t		vsubpd	%%ymm13,%%ymm9 ,%%ymm9 		\n\t"\
		"vaddpd	%%ymm3,%%ymm3,%%ymm3		\n\t		vmovaps	    (%%r10),%%ymm12			\n\t"\
		"vaddpd	%%ymm2,%%ymm2,%%ymm2		\n\t		vmovaps	0x20(%%r10),%%ymm13			\n\t"\
		"vmovaps	%%ymm6,    (%%rax)		\n\t		movq	%[__o6],%%r11				\n\t"\
		"vmovaps	%%ymm7,0x20(%%rdx)		\n\t		movq	%[__o9],%%r14				\n\t"\
		"vaddpd	%%ymm6,%%ymm3,%%ymm3		\n\t		vsubpd	%%ymm11,%%ymm14,%%ymm14		\n\t"\
		"vaddpd	%%ymm7,%%ymm2,%%ymm2		\n\t		vsubpd	%%ymm10,%%ymm15,%%ymm15		\n\t"\
		"vmovaps	%%ymm3,    (%%rdx)		\n\t		vaddpd	%%ymm11,%%ymm11,%%ymm11		\n\t"\
		"vmovaps	%%ymm2,0x20(%%rax)		\n\t		vaddpd	%%ymm10,%%ymm10,%%ymm10		\n\t"\
		"movq	%[__o2],%%rbx				\n\t		vmovaps	%%ymm14,    (%%r11)			\n\t"\
		"movq	%[__o3],%%rcx				\n\t		vmovaps	%%ymm15,0x20(%%r14)			\n\t"\
		"vsubpd	%%ymm1,%%ymm4,%%ymm4		\n\t		vaddpd	%%ymm14,%%ymm11,%%ymm11		\n\t"\
		"vsubpd	%%ymm0,%%ymm5,%%ymm5		\n\t		vaddpd	%%ymm15,%%ymm10,%%ymm10		\n\t"\
		"vaddpd	%%ymm1,%%ymm1,%%ymm1		\n\t		vmovaps	%%ymm11,    (%%r14)			\n\t"\
		"vaddpd	%%ymm0,%%ymm0,%%ymm0		\n\t		vmovaps	%%ymm10,0x20(%%r11)			\n\t"\
		"vmovaps	%%ymm4,    (%%rbx)		\n\t		movq	%[__o7],%%r12				\n\t"\
		"vmovaps	%%ymm5,0x20(%%rcx)		\n\t		movq	%[__o8],%%r13				\n\t"\
		"vaddpd	%%ymm4,%%ymm1,%%ymm1		\n\t		vsubpd	%%ymm9 ,%%ymm12,%%ymm12		\n\t"\
		"vaddpd	%%ymm5,%%ymm0,%%ymm0		\n\t		vsubpd	%%ymm8 ,%%ymm13,%%ymm13		\n\t"\
		"vmovaps	%%ymm1,    (%%rcx)		\n\t		vaddpd	%%ymm9 ,%%ymm9 ,%%ymm9 		\n\t"\
		"vmovaps	%%ymm0,0x20(%%rbx)		\n\t		vaddpd	%%ymm8 ,%%ymm8 ,%%ymm8 		\n\t"\
		"												vmovaps	%%ymm12,    (%%r12)			\n\t"\
		"												vmovaps	%%ymm13,0x20(%%r13)			\n\t"\
		"												vaddpd	%%ymm12,%%ymm9 ,%%ymm9 		\n\t"\
		"												vaddpd	%%ymm13,%%ymm8 ,%%ymm8 		\n\t"\
		"												vmovaps	%%ymm9 ,    (%%r13)			\n\t"\
		"												vmovaps	%%ymm8 ,0x20(%%r12)			\n\t"\
		:					/* outputs: none */\
		: [__cc1] "m" (Xcc1)	/* All inputs from memory addresses here */\
		 ,[__i0] "m" (Xi0)\
		 ,[__i1] "m" (Xi1)\
		 ,[__i2] "m" (Xi2)\
		 ,[__i3] "m" (Xi3)\
		 ,[__i4] "m" (Xi4)\
		 ,[__o0] "m" (Xo0)\
		 ,[__o1] "m" (Xo1)\
		 ,[__o2] "m" (Xo2)\
		 ,[__o3] "m" (Xo3)\
		 ,[__o4] "m" (Xo4)\
		 ,[__i5] "m" (Yi0)\
		 ,[__i6] "m" (Yi1)\
		 ,[__i7] "m" (Yi2)\
		 ,[__i8] "m" (Yi3)\
		 ,[__i9] "m" (Yi4)\
		 ,[__o5] "m" (Yo0)\
		 ,[__o6] "m" (Yo1)\
		 ,[__o7] "m" (Yo2)\
		 ,[__o8] "m" (Yo3)\
		 ,[__o9] "m" (Yo4)\
		: "cc","memory","rax","rbx","rcx","rdx","rsi","rdi","r10","r11","r12","r13","r14","r15","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"		/* Clobbered registers */\
	);\
	}

	/*...Radix-7 DFT: Inputs in memlocs __i0-6, outputs into __o0-6, possibly coincident with inputs:\ */\
	#define SSE2_RADIX_07_DFT(Xi0,Xi1,Xi2,Xi3,Xi4,Xi5,Xi6, Xcc, Xo0,Xo1,Xo2,Xo3,Xo4,Xo5,Xo6)\
	{\
	__asm__ volatile (\
		"movq	%[__i1],%%rax		\n\t"\
		"movq	%[__i2],%%rbx		\n\t"\
		"movq	%[__i3],%%rcx		\n\t"\
		"movq	%[__i4],%%rdx		\n\t"\
		"movq	%[__i5],%%rsi		\n\t"\
		"movq	%[__i6],%%rdi		\n\t	/*** Imaginary Parts: ***/	\n\t"\
		"vmovaps	(%%rax),%%ymm6		\n\t	vmovaps	0x20(%%rax),%%ymm14	\n\t"\
		"vmovaps	(%%rdi),%%ymm1		\n\t	vmovaps	0x20(%%rdi),%%ymm9 	\n\t"\
		"vmovaps	(%%rbx),%%ymm5		\n\t	vmovaps	0x20(%%rbx),%%ymm13	\n\t"\
		"vmovaps	(%%rsi),%%ymm2		\n\t	vmovaps	0x20(%%rsi),%%ymm10	\n\t"\
		"vmovaps	(%%rcx),%%ymm4		\n\t	vmovaps	0x20(%%rcx),%%ymm12	\n\t"\
		"vmovaps	(%%rdx),%%ymm3		\n\t	vmovaps	0x20(%%rdx),%%ymm11	\n\t"\
		"movq	%[__i0],%%rbx		\n\t"\
		"vsubpd	%%ymm1,%%ymm6,%%ymm6	\n\t	vsubpd	%%ymm9 ,%%ymm14,%%ymm14		\n\t"\
		"vaddpd	%%ymm1,%%ymm1,%%ymm1	\n\t	vaddpd	%%ymm9 ,%%ymm9 ,%%ymm9 		\n\t"\
		"vaddpd	%%ymm6,%%ymm1,%%ymm1	\n\t	vaddpd	%%ymm14,%%ymm9 ,%%ymm9  	\n\t"\
		"vsubpd	%%ymm2,%%ymm5,%%ymm5	\n\t	vsubpd	%%ymm10,%%ymm13,%%ymm13		\n\t"\
		"vaddpd	%%ymm2,%%ymm2,%%ymm2	\n\t	vaddpd	%%ymm10,%%ymm10,%%ymm10		\n\t"\
		"vaddpd	%%ymm5,%%ymm2,%%ymm2	\n\t	vaddpd	%%ymm13,%%ymm10,%%ymm10		\n\t"\
		"vmovaps	(%%rbx),%%ymm0		\n\t	vmovaps	0x20(%%rbx),%%ymm8 	\n\t"\
		"vsubpd	%%ymm3,%%ymm4,%%ymm4	\n\t	vsubpd	%%ymm11,%%ymm12,%%ymm12		\n\t"\
		"vaddpd	%%ymm3,%%ymm3,%%ymm3	\n\t	vaddpd	%%ymm11,%%ymm11,%%ymm11		\n\t"\
		"vaddpd	%%ymm4,%%ymm3,%%ymm3	\n\t	vaddpd	%%ymm12,%%ymm11,%%ymm11		\n\t"\
		"\n\t"\
		"movq	%[__o0],%%rcx		\n\t"\
		"movq	%[__cc],%%rsi		\n\t"\
		"vmovaps	%%ymm0,0x100(%%rsi)	\n\t	vmovaps	%%ymm8 ,0x140(%%rsi)	\n\t"\
		"vmovaps	%%ymm6,0x120(%%rsi)	\n\t	vmovaps	%%ymm14,0x160(%%rsi)	\n\t"\
		"vaddpd	%%ymm1,%%ymm0,%%ymm0	\n\t	vaddpd	%%ymm9 ,%%ymm8 ,%%ymm8 	\n\t"\
		"vmovaps	%%ymm5,%%ymm7		\n\t	vmovaps	%%ymm13,%%ymm15		\n\t"\
		"vaddpd	%%ymm2,%%ymm3,%%ymm3	\n\t	vaddpd	%%ymm10,%%ymm11,%%ymm11		\n\t"\
		"vsubpd	%%ymm4,%%ymm5,%%ymm5	\n\t	vsubpd	%%ymm12,%%ymm13,%%ymm13		\n\t"\
		"vsubpd	%%ymm2,%%ymm1,%%ymm1	\n\t	vsubpd	%%ymm10,%%ymm9 ,%%ymm9  	\n\t"\
		"vsubpd	%%ymm7,%%ymm6,%%ymm6	\n\t	vsubpd	%%ymm15,%%ymm14,%%ymm14		\n\t"\
		"vaddpd	%%ymm2,%%ymm2,%%ymm2	\n\t	vaddpd	%%ymm10,%%ymm10,%%ymm10		\n\t"\
		"vaddpd	%%ymm7,%%ymm4,%%ymm4	\n\t	vaddpd	%%ymm15,%%ymm12,%%ymm12		\n\t"\
		"vaddpd	%%ymm3,%%ymm0,%%ymm0	\n\t	vaddpd	%%ymm11,%%ymm8 ,%%ymm8  	\n\t"\
	"vaddpd	0x120(%%rsi),%%ymm5,%%ymm5	\n\t	vaddpd	0x160(%%rsi),%%ymm13,%%ymm13\n\t"\
		"vsubpd	%%ymm2,%%ymm3,%%ymm3	\n\t	vsubpd	%%ymm10,%%ymm11,%%ymm11		\n\t"\
		"vmovaps	%%ymm4,%%ymm7		\n\t	vmovaps	%%ymm12,%%ymm15		\n\t"\
		"vmovaps	%%ymm0,    (%%rcx)	\n\t	vmovaps	%%ymm8 ,0x20(%%rcx)	\n\t"/* B0 */\
		"vsubpd	%%ymm6,%%ymm4,%%ymm4	\n\t	vsubpd		%%ymm14,%%ymm12,%%ymm12		\n\t"\
		"vmovaps	%%ymm1,%%ymm2		\n\t	vmovaps	%%ymm9 ,%%ymm10		\n\t"\
	"vsubpd 0x100(%%rsi),%%ymm0,%%ymm0	\n\t	vsubpd 0x140(%%rsi),%%ymm8 ,%%ymm8  \n\t"\
	"vmulpd	0x20(%%rsi),%%ymm5,%%ymm5	\n\t	vmulpd	0x20(%%rsi),%%ymm13,%%ymm13	\n\t"\
	"vaddpd	%%ymm3,%%ymm2,%%ymm2		\n\t	vaddpd	%%ymm11,%%ymm10,%%ymm10		\n\t"\
	"vmulpd	0x80(%%rsi),%%ymm3,%%ymm3	\n\t	vmulpd	0x80(%%rsi),%%ymm11,%%ymm11	\n\t"\
	"vmulpd	0xe0(%%rsi),%%ymm4,%%ymm4	\n\t	vmulpd	0xe0(%%rsi),%%ymm12,%%ymm12	\n\t"\
	"vmulpd	0x40(%%rsi),%%ymm1,%%ymm1	\n\t	vmulpd	0x40(%%rsi),%%ymm9 ,%%ymm9  \n\t"\
	"vmulpd	0x60(%%rsi),%%ymm6,%%ymm6	\n\t	vmulpd	0x60(%%rsi),%%ymm14,%%ymm14	\n\t"\
	"vmulpd	    (%%rsi),%%ymm0,%%ymm0	\n\t	vmulpd	    (%%rsi),%%ymm8 ,%%ymm8  \n\t"\
	"vmulpd	0xa0(%%rsi),%%ymm7,%%ymm7	\n\t	vmulpd	0xa0(%%rsi),%%ymm15,%%ymm15	\n\t"\
	"vmulpd	0xc0(%%rsi),%%ymm2,%%ymm2	\n\t	vmulpd	0xc0(%%rsi),%%ymm10,%%ymm10	\n\t"\
	"vaddpd	    (%%rcx),%%ymm0,%%ymm0	\n\t	vaddpd	0x20(%%rcx),%%ymm8 ,%%ymm8  \n\t"\
	"vaddpd	%%ymm4,%%ymm6,%%ymm6		\n\t	vaddpd	%%ymm12,%%ymm14,%%ymm14		\n\t"\
	"vsubpd	%%ymm2,%%ymm1,%%ymm1		\n\t	vsubpd	%%ymm10,%%ymm9 ,%%ymm9  	\n\t"\
	"vsubpd	%%ymm7,%%ymm4,%%ymm4		\n\t	vsubpd	%%ymm15,%%ymm12,%%ymm12		\n\t"\
	"vsubpd	%%ymm2,%%ymm3,%%ymm3		\n\t	vsubpd	%%ymm10,%%ymm11,%%ymm11		\n\t"\
		"movq	%[__o1],%%rax		\n\t"\
		"movq	%[__o2],%%rbx		\n\t"\
		"movq	%[__o3],%%rcx		\n\t"\
		"movq	%[__o4],%%rdx		\n\t"\
		"movq	%[__o5],%%rsi		\n\t"\
		"movq	%[__o6],%%rdi		\n\t"\
		"vmovaps	%%ymm0,%%ymm2		\n\t	vmovaps	%%ymm8 ,%%ymm10		\n\t"\
		"vmovaps	%%ymm5,%%ymm7		\n\t	vmovaps	%%ymm13,%%ymm15		\n\t"\
		"vaddpd	%%ymm1,%%ymm0,%%ymm0	\n\t	vaddpd	%%ymm9 ,%%ymm8 ,%%ymm8  	\n\t"\
		"vaddpd	%%ymm6,%%ymm5,%%ymm5	\n\t	vaddpd	%%ymm14,%%ymm13,%%ymm13		\n\t"\
		"vaddpd	%%ymm3,%%ymm1,%%ymm1	\n\t	vaddpd	%%ymm11,%%ymm9 ,%%ymm9  	\n\t"\
		"vaddpd	%%ymm4,%%ymm6,%%ymm6	\n\t	vaddpd	%%ymm12,%%ymm14,%%ymm14		\n\t"\
		"vaddpd	%%ymm2,%%ymm3,%%ymm3	\n\t	vaddpd	%%ymm10,%%ymm11,%%ymm11		\n\t"\
		"vaddpd	%%ymm7,%%ymm4,%%ymm4	\n\t	vaddpd	%%ymm15,%%ymm12,%%ymm12		\n\t"\
		"vsubpd	%%ymm1,%%ymm2,%%ymm2	\n\t	vsubpd	%%ymm9 ,%%ymm10,%%ymm10		\n\t"\
		"vsubpd	%%ymm6,%%ymm7,%%ymm7	\n\t	vsubpd	%%ymm14,%%ymm15,%%ymm15		\n\t"\
		/* ymm1,6,9,14 free ... Note the order reversal on the 3rd pair of outputs: */\
		"vsubpd	%%ymm13,%%ymm0 ,%%ymm0 	\n\t	vsubpd	%%ymm15,%%ymm2 ,%%ymm2 		\n\t	vsubpd	%%ymm12,%%ymm3 ,%%ymm3 		\n\t"\
		"vsubpd	%%ymm5 ,%%ymm8 ,%%ymm8  \n\t	vsubpd	%%ymm7 ,%%ymm10,%%ymm10		\n\t	vsubpd	%%ymm4 ,%%ymm11,%%ymm11		\n\t"\
		"vaddpd	%%ymm13,%%ymm13,%%ymm13	\n\t	vaddpd	%%ymm15,%%ymm15,%%ymm15		\n\t	vaddpd	%%ymm12,%%ymm12,%%ymm12		\n\t"\
		"vaddpd	%%ymm5 ,%%ymm5 ,%%ymm5 	\n\t	vaddpd	%%ymm7 ,%%ymm7 ,%%ymm7 		\n\t	vaddpd	%%ymm4 ,%%ymm4 ,%%ymm4 		\n\t"\
		"vaddpd	%%ymm0 ,%%ymm13,%%ymm13	\n\t	vaddpd	%%ymm2 ,%%ymm15,%%ymm15		\n\t	vaddpd	%%ymm3 ,%%ymm12,%%ymm12		\n\t"\
		"vaddpd	%%ymm8 ,%%ymm5 ,%%ymm5 	\n\t	vaddpd	%%ymm10,%%ymm7 ,%%ymm7 		\n\t	vaddpd	%%ymm11,%%ymm4 ,%%ymm4 		\n\t"\
		"vmovaps	%%ymm0 ,    (%%rax)	\n\t	vmovaps	%%ymm2 ,    (%%rbx)	\n\t	vmovaps	%%ymm3 ,    (%%rdx)	\n\t"/* B124r */\
		"vmovaps	%%ymm8 ,0x20(%%rdi)	\n\t	vmovaps	%%ymm10,0x20(%%rsi)	\n\t	vmovaps	%%ymm11,0x20(%%rcx)	\n\t"/* B653i */\
		"vmovaps	%%ymm13,    (%%rdi)	\n\t	vmovaps	%%ymm15,    (%%rsi)	\n\t	vmovaps	%%ymm12,    (%%rcx)	\n\t"/* B653r */\
		"vmovaps	%%ymm5 ,0x20(%%rax)	\n\t	vmovaps	%%ymm7 ,0x20(%%rbx)	\n\t	vmovaps	%%ymm4 ,0x20(%%rdx)	\n\t"/* B124i */\
		"\n\t"\
		:					/* outputs: none */\
		: [__i0] "m" (Xi0)	/* All inputs from memory addresses here */\
		 ,[__i1] "m" (Xi1)\
		 ,[__i2] "m" (Xi2)\
		 ,[__i3] "m" (Xi3)\
		 ,[__i4] "m" (Xi4)\
		 ,[__i5] "m" (Xi5)\
		 ,[__i6] "m" (Xi6)\
		 ,[__cc] "m" (Xcc)\
		 ,[__o0] "m" (Xo0)\
		 ,[__o1] "m" (Xo1)\
		 ,[__o2] "m" (Xo2)\
		 ,[__o3] "m" (Xo3)\
		 ,[__o4] "m" (Xo4)\
		 ,[__o5] "m" (Xo5)\
		 ,[__o6] "m" (Xo6)\
		: "cc","memory","rax","rbx","rcx","rdx","rsi","rdi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"		/* Clobbered registers */\
	);\
	}

	/* Twiddleless version of SSE2_RADIX8_DIF_TWIDDLE. Inputs enter in memory locations __r0 + [__i1,__i2,__i3,__i4,__i5,__i6,__i7],;
	where r0 is a memory address and the i's are LITERAL [BYTE] OFFSETS. Outputs go into memory locations __o0,__o1,__o2,__o3,__o4,__o5,__o6,__o7, assumed disjoint with inputs:\
	*/
	#define SSE2_RADIX8_DIF_0TWIDDLE(Xr0, Xi1,Xi2,Xi3,Xi4,Xi5,Xi6,Xi7, Xo0,Xo1,Xo2,Xo3,Xo4,Xo5,Xo6,Xo7, Xisrt2)\
	{\
	__asm__ volatile (\
		"/* 1st of 2 radix-4 subtransforms, data in ymm0-7: */\n\t	/* 2nd of 2 radix-4 subtransforms, data in ymm8-15: */\n\t"\
		"movq	%[__r0],%%rax	/* i0 = r00 */	\n\t		movq	%[__i1],%%r10		/* i1 */	\n\t"\
		"movq	%[__i2],%%rbx	/* i2 */		\n\t		movq	%[__i3],%%r11		/* i3 */	\n\t"\
		"movq	%[__i4],%%rcx	/* i4 */		\n\t		movq	%[__i5],%%r12		/* i5 */	\n\t"\
		"movq	%[__i6],%%rdx	/* i6 */		\n\t		movq	%[__i7],%%r13		/* i7 */	\n\t"\
		"addq	%%rax,%%rbx						\n\t		addq	%%rax,%%r10						\n\t"\
		"addq	%%rax,%%rcx						\n\t		addq	%%rax,%%r11						\n\t"\
		"addq	%%rax,%%rdx						\n\t		addq	%%rax,%%r12						\n\t"\
		"movq	%[__isrt2],%%rsi				\n\t		addq	%%rax,%%r13						\n\t"\
		"										\n\t		/* p1,5 combo: x+y into ymm8 /1, x-y in ymm10/3: */	\n\t"\
		"/* p0,4 combo: x+-y into ymm0/1, 2/3, resp: */\n\t	vmovaps	     (%%r12),%%ymm8 			\n\t"\
		"										\n\t		vmovaps	0x020(%%r12),%%ymm9 			\n\t"\
		"vmovaps	     (%%rcx),%%ymm0			\n\t		vmovaps	     (%%r10),%%ymm10			\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm1			\n\t		vmovaps	0x020(%%r10),%%ymm11			\n\t"\
		"vmovaps	     (%%rax),%%ymm2			\n\t		vsubpd	%%ymm8 ,%%ymm10,%%ymm10			\n\t"\
		"vmovaps	0x020(%%rax),%%ymm3			\n\t		vsubpd	%%ymm9 ,%%ymm11,%%ymm11			\n\t"\
		"vsubpd	%%ymm0,%%ymm2,%%ymm2			\n\t		vaddpd	%%ymm8 ,%%ymm8 ,%%ymm8 			\n\t"\
		"vsubpd	%%ymm1,%%ymm3,%%ymm3			\n\t		vaddpd	%%ymm9 ,%%ymm9 ,%%ymm9 			\n\t"\
		"vaddpd	%%ymm0,%%ymm0,%%ymm0			\n\t		vaddpd	%%ymm10,%%ymm8 ,%%ymm8 			\n\t"\
		"vaddpd	%%ymm1,%%ymm1,%%ymm1			\n\t		vaddpd	%%ymm11,%%ymm9 ,%%ymm9 			\n\t"\
		"vaddpd	%%ymm2,%%ymm0,%%ymm0			\n\t		/* p3,7 combo: x+y into ymm14/7, x-y in ymm12/5: */	\n\t"\
		"vaddpd	%%ymm3,%%ymm1,%%ymm1			\n\t		vmovaps	     (%%r11),%%ymm12			\n\t"\
		"										\n\t		vmovaps	0x020(%%r11),%%ymm13			\n\t"\
		"/* p2,6 combo: x+-y into ymm4/5, 6/7, resp: */\n\t	vmovaps	     (%%r13),%%ymm14			\n\t"\
		"										\n\t		vmovaps	0x020(%%r13),%%ymm15			\n\t"\
		"vmovaps	     (%%rdx),%%ymm4			\n\t		vsubpd	%%ymm14,%%ymm12,%%ymm12			\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm5			\n\t		vsubpd	%%ymm15,%%ymm13,%%ymm13			\n\t"\
		"vmovaps	     (%%rbx),%%ymm6			\n\t		vaddpd	%%ymm14,%%ymm14,%%ymm14			\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm7			\n\t		vaddpd	%%ymm15,%%ymm15,%%ymm15			\n\t"\
		"vsubpd	%%ymm4,%%ymm6,%%ymm6			\n\t		vaddpd	%%ymm12,%%ymm14,%%ymm14			\n\t"\
		"vsubpd	%%ymm5,%%ymm7,%%ymm7			\n\t		vaddpd	%%ymm13,%%ymm15,%%ymm15			\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4			\n\t		/* Finish radix-4 butterfly, tmp-store 1st of 4 outputs to free up 2 registers: */\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5			\n\t		vsubpd	%%ymm14,%%ymm8 ,%%ymm8 			\n\t"\
		"vaddpd	%%ymm6,%%ymm4,%%ymm4			\n\t		vsubpd	%%ymm15,%%ymm9 ,%%ymm9 			\n\t"\
		"vaddpd	%%ymm7,%%ymm5,%%ymm5			\n\t		vsubpd	%%ymm13,%%ymm10,%%ymm10			\n\t"\
		"										\n\t		vsubpd	%%ymm12,%%ymm11,%%ymm11			\n\t"\
		"vsubpd	%%ymm4,%%ymm0,%%ymm0			\n\t		vaddpd	%%ymm14,%%ymm14,%%ymm14			\n\t"\
		"vsubpd	%%ymm7,%%ymm2,%%ymm2			\n\t		vaddpd	%%ymm13,%%ymm13,%%ymm13			\n\t"\
		"vsubpd	%%ymm5,%%ymm1,%%ymm1			\n\t		vaddpd	%%ymm15,%%ymm15,%%ymm15			\n\t"\
		"vsubpd	%%ymm6,%%ymm3,%%ymm3			\n\t		vaddpd	%%ymm12,%%ymm12,%%ymm12			\n\t"\
		"													vaddpd	%%ymm8 ,%%ymm14,%%ymm14			\n\t"\
		"													vaddpd	%%ymm10,%%ymm13,%%ymm13			\n\t"\
		"													vaddpd	%%ymm9 ,%%ymm15,%%ymm15			\n\t"\
		"													vaddpd	%%ymm11,%%ymm12,%%ymm12			\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4			\n\t		vmovaps	%%ymm14,     (%%r10)			\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7			\n\t		vmovaps	%%ymm15,0x020(%%r10)			\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5			\n\t		vmovaps	%%ymm10,%%ymm14					\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6			\n\t		vmovaps	%%ymm13,%%ymm15					\n\t"\
		"vaddpd	%%ymm0,%%ymm4,%%ymm4			\n\t		vsubpd	%%ymm12,%%ymm10,%%ymm10			\n\t"\
		"vaddpd	%%ymm2,%%ymm7,%%ymm7			\n\t		vsubpd	%%ymm11,%%ymm13,%%ymm13			\n\t"\
		"vaddpd	%%ymm1,%%ymm5,%%ymm5			\n\t		vaddpd	%%ymm14,%%ymm12,%%ymm12			\n\t"\
		"vaddpd	%%ymm3,%%ymm6,%%ymm6			\n\t		vaddpd	%%ymm15,%%ymm11,%%ymm11			\n\t"\
		"													vmovaps	(%%rsi),%%ymm14	/* isrt2 */		\n\t"\
		"													vmulpd	%%ymm14,%%ymm10,%%ymm10			\n\t"\
		"													vmulpd	%%ymm14,%%ymm13,%%ymm13			\n\t"\
		"													vmulpd	%%ymm14,%%ymm12,%%ymm12			\n\t"\
		"										\n\t		vmulpd	%%ymm14,%%ymm11,%%ymm11			\n\t"\
		"vmovaps      (%%r10),%%ymm14 /* reload spill */\n\t"\
		"vmovaps 0x020(%%r10),%%ymm15 /* reload spill */\n\t"\
		"										\n\t"\
		"/* Inline of SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_A(r0): Combine radix-4 subtransforms and write outputs: */\n\t"\
		"/***** t0,1,2,3,4,5,6,7 in ymm[ 4, 5| 2,6| 0, 1| 7,3] *****/\n\t"\
		"/***** t8,9,a,b,c,d,e,f in ymm[14,15|10,12| 8, 9|13,11] */\n\t"\
		"movq	%[__o4],%%rax					\n\t		vsubpd   %%ymm10,%%ymm2	,%%ymm2			\n\t"\
		"movq	%[__o5],%%rbx					\n\t		vsubpd   %%ymm12,%%ymm6	,%%ymm6			\n\t"\
		"movq	%[__o6],%%rcx					\n\t		vaddpd   %%ymm10,%%ymm10,%%ymm10		\n\t"\
		"movq	%[__o7],%%rdx					\n\t		vaddpd   %%ymm12,%%ymm12,%%ymm12		\n\t"\
		"										\n\t		vaddpd   %%ymm2 ,%%ymm10,%%ymm10		\n\t"\
		"vsubpd   %%ymm11,%%ymm7 ,%%ymm7 		\n\t		vaddpd   %%ymm6 ,%%ymm12,%%ymm12		\n\t"\
		"vsubpd   %%ymm13,%%ymm3 ,%%ymm3 		\n\t"\
		"vaddpd   %%ymm11,%%ymm11,%%ymm11		\n\t		vmovaps	%%ymm2 ,     (%%rbx)	/* o5r */	\n\t"\
		"vaddpd   %%ymm13,%%ymm13,%%ymm13		\n\t		vmovaps	%%ymm6 ,0x020(%%rbx)	/* o5i */	\n\t"\
		"vaddpd   %%ymm7 ,%%ymm11,%%ymm11		\n\t		vmovaps	%%ymm10,     (%%rax)	/* o4r */	\n\t"\
		"vaddpd   %%ymm3 ,%%ymm13,%%ymm13		\n\t		vmovaps	%%ymm12,0x020(%%rax)	/* o4i */	\n\t"\
		"										\n\t"\
		"vmovaps	%%ymm7 ,     (%%rcx)	/* o6r */\n\t"\
		"vmovaps	%%ymm3 ,0x020(%%rdx)	/* o7i */\n\t"\
		"vmovaps	%%ymm11,     (%%rdx)	/* o7r */\n\t"\
		"vmovaps	%%ymm13,0x020(%%rcx)	/* o6i */\n\t"\
		"										\n\t"\
		"movq	%[__o0],%%rax					\n\t"\
		"movq	%[__o1],%%rbx					\n\t"\
		"movq	%[__o2],%%rcx					\n\t"\
		"movq	%[__o3],%%rdx					\n\t"\
		"										\n\t"\
		"vsubpd	%%ymm14,%%ymm4 ,%%ymm4 			\n\t"\
		"vsubpd	%%ymm15,%%ymm5 ,%%ymm5 			\n\t"\
		"vsubpd	%%ymm9 ,%%ymm0 ,%%ymm0 			\n\t"\
		"vsubpd	%%ymm8 ,%%ymm1 ,%%ymm1 			\n\t"\
		"vaddpd	%%ymm14,%%ymm14,%%ymm14			\n\t		vmovaps	%%ymm4 ,     (%%rbx)	/* o1r */	\n\t"\
		"vaddpd	%%ymm15,%%ymm15,%%ymm15			\n\t		vmovaps	%%ymm5 ,0x020(%%rbx)	/* o1i */	\n\t"\
		"vaddpd	%%ymm9 ,%%ymm9 ,%%ymm9 			\n\t		vmovaps	%%ymm0 ,     (%%rcx)	/* o2r */	\n\t"\
		"vaddpd	%%ymm8 ,%%ymm8 ,%%ymm8 			\n\t		vmovaps	%%ymm1 ,0x020(%%rdx)	/* o3i */	\n\t"\
		"vaddpd	%%ymm4 ,%%ymm14,%%ymm14			\n\t"\
		"vaddpd	%%ymm5 ,%%ymm15,%%ymm15			\n\t"\
		"vaddpd	%%ymm0 ,%%ymm9 ,%%ymm9 			\n\t"\
		"vaddpd	%%ymm1 ,%%ymm8 ,%%ymm8 			\n\t"\
		"										\n\t"\
		"vmovaps	%%ymm14,     (%%rax)	/* o0r */\n\t"\
		"vmovaps	%%ymm15,0x020(%%rax)	/* o0r */\n\t"\
		"vmovaps	%%ymm9 ,     (%%rdx)	/* o3r */\n\t"\
		"vmovaps	%%ymm8 ,0x020(%%rcx)	/* o2i */\n\t"\
		"										\n\t"\
		:					/* outputs: none */\
		: [__r0] "m" (Xr0)	/* All inputs from memory addresses here */\
		 ,[__i1] "e" (Xi1)\
		 ,[__i2] "e" (Xi2)\
		 ,[__i3] "e" (Xi3)\
		 ,[__i4] "e" (Xi4)\
		 ,[__i5] "e" (Xi5)\
		 ,[__i6] "e" (Xi6)\
		 ,[__i7] "e" (Xi7)\
		 ,[__o0] "m" (Xo0)\
		 ,[__o1] "m" (Xo1)\
		 ,[__o2] "m" (Xo2)\
		 ,[__o3] "m" (Xo3)\
		 ,[__o4] "m" (Xo4)\
		 ,[__o5] "m" (Xo5)\
		 ,[__o6] "m" (Xo6)\
		 ,[__o7] "m" (Xo7)\
		 ,[__isrt2] "m" (Xisrt2)\
		: "cc","memory","rax","rbx","rcx","rdx","rsi","r10","r11","r12","r13","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"		/* Clobbered registers */\
	);\
	}

	/* Twiddleless version of SSE2_RADIX8_DIT_TWIDDLE. Inputs enter in memory locations __i0,__i1,__i2,__i3,__i4,__i5,__i6,__i7.
	Outputs go into 16 contiguous 32-byte memory locations starting at __out and assumed disjoint with inputs.
	This macro built on the same code template as SSE2_RADIX8_DIF_TWIDDLE0, but with the I/O-location indices mutually bit reversed:
	01234567 <--> 04261537, which can be effected via the pairwise swaps 1 <--> 4 and 3 <--> 6.
	*/
	#define	SSE2_RADIX8_DIT_0TWIDDLE(Xi0,Xi1,Xi2,Xi3,Xi4,Xi5,Xi6,Xi7, Xout, Xisrt2)\
	{\
	__asm__ volatile (\
		"/* 1st of 2 radix-4 subtransforms, data in ymm0-7: */\n\t	/* 2nd of 2 radix-4 subtransforms, data in ymm8-15: */\n\t"\
		"movq	%[__i0],%%rax					\n\t		movq	%[__i4],%%r10					\n\t"\
		"movq	%[__i1],%%rbx					\n\t		movq	%[__i5],%%r11					\n\t"\
		"movq	%[__i2],%%rcx					\n\t		movq	%[__i6],%%r12					\n\t"\
		"movq	%[__i3],%%rdx					\n\t		movq	%[__i7],%%r13					\n\t"\
		"										\n\t		/* p1,5 combo: x+-y into ymm8/1, 10/3, resp: */	\n\t"\
		"/* p0,4 combo: x+-y into ymm0/1, 2/3, resp: */\n\tvmovaps	     (%%r11),%%ymm8 			\n\t"\
		"										\n\t		vmovaps	0x020(%%r11),%%ymm9 			\n\t"\
		"vmovaps	     (%%rbx),%%ymm0			\n\t		vmovaps	     (%%r10),%%ymm10			\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm1			\n\t		vmovaps	0x020(%%r10),%%ymm11			\n\t"\
		"vmovaps	     (%%rax),%%ymm2			\n\t		vsubpd	%%ymm8 ,%%ymm10,%%ymm10			\n\t"\
		"vmovaps	0x020(%%rax),%%ymm3			\n\t		vsubpd	%%ymm9 ,%%ymm11,%%ymm11			\n\t"\
		"vsubpd	%%ymm0,%%ymm2,%%ymm2			\n\t		vaddpd	%%ymm8 ,%%ymm8 ,%%ymm8 			\n\t"\
		"vsubpd	%%ymm1,%%ymm3,%%ymm3			\n\t		vaddpd	%%ymm9 ,%%ymm9 ,%%ymm9 			\n\t"\
		"vaddpd	%%ymm0,%%ymm0,%%ymm0			\n\t		vaddpd	%%ymm10,%%ymm8 ,%%ymm8 			\n\t"\
		"vaddpd	%%ymm1,%%ymm1,%%ymm1			\n\t		vaddpd	%%ymm11,%%ymm9 ,%%ymm9 			\n\t"\
		"vaddpd	%%ymm2,%%ymm0,%%ymm0			\n\t		/* p3,7 combo: x+-y into ymm14/7, 12/5, resp: */	\n\t"\
		"vaddpd	%%ymm3,%%ymm1,%%ymm1			\n\t		vmovaps	     (%%r12),%%ymm12			\n\t"\
		"										\n\t		vmovaps	0x020(%%r12),%%ymm13			\n\t"\
		"/* p2,6 combo: x+-y into ymm4/5, 6/7, resp: */\n\t	vmovaps	     (%%r13),%%ymm14			\n\t"\
		"										\n\t		vmovaps	0x020(%%r13),%%ymm15			\n\t"\
		"vmovaps	     (%%rdx),%%ymm4			\n\t		vsubpd	%%ymm14,%%ymm12,%%ymm12			\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm5			\n\t		vsubpd	%%ymm15,%%ymm13,%%ymm13			\n\t"\
		"vmovaps	     (%%rcx),%%ymm6			\n\t		vaddpd	%%ymm14,%%ymm14,%%ymm14			\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm7			\n\t		vaddpd	%%ymm15,%%ymm15,%%ymm15			\n\t"\
		"vsubpd	%%ymm4,%%ymm6,%%ymm6			\n\t		vaddpd	%%ymm12,%%ymm14,%%ymm14			\n\t"\
		"vsubpd	%%ymm5,%%ymm7,%%ymm7			\n\t		vaddpd	%%ymm13,%%ymm15,%%ymm15			\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4			\n\t		/* Finish radix-4 butterfly, tmp-store 1st of 4 outputs to free up 2 registers: */\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5			\n\t		vsubpd	%%ymm14,%%ymm8 ,%%ymm8 			\n\t"\
		"vaddpd	%%ymm6,%%ymm4,%%ymm4			\n\t		vsubpd	%%ymm15,%%ymm9 ,%%ymm9 			\n\t"\
		"vaddpd	%%ymm7,%%ymm5,%%ymm5			\n\t		vsubpd	%%ymm13,%%ymm10,%%ymm10			\n\t"\
		"										\n\t		vsubpd	%%ymm12,%%ymm11,%%ymm11			\n\t"\
		"vsubpd	%%ymm4,%%ymm0,%%ymm0			\n\t		vaddpd	%%ymm14,%%ymm14,%%ymm14			\n\t"\
		"vsubpd	%%ymm7,%%ymm2,%%ymm2			\n\t		vaddpd	%%ymm13,%%ymm13,%%ymm13			\n\t"\
		"vsubpd	%%ymm5,%%ymm1,%%ymm1			\n\t		vaddpd	%%ymm15,%%ymm15,%%ymm15			\n\t"\
		"vsubpd	%%ymm6,%%ymm3,%%ymm3			\n\t		vaddpd	%%ymm12,%%ymm12,%%ymm12			\n\t"\
		"													vaddpd	%%ymm8 ,%%ymm14,%%ymm14			\n\t"\
		"													vaddpd	%%ymm10,%%ymm13,%%ymm13			\n\t"\
		"													vaddpd	%%ymm9 ,%%ymm15,%%ymm15			\n\t"\
		"													vaddpd	%%ymm11,%%ymm12,%%ymm12			\n\t"\
		"													movq	%[__isrt2],%%rsi	/* isrt2 */	\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4			\n\t		vmovaps	%%ymm14,     (%%rax)	/* spill*/	\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7			\n\t		vmovaps	%%ymm15,0x020(%%rax)	/* spill*/	\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5			\n\t		vmovaps	%%ymm10,%%ymm14					\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6			\n\t		vmovaps	%%ymm13,%%ymm15					\n\t"\
		"vaddpd	%%ymm0,%%ymm4,%%ymm4			\n\t		vsubpd	%%ymm12,%%ymm10,%%ymm10			\n\t"\
		"vaddpd	%%ymm2,%%ymm7,%%ymm7			\n\t		vsubpd	%%ymm11,%%ymm13,%%ymm13			\n\t"\
		"vaddpd	%%ymm1,%%ymm5,%%ymm5			\n\t		vaddpd	%%ymm14,%%ymm12,%%ymm12			\n\t"\
		"vaddpd	%%ymm3,%%ymm6,%%ymm6			\n\t		vaddpd	%%ymm15,%%ymm11,%%ymm11			\n\t"\
		"													vmovaps	(%%rsi),%%ymm14		/* isrt2 */	\n\t"\
		"													vmulpd	%%ymm14,%%ymm10,%%ymm10			\n\t"\
		"													vmulpd	%%ymm14,%%ymm13,%%ymm13			\n\t"\
		"													vmulpd	%%ymm14,%%ymm12,%%ymm12			\n\t"\
		"													vmulpd	%%ymm14,%%ymm11,%%ymm11			\n\t"\
		"/* Combine radix-4 subtransforms and write outputs: */\n\t"\
		"\n\t"\
		"vmovaps      (%%rax),%%ymm14/* reload spill */\n\t	vsubpd   %%ymm10,%%ymm2	,%%ymm2			\n\t"\
		"vmovaps 0x020(%%rax),%%ymm15/* reload spill */\n\t	vsubpd   %%ymm12,%%ymm6	,%%ymm6			\n\t"\
		"													vaddpd   %%ymm10,%%ymm10,%%ymm10		\n\t"\
		"movq	%[__out],%%rax					\n\t		vaddpd   %%ymm12,%%ymm12,%%ymm12		\n\t"\
		"										\n\t		vaddpd   %%ymm2 ,%%ymm10,%%ymm10		\n\t"\
		"vsubpd   %%ymm11,%%ymm7 ,%%ymm7 		\n\t		vaddpd   %%ymm6 ,%%ymm12,%%ymm12		\n\t"\
		"vsubpd   %%ymm13,%%ymm3 ,%%ymm3 		\n\t"\
		"vaddpd   %%ymm11,%%ymm11,%%ymm11		\n\t		vmovaps	%%ymm2 ,0x140(%%rax)	/* o5r */	\n\t"\
		"vaddpd   %%ymm13,%%ymm13,%%ymm13		\n\t		vmovaps	%%ymm6 ,0x160(%%rax)	/* o5i */	\n\t"\
		"vaddpd   %%ymm7 ,%%ymm11,%%ymm11		\n\t		vmovaps	%%ymm10,0x040(%%rax)	/* o1r */	\n\t"\
		"vaddpd   %%ymm3 ,%%ymm13,%%ymm13		\n\t		vmovaps	%%ymm12,0x060(%%rax)	/* o1i */	\n\t"\
		"										\n\t"\
		"vmovaps	%%ymm7 ,0x0c0(%%rax)	/* o3r */\n\t"\
		"vmovaps	%%ymm3 ,0x1e0(%%rax)	/* o7i */\n\t"\
		"vmovaps	%%ymm11,0x1c0(%%rax)	/* o7r */\n\t"\
		"vmovaps	%%ymm13,0x0e0(%%rax)	/* o3i */\n\t"\
		"										\n\t"\
		"vsubpd	%%ymm14,%%ymm4 ,%%ymm4 			\n\t"\
		"vsubpd	%%ymm15,%%ymm5 ,%%ymm5 			\n\t"\
		"vsubpd	%%ymm9 ,%%ymm0 ,%%ymm0 			\n\t"\
		"vsubpd	%%ymm8 ,%%ymm1 ,%%ymm1 			\n\t"\
		"vaddpd	%%ymm14,%%ymm14,%%ymm14			\n\t		vmovaps	%%ymm4 ,0x100(%%rax)	/* o4r */	\n\t"\
		"vaddpd	%%ymm15,%%ymm15,%%ymm15			\n\t		vmovaps	%%ymm5 ,0x120(%%rax)	/* o4i */	\n\t"\
		"vaddpd	%%ymm9 ,%%ymm9 ,%%ymm9 			\n\t		vmovaps	%%ymm0 ,0x080(%%rax)	/* o2r */	\n\t"\
		"vaddpd	%%ymm8 ,%%ymm8 ,%%ymm8 			\n\t		vmovaps	%%ymm1 ,0x1a0(%%rax)	/* o6i */	\n\t"\
		"vaddpd	%%ymm4 ,%%ymm14,%%ymm14			\n\t"\
		"vaddpd	%%ymm5 ,%%ymm15,%%ymm15			\n\t"\
		"vaddpd	%%ymm0 ,%%ymm9 ,%%ymm9 			\n\t"\
		"vaddpd	%%ymm1 ,%%ymm8 ,%%ymm8 			\n\t"\
		"										\n\t"\
		"vmovaps	%%ymm14,     (%%rax)	/* o0r */\n\t"\
		"vmovaps	%%ymm15,0x020(%%rax)	/* o0r */\n\t"\
		"vmovaps	%%ymm9 ,0x180(%%rax)	/* o6r */\n\t"\
		"vmovaps	%%ymm8 ,0x0a0(%%rax)	/* o2i */\n\t"\
		:					/* outputs: none */\
		: [__i0] "m" (Xi0)	/* All iputs from memory addresses here */\
		 ,[__i1] "m" (Xi1)\
		 ,[__i2] "m" (Xi2)\
		 ,[__i3] "m" (Xi3)\
		 ,[__i4] "m" (Xi4)\
		 ,[__i5] "m" (Xi5)\
		 ,[__i6] "m" (Xi6)\
		 ,[__i7] "m" (Xi7)\
		 ,[__out] "m" (Xout)\
		 ,[__isrt2] "m" (Xisrt2)\
		: "cc","memory","rax","rbx","rcx","rdx","rsi","r10","r11","r12","r13","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"		/* Clobbered registers */\
	);\
	}

	// Same as SSE2_RADIX8_DIT_0TWIDDLE but with user-specifiable [i.e. not nec. contiguous] output addresses:
	#define	SSE2_RADIX8_DIT_0TWIDDLE_OOP(Xi0,Xi1,Xi2,Xi3,Xi4,Xi5,Xi6,Xi7, Xo0,Xo1,Xo2,Xo3,Xo4,Xo5,Xo6,Xo7, Xisrt2)\
	{\
	__asm__ volatile (\
		"/* 1st of 2 radix-4 subtransforms, data in ymm0-7: */\n\t	/* 2nd of 2 radix-4 subtransforms, data in ymm8-15: */\n\t"\
		"movq	%[__i0],%%rax					\n\t		movq	%[__i4],%%r10					\n\t"\
		"movq	%[__i1],%%rbx					\n\t		movq	%[__i5],%%r11					\n\t"\
		"movq	%[__i2],%%rcx					\n\t		movq	%[__i6],%%r12					\n\t"\
		"movq	%[__i3],%%rdx					\n\t		movq	%[__i7],%%r13					\n\t"\
		"										\n\t		/* p1,5 combo: x+-y into ymm8/1, 10/3, resp: */	\n\t"\
		"/* p0,4 combo: x+-y into ymm0/1, 2/3, resp: */\n\tvmovaps	     (%%r11),%%ymm8 			\n\t"\
		"										\n\t		vmovaps	0x020(%%r11),%%ymm9 			\n\t"\
		"vmovaps	     (%%rbx),%%ymm0			\n\t		vmovaps	     (%%r10),%%ymm10			\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm1			\n\t		vmovaps	0x020(%%r10),%%ymm11			\n\t"\
		"vmovaps	     (%%rax),%%ymm2			\n\t		vsubpd	%%ymm8 ,%%ymm10,%%ymm10			\n\t"\
		"vmovaps	0x020(%%rax),%%ymm3			\n\t		vsubpd	%%ymm9 ,%%ymm11,%%ymm11			\n\t"\
		"vsubpd	%%ymm0,%%ymm2,%%ymm2			\n\t		vaddpd	%%ymm8 ,%%ymm8 ,%%ymm8 			\n\t"\
		"vsubpd	%%ymm1,%%ymm3,%%ymm3			\n\t		vaddpd	%%ymm9 ,%%ymm9 ,%%ymm9 			\n\t"\
		"vaddpd	%%ymm0,%%ymm0,%%ymm0			\n\t		vaddpd	%%ymm10,%%ymm8 ,%%ymm8 			\n\t"\
		"vaddpd	%%ymm1,%%ymm1,%%ymm1			\n\t		vaddpd	%%ymm11,%%ymm9 ,%%ymm9 			\n\t"\
		"vaddpd	%%ymm2,%%ymm0,%%ymm0			\n\t		/* p3,7 combo: x+-y into ymm14/7, 12/5, resp: */	\n\t"\
		"vaddpd	%%ymm3,%%ymm1,%%ymm1			\n\t		vmovaps	     (%%r12),%%ymm12			\n\t"\
		"										\n\t		vmovaps	0x020(%%r12),%%ymm13			\n\t"\
		"/* p2,6 combo: x+-y into ymm4/5, 6/7, resp: */\n\t	vmovaps	     (%%r13),%%ymm14			\n\t"\
		"										\n\t		vmovaps	0x020(%%r13),%%ymm15			\n\t"\
		"vmovaps	     (%%rdx),%%ymm4			\n\t		vsubpd	%%ymm14,%%ymm12,%%ymm12			\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm5			\n\t		vsubpd	%%ymm15,%%ymm13,%%ymm13			\n\t"\
		"vmovaps	     (%%rcx),%%ymm6			\n\t		vaddpd	%%ymm14,%%ymm14,%%ymm14			\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm7			\n\t		vaddpd	%%ymm15,%%ymm15,%%ymm15			\n\t"\
		"vsubpd	%%ymm4,%%ymm6,%%ymm6			\n\t		vaddpd	%%ymm12,%%ymm14,%%ymm14			\n\t"\
		"vsubpd	%%ymm5,%%ymm7,%%ymm7			\n\t		vaddpd	%%ymm13,%%ymm15,%%ymm15			\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4			\n\t		/* Finish radix-4 butterfly, tmp-store 1st of 4 outputs to free up 2 registers: */\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5			\n\t		vsubpd	%%ymm14,%%ymm8 ,%%ymm8 			\n\t"\
		"vaddpd	%%ymm6,%%ymm4,%%ymm4			\n\t		vsubpd	%%ymm15,%%ymm9 ,%%ymm9 			\n\t"\
		"vaddpd	%%ymm7,%%ymm5,%%ymm5			\n\t		vsubpd	%%ymm13,%%ymm10,%%ymm10			\n\t"\
		"										\n\t		vsubpd	%%ymm12,%%ymm11,%%ymm11			\n\t"\
		"vsubpd	%%ymm4,%%ymm0,%%ymm0			\n\t		vaddpd	%%ymm14,%%ymm14,%%ymm14			\n\t"\
		"vsubpd	%%ymm7,%%ymm2,%%ymm2			\n\t		vaddpd	%%ymm13,%%ymm13,%%ymm13			\n\t"\
		"vsubpd	%%ymm5,%%ymm1,%%ymm1			\n\t		vaddpd	%%ymm15,%%ymm15,%%ymm15			\n\t"\
		"vsubpd	%%ymm6,%%ymm3,%%ymm3			\n\t		vaddpd	%%ymm12,%%ymm12,%%ymm12			\n\t"\
		"													vaddpd	%%ymm8 ,%%ymm14,%%ymm14			\n\t"\
		"													vaddpd	%%ymm10,%%ymm13,%%ymm13			\n\t"\
		"													vaddpd	%%ymm9 ,%%ymm15,%%ymm15			\n\t"\
		"													vaddpd	%%ymm11,%%ymm12,%%ymm12			\n\t"\
		"													movq	%[__isrt2],%%rsi	/* isrt2 */	\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4			\n\t		vmovaps	%%ymm14,     (%%rax)	/* spill*/	\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7			\n\t		vmovaps	%%ymm15,0x020(%%rax)	/* spill*/	\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5			\n\t		vmovaps	%%ymm10,%%ymm14					\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6			\n\t		vmovaps	%%ymm13,%%ymm15					\n\t"\
		"vaddpd	%%ymm0,%%ymm4,%%ymm4			\n\t		vsubpd	%%ymm12,%%ymm10,%%ymm10			\n\t"\
		"vaddpd	%%ymm2,%%ymm7,%%ymm7			\n\t		vsubpd	%%ymm11,%%ymm13,%%ymm13			\n\t"\
		"vaddpd	%%ymm1,%%ymm5,%%ymm5			\n\t		vaddpd	%%ymm14,%%ymm12,%%ymm12			\n\t"\
		"vaddpd	%%ymm3,%%ymm6,%%ymm6			\n\t		vaddpd	%%ymm15,%%ymm11,%%ymm11			\n\t"\
		"													vmovaps	(%%rsi),%%ymm14		/* isrt2 */	\n\t"\
		"													vmulpd	%%ymm14,%%ymm10,%%ymm10			\n\t"\
		"													vmulpd	%%ymm14,%%ymm13,%%ymm13			\n\t"\
		"													vmulpd	%%ymm14,%%ymm12,%%ymm12			\n\t"\
		"													vmulpd	%%ymm14,%%ymm11,%%ymm11			\n\t"\
		"/* Combine radix-4 subtransforms and write outputs: */\n\t"\
		"\n\t"\
		"vmovaps      (%%rax),%%ymm14/* reload spill */\n\t	vsubpd   %%ymm10,%%ymm2	,%%ymm2			\n\t"\
		"vmovaps 0x020(%%rax),%%ymm15/* reload spill */\n\t	vsubpd   %%ymm12,%%ymm6	,%%ymm6			\n\t"\
		"movq	%[__o1],%%rax					\n\t		movq	%[__o5],%%rcx					\n\t"\
		"													vaddpd   %%ymm10,%%ymm10,%%ymm10		\n\t"\
		"													vaddpd   %%ymm12,%%ymm12,%%ymm12		\n\t"\
		"										\n\t		vaddpd   %%ymm2 ,%%ymm10,%%ymm10		\n\t"\
		"vsubpd   %%ymm11,%%ymm7 ,%%ymm7 		\n\t		vaddpd   %%ymm6 ,%%ymm12,%%ymm12		\n\t"\
		"vsubpd   %%ymm13,%%ymm3 ,%%ymm3 		\n\t"\
		"movq	%[__o3],%%rbx					\n\t		movq	%[__o7],%%rdx					\n\t"\
		"vaddpd   %%ymm11,%%ymm11,%%ymm11		\n\t		vmovaps	%%ymm2 ,    (%%rcx)	/* o5r */	\n\t"\
		"vaddpd   %%ymm13,%%ymm13,%%ymm13		\n\t		vmovaps	%%ymm6 ,0x20(%%rcx)	/* o5i */	\n\t"\
		"vaddpd   %%ymm7 ,%%ymm11,%%ymm11		\n\t		vmovaps	%%ymm10,    (%%rax)	/* o1r */	\n\t"\
		"vaddpd   %%ymm3 ,%%ymm13,%%ymm13		\n\t		vmovaps	%%ymm12,0x20(%%rax)	/* o1i */	\n\t"\
		"movq	%[__o0],%%rax					\n\t		movq	%[__o4],%%rcx					\n\t"\
		"										\n\t"\
		"vmovaps	%%ymm7 ,    (%%rbx)	/* o3r */	\n\t"\
		"vmovaps	%%ymm3 ,0x20(%%rdx)	/* o7i */	\n\t"\
		"vmovaps	%%ymm11,    (%%rdx)	/* o7r */	\n\t"\
		"vmovaps	%%ymm13,0x20(%%rbx)	/* o3i */	\n\t"\
		"										\n\t"\
		"movq	%[__o2],%%rbx					\n\t		movq	%[__o6],%%rdx					\n\t"\
		"vsubpd	%%ymm14,%%ymm4 ,%%ymm4 			\n\t"\
		"vsubpd	%%ymm15,%%ymm5 ,%%ymm5 			\n\t"\
		"vsubpd	%%ymm9 ,%%ymm0 ,%%ymm0 			\n\t"\
		"vsubpd	%%ymm8 ,%%ymm1 ,%%ymm1 			\n\t"\
		"vaddpd	%%ymm14,%%ymm14,%%ymm14			\n\t		vmovaps	%%ymm4 ,    (%%rcx)	/* o4r */	\n\t"\
		"vaddpd	%%ymm15,%%ymm15,%%ymm15			\n\t		vmovaps	%%ymm5 ,0x20(%%rcx)	/* o4i */	\n\t"\
		"vaddpd	%%ymm9 ,%%ymm9 ,%%ymm9 			\n\t		vmovaps	%%ymm0 ,    (%%rbx)	/* o2r */	\n\t"\
		"vaddpd	%%ymm8 ,%%ymm8 ,%%ymm8 			\n\t		vmovaps	%%ymm1 ,0x20(%%rdx)	/* o6i */	\n\t"\
		"vaddpd	%%ymm4 ,%%ymm14,%%ymm14			\n\t"\
		"vaddpd	%%ymm5 ,%%ymm15,%%ymm15			\n\t"\
		"vaddpd	%%ymm0 ,%%ymm9 ,%%ymm9 			\n\t"\
		"vaddpd	%%ymm1 ,%%ymm8 ,%%ymm8 			\n\t"\
		"										\n\t"\
		"vmovaps	%%ymm14,    (%%rax)	/* o0r */	\n\t"\
		"vmovaps	%%ymm15,0x20(%%rax)	/* o0i */	\n\t"\
		"vmovaps	%%ymm9 ,    (%%rdx)	/* o6r */	\n\t"\
		"vmovaps	%%ymm8 ,0x20(%%rbx)	/* o2i */	\n\t"\
		:					/* outputs: none */\
		: [__i0] "m" (Xi0)	/* All iputs from memory addresses here */\
		 ,[__i1] "m" (Xi1)\
		 ,[__i2] "m" (Xi2)\
		 ,[__i3] "m" (Xi3)\
		 ,[__i4] "m" (Xi4)\
		 ,[__i5] "m" (Xi5)\
		 ,[__i6] "m" (Xi6)\
		 ,[__i7] "m" (Xi7)\
		 ,[__o0] "m" (Xo0)\
		 ,[__o1] "m" (Xo1)\
		 ,[__o2] "m" (Xo2)\
		 ,[__o3] "m" (Xo3)\
		 ,[__o4] "m" (Xo4)\
		 ,[__o5] "m" (Xo5)\
		 ,[__o6] "m" (Xo6)\
		 ,[__o7] "m" (Xo7)\
		 ,[__isrt2] "m" (Xisrt2)\
		: "cc","memory","rax","rbx","rcx","rdx","rsi","r10","r11","r12","r13","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"		/* Clobbered registers */\
	);\
	}

	// AVX analog of dft_macro.h::RADIX_08_DIF_TWIDDLE_OOP - Result of adding separate I/O addressing to
	// radix8_dif_dit_pass_gcc64.h::SSE2_RADIX8_DIF_TWIDDLE:
	#define SSE2_RADIX8_DIF_TWIDDLE_OOP(Xi0,Xi1,Xi2,Xi3,Xi4,Xi5,Xi6,Xi7, Xo0,Xo1,Xo2,Xo3,Xo4,Xo5,Xo6,Xo7 ,Xc1,Xs1,Xc2,Xs2,Xc3,Xs3,Xc4,Xs4,Xc5,Xs5,Xc6,Xs6,Xc7,Xs7)\
	{\
	__asm__ volatile (\
		"													movq		%[i1]	,%%r10				\n\t"\
		"													movq		%[i5]	,%%r11				\n\t"\
		"											movq	%[c1],%%r12	\n\t	movq	%[c5],%%r14	\n\t"\
		"											movq	%[s1],%%r13	\n\t	movq	%[s5],%%r15	\n\t"\
		"													vmovaps	    (%%r10)	,%%ymm8 			\n\t"\
		"movq		%[i0]	,%%rax				\n\t		vmovaps	0x20(%%r10)	,%%ymm10			\n\t"\
		"movq		%[i4]	,%%rbx				\n\t		vmovaps	    (%%r10)	,%%ymm9 			\n\t"\
		"movq		%[c4]	,%%rcx				\n\t		vmovaps	0x20(%%r10)	,%%ymm11			\n\t"\
		"movq		%[s4]	,%%rsi			\n\t"\
	/* [rsi] (and if needed rdi) points to sine components of each sincos pair, which is not really a pair here in terms of relative addressing: */\
		"vmovaps	    (%%rax)	,%%ymm0			\n\t		vmulpd	    (%%r12)	,%%ymm8 ,%%ymm8 	\n\t"\
		"vmovaps	0x20(%%rax)	,%%ymm1			\n\t		vmulpd	    (%%r13)	,%%ymm10,%%ymm10	\n\t"\
		"vmovaps	    (%%rax)	,%%ymm6			\n\t		vmulpd	    (%%r13)	,%%ymm9 ,%%ymm9 	\n\t"\
		"vmovaps	0x20(%%rax)	,%%ymm7			\n\t		vmulpd	    (%%r12)	,%%ymm11,%%ymm11	\n\t"\
		"vmovaps	    (%%rbx)	,%%ymm2			\n\t		vsubpd	%%ymm10		,%%ymm8 ,%%ymm8 	\n\t"\
		"vmovaps	0x20(%%rbx)	,%%ymm3			\n\t		vaddpd	%%ymm11		,%%ymm9 ,%%ymm9 	\n\t"\
		"vmovaps	    (%%rbx)	,%%ymm4			\n\t		vmovaps	    (%%r11)	,%%ymm10			\n\t"\
		"vmovaps	0x20(%%rbx)	,%%ymm5			\n\t		vmovaps	0x20(%%r11)	,%%ymm11			\n\t"\
		"vmulpd	    (%%rcx)	,%%ymm2,%%ymm2		\n\t		vmovaps	    (%%r11)	,%%ymm12			\n\t"\
		"vmulpd	    (%%rcx)	,%%ymm3,%%ymm3		\n\t		vmovaps	0x20(%%r11)	,%%ymm13			\n\t"\
		"vmulpd	    (%%rsi)	,%%ymm4,%%ymm4		\n\t		vmulpd	    (%%r14)	,%%ymm10,%%ymm10	\n\t"\
		"vmulpd	    (%%rsi)	,%%ymm5,%%ymm5		\n\t		vmulpd	    (%%r15)	,%%ymm11,%%ymm11	\n\t"\
		"vsubpd	%%ymm5		,%%ymm2,%%ymm2		\n\t		vmulpd	    (%%r15)	,%%ymm12,%%ymm12	\n\t"\
		"vaddpd	%%ymm4		,%%ymm3,%%ymm3		\n\t		vmulpd	    (%%r14)	,%%ymm13,%%ymm13	\n\t"\
		"vaddpd	%%ymm2		,%%ymm0,%%ymm0		\n\t		vsubpd	%%ymm11		,%%ymm10,%%ymm10	\n\t"\
		"vaddpd	%%ymm3		,%%ymm1,%%ymm1		\n\t		vaddpd	%%ymm13		,%%ymm12,%%ymm12	\n\t"\
		"vsubpd	%%ymm2		,%%ymm6,%%ymm6		\n\t		vmovaps	%%ymm10		,%%ymm11			\n\t"\
		"vsubpd	%%ymm3		,%%ymm7,%%ymm7		\n\t		vmovaps	%%ymm12		,%%ymm13			\n\t"\
		"vmovaps	%%ymm0		,    (%%rax)	\n\t		vaddpd	%%ymm8 		,%%ymm10,%%ymm10	\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rax)	\n\t		vsubpd	%%ymm11		,%%ymm8 ,%%ymm8 	\n\t"\
		"vmovaps	%%ymm6		,    (%%rbx)	\n\t		vaddpd	%%ymm9 		,%%ymm12,%%ymm12	\n\t"\
		"vmovaps	%%ymm7		,0x20(%%rbx)	\n\t		vsubpd	%%ymm13		,%%ymm9 ,%%ymm9 	\n\t"\
		"movq		%[i2]	,%%rax				\n\t		vmovaps	%%ymm10		,    (%%r10)	\n\t"\
		"movq		%[i6]	,%%rbx				\n\t		vmovaps	%%ymm12		,0x20(%%r10)	\n\t"\
		"movq %[c2],%%rcx \n\t movq %[s2],%%rsi	\n\t		vmovaps	%%ymm8 		,    (%%r11)	\n\t"\
		"movq %[c6],%%rdx \n\t movq %[s6],%%rdi	\n\t		vmovaps	%%ymm9 		,0x20(%%r11)	\n\t"\
		"vmovaps	    (%%rax)	,%%ymm0			\n\t		movq		%[i3]	,%%r10				\n\t"\
		"vmovaps	0x20(%%rax)	,%%ymm2			\n\t		movq		%[i7]	,%%r11				\n\t"\
		"vmovaps	    (%%rax)	,%%ymm1			\n\t		movq %[c3],%%r12 \n\t movq %[s3],%%r13	\n\t"\
		"vmovaps	0x20(%%rax)	,%%ymm3			\n\t		movq %[c7],%%r14 \n\t movq %[s7],%%r15	\n\t"\
		"vmulpd	    (%%rcx)	,%%ymm0,%%ymm0		\n\t		vmovaps	    (%%r10)	,%%ymm8 			\n\t"\
		"vmulpd	    (%%rsi)	,%%ymm2,%%ymm2		\n\t		vmovaps	0x20(%%r10)	,%%ymm10			\n\t"\
		"vmulpd	    (%%rsi)	,%%ymm1,%%ymm1		\n\t		vmovaps	    (%%r10)	,%%ymm9 			\n\t"\
		"vmulpd	    (%%rcx)	,%%ymm3,%%ymm3		\n\t		vmovaps	0x20(%%r10)	,%%ymm11			\n\t"\
		"vsubpd	%%ymm2		,%%ymm0,%%ymm0		\n\t		vmulpd	    (%%r12)	,%%ymm8 ,%%ymm8 	\n\t"\
		"vaddpd	%%ymm3		,%%ymm1,%%ymm1		\n\t		vmulpd	    (%%r13)	,%%ymm10,%%ymm10	\n\t"\
		"vmovaps	    (%%rbx)	,%%ymm2			\n\t		vmulpd	    (%%r13)	,%%ymm9 ,%%ymm9 	\n\t"\
		"vmovaps	0x20(%%rbx)	,%%ymm3			\n\t		vmulpd	    (%%r12)	,%%ymm11,%%ymm11	\n\t"\
		"vmovaps	    (%%rbx)	,%%ymm4			\n\t		vsubpd	%%ymm10		,%%ymm8 ,%%ymm8 	\n\t"\
		"vmovaps	0x20(%%rbx)	,%%ymm5			\n\t		vaddpd	%%ymm11		,%%ymm9 ,%%ymm9 	\n\t"\
		"vmulpd	    (%%rdx)	,%%ymm2,%%ymm2		\n\t		vmovaps	    (%%r11)	,%%ymm10			\n\t"\
		"vmulpd	    (%%rdi)	,%%ymm3,%%ymm3		\n\t		vmovaps	0x20(%%r11)	,%%ymm11			\n\t"\
		"vmulpd	    (%%rdi)	,%%ymm4,%%ymm4		\n\t		vmovaps	    (%%r11)	,%%ymm12			\n\t"\
		"vmulpd	    (%%rdx)	,%%ymm5,%%ymm5		\n\t		vmovaps	0x20(%%r11)	,%%ymm13			\n\t"\
		"vsubpd	%%ymm3		,%%ymm2,%%ymm2		\n\t		vmulpd	    (%%r14)	,%%ymm10,%%ymm10	\n\t"\
		"vaddpd	%%ymm5		,%%ymm4,%%ymm4		\n\t		vmulpd	    (%%r15)	,%%ymm11,%%ymm11	\n\t"\
		"vmovaps	%%ymm2		,%%ymm3			\n\t		vmulpd	    (%%r15)	,%%ymm12,%%ymm12	\n\t"\
		"vmovaps	%%ymm4		,%%ymm5			\n\t		vmulpd	    (%%r14)	,%%ymm13,%%ymm13	\n\t"\
		"vaddpd	%%ymm0		,%%ymm2,%%ymm2		\n\t		vsubpd	%%ymm11		,%%ymm10,%%ymm10	\n\t"\
		"vsubpd	%%ymm3		,%%ymm0,%%ymm0		\n\t		vaddpd	%%ymm13		,%%ymm12,%%ymm12	\n\t"\
		"vaddpd	%%ymm1		,%%ymm4,%%ymm4		\n\t		vmovaps	%%ymm10		,%%ymm11			\n\t"\
		"vsubpd	%%ymm5		,%%ymm1,%%ymm1		\n\t		vmovaps	%%ymm12		,%%ymm13			\n\t"\
		"vmovaps	%%ymm2		,    (%%rax)	\n\t		vaddpd	%%ymm8 		,%%ymm10,%%ymm10	\n\t"\
		"vmovaps	%%ymm4		,0x20(%%rax)	\n\t		vsubpd	%%ymm11		,%%ymm8 ,%%ymm8 	\n\t"\
		"vmovaps	%%ymm0		,    (%%rbx)	\n\t		vaddpd	%%ymm9 		,%%ymm12,%%ymm12	\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rbx)	\n\t		vsubpd	%%ymm13		,%%ymm9 ,%%ymm9 	\n\t"\
		"													vmovaps	%%ymm10		,    (%%r10)	\n\t"\
		"													vmovaps	%%ymm12		,0x20(%%r10)	\n\t"\
		"													vmovaps	%%ymm8 		,    (%%r11)	\n\t"\
		"													vmovaps	%%ymm9 		,0x20(%%r11)	\n\t"\
/* combine to get 2 length-4 output subtransforms... */\
		"movq		%[i0]	,%%rax				\n\t		movq		%[i4]	,%%r10				\n\t"\
		"movq		%[i2]	,%%rbx				\n\t		movq		%[i6]	,%%r11				\n\t"\
		"vmovaps	    (%%rax)	,%%ymm0			\n\t		vmovaps	    (%%r10)	,%%ymm8 			\n\t"\
		"vmovaps	0x20(%%rax)	,%%ymm1			\n\t		vmovaps	0x20(%%r10)	,%%ymm9 			\n\t"\
		"vmovaps	%%ymm0		,%%ymm4			\n\t		vmovaps	%%ymm8 		,%%ymm12			\n\t"\
		"vmovaps	%%ymm1		,%%ymm5			\n\t		vmovaps	%%ymm9 		,%%ymm13			\n\t"\
		"vaddpd	    (%%rbx)	,%%ymm0,%%ymm0		\n\t		vsubpd	0x20(%%r11)	,%%ymm8 ,%%ymm8 	\n\t"\
		"vsubpd	    (%%rbx)	,%%ymm4,%%ymm4		\n\t		vaddpd	0x20(%%r11)	,%%ymm12,%%ymm12	\n\t"\
		"vaddpd	0x20(%%rbx)	,%%ymm1,%%ymm1		\n\t		vaddpd	    (%%r11)	,%%ymm9 ,%%ymm9 	\n\t"\
		"vsubpd	0x20(%%rbx)	,%%ymm5,%%ymm5		\n\t		vsubpd	    (%%r11)	,%%ymm13,%%ymm13	\n\t"\
		"movq		%[o0]	,%%rax				\n\t		movq		%[o4]	,%%r10				\n\t"\
		"movq		%[o2]	,%%rbx				\n\t		movq		%[o6]	,%%r11				\n\t"\
		"vmovaps	%%ymm0		,    (%%rax)	\n\t		vmovaps	%%ymm8 		,    (%%r10)	\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rax)	\n\t		vmovaps	%%ymm9 		,0x20(%%r10)	\n\t"\
		"vmovaps	%%ymm4		,    (%%rbx)	\n\t		vmovaps	%%ymm12		,    (%%r11)	\n\t"\
		"vmovaps	%%ymm5		,0x20(%%rbx)	\n\t		vmovaps	%%ymm13		,0x20(%%r11)	\n\t"\
		"movq		%[i1]	,%%rcx				\n\t		movq		%[i5]	,%%r12				\n\t"\
		"movq		%[i3]	,%%rdx				\n\t		movq		%[i7]	,%%r13				\n\t"\
		"vmovaps	    (%%rcx)	,%%ymm2			\n\t		vmovaps	    (%%r12)	,%%ymm10			\n\t"\
		"vmovaps	0x20(%%rcx)	,%%ymm3			\n\t		vmovaps	0x20(%%r12)	,%%ymm11			\n\t"\
		"vmovaps	%%ymm2		,%%ymm6			\n\t		vmovaps	%%ymm10		,%%ymm14			\n\t"\
		"vmovaps	%%ymm3		,%%ymm7			\n\t		vmovaps	%%ymm11		,%%ymm15			\n\t"\
		"vaddpd	    (%%rdx)	,%%ymm2,%%ymm2		\n\t		vsubpd	0x20(%%r13)	,%%ymm10,%%ymm10	\n\t"\
		"vsubpd	    (%%rdx)	,%%ymm6,%%ymm6		\n\t		vaddpd	0x20(%%r13)	,%%ymm14,%%ymm14	\n\t"\
		"vaddpd	0x20(%%rdx)	,%%ymm3,%%ymm3		\n\t		vaddpd	    (%%r13)	,%%ymm11,%%ymm11	\n\t"\
		"vsubpd	0x20(%%rdx)	,%%ymm7,%%ymm7		\n\t		vsubpd	    (%%r13)	,%%ymm15,%%ymm15	\n\t"\
		"movq		%[o1]	,%%rcx				\n\t		movq		%[o5]	,%%r12				\n\t"\
		"movq		%[o3]	,%%rdx				\n\t		movq		%[o7]	,%%r13				\n\t"\
		"vsubpd	%%ymm2		,%%ymm0,%%ymm0		\n\t		vmovaps	%%ymm12		,    (%%r13)	\n\t"\
		"vsubpd	%%ymm3		,%%ymm1,%%ymm1		\n\t		vmovaps	%%ymm13		,0x20(%%r13)	\n\t"\
	/* Use the cosine term of the [c1,s1] pair, which is the *middle* [4th of 7] of our 7 input pairs, in terms \
	of the input-arg bit-reversal reordering defined in the __X[c,s] --> [c,s] mapping below and happens to \
	always in fact *be* a true cosine term, which is a requirement for our "decr 1 gives isrt2" data-copy scheme: */\
		"												movq	%[c1],%%r14	\n\t"\
		"vsubpd	%%ymm7		,%%ymm4,%%ymm4		\n\t	subq	$0x20,%%r14	\n\t"/* isrt2 in [c1]-1 */\
		"vsubpd	%%ymm6		,%%ymm5,%%ymm5		\n\t		vmovaps	%%ymm10		,%%ymm13			\n\t"\
		"vaddpd	    (%%rax)	,%%ymm2,%%ymm2		\n\t		vsubpd	%%ymm11		,%%ymm10,%%ymm10	\n\t"\
		"vaddpd	0x20(%%rax)	,%%ymm3,%%ymm3		\n\t		vaddpd	%%ymm11		,%%ymm13,%%ymm13	\n\t"\
		"vaddpd	    (%%rbx)	,%%ymm7,%%ymm7		\n\t		vmulpd	    (%%r14)	,%%ymm10,%%ymm10	\n\t"\
		"vaddpd	0x20(%%rbx)	,%%ymm6,%%ymm6		\n\t		vmulpd	    (%%r14)	,%%ymm13,%%ymm13	\n\t"\
		"vmovaps	%%ymm2,    (%%rax)	/* [o0].re */\n\t	vmovaps	0x20(%%r13)	,%%ymm11			\n\t"\
		"vmovaps	%%ymm3,0x20(%%rax)	/* [o0].im */\n\t	vmovaps	%%ymm15		,%%ymm12			\n\t"\
		"vmovaps	%%ymm4,    (%%rbx)	/* [o2].re */\n\t	vaddpd	%%ymm14		,%%ymm12,%%ymm12	\n\t"\
		"vmovaps	%%ymm6,0x20(%%rbx)	/* [o2].im */\n\t	vsubpd	%%ymm14		,%%ymm15,%%ymm15	\n\t"\
		"vmovaps	%%ymm0,    (%%rcx)	/* [o1].re */\n\t	vmulpd	    (%%r14)	,%%ymm12,%%ymm12	\n\t"\
		"vmovaps	%%ymm1,0x20(%%rcx)	/* [o1].im */\n\t	vmulpd	    (%%r14)	,%%ymm15,%%ymm15	\n\t"\
		"vmovaps	%%ymm7,    (%%rdx)	/* [o3].re */\n\t	vmovaps		(%%r13)	,%%ymm14			\n\t"\
		"vmovaps	%%ymm5,0x20(%%rdx)	/* [o3].im */\n\t	vsubpd	%%ymm10		,%%ymm8 ,%%ymm8 	\n\t"\
		"													vsubpd	%%ymm13		,%%ymm9 ,%%ymm9 	\n\t"\
		"													vsubpd	%%ymm12		,%%ymm14,%%ymm14	\n\t"\
		"													vsubpd	%%ymm15		,%%ymm11,%%ymm11	\n\t"\
		"													vaddpd	    (%%r10)	,%%ymm10,%%ymm10	\n\t"\
		"													vaddpd	0x20(%%r10)	,%%ymm13,%%ymm13	\n\t"\
		"													vaddpd	    (%%r11)	,%%ymm12,%%ymm12	\n\t"\
		"													vaddpd	0x20(%%r11)	,%%ymm15,%%ymm15	\n\t"\
		"													vmovaps	%%ymm10,    (%%r10)	\n\t"/* [o4].re */\
		"													vmovaps	%%ymm13,0x20(%%r10)	\n\t"/* [o4].im */\
		"													vmovaps	%%ymm14,    (%%r11)	\n\t"/* [o6].re */\
		"													vmovaps	%%ymm11,0x20(%%r11)	\n\t"/* [o6].im */\
		"													vmovaps	%%ymm8 ,    (%%r12)	\n\t"/* [o5].re */\
		"													vmovaps	%%ymm9 ,0x20(%%r12)	\n\t"/* [o5].im */\
		"													vmovaps	%%ymm12,    (%%r13)	\n\t"/* [o7].re */\
		"													vmovaps	%%ymm15,0x20(%%r13)	\n\t"/* [o7].im */\
		:					/* outputs: none */\
		: [i0] "m" (Xi0)	/* All inputs from memory addresses here */\
		 ,[i1] "m" (Xi1)\
		 ,[i2] "m" (Xi2)\
		 ,[i3] "m" (Xi3)\
		 ,[i4] "m" (Xi4)\
		 ,[i5] "m" (Xi5)\
		 ,[i6] "m" (Xi6)\
		 ,[i7] "m" (Xi7)\
		 ,[o0] "m" (Xo0)\
		 ,[o1] "m" (Xo1)\
		 ,[o2] "m" (Xo2)\
		 ,[o3] "m" (Xo3)\
		 ,[o4] "m" (Xo4)\
		 ,[o5] "m" (Xo5)\
		 ,[o6] "m" (Xo6)\
		 ,[o7] "m" (Xo7)\
		 ,[c4] "m" (Xc1),[s4] "m" (Xs1)\
		 ,[c2] "m" (Xc2),[s2] "m" (Xs2)\
		 ,[c6] "m" (Xc3),[s6] "m" (Xs3)\
		 ,[c1] "m" (Xc4),[s1] "m" (Xs4)\
		 ,[c5] "m" (Xc5),[s5] "m" (Xs5)\
		 ,[c3] "m" (Xc6),[s3] "m" (Xs6)\
		 ,[c7] "m" (Xc7),[s7] "m" (Xs7)\
		: "cc","memory","rax","rbx","rcx","rdx","rsi","rdi","r10","r11","r12","r13","r14","r15","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
	);\
	}

	#define SSE2_RADIX8_DIT_TWIDDLE_OOP(Xi0,Xi1,Xi2,Xi3,Xi4,Xi5,Xi6,Xi7, Xo0,Xo1,Xo2,Xo3,Xo4,Xo5,Xo6,Xo7 ,Xc1,Xs1,Xc2,Xs2,Xc3,Xs3,Xc4,Xs4,Xc5,Xs5,Xc6,Xs6,Xc7,Xs7)\
	{\
	__asm__ volatile (\
	/* Block 0/1 has just one twiddle-CMUL: */\
		"movq		%[i0],%%rax				\n\t"\
		"movq		%[i1],%%rbx				\n\t"\
		"movq		%[c1],%%rdi				\n\t"/* [rdi,rsi] point to [cos,sin] components of each sincos pair, */\
		"movq		%[s1],%%rsi				\n\t"/* which is not really a pair here in terms of relative addressing: */\
		"vmovaps	    (%%rbx),%%ymm4 		\n\t	vmovaps		0x20(%%rbx),%%ymm5 		\n\t"/* _r4  = __tr1;	_r5  = __ti1; */\
		"vmovaps	    (%%rax),%%ymm0 		\n\t	vmovaps		0x20(%%rax),%%ymm1 		\n\t"/* _r0  = __tr0;	_r1  = __ti0; */\
		"vmovaps	%%ymm5 ,%%ymm6 			\n\t	vmovaps		%%ymm4 ,%%ymm7 			\n\t"/* _r6  = _r5;		_r7  = _r4;		** [r4,r5] = CMUL(__t1,__W1): */\
		"vmulpd		(%%rdi),%%ymm4 ,%%ymm4 	\n\t	vmulpd		(%%rdi),%%ymm5 ,%%ymm5 	\n\t"/* _r4 *= __Wr1;	_r5 *= __Wr1; */\
		"vmulpd		(%%rsi),%%ymm6 ,%%ymm6 	\n\t	vmulpd		(%%rsi),%%ymm7 ,%%ymm7 	\n\t"/* _r6 *= __Wi1;	_r7 *= __Wi1; */\
		"vaddpd		%%ymm6 ,%%ymm4 ,%%ymm4 	\n\t	vsubpd		%%ymm7 ,%%ymm5 ,%%ymm5 	\n\t"/* _r4 += _r6;		_r5 -= _r7; */\
		"vmovaps	%%ymm0 ,%%ymm2 			\n\t	vmovaps		%%ymm1 ,%%ymm3 			\n\t"/* _r2  = _r0;		_r3  = _r1; */\
		"vaddpd		%%ymm4 ,%%ymm0 ,%%ymm0 	\n\t	vaddpd		%%ymm5 ,%%ymm1 ,%%ymm1 	\n\t"/* _r0 += _r4;		_r1 += _r5; */\
		"vsubpd		%%ymm4 ,%%ymm2 ,%%ymm2 	\n\t	vsubpd		%%ymm5 ,%%ymm3 ,%%ymm3 	\n\t"/* _r2 -= _r4;		_r3 -= _r5; */\
		"vmovaps	%%ymm0 ,    (%%rax)		\n\t	vmovaps		%%ymm1 ,0x20(%%rax)		\n\t"/* __tr0 = _r0;	__ti0 = _r1; */\
		"vmovaps	%%ymm2 ,    (%%rbx)		\n\t	vmovaps		%%ymm3 ,0x20(%%rbx)		\n\t"/* __tr1 = _r2;	__ti1 = _r3; */\
	/* Blocks 2/3 use separate register subset, can be done overlapped with 0/1: */\
		"movq		%[i2],%%rcx				\n\t"\
		"movq		%[c2],%%r10				\n\t"\
		"movq		%[s2],%%r11				\n\t"/* [r8,r9] = CMUL(__t2,__W2): */\
		"vmovaps		(%%rcx),%%ymm8 		\n\t	vmovaps		0x20(%%rcx),%%ymm9 		\n\t"/* _r8  = __tr2;	_r9  = __ti2; */\
		"vmovaps	%%ymm9 ,%%ymm10			\n\t	vmovaps		%%ymm8 ,%%ymm11			\n\t"/* _ra  = _r9;		_rb  = _r8; */\
		"vmulpd		(%%r10),%%ymm8 ,%%ymm8 	\n\t	vmulpd		(%%r10),%%ymm9 ,%%ymm9 	\n\t"/* _r8 *= __Wr2;	_r9 *= __Wr2; */\
		"vmulpd		(%%r11),%%ymm10,%%ymm10	\n\t	vmulpd		(%%r11),%%ymm11,%%ymm11	\n\t"/* _ra *= __Wi2;	_rb *= __Wi2; */\
		"vaddpd		%%ymm10,%%ymm8 ,%%ymm8 	\n\t	vsubpd		%%ymm11,%%ymm9 ,%%ymm9 	\n\t"/* _r8 += _ra;		_r9 -= _rb; */\
		"movq		%[i3],%%rdx				\n\t"\
		"movq		%[c3],%%r12				\n\t"\
		"movq		%[s3],%%r13				\n\t"/* [rc,rd] = CMUL(__t3,__W3): */\
		"vmovaps		(%%rdx),%%ymm12		\n\t	vmovaps		0x20(%%rdx),%%ymm13		\n\t"/* _rc  = __tr3;	_rd  = __ti3; */\
		"vmovaps	%%ymm13,%%ymm14			\n\t	vmovaps		%%ymm12,%%ymm15			\n\t"/* _re  = _rd;		_rf  = _rc; */\
		"vmulpd		(%%r12),%%ymm12,%%ymm12	\n\t	vmulpd		(%%r12),%%ymm13,%%ymm13	\n\t"/* _rc *= __Wr3;	_rd *= __Wr3; */\
		"vmulpd		(%%r13),%%ymm14,%%ymm14	\n\t	vmulpd		(%%r13),%%ymm15,%%ymm15	\n\t"/* _re *= __Wi3;	_rf *= __Wi3; */\
		"vaddpd		%%ymm14,%%ymm12,%%ymm12	\n\t	vsubpd		%%ymm15,%%ymm13,%%ymm13	\n\t"/* _rc += _re;		_rd -= _rf; */\
		/* Now do radix-2 butterfly: */\
		"vmovaps	%%ymm8 ,%%ymm10			\n\t	vmovaps		%%ymm9 ,%%ymm11			\n\t"/* _ra  = _r8;		_rb  = _r9; */\
		"vaddpd		%%ymm12,%%ymm8 ,%%ymm8 	\n\t	vaddpd		%%ymm13,%%ymm9 ,%%ymm9 	\n\t"/* _r8 += _rc;		_r9 += _rd; */\
		"vsubpd		%%ymm12,%%ymm10,%%ymm10	\n\t	vsubpd		%%ymm13,%%ymm11,%%ymm11	\n\t"/* _ra -= _rc;		_rb -= _rd; */\
		"vmovaps	%%ymm8 ,    (%%rcx)		\n\t	vmovaps		%%ymm9 ,0x20(%%rcx)		\n\t"/* __tr2 = _r8;	__ti2 = _r9; */\
		"vmovaps	%%ymm10,    (%%rdx)		\n\t	vmovaps		%%ymm11,0x20(%%rdx)		\n\t"/* __tr3 = _ra;	__ti3 = _rb; */\
	/* Blocks 4/5: */\
		"movq		%[i4],%%rax				\n\t"\
		"movq		%[c4],%%rdi				\n\t"\
		"movq		%[s4],%%rsi				\n\t"/* [r0,r1] = CMUL(__t4,__W4): */\
		"vmovaps		(%%rax),%%ymm0 		\n\t	vmovaps		0x20(%%rax),%%ymm1 		\n\t"/* _r0  = __tr4;	_r1  = __ti4; */\
		"vmovaps	%%ymm1 ,%%ymm2 			\n\t	vmovaps		%%ymm0 ,%%ymm3 			\n\t"/* _r2  = _r1;		_r3  = _r0; */\
		"vmulpd		(%%rdi),%%ymm0 ,%%ymm0 	\n\t	vmulpd		(%%rdi),%%ymm1 ,%%ymm1 	\n\t"/* _r0 *= __Wr4;	_r1 *= __Wr4; */\
		"vmulpd		(%%rsi),%%ymm2 ,%%ymm2 	\n\t	vmulpd		(%%rsi),%%ymm3 ,%%ymm3 	\n\t"/* _r2 *= __Wi4;	_r3 *= __Wi4; */\
		"vaddpd		%%ymm2 ,%%ymm0 ,%%ymm0 	\n\t	vsubpd		%%ymm3 ,%%ymm1 ,%%ymm1 	\n\t"/* _r0 += _r2;		_r1 -= _r3; */\
		"movq		%[i5],%%rbx				\n\t"\
		"movq		%[c5],%%r8 				\n\t"\
		"movq		%[s5],%%r9 				\n\t"/* [r4,r5] = CMUL(__t5,__W5): */\
		"vmovaps		(%%rbx),%%ymm4 		\n\t	vmovaps		0x20(%%rbx),%%ymm5 		\n\t"/* _r4  = __tr5;	_r5  = __ti5; */\
		"vmovaps	%%ymm5 ,%%ymm6 			\n\t	vmovaps		%%ymm4 ,%%ymm7 			\n\t"/* _r6  = _r5;		_r7  = _r4; */\
		"vmulpd		(%%r8 ),%%ymm4 ,%%ymm4 	\n\t	vmulpd		(%%r8 ),%%ymm5 ,%%ymm5 	\n\t"/* _r4 *= __Wr5;	_r5 *= __Wr5; */\
		"vmulpd		(%%r9 ),%%ymm6 ,%%ymm6 	\n\t	vmulpd		(%%r9 ),%%ymm7 ,%%ymm7 	\n\t"/* _r6 *= __Wi5;	_r7 *= __Wi5; */\
		"vaddpd		%%ymm6 ,%%ymm4 ,%%ymm4 	\n\t	vsubpd		%%ymm7 ,%%ymm5 ,%%ymm5 	\n\t"/* _r4 += _r6;		_r5 -= _r7; */\
		/* Now do radix-2 butterfly: */\
		"vmovaps	%%ymm0 ,%%ymm2 			\n\t	vmovaps		%%ymm1 ,%%ymm3 			\n\t"/* _r2  = _r0;		_r3  = _r1; */\
		"vaddpd		%%ymm4 ,%%ymm0 ,%%ymm0 	\n\t	vaddpd		%%ymm5 ,%%ymm1 ,%%ymm1 	\n\t"/* _r0 += _r4;		_r1 += _r5; */\
		"vsubpd		%%ymm4 ,%%ymm2 ,%%ymm2 	\n\t	vsubpd		%%ymm5 ,%%ymm3 ,%%ymm3 	\n\t"/* _r2 -= _r4;		_r3 -= _r5; */\
	/* Blocks 6/7 use separate register subset, can be done overlapped with 4/5: */\
		"movq		%[i6],%%rcx				\n\t"\
		"movq		%[c6],%%r10				\n\t"\
		"movq		%[s6],%%r11				\n\t"/* [r8,r9] = CMUL(__t6,__W6): */\
		"vmovaps		(%%rcx),%%ymm8 		\n\t	vmovaps		0x20(%%rcx),%%ymm9 		\n\t"/* _r8  = __tr6;	_r9  = __ti6; */\
		"vmovaps	%%ymm9 ,%%ymm10			\n\t	vmovaps		%%ymm8 ,%%ymm11			\n\t"/* _ra  = _r9;		_rb  = _r8; */\
		"vmulpd		(%%r10),%%ymm8 ,%%ymm8 	\n\t	vmulpd		(%%r10),%%ymm9 ,%%ymm9 	\n\t"/* _r8 *= __Wr6;	_r9 *= __Wr6; */\
		"vmulpd		(%%r11),%%ymm10,%%ymm10	\n\t	vmulpd		(%%r11),%%ymm11,%%ymm11	\n\t"/* _ra *= __Wi6;	_rb *= __Wi6; */\
		"vaddpd		%%ymm10,%%ymm8 ,%%ymm8 	\n\t	vsubpd		%%ymm11,%%ymm9 ,%%ymm9 	\n\t"/* _r8 += _ra;		_r9 -= _rb; */\
		"movq		%[i7],%%rdx				\n\t"\
		"movq		%[c7],%%r12				\n\t"\
		"movq		%[s7],%%r13				\n\t"/* [rc,rd] = CMUL(__t7,__W7): */\
		"vmovaps		(%%rdx),%%ymm12		\n\t	vmovaps		0x20(%%rdx),%%ymm13		\n\t"/* _rc  = __tr7;	_rd  = __ti7; */\
		"vmovaps	%%ymm13,%%ymm14			\n\t	vmovaps		%%ymm12,%%ymm15			\n\t"/* _re  = _rd;		_rf  = _rc; */\
		"vmulpd		(%%r12),%%ymm12,%%ymm12	\n\t	vmulpd		(%%r12),%%ymm13,%%ymm13	\n\t"/* _rc *= __Wr7;	_rd *= __Wr7; */\
		"vmulpd		(%%r13),%%ymm14,%%ymm14	\n\t	vmulpd		(%%r13),%%ymm15,%%ymm15	\n\t"/* _re *= __Wi7;	_rf *= __Wi7; */\
		"vaddpd		%%ymm14,%%ymm12,%%ymm12	\n\t	vsubpd		%%ymm15,%%ymm13,%%ymm13	\n\t"/* _rc += _re;		_rd -= _rf; */\
		/* Now do radix-2 butterfly: */\
		"vmovaps	%%ymm8 ,%%ymm10			\n\t	vmovaps		%%ymm9 ,%%ymm11			\n\t"/* _ra  = _r8;		_rb  = _r9; */\
		"vaddpd		%%ymm12,%%ymm8 ,%%ymm8 	\n\t	vaddpd		%%ymm13,%%ymm9 ,%%ymm9 	\n\t"/* _r8 += _rc;		_r9 += _rd; */\
		"vsubpd		%%ymm12,%%ymm10,%%ymm10	\n\t	vsubpd		%%ymm13,%%ymm11,%%ymm11	\n\t"/* _ra -= _rc;		_rb -= _rd; */\
	/* Reload Block 0-3 outputs into r4-7,c-f, combine to get the 2 length-4 subtransform... */\
		"movq		%[i0],%%rax				\n\t"\
		"movq		%[i1],%%rbx				\n\t"\
		"movq		%[i2],%%rcx				\n\t"\
		"movq		%[i3],%%rdx				\n\t"\
		"vmovaps		(%%rax),%%ymm4 		\n\t	vmovaps		0x20(%%rax),%%ymm5 		\n\t"/* _r4 = __tr0;	_r5 = __ti0; */\
		"vmovaps		(%%rbx),%%ymm6 		\n\t	vmovaps		0x20(%%rbx),%%ymm7 		\n\t"/* _r6 = __tr1;	_r7 = __ti1; */\
		"vmovaps		(%%rcx),%%ymm12		\n\t	vmovaps		0x20(%%rcx),%%ymm13		\n\t"/* _rc = __tr2;	_rd = __ti2; */\
		"vmovaps		(%%rdx),%%ymm14		\n\t	vmovaps		0x20(%%rdx),%%ymm15		\n\t"/* _re = __tr3;	_rf = __ti3; */\
		"movq		%[o0],%%rax				\n\t"/* Assumes user stuck a (vec_dbl)2.0 into this output slot prior to macro call. */\
		"vsubpd		%%ymm12,%%ymm4 ,%%ymm4 	\n\t	vsubpd		%%ymm13,%%ymm5 ,%%ymm5 	\n\t"/* _r4 -= _rc;		_r5 -= _rd; */\
		"vsubpd		%%ymm15,%%ymm6 ,%%ymm6 	\n\t	vsubpd		%%ymm14,%%ymm7 ,%%ymm7 	\n\t"/* _r6 -= _rf;		_r7 -= _re; */\
		"vsubpd		%%ymm8 ,%%ymm0 ,%%ymm0 	\n\t	vsubpd		%%ymm9 ,%%ymm1 ,%%ymm1 	\n\t"/* _r0 -= _r8;		_r1 -= _r9; */\
		"vsubpd		%%ymm11,%%ymm2 ,%%ymm2 	\n\t	vsubpd		%%ymm10,%%ymm3 ,%%ymm3 	\n\t"/* _r2 -= _rb;		_r3 -= _ra; */\
		/* We hope the microcode execution engine sticks the datum at (%%rax) into a virtual register and inlines the MULs with the above SUBs: */\
		"vmulpd		(%%rax),%%ymm12,%%ymm12	\n\t	vmulpd		(%%rax),%%ymm13,%%ymm13	\n\t"/* _rc *= _two;	_rd *= _two; */\
		"vmulpd		(%%rax),%%ymm15,%%ymm15	\n\t	vmulpd		(%%rax),%%ymm14,%%ymm14	\n\t"/* _rf *= _two;	_re *= _two; */\
		"vmulpd		(%%rax),%%ymm8 ,%%ymm8 	\n\t	vmulpd		(%%rax),%%ymm9 ,%%ymm9 	\n\t"/* _r8 *= _two;	_r9 *= _two; */\
		"vmulpd		(%%rax),%%ymm11,%%ymm11	\n\t	vmulpd		(%%rax),%%ymm10,%%ymm10	\n\t"/* _rb *= _two;	_ra *= _two; */\
		"vaddpd		%%ymm4 ,%%ymm12,%%ymm12	\n\t	vaddpd		%%ymm5 ,%%ymm13,%%ymm13	\n\t"/* _rc += _r4;		_rd += _r5; */\
		"vaddpd		%%ymm6 ,%%ymm15,%%ymm15	\n\t	vaddpd		%%ymm7 ,%%ymm14,%%ymm14	\n\t"/* _rf += _r6;		_re += _r7; */\
		"vaddpd		%%ymm0 ,%%ymm8 ,%%ymm8 	\n\t	vaddpd		%%ymm1 ,%%ymm9 ,%%ymm9 	\n\t"/* _r8 += _r0;		_r9 += _r1; */\
		"vaddpd		%%ymm2 ,%%ymm11,%%ymm11	\n\t	vaddpd		%%ymm3 ,%%ymm10,%%ymm10	\n\t"/* _rb += _r2;		_ra += _r3; */\
		/* In terms of our original scalar-code prototyping macro, the data are: __tr0 = _r[c,f,4,6,8,b,0,2], __ti0 = _r[d,7,5,e,9,3,1,a]; */\
	/* Now combine the two half-transforms: */\
		/* Need r2/3+- a/b combos for the *ISRT2 preceding the output 4-7 radix-2 butterflies, so start them first: */\
		"vsubpd		%%ymm3 ,%%ymm11,%%ymm11	\n\t	vsubpd		%%ymm10,%%ymm2 ,%%ymm2 	\n\t"/* _rb -= _r3;		_r2 -= _ra; */\
		"vsubpd		%%ymm8 ,%%ymm12,%%ymm12	\n\t	vsubpd		%%ymm9 ,%%ymm13,%%ymm13	\n\t"/* _rc -= _r8;		_rd -= _r9; */\
		"vsubpd		%%ymm1 ,%%ymm4 ,%%ymm4 	\n\t	vsubpd		%%ymm0 ,%%ymm5 ,%%ymm5 	\n\t"/* _r4 -= _r1;		_r5 -= _r0; */\
		"vmulpd		(%%rax),%%ymm3 ,%%ymm3 	\n\t	vmulpd		(%%rax),%%ymm10,%%ymm10	\n\t"/* _r3 *= _two;	_ra *= _two; */\
		"vmulpd		(%%rax),%%ymm8 ,%%ymm8 	\n\t	vmulpd		(%%rax),%%ymm9 ,%%ymm9 	\n\t"/* _r8 *= _two;	_r9 *= _two; */\
		"vmulpd		(%%rax),%%ymm1 ,%%ymm1 	\n\t	vmulpd		(%%rax),%%ymm0 ,%%ymm0 	\n\t"/* _r1 *= _two;	_r0 *= _two; */\
		"vaddpd		%%ymm11,%%ymm3 ,%%ymm3 	\n\t	vaddpd		%%ymm2 ,%%ymm10,%%ymm10	\n\t"/* _r3 += _rb;		_ra += _r2; */\
		"vaddpd		%%ymm12,%%ymm8 ,%%ymm8 	\n\t	vaddpd		%%ymm13,%%ymm9 ,%%ymm9 	\n\t"/* _r8 += _rc;		_r9 += _rd; */\
		"vaddpd		%%ymm4 ,%%ymm1 ,%%ymm1 	\n\t	vaddpd		%%ymm5 ,%%ymm0 ,%%ymm0 	\n\t"/* _r1 += _r4;		_r0 += _r5; */\
		/*movq		%[o0],%%rax		[o0] already in rax */	\
		"movq		%[o1],%%rbx				\n\t"\
		"movq		%[o2],%%rcx				\n\t"\
		"movq		%[o3],%%rdx				\n\t"\
		"vmovaps	%%ymm12,    (%%rbx)		\n\t	vmovaps		%%ymm13,0x20(%%rbx)		\n\t"/* __Br1 = _rc;	__Bi1 = _rd; */\
		/* Use that _rc,d free to stick 2.0 into _rc and that [c4] in rdi to load ISRT2 from c4-1 into _rd: */\
		"vmovaps		(%%rax),%%ymm12		\n\t	vmovaps		-0x20(%%rdi),%%ymm13	\n\t"/* _rc = 2.0;		_rd = ISRT2; */\
		"vmovaps	%%ymm4 ,    (%%rdx)		\n\t	vmovaps		%%ymm0 ,0x20(%%rdx)		\n\t"/* __Br3 = _r4;	__Bi3 = _r0; */\
		"vmovaps	%%ymm8 ,    (%%rax)		\n\t	vmovaps		%%ymm9 ,0x20(%%rax)		\n\t"/* __Br0 = _r8;	__Bi0 = _r9; */\
		"vmovaps	%%ymm1 ,    (%%rcx)		\n\t	vmovaps		%%ymm5 ,0x20(%%rcx)		\n\t"/* __Br2 = _r1;	__Bi2 = _r5; */\
		"vmulpd		%%ymm13,%%ymm3 ,%%ymm3 	\n\t	vmulpd		%%ymm13,%%ymm11,%%ymm11	\n\t"/* _r3 *= ISRT2;	_rb *= ISRT2; */\
		"vmulpd		%%ymm13,%%ymm2 ,%%ymm2 	\n\t	vmulpd		%%ymm13,%%ymm10,%%ymm10	\n\t"/* _r2 *= ISRT2;	_ra *= ISRT2; */\
		"vsubpd		%%ymm3 ,%%ymm15,%%ymm15	\n\t	vsubpd		%%ymm11,%%ymm7 ,%%ymm7 	\n\t"/* _rf -= _r3;		_r7 -= _rb; */\
		"vsubpd		%%ymm2 ,%%ymm6 ,%%ymm6 	\n\t	vsubpd		%%ymm10,%%ymm14,%%ymm14	\n\t"/* _r6 -= _r2;		_re -= _ra; */\
		"vmulpd		%%ymm12,%%ymm3 ,%%ymm3 	\n\t	vmulpd		%%ymm12,%%ymm11,%%ymm11	\n\t"/* _r3 *= _two;	_rb *= _two; */\
		"vmulpd		%%ymm12,%%ymm2 ,%%ymm2 	\n\t	vmulpd		%%ymm12,%%ymm10,%%ymm10	\n\t"/* _r2 *= _two;	_ra *= _two; */\
		"vaddpd		%%ymm15,%%ymm3 ,%%ymm3 	\n\t	vaddpd		%%ymm7 ,%%ymm11,%%ymm11	\n\t"/* _r3 += _rf;		_rb += _r7; */\
		"vaddpd		%%ymm6 ,%%ymm2 ,%%ymm2 	\n\t	vaddpd		%%ymm14,%%ymm10,%%ymm10	\n\t"/* _r2 += _r6;		_ra += _re; */\
		"movq		%[o4],%%rax				\n\t"\
		"movq		%[o5],%%rbx				\n\t"\
		"movq		%[o6],%%rcx				\n\t"\
		"movq		%[o7],%%rdx				\n\t"\
		"vmovaps	%%ymm3 ,    (%%rax)		\n\t	vmovaps		%%ymm7 ,0x20(%%rax)		\n\t"/* __Br4 = _r3;	__Bi4 = _r7; */\
		"vmovaps	%%ymm15,    (%%rbx)		\n\t	vmovaps		%%ymm11,0x20(%%rbx)		\n\t"/* __Br5 = _rf;	__Bi5 = _rb; */\
		"vmovaps	%%ymm6 ,    (%%rcx)		\n\t	vmovaps		%%ymm14,0x20(%%rcx)		\n\t"/* __Br6 = _r6;	__Bi6 = _re; */\
		"vmovaps	%%ymm2 ,    (%%rdx)		\n\t	vmovaps		%%ymm10,0x20(%%rdx)		\n\t"/* __Br7 = _r2;	__Bi7 = _ra; */\
		:					/* outputs: none */\
		: [i0] "m" (Xi0)	/* All inputs from memory addresses here */\
		 ,[i1] "m" (Xi1)\
		 ,[i2] "m" (Xi2)\
		 ,[i3] "m" (Xi3)\
		 ,[i4] "m" (Xi4)\
		 ,[i5] "m" (Xi5)\
		 ,[i6] "m" (Xi6)\
		 ,[i7] "m" (Xi7)\
		 ,[o0] "m" (Xo0)\
		 ,[o1] "m" (Xo1)\
		 ,[o2] "m" (Xo2)\
		 ,[o3] "m" (Xo3)\
		 ,[o4] "m" (Xo4)\
		 ,[o5] "m" (Xo5)\
		 ,[o6] "m" (Xo6)\
		 ,[o7] "m" (Xo7)\
		 ,[c1] "m" (Xc1),[s1] "m" (Xs1)\
		 ,[c2] "m" (Xc2),[s2] "m" (Xs2)\
		 ,[c3] "m" (Xc3),[s3] "m" (Xs3)\
		 ,[c4] "m" (Xc4),[s4] "m" (Xs4)\
		 ,[c5] "m" (Xc5),[s5] "m" (Xs5)\
		 ,[c6] "m" (Xc6),[s6] "m" (Xs6)\
		 ,[c7] "m" (Xc7),[s7] "m" (Xs7)\
		: "cc","memory","rax","rbx","rcx","rdx","rsi","rdi","r8","r9","r10","r11","r12","r13","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
	);\
	}

	// Based on the SSE2_RADIX16_DIT_NOTWIDDLE macro in radix16_ditN_cy_dif1_gcc64.h, but with completely
	// specifiable 16-input addressing required for usage as the power-of-2 component of a twiddleless
	// radix = [odd*2^n] DFT routine.
	// We use just a single output base-pointer plus literal ostrides which are [1,2,3,4]-multiples of
	// __01; this allows us to cut GP-register usage, which is absolutely a must for the 32-bit version
	// of the macro, and is a benefit to the 64-bit versions which code-fold to yield 2 side-by-side
	// streams of independently executable instructions, one for data in xmm0-7, the other using xmm8-15.
	#define SSE2_RADIX16_DIT_0TWIDDLE(Xin0,Xin1,Xin2,Xin3,Xin4,Xin5,Xin6,Xin7,Xin8,Xin9,Xina,Xinb,Xinc,Xind,Xine,Xinf, Xisrt2, Xout0,Xo1,Xo2,Xo3,Xo4)\
	{\
	__asm__ volatile (\
		"movq	%[__in0],%%rax		\n\t"\
		"movq	%[__in1],%%rbx		\n\t"\
		"movq	%[__in2],%%rcx		\n\t"\
		"movq	%[__in3],%%rdx		\n\t"\
		/* SSE2_RADIX4_DIT_0TWIDDLE_B(r0 ): */\
		"vmovaps	    (%%rax),%%ymm2	\n\t"\
		"vmovaps	    (%%rcx),%%ymm6	\n\t"\
		"vmovaps	0x20(%%rax),%%ymm3	\n\t"\
		"vmovaps	0x20(%%rcx),%%ymm7	\n\t"\
		"vmovaps	    (%%rbx),%%ymm0	\n\t"\
		"vmovaps	    (%%rdx),%%ymm4	\n\t"\
		"vmovaps	0x20(%%rbx),%%ymm1	\n\t"\
		"vmovaps	0x20(%%rdx),%%ymm5	\n\t"\
		"vsubpd	%%ymm0,%%ymm2,%%ymm2		\n\t"\
		"vsubpd	%%ymm4,%%ymm6,%%ymm6		\n\t"\
		"vsubpd	%%ymm1,%%ymm3,%%ymm3		\n\t"\
		"vsubpd	%%ymm5,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	%%ymm0,%%ymm0,%%ymm0		\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	%%ymm1,%%ymm1,%%ymm1		\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	%%ymm2,%%ymm0,%%ymm0		\n\t"\
		"vaddpd	%%ymm6,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	%%ymm3,%%ymm1,%%ymm1		\n\t"\
		"vaddpd	%%ymm7,%%ymm5,%%ymm5		\n\t"\
		"movq	%[__out0],%%rsi		\n\t"\
		"vsubpd	%%ymm4,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	%%ymm7,%%ymm2,%%ymm2		\n\t"\
		"vsubpd	%%ymm5,%%ymm1,%%ymm1		\n\t"\
		"vsubpd	%%ymm6,%%ymm3,%%ymm3		\n\t"\
		"leaq	0x20(%%rsi),%%rdi	\n\t"/* Need separate address Im parts of outputs due to literal-offsets below */\
		"vmovaps	%%ymm0,%c[__o2](%%rsi)	\n\t"\
		"vmovaps	%%ymm2,%c[__o3](%%rsi)	\n\t"\
		"vmovaps	%%ymm1,%c[__o2](%%rdi)	\n\t"\
		"vmovaps	%%ymm3,%c[__o1](%%rdi)	\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	%%ymm0,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	%%ymm2,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	%%ymm1,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	%%ymm3,%%ymm6,%%ymm6		\n\t"\
		"vmovaps	%%ymm4,        (%%rsi)	\n\t"\
		"vmovaps	%%ymm7,%c[__o1](%%rsi)	\n\t"\
		"vmovaps	%%ymm5,        (%%rdi)	\n\t"\
		"vmovaps	%%ymm6,%c[__o3](%%rdi)	\n\t"\
		/* SSE2_RADIX4_DIT_0TWIDDLE_B(r8 ): */\
		"movq	%[__in4],%%rax		\n\t	leaq	%c[__o4](%%rsi),%%rsi	\n\t"/* __out0 + 4*ostride */\
		"movq	%[__in5],%%rbx		\n\t"\
		"movq	%[__in6],%%rcx		\n\t"\
		"movq	%[__in7],%%rdx		\n\t"\
		"vmovaps	    (%%rax),%%ymm2	\n\t"\
		"vmovaps	    (%%rcx),%%ymm6	\n\t"\
		"vmovaps	0x20(%%rax),%%ymm3	\n\t"\
		"vmovaps	0x20(%%rcx),%%ymm7	\n\t"\
		"vmovaps	    (%%rbx),%%ymm0	\n\t"\
		"vmovaps	    (%%rdx),%%ymm4	\n\t"\
		"vmovaps	0x20(%%rbx),%%ymm1	\n\t"\
		"vmovaps	0x20(%%rdx),%%ymm5	\n\t"\
		"vsubpd	%%ymm0,%%ymm2,%%ymm2		\n\t"\
		"vsubpd	%%ymm4,%%ymm6,%%ymm6		\n\t"\
		"vsubpd	%%ymm1,%%ymm3,%%ymm3		\n\t"\
		"vsubpd	%%ymm5,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	%%ymm0,%%ymm0,%%ymm0		\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	%%ymm1,%%ymm1,%%ymm1		\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	%%ymm2,%%ymm0,%%ymm0		\n\t"\
		"vaddpd	%%ymm6,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	%%ymm3,%%ymm1,%%ymm1		\n\t"\
		"vaddpd	%%ymm7,%%ymm5,%%ymm5		\n\t"\
		"vsubpd	%%ymm4,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	%%ymm7,%%ymm2,%%ymm2		\n\t"\
		"vsubpd	%%ymm5,%%ymm1,%%ymm1		\n\t"\
		"vsubpd	%%ymm6,%%ymm3,%%ymm3		\n\t"\
		"leaq	0x20(%%rsi),%%rdi	\n\t"\
		"vmovaps	%%ymm0,%c[__o2](%%rsi)	\n\t"\
		"vmovaps	%%ymm2,%c[__o3](%%rsi)	\n\t"\
		"vmovaps	%%ymm1,%c[__o2](%%rdi)	\n\t"\
		"vmovaps	%%ymm3,%c[__o1](%%rdi)	\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	%%ymm0,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	%%ymm2,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	%%ymm1,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	%%ymm3,%%ymm6,%%ymm6		\n\t"\
		"vmovaps	%%ymm4,        (%%rsi)	\n\t"\
		"vmovaps	%%ymm7,%c[__o1](%%rsi)	\n\t"\
		"vmovaps	%%ymm5,        (%%rdi)	\n\t"\
		"vmovaps	%%ymm6,%c[__o3](%%rdi)	\n\t"\
		/* SSE2_RADIX4_DIT_0TWIDDLE_B(r16): */\
		"movq	%[__in8],%%rax			\n\t	leaq	%c[__o4](%%rsi),%%rsi	\n\t"/* __out0 + 8*ostride */\
		"movq	%[__in9],%%rbx			\n\t"\
		"movq	%[__ina],%%rcx			\n\t"\
		"movq	%[__inb],%%rdx			\n\t"\
		"vmovaps	    (%%rax),%%ymm2	\n\t"\
		"vmovaps	    (%%rcx),%%ymm6	\n\t"\
		"vmovaps	0x20(%%rax),%%ymm3	\n\t"\
		"vmovaps	0x20(%%rcx),%%ymm7	\n\t"\
		"vmovaps	    (%%rbx),%%ymm0	\n\t"\
		"vmovaps	    (%%rdx),%%ymm4	\n\t"\
		"vmovaps	0x20(%%rbx),%%ymm1	\n\t"\
		"vmovaps	0x20(%%rdx),%%ymm5	\n\t"\
		"vsubpd	%%ymm0,%%ymm2,%%ymm2		\n\t"\
		"vsubpd	%%ymm4,%%ymm6,%%ymm6		\n\t"\
		"vsubpd	%%ymm1,%%ymm3,%%ymm3		\n\t"\
		"vsubpd	%%ymm5,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	%%ymm0,%%ymm0,%%ymm0		\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	%%ymm1,%%ymm1,%%ymm1		\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	%%ymm2,%%ymm0,%%ymm0		\n\t"\
		"vaddpd	%%ymm6,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	%%ymm3,%%ymm1,%%ymm1		\n\t"\
		"vaddpd	%%ymm7,%%ymm5,%%ymm5		\n\t"\
		"vsubpd	%%ymm4,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	%%ymm7,%%ymm2,%%ymm2		\n\t"\
		"vsubpd	%%ymm5,%%ymm1,%%ymm1		\n\t"\
		"vsubpd	%%ymm6,%%ymm3,%%ymm3		\n\t"\
		"leaq	0x20(%%rsi),%%rdi	\n\t"\
		"vmovaps	%%ymm0,%c[__o2](%%rsi)	\n\t"\
		"vmovaps	%%ymm2,%c[__o3](%%rsi)	\n\t"\
		"vmovaps	%%ymm1,%c[__o2](%%rdi)	\n\t"\
		"vmovaps	%%ymm3,%c[__o1](%%rdi)	\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	%%ymm0,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	%%ymm2,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	%%ymm1,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	%%ymm3,%%ymm6,%%ymm6		\n\t"\
		"vmovaps	%%ymm4,        (%%rsi)	\n\t"\
		"vmovaps	%%ymm7,%c[__o1](%%rsi)	\n\t"\
		"vmovaps	%%ymm5,        (%%rdi)	\n\t"\
		"vmovaps	%%ymm6,%c[__o3](%%rdi)	\n\t"\
		/* SSE2_RADIX4_DIT_0TWIDDLE_B(r24): */\
		"movq	%[__inc],%%rax			\n\t	leaq	%c[__o4](%%rsi),%%rsi	\n\t"/* __out0 + c*ostride */\
		"movq	%[__ind],%%rbx			\n\t"\
		"movq	%[__ine],%%rcx			\n\t"\
		"movq	%[__inf],%%rdx			\n\t"\
		"vmovaps	    (%%rax),%%ymm2	\n\t"\
		"vmovaps	    (%%rcx),%%ymm6	\n\t"\
		"vmovaps	0x20(%%rax),%%ymm3	\n\t"\
		"vmovaps	0x20(%%rcx),%%ymm7	\n\t"\
		"vmovaps	    (%%rbx),%%ymm0	\n\t"\
		"vmovaps	    (%%rdx),%%ymm4	\n\t"\
		"vmovaps	0x20(%%rbx),%%ymm1	\n\t"\
		"vmovaps	0x20(%%rdx),%%ymm5	\n\t"\
		"vsubpd	%%ymm0,%%ymm2,%%ymm2		\n\t"\
		"vsubpd	%%ymm4,%%ymm6,%%ymm6		\n\t"\
		"vsubpd	%%ymm1,%%ymm3,%%ymm3		\n\t"\
		"vsubpd	%%ymm5,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	%%ymm0,%%ymm0,%%ymm0		\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	%%ymm1,%%ymm1,%%ymm1		\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	%%ymm2,%%ymm0,%%ymm0		\n\t"\
		"vaddpd	%%ymm6,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	%%ymm3,%%ymm1,%%ymm1		\n\t"\
		"vaddpd	%%ymm7,%%ymm5,%%ymm5		\n\t"\
		"vsubpd	%%ymm4,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	%%ymm7,%%ymm2,%%ymm2		\n\t"\
		"vsubpd	%%ymm5,%%ymm1,%%ymm1		\n\t"\
		"vsubpd	%%ymm6,%%ymm3,%%ymm3		\n\t"\
		"leaq	0x20(%%rsi),%%rdi	\n\t"\
		"vmovaps	%%ymm0,%c[__o2](%%rsi)	\n\t"\
		"vmovaps	%%ymm2,%c[__o3](%%rsi)	\n\t"\
		"vmovaps	%%ymm1,%c[__o2](%%rdi)	\n\t"\
		"vmovaps	%%ymm3,%c[__o1](%%rdi)	\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	%%ymm0,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	%%ymm2,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	%%ymm1,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	%%ymm3,%%ymm6,%%ymm6		\n\t"\
		"vmovaps	%%ymm4,        (%%rsi)	\n\t"\
		"vmovaps	%%ymm7,%c[__o1](%%rsi)	\n\t"\
		"vmovaps	%%ymm5,        (%%rdi)	\n\t"\
		"vmovaps	%%ymm6,%c[__o3](%%rdi)	\n\t"\
	/*** Now do 4 DFTs with internal twiddles on the 4*stride - separated data: ***/\
		"movq	%[__out0],%%rax		\n\t"\
		"leaq	%c[__o4](%%rax),%%rbx	\n\t"/* __out0 +   [4*ostride] */\
		"leaq	%c[__o4](%%rbx),%%rcx	\n\t"/* __out0 + 2*[4*ostride] */\
		"leaq	%c[__o4](%%rcx),%%rdx	\n\t"/* __out0 + 3*[4*ostride] */\
		/* Block 0: */\
		"vmovaps	    (%%rax),%%ymm0	\n\t"\
		"vmovaps	0x20(%%rax),%%ymm1	\n\t"\
		"vmovaps	    (%%rbx),%%ymm2	\n\t"\
		"vmovaps	0x20(%%rbx),%%ymm3	\n\t"\
		"vsubpd	    (%%rbx),%%ymm0,%%ymm0	\n\t"\
		"vsubpd	0x20(%%rbx),%%ymm1,%%ymm1	\n\t"\
		"vaddpd	    (%%rax),%%ymm2,%%ymm2	\n\t"\
		"vaddpd	0x20(%%rax),%%ymm3,%%ymm3	\n\t"\
		"vmovaps	    (%%rcx),%%ymm4	\n\t"\
		"vmovaps	0x20(%%rcx),%%ymm5	\n\t"\
		"vmovaps	    (%%rdx),%%ymm6	\n\t"\
		"vmovaps	0x20(%%rdx),%%ymm7	\n\t"\
		"vsubpd	    (%%rdx),%%ymm4,%%ymm4	\n\t"\
		"vsubpd	0x20(%%rdx),%%ymm5,%%ymm5	\n\t"\
		"vaddpd	    (%%rcx),%%ymm6,%%ymm6	\n\t"\
		"vaddpd	0x20(%%rcx),%%ymm7,%%ymm7	\n\t"\
		"vsubpd	%%ymm6,%%ymm2,%%ymm2		\n\t"\
		"vsubpd	%%ymm7,%%ymm3,%%ymm3		\n\t"\
		"vmovaps	%%ymm2,    (%%rcx)	\n\t"\
		"vmovaps	%%ymm3,0x20(%%rcx)	\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	%%ymm2,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	%%ymm3,%%ymm7,%%ymm7		\n\t"\
		"vmovaps	%%ymm6,    (%%rax)	\n\t"\
		"vmovaps	%%ymm7,0x20(%%rax)	\n\t"\
		"vsubpd	%%ymm5,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	%%ymm4,%%ymm1,%%ymm1		\n\t"\
		"vmovaps	%%ymm0,    (%%rdx)	\n\t"\
		"vmovaps	%%ymm1,0x20(%%rbx)	\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	%%ymm0,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	%%ymm1,%%ymm4,%%ymm4		\n\t"\
		"vmovaps	%%ymm5,    (%%rbx)	\n\t"\
		"vmovaps	%%ymm4,0x20(%%rdx)	\n\t"\
		/* Block 2: */\
		"leaq	%c[__o2](%%rax),%%rax	\n\t"/* All addresses += 2*ostride */\
		"leaq	%c[__o2](%%rbx),%%rbx	\n\t"\
		"leaq	%c[__o2](%%rcx),%%rcx	\n\t"\
		"leaq	%c[__o2](%%rdx),%%rdx	\n\t"\
		"movq	%[__isrt2],%%rdi	\n\t"\
		"vmovaps	(%%rdi),%%ymm2		\n\t"\
		"vmovaps	    (%%rcx),%%ymm4	\n\t"\
		"vmovaps	0x20(%%rcx),%%ymm5	\n\t"\
		"vmovaps	    (%%rdx),%%ymm0	\n\t"\
		"vmovaps	0x20(%%rdx),%%ymm1	\n\t"\
		"vaddpd	0x20(%%rcx),%%ymm4,%%ymm4	\n\t"\
		"vsubpd	    (%%rcx),%%ymm5,%%ymm5	\n\t"\
		"vsubpd	0x20(%%rdx),%%ymm0,%%ymm0	\n\t"\
		"vaddpd	    (%%rdx),%%ymm1,%%ymm1	\n\t"\
		"vmulpd	%%ymm2,%%ymm4,%%ymm4					\n\t"\
		"vmulpd	%%ymm2,%%ymm5,%%ymm5					\n\t"\
		"vmulpd	%%ymm2,%%ymm0,%%ymm0					\n\t"\
		"vmulpd	%%ymm2,%%ymm1,%%ymm1					\n\t"\
		"vmovaps	%%ymm4,%%ymm6		\n\t"\
		"vmovaps	%%ymm5,%%ymm7		\n\t"\
		"vsubpd	%%ymm0,%%ymm4,%%ymm4					\n\t"\
		"vsubpd	%%ymm1,%%ymm5,%%ymm5					\n\t"\
		"vaddpd	%%ymm0,%%ymm6,%%ymm6					\n\t"\
		"vaddpd	%%ymm1,%%ymm7,%%ymm7					\n\t"\
		"vmovaps	    (%%rax),%%ymm0	\n\t"\
		"vmovaps	0x20(%%rax),%%ymm1	\n\t"\
		"vmovaps	    (%%rbx),%%ymm2	\n\t"\
		"vmovaps	0x20(%%rbx),%%ymm3	\n\t"\
		"vsubpd	0x20(%%rbx),%%ymm0,%%ymm0	\n\t"\
		"vsubpd	    (%%rbx),%%ymm1,%%ymm1	\n\t"\
		"vaddpd	    (%%rax),%%ymm3,%%ymm3	\n\t"\
		"vaddpd	0x20(%%rax),%%ymm2,%%ymm2	\n\t"\
		"vsubpd	%%ymm4,%%ymm3,%%ymm3		\n\t"\
		"vsubpd	%%ymm5,%%ymm1,%%ymm1		\n\t"\
		"vmovaps	%%ymm3,    (%%rcx)	\n\t"\
		"vmovaps	%%ymm1,0x20(%%rcx)	\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	%%ymm3,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	%%ymm1,%%ymm5,%%ymm5		\n\t"\
		"vmovaps	%%ymm4,    (%%rax)	\n\t"\
		"vmovaps	%%ymm5,0x20(%%rax)	\n\t"\
		"vsubpd	%%ymm7,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	%%ymm6,%%ymm2,%%ymm2		\n\t"\
		"vmovaps	%%ymm0,    (%%rdx)	\n\t"\
		"vmovaps	%%ymm2,0x20(%%rbx)	\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	%%ymm0,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	%%ymm2,%%ymm6,%%ymm6		\n\t"\
		"vmovaps	%%ymm7,    (%%rbx)	\n\t"\
		"vmovaps	%%ymm6,0x20(%%rdx)	\n\t"\
		/* Block 1: */\
		"leaq	-%c[__o1](%%rax),%%rax	\n\t"/* All addresses -= 1*ostride */\
		"leaq	-%c[__o1](%%rbx),%%rbx	\n\t"\
		"leaq	-%c[__o1](%%rcx),%%rcx	\n\t"\
		"leaq	-%c[__o1](%%rdx),%%rdx	\n\t"\
		"leaq	0x20(%%rdi),%%rsi	\n\t"/* cc0 */\
		"vmovaps	    (%%rdx),%%ymm0	\n\t"\
		"vmovaps	0x20(%%rdx),%%ymm1	\n\t"\
		"vmovaps	    (%%rdx),%%ymm2	\n\t"\
		"vmovaps	0x20(%%rdx),%%ymm3	\n\t"\
		"vmulpd	0x20(%%rsi),%%ymm0,%%ymm0	\n\t"\
		"vmulpd	0x20(%%rsi),%%ymm1,%%ymm1	\n\t"\
		"vmulpd	    (%%rsi),%%ymm2,%%ymm2	\n\t"\
		"vmulpd	    (%%rsi),%%ymm3,%%ymm3	\n\t"\
		"vsubpd	%%ymm2,%%ymm1,%%ymm1		\n\t"\
		"vaddpd	%%ymm3,%%ymm0,%%ymm0		\n\t"\
		"vmovaps	    (%%rcx),%%ymm4	\n\t"\
		"vmovaps	0x20(%%rcx),%%ymm5	\n\t"\
		"vmovaps	    (%%rcx),%%ymm6	\n\t"\
		"vmovaps	0x20(%%rcx),%%ymm7	\n\t"\
		"vmulpd	    (%%rsi),%%ymm4,%%ymm4	\n\t"\
		"vmulpd	    (%%rsi),%%ymm5,%%ymm5	\n\t"\
		"vmulpd	0x20(%%rsi),%%ymm6,%%ymm6	\n\t"\
		"vmulpd	0x20(%%rsi),%%ymm7,%%ymm7	\n\t"\
		"vsubpd	%%ymm6,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	%%ymm7,%%ymm4,%%ymm4		\n\t"\
		"vmovaps	%%ymm5,%%ymm7		\n\t"\
		"vmovaps	%%ymm4,%%ymm6		\n\t"\
		"vaddpd	%%ymm0,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	%%ymm1,%%ymm5,%%ymm5		\n\t"\
		"vsubpd	%%ymm0,%%ymm6,%%ymm6		\n\t"\
		"vsubpd	%%ymm1,%%ymm7,%%ymm7		\n\t"\
		"vmovaps	    (%%rbx),%%ymm2	\n\t"\
		"vmovaps	0x20(%%rbx),%%ymm3	\n\t"\
		"vmovaps	    (%%rax),%%ymm0	\n\t"\
		"vmovaps	0x20(%%rax),%%ymm1	\n\t"\
		"vaddpd	0x20(%%rbx),%%ymm2,%%ymm2	\n\t"\
		"vsubpd	    (%%rbx),%%ymm3,%%ymm3	\n\t"\
		"vmulpd	    (%%rdi),%%ymm2,%%ymm2	\n\t"\
		"vmulpd	    (%%rdi),%%ymm3,%%ymm3	\n\t"\
		"vsubpd	%%ymm2,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	%%ymm3,%%ymm1,%%ymm1		\n\t"\
		"vaddpd	%%ymm2,%%ymm2,%%ymm2		\n\t"\
		"vaddpd	%%ymm3,%%ymm3,%%ymm3		\n\t"\
		"vaddpd	%%ymm0,%%ymm2,%%ymm2		\n\t"\
		"vaddpd	%%ymm1,%%ymm3,%%ymm3		\n\t"\
		"vsubpd	%%ymm4,%%ymm2,%%ymm2		\n\t"\
		"vsubpd	%%ymm5,%%ymm3,%%ymm3		\n\t"\
		"vmovaps	%%ymm2,    (%%rcx)	\n\t"\
		"vmovaps	%%ymm3,0x20(%%rcx)	\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	%%ymm2,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	%%ymm3,%%ymm5,%%ymm5		\n\t"\
		"vmovaps	%%ymm4,    (%%rax)	\n\t"\
		"vmovaps	%%ymm5,0x20(%%rax)	\n\t"\
		"vsubpd	%%ymm7,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	%%ymm6,%%ymm1,%%ymm1		\n\t"\
		"vmovaps	%%ymm0,    (%%rdx)	\n\t"\
		"vmovaps	%%ymm1,0x20(%%rbx)	\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	%%ymm0,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	%%ymm1,%%ymm6,%%ymm6		\n\t"\
		"vmovaps	%%ymm7,    (%%rbx)	\n\t"\
		"vmovaps	%%ymm6,0x20(%%rdx)	\n\t"\
		/* Block 3: */\
		"leaq	%c[__o2](%%rax),%%rax	\n\t"/* All addresses += 2*ostride */\
		"leaq	%c[__o2](%%rbx),%%rbx	\n\t"\
		"leaq	%c[__o2](%%rcx),%%rcx	\n\t"\
		"leaq	%c[__o2](%%rdx),%%rdx	\n\t"\
		"vmovaps	    (%%rdx),%%ymm0	\n\t"\
		"vmovaps	0x20(%%rdx),%%ymm1	\n\t"\
		"vmovaps	    (%%rdx),%%ymm2	\n\t"\
		"vmovaps	0x20(%%rdx),%%ymm3	\n\t"\
		"vmulpd	    (%%rsi),%%ymm0,%%ymm0	\n\t"\
		"vmulpd	    (%%rsi),%%ymm1,%%ymm1	\n\t"\
		"vmulpd	0x20(%%rsi),%%ymm2,%%ymm2	\n\t"\
		"vmulpd	0x20(%%rsi),%%ymm3,%%ymm3	\n\t"\
		"vsubpd	%%ymm2,%%ymm1,%%ymm1		\n\t"\
		"vaddpd	%%ymm3,%%ymm0,%%ymm0		\n\t"\
		"vmovaps	    (%%rcx),%%ymm4	\n\t"\
		"vmovaps	0x20(%%rcx),%%ymm5	\n\t"\
		"vmovaps	    (%%rcx),%%ymm6	\n\t"\
		"vmovaps	0x20(%%rcx),%%ymm7	\n\t"\
		"vmulpd	0x20(%%rsi),%%ymm4,%%ymm4	\n\t"\
		"vmulpd	0x20(%%rsi),%%ymm5,%%ymm5	\n\t"\
		"vmulpd	    (%%rsi),%%ymm6,%%ymm6	\n\t"\
		"vmulpd	    (%%rsi),%%ymm7,%%ymm7	\n\t"\
		"vsubpd	%%ymm6,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	%%ymm7,%%ymm4,%%ymm4		\n\t"\
		"vmovaps	%%ymm5,%%ymm7		\n\t"\
		"vmovaps	%%ymm4,%%ymm6		\n\t"\
		"vaddpd	%%ymm0,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	%%ymm1,%%ymm5,%%ymm5		\n\t"\
		"vsubpd	%%ymm0,%%ymm6,%%ymm6		\n\t"\
		"vsubpd	%%ymm1,%%ymm7,%%ymm7		\n\t"\
		"vmovaps	    (%%rbx),%%ymm2	\n\t"\
		"vmovaps	0x20(%%rbx),%%ymm3	\n\t"\
		"vmovaps	    (%%rax),%%ymm0	\n\t"\
		"vmovaps	0x20(%%rax),%%ymm1	\n\t"\
		"vsubpd	0x20(%%rbx),%%ymm2,%%ymm2	\n\t"\
		"vaddpd	    (%%rbx),%%ymm3,%%ymm3	\n\t"\
		"vmulpd	    (%%rdi),%%ymm2,%%ymm2	\n\t"\
		"vmulpd	    (%%rdi),%%ymm3,%%ymm3	\n\t"\
		"vsubpd	%%ymm2,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	%%ymm3,%%ymm1,%%ymm1		\n\t"\
		"vaddpd	%%ymm2,%%ymm2,%%ymm2		\n\t"\
		"vaddpd	%%ymm3,%%ymm3,%%ymm3		\n\t"\
		"vaddpd	%%ymm0,%%ymm2,%%ymm2		\n\t"\
		"vaddpd	%%ymm1,%%ymm3,%%ymm3		\n\t"\
		"vsubpd	%%ymm6,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	%%ymm7,%%ymm1,%%ymm1		\n\t"\
		"vmovaps	%%ymm0,    (%%rcx)	\n\t"\
		"vmovaps	%%ymm1,0x20(%%rcx)	\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	%%ymm0,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	%%ymm1,%%ymm7,%%ymm7		\n\t"\
		"vmovaps	%%ymm6,    (%%rax)	\n\t"\
		"vmovaps	%%ymm7,0x20(%%rax)	\n\t"\
		"vsubpd	%%ymm5,%%ymm2,%%ymm2		\n\t"\
		"vsubpd	%%ymm4,%%ymm3,%%ymm3		\n\t"\
		"vmovaps	%%ymm2,    (%%rdx)	\n\t"\
		"vmovaps	%%ymm3,0x20(%%rbx)	\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	%%ymm2,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	%%ymm3,%%ymm4,%%ymm4		\n\t"\
		"vmovaps	%%ymm5,    (%%rbx)	\n\t"\
		"vmovaps	%%ymm4,0x20(%%rdx)	\n\t"\
		:					/* outputs: none */\
		:[__in0] "m" (Xin0)	/* All inputs from memory addresses here */\
		,[__in1] "m" (Xin1)\
		,[__in2] "m" (Xin2)\
		,[__in3] "m" (Xin3)\
		,[__in4] "m" (Xin4)\
		,[__in5] "m" (Xin5)\
		,[__in6] "m" (Xin6)\
		,[__in7] "m" (Xin7)\
		,[__in8] "m" (Xin8)\
		,[__in9] "m" (Xin9)\
		,[__ina] "m" (Xina)\
		,[__inb] "m" (Xinb)\
		,[__inc] "m" (Xinc)\
		,[__ind] "m" (Xind)\
		,[__ine] "m" (Xine)\
		,[__inf] "m" (Xinf)\
		,[__isrt2] "m" (Xisrt2)\
		,[__out0] "m" (Xout0)\
		,[__o1] "e" (Xo1)\
		,[__o2] "e" (Xo2)\
		,[__o3] "e" (Xo3)\
		,[__o4] "e" (Xo4)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
	);\
	}

	// Based on the SSE2_RADIX16_DIF_NOTWIDDLE macro in radix16_ditN_cy_dif1_gcc64.h, but with completely
	// specifiable 16-output addressing required for usage as the power-of-2 component of a twiddleless
	// radix = [odd*2^n] DFT routine.
	#define SSE2_RADIX16_DIF_0TWIDDLE(Xin0,Xi1,Xi2,Xi3,Xi4, Xisrt2, Xout0,Xout1,Xout2,Xout3,Xout4,Xout5,Xout6,Xout7,Xout8,Xout9,Xouta,Xoutb,Xoutc,Xoutd,Xoute,Xoutf)\
	{\
	__asm__ volatile (\
		/* SSE2_RADIX4_DIF_IN_PLACE(r1 , r17, r9 , r25): */\
		"movq	%[__in0],%%rax		\n\t"\
		"leaq	%c[__i4](%%rax),%%rcx	\n\t"/* __in0 +   [4*istride]; note BR of [a,b,c,d]-ptrs, i.e. b/c swap */\
		"leaq	%c[__i4](%%rcx),%%rbx	\n\t"/* __in0 + 2*[4*istride] */\
		"leaq	%c[__i4](%%rbx),%%rdx	\n\t"/* __in0 + 3*[4*istride] */\
		"vmovaps	    (%%rax),%%ymm0	\n\t"\
		"vmovaps	0x20(%%rax),%%ymm1	\n\t"\
		"vmovaps	    (%%rax),%%ymm2	\n\t"\
		"vmovaps	0x20(%%rax),%%ymm3	\n\t"\
		"vaddpd	    (%%rbx),%%ymm0,%%ymm0	\n\t"\
		"vaddpd	0x20(%%rbx),%%ymm1,%%ymm1	\n\t"\
		"vsubpd	    (%%rbx),%%ymm2,%%ymm2	\n\t"\
		"vsubpd	0x20(%%rbx),%%ymm3,%%ymm3	\n\t"\
		"vmovaps	    (%%rcx),%%ymm4	\n\t"\
		"vmovaps	0x20(%%rcx),%%ymm5	\n\t"\
		"vmovaps	    (%%rcx),%%ymm6	\n\t"\
		"vmovaps	0x20(%%rcx),%%ymm7	\n\t"\
		"vaddpd	    (%%rdx),%%ymm4,%%ymm4	\n\t"\
		"vaddpd	0x20(%%rdx),%%ymm5,%%ymm5	\n\t"\
		"vsubpd	    (%%rdx),%%ymm6,%%ymm6	\n\t"\
		"vsubpd	0x20(%%rdx),%%ymm7,%%ymm7	\n\t"\
		"vsubpd	%%ymm4,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	%%ymm5,%%ymm1,%%ymm1		\n\t"\
		"vmovaps	%%ymm0,    (%%rbx)	\n\t"\
		"vmovaps	%%ymm1,0x20(%%rbx)	\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	%%ymm0,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	%%ymm1,%%ymm5,%%ymm5		\n\t"\
		"vmovaps	%%ymm4,    (%%rax)	\n\t"\
		"vmovaps	%%ymm5,0x20(%%rax)	\n\t"\
		"vsubpd	%%ymm7,%%ymm2,%%ymm2		\n\t"\
		"vsubpd	%%ymm6,%%ymm3,%%ymm3		\n\t"\
		"vmovaps	%%ymm2,    (%%rcx)	\n\t"\
		"vmovaps	%%ymm3,0x20(%%rdx)	\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	%%ymm2,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	%%ymm3,%%ymm6,%%ymm6		\n\t"\
		"vmovaps	%%ymm7,    (%%rdx)	\n\t"\
		"vmovaps	%%ymm6,0x20(%%rcx)	\n\t"\
		/* SSE2_RADIX4_DIF_IN_PLACE(r5 , r21, r13, r29): */\
		"leaq	%c[__i2](%%rax),%%rax	\n\t"/* All addresses += 2*ostride */\
		"leaq	%c[__i2](%%rbx),%%rbx	\n\t"\
		"leaq	%c[__i2](%%rcx),%%rcx	\n\t"\
		"leaq	%c[__i2](%%rdx),%%rdx	\n\t"\
		"vmovaps	    (%%rax),%%ymm0	\n\t"\
		"vmovaps	0x20(%%rax),%%ymm1	\n\t"\
		"vmovaps	    (%%rax),%%ymm2	\n\t"\
		"vmovaps	0x20(%%rax),%%ymm3	\n\t"\
		"vaddpd	    (%%rbx),%%ymm0,%%ymm0	\n\t"\
		"vaddpd	0x20(%%rbx),%%ymm1,%%ymm1	\n\t"\
		"vsubpd	    (%%rbx),%%ymm2,%%ymm2	\n\t"\
		"vsubpd	0x20(%%rbx),%%ymm3,%%ymm3	\n\t"\
		"vmovaps	    (%%rcx),%%ymm4	\n\t"\
		"vmovaps	0x20(%%rcx),%%ymm5	\n\t"\
		"vmovaps	    (%%rcx),%%ymm6	\n\t"\
		"vmovaps	0x20(%%rcx),%%ymm7	\n\t"\
		"vaddpd	    (%%rdx),%%ymm4,%%ymm4	\n\t"\
		"vaddpd	0x20(%%rdx),%%ymm5,%%ymm5	\n\t"\
		"vsubpd	    (%%rdx),%%ymm6,%%ymm6	\n\t"\
		"vsubpd	0x20(%%rdx),%%ymm7,%%ymm7	\n\t"\
		"vsubpd	%%ymm4,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	%%ymm5,%%ymm1,%%ymm1		\n\t"\
		"vmovaps	%%ymm0,    (%%rbx)	\n\t"\
		"vmovaps	%%ymm1,0x20(%%rbx)	\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	%%ymm0,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	%%ymm1,%%ymm5,%%ymm5		\n\t"\
		"vmovaps	%%ymm4,    (%%rax)	\n\t"\
		"vmovaps	%%ymm5,0x20(%%rax)	\n\t"\
		"vsubpd	%%ymm7,%%ymm2,%%ymm2		\n\t"\
		"vsubpd	%%ymm6,%%ymm3,%%ymm3		\n\t"\
		"vmovaps	%%ymm2,    (%%rcx)	\n\t"\
		"vmovaps	%%ymm3,0x20(%%rdx)	\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	%%ymm2,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	%%ymm3,%%ymm6,%%ymm6		\n\t"\
		"vmovaps	%%ymm7,    (%%rdx)	\n\t"\
		"vmovaps	%%ymm6,0x20(%%rcx)	\n\t"\
		/* SSE2_RADIX4_DIF_IN_PLACE(r3 , r19, r11, r27): */\
		"leaq	-%c[__i1](%%rax),%%rax	\n\t"/* All addresses -= 1*ostride */\
		"leaq	-%c[__i1](%%rbx),%%rbx	\n\t"\
		"leaq	-%c[__i1](%%rcx),%%rcx	\n\t"\
		"leaq	-%c[__i1](%%rdx),%%rdx	\n\t"\
		"vmovaps	    (%%rax),%%ymm0	\n\t"\
		"vmovaps	0x20(%%rax),%%ymm1	\n\t"\
		"vmovaps	    (%%rax),%%ymm2	\n\t"\
		"vmovaps	0x20(%%rax),%%ymm3	\n\t"\
		"vaddpd	    (%%rbx),%%ymm0,%%ymm0	\n\t"\
		"vaddpd	0x20(%%rbx),%%ymm1,%%ymm1	\n\t"\
		"vsubpd	    (%%rbx),%%ymm2,%%ymm2	\n\t"\
		"vsubpd	0x20(%%rbx),%%ymm3,%%ymm3	\n\t"\
		"vmovaps	    (%%rcx),%%ymm4	\n\t"\
		"vmovaps	0x20(%%rcx),%%ymm5	\n\t"\
		"vmovaps	    (%%rcx),%%ymm6	\n\t"\
		"vmovaps	0x20(%%rcx),%%ymm7	\n\t"\
		"vaddpd	    (%%rdx),%%ymm4,%%ymm4	\n\t"\
		"vaddpd	0x20(%%rdx),%%ymm5,%%ymm5	\n\t"\
		"vsubpd	    (%%rdx),%%ymm6,%%ymm6	\n\t"\
		"vsubpd	0x20(%%rdx),%%ymm7,%%ymm7	\n\t"\
		"vsubpd	%%ymm4,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	%%ymm5,%%ymm1,%%ymm1		\n\t"\
		"vmovaps	%%ymm0,    (%%rbx)	\n\t"\
		"vmovaps	%%ymm1,0x20(%%rbx)	\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	%%ymm0,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	%%ymm1,%%ymm5,%%ymm5		\n\t"\
		"vmovaps	%%ymm4,    (%%rax)	\n\t"\
		"vmovaps	%%ymm5,0x20(%%rax)	\n\t"\
		"vsubpd	%%ymm7,%%ymm2,%%ymm2		\n\t"\
		"vsubpd	%%ymm6,%%ymm3,%%ymm3		\n\t"\
		"vmovaps	%%ymm2,    (%%rcx)	\n\t"\
		"vmovaps	%%ymm3,0x20(%%rdx)	\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	%%ymm2,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	%%ymm3,%%ymm6,%%ymm6		\n\t"\
		"vmovaps	%%ymm7,    (%%rdx)	\n\t"\
		"vmovaps	%%ymm6,0x20(%%rcx)	\n\t"\
		/* SSE2_RADIX4_DIF_IN_PLACE(r7 , r23, r15, r31): */\
		"leaq	%c[__i2](%%rax),%%rax	\n\t"/* All addresses += 2*ostride */\
		"leaq	%c[__i2](%%rbx),%%rbx	\n\t"\
		"leaq	%c[__i2](%%rcx),%%rcx	\n\t"\
		"leaq	%c[__i2](%%rdx),%%rdx	\n\t"\
		"vmovaps	    (%%rax),%%ymm0	\n\t"\
		"vmovaps	0x20(%%rax),%%ymm1	\n\t"\
		"vmovaps	    (%%rax),%%ymm2	\n\t"\
		"vmovaps	0x20(%%rax),%%ymm3	\n\t"\
		"vaddpd	    (%%rbx),%%ymm0,%%ymm0	\n\t"\
		"vaddpd	0x20(%%rbx),%%ymm1,%%ymm1	\n\t"\
		"vsubpd	    (%%rbx),%%ymm2,%%ymm2	\n\t"\
		"vsubpd	0x20(%%rbx),%%ymm3,%%ymm3	\n\t"\
		"vmovaps	    (%%rcx),%%ymm4	\n\t"\
		"vmovaps	0x20(%%rcx),%%ymm5	\n\t"\
		"vmovaps	    (%%rcx),%%ymm6	\n\t"\
		"vmovaps	0x20(%%rcx),%%ymm7	\n\t"\
		"vaddpd	    (%%rdx),%%ymm4,%%ymm4	\n\t"\
		"vaddpd	0x20(%%rdx),%%ymm5,%%ymm5	\n\t"\
		"vsubpd	    (%%rdx),%%ymm6,%%ymm6	\n\t"\
		"vsubpd	0x20(%%rdx),%%ymm7,%%ymm7	\n\t"\
		"vsubpd	%%ymm4,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	%%ymm5,%%ymm1,%%ymm1		\n\t"\
		"vmovaps	%%ymm0,    (%%rbx)	\n\t"\
		"vmovaps	%%ymm1,0x20(%%rbx)	\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	%%ymm0,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	%%ymm1,%%ymm5,%%ymm5		\n\t"\
		"vmovaps	%%ymm4,    (%%rax)	\n\t"\
		"vmovaps	%%ymm5,0x20(%%rax)	\n\t"\
		"vsubpd	%%ymm7,%%ymm2,%%ymm2		\n\t"\
		"vsubpd	%%ymm6,%%ymm3,%%ymm3		\n\t"\
		"vmovaps	%%ymm2,    (%%rcx)	\n\t"\
		"vmovaps	%%ymm3,0x20(%%rdx)	\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	%%ymm2,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	%%ymm3,%%ymm6,%%ymm6		\n\t"\
		"vmovaps	%%ymm7,    (%%rdx)	\n\t"\
		"vmovaps	%%ymm6,0x20(%%rcx)	\n\t"\
	/*** Now do 4 DFTs with internal twiddles on the 1*stride - separated data. Do blocks in order 0,2,1,3 to allow increment-only of rsi-datum from 1 block to the next: ***/\
		/* Block 0: r0-3 */\
		"movq	%[__in0],%%rsi	\n\t"\
		"movq	%[__out0],%%rax		\n\t"\
		"movq	%[__out1],%%rbx		\n\t"\
		"movq	%[__out2],%%rcx		\n\t"\
		"movq	%[__out3],%%rdx		\n\t"\
		"leaq	0x20(%%rsi),%%rdi	\n\t"/* Need separate address Im parts of outputs due to literal-offsets below */\
		"vmovaps	        (%%rsi),%%ymm0	\n\t"\
		"vmovaps	        (%%rdi),%%ymm1	\n\t"\
		"vmovaps	%c[__i2](%%rsi),%%ymm2	\n\t"\
		"vmovaps	%c[__i2](%%rdi),%%ymm3	\n\t"\
		"vsubpd	%c[__i2](%%rsi),%%ymm0,%%ymm0	\n\t"\
		"vsubpd	%c[__i2](%%rdi),%%ymm1,%%ymm1	\n\t"\
		"vaddpd	        (%%rsi),%%ymm2,%%ymm2	\n\t"\
		"vaddpd	        (%%rdi),%%ymm3,%%ymm3	\n\t"\
		"vmovaps	%c[__i1](%%rsi),%%ymm4	\n\t"\
		"vmovaps	%c[__i1](%%rdi),%%ymm5	\n\t"\
		"vmovaps	%c[__i3](%%rsi),%%ymm6	\n\t"\
		"vmovaps	%c[__i3](%%rdi),%%ymm7	\n\t"\
		"vsubpd	%c[__i3](%%rsi),%%ymm4,%%ymm4	\n\t"\
		"vsubpd	%c[__i3](%%rdi),%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%c[__i1](%%rsi),%%ymm6,%%ymm6	\n\t"\
		"vaddpd	%c[__i1](%%rdi),%%ymm7,%%ymm7	\n\t"\
		"vsubpd	%%ymm6,%%ymm2,%%ymm2		\n\t"\
		"vsubpd	%%ymm7,%%ymm3,%%ymm3		\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vmovaps	%%ymm2,    (%%rbx)	\n\t"\
		"vmovaps	%%ymm3,0x20(%%rbx)	\n\t"\
		"vaddpd	%%ymm2,%%ymm6,%%ymm6	\n\t"\
		"vaddpd	%%ymm3,%%ymm7,%%ymm7	\n\t"\
		"vmovaps	%%ymm6,    (%%rax)	\n\t"\
		"vmovaps	%%ymm7,0x20(%%rax)	\n\t"\
		"vsubpd	%%ymm5,%%ymm0,%%ymm0	\n\t"\
		"vsubpd	%%ymm4,%%ymm1,%%ymm1	\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4	\n\t"\
		"vmovaps	%%ymm0,    (%%rcx)	\n\t"\
		"vmovaps	%%ymm1,0x20(%%rdx)	\n\t"\
		"vaddpd	%%ymm0,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm1,%%ymm4,%%ymm4	\n\t"\
		"vmovaps	%%ymm5,    (%%rdx)	\n\t"\
		"vmovaps	%%ymm4,0x20(%%rcx)	\n\t"\
		/* Block 2: */\
		"movq	%[__out8],%%rax		\n\t	leaq	%c[__i4](%%rsi),%%rsi	\n\t"/* __in0 + 4*ostride */\
		"movq	%[__out9],%%rbx		\n\t"\
		"movq	%[__outa],%%rcx		\n\t"\
		"movq	%[__outb],%%rdx		\n\t"\
		"leaq	0x20(%%rsi),%%rdi	\n\t"\
		"vmovaps	%c[__i1](%%rsi),%%ymm4	\n\t"\
		"vmovaps	%c[__i3](%%rsi),%%ymm6	\n\t"\
		"vmovaps	%c[__i1](%%rdi),%%ymm5	\n\t"\
		"vmovaps	%c[__i3](%%rdi),%%ymm7	\n\t"\
		"vmovaps	%%ymm4,%%ymm0		\n\t"\
		"vmovaps	%%ymm6,%%ymm2		\n\t"\
		"vmovaps	%%ymm5,%%ymm1		\n\t"\
		"vmovaps	%%ymm7,%%ymm3		\n\t"\
		"movq	%[__isrt2],%%rdi	\n\t"\
		"addq	$0x20,%%rdi	\n\t"/* cc0 */\
		"vmulpd	    (%%rdi),%%ymm4,%%ymm4	\n\t"\
		"vmulpd	0x20(%%rdi),%%ymm6,%%ymm6	\n\t"\
		"vmulpd	0x20(%%rdi),%%ymm1,%%ymm1	\n\t"\
		"vmulpd	    (%%rdi),%%ymm3,%%ymm3	\n\t"\
		"vmulpd	    (%%rdi),%%ymm5,%%ymm5	\n\t"\
		"vmulpd	0x20(%%rdi),%%ymm7,%%ymm7	\n\t"\
		"vmulpd	0x20(%%rdi),%%ymm0,%%ymm0	\n\t"\
		"vmulpd	    (%%rdi),%%ymm2,%%ymm2	\n\t"\
		"vsubpd	%%ymm1,%%ymm4,%%ymm4		\n\t"\
		"vsubpd	%%ymm3,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	%%ymm0,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	%%ymm2,%%ymm7,%%ymm7		\n\t"\
		"vsubpd	%%ymm6,%%ymm4,%%ymm4		\n\t"\
		"vsubpd	%%ymm7,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	%%ymm4,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	%%ymm5,%%ymm7,%%ymm7		\n\t"\
		"leaq	0x20(%%rsi),%%rdi	\n\t"\
		"vmovaps	%c[__i2](%%rsi),%%ymm2	\n\t"\
		"vmovaps	%c[__i2](%%rdi),%%ymm3	\n\t"\
		"vsubpd	%c[__i2](%%rdi),%%ymm2,%%ymm2	\n\t"\
		"vaddpd	%c[__i2](%%rsi),%%ymm3,%%ymm3	\n\t"\
		"movq	%[__isrt2],%%rdi	\n\t"\
		"vmulpd	(%%rdi),%%ymm2,%%ymm2	\n\t"/* mul by isrt2 */\
		"vmulpd	(%%rdi),%%ymm3,%%ymm3	\n\t"\
		"leaq	0x20(%%rsi),%%rdi	\n\t"\
		"vmovaps	        (%%rsi),%%ymm0	\n\t"\
		"vmovaps	        (%%rdi),%%ymm1	\n\t"\
		"vsubpd	%%ymm2,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	%%ymm3,%%ymm1,%%ymm1		\n\t"\
		"vaddpd	        (%%rsi),%%ymm2,%%ymm2	\n\t"\
		"vaddpd	        (%%rdi),%%ymm3,%%ymm3	\n\t"\
		"vsubpd	%%ymm6,%%ymm2,%%ymm2		\n\t"\
		"vsubpd	%%ymm7,%%ymm3,%%ymm3		\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vmovaps	%%ymm2,    (%%rbx)	\n\t"\
		"vmovaps	%%ymm3,0x20(%%rbx)	\n\t"\
		"vaddpd	%%ymm2,%%ymm6,%%ymm6	\n\t"\
		"vaddpd	%%ymm3,%%ymm7,%%ymm7	\n\t"\
		"vmovaps	%%ymm6,    (%%rax)	\n\t"\
		"vmovaps	%%ymm7,0x20(%%rax)	\n\t"\
		"vsubpd	%%ymm5,%%ymm0,%%ymm0	\n\t"\
		"vsubpd	%%ymm4,%%ymm1,%%ymm1	\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4	\n\t"\
		"vmovaps	%%ymm0,    (%%rcx)	\n\t"\
		"vmovaps	%%ymm1,0x20(%%rdx)	\n\t"\
		"vaddpd	%%ymm0,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm1,%%ymm4,%%ymm4	\n\t"\
		"vmovaps	%%ymm5,    (%%rdx)	\n\t"\
		"vmovaps	%%ymm4,0x20(%%rcx)	\n\t"\
		/* Block 1: r8-b */\
		"movq	%[__out4],%%rax		\n\t	leaq	%c[__i4](%%rsi),%%rsi	\n\t"/* __in0 + 8*ostride */\
		"movq	%[__out5],%%rbx		\n\t"\
		"movq	%[__out6],%%rcx		\n\t"\
		"movq	%[__out7],%%rdx		\n\t"\
		"leaq	0x20(%%rsi),%%rdi	\n\t"\
		"vmovaps	        (%%rsi),%%ymm0	\n\t"\
		"vmovaps	        (%%rdi),%%ymm1	\n\t"\
		"vmovaps	%c[__i2](%%rsi),%%ymm2	\n\t"\
		"vmovaps	%c[__i2](%%rdi),%%ymm3	\n\t"\
		"vsubpd	%c[__i2](%%rdi),%%ymm0,%%ymm0	\n\t"\
		"vsubpd	%c[__i2](%%rsi),%%ymm1,%%ymm1	\n\t"\
		"vaddpd	        (%%rdi),%%ymm2,%%ymm2	\n\t"\
		"vaddpd	        (%%rsi),%%ymm3,%%ymm3	\n\t"\
		"vmovaps	%c[__i1](%%rsi),%%ymm4	\n\t"\
		"vmovaps	%c[__i1](%%rdi),%%ymm5	\n\t"\
		"vmovaps	%c[__i3](%%rsi),%%ymm6	\n\t"\
		"vmovaps	%c[__i3](%%rdi),%%ymm7	\n\t"\
		"vsubpd	%c[__i1](%%rdi),%%ymm4,%%ymm4	\n\t"\
		"vaddpd	%c[__i1](%%rsi),%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%c[__i3](%%rdi),%%ymm6,%%ymm6	\n\t"\
		"vsubpd	%c[__i3](%%rsi),%%ymm7,%%ymm7	\n\t"\
		"movq	%[__isrt2],%%rdi	\n\t"\
		"vmulpd	(%%rdi),%%ymm4,%%ymm4		\n\t"\
		"vmulpd	(%%rdi),%%ymm5,%%ymm5		\n\t"\
		"vmulpd	(%%rdi),%%ymm6,%%ymm6		\n\t"\
		"vmulpd	(%%rdi),%%ymm7,%%ymm7		\n\t"\
		"vsubpd	%%ymm6,%%ymm4,%%ymm4		\n\t"\
		"vsubpd	%%ymm7,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	%%ymm4,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	%%ymm5,%%ymm7,%%ymm7		\n\t"\
		"vsubpd	%%ymm4,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	%%ymm5,%%ymm2,%%ymm2		\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vmovaps	%%ymm0,    (%%rbx)	\n\t"\
		"vmovaps	%%ymm2,0x20(%%rbx)	\n\t"\
		"vaddpd	%%ymm0,%%ymm4,%%ymm4	\n\t"\
		"vaddpd	%%ymm2,%%ymm5,%%ymm5	\n\t"\
		"vmovaps	%%ymm4,    (%%rax)	\n\t"\
		"vmovaps	%%ymm5,0x20(%%rax)	\n\t"\
		"vsubpd	%%ymm7,%%ymm3,%%ymm3	\n\t"\
		"vsubpd	%%ymm6,%%ymm1,%%ymm1	\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7	\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6	\n\t"\
		"vmovaps	%%ymm3,    (%%rcx)	\n\t"\
		"vmovaps	%%ymm1,0x20(%%rdx)	\n\t"\
		"vaddpd	%%ymm3,%%ymm7,%%ymm7	\n\t"\
		"vaddpd	%%ymm1,%%ymm6,%%ymm6	\n\t"\
		"vmovaps	%%ymm7,    (%%rdx)	\n\t"\
		"vmovaps	%%ymm6,0x20(%%rcx)	\n\t"\
		/* Block 3: */\
		"movq	%[__outc],%%rax		\n\t	leaq	%c[__i4](%%rsi),%%rsi	\n\t"/* __in0 + c*ostride */\
		"movq	%[__outd],%%rbx		\n\t"\
		"movq	%[__oute],%%rcx		\n\t"\
		"movq	%[__outf],%%rdx		\n\t"\
		"leaq	0x20(%%rsi),%%rdi	\n\t"\
		"vmovaps	%c[__i1](%%rsi),%%ymm4	\n\t"\
		"vmovaps	%c[__i3](%%rsi),%%ymm6	\n\t"\
		"vmovaps	%c[__i1](%%rdi),%%ymm5	\n\t"\
		"vmovaps	%c[__i3](%%rdi),%%ymm7	\n\t"\
		"vmovaps	%%ymm4,%%ymm0		\n\t"\
		"vmovaps	%%ymm6,%%ymm2		\n\t"\
		"vmovaps	%%ymm5,%%ymm1		\n\t"\
		"vmovaps	%%ymm7,%%ymm3		\n\t"\
		"movq	%[__isrt2],%%rdi	\n\t"\
		"addq	$0x20,%%rdi	\n\t"/* cc0 */\
		"vmulpd	0x20(%%rdi),%%ymm4,%%ymm4	\n\t"\
		"vmulpd	    (%%rdi),%%ymm6,%%ymm6	\n\t"\
		"vmulpd	    (%%rdi),%%ymm1,%%ymm1	\n\t"\
		"vmulpd	0x20(%%rdi),%%ymm3,%%ymm3	\n\t"\
		"vmulpd	0x20(%%rdi),%%ymm5,%%ymm5	\n\t"\
		"vmulpd	    (%%rdi),%%ymm7,%%ymm7	\n\t"\
		"vmulpd	    (%%rdi),%%ymm0,%%ymm0	\n\t"\
		"vmulpd	0x20(%%rdi),%%ymm2,%%ymm2	\n\t"\
		"vsubpd	%%ymm1,%%ymm4,%%ymm4		\n\t"\
		"vsubpd	%%ymm3,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	%%ymm0,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	%%ymm2,%%ymm7,%%ymm7		\n\t"\
		"vsubpd	%%ymm6,%%ymm4,%%ymm4		\n\t"\
		"vsubpd	%%ymm7,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	%%ymm4,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	%%ymm5,%%ymm7,%%ymm7		\n\t"\
		"leaq	0x20(%%rsi),%%rdi	\n\t"\
		"vmovaps	%c[__i2](%%rsi),%%ymm2	\n\t"\
		"vmovaps	%c[__i2](%%rdi),%%ymm3	\n\t"\
		"vaddpd	%c[__i2](%%rdi),%%ymm2,%%ymm2	\n\t"\
		"vsubpd	%c[__i2](%%rsi),%%ymm3,%%ymm3	\n\t"\
		"movq	%[__isrt2],%%rdi	\n\t"\
		"vmulpd	(%%rdi),%%ymm2,%%ymm2	\n\t"/* mul by isrt2 */\
		"vmulpd	(%%rdi),%%ymm3,%%ymm3	\n\t"\
		"leaq	0x20(%%rsi),%%rdi	\n\t"\
		"vmovaps	        (%%rsi),%%ymm0	\n\t"\
		"vmovaps	        (%%rdi),%%ymm1	\n\t"\
		"vsubpd	%%ymm2,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	%%ymm3,%%ymm1,%%ymm1		\n\t"\
		"vaddpd	        (%%rsi),%%ymm2,%%ymm2	\n\t"\
		"vaddpd	        (%%rdi),%%ymm3,%%ymm3	\n\t"\
		"vsubpd	%%ymm4,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	%%ymm5,%%ymm1,%%ymm1		\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vmovaps	%%ymm0,    (%%rbx)	\n\t"\
		"vmovaps	%%ymm1,0x20(%%rbx)	\n\t"\
		"vaddpd	%%ymm0,%%ymm4,%%ymm4	\n\t"\
		"vaddpd	%%ymm1,%%ymm5,%%ymm5	\n\t"\
		"vmovaps	%%ymm4,    (%%rax)	\n\t"\
		"vmovaps	%%ymm5,0x20(%%rax)	\n\t"\
		"vsubpd	%%ymm7,%%ymm2,%%ymm2	\n\t"\
		"vsubpd	%%ymm6,%%ymm3,%%ymm3	\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7	\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6	\n\t"\
		"vmovaps	%%ymm2,    (%%rcx)	\n\t"\
		"vmovaps	%%ymm3,0x20(%%rdx)	\n\t"\
		"vaddpd	%%ymm2,%%ymm7,%%ymm7	\n\t"\
		"vaddpd	%%ymm3,%%ymm6,%%ymm6	\n\t"\
		"vmovaps	%%ymm7,    (%%rdx)	\n\t"\
		"vmovaps	%%ymm6,0x20(%%rcx)	\n\t"\
		:					/* outputs: none */\
		:[__in0] "m" (Xin0)	/* All inputs from memory addresses here */\
		,[__i1] "e" (Xi1)\
		,[__i2] "e" (Xi2)\
		,[__i3] "e" (Xi3)\
		,[__i4] "e" (Xi4)\
		,[__isrt2] "m" (Xisrt2)\
		,[__out0] "m" (Xout0)\
		,[__out1] "m" (Xout1)\
		,[__out2] "m" (Xout2)\
		,[__out3] "m" (Xout3)\
		,[__out4] "m" (Xout4)\
		,[__out5] "m" (Xout5)\
		,[__out6] "m" (Xout6)\
		,[__out7] "m" (Xout7)\
		,[__out8] "m" (Xout8)\
		,[__out9] "m" (Xout9)\
		,[__outa] "m" (Xouta)\
		,[__outb] "m" (Xoutb)\
		,[__outc] "m" (Xoutc)\
		,[__outd] "m" (Xoutd)\
		,[__oute] "m" (Xoute)\
		,[__outf] "m" (Xoutf)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
	);\
	}

#elif defined(USE_SSE2)

	// Based on the SSE2_RADIX16_DIT_NOTWIDDLE macro in radix16_ditN_cy_dif1_gcc64.h, but with completely
	// specifiable 16-input addressing required for usage as the power-of-2 component of a twiddleless
	// radix = [odd*2^n] DFT routine.
	// We use just a single output base-pointer plus literal ostrides which are [1,2,3,4]-multiples of
	// __01; this allows us to cut GP-register usage, which is absolutely a must for the 32-bit version
	// of the macro, and is a benefit to the 64-bit versions which code-fold to yield 2 side-by-side
	// streams of independently executable instructions, one for data in xmm0-7, the other using xmm8-15.
	#define SSE2_RADIX16_DIT_0TWIDDLE(Xin0,Xin1,Xin2,Xin3,Xin4,Xin5,Xin6,Xin7,Xin8,Xin9,Xina,Xinb,Xinc,Xind,Xine,Xinf, Xisrt2, Xout0,Xo1,Xo2,Xo3,Xo4)\
	{\
	__asm__ volatile (\
		"movq	%[__in0],%%rax		\n\t"\
		"movq	%[__in1],%%rbx		\n\t"\
		"movq	%[__in2],%%rcx		\n\t"\
		"movq	%[__in3],%%rdx		\n\t"\
		"/* SSE2_RADIX4_DIT_0TWIDDLE_B(r0 ): */\n\t"\
		"movaps	    (%%rax),%%xmm2	\n\t"\
		"movaps	    (%%rcx),%%xmm6	\n\t"\
		"movaps	0x10(%%rax),%%xmm3	\n\t"\
		"movaps	0x10(%%rcx),%%xmm7	\n\t"\
		"movaps	    (%%rbx),%%xmm0	\n\t"\
		"movaps	    (%%rdx),%%xmm4	\n\t"\
		"movaps	0x10(%%rbx),%%xmm1	\n\t"\
		"movaps	0x10(%%rdx),%%xmm5	\n\t"\
		"subpd	%%xmm0,%%xmm2		\n\t"\
		"subpd	%%xmm4,%%xmm6		\n\t"\
		"subpd	%%xmm1,%%xmm3		\n\t"\
		"subpd	%%xmm5,%%xmm7		\n\t"\
		"addpd	%%xmm0,%%xmm0		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm1,%%xmm1		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm2,%%xmm0		\n\t"\
		"addpd	%%xmm6,%%xmm4		\n\t"\
		"addpd	%%xmm3,%%xmm1		\n\t"\
		"addpd	%%xmm7,%%xmm5		\n\t"\
		"movq	%[__out0],%%rsi		\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm2		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"subpd	%%xmm6,%%xmm3		\n\t"\
		"leaq	0x10(%%rsi),%%rdi	\n\t"/* Need separate address Im parts of outputs due to literal-offsets below */\
		"movaps	%%xmm0,%c[__o2](%%rsi)	\n\t"\
		"movaps	%%xmm2,%c[__o3](%%rsi)	\n\t"\
		"movaps	%%xmm1,%c[__o2](%%rdi)	\n\t"\
		"movaps	%%xmm3,%c[__o1](%%rdi)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"addpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	%%xmm4,        (%%rsi)	\n\t"\
		"movaps	%%xmm7,%c[__o1](%%rsi)	\n\t"\
		"movaps	%%xmm5,        (%%rdi)	\n\t"\
		"movaps	%%xmm6,%c[__o3](%%rdi)	\n\t"\
		"/* SSE2_RADIX4_DIT_0TWIDDLE_B(r8 ): */\n\t"\
		"movq	%[__in4],%%rax		\n\t	leaq	%c[__o4](%%rsi),%%rsi	\n\t"/* __out0 + 4*ostride */\
		"movq	%[__in5],%%rbx		\n\t"\
		"movq	%[__in6],%%rcx		\n\t"\
		"movq	%[__in7],%%rdx		\n\t"\
		"movaps	    (%%rax),%%xmm2	\n\t"\
		"movaps	    (%%rcx),%%xmm6	\n\t"\
		"movaps	0x10(%%rax),%%xmm3	\n\t"\
		"movaps	0x10(%%rcx),%%xmm7	\n\t"\
		"movaps	    (%%rbx),%%xmm0	\n\t"\
		"movaps	    (%%rdx),%%xmm4	\n\t"\
		"movaps	0x10(%%rbx),%%xmm1	\n\t"\
		"movaps	0x10(%%rdx),%%xmm5	\n\t"\
		"subpd	%%xmm0,%%xmm2		\n\t"\
		"subpd	%%xmm4,%%xmm6		\n\t"\
		"subpd	%%xmm1,%%xmm3		\n\t"\
		"subpd	%%xmm5,%%xmm7		\n\t"\
		"addpd	%%xmm0,%%xmm0		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm1,%%xmm1		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm2,%%xmm0		\n\t"\
		"addpd	%%xmm6,%%xmm4		\n\t"\
		"addpd	%%xmm3,%%xmm1		\n\t"\
		"addpd	%%xmm7,%%xmm5		\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm2		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"subpd	%%xmm6,%%xmm3		\n\t"\
		"leaq	0x10(%%rsi),%%rdi	\n\t"\
		"movaps	%%xmm0,%c[__o2](%%rsi)	\n\t"\
		"movaps	%%xmm2,%c[__o3](%%rsi)	\n\t"\
		"movaps	%%xmm1,%c[__o2](%%rdi)	\n\t"\
		"movaps	%%xmm3,%c[__o1](%%rdi)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"addpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	%%xmm4,        (%%rsi)	\n\t"\
		"movaps	%%xmm7,%c[__o1](%%rsi)	\n\t"\
		"movaps	%%xmm5,        (%%rdi)	\n\t"\
		"movaps	%%xmm6,%c[__o3](%%rdi)	\n\t"\
		"/* SSE2_RADIX4_DIT_0TWIDDLE_B(r16): */\n\t"\
		"movq	%[__in8],%%rax			\n\t	leaq	%c[__o4](%%rsi),%%rsi	\n\t"/* __out0 + 8*ostride */\
		"movq	%[__in9],%%rbx			\n\t"\
		"movq	%[__ina],%%rcx			\n\t"\
		"movq	%[__inb],%%rdx			\n\t"\
		"movaps	    (%%rax),%%xmm2	\n\t"\
		"movaps	    (%%rcx),%%xmm6	\n\t"\
		"movaps	0x10(%%rax),%%xmm3	\n\t"\
		"movaps	0x10(%%rcx),%%xmm7	\n\t"\
		"movaps	    (%%rbx),%%xmm0	\n\t"\
		"movaps	    (%%rdx),%%xmm4	\n\t"\
		"movaps	0x10(%%rbx),%%xmm1	\n\t"\
		"movaps	0x10(%%rdx),%%xmm5	\n\t"\
		"subpd	%%xmm0,%%xmm2		\n\t"\
		"subpd	%%xmm4,%%xmm6		\n\t"\
		"subpd	%%xmm1,%%xmm3		\n\t"\
		"subpd	%%xmm5,%%xmm7		\n\t"\
		"addpd	%%xmm0,%%xmm0		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm1,%%xmm1		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm2,%%xmm0		\n\t"\
		"addpd	%%xmm6,%%xmm4		\n\t"\
		"addpd	%%xmm3,%%xmm1		\n\t"\
		"addpd	%%xmm7,%%xmm5		\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm2		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"subpd	%%xmm6,%%xmm3		\n\t"\
		"leaq	0x10(%%rsi),%%rdi	\n\t"\
		"movaps	%%xmm0,%c[__o2](%%rsi)	\n\t"\
		"movaps	%%xmm2,%c[__o3](%%rsi)	\n\t"\
		"movaps	%%xmm1,%c[__o2](%%rdi)	\n\t"\
		"movaps	%%xmm3,%c[__o1](%%rdi)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"addpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	%%xmm4,        (%%rsi)	\n\t"\
		"movaps	%%xmm7,%c[__o1](%%rsi)	\n\t"\
		"movaps	%%xmm5,        (%%rdi)	\n\t"\
		"movaps	%%xmm6,%c[__o3](%%rdi)	\n\t"\
		"/* SSE2_RADIX4_DIT_0TWIDDLE_B(r24): */\n\t"\
		"movq	%[__inc],%%rax			\n\t	leaq	%c[__o4](%%rsi),%%rsi	\n\t"/* __out0 + c*ostride */\
		"movq	%[__ind],%%rbx			\n\t"\
		"movq	%[__ine],%%rcx			\n\t"\
		"movq	%[__inf],%%rdx			\n\t"\
		"movaps	    (%%rax),%%xmm2	\n\t"\
		"movaps	    (%%rcx),%%xmm6	\n\t"\
		"movaps	0x10(%%rax),%%xmm3	\n\t"\
		"movaps	0x10(%%rcx),%%xmm7	\n\t"\
		"movaps	    (%%rbx),%%xmm0	\n\t"\
		"movaps	    (%%rdx),%%xmm4	\n\t"\
		"movaps	0x10(%%rbx),%%xmm1	\n\t"\
		"movaps	0x10(%%rdx),%%xmm5	\n\t"\
		"subpd	%%xmm0,%%xmm2		\n\t"\
		"subpd	%%xmm4,%%xmm6		\n\t"\
		"subpd	%%xmm1,%%xmm3		\n\t"\
		"subpd	%%xmm5,%%xmm7		\n\t"\
		"addpd	%%xmm0,%%xmm0		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm1,%%xmm1		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm2,%%xmm0		\n\t"\
		"addpd	%%xmm6,%%xmm4		\n\t"\
		"addpd	%%xmm3,%%xmm1		\n\t"\
		"addpd	%%xmm7,%%xmm5		\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm2		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"subpd	%%xmm6,%%xmm3		\n\t"\
		"leaq	0x10(%%rsi),%%rdi	\n\t"\
		"movaps	%%xmm0,%c[__o2](%%rsi)	\n\t"\
		"movaps	%%xmm2,%c[__o3](%%rsi)	\n\t"\
		"movaps	%%xmm1,%c[__o2](%%rdi)	\n\t"\
		"movaps	%%xmm3,%c[__o1](%%rdi)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"addpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	%%xmm4,        (%%rsi)	\n\t"\
		"movaps	%%xmm7,%c[__o1](%%rsi)	\n\t"\
		"movaps	%%xmm5,        (%%rdi)	\n\t"\
		"movaps	%%xmm6,%c[__o3](%%rdi)	\n\t"\
	/*** Now do 4 DFTs with internal twiddles on the 4*stride - separated data: ***/\
		"movq	%[__out0],%%rax		\n\t"\
		"leaq	%c[__o4](%%rax),%%rbx	\n\t"/* __out0 +   [4*ostride] */\
		"leaq	%c[__o4](%%rbx),%%rcx	\n\t"/* __out0 + 2*[4*ostride] */\
		"leaq	%c[__o4](%%rcx),%%rdx	\n\t"/* __out0 + 3*[4*ostride] */\
		/* Block 0: */\
		"movaps	    (%%rax),%%xmm0	\n\t"\
		"movaps	0x10(%%rax),%%xmm1	\n\t"\
		"movaps	    (%%rbx),%%xmm2	\n\t"\
		"movaps	0x10(%%rbx),%%xmm3	\n\t"\
		"subpd	    (%%rbx),%%xmm0	\n\t"\
		"subpd	0x10(%%rbx),%%xmm1	\n\t"\
		"addpd	    (%%rax),%%xmm2	\n\t"\
		"addpd	0x10(%%rax),%%xmm3	\n\t"\
		"movaps	    (%%rcx),%%xmm4	\n\t"\
		"movaps	0x10(%%rcx),%%xmm5	\n\t"\
		"movaps	    (%%rdx),%%xmm6	\n\t"\
		"movaps	0x10(%%rdx),%%xmm7	\n\t"\
		"subpd	    (%%rdx),%%xmm4	\n\t"\
		"subpd	0x10(%%rdx),%%xmm5	\n\t"\
		"addpd	    (%%rcx),%%xmm6	\n\t"\
		"addpd	0x10(%%rcx),%%xmm7	\n\t"\
		"							\n\t"\
		"subpd	%%xmm6,%%xmm2		\n\t"\
		"subpd	%%xmm7,%%xmm3		\n\t"\
		"movaps	%%xmm2,    (%%rcx)	\n\t"\
		"movaps	%%xmm3,0x10(%%rcx)	\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm2,%%xmm6		\n\t"\
		"addpd	%%xmm3,%%xmm7		\n\t"\
		"movaps	%%xmm6,    (%%rax)	\n\t"\
		"movaps	%%xmm7,0x10(%%rax)	\n\t"\
		"							\n\t"\
		"subpd	%%xmm5,%%xmm0		\n\t"\
		"subpd	%%xmm4,%%xmm1		\n\t"\
		"movaps	%%xmm0,    (%%rdx)	\n\t"\
		"movaps	%%xmm1,0x10(%%rbx)	\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm0,%%xmm5		\n\t"\
		"addpd	%%xmm1,%%xmm4		\n\t"\
		"movaps	%%xmm5,    (%%rbx)	\n\t"\
		"movaps	%%xmm4,0x10(%%rdx)	\n\t"\
		"							\n\t"\
		/* Block 2: */\
		"leaq	%c[__o2](%%rax),%%rax	\n\t"/* All addresses += 2*ostride */\
		"leaq	%c[__o2](%%rbx),%%rbx	\n\t"\
		"leaq	%c[__o2](%%rcx),%%rcx	\n\t"\
		"leaq	%c[__o2](%%rdx),%%rdx	\n\t"\
		"movq	%[__isrt2],%%rdi	\n\t"\
		"movaps	(%%rdi),%%xmm2		\n\t"\
		"							\n\t"\
		"movaps	    (%%rcx),%%xmm4	\n\t"\
		"movaps	0x10(%%rcx),%%xmm5	\n\t"\
		"movaps	    (%%rdx),%%xmm0	\n\t"\
		"movaps	0x10(%%rdx),%%xmm1	\n\t"\
		"						\n\t"\
		"addpd	0x10(%%rcx),%%xmm4	\n\t"\
		"subpd	    (%%rcx),%%xmm5	\n\t"\
		"subpd	0x10(%%rdx),%%xmm0	\n\t"\
		"addpd	    (%%rdx),%%xmm1	\n\t"\
		"mulpd	%%xmm2,%%xmm4		\n\t"\
		"mulpd	%%xmm2,%%xmm5		\n\t"\
		"mulpd	%%xmm2,%%xmm0		\n\t"\
		"mulpd	%%xmm2,%%xmm1		\n\t"\
		"movaps	%%xmm4,%%xmm6		\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t"\
		"							\n\t"\
		"subpd	%%xmm0,%%xmm4		\n\t"\
		"subpd	%%xmm1,%%xmm5		\n\t"\
		"addpd	%%xmm0,%%xmm6		\n\t"\
		"addpd	%%xmm1,%%xmm7		\n\t"\
		"							\n\t"\
		"movaps	    (%%rax),%%xmm0	\n\t"\
		"movaps	0x10(%%rax),%%xmm1	\n\t"\
		"movaps	    (%%rbx),%%xmm2	\n\t"\
		"movaps	0x10(%%rbx),%%xmm3	\n\t"\
		"							\n\t"\
		"subpd	0x10(%%rbx),%%xmm0	\n\t"\
		"subpd	    (%%rbx),%%xmm1	\n\t"\
		"addpd	    (%%rax),%%xmm3	\n\t"\
		"addpd	0x10(%%rax),%%xmm2	\n\t"\
		"							\n\t"\
		"subpd	%%xmm4,%%xmm3		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"movaps	%%xmm3,    (%%rcx)	\n\t"\
		"movaps	%%xmm1,0x10(%%rcx)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm3,%%xmm4		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"movaps	%%xmm4,    (%%rax)	\n\t"\
		"movaps	%%xmm5,0x10(%%rax)	\n\t"\
		"							\n\t"\
		"subpd	%%xmm7,%%xmm0		\n\t"\
		"subpd	%%xmm6,%%xmm2		\n\t"\
		"movaps	%%xmm0,    (%%rdx)	\n\t"\
		"movaps	%%xmm2,0x10(%%rbx)	\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm7		\n\t"\
		"addpd	%%xmm2,%%xmm6		\n\t"\
		"movaps	%%xmm7,    (%%rbx)	\n\t"\
		"movaps	%%xmm6,0x10(%%rdx)	\n\t"\
		"							\n\t"\
		/* Block 1: */\
		"leaq	-%c[__o1](%%rax),%%rax	\n\t"/* All addresses -= 1*ostride */\
		"leaq	-%c[__o1](%%rbx),%%rbx	\n\t"\
		"leaq	-%c[__o1](%%rcx),%%rcx	\n\t"\
		"leaq	-%c[__o1](%%rdx),%%rdx	\n\t"\
		"leaq	0x10(%%rdi),%%rsi	\n\t"/* cc0 */\
		"							\n\t"\
		"movaps	    (%%rdx),%%xmm0	\n\t"\
		"movaps	0x10(%%rdx),%%xmm1	\n\t"\
		"movaps	    (%%rdx),%%xmm2	\n\t"\
		"movaps	0x10(%%rdx),%%xmm3	\n\t"\
		"							\n\t"\
		"mulpd	0x10(%%rsi),%%xmm0	\n\t"\
		"mulpd	0x10(%%rsi),%%xmm1	\n\t"\
		"mulpd	    (%%rsi),%%xmm2	\n\t"\
		"mulpd	    (%%rsi),%%xmm3	\n\t"\
		"subpd	%%xmm2,%%xmm1		\n\t"\
		"addpd	%%xmm3,%%xmm0		\n\t"\
		"							\n\t"\
		"movaps	    (%%rcx),%%xmm4	\n\t"\
		"movaps	0x10(%%rcx),%%xmm5	\n\t"\
		"movaps	    (%%rcx),%%xmm6	\n\t"\
		"movaps	0x10(%%rcx),%%xmm7	\n\t"\
		"							\n\t"\
		"mulpd	    (%%rsi),%%xmm4	\n\t"\
		"mulpd	    (%%rsi),%%xmm5	\n\t"\
		"mulpd	0x10(%%rsi),%%xmm6	\n\t"\
		"mulpd	0x10(%%rsi),%%xmm7	\n\t"\
		"subpd	%%xmm6,%%xmm5		\n\t"\
		"addpd	%%xmm7,%%xmm4		\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t"\
		"movaps	%%xmm4,%%xmm6		\n\t"\
		"							\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"subpd	%%xmm0,%%xmm6		\n\t"\
		"subpd	%%xmm1,%%xmm7		\n\t"\
		"							\n\t"\
		"movaps	    (%%rbx),%%xmm2	\n\t"\
		"movaps	0x10(%%rbx),%%xmm3	\n\t"\
		"movaps	    (%%rax),%%xmm0	\n\t"\
		"movaps	0x10(%%rax),%%xmm1	\n\t"\
		"addpd	0x10(%%rbx),%%xmm2	\n\t"\
		"subpd	    (%%rbx),%%xmm3	\n\t"\
		"mulpd	    (%%rdi),%%xmm2	\n\t"\
		"mulpd	    (%%rdi),%%xmm3	\n\t"\
		"							\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"							\n\t"\
		"subpd	%%xmm4,%%xmm2		\n\t"\
		"subpd	%%xmm5,%%xmm3		\n\t"\
		"movaps	%%xmm2,    (%%rcx)	\n\t"\
		"movaps	%%xmm3,0x10(%%rcx)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm2,%%xmm4		\n\t"\
		"addpd	%%xmm3,%%xmm5		\n\t"\
		"movaps	%%xmm4,    (%%rax)	\n\t"\
		"movaps	%%xmm5,0x10(%%rax)	\n\t"\
		"							\n\t"\
		"subpd	%%xmm7,%%xmm0		\n\t"\
		"subpd	%%xmm6,%%xmm1		\n\t"\
		"movaps	%%xmm0,    (%%rdx)	\n\t"\
		"movaps	%%xmm1,0x10(%%rbx)	\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm7		\n\t"\
		"addpd	%%xmm1,%%xmm6		\n\t"\
		"movaps	%%xmm7,    (%%rbx)	\n\t"\
		"movaps	%%xmm6,0x10(%%rdx)	\n\t"\
		"							\n\t"\
		/* Block 3: */\
		"leaq	%c[__o2](%%rax),%%rax	\n\t"/* All addresses += 2*ostride */\
		"leaq	%c[__o2](%%rbx),%%rbx	\n\t"\
		"leaq	%c[__o2](%%rcx),%%rcx	\n\t"\
		"leaq	%c[__o2](%%rdx),%%rdx	\n\t"\
		"							\n\t"\
		"movaps	    (%%rdx),%%xmm0	\n\t"\
		"movaps	0x10(%%rdx),%%xmm1	\n\t"\
		"movaps	    (%%rdx),%%xmm2	\n\t"\
		"movaps	0x10(%%rdx),%%xmm3	\n\t"\
		"						\n\t"\
		"mulpd	    (%%rsi),%%xmm0	\n\t"\
		"mulpd	    (%%rsi),%%xmm1	\n\t"\
		"mulpd	0x10(%%rsi),%%xmm2	\n\t"\
		"mulpd	0x10(%%rsi),%%xmm3	\n\t"\
		"subpd	%%xmm2,%%xmm1		\n\t"\
		"addpd	%%xmm3,%%xmm0		\n\t"\
		"							\n\t"\
		"movaps	    (%%rcx),%%xmm4	\n\t"\
		"movaps	0x10(%%rcx),%%xmm5	\n\t"\
		"movaps	    (%%rcx),%%xmm6	\n\t"\
		"movaps	0x10(%%rcx),%%xmm7	\n\t"\
		"						\n\t"\
		"mulpd	0x10(%%rsi),%%xmm4	\n\t"\
		"mulpd	0x10(%%rsi),%%xmm5	\n\t"\
		"mulpd	    (%%rsi),%%xmm6	\n\t"\
		"mulpd	    (%%rsi),%%xmm7	\n\t"\
		"subpd	%%xmm6,%%xmm5		\n\t"\
		"addpd	%%xmm7,%%xmm4		\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t"\
		"movaps	%%xmm4,%%xmm6		\n\t"\
		"							\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"subpd	%%xmm0,%%xmm6		\n\t"\
		"subpd	%%xmm1,%%xmm7		\n\t"\
		"							\n\t"\
		"movaps	    (%%rbx),%%xmm2	\n\t"\
		"movaps	0x10(%%rbx),%%xmm3	\n\t"\
		"movaps	    (%%rax),%%xmm0	\n\t"\
		"movaps	0x10(%%rax),%%xmm1	\n\t"\
		"subpd	0x10(%%rbx),%%xmm2	\n\t"\
		"addpd	    (%%rbx),%%xmm3	\n\t"\
		"mulpd	    (%%rdi),%%xmm2	\n\t"\
		"mulpd	    (%%rdi),%%xmm3	\n\t"\
		"							\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"							\n\t"\
		"subpd	%%xmm6,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm1		\n\t"\
		"movaps	%%xmm0,    (%%rcx)	\n\t"\
		"movaps	%%xmm1,0x10(%%rcx)	\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm0,%%xmm6		\n\t"\
		"addpd	%%xmm1,%%xmm7		\n\t"\
		"movaps	%%xmm6,    (%%rax)	\n\t"\
		"movaps	%%xmm7,0x10(%%rax)	\n\t"\
		"							\n\t"\
		"subpd	%%xmm5,%%xmm2		\n\t"\
		"subpd	%%xmm4,%%xmm3		\n\t"\
		"movaps	%%xmm2,    (%%rdx)	\n\t"\
		"movaps	%%xmm3,0x10(%%rbx)	\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm2,%%xmm5		\n\t"\
		"addpd	%%xmm3,%%xmm4		\n\t"\
		"movaps	%%xmm5,    (%%rbx)	\n\t"\
		"movaps	%%xmm4,0x10(%%rdx)	\n\t"\
		:					/* outputs: none */\
		:[__in0] "m" (Xin0)	/* All inputs from memory addresses here */\
		,[__in1] "m" (Xin1)\
		,[__in2] "m" (Xin2)\
		,[__in3] "m" (Xin3)\
		,[__in4] "m" (Xin4)\
		,[__in5] "m" (Xin5)\
		,[__in6] "m" (Xin6)\
		,[__in7] "m" (Xin7)\
		,[__in8] "m" (Xin8)\
		,[__in9] "m" (Xin9)\
		,[__ina] "m" (Xina)\
		,[__inb] "m" (Xinb)\
		,[__inc] "m" (Xinc)\
		,[__ind] "m" (Xind)\
		,[__ine] "m" (Xine)\
		,[__inf] "m" (Xinf)\
		,[__isrt2] "m" (Xisrt2)\
		,[__out0] "m" (Xout0)\
		,[__o1] "e" (Xo1)\
		,[__o2] "e" (Xo2)\
		,[__o3] "e" (Xo3)\
		,[__o4] "e" (Xo4)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
	);\
	}

	// Based on the SSE2_RADIX16_DIF_NOTWIDDLE macro in radix16_ditN_cy_dif1_gcc64.h, but with completely
	// specifiable 16-output addressing required for usage as the power-of-2 component of a twiddleless
	// radix = [odd*2^n] DFT routine.
	#define SSE2_RADIX16_DIF_0TWIDDLE(Xin0,Xi1,Xi2,Xi3,Xi4, Xisrt2, Xout0,Xout1,Xout2,Xout3,Xout4,Xout5,Xout6,Xout7,Xout8,Xout9,Xouta,Xoutb,Xoutc,Xoutd,Xoute,Xoutf)\
	{\
	__asm__ volatile (\
		"/* SSE2_RADIX4_DIF_IN_PLACE(r1 , r17, r9 , r25): */\n\t"\
		"movq	%[__in0],%%rax		\n\t"\
		"leaq	%c[__i4](%%rax),%%rcx	\n\t"/* __in0 +   [4*istride]; note BR of [a,b,c,d]-ptrs, i.e. b/c swap */\
		"leaq	%c[__i4](%%rcx),%%rbx	\n\t"/* __in0 + 2*[4*istride] */\
		"leaq	%c[__i4](%%rbx),%%rdx	\n\t"/* __in0 + 3*[4*istride] */\
		"							\n\t"\
		"movaps	    (%%rax),%%xmm0	\n\t"\
		"movaps	0x10(%%rax),%%xmm1	\n\t"\
		"movaps	    (%%rax),%%xmm2	\n\t"\
		"movaps	0x10(%%rax),%%xmm3	\n\t"\
		"							\n\t"\
		"addpd	    (%%rbx),%%xmm0	\n\t"\
		"addpd	0x10(%%rbx),%%xmm1	\n\t"\
		"subpd	    (%%rbx),%%xmm2	\n\t"\
		"subpd	0x10(%%rbx),%%xmm3	\n\t"\
		"							\n\t"\
		"movaps	    (%%rcx),%%xmm4	\n\t"\
		"movaps	0x10(%%rcx),%%xmm5	\n\t"\
		"movaps	    (%%rcx),%%xmm6	\n\t"\
		"movaps	0x10(%%rcx),%%xmm7	\n\t"\
		"							\n\t"\
		"addpd	    (%%rdx),%%xmm4	\n\t"\
		"addpd	0x10(%%rdx),%%xmm5	\n\t"\
		"subpd	    (%%rdx),%%xmm6	\n\t"\
		"subpd	0x10(%%rdx),%%xmm7	\n\t"\
		"							\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"movaps	%%xmm0,    (%%rbx)	\n\t"\
		"movaps	%%xmm1,0x10(%%rbx)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"movaps	%%xmm4,    (%%rax)	\n\t"\
		"movaps	%%xmm5,0x10(%%rax)	\n\t"\
		"							\n\t"\
		"subpd	%%xmm7,%%xmm2		\n\t"\
		"subpd	%%xmm6,%%xmm3		\n\t"\
		"movaps	%%xmm2,    (%%rcx)	\n\t"\
		"movaps	%%xmm3,0x10(%%rdx)	\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	%%xmm7,    (%%rdx)	\n\t"\
		"movaps	%%xmm6,0x10(%%rcx)	\n\t"\
		"/* SSE2_RADIX4_DIF_IN_PLACE(r5 , r21, r13, r29): */\n\t"\
		"leaq	%c[__i2](%%rax),%%rax	\n\t"/* All addresses += 2*ostride */\
		"leaq	%c[__i2](%%rbx),%%rbx	\n\t"\
		"leaq	%c[__i2](%%rcx),%%rcx	\n\t"\
		"leaq	%c[__i2](%%rdx),%%rdx	\n\t"\
		"							\n\t"\
		"movaps	    (%%rax),%%xmm0	\n\t"\
		"movaps	0x10(%%rax),%%xmm1	\n\t"\
		"movaps	    (%%rax),%%xmm2	\n\t"\
		"movaps	0x10(%%rax),%%xmm3	\n\t"\
		"							\n\t"\
		"addpd	    (%%rbx),%%xmm0	\n\t"\
		"addpd	0x10(%%rbx),%%xmm1	\n\t"\
		"subpd	    (%%rbx),%%xmm2	\n\t"\
		"subpd	0x10(%%rbx),%%xmm3	\n\t"\
		"							\n\t"\
		"movaps	    (%%rcx),%%xmm4	\n\t"\
		"movaps	0x10(%%rcx),%%xmm5	\n\t"\
		"movaps	    (%%rcx),%%xmm6	\n\t"\
		"movaps	0x10(%%rcx),%%xmm7	\n\t"\
		"							\n\t"\
		"addpd	    (%%rdx),%%xmm4	\n\t"\
		"addpd	0x10(%%rdx),%%xmm5	\n\t"\
		"subpd	    (%%rdx),%%xmm6	\n\t"\
		"subpd	0x10(%%rdx),%%xmm7	\n\t"\
		"							\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"movaps	%%xmm0,    (%%rbx)	\n\t"\
		"movaps	%%xmm1,0x10(%%rbx)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"movaps	%%xmm4,    (%%rax)	\n\t"\
		"movaps	%%xmm5,0x10(%%rax)	\n\t"\
		"							\n\t"\
		"subpd	%%xmm7,%%xmm2		\n\t"\
		"subpd	%%xmm6,%%xmm3		\n\t"\
		"movaps	%%xmm2,    (%%rcx)	\n\t"\
		"movaps	%%xmm3,0x10(%%rdx)	\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	%%xmm7,    (%%rdx)	\n\t"\
		"movaps	%%xmm6,0x10(%%rcx)	\n\t"\
		"/* SSE2_RADIX4_DIF_IN_PLACE(r3 , r19, r11, r27): */\n\t"\
		"leaq	-%c[__i1](%%rax),%%rax	\n\t"/* All addresses -= 1*ostride */\
		"leaq	-%c[__i1](%%rbx),%%rbx	\n\t"\
		"leaq	-%c[__i1](%%rcx),%%rcx	\n\t"\
		"leaq	-%c[__i1](%%rdx),%%rdx	\n\t"\
		"							\n\t"\
		"movaps	    (%%rax),%%xmm0	\n\t"\
		"movaps	0x10(%%rax),%%xmm1	\n\t"\
		"movaps	    (%%rax),%%xmm2	\n\t"\
		"movaps	0x10(%%rax),%%xmm3	\n\t"\
		"							\n\t"\
		"addpd	    (%%rbx),%%xmm0	\n\t"\
		"addpd	0x10(%%rbx),%%xmm1	\n\t"\
		"subpd	    (%%rbx),%%xmm2	\n\t"\
		"subpd	0x10(%%rbx),%%xmm3	\n\t"\
		"							\n\t"\
		"movaps	    (%%rcx),%%xmm4	\n\t"\
		"movaps	0x10(%%rcx),%%xmm5	\n\t"\
		"movaps	    (%%rcx),%%xmm6	\n\t"\
		"movaps	0x10(%%rcx),%%xmm7	\n\t"\
		"							\n\t"\
		"addpd	    (%%rdx),%%xmm4	\n\t"\
		"addpd	0x10(%%rdx),%%xmm5	\n\t"\
		"subpd	    (%%rdx),%%xmm6	\n\t"\
		"subpd	0x10(%%rdx),%%xmm7	\n\t"\
		"							\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"movaps	%%xmm0,    (%%rbx)	\n\t"\
		"movaps	%%xmm1,0x10(%%rbx)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"movaps	%%xmm4,    (%%rax)	\n\t"\
		"movaps	%%xmm5,0x10(%%rax)	\n\t"\
		"							\n\t"\
		"subpd	%%xmm7,%%xmm2		\n\t"\
		"subpd	%%xmm6,%%xmm3		\n\t"\
		"movaps	%%xmm2,    (%%rcx)	\n\t"\
		"movaps	%%xmm3,0x10(%%rdx)	\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	%%xmm7,    (%%rdx)	\n\t"\
		"movaps	%%xmm6,0x10(%%rcx)	\n\t"\
		"/* SSE2_RADIX4_DIF_IN_PLACE(r7 , r23, r15, r31): */\n\t"\
		"leaq	%c[__i2](%%rax),%%rax	\n\t"/* All addresses += 2*ostride */\
		"leaq	%c[__i2](%%rbx),%%rbx	\n\t"\
		"leaq	%c[__i2](%%rcx),%%rcx	\n\t"\
		"leaq	%c[__i2](%%rdx),%%rdx	\n\t"\
		"							\n\t"\
		"movaps	    (%%rax),%%xmm0	\n\t"\
		"movaps	0x10(%%rax),%%xmm1	\n\t"\
		"movaps	    (%%rax),%%xmm2	\n\t"\
		"movaps	0x10(%%rax),%%xmm3	\n\t"\
		"							\n\t"\
		"addpd	    (%%rbx),%%xmm0	\n\t"\
		"addpd	0x10(%%rbx),%%xmm1	\n\t"\
		"subpd	    (%%rbx),%%xmm2	\n\t"\
		"subpd	0x10(%%rbx),%%xmm3	\n\t"\
		"							\n\t"\
		"movaps	    (%%rcx),%%xmm4	\n\t"\
		"movaps	0x10(%%rcx),%%xmm5	\n\t"\
		"movaps	    (%%rcx),%%xmm6	\n\t"\
		"movaps	0x10(%%rcx),%%xmm7	\n\t"\
		"							\n\t"\
		"addpd	    (%%rdx),%%xmm4	\n\t"\
		"addpd	0x10(%%rdx),%%xmm5	\n\t"\
		"subpd	    (%%rdx),%%xmm6	\n\t"\
		"subpd	0x10(%%rdx),%%xmm7	\n\t"\
		"							\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"movaps	%%xmm0,    (%%rbx)	\n\t"\
		"movaps	%%xmm1,0x10(%%rbx)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"movaps	%%xmm4,    (%%rax)	\n\t"\
		"movaps	%%xmm5,0x10(%%rax)	\n\t"\
		"							\n\t"\
		"subpd	%%xmm7,%%xmm2		\n\t"\
		"subpd	%%xmm6,%%xmm3		\n\t"\
		"movaps	%%xmm2,    (%%rcx)	\n\t"\
		"movaps	%%xmm3,0x10(%%rdx)	\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	%%xmm7,    (%%rdx)	\n\t"\
		"movaps	%%xmm6,0x10(%%rcx)	\n\t"\
		"							\n\t"\
	/*** Now do 4 DFTs with internal twiddles on the 1*stride - separated data. Do blocks in order 0,2,1,3 to allow increment-only of rsi-datum from 1 block to the next: ***/\
		/* Block 0: r0-3 */\
		"movq	%[__in0],%%rsi	\n\t"\
		"movq	%[__out0],%%rax		\n\t"\
		"movq	%[__out1],%%rbx		\n\t"\
		"movq	%[__out2],%%rcx		\n\t"\
		"movq	%[__out3],%%rdx		\n\t"\
		"							\n\t"\
		"leaq	0x10(%%rsi),%%rdi	\n\t"/* Need separate address Im parts of outputs due to literal-offsets below */\
		"movaps	        (%%rsi),%%xmm0	\n\t"\
		"movaps	        (%%rdi),%%xmm1	\n\t"\
		"movaps	%c[__i2](%%rsi),%%xmm2	\n\t"\
		"movaps	%c[__i2](%%rdi),%%xmm3	\n\t"\
		"						\n\t"\
		"subpd	%c[__i2](%%rsi),%%xmm0	\n\t"\
		"subpd	%c[__i2](%%rdi),%%xmm1	\n\t"\
		"addpd	        (%%rsi),%%xmm2	\n\t"\
		"addpd	        (%%rdi),%%xmm3	\n\t"\
		"						\n\t"\
		"movaps	%c[__i1](%%rsi),%%xmm4	\n\t"\
		"movaps	%c[__i1](%%rdi),%%xmm5	\n\t"\
		"movaps	%c[__i3](%%rsi),%%xmm6	\n\t"\
		"movaps	%c[__i3](%%rdi),%%xmm7	\n\t"\
		"						\n\t"\
		"subpd	%c[__i3](%%rsi),%%xmm4	\n\t"\
		"subpd	%c[__i3](%%rdi),%%xmm5	\n\t"\
		"addpd	%c[__i1](%%rsi),%%xmm6	\n\t"\
		"addpd	%c[__i1](%%rdi),%%xmm7	\n\t"\
		"							\n\t"\
		"subpd	%%xmm6,%%xmm2		\n\t"\
		"subpd	%%xmm7,%%xmm3		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"movaps	%%xmm2,    (%%rbx)	\n\t"\
		"movaps	%%xmm3,0x10(%%rbx)	\n\t"\
		"addpd	%%xmm2,		%%xmm6	\n\t"\
		"addpd	%%xmm3,		%%xmm7	\n\t"\
		"movaps	%%xmm6,    (%%rax)	\n\t"\
		"movaps	%%xmm7,0x10(%%rax)	\n\t"\
		"							\n\t"\
		"subpd	%%xmm5,		%%xmm0	\n\t"\
		"subpd	%%xmm4,		%%xmm1	\n\t"\
		"addpd	%%xmm5,		%%xmm5	\n\t"\
		"addpd	%%xmm4,		%%xmm4	\n\t"\
		"movaps	%%xmm0,    (%%rcx)	\n\t"\
		"movaps	%%xmm1,0x10(%%rdx)	\n\t"\
		"addpd	%%xmm0,		%%xmm5	\n\t"\
		"addpd	%%xmm1,		%%xmm4	\n\t"\
		"movaps	%%xmm5,    (%%rdx)	\n\t"\
		"movaps	%%xmm4,0x10(%%rcx)	\n\t"\
		"							\n\t"\
		/* Block 2: */\
		"movq	%[__out8],%%rax		\n\t	leaq	%c[__i4](%%rsi),%%rsi	\n\t"/* __in0 + 4*ostride */\
		"movq	%[__out9],%%rbx		\n\t"\
		"movq	%[__outa],%%rcx		\n\t"\
		"movq	%[__outb],%%rdx		\n\t"\
		"leaq	0x10(%%rsi),%%rdi	\n\t"\
		"movaps	%c[__i1](%%rsi),%%xmm4	\n\t"\
		"movaps	%c[__i3](%%rsi),%%xmm6	\n\t"\
		"movaps	%c[__i1](%%rdi),%%xmm5	\n\t"\
		"movaps	%c[__i3](%%rdi),%%xmm7	\n\t"\
		"movaps	%%xmm4,%%xmm0		\n\t"\
		"movaps	%%xmm6,%%xmm2		\n\t"\
		"movaps	%%xmm5,%%xmm1		\n\t"\
		"movaps	%%xmm7,%%xmm3		\n\t"\
		"							\n\t"\
		"movq	%[__isrt2],%%rdi	\n\t"\
		"addq	$0x10,%%rdi	\n\t"/* cc0 */\
		"mulpd	    (%%rdi),%%xmm4	\n\t"\
		"mulpd	0x10(%%rdi),%%xmm6	\n\t"\
		"mulpd	0x10(%%rdi),%%xmm1	\n\t"\
		"mulpd	    (%%rdi),%%xmm3	\n\t"\
		"mulpd	    (%%rdi),%%xmm5	\n\t"\
		"mulpd	0x10(%%rdi),%%xmm7	\n\t"\
		"mulpd	0x10(%%rdi),%%xmm0	\n\t"\
		"mulpd	    (%%rdi),%%xmm2	\n\t"\
		"subpd	%%xmm1,%%xmm4		\n\t"\
		"subpd	%%xmm3,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm5		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"							\n\t"\
		"subpd	%%xmm6,%%xmm4		\n\t"\
		"subpd	%%xmm7,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm6		\n\t"\
		"addpd	%%xmm5,%%xmm7		\n\t"\
		"							\n\t"\
		"leaq	0x10(%%rsi),%%rdi	\n\t"\
		"movaps	%c[__i2](%%rsi),%%xmm2	\n\t"\
		"movaps	%c[__i2](%%rdi),%%xmm3	\n\t"\
		"subpd	%c[__i2](%%rdi),%%xmm2	\n\t"\
		"addpd	%c[__i2](%%rsi),%%xmm3	\n\t"\
		"movq	%[__isrt2],%%rdi	\n\t"\
		"mulpd	(%%rdi),%%xmm2	\n\t"/* mul by isrt2 */\
		"mulpd	(%%rdi),%%xmm3	\n\t"\
		"							\n\t"\
		"leaq	0x10(%%rsi),%%rdi	\n\t"\
		"movaps	        (%%rsi),%%xmm0	\n\t"\
		"movaps	        (%%rdi),%%xmm1	\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"addpd	        (%%rsi),%%xmm2	\n\t"\
		"addpd	        (%%rdi),%%xmm3	\n\t"\
		"							\n\t"\
		"subpd	%%xmm6,%%xmm2		\n\t"\
		"subpd	%%xmm7,%%xmm3		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"movaps	%%xmm2,    (%%rbx)	\n\t"\
		"movaps	%%xmm3,0x10(%%rbx)	\n\t"\
		"addpd	%%xmm2,		%%xmm6	\n\t"\
		"addpd	%%xmm3,		%%xmm7	\n\t"\
		"movaps	%%xmm6,    (%%rax)	\n\t"\
		"movaps	%%xmm7,0x10(%%rax)	\n\t"\
		"							\n\t"\
		"subpd	%%xmm5,		%%xmm0	\n\t"\
		"subpd	%%xmm4,		%%xmm1	\n\t"\
		"addpd	%%xmm5,		%%xmm5	\n\t"\
		"addpd	%%xmm4,		%%xmm4	\n\t"\
		"movaps	%%xmm0,    (%%rcx)	\n\t"\
		"movaps	%%xmm1,0x10(%%rdx)	\n\t"\
		"addpd	%%xmm0,		%%xmm5	\n\t"\
		"addpd	%%xmm1,		%%xmm4	\n\t"\
		"movaps	%%xmm5,    (%%rdx)	\n\t"\
		"movaps	%%xmm4,0x10(%%rcx)	\n\t"\
		"							\n\t"\
		/* Block 1: r8-b */\
		"movq	%[__out4],%%rax		\n\t	leaq	%c[__i4](%%rsi),%%rsi	\n\t"/* __in0 + 8*ostride */\
		"movq	%[__out5],%%rbx		\n\t"\
		"movq	%[__out6],%%rcx		\n\t"\
		"movq	%[__out7],%%rdx		\n\t"\
		"leaq	0x10(%%rsi),%%rdi	\n\t"\
		"movaps	        (%%rsi),%%xmm0	\n\t"\
		"movaps	        (%%rdi),%%xmm1	\n\t"\
		"movaps	%c[__i2](%%rsi),%%xmm2	\n\t"\
		"movaps	%c[__i2](%%rdi),%%xmm3	\n\t"\
		"						\n\t"\
		"subpd	%c[__i2](%%rdi),%%xmm0	\n\t"\
		"subpd	%c[__i2](%%rsi),%%xmm1	\n\t"\
		"addpd	        (%%rdi),%%xmm2	\n\t"\
		"addpd	        (%%rsi),%%xmm3	\n\t"\
		"						\n\t"\
		"movaps	%c[__i1](%%rsi),%%xmm4	\n\t"\
		"movaps	%c[__i1](%%rdi),%%xmm5	\n\t"\
		"movaps	%c[__i3](%%rsi),%%xmm6	\n\t"\
		"movaps	%c[__i3](%%rdi),%%xmm7	\n\t"\
		"						\n\t"\
		"subpd	%c[__i1](%%rdi),%%xmm4	\n\t"\
		"addpd	%c[__i1](%%rsi),%%xmm5	\n\t"\
		"addpd	%c[__i3](%%rdi),%%xmm6	\n\t"\
		"subpd	%c[__i3](%%rsi),%%xmm7	\n\t"\
		"							\n\t"\
		"movq	%[__isrt2],%%rdi	\n\t"\
		"mulpd	(%%rdi),%%xmm4		\n\t"\
		"mulpd	(%%rdi),%%xmm5		\n\t"\
		"mulpd	(%%rdi),%%xmm6		\n\t"\
		"mulpd	(%%rdi),%%xmm7		\n\t"\
		"							\n\t"\
		"subpd	%%xmm6,%%xmm4		\n\t"\
		"subpd	%%xmm7,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm6		\n\t"\
		"addpd	%%xmm5,%%xmm7		\n\t"\
		"							\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm2		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"							\n\t"\
		"movaps	%%xmm0,    (%%rbx)	\n\t"\
		"movaps	%%xmm2,0x10(%%rbx)	\n\t"\
		"addpd	%%xmm0,		%%xmm4	\n\t"\
		"addpd	%%xmm2,		%%xmm5	\n\t"\
		"movaps	%%xmm4,    (%%rax)	\n\t"\
		"movaps	%%xmm5,0x10(%%rax)	\n\t"\
		"							\n\t"\
		"subpd	%%xmm7,		%%xmm3	\n\t"\
		"subpd	%%xmm6,		%%xmm1	\n\t"\
		"addpd	%%xmm7,		%%xmm7	\n\t"\
		"addpd	%%xmm6,		%%xmm6	\n\t"\
		"movaps	%%xmm3,    (%%rcx)	\n\t"\
		"movaps	%%xmm1,0x10(%%rdx)	\n\t"\
		"addpd	%%xmm3,		%%xmm7	\n\t"\
		"addpd	%%xmm1,		%%xmm6	\n\t"\
		"movaps	%%xmm7,    (%%rdx)	\n\t"\
		"movaps	%%xmm6,0x10(%%rcx)	\n\t"\
		"							\n\t"\
		/* Block 3: */\
		"movq	%[__outc],%%rax		\n\t	leaq	%c[__i4](%%rsi),%%rsi	\n\t"/* __in0 + c*ostride */\
		"movq	%[__outd],%%rbx		\n\t"\
		"movq	%[__oute],%%rcx		\n\t"\
		"movq	%[__outf],%%rdx		\n\t"\
		"leaq	0x10(%%rsi),%%rdi	\n\t"\
		"movaps	%c[__i1](%%rsi),%%xmm4	\n\t"\
		"movaps	%c[__i3](%%rsi),%%xmm6	\n\t"\
		"movaps	%c[__i1](%%rdi),%%xmm5	\n\t"\
		"movaps	%c[__i3](%%rdi),%%xmm7	\n\t"\
		"movaps	%%xmm4,%%xmm0		\n\t"\
		"movaps	%%xmm6,%%xmm2		\n\t"\
		"movaps	%%xmm5,%%xmm1		\n\t"\
		"movaps	%%xmm7,%%xmm3		\n\t"\
		"							\n\t"\
		"movq	%[__isrt2],%%rdi	\n\t"\
		"addq	$0x10,%%rdi	\n\t"/* cc0 */\
		"mulpd	0x10(%%rdi),%%xmm4	\n\t"\
		"mulpd	    (%%rdi),%%xmm6	\n\t"\
		"mulpd	    (%%rdi),%%xmm1	\n\t"\
		"mulpd	0x10(%%rdi),%%xmm3	\n\t"\
		"mulpd	0x10(%%rdi),%%xmm5	\n\t"\
		"mulpd	    (%%rdi),%%xmm7	\n\t"\
		"mulpd	    (%%rdi),%%xmm0	\n\t"\
		"mulpd	0x10(%%rdi),%%xmm2	\n\t"\
		"subpd	%%xmm1,%%xmm4		\n\t"\
		"subpd	%%xmm3,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm5		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"							\n\t"\
		"subpd	%%xmm6,%%xmm4		\n\t"\
		"subpd	%%xmm7,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm6		\n\t"\
		"addpd	%%xmm5,%%xmm7		\n\t"\
		"							\n\t"\
		"leaq	0x10(%%rsi),%%rdi	\n\t"\
		"movaps	%c[__i2](%%rsi),%%xmm2	\n\t"\
		"movaps	%c[__i2](%%rdi),%%xmm3	\n\t"\
		"addpd	%c[__i2](%%rdi),%%xmm2	\n\t"\
		"subpd	%c[__i2](%%rsi),%%xmm3	\n\t"\
		"movq	%[__isrt2],%%rdi	\n\t"\
		"mulpd	(%%rdi),%%xmm2	\n\t"/* mul by isrt2 */\
		"mulpd	(%%rdi),%%xmm3	\n\t"\
		"							\n\t"\
		"leaq	0x10(%%rsi),%%rdi	\n\t"\
		"movaps	        (%%rsi),%%xmm0	\n\t"\
		"movaps	        (%%rdi),%%xmm1	\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"addpd	        (%%rsi),%%xmm2	\n\t"\
		"addpd	        (%%rdi),%%xmm3	\n\t"\
		"							\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"movaps	%%xmm0,    (%%rbx)	\n\t"\
		"movaps	%%xmm1,0x10(%%rbx)	\n\t"\
		"addpd	%%xmm0,		%%xmm4	\n\t"\
		"addpd	%%xmm1,		%%xmm5	\n\t"\
		"movaps	%%xmm4,    (%%rax)	\n\t"\
		"movaps	%%xmm5,0x10(%%rax)	\n\t"\
		"							\n\t"\
		"subpd	%%xmm7,		%%xmm2	\n\t"\
		"subpd	%%xmm6,		%%xmm3	\n\t"\
		"addpd	%%xmm7,		%%xmm7	\n\t"\
		"addpd	%%xmm6,		%%xmm6	\n\t"\
		"movaps	%%xmm2,    (%%rcx)	\n\t"\
		"movaps	%%xmm3,0x10(%%rdx)	\n\t"\
		"addpd	%%xmm2,		%%xmm7	\n\t"\
		"addpd	%%xmm3,		%%xmm6	\n\t"\
		"movaps	%%xmm7,    (%%rdx)	\n\t"\
		"movaps	%%xmm6,0x10(%%rcx)	\n\t"\
		:					/* outputs: none */\
		:[__in0] "m" (Xin0)	/* All inputs from memory addresses here */\
		,[__i1] "e" (Xi1)\
		,[__i2] "e" (Xi2)\
		,[__i3] "e" (Xi3)\
		,[__i4] "e" (Xi4)\
		,[__isrt2] "m" (Xisrt2)\
		,[__out0] "m" (Xout0)\
		,[__out1] "m" (Xout1)\
		,[__out2] "m" (Xout2)\
		,[__out3] "m" (Xout3)\
		,[__out4] "m" (Xout4)\
		,[__out5] "m" (Xout5)\
		,[__out6] "m" (Xout6)\
		,[__out7] "m" (Xout7)\
		,[__out8] "m" (Xout8)\
		,[__out9] "m" (Xout9)\
		,[__outa] "m" (Xouta)\
		,[__outb] "m" (Xoutb)\
		,[__outc] "m" (Xoutc)\
		,[__outd] "m" (Xoutd)\
		,[__oute] "m" (Xoute)\
		,[__outf] "m" (Xoutf)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
	);\
	}

	#define SSE2_CMUL_EXPO(XcA,XcB,XcAmB,XcApB)\
	{\
	__asm__ volatile (\
		"movq	%[__cA]		,%%rax\n\t"\
		"movq	%[__cB]		,%%rbx\n\t"\
		"movq	%[__cAmB]	,%%rcx\n\t"\
		"movq	%[__cApB]	,%%rdx\n\t"\
		"\n\t"\
		"movaps	    (%%rax),%%xmm0\n\t"\
		"movaps	0x10(%%rax),%%xmm2\n\t"\
		"movaps	    (%%rbx),%%xmm4\n\t"\
		"movaps	0x10(%%rbx),%%xmm5\n\t"\
		"movaps	%%xmm0,%%xmm1\n\t"\
		"movaps	%%xmm2,%%xmm3\n\t"\
		"\n\t"\
		"mulpd	%%xmm4,%%xmm0\n\t"\
		"mulpd	%%xmm5,%%xmm1\n\t"\
		"mulpd	%%xmm4,%%xmm2\n\t"\
		"mulpd	%%xmm5,%%xmm3\n\t"\
		"movaps	%%xmm0,%%xmm4\n\t"\
		"movaps	%%xmm1,%%xmm5\n\t"\
		"addpd	%%xmm3,%%xmm0\n\t"\
		"subpd	%%xmm2,%%xmm1\n\t"\
		"subpd	%%xmm3,%%xmm4\n\t"\
		"addpd	%%xmm2,%%xmm5\n\t"\
		"movaps	%%xmm0,    (%%rcx)\n\t"\
		"movaps	%%xmm1,0x10(%%rcx)\n\t"\
		"movaps	%%xmm4,    (%%rdx)\n\t"\
		"movaps	%%xmm5,0x10(%%rdx)\n\t"\
		:					/* outputs: none */\
		: [__cA]  "m" (XcA)	/* All inputs from memory addresses here */\
		 ,[__cB]  "m" (XcB)\
		 ,[__cAmB] "m" (XcAmB)\
		 ,[__cApB] "m" (XcApB)\
		: "cc","memory","rax","rbx","rcx","rdx","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5"		/* Clobbered registers */\
	);\
	}

	#define PAIR_SQUARE_4_SSE2(XtAr, XtBr, XtCr, XtDr, Xc, Xs, Xforth)\
	{\
	__asm__ volatile (\
		"/*   calculate cross-product terms...\n\t"\
		"	__rt=__tAr* ~tDr+__tAi* ~tDi; __rt=__rt+__rt;\n\t"\
		"	__it=__tAi* ~tDr-__tAr* ~tDi; __it=__it+__it;\n\t"\
		"*/\n\t"\
		"movq	%[__tDr]	,%%rdx\n\t"\
		"movq	%[__tAr]	,%%rax\n\t"\
		"\n\t"\
		"movaps	    (%%rdx)	,%%xmm6		/* tDr */\n\t"\
		"movaps	0x10(%%rdx)	,%%xmm7		/* tDi */\n\t"\
		"movaps	    (%%rax)	,%%xmm0		/* tAr */\n\t"\
		"movaps	0x10(%%rax)	,%%xmm3		/* tAi */\n\t"\
		"shufpd	$1	,%%xmm6	,%%xmm6	/*~tDr */\n\t"\
		"shufpd	$1	,%%xmm7	,%%xmm7	/*~tDi */\n\t"\
		"movaps	    (%%rax)	,%%xmm2		/* cpy tAr */\n\t"\
		"movaps	0x10(%%rax)	,%%xmm1		/* cpy tAi */\n\t"\
		"\n\t"\
		"mulpd	%%xmm6		,%%xmm0	/* tAr*~tDr */\n\t"\
		"mulpd	%%xmm7		,%%xmm3	/* tAi*~tDi */\n\t"\
		"mulpd	%%xmm6		,%%xmm1	/* tAi*~tDr */\n\t"\
		"mulpd	%%xmm7		,%%xmm2	/* tAr*~tDi */\n\t"\
		"addpd	%%xmm3		,%%xmm0	/* rt */\n\t"\
		"subpd	%%xmm2		,%%xmm1	/* it */\n\t"\
		"addpd	%%xmm0		,%%xmm0	/* rt=rt+rt */\n\t"\
		"addpd	%%xmm1		,%%xmm1	/* it=it+it; xmm2-7 free */\n\t"\
		"/*\n\t"\
		"	__st=__tBr* ~tCr+__tBi* ~tCi; __st=__st+__st;\n\t"\
		"	__jt=__tBi* ~tCr-__tBr* ~tCi; __jt=__jt+__jt;\n\t"\
		"*/\n\t"\
		"movq	%[__tCr]	,%%rcx\n\t"\
		"movq	%[__tBr]	,%%rbx\n\t"\
		"\n\t"\
		"movaps	    (%%rcx)	,%%xmm6		/* tCr */\n\t"\
		"movaps	0x10(%%rcx)	,%%xmm7		/* tCi */\n\t"\
		"movaps	    (%%rbx)	,%%xmm2		/* tBr */\n\t"\
		"movaps	0x10(%%rbx)	,%%xmm5		/* tBi */\n\t"\
		"shufpd	$1	,%%xmm6	,%%xmm6	/*~tCr */\n\t"\
		"shufpd	$1	,%%xmm7	,%%xmm7	/*~tCi */\n\t"\
		"movaps	    (%%rbx)	,%%xmm4		/* cpy tBr */\n\t"\
		"movaps	0x10(%%rbx)	,%%xmm3		/* cpy tBi */\n\t"\
		"\n\t"\
		"mulpd	%%xmm6		,%%xmm2	/* tBr*~tCr */\n\t"\
		"mulpd	%%xmm7		,%%xmm5	/* tBi*~tCi */\n\t"\
		"mulpd	%%xmm6		,%%xmm3	/* tBi*~tCr */\n\t"\
		"mulpd	%%xmm7		,%%xmm4	/* tBr*~tCi */\n\t"\
		"addpd	%%xmm5		,%%xmm2	/* st */\n\t"\
		"subpd	%%xmm4		,%%xmm3	/* jt */\n\t"\
		"addpd	%%xmm2		,%%xmm2	/* st=st+st */\n\t"\
		"addpd	%%xmm3		,%%xmm3	/* jt=jt+jt; xmm4-7 free */\n\t"\
		"\n\t"\
		"/*   now calculate square terms and __store back in the same temporaries:	*/\n\t"\
		"/*	__tmp=(__tAr+__tAi)*(__tAr-__tAi); __tAi=__tAr*__tAi; __tAi=__tAi+__tAi; __tAr=__tmp;	*/\n\t"\
		"\n\t"\
		"movaps	    (%%rax)	,%%xmm4		/* __tAr */\n\t"\
		"movaps	0x10(%%rax)	,%%xmm5		/* __tAi */\n\t"\
		"subpd	%%xmm5		,%%xmm4		/* (__tAr-__tAi) */\n\t"\
		"addpd	%%xmm5		,%%xmm5		/*      2*__tAi  */\n\t"\
		"addpd	%%xmm4		,%%xmm5		/* (__tAr+__tAi) */\n\t"\
		"mulpd	%%xmm5		,%%xmm4		/*>__tAr */\n\t"\
		"\n\t"\
		"movaps	    (%%rax)	,%%xmm5		/* __tAr */\n\t"\
		"mulpd	0x10(%%rax)	,%%xmm5		/* __tAr*__tAi */\n\t"\
		"addpd	%%xmm5		,%%xmm5		/*>__tAi */\n\t"\
		"movaps	%%xmm4	,    (%%rax)	/* tmp store >__tAr */\n\t"\
		"movaps	%%xmm5	,0x10(%%rax)	/* tmp store >__tAi */\n\t"\
		"\n\t"\
		"subpd	%%xmm4		,%%xmm0	/* rt-__tAr */\n\t"\
		"subpd	%%xmm5		,%%xmm1	/* it-__tAi; xmm4-7 free */\n\t"\
		"\n\t"\
		"/*	__tmp=(__tBr+__tBi)*(__tBr-__tBi); __tBi=__tBr*__tBi; __tBi=__tBi+__tBi; __tBr=__tmp;	*/\n\t"\
		"/*** [Can be done in parallel with above segment] ***/\n\t"\
		"\n\t"\
		"movaps	    (%%rbx)	,%%xmm6		/* __tBr */\n\t"\
		"movaps	0x10(%%rbx)	,%%xmm7		/* __tBi */\n\t"\
		"subpd	%%xmm7		,%%xmm6		/* (__tBr-__tBi) */\n\t"\
		"addpd	%%xmm7		,%%xmm7		/*      2*__tBi  */\n\t"\
		"addpd	%%xmm6		,%%xmm7		/* (__tBr+__tBi) */\n\t"\
		"mulpd	%%xmm7		,%%xmm6		/*>__tBr */\n\t"\
		"\n\t"\
		"movaps	    (%%rbx)	,%%xmm7		/* __tBr */\n\t"\
		"mulpd	0x10(%%rbx)	,%%xmm7		/* __tBr*__tBi */\n\t"\
		"addpd	%%xmm7		,%%xmm7		/*>__tBi */\n\t"\
		"movaps	%%xmm6	,    (%%rbx)	/* tmp store >__tBr */\n\t"\
		"movaps	%%xmm7	,0x10(%%rbx)	/* tmp store >__tBi */\n\t"\
		"\n\t"\
		"subpd	%%xmm6		,%%xmm2	/* st-__tBr */\n\t"\
		"subpd	%%xmm7		,%%xmm3	/* jt-__tBi; xmm4-7 free */\n\t"\
		"\n\t"\
		"/*	__tmp=(__tDr+__tDi)*(__tDr-__tDi); __tDi=__tDr*__tDi; __tDi=__tDi+__tDi; __tDr=__tmp;	*/\n\t"\
		"\n\t"\
		"movaps	    (%%rdx)	,%%xmm4		/* __tDr */\n\t"\
		"movaps	0x10(%%rdx)	,%%xmm5		/* __tDi */\n\t"\
		"subpd	%%xmm5		,%%xmm4		/* (__tDr-__tDi) */\n\t"\
		"addpd	%%xmm5		,%%xmm5		/*      2*__tDi  */\n\t"\
		"addpd	%%xmm4		,%%xmm5		/* (__tDr+__tDi) */\n\t"\
		"mulpd	%%xmm5		,%%xmm4		/*>__tDr */\n\t"\
		"\n\t"\
		"movaps	    (%%rdx)	,%%xmm5		/* __tDr */\n\t"\
		"mulpd	0x10(%%rdx)	,%%xmm5		/* __tDr*__tDi */\n\t"\
		"addpd	%%xmm5		,%%xmm5		/*>__tDi */\n\t"\
		"movaps	%%xmm4	,    (%%rdx)	/* tmp store ~tDr */\n\t"\
		"movaps	%%xmm5	,0x10(%%rdx)	/* tmp store ~tDi */\n\t"\
		"shufpd	$1	,%%xmm4	,%%xmm4	/*~tDr */\n\t"\
		"shufpd	$1	,%%xmm5	,%%xmm5	/*~tDi */\n\t"\
		"\n\t"\
		"subpd	%%xmm4		,%%xmm0	/* rt-__tAr- ~tDr */\n\t"\
		"addpd	%%xmm5		,%%xmm1	/* it-__tAi+ ~tDi; xmm4-7 free */\n\t"\
		"\n\t"\
		"/*	__tmp=(__tCr+__tCi)*(__tCr-__tCi); __tCi=__tCr*__tCi; __tCi=__tCi+__tCi; __tCr=__tmp;	*/\n\t"\
		"/*** [Can be done in parallel with above segment] ***/\n\t"\
		"\n\t"\
		"movaps	    (%%rcx)	,%%xmm6		/* __tCr */\n\t"\
		"movaps	0x10(%%rcx)	,%%xmm7		/* __tCi */\n\t"\
		"subpd	%%xmm7		,%%xmm6		/* (__tCr-__tCi) */\n\t"\
		"addpd	%%xmm7		,%%xmm7		/*      2*__tCi  */\n\t"\
		"addpd	%%xmm6		,%%xmm7		/* (__tCr+__tCi) */\n\t"\
		"mulpd	%%xmm7		,%%xmm6		/*>__tCr */\n\t"\
		"\n\t"\
		"movaps	    (%%rcx)	,%%xmm7		/* __tCr */\n\t"\
		"mulpd	0x10(%%rcx)	,%%xmm7		/* __tCr*__tCi */\n\t"\
		"addpd	%%xmm7		,%%xmm7		/*>__tCi */\n\t"\
		"movaps	%%xmm6	,    (%%rcx)	/* tmp store ~tCr */\n\t"\
		"movaps	%%xmm7	,0x10(%%rcx)	/* tmp store ~tCi */\n\t"\
		"shufpd	$1	,%%xmm6	,%%xmm6	/*~tCr */\n\t"\
		"shufpd	$1	,%%xmm7	,%%xmm7	/*~tCi */\n\t"\
		"\n\t"\
		"subpd	%%xmm6		,%%xmm2	/* st-__tBr- ~tCr */\n\t"\
		"addpd	%%xmm7		,%%xmm3	/* jt-__tBi+ ~tCi; xmm4-7 free */\n\t"\
		"/*\n\t"\
		"	__tmp=((1.0+__c)*__rt-__s*__it)*0.25;\n\t"\
		"	__it =((1.0+__c)*__it+__s*__rt)*0.25;	__rt=__tmp;\n\t"\
		"*/\n\t"\
		"/*** [Can be done in parallel with above segment] ***/\n\t"\
		"movq	%[__c]		,%%rax\n\t"\
		"movq	%[__s]		,%%rbx\n\t"\
		"movq	%[__forth]	,%%rdx\n\t"\
		"movaps	%%xmm0		,%%xmm4		/* cpy rt */\n\t"\
		"movaps	%%xmm1		,%%xmm5		/* cpy it */\n\t"\
		"mulpd	(%%rax)		,%%xmm0		/* c*rt */\n\t"\
		"mulpd	(%%rax)		,%%xmm1		/* c*it */\n\t"\
		"addpd	%%xmm4		,%%xmm0		/* (c+1.0)*rt */\n\t"\
		"addpd	%%xmm5		,%%xmm1		/* (c+1.0)*it */\n\t"\
		"mulpd	(%%rbx)		,%%xmm4		/* s*rt */\n\t"\
		"mulpd	(%%rbx)		,%%xmm5		/* s*it */\n\t"\
		"subpd	%%xmm5		,%%xmm0		/* (c+1.0)*rt-s*it */\n\t"\
		"addpd	%%xmm4		,%%xmm1		/* (c+1.0)*it+s*rt; xmm4,5 free */\n\t"\
		"mulpd	(%%rdx)		,%%xmm0	/* -rt Both of these inherit the sign flip [w.r.to the non-SSE2 PAIR_SQUARE_4 macro] */\n\t"\
		"mulpd	(%%rdx)		,%%xmm1	/* -it that resulted from the in-place-friendlier (rt-__tAr- ~tDr) reordering above. */\n\t"\
		"/*\n\t"\
		"	__tmp=((1.0-__s)*__st-__c*__jt)*0.25;\n\t"\
		"	__jt =((1.0-__s)*__jt+__c*__st)*0.25	__st=__tmp;\n\t"\
		"*/\n\t"\
		"/*** [Can be done in parallel wjth above segment] ***/\n\t"\
		"movaps	%%xmm2		,%%xmm6		/* cpy st */\n\t"\
		"movaps	%%xmm3		,%%xmm7		/* cpy jt */\n\t"\
		"mulpd	(%%rbx)		,%%xmm2		/* s*st */\n\t"\
		"mulpd	(%%rbx)		,%%xmm3		/* s*jt */\n\t"\
		"subpd	%%xmm6		,%%xmm2		/* (s-1.0)*st, note sign flip! */\n\t"\
		"subpd	%%xmm7		,%%xmm3		/* (s-1.0)*jt, note sign flip! */\n\t"\
		"mulpd	(%%rax)		,%%xmm6		/* c*st */\n\t"\
		"mulpd	(%%rax)		,%%xmm7		/* c*jt */\n\t"\
		"addpd	%%xmm7		,%%xmm2		/* -[(1.0-s)*st-c*jt] */\n\t"\
		"subpd	%%xmm6		,%%xmm3		/* -[(1.0-s)*jt+c*st]; xmm6,7 free */\n\t"\
		"mulpd	(%%rdx)		,%%xmm2	/* +st Sign flip due to (s-1.0) reordering here */\n\t"\
		"mulpd	(%%rdx)		,%%xmm3	/* +jt cancels earlier one due to in-place-friendlier (st-__tBr- ~tCr) reordering above. */\n\t"\
		"/*...and now complete and store the results. We flip the signs on st and jt here to undo the above -st,-jt negations. */\n\t"\
		"/*	__tAr = (__tAr+__rt);\n\t"\
		"	__tAi = (__tAi+__it);\n\t"\
		"	__tBr = (__tBr-__st);\n\t"\
		"	__tBi = (__tBi-__jt);\n\t"\
		"*/\n\t"\
		"movq	%[__tAr]	,%%rax\n\t"\
		"movq	%[__tBr]	,%%rbx\n\t"\
		"\n\t"\
		"movaps	    (%%rax)	,%%xmm4		/* __tAr */\n\t"\
		"movaps	0x10(%%rax)	,%%xmm5		/* __tAi */\n\t"\
		"movaps	    (%%rbx)	,%%xmm6		/* __tBr */\n\t"\
		"movaps	0x10(%%rbx)	,%%xmm7		/* __tBi */\n\t"\
		"addpd	%%xmm0		,%%xmm4		/* (__tAr+__rt) */\n\t"\
		"addpd	%%xmm1		,%%xmm5		/* (__tAi+__it) */\n\t"\
		"subpd	%%xmm2		,%%xmm6		/* (__tBr-__st) */\n\t"\
		"subpd	%%xmm3		,%%xmm7		/* (__tBi-__jt) */\n\t"\
		"movaps	%%xmm4	,    (%%rax)	/* store >__tAr */\n\t"\
		"movaps	%%xmm5	,0x10(%%rax)	/* store >__tAi */\n\t"\
		"movaps	%%xmm6	,    (%%rbx)	/* store >__tBr */\n\t"\
		"movaps	%%xmm7	,0x10(%%rbx)	/* store >__tBi */\n\t"\
		"/*...N-j terms are as above, but with the replacements: __tAr<--> ~tDr, __tAi<--> ~tDi, __it|-->-__it. */\n\t"\
		"/*	__tDr = (__tDr+ ~rt);\n\t"\
		"	__tDi = (__tDi- ~it);\n\t"\
		"	__tCr = (__tCr- ~st);\n\t"\
		"	__tCi = (__tCi+ ~jt);\n\t"\
		"*/\n\t"\
		"movq	%[__tCr]	,%%rcx\n\t"\
		"movq	%[__tDr]	,%%rdx\n\t"\
		"\n\t"\
		"shufpd	$1	,%%xmm0	,%%xmm0		/* ~rt */\n\t"\
		"shufpd	$1	,%%xmm1	,%%xmm1		/* ~it */\n\t"\
		"shufpd	$1	,%%xmm2	,%%xmm2		/* ~st */\n\t"\
		"shufpd	$1	,%%xmm3	,%%xmm3		/* ~jt */\n\t"\
		"\n\t"\
		"movaps	    (%%rdx)	,%%xmm4		/* __tDr */\n\t"\
		"movaps	0x10(%%rdx)	,%%xmm5		/* __tDi */\n\t"\
		"movaps	    (%%rcx)	,%%xmm6		/* __tCr */\n\t"\
		"movaps	0x10(%%rcx)	,%%xmm7		/* __tCi */\n\t"\
		"addpd	%%xmm0		,%%xmm4		/* (__tDr+ ~rt) */\n\t"\
		"subpd	%%xmm1		,%%xmm5		/* (__tDi- ~it) */\n\t"\
		"subpd	%%xmm2		,%%xmm6		/* (__tCr- ~st) */\n\t"\
		"addpd	%%xmm3		,%%xmm7		/* (__tCi+ ~jt) */\n\t"\
		"movaps	%%xmm4	,    (%%rdx)	/* store >__tDr */\n\t"\
		"movaps	%%xmm5	,0x10(%%rdx)	/* store >__tDi */\n\t"\
		"movaps	%%xmm6	,    (%%rcx)	/* store >__tCr */\n\t"\
		"movaps	%%xmm7	,0x10(%%rcx)	/* store >__tCi */\n\t"\
		:					/* outputs: none */\
		: [__tAr] "m" (XtAr)	/* All inputs from memory addresses here */\
		 ,[__tBr] "m" (XtBr)\
		 ,[__tCr] "m" (XtCr)\
		 ,[__tDr] "m" (XtDr)\
		 ,[__c] "m" (Xc)\
		 ,[__s] "m" (Xs)\
		 ,[__forth] "m" (Xforth)\
		: "cc","memory","rax","rbx","rcx","rdx","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
	);\
	}

	#define SSE2_RADIX_03_DFT(Xi0,Xi1,Xi2, Xcc1, Xo0,Xo1,Xo2)\
	{\
	__asm__ volatile (\
			"movq	%[__i0],%%rax		\n\t"\
			"movq	%[__i1],%%rbx		\n\t"\
			"movq	%[__i2],%%rcx		\n\t"\
			"movq	%[__cc1],%%rdx		\n\t"\
			"\n\t"\
			"movaps	    (%%rbx),%%xmm2	\n\t"\
			"movaps	0x10(%%rbx),%%xmm3	\n\t"\
			"movaps	    (%%rax),%%xmm0	\n\t"\
			"movaps	0x10(%%rax),%%xmm1	\n\t"\
			"movaps	    (%%rcx),%%xmm6	\n\t"\
			"movaps	0x10(%%rcx),%%xmm7	\n\t"\
			"movaps	%%xmm2,%%xmm4		\n\t"\
			"movaps	%%xmm3,%%xmm5		\n\t"\
			"\n\t"\
			"movq	%[__o0],%%rax		\n\t"\
			"movq	%[__o1],%%rbx		\n\t"\
			"movq	%[__o2],%%rcx		\n\t"\
			"addpd	%%xmm6,%%xmm2		\n\t"\
			"addpd	%%xmm7,%%xmm3		\n\t"\
			"subpd	%%xmm6,%%xmm4		\n\t"\
			"subpd	%%xmm7,%%xmm5		\n\t"\
			"addpd	%%xmm2,%%xmm0		\n\t"\
			"addpd	%%xmm3,%%xmm1		\n\t"\
			"movaps	    (%%rdx),%%xmm6	\n\t"\
			"movaps	0x10(%%rdx),%%xmm7	\n\t"\
			"movaps	%%xmm0,    (%%rax)	\n\t"\
			"movaps	%%xmm1,0x10(%%rax)	\n\t"\
			"\n\t"\
			"mulpd	%%xmm6,%%xmm2		\n\t"\
			"mulpd	%%xmm6,%%xmm3		\n\t"\
			"mulpd	%%xmm7,%%xmm4		\n\t"\
			"mulpd	%%xmm7,%%xmm5		\n\t"\
			"addpd	%%xmm0,%%xmm2		\n\t"\
			"addpd	%%xmm1,%%xmm3		\n\t"\
			"\n\t"\
			"movaps	%%xmm2,%%xmm0		\n\t"\
			"movaps	%%xmm3,%%xmm1		\n\t"\
			"\n\t"\
			"subpd	%%xmm5,%%xmm2		\n\t"\
			"addpd	%%xmm4,%%xmm3		\n\t"\
			"addpd	%%xmm5,%%xmm0		\n\t"\
			"subpd	%%xmm4,%%xmm1		\n\t"\
			"\n\t"\
			"movaps	%%xmm2,    (%%rbx)	\n\t"\
			"movaps	%%xmm3,0x10(%%rbx)	\n\t"\
			"movaps	%%xmm0,    (%%rcx)	\n\t"\
			"movaps	%%xmm1,0x10(%%rcx)	\n\t"\
		:					/* outputs: none */\
		: [__i0] "m" (Xi0)	/* All inputs from memory addresses here */\
		 ,[__i1] "m" (Xi1)\
		 ,[__i2] "m" (Xi2)\
		 ,[__cc1] "m" (Xcc1)\
		 ,[__o0] "m" (Xo0)\
		 ,[__o1] "m" (Xo1)\
		 ,[__o2] "m" (Xo2)\
		: "cc","memory","rax","rbx","rcx","rdx","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
	);\
	}

	#define SSE2_RADIX_03_DFT_X2(Xcc0, Xi0,Xi1,Xi2, Xo0,Xo1,Xo2, Xj0,Xj1,Xj2, Xu0,Xu1,Xu2)\
	{\
	__asm__ volatile (\
		"movq	%[__i0],%%rax		\n\t	movq	%[__j0],%%r10		\n\t"\
		"movq	%[__i1],%%rbx		\n\t	movq	%[__j1],%%r11		\n\t"\
		"movq	%[__i2],%%rcx		\n\t	movq	%[__j2],%%r12		\n\t"\
		"movq	%[__cc0],%%rdx		\n\t"\
		"movaps	    (%%rbx),%%xmm2	\n\t	movaps	    (%%r11),%%xmm10	\n\t"\
		"movaps	0x10(%%rbx),%%xmm3	\n\t	movaps	0x10(%%r11),%%xmm11	\n\t"\
		"movaps	    (%%rax),%%xmm0	\n\t	movaps	    (%%r10),%%xmm8 	\n\t"\
		"movaps	0x10(%%rax),%%xmm1	\n\t	movaps	0x10(%%r10),%%xmm9 	\n\t"\
		"movaps	    (%%rcx),%%xmm6	\n\t	movaps	    (%%r12),%%xmm14	\n\t"\
		"movaps	0x10(%%rcx),%%xmm7	\n\t	movaps	0x10(%%r12),%%xmm15	\n\t"\
		"movaps	%%xmm2,%%xmm4		\n\t	movaps	%%xmm10,%%xmm12		\n\t"\
		"movaps	%%xmm3,%%xmm5		\n\t	movaps	%%xmm11,%%xmm13		\n\t"\
		"movq	%[__o0],%%rax		\n\t	movq	%[__u0],%%r10		\n\t"\
		"movq	%[__o1],%%rbx		\n\t	movq	%[__u1],%%r11		\n\t"\
		"movq	%[__o2],%%rcx		\n\t	movq	%[__u2],%%r12		\n\t"\
		"addpd	%%xmm6,%%xmm2		\n\t	addpd	%%xmm14,%%xmm10		\n\t"\
		"addpd	%%xmm7,%%xmm3		\n\t	addpd	%%xmm15,%%xmm11		\n\t"\
		"subpd	%%xmm6,%%xmm4		\n\t	subpd	%%xmm14,%%xmm12		\n\t"\
		"subpd	%%xmm7,%%xmm5		\n\t	subpd	%%xmm15,%%xmm13		\n\t"\
		"addpd	%%xmm2,%%xmm0		\n\t	addpd	%%xmm10,%%xmm8 		\n\t"\
		"addpd	%%xmm3,%%xmm1		\n\t	addpd	%%xmm11,%%xmm9 		\n\t"\
		"movaps	    (%%rdx),%%xmm6	\n\t"\
		"movaps	0x10(%%rdx),%%xmm7	\n\t"\
		"movaps	%%xmm0,     (%%rax)	\n\t	movaps	%%xmm8 ,     (%%r10)\n\t"\
		"movaps	%%xmm1,0x010(%%rax)	\n\t	movaps	%%xmm9 ,0x010(%%r10)\n\t"\
		"mulpd	%%xmm6,%%xmm2		\n\t	mulpd	%%xmm6 ,%%xmm10		\n\t"\
		"mulpd	%%xmm6,%%xmm3		\n\t	mulpd	%%xmm6 ,%%xmm11		\n\t"\
		"mulpd	%%xmm7,%%xmm4		\n\t	mulpd	%%xmm7 ,%%xmm12		\n\t"\
		"mulpd	%%xmm7,%%xmm5		\n\t	mulpd	%%xmm7 ,%%xmm13		\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t	addpd	%%xmm8 ,%%xmm10		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t	addpd	%%xmm9 ,%%xmm11		\n\t"\
		"movaps	%%xmm2,%%xmm0		\n\t	movaps	%%xmm10,%%xmm8 		\n\t"\
		"movaps	%%xmm3,%%xmm1		\n\t	movaps	%%xmm11,%%xmm9 		\n\t"\
		"subpd	%%xmm5,%%xmm2		\n\t	subpd	%%xmm13,%%xmm10		\n\t"\
		"addpd	%%xmm4,%%xmm3		\n\t	addpd	%%xmm12,%%xmm11		\n\t"\
		"addpd	%%xmm5,%%xmm0		\n\t	addpd	%%xmm13,%%xmm8 		\n\t"\
		"subpd	%%xmm4,%%xmm1		\n\t	subpd	%%xmm12,%%xmm9 		\n\t"\
		"movaps	%%xmm2,     (%%rbx)	\n\t	movaps	%%xmm10,     (%%r11)\n\t"\
		"movaps	%%xmm3,0x010(%%rbx)	\n\t	movaps	%%xmm11,0x010(%%r11)\n\t"\
		"movaps	%%xmm0,     (%%rcx)	\n\t	movaps	%%xmm8 ,     (%%r12)\n\t"\
		"movaps	%%xmm1,0x010(%%rcx)	\n\t	movaps	%%xmm9 ,0x010(%%r12)\n\t"\
		:					/* outputs: none */\
		: [__cc0] "m" (Xcc0)	/* All inputs from memory addresses here */\
		 ,[__i0] "m" (Xi0)\
		 ,[__i1] "m" (Xi1)\
		 ,[__i2] "m" (Xi2)\
		 ,[__o0] "m" (Xo0)\
		 ,[__o1] "m" (Xo1)\
		 ,[__o2] "m" (Xo2)\
		 ,[__j0] "m" (Xj0)\
		 ,[__j1] "m" (Xj1)\
		 ,[__j2] "m" (Xj2)\
		 ,[__u0] "m" (Xu0)\
		 ,[__u1] "m" (Xu1)\
		 ,[__u2] "m" (Xu2)\
		: "cc","memory","rax","rbx","rcx","rdx","r10","r11","r12","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"		/* Clobbered registers */\
	);\
	}

	#define SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(Xadd0, Xadd1, Xadd2, Xadd3, Xtmp, Xstride)\
	{\
	__asm__ volatile (\
		"movq	%[__tmp]   ,%%rax	\n\t"\
		"movq	%[__stride],%%rsi	\n\t"\
		"movq	%%rax,%%rbx			\n\t"\
		"addq	%%rsi,%%rbx			/* add_in1  */\n\t"\
		"shlq	$1,%%rsi			/* stride*2 */\n\t"\
		"movaps	    (%%rax),%%xmm0	\n\t"\
		"movaps	    (%%rbx),%%xmm2	\n\t"\
		"movaps	0x10(%%rax),%%xmm1	\n\t"\
		"movaps	0x10(%%rbx),%%xmm3	\n\t"\
		"movaps	    (%%rax),%%xmm4	\n\t"\
		"movaps	    (%%rbx),%%xmm6	\n\t"\
		"movaps	0x10(%%rax),%%xmm5	\n\t"\
		"movaps	0x10(%%rbx),%%xmm7	\n\t"\
		"addq	%%rsi,%%rax			/* add_in2  */\n\t"\
		"addq	%%rsi,%%rbx			/* add_in3  */\n\t"\
		"addpd	    (%%rax),%%xmm0	\n\t"\
		"addpd	    (%%rbx),%%xmm2	\n\t"\
		"addpd	0x10(%%rax),%%xmm1	\n\t"\
		"addpd	0x10(%%rbx),%%xmm3	\n\t"\
		"subpd	    (%%rax),%%xmm4	\n\t"\
		"subpd	    (%%rbx),%%xmm6	\n\t"\
		"subpd	0x10(%%rax),%%xmm5	\n\t"\
		"subpd	0x10(%%rbx),%%xmm7	\n\t"\
		"/* Finish radix-4 butterfly and store results into main-array slots: */\n\t"\
		"movq	%[__add0],%%rax		\n\t"\
		"movq	%[__add1],%%rbx		\n\t"\
		"movq	%[__add2],%%rcx		\n\t"\
		"movq	%[__add3],%%rdx		\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm4		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"subpd	%%xmm6,%%xmm5		\n\t"\
		"movaps	%%xmm0,     (%%rbx)	\n\t"\
		"movaps	%%xmm4,     (%%rcx)	\n\t"\
		"movaps	%%xmm1,0x010(%%rbx)	\n\t"\
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
		: "cc","memory","rax","rbx","rcx","rdx","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
	);\
	}

	/* DIF radix-4 subconvolution, sans twiddles, inputs in __i0-3, outputs in __o0-3, possibly coincident with inputs: */
	#define SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_E(Xi0,Xi1,Xi2,Xi3, Xo0,Xo1,Xo2,Xo3)\
	{\
	__asm__ volatile (\
		"movq	%[__i0],%%rax		\n\t"\
		"movq	%[__i1],%%rbx		\n\t"\
		"movq	%[__i2],%%rcx		\n\t"\
		"movq	%[__i3],%%rdx		\n\t"\
		"movaps	    (%%rax),%%xmm0	\n\t"\
		"movaps	    (%%rbx),%%xmm4	\n\t"\
		"movaps	0x10(%%rax),%%xmm1	\n\t"\
		"movaps	0x10(%%rbx),%%xmm5	\n\t"\
		"movaps	%%xmm0,%%xmm2	\n\t"\
		"movaps	%%xmm4,%%xmm6	\n\t"\
		"movaps	%%xmm1,%%xmm3	\n\t"\
		"movaps	%%xmm5,%%xmm7	\n\t"\
		"addpd	    (%%rcx),%%xmm0	\n\t"\
		"addpd	    (%%rdx),%%xmm4	\n\t"\
		"addpd	0x10(%%rcx),%%xmm1	\n\t"\
		"addpd	0x10(%%rdx),%%xmm5	\n\t"\
		"subpd	    (%%rcx),%%xmm2	\n\t"\
		"subpd	    (%%rdx),%%xmm6	\n\t"\
		"subpd	0x10(%%rcx),%%xmm3	\n\t"\
		"subpd	0x10(%%rdx),%%xmm7	\n\t"\
		"/* Finish radix-4 butterfly and store results into main-array slots: */\n\t"\
		"movq	%[__o0],%%rax		\n\t"\
		"movq	%[__o1],%%rbx		\n\t"\
		"movq	%[__o2],%%rcx		\n\t"\
		"movq	%[__o3],%%rdx		\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm2		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"subpd	%%xmm6,%%xmm3		\n\t"\
		"movaps	%%xmm0,    (%%rbx)	\n\t"\
		"movaps	%%xmm2,    (%%rcx)	\n\t"\
		"movaps	%%xmm1,0x10(%%rbx)	\n\t"\
		"movaps	%%xmm3,0x10(%%rdx)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"addpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	%%xmm4,    (%%rax)	\n\t"\
		"movaps	%%xmm7,    (%%rdx)	\n\t"\
		"movaps	%%xmm5,0x10(%%rax)	\n\t"\
		"movaps	%%xmm6,0x10(%%rcx)	\n\t"\
		:					/* outputs: none */\
		: [__i0] "m" (Xi0)	/* All inputs from memory addresses here */\
		 ,[__i1] "m" (Xi1)\
		 ,[__i2] "m" (Xi2)\
		 ,[__i3] "m" (Xi3)\
		 ,[__o0] "m" (Xo0)\
		 ,[__o1] "m" (Xo1)\
		 ,[__o2] "m" (Xo2)\
		 ,[__o3] "m" (Xo3)\
		: "cc","memory","rax","rbx","rcx","rdx","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
	);\
	}

	#define SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(Xadd0, Xadd1, Xadd2, Xadd3, Xtmp, Xstride)\
	{\
	__asm__ volatile (\
		"movq	%[__add0],%%rax		\n\t"\
		"movq	%[__add1],%%rbx		\n\t"\
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
		"addpd	    (%%rbx),%%xmm0	\n\t"\
		"addpd	    (%%rdx),%%xmm4	\n\t"\
		"addpd	0x10(%%rbx),%%xmm1	\n\t"\
		"addpd	0x10(%%rdx),%%xmm5	\n\t"\
		"subpd	    (%%rbx),%%xmm2	\n\t"\
		"subpd	    (%%rdx),%%xmm6	\n\t"\
		"subpd	0x10(%%rbx),%%xmm3	\n\t"\
		"subpd	0x10(%%rdx),%%xmm7	\n\t"\
		"movq	%%rax,%%rbx			\n\t"\
		"addq	%%rcx,%%rbx			\n\t"\
		"movq	%%rbx,%%rdx			\n\t"\
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
		"movaps	%%xmm3,0x010(%%rbx)	\n\t"\
		"addpd	%%xmm4,%%xmm4			\n\t"\
		"addpd	%%xmm7,%%xmm7			\n\t"\
		"addpd	%%xmm5,%%xmm5			\n\t"\
		"addpd	%%xmm6,%%xmm6			\n\t"\
		"addpd	%%xmm0,%%xmm4			\n\t"\
		"addpd	%%xmm2,%%xmm7			\n\t"\
		"addpd	%%xmm1,%%xmm5			\n\t"\
		"addpd	%%xmm3,%%xmm6			\n\t"\
		"movaps	%%xmm4,     (%%rax)	\n\t"\
		"movaps	%%xmm7,     (%%rbx)	\n\t"\
		"movaps	%%xmm5,0x010(%%rax)	\n\t"\
		"movaps	%%xmm6,0x010(%%rdx)	\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__add1] "m" (Xadd1)\
		 ,[__add2] "m" (Xadd2)\
		 ,[__add3] "m" (Xadd3)\
		 ,[__tmp] "m" (Xtmp)\
		 ,[__stride] "e" (Xstride)\
		: "cc","memory","rax","rbx","rcx","rdx","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
	);\
	}

	/* DIT radix-4 subconvolution, sans twiddles, inputs in __i0-3, outputs in __o0-3, possibly coincident with inputs: */
	#define SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_E(Xi0,Xi1,Xi2,Xi3, Xo0,Xo1,Xo2,Xo3)\
	{\
	__asm__ volatile (\
		"movq	%[__i0],%%rax		\n\t"\
		"movq	%[__i1],%%rbx		\n\t"\
		"movq	%[__i2],%%rcx		\n\t"\
		"movq	%[__i3],%%rdx		\n\t"\
		"movaps	    (%%rax),%%xmm0	\n\t"\
		"movaps	    (%%rcx),%%xmm4	\n\t"\
		"movaps	0x10(%%rax),%%xmm1	\n\t"\
		"movaps	0x10(%%rcx),%%xmm5	\n\t"\
		"movaps	%%xmm0,%%xmm2			\n\t"\
		"movaps	%%xmm4,%%xmm6			\n\t"\
		"movaps	%%xmm1,%%xmm3			\n\t"\
		"movaps	%%xmm5,%%xmm7			\n\t"\
		"addpd	    (%%rbx),%%xmm0	\n\t"\
		"addpd	    (%%rdx),%%xmm4	\n\t"\
		"addpd	0x10(%%rbx),%%xmm1	\n\t"\
		"addpd	0x10(%%rdx),%%xmm5	\n\t"\
		"subpd	    (%%rbx),%%xmm2	\n\t"\
		"subpd	    (%%rdx),%%xmm6	\n\t"\
		"subpd	0x10(%%rbx),%%xmm3	\n\t"\
		"subpd	0x10(%%rdx),%%xmm7	\n\t"\
		"/* Finish radix-4 butterfly and store results into output-array slots: */\n\t"\
		"movq	%[__o0],%%rax		\n\t"\
		"movq	%[__o1],%%rbx		\n\t"\
		"movq	%[__o2],%%rcx		\n\t"\
		"movq	%[__o3],%%rdx		\n\t"\
		"subpd	%%xmm4,%%xmm0			\n\t"\
		"subpd	%%xmm7,%%xmm2			\n\t"\
		"subpd	%%xmm5,%%xmm1			\n\t"\
		"subpd	%%xmm6,%%xmm3			\n\t"\
		"movaps	%%xmm0,     (%%rcx)	\n\t"\
		"movaps	%%xmm2,     (%%rdx)	\n\t"\
		"movaps	%%xmm1,0x010(%%rcx)	\n\t"\
		"movaps	%%xmm3,0x010(%%rbx)	\n\t"\
		"addpd	%%xmm4,%%xmm4			\n\t"\
		"addpd	%%xmm7,%%xmm7			\n\t"\
		"addpd	%%xmm5,%%xmm5			\n\t"\
		"addpd	%%xmm6,%%xmm6			\n\t"\
		"addpd	%%xmm0,%%xmm4			\n\t"\
		"addpd	%%xmm2,%%xmm7			\n\t"\
		"addpd	%%xmm1,%%xmm5			\n\t"\
		"addpd	%%xmm3,%%xmm6			\n\t"\
		"movaps	%%xmm4,     (%%rax)	\n\t"\
		"movaps	%%xmm7,     (%%rbx)	\n\t"\
		"movaps	%%xmm5,0x010(%%rax)	\n\t"\
		"movaps	%%xmm6,0x010(%%rdx)	\n\t"\
		:					/* outputs: none */\
		: [__i0] "m" (Xi0)	/* All inputs from memory addresses here */\
		 ,[__i1] "m" (Xi1)\
		 ,[__i2] "m" (Xi2)\
		 ,[__i3] "m" (Xi3)\
		 ,[__o0] "m" (Xo0)\
		 ,[__o1] "m" (Xo1)\
		 ,[__o2] "m" (Xo2)\
		 ,[__o3] "m" (Xo3)\
		: "cc","memory","rax","rbx","rcx","rdx","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
	);\
	}

	#define SSE2_RADIX_05_DFT_0TWIDDLE(Xi0,Xi1,Xi2,Xi3,Xi4, Xcc1, Xo0,Xo1,Xo2,Xo3,Xo4)\
	{\
	__asm__ volatile (\
		"movq	%[__i0],%%rsi		\n\t"\
		"movq	%[__i1],%%rax		\n\t"\
		"movq	%[__i2],%%rbx		\n\t"\
		"movq	%[__i3],%%rcx		\n\t"\
		"movq	%[__i4],%%rdx		\n\t"\
		"movq	%[__o0],%%rdi		\n\t"\
		"movaps	    (%%rax),%%xmm0	\n\t"\
		"movaps	0x10(%%rax),%%xmm1	\n\t"\
		"movaps	    (%%rbx),%%xmm2	\n\t"\
		"movaps	0x10(%%rbx),%%xmm3	\n\t"\
		"movaps	    (%%rcx),%%xmm4	\n\t"\
		"movaps	0x10(%%rcx),%%xmm5	\n\t"\
		"movaps	    (%%rdx),%%xmm6	\n\t"\
		"movaps	0x10(%%rdx),%%xmm7	\n\t"\
		"subpd	%%xmm6,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm1		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm0,%%xmm6		\n\t"\
		"addpd	%%xmm1,%%xmm7		\n\t"\
		"subpd	%%xmm4,%%xmm2		\n\t"\
		"subpd	%%xmm5,%%xmm3		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm2,%%xmm4		\n\t"\
		"addpd	%%xmm3,%%xmm5		\n\t"\
	"movq	%[__cc1],%%rax		\n\t"\
		"subpd	%%xmm4,%%xmm6		\n\t"\
		"subpd	%%xmm5,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm5		\n\t"\
		"addpd	    (%%rsi),%%xmm4	\n\t"\
		"addpd	0x10(%%rsi),%%xmm5	\n\t"\
		"movaps	%%xmm4,    (%%rdi)	\n\t"\
		"movaps	%%xmm5,0x10(%%rdi)	\n\t"\
		"mulpd	0x10(%%rax),%%xmm6	\n\t"\
		"mulpd	0x10(%%rax),%%xmm7	\n\t"\
		"subpd	     (%%rsi),%%xmm4	\n\t"\
		"subpd	0x010(%%rsi),%%xmm5	\n\t"\
		"mulpd	    (%%rax),%%xmm4	\n\t"\
		"mulpd	    (%%rax),%%xmm5	\n\t"\
		"addpd	     (%%rdi),%%xmm4	\n\t"\
		"addpd	0x010(%%rdi),%%xmm5	\n\t"\
		"subpd	%%xmm6,%%xmm4		\n\t"\
		"subpd	%%xmm7,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm6		\n\t"\
		"addpd	%%xmm5,%%xmm7		\n\t"\
		"movaps	%%xmm4,    (%%rsi)	\n\t"\
		"movaps	%%xmm5,0x10(%%rsi)	\n\t"\
		"movaps	%%xmm0,%%xmm4		\n\t"\
		"movaps	%%xmm1,%%xmm5		\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"mulpd	0x20(%%rax),%%xmm0	\n\t"\
		"mulpd	0x20(%%rax),%%xmm1	\n\t"\
		"mulpd	0x30(%%rax),%%xmm2	\n\t"\
		"mulpd	0x30(%%rax),%%xmm3	\n\t"\
		"mulpd	0x40(%%rax),%%xmm4	\n\t"\
		"mulpd	0x40(%%rax),%%xmm5	\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"movaps	    (%%rsi),%%xmm4	\n\t"\
		"movaps	0x10(%%rsi),%%xmm5	\n\t"\
		"movq	%[__o1],%%rax		\n\t"\
		"movq	%[__o4],%%rdx		\n\t"\
		"subpd	%%xmm3,%%xmm6		\n\t"\
		"subpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"movaps	%%xmm6,    (%%rax)	\n\t"\
		"movaps	%%xmm7,0x10(%%rdx)	\n\t"\
		"addpd	%%xmm6,%%xmm3		\n\t"\
		"addpd	%%xmm7,%%xmm2		\n\t"\
		"movaps	%%xmm3,    (%%rdx)	\n\t"\
		"movaps	%%xmm2,0x10(%%rax)	\n\t"\
		"movq	%[__o2],%%rbx		\n\t"\
		"movq	%[__o3],%%rcx		\n\t"\
		"subpd	%%xmm1,%%xmm4		\n\t"\
		"subpd	%%xmm0,%%xmm5		\n\t"\
		"addpd	%%xmm1,%%xmm1		\n\t"\
		"addpd	%%xmm0,%%xmm0		\n\t"\
		"movaps	%%xmm4,    (%%rbx)	\n\t"\
		"movaps	%%xmm5,0x10(%%rcx)	\n\t"\
		"addpd	%%xmm4,%%xmm1		\n\t"\
		"addpd	%%xmm5,%%xmm0		\n\t"\
		"movaps	%%xmm1,    (%%rcx)	\n\t"\
		"movaps	%%xmm0,0x10(%%rbx)	\n\t"\
		"							\n\t"\
		:					/* outputs: none */\
		: [__i0] "m" (Xi0)	/* All inputs from memory addresses here */\
		 ,[__i1] "m" (Xi1)\
		 ,[__i2] "m" (Xi2)\
		 ,[__i3] "m" (Xi3)\
		 ,[__i4] "m" (Xi4)\
		 ,[__cc1] "m" (Xcc1)\
		 ,[__o0] "m" (Xo0)\
		 ,[__o1] "m" (Xo1)\
		 ,[__o2] "m" (Xo2)\
		 ,[__o3] "m" (Xo3)\
		 ,[__o4] "m" (Xo4)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
	);\
	}

	/* 16-xmm-register version does 2 of the above side-by-side: */
	#define SSE2_RADIX_05_DFT_0TWID_X2(Xcc1, Xi0,Xi1,Xi2,Xi3,Xi4, Xo0,Xo1,Xo2,Xo3,Xo4, Yi0,Yi1,Yi2,Yi3,Yi4, Yo0,Yo1,Yo2,Yo3,Yo4)\
	{\
	__asm__ volatile (\
		"movq	%[__i0],%%rsi		\n\t"\
		"movq	%[__i1],%%rax		\n\t"\
		"movq	%[__i2],%%rbx		\n\t"\
		"movq	%[__i3],%%rcx		\n\t"\
		"movq	%[__i4],%%rdx		\n\t"\
		"movq	%[__o0],%%rdi		\n\t"\
		"movaps	    (%%rax),%%xmm0	\n\t		movq	%[__i5],%%r10		\n\t"\
		"movaps	0x10(%%rax),%%xmm1	\n\t		movq	%[__i6],%%r11		\n\t"\
		"movaps	    (%%rbx),%%xmm2	\n\t		movq	%[__i7],%%r12		\n\t"\
		"movaps	0x10(%%rbx),%%xmm3	\n\t		movq	%[__i8],%%r13		\n\t"\
		"movaps	    (%%rcx),%%xmm4	\n\t		movq	%[__i9],%%r14		\n\t"\
		"movaps	0x10(%%rcx),%%xmm5	\n\t		movq	%[__o5],%%r15		\n\t"\
		"movaps	    (%%rdx),%%xmm6	\n\t		movaps	    (%%r11),%%xmm8 	\n\t"\
		"movaps	0x10(%%rdx),%%xmm7	\n\t		movaps	0x10(%%r11),%%xmm9 	\n\t"\
		"subpd	%%xmm6,%%xmm0		\n\t		movaps	    (%%r12),%%xmm10	\n\t"\
		"subpd	%%xmm7,%%xmm1		\n\t		movaps	0x10(%%r12),%%xmm11	\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t		movaps	    (%%r13),%%xmm12	\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t		movaps	0x10(%%r13),%%xmm13	\n\t"\
		"addpd	%%xmm0,%%xmm6		\n\t		movaps	    (%%r14),%%xmm14	\n\t"\
		"addpd	%%xmm1,%%xmm7		\n\t		movaps	0x10(%%r14),%%xmm15	\n\t"\
		"subpd	%%xmm4,%%xmm2		\n\t		subpd	%%xmm14,%%xmm8 		\n\t"\
		"subpd	%%xmm5,%%xmm3		\n\t		subpd	%%xmm15,%%xmm9 		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t		addpd	%%xmm14,%%xmm14		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t		addpd	%%xmm15,%%xmm15		\n\t"\
		"addpd	%%xmm2,%%xmm4		\n\t		addpd	%%xmm8 ,%%xmm14		\n\t"\
		"addpd	%%xmm3,%%xmm5		\n\t		addpd	%%xmm9 ,%%xmm15		\n\t"\
		"movq	%[__cc1],%%rax		/* Sincos data shared by both streams */\n\t"\
		"subpd	%%xmm4,%%xmm6		\n\t		subpd	%%xmm12,%%xmm10		\n\t"\
		"subpd	%%xmm5,%%xmm7		\n\t		subpd	%%xmm13,%%xmm11		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t		addpd	%%xmm12,%%xmm12		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t		addpd	%%xmm13,%%xmm13		\n\t"\
		"addpd	%%xmm6,%%xmm4		\n\t		addpd	%%xmm10,%%xmm12		\n\t"\
		"addpd	%%xmm7,%%xmm5		\n\t		addpd	%%xmm11,%%xmm13		\n\t"\
		"addpd	    (%%rsi),%%xmm4	\n\t		subpd	%%xmm12,%%xmm14		\n\t"\
		"addpd	0x10(%%rsi),%%xmm5	\n\t		subpd	%%xmm13,%%xmm15		\n\t"\
		"movaps	%%xmm4,    (%%rdi)	\n\t		addpd	%%xmm12,%%xmm12		\n\t"\
		"movaps	%%xmm5,0x10(%%rdi)	\n\t		addpd	%%xmm13,%%xmm13		\n\t"\
		"mulpd	0x10(%%rax),%%xmm6	\n\t		addpd	%%xmm14,%%xmm12		\n\t"\
		"mulpd	0x10(%%rax),%%xmm7	\n\t		addpd	%%xmm15,%%xmm13		\n\t"\
		"subpd	    (%%rsi),%%xmm4	\n\t		addpd	    (%%r10),%%xmm12	\n\t"\
		"subpd	0x10(%%rsi),%%xmm5	\n\t		addpd	0x10(%%r10),%%xmm13	\n\t"\
		"mulpd	    (%%rax),%%xmm4	\n\t		movaps	%%xmm12,    (%%r15)	\n\t"\
		"mulpd	    (%%rax),%%xmm5	\n\t		movaps	%%xmm13,0x10(%%r15)	\n\t"\
		"addpd	    (%%rdi),%%xmm4	\n\t		mulpd	0x10(%%rax),%%xmm14	\n\t"\
		"addpd	0x10(%%rdi),%%xmm5	\n\t		mulpd	0x10(%%rax),%%xmm15	\n\t"\
		"subpd	%%xmm6,%%xmm4		\n\t		subpd	    (%%r10),%%xmm12	\n\t"\
		"subpd	%%xmm7,%%xmm5		\n\t		subpd	0x10(%%r10),%%xmm13	\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t		mulpd	    (%%rax),%%xmm12	\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t		mulpd	    (%%rax),%%xmm13	\n\t"\
		"addpd	%%xmm4,%%xmm6		\n\t		addpd	    (%%r15),%%xmm12	\n\t"\
		"addpd	%%xmm5,%%xmm7		\n\t		addpd	0x10(%%r15),%%xmm13	\n\t"\
		"movaps	%%xmm4,    (%%rsi)	\n\t		subpd	%%xmm14,%%xmm12		\n\t"\
		"movaps	%%xmm5,0x10(%%rsi)	\n\t		subpd	%%xmm15,%%xmm13		\n\t"\
		"movaps	%%xmm0,%%xmm4		\n\t		addpd	%%xmm14,%%xmm14		\n\t"\
		"movaps	%%xmm1,%%xmm5		\n\t		addpd	%%xmm15,%%xmm15		\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t		addpd	%%xmm12,%%xmm14		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t		addpd	%%xmm13,%%xmm15		\n\t"\
		"mulpd	0x20(%%rax),%%xmm0	\n\t		movaps	%%xmm12,    (%%r10)	\n\t"\
		"mulpd	0x20(%%rax),%%xmm1	\n\t		movaps	%%xmm13,0x10(%%r10)	\n\t"\
		"mulpd	0x30(%%rax),%%xmm2	\n\t		movaps	%%xmm8 ,%%xmm12		\n\t"\
		"mulpd	0x30(%%rax),%%xmm3	\n\t		movaps	%%xmm9 ,%%xmm13		\n\t"\
		"mulpd	0x40(%%rax),%%xmm4	\n\t		subpd	%%xmm10,%%xmm8 		\n\t"\
		"mulpd	0x40(%%rax),%%xmm5	\n\t		subpd	%%xmm11,%%xmm9 		\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t		mulpd	0x20(%%rax),%%xmm8 	\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t		mulpd	0x20(%%rax),%%xmm9 	\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t		mulpd	0x30(%%rax),%%xmm10	\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t		mulpd	0x30(%%rax),%%xmm11	\n\t"\
		"movaps	    (%%rsi),%%xmm4	\n\t		mulpd	0x40(%%rax),%%xmm12	\n\t"\
		"movaps	0x10(%%rsi),%%xmm5	\n\t		mulpd	0x40(%%rax),%%xmm13	\n\t"\
		"/*------ End of shared-trig-data accessed via rax register ------*/\n\t"\
		"movq	%[__o1],%%rax		\n\t		addpd	%%xmm8 ,%%xmm10		\n\t"\
		"movq	%[__o4],%%rdx		\n\t		addpd	%%xmm9 ,%%xmm11		\n\t"\
		"subpd	%%xmm3,%%xmm6		\n\t		subpd	%%xmm12,%%xmm8 		\n\t"\
		"subpd	%%xmm2,%%xmm7		\n\t		subpd	%%xmm13,%%xmm9 		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t		movaps	    (%%r10),%%xmm12	\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t		movaps	0x10(%%r10),%%xmm13	\n\t"\
		"movaps	%%xmm6,    (%%rax)	\n\t		movq	%[__o6],%%r11		\n\t"\
		"movaps	%%xmm7,0x10(%%rdx)	\n\t		movq	%[__o9],%%r14		\n\t"\
		"addpd	%%xmm6,%%xmm3		\n\t		subpd	%%xmm11,%%xmm14		\n\t"\
		"addpd	%%xmm7,%%xmm2		\n\t		subpd	%%xmm10,%%xmm15		\n\t"\
		"movaps	%%xmm3,    (%%rdx)	\n\t		addpd	%%xmm11,%%xmm11		\n\t"\
		"movaps	%%xmm2,0x10(%%rax)	\n\t		addpd	%%xmm10,%%xmm10		\n\t"\
		"movq	%[__o2],%%rbx		\n\t		movaps	%%xmm14,    (%%r11)	\n\t"\
		"movq	%[__o3],%%rcx		\n\t		movaps	%%xmm15,0x10(%%r14)	\n\t"\
		"subpd	%%xmm1,%%xmm4		\n\t		addpd	%%xmm14,%%xmm11		\n\t"\
		"subpd	%%xmm0,%%xmm5		\n\t		addpd	%%xmm15,%%xmm10		\n\t"\
		"addpd	%%xmm1,%%xmm1		\n\t		movaps	%%xmm11,    (%%r14)	\n\t"\
		"addpd	%%xmm0,%%xmm0		\n\t		movaps	%%xmm10,0x10(%%r11)	\n\t"\
		"movaps	%%xmm4,    (%%rbx)	\n\t		movq	%[__o7],%%r12		\n\t"\
		"movaps	%%xmm5,0x10(%%rcx)	\n\t		movq	%[__o8],%%r13		\n\t"\
		"addpd	%%xmm4,%%xmm1		\n\t		subpd	%%xmm9 ,%%xmm12		\n\t"\
		"addpd	%%xmm5,%%xmm0		\n\t		subpd	%%xmm8 ,%%xmm13		\n\t"\
		"movaps	%%xmm1,    (%%rcx)	\n\t		addpd	%%xmm9 ,%%xmm9 		\n\t"\
		"movaps	%%xmm0,0x10(%%rbx)	\n\t		addpd	%%xmm8 ,%%xmm8 		\n\t"\
		"										movaps	%%xmm12,    (%%r12)	\n\t"\
		"										movaps	%%xmm13,0x10(%%r13)	\n\t"\
		"										addpd	%%xmm12,%%xmm9 		\n\t"\
		"										addpd	%%xmm13,%%xmm8 		\n\t"\
		"										movaps	%%xmm9 ,    (%%r13)	\n\t"\
		"										movaps	%%xmm8 ,0x10(%%r12)	\n\t"\
		"							\n\t"\
		:					/* outputs: none */\
		: [__cc1] "m" (Xcc1)	/* All inputs from memory addresses here */\
		 ,[__i0] "m" (Xi0)\
		 ,[__i1] "m" (Xi1)\
		 ,[__i2] "m" (Xi2)\
		 ,[__i3] "m" (Xi3)\
		 ,[__i4] "m" (Xi4)\
		 ,[__o0] "m" (Xo0)\
		 ,[__o1] "m" (Xo1)\
		 ,[__o2] "m" (Xo2)\
		 ,[__o3] "m" (Xo3)\
		 ,[__o4] "m" (Xo4)\
		 ,[__i5] "m" (Yi0)\
		 ,[__i6] "m" (Yi1)\
		 ,[__i7] "m" (Yi2)\
		 ,[__i8] "m" (Yi3)\
		 ,[__i9] "m" (Yi4)\
		 ,[__o5] "m" (Yo0)\
		 ,[__o6] "m" (Yo1)\
		 ,[__o7] "m" (Yo2)\
		 ,[__o8] "m" (Yo3)\
		 ,[__o9] "m" (Yo4)\
		: "cc","memory","rax","rbx","rcx","rdx","rsi","rdi","r10","r11","r12","r13","r14","r15","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"		/* Clobbered registers */\
	);\
	}

	/*...Radix-7 DFT: Inputs in memlocs __i0-6, outputs into __o0-6, possibly coincident with inputs:\ */\
	#define SSE2_RADIX_07_DFT(Xi0,Xi1,Xi2,Xi3,Xi4,Xi5,Xi6, Xcc, Xo0,Xo1,Xo2,Xo3,Xo4,Xo5,Xo6)\
	{\
	__asm__ volatile (\
		"movq	%[__i1],%%rax		\n\t"\
		"movq	%[__i2],%%rbx		\n\t"\
		"movq	%[__i3],%%rcx		\n\t"\
		"movq	%[__i4],%%rdx		\n\t"\
		"movq	%[__i5],%%rsi		\n\t"\
		"movq	%[__i6],%%rdi		\n\t	/*** Imaginary Parts: ***/	\n\t"\
		"movaps	(%%rax),%%xmm6		\n\t	movaps	0x10(%%rax),%%xmm14	\n\t"\
		"movaps	(%%rdi),%%xmm1		\n\t	movaps	0x10(%%rdi),%%xmm9 	\n\t"\
		"movaps	(%%rbx),%%xmm5		\n\t	movaps	0x10(%%rbx),%%xmm13	\n\t"\
		"movaps	(%%rsi),%%xmm2		\n\t	movaps	0x10(%%rsi),%%xmm10	\n\t"\
		"movaps	(%%rcx),%%xmm4		\n\t	movaps	0x10(%%rcx),%%xmm12	\n\t"\
		"movaps	(%%rdx),%%xmm3		\n\t	movaps	0x10(%%rdx),%%xmm11	\n\t"\
		"movq	%[__i0],%%rbx		\n\t"\
		"subpd	%%xmm1,%%xmm6		\n\t	subpd	%%xmm9 ,%%xmm14		\n\t"\
		"addpd	%%xmm1,%%xmm1		\n\t	addpd	%%xmm9 ,%%xmm9 		\n\t"\
		"addpd	%%xmm6,%%xmm1		\n\t	addpd	%%xmm14,%%xmm9  	\n\t"\
		"subpd	%%xmm2,%%xmm5		\n\t	subpd	%%xmm10,%%xmm13		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t	addpd	%%xmm10,%%xmm10		\n\t"\
		"addpd	%%xmm5,%%xmm2		\n\t	addpd	%%xmm13,%%xmm10		\n\t"\
		"movaps	(%%rbx),%%xmm0		\n\t	movaps	0x10(%%rbx),%%xmm8 	\n\t"\
		"subpd	%%xmm3,%%xmm4		\n\t	subpd	%%xmm11,%%xmm12		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t	addpd	%%xmm11,%%xmm11		\n\t"\
		"addpd	%%xmm4,%%xmm3		\n\t	addpd	%%xmm12,%%xmm11		\n\t"\
		"\n\t"\
		"movq	%[__o0],%%rcx		\n\t"\
		"movq	%[__cc],%%rsi		\n\t"\
		"movaps	%%xmm0,0x80(%%rsi)	\n\t	movaps	%%xmm8 ,0xa0(%%rsi)	\n\t"\
		"movaps	%%xmm6,0x90(%%rsi)	\n\t	movaps	%%xmm14,0xb0(%%rsi)	\n\t"\
		"addpd	%%xmm1,%%xmm0		\n\t	addpd	%%xmm9 ,%%xmm8  	\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t	movaps	%%xmm13,%%xmm15		\n\t"\
		"addpd	%%xmm2,%%xmm3		\n\t	addpd	%%xmm10,%%xmm11		\n\t"\
		"subpd	%%xmm4,%%xmm5		\n\t	subpd	%%xmm12,%%xmm13		\n\t"\
		"subpd	%%xmm2,%%xmm1		\n\t	subpd	%%xmm10,%%xmm9  	\n\t"\
		"subpd	%%xmm7,%%xmm6		\n\t	subpd	%%xmm15,%%xmm14		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t	addpd	%%xmm10,%%xmm10		\n\t"\
		"addpd	%%xmm7,%%xmm4		\n\t	addpd	%%xmm15,%%xmm12		\n\t"\
		"addpd	%%xmm3,%%xmm0		\n\t	addpd	%%xmm11,%%xmm8  	\n\t"\
		"addpd	0x90(%%rsi),%%xmm5	\n\t	addpd	0xb0(%%rsi),%%xmm13	\n\t"\
		"subpd	%%xmm2,%%xmm3		\n\t	subpd	%%xmm10,%%xmm11		\n\t"\
		"movaps	%%xmm4,%%xmm7		\n\t	movaps	%%xmm12,%%xmm15		\n\t"\
		"movaps	%%xmm0,    (%%rcx)	\n\t	movaps	%%xmm8 ,0x10(%%rcx)	\n\t"/* B0 */\
		"subpd	%%xmm6,%%xmm4		\n\t	subpd	%%xmm14,%%xmm12		\n\t"\
		"movaps	%%xmm1,%%xmm2		\n\t	movaps	%%xmm9 ,%%xmm10		\n\t"\
		"subpd	0x80(%%rsi),%%xmm0	\n\t	subpd	0xa0(%%rsi),%%xmm8  \n\t"\
		"mulpd	0x10(%%rsi),%%xmm5	\n\t	mulpd	0x10(%%rsi),%%xmm13	\n\t"\
		"addpd	%%xmm3,%%xmm2		\n\t	addpd	%%xmm11,%%xmm10		\n\t"\
		"mulpd	0x40(%%rsi),%%xmm3	\n\t	mulpd	0x40(%%rsi),%%xmm11	\n\t"\
		"mulpd	0x70(%%rsi),%%xmm4	\n\t	mulpd	0x70(%%rsi),%%xmm12	\n\t"\
		"mulpd	0x20(%%rsi),%%xmm1	\n\t	mulpd	0x20(%%rsi),%%xmm9  \n\t"\
		"mulpd	0x30(%%rsi),%%xmm6	\n\t	mulpd	0x30(%%rsi),%%xmm14	\n\t"\
		"mulpd	    (%%rsi),%%xmm0	\n\t	mulpd	    (%%rsi),%%xmm8  \n\t"\
		"mulpd	0x50(%%rsi),%%xmm7	\n\t	mulpd	0x50(%%rsi),%%xmm15	\n\t"\
		"mulpd	0x60(%%rsi),%%xmm2	\n\t	mulpd	0x60(%%rsi),%%xmm10	\n\t"\
		"addpd	    (%%rcx),%%xmm0	\n\t	addpd	0x10(%%rcx),%%xmm8  \n\t"\
		"addpd	%%xmm4,%%xmm6		\n\t	addpd	%%xmm12,%%xmm14		\n\t"\
		"subpd	%%xmm2,%%xmm1		\n\t	subpd	%%xmm10,%%xmm9  	\n\t"\
		"subpd	%%xmm7,%%xmm4		\n\t	subpd	%%xmm15,%%xmm12		\n\t"\
		"subpd	%%xmm2,%%xmm3		\n\t	subpd	%%xmm10,%%xmm11		\n\t"\
		"movq	%[__o1],%%rax		\n\t"\
		"movq	%[__o2],%%rbx		\n\t"\
		"movq	%[__o3],%%rcx		\n\t"\
		"movq	%[__o4],%%rdx		\n\t"\
		"movq	%[__o5],%%rsi		\n\t"\
		"movq	%[__o6],%%rdi		\n\t"\
		"movaps	%%xmm0,%%xmm2		\n\t	movaps	%%xmm8 ,%%xmm10		\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t	movaps	%%xmm13,%%xmm15		\n\t"\
		"addpd	%%xmm1,%%xmm0		\n\t	addpd	%%xmm9 ,%%xmm8  	\n\t"\
		"addpd	%%xmm6,%%xmm5		\n\t	addpd	%%xmm14,%%xmm13		\n\t"\
		"addpd	%%xmm3,%%xmm1		\n\t	addpd	%%xmm11,%%xmm9  	\n\t"\
		"addpd	%%xmm4,%%xmm6		\n\t	addpd	%%xmm12,%%xmm14		\n\t"\
		"addpd	%%xmm2,%%xmm3		\n\t	addpd	%%xmm10,%%xmm11		\n\t"\
		"addpd	%%xmm7,%%xmm4		\n\t	addpd	%%xmm15,%%xmm12		\n\t"\
		"subpd	%%xmm1,%%xmm2		\n\t	subpd	%%xmm9 ,%%xmm10		\n\t"\
		"subpd	%%xmm6,%%xmm7		\n\t	subpd	%%xmm14,%%xmm15		\n\t"\
		"/* xmm1,6,9,14 free ... Note the order reversal on the 3rd pair of outputs: */\n\t"\
		"subpd	%%xmm13,%%xmm0		\n\t	subpd	%%xmm15,%%xmm2		\n\t	subpd	%%xmm12,%%xmm3 		\n\t"\
		"subpd	%%xmm5 ,%%xmm8  	\n\t	subpd	%%xmm7 ,%%xmm10		\n\t	subpd	%%xmm4 ,%%xmm11		\n\t"\
		"addpd	%%xmm13,%%xmm13		\n\t	addpd	%%xmm15,%%xmm15		\n\t	addpd	%%xmm12,%%xmm12		\n\t"\
		"addpd	%%xmm5 ,%%xmm5		\n\t	addpd	%%xmm7 ,%%xmm7		\n\t	addpd	%%xmm4 ,%%xmm4 		\n\t"\
		"addpd	%%xmm0 ,%%xmm13		\n\t	addpd	%%xmm2 ,%%xmm15		\n\t	addpd	%%xmm3 ,%%xmm12		\n\t"\
		"addpd	%%xmm8 ,%%xmm5		\n\t	addpd	%%xmm10,%%xmm7		\n\t	addpd	%%xmm11,%%xmm4 		\n\t"\
		"movaps	%%xmm0 ,    (%%rax)	\n\t	movaps	%%xmm2 ,    (%%rbx)	\n\t	movaps	%%xmm3 ,    (%%rdx)	\n\t"/* B124r */\
		"movaps	%%xmm8 ,0x10(%%rdi)	\n\t	movaps	%%xmm10,0x10(%%rsi)	\n\t	movaps	%%xmm11,0x10(%%rcx)	\n\t"/* B653i */\
		"movaps	%%xmm13,    (%%rdi)	\n\t	movaps	%%xmm15,    (%%rsi)	\n\t	movaps	%%xmm12,    (%%rcx)	\n\t"/* B653r */\
		"movaps	%%xmm5 ,0x10(%%rax)	\n\t	movaps	%%xmm7 ,0x10(%%rbx)	\n\t	movaps	%%xmm4 ,0x10(%%rdx)	\n\t"/* B124i */\
		"\n\t"\
		:					/* outputs: none */\
		: [__i0] "m" (Xi0)	/* All inputs from memory addresses here */\
		 ,[__i1] "m" (Xi1)\
		 ,[__i2] "m" (Xi2)\
		 ,[__i3] "m" (Xi3)\
		 ,[__i4] "m" (Xi4)\
		 ,[__i5] "m" (Xi5)\
		 ,[__i6] "m" (Xi6)\
		 ,[__cc] "m" (Xcc)\
		 ,[__o0] "m" (Xo0)\
		 ,[__o1] "m" (Xo1)\
		 ,[__o2] "m" (Xo2)\
		 ,[__o3] "m" (Xo3)\
		 ,[__o4] "m" (Xo4)\
		 ,[__o5] "m" (Xo5)\
		 ,[__o6] "m" (Xo6)\
		: "cc","memory","rax","rbx","rcx","rdx","rsi","rdi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"		/* Clobbered registers */\
	);\
	}

	/* Twiddleless version of SSE2_RADIX8_DIF_TWIDDLE. Inputs enter in memory locations __r0 + [__i1,__i2,__i3,__i4,__i5,__i6,__i7],;
	where r0 is a memory address and the i's are LITERAL [BYTE] OFFSETS. Outputs go into memory locations __o0,__o1,__o2,__o3,__o4,__o5,__o6,__o7, assumed disjoint with inputs:\
	*/
	#define SSE2_RADIX8_DIF_0TWIDDLE(Xr0, Xi1,Xi2,Xi3,Xi4,Xi5,Xi6,Xi7, Xo0,Xo1,Xo2,Xo3,Xo4,Xo5,Xo6,Xo7, Xisrt2)\
	{\
	__asm__ volatile (\
		"/* 1st of 2 radix-4 subtransforms, data in xmm0-7: */\n\t	/* 2nd of 2 radix-4 subtransforms, data in xmm8-15: */\n\t"\
		"movq	%[__r0],%%rax	/* i0 = r00 */	\n\t			movq	%[__i1],%%r10		/* i1 */	\n\t"\
		"movq	%[__i2],%%rbx	/* i2 */		\n\t			movq	%[__i3],%%r11		/* i3 */	\n\t"\
		"movq	%[__i4],%%rcx	/* i4 */		\n\t			movq	%[__i5],%%r12		/* i5 */	\n\t"\
		"movq	%[__i6],%%rdx	/* i6 */		\n\t			movq	%[__i7],%%r13		/* i7 */	\n\t"\
		"addq	%%rax,%%rbx						\n\t			addq	%%rax,%%r10						\n\t"\
		"addq	%%rax,%%rcx						\n\t			addq	%%rax,%%r11						\n\t"\
		"addq	%%rax,%%rdx						\n\t			addq	%%rax,%%r12						\n\t"\
		"movq	%[__isrt2],%%rsi				\n\t			addq	%%rax,%%r13						\n\t"\
		"										\n\t			/* p1,5 combo: x+y into xmm8 /1, x-y in xmm10/3: */	\n\t"\
		"/* p0,4 combo: x+y into xmm0/1, x-y in xmm2/3: */\n\t	movaps	    (%%r12),%%xmm8 				\n\t"\
		"										\n\t			movaps	0x10(%%r12),%%xmm9 				\n\t"\
		"movaps	    (%%rcx),%%xmm0				\n\t			movaps	    (%%r10),%%xmm10				\n\t"\
		"movaps	0x10(%%rcx),%%xmm1				\n\t			movaps	0x10(%%r10),%%xmm11				\n\t"\
		"movaps	    (%%rax),%%xmm2				\n\t			subpd	%%xmm8 ,%%xmm10					\n\t"\
		"movaps	0x10(%%rax),%%xmm3				\n\t			subpd	%%xmm9 ,%%xmm11					\n\t"\
		"subpd	%%xmm0,%%xmm2					\n\t			addpd	%%xmm8 ,%%xmm8 					\n\t"\
		"subpd	%%xmm1,%%xmm3					\n\t			addpd	%%xmm9 ,%%xmm9 					\n\t"\
		"addpd	%%xmm0,%%xmm0					\n\t			addpd	%%xmm10,%%xmm8 					\n\t"\
		"addpd	%%xmm1,%%xmm1					\n\t			addpd	%%xmm11,%%xmm9 					\n\t"\
		"addpd	%%xmm2,%%xmm0					\n\t			/* p3,7 combo: x+y into xmm14/7, x-y in xmm12/5: */	\n\t"\
		"addpd	%%xmm3,%%xmm1					\n\t			movaps	    (%%r11),%%xmm12				\n\t"\
		"										\n\t			movaps	0x10(%%r11),%%xmm13				\n\t"\
		"/* p2,6 combo: x+y into xmm4/5, x-y in xmm6/7: */\n\t	movaps	    (%%r13),%%xmm14				\n\t"\
		"										\n\t			movaps	0x10(%%r13),%%xmm15				\n\t"\
		"movaps	    (%%rdx),%%xmm4				\n\t			subpd	%%xmm14,%%xmm12					\n\t"\
		"movaps	0x10(%%rdx),%%xmm5				\n\t			subpd	%%xmm15,%%xmm13					\n\t"\
		"movaps	    (%%rbx),%%xmm6				\n\t			addpd	%%xmm14,%%xmm14					\n\t"\
		"movaps	0x10(%%rbx),%%xmm7				\n\t			addpd	%%xmm15,%%xmm15					\n\t"\
		"subpd	%%xmm4,%%xmm6					\n\t			addpd	%%xmm12,%%xmm14					\n\t"\
		"subpd	%%xmm5,%%xmm7					\n\t			addpd	%%xmm13,%%xmm15					\n\t"\
		"addpd	%%xmm4,%%xmm4					\n\t			/* Finish radix-4 butterfly, tmp-store 1st of 4 outputs to free up 2 registers: */\n\t"\
		"addpd	%%xmm5,%%xmm5					\n\t			subpd	%%xmm14,%%xmm8 					\n\t"\
		"addpd	%%xmm6,%%xmm4					\n\t			subpd	%%xmm15,%%xmm9 					\n\t"\
		"addpd	%%xmm7,%%xmm5					\n\t			subpd	%%xmm13,%%xmm10					\n\t"\
		"										\n\t			subpd	%%xmm12,%%xmm11					\n\t"\
		"subpd	%%xmm4,%%xmm0					\n\t			addpd	%%xmm14,%%xmm14					\n\t"\
		"subpd	%%xmm7,%%xmm2					\n\t			addpd	%%xmm13,%%xmm13					\n\t"\
		"subpd	%%xmm5,%%xmm1					\n\t			addpd	%%xmm15,%%xmm15					\n\t"\
		"subpd	%%xmm6,%%xmm3					\n\t			addpd	%%xmm12,%%xmm12					\n\t"\
		"														addpd	%%xmm8 ,%%xmm14					\n\t"\
		"														addpd	%%xmm10,%%xmm13					\n\t"\
		"														addpd	%%xmm9 ,%%xmm15					\n\t"\
		"														addpd	%%xmm11,%%xmm12					\n\t"\
		"addpd	%%xmm4,%%xmm4					\n\t			movaps	%%xmm14,    (%%r10)				\n\t"\
		"addpd	%%xmm7,%%xmm7					\n\t			movaps	%%xmm15,0x10(%%r10)				\n\t"\
		"addpd	%%xmm5,%%xmm5					\n\t			movaps	%%xmm10,%%xmm14					\n\t"\
		"addpd	%%xmm6,%%xmm6					\n\t			movaps	%%xmm13,%%xmm15					\n\t"\
		"addpd	%%xmm0,%%xmm4					\n\t			subpd	%%xmm12,%%xmm10					\n\t"\
		"addpd	%%xmm2,%%xmm7					\n\t			subpd	%%xmm11,%%xmm13					\n\t"\
		"addpd	%%xmm1,%%xmm5					\n\t			addpd	%%xmm14,%%xmm12					\n\t"\
		"addpd	%%xmm3,%%xmm6					\n\t			addpd	%%xmm15,%%xmm11					\n\t"\
		"														movaps	(%%rsi),%%xmm14	/* isrt2 */		\n\t"\
		"														mulpd	%%xmm14,%%xmm10					\n\t"\
		"														mulpd	%%xmm14,%%xmm13					\n\t"\
		"														mulpd	%%xmm14,%%xmm12					\n\t"\
		"										\n\t			mulpd	%%xmm14,%%xmm11					\n\t"\
		"movaps	    (%%r10),%%xmm14	/* restore spilled */\n\t"\
		"movaps	0x10(%%r10),%%xmm15	/* restore spilled */\n\t"\
		"										\n\t"\
		"/* Inline of SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_A(r0): Combine radix-4 subtransforms and write outputs: */\n\t"\
		"/***** t0,1,2,3,4,5,6,7 in xmm[ 4, 5| 2, 6| 0, 1| 7, 3] *****/\n\t"\
		"/***** t8,9,a,b,c,d,e,f in xmm[14,15|10,12| 8, 9|13,11] */\n\t"\
		"movq	%[__o4],%%rax					\n\t			subpd   %%xmm10,%%xmm2					\n\t"\
		"movq	%[__o5],%%rbx					\n\t			subpd   %%xmm12,%%xmm6					\n\t"\
		"movq	%[__o6],%%rcx					\n\t			addpd   %%xmm10,%%xmm10					\n\t"\
		"movq	%[__o7],%%rdx					\n\t			addpd   %%xmm12,%%xmm12					\n\t"\
		"										\n\t			addpd   %%xmm2,%%xmm10					\n\t"\
		"subpd   %%xmm11,%%xmm7					\n\t			addpd   %%xmm6,%%xmm12					\n\t"\
		"subpd   %%xmm13,%%xmm3					\n\t													\n\t"\
		"addpd   %%xmm11,%%xmm11				\n\t			movaps	%%xmm2 ,    (%%rbx)	/* o5r */	\n\t"\
		"addpd   %%xmm13,%%xmm13				\n\t			movaps	%%xmm6 ,0x10(%%rbx)	/* o5i */	\n\t"\
		"addpd   %%xmm7,%%xmm11					\n\t			movaps	%%xmm10,    (%%rax)	/* o4r */	\n\t"\
		"addpd   %%xmm3,%%xmm13					\n\t			movaps	%%xmm12,0x10(%%rax)	/* o4i */	\n\t"\
		"										\n\t"\
		"movaps	%%xmm7 ,    (%%rcx)	/* o6r */	\n\t"\
		"movaps	%%xmm3 ,0x10(%%rdx)	/* o7i */	\n\t"\
		"movaps	%%xmm11,    (%%rdx)	/* o7r */	\n\t"\
		"movaps	%%xmm13,0x10(%%rcx)	/* o6i */	\n\t"\
		"										\n\t"\
		"movq	%[__o0],%%rax					\n\t"\
		"movq	%[__o1],%%rbx					\n\t"\
		"movq	%[__o2],%%rcx					\n\t"\
		"movq	%[__o3],%%rdx					\n\t"\
		"										\n\t"\
		"subpd	%%xmm14,%%xmm4 					\n\t"\
		"subpd	%%xmm15,%%xmm5 					\n\t"\
		"subpd	%%xmm9 ,%%xmm0 					\n\t"\
		"subpd	%%xmm8 ,%%xmm1 					\n\t"\
		"addpd	%%xmm14,%%xmm14					\n\t			movaps	%%xmm4 ,    (%%rbx)	/* o1r */	\n\t"\
		"addpd	%%xmm15,%%xmm15					\n\t			movaps	%%xmm5 ,0x10(%%rbx)	/* o1i */	\n\t"\
		"addpd	%%xmm9 ,%%xmm9 					\n\t			movaps	%%xmm0 ,    (%%rcx)	/* o2r */	\n\t"\
		"addpd	%%xmm8 ,%%xmm8 					\n\t			movaps	%%xmm1 ,0x10(%%rdx)	/* o3i */	\n\t"\
		"addpd	%%xmm4 ,%%xmm14					\n\t"\
		"addpd	%%xmm5 ,%%xmm15					\n\t"\
		"addpd	%%xmm0 ,%%xmm9 					\n\t"\
		"addpd	%%xmm1 ,%%xmm8 					\n\t"\
		"										\n\t"\
		"movaps	%%xmm14,    (%%rax)	/* o0r */	\n\t"\
		"movaps	%%xmm15,0x10(%%rax)	/* o0r */	\n\t"\
		"movaps	%%xmm9 ,    (%%rdx)	/* o3r */	\n\t"\
		"movaps	%%xmm8 ,0x10(%%rcx)	/* o2i */	\n\t"\
		"										\n\t"\
		:					/* outputs: none */\
		: [__r0] "m" (Xr0)	/* All inputs from memory addresses here */\
		 ,[__i1] "e" (Xi1)\
		 ,[__i2] "e" (Xi2)\
		 ,[__i3] "e" (Xi3)\
		 ,[__i4] "e" (Xi4)\
		 ,[__i5] "e" (Xi5)\
		 ,[__i6] "e" (Xi6)\
		 ,[__i7] "e" (Xi7)\
		 ,[__o0] "m" (Xo0)\
		 ,[__o1] "m" (Xo1)\
		 ,[__o2] "m" (Xo2)\
		 ,[__o3] "m" (Xo3)\
		 ,[__o4] "m" (Xo4)\
		 ,[__o5] "m" (Xo5)\
		 ,[__o6] "m" (Xo6)\
		 ,[__o7] "m" (Xo7)\
		 ,[__isrt2] "m" (Xisrt2)\
		: "cc","memory","rax","rbx","rcx","rdx","rsi","r10","r11","r12","r13","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"		/* Clobbered registers */\
	);\
	}

	/* Twiddleless version of SSE2_RADIX8_DIT_TWIDDLE. Inputs enter in memory locations __i0,__i1,__i2,__i3,__i4,__i5,__i6,__i7.
	Outputs go into 16 contiguous 32-byte memory locations starting at __out and assumed disjoint with inputs.
	This macro built on the same code template as SSE2_RADIX8_DIF_TWIDDLE0, but with the I/O-location indices mutually bit reversed:
	01234567 <--> 04261537, which can be effected via the pairwise swaps 1 <--> 4 and 3 <--> 6.
	*/
	#define	SSE2_RADIX8_DIT_0TWIDDLE(Xi0,Xi1,Xi2,Xi3,Xi4,Xi5,Xi6,Xi7, Xout, Xisrt2)\
	{\
	__asm__ volatile (\
		"/* 1st of 2 radix-4 subtransforms, data in xmm0-7: */\n\t	/* 2nd of 2 radix-4 subtransforms, data in xmm8-15: */\n\t"\
		"movq	%[__i0],%%rax					\n\t			movq	%[__i4],%%r10					\n\t"\
		"movq	%[__i1],%%rbx					\n\t			movq	%[__i5],%%r11					\n\t"\
		"movq	%[__i2],%%rcx					\n\t			movq	%[__i6],%%r12					\n\t"\
		"movq	%[__i3],%%rdx					\n\t			movq	%[__i7],%%r13					\n\t"\
		"										\n\t			/* p1,5 combo: x+y into xmm8 /1, x-y in xmm10/3: */	\n\t"\
		"/* p0,4 combo: x+y into xmm0/1, x-y in xmm2/3: */\n\t	movaps	    (%%r11),%%xmm8 				\n\t"\
		"										\n\t			movaps	0x10(%%r11),%%xmm9 				\n\t"\
		"movaps	    (%%rbx),%%xmm0				\n\t			movaps	    (%%r10),%%xmm10				\n\t"\
		"movaps	0x10(%%rbx),%%xmm1				\n\t			movaps	0x10(%%r10),%%xmm11				\n\t"\
		"movaps	    (%%rax),%%xmm2				\n\t			subpd	%%xmm8 ,%%xmm10					\n\t"\
		"movaps	0x10(%%rax),%%xmm3				\n\t			subpd	%%xmm9 ,%%xmm11					\n\t"\
		"subpd	%%xmm0,%%xmm2					\n\t			addpd	%%xmm8 ,%%xmm8 					\n\t"\
		"subpd	%%xmm1,%%xmm3					\n\t			addpd	%%xmm9 ,%%xmm9 					\n\t"\
		"addpd	%%xmm0,%%xmm0					\n\t			addpd	%%xmm10,%%xmm8 					\n\t"\
		"addpd	%%xmm1,%%xmm1					\n\t			addpd	%%xmm11,%%xmm9 					\n\t"\
		"addpd	%%xmm2,%%xmm0					\n\t			/* p3,7 combo: x+y into xmm14/7, x-y in xmm12/5: */	\n\t"\
		"addpd	%%xmm3,%%xmm1					\n\t			movaps	    (%%r12),%%xmm12				\n\t"\
		"										\n\t			movaps	0x10(%%r12),%%xmm13				\n\t"\
		"/* p2,6 combo: x+y into xmm4/5, x-y in xmm6/7: */\n\t	movaps	    (%%r13),%%xmm14				\n\t"\
		"										\n\t			movaps	0x10(%%r13),%%xmm15				\n\t"\
		"movaps	    (%%rdx),%%xmm4				\n\t			subpd	%%xmm14,%%xmm12					\n\t"\
		"movaps	0x10(%%rdx),%%xmm5				\n\t			subpd	%%xmm15,%%xmm13					\n\t"\
		"movaps	    (%%rcx),%%xmm6				\n\t			addpd	%%xmm14,%%xmm14					\n\t"\
		"movaps	0x10(%%rcx),%%xmm7				\n\t			addpd	%%xmm15,%%xmm15					\n\t"\
		"subpd	%%xmm4,%%xmm6					\n\t			addpd	%%xmm12,%%xmm14					\n\t"\
		"subpd	%%xmm5,%%xmm7					\n\t			addpd	%%xmm13,%%xmm15					\n\t"\
		"addpd	%%xmm4,%%xmm4					\n\t			/* Finish radix-4 butterfly, tmp-store 1st of 4 outputs to free up 2 registers: */\n\t"\
		"addpd	%%xmm5,%%xmm5					\n\t			subpd	%%xmm14,%%xmm8 					\n\t"\
		"addpd	%%xmm6,%%xmm4					\n\t			subpd	%%xmm15,%%xmm9 					\n\t"\
		"addpd	%%xmm7,%%xmm5					\n\t			subpd	%%xmm13,%%xmm10					\n\t"\
		"										\n\t			subpd	%%xmm12,%%xmm11					\n\t"\
		"subpd	%%xmm4,%%xmm0					\n\t			addpd	%%xmm14,%%xmm14					\n\t"\
		"subpd	%%xmm7,%%xmm2					\n\t			addpd	%%xmm13,%%xmm13					\n\t"\
		"subpd	%%xmm5,%%xmm1					\n\t			addpd	%%xmm15,%%xmm15					\n\t"\
		"subpd	%%xmm6,%%xmm3					\n\t			addpd	%%xmm12,%%xmm12					\n\t"\
		"														addpd	%%xmm8 ,%%xmm14					\n\t"\
		"														addpd	%%xmm10,%%xmm13					\n\t"\
		"														addpd	%%xmm9 ,%%xmm15					\n\t"\
		"														addpd	%%xmm11,%%xmm12					\n\t"\
		"														movq	%[__isrt2],%%rsi	/* isrt2 */	\n\t"\
		"addpd	%%xmm4,%%xmm4					\n\t			movaps	%%xmm14,    (%%rax)	/* spill */	\n\t"\
		"addpd	%%xmm7,%%xmm7					\n\t			movaps	%%xmm15,0x10(%%rax)	/* spill */	\n\t"\
		"addpd	%%xmm5,%%xmm5					\n\t			movaps	%%xmm10,%%xmm14					\n\t"\
		"addpd	%%xmm6,%%xmm6					\n\t			movaps	%%xmm13,%%xmm15					\n\t"\
		"addpd	%%xmm0,%%xmm4					\n\t			subpd	%%xmm12,%%xmm10					\n\t"\
		"addpd	%%xmm2,%%xmm7					\n\t			subpd	%%xmm11,%%xmm13					\n\t"\
		"addpd	%%xmm1,%%xmm5					\n\t			addpd	%%xmm14,%%xmm12					\n\t"\
		"addpd	%%xmm3,%%xmm6					\n\t			addpd	%%xmm15,%%xmm11					\n\t"\
		"														movaps	(%%rsi),%%xmm14		/* isrt2 */	\n\t"\
		"														mulpd	%%xmm14,%%xmm10					\n\t"\
		"														mulpd	%%xmm14,%%xmm13					\n\t"\
		"														mulpd	%%xmm14,%%xmm12					\n\t"\
		"														mulpd	%%xmm14,%%xmm11					\n\t"\
		"/* Combine radix-4 subtransforms and write outputs: */\n\t"\
		"\n\t"\
		"movaps	    (%%rax),%%xmm14	/* restore spilled */\n\t	subpd   %%xmm10,%%xmm2					\n\t"\
		"movaps	0x10(%%rax),%%xmm15	/* restore spilled */\n\t	subpd   %%xmm12,%%xmm6					\n\t"\
		"														addpd   %%xmm10,%%xmm10					\n\t"\
		"movq	%[__out],%%rax					\n\t			addpd   %%xmm12,%%xmm12					\n\t"\
		"										\n\t			addpd   %%xmm2,%%xmm10					\n\t"\
		"subpd   %%xmm11,%%xmm7					\n\t			addpd   %%xmm6,%%xmm12					\n\t"\
		"subpd   %%xmm13,%%xmm3					\n\t													\n\t"\
		"addpd   %%xmm11,%%xmm11				\n\t			movaps	%%xmm2 ,0xa0(%%rax)	/* o5r */	\n\t"\
		"addpd   %%xmm13,%%xmm13				\n\t			movaps	%%xmm6 ,0xb0(%%rax)	/* o5i */	\n\t"\
		"addpd   %%xmm7,%%xmm11					\n\t			movaps	%%xmm10,0x20(%%rax)	/* o1r */	\n\t"\
		"addpd   %%xmm3,%%xmm13					\n\t			movaps	%%xmm12,0x30(%%rax)	/* o1i */	\n\t"\
		"										\n\t"\
		"movaps	%%xmm7 ,0x60(%%rax)	/* o3r */	\n\t"\
		"movaps	%%xmm3 ,0xf0(%%rax)	/* o7i */	\n\t"\
		"movaps	%%xmm11,0xe0(%%rax)	/* o7r */	\n\t"\
		"movaps	%%xmm13,0x70(%%rax)	/* o3i */	\n\t"\
		"										\n\t"\
		"subpd	%%xmm14,%%xmm4 					\n\t"\
		"subpd	%%xmm15,%%xmm5 					\n\t"\
		"subpd	%%xmm9 ,%%xmm0 					\n\t"\
		"subpd	%%xmm8 ,%%xmm1 					\n\t"\
		"addpd	%%xmm14,%%xmm14					\n\t			movaps	%%xmm4 ,0x80(%%rax)	/* o4r */	\n\t"\
		"addpd	%%xmm15,%%xmm15					\n\t			movaps	%%xmm5 ,0x90(%%rax)	/* o4i */	\n\t"\
		"addpd	%%xmm9 ,%%xmm9 					\n\t			movaps	%%xmm0 ,0x40(%%rax)	/* o2r */	\n\t"\
		"addpd	%%xmm8 ,%%xmm8 					\n\t			movaps	%%xmm1 ,0xd0(%%rax)	/* o6i */	\n\t"\
		"addpd	%%xmm4 ,%%xmm14					\n\t"\
		"addpd	%%xmm5 ,%%xmm15					\n\t"\
		"addpd	%%xmm0 ,%%xmm9 					\n\t"\
		"addpd	%%xmm1 ,%%xmm8 					\n\t"\
		"										\n\t"\
		"movaps	%%xmm14,    (%%rax)	/* o0r */	\n\t"\
		"movaps	%%xmm15,0x10(%%rax)	/* o0r */	\n\t"\
		"movaps	%%xmm9 ,0xc0(%%rax)	/* o6r */	\n\t"\
		"movaps	%%xmm8 ,0x50(%%rax)	/* o2i */	\n\t"\
		"										\n\t"\
		:					/* outputs: none */\
		: [__i0] "m" (Xi0)	/* All iputs from memory addresses here */\
		 ,[__i1] "m" (Xi1)\
		 ,[__i2] "m" (Xi2)\
		 ,[__i3] "m" (Xi3)\
		 ,[__i4] "m" (Xi4)\
		 ,[__i5] "m" (Xi5)\
		 ,[__i6] "m" (Xi6)\
		 ,[__i7] "m" (Xi7)\
		 ,[__out] "m" (Xout)\
		 ,[__isrt2] "m" (Xisrt2)\
		: "cc","memory","rax","rbx","rcx","rdx","rsi","r10","r11","r12","r13","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"		/* Clobbered registers */\
	);\
	}

	// Same as SSE2_RADIX8_DIT_0TWIDDLE but with user-specifiable [i.e. not nec. contiguous] output addresses:
	#define	SSE2_RADIX8_DIT_0TWIDDLE_OOP(Xi0,Xi1,Xi2,Xi3,Xi4,Xi5,Xi6,Xi7, Xo0,Xo1,Xo2,Xo3,Xo4,Xo5,Xo6,Xo7, Xisrt2)\
	{\
	__asm__ volatile (\
		"/* 1st of 2 radix-4 subtransforms, data in xmm0-7: */\n\t	/* 2nd of 2 radix-4 subtransforms, data in xmm8-15: */\n\t"\
		"movq	%[__i0],%%rax					\n\t			movq	%[__i4],%%r10					\n\t"\
		"movq	%[__i1],%%rbx					\n\t			movq	%[__i5],%%r11					\n\t"\
		"movq	%[__i2],%%rcx					\n\t			movq	%[__i6],%%r12					\n\t"\
		"movq	%[__i3],%%rdx					\n\t			movq	%[__i7],%%r13					\n\t"\
		"										\n\t			/* p1,5 combo: x+y into xmm8 /1, x-y in xmm10/3: */	\n\t"\
		"/* p0,4 combo: x+y into xmm0/1, x-y in xmm2/3: */\n\t	movaps	    (%%r11),%%xmm8 				\n\t"\
		"										\n\t			movaps	0x10(%%r11),%%xmm9 				\n\t"\
		"movaps	    (%%rbx),%%xmm0				\n\t			movaps	    (%%r10),%%xmm10				\n\t"\
		"movaps	0x10(%%rbx),%%xmm1				\n\t			movaps	0x10(%%r10),%%xmm11				\n\t"\
		"movaps	    (%%rax),%%xmm2				\n\t			subpd	%%xmm8 ,%%xmm10					\n\t"\
		"movaps	0x10(%%rax),%%xmm3				\n\t			subpd	%%xmm9 ,%%xmm11					\n\t"\
		"subpd	%%xmm0,%%xmm2					\n\t			addpd	%%xmm8 ,%%xmm8 					\n\t"\
		"subpd	%%xmm1,%%xmm3					\n\t			addpd	%%xmm9 ,%%xmm9 					\n\t"\
		"addpd	%%xmm0,%%xmm0					\n\t			addpd	%%xmm10,%%xmm8 					\n\t"\
		"addpd	%%xmm1,%%xmm1					\n\t			addpd	%%xmm11,%%xmm9 					\n\t"\
		"addpd	%%xmm2,%%xmm0					\n\t			/* p3,7 combo: x+y into xmm14/7, x-y in xmm12/5: */	\n\t"\
		"addpd	%%xmm3,%%xmm1					\n\t			movaps	    (%%r12),%%xmm12				\n\t"\
		"										\n\t			movaps	0x10(%%r12),%%xmm13				\n\t"\
		"/* p2,6 combo: x+y into xmm4/5, x-y in xmm6/7: */\n\t	movaps	    (%%r13),%%xmm14				\n\t"\
		"										\n\t			movaps	0x10(%%r13),%%xmm15				\n\t"\
		"movaps	    (%%rdx),%%xmm4				\n\t			subpd	%%xmm14,%%xmm12					\n\t"\
		"movaps	0x10(%%rdx),%%xmm5				\n\t			subpd	%%xmm15,%%xmm13					\n\t"\
		"movaps	    (%%rcx),%%xmm6				\n\t			addpd	%%xmm14,%%xmm14					\n\t"\
		"movaps	0x10(%%rcx),%%xmm7				\n\t			addpd	%%xmm15,%%xmm15					\n\t"\
		"subpd	%%xmm4,%%xmm6					\n\t			addpd	%%xmm12,%%xmm14					\n\t"\
		"subpd	%%xmm5,%%xmm7					\n\t			addpd	%%xmm13,%%xmm15					\n\t"\
		"addpd	%%xmm4,%%xmm4					\n\t			/* Finish radix-4 butterfly, tmp-store 1st of 4 outputs to free up 2 registers: */\n\t"\
		"addpd	%%xmm5,%%xmm5					\n\t			subpd	%%xmm14,%%xmm8 					\n\t"\
		"addpd	%%xmm6,%%xmm4					\n\t			subpd	%%xmm15,%%xmm9 					\n\t"\
		"addpd	%%xmm7,%%xmm5					\n\t			subpd	%%xmm13,%%xmm10					\n\t"\
		"										\n\t			subpd	%%xmm12,%%xmm11					\n\t"\
		"subpd	%%xmm4,%%xmm0					\n\t			addpd	%%xmm14,%%xmm14					\n\t"\
		"subpd	%%xmm7,%%xmm2					\n\t			addpd	%%xmm13,%%xmm13					\n\t"\
		"subpd	%%xmm5,%%xmm1					\n\t			addpd	%%xmm15,%%xmm15					\n\t"\
		"subpd	%%xmm6,%%xmm3					\n\t			addpd	%%xmm12,%%xmm12					\n\t"\
		"														addpd	%%xmm8 ,%%xmm14					\n\t"\
		"														addpd	%%xmm10,%%xmm13					\n\t"\
		"														addpd	%%xmm9 ,%%xmm15					\n\t"\
		"														addpd	%%xmm11,%%xmm12					\n\t"\
		"														movq	%[__isrt2],%%rsi	/* isrt2 */	\n\t"\
		"addpd	%%xmm4,%%xmm4					\n\t			movaps	%%xmm14,    (%%rax)	/* spill */	\n\t"\
		"addpd	%%xmm7,%%xmm7					\n\t			movaps	%%xmm15,0x10(%%rax)	/* spill */	\n\t"\
		"addpd	%%xmm5,%%xmm5					\n\t			movaps	%%xmm10,%%xmm14					\n\t"\
		"addpd	%%xmm6,%%xmm6					\n\t			movaps	%%xmm13,%%xmm15					\n\t"\
		"addpd	%%xmm0,%%xmm4					\n\t			subpd	%%xmm12,%%xmm10					\n\t"\
		"addpd	%%xmm2,%%xmm7					\n\t			subpd	%%xmm11,%%xmm13					\n\t"\
		"addpd	%%xmm1,%%xmm5					\n\t			addpd	%%xmm14,%%xmm12					\n\t"\
		"addpd	%%xmm3,%%xmm6					\n\t			addpd	%%xmm15,%%xmm11					\n\t"\
		"														movaps	(%%rsi),%%xmm14		/* isrt2 */	\n\t"\
		"														mulpd	%%xmm14,%%xmm10					\n\t"\
		"														mulpd	%%xmm14,%%xmm13					\n\t"\
		"														mulpd	%%xmm14,%%xmm12					\n\t"\
		"														mulpd	%%xmm14,%%xmm11					\n\t"\
		"/* Combine radix-4 subtransforms and write outputs: */\n\t"\
		"\n\t"\
		"movaps	    (%%rax),%%xmm14	/* restore spilled */\n\t	subpd   %%xmm10,%%xmm2					\n\t"\
		"movaps	0x10(%%rax),%%xmm15	/* restore spilled */\n\t	subpd   %%xmm12,%%xmm6					\n\t"\
		"movq	%[__o1],%%rax					\n\t			movq	%[__o5],%%rcx					\n\t"\
		"														addpd   %%xmm10,%%xmm10					\n\t"\
		"										\n\t			addpd   %%xmm12,%%xmm12					\n\t"\
		"										\n\t			addpd   %%xmm2,%%xmm10					\n\t"\
		"subpd   %%xmm11,%%xmm7					\n\t			addpd   %%xmm6,%%xmm12					\n\t"\
		"subpd   %%xmm13,%%xmm3					\n\t													\n\t"\
		"movq	%[__o3],%%rbx					\n\t			movq	%[__o7],%%rdx					\n\t"\
		"addpd   %%xmm11,%%xmm11				\n\t			movaps	%%xmm2 ,    (%%rcx)	/* o5r */	\n\t"\
		"addpd   %%xmm13,%%xmm13				\n\t			movaps	%%xmm6 ,0x10(%%rcx)	/* o5i */	\n\t"\
		"addpd   %%xmm7,%%xmm11					\n\t			movaps	%%xmm10,    (%%rax)	/* o1r */	\n\t"\
		"addpd   %%xmm3,%%xmm13					\n\t			movaps	%%xmm12,0x10(%%rax)	/* o1i */	\n\t"\
		"movq	%[__o0],%%rax					\n\t			movq	%[__o4],%%rcx					\n\t"\
		"										\n\t"\
		"movaps	%%xmm7 ,    (%%rbx)	/* o3r */	\n\t"\
		"movaps	%%xmm3 ,0x10(%%rdx)	/* o7i */	\n\t"\
		"movaps	%%xmm11,    (%%rdx)	/* o7r */	\n\t"\
		"movaps	%%xmm13,0x10(%%rbx)	/* o3i */	\n\t"\
		"										\n\t"\
		"movq	%[__o2],%%rbx					\n\t			movq	%[__o6],%%rdx					\n\t"\
		"subpd	%%xmm14,%%xmm4 					\n\t"\
		"subpd	%%xmm15,%%xmm5 					\n\t"\
		"subpd	%%xmm9 ,%%xmm0 					\n\t"\
		"subpd	%%xmm8 ,%%xmm1 					\n\t"\
		"addpd	%%xmm14,%%xmm14					\n\t			movaps	%%xmm4 ,    (%%rcx)	/* o4r */	\n\t"\
		"addpd	%%xmm15,%%xmm15					\n\t			movaps	%%xmm5 ,0x10(%%rcx)	/* o4i */	\n\t"\
		"addpd	%%xmm9 ,%%xmm9 					\n\t			movaps	%%xmm0 ,    (%%rbx)	/* o2r */	\n\t"\
		"addpd	%%xmm8 ,%%xmm8 					\n\t			movaps	%%xmm1 ,0x10(%%rdx)	/* o6i */	\n\t"\
		"addpd	%%xmm4 ,%%xmm14					\n\t"\
		"addpd	%%xmm5 ,%%xmm15					\n\t"\
		"addpd	%%xmm0 ,%%xmm9 					\n\t"\
		"addpd	%%xmm1 ,%%xmm8 					\n\t"\
		"										\n\t"\
		"movaps	%%xmm14,    (%%rax)	/* o0r */	\n\t"\
		"movaps	%%xmm15,0x10(%%rax)	/* o0i */	\n\t"\
		"movaps	%%xmm9 ,    (%%rdx)	/* o6r */	\n\t"\
		"movaps	%%xmm8 ,0x10(%%rbx)	/* o2i */	\n\t"\
		"										\n\t"\
		:					/* outputs: none */\
		: [__i0] "m" (Xi0)	/* All iputs from memory addresses here */\
		 ,[__i1] "m" (Xi1)\
		 ,[__i2] "m" (Xi2)\
		 ,[__i3] "m" (Xi3)\
		 ,[__i4] "m" (Xi4)\
		 ,[__i5] "m" (Xi5)\
		 ,[__i6] "m" (Xi6)\
		 ,[__i7] "m" (Xi7)\
		 ,[__o0] "m" (Xo0)\
		 ,[__o1] "m" (Xo1)\
		 ,[__o2] "m" (Xo2)\
		 ,[__o3] "m" (Xo3)\
		 ,[__o4] "m" (Xo4)\
		 ,[__o5] "m" (Xo5)\
		 ,[__o6] "m" (Xo6)\
		 ,[__o7] "m" (Xo7)\
		 ,[__isrt2] "m" (Xisrt2)\
		: "cc","memory","rax","rbx","rcx","rdx","rsi","r10","r11","r12","r13","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"		/* Clobbered registers */\
	);\
	}

	// SSE2 analog of dft_macro.h::RADIX_08_DIF_TWIDDLE_OOP - Result of adding separate I/O addressing to
	// radix8_dif_dit_pass_gcc64.h::SSE2_RADIX8_DIF_TWIDDLE:
	#define SSE2_RADIX8_DIF_TWIDDLE_OOP(Xi0,Xi1,Xi2,Xi3,Xi4,Xi5,Xi6,Xi7, Xo0,Xo1,Xo2,Xo3,Xo4,Xo5,Xo6,Xo7 ,Xc1,Xs1,Xc2,Xs2,Xc3,Xs3,Xc4,Xs4,Xc5,Xs5,Xc6,Xs6,Xc7,Xs7)\
	{\
	__asm__ volatile (\
		/* Block 0,4: */\
		"movq		%[i0]	,%%rax			\n\t"\
		"movq		%[i4]	,%%rbx			\n\t"\
		"movq		%[c4]	,%%rcx			\n\t"\
		"movq		%[s4]	,%%rsi			\n\t"\
	/* [rsi] (and if needed rdi) points to sine components of each sincos pair, which is not really a pair here in terms of relative addressing: */\
		"movaps		    (%%rbx)	,%%xmm2		\n\t"\
		"movaps		0x10(%%rbx)	,%%xmm3		\n\t"\
		"movaps		(%%rcx)	,%%xmm0		\n\t"\
		"movaps		(%%rsi)	,%%xmm1		\n\t"\
		"movaps		%%xmm2	,%%xmm4		\n\t"\
		"movaps		%%xmm3	,%%xmm5		\n\t"\
		"mulpd		%%xmm0	,%%xmm2		\n\t"\
		"mulpd		%%xmm0	,%%xmm3		\n\t"\
		"mulpd		%%xmm1	,%%xmm4		\n\t"\
		"mulpd		%%xmm1	,%%xmm5		\n\t"\
		"movaps		    (%%rax)	,%%xmm0		\n\t"\
		"movaps		0x10(%%rax)	,%%xmm1		\n\t"\
		"subpd		%%xmm5		,%%xmm2		\n\t"\
		"addpd		%%xmm4		,%%xmm3		\n\t"\
		"movaps		%%xmm0	,%%xmm6		\n\t"\
		"movaps		%%xmm1	,%%xmm7		\n\t"\
		"addpd		%%xmm2		,%%xmm0		\n\t"\
		"addpd		%%xmm3		,%%xmm1		\n\t"\
		"subpd		%%xmm2		,%%xmm6		\n\t"\
		"subpd		%%xmm3		,%%xmm7		\n\t"\
		"movaps		%%xmm0		,    (%%rax)	\n\t"\
		"movaps		%%xmm1		,0x10(%%rax)	\n\t"\
		"movaps		%%xmm6		,    (%%rbx)	\n\t"\
		"movaps		%%xmm7		,0x10(%%rbx)	\n\t"\
		/* Block 2,6: */\
		"movq		%[i2]	,%%rax			\n\t"\
		"movq		%[i6]	,%%rbx			\n\t"\
		"movq		%[c2]	,%%rcx			\n\t"\
		"movq		%[c6]	,%%rdx			\n\t"\
		"movq		%[s2]	,%%rsi			\n\t"\
		"movq		%[s6]	,%%rdi			\n\t"\
		"movaps		    (%%rax)	,%%xmm0		\n\t"\
		"movaps		0x10(%%rax)	,%%xmm2		\n\t"\
		"movaps		(%%rcx)	,%%xmm6		\n\t"\
		"movaps		(%%rsi)	,%%xmm7		\n\t"\
		"movaps		%%xmm0	,%%xmm1		\n\t"\
		"movaps		%%xmm2	,%%xmm3		\n\t"\
		"mulpd		%%xmm6	,%%xmm0		\n\t"\
		"mulpd		%%xmm7	,%%xmm2		\n\t"\
		"mulpd		%%xmm7	,%%xmm1		\n\t"\
		"mulpd		%%xmm6	,%%xmm3		\n\t"\
		"subpd		%%xmm2		,%%xmm0		\n\t"\
		"addpd		%%xmm3		,%%xmm1		\n\t"\
		"movaps		    (%%rbx)	,%%xmm2		\n\t"\
		"movaps		0x10(%%rbx)	,%%xmm3		\n\t"\
		"movaps		(%%rdx)	,%%xmm6		\n\t"\
		"movaps		(%%rdi)	,%%xmm7		\n\t"\
		"movaps		%%xmm2	,%%xmm4		\n\t"\
		"movaps		%%xmm3	,%%xmm5		\n\t"\
		"mulpd		%%xmm6	,%%xmm2		\n\t"\
		"mulpd		%%xmm7	,%%xmm3		\n\t"\
		"mulpd		%%xmm7	,%%xmm4		\n\t"\
		"mulpd		%%xmm6	,%%xmm5		\n\t"\
		"subpd		%%xmm3		,%%xmm2		\n\t"\
		"addpd		%%xmm5		,%%xmm4		\n\t"\
		"movaps		%%xmm2		,%%xmm3		\n\t"\
		"movaps		%%xmm4		,%%xmm5		\n\t"\
		"addpd		%%xmm0		,%%xmm2		\n\t"\
		"subpd		%%xmm3		,%%xmm0		\n\t"\
		"addpd		%%xmm1		,%%xmm4		\n\t"\
		"subpd		%%xmm5		,%%xmm1		\n\t"\
		"movaps		%%xmm2		,    (%%rax)	\n\t"\
		"movaps		%%xmm4		,0x10(%%rax)	\n\t"\
		"movaps		%%xmm0		,    (%%rbx)	\n\t"\
		"movaps		%%xmm1		,0x10(%%rbx)	\n\t"\
		/* Block 1,5: */\
		"movq		%[i1]	,%%rax			\n\t"\
		"movq		%[i5]	,%%rbx			\n\t"\
		"movq		%[c1]	,%%rcx			\n\t"\
		"movq		%[c5]	,%%rdx			\n\t"\
		"movq		%[s1]	,%%rsi			\n\t"\
		"movq		%[s5]	,%%rdi			\n\t"\
		"movaps		    (%%rax)	,%%xmm0		\n\t"\
		"movaps		0x10(%%rax)	,%%xmm2		\n\t"\
		"movaps		(%%rcx)	,%%xmm6		\n\t"\
		"movaps		(%%rsi)	,%%xmm7		\n\t"\
		"movaps		%%xmm0	,%%xmm1		\n\t"\
		"movaps		%%xmm2	,%%xmm3		\n\t"\
		"mulpd		%%xmm6	,%%xmm0		\n\t"\
		"mulpd		%%xmm7	,%%xmm2		\n\t"\
		"mulpd		%%xmm7	,%%xmm1		\n\t"\
		"mulpd		%%xmm6	,%%xmm3		\n\t"\
		"subpd		%%xmm2		,%%xmm0		\n\t"\
		"addpd		%%xmm3		,%%xmm1		\n\t"\
		"movaps		    (%%rbx)	,%%xmm2		\n\t"\
		"movaps		0x10(%%rbx)	,%%xmm3		\n\t"\
		"movaps		(%%rdx)	,%%xmm6		\n\t"\
		"movaps		(%%rdi)	,%%xmm7		\n\t"\
		"movaps		%%xmm2	,%%xmm4		\n\t"\
		"movaps		%%xmm3	,%%xmm5		\n\t"\
		"mulpd		%%xmm6	,%%xmm2		\n\t"\
		"mulpd		%%xmm7	,%%xmm3		\n\t"\
		"mulpd		%%xmm7	,%%xmm4		\n\t"\
		"mulpd		%%xmm6	,%%xmm5		\n\t"\
		"subpd		%%xmm3		,%%xmm2		\n\t"\
		"addpd		%%xmm5		,%%xmm4		\n\t"\
		"movaps		%%xmm2		,%%xmm3		\n\t"\
		"movaps		%%xmm4		,%%xmm5		\n\t"\
		"addpd		%%xmm0		,%%xmm2		\n\t"\
		"subpd		%%xmm3		,%%xmm0		\n\t"\
		"addpd		%%xmm1		,%%xmm4		\n\t"\
		"subpd		%%xmm5		,%%xmm1		\n\t"\
		"movaps		%%xmm2		,    (%%rax)	\n\t"\
		"movaps		%%xmm4		,0x10(%%rax)	\n\t"\
		"movaps		%%xmm0		,    (%%rbx)	\n\t"\
		"movaps		%%xmm1		,0x10(%%rbx)	\n\t"\
		/* Block 3,7: */\
		"movq		%[i3]	,%%rax			\n\t"\
		"movq		%[i7]	,%%rbx			\n\t"\
		"movq		%[c3]	,%%rcx			\n\t"\
		"movq		%[c7]	,%%rdx			\n\t"\
		"movq		%[s3]	,%%rsi			\n\t"\
		"movq		%[s7]	,%%rdi			\n\t"\
		"movaps		    (%%rax)	,%%xmm0		\n\t"\
		"movaps		0x10(%%rax)	,%%xmm2		\n\t"\
		"movaps		(%%rcx)	,%%xmm6		\n\t"\
		"movaps		(%%rsi)	,%%xmm7		\n\t"\
		"movaps		%%xmm0	,%%xmm1		\n\t"\
		"movaps		%%xmm2	,%%xmm3		\n\t"\
		"mulpd		%%xmm6	,%%xmm0		\n\t"\
		"mulpd		%%xmm7	,%%xmm2		\n\t"\
		"mulpd		%%xmm7	,%%xmm1		\n\t"\
		"mulpd		%%xmm6	,%%xmm3		\n\t"\
		"subpd		%%xmm2		,%%xmm0		\n\t"\
		"addpd		%%xmm3		,%%xmm1		\n\t"\
		"movaps		    (%%rbx)	,%%xmm2		\n\t"\
		"movaps		0x10(%%rbx)	,%%xmm3		\n\t"\
		"movaps		(%%rdx)	,%%xmm6		\n\t"\
		"movaps		(%%rdi)	,%%xmm7		\n\t"\
		"movaps		%%xmm2	,%%xmm4		\n\t"\
		"movaps		%%xmm3	,%%xmm5		\n\t"\
		"mulpd		%%xmm6	,%%xmm2		\n\t"\
		"mulpd		%%xmm7	,%%xmm3		\n\t"\
		"mulpd		%%xmm7	,%%xmm4		\n\t"\
		"mulpd		%%xmm6	,%%xmm5		\n\t"\
		"subpd		%%xmm3		,%%xmm2		\n\t"\
		"addpd		%%xmm5		,%%xmm4		\n\t"\
		"movaps		%%xmm2		,%%xmm3		\n\t"\
		"movaps		%%xmm4		,%%xmm5		\n\t"\
		"addpd		%%xmm0		,%%xmm2		\n\t"\
		"subpd		%%xmm3		,%%xmm0		\n\t"\
		"addpd		%%xmm1		,%%xmm4		\n\t"\
		"subpd		%%xmm5		,%%xmm1		\n\t"\
		"movaps		%%xmm2		,    (%%rax)	\n\t"\
		"movaps		%%xmm4		,0x10(%%rax)	\n\t"\
		"movaps		%%xmm0		,    (%%rbx)	\n\t"\
		"movaps		%%xmm1		,0x10(%%rbx)	\n\t"\
		/* Now combine the 2 half-transforms: */\
		/* Combo 1: */\
		"movq		%[i0]	,%%rax			\n\t"\
		"movq		%[i2]	,%%rbx			\n\t"\
		"movaps		    (%%rax)	,%%xmm0		\n\t"\
		"movaps		0x10(%%rax)	,%%xmm1		\n\t"\
		"movaps		    (%%rbx)	,%%xmm2		\n\t"\
		"movaps		0x10(%%rbx)	,%%xmm3		\n\t"\
		"movaps		%%xmm0		,%%xmm4		\n\t"\
		"movaps		%%xmm1		,%%xmm5		\n\t"\
		"addpd		%%xmm2	,%%xmm0		\n\t"\
		"subpd		%%xmm2	,%%xmm4		\n\t"\
		"addpd		%%xmm3	,%%xmm1		\n\t"\
		"subpd		%%xmm3	,%%xmm5		\n\t"\
		"movaps		%%xmm0		,    (%%rax)	\n\t"\
		"movaps		%%xmm1		,0x10(%%rax)	\n\t"\
		"movaps		%%xmm4		,    (%%rbx)	\n\t"\
		"movaps		%%xmm5		,0x10(%%rbx)	\n\t"\
		"movq		%[i1]	,%%rcx			\n\t"\
		"movq		%[i3]	,%%rdx			\n\t"\
		"movaps		    (%%rcx)	,%%xmm2		\n\t"\
		"movaps		0x10(%%rcx)	,%%xmm3		\n\t"\
		"movaps		    (%%rdx)	,%%xmm8		\n\t"/* Use xtra 64-bit reg xmm8 here */\
		"movaps		0x10(%%rdx)	,%%xmm9		\n\t"/* Use xtra 64-bit reg xmm9 here */\
		"movq		%[o1]	,%%rcx			\n\t"\
		"movq		%[o3]	,%%rdx			\n\t"\
		"movaps		%%xmm2		,%%xmm6		\n\t"\
		"movaps		%%xmm3		,%%xmm7		\n\t"\
		"addpd		%%xmm8	,%%xmm2		\n\t"\
		"subpd		%%xmm8	,%%xmm6		\n\t"\
		"addpd		%%xmm9	,%%xmm3		\n\t"\
		"subpd		%%xmm9	,%%xmm7		\n\t"\
		"subpd		%%xmm2		,%%xmm0		\n\t"\
		"subpd		%%xmm3		,%%xmm1		\n\t"\
		"subpd		%%xmm7		,%%xmm4		\n\t"\
		"subpd		%%xmm6		,%%xmm5		\n\t"\
		"addpd			(%%rax)	,%%xmm2		\n\t"\
		"addpd		0x10(%%rax)	,%%xmm3		\n\t"\
		"addpd			(%%rbx)	,%%xmm7		\n\t"\
		"addpd		0x10(%%rbx)	,%%xmm6		\n\t"\
		"movq		%[o0]	,%%rax			\n\t"\
		"movq		%[o2]	,%%rbx			\n\t"\
		"movaps		%%xmm2		,    (%%rax)	\n\t"/* [o0].re */\
		"movaps		%%xmm3		,0x10(%%rax)	\n\t"/* [o0].im */\
		"movaps		%%xmm4		,    (%%rbx)	\n\t"/* [o2].re */\
		"movaps		%%xmm6		,0x10(%%rbx)	\n\t"/* [o2].im */\
		"movaps		%%xmm0		,    (%%rcx)	\n\t"/* [o1].re */\
		"movaps		%%xmm1		,0x10(%%rcx)	\n\t"/* [o1].im */\
		"movaps		%%xmm7		,    (%%rdx)	\n\t"/* [o3].re */\
		"movaps		%%xmm5		,0x10(%%rdx)	\n\t"/* [o3].im */\
		/* Combo 2: */\
		"movq		%[i4]	,%%rax			\n\t"\
		"movq		%[i6]	,%%rbx			\n\t"\
		"movaps		    (%%rax)	,%%xmm0		\n\t"\
		"movaps		0x10(%%rax)	,%%xmm1		\n\t"\
		"movaps		    (%%rbx)	,%%xmm2		\n\t"\
		"movaps		0x10(%%rbx)	,%%xmm3		\n\t"\
		"movaps		%%xmm0		,%%xmm4		\n\t"\
		"movaps		%%xmm1		,%%xmm5		\n\t"\
		"subpd		%%xmm3	,%%xmm0		\n\t"\
		"addpd		%%xmm3	,%%xmm4		\n\t"\
		"addpd		%%xmm2	,%%xmm1		\n\t"\
		"subpd		%%xmm2	,%%xmm5		\n\t"\
		"movaps		%%xmm0		,    (%%rax)	\n\t"\
		"movaps		%%xmm1		,0x10(%%rax)	\n\t"\
		"movaps		%%xmm4		,    (%%rbx)	\n\t"\
		"movaps		%%xmm5		,0x10(%%rbx)	\n\t"\
		"movq		%[i5]	,%%rcx			\n\t"\
		"movq		%[i7]	,%%rdx			\n\t"\
		"movaps		    (%%rcx)	,%%xmm2		\n\t"\
		"movaps		0x10(%%rcx)	,%%xmm3		\n\t"\
		"movaps		    (%%rdx)	,%%xmm8		\n\t"/* Use xtra 64-bit reg xmm8 here */\
		"movaps		0x10(%%rdx)	,%%xmm9		\n\t"/* Use xtra 64-bit reg xmm9 here */\
		"movaps		%%xmm2		,%%xmm6		\n\t"\
		"movaps		%%xmm3		,%%xmm7		\n\t"\
		"subpd		%%xmm9	,%%xmm2		\n\t"\
		"addpd		%%xmm9	,%%xmm6		\n\t"\
		"addpd		%%xmm8	,%%xmm3		\n\t"\
		"subpd		%%xmm8	,%%xmm7		\n\t"\
		"movaps		%%xmm4		,    (%%rdx)	\n\t"\
		"movaps		%%xmm5		,0x10(%%rdx)	\n\t"\
	/* Use the cosine term of the [c1,s1] pair, which is the *middle* [4th of 7] of our 7 input pairs, in terms \
	of the input-arg bit-reversal reordering defined in the __X[c,s] --> [c,s] mapping below and happens to \
	always in fact *be* a true cosine term, which is a requirement for our "decr 1 gives isrt2" data-copy scheme: */\
		"movq		%[c1],%%rsi			\n\t"/* isrt2 in [c1]-1 */\
		"movaps	-0x10(%%rsi),%%xmm8		\n\t"/* Use xtra 64-bit reg xmm8 here */\
		"movaps		%%xmm2		,%%xmm5		\n\t"\
		"subpd		%%xmm3		,%%xmm2		\n\t"\
		"addpd		%%xmm3		,%%xmm5		\n\t"\
		"mulpd			%%xmm8	,%%xmm2		\n\t"\
		"mulpd			%%xmm8	,%%xmm5		\n\t"\
		"movaps		0x10(%%rdx)	,%%xmm3		\n\t"\
		"movaps		%%xmm7		,%%xmm4		\n\t"\
		"addpd		%%xmm6		,%%xmm4		\n\t"\
		"subpd		%%xmm6		,%%xmm7		\n\t"\
		"mulpd			%%xmm8	,%%xmm4		\n\t"\
		"mulpd			%%xmm8	,%%xmm7		\n\t"\
		"movaps			(%%rdx)	,%%xmm6		\n\t"/* last ref to rdx-for-in-address */\
		"movq		%[o5]	,%%rcx			\n\t"\
		"movq		%[o7]	,%%rdx			\n\t"\
		"subpd		%%xmm2		,%%xmm0		\n\t"\
		"subpd		%%xmm5		,%%xmm1		\n\t"\
		"subpd		%%xmm4		,%%xmm6		\n\t"\
		"subpd		%%xmm7		,%%xmm3		\n\t"\
		"addpd			(%%rax)	,%%xmm2		\n\t"\
		"addpd		0x10(%%rax)	,%%xmm5		\n\t"\
		"addpd			(%%rbx)	,%%xmm4		\n\t"\
		"addpd		0x10(%%rbx)	,%%xmm7		\n\t"\
		"movq		%[o4]	,%%rax			\n\t"\
		"movq		%[o6]	,%%rbx			\n\t"\
		"movaps		%%xmm2		,    (%%rax)	\n\t"/* [o4].re */\
		"movaps		%%xmm5		,0x10(%%rax)	\n\t"/* [o4].im */\
		"movaps		%%xmm6		,    (%%rbx)	\n\t"/* [o6].re */\
		"movaps		%%xmm3		,0x10(%%rbx)	\n\t"/* [o6].im */\
		"movaps		%%xmm0		,    (%%rcx)	\n\t"/* [o5].re */\
		"movaps		%%xmm1		,0x10(%%rcx)	\n\t"/* [o5].im */\
		"movaps		%%xmm4		,    (%%rdx)	\n\t"/* [o7].re */\
		"movaps		%%xmm7		,0x10(%%rdx)	\n\t"/* [o7].im */\
		:					/* outputs: none */\
		: [i0] "m" (Xi0)	/* All inputs from memory addresses here */\
		 ,[i1] "m" (Xi1)\
		 ,[i2] "m" (Xi2)\
		 ,[i3] "m" (Xi3)\
		 ,[i4] "m" (Xi4)\
		 ,[i5] "m" (Xi5)\
		 ,[i6] "m" (Xi6)\
		 ,[i7] "m" (Xi7)\
		 ,[o0] "m" (Xo0)\
		 ,[o1] "m" (Xo1)\
		 ,[o2] "m" (Xo2)\
		 ,[o3] "m" (Xo3)\
		 ,[o4] "m" (Xo4)\
		 ,[o5] "m" (Xo5)\
		 ,[o6] "m" (Xo6)\
		 ,[o7] "m" (Xo7)\
		 ,[c4] "m" (Xc1),[s4] "m" (Xs1)\
		 ,[c2] "m" (Xc2),[s2] "m" (Xs2)\
		 ,[c6] "m" (Xc3),[s6] "m" (Xs3)\
		 ,[c1] "m" (Xc4),[s1] "m" (Xs4)\
		 ,[c5] "m" (Xc5),[s5] "m" (Xs5)\
		 ,[c3] "m" (Xc6),[s3] "m" (Xs6)\
		 ,[c7] "m" (Xc7),[s7] "m" (Xs7)\
		: "cc","memory","rax","rbx","rcx","rdx","rsi","rdi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
	);\
	}

	// SSE2 analog of dft_macro.h::RADIX_08_DIT_TWIDDLE_OOP - Result of sign-flippage and adding separate I/O addressing to
	// radix8_dif_dit_pass_gcc64.h::SSE2_RADIX8_DIF_TWIDDLE. We begin with the DIF macro here because we need a pre-twiddles
	// implementation for our purposes, whereas SSE2_RADIX8_DIT_TWIDDLE is post-twiddles.
	//
	// SIMD Opcount: 102 load/store [30 implicit], 66 add/sub, 50 mul. Compare to DFT macros used for radix-8-pass-with-twiddles:
	// DIF opcount : 140 load/store [56 implicit], 66 add/sub, 32 mul
	// DIT opcount :  85 load/store [36 implicit], 68 add/sub, 32 mul . 
	//
	#define SSE2_RADIX8_DIT_TWIDDLE_OOP(Xi0,Xi1,Xi2,Xi3,Xi4,Xi5,Xi6,Xi7, Xo0,Xo1,Xo2,Xo3,Xo4,Xo5,Xo6,Xo7 ,Xc1,Xs1,Xc2,Xs2,Xc3,Xs3,Xc4,Xs4,Xc5,Xs5,Xc6,Xs6,Xc7,Xs7)\
	{\
	__asm__ volatile (\
	/* Block 0/1 has just one twiddle-CMUL: */\
		"movq		%[i0],%%rax			\n\t"\
		"movq		%[i1],%%rbx			\n\t"\
		"movq		%[c1],%%rdi			\n\t"/* [rdi,rsi] point to [cos,sin] components of each sincos pair, */\
		"movq		%[s1],%%rsi			\n\t"/* which is not really a pair here in terms of relative addressing: */\
		"movaps		    (%%rbx),%%xmm4 	\n\t	movaps		0x10(%%rbx),%%xmm5 	\n\t"/* _r4  = __tr1;	_r5  = __ti1; */\
		"movaps		    (%%rax),%%xmm0 	\n\t	movaps		0x10(%%rax),%%xmm1 	\n\t"/* _r0  = __tr0;	_r1  = __ti0; */\
		"movaps		%%xmm5 ,%%xmm6 		\n\t	movaps		%%xmm4 ,%%xmm7 		\n\t"/* _r6  = _r5;		_r7  = _r4;		** [r4,r5] = CMUL(__t1,__W1): */\
		"mulpd		(%%rdi),%%xmm4 		\n\t	mulpd		(%%rdi),%%xmm5 		\n\t"/* _r4 *= __Wr1;	_r5 *= __Wr1; */\
		"mulpd		(%%rsi),%%xmm6 		\n\t	mulpd		(%%rsi),%%xmm7 		\n\t"/* _r6 *= __Wi1;	_r7 *= __Wi1; */\
		"addpd		%%xmm6 ,%%xmm4 		\n\t	subpd		%%xmm7 ,%%xmm5 		\n\t"/* _r4 += _r6;		_r5 -= _r7; */\
		"movaps		%%xmm0 ,%%xmm2 		\n\t	movaps		%%xmm1 ,%%xmm3 		\n\t"/* _r2  = _r0;		_r3  = _r1; */\
		"addpd		%%xmm4 ,%%xmm0 		\n\t	addpd		%%xmm5 ,%%xmm1 		\n\t"/* _r0 += _r4;		_r1 += _r5; */\
		"subpd		%%xmm4 ,%%xmm2 		\n\t	subpd		%%xmm5 ,%%xmm3 		\n\t"/* _r2 -= _r4;		_r3 -= _r5; */\
		"movaps		%%xmm0 ,    (%%rax)	\n\t	movaps		%%xmm1 ,0x10(%%rax)	\n\t"/* __tr0 = _r0;	__ti0 = _r1; */\
		"movaps		%%xmm2 ,    (%%rbx)	\n\t	movaps		%%xmm3 ,0x10(%%rbx)	\n\t"/* __tr1 = _r2;	__ti1 = _r3; */\
	/* Blocks 2/3 use separate register subset, can be done overlapped with 0/1: */\
		"movq		%[i2],%%rcx			\n\t"\
		"movq		%[c2],%%r10			\n\t"\
		"movq		%[s2],%%r11			\n\t"/* [r8,r9] = CMUL(__t2,__W2): */\
		"movaps		    (%%rcx),%%xmm8 	\n\t	movaps		0x10(%%rcx),%%xmm9 	\n\t"/* _r8  = __tr2;	_r9  = __ti2; */\
		"movaps		%%xmm9 ,%%xmm10		\n\t	movaps		%%xmm8 ,%%xmm11		\n\t"/* _ra  = _r9;		_rb  = _r8; */\
		"mulpd		(%%r10),%%xmm8 		\n\t	mulpd		(%%r10),%%xmm9 		\n\t"/* _r8 *= __Wr2;	_r9 *= __Wr2; */\
		"mulpd		(%%r11),%%xmm10		\n\t	mulpd		(%%r11),%%xmm11		\n\t"/* _ra *= __Wi2;	_rb *= __Wi2; */\
		"addpd		%%xmm10,%%xmm8 		\n\t	subpd		%%xmm11,%%xmm9 		\n\t"/* _r8 += _ra;		_r9 -= _rb; */\
		"movq		%[i3],%%rdx			\n\t"\
		"movq		%[c3],%%r12			\n\t"\
		"movq		%[s3],%%r13			\n\t"/* [rc,rd] = CMUL(__t3,__W3): */\
		"movaps		    (%%rdx),%%xmm12	\n\t	movaps		0x10(%%rdx),%%xmm13	\n\t"/* _rc  = __tr3;	_rd  = __ti3; */\
		"movaps		%%xmm13,%%xmm14		\n\t	movaps		%%xmm12,%%xmm15		\n\t"/* _re  = _rd;		_rf  = _rc; */\
		"mulpd		(%%r12),%%xmm12		\n\t	mulpd		(%%r12),%%xmm13		\n\t"/* _rc *= __Wr3;	_rd *= __Wr3; */\
		"mulpd		(%%r13),%%xmm14		\n\t	mulpd		(%%r13),%%xmm15		\n\t"/* _re *= __Wi3;	_rf *= __Wi3; */\
		"addpd		%%xmm14,%%xmm12		\n\t	subpd		%%xmm15,%%xmm13		\n\t"/* _rc += _re;		_rd -= _rf; */\
		/* Now do radix-2 butterfly: */\
		"movaps		%%xmm8 ,%%xmm10		\n\t	movaps		%%xmm9 ,%%xmm11		\n\t"/* _ra  = _r8;		_rb  = _r9; */\
		"addpd		%%xmm12,%%xmm8 		\n\t	addpd		%%xmm13,%%xmm9 		\n\t"/* _r8 += _rc;		_r9 += _rd; */\
		"subpd		%%xmm12,%%xmm10		\n\t	subpd		%%xmm13,%%xmm11		\n\t"/* _ra -= _rc;		_rb -= _rd; */\
		"movaps		%%xmm8 ,    (%%rcx)	\n\t	movaps		%%xmm9 ,0x10(%%rcx)	\n\t"/* __tr2 = _r8;	__ti2 = _r9; */\
		"movaps		%%xmm10,    (%%rdx)	\n\t	movaps		%%xmm11,0x10(%%rdx)	\n\t"/* __tr3 = _ra;	__ti3 = _rb; */\
	/* Blocks 4/5: */\
		"movq		%[i4],%%rax			\n\t"\
		"movq		%[c4],%%rdi			\n\t"\
		"movq		%[s4],%%rsi			\n\t"/* [r0,r1] = CMUL(__t4,__W4): */\
		"movaps		    (%%rax),%%xmm0 	\n\t	movaps		0x10(%%rax),%%xmm1 	\n\t"/* _r0  = __tr4;	_r1  = __ti4; */\
		"movaps		%%xmm1 ,%%xmm2 		\n\t	movaps		%%xmm0 ,%%xmm3 		\n\t"/* _r2  = _r1;		_r3  = _r0; */\
		"mulpd		(%%rdi),%%xmm0 		\n\t	mulpd		(%%rdi),%%xmm1 		\n\t"/* _r0 *= __Wr4;	_r1 *= __Wr4; */\
		"mulpd		(%%rsi),%%xmm2 		\n\t	mulpd		(%%rsi),%%xmm3 		\n\t"/* _r2 *= __Wi4;	_r3 *= __Wi4; */\
		"addpd		%%xmm2 ,%%xmm0 		\n\t	subpd		%%xmm3 ,%%xmm1 		\n\t"/* _r0 += _r2;		_r1 -= _r3; */\
		"movq		%[i5],%%rbx			\n\t"\
		"movq		%[c5],%%r8 			\n\t"\
		"movq		%[s5],%%r9 			\n\t"/* [r4,r5] = CMUL(__t5,__W5): */\
		"movaps		    (%%rbx),%%xmm4 	\n\t	movaps		0x10(%%rbx),%%xmm5 	\n\t"/* _r4  = __tr5;	_r5  = __ti5; */\
		"movaps		%%xmm5 ,%%xmm6 		\n\t	movaps		%%xmm4 ,%%xmm7 		\n\t"/* _r6  = _r5;		_r7  = _r4; */\
		"mulpd		(%%r8 ),%%xmm4 		\n\t	mulpd		(%%r8 ),%%xmm5 		\n\t"/* _r4 *= __Wr5;	_r5 *= __Wr5; */\
		"mulpd		(%%r9 ),%%xmm6 		\n\t	mulpd		(%%r9 ),%%xmm7 		\n\t"/* _r6 *= __Wi5;	_r7 *= __Wi5; */\
		"addpd		%%xmm6 ,%%xmm4 		\n\t	subpd		%%xmm7 ,%%xmm5 		\n\t"/* _r4 += _r6;		_r5 -= _r7; */\
		/* Now do radix-2 butterfly: */\
		"movaps		%%xmm0 ,%%xmm2 		\n\t	movaps		%%xmm1 ,%%xmm3 		\n\t"/* _r2  = _r0;		_r3  = _r1; */\
		"addpd		%%xmm4 ,%%xmm0 		\n\t	addpd		%%xmm5 ,%%xmm1 		\n\t"/* _r0 += _r4;		_r1 += _r5; */\
		"subpd		%%xmm4 ,%%xmm2 		\n\t	subpd		%%xmm5 ,%%xmm3 		\n\t"/* _r2 -= _r4;		_r3 -= _r5; */\
	/* Blocks 6/7 use separate register subset, can be done overlapped with 4/5: */\
		"movq		%[i6],%%rcx			\n\t"\
		"movq		%[c6],%%r10			\n\t"\
		"movq		%[s6],%%r11			\n\t"/* [r8,r9] = CMUL(__t6,__W6): */\
		"movaps		    (%%rcx),%%xmm8 	\n\t	movaps		0x10(%%rcx),%%xmm9 	\n\t"/* _r8  = __tr6;	_r9  = __ti6; */\
		"movaps		%%xmm9 ,%%xmm10		\n\t	movaps		%%xmm8 ,%%xmm11		\n\t"/* _ra  = _r9;		_rb  = _r8; */\
		"mulpd		(%%r10),%%xmm8 		\n\t	mulpd		(%%r10),%%xmm9 		\n\t"/* _r8 *= __Wr6;	_r9 *= __Wr6; */\
		"mulpd		(%%r11),%%xmm10		\n\t	mulpd		(%%r11),%%xmm11		\n\t"/* _ra *= __Wi6;	_rb *= __Wi6; */\
		"addpd		%%xmm10,%%xmm8 		\n\t	subpd		%%xmm11,%%xmm9 		\n\t"/* _r8 += _ra;		_r9 -= _rb; */\
		"movq		%[i7],%%rdx			\n\t"\
		"movq		%[c7],%%r12			\n\t"\
		"movq		%[s7],%%r13			\n\t"/* [rc,rd] = CMUL(__t7,__W7): */\
		"movaps		    (%%rdx),%%xmm12	\n\t	movaps		0x10(%%rdx),%%xmm13	\n\t"/* _rc  = __tr7;	_rd  = __ti7; */\
		"movaps		%%xmm13,%%xmm14		\n\t	movaps		%%xmm12,%%xmm15		\n\t"/* _re  = _rd;		_rf  = _rc; */\
		"mulpd		(%%r12),%%xmm12		\n\t	mulpd		(%%r12),%%xmm13		\n\t"/* _rc *= __Wr7;	_rd *= __Wr7; */\
		"mulpd		(%%r13),%%xmm14		\n\t	mulpd		(%%r13),%%xmm15		\n\t"/* _re *= __Wi7;	_rf *= __Wi7; */\
		"addpd		%%xmm14,%%xmm12		\n\t	subpd		%%xmm15,%%xmm13		\n\t"/* _rc += _re;		_rd -= _rf; */\
		/* Now do radix-2 butterfly: */\
		"movaps		%%xmm8 ,%%xmm10		\n\t	movaps		%%xmm9 ,%%xmm11		\n\t"/* _ra  = _r8;		_rb  = _r9; */\
		"addpd		%%xmm12,%%xmm8 		\n\t	addpd		%%xmm13,%%xmm9 		\n\t"/* _r8 += _rc;		_r9 += _rd; */\
		"subpd		%%xmm12,%%xmm10		\n\t	subpd		%%xmm13,%%xmm11		\n\t"/* _ra -= _rc;		_rb -= _rd; */\
	/* Reload Block 0-3 outputs into r4-7,c-f, combine to get the 2 length-4 subtransform... */\
		"movq		%[i0],%%rax			\n\t"\
		"movq		%[i1],%%rbx			\n\t"\
		"movq		%[i2],%%rcx			\n\t"\
		"movq		%[i3],%%rdx			\n\t"\
		"movaps		    (%%rax),%%xmm4 	\n\t	movaps		0x10(%%rax),%%xmm5 	\n\t"/* _r4 = __tr0;	_r5 = __ti0; */\
		"movaps		    (%%rbx),%%xmm6 	\n\t	movaps		0x10(%%rbx),%%xmm7 	\n\t"/* _r6 = __tr1;	_r7 = __ti1; */\
		"movaps		    (%%rcx),%%xmm12	\n\t	movaps		0x10(%%rcx),%%xmm13	\n\t"/* _rc = __tr2;	_rd = __ti2; */\
		"movaps		    (%%rdx),%%xmm14	\n\t	movaps		0x10(%%rdx),%%xmm15	\n\t"/* _re = __tr3;	_rf = __ti3; */\
		"movq		%[o0],%%rax			\n\t"/* Assumes user stuck a (vec_dbl)2.0 into this output slot prior to macro call. */\
		"subpd		%%xmm12,%%xmm4 		\n\t	subpd		%%xmm13,%%xmm5 		\n\t"/* _r4 -= _rc;		_r5 -= _rd; */\
		"subpd		%%xmm15,%%xmm6 		\n\t	subpd		%%xmm14,%%xmm7 		\n\t"/* _r6 -= _rf;		_r7 -= _re; */\
		"subpd		%%xmm8 ,%%xmm0 		\n\t	subpd		%%xmm9 ,%%xmm1 		\n\t"/* _r0 -= _r8;		_r1 -= _r9; */\
		"subpd		%%xmm11,%%xmm2 		\n\t	subpd		%%xmm10,%%xmm3 		\n\t"/* _r2 -= _rb;		_r3 -= _ra; */\
		/* We hope the microcode execution engine sticks the datum at (%%rax) into a virtual register and inlines the MULs with the above SUBs: */\
		"mulpd		(%%rax),%%xmm12		\n\t	mulpd		(%%rax),%%xmm13		\n\t"/* _rc *= _two;	_rd *= _two; */\
		"mulpd		(%%rax),%%xmm15		\n\t	mulpd		(%%rax),%%xmm14		\n\t"/* _rf *= _two;	_re *= _two; */\
		"mulpd		(%%rax),%%xmm8 		\n\t	mulpd		(%%rax),%%xmm9 		\n\t"/* _r8 *= _two;	_r9 *= _two; */\
		"mulpd		(%%rax),%%xmm11		\n\t	mulpd		(%%rax),%%xmm10		\n\t"/* _rb *= _two;	_ra *= _two; */\
		"addpd		%%xmm4 ,%%xmm12		\n\t	addpd		%%xmm5 ,%%xmm13		\n\t"/* _rc += _r4;		_rd += _r5; */\
		"addpd		%%xmm6 ,%%xmm15		\n\t	addpd		%%xmm7 ,%%xmm14		\n\t"/* _rf += _r6;		_re += _r7; */\
		"addpd		%%xmm0 ,%%xmm8 		\n\t	addpd		%%xmm1 ,%%xmm9 		\n\t"/* _r8 += _r0;		_r9 += _r1; */\
		"addpd		%%xmm2 ,%%xmm11		\n\t	addpd		%%xmm3 ,%%xmm10		\n\t"/* _rb += _r2;		_ra += _r3; */\
		/* In terms of our original scalar-code prototyping macro, the data are: __tr0 = _r[c,f,4,6,8,b,0,2], __ti0 = _r[d,7,5,e,9,3,1,a]; */\
	/* Now combine the two half-transforms: */\
		/* Need r2/3 +- a/b combos for the *ISRT2 preceding the output 4-7 radix-2 butterflies, so start them first: */\
		"subpd		%%xmm3 ,%%xmm11		\n\t	subpd		%%xmm10,%%xmm2 		\n\t"/* _rb -= _r3;		_r2 -= _ra; */\
		"subpd		%%xmm8 ,%%xmm12		\n\t	subpd		%%xmm9 ,%%xmm13		\n\t"/* _rc -= _r8;		_rd -= _r9; */\
		"subpd		%%xmm1 ,%%xmm4 		\n\t	subpd		%%xmm0 ,%%xmm5 		\n\t"/* _r4 -= _r1;		_r5 -= _r0; */\
		"mulpd		(%%rax),%%xmm3 		\n\t	mulpd		(%%rax),%%xmm10		\n\t"/* _r3 *= _two;	_ra *= _two; */\
		"mulpd		(%%rax),%%xmm8 		\n\t	mulpd		(%%rax),%%xmm9 		\n\t"/* _r8 *= _two;	_r9 *= _two; */\
		"mulpd		(%%rax),%%xmm1 		\n\t	mulpd		(%%rax),%%xmm0 		\n\t"/* _r1 *= _two;	_r0 *= _two; */\
		"addpd		%%xmm11,%%xmm3 		\n\t	addpd		%%xmm2 ,%%xmm10		\n\t"/* _r3 += _rb;		_ra += _r2; */\
		"addpd		%%xmm12,%%xmm8 		\n\t	addpd		%%xmm13,%%xmm9 		\n\t"/* _r8 += _rc;		_r9 += _rd; */\
		"addpd		%%xmm4 ,%%xmm1 		\n\t	addpd		%%xmm5 ,%%xmm0 		\n\t"/* _r1 += _r4;		_r0 += _r5; */\
		/*movq		%[o0],%%rax		[o0] already in rax */	\
		"movq		%[o1],%%rbx			\n\t"\
		"movq		%[o2],%%rcx			\n\t"\
		"movq		%[o3],%%rdx			\n\t"\
		"movaps		%%xmm12,    (%%rbx)	\n\t	movaps		%%xmm13,0x10(%%rbx)	\n\t"/* __Br1 = _rc;	__Bi1 = _rd; */\
		/* Use that _rc,d free to stick 2.0 into _rc and that [c4] in rdi to load ISRT2 from c4-1 into _rd: */\
		"movaps		    (%%rax),%%xmm12	\n\t	movaps		-0x10(%%rdi),%%xmm13\n\t"/* _rc = 2.0;		_rd = ISRT2; */\
		"movaps		%%xmm4 ,    (%%rdx)	\n\t	movaps		%%xmm0 ,0x10(%%rdx)	\n\t"/* __Br3 = _r4;	__Bi3 = _r0; */\
		"movaps		%%xmm8 ,    (%%rax)	\n\t	movaps		%%xmm9 ,0x10(%%rax)	\n\t"/* __Br0 = _r8;	__Bi0 = _r9; */\
		"movaps		%%xmm1 ,    (%%rcx)	\n\t	movaps		%%xmm5 ,0x10(%%rcx)	\n\t"/* __Br2 = _r1;	__Bi2 = _r5; */\
		"mulpd		%%xmm13,%%xmm3 		\n\t	mulpd		%%xmm13,%%xmm11		\n\t"/* _r3 *= ISRT2;	_rb *= ISRT2; */\
		"mulpd		%%xmm13,%%xmm2 		\n\t	mulpd		%%xmm13,%%xmm10		\n\t"/* _r2 *= ISRT2;	_ra *= ISRT2; */\
		"subpd		%%xmm3 ,%%xmm15		\n\t	subpd		%%xmm11,%%xmm7 		\n\t"/* _rf -= _r3;		_r7 -= _rb; */\
		"subpd		%%xmm2 ,%%xmm6 		\n\t	subpd		%%xmm10,%%xmm14		\n\t"/* _r6 -= _r2;		_re -= _ra; */\
		"mulpd		%%xmm12,%%xmm3 		\n\t	mulpd		%%xmm12,%%xmm11		\n\t"/* _r3 *= _two;	_rb *= _two; */\
		"mulpd		%%xmm12,%%xmm2 		\n\t	mulpd		%%xmm12,%%xmm10		\n\t"/* _r2 *= _two;	_ra *= _two; */\
		"addpd		%%xmm15,%%xmm3 		\n\t	addpd		%%xmm7 ,%%xmm11		\n\t"/* _r3 += _rf;		_rb += _r7; */\
		"addpd		%%xmm6 ,%%xmm2 		\n\t	addpd		%%xmm14,%%xmm10		\n\t"/* _r2 += _r6;		_ra += _re; */\
		"movq		%[o4],%%rax			\n\t"\
		"movq		%[o5],%%rbx			\n\t"\
		"movq		%[o6],%%rcx			\n\t"\
		"movq		%[o7],%%rdx			\n\t"\
		"movaps		%%xmm3 ,    (%%rax)	\n\t	movaps		%%xmm7 ,0x10(%%rax)	\n\t"/* __Br4 = _r3;	__Bi4 = _r7; */\
		"movaps		%%xmm15,    (%%rbx)	\n\t	movaps		%%xmm11,0x10(%%rbx)	\n\t"/* __Br5 = _rf;	__Bi5 = _rb; */\
		"movaps		%%xmm6 ,    (%%rcx)	\n\t	movaps		%%xmm14,0x10(%%rcx)	\n\t"/* __Br6 = _r6;	__Bi6 = _re; */\
		"movaps		%%xmm2 ,    (%%rdx)	\n\t	movaps		%%xmm10,0x10(%%rdx)	\n\t"/* __Br7 = _r2;	__Bi7 = _ra; */\
		:					/* outputs: none */\
		: [i0] "m" (Xi0)	/* All inputs from memory addresses here */\
		 ,[i1] "m" (Xi1)\
		 ,[i2] "m" (Xi2)\
		 ,[i3] "m" (Xi3)\
		 ,[i4] "m" (Xi4)\
		 ,[i5] "m" (Xi5)\
		 ,[i6] "m" (Xi6)\
		 ,[i7] "m" (Xi7)\
		 ,[o0] "m" (Xo0)\
		 ,[o1] "m" (Xo1)\
		 ,[o2] "m" (Xo2)\
		 ,[o3] "m" (Xo3)\
		 ,[o4] "m" (Xo4)\
		 ,[o5] "m" (Xo5)\
		 ,[o6] "m" (Xo6)\
		 ,[o7] "m" (Xo7)\
		 ,[c1] "m" (Xc1),[s1] "m" (Xs1)\
		 ,[c2] "m" (Xc2),[s2] "m" (Xs2)\
		 ,[c3] "m" (Xc3),[s3] "m" (Xs3)\
		 ,[c4] "m" (Xc4),[s4] "m" (Xs4)\
		 ,[c5] "m" (Xc5),[s5] "m" (Xs5)\
		 ,[c6] "m" (Xc6),[s6] "m" (Xs6)\
		 ,[c7] "m" (Xc7),[s7] "m" (Xs7)\
		: "cc","memory","rax","rbx","rcx","rdx","rsi","rdi","r8","r9","r10","r11","r12","r13","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
	);\
	}

#endif	// AVX / SSE2 toggle

#endif	/* sse2_macro_gcc_h_included */

