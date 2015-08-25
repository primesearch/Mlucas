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
#ifndef radix52_ditN_cy_dif1_gcc_h_included
#define radix52_ditN_cy_dif1_gcc_h_included

	#define	SSE2_RADIX52_DIT_NOTWIDDLE(add0,p01,p02,p03,p04,r00)\
	{\
	__asm__ volatile (\
		"/*	Block 1: add0,1,3,2 = &a[j1+p00]+p0,1,2,3: SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_C(eax,ebx,edx,ecx, r00, 0x1a0) */\n\t"\
			"movq	%[__r00],%%rsi		\n\t"\
			"movq	%[__add0],%%rax		\n\t"\
			"movslq	%[__p01],%%rbx		\n\t"\
			"movslq	%[__p02],%%rcx		\n\t"\
			"movslq	%[__p03],%%rdx		\n\t"\
			"movslq	%[__p04],%%rdi							\n\t"\
			"shlq	$3,%%rbx								\n\t"\
			"shlq	$3,%%rcx								\n\t"\
			"shlq	$3,%%rdx								\n\t"\
			"shlq	$3,%%rdi /* p04, shifted-for-doubles */	\n\t"\
			"addq	%%rax,%%rbx								\n\t		movq	%%rdi,%%r8	/* copy p04 */\n\t"\
			"addq	%%rax,%%rcx								\n\t"\
			"addq	%%rax,%%rdx								\n\t		shlq	$2,%%r8		/* r8 = p16 */\n\t"\
			"/* ecx <-> edx */								\n\t"\
			"movaps	     (%%rax),%%xmm2						\n\t		addq	%%rdi,%%r8	/* r8  = p20 */\n\t"\
			"movaps	     (%%rdx),%%xmm6						\n\t"\
			"movaps	0x010(%%rax),%%xmm3						\n\t		shlq	$1,%%rdi	/* rdi = p08 */\n\t"\
			"movaps	0x010(%%rdx),%%xmm7						\n\t"\
			"movaps	     (%%rbx),%%xmm0						\n\t		addq	%%rdi,%%r8	/* r8  = p28. done with this one. */\n\t"\
			"movaps	     (%%rcx),%%xmm4						\n\t"\
			"movaps	0x010(%%rbx),%%xmm1						\n\t	leaq (%%rdi,%%rdi,2),%%rdi	/* rdi = p24, for alternating +p28/-p24 in right/left cols */\n\t"\
			"movaps	0x010(%%rcx),%%xmm5						\n\t"\
			"subpd	%%xmm0,%%xmm2							\n\t"\
			"subpd	%%xmm4,%%xmm6							\n\t		/*	Block 9: add2,3,0,1 = &a[j1+p28]+p0,1,2,3: */	\n\t"\
			"subpd	%%xmm1,%%xmm3							\n\t			leaq	0x100(%%rsi),%%r9	/* r08 */			\n\t"\
			"subpd	%%xmm5,%%xmm7							\n\t			addq	%%r8,%%rax	/* add44 = add00+p28 */		\n\t"\
			"addpd	%%xmm0,%%xmm0							\n\t			addq	%%r8,%%rbx								\n\t"\
			"addpd	%%xmm4,%%xmm4							\n\t			addq	%%r8,%%rcx								\n\t"\
			"addpd	%%xmm1,%%xmm1							\n\t			addq	%%r8,%%rdx								\n\t"\
			"addpd	%%xmm5,%%xmm5							\n\t			/* e[ab]x <-> e[cd]x */							\n\t"\
			"addpd	%%xmm2,%%xmm0							\n\t			movaps	     (%%rcx),%%xmm10					\n\t"\
			"addpd	%%xmm6,%%xmm4							\n\t			movaps	     (%%rax),%%xmm14					\n\t"\
			"addpd	%%xmm3,%%xmm1							\n\t			movaps	0x010(%%rcx),%%xmm11					\n\t"\
			"addpd	%%xmm7,%%xmm5							\n\t			movaps	0x010(%%rax),%%xmm15					\n\t"\
			"subpd	%%xmm4,%%xmm0							\n\t			movaps	     (%%rdx),%%xmm8 					\n\t"\
			"subpd	%%xmm7,%%xmm2							\n\t			movaps	     (%%rbx),%%xmm12					\n\t"\
			"subpd	%%xmm5,%%xmm1							\n\t			movaps	0x010(%%rdx),%%xmm9 					\n\t"\
			"subpd	%%xmm6,%%xmm3							\n\t			movaps	0x010(%%rbx),%%xmm13					\n\t"\
			"movaps	%%xmm0,0x340(%%rsi)						\n\t			subpd	%%xmm8 ,%%xmm10							\n\t"\
			"movaps	%%xmm2,0x4e0(%%rsi)						\n\t			subpd	%%xmm12,%%xmm14							\n\t"\
			"movaps	%%xmm1,0x350(%%rsi)						\n\t			subpd	%%xmm9 ,%%xmm11							\n\t"\
			"movaps	%%xmm3,0x1b0(%%rsi)						\n\t			subpd	%%xmm13,%%xmm15							\n\t"\
			"addpd	%%xmm4,%%xmm4							\n\t			addpd	%%xmm8 ,%%xmm8 							\n\t"\
			"addpd	%%xmm7,%%xmm7							\n\t			addpd	%%xmm12,%%xmm12							\n\t"\
			"addpd	%%xmm5,%%xmm5							\n\t			addpd	%%xmm9 ,%%xmm9 							\n\t"\
			"addpd	%%xmm6,%%xmm6							\n\t			addpd	%%xmm13,%%xmm13							\n\t"\
			"addpd	%%xmm0,%%xmm4							\n\t			addpd	%%xmm10,%%xmm8 							\n\t"\
			"addpd	%%xmm2,%%xmm7							\n\t			addpd	%%xmm14,%%xmm12							\n\t"\
			"addpd	%%xmm1,%%xmm5							\n\t			addpd	%%xmm11,%%xmm9 							\n\t"\
			"addpd	%%xmm3,%%xmm6							\n\t			addpd	%%xmm15,%%xmm13							\n\t"\
			"movaps	%%xmm4,     (%%rsi)						\n\t			subpd	%%xmm12,%%xmm8 							\n\t"\
			"movaps	%%xmm7,0x1a0(%%rsi)						\n\t			subpd	%%xmm15,%%xmm10							\n\t"\
			"movaps	%%xmm5,0x010(%%rsi)						\n\t			subpd	%%xmm13,%%xmm9 							\n\t"\
			"movaps	%%xmm6,0x4f0(%%rsi)						\n\t			subpd	%%xmm14,%%xmm11							\n\t"\
			"												\n\t			movaps	%%xmm8 ,0x340(%%r9)						\n\t"\
		"/*	Block 4: add3,2,1,0 = &a[j1+p04]+p0,1,2,3: */	\n\t			movaps	%%xmm10,0x4e0(%%r9)						\n\t"\
			"addq	$0x60,%%rsi	/* r03 */					\n\t			movaps	%%xmm9 ,0x350(%%r9)						\n\t"\
			"subq	%%rdi,%%rax	/* add04 = add24-p24 */		\n\t			movaps	%%xmm11,0x1b0(%%r9)						\n\t"\
			"subq	%%rdi,%%rbx								\n\t			addpd	%%xmm12,%%xmm12							\n\t"\
			"subq	%%rdi,%%rcx								\n\t			addpd	%%xmm15,%%xmm15							\n\t"\
			"subq	%%rdi,%%rdx								\n\t			addpd	%%xmm13,%%xmm13							\n\t"\
			"/* eax <-> edx, ebx <-> ecx */					\n\t			addpd	%%xmm14,%%xmm14							\n\t"\
			"movaps	     (%%rdx),%%xmm2						\n\t			addpd	%%xmm8 ,%%xmm12							\n\t"\
			"movaps	     (%%rbx),%%xmm6						\n\t			addpd	%%xmm10,%%xmm15							\n\t"\
			"movaps	0x010(%%rdx),%%xmm3						\n\t			addpd	%%xmm9 ,%%xmm13							\n\t"\
			"movaps	0x010(%%rbx),%%xmm7						\n\t			addpd	%%xmm11,%%xmm14							\n\t"\
			"movaps	     (%%rcx),%%xmm0						\n\t			movaps	%%xmm12,     (%%r9)						\n\t"\
			"movaps	     (%%rax),%%xmm4						\n\t			movaps	%%xmm15,0x1a0(%%r9)						\n\t"\
			"movaps	0x010(%%rcx),%%xmm1						\n\t			movaps	%%xmm13,0x010(%%r9)						\n\t"\
			"movaps	0x010(%%rax),%%xmm5						\n\t			movaps	%%xmm14,0x4f0(%%r9)						\n\t"\
			"subpd	%%xmm0,%%xmm2							\n\t		/*	Block 12: add0,1,3,2 = &a[j1+p32]+p0,1,2,3: */	\n\t"\
			"subpd	%%xmm4,%%xmm6							\n\t			addq	$0x60,%%r9	/* r0b */					\n\t"\
			"subpd	%%xmm1,%%xmm3							\n\t			addq	%%r8,%%rax	/* add32 = add04+p28 */		\n\t"\
			"subpd	%%xmm5,%%xmm7							\n\t			addq	%%r8,%%rbx								\n\t"\
			"addpd	%%xmm0,%%xmm0							\n\t			addq	%%r8,%%rcx								\n\t"\
			"addpd	%%xmm4,%%xmm4							\n\t			addq	%%r8,%%rdx								\n\t"\
			"addpd	%%xmm1,%%xmm1							\n\t			/* ecx <-> edx */								\n\t"\
			"addpd	%%xmm5,%%xmm5							\n\t			movaps	     (%%rax),%%xmm10					\n\t"\
			"addpd	%%xmm2,%%xmm0							\n\t			movaps	     (%%rdx),%%xmm14					\n\t"\
			"addpd	%%xmm6,%%xmm4							\n\t			movaps	0x010(%%rax),%%xmm11					\n\t"\
			"addpd	%%xmm3,%%xmm1							\n\t			movaps	0x010(%%rdx),%%xmm15					\n\t"\
			"addpd	%%xmm7,%%xmm5							\n\t			movaps	     (%%rbx),%%xmm8 					\n\t"\
			"subpd	%%xmm4,%%xmm0							\n\t			movaps	     (%%rcx),%%xmm12					\n\t"\
			"subpd	%%xmm7,%%xmm2							\n\t			movaps	0x010(%%rbx),%%xmm9 					\n\t"\
			"subpd	%%xmm5,%%xmm1							\n\t			movaps	0x010(%%rcx),%%xmm13					\n\t"\
			"subpd	%%xmm6,%%xmm3							\n\t			subpd	%%xmm8 ,%%xmm10							\n\t"\
			"movaps	%%xmm0,0x340(%%rsi)						\n\t			subpd	%%xmm12,%%xmm14							\n\t"\
			"movaps	%%xmm2,0x4e0(%%rsi)						\n\t			subpd	%%xmm9 ,%%xmm11							\n\t"\
			"movaps	%%xmm1,0x350(%%rsi)						\n\t			subpd	%%xmm13,%%xmm15							\n\t"\
			"movaps	%%xmm3,0x1b0(%%rsi)						\n\t			addpd	%%xmm8 ,%%xmm8 							\n\t"\
			"addpd	%%xmm4,%%xmm4							\n\t			addpd	%%xmm12,%%xmm12							\n\t"\
			"addpd	%%xmm7,%%xmm7							\n\t			addpd	%%xmm9 ,%%xmm9 							\n\t"\
			"addpd	%%xmm5,%%xmm5							\n\t			addpd	%%xmm13,%%xmm13							\n\t"\
			"addpd	%%xmm6,%%xmm6							\n\t			addpd	%%xmm10,%%xmm8 							\n\t"\
			"addpd	%%xmm0,%%xmm4							\n\t			addpd	%%xmm14,%%xmm12							\n\t"\
			"addpd	%%xmm2,%%xmm7							\n\t			addpd	%%xmm11,%%xmm9 							\n\t"\
			"addpd	%%xmm1,%%xmm5							\n\t			addpd	%%xmm15,%%xmm13							\n\t"\
			"addpd	%%xmm3,%%xmm6							\n\t			subpd	%%xmm12,%%xmm8 							\n\t"\
			"movaps	%%xmm4,     (%%rsi)						\n\t			subpd	%%xmm15,%%xmm10							\n\t"\
			"movaps	%%xmm7,0x1a0(%%rsi)						\n\t			subpd	%%xmm13,%%xmm9 							\n\t"\
			"movaps	%%xmm5,0x010(%%rsi)						\n\t			subpd	%%xmm14,%%xmm11							\n\t"\
			"movaps	%%xmm6,0x4f0(%%rsi)						\n\t			movaps	%%xmm8 ,0x340(%%r9)						\n\t"\
		"/*	Block 7: add1,0,2,3 = &a[j1+p08]+p0,1,2,3: */	\n\t			movaps	%%xmm10,0x4e0(%%r9)						\n\t"\
			"addq	$0x60,%%rsi	/* r06 */					\n\t			movaps	%%xmm9 ,0x350(%%r9)						\n\t"\
			"subq	%%rdi,%%rax	/* add08 = add32-p24 */		\n\t			movaps	%%xmm11,0x1b0(%%r9)						\n\t"\
			"subq	%%rdi,%%rbx								\n\t			addpd	%%xmm12,%%xmm12							\n\t"\
			"subq	%%rdi,%%rcx								\n\t			addpd	%%xmm15,%%xmm15							\n\t"\
			"subq	%%rdi,%%rdx								\n\t			addpd	%%xmm13,%%xmm13							\n\t"\
			"/* eax <-> ebx */								\n\t			addpd	%%xmm14,%%xmm14							\n\t"\
			"movaps	     (%%rbx),%%xmm2						\n\t			addpd	%%xmm8 ,%%xmm12							\n\t"\
			"movaps	     (%%rcx),%%xmm6						\n\t			addpd	%%xmm10,%%xmm15							\n\t"\
			"movaps	0x010(%%rbx),%%xmm3						\n\t			addpd	%%xmm9 ,%%xmm13							\n\t"\
			"movaps	0x010(%%rcx),%%xmm7						\n\t			addpd	%%xmm11,%%xmm14							\n\t"\
			"movaps	     (%%rax),%%xmm0						\n\t			movaps	%%xmm12,     (%%r9)						\n\t"\
			"movaps	     (%%rdx),%%xmm4						\n\t			movaps	%%xmm15,0x1a0(%%r9)						\n\t"\
			"movaps	0x010(%%rax),%%xmm1						\n\t			movaps	%%xmm13,0x010(%%r9)						\n\t"\
			"movaps	0x010(%%rdx),%%xmm5						\n\t			movaps	%%xmm14,0x4f0(%%r9)						\n\t"\
			"subpd	%%xmm0,%%xmm2							\n\t															\n\t"\
			"subpd	%%xmm4,%%xmm6							\n\t		/* Block 2: add3,2,1,0 = &a[j1+p36]+p0,1,2,3: */	\n\t"\
			"subpd	%%xmm1,%%xmm3							\n\t			subq	$0x140,%%r9	/* r01 */					\n\t"\
			"subpd	%%xmm5,%%xmm7							\n\t			addq	%%r8,%%rax	/* add36 = add08+p28 */		\n\t"\
			"addpd	%%xmm0,%%xmm0							\n\t			addq	%%r8,%%rbx								\n\t"\
			"addpd	%%xmm4,%%xmm4							\n\t			addq	%%r8,%%rcx								\n\t"\
			"addpd	%%xmm1,%%xmm1							\n\t			addq	%%r8,%%rdx								\n\t"\
			"addpd	%%xmm5,%%xmm5							\n\t			/* eax <-> edx, ebx <-> ecx */					\n\t"\
			"addpd	%%xmm2,%%xmm0							\n\t			movaps	     (%%rdx),%%xmm10					\n\t"\
			"addpd	%%xmm6,%%xmm4							\n\t			movaps	     (%%rbx),%%xmm14					\n\t"\
			"addpd	%%xmm3,%%xmm1							\n\t			movaps	0x010(%%rdx),%%xmm11					\n\t"\
			"addpd	%%xmm7,%%xmm5							\n\t			movaps	0x010(%%rbx),%%xmm15					\n\t"\
			"subpd	%%xmm4,%%xmm0							\n\t			movaps	     (%%rcx),%%xmm8 					\n\t"\
			"subpd	%%xmm7,%%xmm2							\n\t			movaps	     (%%rax),%%xmm12					\n\t"\
			"subpd	%%xmm5,%%xmm1							\n\t			movaps	0x010(%%rcx),%%xmm9 					\n\t"\
			"subpd	%%xmm6,%%xmm3							\n\t			movaps	0x010(%%rax),%%xmm13					\n\t"\
			"movaps	%%xmm0,0x340(%%rsi)						\n\t			subpd	%%xmm8 ,%%xmm10							\n\t"\
			"movaps	%%xmm2,0x4e0(%%rsi)						\n\t			subpd	%%xmm12,%%xmm14							\n\t"\
			"movaps	%%xmm1,0x350(%%rsi)						\n\t			subpd	%%xmm9 ,%%xmm11							\n\t"\
			"movaps	%%xmm3,0x1b0(%%rsi)						\n\t			subpd	%%xmm13,%%xmm15							\n\t"\
			"addpd	%%xmm4,%%xmm4							\n\t			addpd	%%xmm8 ,%%xmm8 							\n\t"\
			"addpd	%%xmm7,%%xmm7							\n\t			addpd	%%xmm12,%%xmm12							\n\t"\
			"addpd	%%xmm5,%%xmm5							\n\t			addpd	%%xmm9 ,%%xmm9 							\n\t"\
			"addpd	%%xmm6,%%xmm6							\n\t			addpd	%%xmm13,%%xmm13							\n\t"\
			"addpd	%%xmm0,%%xmm4							\n\t			addpd	%%xmm10,%%xmm8 							\n\t"\
			"addpd	%%xmm2,%%xmm7							\n\t			addpd	%%xmm14,%%xmm12							\n\t"\
			"addpd	%%xmm1,%%xmm5							\n\t			addpd	%%xmm11,%%xmm9 							\n\t"\
			"addpd	%%xmm3,%%xmm6							\n\t			addpd	%%xmm15,%%xmm13							\n\t"\
			"movaps	%%xmm4,     (%%rsi)						\n\t			subpd	%%xmm12,%%xmm8 							\n\t"\
			"movaps	%%xmm7,0x1a0(%%rsi)						\n\t			subpd	%%xmm15,%%xmm10							\n\t"\
			"movaps	%%xmm5,0x010(%%rsi)						\n\t			subpd	%%xmm13,%%xmm9 							\n\t"\
			"movaps	%%xmm6,0x4f0(%%rsi)						\n\t			subpd	%%xmm14,%%xmm11							\n\t"\
			"												\n\t			movaps	%%xmm8 ,0x340(%%r9)						\n\t"\
		"/*	Block 10: add2,3,0,1 = &a[j1+p12]+p0,1,2,3: */	\n\t			movaps	%%xmm10,0x4e0(%%r9)						\n\t"\
			"addq	$0x60,%%rsi	/* r09 */					\n\t			movaps	%%xmm9 ,0x350(%%r9)						\n\t"\
			"subq	%%rdi,%%rax	/* add12 = add36-p24 */		\n\t			movaps	%%xmm11,0x1b0(%%r9)						\n\t"\
			"subq	%%rdi,%%rbx								\n\t			addpd	%%xmm12,%%xmm12							\n\t"\
			"subq	%%rdi,%%rcx								\n\t			addpd	%%xmm15,%%xmm15							\n\t"\
			"subq	%%rdi,%%rdx								\n\t			addpd	%%xmm13,%%xmm13							\n\t"\
			"/* e[ab]x <-> e[cd]x */						\n\t			addpd	%%xmm14,%%xmm14							\n\t"\
			"movaps	     (%%rcx),%%xmm2						\n\t			addpd	%%xmm8 ,%%xmm12							\n\t"\
			"movaps	     (%%rax),%%xmm6						\n\t			addpd	%%xmm10,%%xmm15							\n\t"\
			"movaps	0x010(%%rcx),%%xmm3						\n\t			addpd	%%xmm9 ,%%xmm13							\n\t"\
			"movaps	0x010(%%rax),%%xmm7						\n\t			addpd	%%xmm11,%%xmm14							\n\t"\
			"movaps	     (%%rdx),%%xmm0						\n\t			movaps	%%xmm12,     (%%r9)						\n\t"\
			"movaps	     (%%rbx),%%xmm4						\n\t			movaps	%%xmm15,0x1a0(%%r9)						\n\t"\
			"movaps	0x010(%%rdx),%%xmm1						\n\t			movaps	%%xmm13,0x010(%%r9)						\n\t"\
			"movaps	0x010(%%rbx),%%xmm5						\n\t			movaps	%%xmm14,0x4f0(%%r9)						\n\t"\
			"subpd	%%xmm0,%%xmm2							\n\t		/*	Block 5: add1,0,2,3 = &a[j1+p40]+p0,1,2,3: */	\n\t"\
			"subpd	%%xmm4,%%xmm6							\n\t			addq	$0x60,%%r9	/* r04 */					\n\t"\
			"subpd	%%xmm1,%%xmm3							\n\t			addq	%%r8,%%rax	/* add40 = add12+p28 */		\n\t"\
			"subpd	%%xmm5,%%xmm7							\n\t			addq	%%r8,%%rbx								\n\t"\
			"addpd	%%xmm0,%%xmm0							\n\t			addq	%%r8,%%rcx								\n\t"\
			"addpd	%%xmm4,%%xmm4							\n\t			addq	%%r8,%%rdx								\n\t"\
			"addpd	%%xmm1,%%xmm1							\n\t			/* eax <-> ebx */								\n\t"\
			"addpd	%%xmm5,%%xmm5							\n\t			movaps	     (%%rbx),%%xmm10					\n\t"\
			"addpd	%%xmm2,%%xmm0							\n\t			movaps	     (%%rcx),%%xmm14					\n\t"\
			"addpd	%%xmm6,%%xmm4							\n\t			movaps	0x010(%%rbx),%%xmm11					\n\t"\
			"addpd	%%xmm3,%%xmm1							\n\t			movaps	0x010(%%rcx),%%xmm15					\n\t"\
			"addpd	%%xmm7,%%xmm5							\n\t			movaps	     (%%rax),%%xmm8 					\n\t"\
			"subpd	%%xmm4,%%xmm0							\n\t			movaps	     (%%rdx),%%xmm12					\n\t"\
			"subpd	%%xmm7,%%xmm2							\n\t			movaps	0x010(%%rax),%%xmm9 					\n\t"\
			"subpd	%%xmm5,%%xmm1							\n\t			movaps	0x010(%%rdx),%%xmm13					\n\t"\
			"subpd	%%xmm6,%%xmm3							\n\t			subpd	%%xmm8 ,%%xmm10							\n\t"\
			"movaps	%%xmm0,0x340(%%rsi)						\n\t			subpd	%%xmm12,%%xmm14							\n\t"\
			"movaps	%%xmm2,0x4e0(%%rsi)						\n\t			subpd	%%xmm9 ,%%xmm11							\n\t"\
			"movaps	%%xmm1,0x350(%%rsi)						\n\t			subpd	%%xmm13,%%xmm15							\n\t"\
			"movaps	%%xmm3,0x1b0(%%rsi)						\n\t			addpd	%%xmm8 ,%%xmm8 							\n\t"\
			"addpd	%%xmm4,%%xmm4							\n\t			addpd	%%xmm12,%%xmm12							\n\t"\
			"addpd	%%xmm7,%%xmm7							\n\t			addpd	%%xmm9 ,%%xmm9 							\n\t"\
			"addpd	%%xmm5,%%xmm5							\n\t			addpd	%%xmm13,%%xmm13							\n\t"\
			"addpd	%%xmm6,%%xmm6							\n\t			addpd	%%xmm10,%%xmm8 							\n\t"\
			"addpd	%%xmm0,%%xmm4							\n\t			addpd	%%xmm14,%%xmm12							\n\t"\
			"addpd	%%xmm2,%%xmm7							\n\t			addpd	%%xmm11,%%xmm9 							\n\t"\
			"addpd	%%xmm1,%%xmm5							\n\t			addpd	%%xmm15,%%xmm13							\n\t"\
			"addpd	%%xmm3,%%xmm6							\n\t			subpd	%%xmm12,%%xmm8 							\n\t"\
			"movaps	%%xmm4,     (%%rsi)						\n\t			subpd	%%xmm15,%%xmm10							\n\t"\
			"movaps	%%xmm7,0x1a0(%%rsi)						\n\t			subpd	%%xmm13,%%xmm9 							\n\t"\
			"movaps	%%xmm5,0x010(%%rsi)						\n\t			subpd	%%xmm14,%%xmm11							\n\t"\
			"movaps	%%xmm6,0x4f0(%%rsi)						\n\t			movaps	%%xmm8 ,0x340(%%r9)						\n\t"\
		"/*	Block 13: add0,1,3,2 = &a[j1+p16]+p0,1,2,3: */	\n\t			movaps	%%xmm10,0x4e0(%%r9)						\n\t"\
			"addq	$0x60,%%rsi	/* r0c */					\n\t			movaps	%%xmm9 ,0x350(%%r9)						\n\t"\
			"subq	%%rdi,%%rax	/* add16 = add40-p24 */		\n\t			movaps	%%xmm11,0x1b0(%%r9)						\n\t"\
			"subq	%%rdi,%%rbx								\n\t			addpd	%%xmm12,%%xmm12							\n\t"\
			"subq	%%rdi,%%rcx								\n\t			addpd	%%xmm15,%%xmm15							\n\t"\
			"subq	%%rdi,%%rdx								\n\t			addpd	%%xmm13,%%xmm13							\n\t"\
			"/* ecx <-> edx */								\n\t			addpd	%%xmm14,%%xmm14							\n\t"\
			"movaps	     (%%rax),%%xmm2						\n\t			addpd	%%xmm8 ,%%xmm12							\n\t"\
			"movaps	     (%%rdx),%%xmm6						\n\t			addpd	%%xmm10,%%xmm15							\n\t"\
			"movaps	0x010(%%rax),%%xmm3						\n\t			addpd	%%xmm9 ,%%xmm13							\n\t"\
			"movaps	0x010(%%rdx),%%xmm7						\n\t			addpd	%%xmm11,%%xmm14							\n\t"\
			"movaps	     (%%rbx),%%xmm0						\n\t			movaps	%%xmm12,     (%%r9)						\n\t"\
			"movaps	     (%%rcx),%%xmm4						\n\t			movaps	%%xmm15,0x1a0(%%r9)						\n\t"\
			"movaps	0x010(%%rbx),%%xmm1						\n\t			movaps	%%xmm13,0x010(%%r9)						\n\t"\
			"movaps	0x010(%%rcx),%%xmm5						\n\t			movaps	%%xmm14,0x4f0(%%r9)						\n\t"\
			"subpd	%%xmm0,%%xmm2							\n\t															\n\t"\
			"subpd	%%xmm4,%%xmm6							\n\t		/* Block 8: add2,3,0,1 = &a[j1+p44]+p0,1,2,3: */	\n\t"\
			"subpd	%%xmm1,%%xmm3							\n\t			addq	$0x60,%%r9	/* r07 */					\n\t"\
			"subpd	%%xmm5,%%xmm7							\n\t			addq	%%r8,%%rax	/* add44 = add16+p28 */		\n\t"\
			"addpd	%%xmm0,%%xmm0							\n\t			addq	%%r8,%%rbx								\n\t"\
			"addpd	%%xmm4,%%xmm4							\n\t			addq	%%r8,%%rcx								\n\t"\
			"addpd	%%xmm1,%%xmm1							\n\t			addq	%%r8,%%rdx								\n\t"\
			"addpd	%%xmm5,%%xmm5							\n\t			/* e[ab]x <-> e[cd]x */							\n\t"\
			"addpd	%%xmm2,%%xmm0							\n\t			movaps	     (%%rcx),%%xmm10					\n\t"\
			"addpd	%%xmm6,%%xmm4							\n\t			movaps	     (%%rax),%%xmm14					\n\t"\
			"addpd	%%xmm3,%%xmm1							\n\t			movaps	0x010(%%rcx),%%xmm11					\n\t"\
			"addpd	%%xmm7,%%xmm5							\n\t			movaps	0x010(%%rax),%%xmm15					\n\t"\
			"subpd	%%xmm4,%%xmm0							\n\t			movaps	     (%%rdx),%%xmm8 					\n\t"\
			"subpd	%%xmm7,%%xmm2							\n\t			movaps	     (%%rbx),%%xmm12					\n\t"\
			"subpd	%%xmm5,%%xmm1							\n\t			movaps	0x010(%%rdx),%%xmm9 					\n\t"\
			"subpd	%%xmm6,%%xmm3							\n\t			movaps	0x010(%%rbx),%%xmm13					\n\t"\
			"movaps	%%xmm0,0x340(%%rsi)						\n\t			subpd	%%xmm8 ,%%xmm10							\n\t"\
			"movaps	%%xmm2,0x4e0(%%rsi)						\n\t			subpd	%%xmm12,%%xmm14							\n\t"\
			"movaps	%%xmm1,0x350(%%rsi)						\n\t			subpd	%%xmm9 ,%%xmm11							\n\t"\
			"movaps	%%xmm3,0x1b0(%%rsi)						\n\t			subpd	%%xmm13,%%xmm15							\n\t"\
			"addpd	%%xmm4,%%xmm4							\n\t			addpd	%%xmm8 ,%%xmm8 							\n\t"\
			"addpd	%%xmm7,%%xmm7							\n\t			addpd	%%xmm12,%%xmm12							\n\t"\
			"addpd	%%xmm5,%%xmm5							\n\t			addpd	%%xmm9 ,%%xmm9 							\n\t"\
			"addpd	%%xmm6,%%xmm6							\n\t			addpd	%%xmm13,%%xmm13							\n\t"\
			"addpd	%%xmm0,%%xmm4							\n\t			addpd	%%xmm10,%%xmm8 							\n\t"\
			"addpd	%%xmm2,%%xmm7							\n\t			addpd	%%xmm14,%%xmm12							\n\t"\
			"addpd	%%xmm1,%%xmm5							\n\t			addpd	%%xmm11,%%xmm9 							\n\t"\
			"addpd	%%xmm3,%%xmm6							\n\t			addpd	%%xmm15,%%xmm13							\n\t"\
			"movaps	%%xmm4,     (%%rsi)						\n\t			subpd	%%xmm12,%%xmm8 							\n\t"\
			"movaps	%%xmm7,0x1a0(%%rsi)						\n\t			subpd	%%xmm15,%%xmm10							\n\t"\
			"movaps	%%xmm5,0x010(%%rsi)						\n\t			subpd	%%xmm13,%%xmm9 							\n\t"\
			"movaps	%%xmm6,0x4f0(%%rsi)						\n\t			subpd	%%xmm14,%%xmm11							\n\t"\
		"/*	Block 3: add3,2,1,0 = &a[j1+p20]+p0,1,2,3: */	\n\t			movaps	%%xmm8 ,0x340(%%r9)						\n\t"\
			"subq	$0x140,%%rsi	/* r02 */				\n\t			movaps	%%xmm10,0x4e0(%%r9)						\n\t"\
			"subq	%%rdi,%%rax	/* add20 = add44-p24 */		\n\t			movaps	%%xmm9 ,0x350(%%r9)						\n\t"\
			"subq	%%rdi,%%rbx								\n\t			movaps	%%xmm11,0x1b0(%%r9)						\n\t"\
			"subq	%%rdi,%%rcx								\n\t			addpd	%%xmm12,%%xmm12							\n\t"\
			"subq	%%rdi,%%rdx								\n\t			addpd	%%xmm15,%%xmm15							\n\t"\
			"/* eax <-> edx, ebx <-> ecx */					\n\t			addpd	%%xmm13,%%xmm13							\n\t"\
			"movaps	     (%%rdx),%%xmm2						\n\t			addpd	%%xmm14,%%xmm14							\n\t"\
			"movaps	     (%%rbx),%%xmm6						\n\t			addpd	%%xmm8 ,%%xmm12							\n\t"\
			"movaps	0x010(%%rdx),%%xmm3						\n\t			addpd	%%xmm10,%%xmm15							\n\t"\
			"movaps	0x010(%%rbx),%%xmm7						\n\t			addpd	%%xmm9 ,%%xmm13							\n\t"\
			"movaps	     (%%rcx),%%xmm0						\n\t			addpd	%%xmm11,%%xmm14							\n\t"\
			"movaps	     (%%rax),%%xmm4						\n\t			movaps	%%xmm12,     (%%r9)						\n\t"\
			"movaps	0x010(%%rcx),%%xmm1						\n\t			movaps	%%xmm15,0x1a0(%%r9)						\n\t"\
			"movaps	0x010(%%rax),%%xmm5						\n\t			movaps	%%xmm13,0x010(%%r9)						\n\t"\
			"subpd	%%xmm0,%%xmm2							\n\t			movaps	%%xmm14,0x4f0(%%r9)						\n\t"\
			"subpd	%%xmm4,%%xmm6							\n\t		/*	Block 11: add0,1,3,2 = &a[j1+p48]+p0,1,2,3: */	\n\t"\
			"subpd	%%xmm1,%%xmm3							\n\t			addq	$0x60,%%r9	/* r10 */					\n\t"\
			"subpd	%%xmm5,%%xmm7							\n\t			addq	%%r8,%%rax	/* add48 = add20+p28 */		\n\t"\
			"addpd	%%xmm0,%%xmm0							\n\t			addq	%%r8,%%rbx								\n\t"\
			"addpd	%%xmm4,%%xmm4							\n\t			addq	%%r8,%%rcx								\n\t"\
			"addpd	%%xmm1,%%xmm1							\n\t			addq	%%r8,%%rdx								\n\t"\
			"addpd	%%xmm5,%%xmm5							\n\t			/* ecx <-> edx */								\n\t"\
			"addpd	%%xmm2,%%xmm0							\n\t			movaps	     (%%rax),%%xmm10					\n\t"\
			"addpd	%%xmm6,%%xmm4							\n\t			movaps	     (%%rdx),%%xmm14					\n\t"\
			"addpd	%%xmm3,%%xmm1							\n\t			movaps	0x010(%%rax),%%xmm11					\n\t"\
			"addpd	%%xmm7,%%xmm5							\n\t			movaps	0x010(%%rdx),%%xmm15					\n\t"\
			"subpd	%%xmm4,%%xmm0							\n\t			movaps	     (%%rbx),%%xmm8 					\n\t"\
			"subpd	%%xmm7,%%xmm2							\n\t			movaps	     (%%rcx),%%xmm12					\n\t"\
			"subpd	%%xmm5,%%xmm1							\n\t			movaps	0x010(%%rbx),%%xmm9 					\n\t"\
			"subpd	%%xmm6,%%xmm3							\n\t			movaps	0x010(%%rcx),%%xmm13					\n\t"\
			"movaps	%%xmm0,0x340(%%rsi)						\n\t			subpd	%%xmm8 ,%%xmm10							\n\t"\
			"movaps	%%xmm2,0x4e0(%%rsi)						\n\t			subpd	%%xmm12,%%xmm14							\n\t"\
			"movaps	%%xmm1,0x350(%%rsi)						\n\t			subpd	%%xmm9 ,%%xmm11							\n\t"\
			"movaps	%%xmm3,0x1b0(%%rsi)						\n\t			subpd	%%xmm13,%%xmm15							\n\t"\
			"addpd	%%xmm4,%%xmm4							\n\t			addpd	%%xmm8 ,%%xmm8 							\n\t"\
			"addpd	%%xmm7,%%xmm7							\n\t			addpd	%%xmm12,%%xmm12							\n\t"\
			"addpd	%%xmm5,%%xmm5							\n\t			addpd	%%xmm9 ,%%xmm9 							\n\t"\
			"addpd	%%xmm6,%%xmm6							\n\t			addpd	%%xmm13,%%xmm13							\n\t"\
			"addpd	%%xmm0,%%xmm4							\n\t			addpd	%%xmm10,%%xmm8 							\n\t"\
			"addpd	%%xmm2,%%xmm7							\n\t			addpd	%%xmm14,%%xmm12							\n\t"\
			"addpd	%%xmm1,%%xmm5							\n\t			addpd	%%xmm11,%%xmm9 							\n\t"\
			"addpd	%%xmm3,%%xmm6							\n\t			addpd	%%xmm15,%%xmm13							\n\t"\
			"movaps	%%xmm4,     (%%rsi)						\n\t			subpd	%%xmm12,%%xmm8 							\n\t"\
			"movaps	%%xmm7,0x1a0(%%rsi)						\n\t			subpd	%%xmm15,%%xmm10							\n\t"\
			"movaps	%%xmm5,0x010(%%rsi)						\n\t			subpd	%%xmm13,%%xmm9 							\n\t"\
			"movaps	%%xmm6,0x4f0(%%rsi)						\n\t			subpd	%%xmm14,%%xmm11							\n\t"\
			"												\n\t			movaps	%%xmm8 ,0x340(%%r9)						\n\t"\
		"/*	Block 6: add1,0,2,3 = &a[j1+p24]+p0,1,2,3: */	\n\t			movaps	%%xmm10,0x4e0(%%r9)						\n\t"\
			"addq	$0x60,%%rsi	/* r05 */					\n\t			movaps	%%xmm9 ,0x350(%%r9)						\n\t"\
			"subq	%%rdi,%%rax	/* add24 = add48-p24 */		\n\t			movaps	%%xmm11,0x1b0(%%r9)						\n\t"\
			"subq	%%rdi,%%rbx								\n\t			addpd	%%xmm12,%%xmm12							\n\t"\
			"subq	%%rdi,%%rcx								\n\t			addpd	%%xmm15,%%xmm15							\n\t"\
			"subq	%%rdi,%%rdx								\n\t			addpd	%%xmm13,%%xmm13							\n\t"\
			"/* eax <-> ebx */								\n\t			addpd	%%xmm14,%%xmm14							\n\t"\
			"movaps	     (%%rbx),%%xmm2						\n\t			addpd	%%xmm8 ,%%xmm12							\n\t"\
			"movaps	     (%%rcx),%%xmm6						\n\t			addpd	%%xmm10,%%xmm15							\n\t"\
			"movaps	0x010(%%rbx),%%xmm3						\n\t			addpd	%%xmm9 ,%%xmm13							\n\t"\
			"movaps	0x010(%%rcx),%%xmm7						\n\t			addpd	%%xmm11,%%xmm14							\n\t"\
			"movaps	     (%%rax),%%xmm0						\n\t			movaps	%%xmm12,     (%%r9)						\n\t"\
			"movaps	     (%%rdx),%%xmm4						\n\t			movaps	%%xmm15,0x1a0(%%r9)						\n\t"\
			"movaps	0x010(%%rax),%%xmm1						\n\t			movaps	%%xmm13,0x010(%%r9)						\n\t"\
			"movaps	0x010(%%rdx),%%xmm5						\n\t			movaps	%%xmm14,0x4f0(%%r9)						\n\t"\
			"subpd	%%xmm0,%%xmm2							\n\t"\
			"subpd	%%xmm4,%%xmm6							\n\t"\
			"subpd	%%xmm1,%%xmm3							\n\t"\
			"subpd	%%xmm5,%%xmm7							\n\t"\
			"addpd	%%xmm0,%%xmm0							\n\t"\
			"addpd	%%xmm4,%%xmm4							\n\t"\
			"addpd	%%xmm1,%%xmm1							\n\t"\
			"addpd	%%xmm5,%%xmm5							\n\t"\
			"addpd	%%xmm2,%%xmm0							\n\t"\
			"addpd	%%xmm6,%%xmm4							\n\t"\
			"addpd	%%xmm3,%%xmm1							\n\t"\
			"addpd	%%xmm7,%%xmm5							\n\t"\
			"subpd	%%xmm4,%%xmm0							\n\t"\
			"subpd	%%xmm7,%%xmm2							\n\t"\
			"subpd	%%xmm5,%%xmm1							\n\t"\
			"subpd	%%xmm6,%%xmm3							\n\t"\
			"movaps	%%xmm0,0x340(%%rsi)						\n\t"\
			"movaps	%%xmm2,0x4e0(%%rsi)						\n\t"\
			"movaps	%%xmm1,0x350(%%rsi)						\n\t"\
			"movaps	%%xmm3,0x1b0(%%rsi)						\n\t"\
			"addpd	%%xmm4,%%xmm4							\n\t"\
			"addpd	%%xmm7,%%xmm7							\n\t"\
			"addpd	%%xmm5,%%xmm5							\n\t"\
			"addpd	%%xmm6,%%xmm6							\n\t"\
			"addpd	%%xmm0,%%xmm4							\n\t"\
			"addpd	%%xmm2,%%xmm7							\n\t"\
			"addpd	%%xmm1,%%xmm5							\n\t"\
			"addpd	%%xmm3,%%xmm6							\n\t"\
			"movaps	%%xmm4,     (%%rsi)						\n\t"\
			"movaps	%%xmm7,0x1a0(%%rsi)						\n\t"\
			"movaps	%%xmm5,0x010(%%rsi)						\n\t"\
			"movaps	%%xmm6,0x4f0(%%rsi)						\n\t"\
			:					/* outputs: none */\
			: [__add0] "m" (add0)	/* All inputs from memory addresses here */\
			 ,[__p01] "m" (p01)\
			 ,[__p02] "m" (p02)\
			 ,[__p03] "m" (p03)\
			 ,[__p04] "m" (p04)\
			 ,[__r00] "m" (r00)\
			: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r8","r9","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"		/* Clobbered registers */\
	);\
	}

	#define	SSE2_RADIX52_DIF_NOTWIDDLE(add0,p01,p02,p03,p04,r00)\
	{\
	__asm__ volatile (\
		"/*	Block 1: add0,1,3,2 = &a[j1+p00]+p0,1,2,3: SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r00, 0x1a0) */\n\t"\
			"movq	%[__r00],%%rsi						\n\t"\
			"movq	%[__add0],%%rax						\n\t"\
			"movslq	%[__p01],%%rbx						\n\t"\
			"movslq	%[__p02],%%rcx						\n\t"\
			"movslq	%[__p03],%%rdx						\n\t"\
			"movslq	%[__p04],%%rdi						\n\t"\
			"shlq	$3,%%rbx							\n\t"\
			"shlq	$3,%%rcx							\n\t"\
			"shlq	$3,%%rdx							\n\t"\
			"shlq	$3,%%rdi /* p04, shift-for-double */\n\t"\
			"addq	%%rax,%%rbx							\n\t"\
			"addq	%%rax,%%rcx							\n\t"\
			"addq	%%rax,%%rdx							\n\t"\
			"movaps	     (%%rsi),%%xmm4					\n\t"\
			"movaps	0x1a0(%%rsi),%%xmm6					\n\t"\
			"movaps	0x010(%%rsi),%%xmm5					\n\t"\
			"movaps	0x1b0(%%rsi),%%xmm7					\n\t"\
			"movaps	0x340(%%rsi),%%xmm0					\n\t		movq	%%rdi,%%r8	/* copy p04 */\n\t"\
			"movaps	0x4e0(%%rsi),%%xmm2					\n\t"\
			"movaps	0x350(%%rsi),%%xmm1					\n\t		shlq	$2,%%r8		/* r8 = p16 */\n\t"\
			"movaps	0x4f0(%%rsi),%%xmm3					\n\t"\
			"subpd	%%xmm0,%%xmm4						\n\t		addq	%%rdi,%%r8	/* r8  = p20 */\n\t"\
			"subpd	%%xmm2,%%xmm6						\n\t"\
			"subpd	%%xmm1,%%xmm5						\n\t		shlq	$1,%%rdi	/* rdi = p08 */\n\t"\
			"subpd	%%xmm3,%%xmm7						\n\t"\
			"addpd	%%xmm0,%%xmm0						\n\t		addq	%%rdi,%%r8	/* r8  = p28. done with this one. */\n\t"\
			"addpd	%%xmm2,%%xmm2						\n\t"\
			"addpd	%%xmm1,%%xmm1						\n\t	leaq (%%rdi,%%rdi,2),%%rdi	/* rdi = p24, for alternating +p28/-p24 in right/left cols */\n\t"\
			"addpd	%%xmm3,%%xmm3						\n\t"\
			"addpd	%%xmm4,%%xmm0						\n\t"\
			"addpd	%%xmm6,%%xmm2						\n\t	/*	Block 7: add2,3,0,1 = a[j1+p28]+p0-3 -> r06: */\n\t"\
			"addpd	%%xmm5,%%xmm1						\n\t		leaq	0xc0(%%rsi),%%r9	/* r06 */		\n\t"\
			"addpd	%%xmm7,%%xmm3						\n\t		movaps	     (%%r9),%%xmm12					\n\t"\
		"/* Unpermuted output addresses: bcbd|adac*/	\n\t		movaps	0x1a0(%%r9),%%xmm14					\n\t"\
		"/*	Swap output regs c / d: */					\n\t		movaps	0x010(%%r9),%%xmm13					\n\t"\
			"subpd	%%xmm2,%%xmm0						\n\t		movaps	0x1b0(%%r9),%%xmm15					\n\t"\
			"subpd	%%xmm7,%%xmm4						\n\t		movaps	0x340(%%r9),%%xmm8 					\n\t"\
			"subpd	%%xmm3,%%xmm1						\n\t		movaps	0x4e0(%%r9),%%xmm10					\n\t"\
			"subpd	%%xmm6,%%xmm5						\n\t		movaps	0x350(%%r9),%%xmm9 					\n\t"\
			"movaps	%%xmm0,    (%%rbx)					\n\t		movaps	0x4f0(%%r9),%%xmm11					\n\t"\
			"movaps	%%xmm4,    (%%rdx)					\n\t		subpd	%%xmm8 ,%%xmm12						\n\t"\
			"movaps	%%xmm1,0x10(%%rbx)					\n\t		subpd	%%xmm10,%%xmm14						\n\t"\
			"movaps	%%xmm5,0x10(%%rcx)					\n\t		subpd	%%xmm9 ,%%xmm13						\n\t"\
			"addpd	%%xmm2,%%xmm2						\n\t		subpd	%%xmm11,%%xmm15						\n\t"\
			"addpd	%%xmm7,%%xmm7						\n\t		addpd	%%xmm8 ,%%xmm8 						\n\t"\
			"addpd	%%xmm3,%%xmm3						\n\t		addpd	%%xmm10,%%xmm10						\n\t"\
			"addpd	%%xmm6,%%xmm6						\n\t		addpd	%%xmm9 ,%%xmm9 						\n\t"\
			"addpd	%%xmm0,%%xmm2						\n\t		addpd	%%xmm11,%%xmm11						\n\t"\
			"addpd	%%xmm4,%%xmm7						\n\t		addpd	%%xmm12,%%xmm8 						\n\t"\
			"addpd	%%xmm1,%%xmm3						\n\t		addpd	%%xmm14,%%xmm10						\n\t"\
			"addpd	%%xmm5,%%xmm6						\n\t		addpd	%%xmm13,%%xmm9 						\n\t"\
			"movaps	%%xmm2,    (%%rax)					\n\t		addpd	%%xmm15,%%xmm11						\n\t"\
			"movaps	%%xmm7,    (%%rcx)					\n\t	/* Swap outputs 2301[cdab]->0123[abcd], or abcd -> cdab: */\n\t"\
			"movaps	%%xmm3,0x10(%%rax)					\n\t		addq	%%r8,%%rax /* add28 */				\n\t"\
			"movaps	%%xmm6,0x10(%%rdx)					\n\t		addq	%%r8,%%rbx							\n\t"\
	"/*	Block 13: add3,2,1,0 = a[j1+p04]+p0-3 -> r0c: */\n\t		addq	%%r8,%%rcx							\n\t"\
	"/* out-regs abcd -> dcba (simple reverse-order). */\n\t		addq	%%r8,%%rdx							\n\t"\
			"addq	$0x180,%%rsi	/* r0c */			\n\t		subpd	%%xmm10,%%xmm8 						\n\t"\
			"movaps	     (%%rsi),%%xmm4					\n\t		subpd	%%xmm15,%%xmm12						\n\t"\
			"movaps	0x1a0(%%rsi),%%xmm6					\n\t		subpd	%%xmm11,%%xmm9 						\n\t"\
			"movaps	0x010(%%rsi),%%xmm5					\n\t		subpd	%%xmm14,%%xmm13						\n\t"\
			"movaps	0x1b0(%%rsi),%%xmm7					\n\t		movaps	%%xmm8 ,    (%%rdx)					\n\t"\
			"movaps	0x340(%%rsi),%%xmm0					\n\t		movaps	%%xmm12,    (%%rax)					\n\t"\
			"movaps	0x4e0(%%rsi),%%xmm2					\n\t		movaps	%%xmm9 ,0x10(%%rdx)					\n\t"\
			"movaps	0x350(%%rsi),%%xmm1					\n\t		movaps	%%xmm13,0x10(%%rbx)					\n\t"\
			"movaps	0x4f0(%%rsi),%%xmm3					\n\t		addpd	%%xmm10,%%xmm10						\n\t"\
			"subpd	%%xmm0,%%xmm4						\n\t		addpd	%%xmm15,%%xmm15						\n\t"\
			"subpd	%%xmm2,%%xmm6						\n\t		addpd	%%xmm11,%%xmm11						\n\t"\
			"subpd	%%xmm1,%%xmm5						\n\t		addpd	%%xmm14,%%xmm14						\n\t"\
			"subpd	%%xmm3,%%xmm7						\n\t		addpd	%%xmm8 ,%%xmm10						\n\t"\
			"addpd	%%xmm0,%%xmm0						\n\t		addpd	%%xmm12,%%xmm15						\n\t"\
			"addpd	%%xmm2,%%xmm2						\n\t		addpd	%%xmm9 ,%%xmm11						\n\t"\
			"addpd	%%xmm1,%%xmm1						\n\t		addpd	%%xmm13,%%xmm14						\n\t"\
			"addpd	%%xmm3,%%xmm3						\n\t		movaps	%%xmm10,    (%%rcx)					\n\t"\
			"addpd	%%xmm4,%%xmm0						\n\t		movaps	%%xmm15,    (%%rbx)					\n\t"\
			"addpd	%%xmm6,%%xmm2						\n\t		movaps	%%xmm11,0x10(%%rcx)					\n\t"\
			"addpd	%%xmm5,%%xmm1						\n\t		movaps	%%xmm14,0x10(%%rax)					\n\t"\
			"addpd	%%xmm7,%%xmm3						\n\t"\
			"/* Output addresses: swap abcd -> cdba: */	\n\t"\
			"subq	%%rdi,%%rax /* add04 */				\n\t	/*	Block 6: add0,1,3,2 = a[j1+p36]+p0-3 -> r05: */\n\t"\
			"subq	%%rdi,%%rbx							\n\t		subq	$0x20,%%r9	/* r05 */				\n\t"\
			"subq	%%rdi,%%rcx							\n\t		movaps	     (%%r9),%%xmm12					\n\t"\
			"subq	%%rdi,%%rdx							\n\t		movaps	0x1a0(%%r9),%%xmm14					\n\t"\
			"subpd	%%xmm2,%%xmm0						\n\t		movaps	0x010(%%r9),%%xmm13					\n\t"\
			"subpd	%%xmm7,%%xmm4						\n\t		movaps	0x1b0(%%r9),%%xmm15					\n\t"\
			"subpd	%%xmm3,%%xmm1						\n\t		movaps	0x340(%%r9),%%xmm8 					\n\t"\
			"subpd	%%xmm6,%%xmm5						\n\t		movaps	0x4e0(%%r9),%%xmm10					\n\t"\
			"movaps	%%xmm0,    (%%rcx)					\n\t		movaps	0x350(%%r9),%%xmm9 					\n\t"\
			"movaps	%%xmm4,    (%%rbx)					\n\t		movaps	0x4f0(%%r9),%%xmm11					\n\t"\
			"movaps	%%xmm1,0x10(%%rcx)					\n\t		subpd	%%xmm8 ,%%xmm12						\n\t"\
			"movaps	%%xmm5,0x10(%%rax)					\n\t		subpd	%%xmm10,%%xmm14						\n\t"\
			"addpd	%%xmm2,%%xmm2						\n\t		subpd	%%xmm9 ,%%xmm13						\n\t"\
			"addpd	%%xmm7,%%xmm7						\n\t		subpd	%%xmm11,%%xmm15						\n\t"\
			"addpd	%%xmm3,%%xmm3						\n\t		addpd	%%xmm8 ,%%xmm8 						\n\t"\
			"addpd	%%xmm6,%%xmm6						\n\t		addpd	%%xmm10,%%xmm10						\n\t"\
			"addpd	%%xmm0,%%xmm2						\n\t		addpd	%%xmm9 ,%%xmm9 						\n\t"\
			"addpd	%%xmm4,%%xmm7						\n\t		addpd	%%xmm11,%%xmm11						\n\t"\
			"addpd	%%xmm1,%%xmm3						\n\t		addpd	%%xmm12,%%xmm8 						\n\t"\
			"addpd	%%xmm5,%%xmm6						\n\t		addpd	%%xmm14,%%xmm10						\n\t"\
			"movaps	%%xmm2,    (%%rdx)					\n\t		addpd	%%xmm13,%%xmm9 						\n\t"\
			"movaps	%%xmm7,    (%%rax)					\n\t		addpd	%%xmm15,%%xmm11						\n\t"\
			"movaps	%%xmm3,0x10(%%rdx)					\n\t		/* Output addresses: swap c/d: */			\n\t"\
			"movaps	%%xmm6,0x10(%%rbx)					\n\t		addq	%%r8,%%rax /* add32 */				\n\t"\
		"/*	Block 12: add1,0,2,3 = a[j1+p08]+p0-3 -> r0b: */\n\t	addq	%%r8,%%rbx							\n\t"\
			"subq	$0x20,%%rsi	/* r0b */				\n\t		addq	%%r8,%%rcx							\n\t"\
			"movaps	     (%%rsi),%%xmm4					\n\t		addq	%%r8,%%rdx							\n\t"\
			"movaps	0x1a0(%%rsi),%%xmm6					\n\t		subpd	%%xmm10,%%xmm8 						\n\t"\
			"movaps	0x010(%%rsi),%%xmm5					\n\t		subpd	%%xmm15,%%xmm12						\n\t"\
			"movaps	0x1b0(%%rsi),%%xmm7					\n\t		subpd	%%xmm11,%%xmm9 						\n\t"\
			"movaps	0x340(%%rsi),%%xmm0					\n\t		subpd	%%xmm14,%%xmm13						\n\t"\
			"movaps	0x4e0(%%rsi),%%xmm2					\n\t		movaps	%%xmm8 ,    (%%rbx)					\n\t"\
			"movaps	0x350(%%rsi),%%xmm1					\n\t		movaps	%%xmm12,    (%%rdx)					\n\t"\
			"movaps	0x4f0(%%rsi),%%xmm3					\n\t		movaps	%%xmm9 ,0x10(%%rbx)					\n\t"\
			"subpd	%%xmm0,%%xmm4						\n\t		movaps	%%xmm13,0x10(%%rcx)					\n\t"\
			"subpd	%%xmm2,%%xmm6						\n\t		addpd	%%xmm10,%%xmm10						\n\t"\
			"subpd	%%xmm1,%%xmm5						\n\t		addpd	%%xmm15,%%xmm15						\n\t"\
			"subpd	%%xmm3,%%xmm7						\n\t		addpd	%%xmm11,%%xmm11						\n\t"\
			"addpd	%%xmm0,%%xmm0						\n\t		addpd	%%xmm14,%%xmm14						\n\t"\
			"addpd	%%xmm2,%%xmm2						\n\t		addpd	%%xmm8 ,%%xmm10						\n\t"\
			"addpd	%%xmm1,%%xmm1						\n\t		addpd	%%xmm12,%%xmm15						\n\t"\
			"addpd	%%xmm3,%%xmm3						\n\t		addpd	%%xmm9 ,%%xmm11						\n\t"\
			"addpd	%%xmm4,%%xmm0						\n\t		addpd	%%xmm13,%%xmm14						\n\t"\
			"addpd	%%xmm6,%%xmm2						\n\t		movaps	%%xmm10,    (%%rax)					\n\t"\
			"addpd	%%xmm5,%%xmm1						\n\t		movaps	%%xmm15,    (%%rcx)					\n\t"\
			"addpd	%%xmm7,%%xmm3						\n\t		movaps	%%xmm11,0x10(%%rax)					\n\t"\
			"/* Output addresses: swap a/b: */			\n\t		movaps	%%xmm14,0x10(%%rdx)					\n\t"\
			"subq	%%rdi,%%rax /* add08 */				\n\t	/*	Block 5: add3,2,1,0 = a[j1+p36]+p0-3 -> r04: */\n\t"\
			"subq	%%rdi,%%rbx							\n\t		subq	$0x20,%%r9	/* r04 */				\n\t"\
			"subq	%%rdi,%%rcx							\n\t		movaps	     (%%r9),%%xmm12					\n\t"\
			"subq	%%rdi,%%rdx							\n\t		movaps	0x1a0(%%r9),%%xmm14					\n\t"\
			"subpd	%%xmm2,%%xmm0						\n\t		movaps	0x010(%%r9),%%xmm13					\n\t"\
			"subpd	%%xmm7,%%xmm4						\n\t		movaps	0x1b0(%%r9),%%xmm15					\n\t"\
			"subpd	%%xmm3,%%xmm1						\n\t		movaps	0x340(%%r9),%%xmm8 					\n\t"\
			"subpd	%%xmm6,%%xmm5						\n\t		movaps	0x4e0(%%r9),%%xmm10					\n\t"\
			"movaps	%%xmm0,    (%%rax)					\n\t		movaps	0x350(%%r9),%%xmm9 					\n\t"\
			"movaps	%%xmm4,    (%%rcx)					\n\t		movaps	0x4f0(%%r9),%%xmm11					\n\t"\
			"movaps	%%xmm1,0x10(%%rax)					\n\t		subpd	%%xmm8 ,%%xmm12						\n\t"\
			"movaps	%%xmm5,0x10(%%rdx)					\n\t		subpd	%%xmm10,%%xmm14						\n\t"\
			"addpd	%%xmm2,%%xmm2						\n\t		subpd	%%xmm9 ,%%xmm13						\n\t"\
			"addpd	%%xmm7,%%xmm7						\n\t		subpd	%%xmm11,%%xmm15						\n\t"\
			"addpd	%%xmm3,%%xmm3						\n\t		addpd	%%xmm8 ,%%xmm8 						\n\t"\
			"addpd	%%xmm6,%%xmm6						\n\t		addpd	%%xmm10,%%xmm10						\n\t"\
			"addpd	%%xmm0,%%xmm2						\n\t		addpd	%%xmm9 ,%%xmm9 						\n\t"\
			"addpd	%%xmm4,%%xmm7						\n\t		addpd	%%xmm11,%%xmm11						\n\t"\
			"addpd	%%xmm1,%%xmm3						\n\t		addpd	%%xmm12,%%xmm8 						\n\t"\
			"addpd	%%xmm5,%%xmm6						\n\t		addpd	%%xmm14,%%xmm10						\n\t"\
			"movaps	%%xmm2,    (%%rbx)					\n\t		addpd	%%xmm13,%%xmm9 						\n\t"\
			"movaps	%%xmm7,    (%%rdx)					\n\t		addpd	%%xmm15,%%xmm11						\n\t"\
			"movaps	%%xmm3,0x10(%%rbx)					\n\t		/* Swap output addresses abcd -> dcba: */	\n\t"\
			"movaps	%%xmm6,0x10(%%rcx)					\n\t		addq	%%r8,%%rax /* add36 */				\n\t"\
		"/*	Block 11: add2,3,0,1 = a[j1+p12]+p0-3 -> r0a: */\n\t	addq	%%r8,%%rbx							\n\t"\
			"subq	$0x20,%%rsi	/* r0b */				\n\t		addq	%%r8,%%rcx							\n\t"\
			"movaps	     (%%rsi),%%xmm4					\n\t		addq	%%r8,%%rdx							\n\t"\
			"movaps	0x1a0(%%rsi),%%xmm6					\n\t		subpd	%%xmm10,%%xmm8 						\n\t"\
			"movaps	0x010(%%rsi),%%xmm5					\n\t		subpd	%%xmm15,%%xmm12						\n\t"\
			"movaps	0x1b0(%%rsi),%%xmm7					\n\t		subpd	%%xmm11,%%xmm9 						\n\t"\
			"movaps	0x340(%%rsi),%%xmm0					\n\t		subpd	%%xmm14,%%xmm13						\n\t"\
			"movaps	0x4e0(%%rsi),%%xmm2					\n\t		movaps	%%xmm8 ,    (%%rcx)					\n\t"\
			"movaps	0x350(%%rsi),%%xmm1					\n\t		movaps	%%xmm12,    (%%rbx)					\n\t"\
			"movaps	0x4f0(%%rsi),%%xmm3					\n\t		movaps	%%xmm9 ,0x10(%%rcx)					\n\t"\
			"subpd	%%xmm0,%%xmm4						\n\t		movaps	%%xmm13,0x10(%%rax)					\n\t"\
			"subpd	%%xmm2,%%xmm6						\n\t		addpd	%%xmm10,%%xmm10						\n\t"\
			"subpd	%%xmm1,%%xmm5						\n\t		addpd	%%xmm15,%%xmm15						\n\t"\
			"subpd	%%xmm3,%%xmm7						\n\t		addpd	%%xmm11,%%xmm11						\n\t"\
			"addpd	%%xmm0,%%xmm0						\n\t		addpd	%%xmm14,%%xmm14						\n\t"\
			"addpd	%%xmm2,%%xmm2						\n\t		addpd	%%xmm8 ,%%xmm10						\n\t"\
			"addpd	%%xmm1,%%xmm1						\n\t		addpd	%%xmm12,%%xmm15						\n\t"\
			"addpd	%%xmm3,%%xmm3						\n\t		addpd	%%xmm9 ,%%xmm11						\n\t"\
			"addpd	%%xmm4,%%xmm0						\n\t		addpd	%%xmm13,%%xmm14						\n\t"\
			"addpd	%%xmm6,%%xmm2						\n\t		movaps	%%xmm10,    (%%rdx)					\n\t"\
			"addpd	%%xmm5,%%xmm1						\n\t		movaps	%%xmm15,    (%%rax)					\n\t"\
			"addpd	%%xmm7,%%xmm3						\n\t		movaps	%%xmm11,0x10(%%rdx)					\n\t"\
			"/* Output addresses: swap abcd -> cdab: */	\n\t		movaps	%%xmm14,0x10(%%rbx)					\n\t"\
			"subq	%%rdi,%%rax /* add12 */				\n\t	/*	Block 4: add1,0,2,3 = a[j1+p40]+p0-3 -> r03: */\n\t"\
			"subq	%%rdi,%%rbx							\n\t		subq	$0x20,%%r9	/* r03 */				\n\t"\
			"subq	%%rdi,%%rcx							\n\t		movaps	     (%%r9),%%xmm12					\n\t"\
			"subq	%%rdi,%%rdx							\n\t		movaps	0x1a0(%%r9),%%xmm14					\n\t"\
			"subpd	%%xmm2,%%xmm0						\n\t		movaps	0x010(%%r9),%%xmm13					\n\t"\
			"subpd	%%xmm7,%%xmm4						\n\t		movaps	0x1b0(%%r9),%%xmm15					\n\t"\
			"subpd	%%xmm3,%%xmm1						\n\t		movaps	0x340(%%r9),%%xmm8 					\n\t"\
			"subpd	%%xmm6,%%xmm5						\n\t		movaps	0x4e0(%%r9),%%xmm10					\n\t"\
			"movaps	%%xmm0,    (%%rdx)					\n\t		movaps	0x350(%%r9),%%xmm9 					\n\t"\
			"movaps	%%xmm4,    (%%rax)					\n\t		movaps	0x4f0(%%r9),%%xmm11					\n\t"\
			"movaps	%%xmm1,0x10(%%rdx)					\n\t		subpd	%%xmm8 ,%%xmm12						\n\t"\
			"movaps	%%xmm5,0x10(%%rbx)					\n\t		subpd	%%xmm10,%%xmm14						\n\t"\
			"addpd	%%xmm2,%%xmm2						\n\t		subpd	%%xmm9 ,%%xmm13						\n\t"\
			"addpd	%%xmm7,%%xmm7						\n\t		subpd	%%xmm11,%%xmm15						\n\t"\
			"addpd	%%xmm3,%%xmm3						\n\t		addpd	%%xmm8 ,%%xmm8 						\n\t"\
			"addpd	%%xmm6,%%xmm6						\n\t		addpd	%%xmm10,%%xmm10						\n\t"\
			"addpd	%%xmm0,%%xmm2						\n\t		addpd	%%xmm9 ,%%xmm9 						\n\t"\
			"addpd	%%xmm4,%%xmm7						\n\t		addpd	%%xmm11,%%xmm11						\n\t"\
			"addpd	%%xmm1,%%xmm3						\n\t		addpd	%%xmm12,%%xmm8 						\n\t"\
			"addpd	%%xmm5,%%xmm6						\n\t		addpd	%%xmm14,%%xmm10						\n\t"\
			"movaps	%%xmm2,    (%%rcx)					\n\t		addpd	%%xmm13,%%xmm9 						\n\t"\
			"movaps	%%xmm7,    (%%rbx)					\n\t		addpd	%%xmm15,%%xmm11						\n\t"\
			"movaps	%%xmm3,0x10(%%rcx)					\n\t		/* Output addresses: swap a/b: */			\n\t"\
			"movaps	%%xmm6,0x10(%%rax)					\n\t		addq	%%r8,%%rax /* add40 */				\n\t"\
		"/*	Block 10: add0,1,3,2 = a[j1+p16]+p0-3 -> r09: */\n\t	addq	%%r8,%%rbx							\n\t"\
			"subq	$0x20,%%rsi	/* r09 */				\n\t		addq	%%r8,%%rcx							\n\t"\
			"movaps	     (%%rsi),%%xmm4					\n\t		addq	%%r8,%%rdx							\n\t"\
			"movaps	0x1a0(%%rsi),%%xmm6					\n\t		subpd	%%xmm10,%%xmm8 						\n\t"\
			"movaps	0x010(%%rsi),%%xmm5					\n\t		subpd	%%xmm15,%%xmm12						\n\t"\
			"movaps	0x1b0(%%rsi),%%xmm7					\n\t		subpd	%%xmm11,%%xmm9 						\n\t"\
			"movaps	0x340(%%rsi),%%xmm0					\n\t		subpd	%%xmm14,%%xmm13						\n\t"\
			"movaps	0x4e0(%%rsi),%%xmm2					\n\t		movaps	%%xmm8 ,    (%%rax)					\n\t"\
			"movaps	0x350(%%rsi),%%xmm1					\n\t		movaps	%%xmm12,    (%%rcx)					\n\t"\
			"movaps	0x4f0(%%rsi),%%xmm3					\n\t		movaps	%%xmm9 ,0x10(%%rax)					\n\t"\
			"subpd	%%xmm0,%%xmm4						\n\t		movaps	%%xmm13,0x10(%%rdx)					\n\t"\
			"subpd	%%xmm2,%%xmm6						\n\t		addpd	%%xmm10,%%xmm10						\n\t"\
			"subpd	%%xmm1,%%xmm5						\n\t		addpd	%%xmm15,%%xmm15						\n\t"\
			"subpd	%%xmm3,%%xmm7						\n\t		addpd	%%xmm11,%%xmm11						\n\t"\
			"addpd	%%xmm0,%%xmm0						\n\t		addpd	%%xmm14,%%xmm14						\n\t"\
			"addpd	%%xmm2,%%xmm2						\n\t		addpd	%%xmm8 ,%%xmm10						\n\t"\
			"addpd	%%xmm1,%%xmm1						\n\t		addpd	%%xmm12,%%xmm15						\n\t"\
			"addpd	%%xmm3,%%xmm3						\n\t		addpd	%%xmm9 ,%%xmm11						\n\t"\
			"addpd	%%xmm4,%%xmm0						\n\t		addpd	%%xmm13,%%xmm14						\n\t"\
			"addpd	%%xmm6,%%xmm2						\n\t		movaps	%%xmm10,    (%%rbx)					\n\t"\
			"addpd	%%xmm5,%%xmm1						\n\t		movaps	%%xmm15,    (%%rdx)					\n\t"\
			"addpd	%%xmm7,%%xmm3						\n\t		movaps	%%xmm11,0x10(%%rbx)					\n\t"\
			"/* Output addresses: swap c/d: */			\n\t		movaps	%%xmm14,0x10(%%rcx)					\n\t"\
			"subq	%%rdi,%%rax /* add16 */				\n\t	/*	Block 3: add2,3,0,1 = a[j1+p44]+p0-3 -> r05: */\n\t"\
			"subq	%%rdi,%%rbx							\n\t		subq	$0x20,%%r9	/* r02 */				\n\t"\
			"subq	%%rdi,%%rcx							\n\t		movaps	     (%%r9),%%xmm12					\n\t"\
			"subq	%%rdi,%%rdx							\n\t		movaps	0x1a0(%%r9),%%xmm14					\n\t"\
			"subpd	%%xmm2,%%xmm0						\n\t		movaps	0x010(%%r9),%%xmm13					\n\t"\
			"subpd	%%xmm7,%%xmm4						\n\t		movaps	0x1b0(%%r9),%%xmm15					\n\t"\
			"subpd	%%xmm3,%%xmm1						\n\t		movaps	0x340(%%r9),%%xmm8 					\n\t"\
			"subpd	%%xmm6,%%xmm5						\n\t		movaps	0x4e0(%%r9),%%xmm10					\n\t"\
			"movaps	%%xmm0,    (%%rbx)					\n\t		movaps	0x350(%%r9),%%xmm9 					\n\t"\
			"movaps	%%xmm4,    (%%rdx)					\n\t		movaps	0x4f0(%%r9),%%xmm11					\n\t"\
			"movaps	%%xmm1,0x10(%%rbx)					\n\t		subpd	%%xmm8 ,%%xmm12						\n\t"\
			"movaps	%%xmm5,0x10(%%rcx)					\n\t		subpd	%%xmm10,%%xmm14						\n\t"\
			"addpd	%%xmm2,%%xmm2						\n\t		subpd	%%xmm9 ,%%xmm13						\n\t"\
			"addpd	%%xmm7,%%xmm7						\n\t		subpd	%%xmm11,%%xmm15						\n\t"\
			"addpd	%%xmm3,%%xmm3						\n\t		addpd	%%xmm8 ,%%xmm8 						\n\t"\
			"addpd	%%xmm6,%%xmm6						\n\t		addpd	%%xmm10,%%xmm10						\n\t"\
			"addpd	%%xmm0,%%xmm2						\n\t		addpd	%%xmm9 ,%%xmm9 						\n\t"\
			"addpd	%%xmm4,%%xmm7						\n\t		addpd	%%xmm11,%%xmm11						\n\t"\
			"addpd	%%xmm1,%%xmm3						\n\t		addpd	%%xmm12,%%xmm8 						\n\t"\
			"addpd	%%xmm5,%%xmm6						\n\t		addpd	%%xmm14,%%xmm10						\n\t"\
			"movaps	%%xmm2,    (%%rax)					\n\t		addpd	%%xmm13,%%xmm9 						\n\t"\
			"movaps	%%xmm7,    (%%rcx)					\n\t		addpd	%%xmm15,%%xmm11						\n\t"\
			"movaps	%%xmm3,0x10(%%rax)					\n\t		/* Output addresses: swap abcd -> cdab: */	\n\t"\
			"movaps	%%xmm6,0x10(%%rdx)					\n\t		addq	%%r8,%%rax /* add44 */				\n\t"\
		"/*	Block 9: add3,2,1,0 = a[j1+p20]+p0-3 -> r08: */\n\t		addq	%%r8,%%rbx							\n\t"\
			"subq	$0x20,%%rsi	/* r08 */				\n\t		addq	%%r8,%%rcx							\n\t"\
			"movaps	     (%%rsi),%%xmm4					\n\t		addq	%%r8,%%rdx							\n\t"\
			"movaps	0x1a0(%%rsi),%%xmm6					\n\t		subpd	%%xmm10,%%xmm8 						\n\t"\
			"movaps	0x010(%%rsi),%%xmm5					\n\t		subpd	%%xmm15,%%xmm12						\n\t"\
			"movaps	0x1b0(%%rsi),%%xmm7					\n\t		subpd	%%xmm11,%%xmm9 						\n\t"\
			"movaps	0x340(%%rsi),%%xmm0					\n\t		subpd	%%xmm14,%%xmm13						\n\t"\
			"movaps	0x4e0(%%rsi),%%xmm2					\n\t		movaps	%%xmm8 ,    (%%rdx)					\n\t"\
			"movaps	0x350(%%rsi),%%xmm1					\n\t		movaps	%%xmm12,    (%%rax)					\n\t"\
			"movaps	0x4f0(%%rsi),%%xmm3					\n\t		movaps	%%xmm9 ,0x10(%%rdx)					\n\t"\
			"subpd	%%xmm0,%%xmm4						\n\t		movaps	%%xmm13,0x10(%%rbx)					\n\t"\
			"subpd	%%xmm2,%%xmm6						\n\t		addpd	%%xmm10,%%xmm10						\n\t"\
			"subpd	%%xmm1,%%xmm5						\n\t		addpd	%%xmm15,%%xmm15						\n\t"\
			"subpd	%%xmm3,%%xmm7						\n\t		addpd	%%xmm11,%%xmm11						\n\t"\
			"addpd	%%xmm0,%%xmm0						\n\t		addpd	%%xmm14,%%xmm14						\n\t"\
			"addpd	%%xmm2,%%xmm2						\n\t		addpd	%%xmm8 ,%%xmm10						\n\t"\
			"addpd	%%xmm1,%%xmm1						\n\t		addpd	%%xmm12,%%xmm15						\n\t"\
			"addpd	%%xmm3,%%xmm3						\n\t		addpd	%%xmm9 ,%%xmm11						\n\t"\
			"addpd	%%xmm4,%%xmm0						\n\t		addpd	%%xmm13,%%xmm14						\n\t"\
			"addpd	%%xmm6,%%xmm2						\n\t		movaps	%%xmm10,    (%%rcx)					\n\t"\
			"addpd	%%xmm5,%%xmm1						\n\t		movaps	%%xmm15,    (%%rbx)					\n\t"\
			"addpd	%%xmm7,%%xmm3						\n\t		movaps	%%xmm11,0x10(%%rcx)					\n\t"\
			"/* Output addresses: swap abcd -> dcba: */	\n\t		movaps	%%xmm14,0x10(%%rax)					\n\t"\
			"subq	%%rdi,%%rax /* add20 */				\n\t	/*	Block 2: add0,1,3,2 = a[j1+p48]+p0-3 -> r04: */\n\t"\
			"subq	%%rdi,%%rbx							\n\t		subq	$0x20,%%r9	/* r01 */				\n\t"\
			"subq	%%rdi,%%rcx							\n\t		movaps	     (%%r9),%%xmm12					\n\t"\
			"subq	%%rdi,%%rdx							\n\t		movaps	0x1a0(%%r9),%%xmm14					\n\t"\
			"subpd	%%xmm2,%%xmm0						\n\t		movaps	0x010(%%r9),%%xmm13					\n\t"\
			"subpd	%%xmm7,%%xmm4						\n\t		movaps	0x1b0(%%r9),%%xmm15					\n\t"\
			"subpd	%%xmm3,%%xmm1						\n\t		movaps	0x340(%%r9),%%xmm8 					\n\t"\
			"subpd	%%xmm6,%%xmm5						\n\t		movaps	0x4e0(%%r9),%%xmm10					\n\t"\
			"movaps	%%xmm0,    (%%rcx)					\n\t		movaps	0x350(%%r9),%%xmm9 					\n\t"\
			"movaps	%%xmm4,    (%%rbx)					\n\t		movaps	0x4f0(%%r9),%%xmm11					\n\t"\
			"movaps	%%xmm1,0x10(%%rcx)					\n\t		subpd	%%xmm8 ,%%xmm12						\n\t"\
			"movaps	%%xmm5,0x10(%%rax)					\n\t		subpd	%%xmm10,%%xmm14						\n\t"\
			"addpd	%%xmm2,%%xmm2						\n\t		subpd	%%xmm9 ,%%xmm13						\n\t"\
			"addpd	%%xmm7,%%xmm7						\n\t		subpd	%%xmm11,%%xmm15						\n\t"\
			"addpd	%%xmm3,%%xmm3						\n\t		addpd	%%xmm8 ,%%xmm8 						\n\t"\
			"addpd	%%xmm6,%%xmm6						\n\t		addpd	%%xmm10,%%xmm10						\n\t"\
			"addpd	%%xmm0,%%xmm2						\n\t		addpd	%%xmm9 ,%%xmm9 						\n\t"\
			"addpd	%%xmm4,%%xmm7						\n\t		addpd	%%xmm11,%%xmm11						\n\t"\
			"addpd	%%xmm1,%%xmm3						\n\t		addpd	%%xmm12,%%xmm8 						\n\t"\
			"addpd	%%xmm5,%%xmm6						\n\t		addpd	%%xmm14,%%xmm10						\n\t"\
			"movaps	%%xmm2,    (%%rdx)					\n\t		addpd	%%xmm13,%%xmm9 						\n\t"\
			"movaps	%%xmm7,    (%%rax)					\n\t		addpd	%%xmm15,%%xmm11						\n\t"\
			"movaps	%%xmm3,0x10(%%rdx)					\n\t		/* Output addresses: swap c/d: */			\n\t"\
			"movaps	%%xmm6,0x10(%%rbx)					\n\t		addq	%%r8,%%rax /* add48 */				\n\t"\
		"/*	Block 8: add1,0,2,3 = a[j1+p24]+p0-3 -> r07: */\n\t		addq	%%r8,%%rbx							\n\t"\
			"subq	$0x20,%%rsi	/* r07 */				\n\t		addq	%%r8,%%rcx							\n\t"\
			"movaps	     (%%rsi),%%xmm4					\n\t		addq	%%r8,%%rdx							\n\t"\
			"movaps	0x1a0(%%rsi),%%xmm6					\n\t		subpd	%%xmm10,%%xmm8 						\n\t"\
			"movaps	0x010(%%rsi),%%xmm5					\n\t		subpd	%%xmm15,%%xmm12						\n\t"\
			"movaps	0x1b0(%%rsi),%%xmm7					\n\t		subpd	%%xmm11,%%xmm9 						\n\t"\
			"movaps	0x340(%%rsi),%%xmm0					\n\t		subpd	%%xmm14,%%xmm13						\n\t"\
			"movaps	0x4e0(%%rsi),%%xmm2					\n\t		movaps	%%xmm8 ,    (%%rbx)					\n\t"\
			"movaps	0x350(%%rsi),%%xmm1					\n\t		movaps	%%xmm12,    (%%rdx)					\n\t"\
			"movaps	0x4f0(%%rsi),%%xmm3					\n\t		movaps	%%xmm9 ,0x10(%%rbx)					\n\t"\
			"subpd	%%xmm0,%%xmm4						\n\t		movaps	%%xmm13,0x10(%%rcx)					\n\t"\
			"subpd	%%xmm2,%%xmm6						\n\t		addpd	%%xmm10,%%xmm10						\n\t"\
			"subpd	%%xmm1,%%xmm5						\n\t		addpd	%%xmm15,%%xmm15						\n\t"\
			"subpd	%%xmm3,%%xmm7						\n\t		addpd	%%xmm11,%%xmm11						\n\t"\
			"addpd	%%xmm0,%%xmm0						\n\t		addpd	%%xmm14,%%xmm14						\n\t"\
			"addpd	%%xmm2,%%xmm2						\n\t		addpd	%%xmm8 ,%%xmm10						\n\t"\
			"addpd	%%xmm1,%%xmm1						\n\t		addpd	%%xmm12,%%xmm15						\n\t"\
			"addpd	%%xmm3,%%xmm3						\n\t		addpd	%%xmm9 ,%%xmm11						\n\t"\
			"addpd	%%xmm4,%%xmm0						\n\t		addpd	%%xmm13,%%xmm14						\n\t"\
			"addpd	%%xmm6,%%xmm2						\n\t		movaps	%%xmm10,    (%%rax)					\n\t"\
			"addpd	%%xmm5,%%xmm1						\n\t		movaps	%%xmm15,    (%%rcx)					\n\t"\
			"addpd	%%xmm7,%%xmm3						\n\t		movaps	%%xmm11,0x10(%%rax)					\n\t"\
			"/* Output addresses: swap a/b: */			\n\t		movaps	%%xmm14,0x10(%%rdx)					\n\t"\
			"subq	%%rdi,%%rax /* add24 */				\n\t"\
			"subq	%%rdi,%%rbx							\n\t"\
			"subq	%%rdi,%%rcx							\n\t"\
			"subq	%%rdi,%%rdx							\n\t"\
			"subpd	%%xmm2,%%xmm0						\n\t"\
			"subpd	%%xmm7,%%xmm4						\n\t"\
			"subpd	%%xmm3,%%xmm1						\n\t"\
			"subpd	%%xmm6,%%xmm5						\n\t"\
			"movaps	%%xmm0,    (%%rax)					\n\t"\
			"movaps	%%xmm4,    (%%rcx)					\n\t"\
			"movaps	%%xmm1,0x10(%%rax)					\n\t"\
			"movaps	%%xmm5,0x10(%%rdx)					\n\t"\
			"addpd	%%xmm2,%%xmm2						\n\t"\
			"addpd	%%xmm7,%%xmm7						\n\t"\
			"addpd	%%xmm3,%%xmm3						\n\t"\
			"addpd	%%xmm6,%%xmm6						\n\t"\
			"addpd	%%xmm0,%%xmm2						\n\t"\
			"addpd	%%xmm4,%%xmm7						\n\t"\
			"addpd	%%xmm1,%%xmm3						\n\t"\
			"addpd	%%xmm5,%%xmm6						\n\t"\
			"movaps	%%xmm2,    (%%rbx)					\n\t"\
			"movaps	%%xmm7,    (%%rdx)					\n\t"\
			"movaps	%%xmm3,0x10(%%rbx)					\n\t"\
			"movaps	%%xmm6,0x10(%%rcx)					\n\t"\
			:					/* outputs: none */\
			: [__add0] "m" (add0)	/* All inputs from memory addresses here */\
			 ,[__p01] "m" (p01)\
			 ,[__p02] "m" (p02)\
			 ,[__p03] "m" (p03)\
			 ,[__p04] "m" (p04)\
			 ,[__r00] "m" (r00)\
			: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r8","r9","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"		/* Clobbered registers */\
	);\
	}

#endif	/* radix52_ditN_cy_dif1_gcc_h_included */

