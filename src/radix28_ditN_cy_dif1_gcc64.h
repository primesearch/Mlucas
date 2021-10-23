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

/*******************************************************************************
   We now include this header file if it was not included before.
*******************************************************************************/
#ifndef radix28_ditN_cy_dif1_gcc_h_included
#define radix28_ditN_cy_dif1_gcc_h_included

#ifdef USE_AVX	// AVX and AVX2 both use 256-bit registers

	#define	SSE2_RADIX28_DIT_NOTWIDDLE(Xadd,Xp01,Xp02,Xp03,Xp04,Xp08,Xp12,Xp16,Xp20,Xp24,Xout,Xcc0)\
	{\
	__asm__ volatile (\
	"/* 1:SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add00+p0132, s1p00,03,02,01): */	\n\t"\
		"movq	%[__out],%%rsi												\n\t"\
		"movq	%[__add],%%rax												\n\t"\
		"movslq	%[__p01],%%rbx												\n\t"\
		"movslq	%[__p02],%%rcx												\n\t"\
		"movslq	%[__p03],%%rdx												\n\t"\
		"shlq	$3,%%rbx													\n\t"\
		"shlq	$3,%%rcx													\n\t"\
		"shlq	$3,%%rdx													\n\t"\
		"addq	%%rax,%%rbx													\n\t"\
		"addq	%%rax,%%rcx													\n\t"\
		"addq	%%rax,%%rdx													\n\t"\
		"/* ecx <-> edx */													\n\t"\
		"vmovaps	     (%%rax),%%ymm0										\n\t"\
		"vmovaps	     (%%rdx),%%ymm4										\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1										\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm5										\n\t"\
		"vmovaps	%%ymm0,%%ymm2											\n\t"\
		"vmovaps	%%ymm4,%%ymm6											\n\t"\
		"vmovaps	%%ymm1,%%ymm3											\n\t"\
		"vmovaps	%%ymm5,%%ymm7											\n\t"\
		"vaddpd	     (%%rbx),%%ymm0,%%ymm0									\n\t"\
		"vaddpd	     (%%rcx),%%ymm4,%%ymm4									\n\t"\
		"vaddpd	0x020(%%rbx),%%ymm1,%%ymm1									\n\t"\
		"vaddpd	0x020(%%rcx),%%ymm5,%%ymm5									\n\t"\
		"vsubpd	     (%%rbx),%%ymm2,%%ymm2									\n\t"\
		"vsubpd	     (%%rcx),%%ymm6,%%ymm6									\n\t/* 2:SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add12+p3210, s1p04,07,06,05r): */	\n\t"\
		"vsubpd	0x020(%%rbx),%%ymm3,%%ymm3									\n\t	movslq	%[__p12],%%rdi												\n\t"\
		"vsubpd	0x020(%%rcx),%%ymm7,%%ymm7									\n\t	shlq	$3,%%rdi													\n\t"\
		"/* Finish radix-4 butterfly and store results: */					\n\t	addq	%%rdi,%%rax		/* add0 + p12 */							\n\t"\
		"/* Swap _o1 <-> _o3: 0x040,0x060 <-> 0x0c0,0x0e0: */				\n\t	addq	%%rdi,%%rbx													\n\t"\
		"vsubpd	%%ymm4,%%ymm0,%%ymm0										\n\t	addq	%%rdi,%%rcx													\n\t"\
		"vsubpd	%%ymm7,%%ymm2,%%ymm2										\n\t	addq	%%rdi,%%rdx													\n\t"\
		"vsubpd	%%ymm5,%%ymm1,%%ymm1										\n\t	/* eax <-> edx, ebx <-> ecx */										\n\t"\
		"vsubpd	%%ymm6,%%ymm3,%%ymm3										\n\t	vmovaps	     (%%rdx),%%ymm8 										\n\t"\
		"vmovaps	%%ymm0,0x080(%%rsi)										\n\t	vmovaps	     (%%rbx),%%ymm12										\n\t"\
		"vmovaps	%%ymm2,0x040(%%rsi)										\n\t	vmovaps	0x020(%%rdx),%%ymm9 										\n\t"\
		"vmovaps	%%ymm1,0x0a0(%%rsi)										\n\t	vmovaps	0x020(%%rbx),%%ymm13										\n\t"\
		"vmovaps	%%ymm3,0x0e0(%%rsi)										\n\t	vmovaps	%%ymm8 ,%%ymm10											\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4										\n\t	vmovaps	%%ymm12,%%ymm14											\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7										\n\t	vmovaps	%%ymm9 ,%%ymm11											\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5										\n\t	vmovaps	%%ymm13,%%ymm15											\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6										\n\t	vaddpd	     (%%rcx),%%ymm8 ,%%ymm8 									\n\t"\
		"vaddpd	%%ymm0,%%ymm4,%%ymm4										\n\t	vaddpd	     (%%rax),%%ymm12,%%ymm12									\n\t"\
		"vaddpd	%%ymm2,%%ymm7,%%ymm7										\n\t	vaddpd	0x020(%%rcx),%%ymm9 ,%%ymm9 									\n\t"\
		"vaddpd	%%ymm1,%%ymm5,%%ymm5										\n\t	vaddpd	0x020(%%rax),%%ymm13,%%ymm13									\n\t"\
		"vaddpd	%%ymm3,%%ymm6,%%ymm6										\n\t	vsubpd	     (%%rcx),%%ymm10,%%ymm10									\n\t"\
		"vmovaps	%%ymm4,     (%%rsi)										\n\t	vsubpd	     (%%rax),%%ymm14,%%ymm14									\n\t"\
		"vmovaps	%%ymm7,0x0c0(%%rsi)										\n\t	vsubpd	0x020(%%rcx),%%ymm11,%%ymm11									\n\t"\
		"vmovaps	%%ymm5,0x020(%%rsi)										\n\t	vsubpd	0x020(%%rax),%%ymm15,%%ymm15									\n\t"\
		"vmovaps	%%ymm6,0x060(%%rsi)										\n\t	/* Finish radix-4 butterfly and store results: */					\n\t"\
	"/* 3:SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add24+p1023, s1p08,11,10,09r): */	\n\t	addq	$0x100,%%rsi												\n\t"\
		"addq	%%rdi,%%rax		/* add0 + p24 */							\n\t	vsubpd	%%ymm12,%%ymm8 ,%%ymm8 										\n\t"\
		"addq	%%rdi,%%rbx													\n\t	vsubpd	%%ymm15,%%ymm10,%%ymm10										\n\t"\
		"addq	%%rdi,%%rcx													\n\t	vsubpd	%%ymm13,%%ymm9 ,%%ymm9 										\n\t"\
		"addq	%%rdi,%%rdx													\n\t	vsubpd	%%ymm14,%%ymm11,%%ymm11										\n\t"\
		"/* eax <-> ebx */													\n\t	vmovaps	%%ymm8 ,0x080(%%rsi)										\n\t"\
		"vmovaps	     (%%rbx),%%ymm0										\n\t	vmovaps	%%ymm10,0x040(%%rsi)										\n\t"\
		"vmovaps	     (%%rcx),%%ymm4										\n\t	vmovaps	%%ymm9 ,0x0a0(%%rsi)										\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm1										\n\t	vmovaps	%%ymm11,0x0e0(%%rsi)										\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5										\n\t	vaddpd	%%ymm12,%%ymm12,%%ymm12										\n\t"\
		"vmovaps	%%ymm0,%%ymm2											\n\t	vaddpd	%%ymm15,%%ymm15,%%ymm15										\n\t"\
		"vmovaps	%%ymm4,%%ymm6											\n\t	vaddpd	%%ymm13,%%ymm13,%%ymm13										\n\t"\
		"vmovaps	%%ymm1,%%ymm3											\n\t	vaddpd	%%ymm14,%%ymm14,%%ymm14										\n\t"\
		"vmovaps	%%ymm5,%%ymm7											\n\t	vaddpd	%%ymm8 ,%%ymm12,%%ymm12										\n\t"\
		"vaddpd	     (%%rax),%%ymm0,%%ymm0									\n\t	vaddpd	%%ymm10,%%ymm15,%%ymm15										\n\t"\
		"vaddpd	     (%%rdx),%%ymm4,%%ymm4									\n\t	vaddpd	%%ymm9 ,%%ymm13,%%ymm13										\n\t"\
		"vaddpd	0x020(%%rax),%%ymm1,%%ymm1									\n\t	vaddpd	%%ymm11,%%ymm14,%%ymm14										\n\t"\
		"vaddpd	0x020(%%rdx),%%ymm5,%%ymm5									\n\t	vmovaps	%%ymm12,     (%%rsi)										\n\t"\
		"vsubpd	     (%%rax),%%ymm2,%%ymm2									\n\t	vmovaps	%%ymm15,0x0c0(%%rsi)										\n\t"\
		"vsubpd	     (%%rdx),%%ymm6,%%ymm6									\n\t	vmovaps	%%ymm13,0x020(%%rsi)										\n\t"\
		"vsubpd	0x020(%%rax),%%ymm3,%%ymm3									\n\t	vmovaps	%%ymm14,0x060(%%rsi)										\n\t"\
		"vsubpd	0x020(%%rdx),%%ymm7,%%ymm7									\n\t/* 4:SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add08+p1023, s1p12,15,14,13r): */	\n\t"\
		"/* Finish radix-4 butterfly and store results: */					\n\t	movslq	%[__p16],%%rdi												\n\t"\
		"addq	$0x100,%%rsi												\n\t	shlq	$3,%%rdi													\n\t"\
		"vsubpd	%%ymm4,%%ymm0,%%ymm0										\n\t	subq	%%rdi,%%rax		/* add0 + p08 */							\n\t"\
		"vsubpd	%%ymm7,%%ymm2,%%ymm2										\n\t	subq	%%rdi,%%rbx													\n\t"\
		"vsubpd	%%ymm5,%%ymm1,%%ymm1										\n\t	subq	%%rdi,%%rcx													\n\t"\
		"vsubpd	%%ymm6,%%ymm3,%%ymm3										\n\t	subq	%%rdi,%%rdx													\n\t"\
		"vmovaps	%%ymm0,0x080(%%rsi)										\n\t	/* eax <-> %%rbx */												\n\t"\
		"vmovaps	%%ymm2,0x040(%%rsi)										\n\t	vmovaps	     (%%rbx),%%ymm8 										\n\t"\
		"vmovaps	%%ymm1,0x0a0(%%rsi)										\n\t	vmovaps	     (%%rcx),%%ymm12										\n\t"\
		"vmovaps	%%ymm3,0x0e0(%%rsi)										\n\t	vmovaps	0x020(%%rbx),%%ymm9 										\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4										\n\t	vmovaps	0x020(%%rcx),%%ymm13										\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7										\n\t	vmovaps	%%ymm8 ,%%ymm10											\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5										\n\t	vmovaps	%%ymm12,%%ymm14											\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6										\n\t	vmovaps	%%ymm9 ,%%ymm11											\n\t"\
		"vaddpd	%%ymm0,%%ymm4,%%ymm4										\n\t	vmovaps	%%ymm13,%%ymm15											\n\t"\
		"vaddpd	%%ymm2,%%ymm7,%%ymm7										\n\t	vaddpd	     (%%rax),%%ymm8 ,%%ymm8 									\n\t"\
		"vaddpd	%%ymm1,%%ymm5,%%ymm5										\n\t	vaddpd	     (%%rdx),%%ymm12,%%ymm12									\n\t"\
		"vaddpd	%%ymm3,%%ymm6,%%ymm6										\n\t	vaddpd	0x020(%%rax),%%ymm9 ,%%ymm9 									\n\t"\
		"vmovaps	%%ymm4,     (%%rsi)										\n\t	vaddpd	0x020(%%rdx),%%ymm13,%%ymm13									\n\t"\
		"vmovaps	%%ymm7,0x0c0(%%rsi)										\n\t	vsubpd	     (%%rax),%%ymm10,%%ymm10									\n\t"\
		"vmovaps	%%ymm5,0x020(%%rsi)										\n\t	vsubpd	     (%%rdx),%%ymm14,%%ymm14									\n\t"\
		"vmovaps	%%ymm6,0x060(%%rsi)										\n\t	vsubpd	0x020(%%rax),%%ymm11,%%ymm11									\n\t"\
	"/* 5:SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add20+p2301, s1p16,19,18,17r): */	\n\t	vsubpd	0x020(%%rdx),%%ymm15,%%ymm15									\n\t"\
		"movslq	%[__p12],%%rdi												\n\t	/* Finish radix-4 butterfly and store results: */					\n\t"\
		"shlq	$3,%%rdi													\n\t	addq	$0x100,%%rsi												\n\t"\
		"addq	%%rdi,%%rax		/* add0 + p20 */							\n\t	vsubpd	%%ymm12,%%ymm8 ,%%ymm8 										\n\t"\
		"addq	%%rdi,%%rbx													\n\t	vsubpd	%%ymm15,%%ymm10,%%ymm10										\n\t"\
		"addq	%%rdi,%%rcx													\n\t	vsubpd	%%ymm13,%%ymm9 ,%%ymm9 										\n\t"\
		"addq	%%rdi,%%rdx													\n\t	vsubpd	%%ymm14,%%ymm11,%%ymm11										\n\t"\
		"/* eax <-> %%rcx, %%rbx <-> edx */									\n\t	vmovaps	%%ymm8 ,0x080(%%rsi)										\n\t"\
		"vmovaps	     (%%rcx),%%ymm0										\n\t	vmovaps	%%ymm10,0x040(%%rsi)										\n\t"\
		"vmovaps	     (%%rax),%%ymm4										\n\t	vmovaps	%%ymm9 ,0x0a0(%%rsi)										\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm1										\n\t	vmovaps	%%ymm11,0x0e0(%%rsi)										\n\t"\
		"vmovaps	0x020(%%rax),%%ymm5										\n\t	vaddpd	%%ymm12,%%ymm12,%%ymm12										\n\t"\
		"vmovaps	%%ymm0,%%ymm2											\n\t	vaddpd	%%ymm15,%%ymm15,%%ymm15										\n\t"\
		"vmovaps	%%ymm4,%%ymm6											\n\t	vaddpd	%%ymm13,%%ymm13,%%ymm13										\n\t"\
		"vmovaps	%%ymm1,%%ymm3											\n\t	vaddpd	%%ymm14,%%ymm14,%%ymm14										\n\t"\
		"vmovaps	%%ymm5,%%ymm7											\n\t	vaddpd	%%ymm8 ,%%ymm12,%%ymm12										\n\t"\
		"vaddpd	     (%%rdx),%%ymm0,%%ymm0									\n\t	vaddpd	%%ymm10,%%ymm15,%%ymm15										\n\t"\
		"vaddpd	     (%%rbx),%%ymm4,%%ymm4									\n\t	vaddpd	%%ymm9 ,%%ymm13,%%ymm13										\n\t"\
		"vaddpd	0x020(%%rdx),%%ymm1,%%ymm1									\n\t	vaddpd	%%ymm11,%%ymm14,%%ymm14										\n\t"\
		"vaddpd	0x020(%%rbx),%%ymm5,%%ymm5									\n\t	vmovaps	%%ymm12,     (%%rsi)										\n\t"\
		"vsubpd	     (%%rdx),%%ymm2,%%ymm2									\n\t	vmovaps	%%ymm15,0x0c0(%%rsi)										\n\t"\
		"vsubpd	     (%%rbx),%%ymm6,%%ymm6									\n\t	vmovaps	%%ymm13,0x020(%%rsi)										\n\t"\
		"vsubpd	0x020(%%rdx),%%ymm3,%%ymm3									\n\t	vmovaps	%%ymm14,0x060(%%rsi)										\n\t"\
		"vsubpd	0x020(%%rbx),%%ymm7,%%ymm7									\n\t/* 6:SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add04+p2301, s1p20,23,22,21r): */	\n\t"\
		"/* Finish radix-4 butterfly and store results: */					\n\t	movslq	%[__p16],%%rdi												\n\t"\
		"addq	$0x100,%%rsi												\n\t	shlq	$3,%%rdi													\n\t"\
		"vsubpd	%%ymm4,%%ymm0,%%ymm0										\n\t	subq	%%rdi,%%rax		/* add0 + p04 */							\n\t"\
		"vsubpd	%%ymm7,%%ymm2,%%ymm2										\n\t	subq	%%rdi,%%rbx													\n\t"\
		"vsubpd	%%ymm5,%%ymm1,%%ymm1										\n\t	subq	%%rdi,%%rcx													\n\t"\
		"vsubpd	%%ymm6,%%ymm3,%%ymm3										\n\t	subq	%%rdi,%%rdx													\n\t"\
		"vmovaps	%%ymm0,0x080(%%rsi)										\n\t	/* eax <-> %%rcx, %%rbx <-> edx */									\n\t"\
		"vmovaps	%%ymm2,0x040(%%rsi)										\n\t	vmovaps	     (%%rcx),%%ymm8 										\n\t"\
		"vmovaps	%%ymm1,0x0a0(%%rsi)										\n\t	vmovaps	     (%%rax),%%ymm12										\n\t"\
		"vmovaps	%%ymm3,0x0e0(%%rsi)										\n\t	vmovaps	0x020(%%rcx),%%ymm9 										\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4										\n\t	vmovaps	0x020(%%rax),%%ymm13										\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7										\n\t	vmovaps	%%ymm8 ,%%ymm10											\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5										\n\t	vmovaps	%%ymm12,%%ymm14											\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6										\n\t	vmovaps	%%ymm9 ,%%ymm11											\n\t"\
		"vaddpd	%%ymm0,%%ymm4,%%ymm4										\n\t	vmovaps	%%ymm13,%%ymm15											\n\t"\
		"vaddpd	%%ymm2,%%ymm7,%%ymm7										\n\t	vaddpd	     (%%rdx),%%ymm8 ,%%ymm8 									\n\t"\
		"vaddpd	%%ymm1,%%ymm5,%%ymm5										\n\t	vaddpd	     (%%rbx),%%ymm12,%%ymm12									\n\t"\
		"vaddpd	%%ymm3,%%ymm6,%%ymm6										\n\t	vaddpd	0x020(%%rdx),%%ymm9 ,%%ymm9 									\n\t"\
		"vmovaps	%%ymm4,     (%%rsi)										\n\t	vaddpd	0x020(%%rbx),%%ymm13,%%ymm13									\n\t"\
		"vmovaps	%%ymm7,0x0c0(%%rsi)										\n\t	vsubpd	     (%%rdx),%%ymm10,%%ymm10									\n\t"\
		"vmovaps	%%ymm5,0x020(%%rsi)										\n\t	vsubpd	     (%%rbx),%%ymm14,%%ymm14									\n\t"\
		"vmovaps	%%ymm6,0x060(%%rsi)										\n\t	vsubpd	0x020(%%rdx),%%ymm11,%%ymm11									\n\t"\
	"/* 7:SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add16+p0132, s1p24,27,26,25r): */	\n\t	vsubpd	0x020(%%rbx),%%ymm15,%%ymm15									\n\t"\
		"movslq	%[__p12],%%rdi												\n\t	/* Finish radix-4 butterfly and store results: */					\n\t"\
		"shlq	$3,%%rdi													\n\t	addq	$0x100,%%rsi												\n\t"\
		"addq	%%rdi,%%rax		/* add0 + p16 */							\n\t	vsubpd	%%ymm12,%%ymm8 ,%%ymm8 										\n\t"\
		"addq	%%rdi,%%rbx													\n\t	vsubpd	%%ymm15,%%ymm10,%%ymm10										\n\t"\
		"addq	%%rdi,%%rcx													\n\t	vsubpd	%%ymm13,%%ymm9 ,%%ymm9 										\n\t"\
		"addq	%%rdi,%%rdx													\n\t	vsubpd	%%ymm14,%%ymm11,%%ymm11										\n\t"\
		"/* ecx <-> %%rdx */												\n\t	vmovaps	%%ymm8 ,0x080(%%rsi)										\n\t"\
		"vmovaps	     (%%rax),%%ymm0										\n\t	vmovaps	%%ymm10,0x040(%%rsi)										\n\t"\
		"vmovaps	     (%%rdx),%%ymm4										\n\t	vmovaps	%%ymm9 ,0x0a0(%%rsi)										\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1										\n\t	vmovaps	%%ymm11,0x0e0(%%rsi)										\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm5										\n\t	vaddpd	%%ymm12,%%ymm12,%%ymm12										\n\t"\
		"vmovaps	%%ymm0,%%ymm2											\n\t	vaddpd	%%ymm15,%%ymm15,%%ymm15										\n\t"\
		"vmovaps	%%ymm4,%%ymm6											\n\t	vaddpd	%%ymm13,%%ymm13,%%ymm13										\n\t"\
		"vmovaps	%%ymm1,%%ymm3											\n\t	vaddpd	%%ymm14,%%ymm14,%%ymm14										\n\t"\
		"vmovaps	%%ymm5,%%ymm7											\n\t	vaddpd	%%ymm8 ,%%ymm12,%%ymm12										\n\t"\
		"vaddpd	     (%%rbx),%%ymm0,%%ymm0									\n\t	vaddpd	%%ymm10,%%ymm15,%%ymm15										\n\t"\
		"vaddpd	     (%%rcx),%%ymm4,%%ymm4									\n\t	vaddpd	%%ymm9 ,%%ymm13,%%ymm13										\n\t"\
		"vaddpd	0x020(%%rbx),%%ymm1,%%ymm1									\n\t	vaddpd	%%ymm11,%%ymm14,%%ymm14										\n\t"\
		"vaddpd	0x020(%%rcx),%%ymm5,%%ymm5									\n\t	vmovaps	%%ymm12,     (%%rsi)										\n\t"\
		"vsubpd	     (%%rbx),%%ymm2,%%ymm2									\n\t	vmovaps	%%ymm15,0x0c0(%%rsi)										\n\t"\
		"vsubpd	     (%%rcx),%%ymm6,%%ymm6									\n\t	vmovaps	%%ymm13,0x020(%%rsi)										\n\t"\
		"vsubpd	0x020(%%rbx),%%ymm3,%%ymm3									\n\t	vmovaps	%%ymm14,0x060(%%rsi)										\n\t"\
		"vsubpd	0x020(%%rcx),%%ymm7,%%ymm7									\n\t"\
		"/* Finish radix-4 butterfly and store results: */					\n\t"\
		"addq	$0x100,%%rsi												\n\t"\
		"vsubpd	%%ymm4,%%ymm0,%%ymm0										\n\t"\
		"vsubpd	%%ymm7,%%ymm2,%%ymm2										\n\t"\
		"vsubpd	%%ymm5,%%ymm1,%%ymm1										\n\t"\
		"vsubpd	%%ymm6,%%ymm3,%%ymm3										\n\t"\
		"vmovaps	%%ymm0,0x080(%%rsi)										\n\t"\
		"vmovaps	%%ymm2,0x040(%%rsi)										\n\t"\
		"vmovaps	%%ymm1,0x0a0(%%rsi)										\n\t"\
		"vmovaps	%%ymm3,0x0e0(%%rsi)										\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4										\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7										\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5										\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6										\n\t"\
		"vaddpd	%%ymm0,%%ymm4,%%ymm4										\n\t"\
		"vaddpd	%%ymm2,%%ymm7,%%ymm7										\n\t"\
		"vaddpd	%%ymm1,%%ymm5,%%ymm5										\n\t"\
		"vaddpd	%%ymm3,%%ymm6,%%ymm6										\n\t"\
		"vmovaps	%%ymm4,     (%%rsi)										\n\t"\
		"vmovaps	%%ymm7,0x0c0(%%rsi)										\n\t"\
		"vmovaps	%%ymm5,0x020(%%rsi)										\n\t"\
		"vmovaps	%%ymm6,0x060(%%rsi)										\n\t"\
		"/***************************************/							\n\t"\
		"/*...and now do 4 radix-7 transforms...*/							\n\t"\
		"/***************************************/							\n\t"\
	"/* SSE2_RADIX_07_DFT(00,04,08,12,16,20,24 -> 00,08,16,24,04,12,20): */	\n\t"\
		"movq	%[__out],%%rdi					\n\t"\
		"movq	%[__cc0],%%rax					\n\t"\
		"movq	$0x040,%%rbx					\n\t"\
		"movq	$0x080,%%rcx					\n\t"\
		"movq	$0x0c0,%%rdx					\n\t"\
		"movq	$0x100,%%rsi					\n\t"\
		"addq	%%rax,%%rbx						\n\t"\
		"addq	%%rax,%%rcx						\n\t"\
		"addq	%%rax,%%rdx						\n\t"\
		"addq	%%rax,%%rsi						\n\t	/*** Imaginary Parts: ***/					\n\t"\
		"vmovaps	0x100(%%rdi),%%ymm6			\n\t	vmovaps	0x120(%%rdi),%%ymm14				\n\t"\
		"vmovaps	0x600(%%rdi),%%ymm1			\n\t	vmovaps	0x620(%%rdi),%%ymm9 				\n\t"\
		"vmovaps	0x200(%%rdi),%%ymm5			\n\t	vmovaps	0x220(%%rdi),%%ymm13				\n\t"\
		"vmovaps	0x500(%%rdi),%%ymm2			\n\t	vmovaps	0x520(%%rdi),%%ymm10				\n\t"\
		"vmovaps	0x300(%%rdi),%%ymm4			\n\t	vmovaps	0x320(%%rdi),%%ymm12				\n\t"\
		"vmovaps	0x400(%%rdi),%%ymm3			\n\t	vmovaps	0x420(%%rdi),%%ymm11				\n\t"\
		"vsubpd	%%ymm1,%%ymm6,%%ymm6			\n\t	vsubpd	%%ymm9 ,%%ymm14,%%ymm14				\n\t"\
		"vaddpd	%%ymm1,%%ymm1,%%ymm1			\n\t	vaddpd	%%ymm9 ,%%ymm9 ,%%ymm9 				\n\t"\
		"vaddpd	%%ymm6,%%ymm1,%%ymm1			\n\t	vaddpd	%%ymm14,%%ymm9 ,%%ymm9 				\n\t"\
		"vsubpd	%%ymm2,%%ymm5,%%ymm5			\n\t	vsubpd	%%ymm10,%%ymm13,%%ymm13				\n\t"\
		"vaddpd	%%ymm2,%%ymm2,%%ymm2			\n\t	vaddpd	%%ymm10,%%ymm10,%%ymm10				\n\t"\
		"vaddpd	%%ymm5,%%ymm2,%%ymm2			\n\t	vaddpd	%%ymm13,%%ymm10,%%ymm10				\n\t"\
		"vmovaps	     (%%rdi),%%ymm0			\n\t	vmovaps	0x020(%%rdi),%%ymm8 				\n\t"\
		"vsubpd	%%ymm3,%%ymm4,%%ymm4			\n\t	vsubpd	%%ymm11,%%ymm12,%%ymm12				\n\t"\
		"vaddpd	%%ymm3,%%ymm3,%%ymm3			\n\t	vaddpd	%%ymm11,%%ymm11,%%ymm11				\n\t"\
		"vaddpd	%%ymm4,%%ymm3,%%ymm3			\n\t	vaddpd	%%ymm12,%%ymm11,%%ymm11				\n\t"\
		"/* Only alloc 2 scratch storage slots in calling code, so 16-reg version use 2 of the eventual outputs for Im parts */\n\t"\
		"vmovaps	%%ymm0,     (%%rsi)			\n\t	vmovaps	%%ymm8 ,0x100(%%rdi)				\n\t"\
		"vmovaps	%%ymm6,0x020(%%rsi)			\n\t	vmovaps	%%ymm14,0x200(%%rdi)				\n\t"\
		"vaddpd	      %%ymm1,%%ymm0,%%ymm0		\n\t	vaddpd	      %%ymm9 ,%%ymm8 ,%%ymm8 		\n\t"\
		"vmovaps	%%ymm5,%%ymm7				\n\t	vmovaps	%%ymm13,%%ymm15						\n\t"\
		"vaddpd	      %%ymm2,%%ymm3,%%ymm3		\n\t	vaddpd	      %%ymm10,%%ymm11,%%ymm11		\n\t"\
		"vsubpd	      %%ymm4,%%ymm5,%%ymm5		\n\t	vsubpd	      %%ymm12,%%ymm13,%%ymm13		\n\t"\
		"vsubpd	      %%ymm2,%%ymm1,%%ymm1		\n\t	vsubpd	      %%ymm10,%%ymm9 ,%%ymm9 		\n\t"\
		"vsubpd	      %%ymm7,%%ymm6,%%ymm6		\n\t	vsubpd	      %%ymm15,%%ymm14,%%ymm14		\n\t"\
		"vaddpd	      %%ymm2,%%ymm2,%%ymm2		\n\t	vaddpd	      %%ymm10,%%ymm10,%%ymm10		\n\t"\
		"vaddpd	      %%ymm7,%%ymm4,%%ymm4		\n\t	vaddpd	      %%ymm15,%%ymm12,%%ymm12		\n\t"\
		"vaddpd	      %%ymm3,%%ymm0,%%ymm0		\n\t	vaddpd	      %%ymm11,%%ymm8 ,%%ymm8 		\n\t"\
		"vaddpd	0x020(%%rsi),%%ymm5,%%ymm5		\n\t	vaddpd	0x200(%%rdi),%%ymm13,%%ymm13		\n\t"\
		"vsubpd	      %%ymm2,%%ymm3,%%ymm3		\n\t	vsubpd	      %%ymm10,%%ymm11,%%ymm11		\n\t"\
		"vmovaps	%%ymm4,%%ymm7				\n\t	vmovaps	%%ymm12,%%ymm15						\n\t"\
		"vmovaps	%%ymm0,      (%%rdi)		\n\t	vmovaps	%%ymm8 ,0x020(%%rdi)				\n\t"\
		"vsubpd	      %%ymm6,%%ymm4,%%ymm4		\n\t	vsubpd	      %%ymm14,%%ymm12,%%ymm12		\n\t"\
		"vmovaps	%%ymm1,%%ymm2				\n\t	vmovaps	%%ymm9 ,%%ymm10						\n\t"\
		"vsubpd	     (%%rsi),%%ymm0,%%ymm0		\n\t	vsubpd	0x100(%%rdi),%%ymm8 ,%%ymm8 		\n\t"\
		"vmulpd	0x020(%%rax),%%ymm5,%%ymm5		\n\t	vmulpd	0x020(%%rax),%%ymm13,%%ymm13		\n\t"\
		"vaddpd	      %%ymm3,%%ymm2,%%ymm2		\n\t	vaddpd	      %%ymm11,%%ymm10,%%ymm10		\n\t"\
		"vmulpd	     (%%rcx),%%ymm3,%%ymm3		\n\t	vmulpd	     (%%rcx),%%ymm11,%%ymm11		\n\t"\
		"vmulpd	0x020(%%rdx),%%ymm4,%%ymm4		\n\t	vmulpd	0x020(%%rdx),%%ymm12,%%ymm12		\n\t"\
		"vmulpd	     (%%rbx),%%ymm1,%%ymm1		\n\t	vmulpd	     (%%rbx),%%ymm9 ,%%ymm9 		\n\t"\
		"vmulpd	0x020(%%rbx),%%ymm6,%%ymm6		\n\t	vmulpd	0x020(%%rbx),%%ymm14,%%ymm14		\n\t"\
		"vmulpd	     (%%rax),%%ymm0,%%ymm0		\n\t	vmulpd	     (%%rax),%%ymm8 ,%%ymm8 		\n\t"\
		"vmulpd	0x020(%%rcx),%%ymm7,%%ymm7		\n\t	vmulpd	0x020(%%rcx),%%ymm15,%%ymm15		\n\t"\
		"vmulpd	     (%%rdx),%%ymm2,%%ymm2		\n\t	vmulpd	     (%%rdx),%%ymm10,%%ymm10		\n\t"\
		"vaddpd	     (%%rdi),%%ymm0,%%ymm0		\n\t	vaddpd	0x020(%%rdi),%%ymm8 ,%%ymm8 		\n\t"\
		"vaddpd	      %%ymm4,%%ymm6,%%ymm6		\n\t	vaddpd	      %%ymm12,%%ymm14,%%ymm14		\n\t"\
		"vsubpd	      %%ymm2,%%ymm1,%%ymm1		\n\t	vsubpd	      %%ymm10,%%ymm9 ,%%ymm9 		\n\t"\
		"vsubpd	      %%ymm7,%%ymm4,%%ymm4		\n\t	vsubpd	      %%ymm15,%%ymm12,%%ymm12		\n\t"\
		"vsubpd	      %%ymm2,%%ymm3,%%ymm3		\n\t	vsubpd	      %%ymm10,%%ymm11,%%ymm11		\n\t"\
		"vmovaps	%%ymm0,%%ymm2				\n\t	vmovaps	%%ymm8 ,%%ymm10						\n\t"\
		"vmovaps	%%ymm5,%%ymm7				\n\t	vmovaps	%%ymm13,%%ymm15						\n\t"\
		"vaddpd	      %%ymm1,%%ymm0,%%ymm0		\n\t	vaddpd	      %%ymm9 ,%%ymm8 ,%%ymm8 		\n\t"\
		"vaddpd	      %%ymm6,%%ymm5,%%ymm5		\n\t	vaddpd	      %%ymm14,%%ymm13,%%ymm13		\n\t"\
		"vaddpd	      %%ymm3,%%ymm1,%%ymm1		\n\t	vaddpd	      %%ymm11,%%ymm9 ,%%ymm9 		\n\t"\
		"vaddpd	      %%ymm4,%%ymm6,%%ymm6		\n\t	vaddpd	      %%ymm12,%%ymm14,%%ymm14		\n\t"\
		"vaddpd	      %%ymm2,%%ymm3,%%ymm3		\n\t	vaddpd	      %%ymm10,%%ymm11,%%ymm11		\n\t"\
		"vaddpd	      %%ymm7,%%ymm4,%%ymm4		\n\t	vaddpd	      %%ymm15,%%ymm12,%%ymm12		\n\t"\
		"vsubpd	      %%ymm1,%%ymm2,%%ymm2		\n\t	vsubpd	      %%ymm9 ,%%ymm10,%%ymm10		\n\t"\
		"vsubpd	      %%ymm6,%%ymm7,%%ymm7		\n\t	vsubpd	      %%ymm14,%%ymm15,%%ymm15		\n\t"\
		"/* ymm1,6,9,14 free ... Note the order reversal on the 3rd pair of outputs: */\n\t"\
		"vsubpd	%%ymm13,%%ymm0 ,%%ymm0		\n\t	vsubpd	%%ymm15,%%ymm2 ,%%ymm2	\n\t	vsubpd	%%ymm12,%%ymm3 ,%%ymm3 		\n\t"\
		"vsubpd	%%ymm5 ,%%ymm8 ,%%ymm8 		\n\t	vsubpd	%%ymm7 ,%%ymm10,%%ymm10	\n\t	vsubpd	%%ymm4 ,%%ymm11,%%ymm11		\n\t"\
		"vaddpd	%%ymm13,%%ymm13,%%ymm13		\n\t	vaddpd	%%ymm15,%%ymm15,%%ymm15	\n\t	vaddpd	%%ymm12,%%ymm12,%%ymm12		\n\t"\
		"vaddpd	%%ymm5 ,%%ymm5 ,%%ymm5		\n\t	vaddpd	%%ymm7 ,%%ymm7 ,%%ymm7	\n\t	vaddpd	%%ymm4 ,%%ymm4 ,%%ymm4 		\n\t"\
		"vaddpd	%%ymm0 ,%%ymm13,%%ymm13		\n\t	vaddpd	%%ymm2 ,%%ymm15,%%ymm15	\n\t	vaddpd	%%ymm3 ,%%ymm12,%%ymm12		\n\t"\
		"vaddpd	%%ymm8 ,%%ymm5 ,%%ymm5		\n\t	vaddpd	%%ymm10,%%ymm7 ,%%ymm7	\n\t	vaddpd	%%ymm11,%%ymm4 ,%%ymm4 		\n\t"\
		"vmovaps	%%ymm0 ,0x200(%%rdi)	\n\t	vmovaps	%%ymm2 ,0x400(%%rdi)	\n\t	vmovaps	%%ymm3 ,0x100(%%rdi)		\n\t"\
		"vmovaps	%%ymm8 ,0x520(%%rdi)	\n\t	vmovaps	%%ymm10,0x320(%%rdi)	\n\t	vmovaps	%%ymm11,0x620(%%rdi)		\n\t"\
		"vmovaps	%%ymm13,0x500(%%rdi)	\n\t	vmovaps	%%ymm15,0x300(%%rdi)	\n\t	vmovaps	%%ymm12,0x600(%%rdi)		\n\t"\
		"vmovaps	%%ymm5 ,0x220(%%rdi)	\n\t	vmovaps	%%ymm7 ,0x420(%%rdi)	\n\t	vmovaps	%%ymm4 ,0x120(%%rdi)		\n\t"\
	"/* SSE2_RADIX_07_DFT(01,05,09,13,17,21,25 -> 21,01,09,17,25,05,13): */	\n\t"\
		"addq	$0x040,%%rdi					\n\t	/*** Imaginary Parts: ***/					\n\t"\
		"vmovaps	0x100(%%rdi),%%ymm6			\n\t	vmovaps	0x120(%%rdi),%%ymm14				\n\t"\
		"vmovaps	0x600(%%rdi),%%ymm1			\n\t	vmovaps	0x620(%%rdi),%%ymm9 				\n\t"\
		"vmovaps	0x200(%%rdi),%%ymm5			\n\t	vmovaps	0x220(%%rdi),%%ymm13				\n\t"\
		"vmovaps	0x500(%%rdi),%%ymm2			\n\t	vmovaps	0x520(%%rdi),%%ymm10				\n\t"\
		"vmovaps	0x300(%%rdi),%%ymm4			\n\t	vmovaps	0x320(%%rdi),%%ymm12				\n\t"\
		"vmovaps	0x400(%%rdi),%%ymm3			\n\t	vmovaps	0x420(%%rdi),%%ymm11				\n\t"\
		"vsubpd	%%ymm1,%%ymm6,%%ymm6			\n\t	vsubpd	%%ymm9 ,%%ymm14,%%ymm14				\n\t"\
		"vaddpd	%%ymm1,%%ymm1,%%ymm1			\n\t	vaddpd	%%ymm9 ,%%ymm9 ,%%ymm9 				\n\t"\
		"vaddpd	%%ymm6,%%ymm1,%%ymm1			\n\t	vaddpd	%%ymm14,%%ymm9 ,%%ymm9 				\n\t"\
		"vsubpd	%%ymm2,%%ymm5,%%ymm5			\n\t	vsubpd	%%ymm10,%%ymm13,%%ymm13				\n\t"\
		"vaddpd	%%ymm2,%%ymm2,%%ymm2			\n\t	vaddpd	%%ymm10,%%ymm10,%%ymm10				\n\t"\
		"vaddpd	%%ymm5,%%ymm2,%%ymm2			\n\t	vaddpd	%%ymm13,%%ymm10,%%ymm10				\n\t"\
		"vmovaps	     (%%rdi),%%ymm0			\n\t	vmovaps	0x020(%%rdi),%%ymm8 				\n\t"\
		"vsubpd	%%ymm3,%%ymm4,%%ymm4			\n\t	vsubpd	%%ymm11,%%ymm12,%%ymm12				\n\t"\
		"vaddpd	%%ymm3,%%ymm3,%%ymm3			\n\t	vaddpd	%%ymm11,%%ymm11,%%ymm11				\n\t"\
		"vaddpd	%%ymm4,%%ymm3,%%ymm3			\n\t	vaddpd	%%ymm12,%%ymm11,%%ymm11				\n\t"\
		"/* Only alloc 2 scratch storage slots in calling code, so 16-reg version use 2 of the eventual outputs for Im parts */\n\t"\
		"vmovaps	%%ymm0,     (%%rsi)			\n\t	vmovaps	%%ymm8 ,0x100(%%rdi)				\n\t"\
		"vmovaps	%%ymm6,0x020(%%rsi)			\n\t	vmovaps	%%ymm14,0x200(%%rdi)				\n\t"\
		"vaddpd	      %%ymm1,%%ymm0,%%ymm0		\n\t	vaddpd	      %%ymm9 ,%%ymm8 ,%%ymm8 		\n\t"\
		"vmovaps	%%ymm5,%%ymm7				\n\t	vmovaps	%%ymm13,%%ymm15						\n\t"\
		"vaddpd	      %%ymm2,%%ymm3,%%ymm3		\n\t	vaddpd	      %%ymm10,%%ymm11,%%ymm11		\n\t"\
		"vsubpd	      %%ymm4,%%ymm5,%%ymm5		\n\t	vsubpd	      %%ymm12,%%ymm13,%%ymm13		\n\t"\
		"vsubpd	      %%ymm2,%%ymm1,%%ymm1		\n\t	vsubpd	      %%ymm10,%%ymm9 ,%%ymm9 		\n\t"\
		"vsubpd	      %%ymm7,%%ymm6,%%ymm6		\n\t	vsubpd	      %%ymm15,%%ymm14,%%ymm14		\n\t"\
		"vaddpd	      %%ymm2,%%ymm2,%%ymm2		\n\t	vaddpd	      %%ymm10,%%ymm10,%%ymm10		\n\t"\
		"vaddpd	      %%ymm7,%%ymm4,%%ymm4		\n\t	vaddpd	      %%ymm15,%%ymm12,%%ymm12		\n\t"\
		"vaddpd	      %%ymm3,%%ymm0,%%ymm0		\n\t	vaddpd	      %%ymm11,%%ymm8 ,%%ymm8 		\n\t"\
		"vaddpd	0x020(%%rsi),%%ymm5,%%ymm5		\n\t	vaddpd	0x200(%%rdi),%%ymm13,%%ymm13		\n\t"\
		"vsubpd	      %%ymm2,%%ymm3,%%ymm3		\n\t	vsubpd	      %%ymm10,%%ymm11,%%ymm11		\n\t"\
		"vmovaps	%%ymm4,%%ymm7				\n\t	vmovaps	%%ymm12,%%ymm15						\n\t"\
		"vmovaps	%%ymm0, 0x500(%%rdi)		\n\t	vmovaps	%%ymm8 ,0x520(%%rdi)				\n\t"\
		"vsubpd	      %%ymm6,%%ymm4,%%ymm4		\n\t	vsubpd	      %%ymm14,%%ymm12,%%ymm12		\n\t"\
		"vmovaps	%%ymm1,%%ymm2				\n\t	vmovaps	%%ymm9 ,%%ymm10						\n\t"\
		"vsubpd	     (%%rsi),%%ymm0,%%ymm0		\n\t	vsubpd	0x100(%%rdi),%%ymm8 ,%%ymm8 		\n\t"\
		"vmulpd	0x020(%%rax),%%ymm5,%%ymm5		\n\t	vmulpd	0x020(%%rax),%%ymm13,%%ymm13		\n\t"\
		"vaddpd	      %%ymm3,%%ymm2,%%ymm2		\n\t	vaddpd	      %%ymm11,%%ymm10,%%ymm10		\n\t"\
		"vmulpd	     (%%rcx),%%ymm3,%%ymm3		\n\t	vmulpd	     (%%rcx),%%ymm11,%%ymm11		\n\t"\
		"vmulpd	0x020(%%rdx),%%ymm4,%%ymm4		\n\t	vmulpd	0x020(%%rdx),%%ymm12,%%ymm12		\n\t"\
		"vmulpd	     (%%rbx),%%ymm1,%%ymm1		\n\t	vmulpd	     (%%rbx),%%ymm9 ,%%ymm9 		\n\t"\
		"vmulpd	0x020(%%rbx),%%ymm6,%%ymm6		\n\t	vmulpd	0x020(%%rbx),%%ymm14,%%ymm14		\n\t"\
		"vmulpd	     (%%rax),%%ymm0,%%ymm0		\n\t	vmulpd	     (%%rax),%%ymm8 ,%%ymm8 		\n\t"\
		"vmulpd	0x020(%%rcx),%%ymm7,%%ymm7		\n\t	vmulpd	0x020(%%rcx),%%ymm15,%%ymm15		\n\t"\
		"vmulpd	     (%%rdx),%%ymm2,%%ymm2		\n\t	vmulpd	     (%%rdx),%%ymm10,%%ymm10		\n\t"\
		"vaddpd	0x500(%%rdi),%%ymm0,%%ymm0		\n\t	vaddpd	0x520(%%rdi),%%ymm8 ,%%ymm8 		\n\t"\
		"vaddpd	      %%ymm4,%%ymm6,%%ymm6		\n\t	vaddpd	      %%ymm12,%%ymm14,%%ymm14		\n\t"\
		"vsubpd	      %%ymm2,%%ymm1,%%ymm1		\n\t	vsubpd	      %%ymm10,%%ymm9 ,%%ymm9 		\n\t"\
		"vsubpd	      %%ymm7,%%ymm4,%%ymm4		\n\t	vsubpd	      %%ymm15,%%ymm12,%%ymm12		\n\t"\
		"vsubpd	      %%ymm2,%%ymm3,%%ymm3		\n\t	vsubpd	      %%ymm10,%%ymm11,%%ymm11		\n\t"\
		"vmovaps	%%ymm0,%%ymm2				\n\t	vmovaps	%%ymm8 ,%%ymm10						\n\t"\
		"vmovaps	%%ymm5,%%ymm7				\n\t	vmovaps	%%ymm13,%%ymm15						\n\t"\
		"vaddpd	      %%ymm1,%%ymm0,%%ymm0		\n\t	vaddpd	      %%ymm9 ,%%ymm8 ,%%ymm8 		\n\t"\
		"vaddpd	      %%ymm6,%%ymm5,%%ymm5		\n\t	vaddpd	      %%ymm14,%%ymm13,%%ymm13		\n\t"\
		"vaddpd	      %%ymm3,%%ymm1,%%ymm1		\n\t	vaddpd	      %%ymm11,%%ymm9 ,%%ymm9 		\n\t"\
		"vaddpd	      %%ymm4,%%ymm6,%%ymm6		\n\t	vaddpd	      %%ymm12,%%ymm14,%%ymm14		\n\t"\
		"vaddpd	      %%ymm2,%%ymm3,%%ymm3		\n\t	vaddpd	      %%ymm10,%%ymm11,%%ymm11		\n\t"\
		"vaddpd	      %%ymm7,%%ymm4,%%ymm4		\n\t	vaddpd	      %%ymm15,%%ymm12,%%ymm12		\n\t"\
		"vsubpd	      %%ymm1,%%ymm2,%%ymm2		\n\t	vsubpd	      %%ymm9 ,%%ymm10,%%ymm10		\n\t"\
		"vsubpd	      %%ymm6,%%ymm7,%%ymm7		\n\t	vsubpd	      %%ymm14,%%ymm15,%%ymm15		\n\t"\
		"/* ymm1,6,9,14 free ... Note the order reversal on the 3rd pair of outputs: */\n\t"\
		"vsubpd	%%ymm13,%%ymm0 ,%%ymm0		\n\t	vsubpd	%%ymm15,%%ymm2 ,%%ymm2	\n\t	vsubpd	%%ymm12,%%ymm3 ,%%ymm3 		\n\t"\
		"vsubpd	%%ymm5 ,%%ymm8 ,%%ymm8 		\n\t	vsubpd	%%ymm7 ,%%ymm10,%%ymm10	\n\t	vsubpd	%%ymm4 ,%%ymm11,%%ymm11		\n\t"\
		"vaddpd	%%ymm13,%%ymm13,%%ymm13		\n\t	vaddpd	%%ymm15,%%ymm15,%%ymm15	\n\t	vaddpd	%%ymm12,%%ymm12,%%ymm12		\n\t"\
		"vaddpd	%%ymm5 ,%%ymm5 ,%%ymm5		\n\t	vaddpd	%%ymm7 ,%%ymm7 ,%%ymm7	\n\t	vaddpd	%%ymm4 ,%%ymm4 ,%%ymm4 		\n\t"\
		"vaddpd	%%ymm0 ,%%ymm13,%%ymm13		\n\t	vaddpd	%%ymm2 ,%%ymm15,%%ymm15	\n\t	vaddpd	%%ymm3 ,%%ymm12,%%ymm12		\n\t"\
		"vaddpd	%%ymm8 ,%%ymm5 ,%%ymm5		\n\t	vaddpd	%%ymm10,%%ymm7 ,%%ymm7	\n\t	vaddpd	%%ymm11,%%ymm4 ,%%ymm4 		\n\t"\
		"vmovaps	%%ymm0 ,     (%%rdi)	\n\t	vmovaps	%%ymm2 ,0x200(%%rdi)	\n\t	vmovaps	%%ymm3 ,0x600(%%rdi)		\n\t"\
		"vmovaps	%%ymm8 ,0x320(%%rdi)	\n\t	vmovaps	%%ymm10,0x120(%%rdi)	\n\t	vmovaps	%%ymm11,0x420(%%rdi)		\n\t"\
		"vmovaps	%%ymm13,0x300(%%rdi)	\n\t	vmovaps	%%ymm15,0x100(%%rdi)	\n\t	vmovaps	%%ymm12,0x400(%%rdi)		\n\t"\
		"vmovaps	%%ymm5 ,0x020(%%rdi)	\n\t	vmovaps	%%ymm7 ,0x220(%%rdi)	\n\t	vmovaps	%%ymm4 ,0x620(%%rdi)		\n\t"\
	"/* SSE2_RADIX_07_DFT(02,06,10,14,18,22,26 -> 14,22,02,10,18,26,06): */	\n\t"\
		"addq	$0x040,%%rdi					\n\t	/*** Imaginary Parts: ***/					\n\t"\
		"vmovaps	0x100(%%rdi),%%ymm6			\n\t	vmovaps	0x120(%%rdi),%%ymm14				\n\t"\
		"vmovaps	0x600(%%rdi),%%ymm1			\n\t	vmovaps	0x620(%%rdi),%%ymm9 				\n\t"\
		"vmovaps	0x200(%%rdi),%%ymm5			\n\t	vmovaps	0x220(%%rdi),%%ymm13				\n\t"\
		"vmovaps	0x500(%%rdi),%%ymm2			\n\t	vmovaps	0x520(%%rdi),%%ymm10				\n\t"\
		"vmovaps	0x300(%%rdi),%%ymm4			\n\t	vmovaps	0x320(%%rdi),%%ymm12				\n\t"\
		"vmovaps	0x400(%%rdi),%%ymm3			\n\t	vmovaps	0x420(%%rdi),%%ymm11				\n\t"\
		"vsubpd	%%ymm1,%%ymm6,%%ymm6			\n\t	vsubpd	%%ymm9 ,%%ymm14,%%ymm14				\n\t"\
		"vaddpd	%%ymm1,%%ymm1,%%ymm1			\n\t	vaddpd	%%ymm9 ,%%ymm9 ,%%ymm9 				\n\t"\
		"vaddpd	%%ymm6,%%ymm1,%%ymm1			\n\t	vaddpd	%%ymm14,%%ymm9 ,%%ymm9 				\n\t"\
		"vsubpd	%%ymm2,%%ymm5,%%ymm5			\n\t	vsubpd	%%ymm10,%%ymm13,%%ymm13				\n\t"\
		"vaddpd	%%ymm2,%%ymm2,%%ymm2			\n\t	vaddpd	%%ymm10,%%ymm10,%%ymm10				\n\t"\
		"vaddpd	%%ymm5,%%ymm2,%%ymm2			\n\t	vaddpd	%%ymm13,%%ymm10,%%ymm10				\n\t"\
		"vmovaps	     (%%rdi),%%ymm0			\n\t	vmovaps	0x020(%%rdi),%%ymm8 				\n\t"\
		"vsubpd	%%ymm3,%%ymm4,%%ymm4			\n\t	vsubpd	%%ymm11,%%ymm12,%%ymm12				\n\t"\
		"vaddpd	%%ymm3,%%ymm3,%%ymm3			\n\t	vaddpd	%%ymm11,%%ymm11,%%ymm11				\n\t"\
		"vaddpd	%%ymm4,%%ymm3,%%ymm3			\n\t	vaddpd	%%ymm12,%%ymm11,%%ymm11				\n\t"\
		"/* Only alloc 2 scratch storage slots in calling code, so 16-reg version use 2 of the eventual outputs for Im parts */\n\t"\
		"vmovaps	%%ymm0,     (%%rsi)			\n\t	vmovaps	%%ymm8 ,0x100(%%rdi)				\n\t"\
		"vmovaps	%%ymm6,0x020(%%rsi)			\n\t	vmovaps	%%ymm14,0x200(%%rdi)				\n\t"\
		"vaddpd	      %%ymm1,%%ymm0,%%ymm0		\n\t	vaddpd	      %%ymm9 ,%%ymm8 ,%%ymm8 		\n\t"\
		"vmovaps	%%ymm5,%%ymm7				\n\t	vmovaps	%%ymm13,%%ymm15						\n\t"\
		"vaddpd	      %%ymm2,%%ymm3,%%ymm3		\n\t	vaddpd	      %%ymm10,%%ymm11,%%ymm11		\n\t"\
		"vsubpd	      %%ymm4,%%ymm5,%%ymm5		\n\t	vsubpd	      %%ymm12,%%ymm13,%%ymm13		\n\t"\
		"vsubpd	      %%ymm2,%%ymm1,%%ymm1		\n\t	vsubpd	      %%ymm10,%%ymm9 ,%%ymm9 		\n\t"\
		"vsubpd	      %%ymm7,%%ymm6,%%ymm6		\n\t	vsubpd	      %%ymm15,%%ymm14,%%ymm14		\n\t"\
		"vaddpd	      %%ymm2,%%ymm2,%%ymm2		\n\t	vaddpd	      %%ymm10,%%ymm10,%%ymm10		\n\t"\
		"vaddpd	      %%ymm7,%%ymm4,%%ymm4		\n\t	vaddpd	      %%ymm15,%%ymm12,%%ymm12		\n\t"\
		"vaddpd	      %%ymm3,%%ymm0,%%ymm0		\n\t	vaddpd	      %%ymm11,%%ymm8 ,%%ymm8 		\n\t"\
		"vaddpd	0x020(%%rsi),%%ymm5,%%ymm5		\n\t	vaddpd	0x200(%%rdi),%%ymm13,%%ymm13		\n\t"\
		"vsubpd	      %%ymm2,%%ymm3,%%ymm3		\n\t	vsubpd	      %%ymm10,%%ymm11,%%ymm11		\n\t"\
		"vmovaps	%%ymm4,%%ymm7				\n\t	vmovaps	%%ymm12,%%ymm15						\n\t"\
		"vmovaps	%%ymm0, 0x300(%%rdi)		\n\t	vmovaps	%%ymm8 ,0x320(%%rdi)				\n\t"\
		"vsubpd	      %%ymm6,%%ymm4,%%ymm4		\n\t	vsubpd	      %%ymm14,%%ymm12,%%ymm12		\n\t"\
		"vmovaps	%%ymm1,%%ymm2				\n\t	vmovaps	%%ymm9 ,%%ymm10						\n\t"\
		"vsubpd	     (%%rsi),%%ymm0,%%ymm0		\n\t	vsubpd	0x100(%%rdi),%%ymm8 ,%%ymm8 		\n\t"\
		"vmulpd	0x020(%%rax),%%ymm5,%%ymm5		\n\t	vmulpd	0x020(%%rax),%%ymm13,%%ymm13		\n\t"\
		"vaddpd	      %%ymm3,%%ymm2,%%ymm2		\n\t	vaddpd	      %%ymm11,%%ymm10,%%ymm10		\n\t"\
		"vmulpd	     (%%rcx),%%ymm3,%%ymm3		\n\t	vmulpd	     (%%rcx),%%ymm11,%%ymm11		\n\t"\
		"vmulpd	0x020(%%rdx),%%ymm4,%%ymm4		\n\t	vmulpd	0x020(%%rdx),%%ymm12,%%ymm12		\n\t"\
		"vmulpd	     (%%rbx),%%ymm1,%%ymm1		\n\t	vmulpd	     (%%rbx),%%ymm9 ,%%ymm9 		\n\t"\
		"vmulpd	0x020(%%rbx),%%ymm6,%%ymm6		\n\t	vmulpd	0x020(%%rbx),%%ymm14,%%ymm14		\n\t"\
		"vmulpd	     (%%rax),%%ymm0,%%ymm0		\n\t	vmulpd	     (%%rax),%%ymm8 ,%%ymm8 		\n\t"\
		"vmulpd	0x020(%%rcx),%%ymm7,%%ymm7		\n\t	vmulpd	0x020(%%rcx),%%ymm15,%%ymm15		\n\t"\
		"vmulpd	     (%%rdx),%%ymm2,%%ymm2		\n\t	vmulpd	     (%%rdx),%%ymm10,%%ymm10		\n\t"\
		"vaddpd	0x300(%%rdi),%%ymm0,%%ymm0		\n\t	vaddpd	0x320(%%rdi),%%ymm8 ,%%ymm8 		\n\t"\
		"vaddpd	      %%ymm4,%%ymm6,%%ymm6		\n\t	vaddpd	      %%ymm12,%%ymm14,%%ymm14		\n\t"\
		"vsubpd	      %%ymm2,%%ymm1,%%ymm1		\n\t	vsubpd	      %%ymm10,%%ymm9 ,%%ymm9 		\n\t"\
		"vsubpd	      %%ymm7,%%ymm4,%%ymm4		\n\t	vsubpd	      %%ymm15,%%ymm12,%%ymm12		\n\t"\
		"vsubpd	      %%ymm2,%%ymm3,%%ymm3		\n\t	vsubpd	      %%ymm10,%%ymm11,%%ymm11		\n\t"\
		"vmovaps	%%ymm0,%%ymm2				\n\t	vmovaps	%%ymm8 ,%%ymm10						\n\t"\
		"vmovaps	%%ymm5,%%ymm7				\n\t	vmovaps	%%ymm13,%%ymm15						\n\t"\
		"vaddpd	      %%ymm1,%%ymm0,%%ymm0		\n\t	vaddpd	      %%ymm9 ,%%ymm8 ,%%ymm8 		\n\t"\
		"vaddpd	      %%ymm6,%%ymm5,%%ymm5		\n\t	vaddpd	      %%ymm14,%%ymm13,%%ymm13		\n\t"\
		"vaddpd	      %%ymm3,%%ymm1,%%ymm1		\n\t	vaddpd	      %%ymm11,%%ymm9 ,%%ymm9 		\n\t"\
		"vaddpd	      %%ymm4,%%ymm6,%%ymm6		\n\t	vaddpd	      %%ymm12,%%ymm14,%%ymm14		\n\t"\
		"vaddpd	      %%ymm2,%%ymm3,%%ymm3		\n\t	vaddpd	      %%ymm10,%%ymm11,%%ymm11		\n\t"\
		"vaddpd	      %%ymm7,%%ymm4,%%ymm4		\n\t	vaddpd	      %%ymm15,%%ymm12,%%ymm12		\n\t"\
		"vsubpd	      %%ymm1,%%ymm2,%%ymm2		\n\t	vsubpd	      %%ymm9 ,%%ymm10,%%ymm10		\n\t"\
		"vsubpd	      %%ymm6,%%ymm7,%%ymm7		\n\t	vsubpd	      %%ymm14,%%ymm15,%%ymm15		\n\t"\
		"/* ymm1,6,9,14 free ... Note the order reversal on the 3rd pair of outputs: */\n\t"\
		"vsubpd	%%ymm13,%%ymm0 ,%%ymm0		\n\t	vsubpd	%%ymm15,%%ymm2 ,%%ymm2	\n\t	vsubpd	%%ymm12,%%ymm3 ,%%ymm3 		\n\t"\
		"vsubpd	%%ymm5 ,%%ymm8 ,%%ymm8 		\n\t	vsubpd	%%ymm7 ,%%ymm10,%%ymm10	\n\t	vsubpd	%%ymm4 ,%%ymm11,%%ymm11		\n\t"\
		"vaddpd	%%ymm13,%%ymm13,%%ymm13		\n\t	vaddpd	%%ymm15,%%ymm15,%%ymm15	\n\t	vaddpd	%%ymm12,%%ymm12,%%ymm12		\n\t"\
		"vaddpd	%%ymm5 ,%%ymm5 ,%%ymm5		\n\t	vaddpd	%%ymm7 ,%%ymm7 ,%%ymm7	\n\t	vaddpd	%%ymm4 ,%%ymm4 ,%%ymm4 		\n\t"\
		"vaddpd	%%ymm0 ,%%ymm13,%%ymm13		\n\t	vaddpd	%%ymm2 ,%%ymm15,%%ymm15	\n\t	vaddpd	%%ymm3 ,%%ymm12,%%ymm12		\n\t"\
		"vaddpd	%%ymm8 ,%%ymm5 ,%%ymm5		\n\t	vaddpd	%%ymm10,%%ymm7 ,%%ymm7	\n\t	vaddpd	%%ymm11,%%ymm4 ,%%ymm4 		\n\t"\
		"vmovaps	%%ymm0 ,0x500(%%rdi)	\n\t	vmovaps	%%ymm2 ,     (%%rdi)	\n\t	vmovaps	%%ymm3 ,0x400(%%rdi)		\n\t"\
		"vmovaps	%%ymm8 ,0x120(%%rdi)	\n\t	vmovaps	%%ymm10,0x620(%%rdi)	\n\t	vmovaps	%%ymm11,0x220(%%rdi)		\n\t"\
		"vmovaps	%%ymm13,0x100(%%rdi)	\n\t	vmovaps	%%ymm15,0x600(%%rdi)	\n\t	vmovaps	%%ymm12,0x200(%%rdi)		\n\t"\
		"vmovaps	%%ymm5 ,0x520(%%rdi)	\n\t	vmovaps	%%ymm7 ,0x020(%%rdi)	\n\t	vmovaps	%%ymm4 ,0x420(%%rdi)		\n\t"\
	"/* SSE2_RADIX_07_DFT(03,07,11,15,19,23,27 -> 07,15,23,03,11,19,27): */	\n\t"\
		"addq	$0x040,%%rdi					\n\t	/*** Imaginary Parts: ***/					\n\t"\
		"vmovaps	0x100(%%rdi),%%ymm6			\n\t	vmovaps	0x120(%%rdi),%%ymm14				\n\t"\
		"vmovaps	0x600(%%rdi),%%ymm1			\n\t	vmovaps	0x620(%%rdi),%%ymm9 				\n\t"\
		"vmovaps	0x200(%%rdi),%%ymm5			\n\t	vmovaps	0x220(%%rdi),%%ymm13				\n\t"\
		"vmovaps	0x500(%%rdi),%%ymm2			\n\t	vmovaps	0x520(%%rdi),%%ymm10				\n\t"\
		"vmovaps	0x300(%%rdi),%%ymm4			\n\t	vmovaps	0x320(%%rdi),%%ymm12				\n\t"\
		"vmovaps	0x400(%%rdi),%%ymm3			\n\t	vmovaps	0x420(%%rdi),%%ymm11				\n\t"\
		"vsubpd	%%ymm1,%%ymm6,%%ymm6			\n\t	vsubpd	%%ymm9 ,%%ymm14,%%ymm14				\n\t"\
		"vaddpd	%%ymm1,%%ymm1,%%ymm1			\n\t	vaddpd	%%ymm9 ,%%ymm9 ,%%ymm9 				\n\t"\
		"vaddpd	%%ymm6,%%ymm1,%%ymm1			\n\t	vaddpd	%%ymm14,%%ymm9 ,%%ymm9 				\n\t"\
		"vsubpd	%%ymm2,%%ymm5,%%ymm5			\n\t	vsubpd	%%ymm10,%%ymm13,%%ymm13				\n\t"\
		"vaddpd	%%ymm2,%%ymm2,%%ymm2			\n\t	vaddpd	%%ymm10,%%ymm10,%%ymm10				\n\t"\
		"vaddpd	%%ymm5,%%ymm2,%%ymm2			\n\t	vaddpd	%%ymm13,%%ymm10,%%ymm10				\n\t"\
		"vmovaps	     (%%rdi),%%ymm0			\n\t	vmovaps	0x020(%%rdi),%%ymm8 				\n\t"\
		"vsubpd	%%ymm3,%%ymm4,%%ymm4			\n\t	vsubpd	%%ymm11,%%ymm12,%%ymm12				\n\t"\
		"vaddpd	%%ymm3,%%ymm3,%%ymm3			\n\t	vaddpd	%%ymm11,%%ymm11,%%ymm11				\n\t"\
		"vaddpd	%%ymm4,%%ymm3,%%ymm3			\n\t	vaddpd	%%ymm12,%%ymm11,%%ymm11				\n\t"\
		"/* Only alloc 2 scratch storage slots in calling code, so 16-reg version use 2 of the eventual outputs for Im parts */\n\t"\
		"vmovaps	%%ymm0,     (%%rsi)			\n\t	vmovaps	%%ymm8 ,0x300(%%rdi)				\n\t"\
		"vmovaps	%%ymm6,0x020(%%rsi)			\n\t	vmovaps	%%ymm14,0x200(%%rdi)				\n\t"\
		"vaddpd	      %%ymm1,%%ymm0,%%ymm0		\n\t	vaddpd	      %%ymm9 ,%%ymm8 ,%%ymm8 		\n\t"\
		"vmovaps	%%ymm5,%%ymm7				\n\t	vmovaps	%%ymm13,%%ymm15						\n\t"\
		"vaddpd	      %%ymm2,%%ymm3,%%ymm3		\n\t	vaddpd	      %%ymm10,%%ymm11,%%ymm11		\n\t"\
		"vsubpd	      %%ymm4,%%ymm5,%%ymm5		\n\t	vsubpd	      %%ymm12,%%ymm13,%%ymm13		\n\t"\
		"vsubpd	      %%ymm2,%%ymm1,%%ymm1		\n\t	vsubpd	      %%ymm10,%%ymm9 ,%%ymm9 		\n\t"\
		"vsubpd	      %%ymm7,%%ymm6,%%ymm6		\n\t	vsubpd	      %%ymm15,%%ymm14,%%ymm14		\n\t"\
		"vaddpd	      %%ymm2,%%ymm2,%%ymm2		\n\t	vaddpd	      %%ymm10,%%ymm10,%%ymm10		\n\t"\
		"vaddpd	      %%ymm7,%%ymm4,%%ymm4		\n\t	vaddpd	      %%ymm15,%%ymm12,%%ymm12		\n\t"\
		"vaddpd	      %%ymm3,%%ymm0,%%ymm0		\n\t	vaddpd	      %%ymm11,%%ymm8 ,%%ymm8 		\n\t"\
		"vaddpd	0x020(%%rsi),%%ymm5,%%ymm5		\n\t	vaddpd	0x200(%%rdi),%%ymm13,%%ymm13		\n\t"\
		"vsubpd	      %%ymm2,%%ymm3,%%ymm3		\n\t	vsubpd	      %%ymm10,%%ymm11,%%ymm11		\n\t"\
		"vmovaps	%%ymm4,%%ymm7				\n\t	vmovaps	%%ymm12,%%ymm15						\n\t"\
		"vmovaps	%%ymm0, 0x100(%%rdi)		\n\t	vmovaps	%%ymm8 ,0x120(%%rdi)				\n\t"\
		"vsubpd	      %%ymm6,%%ymm4,%%ymm4		\n\t	vsubpd	      %%ymm14,%%ymm12,%%ymm12		\n\t"\
		"vmovaps	%%ymm1,%%ymm2				\n\t	vmovaps	%%ymm9 ,%%ymm10						\n\t"\
		"vsubpd	     (%%rsi),%%ymm0,%%ymm0		\n\t	vsubpd	0x300(%%rdi),%%ymm8 ,%%ymm8 		\n\t"\
		"vmulpd	0x020(%%rax),%%ymm5,%%ymm5		\n\t	vmulpd	0x020(%%rax),%%ymm13,%%ymm13		\n\t"\
		"vaddpd	      %%ymm3,%%ymm2,%%ymm2		\n\t	vaddpd	      %%ymm11,%%ymm10,%%ymm10		\n\t"\
		"vmulpd	     (%%rcx),%%ymm3,%%ymm3		\n\t	vmulpd	     (%%rcx),%%ymm11,%%ymm11		\n\t"\
		"vmulpd	0x020(%%rdx),%%ymm4,%%ymm4		\n\t	vmulpd	0x020(%%rdx),%%ymm12,%%ymm12		\n\t"\
		"vmulpd	     (%%rbx),%%ymm1,%%ymm1		\n\t	vmulpd	     (%%rbx),%%ymm9 ,%%ymm9 		\n\t"\
		"vmulpd	0x020(%%rbx),%%ymm6,%%ymm6		\n\t	vmulpd	0x020(%%rbx),%%ymm14,%%ymm14		\n\t"\
		"vmulpd	     (%%rax),%%ymm0,%%ymm0		\n\t	vmulpd	     (%%rax),%%ymm8 ,%%ymm8 		\n\t"\
		"vmulpd	0x020(%%rcx),%%ymm7,%%ymm7		\n\t	vmulpd	0x020(%%rcx),%%ymm15,%%ymm15		\n\t"\
		"vmulpd	     (%%rdx),%%ymm2,%%ymm2		\n\t	vmulpd	     (%%rdx),%%ymm10,%%ymm10		\n\t"\
		"vaddpd	0x100(%%rdi),%%ymm0,%%ymm0		\n\t	vaddpd	0x120(%%rdi),%%ymm8 ,%%ymm8 		\n\t"\
		"vaddpd	      %%ymm4,%%ymm6,%%ymm6		\n\t	vaddpd	      %%ymm12,%%ymm14,%%ymm14		\n\t"\
		"vsubpd	      %%ymm2,%%ymm1,%%ymm1		\n\t	vsubpd	      %%ymm10,%%ymm9 ,%%ymm9 		\n\t"\
		"vsubpd	      %%ymm7,%%ymm4,%%ymm4		\n\t	vsubpd	      %%ymm15,%%ymm12,%%ymm12		\n\t"\
		"vsubpd	      %%ymm2,%%ymm3,%%ymm3		\n\t	vsubpd	      %%ymm10,%%ymm11,%%ymm11		\n\t"\
		"vmovaps	%%ymm0,%%ymm2				\n\t	vmovaps	%%ymm8 ,%%ymm10						\n\t"\
		"vmovaps	%%ymm5,%%ymm7				\n\t	vmovaps	%%ymm13,%%ymm15						\n\t"\
		"vaddpd	      %%ymm1,%%ymm0,%%ymm0		\n\t	vaddpd	      %%ymm9 ,%%ymm8 ,%%ymm8 		\n\t"\
		"vaddpd	      %%ymm6,%%ymm5,%%ymm5		\n\t	vaddpd	      %%ymm14,%%ymm13,%%ymm13		\n\t"\
		"vaddpd	      %%ymm3,%%ymm1,%%ymm1		\n\t	vaddpd	      %%ymm11,%%ymm9 ,%%ymm9 		\n\t"\
		"vaddpd	      %%ymm4,%%ymm6,%%ymm6		\n\t	vaddpd	      %%ymm12,%%ymm14,%%ymm14		\n\t"\
		"vaddpd	      %%ymm2,%%ymm3,%%ymm3		\n\t	vaddpd	      %%ymm10,%%ymm11,%%ymm11		\n\t"\
		"vaddpd	      %%ymm7,%%ymm4,%%ymm4		\n\t	vaddpd	      %%ymm15,%%ymm12,%%ymm12		\n\t"\
		"vsubpd	      %%ymm1,%%ymm2,%%ymm2		\n\t	vsubpd	      %%ymm9 ,%%ymm10,%%ymm10		\n\t"\
		"vsubpd	      %%ymm6,%%ymm7,%%ymm7		\n\t	vsubpd	      %%ymm14,%%ymm15,%%ymm15		\n\t"\
		"/* ymm1,6,9,14 free ... Note the order reversal on the 3rd pair of outputs: */\n\t"\
		"vsubpd	%%ymm13,%%ymm0 ,%%ymm0		\n\t	vsubpd	%%ymm15,%%ymm2 ,%%ymm2	\n\t	vsubpd	%%ymm12,%%ymm3 ,%%ymm3 		\n\t"\
		"vsubpd	%%ymm5 ,%%ymm8 ,%%ymm8 		\n\t	vsubpd	%%ymm7 ,%%ymm10,%%ymm10	\n\t	vsubpd	%%ymm4 ,%%ymm11,%%ymm11		\n\t"\
		"vaddpd	%%ymm13,%%ymm13,%%ymm13		\n\t	vaddpd	%%ymm15,%%ymm15,%%ymm15	\n\t	vaddpd	%%ymm12,%%ymm12,%%ymm12		\n\t"\
		"vaddpd	%%ymm5 ,%%ymm5 ,%%ymm5		\n\t	vaddpd	%%ymm7 ,%%ymm7 ,%%ymm7	\n\t	vaddpd	%%ymm4 ,%%ymm4 ,%%ymm4 		\n\t"\
		"vaddpd	%%ymm0 ,%%ymm13,%%ymm13		\n\t	vaddpd	%%ymm2 ,%%ymm15,%%ymm15	\n\t	vaddpd	%%ymm3 ,%%ymm12,%%ymm12		\n\t"\
		"vaddpd	%%ymm8 ,%%ymm5 ,%%ymm5		\n\t	vaddpd	%%ymm10,%%ymm7 ,%%ymm7	\n\t	vaddpd	%%ymm11,%%ymm4 ,%%ymm4 		\n\t"\
		"vmovaps	%%ymm0 ,0x300(%%rdi)	\n\t	vmovaps	%%ymm2 ,0x500(%%rdi)	\n\t	vmovaps	%%ymm3 ,0x200(%%rdi)		\n\t"\
		"vmovaps	%%ymm8 ,0x620(%%rdi)	\n\t	vmovaps	%%ymm10,0x420(%%rdi)	\n\t	vmovaps	%%ymm11,0x020(%%rdi)		\n\t"\
		"vmovaps	%%ymm13,0x600(%%rdi)	\n\t	vmovaps	%%ymm15,0x400(%%rdi)	\n\t	vmovaps	%%ymm12,     (%%rdi)		\n\t"\
		"vmovaps	%%ymm5 ,0x320(%%rdi)	\n\t	vmovaps	%%ymm7 ,0x520(%%rdi)	\n\t	vmovaps	%%ymm4 ,0x220(%%rdi)		\n\t"\
		:					/* outputs: none */\
		: [__add] "m" (Xadd)	/* All inputs from memory addresses here */\
		 ,[__p01] "m" (Xp01)\
		 ,[__p02] "m" (Xp02)\
		 ,[__p03] "m" (Xp03)\
		 ,[__p04] "m" (Xp04)\
		 ,[__p08] "m" (Xp08)\
		 ,[__p12] "m" (Xp12)\
		 ,[__p16] "m" (Xp16)\
		 ,[__p20] "m" (Xp20)\
		 ,[__p24] "m" (Xp24)\
		 ,[__out] "m" (Xout)\
		 ,[__cc0] "m" (Xcc0)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"		/* Clobbered registers */\
	);\
	}

	#define	SSE2_RADIX28_DIF_NOTWIDDLE(Xadd,Xp01,Xp02,Xp03,Xp04,Xp08,Xp12,Xp16,Xp20,Xp24,Xout,Xcc0)\
	{\
	__asm__ volatile (\
	"/* SSE2_RADIX_07_DFT(00,24,20,16,12,08,04r -> 00,04,08,12,16,20,24): */\n\t"\
		"movq	%[__out],%%rdi					\n\t"\
		"movq	%[__cc0],%%rax					\n\t"\
		"movq	$0x040,%%rbx					\n\t"\
		"movq	$0x080,%%rcx					\n\t"\
		"movq	$0x0c0,%%rdx					\n\t"\
		"movq	$0x100,%%rsi					\n\t"\
		"addq	%%rax,%%rbx						\n\t"\
		"addq	%%rax,%%rcx						\n\t"\
		"addq	%%rax,%%rdx						\n\t"\
		"addq	%%rax,%%rsi						\n\t	/*** Imaginary Parts: ***/					\n\t"\
		"vmovaps	0x600(%%rdi),%%ymm6			\n\t	vmovaps	0x620(%%rdi),%%ymm14				\n\t"\
		"vmovaps	0x100(%%rdi),%%ymm1			\n\t	vmovaps	0x120(%%rdi),%%ymm9 				\n\t"\
		"vmovaps	0x500(%%rdi),%%ymm5			\n\t	vmovaps	0x520(%%rdi),%%ymm13				\n\t"\
		"vmovaps	0x200(%%rdi),%%ymm2			\n\t	vmovaps	0x220(%%rdi),%%ymm10				\n\t"\
		"vmovaps	0x400(%%rdi),%%ymm4			\n\t	vmovaps	0x420(%%rdi),%%ymm12				\n\t"\
		"vmovaps	0x300(%%rdi),%%ymm3			\n\t	vmovaps	0x320(%%rdi),%%ymm11				\n\t"\
		"vsubpd	%%ymm1,%%ymm6,%%ymm6			\n\t	vsubpd	%%ymm9 ,%%ymm14,%%ymm14				\n\t"\
		"vaddpd	%%ymm1,%%ymm1,%%ymm1			\n\t	vaddpd	%%ymm9 ,%%ymm9 ,%%ymm9 				\n\t"\
		"vaddpd	%%ymm6,%%ymm1,%%ymm1			\n\t	vaddpd	%%ymm14,%%ymm9 ,%%ymm9 				\n\t"\
		"vsubpd	%%ymm2,%%ymm5,%%ymm5			\n\t	vsubpd	%%ymm10,%%ymm13,%%ymm13				\n\t"\
		"vaddpd	%%ymm2,%%ymm2,%%ymm2			\n\t	vaddpd	%%ymm10,%%ymm10,%%ymm10				\n\t"\
		"vaddpd	%%ymm5,%%ymm2,%%ymm2			\n\t	vaddpd	%%ymm13,%%ymm10,%%ymm10				\n\t"\
		"vmovaps	     (%%rdi),%%ymm0			\n\t	vmovaps	0x020(%%rdi),%%ymm8 				\n\t"\
		"vsubpd	%%ymm3,%%ymm4,%%ymm4			\n\t	vsubpd	%%ymm11,%%ymm12,%%ymm12				\n\t"\
		"vaddpd	%%ymm3,%%ymm3,%%ymm3			\n\t	vaddpd	%%ymm11,%%ymm11,%%ymm11				\n\t"\
		"vaddpd	%%ymm4,%%ymm3,%%ymm3			\n\t	vaddpd	%%ymm12,%%ymm11,%%ymm11				\n\t"\
		"/* Only alloc 2 scratch storage slots in calling code, so 16-reg version use 2 of the eventual outputs for Im parts */\n\t"\
		"vmovaps	%%ymm0,     (%%rsi)			\n\t	vmovaps	%%ymm8 ,0x100(%%rdi)				\n\t"\
		"vmovaps	%%ymm6,0x020(%%rsi)			\n\t	vmovaps	%%ymm14,0x200(%%rdi)				\n\t"\
		"vaddpd	      %%ymm1,%%ymm0,%%ymm0		\n\t	vaddpd	      %%ymm9 ,%%ymm8 ,%%ymm8 		\n\t"\
		"vmovaps	%%ymm5,%%ymm7				\n\t	vmovaps	%%ymm13,%%ymm15						\n\t"\
		"vaddpd	      %%ymm2,%%ymm3,%%ymm3		\n\t	vaddpd	      %%ymm10,%%ymm11,%%ymm11		\n\t"\
		"vsubpd	      %%ymm4,%%ymm5,%%ymm5		\n\t	vsubpd	      %%ymm12,%%ymm13,%%ymm13		\n\t"\
		"vsubpd	      %%ymm2,%%ymm1,%%ymm1		\n\t	vsubpd	      %%ymm10,%%ymm9 ,%%ymm9 		\n\t"\
		"vsubpd	      %%ymm7,%%ymm6,%%ymm6		\n\t	vsubpd	      %%ymm15,%%ymm14,%%ymm14		\n\t"\
		"vaddpd	      %%ymm2,%%ymm2,%%ymm2		\n\t	vaddpd	      %%ymm10,%%ymm10,%%ymm10		\n\t"\
		"vaddpd	      %%ymm7,%%ymm4,%%ymm4		\n\t	vaddpd	      %%ymm15,%%ymm12,%%ymm12		\n\t"\
		"vaddpd	      %%ymm3,%%ymm0,%%ymm0		\n\t	vaddpd	      %%ymm11,%%ymm8 ,%%ymm8 		\n\t"\
		"vaddpd	0x020(%%rsi),%%ymm5,%%ymm5		\n\t	vaddpd	0x200(%%rdi),%%ymm13,%%ymm13		\n\t"\
		"vsubpd	      %%ymm2,%%ymm3,%%ymm3		\n\t	vsubpd	      %%ymm10,%%ymm11,%%ymm11		\n\t"\
		"vmovaps	%%ymm4,%%ymm7				\n\t	vmovaps	%%ymm12,%%ymm15						\n\t"\
		"vmovaps	%%ymm0,      (%%rdi)		\n\t	vmovaps	%%ymm8 ,0x020(%%rdi)				\n\t"\
		"vsubpd	      %%ymm6,%%ymm4,%%ymm4		\n\t	vsubpd	      %%ymm14,%%ymm12,%%ymm12		\n\t"\
		"vmovaps	%%ymm1,%%ymm2				\n\t	vmovaps	%%ymm9 ,%%ymm10						\n\t"\
		"vsubpd	     (%%rsi),%%ymm0,%%ymm0		\n\t	vsubpd	0x100(%%rdi),%%ymm8 ,%%ymm8 		\n\t"\
		"vmulpd	0x020(%%rax),%%ymm5,%%ymm5		\n\t	vmulpd	0x020(%%rax),%%ymm13,%%ymm13		\n\t"\
		"vaddpd	      %%ymm3,%%ymm2,%%ymm2		\n\t	vaddpd	      %%ymm11,%%ymm10,%%ymm10		\n\t"\
		"vmulpd	     (%%rcx),%%ymm3,%%ymm3		\n\t	vmulpd	     (%%rcx),%%ymm11,%%ymm11		\n\t"\
		"vmulpd	0x020(%%rdx),%%ymm4,%%ymm4		\n\t	vmulpd	0x020(%%rdx),%%ymm12,%%ymm12		\n\t"\
		"vmulpd	     (%%rbx),%%ymm1,%%ymm1		\n\t	vmulpd	     (%%rbx),%%ymm9 ,%%ymm9 		\n\t"\
		"vmulpd	0x020(%%rbx),%%ymm6,%%ymm6		\n\t	vmulpd	0x020(%%rbx),%%ymm14,%%ymm14		\n\t"\
		"vmulpd	     (%%rax),%%ymm0,%%ymm0		\n\t	vmulpd	     (%%rax),%%ymm8 ,%%ymm8 		\n\t"\
		"vmulpd	0x020(%%rcx),%%ymm7,%%ymm7		\n\t	vmulpd	0x020(%%rcx),%%ymm15,%%ymm15		\n\t"\
		"vmulpd	     (%%rdx),%%ymm2,%%ymm2		\n\t	vmulpd	     (%%rdx),%%ymm10,%%ymm10		\n\t"\
		"vaddpd	     (%%rdi),%%ymm0,%%ymm0		\n\t	vaddpd	0x020(%%rdi),%%ymm8 ,%%ymm8 		\n\t"\
		"vaddpd	      %%ymm4,%%ymm6,%%ymm6		\n\t	vaddpd	      %%ymm12,%%ymm14,%%ymm14		\n\t"\
		"vsubpd	      %%ymm2,%%ymm1,%%ymm1		\n\t	vsubpd	      %%ymm10,%%ymm9 ,%%ymm9 		\n\t"\
		"vsubpd	      %%ymm7,%%ymm4,%%ymm4		\n\t	vsubpd	      %%ymm15,%%ymm12,%%ymm12		\n\t"\
		"vsubpd	      %%ymm2,%%ymm3,%%ymm3		\n\t	vsubpd	      %%ymm10,%%ymm11,%%ymm11		\n\t"\
		"vmovaps	%%ymm0,%%ymm2				\n\t	vmovaps	%%ymm8 ,%%ymm10						\n\t"\
		"vmovaps	%%ymm5,%%ymm7				\n\t	vmovaps	%%ymm13,%%ymm15						\n\t"\
		"vaddpd	      %%ymm1,%%ymm0,%%ymm0		\n\t	vaddpd	      %%ymm9 ,%%ymm8 ,%%ymm8 		\n\t"\
		"vaddpd	      %%ymm6,%%ymm5,%%ymm5		\n\t	vaddpd	      %%ymm14,%%ymm13,%%ymm13		\n\t"\
		"vaddpd	      %%ymm3,%%ymm1,%%ymm1		\n\t	vaddpd	      %%ymm11,%%ymm9 ,%%ymm9 		\n\t"\
		"vaddpd	      %%ymm4,%%ymm6,%%ymm6		\n\t	vaddpd	      %%ymm12,%%ymm14,%%ymm14		\n\t"\
		"vaddpd	      %%ymm2,%%ymm3,%%ymm3		\n\t	vaddpd	      %%ymm10,%%ymm11,%%ymm11		\n\t"\
		"vaddpd	      %%ymm7,%%ymm4,%%ymm4		\n\t	vaddpd	      %%ymm15,%%ymm12,%%ymm12		\n\t"\
		"vsubpd	      %%ymm1,%%ymm2,%%ymm2		\n\t	vsubpd	      %%ymm9 ,%%ymm10,%%ymm10		\n\t"\
		"vsubpd	      %%ymm6,%%ymm7,%%ymm7		\n\t	vsubpd	      %%ymm14,%%ymm15,%%ymm15		\n\t"\
		"/* ymm1,6,9,14 free ... Note the order reversal on the 3rd pair of outputs: */\n\t"\
		"vsubpd	%%ymm13,%%ymm0 ,%%ymm0		\n\t	vsubpd	%%ymm15,%%ymm2 ,%%ymm2	\n\t	vsubpd	%%ymm12,%%ymm3 ,%%ymm3 		\n\t"\
		"vsubpd	%%ymm5 ,%%ymm8 ,%%ymm8 		\n\t	vsubpd	%%ymm7 ,%%ymm10,%%ymm10	\n\t	vsubpd	%%ymm4 ,%%ymm11,%%ymm11		\n\t"\
		"vaddpd	%%ymm13,%%ymm13,%%ymm13		\n\t	vaddpd	%%ymm15,%%ymm15,%%ymm15	\n\t	vaddpd	%%ymm12,%%ymm12,%%ymm12		\n\t"\
		"vaddpd	%%ymm5 ,%%ymm5 ,%%ymm5		\n\t	vaddpd	%%ymm7 ,%%ymm7 ,%%ymm7	\n\t	vaddpd	%%ymm4 ,%%ymm4 ,%%ymm4 		\n\t"\
		"vaddpd	%%ymm0 ,%%ymm13,%%ymm13		\n\t	vaddpd	%%ymm2 ,%%ymm15,%%ymm15	\n\t	vaddpd	%%ymm3 ,%%ymm12,%%ymm12		\n\t"\
		"vaddpd	%%ymm8 ,%%ymm5 ,%%ymm5		\n\t	vaddpd	%%ymm10,%%ymm7 ,%%ymm7	\n\t	vaddpd	%%ymm11,%%ymm4 ,%%ymm4 		\n\t"\
		"vmovaps	%%ymm0 ,0x100(%%rdi)	\n\t	vmovaps	%%ymm2 ,0x200(%%rdi)	\n\t	vmovaps	%%ymm3 ,0x400(%%rdi)		\n\t"\
		"vmovaps	%%ymm8 ,0x620(%%rdi)	\n\t	vmovaps	%%ymm10,0x520(%%rdi)	\n\t	vmovaps	%%ymm11,0x320(%%rdi)		\n\t"\
		"vmovaps	%%ymm13,0x600(%%rdi)	\n\t	vmovaps	%%ymm15,0x500(%%rdi)	\n\t	vmovaps	%%ymm12,0x300(%%rdi)		\n\t"\
		"vmovaps	%%ymm5 ,0x120(%%rdi)	\n\t	vmovaps	%%ymm7 ,0x220(%%rdi)	\n\t	vmovaps	%%ymm4 ,0x420(%%rdi)		\n\t"\
	"/* SSE2_RADIX_07_DFT(21,17,13,09,05,01,25 -> 01,05,09,13,17,21,25): */\n\t"\
		"addq	$0x040,%%rdi	/* s1p01r */	\n\t	/*** Imaginary Parts: ***/					\n\t"\
		"vmovaps	0x400(%%rdi),%%ymm6			\n\t	vmovaps	0x420(%%rdi),%%ymm14				\n\t"\
		"vmovaps	0x600(%%rdi),%%ymm1			\n\t	vmovaps	0x620(%%rdi),%%ymm9 				\n\t"\
		"vmovaps	0x300(%%rdi),%%ymm5			\n\t	vmovaps	0x320(%%rdi),%%ymm13				\n\t"\
		"vmovaps	     (%%rdi),%%ymm2			\n\t	vmovaps	0x020(%%rdi),%%ymm10				\n\t"\
		"vmovaps	0x200(%%rdi),%%ymm4			\n\t	vmovaps	0x220(%%rdi),%%ymm12				\n\t"\
		"vmovaps	0x100(%%rdi),%%ymm3			\n\t	vmovaps	0x120(%%rdi),%%ymm11				\n\t"\
		"vsubpd	%%ymm1,%%ymm6,%%ymm6			\n\t	vsubpd	%%ymm9 ,%%ymm14,%%ymm14				\n\t"\
		"vaddpd	%%ymm1,%%ymm1,%%ymm1			\n\t	vaddpd	%%ymm9 ,%%ymm9 ,%%ymm9 				\n\t"\
		"vaddpd	%%ymm6,%%ymm1,%%ymm1			\n\t	vaddpd	%%ymm14,%%ymm9 ,%%ymm9 				\n\t"\
		"vsubpd	%%ymm2,%%ymm5,%%ymm5			\n\t	vsubpd	%%ymm10,%%ymm13,%%ymm13				\n\t"\
		"vaddpd	%%ymm2,%%ymm2,%%ymm2			\n\t	vaddpd	%%ymm10,%%ymm10,%%ymm10				\n\t"\
		"vaddpd	%%ymm5,%%ymm2,%%ymm2			\n\t	vaddpd	%%ymm13,%%ymm10,%%ymm10				\n\t"\
		"vmovaps	0x500(%%rdi),%%ymm0			\n\t	vmovaps	0x520(%%rdi),%%ymm8 				\n\t"\
		"vsubpd	%%ymm3,%%ymm4,%%ymm4			\n\t	vsubpd	%%ymm11,%%ymm12,%%ymm12				\n\t"\
		"vaddpd	%%ymm3,%%ymm3,%%ymm3			\n\t	vaddpd	%%ymm11,%%ymm11,%%ymm11				\n\t"\
		"vaddpd	%%ymm4,%%ymm3,%%ymm3			\n\t	vaddpd	%%ymm12,%%ymm11,%%ymm11				\n\t"\
		"/* Only alloc 2 scratch storage slots in calling code, so 16-reg version use 2 of the eventual outputs for Im parts */\n\t"\
		"vmovaps	%%ymm0,     (%%rsi)			\n\t	vmovaps	%%ymm8 ,0x100(%%rdi)				\n\t"\
		"vmovaps	%%ymm6,0x020(%%rsi)			\n\t	vmovaps	%%ymm14,0x200(%%rdi)				\n\t"\
		"vaddpd	      %%ymm1,%%ymm0,%%ymm0		\n\t	vaddpd	      %%ymm9 ,%%ymm8 ,%%ymm8 		\n\t"\
		"vmovaps	%%ymm5,%%ymm7				\n\t	vmovaps	%%ymm13,%%ymm15						\n\t"\
		"vaddpd	      %%ymm2,%%ymm3,%%ymm3		\n\t	vaddpd	      %%ymm10,%%ymm11,%%ymm11		\n\t"\
		"vsubpd	      %%ymm4,%%ymm5,%%ymm5		\n\t	vsubpd	      %%ymm12,%%ymm13,%%ymm13		\n\t"\
		"vsubpd	      %%ymm2,%%ymm1,%%ymm1		\n\t	vsubpd	      %%ymm10,%%ymm9 ,%%ymm9 		\n\t"\
		"vsubpd	      %%ymm7,%%ymm6,%%ymm6		\n\t	vsubpd	      %%ymm15,%%ymm14,%%ymm14		\n\t"\
		"vaddpd	      %%ymm2,%%ymm2,%%ymm2		\n\t	vaddpd	      %%ymm10,%%ymm10,%%ymm10		\n\t"\
		"vaddpd	      %%ymm7,%%ymm4,%%ymm4		\n\t	vaddpd	      %%ymm15,%%ymm12,%%ymm12		\n\t"\
		"vaddpd	      %%ymm3,%%ymm0,%%ymm0		\n\t	vaddpd	      %%ymm11,%%ymm8 ,%%ymm8 		\n\t"\
		"vaddpd	0x020(%%rsi),%%ymm5,%%ymm5		\n\t	vaddpd	0x200(%%rdi),%%ymm13,%%ymm13		\n\t"\
		"vsubpd	      %%ymm2,%%ymm3,%%ymm3		\n\t	vsubpd	      %%ymm10,%%ymm11,%%ymm11		\n\t"\
		"vmovaps	%%ymm4,%%ymm7				\n\t	vmovaps	%%ymm12,%%ymm15						\n\t"\
		"vmovaps	%%ymm0,      (%%rdi)		\n\t	vmovaps	%%ymm8 ,0x020(%%rdi)				\n\t"\
		"vsubpd	      %%ymm6,%%ymm4,%%ymm4		\n\t	vsubpd	      %%ymm14,%%ymm12,%%ymm12		\n\t"\
		"vmovaps	%%ymm1,%%ymm2				\n\t	vmovaps	%%ymm9 ,%%ymm10						\n\t"\
		"vsubpd	     (%%rsi),%%ymm0,%%ymm0		\n\t	vsubpd	0x100(%%rdi),%%ymm8 ,%%ymm8 		\n\t"\
		"vmulpd	0x020(%%rax),%%ymm5,%%ymm5		\n\t	vmulpd	0x020(%%rax),%%ymm13,%%ymm13		\n\t"\
		"vaddpd	      %%ymm3,%%ymm2,%%ymm2		\n\t	vaddpd	      %%ymm11,%%ymm10,%%ymm10		\n\t"\
		"vmulpd	     (%%rcx),%%ymm3,%%ymm3		\n\t	vmulpd	     (%%rcx),%%ymm11,%%ymm11		\n\t"\
		"vmulpd	0x020(%%rdx),%%ymm4,%%ymm4		\n\t	vmulpd	0x020(%%rdx),%%ymm12,%%ymm12		\n\t"\
		"vmulpd	     (%%rbx),%%ymm1,%%ymm1		\n\t	vmulpd	     (%%rbx),%%ymm9 ,%%ymm9 		\n\t"\
		"vmulpd	0x020(%%rbx),%%ymm6,%%ymm6		\n\t	vmulpd	0x020(%%rbx),%%ymm14,%%ymm14		\n\t"\
		"vmulpd	     (%%rax),%%ymm0,%%ymm0		\n\t	vmulpd	     (%%rax),%%ymm8 ,%%ymm8 		\n\t"\
		"vmulpd	0x020(%%rcx),%%ymm7,%%ymm7		\n\t	vmulpd	0x020(%%rcx),%%ymm15,%%ymm15		\n\t"\
		"vmulpd	     (%%rdx),%%ymm2,%%ymm2		\n\t	vmulpd	     (%%rdx),%%ymm10,%%ymm10		\n\t"\
		"vaddpd	     (%%rdi),%%ymm0,%%ymm0		\n\t	vaddpd	0x020(%%rdi),%%ymm8 ,%%ymm8 		\n\t"\
		"vaddpd	      %%ymm4,%%ymm6,%%ymm6		\n\t	vaddpd	      %%ymm12,%%ymm14,%%ymm14		\n\t"\
		"vsubpd	      %%ymm2,%%ymm1,%%ymm1		\n\t	vsubpd	      %%ymm10,%%ymm9 ,%%ymm9 		\n\t"\
		"vsubpd	      %%ymm7,%%ymm4,%%ymm4		\n\t	vsubpd	      %%ymm15,%%ymm12,%%ymm12		\n\t"\
		"vsubpd	      %%ymm2,%%ymm3,%%ymm3		\n\t	vsubpd	      %%ymm10,%%ymm11,%%ymm11		\n\t"\
		"vmovaps	%%ymm0,%%ymm2				\n\t	vmovaps	%%ymm8 ,%%ymm10						\n\t"\
		"vmovaps	%%ymm5,%%ymm7				\n\t	vmovaps	%%ymm13,%%ymm15						\n\t"\
		"vaddpd	      %%ymm1,%%ymm0,%%ymm0		\n\t	vaddpd	      %%ymm9 ,%%ymm8 ,%%ymm8 		\n\t"\
		"vaddpd	      %%ymm6,%%ymm5,%%ymm5		\n\t	vaddpd	      %%ymm14,%%ymm13,%%ymm13		\n\t"\
		"vaddpd	      %%ymm3,%%ymm1,%%ymm1		\n\t	vaddpd	      %%ymm11,%%ymm9 ,%%ymm9 		\n\t"\
		"vaddpd	      %%ymm4,%%ymm6,%%ymm6		\n\t	vaddpd	      %%ymm12,%%ymm14,%%ymm14		\n\t"\
		"vaddpd	      %%ymm2,%%ymm3,%%ymm3		\n\t	vaddpd	      %%ymm10,%%ymm11,%%ymm11		\n\t"\
		"vaddpd	      %%ymm7,%%ymm4,%%ymm4		\n\t	vaddpd	      %%ymm15,%%ymm12,%%ymm12		\n\t"\
		"vsubpd	      %%ymm1,%%ymm2,%%ymm2		\n\t	vsubpd	      %%ymm9 ,%%ymm10,%%ymm10		\n\t"\
		"vsubpd	      %%ymm6,%%ymm7,%%ymm7		\n\t	vsubpd	      %%ymm14,%%ymm15,%%ymm15		\n\t"\
		"/* ymm1,6,9,14 free ... Note the order reversal on the 3rd pair of outputs: */\n\t"\
		"vsubpd	%%ymm13,%%ymm0 ,%%ymm0		\n\t	vsubpd	%%ymm15,%%ymm2 ,%%ymm2	\n\t	vsubpd	%%ymm12,%%ymm3 ,%%ymm3 		\n\t"\
		"vsubpd	%%ymm5 ,%%ymm8 ,%%ymm8 		\n\t	vsubpd	%%ymm7 ,%%ymm10,%%ymm10	\n\t	vsubpd	%%ymm4 ,%%ymm11,%%ymm11		\n\t"\
		"vaddpd	%%ymm13,%%ymm13,%%ymm13		\n\t	vaddpd	%%ymm15,%%ymm15,%%ymm15	\n\t	vaddpd	%%ymm12,%%ymm12,%%ymm12		\n\t"\
		"vaddpd	%%ymm5 ,%%ymm5 ,%%ymm5		\n\t	vaddpd	%%ymm7 ,%%ymm7 ,%%ymm7	\n\t	vaddpd	%%ymm4 ,%%ymm4 ,%%ymm4 		\n\t"\
		"vaddpd	%%ymm0 ,%%ymm13,%%ymm13		\n\t	vaddpd	%%ymm2 ,%%ymm15,%%ymm15	\n\t	vaddpd	%%ymm3 ,%%ymm12,%%ymm12		\n\t"\
		"vaddpd	%%ymm8 ,%%ymm5 ,%%ymm5		\n\t	vaddpd	%%ymm10,%%ymm7 ,%%ymm7	\n\t	vaddpd	%%ymm11,%%ymm4 ,%%ymm4 		\n\t"\
		"vmovaps	%%ymm0 ,0x100(%%rdi)	\n\t	vmovaps	%%ymm2 ,0x200(%%rdi)	\n\t	vmovaps	%%ymm3 ,0x400(%%rdi)		\n\t"\
		"vmovaps	%%ymm8 ,0x620(%%rdi)	\n\t	vmovaps	%%ymm10,0x520(%%rdi)	\n\t	vmovaps	%%ymm11,0x320(%%rdi)		\n\t"\
		"vmovaps	%%ymm13,0x600(%%rdi)	\n\t	vmovaps	%%ymm15,0x500(%%rdi)	\n\t	vmovaps	%%ymm12,0x300(%%rdi)		\n\t"\
		"vmovaps	%%ymm5 ,0x120(%%rdi)	\n\t	vmovaps	%%ymm7 ,0x220(%%rdi)	\n\t	vmovaps	%%ymm4 ,0x420(%%rdi)		\n\t"\
	"/* SSE2_RADIX_07_DFT(14,10,06,02,26,22,18 -> 02,06,10,14,18,22,26): */\n\t"\
		"addq	$0x040,%%rdi	/* s1p02r */	\n\t	/*** Imaginary Parts: ***/					\n\t"\
		"vmovaps	0x200(%%rdi),%%ymm6			\n\t	vmovaps	0x220(%%rdi),%%ymm14				\n\t"\
		"vmovaps	0x400(%%rdi),%%ymm1			\n\t	vmovaps	0x420(%%rdi),%%ymm9 				\n\t"\
		"vmovaps	0x100(%%rdi),%%ymm5			\n\t	vmovaps	0x120(%%rdi),%%ymm13				\n\t"\
		"vmovaps	0x500(%%rdi),%%ymm2			\n\t	vmovaps	0x520(%%rdi),%%ymm10				\n\t"\
		"vmovaps	     (%%rdi),%%ymm4			\n\t	vmovaps	0x020(%%rdi),%%ymm12				\n\t"\
		"vmovaps	0x600(%%rdi),%%ymm3			\n\t	vmovaps	0x620(%%rdi),%%ymm11				\n\t"\
		"vsubpd	%%ymm1,%%ymm6,%%ymm6			\n\t	vsubpd	%%ymm9 ,%%ymm14,%%ymm14				\n\t"\
		"vaddpd	%%ymm1,%%ymm1,%%ymm1			\n\t	vaddpd	%%ymm9 ,%%ymm9 ,%%ymm9 				\n\t"\
		"vaddpd	%%ymm6,%%ymm1,%%ymm1			\n\t	vaddpd	%%ymm14,%%ymm9 ,%%ymm9 				\n\t"\
		"vsubpd	%%ymm2,%%ymm5,%%ymm5			\n\t	vsubpd	%%ymm10,%%ymm13,%%ymm13				\n\t"\
		"vaddpd	%%ymm2,%%ymm2,%%ymm2			\n\t	vaddpd	%%ymm10,%%ymm10,%%ymm10				\n\t"\
		"vaddpd	%%ymm5,%%ymm2,%%ymm2			\n\t	vaddpd	%%ymm13,%%ymm10,%%ymm10				\n\t"\
		"vmovaps	0x300(%%rdi),%%ymm0			\n\t	vmovaps	0x320(%%rdi),%%ymm8 				\n\t"\
		"vsubpd	%%ymm3,%%ymm4,%%ymm4			\n\t	vsubpd	%%ymm11,%%ymm12,%%ymm12				\n\t"\
		"vaddpd	%%ymm3,%%ymm3,%%ymm3			\n\t	vaddpd	%%ymm11,%%ymm11,%%ymm11				\n\t"\
		"vaddpd	%%ymm4,%%ymm3,%%ymm3			\n\t	vaddpd	%%ymm12,%%ymm11,%%ymm11				\n\t"\
		"/* Only alloc 2 scratch storage slots in calling code, so 16-reg version use 2 of the eventual outputs for Im parts */\n\t"\
		"vmovaps	%%ymm0,     (%%rsi)			\n\t	vmovaps	%%ymm8 ,0x100(%%rdi)				\n\t"\
		"vmovaps	%%ymm6,0x020(%%rsi)			\n\t	vmovaps	%%ymm14,0x200(%%rdi)				\n\t"\
		"vaddpd	      %%ymm1,%%ymm0,%%ymm0		\n\t	vaddpd	      %%ymm9 ,%%ymm8 ,%%ymm8 		\n\t"\
		"vmovaps	%%ymm5,%%ymm7				\n\t	vmovaps	%%ymm13,%%ymm15						\n\t"\
		"vaddpd	      %%ymm2,%%ymm3,%%ymm3		\n\t	vaddpd	      %%ymm10,%%ymm11,%%ymm11		\n\t"\
		"vsubpd	      %%ymm4,%%ymm5,%%ymm5		\n\t	vsubpd	      %%ymm12,%%ymm13,%%ymm13		\n\t"\
		"vsubpd	      %%ymm2,%%ymm1,%%ymm1		\n\t	vsubpd	      %%ymm10,%%ymm9 ,%%ymm9 		\n\t"\
		"vsubpd	      %%ymm7,%%ymm6,%%ymm6		\n\t	vsubpd	      %%ymm15,%%ymm14,%%ymm14		\n\t"\
		"vaddpd	      %%ymm2,%%ymm2,%%ymm2		\n\t	vaddpd	      %%ymm10,%%ymm10,%%ymm10		\n\t"\
		"vaddpd	      %%ymm7,%%ymm4,%%ymm4		\n\t	vaddpd	      %%ymm15,%%ymm12,%%ymm12		\n\t"\
		"vaddpd	      %%ymm3,%%ymm0,%%ymm0		\n\t	vaddpd	      %%ymm11,%%ymm8 ,%%ymm8 		\n\t"\
		"vaddpd	0x020(%%rsi),%%ymm5,%%ymm5		\n\t	vaddpd	0x200(%%rdi),%%ymm13,%%ymm13		\n\t"\
		"vsubpd	      %%ymm2,%%ymm3,%%ymm3		\n\t	vsubpd	      %%ymm10,%%ymm11,%%ymm11		\n\t"\
		"vmovaps	%%ymm4,%%ymm7				\n\t	vmovaps	%%ymm12,%%ymm15						\n\t"\
		"vmovaps	%%ymm0,      (%%rdi)		\n\t	vmovaps	%%ymm8 ,0x020(%%rdi)				\n\t"\
		"vsubpd	      %%ymm6,%%ymm4,%%ymm4		\n\t	vsubpd	      %%ymm14,%%ymm12,%%ymm12		\n\t"\
		"vmovaps	%%ymm1,%%ymm2				\n\t	vmovaps	%%ymm9 ,%%ymm10						\n\t"\
		"vsubpd	     (%%rsi),%%ymm0,%%ymm0		\n\t	vsubpd	0x100(%%rdi),%%ymm8 ,%%ymm8 		\n\t"\
		"vmulpd	0x020(%%rax),%%ymm5,%%ymm5		\n\t	vmulpd	0x020(%%rax),%%ymm13,%%ymm13		\n\t"\
		"vaddpd	      %%ymm3,%%ymm2,%%ymm2		\n\t	vaddpd	      %%ymm11,%%ymm10,%%ymm10		\n\t"\
		"vmulpd	     (%%rcx),%%ymm3,%%ymm3		\n\t	vmulpd	     (%%rcx),%%ymm11,%%ymm11		\n\t"\
		"vmulpd	0x020(%%rdx),%%ymm4,%%ymm4		\n\t	vmulpd	0x020(%%rdx),%%ymm12,%%ymm12		\n\t"\
		"vmulpd	     (%%rbx),%%ymm1,%%ymm1		\n\t	vmulpd	     (%%rbx),%%ymm9 ,%%ymm9 		\n\t"\
		"vmulpd	0x020(%%rbx),%%ymm6,%%ymm6		\n\t	vmulpd	0x020(%%rbx),%%ymm14,%%ymm14		\n\t"\
		"vmulpd	     (%%rax),%%ymm0,%%ymm0		\n\t	vmulpd	     (%%rax),%%ymm8 ,%%ymm8 		\n\t"\
		"vmulpd	0x020(%%rcx),%%ymm7,%%ymm7		\n\t	vmulpd	0x020(%%rcx),%%ymm15,%%ymm15		\n\t"\
		"vmulpd	     (%%rdx),%%ymm2,%%ymm2		\n\t	vmulpd	     (%%rdx),%%ymm10,%%ymm10		\n\t"\
		"vaddpd	     (%%rdi),%%ymm0,%%ymm0		\n\t	vaddpd	0x020(%%rdi),%%ymm8 ,%%ymm8 		\n\t"\
		"vaddpd	      %%ymm4,%%ymm6,%%ymm6		\n\t	vaddpd	      %%ymm12,%%ymm14,%%ymm14		\n\t"\
		"vsubpd	      %%ymm2,%%ymm1,%%ymm1		\n\t	vsubpd	      %%ymm10,%%ymm9 ,%%ymm9 		\n\t"\
		"vsubpd	      %%ymm7,%%ymm4,%%ymm4		\n\t	vsubpd	      %%ymm15,%%ymm12,%%ymm12		\n\t"\
		"vsubpd	      %%ymm2,%%ymm3,%%ymm3		\n\t	vsubpd	      %%ymm10,%%ymm11,%%ymm11		\n\t"\
		"vmovaps	%%ymm0,%%ymm2				\n\t	vmovaps	%%ymm8 ,%%ymm10						\n\t"\
		"vmovaps	%%ymm5,%%ymm7				\n\t	vmovaps	%%ymm13,%%ymm15						\n\t"\
		"vaddpd	      %%ymm1,%%ymm0,%%ymm0		\n\t	vaddpd	      %%ymm9 ,%%ymm8 ,%%ymm8 		\n\t"\
		"vaddpd	      %%ymm6,%%ymm5,%%ymm5		\n\t	vaddpd	      %%ymm14,%%ymm13,%%ymm13		\n\t"\
		"vaddpd	      %%ymm3,%%ymm1,%%ymm1		\n\t	vaddpd	      %%ymm11,%%ymm9 ,%%ymm9 		\n\t"\
		"vaddpd	      %%ymm4,%%ymm6,%%ymm6		\n\t	vaddpd	      %%ymm12,%%ymm14,%%ymm14		\n\t"\
		"vaddpd	      %%ymm2,%%ymm3,%%ymm3		\n\t	vaddpd	      %%ymm10,%%ymm11,%%ymm11		\n\t"\
		"vaddpd	      %%ymm7,%%ymm4,%%ymm4		\n\t	vaddpd	      %%ymm15,%%ymm12,%%ymm12		\n\t"\
		"vsubpd	      %%ymm1,%%ymm2,%%ymm2		\n\t	vsubpd	      %%ymm9 ,%%ymm10,%%ymm10		\n\t"\
		"vsubpd	      %%ymm6,%%ymm7,%%ymm7		\n\t	vsubpd	      %%ymm14,%%ymm15,%%ymm15		\n\t"\
		"/* ymm1,6,9,14 free ... Note the order reversal on the 3rd pair of outputs: */\n\t"\
		"vsubpd	%%ymm13,%%ymm0 ,%%ymm0		\n\t	vsubpd	%%ymm15,%%ymm2 ,%%ymm2	\n\t	vsubpd	%%ymm12,%%ymm3 ,%%ymm3 		\n\t"\
		"vsubpd	%%ymm5 ,%%ymm8 ,%%ymm8 		\n\t	vsubpd	%%ymm7 ,%%ymm10,%%ymm10	\n\t	vsubpd	%%ymm4 ,%%ymm11,%%ymm11		\n\t"\
		"vaddpd	%%ymm13,%%ymm13,%%ymm13		\n\t	vaddpd	%%ymm15,%%ymm15,%%ymm15	\n\t	vaddpd	%%ymm12,%%ymm12,%%ymm12		\n\t"\
		"vaddpd	%%ymm5 ,%%ymm5 ,%%ymm5		\n\t	vaddpd	%%ymm7 ,%%ymm7 ,%%ymm7	\n\t	vaddpd	%%ymm4 ,%%ymm4 ,%%ymm4 		\n\t"\
		"vaddpd	%%ymm0 ,%%ymm13,%%ymm13		\n\t	vaddpd	%%ymm2 ,%%ymm15,%%ymm15	\n\t	vaddpd	%%ymm3 ,%%ymm12,%%ymm12		\n\t"\
		"vaddpd	%%ymm8 ,%%ymm5 ,%%ymm5		\n\t	vaddpd	%%ymm10,%%ymm7 ,%%ymm7	\n\t	vaddpd	%%ymm11,%%ymm4 ,%%ymm4 		\n\t"\
		"vmovaps	%%ymm0 ,0x100(%%rdi)	\n\t	vmovaps	%%ymm2 ,0x200(%%rdi)	\n\t	vmovaps	%%ymm3 ,0x400(%%rdi)		\n\t"\
		"vmovaps	%%ymm8 ,0x620(%%rdi)	\n\t	vmovaps	%%ymm10,0x520(%%rdi)	\n\t	vmovaps	%%ymm11,0x320(%%rdi)		\n\t"\
		"vmovaps	%%ymm13,0x600(%%rdi)	\n\t	vmovaps	%%ymm15,0x500(%%rdi)	\n\t	vmovaps	%%ymm12,0x300(%%rdi)		\n\t"\
		"vmovaps	%%ymm5 ,0x120(%%rdi)	\n\t	vmovaps	%%ymm7 ,0x220(%%rdi)	\n\t	vmovaps	%%ymm4 ,0x420(%%rdi)		\n\t"\
	"/* SSE2_RADIX_07_DFT(07,03,27,23,19,15,11 -> 03,07,11,15,19,23,27): */\n\t"\
		"addq	$0x040,%%rdi	/* s1p03r */	\n\t	/*** Imaginary Parts: ***/					\n\t"\
		"vmovaps	     (%%rdi),%%ymm6			\n\t	vmovaps	0x020(%%rdi),%%ymm14				\n\t"\
		"vmovaps	0x200(%%rdi),%%ymm1			\n\t	vmovaps	0x220(%%rdi),%%ymm9 				\n\t"\
		"vmovaps	0x600(%%rdi),%%ymm5			\n\t	vmovaps	0x620(%%rdi),%%ymm13				\n\t"\
		"vmovaps	0x300(%%rdi),%%ymm2			\n\t	vmovaps	0x320(%%rdi),%%ymm10				\n\t"\
		"vmovaps	0x500(%%rdi),%%ymm4			\n\t	vmovaps	0x520(%%rdi),%%ymm12				\n\t"\
		"vmovaps	0x400(%%rdi),%%ymm3			\n\t	vmovaps	0x420(%%rdi),%%ymm11				\n\t"\
		"vsubpd	%%ymm1,%%ymm6,%%ymm6			\n\t	vsubpd	%%ymm9 ,%%ymm14,%%ymm14				\n\t"\
		"vaddpd	%%ymm1,%%ymm1,%%ymm1			\n\t	vaddpd	%%ymm9 ,%%ymm9 ,%%ymm9 				\n\t"\
		"vaddpd	%%ymm6,%%ymm1,%%ymm1			\n\t	vaddpd	%%ymm14,%%ymm9 ,%%ymm9 				\n\t"\
		"vsubpd	%%ymm2,%%ymm5,%%ymm5			\n\t	vsubpd	%%ymm10,%%ymm13,%%ymm13				\n\t"\
		"vaddpd	%%ymm2,%%ymm2,%%ymm2			\n\t	vaddpd	%%ymm10,%%ymm10,%%ymm10				\n\t"\
		"vaddpd	%%ymm5,%%ymm2,%%ymm2			\n\t	vaddpd	%%ymm13,%%ymm10,%%ymm10				\n\t"\
		"vmovaps	0x100(%%rdi),%%ymm0			\n\t	vmovaps	0x120(%%rdi),%%ymm8 				\n\t"\
		"vsubpd	%%ymm3,%%ymm4,%%ymm4			\n\t	vsubpd	%%ymm11,%%ymm12,%%ymm12				\n\t"\
		"vaddpd	%%ymm3,%%ymm3,%%ymm3			\n\t	vaddpd	%%ymm11,%%ymm11,%%ymm11				\n\t"\
		"vaddpd	%%ymm4,%%ymm3,%%ymm3			\n\t	vaddpd	%%ymm12,%%ymm11,%%ymm11				\n\t"\
		"/* Only alloc 2 scratch storage slots in calling code, so 16-reg version use 2 of the eventual outputs for Im parts */\n\t"\
		"vmovaps	%%ymm0,     (%%rsi)			\n\t	vmovaps	%%ymm8 ,0x100(%%rdi)				\n\t"\
		"vmovaps	%%ymm6,0x020(%%rsi)			\n\t	vmovaps	%%ymm14,0x200(%%rdi)				\n\t"\
		"vaddpd	      %%ymm1,%%ymm0,%%ymm0		\n\t	vaddpd	      %%ymm9 ,%%ymm8 ,%%ymm8 		\n\t"\
		"vmovaps	%%ymm5,%%ymm7				\n\t	vmovaps	%%ymm13,%%ymm15						\n\t"\
		"vaddpd	      %%ymm2,%%ymm3,%%ymm3		\n\t	vaddpd	      %%ymm10,%%ymm11,%%ymm11		\n\t"\
		"vsubpd	      %%ymm4,%%ymm5,%%ymm5		\n\t	vsubpd	      %%ymm12,%%ymm13,%%ymm13		\n\t"\
		"vsubpd	      %%ymm2,%%ymm1,%%ymm1		\n\t	vsubpd	      %%ymm10,%%ymm9 ,%%ymm9 		\n\t"\
		"vsubpd	      %%ymm7,%%ymm6,%%ymm6		\n\t	vsubpd	      %%ymm15,%%ymm14,%%ymm14		\n\t"\
		"vaddpd	      %%ymm2,%%ymm2,%%ymm2		\n\t	vaddpd	      %%ymm10,%%ymm10,%%ymm10		\n\t"\
		"vaddpd	      %%ymm7,%%ymm4,%%ymm4		\n\t	vaddpd	      %%ymm15,%%ymm12,%%ymm12		\n\t"\
		"vaddpd	      %%ymm3,%%ymm0,%%ymm0		\n\t	vaddpd	      %%ymm11,%%ymm8 ,%%ymm8 		\n\t"\
		"vaddpd	0x020(%%rsi),%%ymm5,%%ymm5		\n\t	vaddpd	0x200(%%rdi),%%ymm13,%%ymm13		\n\t"\
		"vsubpd	      %%ymm2,%%ymm3,%%ymm3		\n\t	vsubpd	      %%ymm10,%%ymm11,%%ymm11		\n\t"\
		"vmovaps	%%ymm4,%%ymm7				\n\t	vmovaps	%%ymm12,%%ymm15						\n\t"\
		"vmovaps	%%ymm0,      (%%rdi)		\n\t	vmovaps	%%ymm8 ,0x020(%%rdi)				\n\t"\
		"vsubpd	      %%ymm6,%%ymm4,%%ymm4		\n\t	vsubpd	      %%ymm14,%%ymm12,%%ymm12		\n\t"\
		"vmovaps	%%ymm1,%%ymm2				\n\t	vmovaps	%%ymm9 ,%%ymm10						\n\t"\
		"vsubpd	     (%%rsi),%%ymm0,%%ymm0		\n\t	vsubpd	0x100(%%rdi),%%ymm8 ,%%ymm8 		\n\t"\
		"vmulpd	0x020(%%rax),%%ymm5,%%ymm5		\n\t	vmulpd	0x020(%%rax),%%ymm13,%%ymm13		\n\t"\
		"vaddpd	      %%ymm3,%%ymm2,%%ymm2		\n\t	vaddpd	      %%ymm11,%%ymm10,%%ymm10		\n\t"\
		"vmulpd	     (%%rcx),%%ymm3,%%ymm3		\n\t	vmulpd	     (%%rcx),%%ymm11,%%ymm11		\n\t"\
		"vmulpd	0x020(%%rdx),%%ymm4,%%ymm4		\n\t	vmulpd	0x020(%%rdx),%%ymm12,%%ymm12		\n\t"\
		"vmulpd	     (%%rbx),%%ymm1,%%ymm1		\n\t	vmulpd	     (%%rbx),%%ymm9 ,%%ymm9 		\n\t"\
		"vmulpd	0x020(%%rbx),%%ymm6,%%ymm6		\n\t	vmulpd	0x020(%%rbx),%%ymm14,%%ymm14		\n\t"\
		"vmulpd	     (%%rax),%%ymm0,%%ymm0		\n\t	vmulpd	     (%%rax),%%ymm8 ,%%ymm8 		\n\t"\
		"vmulpd	0x020(%%rcx),%%ymm7,%%ymm7		\n\t	vmulpd	0x020(%%rcx),%%ymm15,%%ymm15		\n\t"\
		"vmulpd	     (%%rdx),%%ymm2,%%ymm2		\n\t	vmulpd	     (%%rdx),%%ymm10,%%ymm10		\n\t"\
		"vaddpd	     (%%rdi),%%ymm0,%%ymm0		\n\t	vaddpd	0x020(%%rdi),%%ymm8 ,%%ymm8 		\n\t"\
		"vaddpd	      %%ymm4,%%ymm6,%%ymm6		\n\t	vaddpd	      %%ymm12,%%ymm14,%%ymm14		\n\t"\
		"vsubpd	      %%ymm2,%%ymm1,%%ymm1		\n\t	vsubpd	      %%ymm10,%%ymm9 ,%%ymm9 		\n\t"\
		"vsubpd	      %%ymm7,%%ymm4,%%ymm4		\n\t	vsubpd	      %%ymm15,%%ymm12,%%ymm12		\n\t"\
		"vsubpd	      %%ymm2,%%ymm3,%%ymm3		\n\t	vsubpd	      %%ymm10,%%ymm11,%%ymm11		\n\t"\
		"vmovaps	%%ymm0,%%ymm2				\n\t	vmovaps	%%ymm8 ,%%ymm10						\n\t"\
		"vmovaps	%%ymm5,%%ymm7				\n\t	vmovaps	%%ymm13,%%ymm15						\n\t"\
		"vaddpd	      %%ymm1,%%ymm0,%%ymm0		\n\t	vaddpd	      %%ymm9 ,%%ymm8 ,%%ymm8 		\n\t"\
		"vaddpd	      %%ymm6,%%ymm5,%%ymm5		\n\t	vaddpd	      %%ymm14,%%ymm13,%%ymm13		\n\t"\
		"vaddpd	      %%ymm3,%%ymm1,%%ymm1		\n\t	vaddpd	      %%ymm11,%%ymm9 ,%%ymm9 		\n\t"\
		"vaddpd	      %%ymm4,%%ymm6,%%ymm6		\n\t	vaddpd	      %%ymm12,%%ymm14,%%ymm14		\n\t"\
		"vaddpd	      %%ymm2,%%ymm3,%%ymm3		\n\t	vaddpd	      %%ymm10,%%ymm11,%%ymm11		\n\t"\
		"vaddpd	      %%ymm7,%%ymm4,%%ymm4		\n\t	vaddpd	      %%ymm15,%%ymm12,%%ymm12		\n\t"\
		"vsubpd	      %%ymm1,%%ymm2,%%ymm2		\n\t	vsubpd	      %%ymm9 ,%%ymm10,%%ymm10		\n\t"\
		"vsubpd	      %%ymm6,%%ymm7,%%ymm7		\n\t	vsubpd	      %%ymm14,%%ymm15,%%ymm15		\n\t"\
		"/* ymm1,6,9,14 free ... Note the order reversal on the 3rd pair of outputs: */\n\t"\
		"vsubpd	%%ymm13,%%ymm0 ,%%ymm0		\n\t	vsubpd	%%ymm15,%%ymm2 ,%%ymm2	\n\t	vsubpd	%%ymm12,%%ymm3 ,%%ymm3 		\n\t"\
		"vsubpd	%%ymm5 ,%%ymm8 ,%%ymm8 		\n\t	vsubpd	%%ymm7 ,%%ymm10,%%ymm10	\n\t	vsubpd	%%ymm4 ,%%ymm11,%%ymm11		\n\t"\
		"vaddpd	%%ymm13,%%ymm13,%%ymm13		\n\t	vaddpd	%%ymm15,%%ymm15,%%ymm15	\n\t	vaddpd	%%ymm12,%%ymm12,%%ymm12		\n\t"\
		"vaddpd	%%ymm5 ,%%ymm5 ,%%ymm5		\n\t	vaddpd	%%ymm7 ,%%ymm7 ,%%ymm7	\n\t	vaddpd	%%ymm4 ,%%ymm4 ,%%ymm4 		\n\t"\
		"vaddpd	%%ymm0 ,%%ymm13,%%ymm13		\n\t	vaddpd	%%ymm2 ,%%ymm15,%%ymm15	\n\t	vaddpd	%%ymm3 ,%%ymm12,%%ymm12		\n\t"\
		"vaddpd	%%ymm8 ,%%ymm5 ,%%ymm5		\n\t	vaddpd	%%ymm10,%%ymm7 ,%%ymm7	\n\t	vaddpd	%%ymm11,%%ymm4 ,%%ymm4 		\n\t"\
		"vmovaps	%%ymm0 ,0x100(%%rdi)	\n\t	vmovaps	%%ymm2 ,0x200(%%rdi)	\n\t	vmovaps	%%ymm3 ,0x400(%%rdi)		\n\t"\
		"vmovaps	%%ymm8 ,0x620(%%rdi)	\n\t	vmovaps	%%ymm10,0x520(%%rdi)	\n\t	vmovaps	%%ymm11,0x320(%%rdi)		\n\t"\
		"vmovaps	%%ymm13,0x600(%%rdi)	\n\t	vmovaps	%%ymm15,0x500(%%rdi)	\n\t	vmovaps	%%ymm12,0x300(%%rdi)		\n\t"\
		"vmovaps	%%ymm5 ,0x120(%%rdi)	\n\t	vmovaps	%%ymm7 ,0x220(%%rdi)	\n\t	vmovaps	%%ymm4 ,0x420(%%rdi)		\n\t"\
		"\n\t"\
	"/* Since doing radix-7 in-place here, inputs of radix-4 are in consecutive memory locs, i.e. 0x40 bytes apart, e.g. the distance between s1p00r and s1p01r: */\n\t"\
		"\n\t"\
	"/* 1:SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_E(s1p00,01,02,03 -> add00+p0123): */	\n\t"\
		"movq	%[__out],%%rsi	/* s1p00r */									\n\t"\
		"movq	%[__add],%%rax													\n\t"\
		"movslq	%[__p01],%%rbx													\n\t"\
		"movslq	%[__p02],%%rcx													\n\t"\
		"movslq	%[__p03],%%rdx													\n\t"\
		"shlq	$3,%%rbx		/* Pointer offset for floating doubles */		\n\t"\
		"shlq	$3,%%rcx														\n\t"\
		"shlq	$3,%%rdx														\n\t"\
		"addq	%%rax,%%rbx														\n\t"\
		"addq	%%rax,%%rcx														\n\t"\
		"addq	%%rax,%%rdx														\n\t"\
		"vmovaps	     (%%rsi),%%ymm0											\n\t"\
		"vmovaps	0x040(%%rsi),%%ymm4											\n\t"\
		"vmovaps	0x020(%%rsi),%%ymm1											\n\t"\
		"vmovaps	0x060(%%rsi),%%ymm5											\n\t"\
		"vmovaps	%%ymm0,%%ymm2												\n\t"\
		"vmovaps	%%ymm4,%%ymm6												\n\t"\
		"vmovaps	%%ymm1,%%ymm3												\n\t"\
		"vmovaps	%%ymm5,%%ymm7												\n\t"\
		"vaddpd	0x080(%%rsi),%%ymm0,%%ymm0										\n\t"\
		"vaddpd	0x0c0(%%rsi),%%ymm4,%%ymm4										\n\t"\
		"vaddpd	0x0a0(%%rsi),%%ymm1,%%ymm1										\n\t"\
		"vaddpd	0x0e0(%%rsi),%%ymm5,%%ymm5										\n\t"\
		"vsubpd	0x080(%%rsi),%%ymm2,%%ymm2										\n\t	/* 2:SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_E(s1p04,05,06,07 -> add24+p1032): */\n\t"\
		"vsubpd	0x0c0(%%rsi),%%ymm6,%%ymm6										\n\t		movslq	%[__p24],%%rdi											\n\t"\
		"vsubpd	0x0a0(%%rsi),%%ymm3,%%ymm3										\n\t		shlq	$3,%%rdi												\n\t"\
		"vsubpd	0x0e0(%%rsi),%%ymm7,%%ymm7										\n\t		addq	$0x100,%%rsi											\n\t"\
		"/* Finish radix-4 butterfly and store: */								\n\t		vmovaps	     (%%rsi),%%ymm8 									\n\t"\
		"vsubpd	%%ymm4,%%ymm0,%%ymm0											\n\t		vmovaps	0x040(%%rsi),%%ymm12									\n\t"\
		"vsubpd	%%ymm7,%%ymm2,%%ymm2											\n\t		vmovaps	0x020(%%rsi),%%ymm9 									\n\t"\
		"vsubpd	%%ymm5,%%ymm1,%%ymm1											\n\t		vmovaps	0x060(%%rsi),%%ymm13									\n\t"\
		"vsubpd	%%ymm6,%%ymm3,%%ymm3											\n\t		vmovaps	%%ymm8 ,%%ymm10											\n\t"\
		"vmovaps	%%ymm0,     (%%rbx)											\n\t		vmovaps	%%ymm12,%%ymm14											\n\t"\
		"vmovaps	%%ymm2,     (%%rcx)											\n\t		vmovaps	%%ymm9 ,%%ymm11											\n\t"\
		"vmovaps	%%ymm1,0x020(%%rbx)											\n\t		vmovaps	%%ymm13,%%ymm15											\n\t"\
		"vmovaps	%%ymm3,0x020(%%rdx)											\n\t		vaddpd	0x080(%%rsi),%%ymm8 ,%%ymm8 							\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4											\n\t		vaddpd	0x0c0(%%rsi),%%ymm12,%%ymm12							\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7											\n\t		vaddpd	0x0a0(%%rsi),%%ymm9 ,%%ymm9 							\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5											\n\t		vaddpd	0x0e0(%%rsi),%%ymm13,%%ymm13							\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6											\n\t		vsubpd	0x080(%%rsi),%%ymm10,%%ymm10							\n\t"\
		"vaddpd	%%ymm0,%%ymm4,%%ymm4											\n\t		vsubpd	0x0c0(%%rsi),%%ymm14,%%ymm14							\n\t"\
		"vaddpd	%%ymm2,%%ymm7,%%ymm7											\n\t		vsubpd	0x0a0(%%rsi),%%ymm11,%%ymm11							\n\t"\
		"vaddpd	%%ymm1,%%ymm5,%%ymm5											\n\t		vsubpd	0x0e0(%%rsi),%%ymm15,%%ymm15							\n\t"\
		"vaddpd	%%ymm3,%%ymm6,%%ymm6											\n\t		/* Finish radix-4 butterfly and store: */						\n\t"\
		"vmovaps	%%ymm4,     (%%rax)											\n\t		addq	%%rdi,%%rbx		/* reorder to not conflict with left-col store addresses */\n\t"\
		"vmovaps	%%ymm7,     (%%rdx)											\n\t		addq	%%rdi,%%rdx												\n\t"\
		"vmovaps	%%ymm5,0x020(%%rax)											\n\t		addq	%%rdi,%%rax												\n\t"\
		"vmovaps	%%ymm6,0x020(%%rcx)											\n\t		addq	%%rdi,%%rcx												\n\t"\
	"/* 3:SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_E(s1p08,09,10,11 -> add20+p2310): */	\n\t		/* Swap _o1 <-> _o2, _o3 <-> _o4: eax <-> ebx, ecx <-> edx */	\n\t"\
		"addq	$0x100,%%rsi													\n\t		vsubpd	%%ymm12,%%ymm8 ,%%ymm8 									\n\t"\
		"movslq	%[__p04],%%rdi													\n\t		vsubpd	%%ymm15,%%ymm10,%%ymm10									\n\t"\
		"shlq	$3,%%rdi														\n\t		vsubpd	%%ymm13,%%ymm9 ,%%ymm9 									\n\t"\
		"vmovaps	     (%%rsi),%%ymm0											\n\t		vsubpd	%%ymm14,%%ymm11,%%ymm11									\n\t"\
		"vmovaps	0x040(%%rsi),%%ymm4											\n\t		vmovaps	%%ymm8 ,     (%%rax)									\n\t"\
		"vmovaps	0x020(%%rsi),%%ymm1											\n\t		vmovaps	%%ymm10,     (%%rdx)									\n\t"\
		"vmovaps	0x060(%%rsi),%%ymm5											\n\t		vmovaps	%%ymm9 ,0x020(%%rax)									\n\t"\
		"vmovaps	%%ymm0,%%ymm2												\n\t		vmovaps	%%ymm11,0x020(%%rcx)									\n\t"\
		"vmovaps	%%ymm4,%%ymm6												\n\t		vaddpd	%%ymm12,%%ymm12,%%ymm12									\n\t"\
		"vmovaps	%%ymm1,%%ymm3												\n\t		vaddpd	%%ymm15,%%ymm15,%%ymm15									\n\t"\
		"vmovaps	%%ymm5,%%ymm7												\n\t		vaddpd	%%ymm13,%%ymm13,%%ymm13									\n\t"\
		"vaddpd	0x080(%%rsi),%%ymm0,%%ymm0										\n\t		vaddpd	%%ymm14,%%ymm14,%%ymm14									\n\t"\
		"vaddpd	0x0c0(%%rsi),%%ymm4,%%ymm4										\n\t		vaddpd	%%ymm8 ,%%ymm12,%%ymm12									\n\t"\
		"vaddpd	0x0a0(%%rsi),%%ymm1,%%ymm1										\n\t		vaddpd	%%ymm10,%%ymm15,%%ymm15									\n\t"\
		"vaddpd	0x0e0(%%rsi),%%ymm5,%%ymm5										\n\t		vaddpd	%%ymm9 ,%%ymm13,%%ymm13									\n\t"\
		"vsubpd	0x080(%%rsi),%%ymm2,%%ymm2										\n\t		vaddpd	%%ymm11,%%ymm14,%%ymm14									\n\t"\
		"vsubpd	0x0c0(%%rsi),%%ymm6,%%ymm6										\n\t		vmovaps	%%ymm12,     (%%rbx)									\n\t"\
		"vsubpd	0x0a0(%%rsi),%%ymm3,%%ymm3										\n\t		vmovaps	%%ymm15,     (%%rcx)									\n\t"\
		"vsubpd	0x0e0(%%rsi),%%ymm7,%%ymm7										\n\t		vmovaps	%%ymm13,0x020(%%rbx)									\n\t"\
		"/* Finish radix-4 butterfly and store: */								\n\t		vmovaps	%%ymm14,0x020(%%rdx)									\n\t"\
		"subq	%%rdi,%%rax														\n\t"\
		"subq	%%rdi,%%rbx														\n\t"\
		"subq	%%rdi,%%rcx														\n\t	/* 4:SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_E(s1p12,13,14,15 -> add16+p0123): */\n\t"\
		"subq	%%rdi,%%rdx														\n\t		addq	$0x100,%%rsi											\n\t"\
		"/* Swap _o0123 <-> _o2310: e[abcd]x <-> e[cdba]x */					\n\t		vmovaps	     (%%rsi),%%ymm8 									\n\t"\
		"vsubpd	%%ymm4,%%ymm0,%%ymm0											\n\t		vmovaps	0x040(%%rsi),%%ymm12									\n\t"\
		"vsubpd	%%ymm7,%%ymm2,%%ymm2											\n\t		vmovaps	0x020(%%rsi),%%ymm9 									\n\t"\
		"vsubpd	%%ymm5,%%ymm1,%%ymm1											\n\t		vmovaps	0x060(%%rsi),%%ymm13									\n\t"\
		"vsubpd	%%ymm6,%%ymm3,%%ymm3											\n\t		vmovaps	%%ymm8 ,%%ymm10											\n\t"\
		"vmovaps	%%ymm0,     (%%rdx)											\n\t		vmovaps	%%ymm12,%%ymm14											\n\t"\
		"vmovaps	%%ymm2,     (%%rbx)											\n\t		vmovaps	%%ymm9 ,%%ymm11											\n\t"\
		"vmovaps	%%ymm1,0x020(%%rdx)											\n\t		vmovaps	%%ymm13,%%ymm15											\n\t"\
		"vmovaps	%%ymm3,0x020(%%rax)											\n\t		vaddpd	0x080(%%rsi),%%ymm8 ,%%ymm8 							\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4											\n\t		vaddpd	0x0c0(%%rsi),%%ymm12,%%ymm12							\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7											\n\t		vaddpd	0x0a0(%%rsi),%%ymm9 ,%%ymm9 							\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5											\n\t		vaddpd	0x0e0(%%rsi),%%ymm13,%%ymm13							\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6											\n\t		vsubpd	0x080(%%rsi),%%ymm10,%%ymm10							\n\t"\
		"vaddpd	%%ymm0,%%ymm4,%%ymm4											\n\t		vsubpd	0x0c0(%%rsi),%%ymm14,%%ymm14							\n\t"\
		"vaddpd	%%ymm2,%%ymm7,%%ymm7											\n\t		vsubpd	0x0a0(%%rsi),%%ymm11,%%ymm11							\n\t"\
		"vaddpd	%%ymm1,%%ymm5,%%ymm5											\n\t		vsubpd	0x0e0(%%rsi),%%ymm15,%%ymm15							\n\t"\
		"vaddpd	%%ymm3,%%ymm6,%%ymm6											\n\t		/* Finish radix-4 butterfly and store: */						\n\t"\
		"vmovaps	%%ymm4,     (%%rcx)											\n\t		subq	%%rdi,%%rdx												\n\t"\
		"vmovaps	%%ymm7,     (%%rax)											\n\t		subq	%%rdi,%%rax												\n\t"\
		"vmovaps	%%ymm5,0x020(%%rcx)											\n\t		subq	%%rdi,%%rcx												\n\t"\
		"vmovaps	%%ymm6,0x020(%%rbx)											\n\t		subq	%%rdi,%%rbx												\n\t"\
	"/* 5:SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_E(s1p16,17,18,19 -> add12+p3201): */	\n\t		vsubpd	%%ymm12,%%ymm8 ,%%ymm8 									\n\t"\
		"addq	$0x100,%%rsi													\n\t		vsubpd	%%ymm15,%%ymm10,%%ymm10									\n\t"\
		"vmovaps	     (%%rsi),%%ymm0											\n\t		vsubpd	%%ymm13,%%ymm9 ,%%ymm9 									\n\t"\
		"vmovaps	0x040(%%rsi),%%ymm4											\n\t		vsubpd	%%ymm14,%%ymm11,%%ymm11									\n\t"\
		"vmovaps	0x020(%%rsi),%%ymm1											\n\t		vmovaps	%%ymm8 ,     (%%rbx)									\n\t"\
		"vmovaps	0x060(%%rsi),%%ymm5											\n\t		vmovaps	%%ymm10,     (%%rcx)									\n\t"\
		"vmovaps	%%ymm0,%%ymm2												\n\t		vmovaps	%%ymm9 ,0x020(%%rbx)									\n\t"\
		"vmovaps	%%ymm4,%%ymm6												\n\t		vmovaps	%%ymm11,0x020(%%rdx)									\n\t"\
		"vmovaps	%%ymm1,%%ymm3												\n\t		vaddpd	%%ymm12,%%ymm12,%%ymm12									\n\t"\
		"vmovaps	%%ymm5,%%ymm7												\n\t		vaddpd	%%ymm15,%%ymm15,%%ymm15									\n\t"\
		"vaddpd	0x080(%%rsi),%%ymm0,%%ymm0										\n\t		vaddpd	%%ymm13,%%ymm13,%%ymm13									\n\t"\
		"vaddpd	0x0c0(%%rsi),%%ymm4,%%ymm4										\n\t		vaddpd	%%ymm14,%%ymm14,%%ymm14									\n\t"\
		"vaddpd	0x0a0(%%rsi),%%ymm1,%%ymm1										\n\t		vaddpd	%%ymm8 ,%%ymm12,%%ymm12									\n\t"\
		"vaddpd	0x0e0(%%rsi),%%ymm5,%%ymm5										\n\t		vaddpd	%%ymm10,%%ymm15,%%ymm15									\n\t"\
		"vsubpd	0x080(%%rsi),%%ymm2,%%ymm2										\n\t		vaddpd	%%ymm9 ,%%ymm13,%%ymm13									\n\t"\
		"vsubpd	0x0c0(%%rsi),%%ymm6,%%ymm6										\n\t		vaddpd	%%ymm11,%%ymm14,%%ymm14									\n\t"\
		"vsubpd	0x0a0(%%rsi),%%ymm3,%%ymm3										\n\t		vmovaps	%%ymm12,     (%%rax)									\n\t"\
		"vsubpd	0x0e0(%%rsi),%%ymm7,%%ymm7										\n\t		vmovaps	%%ymm15,     (%%rdx)									\n\t"\
		"/* Finish radix-4 butterfly and store: */								\n\t		vmovaps	%%ymm13,0x020(%%rax)									\n\t"\
		"subq	%%rdi,%%rax														\n\t		vmovaps	%%ymm14,0x020(%%rcx)									\n\t"\
		"subq	%%rdi,%%rbx														\n\t"\
		"subq	%%rdi,%%rcx														\n\t	/* 6:SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_E(s1p20,21,22,23 -> add08+p1032): */\n\t"\
		"subq	%%rdi,%%rdx														\n\t		addq	$0x100,%%rsi											\n\t"\
		"/* Swap _o0123 <-> _o3201: e[abcd]x <-> e[dcab]x */					\n\t		vmovaps	     (%%rsi),%%ymm8 									\n\t"\
		"vsubpd	%%ymm4,%%ymm0,%%ymm0											\n\t		vmovaps	0x040(%%rsi),%%ymm12									\n\t"\
		"vsubpd	%%ymm7,%%ymm2,%%ymm2											\n\t		vmovaps	0x020(%%rsi),%%ymm9 									\n\t"\
		"vsubpd	%%ymm5,%%ymm1,%%ymm1											\n\t		vmovaps	0x060(%%rsi),%%ymm13									\n\t"\
		"vsubpd	%%ymm6,%%ymm3,%%ymm3											\n\t		vmovaps	%%ymm8 ,%%ymm10											\n\t"\
		"vmovaps	%%ymm0,     (%%rcx)											\n\t		vmovaps	%%ymm12,%%ymm14											\n\t"\
		"vmovaps	%%ymm2,     (%%rax)											\n\t		vmovaps	%%ymm9 ,%%ymm11											\n\t"\
		"vmovaps	%%ymm1,0x020(%%rcx)											\n\t		vmovaps	%%ymm13,%%ymm15											\n\t"\
		"vmovaps	%%ymm3,0x020(%%rbx)											\n\t		vaddpd	0x080(%%rsi),%%ymm8 ,%%ymm8 							\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4											\n\t		vaddpd	0x0c0(%%rsi),%%ymm12,%%ymm12							\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7											\n\t		vaddpd	0x0a0(%%rsi),%%ymm9 ,%%ymm9 							\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5											\n\t		vaddpd	0x0e0(%%rsi),%%ymm13,%%ymm13							\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6											\n\t		vsubpd	0x080(%%rsi),%%ymm10,%%ymm10							\n\t"\
		"vaddpd	%%ymm0,%%ymm4,%%ymm4											\n\t		vsubpd	0x0c0(%%rsi),%%ymm14,%%ymm14							\n\t"\
		"vaddpd	%%ymm2,%%ymm7,%%ymm7											\n\t		vsubpd	0x0a0(%%rsi),%%ymm11,%%ymm11							\n\t"\
		"vaddpd	%%ymm1,%%ymm5,%%ymm5											\n\t		vsubpd	0x0e0(%%rsi),%%ymm15,%%ymm15							\n\t"\
		"vaddpd	%%ymm3,%%ymm6,%%ymm6											\n\t		/* Finish radix-4 butterfly and store: */						\n\t"\
		"vmovaps	%%ymm4,     (%%rdx)											\n\t		subq	%%rdi,%%rcx												\n\t"\
		"vmovaps	%%ymm7,     (%%rbx)											\n\t		subq	%%rdi,%%rbx												\n\t"\
		"vmovaps	%%ymm5,0x020(%%rdx)											\n\t		subq	%%rdi,%%rdx												\n\t"\
		"vmovaps	%%ymm6,0x020(%%rax)											\n\t		subq	%%rdi,%%rax												\n\t"\
	"/* 7:SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_E(s1p24,25,26,27 -> add04+p2310): */	\n\t		/* Swap _o0 <-> _o1, _o2 <-> _o3: eax <-> ebx, ecx <-> edx */	\n\t"\
		"addq	$0x100,%%rsi													\n\t		vsubpd	%%ymm12,%%ymm8 ,%%ymm8 									\n\t"\
		"vmovaps	     (%%rsi),%%ymm0											\n\t		vsubpd	%%ymm15,%%ymm10,%%ymm10									\n\t"\
		"vmovaps	0x040(%%rsi),%%ymm4											\n\t		vsubpd	%%ymm13,%%ymm9 ,%%ymm9 									\n\t"\
		"vmovaps	0x020(%%rsi),%%ymm1											\n\t		vsubpd	%%ymm14,%%ymm11,%%ymm11									\n\t"\
		"vmovaps	0x060(%%rsi),%%ymm5											\n\t		vmovaps	%%ymm8 ,     (%%rax)									\n\t"\
		"vmovaps	%%ymm0,%%ymm2												\n\t		vmovaps	%%ymm10,     (%%rdx)									\n\t"\
		"vmovaps	%%ymm4,%%ymm6												\n\t		vmovaps	%%ymm9 ,0x020(%%rax)									\n\t"\
		"vmovaps	%%ymm1,%%ymm3												\n\t		vmovaps	%%ymm11,0x020(%%rcx)									\n\t"\
		"vmovaps	%%ymm5,%%ymm7												\n\t		vaddpd	%%ymm12,%%ymm12,%%ymm12									\n\t"\
		"vaddpd	0x080(%%rsi),%%ymm0,%%ymm0										\n\t		vaddpd	%%ymm15,%%ymm15,%%ymm15									\n\t"\
		"vaddpd	0x0c0(%%rsi),%%ymm4,%%ymm4										\n\t		vaddpd	%%ymm13,%%ymm13,%%ymm13									\n\t"\
		"vaddpd	0x0a0(%%rsi),%%ymm1,%%ymm1										\n\t		vaddpd	%%ymm14,%%ymm14,%%ymm14									\n\t"\
		"vaddpd	0x0e0(%%rsi),%%ymm5,%%ymm5										\n\t		vaddpd	%%ymm8 ,%%ymm12,%%ymm12									\n\t"\
		"vsubpd	0x080(%%rsi),%%ymm2,%%ymm2										\n\t		vaddpd	%%ymm10,%%ymm15,%%ymm15									\n\t"\
		"vsubpd	0x0c0(%%rsi),%%ymm6,%%ymm6										\n\t		vaddpd	%%ymm9 ,%%ymm13,%%ymm13									\n\t"\
		"vsubpd	0x0a0(%%rsi),%%ymm3,%%ymm3										\n\t		vaddpd	%%ymm11,%%ymm14,%%ymm14									\n\t"\
		"vsubpd	0x0e0(%%rsi),%%ymm7,%%ymm7										\n\t		vmovaps	%%ymm12,     (%%rbx)									\n\t"\
		"/* Finish radix-4 butterfly and store: */								\n\t		vmovaps	%%ymm15,     (%%rcx)									\n\t"\
		"subq	%%rdi,%%rax														\n\t		vmovaps	%%ymm13,0x020(%%rbx)									\n\t"\
		"subq	%%rdi,%%rcx														\n\t		vmovaps	%%ymm14,0x020(%%rdx)									\n\t"\
		"subq	%%rdi,%%rbx														\n\t"\
		"subq	%%rdi,%%rdx														\n\t"\
		"/* Swap _o0123 <-> _o2310: e[abcd]x <-> e[cdba]x */					\n\t"\
		"vsubpd	%%ymm4,%%ymm0,%%ymm0											\n\t"\
		"vsubpd	%%ymm7,%%ymm2,%%ymm2											\n\t"\
		"vsubpd	%%ymm5,%%ymm1,%%ymm1											\n\t"\
		"vsubpd	%%ymm6,%%ymm3,%%ymm3											\n\t"\
		"vmovaps	%%ymm0,     (%%rdx)											\n\t"\
		"vmovaps	%%ymm2,     (%%rbx)											\n\t"\
		"vmovaps	%%ymm1,0x020(%%rdx)											\n\t"\
		"vmovaps	%%ymm3,0x020(%%rax)											\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4											\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7											\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5											\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6											\n\t"\
		"vaddpd	%%ymm0,%%ymm4,%%ymm4											\n\t"\
		"vaddpd	%%ymm2,%%ymm7,%%ymm7											\n\t"\
		"vaddpd	%%ymm1,%%ymm5,%%ymm5											\n\t"\
		"vaddpd	%%ymm3,%%ymm6,%%ymm6											\n\t"\
		"vmovaps	%%ymm4,     (%%rcx)											\n\t"\
		"vmovaps	%%ymm7,     (%%rax)											\n\t"\
		"vmovaps	%%ymm5,0x020(%%rcx)											\n\t"\
		"vmovaps	%%ymm6,0x020(%%rbx)											\n\t"\
		:					/* outputs: none */\
		: [__add] "m" (Xadd)	/* All inputs from memory addresses here */\
		 ,[__p01] "m" (Xp01)\
		 ,[__p02] "m" (Xp02)\
		 ,[__p03] "m" (Xp03)\
		 ,[__p04] "m" (Xp04)\
		 ,[__p08] "m" (Xp08)\
		 ,[__p12] "m" (Xp12)\
		 ,[__p16] "m" (Xp16)\
		 ,[__p20] "m" (Xp20)\
		 ,[__p24] "m" (Xp24)\
		 ,[__out] "m" (Xout)\
		 ,[__cc0] "m" (Xcc0)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"		/* Clobbered registers */\
	);\
	}

#elif defined(USE_SSE2)

	#define	SSE2_RADIX28_DIT_NOTWIDDLE(Xadd,Xp01,Xp02,Xp03,Xp04,Xp08,Xp12,Xp16,Xp20,Xp24,Xout,Xcc0)\
	{\
	__asm__ volatile (\
	"/* 1:SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add00+p0132, s1p00,03,02,01): */	\n\t"\
		"movq	%[__out],%%rsi										\n\t"\
		"movq	%[__add],%%rax										\n\t"\
		"movslq	%[__p01],%%rbx										\n\t"\
		"movslq	%[__p02],%%rcx										\n\t"\
		"movslq	%[__p03],%%rdx										\n\t"\
		"shlq	$3,%%rbx											\n\t"\
		"shlq	$3,%%rcx											\n\t"\
		"shlq	$3,%%rdx											\n\t"\
		"addq	%%rax,%%rbx											\n\t"\
		"addq	%%rax,%%rcx											\n\t"\
		"addq	%%rax,%%rdx											\n\t"\
		"/* ecx <-> edx */											\n\t"\
		"movaps	     (%%rax),%%xmm0									\n\t"\
		"movaps	     (%%rdx),%%xmm4									\n\t"\
		"movaps	0x010(%%rax),%%xmm1									\n\t"\
		"movaps	0x010(%%rdx),%%xmm5									\n\t"\
		"movaps	%%xmm0,%%xmm2										\n\t"\
		"movaps	%%xmm4,%%xmm6										\n\t"\
		"movaps	%%xmm1,%%xmm3										\n\t"\
		"movaps	%%xmm5,%%xmm7										\n\t"\
		"addpd	     (%%rbx),%%xmm0									\n\t"\
		"addpd	     (%%rcx),%%xmm4									\n\t"\
		"addpd	0x010(%%rbx),%%xmm1									\n\t"\
		"addpd	0x010(%%rcx),%%xmm5									\n\t"\
		"subpd	     (%%rbx),%%xmm2									\n\t"\
		"subpd	     (%%rcx),%%xmm6									\n\t/* 2:SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add12+p3210, s1p04,07,06,05r): */	\n\t"\
		"subpd	0x010(%%rbx),%%xmm3									\n\t	movslq	%[__p12],%%rdi								\n\t"\
		"subpd	0x010(%%rcx),%%xmm7									\n\t	shlq	$3,%%rdi									\n\t"\
		"/* Finish radix-4 butterfly and store results: */			\n\t	addq	%%rdi,%%rax		/* add0 + p12 */			\n\t"\
		"/* Swap _o1 <-> _o3: 0x020,0x030 <-> 0x060,0x070: */		\n\t	addq	%%rdi,%%rbx									\n\t"\
		"subpd	%%xmm4,%%xmm0										\n\t	addq	%%rdi,%%rcx									\n\t"\
		"subpd	%%xmm7,%%xmm2										\n\t	addq	%%rdi,%%rdx									\n\t"\
		"subpd	%%xmm5,%%xmm1										\n\t	/* eax <-> edx, ebx <-> ecx */						\n\t"\
		"subpd	%%xmm6,%%xmm3										\n\t	movaps	     (%%rdx),%%xmm8 						\n\t"\
		"movaps	%%xmm0,0x040(%%rsi)									\n\t	movaps	     (%%rbx),%%xmm12						\n\t"\
		"movaps	%%xmm2,0x020(%%rsi)									\n\t	movaps	0x010(%%rdx),%%xmm9 						\n\t"\
		"movaps	%%xmm1,0x050(%%rsi)									\n\t	movaps	0x010(%%rbx),%%xmm13						\n\t"\
		"movaps	%%xmm3,0x070(%%rsi)									\n\t	movaps	%%xmm8 ,%%xmm10								\n\t"\
		"addpd	%%xmm4,%%xmm4										\n\t	movaps	%%xmm12,%%xmm14								\n\t"\
		"addpd	%%xmm7,%%xmm7										\n\t	movaps	%%xmm9 ,%%xmm11								\n\t"\
		"addpd	%%xmm5,%%xmm5										\n\t	movaps	%%xmm13,%%xmm15								\n\t"\
		"addpd	%%xmm6,%%xmm6										\n\t	addpd	     (%%rcx),%%xmm8  						\n\t"\
		"addpd	%%xmm0,%%xmm4										\n\t	addpd	     (%%rax),%%xmm12						\n\t"\
		"addpd	%%xmm2,%%xmm7										\n\t	addpd	0x010(%%rcx),%%xmm9  						\n\t"\
		"addpd	%%xmm1,%%xmm5										\n\t	addpd	0x010(%%rax),%%xmm13						\n\t"\
		"addpd	%%xmm3,%%xmm6										\n\t	subpd	     (%%rcx),%%xmm10						\n\t"\
		"movaps	%%xmm4,     (%%rsi)									\n\t	subpd	     (%%rax),%%xmm14						\n\t"\
		"movaps	%%xmm7,0x060(%%rsi)									\n\t	subpd	0x010(%%rcx),%%xmm11						\n\t"\
		"movaps	%%xmm5,0x010(%%rsi)									\n\t	subpd	0x010(%%rax),%%xmm15						\n\t"\
		"movaps	%%xmm6,0x030(%%rsi)									\n\t	/* Finish radix-4 butterfly and store results: */	\n\t"\
	"/* 3:SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add24+p1023, s1p08,11,10,09r): */	\n\t	addq	$0x080,%%rsi						\n\t"\
		"addq	%%rdi,%%rax		/* add0 + p24 */					\n\t	subpd	%%xmm12,%%xmm8  							\n\t"\
		"addq	%%rdi,%%rbx											\n\t	subpd	%%xmm15,%%xmm10								\n\t"\
		"addq	%%rdi,%%rcx											\n\t	subpd	%%xmm13,%%xmm9  							\n\t"\
		"addq	%%rdi,%%rdx											\n\t	subpd	%%xmm14,%%xmm11								\n\t"\
		"/* eax <-> ebx */											\n\t	movaps	%%xmm8 ,0x040(%%rsi)						\n\t"\
		"movaps	     (%%rbx),%%xmm0									\n\t	movaps	%%xmm10,0x020(%%rsi)						\n\t"\
		"movaps	     (%%rcx),%%xmm4									\n\t	movaps	%%xmm9 ,0x050(%%rsi)						\n\t"\
		"movaps	0x010(%%rbx),%%xmm1									\n\t	movaps	%%xmm11,0x070(%%rsi)						\n\t"\
		"movaps	0x010(%%rcx),%%xmm5									\n\t	addpd	%%xmm12,%%xmm12								\n\t"\
		"movaps	%%xmm0,%%xmm2										\n\t	addpd	%%xmm15,%%xmm15								\n\t"\
		"movaps	%%xmm4,%%xmm6										\n\t	addpd	%%xmm13,%%xmm13								\n\t"\
		"movaps	%%xmm1,%%xmm3										\n\t	addpd	%%xmm14,%%xmm14								\n\t"\
		"movaps	%%xmm5,%%xmm7										\n\t	addpd	%%xmm8 ,%%xmm12								\n\t"\
		"addpd	     (%%rax),%%xmm0									\n\t	addpd	%%xmm10,%%xmm15								\n\t"\
		"addpd	     (%%rdx),%%xmm4									\n\t	addpd	%%xmm9 ,%%xmm13								\n\t"\
		"addpd	0x010(%%rax),%%xmm1									\n\t	addpd	%%xmm11,%%xmm14								\n\t"\
		"addpd	0x010(%%rdx),%%xmm5									\n\t	movaps	%%xmm12,     (%%rsi)						\n\t"\
		"subpd	     (%%rax),%%xmm2									\n\t	movaps	%%xmm15,0x060(%%rsi)						\n\t"\
		"subpd	     (%%rdx),%%xmm6									\n\t	movaps	%%xmm13,0x010(%%rsi)						\n\t"\
		"subpd	0x010(%%rax),%%xmm3									\n\t	movaps	%%xmm14,0x030(%%rsi)						\n\t"\
		"subpd	0x010(%%rdx),%%xmm7									\n\t/* 4:SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add08+p1023, s1p12,15,14,13r): */	\n\t"\
		"/* Finish radix-4 butterfly and store results: */			\n\t	movslq	%[__p16],%%rdi								\n\t"\
		"addq	$0x080,%%rsi										\n\t	shlq	$3,%%rdi									\n\t"\
		"subpd	%%xmm4,%%xmm0										\n\t	subq	%%rdi,%%rax		/* add0 + p08 */			\n\t"\
		"subpd	%%xmm7,%%xmm2										\n\t	subq	%%rdi,%%rbx									\n\t"\
		"subpd	%%xmm5,%%xmm1										\n\t	subq	%%rdi,%%rcx									\n\t"\
		"subpd	%%xmm6,%%xmm3										\n\t	subq	%%rdi,%%rdx									\n\t"\
		"movaps	%%xmm0,0x040(%%rsi)									\n\t	/* eax <-> %%rbx */									\n\t"\
		"movaps	%%xmm2,0x020(%%rsi)									\n\t	movaps	     (%%rbx),%%xmm8 						\n\t"\
		"movaps	%%xmm1,0x050(%%rsi)									\n\t	movaps	     (%%rcx),%%xmm12						\n\t"\
		"movaps	%%xmm3,0x070(%%rsi)									\n\t	movaps	0x010(%%rbx),%%xmm9 						\n\t"\
		"addpd	%%xmm4,%%xmm4										\n\t	movaps	0x010(%%rcx),%%xmm13						\n\t"\
		"addpd	%%xmm7,%%xmm7										\n\t	movaps	%%xmm8 ,%%xmm10								\n\t"\
		"addpd	%%xmm5,%%xmm5										\n\t	movaps	%%xmm12,%%xmm14								\n\t"\
		"addpd	%%xmm6,%%xmm6										\n\t	movaps	%%xmm9 ,%%xmm11								\n\t"\
		"addpd	%%xmm0,%%xmm4										\n\t	movaps	%%xmm13,%%xmm15								\n\t"\
		"addpd	%%xmm2,%%xmm7										\n\t	addpd	     (%%rax),%%xmm8  						\n\t"\
		"addpd	%%xmm1,%%xmm5										\n\t	addpd	     (%%rdx),%%xmm12						\n\t"\
		"addpd	%%xmm3,%%xmm6										\n\t	addpd	0x010(%%rax),%%xmm9  						\n\t"\
		"movaps	%%xmm4,     (%%rsi)									\n\t	addpd	0x010(%%rdx),%%xmm13						\n\t"\
		"movaps	%%xmm7,0x060(%%rsi)									\n\t	subpd	     (%%rax),%%xmm10						\n\t"\
		"movaps	%%xmm5,0x010(%%rsi)									\n\t	subpd	     (%%rdx),%%xmm14						\n\t"\
		"movaps	%%xmm6,0x030(%%rsi)									\n\t	subpd	0x010(%%rax),%%xmm11						\n\t"\
	"/* 5:SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add20+p2301, s1p16,19,18,17r): */	\n\t	subpd	0x010(%%rdx),%%xmm15				\n\t"\
		"movslq	%[__p12],%%rdi										\n\t	/* Finish radix-4 butterfly and store results: */	\n\t"\
		"shlq	$3,%%rdi											\n\t	addq	$0x080,%%rsi								\n\t"\
		"addq	%%rdi,%%rax		/* add0 + p20 */					\n\t	subpd	%%xmm12,%%xmm8  							\n\t"\
		"addq	%%rdi,%%rbx											\n\t	subpd	%%xmm15,%%xmm10								\n\t"\
		"addq	%%rdi,%%rcx											\n\t	subpd	%%xmm13,%%xmm9  							\n\t"\
		"addq	%%rdi,%%rdx											\n\t	subpd	%%xmm14,%%xmm11								\n\t"\
		"/* eax <-> %%rcx, %%rbx <-> edx */							\n\t	movaps	%%xmm8 ,0x040(%%rsi)						\n\t"\
		"movaps	     (%%rcx),%%xmm0									\n\t	movaps	%%xmm10,0x020(%%rsi)						\n\t"\
		"movaps	     (%%rax),%%xmm4									\n\t	movaps	%%xmm9 ,0x050(%%rsi)						\n\t"\
		"movaps	0x010(%%rcx),%%xmm1									\n\t	movaps	%%xmm11,0x070(%%rsi)						\n\t"\
		"movaps	0x010(%%rax),%%xmm5									\n\t	addpd	%%xmm12,%%xmm12								\n\t"\
		"movaps	%%xmm0,%%xmm2										\n\t	addpd	%%xmm15,%%xmm15								\n\t"\
		"movaps	%%xmm4,%%xmm6										\n\t	addpd	%%xmm13,%%xmm13								\n\t"\
		"movaps	%%xmm1,%%xmm3										\n\t	addpd	%%xmm14,%%xmm14								\n\t"\
		"movaps	%%xmm5,%%xmm7										\n\t	addpd	%%xmm8 ,%%xmm12								\n\t"\
		"addpd	     (%%rdx),%%xmm0									\n\t	addpd	%%xmm10,%%xmm15								\n\t"\
		"addpd	     (%%rbx),%%xmm4									\n\t	addpd	%%xmm9 ,%%xmm13								\n\t"\
		"addpd	0x010(%%rdx),%%xmm1									\n\t	addpd	%%xmm11,%%xmm14								\n\t"\
		"addpd	0x010(%%rbx),%%xmm5									\n\t	movaps	%%xmm12,     (%%rsi)						\n\t"\
		"subpd	     (%%rdx),%%xmm2									\n\t	movaps	%%xmm15,0x060(%%rsi)						\n\t"\
		"subpd	     (%%rbx),%%xmm6									\n\t	movaps	%%xmm13,0x010(%%rsi)						\n\t"\
		"subpd	0x010(%%rdx),%%xmm3									\n\t	movaps	%%xmm14,0x030(%%rsi)						\n\t"\
		"subpd	0x010(%%rbx),%%xmm7									\n\t/* 6:SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add04+p2301, s1p20,23,22,21r): */	\n\t"\
		"/* Finish radix-4 butterfly and store results: */			\n\t	movslq	%[__p16],%%rdi								\n\t"\
		"addq	$0x080,%%rsi										\n\t	shlq	$3,%%rdi									\n\t"\
		"subpd	%%xmm4,%%xmm0										\n\t	subq	%%rdi,%%rax		/* add0 + p04 */			\n\t"\
		"subpd	%%xmm7,%%xmm2										\n\t	subq	%%rdi,%%rbx									\n\t"\
		"subpd	%%xmm5,%%xmm1										\n\t	subq	%%rdi,%%rcx									\n\t"\
		"subpd	%%xmm6,%%xmm3										\n\t	subq	%%rdi,%%rdx									\n\t"\
		"movaps	%%xmm0,0x040(%%rsi)									\n\t	/* eax <-> %%rcx, %%rbx <-> edx */					\n\t"\
		"movaps	%%xmm2,0x020(%%rsi)									\n\t	movaps	     (%%rcx),%%xmm8 						\n\t"\
		"movaps	%%xmm1,0x050(%%rsi)									\n\t	movaps	     (%%rax),%%xmm12						\n\t"\
		"movaps	%%xmm3,0x070(%%rsi)									\n\t	movaps	0x010(%%rcx),%%xmm9 						\n\t"\
		"addpd	%%xmm4,%%xmm4										\n\t	movaps	0x010(%%rax),%%xmm13						\n\t"\
		"addpd	%%xmm7,%%xmm7										\n\t	movaps	%%xmm8 ,%%xmm10								\n\t"\
		"addpd	%%xmm5,%%xmm5										\n\t	movaps	%%xmm12,%%xmm14								\n\t"\
		"addpd	%%xmm6,%%xmm6										\n\t	movaps	%%xmm9 ,%%xmm11								\n\t"\
		"addpd	%%xmm0,%%xmm4										\n\t	movaps	%%xmm13,%%xmm15								\n\t"\
		"addpd	%%xmm2,%%xmm7										\n\t	addpd	     (%%rdx),%%xmm8  						\n\t"\
		"addpd	%%xmm1,%%xmm5										\n\t	addpd	     (%%rbx),%%xmm12						\n\t"\
		"addpd	%%xmm3,%%xmm6										\n\t	addpd	0x010(%%rdx),%%xmm9  						\n\t"\
		"movaps	%%xmm4,     (%%rsi)									\n\t	addpd	0x010(%%rbx),%%xmm13						\n\t"\
		"movaps	%%xmm7,0x060(%%rsi)									\n\t	subpd	     (%%rdx),%%xmm10						\n\t"\
		"movaps	%%xmm5,0x010(%%rsi)									\n\t	subpd	     (%%rbx),%%xmm14						\n\t"\
		"movaps	%%xmm6,0x030(%%rsi)									\n\t	subpd	0x010(%%rdx),%%xmm11						\n\t"\
	"/* 7:SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add16+p0132, s1p24,27,26,25r): */	\n\t	subpd	0x010(%%rbx),%%xmm15				\n\t"\
		"movslq	%[__p12],%%rdi										\n\t	/* Finish radix-4 butterfly and store results: */	\n\t"\
		"shlq	$3,%%rdi											\n\t	addq	$0x080,%%rsi								\n\t"\
		"addq	%%rdi,%%rax		/* add0 + p16 */					\n\t	subpd	%%xmm12,%%xmm8  							\n\t"\
		"addq	%%rdi,%%rbx											\n\t	subpd	%%xmm15,%%xmm10								\n\t"\
		"addq	%%rdi,%%rcx											\n\t	subpd	%%xmm13,%%xmm9  							\n\t"\
		"addq	%%rdi,%%rdx											\n\t	subpd	%%xmm14,%%xmm11								\n\t"\
		"/* ecx <-> %%rdx */										\n\t	movaps	%%xmm8 ,0x040(%%rsi)						\n\t"\
		"movaps	     (%%rax),%%xmm0									\n\t	movaps	%%xmm10,0x020(%%rsi)						\n\t"\
		"movaps	     (%%rdx),%%xmm4									\n\t	movaps	%%xmm9 ,0x050(%%rsi)						\n\t"\
		"movaps	0x010(%%rax),%%xmm1									\n\t	movaps	%%xmm11,0x070(%%rsi)						\n\t"\
		"movaps	0x010(%%rdx),%%xmm5									\n\t	addpd	%%xmm12,%%xmm12								\n\t"\
		"movaps	%%xmm0,%%xmm2										\n\t	addpd	%%xmm15,%%xmm15								\n\t"\
		"movaps	%%xmm4,%%xmm6										\n\t	addpd	%%xmm13,%%xmm13								\n\t"\
		"movaps	%%xmm1,%%xmm3										\n\t	addpd	%%xmm14,%%xmm14								\n\t"\
		"movaps	%%xmm5,%%xmm7										\n\t	addpd	%%xmm8 ,%%xmm12								\n\t"\
		"addpd	     (%%rbx),%%xmm0									\n\t	addpd	%%xmm10,%%xmm15								\n\t"\
		"addpd	     (%%rcx),%%xmm4									\n\t	addpd	%%xmm9 ,%%xmm13								\n\t"\
		"addpd	0x010(%%rbx),%%xmm1									\n\t	addpd	%%xmm11,%%xmm14								\n\t"\
		"addpd	0x010(%%rcx),%%xmm5									\n\t	movaps	%%xmm12,     (%%rsi)						\n\t"\
		"subpd	     (%%rbx),%%xmm2									\n\t	movaps	%%xmm15,0x060(%%rsi)						\n\t"\
		"subpd	     (%%rcx),%%xmm6									\n\t	movaps	%%xmm13,0x010(%%rsi)						\n\t"\
		"subpd	0x010(%%rbx),%%xmm3									\n\t	movaps	%%xmm14,0x030(%%rsi)						\n\t"\
		"subpd	0x010(%%rcx),%%xmm7									\n\t"\
		"/* Finish radix-4 butterfly and store results: */			\n\t"\
		"addq	$0x080,%%rsi										\n\t"\
		"subpd	%%xmm4,%%xmm0										\n\t"\
		"subpd	%%xmm7,%%xmm2										\n\t"\
		"subpd	%%xmm5,%%xmm1										\n\t"\
		"subpd	%%xmm6,%%xmm3										\n\t"\
		"movaps	%%xmm0,0x040(%%rsi)									\n\t"\
		"movaps	%%xmm2,0x020(%%rsi)									\n\t"\
		"movaps	%%xmm1,0x050(%%rsi)									\n\t"\
		"movaps	%%xmm3,0x070(%%rsi)									\n\t"\
		"addpd	%%xmm4,%%xmm4										\n\t"\
		"addpd	%%xmm7,%%xmm7										\n\t"\
		"addpd	%%xmm5,%%xmm5										\n\t"\
		"addpd	%%xmm6,%%xmm6										\n\t"\
		"addpd	%%xmm0,%%xmm4										\n\t"\
		"addpd	%%xmm2,%%xmm7										\n\t"\
		"addpd	%%xmm1,%%xmm5										\n\t"\
		"addpd	%%xmm3,%%xmm6										\n\t"\
		"movaps	%%xmm4,     (%%rsi)									\n\t"\
		"movaps	%%xmm7,0x060(%%rsi)									\n\t"\
		"movaps	%%xmm5,0x010(%%rsi)									\n\t"\
		"movaps	%%xmm6,0x030(%%rsi)									\n\t"\
		"/***************************************/					\n\t"\
		"/*...and now do 4 radix-7 transforms...*/					\n\t"\
		"/***************************************/					\n\t"\
	"/* SSE2_RADIX_07_DFT(00,04,08,12,16,20,24 -> 00,08,16,24,04,12,20): */	\n\t"\
		"movq	%[__out],%%rdi			\n\t"\
		"movq	%[__cc0],%%rax			\n\t"\
		"movq	$0x020,%%rbx			\n\t"\
		"movq	$0x040,%%rcx			\n\t"\
		"movq	$0x060,%%rdx			\n\t"\
		"movq	$0x080,%%rsi			\n\t"\
		"addq	%%rax,%%rbx				\n\t"\
		"addq	%%rax,%%rcx				\n\t"\
		"addq	%%rax,%%rdx				\n\t"\
		"addq	%%rax,%%rsi				\n\t	/*** Imaginary Parts: ***/			\n\t"\
		"movaps	0x080(%%rdi),%%xmm6		\n\t	movaps	0x090(%%rdi),%%xmm14		\n\t"\
		"movaps	0x300(%%rdi),%%xmm1		\n\t	movaps	0x310(%%rdi),%%xmm9 		\n\t"\
		"movaps	0x100(%%rdi),%%xmm5		\n\t	movaps	0x110(%%rdi),%%xmm13		\n\t"\
		"movaps	0x280(%%rdi),%%xmm2		\n\t	movaps	0x290(%%rdi),%%xmm10		\n\t"\
		"movaps	0x180(%%rdi),%%xmm4		\n\t	movaps	0x190(%%rdi),%%xmm12		\n\t"\
		"movaps	0x200(%%rdi),%%xmm3		\n\t	movaps	0x210(%%rdi),%%xmm11		\n\t"\
		"subpd	%%xmm1,%%xmm6			\n\t	subpd	%%xmm9 ,%%xmm14				\n\t"\
		"addpd	%%xmm1,%%xmm1			\n\t	addpd	%%xmm9  ,%%xmm9 			\n\t"\
		"addpd	%%xmm6,%%xmm1			\n\t	addpd	%%xmm14,%%xmm9  			\n\t"\
		"subpd	%%xmm2,%%xmm5			\n\t	subpd	%%xmm10,%%xmm13				\n\t"\
		"addpd	%%xmm2,%%xmm2			\n\t	addpd	%%xmm10,%%xmm10				\n\t"\
		"addpd	%%xmm5,%%xmm2			\n\t	addpd	%%xmm13,%%xmm10				\n\t"\
		"movaps	     (%%rdi),%%xmm0		\n\t	movaps	0x010(%%rdi),%%xmm8 		\n\t"\
		"subpd	%%xmm3,%%xmm4			\n\t	subpd	%%xmm11,%%xmm12				\n\t"\
		"addpd	%%xmm3,%%xmm3			\n\t	addpd	%%xmm11,%%xmm11				\n\t"\
		"addpd	%%xmm4,%%xmm3			\n\t	addpd	%%xmm12,%%xmm11				\n\t"\
		"/* Only alloc 2 scratch storageslots in calling code, so 16-reg version use 2 of the eventual outputs for Im parts */\n\t"\
		"movaps	%%xmm0,     (%%rsi)		\n\t	movaps	%%xmm8 ,0x080(%%rdi)		\n\t"\
		"movaps	%%xmm6,0x010(%%rsi)		\n\t	movaps	%%xmm14,0x100(%%rdi)		\n\t"\
		"addpd	      %%xmm1,%%xmm0		\n\t	addpd	      %%xmm9 ,%%xmm8  		\n\t"\
		"movaps	%%xmm5,%%xmm7			\n\t	movaps	%%xmm13,%%xmm15				\n\t"\
		"addpd	      %%xmm2,%%xmm3		\n\t	addpd	      %%xmm10,%%xmm11		\n\t"\
		"subpd	      %%xmm4,%%xmm5		\n\t	subpd	      %%xmm12,%%xmm13		\n\t"\
		"subpd	      %%xmm2,%%xmm1		\n\t	subpd	      %%xmm10,%%xmm9  		\n\t"\
		"subpd	      %%xmm7,%%xmm6		\n\t	subpd	      %%xmm15,%%xmm14		\n\t"\
		"addpd	      %%xmm2,%%xmm2		\n\t	addpd	      %%xmm10,%%xmm10		\n\t"\
		"addpd	      %%xmm7,%%xmm4		\n\t	addpd	      %%xmm15,%%xmm12		\n\t"\
		"addpd	      %%xmm3,%%xmm0		\n\t	addpd	      %%xmm11,%%xmm8  		\n\t"\
		"addpd	0x010(%%rsi),%%xmm5		\n\t	addpd	0x100(%%rdi),%%xmm13		\n\t"\
		"subpd	      %%xmm2,%%xmm3		\n\t	subpd	      %%xmm10,%%xmm11		\n\t"\
		"movaps	%%xmm4,%%xmm7			\n\t	movaps	%%xmm12,%%xmm15				\n\t"\
		"movaps	%%xmm0,      (%%rdi)	\n\t	movaps	%%xmm8 ,0x010(%%rdi)		\n\t"\
		"subpd	      %%xmm6,%%xmm4		\n\t	subpd	      %%xmm14,%%xmm12		\n\t"\
		"movaps	%%xmm1,%%xmm2			\n\t	movaps	%%xmm9 ,%%xmm10				\n\t"\
		"subpd	     (%%rsi),%%xmm0		\n\t	subpd	0x080(%%rdi),%%xmm8  		\n\t"\
		"mulpd	0x010(%%rax),%%xmm5		\n\t	mulpd	0x010(%%rax),%%xmm13		\n\t"\
		"addpd	      %%xmm3,%%xmm2		\n\t	addpd	      %%xmm11,%%xmm10		\n\t"\
		"mulpd	     (%%rcx),%%xmm3		\n\t	mulpd	     (%%rcx),%%xmm11		\n\t"\
		"mulpd	0x010(%%rdx),%%xmm4		\n\t	mulpd	0x010(%%rdx),%%xmm12		\n\t"\
		"mulpd	     (%%rbx),%%xmm1		\n\t	mulpd	     (%%rbx),%%xmm9  		\n\t"\
		"mulpd	0x010(%%rbx),%%xmm6		\n\t	mulpd	0x010(%%rbx),%%xmm14		\n\t"\
		"mulpd	     (%%rax),%%xmm0		\n\t	mulpd	     (%%rax),%%xmm8  		\n\t"\
		"mulpd	0x010(%%rcx),%%xmm7		\n\t	mulpd	0x010(%%rcx),%%xmm15		\n\t"\
		"mulpd	     (%%rdx),%%xmm2		\n\t	mulpd	     (%%rdx),%%xmm10		\n\t"\
		"addpd	     (%%rdi),%%xmm0		\n\t	addpd	0x010(%%rdi),%%xmm8  		\n\t"\
		"addpd	      %%xmm4,%%xmm6		\n\t	addpd	      %%xmm12,%%xmm14		\n\t"\
		"subpd	      %%xmm2,%%xmm1		\n\t	subpd	      %%xmm10,%%xmm9  		\n\t"\
		"subpd	      %%xmm7,%%xmm4		\n\t	subpd	      %%xmm15,%%xmm12		\n\t"\
		"subpd	      %%xmm2,%%xmm3		\n\t	subpd	      %%xmm10,%%xmm11		\n\t"\
		"movaps	%%xmm0,%%xmm2			\n\t	movaps	%%xmm8 ,%%xmm10				\n\t"\
		"movaps	%%xmm5,%%xmm7			\n\t	movaps	%%xmm13,%%xmm15				\n\t"\
		"addpd	      %%xmm1,%%xmm0		\n\t	addpd	      %%xmm9 ,%%xmm8  		\n\t"\
		"addpd	      %%xmm6,%%xmm5		\n\t	addpd	      %%xmm14,%%xmm13		\n\t"\
		"addpd	      %%xmm3,%%xmm1		\n\t	addpd	      %%xmm11,%%xmm9  		\n\t"\
		"addpd	      %%xmm4,%%xmm6		\n\t	addpd	      %%xmm12,%%xmm14		\n\t"\
		"addpd	      %%xmm2,%%xmm3		\n\t	addpd	      %%xmm10,%%xmm11		\n\t"\
		"addpd	      %%xmm7,%%xmm4		\n\t	addpd	      %%xmm15,%%xmm12		\n\t"\
		"subpd	      %%xmm1,%%xmm2		\n\t	subpd	      %%xmm9 ,%%xmm10		\n\t"\
		"subpd	      %%xmm6,%%xmm7		\n\t	subpd	      %%xmm14,%%xmm15		\n\t"\
		"/* xmm1,6,9,14 free ... Note the order reversal on the 3rd pair of outputs: */\n\t"\
		"subpd	%%xmm13,%%xmm0			\n\t	subpd	%%xmm15,%%xmm2			\n\t	subpd	%%xmm12,%%xmm3 			\n\t"\
		"subpd	%%xmm5 ,%%xmm8  		\n\t	subpd	%%xmm7 ,%%xmm10			\n\t	subpd	%%xmm4 ,%%xmm11			\n\t"\
		"addpd	%%xmm13,%%xmm13			\n\t	addpd	%%xmm15,%%xmm15			\n\t	addpd	%%xmm12,%%xmm12			\n\t"\
		"addpd	%%xmm5 ,%%xmm5			\n\t	addpd	%%xmm7 ,%%xmm7			\n\t	addpd	%%xmm4 ,%%xmm4 			\n\t"\
		"addpd	%%xmm0 ,%%xmm13			\n\t	addpd	%%xmm2 ,%%xmm15			\n\t	addpd	%%xmm3 ,%%xmm12			\n\t"\
		"addpd	%%xmm8 ,%%xmm5			\n\t	addpd	%%xmm10,%%xmm7			\n\t	addpd	%%xmm11,%%xmm4 			\n\t"\
		"movaps	%%xmm0 ,0x100(%%rdi)	\n\t	movaps	%%xmm2 ,0x200(%%rdi)	\n\t	movaps	%%xmm3 ,0x080(%%rdi)	\n\t"\
		"movaps	%%xmm8 ,0x290(%%rdi)	\n\t	movaps	%%xmm10,0x190(%%rdi)	\n\t	movaps	%%xmm11,0x310(%%rdi)	\n\t"\
		"movaps	%%xmm13,0x280(%%rdi)	\n\t	movaps	%%xmm15,0x180(%%rdi)	\n\t	movaps	%%xmm12,0x300(%%rdi)	\n\t"\
		"movaps	%%xmm5 ,0x110(%%rdi)	\n\t	movaps	%%xmm7 ,0x210(%%rdi)	\n\t	movaps	%%xmm4 ,0x090(%%rdi)	\n\t"\
	"/* SSE2_RADIX_07_DFT(01,05,09,13,17,21,25 -> 21,01,09,17,25,05,13): */	\n\t"\
		"addq	$0x020,%%rdi			\n\t	/*** Imaginary Parts: ***/			\n\t"\
		"movaps	0x080(%%rdi),%%xmm6		\n\t	movaps	0x090(%%rdi),%%xmm14		\n\t"\
		"movaps	0x300(%%rdi),%%xmm1		\n\t	movaps	0x310(%%rdi),%%xmm9 		\n\t"\
		"movaps	0x100(%%rdi),%%xmm5		\n\t	movaps	0x110(%%rdi),%%xmm13		\n\t"\
		"movaps	0x280(%%rdi),%%xmm2		\n\t	movaps	0x290(%%rdi),%%xmm10		\n\t"\
		"movaps	0x180(%%rdi),%%xmm4		\n\t	movaps	0x190(%%rdi),%%xmm12		\n\t"\
		"movaps	0x200(%%rdi),%%xmm3		\n\t	movaps	0x210(%%rdi),%%xmm11		\n\t"\
		"subpd	%%xmm1,%%xmm6			\n\t	subpd	%%xmm9 ,%%xmm14				\n\t"\
		"addpd	%%xmm1,%%xmm1			\n\t	addpd	%%xmm9  ,%%xmm9 			\n\t"\
		"addpd	%%xmm6,%%xmm1			\n\t	addpd	%%xmm14,%%xmm9  			\n\t"\
		"subpd	%%xmm2,%%xmm5			\n\t	subpd	%%xmm10,%%xmm13				\n\t"\
		"addpd	%%xmm2,%%xmm2			\n\t	addpd	%%xmm10,%%xmm10				\n\t"\
		"addpd	%%xmm5,%%xmm2			\n\t	addpd	%%xmm13,%%xmm10				\n\t"\
		"movaps	     (%%rdi),%%xmm0		\n\t	movaps	0x010(%%rdi),%%xmm8 		\n\t"\
		"subpd	%%xmm3,%%xmm4			\n\t	subpd	%%xmm11,%%xmm12				\n\t"\
		"addpd	%%xmm3,%%xmm3			\n\t	addpd	%%xmm11,%%xmm11				\n\t"\
		"addpd	%%xmm4,%%xmm3			\n\t	addpd	%%xmm12,%%xmm11				\n\t"\
		"/* Only alloc 2 scratch storage slots in calling code, so 16-reg version use 2 of the eventual outputs for Im parts */\n\t"\
		"movaps	%%xmm0,     (%%rsi)		\n\t	movaps	%%xmm8 ,0x080(%%rdi)		\n\t"\
		"movaps	%%xmm6,0x010(%%rsi)		\n\t	movaps	%%xmm14,0x100(%%rdi)		\n\t"\
		"addpd	      %%xmm1,%%xmm0		\n\t	addpd	      %%xmm9 ,%%xmm8  		\n\t"\
		"movaps	%%xmm5,%%xmm7			\n\t	movaps	%%xmm13,%%xmm15				\n\t"\
		"addpd	      %%xmm2,%%xmm3		\n\t	addpd	      %%xmm10,%%xmm11		\n\t"\
		"subpd	      %%xmm4,%%xmm5		\n\t	subpd	      %%xmm12,%%xmm13		\n\t"\
		"subpd	      %%xmm2,%%xmm1		\n\t	subpd	      %%xmm10,%%xmm9  		\n\t"\
		"subpd	      %%xmm7,%%xmm6		\n\t	subpd	      %%xmm15,%%xmm14		\n\t"\
		"addpd	      %%xmm2,%%xmm2		\n\t	addpd	      %%xmm10,%%xmm10		\n\t"\
		"addpd	      %%xmm7,%%xmm4		\n\t	addpd	      %%xmm15,%%xmm12		\n\t"\
		"addpd	      %%xmm3,%%xmm0		\n\t	addpd	      %%xmm11,%%xmm8  		\n\t"\
		"addpd	0x010(%%rsi),%%xmm5		\n\t	addpd	0x100(%%rdi),%%xmm13		\n\t"\
		"subpd	      %%xmm2,%%xmm3		\n\t	subpd	      %%xmm10,%%xmm11		\n\t"\
		"movaps	%%xmm4,%%xmm7			\n\t	movaps	%%xmm12,%%xmm15				\n\t"\
		"movaps	%%xmm0, 0x280(%%rdi)	\n\t	movaps	%%xmm8 ,0x290(%%rdi)		\n\t"\
		"subpd	      %%xmm6,%%xmm4		\n\t	subpd	      %%xmm14,%%xmm12		\n\t"\
		"movaps	%%xmm1,%%xmm2			\n\t	movaps	%%xmm9 ,%%xmm10				\n\t"\
		"subpd	     (%%rsi),%%xmm0		\n\t	subpd	0x080(%%rdi),%%xmm8  		\n\t"\
		"mulpd	0x010(%%rax),%%xmm5		\n\t	mulpd	0x010(%%rax),%%xmm13		\n\t"\
		"addpd	      %%xmm3,%%xmm2		\n\t	addpd	      %%xmm11,%%xmm10		\n\t"\
		"mulpd	     (%%rcx),%%xmm3		\n\t	mulpd	     (%%rcx),%%xmm11		\n\t"\
		"mulpd	0x010(%%rdx),%%xmm4		\n\t	mulpd	0x010(%%rdx),%%xmm12		\n\t"\
		"mulpd	     (%%rbx),%%xmm1		\n\t	mulpd	     (%%rbx),%%xmm9  		\n\t"\
		"mulpd	0x010(%%rbx),%%xmm6		\n\t	mulpd	0x010(%%rbx),%%xmm14		\n\t"\
		"mulpd	     (%%rax),%%xmm0		\n\t	mulpd	     (%%rax),%%xmm8  		\n\t"\
		"mulpd	0x010(%%rcx),%%xmm7		\n\t	mulpd	0x010(%%rcx),%%xmm15		\n\t"\
		"mulpd	     (%%rdx),%%xmm2		\n\t	mulpd	     (%%rdx),%%xmm10		\n\t"\
		"addpd	0x280(%%rdi),%%xmm0		\n\t	addpd	0x290(%%rdi),%%xmm8  		\n\t"\
		"addpd	      %%xmm4,%%xmm6		\n\t	addpd	      %%xmm12,%%xmm14		\n\t"\
		"subpd	      %%xmm2,%%xmm1		\n\t	subpd	      %%xmm10,%%xmm9  		\n\t"\
		"subpd	      %%xmm7,%%xmm4		\n\t	subpd	      %%xmm15,%%xmm12		\n\t"\
		"subpd	      %%xmm2,%%xmm3		\n\t	subpd	      %%xmm10,%%xmm11		\n\t"\
		"movaps	%%xmm0,%%xmm2			\n\t	movaps	%%xmm8 ,%%xmm10				\n\t"\
		"movaps	%%xmm5,%%xmm7			\n\t	movaps	%%xmm13,%%xmm15				\n\t"\
		"addpd	      %%xmm1,%%xmm0		\n\t	addpd	      %%xmm9 ,%%xmm8  		\n\t"\
		"addpd	      %%xmm6,%%xmm5		\n\t	addpd	      %%xmm14,%%xmm13		\n\t"\
		"addpd	      %%xmm3,%%xmm1		\n\t	addpd	      %%xmm11,%%xmm9  		\n\t"\
		"addpd	      %%xmm4,%%xmm6		\n\t	addpd	      %%xmm12,%%xmm14		\n\t"\
		"addpd	      %%xmm2,%%xmm3		\n\t	addpd	      %%xmm10,%%xmm11		\n\t"\
		"addpd	      %%xmm7,%%xmm4		\n\t	addpd	      %%xmm15,%%xmm12		\n\t"\
		"subpd	      %%xmm1,%%xmm2		\n\t	subpd	      %%xmm9 ,%%xmm10		\n\t"\
		"subpd	      %%xmm6,%%xmm7		\n\t	subpd	      %%xmm14,%%xmm15		\n\t"\
		"/* xmm1,6,9,14 free ... Note the order reversal on the 3rd pair of outputs: */\n\t"\
		"subpd	%%xmm13,%%xmm0			\n\t	subpd	%%xmm15,%%xmm2			\n\t	subpd	%%xmm12,%%xmm3 			\n\t"\
		"subpd	%%xmm5 ,%%xmm8  		\n\t	subpd	%%xmm7 ,%%xmm10			\n\t	subpd	%%xmm4 ,%%xmm11			\n\t"\
		"addpd	%%xmm13,%%xmm13			\n\t	addpd	%%xmm15,%%xmm15			\n\t	addpd	%%xmm12,%%xmm12			\n\t"\
		"addpd	%%xmm5 ,%%xmm5			\n\t	addpd	%%xmm7 ,%%xmm7			\n\t	addpd	%%xmm4 ,%%xmm4 			\n\t"\
		"addpd	%%xmm0 ,%%xmm13			\n\t	addpd	%%xmm2 ,%%xmm15			\n\t	addpd	%%xmm3 ,%%xmm12			\n\t"\
		"addpd	%%xmm8 ,%%xmm5			\n\t	addpd	%%xmm10,%%xmm7			\n\t	addpd	%%xmm11,%%xmm4 			\n\t"\
		"movaps	%%xmm0 ,     (%%rdi)	\n\t	movaps	%%xmm2 ,0x100(%%rdi)	\n\t	movaps	%%xmm3 ,0x300(%%rdi)	\n\t"\
		"movaps	%%xmm8 ,0x190(%%rdi)	\n\t	movaps	%%xmm10,0x090(%%rdi)	\n\t	movaps	%%xmm11,0x210(%%rdi)	\n\t"\
		"movaps	%%xmm13,0x180(%%rdi)	\n\t	movaps	%%xmm15,0x080(%%rdi)	\n\t	movaps	%%xmm12,0x200(%%rdi)	\n\t"\
		"movaps	%%xmm5 ,0x010(%%rdi)	\n\t	movaps	%%xmm7 ,0x110(%%rdi)	\n\t	movaps	%%xmm4 ,0x310(%%rdi)	\n\t"\
	"/* SSE2_RADIX_07_DFT(02,06,10,14,18,22,26 -> 14,22,02,10,18,26,06): */	\n\t"\
		"addq	$0x020,%%rdi			\n\t	/*** Imaginary Parts: ***/			\n\t"\
		"movaps	0x080(%%rdi),%%xmm6		\n\t	movaps	0x090(%%rdi),%%xmm14		\n\t"\
		"movaps	0x300(%%rdi),%%xmm1		\n\t	movaps	0x310(%%rdi),%%xmm9 		\n\t"\
		"movaps	0x100(%%rdi),%%xmm5		\n\t	movaps	0x110(%%rdi),%%xmm13		\n\t"\
		"movaps	0x280(%%rdi),%%xmm2		\n\t	movaps	0x290(%%rdi),%%xmm10		\n\t"\
		"movaps	0x180(%%rdi),%%xmm4		\n\t	movaps	0x190(%%rdi),%%xmm12		\n\t"\
		"movaps	0x200(%%rdi),%%xmm3		\n\t	movaps	0x210(%%rdi),%%xmm11		\n\t"\
		"subpd	%%xmm1,%%xmm6			\n\t	subpd	%%xmm9 ,%%xmm14				\n\t"\
		"addpd	%%xmm1,%%xmm1			\n\t	addpd	%%xmm9  ,%%xmm9 			\n\t"\
		"addpd	%%xmm6,%%xmm1			\n\t	addpd	%%xmm14,%%xmm9  			\n\t"\
		"subpd	%%xmm2,%%xmm5			\n\t	subpd	%%xmm10,%%xmm13				\n\t"\
		"addpd	%%xmm2,%%xmm2			\n\t	addpd	%%xmm10,%%xmm10				\n\t"\
		"addpd	%%xmm5,%%xmm2			\n\t	addpd	%%xmm13,%%xmm10				\n\t"\
		"movaps	     (%%rdi),%%xmm0		\n\t	movaps	0x010(%%rdi),%%xmm8 		\n\t"\
		"subpd	%%xmm3,%%xmm4			\n\t	subpd	%%xmm11,%%xmm12				\n\t"\
		"addpd	%%xmm3,%%xmm3			\n\t	addpd	%%xmm11,%%xmm11				\n\t"\
		"addpd	%%xmm4,%%xmm3			\n\t	addpd	%%xmm12,%%xmm11				\n\t"\
		"/* Only alloc 2 scratch storage slots in calling code, so 16-reg version use 2 of the eventual outputs for Im parts */\n\t"\
		"movaps	%%xmm0,     (%%rsi)		\n\t	movaps	%%xmm8 ,0x080(%%rdi)		\n\t"\
		"movaps	%%xmm6,0x010(%%rsi)		\n\t	movaps	%%xmm14,0x100(%%rdi)		\n\t"\
		"addpd	      %%xmm1,%%xmm0		\n\t	addpd	      %%xmm9 ,%%xmm8  		\n\t"\
		"movaps	%%xmm5,%%xmm7			\n\t	movaps	%%xmm13,%%xmm15				\n\t"\
		"addpd	      %%xmm2,%%xmm3		\n\t	addpd	      %%xmm10,%%xmm11		\n\t"\
		"subpd	      %%xmm4,%%xmm5		\n\t	subpd	      %%xmm12,%%xmm13		\n\t"\
		"subpd	      %%xmm2,%%xmm1		\n\t	subpd	      %%xmm10,%%xmm9  		\n\t"\
		"subpd	      %%xmm7,%%xmm6		\n\t	subpd	      %%xmm15,%%xmm14		\n\t"\
		"addpd	      %%xmm2,%%xmm2		\n\t	addpd	      %%xmm10,%%xmm10		\n\t"\
		"addpd	      %%xmm7,%%xmm4		\n\t	addpd	      %%xmm15,%%xmm12		\n\t"\
		"addpd	      %%xmm3,%%xmm0		\n\t	addpd	      %%xmm11,%%xmm8  		\n\t"\
		"addpd	0x010(%%rsi),%%xmm5		\n\t	addpd	0x100(%%rdi),%%xmm13		\n\t"\
		"subpd	      %%xmm2,%%xmm3		\n\t	subpd	      %%xmm10,%%xmm11		\n\t"\
		"movaps	%%xmm4,%%xmm7			\n\t	movaps	%%xmm12,%%xmm15				\n\t"\
		"movaps	%%xmm0, 0x180(%%rdi)	\n\t	movaps	%%xmm8 ,0x190(%%rdi)		\n\t"\
		"subpd	      %%xmm6,%%xmm4		\n\t	subpd	      %%xmm14,%%xmm12		\n\t"\
		"movaps	%%xmm1,%%xmm2			\n\t	movaps	%%xmm9 ,%%xmm10				\n\t"\
		"subpd	     (%%rsi),%%xmm0		\n\t	subpd	0x080(%%rdi),%%xmm8  		\n\t"\
		"mulpd	0x010(%%rax),%%xmm5		\n\t	mulpd	0x010(%%rax),%%xmm13		\n\t"\
		"addpd	      %%xmm3,%%xmm2		\n\t	addpd	      %%xmm11,%%xmm10		\n\t"\
		"mulpd	     (%%rcx),%%xmm3		\n\t	mulpd	     (%%rcx),%%xmm11		\n\t"\
		"mulpd	0x010(%%rdx),%%xmm4		\n\t	mulpd	0x010(%%rdx),%%xmm12		\n\t"\
		"mulpd	     (%%rbx),%%xmm1		\n\t	mulpd	     (%%rbx),%%xmm9  		\n\t"\
		"mulpd	0x010(%%rbx),%%xmm6		\n\t	mulpd	0x010(%%rbx),%%xmm14		\n\t"\
		"mulpd	     (%%rax),%%xmm0		\n\t	mulpd	     (%%rax),%%xmm8  		\n\t"\
		"mulpd	0x010(%%rcx),%%xmm7		\n\t	mulpd	0x010(%%rcx),%%xmm15		\n\t"\
		"mulpd	     (%%rdx),%%xmm2		\n\t	mulpd	     (%%rdx),%%xmm10		\n\t"\
		"addpd	0x180(%%rdi),%%xmm0		\n\t	addpd	0x190(%%rdi),%%xmm8  		\n\t"\
		"addpd	      %%xmm4,%%xmm6		\n\t	addpd	      %%xmm12,%%xmm14		\n\t"\
		"subpd	      %%xmm2,%%xmm1		\n\t	subpd	      %%xmm10,%%xmm9  		\n\t"\
		"subpd	      %%xmm7,%%xmm4		\n\t	subpd	      %%xmm15,%%xmm12		\n\t"\
		"subpd	      %%xmm2,%%xmm3		\n\t	subpd	      %%xmm10,%%xmm11		\n\t"\
		"movaps	%%xmm0,%%xmm2			\n\t	movaps	%%xmm8 ,%%xmm10				\n\t"\
		"movaps	%%xmm5,%%xmm7			\n\t	movaps	%%xmm13,%%xmm15				\n\t"\
		"addpd	      %%xmm1,%%xmm0		\n\t	addpd	      %%xmm9 ,%%xmm8  		\n\t"\
		"addpd	      %%xmm6,%%xmm5		\n\t	addpd	      %%xmm14,%%xmm13		\n\t"\
		"addpd	      %%xmm3,%%xmm1		\n\t	addpd	      %%xmm11,%%xmm9  		\n\t"\
		"addpd	      %%xmm4,%%xmm6		\n\t	addpd	      %%xmm12,%%xmm14		\n\t"\
		"addpd	      %%xmm2,%%xmm3		\n\t	addpd	      %%xmm10,%%xmm11		\n\t"\
		"addpd	      %%xmm7,%%xmm4		\n\t	addpd	      %%xmm15,%%xmm12		\n\t"\
		"subpd	      %%xmm1,%%xmm2		\n\t	subpd	      %%xmm9 ,%%xmm10		\n\t"\
		"subpd	      %%xmm6,%%xmm7		\n\t	subpd	      %%xmm14,%%xmm15		\n\t"\
		"/* xmm1,6,9,14 free ... Note the order reversal on the 3rd pair of outputs: */\n\t"\
		"subpd	%%xmm13,%%xmm0			\n\t	subpd	%%xmm15,%%xmm2			\n\t	subpd	%%xmm12,%%xmm3 			\n\t"\
		"subpd	%%xmm5 ,%%xmm8  		\n\t	subpd	%%xmm7 ,%%xmm10			\n\t	subpd	%%xmm4 ,%%xmm11			\n\t"\
		"addpd	%%xmm13,%%xmm13			\n\t	addpd	%%xmm15,%%xmm15			\n\t	addpd	%%xmm12,%%xmm12			\n\t"\
		"addpd	%%xmm5 ,%%xmm5			\n\t	addpd	%%xmm7 ,%%xmm7			\n\t	addpd	%%xmm4 ,%%xmm4 			\n\t"\
		"addpd	%%xmm0 ,%%xmm13			\n\t	addpd	%%xmm2 ,%%xmm15			\n\t	addpd	%%xmm3 ,%%xmm12			\n\t"\
		"addpd	%%xmm8 ,%%xmm5			\n\t	addpd	%%xmm10,%%xmm7			\n\t	addpd	%%xmm11,%%xmm4 			\n\t"\
		"movaps	%%xmm0 ,0x280(%%rdi)	\n\t	movaps	%%xmm2 ,     (%%rdi)	\n\t	movaps	%%xmm3 ,0x200(%%rdi)	\n\t"\
		"movaps	%%xmm8 ,0x090(%%rdi)	\n\t	movaps	%%xmm10,0x310(%%rdi)	\n\t	movaps	%%xmm11,0x110(%%rdi)	\n\t"\
		"movaps	%%xmm13,0x080(%%rdi)	\n\t	movaps	%%xmm15,0x300(%%rdi)	\n\t	movaps	%%xmm12,0x100(%%rdi)	\n\t"\
		"movaps	%%xmm5 ,0x290(%%rdi)	\n\t	movaps	%%xmm7 ,0x010(%%rdi)	\n\t	movaps	%%xmm4 ,0x210(%%rdi)	\n\t"\
	"/* SSE2_RADIX_07_DFT(03,07,11,15,19,23,27 -> 07,15,23,03,11,19,27): */	\n\t"\
		"addq	$0x020,%%rdi			\n\t	/*** Imaginary Parts: ***/			\n\t"\
		"movaps	0x080(%%rdi),%%xmm6		\n\t	movaps	0x090(%%rdi),%%xmm14		\n\t"\
		"movaps	0x300(%%rdi),%%xmm1		\n\t	movaps	0x310(%%rdi),%%xmm9 		\n\t"\
		"movaps	0x100(%%rdi),%%xmm5		\n\t	movaps	0x110(%%rdi),%%xmm13		\n\t"\
		"movaps	0x280(%%rdi),%%xmm2		\n\t	movaps	0x290(%%rdi),%%xmm10		\n\t"\
		"movaps	0x180(%%rdi),%%xmm4		\n\t	movaps	0x190(%%rdi),%%xmm12		\n\t"\
		"movaps	0x200(%%rdi),%%xmm3		\n\t	movaps	0x210(%%rdi),%%xmm11		\n\t"\
		"subpd	%%xmm1,%%xmm6			\n\t	subpd	%%xmm9 ,%%xmm14				\n\t"\
		"addpd	%%xmm1,%%xmm1			\n\t	addpd	%%xmm9  ,%%xmm9 			\n\t"\
		"addpd	%%xmm6,%%xmm1			\n\t	addpd	%%xmm14,%%xmm9  			\n\t"\
		"subpd	%%xmm2,%%xmm5			\n\t	subpd	%%xmm10,%%xmm13				\n\t"\
		"addpd	%%xmm2,%%xmm2			\n\t	addpd	%%xmm10,%%xmm10				\n\t"\
		"addpd	%%xmm5,%%xmm2			\n\t	addpd	%%xmm13,%%xmm10				\n\t"\
		"movaps	     (%%rdi),%%xmm0		\n\t	movaps	0x010(%%rdi),%%xmm8 		\n\t"\
		"subpd	%%xmm3,%%xmm4			\n\t	subpd	%%xmm11,%%xmm12				\n\t"\
		"addpd	%%xmm3,%%xmm3			\n\t	addpd	%%xmm11,%%xmm11				\n\t"\
		"addpd	%%xmm4,%%xmm3			\n\t	addpd	%%xmm12,%%xmm11				\n\t"\
		"/* Only alloc 2 scratch storage slots in calling code, so 16-reg version use 2 of the eventual outputs for Im parts */\n\t"\
		"movaps	%%xmm0,     (%%rsi)		\n\t	movaps	%%xmm8 ,0x180(%%rdi)		\n\t"\
		"movaps	%%xmm6,0x010(%%rsi)		\n\t	movaps	%%xmm14,0x100(%%rdi)		\n\t"\
		"addpd	      %%xmm1,%%xmm0		\n\t	addpd	      %%xmm9 ,%%xmm8  		\n\t"\
		"movaps	%%xmm5,%%xmm7			\n\t	movaps	%%xmm13,%%xmm15				\n\t"\
		"addpd	      %%xmm2,%%xmm3		\n\t	addpd	      %%xmm10,%%xmm11		\n\t"\
		"subpd	      %%xmm4,%%xmm5		\n\t	subpd	      %%xmm12,%%xmm13		\n\t"\
		"subpd	      %%xmm2,%%xmm1		\n\t	subpd	      %%xmm10,%%xmm9  		\n\t"\
		"subpd	      %%xmm7,%%xmm6		\n\t	subpd	      %%xmm15,%%xmm14		\n\t"\
		"addpd	      %%xmm2,%%xmm2		\n\t	addpd	      %%xmm10,%%xmm10		\n\t"\
		"addpd	      %%xmm7,%%xmm4		\n\t	addpd	      %%xmm15,%%xmm12		\n\t"\
		"addpd	      %%xmm3,%%xmm0		\n\t	addpd	      %%xmm11,%%xmm8  		\n\t"\
		"addpd	0x010(%%rsi),%%xmm5		\n\t	addpd	0x100(%%rdi),%%xmm13		\n\t"\
		"subpd	      %%xmm2,%%xmm3		\n\t	subpd	      %%xmm10,%%xmm11		\n\t"\
		"movaps	%%xmm4,%%xmm7			\n\t	movaps	%%xmm12,%%xmm15				\n\t"\
		"movaps	%%xmm0, 0x080(%%rdi)	\n\t	movaps	%%xmm8 ,0x090(%%rdi)		\n\t"\
		"subpd	      %%xmm6,%%xmm4		\n\t	subpd	      %%xmm14,%%xmm12		\n\t"\
		"movaps	%%xmm1,%%xmm2			\n\t	movaps	%%xmm9 ,%%xmm10				\n\t"\
		"subpd	     (%%rsi),%%xmm0		\n\t	subpd	0x180(%%rdi),%%xmm8  		\n\t"\
		"mulpd	0x010(%%rax),%%xmm5		\n\t	mulpd	0x010(%%rax),%%xmm13		\n\t"\
		"addpd	      %%xmm3,%%xmm2		\n\t	addpd	      %%xmm11,%%xmm10		\n\t"\
		"mulpd	     (%%rcx),%%xmm3		\n\t	mulpd	     (%%rcx),%%xmm11		\n\t"\
		"mulpd	0x010(%%rdx),%%xmm4		\n\t	mulpd	0x010(%%rdx),%%xmm12		\n\t"\
		"mulpd	     (%%rbx),%%xmm1		\n\t	mulpd	     (%%rbx),%%xmm9  		\n\t"\
		"mulpd	0x010(%%rbx),%%xmm6		\n\t	mulpd	0x010(%%rbx),%%xmm14		\n\t"\
		"mulpd	     (%%rax),%%xmm0		\n\t	mulpd	     (%%rax),%%xmm8  		\n\t"\
		"mulpd	0x010(%%rcx),%%xmm7		\n\t	mulpd	0x010(%%rcx),%%xmm15		\n\t"\
		"mulpd	     (%%rdx),%%xmm2		\n\t	mulpd	     (%%rdx),%%xmm10		\n\t"\
		"addpd	0x080(%%rdi),%%xmm0		\n\t	addpd	0x090(%%rdi),%%xmm8  		\n\t"\
		"addpd	      %%xmm4,%%xmm6		\n\t	addpd	      %%xmm12,%%xmm14		\n\t"\
		"subpd	      %%xmm2,%%xmm1		\n\t	subpd	      %%xmm10,%%xmm9  		\n\t"\
		"subpd	      %%xmm7,%%xmm4		\n\t	subpd	      %%xmm15,%%xmm12		\n\t"\
		"subpd	      %%xmm2,%%xmm3		\n\t	subpd	      %%xmm10,%%xmm11		\n\t"\
		"movaps	%%xmm0,%%xmm2			\n\t	movaps	%%xmm8 ,%%xmm10				\n\t"\
		"movaps	%%xmm5,%%xmm7			\n\t	movaps	%%xmm13,%%xmm15				\n\t"\
		"addpd	      %%xmm1,%%xmm0		\n\t	addpd	      %%xmm9 ,%%xmm8  		\n\t"\
		"addpd	      %%xmm6,%%xmm5		\n\t	addpd	      %%xmm14,%%xmm13		\n\t"\
		"addpd	      %%xmm3,%%xmm1		\n\t	addpd	      %%xmm11,%%xmm9  		\n\t"\
		"addpd	      %%xmm4,%%xmm6		\n\t	addpd	      %%xmm12,%%xmm14		\n\t"\
		"addpd	      %%xmm2,%%xmm3		\n\t	addpd	      %%xmm10,%%xmm11		\n\t"\
		"addpd	      %%xmm7,%%xmm4		\n\t	addpd	      %%xmm15,%%xmm12		\n\t"\
		"subpd	      %%xmm1,%%xmm2		\n\t	subpd	      %%xmm9 ,%%xmm10		\n\t"\
		"subpd	      %%xmm6,%%xmm7		\n\t	subpd	      %%xmm14,%%xmm15		\n\t"\
		"/* xmm1,6,9,14 free ... Note the order reversal on the 3rd pair of outputs: */\n\t"\
		"subpd	%%xmm13,%%xmm0			\n\t	subpd	%%xmm15,%%xmm2			\n\t	subpd	%%xmm12,%%xmm3 			\n\t"\
		"subpd	%%xmm5 ,%%xmm8  		\n\t	subpd	%%xmm7 ,%%xmm10			\n\t	subpd	%%xmm4 ,%%xmm11			\n\t"\
		"addpd	%%xmm13,%%xmm13			\n\t	addpd	%%xmm15,%%xmm15			\n\t	addpd	%%xmm12,%%xmm12			\n\t"\
		"addpd	%%xmm5 ,%%xmm5			\n\t	addpd	%%xmm7 ,%%xmm7			\n\t	addpd	%%xmm4 ,%%xmm4 			\n\t"\
		"addpd	%%xmm0 ,%%xmm13			\n\t	addpd	%%xmm2 ,%%xmm15			\n\t	addpd	%%xmm3 ,%%xmm12			\n\t"\
		"addpd	%%xmm8 ,%%xmm5			\n\t	addpd	%%xmm10,%%xmm7			\n\t	addpd	%%xmm11,%%xmm4 			\n\t"\
		"movaps	%%xmm0 ,0x180(%%rdi)	\n\t	movaps	%%xmm2 ,0x280(%%rdi)	\n\t	movaps	%%xmm3 ,0x100(%%rdi)	\n\t"\
		"movaps	%%xmm8 ,0x310(%%rdi)	\n\t	movaps	%%xmm10,0x210(%%rdi)	\n\t	movaps	%%xmm11,0x010(%%rdi)	\n\t"\
		"movaps	%%xmm13,0x300(%%rdi)	\n\t	movaps	%%xmm15,0x200(%%rdi)	\n\t	movaps	%%xmm12,     (%%rdi)	\n\t"\
		"movaps	%%xmm5 ,0x190(%%rdi)	\n\t	movaps	%%xmm7 ,0x290(%%rdi)	\n\t	movaps	%%xmm4 ,0x110(%%rdi)	\n\t"\
		:					/* outputs: none */\
		: [__add] "m" (Xadd)	/* All inputs from memory addresses here */\
		 ,[__p01] "m" (Xp01)\
		 ,[__p02] "m" (Xp02)\
		 ,[__p03] "m" (Xp03)\
		 ,[__p04] "m" (Xp04)\
		 ,[__p08] "m" (Xp08)\
		 ,[__p12] "m" (Xp12)\
		 ,[__p16] "m" (Xp16)\
		 ,[__p20] "m" (Xp20)\
		 ,[__p24] "m" (Xp24)\
		 ,[__out] "m" (Xout)\
		 ,[__cc0] "m" (Xcc0)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"		/* Clobbered registers */\
	);\
	}

	#define	SSE2_RADIX28_DIF_NOTWIDDLE(Xadd,Xp01,Xp02,Xp03,Xp04,Xp08,Xp12,Xp16,Xp20,Xp24,Xout,Xcc0)\
	{\
	__asm__ volatile (\
	"/* SSE2_RADIX_07_DFT(00,24,20,16,12,08,04r -> 00,04,08,12,16,20,24): */\n\t"\
		"movq	%[__out],%%rdi			\n\t"\
		"movq	%[__cc0],%%rax			\n\t"\
		"movq	$0x020,%%rbx			\n\t"\
		"movq	$0x040,%%rcx			\n\t"\
		"movq	$0x060,%%rdx			\n\t"\
		"movq	$0x080,%%rsi			\n\t"\
		"addq	%%rax,%%rbx				\n\t"\
		"addq	%%rax,%%rcx				\n\t"\
		"addq	%%rax,%%rdx				\n\t"\
		"addq	%%rax,%%rsi				\n\t	/*** Imaginary Parts: ***/			\n\t"\
		"movaps	0x300(%%rdi),%%xmm6		\n\t	movaps	0x310(%%rdi),%%xmm14		\n\t"\
		"movaps	0x080(%%rdi),%%xmm1		\n\t	movaps	0x090(%%rdi),%%xmm9 		\n\t"\
		"movaps	0x280(%%rdi),%%xmm5		\n\t	movaps	0x290(%%rdi),%%xmm13		\n\t"\
		"movaps	0x100(%%rdi),%%xmm2		\n\t	movaps	0x110(%%rdi),%%xmm10		\n\t"\
		"movaps	0x200(%%rdi),%%xmm4		\n\t	movaps	0x210(%%rdi),%%xmm12		\n\t"\
		"movaps	0x180(%%rdi),%%xmm3		\n\t	movaps	0x190(%%rdi),%%xmm11		\n\t"\
		"subpd	%%xmm1,%%xmm6			\n\t	subpd	%%xmm9 ,%%xmm14				\n\t"\
		"addpd	%%xmm1,%%xmm1			\n\t	addpd	%%xmm9  ,%%xmm9 			\n\t"\
		"addpd	%%xmm6,%%xmm1			\n\t	addpd	%%xmm14,%%xmm9  			\n\t"\
		"subpd	%%xmm2,%%xmm5			\n\t	subpd	%%xmm10,%%xmm13				\n\t"\
		"addpd	%%xmm2,%%xmm2			\n\t	addpd	%%xmm10,%%xmm10				\n\t"\
		"addpd	%%xmm5,%%xmm2			\n\t	addpd	%%xmm13,%%xmm10				\n\t"\
		"movaps	     (%%rdi),%%xmm0		\n\t	movaps	0x010(%%rdi),%%xmm8 		\n\t"\
		"subpd	%%xmm3,%%xmm4			\n\t	subpd	%%xmm11,%%xmm12				\n\t"\
		"addpd	%%xmm3,%%xmm3			\n\t	addpd	%%xmm11,%%xmm11				\n\t"\
		"addpd	%%xmm4,%%xmm3			\n\t	addpd	%%xmm12,%%xmm11				\n\t"\
		"/* Only alloc 2 scratch storageslots in calling code, so 16-reg version use 2 of the eventual outputs for Im parts */\n\t"\
		"movaps	%%xmm0,     (%%rsi)		\n\t	movaps	%%xmm8 ,0x080(%%rdi)		\n\t"\
		"movaps	%%xmm6,0x010(%%rsi)		\n\t	movaps	%%xmm14,0x100(%%rdi)		\n\t"\
		"addpd	      %%xmm1,%%xmm0		\n\t	addpd	      %%xmm9 ,%%xmm8  		\n\t"\
		"movaps	%%xmm5,%%xmm7			\n\t	movaps	%%xmm13,%%xmm15				\n\t"\
		"addpd	      %%xmm2,%%xmm3		\n\t	addpd	      %%xmm10,%%xmm11		\n\t"\
		"subpd	      %%xmm4,%%xmm5		\n\t	subpd	      %%xmm12,%%xmm13		\n\t"\
		"subpd	      %%xmm2,%%xmm1		\n\t	subpd	      %%xmm10,%%xmm9  		\n\t"\
		"subpd	      %%xmm7,%%xmm6		\n\t	subpd	      %%xmm15,%%xmm14		\n\t"\
		"addpd	      %%xmm2,%%xmm2		\n\t	addpd	      %%xmm10,%%xmm10		\n\t"\
		"addpd	      %%xmm7,%%xmm4		\n\t	addpd	      %%xmm15,%%xmm12		\n\t"\
		"addpd	      %%xmm3,%%xmm0		\n\t	addpd	      %%xmm11,%%xmm8  		\n\t"\
		"addpd	0x010(%%rsi),%%xmm5		\n\t	addpd	0x100(%%rdi),%%xmm13		\n\t"\
		"subpd	      %%xmm2,%%xmm3		\n\t	subpd	      %%xmm10,%%xmm11		\n\t"\
		"movaps	%%xmm4,%%xmm7			\n\t	movaps	%%xmm12,%%xmm15				\n\t"\
		"movaps	%%xmm0,      (%%rdi)	\n\t	movaps	%%xmm8 ,0x010(%%rdi)		\n\t"\
		"subpd	      %%xmm6,%%xmm4		\n\t	subpd	      %%xmm14,%%xmm12		\n\t"\
		"movaps	%%xmm1,%%xmm2			\n\t	movaps	%%xmm9 ,%%xmm10				\n\t"\
		"subpd	     (%%rsi),%%xmm0		\n\t	subpd	0x080(%%rdi),%%xmm8  		\n\t"\
		"mulpd	0x010(%%rax),%%xmm5		\n\t	mulpd	0x010(%%rax),%%xmm13		\n\t"\
		"addpd	      %%xmm3,%%xmm2		\n\t	addpd	      %%xmm11,%%xmm10		\n\t"\
		"mulpd	     (%%rcx),%%xmm3		\n\t	mulpd	     (%%rcx),%%xmm11		\n\t"\
		"mulpd	0x010(%%rdx),%%xmm4		\n\t	mulpd	0x010(%%rdx),%%xmm12		\n\t"\
		"mulpd	     (%%rbx),%%xmm1		\n\t	mulpd	     (%%rbx),%%xmm9  		\n\t"\
		"mulpd	0x010(%%rbx),%%xmm6		\n\t	mulpd	0x010(%%rbx),%%xmm14		\n\t"\
		"mulpd	     (%%rax),%%xmm0		\n\t	mulpd	     (%%rax),%%xmm8  		\n\t"\
		"mulpd	0x010(%%rcx),%%xmm7		\n\t	mulpd	0x010(%%rcx),%%xmm15		\n\t"\
		"mulpd	     (%%rdx),%%xmm2		\n\t	mulpd	     (%%rdx),%%xmm10		\n\t"\
		"addpd	     (%%rdi),%%xmm0		\n\t	addpd	0x010(%%rdi),%%xmm8  		\n\t"\
		"addpd	      %%xmm4,%%xmm6		\n\t	addpd	      %%xmm12,%%xmm14		\n\t"\
		"subpd	      %%xmm2,%%xmm1		\n\t	subpd	      %%xmm10,%%xmm9  		\n\t"\
		"subpd	      %%xmm7,%%xmm4		\n\t	subpd	      %%xmm15,%%xmm12		\n\t"\
		"subpd	      %%xmm2,%%xmm3		\n\t	subpd	      %%xmm10,%%xmm11		\n\t"\
		"movaps	%%xmm0,%%xmm2			\n\t	movaps	%%xmm8 ,%%xmm10				\n\t"\
		"movaps	%%xmm5,%%xmm7			\n\t	movaps	%%xmm13,%%xmm15				\n\t"\
		"addpd	      %%xmm1,%%xmm0		\n\t	addpd	      %%xmm9 ,%%xmm8  		\n\t"\
		"addpd	      %%xmm6,%%xmm5		\n\t	addpd	      %%xmm14,%%xmm13		\n\t"\
		"addpd	      %%xmm3,%%xmm1		\n\t	addpd	      %%xmm11,%%xmm9  		\n\t"\
		"addpd	      %%xmm4,%%xmm6		\n\t	addpd	      %%xmm12,%%xmm14		\n\t"\
		"addpd	      %%xmm2,%%xmm3		\n\t	addpd	      %%xmm10,%%xmm11		\n\t"\
		"addpd	      %%xmm7,%%xmm4		\n\t	addpd	      %%xmm15,%%xmm12		\n\t"\
		"subpd	      %%xmm1,%%xmm2		\n\t	subpd	      %%xmm9 ,%%xmm10		\n\t"\
		"subpd	      %%xmm6,%%xmm7		\n\t	subpd	      %%xmm14,%%xmm15		\n\t"\
		"/* xmm1,6,9,14 free ... Note the order reversal on the 3rd pair of outputs: */\n\t"\
		"subpd	%%xmm13,%%xmm0			\n\t	subpd	%%xmm15,%%xmm2			\n\t	subpd	%%xmm12,%%xmm3 			\n\t"\
		"subpd	%%xmm5 ,%%xmm8  		\n\t	subpd	%%xmm7 ,%%xmm10			\n\t	subpd	%%xmm4 ,%%xmm11			\n\t"\
		"addpd	%%xmm13,%%xmm13			\n\t	addpd	%%xmm15,%%xmm15			\n\t	addpd	%%xmm12,%%xmm12			\n\t"\
		"addpd	%%xmm5 ,%%xmm5			\n\t	addpd	%%xmm7 ,%%xmm7			\n\t	addpd	%%xmm4 ,%%xmm4 			\n\t"\
		"addpd	%%xmm0 ,%%xmm13			\n\t	addpd	%%xmm2 ,%%xmm15			\n\t	addpd	%%xmm3 ,%%xmm12			\n\t"\
		"addpd	%%xmm8 ,%%xmm5			\n\t	addpd	%%xmm10,%%xmm7			\n\t	addpd	%%xmm11,%%xmm4 			\n\t"\
		"movaps	%%xmm0 ,0x080(%%rdi)	\n\t	movaps	%%xmm2 ,0x100(%%rdi)	\n\t	movaps	%%xmm3 ,0x200(%%rdi)	\n\t"\
		"movaps	%%xmm8 ,0x310(%%rdi)	\n\t	movaps	%%xmm10,0x290(%%rdi)	\n\t	movaps	%%xmm11,0x190(%%rdi)	\n\t"\
		"movaps	%%xmm13,0x300(%%rdi)	\n\t	movaps	%%xmm15,0x280(%%rdi)	\n\t	movaps	%%xmm12,0x180(%%rdi)	\n\t"\
		"movaps	%%xmm5 ,0x090(%%rdi)	\n\t	movaps	%%xmm7 ,0x110(%%rdi)	\n\t	movaps	%%xmm4 ,0x210(%%rdi)	\n\t"\
	"/* SSE2_RADIX_07_DFT(21,17,13,09,05,01,25 -> 01,05,09,13,17,21,25): */\n\t"\
		"addq	$0x20,%%rdi	/* s1p01r */\n\t	/*** Imaginary Parts: ***/			\n\t"\
		"movaps	0x200(%%rdi),%%xmm6		\n\t	movaps	0x210(%%rdi),%%xmm14		\n\t"\
		"movaps	0x300(%%rdi),%%xmm1		\n\t	movaps	0x310(%%rdi),%%xmm9 		\n\t"\
		"movaps	0x180(%%rdi),%%xmm5		\n\t	movaps	0x190(%%rdi),%%xmm13		\n\t"\
		"movaps	     (%%rdi),%%xmm2		\n\t	movaps	0x010(%%rdi),%%xmm10		\n\t"\
		"movaps	0x100(%%rdi),%%xmm4		\n\t	movaps	0x110(%%rdi),%%xmm12		\n\t"\
		"movaps	0x080(%%rdi),%%xmm3		\n\t	movaps	0x090(%%rdi),%%xmm11		\n\t"\
		"subpd	%%xmm1,%%xmm6			\n\t	subpd	%%xmm9 ,%%xmm14				\n\t"\
		"addpd	%%xmm1,%%xmm1			\n\t	addpd	%%xmm9  ,%%xmm9 			\n\t"\
		"addpd	%%xmm6,%%xmm1			\n\t	addpd	%%xmm14,%%xmm9  			\n\t"\
		"subpd	%%xmm2,%%xmm5			\n\t	subpd	%%xmm10,%%xmm13				\n\t"\
		"addpd	%%xmm2,%%xmm2			\n\t	addpd	%%xmm10,%%xmm10				\n\t"\
		"addpd	%%xmm5,%%xmm2			\n\t	addpd	%%xmm13,%%xmm10				\n\t"\
		"movaps	0x280(%%rdi),%%xmm0		\n\t	movaps	0x290(%%rdi),%%xmm8 		\n\t"\
		"subpd	%%xmm3,%%xmm4			\n\t	subpd	%%xmm11,%%xmm12				\n\t"\
		"addpd	%%xmm3,%%xmm3			\n\t	addpd	%%xmm11,%%xmm11				\n\t"\
		"addpd	%%xmm4,%%xmm3			\n\t	addpd	%%xmm12,%%xmm11				\n\t"\
		"/* Only alloc 2 scratch storageslots in calling code, so 16-reg version use 2 of the eventual outputs for Im parts */\n\t"\
		"movaps	%%xmm0,     (%%rsi)		\n\t	movaps	%%xmm8 ,0x080(%%rdi)		\n\t"\
		"movaps	%%xmm6,0x010(%%rsi)		\n\t	movaps	%%xmm14,0x100(%%rdi)		\n\t"\
		"addpd	      %%xmm1,%%xmm0		\n\t	addpd	      %%xmm9 ,%%xmm8  		\n\t"\
		"movaps	%%xmm5,%%xmm7			\n\t	movaps	%%xmm13,%%xmm15				\n\t"\
		"addpd	      %%xmm2,%%xmm3		\n\t	addpd	      %%xmm10,%%xmm11		\n\t"\
		"subpd	      %%xmm4,%%xmm5		\n\t	subpd	      %%xmm12,%%xmm13		\n\t"\
		"subpd	      %%xmm2,%%xmm1		\n\t	subpd	      %%xmm10,%%xmm9  		\n\t"\
		"subpd	      %%xmm7,%%xmm6		\n\t	subpd	      %%xmm15,%%xmm14		\n\t"\
		"addpd	      %%xmm2,%%xmm2		\n\t	addpd	      %%xmm10,%%xmm10		\n\t"\
		"addpd	      %%xmm7,%%xmm4		\n\t	addpd	      %%xmm15,%%xmm12		\n\t"\
		"addpd	      %%xmm3,%%xmm0		\n\t	addpd	      %%xmm11,%%xmm8  		\n\t"\
		"addpd	0x010(%%rsi),%%xmm5		\n\t	addpd	0x100(%%rdi),%%xmm13		\n\t"\
		"subpd	      %%xmm2,%%xmm3		\n\t	subpd	      %%xmm10,%%xmm11		\n\t"\
		"movaps	%%xmm4,%%xmm7			\n\t	movaps	%%xmm12,%%xmm15				\n\t"\
		"movaps	%%xmm0,      (%%rdi)	\n\t	movaps	%%xmm8 ,0x010(%%rdi)		\n\t"\
		"subpd	      %%xmm6,%%xmm4		\n\t	subpd	      %%xmm14,%%xmm12		\n\t"\
		"movaps	%%xmm1,%%xmm2			\n\t	movaps	%%xmm9 ,%%xmm10				\n\t"\
		"subpd	     (%%rsi),%%xmm0		\n\t	subpd	0x080(%%rdi),%%xmm8  		\n\t"\
		"mulpd	0x010(%%rax),%%xmm5		\n\t	mulpd	0x010(%%rax),%%xmm13		\n\t"\
		"addpd	      %%xmm3,%%xmm2		\n\t	addpd	      %%xmm11,%%xmm10		\n\t"\
		"mulpd	     (%%rcx),%%xmm3		\n\t	mulpd	     (%%rcx),%%xmm11		\n\t"\
		"mulpd	0x010(%%rdx),%%xmm4		\n\t	mulpd	0x010(%%rdx),%%xmm12		\n\t"\
		"mulpd	     (%%rbx),%%xmm1		\n\t	mulpd	     (%%rbx),%%xmm9  		\n\t"\
		"mulpd	0x010(%%rbx),%%xmm6		\n\t	mulpd	0x010(%%rbx),%%xmm14		\n\t"\
		"mulpd	     (%%rax),%%xmm0		\n\t	mulpd	     (%%rax),%%xmm8  		\n\t"\
		"mulpd	0x010(%%rcx),%%xmm7		\n\t	mulpd	0x010(%%rcx),%%xmm15		\n\t"\
		"mulpd	     (%%rdx),%%xmm2		\n\t	mulpd	     (%%rdx),%%xmm10		\n\t"\
		"addpd	     (%%rdi),%%xmm0		\n\t	addpd	0x010(%%rdi),%%xmm8  		\n\t"\
		"addpd	      %%xmm4,%%xmm6		\n\t	addpd	      %%xmm12,%%xmm14		\n\t"\
		"subpd	      %%xmm2,%%xmm1		\n\t	subpd	      %%xmm10,%%xmm9  		\n\t"\
		"subpd	      %%xmm7,%%xmm4		\n\t	subpd	      %%xmm15,%%xmm12		\n\t"\
		"subpd	      %%xmm2,%%xmm3		\n\t	subpd	      %%xmm10,%%xmm11		\n\t"\
		"movaps	%%xmm0,%%xmm2			\n\t	movaps	%%xmm8 ,%%xmm10				\n\t"\
		"movaps	%%xmm5,%%xmm7			\n\t	movaps	%%xmm13,%%xmm15				\n\t"\
		"addpd	      %%xmm1,%%xmm0		\n\t	addpd	      %%xmm9 ,%%xmm8  		\n\t"\
		"addpd	      %%xmm6,%%xmm5		\n\t	addpd	      %%xmm14,%%xmm13		\n\t"\
		"addpd	      %%xmm3,%%xmm1		\n\t	addpd	      %%xmm11,%%xmm9  		\n\t"\
		"addpd	      %%xmm4,%%xmm6		\n\t	addpd	      %%xmm12,%%xmm14		\n\t"\
		"addpd	      %%xmm2,%%xmm3		\n\t	addpd	      %%xmm10,%%xmm11		\n\t"\
		"addpd	      %%xmm7,%%xmm4		\n\t	addpd	      %%xmm15,%%xmm12		\n\t"\
		"subpd	      %%xmm1,%%xmm2		\n\t	subpd	      %%xmm9 ,%%xmm10		\n\t"\
		"subpd	      %%xmm6,%%xmm7		\n\t	subpd	      %%xmm14,%%xmm15		\n\t"\
		"/* xmm1,6,9,14 free ... Note the order reversal on the 3rd pair of outputs: */\n\t"\
		"subpd	%%xmm13,%%xmm0			\n\t	subpd	%%xmm15,%%xmm2			\n\t	subpd	%%xmm12,%%xmm3 			\n\t"\
		"subpd	%%xmm5 ,%%xmm8  		\n\t	subpd	%%xmm7 ,%%xmm10			\n\t	subpd	%%xmm4 ,%%xmm11			\n\t"\
		"addpd	%%xmm13,%%xmm13			\n\t	addpd	%%xmm15,%%xmm15			\n\t	addpd	%%xmm12,%%xmm12			\n\t"\
		"addpd	%%xmm5 ,%%xmm5			\n\t	addpd	%%xmm7 ,%%xmm7			\n\t	addpd	%%xmm4 ,%%xmm4 			\n\t"\
		"addpd	%%xmm0 ,%%xmm13			\n\t	addpd	%%xmm2 ,%%xmm15			\n\t	addpd	%%xmm3 ,%%xmm12			\n\t"\
		"addpd	%%xmm8 ,%%xmm5			\n\t	addpd	%%xmm10,%%xmm7			\n\t	addpd	%%xmm11,%%xmm4 			\n\t"\
		"movaps	%%xmm0 ,0x080(%%rdi)	\n\t	movaps	%%xmm2 ,0x100(%%rdi)	\n\t	movaps	%%xmm3 ,0x200(%%rdi)	\n\t"\
		"movaps	%%xmm8 ,0x310(%%rdi)	\n\t	movaps	%%xmm10,0x290(%%rdi)	\n\t	movaps	%%xmm11,0x190(%%rdi)	\n\t"\
		"movaps	%%xmm13,0x300(%%rdi)	\n\t	movaps	%%xmm15,0x280(%%rdi)	\n\t	movaps	%%xmm12,0x180(%%rdi)	\n\t"\
		"movaps	%%xmm5 ,0x090(%%rdi)	\n\t	movaps	%%xmm7 ,0x110(%%rdi)	\n\t	movaps	%%xmm4 ,0x210(%%rdi)	\n\t"\
	"/* SSE2_RADIX_07_DFT(14,10,06,02,26,22,18 -> 02,06,10,14,18,22,26): */\n\t"\
		"addq	$0x20,%%rdi	/* s1p02r */\n\t	/*** Imaginary Parts: ***/			\n\t"\
		"movaps	0x100(%%rdi),%%xmm6		\n\t	movaps	0x110(%%rdi),%%xmm14		\n\t"\
		"movaps	0x200(%%rdi),%%xmm1		\n\t	movaps	0x210(%%rdi),%%xmm9 		\n\t"\
		"movaps	0x080(%%rdi),%%xmm5		\n\t	movaps	0x090(%%rdi),%%xmm13		\n\t"\
		"movaps	0x280(%%rdi),%%xmm2		\n\t	movaps	0x290(%%rdi),%%xmm10		\n\t"\
		"movaps	     (%%rdi),%%xmm4		\n\t	movaps	0x010(%%rdi),%%xmm12		\n\t"\
		"movaps	0x300(%%rdi),%%xmm3		\n\t	movaps	0x310(%%rdi),%%xmm11		\n\t"\
		"subpd	%%xmm1,%%xmm6			\n\t	subpd	%%xmm9 ,%%xmm14				\n\t"\
		"addpd	%%xmm1,%%xmm1			\n\t	addpd	%%xmm9  ,%%xmm9 			\n\t"\
		"addpd	%%xmm6,%%xmm1			\n\t	addpd	%%xmm14,%%xmm9  			\n\t"\
		"subpd	%%xmm2,%%xmm5			\n\t	subpd	%%xmm10,%%xmm13				\n\t"\
		"addpd	%%xmm2,%%xmm2			\n\t	addpd	%%xmm10,%%xmm10				\n\t"\
		"addpd	%%xmm5,%%xmm2			\n\t	addpd	%%xmm13,%%xmm10				\n\t"\
		"movaps	0x180(%%rdi),%%xmm0		\n\t	movaps	0x190(%%rdi),%%xmm8 		\n\t"\
		"subpd	%%xmm3,%%xmm4			\n\t	subpd	%%xmm11,%%xmm12				\n\t"\
		"addpd	%%xmm3,%%xmm3			\n\t	addpd	%%xmm11,%%xmm11				\n\t"\
		"addpd	%%xmm4,%%xmm3			\n\t	addpd	%%xmm12,%%xmm11				\n\t"\
		"/* Only alloc 2 scratch storageslots in calling code, so 16-reg version use 2 of the eventual outputs for Im parts */\n\t"\
		"movaps	%%xmm0,     (%%rsi)		\n\t	movaps	%%xmm8 ,0x080(%%rdi)		\n\t"\
		"movaps	%%xmm6,0x010(%%rsi)		\n\t	movaps	%%xmm14,0x100(%%rdi)		\n\t"\
		"addpd	      %%xmm1,%%xmm0		\n\t	addpd	      %%xmm9 ,%%xmm8  		\n\t"\
		"movaps	%%xmm5,%%xmm7			\n\t	movaps	%%xmm13,%%xmm15				\n\t"\
		"addpd	      %%xmm2,%%xmm3		\n\t	addpd	      %%xmm10,%%xmm11		\n\t"\
		"subpd	      %%xmm4,%%xmm5		\n\t	subpd	      %%xmm12,%%xmm13		\n\t"\
		"subpd	      %%xmm2,%%xmm1		\n\t	subpd	      %%xmm10,%%xmm9  		\n\t"\
		"subpd	      %%xmm7,%%xmm6		\n\t	subpd	      %%xmm15,%%xmm14		\n\t"\
		"addpd	      %%xmm2,%%xmm2		\n\t	addpd	      %%xmm10,%%xmm10		\n\t"\
		"addpd	      %%xmm7,%%xmm4		\n\t	addpd	      %%xmm15,%%xmm12		\n\t"\
		"addpd	      %%xmm3,%%xmm0		\n\t	addpd	      %%xmm11,%%xmm8  		\n\t"\
		"addpd	0x010(%%rsi),%%xmm5		\n\t	addpd	0x100(%%rdi),%%xmm13		\n\t"\
		"subpd	      %%xmm2,%%xmm3		\n\t	subpd	      %%xmm10,%%xmm11		\n\t"\
		"movaps	%%xmm4,%%xmm7			\n\t	movaps	%%xmm12,%%xmm15				\n\t"\
		"movaps	%%xmm0,      (%%rdi)	\n\t	movaps	%%xmm8 ,0x010(%%rdi)		\n\t"\
		"subpd	      %%xmm6,%%xmm4		\n\t	subpd	      %%xmm14,%%xmm12		\n\t"\
		"movaps	%%xmm1,%%xmm2			\n\t	movaps	%%xmm9 ,%%xmm10				\n\t"\
		"subpd	     (%%rsi),%%xmm0		\n\t	subpd	0x080(%%rdi),%%xmm8  		\n\t"\
		"mulpd	0x010(%%rax),%%xmm5		\n\t	mulpd	0x010(%%rax),%%xmm13		\n\t"\
		"addpd	      %%xmm3,%%xmm2		\n\t	addpd	      %%xmm11,%%xmm10		\n\t"\
		"mulpd	     (%%rcx),%%xmm3		\n\t	mulpd	     (%%rcx),%%xmm11		\n\t"\
		"mulpd	0x010(%%rdx),%%xmm4		\n\t	mulpd	0x010(%%rdx),%%xmm12		\n\t"\
		"mulpd	     (%%rbx),%%xmm1		\n\t	mulpd	     (%%rbx),%%xmm9  		\n\t"\
		"mulpd	0x010(%%rbx),%%xmm6		\n\t	mulpd	0x010(%%rbx),%%xmm14		\n\t"\
		"mulpd	     (%%rax),%%xmm0		\n\t	mulpd	     (%%rax),%%xmm8  		\n\t"\
		"mulpd	0x010(%%rcx),%%xmm7		\n\t	mulpd	0x010(%%rcx),%%xmm15		\n\t"\
		"mulpd	     (%%rdx),%%xmm2		\n\t	mulpd	     (%%rdx),%%xmm10		\n\t"\
		"addpd	     (%%rdi),%%xmm0		\n\t	addpd	0x010(%%rdi),%%xmm8  		\n\t"\
		"addpd	      %%xmm4,%%xmm6		\n\t	addpd	      %%xmm12,%%xmm14		\n\t"\
		"subpd	      %%xmm2,%%xmm1		\n\t	subpd	      %%xmm10,%%xmm9  		\n\t"\
		"subpd	      %%xmm7,%%xmm4		\n\t	subpd	      %%xmm15,%%xmm12		\n\t"\
		"subpd	      %%xmm2,%%xmm3		\n\t	subpd	      %%xmm10,%%xmm11		\n\t"\
		"movaps	%%xmm0,%%xmm2			\n\t	movaps	%%xmm8 ,%%xmm10				\n\t"\
		"movaps	%%xmm5,%%xmm7			\n\t	movaps	%%xmm13,%%xmm15				\n\t"\
		"addpd	      %%xmm1,%%xmm0		\n\t	addpd	      %%xmm9 ,%%xmm8  		\n\t"\
		"addpd	      %%xmm6,%%xmm5		\n\t	addpd	      %%xmm14,%%xmm13		\n\t"\
		"addpd	      %%xmm3,%%xmm1		\n\t	addpd	      %%xmm11,%%xmm9  		\n\t"\
		"addpd	      %%xmm4,%%xmm6		\n\t	addpd	      %%xmm12,%%xmm14		\n\t"\
		"addpd	      %%xmm2,%%xmm3		\n\t	addpd	      %%xmm10,%%xmm11		\n\t"\
		"addpd	      %%xmm7,%%xmm4		\n\t	addpd	      %%xmm15,%%xmm12		\n\t"\
		"subpd	      %%xmm1,%%xmm2		\n\t	subpd	      %%xmm9 ,%%xmm10		\n\t"\
		"subpd	      %%xmm6,%%xmm7		\n\t	subpd	      %%xmm14,%%xmm15		\n\t"\
		"/* xmm1,6,9,14 free ... Note the order reversal on the 3rd pair of outputs: */\n\t"\
		"subpd	%%xmm13,%%xmm0			\n\t	subpd	%%xmm15,%%xmm2			\n\t	subpd	%%xmm12,%%xmm3 			\n\t"\
		"subpd	%%xmm5 ,%%xmm8  		\n\t	subpd	%%xmm7 ,%%xmm10			\n\t	subpd	%%xmm4 ,%%xmm11			\n\t"\
		"addpd	%%xmm13,%%xmm13			\n\t	addpd	%%xmm15,%%xmm15			\n\t	addpd	%%xmm12,%%xmm12			\n\t"\
		"addpd	%%xmm5 ,%%xmm5			\n\t	addpd	%%xmm7 ,%%xmm7			\n\t	addpd	%%xmm4 ,%%xmm4 			\n\t"\
		"addpd	%%xmm0 ,%%xmm13			\n\t	addpd	%%xmm2 ,%%xmm15			\n\t	addpd	%%xmm3 ,%%xmm12			\n\t"\
		"addpd	%%xmm8 ,%%xmm5			\n\t	addpd	%%xmm10,%%xmm7			\n\t	addpd	%%xmm11,%%xmm4 			\n\t"\
		"movaps	%%xmm0 ,0x080(%%rdi)	\n\t	movaps	%%xmm2 ,0x100(%%rdi)	\n\t	movaps	%%xmm3 ,0x200(%%rdi)	\n\t"\
		"movaps	%%xmm8 ,0x310(%%rdi)	\n\t	movaps	%%xmm10,0x290(%%rdi)	\n\t	movaps	%%xmm11,0x190(%%rdi)	\n\t"\
		"movaps	%%xmm13,0x300(%%rdi)	\n\t	movaps	%%xmm15,0x280(%%rdi)	\n\t	movaps	%%xmm12,0x180(%%rdi)	\n\t"\
		"movaps	%%xmm5 ,0x090(%%rdi)	\n\t	movaps	%%xmm7 ,0x110(%%rdi)	\n\t	movaps	%%xmm4 ,0x210(%%rdi)	\n\t"\
	"/* SSE2_RADIX_07_DFT(07,03,27,23,19,15,11 -> 03,07,11,15,19,23,27): */\n\t"\
		"addq	$0x20,%%rdi	/* s1p03r */\n\t	/*** Imaginary Parts: ***/			\n\t"\
		"movaps	     (%%rdi),%%xmm6		\n\t	movaps	0x010(%%rdi),%%xmm14		\n\t"\
		"movaps	0x100(%%rdi),%%xmm1		\n\t	movaps	0x110(%%rdi),%%xmm9 		\n\t"\
		"movaps	0x300(%%rdi),%%xmm5		\n\t	movaps	0x310(%%rdi),%%xmm13		\n\t"\
		"movaps	0x180(%%rdi),%%xmm2		\n\t	movaps	0x190(%%rdi),%%xmm10		\n\t"\
		"movaps	0x280(%%rdi),%%xmm4		\n\t	movaps	0x290(%%rdi),%%xmm12		\n\t"\
		"movaps	0x200(%%rdi),%%xmm3		\n\t	movaps	0x210(%%rdi),%%xmm11		\n\t"\
		"subpd	%%xmm1,%%xmm6			\n\t	subpd	%%xmm9 ,%%xmm14				\n\t"\
		"addpd	%%xmm1,%%xmm1			\n\t	addpd	%%xmm9  ,%%xmm9 			\n\t"\
		"addpd	%%xmm6,%%xmm1			\n\t	addpd	%%xmm14,%%xmm9  			\n\t"\
		"subpd	%%xmm2,%%xmm5			\n\t	subpd	%%xmm10,%%xmm13				\n\t"\
		"addpd	%%xmm2,%%xmm2			\n\t	addpd	%%xmm10,%%xmm10				\n\t"\
		"addpd	%%xmm5,%%xmm2			\n\t	addpd	%%xmm13,%%xmm10				\n\t"\
		"movaps	0x080(%%rdi),%%xmm0		\n\t	movaps	0x090(%%rdi),%%xmm8 		\n\t"\
		"subpd	%%xmm3,%%xmm4			\n\t	subpd	%%xmm11,%%xmm12				\n\t"\
		"addpd	%%xmm3,%%xmm3			\n\t	addpd	%%xmm11,%%xmm11				\n\t"\
		"addpd	%%xmm4,%%xmm3			\n\t	addpd	%%xmm12,%%xmm11				\n\t"\
		"/* Only alloc 2 scratch storageslots in calling code, so 16-reg version use 2 of the eventual outputs for Im parts */\n\t"\
		"movaps	%%xmm0,     (%%rsi)		\n\t	movaps	%%xmm8 ,0x080(%%rdi)		\n\t"\
		"movaps	%%xmm6,0x010(%%rsi)		\n\t	movaps	%%xmm14,0x100(%%rdi)		\n\t"\
		"addpd	      %%xmm1,%%xmm0		\n\t	addpd	      %%xmm9 ,%%xmm8  		\n\t"\
		"movaps	%%xmm5,%%xmm7			\n\t	movaps	%%xmm13,%%xmm15				\n\t"\
		"addpd	      %%xmm2,%%xmm3		\n\t	addpd	      %%xmm10,%%xmm11		\n\t"\
		"subpd	      %%xmm4,%%xmm5		\n\t	subpd	      %%xmm12,%%xmm13		\n\t"\
		"subpd	      %%xmm2,%%xmm1		\n\t	subpd	      %%xmm10,%%xmm9  		\n\t"\
		"subpd	      %%xmm7,%%xmm6		\n\t	subpd	      %%xmm15,%%xmm14		\n\t"\
		"addpd	      %%xmm2,%%xmm2		\n\t	addpd	      %%xmm10,%%xmm10		\n\t"\
		"addpd	      %%xmm7,%%xmm4		\n\t	addpd	      %%xmm15,%%xmm12		\n\t"\
		"addpd	      %%xmm3,%%xmm0		\n\t	addpd	      %%xmm11,%%xmm8  		\n\t"\
		"addpd	0x010(%%rsi),%%xmm5		\n\t	addpd	0x100(%%rdi),%%xmm13		\n\t"\
		"subpd	      %%xmm2,%%xmm3		\n\t	subpd	      %%xmm10,%%xmm11		\n\t"\
		"movaps	%%xmm4,%%xmm7			\n\t	movaps	%%xmm12,%%xmm15				\n\t"\
		"movaps	%%xmm0,      (%%rdi)	\n\t	movaps	%%xmm8 ,0x010(%%rdi)		\n\t"\
		"subpd	      %%xmm6,%%xmm4		\n\t	subpd	      %%xmm14,%%xmm12		\n\t"\
		"movaps	%%xmm1,%%xmm2			\n\t	movaps	%%xmm9 ,%%xmm10				\n\t"\
		"subpd	     (%%rsi),%%xmm0		\n\t	subpd	0x080(%%rdi),%%xmm8  		\n\t"\
		"mulpd	0x010(%%rax),%%xmm5		\n\t	mulpd	0x010(%%rax),%%xmm13		\n\t"\
		"addpd	      %%xmm3,%%xmm2		\n\t	addpd	      %%xmm11,%%xmm10		\n\t"\
		"mulpd	     (%%rcx),%%xmm3		\n\t	mulpd	     (%%rcx),%%xmm11		\n\t"\
		"mulpd	0x010(%%rdx),%%xmm4		\n\t	mulpd	0x010(%%rdx),%%xmm12		\n\t"\
		"mulpd	     (%%rbx),%%xmm1		\n\t	mulpd	     (%%rbx),%%xmm9  		\n\t"\
		"mulpd	0x010(%%rbx),%%xmm6		\n\t	mulpd	0x010(%%rbx),%%xmm14		\n\t"\
		"mulpd	     (%%rax),%%xmm0		\n\t	mulpd	     (%%rax),%%xmm8  		\n\t"\
		"mulpd	0x010(%%rcx),%%xmm7		\n\t	mulpd	0x010(%%rcx),%%xmm15		\n\t"\
		"mulpd	     (%%rdx),%%xmm2		\n\t	mulpd	     (%%rdx),%%xmm10		\n\t"\
		"addpd	     (%%rdi),%%xmm0		\n\t	addpd	0x010(%%rdi),%%xmm8  		\n\t"\
		"addpd	      %%xmm4,%%xmm6		\n\t	addpd	      %%xmm12,%%xmm14		\n\t"\
		"subpd	      %%xmm2,%%xmm1		\n\t	subpd	      %%xmm10,%%xmm9  		\n\t"\
		"subpd	      %%xmm7,%%xmm4		\n\t	subpd	      %%xmm15,%%xmm12		\n\t"\
		"subpd	      %%xmm2,%%xmm3		\n\t	subpd	      %%xmm10,%%xmm11		\n\t"\
		"movaps	%%xmm0,%%xmm2			\n\t	movaps	%%xmm8 ,%%xmm10				\n\t"\
		"movaps	%%xmm5,%%xmm7			\n\t	movaps	%%xmm13,%%xmm15				\n\t"\
		"addpd	      %%xmm1,%%xmm0		\n\t	addpd	      %%xmm9 ,%%xmm8  		\n\t"\
		"addpd	      %%xmm6,%%xmm5		\n\t	addpd	      %%xmm14,%%xmm13		\n\t"\
		"addpd	      %%xmm3,%%xmm1		\n\t	addpd	      %%xmm11,%%xmm9  		\n\t"\
		"addpd	      %%xmm4,%%xmm6		\n\t	addpd	      %%xmm12,%%xmm14		\n\t"\
		"addpd	      %%xmm2,%%xmm3		\n\t	addpd	      %%xmm10,%%xmm11		\n\t"\
		"addpd	      %%xmm7,%%xmm4		\n\t	addpd	      %%xmm15,%%xmm12		\n\t"\
		"subpd	      %%xmm1,%%xmm2		\n\t	subpd	      %%xmm9 ,%%xmm10		\n\t"\
		"subpd	      %%xmm6,%%xmm7		\n\t	subpd	      %%xmm14,%%xmm15		\n\t"\
		"/* xmm1,6,9,14 free ... Note the order reversal on the 3rd pair of outputs: */\n\t"\
		"subpd	%%xmm13,%%xmm0			\n\t	subpd	%%xmm15,%%xmm2			\n\t	subpd	%%xmm12,%%xmm3 			\n\t"\
		"subpd	%%xmm5 ,%%xmm8  		\n\t	subpd	%%xmm7 ,%%xmm10			\n\t	subpd	%%xmm4 ,%%xmm11			\n\t"\
		"addpd	%%xmm13,%%xmm13			\n\t	addpd	%%xmm15,%%xmm15			\n\t	addpd	%%xmm12,%%xmm12			\n\t"\
		"addpd	%%xmm5 ,%%xmm5			\n\t	addpd	%%xmm7 ,%%xmm7			\n\t	addpd	%%xmm4 ,%%xmm4 			\n\t"\
		"addpd	%%xmm0 ,%%xmm13			\n\t	addpd	%%xmm2 ,%%xmm15			\n\t	addpd	%%xmm3 ,%%xmm12			\n\t"\
		"addpd	%%xmm8 ,%%xmm5			\n\t	addpd	%%xmm10,%%xmm7			\n\t	addpd	%%xmm11,%%xmm4 			\n\t"\
		"movaps	%%xmm0 ,0x080(%%rdi)	\n\t	movaps	%%xmm2 ,0x100(%%rdi)	\n\t	movaps	%%xmm3 ,0x200(%%rdi)	\n\t"\
		"movaps	%%xmm8 ,0x310(%%rdi)	\n\t	movaps	%%xmm10,0x290(%%rdi)	\n\t	movaps	%%xmm11,0x190(%%rdi)	\n\t"\
		"movaps	%%xmm13,0x300(%%rdi)	\n\t	movaps	%%xmm15,0x280(%%rdi)	\n\t	movaps	%%xmm12,0x180(%%rdi)	\n\t"\
		"movaps	%%xmm5 ,0x090(%%rdi)	\n\t	movaps	%%xmm7 ,0x110(%%rdi)	\n\t	movaps	%%xmm4 ,0x210(%%rdi)	\n\t"\
		"\n\t"\
	"/* Since doing radix-7 in-place here, inputs of radix-4 are in consecutive memory locs, i.e. 0x20 bytes apart, e.g. the distance between s1p00r and s1p01r: */\n\t"\
		"\n\t"\
	"/* 1:SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_E(s1p00,01,02,03 -> add00+p0123): */	\n\t"\
		"movq	%[__out],%%rsi	/* s1p00r */						\n\t"\
		"movq	%[__add],%%rax										\n\t"\
		"movslq	%[__p01],%%rbx										\n\t"\
		"movslq	%[__p02],%%rcx										\n\t"\
		"movslq	%[__p03],%%rdx										\n\t"\
		"shlq	$3,%%rbx	/* Ptr offset for floating doubles */	\n\t"\
		"shlq	$3,%%rcx											\n\t"\
		"shlq	$3,%%rdx											\n\t"\
		"addq	%%rax,%%rbx											\n\t"\
		"addq	%%rax,%%rcx											\n\t"\
		"addq	%%rax,%%rdx											\n\t"\
		"movaps	     (%%rsi),%%xmm0									\n\t"\
		"movaps	0x020(%%rsi),%%xmm4									\n\t"\
		"movaps	0x010(%%rsi),%%xmm1									\n\t"\
		"movaps	0x030(%%rsi),%%xmm5									\n\t"\
		"movaps	%%xmm0,%%xmm2										\n\t"\
		"movaps	%%xmm4,%%xmm6										\n\t"\
		"movaps	%%xmm1,%%xmm3										\n\t"\
		"movaps	%%xmm5,%%xmm7										\n\t"\
		"addpd	0x040(%%rsi),%%xmm0									\n\t"\
		"addpd	0x060(%%rsi),%%xmm4									\n\t"\
		"addpd	0x050(%%rsi),%%xmm1									\n\t"\
		"addpd	0x070(%%rsi),%%xmm5									\n\t"\
		"subpd	0x040(%%rsi),%%xmm2									\n\t	/* 2:SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_E(s1p04,05,06,07 -> add24+p1032): */\n\t"\
		"subpd	0x060(%%rsi),%%xmm6									\n\t		movslq	%[__p24],%%rdi										\n\t"\
		"subpd	0x050(%%rsi),%%xmm3									\n\t		shlq	$3,%%rdi											\n\t"\
		"subpd	0x070(%%rsi),%%xmm7									\n\t		addq	$0x080,%%rsi										\n\t"\
		"/* Finish radix-4 butterfly and store: */					\n\t		movaps	     (%%rsi),%%xmm8 								\n\t"\
		"subpd	%%xmm4,%%xmm0										\n\t		movaps	0x020(%%rsi),%%xmm12								\n\t"\
		"subpd	%%xmm7,%%xmm2										\n\t		movaps	0x010(%%rsi),%%xmm9 								\n\t"\
		"subpd	%%xmm5,%%xmm1										\n\t		movaps	0x030(%%rsi),%%xmm13								\n\t"\
		"subpd	%%xmm6,%%xmm3										\n\t		movaps	%%xmm8 ,%%xmm10										\n\t"\
		"movaps	%%xmm0,     (%%rbx)									\n\t		movaps	%%xmm12,%%xmm14										\n\t"\
		"movaps	%%xmm2,     (%%rcx)									\n\t		movaps	%%xmm9 ,%%xmm11										\n\t"\
		"movaps	%%xmm1,0x010(%%rbx)									\n\t		movaps	%%xmm13,%%xmm15										\n\t"\
		"movaps	%%xmm3,0x010(%%rdx)									\n\t		addpd	0x040(%%rsi),%%xmm8 								\n\t"\
		"addpd	%%xmm4,%%xmm4										\n\t		addpd	0x060(%%rsi),%%xmm12								\n\t"\
		"addpd	%%xmm7,%%xmm7										\n\t		addpd	0x050(%%rsi),%%xmm9 								\n\t"\
		"addpd	%%xmm5,%%xmm5										\n\t		addpd	0x070(%%rsi),%%xmm13								\n\t"\
		"addpd	%%xmm6,%%xmm6										\n\t		subpd	0x040(%%rsi),%%xmm10								\n\t"\
		"addpd	%%xmm0,%%xmm4										\n\t		subpd	0x060(%%rsi),%%xmm14								\n\t"\
		"addpd	%%xmm2,%%xmm7										\n\t		subpd	0x050(%%rsi),%%xmm11								\n\t"\
		"addpd	%%xmm1,%%xmm5										\n\t		subpd	0x070(%%rsi),%%xmm15								\n\t"\
		"addpd	%%xmm3,%%xmm6										\n\t		/* Finish radix-4 butterfly and store: */					\n\t"\
		"movaps	%%xmm4,     (%%rax)									\n\t		addq	%%rdi,%%rbx		/* reorder to not conflict with left-col store addresses */\n\t"\
		"movaps	%%xmm7,     (%%rdx)									\n\t		addq	%%rdi,%%rdx											\n\t"\
		"movaps	%%xmm5,0x010(%%rax)									\n\t		addq	%%rdi,%%rax											\n\t"\
		"movaps	%%xmm6,0x010(%%rcx)									\n\t		addq	%%rdi,%%rcx											\n\t"\
	"/* 3:SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_E(s1p08,09,10,11 -> add20+p2310): */\n\t	/* Swap _o1 <-> _o2, _o3 <-> _o4: eax <-> ebx, ecx <-> edx */\n\t"\
		"addq	$0x080,%%rsi										\n\t		subpd	%%xmm12,%%xmm8 										\n\t"\
		"movslq	%[__p04],%%rdi										\n\t		subpd	%%xmm15,%%xmm10										\n\t"\
		"shlq	$3,%%rdi											\n\t		subpd	%%xmm13,%%xmm9 										\n\t"\
		"movaps	     (%%rsi),%%xmm0									\n\t		subpd	%%xmm14,%%xmm11										\n\t"\
		"movaps	0x020(%%rsi),%%xmm4									\n\t		movaps	%%xmm8 ,     (%%rax)								\n\t"\
		"movaps	0x010(%%rsi),%%xmm1									\n\t		movaps	%%xmm10,     (%%rdx)								\n\t"\
		"movaps	0x030(%%rsi),%%xmm5									\n\t		movaps	%%xmm9 ,0x010(%%rax)								\n\t"\
		"movaps	%%xmm0,%%xmm2										\n\t		movaps	%%xmm11,0x010(%%rcx)								\n\t"\
		"movaps	%%xmm4,%%xmm6										\n\t		addpd	%%xmm12,%%xmm12										\n\t"\
		"movaps	%%xmm1,%%xmm3										\n\t		addpd	%%xmm15,%%xmm15										\n\t"\
		"movaps	%%xmm5,%%xmm7										\n\t		addpd	%%xmm13,%%xmm13										\n\t"\
		"addpd	0x040(%%rsi),%%xmm0									\n\t		addpd	%%xmm14,%%xmm14										\n\t"\
		"addpd	0x060(%%rsi),%%xmm4									\n\t		addpd	%%xmm8 ,%%xmm12										\n\t"\
		"addpd	0x050(%%rsi),%%xmm1									\n\t		addpd	%%xmm10,%%xmm15										\n\t"\
		"addpd	0x070(%%rsi),%%xmm5									\n\t		addpd	%%xmm9 ,%%xmm13										\n\t"\
		"subpd	0x040(%%rsi),%%xmm2									\n\t		addpd	%%xmm11,%%xmm14										\n\t"\
		"subpd	0x060(%%rsi),%%xmm6									\n\t		movaps	%%xmm12,     (%%rbx)								\n\t"\
		"subpd	0x050(%%rsi),%%xmm3									\n\t		movaps	%%xmm15,     (%%rcx)								\n\t"\
		"subpd	0x070(%%rsi),%%xmm7									\n\t		movaps	%%xmm13,0x010(%%rbx)								\n\t"\
		"/* Finish radix-4 butterfly and store: */					\n\t		movaps	%%xmm14,0x010(%%rdx)								\n\t"\
		"subq	%%rdi,%%rax											\n\t"\
		"subq	%%rdi,%%rbx											\n\t"\
		"subq	%%rdi,%%rcx											\n\t	/* 4:SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_E(s1p12,13,14,15 -> add16+p0123): */\n\t"\
		"subq	%%rdi,%%rdx											\n\t		addq	$0x080,%%rsi										\n\t"\
		"/* Swap _o0123 <-> _o2310: e[abcd]x <-> e[cdba]x */		\n\t		movaps	     (%%rsi),%%xmm8 								\n\t"\
		"subpd	%%xmm4,%%xmm0										\n\t		movaps	0x020(%%rsi),%%xmm12								\n\t"\
		"subpd	%%xmm7,%%xmm2										\n\t		movaps	0x010(%%rsi),%%xmm9 								\n\t"\
		"subpd	%%xmm5,%%xmm1										\n\t		movaps	0x030(%%rsi),%%xmm13								\n\t"\
		"subpd	%%xmm6,%%xmm3										\n\t		movaps	%%xmm8 ,%%xmm10										\n\t"\
		"movaps	%%xmm0,     (%%rdx)									\n\t		movaps	%%xmm12,%%xmm14										\n\t"\
		"movaps	%%xmm2,     (%%rbx)									\n\t		movaps	%%xmm9 ,%%xmm11										\n\t"\
		"movaps	%%xmm1,0x010(%%rdx)									\n\t		movaps	%%xmm13,%%xmm15										\n\t"\
		"movaps	%%xmm3,0x010(%%rax)									\n\t		addpd	0x040(%%rsi),%%xmm8 								\n\t"\
		"addpd	%%xmm4,%%xmm4										\n\t		addpd	0x060(%%rsi),%%xmm12								\n\t"\
		"addpd	%%xmm7,%%xmm7										\n\t		addpd	0x050(%%rsi),%%xmm9 								\n\t"\
		"addpd	%%xmm5,%%xmm5										\n\t		addpd	0x070(%%rsi),%%xmm13								\n\t"\
		"addpd	%%xmm6,%%xmm6										\n\t		subpd	0x040(%%rsi),%%xmm10								\n\t"\
		"addpd	%%xmm0,%%xmm4										\n\t		subpd	0x060(%%rsi),%%xmm14								\n\t"\
		"addpd	%%xmm2,%%xmm7										\n\t		subpd	0x050(%%rsi),%%xmm11								\n\t"\
		"addpd	%%xmm1,%%xmm5										\n\t		subpd	0x070(%%rsi),%%xmm15								\n\t"\
		"addpd	%%xmm3,%%xmm6										\n\t		/* Finish radix-4 butterfly and store: */					\n\t"\
		"movaps	%%xmm4,     (%%rcx)									\n\t		subq	%%rdi,%%rdx											\n\t"\
		"movaps	%%xmm7,     (%%rax)									\n\t		subq	%%rdi,%%rax											\n\t"\
		"movaps	%%xmm5,0x010(%%rcx)									\n\t		subq	%%rdi,%%rcx											\n\t"\
		"movaps	%%xmm6,0x010(%%rbx)									\n\t		subq	%%rdi,%%rbx											\n\t"\
	"/* 5:SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_E(s1p16,17,18,19 -> add12+p3201): */\n\t	subpd	%%xmm12,%%xmm8 									\n\t"\
		"addq	$0x080,%%rsi										\n\t		subpd	%%xmm15,%%xmm10										\n\t"\
		"movaps	     (%%rsi),%%xmm0									\n\t		subpd	%%xmm13,%%xmm9 										\n\t"\
		"movaps	0x020(%%rsi),%%xmm4									\n\t		subpd	%%xmm14,%%xmm11										\n\t"\
		"movaps	0x010(%%rsi),%%xmm1									\n\t		movaps	%%xmm8 ,     (%%rbx)								\n\t"\
		"movaps	0x030(%%rsi),%%xmm5									\n\t		movaps	%%xmm10,     (%%rcx)								\n\t"\
		"movaps	%%xmm0,%%xmm2										\n\t		movaps	%%xmm9 ,0x010(%%rbx)								\n\t"\
		"movaps	%%xmm4,%%xmm6										\n\t		movaps	%%xmm11,0x010(%%rdx)								\n\t"\
		"movaps	%%xmm1,%%xmm3										\n\t		addpd	%%xmm12,%%xmm12										\n\t"\
		"movaps	%%xmm5,%%xmm7										\n\t		addpd	%%xmm15,%%xmm15										\n\t"\
		"addpd	0x040(%%rsi),%%xmm0									\n\t		addpd	%%xmm13,%%xmm13										\n\t"\
		"addpd	0x060(%%rsi),%%xmm4									\n\t		addpd	%%xmm14,%%xmm14										\n\t"\
		"addpd	0x050(%%rsi),%%xmm1									\n\t		addpd	%%xmm8 ,%%xmm12										\n\t"\
		"addpd	0x070(%%rsi),%%xmm5									\n\t		addpd	%%xmm10,%%xmm15										\n\t"\
		"subpd	0x040(%%rsi),%%xmm2									\n\t		addpd	%%xmm9 ,%%xmm13										\n\t"\
		"subpd	0x060(%%rsi),%%xmm6									\n\t		addpd	%%xmm11,%%xmm14										\n\t"\
		"subpd	0x050(%%rsi),%%xmm3									\n\t		movaps	%%xmm12,     (%%rax)								\n\t"\
		"subpd	0x070(%%rsi),%%xmm7									\n\t		movaps	%%xmm15,     (%%rdx)								\n\t"\
		"/* Finish radix-4 butterfly and store: */					\n\t		movaps	%%xmm13,0x010(%%rax)								\n\t"\
		"subq	%%rdi,%%rax											\n\t		movaps	%%xmm14,0x010(%%rcx)								\n\t"\
		"subq	%%rdi,%%rbx											\n\t"\
		"subq	%%rdi,%%rcx											\n\t	/* 6:SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_E(s1p20,21,22,23 -> add08+p1032): */\n\t"\
		"subq	%%rdi,%%rdx											\n\t		addq	$0x080,%%rsi										\n\t"\
		"/* Swap _o0123 <-> _o3201: e[abcd]x <-> e[dcab]x */		\n\t		movaps	     (%%rsi),%%xmm8 								\n\t"\
		"subpd	%%xmm4,%%xmm0										\n\t		movaps	0x020(%%rsi),%%xmm12								\n\t"\
		"subpd	%%xmm7,%%xmm2										\n\t		movaps	0x010(%%rsi),%%xmm9 								\n\t"\
		"subpd	%%xmm5,%%xmm1										\n\t		movaps	0x030(%%rsi),%%xmm13								\n\t"\
		"subpd	%%xmm6,%%xmm3										\n\t		movaps	%%xmm8 ,%%xmm10										\n\t"\
		"movaps	%%xmm0,     (%%rcx)									\n\t		movaps	%%xmm12,%%xmm14										\n\t"\
		"movaps	%%xmm2,     (%%rax)									\n\t		movaps	%%xmm9 ,%%xmm11										\n\t"\
		"movaps	%%xmm1,0x010(%%rcx)									\n\t		movaps	%%xmm13,%%xmm15										\n\t"\
		"movaps	%%xmm3,0x010(%%rbx)									\n\t		addpd	0x040(%%rsi),%%xmm8 								\n\t"\
		"addpd	%%xmm4,%%xmm4										\n\t		addpd	0x060(%%rsi),%%xmm12								\n\t"\
		"addpd	%%xmm7,%%xmm7										\n\t		addpd	0x050(%%rsi),%%xmm9 								\n\t"\
		"addpd	%%xmm5,%%xmm5										\n\t		addpd	0x070(%%rsi),%%xmm13								\n\t"\
		"addpd	%%xmm6,%%xmm6										\n\t		subpd	0x040(%%rsi),%%xmm10								\n\t"\
		"addpd	%%xmm0,%%xmm4										\n\t		subpd	0x060(%%rsi),%%xmm14								\n\t"\
		"addpd	%%xmm2,%%xmm7										\n\t		subpd	0x050(%%rsi),%%xmm11								\n\t"\
		"addpd	%%xmm1,%%xmm5										\n\t		subpd	0x070(%%rsi),%%xmm15								\n\t"\
		"addpd	%%xmm3,%%xmm6										\n\t		/* Finish radix-4 butterfly and store: */					\n\t"\
		"movaps	%%xmm4,     (%%rdx)									\n\t		subq	%%rdi,%%rcx											\n\t"\
		"movaps	%%xmm7,     (%%rbx)									\n\t		subq	%%rdi,%%rbx											\n\t"\
		"movaps	%%xmm5,0x010(%%rdx)									\n\t		subq	%%rdi,%%rdx											\n\t"\
		"movaps	%%xmm6,0x010(%%rax)									\n\t		subq	%%rdi,%%rax											\n\t"\
	"/* 7:SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_E(s1p24,25,26,27 -> add04+p2310): */\n\t	/* Swap _o0 <-> _o1, _o2 <-> _o3: eax <-> ebx, ecx <-> edx */\n\t"\
		"addq	$0x080,%%rsi										\n\t		subpd	%%xmm12,%%xmm8 										\n\t"\
		"movaps	     (%%rsi),%%xmm0									\n\t		subpd	%%xmm15,%%xmm10										\n\t"\
		"movaps	0x020(%%rsi),%%xmm4									\n\t		subpd	%%xmm13,%%xmm9 										\n\t"\
		"movaps	0x010(%%rsi),%%xmm1									\n\t		subpd	%%xmm14,%%xmm11										\n\t"\
		"movaps	0x030(%%rsi),%%xmm5									\n\t		movaps	%%xmm8 ,     (%%rax)								\n\t"\
		"movaps	%%xmm0,%%xmm2										\n\t		movaps	%%xmm10,     (%%rdx)								\n\t"\
		"movaps	%%xmm4,%%xmm6										\n\t		movaps	%%xmm9 ,0x010(%%rax)								\n\t"\
		"movaps	%%xmm1,%%xmm3										\n\t		movaps	%%xmm11,0x010(%%rcx)								\n\t"\
		"movaps	%%xmm5,%%xmm7										\n\t		addpd	%%xmm12,%%xmm12										\n\t"\
		"addpd	0x040(%%rsi),%%xmm0									\n\t		addpd	%%xmm15,%%xmm15										\n\t"\
		"addpd	0x060(%%rsi),%%xmm4									\n\t		addpd	%%xmm13,%%xmm13										\n\t"\
		"addpd	0x050(%%rsi),%%xmm1									\n\t		addpd	%%xmm14,%%xmm14										\n\t"\
		"addpd	0x070(%%rsi),%%xmm5									\n\t		addpd	%%xmm8 ,%%xmm12										\n\t"\
		"subpd	0x040(%%rsi),%%xmm2									\n\t		addpd	%%xmm10,%%xmm15										\n\t"\
		"subpd	0x060(%%rsi),%%xmm6									\n\t		addpd	%%xmm9 ,%%xmm13										\n\t"\
		"subpd	0x050(%%rsi),%%xmm3									\n\t		addpd	%%xmm11,%%xmm14										\n\t"\
		"subpd	0x070(%%rsi),%%xmm7									\n\t		movaps	%%xmm12,     (%%rbx)								\n\t"\
		"/* Finish radix-4 butterfly and store: */					\n\t		movaps	%%xmm15,     (%%rcx)								\n\t"\
		"subq	%%rdi,%%rax											\n\t		movaps	%%xmm13,0x010(%%rbx)								\n\t"\
		"subq	%%rdi,%%rcx											\n\t		movaps	%%xmm14,0x010(%%rdx)								\n\t"\
		"subq	%%rdi,%%rbx											\n\t"\
		"subq	%%rdi,%%rdx											\n\t"\
		"/* Swap _o0123 <-> _o2310: e[abcd]x <-> e[cdba]x */		\n\t"\
		"subpd	%%xmm4,%%xmm0										\n\t"\
		"subpd	%%xmm7,%%xmm2										\n\t"\
		"subpd	%%xmm5,%%xmm1										\n\t"\
		"subpd	%%xmm6,%%xmm3										\n\t"\
		"movaps	%%xmm0,     (%%rdx)									\n\t"\
		"movaps	%%xmm2,     (%%rbx)									\n\t"\
		"movaps	%%xmm1,0x010(%%rdx)									\n\t"\
		"movaps	%%xmm3,0x010(%%rax)									\n\t"\
		"addpd	%%xmm4,%%xmm4										\n\t"\
		"addpd	%%xmm7,%%xmm7										\n\t"\
		"addpd	%%xmm5,%%xmm5										\n\t"\
		"addpd	%%xmm6,%%xmm6										\n\t"\
		"addpd	%%xmm0,%%xmm4										\n\t"\
		"addpd	%%xmm2,%%xmm7										\n\t"\
		"addpd	%%xmm1,%%xmm5										\n\t"\
		"addpd	%%xmm3,%%xmm6										\n\t"\
		"movaps	%%xmm4,     (%%rcx)									\n\t"\
		"movaps	%%xmm7,     (%%rax)									\n\t"\
		"movaps	%%xmm5,0x010(%%rcx)									\n\t"\
		"movaps	%%xmm6,0x010(%%rbx)									\n\t"\
		:					/* outputs: none */\
		: [__add] "m" (Xadd)	/* All inputs from memory addresses here */\
		 ,[__p01] "m" (Xp01)\
		 ,[__p02] "m" (Xp02)\
		 ,[__p03] "m" (Xp03)\
		 ,[__p04] "m" (Xp04)\
		 ,[__p08] "m" (Xp08)\
		 ,[__p12] "m" (Xp12)\
		 ,[__p16] "m" (Xp16)\
		 ,[__p20] "m" (Xp20)\
		 ,[__p24] "m" (Xp24)\
		 ,[__out] "m" (Xout)\
		 ,[__cc0] "m" (Xcc0)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"		/* Clobbered registers */\
	);\
	}

#endif

#endif	/* radix28_ditN_cy_dif1_gcc_h_included */
