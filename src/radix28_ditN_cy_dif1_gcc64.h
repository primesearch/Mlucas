/*******************************************************************************
*                                                                              *
*   (C) 1997-2019 by Ernst W. Mayer.                                           *
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

  #define USE_64BIT_ASM_STYLE	1	// My x86 timings indicate fancier versions below using ymm0-15 for the radix-16 DFT is faster.

  #if USE_64BIT_ASM_STYLE	// True: Deeper 64-bit-ified version of the 32-bit version of the ASM macros, using all of ymm0-15

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

  #else // USE_64BIT_ASM_STYLE = False: Use simple 64-bit-ified version of the analogous 32-bit ASM macros, i.e. using just ymm0-7.

	#define	SSE2_RADIX28_DIT_NOTWIDDLE(Xadd,Xp01,Xp02,Xp03,Xp04,Xp08,Xp12,Xp16,Xp20,Xp24,Xout,Xcc0)\
	{\
	__asm__ volatile (\
		"/* SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_E(add00+p[0,1,3,2], s1p00r,s1p03r,s1p02r,s1p01r): */\n\t"\
			"movq	%[__out],%%rsi	/* s1p00r */\n\t"\
			"movq	%[__add],%%rax\n\t"\
			"movslq	%[__p01],%%rbx\n\t"\
			"movslq	%[__p02],%%rcx\n\t"\
			"movslq	%[__p03],%%rdx\n\t"\
			"shlq	$3,%%rbx		/* Pointer offset for floating doubles */\n\t"\
			"shlq	$3,%%rcx\n\t"\
			"shlq	$3,%%rdx\n\t"\
			"addq	%%rax,%%rbx\n\t"\
			"addq	%%rax,%%rcx\n\t"\
			"addq	%%rax,%%rdx\n\t"\
			"/* ecx <-> edx */\n\t"\
			"vmovaps	     (%%rax),%%ymm0		/* a(jt   ) */\n\t"\
			"vmovaps	     (%%rdx),%%ymm4		/* a(jt+p2) */\n\t"\
			"vmovaps	0x020(%%rax),%%ymm1		/* a(jp   ) */\n\t"\
			"vmovaps	0x020(%%rdx),%%ymm5		/* a(jp+p2) */\n\t"\
			"vmovaps	%%ymm0,%%ymm2			/* ymm2 <- cpy a(jt   ) */\n\t"\
			"vmovaps	%%ymm4,%%ymm6			/* ymm4 <- cpy a(jt+p2) */\n\t"\
			"vmovaps	%%ymm1,%%ymm3			/* ymm3 <- cpy a(jp   ) */\n\t"\
			"vmovaps	%%ymm5,%%ymm7			/* ymm5 <- cpy a(jp+p2) */\n\t"\
			"vaddpd	     (%%rbx),%%ymm0,%%ymm0		/* t1 */\n\t"\
			"vaddpd	     (%%rcx),%%ymm4,%%ymm4		/* t5 */\n\t"\
			"vaddpd	0x020(%%rbx),%%ymm1,%%ymm1		/* t2 */\n\t"\
			"vaddpd	0x020(%%rcx),%%ymm5,%%ymm5		/* t6 */\n\t"\
			"vsubpd	     (%%rbx),%%ymm2,%%ymm2		/* t3 */\n\t"\
			"vsubpd	     (%%rcx),%%ymm6,%%ymm6		/* t7 */\n\t"\
			"vsubpd	0x020(%%rbx),%%ymm3,%%ymm3		/* t4 */\n\t"\
			"vsubpd	0x020(%%rcx),%%ymm7,%%ymm7		/* t8 */\n\t"\
			"/* Finish radix-4 butterfly and store results into output-array slots: */\n\t"\
			"/* Swap _o1 <-> _o3: 0x040,0x060 <-> 0x0c0,0x0e0: */\n\t"\
			"vsubpd	%%ymm4,%%ymm0,%%ymm0			/* ~t5 <- t1 -t5 */\n\t"\
			"vsubpd	%%ymm7,%%ymm2,%%ymm2			/* ~t7 <- t3 -t8 */\n\t"\
			"vsubpd	%%ymm5,%%ymm1,%%ymm1			/* ~t6 <- t2 -t6 */\n\t"\
			"vsubpd	%%ymm6,%%ymm3,%%ymm3			/* ~t4 <- t4 -t7 */\n\t"\
			"vmovaps	%%ymm0,0x080(%%rsi)		/* <- ~t5 */\n\t"\
			"vmovaps	%%ymm2,0x040(%%rsi)		/* <- ~t7 */\n\t"\
			"vmovaps	%%ymm1,0x0a0(%%rsi)		/* <- ~t6 */\n\t"\
			"vmovaps	%%ymm3,0x0e0(%%rsi)		/* <- ~t4 */\n\t"\
			"vaddpd	%%ymm4,%%ymm4,%%ymm4			/*          2*t5 */\n\t"\
			"vaddpd	%%ymm7,%%ymm7,%%ymm7			/*          2*t8 */\n\t"\
			"vaddpd	%%ymm5,%%ymm5,%%ymm5			/*          2*t6 */\n\t"\
			"vaddpd	%%ymm6,%%ymm6,%%ymm6			/*          2*t7 */\n\t"\
			"vaddpd	%%ymm0,%%ymm4,%%ymm4			/* ~t1 <- t1 +t5 */\n\t"\
			"vaddpd	%%ymm2,%%ymm7,%%ymm7			/* ~t3 <- t3 +t8 */\n\t"\
			"vaddpd	%%ymm1,%%ymm5,%%ymm5			/* ~t2 <- t2 +t6 */\n\t"\
			"vaddpd	%%ymm3,%%ymm6,%%ymm6			/* ~t8 <- t4 +t7 */\n\t"\
			"vmovaps	%%ymm4,      (%%rsi)		/* <- ~t1 */\n\t"\
			"vmovaps	%%ymm7,0x0c0(%%rsi)		/* <- ~t3 */\n\t"\
			"vmovaps	%%ymm5,0x020(%%rsi)		/* <- ~t2 */\n\t"\
			"vmovaps	%%ymm6,0x060(%%rsi)		/* <- ~t8 */\n\t"\
		"/* SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_E(add12+p[3,2,1,0], s1p04r,s1p07r,s1p06r,s1p05r): */\n\t"\
			"addq	$0x100,%%rsi	/* s1p04r */\n\t"\
			"subq	%%rax,%%rbx		/* p01 */\n\t"\
			"subq	%%rax,%%rcx		/* p02 */\n\t"\
			"subq	%%rax,%%rdx		/* p03 */\n\t"\
			"movslq	%[__p12],%%rdi\n\t"\
			"shlq	$3,%%rdi		/* Pointer offset for floating doubles */\n\t"\
			"addq	%%rdi,%%rax		/* add12 */\n\t"\
			"addq	%%rax,%%rbx\n\t"\
			"addq	%%rax,%%rcx\n\t"\
			"addq	%%rax,%%rdx\n\t"\
			"/* eax <-> edx, ebx <-> ecx */\n\t"\
			"vmovaps	     (%%rdx),%%ymm0		/* a(jt   ) */\n\t"\
			"vmovaps	     (%%rbx),%%ymm4		/* a(jt+p2) */\n\t"\
			"vmovaps	0x020(%%rdx),%%ymm1		/* a(jp   ) */\n\t"\
			"vmovaps	0x020(%%rbx),%%ymm5		/* a(jp+p2) */\n\t"\
			"vmovaps	%%ymm0,%%ymm2			/* ymm2 <- cpy a(jt   ) */\n\t"\
			"vmovaps	%%ymm4,%%ymm6			/* ymm4 <- cpy a(jt+p2) */\n\t"\
			"vmovaps	%%ymm1,%%ymm3			/* ymm3 <- cpy a(jp   ) */\n\t"\
			"vmovaps	%%ymm5,%%ymm7			/* ymm5 <- cpy a(jp+p2) */\n\t"\
			"vaddpd	     (%%rcx),%%ymm0,%%ymm0		/* t1 */\n\t"\
			"vaddpd	     (%%rax),%%ymm4,%%ymm4		/* t5 */\n\t"\
			"vaddpd	0x020(%%rcx),%%ymm1,%%ymm1		/* t2 */\n\t"\
			"vaddpd	0x020(%%rax),%%ymm5,%%ymm5		/* t6 */\n\t"\
			"vsubpd	     (%%rcx),%%ymm2,%%ymm2		/* t3 */\n\t"\
			"vsubpd	     (%%rax),%%ymm6,%%ymm6		/* t7 */\n\t"\
			"vsubpd	0x020(%%rcx),%%ymm3,%%ymm3		/* t4 */\n\t"\
			"vsubpd	0x020(%%rax),%%ymm7,%%ymm7		/* t8 */\n\t"\
			"/* Finish radix-4 butterfly and store results into output-array slots: */\n\t"\
			"/* Swap _o1 <-> _o3: 0x040,0x060 <-> 0x0c0,0x0e0: */\n\t"\
			"vsubpd	%%ymm4,%%ymm0,%%ymm0			/* ~t5 <- t1 -t5 */\n\t"\
			"vsubpd	%%ymm7,%%ymm2,%%ymm2			/* ~t7 <- t3 -t8 */\n\t"\
			"vsubpd	%%ymm5,%%ymm1,%%ymm1			/* ~t6 <- t2 -t6 */\n\t"\
			"vsubpd	%%ymm6,%%ymm3,%%ymm3			/* ~t4 <- t4 -t7 */\n\t"\
			"vmovaps	%%ymm0,0x080(%%rsi)		/* <- ~t5 */\n\t"\
			"vmovaps	%%ymm2,0x040(%%rsi)		/* <- ~t7 */\n\t"\
			"vmovaps	%%ymm1,0x0a0(%%rsi)		/* <- ~t6 */\n\t"\
			"vmovaps	%%ymm3,0x0e0(%%rsi)		/* <- ~t4 */\n\t"\
			"vaddpd	%%ymm4,%%ymm4,%%ymm4			/*          2*t5 */\n\t"\
			"vaddpd	%%ymm7,%%ymm7,%%ymm7			/*          2*t8 */\n\t"\
			"vaddpd	%%ymm5,%%ymm5,%%ymm5			/*          2*t6 */\n\t"\
			"vaddpd	%%ymm6,%%ymm6,%%ymm6			/*          2*t7 */\n\t"\
			"vaddpd	%%ymm0,%%ymm4,%%ymm4			/* ~t1 <- t1 +t5 */\n\t"\
			"vaddpd	%%ymm2,%%ymm7,%%ymm7			/* ~t3 <- t3 +t8 */\n\t"\
			"vaddpd	%%ymm1,%%ymm5,%%ymm5			/* ~t2 <- t2 +t6 */\n\t"\
			"vaddpd	%%ymm3,%%ymm6,%%ymm6			/* ~t8 <- t4 +t7 */\n\t"\
			"vmovaps	%%ymm4,      (%%rsi)		/* <- ~t1 */\n\t"\
			"vmovaps	%%ymm7,0x0c0(%%rsi)		/* <- ~t3 */\n\t"\
			"vmovaps	%%ymm5,0x020(%%rsi)		/* <- ~t2 */\n\t"\
			"vmovaps	%%ymm6,0x060(%%rsi)		/* <- ~t8 */\n\t"\
		"/* SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_E(add24+p[1,0,2,3], s1p08r,s1p11r,s1p10r,s1p09r): */\n\t"\
			"addq	$0x100,%%rsi	/* s1p08r */\n\t"\
			"subq	%%rax,%%rbx		/* p01 */\n\t"\
			"subq	%%rax,%%rcx		/* p02 */\n\t"\
			"subq	%%rax,%%rdx		/* p03 */\n\t"\
			"subq	%%rdi,%%rax		/* add0 */\n\t"\
			"movslq	%[__p24],%%rdi\n\t"\
			"shlq	$3,%%rdi		/* Pointer offset for floating doubles */\n\t"\
			"addq	%%rdi,%%rax		/* add24 */\n\t"\
			"addq	%%rax,%%rbx\n\t"\
			"addq	%%rax,%%rcx\n\t"\
			"addq	%%rax,%%rdx\n\t"\
			"/* eax <-> ebx */\n\t"\
			"vmovaps	     (%%rbx),%%ymm0		/* a(jt   ) */\n\t"\
			"vmovaps	     (%%rcx),%%ymm4		/* a(jt+p2) */\n\t"\
			"vmovaps	0x020(%%rbx),%%ymm1		/* a(jp   ) */\n\t"\
			"vmovaps	0x020(%%rcx),%%ymm5		/* a(jp+p2) */\n\t"\
			"vmovaps	%%ymm0,%%ymm2			/* ymm2 <- cpy a(jt   ) */\n\t"\
			"vmovaps	%%ymm4,%%ymm6			/* ymm4 <- cpy a(jt+p2) */\n\t"\
			"vmovaps	%%ymm1,%%ymm3			/* ymm3 <- cpy a(jp   ) */\n\t"\
			"vmovaps	%%ymm5,%%ymm7			/* ymm5 <- cpy a(jp+p2) */\n\t"\
			"vaddpd	     (%%rax),%%ymm0,%%ymm0		/* t1 */\n\t"\
			"vaddpd	     (%%rdx),%%ymm4,%%ymm4		/* t5 */\n\t"\
			"vaddpd	0x020(%%rax),%%ymm1,%%ymm1		/* t2 */\n\t"\
			"vaddpd	0x020(%%rdx),%%ymm5,%%ymm5		/* t6 */\n\t"\
			"vsubpd	     (%%rax),%%ymm2,%%ymm2		/* t3 */\n\t"\
			"vsubpd	     (%%rdx),%%ymm6,%%ymm6		/* t7 */\n\t"\
			"vsubpd	0x020(%%rax),%%ymm3,%%ymm3		/* t4 */\n\t"\
			"vsubpd	0x020(%%rdx),%%ymm7,%%ymm7		/* t8 */\n\t"\
			"/* Finish radix-4 butterfly and store results into output-array slots: */\n\t"\
			"/* Swap _o1 <-> _o3: 0x040,0x060 <-> 0x0c0,0x0e0: */\n\t"\
			"vsubpd	%%ymm4,%%ymm0,%%ymm0			/* ~t5 <- t1 -t5 */\n\t"\
			"vsubpd	%%ymm7,%%ymm2,%%ymm2			/* ~t7 <- t3 -t8 */\n\t"\
			"vsubpd	%%ymm5,%%ymm1,%%ymm1			/* ~t6 <- t2 -t6 */\n\t"\
			"vsubpd	%%ymm6,%%ymm3,%%ymm3			/* ~t4 <- t4 -t7 */\n\t"\
			"vmovaps	%%ymm0,0x080(%%rsi)		/* <- ~t5 */\n\t"\
			"vmovaps	%%ymm2,0x040(%%rsi)		/* <- ~t7 */\n\t"\
			"vmovaps	%%ymm1,0x0a0(%%rsi)		/* <- ~t6 */\n\t"\
			"vmovaps	%%ymm3,0x0e0(%%rsi)		/* <- ~t4 */\n\t"\
			"vaddpd	%%ymm4,%%ymm4,%%ymm4			/*          2*t5 */\n\t"\
			"vaddpd	%%ymm7,%%ymm7,%%ymm7			/*          2*t8 */\n\t"\
			"vaddpd	%%ymm5,%%ymm5,%%ymm5			/*          2*t6 */\n\t"\
			"vaddpd	%%ymm6,%%ymm6,%%ymm6			/*          2*t7 */\n\t"\
			"vaddpd	%%ymm0,%%ymm4,%%ymm4			/* ~t1 <- t1 +t5 */\n\t"\
			"vaddpd	%%ymm2,%%ymm7,%%ymm7			/* ~t3 <- t3 +t8 */\n\t"\
			"vaddpd	%%ymm1,%%ymm5,%%ymm5			/* ~t2 <- t2 +t6 */\n\t"\
			"vaddpd	%%ymm3,%%ymm6,%%ymm6			/* ~t8 <- t4 +t7 */\n\t"\
			"vmovaps	%%ymm4,      (%%rsi)		/* <- ~t1 */\n\t"\
			"vmovaps	%%ymm7,0x0c0(%%rsi)		/* <- ~t3 */\n\t"\
			"vmovaps	%%ymm5,0x020(%%rsi)		/* <- ~t2 */\n\t"\
			"vmovaps	%%ymm6,0x060(%%rsi)		/* <- ~t8 */\n\t"\
		"/* SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_E(add08+p[1,0,2,3], s1p12r,s1p15r,s1p14r,s1p13r): */\n\t"\
			"addq	$0x100,%%rsi	/* s1p12r */\n\t"\
			"subq	%%rax,%%rbx		/* p01 */\n\t"\
			"subq	%%rax,%%rcx		/* p02 */\n\t"\
			"subq	%%rax,%%rdx		/* p03 */\n\t"\
			"subq	%%rdi,%%rax		/* add0 */\n\t"\
			"movslq	%[__p08],%%rdi\n\t"\
			"shlq	$3,%%rdi		/* Pointer offset for floating doubles */\n\t"\
			"addq	%%rdi,%%rax		/* add08 */\n\t"\
			"addq	%%rax,%%rbx\n\t"\
			"addq	%%rax,%%rcx\n\t"\
			"addq	%%rax,%%rdx\n\t"\
			"/* eax <-> %%rbx */\n\t"\
			"vmovaps	     (%%rbx),%%ymm0		/* a(jt   ) */\n\t"\
			"vmovaps	     (%%rcx),%%ymm4		/* a(jt+p2) */\n\t"\
			"vmovaps	0x020(%%rbx),%%ymm1		/* a(jp   ) */\n\t"\
			"vmovaps	0x020(%%rcx),%%ymm5		/* a(jp+p2) */\n\t"\
			"vmovaps	%%ymm0,%%ymm2			/* ymm2 <- cpy a(jt   ) */\n\t"\
			"vmovaps	%%ymm4,%%ymm6			/* ymm4 <- cpy a(jt+p2) */\n\t"\
			"vmovaps	%%ymm1,%%ymm3			/* ymm3 <- cpy a(jp   ) */\n\t"\
			"vmovaps	%%ymm5,%%ymm7			/* ymm5 <- cpy a(jp+p2) */\n\t"\
			"vaddpd	     (%%rax),%%ymm0,%%ymm0		/* t1 */\n\t"\
			"vaddpd	     (%%rdx),%%ymm4,%%ymm4		/* t5 */\n\t"\
			"vaddpd	0x020(%%rax),%%ymm1,%%ymm1		/* t2 */\n\t"\
			"vaddpd	0x020(%%rdx),%%ymm5,%%ymm5		/* t6 */\n\t"\
			"vsubpd	     (%%rax),%%ymm2,%%ymm2		/* t3 */\n\t"\
			"vsubpd	     (%%rdx),%%ymm6,%%ymm6		/* t7 */\n\t"\
			"vsubpd	0x020(%%rax),%%ymm3,%%ymm3		/* t4 */\n\t"\
			"vsubpd	0x020(%%rdx),%%ymm7,%%ymm7		/* t8 */\n\t"\
			"/* Finish radix-4 butterfly and store results into output-array slots: */\n\t"\
			"/* Swap _o1 <-> _o3: 0x040,0x060 <-> 0x0c0,0x0e0: */\n\t"\
			"vsubpd	%%ymm4,%%ymm0,%%ymm0			/* ~t5 <- t1 -t5 */\n\t"\
			"vsubpd	%%ymm7,%%ymm2,%%ymm2			/* ~t7 <- t3 -t8 */\n\t"\
			"vsubpd	%%ymm5,%%ymm1,%%ymm1			/* ~t6 <- t2 -t6 */\n\t"\
			"vsubpd	%%ymm6,%%ymm3,%%ymm3			/* ~t4 <- t4 -t7 */\n\t"\
			"vmovaps	%%ymm0,0x080(%%rsi)		/* <- ~t5 */\n\t"\
			"vmovaps	%%ymm2,0x040(%%rsi)		/* <- ~t7 */\n\t"\
			"vmovaps	%%ymm1,0x0a0(%%rsi)		/* <- ~t6 */\n\t"\
			"vmovaps	%%ymm3,0x0e0(%%rsi)		/* <- ~t4 */\n\t"\
			"vaddpd	%%ymm4,%%ymm4,%%ymm4			/*          2*t5 */\n\t"\
			"vaddpd	%%ymm7,%%ymm7,%%ymm7			/*          2*t8 */\n\t"\
			"vaddpd	%%ymm5,%%ymm5,%%ymm5			/*          2*t6 */\n\t"\
			"vaddpd	%%ymm6,%%ymm6,%%ymm6			/*          2*t7 */\n\t"\
			"vaddpd	%%ymm0,%%ymm4,%%ymm4			/* ~t1 <- t1 +t5 */\n\t"\
			"vaddpd	%%ymm2,%%ymm7,%%ymm7			/* ~t3 <- t3 +t8 */\n\t"\
			"vaddpd	%%ymm1,%%ymm5,%%ymm5			/* ~t2 <- t2 +t6 */\n\t"\
			"vaddpd	%%ymm3,%%ymm6,%%ymm6			/* ~t8 <- t4 +t7 */\n\t"\
			"vmovaps	%%ymm4,      (%%rsi)		/* <- ~t1 */\n\t"\
			"vmovaps	%%ymm7,0x0c0(%%rsi)		/* <- ~t3 */\n\t"\
			"vmovaps	%%ymm5,0x020(%%rsi)		/* <- ~t2 */\n\t"\
			"vmovaps	%%ymm6,0x060(%%rsi)		/* <- ~t8 */\n\t"\
		"/* SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_E(add20+p[2,3,0,1], s1p16r,s1p19r,s1p18r,s1p17r): */\n\t"\
			"addq	$0x100,%%rsi	/* s1p16r */\n\t"\
			"subq	%%rax,%%rbx		/* p01 */\n\t"\
			"subq	%%rax,%%rcx		/* p02 */\n\t"\
			"subq	%%rax,%%rdx		/* p03 */\n\t"\
			"subq	%%rdi,%%rax		/* add0 */\n\t"\
			"movslq	%[__p20],%%rdi\n\t"\
			"shlq	$3,%%rdi		/* Pointer offset for floating doubles */\n\t"\
			"addq	%%rdi,%%rax		/* add20 */\n\t"\
			"addq	%%rax,%%rbx\n\t"\
			"addq	%%rax,%%rcx\n\t"\
			"addq	%%rax,%%rdx\n\t"\
			"/* eax <-> %%rcx, %%rbx <-> edx */\n\t"\
			"vmovaps	     (%%rcx),%%ymm0		/* a(jt   ) */\n\t"\
			"vmovaps	     (%%rax),%%ymm4		/* a(jt+p2) */\n\t"\
			"vmovaps	0x020(%%rcx),%%ymm1		/* a(jp   ) */\n\t"\
			"vmovaps	0x020(%%rax),%%ymm5		/* a(jp+p2) */\n\t"\
			"vmovaps	%%ymm0,%%ymm2			/* ymm2 <- cpy a(jt   ) */\n\t"\
			"vmovaps	%%ymm4,%%ymm6			/* ymm4 <- cpy a(jt+p2) */\n\t"\
			"vmovaps	%%ymm1,%%ymm3			/* ymm3 <- cpy a(jp   ) */\n\t"\
			"vmovaps	%%ymm5,%%ymm7			/* ymm5 <- cpy a(jp+p2) */\n\t"\
			"vaddpd	     (%%rdx),%%ymm0,%%ymm0		/* t1 */\n\t"\
			"vaddpd	     (%%rbx),%%ymm4,%%ymm4		/* t5 */\n\t"\
			"vaddpd	0x020(%%rdx),%%ymm1,%%ymm1		/* t2 */\n\t"\
			"vaddpd	0x020(%%rbx),%%ymm5,%%ymm5		/* t6 */\n\t"\
			"vsubpd	     (%%rdx),%%ymm2,%%ymm2		/* t3 */\n\t"\
			"vsubpd	     (%%rbx),%%ymm6,%%ymm6		/* t7 */\n\t"\
			"vsubpd	0x020(%%rdx),%%ymm3,%%ymm3		/* t4 */\n\t"\
			"vsubpd	0x020(%%rbx),%%ymm7,%%ymm7		/* t8 */\n\t"\
			"/* Finish radix-4 butterfly and store results into output-array slots: */\n\t"\
			"/* Swap _o1 <-> _o3: 0x040,0x060 <-> 0x0c0,0x0e0: */\n\t"\
			"vsubpd	%%ymm4,%%ymm0,%%ymm0			/* ~t5 <- t1 -t5 */\n\t"\
			"vsubpd	%%ymm7,%%ymm2,%%ymm2			/* ~t7 <- t3 -t8 */\n\t"\
			"vsubpd	%%ymm5,%%ymm1,%%ymm1			/* ~t6 <- t2 -t6 */\n\t"\
			"vsubpd	%%ymm6,%%ymm3,%%ymm3			/* ~t4 <- t4 -t7 */\n\t"\
			"vmovaps	%%ymm0,0x080(%%rsi)		/* <- ~t5 */\n\t"\
			"vmovaps	%%ymm2,0x040(%%rsi)		/* <- ~t7 */\n\t"\
			"vmovaps	%%ymm1,0x0a0(%%rsi)		/* <- ~t6 */\n\t"\
			"vmovaps	%%ymm3,0x0e0(%%rsi)		/* <- ~t4 */\n\t"\
			"vaddpd	%%ymm4,%%ymm4,%%ymm4			/*          2*t5 */\n\t"\
			"vaddpd	%%ymm7,%%ymm7,%%ymm7			/*          2*t8 */\n\t"\
			"vaddpd	%%ymm5,%%ymm5,%%ymm5			/*          2*t6 */\n\t"\
			"vaddpd	%%ymm6,%%ymm6,%%ymm6			/*          2*t7 */\n\t"\
			"vaddpd	%%ymm0,%%ymm4,%%ymm4			/* ~t1 <- t1 +t5 */\n\t"\
			"vaddpd	%%ymm2,%%ymm7,%%ymm7			/* ~t3 <- t3 +t8 */\n\t"\
			"vaddpd	%%ymm1,%%ymm5,%%ymm5			/* ~t2 <- t2 +t6 */\n\t"\
			"vaddpd	%%ymm3,%%ymm6,%%ymm6			/* ~t8 <- t4 +t7 */\n\t"\
			"vmovaps	%%ymm4,      (%%rsi)		/* <- ~t1 */\n\t"\
			"vmovaps	%%ymm7,0x0c0(%%rsi)		/* <- ~t3 */\n\t"\
			"vmovaps	%%ymm5,0x020(%%rsi)		/* <- ~t2 */\n\t"\
			"vmovaps	%%ymm6,0x060(%%rsi)		/* <- ~t8 */\n\t"\
		"/* SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_E(add04+p[2,3,0,1], s1p20r,s1p23r,s1p22r,s1p21r): */\n\t"\
			"addq	$0x100,%%rsi	/* s1p20r */\n\t"\
			"subq	%%rax,%%rbx		/* p01 */\n\t"\
			"subq	%%rax,%%rcx		/* p02 */\n\t"\
			"subq	%%rax,%%rdx		/* p03 */\n\t"\
			"subq	%%rdi,%%rax		/* add0 */\n\t"\
			"movslq	%[__p04],%%rdi\n\t"\
			"shlq	$3,%%rdi		/* Pointer offset for floating doubles */\n\t"\
			"addq	%%rdi,%%rax		/* add04 */\n\t"\
			"addq	%%rax,%%rbx\n\t"\
			"addq	%%rax,%%rcx\n\t"\
			"addq	%%rax,%%rdx\n\t"\
			"/* eax <-> %%rcx, %%rbx <-> edx */\n\t"\
			"vmovaps	     (%%rcx),%%ymm0		/* a(jt   ) */\n\t"\
			"vmovaps	     (%%rax),%%ymm4		/* a(jt+p2) */\n\t"\
			"vmovaps	0x020(%%rcx),%%ymm1		/* a(jp   ) */\n\t"\
			"vmovaps	0x020(%%rax),%%ymm5		/* a(jp+p2) */\n\t"\
			"vmovaps	%%ymm0,%%ymm2			/* ymm2 <- cpy a(jt   ) */\n\t"\
			"vmovaps	%%ymm4,%%ymm6			/* ymm4 <- cpy a(jt+p2) */\n\t"\
			"vmovaps	%%ymm1,%%ymm3			/* ymm3 <- cpy a(jp   ) */\n\t"\
			"vmovaps	%%ymm5,%%ymm7			/* ymm5 <- cpy a(jp+p2) */\n\t"\
			"vaddpd	     (%%rdx),%%ymm0,%%ymm0		/* t1 */\n\t"\
			"vaddpd	     (%%rbx),%%ymm4,%%ymm4		/* t5 */\n\t"\
			"vaddpd	0x020(%%rdx),%%ymm1,%%ymm1		/* t2 */\n\t"\
			"vaddpd	0x020(%%rbx),%%ymm5,%%ymm5		/* t6 */\n\t"\
			"vsubpd	     (%%rdx),%%ymm2,%%ymm2		/* t3 */\n\t"\
			"vsubpd	     (%%rbx),%%ymm6,%%ymm6		/* t7 */\n\t"\
			"vsubpd	0x020(%%rdx),%%ymm3,%%ymm3		/* t4 */\n\t"\
			"vsubpd	0x020(%%rbx),%%ymm7,%%ymm7		/* t8 */\n\t"\
			"/* Finish radix-4 butterfly and store results into output-array slots: */\n\t"\
			"/* Swap _o1 <-> _o3: 0x040,0x060 <-> 0x0c0,0x0e0: */\n\t"\
			"vsubpd	%%ymm4,%%ymm0,%%ymm0			/* ~t5 <- t1 -t5 */\n\t"\
			"vsubpd	%%ymm7,%%ymm2,%%ymm2			/* ~t7 <- t3 -t8 */\n\t"\
			"vsubpd	%%ymm5,%%ymm1,%%ymm1			/* ~t6 <- t2 -t6 */\n\t"\
			"vsubpd	%%ymm6,%%ymm3,%%ymm3			/* ~t4 <- t4 -t7 */\n\t"\
			"vmovaps	%%ymm0,0x080(%%rsi)		/* <- ~t5 */\n\t"\
			"vmovaps	%%ymm2,0x040(%%rsi)		/* <- ~t7 */\n\t"\
			"vmovaps	%%ymm1,0x0a0(%%rsi)		/* <- ~t6 */\n\t"\
			"vmovaps	%%ymm3,0x0e0(%%rsi)		/* <- ~t4 */\n\t"\
			"vaddpd	%%ymm4,%%ymm4,%%ymm4			/*          2*t5 */\n\t"\
			"vaddpd	%%ymm7,%%ymm7,%%ymm7			/*          2*t8 */\n\t"\
			"vaddpd	%%ymm5,%%ymm5,%%ymm5			/*          2*t6 */\n\t"\
			"vaddpd	%%ymm6,%%ymm6,%%ymm6			/*          2*t7 */\n\t"\
			"vaddpd	%%ymm0,%%ymm4,%%ymm4			/* ~t1 <- t1 +t5 */\n\t"\
			"vaddpd	%%ymm2,%%ymm7,%%ymm7			/* ~t3 <- t3 +t8 */\n\t"\
			"vaddpd	%%ymm1,%%ymm5,%%ymm5			/* ~t2 <- t2 +t6 */\n\t"\
			"vaddpd	%%ymm3,%%ymm6,%%ymm6			/* ~t8 <- t4 +t7 */\n\t"\
			"vmovaps	%%ymm4,      (%%rsi)		/* <- ~t1 */\n\t"\
			"vmovaps	%%ymm7,0x0c0(%%rsi)		/* <- ~t3 */\n\t"\
			"vmovaps	%%ymm5,0x020(%%rsi)		/* <- ~t2 */\n\t"\
			"vmovaps	%%ymm6,0x060(%%rsi)		/* <- ~t8 */\n\t"\
		"/* SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_E(add16+p[0,1,3,2], s1p24r,s1p27r,s1p26r,s1p25r): */\n\t"\
			"addq	$0x100,%%rsi	/* s1p24r */\n\t"\
			"subq	%%rax,%%rbx		/* p01 */\n\t"\
			"subq	%%rax,%%rcx		/* p02 */\n\t"\
			"subq	%%rax,%%rdx		/* p03 */\n\t"\
			"subq	%%rdi,%%rax		/* add0 */\n\t"\
			"movslq	%[__p16],%%rdi\n\t"\
			"shlq	$3,%%rdi		/* Pointer offset for floating doubles */\n\t"\
			"addq	%%rdi,%%rax		/* add16 */\n\t"\
			"addq	%%rax,%%rbx\n\t"\
			"addq	%%rax,%%rcx\n\t"\
			"addq	%%rax,%%rdx\n\t"\
			"/* ecx <-> %%rdx */\n\t"\
			"vmovaps	     (%%rax),%%ymm0		/* a(jt   ) */\n\t"\
			"vmovaps	     (%%rdx),%%ymm4		/* a(jt+p2) */\n\t"\
			"vmovaps	0x020(%%rax),%%ymm1		/* a(jp   ) */\n\t"\
			"vmovaps	0x020(%%rdx),%%ymm5		/* a(jp+p2) */\n\t"\
			"vmovaps	%%ymm0,%%ymm2			/* ymm2 <- cpy a(jt   ) */\n\t"\
			"vmovaps	%%ymm4,%%ymm6			/* ymm4 <- cpy a(jt+p2) */\n\t"\
			"vmovaps	%%ymm1,%%ymm3			/* ymm3 <- cpy a(jp   ) */\n\t"\
			"vmovaps	%%ymm5,%%ymm7			/* ymm5 <- cpy a(jp+p2) */\n\t"\
			"vaddpd	     (%%rbx),%%ymm0,%%ymm0		/* t1 */\n\t"\
			"vaddpd	     (%%rcx),%%ymm4,%%ymm4		/* t5 */\n\t"\
			"vaddpd	0x020(%%rbx),%%ymm1,%%ymm1		/* t2 */\n\t"\
			"vaddpd	0x020(%%rcx),%%ymm5,%%ymm5		/* t6 */\n\t"\
			"vsubpd	     (%%rbx),%%ymm2,%%ymm2		/* t3 */\n\t"\
			"vsubpd	     (%%rcx),%%ymm6,%%ymm6		/* t7 */\n\t"\
			"vsubpd	0x020(%%rbx),%%ymm3,%%ymm3		/* t4 */\n\t"\
			"vsubpd	0x020(%%rcx),%%ymm7,%%ymm7		/* t8 */\n\t"\
			"/* Finish radix-4 butterfly and store results into output-array slots: */\n\t"\
			"/* Swap _o1 <-> _o3: 0x040,0x060 <-> 0x0c0,0x0e0: */\n\t"\
			"vsubpd	%%ymm4,%%ymm0,%%ymm0			/* ~t5 <- t1 -t5 */\n\t"\
			"vsubpd	%%ymm7,%%ymm2,%%ymm2			/* ~t7 <- t3 -t8 */\n\t"\
			"vsubpd	%%ymm5,%%ymm1,%%ymm1			/* ~t6 <- t2 -t6 */\n\t"\
			"vsubpd	%%ymm6,%%ymm3,%%ymm3			/* ~t4 <- t4 -t7 */\n\t"\
			"vmovaps	%%ymm0,0x080(%%rsi)		/* <- ~t5 */\n\t"\
			"vmovaps	%%ymm2,0x040(%%rsi)		/* <- ~t7 */\n\t"\
			"vmovaps	%%ymm1,0x0a0(%%rsi)		/* <- ~t6 */\n\t"\
			"vmovaps	%%ymm3,0x0e0(%%rsi)		/* <- ~t4 */\n\t"\
			"vaddpd	%%ymm4,%%ymm4,%%ymm4			/*          2*t5 */\n\t"\
			"vaddpd	%%ymm7,%%ymm7,%%ymm7			/*          2*t8 */\n\t"\
			"vaddpd	%%ymm5,%%ymm5,%%ymm5			/*          2*t6 */\n\t"\
			"vaddpd	%%ymm6,%%ymm6,%%ymm6			/*          2*t7 */\n\t"\
			"vaddpd	%%ymm0,%%ymm4,%%ymm4			/* ~t1 <- t1 +t5 */\n\t"\
			"vaddpd	%%ymm2,%%ymm7,%%ymm7			/* ~t3 <- t3 +t8 */\n\t"\
			"vaddpd	%%ymm1,%%ymm5,%%ymm5			/* ~t2 <- t2 +t6 */\n\t"\
			"vaddpd	%%ymm3,%%ymm6,%%ymm6			/* ~t8 <- t4 +t7 */\n\t"\
			"vmovaps	%%ymm4,      (%%rsi)		/* <- ~t1 */\n\t"\
			"vmovaps	%%ymm7,0x0c0(%%rsi)		/* <- ~t3 */\n\t"\
			"vmovaps	%%ymm5,0x020(%%rsi)		/* <- ~t2 */\n\t"\
			"vmovaps	%%ymm6,0x060(%%rsi)		/* <- ~t8 */\n\t"\
			"\n\t"\
			"/*...and now do 4 radix-7 transforms...*/\n\t"\
			"\n\t"\
		"/* SSE2_RADIX_07_DFT(s1p00r,s1p04r,s1p08r,s1p12r,s1p16r,s1p20r,s1p24r,cc0,s1p00r,s1p08r,s1p16r,s1p24r,s1p04r,s1p12r,s1p20r): */\n\t"\
			"/*\n\t"\
			"t1r=A1r+A6r;	t2r=A2r+A5r;	t3r=A3r+A4r;\n\t"\
			"t6r=A1r-A6r;	t5r=A2r-A5r;	t4r=A3r-A4r;\n\t"\
			"*/\n\t"\
			"movq	%[__out],%%rdi	/* s1p00r */\n\t"\
			"movq	%[__cc0],%%rax	/* cc0 */\n\t"\
			"movq	$0x040,%%rbx\n\t"\
			"movq	$0x080,%%rcx\n\t"\
			"movq	$0x0c0,%%rdx\n\t"\
			"movq	$0x100,%%rsi\n\t"\
			"addq	%%rax,%%rbx		/* cc1 */\n\t"\
			"addq	%%rax,%%rcx		/* cc2 */\n\t"\
			"addq	%%rax,%%rdx		/* cc3 */\n\t"\
			"addq	%%rax,%%rsi		/* scratch storage */\n\t"\
			"vmovaps	0x100(%%rdi),%%ymm6			/* A1r = s1p04r */\n\t"\
			"vmovaps	0x600(%%rdi),%%ymm1			/* A6r = s1p24r */\n\t"\
			"vmovaps	0x200(%%rdi),%%ymm5			/* A2r = s1p08r */\n\t"\
			"vmovaps	0x500(%%rdi),%%ymm2			/* A5r = s1p20r */\n\t"\
			"vmovaps	0x300(%%rdi),%%ymm4			/* A3r = s1p12r */\n\t"\
			"vmovaps	0x400(%%rdi),%%ymm3			/* A4r = s1p16r */\n\t"\
			"\n\t"\
			"vsubpd	%%ymm1,%%ymm6,%%ymm6			/* t6r = A1r-A6r */\n\t"\
			"vaddpd	%%ymm1,%%ymm1,%%ymm1			/*         2*A6r */\n\t"\
			"vaddpd	%%ymm6,%%ymm1,%%ymm1			/* t1r = A1r+A6r */\n\t"\
			"\n\t"\
			"vsubpd	%%ymm2,%%ymm5,%%ymm5			/* t5r = A2r-A5r */\n\t"\
			"vaddpd	%%ymm2,%%ymm2,%%ymm2			/*         2*A5r */\n\t"\
			"vaddpd	%%ymm5,%%ymm2,%%ymm2			/* t2r = A2r+A5r */\n\t"\
			"\n\t"\
			"vmovaps	     (%%rdi),%%ymm0		/* Ar0 */\n\t"\
			"vsubpd	%%ymm3,%%ymm4,%%ymm4			/* t4r = A3r-A4r */\n\t"\
			"vaddpd	%%ymm3,%%ymm3,%%ymm3			/*         2*A4r */\n\t"\
			"vaddpd	%%ymm4,%%ymm3,%%ymm3			/* t3r = A3r+A4r */\n\t"\
			"/*\n\t"\
			"	rt  = t1r+t2r+t3r;			\n\t"\
			"	B0r = rt + A0r;				\n\t"\
			"	t0r = rt*cx0 + A0r;			t3r=(t6r-t4r+t5r)*sx0;\n\t"\
			"	t1r = t1r-t2r;				t6r= t6r-t5r;\n\t"\
			"	t2r = t3r-t2r;				t5r= t4r+t5r;\n\t"\
			"	t3r =(t1r+t2r)*cx3;			t4r=(t5r-t6r)*sx3;\n\t"\
			"	t1r = t1r*cx1;				t6r= t6r*sx1;\n\t"\
			"	t2r = t2r*cx2;				t5r= t5r*sx2;\n\t"\
			"	tt  = t1r-t3r;				t6r= t4r+t6r;\n\t"\
			"	t2r = t2r-t3r;				t5r= t4r-t5r;\n\t"\
			"	\n\t"\
			"	t1r= t0r- tt-t2r;			t4r= t3r-t6r-t5r;\n\t"\
			"	t2r= t0r+t2r;				t5r= t3r+t5r;\n\t"\
			"	t0r= t0r+ tt;				t3r= t3r+t6r;\n\t"\
			"*/\n\t"\
			"vmovaps	%%ymm0,     (%%rsi)		/* cpy t0 into scratch sincos slot */\n\t"\
			"vmovaps	%%ymm6,0x020(%%rsi)		/* cpy t6 into scratch sincos slot */\n\t"\
			"vaddpd	      %%ymm1,%%ymm0,%%ymm0		/*~A0 = A0+t1 */\n\t"\
			"vmovaps	%%ymm5,%%ymm7			/* cpy t5 */\n\t"\
			"vaddpd	      %%ymm2,%%ymm3,%%ymm3		/*~t3 = t3+t2 */\n\t"\
			"vsubpd	      %%ymm4,%%ymm5,%%ymm5		/*~t5 = t5-t4 */\n\t"\
			"vsubpd	      %%ymm2,%%ymm1,%%ymm1		/*~t1 = t1-t2 */\n\t"\
			"vsubpd	      %%ymm7,%%ymm6,%%ymm6		/*~t6 = t6-t5 */\n\t"\
			"vaddpd	      %%ymm2,%%ymm2,%%ymm2		/* 2*t2 */\n\t"\
			"vaddpd	      %%ymm7,%%ymm4,%%ymm4		/*~t5 = t4+t5 */\n\t"\
			"vaddpd	      %%ymm3,%%ymm0,%%ymm0		/* B0 */\n\t"\
			"vaddpd	0x020(%%rsi),%%ymm5,%%ymm5		/* t3 = (t5-t4)+t6 */\n\t"\
			"vsubpd	      %%ymm2,%%ymm3,%%ymm3		/*~t2 =  (t2+t3) - 2*t2 = t3-t2 */\n\t"\
			"vmovaps	%%ymm4,%%ymm7			/* cpy t5 */\n\t"\
			"vmovaps	%%ymm0,      (%%rdi)		/* <-B0 = s1p00r, ymm0 FREE */\n\t"\
			"vsubpd	      %%ymm6,%%ymm4,%%ymm4		/* t4 = ~t5-~t6 */\n\t"\
			"vmovaps	%%ymm1,%%ymm2			/* cpy ~t1 */\n\t"\
			"vsubpd	     (%%rsi),%%ymm0,%%ymm0		/* r = B0 - t0 */\n\t"\
			"vmulpd	0x020(%%rax),%%ymm5,%%ymm5		/*~t3 = t3*sx0 */\n\t"\
			"vaddpd	      %%ymm3,%%ymm2,%%ymm2		/* ~t1+~t2 */\n\t"\
			"vmulpd	     (%%rcx),%%ymm3,%%ymm3		/* t2 = t2*cx2 */\n\t"\
			"vmulpd	0x020(%%rdx),%%ymm4,%%ymm4		/*~t4 = t4*sx3 */\n\t"\
			"vmulpd	     (%%rbx),%%ymm1,%%ymm1		/* t1 = t1*cx1 */\n\t"\
			"vmulpd	0x020(%%rbx),%%ymm6,%%ymm6		/*~t6 = t6*sx1 */\n\t"\
			"vmulpd	     (%%rax),%%ymm0,%%ymm0		/* ~r = r*(cx0-1) */\n\t"\
			"vmulpd	0x020(%%rcx),%%ymm7,%%ymm7		/*~t5 = t5*sx2 */\n\t"\
			"vmulpd	     (%%rdx),%%ymm2,%%ymm2		/* t3 =(t1+t2)*cx3 */\n\t"\
			"vaddpd	     (%%rdi),%%ymm0,%%ymm0		/* t0 =~r + B0 */\n\t"\
			"vaddpd	      %%ymm4,%%ymm6,%%ymm6		/*~t6 = t4+t6 */\n\t"\
			"vsubpd	      %%ymm2,%%ymm1,%%ymm1		/* tt = t1-t3 */\n\t"\
			"vsubpd	      %%ymm7,%%ymm4,%%ymm4		/*~t5 = t4-t5, ymm7 FREE */\n\t"\
			"vsubpd	      %%ymm2,%%ymm3,%%ymm3		/* t2 = t2-t3, ymm2 FREE */\n\t"\
			"\n\t"\
			"vmovaps	%%ymm0,%%ymm2			/* cpy t0 */\n\t"\
			"vmovaps	%%ymm5,%%ymm7			/* cpy t3 */\n\t"\
			"vaddpd	      %%ymm1,%%ymm0,%%ymm0			/*~t0 = t0+tt */\n\t"\
			"vaddpd	      %%ymm6,%%ymm5,%%ymm5			/*~t3 = t3+t6 */\n\t"\
			"vaddpd	      %%ymm3,%%ymm1,%%ymm1			/*~tt = tt+t2 */\n\t"\
			"vaddpd	      %%ymm4,%%ymm6,%%ymm6			/*      t6+t5 */\n\t"\
			"vaddpd	      %%ymm2,%%ymm3,%%ymm3			/*~t2 = t2+t0 */\n\t"\
			"vaddpd	      %%ymm7,%%ymm4,%%ymm4			/*~t5 = t5+t3 */\n\t"\
			"vsubpd	      %%ymm1,%%ymm2,%%ymm2			/*~t1 = t0-tt-t2 */\n\t"\
			"vsubpd	      %%ymm6,%%ymm7,%%ymm7			/*~t4 = t3-t6-t5 */\n\t"\
			"vmovaps	%%ymm0,0x200(%%rdi)			/* B1 <- t0 = s1p08r */\n\t"\
			"vmovaps	%%ymm5,0x500(%%rdi)			/* B6 <- t3 = s1p20r */\n\t"\
			"vmovaps	%%ymm2,0x400(%%rdi)			/* B2 <- t1 = s1p16r */\n\t"\
			"vmovaps	%%ymm7,0x300(%%rdi)			/* B5 <- t4 = s1p12r */\n\t"\
			"vmovaps	%%ymm3,0x600(%%rdi)			/* B3 <- t2 = s1p24r */\n\t"\
			"vmovaps	%%ymm4,0x100(%%rdi)			/* B4 <- t5 = s1p04r */\n\t"\
			"/************************** Imaginary Parts: ******************************************/\n\t"\
			"vmovaps	0x120(%%rdi),%%ymm6		/* A1i = s1p04i */\n\t"\
			"vmovaps	0x620(%%rdi),%%ymm1		/* A6i = s1p24i */\n\t"\
			"vmovaps	0x220(%%rdi),%%ymm5		/* A2i = s1p08i */\n\t"\
			"vmovaps	0x520(%%rdi),%%ymm2		/* A5i = s1p20i */\n\t"\
			"vmovaps	0x320(%%rdi),%%ymm4		/* A3i = s1p12i */\n\t"\
			"vmovaps	0x420(%%rdi),%%ymm3		/* A4i = s1p16i */\n\t"\
			"\n\t"\
			"vsubpd	%%ymm1,%%ymm6,%%ymm6			/* t6i = A1i-A6i */\n\t"\
			"vaddpd	%%ymm1,%%ymm1,%%ymm1			/*         2*A6i */\n\t"\
			"vaddpd	%%ymm6,%%ymm1,%%ymm1			/* t1i = A1i+A6i */\n\t"\
			"\n\t"\
			"vsubpd	%%ymm2,%%ymm5,%%ymm5			/* t5i = A2i-A5i */\n\t"\
			"vaddpd	%%ymm2,%%ymm2,%%ymm2			/*         2*A5i */\n\t"\
			"vaddpd	%%ymm5,%%ymm2,%%ymm2			/* t2i = A2i+A5i */\n\t"\
			"\n\t"\
			"vmovaps	0x020(%%rdi),%%ymm0		/* Ai0 */\n\t"\
			"vsubpd	%%ymm3,%%ymm4,%%ymm4			/* t4i = A3i-A4i */\n\t"\
			"vaddpd	%%ymm3,%%ymm3,%%ymm3			/*         2*A4i */\n\t"\
			"vaddpd	%%ymm4,%%ymm3,%%ymm3			/* t3i = A3i+A4i */\n\t"\
			"/*\n\t"\
			"	it  = t1i+t2i+t3i;			\n\t"\
			"	B0i = it + A0i;				\n\t"\
			"	t0i = it*cx0 + A0i;			t3i=(t6i-t4i+t5i)*sx0;\n\t"\
			"	t1i = t1i-t2i;				t6i= t6i-t5i;\n\t"\
			"	t2i = t2i-t3i;				t5i= t4i+t5i;\n\t"\
			"	t3i =(t1i-t2i)*cx3;			t4i=(t5i-t6i)*sx3;\n\t"\
			"	t1i = t1i*cx1;				t6i= t6i*sx1;\n\t"\
			"	t2i = t2i*cx2;				t5i= t5i*sx2;\n\t"\
			"	it  = t1i-t3i;				t6i= t4i+t6i;\n\t"\
			"	t2i = t2i-t3i;				t5i= t4i-t5i;\n\t"\
			"	\n\t"\
			"	t1i= t0i- it-t2i;			t4i= t3i-t6i-t5i;\n\t"\
			"	t2i= t0i+t2i;				t5i= t3i+t5i;\n\t"\
			"	t0i= t0i+ it;				t3i= t3i+t6i;\n\t"\
			"*/\n\t"\
			"vmovaps	%%ymm0,     (%%rsi)		/* cpy t0 into scratch sincos slot */\n\t"\
			"vmovaps	%%ymm6,0x020(%%rsi)		/* cpy t6 into scratch sincos slot */\n\t"\
			"vaddpd	      %%ymm1,%%ymm0,%%ymm0			/*~A0 = A0+t1 */\n\t"\
			"vmovaps	%%ymm5,%%ymm7			/* cpy t5 */\n\t"\
			"vaddpd	      %%ymm2,%%ymm3,%%ymm3			/*~t3 = t3+t2 */\n\t"\
			"vsubpd	      %%ymm4,%%ymm5,%%ymm5			/*~t5 = t5-t4 */\n\t"\
			"vsubpd	      %%ymm2,%%ymm1,%%ymm1			/*~t1 = t1-t2 */\n\t"\
			"vsubpd	      %%ymm7,%%ymm6,%%ymm6			/*~t6 = t6-t5 */\n\t"\
			"vaddpd	      %%ymm2,%%ymm2,%%ymm2			/* 2*t2 */\n\t"\
			"vaddpd	      %%ymm7,%%ymm4,%%ymm4			/*~t5 = t4+t5 */\n\t"\
			"vaddpd	      %%ymm3,%%ymm0,%%ymm0			/* B0 */\n\t"\
			"vaddpd	0x020(%%rsi),%%ymm5,%%ymm5		/* t3 = (t5-t4)+t6 */\n\t"\
			"vsubpd	      %%ymm2,%%ymm3,%%ymm3			/*~t2 =  (t2+t3) - 2*t2 = t3-t2 */\n\t"\
			"vmovaps	%%ymm4,%%ymm7			/* cpy t5 */\n\t"\
			"vmovaps	%%ymm0,0x020(%%rdi)		/* <-B0 = s1p00i, ymm0 FREE */\n\t"\
			"vsubpd	      %%ymm6,%%ymm4,%%ymm4			/* t4 = ~t5-~t6 */\n\t"\
			"vmovaps	%%ymm1,%%ymm2			/* cpy ~t1 */\n\t"\
			"vsubpd	     (%%rsi),%%ymm0,%%ymm0		/* r = B0 - t0 */\n\t"\
			"vmulpd	0x020(%%rax),%%ymm5,%%ymm5		/*~t3 = t3*sx0 */\n\t"\
			"vaddpd	      %%ymm3,%%ymm2,%%ymm2			/* ~t1+~t2 */\n\t"\
			"vmulpd	     (%%rcx),%%ymm3,%%ymm3		/* t2 = t2*cx2 */\n\t"\
			"vmulpd	0x020(%%rdx),%%ymm4,%%ymm4		/*~t4 = t4*sx3 */\n\t"\
			"vmulpd	     (%%rbx),%%ymm1,%%ymm1		/* t1 = t1*cx1 */\n\t"\
			"vmulpd	0x020(%%rbx),%%ymm6,%%ymm6		/*~t6 = t6*sx1 */\n\t"\
			"vmulpd	     (%%rax),%%ymm0,%%ymm0		/* ~r = r*(cx0-1) */\n\t"\
			"vmulpd	0x020(%%rcx),%%ymm7,%%ymm7		/*~t5 = t5*sx2 */\n\t"\
			"vmulpd	     (%%rdx),%%ymm2,%%ymm2		/* t3 =(t1+t2)*cx3 */\n\t"\
			"vaddpd	0x020(%%rdi),%%ymm0,%%ymm0		/* t0 =~r + B0 */\n\t"\
			"vaddpd	      %%ymm4,%%ymm6,%%ymm6			/*~t6 = t4+t6 */\n\t"\
			"vsubpd	      %%ymm2,%%ymm1,%%ymm1			/* tt = t1-t3 */\n\t"\
			"vsubpd	      %%ymm7,%%ymm4,%%ymm4			/*~t5 = t4-t5, ymm7 FREE */\n\t"\
			"vsubpd	      %%ymm2,%%ymm3,%%ymm3			/* t2 = t2-t3, ymm2 FREE */\n\t"\
			"vmovaps	%%ymm0,%%ymm2			/* cpy t0 */\n\t"\
			"vmovaps	%%ymm5,%%ymm7			/* cpy t3 */\n\t"\
			"vaddpd	      %%ymm1,%%ymm0,%%ymm0			/*~t0 = t0+tt */\n\t"\
			"vaddpd	      %%ymm6,%%ymm5,%%ymm5			/*~t3 = t3+t6 */\n\t"\
			"vaddpd	      %%ymm3,%%ymm1,%%ymm1			/*~tt = tt+t2 */\n\t"\
			"vaddpd	      %%ymm4,%%ymm6,%%ymm6			/*      t6+t5 */\n\t"\
			"vaddpd	      %%ymm2,%%ymm3,%%ymm3			/*~t2 = t2+t0 */\n\t"\
			"vaddpd	      %%ymm7,%%ymm4,%%ymm4			/*~t5 = t5+t3 */\n\t"\
			"vsubpd	      %%ymm1,%%ymm2,%%ymm2			/*~t1 = t0-tt-t2, ymm1 FREE */\n\t"\
			"vsubpd	      %%ymm6,%%ymm7,%%ymm7			/*~t4 = t3-t6-t5, ymm6 FREE */\n\t"\
			"/*\n\t"\
			"	B1r =t0r-t3i;				B1i*=t0i+t3r;\n\t"\
			"	B2r =t1r-t4i;				B2i*=t1i+t4r;\n\t"\
			"	B3r*=t2r+t5i;				B3i =t2i-t5r;\n\t"\
			"	B4r*=t2r-t5i;				B4i =t2i+t5r;\n\t"\
			"	B5r =t1r+t4i;				B5i*=t1i-t4r;\n\t"\
			"	B6r =t0r+t3i;				B6i*=t0i-t3r;\n\t"\
			"*/\n\t"\
			"/* ymm1,6 FREE */\n\t"\
			"vmovaps	0x200(%%rdi),%%ymm1		/* t0r */\n\t"\
			"vmovaps	0x500(%%rdi),%%ymm6		/* t3r */\n\t"\
			"vsubpd	%%ymm5,%%ymm1,%%ymm1			/* B1r =t0r-t3i */\n\t"\
			"vsubpd	%%ymm6,%%ymm0,%%ymm0			/* B6i =t0i-t3r */\n\t"\
			"vaddpd	%%ymm5,%%ymm5,%%ymm5			/*        2*t3i */\n\t"\
			"vaddpd	%%ymm6,%%ymm6,%%ymm6			/*        2*t3r */\n\t"\
			"vaddpd	%%ymm1,%%ymm5,%%ymm5			/* B6r =t0r+t3i */\n\t"\
			"vaddpd	%%ymm0,%%ymm6,%%ymm6			/* B1i =t0i+t3r */\n\t"\
			"vmovaps	%%ymm1,0x200(%%rdi)		/* <-B1r = s1p08r */\n\t"\
			"vmovaps	%%ymm0,0x520(%%rdi)		/* <-B6i = s1p20r */\n\t"\
			"vmovaps	%%ymm5,0x500(%%rdi)		/* <-B6r = s1p20r */\n\t"\
			"vmovaps	%%ymm6,0x220(%%rdi)		/* <-B1i = s1p08i */\n\t"\
			"\n\t"\
			"vmovaps	0x400(%%rdi),%%ymm1		/* t1r */\n\t"\
			"vmovaps	0x300(%%rdi),%%ymm6		/* t4r */\n\t"\
			"vsubpd	%%ymm7,%%ymm1,%%ymm1			/* B2r =t1r-t4i */\n\t"\
			"vsubpd	%%ymm6,%%ymm2,%%ymm2			/* B5i =t1i-t4r */\n\t"\
			"vaddpd	%%ymm7,%%ymm7,%%ymm7			/*        2*t4i */\n\t"\
			"vaddpd	%%ymm6,%%ymm6,%%ymm6			/*        2*t4r */\n\t"\
			"vaddpd	%%ymm1,%%ymm7,%%ymm7			/* B5r =t1r+t4i */\n\t"\
			"vaddpd	%%ymm2,%%ymm6,%%ymm6			/* B2i =t1i+t4r */\n\t"\
			"vmovaps	%%ymm1,0x400(%%rdi)		/* <-B2r = s1p16r */\n\t"\
			"vmovaps	%%ymm2,0x320(%%rdi)		/* <-B5i = s1p12r */\n\t"\
			"vmovaps	%%ymm7,0x300(%%rdi)		/* <-B5r = s1p12r */\n\t"\
			"vmovaps	%%ymm6,0x420(%%rdi)		/* <-B2i = s1p16i */\n\t"\
			"\n\t"\
			"/* Note the order reversal on this pair of outputs: */\n\t"\
			"vmovaps	0x600(%%rdi),%%ymm0		/* t2r */\n\t"\
			"vmovaps	0x100(%%rdi),%%ymm5		/* t5r */\n\t"\
			"vsubpd	%%ymm4,%%ymm0,%%ymm0			/* B4r =t2r-t5i */\n\t"\
			"vsubpd	%%ymm5,%%ymm3,%%ymm3			/* B3i =t2i-t5r */\n\t"\
			"vaddpd	%%ymm4,%%ymm4,%%ymm4			/*        2*t5i */\n\t"\
			"vaddpd	%%ymm5,%%ymm5,%%ymm5			/*        2*t5r */\n\t"\
			"vaddpd	%%ymm0,%%ymm4,%%ymm4			/* B3r =t2r+t5i */\n\t"\
			"vaddpd	%%ymm3,%%ymm5,%%ymm5			/* B4i =t2i+t5r */\n\t"\
			"vmovaps	%%ymm0,0x100(%%rdi)		/* <-B4r = s1p04r */\n\t"\
			"vmovaps	%%ymm3,0x620(%%rdi)		/* <-B3i = s1p24r */\n\t"\
			"vmovaps	%%ymm4,0x600(%%rdi)		/* <-B3r = s1p24r */\n\t"\
			"vmovaps	%%ymm5,0x120(%%rdi)		/* <-B4i = s1p04i */\n\t"\
			"\n\t"\
		"/* SSE2_RADIX_07_DFT(s1p01r,s1p05r,s1p09r,s1p13r,s1p17r,s1p21r,s1p25r,cc0,s1p21r,s1p01r,s1p09r,s1p17r,s1p25r,s1p05r,s1p13r): */\n\t"\
			"addq	$0x040,%%rdi	/* s1p01r */\n\t"\
			"vmovaps	0x100(%%rdi),%%ymm6			/* A1r = s1p05r */\n\t"\
			"vmovaps	0x600(%%rdi),%%ymm1			/* A6r = s1p25r */\n\t"\
			"vmovaps	0x200(%%rdi),%%ymm5			/* A2r = s1p09r */\n\t"\
			"vmovaps	0x500(%%rdi),%%ymm2			/* A5r = s1p21r */\n\t"\
			"vmovaps	0x300(%%rdi),%%ymm4			/* A3r = s1p13r */\n\t"\
			"vmovaps	0x400(%%rdi),%%ymm3			/* A4r = s1p17r */\n\t"\
			"\n\t"\
			"vsubpd	%%ymm1,%%ymm6,%%ymm6			/* t6r = A1r-A6r */\n\t"\
			"vaddpd	%%ymm1,%%ymm1,%%ymm1			/*         2*A6r */\n\t"\
			"vaddpd	%%ymm6,%%ymm1,%%ymm1			/* t1r = A1r+A6r */\n\t"\
			"\n\t"\
			"vsubpd	%%ymm2,%%ymm5,%%ymm5			/* t5r = A2r-A5r */\n\t"\
			"vaddpd	%%ymm2,%%ymm2,%%ymm2			/*         2*A5r */\n\t"\
			"vaddpd	%%ymm5,%%ymm2,%%ymm2			/* t2r = A2r+A5r */\n\t"\
			"\n\t"\
			"vmovaps	     (%%rdi),%%ymm0		/* Ar0 */\n\t"\
			"vsubpd	%%ymm3,%%ymm4,%%ymm4			/* t4r = A3r-A4r */\n\t"\
			"vaddpd	%%ymm3,%%ymm3,%%ymm3			/*         2*A4r */\n\t"\
			"vaddpd	%%ymm4,%%ymm3,%%ymm3			/* t3r = A3r+A4r */\n\t"\
			"\n\t"\
			"vmovaps	%%ymm0,     (%%rsi)		/* cpy t0 into scratch sincos slot */\n\t"\
			"vmovaps	%%ymm6,0x020(%%rsi)		/* cpy t6 into scratch sincos slot */\n\t"\
			"vaddpd	      %%ymm1,%%ymm0,%%ymm0		/*~A0 = A0+t1 */\n\t"\
			"vmovaps	%%ymm5,%%ymm7			/* cpy t5 */\n\t"\
			"vaddpd	      %%ymm2,%%ymm3,%%ymm3		/*~t3 = t3+t2 */\n\t"\
			"vsubpd	      %%ymm4,%%ymm5,%%ymm5		/*~t5 = t5-t4 */\n\t"\
			"vsubpd	      %%ymm2,%%ymm1,%%ymm1		/*~t1 = t1-t2 */\n\t"\
			"vsubpd	      %%ymm7,%%ymm6,%%ymm6		/*~t6 = t6-t5 */\n\t"\
			"vaddpd	      %%ymm2,%%ymm2,%%ymm2		/* 2*t2 */\n\t"\
			"vaddpd	      %%ymm7,%%ymm4,%%ymm4		/*~t5 = t4+t5 */\n\t"\
			"vaddpd	      %%ymm3,%%ymm0,%%ymm0		/* B0 */\n\t"\
			"vaddpd	0x020(%%rsi),%%ymm5,%%ymm5		/* t3 = (t5-t4)+t6 */\n\t"\
			"vsubpd	      %%ymm2,%%ymm3,%%ymm3		/*~t2 =  (t2+t3) - 2*t2 = t3-t2 */\n\t"\
			"vmovaps	%%ymm4,%%ymm7			/* cpy t5 */\n\t"\
			"vmovaps	%%ymm0,0x500(%%rdi)		/* <-B0 = s1p00r, ymm0 FREE */\n\t"\
			"vsubpd	      %%ymm6,%%ymm4,%%ymm4		/* t4 = ~t5-~t6 */\n\t"\
			"vmovaps	%%ymm1,%%ymm2			/* cpy ~t1 */\n\t"\
			"vsubpd	     (%%rsi),%%ymm0,%%ymm0		/* r = B0 - t0 */\n\t"\
			"vmulpd	0x020(%%rax),%%ymm5,%%ymm5		/*~t3 = t3*sx0 */\n\t"\
			"vaddpd	      %%ymm3,%%ymm2,%%ymm2		/* ~t1+~t2 */\n\t"\
			"vmulpd	     (%%rcx),%%ymm3,%%ymm3		/* t2 = t2*cx2 */\n\t"\
			"vmulpd	0x020(%%rdx),%%ymm4,%%ymm4		/*~t4 = t4*sx3 */\n\t"\
			"vmulpd	     (%%rbx),%%ymm1,%%ymm1		/* t1 = t1*cx1 */\n\t"\
			"vmulpd	0x020(%%rbx),%%ymm6,%%ymm6		/*~t6 = t6*sx1 */\n\t"\
			"vmulpd	     (%%rax),%%ymm0,%%ymm0		/* ~r = r*(cx0-1) */\n\t"\
			"vmulpd	0x020(%%rcx),%%ymm7,%%ymm7		/*~t5 = t5*sx2 */\n\t"\
			"vmulpd	     (%%rdx),%%ymm2,%%ymm2		/* t3 =(t1+t2)*cx3 */\n\t"\
			"vaddpd	0x500(%%rdi),%%ymm0,%%ymm0		/* t0 =~r + B0 */\n\t"\
			"vaddpd	      %%ymm4,%%ymm6,%%ymm6		/*~t6 = t4+t6 */\n\t"\
			"vsubpd	      %%ymm2,%%ymm1,%%ymm1		/* tt = t1-t3 */\n\t"\
			"vsubpd	      %%ymm7,%%ymm4,%%ymm4		/*~t5 = t4-t5, ymm7 FREE */\n\t"\
			"vsubpd	      %%ymm2,%%ymm3,%%ymm3		/* t2 = t2-t3, ymm2 FREE */\n\t"\
			"\n\t"\
			"vmovaps	%%ymm0,%%ymm2			/* cpy t0 */\n\t"\
			"vmovaps	%%ymm5,%%ymm7			/* cpy t3 */\n\t"\
			"vaddpd	      %%ymm1,%%ymm0,%%ymm0			/*~t0 = t0+tt */\n\t"\
			"vaddpd	      %%ymm6,%%ymm5,%%ymm5			/*~t3 = t3+t6 */\n\t"\
			"vaddpd	      %%ymm3,%%ymm1,%%ymm1			/*~tt = tt+t2 */\n\t"\
			"vaddpd	      %%ymm4,%%ymm6,%%ymm6			/*      t6+t5 */\n\t"\
			"vaddpd	      %%ymm2,%%ymm3,%%ymm3			/*~t2 = t2+t0 */\n\t"\
			"vaddpd	      %%ymm7,%%ymm4,%%ymm4			/*~t5 = t5+t3 */\n\t"\
			"vsubpd	      %%ymm1,%%ymm2,%%ymm2			/*~t1 = t0-tt-t2 */\n\t"\
			"vsubpd	      %%ymm6,%%ymm7,%%ymm7			/*~t4 = t3-t6-t5 */\n\t"\
			"vmovaps	%%ymm0,     (%%rdi)			/* B1 <- t0 = s1p01r */\n\t"\
			"vmovaps	%%ymm5,0x300(%%rdi)			/* B6 <- t3 = s1p13r */\n\t"\
			"vmovaps	%%ymm2,0x200(%%rdi)			/* B2 <- t1 = s1p09r */\n\t"\
			"vmovaps	%%ymm7,0x100(%%rdi)			/* B5 <- t4 = s1p05r */\n\t"\
			"vmovaps	%%ymm3,0x400(%%rdi)			/* B3 <- t2 = s1p17r */\n\t"\
			"vmovaps	%%ymm4,0x600(%%rdi)			/* B4 <- t5 = s1p25r */\n\t"\
			"/************************** Imaginary Parts: ******************************************/\n\t"\
			"vmovaps	0x120(%%rdi),%%ymm6		/* A1i = s1p05i */\n\t"\
			"vmovaps	0x620(%%rdi),%%ymm1		/* A6i = s1p25i */\n\t"\
			"vmovaps	0x220(%%rdi),%%ymm5		/* A2i = s1p09i */\n\t"\
			"vmovaps	0x520(%%rdi),%%ymm2		/* A5i = s1p21i */\n\t"\
			"vmovaps	0x320(%%rdi),%%ymm4		/* A3i = s1p13i */\n\t"\
			"vmovaps	0x420(%%rdi),%%ymm3		/* A4i = s1p17i */\n\t"\
			"\n\t"\
			"vsubpd	%%ymm1,%%ymm6,%%ymm6			/* t6i = A1i-A6i */\n\t"\
			"vaddpd	%%ymm1,%%ymm1,%%ymm1			/*         2*A6i */\n\t"\
			"vaddpd	%%ymm6,%%ymm1,%%ymm1			/* t1i = A1i+A6i */\n\t"\
			"\n\t"\
			"vsubpd	%%ymm2,%%ymm5,%%ymm5			/* t5i = A2i-A5i */\n\t"\
			"vaddpd	%%ymm2,%%ymm2,%%ymm2			/*         2*A5i */\n\t"\
			"vaddpd	%%ymm5,%%ymm2,%%ymm2			/* t2i = A2i+A5i */\n\t"\
			"\n\t"\
			"vmovaps	0x020(%%rdi),%%ymm0		/* Ai0 */\n\t"\
			"vsubpd	%%ymm3,%%ymm4,%%ymm4			/* t4i = A3i-A4i */\n\t"\
			"vaddpd	%%ymm3,%%ymm3,%%ymm3			/*         2*A4i */\n\t"\
			"vaddpd	%%ymm4,%%ymm3,%%ymm3			/* t3i = A3i+A4i */\n\t"\
			"\n\t"\
			"vmovaps	%%ymm0,     (%%rsi)		/* cpy t0 into scratch sincos slot */\n\t"\
			"vmovaps	%%ymm6,0x020(%%rsi)		/* cpy t6 into scratch sincos slot */\n\t"\
			"vaddpd	      %%ymm1,%%ymm0,%%ymm0			/*~A0 = A0+t1 */\n\t"\
			"vmovaps	%%ymm5,%%ymm7					/* cpy t5 */\n\t"\
			"vaddpd	      %%ymm2,%%ymm3,%%ymm3			/*~t3 = t3+t2 */\n\t"\
			"vsubpd	      %%ymm4,%%ymm5,%%ymm5			/*~t5 = t5-t4 */\n\t"\
			"vsubpd	      %%ymm2,%%ymm1,%%ymm1			/*~t1 = t1-t2 */\n\t"\
			"vsubpd	      %%ymm7,%%ymm6,%%ymm6			/*~t6 = t6-t5 */\n\t"\
			"vaddpd	      %%ymm2,%%ymm2,%%ymm2			/* 2*t2 */\n\t"\
			"vaddpd	      %%ymm7,%%ymm4,%%ymm4			/*~t5 = t4+t5 */\n\t"\
			"vaddpd	      %%ymm3,%%ymm0,%%ymm0			/* B0 */\n\t"\
			"vaddpd	0x020(%%rsi),%%ymm5,%%ymm5		/* t3 = (t5-t4)+t6 */\n\t"\
			"vsubpd	      %%ymm2,%%ymm3,%%ymm3			/*~t2 =  (t2+t3) - 2*t2 = t3-t2 */\n\t"\
			"vmovaps	%%ymm4,%%ymm7					/* cpy t5 */\n\t"\
			"vmovaps	%%ymm0,0x520(%%rdi)		/* <-B0 = s1p21i, ymm0 FREE */\n\t"\
			"vsubpd	      %%ymm6,%%ymm4,%%ymm4			/* t4 = ~t5-~t6 */\n\t"\
			"vmovaps	%%ymm1,%%ymm2					/* cpy ~t1 */\n\t"\
			"vsubpd	     (%%rsi),%%ymm0,%%ymm0		/* r = B0 - t0 */\n\t"\
			"vmulpd	0x020(%%rax),%%ymm5,%%ymm5		/*~t3 = t3*sx0 */\n\t"\
			"vaddpd	      %%ymm3,%%ymm2,%%ymm2			/* ~t1+~t2 */\n\t"\
			"vmulpd	     (%%rcx),%%ymm3,%%ymm3		/* t2 = t2*cx2 */\n\t"\
			"vmulpd	0x020(%%rdx),%%ymm4,%%ymm4		/*~t4 = t4*sx3 */\n\t"\
			"vmulpd	     (%%rbx),%%ymm1,%%ymm1		/* t1 = t1*cx1 */\n\t"\
			"vmulpd	0x020(%%rbx),%%ymm6,%%ymm6		/*~t6 = t6*sx1 */\n\t"\
			"vmulpd	     (%%rax),%%ymm0,%%ymm0		/* ~r = r*(cx0-1) */\n\t"\
			"vmulpd	0x020(%%rcx),%%ymm7,%%ymm7		/*~t5 = t5*sx2 */\n\t"\
			"vmulpd	     (%%rdx),%%ymm2,%%ymm2		/* t3 =(t1+t2)*cx3 */\n\t"\
			"vaddpd	0x520(%%rdi),%%ymm0,%%ymm0		/* t0 =~r + B0 */\n\t"\
			"vaddpd	      %%ymm4,%%ymm6,%%ymm6			/*~t6 = t4+t6 */\n\t"\
			"vsubpd	      %%ymm2,%%ymm1,%%ymm1			/* tt = t1-t3 */\n\t"\
			"vsubpd	      %%ymm7,%%ymm4,%%ymm4			/*~t5 = t4-t5, ymm7 FREE */\n\t"\
			"vsubpd	      %%ymm2,%%ymm3,%%ymm3			/* t2 = t2-t3, ymm2 FREE */\n\t"\
			"vmovaps	%%ymm0,%%ymm2					/* cpy t0 */\n\t"\
			"vmovaps	%%ymm5,%%ymm7					/* cpy t3 */\n\t"\
			"vaddpd	      %%ymm1,%%ymm0,%%ymm0			/*~t0 = t0+tt */\n\t"\
			"vaddpd	      %%ymm6,%%ymm5,%%ymm5			/*~t3 = t3+t6 */\n\t"\
			"vaddpd	      %%ymm3,%%ymm1,%%ymm1			/*~tt = tt+t2 */\n\t"\
			"vaddpd	      %%ymm4,%%ymm6,%%ymm6			/*      t6+t5 */\n\t"\
			"vaddpd	      %%ymm2,%%ymm3,%%ymm3			/*~t2 = t2+t0 */\n\t"\
			"vaddpd	      %%ymm7,%%ymm4,%%ymm4			/*~t5 = t5+t3 */\n\t"\
			"vsubpd	      %%ymm1,%%ymm2,%%ymm2			/*~t1 = t0-tt-t2, ymm1 FREE */\n\t"\
			"vsubpd	      %%ymm6,%%ymm7,%%ymm7			/*~t4 = t3-t6-t5, ymm6 FREE */\n\t"\
			"\n\t"\
			"/* ymm1,6 FREE */\n\t"\
			"vmovaps	      (%%rdi),%%ymm1		/* t0r */\n\t"\
			"vmovaps	0x300(%%rdi),%%ymm6		/* t3r */\n\t"\
			"vsubpd	%%ymm5,%%ymm1,%%ymm1			/* B1r =t0r-t3i */\n\t"\
			"vsubpd	%%ymm6,%%ymm0,%%ymm0			/* B6i =t0i-t3r */\n\t"\
			"vaddpd	%%ymm5,%%ymm5,%%ymm5			/*        2*t3i */\n\t"\
			"vaddpd	%%ymm6,%%ymm6,%%ymm6			/*        2*t3r */\n\t"\
			"vaddpd	%%ymm1,%%ymm5,%%ymm5			/* B6r =t0r+t3i */\n\t"\
			"vaddpd	%%ymm0,%%ymm6,%%ymm6			/* B1i =t0i+t3r */\n\t"\
			"vmovaps	%%ymm1,      (%%rdi)		/* <-B1r = s1p08r */\n\t"\
			"vmovaps	%%ymm0,0x320(%%rdi)		/* <-B6i = s1p20r */\n\t"\
			"vmovaps	%%ymm5,0x300(%%rdi)		/* <-B6r = s1p20r */\n\t"\
			"vmovaps	%%ymm6,0x020(%%rdi)		/* <-B1i = s1p08i */\n\t"\
			"\n\t"\
			"vmovaps	0x200(%%rdi),%%ymm1		/* t1r */\n\t"\
			"vmovaps	0x100(%%rdi),%%ymm6		/* t4r */\n\t"\
			"vsubpd	%%ymm7,%%ymm1,%%ymm1			/* B2r =t1r-t4i */\n\t"\
			"vsubpd	%%ymm6,%%ymm2,%%ymm2			/* B5i =t1i-t4r */\n\t"\
			"vaddpd	%%ymm7,%%ymm7,%%ymm7			/*        2*t4i */\n\t"\
			"vaddpd	%%ymm6,%%ymm6,%%ymm6			/*        2*t4r */\n\t"\
			"vaddpd	%%ymm1,%%ymm7,%%ymm7			/* B5r =t1r+t4i */\n\t"\
			"vaddpd	%%ymm2,%%ymm6,%%ymm6			/* B2i =t1i+t4r */\n\t"\
			"vmovaps	%%ymm1,0x200(%%rdi)		/* <-B2r = s1p16r */\n\t"\
			"vmovaps	%%ymm2,0x120(%%rdi)		/* <-B5i = s1p12r */\n\t"\
			"vmovaps	%%ymm7,0x100(%%rdi)		/* <-B5r = s1p12r */\n\t"\
			"vmovaps	%%ymm6,0x220(%%rdi)		/* <-B2i = s1p16i */\n\t"\
			"\n\t"\
			"/* Note the order reversal on this pair of outputs: */\n\t"\
			"vmovaps	0x400(%%rdi),%%ymm0		/* t2r */\n\t"\
			"vmovaps	0x600(%%rdi),%%ymm5		/* t5r */\n\t"\
			"vsubpd	%%ymm4,%%ymm0,%%ymm0			/* B4r =t2r-t5i */\n\t"\
			"vsubpd	%%ymm5,%%ymm3,%%ymm3			/* B3i =t2i-t5r */\n\t"\
			"vaddpd	%%ymm4,%%ymm4,%%ymm4			/*        2*t5i */\n\t"\
			"vaddpd	%%ymm5,%%ymm5,%%ymm5			/*        2*t5r */\n\t"\
			"vaddpd	%%ymm0,%%ymm4,%%ymm4			/* B3r =t2r+t5i */\n\t"\
			"vaddpd	%%ymm3,%%ymm5,%%ymm5			/* B4i =t2i+t5r */\n\t"\
			"vmovaps	%%ymm0,0x600(%%rdi)		/* <-B4r = s1p04r */\n\t"\
			"vmovaps	%%ymm3,0x420(%%rdi)		/* <-B3i = s1p24r */\n\t"\
			"vmovaps	%%ymm4,0x400(%%rdi)		/* <-B3r = s1p24r */\n\t"\
			"vmovaps	%%ymm5,0x620(%%rdi)		/* <-B4i = s1p04i */\n\t"\
			"\n\t"\
		"/* SSE2_RADIX_07_DFT(s1p02r,s1p06r,s1p10r,s1p14r,s1p18r,s1p22r,s1p26r,cc0,s1p14r,s1p22r,s1p02r,s1p10r,s1p18r,s1p26r,s1p06r): */\n\t"\
			"addq	$0x040,%%rdi	/* s1p02r */\n\t"\
			"vmovaps	0x100(%%rdi),%%ymm6			/* A1r = s1p06r */\n\t"\
			"vmovaps	0x600(%%rdi),%%ymm1			/* A6r = s1p26r */\n\t"\
			"vmovaps	0x200(%%rdi),%%ymm5			/* A2r = s1p10r */\n\t"\
			"vmovaps	0x500(%%rdi),%%ymm2			/* A5r = s1p22r */\n\t"\
			"vmovaps	0x300(%%rdi),%%ymm4			/* A3r = s1p14r */\n\t"\
			"vmovaps	0x400(%%rdi),%%ymm3			/* A4r = s1p18r */\n\t"\
			"\n\t"\
			"vsubpd	%%ymm1,%%ymm6,%%ymm6			/* t6r = A1r-A6r */\n\t"\
			"vaddpd	%%ymm1,%%ymm1,%%ymm1			/*         2*A6r */\n\t"\
			"vaddpd	%%ymm6,%%ymm1,%%ymm1			/* t1r = A1r+A6r */\n\t"\
			"\n\t"\
			"vsubpd	%%ymm2,%%ymm5,%%ymm5			/* t5r = A2r-A5r */\n\t"\
			"vaddpd	%%ymm2,%%ymm2,%%ymm2			/*         2*A5r */\n\t"\
			"vaddpd	%%ymm5,%%ymm2,%%ymm2			/* t2r = A2r+A5r */\n\t"\
			"\n\t"\
			"vmovaps	     (%%rdi),%%ymm0		/* Ar0 */\n\t"\
			"vsubpd	%%ymm3,%%ymm4,%%ymm4			/* t4r = A3r-A4r */\n\t"\
			"vaddpd	%%ymm3,%%ymm3,%%ymm3			/*         2*A4r */\n\t"\
			"vaddpd	%%ymm4,%%ymm3,%%ymm3			/* t3r = A3r+A4r */\n\t"\
			"\n\t"\
			"vmovaps	%%ymm0,     (%%rsi)		/* cpy t0 into scratch sincos slot */\n\t"\
			"vmovaps	%%ymm6,0x020(%%rsi)		/* cpy t6 into scratch sincos slot */\n\t"\
			"vaddpd	      %%ymm1,%%ymm0,%%ymm0		/*~A0 = A0+t1 */\n\t"\
			"vmovaps	%%ymm5,%%ymm7			/* cpy t5 */\n\t"\
			"vaddpd	      %%ymm2,%%ymm3,%%ymm3		/*~t3 = t3+t2 */\n\t"\
			"vsubpd	      %%ymm4,%%ymm5,%%ymm5		/*~t5 = t5-t4 */\n\t"\
			"vsubpd	      %%ymm2,%%ymm1,%%ymm1		/*~t1 = t1-t2 */\n\t"\
			"vsubpd	      %%ymm7,%%ymm6,%%ymm6		/*~t6 = t6-t5 */\n\t"\
			"vaddpd	      %%ymm2,%%ymm2,%%ymm2		/* 2*t2 */\n\t"\
			"vaddpd	      %%ymm7,%%ymm4,%%ymm4		/*~t5 = t4+t5 */\n\t"\
			"vaddpd	      %%ymm3,%%ymm0,%%ymm0		/* B0 */\n\t"\
			"vaddpd	0x020(%%rsi),%%ymm5,%%ymm5		/* t3 = (t5-t4)+t6 */\n\t"\
			"vsubpd	      %%ymm2,%%ymm3,%%ymm3		/*~t2 =  (t2+t3) - 2*t2 = t3-t2 */\n\t"\
			"vmovaps	%%ymm4,%%ymm7			/* cpy t5 */\n\t"\
			"vmovaps	%%ymm0,0x300(%%rdi)		/* <-B0 = s1p00r, ymm0 FREE */\n\t"\
			"vsubpd	      %%ymm6,%%ymm4,%%ymm4		/* t4 = ~t5-~t6 */\n\t"\
			"vmovaps	%%ymm1,%%ymm2			/* cpy ~t1 */\n\t"\
			"vsubpd	     (%%rsi),%%ymm0,%%ymm0		/* r = B0 - t0 */\n\t"\
			"vmulpd	0x020(%%rax),%%ymm5,%%ymm5		/*~t3 = t3*sx0 */\n\t"\
			"vaddpd	      %%ymm3,%%ymm2,%%ymm2		/* ~t1+~t2 */\n\t"\
			"vmulpd	     (%%rcx),%%ymm3,%%ymm3		/* t2 = t2*cx2 */\n\t"\
			"vmulpd	0x020(%%rdx),%%ymm4,%%ymm4		/*~t4 = t4*sx3 */\n\t"\
			"vmulpd	     (%%rbx),%%ymm1,%%ymm1		/* t1 = t1*cx1 */\n\t"\
			"vmulpd	0x020(%%rbx),%%ymm6,%%ymm6		/*~t6 = t6*sx1 */\n\t"\
			"vmulpd	     (%%rax),%%ymm0,%%ymm0		/* ~r = r*(cx0-1) */\n\t"\
			"vmulpd	0x020(%%rcx),%%ymm7,%%ymm7		/*~t5 = t5*sx2 */\n\t"\
			"vmulpd	     (%%rdx),%%ymm2,%%ymm2		/* t3 =(t1+t2)*cx3 */\n\t"\
			"vaddpd	0x300(%%rdi),%%ymm0,%%ymm0		/* t0 =~r + B0 */\n\t"\
			"vaddpd	      %%ymm4,%%ymm6,%%ymm6		/*~t6 = t4+t6 */\n\t"\
			"vsubpd	      %%ymm2,%%ymm1,%%ymm1		/* tt = t1-t3 */\n\t"\
			"vsubpd	      %%ymm7,%%ymm4,%%ymm4		/*~t5 = t4-t5, ymm7 FREE */\n\t"\
			"vsubpd	      %%ymm2,%%ymm3,%%ymm3		/* t2 = t2-t3, ymm2 FREE */\n\t"\
			"\n\t"\
			"vmovaps	%%ymm0,%%ymm2			/* cpy t0 */\n\t"\
			"vmovaps	%%ymm5,%%ymm7			/* cpy t3 */\n\t"\
			"vaddpd	      %%ymm1,%%ymm0,%%ymm0			/*~t0 = t0+tt */\n\t"\
			"vaddpd	      %%ymm6,%%ymm5,%%ymm5			/*~t3 = t3+t6 */\n\t"\
			"vaddpd	      %%ymm3,%%ymm1,%%ymm1			/*~tt = tt+t2 */\n\t"\
			"vaddpd	      %%ymm4,%%ymm6,%%ymm6			/*      t6+t5 */\n\t"\
			"vaddpd	      %%ymm2,%%ymm3,%%ymm3			/*~t2 = t2+t0 */\n\t"\
			"vaddpd	      %%ymm7,%%ymm4,%%ymm4			/*~t5 = t5+t3 */\n\t"\
			"vsubpd	      %%ymm1,%%ymm2,%%ymm2			/*~t1 = t0-tt-t2 */\n\t"\
			"vsubpd	      %%ymm6,%%ymm7,%%ymm7			/*~t4 = t3-t6-t5 */\n\t"\
			"vmovaps	%%ymm0,0x500(%%rdi)			/* B1 <- t0 = s1p22r */\n\t"\
			"vmovaps	%%ymm5,0x100(%%rdi)			/* B6 <- t3 = s1p06r */\n\t"\
			"vmovaps	%%ymm2,      (%%rdi)			/* B2 <- t1 = s1p02r */\n\t"\
			"vmovaps	%%ymm7,0x600(%%rdi)			/* B5 <- t4 = s1p26r */\n\t"\
			"vmovaps	%%ymm3,0x200(%%rdi)			/* B3 <- t2 = s1p10r */\n\t"\
			"vmovaps	%%ymm4,0x400(%%rdi)			/* B4 <- t5 = s1p18r */\n\t"\
			"/************************** Imaginary Parts: ******************************************/\n\t"\
			"vmovaps	0x120(%%rdi),%%ymm6		/* A1i = s1p06i */\n\t"\
			"vmovaps	0x620(%%rdi),%%ymm1		/* A6i = s1p26i */\n\t"\
			"vmovaps	0x220(%%rdi),%%ymm5		/* A2i = s1p10i */\n\t"\
			"vmovaps	0x520(%%rdi),%%ymm2		/* A5i = s1p22i */\n\t"\
			"vmovaps	0x320(%%rdi),%%ymm4		/* A3i = s1p14i */\n\t"\
			"vmovaps	0x420(%%rdi),%%ymm3		/* A4i = s1p18i */\n\t"\
			"\n\t"\
			"vsubpd	%%ymm1,%%ymm6,%%ymm6			/* t6i = A1i-A6i */\n\t"\
			"vaddpd	%%ymm1,%%ymm1,%%ymm1			/*         2*A6i */\n\t"\
			"vaddpd	%%ymm6,%%ymm1,%%ymm1			/* t1i = A1i+A6i */\n\t"\
			"\n\t"\
			"vsubpd	%%ymm2,%%ymm5,%%ymm5			/* t5i = A2i-A5i */\n\t"\
			"vaddpd	%%ymm2,%%ymm2,%%ymm2			/*         2*A5i */\n\t"\
			"vaddpd	%%ymm5,%%ymm2,%%ymm2			/* t2i = A2i+A5i */\n\t"\
			"\n\t"\
			"vmovaps	0x020(%%rdi),%%ymm0		/* Ai0 */\n\t"\
			"vsubpd	%%ymm3,%%ymm4,%%ymm4			/* t4i = A3i-A4i */\n\t"\
			"vaddpd	%%ymm3,%%ymm3,%%ymm3			/*         2*A4i */\n\t"\
			"vaddpd	%%ymm4,%%ymm3,%%ymm3			/* t3i = A3i+A4i */\n\t"\
			"\n\t"\
			"vmovaps	%%ymm0,     (%%rsi)		/* cpy t0 into scratch sincos slot */\n\t"\
			"vmovaps	%%ymm6,0x020(%%rsi)		/* cpy t6 into scratch sincos slot */\n\t"\
			"vaddpd	      %%ymm1,%%ymm0,%%ymm0			/*~A0 = A0+t1 */\n\t"\
			"vmovaps	%%ymm5,%%ymm7					/* cpy t5 */\n\t"\
			"vaddpd	      %%ymm2,%%ymm3,%%ymm3			/*~t3 = t3+t2 */\n\t"\
			"vsubpd	      %%ymm4,%%ymm5,%%ymm5			/*~t5 = t5-t4 */\n\t"\
			"vsubpd	      %%ymm2,%%ymm1,%%ymm1			/*~t1 = t1-t2 */\n\t"\
			"vsubpd	      %%ymm7,%%ymm6,%%ymm6			/*~t6 = t6-t5 */\n\t"\
			"vaddpd	      %%ymm2,%%ymm2,%%ymm2			/* 2*t2 */\n\t"\
			"vaddpd	      %%ymm7,%%ymm4,%%ymm4			/*~t5 = t4+t5 */\n\t"\
			"vaddpd	      %%ymm3,%%ymm0,%%ymm0			/* B0 */\n\t"\
			"vaddpd	0x020(%%rsi),%%ymm5,%%ymm5		/* t3 = (t5-t4)+t6 */\n\t"\
			"vsubpd	      %%ymm2,%%ymm3,%%ymm3			/*~t2 =  (t2+t3) - 2*t2 = t3-t2 */\n\t"\
			"vmovaps	%%ymm4,%%ymm7					/* cpy t5 */\n\t"\
			"vmovaps	%%ymm0,0x320(%%rdi)		/* <-B0 = s1p14i, ymm0 FREE */\n\t"\
			"vsubpd	      %%ymm6,%%ymm4,%%ymm4			/* t4 = ~t5-~t6 */\n\t"\
			"vmovaps	%%ymm1,%%ymm2					/* cpy ~t1 */\n\t"\
			"vsubpd	     (%%rsi),%%ymm0,%%ymm0		/* r = B0 - t0 */\n\t"\
			"vmulpd	0x020(%%rax),%%ymm5,%%ymm5		/*~t3 = t3*sx0 */\n\t"\
			"vaddpd	      %%ymm3,%%ymm2,%%ymm2			/* ~t1+~t2 */\n\t"\
			"vmulpd	     (%%rcx),%%ymm3,%%ymm3		/* t2 = t2*cx2 */\n\t"\
			"vmulpd	0x020(%%rdx),%%ymm4,%%ymm4		/*~t4 = t4*sx3 */\n\t"\
			"vmulpd	     (%%rbx),%%ymm1,%%ymm1		/* t1 = t1*cx1 */\n\t"\
			"vmulpd	0x020(%%rbx),%%ymm6,%%ymm6		/*~t6 = t6*sx1 */\n\t"\
			"vmulpd	     (%%rax),%%ymm0,%%ymm0		/* ~r = r*(cx0-1) */\n\t"\
			"vmulpd	0x020(%%rcx),%%ymm7,%%ymm7		/*~t5 = t5*sx2 */\n\t"\
			"vmulpd	     (%%rdx),%%ymm2,%%ymm2		/* t3 =(t1+t2)*cx3 */\n\t"\
			"vaddpd	0x320(%%rdi),%%ymm0,%%ymm0		/* t0 =~r + B0 */\n\t"\
			"vaddpd	      %%ymm4,%%ymm6,%%ymm6			/*~t6 = t4+t6 */\n\t"\
			"vsubpd	      %%ymm2,%%ymm1,%%ymm1			/* tt = t1-t3 */\n\t"\
			"vsubpd	      %%ymm7,%%ymm4,%%ymm4			/*~t5 = t4-t5, ymm7 FREE */\n\t"\
			"vsubpd	      %%ymm2,%%ymm3,%%ymm3			/* t2 = t2-t3, ymm2 FREE */\n\t"\
			"vmovaps	%%ymm0,%%ymm2					/* cpy t0 */\n\t"\
			"vmovaps	%%ymm5,%%ymm7					/* cpy t3 */\n\t"\
			"vaddpd	      %%ymm1,%%ymm0,%%ymm0			/*~t0 = t0+tt */\n\t"\
			"vaddpd	      %%ymm6,%%ymm5,%%ymm5			/*~t3 = t3+t6 */\n\t"\
			"vaddpd	      %%ymm3,%%ymm1,%%ymm1			/*~tt = tt+t2 */\n\t"\
			"vaddpd	      %%ymm4,%%ymm6,%%ymm6			/*      t6+t5 */\n\t"\
			"vaddpd	      %%ymm2,%%ymm3,%%ymm3			/*~t2 = t2+t0 */\n\t"\
			"vaddpd	      %%ymm7,%%ymm4,%%ymm4			/*~t5 = t5+t3 */\n\t"\
			"vsubpd	      %%ymm1,%%ymm2,%%ymm2			/*~t1 = t0-tt-t2, ymm1 FREE */\n\t"\
			"vsubpd	      %%ymm6,%%ymm7,%%ymm7			/*~t4 = t3-t6-t5, ymm6 FREE */\n\t"\
			"\n\t"\
			"/* ymm1,6 FREE */\n\t"\
			"vmovaps	0x500(%%rdi),%%ymm1		/* t0r */\n\t"\
			"vmovaps	0x100(%%rdi),%%ymm6		/* t3r */\n\t"\
			"vsubpd	%%ymm5,%%ymm1,%%ymm1			/* B1r =t0r-t3i */\n\t"\
			"vsubpd	%%ymm6,%%ymm0,%%ymm0			/* B6i =t0i-t3r */\n\t"\
			"vaddpd	%%ymm5,%%ymm5,%%ymm5			/*        2*t3i */\n\t"\
			"vaddpd	%%ymm6,%%ymm6,%%ymm6			/*        2*t3r */\n\t"\
			"vaddpd	%%ymm1,%%ymm5,%%ymm5			/* B6r =t0r+t3i */\n\t"\
			"vaddpd	%%ymm0,%%ymm6,%%ymm6			/* B1i =t0i+t3r */\n\t"\
			"vmovaps	%%ymm1,0x500(%%rdi)		/* <-B1r = s1p08r */\n\t"\
			"vmovaps	%%ymm0,0x120(%%rdi)		/* <-B6i = s1p20r */\n\t"\
			"vmovaps	%%ymm5,0x100(%%rdi)		/* <-B6r = s1p20r */\n\t"\
			"vmovaps	%%ymm6,0x520(%%rdi)		/* <-B1i = s1p08i */\n\t"\
			"\n\t"\
			"vmovaps	      (%%rdi),%%ymm1		/* t1r */\n\t"\
			"vmovaps	0x600(%%rdi),%%ymm6		/* t4r */\n\t"\
			"vsubpd	%%ymm7,%%ymm1,%%ymm1			/* B2r =t1r-t4i */\n\t"\
			"vsubpd	%%ymm6,%%ymm2,%%ymm2			/* B5i =t1i-t4r */\n\t"\
			"vaddpd	%%ymm7,%%ymm7,%%ymm7			/*        2*t4i */\n\t"\
			"vaddpd	%%ymm6,%%ymm6,%%ymm6			/*        2*t4r */\n\t"\
			"vaddpd	%%ymm1,%%ymm7,%%ymm7			/* B5r =t1r+t4i */\n\t"\
			"vaddpd	%%ymm2,%%ymm6,%%ymm6			/* B2i =t1i+t4r */\n\t"\
			"vmovaps	%%ymm1,      (%%rdi)		/* <-B2r = s1p16r */\n\t"\
			"vmovaps	%%ymm2,0x620(%%rdi)		/* <-B5i = s1p12r */\n\t"\
			"vmovaps	%%ymm7,0x600(%%rdi)		/* <-B5r = s1p12r */\n\t"\
			"vmovaps	%%ymm6,0x020(%%rdi)		/* <-B2i = s1p16i */\n\t"\
			"\n\t"\
			"/* Note the order reversal on this pair of outputs: */\n\t"\
			"vmovaps	0x200(%%rdi),%%ymm0		/* t2r */\n\t"\
			"vmovaps	0x400(%%rdi),%%ymm5		/* t5r */\n\t"\
			"vsubpd	%%ymm4,%%ymm0,%%ymm0			/* B4r =t2r-t5i */\n\t"\
			"vsubpd	%%ymm5,%%ymm3,%%ymm3			/* B3i =t2i-t5r */\n\t"\
			"vaddpd	%%ymm4,%%ymm4,%%ymm4			/*        2*t5i */\n\t"\
			"vaddpd	%%ymm5,%%ymm5,%%ymm5			/*        2*t5r */\n\t"\
			"vaddpd	%%ymm0,%%ymm4,%%ymm4			/* B3r =t2r+t5i */\n\t"\
			"vaddpd	%%ymm3,%%ymm5,%%ymm5			/* B4i =t2i+t5r */\n\t"\
			"vmovaps	%%ymm0,0x400(%%rdi)		/* <-B4r = s1p04r */\n\t"\
			"vmovaps	%%ymm3,0x220(%%rdi)		/* <-B3i = s1p24r */\n\t"\
			"vmovaps	%%ymm4,0x200(%%rdi)		/* <-B3r = s1p24r */\n\t"\
			"vmovaps	%%ymm5,0x420(%%rdi)		/* <-B4i = s1p04i */\n\t"\
			"\n\t"\
		"/* SSE2_RADIX_07_DFT(s1p03r,s1p07r,s1p11r,s1p15r,s1p19r,s1p23r,s1p27r,cc0,s1p07r,s1p15r,s1p23r,s1p03r,s1p11r,s1p19r,s1p27r): */\n\t"\
			"addq	$0x040,%%rdi	/* s1p03r */\n\t"\
			"vmovaps	0x100(%%rdi),%%ymm6			/* A1r = s1p07r */\n\t"\
			"vmovaps	0x600(%%rdi),%%ymm1			/* A6r = s1p27r */\n\t"\
			"vmovaps	0x200(%%rdi),%%ymm5			/* A2r = s1p11r */\n\t"\
			"vmovaps	0x500(%%rdi),%%ymm2			/* A5r = s1p23r */\n\t"\
			"vmovaps	0x300(%%rdi),%%ymm4			/* A3r = s1p15r */\n\t"\
			"vmovaps	0x400(%%rdi),%%ymm3			/* A4r = s1p19r */\n\t"\
			"\n\t"\
			"vsubpd	%%ymm1,%%ymm6,%%ymm6			/* t6r = A1r-A6r */\n\t"\
			"vaddpd	%%ymm1,%%ymm1,%%ymm1			/*         2*A6r */\n\t"\
			"vaddpd	%%ymm6,%%ymm1,%%ymm1			/* t1r = A1r+A6r */\n\t"\
			"\n\t"\
			"vsubpd	%%ymm2,%%ymm5,%%ymm5			/* t5r = A2r-A5r */\n\t"\
			"vaddpd	%%ymm2,%%ymm2,%%ymm2			/*         2*A5r */\n\t"\
			"vaddpd	%%ymm5,%%ymm2,%%ymm2			/* t2r = A2r+A5r */\n\t"\
			"\n\t"\
			"vmovaps	     (%%rdi),%%ymm0		/* Ar0 */\n\t"\
			"vsubpd	%%ymm3,%%ymm4,%%ymm4			/* t4r = A3r-A4r */\n\t"\
			"vaddpd	%%ymm3,%%ymm3,%%ymm3			/*         2*A4r */\n\t"\
			"vaddpd	%%ymm4,%%ymm3,%%ymm3			/* t3r = A3r+A4r */\n\t"\
			"\n\t"\
			"vmovaps	%%ymm0,     (%%rsi)		/* cpy t0 into scratch sincos slot */\n\t"\
			"vmovaps	%%ymm6,0x020(%%rsi)		/* cpy t6 into scratch sincos slot */\n\t"\
			"vaddpd	      %%ymm1,%%ymm0,%%ymm0		/*~A0 = A0+t1 */\n\t"\
			"vmovaps	%%ymm5,%%ymm7			/* cpy t5 */\n\t"\
			"vaddpd	      %%ymm2,%%ymm3,%%ymm3		/*~t3 = t3+t2 */\n\t"\
			"vsubpd	      %%ymm4,%%ymm5,%%ymm5		/*~t5 = t5-t4 */\n\t"\
			"vsubpd	      %%ymm2,%%ymm1,%%ymm1		/*~t1 = t1-t2 */\n\t"\
			"vsubpd	      %%ymm7,%%ymm6,%%ymm6		/*~t6 = t6-t5 */\n\t"\
			"vaddpd	      %%ymm2,%%ymm2,%%ymm2		/* 2*t2 */\n\t"\
			"vaddpd	      %%ymm7,%%ymm4,%%ymm4		/*~t5 = t4+t5 */\n\t"\
			"vaddpd	      %%ymm3,%%ymm0,%%ymm0		/* B0 */\n\t"\
			"vaddpd	0x020(%%rsi),%%ymm5,%%ymm5		/* t3 = (t5-t4)+t6 */\n\t"\
			"vsubpd	      %%ymm2,%%ymm3,%%ymm3		/*~t2 =  (t2+t3) - 2*t2 = t3-t2 */\n\t"\
			"vmovaps	%%ymm4,%%ymm7			/* cpy t5 */\n\t"\
			"vmovaps	%%ymm0,0x100(%%rdi)		/* <-B0 = s1p00r, ymm0 FREE */\n\t"\
			"vsubpd	      %%ymm6,%%ymm4,%%ymm4		/* t4 = ~t5-~t6 */\n\t"\
			"vmovaps	%%ymm1,%%ymm2			/* cpy ~t1 */\n\t"\
			"vsubpd	     (%%rsi),%%ymm0,%%ymm0		/* r = B0 - t0 */\n\t"\
			"vmulpd	0x020(%%rax),%%ymm5,%%ymm5		/*~t3 = t3*sx0 */\n\t"\
			"vaddpd	      %%ymm3,%%ymm2,%%ymm2		/* ~t1+~t2 */\n\t"\
			"vmulpd	     (%%rcx),%%ymm3,%%ymm3		/* t2 = t2*cx2 */\n\t"\
			"vmulpd	0x020(%%rdx),%%ymm4,%%ymm4		/*~t4 = t4*sx3 */\n\t"\
			"vmulpd	     (%%rbx),%%ymm1,%%ymm1		/* t1 = t1*cx1 */\n\t"\
			"vmulpd	0x020(%%rbx),%%ymm6,%%ymm6		/*~t6 = t6*sx1 */\n\t"\
			"vmulpd	     (%%rax),%%ymm0,%%ymm0		/* ~r = r*(cx0-1) */\n\t"\
			"vmulpd	0x020(%%rcx),%%ymm7,%%ymm7		/*~t5 = t5*sx2 */\n\t"\
			"vmulpd	     (%%rdx),%%ymm2,%%ymm2		/* t3 =(t1+t2)*cx3 */\n\t"\
			"vaddpd	0x100(%%rdi),%%ymm0,%%ymm0		/* t0 =~r + B0 */\n\t"\
			"vaddpd	      %%ymm4,%%ymm6,%%ymm6		/*~t6 = t4+t6 */\n\t"\
			"vsubpd	      %%ymm2,%%ymm1,%%ymm1		/* tt = t1-t3 */\n\t"\
			"vsubpd	      %%ymm7,%%ymm4,%%ymm4		/*~t5 = t4-t5, ymm7 FREE */\n\t"\
			"vsubpd	      %%ymm2,%%ymm3,%%ymm3		/* t2 = t2-t3, ymm2 FREE */\n\t"\
			"\n\t"\
			"vmovaps	%%ymm0,%%ymm2			/* cpy t0 */\n\t"\
			"vmovaps	%%ymm5,%%ymm7			/* cpy t3 */\n\t"\
			"vaddpd	      %%ymm1,%%ymm0,%%ymm0			/*~t0 = t0+tt */\n\t"\
			"vaddpd	      %%ymm6,%%ymm5,%%ymm5			/*~t3 = t3+t6 */\n\t"\
			"vaddpd	      %%ymm3,%%ymm1,%%ymm1			/*~tt = tt+t2 */\n\t"\
			"vaddpd	      %%ymm4,%%ymm6,%%ymm6			/*      t6+t5 */\n\t"\
			"vaddpd	      %%ymm2,%%ymm3,%%ymm3			/*~t2 = t2+t0 */\n\t"\
			"vaddpd	      %%ymm7,%%ymm4,%%ymm4			/*~t5 = t5+t3 */\n\t"\
			"vsubpd	      %%ymm1,%%ymm2,%%ymm2			/*~t1 = t0-tt-t2 */\n\t"\
			"vsubpd	      %%ymm6,%%ymm7,%%ymm7			/*~t4 = t3-t6-t5 */\n\t"\
			"vmovaps	%%ymm0,0x300(%%rdi)			/* B1 <- t0 = s1p15r */\n\t"\
			"vmovaps	%%ymm5,0x600(%%rdi)			/* B6 <- t3 = s1p27r */\n\t"\
			"vmovaps	%%ymm2,0x500(%%rdi)			/* B2 <- t1 = s1p23r */\n\t"\
			"vmovaps	%%ymm7,0x400(%%rdi)			/* B5 <- t4 = s1p19r */\n\t"\
			"vmovaps	%%ymm3,      (%%rdi)			/* B3 <- t2 = s1p03r */\n\t"\
			"vmovaps	%%ymm4,0x200(%%rdi)			/* B4 <- t5 = s1p11r */\n\t"\
			"/************************** Imaginary Parts: ******************************************/\n\t"\
			"vmovaps	0x120(%%rdi),%%ymm6		/* A1i = s1p07i */\n\t"\
			"vmovaps	0x620(%%rdi),%%ymm1		/* A6i = s1p27i */\n\t"\
			"vmovaps	0x220(%%rdi),%%ymm5		/* A2i = s1p11i */\n\t"\
			"vmovaps	0x520(%%rdi),%%ymm2		/* A5i = s1p23i */\n\t"\
			"vmovaps	0x320(%%rdi),%%ymm4		/* A3i = s1p15i */\n\t"\
			"vmovaps	0x420(%%rdi),%%ymm3		/* A4i = s1p19i */\n\t"\
			"\n\t"\
			"vsubpd	%%ymm1,%%ymm6,%%ymm6			/* t6i = A1i-A6i */\n\t"\
			"vaddpd	%%ymm1,%%ymm1,%%ymm1			/*         2*A6i */\n\t"\
			"vaddpd	%%ymm6,%%ymm1,%%ymm1			/* t1i = A1i+A6i */\n\t"\
			"\n\t"\
			"vsubpd	%%ymm2,%%ymm5,%%ymm5			/* t5i = A2i-A5i */\n\t"\
			"vaddpd	%%ymm2,%%ymm2,%%ymm2			/*         2*A5i */\n\t"\
			"vaddpd	%%ymm5,%%ymm2,%%ymm2			/* t2i = A2i+A5i */\n\t"\
			"\n\t"\
			"vmovaps	0x020(%%rdi),%%ymm0		/* Ai0 */\n\t"\
			"vsubpd	%%ymm3,%%ymm4,%%ymm4			/* t4i = A3i-A4i */\n\t"\
			"vaddpd	%%ymm3,%%ymm3,%%ymm3			/*         2*A4i */\n\t"\
			"vaddpd	%%ymm4,%%ymm3,%%ymm3			/* t3i = A3i+A4i */\n\t"\
			"\n\t"\
			"vmovaps	%%ymm0,     (%%rsi)		/* cpy t0 into scratch sincos slot */\n\t"\
			"vmovaps	%%ymm6,0x020(%%rsi)		/* cpy t6 into scratch sincos slot */\n\t"\
			"vaddpd	      %%ymm1,%%ymm0,%%ymm0			/*~A0 = A0+t1 */\n\t"\
			"vmovaps	%%ymm5,%%ymm7					/* cpy t5 */\n\t"\
			"vaddpd	      %%ymm2,%%ymm3,%%ymm3			/*~t3 = t3+t2 */\n\t"\
			"vsubpd	      %%ymm4,%%ymm5,%%ymm5			/*~t5 = t5-t4 */\n\t"\
			"vsubpd	      %%ymm2,%%ymm1,%%ymm1			/*~t1 = t1-t2 */\n\t"\
			"vsubpd	      %%ymm7,%%ymm6,%%ymm6			/*~t6 = t6-t5 */\n\t"\
			"vaddpd	      %%ymm2,%%ymm2,%%ymm2			/* 2*t2 */\n\t"\
			"vaddpd	      %%ymm7,%%ymm4,%%ymm4			/*~t5 = t4+t5 */\n\t"\
			"vaddpd	      %%ymm3,%%ymm0,%%ymm0			/* B0 */\n\t"\
			"vaddpd	0x020(%%rsi),%%ymm5,%%ymm5		/* t3 = (t5-t4)+t6 */\n\t"\
			"vsubpd	      %%ymm2,%%ymm3,%%ymm3			/*~t2 =  (t2+t3) - 2*t2 = t3-t2 */\n\t"\
			"vmovaps	%%ymm4,%%ymm7					/* cpy t5 */\n\t"\
			"vmovaps	%%ymm0,0x120(%%rdi)		/* <-B0 = s1p07i, ymm0 FREE */\n\t"\
			"vsubpd	      %%ymm6,%%ymm4,%%ymm4			/* t4 = ~t5-~t6 */\n\t"\
			"vmovaps	%%ymm1,%%ymm2					/* cpy ~t1 */\n\t"\
			"vsubpd	     (%%rsi),%%ymm0,%%ymm0		/* r = B0 - t0 */\n\t"\
			"vmulpd	0x020(%%rax),%%ymm5,%%ymm5		/*~t3 = t3*sx0 */\n\t"\
			"vaddpd	      %%ymm3,%%ymm2,%%ymm2			/* ~t1+~t2 */\n\t"\
			"vmulpd	     (%%rcx),%%ymm3,%%ymm3		/* t2 = t2*cx2 */\n\t"\
			"vmulpd	0x020(%%rdx),%%ymm4,%%ymm4		/*~t4 = t4*sx3 */\n\t"\
			"vmulpd	     (%%rbx),%%ymm1,%%ymm1		/* t1 = t1*cx1 */\n\t"\
			"vmulpd	0x020(%%rbx),%%ymm6,%%ymm6		/*~t6 = t6*sx1 */\n\t"\
			"vmulpd	     (%%rax),%%ymm0,%%ymm0		/* ~r = r*(cx0-1) */\n\t"\
			"vmulpd	0x020(%%rcx),%%ymm7,%%ymm7		/*~t5 = t5*sx2 */\n\t"\
			"vmulpd	     (%%rdx),%%ymm2,%%ymm2		/* t3 =(t1+t2)*cx3 */\n\t"\
			"vaddpd	0x120(%%rdi),%%ymm0,%%ymm0		/* t0 =~r + B0 */\n\t"\
			"vaddpd	      %%ymm4,%%ymm6,%%ymm6			/*~t6 = t4+t6 */\n\t"\
			"vsubpd	      %%ymm2,%%ymm1,%%ymm1			/* tt = t1-t3 */\n\t"\
			"vsubpd	      %%ymm7,%%ymm4,%%ymm4			/*~t5 = t4-t5, ymm7 FREE */\n\t"\
			"vsubpd	      %%ymm2,%%ymm3,%%ymm3			/* t2 = t2-t3, ymm2 FREE */\n\t"\
			"vmovaps	%%ymm0,%%ymm2					/* cpy t0 */\n\t"\
			"vmovaps	%%ymm5,%%ymm7					/* cpy t3 */\n\t"\
			"vaddpd	      %%ymm1,%%ymm0,%%ymm0			/*~t0 = t0+tt */\n\t"\
			"vaddpd	      %%ymm6,%%ymm5,%%ymm5			/*~t3 = t3+t6 */\n\t"\
			"vaddpd	      %%ymm3,%%ymm1,%%ymm1			/*~tt = tt+t2 */\n\t"\
			"vaddpd	      %%ymm4,%%ymm6,%%ymm6			/*      t6+t5 */\n\t"\
			"vaddpd	      %%ymm2,%%ymm3,%%ymm3			/*~t2 = t2+t0 */\n\t"\
			"vaddpd	      %%ymm7,%%ymm4,%%ymm4			/*~t5 = t5+t3 */\n\t"\
			"vsubpd	      %%ymm1,%%ymm2,%%ymm2			/*~t1 = t0-tt-t2, ymm1 FREE */\n\t"\
			"vsubpd	      %%ymm6,%%ymm7,%%ymm7			/*~t4 = t3-t6-t5, ymm6 FREE */\n\t"\
			"\n\t"\
			"/* ymm1,6 FREE */\n\t"\
			"vmovaps	0x300(%%rdi),%%ymm1		/* t0r */\n\t"\
			"vmovaps	0x600(%%rdi),%%ymm6		/* t3r */\n\t"\
			"vsubpd	%%ymm5,%%ymm1,%%ymm1			/* B1r =t0r-t3i */\n\t"\
			"vsubpd	%%ymm6,%%ymm0,%%ymm0			/* B6i =t0i-t3r */\n\t"\
			"vaddpd	%%ymm5,%%ymm5,%%ymm5			/*        2*t3i */\n\t"\
			"vaddpd	%%ymm6,%%ymm6,%%ymm6			/*        2*t3r */\n\t"\
			"vaddpd	%%ymm1,%%ymm5,%%ymm5			/* B6r =t0r+t3i */\n\t"\
			"vaddpd	%%ymm0,%%ymm6,%%ymm6			/* B1i =t0i+t3r */\n\t"\
			"vmovaps	%%ymm1,0x300(%%rdi)		/* <-B1r = s1p08r */\n\t"\
			"vmovaps	%%ymm0,0x620(%%rdi)		/* <-B6i = s1p20r */\n\t"\
			"vmovaps	%%ymm5,0x600(%%rdi)		/* <-B6r = s1p20r */\n\t"\
			"vmovaps	%%ymm6,0x320(%%rdi)		/* <-B1i = s1p08i */\n\t"\
			"\n\t"\
			"vmovaps	0x500(%%rdi),%%ymm1		/* t1r */\n\t"\
			"vmovaps	0x400(%%rdi),%%ymm6		/* t4r */\n\t"\
			"vsubpd	%%ymm7,%%ymm1,%%ymm1			/* B2r =t1r-t4i */\n\t"\
			"vsubpd	%%ymm6,%%ymm2,%%ymm2			/* B5i =t1i-t4r */\n\t"\
			"vaddpd	%%ymm7,%%ymm7,%%ymm7			/*        2*t4i */\n\t"\
			"vaddpd	%%ymm6,%%ymm6,%%ymm6			/*        2*t4r */\n\t"\
			"vaddpd	%%ymm1,%%ymm7,%%ymm7			/* B5r =t1r+t4i */\n\t"\
			"vaddpd	%%ymm2,%%ymm6,%%ymm6			/* B2i =t1i+t4r */\n\t"\
			"vmovaps	%%ymm1,0x500(%%rdi)		/* <-B2r = s1p16r */\n\t"\
			"vmovaps	%%ymm2,0x420(%%rdi)		/* <-B5i = s1p12r */\n\t"\
			"vmovaps	%%ymm7,0x400(%%rdi)		/* <-B5r = s1p12r */\n\t"\
			"vmovaps	%%ymm6,0x520(%%rdi)		/* <-B2i = s1p16i */\n\t"\
			"\n\t"\
			"/* Note the order reversal on this pair of outputs: */\n\t"\
			"vmovaps	      (%%rdi),%%ymm0		/* t2r */\n\t"\
			"vmovaps	0x200(%%rdi),%%ymm5		/* t5r */\n\t"\
			"vsubpd	%%ymm4,%%ymm0,%%ymm0			/* B4r =t2r-t5i */\n\t"\
			"vsubpd	%%ymm5,%%ymm3,%%ymm3			/* B3i =t2i-t5r */\n\t"\
			"vaddpd	%%ymm4,%%ymm4,%%ymm4			/*        2*t5i */\n\t"\
			"vaddpd	%%ymm5,%%ymm5,%%ymm5			/*        2*t5r */\n\t"\
			"vaddpd	%%ymm0,%%ymm4,%%ymm4			/* B3r =t2r+t5i */\n\t"\
			"vaddpd	%%ymm3,%%ymm5,%%ymm5			/* B4i =t2i+t5r */\n\t"\
			"vmovaps	%%ymm0,0x200(%%rdi)		/* <-B4r = s1p04r */\n\t"\
			"vmovaps	%%ymm3,0x020(%%rdi)		/* <-B3i = s1p24r */\n\t"\
			"vmovaps	%%ymm4,      (%%rdi)		/* <-B3r = s1p24r */\n\t"\
			"vmovaps	%%ymm5,0x220(%%rdi)		/* <-B4i = s1p04i */\n\t"\
			"\n\t"\
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
			: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
		);\
	}


	#define	SSE2_RADIX28_DIF_NOTWIDDLE(Xadd,Xp01,Xp02,Xp03,Xp04,Xp08,Xp12,Xp16,Xp20,Xp24,Xout,Xcc0)\
	{\
	__asm__ volatile (\
		"/* SSE2_RADIX_07_DFT(s1p00r,s1p24r,s1p20r,s1p16r,s1p12r,s1p08r,s1p04r,cc0,s1p00r,s1p04r,s1p08r,s1p12r,s1p16r,s1p20r,s1p24r): */\n\t"\
			"movq	%[__out],%%rdi	/* s1p00r */\n\t"\
			"movq	%[__cc0],%%rax	/* cc0 */\n\t"\
			"movq	$0x040,%%rbx\n\t"\
			"movq	$0x080,%%rcx\n\t"\
			"movq	$0x0c0,%%rdx\n\t"\
			"movq	$0x100,%%rsi\n\t"\
			"addq	%%rax,%%rbx		/* cc1 */\n\t"\
			"addq	%%rax,%%rcx		/* cc2 */\n\t"\
			"addq	%%rax,%%rdx		/* cc3 */\n\t"\
			"addq	%%rax,%%rsi		/* scratch storage */\n\t"\
			"vmovaps	0x600(%%rdi),%%ymm6			/* A1r = s1p24r */\n\t"\
			"vmovaps	0x100(%%rdi),%%ymm1			/* A6r = s1p04r */\n\t"\
			"vmovaps	0x500(%%rdi),%%ymm5			/* A2r = s1p20r */\n\t"\
			"vmovaps	0x200(%%rdi),%%ymm2			/* A5r = s1p08r */\n\t"\
			"vmovaps	0x400(%%rdi),%%ymm4			/* A3r = s1p16r */\n\t"\
			"vmovaps	0x300(%%rdi),%%ymm3			/* A4r = s1p12r */\n\t"\
			"\n\t"\
			"vsubpd	%%ymm1,%%ymm6,%%ymm6			/* t6r = A1r-A6r */\n\t"\
			"vaddpd	%%ymm1,%%ymm1,%%ymm1			/*         2*A6r */\n\t"\
			"vaddpd	%%ymm6,%%ymm1,%%ymm1			/* t1r = A1r+A6r */\n\t"\
			"\n\t"\
			"vsubpd	%%ymm2,%%ymm5,%%ymm5			/* t5r = A2r-A5r */\n\t"\
			"vaddpd	%%ymm2,%%ymm2,%%ymm2			/*         2*A5r */\n\t"\
			"vaddpd	%%ymm5,%%ymm2,%%ymm2			/* t2r = A2r+A5r */\n\t"\
			"\n\t"\
			"vmovaps	     (%%rdi),%%ymm0		/* Ar0 */\n\t"\
			"vsubpd	%%ymm3,%%ymm4,%%ymm4			/* t4r = A3r-A4r */\n\t"\
			"vaddpd	%%ymm3,%%ymm3,%%ymm3			/*         2*A4r */\n\t"\
			"vaddpd	%%ymm4,%%ymm3,%%ymm3			/* t3r = A3r+A4r */\n\t"\
			"\n\t"\
			"vmovaps	%%ymm0,     (%%rsi)		/* cpy t0 into scratch sincos slot */\n\t"\
			"vmovaps	%%ymm6,0x020(%%rsi)		/* cpy t6 into scratch sincos slot */\n\t"\
			"vaddpd	      %%ymm1,%%ymm0,%%ymm0		/*~A0 = A0+t1 */\n\t"\
			"vmovaps	%%ymm5,%%ymm7			/* cpy t5 */\n\t"\
			"vaddpd	      %%ymm2,%%ymm3,%%ymm3		/*~t3 = t3+t2 */\n\t"\
			"vsubpd	      %%ymm4,%%ymm5,%%ymm5		/*~t5 = t5-t4 */\n\t"\
			"vsubpd	      %%ymm2,%%ymm1,%%ymm1		/*~t1 = t1-t2 */\n\t"\
			"vsubpd	      %%ymm7,%%ymm6,%%ymm6		/*~t6 = t6-t5 */\n\t"\
			"vaddpd	      %%ymm2,%%ymm2,%%ymm2		/* 2*t2 */\n\t"\
			"vaddpd	      %%ymm7,%%ymm4,%%ymm4		/*~t5 = t4+t5 */\n\t"\
			"vaddpd	      %%ymm3,%%ymm0,%%ymm0		/* B0 */\n\t"\
			"vaddpd	0x020(%%rsi),%%ymm5,%%ymm5		/* t3 = (t5-t4)+t6 */\n\t"\
			"vsubpd	      %%ymm2,%%ymm3,%%ymm3		/*~t2 =  (t2+t3) - 2*t2 = t3-t2 */\n\t"\
			"vmovaps	%%ymm4,%%ymm7			/* cpy t5 */\n\t"\
			"vmovaps	%%ymm0,      (%%rdi)		/* <-B0 = s1p00r, ymm0 FREE */\n\t"\
			"vsubpd	      %%ymm6,%%ymm4,%%ymm4		/* t4 = ~t5-~t6 */\n\t"\
			"vmovaps	%%ymm1,%%ymm2			/* cpy ~t1 */\n\t"\
			"vsubpd	     (%%rsi),%%ymm0,%%ymm0		/* r = B0 - t0 */\n\t"\
			"vmulpd	0x020(%%rax),%%ymm5,%%ymm5		/*~t3 = t3*sx0 */\n\t"\
			"vaddpd	      %%ymm3,%%ymm2,%%ymm2		/* ~t1+~t2 */\n\t"\
			"vmulpd	     (%%rcx),%%ymm3,%%ymm3		/* t2 = t2*cx2 */\n\t"\
			"vmulpd	0x020(%%rdx),%%ymm4,%%ymm4		/*~t4 = t4*sx3 */\n\t"\
			"vmulpd	     (%%rbx),%%ymm1,%%ymm1		/* t1 = t1*cx1 */\n\t"\
			"vmulpd	0x020(%%rbx),%%ymm6,%%ymm6		/*~t6 = t6*sx1 */\n\t"\
			"vmulpd	     (%%rax),%%ymm0,%%ymm0		/* ~r = r*(cx0-1) */\n\t"\
			"vmulpd	0x020(%%rcx),%%ymm7,%%ymm7		/*~t5 = t5*sx2 */\n\t"\
			"vmulpd	     (%%rdx),%%ymm2,%%ymm2		/* t3 =(t1+t2)*cx3 */\n\t"\
			"vaddpd	     (%%rdi),%%ymm0,%%ymm0		/* t0 =~r + B0 */\n\t"\
			"vaddpd	      %%ymm4,%%ymm6,%%ymm6		/*~t6 = t4+t6 */\n\t"\
			"vsubpd	      %%ymm2,%%ymm1,%%ymm1		/* tt = t1-t3 */\n\t"\
			"vsubpd	      %%ymm7,%%ymm4,%%ymm4		/*~t5 = t4-t5, ymm7 FREE */\n\t"\
			"vsubpd	      %%ymm2,%%ymm3,%%ymm3		/* t2 = t2-t3, ymm2 FREE */\n\t"\
			"\n\t"\
			"vmovaps	%%ymm0,%%ymm2			/* cpy t0 */\n\t"\
			"vmovaps	%%ymm5,%%ymm7			/* cpy t3 */\n\t"\
			"vaddpd	      %%ymm1,%%ymm0,%%ymm0			/*~t0 = t0+tt */\n\t"\
			"vaddpd	      %%ymm6,%%ymm5,%%ymm5			/*~t3 = t3+t6 */\n\t"\
			"vaddpd	      %%ymm3,%%ymm1,%%ymm1			/*~tt = tt+t2 */\n\t"\
			"vaddpd	      %%ymm4,%%ymm6,%%ymm6			/*      t6+t5 */\n\t"\
			"vaddpd	      %%ymm2,%%ymm3,%%ymm3			/*~t2 = t2+t0 */\n\t"\
			"vaddpd	      %%ymm7,%%ymm4,%%ymm4			/*~t5 = t5+t3 */\n\t"\
			"vsubpd	      %%ymm1,%%ymm2,%%ymm2			/*~t1 = t0-tt-t2 */\n\t"\
			"vsubpd	      %%ymm6,%%ymm7,%%ymm7			/*~t4 = t3-t6-t5 */\n\t"\
			"vmovaps	%%ymm0,0x100(%%rdi)			/* B1 <- t0 = s1p04r */\n\t"\
			"vmovaps	%%ymm5,0x600(%%rdi)			/* B6 <- t3 = s1p24r */\n\t"\
			"vmovaps	%%ymm2,0x200(%%rdi)			/* B2 <- t1 = s1p08r */\n\t"\
			"vmovaps	%%ymm7,0x500(%%rdi)			/* B5 <- t4 = s1p20r */\n\t"\
			"vmovaps	%%ymm3,0x300(%%rdi)			/* B3 <- t2 = s1p12r */\n\t"\
			"vmovaps	%%ymm4,0x400(%%rdi)			/* B4 <- t5 = s1p16r */\n\t"\
			"/************************** Imaginary Parts: ******************************************/\n\t"\
			"vmovaps	0x620(%%rdi),%%ymm6		/* A1i = s1p24i */\n\t"\
			"vmovaps	0x120(%%rdi),%%ymm1		/* A6i = s1p04i */\n\t"\
			"vmovaps	0x520(%%rdi),%%ymm5		/* A2i = s1p20i */\n\t"\
			"vmovaps	0x220(%%rdi),%%ymm2		/* A5i = s1p08i */\n\t"\
			"vmovaps	0x420(%%rdi),%%ymm4		/* A3i = s1p16i */\n\t"\
			"vmovaps	0x320(%%rdi),%%ymm3		/* A4i = s1p12i */\n\t"\
			"\n\t"\
			"vsubpd	%%ymm1,%%ymm6,%%ymm6			/* t6i = A1i-A6i */\n\t"\
			"vaddpd	%%ymm1,%%ymm1,%%ymm1			/*         2*A6i */\n\t"\
			"vaddpd	%%ymm6,%%ymm1,%%ymm1			/* t1i = A1i+A6i */\n\t"\
			"\n\t"\
			"vsubpd	%%ymm2,%%ymm5,%%ymm5			/* t5i = A2i-A5i */\n\t"\
			"vaddpd	%%ymm2,%%ymm2,%%ymm2			/*         2*A5i */\n\t"\
			"vaddpd	%%ymm5,%%ymm2,%%ymm2			/* t2i = A2i+A5i */\n\t"\
			"\n\t"\
			"vmovaps	0x020(%%rdi),%%ymm0		/* Ai0 */\n\t"\
			"vsubpd	%%ymm3,%%ymm4,%%ymm4			/* t4i = A3i-A4i */\n\t"\
			"vaddpd	%%ymm3,%%ymm3,%%ymm3			/*         2*A4i */\n\t"\
			"vaddpd	%%ymm4,%%ymm3,%%ymm3			/* t3i = A3i+A4i */\n\t"\
			"\n\t"\
			"vmovaps	%%ymm0,     (%%rsi)		/* cpy t0 into scratch sincos slot */\n\t"\
			"vmovaps	%%ymm6,0x020(%%rsi)		/* cpy t6 into scratch sincos slot */\n\t"\
			"vaddpd	      %%ymm1,%%ymm0,%%ymm0			/*~A0 = A0+t1 */\n\t"\
			"vmovaps	%%ymm5,%%ymm7					/* cpy t5 */\n\t"\
			"vaddpd	      %%ymm2,%%ymm3,%%ymm3			/*~t3 = t3+t2 */\n\t"\
			"vsubpd	      %%ymm4,%%ymm5,%%ymm5			/*~t5 = t5-t4 */\n\t"\
			"vsubpd	      %%ymm2,%%ymm1,%%ymm1			/*~t1 = t1-t2 */\n\t"\
			"vsubpd	      %%ymm7,%%ymm6,%%ymm6			/*~t6 = t6-t5 */\n\t"\
			"vaddpd	      %%ymm2,%%ymm2,%%ymm2			/* 2*t2 */\n\t"\
			"vaddpd	      %%ymm7,%%ymm4,%%ymm4			/*~t5 = t4+t5 */\n\t"\
			"vaddpd	      %%ymm3,%%ymm0,%%ymm0			/* B0 */\n\t"\
			"vaddpd	0x020(%%rsi),%%ymm5,%%ymm5		/* t3 = (t5-t4)+t6 */\n\t"\
			"vsubpd	      %%ymm2,%%ymm3,%%ymm3			/*~t2 =  (t2+t3) - 2*t2 = t3-t2 */\n\t"\
			"vmovaps	%%ymm4,%%ymm7					/* cpy t5 */\n\t"\
			"vmovaps	%%ymm0,0x020(%%rdi)		/* <-B0 = s1p00i, ymm0 FREE */\n\t"\
			"vsubpd	      %%ymm6,%%ymm4,%%ymm4			/* t4 = ~t5-~t6 */\n\t"\
			"vmovaps	%%ymm1,%%ymm2					/* cpy ~t1 */\n\t"\
			"vsubpd	     (%%rsi),%%ymm0,%%ymm0		/* r = B0 - t0 */\n\t"\
			"vmulpd	0x020(%%rax),%%ymm5,%%ymm5		/*~t3 = t3*sx0 */\n\t"\
			"vaddpd	      %%ymm3,%%ymm2,%%ymm2			/* ~t1+~t2 */\n\t"\
			"vmulpd	     (%%rcx),%%ymm3,%%ymm3		/* t2 = t2*cx2 */\n\t"\
			"vmulpd	0x020(%%rdx),%%ymm4,%%ymm4		/*~t4 = t4*sx3 */\n\t"\
			"vmulpd	     (%%rbx),%%ymm1,%%ymm1		/* t1 = t1*cx1 */\n\t"\
			"vmulpd	0x020(%%rbx),%%ymm6,%%ymm6		/*~t6 = t6*sx1 */\n\t"\
			"vmulpd	     (%%rax),%%ymm0,%%ymm0		/* ~r = r*(cx0-1) */\n\t"\
			"vmulpd	0x020(%%rcx),%%ymm7,%%ymm7		/*~t5 = t5*sx2 */\n\t"\
			"vmulpd	     (%%rdx),%%ymm2,%%ymm2		/* t3 =(t1+t2)*cx3 */\n\t"\
			"vaddpd	0x020(%%rdi),%%ymm0,%%ymm0		/* t0 =~r + B0 */\n\t"\
			"vaddpd	      %%ymm4,%%ymm6,%%ymm6			/*~t6 = t4+t6 */\n\t"\
			"vsubpd	      %%ymm2,%%ymm1,%%ymm1			/* tt = t1-t3 */\n\t"\
			"vsubpd	      %%ymm7,%%ymm4,%%ymm4			/*~t5 = t4-t5, ymm7 FREE */\n\t"\
			"vsubpd	      %%ymm2,%%ymm3,%%ymm3			/* t2 = t2-t3, ymm2 FREE */\n\t"\
			"vmovaps	%%ymm0,%%ymm2					/* cpy t0 */\n\t"\
			"vmovaps	%%ymm5,%%ymm7					/* cpy t3 */\n\t"\
			"vaddpd	      %%ymm1,%%ymm0,%%ymm0			/*~t0 = t0+tt */\n\t"\
			"vaddpd	      %%ymm6,%%ymm5,%%ymm5			/*~t3 = t3+t6 */\n\t"\
			"vaddpd	      %%ymm3,%%ymm1,%%ymm1			/*~tt = tt+t2 */\n\t"\
			"vaddpd	      %%ymm4,%%ymm6,%%ymm6			/*      t6+t5 */\n\t"\
			"vaddpd	      %%ymm2,%%ymm3,%%ymm3			/*~t2 = t2+t0 */\n\t"\
			"vaddpd	      %%ymm7,%%ymm4,%%ymm4			/*~t5 = t5+t3 */\n\t"\
			"vsubpd	      %%ymm1,%%ymm2,%%ymm2			/*~t1 = t0-tt-t2, ymm1 FREE */\n\t"\
			"vsubpd	      %%ymm6,%%ymm7,%%ymm7			/*~t4 = t3-t6-t5, ymm6 FREE */\n\t"\
			"\n\t"\
			"/* ymm1,6 FREE */\n\t"\
			"vmovaps	0x100(%%rdi),%%ymm1		/* t0r */\n\t"\
			"vmovaps	0x600(%%rdi),%%ymm6		/* t3r */\n\t"\
			"vsubpd	%%ymm5,%%ymm1,%%ymm1			/* B1r =t0r-t3i */\n\t"\
			"vsubpd	%%ymm6,%%ymm0,%%ymm0			/* B6i =t0i-t3r */\n\t"\
			"vaddpd	%%ymm5,%%ymm5,%%ymm5			/*        2*t3i */\n\t"\
			"vaddpd	%%ymm6,%%ymm6,%%ymm6			/*        2*t3r */\n\t"\
			"vaddpd	%%ymm1,%%ymm5,%%ymm5			/* B6r =t0r+t3i */\n\t"\
			"vaddpd	%%ymm0,%%ymm6,%%ymm6			/* B1i =t0i+t3r */\n\t"\
			"vmovaps	%%ymm1,0x100(%%rdi)		/* <-B1r = s1p08r */\n\t"\
			"vmovaps	%%ymm0,0x620(%%rdi)		/* <-B6i = s1p20r */\n\t"\
			"vmovaps	%%ymm5,0x600(%%rdi)		/* <-B6r = s1p20r */\n\t"\
			"vmovaps	%%ymm6,0x120(%%rdi)		/* <-B1i = s1p08i */\n\t"\
			"\n\t"\
			"vmovaps	0x200(%%rdi),%%ymm1		/* t1r */\n\t"\
			"vmovaps	0x500(%%rdi),%%ymm6		/* t4r */\n\t"\
			"vsubpd	%%ymm7,%%ymm1,%%ymm1			/* B2r =t1r-t4i */\n\t"\
			"vsubpd	%%ymm6,%%ymm2,%%ymm2			/* B5i =t1i-t4r */\n\t"\
			"vaddpd	%%ymm7,%%ymm7,%%ymm7			/*        2*t4i */\n\t"\
			"vaddpd	%%ymm6,%%ymm6,%%ymm6			/*        2*t4r */\n\t"\
			"vaddpd	%%ymm1,%%ymm7,%%ymm7			/* B5r =t1r+t4i */\n\t"\
			"vaddpd	%%ymm2,%%ymm6,%%ymm6			/* B2i =t1i+t4r */\n\t"\
			"vmovaps	%%ymm1,0x200(%%rdi)		/* <-B2r = s1p16r */\n\t"\
			"vmovaps	%%ymm2,0x520(%%rdi)		/* <-B5i = s1p12r */\n\t"\
			"vmovaps	%%ymm7,0x500(%%rdi)		/* <-B5r = s1p12r */\n\t"\
			"vmovaps	%%ymm6,0x220(%%rdi)		/* <-B2i = s1p16i */\n\t"\
			"\n\t"\
			"/* Note the order reversal on this pair of outputs: */\n\t"\
			"vmovaps	0x300(%%rdi),%%ymm0		/* t2r */\n\t"\
			"vmovaps	0x400(%%rdi),%%ymm5		/* t5r */\n\t"\
			"vsubpd	%%ymm4,%%ymm0,%%ymm0			/* B4r =t2r-t5i */\n\t"\
			"vsubpd	%%ymm5,%%ymm3,%%ymm3			/* B3i =t2i-t5r */\n\t"\
			"vaddpd	%%ymm4,%%ymm4,%%ymm4			/*        2*t5i */\n\t"\
			"vaddpd	%%ymm5,%%ymm5,%%ymm5			/*        2*t5r */\n\t"\
			"vaddpd	%%ymm0,%%ymm4,%%ymm4			/* B3r =t2r+t5i */\n\t"\
			"vaddpd	%%ymm3,%%ymm5,%%ymm5			/* B4i =t2i+t5r */\n\t"\
			"vmovaps	%%ymm0,0x400(%%rdi)		/* <-B4r = s1p04r */\n\t"\
			"vmovaps	%%ymm3,0x320(%%rdi)		/* <-B3i = s1p24r */\n\t"\
			"vmovaps	%%ymm4,0x300(%%rdi)		/* <-B3r = s1p24r */\n\t"\
			"vmovaps	%%ymm5,0x420(%%rdi)		/* <-B4i = s1p04i */\n\t"\
			"\n\t"\
		"/* SSE2_RADIX_07_DFT(s1p21r,s1p17r,s1p13r,s1p09r,s1p05r,s1p01r,s1p25r,cc0,s1p01r,s1p05r,s1p09r,s1p13r,s1p17r,s1p21r,s1p25r): */\n\t"\
			"addq	$0x040,%%rdi	/* s1p01r */\n\t"\
			"vmovaps	0x400(%%rdi),%%ymm6			/* A1r = s1p17r */\n\t"\
			"vmovaps	0x600(%%rdi),%%ymm1			/* A6r = s1p25r */\n\t"\
			"vmovaps	0x300(%%rdi),%%ymm5			/* A2r = s1p13r */\n\t"\
			"vmovaps	      (%%rdi),%%ymm2			/* A5r = s1p01r */\n\t"\
			"vmovaps	0x200(%%rdi),%%ymm4			/* A3r = s1p09r */\n\t"\
			"vmovaps	0x100(%%rdi),%%ymm3			/* A4r = s1p05r */\n\t"\
			"\n\t"\
			"vsubpd	%%ymm1,%%ymm6,%%ymm6			/* t6r = A1r-A6r */\n\t"\
			"vaddpd	%%ymm1,%%ymm1,%%ymm1			/*         2*A6r */\n\t"\
			"vaddpd	%%ymm6,%%ymm1,%%ymm1			/* t1r = A1r+A6r */\n\t"\
			"\n\t"\
			"vsubpd	%%ymm2,%%ymm5,%%ymm5			/* t5r = A2r-A5r */\n\t"\
			"vaddpd	%%ymm2,%%ymm2,%%ymm2			/*         2*A5r */\n\t"\
			"vaddpd	%%ymm5,%%ymm2,%%ymm2			/* t2r = A2r+A5r */\n\t"\
			"\n\t"\
			"vmovaps	0x500(%%rdi),%%ymm0		/* Ar0 = s1p21r */\n\t"\
			"vsubpd	%%ymm3,%%ymm4,%%ymm4			/* t4r = A3r-A4r */\n\t"\
			"vaddpd	%%ymm3,%%ymm3,%%ymm3			/*         2*A4r */\n\t"\
			"vaddpd	%%ymm4,%%ymm3,%%ymm3			/* t3r = A3r+A4r */\n\t"\
			"\n\t"\
			"vmovaps	%%ymm0,     (%%rsi)		/* cpy t0 into scratch sincos slot */\n\t"\
			"vmovaps	%%ymm6,0x020(%%rsi)		/* cpy t6 into scratch sincos slot */\n\t"\
			"vaddpd	      %%ymm1,%%ymm0,%%ymm0		/*~A0 = A0+t1 */\n\t"\
			"vmovaps	%%ymm5,%%ymm7			/* cpy t5 */\n\t"\
			"vaddpd	      %%ymm2,%%ymm3,%%ymm3		/*~t3 = t3+t2 */\n\t"\
			"vsubpd	      %%ymm4,%%ymm5,%%ymm5		/*~t5 = t5-t4 */\n\t"\
			"vsubpd	      %%ymm2,%%ymm1,%%ymm1		/*~t1 = t1-t2 */\n\t"\
			"vsubpd	      %%ymm7,%%ymm6,%%ymm6		/*~t6 = t6-t5 */\n\t"\
			"vaddpd	      %%ymm2,%%ymm2,%%ymm2		/* 2*t2 */\n\t"\
			"vaddpd	      %%ymm7,%%ymm4,%%ymm4		/*~t5 = t4+t5 */\n\t"\
			"vaddpd	      %%ymm3,%%ymm0,%%ymm0		/* B0 */\n\t"\
			"vaddpd	0x020(%%rsi),%%ymm5,%%ymm5		/* t3 = (t5-t4)+t6 */\n\t"\
			"vsubpd	      %%ymm2,%%ymm3,%%ymm3		/*~t2 =  (t2+t3) - 2*t2 = t3-t2 */\n\t"\
			"vmovaps	%%ymm4,%%ymm7			/* cpy t5 */\n\t"\
			"vmovaps	%%ymm0,      (%%rdi)		/* <-B0 = s1p01r, ymm0 FREE */\n\t"\
			"vsubpd	      %%ymm6,%%ymm4,%%ymm4		/* t4 = ~t5-~t6 */\n\t"\
			"vmovaps	%%ymm1,%%ymm2			/* cpy ~t1 */\n\t"\
			"vsubpd	     (%%rsi),%%ymm0,%%ymm0		/* r = B0 - t0 */\n\t"\
			"vmulpd	0x020(%%rax),%%ymm5,%%ymm5		/*~t3 = t3*sx0 */\n\t"\
			"vaddpd	      %%ymm3,%%ymm2,%%ymm2		/* ~t1+~t2 */\n\t"\
			"vmulpd	     (%%rcx),%%ymm3,%%ymm3		/* t2 = t2*cx2 */\n\t"\
			"vmulpd	0x020(%%rdx),%%ymm4,%%ymm4		/*~t4 = t4*sx3 */\n\t"\
			"vmulpd	     (%%rbx),%%ymm1,%%ymm1		/* t1 = t1*cx1 */\n\t"\
			"vmulpd	0x020(%%rbx),%%ymm6,%%ymm6		/*~t6 = t6*sx1 */\n\t"\
			"vmulpd	     (%%rax),%%ymm0,%%ymm0		/* ~r = r*(cx0-1) */\n\t"\
			"vmulpd	0x020(%%rcx),%%ymm7,%%ymm7		/*~t5 = t5*sx2 */\n\t"\
			"vmulpd	     (%%rdx),%%ymm2,%%ymm2		/* t3 =(t1+t2)*cx3 */\n\t"\
			"vaddpd	     (%%rdi),%%ymm0,%%ymm0		/* t0 =~r + B0 */\n\t"\
			"vaddpd	      %%ymm4,%%ymm6,%%ymm6		/*~t6 = t4+t6 */\n\t"\
			"vsubpd	      %%ymm2,%%ymm1,%%ymm1		/* tt = t1-t3 */\n\t"\
			"vsubpd	      %%ymm7,%%ymm4,%%ymm4		/*~t5 = t4-t5, ymm7 FREE */\n\t"\
			"vsubpd	      %%ymm2,%%ymm3,%%ymm3		/* t2 = t2-t3, ymm2 FREE */\n\t"\
			"\n\t"\
			"vmovaps	%%ymm0,%%ymm2			/* cpy t0 */\n\t"\
			"vmovaps	%%ymm5,%%ymm7			/* cpy t3 */\n\t"\
			"vaddpd	      %%ymm1,%%ymm0,%%ymm0			/*~t0 = t0+tt */\n\t"\
			"vaddpd	      %%ymm6,%%ymm5,%%ymm5			/*~t3 = t3+t6 */\n\t"\
			"vaddpd	      %%ymm3,%%ymm1,%%ymm1			/*~tt = tt+t2 */\n\t"\
			"vaddpd	      %%ymm4,%%ymm6,%%ymm6			/*      t6+t5 */\n\t"\
			"vaddpd	      %%ymm2,%%ymm3,%%ymm3			/*~t2 = t2+t0 */\n\t"\
			"vaddpd	      %%ymm7,%%ymm4,%%ymm4			/*~t5 = t5+t3 */\n\t"\
			"vsubpd	      %%ymm1,%%ymm2,%%ymm2			/*~t1 = t0-tt-t2 */\n\t"\
			"vsubpd	      %%ymm6,%%ymm7,%%ymm7			/*~t4 = t3-t6-t5 */\n\t"\
			"vmovaps	%%ymm0,0x100(%%rdi)			/* B1 <- t0 = s1p05r */\n\t"\
			"vmovaps	%%ymm5,0x600(%%rdi)			/* B6 <- t3 = s1p25r */\n\t"\
			"vmovaps	%%ymm2,0x200(%%rdi)			/* B2 <- t1 = s1p09r */\n\t"\
			"vmovaps	%%ymm7,0x500(%%rdi)			/* B5 <- t4 = s1p21r */\n\t"\
			"vmovaps	%%ymm3,0x300(%%rdi)			/* B3 <- t2 = s1p13r */\n\t"\
			"vmovaps	%%ymm4,0x400(%%rdi)			/* B4 <- t5 = s1p17r */\n\t"\
			"/************************** Imaginary Parts: ******************************************/\n\t"\
			"vmovaps	0x420(%%rdi),%%ymm6		/* A1i = s1p17i */\n\t"\
			"vmovaps	0x620(%%rdi),%%ymm1		/* A6i = s1p25i */\n\t"\
			"vmovaps	0x320(%%rdi),%%ymm5		/* A2i = s1p13i */\n\t"\
			"vmovaps	0x020(%%rdi),%%ymm2		/* A5i = s1p01i */\n\t"\
			"vmovaps	0x220(%%rdi),%%ymm4		/* A3i = s1p09i */\n\t"\
			"vmovaps	0x120(%%rdi),%%ymm3		/* A4i = s1p05i */\n\t"\
			"\n\t"\
			"vsubpd	%%ymm1,%%ymm6,%%ymm6			/* t6i = A1i-A6i */\n\t"\
			"vaddpd	%%ymm1,%%ymm1,%%ymm1			/*         2*A6i */\n\t"\
			"vaddpd	%%ymm6,%%ymm1,%%ymm1			/* t1i = A1i+A6i */\n\t"\
			"\n\t"\
			"vsubpd	%%ymm2,%%ymm5,%%ymm5			/* t5i = A2i-A5i */\n\t"\
			"vaddpd	%%ymm2,%%ymm2,%%ymm2			/*         2*A5i */\n\t"\
			"vaddpd	%%ymm5,%%ymm2,%%ymm2			/* t2i = A2i+A5i */\n\t"\
			"\n\t"\
			"vmovaps	0x520(%%rdi),%%ymm0		/* Ai0 = s1p21i */\n\t"\
			"vsubpd	%%ymm3,%%ymm4,%%ymm4			/* t4i = A3i-A4i */\n\t"\
			"vaddpd	%%ymm3,%%ymm3,%%ymm3			/*         2*A4i */\n\t"\
			"vaddpd	%%ymm4,%%ymm3,%%ymm3			/* t3i = A3i+A4i */\n\t"\
			"\n\t"\
			"vmovaps	%%ymm0,     (%%rsi)		/* cpy t0 into scratch sincos slot */\n\t"\
			"vmovaps	%%ymm6,0x020(%%rsi)		/* cpy t6 into scratch sincos slot */\n\t"\
			"vaddpd	      %%ymm1,%%ymm0,%%ymm0			/*~A0 = A0+t1 */\n\t"\
			"vmovaps	%%ymm5,%%ymm7					/* cpy t5 */\n\t"\
			"vaddpd	      %%ymm2,%%ymm3,%%ymm3			/*~t3 = t3+t2 */\n\t"\
			"vsubpd	      %%ymm4,%%ymm5,%%ymm5			/*~t5 = t5-t4 */\n\t"\
			"vsubpd	      %%ymm2,%%ymm1,%%ymm1			/*~t1 = t1-t2 */\n\t"\
			"vsubpd	      %%ymm7,%%ymm6,%%ymm6			/*~t6 = t6-t5 */\n\t"\
			"vaddpd	      %%ymm2,%%ymm2,%%ymm2			/* 2*t2 */\n\t"\
			"vaddpd	      %%ymm7,%%ymm4,%%ymm4			/*~t5 = t4+t5 */\n\t"\
			"vaddpd	      %%ymm3,%%ymm0,%%ymm0			/* B0 */\n\t"\
			"vaddpd	0x020(%%rsi),%%ymm5,%%ymm5		/* t3 = (t5-t4)+t6 */\n\t"\
			"vsubpd	      %%ymm2,%%ymm3,%%ymm3			/*~t2 =  (t2+t3) - 2*t2 = t3-t2 */\n\t"\
			"vmovaps	%%ymm4,%%ymm7					/* cpy t5 */\n\t"\
			"vmovaps	%%ymm0,0x020(%%rdi)		/* <-B0 = s1p00i, ymm0 FREE */\n\t"\
			"vsubpd	      %%ymm6,%%ymm4,%%ymm4			/* t4 = ~t5-~t6 */\n\t"\
			"vmovaps	%%ymm1,%%ymm2					/* cpy ~t1 */\n\t"\
			"vsubpd	     (%%rsi),%%ymm0,%%ymm0		/* r = B0 - t0 */\n\t"\
			"vmulpd	0x020(%%rax),%%ymm5,%%ymm5		/*~t3 = t3*sx0 */\n\t"\
			"vaddpd	      %%ymm3,%%ymm2,%%ymm2			/* ~t1+~t2 */\n\t"\
			"vmulpd	     (%%rcx),%%ymm3,%%ymm3		/* t2 = t2*cx2 */\n\t"\
			"vmulpd	0x020(%%rdx),%%ymm4,%%ymm4		/*~t4 = t4*sx3 */\n\t"\
			"vmulpd	     (%%rbx),%%ymm1,%%ymm1		/* t1 = t1*cx1 */\n\t"\
			"vmulpd	0x020(%%rbx),%%ymm6,%%ymm6		/*~t6 = t6*sx1 */\n\t"\
			"vmulpd	     (%%rax),%%ymm0,%%ymm0		/* ~r = r*(cx0-1) */\n\t"\
			"vmulpd	0x020(%%rcx),%%ymm7,%%ymm7		/*~t5 = t5*sx2 */\n\t"\
			"vmulpd	     (%%rdx),%%ymm2,%%ymm2		/* t3 =(t1+t2)*cx3 */\n\t"\
			"vaddpd	0x020(%%rdi),%%ymm0,%%ymm0		/* t0 =~r + B0 */\n\t"\
			"vaddpd	      %%ymm4,%%ymm6,%%ymm6			/*~t6 = t4+t6 */\n\t"\
			"vsubpd	      %%ymm2,%%ymm1,%%ymm1			/* tt = t1-t3 */\n\t"\
			"vsubpd	      %%ymm7,%%ymm4,%%ymm4			/*~t5 = t4-t5, ymm7 FREE */\n\t"\
			"vsubpd	      %%ymm2,%%ymm3,%%ymm3			/* t2 = t2-t3, ymm2 FREE */\n\t"\
			"vmovaps	%%ymm0,%%ymm2					/* cpy t0 */\n\t"\
			"vmovaps	%%ymm5,%%ymm7					/* cpy t3 */\n\t"\
			"vaddpd	      %%ymm1,%%ymm0,%%ymm0			/*~t0 = t0+tt */\n\t"\
			"vaddpd	      %%ymm6,%%ymm5,%%ymm5			/*~t3 = t3+t6 */\n\t"\
			"vaddpd	      %%ymm3,%%ymm1,%%ymm1			/*~tt = tt+t2 */\n\t"\
			"vaddpd	      %%ymm4,%%ymm6,%%ymm6			/*      t6+t5 */\n\t"\
			"vaddpd	      %%ymm2,%%ymm3,%%ymm3			/*~t2 = t2+t0 */\n\t"\
			"vaddpd	      %%ymm7,%%ymm4,%%ymm4			/*~t5 = t5+t3 */\n\t"\
			"vsubpd	      %%ymm1,%%ymm2,%%ymm2			/*~t1 = t0-tt-t2, ymm1 FREE */\n\t"\
			"vsubpd	      %%ymm6,%%ymm7,%%ymm7			/*~t4 = t3-t6-t5, ymm6 FREE */\n\t"\
			"\n\t"\
			"/* ymm1,6 FREE */\n\t"\
			"vmovaps	0x100(%%rdi),%%ymm1		/* t0r */\n\t"\
			"vmovaps	0x600(%%rdi),%%ymm6		/* t3r */\n\t"\
			"vsubpd	%%ymm5,%%ymm1,%%ymm1			/* B1r =t0r-t3i */\n\t"\
			"vsubpd	%%ymm6,%%ymm0,%%ymm0			/* B6i =t0i-t3r */\n\t"\
			"vaddpd	%%ymm5,%%ymm5,%%ymm5			/*        2*t3i */\n\t"\
			"vaddpd	%%ymm6,%%ymm6,%%ymm6			/*        2*t3r */\n\t"\
			"vaddpd	%%ymm1,%%ymm5,%%ymm5			/* B6r =t0r+t3i */\n\t"\
			"vaddpd	%%ymm0,%%ymm6,%%ymm6			/* B1i =t0i+t3r */\n\t"\
			"vmovaps	%%ymm1,0x100(%%rdi)		/* <-B1r = s1p08r */\n\t"\
			"vmovaps	%%ymm0,0x620(%%rdi)		/* <-B6i = s1p20r */\n\t"\
			"vmovaps	%%ymm5,0x600(%%rdi)		/* <-B6r = s1p20r */\n\t"\
			"vmovaps	%%ymm6,0x120(%%rdi)		/* <-B1i = s1p08i */\n\t"\
			"\n\t"\
			"vmovaps	0x200(%%rdi),%%ymm1		/* t1r */\n\t"\
			"vmovaps	0x500(%%rdi),%%ymm6		/* t4r */\n\t"\
			"vsubpd	%%ymm7,%%ymm1,%%ymm1			/* B2r =t1r-t4i */\n\t"\
			"vsubpd	%%ymm6,%%ymm2,%%ymm2			/* B5i =t1i-t4r */\n\t"\
			"vaddpd	%%ymm7,%%ymm7,%%ymm7			/*        2*t4i */\n\t"\
			"vaddpd	%%ymm6,%%ymm6,%%ymm6			/*        2*t4r */\n\t"\
			"vaddpd	%%ymm1,%%ymm7,%%ymm7			/* B5r =t1r+t4i */\n\t"\
			"vaddpd	%%ymm2,%%ymm6,%%ymm6			/* B2i =t1i+t4r */\n\t"\
			"vmovaps	%%ymm1,0x200(%%rdi)		/* <-B2r = s1p16r */\n\t"\
			"vmovaps	%%ymm2,0x520(%%rdi)		/* <-B5i = s1p12r */\n\t"\
			"vmovaps	%%ymm7,0x500(%%rdi)		/* <-B5r = s1p12r */\n\t"\
			"vmovaps	%%ymm6,0x220(%%rdi)		/* <-B2i = s1p16i */\n\t"\
			"\n\t"\
			"/* Note the order reversal on this pair of outputs: */\n\t"\
			"vmovaps	0x300(%%rdi),%%ymm0		/* t2r */\n\t"\
			"vmovaps	0x400(%%rdi),%%ymm5		/* t5r */\n\t"\
			"vsubpd	%%ymm4,%%ymm0,%%ymm0			/* B4r =t2r-t5i */\n\t"\
			"vsubpd	%%ymm5,%%ymm3,%%ymm3			/* B3i =t2i-t5r */\n\t"\
			"vaddpd	%%ymm4,%%ymm4,%%ymm4			/*        2*t5i */\n\t"\
			"vaddpd	%%ymm5,%%ymm5,%%ymm5			/*        2*t5r */\n\t"\
			"vaddpd	%%ymm0,%%ymm4,%%ymm4			/* B3r =t2r+t5i */\n\t"\
			"vaddpd	%%ymm3,%%ymm5,%%ymm5			/* B4i =t2i+t5r */\n\t"\
			"vmovaps	%%ymm0,0x400(%%rdi)		/* <-B4r = s1p04r */\n\t"\
			"vmovaps	%%ymm3,0x320(%%rdi)		/* <-B3i = s1p24r */\n\t"\
			"vmovaps	%%ymm4,0x300(%%rdi)		/* <-B3r = s1p24r */\n\t"\
			"vmovaps	%%ymm5,0x420(%%rdi)		/* <-B4i = s1p04i */\n\t"\
			"\n\t"\
		"/* SSE2_RADIX_07_DFT(s1p14r,s1p10r,s1p06r,s1p02r,s1p26r,s1p22r,s1p18r,cc0,s1p02r,s1p06r,s1p10r,s1p14r,s1p18r,s1p22r,s1p26r): */\n\t"\
			"addq	$0x040,%%rdi	/* s1p02r */\n\t"\
			"vmovaps	0x200(%%rdi),%%ymm6			/* A1r = s1p10r */\n\t"\
			"vmovaps	0x400(%%rdi),%%ymm1			/* A6r = s1p18r */\n\t"\
			"vmovaps	0x100(%%rdi),%%ymm5			/* A2r = s1p06r */\n\t"\
			"vmovaps	0x500(%%rdi),%%ymm2			/* A5r = s1p22r */\n\t"\
			"vmovaps	      (%%rdi),%%ymm4			/* A3r = s1p02r */\n\t"\
			"vmovaps	0x600(%%rdi),%%ymm3			/* A4r = s1p26r */\n\t"\
			"\n\t"\
			"vsubpd	%%ymm1,%%ymm6,%%ymm6			/* t6r = A1r-A6r */\n\t"\
			"vaddpd	%%ymm1,%%ymm1,%%ymm1			/*         2*A6r */\n\t"\
			"vaddpd	%%ymm6,%%ymm1,%%ymm1			/* t1r = A1r+A6r */\n\t"\
			"\n\t"\
			"vsubpd	%%ymm2,%%ymm5,%%ymm5			/* t5r = A2r-A5r */\n\t"\
			"vaddpd	%%ymm2,%%ymm2,%%ymm2			/*         2*A5r */\n\t"\
			"vaddpd	%%ymm5,%%ymm2,%%ymm2			/* t2r = A2r+A5r */\n\t"\
			"\n\t"\
			"vmovaps	0x300(%%rdi),%%ymm0		/* Ar0 = s1p14r */\n\t"\
			"vsubpd	%%ymm3,%%ymm4,%%ymm4			/* t4r = A3r-A4r */\n\t"\
			"vaddpd	%%ymm3,%%ymm3,%%ymm3			/*         2*A4r */\n\t"\
			"vaddpd	%%ymm4,%%ymm3,%%ymm3			/* t3r = A3r+A4r */\n\t"\
			"\n\t"\
			"vmovaps	%%ymm0,     (%%rsi)		/* cpy t0 into scratch sincos slot */\n\t"\
			"vmovaps	%%ymm6,0x020(%%rsi)		/* cpy t6 into scratch sincos slot */\n\t"\
			"vaddpd	      %%ymm1,%%ymm0,%%ymm0		/*~A0 = A0+t1 */\n\t"\
			"vmovaps	%%ymm5,%%ymm7			/* cpy t5 */\n\t"\
			"vaddpd	      %%ymm2,%%ymm3,%%ymm3		/*~t3 = t3+t2 */\n\t"\
			"vsubpd	      %%ymm4,%%ymm5,%%ymm5		/*~t5 = t5-t4 */\n\t"\
			"vsubpd	      %%ymm2,%%ymm1,%%ymm1		/*~t1 = t1-t2 */\n\t"\
			"vsubpd	      %%ymm7,%%ymm6,%%ymm6		/*~t6 = t6-t5 */\n\t"\
			"vaddpd	      %%ymm2,%%ymm2,%%ymm2		/* 2*t2 */\n\t"\
			"vaddpd	      %%ymm7,%%ymm4,%%ymm4		/*~t5 = t4+t5 */\n\t"\
			"vaddpd	      %%ymm3,%%ymm0,%%ymm0		/* B0 */\n\t"\
			"vaddpd	0x020(%%rsi),%%ymm5,%%ymm5		/* t3 = (t5-t4)+t6 */\n\t"\
			"vsubpd	      %%ymm2,%%ymm3,%%ymm3		/*~t2 =  (t2+t3) - 2*t2 = t3-t2 */\n\t"\
			"vmovaps	%%ymm4,%%ymm7			/* cpy t5 */\n\t"\
			"vmovaps	%%ymm0,      (%%rdi)		/* <-B0 = s1p01r, ymm0 FREE */\n\t"\
			"vsubpd	      %%ymm6,%%ymm4,%%ymm4		/* t4 = ~t5-~t6 */\n\t"\
			"vmovaps	%%ymm1,%%ymm2			/* cpy ~t1 */\n\t"\
			"vsubpd	     (%%rsi),%%ymm0,%%ymm0		/* r = B0 - t0 */\n\t"\
			"vmulpd	0x020(%%rax),%%ymm5,%%ymm5		/*~t3 = t3*sx0 */\n\t"\
			"vaddpd	      %%ymm3,%%ymm2,%%ymm2		/* ~t1+~t2 */\n\t"\
			"vmulpd	     (%%rcx),%%ymm3,%%ymm3		/* t2 = t2*cx2 */\n\t"\
			"vmulpd	0x020(%%rdx),%%ymm4,%%ymm4		/*~t4 = t4*sx3 */\n\t"\
			"vmulpd	     (%%rbx),%%ymm1,%%ymm1		/* t1 = t1*cx1 */\n\t"\
			"vmulpd	0x020(%%rbx),%%ymm6,%%ymm6		/*~t6 = t6*sx1 */\n\t"\
			"vmulpd	     (%%rax),%%ymm0,%%ymm0		/* ~r = r*(cx0-1) */\n\t"\
			"vmulpd	0x020(%%rcx),%%ymm7,%%ymm7		/*~t5 = t5*sx2 */\n\t"\
			"vmulpd	     (%%rdx),%%ymm2,%%ymm2		/* t3 =(t1+t2)*cx3 */\n\t"\
			"vaddpd	     (%%rdi),%%ymm0,%%ymm0		/* t0 =~r + B0 */\n\t"\
			"vaddpd	      %%ymm4,%%ymm6,%%ymm6		/*~t6 = t4+t6 */\n\t"\
			"vsubpd	      %%ymm2,%%ymm1,%%ymm1		/* tt = t1-t3 */\n\t"\
			"vsubpd	      %%ymm7,%%ymm4,%%ymm4		/*~t5 = t4-t5, ymm7 FREE */\n\t"\
			"vsubpd	      %%ymm2,%%ymm3,%%ymm3		/* t2 = t2-t3, ymm2 FREE */\n\t"\
			"\n\t"\
			"vmovaps	%%ymm0,%%ymm2			/* cpy t0 */\n\t"\
			"vmovaps	%%ymm5,%%ymm7			/* cpy t3 */\n\t"\
			"vaddpd	      %%ymm1,%%ymm0,%%ymm0			/*~t0 = t0+tt */\n\t"\
			"vaddpd	      %%ymm6,%%ymm5,%%ymm5			/*~t3 = t3+t6 */\n\t"\
			"vaddpd	      %%ymm3,%%ymm1,%%ymm1			/*~tt = tt+t2 */\n\t"\
			"vaddpd	      %%ymm4,%%ymm6,%%ymm6			/*      t6+t5 */\n\t"\
			"vaddpd	      %%ymm2,%%ymm3,%%ymm3			/*~t2 = t2+t0 */\n\t"\
			"vaddpd	      %%ymm7,%%ymm4,%%ymm4			/*~t5 = t5+t3 */\n\t"\
			"vsubpd	      %%ymm1,%%ymm2,%%ymm2			/*~t1 = t0-tt-t2 */\n\t"\
			"vsubpd	      %%ymm6,%%ymm7,%%ymm7			/*~t4 = t3-t6-t5 */\n\t"\
			"vmovaps	%%ymm0,0x100(%%rdi)			/* B1 <- t0 = s1p05r */\n\t"\
			"vmovaps	%%ymm5,0x600(%%rdi)			/* B6 <- t3 = s1p25r */\n\t"\
			"vmovaps	%%ymm2,0x200(%%rdi)			/* B2 <- t1 = s1p09r */\n\t"\
			"vmovaps	%%ymm7,0x500(%%rdi)			/* B5 <- t4 = s1p21r */\n\t"\
			"vmovaps	%%ymm3,0x300(%%rdi)			/* B3 <- t2 = s1p13r */\n\t"\
			"vmovaps	%%ymm4,0x400(%%rdi)			/* B4 <- t5 = s1p17r */\n\t"\
			"/************************** Imaginary Parts: ******************************************/\n\t"\
			"vmovaps	0x220(%%rdi),%%ymm6		/* A1i = s1p10i */\n\t"\
			"vmovaps	0x420(%%rdi),%%ymm1		/* A6i = s1p18i */\n\t"\
			"vmovaps	0x120(%%rdi),%%ymm5		/* A2i = s1p06i */\n\t"\
			"vmovaps	0x520(%%rdi),%%ymm2		/* A5i = s1p22i */\n\t"\
			"vmovaps	0x020(%%rdi),%%ymm4		/* A3i = s1p02i */\n\t"\
			"vmovaps	0x620(%%rdi),%%ymm3		/* A4i = s1p26i */\n\t"\
			"\n\t"\
			"vsubpd	%%ymm1,%%ymm6,%%ymm6			/* t6i = A1i-A6i */\n\t"\
			"vaddpd	%%ymm1,%%ymm1,%%ymm1			/*         2*A6i */\n\t"\
			"vaddpd	%%ymm6,%%ymm1,%%ymm1			/* t1i = A1i+A6i */\n\t"\
			"\n\t"\
			"vsubpd	%%ymm2,%%ymm5,%%ymm5			/* t5i = A2i-A5i */\n\t"\
			"vaddpd	%%ymm2,%%ymm2,%%ymm2			/*         2*A5i */\n\t"\
			"vaddpd	%%ymm5,%%ymm2,%%ymm2			/* t2i = A2i+A5i */\n\t"\
			"\n\t"\
			"vmovaps	0x320(%%rdi),%%ymm0		/* Ai0 = s1p14i */\n\t"\
			"vsubpd	%%ymm3,%%ymm4,%%ymm4			/* t4i = A3i-A4i */\n\t"\
			"vaddpd	%%ymm3,%%ymm3,%%ymm3			/*         2*A4i */\n\t"\
			"vaddpd	%%ymm4,%%ymm3,%%ymm3			/* t3i = A3i+A4i */\n\t"\
			"\n\t"\
			"vmovaps	%%ymm0,     (%%rsi)		/* cpy t0 into scratch sincos slot */\n\t"\
			"vmovaps	%%ymm6,0x020(%%rsi)		/* cpy t6 into scratch sincos slot */\n\t"\
			"vaddpd	      %%ymm1,%%ymm0,%%ymm0			/*~A0 = A0+t1 */\n\t"\
			"vmovaps	%%ymm5,%%ymm7					/* cpy t5 */\n\t"\
			"vaddpd	      %%ymm2,%%ymm3,%%ymm3			/*~t3 = t3+t2 */\n\t"\
			"vsubpd	      %%ymm4,%%ymm5,%%ymm5			/*~t5 = t5-t4 */\n\t"\
			"vsubpd	      %%ymm2,%%ymm1,%%ymm1			/*~t1 = t1-t2 */\n\t"\
			"vsubpd	      %%ymm7,%%ymm6,%%ymm6			/*~t6 = t6-t5 */\n\t"\
			"vaddpd	      %%ymm2,%%ymm2,%%ymm2			/* 2*t2 */\n\t"\
			"vaddpd	      %%ymm7,%%ymm4,%%ymm4			/*~t5 = t4+t5 */\n\t"\
			"vaddpd	      %%ymm3,%%ymm0,%%ymm0			/* B0 */\n\t"\
			"vaddpd	0x020(%%rsi),%%ymm5,%%ymm5		/* t3 = (t5-t4)+t6 */\n\t"\
			"vsubpd	      %%ymm2,%%ymm3,%%ymm3			/*~t2 =  (t2+t3) - 2*t2 = t3-t2 */\n\t"\
			"vmovaps	%%ymm4,%%ymm7					/* cpy t5 */\n\t"\
			"vmovaps	%%ymm0,0x020(%%rdi)		/* <-B0 = s1p00i, ymm0 FREE */\n\t"\
			"vsubpd	      %%ymm6,%%ymm4,%%ymm4			/* t4 = ~t5-~t6 */\n\t"\
			"vmovaps	%%ymm1,%%ymm2					/* cpy ~t1 */\n\t"\
			"vsubpd	     (%%rsi),%%ymm0,%%ymm0		/* r = B0 - t0 */\n\t"\
			"vmulpd	0x020(%%rax),%%ymm5,%%ymm5		/*~t3 = t3*sx0 */\n\t"\
			"vaddpd	      %%ymm3,%%ymm2,%%ymm2			/* ~t1+~t2 */\n\t"\
			"vmulpd	     (%%rcx),%%ymm3,%%ymm3		/* t2 = t2*cx2 */\n\t"\
			"vmulpd	0x020(%%rdx),%%ymm4,%%ymm4		/*~t4 = t4*sx3 */\n\t"\
			"vmulpd	     (%%rbx),%%ymm1,%%ymm1		/* t1 = t1*cx1 */\n\t"\
			"vmulpd	0x020(%%rbx),%%ymm6,%%ymm6		/*~t6 = t6*sx1 */\n\t"\
			"vmulpd	     (%%rax),%%ymm0,%%ymm0		/* ~r = r*(cx0-1) */\n\t"\
			"vmulpd	0x020(%%rcx),%%ymm7,%%ymm7		/*~t5 = t5*sx2 */\n\t"\
			"vmulpd	     (%%rdx),%%ymm2,%%ymm2		/* t3 =(t1+t2)*cx3 */\n\t"\
			"vaddpd	0x020(%%rdi),%%ymm0,%%ymm0		/* t0 =~r + B0 */\n\t"\
			"vaddpd	      %%ymm4,%%ymm6,%%ymm6			/*~t6 = t4+t6 */\n\t"\
			"vsubpd	      %%ymm2,%%ymm1,%%ymm1			/* tt = t1-t3 */\n\t"\
			"vsubpd	      %%ymm7,%%ymm4,%%ymm4			/*~t5 = t4-t5, ymm7 FREE */\n\t"\
			"vsubpd	      %%ymm2,%%ymm3,%%ymm3			/* t2 = t2-t3, ymm2 FREE */\n\t"\
			"vmovaps	%%ymm0,%%ymm2					/* cpy t0 */\n\t"\
			"vmovaps	%%ymm5,%%ymm7					/* cpy t3 */\n\t"\
			"vaddpd	      %%ymm1,%%ymm0,%%ymm0			/*~t0 = t0+tt */\n\t"\
			"vaddpd	      %%ymm6,%%ymm5,%%ymm5			/*~t3 = t3+t6 */\n\t"\
			"vaddpd	      %%ymm3,%%ymm1,%%ymm1			/*~tt = tt+t2 */\n\t"\
			"vaddpd	      %%ymm4,%%ymm6,%%ymm6			/*      t6+t5 */\n\t"\
			"vaddpd	      %%ymm2,%%ymm3,%%ymm3			/*~t2 = t2+t0 */\n\t"\
			"vaddpd	      %%ymm7,%%ymm4,%%ymm4			/*~t5 = t5+t3 */\n\t"\
			"vsubpd	      %%ymm1,%%ymm2,%%ymm2			/*~t1 = t0-tt-t2, ymm1 FREE */\n\t"\
			"vsubpd	      %%ymm6,%%ymm7,%%ymm7			/*~t4 = t3-t6-t5, ymm6 FREE */\n\t"\
			"\n\t"\
			"/* ymm1,6 FREE */\n\t"\
			"vmovaps	0x100(%%rdi),%%ymm1		/* t0r */\n\t"\
			"vmovaps	0x600(%%rdi),%%ymm6		/* t3r */\n\t"\
			"vsubpd	%%ymm5,%%ymm1,%%ymm1			/* B1r =t0r-t3i */\n\t"\
			"vsubpd	%%ymm6,%%ymm0,%%ymm0			/* B6i =t0i-t3r */\n\t"\
			"vaddpd	%%ymm5,%%ymm5,%%ymm5			/*        2*t3i */\n\t"\
			"vaddpd	%%ymm6,%%ymm6,%%ymm6			/*        2*t3r */\n\t"\
			"vaddpd	%%ymm1,%%ymm5,%%ymm5			/* B6r =t0r+t3i */\n\t"\
			"vaddpd	%%ymm0,%%ymm6,%%ymm6			/* B1i =t0i+t3r */\n\t"\
			"vmovaps	%%ymm1,0x100(%%rdi)		/* <-B1r = s1p08r */\n\t"\
			"vmovaps	%%ymm0,0x620(%%rdi)		/* <-B6i = s1p20r */\n\t"\
			"vmovaps	%%ymm5,0x600(%%rdi)		/* <-B6r = s1p20r */\n\t"\
			"vmovaps	%%ymm6,0x120(%%rdi)		/* <-B1i = s1p08i */\n\t"\
			"\n\t"\
			"vmovaps	0x200(%%rdi),%%ymm1		/* t1r */\n\t"\
			"vmovaps	0x500(%%rdi),%%ymm6		/* t4r */\n\t"\
			"vsubpd	%%ymm7,%%ymm1,%%ymm1			/* B2r =t1r-t4i */\n\t"\
			"vsubpd	%%ymm6,%%ymm2,%%ymm2			/* B5i =t1i-t4r */\n\t"\
			"vaddpd	%%ymm7,%%ymm7,%%ymm7			/*        2*t4i */\n\t"\
			"vaddpd	%%ymm6,%%ymm6,%%ymm6			/*        2*t4r */\n\t"\
			"vaddpd	%%ymm1,%%ymm7,%%ymm7			/* B5r =t1r+t4i */\n\t"\
			"vaddpd	%%ymm2,%%ymm6,%%ymm6			/* B2i =t1i+t4r */\n\t"\
			"vmovaps	%%ymm1,0x200(%%rdi)		/* <-B2r = s1p16r */\n\t"\
			"vmovaps	%%ymm2,0x520(%%rdi)		/* <-B5i = s1p12r */\n\t"\
			"vmovaps	%%ymm7,0x500(%%rdi)		/* <-B5r = s1p12r */\n\t"\
			"vmovaps	%%ymm6,0x220(%%rdi)		/* <-B2i = s1p16i */\n\t"\
			"\n\t"\
			"/* Note the order reversal on this pair of outputs: */\n\t"\
			"vmovaps	0x300(%%rdi),%%ymm0		/* t2r */\n\t"\
			"vmovaps	0x400(%%rdi),%%ymm5		/* t5r */\n\t"\
			"vsubpd	%%ymm4,%%ymm0,%%ymm0			/* B4r =t2r-t5i */\n\t"\
			"vsubpd	%%ymm5,%%ymm3,%%ymm3			/* B3i =t2i-t5r */\n\t"\
			"vaddpd	%%ymm4,%%ymm4,%%ymm4			/*        2*t5i */\n\t"\
			"vaddpd	%%ymm5,%%ymm5,%%ymm5			/*        2*t5r */\n\t"\
			"vaddpd	%%ymm0,%%ymm4,%%ymm4			/* B3r =t2r+t5i */\n\t"\
			"vaddpd	%%ymm3,%%ymm5,%%ymm5			/* B4i =t2i+t5r */\n\t"\
			"vmovaps	%%ymm0,0x400(%%rdi)		/* <-B4r = s1p04r */\n\t"\
			"vmovaps	%%ymm3,0x320(%%rdi)		/* <-B3i = s1p24r */\n\t"\
			"vmovaps	%%ymm4,0x300(%%rdi)		/* <-B3r = s1p24r */\n\t"\
			"vmovaps	%%ymm5,0x420(%%rdi)		/* <-B4i = s1p04i */\n\t"\
			"\n\t"\
		"/* SSE2_RADIX_07_DFT(s1p07r,s1p03r,s1p27r,s1p23r,s1p19r,s1p15r,s1p11r,cc0,s1p03r,s1p07r,s1p11r,s1p15r,s1p19r,s1p23r,s1p27r): */\n\t"\
			"addq	$0x040,%%rdi	/* s1p03r */\n\t"\
			"vmovaps	      (%%rdi),%%ymm6			/* A1r = s1p03r */\n\t"\
			"vmovaps	0x200(%%rdi),%%ymm1			/* A6r = s1p11r */\n\t"\
			"vmovaps	0x600(%%rdi),%%ymm5			/* A2r = s1p27r */\n\t"\
			"vmovaps	0x300(%%rdi),%%ymm2			/* A5r = s1p15r */\n\t"\
			"vmovaps	0x500(%%rdi),%%ymm4			/* A3r = s1p23r */\n\t"\
			"vmovaps	0x400(%%rdi),%%ymm3			/* A4r = s1p19r */\n\t"\
			"\n\t"\
			"vsubpd	%%ymm1,%%ymm6,%%ymm6			/* t6r = A1r-A6r */\n\t"\
			"vaddpd	%%ymm1,%%ymm1,%%ymm1			/*         2*A6r */\n\t"\
			"vaddpd	%%ymm6,%%ymm1,%%ymm1			/* t1r = A1r+A6r */\n\t"\
			"\n\t"\
			"vsubpd	%%ymm2,%%ymm5,%%ymm5			/* t5r = A2r-A5r */\n\t"\
			"vaddpd	%%ymm2,%%ymm2,%%ymm2			/*         2*A5r */\n\t"\
			"vaddpd	%%ymm5,%%ymm2,%%ymm2			/* t2r = A2r+A5r */\n\t"\
			"\n\t"\
			"vmovaps	0x100(%%rdi),%%ymm0		/* Ar0 = s1p07r */\n\t"\
			"vsubpd	%%ymm3,%%ymm4,%%ymm4			/* t4r = A3r-A4r */\n\t"\
			"vaddpd	%%ymm3,%%ymm3,%%ymm3			/*         2*A4r */\n\t"\
			"vaddpd	%%ymm4,%%ymm3,%%ymm3			/* t3r = A3r+A4r */\n\t"\
			"\n\t"\
			"vmovaps	%%ymm0,     (%%rsi)		/* cpy t0 into scratch sincos slot */\n\t"\
			"vmovaps	%%ymm6,0x020(%%rsi)		/* cpy t6 into scratch sincos slot */\n\t"\
			"vaddpd	      %%ymm1,%%ymm0,%%ymm0		/*~A0 = A0+t1 */\n\t"\
			"vmovaps	%%ymm5,%%ymm7			/* cpy t5 */\n\t"\
			"vaddpd	      %%ymm2,%%ymm3,%%ymm3		/*~t3 = t3+t2 */\n\t"\
			"vsubpd	      %%ymm4,%%ymm5,%%ymm5		/*~t5 = t5-t4 */\n\t"\
			"vsubpd	      %%ymm2,%%ymm1,%%ymm1		/*~t1 = t1-t2 */\n\t"\
			"vsubpd	      %%ymm7,%%ymm6,%%ymm6		/*~t6 = t6-t5 */\n\t"\
			"vaddpd	      %%ymm2,%%ymm2,%%ymm2		/* 2*t2 */\n\t"\
			"vaddpd	      %%ymm7,%%ymm4,%%ymm4		/*~t5 = t4+t5 */\n\t"\
			"vaddpd	      %%ymm3,%%ymm0,%%ymm0		/* B0 */\n\t"\
			"vaddpd	0x020(%%rsi),%%ymm5,%%ymm5		/* t3 = (t5-t4)+t6 */\n\t"\
			"vsubpd	      %%ymm2,%%ymm3,%%ymm3		/*~t2 =  (t2+t3) - 2*t2 = t3-t2 */\n\t"\
			"vmovaps	%%ymm4,%%ymm7			/* cpy t5 */\n\t"\
			"vmovaps	%%ymm0,      (%%rdi)		/* <-B0 = s1p01r, ymm0 FREE */\n\t"\
			"vsubpd	      %%ymm6,%%ymm4,%%ymm4		/* t4 = ~t5-~t6 */\n\t"\
			"vmovaps	%%ymm1,%%ymm2			/* cpy ~t1 */\n\t"\
			"vsubpd	     (%%rsi),%%ymm0,%%ymm0		/* r = B0 - t0 */\n\t"\
			"vmulpd	0x020(%%rax),%%ymm5,%%ymm5		/*~t3 = t3*sx0 */\n\t"\
			"vaddpd	      %%ymm3,%%ymm2,%%ymm2		/* ~t1+~t2 */\n\t"\
			"vmulpd	     (%%rcx),%%ymm3,%%ymm3		/* t2 = t2*cx2 */\n\t"\
			"vmulpd	0x020(%%rdx),%%ymm4,%%ymm4		/*~t4 = t4*sx3 */\n\t"\
			"vmulpd	     (%%rbx),%%ymm1,%%ymm1		/* t1 = t1*cx1 */\n\t"\
			"vmulpd	0x020(%%rbx),%%ymm6,%%ymm6		/*~t6 = t6*sx1 */\n\t"\
			"vmulpd	     (%%rax),%%ymm0,%%ymm0		/* ~r = r*(cx0-1) */\n\t"\
			"vmulpd	0x020(%%rcx),%%ymm7,%%ymm7		/*~t5 = t5*sx2 */\n\t"\
			"vmulpd	     (%%rdx),%%ymm2,%%ymm2		/* t3 =(t1+t2)*cx3 */\n\t"\
			"vaddpd	     (%%rdi),%%ymm0,%%ymm0		/* t0 =~r + B0 */\n\t"\
			"vaddpd	      %%ymm4,%%ymm6,%%ymm6		/*~t6 = t4+t6 */\n\t"\
			"vsubpd	      %%ymm2,%%ymm1,%%ymm1		/* tt = t1-t3 */\n\t"\
			"vsubpd	      %%ymm7,%%ymm4,%%ymm4		/*~t5 = t4-t5, ymm7 FREE */\n\t"\
			"vsubpd	      %%ymm2,%%ymm3,%%ymm3		/* t2 = t2-t3, ymm2 FREE */\n\t"\
			"\n\t"\
			"vmovaps	%%ymm0,%%ymm2			/* cpy t0 */\n\t"\
			"vmovaps	%%ymm5,%%ymm7			/* cpy t3 */\n\t"\
			"vaddpd	      %%ymm1,%%ymm0,%%ymm0			/*~t0 = t0+tt */\n\t"\
			"vaddpd	      %%ymm6,%%ymm5,%%ymm5			/*~t3 = t3+t6 */\n\t"\
			"vaddpd	      %%ymm3,%%ymm1,%%ymm1			/*~tt = tt+t2 */\n\t"\
			"vaddpd	      %%ymm4,%%ymm6,%%ymm6			/*      t6+t5 */\n\t"\
			"vaddpd	      %%ymm2,%%ymm3,%%ymm3			/*~t2 = t2+t0 */\n\t"\
			"vaddpd	      %%ymm7,%%ymm4,%%ymm4			/*~t5 = t5+t3 */\n\t"\
			"vsubpd	      %%ymm1,%%ymm2,%%ymm2			/*~t1 = t0-tt-t2 */\n\t"\
			"vsubpd	      %%ymm6,%%ymm7,%%ymm7			/*~t4 = t3-t6-t5 */\n\t"\
			"vmovaps	%%ymm0,0x100(%%rdi)			/* B1 <- t0 = s1p05r */\n\t"\
			"vmovaps	%%ymm5,0x600(%%rdi)			/* B6 <- t3 = s1p25r */\n\t"\
			"vmovaps	%%ymm2,0x200(%%rdi)			/* B2 <- t1 = s1p09r */\n\t"\
			"vmovaps	%%ymm7,0x500(%%rdi)			/* B5 <- t4 = s1p21r */\n\t"\
			"vmovaps	%%ymm3,0x300(%%rdi)			/* B3 <- t2 = s1p13r */\n\t"\
			"vmovaps	%%ymm4,0x400(%%rdi)			/* B4 <- t5 = s1p17r */\n\t"\
			"/************************** Imaginary Parts: ******************************************/\n\t"\
			"vmovaps	0x020(%%rdi),%%ymm6			/* A1i = s1p03i */\n\t"\
			"vmovaps	0x220(%%rdi),%%ymm1			/* A6i = s1p11i */\n\t"\
			"vmovaps	0x620(%%rdi),%%ymm5			/* A2i = s1p27i */\n\t"\
			"vmovaps	0x320(%%rdi),%%ymm2			/* A5i = s1p15i */\n\t"\
			"vmovaps	0x520(%%rdi),%%ymm4			/* A3i = s1p23i */\n\t"\
			"vmovaps	0x420(%%rdi),%%ymm3			/* A4i = s1p19i */\n\t"\
			"\n\t"\
			"vsubpd	%%ymm1,%%ymm6,%%ymm6			/* t6i = A1i-A6i */\n\t"\
			"vaddpd	%%ymm1,%%ymm1,%%ymm1			/*         2*A6i */\n\t"\
			"vaddpd	%%ymm6,%%ymm1,%%ymm1			/* t1i = A1i+A6i */\n\t"\
			"\n\t"\
			"vsubpd	%%ymm2,%%ymm5,%%ymm5			/* t5i = A2i-A5i */\n\t"\
			"vaddpd	%%ymm2,%%ymm2,%%ymm2			/*         2*A5i */\n\t"\
			"vaddpd	%%ymm5,%%ymm2,%%ymm2			/* t2i = A2i+A5i */\n\t"\
			"\n\t"\
			"vmovaps	0x120(%%rdi),%%ymm0		/* Ai0 = s1p07i */\n\t"\
			"vsubpd	%%ymm3,%%ymm4,%%ymm4			/* t4i = A3i-A4i */\n\t"\
			"vaddpd	%%ymm3,%%ymm3,%%ymm3			/*         2*A4i */\n\t"\
			"vaddpd	%%ymm4,%%ymm3,%%ymm3			/* t3i = A3i+A4i */\n\t"\
			"\n\t"\
			"vmovaps	%%ymm0,     (%%rsi)		/* cpy t0 into scratch sincos slot */\n\t"\
			"vmovaps	%%ymm6,0x020(%%rsi)		/* cpy t6 into scratch sincos slot */\n\t"\
			"vaddpd	      %%ymm1,%%ymm0,%%ymm0			/*~A0 = A0+t1 */\n\t"\
			"vmovaps	%%ymm5,%%ymm7					/* cpy t5 */\n\t"\
			"vaddpd	      %%ymm2,%%ymm3,%%ymm3			/*~t3 = t3+t2 */\n\t"\
			"vsubpd	      %%ymm4,%%ymm5,%%ymm5			/*~t5 = t5-t4 */\n\t"\
			"vsubpd	      %%ymm2,%%ymm1,%%ymm1			/*~t1 = t1-t2 */\n\t"\
			"vsubpd	      %%ymm7,%%ymm6,%%ymm6			/*~t6 = t6-t5 */\n\t"\
			"vaddpd	      %%ymm2,%%ymm2,%%ymm2			/* 2*t2 */\n\t"\
			"vaddpd	      %%ymm7,%%ymm4,%%ymm4			/*~t5 = t4+t5 */\n\t"\
			"vaddpd	      %%ymm3,%%ymm0,%%ymm0			/* B0 */\n\t"\
			"vaddpd	0x020(%%rsi),%%ymm5,%%ymm5		/* t3 = (t5-t4)+t6 */\n\t"\
			"vsubpd	      %%ymm2,%%ymm3,%%ymm3			/*~t2 =  (t2+t3) - 2*t2 = t3-t2 */\n\t"\
			"vmovaps	%%ymm4,%%ymm7					/* cpy t5 */\n\t"\
			"vmovaps	%%ymm0,0x020(%%rdi)		/* <-B0 = s1p00i, ymm0 FREE */\n\t"\
			"vsubpd	      %%ymm6,%%ymm4,%%ymm4			/* t4 = ~t5-~t6 */\n\t"\
			"vmovaps	%%ymm1,%%ymm2					/* cpy ~t1 */\n\t"\
			"vsubpd	     (%%rsi),%%ymm0,%%ymm0		/* r = B0 - t0 */\n\t"\
			"vmulpd	0x020(%%rax),%%ymm5,%%ymm5		/*~t3 = t3*sx0 */\n\t"\
			"vaddpd	      %%ymm3,%%ymm2,%%ymm2			/* ~t1+~t2 */\n\t"\
			"vmulpd	     (%%rcx),%%ymm3,%%ymm3		/* t2 = t2*cx2 */\n\t"\
			"vmulpd	0x020(%%rdx),%%ymm4,%%ymm4		/*~t4 = t4*sx3 */\n\t"\
			"vmulpd	     (%%rbx),%%ymm1,%%ymm1		/* t1 = t1*cx1 */\n\t"\
			"vmulpd	0x020(%%rbx),%%ymm6,%%ymm6		/*~t6 = t6*sx1 */\n\t"\
			"vmulpd	     (%%rax),%%ymm0,%%ymm0		/* ~r = r*(cx0-1) */\n\t"\
			"vmulpd	0x020(%%rcx),%%ymm7,%%ymm7		/*~t5 = t5*sx2 */\n\t"\
			"vmulpd	     (%%rdx),%%ymm2,%%ymm2		/* t3 =(t1+t2)*cx3 */\n\t"\
			"vaddpd	0x020(%%rdi),%%ymm0,%%ymm0		/* t0 =~r + B0 */\n\t"\
			"vaddpd	      %%ymm4,%%ymm6,%%ymm6			/*~t6 = t4+t6 */\n\t"\
			"vsubpd	      %%ymm2,%%ymm1,%%ymm1			/* tt = t1-t3 */\n\t"\
			"vsubpd	      %%ymm7,%%ymm4,%%ymm4			/*~t5 = t4-t5, ymm7 FREE */\n\t"\
			"vsubpd	      %%ymm2,%%ymm3,%%ymm3			/* t2 = t2-t3, ymm2 FREE */\n\t"\
			"vmovaps	%%ymm0,%%ymm2					/* cpy t0 */\n\t"\
			"vmovaps	%%ymm5,%%ymm7					/* cpy t3 */\n\t"\
			"vaddpd	      %%ymm1,%%ymm0,%%ymm0			/*~t0 = t0+tt */\n\t"\
			"vaddpd	      %%ymm6,%%ymm5,%%ymm5			/*~t3 = t3+t6 */\n\t"\
			"vaddpd	      %%ymm3,%%ymm1,%%ymm1			/*~tt = tt+t2 */\n\t"\
			"vaddpd	      %%ymm4,%%ymm6,%%ymm6			/*      t6+t5 */\n\t"\
			"vaddpd	      %%ymm2,%%ymm3,%%ymm3			/*~t2 = t2+t0 */\n\t"\
			"vaddpd	      %%ymm7,%%ymm4,%%ymm4			/*~t5 = t5+t3 */\n\t"\
			"vsubpd	      %%ymm1,%%ymm2,%%ymm2			/*~t1 = t0-tt-t2, ymm1 FREE */\n\t"\
			"vsubpd	      %%ymm6,%%ymm7,%%ymm7			/*~t4 = t3-t6-t5, ymm6 FREE */\n\t"\
			"\n\t"\
			"/* ymm1,6 FREE */\n\t"\
			"vmovaps	0x100(%%rdi),%%ymm1		/* t0r */\n\t"\
			"vmovaps	0x600(%%rdi),%%ymm6		/* t3r */\n\t"\
			"vsubpd	%%ymm5,%%ymm1,%%ymm1			/* B1r =t0r-t3i */\n\t"\
			"vsubpd	%%ymm6,%%ymm0,%%ymm0			/* B6i =t0i-t3r */\n\t"\
			"vaddpd	%%ymm5,%%ymm5,%%ymm5			/*        2*t3i */\n\t"\
			"vaddpd	%%ymm6,%%ymm6,%%ymm6			/*        2*t3r */\n\t"\
			"vaddpd	%%ymm1,%%ymm5,%%ymm5			/* B6r =t0r+t3i */\n\t"\
			"vaddpd	%%ymm0,%%ymm6,%%ymm6			/* B1i =t0i+t3r */\n\t"\
			"vmovaps	%%ymm1,0x100(%%rdi)		/* <-B1r = s1p08r */\n\t"\
			"vmovaps	%%ymm0,0x620(%%rdi)		/* <-B6i = s1p20r */\n\t"\
			"vmovaps	%%ymm5,0x600(%%rdi)		/* <-B6r = s1p20r */\n\t"\
			"vmovaps	%%ymm6,0x120(%%rdi)		/* <-B1i = s1p08i */\n\t"\
			"\n\t"\
			"vmovaps	0x200(%%rdi),%%ymm1		/* t1r */\n\t"\
			"vmovaps	0x500(%%rdi),%%ymm6		/* t4r */\n\t"\
			"vsubpd	%%ymm7,%%ymm1,%%ymm1			/* B2r =t1r-t4i */\n\t"\
			"vsubpd	%%ymm6,%%ymm2,%%ymm2			/* B5i =t1i-t4r */\n\t"\
			"vaddpd	%%ymm7,%%ymm7,%%ymm7			/*        2*t4i */\n\t"\
			"vaddpd	%%ymm6,%%ymm6,%%ymm6			/*        2*t4r */\n\t"\
			"vaddpd	%%ymm1,%%ymm7,%%ymm7			/* B5r =t1r+t4i */\n\t"\
			"vaddpd	%%ymm2,%%ymm6,%%ymm6			/* B2i =t1i+t4r */\n\t"\
			"vmovaps	%%ymm1,0x200(%%rdi)		/* <-B2r = s1p16r */\n\t"\
			"vmovaps	%%ymm2,0x520(%%rdi)		/* <-B5i = s1p12r */\n\t"\
			"vmovaps	%%ymm7,0x500(%%rdi)		/* <-B5r = s1p12r */\n\t"\
			"vmovaps	%%ymm6,0x220(%%rdi)		/* <-B2i = s1p16i */\n\t"\
			"\n\t"\
			"/* Note the order reversal on this pair of outputs: */\n\t"\
			"vmovaps	0x300(%%rdi),%%ymm0		/* t2r */\n\t"\
			"vmovaps	0x400(%%rdi),%%ymm5		/* t5r */\n\t"\
			"vsubpd	%%ymm4,%%ymm0,%%ymm0			/* B4r =t2r-t5i */\n\t"\
			"vsubpd	%%ymm5,%%ymm3,%%ymm3			/* B3i =t2i-t5r */\n\t"\
			"vaddpd	%%ymm4,%%ymm4,%%ymm4			/*        2*t5i */\n\t"\
			"vaddpd	%%ymm5,%%ymm5,%%ymm5			/*        2*t5r */\n\t"\
			"vaddpd	%%ymm0,%%ymm4,%%ymm4			/* B3r =t2r+t5i */\n\t"\
			"vaddpd	%%ymm3,%%ymm5,%%ymm5			/* B4i =t2i+t5r */\n\t"\
			"vmovaps	%%ymm0,0x400(%%rdi)		/* <-B4r = s1p04r */\n\t"\
			"vmovaps	%%ymm3,0x320(%%rdi)		/* <-B3i = s1p24r */\n\t"\
			"vmovaps	%%ymm4,0x300(%%rdi)		/* <-B3r = s1p24r */\n\t"\
			"vmovaps	%%ymm5,0x420(%%rdi)		/* <-B4i = s1p04i */\n\t"\
			"\n\t"\
		"/* Since doing radix-7 in-place here, inputs of radix-4 are in consecutive memory locs, i.e. 0x40 bytes apart, e.g. the distance between s1p00r and s1p01r: */\n\t"\
			"\n\t"\
		"/* SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_E(s1p00r,s1p01r,s1p02r,s1p03r, add00+p[0,1,2,3]): */\n\t"\
			"movq	%[__out],%%rsi	/* s1p00r */\n\t"\
			"movq	%[__add],%%rax\n\t"\
			"movslq	%[__p01],%%rbx\n\t"\
			"movslq	%[__p02],%%rcx\n\t"\
			"movslq	%[__p03],%%rdx\n\t"\
			"shlq	$3,%%rbx		/* Pointer offset for floating doubles */\n\t"\
			"shlq	$3,%%rcx\n\t"\
			"shlq	$3,%%rdx\n\t"\
			"addq	%%rax,%%rbx\n\t"\
			"addq	%%rax,%%rcx\n\t"\
			"addq	%%rax,%%rdx\n\t"\
			"\n\t"\
			"vmovaps	     (%%rsi),%%ymm0		/* in0r = s1p00r */\n\t"\
			"vmovaps	0x040(%%rsi),%%ymm4		/* in1r = s1p01r */\n\t"\
			"vmovaps	0x020(%%rsi),%%ymm1		/* in0i = s1p00i */\n\t"\
			"vmovaps	0x060(%%rsi),%%ymm5		/* in1i = s1p01i */\n\t"\
			"vmovaps	%%ymm0,%%ymm2			/* ymm2 <- in0r */\n\t"\
			"vmovaps	%%ymm4,%%ymm6			/* ymm4 <- in1r */\n\t"\
			"vmovaps	%%ymm1,%%ymm3			/* ymm3 <- in0i */\n\t"\
			"vmovaps	%%ymm5,%%ymm7			/* ymm5 <- in1i */\n\t"\
			"\n\t"\
			"vaddpd	0x080(%%rsi),%%ymm0,%%ymm0		/* t1 */\n\t"\
			"vaddpd	0x0c0(%%rsi),%%ymm4,%%ymm4		/* t5 */\n\t"\
			"vaddpd	0x0a0(%%rsi),%%ymm1,%%ymm1		/* t2 */\n\t"\
			"vaddpd	0x0e0(%%rsi),%%ymm5,%%ymm5		/* t6 */\n\t"\
			"vsubpd	0x080(%%rsi),%%ymm2,%%ymm2		/* t3 */\n\t"\
			"vsubpd	0x0c0(%%rsi),%%ymm6,%%ymm6		/* t7 */\n\t"\
			"vsubpd	0x0a0(%%rsi),%%ymm3,%%ymm3		/* t4 */\n\t"\
			"vsubpd	0x0e0(%%rsi),%%ymm7,%%ymm7		/* t8 */\n\t"\
			"/* Finish radix-4 butterfly and store results into output-array slots: */\n\t"\
			"vsubpd	%%ymm4,%%ymm0,%%ymm0			/* ~t5 <- t1 -t5 */\n\t"\
			"vsubpd	%%ymm7,%%ymm2,%%ymm2			/* ~t7 <- t3 -t8 */\n\t"\
			"vsubpd	%%ymm5,%%ymm1,%%ymm1			/* ~t6 <- t2 -t6 */\n\t"\
			"vsubpd	%%ymm6,%%ymm3,%%ymm3			/* ~t4 <- t4 -t7 */\n\t"\
			"vmovaps	%%ymm0,      (%%rbx)		/* <- ~t5 */\n\t"\
			"vmovaps	%%ymm2,      (%%rcx)		/* <- ~t7 */\n\t"\
			"vmovaps	%%ymm1,0x020(%%rbx)		/* <- ~t6 */\n\t"\
			"vmovaps	%%ymm3,0x020(%%rdx)		/* <- ~t4 */\n\t"\
			"vaddpd	%%ymm4,%%ymm4,%%ymm4			/*          2*t5 */\n\t"\
			"vaddpd	%%ymm7,%%ymm7,%%ymm7			/*          2*t8 */\n\t"\
			"vaddpd	%%ymm5,%%ymm5,%%ymm5			/*          2*t6 */\n\t"\
			"vaddpd	%%ymm6,%%ymm6,%%ymm6			/*          2*t7 */\n\t"\
			"vaddpd	%%ymm0,%%ymm4,%%ymm4			/* ~t1 <- t1 +t5 */\n\t"\
			"vaddpd	%%ymm2,%%ymm7,%%ymm7			/* ~t3 <- t3 +t8 */\n\t"\
			"vaddpd	%%ymm1,%%ymm5,%%ymm5			/* ~t2 <- t2 +t6 */\n\t"\
			"vaddpd	%%ymm3,%%ymm6,%%ymm6			/* ~t8 <- t4 +t7 */\n\t"\
			"vmovaps	%%ymm4,      (%%rax)		/* <- ~t1 */\n\t"\
			"vmovaps	%%ymm7,      (%%rdx)		/* <- ~t3 */\n\t"\
			"vmovaps	%%ymm5,0x020(%%rax)		/* <- ~t2 */\n\t"\
			"vmovaps	%%ymm6,0x020(%%rcx)		/* <- ~t8 */\n\t"\
		"/* SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_E(s1p04r,s1p05r,s1p06r,s1p07r, add24+p[1,0,3,2]): */\n\t"\
			"addq	$0x100,%%rsi		/* s1p04r */\n\t"\
			"movslq	%[__p24],%%rdi\n\t"\
			"shlq	$3,%%rdi		/* Pointer offset for floating doubles */\n\t"\
			"subq	%%rax,%%rbx\n\t"\
			"subq	%%rax,%%rcx\n\t"\
			"subq	%%rax,%%rdx\n\t"\
			"addq	%%rdi,%%rax		/* add24 */\n\t"\
			"addq	%%rax,%%rbx\n\t"\
			"addq	%%rax,%%rcx\n\t"\
			"addq	%%rax,%%rdx\n\t"\
			"\n\t"\
			"vmovaps	     (%%rsi),%%ymm0		/* in0r = s1p00r */\n\t"\
			"vmovaps	0x040(%%rsi),%%ymm4		/* in1r = s1p01r */\n\t"\
			"vmovaps	0x020(%%rsi),%%ymm1		/* in0i = s1p00i */\n\t"\
			"vmovaps	0x060(%%rsi),%%ymm5		/* in1i = s1p01i */\n\t"\
			"vmovaps	%%ymm0,%%ymm2			/* ymm2 <- in0r */\n\t"\
			"vmovaps	%%ymm4,%%ymm6			/* ymm4 <- in1r */\n\t"\
			"vmovaps	%%ymm1,%%ymm3			/* ymm3 <- in0i */\n\t"\
			"vmovaps	%%ymm5,%%ymm7			/* ymm5 <- in1i */\n\t"\
			"\n\t"\
			"vaddpd	0x080(%%rsi),%%ymm0,%%ymm0		/* t1 */\n\t"\
			"vaddpd	0x0c0(%%rsi),%%ymm4,%%ymm4		/* t5 */\n\t"\
			"vaddpd	0x0a0(%%rsi),%%ymm1,%%ymm1		/* t2 */\n\t"\
			"vaddpd	0x0e0(%%rsi),%%ymm5,%%ymm5		/* t6 */\n\t"\
			"vsubpd	0x080(%%rsi),%%ymm2,%%ymm2		/* t3 */\n\t"\
			"vsubpd	0x0c0(%%rsi),%%ymm6,%%ymm6		/* t7 */\n\t"\
			"vsubpd	0x0a0(%%rsi),%%ymm3,%%ymm3		/* t4 */\n\t"\
			"vsubpd	0x0e0(%%rsi),%%ymm7,%%ymm7		/* t8 */\n\t"\
			"/* Finish radix-4 butterfly and store results into output-array slots: */\n\t"\
			"/* Swap _o1 <-> _o2, _o3 <-> _o4: eax <-> ebx, ecx <-> edx */\n\t"\
			"vsubpd	%%ymm4,%%ymm0,%%ymm0			/* ~t5 <- t1 -t5 */\n\t"\
			"vsubpd	%%ymm7,%%ymm2,%%ymm2			/* ~t7 <- t3 -t8 */\n\t"\
			"vsubpd	%%ymm5,%%ymm1,%%ymm1			/* ~t6 <- t2 -t6 */\n\t"\
			"vsubpd	%%ymm6,%%ymm3,%%ymm3			/* ~t4 <- t4 -t7 */\n\t"\
			"vmovaps	%%ymm0,      (%%rax)		/* <- ~t5 */\n\t"\
			"vmovaps	%%ymm2,      (%%rdx)		/* <- ~t7 */\n\t"\
			"vmovaps	%%ymm1,0x020(%%rax)		/* <- ~t6 */\n\t"\
			"vmovaps	%%ymm3,0x020(%%rcx)		/* <- ~t4 */\n\t"\
			"vaddpd	%%ymm4,%%ymm4,%%ymm4			/*          2*t5 */\n\t"\
			"vaddpd	%%ymm7,%%ymm7,%%ymm7			/*          2*t8 */\n\t"\
			"vaddpd	%%ymm5,%%ymm5,%%ymm5			/*          2*t6 */\n\t"\
			"vaddpd	%%ymm6,%%ymm6,%%ymm6			/*          2*t7 */\n\t"\
			"vaddpd	%%ymm0,%%ymm4,%%ymm4			/* ~t1 <- t1 +t5 */\n\t"\
			"vaddpd	%%ymm2,%%ymm7,%%ymm7			/* ~t3 <- t3 +t8 */\n\t"\
			"vaddpd	%%ymm1,%%ymm5,%%ymm5			/* ~t2 <- t2 +t6 */\n\t"\
			"vaddpd	%%ymm3,%%ymm6,%%ymm6			/* ~t8 <- t4 +t7 */\n\t"\
			"vmovaps	%%ymm4,      (%%rbx)		/* <- ~t1 */\n\t"\
			"vmovaps	%%ymm7,      (%%rcx)		/* <- ~t3 */\n\t"\
			"vmovaps	%%ymm5,0x020(%%rbx)		/* <- ~t2 */\n\t"\
			"vmovaps	%%ymm6,0x020(%%rdx)		/* <- ~t8 */\n\t"\
		"/* SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_E(s1p08r,s1p09r,s1p10r,s1p11r, add20+p[2,3,1,0]): */\n\t"\
			"addq	$0x100,%%rsi		/* s1p08r */\n\t"\
			"subq	%%rax,%%rbx\n\t"\
			"subq	%%rax,%%rcx\n\t"\
			"subq	%%rax,%%rdx\n\t"\
			"subq	%%rdi,%%rax		/* add0  */\n\t"\
			"movslq	%[__p20],%%rdi\n\t"\
			"shlq	$3,%%rdi		/* Pointer offset for floating doubles */\n\t"\
			"addq	%%rdi,%%rax		/* add20 */\n\t"\
			"addq	%%rax,%%rbx\n\t"\
			"addq	%%rax,%%rcx\n\t"\
			"addq	%%rax,%%rdx\n\t"\
			"\n\t"\
			"vmovaps	     (%%rsi),%%ymm0		/* in0r = s1p00r */\n\t"\
			"vmovaps	0x040(%%rsi),%%ymm4		/* in1r = s1p01r */\n\t"\
			"vmovaps	0x020(%%rsi),%%ymm1		/* in0i = s1p00i */\n\t"\
			"vmovaps	0x060(%%rsi),%%ymm5		/* in1i = s1p01i */\n\t"\
			"vmovaps	%%ymm0,%%ymm2			/* ymm2 <- in0r */\n\t"\
			"vmovaps	%%ymm4,%%ymm6			/* ymm4 <- in1r */\n\t"\
			"vmovaps	%%ymm1,%%ymm3			/* ymm3 <- in0i */\n\t"\
			"vmovaps	%%ymm5,%%ymm7			/* ymm5 <- in1i */\n\t"\
			"\n\t"\
			"vaddpd	0x080(%%rsi),%%ymm0,%%ymm0		/* t1 */\n\t"\
			"vaddpd	0x0c0(%%rsi),%%ymm4,%%ymm4		/* t5 */\n\t"\
			"vaddpd	0x0a0(%%rsi),%%ymm1,%%ymm1		/* t2 */\n\t"\
			"vaddpd	0x0e0(%%rsi),%%ymm5,%%ymm5		/* t6 */\n\t"\
			"vsubpd	0x080(%%rsi),%%ymm2,%%ymm2		/* t3 */\n\t"\
			"vsubpd	0x0c0(%%rsi),%%ymm6,%%ymm6		/* t7 */\n\t"\
			"vsubpd	0x0a0(%%rsi),%%ymm3,%%ymm3		/* t4 */\n\t"\
			"vsubpd	0x0e0(%%rsi),%%ymm7,%%ymm7		/* t8 */\n\t"\
			"/* Finish radix-4 butterfly and store results into output-array slots: */\n\t"\
			"/* Swap _o0123 <-> _o2310: e[abcd]x <-> e[cdba]x */\n\t"\
			"vsubpd	%%ymm4,%%ymm0,%%ymm0			/* ~t5 <- t1 -t5 */\n\t"\
			"vsubpd	%%ymm7,%%ymm2,%%ymm2			/* ~t7 <- t3 -t8 */\n\t"\
			"vsubpd	%%ymm5,%%ymm1,%%ymm1			/* ~t6 <- t2 -t6 */\n\t"\
			"vsubpd	%%ymm6,%%ymm3,%%ymm3			/* ~t4 <- t4 -t7 */\n\t"\
			"vmovaps	%%ymm0,      (%%rdx)		/* <- ~t5 */\n\t"\
			"vmovaps	%%ymm2,      (%%rbx)		/* <- ~t7 */\n\t"\
			"vmovaps	%%ymm1,0x020(%%rdx)		/* <- ~t6 */\n\t"\
			"vmovaps	%%ymm3,0x020(%%rax)		/* <- ~t4 */\n\t"\
			"vaddpd	%%ymm4,%%ymm4,%%ymm4			/*          2*t5 */\n\t"\
			"vaddpd	%%ymm7,%%ymm7,%%ymm7			/*          2*t8 */\n\t"\
			"vaddpd	%%ymm5,%%ymm5,%%ymm5			/*          2*t6 */\n\t"\
			"vaddpd	%%ymm6,%%ymm6,%%ymm6			/*          2*t7 */\n\t"\
			"vaddpd	%%ymm0,%%ymm4,%%ymm4			/* ~t1 <- t1 +t5 */\n\t"\
			"vaddpd	%%ymm2,%%ymm7,%%ymm7			/* ~t3 <- t3 +t8 */\n\t"\
			"vaddpd	%%ymm1,%%ymm5,%%ymm5			/* ~t2 <- t2 +t6 */\n\t"\
			"vaddpd	%%ymm3,%%ymm6,%%ymm6			/* ~t8 <- t4 +t7 */\n\t"\
			"vmovaps	%%ymm4,      (%%rcx)		/* <- ~t1 */\n\t"\
			"vmovaps	%%ymm7,      (%%rax)		/* <- ~t3 */\n\t"\
			"vmovaps	%%ymm5,0x020(%%rcx)		/* <- ~t2 */\n\t"\
			"vmovaps	%%ymm6,0x020(%%rbx)		/* <- ~t8 */\n\t"\
		"/* SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_E(s1p12r,s1p13r,s1p14r,s1p15r, add16+p[0,1,2,3]): */\n\t"\
			"addq	$0x100,%%rsi		/* s1p12r */\n\t"\
			"subq	%%rax,%%rbx\n\t"\
			"subq	%%rax,%%rcx\n\t"\
			"subq	%%rax,%%rdx\n\t"\
			"subq	%%rdi,%%rax		/* add0  */\n\t"\
			"movslq	%[__p16],%%rdi\n\t"\
			"shlq	$3,%%rdi		/* Pointer offset for floating doubles */\n\t"\
			"addq	%%rdi,%%rax		/* add16 */\n\t"\
			"addq	%%rax,%%rbx\n\t"\
			"addq	%%rax,%%rcx\n\t"\
			"addq	%%rax,%%rdx\n\t"\
			"\n\t"\
			"vmovaps	     (%%rsi),%%ymm0		/* in0r = s1p00r */\n\t"\
			"vmovaps	0x040(%%rsi),%%ymm4		/* in1r = s1p01r */\n\t"\
			"vmovaps	0x020(%%rsi),%%ymm1		/* in0i = s1p00i */\n\t"\
			"vmovaps	0x060(%%rsi),%%ymm5		/* in1i = s1p01i */\n\t"\
			"vmovaps	%%ymm0,%%ymm2			/* ymm2 <- in0r */\n\t"\
			"vmovaps	%%ymm4,%%ymm6			/* ymm4 <- in1r */\n\t"\
			"vmovaps	%%ymm1,%%ymm3			/* ymm3 <- in0i */\n\t"\
			"vmovaps	%%ymm5,%%ymm7			/* ymm5 <- in1i */\n\t"\
			"\n\t"\
			"vaddpd	0x080(%%rsi),%%ymm0,%%ymm0		/* t1 */\n\t"\
			"vaddpd	0x0c0(%%rsi),%%ymm4,%%ymm4		/* t5 */\n\t"\
			"vaddpd	0x0a0(%%rsi),%%ymm1,%%ymm1		/* t2 */\n\t"\
			"vaddpd	0x0e0(%%rsi),%%ymm5,%%ymm5		/* t6 */\n\t"\
			"vsubpd	0x080(%%rsi),%%ymm2,%%ymm2		/* t3 */\n\t"\
			"vsubpd	0x0c0(%%rsi),%%ymm6,%%ymm6		/* t7 */\n\t"\
			"vsubpd	0x0a0(%%rsi),%%ymm3,%%ymm3		/* t4 */\n\t"\
			"vsubpd	0x0e0(%%rsi),%%ymm7,%%ymm7		/* t8 */\n\t"\
			"/* Finish radix-4 butterfly and store results into output-array slots: */\n\t"\
			"vsubpd	%%ymm4,%%ymm0,%%ymm0			/* ~t5 <- t1 -t5 */\n\t"\
			"vsubpd	%%ymm7,%%ymm2,%%ymm2			/* ~t7 <- t3 -t8 */\n\t"\
			"vsubpd	%%ymm5,%%ymm1,%%ymm1			/* ~t6 <- t2 -t6 */\n\t"\
			"vsubpd	%%ymm6,%%ymm3,%%ymm3			/* ~t4 <- t4 -t7 */\n\t"\
			"vmovaps	%%ymm0,      (%%rbx)		/* <- ~t5 */\n\t"\
			"vmovaps	%%ymm2,      (%%rcx)		/* <- ~t7 */\n\t"\
			"vmovaps	%%ymm1,0x020(%%rbx)		/* <- ~t6 */\n\t"\
			"vmovaps	%%ymm3,0x020(%%rdx)		/* <- ~t4 */\n\t"\
			"vaddpd	%%ymm4,%%ymm4,%%ymm4			/*          2*t5 */\n\t"\
			"vaddpd	%%ymm7,%%ymm7,%%ymm7			/*          2*t8 */\n\t"\
			"vaddpd	%%ymm5,%%ymm5,%%ymm5			/*          2*t6 */\n\t"\
			"vaddpd	%%ymm6,%%ymm6,%%ymm6			/*          2*t7 */\n\t"\
			"vaddpd	%%ymm0,%%ymm4,%%ymm4			/* ~t1 <- t1 +t5 */\n\t"\
			"vaddpd	%%ymm2,%%ymm7,%%ymm7			/* ~t3 <- t3 +t8 */\n\t"\
			"vaddpd	%%ymm1,%%ymm5,%%ymm5			/* ~t2 <- t2 +t6 */\n\t"\
			"vaddpd	%%ymm3,%%ymm6,%%ymm6			/* ~t8 <- t4 +t7 */\n\t"\
			"vmovaps	%%ymm4,      (%%rax)		/* <- ~t1 */\n\t"\
			"vmovaps	%%ymm7,      (%%rdx)		/* <- ~t3 */\n\t"\
			"vmovaps	%%ymm5,0x020(%%rax)		/* <- ~t2 */\n\t"\
			"vmovaps	%%ymm6,0x020(%%rcx)		/* <- ~t8 */\n\t"\
		"/* SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_E(s1p16r,s1p17r,s1p18r,s1p19r, add12+p[3,2,0,1]): */\n\t"\
			"addq	$0x100,%%rsi		/* s1p16r */\n\t"\
			"subq	%%rax,%%rbx\n\t"\
			"subq	%%rax,%%rcx\n\t"\
			"subq	%%rax,%%rdx\n\t"\
			"subq	%%rdi,%%rax		/* add0  */\n\t"\
			"movslq	%[__p12],%%rdi\n\t"\
			"shlq	$3,%%rdi		/* Pointer offset for floating doubles */\n\t"\
			"addq	%%rdi,%%rax		/* add12 */\n\t"\
			"addq	%%rax,%%rbx\n\t"\
			"addq	%%rax,%%rcx\n\t"\
			"addq	%%rax,%%rdx\n\t"\
			"\n\t"\
			"vmovaps	     (%%rsi),%%ymm0		/* in0r = s1p00r */\n\t"\
			"vmovaps	0x040(%%rsi),%%ymm4		/* in1r = s1p01r */\n\t"\
			"vmovaps	0x020(%%rsi),%%ymm1		/* in0i = s1p00i */\n\t"\
			"vmovaps	0x060(%%rsi),%%ymm5		/* in1i = s1p01i */\n\t"\
			"vmovaps	%%ymm0,%%ymm2			/* ymm2 <- in0r */\n\t"\
			"vmovaps	%%ymm4,%%ymm6			/* ymm4 <- in1r */\n\t"\
			"vmovaps	%%ymm1,%%ymm3			/* ymm3 <- in0i */\n\t"\
			"vmovaps	%%ymm5,%%ymm7			/* ymm5 <- in1i */\n\t"\
			"\n\t"\
			"vaddpd	0x080(%%rsi),%%ymm0,%%ymm0		/* t1 */\n\t"\
			"vaddpd	0x0c0(%%rsi),%%ymm4,%%ymm4		/* t5 */\n\t"\
			"vaddpd	0x0a0(%%rsi),%%ymm1,%%ymm1		/* t2 */\n\t"\
			"vaddpd	0x0e0(%%rsi),%%ymm5,%%ymm5		/* t6 */\n\t"\
			"vsubpd	0x080(%%rsi),%%ymm2,%%ymm2		/* t3 */\n\t"\
			"vsubpd	0x0c0(%%rsi),%%ymm6,%%ymm6		/* t7 */\n\t"\
			"vsubpd	0x0a0(%%rsi),%%ymm3,%%ymm3		/* t4 */\n\t"\
			"vsubpd	0x0e0(%%rsi),%%ymm7,%%ymm7		/* t8 */\n\t"\
			"/* Finish radix-4 butterfly and store results into output-array slots: */\n\t"\
			"/* Swap _o0123 <-> _o3201: e[abcd]x <-> e[dcab]x */\n\t"\
			"vsubpd	%%ymm4,%%ymm0,%%ymm0			/* ~t5 <- t1 -t5 */\n\t"\
			"vsubpd	%%ymm7,%%ymm2,%%ymm2			/* ~t7 <- t3 -t8 */\n\t"\
			"vsubpd	%%ymm5,%%ymm1,%%ymm1			/* ~t6 <- t2 -t6 */\n\t"\
			"vsubpd	%%ymm6,%%ymm3,%%ymm3			/* ~t4 <- t4 -t7 */\n\t"\
			"vmovaps	%%ymm0,      (%%rcx)		/* <- ~t5 */\n\t"\
			"vmovaps	%%ymm2,      (%%rax)		/* <- ~t7 */\n\t"\
			"vmovaps	%%ymm1,0x020(%%rcx)		/* <- ~t6 */\n\t"\
			"vmovaps	%%ymm3,0x020(%%rbx)		/* <- ~t4 */\n\t"\
			"vaddpd	%%ymm4,%%ymm4,%%ymm4			/*          2*t5 */\n\t"\
			"vaddpd	%%ymm7,%%ymm7,%%ymm7			/*          2*t8 */\n\t"\
			"vaddpd	%%ymm5,%%ymm5,%%ymm5			/*          2*t6 */\n\t"\
			"vaddpd	%%ymm6,%%ymm6,%%ymm6			/*          2*t7 */\n\t"\
			"vaddpd	%%ymm0,%%ymm4,%%ymm4			/* ~t1 <- t1 +t5 */\n\t"\
			"vaddpd	%%ymm2,%%ymm7,%%ymm7			/* ~t3 <- t3 +t8 */\n\t"\
			"vaddpd	%%ymm1,%%ymm5,%%ymm5			/* ~t2 <- t2 +t6 */\n\t"\
			"vaddpd	%%ymm3,%%ymm6,%%ymm6			/* ~t8 <- t4 +t7 */\n\t"\
			"vmovaps	%%ymm4,      (%%rdx)		/* <- ~t1 */\n\t"\
			"vmovaps	%%ymm7,      (%%rbx)		/* <- ~t3 */\n\t"\
			"vmovaps	%%ymm5,0x020(%%rdx)		/* <- ~t2 */\n\t"\
			"vmovaps	%%ymm6,0x020(%%rax)		/* <- ~t8 */\n\t"\
		"/* SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_E(s1p20r,s1p21r,s1p22r,s1p23r, add08+p[1,0,3,2]): */\n\t"\
			"addq	$0x100,%%rsi		/* s1p20r */\n\t"\
			"subq	%%rax,%%rbx\n\t"\
			"subq	%%rax,%%rcx\n\t"\
			"subq	%%rax,%%rdx\n\t"\
			"subq	%%rdi,%%rax		/* add0  */\n\t"\
			"movslq	%[__p08],%%rdi\n\t"\
			"shlq	$3,%%rdi		/* Pointer offset for floating doubles */\n\t"\
			"addq	%%rdi,%%rax		/* add08 */\n\t"\
			"addq	%%rax,%%rbx\n\t"\
			"addq	%%rax,%%rcx\n\t"\
			"addq	%%rax,%%rdx\n\t"\
			"\n\t"\
			"vmovaps	     (%%rsi),%%ymm0		/* in0r = s1p00r */\n\t"\
			"vmovaps	0x040(%%rsi),%%ymm4		/* in1r = s1p01r */\n\t"\
			"vmovaps	0x020(%%rsi),%%ymm1		/* in0i = s1p00i */\n\t"\
			"vmovaps	0x060(%%rsi),%%ymm5		/* in1i = s1p01i */\n\t"\
			"vmovaps	%%ymm0,%%ymm2			/* ymm2 <- in0r */\n\t"\
			"vmovaps	%%ymm4,%%ymm6			/* ymm4 <- in1r */\n\t"\
			"vmovaps	%%ymm1,%%ymm3			/* ymm3 <- in0i */\n\t"\
			"vmovaps	%%ymm5,%%ymm7			/* ymm5 <- in1i */\n\t"\
			"\n\t"\
			"vaddpd	0x080(%%rsi),%%ymm0,%%ymm0		/* t1 */\n\t"\
			"vaddpd	0x0c0(%%rsi),%%ymm4,%%ymm4		/* t5 */\n\t"\
			"vaddpd	0x0a0(%%rsi),%%ymm1,%%ymm1		/* t2 */\n\t"\
			"vaddpd	0x0e0(%%rsi),%%ymm5,%%ymm5		/* t6 */\n\t"\
			"vsubpd	0x080(%%rsi),%%ymm2,%%ymm2		/* t3 */\n\t"\
			"vsubpd	0x0c0(%%rsi),%%ymm6,%%ymm6		/* t7 */\n\t"\
			"vsubpd	0x0a0(%%rsi),%%ymm3,%%ymm3		/* t4 */\n\t"\
			"vsubpd	0x0e0(%%rsi),%%ymm7,%%ymm7		/* t8 */\n\t"\
			"/* Finish radix-4 butterfly and store results into output-array slots: */\n\t"\
			"/* Swap _o0 <-> _o1, _o2 <-> _o3: eax <-> ebx, ecx <-> edx */\n\t"\
			"vsubpd	%%ymm4,%%ymm0,%%ymm0			/* ~t5 <- t1 -t5 */\n\t"\
			"vsubpd	%%ymm7,%%ymm2,%%ymm2			/* ~t7 <- t3 -t8 */\n\t"\
			"vsubpd	%%ymm5,%%ymm1,%%ymm1			/* ~t6 <- t2 -t6 */\n\t"\
			"vsubpd	%%ymm6,%%ymm3,%%ymm3			/* ~t4 <- t4 -t7 */\n\t"\
			"vmovaps	%%ymm0,      (%%rax)		/* <- ~t5 */\n\t"\
			"vmovaps	%%ymm2,      (%%rdx)		/* <- ~t7 */\n\t"\
			"vmovaps	%%ymm1,0x020(%%rax)		/* <- ~t6 */\n\t"\
			"vmovaps	%%ymm3,0x020(%%rcx)		/* <- ~t4 */\n\t"\
			"vaddpd	%%ymm4,%%ymm4,%%ymm4			/*          2*t5 */\n\t"\
			"vaddpd	%%ymm7,%%ymm7,%%ymm7			/*          2*t8 */\n\t"\
			"vaddpd	%%ymm5,%%ymm5,%%ymm5			/*          2*t6 */\n\t"\
			"vaddpd	%%ymm6,%%ymm6,%%ymm6			/*          2*t7 */\n\t"\
			"vaddpd	%%ymm0,%%ymm4,%%ymm4			/* ~t1 <- t1 +t5 */\n\t"\
			"vaddpd	%%ymm2,%%ymm7,%%ymm7			/* ~t3 <- t3 +t8 */\n\t"\
			"vaddpd	%%ymm1,%%ymm5,%%ymm5			/* ~t2 <- t2 +t6 */\n\t"\
			"vaddpd	%%ymm3,%%ymm6,%%ymm6			/* ~t8 <- t4 +t7 */\n\t"\
			"vmovaps	%%ymm4,      (%%rbx)		/* <- ~t1 */\n\t"\
			"vmovaps	%%ymm7,      (%%rcx)		/* <- ~t3 */\n\t"\
			"vmovaps	%%ymm5,0x020(%%rbx)		/* <- ~t2 */\n\t"\
			"vmovaps	%%ymm6,0x020(%%rdx)		/* <- ~t8 */\n\t"\
		"/* SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_E(s1p24r,s1p25r,s1p26r,s1p27r, add04+p[2,3,1,0]): */\n\t"\
			"addq	$0x100,%%rsi		/* s1p24r */\n\t"\
			"subq	%%rax,%%rbx\n\t"\
			"subq	%%rax,%%rcx\n\t"\
			"subq	%%rax,%%rdx\n\t"\
			"subq	%%rdi,%%rax		/* add0  */\n\t"\
			"movslq	%[__p04],%%rdi\n\t"\
			"shlq	$3,%%rdi		/* Pointer offset for floating doubles */\n\t"\
			"addq	%%rdi,%%rax		/* add04 */\n\t"\
			"addq	%%rax,%%rbx\n\t"\
			"addq	%%rax,%%rcx\n\t"\
			"addq	%%rax,%%rdx\n\t"\
			"\n\t"\
			"vmovaps	     (%%rsi),%%ymm0		/* in0r = s1p00r */\n\t"\
			"vmovaps	0x040(%%rsi),%%ymm4		/* in1r = s1p01r */\n\t"\
			"vmovaps	0x020(%%rsi),%%ymm1		/* in0i = s1p00i */\n\t"\
			"vmovaps	0x060(%%rsi),%%ymm5		/* in1i = s1p01i */\n\t"\
			"vmovaps	%%ymm0,%%ymm2			/* ymm2 <- in0r */\n\t"\
			"vmovaps	%%ymm4,%%ymm6			/* ymm4 <- in1r */\n\t"\
			"vmovaps	%%ymm1,%%ymm3			/* ymm3 <- in0i */\n\t"\
			"vmovaps	%%ymm5,%%ymm7			/* ymm5 <- in1i */\n\t"\
			"\n\t"\
			"vaddpd	0x080(%%rsi),%%ymm0,%%ymm0		/* t1 */\n\t"\
			"vaddpd	0x0c0(%%rsi),%%ymm4,%%ymm4		/* t5 */\n\t"\
			"vaddpd	0x0a0(%%rsi),%%ymm1,%%ymm1		/* t2 */\n\t"\
			"vaddpd	0x0e0(%%rsi),%%ymm5,%%ymm5		/* t6 */\n\t"\
			"vsubpd	0x080(%%rsi),%%ymm2,%%ymm2		/* t3 */\n\t"\
			"vsubpd	0x0c0(%%rsi),%%ymm6,%%ymm6		/* t7 */\n\t"\
			"vsubpd	0x0a0(%%rsi),%%ymm3,%%ymm3		/* t4 */\n\t"\
			"vsubpd	0x0e0(%%rsi),%%ymm7,%%ymm7		/* t8 */\n\t"\
			"/* Finish radix-4 butterfly and store results into output-array slots: */\n\t"\
			"/* Swap _o0123 <-> _o2310: e[abcd]x <-> e[cdba]x */\n\t"\
			"vsubpd	%%ymm4,%%ymm0,%%ymm0			/* ~t5 <- t1 -t5 */\n\t"\
			"vsubpd	%%ymm7,%%ymm2,%%ymm2			/* ~t7 <- t3 -t8 */\n\t"\
			"vsubpd	%%ymm5,%%ymm1,%%ymm1			/* ~t6 <- t2 -t6 */\n\t"\
			"vsubpd	%%ymm6,%%ymm3,%%ymm3			/* ~t4 <- t4 -t7 */\n\t"\
			"vmovaps	%%ymm0,      (%%rdx)		/* <- ~t5 */\n\t"\
			"vmovaps	%%ymm2,      (%%rbx)		/* <- ~t7 */\n\t"\
			"vmovaps	%%ymm1,0x020(%%rdx)		/* <- ~t6 */\n\t"\
			"vmovaps	%%ymm3,0x020(%%rax)		/* <- ~t4 */\n\t"\
			"vaddpd	%%ymm4,%%ymm4,%%ymm4			/*          2*t5 */\n\t"\
			"vaddpd	%%ymm7,%%ymm7,%%ymm7			/*          2*t8 */\n\t"\
			"vaddpd	%%ymm5,%%ymm5,%%ymm5			/*          2*t6 */\n\t"\
			"vaddpd	%%ymm6,%%ymm6,%%ymm6			/*          2*t7 */\n\t"\
			"vaddpd	%%ymm0,%%ymm4,%%ymm4			/* ~t1 <- t1 +t5 */\n\t"\
			"vaddpd	%%ymm2,%%ymm7,%%ymm7			/* ~t3 <- t3 +t8 */\n\t"\
			"vaddpd	%%ymm1,%%ymm5,%%ymm5			/* ~t2 <- t2 +t6 */\n\t"\
			"vaddpd	%%ymm3,%%ymm6,%%ymm6			/* ~t8 <- t4 +t7 */\n\t"\
			"vmovaps	%%ymm4,      (%%rcx)		/* <- ~t1 */\n\t"\
			"vmovaps	%%ymm7,      (%%rax)		/* <- ~t3 */\n\t"\
			"vmovaps	%%ymm5,0x020(%%rcx)		/* <- ~t2 */\n\t"\
			"vmovaps	%%ymm6,0x020(%%rbx)		/* <- ~t8 */\n\t"\
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
			: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
		);\
	}

  #endif	// USE_64BIT_ASM_STYLE ?

#elif defined(USE_SSE2)

  #define USE_64BIT_ASM_STYLE	1	// My x86 timings indicate fancier versions below using ymm0-15 for the radix-16 DFT is faster.

  #if USE_64BIT_ASM_STYLE	// True: Deeper 64-bit-ified version of the 32-bit version of the ASM macros, using all of ymm0-15

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

  #else // USE_64BIT_ASM_STYLE = False: Use simple 64-bit-ified version of the analogous 32-bit ASM macros, i.e. using just ymm0-7.

	#define	SSE2_RADIX28_DIT_NOTWIDDLE(Xadd,Xp01,Xp02,Xp03,Xp04,Xp08,Xp12,Xp16,Xp20,Xp24,Xout,Xcc0)\
	{\
	__asm__ volatile (\
		"/* SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_E(add00+p[0,1,3,2], s1p00r,s1p03r,s1p02r,s1p01r): */\n\t"\
			"movq	%[__out],%%rsi	/* s1p00r */\n\t"\
			"movq	%[__add],%%rax\n\t"\
			"movslq	%[__p01],%%rbx\n\t"\
			"movslq	%[__p02],%%rcx\n\t"\
			"movslq	%[__p03],%%rdx\n\t"\
			"shlq	$3,%%rbx		/* Pointer offset for floating doubles */\n\t"\
			"shlq	$3,%%rcx\n\t"\
			"shlq	$3,%%rdx\n\t"\
			"addq	%%rax,%%rbx\n\t"\
			"addq	%%rax,%%rcx\n\t"\
			"addq	%%rax,%%rdx\n\t"\
			"/* ecx <-> edx */\n\t"\
			"movaps	    (%%rax),%%xmm0		/* a(jt   ) */\n\t"\
			"movaps	    (%%rdx),%%xmm4		/* a(jt+p2) */\n\t"\
			"movaps	0x10(%%rax),%%xmm1		/* a(jp   ) */\n\t"\
			"movaps	0x10(%%rdx),%%xmm5		/* a(jp+p2) */\n\t"\
			"movaps	%%xmm0,%%xmm2			/* xmm2 <- cpy a(jt   ) */\n\t"\
			"movaps	%%xmm4,%%xmm6			/* xmm4 <- cpy a(jt+p2) */\n\t"\
			"movaps	%%xmm1,%%xmm3			/* xmm3 <- cpy a(jp   ) */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* xmm5 <- cpy a(jp+p2) */\n\t"\
			"addpd	    (%%rbx),%%xmm0		/* t1 */\n\t"\
			"addpd	    (%%rcx),%%xmm4		/* t5 */\n\t"\
			"addpd	0x10(%%rbx),%%xmm1		/* t2 */\n\t"\
			"addpd	0x10(%%rcx),%%xmm5		/* t6 */\n\t"\
			"subpd	    (%%rbx),%%xmm2		/* t3 */\n\t"\
			"subpd	    (%%rcx),%%xmm6		/* t7 */\n\t"\
			"subpd	0x10(%%rbx),%%xmm3		/* t4 */\n\t"\
			"subpd	0x10(%%rcx),%%xmm7		/* t8 */\n\t"\
			"/* Finish radix-4 butterfly and store results into output-array slots: */\n\t"\
			"/* Swap _o1 <-> _o3: 0x20,0x30 <-> 0x60,0x70: */\n\t"\
			"subpd	%%xmm4,%%xmm0			/* ~t5 <- t1 -t5 */\n\t"\
			"subpd	%%xmm7,%%xmm2			/* ~t7 <- t3 -t8 */\n\t"\
			"subpd	%%xmm5,%%xmm1			/* ~t6 <- t2 -t6 */\n\t"\
			"subpd	%%xmm6,%%xmm3			/* ~t4 <- t4 -t7 */\n\t"\
			"movaps	%%xmm0,0x040(%%rsi)		/* <- ~t5 */\n\t"\
			"movaps	%%xmm2,0x020(%%rsi)		/* <- ~t7 */\n\t"\
			"movaps	%%xmm1,0x050(%%rsi)		/* <- ~t6 */\n\t"\
			"movaps	%%xmm3,0x070(%%rsi)		/* <- ~t4 */\n\t"\
			"addpd	%%xmm4,%%xmm4			/*          2*t5 */\n\t"\
			"addpd	%%xmm7,%%xmm7			/*          2*t8 */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*          2*t6 */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*          2*t7 */\n\t"\
			"addpd	%%xmm0,%%xmm4			/* ~t1 <- t1 +t5 */\n\t"\
			"addpd	%%xmm2,%%xmm7			/* ~t3 <- t3 +t8 */\n\t"\
			"addpd	%%xmm1,%%xmm5			/* ~t2 <- t2 +t6 */\n\t"\
			"addpd	%%xmm3,%%xmm6			/* ~t8 <- t4 +t7 */\n\t"\
			"movaps	%%xmm4,     (%%rsi)		/* <- ~t1 */\n\t"\
			"movaps	%%xmm7,0x060(%%rsi)		/* <- ~t3 */\n\t"\
			"movaps	%%xmm5,0x010(%%rsi)		/* <- ~t2 */\n\t"\
			"movaps	%%xmm6,0x030(%%rsi)		/* <- ~t8 */\n\t"\
		"/* SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_E(add12+p[3,2,1,0], s1p04r,s1p07r,s1p06r,s1p05r): */\n\t"\
			"addq	$0x80,%%rsi	/* s1p04r */\n\t"\
			"subq	%%rax,%%rbx		/* p01 */\n\t"\
			"subq	%%rax,%%rcx		/* p02 */\n\t"\
			"subq	%%rax,%%rdx		/* p03 */\n\t"\
			"movslq	%[__p12],%%rdi\n\t"\
			"shlq	$3,%%rdi		/* Pointer offset for floating doubles */\n\t"\
			"addq	%%rdi,%%rax		/* add12 */\n\t"\
			"addq	%%rax,%%rbx\n\t"\
			"addq	%%rax,%%rcx\n\t"\
			"addq	%%rax,%%rdx\n\t"\
			"/* eax <-> edx, ebx <-> ecx */\n\t"\
			"movaps	    (%%rdx),%%xmm0		/* a(jt   ) */\n\t"\
			"movaps	    (%%rbx),%%xmm4		/* a(jt+p2) */\n\t"\
			"movaps	0x10(%%rdx),%%xmm1		/* a(jp   ) */\n\t"\
			"movaps	0x10(%%rbx),%%xmm5		/* a(jp+p2) */\n\t"\
			"movaps	%%xmm0,%%xmm2			/* xmm2 <- cpy a(jt   ) */\n\t"\
			"movaps	%%xmm4,%%xmm6			/* xmm4 <- cpy a(jt+p2) */\n\t"\
			"movaps	%%xmm1,%%xmm3			/* xmm3 <- cpy a(jp   ) */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* xmm5 <- cpy a(jp+p2) */\n\t"\
			"addpd	    (%%rcx),%%xmm0		/* t1 */\n\t"\
			"addpd	    (%%rax),%%xmm4		/* t5 */\n\t"\
			"addpd	0x10(%%rcx),%%xmm1		/* t2 */\n\t"\
			"addpd	0x10(%%rax),%%xmm5		/* t6 */\n\t"\
			"subpd	    (%%rcx),%%xmm2		/* t3 */\n\t"\
			"subpd	    (%%rax),%%xmm6		/* t7 */\n\t"\
			"subpd	0x10(%%rcx),%%xmm3		/* t4 */\n\t"\
			"subpd	0x10(%%rax),%%xmm7		/* t8 */\n\t"\
			"/* Finish radix-4 butterfly and store results into output-array slots: */\n\t"\
			"/* Swap _o1 <-> _o3: 0x20,0x30 <-> 0x60,0x70: */\n\t"\
			"subpd	%%xmm4,%%xmm0			/* ~t5 <- t1 -t5 */\n\t"\
			"subpd	%%xmm7,%%xmm2			/* ~t7 <- t3 -t8 */\n\t"\
			"subpd	%%xmm5,%%xmm1			/* ~t6 <- t2 -t6 */\n\t"\
			"subpd	%%xmm6,%%xmm3			/* ~t4 <- t4 -t7 */\n\t"\
			"movaps	%%xmm0,0x040(%%rsi)		/* <- ~t5 */\n\t"\
			"movaps	%%xmm2,0x020(%%rsi)		/* <- ~t7 */\n\t"\
			"movaps	%%xmm1,0x050(%%rsi)		/* <- ~t6 */\n\t"\
			"movaps	%%xmm3,0x070(%%rsi)		/* <- ~t4 */\n\t"\
			"addpd	%%xmm4,%%xmm4			/*          2*t5 */\n\t"\
			"addpd	%%xmm7,%%xmm7			/*          2*t8 */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*          2*t6 */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*          2*t7 */\n\t"\
			"addpd	%%xmm0,%%xmm4			/* ~t1 <- t1 +t5 */\n\t"\
			"addpd	%%xmm2,%%xmm7			/* ~t3 <- t3 +t8 */\n\t"\
			"addpd	%%xmm1,%%xmm5			/* ~t2 <- t2 +t6 */\n\t"\
			"addpd	%%xmm3,%%xmm6			/* ~t8 <- t4 +t7 */\n\t"\
			"movaps	%%xmm4,     (%%rsi)		/* <- ~t1 */\n\t"\
			"movaps	%%xmm7,0x060(%%rsi)		/* <- ~t3 */\n\t"\
			"movaps	%%xmm5,0x010(%%rsi)		/* <- ~t2 */\n\t"\
			"movaps	%%xmm6,0x030(%%rsi)		/* <- ~t8 */\n\t"\
		"/* SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_E(add24+p[1,0,2,3], s1p08r,s1p11r,s1p10r,s1p09r): */\n\t"\
			"addq	$0x80,%%rsi	/* s1p08r */\n\t"\
			"subq	%%rax,%%rbx		/* p01 */\n\t"\
			"subq	%%rax,%%rcx		/* p02 */\n\t"\
			"subq	%%rax,%%rdx		/* p03 */\n\t"\
			"subq	%%rdi,%%rax		/* add0 */\n\t"\
			"movslq	%[__p24],%%rdi\n\t"\
			"shlq	$3,%%rdi		/* Pointer offset for floating doubles */\n\t"\
			"addq	%%rdi,%%rax		/* add24 */\n\t"\
			"addq	%%rax,%%rbx\n\t"\
			"addq	%%rax,%%rcx\n\t"\
			"addq	%%rax,%%rdx\n\t"\
			"/* eax <-> ebx */\n\t"\
			"movaps	    (%%rbx),%%xmm0		/* a(jt   ) */\n\t"\
			"movaps	    (%%rcx),%%xmm4		/* a(jt+p2) */\n\t"\
			"movaps	0x10(%%rbx),%%xmm1		/* a(jp   ) */\n\t"\
			"movaps	0x10(%%rcx),%%xmm5		/* a(jp+p2) */\n\t"\
			"movaps	%%xmm0,%%xmm2			/* xmm2 <- cpy a(jt   ) */\n\t"\
			"movaps	%%xmm4,%%xmm6			/* xmm4 <- cpy a(jt+p2) */\n\t"\
			"movaps	%%xmm1,%%xmm3			/* xmm3 <- cpy a(jp   ) */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* xmm5 <- cpy a(jp+p2) */\n\t"\
			"addpd	    (%%rax),%%xmm0		/* t1 */\n\t"\
			"addpd	    (%%rdx),%%xmm4		/* t5 */\n\t"\
			"addpd	0x10(%%rax),%%xmm1		/* t2 */\n\t"\
			"addpd	0x10(%%rdx),%%xmm5		/* t6 */\n\t"\
			"subpd	    (%%rax),%%xmm2		/* t3 */\n\t"\
			"subpd	    (%%rdx),%%xmm6		/* t7 */\n\t"\
			"subpd	0x10(%%rax),%%xmm3		/* t4 */\n\t"\
			"subpd	0x10(%%rdx),%%xmm7		/* t8 */\n\t"\
			"/* Finish radix-4 butterfly and store results into output-array slots: */\n\t"\
			"/* Swap _o1 <-> _o3: 0x20,0x30 <-> 0x60,0x70: */\n\t"\
			"subpd	%%xmm4,%%xmm0			/* ~t5 <- t1 -t5 */\n\t"\
			"subpd	%%xmm7,%%xmm2			/* ~t7 <- t3 -t8 */\n\t"\
			"subpd	%%xmm5,%%xmm1			/* ~t6 <- t2 -t6 */\n\t"\
			"subpd	%%xmm6,%%xmm3			/* ~t4 <- t4 -t7 */\n\t"\
			"movaps	%%xmm0,0x040(%%rsi)		/* <- ~t5 */\n\t"\
			"movaps	%%xmm2,0x020(%%rsi)		/* <- ~t7 */\n\t"\
			"movaps	%%xmm1,0x050(%%rsi)		/* <- ~t6 */\n\t"\
			"movaps	%%xmm3,0x070(%%rsi)		/* <- ~t4 */\n\t"\
			"addpd	%%xmm4,%%xmm4			/*          2*t5 */\n\t"\
			"addpd	%%xmm7,%%xmm7			/*          2*t8 */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*          2*t6 */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*          2*t7 */\n\t"\
			"addpd	%%xmm0,%%xmm4			/* ~t1 <- t1 +t5 */\n\t"\
			"addpd	%%xmm2,%%xmm7			/* ~t3 <- t3 +t8 */\n\t"\
			"addpd	%%xmm1,%%xmm5			/* ~t2 <- t2 +t6 */\n\t"\
			"addpd	%%xmm3,%%xmm6			/* ~t8 <- t4 +t7 */\n\t"\
			"movaps	%%xmm4,     (%%rsi)		/* <- ~t1 */\n\t"\
			"movaps	%%xmm7,0x060(%%rsi)		/* <- ~t3 */\n\t"\
			"movaps	%%xmm5,0x010(%%rsi)		/* <- ~t2 */\n\t"\
			"movaps	%%xmm6,0x030(%%rsi)		/* <- ~t8 */\n\t"\
		"/* SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_E(add08+p[1,0,2,3], s1p12r,s1p15r,s1p14r,s1p13r): */\n\t"\
			"addq	$0x80,%%rsi	/* s1p12r */\n\t"\
			"subq	%%rax,%%rbx		/* p01 */\n\t"\
			"subq	%%rax,%%rcx		/* p02 */\n\t"\
			"subq	%%rax,%%rdx		/* p03 */\n\t"\
			"subq	%%rdi,%%rax		/* add0 */\n\t"\
			"movslq	%[__p08],%%rdi\n\t"\
			"shlq	$3,%%rdi		/* Pointer offset for floating doubles */\n\t"\
			"addq	%%rdi,%%rax		/* add08 */\n\t"\
			"addq	%%rax,%%rbx\n\t"\
			"addq	%%rax,%%rcx\n\t"\
			"addq	%%rax,%%rdx\n\t"\
			"/* eax <-> %%rbx */\n\t"\
			"movaps	    (%%rbx),%%xmm0		/* a(jt   ) */\n\t"\
			"movaps	    (%%rcx),%%xmm4		/* a(jt+p2) */\n\t"\
			"movaps	0x10(%%rbx),%%xmm1		/* a(jp   ) */\n\t"\
			"movaps	0x10(%%rcx),%%xmm5		/* a(jp+p2) */\n\t"\
			"movaps	%%xmm0,%%xmm2			/* xmm2 <- cpy a(jt   ) */\n\t"\
			"movaps	%%xmm4,%%xmm6			/* xmm4 <- cpy a(jt+p2) */\n\t"\
			"movaps	%%xmm1,%%xmm3			/* xmm3 <- cpy a(jp   ) */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* xmm5 <- cpy a(jp+p2) */\n\t"\
			"addpd	    (%%rax),%%xmm0		/* t1 */\n\t"\
			"addpd	    (%%rdx),%%xmm4		/* t5 */\n\t"\
			"addpd	0x10(%%rax),%%xmm1		/* t2 */\n\t"\
			"addpd	0x10(%%rdx),%%xmm5		/* t6 */\n\t"\
			"subpd	    (%%rax),%%xmm2		/* t3 */\n\t"\
			"subpd	    (%%rdx),%%xmm6		/* t7 */\n\t"\
			"subpd	0x10(%%rax),%%xmm3		/* t4 */\n\t"\
			"subpd	0x10(%%rdx),%%xmm7		/* t8 */\n\t"\
			"/* Finish radix-4 butterfly and store results into output-array slots: */\n\t"\
			"/* Swap _o1 <-> _o3: 0x20,0x30 <-> 0x60,0x70: */\n\t"\
			"subpd	%%xmm4,%%xmm0			/* ~t5 <- t1 -t5 */\n\t"\
			"subpd	%%xmm7,%%xmm2			/* ~t7 <- t3 -t8 */\n\t"\
			"subpd	%%xmm5,%%xmm1			/* ~t6 <- t2 -t6 */\n\t"\
			"subpd	%%xmm6,%%xmm3			/* ~t4 <- t4 -t7 */\n\t"\
			"movaps	%%xmm0,0x040(%%rsi)		/* <- ~t5 */\n\t"\
			"movaps	%%xmm2,0x020(%%rsi)		/* <- ~t7 */\n\t"\
			"movaps	%%xmm1,0x050(%%rsi)		/* <- ~t6 */\n\t"\
			"movaps	%%xmm3,0x070(%%rsi)		/* <- ~t4 */\n\t"\
			"addpd	%%xmm4,%%xmm4			/*          2*t5 */\n\t"\
			"addpd	%%xmm7,%%xmm7			/*          2*t8 */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*          2*t6 */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*          2*t7 */\n\t"\
			"addpd	%%xmm0,%%xmm4			/* ~t1 <- t1 +t5 */\n\t"\
			"addpd	%%xmm2,%%xmm7			/* ~t3 <- t3 +t8 */\n\t"\
			"addpd	%%xmm1,%%xmm5			/* ~t2 <- t2 +t6 */\n\t"\
			"addpd	%%xmm3,%%xmm6			/* ~t8 <- t4 +t7 */\n\t"\
			"movaps	%%xmm4,     (%%rsi)		/* <- ~t1 */\n\t"\
			"movaps	%%xmm7,0x060(%%rsi)		/* <- ~t3 */\n\t"\
			"movaps	%%xmm5,0x010(%%rsi)		/* <- ~t2 */\n\t"\
			"movaps	%%xmm6,0x030(%%rsi)		/* <- ~t8 */\n\t"\
		"/* SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_E(add20+p[2,3,0,1], s1p16r,s1p19r,s1p18r,s1p17r): */\n\t"\
			"addq	$0x80,%%rsi	/* s1p16r */\n\t"\
			"subq	%%rax,%%rbx		/* p01 */\n\t"\
			"subq	%%rax,%%rcx		/* p02 */\n\t"\
			"subq	%%rax,%%rdx		/* p03 */\n\t"\
			"subq	%%rdi,%%rax		/* add0 */\n\t"\
			"movslq	%[__p20],%%rdi\n\t"\
			"shlq	$3,%%rdi		/* Pointer offset for floating doubles */\n\t"\
			"addq	%%rdi,%%rax		/* add20 */\n\t"\
			"addq	%%rax,%%rbx\n\t"\
			"addq	%%rax,%%rcx\n\t"\
			"addq	%%rax,%%rdx\n\t"\
			"/* eax <-> %%rcx, %%rbx <-> edx */\n\t"\
			"movaps	    (%%rcx),%%xmm0		/* a(jt   ) */\n\t"\
			"movaps	    (%%rax),%%xmm4		/* a(jt+p2) */\n\t"\
			"movaps	0x10(%%rcx),%%xmm1		/* a(jp   ) */\n\t"\
			"movaps	0x10(%%rax),%%xmm5		/* a(jp+p2) */\n\t"\
			"movaps	%%xmm0,%%xmm2			/* xmm2 <- cpy a(jt   ) */\n\t"\
			"movaps	%%xmm4,%%xmm6			/* xmm4 <- cpy a(jt+p2) */\n\t"\
			"movaps	%%xmm1,%%xmm3			/* xmm3 <- cpy a(jp   ) */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* xmm5 <- cpy a(jp+p2) */\n\t"\
			"addpd	    (%%rdx),%%xmm0		/* t1 */\n\t"\
			"addpd	    (%%rbx),%%xmm4		/* t5 */\n\t"\
			"addpd	0x10(%%rdx),%%xmm1		/* t2 */\n\t"\
			"addpd	0x10(%%rbx),%%xmm5		/* t6 */\n\t"\
			"subpd	    (%%rdx),%%xmm2		/* t3 */\n\t"\
			"subpd	    (%%rbx),%%xmm6		/* t7 */\n\t"\
			"subpd	0x10(%%rdx),%%xmm3		/* t4 */\n\t"\
			"subpd	0x10(%%rbx),%%xmm7		/* t8 */\n\t"\
			"/* Finish radix-4 butterfly and store results into output-array slots: */\n\t"\
			"/* Swap _o1 <-> _o3: 0x20,0x30 <-> 0x60,0x70: */\n\t"\
			"subpd	%%xmm4,%%xmm0			/* ~t5 <- t1 -t5 */\n\t"\
			"subpd	%%xmm7,%%xmm2			/* ~t7 <- t3 -t8 */\n\t"\
			"subpd	%%xmm5,%%xmm1			/* ~t6 <- t2 -t6 */\n\t"\
			"subpd	%%xmm6,%%xmm3			/* ~t4 <- t4 -t7 */\n\t"\
			"movaps	%%xmm0,0x040(%%rsi)		/* <- ~t5 */\n\t"\
			"movaps	%%xmm2,0x020(%%rsi)		/* <- ~t7 */\n\t"\
			"movaps	%%xmm1,0x050(%%rsi)		/* <- ~t6 */\n\t"\
			"movaps	%%xmm3,0x070(%%rsi)		/* <- ~t4 */\n\t"\
			"addpd	%%xmm4,%%xmm4			/*          2*t5 */\n\t"\
			"addpd	%%xmm7,%%xmm7			/*          2*t8 */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*          2*t6 */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*          2*t7 */\n\t"\
			"addpd	%%xmm0,%%xmm4			/* ~t1 <- t1 +t5 */\n\t"\
			"addpd	%%xmm2,%%xmm7			/* ~t3 <- t3 +t8 */\n\t"\
			"addpd	%%xmm1,%%xmm5			/* ~t2 <- t2 +t6 */\n\t"\
			"addpd	%%xmm3,%%xmm6			/* ~t8 <- t4 +t7 */\n\t"\
			"movaps	%%xmm4,     (%%rsi)		/* <- ~t1 */\n\t"\
			"movaps	%%xmm7,0x060(%%rsi)		/* <- ~t3 */\n\t"\
			"movaps	%%xmm5,0x010(%%rsi)		/* <- ~t2 */\n\t"\
			"movaps	%%xmm6,0x030(%%rsi)		/* <- ~t8 */\n\t"\
		"/* SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_E(add04+p[2,3,0,1], s1p20r,s1p23r,s1p22r,s1p21r): */\n\t"\
			"addq	$0x80,%%rsi	/* s1p20r */\n\t"\
			"subq	%%rax,%%rbx		/* p01 */\n\t"\
			"subq	%%rax,%%rcx		/* p02 */\n\t"\
			"subq	%%rax,%%rdx		/* p03 */\n\t"\
			"subq	%%rdi,%%rax		/* add0 */\n\t"\
			"movslq	%[__p04],%%rdi\n\t"\
			"shlq	$3,%%rdi		/* Pointer offset for floating doubles */\n\t"\
			"addq	%%rdi,%%rax		/* add04 */\n\t"\
			"addq	%%rax,%%rbx\n\t"\
			"addq	%%rax,%%rcx\n\t"\
			"addq	%%rax,%%rdx\n\t"\
			"/* eax <-> %%rcx, %%rbx <-> edx */\n\t"\
			"movaps	    (%%rcx),%%xmm0		/* a(jt   ) */\n\t"\
			"movaps	    (%%rax),%%xmm4		/* a(jt+p2) */\n\t"\
			"movaps	0x10(%%rcx),%%xmm1		/* a(jp   ) */\n\t"\
			"movaps	0x10(%%rax),%%xmm5		/* a(jp+p2) */\n\t"\
			"movaps	%%xmm0,%%xmm2			/* xmm2 <- cpy a(jt   ) */\n\t"\
			"movaps	%%xmm4,%%xmm6			/* xmm4 <- cpy a(jt+p2) */\n\t"\
			"movaps	%%xmm1,%%xmm3			/* xmm3 <- cpy a(jp   ) */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* xmm5 <- cpy a(jp+p2) */\n\t"\
			"addpd	    (%%rdx),%%xmm0		/* t1 */\n\t"\
			"addpd	    (%%rbx),%%xmm4		/* t5 */\n\t"\
			"addpd	0x10(%%rdx),%%xmm1		/* t2 */\n\t"\
			"addpd	0x10(%%rbx),%%xmm5		/* t6 */\n\t"\
			"subpd	    (%%rdx),%%xmm2		/* t3 */\n\t"\
			"subpd	    (%%rbx),%%xmm6		/* t7 */\n\t"\
			"subpd	0x10(%%rdx),%%xmm3		/* t4 */\n\t"\
			"subpd	0x10(%%rbx),%%xmm7		/* t8 */\n\t"\
			"/* Finish radix-4 butterfly and store results into output-array slots: */\n\t"\
			"/* Swap _o1 <-> _o3: 0x20,0x30 <-> 0x60,0x70: */\n\t"\
			"subpd	%%xmm4,%%xmm0			/* ~t5 <- t1 -t5 */\n\t"\
			"subpd	%%xmm7,%%xmm2			/* ~t7 <- t3 -t8 */\n\t"\
			"subpd	%%xmm5,%%xmm1			/* ~t6 <- t2 -t6 */\n\t"\
			"subpd	%%xmm6,%%xmm3			/* ~t4 <- t4 -t7 */\n\t"\
			"movaps	%%xmm0,0x040(%%rsi)		/* <- ~t5 */\n\t"\
			"movaps	%%xmm2,0x020(%%rsi)		/* <- ~t7 */\n\t"\
			"movaps	%%xmm1,0x050(%%rsi)		/* <- ~t6 */\n\t"\
			"movaps	%%xmm3,0x070(%%rsi)		/* <- ~t4 */\n\t"\
			"addpd	%%xmm4,%%xmm4			/*          2*t5 */\n\t"\
			"addpd	%%xmm7,%%xmm7			/*          2*t8 */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*          2*t6 */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*          2*t7 */\n\t"\
			"addpd	%%xmm0,%%xmm4			/* ~t1 <- t1 +t5 */\n\t"\
			"addpd	%%xmm2,%%xmm7			/* ~t3 <- t3 +t8 */\n\t"\
			"addpd	%%xmm1,%%xmm5			/* ~t2 <- t2 +t6 */\n\t"\
			"addpd	%%xmm3,%%xmm6			/* ~t8 <- t4 +t7 */\n\t"\
			"movaps	%%xmm4,     (%%rsi)		/* <- ~t1 */\n\t"\
			"movaps	%%xmm7,0x060(%%rsi)		/* <- ~t3 */\n\t"\
			"movaps	%%xmm5,0x010(%%rsi)		/* <- ~t2 */\n\t"\
			"movaps	%%xmm6,0x030(%%rsi)		/* <- ~t8 */\n\t"\
		"/* SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_E(add16+p[0,1,3,2], s1p24r,s1p27r,s1p26r,s1p25r): */\n\t"\
			"addq	$0x80,%%rsi	/* s1p24r */\n\t"\
			"subq	%%rax,%%rbx		/* p01 */\n\t"\
			"subq	%%rax,%%rcx		/* p02 */\n\t"\
			"subq	%%rax,%%rdx		/* p03 */\n\t"\
			"subq	%%rdi,%%rax		/* add0 */\n\t"\
			"movslq	%[__p16],%%rdi\n\t"\
			"shlq	$3,%%rdi		/* Pointer offset for floating doubles */\n\t"\
			"addq	%%rdi,%%rax		/* add16 */\n\t"\
			"addq	%%rax,%%rbx\n\t"\
			"addq	%%rax,%%rcx\n\t"\
			"addq	%%rax,%%rdx\n\t"\
			"/* ecx <-> %%rdx */\n\t"\
			"movaps	    (%%rax),%%xmm0		/* a(jt   ) */\n\t"\
			"movaps	    (%%rdx),%%xmm4		/* a(jt+p2) */\n\t"\
			"movaps	0x10(%%rax),%%xmm1		/* a(jp   ) */\n\t"\
			"movaps	0x10(%%rdx),%%xmm5		/* a(jp+p2) */\n\t"\
			"movaps	%%xmm0,%%xmm2			/* xmm2 <- cpy a(jt   ) */\n\t"\
			"movaps	%%xmm4,%%xmm6			/* xmm4 <- cpy a(jt+p2) */\n\t"\
			"movaps	%%xmm1,%%xmm3			/* xmm3 <- cpy a(jp   ) */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* xmm5 <- cpy a(jp+p2) */\n\t"\
			"addpd	    (%%rbx),%%xmm0		/* t1 */\n\t"\
			"addpd	    (%%rcx),%%xmm4		/* t5 */\n\t"\
			"addpd	0x10(%%rbx),%%xmm1		/* t2 */\n\t"\
			"addpd	0x10(%%rcx),%%xmm5		/* t6 */\n\t"\
			"subpd	    (%%rbx),%%xmm2		/* t3 */\n\t"\
			"subpd	    (%%rcx),%%xmm6		/* t7 */\n\t"\
			"subpd	0x10(%%rbx),%%xmm3		/* t4 */\n\t"\
			"subpd	0x10(%%rcx),%%xmm7		/* t8 */\n\t"\
			"/* Finish radix-4 butterfly and store results into output-array slots: */\n\t"\
			"/* Swap _o1 <-> _o3: 0x20,0x30 <-> 0x60,0x70: */\n\t"\
			"subpd	%%xmm4,%%xmm0			/* ~t5 <- t1 -t5 */\n\t"\
			"subpd	%%xmm7,%%xmm2			/* ~t7 <- t3 -t8 */\n\t"\
			"subpd	%%xmm5,%%xmm1			/* ~t6 <- t2 -t6 */\n\t"\
			"subpd	%%xmm6,%%xmm3			/* ~t4 <- t4 -t7 */\n\t"\
			"movaps	%%xmm0,0x040(%%rsi)		/* <- ~t5 */\n\t"\
			"movaps	%%xmm2,0x020(%%rsi)		/* <- ~t7 */\n\t"\
			"movaps	%%xmm1,0x050(%%rsi)		/* <- ~t6 */\n\t"\
			"movaps	%%xmm3,0x070(%%rsi)		/* <- ~t4 */\n\t"\
			"addpd	%%xmm4,%%xmm4			/*          2*t5 */\n\t"\
			"addpd	%%xmm7,%%xmm7			/*          2*t8 */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*          2*t6 */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*          2*t7 */\n\t"\
			"addpd	%%xmm0,%%xmm4			/* ~t1 <- t1 +t5 */\n\t"\
			"addpd	%%xmm2,%%xmm7			/* ~t3 <- t3 +t8 */\n\t"\
			"addpd	%%xmm1,%%xmm5			/* ~t2 <- t2 +t6 */\n\t"\
			"addpd	%%xmm3,%%xmm6			/* ~t8 <- t4 +t7 */\n\t"\
			"movaps	%%xmm4,     (%%rsi)		/* <- ~t1 */\n\t"\
			"movaps	%%xmm7,0x060(%%rsi)		/* <- ~t3 */\n\t"\
			"movaps	%%xmm5,0x010(%%rsi)		/* <- ~t2 */\n\t"\
			"movaps	%%xmm6,0x030(%%rsi)		/* <- ~t8 */\n\t"\
			"\n\t"\
			"/*...and now do 4 radix-7 transforms...*/\n\t"\
			"\n\t"\
		"/* SSE2_RADIX_07_DFT(s1p00r,s1p04r,s1p08r,s1p12r,s1p16r,s1p20r,s1p24r,cc0,s1p00r,s1p08r,s1p16r,s1p24r,s1p04r,s1p12r,s1p20r): */\n\t"\
			"/*\n\t"\
			"t1r=A1r+A6r;	t2r=A2r+A5r;	t3r=A3r+A4r;\n\t"\
			"t6r=A1r-A6r;	t5r=A2r-A5r;	t4r=A3r-A4r;\n\t"\
			"*/\n\t"\
			"movq	%[__out],%%rdi	/* s1p00r */\n\t"\
			"movq	%[__cc0],%%rax	/* cc0 */\n\t"\
			"movq	$0x20,%%rbx\n\t"\
			"movq	$0x40,%%rcx\n\t"\
			"movq	$0x60,%%rdx\n\t"\
			"movq	$0x80,%%rsi\n\t"\
			"addq	%%rax,%%rbx		/* cc1 */\n\t"\
			"addq	%%rax,%%rcx		/* cc2 */\n\t"\
			"addq	%%rax,%%rdx		/* cc3 */\n\t"\
			"addq	%%rax,%%rsi		/* scratch storage */\n\t"\
			"movaps	0x080(%%rdi),%%xmm6			/* A1r = s1p04r */\n\t"\
			"movaps	0x300(%%rdi),%%xmm1			/* A6r = s1p24r */\n\t"\
			"movaps	0x100(%%rdi),%%xmm5			/* A2r = s1p08r */\n\t"\
			"movaps	0x280(%%rdi),%%xmm2			/* A5r = s1p20r */\n\t"\
			"movaps	0x180(%%rdi),%%xmm4			/* A3r = s1p12r */\n\t"\
			"movaps	0x200(%%rdi),%%xmm3			/* A4r = s1p16r */\n\t"\
			"\n\t"\
			"subpd	%%xmm1,%%xmm6			/* t6r = A1r-A6r */\n\t"\
			"addpd	%%xmm1,%%xmm1			/*         2*A6r */\n\t"\
			"addpd	%%xmm6,%%xmm1			/* t1r = A1r+A6r */\n\t"\
			"\n\t"\
			"subpd	%%xmm2,%%xmm5			/* t5r = A2r-A5r */\n\t"\
			"addpd	%%xmm2,%%xmm2			/*         2*A5r */\n\t"\
			"addpd	%%xmm5,%%xmm2			/* t2r = A2r+A5r */\n\t"\
			"\n\t"\
			"movaps	    (%%rdi),%%xmm0		/* Ar0 */\n\t"\
			"subpd	%%xmm3,%%xmm4			/* t4r = A3r-A4r */\n\t"\
			"addpd	%%xmm3,%%xmm3			/*         2*A4r */\n\t"\
			"addpd	%%xmm4,%%xmm3			/* t3r = A3r+A4r */\n\t"\
			"/*\n\t"\
			"	rt  = t1r+t2r+t3r;			\n\t"\
			"	B0r = rt + A0r;				\n\t"\
			"	t0r = rt*cx0 + A0r;			t3r=(t6r-t4r+t5r)*sx0;\n\t"\
			"	t1r = t1r-t2r;				t6r= t6r-t5r;\n\t"\
			"	t2r = t3r-t2r;				t5r= t4r+t5r;\n\t"\
			"	t3r =(t1r+t2r)*cx3;			t4r=(t5r-t6r)*sx3;\n\t"\
			"	t1r = t1r*cx1;				t6r= t6r*sx1;\n\t"\
			"	t2r = t2r*cx2;				t5r= t5r*sx2;\n\t"\
			"	tt  = t1r-t3r;				t6r= t4r+t6r;\n\t"\
			"	t2r = t2r-t3r;				t5r= t4r-t5r;\n\t"\
			"	\n\t"\
			"	t1r= t0r- tt-t2r;			t4r= t3r-t6r-t5r;\n\t"\
			"	t2r= t0r+t2r;				t5r= t3r+t5r;\n\t"\
			"	t0r= t0r+ tt;				t3r= t3r+t6r;\n\t"\
			"*/\n\t"\
			"movaps	%%xmm0,    (%%rsi)		/* cpy t0 into scratch sincos slot */\n\t"\
			"movaps	%%xmm6,0x10(%%rsi)		/* cpy t6 into scratch sincos slot */\n\t"\
			"addpd	%%xmm1,%%xmm0			/*~A0 = A0+t1 */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* cpy t5 */\n\t"\
			"addpd	%%xmm2,%%xmm3			/*~t3 = t3+t2 */\n\t"\
			"subpd	%%xmm4,%%xmm5			/*~t5 = t5-t4 */\n\t"\
			"subpd	%%xmm2,%%xmm1			/*~t1 = t1-t2 */\n\t"\
			"subpd	%%xmm7,%%xmm6			/*~t6 = t6-t5 */\n\t"\
			"addpd	%%xmm2,%%xmm2			/* 2*t2 */\n\t"\
			"addpd	%%xmm7,%%xmm4			/*~t5 = t4+t5 */\n\t"\
			"addpd	%%xmm3,%%xmm0			/* B0 */\n\t"\
			"addpd	0x10(%%rsi),%%xmm5		/* t3 = (t5-t4)+t6 */\n\t"\
			"subpd	%%xmm2,%%xmm3			/*~t2 =  (t2+t3) - 2*t2 = t3-t2 */\n\t"\
			"movaps	%%xmm4,%%xmm7			/* cpy t5 */\n\t"\
			"movaps	%%xmm0,     (%%rdi)		/* <-B0 = s1p00r, xmm0 FREE */\n\t"\
			"subpd	%%xmm6,%%xmm4			/* t4 = ~t5-~t6 */\n\t"\
			"movaps	%%xmm1,%%xmm2			/* cpy ~t1 */\n\t"\
			"subpd	    (%%rsi),%%xmm0		/* r = B0 - t0 */\n\t"\
			"mulpd	0x10(%%rax),%%xmm5		/*~t3 = t3*sx0 */\n\t"\
			"addpd	%%xmm3,%%xmm2			/* ~t1+~t2 */\n\t"\
			"mulpd	    (%%rcx),%%xmm3		/* t2 = t2*cx2 */\n\t"\
			"mulpd	0x10(%%rdx),%%xmm4		/*~t4 = t4*sx3 */\n\t"\
			"mulpd	    (%%rbx),%%xmm1		/* t1 = t1*cx1 */\n\t"\
			"mulpd	0x10(%%rbx),%%xmm6		/*~t6 = t6*sx1 */\n\t"\
			"mulpd	    (%%rax),%%xmm0		/* ~r = r*(cx0-1) */\n\t"\
			"mulpd	0x10(%%rcx),%%xmm7		/*~t5 = t5*sx2 */\n\t"\
			"mulpd	    (%%rdx),%%xmm2		/* t3 =(t1+t2)*cx3 */\n\t"\
			"addpd	    (%%rdi),%%xmm0		/* t0 =~r + B0 */\n\t"\
			"addpd	%%xmm4,%%xmm6			/*~t6 = t4+t6 */\n\t"\
			"subpd	%%xmm2,%%xmm1			/* tt = t1-t3 */\n\t"\
			"subpd	%%xmm7,%%xmm4			/*~t5 = t4-t5, xmm7 FREE */\n\t"\
			"subpd	%%xmm2,%%xmm3			/* t2 = t2-t3, xmm2 FREE */\n\t"\
			"\n\t"\
			"movaps	%%xmm0,%%xmm2			/* cpy t0 */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* cpy t3 */\n\t"\
			"addpd	%%xmm1,%%xmm0			/*~t0 = t0+tt */\n\t"\
			"addpd	%%xmm6,%%xmm5			/*~t3 = t3+t6 */\n\t"\
			"addpd	%%xmm3,%%xmm1			/*~tt = tt+t2 */\n\t"\
			"addpd	%%xmm4,%%xmm6			/*      t6+t5 */\n\t"\
			"addpd	%%xmm2,%%xmm3			/*~t2 = t2+t0 */\n\t"\
			"addpd	%%xmm7,%%xmm4			/*~t5 = t5+t3 */\n\t"\
			"subpd	%%xmm1,%%xmm2			/*~t1 = t0-tt-t2 */\n\t"\
			"subpd	%%xmm6,%%xmm7			/*~t4 = t3-t6-t5 */\n\t"\
			"movaps	%%xmm0,0x100(%%rdi)			/* B1 <- t0 = s1p08r */\n\t"\
			"movaps	%%xmm5,0x280(%%rdi)			/* B6 <- t3 = s1p20r */\n\t"\
			"movaps	%%xmm2,0x200(%%rdi)			/* B2 <- t1 = s1p16r */\n\t"\
			"movaps	%%xmm7,0x180(%%rdi)			/* B5 <- t4 = s1p12r */\n\t"\
			"movaps	%%xmm3,0x300(%%rdi)			/* B3 <- t2 = s1p24r */\n\t"\
			"movaps	%%xmm4,0x080(%%rdi)			/* B4 <- t5 = s1p04r */\n\t"\
			"/************************** Imaginary Parts: ******************************************/\n\t"\
			"movaps	0x090(%%rdi),%%xmm6		/* A1i = s1p04i */\n\t"\
			"movaps	0x310(%%rdi),%%xmm1		/* A6i = s1p24i */\n\t"\
			"movaps	0x110(%%rdi),%%xmm5		/* A2i = s1p08i */\n\t"\
			"movaps	0x290(%%rdi),%%xmm2		/* A5i = s1p20i */\n\t"\
			"movaps	0x190(%%rdi),%%xmm4		/* A3i = s1p12i */\n\t"\
			"movaps	0x210(%%rdi),%%xmm3		/* A4i = s1p16i */\n\t"\
			"\n\t"\
			"subpd	%%xmm1,%%xmm6			/* t6i = A1i-A6i */\n\t"\
			"addpd	%%xmm1,%%xmm1			/*         2*A6i */\n\t"\
			"addpd	%%xmm6,%%xmm1			/* t1i = A1i+A6i */\n\t"\
			"\n\t"\
			"subpd	%%xmm2,%%xmm5			/* t5i = A2i-A5i */\n\t"\
			"addpd	%%xmm2,%%xmm2			/*         2*A5i */\n\t"\
			"addpd	%%xmm5,%%xmm2			/* t2i = A2i+A5i */\n\t"\
			"\n\t"\
			"movaps	0x10(%%rdi),%%xmm0		/* Ai0 */\n\t"\
			"subpd	%%xmm3,%%xmm4			/* t4i = A3i-A4i */\n\t"\
			"addpd	%%xmm3,%%xmm3			/*         2*A4i */\n\t"\
			"addpd	%%xmm4,%%xmm3			/* t3i = A3i+A4i */\n\t"\
			"/*\n\t"\
			"	it  = t1i+t2i+t3i;			\n\t"\
			"	B0i = it + A0i;				\n\t"\
			"	t0i = it*cx0 + A0i;			t3i=(t6i-t4i+t5i)*sx0;\n\t"\
			"	t1i = t1i-t2i;				t6i= t6i-t5i;\n\t"\
			"	t2i = t2i-t3i;				t5i= t4i+t5i;\n\t"\
			"	t3i =(t1i-t2i)*cx3;			t4i=(t5i-t6i)*sx3;\n\t"\
			"	t1i = t1i*cx1;				t6i= t6i*sx1;\n\t"\
			"	t2i = t2i*cx2;				t5i= t5i*sx2;\n\t"\
			"	it  = t1i-t3i;				t6i= t4i+t6i;\n\t"\
			"	t2i = t2i-t3i;				t5i= t4i-t5i;\n\t"\
			"	\n\t"\
			"	t1i= t0i- it-t2i;			t4i= t3i-t6i-t5i;\n\t"\
			"	t2i= t0i+t2i;				t5i= t3i+t5i;\n\t"\
			"	t0i= t0i+ it;				t3i= t3i+t6i;\n\t"\
			"*/\n\t"\
			"movaps	%%xmm0,    (%%rsi)		/* cpy t0 into scratch sincos slot */\n\t"\
			"movaps	%%xmm6,0x10(%%rsi)		/* cpy t6 into scratch sincos slot */\n\t"\
			"addpd	%%xmm1,%%xmm0			/*~A0 = A0+t1 */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* cpy t5 */\n\t"\
			"addpd	%%xmm2,%%xmm3			/*~t3 = t3+t2 */\n\t"\
			"subpd	%%xmm4,%%xmm5			/*~t5 = t5-t4 */\n\t"\
			"subpd	%%xmm2,%%xmm1			/*~t1 = t1-t2 */\n\t"\
			"subpd	%%xmm7,%%xmm6			/*~t6 = t6-t5 */\n\t"\
			"addpd	%%xmm2,%%xmm2			/* 2*t2 */\n\t"\
			"addpd	%%xmm7,%%xmm4			/*~t5 = t4+t5 */\n\t"\
			"addpd	%%xmm3,%%xmm0			/* B0 */\n\t"\
			"addpd	0x10(%%rsi),%%xmm5		/* t3 = (t5-t4)+t6 */\n\t"\
			"subpd	%%xmm2,%%xmm3			/*~t2 =  (t2+t3) - 2*t2 = t3-t2 */\n\t"\
			"movaps	%%xmm4,%%xmm7			/* cpy t5 */\n\t"\
			"movaps	%%xmm0,0x010(%%rdi)		/* <-B0 = s1p00i, xmm0 FREE */\n\t"\
			"subpd	%%xmm6,%%xmm4			/* t4 = ~t5-~t6 */\n\t"\
			"movaps	%%xmm1,%%xmm2			/* cpy ~t1 */\n\t"\
			"subpd	    (%%rsi),%%xmm0		/* r = B0 - t0 */\n\t"\
			"mulpd	0x10(%%rax),%%xmm5		/*~t3 = t3*sx0 */\n\t"\
			"addpd	%%xmm3,%%xmm2			/* ~t1+~t2 */\n\t"\
			"mulpd	    (%%rcx),%%xmm3		/* t2 = t2*cx2 */\n\t"\
			"mulpd	0x10(%%rdx),%%xmm4		/*~t4 = t4*sx3 */\n\t"\
			"mulpd	    (%%rbx),%%xmm1		/* t1 = t1*cx1 */\n\t"\
			"mulpd	0x10(%%rbx),%%xmm6		/*~t6 = t6*sx1 */\n\t"\
			"mulpd	    (%%rax),%%xmm0		/* ~r = r*(cx0-1) */\n\t"\
			"mulpd	0x10(%%rcx),%%xmm7		/*~t5 = t5*sx2 */\n\t"\
			"mulpd	    (%%rdx),%%xmm2		/* t3 =(t1+t2)*cx3 */\n\t"\
			"addpd	0x10(%%rdi),%%xmm0		/* t0 =~r + B0 */\n\t"\
			"addpd	%%xmm4,%%xmm6			/*~t6 = t4+t6 */\n\t"\
			"subpd	%%xmm2,%%xmm1			/* tt = t1-t3 */\n\t"\
			"subpd	%%xmm7,%%xmm4			/*~t5 = t4-t5, xmm7 FREE */\n\t"\
			"subpd	%%xmm2,%%xmm3			/* t2 = t2-t3, xmm2 FREE */\n\t"\
			"movaps	%%xmm0,%%xmm2			/* cpy t0 */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* cpy t3 */\n\t"\
			"addpd	%%xmm1,%%xmm0			/*~t0 = t0+tt */\n\t"\
			"addpd	%%xmm6,%%xmm5			/*~t3 = t3+t6 */\n\t"\
			"addpd	%%xmm3,%%xmm1			/*~tt = tt+t2 */\n\t"\
			"addpd	%%xmm4,%%xmm6			/*      t6+t5 */\n\t"\
			"addpd	%%xmm2,%%xmm3			/*~t2 = t2+t0 */\n\t"\
			"addpd	%%xmm7,%%xmm4			/*~t5 = t5+t3 */\n\t"\
			"subpd	%%xmm1,%%xmm2			/*~t1 = t0-tt-t2, xmm1 FREE */\n\t"\
			"subpd	%%xmm6,%%xmm7			/*~t4 = t3-t6-t5, xmm6 FREE */\n\t"\
			"/*\n\t"\
			"	B1r =t0r-t3i;				B1i*=t0i+t3r;\n\t"\
			"	B2r =t1r-t4i;				B2i*=t1i+t4r;\n\t"\
			"	B3r*=t2r+t5i;				B3i =t2i-t5r;\n\t"\
			"	B4r*=t2r-t5i;				B4i =t2i+t5r;\n\t"\
			"	B5r =t1r+t4i;				B5i*=t1i-t4r;\n\t"\
			"	B6r =t0r+t3i;				B6i*=t0i-t3r;\n\t"\
			"*/\n\t"\
			"/* xmm1,6 FREE */\n\t"\
			"movaps	0x100(%%rdi),%%xmm1		/* t0r */\n\t"\
			"movaps	0x280(%%rdi),%%xmm6		/* t3r */\n\t"\
			"subpd	%%xmm5,%%xmm1			/* B1r =t0r-t3i */\n\t"\
			"subpd	%%xmm6,%%xmm0			/* B6i =t0i-t3r */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*        2*t3i */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*        2*t3r */\n\t"\
			"addpd	%%xmm1,%%xmm5			/* B6r =t0r+t3i */\n\t"\
			"addpd	%%xmm0,%%xmm6			/* B1i =t0i+t3r */\n\t"\
			"movaps	%%xmm1,0x100(%%rdi)		/* <-B1r = s1p08r */\n\t"\
			"movaps	%%xmm0,0x290(%%rdi)		/* <-B6i = s1p20r */\n\t"\
			"movaps	%%xmm5,0x280(%%rdi)		/* <-B6r = s1p20r */\n\t"\
			"movaps	%%xmm6,0x110(%%rdi)		/* <-B1i = s1p08i */\n\t"\
			"\n\t"\
			"movaps	0x200(%%rdi),%%xmm1		/* t1r */\n\t"\
			"movaps	0x180(%%rdi),%%xmm6		/* t4r */\n\t"\
			"subpd	%%xmm7,%%xmm1			/* B2r =t1r-t4i */\n\t"\
			"subpd	%%xmm6,%%xmm2			/* B5i =t1i-t4r */\n\t"\
			"addpd	%%xmm7,%%xmm7			/*        2*t4i */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*        2*t4r */\n\t"\
			"addpd	%%xmm1,%%xmm7			/* B5r =t1r+t4i */\n\t"\
			"addpd	%%xmm2,%%xmm6			/* B2i =t1i+t4r */\n\t"\
			"movaps	%%xmm1,0x200(%%rdi)		/* <-B2r = s1p16r */\n\t"\
			"movaps	%%xmm2,0x190(%%rdi)		/* <-B5i = s1p12r */\n\t"\
			"movaps	%%xmm7,0x180(%%rdi)		/* <-B5r = s1p12r */\n\t"\
			"movaps	%%xmm6,0x210(%%rdi)		/* <-B2i = s1p16i */\n\t"\
			"\n\t"\
			"/* Note the order reversal on this pair of outputs: */\n\t"\
			"movaps	0x300(%%rdi),%%xmm0		/* t2r */\n\t"\
			"movaps	0x080(%%rdi),%%xmm5		/* t5r */\n\t"\
			"subpd	%%xmm4,%%xmm0			/* B4r =t2r-t5i */\n\t"\
			"subpd	%%xmm5,%%xmm3			/* B3i =t2i-t5r */\n\t"\
			"addpd	%%xmm4,%%xmm4			/*        2*t5i */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*        2*t5r */\n\t"\
			"addpd	%%xmm0,%%xmm4			/* B3r =t2r+t5i */\n\t"\
			"addpd	%%xmm3,%%xmm5			/* B4i =t2i+t5r */\n\t"\
			"movaps	%%xmm0,0x080(%%rdi)		/* <-B4r = s1p04r */\n\t"\
			"movaps	%%xmm3,0x310(%%rdi)		/* <-B3i = s1p24r */\n\t"\
			"movaps	%%xmm4,0x300(%%rdi)		/* <-B3r = s1p24r */\n\t"\
			"movaps	%%xmm5,0x090(%%rdi)		/* <-B4i = s1p04i */\n\t"\
			"\n\t"\
		"/* SSE2_RADIX_07_DFT(s1p01r,s1p05r,s1p09r,s1p13r,s1p17r,s1p21r,s1p25r,cc0,s1p21r,s1p01r,s1p09r,s1p17r,s1p25r,s1p05r,s1p13r): */\n\t"\
			"addq	$0x20,%%rdi	/* s1p01r */\n\t"\
			"movaps	0x080(%%rdi),%%xmm6			/* A1r = s1p05r */\n\t"\
			"movaps	0x300(%%rdi),%%xmm1			/* A6r = s1p25r */\n\t"\
			"movaps	0x100(%%rdi),%%xmm5			/* A2r = s1p09r */\n\t"\
			"movaps	0x280(%%rdi),%%xmm2			/* A5r = s1p21r */\n\t"\
			"movaps	0x180(%%rdi),%%xmm4			/* A3r = s1p13r */\n\t"\
			"movaps	0x200(%%rdi),%%xmm3			/* A4r = s1p17r */\n\t"\
			"\n\t"\
			"subpd	%%xmm1,%%xmm6			/* t6r = A1r-A6r */\n\t"\
			"addpd	%%xmm1,%%xmm1			/*         2*A6r */\n\t"\
			"addpd	%%xmm6,%%xmm1			/* t1r = A1r+A6r */\n\t"\
			"\n\t"\
			"subpd	%%xmm2,%%xmm5			/* t5r = A2r-A5r */\n\t"\
			"addpd	%%xmm2,%%xmm2			/*         2*A5r */\n\t"\
			"addpd	%%xmm5,%%xmm2			/* t2r = A2r+A5r */\n\t"\
			"\n\t"\
			"movaps	    (%%rdi),%%xmm0		/* Ar0 */\n\t"\
			"subpd	%%xmm3,%%xmm4			/* t4r = A3r-A4r */\n\t"\
			"addpd	%%xmm3,%%xmm3			/*         2*A4r */\n\t"\
			"addpd	%%xmm4,%%xmm3			/* t3r = A3r+A4r */\n\t"\
			"\n\t"\
			"movaps	%%xmm0,    (%%rsi)		/* cpy t0 into scratch sincos slot */\n\t"\
			"movaps	%%xmm6,0x10(%%rsi)		/* cpy t6 into scratch sincos slot */\n\t"\
			"addpd	%%xmm1,%%xmm0			/*~A0 = A0+t1 */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* cpy t5 */\n\t"\
			"addpd	%%xmm2,%%xmm3			/*~t3 = t3+t2 */\n\t"\
			"subpd	%%xmm4,%%xmm5			/*~t5 = t5-t4 */\n\t"\
			"subpd	%%xmm2,%%xmm1			/*~t1 = t1-t2 */\n\t"\
			"subpd	%%xmm7,%%xmm6			/*~t6 = t6-t5 */\n\t"\
			"addpd	%%xmm2,%%xmm2			/* 2*t2 */\n\t"\
			"addpd	%%xmm7,%%xmm4			/*~t5 = t4+t5 */\n\t"\
			"addpd	%%xmm3,%%xmm0			/* B0 */\n\t"\
			"addpd	0x10(%%rsi),%%xmm5		/* t3 = (t5-t4)+t6 */\n\t"\
			"subpd	%%xmm2,%%xmm3			/*~t2 =  (t2+t3) - 2*t2 = t3-t2 */\n\t"\
			"movaps	%%xmm4,%%xmm7			/* cpy t5 */\n\t"\
			"movaps	%%xmm0,0x280(%%rdi)		/* <-B0 = s1p21r, xmm0 FREE */\n\t"\
			"subpd	%%xmm6,%%xmm4			/* t4 = ~t5-~t6 */\n\t"\
			"movaps	%%xmm1,%%xmm2			/* cpy ~t1 */\n\t"\
			"subpd	    (%%rsi),%%xmm0		/* r = B0 - t0 */\n\t"\
			"mulpd	0x10(%%rax),%%xmm5		/*~t3 = t3*sx0 */\n\t"\
			"addpd	%%xmm3,%%xmm2			/* ~t1+~t2 */\n\t"\
			"mulpd	    (%%rcx),%%xmm3		/* t2 = t2*cx2 */\n\t"\
			"mulpd	0x10(%%rdx),%%xmm4		/*~t4 = t4*sx3 */\n\t"\
			"mulpd	    (%%rbx),%%xmm1		/* t1 = t1*cx1 */\n\t"\
			"mulpd	0x10(%%rbx),%%xmm6		/*~t6 = t6*sx1 */\n\t"\
			"mulpd	    (%%rax),%%xmm0		/* ~r = r*(cx0-1) */\n\t"\
			"mulpd	0x10(%%rcx),%%xmm7		/*~t5 = t5*sx2 */\n\t"\
			"mulpd	    (%%rdx),%%xmm2		/* t3 =(t1+t2)*cx3 */\n\t"\
			"addpd	0x280(%%rdi),%%xmm0		/* t0 =~r + B0 */\n\t"\
			"addpd	%%xmm4,%%xmm6			/*~t6 = t4+t6 */\n\t"\
			"subpd	%%xmm2,%%xmm1			/* tt = t1-t3 */\n\t"\
			"subpd	%%xmm7,%%xmm4			/*~t5 = t4-t5, xmm7 FREE */\n\t"\
			"subpd	%%xmm2,%%xmm3			/* t2 = t2-t3, xmm2 FREE */\n\t"\
			"\n\t"\
			"movaps	%%xmm0,%%xmm2			/* cpy t0 */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* cpy t3 */\n\t"\
			"addpd	%%xmm1,%%xmm0			/*~t0 = t0+tt */\n\t"\
			"addpd	%%xmm6,%%xmm5			/*~t3 = t3+t6 */\n\t"\
			"addpd	%%xmm3,%%xmm1			/*~tt = tt+t2 */\n\t"\
			"addpd	%%xmm4,%%xmm6			/*      t6+t5 */\n\t"\
			"addpd	%%xmm2,%%xmm3			/*~t2 = t2+t0 */\n\t"\
			"addpd	%%xmm7,%%xmm4			/*~t5 = t5+t3 */\n\t"\
			"subpd	%%xmm1,%%xmm2			/*~t1 = t0-tt-t2 */\n\t"\
			"subpd	%%xmm6,%%xmm7			/*~t4 = t3-t6-t5 */\n\t"\
			"movaps	%%xmm0,     (%%rdi)			/* B1 <- t0 = s1p01r */\n\t"\
			"movaps	%%xmm5,0x180(%%rdi)			/* B6 <- t3 = s1p13r */\n\t"\
			"movaps	%%xmm2,0x100(%%rdi)			/* B2 <- t1 = s1p09r */\n\t"\
			"movaps	%%xmm7,0x080(%%rdi)			/* B5 <- t4 = s1p05r */\n\t"\
			"movaps	%%xmm3,0x200(%%rdi)			/* B3 <- t2 = s1p17r */\n\t"\
			"movaps	%%xmm4,0x300(%%rdi)			/* B4 <- t5 = s1p25r */\n\t"\
			"/************************** Imaginary Parts: ******************************************/\n\t"\
			"movaps	0x090(%%rdi),%%xmm6		/* A1i = s1p05i */\n\t"\
			"movaps	0x310(%%rdi),%%xmm1		/* A6i = s1p25i */\n\t"\
			"movaps	0x110(%%rdi),%%xmm5		/* A2i = s1p09i */\n\t"\
			"movaps	0x290(%%rdi),%%xmm2		/* A5i = s1p21i */\n\t"\
			"movaps	0x190(%%rdi),%%xmm4		/* A3i = s1p13i */\n\t"\
			"movaps	0x210(%%rdi),%%xmm3		/* A4i = s1p17i */\n\t"\
			"\n\t"\
			"subpd	%%xmm1,%%xmm6			/* t6i = A1i-A6i */\n\t"\
			"addpd	%%xmm1,%%xmm1			/*         2*A6i */\n\t"\
			"addpd	%%xmm6,%%xmm1			/* t1i = A1i+A6i */\n\t"\
			"\n\t"\
			"subpd	%%xmm2,%%xmm5			/* t5i = A2i-A5i */\n\t"\
			"addpd	%%xmm2,%%xmm2			/*         2*A5i */\n\t"\
			"addpd	%%xmm5,%%xmm2			/* t2i = A2i+A5i */\n\t"\
			"\n\t"\
			"movaps	0x10(%%rdi),%%xmm0		/* Ai0 */\n\t"\
			"subpd	%%xmm3,%%xmm4			/* t4i = A3i-A4i */\n\t"\
			"addpd	%%xmm3,%%xmm3			/*         2*A4i */\n\t"\
			"addpd	%%xmm4,%%xmm3			/* t3i = A3i+A4i */\n\t"\
			"\n\t"\
			"movaps	%%xmm0,    (%%rsi)		/* cpy t0 into scratch sincos slot */\n\t"\
			"movaps	%%xmm6,0x10(%%rsi)		/* cpy t6 into scratch sincos slot */\n\t"\
			"addpd	%%xmm1,%%xmm0			/*~A0 = A0+t1 */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* cpy t5 */\n\t"\
			"addpd	%%xmm2,%%xmm3			/*~t3 = t3+t2 */\n\t"\
			"subpd	%%xmm4,%%xmm5			/*~t5 = t5-t4 */\n\t"\
			"subpd	%%xmm2,%%xmm1			/*~t1 = t1-t2 */\n\t"\
			"subpd	%%xmm7,%%xmm6			/*~t6 = t6-t5 */\n\t"\
			"addpd	%%xmm2,%%xmm2			/* 2*t2 */\n\t"\
			"addpd	%%xmm7,%%xmm4			/*~t5 = t4+t5 */\n\t"\
			"addpd	%%xmm3,%%xmm0			/* B0 */\n\t"\
			"addpd	0x10(%%rsi),%%xmm5		/* t3 = (t5-t4)+t6 */\n\t"\
			"subpd	%%xmm2,%%xmm3			/*~t2 =  (t2+t3) - 2*t2 = t3-t2 */\n\t"\
			"movaps	%%xmm4,%%xmm7			/* cpy t5 */\n\t"\
			"movaps	%%xmm0,0x290(%%rdi)		/* <-B0 = s1p21i, xmm0 FREE */\n\t"\
			"subpd	%%xmm6,%%xmm4			/* t4 = ~t5-~t6 */\n\t"\
			"movaps	%%xmm1,%%xmm2			/* cpy ~t1 */\n\t"\
			"subpd	    (%%rsi),%%xmm0		/* r = B0 - t0 */\n\t"\
			"mulpd	0x10(%%rax),%%xmm5		/*~t3 = t3*sx0 */\n\t"\
			"addpd	%%xmm3,%%xmm2			/* ~t1+~t2 */\n\t"\
			"mulpd	    (%%rcx),%%xmm3		/* t2 = t2*cx2 */\n\t"\
			"mulpd	0x10(%%rdx),%%xmm4		/*~t4 = t4*sx3 */\n\t"\
			"mulpd	    (%%rbx),%%xmm1		/* t1 = t1*cx1 */\n\t"\
			"mulpd	0x10(%%rbx),%%xmm6		/*~t6 = t6*sx1 */\n\t"\
			"mulpd	    (%%rax),%%xmm0		/* ~r = r*(cx0-1) */\n\t"\
			"mulpd	0x10(%%rcx),%%xmm7		/*~t5 = t5*sx2 */\n\t"\
			"mulpd	    (%%rdx),%%xmm2		/* t3 =(t1+t2)*cx3 */\n\t"\
			"addpd	0x290(%%rdi),%%xmm0		/* t0 =~r + B0 */\n\t"\
			"addpd	%%xmm4,%%xmm6			/*~t6 = t4+t6 */\n\t"\
			"subpd	%%xmm2,%%xmm1			/* tt = t1-t3 */\n\t"\
			"subpd	%%xmm7,%%xmm4			/*~t5 = t4-t5, xmm7 FREE */\n\t"\
			"subpd	%%xmm2,%%xmm3			/* t2 = t2-t3, xmm2 FREE */\n\t"\
			"movaps	%%xmm0,%%xmm2			/* cpy t0 */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* cpy t3 */\n\t"\
			"addpd	%%xmm1,%%xmm0			/*~t0 = t0+tt */\n\t"\
			"addpd	%%xmm6,%%xmm5			/*~t3 = t3+t6 */\n\t"\
			"addpd	%%xmm3,%%xmm1			/*~tt = tt+t2 */\n\t"\
			"addpd	%%xmm4,%%xmm6			/*      t6+t5 */\n\t"\
			"addpd	%%xmm2,%%xmm3			/*~t2 = t2+t0 */\n\t"\
			"addpd	%%xmm7,%%xmm4			/*~t5 = t5+t3 */\n\t"\
			"subpd	%%xmm1,%%xmm2			/*~t1 = t0-tt-t2, xmm1 FREE */\n\t"\
			"subpd	%%xmm6,%%xmm7			/*~t4 = t3-t6-t5, xmm6 FREE */\n\t"\
			"\n\t"\
			"/* xmm1,6 FREE */\n\t"\
			"movaps	     (%%rdi),%%xmm1		/* t0r */\n\t"\
			"movaps	0x180(%%rdi),%%xmm6		/* t3r */\n\t"\
			"subpd	%%xmm5,%%xmm1			/* B1r =t0r-t3i */\n\t"\
			"subpd	%%xmm6,%%xmm0			/* B6i =t0i-t3r */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*        2*t3i */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*        2*t3r */\n\t"\
			"addpd	%%xmm1,%%xmm5			/* B6r =t0r+t3i */\n\t"\
			"addpd	%%xmm0,%%xmm6			/* B1i =t0i+t3r */\n\t"\
			"movaps	%%xmm1,     (%%rdi)		/* <-B1r = s1p08r */\n\t"\
			"movaps	%%xmm0,0x190(%%rdi)		/* <-B6i = s1p20r */\n\t"\
			"movaps	%%xmm5,0x180(%%rdi)		/* <-B6r = s1p20r */\n\t"\
			"movaps	%%xmm6,0x010(%%rdi)		/* <-B1i = s1p08i */\n\t"\
			"\n\t"\
			"movaps	0x100(%%rdi),%%xmm1		/* t1r */\n\t"\
			"movaps	0x080(%%rdi),%%xmm6		/* t4r */\n\t"\
			"subpd	%%xmm7,%%xmm1			/* B2r =t1r-t4i */\n\t"\
			"subpd	%%xmm6,%%xmm2			/* B5i =t1i-t4r */\n\t"\
			"addpd	%%xmm7,%%xmm7			/*        2*t4i */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*        2*t4r */\n\t"\
			"addpd	%%xmm1,%%xmm7			/* B5r =t1r+t4i */\n\t"\
			"addpd	%%xmm2,%%xmm6			/* B2i =t1i+t4r */\n\t"\
			"movaps	%%xmm1,0x100(%%rdi)		/* <-B2r = s1p16r */\n\t"\
			"movaps	%%xmm2,0x090(%%rdi)		/* <-B5i = s1p12r */\n\t"\
			"movaps	%%xmm7,0x080(%%rdi)		/* <-B5r = s1p12r */\n\t"\
			"movaps	%%xmm6,0x110(%%rdi)		/* <-B2i = s1p16i */\n\t"\
			"\n\t"\
			"/* Note the order reversal on this pair of outputs: */\n\t"\
			"movaps	0x200(%%rdi),%%xmm0		/* t2r */\n\t"\
			"movaps	0x300(%%rdi),%%xmm5		/* t5r */\n\t"\
			"subpd	%%xmm4,%%xmm0			/* B4r =t2r-t5i */\n\t"\
			"subpd	%%xmm5,%%xmm3			/* B3i =t2i-t5r */\n\t"\
			"addpd	%%xmm4,%%xmm4			/*        2*t5i */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*        2*t5r */\n\t"\
			"addpd	%%xmm0,%%xmm4			/* B3r =t2r+t5i */\n\t"\
			"addpd	%%xmm3,%%xmm5			/* B4i =t2i+t5r */\n\t"\
			"movaps	%%xmm0,0x300(%%rdi)		/* <-B4r = s1p04r */\n\t"\
			"movaps	%%xmm3,0x210(%%rdi)		/* <-B3i = s1p24r */\n\t"\
			"movaps	%%xmm4,0x200(%%rdi)		/* <-B3r = s1p24r */\n\t"\
			"movaps	%%xmm5,0x310(%%rdi)		/* <-B4i = s1p04i */\n\t"\
			"\n\t"\
		"/* SSE2_RADIX_07_DFT(s1p02r,s1p06r,s1p10r,s1p14r,s1p18r,s1p22r,s1p26r,cc0,s1p14r,s1p22r,s1p02r,s1p10r,s1p18r,s1p26r,s1p06r): */\n\t"\
			"addq	$0x20,%%rdi	/* s1p02r */\n\t"\
			"movaps	0x080(%%rdi),%%xmm6			/* A1r = s1p06r */\n\t"\
			"movaps	0x300(%%rdi),%%xmm1			/* A6r = s1p26r */\n\t"\
			"movaps	0x100(%%rdi),%%xmm5			/* A2r = s1p10r */\n\t"\
			"movaps	0x280(%%rdi),%%xmm2			/* A5r = s1p22r */\n\t"\
			"movaps	0x180(%%rdi),%%xmm4			/* A3r = s1p14r */\n\t"\
			"movaps	0x200(%%rdi),%%xmm3			/* A4r = s1p18r */\n\t"\
			"\n\t"\
			"subpd	%%xmm1,%%xmm6			/* t6r = A1r-A6r */\n\t"\
			"addpd	%%xmm1,%%xmm1			/*         2*A6r */\n\t"\
			"addpd	%%xmm6,%%xmm1			/* t1r = A1r+A6r */\n\t"\
			"\n\t"\
			"subpd	%%xmm2,%%xmm5			/* t5r = A2r-A5r */\n\t"\
			"addpd	%%xmm2,%%xmm2			/*         2*A5r */\n\t"\
			"addpd	%%xmm5,%%xmm2			/* t2r = A2r+A5r */\n\t"\
			"\n\t"\
			"movaps	    (%%rdi),%%xmm0		/* Ar0 */\n\t"\
			"subpd	%%xmm3,%%xmm4			/* t4r = A3r-A4r */\n\t"\
			"addpd	%%xmm3,%%xmm3			/*         2*A4r */\n\t"\
			"addpd	%%xmm4,%%xmm3			/* t3r = A3r+A4r */\n\t"\
			"\n\t"\
			"movaps	%%xmm0,    (%%rsi)		/* cpy t0 into scratch sincos slot */\n\t"\
			"movaps	%%xmm6,0x10(%%rsi)		/* cpy t6 into scratch sincos slot */\n\t"\
			"addpd	%%xmm1,%%xmm0			/*~A0 = A0+t1 */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* cpy t5 */\n\t"\
			"addpd	%%xmm2,%%xmm3			/*~t3 = t3+t2 */\n\t"\
			"subpd	%%xmm4,%%xmm5			/*~t5 = t5-t4 */\n\t"\
			"subpd	%%xmm2,%%xmm1			/*~t1 = t1-t2 */\n\t"\
			"subpd	%%xmm7,%%xmm6			/*~t6 = t6-t5 */\n\t"\
			"addpd	%%xmm2,%%xmm2			/* 2*t2 */\n\t"\
			"addpd	%%xmm7,%%xmm4			/*~t5 = t4+t5 */\n\t"\
			"addpd	%%xmm3,%%xmm0			/* B0 */\n\t"\
			"addpd	0x10(%%rsi),%%xmm5		/* t3 = (t5-t4)+t6 */\n\t"\
			"subpd	%%xmm2,%%xmm3			/*~t2 =  (t2+t3) - 2*t2 = t3-t2 */\n\t"\
			"movaps	%%xmm4,%%xmm7			/* cpy t5 */\n\t"\
			"movaps	%%xmm0,0x180(%%rdi)		/* <-B0 = s1p14r, xmm0 FREE */\n\t"\
			"subpd	%%xmm6,%%xmm4			/* t4 = ~t5-~t6 */\n\t"\
			"movaps	%%xmm1,%%xmm2			/* cpy ~t1 */\n\t"\
			"subpd	    (%%rsi),%%xmm0		/* r = B0 - t0 */\n\t"\
			"mulpd	0x10(%%rax),%%xmm5		/*~t3 = t3*sx0 */\n\t"\
			"addpd	%%xmm3,%%xmm2			/* ~t1+~t2 */\n\t"\
			"mulpd	    (%%rcx),%%xmm3		/* t2 = t2*cx2 */\n\t"\
			"mulpd	0x10(%%rdx),%%xmm4		/*~t4 = t4*sx3 */\n\t"\
			"mulpd	    (%%rbx),%%xmm1		/* t1 = t1*cx1 */\n\t"\
			"mulpd	0x10(%%rbx),%%xmm6		/*~t6 = t6*sx1 */\n\t"\
			"mulpd	    (%%rax),%%xmm0		/* ~r = r*(cx0-1) */\n\t"\
			"mulpd	0x10(%%rcx),%%xmm7		/*~t5 = t5*sx2 */\n\t"\
			"mulpd	    (%%rdx),%%xmm2		/* t3 =(t1+t2)*cx3 */\n\t"\
			"addpd	0x180(%%rdi),%%xmm0		/* t0 =~r + B0 */\n\t"\
			"addpd	%%xmm4,%%xmm6			/*~t6 = t4+t6 */\n\t"\
			"subpd	%%xmm2,%%xmm1			/* tt = t1-t3 */\n\t"\
			"subpd	%%xmm7,%%xmm4			/*~t5 = t4-t5, xmm7 FREE */\n\t"\
			"subpd	%%xmm2,%%xmm3			/* t2 = t2-t3, xmm2 FREE */\n\t"\
			"\n\t"\
			"movaps	%%xmm0,%%xmm2			/* cpy t0 */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* cpy t3 */\n\t"\
			"addpd	%%xmm1,%%xmm0			/*~t0 = t0+tt */\n\t"\
			"addpd	%%xmm6,%%xmm5			/*~t3 = t3+t6 */\n\t"\
			"addpd	%%xmm3,%%xmm1			/*~tt = tt+t2 */\n\t"\
			"addpd	%%xmm4,%%xmm6			/*      t6+t5 */\n\t"\
			"addpd	%%xmm2,%%xmm3			/*~t2 = t2+t0 */\n\t"\
			"addpd	%%xmm7,%%xmm4			/*~t5 = t5+t3 */\n\t"\
			"subpd	%%xmm1,%%xmm2			/*~t1 = t0-tt-t2 */\n\t"\
			"subpd	%%xmm6,%%xmm7			/*~t4 = t3-t6-t5 */\n\t"\
			"movaps	%%xmm0,0x280(%%rdi)			/* B1 <- t0 = s1p22r */\n\t"\
			"movaps	%%xmm5,0x080(%%rdi)			/* B6 <- t3 = s1p06r */\n\t"\
			"movaps	%%xmm2,     (%%rdi)			/* B2 <- t1 = s1p02r */\n\t"\
			"movaps	%%xmm7,0x300(%%rdi)			/* B5 <- t4 = s1p26r */\n\t"\
			"movaps	%%xmm3,0x100(%%rdi)			/* B3 <- t2 = s1p10r */\n\t"\
			"movaps	%%xmm4,0x200(%%rdi)			/* B4 <- t5 = s1p18r */\n\t"\
			"/************************** Imaginary Parts: ******************************************/\n\t"\
			"movaps	0x090(%%rdi),%%xmm6		/* A1i = s1p06i */\n\t"\
			"movaps	0x310(%%rdi),%%xmm1		/* A6i = s1p26i */\n\t"\
			"movaps	0x110(%%rdi),%%xmm5		/* A2i = s1p10i */\n\t"\
			"movaps	0x290(%%rdi),%%xmm2		/* A5i = s1p22i */\n\t"\
			"movaps	0x190(%%rdi),%%xmm4		/* A3i = s1p14i */\n\t"\
			"movaps	0x210(%%rdi),%%xmm3		/* A4i = s1p18i */\n\t"\
			"\n\t"\
			"subpd	%%xmm1,%%xmm6			/* t6i = A1i-A6i */\n\t"\
			"addpd	%%xmm1,%%xmm1			/*         2*A6i */\n\t"\
			"addpd	%%xmm6,%%xmm1			/* t1i = A1i+A6i */\n\t"\
			"\n\t"\
			"subpd	%%xmm2,%%xmm5			/* t5i = A2i-A5i */\n\t"\
			"addpd	%%xmm2,%%xmm2			/*         2*A5i */\n\t"\
			"addpd	%%xmm5,%%xmm2			/* t2i = A2i+A5i */\n\t"\
			"\n\t"\
			"movaps	0x10(%%rdi),%%xmm0		/* Ai0 */\n\t"\
			"subpd	%%xmm3,%%xmm4			/* t4i = A3i-A4i */\n\t"\
			"addpd	%%xmm3,%%xmm3			/*         2*A4i */\n\t"\
			"addpd	%%xmm4,%%xmm3			/* t3i = A3i+A4i */\n\t"\
			"\n\t"\
			"movaps	%%xmm0,    (%%rsi)		/* cpy t0 into scratch sincos slot */\n\t"\
			"movaps	%%xmm6,0x10(%%rsi)		/* cpy t6 into scratch sincos slot */\n\t"\
			"addpd	%%xmm1,%%xmm0			/*~A0 = A0+t1 */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* cpy t5 */\n\t"\
			"addpd	%%xmm2,%%xmm3			/*~t3 = t3+t2 */\n\t"\
			"subpd	%%xmm4,%%xmm5			/*~t5 = t5-t4 */\n\t"\
			"subpd	%%xmm2,%%xmm1			/*~t1 = t1-t2 */\n\t"\
			"subpd	%%xmm7,%%xmm6			/*~t6 = t6-t5 */\n\t"\
			"addpd	%%xmm2,%%xmm2			/* 2*t2 */\n\t"\
			"addpd	%%xmm7,%%xmm4			/*~t5 = t4+t5 */\n\t"\
			"addpd	%%xmm3,%%xmm0			/* B0 */\n\t"\
			"addpd	0x10(%%rsi),%%xmm5		/* t3 = (t5-t4)+t6 */\n\t"\
			"subpd	%%xmm2,%%xmm3			/*~t2 =  (t2+t3) - 2*t2 = t3-t2 */\n\t"\
			"movaps	%%xmm4,%%xmm7			/* cpy t5 */\n\t"\
			"movaps	%%xmm0,0x190(%%rdi)		/* <-B0 = s1p14i, xmm0 FREE */\n\t"\
			"subpd	%%xmm6,%%xmm4			/* t4 = ~t5-~t6 */\n\t"\
			"movaps	%%xmm1,%%xmm2			/* cpy ~t1 */\n\t"\
			"subpd	    (%%rsi),%%xmm0		/* r = B0 - t0 */\n\t"\
			"mulpd	0x10(%%rax),%%xmm5		/*~t3 = t3*sx0 */\n\t"\
			"addpd	%%xmm3,%%xmm2			/* ~t1+~t2 */\n\t"\
			"mulpd	    (%%rcx),%%xmm3		/* t2 = t2*cx2 */\n\t"\
			"mulpd	0x10(%%rdx),%%xmm4		/*~t4 = t4*sx3 */\n\t"\
			"mulpd	    (%%rbx),%%xmm1		/* t1 = t1*cx1 */\n\t"\
			"mulpd	0x10(%%rbx),%%xmm6		/*~t6 = t6*sx1 */\n\t"\
			"mulpd	    (%%rax),%%xmm0		/* ~r = r*(cx0-1) */\n\t"\
			"mulpd	0x10(%%rcx),%%xmm7		/*~t5 = t5*sx2 */\n\t"\
			"mulpd	    (%%rdx),%%xmm2		/* t3 =(t1+t2)*cx3 */\n\t"\
			"addpd	0x190(%%rdi),%%xmm0		/* t0 =~r + B0 */\n\t"\
			"addpd	%%xmm4,%%xmm6			/*~t6 = t4+t6 */\n\t"\
			"subpd	%%xmm2,%%xmm1			/* tt = t1-t3 */\n\t"\
			"subpd	%%xmm7,%%xmm4			/*~t5 = t4-t5, xmm7 FREE */\n\t"\
			"subpd	%%xmm2,%%xmm3			/* t2 = t2-t3, xmm2 FREE */\n\t"\
			"movaps	%%xmm0,%%xmm2			/* cpy t0 */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* cpy t3 */\n\t"\
			"addpd	%%xmm1,%%xmm0			/*~t0 = t0+tt */\n\t"\
			"addpd	%%xmm6,%%xmm5			/*~t3 = t3+t6 */\n\t"\
			"addpd	%%xmm3,%%xmm1			/*~tt = tt+t2 */\n\t"\
			"addpd	%%xmm4,%%xmm6			/*      t6+t5 */\n\t"\
			"addpd	%%xmm2,%%xmm3			/*~t2 = t2+t0 */\n\t"\
			"addpd	%%xmm7,%%xmm4			/*~t5 = t5+t3 */\n\t"\
			"subpd	%%xmm1,%%xmm2			/*~t1 = t0-tt-t2, xmm1 FREE */\n\t"\
			"subpd	%%xmm6,%%xmm7			/*~t4 = t3-t6-t5, xmm6 FREE */\n\t"\
			"\n\t"\
			"/* xmm1,6 FREE */\n\t"\
			"movaps	0x280(%%rdi),%%xmm1		/* t0r */\n\t"\
			"movaps	0x080(%%rdi),%%xmm6		/* t3r */\n\t"\
			"subpd	%%xmm5,%%xmm1			/* B1r =t0r-t3i */\n\t"\
			"subpd	%%xmm6,%%xmm0			/* B6i =t0i-t3r */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*        2*t3i */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*        2*t3r */\n\t"\
			"addpd	%%xmm1,%%xmm5			/* B6r =t0r+t3i */\n\t"\
			"addpd	%%xmm0,%%xmm6			/* B1i =t0i+t3r */\n\t"\
			"movaps	%%xmm1,0x280(%%rdi)		/* <-B1r = s1p08r */\n\t"\
			"movaps	%%xmm0,0x090(%%rdi)		/* <-B6i = s1p20r */\n\t"\
			"movaps	%%xmm5,0x080(%%rdi)		/* <-B6r = s1p20r */\n\t"\
			"movaps	%%xmm6,0x290(%%rdi)		/* <-B1i = s1p08i */\n\t"\
			"\n\t"\
			"movaps	     (%%rdi),%%xmm1		/* t1r */\n\t"\
			"movaps	0x300(%%rdi),%%xmm6		/* t4r */\n\t"\
			"subpd	%%xmm7,%%xmm1			/* B2r =t1r-t4i */\n\t"\
			"subpd	%%xmm6,%%xmm2			/* B5i =t1i-t4r */\n\t"\
			"addpd	%%xmm7,%%xmm7			/*        2*t4i */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*        2*t4r */\n\t"\
			"addpd	%%xmm1,%%xmm7			/* B5r =t1r+t4i */\n\t"\
			"addpd	%%xmm2,%%xmm6			/* B2i =t1i+t4r */\n\t"\
			"movaps	%%xmm1,     (%%rdi)		/* <-B2r = s1p16r */\n\t"\
			"movaps	%%xmm2,0x310(%%rdi)		/* <-B5i = s1p12r */\n\t"\
			"movaps	%%xmm7,0x300(%%rdi)		/* <-B5r = s1p12r */\n\t"\
			"movaps	%%xmm6,0x010(%%rdi)		/* <-B2i = s1p16i */\n\t"\
			"\n\t"\
			"/* Note the order reversal on this pair of outputs: */\n\t"\
			"movaps	0x100(%%rdi),%%xmm0		/* t2r */\n\t"\
			"movaps	0x200(%%rdi),%%xmm5		/* t5r */\n\t"\
			"subpd	%%xmm4,%%xmm0			/* B4r =t2r-t5i */\n\t"\
			"subpd	%%xmm5,%%xmm3			/* B3i =t2i-t5r */\n\t"\
			"addpd	%%xmm4,%%xmm4			/*        2*t5i */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*        2*t5r */\n\t"\
			"addpd	%%xmm0,%%xmm4			/* B3r =t2r+t5i */\n\t"\
			"addpd	%%xmm3,%%xmm5			/* B4i =t2i+t5r */\n\t"\
			"movaps	%%xmm0,0x200(%%rdi)		/* <-B4r = s1p04r */\n\t"\
			"movaps	%%xmm3,0x110(%%rdi)		/* <-B3i = s1p24r */\n\t"\
			"movaps	%%xmm4,0x100(%%rdi)		/* <-B3r = s1p24r */\n\t"\
			"movaps	%%xmm5,0x210(%%rdi)		/* <-B4i = s1p04i */\n\t"\
			"\n\t"\
		"/* SSE2_RADIX_07_DFT(s1p03r,s1p07r,s1p11r,s1p15r,s1p19r,s1p23r,s1p27r,cc0,s1p07r,s1p15r,s1p23r,s1p03r,s1p11r,s1p19r,s1p27r): */\n\t"\
			"addq	$0x20,%%rdi	/* s1p03r */\n\t"\
			"movaps	0x080(%%rdi),%%xmm6			/* A1r = s1p07r */\n\t"\
			"movaps	0x300(%%rdi),%%xmm1			/* A6r = s1p27r */\n\t"\
			"movaps	0x100(%%rdi),%%xmm5			/* A2r = s1p11r */\n\t"\
			"movaps	0x280(%%rdi),%%xmm2			/* A5r = s1p23r */\n\t"\
			"movaps	0x180(%%rdi),%%xmm4			/* A3r = s1p15r */\n\t"\
			"movaps	0x200(%%rdi),%%xmm3			/* A4r = s1p19r */\n\t"\
			"\n\t"\
			"subpd	%%xmm1,%%xmm6			/* t6r = A1r-A6r */\n\t"\
			"addpd	%%xmm1,%%xmm1			/*         2*A6r */\n\t"\
			"addpd	%%xmm6,%%xmm1			/* t1r = A1r+A6r */\n\t"\
			"\n\t"\
			"subpd	%%xmm2,%%xmm5			/* t5r = A2r-A5r */\n\t"\
			"addpd	%%xmm2,%%xmm2			/*         2*A5r */\n\t"\
			"addpd	%%xmm5,%%xmm2			/* t2r = A2r+A5r */\n\t"\
			"\n\t"\
			"movaps	    (%%rdi),%%xmm0		/* Ar0 */\n\t"\
			"subpd	%%xmm3,%%xmm4			/* t4r = A3r-A4r */\n\t"\
			"addpd	%%xmm3,%%xmm3			/*         2*A4r */\n\t"\
			"addpd	%%xmm4,%%xmm3			/* t3r = A3r+A4r */\n\t"\
			"\n\t"\
			"movaps	%%xmm0,    (%%rsi)		/* cpy t0 into scratch sincos slot */\n\t"\
			"movaps	%%xmm6,0x10(%%rsi)		/* cpy t6 into scratch sincos slot */\n\t"\
			"addpd	%%xmm1,%%xmm0			/*~A0 = A0+t1 */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* cpy t5 */\n\t"\
			"addpd	%%xmm2,%%xmm3			/*~t3 = t3+t2 */\n\t"\
			"subpd	%%xmm4,%%xmm5			/*~t5 = t5-t4 */\n\t"\
			"subpd	%%xmm2,%%xmm1			/*~t1 = t1-t2 */\n\t"\
			"subpd	%%xmm7,%%xmm6			/*~t6 = t6-t5 */\n\t"\
			"addpd	%%xmm2,%%xmm2			/* 2*t2 */\n\t"\
			"addpd	%%xmm7,%%xmm4			/*~t5 = t4+t5 */\n\t"\
			"addpd	%%xmm3,%%xmm0			/* B0 */\n\t"\
			"addpd	0x10(%%rsi),%%xmm5		/* t3 = (t5-t4)+t6 */\n\t"\
			"subpd	%%xmm2,%%xmm3			/*~t2 =  (t2+t3) - 2*t2 = t3-t2 */\n\t"\
			"movaps	%%xmm4,%%xmm7			/* cpy t5 */\n\t"\
			"movaps	%%xmm0,0x080(%%rdi)		/* <-B0 = s1p07r, xmm0 FREE */\n\t"\
			"subpd	%%xmm6,%%xmm4			/* t4 = ~t5-~t6 */\n\t"\
			"movaps	%%xmm1,%%xmm2			/* cpy ~t1 */\n\t"\
			"subpd	    (%%rsi),%%xmm0		/* r = B0 - t0 */\n\t"\
			"mulpd	0x10(%%rax),%%xmm5		/*~t3 = t3*sx0 */\n\t"\
			"addpd	%%xmm3,%%xmm2			/* ~t1+~t2 */\n\t"\
			"mulpd	    (%%rcx),%%xmm3		/* t2 = t2*cx2 */\n\t"\
			"mulpd	0x10(%%rdx),%%xmm4		/*~t4 = t4*sx3 */\n\t"\
			"mulpd	    (%%rbx),%%xmm1		/* t1 = t1*cx1 */\n\t"\
			"mulpd	0x10(%%rbx),%%xmm6		/*~t6 = t6*sx1 */\n\t"\
			"mulpd	    (%%rax),%%xmm0		/* ~r = r*(cx0-1) */\n\t"\
			"mulpd	0x10(%%rcx),%%xmm7		/*~t5 = t5*sx2 */\n\t"\
			"mulpd	    (%%rdx),%%xmm2		/* t3 =(t1+t2)*cx3 */\n\t"\
			"addpd	0x080(%%rdi),%%xmm0		/* t0 =~r + B0 */\n\t"\
			"addpd	%%xmm4,%%xmm6			/*~t6 = t4+t6 */\n\t"\
			"subpd	%%xmm2,%%xmm1			/* tt = t1-t3 */\n\t"\
			"subpd	%%xmm7,%%xmm4			/*~t5 = t4-t5, xmm7 FREE */\n\t"\
			"subpd	%%xmm2,%%xmm3			/* t2 = t2-t3, xmm2 FREE */\n\t"\
			"\n\t"\
			"movaps	%%xmm0,%%xmm2			/* cpy t0 */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* cpy t3 */\n\t"\
			"addpd	%%xmm1,%%xmm0			/*~t0 = t0+tt */\n\t"\
			"addpd	%%xmm6,%%xmm5			/*~t3 = t3+t6 */\n\t"\
			"addpd	%%xmm3,%%xmm1			/*~tt = tt+t2 */\n\t"\
			"addpd	%%xmm4,%%xmm6			/*      t6+t5 */\n\t"\
			"addpd	%%xmm2,%%xmm3			/*~t2 = t2+t0 */\n\t"\
			"addpd	%%xmm7,%%xmm4			/*~t5 = t5+t3 */\n\t"\
			"subpd	%%xmm1,%%xmm2			/*~t1 = t0-tt-t2 */\n\t"\
			"subpd	%%xmm6,%%xmm7			/*~t4 = t3-t6-t5 */\n\t"\
			"movaps	%%xmm0,0x180(%%rdi)			/* B1 <- t0 = s1p15r */\n\t"\
			"movaps	%%xmm5,0x300(%%rdi)			/* B6 <- t3 = s1p27r */\n\t"\
			"movaps	%%xmm2,0x280(%%rdi)			/* B2 <- t1 = s1p23r */\n\t"\
			"movaps	%%xmm7,0x200(%%rdi)			/* B5 <- t4 = s1p19r */\n\t"\
			"movaps	%%xmm3,     (%%rdi)			/* B3 <- t2 = s1p03r */\n\t"\
			"movaps	%%xmm4,0x100(%%rdi)			/* B4 <- t5 = s1p11r */\n\t"\
			"/************************** Imaginary Parts: ******************************************/\n\t"\
			"movaps	0x090(%%rdi),%%xmm6		/* A1i = s1p07i */\n\t"\
			"movaps	0x310(%%rdi),%%xmm1		/* A6i = s1p27i */\n\t"\
			"movaps	0x110(%%rdi),%%xmm5		/* A2i = s1p11i */\n\t"\
			"movaps	0x290(%%rdi),%%xmm2		/* A5i = s1p23i */\n\t"\
			"movaps	0x190(%%rdi),%%xmm4		/* A3i = s1p15i */\n\t"\
			"movaps	0x210(%%rdi),%%xmm3		/* A4i = s1p19i */\n\t"\
			"\n\t"\
			"subpd	%%xmm1,%%xmm6			/* t6i = A1i-A6i */\n\t"\
			"addpd	%%xmm1,%%xmm1			/*         2*A6i */\n\t"\
			"addpd	%%xmm6,%%xmm1			/* t1i = A1i+A6i */\n\t"\
			"\n\t"\
			"subpd	%%xmm2,%%xmm5			/* t5i = A2i-A5i */\n\t"\
			"addpd	%%xmm2,%%xmm2			/*         2*A5i */\n\t"\
			"addpd	%%xmm5,%%xmm2			/* t2i = A2i+A5i */\n\t"\
			"\n\t"\
			"movaps	0x10(%%rdi),%%xmm0		/* Ai0 */\n\t"\
			"subpd	%%xmm3,%%xmm4			/* t4i = A3i-A4i */\n\t"\
			"addpd	%%xmm3,%%xmm3			/*         2*A4i */\n\t"\
			"addpd	%%xmm4,%%xmm3			/* t3i = A3i+A4i */\n\t"\
			"\n\t"\
			"movaps	%%xmm0,    (%%rsi)		/* cpy t0 into scratch sincos slot */\n\t"\
			"movaps	%%xmm6,0x10(%%rsi)		/* cpy t6 into scratch sincos slot */\n\t"\
			"addpd	%%xmm1,%%xmm0			/*~A0 = A0+t1 */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* cpy t5 */\n\t"\
			"addpd	%%xmm2,%%xmm3			/*~t3 = t3+t2 */\n\t"\
			"subpd	%%xmm4,%%xmm5			/*~t5 = t5-t4 */\n\t"\
			"subpd	%%xmm2,%%xmm1			/*~t1 = t1-t2 */\n\t"\
			"subpd	%%xmm7,%%xmm6			/*~t6 = t6-t5 */\n\t"\
			"addpd	%%xmm2,%%xmm2			/* 2*t2 */\n\t"\
			"addpd	%%xmm7,%%xmm4			/*~t5 = t4+t5 */\n\t"\
			"addpd	%%xmm3,%%xmm0			/* B0 */\n\t"\
			"addpd	0x10(%%rsi),%%xmm5		/* t3 = (t5-t4)+t6 */\n\t"\
			"subpd	%%xmm2,%%xmm3			/*~t2 =  (t2+t3) - 2*t2 = t3-t2 */\n\t"\
			"movaps	%%xmm4,%%xmm7			/* cpy t5 */\n\t"\
			"movaps	%%xmm0,0x090(%%rdi)		/* <-B0 = s1p07i, xmm0 FREE */\n\t"\
			"subpd	%%xmm6,%%xmm4			/* t4 = ~t5-~t6 */\n\t"\
			"movaps	%%xmm1,%%xmm2			/* cpy ~t1 */\n\t"\
			"subpd	    (%%rsi),%%xmm0		/* r = B0 - t0 */\n\t"\
			"mulpd	0x10(%%rax),%%xmm5		/*~t3 = t3*sx0 */\n\t"\
			"addpd	%%xmm3,%%xmm2			/* ~t1+~t2 */\n\t"\
			"mulpd	    (%%rcx),%%xmm3		/* t2 = t2*cx2 */\n\t"\
			"mulpd	0x10(%%rdx),%%xmm4		/*~t4 = t4*sx3 */\n\t"\
			"mulpd	    (%%rbx),%%xmm1		/* t1 = t1*cx1 */\n\t"\
			"mulpd	0x10(%%rbx),%%xmm6		/*~t6 = t6*sx1 */\n\t"\
			"mulpd	    (%%rax),%%xmm0		/* ~r = r*(cx0-1) */\n\t"\
			"mulpd	0x10(%%rcx),%%xmm7		/*~t5 = t5*sx2 */\n\t"\
			"mulpd	    (%%rdx),%%xmm2		/* t3 =(t1+t2)*cx3 */\n\t"\
			"addpd	0x090(%%rdi),%%xmm0		/* t0 =~r + B0 */\n\t"\
			"addpd	%%xmm4,%%xmm6			/*~t6 = t4+t6 */\n\t"\
			"subpd	%%xmm2,%%xmm1			/* tt = t1-t3 */\n\t"\
			"subpd	%%xmm7,%%xmm4			/*~t5 = t4-t5, xmm7 FREE */\n\t"\
			"subpd	%%xmm2,%%xmm3			/* t2 = t2-t3, xmm2 FREE */\n\t"\
			"movaps	%%xmm0,%%xmm2			/* cpy t0 */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* cpy t3 */\n\t"\
			"addpd	%%xmm1,%%xmm0			/*~t0 = t0+tt */\n\t"\
			"addpd	%%xmm6,%%xmm5			/*~t3 = t3+t6 */\n\t"\
			"addpd	%%xmm3,%%xmm1			/*~tt = tt+t2 */\n\t"\
			"addpd	%%xmm4,%%xmm6			/*      t6+t5 */\n\t"\
			"addpd	%%xmm2,%%xmm3			/*~t2 = t2+t0 */\n\t"\
			"addpd	%%xmm7,%%xmm4			/*~t5 = t5+t3 */\n\t"\
			"subpd	%%xmm1,%%xmm2			/*~t1 = t0-tt-t2, xmm1 FREE */\n\t"\
			"subpd	%%xmm6,%%xmm7			/*~t4 = t3-t6-t5, xmm6 FREE */\n\t"\
			"\n\t"\
			"/* xmm1,6 FREE */\n\t"\
			"movaps	0x180(%%rdi),%%xmm1		/* t0r */\n\t"\
			"movaps	0x300(%%rdi),%%xmm6		/* t3r */\n\t"\
			"subpd	%%xmm5,%%xmm1			/* B1r =t0r-t3i */\n\t"\
			"subpd	%%xmm6,%%xmm0			/* B6i =t0i-t3r */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*        2*t3i */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*        2*t3r */\n\t"\
			"addpd	%%xmm1,%%xmm5			/* B6r =t0r+t3i */\n\t"\
			"addpd	%%xmm0,%%xmm6			/* B1i =t0i+t3r */\n\t"\
			"movaps	%%xmm1,0x180(%%rdi)		/* <-B1r = s1p08r */\n\t"\
			"movaps	%%xmm0,0x310(%%rdi)		/* <-B6i = s1p20r */\n\t"\
			"movaps	%%xmm5,0x300(%%rdi)		/* <-B6r = s1p20r */\n\t"\
			"movaps	%%xmm6,0x190(%%rdi)		/* <-B1i = s1p08i */\n\t"\
			"\n\t"\
			"movaps	0x280(%%rdi),%%xmm1		/* t1r */\n\t"\
			"movaps	0x200(%%rdi),%%xmm6		/* t4r */\n\t"\
			"subpd	%%xmm7,%%xmm1			/* B2r =t1r-t4i */\n\t"\
			"subpd	%%xmm6,%%xmm2			/* B5i =t1i-t4r */\n\t"\
			"addpd	%%xmm7,%%xmm7			/*        2*t4i */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*        2*t4r */\n\t"\
			"addpd	%%xmm1,%%xmm7			/* B5r =t1r+t4i */\n\t"\
			"addpd	%%xmm2,%%xmm6			/* B2i =t1i+t4r */\n\t"\
			"movaps	%%xmm1,0x280(%%rdi)		/* <-B2r = s1p16r */\n\t"\
			"movaps	%%xmm2,0x210(%%rdi)		/* <-B5i = s1p12r */\n\t"\
			"movaps	%%xmm7,0x200(%%rdi)		/* <-B5r = s1p12r */\n\t"\
			"movaps	%%xmm6,0x290(%%rdi)		/* <-B2i = s1p16i */\n\t"\
			"\n\t"\
			"/* Note the order reversal on this pair of outputs: */\n\t"\
			"movaps	     (%%rdi),%%xmm0		/* t2r */\n\t"\
			"movaps	0x100(%%rdi),%%xmm5		/* t5r */\n\t"\
			"subpd	%%xmm4,%%xmm0			/* B4r =t2r-t5i */\n\t"\
			"subpd	%%xmm5,%%xmm3			/* B3i =t2i-t5r */\n\t"\
			"addpd	%%xmm4,%%xmm4			/*        2*t5i */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*        2*t5r */\n\t"\
			"addpd	%%xmm0,%%xmm4			/* B3r =t2r+t5i */\n\t"\
			"addpd	%%xmm3,%%xmm5			/* B4i =t2i+t5r */\n\t"\
			"movaps	%%xmm0,0x100(%%rdi)		/* <-B4r = s1p04r */\n\t"\
			"movaps	%%xmm3,0x010(%%rdi)		/* <-B3i = s1p24r */\n\t"\
			"movaps	%%xmm4,     (%%rdi)		/* <-B3r = s1p24r */\n\t"\
			"movaps	%%xmm5,0x110(%%rdi)		/* <-B4i = s1p04i */\n\t"\
			"\n\t"\
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
			: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
		);\
	}

	#define	SSE2_RADIX28_DIF_NOTWIDDLE(Xadd,Xp01,Xp02,Xp03,Xp04,Xp08,Xp12,Xp16,Xp20,Xp24,Xout,Xcc0)\
	{\
	__asm__ volatile (\
		"/* SSE2_RADIX_07_DFT(s1p00r,s1p24r,s1p20r,s1p16r,s1p12r,s1p08r,s1p04r,cc0,s1p00r,s1p04r,s1p08r,s1p12r,s1p16r,s1p20r,s1p24r): */\n\t"\
			"movq	%[__out],%%rdi	/* s1p00r */\n\t"\
			"movq	%[__cc0],%%rax	/* cc0 */\n\t"\
			"movq	$0x20,%%rbx\n\t"\
			"movq	$0x40,%%rcx\n\t"\
			"movq	$0x60,%%rdx\n\t"\
			"movq	$0x80,%%rsi\n\t"\
			"addq	%%rax,%%rbx		/* cc1 */\n\t"\
			"addq	%%rax,%%rcx		/* cc2 */\n\t"\
			"addq	%%rax,%%rdx		/* cc3 */\n\t"\
			"addq	%%rax,%%rsi		/* scratch storage */\n\t"\
			"movaps	0x300(%%rdi),%%xmm6			/* A1r = s1p24r */\n\t"\
			"movaps	0x080(%%rdi),%%xmm1			/* A6r = s1p04r */\n\t"\
			"movaps	0x280(%%rdi),%%xmm5			/* A2r = s1p20r */\n\t"\
			"movaps	0x100(%%rdi),%%xmm2			/* A5r = s1p08r */\n\t"\
			"movaps	0x200(%%rdi),%%xmm4			/* A3r = s1p16r */\n\t"\
			"movaps	0x180(%%rdi),%%xmm3			/* A4r = s1p12r */\n\t"\
			"\n\t"\
			"subpd	%%xmm1,%%xmm6			/* t6r = A1r-A6r */\n\t"\
			"addpd	%%xmm1,%%xmm1			/*         2*A6r */\n\t"\
			"addpd	%%xmm6,%%xmm1			/* t1r = A1r+A6r */\n\t"\
			"\n\t"\
			"subpd	%%xmm2,%%xmm5			/* t5r = A2r-A5r */\n\t"\
			"addpd	%%xmm2,%%xmm2			/*         2*A5r */\n\t"\
			"addpd	%%xmm5,%%xmm2			/* t2r = A2r+A5r */\n\t"\
			"\n\t"\
			"movaps	    (%%rdi),%%xmm0		/* Ar0 */\n\t"\
			"subpd	%%xmm3,%%xmm4			/* t4r = A3r-A4r */\n\t"\
			"addpd	%%xmm3,%%xmm3			/*         2*A4r */\n\t"\
			"addpd	%%xmm4,%%xmm3			/* t3r = A3r+A4r */\n\t"\
			"\n\t"\
			"movaps	%%xmm0,    (%%rsi)		/* cpy t0 into scratch sincos slot */\n\t"\
			"movaps	%%xmm6,0x10(%%rsi)		/* cpy t6 into scratch sincos slot */\n\t"\
			"addpd	%%xmm1,%%xmm0			/*~A0 = A0+t1 */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* cpy t5 */\n\t"\
			"addpd	%%xmm2,%%xmm3			/*~t3 = t3+t2 */\n\t"\
			"subpd	%%xmm4,%%xmm5			/*~t5 = t5-t4 */\n\t"\
			"subpd	%%xmm2,%%xmm1			/*~t1 = t1-t2 */\n\t"\
			"subpd	%%xmm7,%%xmm6			/*~t6 = t6-t5 */\n\t"\
			"addpd	%%xmm2,%%xmm2			/* 2*t2 */\n\t"\
			"addpd	%%xmm7,%%xmm4			/*~t5 = t4+t5 */\n\t"\
			"addpd	%%xmm3,%%xmm0			/* B0 */\n\t"\
			"addpd	0x10(%%rsi),%%xmm5		/* t3 = (t5-t4)+t6 */\n\t"\
			"subpd	%%xmm2,%%xmm3			/*~t2 =  (t2+t3) - 2*t2 = t3-t2 */\n\t"\
			"movaps	%%xmm4,%%xmm7			/* cpy t5 */\n\t"\
			"movaps	%%xmm0,     (%%rdi)		/* <-B0 = s1p00r, xmm0 FREE */\n\t"\
			"subpd	%%xmm6,%%xmm4			/* t4 = ~t5-~t6 */\n\t"\
			"movaps	%%xmm1,%%xmm2			/* cpy ~t1 */\n\t"\
			"subpd	    (%%rsi),%%xmm0		/* r = B0 - t0 */\n\t"\
			"mulpd	0x10(%%rax),%%xmm5		/*~t3 = t3*sx0 */\n\t"\
			"addpd	%%xmm3,%%xmm2			/* ~t1+~t2 */\n\t"\
			"mulpd	    (%%rcx),%%xmm3		/* t2 = t2*cx2 */\n\t"\
			"mulpd	0x10(%%rdx),%%xmm4		/*~t4 = t4*sx3 */\n\t"\
			"mulpd	    (%%rbx),%%xmm1		/* t1 = t1*cx1 */\n\t"\
			"mulpd	0x10(%%rbx),%%xmm6		/*~t6 = t6*sx1 */\n\t"\
			"mulpd	    (%%rax),%%xmm0		/* ~r = r*(cx0-1) */\n\t"\
			"mulpd	0x10(%%rcx),%%xmm7		/*~t5 = t5*sx2 */\n\t"\
			"mulpd	    (%%rdx),%%xmm2		/* t3 =(t1+t2)*cx3 */\n\t"\
			"addpd	    (%%rdi),%%xmm0		/* t0 =~r + B0 */\n\t"\
			"addpd	%%xmm4,%%xmm6			/*~t6 = t4+t6 */\n\t"\
			"subpd	%%xmm2,%%xmm1			/* tt = t1-t3 */\n\t"\
			"subpd	%%xmm7,%%xmm4			/*~t5 = t4-t5, xmm7 FREE */\n\t"\
			"subpd	%%xmm2,%%xmm3			/* t2 = t2-t3, xmm2 FREE */\n\t"\
			"\n\t"\
			"movaps	%%xmm0,%%xmm2			/* cpy t0 */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* cpy t3 */\n\t"\
			"addpd	%%xmm1,%%xmm0			/*~t0 = t0+tt */\n\t"\
			"addpd	%%xmm6,%%xmm5			/*~t3 = t3+t6 */\n\t"\
			"addpd	%%xmm3,%%xmm1			/*~tt = tt+t2 */\n\t"\
			"addpd	%%xmm4,%%xmm6			/*      t6+t5 */\n\t"\
			"addpd	%%xmm2,%%xmm3			/*~t2 = t2+t0 */\n\t"\
			"addpd	%%xmm7,%%xmm4			/*~t5 = t5+t3 */\n\t"\
			"subpd	%%xmm1,%%xmm2			/*~t1 = t0-tt-t2 */\n\t"\
			"subpd	%%xmm6,%%xmm7			/*~t4 = t3-t6-t5 */\n\t"\
			"movaps	%%xmm0,0x080(%%rdi)			/* B1 <- t0 = s1p04r */\n\t"\
			"movaps	%%xmm5,0x300(%%rdi)			/* B6 <- t3 = s1p24r */\n\t"\
			"movaps	%%xmm2,0x100(%%rdi)			/* B2 <- t1 = s1p08r */\n\t"\
			"movaps	%%xmm7,0x280(%%rdi)			/* B5 <- t4 = s1p20r */\n\t"\
			"movaps	%%xmm3,0x180(%%rdi)			/* B3 <- t2 = s1p12r */\n\t"\
			"movaps	%%xmm4,0x200(%%rdi)			/* B4 <- t5 = s1p16r */\n\t"\
			"/************************** Imaginary Parts: ******************************************/\n\t"\
			"movaps	0x310(%%rdi),%%xmm6		/* A1i = s1p24i */\n\t"\
			"movaps	0x090(%%rdi),%%xmm1		/* A6i = s1p04i */\n\t"\
			"movaps	0x290(%%rdi),%%xmm5		/* A2i = s1p20i */\n\t"\
			"movaps	0x110(%%rdi),%%xmm2		/* A5i = s1p08i */\n\t"\
			"movaps	0x210(%%rdi),%%xmm4		/* A3i = s1p16i */\n\t"\
			"movaps	0x190(%%rdi),%%xmm3		/* A4i = s1p12i */\n\t"\
			"\n\t"\
			"subpd	%%xmm1,%%xmm6			/* t6i = A1i-A6i */\n\t"\
			"addpd	%%xmm1,%%xmm1			/*         2*A6i */\n\t"\
			"addpd	%%xmm6,%%xmm1			/* t1i = A1i+A6i */\n\t"\
			"\n\t"\
			"subpd	%%xmm2,%%xmm5			/* t5i = A2i-A5i */\n\t"\
			"addpd	%%xmm2,%%xmm2			/*         2*A5i */\n\t"\
			"addpd	%%xmm5,%%xmm2			/* t2i = A2i+A5i */\n\t"\
			"\n\t"\
			"movaps	0x10(%%rdi),%%xmm0		/* Ai0 */\n\t"\
			"subpd	%%xmm3,%%xmm4			/* t4i = A3i-A4i */\n\t"\
			"addpd	%%xmm3,%%xmm3			/*         2*A4i */\n\t"\
			"addpd	%%xmm4,%%xmm3			/* t3i = A3i+A4i */\n\t"\
			"\n\t"\
			"movaps	%%xmm0,    (%%rsi)		/* cpy t0 into scratch sincos slot */\n\t"\
			"movaps	%%xmm6,0x10(%%rsi)		/* cpy t6 into scratch sincos slot */\n\t"\
			"addpd	%%xmm1,%%xmm0			/*~A0 = A0+t1 */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* cpy t5 */\n\t"\
			"addpd	%%xmm2,%%xmm3			/*~t3 = t3+t2 */\n\t"\
			"subpd	%%xmm4,%%xmm5			/*~t5 = t5-t4 */\n\t"\
			"subpd	%%xmm2,%%xmm1			/*~t1 = t1-t2 */\n\t"\
			"subpd	%%xmm7,%%xmm6			/*~t6 = t6-t5 */\n\t"\
			"addpd	%%xmm2,%%xmm2			/* 2*t2 */\n\t"\
			"addpd	%%xmm7,%%xmm4			/*~t5 = t4+t5 */\n\t"\
			"addpd	%%xmm3,%%xmm0			/* B0 */\n\t"\
			"addpd	0x10(%%rsi),%%xmm5		/* t3 = (t5-t4)+t6 */\n\t"\
			"subpd	%%xmm2,%%xmm3			/*~t2 =  (t2+t3) - 2*t2 = t3-t2 */\n\t"\
			"movaps	%%xmm4,%%xmm7			/* cpy t5 */\n\t"\
			"movaps	%%xmm0,0x010(%%rdi)		/* <-B0 = s1p00i, xmm0 FREE */\n\t"\
			"subpd	%%xmm6,%%xmm4			/* t4 = ~t5-~t6 */\n\t"\
			"movaps	%%xmm1,%%xmm2			/* cpy ~t1 */\n\t"\
			"subpd	    (%%rsi),%%xmm0		/* r = B0 - t0 */\n\t"\
			"mulpd	0x10(%%rax),%%xmm5		/*~t3 = t3*sx0 */\n\t"\
			"addpd	%%xmm3,%%xmm2			/* ~t1+~t2 */\n\t"\
			"mulpd	    (%%rcx),%%xmm3		/* t2 = t2*cx2 */\n\t"\
			"mulpd	0x10(%%rdx),%%xmm4		/*~t4 = t4*sx3 */\n\t"\
			"mulpd	    (%%rbx),%%xmm1		/* t1 = t1*cx1 */\n\t"\
			"mulpd	0x10(%%rbx),%%xmm6		/*~t6 = t6*sx1 */\n\t"\
			"mulpd	    (%%rax),%%xmm0		/* ~r = r*(cx0-1) */\n\t"\
			"mulpd	0x10(%%rcx),%%xmm7		/*~t5 = t5*sx2 */\n\t"\
			"mulpd	    (%%rdx),%%xmm2		/* t3 =(t1+t2)*cx3 */\n\t"\
			"addpd	0x10(%%rdi),%%xmm0		/* t0 =~r + B0 */\n\t"\
			"addpd	%%xmm4,%%xmm6			/*~t6 = t4+t6 */\n\t"\
			"subpd	%%xmm2,%%xmm1			/* tt = t1-t3 */\n\t"\
			"subpd	%%xmm7,%%xmm4			/*~t5 = t4-t5, xmm7 FREE */\n\t"\
			"subpd	%%xmm2,%%xmm3			/* t2 = t2-t3, xmm2 FREE */\n\t"\
			"movaps	%%xmm0,%%xmm2			/* cpy t0 */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* cpy t3 */\n\t"\
			"addpd	%%xmm1,%%xmm0			/*~t0 = t0+tt */\n\t"\
			"addpd	%%xmm6,%%xmm5			/*~t3 = t3+t6 */\n\t"\
			"addpd	%%xmm3,%%xmm1			/*~tt = tt+t2 */\n\t"\
			"addpd	%%xmm4,%%xmm6			/*      t6+t5 */\n\t"\
			"addpd	%%xmm2,%%xmm3			/*~t2 = t2+t0 */\n\t"\
			"addpd	%%xmm7,%%xmm4			/*~t5 = t5+t3 */\n\t"\
			"subpd	%%xmm1,%%xmm2			/*~t1 = t0-tt-t2, xmm1 FREE */\n\t"\
			"subpd	%%xmm6,%%xmm7			/*~t4 = t3-t6-t5, xmm6 FREE */\n\t"\
			"\n\t"\
			"/* xmm1,6 FREE */\n\t"\
			"movaps	0x080(%%rdi),%%xmm1		/* t0r */\n\t"\
			"movaps	0x300(%%rdi),%%xmm6		/* t3r */\n\t"\
			"subpd	%%xmm5,%%xmm1			/* B1r =t0r-t3i */\n\t"\
			"subpd	%%xmm6,%%xmm0			/* B6i =t0i-t3r */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*        2*t3i */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*        2*t3r */\n\t"\
			"addpd	%%xmm1,%%xmm5			/* B6r =t0r+t3i */\n\t"\
			"addpd	%%xmm0,%%xmm6			/* B1i =t0i+t3r */\n\t"\
			"movaps	%%xmm1,0x080(%%rdi)		/* <-B1r = s1p08r */\n\t"\
			"movaps	%%xmm0,0x310(%%rdi)		/* <-B6i = s1p20r */\n\t"\
			"movaps	%%xmm5,0x300(%%rdi)		/* <-B6r = s1p20r */\n\t"\
			"movaps	%%xmm6,0x090(%%rdi)		/* <-B1i = s1p08i */\n\t"\
			"\n\t"\
			"movaps	0x100(%%rdi),%%xmm1		/* t1r */\n\t"\
			"movaps	0x280(%%rdi),%%xmm6		/* t4r */\n\t"\
			"subpd	%%xmm7,%%xmm1			/* B2r =t1r-t4i */\n\t"\
			"subpd	%%xmm6,%%xmm2			/* B5i =t1i-t4r */\n\t"\
			"addpd	%%xmm7,%%xmm7			/*        2*t4i */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*        2*t4r */\n\t"\
			"addpd	%%xmm1,%%xmm7			/* B5r =t1r+t4i */\n\t"\
			"addpd	%%xmm2,%%xmm6			/* B2i =t1i+t4r */\n\t"\
			"movaps	%%xmm1,0x100(%%rdi)		/* <-B2r = s1p16r */\n\t"\
			"movaps	%%xmm2,0x290(%%rdi)		/* <-B5i = s1p12r */\n\t"\
			"movaps	%%xmm7,0x280(%%rdi)		/* <-B5r = s1p12r */\n\t"\
			"movaps	%%xmm6,0x110(%%rdi)		/* <-B2i = s1p16i */\n\t"\
			"\n\t"\
			"/* Note the order reversal on this pair of outputs: */\n\t"\
			"movaps	0x180(%%rdi),%%xmm0		/* t2r */\n\t"\
			"movaps	0x200(%%rdi),%%xmm5		/* t5r */\n\t"\
			"subpd	%%xmm4,%%xmm0			/* B4r =t2r-t5i */\n\t"\
			"subpd	%%xmm5,%%xmm3			/* B3i =t2i-t5r */\n\t"\
			"addpd	%%xmm4,%%xmm4			/*        2*t5i */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*        2*t5r */\n\t"\
			"addpd	%%xmm0,%%xmm4			/* B3r =t2r+t5i */\n\t"\
			"addpd	%%xmm3,%%xmm5			/* B4i =t2i+t5r */\n\t"\
			"movaps	%%xmm0,0x200(%%rdi)		/* <-B4r = s1p04r */\n\t"\
			"movaps	%%xmm3,0x190(%%rdi)		/* <-B3i = s1p24r */\n\t"\
			"movaps	%%xmm4,0x180(%%rdi)		/* <-B3r = s1p24r */\n\t"\
			"movaps	%%xmm5,0x210(%%rdi)		/* <-B4i = s1p04i */\n\t"\
			"\n\t"\
		"/* SSE2_RADIX_07_DFT(s1p21r,s1p17r,s1p13r,s1p09r,s1p05r,s1p01r,s1p25r,cc0,s1p01r,s1p05r,s1p09r,s1p13r,s1p17r,s1p21r,s1p25r): */\n\t"\
			"addq	$0x20,%%rdi	/* s1p01r */\n\t"\
			"movaps	0x200(%%rdi),%%xmm6			/* A1r = s1p17r */\n\t"\
			"movaps	0x300(%%rdi),%%xmm1			/* A6r = s1p25r */\n\t"\
			"movaps	0x180(%%rdi),%%xmm5			/* A2r = s1p13r */\n\t"\
			"movaps	     (%%rdi),%%xmm2			/* A5r = s1p01r */\n\t"\
			"movaps	0x100(%%rdi),%%xmm4			/* A3r = s1p09r */\n\t"\
			"movaps	0x080(%%rdi),%%xmm3			/* A4r = s1p05r */\n\t"\
			"\n\t"\
			"subpd	%%xmm1,%%xmm6			/* t6r = A1r-A6r */\n\t"\
			"addpd	%%xmm1,%%xmm1			/*         2*A6r */\n\t"\
			"addpd	%%xmm6,%%xmm1			/* t1r = A1r+A6r */\n\t"\
			"\n\t"\
			"subpd	%%xmm2,%%xmm5			/* t5r = A2r-A5r */\n\t"\
			"addpd	%%xmm2,%%xmm2			/*         2*A5r */\n\t"\
			"addpd	%%xmm5,%%xmm2			/* t2r = A2r+A5r */\n\t"\
			"\n\t"\
			"movaps	0x280(%%rdi),%%xmm0		/* Ar0 = s1p21r */\n\t"\
			"subpd	%%xmm3,%%xmm4			/* t4r = A3r-A4r */\n\t"\
			"addpd	%%xmm3,%%xmm3			/*         2*A4r */\n\t"\
			"addpd	%%xmm4,%%xmm3			/* t3r = A3r+A4r */\n\t"\
			"\n\t"\
			"movaps	%%xmm0,    (%%rsi)		/* cpy t0 into scratch sincos slot */\n\t"\
			"movaps	%%xmm6,0x10(%%rsi)		/* cpy t6 into scratch sincos slot */\n\t"\
			"addpd	%%xmm1,%%xmm0			/*~A0 = A0+t1 */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* cpy t5 */\n\t"\
			"addpd	%%xmm2,%%xmm3			/*~t3 = t3+t2 */\n\t"\
			"subpd	%%xmm4,%%xmm5			/*~t5 = t5-t4 */\n\t"\
			"subpd	%%xmm2,%%xmm1			/*~t1 = t1-t2 */\n\t"\
			"subpd	%%xmm7,%%xmm6			/*~t6 = t6-t5 */\n\t"\
			"addpd	%%xmm2,%%xmm2			/* 2*t2 */\n\t"\
			"addpd	%%xmm7,%%xmm4			/*~t5 = t4+t5 */\n\t"\
			"addpd	%%xmm3,%%xmm0			/* B0 */\n\t"\
			"addpd	0x10(%%rsi),%%xmm5		/* t3 = (t5-t4)+t6 */\n\t"\
			"subpd	%%xmm2,%%xmm3			/*~t2 =  (t2+t3) - 2*t2 = t3-t2 */\n\t"\
			"movaps	%%xmm4,%%xmm7			/* cpy t5 */\n\t"\
			"movaps	%%xmm0,     (%%rdi)		/* <-B0 = s1p01r, xmm0 FREE */\n\t"\
			"subpd	%%xmm6,%%xmm4			/* t4 = ~t5-~t6 */\n\t"\
			"movaps	%%xmm1,%%xmm2			/* cpy ~t1 */\n\t"\
			"subpd	    (%%rsi),%%xmm0		/* r = B0 - t0 */\n\t"\
			"mulpd	0x10(%%rax),%%xmm5		/*~t3 = t3*sx0 */\n\t"\
			"addpd	%%xmm3,%%xmm2			/* ~t1+~t2 */\n\t"\
			"mulpd	    (%%rcx),%%xmm3		/* t2 = t2*cx2 */\n\t"\
			"mulpd	0x10(%%rdx),%%xmm4		/*~t4 = t4*sx3 */\n\t"\
			"mulpd	    (%%rbx),%%xmm1		/* t1 = t1*cx1 */\n\t"\
			"mulpd	0x10(%%rbx),%%xmm6		/*~t6 = t6*sx1 */\n\t"\
			"mulpd	    (%%rax),%%xmm0		/* ~r = r*(cx0-1) */\n\t"\
			"mulpd	0x10(%%rcx),%%xmm7		/*~t5 = t5*sx2 */\n\t"\
			"mulpd	    (%%rdx),%%xmm2		/* t3 =(t1+t2)*cx3 */\n\t"\
			"addpd	    (%%rdi),%%xmm0		/* t0 =~r + B0 */\n\t"\
			"addpd	%%xmm4,%%xmm6			/*~t6 = t4+t6 */\n\t"\
			"subpd	%%xmm2,%%xmm1			/* tt = t1-t3 */\n\t"\
			"subpd	%%xmm7,%%xmm4			/*~t5 = t4-t5, xmm7 FREE */\n\t"\
			"subpd	%%xmm2,%%xmm3			/* t2 = t2-t3, xmm2 FREE */\n\t"\
			"\n\t"\
			"movaps	%%xmm0,%%xmm2			/* cpy t0 */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* cpy t3 */\n\t"\
			"addpd	%%xmm1,%%xmm0			/*~t0 = t0+tt */\n\t"\
			"addpd	%%xmm6,%%xmm5			/*~t3 = t3+t6 */\n\t"\
			"addpd	%%xmm3,%%xmm1			/*~tt = tt+t2 */\n\t"\
			"addpd	%%xmm4,%%xmm6			/*      t6+t5 */\n\t"\
			"addpd	%%xmm2,%%xmm3			/*~t2 = t2+t0 */\n\t"\
			"addpd	%%xmm7,%%xmm4			/*~t5 = t5+t3 */\n\t"\
			"subpd	%%xmm1,%%xmm2			/*~t1 = t0-tt-t2 */\n\t"\
			"subpd	%%xmm6,%%xmm7			/*~t4 = t3-t6-t5 */\n\t"\
			"movaps	%%xmm0,0x080(%%rdi)			/* B1 <- t0 = s1p05r */\n\t"\
			"movaps	%%xmm5,0x300(%%rdi)			/* B6 <- t3 = s1p25r */\n\t"\
			"movaps	%%xmm2,0x100(%%rdi)			/* B2 <- t1 = s1p09r */\n\t"\
			"movaps	%%xmm7,0x280(%%rdi)			/* B5 <- t4 = s1p21r */\n\t"\
			"movaps	%%xmm3,0x180(%%rdi)			/* B3 <- t2 = s1p13r */\n\t"\
			"movaps	%%xmm4,0x200(%%rdi)			/* B4 <- t5 = s1p17r */\n\t"\
			"/************************** Imaginary Parts: ******************************************/\n\t"\
			"movaps	0x210(%%rdi),%%xmm6		/* A1i = s1p17i */\n\t"\
			"movaps	0x310(%%rdi),%%xmm1		/* A6i = s1p25i */\n\t"\
			"movaps	0x190(%%rdi),%%xmm5		/* A2i = s1p13i */\n\t"\
			"movaps	0x010(%%rdi),%%xmm2		/* A5i = s1p01i */\n\t"\
			"movaps	0x110(%%rdi),%%xmm4		/* A3i = s1p09i */\n\t"\
			"movaps	0x090(%%rdi),%%xmm3		/* A4i = s1p05i */\n\t"\
			"\n\t"\
			"subpd	%%xmm1,%%xmm6			/* t6i = A1i-A6i */\n\t"\
			"addpd	%%xmm1,%%xmm1			/*         2*A6i */\n\t"\
			"addpd	%%xmm6,%%xmm1			/* t1i = A1i+A6i */\n\t"\
			"\n\t"\
			"subpd	%%xmm2,%%xmm5			/* t5i = A2i-A5i */\n\t"\
			"addpd	%%xmm2,%%xmm2			/*         2*A5i */\n\t"\
			"addpd	%%xmm5,%%xmm2			/* t2i = A2i+A5i */\n\t"\
			"\n\t"\
			"movaps	0x290(%%rdi),%%xmm0		/* Ai0 = s1p21i */\n\t"\
			"subpd	%%xmm3,%%xmm4			/* t4i = A3i-A4i */\n\t"\
			"addpd	%%xmm3,%%xmm3			/*         2*A4i */\n\t"\
			"addpd	%%xmm4,%%xmm3			/* t3i = A3i+A4i */\n\t"\
			"\n\t"\
			"movaps	%%xmm0,    (%%rsi)		/* cpy t0 into scratch sincos slot */\n\t"\
			"movaps	%%xmm6,0x10(%%rsi)		/* cpy t6 into scratch sincos slot */\n\t"\
			"addpd	%%xmm1,%%xmm0			/*~A0 = A0+t1 */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* cpy t5 */\n\t"\
			"addpd	%%xmm2,%%xmm3			/*~t3 = t3+t2 */\n\t"\
			"subpd	%%xmm4,%%xmm5			/*~t5 = t5-t4 */\n\t"\
			"subpd	%%xmm2,%%xmm1			/*~t1 = t1-t2 */\n\t"\
			"subpd	%%xmm7,%%xmm6			/*~t6 = t6-t5 */\n\t"\
			"addpd	%%xmm2,%%xmm2			/* 2*t2 */\n\t"\
			"addpd	%%xmm7,%%xmm4			/*~t5 = t4+t5 */\n\t"\
			"addpd	%%xmm3,%%xmm0			/* B0 */\n\t"\
			"addpd	0x10(%%rsi),%%xmm5		/* t3 = (t5-t4)+t6 */\n\t"\
			"subpd	%%xmm2,%%xmm3			/*~t2 =  (t2+t3) - 2*t2 = t3-t2 */\n\t"\
			"movaps	%%xmm4,%%xmm7			/* cpy t5 */\n\t"\
			"movaps	%%xmm0,0x010(%%rdi)		/* <-B0 = s1p00i, xmm0 FREE */\n\t"\
			"subpd	%%xmm6,%%xmm4			/* t4 = ~t5-~t6 */\n\t"\
			"movaps	%%xmm1,%%xmm2			/* cpy ~t1 */\n\t"\
			"subpd	    (%%rsi),%%xmm0		/* r = B0 - t0 */\n\t"\
			"mulpd	0x10(%%rax),%%xmm5		/*~t3 = t3*sx0 */\n\t"\
			"addpd	%%xmm3,%%xmm2			/* ~t1+~t2 */\n\t"\
			"mulpd	    (%%rcx),%%xmm3		/* t2 = t2*cx2 */\n\t"\
			"mulpd	0x10(%%rdx),%%xmm4		/*~t4 = t4*sx3 */\n\t"\
			"mulpd	    (%%rbx),%%xmm1		/* t1 = t1*cx1 */\n\t"\
			"mulpd	0x10(%%rbx),%%xmm6		/*~t6 = t6*sx1 */\n\t"\
			"mulpd	    (%%rax),%%xmm0		/* ~r = r*(cx0-1) */\n\t"\
			"mulpd	0x10(%%rcx),%%xmm7		/*~t5 = t5*sx2 */\n\t"\
			"mulpd	    (%%rdx),%%xmm2		/* t3 =(t1+t2)*cx3 */\n\t"\
			"addpd	0x10(%%rdi),%%xmm0		/* t0 =~r + B0 */\n\t"\
			"addpd	%%xmm4,%%xmm6			/*~t6 = t4+t6 */\n\t"\
			"subpd	%%xmm2,%%xmm1			/* tt = t1-t3 */\n\t"\
			"subpd	%%xmm7,%%xmm4			/*~t5 = t4-t5, xmm7 FREE */\n\t"\
			"subpd	%%xmm2,%%xmm3			/* t2 = t2-t3, xmm2 FREE */\n\t"\
			"movaps	%%xmm0,%%xmm2			/* cpy t0 */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* cpy t3 */\n\t"\
			"addpd	%%xmm1,%%xmm0			/*~t0 = t0+tt */\n\t"\
			"addpd	%%xmm6,%%xmm5			/*~t3 = t3+t6 */\n\t"\
			"addpd	%%xmm3,%%xmm1			/*~tt = tt+t2 */\n\t"\
			"addpd	%%xmm4,%%xmm6			/*      t6+t5 */\n\t"\
			"addpd	%%xmm2,%%xmm3			/*~t2 = t2+t0 */\n\t"\
			"addpd	%%xmm7,%%xmm4			/*~t5 = t5+t3 */\n\t"\
			"subpd	%%xmm1,%%xmm2			/*~t1 = t0-tt-t2, xmm1 FREE */\n\t"\
			"subpd	%%xmm6,%%xmm7			/*~t4 = t3-t6-t5, xmm6 FREE */\n\t"\
			"\n\t"\
			"/* xmm1,6 FREE */\n\t"\
			"movaps	0x080(%%rdi),%%xmm1		/* t0r */\n\t"\
			"movaps	0x300(%%rdi),%%xmm6		/* t3r */\n\t"\
			"subpd	%%xmm5,%%xmm1			/* B1r =t0r-t3i */\n\t"\
			"subpd	%%xmm6,%%xmm0			/* B6i =t0i-t3r */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*        2*t3i */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*        2*t3r */\n\t"\
			"addpd	%%xmm1,%%xmm5			/* B6r =t0r+t3i */\n\t"\
			"addpd	%%xmm0,%%xmm6			/* B1i =t0i+t3r */\n\t"\
			"movaps	%%xmm1,0x080(%%rdi)		/* <-B1r = s1p08r */\n\t"\
			"movaps	%%xmm0,0x310(%%rdi)		/* <-B6i = s1p20r */\n\t"\
			"movaps	%%xmm5,0x300(%%rdi)		/* <-B6r = s1p20r */\n\t"\
			"movaps	%%xmm6,0x090(%%rdi)		/* <-B1i = s1p08i */\n\t"\
			"\n\t"\
			"movaps	0x100(%%rdi),%%xmm1		/* t1r */\n\t"\
			"movaps	0x280(%%rdi),%%xmm6		/* t4r */\n\t"\
			"subpd	%%xmm7,%%xmm1			/* B2r =t1r-t4i */\n\t"\
			"subpd	%%xmm6,%%xmm2			/* B5i =t1i-t4r */\n\t"\
			"addpd	%%xmm7,%%xmm7			/*        2*t4i */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*        2*t4r */\n\t"\
			"addpd	%%xmm1,%%xmm7			/* B5r =t1r+t4i */\n\t"\
			"addpd	%%xmm2,%%xmm6			/* B2i =t1i+t4r */\n\t"\
			"movaps	%%xmm1,0x100(%%rdi)		/* <-B2r = s1p16r */\n\t"\
			"movaps	%%xmm2,0x290(%%rdi)		/* <-B5i = s1p12r */\n\t"\
			"movaps	%%xmm7,0x280(%%rdi)		/* <-B5r = s1p12r */\n\t"\
			"movaps	%%xmm6,0x110(%%rdi)		/* <-B2i = s1p16i */\n\t"\
			"\n\t"\
			"/* Note the order reversal on this pair of outputs: */\n\t"\
			"movaps	0x180(%%rdi),%%xmm0		/* t2r */\n\t"\
			"movaps	0x200(%%rdi),%%xmm5		/* t5r */\n\t"\
			"subpd	%%xmm4,%%xmm0			/* B4r =t2r-t5i */\n\t"\
			"subpd	%%xmm5,%%xmm3			/* B3i =t2i-t5r */\n\t"\
			"addpd	%%xmm4,%%xmm4			/*        2*t5i */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*        2*t5r */\n\t"\
			"addpd	%%xmm0,%%xmm4			/* B3r =t2r+t5i */\n\t"\
			"addpd	%%xmm3,%%xmm5			/* B4i =t2i+t5r */\n\t"\
			"movaps	%%xmm0,0x200(%%rdi)		/* <-B4r = s1p04r */\n\t"\
			"movaps	%%xmm3,0x190(%%rdi)		/* <-B3i = s1p24r */\n\t"\
			"movaps	%%xmm4,0x180(%%rdi)		/* <-B3r = s1p24r */\n\t"\
			"movaps	%%xmm5,0x210(%%rdi)		/* <-B4i = s1p04i */\n\t"\
			"\n\t"\
		"/* SSE2_RADIX_07_DFT(s1p14r,s1p10r,s1p06r,s1p02r,s1p26r,s1p22r,s1p18r,cc0,s1p02r,s1p06r,s1p10r,s1p14r,s1p18r,s1p22r,s1p26r): */\n\t"\
			"addq	$0x20,%%rdi	/* s1p02r */\n\t"\
			"movaps	0x100(%%rdi),%%xmm6			/* A1r = s1p10r */\n\t"\
			"movaps	0x200(%%rdi),%%xmm1			/* A6r = s1p18r */\n\t"\
			"movaps	0x080(%%rdi),%%xmm5			/* A2r = s1p06r */\n\t"\
			"movaps	0x280(%%rdi),%%xmm2			/* A5r = s1p22r */\n\t"\
			"movaps	     (%%rdi),%%xmm4			/* A3r = s1p02r */\n\t"\
			"movaps	0x300(%%rdi),%%xmm3			/* A4r = s1p26r */\n\t"\
			"\n\t"\
			"subpd	%%xmm1,%%xmm6			/* t6r = A1r-A6r */\n\t"\
			"addpd	%%xmm1,%%xmm1			/*         2*A6r */\n\t"\
			"addpd	%%xmm6,%%xmm1			/* t1r = A1r+A6r */\n\t"\
			"\n\t"\
			"subpd	%%xmm2,%%xmm5			/* t5r = A2r-A5r */\n\t"\
			"addpd	%%xmm2,%%xmm2			/*         2*A5r */\n\t"\
			"addpd	%%xmm5,%%xmm2			/* t2r = A2r+A5r */\n\t"\
			"\n\t"\
			"movaps	0x180(%%rdi),%%xmm0		/* Ar0 = s1p14r */\n\t"\
			"subpd	%%xmm3,%%xmm4			/* t4r = A3r-A4r */\n\t"\
			"addpd	%%xmm3,%%xmm3			/*         2*A4r */\n\t"\
			"addpd	%%xmm4,%%xmm3			/* t3r = A3r+A4r */\n\t"\
			"\n\t"\
			"movaps	%%xmm0,    (%%rsi)		/* cpy t0 into scratch sincos slot */\n\t"\
			"movaps	%%xmm6,0x10(%%rsi)		/* cpy t6 into scratch sincos slot */\n\t"\
			"addpd	%%xmm1,%%xmm0			/*~A0 = A0+t1 */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* cpy t5 */\n\t"\
			"addpd	%%xmm2,%%xmm3			/*~t3 = t3+t2 */\n\t"\
			"subpd	%%xmm4,%%xmm5			/*~t5 = t5-t4 */\n\t"\
			"subpd	%%xmm2,%%xmm1			/*~t1 = t1-t2 */\n\t"\
			"subpd	%%xmm7,%%xmm6			/*~t6 = t6-t5 */\n\t"\
			"addpd	%%xmm2,%%xmm2			/* 2*t2 */\n\t"\
			"addpd	%%xmm7,%%xmm4			/*~t5 = t4+t5 */\n\t"\
			"addpd	%%xmm3,%%xmm0			/* B0 */\n\t"\
			"addpd	0x10(%%rsi),%%xmm5		/* t3 = (t5-t4)+t6 */\n\t"\
			"subpd	%%xmm2,%%xmm3			/*~t2 =  (t2+t3) - 2*t2 = t3-t2 */\n\t"\
			"movaps	%%xmm4,%%xmm7			/* cpy t5 */\n\t"\
			"movaps	%%xmm0,     (%%rdi)		/* <-B0 = s1p01r, xmm0 FREE */\n\t"\
			"subpd	%%xmm6,%%xmm4			/* t4 = ~t5-~t6 */\n\t"\
			"movaps	%%xmm1,%%xmm2			/* cpy ~t1 */\n\t"\
			"subpd	    (%%rsi),%%xmm0		/* r = B0 - t0 */\n\t"\
			"mulpd	0x10(%%rax),%%xmm5		/*~t3 = t3*sx0 */\n\t"\
			"addpd	%%xmm3,%%xmm2			/* ~t1+~t2 */\n\t"\
			"mulpd	    (%%rcx),%%xmm3		/* t2 = t2*cx2 */\n\t"\
			"mulpd	0x10(%%rdx),%%xmm4		/*~t4 = t4*sx3 */\n\t"\
			"mulpd	    (%%rbx),%%xmm1		/* t1 = t1*cx1 */\n\t"\
			"mulpd	0x10(%%rbx),%%xmm6		/*~t6 = t6*sx1 */\n\t"\
			"mulpd	    (%%rax),%%xmm0		/* ~r = r*(cx0-1) */\n\t"\
			"mulpd	0x10(%%rcx),%%xmm7		/*~t5 = t5*sx2 */\n\t"\
			"mulpd	    (%%rdx),%%xmm2		/* t3 =(t1+t2)*cx3 */\n\t"\
			"addpd	    (%%rdi),%%xmm0		/* t0 =~r + B0 */\n\t"\
			"addpd	%%xmm4,%%xmm6			/*~t6 = t4+t6 */\n\t"\
			"subpd	%%xmm2,%%xmm1			/* tt = t1-t3 */\n\t"\
			"subpd	%%xmm7,%%xmm4			/*~t5 = t4-t5, xmm7 FREE */\n\t"\
			"subpd	%%xmm2,%%xmm3			/* t2 = t2-t3, xmm2 FREE */\n\t"\
			"\n\t"\
			"movaps	%%xmm0,%%xmm2			/* cpy t0 */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* cpy t3 */\n\t"\
			"addpd	%%xmm1,%%xmm0			/*~t0 = t0+tt */\n\t"\
			"addpd	%%xmm6,%%xmm5			/*~t3 = t3+t6 */\n\t"\
			"addpd	%%xmm3,%%xmm1			/*~tt = tt+t2 */\n\t"\
			"addpd	%%xmm4,%%xmm6			/*      t6+t5 */\n\t"\
			"addpd	%%xmm2,%%xmm3			/*~t2 = t2+t0 */\n\t"\
			"addpd	%%xmm7,%%xmm4			/*~t5 = t5+t3 */\n\t"\
			"subpd	%%xmm1,%%xmm2			/*~t1 = t0-tt-t2 */\n\t"\
			"subpd	%%xmm6,%%xmm7			/*~t4 = t3-t6-t5 */\n\t"\
			"movaps	%%xmm0,0x080(%%rdi)			/* B1 <- t0 = s1p05r */\n\t"\
			"movaps	%%xmm5,0x300(%%rdi)			/* B6 <- t3 = s1p25r */\n\t"\
			"movaps	%%xmm2,0x100(%%rdi)			/* B2 <- t1 = s1p09r */\n\t"\
			"movaps	%%xmm7,0x280(%%rdi)			/* B5 <- t4 = s1p21r */\n\t"\
			"movaps	%%xmm3,0x180(%%rdi)			/* B3 <- t2 = s1p13r */\n\t"\
			"movaps	%%xmm4,0x200(%%rdi)			/* B4 <- t5 = s1p17r */\n\t"\
			"/************************** Imaginary Parts: ******************************************/\n\t"\
			"movaps	0x110(%%rdi),%%xmm6		/* A1i = s1p10i */\n\t"\
			"movaps	0x210(%%rdi),%%xmm1		/* A6i = s1p18i */\n\t"\
			"movaps	0x090(%%rdi),%%xmm5		/* A2i = s1p06i */\n\t"\
			"movaps	0x290(%%rdi),%%xmm2		/* A5i = s1p22i */\n\t"\
			"movaps	0x010(%%rdi),%%xmm4		/* A3i = s1p02i */\n\t"\
			"movaps	0x310(%%rdi),%%xmm3		/* A4i = s1p26i */\n\t"\
			"\n\t"\
			"subpd	%%xmm1,%%xmm6			/* t6i = A1i-A6i */\n\t"\
			"addpd	%%xmm1,%%xmm1			/*         2*A6i */\n\t"\
			"addpd	%%xmm6,%%xmm1			/* t1i = A1i+A6i */\n\t"\
			"\n\t"\
			"subpd	%%xmm2,%%xmm5			/* t5i = A2i-A5i */\n\t"\
			"addpd	%%xmm2,%%xmm2			/*         2*A5i */\n\t"\
			"addpd	%%xmm5,%%xmm2			/* t2i = A2i+A5i */\n\t"\
			"\n\t"\
			"movaps	0x190(%%rdi),%%xmm0		/* Ai0 = s1p14i */\n\t"\
			"subpd	%%xmm3,%%xmm4			/* t4i = A3i-A4i */\n\t"\
			"addpd	%%xmm3,%%xmm3			/*         2*A4i */\n\t"\
			"addpd	%%xmm4,%%xmm3			/* t3i = A3i+A4i */\n\t"\
			"\n\t"\
			"movaps	%%xmm0,    (%%rsi)		/* cpy t0 into scratch sincos slot */\n\t"\
			"movaps	%%xmm6,0x10(%%rsi)		/* cpy t6 into scratch sincos slot */\n\t"\
			"addpd	%%xmm1,%%xmm0			/*~A0 = A0+t1 */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* cpy t5 */\n\t"\
			"addpd	%%xmm2,%%xmm3			/*~t3 = t3+t2 */\n\t"\
			"subpd	%%xmm4,%%xmm5			/*~t5 = t5-t4 */\n\t"\
			"subpd	%%xmm2,%%xmm1			/*~t1 = t1-t2 */\n\t"\
			"subpd	%%xmm7,%%xmm6			/*~t6 = t6-t5 */\n\t"\
			"addpd	%%xmm2,%%xmm2			/* 2*t2 */\n\t"\
			"addpd	%%xmm7,%%xmm4			/*~t5 = t4+t5 */\n\t"\
			"addpd	%%xmm3,%%xmm0			/* B0 */\n\t"\
			"addpd	0x10(%%rsi),%%xmm5		/* t3 = (t5-t4)+t6 */\n\t"\
			"subpd	%%xmm2,%%xmm3			/*~t2 =  (t2+t3) - 2*t2 = t3-t2 */\n\t"\
			"movaps	%%xmm4,%%xmm7			/* cpy t5 */\n\t"\
			"movaps	%%xmm0,0x010(%%rdi)		/* <-B0 = s1p00i, xmm0 FREE */\n\t"\
			"subpd	%%xmm6,%%xmm4			/* t4 = ~t5-~t6 */\n\t"\
			"movaps	%%xmm1,%%xmm2			/* cpy ~t1 */\n\t"\
			"subpd	    (%%rsi),%%xmm0		/* r = B0 - t0 */\n\t"\
			"mulpd	0x10(%%rax),%%xmm5		/*~t3 = t3*sx0 */\n\t"\
			"addpd	%%xmm3,%%xmm2			/* ~t1+~t2 */\n\t"\
			"mulpd	    (%%rcx),%%xmm3		/* t2 = t2*cx2 */\n\t"\
			"mulpd	0x10(%%rdx),%%xmm4		/*~t4 = t4*sx3 */\n\t"\
			"mulpd	    (%%rbx),%%xmm1		/* t1 = t1*cx1 */\n\t"\
			"mulpd	0x10(%%rbx),%%xmm6		/*~t6 = t6*sx1 */\n\t"\
			"mulpd	    (%%rax),%%xmm0		/* ~r = r*(cx0-1) */\n\t"\
			"mulpd	0x10(%%rcx),%%xmm7		/*~t5 = t5*sx2 */\n\t"\
			"mulpd	    (%%rdx),%%xmm2		/* t3 =(t1+t2)*cx3 */\n\t"\
			"addpd	0x10(%%rdi),%%xmm0		/* t0 =~r + B0 */\n\t"\
			"addpd	%%xmm4,%%xmm6			/*~t6 = t4+t6 */\n\t"\
			"subpd	%%xmm2,%%xmm1			/* tt = t1-t3 */\n\t"\
			"subpd	%%xmm7,%%xmm4			/*~t5 = t4-t5, xmm7 FREE */\n\t"\
			"subpd	%%xmm2,%%xmm3			/* t2 = t2-t3, xmm2 FREE */\n\t"\
			"movaps	%%xmm0,%%xmm2			/* cpy t0 */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* cpy t3 */\n\t"\
			"addpd	%%xmm1,%%xmm0			/*~t0 = t0+tt */\n\t"\
			"addpd	%%xmm6,%%xmm5			/*~t3 = t3+t6 */\n\t"\
			"addpd	%%xmm3,%%xmm1			/*~tt = tt+t2 */\n\t"\
			"addpd	%%xmm4,%%xmm6			/*      t6+t5 */\n\t"\
			"addpd	%%xmm2,%%xmm3			/*~t2 = t2+t0 */\n\t"\
			"addpd	%%xmm7,%%xmm4			/*~t5 = t5+t3 */\n\t"\
			"subpd	%%xmm1,%%xmm2			/*~t1 = t0-tt-t2, xmm1 FREE */\n\t"\
			"subpd	%%xmm6,%%xmm7			/*~t4 = t3-t6-t5, xmm6 FREE */\n\t"\
			"\n\t"\
			"/* xmm1,6 FREE */\n\t"\
			"movaps	0x080(%%rdi),%%xmm1		/* t0r */\n\t"\
			"movaps	0x300(%%rdi),%%xmm6		/* t3r */\n\t"\
			"subpd	%%xmm5,%%xmm1			/* B1r =t0r-t3i */\n\t"\
			"subpd	%%xmm6,%%xmm0			/* B6i =t0i-t3r */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*        2*t3i */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*        2*t3r */\n\t"\
			"addpd	%%xmm1,%%xmm5			/* B6r =t0r+t3i */\n\t"\
			"addpd	%%xmm0,%%xmm6			/* B1i =t0i+t3r */\n\t"\
			"movaps	%%xmm1,0x080(%%rdi)		/* <-B1r = s1p08r */\n\t"\
			"movaps	%%xmm0,0x310(%%rdi)		/* <-B6i = s1p20r */\n\t"\
			"movaps	%%xmm5,0x300(%%rdi)		/* <-B6r = s1p20r */\n\t"\
			"movaps	%%xmm6,0x090(%%rdi)		/* <-B1i = s1p08i */\n\t"\
			"\n\t"\
			"movaps	0x100(%%rdi),%%xmm1		/* t1r */\n\t"\
			"movaps	0x280(%%rdi),%%xmm6		/* t4r */\n\t"\
			"subpd	%%xmm7,%%xmm1			/* B2r =t1r-t4i */\n\t"\
			"subpd	%%xmm6,%%xmm2			/* B5i =t1i-t4r */\n\t"\
			"addpd	%%xmm7,%%xmm7			/*        2*t4i */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*        2*t4r */\n\t"\
			"addpd	%%xmm1,%%xmm7			/* B5r =t1r+t4i */\n\t"\
			"addpd	%%xmm2,%%xmm6			/* B2i =t1i+t4r */\n\t"\
			"movaps	%%xmm1,0x100(%%rdi)		/* <-B2r = s1p16r */\n\t"\
			"movaps	%%xmm2,0x290(%%rdi)		/* <-B5i = s1p12r */\n\t"\
			"movaps	%%xmm7,0x280(%%rdi)		/* <-B5r = s1p12r */\n\t"\
			"movaps	%%xmm6,0x110(%%rdi)		/* <-B2i = s1p16i */\n\t"\
			"\n\t"\
			"/* Note the order reversal on this pair of outputs: */\n\t"\
			"movaps	0x180(%%rdi),%%xmm0		/* t2r */\n\t"\
			"movaps	0x200(%%rdi),%%xmm5		/* t5r */\n\t"\
			"subpd	%%xmm4,%%xmm0			/* B4r =t2r-t5i */\n\t"\
			"subpd	%%xmm5,%%xmm3			/* B3i =t2i-t5r */\n\t"\
			"addpd	%%xmm4,%%xmm4			/*        2*t5i */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*        2*t5r */\n\t"\
			"addpd	%%xmm0,%%xmm4			/* B3r =t2r+t5i */\n\t"\
			"addpd	%%xmm3,%%xmm5			/* B4i =t2i+t5r */\n\t"\
			"movaps	%%xmm0,0x200(%%rdi)		/* <-B4r = s1p04r */\n\t"\
			"movaps	%%xmm3,0x190(%%rdi)		/* <-B3i = s1p24r */\n\t"\
			"movaps	%%xmm4,0x180(%%rdi)		/* <-B3r = s1p24r */\n\t"\
			"movaps	%%xmm5,0x210(%%rdi)		/* <-B4i = s1p04i */\n\t"\
			"\n\t"\
		"/* SSE2_RADIX_07_DFT(s1p07r,s1p03r,s1p27r,s1p23r,s1p19r,s1p15r,s1p11r,cc0,s1p03r,s1p07r,s1p11r,s1p15r,s1p19r,s1p23r,s1p27r): */\n\t"\
			"addq	$0x20,%%rdi	/* s1p03r */\n\t"\
			"movaps	     (%%rdi),%%xmm6			/* A1r = s1p03r */\n\t"\
			"movaps	0x100(%%rdi),%%xmm1			/* A6r = s1p11r */\n\t"\
			"movaps	0x300(%%rdi),%%xmm5			/* A2r = s1p27r */\n\t"\
			"movaps	0x180(%%rdi),%%xmm2			/* A5r = s1p15r */\n\t"\
			"movaps	0x280(%%rdi),%%xmm4			/* A3r = s1p23r */\n\t"\
			"movaps	0x200(%%rdi),%%xmm3			/* A4r = s1p19r */\n\t"\
			"\n\t"\
			"subpd	%%xmm1,%%xmm6			/* t6r = A1r-A6r */\n\t"\
			"addpd	%%xmm1,%%xmm1			/*         2*A6r */\n\t"\
			"addpd	%%xmm6,%%xmm1			/* t1r = A1r+A6r */\n\t"\
			"\n\t"\
			"subpd	%%xmm2,%%xmm5			/* t5r = A2r-A5r */\n\t"\
			"addpd	%%xmm2,%%xmm2			/*         2*A5r */\n\t"\
			"addpd	%%xmm5,%%xmm2			/* t2r = A2r+A5r */\n\t"\
			"\n\t"\
			"movaps	0x080(%%rdi),%%xmm0		/* Ar0 = s1p07r */\n\t"\
			"subpd	%%xmm3,%%xmm4			/* t4r = A3r-A4r */\n\t"\
			"addpd	%%xmm3,%%xmm3			/*         2*A4r */\n\t"\
			"addpd	%%xmm4,%%xmm3			/* t3r = A3r+A4r */\n\t"\
			"\n\t"\
			"movaps	%%xmm0,    (%%rsi)		/* cpy t0 into scratch sincos slot */\n\t"\
			"movaps	%%xmm6,0x10(%%rsi)		/* cpy t6 into scratch sincos slot */\n\t"\
			"addpd	%%xmm1,%%xmm0			/*~A0 = A0+t1 */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* cpy t5 */\n\t"\
			"addpd	%%xmm2,%%xmm3			/*~t3 = t3+t2 */\n\t"\
			"subpd	%%xmm4,%%xmm5			/*~t5 = t5-t4 */\n\t"\
			"subpd	%%xmm2,%%xmm1			/*~t1 = t1-t2 */\n\t"\
			"subpd	%%xmm7,%%xmm6			/*~t6 = t6-t5 */\n\t"\
			"addpd	%%xmm2,%%xmm2			/* 2*t2 */\n\t"\
			"addpd	%%xmm7,%%xmm4			/*~t5 = t4+t5 */\n\t"\
			"addpd	%%xmm3,%%xmm0			/* B0 */\n\t"\
			"addpd	0x10(%%rsi),%%xmm5		/* t3 = (t5-t4)+t6 */\n\t"\
			"subpd	%%xmm2,%%xmm3			/*~t2 =  (t2+t3) - 2*t2 = t3-t2 */\n\t"\
			"movaps	%%xmm4,%%xmm7			/* cpy t5 */\n\t"\
			"movaps	%%xmm0,     (%%rdi)		/* <-B0 = s1p01r, xmm0 FREE */\n\t"\
			"subpd	%%xmm6,%%xmm4			/* t4 = ~t5-~t6 */\n\t"\
			"movaps	%%xmm1,%%xmm2			/* cpy ~t1 */\n\t"\
			"subpd	    (%%rsi),%%xmm0		/* r = B0 - t0 */\n\t"\
			"mulpd	0x10(%%rax),%%xmm5		/*~t3 = t3*sx0 */\n\t"\
			"addpd	%%xmm3,%%xmm2			/* ~t1+~t2 */\n\t"\
			"mulpd	    (%%rcx),%%xmm3		/* t2 = t2*cx2 */\n\t"\
			"mulpd	0x10(%%rdx),%%xmm4		/*~t4 = t4*sx3 */\n\t"\
			"mulpd	    (%%rbx),%%xmm1		/* t1 = t1*cx1 */\n\t"\
			"mulpd	0x10(%%rbx),%%xmm6		/*~t6 = t6*sx1 */\n\t"\
			"mulpd	    (%%rax),%%xmm0		/* ~r = r*(cx0-1) */\n\t"\
			"mulpd	0x10(%%rcx),%%xmm7		/*~t5 = t5*sx2 */\n\t"\
			"mulpd	    (%%rdx),%%xmm2		/* t3 =(t1+t2)*cx3 */\n\t"\
			"addpd	    (%%rdi),%%xmm0		/* t0 =~r + B0 */\n\t"\
			"addpd	%%xmm4,%%xmm6			/*~t6 = t4+t6 */\n\t"\
			"subpd	%%xmm2,%%xmm1			/* tt = t1-t3 */\n\t"\
			"subpd	%%xmm7,%%xmm4			/*~t5 = t4-t5, xmm7 FREE */\n\t"\
			"subpd	%%xmm2,%%xmm3			/* t2 = t2-t3, xmm2 FREE */\n\t"\
			"\n\t"\
			"movaps	%%xmm0,%%xmm2			/* cpy t0 */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* cpy t3 */\n\t"\
			"addpd	%%xmm1,%%xmm0			/*~t0 = t0+tt */\n\t"\
			"addpd	%%xmm6,%%xmm5			/*~t3 = t3+t6 */\n\t"\
			"addpd	%%xmm3,%%xmm1			/*~tt = tt+t2 */\n\t"\
			"addpd	%%xmm4,%%xmm6			/*      t6+t5 */\n\t"\
			"addpd	%%xmm2,%%xmm3			/*~t2 = t2+t0 */\n\t"\
			"addpd	%%xmm7,%%xmm4			/*~t5 = t5+t3 */\n\t"\
			"subpd	%%xmm1,%%xmm2			/*~t1 = t0-tt-t2 */\n\t"\
			"subpd	%%xmm6,%%xmm7			/*~t4 = t3-t6-t5 */\n\t"\
			"movaps	%%xmm0,0x080(%%rdi)			/* B1 <- t0 = s1p05r */\n\t"\
			"movaps	%%xmm5,0x300(%%rdi)			/* B6 <- t3 = s1p25r */\n\t"\
			"movaps	%%xmm2,0x100(%%rdi)			/* B2 <- t1 = s1p09r */\n\t"\
			"movaps	%%xmm7,0x280(%%rdi)			/* B5 <- t4 = s1p21r */\n\t"\
			"movaps	%%xmm3,0x180(%%rdi)			/* B3 <- t2 = s1p13r */\n\t"\
			"movaps	%%xmm4,0x200(%%rdi)			/* B4 <- t5 = s1p17r */\n\t"\
			"/************************** Imaginary Parts: ******************************************/\n\t"\
			"movaps	0x010(%%rdi),%%xmm6			/* A1i = s1p03i */\n\t"\
			"movaps	0x110(%%rdi),%%xmm1			/* A6i = s1p11i */\n\t"\
			"movaps	0x310(%%rdi),%%xmm5			/* A2i = s1p27i */\n\t"\
			"movaps	0x190(%%rdi),%%xmm2			/* A5i = s1p15i */\n\t"\
			"movaps	0x290(%%rdi),%%xmm4			/* A3i = s1p23i */\n\t"\
			"movaps	0x210(%%rdi),%%xmm3			/* A4i = s1p19i */\n\t"\
			"\n\t"\
			"subpd	%%xmm1,%%xmm6			/* t6i = A1i-A6i */\n\t"\
			"addpd	%%xmm1,%%xmm1			/*         2*A6i */\n\t"\
			"addpd	%%xmm6,%%xmm1			/* t1i = A1i+A6i */\n\t"\
			"\n\t"\
			"subpd	%%xmm2,%%xmm5			/* t5i = A2i-A5i */\n\t"\
			"addpd	%%xmm2,%%xmm2			/*         2*A5i */\n\t"\
			"addpd	%%xmm5,%%xmm2			/* t2i = A2i+A5i */\n\t"\
			"\n\t"\
			"movaps	0x090(%%rdi),%%xmm0		/* Ai0 = s1p07i */\n\t"\
			"subpd	%%xmm3,%%xmm4			/* t4i = A3i-A4i */\n\t"\
			"addpd	%%xmm3,%%xmm3			/*         2*A4i */\n\t"\
			"addpd	%%xmm4,%%xmm3			/* t3i = A3i+A4i */\n\t"\
			"\n\t"\
			"movaps	%%xmm0,    (%%rsi)		/* cpy t0 into scratch sincos slot */\n\t"\
			"movaps	%%xmm6,0x10(%%rsi)		/* cpy t6 into scratch sincos slot */\n\t"\
			"addpd	%%xmm1,%%xmm0			/*~A0 = A0+t1 */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* cpy t5 */\n\t"\
			"addpd	%%xmm2,%%xmm3			/*~t3 = t3+t2 */\n\t"\
			"subpd	%%xmm4,%%xmm5			/*~t5 = t5-t4 */\n\t"\
			"subpd	%%xmm2,%%xmm1			/*~t1 = t1-t2 */\n\t"\
			"subpd	%%xmm7,%%xmm6			/*~t6 = t6-t5 */\n\t"\
			"addpd	%%xmm2,%%xmm2			/* 2*t2 */\n\t"\
			"addpd	%%xmm7,%%xmm4			/*~t5 = t4+t5 */\n\t"\
			"addpd	%%xmm3,%%xmm0			/* B0 */\n\t"\
			"addpd	0x10(%%rsi),%%xmm5		/* t3 = (t5-t4)+t6 */\n\t"\
			"subpd	%%xmm2,%%xmm3			/*~t2 =  (t2+t3) - 2*t2 = t3-t2 */\n\t"\
			"movaps	%%xmm4,%%xmm7			/* cpy t5 */\n\t"\
			"movaps	%%xmm0,0x010(%%rdi)		/* <-B0 = s1p00i, xmm0 FREE */\n\t"\
			"subpd	%%xmm6,%%xmm4			/* t4 = ~t5-~t6 */\n\t"\
			"movaps	%%xmm1,%%xmm2			/* cpy ~t1 */\n\t"\
			"subpd	    (%%rsi),%%xmm0		/* r = B0 - t0 */\n\t"\
			"mulpd	0x10(%%rax),%%xmm5		/*~t3 = t3*sx0 */\n\t"\
			"addpd	%%xmm3,%%xmm2			/* ~t1+~t2 */\n\t"\
			"mulpd	    (%%rcx),%%xmm3		/* t2 = t2*cx2 */\n\t"\
			"mulpd	0x10(%%rdx),%%xmm4		/*~t4 = t4*sx3 */\n\t"\
			"mulpd	    (%%rbx),%%xmm1		/* t1 = t1*cx1 */\n\t"\
			"mulpd	0x10(%%rbx),%%xmm6		/*~t6 = t6*sx1 */\n\t"\
			"mulpd	    (%%rax),%%xmm0		/* ~r = r*(cx0-1) */\n\t"\
			"mulpd	0x10(%%rcx),%%xmm7		/*~t5 = t5*sx2 */\n\t"\
			"mulpd	    (%%rdx),%%xmm2		/* t3 =(t1+t2)*cx3 */\n\t"\
			"addpd	0x10(%%rdi),%%xmm0		/* t0 =~r + B0 */\n\t"\
			"addpd	%%xmm4,%%xmm6			/*~t6 = t4+t6 */\n\t"\
			"subpd	%%xmm2,%%xmm1			/* tt = t1-t3 */\n\t"\
			"subpd	%%xmm7,%%xmm4			/*~t5 = t4-t5, xmm7 FREE */\n\t"\
			"subpd	%%xmm2,%%xmm3			/* t2 = t2-t3, xmm2 FREE */\n\t"\
			"movaps	%%xmm0,%%xmm2			/* cpy t0 */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* cpy t3 */\n\t"\
			"addpd	%%xmm1,%%xmm0			/*~t0 = t0+tt */\n\t"\
			"addpd	%%xmm6,%%xmm5			/*~t3 = t3+t6 */\n\t"\
			"addpd	%%xmm3,%%xmm1			/*~tt = tt+t2 */\n\t"\
			"addpd	%%xmm4,%%xmm6			/*      t6+t5 */\n\t"\
			"addpd	%%xmm2,%%xmm3			/*~t2 = t2+t0 */\n\t"\
			"addpd	%%xmm7,%%xmm4			/*~t5 = t5+t3 */\n\t"\
			"subpd	%%xmm1,%%xmm2			/*~t1 = t0-tt-t2, xmm1 FREE */\n\t"\
			"subpd	%%xmm6,%%xmm7			/*~t4 = t3-t6-t5, xmm6 FREE */\n\t"\
			"\n\t"\
			"/* xmm1,6 FREE */\n\t"\
			"movaps	0x080(%%rdi),%%xmm1		/* t0r */\n\t"\
			"movaps	0x300(%%rdi),%%xmm6		/* t3r */\n\t"\
			"subpd	%%xmm5,%%xmm1			/* B1r =t0r-t3i */\n\t"\
			"subpd	%%xmm6,%%xmm0			/* B6i =t0i-t3r */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*        2*t3i */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*        2*t3r */\n\t"\
			"addpd	%%xmm1,%%xmm5			/* B6r =t0r+t3i */\n\t"\
			"addpd	%%xmm0,%%xmm6			/* B1i =t0i+t3r */\n\t"\
			"movaps	%%xmm1,0x080(%%rdi)		/* <-B1r = s1p08r */\n\t"\
			"movaps	%%xmm0,0x310(%%rdi)		/* <-B6i = s1p20r */\n\t"\
			"movaps	%%xmm5,0x300(%%rdi)		/* <-B6r = s1p20r */\n\t"\
			"movaps	%%xmm6,0x090(%%rdi)		/* <-B1i = s1p08i */\n\t"\
			"\n\t"\
			"movaps	0x100(%%rdi),%%xmm1		/* t1r */\n\t"\
			"movaps	0x280(%%rdi),%%xmm6		/* t4r */\n\t"\
			"subpd	%%xmm7,%%xmm1			/* B2r =t1r-t4i */\n\t"\
			"subpd	%%xmm6,%%xmm2			/* B5i =t1i-t4r */\n\t"\
			"addpd	%%xmm7,%%xmm7			/*        2*t4i */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*        2*t4r */\n\t"\
			"addpd	%%xmm1,%%xmm7			/* B5r =t1r+t4i */\n\t"\
			"addpd	%%xmm2,%%xmm6			/* B2i =t1i+t4r */\n\t"\
			"movaps	%%xmm1,0x100(%%rdi)		/* <-B2r = s1p16r */\n\t"\
			"movaps	%%xmm2,0x290(%%rdi)		/* <-B5i = s1p12r */\n\t"\
			"movaps	%%xmm7,0x280(%%rdi)		/* <-B5r = s1p12r */\n\t"\
			"movaps	%%xmm6,0x110(%%rdi)		/* <-B2i = s1p16i */\n\t"\
			"\n\t"\
			"/* Note the order reversal on this pair of outputs: */\n\t"\
			"movaps	0x180(%%rdi),%%xmm0		/* t2r */\n\t"\
			"movaps	0x200(%%rdi),%%xmm5		/* t5r */\n\t"\
			"subpd	%%xmm4,%%xmm0			/* B4r =t2r-t5i */\n\t"\
			"subpd	%%xmm5,%%xmm3			/* B3i =t2i-t5r */\n\t"\
			"addpd	%%xmm4,%%xmm4			/*        2*t5i */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*        2*t5r */\n\t"\
			"addpd	%%xmm0,%%xmm4			/* B3r =t2r+t5i */\n\t"\
			"addpd	%%xmm3,%%xmm5			/* B4i =t2i+t5r */\n\t"\
			"movaps	%%xmm0,0x200(%%rdi)		/* <-B4r = s1p04r */\n\t"\
			"movaps	%%xmm3,0x190(%%rdi)		/* <-B3i = s1p24r */\n\t"\
			"movaps	%%xmm4,0x180(%%rdi)		/* <-B3r = s1p24r */\n\t"\
			"movaps	%%xmm5,0x210(%%rdi)		/* <-B4i = s1p04i */\n\t"\
			"\n\t"\
		"/* Since doing radix-7 in-place here, inputs of radix-4 are in consecutive memory locs, i.e. 0x20 bytes apart, e.g. the distance between s1p00r and s1p01r: */\n\t"\
			"\n\t"\
		"/* SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_E(s1p00r,s1p01r,s1p02r,s1p03r, add00+p[0,1,2,3]): */\n\t"\
			"movq	%[__out],%%rsi	/* s1p00r */\n\t"\
			"movq	%[__add],%%rax\n\t"\
			"movslq	%[__p01],%%rbx\n\t"\
			"movslq	%[__p02],%%rcx\n\t"\
			"movslq	%[__p03],%%rdx\n\t"\
			"shlq	$3,%%rbx		/* Pointer offset for floating doubles */\n\t"\
			"shlq	$3,%%rcx\n\t"\
			"shlq	$3,%%rdx\n\t"\
			"addq	%%rax,%%rbx\n\t"\
			"addq	%%rax,%%rcx\n\t"\
			"addq	%%rax,%%rdx\n\t"\
			"\n\t"\
			"movaps	    (%%rsi),%%xmm0		/* in0r = s1p00r */\n\t"\
			"movaps	0x20(%%rsi),%%xmm4		/* in1r = s1p01r */\n\t"\
			"movaps	0x10(%%rsi),%%xmm1		/* in0i = s1p00i */\n\t"\
			"movaps	0x30(%%rsi),%%xmm5		/* in1i = s1p01i */\n\t"\
			"movaps	%%xmm0,%%xmm2			/* xmm2 <- in0r */\n\t"\
			"movaps	%%xmm4,%%xmm6			/* xmm4 <- in1r */\n\t"\
			"movaps	%%xmm1,%%xmm3			/* xmm3 <- in0i */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* xmm5 <- in1i */\n\t"\
			"\n\t"\
			"addpd	0x40(%%rsi),%%xmm0		/* t1 */\n\t"\
			"addpd	0x60(%%rsi),%%xmm4		/* t5 */\n\t"\
			"addpd	0x50(%%rsi),%%xmm1		/* t2 */\n\t"\
			"addpd	0x70(%%rsi),%%xmm5		/* t6 */\n\t"\
			"subpd	0x40(%%rsi),%%xmm2		/* t3 */\n\t"\
			"subpd	0x60(%%rsi),%%xmm6		/* t7 */\n\t"\
			"subpd	0x50(%%rsi),%%xmm3		/* t4 */\n\t"\
			"subpd	0x70(%%rsi),%%xmm7		/* t8 */\n\t"\
			"/* Finish radix-4 butterfly and store results into output-array slots: */\n\t"\
			"subpd	%%xmm4,%%xmm0			/* ~t5 <- t1 -t5 */\n\t"\
			"subpd	%%xmm7,%%xmm2			/* ~t7 <- t3 -t8 */\n\t"\
			"subpd	%%xmm5,%%xmm1			/* ~t6 <- t2 -t6 */\n\t"\
			"subpd	%%xmm6,%%xmm3			/* ~t4 <- t4 -t7 */\n\t"\
			"movaps	%%xmm0,     (%%rbx)		/* <- ~t5 */\n\t"\
			"movaps	%%xmm2,     (%%rcx)		/* <- ~t7 */\n\t"\
			"movaps	%%xmm1,0x010(%%rbx)		/* <- ~t6 */\n\t"\
			"movaps	%%xmm3,0x010(%%rdx)		/* <- ~t4 */\n\t"\
			"addpd	%%xmm4,%%xmm4			/*          2*t5 */\n\t"\
			"addpd	%%xmm7,%%xmm7			/*          2*t8 */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*          2*t6 */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*          2*t7 */\n\t"\
			"addpd	%%xmm0,%%xmm4			/* ~t1 <- t1 +t5 */\n\t"\
			"addpd	%%xmm2,%%xmm7			/* ~t3 <- t3 +t8 */\n\t"\
			"addpd	%%xmm1,%%xmm5			/* ~t2 <- t2 +t6 */\n\t"\
			"addpd	%%xmm3,%%xmm6			/* ~t8 <- t4 +t7 */\n\t"\
			"movaps	%%xmm4,     (%%rax)		/* <- ~t1 */\n\t"\
			"movaps	%%xmm7,     (%%rdx)		/* <- ~t3 */\n\t"\
			"movaps	%%xmm5,0x010(%%rax)		/* <- ~t2 */\n\t"\
			"movaps	%%xmm6,0x010(%%rcx)		/* <- ~t8 */\n\t"\
		"/* SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_E(s1p04r,s1p05r,s1p06r,s1p07r, add24+p[1,0,3,2]): */\n\t"\
			"addq	$0x80,%%rsi		/* s1p04r */\n\t"\
			"movslq	%[__p24],%%rdi\n\t"\
			"shlq	$3,%%rdi		/* Pointer offset for floating doubles */\n\t"\
			"subq	%%rax,%%rbx\n\t"\
			"subq	%%rax,%%rcx\n\t"\
			"subq	%%rax,%%rdx\n\t"\
			"addq	%%rdi,%%rax		/* add24 */\n\t"\
			"addq	%%rax,%%rbx\n\t"\
			"addq	%%rax,%%rcx\n\t"\
			"addq	%%rax,%%rdx\n\t"\
			"\n\t"\
			"movaps	    (%%rsi),%%xmm0		/* in0r = s1p00r */\n\t"\
			"movaps	0x20(%%rsi),%%xmm4		/* in1r = s1p01r */\n\t"\
			"movaps	0x10(%%rsi),%%xmm1		/* in0i = s1p00i */\n\t"\
			"movaps	0x30(%%rsi),%%xmm5		/* in1i = s1p01i */\n\t"\
			"movaps	%%xmm0,%%xmm2			/* xmm2 <- in0r */\n\t"\
			"movaps	%%xmm4,%%xmm6			/* xmm4 <- in1r */\n\t"\
			"movaps	%%xmm1,%%xmm3			/* xmm3 <- in0i */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* xmm5 <- in1i */\n\t"\
			"\n\t"\
			"addpd	0x40(%%rsi),%%xmm0		/* t1 */\n\t"\
			"addpd	0x60(%%rsi),%%xmm4		/* t5 */\n\t"\
			"addpd	0x50(%%rsi),%%xmm1		/* t2 */\n\t"\
			"addpd	0x70(%%rsi),%%xmm5		/* t6 */\n\t"\
			"subpd	0x40(%%rsi),%%xmm2		/* t3 */\n\t"\
			"subpd	0x60(%%rsi),%%xmm6		/* t7 */\n\t"\
			"subpd	0x50(%%rsi),%%xmm3		/* t4 */\n\t"\
			"subpd	0x70(%%rsi),%%xmm7		/* t8 */\n\t"\
			"/* Finish radix-4 butterfly and store results into output-array slots: */\n\t"\
			"/* Swap _o1 <-> _o2, _o3 <-> _o4: eax <-> ebx, ecx <-> edx */\n\t"\
			"subpd	%%xmm4,%%xmm0			/* ~t5 <- t1 -t5 */\n\t"\
			"subpd	%%xmm7,%%xmm2			/* ~t7 <- t3 -t8 */\n\t"\
			"subpd	%%xmm5,%%xmm1			/* ~t6 <- t2 -t6 */\n\t"\
			"subpd	%%xmm6,%%xmm3			/* ~t4 <- t4 -t7 */\n\t"\
			"movaps	%%xmm0,     (%%rax)		/* <- ~t5 */\n\t"\
			"movaps	%%xmm2,     (%%rdx)		/* <- ~t7 */\n\t"\
			"movaps	%%xmm1,0x010(%%rax)		/* <- ~t6 */\n\t"\
			"movaps	%%xmm3,0x010(%%rcx)		/* <- ~t4 */\n\t"\
			"addpd	%%xmm4,%%xmm4			/*          2*t5 */\n\t"\
			"addpd	%%xmm7,%%xmm7			/*          2*t8 */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*          2*t6 */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*          2*t7 */\n\t"\
			"addpd	%%xmm0,%%xmm4			/* ~t1 <- t1 +t5 */\n\t"\
			"addpd	%%xmm2,%%xmm7			/* ~t3 <- t3 +t8 */\n\t"\
			"addpd	%%xmm1,%%xmm5			/* ~t2 <- t2 +t6 */\n\t"\
			"addpd	%%xmm3,%%xmm6			/* ~t8 <- t4 +t7 */\n\t"\
			"movaps	%%xmm4,     (%%rbx)		/* <- ~t1 */\n\t"\
			"movaps	%%xmm7,     (%%rcx)		/* <- ~t3 */\n\t"\
			"movaps	%%xmm5,0x010(%%rbx)		/* <- ~t2 */\n\t"\
			"movaps	%%xmm6,0x010(%%rdx)		/* <- ~t8 */\n\t"\
		"/* SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_E(s1p08r,s1p09r,s1p10r,s1p11r, add20+p[2,3,1,0]): */\n\t"\
			"addq	$0x80,%%rsi		/* s1p08r */\n\t"\
			"subq	%%rax,%%rbx\n\t"\
			"subq	%%rax,%%rcx\n\t"\
			"subq	%%rax,%%rdx\n\t"\
			"subq	%%rdi,%%rax		/* add0  */\n\t"\
			"movslq	%[__p20],%%rdi\n\t"\
			"shlq	$3,%%rdi		/* Pointer offset for floating doubles */\n\t"\
			"addq	%%rdi,%%rax		/* add20 */\n\t"\
			"addq	%%rax,%%rbx\n\t"\
			"addq	%%rax,%%rcx\n\t"\
			"addq	%%rax,%%rdx\n\t"\
			"\n\t"\
			"movaps	    (%%rsi),%%xmm0		/* in0r = s1p00r */\n\t"\
			"movaps	0x20(%%rsi),%%xmm4		/* in1r = s1p01r */\n\t"\
			"movaps	0x10(%%rsi),%%xmm1		/* in0i = s1p00i */\n\t"\
			"movaps	0x30(%%rsi),%%xmm5		/* in1i = s1p01i */\n\t"\
			"movaps	%%xmm0,%%xmm2			/* xmm2 <- in0r */\n\t"\
			"movaps	%%xmm4,%%xmm6			/* xmm4 <- in1r */\n\t"\
			"movaps	%%xmm1,%%xmm3			/* xmm3 <- in0i */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* xmm5 <- in1i */\n\t"\
			"\n\t"\
			"addpd	0x40(%%rsi),%%xmm0		/* t1 */\n\t"\
			"addpd	0x60(%%rsi),%%xmm4		/* t5 */\n\t"\
			"addpd	0x50(%%rsi),%%xmm1		/* t2 */\n\t"\
			"addpd	0x70(%%rsi),%%xmm5		/* t6 */\n\t"\
			"subpd	0x40(%%rsi),%%xmm2		/* t3 */\n\t"\
			"subpd	0x60(%%rsi),%%xmm6		/* t7 */\n\t"\
			"subpd	0x50(%%rsi),%%xmm3		/* t4 */\n\t"\
			"subpd	0x70(%%rsi),%%xmm7		/* t8 */\n\t"\
			"/* Finish radix-4 butterfly and store results into output-array slots: */\n\t"\
			"/* Swap _o0123 <-> _o2310: e[abcd]x <-> e[cdba]x */\n\t"\
			"subpd	%%xmm4,%%xmm0			/* ~t5 <- t1 -t5 */\n\t"\
			"subpd	%%xmm7,%%xmm2			/* ~t7 <- t3 -t8 */\n\t"\
			"subpd	%%xmm5,%%xmm1			/* ~t6 <- t2 -t6 */\n\t"\
			"subpd	%%xmm6,%%xmm3			/* ~t4 <- t4 -t7 */\n\t"\
			"movaps	%%xmm0,     (%%rdx)		/* <- ~t5 */\n\t"\
			"movaps	%%xmm2,     (%%rbx)		/* <- ~t7 */\n\t"\
			"movaps	%%xmm1,0x010(%%rdx)		/* <- ~t6 */\n\t"\
			"movaps	%%xmm3,0x010(%%rax)		/* <- ~t4 */\n\t"\
			"addpd	%%xmm4,%%xmm4			/*          2*t5 */\n\t"\
			"addpd	%%xmm7,%%xmm7			/*          2*t8 */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*          2*t6 */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*          2*t7 */\n\t"\
			"addpd	%%xmm0,%%xmm4			/* ~t1 <- t1 +t5 */\n\t"\
			"addpd	%%xmm2,%%xmm7			/* ~t3 <- t3 +t8 */\n\t"\
			"addpd	%%xmm1,%%xmm5			/* ~t2 <- t2 +t6 */\n\t"\
			"addpd	%%xmm3,%%xmm6			/* ~t8 <- t4 +t7 */\n\t"\
			"movaps	%%xmm4,     (%%rcx)		/* <- ~t1 */\n\t"\
			"movaps	%%xmm7,     (%%rax)		/* <- ~t3 */\n\t"\
			"movaps	%%xmm5,0x010(%%rcx)		/* <- ~t2 */\n\t"\
			"movaps	%%xmm6,0x010(%%rbx)		/* <- ~t8 */\n\t"\
		"/* SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_E(s1p12r,s1p13r,s1p14r,s1p15r, add16+p[0,1,2,3]): */\n\t"\
			"addq	$0x80,%%rsi		/* s1p12r */\n\t"\
			"subq	%%rax,%%rbx\n\t"\
			"subq	%%rax,%%rcx\n\t"\
			"subq	%%rax,%%rdx\n\t"\
			"subq	%%rdi,%%rax		/* add0  */\n\t"\
			"movslq	%[__p16],%%rdi\n\t"\
			"shlq	$3,%%rdi		/* Pointer offset for floating doubles */\n\t"\
			"addq	%%rdi,%%rax		/* add16 */\n\t"\
			"addq	%%rax,%%rbx\n\t"\
			"addq	%%rax,%%rcx\n\t"\
			"addq	%%rax,%%rdx\n\t"\
			"\n\t"\
			"movaps	    (%%rsi),%%xmm0		/* in0r = s1p00r */\n\t"\
			"movaps	0x20(%%rsi),%%xmm4		/* in1r = s1p01r */\n\t"\
			"movaps	0x10(%%rsi),%%xmm1		/* in0i = s1p00i */\n\t"\
			"movaps	0x30(%%rsi),%%xmm5		/* in1i = s1p01i */\n\t"\
			"movaps	%%xmm0,%%xmm2			/* xmm2 <- in0r */\n\t"\
			"movaps	%%xmm4,%%xmm6			/* xmm4 <- in1r */\n\t"\
			"movaps	%%xmm1,%%xmm3			/* xmm3 <- in0i */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* xmm5 <- in1i */\n\t"\
			"\n\t"\
			"addpd	0x40(%%rsi),%%xmm0		/* t1 */\n\t"\
			"addpd	0x60(%%rsi),%%xmm4		/* t5 */\n\t"\
			"addpd	0x50(%%rsi),%%xmm1		/* t2 */\n\t"\
			"addpd	0x70(%%rsi),%%xmm5		/* t6 */\n\t"\
			"subpd	0x40(%%rsi),%%xmm2		/* t3 */\n\t"\
			"subpd	0x60(%%rsi),%%xmm6		/* t7 */\n\t"\
			"subpd	0x50(%%rsi),%%xmm3		/* t4 */\n\t"\
			"subpd	0x70(%%rsi),%%xmm7		/* t8 */\n\t"\
			"/* Finish radix-4 butterfly and store results into output-array slots: */\n\t"\
			"subpd	%%xmm4,%%xmm0			/* ~t5 <- t1 -t5 */\n\t"\
			"subpd	%%xmm7,%%xmm2			/* ~t7 <- t3 -t8 */\n\t"\
			"subpd	%%xmm5,%%xmm1			/* ~t6 <- t2 -t6 */\n\t"\
			"subpd	%%xmm6,%%xmm3			/* ~t4 <- t4 -t7 */\n\t"\
			"movaps	%%xmm0,     (%%rbx)		/* <- ~t5 */\n\t"\
			"movaps	%%xmm2,     (%%rcx)		/* <- ~t7 */\n\t"\
			"movaps	%%xmm1,0x010(%%rbx)		/* <- ~t6 */\n\t"\
			"movaps	%%xmm3,0x010(%%rdx)		/* <- ~t4 */\n\t"\
			"addpd	%%xmm4,%%xmm4			/*          2*t5 */\n\t"\
			"addpd	%%xmm7,%%xmm7			/*          2*t8 */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*          2*t6 */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*          2*t7 */\n\t"\
			"addpd	%%xmm0,%%xmm4			/* ~t1 <- t1 +t5 */\n\t"\
			"addpd	%%xmm2,%%xmm7			/* ~t3 <- t3 +t8 */\n\t"\
			"addpd	%%xmm1,%%xmm5			/* ~t2 <- t2 +t6 */\n\t"\
			"addpd	%%xmm3,%%xmm6			/* ~t8 <- t4 +t7 */\n\t"\
			"movaps	%%xmm4,     (%%rax)		/* <- ~t1 */\n\t"\
			"movaps	%%xmm7,     (%%rdx)		/* <- ~t3 */\n\t"\
			"movaps	%%xmm5,0x010(%%rax)		/* <- ~t2 */\n\t"\
			"movaps	%%xmm6,0x010(%%rcx)		/* <- ~t8 */\n\t"\
		"/* SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_E(s1p16r,s1p17r,s1p18r,s1p19r, add12+p[3,2,0,1]): */\n\t"\
			"addq	$0x80,%%rsi		/* s1p16r */\n\t"\
			"subq	%%rax,%%rbx\n\t"\
			"subq	%%rax,%%rcx\n\t"\
			"subq	%%rax,%%rdx\n\t"\
			"subq	%%rdi,%%rax		/* add0  */\n\t"\
			"movslq	%[__p12],%%rdi\n\t"\
			"shlq	$3,%%rdi		/* Pointer offset for floating doubles */\n\t"\
			"addq	%%rdi,%%rax		/* add12 */\n\t"\
			"addq	%%rax,%%rbx\n\t"\
			"addq	%%rax,%%rcx\n\t"\
			"addq	%%rax,%%rdx\n\t"\
			"\n\t"\
			"movaps	    (%%rsi),%%xmm0		/* in0r = s1p00r */\n\t"\
			"movaps	0x20(%%rsi),%%xmm4		/* in1r = s1p01r */\n\t"\
			"movaps	0x10(%%rsi),%%xmm1		/* in0i = s1p00i */\n\t"\
			"movaps	0x30(%%rsi),%%xmm5		/* in1i = s1p01i */\n\t"\
			"movaps	%%xmm0,%%xmm2			/* xmm2 <- in0r */\n\t"\
			"movaps	%%xmm4,%%xmm6			/* xmm4 <- in1r */\n\t"\
			"movaps	%%xmm1,%%xmm3			/* xmm3 <- in0i */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* xmm5 <- in1i */\n\t"\
			"\n\t"\
			"addpd	0x40(%%rsi),%%xmm0		/* t1 */\n\t"\
			"addpd	0x60(%%rsi),%%xmm4		/* t5 */\n\t"\
			"addpd	0x50(%%rsi),%%xmm1		/* t2 */\n\t"\
			"addpd	0x70(%%rsi),%%xmm5		/* t6 */\n\t"\
			"subpd	0x40(%%rsi),%%xmm2		/* t3 */\n\t"\
			"subpd	0x60(%%rsi),%%xmm6		/* t7 */\n\t"\
			"subpd	0x50(%%rsi),%%xmm3		/* t4 */\n\t"\
			"subpd	0x70(%%rsi),%%xmm7		/* t8 */\n\t"\
			"/* Finish radix-4 butterfly and store results into output-array slots: */\n\t"\
			"/* Swap _o0123 <-> _o3201: e[abcd]x <-> e[dcab]x */\n\t"\
			"subpd	%%xmm4,%%xmm0			/* ~t5 <- t1 -t5 */\n\t"\
			"subpd	%%xmm7,%%xmm2			/* ~t7 <- t3 -t8 */\n\t"\
			"subpd	%%xmm5,%%xmm1			/* ~t6 <- t2 -t6 */\n\t"\
			"subpd	%%xmm6,%%xmm3			/* ~t4 <- t4 -t7 */\n\t"\
			"movaps	%%xmm0,     (%%rcx)		/* <- ~t5 */\n\t"\
			"movaps	%%xmm2,     (%%rax)		/* <- ~t7 */\n\t"\
			"movaps	%%xmm1,0x010(%%rcx)		/* <- ~t6 */\n\t"\
			"movaps	%%xmm3,0x010(%%rbx)		/* <- ~t4 */\n\t"\
			"addpd	%%xmm4,%%xmm4			/*          2*t5 */\n\t"\
			"addpd	%%xmm7,%%xmm7			/*          2*t8 */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*          2*t6 */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*          2*t7 */\n\t"\
			"addpd	%%xmm0,%%xmm4			/* ~t1 <- t1 +t5 */\n\t"\
			"addpd	%%xmm2,%%xmm7			/* ~t3 <- t3 +t8 */\n\t"\
			"addpd	%%xmm1,%%xmm5			/* ~t2 <- t2 +t6 */\n\t"\
			"addpd	%%xmm3,%%xmm6			/* ~t8 <- t4 +t7 */\n\t"\
			"movaps	%%xmm4,     (%%rdx)		/* <- ~t1 */\n\t"\
			"movaps	%%xmm7,     (%%rbx)		/* <- ~t3 */\n\t"\
			"movaps	%%xmm5,0x010(%%rdx)		/* <- ~t2 */\n\t"\
			"movaps	%%xmm6,0x010(%%rax)		/* <- ~t8 */\n\t"\
		"/* SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_E(s1p20r,s1p21r,s1p22r,s1p23r, add08+p[1,0,3,2]): */\n\t"\
			"addq	$0x80,%%rsi		/* s1p20r */\n\t"\
			"subq	%%rax,%%rbx\n\t"\
			"subq	%%rax,%%rcx\n\t"\
			"subq	%%rax,%%rdx\n\t"\
			"subq	%%rdi,%%rax		/* add0  */\n\t"\
			"movslq	%[__p08],%%rdi\n\t"\
			"shlq	$3,%%rdi		/* Pointer offset for floating doubles */\n\t"\
			"addq	%%rdi,%%rax		/* add08 */\n\t"\
			"addq	%%rax,%%rbx\n\t"\
			"addq	%%rax,%%rcx\n\t"\
			"addq	%%rax,%%rdx\n\t"\
			"\n\t"\
			"movaps	    (%%rsi),%%xmm0		/* in0r = s1p00r */\n\t"\
			"movaps	0x20(%%rsi),%%xmm4		/* in1r = s1p01r */\n\t"\
			"movaps	0x10(%%rsi),%%xmm1		/* in0i = s1p00i */\n\t"\
			"movaps	0x30(%%rsi),%%xmm5		/* in1i = s1p01i */\n\t"\
			"movaps	%%xmm0,%%xmm2			/* xmm2 <- in0r */\n\t"\
			"movaps	%%xmm4,%%xmm6			/* xmm4 <- in1r */\n\t"\
			"movaps	%%xmm1,%%xmm3			/* xmm3 <- in0i */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* xmm5 <- in1i */\n\t"\
			"\n\t"\
			"addpd	0x40(%%rsi),%%xmm0		/* t1 */\n\t"\
			"addpd	0x60(%%rsi),%%xmm4		/* t5 */\n\t"\
			"addpd	0x50(%%rsi),%%xmm1		/* t2 */\n\t"\
			"addpd	0x70(%%rsi),%%xmm5		/* t6 */\n\t"\
			"subpd	0x40(%%rsi),%%xmm2		/* t3 */\n\t"\
			"subpd	0x60(%%rsi),%%xmm6		/* t7 */\n\t"\
			"subpd	0x50(%%rsi),%%xmm3		/* t4 */\n\t"\
			"subpd	0x70(%%rsi),%%xmm7		/* t8 */\n\t"\
			"/* Finish radix-4 butterfly and store results into output-array slots: */\n\t"\
			"/* Swap _o0 <-> _o1, _o2 <-> _o3: eax <-> ebx, ecx <-> edx */\n\t"\
			"subpd	%%xmm4,%%xmm0			/* ~t5 <- t1 -t5 */\n\t"\
			"subpd	%%xmm7,%%xmm2			/* ~t7 <- t3 -t8 */\n\t"\
			"subpd	%%xmm5,%%xmm1			/* ~t6 <- t2 -t6 */\n\t"\
			"subpd	%%xmm6,%%xmm3			/* ~t4 <- t4 -t7 */\n\t"\
			"movaps	%%xmm0,     (%%rax)		/* <- ~t5 */\n\t"\
			"movaps	%%xmm2,     (%%rdx)		/* <- ~t7 */\n\t"\
			"movaps	%%xmm1,0x010(%%rax)		/* <- ~t6 */\n\t"\
			"movaps	%%xmm3,0x010(%%rcx)		/* <- ~t4 */\n\t"\
			"addpd	%%xmm4,%%xmm4			/*          2*t5 */\n\t"\
			"addpd	%%xmm7,%%xmm7			/*          2*t8 */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*          2*t6 */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*          2*t7 */\n\t"\
			"addpd	%%xmm0,%%xmm4			/* ~t1 <- t1 +t5 */\n\t"\
			"addpd	%%xmm2,%%xmm7			/* ~t3 <- t3 +t8 */\n\t"\
			"addpd	%%xmm1,%%xmm5			/* ~t2 <- t2 +t6 */\n\t"\
			"addpd	%%xmm3,%%xmm6			/* ~t8 <- t4 +t7 */\n\t"\
			"movaps	%%xmm4,     (%%rbx)		/* <- ~t1 */\n\t"\
			"movaps	%%xmm7,     (%%rcx)		/* <- ~t3 */\n\t"\
			"movaps	%%xmm5,0x010(%%rbx)		/* <- ~t2 */\n\t"\
			"movaps	%%xmm6,0x010(%%rdx)		/* <- ~t8 */\n\t"\
		"/* SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_E(s1p24r,s1p25r,s1p26r,s1p27r, add04+p[2,3,1,0]): */\n\t"\
			"addq	$0x80,%%rsi		/* s1p24r */\n\t"\
			"subq	%%rax,%%rbx\n\t"\
			"subq	%%rax,%%rcx\n\t"\
			"subq	%%rax,%%rdx\n\t"\
			"subq	%%rdi,%%rax		/* add0  */\n\t"\
			"movslq	%[__p04],%%rdi\n\t"\
			"shlq	$3,%%rdi		/* Pointer offset for floating doubles */\n\t"\
			"addq	%%rdi,%%rax		/* add04 */\n\t"\
			"addq	%%rax,%%rbx\n\t"\
			"addq	%%rax,%%rcx\n\t"\
			"addq	%%rax,%%rdx\n\t"\
			"\n\t"\
			"movaps	    (%%rsi),%%xmm0		/* in0r = s1p00r */\n\t"\
			"movaps	0x20(%%rsi),%%xmm4		/* in1r = s1p01r */\n\t"\
			"movaps	0x10(%%rsi),%%xmm1		/* in0i = s1p00i */\n\t"\
			"movaps	0x30(%%rsi),%%xmm5		/* in1i = s1p01i */\n\t"\
			"movaps	%%xmm0,%%xmm2			/* xmm2 <- in0r */\n\t"\
			"movaps	%%xmm4,%%xmm6			/* xmm4 <- in1r */\n\t"\
			"movaps	%%xmm1,%%xmm3			/* xmm3 <- in0i */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* xmm5 <- in1i */\n\t"\
			"\n\t"\
			"addpd	0x40(%%rsi),%%xmm0		/* t1 */\n\t"\
			"addpd	0x60(%%rsi),%%xmm4		/* t5 */\n\t"\
			"addpd	0x50(%%rsi),%%xmm1		/* t2 */\n\t"\
			"addpd	0x70(%%rsi),%%xmm5		/* t6 */\n\t"\
			"subpd	0x40(%%rsi),%%xmm2		/* t3 */\n\t"\
			"subpd	0x60(%%rsi),%%xmm6		/* t7 */\n\t"\
			"subpd	0x50(%%rsi),%%xmm3		/* t4 */\n\t"\
			"subpd	0x70(%%rsi),%%xmm7		/* t8 */\n\t"\
			"/* Finish radix-4 butterfly and store results into output-array slots: */\n\t"\
			"/* Swap _o0123 <-> _o2310: e[abcd]x <-> e[cdba]x */\n\t"\
			"subpd	%%xmm4,%%xmm0			/* ~t5 <- t1 -t5 */\n\t"\
			"subpd	%%xmm7,%%xmm2			/* ~t7 <- t3 -t8 */\n\t"\
			"subpd	%%xmm5,%%xmm1			/* ~t6 <- t2 -t6 */\n\t"\
			"subpd	%%xmm6,%%xmm3			/* ~t4 <- t4 -t7 */\n\t"\
			"movaps	%%xmm0,     (%%rdx)		/* <- ~t5 */\n\t"\
			"movaps	%%xmm2,     (%%rbx)		/* <- ~t7 */\n\t"\
			"movaps	%%xmm1,0x010(%%rdx)		/* <- ~t6 */\n\t"\
			"movaps	%%xmm3,0x010(%%rax)		/* <- ~t4 */\n\t"\
			"addpd	%%xmm4,%%xmm4			/*          2*t5 */\n\t"\
			"addpd	%%xmm7,%%xmm7			/*          2*t8 */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*          2*t6 */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*          2*t7 */\n\t"\
			"addpd	%%xmm0,%%xmm4			/* ~t1 <- t1 +t5 */\n\t"\
			"addpd	%%xmm2,%%xmm7			/* ~t3 <- t3 +t8 */\n\t"\
			"addpd	%%xmm1,%%xmm5			/* ~t2 <- t2 +t6 */\n\t"\
			"addpd	%%xmm3,%%xmm6			/* ~t8 <- t4 +t7 */\n\t"\
			"movaps	%%xmm4,     (%%rcx)		/* <- ~t1 */\n\t"\
			"movaps	%%xmm7,     (%%rax)		/* <- ~t3 */\n\t"\
			"movaps	%%xmm5,0x010(%%rcx)		/* <- ~t2 */\n\t"\
			"movaps	%%xmm6,0x010(%%rbx)		/* <- ~t8 */\n\t"\
			"\n\t"\
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
			: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
		);\
	}

  #endif	// USE_64BIT_ASM_STYLE ?

#endif

#endif	/* radix28_ditN_cy_dif1_gcc_h_included */
