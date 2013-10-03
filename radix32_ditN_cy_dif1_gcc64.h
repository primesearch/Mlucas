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
#ifndef radix32_ditN_cy_dif1_gcc_h_included
#define radix32_ditN_cy_dif1_gcc_h_included

#ifdef USE_AVX	// AVX and AVX2 both use 256-bit registers

  #define USE_64BIT_ASM_STYLE	0

  #if !USE_64BIT_ASM_STYLE	// USE_64BIT_ASM_STYLE = False: Use simple 64-bit-ified version of the analogous 32-bit ASM macros, i.e. using just xmm0-7.

	#define SSE2_RADIX32_DIT_NOTWIDDLE(Xadd,Xp01,Xp02,Xp03,Xp04,Xp08,Xp10,Xp18,Xr00,Xisrt2,Xcc0)\
	{\
	__asm__ volatile (\
		"\n\t"/*...Block 1: */\
		"movq	%[__r00],%%rsi	\n\t"\
		"movq	%[__add],%%rax	\n\t"\
		"movslq	%[__p01],%%rbx	\n\t"\
		"movslq	%[__p02],%%rcx	\n\t"\
		"movslq	%[__p03],%%rdx	\n\t"\
		"shlq	$3,%%rbx		\n\t"\
		"shlq	$3,%%rcx		\n\t"\
		"shlq	$3,%%rdx		\n\t"\
		"addq	%%rax,%%rbx		\n\t"\
		"addq	%%rax,%%rcx		\n\t"\
		"addq	%%rax,%%rdx		\n\t"\
		"\n\t"/* SSE2_RADIX4_DIT_0TWIDDLE(add0, add1, add2, add3, r00): */\
		"vmovaps	     (%%rax),%%ymm2	\n\t"\
		"vmovaps	     (%%rcx),%%ymm6	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm3	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm7	\n\t"\
		"vmovaps	     (%%rbx),%%ymm0	\n\t"\
		"vmovaps	     (%%rdx),%%ymm4	\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm1	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm5	\n\t"\
		"vsubpd	      %%ymm0,%%ymm2,%%ymm2		\n\t"\
		"vsubpd	      %%ymm4,%%ymm6,%%ymm6		\n\t"\
		"vsubpd	      %%ymm1,%%ymm3,%%ymm3		\n\t"\
		"vsubpd	      %%ymm5,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm0,%%ymm0,%%ymm0		\n\t"\
		"vaddpd	      %%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	      %%ymm1,%%ymm1,%%ymm1		\n\t"\
		"vaddpd	      %%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm2,%%ymm0,%%ymm0		\n\t"\
		"vaddpd	      %%ymm6,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	      %%ymm3,%%ymm1,%%ymm1		\n\t"\
		"vaddpd	      %%ymm7,%%ymm5,%%ymm5		\n\t"\
		"vsubpd	      %%ymm4,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	      %%ymm7,%%ymm2,%%ymm2		\n\t"\
		"vsubpd	      %%ymm5,%%ymm1,%%ymm1		\n\t"\
		"vsubpd	      %%ymm6,%%ymm3,%%ymm3		\n\t"\
		"vmovaps	%%ymm0,0x080(%%rsi)	\n\t"\
		"vmovaps	%%ymm2,0x0c0(%%rsi)	\n\t"\
		"vmovaps	%%ymm1,0x0a0(%%rsi)	\n\t"\
		"vmovaps	%%ymm3,0x060(%%rsi)	\n\t"\
		"vaddpd	      %%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm0,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	      %%ymm2,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm1,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm3,%%ymm6,%%ymm6		\n\t"\
		"vmovaps	%%ymm4,     (%%rsi)	\n\t"\
		"vmovaps	%%ymm7,0x040(%%rsi)	\n\t"\
		"vmovaps	%%ymm5,0x020(%%rsi)	\n\t"\
		"vmovaps	%%ymm6,0x0e0(%%rsi)	\n\t"\
		"\n\t"/* SSE2_RADIX4_DIT_0TWIDDLE_2NDOFTWO(add4, add5, add6, add7, r08): */\
		"addq	$0x100,%%rsi			\n\t"\
		"movslq	%[__p04],%%rdi		\n\t"\
		"shlq	$3,%%rdi			\n\t"\
		"addq	%%rdi,%%rax			\n\t"\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	%%rdi,%%rcx			\n\t"\
		"addq	%%rdi,%%rdx			\n\t"\
		"\n\t"\
		"vmovaps	     (%%rax),%%ymm2	\n\t"\
		"vmovaps	     (%%rcx),%%ymm6	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm3	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm7	\n\t"\
		"vmovaps	     (%%rbx),%%ymm0	\n\t"\
		"vmovaps	     (%%rdx),%%ymm4	\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm1	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm5	\n\t"\
		"\n\t"\
		"movq	%[__isrt2],%%rdi	\n\t"\
		"vsubpd	      %%ymm0,%%ymm2,%%ymm2		\n\t"\
		"vsubpd	      %%ymm4,%%ymm6,%%ymm6		\n\t"\
		"vsubpd	      %%ymm1,%%ymm3,%%ymm3		\n\t"\
		"vsubpd	      %%ymm5,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm0,%%ymm0,%%ymm0		\n\t"\
		"vaddpd	      %%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	      %%ymm1,%%ymm1,%%ymm1		\n\t"\
		"vaddpd	      %%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm2,%%ymm0,%%ymm0		\n\t"\
		"vaddpd	      %%ymm6,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	      %%ymm3,%%ymm1,%%ymm1		\n\t"\
		"vaddpd	      %%ymm7,%%ymm5,%%ymm5		\n\t"\
		"vsubpd	      %%ymm4,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	      %%ymm5,%%ymm1,%%ymm1		\n\t"\
		"vmovaps	%%ymm0,0x080(%%rsi)	\n\t"\
		"vmovaps	%%ymm1,0x0a0(%%rsi)	\n\t"\
		"vsubpd	      %%ymm7,%%ymm2,%%ymm2		\n\t"\
		"vaddpd	      %%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vsubpd	      %%ymm6,%%ymm3,%%ymm3		\n\t"\
		"vaddpd	      %%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm0,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm1,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm2,%%ymm7,%%ymm7		\n\t"\
		"vmovaps	%%ymm4,     (%%rsi)	\n\t"\
		"vaddpd	      %%ymm3,%%ymm6,%%ymm6		\n\t"\
		"vmovaps	%%ymm5,0x020(%%rsi)	\n\t"\
		"vmovaps	(%%rdi),%%ymm5	\n\t"/* ISRT2 */\
		"vmovaps	%%ymm3,%%ymm0		\n\t"\
		"vmovaps	%%ymm6,%%ymm1		\n\t"\
		"vsubpd	      %%ymm7,%%ymm3,%%ymm3		\n\t"\
		"vsubpd	      %%ymm2,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm7,%%ymm0,%%ymm0		\n\t"\
		"vaddpd	      %%ymm2,%%ymm1,%%ymm1		\n\t"\
		"vmulpd	      %%ymm5,%%ymm3,%%ymm3		\n\t"\
		"vmulpd	      %%ymm5,%%ymm6,%%ymm6		\n\t"\
		"vmulpd	      %%ymm5,%%ymm0,%%ymm0		\n\t"\
		"vmulpd	      %%ymm5,%%ymm1,%%ymm1		\n\t"\
		"vmovaps	%%ymm3,0x060(%%rsi)	\n\t"\
		"vmovaps	%%ymm6,0x0e0(%%rsi)	\n\t"\
		"vmovaps	%%ymm0,0x040(%%rsi)	\n\t"\
		"vmovaps	%%ymm1,0x0c0(%%rsi)	\n\t"\
		"\n\t"/* SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS(r00,r02,r04,r06,r08,r0A,r0C,r0E): */\
		"vmovaps		-0x100(%%rsi),%%ymm0	\n\t"/* r00 */\
		"vmovaps		-0x080(%%rsi),%%ymm4	\n\t"/* r04 */\
		"vmovaps		-0x0e0(%%rsi),%%ymm1	\n\t"/* r01 */\
		"vmovaps		-0x060(%%rsi),%%ymm5	\n\t"/* r05 */\
		"vmovaps		      (%%rsi),%%ymm2	\n\t"/* r08 */\
		"vmovaps		 0x0a0(%%rsi),%%ymm7	\n\t"/* r0D */\
		"vmovaps		 0x020(%%rsi),%%ymm3	\n\t"/* r09 */\
		"vmovaps		 0x080(%%rsi),%%ymm6	\n\t"/* r0C */\
		"vsubpd	      %%ymm2,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	      %%ymm7,%%ymm4,%%ymm4		\n\t"\
		"vsubpd	      %%ymm3,%%ymm1,%%ymm1		\n\t"\
		"vsubpd	      %%ymm6,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm2,%%ymm2,%%ymm2		\n\t"\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm3,%%ymm3,%%ymm3		\n\t"\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm0,%%ymm2,%%ymm2		\n\t"\
		"vaddpd	      %%ymm4,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm1,%%ymm3,%%ymm3		\n\t"\
		"vaddpd	      %%ymm5,%%ymm6,%%ymm6		\n\t"\
		"vmovaps		%%ymm0,      (%%rsi)\n\t"\
		"vmovaps		%%ymm4, 0x080(%%rsi)\n\t"\
		"vmovaps		%%ymm1, 0x020(%%rsi)\n\t"\
		"vmovaps		%%ymm5,-0x060(%%rsi)\n\t"\
		"vmovaps		%%ymm2,-0x100(%%rsi)\n\t"\
		"vmovaps		%%ymm7,-0x080(%%rsi)\n\t"\
		"vmovaps		%%ymm3,-0x0e0(%%rsi)\n\t"\
		"vmovaps		%%ymm6, 0x0a0(%%rsi)\n\t"\
		"vmovaps		-0x0c0(%%rsi),%%ymm0	\n\t"/* r02 */\
		"vmovaps		-0x040(%%rsi),%%ymm4	\n\t"/* r06 */\
		"vmovaps		-0x0a0(%%rsi),%%ymm1	\n\t"/* r03 */\
		"vmovaps		-0x020(%%rsi),%%ymm5	\n\t"/* r07 */\
		"vmovaps		 0x040(%%rsi),%%ymm2	\n\t"/* r0A */\
		"vmovaps		 0x0e0(%%rsi),%%ymm7	\n\t"/* r0F */\
		"vmovaps		 0x060(%%rsi),%%ymm3	\n\t"/* r0B */\
		"vmovaps		 0x0c0(%%rsi),%%ymm6	\n\t"/* r0E */\
		"vsubpd	      %%ymm2,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	      %%ymm7,%%ymm4,%%ymm4		\n\t"\
		"vsubpd	      %%ymm3,%%ymm1,%%ymm1		\n\t"\
		"vsubpd	      %%ymm6,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm2,%%ymm2,%%ymm2		\n\t"\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm3,%%ymm3,%%ymm3		\n\t"\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm0,%%ymm2,%%ymm2		\n\t"\
		"vaddpd	      %%ymm4,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm1,%%ymm3,%%ymm3		\n\t"\
		"vaddpd	      %%ymm5,%%ymm6,%%ymm6		\n\t"\
		"vmovaps		%%ymm0, 0x040(%%rsi)\n\t"\
		"vmovaps		%%ymm4, 0x0c0(%%rsi)\n\t"\
		"vmovaps		%%ymm1, 0x060(%%rsi)\n\t"\
		"vmovaps		%%ymm5,-0x020(%%rsi)\n\t"\
		"vmovaps		%%ymm2,-0x0c0(%%rsi)\n\t"\
		"vmovaps		%%ymm7,-0x040(%%rsi)\n\t"\
		"vmovaps		%%ymm3,-0x0a0(%%rsi)\n\t"\
		"vmovaps		%%ymm6, 0x0e0(%%rsi)\n\t"\
		"\n\t"\
		"\n\t"/*...Block 2: */\
		"addq	$0x100,%%rsi			\n\t"\
		"movslq	%[__p04],%%rdi		\n\t"\
		"shlq	$3,%%rdi			\n\t"\
		"addq	%%rdi,%%rax			\n\t"\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	%%rdi,%%rcx			\n\t"\
		"addq	%%rdi,%%rdx			\n\t"\
		"\n\t"\
		"\n\t"/* SSE2_RADIX4_DIT_0TWIDDLE(add0, add1, add2, add3, r10): */\
		"vmovaps	     (%%rax),%%ymm2	\n\t"\
		"vmovaps	     (%%rcx),%%ymm6	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm3	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm7	\n\t"\
		"vmovaps	     (%%rbx),%%ymm0	\n\t"\
		"vmovaps	     (%%rdx),%%ymm4	\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm1	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm5	\n\t"\
		"vsubpd	      %%ymm0,%%ymm2,%%ymm2		\n\t"\
		"vsubpd	      %%ymm4,%%ymm6,%%ymm6		\n\t"\
		"vsubpd	      %%ymm1,%%ymm3,%%ymm3		\n\t"\
		"vsubpd	      %%ymm5,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm0,%%ymm0,%%ymm0		\n\t"\
		"vaddpd	      %%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	      %%ymm1,%%ymm1,%%ymm1		\n\t"\
		"vaddpd	      %%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm2,%%ymm0,%%ymm0		\n\t"\
		"vaddpd	      %%ymm6,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	      %%ymm3,%%ymm1,%%ymm1		\n\t"\
		"vaddpd	      %%ymm7,%%ymm5,%%ymm5		\n\t"\
		"vsubpd	      %%ymm4,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	      %%ymm7,%%ymm2,%%ymm2		\n\t"\
		"vsubpd	      %%ymm5,%%ymm1,%%ymm1		\n\t"\
		"vsubpd	      %%ymm6,%%ymm3,%%ymm3		\n\t"\
		"vmovaps	%%ymm0,0x080(%%rsi)	\n\t"\
		"vmovaps	%%ymm2,0x0c0(%%rsi)	\n\t"\
		"vmovaps	%%ymm1,0x0a0(%%rsi)	\n\t"\
		"vmovaps	%%ymm3,0x060(%%rsi)	\n\t"\
		"vaddpd	      %%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm0,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	      %%ymm2,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm1,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm3,%%ymm6,%%ymm6		\n\t"\
		"vmovaps	%%ymm4,     (%%rsi)	\n\t"\
		"vmovaps	%%ymm7,0x040(%%rsi)	\n\t"\
		"vmovaps	%%ymm5,0x020(%%rsi)	\n\t"\
		"vmovaps	%%ymm6,0x0e0(%%rsi)	\n\t"\
		"\n\t"/* SSE2_RADIX4_DIT_0TWIDDLE_2NDOFTWO(add4, add5, add6, add7, r18): */\
		"addq	$0x100,%%rsi			\n\t"\
		"movslq	%[__p04],%%rdi		\n\t"\
		"shlq	$3,%%rdi			\n\t"\
		"addq	%%rdi,%%rax			\n\t"\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	%%rdi,%%rcx			\n\t"\
		"addq	%%rdi,%%rdx			\n\t"\
		"\n\t"\
		"vmovaps	     (%%rax),%%ymm2	\n\t"\
		"vmovaps	     (%%rcx),%%ymm6	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm3	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm7	\n\t"\
		"vmovaps	     (%%rbx),%%ymm0	\n\t"\
		"vmovaps	     (%%rdx),%%ymm4	\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm1	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm5	\n\t"\
		"movq	%[__isrt2],%%rdi	\n\t"\
		"vsubpd	      %%ymm0,%%ymm2,%%ymm2		\n\t"\
		"vsubpd	      %%ymm4,%%ymm6,%%ymm6		\n\t"\
		"vsubpd	      %%ymm1,%%ymm3,%%ymm3		\n\t"\
		"vsubpd	      %%ymm5,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm0,%%ymm0,%%ymm0		\n\t"\
		"vaddpd	      %%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	      %%ymm1,%%ymm1,%%ymm1		\n\t"\
		"vaddpd	      %%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm2,%%ymm0,%%ymm0		\n\t"\
		"vaddpd	      %%ymm6,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	      %%ymm3,%%ymm1,%%ymm1		\n\t"\
		"vaddpd	      %%ymm7,%%ymm5,%%ymm5		\n\t"\
		"vsubpd	      %%ymm4,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	      %%ymm5,%%ymm1,%%ymm1		\n\t"\
		"vmovaps	%%ymm0,0x080(%%rsi)	\n\t"\
		"vmovaps	%%ymm1,0x0a0(%%rsi)	\n\t"\
		"vsubpd	      %%ymm7,%%ymm2,%%ymm2		\n\t"\
		"vaddpd	      %%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vsubpd	      %%ymm6,%%ymm3,%%ymm3		\n\t"\
		"vaddpd	      %%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm0,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm1,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm2,%%ymm7,%%ymm7		\n\t"\
		"vmovaps	%%ymm4,     (%%rsi)	\n\t"\
		"vaddpd	      %%ymm3,%%ymm6,%%ymm6		\n\t"\
		"vmovaps	%%ymm5,0x020(%%rsi)	\n\t"\
		"vmovaps	(%%rdi),%%ymm5	\n\t"/* ISRT2 */\
		"vmovaps	%%ymm3,%%ymm0		\n\t"\
		"vmovaps	%%ymm6,%%ymm1		\n\t"\
		"vsubpd	      %%ymm7,%%ymm3,%%ymm3		\n\t"\
		"vsubpd	      %%ymm2,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm7,%%ymm0,%%ymm0		\n\t"\
		"vaddpd	      %%ymm2,%%ymm1,%%ymm1		\n\t"\
		"vmulpd	      %%ymm5,%%ymm3,%%ymm3		\n\t"\
		"vmulpd	      %%ymm5,%%ymm6,%%ymm6		\n\t"\
		"vmulpd	      %%ymm5,%%ymm0,%%ymm0		\n\t"\
		"vmulpd	      %%ymm5,%%ymm1,%%ymm1		\n\t"\
		"vmovaps	%%ymm3,0x060(%%rsi)	\n\t"\
		"vmovaps	%%ymm6,0x0e0(%%rsi)	\n\t"\
		"vmovaps	%%ymm0,0x040(%%rsi)	\n\t"\
		"vmovaps	%%ymm1,0x0c0(%%rsi)	\n\t"\
		"\n\t"/* SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS(r10,r12,r14,r16,r18,r1A,r1C,r1E): */\
		"vmovaps		-0x100(%%rsi),%%ymm0	\n\t"/* r10 */\
		"vmovaps		-0x080(%%rsi),%%ymm4	\n\t"/* r14 */\
		"vmovaps		-0x0e0(%%rsi),%%ymm1	\n\t"/* r11 */\
		"vmovaps		-0x060(%%rsi),%%ymm5	\n\t"/* r15 */\
		"vmovaps		      (%%rsi),%%ymm2	\n\t"/* r18 */\
		"vmovaps		 0x0a0(%%rsi),%%ymm7	\n\t"/* r1D */\
		"vmovaps		 0x020(%%rsi),%%ymm3	\n\t"/* r19 */\
		"vmovaps		 0x080(%%rsi),%%ymm6	\n\t"/* r1C */\
		"vsubpd	      %%ymm2,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	      %%ymm7,%%ymm4,%%ymm4		\n\t"\
		"vsubpd	      %%ymm3,%%ymm1,%%ymm1		\n\t"\
		"vsubpd	      %%ymm6,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm2,%%ymm2,%%ymm2		\n\t"\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm3,%%ymm3,%%ymm3		\n\t"\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm0,%%ymm2,%%ymm2		\n\t"\
		"vaddpd	      %%ymm4,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm1,%%ymm3,%%ymm3		\n\t"\
		"vaddpd	      %%ymm5,%%ymm6,%%ymm6		\n\t"\
		"vmovaps		%%ymm0,      (%%rsi)\n\t"\
		"vmovaps		%%ymm4, 0x080(%%rsi)\n\t"\
		"vmovaps		%%ymm1, 0x020(%%rsi)\n\t"\
		"vmovaps		%%ymm5,-0x060(%%rsi)\n\t"\
		"vmovaps		%%ymm2,-0x100(%%rsi)\n\t"\
		"vmovaps		%%ymm7,-0x080(%%rsi)\n\t"\
		"vmovaps		%%ymm3,-0x0e0(%%rsi)\n\t"\
		"vmovaps		%%ymm6, 0x0a0(%%rsi)\n\t"\
		"vmovaps		-0x0c0(%%rsi),%%ymm0	\n\t"/* r12 */\
		"vmovaps		-0x040(%%rsi),%%ymm4	\n\t"/* r16 */\
		"vmovaps		-0x0a0(%%rsi),%%ymm1	\n\t"/* r13 */\
		"vmovaps		-0x020(%%rsi),%%ymm5	\n\t"/* r17 */\
		"vmovaps		 0x040(%%rsi),%%ymm2	\n\t"/* r1A */\
		"vmovaps		 0x0e0(%%rsi),%%ymm7	\n\t"/* r1F */\
		"vmovaps		 0x060(%%rsi),%%ymm3	\n\t"/* r1B */\
		"vmovaps		 0x0c0(%%rsi),%%ymm6	\n\t"/* r1E */\
		"vsubpd	      %%ymm2,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	      %%ymm7,%%ymm4,%%ymm4		\n\t"\
		"vsubpd	      %%ymm3,%%ymm1,%%ymm1		\n\t"\
		"vsubpd	      %%ymm6,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm2,%%ymm2,%%ymm2		\n\t"\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm3,%%ymm3,%%ymm3		\n\t"\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm0,%%ymm2,%%ymm2		\n\t"\
		"vaddpd	      %%ymm4,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm1,%%ymm3,%%ymm3		\n\t"\
		"vaddpd	      %%ymm5,%%ymm6,%%ymm6		\n\t"\
		"vmovaps		%%ymm0, 0x040(%%rsi)\n\t"\
		"vmovaps		%%ymm4, 0x0c0(%%rsi)\n\t"\
		"vmovaps		%%ymm1, 0x060(%%rsi)\n\t"\
		"vmovaps		%%ymm5,-0x020(%%rsi)\n\t"\
		"vmovaps		%%ymm2,-0x0c0(%%rsi)\n\t"\
		"vmovaps		%%ymm7,-0x040(%%rsi)\n\t"\
		"vmovaps		%%ymm3,-0x0a0(%%rsi)\n\t"\
		"vmovaps		%%ymm6, 0x0e0(%%rsi)\n\t"\
		"\n\t"\
		"\n\t"/*...Block 3: */\
		"addq	$0x100,%%rsi			\n\t"\
		"movslq	%[__p04],%%rdi		\n\t"\
		"shlq	$3,%%rdi			\n\t"\
		"addq	%%rdi,%%rax			\n\t"\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	%%rdi,%%rcx			\n\t"\
		"addq	%%rdi,%%rdx			\n\t"\
		"\n\t"\
		"\n\t"/* SSE2_RADIX4_DIT_0TWIDDLE(add0, add1, add2, add3, r20): */\
		"vmovaps	     (%%rax),%%ymm2	\n\t"\
		"vmovaps	     (%%rcx),%%ymm6	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm3	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm7	\n\t"\
		"vmovaps	     (%%rbx),%%ymm0	\n\t"\
		"vmovaps	     (%%rdx),%%ymm4	\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm1	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm5	\n\t"\
		"vsubpd	      %%ymm0,%%ymm2,%%ymm2		\n\t"\
		"vsubpd	      %%ymm4,%%ymm6,%%ymm6		\n\t"\
		"vsubpd	      %%ymm1,%%ymm3,%%ymm3		\n\t"\
		"vsubpd	      %%ymm5,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm0,%%ymm0,%%ymm0		\n\t"\
		"vaddpd	      %%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	      %%ymm1,%%ymm1,%%ymm1		\n\t"\
		"vaddpd	      %%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm2,%%ymm0,%%ymm0		\n\t"\
		"vaddpd	      %%ymm6,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	      %%ymm3,%%ymm1,%%ymm1		\n\t"\
		"vaddpd	      %%ymm7,%%ymm5,%%ymm5		\n\t"\
		"vsubpd	      %%ymm4,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	      %%ymm7,%%ymm2,%%ymm2		\n\t"\
		"vsubpd	      %%ymm5,%%ymm1,%%ymm1		\n\t"\
		"vsubpd	      %%ymm6,%%ymm3,%%ymm3		\n\t"\
		"vmovaps	%%ymm0,0x080(%%rsi)	\n\t"\
		"vmovaps	%%ymm2,0x0c0(%%rsi)	\n\t"\
		"vmovaps	%%ymm1,0x0a0(%%rsi)	\n\t"\
		"vmovaps	%%ymm3,0x060(%%rsi)	\n\t"\
		"vaddpd	      %%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm0,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	      %%ymm2,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm1,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm3,%%ymm6,%%ymm6		\n\t"\
		"vmovaps	%%ymm4,     (%%rsi)	\n\t"\
		"vmovaps	%%ymm7,0x040(%%rsi)	\n\t"\
		"vmovaps	%%ymm5,0x020(%%rsi)	\n\t"\
		"vmovaps	%%ymm6,0x0e0(%%rsi)	\n\t"\
		"\n\t"/* SSE2_RADIX4_DIT_0TWIDDLE_2NDOFTWO(add4, add5, add6, add7, r28): */\
		"addq	$0x100,%%rsi			\n\t"\
		"movslq	%[__p04],%%rdi		\n\t"\
		"shlq	$3,%%rdi			\n\t"\
		"addq	%%rdi,%%rax			\n\t"\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	%%rdi,%%rcx			\n\t"\
		"addq	%%rdi,%%rdx			\n\t"\
		"\n\t"\
		"vmovaps	     (%%rax),%%ymm2	\n\t"\
		"vmovaps	     (%%rcx),%%ymm6	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm3	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm7	\n\t"\
		"vmovaps	     (%%rbx),%%ymm0	\n\t"\
		"vmovaps	     (%%rdx),%%ymm4	\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm1	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm5	\n\t"\
		"movq	%[__isrt2],%%rdi	\n\t"\
		"vsubpd	      %%ymm0,%%ymm2,%%ymm2		\n\t"\
		"vsubpd	      %%ymm4,%%ymm6,%%ymm6		\n\t"\
		"vsubpd	      %%ymm1,%%ymm3,%%ymm3		\n\t"\
		"vsubpd	      %%ymm5,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm0,%%ymm0,%%ymm0		\n\t"\
		"vaddpd	      %%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	      %%ymm1,%%ymm1,%%ymm1		\n\t"\
		"vaddpd	      %%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm2,%%ymm0,%%ymm0		\n\t"\
		"vaddpd	      %%ymm6,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	      %%ymm3,%%ymm1,%%ymm1		\n\t"\
		"vaddpd	      %%ymm7,%%ymm5,%%ymm5		\n\t"\
		"vsubpd	      %%ymm4,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	      %%ymm5,%%ymm1,%%ymm1		\n\t"\
		"vmovaps	%%ymm0,0x080(%%rsi)	\n\t"\
		"vmovaps	%%ymm1,0x0a0(%%rsi)	\n\t"\
		"vsubpd	      %%ymm7,%%ymm2,%%ymm2		\n\t"\
		"vaddpd	      %%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vsubpd	      %%ymm6,%%ymm3,%%ymm3		\n\t"\
		"vaddpd	      %%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm0,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm1,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm2,%%ymm7,%%ymm7		\n\t"\
		"vmovaps	%%ymm4,     (%%rsi)	\n\t"\
		"vaddpd	      %%ymm3,%%ymm6,%%ymm6		\n\t"\
		"vmovaps	%%ymm5,0x020(%%rsi)	\n\t"\
		"vmovaps	(%%rdi),%%ymm5	\n\t"/* ISRT2 */\
		"vmovaps	%%ymm3,%%ymm0		\n\t"\
		"vmovaps	%%ymm6,%%ymm1		\n\t"\
		"vsubpd	      %%ymm7,%%ymm3,%%ymm3		\n\t"\
		"vsubpd	      %%ymm2,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm7,%%ymm0,%%ymm0		\n\t"\
		"vaddpd	      %%ymm2,%%ymm1,%%ymm1		\n\t"\
		"vmulpd	      %%ymm5,%%ymm3,%%ymm3		\n\t"\
		"vmulpd	      %%ymm5,%%ymm6,%%ymm6		\n\t"\
		"vmulpd	      %%ymm5,%%ymm0,%%ymm0		\n\t"\
		"vmulpd	      %%ymm5,%%ymm1,%%ymm1		\n\t"\
		"vmovaps	%%ymm3,0x060(%%rsi)	\n\t"\
		"vmovaps	%%ymm6,0x0e0(%%rsi)	\n\t"\
		"vmovaps	%%ymm0,0x040(%%rsi)	\n\t"\
		"vmovaps	%%ymm1,0x0c0(%%rsi)	\n\t"\
		"\n\t"/* SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS(r20,r22,r24,r26,r28,r2A,r2C,r2E): */\
		"vmovaps		-0x100(%%rsi),%%ymm0	\n\t"/* r20 */\
		"vmovaps		-0x080(%%rsi),%%ymm4	\n\t"/* r24 */\
		"vmovaps		-0x0e0(%%rsi),%%ymm1	\n\t"/* r21 */\
		"vmovaps		-0x060(%%rsi),%%ymm5	\n\t"/* r25 */\
		"vmovaps		      (%%rsi),%%ymm2	\n\t"/* r28 */\
		"vmovaps		 0x0a0(%%rsi),%%ymm7	\n\t"/* r2D */\
		"vmovaps		 0x020(%%rsi),%%ymm3	\n\t"/* r29 */\
		"vmovaps		 0x080(%%rsi),%%ymm6	\n\t"/* r2C */\
		"vsubpd	      %%ymm2,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	      %%ymm7,%%ymm4,%%ymm4		\n\t"\
		"vsubpd	      %%ymm3,%%ymm1,%%ymm1		\n\t"\
		"vsubpd	      %%ymm6,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm2,%%ymm2,%%ymm2		\n\t"\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm3,%%ymm3,%%ymm3		\n\t"\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm0,%%ymm2,%%ymm2		\n\t"\
		"vaddpd	      %%ymm4,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm1,%%ymm3,%%ymm3		\n\t"\
		"vaddpd	      %%ymm5,%%ymm6,%%ymm6		\n\t"\
		"vmovaps		%%ymm0,      (%%rsi)\n\t"\
		"vmovaps		%%ymm4, 0x080(%%rsi)\n\t"\
		"vmovaps		%%ymm1, 0x020(%%rsi)\n\t"\
		"vmovaps		%%ymm5,-0x060(%%rsi)\n\t"\
		"vmovaps		%%ymm2,-0x100(%%rsi)\n\t"\
		"vmovaps		%%ymm7,-0x080(%%rsi)\n\t"\
		"vmovaps		%%ymm3,-0x0e0(%%rsi)\n\t"\
		"vmovaps		%%ymm6, 0x0a0(%%rsi)\n\t"\
		"vmovaps		-0x0c0(%%rsi),%%ymm0	\n\t"/* r22 */\
		"vmovaps		-0x040(%%rsi),%%ymm4	\n\t"/* r26 */\
		"vmovaps		-0x0a0(%%rsi),%%ymm1	\n\t"/* r23 */\
		"vmovaps		-0x020(%%rsi),%%ymm5	\n\t"/* r27 */\
		"vmovaps		 0x040(%%rsi),%%ymm2	\n\t"/* r2A */\
		"vmovaps		 0x0e0(%%rsi),%%ymm7	\n\t"/* r2F */\
		"vmovaps		 0x060(%%rsi),%%ymm3	\n\t"/* r2B */\
		"vmovaps		 0x0c0(%%rsi),%%ymm6	\n\t"/* r2E */\
		"vsubpd	      %%ymm2,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	      %%ymm7,%%ymm4,%%ymm4		\n\t"\
		"vsubpd	      %%ymm3,%%ymm1,%%ymm1		\n\t"\
		"vsubpd	      %%ymm6,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm2,%%ymm2,%%ymm2		\n\t"\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm3,%%ymm3,%%ymm3		\n\t"\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm0,%%ymm2,%%ymm2		\n\t"\
		"vaddpd	      %%ymm4,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm1,%%ymm3,%%ymm3		\n\t"\
		"vaddpd	      %%ymm5,%%ymm6,%%ymm6		\n\t"\
		"vmovaps		%%ymm0, 0x040(%%rsi)\n\t"\
		"vmovaps		%%ymm4, 0x0c0(%%rsi)\n\t"\
		"vmovaps		%%ymm1, 0x060(%%rsi)\n\t"\
		"vmovaps		%%ymm5,-0x020(%%rsi)\n\t"\
		"vmovaps		%%ymm2,-0x0c0(%%rsi)\n\t"\
		"vmovaps		%%ymm7,-0x040(%%rsi)\n\t"\
		"vmovaps		%%ymm3,-0x0a0(%%rsi)\n\t"\
		"vmovaps		%%ymm6, 0x0e0(%%rsi)\n\t"\
		"\n\t"\
		"\n\t"/*...Block 4: */\
		"addq	$0x100,%%rsi			\n\t"\
		"movslq	%[__p04],%%rdi		\n\t"\
		"shlq	$3,%%rdi			\n\t"\
		"addq	%%rdi,%%rax			\n\t"\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	%%rdi,%%rcx			\n\t"\
		"addq	%%rdi,%%rdx			\n\t"\
		"\n\t"\
		"\n\t"/* SSE2_RADIX4_DIT_0TWIDDLE(add0, add1, add2, add3, r30): */\
		"vmovaps	     (%%rax),%%ymm2	\n\t"\
		"vmovaps	     (%%rcx),%%ymm6	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm3	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm7	\n\t"\
		"vmovaps	     (%%rbx),%%ymm0	\n\t"\
		"vmovaps	     (%%rdx),%%ymm4	\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm1	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm5	\n\t"\
		"vsubpd	      %%ymm0,%%ymm2,%%ymm2		\n\t"\
		"vsubpd	      %%ymm4,%%ymm6,%%ymm6		\n\t"\
		"vsubpd	      %%ymm1,%%ymm3,%%ymm3		\n\t"\
		"vsubpd	      %%ymm5,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm0,%%ymm0,%%ymm0		\n\t"\
		"vaddpd	      %%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	      %%ymm1,%%ymm1,%%ymm1		\n\t"\
		"vaddpd	      %%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm2,%%ymm0,%%ymm0		\n\t"\
		"vaddpd	      %%ymm6,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	      %%ymm3,%%ymm1,%%ymm1		\n\t"\
		"vaddpd	      %%ymm7,%%ymm5,%%ymm5		\n\t"\
		"vsubpd	      %%ymm4,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	      %%ymm7,%%ymm2,%%ymm2		\n\t"\
		"vsubpd	      %%ymm5,%%ymm1,%%ymm1		\n\t"\
		"vsubpd	      %%ymm6,%%ymm3,%%ymm3		\n\t"\
		"vmovaps	%%ymm0,0x080(%%rsi)	\n\t"\
		"vmovaps	%%ymm2,0x0c0(%%rsi)	\n\t"\
		"vmovaps	%%ymm1,0x0a0(%%rsi)	\n\t"\
		"vmovaps	%%ymm3,0x060(%%rsi)	\n\t"\
		"vaddpd	      %%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm0,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	      %%ymm2,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm1,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm3,%%ymm6,%%ymm6		\n\t"\
		"vmovaps	%%ymm4,     (%%rsi)	\n\t"\
		"vmovaps	%%ymm7,0x040(%%rsi)	\n\t"\
		"vmovaps	%%ymm5,0x020(%%rsi)	\n\t"\
		"vmovaps	%%ymm6,0x0e0(%%rsi)	\n\t"\
		"\n\t"/* SSE2_RADIX4_DIT_0TWIDDLE_2NDOFTWO(add4, add5, add6, add7, r38): */\
		"addq	$0x100,%%rsi			\n\t"\
		"movslq	%[__p04],%%rdi		\n\t"\
		"shlq	$3,%%rdi			\n\t"\
		"addq	%%rdi,%%rax			\n\t"\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	%%rdi,%%rcx			\n\t"\
		"addq	%%rdi,%%rdx			\n\t"\
		"\n\t"\
		"vmovaps	     (%%rax),%%ymm2	\n\t"\
		"vmovaps	     (%%rcx),%%ymm6	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm3	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm7	\n\t"\
		"vmovaps	     (%%rbx),%%ymm0	\n\t"\
		"vmovaps	     (%%rdx),%%ymm4	\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm1	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm5	\n\t"\
		"movq	%[__isrt2],%%rdi	\n\t"\
		"vsubpd	      %%ymm0,%%ymm2,%%ymm2		\n\t"\
		"vsubpd	      %%ymm4,%%ymm6,%%ymm6		\n\t"\
		"vsubpd	      %%ymm1,%%ymm3,%%ymm3		\n\t"\
		"vsubpd	      %%ymm5,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm0,%%ymm0,%%ymm0		\n\t"\
		"vaddpd	      %%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	      %%ymm1,%%ymm1,%%ymm1		\n\t"\
		"vaddpd	      %%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm2,%%ymm0,%%ymm0		\n\t"\
		"vaddpd	      %%ymm6,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	      %%ymm3,%%ymm1,%%ymm1		\n\t"\
		"vaddpd	      %%ymm7,%%ymm5,%%ymm5		\n\t"\
		"vsubpd	      %%ymm4,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	      %%ymm5,%%ymm1,%%ymm1		\n\t"\
		"vmovaps	%%ymm0,0x080(%%rsi)	\n\t"\
		"vmovaps	%%ymm1,0x0a0(%%rsi)	\n\t"\
		"vsubpd	      %%ymm7,%%ymm2,%%ymm2		\n\t"\
		"vaddpd	      %%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vsubpd	      %%ymm6,%%ymm3,%%ymm3		\n\t"\
		"vaddpd	      %%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm0,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm1,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm2,%%ymm7,%%ymm7		\n\t"\
		"vmovaps	%%ymm4,     (%%rsi)	\n\t"\
		"vaddpd	      %%ymm3,%%ymm6,%%ymm6		\n\t"\
		"vmovaps	%%ymm5,0x020(%%rsi)	\n\t"\
		"vmovaps	(%%rdi),%%ymm5	\n\t"/* ISRT2 */\
		"vmovaps	%%ymm3,%%ymm0		\n\t"\
		"vmovaps	%%ymm6,%%ymm1		\n\t"\
		"vsubpd	      %%ymm7,%%ymm3,%%ymm3		\n\t"\
		"vsubpd	      %%ymm2,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm7,%%ymm0,%%ymm0		\n\t"\
		"vaddpd	      %%ymm2,%%ymm1,%%ymm1		\n\t"\
		"vmulpd	      %%ymm5,%%ymm3,%%ymm3		\n\t"\
		"vmulpd	      %%ymm5,%%ymm6,%%ymm6		\n\t"\
		"vmulpd	      %%ymm5,%%ymm0,%%ymm0		\n\t"\
		"vmulpd	      %%ymm5,%%ymm1,%%ymm1		\n\t"\
		"vmovaps	%%ymm3,0x060(%%rsi)	\n\t"\
		"vmovaps	%%ymm6,0x0e0(%%rsi)	\n\t"\
		"vmovaps	%%ymm0,0x040(%%rsi)	\n\t"\
		"vmovaps	%%ymm1,0x0c0(%%rsi)	\n\t"\
		"\n\t"/* SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS(r30,r32,r34,r36,r38,r3A,r3C,r3E): */\
		"vmovaps		-0x100(%%rsi),%%ymm0	\n\t"/* r30 */\
		"vmovaps		-0x080(%%rsi),%%ymm4	\n\t"/* r34 */\
		"vmovaps		-0x0e0(%%rsi),%%ymm1	\n\t"/* r31 */\
		"vmovaps		-0x060(%%rsi),%%ymm5	\n\t"/* r35 */\
		"vmovaps		      (%%rsi),%%ymm2	\n\t"/* r38 */\
		"vmovaps		 0x0a0(%%rsi),%%ymm7	\n\t"/* r3D */\
		"vmovaps		 0x020(%%rsi),%%ymm3	\n\t"/* r39 */\
		"vmovaps		 0x080(%%rsi),%%ymm6	\n\t"/* r3C */\
		"vsubpd	      %%ymm2,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	      %%ymm7,%%ymm4,%%ymm4		\n\t"\
		"vsubpd	      %%ymm3,%%ymm1,%%ymm1		\n\t"\
		"vsubpd	      %%ymm6,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm2,%%ymm2,%%ymm2		\n\t"\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm3,%%ymm3,%%ymm3		\n\t"\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm0,%%ymm2,%%ymm2		\n\t"\
		"vaddpd	      %%ymm4,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm1,%%ymm3,%%ymm3		\n\t"\
		"vaddpd	      %%ymm5,%%ymm6,%%ymm6		\n\t"\
		"vmovaps		%%ymm0,      (%%rsi)\n\t"\
		"vmovaps		%%ymm4, 0x080(%%rsi)\n\t"\
		"vmovaps		%%ymm1, 0x020(%%rsi)\n\t"\
		"vmovaps		%%ymm5,-0x060(%%rsi)\n\t"\
		"vmovaps		%%ymm2,-0x100(%%rsi)\n\t"\
		"vmovaps		%%ymm7,-0x080(%%rsi)\n\t"\
		"vmovaps		%%ymm3,-0x0e0(%%rsi)\n\t"\
		"vmovaps		%%ymm6, 0x0a0(%%rsi)\n\t"\
		"vmovaps		-0x0c0(%%rsi),%%ymm0	\n\t"/* r32 */\
		"vmovaps		-0x040(%%rsi),%%ymm4	\n\t"/* r36 */\
		"vmovaps		-0x0a0(%%rsi),%%ymm1	\n\t"/* r33 */\
		"vmovaps		-0x020(%%rsi),%%ymm5	\n\t"/* r37 */\
		"vmovaps		 0x040(%%rsi),%%ymm2	\n\t"/* r3A */\
		"vmovaps		 0x0e0(%%rsi),%%ymm7	\n\t"/* r3F */\
		"vmovaps		 0x060(%%rsi),%%ymm3	\n\t"/* r3B */\
		"vmovaps		 0x0c0(%%rsi),%%ymm6	\n\t"/* r3E */\
		"vsubpd	      %%ymm2,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	      %%ymm7,%%ymm4,%%ymm4		\n\t"\
		"vsubpd	      %%ymm3,%%ymm1,%%ymm1		\n\t"\
		"vsubpd	      %%ymm6,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm2,%%ymm2,%%ymm2		\n\t"\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm3,%%ymm3,%%ymm3		\n\t"\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm0,%%ymm2,%%ymm2		\n\t"\
		"vaddpd	      %%ymm4,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm1,%%ymm3,%%ymm3		\n\t"\
		"vaddpd	      %%ymm5,%%ymm6,%%ymm6		\n\t"\
		"vmovaps		%%ymm0, 0x040(%%rsi)\n\t"\
		"vmovaps		%%ymm4, 0x0c0(%%rsi)\n\t"\
		"vmovaps		%%ymm1, 0x060(%%rsi)\n\t"\
		"vmovaps		%%ymm5,-0x020(%%rsi)\n\t"\
		"vmovaps		%%ymm2,-0x0c0(%%rsi)\n\t"\
		"vmovaps		%%ymm7,-0x040(%%rsi)\n\t"\
		"vmovaps		%%ymm3,-0x0a0(%%rsi)\n\t"\
		"vmovaps		%%ymm6, 0x0e0(%%rsi)\n\t"\
		"\n\t"\
	"\n\t"/*...and now do eight radix-4 transforms, including the internal twiddle factors:	*/\
		"\n\t"/*...Block 1: r00,r10,r20,r30	*/\
		"\n\t"\
		"movq	%[__isrt2],%%rdi	\n\t"\
		"movq	%[__cc0],%%rsi	\n\t"\
		"movq	%[__r00],%%rax	\n\t"\
		"movq	%%rax,%%rbx		\n\t"\
		"movq	%%rax,%%rcx		\n\t"\
		"movq	%%rax,%%rdx		\n\t"\
		"addq	$0x200,%%rbx	\n\t"\
		"addq	$0x400,%%rcx	\n\t"\
		"addq	$0x600,%%rdx	\n\t"\
		"\n\t"\
		"vmovaps	     (%%rax),%%ymm0		\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1		\n\t"\
		"vmovaps	     (%%rbx),%%ymm2		\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm3		\n\t"\
		"vsubpd	     (%%rbx),%%ymm0,%%ymm0		\n\t"\
		"vsubpd	0x020(%%rbx),%%ymm1,%%ymm1		\n\t"\
		"vaddpd	     (%%rax),%%ymm2,%%ymm2		\n\t"\
		"vaddpd	0x020(%%rax),%%ymm3,%%ymm3		\n\t"\
		"vmovaps	     (%%rcx),%%ymm4		\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5		\n\t"\
		"vmovaps	     (%%rdx),%%ymm6		\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm7		\n\t"\
		"vsubpd	     (%%rdx),%%ymm4,%%ymm4		\n\t"\
		"vsubpd	0x020(%%rdx),%%ymm5,%%ymm5		\n\t"\
		"vaddpd	     (%%rcx),%%ymm6,%%ymm6		\n\t"\
		"vaddpd	0x020(%%rcx),%%ymm7,%%ymm7		\n\t"\
		"\n\t"/* SSE2_RADIX4_DIT_LASTPAIR_IN_PLACE_STRIDE(ymm2,ymm3,ymm0,ymm1,ymm6,ymm7,ymm4,ymm5): swap ymm[01]<->[23], ymm[45]<->[67] */\
		"vaddpd	      %%ymm6,%%ymm2,%%ymm2		\n\t"\
		"vaddpd	      %%ymm7,%%ymm3,%%ymm3		\n\t"\
		"vmovaps	%%ymm2,     (%%rax)	\n\t"\
		"vmovaps	%%ymm3,0x020(%%rax)	\n\t"\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vsubpd	      %%ymm6,%%ymm2,%%ymm2		\n\t"\
		"vsubpd	      %%ymm7,%%ymm3,%%ymm3		\n\t"\
		"vmovaps	%%ymm2,     (%%rcx)	\n\t"\
		"vmovaps	%%ymm3,0x020(%%rcx)	\n\t"\
		"vaddpd	      %%ymm5,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	      %%ymm4,%%ymm1,%%ymm1		\n\t"\
		"vmovaps	%%ymm0,     (%%rbx)	\n\t"\
		"vmovaps	%%ymm1,0x020(%%rbx)	\n\t"\
		"vaddpd	      %%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vsubpd	      %%ymm5,%%ymm0,%%ymm0		\n\t"\
		"vaddpd	      %%ymm4,%%ymm1,%%ymm1		\n\t"\
		"vmovaps	%%ymm0,     (%%rdx)	\n\t"\
		"vmovaps	%%ymm1,0x020(%%rdx)	\n\t"\
		"\n\t"\
		"\n\t"/*...Block 5: r08,r18,r28,r38	*/\
		"\n\t"\
		"addq	$0x100,%%rax		\n\t"/* r08 */\
		"addq	$0x100,%%rbx		\n\t"\
		"addq	$0x100,%%rcx		\n\t"\
		"addq	$0x100,%%rdx		\n\t"\
		"vmovaps	(%%rdi),%%ymm2	\n\t"/* isrt2 */\
		"\n\t"\
		"vmovaps	     (%%rcx),%%ymm4	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5	\n\t"\
		"vmovaps	     (%%rdx),%%ymm0	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm1	\n\t"\
		"vaddpd	0x020(%%rcx),%%ymm4,%%ymm4		\n\t"\
		"vsubpd	     (%%rcx),%%ymm5,%%ymm5		\n\t"\
		"vsubpd	0x020(%%rdx),%%ymm0,%%ymm0		\n\t"\
		"vaddpd	     (%%rdx),%%ymm1,%%ymm1		\n\t"\
		"vmulpd	      %%ymm2,%%ymm4,%%ymm4		\n\t"\
		"vmulpd	      %%ymm2,%%ymm5,%%ymm5		\n\t"\
		"vmulpd	      %%ymm2,%%ymm0,%%ymm0		\n\t"\
		"vmulpd	      %%ymm2,%%ymm1,%%ymm1		\n\t"\
		"vmovaps	%%ymm4,%%ymm6		\n\t"\
		"vmovaps	%%ymm5,%%ymm7		\n\t"\
		"vsubpd	      %%ymm0,%%ymm4,%%ymm4		\n\t"\
		"vsubpd	      %%ymm1,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm0,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm1,%%ymm7,%%ymm7		\n\t"\
		"vmovaps	     (%%rax),%%ymm0	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1	\n\t"\
		"vmovaps	     (%%rbx),%%ymm2	\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm3	\n\t"\
		"vsubpd	0x020(%%rbx),%%ymm0,%%ymm0		\n\t"\
		"vsubpd	     (%%rbx),%%ymm1,%%ymm1		\n\t"\
		"vaddpd	     (%%rax),%%ymm3,%%ymm3		\n\t"\
		"vaddpd	0x020(%%rax),%%ymm2,%%ymm2		\n\t"\
		"\n\t"\
		"\n\t"/* SSE2_RADIX4_DIT_LASTPAIR_IN_PLACE_STRIDE(ymm3,ymm1,ymm0,ymm2,ymm4,ymm5,ymm6,ymm7): swap ymm0123<->3102 */\
		"vaddpd	      %%ymm4,%%ymm3,%%ymm3		\n\t"\
		"vaddpd	      %%ymm5,%%ymm1,%%ymm1		\n\t"\
		"vmovaps	%%ymm3,     (%%rax)	\n\t"\
		"vmovaps	%%ymm1,0x020(%%rax)	\n\t"\
		"vaddpd	      %%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	      %%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vsubpd	      %%ymm4,%%ymm3,%%ymm3		\n\t"\
		"vsubpd	      %%ymm5,%%ymm1,%%ymm1		\n\t"\
		"vmovaps	%%ymm3,     (%%rcx)	\n\t"\
		"vmovaps	%%ymm1,0x020(%%rcx)	\n\t"\
		"vaddpd	      %%ymm7,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	      %%ymm6,%%ymm2,%%ymm2		\n\t"\
		"vmovaps	%%ymm0,     (%%rbx)	\n\t"\
		"vmovaps	%%ymm2,0x020(%%rbx)	\n\t"\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vsubpd	      %%ymm7,%%ymm0,%%ymm0		\n\t"\
		"vaddpd	      %%ymm6,%%ymm2,%%ymm2		\n\t"\
		"vmovaps	%%ymm0,     (%%rdx)	\n\t"\
		"vmovaps	%%ymm2,0x020(%%rdx)	\n\t"\
		"\n\t"\
		"\n\t"/*...Block 3: r04,r14,r24,r34	*/\
		"\n\t"\
		"subq	$0x080,%%rax 	\n\t"/* r04 */\
		"subq	$0x080,%%rbx		\n\t"\
		"subq	$0x080,%%rcx		\n\t"\
		"subq	$0x080,%%rdx		\n\t"\
		"\n\t"\
		"vmovaps	     (%%rcx),%%ymm4	\n\t"\
		"vmovaps	     (%%rdx),%%ymm0	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm1	\n\t"\
		"vmovaps	     (%%rcx),%%ymm6	\n\t"\
		"vmovaps	     (%%rdx),%%ymm2	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm7	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm3	\n\t"\
		"\n\t"\
		"vmulpd	     (%%rsi),%%ymm4,%%ymm4		\n\t"\
		"vmulpd	0x020(%%rsi),%%ymm0,%%ymm0		\n\t"\
		"vmulpd	     (%%rsi),%%ymm5,%%ymm5		\n\t"\
		"vmulpd	0x020(%%rsi),%%ymm1,%%ymm1		\n\t"\
		"vmulpd	0x020(%%rsi),%%ymm6,%%ymm6		\n\t"\
		"vmulpd	     (%%rsi),%%ymm2,%%ymm2		\n\t"\
		"vmulpd	0x020(%%rsi),%%ymm7,%%ymm7		\n\t"\
		"vmulpd	     (%%rsi),%%ymm3,%%ymm3		\n\t"\
		"vsubpd	      %%ymm6,%%ymm5,%%ymm5		\n\t"\
		"vsubpd	      %%ymm2,%%ymm1,%%ymm1		\n\t"\
		"vaddpd	      %%ymm7,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	      %%ymm3,%%ymm0,%%ymm0		\n\t"\
		"vmovaps	%%ymm5,%%ymm7		\n\t"\
		"vmovaps	%%ymm4,%%ymm6		\n\t"\
		"vaddpd	      %%ymm0,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	      %%ymm1,%%ymm5,%%ymm5		\n\t"\
		"vsubpd	      %%ymm0,%%ymm6,%%ymm6		\n\t"\
		"vsubpd	      %%ymm1,%%ymm7,%%ymm7		\n\t"\
		"vmovaps	     (%%rbx),%%ymm2	\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm3	\n\t"\
		"vmovaps	     (%%rax),%%ymm0	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1	\n\t"\
		"vaddpd	0x020(%%rbx),%%ymm2,%%ymm2		\n\t"\
		"vsubpd	     (%%rbx),%%ymm3,%%ymm3		\n\t"\
		"vmulpd	     (%%rdi),%%ymm2,%%ymm2		\n\t"\
		"vmulpd	     (%%rdi),%%ymm3,%%ymm3		\n\t"\
		"vsubpd	      %%ymm2,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	      %%ymm3,%%ymm1,%%ymm1		\n\t"\
		"vaddpd	      %%ymm2,%%ymm2,%%ymm2		\n\t"\
		"vaddpd	      %%ymm3,%%ymm3,%%ymm3		\n\t"\
		"vaddpd	      %%ymm0,%%ymm2,%%ymm2		\n\t"\
		"vaddpd	      %%ymm1,%%ymm3,%%ymm3		\n\t"\
		"\n\t"\
		"\n\t"/* SSE2_RADIX4_DIT_LASTPAIR_IN_PLACE_STRIDE(ymm2,ymm3,ymm0,ymm1,ymm4,ymm5,ymm6,ymm7): swap ymm[01]<->[23] */\
		"vaddpd	      %%ymm4,%%ymm2,%%ymm2		\n\t"\
		"vaddpd	      %%ymm5,%%ymm3,%%ymm3		\n\t"\
		"vmovaps	%%ymm2,     (%%rax)	\n\t"\
		"vmovaps	%%ymm3,0x020(%%rax)	\n\t"\
		"vaddpd	      %%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	      %%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vsubpd	      %%ymm4,%%ymm2,%%ymm2		\n\t"\
		"vsubpd	      %%ymm5,%%ymm3,%%ymm3		\n\t"\
		"vmovaps	%%ymm2,     (%%rcx)	\n\t"\
		"vmovaps	%%ymm3,0x020(%%rcx)	\n\t"\
		"vaddpd	      %%ymm7,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	      %%ymm6,%%ymm1,%%ymm1		\n\t"\
		"vmovaps	%%ymm0,     (%%rbx)	\n\t"\
		"vmovaps	%%ymm1,0x020(%%rbx)	\n\t"\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vsubpd	      %%ymm7,%%ymm0,%%ymm0		\n\t"\
		"vaddpd	      %%ymm6,%%ymm1,%%ymm1		\n\t"\
		"vmovaps	%%ymm0,     (%%rdx)	\n\t"\
		"vmovaps	%%ymm1,0x020(%%rdx)	\n\t"\
		"\n\t"\
		"\n\t"/*...Block 7: r0C,r1C,r2C,r3C	*/\
		"\n\t"\
		"addq	$0x100,%%rax 	\n\t"/* r0C */\
		"addq	$0x100,%%rbx		\n\t"\
		"addq	$0x100,%%rcx		\n\t"\
		"addq	$0x100,%%rdx		\n\t"\
		"\n\t"\
		"vmovaps	     (%%rcx),%%ymm4	\n\t"\
		"vmovaps	     (%%rdx),%%ymm0	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm1	\n\t"\
		"vmovaps	     (%%rcx),%%ymm6	\n\t"\
		"vmovaps	     (%%rdx),%%ymm2	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm7	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm3	\n\t"\
		"\n\t"\
		"vmulpd	0x020(%%rsi),%%ymm4,%%ymm4		\n\t"\
		"vmulpd	     (%%rsi),%%ymm0,%%ymm0		\n\t"\
		"vmulpd	0x020(%%rsi),%%ymm5,%%ymm5		\n\t"\
		"vmulpd	     (%%rsi),%%ymm1,%%ymm1		\n\t"\
		"vmulpd	     (%%rsi),%%ymm6,%%ymm6		\n\t"\
		"vmulpd	0x020(%%rsi),%%ymm2,%%ymm2		\n\t"\
		"vmulpd	     (%%rsi),%%ymm7,%%ymm7		\n\t"\
		"vmulpd	0x020(%%rsi),%%ymm3,%%ymm3		\n\t"\
		"vsubpd	      %%ymm6,%%ymm5,%%ymm5		\n\t"\
		"vsubpd	      %%ymm2,%%ymm1,%%ymm1		\n\t"\
		"vaddpd	      %%ymm7,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	      %%ymm3,%%ymm0,%%ymm0		\n\t"\
		"vmovaps	%%ymm5,%%ymm7		\n\t"\
		"vmovaps	%%ymm4,%%ymm6		\n\t"\
		"vaddpd	      %%ymm0,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	      %%ymm1,%%ymm5,%%ymm5		\n\t"\
		"vsubpd	      %%ymm0,%%ymm6,%%ymm6		\n\t"\
		"vsubpd	      %%ymm1,%%ymm7,%%ymm7		\n\t"\
		"vmovaps	     (%%rbx),%%ymm2	\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm3	\n\t"\
		"vmovaps	     (%%rax),%%ymm0	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1	\n\t"\
		"vsubpd	0x020(%%rbx),%%ymm2,%%ymm2		\n\t"\
		"vaddpd	     (%%rbx),%%ymm3,%%ymm3		\n\t"\
		"vmulpd	     (%%rdi),%%ymm2,%%ymm2		\n\t"\
		"vmulpd	     (%%rdi),%%ymm3,%%ymm3		\n\t"\
		"vsubpd	      %%ymm2,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	      %%ymm3,%%ymm1,%%ymm1		\n\t"\
		"vaddpd	      %%ymm2,%%ymm2,%%ymm2		\n\t"\
		"vaddpd	      %%ymm3,%%ymm3,%%ymm3		\n\t"\
		"vaddpd	      %%ymm0,%%ymm2,%%ymm2		\n\t"\
		"vaddpd	      %%ymm1,%%ymm3,%%ymm3		\n\t"\
		"\n\t"\
		"\n\t"/* SSE2_RADIX4_DIT_LASTPAIR_IN_PLACE_STRIDE(ymm0,ymm1,ymm2,ymm3,ymm6,ymm7,ymm4,ymm5): swap ymm[45]<->[67] */\
		"vaddpd	      %%ymm6,%%ymm0,%%ymm0		\n\t"\
		"vaddpd	      %%ymm7,%%ymm1,%%ymm1		\n\t"\
		"vmovaps	%%ymm0,     (%%rax)	\n\t"\
		"vmovaps	%%ymm1,0x020(%%rax)	\n\t"\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vsubpd	      %%ymm6,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	      %%ymm7,%%ymm1,%%ymm1		\n\t"\
		"vmovaps	%%ymm0,     (%%rcx)	\n\t"\
		"vmovaps	%%ymm1,0x020(%%rcx)	\n\t"\
		"vaddpd	      %%ymm5,%%ymm2,%%ymm2		\n\t"\
		"vsubpd	      %%ymm4,%%ymm3,%%ymm3		\n\t"\
		"vmovaps	%%ymm2,     (%%rbx)	\n\t"\
		"vmovaps	%%ymm3,0x020(%%rbx)	\n\t"\
		"vaddpd	      %%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vsubpd	      %%ymm5,%%ymm2,%%ymm2		\n\t"\
		"vaddpd	      %%ymm4,%%ymm3,%%ymm3		\n\t"\
		"vmovaps	%%ymm2,     (%%rdx)	\n\t"\
		"vmovaps	%%ymm3,0x020(%%rdx)	\n\t"\
		"\n\t"\
		"\n\t"/*...Block 2: r02,r12,r22,r32	*/\
		"\n\t"\
		"subq	$0x140,%%rax 	\n\t"/* r02 */\
		"subq	$0x140,%%rbx			\n\t"\
		"subq	$0x140,%%rcx			\n\t"\
		"subq	$0x140,%%rdx			\n\t"\
		"addq	$0x060,%%rdi \n\t"/* cc1 */\
		"addq	$0x080,%%rsi \n\t"/* cc3 */\
		"\n\t"\
		"vmovaps	     (%%rcx),%%ymm4	\n\t"\
		"vmovaps	     (%%rdx),%%ymm0	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm1	\n\t"\
		"vmovaps	     (%%rcx),%%ymm6	\n\t"\
		"vmovaps	     (%%rdx),%%ymm2	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm7	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm3	\n\t"\
		"\n\t"\
		"vmulpd	     (%%rdi),%%ymm4,%%ymm4		\n\t"\
		"vmulpd	     (%%rsi),%%ymm0,%%ymm0		\n\t"\
		"vmulpd	     (%%rdi),%%ymm5,%%ymm5		\n\t"\
		"vmulpd	     (%%rsi),%%ymm1,%%ymm1		\n\t"\
		"vmulpd	0x020(%%rdi),%%ymm6,%%ymm6		\n\t"\
		"vmulpd	0x020(%%rsi),%%ymm2,%%ymm2		\n\t"\
		"vmulpd	0x020(%%rdi),%%ymm7,%%ymm7		\n\t"\
		"vmulpd	0x020(%%rsi),%%ymm3,%%ymm3		\n\t"\
		"vsubpd	      %%ymm6,%%ymm5,%%ymm5		\n\t"\
		"vsubpd	      %%ymm2,%%ymm1,%%ymm1		\n\t"\
		"vaddpd	      %%ymm7,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	      %%ymm3,%%ymm0,%%ymm0		\n\t"\
		"vmovaps	%%ymm5,%%ymm7		\n\t"\
		"vmovaps	%%ymm4,%%ymm6		\n\t"\
		"vsubpd	      %%ymm0,%%ymm4,%%ymm4		\n\t"\
		"vsubpd	      %%ymm1,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm0,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm1,%%ymm7,%%ymm7		\n\t"\
		"\n\t"\
		"subq	$0x080,%%rsi \n\t"/* cc0 */\
		"vmovaps	     (%%rbx),%%ymm2	\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm3	\n\t"\
		"vmovaps	     (%%rbx),%%ymm0	\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm1	\n\t"\
		"vmulpd	     (%%rsi),%%ymm2,%%ymm2		\n\t"\
		"vmulpd	0x020(%%rsi),%%ymm1,%%ymm1		\n\t"\
		"vmulpd	     (%%rsi),%%ymm3,%%ymm3		\n\t"\
		"vmulpd	0x020(%%rsi),%%ymm0,%%ymm0		\n\t"\
		"vaddpd	      %%ymm1,%%ymm2,%%ymm2		\n\t"\
		"vsubpd	      %%ymm0,%%ymm3,%%ymm3		\n\t"\
		"vmovaps	     (%%rax),%%ymm0	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1	\n\t"\
		"vsubpd	      %%ymm2,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	      %%ymm3,%%ymm1,%%ymm1		\n\t"\
		"vaddpd	      %%ymm2,%%ymm2,%%ymm2		\n\t"\
		"vaddpd	      %%ymm3,%%ymm3,%%ymm3		\n\t"\
		"vaddpd	      %%ymm0,%%ymm2,%%ymm2		\n\t"\
		"vaddpd	      %%ymm1,%%ymm3,%%ymm3		\n\t"\
		"\n\t"\
		"\n\t"/* SSE2_RADIX4_DIT_LASTPAIR_IN_PLACE_STRIDE(ymm2,ymm3,ymm0,ymm1,ymm6,ymm7,ymm4,ymm5): swap ymm[01]<->[23], ymm[45]<->[67] */\
		"vaddpd	      %%ymm6,%%ymm2,%%ymm2		\n\t"\
		"vaddpd	      %%ymm7,%%ymm3,%%ymm3		\n\t"\
		"vmovaps	%%ymm2,     (%%rax)	\n\t"\
		"vmovaps	%%ymm3,0x020(%%rax)	\n\t"\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vsubpd	      %%ymm6,%%ymm2,%%ymm2		\n\t"\
		"vsubpd	      %%ymm7,%%ymm3,%%ymm3		\n\t"\
		"vmovaps	%%ymm2,     (%%rcx)	\n\t"\
		"vmovaps	%%ymm3,0x020(%%rcx)	\n\t"\
		"vaddpd	      %%ymm5,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	      %%ymm4,%%ymm1,%%ymm1		\n\t"\
		"vmovaps	%%ymm0,     (%%rbx)	\n\t"\
		"vmovaps	%%ymm1,0x020(%%rbx)	\n\t"\
		"vaddpd	      %%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vsubpd	      %%ymm5,%%ymm0,%%ymm0		\n\t"\
		"vaddpd	      %%ymm4,%%ymm1,%%ymm1		\n\t"\
		"vmovaps	%%ymm0,     (%%rdx)	\n\t"\
		"vmovaps	%%ymm1,0x020(%%rdx)	\n\t"\
		"\n\t"\
		"\n\t"/*...Block 6: r0A,r1A,r2A,r3A	*/\
		"\n\t"\
		"addq	$0x100,%%rax 	\n\t"/* r0A */\
		"addq	$0x100,%%rbx			\n\t"\
		"addq	$0x100,%%rcx			\n\t"\
		"addq	$0x100,%%rdx			\n\t"\
		"addq	$0x080,%%rsi \n\t"/* cc3 */\
		"\n\t"\
		"vmovaps	     (%%rcx),%%ymm4	\n\t"\
		"vmovaps	     (%%rdx),%%ymm0	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm1	\n\t"\
		"vmovaps	     (%%rcx),%%ymm6	\n\t"\
		"vmovaps	     (%%rdx),%%ymm2	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm7	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm3	\n\t"\
		"\n\t"\
		"vmulpd	0x020(%%rsi),%%ymm4,%%ymm4		\n\t"\
		"vmulpd	     (%%rdi),%%ymm0,%%ymm0		\n\t"\
		"vmulpd	0x020(%%rsi),%%ymm5,%%ymm5		\n\t"\
		"vmulpd	     (%%rdi),%%ymm1,%%ymm1		\n\t"\
		"vmulpd	     (%%rsi),%%ymm6,%%ymm6		\n\t"\
		"vmulpd	0x020(%%rdi),%%ymm2,%%ymm2		\n\t"\
		"vmulpd	     (%%rsi),%%ymm7,%%ymm7		\n\t"\
		"vmulpd	0x020(%%rdi),%%ymm3,%%ymm3		\n\t"\
		"vsubpd	      %%ymm6,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm2,%%ymm1,%%ymm1		\n\t"\
		"vaddpd	      %%ymm7,%%ymm4,%%ymm4		\n\t"\
		"vsubpd	      %%ymm3,%%ymm0,%%ymm0		\n\t"\
		"vmovaps	%%ymm5,%%ymm7		\n\t"\
		"vmovaps	%%ymm4,%%ymm6		\n\t"\
		"vaddpd	      %%ymm0,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	      %%ymm1,%%ymm5,%%ymm5		\n\t"\
		"vsubpd	      %%ymm0,%%ymm6,%%ymm6		\n\t"\
		"vsubpd	      %%ymm1,%%ymm7,%%ymm7		\n\t"\
		"subq	$0x080,%%rsi \n\t"/* cc0 */\
		"vmovaps	     (%%rbx),%%ymm2	\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm3	\n\t"\
		"vmovaps	     (%%rbx),%%ymm0	\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm1	\n\t"\
		"vmulpd	0x020(%%rsi),%%ymm2,%%ymm2		\n\t"\
		"vmulpd	     (%%rsi),%%ymm1,%%ymm1		\n\t"\
		"vmulpd	0x020(%%rsi),%%ymm3,%%ymm3		\n\t"\
		"vmulpd	     (%%rsi),%%ymm0,%%ymm0		\n\t"\
		"vsubpd	      %%ymm1,%%ymm2,%%ymm2		\n\t"\
		"vaddpd	      %%ymm0,%%ymm3,%%ymm3		\n\t"\
		"vmovaps	     (%%rax),%%ymm0	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1	\n\t"\
		"vsubpd	      %%ymm2,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	      %%ymm3,%%ymm1,%%ymm1		\n\t"\
		"vaddpd	      %%ymm2,%%ymm2,%%ymm2		\n\t"\
		"vaddpd	      %%ymm3,%%ymm3,%%ymm3		\n\t"\
		"vaddpd	      %%ymm0,%%ymm2,%%ymm2		\n\t"\
		"vaddpd	      %%ymm1,%%ymm3,%%ymm3		\n\t"\
		"\n\t"\
		"\n\t"/* SSE2_RADIX4_DIT_LASTPAIR_IN_PLACE_STRIDE(ymm0,ymm1,ymm2,ymm3,ymm6,ymm7,ymm4,ymm5): swap ymm[45]<->[67] */\
		"vaddpd	      %%ymm6,%%ymm0,%%ymm0		\n\t"\
		"vaddpd	      %%ymm7,%%ymm1,%%ymm1		\n\t"\
		"vmovaps	%%ymm0,     (%%rax)	\n\t"\
		"vmovaps	%%ymm1,0x020(%%rax)	\n\t"\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vsubpd	      %%ymm6,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	      %%ymm7,%%ymm1,%%ymm1		\n\t"\
		"vmovaps	%%ymm0,     (%%rcx)	\n\t"\
		"vmovaps	%%ymm1,0x020(%%rcx)	\n\t"\
		"vaddpd	      %%ymm5,%%ymm2,%%ymm2		\n\t"\
		"vsubpd	      %%ymm4,%%ymm3,%%ymm3		\n\t"\
		"vmovaps	%%ymm2,     (%%rbx)	\n\t"\
		"vmovaps	%%ymm3,0x020(%%rbx)	\n\t"\
		"vaddpd	      %%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vsubpd	      %%ymm5,%%ymm2,%%ymm2		\n\t"\
		"vaddpd	      %%ymm4,%%ymm3,%%ymm3		\n\t"\
		"vmovaps	%%ymm2,     (%%rdx)	\n\t"\
		"vmovaps	%%ymm3,0x020(%%rdx)	\n\t"\
		"\n\t"\
		"\n\t"/*...Block 4: r06,r16,r26,r36	*/\
		"\n\t"\
		"subq	$0x080,%%rax 	\n\t"/* r06 */\
		"subq	$0x080,%%rbx			\n\t"\
		"subq	$0x080,%%rcx			\n\t"\
		"subq	$0x080,%%rdx			\n\t"\
		"addq	$0x080,%%rsi \n\t"/* cc3 */\
		"\n\t"\
		"vmovaps	     (%%rcx),%%ymm4	\n\t"\
		"vmovaps	     (%%rdx),%%ymm0	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm1	\n\t"\
		"vmovaps	     (%%rcx),%%ymm6	\n\t"\
		"vmovaps	     (%%rdx),%%ymm2	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm7	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm3	\n\t"\
		"\n\t"\
		"vmulpd	     (%%rsi),%%ymm4,%%ymm4		\n\t"\
		"vmulpd	0x020(%%rdi),%%ymm0,%%ymm0		\n\t"\
		"vmulpd	     (%%rsi),%%ymm5,%%ymm5		\n\t"\
		"vmulpd	0x020(%%rdi),%%ymm1,%%ymm1		\n\t"\
		"vmulpd	0x020(%%rsi),%%ymm6,%%ymm6		\n\t"\
		"vmulpd	     (%%rdi),%%ymm2,%%ymm2		\n\t"\
		"vmulpd	0x020(%%rsi),%%ymm7,%%ymm7		\n\t"\
		"vmulpd	     (%%rdi),%%ymm3,%%ymm3		\n\t"\
		"vsubpd	      %%ymm6,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm2,%%ymm1,%%ymm1		\n\t"\
		"vaddpd	      %%ymm7,%%ymm4,%%ymm4		\n\t"\
		"vsubpd	      %%ymm3,%%ymm0,%%ymm0		\n\t"\
		"vmovaps	%%ymm5,%%ymm7		\n\t"\
		"vmovaps	%%ymm4,%%ymm6		\n\t"\
		"vaddpd	      %%ymm0,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	      %%ymm1,%%ymm5,%%ymm5		\n\t"\
		"vsubpd	      %%ymm0,%%ymm6,%%ymm6		\n\t"\
		"vsubpd	      %%ymm1,%%ymm7,%%ymm7		\n\t"\
		"subq	$0x080,%%rsi \n\t"/* cc0 */\
		"vmovaps	     (%%rbx),%%ymm2	\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm3	\n\t"\
		"vmovaps	     (%%rbx),%%ymm0	\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm1	\n\t"\
		"vmulpd	0x020(%%rsi),%%ymm2,%%ymm2		\n\t"\
		"vmulpd	     (%%rsi),%%ymm1,%%ymm1		\n\t"\
		"vmulpd	0x020(%%rsi),%%ymm3,%%ymm3		\n\t"\
		"vmulpd	     (%%rsi),%%ymm0,%%ymm0		\n\t"\
		"vaddpd	      %%ymm1,%%ymm2,%%ymm2		\n\t"\
		"vsubpd	      %%ymm0,%%ymm3,%%ymm3		\n\t"\
		"vmovaps	     (%%rax),%%ymm0	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1	\n\t"\
		"vsubpd	      %%ymm2,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	      %%ymm3,%%ymm1,%%ymm1		\n\t"\
		"vaddpd	      %%ymm2,%%ymm2,%%ymm2		\n\t"\
		"vaddpd	      %%ymm3,%%ymm3,%%ymm3		\n\t"\
		"vaddpd	      %%ymm0,%%ymm2,%%ymm2		\n\t"\
		"vaddpd	      %%ymm1,%%ymm3,%%ymm3		\n\t"\
		"\n\t"\
		"\n\t"/* SSE2_RADIX4_DIT_LASTPAIR_IN_PLACE_STRIDE(ymm2,ymm3,ymm0,ymm1,ymm6,ymm7,ymm4,ymm5): swap ymm[01]<->[23], ymm[45]<->[67] */\
		"vaddpd	      %%ymm6,%%ymm2,%%ymm2		\n\t"\
		"vaddpd	      %%ymm7,%%ymm3,%%ymm3		\n\t"\
		"vmovaps	%%ymm2,     (%%rax)	\n\t"\
		"vmovaps	%%ymm3,0x020(%%rax)	\n\t"\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vsubpd	      %%ymm6,%%ymm2,%%ymm2		\n\t"\
		"vsubpd	      %%ymm7,%%ymm3,%%ymm3		\n\t"\
		"vmovaps	%%ymm2,     (%%rcx)	\n\t"\
		"vmovaps	%%ymm3,0x020(%%rcx)	\n\t"\
		"vaddpd	      %%ymm5,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	      %%ymm4,%%ymm1,%%ymm1		\n\t"\
		"vmovaps	%%ymm0,     (%%rbx)	\n\t"\
		"vmovaps	%%ymm1,0x020(%%rbx)	\n\t"\
		"vaddpd	      %%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vsubpd	      %%ymm5,%%ymm0,%%ymm0		\n\t"\
		"vaddpd	      %%ymm4,%%ymm1,%%ymm1		\n\t"\
		"vmovaps	%%ymm0,     (%%rdx)	\n\t"\
		"vmovaps	%%ymm1,0x020(%%rdx)	\n\t"\
		"\n\t"\
		"\n\t"/*...Block 8: r0E,r1E,r2E,r3E	*/\
		"\n\t"\
		"addq	$0x100,%%rax 	\n\t"/* r0E */\
		"addq	$0x100,%%rbx		\n\t"\
		"addq	$0x100,%%rcx		\n\t"\
		"addq	$0x100,%%rdx		\n\t"\
		"addq	$0x080,%%rsi \n\t"/* cc3 */\
		"\n\t"\
		"vmovaps	     (%%rcx),%%ymm4	\n\t"\
		"vmovaps	     (%%rdx),%%ymm0	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm1	\n\t"\
		"vmovaps	     (%%rcx),%%ymm6	\n\t"\
		"vmovaps	     (%%rdx),%%ymm2	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm7	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm3	\n\t"\
		"\n\t"\
		"vmulpd	0x020(%%rdi),%%ymm4,%%ymm4		\n\t"\
		"vmulpd	0x020(%%rsi),%%ymm0,%%ymm0		\n\t"\
		"vmulpd	0x020(%%rdi),%%ymm5,%%ymm5		\n\t"\
		"vmulpd	0x020(%%rsi),%%ymm1,%%ymm1		\n\t"\
		"vmulpd	     (%%rdi),%%ymm6,%%ymm6		\n\t"\
		"vmulpd	     (%%rsi),%%ymm2,%%ymm2		\n\t"\
		"vmulpd	     (%%rdi),%%ymm7,%%ymm7		\n\t"\
		"vmulpd	     (%%rsi),%%ymm3,%%ymm3		\n\t"\
		"vsubpd	      %%ymm6,%%ymm5,%%ymm5		\n\t"\
		"vsubpd	      %%ymm2,%%ymm1,%%ymm1		\n\t"\
		"vaddpd	      %%ymm7,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	      %%ymm3,%%ymm0,%%ymm0		\n\t"\
		"vmovaps	%%ymm5,%%ymm7		\n\t"\
		"vmovaps	%%ymm4,%%ymm6		\n\t"\
		"vaddpd	      %%ymm0,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	      %%ymm1,%%ymm5,%%ymm5		\n\t"\
		"vsubpd	      %%ymm0,%%ymm6,%%ymm6		\n\t"\
		"vsubpd	      %%ymm1,%%ymm7,%%ymm7		\n\t"\
		"subq	$0x080,%%rsi \n\t"/* cc0 */\
		"vmovaps	     (%%rbx),%%ymm2	\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm3	\n\t"\
		"vmovaps	     (%%rbx),%%ymm0	\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm1	\n\t"\
		"vmulpd	     (%%rsi),%%ymm2,%%ymm2		\n\t"\
		"vmulpd	0x020(%%rsi),%%ymm1,%%ymm1		\n\t"\
		"vmulpd	     (%%rsi),%%ymm3,%%ymm3		\n\t"\
		"vmulpd	0x020(%%rsi),%%ymm0,%%ymm0		\n\t"\
		"vsubpd	      %%ymm1,%%ymm2,%%ymm2		\n\t"\
		"vaddpd	      %%ymm0,%%ymm3,%%ymm3		\n\t"\
		"vmovaps	     (%%rax),%%ymm0	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1	\n\t"\
		"vsubpd	      %%ymm2,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	      %%ymm3,%%ymm1,%%ymm1		\n\t"\
		"vaddpd	      %%ymm2,%%ymm2,%%ymm2		\n\t"\
		"vaddpd	      %%ymm3,%%ymm3,%%ymm3		\n\t"\
		"vaddpd	      %%ymm0,%%ymm2,%%ymm2		\n\t"\
		"vaddpd	      %%ymm1,%%ymm3,%%ymm3		\n\t"\
		"\n\t"\
		"\n\t"/* SSE2_RADIX4_DIT_LASTPAIR_IN_PLACE_STRIDE(ymm0,ymm1,ymm2,ymm3,ymm6,ymm7,ymm4,ymm5): swap ymm[45]<->[67] */\
		"vaddpd	      %%ymm6,%%ymm0,%%ymm0		\n\t"\
		"vaddpd	      %%ymm7,%%ymm1,%%ymm1		\n\t"\
		"vmovaps	%%ymm0,     (%%rax)	\n\t"\
		"vmovaps	%%ymm1,0x020(%%rax)	\n\t"\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vsubpd	      %%ymm6,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	      %%ymm7,%%ymm1,%%ymm1		\n\t"\
		"vmovaps	%%ymm0,     (%%rcx)	\n\t"\
		"vmovaps	%%ymm1,0x020(%%rcx)	\n\t"\
		"vaddpd	      %%ymm5,%%ymm2,%%ymm2		\n\t"\
		"vsubpd	      %%ymm4,%%ymm3,%%ymm3		\n\t"\
		"vmovaps	%%ymm2,     (%%rbx)	\n\t"\
		"vmovaps	%%ymm3,0x020(%%rbx)	\n\t"\
		"vaddpd	      %%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vsubpd	      %%ymm5,%%ymm2,%%ymm2		\n\t"\
		"vaddpd	      %%ymm4,%%ymm3,%%ymm3		\n\t"\
		"vmovaps	%%ymm2,     (%%rdx)	\n\t"\
		"vmovaps	%%ymm3,0x020(%%rdx)	\n\t"\
		:					/* outputs: none */\
		: [__add] "m" (Xadd)	/* All inputs from memory addresses here */\
		 ,[__p01] "m" (Xp01)\
		 ,[__p02] "m" (Xp02)\
		 ,[__p03] "m" (Xp03)\
		 ,[__p04] "m" (Xp04)\
		 ,[__p08] "m" (Xp08)\
		 ,[__p10] "m" (Xp10)\
		 ,[__p18] "m" (Xp18)\
		 ,[__r00] "m" (Xr00)\
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__cc0] "m" (Xcc0)\
		: "rax","rbx","rcx","rdx","rdi","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","memory"		/* Clobbered registers */\
	);\
	}


	#define SSE2_RADIX32_DIF_NOTWIDDLE(Xadd,Xp01,Xp02,Xp03,Xp04,Xp08,Xp10,Xp18,Xr00,Xisrt2,Xcc0)\
	{\
	__asm__ volatile (\
		"\n\t"/* SSE2_RADIX4_DIF_IN_PLACE(r00,r20,r10,r30): */\
		"movq	%[__isrt2],%%rdi	\n\t"\
		"movq	%[__r00],%%rax	\n\t"\
		"movq	$0x400,%%rbx	\n\t"\
		"movq	$0x200,%%rcx	\n\t"\
		"movq	$0x600,%%rdx	\n\t"\
		"addq	%%rax,%%rbx		\n\t"\
		"addq	%%rax,%%rcx		\n\t"\
		"addq	%%rax,%%rdx		\n\t"\
		"\n\t"\
		"vmovaps	     (%%rax),%%ymm0	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1	\n\t"\
		"vmovaps	     (%%rax),%%ymm2	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm3	\n\t"\
		"vaddpd	     (%%rbx),%%ymm0,%%ymm0		\n\t"\
		"vaddpd	0x020(%%rbx),%%ymm1,%%ymm1		\n\t"\
		"vsubpd	     (%%rbx),%%ymm2,%%ymm2		\n\t"\
		"vsubpd	0x020(%%rbx),%%ymm3,%%ymm3		\n\t"\
		"vmovaps	     (%%rcx),%%ymm4	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5	\n\t"\
		"vmovaps	     (%%rcx),%%ymm6	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm7	\n\t"\
		"vaddpd	     (%%rdx),%%ymm4,%%ymm4		\n\t"\
		"vaddpd	0x020(%%rdx),%%ymm5,%%ymm5		\n\t"\
		"vsubpd	     (%%rdx),%%ymm6,%%ymm6		\n\t"\
		"vsubpd	0x020(%%rdx),%%ymm7,%%ymm7		\n\t"\
		"\n\t"/* Finish radix-4 butterfly and store results into temporary-array slots: */\
		"vsubpd	      %%ymm4,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	      %%ymm5,%%ymm1,%%ymm1		\n\t"\
		"vmovaps	%%ymm0,     (%%rbx)	\n\t"\
		"vmovaps	%%ymm1,0x020(%%rbx)	\n\t"\
		"vaddpd	      %%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	      %%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm0,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	      %%ymm1,%%ymm5,%%ymm5		\n\t"\
		"vmovaps	%%ymm4,     (%%rax)	\n\t"\
		"vmovaps	%%ymm5,0x020(%%rax)	\n\t"\
		"vsubpd	      %%ymm7,%%ymm2,%%ymm2		\n\t"\
		"vsubpd	      %%ymm6,%%ymm3,%%ymm3		\n\t"\
		"vmovaps	%%ymm2,     (%%rcx)	\n\t"\
		"vmovaps	%%ymm3,0x020(%%rdx)	\n\t"\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm2,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm3,%%ymm6,%%ymm6		\n\t"\
		"vmovaps	%%ymm7,     (%%rdx)	\n\t"\
		"vmovaps	%%ymm6,0x020(%%rcx)	\n\t"\
		"\n\t"\
		"\n\t"/* SSE2_RADIX4_DIF_IN_PLACE_2NDOFTWO(r08,r28,r18,r38): */\
		"addq	$0x100,%%rax		\n\t"\
		"addq	$0x100,%%rbx		\n\t"\
		"addq	$0x100,%%rcx		\n\t"\
		"addq	$0x100,%%rdx		\n\t"\
		"vmovaps	      (%%rax),%%ymm0	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1	\n\t"\
		"vmovaps	      (%%rax),%%ymm2	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm3	\n\t"\
		"vaddpd	      (%%rbx),%%ymm0,%%ymm0	\n\t"\
		"vaddpd	0x020(%%rbx),%%ymm1,%%ymm1	\n\t"\
		"vsubpd	      (%%rbx),%%ymm2,%%ymm2	\n\t"\
		"vsubpd	0x020(%%rbx),%%ymm3,%%ymm3	\n\t"\
		"vmovaps	      (%%rcx),%%ymm4	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5	\n\t"\
		"vmovaps	      (%%rcx),%%ymm6	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm7	\n\t"\
		"vaddpd	      (%%rdx),%%ymm6,%%ymm6	\n\t"\
		"vaddpd	0x020(%%rdx),%%ymm7,%%ymm7	\n\t"\
		"vsubpd	      (%%rdx),%%ymm4,%%ymm4	\n\t"\
		"vsubpd	0x020(%%rdx),%%ymm5,%%ymm5	\n\t"\
		"\n\t"/* Finish radix-4 butterfly and store results into temporary-array slots: */\
		"vsubpd	      %%ymm6,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	      %%ymm7,%%ymm1,%%ymm1		\n\t"\
		"vsubpd	      %%ymm5,%%ymm2,%%ymm2		\n\t"\
		"vmovaps	%%ymm0,     (%%rbx)	\n\t"\
		"vmovaps	%%ymm1,0x020(%%rbx)	\n\t"\
		"vsubpd	      %%ymm4,%%ymm3,%%ymm3		\n\t"\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	      %%ymm0,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm2,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm1,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm3,%%ymm4,%%ymm4		\n\t"\
		"vmovaps	%%ymm6,     (%%rax)	\n\t"\
		"vmovaps	%%ymm7,0x020(%%rax)	\n\t"\
		"vmovaps	(%%rdi),%%ymm0	\n\t"/* isrt2 */\
		"vmovaps	%%ymm2,%%ymm6		\n\t"\
		"vmovaps	%%ymm5,%%ymm7		\n\t"\
		"vsubpd	      %%ymm4,%%ymm2,%%ymm2		\n\t"\
		"vsubpd	      %%ymm3,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm4,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm3,%%ymm7,%%ymm7		\n\t"\
		"vmulpd	      %%ymm0,%%ymm2,%%ymm2		\n\t"\
		"vmulpd	      %%ymm0,%%ymm5,%%ymm5		\n\t"\
		"vmulpd	      %%ymm0,%%ymm6,%%ymm6		\n\t"\
		"vmulpd	      %%ymm0,%%ymm7,%%ymm7		\n\t"\
		"vmovaps	%%ymm2,     (%%rcx)	\n\t"\
		"vmovaps	%%ymm5,     (%%rdx)	\n\t"\
		"vmovaps	%%ymm6,0x020(%%rcx)	\n\t"\
		"vmovaps	%%ymm7,0x020(%%rdx)	\n\t"\
		"\n\t"/* SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS(r00,r10,r20,r30,r08,r18,r28,r38): */\
		"subq	$0x100,%%rax	\n\t"/* r00 */\
		"subq	$0x400,%%rbx	\n\t"/* r08 */\
		"addq	$0x100,%%rcx	\n\t"/* r20 */\
		"subq	$0x200,%%rdx	\n\t"/* r28 */\
		"vmovaps	     (%%rax),%%ymm0	\n\t"\
		"vmovaps	     (%%rcx),%%ymm4	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5	\n\t"\
		"vmovaps	     (%%rbx),%%ymm2	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm7	\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm3	\n\t"\
		"vmovaps	     (%%rdx),%%ymm6	\n\t"\
		"vsubpd	      %%ymm2,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	      %%ymm7,%%ymm4,%%ymm4		\n\t"\
		"vsubpd	      %%ymm3,%%ymm1,%%ymm1		\n\t"\
		"vsubpd	      %%ymm6,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm2,%%ymm2,%%ymm2		\n\t"\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm3,%%ymm3,%%ymm3		\n\t"\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm0,%%ymm2,%%ymm2		\n\t"\
		"vaddpd	      %%ymm4,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm1,%%ymm3,%%ymm3		\n\t"\
		"vaddpd	      %%ymm5,%%ymm6,%%ymm6		\n\t"\
		"vmovaps	%%ymm0,     (%%rbx)	\n\t"\
		"vmovaps	%%ymm4,     (%%rcx)	\n\t"\
		"vmovaps	%%ymm1,0x020(%%rbx)	\n\t"\
		"vmovaps	%%ymm5,0x020(%%rdx)	\n\t"\
		"vmovaps	%%ymm2,     (%%rax)	\n\t"\
		"vmovaps	%%ymm7,     (%%rdx)	\n\t"\
		"vmovaps	%%ymm3,0x020(%%rax)	\n\t"\
		"vmovaps	%%ymm6,0x020(%%rcx)	\n\t"\
		"addq	$0x200,%%rax	\n\t"/* r10 */\
		"addq	$0x200,%%rbx	\n\t"/* r18 */\
		"addq	$0x200,%%rcx	\n\t"/* r30 */\
		"addq	$0x200,%%rdx	\n\t"/* r38 */\
		"vmovaps	     (%%rax),%%ymm0	\n\t"\
		"vmovaps	     (%%rcx),%%ymm4	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5	\n\t"\
		"vmovaps	     (%%rbx),%%ymm2	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm7	\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm3	\n\t"\
		"vmovaps	     (%%rdx),%%ymm6	\n\t"\
		"vsubpd	      %%ymm2,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	      %%ymm7,%%ymm4,%%ymm4		\n\t"\
		"vsubpd	      %%ymm3,%%ymm1,%%ymm1		\n\t"\
		"vsubpd	      %%ymm6,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm2,%%ymm2,%%ymm2		\n\t"\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm3,%%ymm3,%%ymm3		\n\t"\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm0,%%ymm2,%%ymm2		\n\t"\
		"vaddpd	      %%ymm4,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm1,%%ymm3,%%ymm3		\n\t"\
		"vaddpd	      %%ymm5,%%ymm6,%%ymm6		\n\t"\
		"\n\t"\
		"vmovaps	%%ymm0,     (%%rbx)	\n\t"\
		"vmovaps	%%ymm4,     (%%rcx)	\n\t"\
		"vmovaps	%%ymm1,0x020(%%rbx)	\n\t"\
		"vmovaps	%%ymm5,0x020(%%rdx)	\n\t"\
		"vmovaps	%%ymm2,     (%%rax)	\n\t"\
		"vmovaps	%%ymm7,     (%%rdx)	\n\t"\
		"vmovaps	%%ymm3,0x020(%%rax)	\n\t"\
		"vmovaps	%%ymm6,0x020(%%rcx)	\n\t"\
		"\n\t"\
		"\n\t"/* SSE2_RADIX4_DIF_IN_PLACE(r04,r24,r14,r34): */\
		"subq	$0x180,%%rax	\n\t"/* r04 */\
		"addq	$0x180,%%rbx	\n\t"/* r24 */\
		"subq	$0x380,%%rcx	\n\t"/* r14 */\
		"subq	$0x080,%%rdx	\n\t"/* r34 */\
		"\n\t"\
		"vmovaps	     (%%rax),%%ymm0	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1	\n\t"\
		"vmovaps	     (%%rax),%%ymm2	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm3	\n\t"\
		"vaddpd	     (%%rbx),%%ymm0,%%ymm0		\n\t"\
		"vaddpd	0x020(%%rbx),%%ymm1,%%ymm1		\n\t"\
		"vsubpd	     (%%rbx),%%ymm2,%%ymm2		\n\t"\
		"vsubpd	0x020(%%rbx),%%ymm3,%%ymm3		\n\t"\
		"vmovaps	     (%%rcx),%%ymm4	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5	\n\t"\
		"vmovaps	     (%%rcx),%%ymm6	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm7	\n\t"\
		"vaddpd	     (%%rdx),%%ymm4,%%ymm4		\n\t"\
		"vaddpd	0x020(%%rdx),%%ymm5,%%ymm5		\n\t"\
		"vsubpd	     (%%rdx),%%ymm6,%%ymm6		\n\t"\
		"vsubpd	0x020(%%rdx),%%ymm7,%%ymm7		\n\t"\
		"\n\t"/* Finish radix-4 butterfly and store results into temporary-array slots: */\
		"vsubpd	      %%ymm4,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	      %%ymm5,%%ymm1,%%ymm1		\n\t"\
		"vmovaps	%%ymm0,     (%%rbx)	\n\t"\
		"vmovaps	%%ymm1,0x020(%%rbx)	\n\t"\
		"vaddpd	      %%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	      %%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm0,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	      %%ymm1,%%ymm5,%%ymm5		\n\t"\
		"vmovaps	%%ymm4,     (%%rax)	\n\t"\
		"vmovaps	%%ymm5,0x020(%%rax)	\n\t"\
		"vsubpd	      %%ymm7,%%ymm2,%%ymm2		\n\t"\
		"vsubpd	      %%ymm6,%%ymm3,%%ymm3		\n\t"\
		"vmovaps	%%ymm2,     (%%rcx)	\n\t"\
		"vmovaps	%%ymm3,0x020(%%rdx)	\n\t"\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm2,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm3,%%ymm6,%%ymm6		\n\t"\
		"vmovaps	%%ymm7,     (%%rdx)	\n\t"\
		"vmovaps	%%ymm6,0x020(%%rcx)	\n\t"\
		"\n\t"\
		"\n\t"/* SSE2_RADIX4_DIF_IN_PLACE_2NDOFTWO(r0C,r2C,r1C,r3C): */\
		"addq	$0x100,%%rax		\n\t"\
		"addq	$0x100,%%rbx		\n\t"\
		"addq	$0x100,%%rcx		\n\t"\
		"addq	$0x100,%%rdx		\n\t"\
		"vmovaps	     (%%rax),%%ymm0	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1	\n\t"\
		"vmovaps	     (%%rax),%%ymm2	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm3	\n\t"\
		"vaddpd	     (%%rbx),%%ymm0,%%ymm0		\n\t"\
		"vaddpd	0x020(%%rbx),%%ymm1,%%ymm1		\n\t"\
		"vsubpd	     (%%rbx),%%ymm2,%%ymm2		\n\t"\
		"vsubpd	0x020(%%rbx),%%ymm3,%%ymm3		\n\t"\
		"vmovaps	     (%%rcx),%%ymm4	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5	\n\t"\
		"vmovaps	     (%%rcx),%%ymm6	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm7	\n\t"\
		"vaddpd	     (%%rdx),%%ymm6,%%ymm6		\n\t"\
		"vaddpd	0x020(%%rdx),%%ymm7,%%ymm7		\n\t"\
		"vsubpd	     (%%rdx),%%ymm4,%%ymm4		\n\t"\
		"vsubpd	0x020(%%rdx),%%ymm5,%%ymm5		\n\t"\
		"\n\t"/* Finish radix-4 butterfly and store results into temporary-array slots: */\
		"vsubpd	      %%ymm6,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	      %%ymm7,%%ymm1,%%ymm1		\n\t"\
		"vsubpd	      %%ymm5,%%ymm2,%%ymm2		\n\t"\
		"vmovaps	%%ymm0,     (%%rbx)	\n\t"\
		"vmovaps	%%ymm1,0x020(%%rbx)	\n\t"\
		"vsubpd	      %%ymm4,%%ymm3,%%ymm3		\n\t"\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	      %%ymm0,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm2,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm1,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm3,%%ymm4,%%ymm4		\n\t"\
		"vmovaps	%%ymm6,     (%%rax)	\n\t"\
		"vmovaps	%%ymm7,0x020(%%rax)	\n\t"\
		"vmovaps	(%%rdi),%%ymm0	\n\t"/* isrt2 */\
		"vmovaps	%%ymm2,%%ymm6		\n\t"\
		"vmovaps	%%ymm5,%%ymm7		\n\t"\
		"vsubpd	      %%ymm4,%%ymm2,%%ymm2		\n\t"\
		"vsubpd	      %%ymm3,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm4,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm3,%%ymm7,%%ymm7		\n\t"\
		"vmulpd	      %%ymm0,%%ymm2,%%ymm2		\n\t"\
		"vmulpd	      %%ymm0,%%ymm5,%%ymm5		\n\t"\
		"vmulpd	      %%ymm0,%%ymm6,%%ymm6		\n\t"\
		"vmulpd	      %%ymm0,%%ymm7,%%ymm7		\n\t"\
		"vmovaps	%%ymm2,     (%%rcx)	\n\t"\
		"vmovaps	%%ymm5,     (%%rdx)	\n\t"\
		"vmovaps	%%ymm6,0x020(%%rcx)	\n\t"\
		"vmovaps	%%ymm7,0x020(%%rdx)	\n\t"\
		"\n\t"/* SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS(r04,r14,r24,r34,r0C,r1C,r2C,r3C): */\
		"subq	$0x100,%%rax	\n\t"/* r04 */\
		"subq	$0x400,%%rbx	\n\t"/* r0C */\
		"addq	$0x100,%%rcx	\n\t"/* r24 */\
		"subq	$0x200,%%rdx	\n\t"/* r2C */\
		"vmovaps	     (%%rax),%%ymm0	\n\t"\
		"vmovaps	     (%%rcx),%%ymm4	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5	\n\t"\
		"vmovaps	     (%%rbx),%%ymm2	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm7	\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm3	\n\t"\
		"vmovaps	     (%%rdx),%%ymm6	\n\t"\
		"vsubpd	      %%ymm2,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	      %%ymm7,%%ymm4,%%ymm4		\n\t"\
		"vsubpd	      %%ymm3,%%ymm1,%%ymm1		\n\t"\
		"vsubpd	      %%ymm6,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm2,%%ymm2,%%ymm2		\n\t"\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm3,%%ymm3,%%ymm3		\n\t"\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm0,%%ymm2,%%ymm2		\n\t"\
		"vaddpd	      %%ymm4,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm1,%%ymm3,%%ymm3		\n\t"\
		"vaddpd	      %%ymm5,%%ymm6,%%ymm6		\n\t"\
		"vmovaps	%%ymm0,     (%%rbx)	\n\t"\
		"vmovaps	%%ymm4,     (%%rcx)	\n\t"\
		"vmovaps	%%ymm1,0x020(%%rbx)	\n\t"\
		"vmovaps	%%ymm5,0x020(%%rdx)	\n\t"\
		"vmovaps	%%ymm2,     (%%rax)	\n\t"\
		"vmovaps	%%ymm7,     (%%rdx)	\n\t"\
		"vmovaps	%%ymm3,0x020(%%rax)	\n\t"\
		"vmovaps	%%ymm6,0x020(%%rcx)	\n\t"\
		"addq	$0x200,%%rax	\n\t"/* r14 */\
		"addq	$0x200,%%rbx	\n\t"/* r1C */\
		"addq	$0x200,%%rcx	\n\t"/* r34 */\
		"addq	$0x200,%%rdx	\n\t"/* r3C */\
		"vmovaps	     (%%rax),%%ymm0	\n\t"\
		"vmovaps	     (%%rcx),%%ymm4	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5	\n\t"\
		"vmovaps	     (%%rbx),%%ymm2	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm7	\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm3	\n\t"\
		"vmovaps	     (%%rdx),%%ymm6	\n\t"\
		"vsubpd	      %%ymm2,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	      %%ymm7,%%ymm4,%%ymm4		\n\t"\
		"vsubpd	      %%ymm3,%%ymm1,%%ymm1		\n\t"\
		"vsubpd	      %%ymm6,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm2,%%ymm2,%%ymm2		\n\t"\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm3,%%ymm3,%%ymm3		\n\t"\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm0,%%ymm2,%%ymm2		\n\t"\
		"vaddpd	      %%ymm4,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm1,%%ymm3,%%ymm3		\n\t"\
		"vaddpd	      %%ymm5,%%ymm6,%%ymm6		\n\t"\
		"\n\t"\
		"vmovaps	%%ymm0,     (%%rbx)	\n\t"\
		"vmovaps	%%ymm4,     (%%rcx)	\n\t"\
		"vmovaps	%%ymm1,0x020(%%rbx)	\n\t"\
		"vmovaps	%%ymm5,0x020(%%rdx)	\n\t"\
		"vmovaps	%%ymm2,     (%%rax)	\n\t"\
		"vmovaps	%%ymm7,     (%%rdx)	\n\t"\
		"vmovaps	%%ymm3,0x020(%%rax)	\n\t"\
		"vmovaps	%%ymm6,0x020(%%rcx)	\n\t"\
		"\n\t"\
		"\n\t"/* SSE2_RADIX4_DIF_IN_PLACE(r02,r22,r12,r32): */\
		"\n\t"\
		"subq	$0x240,%%rax	\n\t"/* r02 */\
		"addq	$0x0c0,%%rbx	\n\t"/* r22 */\
		"subq	$0x440,%%rcx	\n\t"/* r12 */\
		"subq	$0x140,%%rdx	\n\t"/* r32 */\
		"\n\t"\
		"vmovaps	     (%%rax),%%ymm0	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1	\n\t"\
		"vmovaps	     (%%rax),%%ymm2	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm3	\n\t"\
		"vaddpd	     (%%rbx),%%ymm0,%%ymm0		\n\t"\
		"vaddpd	0x020(%%rbx),%%ymm1,%%ymm1		\n\t"\
		"vsubpd	     (%%rbx),%%ymm2,%%ymm2		\n\t"\
		"vsubpd	0x020(%%rbx),%%ymm3,%%ymm3		\n\t"\
		"vmovaps	     (%%rcx),%%ymm4	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5	\n\t"\
		"vmovaps	     (%%rcx),%%ymm6	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm7	\n\t"\
		"vaddpd	     (%%rdx),%%ymm4,%%ymm4		\n\t"\
		"vaddpd	0x020(%%rdx),%%ymm5,%%ymm5		\n\t"\
		"vsubpd	     (%%rdx),%%ymm6,%%ymm6		\n\t"\
		"vsubpd	0x020(%%rdx),%%ymm7,%%ymm7		\n\t"\
		"\n\t"/* Finish radix-4 butterfly and store results into temporary-array slots: */\
		"vsubpd	      %%ymm4,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	      %%ymm5,%%ymm1,%%ymm1		\n\t"\
		"vmovaps	%%ymm0,     (%%rbx)	\n\t"\
		"vmovaps	%%ymm1,0x020(%%rbx)	\n\t"\
		"vaddpd	      %%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	      %%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm0,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	      %%ymm1,%%ymm5,%%ymm5		\n\t"\
		"vmovaps	%%ymm4,     (%%rax)	\n\t"\
		"vmovaps	%%ymm5,0x020(%%rax)	\n\t"\
		"vsubpd	      %%ymm7,%%ymm2,%%ymm2		\n\t"\
		"vsubpd	      %%ymm6,%%ymm3,%%ymm3		\n\t"\
		"vmovaps	%%ymm2,     (%%rcx)	\n\t"\
		"vmovaps	%%ymm3,0x020(%%rdx)	\n\t"\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm2,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm3,%%ymm6,%%ymm6		\n\t"\
		"vmovaps	%%ymm7,     (%%rdx)	\n\t"\
		"vmovaps	%%ymm6,0x020(%%rcx)	\n\t"\
		"\n\t"\
		"\n\t"/* SSE2_RADIX4_DIF_IN_PLACE_2NDOFTWO(r0A,r2A,r1A,r3A): */\
		"addq	$0x100,%%rax		\n\t"\
		"addq	$0x100,%%rbx		\n\t"\
		"addq	$0x100,%%rcx		\n\t"\
		"addq	$0x100,%%rdx		\n\t"\
		"vmovaps	     (%%rax),%%ymm0	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1	\n\t"\
		"vmovaps	     (%%rax),%%ymm2	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm3	\n\t"\
		"vaddpd	     (%%rbx),%%ymm0,%%ymm0		\n\t"\
		"vaddpd	0x020(%%rbx),%%ymm1,%%ymm1		\n\t"\
		"vsubpd	     (%%rbx),%%ymm2,%%ymm2		\n\t"\
		"vsubpd	0x020(%%rbx),%%ymm3,%%ymm3		\n\t"\
		"vmovaps	     (%%rcx),%%ymm4	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5	\n\t"\
		"vmovaps	     (%%rcx),%%ymm6	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm7	\n\t"\
		"vaddpd	     (%%rdx),%%ymm6,%%ymm6		\n\t"\
		"vaddpd	0x020(%%rdx),%%ymm7,%%ymm7		\n\t"\
		"vsubpd	     (%%rdx),%%ymm4,%%ymm4		\n\t"\
		"vsubpd	0x020(%%rdx),%%ymm5,%%ymm5		\n\t"\
		"\n\t"/* Finish radix-4 butterfly and store results into temporary-array slots: */\
		"vsubpd	      %%ymm6,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	      %%ymm7,%%ymm1,%%ymm1		\n\t"\
		"vsubpd	      %%ymm5,%%ymm2,%%ymm2		\n\t"\
		"vmovaps	%%ymm0,     (%%rbx)	\n\t"\
		"vmovaps	%%ymm1,0x020(%%rbx)	\n\t"\
		"vsubpd	      %%ymm4,%%ymm3,%%ymm3		\n\t"\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	      %%ymm0,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm2,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm1,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm3,%%ymm4,%%ymm4		\n\t"\
		"vmovaps	%%ymm6,     (%%rax)	\n\t"\
		"vmovaps	%%ymm7,0x020(%%rax)	\n\t"\
		"vmovaps	(%%rdi),%%ymm0	\n\t"/* isrt2 */\
		"vmovaps	%%ymm2,%%ymm6		\n\t"\
		"vmovaps	%%ymm5,%%ymm7		\n\t"\
		"vsubpd	      %%ymm4,%%ymm2,%%ymm2		\n\t"\
		"vsubpd	      %%ymm3,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm4,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm3,%%ymm7,%%ymm7		\n\t"\
		"vmulpd	      %%ymm0,%%ymm2,%%ymm2		\n\t"\
		"vmulpd	      %%ymm0,%%ymm5,%%ymm5		\n\t"\
		"vmulpd	      %%ymm0,%%ymm6,%%ymm6		\n\t"\
		"vmulpd	      %%ymm0,%%ymm7,%%ymm7		\n\t"\
		"vmovaps	%%ymm2,     (%%rcx)	\n\t"\
		"vmovaps	%%ymm5,     (%%rdx)	\n\t"\
		"vmovaps	%%ymm6,0x020(%%rcx)	\n\t"\
		"vmovaps	%%ymm7,0x020(%%rdx)	\n\t"\
		"\n\t"/* SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS(r02,r12,r22,r32,r0A,r1A,r2A,r3A): */\
		"subq	$0x100,%%rax	\n\t"/* r02 */\
		"subq	$0x400,%%rbx	\n\t"/* r0A */\
		"addq	$0x100,%%rcx	\n\t"/* r22 */\
		"subq	$0x200,%%rdx	\n\t"/* r2A */\
		"vmovaps	     (%%rax),%%ymm0	\n\t"\
		"vmovaps	     (%%rcx),%%ymm4	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5	\n\t"\
		"vmovaps	     (%%rbx),%%ymm2	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm7	\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm3	\n\t"\
		"vmovaps	     (%%rdx),%%ymm6	\n\t"\
		"vsubpd	      %%ymm2,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	      %%ymm7,%%ymm4,%%ymm4		\n\t"\
		"vsubpd	      %%ymm3,%%ymm1,%%ymm1		\n\t"\
		"vsubpd	      %%ymm6,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm2,%%ymm2,%%ymm2		\n\t"\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm3,%%ymm3,%%ymm3		\n\t"\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm0,%%ymm2,%%ymm2		\n\t"\
		"vaddpd	      %%ymm4,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm1,%%ymm3,%%ymm3		\n\t"\
		"vaddpd	      %%ymm5,%%ymm6,%%ymm6		\n\t"\
		"vmovaps	%%ymm0,     (%%rbx)	\n\t"\
		"vmovaps	%%ymm4,     (%%rcx)	\n\t"\
		"vmovaps	%%ymm1,0x020(%%rbx)	\n\t"\
		"vmovaps	%%ymm5,0x020(%%rdx)	\n\t"\
		"vmovaps	%%ymm2,     (%%rax)	\n\t"\
		"vmovaps	%%ymm7,     (%%rdx)	\n\t"\
		"vmovaps	%%ymm3,0x020(%%rax)	\n\t"\
		"vmovaps	%%ymm6,0x020(%%rcx)	\n\t"\
		"addq	$0x200,%%rax	\n\t"/* r12 */\
		"addq	$0x200,%%rbx	\n\t"/* r1A */\
		"addq	$0x200,%%rcx	\n\t"/* r32 */\
		"addq	$0x200,%%rdx	\n\t"/* r3A */\
		"vmovaps	     (%%rax),%%ymm0	\n\t"\
		"vmovaps	     (%%rcx),%%ymm4	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5	\n\t"\
		"vmovaps	     (%%rbx),%%ymm2	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm7	\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm3	\n\t"\
		"vmovaps	     (%%rdx),%%ymm6	\n\t"\
		"vsubpd	      %%ymm2,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	      %%ymm7,%%ymm4,%%ymm4		\n\t"\
		"vsubpd	      %%ymm3,%%ymm1,%%ymm1		\n\t"\
		"vsubpd	      %%ymm6,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm2,%%ymm2,%%ymm2		\n\t"\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm3,%%ymm3,%%ymm3		\n\t"\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm0,%%ymm2,%%ymm2		\n\t"\
		"vaddpd	      %%ymm4,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm1,%%ymm3,%%ymm3		\n\t"\
		"vaddpd	      %%ymm5,%%ymm6,%%ymm6		\n\t"\
		"\n\t"\
		"vmovaps	%%ymm0,     (%%rbx)	\n\t"\
		"vmovaps	%%ymm4,     (%%rcx)	\n\t"\
		"vmovaps	%%ymm1,0x020(%%rbx)	\n\t"\
		"vmovaps	%%ymm5,0x020(%%rdx)	\n\t"\
		"vmovaps	%%ymm2,     (%%rax)	\n\t"\
		"vmovaps	%%ymm7,     (%%rdx)	\n\t"\
		"vmovaps	%%ymm3,0x020(%%rax)	\n\t"\
		"vmovaps	%%ymm6,0x020(%%rcx)	\n\t"\
		"\n\t"\
		"\n\t"/* SSE2_RADIX4_DIF_IN_PLACE(r06,r26,r16,r36): */\
		"\n\t"\
		"subq	$0x180,%%rax	\n\t"/* r06 */\
		"addq	$0x180,%%rbx	\n\t"/* r26 */\
		"subq	$0x380,%%rcx	\n\t"/* r16 */\
		"subq	$0x080,%%rdx	\n\t"/* r36 */\
		"\n\t"\
		"vmovaps	     (%%rax),%%ymm0	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1	\n\t"\
		"vmovaps	     (%%rax),%%ymm2	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm3	\n\t"\
		"vaddpd	     (%%rbx),%%ymm0,%%ymm0		\n\t"\
		"vaddpd	0x020(%%rbx),%%ymm1,%%ymm1		\n\t"\
		"vsubpd	     (%%rbx),%%ymm2,%%ymm2		\n\t"\
		"vsubpd	0x020(%%rbx),%%ymm3,%%ymm3		\n\t"\
		"vmovaps	     (%%rcx),%%ymm4	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5	\n\t"\
		"vmovaps	     (%%rcx),%%ymm6	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm7	\n\t"\
		"vaddpd	     (%%rdx),%%ymm4,%%ymm4		\n\t"\
		"vaddpd	0x020(%%rdx),%%ymm5,%%ymm5		\n\t"\
		"vsubpd	     (%%rdx),%%ymm6,%%ymm6		\n\t"\
		"vsubpd	0x020(%%rdx),%%ymm7,%%ymm7		\n\t"\
		"\n\t"/* Finish radix-4 butterfly and store results into temporary-array slots: */\
		"vsubpd	      %%ymm4,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	      %%ymm5,%%ymm1,%%ymm1		\n\t"\
		"vmovaps	%%ymm0,     (%%rbx)	\n\t"\
		"vmovaps	%%ymm1,0x020(%%rbx)	\n\t"\
		"vaddpd	      %%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	      %%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm0,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	      %%ymm1,%%ymm5,%%ymm5		\n\t"\
		"vmovaps	%%ymm4,     (%%rax)	\n\t"\
		"vmovaps	%%ymm5,0x020(%%rax)	\n\t"\
		"vsubpd	      %%ymm7,%%ymm2,%%ymm2		\n\t"\
		"vsubpd	      %%ymm6,%%ymm3,%%ymm3		\n\t"\
		"vmovaps	%%ymm2,     (%%rcx)	\n\t"\
		"vmovaps	%%ymm3,0x020(%%rdx)	\n\t"\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm2,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm3,%%ymm6,%%ymm6		\n\t"\
		"vmovaps	%%ymm7,     (%%rdx)	\n\t"\
		"vmovaps	%%ymm6,0x020(%%rcx)	\n\t"\
		"\n\t"\
		"\n\t"/* SSE2_RADIX4_DIF_IN_PLACE_2NDOFTWO(r0E,r2E,r1E,r3E): */\
		"addq	$0x100,%%rax		\n\t"\
		"addq	$0x100,%%rbx		\n\t"\
		"addq	$0x100,%%rcx		\n\t"\
		"addq	$0x100,%%rdx		\n\t"\
		"vmovaps	     (%%rax),%%ymm0	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1	\n\t"\
		"vmovaps	     (%%rax),%%ymm2	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm3	\n\t"\
		"vaddpd	     (%%rbx),%%ymm0,%%ymm0		\n\t"\
		"vaddpd	0x020(%%rbx),%%ymm1,%%ymm1		\n\t"\
		"vsubpd	     (%%rbx),%%ymm2,%%ymm2		\n\t"\
		"vsubpd	0x020(%%rbx),%%ymm3,%%ymm3		\n\t"\
		"vmovaps	     (%%rcx),%%ymm4	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5	\n\t"\
		"vmovaps	     (%%rcx),%%ymm6	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm7	\n\t"\
		"vaddpd	     (%%rdx),%%ymm6,%%ymm6		\n\t"\
		"vaddpd	0x020(%%rdx),%%ymm7,%%ymm7		\n\t"\
		"vsubpd	     (%%rdx),%%ymm4,%%ymm4		\n\t"\
		"vsubpd	0x020(%%rdx),%%ymm5,%%ymm5		\n\t"\
		"\n\t"/* Finish radix-4 butterfly and store results into temporary-array slots: */\
		"vsubpd	      %%ymm6,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	      %%ymm7,%%ymm1,%%ymm1		\n\t"\
		"vsubpd	      %%ymm5,%%ymm2,%%ymm2		\n\t"\
		"vmovaps	%%ymm0,     (%%rbx)	\n\t"\
		"vmovaps	%%ymm1,0x020(%%rbx)	\n\t"\
		"vsubpd	      %%ymm4,%%ymm3,%%ymm3		\n\t"\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	      %%ymm0,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm2,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm1,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm3,%%ymm4,%%ymm4		\n\t"\
		"vmovaps	%%ymm6,     (%%rax)	\n\t"\
		"vmovaps	%%ymm7,0x020(%%rax)	\n\t"\
		"vmovaps	(%%rdi),%%ymm0	\n\t"/* isrt2 */\
		"vmovaps	%%ymm2,%%ymm6		\n\t"\
		"vmovaps	%%ymm5,%%ymm7		\n\t"\
		"vsubpd	      %%ymm4,%%ymm2,%%ymm2		\n\t"\
		"vsubpd	      %%ymm3,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm4,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm3,%%ymm7,%%ymm7		\n\t"\
		"vmulpd	      %%ymm0,%%ymm2,%%ymm2		\n\t"\
		"vmulpd	      %%ymm0,%%ymm5,%%ymm5		\n\t"\
		"vmulpd	      %%ymm0,%%ymm6,%%ymm6		\n\t"\
		"vmulpd	      %%ymm0,%%ymm7,%%ymm7		\n\t"\
		"vmovaps	%%ymm2,     (%%rcx)	\n\t"\
		"vmovaps	%%ymm5,     (%%rdx)	\n\t"\
		"vmovaps	%%ymm6,0x020(%%rcx)	\n\t"\
		"vmovaps	%%ymm7,0x020(%%rdx)	\n\t"\
		"\n\t"/* SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS(r06,r16,r26,r36,r0E,r1E,r2E,r3E): */\
		"subq	$0x100,%%rax	\n\t"/* r02 */\
		"subq	$0x400,%%rbx	\n\t"/* r0A */\
		"addq	$0x100,%%rcx	\n\t"/* r22 */\
		"subq	$0x200,%%rdx	\n\t"/* r2A */\
		"vmovaps	     (%%rax),%%ymm0	\n\t"\
		"vmovaps	     (%%rcx),%%ymm4	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5	\n\t"\
		"vmovaps	     (%%rbx),%%ymm2	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm7	\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm3	\n\t"\
		"vmovaps	     (%%rdx),%%ymm6	\n\t"\
		"vsubpd	      %%ymm2,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	      %%ymm7,%%ymm4,%%ymm4		\n\t"\
		"vsubpd	      %%ymm3,%%ymm1,%%ymm1		\n\t"\
		"vsubpd	      %%ymm6,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm2,%%ymm2,%%ymm2		\n\t"\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm3,%%ymm3,%%ymm3		\n\t"\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm0,%%ymm2,%%ymm2		\n\t"\
		"vaddpd	      %%ymm4,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm1,%%ymm3,%%ymm3		\n\t"\
		"vaddpd	      %%ymm5,%%ymm6,%%ymm6		\n\t"\
		"vmovaps	%%ymm0,     (%%rbx)	\n\t"\
		"vmovaps	%%ymm4,     (%%rcx)	\n\t"\
		"vmovaps	%%ymm1,0x020(%%rbx)	\n\t"\
		"vmovaps	%%ymm5,0x020(%%rdx)	\n\t"\
		"vmovaps	%%ymm2,     (%%rax)	\n\t"\
		"vmovaps	%%ymm7,     (%%rdx)	\n\t"\
		"vmovaps	%%ymm3,0x020(%%rax)	\n\t"\
		"vmovaps	%%ymm6,0x020(%%rcx)	\n\t"\
		"addq	$0x200,%%rax	\n\t"/* r12 */\
		"addq	$0x200,%%rbx	\n\t"/* r1A */\
		"addq	$0x200,%%rcx	\n\t"/* r32 */\
		"addq	$0x200,%%rdx	\n\t"/* r3A */\
		"vmovaps	     (%%rax),%%ymm0	\n\t"\
		"vmovaps	     (%%rcx),%%ymm4	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5	\n\t"\
		"vmovaps	     (%%rbx),%%ymm2	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm7	\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm3	\n\t"\
		"vmovaps	     (%%rdx),%%ymm6	\n\t"\
		"vsubpd	      %%ymm2,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	      %%ymm7,%%ymm4,%%ymm4		\n\t"\
		"vsubpd	      %%ymm3,%%ymm1,%%ymm1		\n\t"\
		"vsubpd	      %%ymm6,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm2,%%ymm2,%%ymm2		\n\t"\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm3,%%ymm3,%%ymm3		\n\t"\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm0,%%ymm2,%%ymm2		\n\t"\
		"vaddpd	      %%ymm4,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm1,%%ymm3,%%ymm3		\n\t"\
		"vaddpd	      %%ymm5,%%ymm6,%%ymm6		\n\t"\
		"\n\t"\
		"vmovaps	%%ymm0,     (%%rbx)	\n\t"\
		"vmovaps	%%ymm4,     (%%rcx)	\n\t"\
		"vmovaps	%%ymm1,0x020(%%rbx)	\n\t"\
		"vmovaps	%%ymm5,0x020(%%rdx)	\n\t"\
		"vmovaps	%%ymm2,     (%%rax)	\n\t"\
		"vmovaps	%%ymm7,     (%%rdx)	\n\t"\
		"vmovaps	%%ymm3,0x020(%%rax)	\n\t"\
		"vmovaps	%%ymm6,0x020(%%rcx)	\n\t"\
		"\n\t"\
		"\n\t"/*...and now do eight radix-4 transforms, including the internal twiddle factors: */\
		"\n\t"/*...Block 1: t00,t10,t20,t30 in r00,04,02,06 - note swapped middle 2 indices! */\
		"\n\t"\
		"movq	%[__r00],%%rsi	\n\t"\
		"movq	%[__add],%%rax	\n\t"\
		"movslq	%[__p01],%%rbx	\n\t"\
		"movslq	%[__p02],%%rcx	\n\t"\
		"movslq	%[__p03],%%rdx	\n\t"\
		"shlq	$3,%%rbx		\n\t"\
		"shlq	$3,%%rcx		\n\t"\
		"shlq	$3,%%rdx		\n\t"\
		"addq	%%rax,%%rbx		\n\t"\
		"addq	%%rax,%%rcx		\n\t"\
		"addq	%%rax,%%rdx		\n\t"\
		"vmovaps	     (%%rsi),%%ymm0	\n\t"/* t00 */\
		"vmovaps	0x040(%%rsi),%%ymm4	\n\t"/* t20 */\
		"vmovaps	0x020(%%rsi),%%ymm1	\n\t"/* t01 */\
		"vmovaps	0x060(%%rsi),%%ymm5	\n\t"/* t21 */\
		"vmovaps	0x080(%%rsi),%%ymm2	\n\t"/* t10 */\
		"vmovaps	0x0c0(%%rsi),%%ymm6	\n\t"/* t30 */\
		"vmovaps	0x0a0(%%rsi),%%ymm3	\n\t"/* t11 */\
		"vmovaps	0x0e0(%%rsi),%%ymm7	\n\t"/* t31 */\
		"vsubpd	0x080(%%rsi),%%ymm0,%%ymm0	\n\t"/* t10=t00-rt */\
		"vsubpd	0x0c0(%%rsi),%%ymm4,%%ymm4	\n\t"/* t30=t20-rt */\
		"vsubpd	0x0a0(%%rsi),%%ymm1,%%ymm1	\n\t"/* t11=t01-it */\
		"vsubpd	0x0e0(%%rsi),%%ymm5,%%ymm5	\n\t"/* t31=t21-it */\
		"vaddpd	     (%%rsi),%%ymm2,%%ymm2	\n\t"/* t00=t00+rt */\
		"vaddpd	0x040(%%rsi),%%ymm6,%%ymm6	\n\t"/* t20=t20+rt */\
		"vaddpd	0x020(%%rsi),%%ymm3,%%ymm3	\n\t"/* t01=t01+it */\
		"vaddpd	0x060(%%rsi),%%ymm7,%%ymm7	\n\t"/* t21=t21+it */\
		"vsubpd	      %%ymm6,%%ymm2,%%ymm2		\n\t"/* t00 <- t00-t20 */\
		"vsubpd	      %%ymm5,%%ymm0,%%ymm0		\n\t"/* t10 <- t10-t31 */\
		"vsubpd	      %%ymm7,%%ymm3,%%ymm3		\n\t"/* t01 <- t01-t21 */\
		"vsubpd	      %%ymm4,%%ymm1,%%ymm1		\n\t"/* t11 <- t11-t30 */\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6		\n\t"/*          2*t20 */\
		"vaddpd	      %%ymm5,%%ymm5,%%ymm5		\n\t"/*          2*t31 */\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7		\n\t"/*          2*t21 */\
		"vaddpd	      %%ymm4,%%ymm4,%%ymm4		\n\t"/*          2*t30 */\
		"vmovaps	%%ymm2,     (%%rbx)	\n\t"/* a(jt+p1 ) */\
		"vmovaps	%%ymm0,     (%%rcx)	\n\t"/* a(jt+p2 ) */\
		"vmovaps	%%ymm3,0x020(%%rbx)	\n\t"/* a(jp+p1 ) */\
		"vmovaps	%%ymm1,0x020(%%rdx)	\n\t"/* a(jp+p3 ) */\
		"vaddpd	      %%ymm2,%%ymm6,%%ymm6		\n\t"/* t20 <- t00+t20 */\
		"vaddpd	      %%ymm0,%%ymm5,%%ymm5		\n\t"/* t31 <- t10+t31 */\
		"vaddpd	      %%ymm3,%%ymm7,%%ymm7		\n\t"/* t21 <- t01+t21 */\
		"vaddpd	      %%ymm1,%%ymm4,%%ymm4		\n\t"/* t30 <- t11+t30 */\
		"vmovaps	%%ymm6,     (%%rax)	\n\t"/* a(jt+p0 ) */\
		"vmovaps	%%ymm5,     (%%rdx)	\n\t"/* a(jt+p3 ) */\
		"vmovaps	%%ymm7,0x020(%%rax)	\n\t"/* a(jp+p0 ) */\
		"vmovaps	%%ymm4,0x020(%%rcx)	\n\t"/* a(jp+p2 ) */\
		"\n\t"/*...Block 5: t08,t18,t28,t38	*/\
		"addq	$0x100,%%rsi			\n\t"/* r08 */\
		"movslq	%[__p04],%%rdi		\n\t"\
		"shlq	$3,%%rdi			\n\t"\
		"addq	%%rdi,%%rax			\n\t"/* ap04 */\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	%%rdi,%%rcx			\n\t"\
		"addq	%%rdi,%%rdx			\n\t"\
		"movq	%[__isrt2],%%rdi	\n\t"\
		"\n\t"\
		"vmovaps	0x040(%%rsi),%%ymm4	\n\t"/* t28 */\
		"vmovaps	0x060(%%rsi),%%ymm5	\n\t"/* t29 */\
		"vmovaps	0x0c0(%%rsi),%%ymm6	\n\t"/* t38 */\
		"vmovaps	0x0e0(%%rsi),%%ymm7	\n\t"/* t39 */\
		"vmovaps	     (%%rsi),%%ymm0	\n\t"/* t08 */\
		"vsubpd	0x060(%%rsi),%%ymm4,%%ymm4	\n\t"/* t28-t29 */\
		"vmovaps	0x020(%%rsi),%%ymm1	\n\t"/* t09 */\
		"vaddpd	0x040(%%rsi),%%ymm5,%%ymm5	\n\t"/* t29+t28 */\
		"vmovaps	0x080(%%rsi),%%ymm2	\n\t"/* t18 */\
		"vmulpd	     (%%rdi),%%ymm4,%%ymm4		\n\t"/* t28 = (t28-t29)*ISRT2 */\
		"vmovaps	0x0a0(%%rsi),%%ymm3	\n\t"/* t19 */\
		"vmulpd	     (%%rdi),%%ymm5,%%ymm5		\n\t"/* t29 = (t29+t28)*ISRT2 */\
		"vsubpd	0x0a0(%%rsi),%%ymm0,%%ymm0	\n\t"/* t08=t08-t19*/\
		"vaddpd	0x0e0(%%rsi),%%ymm6,%%ymm6	\n\t"/* t38+t39 */\
		"vsubpd	0x080(%%rsi),%%ymm1,%%ymm1	\n\t"/* t19=t09-t18*/\
		"vsubpd	0x0c0(%%rsi),%%ymm7,%%ymm7	\n\t"/* t39-t38 */\
		"vaddpd	0x020(%%rsi),%%ymm2,%%ymm2	\n\t"/* t09=t18+t09*/\
		"vmulpd	     (%%rdi),%%ymm6,%%ymm6		\n\t"/*  rt = (t38+t39)*ISRT2 */\
		"vaddpd	     (%%rsi),%%ymm3,%%ymm3	\n\t"/* t18=t19+t08*/\
		"vmulpd	     (%%rdi),%%ymm7,%%ymm7		\n\t"/*  it = (t39-t38)*ISRT2 */\
		"vsubpd	      %%ymm6,%%ymm4,%%ymm4		\n\t"/* t28=t28-rt */\
		"vsubpd	      %%ymm7,%%ymm5,%%ymm5		\n\t"/* t29=t29-it */\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6		\n\t"/*      2* rt */\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7		\n\t"/*      2* it */\
		"vaddpd	      %%ymm4,%%ymm6,%%ymm6		\n\t"/* t38=t28+rt */\
		"vaddpd	      %%ymm5,%%ymm7,%%ymm7		\n\t"/* t39=t29+it */\
		"vsubpd	      %%ymm4,%%ymm0,%%ymm0		\n\t"/* t08-t28 */\
		"vsubpd	      %%ymm5,%%ymm2,%%ymm2		\n\t"/* t09-t29 */\
		"vaddpd	      %%ymm4,%%ymm4,%%ymm4		\n\t"/*   2*t28 */\
		"vaddpd	      %%ymm5,%%ymm5,%%ymm5		\n\t"/*   2*t29 */\
		"vsubpd	      %%ymm7,%%ymm3,%%ymm3		\n\t"/* t18-t39 */\
		"vsubpd	      %%ymm6,%%ymm1,%%ymm1		\n\t"/* t19-t38 */\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7		\n\t"/*   2*t39 */\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6		\n\t"/*   2*t38 */\
		"vmovaps	%%ymm0,     (%%rbx)	\n\t"/* a(jt+p1 ) */\
		"vmovaps	%%ymm2,0x020(%%rbx)	\n\t"/* a(jp+p1 ) */\
		"vmovaps	%%ymm3,     (%%rcx)	\n\t"/* a(jt+p2 ) */\
		"vmovaps	%%ymm1,0x020(%%rdx)	\n\t"/* a(jp+p3 ) */\
		"vaddpd	      %%ymm0,%%ymm4,%%ymm4		\n\t"/* t08+t28 */\
		"vaddpd	      %%ymm2,%%ymm5,%%ymm5		\n\t"/* t09+t29 */\
		"vaddpd	      %%ymm3,%%ymm7,%%ymm7		\n\t"/* t18+t39 */\
		"vaddpd	      %%ymm1,%%ymm6,%%ymm6		\n\t"/* t19+t38 */\
		"vmovaps	%%ymm4,     (%%rax)	\n\t"/* a(jt+p0 ) */\
		"vmovaps	%%ymm5,0x020(%%rax)	\n\t"/* a(jp+p0 ) */\
		"vmovaps	%%ymm7,     (%%rdx)	\n\t"/* a(jt+p3 ) */\
		"vmovaps	%%ymm6,0x020(%%rcx)	\n\t"/* a(jp+p2 ) */\
		"\n\t"/*...Block 3: t04,t14,t24,t34	*/\
		"addq	$0x300,%%rsi		\n\t"/* r20 */\
		"movslq	%[__p04],%%rdi		\n\t"\
		"shlq	$3,%%rdi			\n\t"\
		"addq	%%rdi,%%rax			\n\t"/* ap08 */\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	%%rdi,%%rcx			\n\t"\
		"addq	%%rdi,%%rdx			\n\t"\
		"movq	%[__cc0],%%rdi	\n\t"\
		"\n\t"\
		"vmovaps	0x040(%%rsi),%%ymm4	\n\t"/* t24 */\
		"vmovaps	0x0c0(%%rsi),%%ymm6	\n\t"/* t34 */\
		"vmovaps	0x060(%%rsi),%%ymm5	\n\t"/* t25 */\
		"vmovaps	0x0e0(%%rsi),%%ymm7	\n\t"/* t35 */\
		"vmovaps	0x040(%%rsi),%%ymm0	\n\t"/* copy t24 */\
		"vmovaps	0x0c0(%%rsi),%%ymm2	\n\t"/* copy t34 */\
		"vmovaps	0x060(%%rsi),%%ymm1	\n\t"/* copy t25 */\
		"vmovaps	0x0e0(%%rsi),%%ymm3	\n\t"/* copy t35 */\
		"vmulpd	     (%%rdi),%%ymm4,%%ymm4	\n\t"/* t24*c */\
		"vmulpd	0x020(%%rdi),%%ymm6,%%ymm6	\n\t"/* t34*s */\
		"vmulpd	0x020(%%rdi),%%ymm1,%%ymm1	\n\t"/* t25*s */\
		"vmulpd	     (%%rdi),%%ymm3,%%ymm3	\n\t"/* t35*c */\
		"vmulpd	     (%%rdi),%%ymm5,%%ymm5	\n\t"/* t25*c */\
		"vmulpd	0x020(%%rdi),%%ymm7,%%ymm7	\n\t"/* t35*s */\
		"vmulpd	0x020(%%rdi),%%ymm0,%%ymm0	\n\t"/* t24*s */\
		"vmulpd	     (%%rdi),%%ymm2,%%ymm2	\n\t"/* t34*c */\
		"vsubpd	      %%ymm1,%%ymm4,%%ymm4		\n\t"/* ~t24 */\
		"vsubpd	      %%ymm3,%%ymm6,%%ymm6		\n\t"/* rt */\
		"vaddpd	      %%ymm0,%%ymm5,%%ymm5		\n\t"/* ~t25 */\
		"vaddpd	      %%ymm2,%%ymm7,%%ymm7		\n\t"/* it */\
		"subq	$0x020,%%rdi			\n\t"/* isrt2 */\
		"vmovaps	0x080(%%rsi),%%ymm2	\n\t"/* t14 */\
		"vsubpd	      %%ymm6,%%ymm4,%%ymm4		\n\t"/*~t34=t24-rt */\
		"vmovaps	0x0a0(%%rsi),%%ymm3	\n\t"/* t15 */\
		"vsubpd	      %%ymm7,%%ymm5,%%ymm5		\n\t"/*~t35=t25-it */\
		"vsubpd	0x0a0(%%rsi),%%ymm2,%%ymm2	\n\t"/* t14-t15 */\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6		\n\t"/*      2* rt */\
		"vaddpd	0x080(%%rsi),%%ymm3,%%ymm3	\n\t"/* t15+t14 */\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7		\n\t"/*      2* it */\
		"vmulpd	     (%%rdi),%%ymm2,%%ymm2		\n\t"/* rt = (t14-t15)*ISRT2 */\
		"vaddpd	      %%ymm4,%%ymm6,%%ymm6		\n\t"/*~t24=t24+rt */\
		"vmulpd	     (%%rdi),%%ymm3,%%ymm3		\n\t"/* it = (t15+t14)*ISRT2 */\
		"vaddpd	      %%ymm5,%%ymm7,%%ymm7		\n\t"/*~t25=t25+it */\
		"vmovaps	     (%%rsi),%%ymm0	\n\t"/* t04 */\
		"vmovaps	0x020(%%rsi),%%ymm1	\n\t"/* t05 */\
		"vsubpd	      %%ymm2,%%ymm0,%%ymm0		\n\t"/*~t14=t04-rt */\
		"vsubpd	      %%ymm3,%%ymm1,%%ymm1		\n\t"/*~t15=t05-it */\
		"vaddpd	     (%%rsi),%%ymm2,%%ymm2	\n\t"/*~t04=rt +t04*/\
		"vaddpd	0x020(%%rsi),%%ymm3,%%ymm3	\n\t"/*~t05=it +t05*/\
		"vsubpd	      %%ymm6,%%ymm2,%%ymm2		\n\t"/* t04-t24 */\
		"vsubpd	      %%ymm5,%%ymm0,%%ymm0		\n\t"/* t14-t35 */\
		"vsubpd	      %%ymm7,%%ymm3,%%ymm3		\n\t"/* t05-t25 */\
		"vsubpd	      %%ymm4,%%ymm1,%%ymm1		\n\t"/* t15-t34 */\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6		\n\t"/*   2*t24 */\
		"vaddpd	      %%ymm5,%%ymm5,%%ymm5		\n\t"/*          2*t35 */\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7		\n\t"/*   2*t25 */\
		"vaddpd	      %%ymm4,%%ymm4,%%ymm4		\n\t"/*          2*t34 */\
		"vmovaps	%%ymm2,     (%%rbx)	\n\t"/* a(jt+p1 ) */\
		"vmovaps	%%ymm0,     (%%rcx)	\n\t"/* a(jt+p2 ) */\
		"vmovaps	%%ymm3,0x020(%%rbx)	\n\t"/* a(jp+p1 ) */\
		"vmovaps	%%ymm1,0x020(%%rdx)	\n\t"/* a(jp+p3 ) */\
		"vaddpd	      %%ymm2,%%ymm6,%%ymm6		\n\t"/* t04+t24 */\
		"vaddpd	      %%ymm0,%%ymm5,%%ymm5		\n\t"/* t14+t35 */\
		"vaddpd	      %%ymm3,%%ymm7,%%ymm7		\n\t"/* t05+t25 */\
		"vaddpd	      %%ymm1,%%ymm4,%%ymm4		\n\t"/* t15+t34 */\
		"vmovaps	%%ymm6,     (%%rax)	\n\t"/* a(jt+p0 ) */\
		"vmovaps	%%ymm5,     (%%rdx)	\n\t"/* a(jt+p3 ) */\
		"vmovaps	%%ymm7,0x020(%%rax)	\n\t"/* a(jp+p0 ) */\
		"vmovaps	%%ymm4,0x020(%%rcx)	\n\t"/* a(jp+p2 ) */\
		"\n\t"/*...Block 7: t0C,t1C,t2C,t3C	*/\
		"addq	$0x100,%%rsi			\n\t"/* r28 */\
		"movslq	%[__p04],%%rdi		\n\t"\
		"shlq	$3,%%rdi			\n\t"\
		"addq	%%rdi,%%rax			\n\t"/* ap0C */\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	%%rdi,%%rcx			\n\t"\
		"addq	%%rdi,%%rdx			\n\t"\
		"movq	%[__cc0],%%rdi	\n\t"\
		"\n\t"\
		"vmovaps	0x040(%%rsi),%%ymm4	\n\t"/* t2C */\
		"vmovaps	0x0c0(%%rsi),%%ymm6	\n\t"/* t3C */\
		"vmovaps	0x060(%%rsi),%%ymm5	\n\t"/* t2D */\
		"vmovaps	0x0e0(%%rsi),%%ymm7	\n\t"/* t3D */\
		"vmovaps	0x040(%%rsi),%%ymm0	\n\t"/* copy t2C */\
		"vmovaps	0x0c0(%%rsi),%%ymm2	\n\t"/* copy t3C */\
		"vmovaps	0x060(%%rsi),%%ymm1	\n\t"/* copy t2D */\
		"vmovaps	0x0e0(%%rsi),%%ymm3	\n\t"/* copy t3D */\
		"vmulpd	0x020(%%rdi),%%ymm4,%%ymm4	\n\t"/* t2C*s */\
		"vmulpd	     (%%rdi),%%ymm6,%%ymm6	\n\t"/* t3C*c */\
		"vmulpd	     (%%rdi),%%ymm1,%%ymm1	\n\t"/* t2D*c */\
		"vmulpd	0x020(%%rdi),%%ymm3,%%ymm3	\n\t"/* t3D*s */\
		"vmulpd	0x020(%%rdi),%%ymm5,%%ymm5	\n\t"/* t2D*s */\
		"vmulpd	     (%%rdi),%%ymm7,%%ymm7	\n\t"/* t3D*c */\
		"vmulpd	     (%%rdi),%%ymm0,%%ymm0	\n\t"/* t2C*c */\
		"vmulpd	0x020(%%rdi),%%ymm2,%%ymm2	\n\t"/* t3C*s */\
		"vsubpd	      %%ymm1,%%ymm4,%%ymm4		\n\t"/* ~t24 */\
		"vsubpd	      %%ymm3,%%ymm6,%%ymm6		\n\t"/* rt */\
		"vaddpd	      %%ymm0,%%ymm5,%%ymm5		\n\t"/* ~t25 */\
		"vaddpd	      %%ymm2,%%ymm7,%%ymm7		\n\t"/* it */\
		"subq	$0x020,%%rdi			\n\t"/* isrt2 */\
		"vmovaps	0x080(%%rsi),%%ymm2	\n\t"/* t14 */\
		"vsubpd	      %%ymm6,%%ymm4,%%ymm4		\n\t"/*~t2C=t2C-rt */\
		"vmovaps	0x0a0(%%rsi),%%ymm3	\n\t"/* t1D */\
		"vsubpd	      %%ymm7,%%ymm5,%%ymm5		\n\t"/*~t2D=t2D-it */\
		"vaddpd	0x0a0(%%rsi),%%ymm2,%%ymm2	\n\t"/* t1C+t1D */\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6		\n\t"/*      2* rt */\
		"vsubpd	0x080(%%rsi),%%ymm3,%%ymm3	\n\t"/* t1D-t1C */\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7		\n\t"/*      2* it */\
		"vmulpd	     (%%rdi),%%ymm2,%%ymm2		\n\t"/* rt = (t1C+t1D)*ISRT2 */\
		"vaddpd	      %%ymm4,%%ymm6,%%ymm6		\n\t"/*~t3C=t2C+rt */\
		"vmulpd	     (%%rdi),%%ymm3,%%ymm3		\n\t"/* it = (t1D-t1C)*ISRT2 */\
		"vaddpd	      %%ymm5,%%ymm7,%%ymm7		\n\t"/*~t3D=t2D+it */\
		"vmovaps	     (%%rsi),%%ymm0	\n\t"/* t0C */\
		"vmovaps	0x020(%%rsi),%%ymm1	\n\t"/* t0D */\
		"vsubpd	      %%ymm2,%%ymm0,%%ymm0		\n\t"/*~t0C=t0C-rt */\
		"vsubpd	      %%ymm3,%%ymm1,%%ymm1		\n\t"/*~t0D=t0D-it */\
		"vaddpd	     (%%rsi),%%ymm2,%%ymm2	\n\t"/*~t1C=rt +t0C*/\
		"vaddpd	0x020(%%rsi),%%ymm3,%%ymm3	\n\t"/*~t1D=it +t0D*/\
		"vsubpd	      %%ymm4,%%ymm0,%%ymm0		\n\t"/* t0C-t2C */\
		"vsubpd	      %%ymm7,%%ymm2,%%ymm2		\n\t"/* t1C-t3D */\
		"vsubpd	      %%ymm5,%%ymm1,%%ymm1		\n\t"/* t0D-t2D */\
		"vsubpd	      %%ymm6,%%ymm3,%%ymm3		\n\t"/* t1D-t3C */\
		"vaddpd	      %%ymm4,%%ymm4,%%ymm4		\n\t"/*   2*t2C */\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7		\n\t"/*   2*t3D */\
		"vaddpd	      %%ymm5,%%ymm5,%%ymm5		\n\t"/*   2*t2D */\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6		\n\t"/*   2*t3C */\
		"vmovaps	%%ymm0,     (%%rbx)	\n\t"/* a(jt+p1 ) */\
		"vmovaps	%%ymm2,     (%%rcx)	\n\t"/* a(jt+p2 ) */\
		"vmovaps	%%ymm1,0x020(%%rbx)	\n\t"/* a(jp+p1 ) */\
		"vmovaps	%%ymm3,0x020(%%rdx)	\n\t"/* a(jp+p3 ) */\
		"vaddpd	      %%ymm0,%%ymm4,%%ymm4		\n\t"/* t0C+t2C */\
		"vaddpd	      %%ymm2,%%ymm7,%%ymm7		\n\t"/* t1C+t3D */\
		"vaddpd	      %%ymm1,%%ymm5,%%ymm5		\n\t"/* t0D+t2D */\
		"vaddpd	      %%ymm3,%%ymm6,%%ymm6		\n\t"/* t1D+t3C */\
		"vmovaps	%%ymm4,     (%%rax)	\n\t"/* a(jt+p0 ) */\
		"vmovaps	%%ymm7,     (%%rdx)	\n\t"/* a(jt+p3 ) */\
		"vmovaps	%%ymm5,0x020(%%rax)	\n\t"/* a(jp+p0 ) */\
		"vmovaps	%%ymm6,0x020(%%rcx)	\n\t"/* a(jp+p2 ) */\
		"\n\t"/*...Block 2: t02,t12,t22,t32	*/\
		"subq	$0x300,%%rsi		\n\t"/* r10 */\
		"movslq	%[__p04],%%rdi		\n\t"\
		"shlq	$3,%%rdi			\n\t"\
		"addq	%%rdi,%%rax			\n\t"/* ap10 */\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	%%rdi,%%rcx			\n\t"\
		"addq	%%rdi,%%rdx			\n\t"\
		"movq	%[__cc0],%%rdi	\n\t"\
		"\n\t"\
		"vmovaps	0x040(%%rsi),%%ymm4	\n\t"/* t22 */\
		"vmovaps	0x0c0(%%rsi),%%ymm6	\n\t"/* t32 */\
		"vmovaps	0x060(%%rsi),%%ymm5	\n\t"/* t23 */\
		"vmovaps	0x0e0(%%rsi),%%ymm7	\n\t"/* t33 */\
		"vmovaps	0x040(%%rsi),%%ymm0	\n\t"/* copy t22 */\
		"vmovaps	0x0c0(%%rsi),%%ymm2	\n\t"/* copy t32 */\
		"vmovaps	0x060(%%rsi),%%ymm1	\n\t"/* copy t23 */\
		"vmovaps	0x0e0(%%rsi),%%ymm3	\n\t"/* copy t33 */\
		"vmulpd	0x040(%%rdi),%%ymm4,%%ymm4	\n\t"/* t22*c32_1 */\
		"vmulpd	0x080(%%rdi),%%ymm6,%%ymm6	\n\t"/* t32*c32_3 */\
		"vmulpd	0x060(%%rdi),%%ymm1,%%ymm1	\n\t"/* t23*s32_1 */\
		"vmulpd	0x0a0(%%rdi),%%ymm3,%%ymm3	\n\t"/* t33*s32_3 */\
		"vmulpd	0x040(%%rdi),%%ymm5,%%ymm5	\n\t"/* t23*c32_1 */\
		"vmulpd	0x080(%%rdi),%%ymm7,%%ymm7	\n\t"/* t33*c32_3 */\
		"vmulpd	0x060(%%rdi),%%ymm0,%%ymm0	\n\t"/* t22*s32_1 */\
		"vmulpd	0x0a0(%%rdi),%%ymm2,%%ymm2	\n\t"/* t32*s32_3 */\
		"vsubpd	      %%ymm1,%%ymm4,%%ymm4		\n\t"/* ~t22 */\
		"vsubpd	      %%ymm3,%%ymm6,%%ymm6		\n\t"/* rt */\
		"vaddpd	      %%ymm0,%%ymm5,%%ymm5		\n\t"/* ~t23 */\
		"vaddpd	      %%ymm2,%%ymm7,%%ymm7		\n\t"/* it */\
		"vmovaps	0x080(%%rsi),%%ymm2	\n\t"/* t12 */\
		"vmovaps	0x0a0(%%rsi),%%ymm0	\n\t"/* t13 */\
		"vmovaps	0x080(%%rsi),%%ymm1	\n\t"/* copy t12 */\
		"vmovaps	0x0a0(%%rsi),%%ymm3	\n\t"/* copy t13 */\
		"vsubpd	      %%ymm6,%%ymm4,%%ymm4		\n\t"/*~t32=t22-rt */\
		"vmulpd	     (%%rdi),%%ymm2,%%ymm2	\n\t"/* t12*c */\
		"vsubpd	      %%ymm7,%%ymm5,%%ymm5		\n\t"/*~t33=t23-it */\
		"vmulpd	0x020(%%rdi),%%ymm0,%%ymm0	\n\t"/* t13*s */\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6		\n\t"/*      2* rt */\
		"vmulpd	     (%%rdi),%%ymm3,%%ymm3	\n\t"/* t13*c */\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7		\n\t"/*      2* it */\
		"vmulpd	0x020(%%rdi),%%ymm1,%%ymm1	\n\t"/* t12*s */\
		"vaddpd	      %%ymm4,%%ymm6,%%ymm6		\n\t"/*~t22=t22+rt */\
		"vsubpd	      %%ymm0,%%ymm2,%%ymm2		\n\t"/* rt */\
		"vaddpd	      %%ymm5,%%ymm7,%%ymm7		\n\t"/*~t23=t23+it */\
		"vaddpd	      %%ymm1,%%ymm3,%%ymm3		\n\t"/* it */\
		"vmovaps	     (%%rsi),%%ymm0	\n\t"/* t02 */\
		"vmovaps	0x020(%%rsi),%%ymm1	\n\t"/* t03 */\
		"vsubpd	      %%ymm2,%%ymm0,%%ymm0		\n\t"/*~t12=t02-rt */\
		"vsubpd	      %%ymm3,%%ymm1,%%ymm1		\n\t"/*~t13=t03-it */\
		"vaddpd	     (%%rsi),%%ymm2,%%ymm2	\n\t"/*~t02=rt+t02 */\
		"vaddpd	0x020(%%rsi),%%ymm3,%%ymm3	\n\t"/*~t03=it+t03 */\
		"vsubpd	      %%ymm6,%%ymm2,%%ymm2		\n\t"/* t02-t22 */\
		"vsubpd	      %%ymm5,%%ymm0,%%ymm0		\n\t"/* t12-t33 */\
		"vsubpd	      %%ymm7,%%ymm3,%%ymm3		\n\t"/* t03-t23 */\
		"vsubpd	      %%ymm4,%%ymm1,%%ymm1		\n\t"/* t13-t32 */\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6		\n\t"/*   2*t22 */\
		"vaddpd	      %%ymm5,%%ymm5,%%ymm5		\n\t"/*   2*t33 */\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7		\n\t"/*   2*t23 */\
		"vaddpd	      %%ymm4,%%ymm4,%%ymm4		\n\t"/*   2*t32 */\
		"vmovaps	%%ymm2,     (%%rbx)	\n\t"/* a(jt+p1 ) */\
		"vmovaps	%%ymm0,     (%%rcx)	\n\t"/* a(jt+p2 ) */\
		"vmovaps	%%ymm3,0x020(%%rbx)	\n\t"/* a(jp+p1 ) */\
		"vmovaps	%%ymm1,0x020(%%rdx)	\n\t"/* a(jp+p3 ) */\
		"vaddpd	      %%ymm2,%%ymm6,%%ymm6		\n\t"/* t02+t22 */\
		"vaddpd	      %%ymm0,%%ymm5,%%ymm5		\n\t"/* t12+t33 */\
		"vaddpd	      %%ymm3,%%ymm7,%%ymm7		\n\t"/* t03+t23 */\
		"vaddpd	      %%ymm1,%%ymm4,%%ymm4		\n\t"/* t13+t32 */\
		"vmovaps	%%ymm6,     (%%rax)	\n\t"/* a(jt+p0 ) */\
		"vmovaps	%%ymm5,     (%%rdx)	\n\t"/* a(jt+p3 ) */\
		"vmovaps	%%ymm7,0x020(%%rax)	\n\t"/* a(jp+p0 ) */\
		"vmovaps	%%ymm4,0x020(%%rcx)	\n\t"/* a(jp+p2 ) */\
		"\n\t"/*...Block 6: t0A,t1A,t2A,t3A	*/\
		"addq	$0x100,%%rsi			\n\t"/* r18 */\
		"movslq	%[__p04],%%rdi		\n\t"\
		"shlq	$3,%%rdi			\n\t"\
		"addq	%%rdi,%%rax			\n\t"/* ap14 */\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	%%rdi,%%rcx			\n\t"\
		"addq	%%rdi,%%rdx			\n\t"\
		"movq	%[__cc0],%%rdi	\n\t"\
		"\n\t"\
		"vmovaps	0x040(%%rsi),%%ymm4	\n\t"/* t2A */\
		"vmovaps	0x0c0(%%rsi),%%ymm6	\n\t"/* t3A */\
		"vmovaps	0x060(%%rsi),%%ymm5	\n\t"/* t2B */\
		"vmovaps	0x0e0(%%rsi),%%ymm7	\n\t"/* t3B */\
		"vmovaps	0x040(%%rsi),%%ymm0	\n\t"/* copy t2A */\
		"vmovaps	0x0c0(%%rsi),%%ymm2	\n\t"/* copy t3A */\
		"vmovaps	0x060(%%rsi),%%ymm1	\n\t"/* copy t2B */\
		"vmovaps	0x0e0(%%rsi),%%ymm3	\n\t"/* copy t3B */\
		"vmulpd	0x0a0(%%rdi),%%ymm4,%%ymm4	\n\t"/* t2A*s32_3 */\
		"vmulpd	0x040(%%rdi),%%ymm6,%%ymm6	\n\t"/* t3A*c32_1 */\
		"vmulpd	0x080(%%rdi),%%ymm1,%%ymm1	\n\t"/* t2B*c32_3 */\
		"vmulpd	0x060(%%rdi),%%ymm3,%%ymm3	\n\t"/* t3B*s32_1 */\
		"vmulpd	0x0a0(%%rdi),%%ymm5,%%ymm5	\n\t"/* t2B*s32_3 */\
		"vmulpd	0x040(%%rdi),%%ymm7,%%ymm7	\n\t"/* t3B*c32_1 */\
		"vmulpd	0x080(%%rdi),%%ymm0,%%ymm0	\n\t"/* t2A*c32_3 */\
		"vmulpd	0x060(%%rdi),%%ymm2,%%ymm2	\n\t"/* t3A*s32_1 */\
		"vsubpd	      %%ymm1,%%ymm4,%%ymm4		\n\t"/* ~t2A */\
		"vaddpd	      %%ymm3,%%ymm6,%%ymm6		\n\t"/* rt */\
		"vaddpd	      %%ymm0,%%ymm5,%%ymm5		\n\t"/* ~t2B */\
		"vsubpd	      %%ymm2,%%ymm7,%%ymm7		\n\t"/* it */\
		"vmovaps	0x080(%%rsi),%%ymm2	\n\t"/* t1A */\
		"vmovaps	0x0a0(%%rsi),%%ymm0	\n\t"/* t1B */\
		"vmovaps	0x080(%%rsi),%%ymm1	\n\t"/* copy t1A */\
		"vmovaps	0x0a0(%%rsi),%%ymm3	\n\t"/* copy t1B */\
		"vsubpd	      %%ymm6,%%ymm4,%%ymm4		\n\t"/*~t2A=t2A-rt */\
		"vmulpd	0x020(%%rdi),%%ymm2,%%ymm2	\n\t"/* t1A*s */\
		"vsubpd	      %%ymm7,%%ymm5,%%ymm5		\n\t"/*~t2B=t2B-it */\
		"vmulpd	     (%%rdi),%%ymm0,%%ymm0	\n\t"/* t1B*c */\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6		\n\t"/*      2* rt */\
		"vmulpd	0x020(%%rdi),%%ymm3,%%ymm3	\n\t"/* t1B*s */\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7		\n\t"/*      2* it */\
		"vmulpd	     (%%rdi),%%ymm1,%%ymm1	\n\t"/* t1A*c */\
		"vaddpd	      %%ymm4,%%ymm6,%%ymm6		\n\t"/*~t3A=t2A+rt */\
		"vaddpd	      %%ymm0,%%ymm2,%%ymm2		\n\t"/* rt */\
		"vaddpd	      %%ymm5,%%ymm7,%%ymm7		\n\t"/*~t3B=t2B+it */\
		"vsubpd	      %%ymm1,%%ymm3,%%ymm3		\n\t"/* it */\
		"vmovaps	     (%%rsi),%%ymm0	\n\t"/* t0A */\
		"vmovaps	0x020(%%rsi),%%ymm1	\n\t"/* t0B */\
		"vsubpd	      %%ymm2,%%ymm0,%%ymm0		\n\t"/*~t0A=t0A-rt */\
		"vsubpd	      %%ymm3,%%ymm1,%%ymm1		\n\t"/*~t0B=t0B-it */\
		"vaddpd	     (%%rsi),%%ymm2,%%ymm2	\n\t"/*~t1A=rt+t0A */\
		"vaddpd	0x020(%%rsi),%%ymm3,%%ymm3	\n\t"/*~t1B=it+t0B */\
		"vsubpd	      %%ymm4,%%ymm0,%%ymm0		\n\t"/* t0A-t2A */\
		"vsubpd	      %%ymm7,%%ymm2,%%ymm2		\n\t"/* t1A-t3B */\
		"vsubpd	      %%ymm5,%%ymm1,%%ymm1		\n\t"/* t0B-t2B */\
		"vsubpd	      %%ymm6,%%ymm3,%%ymm3		\n\t"/* t1B-t3A */\
		"vaddpd	      %%ymm4,%%ymm4,%%ymm4		\n\t"/*   2*t2A */\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7		\n\t"/*   2*t3B */\
		"vaddpd	      %%ymm5,%%ymm5,%%ymm5		\n\t"/*   2*t2B */\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6		\n\t"/*   2*t3A */\
		"vmovaps	%%ymm0,     (%%rbx)	\n\t"/* a(jt+p1 ) */\
		"vmovaps	%%ymm2,     (%%rcx)	\n\t"/* a(jt+p2 ) */\
		"vmovaps	%%ymm1,0x020(%%rbx)	\n\t"/* a(jp+p1 ) */\
		"vmovaps	%%ymm3,0x020(%%rdx)	\n\t"/* a(jp+p3 ) */\
		"vaddpd	      %%ymm0,%%ymm4,%%ymm4		\n\t"/* t0A+t2A */\
		"vaddpd	      %%ymm2,%%ymm7,%%ymm7		\n\t"/* t1A+t3B */\
		"vaddpd	      %%ymm1,%%ymm5,%%ymm5		\n\t"/* t0B+t2B */\
		"vaddpd	      %%ymm3,%%ymm6,%%ymm6		\n\t"/* t1B+t3A */\
		"vmovaps	%%ymm4,     (%%rax)	\n\t"/* a(jt+p0 ) */\
		"vmovaps	%%ymm7,     (%%rdx)	\n\t"/* a(jt+p3 ) */\
		"vmovaps	%%ymm5,0x020(%%rax)	\n\t"/* a(jp+p0 ) */\
		"vmovaps	%%ymm6,0x020(%%rcx)	\n\t"/* a(jp+p2 ) */\
		"\n\t"/*...Block 4: t06,t16,t26,t36	*/\
		"addq	$0x300,%%rsi		\n\t"/* r30 */\
		"movslq	%[__p04],%%rdi		\n\t"\
		"shlq	$3,%%rdi			\n\t"\
		"addq	%%rdi,%%rax			\n\t"/* ap18 */\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	%%rdi,%%rcx			\n\t"\
		"addq	%%rdi,%%rdx			\n\t"\
		"movq	%[__cc0],%%rdi	\n\t"\
		"\n\t"\
		"vmovaps	0x040(%%rsi),%%ymm4	\n\t"/* t26 */\
		"vmovaps	0x0c0(%%rsi),%%ymm6	\n\t"/* t36 */\
		"vmovaps	0x060(%%rsi),%%ymm5	\n\t"/* t27 */\
		"vmovaps	0x0e0(%%rsi),%%ymm7	\n\t"/* t37 */\
		"vmovaps	0x040(%%rsi),%%ymm0	\n\t"/* copy t26 */\
		"vmovaps	0x0c0(%%rsi),%%ymm2	\n\t"/* copy t36 */\
		"vmovaps	0x060(%%rsi),%%ymm1	\n\t"/* copy t27 */\
		"vmovaps	0x0e0(%%rsi),%%ymm3	\n\t"/* copy t37 */\
		"vmulpd	0x080(%%rdi),%%ymm4,%%ymm4	\n\t"/* t26*s32_3 */\
		"vmulpd	0x060(%%rdi),%%ymm6,%%ymm6	\n\t"/* t36*s32_1 */\
		"vmulpd	0x0a0(%%rdi),%%ymm1,%%ymm1	\n\t"/* t27*s32_3 */\
		"vmulpd	0x040(%%rdi),%%ymm3,%%ymm3	\n\t"/* t37*c32_1 */\
		"vmulpd	0x080(%%rdi),%%ymm5,%%ymm5	\n\t"/* t27*c32_3 */\
		"vmulpd	0x060(%%rdi),%%ymm7,%%ymm7	\n\t"/* t37*s32_1 */\
		"vmulpd	0x0a0(%%rdi),%%ymm0,%%ymm0	\n\t"/* t26*s32_3 */\
		"vmulpd	0x040(%%rdi),%%ymm2,%%ymm2	\n\t"/* t36*c32_1 */\
		"vsubpd	      %%ymm1,%%ymm4,%%ymm4		\n\t"/* ~t26 */\
		"vaddpd	      %%ymm3,%%ymm6,%%ymm6		\n\t"/* rt */\
		"vaddpd	      %%ymm0,%%ymm5,%%ymm5		\n\t"/* ~t27 */\
		"vsubpd	      %%ymm2,%%ymm7,%%ymm7		\n\t"/* it */\
		"vmovaps	0x080(%%rsi),%%ymm2	\n\t"/* t16 */\
		"vmovaps	0x0a0(%%rsi),%%ymm0	\n\t"/* t17 */\
		"vmovaps	0x080(%%rsi),%%ymm1	\n\t"/* copy t16 */\
		"vmovaps	0x0a0(%%rsi),%%ymm3	\n\t"/* copy t17 */\
		"vsubpd	      %%ymm6,%%ymm4,%%ymm4		\n\t"/*~t26=t26-rt */\
		"vmulpd	0x020(%%rdi),%%ymm2,%%ymm2	\n\t"/* t16*s */\
		"vsubpd	      %%ymm7,%%ymm5,%%ymm5		\n\t"/*~t27=t27-it */\
		"vmulpd	     (%%rdi),%%ymm0,%%ymm0	\n\t"/* t17*c */\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6		\n\t"/*      2* rt */\
		"vmulpd	0x020(%%rdi),%%ymm3,%%ymm3	\n\t"/* t17*s */\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7		\n\t"/*      2* it */\
		"vmulpd	     (%%rdi),%%ymm1,%%ymm1	\n\t"/* t16*c */\
		"vaddpd	      %%ymm4,%%ymm6,%%ymm6		\n\t"/*~t36=t26+rt */\
		"vsubpd	      %%ymm0,%%ymm2,%%ymm2		\n\t"/* rt */\
		"vaddpd	      %%ymm5,%%ymm7,%%ymm7		\n\t"/*~t37=t27+it */\
		"vaddpd	      %%ymm1,%%ymm3,%%ymm3		\n\t"/* it */\
		"vmovaps	     (%%rsi),%%ymm0	\n\t"/* t06 */\
		"vmovaps	0x020(%%rsi),%%ymm1	\n\t"/* t07 */\
		"vsubpd	      %%ymm2,%%ymm0,%%ymm0		\n\t"/*~t16=t06-rt */\
		"vsubpd	      %%ymm3,%%ymm1,%%ymm1		\n\t"/*~t17=t07-it */\
		"vaddpd	     (%%rsi),%%ymm2,%%ymm2	\n\t"/*~t06=rt+t06 */\
		"vaddpd	0x020(%%rsi),%%ymm3,%%ymm3	\n\t"/*~t07=it+t07 */\
		"vsubpd	      %%ymm4,%%ymm2,%%ymm2		\n\t"/* t06-t26 */\
		"vsubpd	      %%ymm7,%%ymm0,%%ymm0		\n\t"/* t16-t37 */\
		"vsubpd	      %%ymm5,%%ymm3,%%ymm3		\n\t"/* t07-t27 */\
		"vsubpd	      %%ymm6,%%ymm1,%%ymm1		\n\t"/* t17-t36 */\
		"vaddpd	      %%ymm4,%%ymm4,%%ymm4		\n\t"/*   2*t26 */\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7		\n\t"/*   2*t37 */\
		"vaddpd	      %%ymm5,%%ymm5,%%ymm5		\n\t"/*   2*t27 */\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6		\n\t"/*   2*t36 */\
		"vmovaps	%%ymm2,     (%%rbx)	\n\t"/* a(jt+p1 ) */\
		"vmovaps	%%ymm0,     (%%rcx)	\n\t"/* a(jt+p2 ) */\
		"vmovaps	%%ymm3,0x020(%%rbx)	\n\t"/* a(jp+p1 ) */\
		"vmovaps	%%ymm1,0x020(%%rdx)	\n\t"/* a(jp+p3 ) */\
		"vaddpd	      %%ymm2,%%ymm4,%%ymm4		\n\t"/* t06+t26 */\
		"vaddpd	      %%ymm0,%%ymm7,%%ymm7		\n\t"/* t16+t37 */\
		"vaddpd	      %%ymm3,%%ymm5,%%ymm5		\n\t"/* t07+t27 */\
		"vaddpd	      %%ymm1,%%ymm6,%%ymm6		\n\t"/* t17+t36 */\
		"vmovaps	%%ymm4,     (%%rax)	\n\t"/* a(jt+p0 ) */\
		"vmovaps	%%ymm7,     (%%rdx)	\n\t"/* a(jt+p3 ) */\
		"vmovaps	%%ymm5,0x020(%%rax)	\n\t"/* a(jp+p0 ) */\
		"vmovaps	%%ymm6,0x020(%%rcx)	\n\t"/* a(jp+p2 ) */\
		"\n\t"/*...Block 8: t0E,t1E,t2E,t3E	*/\
		"addq	$0x100,%%rsi			\n\t"/* r38 */\
		"movslq	%[__p04],%%rdi		\n\t"\
		"shlq	$3,%%rdi			\n\t"\
		"addq	%%rdi,%%rax			\n\t"/* ap1C */\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	%%rdi,%%rcx			\n\t"\
		"addq	%%rdi,%%rdx			\n\t"\
		"movq	%[__cc0],%%rdi	\n\t"\
		"\n\t"\
		"vmovaps	0x040(%%rsi),%%ymm4	\n\t"/* t2E */\
		"vmovaps	0x0c0(%%rsi),%%ymm6	\n\t"/* t3E */\
		"vmovaps	0x060(%%rsi),%%ymm5	\n\t"/* t2F */\
		"vmovaps	0x0e0(%%rsi),%%ymm7	\n\t"/* t3F */\
		"vmovaps	0x040(%%rsi),%%ymm0	\n\t"/* copy t2E */\
		"vmovaps	0x0c0(%%rsi),%%ymm2	\n\t"/* copy t3E */\
		"vmovaps	0x060(%%rsi),%%ymm1	\n\t"/* copy t2F */\
		"vmovaps	0x0e0(%%rsi),%%ymm3	\n\t"/* copy t3F */\
		"vmulpd	0x060(%%rdi),%%ymm4,%%ymm4	\n\t"/* t2E*s32_1 */\
		"vmulpd	0x0a0(%%rdi),%%ymm6,%%ymm6	\n\t"/* t3E*c32_3 */\
		"vmulpd	0x040(%%rdi),%%ymm1,%%ymm1	\n\t"/* t2F*c32_1 */\
		"vmulpd	0x080(%%rdi),%%ymm3,%%ymm3	\n\t"/* t3F*s32_3 */\
		"vmulpd	0x060(%%rdi),%%ymm5,%%ymm5	\n\t"/* t2F*s32_1 */\
		"vmulpd	0x0a0(%%rdi),%%ymm7,%%ymm7	\n\t"/* t3F*c32_3 */\
		"vmulpd	0x040(%%rdi),%%ymm0,%%ymm0	\n\t"/* t2E*c32_1 */\
		"vmulpd	0x080(%%rdi),%%ymm2,%%ymm2	\n\t"/* t3E*s32_3 */\
		"vsubpd	      %%ymm1,%%ymm4,%%ymm4		\n\t"/* ~t2E */\
		"vsubpd	      %%ymm3,%%ymm6,%%ymm6		\n\t"/* rt */\
		"vaddpd	      %%ymm0,%%ymm5,%%ymm5		\n\t"/* ~t2F */\
		"vaddpd	      %%ymm2,%%ymm7,%%ymm7		\n\t"/* it */\
		"vmovaps	0x080(%%rsi),%%ymm2	\n\t"/* t1E */\
		"vmovaps	0x0a0(%%rsi),%%ymm0	\n\t"/* t1F */\
		"vmovaps	0x080(%%rsi),%%ymm1	\n\t"/* copy t1E */\
		"vmovaps	0x0a0(%%rsi),%%ymm3	\n\t"/* copy t1F */\
		"vsubpd	      %%ymm6,%%ymm4,%%ymm4		\n\t"/*~t2E=t2E-rt */\
		"vmulpd	     (%%rdi),%%ymm2,%%ymm2	\n\t"/* t1E*c */\
		"vsubpd	      %%ymm7,%%ymm5,%%ymm5		\n\t"/*~t2F=t2F-it */\
		"vmulpd	0x020(%%rdi),%%ymm0,%%ymm0	\n\t"/* t1F*s */\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6		\n\t"/*      2* rt */\
		"vmulpd	     (%%rdi),%%ymm3,%%ymm3	\n\t"/* t1F*c */\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7		\n\t"/*      2* it */\
		"vmulpd	0x020(%%rdi),%%ymm1,%%ymm1	\n\t"/* t1E*s */\
		"vaddpd	      %%ymm4,%%ymm6,%%ymm6		\n\t"/*~t3E=t2E+rt */\
		"vaddpd	      %%ymm0,%%ymm2,%%ymm2		\n\t"/* rt */\
		"vaddpd	      %%ymm5,%%ymm7,%%ymm7		\n\t"/*~t3F=t2F+it */\
		"vsubpd	      %%ymm1,%%ymm3,%%ymm3		\n\t"/* it */\
		"vmovaps	     (%%rsi),%%ymm0	\n\t"/* t0E */\
		"vmovaps	0x020(%%rsi),%%ymm1	\n\t"/* t0F */\
		"vsubpd	      %%ymm2,%%ymm0,%%ymm0		\n\t"/*~t0E=t0E-rt */\
		"vsubpd	      %%ymm3,%%ymm1,%%ymm1		\n\t"/*~t0F=t0F-it */\
		"vaddpd	     (%%rsi),%%ymm2,%%ymm2	\n\t"/*~t1E=rt+t0E */\
		"vaddpd	0x020(%%rsi),%%ymm3,%%ymm3	\n\t"/*~t1F=it+t0F */\
		"vsubpd	      %%ymm4,%%ymm0,%%ymm0		\n\t"/* t0E-t2E */\
		"vsubpd	      %%ymm7,%%ymm2,%%ymm2		\n\t"/* t1E-t3F */\
		"vsubpd	      %%ymm5,%%ymm1,%%ymm1		\n\t"/* t0F-t2F */\
		"vsubpd	      %%ymm6,%%ymm3,%%ymm3		\n\t"/* t1F-t3E */\
		"vaddpd	      %%ymm4,%%ymm4,%%ymm4		\n\t"/*   2*t2E */\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7		\n\t"/*   2*t3F */\
		"vaddpd	      %%ymm5,%%ymm5,%%ymm5		\n\t"/*   2*t2F */\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6		\n\t"/*   2*t3E */\
		"vmovaps	%%ymm0,     (%%rbx)	\n\t"/* a(jt+p1 ) */\
		"vmovaps	%%ymm2,     (%%rcx)	\n\t"/* a(jt+p2 ) */\
		"vmovaps	%%ymm1,0x020(%%rbx)	\n\t"/* a(jp+p1 ) */\
		"vmovaps	%%ymm3,0x020(%%rdx)	\n\t"/* a(jp+p3 ) */\
		"vaddpd	      %%ymm0,%%ymm4,%%ymm4		\n\t"/* t0E+t2E */\
		"vaddpd	      %%ymm2,%%ymm7,%%ymm7		\n\t"/* t1E+t3F */\
		"vaddpd	      %%ymm1,%%ymm5,%%ymm5		\n\t"/* t0F+t2F */\
		"vaddpd	      %%ymm3,%%ymm6,%%ymm6		\n\t"/* t1F+t3E */\
		"vmovaps	%%ymm4,     (%%rax)	\n\t"/* a(jt+p0 ) */\
		"vmovaps	%%ymm7,     (%%rdx)	\n\t"/* a(jt+p3 ) */\
		"vmovaps	%%ymm5,0x020(%%rax)	\n\t"/* a(jp+p0 ) */\
		"vmovaps	%%ymm6,0x020(%%rcx)	\n\t"/* a(jp+p2 ) */\
		:					/* outputs: none */\
		: [__add] "m" (Xadd)	/* All inputs from memory addresses here */\
		 ,[__p01] "m" (Xp01)\
		 ,[__p02] "m" (Xp02)\
		 ,[__p03] "m" (Xp03)\
		 ,[__p04] "m" (Xp04)\
		 ,[__p08] "m" (Xp08)\
		 ,[__p10] "m" (Xp10)\
		 ,[__p18] "m" (Xp18)\
		 ,[__r00] "m" (Xr00)\
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__cc0] "m" (Xcc0)\
		: "rax","rbx","rcx","rdx","rdi","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","memory"		/* Clobbered registers */\
	);\
	}

  #else // USE_64BIT_ASM_STYLE = True: Deeper 64-bit-ified version of the above 32-bit ASM macros, using all of ymm0-15

	#error USE_64BIT_ASM_STYLE option currently unavailable for this include file!

  #endif	// USE_64BIT_ASM_STYLE ?

#elif defined(USE_SSE2)

  #define USE_64BIT_ASM_STYLE	0	// My x86 timings indicate fancier versions below using xmm0-15 for the radix-32 DFT is faster.

  #if !USE_64BIT_ASM_STYLE	// False: Use simple 64-bit-ified version of the analogous 32-bit ASM macros, i.e. using just xmm0-7.

	#define SSE2_RADIX32_DIT_NOTWIDDLE(Xadd,Xp01,Xp02,Xp03,Xp04,Xp08,Xp10,Xp18,Xr00,Xisrt2,Xcc0)\
	{\
	__asm__ volatile (\
		"\n\t"/*...Block 1: */\
		"movq	%[__r00],%%rsi	\n\t"\
		"movq	%[__add],%%rax	\n\t"\
		"movslq	%[__p01],%%rbx	\n\t"\
		"movslq	%[__p02],%%rcx	\n\t"\
		"movslq	%[__p03],%%rdx	\n\t"\
		"shlq	$3,%%rbx		\n\t"\
		"shlq	$3,%%rcx		\n\t"\
		"shlq	$3,%%rdx		\n\t"\
		"addq	%%rax,%%rbx		\n\t"\
		"addq	%%rax,%%rcx		\n\t"\
		"addq	%%rax,%%rdx		\n\t"\
		"\n\t"/* SSE2_RADIX4_DIT_0TWIDDLE(add0, add1, add2, add3, r00): */\
		"movaps	     (%%rax),%%xmm2	\n\t"\
		"movaps	     (%%rcx),%%xmm6	\n\t"\
		"movaps	0x010(%%rax),%%xmm3	\n\t"\
		"movaps	0x010(%%rcx),%%xmm7	\n\t"\
		"movaps	     (%%rbx),%%xmm0	\n\t"\
		"movaps	     (%%rdx),%%xmm4	\n\t"\
		"movaps	0x010(%%rbx),%%xmm1	\n\t"\
		"movaps	0x010(%%rdx),%%xmm5	\n\t"\
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
		"\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm2		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"subpd	%%xmm6,%%xmm3		\n\t"\
		"movaps	%%xmm0,0x040(%%rsi)	\n\t"\
		"movaps	%%xmm2,0x060(%%rsi)	\n\t"\
		"movaps	%%xmm1,0x050(%%rsi)	\n\t"\
		"movaps	%%xmm3,0x030(%%rsi)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"addpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	%%xmm4,     (%%rsi)	\n\t"\
		"movaps	%%xmm7,0x020(%%rsi)	\n\t"\
		"movaps	%%xmm5,0x010(%%rsi)	\n\t"\
		"movaps	%%xmm6,0x070(%%rsi)	\n\t"\
		"\n\t"/* SSE2_RADIX4_DIT_0TWIDDLE_2NDOFTWO(add4, add5, add6, add7, r08): */\
		"addq	$0x080,%%rsi			\n\t"\
		"movslq	%[__p04],%%rdi		\n\t"\
		"shlq	$3,%%rdi			\n\t"\
		"addq	%%rdi,%%rax			\n\t"\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	%%rdi,%%rcx			\n\t"\
		"addq	%%rdi,%%rdx			\n\t"\
		"\n\t"\
		"movaps	     (%%rax),%%xmm2	\n\t"\
		"movaps	     (%%rcx),%%xmm6	\n\t"\
		"movaps	0x010(%%rax),%%xmm3	\n\t"\
		"movaps	0x010(%%rcx),%%xmm7	\n\t"\
		"movaps	     (%%rbx),%%xmm0	\n\t"\
		"movaps	     (%%rdx),%%xmm4	\n\t"\
		"movaps	0x010(%%rbx),%%xmm1	\n\t"\
		"movaps	0x010(%%rdx),%%xmm5	\n\t"\
		"\n\t"\
		"movq	%[__isrt2],%%rdi	\n\t"\
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
		"\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"movaps	%%xmm0,0x040(%%rsi)	\n\t"\
		"movaps	%%xmm1,0x050(%%rsi)	\n\t"\
		"subpd	%%xmm7,%%xmm2		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"subpd	%%xmm6,%%xmm3		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"movaps	%%xmm4,     (%%rsi)	\n\t"\
		"addpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	%%xmm5,0x010(%%rsi)	\n\t"\
		"\n\t"\
		"movaps	(%%rdi),%%xmm5	\n\t"/* ISRT2 */\
		"movaps	%%xmm3,%%xmm0		\n\t"\
		"movaps	%%xmm6,%%xmm1		\n\t"\
		"subpd	%%xmm7,%%xmm3		\n\t"\
		"subpd	%%xmm2,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm0		\n\t"\
		"addpd	%%xmm2,%%xmm1		\n\t"\
		"mulpd	%%xmm5,%%xmm3		\n\t"\
		"mulpd	%%xmm5,%%xmm6		\n\t"\
		"mulpd	%%xmm5,%%xmm0		\n\t"\
		"mulpd	%%xmm5,%%xmm1		\n\t"\
		"movaps	%%xmm3,0x030(%%rsi)	\n\t"\
		"movaps	%%xmm6,0x070(%%rsi)	\n\t"\
		"movaps	%%xmm0,0x020(%%rsi)	\n\t"\
		"movaps	%%xmm1,0x060(%%rsi)	\n\t"\
		"\n\t"/* SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS(r00,r02,r04,r06,r08,r0A,r0C,r0E): */\
		"movaps		-0x080(%%rsi),%%xmm0	\n\t"/* r00 */\
		"movaps		-0x040(%%rsi),%%xmm4	\n\t"/* r04 */\
		"movaps		-0x070(%%rsi),%%xmm1	\n\t"/* r01 */\
		"movaps		-0x030(%%rsi),%%xmm5	\n\t"/* r05 */\
		"movaps		      (%%rsi),%%xmm2	\n\t"/* r08 */\
		"movaps		 0x050(%%rsi),%%xmm7	\n\t"/* r0D */\
		"movaps		 0x010(%%rsi),%%xmm3	\n\t"/* r09 */\
		"movaps		 0x040(%%rsi),%%xmm6	\n\t"/* r0C */\
		"subpd		%%xmm2,%%xmm0\n\t"\
		"subpd		%%xmm7,%%xmm4\n\t"\
		"subpd		%%xmm3,%%xmm1\n\t"\
		"subpd		%%xmm6,%%xmm5\n\t"\
		"addpd		%%xmm2,%%xmm2\n\t"\
		"addpd		%%xmm7,%%xmm7\n\t"\
		"addpd		%%xmm3,%%xmm3\n\t"\
		"addpd		%%xmm6,%%xmm6\n\t"\
		"addpd		%%xmm0,%%xmm2\n\t"\
		"addpd		%%xmm4,%%xmm7\n\t"\
		"addpd		%%xmm1,%%xmm3\n\t"\
		"addpd		%%xmm5,%%xmm6\n\t"\
		"movaps		%%xmm0,      (%%rsi)\n\t"\
		"movaps		%%xmm4, 0x040(%%rsi)\n\t"\
		"movaps		%%xmm1, 0x010(%%rsi)\n\t"\
		"movaps		%%xmm5,-0x030(%%rsi)\n\t"\
		"movaps		%%xmm2,-0x080(%%rsi)\n\t"\
		"movaps		%%xmm7,-0x040(%%rsi)\n\t"\
		"movaps		%%xmm3,-0x070(%%rsi)\n\t"\
		"movaps		%%xmm6, 0x050(%%rsi)\n\t"\
		"\n\t"\
		"movaps		-0x060(%%rsi),%%xmm0	\n\t"/* r02 */\
		"movaps		-0x020(%%rsi),%%xmm4	\n\t"/* r06 */\
		"movaps		-0x050(%%rsi),%%xmm1	\n\t"/* r03 */\
		"movaps		-0x010(%%rsi),%%xmm5	\n\t"/* r07 */\
		"movaps		 0x020(%%rsi),%%xmm2	\n\t"/* r0A */\
		"movaps		 0x070(%%rsi),%%xmm7	\n\t"/* r0F */\
		"movaps		 0x030(%%rsi),%%xmm3	\n\t"/* r0B */\
		"movaps		 0x060(%%rsi),%%xmm6	\n\t"/* r0E */\
		"subpd		%%xmm2,%%xmm0\n\t"\
		"subpd		%%xmm7,%%xmm4\n\t"\
		"subpd		%%xmm3,%%xmm1\n\t"\
		"subpd		%%xmm6,%%xmm5\n\t"\
		"addpd		%%xmm2,%%xmm2\n\t"\
		"addpd		%%xmm7,%%xmm7\n\t"\
		"addpd		%%xmm3,%%xmm3\n\t"\
		"addpd		%%xmm6,%%xmm6\n\t"\
		"addpd		%%xmm0,%%xmm2\n\t"\
		"addpd		%%xmm4,%%xmm7\n\t"\
		"addpd		%%xmm1,%%xmm3\n\t"\
		"addpd		%%xmm5,%%xmm6\n\t"\
		"movaps		%%xmm0, 0x020(%%rsi)\n\t"\
		"movaps		%%xmm4, 0x060(%%rsi)\n\t"\
		"movaps		%%xmm1, 0x030(%%rsi)\n\t"\
		"movaps		%%xmm5,-0x010(%%rsi)\n\t"\
		"movaps		%%xmm2,-0x060(%%rsi)\n\t"\
		"movaps		%%xmm7,-0x020(%%rsi)\n\t"\
		"movaps		%%xmm3,-0x050(%%rsi)\n\t"\
		"movaps		%%xmm6, 0x070(%%rsi)\n\t"\
		"\n\t"\
		"\n\t"/*...Block 2: */\
		"addq	$0x080,%%rsi			\n\t"\
		"movslq	%[__p04],%%rdi		\n\t"\
		"shlq	$3,%%rdi			\n\t"\
		"addq	%%rdi,%%rax			\n\t"\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	%%rdi,%%rcx			\n\t"\
		"addq	%%rdi,%%rdx			\n\t"\
		"\n\t"\
		"\n\t"/* SSE2_RADIX4_DIT_0TWIDDLE(add0, add1, add2, add3, r10): */\
		"movaps	     (%%rax),%%xmm2	\n\t"\
		"movaps	     (%%rcx),%%xmm6	\n\t"\
		"movaps	0x010(%%rax),%%xmm3	\n\t"\
		"movaps	0x010(%%rcx),%%xmm7	\n\t"\
		"movaps	     (%%rbx),%%xmm0	\n\t"\
		"movaps	     (%%rdx),%%xmm4	\n\t"\
		"movaps	0x010(%%rbx),%%xmm1	\n\t"\
		"movaps	0x010(%%rdx),%%xmm5	\n\t"\
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
		"\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm2		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"subpd	%%xmm6,%%xmm3		\n\t"\
		"movaps	%%xmm0,0x040(%%rsi)	\n\t"\
		"movaps	%%xmm2,0x060(%%rsi)	\n\t"\
		"movaps	%%xmm1,0x050(%%rsi)	\n\t"\
		"movaps	%%xmm3,0x030(%%rsi)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"addpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	%%xmm4,     (%%rsi)	\n\t"\
		"movaps	%%xmm7,0x020(%%rsi)	\n\t"\
		"movaps	%%xmm5,0x010(%%rsi)	\n\t"\
		"movaps	%%xmm6,0x070(%%rsi)	\n\t"\
		"\n\t"/* SSE2_RADIX4_DIT_0TWIDDLE_2NDOFTWO(add4, add5, add6, add7, r18): */\
		"addq	$0x080,%%rsi			\n\t"\
		"movslq	%[__p04],%%rdi		\n\t"\
		"shlq	$3,%%rdi			\n\t"\
		"addq	%%rdi,%%rax			\n\t"\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	%%rdi,%%rcx			\n\t"\
		"addq	%%rdi,%%rdx			\n\t"\
		"\n\t"\
		"movaps	     (%%rax),%%xmm2	\n\t"\
		"movaps	     (%%rcx),%%xmm6	\n\t"\
		"movaps	0x010(%%rax),%%xmm3	\n\t"\
		"movaps	0x010(%%rcx),%%xmm7	\n\t"\
		"movaps	     (%%rbx),%%xmm0	\n\t"\
		"movaps	     (%%rdx),%%xmm4	\n\t"\
		"movaps	0x010(%%rbx),%%xmm1	\n\t"\
		"movaps	0x010(%%rdx),%%xmm5	\n\t"\
		"\n\t"\
		"movq	%[__isrt2],%%rdi	\n\t"\
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
		"\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"movaps	%%xmm0,0x040(%%rsi)	\n\t"\
		"movaps	%%xmm1,0x050(%%rsi)	\n\t"\
		"subpd	%%xmm7,%%xmm2		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"subpd	%%xmm6,%%xmm3		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"movaps	%%xmm4,     (%%rsi)	\n\t"\
		"addpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	%%xmm5,0x010(%%rsi)	\n\t"\
		"\n\t"\
		"movaps	(%%rdi),%%xmm5	\n\t"/* ISRT2 */\
		"movaps	%%xmm3,%%xmm0		\n\t"\
		"movaps	%%xmm6,%%xmm1		\n\t"\
		"subpd	%%xmm7,%%xmm3		\n\t"\
		"subpd	%%xmm2,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm0		\n\t"\
		"addpd	%%xmm2,%%xmm1		\n\t"\
		"mulpd	%%xmm5,%%xmm3		\n\t"\
		"mulpd	%%xmm5,%%xmm6		\n\t"\
		"mulpd	%%xmm5,%%xmm0		\n\t"\
		"mulpd	%%xmm5,%%xmm1		\n\t"\
		"movaps	%%xmm3,0x030(%%rsi)	\n\t"\
		"movaps	%%xmm6,0x070(%%rsi)	\n\t"\
		"movaps	%%xmm0,0x020(%%rsi)	\n\t"\
		"movaps	%%xmm1,0x060(%%rsi)	\n\t"\
		"\n\t"/* SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS(r10,r12,r14,r16,r18,r1A,r1C,r1E): */\
		"movaps		-0x080(%%rsi),%%xmm0	\n\t"/* r10 */\
		"movaps		-0x040(%%rsi),%%xmm4	\n\t"/* r14 */\
		"movaps		-0x070(%%rsi),%%xmm1	\n\t"/* r11 */\
		"movaps		-0x030(%%rsi),%%xmm5	\n\t"/* r15 */\
		"movaps		      (%%rsi),%%xmm2	\n\t"/* r18 */\
		"movaps		 0x050(%%rsi),%%xmm7	\n\t"/* r1D */\
		"movaps		 0x010(%%rsi),%%xmm3	\n\t"/* r19 */\
		"movaps		 0x040(%%rsi),%%xmm6	\n\t"/* r1C */\
		"subpd		%%xmm2,%%xmm0\n\t"\
		"subpd		%%xmm7,%%xmm4\n\t"\
		"subpd		%%xmm3,%%xmm1\n\t"\
		"subpd		%%xmm6,%%xmm5\n\t"\
		"addpd		%%xmm2,%%xmm2\n\t"\
		"addpd		%%xmm7,%%xmm7\n\t"\
		"addpd		%%xmm3,%%xmm3\n\t"\
		"addpd		%%xmm6,%%xmm6\n\t"\
		"addpd		%%xmm0,%%xmm2\n\t"\
		"addpd		%%xmm4,%%xmm7\n\t"\
		"addpd		%%xmm1,%%xmm3\n\t"\
		"addpd		%%xmm5,%%xmm6\n\t"\
		"movaps		%%xmm0,      (%%rsi)\n\t"\
		"movaps		%%xmm4, 0x040(%%rsi)\n\t"\
		"movaps		%%xmm1, 0x010(%%rsi)\n\t"\
		"movaps		%%xmm5,-0x030(%%rsi)\n\t"\
		"movaps		%%xmm2,-0x080(%%rsi)\n\t"\
		"movaps		%%xmm7,-0x040(%%rsi)\n\t"\
		"movaps		%%xmm3,-0x070(%%rsi)\n\t"\
		"movaps		%%xmm6, 0x050(%%rsi)\n\t"\
		"\n\t"\
		"movaps		-0x060(%%rsi),%%xmm0	\n\t"/* r12 */\
		"movaps		-0x020(%%rsi),%%xmm4	\n\t"/* r16 */\
		"movaps		-0x050(%%rsi),%%xmm1	\n\t"/* r13 */\
		"movaps		-0x010(%%rsi),%%xmm5	\n\t"/* r17 */\
		"movaps		 0x020(%%rsi),%%xmm2	\n\t"/* r1A */\
		"movaps		 0x070(%%rsi),%%xmm7	\n\t"/* r1F */\
		"movaps		 0x030(%%rsi),%%xmm3	\n\t"/* r1B */\
		"movaps		 0x060(%%rsi),%%xmm6	\n\t"/* r1E */\
		"subpd		%%xmm2,%%xmm0\n\t"\
		"subpd		%%xmm7,%%xmm4\n\t"\
		"subpd		%%xmm3,%%xmm1\n\t"\
		"subpd		%%xmm6,%%xmm5\n\t"\
		"addpd		%%xmm2,%%xmm2\n\t"\
		"addpd		%%xmm7,%%xmm7\n\t"\
		"addpd		%%xmm3,%%xmm3\n\t"\
		"addpd		%%xmm6,%%xmm6\n\t"\
		"addpd		%%xmm0,%%xmm2\n\t"\
		"addpd		%%xmm4,%%xmm7\n\t"\
		"addpd		%%xmm1,%%xmm3\n\t"\
		"addpd		%%xmm5,%%xmm6\n\t"\
		"movaps		%%xmm0, 0x020(%%rsi)\n\t"\
		"movaps		%%xmm4, 0x060(%%rsi)\n\t"\
		"movaps		%%xmm1, 0x030(%%rsi)\n\t"\
		"movaps		%%xmm5,-0x010(%%rsi)\n\t"\
		"movaps		%%xmm2,-0x060(%%rsi)\n\t"\
		"movaps		%%xmm7,-0x020(%%rsi)\n\t"\
		"movaps		%%xmm3,-0x050(%%rsi)\n\t"\
		"movaps		%%xmm6, 0x070(%%rsi)\n\t"\
		"\n\t"\
		"\n\t"/*...Block 3: */\
		"addq	$0x080,%%rsi			\n\t"\
		"movslq	%[__p04],%%rdi		\n\t"\
		"shlq	$3,%%rdi			\n\t"\
		"addq	%%rdi,%%rax			\n\t"\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	%%rdi,%%rcx			\n\t"\
		"addq	%%rdi,%%rdx			\n\t"\
		"\n\t"\
		"\n\t"/* SSE2_RADIX4_DIT_0TWIDDLE(add0, add1, add2, add3, r20): */\
		"movaps	     (%%rax),%%xmm2	\n\t"\
		"movaps	     (%%rcx),%%xmm6	\n\t"\
		"movaps	0x010(%%rax),%%xmm3	\n\t"\
		"movaps	0x010(%%rcx),%%xmm7	\n\t"\
		"movaps	     (%%rbx),%%xmm0	\n\t"\
		"movaps	     (%%rdx),%%xmm4	\n\t"\
		"movaps	0x010(%%rbx),%%xmm1	\n\t"\
		"movaps	0x010(%%rdx),%%xmm5	\n\t"\
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
		"\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm2		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"subpd	%%xmm6,%%xmm3		\n\t"\
		"movaps	%%xmm0,0x040(%%rsi)	\n\t"\
		"movaps	%%xmm2,0x060(%%rsi)	\n\t"\
		"movaps	%%xmm1,0x050(%%rsi)	\n\t"\
		"movaps	%%xmm3,0x030(%%rsi)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"addpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	%%xmm4,     (%%rsi)	\n\t"\
		"movaps	%%xmm7,0x020(%%rsi)	\n\t"\
		"movaps	%%xmm5,0x010(%%rsi)	\n\t"\
		"movaps	%%xmm6,0x070(%%rsi)	\n\t"\
		"\n\t"/* SSE2_RADIX4_DIT_0TWIDDLE_2NDOFTWO(add4, add5, add6, add7, r28): */\
		"addq	$0x080,%%rsi			\n\t"\
		"movslq	%[__p04],%%rdi		\n\t"\
		"shlq	$3,%%rdi			\n\t"\
		"addq	%%rdi,%%rax			\n\t"\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	%%rdi,%%rcx			\n\t"\
		"addq	%%rdi,%%rdx			\n\t"\
		"\n\t"\
		"movaps	     (%%rax),%%xmm2	\n\t"\
		"movaps	     (%%rcx),%%xmm6	\n\t"\
		"movaps	0x010(%%rax),%%xmm3	\n\t"\
		"movaps	0x010(%%rcx),%%xmm7	\n\t"\
		"movaps	     (%%rbx),%%xmm0	\n\t"\
		"movaps	     (%%rdx),%%xmm4	\n\t"\
		"movaps	0x010(%%rbx),%%xmm1	\n\t"\
		"movaps	0x010(%%rdx),%%xmm5	\n\t"\
		"\n\t"\
		"movq	%[__isrt2],%%rdi	\n\t"\
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
		"\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"movaps	%%xmm0,0x040(%%rsi)	\n\t"\
		"movaps	%%xmm1,0x050(%%rsi)	\n\t"\
		"subpd	%%xmm7,%%xmm2		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"subpd	%%xmm6,%%xmm3		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"movaps	%%xmm4,     (%%rsi)	\n\t"\
		"addpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	%%xmm5,0x010(%%rsi)	\n\t"\
		"\n\t"\
		"movaps	(%%rdi),%%xmm5	\n\t"/* ISRT2 */\
		"movaps	%%xmm3,%%xmm0		\n\t"\
		"movaps	%%xmm6,%%xmm1		\n\t"\
		"subpd	%%xmm7,%%xmm3		\n\t"\
		"subpd	%%xmm2,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm0		\n\t"\
		"addpd	%%xmm2,%%xmm1		\n\t"\
		"mulpd	%%xmm5,%%xmm3		\n\t"\
		"mulpd	%%xmm5,%%xmm6		\n\t"\
		"mulpd	%%xmm5,%%xmm0		\n\t"\
		"mulpd	%%xmm5,%%xmm1		\n\t"\
		"movaps	%%xmm3,0x030(%%rsi)	\n\t"\
		"movaps	%%xmm6,0x070(%%rsi)	\n\t"\
		"movaps	%%xmm0,0x020(%%rsi)	\n\t"\
		"movaps	%%xmm1,0x060(%%rsi)	\n\t"\
		"\n\t"/* SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS(r20,r22,r24,r26,r28,r2A,r2C,r2E): */\
		"movaps		-0x080(%%rsi),%%xmm0	\n\t"/* r20 */\
		"movaps		-0x040(%%rsi),%%xmm4	\n\t"/* r24 */\
		"movaps		-0x070(%%rsi),%%xmm1	\n\t"/* r21 */\
		"movaps		-0x030(%%rsi),%%xmm5	\n\t"/* r25 */\
		"movaps		      (%%rsi),%%xmm2	\n\t"/* r28 */\
		"movaps		 0x050(%%rsi),%%xmm7	\n\t"/* r2D */\
		"movaps		 0x010(%%rsi),%%xmm3	\n\t"/* r29 */\
		"movaps		 0x040(%%rsi),%%xmm6	\n\t"/* r2C */\
		"subpd		%%xmm2,%%xmm0\n\t"\
		"subpd		%%xmm7,%%xmm4\n\t"\
		"subpd		%%xmm3,%%xmm1\n\t"\
		"subpd		%%xmm6,%%xmm5\n\t"\
		"addpd		%%xmm2,%%xmm2\n\t"\
		"addpd		%%xmm7,%%xmm7\n\t"\
		"addpd		%%xmm3,%%xmm3\n\t"\
		"addpd		%%xmm6,%%xmm6\n\t"\
		"addpd		%%xmm0,%%xmm2\n\t"\
		"addpd		%%xmm4,%%xmm7\n\t"\
		"addpd		%%xmm1,%%xmm3\n\t"\
		"addpd		%%xmm5,%%xmm6\n\t"\
		"movaps		%%xmm0,      (%%rsi)\n\t"\
		"movaps		%%xmm4, 0x040(%%rsi)\n\t"\
		"movaps		%%xmm1, 0x010(%%rsi)\n\t"\
		"movaps		%%xmm5,-0x030(%%rsi)\n\t"\
		"movaps		%%xmm2,-0x080(%%rsi)\n\t"\
		"movaps		%%xmm7,-0x040(%%rsi)\n\t"\
		"movaps		%%xmm3,-0x070(%%rsi)\n\t"\
		"movaps		%%xmm6, 0x050(%%rsi)\n\t"\
		"\n\t"\
		"movaps		-0x060(%%rsi),%%xmm0	\n\t"/* r22 */\
		"movaps		-0x020(%%rsi),%%xmm4	\n\t"/* r26 */\
		"movaps		-0x050(%%rsi),%%xmm1	\n\t"/* r23 */\
		"movaps		-0x010(%%rsi),%%xmm5	\n\t"/* r27 */\
		"movaps		 0x020(%%rsi),%%xmm2	\n\t"/* r2A */\
		"movaps		 0x070(%%rsi),%%xmm7	\n\t"/* r2F */\
		"movaps		 0x030(%%rsi),%%xmm3	\n\t"/* r2B */\
		"movaps		 0x060(%%rsi),%%xmm6	\n\t"/* r2E */\
		"subpd		%%xmm2,%%xmm0\n\t"\
		"subpd		%%xmm7,%%xmm4\n\t"\
		"subpd		%%xmm3,%%xmm1\n\t"\
		"subpd		%%xmm6,%%xmm5\n\t"\
		"addpd		%%xmm2,%%xmm2\n\t"\
		"addpd		%%xmm7,%%xmm7\n\t"\
		"addpd		%%xmm3,%%xmm3\n\t"\
		"addpd		%%xmm6,%%xmm6\n\t"\
		"addpd		%%xmm0,%%xmm2\n\t"\
		"addpd		%%xmm4,%%xmm7\n\t"\
		"addpd		%%xmm1,%%xmm3\n\t"\
		"addpd		%%xmm5,%%xmm6\n\t"\
		"movaps		%%xmm0, 0x020(%%rsi)\n\t"\
		"movaps		%%xmm4, 0x060(%%rsi)\n\t"\
		"movaps		%%xmm1, 0x030(%%rsi)\n\t"\
		"movaps		%%xmm5,-0x010(%%rsi)\n\t"\
		"movaps		%%xmm2,-0x060(%%rsi)\n\t"\
		"movaps		%%xmm7,-0x020(%%rsi)\n\t"\
		"movaps		%%xmm3,-0x050(%%rsi)\n\t"\
		"movaps		%%xmm6, 0x070(%%rsi)\n\t"\
		"\n\t"\
		"\n\t"/*...Block 4: */\
		"addq	$0x080,%%rsi			\n\t"\
		"movslq	%[__p04],%%rdi		\n\t"\
		"shlq	$3,%%rdi			\n\t"\
		"addq	%%rdi,%%rax			\n\t"\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	%%rdi,%%rcx			\n\t"\
		"addq	%%rdi,%%rdx			\n\t"\
		"\n\t"\
		"\n\t"/* SSE2_RADIX4_DIT_0TWIDDLE(add0, add1, add2, add3, r30): */\
		"movaps	     (%%rax),%%xmm2	\n\t"\
		"movaps	     (%%rcx),%%xmm6	\n\t"\
		"movaps	0x010(%%rax),%%xmm3	\n\t"\
		"movaps	0x010(%%rcx),%%xmm7	\n\t"\
		"movaps	     (%%rbx),%%xmm0	\n\t"\
		"movaps	     (%%rdx),%%xmm4	\n\t"\
		"movaps	0x010(%%rbx),%%xmm1	\n\t"\
		"movaps	0x010(%%rdx),%%xmm5	\n\t"\
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
		"\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm2		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"subpd	%%xmm6,%%xmm3		\n\t"\
		"movaps	%%xmm0,0x040(%%rsi)	\n\t"\
		"movaps	%%xmm2,0x060(%%rsi)	\n\t"\
		"movaps	%%xmm1,0x050(%%rsi)	\n\t"\
		"movaps	%%xmm3,0x030(%%rsi)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"addpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	%%xmm4,     (%%rsi)	\n\t"\
		"movaps	%%xmm7,0x020(%%rsi)	\n\t"\
		"movaps	%%xmm5,0x010(%%rsi)	\n\t"\
		"movaps	%%xmm6,0x070(%%rsi)	\n\t"\
		"\n\t"/* SSE2_RADIX4_DIT_0TWIDDLE_2NDOFTWO(add4, add5, add6, add7, r38): */\
		"addq	$0x080,%%rsi			\n\t"\
		"movslq	%[__p04],%%rdi		\n\t"\
		"shlq	$3,%%rdi			\n\t"\
		"addq	%%rdi,%%rax			\n\t"\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	%%rdi,%%rcx			\n\t"\
		"addq	%%rdi,%%rdx			\n\t"\
		"\n\t"\
		"movaps	     (%%rax),%%xmm2	\n\t"\
		"movaps	     (%%rcx),%%xmm6	\n\t"\
		"movaps	0x010(%%rax),%%xmm3	\n\t"\
		"movaps	0x010(%%rcx),%%xmm7	\n\t"\
		"movaps	     (%%rbx),%%xmm0	\n\t"\
		"movaps	     (%%rdx),%%xmm4	\n\t"\
		"movaps	0x010(%%rbx),%%xmm1	\n\t"\
		"movaps	0x010(%%rdx),%%xmm5	\n\t"\
		"\n\t"\
		"movq	%[__isrt2],%%rdi	\n\t"\
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
		"\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"movaps	%%xmm0,0x040(%%rsi)	\n\t"\
		"movaps	%%xmm1,0x050(%%rsi)	\n\t"\
		"subpd	%%xmm7,%%xmm2		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"subpd	%%xmm6,%%xmm3		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"movaps	%%xmm4,     (%%rsi)	\n\t"\
		"addpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	%%xmm5,0x010(%%rsi)	\n\t"\
		"\n\t"\
		"movaps	(%%rdi),%%xmm5	\n\t"/* ISRT2 */\
		"movaps	%%xmm3,%%xmm0		\n\t"\
		"movaps	%%xmm6,%%xmm1		\n\t"\
		"subpd	%%xmm7,%%xmm3		\n\t"\
		"subpd	%%xmm2,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm0		\n\t"\
		"addpd	%%xmm2,%%xmm1		\n\t"\
		"mulpd	%%xmm5,%%xmm3		\n\t"\
		"mulpd	%%xmm5,%%xmm6		\n\t"\
		"mulpd	%%xmm5,%%xmm0		\n\t"\
		"mulpd	%%xmm5,%%xmm1		\n\t"\
		"movaps	%%xmm3,0x030(%%rsi)	\n\t"\
		"movaps	%%xmm6,0x070(%%rsi)	\n\t"\
		"movaps	%%xmm0,0x020(%%rsi)	\n\t"\
		"movaps	%%xmm1,0x060(%%rsi)	\n\t"\
		"\n\t"/* SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS(r30,r32,r34,r36,r38,r3A,r3C,r3E): */\
		"movaps		-0x080(%%rsi),%%xmm0	\n\t"/* r30 */\
		"movaps		-0x040(%%rsi),%%xmm4	\n\t"/* r34 */\
		"movaps		-0x070(%%rsi),%%xmm1	\n\t"/* r31 */\
		"movaps		-0x030(%%rsi),%%xmm5	\n\t"/* r35 */\
		"movaps		      (%%rsi),%%xmm2	\n\t"/* r38 */\
		"movaps		 0x050(%%rsi),%%xmm7	\n\t"/* r3D */\
		"movaps		 0x010(%%rsi),%%xmm3	\n\t"/* r39 */\
		"movaps		 0x040(%%rsi),%%xmm6	\n\t"/* r3C */\
		"subpd		%%xmm2,%%xmm0\n\t"\
		"subpd		%%xmm7,%%xmm4\n\t"\
		"subpd		%%xmm3,%%xmm1\n\t"\
		"subpd		%%xmm6,%%xmm5\n\t"\
		"addpd		%%xmm2,%%xmm2\n\t"\
		"addpd		%%xmm7,%%xmm7\n\t"\
		"addpd		%%xmm3,%%xmm3\n\t"\
		"addpd		%%xmm6,%%xmm6\n\t"\
		"addpd		%%xmm0,%%xmm2\n\t"\
		"addpd		%%xmm4,%%xmm7\n\t"\
		"addpd		%%xmm1,%%xmm3\n\t"\
		"addpd		%%xmm5,%%xmm6\n\t"\
		"movaps		%%xmm0,      (%%rsi)\n\t"\
		"movaps		%%xmm4, 0x040(%%rsi)\n\t"\
		"movaps		%%xmm1, 0x010(%%rsi)\n\t"\
		"movaps		%%xmm5,-0x030(%%rsi)\n\t"\
		"movaps		%%xmm2,-0x080(%%rsi)\n\t"\
		"movaps		%%xmm7,-0x040(%%rsi)\n\t"\
		"movaps		%%xmm3,-0x070(%%rsi)\n\t"\
		"movaps		%%xmm6, 0x050(%%rsi)\n\t"\
		"\n\t"\
		"movaps		-0x060(%%rsi),%%xmm0	\n\t"/* r32 */\
		"movaps		-0x020(%%rsi),%%xmm4	\n\t"/* r36 */\
		"movaps		-0x050(%%rsi),%%xmm1	\n\t"/* r33 */\
		"movaps		-0x010(%%rsi),%%xmm5	\n\t"/* r37 */\
		"movaps		 0x020(%%rsi),%%xmm2	\n\t"/* r3A */\
		"movaps		 0x070(%%rsi),%%xmm7	\n\t"/* r3F */\
		"movaps		 0x030(%%rsi),%%xmm3	\n\t"/* r3B */\
		"movaps		 0x060(%%rsi),%%xmm6	\n\t"/* r3E */\
		"subpd		%%xmm2,%%xmm0\n\t"\
		"subpd		%%xmm7,%%xmm4\n\t"\
		"subpd		%%xmm3,%%xmm1\n\t"\
		"subpd		%%xmm6,%%xmm5\n\t"\
		"addpd		%%xmm2,%%xmm2\n\t"\
		"addpd		%%xmm7,%%xmm7\n\t"\
		"addpd		%%xmm3,%%xmm3\n\t"\
		"addpd		%%xmm6,%%xmm6\n\t"\
		"addpd		%%xmm0,%%xmm2\n\t"\
		"addpd		%%xmm4,%%xmm7\n\t"\
		"addpd		%%xmm1,%%xmm3\n\t"\
		"addpd		%%xmm5,%%xmm6\n\t"\
		"movaps		%%xmm0, 0x020(%%rsi)\n\t"\
		"movaps		%%xmm4, 0x060(%%rsi)\n\t"\
		"movaps		%%xmm1, 0x030(%%rsi)\n\t"\
		"movaps		%%xmm5,-0x010(%%rsi)\n\t"\
		"movaps		%%xmm2,-0x060(%%rsi)\n\t"\
		"movaps		%%xmm7,-0x020(%%rsi)\n\t"\
		"movaps		%%xmm3,-0x050(%%rsi)\n\t"\
		"movaps		%%xmm6, 0x070(%%rsi)\n\t"\
		"\n\t"\
	"\n\t"/*...and now do eight radix-4 transforms, including the internal twiddle factors:	*/\
		"\n\t"/*...Block 1: r00,r10,r20,r30	*/\
		"\n\t"\
		"movq	%[__isrt2],%%rdi	\n\t"\
		"movq	%[__cc0],%%rsi	\n\t"\
		"movq	%[__r00],%%rax	\n\t"\
		"movq	%%rax,%%rbx		\n\t"\
		"movq	%%rax,%%rcx		\n\t"\
		"movq	%%rax,%%rdx		\n\t"\
		"addq	$0x100,%%rbx	\n\t"\
		"addq	$0x200,%%rcx	\n\t"\
		"addq	$0x300,%%rdx	\n\t"\
		"\n\t"\
		"movaps	     (%%rax),%%xmm0		\n\t"\
		"movaps	0x010(%%rax),%%xmm1		\n\t"\
		"movaps	     (%%rbx),%%xmm2		\n\t"\
		"movaps	0x010(%%rbx),%%xmm3		\n\t"\
		"subpd	     (%%rbx),%%xmm0		\n\t"\
		"subpd	0x010(%%rbx),%%xmm1		\n\t"\
		"addpd	     (%%rax),%%xmm2		\n\t"\
		"addpd	0x010(%%rax),%%xmm3		\n\t"\
		"movaps	     (%%rcx),%%xmm4		\n\t"\
		"movaps	0x010(%%rcx),%%xmm5		\n\t"\
		"movaps	     (%%rdx),%%xmm6		\n\t"\
		"movaps	0x010(%%rdx),%%xmm7		\n\t"\
		"subpd	     (%%rdx),%%xmm4		\n\t"\
		"subpd	0x010(%%rdx),%%xmm5		\n\t"\
		"addpd	     (%%rcx),%%xmm6		\n\t"\
		"addpd	0x010(%%rcx),%%xmm7		\n\t"\
		"\n\t"\
		"\n\t"/* SSE2_RADIX4_DIT_LASTPAIR_IN_PLACE_STRIDE(xmm2,xmm3,xmm0,xmm1,xmm6,xmm7,xmm4,xmm5): swap xmm[01]<->[23], xmm[45]<->[67] */\
		"addpd	%%xmm6,%%xmm2		\n\t"\
		"addpd	%%xmm7,%%xmm3		\n\t"\
		"movaps	%%xmm2,     (%%rax)	\n\t"\
		"movaps	%%xmm3,0x010(%%rax)	\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"subpd	%%xmm6,%%xmm2		\n\t"\
		"subpd	%%xmm7,%%xmm3		\n\t"\
		"movaps	%%xmm2,     (%%rcx)	\n\t"\
		"movaps	%%xmm3,0x010(%%rcx)	\n\t"\
		"\n\t"\
		"addpd	%%xmm5,%%xmm0		\n\t"\
		"subpd	%%xmm4,%%xmm1		\n\t"\
		"movaps	%%xmm0,     (%%rbx)	\n\t"\
		"movaps	%%xmm1,0x010(%%rbx)	\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"subpd	%%xmm5,%%xmm0		\n\t"\
		"addpd	%%xmm4,%%xmm1		\n\t"\
		"movaps	%%xmm0,     (%%rdx)	\n\t"\
		"movaps	%%xmm1,0x010(%%rdx)	\n\t"\
		"\n\t"\
		"\n\t"/*...Block 5: r08,r18,r28,r38	*/\
		"\n\t"\
		"addq	$0x080,%%rax		\n\t"/* r08 */\
		"addq	$0x080,%%rbx		\n\t"\
		"addq	$0x080,%%rcx		\n\t"\
		"addq	$0x080,%%rdx		\n\t"\
		"movaps	(%%rdi),%%xmm2	\n\t"/* isrt2 */\
		"\n\t"\
		"movaps	     (%%rcx),%%xmm4	\n\t"\
		"movaps	0x010(%%rcx),%%xmm5	\n\t"\
		"movaps	     (%%rdx),%%xmm0	\n\t"\
		"movaps	0x010(%%rdx),%%xmm1	\n\t"\
		"addpd	0x010(%%rcx),%%xmm4	\n\t"\
		"subpd	     (%%rcx),%%xmm5	\n\t"\
		"subpd	0x010(%%rdx),%%xmm0	\n\t"\
		"addpd	     (%%rdx),%%xmm1	\n\t"\
		"mulpd	%%xmm2,%%xmm4		\n\t"\
		"mulpd	%%xmm2,%%xmm5		\n\t"\
		"mulpd	%%xmm2,%%xmm0		\n\t"\
		"mulpd	%%xmm2,%%xmm1		\n\t"\
		"movaps	%%xmm4,%%xmm6		\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t"\
		"\n\t"\
		"subpd	%%xmm0,%%xmm4		\n\t"\
		"subpd	%%xmm1,%%xmm5		\n\t"\
		"addpd	%%xmm0,%%xmm6		\n\t"\
		"addpd	%%xmm1,%%xmm7		\n\t"\
		"\n\t"\
		"movaps	     (%%rax),%%xmm0	\n\t"\
		"movaps	0x010(%%rax),%%xmm1	\n\t"\
		"movaps	     (%%rbx),%%xmm2	\n\t"\
		"movaps	0x010(%%rbx),%%xmm3	\n\t"\
		"subpd	0x010(%%rbx),%%xmm0	\n\t"\
		"subpd	     (%%rbx),%%xmm1	\n\t"\
		"addpd	     (%%rax),%%xmm3	\n\t"\
		"addpd	0x010(%%rax),%%xmm2	\n\t"\
		"\n\t"\
		"\n\t"/* SSE2_RADIX4_DIT_LASTPAIR_IN_PLACE_STRIDE(xmm3,xmm1,xmm0,xmm2,xmm4,xmm5,xmm6,xmm7): swap xmm0123<->3102 */\
		"addpd	%%xmm4,%%xmm3		\n\t"\
		"addpd	%%xmm5,%%xmm1		\n\t"\
		"movaps	%%xmm3,     (%%rax)	\n\t"\
		"movaps	%%xmm1,0x010(%%rax)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"subpd	%%xmm4,%%xmm3		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"movaps	%%xmm3,     (%%rcx)	\n\t"\
		"movaps	%%xmm1,0x010(%%rcx)	\n\t"\
		"\n\t"\
		"addpd	%%xmm7,%%xmm0		\n\t"\
		"subpd	%%xmm6,%%xmm2		\n\t"\
		"movaps	%%xmm0,     (%%rbx)	\n\t"\
		"movaps	%%xmm2,0x010(%%rbx)	\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"subpd	%%xmm7,%%xmm0		\n\t"\
		"addpd	%%xmm6,%%xmm2		\n\t"\
		"movaps	%%xmm0,     (%%rdx)	\n\t"\
		"movaps	%%xmm2,0x010(%%rdx)	\n\t"\
		"\n\t"\
		"\n\t"/*...Block 3: r04,r14,r24,r34	*/\
		"\n\t"\
		"subq	$0x040,%%rax 	\n\t"/* r04 */\
		"subq	$0x040,%%rbx		\n\t"\
		"subq	$0x040,%%rcx		\n\t"\
		"subq	$0x040,%%rdx		\n\t"\
		"\n\t"\
		"movaps	     (%%rcx),%%xmm4	\n\t"\
		"movaps	     (%%rdx),%%xmm0	\n\t"\
		"movaps	0x010(%%rcx),%%xmm5	\n\t"\
		"movaps	0x010(%%rdx),%%xmm1	\n\t"\
		"movaps	     (%%rcx),%%xmm6	\n\t"\
		"movaps	     (%%rdx),%%xmm2	\n\t"\
		"movaps	0x010(%%rcx),%%xmm7	\n\t"\
		"movaps	0x010(%%rdx),%%xmm3	\n\t"\
		"\n\t"\
		"mulpd	     (%%rsi),%%xmm4	\n\t"\
		"mulpd	0x010(%%rsi),%%xmm0	\n\t"\
		"mulpd	     (%%rsi),%%xmm5	\n\t"\
		"mulpd	0x010(%%rsi),%%xmm1	\n\t"\
		"mulpd	0x010(%%rsi),%%xmm6	\n\t"\
		"mulpd	     (%%rsi),%%xmm2	\n\t"\
		"mulpd	0x010(%%rsi),%%xmm7	\n\t"\
		"mulpd	     (%%rsi),%%xmm3	\n\t"\
		"subpd	%%xmm6,%%xmm5		\n\t"\
		"subpd	%%xmm2,%%xmm1		\n\t"\
		"addpd	%%xmm7,%%xmm4		\n\t"\
		"addpd	%%xmm3,%%xmm0		\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t"\
		"movaps	%%xmm4,%%xmm6		\n\t"\
		"\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"subpd	%%xmm0,%%xmm6		\n\t"\
		"subpd	%%xmm1,%%xmm7		\n\t"\
		"\n\t"\
		"movaps	     (%%rbx),%%xmm2	\n\t"\
		"movaps	0x010(%%rbx),%%xmm3	\n\t"\
		"movaps	     (%%rax),%%xmm0	\n\t"\
		"movaps	0x010(%%rax),%%xmm1	\n\t"\
		"addpd	0x010(%%rbx),%%xmm2	\n\t"\
		"subpd	     (%%rbx),%%xmm3	\n\t"\
		"mulpd	     (%%rdi),%%xmm2		\n\t"\
		"mulpd	     (%%rdi),%%xmm3		\n\t"\
		"\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"\n\t"\
		"\n\t"/* SSE2_RADIX4_DIT_LASTPAIR_IN_PLACE_STRIDE(xmm2,xmm3,xmm0,xmm1,xmm4,xmm5,xmm6,xmm7): swap xmm[01]<->[23] */\
		"addpd	%%xmm4,%%xmm2		\n\t"\
		"addpd	%%xmm5,%%xmm3		\n\t"\
		"movaps	%%xmm2,     (%%rax)	\n\t"\
		"movaps	%%xmm3,0x010(%%rax)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"subpd	%%xmm4,%%xmm2		\n\t"\
		"subpd	%%xmm5,%%xmm3		\n\t"\
		"movaps	%%xmm2,     (%%rcx)	\n\t"\
		"movaps	%%xmm3,0x010(%%rcx)	\n\t"\
		"\n\t"\
		"addpd	%%xmm7,%%xmm0		\n\t"\
		"subpd	%%xmm6,%%xmm1		\n\t"\
		"movaps	%%xmm0,     (%%rbx)	\n\t"\
		"movaps	%%xmm1,0x010(%%rbx)	\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"subpd	%%xmm7,%%xmm0		\n\t"\
		"addpd	%%xmm6,%%xmm1		\n\t"\
		"movaps	%%xmm0,     (%%rdx)	\n\t"\
		"movaps	%%xmm1,0x010(%%rdx)	\n\t"\
		"\n\t"\
		"\n\t"/*...Block 7: r0C,r1C,r2C,r3C	*/\
		"\n\t"\
		"addq	$0x080,%%rax 	\n\t"/* r0C */\
		"addq	$0x080,%%rbx		\n\t"\
		"addq	$0x080,%%rcx		\n\t"\
		"addq	$0x080,%%rdx		\n\t"\
		"\n\t"\
		"movaps	     (%%rcx),%%xmm4	\n\t"\
		"movaps	     (%%rdx),%%xmm0	\n\t"\
		"movaps	0x010(%%rcx),%%xmm5	\n\t"\
		"movaps	0x010(%%rdx),%%xmm1	\n\t"\
		"movaps	     (%%rcx),%%xmm6	\n\t"\
		"movaps	     (%%rdx),%%xmm2	\n\t"\
		"movaps	0x010(%%rcx),%%xmm7	\n\t"\
		"movaps	0x010(%%rdx),%%xmm3	\n\t"\
		"\n\t"\
		"mulpd	0x010(%%rsi),%%xmm4	\n\t"\
		"mulpd	     (%%rsi),%%xmm0	\n\t"\
		"mulpd	0x010(%%rsi),%%xmm5	\n\t"\
		"mulpd	     (%%rsi),%%xmm1	\n\t"\
		"mulpd	     (%%rsi),%%xmm6	\n\t"\
		"mulpd	0x010(%%rsi),%%xmm2	\n\t"\
		"mulpd	     (%%rsi),%%xmm7	\n\t"\
		"mulpd	0x010(%%rsi),%%xmm3	\n\t"\
		"subpd	%%xmm6,%%xmm5		\n\t"\
		"subpd	%%xmm2,%%xmm1		\n\t"\
		"addpd	%%xmm7,%%xmm4		\n\t"\
		"addpd	%%xmm3,%%xmm0		\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t"\
		"movaps	%%xmm4,%%xmm6		\n\t"\
		"\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"subpd	%%xmm0,%%xmm6		\n\t"\
		"subpd	%%xmm1,%%xmm7		\n\t"\
		"\n\t"\
		"movaps	     (%%rbx),%%xmm2	\n\t"\
		"movaps	0x010(%%rbx),%%xmm3	\n\t"\
		"movaps	     (%%rax),%%xmm0	\n\t"\
		"movaps	0x010(%%rax),%%xmm1	\n\t"\
		"subpd	0x010(%%rbx),%%xmm2	\n\t"\
		"addpd	     (%%rbx),%%xmm3	\n\t"\
		"mulpd	     (%%rdi),%%xmm2		\n\t"\
		"mulpd	     (%%rdi),%%xmm3		\n\t"\
		"\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"\n\t"\
		"\n\t"/* SSE2_RADIX4_DIT_LASTPAIR_IN_PLACE_STRIDE(xmm0,xmm1,xmm2,xmm3,xmm6,xmm7,xmm4,xmm5): swap xmm[45]<->[67] */\
		"addpd	%%xmm6,%%xmm0		\n\t"\
		"addpd	%%xmm7,%%xmm1		\n\t"\
		"movaps	%%xmm0,     (%%rax)	\n\t"\
		"movaps	%%xmm1,0x010(%%rax)	\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"subpd	%%xmm6,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm1		\n\t"\
		"movaps	%%xmm0,     (%%rcx)	\n\t"\
		"movaps	%%xmm1,0x010(%%rcx)	\n\t"\
		"\n\t"\
		"addpd	%%xmm5,%%xmm2		\n\t"\
		"subpd	%%xmm4,%%xmm3		\n\t"\
		"movaps	%%xmm2,     (%%rbx)	\n\t"\
		"movaps	%%xmm3,0x010(%%rbx)	\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"subpd	%%xmm5,%%xmm2		\n\t"\
		"addpd	%%xmm4,%%xmm3		\n\t"\
		"movaps	%%xmm2,     (%%rdx)	\n\t"\
		"movaps	%%xmm3,0x010(%%rdx)	\n\t"\
		"\n\t"\
		"\n\t"/*...Block 2: r02,r12,r22,r32	*/\
		"\n\t"\
		"subq	$0x0a0,%%rax 	\n\t"/* r02 */\
		"subq	$0x0a0,%%rbx			\n\t"\
		"subq	$0x0a0,%%rcx			\n\t"\
		"subq	$0x0a0,%%rdx			\n\t"\
		"addq	$0x030,%%rdi \n\t"/* cc1 */\
		"addq	$0x040,%%rsi \n\t"/* cc3 */\
		"\n\t"\
		"movaps	     (%%rcx),%%xmm4	\n\t"\
		"movaps	     (%%rdx),%%xmm0	\n\t"\
		"movaps	0x010(%%rcx),%%xmm5	\n\t"\
		"movaps	0x010(%%rdx),%%xmm1	\n\t"\
		"movaps	     (%%rcx),%%xmm6	\n\t"\
		"movaps	     (%%rdx),%%xmm2	\n\t"\
		"movaps	0x010(%%rcx),%%xmm7	\n\t"\
		"movaps	0x010(%%rdx),%%xmm3	\n\t"\
		"\n\t"\
		"mulpd	     (%%rdi),%%xmm4	\n\t"\
		"mulpd	     (%%rsi),%%xmm0	\n\t"\
		"mulpd	     (%%rdi),%%xmm5	\n\t"\
		"mulpd	     (%%rsi),%%xmm1	\n\t"\
		"mulpd	0x010(%%rdi),%%xmm6	\n\t"\
		"mulpd	0x010(%%rsi),%%xmm2	\n\t"\
		"mulpd	0x010(%%rdi),%%xmm7	\n\t"\
		"mulpd	0x010(%%rsi),%%xmm3	\n\t"\
		"subpd	%%xmm6,%%xmm5		\n\t"\
		"subpd	%%xmm2,%%xmm1		\n\t"\
		"addpd	%%xmm7,%%xmm4		\n\t"\
		"addpd	%%xmm3,%%xmm0		\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t"\
		"movaps	%%xmm4,%%xmm6		\n\t"\
		"\n\t"\
		"subpd	%%xmm0,%%xmm4		\n\t"\
		"subpd	%%xmm1,%%xmm5		\n\t"\
		"addpd	%%xmm0,%%xmm6		\n\t"\
		"addpd	%%xmm1,%%xmm7		\n\t"\
		"\n\t"\
		"subq	$0x040,%%rsi \n\t"/* cc0 */\
		"movaps	     (%%rbx),%%xmm2	\n\t"\
		"movaps	0x010(%%rbx),%%xmm3	\n\t"\
		"movaps	     (%%rbx),%%xmm0	\n\t"\
		"movaps	0x010(%%rbx),%%xmm1	\n\t"\
		"mulpd	     (%%rsi),%%xmm2	\n\t"\
		"mulpd	0x010(%%rsi),%%xmm1	\n\t"\
		"mulpd	     (%%rsi),%%xmm3	\n\t"\
		"mulpd	0x010(%%rsi),%%xmm0	\n\t"\
		"addpd	%%xmm1,%%xmm2		\n\t"\
		"subpd	%%xmm0,%%xmm3		\n\t"\
		"\n\t"\
		"movaps	     (%%rax),%%xmm0	\n\t"\
		"movaps	0x010(%%rax),%%xmm1	\n\t"\
		"\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"\n\t"\
		"\n\t"/* SSE2_RADIX4_DIT_LASTPAIR_IN_PLACE_STRIDE(xmm2,xmm3,xmm0,xmm1,xmm6,xmm7,xmm4,xmm5): swap xmm[01]<->[23], xmm[45]<->[67] */\
		"addpd	%%xmm6,%%xmm2		\n\t"\
		"addpd	%%xmm7,%%xmm3		\n\t"\
		"movaps	%%xmm2,     (%%rax)	\n\t"\
		"movaps	%%xmm3,0x010(%%rax)	\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"subpd	%%xmm6,%%xmm2		\n\t"\
		"subpd	%%xmm7,%%xmm3		\n\t"\
		"movaps	%%xmm2,     (%%rcx)	\n\t"\
		"movaps	%%xmm3,0x010(%%rcx)	\n\t"\
		"\n\t"\
		"addpd	%%xmm5,%%xmm0		\n\t"\
		"subpd	%%xmm4,%%xmm1		\n\t"\
		"movaps	%%xmm0,     (%%rbx)	\n\t"\
		"movaps	%%xmm1,0x010(%%rbx)	\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"subpd	%%xmm5,%%xmm0		\n\t"\
		"addpd	%%xmm4,%%xmm1		\n\t"\
		"movaps	%%xmm0,     (%%rdx)	\n\t"\
		"movaps	%%xmm1,0x010(%%rdx)	\n\t"\
		"\n\t"\
		"\n\t"/*...Block 6: r0A,r1A,r2A,r3A	*/\
		"\n\t"\
		"addq	$0x080,%%rax 	\n\t"/* r0A */\
		"addq	$0x080,%%rbx			\n\t"\
		"addq	$0x080,%%rcx			\n\t"\
		"addq	$0x080,%%rdx			\n\t"\
		"addq	$0x040,%%rsi \n\t"/* cc3 */\
		"\n\t"\
		"movaps	     (%%rcx),%%xmm4	\n\t"\
		"movaps	     (%%rdx),%%xmm0	\n\t"\
		"movaps	0x010(%%rcx),%%xmm5	\n\t"\
		"movaps	0x010(%%rdx),%%xmm1	\n\t"\
		"movaps	     (%%rcx),%%xmm6	\n\t"\
		"movaps	     (%%rdx),%%xmm2	\n\t"\
		"movaps	0x010(%%rcx),%%xmm7	\n\t"\
		"movaps	0x010(%%rdx),%%xmm3	\n\t"\
		"\n\t"\
		"mulpd	0x010(%%rsi),%%xmm4	\n\t"\
		"mulpd	     (%%rdi),%%xmm0	\n\t"\
		"mulpd	0x010(%%rsi),%%xmm5	\n\t"\
		"mulpd	     (%%rdi),%%xmm1	\n\t"\
		"mulpd	     (%%rsi),%%xmm6	\n\t"\
		"mulpd	0x010(%%rdi),%%xmm2	\n\t"\
		"mulpd	     (%%rsi),%%xmm7	\n\t"\
		"mulpd	0x010(%%rdi),%%xmm3	\n\t"\
		"subpd	%%xmm6,%%xmm5		\n\t"\
		"addpd	%%xmm2,%%xmm1		\n\t"\
		"addpd	%%xmm7,%%xmm4		\n\t"\
		"subpd	%%xmm3,%%xmm0		\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t"\
		"movaps	%%xmm4,%%xmm6		\n\t"\
		"\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"subpd	%%xmm0,%%xmm6		\n\t"\
		"subpd	%%xmm1,%%xmm7		\n\t"\
		"\n\t"\
		"subq	$0x040,%%rsi \n\t"/* cc0 */\
		"movaps	     (%%rbx),%%xmm2	\n\t"\
		"movaps	0x010(%%rbx),%%xmm3	\n\t"\
		"movaps	     (%%rbx),%%xmm0	\n\t"\
		"movaps	0x010(%%rbx),%%xmm1	\n\t"\
		"mulpd	0x010(%%rsi),%%xmm2	\n\t"\
		"mulpd	     (%%rsi),%%xmm1	\n\t"\
		"mulpd	0x010(%%rsi),%%xmm3	\n\t"\
		"mulpd	     (%%rsi),%%xmm0	\n\t"\
		"subpd	%%xmm1,%%xmm2		\n\t"\
		"addpd	%%xmm0,%%xmm3		\n\t"\
		"\n\t"\
		"movaps	     (%%rax),%%xmm0	\n\t"\
		"movaps	0x010(%%rax),%%xmm1	\n\t"\
		"\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"\n\t"\
		"\n\t"/* SSE2_RADIX4_DIT_LASTPAIR_IN_PLACE_STRIDE(xmm0,xmm1,xmm2,xmm3,xmm6,xmm7,xmm4,xmm5): swap xmm[45]<->[67] */\
		"addpd	%%xmm6,%%xmm0		\n\t"\
		"addpd	%%xmm7,%%xmm1		\n\t"\
		"movaps	%%xmm0,     (%%rax)	\n\t"\
		"movaps	%%xmm1,0x010(%%rax)	\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"subpd	%%xmm6,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm1		\n\t"\
		"movaps	%%xmm0,     (%%rcx)	\n\t"\
		"movaps	%%xmm1,0x010(%%rcx)	\n\t"\
		"\n\t"\
		"addpd	%%xmm5,%%xmm2		\n\t"\
		"subpd	%%xmm4,%%xmm3		\n\t"\
		"movaps	%%xmm2,     (%%rbx)	\n\t"\
		"movaps	%%xmm3,0x010(%%rbx)	\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"subpd	%%xmm5,%%xmm2		\n\t"\
		"addpd	%%xmm4,%%xmm3		\n\t"\
		"movaps	%%xmm2,     (%%rdx)	\n\t"\
		"movaps	%%xmm3,0x010(%%rdx)	\n\t"\
		"\n\t"\
		"\n\t"/*...Block 4: r06,r16,r26,r36	*/\
		"\n\t"\
		"subq	$0x040,%%rax 	\n\t"/* r06 */\
		"subq	$0x040,%%rbx			\n\t"\
		"subq	$0x040,%%rcx			\n\t"\
		"subq	$0x040,%%rdx			\n\t"\
		"addq	$0x040,%%rsi \n\t"/* cc3 */\
		"\n\t"\
		"movaps	     (%%rcx),%%xmm4	\n\t"\
		"movaps	     (%%rdx),%%xmm0	\n\t"\
		"movaps	0x010(%%rcx),%%xmm5	\n\t"\
		"movaps	0x010(%%rdx),%%xmm1	\n\t"\
		"movaps	     (%%rcx),%%xmm6	\n\t"\
		"movaps	     (%%rdx),%%xmm2	\n\t"\
		"movaps	0x010(%%rcx),%%xmm7	\n\t"\
		"movaps	0x010(%%rdx),%%xmm3	\n\t"\
		"\n\t"\
		"mulpd	     (%%rsi),%%xmm4	\n\t"\
		"mulpd	0x010(%%rdi),%%xmm0	\n\t"\
		"mulpd	     (%%rsi),%%xmm5	\n\t"\
		"mulpd	0x010(%%rdi),%%xmm1	\n\t"\
		"mulpd	0x010(%%rsi),%%xmm6	\n\t"\
		"mulpd	     (%%rdi),%%xmm2	\n\t"\
		"mulpd	0x010(%%rsi),%%xmm7	\n\t"\
		"mulpd	     (%%rdi),%%xmm3	\n\t"\
		"subpd	%%xmm6,%%xmm5		\n\t"\
		"addpd	%%xmm2,%%xmm1		\n\t"\
		"addpd	%%xmm7,%%xmm4		\n\t"\
		"subpd	%%xmm3,%%xmm0		\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t"\
		"movaps	%%xmm4,%%xmm6		\n\t"\
		"\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"subpd	%%xmm0,%%xmm6		\n\t"\
		"subpd	%%xmm1,%%xmm7		\n\t"\
		"\n\t"\
		"subq	$0x040,%%rsi \n\t"/* cc0 */\
		"movaps	     (%%rbx),%%xmm2	\n\t"\
		"movaps	0x010(%%rbx),%%xmm3	\n\t"\
		"movaps	     (%%rbx),%%xmm0	\n\t"\
		"movaps	0x010(%%rbx),%%xmm1	\n\t"\
		"mulpd	0x010(%%rsi),%%xmm2	\n\t"\
		"mulpd	     (%%rsi),%%xmm1	\n\t"\
		"mulpd	0x010(%%rsi),%%xmm3	\n\t"\
		"mulpd	     (%%rsi),%%xmm0	\n\t"\
		"addpd	%%xmm1,%%xmm2		\n\t"\
		"subpd	%%xmm0,%%xmm3		\n\t"\
		"\n\t"\
		"movaps	     (%%rax),%%xmm0	\n\t"\
		"movaps	0x010(%%rax),%%xmm1	\n\t"\
		"\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"\n\t"\
		"\n\t"/* SSE2_RADIX4_DIT_LASTPAIR_IN_PLACE_STRIDE(xmm2,xmm3,xmm0,xmm1,xmm6,xmm7,xmm4,xmm5): swap xmm[01]<->[23], xmm[45]<->[67] */\
		"addpd	%%xmm6,%%xmm2		\n\t"\
		"addpd	%%xmm7,%%xmm3		\n\t"\
		"movaps	%%xmm2,     (%%rax)	\n\t"\
		"movaps	%%xmm3,0x010(%%rax)	\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"subpd	%%xmm6,%%xmm2		\n\t"\
		"subpd	%%xmm7,%%xmm3		\n\t"\
		"movaps	%%xmm2,     (%%rcx)	\n\t"\
		"movaps	%%xmm3,0x010(%%rcx)	\n\t"\
		"\n\t"\
		"addpd	%%xmm5,%%xmm0		\n\t"\
		"subpd	%%xmm4,%%xmm1		\n\t"\
		"movaps	%%xmm0,     (%%rbx)	\n\t"\
		"movaps	%%xmm1,0x010(%%rbx)	\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"subpd	%%xmm5,%%xmm0		\n\t"\
		"addpd	%%xmm4,%%xmm1		\n\t"\
		"movaps	%%xmm0,     (%%rdx)	\n\t"\
		"movaps	%%xmm1,0x010(%%rdx)	\n\t"\
		"\n\t"\
		"\n\t"/*...Block 8: r0E,r1E,r2E,r3E	*/\
		"\n\t"\
		"addq	$0x080,%%rax 	\n\t"/* r0E */\
		"addq	$0x080,%%rbx		\n\t"\
		"addq	$0x080,%%rcx		\n\t"\
		"addq	$0x080,%%rdx		\n\t"\
		"addq	$0x040,%%rsi \n\t"/* cc3 */\
		"\n\t"\
		"movaps	     (%%rcx),%%xmm4	\n\t"\
		"movaps	     (%%rdx),%%xmm0	\n\t"\
		"movaps	0x010(%%rcx),%%xmm5	\n\t"\
		"movaps	0x010(%%rdx),%%xmm1	\n\t"\
		"movaps	     (%%rcx),%%xmm6	\n\t"\
		"movaps	     (%%rdx),%%xmm2	\n\t"\
		"movaps	0x010(%%rcx),%%xmm7	\n\t"\
		"movaps	0x010(%%rdx),%%xmm3	\n\t"\
		"\n\t"\
		"mulpd	0x010(%%rdi),%%xmm4	\n\t"\
		"mulpd	0x010(%%rsi),%%xmm0	\n\t"\
		"mulpd	0x010(%%rdi),%%xmm5	\n\t"\
		"mulpd	0x010(%%rsi),%%xmm1	\n\t"\
		"mulpd	     (%%rdi),%%xmm6	\n\t"\
		"mulpd	     (%%rsi),%%xmm2	\n\t"\
		"mulpd	     (%%rdi),%%xmm7	\n\t"\
		"mulpd	     (%%rsi),%%xmm3	\n\t"\
		"subpd	%%xmm6,%%xmm5		\n\t"\
		"subpd	%%xmm2,%%xmm1		\n\t"\
		"addpd	%%xmm7,%%xmm4		\n\t"\
		"addpd	%%xmm3,%%xmm0		\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t"\
		"movaps	%%xmm4,%%xmm6		\n\t"\
		"\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"subpd	%%xmm0,%%xmm6		\n\t"\
		"subpd	%%xmm1,%%xmm7		\n\t"\
		"\n\t"\
		"subq	$0x040,%%rsi \n\t"/* cc0 */\
		"movaps	     (%%rbx),%%xmm2	\n\t"\
		"movaps	0x010(%%rbx),%%xmm3	\n\t"\
		"movaps	     (%%rbx),%%xmm0	\n\t"\
		"movaps	0x010(%%rbx),%%xmm1	\n\t"\
		"mulpd	     (%%rsi),%%xmm2	\n\t"\
		"mulpd	0x010(%%rsi),%%xmm1	\n\t"\
		"mulpd	     (%%rsi),%%xmm3	\n\t"\
		"mulpd	0x010(%%rsi),%%xmm0	\n\t"\
		"subpd	%%xmm1,%%xmm2		\n\t"\
		"addpd	%%xmm0,%%xmm3		\n\t"\
		"\n\t"\
		"movaps	     (%%rax),%%xmm0	\n\t"\
		"movaps	0x010(%%rax),%%xmm1	\n\t"\
		"\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"\n\t"\
		"\n\t"/* SSE2_RADIX4_DIT_LASTPAIR_IN_PLACE_STRIDE(xmm0,xmm1,xmm2,xmm3,xmm6,xmm7,xmm4,xmm5): swap xmm[45]<->[67] */\
		"addpd	%%xmm6,%%xmm0		\n\t"\
		"addpd	%%xmm7,%%xmm1		\n\t"\
		"movaps	%%xmm0,     (%%rax)	\n\t"\
		"movaps	%%xmm1,0x010(%%rax)	\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"subpd	%%xmm6,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm1		\n\t"\
		"movaps	%%xmm0,     (%%rcx)	\n\t"\
		"movaps	%%xmm1,0x010(%%rcx)	\n\t"\
		"\n\t"\
		"addpd	%%xmm5,%%xmm2		\n\t"\
		"subpd	%%xmm4,%%xmm3		\n\t"\
		"movaps	%%xmm2,     (%%rbx)	\n\t"\
		"movaps	%%xmm3,0x010(%%rbx)	\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"subpd	%%xmm5,%%xmm2		\n\t"\
		"addpd	%%xmm4,%%xmm3		\n\t"\
		"movaps	%%xmm2,     (%%rdx)	\n\t"\
		"movaps	%%xmm3,0x010(%%rdx)	\n\t"\
		:					/* outputs: none */\
		: [__add] "m" (Xadd)	/* All inputs from memory addresses here */\
		 ,[__p01] "m" (Xp01)\
		 ,[__p02] "m" (Xp02)\
		 ,[__p03] "m" (Xp03)\
		 ,[__p04] "m" (Xp04)\
		 ,[__p08] "m" (Xp08)\
		 ,[__p10] "m" (Xp10)\
		 ,[__p18] "m" (Xp18)\
		 ,[__r00] "m" (Xr00)\
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__cc0] "m" (Xcc0)\
		: "rax","rbx","rcx","rdx","rdi","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","memory"		/* Clobbered registers */\
	);\
	}


	#define SSE2_RADIX32_DIF_NOTWIDDLE(Xadd,Xp01,Xp02,Xp03,Xp04,Xp08,Xp10,Xp18,Xr00,Xisrt2,Xcc0)\
	{\
	__asm__ volatile (\
		"\n\t"/* SSE2_RADIX4_DIF_IN_PLACE(r00,r20,r10,r30): */\
		"\n\t"\
		"movq	%[__isrt2],%%rdi	\n\t"\
		"movq	%[__r00],%%rax	\n\t"\
		"movq	$0x200,%%rbx	\n\t"\
		"movq	$0x100,%%rcx	\n\t"\
		"movq	$0x300,%%rdx	\n\t"\
		"addq	%%rax,%%rbx		\n\t"\
		"addq	%%rax,%%rcx		\n\t"\
		"addq	%%rax,%%rdx		\n\t"\
		"\n\t"\
		"movaps	     (%%rax),%%xmm0	\n\t"\
		"movaps	0x010(%%rax),%%xmm1	\n\t"\
		"movaps	     (%%rax),%%xmm2	\n\t"\
		"movaps	0x010(%%rax),%%xmm3	\n\t"\
		"addpd	     (%%rbx),%%xmm0	\n\t"\
		"addpd	0x010(%%rbx),%%xmm1	\n\t"\
		"subpd	     (%%rbx),%%xmm2	\n\t"\
		"subpd	0x010(%%rbx),%%xmm3	\n\t"\
		"\n\t"\
		"movaps	     (%%rcx),%%xmm4	\n\t"\
		"movaps	0x010(%%rcx),%%xmm5	\n\t"\
		"movaps	     (%%rcx),%%xmm6	\n\t"\
		"movaps	0x010(%%rcx),%%xmm7	\n\t"\
		"addpd	     (%%rdx),%%xmm4	\n\t"\
		"addpd	0x010(%%rdx),%%xmm5	\n\t"\
		"subpd	     (%%rdx),%%xmm6	\n\t"\
		"subpd	0x010(%%rdx),%%xmm7	\n\t"\
		"\n\t"/* Finish radix-4 butterfly and store results into temporary-array slots: */\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"movaps	%%xmm0,     (%%rbx)	\n\t"\
		"movaps	%%xmm1,0x010(%%rbx)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"movaps	%%xmm4,     (%%rax)	\n\t"\
		"movaps	%%xmm5,0x010(%%rax)	\n\t"\
		"\n\t"\
		"subpd	%%xmm7,%%xmm2		\n\t"\
		"subpd	%%xmm6,%%xmm3		\n\t"\
		"movaps	%%xmm2,     (%%rcx)	\n\t"\
		"movaps	%%xmm3,0x010(%%rdx)	\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	%%xmm7,     (%%rdx)	\n\t"\
		"movaps	%%xmm6,0x010(%%rcx)	\n\t"\
		"\n\t"\
		"\n\t"/* SSE2_RADIX4_DIF_IN_PLACE_2NDOFTWO(r08,r28,r18,r38): */\
		"addq	$0x080,%%rax		\n\t"\
		"addq	$0x080,%%rbx		\n\t"\
		"addq	$0x080,%%rcx		\n\t"\
		"addq	$0x080,%%rdx		\n\t"\
		"movaps	     (%%rax),%%xmm0	\n\t"\
		"movaps	0x010(%%rax),%%xmm1	\n\t"\
		"movaps	     (%%rax),%%xmm2	\n\t"\
		"movaps	0x010(%%rax),%%xmm3	\n\t"\
		"addpd	     (%%rbx),%%xmm0	\n\t"\
		"addpd	0x010(%%rbx),%%xmm1	\n\t"\
		"subpd	     (%%rbx),%%xmm2	\n\t"\
		"subpd	0x010(%%rbx),%%xmm3	\n\t"\
		"movaps	     (%%rcx),%%xmm4	\n\t"\
		"movaps	0x010(%%rcx),%%xmm5	\n\t"\
		"movaps	     (%%rcx),%%xmm6	\n\t"\
		"movaps	0x010(%%rcx),%%xmm7	\n\t"\
		"addpd	     (%%rdx),%%xmm6	\n\t"\
		"addpd	0x010(%%rdx),%%xmm7	\n\t"\
		"subpd	     (%%rdx),%%xmm4	\n\t"\
		"subpd	0x010(%%rdx),%%xmm5	\n\t"\
		"\n\t"/* Finish radix-4 butterfly and store results into temporary-array slots: */\
		"subpd	%%xmm6,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm1		\n\t"\
		"subpd	%%xmm5,%%xmm2		\n\t"\
		"movaps	%%xmm0,     (%%rbx)	\n\t"\
		"movaps	%%xmm1,0x010(%%rbx)	\n\t"\
		"subpd	%%xmm4,%%xmm3		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm0,%%xmm6		\n\t"\
		"addpd	%%xmm2,%%xmm5		\n\t"\
		"addpd	%%xmm1,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm4		\n\t"\
		"movaps	%%xmm6,     (%%rax)	\n\t"\
		"movaps	%%xmm7,0x010(%%rax)	\n\t"\
		"\n\t"\
		"movaps	(%%rdi),%%xmm0	\n\t"/* isrt2 */\
		"movaps	%%xmm2,%%xmm6		\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t"\
		"subpd	%%xmm4,%%xmm2		\n\t"\
		"subpd	%%xmm3,%%xmm5		\n\t"\
		"addpd	%%xmm4,%%xmm6		\n\t"\
		"addpd	%%xmm3,%%xmm7		\n\t"\
		"mulpd	%%xmm0,%%xmm2		\n\t"\
		"mulpd	%%xmm0,%%xmm5		\n\t"\
		"mulpd	%%xmm0,%%xmm6		\n\t"\
		"mulpd	%%xmm0,%%xmm7		\n\t"\
		"movaps	%%xmm2,     (%%rcx)	\n\t"\
		"movaps	%%xmm5,     (%%rdx)	\n\t"\
		"movaps	%%xmm6,0x010(%%rcx)	\n\t"\
		"movaps	%%xmm7,0x010(%%rdx)	\n\t"\
		"\n\t"/* SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS(r00,r10,r20,r30,r08,r18,r28,r38): */\
		"subq	$0x080,%%rax	\n\t"/* r00 */\
		"subq	$0x200,%%rbx	\n\t"/* r08 */\
		"addq	$0x080,%%rcx	\n\t"/* r20 */\
		"subq	$0x100,%%rdx	\n\t"/* r28 */\
		"movaps	     (%%rax),%%xmm0	\n\t"\
		"movaps	     (%%rcx),%%xmm4	\n\t"\
		"movaps	0x010(%%rax),%%xmm1	\n\t"\
		"movaps	0x010(%%rcx),%%xmm5	\n\t"\
		"movaps	     (%%rbx),%%xmm2	\n\t"\
		"movaps	0x010(%%rdx),%%xmm7	\n\t"\
		"movaps	0x010(%%rbx),%%xmm3	\n\t"\
		"movaps	     (%%rdx),%%xmm6	\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm4		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"subpd	%%xmm6,%%xmm5		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm4,%%xmm7		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"addpd	%%xmm5,%%xmm6		\n\t"\
		"\n\t"\
		"movaps	%%xmm0,     (%%rbx)	\n\t"\
		"movaps	%%xmm4,     (%%rcx)	\n\t"\
		"movaps	%%xmm1,0x010(%%rbx)	\n\t"\
		"movaps	%%xmm5,0x010(%%rdx)	\n\t"\
		"movaps	%%xmm2,     (%%rax)	\n\t"\
		"movaps	%%xmm7,     (%%rdx)	\n\t"\
		"movaps	%%xmm3,0x010(%%rax)	\n\t"\
		"movaps	%%xmm6,0x010(%%rcx)	\n\t"\
		"\n\t"\
		"addq	$0x100,%%rax	\n\t"/* r10 */\
		"addq	$0x100,%%rbx	\n\t"/* r18 */\
		"addq	$0x100,%%rcx	\n\t"/* r30 */\
		"addq	$0x100,%%rdx	\n\t"/* r38 */\
		"\n\t"\
		"movaps	     (%%rax),%%xmm0	\n\t"\
		"movaps	     (%%rcx),%%xmm4	\n\t"\
		"movaps	0x010(%%rax),%%xmm1	\n\t"\
		"movaps	0x010(%%rcx),%%xmm5	\n\t"\
		"movaps	     (%%rbx),%%xmm2	\n\t"\
		"movaps	0x010(%%rdx),%%xmm7	\n\t"\
		"movaps	0x010(%%rbx),%%xmm3	\n\t"\
		"movaps	     (%%rdx),%%xmm6	\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm4		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"subpd	%%xmm6,%%xmm5		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm4,%%xmm7		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"addpd	%%xmm5,%%xmm6		\n\t"\
		"\n\t"\
		"movaps	%%xmm0,     (%%rbx)	\n\t"\
		"movaps	%%xmm4,     (%%rcx)	\n\t"\
		"movaps	%%xmm1,0x010(%%rbx)	\n\t"\
		"movaps	%%xmm5,0x010(%%rdx)	\n\t"\
		"movaps	%%xmm2,     (%%rax)	\n\t"\
		"movaps	%%xmm7,     (%%rdx)	\n\t"\
		"movaps	%%xmm3,0x010(%%rax)	\n\t"\
		"movaps	%%xmm6,0x010(%%rcx)	\n\t"\
		"\n\t"\
		"\n\t"/* SSE2_RADIX4_DIF_IN_PLACE(r04,r24,r14,r34): */\
		"\n\t"\
		"subq	$0x0C0,%%rax	\n\t"/* r04 */\
		"addq	$0x0C0,%%rbx	\n\t"/* r24 */\
		"subq	$0x1C0,%%rcx	\n\t"/* r14 */\
		"subq	$0x040,%%rdx	\n\t"/* r34 */\
		"\n\t"\
		"movaps	     (%%rax),%%xmm0	\n\t"\
		"movaps	0x010(%%rax),%%xmm1	\n\t"\
		"movaps	     (%%rax),%%xmm2	\n\t"\
		"movaps	0x010(%%rax),%%xmm3	\n\t"\
		"addpd	     (%%rbx),%%xmm0	\n\t"\
		"addpd	0x010(%%rbx),%%xmm1	\n\t"\
		"subpd	     (%%rbx),%%xmm2	\n\t"\
		"subpd	0x010(%%rbx),%%xmm3	\n\t"\
		"\n\t"\
		"movaps	     (%%rcx),%%xmm4	\n\t"\
		"movaps	0x010(%%rcx),%%xmm5	\n\t"\
		"movaps	     (%%rcx),%%xmm6	\n\t"\
		"movaps	0x010(%%rcx),%%xmm7	\n\t"\
		"addpd	     (%%rdx),%%xmm4	\n\t"\
		"addpd	0x010(%%rdx),%%xmm5	\n\t"\
		"subpd	     (%%rdx),%%xmm6	\n\t"\
		"subpd	0x010(%%rdx),%%xmm7	\n\t"\
		"\n\t"/* Finish radix-4 butterfly and store results into temporary-array slots: */\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"movaps	%%xmm0,     (%%rbx)	\n\t"\
		"movaps	%%xmm1,0x010(%%rbx)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"movaps	%%xmm4,     (%%rax)	\n\t"\
		"movaps	%%xmm5,0x010(%%rax)	\n\t"\
		"\n\t"\
		"subpd	%%xmm7,%%xmm2		\n\t"\
		"subpd	%%xmm6,%%xmm3		\n\t"\
		"movaps	%%xmm2,     (%%rcx)	\n\t"\
		"movaps	%%xmm3,0x010(%%rdx)	\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	%%xmm7,     (%%rdx)	\n\t"\
		"movaps	%%xmm6,0x010(%%rcx)	\n\t"\
		"\n\t"\
		"\n\t"/* SSE2_RADIX4_DIF_IN_PLACE_2NDOFTWO(r0C,r2C,r1C,r3C): */\
		"addq	$0x080,%%rax		\n\t"\
		"addq	$0x080,%%rbx		\n\t"\
		"addq	$0x080,%%rcx		\n\t"\
		"addq	$0x080,%%rdx		\n\t"\
		"movaps	     (%%rax),%%xmm0	\n\t"\
		"movaps	0x010(%%rax),%%xmm1	\n\t"\
		"movaps	     (%%rax),%%xmm2	\n\t"\
		"movaps	0x010(%%rax),%%xmm3	\n\t"\
		"addpd	     (%%rbx),%%xmm0	\n\t"\
		"addpd	0x010(%%rbx),%%xmm1	\n\t"\
		"subpd	     (%%rbx),%%xmm2	\n\t"\
		"subpd	0x010(%%rbx),%%xmm3	\n\t"\
		"movaps	     (%%rcx),%%xmm4	\n\t"\
		"movaps	0x010(%%rcx),%%xmm5	\n\t"\
		"movaps	     (%%rcx),%%xmm6	\n\t"\
		"movaps	0x010(%%rcx),%%xmm7	\n\t"\
		"addpd	     (%%rdx),%%xmm6	\n\t"\
		"addpd	0x010(%%rdx),%%xmm7	\n\t"\
		"subpd	     (%%rdx),%%xmm4	\n\t"\
		"subpd	0x010(%%rdx),%%xmm5	\n\t"\
		"\n\t"/* Finish radix-4 butterfly and store results into temporary-array slots: */\
		"subpd	%%xmm6,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm1		\n\t"\
		"subpd	%%xmm5,%%xmm2		\n\t"\
		"movaps	%%xmm0,     (%%rbx)	\n\t"\
		"movaps	%%xmm1,0x010(%%rbx)	\n\t"\
		"subpd	%%xmm4,%%xmm3		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm0,%%xmm6		\n\t"\
		"addpd	%%xmm2,%%xmm5		\n\t"\
		"addpd	%%xmm1,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm4		\n\t"\
		"movaps	%%xmm6,     (%%rax)	\n\t"\
		"movaps	%%xmm7,0x010(%%rax)	\n\t"\
		"\n\t"\
		"movaps	(%%rdi),%%xmm0	\n\t"/* isrt2 */\
		"movaps	%%xmm2,%%xmm6		\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t"\
		"subpd	%%xmm4,%%xmm2		\n\t"\
		"subpd	%%xmm3,%%xmm5		\n\t"\
		"addpd	%%xmm4,%%xmm6		\n\t"\
		"addpd	%%xmm3,%%xmm7		\n\t"\
		"mulpd	%%xmm0,%%xmm2		\n\t"\
		"mulpd	%%xmm0,%%xmm5		\n\t"\
		"mulpd	%%xmm0,%%xmm6		\n\t"\
		"mulpd	%%xmm0,%%xmm7		\n\t"\
		"movaps	%%xmm2,     (%%rcx)	\n\t"\
		"movaps	%%xmm5,     (%%rdx)	\n\t"\
		"movaps	%%xmm6,0x010(%%rcx)	\n\t"\
		"movaps	%%xmm7,0x010(%%rdx)	\n\t"\
		"\n\t"/* SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS(r04,r14,r24,r34,r0C,r1C,r2C,r3C): */\
		"subq	$0x080,%%rax	\n\t"/* r04 */\
		"subq	$0x200,%%rbx	\n\t"/* r0C */\
		"addq	$0x080,%%rcx	\n\t"/* r24 */\
		"subq	$0x100,%%rdx	\n\t"/* r2C */\
		"movaps	     (%%rax),%%xmm0	\n\t"\
		"movaps	     (%%rcx),%%xmm4	\n\t"\
		"movaps	0x010(%%rax),%%xmm1	\n\t"\
		"movaps	0x010(%%rcx),%%xmm5	\n\t"\
		"movaps	     (%%rbx),%%xmm2	\n\t"\
		"movaps	0x010(%%rdx),%%xmm7	\n\t"\
		"movaps	0x010(%%rbx),%%xmm3	\n\t"\
		"movaps	     (%%rdx),%%xmm6	\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm4		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"subpd	%%xmm6,%%xmm5		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm4,%%xmm7		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"addpd	%%xmm5,%%xmm6		\n\t"\
		"\n\t"\
		"movaps	%%xmm0,     (%%rbx)	\n\t"\
		"movaps	%%xmm4,     (%%rcx)	\n\t"\
		"movaps	%%xmm1,0x010(%%rbx)	\n\t"\
		"movaps	%%xmm5,0x010(%%rdx)	\n\t"\
		"movaps	%%xmm2,     (%%rax)	\n\t"\
		"movaps	%%xmm7,     (%%rdx)	\n\t"\
		"movaps	%%xmm3,0x010(%%rax)	\n\t"\
		"movaps	%%xmm6,0x010(%%rcx)	\n\t"\
		"\n\t"\
		"addq	$0x100,%%rax	\n\t"/* r14 */\
		"addq	$0x100,%%rbx	\n\t"/* r1C */\
		"addq	$0x100,%%rcx	\n\t"/* r34 */\
		"addq	$0x100,%%rdx	\n\t"/* r3C */\
		"\n\t"\
		"movaps	     (%%rax),%%xmm0	\n\t"\
		"movaps	     (%%rcx),%%xmm4	\n\t"\
		"movaps	0x010(%%rax),%%xmm1	\n\t"\
		"movaps	0x010(%%rcx),%%xmm5	\n\t"\
		"movaps	     (%%rbx),%%xmm2	\n\t"\
		"movaps	0x010(%%rdx),%%xmm7	\n\t"\
		"movaps	0x010(%%rbx),%%xmm3	\n\t"\
		"movaps	     (%%rdx),%%xmm6	\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm4		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"subpd	%%xmm6,%%xmm5		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm4,%%xmm7		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"addpd	%%xmm5,%%xmm6		\n\t"\
		"\n\t"\
		"movaps	%%xmm0,     (%%rbx)	\n\t"\
		"movaps	%%xmm4,     (%%rcx)	\n\t"\
		"movaps	%%xmm1,0x010(%%rbx)	\n\t"\
		"movaps	%%xmm5,0x010(%%rdx)	\n\t"\
		"movaps	%%xmm2,     (%%rax)	\n\t"\
		"movaps	%%xmm7,     (%%rdx)	\n\t"\
		"movaps	%%xmm3,0x010(%%rax)	\n\t"\
		"movaps	%%xmm6,0x010(%%rcx)	\n\t"\
		"\n\t"\
		"\n\t"/* SSE2_RADIX4_DIF_IN_PLACE(r02,r22,r12,r32): */\
		"\n\t"\
		"subq	$0x120,%%rax	\n\t"/* r02 */\
		"addq	$0x060,%%rbx	\n\t"/* r22 */\
		"subq	$0x220,%%rcx	\n\t"/* r12 */\
		"subq	$0x0a0,%%rdx	\n\t"/* r32 */\
		"\n\t"\
		"movaps	     (%%rax),%%xmm0	\n\t"\
		"movaps	0x010(%%rax),%%xmm1	\n\t"\
		"movaps	     (%%rax),%%xmm2	\n\t"\
		"movaps	0x010(%%rax),%%xmm3	\n\t"\
		"addpd	     (%%rbx),%%xmm0	\n\t"\
		"addpd	0x010(%%rbx),%%xmm1	\n\t"\
		"subpd	     (%%rbx),%%xmm2	\n\t"\
		"subpd	0x010(%%rbx),%%xmm3	\n\t"\
		"\n\t"\
		"movaps	     (%%rcx),%%xmm4	\n\t"\
		"movaps	0x010(%%rcx),%%xmm5	\n\t"\
		"movaps	     (%%rcx),%%xmm6	\n\t"\
		"movaps	0x010(%%rcx),%%xmm7	\n\t"\
		"addpd	     (%%rdx),%%xmm4	\n\t"\
		"addpd	0x010(%%rdx),%%xmm5	\n\t"\
		"subpd	     (%%rdx),%%xmm6	\n\t"\
		"subpd	0x010(%%rdx),%%xmm7	\n\t"\
		"\n\t"/* Finish radix-4 butterfly and store results into temporary-array slots: */\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"movaps	%%xmm0,     (%%rbx)	\n\t"\
		"movaps	%%xmm1,0x010(%%rbx)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"movaps	%%xmm4,     (%%rax)	\n\t"\
		"movaps	%%xmm5,0x010(%%rax)	\n\t"\
		"\n\t"\
		"subpd	%%xmm7,%%xmm2		\n\t"\
		"subpd	%%xmm6,%%xmm3		\n\t"\
		"movaps	%%xmm2,     (%%rcx)	\n\t"\
		"movaps	%%xmm3,0x010(%%rdx)	\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	%%xmm7,     (%%rdx)	\n\t"\
		"movaps	%%xmm6,0x010(%%rcx)	\n\t"\
		"\n\t"\
		"\n\t"/* SSE2_RADIX4_DIF_IN_PLACE_2NDOFTWO(r0A,r2A,r1A,r3A): */\
		"addq	$0x080,%%rax		\n\t"\
		"addq	$0x080,%%rbx		\n\t"\
		"addq	$0x080,%%rcx		\n\t"\
		"addq	$0x080,%%rdx		\n\t"\
		"movaps	     (%%rax),%%xmm0	\n\t"\
		"movaps	0x010(%%rax),%%xmm1	\n\t"\
		"movaps	     (%%rax),%%xmm2	\n\t"\
		"movaps	0x010(%%rax),%%xmm3	\n\t"\
		"addpd	     (%%rbx),%%xmm0	\n\t"\
		"addpd	0x010(%%rbx),%%xmm1	\n\t"\
		"subpd	     (%%rbx),%%xmm2	\n\t"\
		"subpd	0x010(%%rbx),%%xmm3	\n\t"\
		"movaps	     (%%rcx),%%xmm4	\n\t"\
		"movaps	0x010(%%rcx),%%xmm5	\n\t"\
		"movaps	     (%%rcx),%%xmm6	\n\t"\
		"movaps	0x010(%%rcx),%%xmm7	\n\t"\
		"addpd	     (%%rdx),%%xmm6	\n\t"\
		"addpd	0x010(%%rdx),%%xmm7	\n\t"\
		"subpd	     (%%rdx),%%xmm4	\n\t"\
		"subpd	0x010(%%rdx),%%xmm5	\n\t"\
		"\n\t"/* Finish radix-4 butterfly and store results into temporary-array slots: */\
		"subpd	%%xmm6,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm1		\n\t"\
		"subpd	%%xmm5,%%xmm2		\n\t"\
		"movaps	%%xmm0,     (%%rbx)	\n\t"\
		"movaps	%%xmm1,0x010(%%rbx)	\n\t"\
		"subpd	%%xmm4,%%xmm3		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm0,%%xmm6		\n\t"\
		"addpd	%%xmm2,%%xmm5		\n\t"\
		"addpd	%%xmm1,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm4		\n\t"\
		"movaps	%%xmm6,     (%%rax)	\n\t"\
		"movaps	%%xmm7,0x010(%%rax)	\n\t"\
		"\n\t"\
		"movaps	(%%rdi),%%xmm0	\n\t"/* isrt2 */\
		"movaps	%%xmm2,%%xmm6		\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t"\
		"subpd	%%xmm4,%%xmm2		\n\t"\
		"subpd	%%xmm3,%%xmm5		\n\t"\
		"addpd	%%xmm4,%%xmm6		\n\t"\
		"addpd	%%xmm3,%%xmm7		\n\t"\
		"mulpd	%%xmm0,%%xmm2		\n\t"\
		"mulpd	%%xmm0,%%xmm5		\n\t"\
		"mulpd	%%xmm0,%%xmm6		\n\t"\
		"mulpd	%%xmm0,%%xmm7		\n\t"\
		"movaps	%%xmm2,     (%%rcx)	\n\t"\
		"movaps	%%xmm5,     (%%rdx)	\n\t"\
		"movaps	%%xmm6,0x010(%%rcx)	\n\t"\
		"movaps	%%xmm7,0x010(%%rdx)	\n\t"\
		"\n\t"/* SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS(r02,r12,r22,r32,r0A,r1A,r2A,r3A): */\
		"subq	$0x080,%%rax	\n\t"/* r02 */\
		"subq	$0x200,%%rbx	\n\t"/* r0A */\
		"addq	$0x080,%%rcx	\n\t"/* r22 */\
		"subq	$0x100,%%rdx	\n\t"/* r2A */\
		"movaps	     (%%rax),%%xmm0	\n\t"\
		"movaps	     (%%rcx),%%xmm4	\n\t"\
		"movaps	0x010(%%rax),%%xmm1	\n\t"\
		"movaps	0x010(%%rcx),%%xmm5	\n\t"\
		"movaps	     (%%rbx),%%xmm2	\n\t"\
		"movaps	0x010(%%rdx),%%xmm7	\n\t"\
		"movaps	0x010(%%rbx),%%xmm3	\n\t"\
		"movaps	     (%%rdx),%%xmm6	\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm4		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"subpd	%%xmm6,%%xmm5		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm4,%%xmm7		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"addpd	%%xmm5,%%xmm6		\n\t"\
		"\n\t"\
		"movaps	%%xmm0,     (%%rbx)	\n\t"\
		"movaps	%%xmm4,     (%%rcx)	\n\t"\
		"movaps	%%xmm1,0x010(%%rbx)	\n\t"\
		"movaps	%%xmm5,0x010(%%rdx)	\n\t"\
		"movaps	%%xmm2,     (%%rax)	\n\t"\
		"movaps	%%xmm7,     (%%rdx)	\n\t"\
		"movaps	%%xmm3,0x010(%%rax)	\n\t"\
		"movaps	%%xmm6,0x010(%%rcx)	\n\t"\
		"\n\t"\
		"addq	$0x100,%%rax	\n\t"/* r12 */\
		"addq	$0x100,%%rbx	\n\t"/* r1A */\
		"addq	$0x100,%%rcx	\n\t"/* r32 */\
		"addq	$0x100,%%rdx	\n\t"/* r3A */\
		"\n\t"\
		"movaps	     (%%rax),%%xmm0	\n\t"\
		"movaps	     (%%rcx),%%xmm4	\n\t"\
		"movaps	0x010(%%rax),%%xmm1	\n\t"\
		"movaps	0x010(%%rcx),%%xmm5	\n\t"\
		"movaps	     (%%rbx),%%xmm2	\n\t"\
		"movaps	0x010(%%rdx),%%xmm7	\n\t"\
		"movaps	0x010(%%rbx),%%xmm3	\n\t"\
		"movaps	     (%%rdx),%%xmm6	\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm4		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"subpd	%%xmm6,%%xmm5		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm4,%%xmm7		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"addpd	%%xmm5,%%xmm6		\n\t"\
		"\n\t"\
		"movaps	%%xmm0,     (%%rbx)	\n\t"\
		"movaps	%%xmm4,     (%%rcx)	\n\t"\
		"movaps	%%xmm1,0x010(%%rbx)	\n\t"\
		"movaps	%%xmm5,0x010(%%rdx)	\n\t"\
		"movaps	%%xmm2,     (%%rax)	\n\t"\
		"movaps	%%xmm7,     (%%rdx)	\n\t"\
		"movaps	%%xmm3,0x010(%%rax)	\n\t"\
		"movaps	%%xmm6,0x010(%%rcx)	\n\t"\
		"\n\t"\
		"\n\t"/* SSE2_RADIX4_DIF_IN_PLACE(r06,r26,r16,r36): */\
		"\n\t"\
		"subq	$0x0C0,%%rax	\n\t"/* r06 */\
		"addq	$0x0C0,%%rbx	\n\t"/* r26 */\
		"subq	$0x1C0,%%rcx	\n\t"/* r16 */\
		"subq	$0x040,%%rdx	\n\t"/* r36 */\
		"\n\t"\
		"movaps	     (%%rax),%%xmm0	\n\t"\
		"movaps	0x010(%%rax),%%xmm1	\n\t"\
		"movaps	     (%%rax),%%xmm2	\n\t"\
		"movaps	0x010(%%rax),%%xmm3	\n\t"\
		"addpd	     (%%rbx),%%xmm0	\n\t"\
		"addpd	0x010(%%rbx),%%xmm1	\n\t"\
		"subpd	     (%%rbx),%%xmm2	\n\t"\
		"subpd	0x010(%%rbx),%%xmm3	\n\t"\
		"\n\t"\
		"movaps	     (%%rcx),%%xmm4	\n\t"\
		"movaps	0x010(%%rcx),%%xmm5	\n\t"\
		"movaps	     (%%rcx),%%xmm6	\n\t"\
		"movaps	0x010(%%rcx),%%xmm7	\n\t"\
		"addpd	     (%%rdx),%%xmm4	\n\t"\
		"addpd	0x010(%%rdx),%%xmm5	\n\t"\
		"subpd	     (%%rdx),%%xmm6	\n\t"\
		"subpd	0x010(%%rdx),%%xmm7	\n\t"\
		"\n\t"/* Finish radix-4 butterfly and store results into temporary-array slots: */\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"movaps	%%xmm0,     (%%rbx)	\n\t"\
		"movaps	%%xmm1,0x010(%%rbx)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"movaps	%%xmm4,     (%%rax)	\n\t"\
		"movaps	%%xmm5,0x010(%%rax)	\n\t"\
		"\n\t"\
		"subpd	%%xmm7,%%xmm2		\n\t"\
		"subpd	%%xmm6,%%xmm3		\n\t"\
		"movaps	%%xmm2,     (%%rcx)	\n\t"\
		"movaps	%%xmm3,0x010(%%rdx)	\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	%%xmm7,     (%%rdx)	\n\t"\
		"movaps	%%xmm6,0x010(%%rcx)	\n\t"\
		"\n\t"\
		"\n\t"/* SSE2_RADIX4_DIF_IN_PLACE_2NDOFTWO(r0E,r2E,r1E,r3E): */\
		"addq	$0x080,%%rax		\n\t"\
		"addq	$0x080,%%rbx		\n\t"\
		"addq	$0x080,%%rcx		\n\t"\
		"addq	$0x080,%%rdx		\n\t"\
		"movaps	     (%%rax),%%xmm0	\n\t"\
		"movaps	0x010(%%rax),%%xmm1	\n\t"\
		"movaps	     (%%rax),%%xmm2	\n\t"\
		"movaps	0x010(%%rax),%%xmm3	\n\t"\
		"addpd	     (%%rbx),%%xmm0	\n\t"\
		"addpd	0x010(%%rbx),%%xmm1	\n\t"\
		"subpd	     (%%rbx),%%xmm2	\n\t"\
		"subpd	0x010(%%rbx),%%xmm3	\n\t"\
		"movaps	     (%%rcx),%%xmm4	\n\t"\
		"movaps	0x010(%%rcx),%%xmm5	\n\t"\
		"movaps	     (%%rcx),%%xmm6	\n\t"\
		"movaps	0x010(%%rcx),%%xmm7	\n\t"\
		"addpd	     (%%rdx),%%xmm6	\n\t"\
		"addpd	0x010(%%rdx),%%xmm7	\n\t"\
		"subpd	     (%%rdx),%%xmm4	\n\t"\
		"subpd	0x010(%%rdx),%%xmm5	\n\t"\
		"\n\t"/* Finish radix-4 butterfly and store results into temporary-array slots: */\
		"subpd	%%xmm6,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm1		\n\t"\
		"subpd	%%xmm5,%%xmm2		\n\t"\
		"movaps	%%xmm0,     (%%rbx)	\n\t"\
		"movaps	%%xmm1,0x010(%%rbx)	\n\t"\
		"subpd	%%xmm4,%%xmm3		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm0,%%xmm6		\n\t"\
		"addpd	%%xmm2,%%xmm5		\n\t"\
		"addpd	%%xmm1,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm4		\n\t"\
		"movaps	%%xmm6,     (%%rax)	\n\t"\
		"movaps	%%xmm7,0x010(%%rax)	\n\t"\
		"\n\t"\
		"movaps	(%%rdi),%%xmm0	\n\t"/* isrt2 */\
		"movaps	%%xmm2,%%xmm6		\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t"\
		"subpd	%%xmm4,%%xmm2		\n\t"\
		"subpd	%%xmm3,%%xmm5		\n\t"\
		"addpd	%%xmm4,%%xmm6		\n\t"\
		"addpd	%%xmm3,%%xmm7		\n\t"\
		"mulpd	%%xmm0,%%xmm2		\n\t"\
		"mulpd	%%xmm0,%%xmm5		\n\t"\
		"mulpd	%%xmm0,%%xmm6		\n\t"\
		"mulpd	%%xmm0,%%xmm7		\n\t"\
		"movaps	%%xmm2,     (%%rcx)	\n\t"\
		"movaps	%%xmm5,     (%%rdx)	\n\t"\
		"movaps	%%xmm6,0x010(%%rcx)	\n\t"\
		"movaps	%%xmm7,0x010(%%rdx)	\n\t"\
		"\n\t"/* SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS(r06,r16,r26,r36,r0E,r1E,r2E,r3E): */\
		"subq	$0x080,%%rax	\n\t"/* r02 */\
		"subq	$0x200,%%rbx	\n\t"/* r0A */\
		"addq	$0x080,%%rcx	\n\t"/* r22 */\
		"subq	$0x100,%%rdx	\n\t"/* r2A */\
		"movaps	     (%%rax),%%xmm0	\n\t"\
		"movaps	     (%%rcx),%%xmm4	\n\t"\
		"movaps	0x010(%%rax),%%xmm1	\n\t"\
		"movaps	0x010(%%rcx),%%xmm5	\n\t"\
		"movaps	     (%%rbx),%%xmm2	\n\t"\
		"movaps	0x010(%%rdx),%%xmm7	\n\t"\
		"movaps	0x010(%%rbx),%%xmm3	\n\t"\
		"movaps	     (%%rdx),%%xmm6	\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm4		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"subpd	%%xmm6,%%xmm5		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm4,%%xmm7		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"addpd	%%xmm5,%%xmm6		\n\t"\
		"\n\t"\
		"movaps	%%xmm0,     (%%rbx)	\n\t"\
		"movaps	%%xmm4,     (%%rcx)	\n\t"\
		"movaps	%%xmm1,0x010(%%rbx)	\n\t"\
		"movaps	%%xmm5,0x010(%%rdx)	\n\t"\
		"movaps	%%xmm2,     (%%rax)	\n\t"\
		"movaps	%%xmm7,     (%%rdx)	\n\t"\
		"movaps	%%xmm3,0x010(%%rax)	\n\t"\
		"movaps	%%xmm6,0x010(%%rcx)	\n\t"\
		"\n\t"\
		"addq	$0x100,%%rax	\n\t"/* r12 */\
		"addq	$0x100,%%rbx	\n\t"/* r1A */\
		"addq	$0x100,%%rcx	\n\t"/* r32 */\
		"addq	$0x100,%%rdx	\n\t"/* r3A */\
		"\n\t"\
		"movaps	     (%%rax),%%xmm0	\n\t"\
		"movaps	     (%%rcx),%%xmm4	\n\t"\
		"movaps	0x010(%%rax),%%xmm1	\n\t"\
		"movaps	0x010(%%rcx),%%xmm5	\n\t"\
		"movaps	     (%%rbx),%%xmm2	\n\t"\
		"movaps	0x010(%%rdx),%%xmm7	\n\t"\
		"movaps	0x010(%%rbx),%%xmm3	\n\t"\
		"movaps	     (%%rdx),%%xmm6	\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm4		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"subpd	%%xmm6,%%xmm5		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm4,%%xmm7		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"addpd	%%xmm5,%%xmm6		\n\t"\
		"\n\t"\
		"movaps	%%xmm0,     (%%rbx)	\n\t"\
		"movaps	%%xmm4,     (%%rcx)	\n\t"\
		"movaps	%%xmm1,0x010(%%rbx)	\n\t"\
		"movaps	%%xmm5,0x010(%%rdx)	\n\t"\
		"movaps	%%xmm2,     (%%rax)	\n\t"\
		"movaps	%%xmm7,     (%%rdx)	\n\t"\
		"movaps	%%xmm3,0x010(%%rax)	\n\t"\
		"movaps	%%xmm6,0x010(%%rcx)	\n\t"\
		"\n\t"\
		"\n\t"/*...and now do eight radix-4 transforms, including the internal twiddle factors: */\
		"\n\t"/*...Block 1: t00,t10,t20,t30 in r00,04,02,06 - note swapped middle 2 indices! */\
		"\n\t"\
		"movq	%[__r00],%%rsi	\n\t"\
		"movq	%[__add],%%rax	\n\t"\
		"movslq	%[__p01],%%rbx	\n\t"\
		"movslq	%[__p02],%%rcx	\n\t"\
		"movslq	%[__p03],%%rdx	\n\t"\
		"shlq	$3,%%rbx		\n\t"\
		"shlq	$3,%%rcx		\n\t"\
		"shlq	$3,%%rdx		\n\t"\
		"addq	%%rax,%%rbx		\n\t"\
		"addq	%%rax,%%rcx		\n\t"\
		"addq	%%rax,%%rdx		\n\t"\
		"movaps	     (%%rsi),%%xmm0	\n\t"/* t00 */\
		"movaps	0x020(%%rsi),%%xmm4	\n\t"/* t20 */\
		"movaps	0x010(%%rsi),%%xmm1	\n\t"/* t01 */\
		"movaps	0x030(%%rsi),%%xmm5	\n\t"/* t21 */\
		"movaps	0x040(%%rsi),%%xmm2	\n\t"/* t10 */\
		"movaps	0x060(%%rsi),%%xmm6	\n\t"/* t30 */\
		"movaps	0x050(%%rsi),%%xmm3	\n\t"/* t11 */\
		"movaps	0x070(%%rsi),%%xmm7	\n\t"/* t31 */\
		"subpd	0x040(%%rsi),%%xmm0	\n\t"/* t10=t00-rt */\
		"subpd	0x060(%%rsi),%%xmm4	\n\t"/* t30=t20-rt */\
		"subpd	0x050(%%rsi),%%xmm1	\n\t"/* t11=t01-it */\
		"subpd	0x070(%%rsi),%%xmm5	\n\t"/* t31=t21-it */\
		"addpd	     (%%rsi),%%xmm2	\n\t"/* t00=t00+rt */\
		"addpd	0x020(%%rsi),%%xmm6	\n\t"/* t20=t20+rt */\
		"addpd	0x010(%%rsi),%%xmm3	\n\t"/* t01=t01+it */\
		"addpd	0x030(%%rsi),%%xmm7	\n\t"/* t21=t21+it */\
		"subpd	%%xmm6,%%xmm2		\n\t"/* t00 <- t00-t20 */\
		"subpd	%%xmm5,%%xmm0		\n\t"/* t10 <- t10-t31 */\
		"subpd	%%xmm7,%%xmm3		\n\t"/* t01 <- t01-t21 */\
		"subpd	%%xmm4,%%xmm1		\n\t"/* t11 <- t11-t30 */\
		"addpd	%%xmm6,%%xmm6		\n\t"/*          2*t20 */\
		"addpd	%%xmm5,%%xmm5		\n\t"/*          2*t31 */\
		"addpd	%%xmm7,%%xmm7		\n\t"/*          2*t21 */\
		"addpd	%%xmm4,%%xmm4		\n\t"/*          2*t30 */\
		"movaps	%%xmm2,     (%%rbx)	\n\t"/* a(jt+p1 ) */\
		"movaps	%%xmm0,     (%%rcx)	\n\t"/* a(jt+p2 ) */\
		"movaps	%%xmm3,0x010(%%rbx)	\n\t"/* a(jp+p1 ) */\
		"movaps	%%xmm1,0x010(%%rdx)	\n\t"/* a(jp+p3 ) */\
		"addpd	%%xmm2,%%xmm6		\n\t"/* t20 <- t00+t20 */\
		"addpd	%%xmm0,%%xmm5		\n\t"/* t31 <- t10+t31 */\
		"addpd	%%xmm3,%%xmm7		\n\t"/* t21 <- t01+t21 */\
		"addpd	%%xmm1,%%xmm4		\n\t"/* t30 <- t11+t30 */\
		"movaps	%%xmm6,     (%%rax)	\n\t"/* a(jt+p0 ) */\
		"movaps	%%xmm5,     (%%rdx)	\n\t"/* a(jt+p3 ) */\
		"movaps	%%xmm7,0x010(%%rax)	\n\t"/* a(jp+p0 ) */\
		"movaps	%%xmm4,0x010(%%rcx)	\n\t"/* a(jp+p2 ) */\
		"\n\t"/*...Block 5: t08,t18,t28,t38	*/\
		"addq	$0x080,%%rsi			\n\t"/* r08 */\
		"movslq	%[__p04],%%rdi		\n\t"\
		"shlq	$3,%%rdi			\n\t"\
		"addq	%%rdi,%%rax			\n\t"/* ap04 */\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	%%rdi,%%rcx			\n\t"\
		"addq	%%rdi,%%rdx			\n\t"\
		"movq	%[__isrt2],%%rdi	\n\t"\
		"\n\t"\
		"movaps	0x020(%%rsi),%%xmm4	\n\t"/* t28 */\
		"movaps	0x030(%%rsi),%%xmm5	\n\t"/* t29 */\
		"movaps	0x060(%%rsi),%%xmm6	\n\t"/* t38 */\
		"movaps	0x070(%%rsi),%%xmm7	\n\t"/* t39 */\
		"movaps	     (%%rsi),%%xmm0	\n\t"/* t08 */\
		"subpd	0x030(%%rsi),%%xmm4	\n\t"/* t28-t29 */\
		"movaps	0x010(%%rsi),%%xmm1	\n\t"/* t09 */\
		"addpd	0x020(%%rsi),%%xmm5	\n\t"/* t29+t28 */\
		"movaps	0x040(%%rsi),%%xmm2	\n\t"/* t18 */\
		"mulpd	     (%%rdi),%%xmm4		\n\t"/* t28 = (t28-t29)*ISRT2 */\
		"movaps	0x050(%%rsi),%%xmm3	\n\t"/* t19 */\
		"mulpd	     (%%rdi),%%xmm5		\n\t"/* t29 = (t29+t28)*ISRT2 */\
		"subpd	0x050(%%rsi),%%xmm0	\n\t"/* t08=t08-t19*/\
		"addpd	0x070(%%rsi),%%xmm6	\n\t"/* t38+t39 */\
		"subpd	0x040(%%rsi),%%xmm1	\n\t"/* t19=t09-t18*/\
		"subpd	0x060(%%rsi),%%xmm7	\n\t"/* t39-t38 */\
		"addpd	0x010(%%rsi),%%xmm2	\n\t"/* t09=t18+t09*/\
		"mulpd	     (%%rdi),%%xmm6		\n\t"/*  rt = (t38+t39)*ISRT2 */\
		"addpd	     (%%rsi),%%xmm3	\n\t"/* t18=t19+t08*/\
		"mulpd	     (%%rdi),%%xmm7		\n\t"/*  it = (t39-t38)*ISRT2 */\
		"\n\t"\
		"subpd	%%xmm6,%%xmm4		\n\t"/* t28=t28-rt */\
		"subpd	%%xmm7,%%xmm5		\n\t"/* t29=t29-it */\
		"addpd	%%xmm6,%%xmm6		\n\t"/*      2* rt */\
		"addpd	%%xmm7,%%xmm7		\n\t"/*      2* it */\
		"addpd	%%xmm4,%%xmm6		\n\t"/* t38=t28+rt */\
		"addpd	%%xmm5,%%xmm7		\n\t"/* t39=t29+it */\
		"subpd	%%xmm4,%%xmm0		\n\t"/* t08-t28 */\
		"subpd	%%xmm5,%%xmm2		\n\t"/* t09-t29 */\
		"addpd	%%xmm4,%%xmm4		\n\t"/*   2*t28 */\
		"addpd	%%xmm5,%%xmm5		\n\t"/*   2*t29 */\
		"\n\t"\
		"subpd	%%xmm7,%%xmm3		\n\t"/* t18-t39 */\
		"subpd	%%xmm6,%%xmm1		\n\t"/* t19-t38 */\
		"addpd	%%xmm7,%%xmm7		\n\t"/*   2*t39 */\
		"addpd	%%xmm6,%%xmm6		\n\t"/*   2*t38 */\
		"movaps	%%xmm0,     (%%rbx)	\n\t"/* a(jt+p1 ) */\
		"movaps	%%xmm2,0x010(%%rbx)	\n\t"/* a(jp+p1 ) */\
		"movaps	%%xmm3,     (%%rcx)	\n\t"/* a(jt+p2 ) */\
		"movaps	%%xmm1,0x010(%%rdx)	\n\t"/* a(jp+p3 ) */\
		"addpd	%%xmm0,%%xmm4		\n\t"/* t08+t28 */\
		"addpd	%%xmm2,%%xmm5		\n\t"/* t09+t29 */\
		"addpd	%%xmm3,%%xmm7		\n\t"/* t18+t39 */\
		"addpd	%%xmm1,%%xmm6		\n\t"/* t19+t38 */\
		"movaps	%%xmm4,     (%%rax)	\n\t"/* a(jt+p0 ) */\
		"movaps	%%xmm5,0x010(%%rax)	\n\t"/* a(jp+p0 ) */\
		"movaps	%%xmm7,     (%%rdx)	\n\t"/* a(jt+p3 ) */\
		"movaps	%%xmm6,0x010(%%rcx)	\n\t"/* a(jp+p2 ) */\
		"\n\t"/*...Block 3: t04,t14,t24,t34	*/\
		"addq	$0x180,%%rsi		\n\t"/* r20 */\
		"movslq	%[__p04],%%rdi		\n\t"\
		"shlq	$3,%%rdi			\n\t"\
		"addq	%%rdi,%%rax			\n\t"/* ap08 */\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	%%rdi,%%rcx			\n\t"\
		"addq	%%rdi,%%rdx			\n\t"\
		"movq	%[__cc0],%%rdi	\n\t"\
		"\n\t"\
		"movaps	0x020(%%rsi),%%xmm4	\n\t"/* t24 */\
		"movaps	0x060(%%rsi),%%xmm6	\n\t"/* t34 */\
		"movaps	0x030(%%rsi),%%xmm5	\n\t"/* t25 */\
		"movaps	0x070(%%rsi),%%xmm7	\n\t"/* t35 */\
		"movaps	0x020(%%rsi),%%xmm0	\n\t"/* copy t24 */\
		"movaps	0x060(%%rsi),%%xmm2	\n\t"/* copy t34 */\
		"movaps	0x030(%%rsi),%%xmm1	\n\t"/* copy t25 */\
		"movaps	0x070(%%rsi),%%xmm3	\n\t"/* copy t35 */\
		"mulpd	     (%%rdi),%%xmm4	\n\t"/* t24*c */\
		"mulpd	0x010(%%rdi),%%xmm6	\n\t"/* t34*s */\
		"mulpd	0x010(%%rdi),%%xmm1	\n\t"/* t25*s */\
		"mulpd	     (%%rdi),%%xmm3	\n\t"/* t35*c */\
		"mulpd	     (%%rdi),%%xmm5	\n\t"/* t25*c */\
		"mulpd	0x010(%%rdi),%%xmm7	\n\t"/* t35*s */\
		"mulpd	0x010(%%rdi),%%xmm0	\n\t"/* t24*s */\
		"mulpd	     (%%rdi),%%xmm2	\n\t"/* t34*c */\
		"subpd	%%xmm1,%%xmm4		\n\t"/* ~t24 */\
		"subpd	%%xmm3,%%xmm6		\n\t"/* rt */\
		"addpd	%%xmm0,%%xmm5		\n\t"/* ~t25 */\
		"addpd	%%xmm2,%%xmm7		\n\t"/* it */\
		"subq	$0x010,%%rdi			\n\t"/* isrt2 */\
		"\n\t"\
		"movaps	0x040(%%rsi),%%xmm2	\n\t"/* t14 */\
		"subpd	%%xmm6,%%xmm4		\n\t"/*~t34=t24-rt */\
		"movaps	0x050(%%rsi),%%xmm3	\n\t"/* t15 */\
		"subpd	%%xmm7,%%xmm5		\n\t"/*~t35=t25-it */\
		"subpd	0x050(%%rsi),%%xmm2	\n\t"/* t14-t15 */\
		"addpd	%%xmm6,%%xmm6		\n\t"/*      2* rt */\
		"addpd	0x040(%%rsi),%%xmm3	\n\t"/* t15+t14 */\
		"addpd	%%xmm7,%%xmm7		\n\t"/*      2* it */\
		"mulpd	     (%%rdi),%%xmm2		\n\t"/* rt = (t14-t15)*ISRT2 */\
		"addpd	%%xmm4,%%xmm6		\n\t"/*~t24=t24+rt */\
		"mulpd	     (%%rdi),%%xmm3		\n\t"/* it = (t15+t14)*ISRT2 */\
		"addpd	%%xmm5,%%xmm7		\n\t"/*~t25=t25+it */\
		"movaps	     (%%rsi),%%xmm0	\n\t"/* t04 */\
		"movaps	0x010(%%rsi),%%xmm1	\n\t"/* t05 */\
		"subpd	%%xmm2,%%xmm0		\n\t"/*~t14=t04-rt */\
		"subpd	%%xmm3,%%xmm1		\n\t"/*~t15=t05-it */\
		"addpd	     (%%rsi),%%xmm2	\n\t"/*~t04=rt +t04*/\
		"addpd	0x010(%%rsi),%%xmm3	\n\t"/*~t05=it +t05*/\
		"subpd	%%xmm6,%%xmm2		\n\t"/* t04-t24 */\
		"subpd	%%xmm5,%%xmm0		\n\t"/* t14-t35 */\
		"subpd	%%xmm7,%%xmm3		\n\t"/* t05-t25 */\
		"subpd	%%xmm4,%%xmm1		\n\t"/* t15-t34 */\
		"addpd	%%xmm6,%%xmm6		\n\t"/*   2*t24 */\
		"addpd	%%xmm5,%%xmm5		\n\t"/*          2*t35 */\
		"addpd	%%xmm7,%%xmm7		\n\t"/*   2*t25 */\
		"addpd	%%xmm4,%%xmm4		\n\t"/*          2*t34 */\
		"movaps	%%xmm2,     (%%rbx)	\n\t"/* a(jt+p1 ) */\
		"movaps	%%xmm0,     (%%rcx)	\n\t"/* a(jt+p2 ) */\
		"movaps	%%xmm3,0x010(%%rbx)	\n\t"/* a(jp+p1 ) */\
		"movaps	%%xmm1,0x010(%%rdx)	\n\t"/* a(jp+p3 ) */\
		"addpd	%%xmm2,%%xmm6		\n\t"/* t04+t24 */\
		"addpd	%%xmm0,%%xmm5		\n\t"/* t14+t35 */\
		"addpd	%%xmm3,%%xmm7		\n\t"/* t05+t25 */\
		"addpd	%%xmm1,%%xmm4		\n\t"/* t15+t34 */\
		"movaps	%%xmm6,     (%%rax)	\n\t"/* a(jt+p0 ) */\
		"movaps	%%xmm5,     (%%rdx)	\n\t"/* a(jt+p3 ) */\
		"movaps	%%xmm7,0x010(%%rax)	\n\t"/* a(jp+p0 ) */\
		"movaps	%%xmm4,0x010(%%rcx)	\n\t"/* a(jp+p2 ) */\
		"\n\t"/*...Block 7: t0C,t1C,t2C,t3C	*/\
		"addq	$0x080,%%rsi			\n\t"/* r28 */\
		"movslq	%[__p04],%%rdi		\n\t"\
		"shlq	$3,%%rdi			\n\t"\
		"addq	%%rdi,%%rax			\n\t"/* ap0C */\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	%%rdi,%%rcx			\n\t"\
		"addq	%%rdi,%%rdx			\n\t"\
		"movq	%[__cc0],%%rdi	\n\t"\
		"\n\t"\
		"movaps	0x020(%%rsi),%%xmm4	\n\t"/* t2C */\
		"movaps	0x060(%%rsi),%%xmm6	\n\t"/* t3C */\
		"movaps	0x030(%%rsi),%%xmm5	\n\t"/* t2D */\
		"movaps	0x070(%%rsi),%%xmm7	\n\t"/* t3D */\
		"movaps	0x020(%%rsi),%%xmm0	\n\t"/* copy t2C */\
		"movaps	0x060(%%rsi),%%xmm2	\n\t"/* copy t3C */\
		"movaps	0x030(%%rsi),%%xmm1	\n\t"/* copy t2D */\
		"movaps	0x070(%%rsi),%%xmm3	\n\t"/* copy t3D */\
		"mulpd	0x010(%%rdi),%%xmm4	\n\t"/* t2C*s */\
		"mulpd	     (%%rdi),%%xmm6	\n\t"/* t3C*c */\
		"mulpd	     (%%rdi),%%xmm1	\n\t"/* t2D*c */\
		"mulpd	0x010(%%rdi),%%xmm3	\n\t"/* t3D*s */\
		"mulpd	0x010(%%rdi),%%xmm5	\n\t"/* t2D*s */\
		"mulpd	     (%%rdi),%%xmm7	\n\t"/* t3D*c */\
		"mulpd	     (%%rdi),%%xmm0	\n\t"/* t2C*c */\
		"mulpd	0x010(%%rdi),%%xmm2	\n\t"/* t3C*s */\
		"subpd	%%xmm1,%%xmm4		\n\t"/* ~t24 */\
		"subpd	%%xmm3,%%xmm6		\n\t"/* rt */\
		"addpd	%%xmm0,%%xmm5		\n\t"/* ~t25 */\
		"addpd	%%xmm2,%%xmm7		\n\t"/* it */\
		"subq	$0x010,%%rdi			\n\t"/* isrt2 */\
		"\n\t"\
		"movaps	0x040(%%rsi),%%xmm2	\n\t"/* t14 */\
		"subpd	%%xmm6,%%xmm4		\n\t"/*~t2C=t2C-rt */\
		"movaps	0x050(%%rsi),%%xmm3	\n\t"/* t1D */\
		"subpd	%%xmm7,%%xmm5		\n\t"/*~t2D=t2D-it */\
		"addpd	0x050(%%rsi),%%xmm2	\n\t"/* t1C+t1D */\
		"addpd	%%xmm6,%%xmm6		\n\t"/*      2* rt */\
		"subpd	0x040(%%rsi),%%xmm3	\n\t"/* t1D-t1C */\
		"addpd	%%xmm7,%%xmm7		\n\t"/*      2* it */\
		"mulpd	     (%%rdi),%%xmm2		\n\t"/* rt = (t1C+t1D)*ISRT2 */\
		"addpd	%%xmm4,%%xmm6		\n\t"/*~t3C=t2C+rt */\
		"mulpd	     (%%rdi),%%xmm3		\n\t"/* it = (t1D-t1C)*ISRT2 */\
		"addpd	%%xmm5,%%xmm7		\n\t"/*~t3D=t2D+it */\
		"movaps	     (%%rsi),%%xmm0	\n\t"/* t0C */\
		"movaps	0x010(%%rsi),%%xmm1	\n\t"/* t0D */\
		"subpd	%%xmm2,%%xmm0		\n\t"/*~t0C=t0C-rt */\
		"subpd	%%xmm3,%%xmm1		\n\t"/*~t0D=t0D-it */\
		"addpd	     (%%rsi),%%xmm2	\n\t"/*~t1C=rt +t0C*/\
		"addpd	0x010(%%rsi),%%xmm3	\n\t"/*~t1D=it +t0D*/\
		"subpd	%%xmm4,%%xmm0		\n\t"/* t0C-t2C */\
		"subpd	%%xmm7,%%xmm2		\n\t"/* t1C-t3D */\
		"subpd	%%xmm5,%%xmm1		\n\t"/* t0D-t2D */\
		"subpd	%%xmm6,%%xmm3		\n\t"/* t1D-t3C */\
		"addpd	%%xmm4,%%xmm4		\n\t"/*   2*t2C */\
		"addpd	%%xmm7,%%xmm7		\n\t"/*   2*t3D */\
		"addpd	%%xmm5,%%xmm5		\n\t"/*   2*t2D */\
		"addpd	%%xmm6,%%xmm6		\n\t"/*   2*t3C */\
		"movaps	%%xmm0,     (%%rbx)	\n\t"/* a(jt+p1 ) */\
		"movaps	%%xmm2,     (%%rcx)	\n\t"/* a(jt+p2 ) */\
		"movaps	%%xmm1,0x010(%%rbx)	\n\t"/* a(jp+p1 ) */\
		"movaps	%%xmm3,0x010(%%rdx)	\n\t"/* a(jp+p3 ) */\
		"addpd	%%xmm0,%%xmm4		\n\t"/* t0C+t2C */\
		"addpd	%%xmm2,%%xmm7		\n\t"/* t1C+t3D */\
		"addpd	%%xmm1,%%xmm5		\n\t"/* t0D+t2D */\
		"addpd	%%xmm3,%%xmm6		\n\t"/* t1D+t3C */\
		"movaps	%%xmm4,     (%%rax)	\n\t"/* a(jt+p0 ) */\
		"movaps	%%xmm7,     (%%rdx)	\n\t"/* a(jt+p3 ) */\
		"movaps	%%xmm5,0x010(%%rax)	\n\t"/* a(jp+p0 ) */\
		"movaps	%%xmm6,0x010(%%rcx)	\n\t"/* a(jp+p2 ) */\
		"\n\t"/*...Block 2: t02,t12,t22,t32	*/\
		"subq	$0x180,%%rsi		\n\t"/* r10 */\
		"movslq	%[__p04],%%rdi		\n\t"\
		"shlq	$3,%%rdi			\n\t"\
		"addq	%%rdi,%%rax			\n\t"/* ap10 */\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	%%rdi,%%rcx			\n\t"\
		"addq	%%rdi,%%rdx			\n\t"\
		"movq	%[__cc0],%%rdi	\n\t"\
		"\n\t"\
		"movaps	0x020(%%rsi),%%xmm4	\n\t"/* t22 */\
		"movaps	0x060(%%rsi),%%xmm6	\n\t"/* t32 */\
		"movaps	0x030(%%rsi),%%xmm5	\n\t"/* t23 */\
		"movaps	0x070(%%rsi),%%xmm7	\n\t"/* t33 */\
		"movaps	0x020(%%rsi),%%xmm0	\n\t"/* copy t22 */\
		"movaps	0x060(%%rsi),%%xmm2	\n\t"/* copy t32 */\
		"movaps	0x030(%%rsi),%%xmm1	\n\t"/* copy t23 */\
		"movaps	0x070(%%rsi),%%xmm3	\n\t"/* copy t33 */\
		"mulpd	0x020(%%rdi),%%xmm4	\n\t"/* t22*c32_1 */\
		"mulpd	0x040(%%rdi),%%xmm6	\n\t"/* t32*c32_3 */\
		"mulpd	0x030(%%rdi),%%xmm1	\n\t"/* t23*s32_1 */\
		"mulpd	0x050(%%rdi),%%xmm3	\n\t"/* t33*s32_3 */\
		"mulpd	0x020(%%rdi),%%xmm5	\n\t"/* t23*c32_1 */\
		"mulpd	0x040(%%rdi),%%xmm7	\n\t"/* t33*c32_3 */\
		"mulpd	0x030(%%rdi),%%xmm0	\n\t"/* t22*s32_1 */\
		"mulpd	0x050(%%rdi),%%xmm2	\n\t"/* t32*s32_3 */\
		"subpd	%%xmm1,%%xmm4		\n\t"/* ~t22 */\
		"subpd	%%xmm3,%%xmm6		\n\t"/* rt */\
		"addpd	%%xmm0,%%xmm5		\n\t"/* ~t23 */\
		"addpd	%%xmm2,%%xmm7		\n\t"/* it */\
		"\n\t"\
		"movaps	0x040(%%rsi),%%xmm2	\n\t"/* t12 */\
		"movaps	0x050(%%rsi),%%xmm0	\n\t"/* t13 */\
		"movaps	0x040(%%rsi),%%xmm1	\n\t"/* copy t12 */\
		"movaps	0x050(%%rsi),%%xmm3	\n\t"/* copy t13 */\
		"subpd	%%xmm6,%%xmm4		\n\t"/*~t32=t22-rt */\
		"mulpd	     (%%rdi),%%xmm2	\n\t"/* t12*c */\
		"subpd	%%xmm7,%%xmm5		\n\t"/*~t33=t23-it */\
		"mulpd	0x010(%%rdi),%%xmm0	\n\t"/* t13*s */\
		"addpd	%%xmm6,%%xmm6		\n\t"/*      2* rt */\
		"mulpd	     (%%rdi),%%xmm3	\n\t"/* t13*c */\
		"addpd	%%xmm7,%%xmm7		\n\t"/*      2* it */\
		"mulpd	0x010(%%rdi),%%xmm1	\n\t"/* t12*s */\
		"addpd	%%xmm4,%%xmm6		\n\t"/*~t22=t22+rt */\
		"subpd	%%xmm0,%%xmm2		\n\t"/* rt */\
		"addpd	%%xmm5,%%xmm7		\n\t"/*~t23=t23+it */\
		"addpd	%%xmm1,%%xmm3		\n\t"/* it */\
		"\n\t"\
		"movaps	     (%%rsi),%%xmm0	\n\t"/* t02 */\
		"movaps	0x010(%%rsi),%%xmm1	\n\t"/* t03 */\
		"subpd	%%xmm2,%%xmm0		\n\t"/*~t12=t02-rt */\
		"subpd	%%xmm3,%%xmm1		\n\t"/*~t13=t03-it */\
		"addpd	     (%%rsi),%%xmm2	\n\t"/*~t02=rt+t02 */\
		"addpd	0x010(%%rsi),%%xmm3	\n\t"/*~t03=it+t03 */\
		"subpd	%%xmm6,%%xmm2		\n\t"/* t02-t22 */\
		"subpd	%%xmm5,%%xmm0		\n\t"/* t12-t33 */\
		"subpd	%%xmm7,%%xmm3		\n\t"/* t03-t23 */\
		"subpd	%%xmm4,%%xmm1		\n\t"/* t13-t32 */\
		"addpd	%%xmm6,%%xmm6		\n\t"/*   2*t22 */\
		"addpd	%%xmm5,%%xmm5		\n\t"/*   2*t33 */\
		"addpd	%%xmm7,%%xmm7		\n\t"/*   2*t23 */\
		"addpd	%%xmm4,%%xmm4		\n\t"/*   2*t32 */\
		"movaps	%%xmm2,     (%%rbx)	\n\t"/* a(jt+p1 ) */\
		"movaps	%%xmm0,     (%%rcx)	\n\t"/* a(jt+p2 ) */\
		"movaps	%%xmm3,0x010(%%rbx)	\n\t"/* a(jp+p1 ) */\
		"movaps	%%xmm1,0x010(%%rdx)	\n\t"/* a(jp+p3 ) */\
		"addpd	%%xmm2,%%xmm6		\n\t"/* t02+t22 */\
		"addpd	%%xmm0,%%xmm5		\n\t"/* t12+t33 */\
		"addpd	%%xmm3,%%xmm7		\n\t"/* t03+t23 */\
		"addpd	%%xmm1,%%xmm4		\n\t"/* t13+t32 */\
		"movaps	%%xmm6,     (%%rax)	\n\t"/* a(jt+p0 ) */\
		"movaps	%%xmm5,     (%%rdx)	\n\t"/* a(jt+p3 ) */\
		"movaps	%%xmm7,0x010(%%rax)	\n\t"/* a(jp+p0 ) */\
		"movaps	%%xmm4,0x010(%%rcx)	\n\t"/* a(jp+p2 ) */\
		"\n\t"/*...Block 6: t0A,t1A,t2A,t3A	*/\
		"addq	$0x080,%%rsi			\n\t"/* r18 */\
		"movslq	%[__p04],%%rdi		\n\t"\
		"shlq	$3,%%rdi			\n\t"\
		"addq	%%rdi,%%rax			\n\t"/* ap14 */\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	%%rdi,%%rcx			\n\t"\
		"addq	%%rdi,%%rdx			\n\t"\
		"movq	%[__cc0],%%rdi	\n\t"\
		"\n\t"\
		"movaps	0x020(%%rsi),%%xmm4	\n\t"/* t2A */\
		"movaps	0x060(%%rsi),%%xmm6	\n\t"/* t3A */\
		"movaps	0x030(%%rsi),%%xmm5	\n\t"/* t2B */\
		"movaps	0x070(%%rsi),%%xmm7	\n\t"/* t3B */\
		"movaps	0x020(%%rsi),%%xmm0	\n\t"/* copy t2A */\
		"movaps	0x060(%%rsi),%%xmm2	\n\t"/* copy t3A */\
		"movaps	0x030(%%rsi),%%xmm1	\n\t"/* copy t2B */\
		"movaps	0x070(%%rsi),%%xmm3	\n\t"/* copy t3B */\
		"mulpd	0x050(%%rdi),%%xmm4	\n\t"/* t2A*s32_3 */\
		"mulpd	0x020(%%rdi),%%xmm6	\n\t"/* t3A*c32_1 */\
		"mulpd	0x040(%%rdi),%%xmm1	\n\t"/* t2B*c32_3 */\
		"mulpd	0x030(%%rdi),%%xmm3	\n\t"/* t3B*s32_1 */\
		"mulpd	0x050(%%rdi),%%xmm5	\n\t"/* t2B*s32_3 */\
		"mulpd	0x020(%%rdi),%%xmm7	\n\t"/* t3B*c32_1 */\
		"mulpd	0x040(%%rdi),%%xmm0	\n\t"/* t2A*c32_3 */\
		"mulpd	0x030(%%rdi),%%xmm2	\n\t"/* t3A*s32_1 */\
		"subpd	%%xmm1,%%xmm4		\n\t"/* ~t2A */\
		"addpd	%%xmm3,%%xmm6		\n\t"/* rt */\
		"addpd	%%xmm0,%%xmm5		\n\t"/* ~t2B */\
		"subpd	%%xmm2,%%xmm7		\n\t"/* it */\
		"\n\t"\
		"movaps	0x040(%%rsi),%%xmm2	\n\t"/* t1A */\
		"movaps	0x050(%%rsi),%%xmm0	\n\t"/* t1B */\
		"movaps	0x040(%%rsi),%%xmm1	\n\t"/* copy t1A */\
		"movaps	0x050(%%rsi),%%xmm3	\n\t"/* copy t1B */\
		"subpd	%%xmm6,%%xmm4		\n\t"/*~t2A=t2A-rt */\
		"mulpd	0x010(%%rdi),%%xmm2	\n\t"/* t1A*s */\
		"subpd	%%xmm7,%%xmm5		\n\t"/*~t2B=t2B-it */\
		"mulpd	     (%%rdi),%%xmm0	\n\t"/* t1B*c */\
		"addpd	%%xmm6,%%xmm6		\n\t"/*      2* rt */\
		"mulpd	0x010(%%rdi),%%xmm3	\n\t"/* t1B*s */\
		"addpd	%%xmm7,%%xmm7		\n\t"/*      2* it */\
		"mulpd	     (%%rdi),%%xmm1	\n\t"/* t1A*c */\
		"addpd	%%xmm4,%%xmm6		\n\t"/*~t3A=t2A+rt */\
		"addpd	%%xmm0,%%xmm2		\n\t"/* rt */\
		"addpd	%%xmm5,%%xmm7		\n\t"/*~t3B=t2B+it */\
		"subpd	%%xmm1,%%xmm3		\n\t"/* it */\
		"\n\t"\
		"movaps	     (%%rsi),%%xmm0	\n\t"/* t0A */\
		"movaps	0x010(%%rsi),%%xmm1	\n\t"/* t0B */\
		"subpd	%%xmm2,%%xmm0		\n\t"/*~t0A=t0A-rt */\
		"subpd	%%xmm3,%%xmm1		\n\t"/*~t0B=t0B-it */\
		"addpd	     (%%rsi),%%xmm2	\n\t"/*~t1A=rt+t0A */\
		"addpd	0x010(%%rsi),%%xmm3	\n\t"/*~t1B=it+t0B */\
		"subpd	%%xmm4,%%xmm0		\n\t"/* t0A-t2A */\
		"subpd	%%xmm7,%%xmm2		\n\t"/* t1A-t3B */\
		"subpd	%%xmm5,%%xmm1		\n\t"/* t0B-t2B */\
		"subpd	%%xmm6,%%xmm3		\n\t"/* t1B-t3A */\
		"addpd	%%xmm4,%%xmm4		\n\t"/*   2*t2A */\
		"addpd	%%xmm7,%%xmm7		\n\t"/*   2*t3B */\
		"addpd	%%xmm5,%%xmm5		\n\t"/*   2*t2B */\
		"addpd	%%xmm6,%%xmm6		\n\t"/*   2*t3A */\
		"movaps	%%xmm0,     (%%rbx)	\n\t"/* a(jt+p1 ) */\
		"movaps	%%xmm2,     (%%rcx)	\n\t"/* a(jt+p2 ) */\
		"movaps	%%xmm1,0x010(%%rbx)	\n\t"/* a(jp+p1 ) */\
		"movaps	%%xmm3,0x010(%%rdx)	\n\t"/* a(jp+p3 ) */\
		"addpd	%%xmm0,%%xmm4		\n\t"/* t0A+t2A */\
		"addpd	%%xmm2,%%xmm7		\n\t"/* t1A+t3B */\
		"addpd	%%xmm1,%%xmm5		\n\t"/* t0B+t2B */\
		"addpd	%%xmm3,%%xmm6		\n\t"/* t1B+t3A */\
		"movaps	%%xmm4,     (%%rax)	\n\t"/* a(jt+p0 ) */\
		"movaps	%%xmm7,     (%%rdx)	\n\t"/* a(jt+p3 ) */\
		"movaps	%%xmm5,0x010(%%rax)	\n\t"/* a(jp+p0 ) */\
		"movaps	%%xmm6,0x010(%%rcx)	\n\t"/* a(jp+p2 ) */\
		"\n\t"/*...Block 4: t06,t16,t26,t36	*/\
		"addq	$0x180,%%rsi		\n\t"/* r30 */\
		"movslq	%[__p04],%%rdi		\n\t"\
		"shlq	$3,%%rdi			\n\t"\
		"addq	%%rdi,%%rax			\n\t"/* ap18 */\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	%%rdi,%%rcx			\n\t"\
		"addq	%%rdi,%%rdx			\n\t"\
		"movq	%[__cc0],%%rdi	\n\t"\
		"\n\t"\
		"movaps	0x020(%%rsi),%%xmm4	\n\t"/* t26 */\
		"movaps	0x060(%%rsi),%%xmm6	\n\t"/* t36 */\
		"movaps	0x030(%%rsi),%%xmm5	\n\t"/* t27 */\
		"movaps	0x070(%%rsi),%%xmm7	\n\t"/* t37 */\
		"movaps	0x020(%%rsi),%%xmm0	\n\t"/* copy t26 */\
		"movaps	0x060(%%rsi),%%xmm2	\n\t"/* copy t36 */\
		"movaps	0x030(%%rsi),%%xmm1	\n\t"/* copy t27 */\
		"movaps	0x070(%%rsi),%%xmm3	\n\t"/* copy t37 */\
		"mulpd	0x040(%%rdi),%%xmm4	\n\t"/* t26*s32_3 */\
		"mulpd	0x030(%%rdi),%%xmm6	\n\t"/* t36*s32_1 */\
		"mulpd	0x050(%%rdi),%%xmm1	\n\t"/* t27*s32_3 */\
		"mulpd	0x020(%%rdi),%%xmm3	\n\t"/* t37*c32_1 */\
		"mulpd	0x040(%%rdi),%%xmm5	\n\t"/* t27*c32_3 */\
		"mulpd	0x030(%%rdi),%%xmm7	\n\t"/* t37*s32_1 */\
		"mulpd	0x050(%%rdi),%%xmm0	\n\t"/* t26*s32_3 */\
		"mulpd	0x020(%%rdi),%%xmm2	\n\t"/* t36*c32_1 */\
		"subpd	%%xmm1,%%xmm4		\n\t"/* ~t26 */\
		"addpd	%%xmm3,%%xmm6		\n\t"/* rt */\
		"addpd	%%xmm0,%%xmm5		\n\t"/* ~t27 */\
		"subpd	%%xmm2,%%xmm7		\n\t"/* it */\
		"\n\t"\
		"movaps	0x040(%%rsi),%%xmm2	\n\t"/* t16 */\
		"movaps	0x050(%%rsi),%%xmm0	\n\t"/* t17 */\
		"movaps	0x040(%%rsi),%%xmm1	\n\t"/* copy t16 */\
		"movaps	0x050(%%rsi),%%xmm3	\n\t"/* copy t17 */\
		"subpd	%%xmm6,%%xmm4		\n\t"/*~t26=t26-rt */\
		"mulpd	0x010(%%rdi),%%xmm2	\n\t"/* t16*s */\
		"subpd	%%xmm7,%%xmm5		\n\t"/*~t27=t27-it */\
		"mulpd	     (%%rdi),%%xmm0	\n\t"/* t17*c */\
		"addpd	%%xmm6,%%xmm6		\n\t"/*      2* rt */\
		"mulpd	0x010(%%rdi),%%xmm3	\n\t"/* t17*s */\
		"addpd	%%xmm7,%%xmm7		\n\t"/*      2* it */\
		"mulpd	     (%%rdi),%%xmm1	\n\t"/* t16*c */\
		"addpd	%%xmm4,%%xmm6		\n\t"/*~t36=t26+rt */\
		"subpd	%%xmm0,%%xmm2		\n\t"/* rt */\
		"addpd	%%xmm5,%%xmm7		\n\t"/*~t37=t27+it */\
		"addpd	%%xmm1,%%xmm3		\n\t"/* it */\
		"\n\t"\
		"movaps	     (%%rsi),%%xmm0	\n\t"/* t06 */\
		"movaps	0x010(%%rsi),%%xmm1	\n\t"/* t07 */\
		"subpd	%%xmm2,%%xmm0		\n\t"/*~t16=t06-rt */\
		"subpd	%%xmm3,%%xmm1		\n\t"/*~t17=t07-it */\
		"addpd	     (%%rsi),%%xmm2	\n\t"/*~t06=rt+t06 */\
		"addpd	0x010(%%rsi),%%xmm3	\n\t"/*~t07=it+t07 */\
		"subpd	%%xmm4,%%xmm2		\n\t"/* t06-t26 */\
		"subpd	%%xmm7,%%xmm0		\n\t"/* t16-t37 */\
		"subpd	%%xmm5,%%xmm3		\n\t"/* t07-t27 */\
		"subpd	%%xmm6,%%xmm1		\n\t"/* t17-t36 */\
		"addpd	%%xmm4,%%xmm4		\n\t"/*   2*t26 */\
		"addpd	%%xmm7,%%xmm7		\n\t"/*   2*t37 */\
		"addpd	%%xmm5,%%xmm5		\n\t"/*   2*t27 */\
		"addpd	%%xmm6,%%xmm6		\n\t"/*   2*t36 */\
		"movaps	%%xmm2,     (%%rbx)	\n\t"/* a(jt+p1 ) */\
		"movaps	%%xmm0,     (%%rcx)	\n\t"/* a(jt+p2 ) */\
		"movaps	%%xmm3,0x010(%%rbx)	\n\t"/* a(jp+p1 ) */\
		"movaps	%%xmm1,0x010(%%rdx)	\n\t"/* a(jp+p3 ) */\
		"addpd	%%xmm2,%%xmm4		\n\t"/* t06+t26 */\
		"addpd	%%xmm0,%%xmm7		\n\t"/* t16+t37 */\
		"addpd	%%xmm3,%%xmm5		\n\t"/* t07+t27 */\
		"addpd	%%xmm1,%%xmm6		\n\t"/* t17+t36 */\
		"movaps	%%xmm4,     (%%rax)	\n\t"/* a(jt+p0 ) */\
		"movaps	%%xmm7,     (%%rdx)	\n\t"/* a(jt+p3 ) */\
		"movaps	%%xmm5,0x010(%%rax)	\n\t"/* a(jp+p0 ) */\
		"movaps	%%xmm6,0x010(%%rcx)	\n\t"/* a(jp+p2 ) */\
		"\n\t"/*...Block 8: t0E,t1E,t2E,t3E	*/\
		"addq	$0x080,%%rsi			\n\t"/* r38 */\
		"movslq	%[__p04],%%rdi		\n\t"\
		"shlq	$3,%%rdi			\n\t"\
		"addq	%%rdi,%%rax			\n\t"/* ap1C */\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	%%rdi,%%rcx			\n\t"\
		"addq	%%rdi,%%rdx			\n\t"\
		"movq	%[__cc0],%%rdi	\n\t"\
		"\n\t"\
		"movaps	0x020(%%rsi),%%xmm4	\n\t"/* t2E */\
		"movaps	0x060(%%rsi),%%xmm6	\n\t"/* t3E */\
		"movaps	0x030(%%rsi),%%xmm5	\n\t"/* t2F */\
		"movaps	0x070(%%rsi),%%xmm7	\n\t"/* t3F */\
		"movaps	0x020(%%rsi),%%xmm0	\n\t"/* copy t2E */\
		"movaps	0x060(%%rsi),%%xmm2	\n\t"/* copy t3E */\
		"movaps	0x030(%%rsi),%%xmm1	\n\t"/* copy t2F */\
		"movaps	0x070(%%rsi),%%xmm3	\n\t"/* copy t3F */\
		"mulpd	0x030(%%rdi),%%xmm4	\n\t"/* t2E*s32_1 */\
		"mulpd	0x050(%%rdi),%%xmm6	\n\t"/* t3E*c32_3 */\
		"mulpd	0x020(%%rdi),%%xmm1	\n\t"/* t2F*c32_1 */\
		"mulpd	0x040(%%rdi),%%xmm3	\n\t"/* t3F*s32_3 */\
		"mulpd	0x030(%%rdi),%%xmm5	\n\t"/* t2F*s32_1 */\
		"mulpd	0x050(%%rdi),%%xmm7	\n\t"/* t3F*c32_3 */\
		"mulpd	0x020(%%rdi),%%xmm0	\n\t"/* t2E*c32_1 */\
		"mulpd	0x040(%%rdi),%%xmm2	\n\t"/* t3E*s32_3 */\
		"subpd	%%xmm1,%%xmm4		\n\t"/* ~t2E */\
		"subpd	%%xmm3,%%xmm6		\n\t"/* rt */\
		"addpd	%%xmm0,%%xmm5		\n\t"/* ~t2F */\
		"addpd	%%xmm2,%%xmm7		\n\t"/* it */\
		"\n\t"\
		"movaps	0x040(%%rsi),%%xmm2	\n\t"/* t1E */\
		"movaps	0x050(%%rsi),%%xmm0	\n\t"/* t1F */\
		"movaps	0x040(%%rsi),%%xmm1	\n\t"/* copy t1E */\
		"movaps	0x050(%%rsi),%%xmm3	\n\t"/* copy t1F */\
		"subpd	%%xmm6,%%xmm4		\n\t"/*~t2E=t2E-rt */\
		"mulpd	     (%%rdi),%%xmm2	\n\t"/* t1E*c */\
		"subpd	%%xmm7,%%xmm5		\n\t"/*~t2F=t2F-it */\
		"mulpd	0x010(%%rdi),%%xmm0	\n\t"/* t1F*s */\
		"addpd	%%xmm6,%%xmm6		\n\t"/*      2* rt */\
		"mulpd	     (%%rdi),%%xmm3	\n\t"/* t1F*c */\
		"addpd	%%xmm7,%%xmm7		\n\t"/*      2* it */\
		"mulpd	0x010(%%rdi),%%xmm1	\n\t"/* t1E*s */\
		"addpd	%%xmm4,%%xmm6		\n\t"/*~t3E=t2E+rt */\
		"addpd	%%xmm0,%%xmm2		\n\t"/* rt */\
		"addpd	%%xmm5,%%xmm7		\n\t"/*~t3F=t2F+it */\
		"subpd	%%xmm1,%%xmm3		\n\t"/* it */\
		"\n\t"\
		"movaps	     (%%rsi),%%xmm0	\n\t"/* t0E */\
		"movaps	0x010(%%rsi),%%xmm1	\n\t"/* t0F */\
		"subpd	%%xmm2,%%xmm0		\n\t"/*~t0E=t0E-rt */\
		"subpd	%%xmm3,%%xmm1		\n\t"/*~t0F=t0F-it */\
		"addpd	     (%%rsi),%%xmm2	\n\t"/*~t1E=rt+t0E */\
		"addpd	0x010(%%rsi),%%xmm3	\n\t"/*~t1F=it+t0F */\
		"subpd	%%xmm4,%%xmm0		\n\t"/* t0E-t2E */\
		"subpd	%%xmm7,%%xmm2		\n\t"/* t1E-t3F */\
		"subpd	%%xmm5,%%xmm1		\n\t"/* t0F-t2F */\
		"subpd	%%xmm6,%%xmm3		\n\t"/* t1F-t3E */\
		"addpd	%%xmm4,%%xmm4		\n\t"/*   2*t2E */\
		"addpd	%%xmm7,%%xmm7		\n\t"/*   2*t3F */\
		"addpd	%%xmm5,%%xmm5		\n\t"/*   2*t2F */\
		"addpd	%%xmm6,%%xmm6		\n\t"/*   2*t3E */\
		"movaps	%%xmm0,     (%%rbx)	\n\t"/* a(jt+p1 ) */\
		"movaps	%%xmm2,     (%%rcx)	\n\t"/* a(jt+p2 ) */\
		"movaps	%%xmm1,0x010(%%rbx)	\n\t"/* a(jp+p1 ) */\
		"movaps	%%xmm3,0x010(%%rdx)	\n\t"/* a(jp+p3 ) */\
		"addpd	%%xmm0,%%xmm4		\n\t"/* t0E+t2E */\
		"addpd	%%xmm2,%%xmm7		\n\t"/* t1E+t3F */\
		"addpd	%%xmm1,%%xmm5		\n\t"/* t0F+t2F */\
		"addpd	%%xmm3,%%xmm6		\n\t"/* t1F+t3E */\
		"movaps	%%xmm4,     (%%rax)	\n\t"/* a(jt+p0 ) */\
		"movaps	%%xmm7,     (%%rdx)	\n\t"/* a(jt+p3 ) */\
		"movaps	%%xmm5,0x010(%%rax)	\n\t"/* a(jp+p0 ) */\
		"movaps	%%xmm6,0x010(%%rcx)	\n\t"/* a(jp+p2 ) */\
		:					/* outputs: none */\
		: [__add] "m" (Xadd)	/* All inputs from memory addresses here */\
		 ,[__p01] "m" (Xp01)\
		 ,[__p02] "m" (Xp02)\
		 ,[__p03] "m" (Xp03)\
		 ,[__p04] "m" (Xp04)\
		 ,[__p08] "m" (Xp08)\
		 ,[__p10] "m" (Xp10)\
		 ,[__p18] "m" (Xp18)\
		 ,[__r00] "m" (Xr00)\
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__cc0] "m" (Xcc0)\
		: "rax","rbx","rcx","rdx","rdi","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","memory"		/* Clobbered registers */\
	);\
	}

  #else // USE_64BIT_ASM_STYLE = True: Deeper 64-bit-ified version of the above 32-bit ASM macros, using all of ymm0-15

	#error USE_64BIT_ASM_STYLE option currently unavailable for this include file!

  #endif	// USE_64BIT_ASM_STYLE ?

#endif	// AVX / SSE2 toggle

#endif	// radix32_ditN_cy_dif1_gcc_h_included

