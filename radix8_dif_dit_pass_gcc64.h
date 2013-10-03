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
#ifndef radix8_dif_dit_pass_gcc_h_included
#define radix8_dif_dit_pass_gcc_h_included

#ifdef USE_AVX	// AVX and AVX2 both use 256-bit registers

  #define USE_64BIT_ASM_STYLE	1	// My x86 timings indicate fancier versions below using ymm0-15 for the radix-16 DFT is faster.

  #if !USE_64BIT_ASM_STYLE	// False: Use simple 64-bit-ified version of the analogous 32-bit ASM macros, i.e. using just ymm0-7.

	#define SSE2_RADIX8_DIF_TWIDDLE(Xadd0,Xadd1,Xadd2,Xadd3,Xadd4,Xadd5,Xadd6,Xadd7,Xisrt2,Xc1,Xc2,Xc3,Xc4,Xc5,Xc6,Xc7)\
	{\
	__asm__ volatile (\
		"movq		%[add0]	,%%rax			\n\t"\
		"movq		%[add4]	,%%rbx			\n\t"\
		"movq		%[c4]	,%%rcx			\n\t"\
		"vmovaps	    (%%rax)	,%%ymm0		\n\t"\
		"vmovaps	0x20(%%rax)	,%%ymm1		\n\t"\
		"vmovaps	    (%%rax)	,%%ymm6		\n\t"\
		"vmovaps	0x20(%%rax)	,%%ymm7		\n\t"\
		"vmovaps	    (%%rbx)	,%%ymm2		\n\t"\
		"vmovaps	0x20(%%rbx)	,%%ymm3		\n\t"\
		"vmovaps	    (%%rbx)	,%%ymm4		\n\t"\
		"vmovaps	0x20(%%rbx)	,%%ymm5		\n\t"\
		"vmulpd			(%%rcx)	,%%ymm2,%%ymm2		\n\t"\
		"vmulpd			(%%rcx)	,%%ymm3,%%ymm3		\n\t"\
		"vmulpd		0x20(%%rcx)	,%%ymm4,%%ymm4		\n\t"\
		"vmulpd		0x20(%%rcx)	,%%ymm5,%%ymm5		\n\t"\
		"vsubpd		%%ymm5		,%%ymm2,%%ymm2		\n\t"\
		"vaddpd		%%ymm4		,%%ymm3,%%ymm3		\n\t"\
		"vaddpd		%%ymm2		,%%ymm0,%%ymm0		\n\t"\
		"vaddpd		%%ymm3		,%%ymm1,%%ymm1		\n\t"\
		"vsubpd		%%ymm2		,%%ymm6,%%ymm6		\n\t"\
		"vsubpd		%%ymm3		,%%ymm7,%%ymm7		\n\t"\
		"vmovaps	%%ymm0		,    (%%rax)	\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rax)	\n\t"\
		"vmovaps	%%ymm6		,    (%%rbx)	\n\t"\
		"vmovaps	%%ymm7		,0x20(%%rbx)	\n\t"\
		"movq		%[add2]	,%%rax			\n\t"\
		"movq		%[add6]	,%%rbx			\n\t"\
		"movq		%[c2]	,%%rcx			\n\t"\
		"movq		%[c6]	,%%rdx			\n\t"\
		"vmovaps	    (%%rax)	,%%ymm0		\n\t"\
		"vmovaps	0x20(%%rax)	,%%ymm2		\n\t"\
		"vmovaps	    (%%rax)	,%%ymm1		\n\t"\
		"vmovaps	0x20(%%rax)	,%%ymm3		\n\t"\
		"vmulpd			(%%rcx)	,%%ymm0,%%ymm0		\n\t"\
		"vmulpd		0x20(%%rcx)	,%%ymm2,%%ymm2		\n\t"\
		"vmulpd		0x20(%%rcx)	,%%ymm1,%%ymm1		\n\t"\
		"vmulpd			(%%rcx)	,%%ymm3,%%ymm3		\n\t"\
		"vsubpd		%%ymm2		,%%ymm0,%%ymm0		\n\t"\
		"vaddpd		%%ymm3		,%%ymm1,%%ymm1		\n\t"\
		"vmovaps	    (%%rbx)	,%%ymm2		\n\t"\
		"vmovaps	0x20(%%rbx)	,%%ymm3		\n\t"\
		"vmovaps	    (%%rbx)	,%%ymm4		\n\t"\
		"vmovaps	0x20(%%rbx)	,%%ymm5		\n\t"\
		"vmulpd			(%%rdx)	,%%ymm2,%%ymm2		\n\t"\
		"vmulpd		0x20(%%rdx)	,%%ymm3,%%ymm3		\n\t"\
		"vmulpd		0x20(%%rdx)	,%%ymm4,%%ymm4		\n\t"\
		"vmulpd			(%%rdx)	,%%ymm5,%%ymm5		\n\t"\
		"vsubpd		%%ymm3		,%%ymm2,%%ymm2		\n\t"\
		"vaddpd		%%ymm5		,%%ymm4,%%ymm4		\n\t"\
		"vmovaps	%%ymm2		,%%ymm3		\n\t"\
		"vmovaps	%%ymm4		,%%ymm5		\n\t"\
		"vaddpd		%%ymm0		,%%ymm2,%%ymm2		\n\t"\
		"vsubpd		%%ymm3		,%%ymm0,%%ymm0		\n\t"\
		"vaddpd		%%ymm1		,%%ymm4,%%ymm4		\n\t"\
		"vsubpd		%%ymm5		,%%ymm1,%%ymm1		\n\t"\
		"vmovaps	%%ymm2		,    (%%rax)	\n\t"\
		"vmovaps	%%ymm4		,0x20(%%rax)	\n\t"\
		"vmovaps	%%ymm0		,    (%%rbx)	\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rbx)	\n\t"\
		"movq		%[add1]	,%%rax			\n\t"\
		"movq		%[add5]	,%%rbx			\n\t"\
		"movq		%[c1]	,%%rcx			\n\t"\
		"movq		%[c5]	,%%rdx			\n\t"\
		"vmovaps	    (%%rax)	,%%ymm0		\n\t"\
		"vmovaps	0x20(%%rax)	,%%ymm2		\n\t"\
		"vmovaps	    (%%rax)	,%%ymm1		\n\t"\
		"vmovaps	0x20(%%rax)	,%%ymm3		\n\t"\
		"vmulpd			(%%rcx)	,%%ymm0,%%ymm0		\n\t"\
		"vmulpd		0x20(%%rcx)	,%%ymm2,%%ymm2		\n\t"\
		"vmulpd		0x20(%%rcx)	,%%ymm1,%%ymm1		\n\t"\
		"vmulpd			(%%rcx)	,%%ymm3,%%ymm3		\n\t"\
		"vsubpd		%%ymm2		,%%ymm0,%%ymm0		\n\t"\
		"vaddpd		%%ymm3		,%%ymm1,%%ymm1		\n\t"\
		"vmovaps	    (%%rbx)	,%%ymm2		\n\t"\
		"vmovaps	0x20(%%rbx)	,%%ymm3		\n\t"\
		"vmovaps	    (%%rbx)	,%%ymm4		\n\t"\
		"vmovaps	0x20(%%rbx)	,%%ymm5		\n\t"\
		"vmulpd			(%%rdx)	,%%ymm2,%%ymm2		\n\t"\
		"vmulpd		0x20(%%rdx)	,%%ymm3,%%ymm3		\n\t"\
		"vmulpd		0x20(%%rdx)	,%%ymm4,%%ymm4		\n\t"\
		"vmulpd			(%%rdx)	,%%ymm5,%%ymm5		\n\t"\
		"vsubpd		%%ymm3		,%%ymm2,%%ymm2		\n\t"\
		"vaddpd		%%ymm5		,%%ymm4,%%ymm4		\n\t"\
		"vmovaps	%%ymm2		,%%ymm3		\n\t"\
		"vmovaps	%%ymm4		,%%ymm5		\n\t"\
		"vaddpd		%%ymm0		,%%ymm2,%%ymm2		\n\t"\
		"vsubpd		%%ymm3		,%%ymm0,%%ymm0		\n\t"\
		"vaddpd		%%ymm1		,%%ymm4,%%ymm4		\n\t"\
		"vsubpd		%%ymm5		,%%ymm1,%%ymm1		\n\t"\
		"vmovaps	%%ymm2		,    (%%rax)	\n\t"\
		"vmovaps	%%ymm4		,0x20(%%rax)	\n\t"\
		"vmovaps	%%ymm0		,    (%%rbx)	\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rbx)	\n\t"\
		"movq		%[add3]	,%%rax			\n\t"\
		"movq		%[add7]	,%%rbx			\n\t"\
		"movq		%[c3]	,%%rcx			\n\t"\
		"movq		%[c7]	,%%rdx			\n\t"\
		"vmovaps	    (%%rax)	,%%ymm0		\n\t"\
		"vmovaps	0x20(%%rax)	,%%ymm2		\n\t"\
		"vmovaps	    (%%rax)	,%%ymm1		\n\t"\
		"vmovaps	0x20(%%rax)	,%%ymm3		\n\t"\
		"vmulpd			(%%rcx)	,%%ymm0,%%ymm0		\n\t"\
		"vmulpd		0x20(%%rcx)	,%%ymm2,%%ymm2		\n\t"\
		"vmulpd		0x20(%%rcx)	,%%ymm1,%%ymm1		\n\t"\
		"vmulpd			(%%rcx)	,%%ymm3,%%ymm3		\n\t"\
		"vsubpd		%%ymm2		,%%ymm0,%%ymm0		\n\t"\
		"vaddpd		%%ymm3		,%%ymm1,%%ymm1		\n\t"\
		"vmovaps	    (%%rbx)	,%%ymm2		\n\t"\
		"vmovaps	0x20(%%rbx)	,%%ymm3		\n\t"\
		"vmovaps	    (%%rbx)	,%%ymm4		\n\t"\
		"vmovaps	0x20(%%rbx)	,%%ymm5		\n\t"\
		"vmulpd			(%%rdx)	,%%ymm2,%%ymm2		\n\t"\
		"vmulpd		0x20(%%rdx)	,%%ymm3,%%ymm3		\n\t"\
		"vmulpd		0x20(%%rdx)	,%%ymm4,%%ymm4		\n\t"\
		"vmulpd			(%%rdx)	,%%ymm5,%%ymm5		\n\t"\
		"vsubpd		%%ymm3		,%%ymm2,%%ymm2		\n\t"\
		"vaddpd		%%ymm5		,%%ymm4,%%ymm4		\n\t"\
		"vmovaps	%%ymm2		,%%ymm3		\n\t"\
		"vmovaps	%%ymm4		,%%ymm5		\n\t"\
		"vaddpd		%%ymm0		,%%ymm2,%%ymm2		\n\t"\
		"vsubpd		%%ymm3		,%%ymm0,%%ymm0		\n\t"\
		"vaddpd		%%ymm1		,%%ymm4,%%ymm4		\n\t"\
		"vsubpd		%%ymm5		,%%ymm1,%%ymm1		\n\t"\
		"vmovaps	%%ymm2		,    (%%rax)	\n\t"\
		"vmovaps	%%ymm4		,0x20(%%rax)	\n\t"\
		"vmovaps	%%ymm0		,    (%%rbx)	\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rbx)	\n\t"\
		"movq		%[add0]	,%%rax			\n\t"\
		"movq		%[add2]	,%%rbx			\n\t"\
		"vmovaps	    (%%rax)	,%%ymm0		\n\t"\
		"vmovaps	0x20(%%rax)	,%%ymm1		\n\t"\
		"vmovaps	%%ymm0		,%%ymm4		\n\t"\
		"vmovaps	%%ymm1		,%%ymm5		\n\t"\
		"vaddpd			(%%rbx)	,%%ymm0,%%ymm0		\n\t"\
		"vsubpd			(%%rbx)	,%%ymm4,%%ymm4		\n\t"\
		"vaddpd		0x20(%%rbx)	,%%ymm1,%%ymm1		\n\t"\
		"vsubpd		0x20(%%rbx)	,%%ymm5,%%ymm5		\n\t"\
		"vmovaps	%%ymm0		,    (%%rax)	\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rax)	\n\t"\
		"vmovaps	%%ymm4		,    (%%rbx)	\n\t"\
		"vmovaps	%%ymm5		,0x20(%%rbx)	\n\t"\
		"movq		%[add1]	,%%rcx			\n\t"\
		"movq		%[add3]	,%%rdx			\n\t"\
		"vmovaps	    (%%rcx)	,%%ymm2		\n\t"\
		"vmovaps	0x20(%%rcx)	,%%ymm3		\n\t"\
		"vmovaps	%%ymm2		,%%ymm6		\n\t"\
		"vmovaps	%%ymm3		,%%ymm7		\n\t"\
		"vaddpd			(%%rdx)	,%%ymm2,%%ymm2		\n\t"\
		"vsubpd			(%%rdx)	,%%ymm6,%%ymm6		\n\t"\
		"vaddpd		0x20(%%rdx)	,%%ymm3,%%ymm3		\n\t"\
		"vsubpd		0x20(%%rdx)	,%%ymm7,%%ymm7		\n\t"\
		"vsubpd		%%ymm2		,%%ymm0,%%ymm0		\n\t"\
		"vsubpd		%%ymm3		,%%ymm1,%%ymm1		\n\t"\
		"vsubpd		%%ymm7		,%%ymm4,%%ymm4		\n\t"\
		"vsubpd		%%ymm6		,%%ymm5,%%ymm5		\n\t"\
		"vaddpd			(%%rax)	,%%ymm2,%%ymm2		\n\t"\
		"vaddpd		0x20(%%rax)	,%%ymm3,%%ymm3		\n\t"\
		"vaddpd			(%%rbx)	,%%ymm7,%%ymm7		\n\t"\
		"vaddpd		0x20(%%rbx)	,%%ymm6,%%ymm6		\n\t"\
		"vmovaps	%%ymm2		,    (%%rax)	\n\t"\
		"vmovaps	%%ymm3		,0x20(%%rax)	\n\t"\
		"vmovaps	%%ymm4		,    (%%rbx)	\n\t"\
		"vmovaps	%%ymm6		,0x20(%%rbx)	\n\t"\
		"vmovaps	%%ymm0		,    (%%rcx)	\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rcx)	\n\t"\
		"vmovaps	%%ymm7		,    (%%rdx)	\n\t"\
		"vmovaps	%%ymm5		,0x20(%%rdx)	\n\t"\
		"movq		%[add4]	,%%rax			\n\t"\
		"movq		%[add6]	,%%rbx			\n\t"\
		"vmovaps	    (%%rax)	,%%ymm0		\n\t"\
		"vmovaps	0x20(%%rax)	,%%ymm1		\n\t"\
		"vmovaps	%%ymm0		,%%ymm4		\n\t"\
		"vmovaps	%%ymm1		,%%ymm5		\n\t"\
		"vsubpd		0x20(%%rbx)	,%%ymm0,%%ymm0		\n\t"\
		"vaddpd		0x20(%%rbx)	,%%ymm4,%%ymm4		\n\t"\
		"vaddpd			(%%rbx)	,%%ymm1,%%ymm1		\n\t"\
		"vsubpd			(%%rbx)	,%%ymm5,%%ymm5		\n\t"\
		"vmovaps	%%ymm0		,    (%%rax)	\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rax)	\n\t"\
		"vmovaps	%%ymm4		,    (%%rbx)	\n\t"\
		"vmovaps	%%ymm5		,0x20(%%rbx)	\n\t"\
		"movq		%[add5]	,%%rcx			\n\t"\
		"movq		%[add7]	,%%rdx			\n\t"\
		"vmovaps	    (%%rcx)	,%%ymm2		\n\t"\
		"vmovaps	0x20(%%rcx)	,%%ymm3		\n\t"\
		"vmovaps	%%ymm2		,%%ymm6		\n\t"\
		"vmovaps	%%ymm3		,%%ymm7		\n\t"\
		"vsubpd		0x20(%%rdx)	,%%ymm2,%%ymm2		\n\t"\
		"vaddpd		0x20(%%rdx)	,%%ymm6,%%ymm6		\n\t"\
		"vaddpd			(%%rdx)	,%%ymm3,%%ymm3		\n\t"\
		"vsubpd			(%%rdx)	,%%ymm7,%%ymm7		\n\t"\
		"vmovaps	%%ymm4		,    (%%rdx)	\n\t"\
		"vmovaps	%%ymm5		,0x20(%%rdx)	\n\t"\
		"movq		%[isrt2]	,%%rax		\n\t"\
		"vmovaps	%%ymm2		,%%ymm5		\n\t"\
		"vsubpd		%%ymm3		,%%ymm2,%%ymm2		\n\t"\
		"vaddpd		%%ymm3		,%%ymm5,%%ymm5		\n\t"\
		"vmulpd			(%%rax)	,%%ymm2,%%ymm2		\n\t"\
		"vmulpd			(%%rax)	,%%ymm5,%%ymm5		\n\t"\
		"vmovaps	0x20(%%rdx)	,%%ymm3		\n\t"\
		"vmovaps	%%ymm7		,%%ymm4		\n\t"\
		"vaddpd		%%ymm6		,%%ymm4,%%ymm4		\n\t"\
		"vsubpd		%%ymm6		,%%ymm7,%%ymm7		\n\t"\
		"vmulpd		(%%rax)		,%%ymm4,%%ymm4		\n\t"\
		"vmulpd		(%%rax)		,%%ymm7,%%ymm7		\n\t"\
		"vmovaps		(%%rdx)	,%%ymm6		\n\t"\
		"movq		%[add4]	,%%rax			\n\t"\
		"vsubpd		%%ymm2		,%%ymm0,%%ymm0		\n\t"\
		"vsubpd		%%ymm5		,%%ymm1,%%ymm1		\n\t"\
		"vsubpd		%%ymm4		,%%ymm6,%%ymm6		\n\t"\
		"vsubpd		%%ymm7		,%%ymm3,%%ymm3		\n\t"\
		"vaddpd			(%%rax)	,%%ymm2,%%ymm2		\n\t"\
		"vaddpd		0x20(%%rax)	,%%ymm5,%%ymm5		\n\t"\
		"vaddpd			(%%rbx)	,%%ymm4,%%ymm4		\n\t"\
		"vaddpd		0x20(%%rbx)	,%%ymm7,%%ymm7		\n\t"\
		"vmovaps	%%ymm2		,    (%%rax)	\n\t"\
		"vmovaps	%%ymm5		,0x20(%%rax)	\n\t"\
		"vmovaps	%%ymm6		,    (%%rbx)	\n\t"\
		"vmovaps	%%ymm3		,0x20(%%rbx)	\n\t"\
		"vmovaps	%%ymm0		,    (%%rcx)	\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rcx)	\n\t"\
		"vmovaps	%%ymm4		,    (%%rdx)	\n\t"\
		"vmovaps	%%ymm7		,0x20(%%rdx)	\n\t"\
		:					/* outputs: none */\
		: [add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[add1] "m" (Xadd1)\
		 ,[add2] "m" (Xadd2)\
		 ,[add3] "m" (Xadd3)\
		 ,[add4] "m" (Xadd4)\
		 ,[add5] "m" (Xadd5)\
		 ,[add6] "m" (Xadd6)\
		 ,[add7] "m" (Xadd7)\
		 ,[isrt2] "m" (Xisrt2)\
		 ,[c1] "m" (Xc1)\
		 ,[c2] "m" (Xc2)\
		 ,[c3] "m" (Xc3)\
		 ,[c4] "m" (Xc4)\
		 ,[c5] "m" (Xc5)\
		 ,[c6] "m" (Xc6)\
		 ,[c7] "m" (Xc7)\
		: "cc","memory","rax","rbx","rcx","rdx","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
	);\
	}

	#define SSE2_RADIX8_DIT_TWIDDLE(Xadd0,Xadd1,Xadd2,Xadd3,Xadd4,Xadd5,Xadd6,Xadd7,Xisrt2,Xc1,Xc2,Xc3,Xc4,Xc5,Xc6,Xc7)\
	{\
	__asm__ volatile (\
		"movq		%[add4]		,%%rax		\n\t"\
		"movq		%[add5]		,%%rbx		\n\t"\
		"vmovaps		(%%rax)	,%%ymm0		\n\t"\
		"vmovaps	0x20(%%rax)	,%%ymm1		\n\t"\
		"vmovaps	%%ymm0		,%%ymm2		\n\t"\
		"vmovaps	%%ymm1		,%%ymm3		\n\t"\
		"vaddpd			(%%rbx)	,%%ymm2,%%ymm2		\n\t"\
		"vaddpd		0x20(%%rbx)	,%%ymm3,%%ymm3		\n\t"\
		"vsubpd			(%%rbx)	,%%ymm0,%%ymm0		\n\t"\
		"vsubpd		0x20(%%rbx)	,%%ymm1,%%ymm1		\n\t"\
		"movq		%[add6]		,%%rcx		\n\t"\
		"movq		%[add7]		,%%rdx		\n\t"\
		"vmovaps		(%%rcx)	,%%ymm4		\n\t"\
		"vmovaps	0x20(%%rcx)	,%%ymm5		\n\t"\
		"vmovaps	%%ymm4		,%%ymm6		\n\t"\
		"vmovaps	%%ymm5		,%%ymm7		\n\t"\
		"vaddpd			(%%rdx)	,%%ymm6,%%ymm6		\n\t"\
		"vaddpd		0x20(%%rdx)	,%%ymm7,%%ymm7		\n\t"\
		"vsubpd			(%%rdx)	,%%ymm4,%%ymm4		\n\t"\
		"vsubpd		0x20(%%rdx)	,%%ymm5,%%ymm5		\n\t"\
		"vmovaps	%%ymm6	,    (%%rcx)	\n\t"\
		"vmovaps	%%ymm7	,0x20(%%rcx)	\n\t"\
		"vmovaps	%%ymm4	,    (%%rdx)	\n\t"\
		"vmovaps	%%ymm5	,0x20(%%rdx)	\n\t"\
		"vaddpd		%%ymm2		,%%ymm6,%%ymm6		\n\t"\
		"vaddpd		%%ymm3		,%%ymm7,%%ymm7		\n\t"\
		"vsubpd			(%%rcx)	,%%ymm2,%%ymm2		\n\t"\
		"vsubpd		0x20(%%rcx)	,%%ymm3,%%ymm3		\n\t"\
		"vmovaps	%%ymm2	,    (%%rcx)	\n\t"\
		"vmovaps	%%ymm3	,0x20(%%rcx)	\n\t"\
		"vmovaps	%%ymm4	,%%ymm2			\n\t"\
		"vmovaps	%%ymm5	,%%ymm3			\n\t"\
		"vaddpd		%%ymm0		,%%ymm5,%%ymm5		\n\t"\
		"vsubpd		%%ymm3		,%%ymm0,%%ymm0		\n\t"\
		"vaddpd		%%ymm1		,%%ymm4,%%ymm4		\n\t"\
		"vsubpd		%%ymm2		,%%ymm1,%%ymm1		\n\t"\
		"vmovaps	%%ymm5		,%%ymm2		\n\t"\
		"vmovaps	%%ymm1		,%%ymm3		\n\t"\
		"vaddpd		%%ymm1		,%%ymm5,%%ymm5		\n\t"\
		"movq		%[isrt2]	,%%rcx		\n\t"\
		"vmovaps	(%%rcx)		,%%ymm1		\n\t"\
		"vsubpd		%%ymm3		,%%ymm2,%%ymm2		\n\t"\
		"vmulpd		%%ymm1		,%%ymm5,%%ymm5		\n\t"\
		"vmulpd		%%ymm1		,%%ymm2,%%ymm2		\n\t"\
		"vmovaps	%%ymm0		,%%ymm3		\n\t"\
		"vmovaps	%%ymm5		,    (%%rbx)\n\t"\
		"vmovaps	%%ymm4		,%%ymm5		\n\t"\
		"vaddpd		%%ymm4		,%%ymm0,%%ymm0		\n\t"\
		"vmovaps	%%ymm2	,0x20(%%rbx)	\n\t"\
		"vsubpd		%%ymm5		,%%ymm3,%%ymm3		\n\t"\
		"vmulpd		%%ymm1		,%%ymm0,%%ymm0		\n\t"\
		"vmulpd		%%ymm1		,%%ymm3,%%ymm3		\n\t"\
		"vmovaps	%%ymm0	,    (%%rdx)	\n\t"\
		"vmovaps	%%ymm3	,0x20(%%rdx)	\n\t"\
		"movq		%[add0]		,%%rax		\n\t"\
		"movq		%[add1]		,%%rbx		\n\t"\
		"vmovaps		(%%rax)	,%%ymm0		\n\t"\
		"vmovaps	0x20(%%rax)	,%%ymm1		\n\t"\
		"vmovaps	%%ymm0		,%%ymm2		\n\t"\
		"vmovaps	%%ymm1		,%%ymm3		\n\t"\
		"vaddpd			(%%rbx)	,%%ymm2,%%ymm2		\n\t"\
		"vaddpd		0x20(%%rbx)	,%%ymm3,%%ymm3		\n\t"\
		"vsubpd			(%%rbx)	,%%ymm0,%%ymm0		\n\t"\
		"vsubpd		0x20(%%rbx)	,%%ymm1,%%ymm1		\n\t"\
		"vmovaps	%%ymm6	,    (%%rax)	\n\t"\
		"vmovaps	%%ymm7	,0x20(%%rax)	\n\t"\
		"vaddpd		%%ymm6		,%%ymm6,%%ymm6		\n\t"\
		"vaddpd		%%ymm7		,%%ymm7,%%ymm7		\n\t"\
		"movq		%[add4]		,%%rcx		\n\t"\
		"vmovaps	%%ymm6	,    (%%rcx)	\n\t"\
		"vmovaps	%%ymm7	,0x20(%%rcx)	\n\t"\
		"movq		%[add2]		,%%rcx		\n\t"\
		"movq		%[add3]		,%%rdx		\n\t"\
		"vmovaps		(%%rcx)	,%%ymm4		\n\t"\
		"vmovaps	0x20(%%rcx)	,%%ymm5		\n\t"\
		"vmovaps	%%ymm4		,%%ymm6		\n\t"\
		"vmovaps	%%ymm5		,%%ymm7		\n\t"\
		"vaddpd			(%%rdx)	,%%ymm6,%%ymm6		\n\t"\
		"vaddpd		0x20(%%rdx)	,%%ymm7,%%ymm7		\n\t"\
		"vsubpd			(%%rdx)	,%%ymm4,%%ymm4		\n\t"\
		"vsubpd		0x20(%%rdx)	,%%ymm5,%%ymm5		\n\t"\
		"vmovaps	%%ymm6	,    (%%rcx)	\n\t"\
		"vmovaps	%%ymm7	,0x20(%%rcx)	\n\t"\
		"vaddpd		%%ymm2		,%%ymm6,%%ymm6		\n\t"\
		"vaddpd		%%ymm3		,%%ymm7,%%ymm7		\n\t"\
		"vsubpd			(%%rcx)	,%%ymm2,%%ymm2		\n\t"\
		"vsubpd		0x20(%%rcx)	,%%ymm3,%%ymm3		\n\t"\
		"vaddpd			(%%rax)	,%%ymm6,%%ymm6		\n\t"\
		"vaddpd		0x20(%%rax)	,%%ymm7,%%ymm7		\n\t"\
		"vmovaps	%%ymm6	,    (%%rax)	\n\t"\
		"vmovaps	%%ymm7	,0x20(%%rax)	\n\t"\
		"movq		%[add4]		,%%rbx		\n\t"\
		"vsubpd			(%%rbx)	,%%ymm6,%%ymm6		\n\t"\
		"vsubpd		0x20(%%rbx)	,%%ymm7,%%ymm7		\n\t"\
		"vmovaps	%%ymm6	,    (%%rbx)	\n\t"\
		"vmovaps	%%ymm7	,0x20(%%rbx)	\n\t"\
		"vmovaps	%%ymm6	,    (%%rdx)	\n\t"\
		"vmovaps	%%ymm7	,0x20(%%rdx)	\n\t"\
		"vmovaps	%%ymm4		,%%ymm6		\n\t"\
		"vmovaps	%%ymm5		,%%ymm7		\n\t"\
		"vaddpd		%%ymm0		,%%ymm5,%%ymm5		\n\t"\
		"vsubpd		%%ymm7		,%%ymm0,%%ymm0		\n\t"\
		"vaddpd		%%ymm1		,%%ymm4,%%ymm4		\n\t"\
		"vsubpd		%%ymm6		,%%ymm1,%%ymm1		\n\t"\
		"movq		%[c4]		,%%rcx		\n\t"\
		"vmovaps	    (%%rbx)	,%%ymm6		\n\t"\
		"vmovaps	0x20(%%rbx)	,%%ymm7		\n\t"\
		"vmulpd		0x20(%%rcx)	,%%ymm6,%%ymm6		\n\t"\
		"vmulpd		0x20(%%rcx)	,%%ymm7,%%ymm7		\n\t"\
		"vmovaps	%%ymm6	,    (%%rbx)	\n\t"\
		"vmovaps	%%ymm7	,0x20(%%rbx)	\n\t"\
		"vmovaps	    (%%rdx)	,%%ymm6		\n\t"\
		"vmovaps	0x20(%%rdx)	,%%ymm7		\n\t"\
		"vmulpd			(%%rcx)	,%%ymm6,%%ymm6		\n\t"\
		"vmulpd			(%%rcx)	,%%ymm7,%%ymm7		\n\t"\
		"vaddpd		0x20(%%rbx)	,%%ymm6,%%ymm6		\n\t"\
		"vsubpd			(%%rbx)	,%%ymm7,%%ymm7		\n\t"\
		"vmovaps	%%ymm6	,    (%%rbx)	\n\t"\
		"vmovaps	%%ymm7	,0x20(%%rbx)	\n\t"\
		"movq		%[add1]		,%%rax		\n\t"\
		"movq		%[add5]		,%%rbx		\n\t"\
		"movq		%[c1]		,%%rcx		\n\t"\
		"movq		%[c5]		,%%rdx		\n\t"\
		"vmovaps	%%ymm5		,%%ymm6		\n\t"\
		"vmovaps	%%ymm1		,%%ymm7		\n\t"\
		"vaddpd			(%%rbx)	,%%ymm5,%%ymm5		\n\t"\
		"vsubpd		0x20(%%rbx)	,%%ymm1,%%ymm1		\n\t"\
		"vsubpd			(%%rbx)	,%%ymm6,%%ymm6		\n\t"\
		"vaddpd		0x20(%%rbx)	,%%ymm7,%%ymm7		\n\t"\
		"vmovaps	%%ymm6	,    (%%rbx)	\n\t"\
		"vmovaps	%%ymm7	,0x20(%%rbx)	\n\t"\
		"vmovaps	%%ymm5		,%%ymm6		\n\t"\
		"vmovaps	%%ymm1		,%%ymm7		\n\t"\
		"vmulpd			(%%rcx)	,%%ymm5,%%ymm5		\n\t"\
		"vmulpd			(%%rcx)	,%%ymm1,%%ymm1		\n\t"\
		"vmulpd		0x20(%%rcx)	,%%ymm6,%%ymm6		\n\t"\
		"vmulpd		0x20(%%rcx)	,%%ymm7,%%ymm7		\n\t"\
		"vsubpd		%%ymm6		,%%ymm1,%%ymm1		\n\t"\
		"vaddpd		%%ymm7		,%%ymm5,%%ymm5		\n\t"\
		"vmovaps	%%ymm1	,0x20(%%rax)	\n\t"\
		"vmovaps	%%ymm5	,    (%%rax)	\n\t"\
		"vmovaps	    (%%rbx)	,%%ymm5		\n\t"\
		"vmovaps	0x20(%%rbx)	,%%ymm1		\n\t"\
		"vmovaps	%%ymm5		,%%ymm6		\n\t"\
		"vmovaps	%%ymm1		,%%ymm7		\n\t"\
		"vmulpd			(%%rdx)	,%%ymm5,%%ymm5		\n\t"\
		"vmulpd			(%%rdx)	,%%ymm1,%%ymm1		\n\t"\
		"vmulpd		0x20(%%rdx)	,%%ymm6,%%ymm6		\n\t"\
		"vmulpd		0x20(%%rdx)	,%%ymm7,%%ymm7		\n\t"\
		"vsubpd		%%ymm6		,%%ymm1,%%ymm1		\n\t"\
		"vaddpd		%%ymm7		,%%ymm5,%%ymm5		\n\t"\
		"vmovaps	%%ymm1	,0x20(%%rbx)	\n\t"\
		"vmovaps	%%ymm5	,    (%%rbx)	\n\t"\
		"movq		%[add2]		,%%rax		\n\t"\
		"movq		%[add6]		,%%rbx		\n\t"\
		"movq		%[c2]		,%%rcx		\n\t"\
		"movq		%[c6]		,%%rdx		\n\t"\
		"vmovaps	%%ymm2		,%%ymm6		\n\t"\
		"vmovaps	%%ymm3		,%%ymm7		\n\t"\
		"vaddpd		0x20(%%rbx)	,%%ymm2,%%ymm2		\n\t"\
		"vsubpd			(%%rbx)	,%%ymm3,%%ymm3		\n\t"\
		"vsubpd		0x20(%%rbx)	,%%ymm6,%%ymm6		\n\t"\
		"vaddpd			(%%rbx)	,%%ymm7,%%ymm7		\n\t"\
		"vmovaps	%%ymm2		,%%ymm1		\n\t"\
		"vmovaps	%%ymm3		,%%ymm5		\n\t"\
		"vmulpd			(%%rcx)	,%%ymm2,%%ymm2		\n\t"\
		"vmulpd			(%%rcx)	,%%ymm3,%%ymm3		\n\t"\
		"vmulpd		0x20(%%rcx)	,%%ymm1,%%ymm1		\n\t"\
		"vmulpd		0x20(%%rcx)	,%%ymm5,%%ymm5		\n\t"\
		"vsubpd		%%ymm1		,%%ymm3,%%ymm3		\n\t"\
		"vaddpd		%%ymm5		,%%ymm2,%%ymm2		\n\t"\
		"vmovaps	%%ymm3	,0x20(%%rax)	\n\t"\
		"vmovaps	%%ymm2	,    (%%rax)	\n\t"\
		"vmovaps	%%ymm6		,%%ymm1		\n\t"\
		"vmovaps	%%ymm7		,%%ymm5		\n\t"\
		"vmulpd			(%%rdx)	,%%ymm6,%%ymm6		\n\t"\
		"vmulpd			(%%rdx)	,%%ymm7,%%ymm7		\n\t"\
		"vmulpd		0x20(%%rdx)	,%%ymm1,%%ymm1		\n\t"\
		"vmulpd		0x20(%%rdx)	,%%ymm5,%%ymm5		\n\t"\
		"vsubpd		%%ymm1		,%%ymm7,%%ymm7		\n\t"\
		"vaddpd		%%ymm5		,%%ymm6,%%ymm6		\n\t"\
		"vmovaps	%%ymm7	,0x20(%%rbx)	\n\t"\
		"vmovaps	%%ymm6	,    (%%rbx)	\n\t"\
		"movq		%[add3]		,%%rax		\n\t"\
		"movq		%[add7]		,%%rbx		\n\t"\
		"movq		%[c3]		,%%rcx		\n\t"\
		"movq		%[c7]		,%%rdx		\n\t"\
		"vmovaps	%%ymm0		,%%ymm6		\n\t"\
		"vmovaps	%%ymm4		,%%ymm7		\n\t"\
		"vsubpd		0x20(%%rbx)	,%%ymm0,%%ymm0		\n\t"\
		"vsubpd			(%%rbx)	,%%ymm4,%%ymm4		\n\t"\
		"vaddpd		0x20(%%rbx)	,%%ymm6,%%ymm6		\n\t"\
		"vaddpd			(%%rbx)	,%%ymm7,%%ymm7		\n\t"\
		"vmovaps	%%ymm0		,%%ymm1		\n\t"\
		"vmovaps	%%ymm4		,%%ymm5		\n\t"\
		"vmulpd			(%%rcx)	,%%ymm0,%%ymm0		\n\t"\
		"vmulpd			(%%rcx)	,%%ymm4,%%ymm4		\n\t"\
		"vmulpd		0x20(%%rcx)	,%%ymm1,%%ymm1		\n\t"\
		"vmulpd		0x20(%%rcx)	,%%ymm5,%%ymm5		\n\t"\
		"vsubpd		%%ymm1		,%%ymm4,%%ymm4		\n\t"\
		"vaddpd		%%ymm5		,%%ymm0,%%ymm0		\n\t"\
		"vmovaps	%%ymm4	,0x20(%%rax)	\n\t"\
		"vmovaps	%%ymm0	,    (%%rax)	\n\t"\
		"vmovaps	%%ymm6		,%%ymm1		\n\t"\
		"vmovaps	%%ymm7		,%%ymm5		\n\t"\
		"vmulpd			(%%rdx)	,%%ymm6,%%ymm6		\n\t"\
		"vmulpd			(%%rdx)	,%%ymm7,%%ymm7		\n\t"\
		"vmulpd		0x20(%%rdx)	,%%ymm1,%%ymm1		\n\t"\
		"vmulpd		0x20(%%rdx)	,%%ymm5,%%ymm5		\n\t"\
		"vsubpd		%%ymm1		,%%ymm7,%%ymm7		\n\t"\
		"vaddpd		%%ymm5		,%%ymm6,%%ymm6		\n\t"\
		"vmovaps	%%ymm7	,0x20(%%rbx)	\n\t"\
		"vmovaps	%%ymm6	,    (%%rbx)	\n\t"\
		:					/* outputs: none */\
		: [add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[add1] "m" (Xadd1)\
		 ,[add2] "m" (Xadd2)\
		 ,[add3] "m" (Xadd3)\
		 ,[add4] "m" (Xadd4)\
		 ,[add5] "m" (Xadd5)\
		 ,[add6] "m" (Xadd6)\
		 ,[add7] "m" (Xadd7)\
		 ,[isrt2] "m" (Xisrt2)\
		 ,[c1] "m" (Xc1)\
		 ,[c2] "m" (Xc2)\
		 ,[c3] "m" (Xc3)\
		 ,[c4] "m" (Xc4)\
		 ,[c5] "m" (Xc5)\
		 ,[c6] "m" (Xc6)\
		 ,[c7] "m" (Xc7)\
		: "cc","memory","rax","rbx","rcx","rdx","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
	);\
	}

  #else // USE_64BIT_ASM_STYLE = True: Deeper 64-bit-ified version of the above 32-bit ASM macros, using all of ymm0-15

	#define SSE2_RADIX8_DIF_TWIDDLE(Xadd0,Xadd1,Xadd2,Xadd3,Xadd4,Xadd5,Xadd6,Xadd7,Xisrt2,Xc1,Xc2,Xc3,Xc4,Xc5,Xc6,Xc7)\
	{\
	__asm__ volatile (\
		"													movq		%[add1]	,%%r10				\n\t"\
		"													movq		%[add5]	,%%r11				\n\t"\
		"													movq		%[c1]	,%%r12				\n\t"\
		"													movq		%[c5]	,%%r13				\n\t"\
		"													vmovaps	    (%%r10)	,%%ymm8 			\n\t"\
		"movq		%[add0]	,%%rax				\n\t		vmovaps	0x20(%%r10)	,%%ymm10			\n\t"\
		"movq		%[add4]	,%%rbx				\n\t		vmovaps	    (%%r10)	,%%ymm9 			\n\t"\
		"movq		%[c4]	,%%rcx				\n\t		vmovaps	0x20(%%r10)	,%%ymm11			\n\t"\
		"vmovaps	    (%%rax)	,%%ymm0			\n\t		vmulpd	    (%%r12)	,%%ymm8 ,%%ymm8 	\n\t"\
		"vmovaps	0x20(%%rax)	,%%ymm1			\n\t		vmulpd	0x20(%%r12)	,%%ymm10,%%ymm10	\n\t"\
		"vmovaps	    (%%rax)	,%%ymm6			\n\t		vmulpd	0x20(%%r12)	,%%ymm9 ,%%ymm9 	\n\t"\
		"vmovaps	0x20(%%rax)	,%%ymm7			\n\t		vmulpd	    (%%r12)	,%%ymm11,%%ymm11	\n\t"\
		"vmovaps	    (%%rbx)	,%%ymm2			\n\t		vsubpd	%%ymm10		,%%ymm8 ,%%ymm8 	\n\t"\
		"vmovaps	0x20(%%rbx)	,%%ymm3			\n\t		vaddpd	%%ymm11		,%%ymm9 ,%%ymm9 	\n\t"\
		"vmovaps	    (%%rbx)	,%%ymm4			\n\t		vmovaps	    (%%r11)	,%%ymm10			\n\t"\
		"vmovaps	0x20(%%rbx)	,%%ymm5			\n\t		vmovaps	0x20(%%r11)	,%%ymm11			\n\t"\
		"vmulpd	    (%%rcx)	,%%ymm2,%%ymm2		\n\t		vmovaps	    (%%r11)	,%%ymm12			\n\t"\
		"vmulpd	    (%%rcx)	,%%ymm3,%%ymm3		\n\t		vmovaps	0x20(%%r11)	,%%ymm13			\n\t"\
		"vmulpd	0x20(%%rcx)	,%%ymm4,%%ymm4		\n\t		vmulpd	    (%%r13)	,%%ymm10,%%ymm10	\n\t"\
		"vmulpd	0x20(%%rcx)	,%%ymm5,%%ymm5		\n\t		vmulpd	0x20(%%r13)	,%%ymm11,%%ymm11	\n\t"\
		"vsubpd	%%ymm5		,%%ymm2,%%ymm2		\n\t		vmulpd	0x20(%%r13)	,%%ymm12,%%ymm12	\n\t"\
		"vaddpd	%%ymm4		,%%ymm3,%%ymm3		\n\t		vmulpd	    (%%r13)	,%%ymm13,%%ymm13	\n\t"\
		"vaddpd	%%ymm2		,%%ymm0,%%ymm0		\n\t		vsubpd	%%ymm11		,%%ymm10,%%ymm10	\n\t"\
		"vaddpd	%%ymm3		,%%ymm1,%%ymm1		\n\t		vaddpd	%%ymm13		,%%ymm12,%%ymm12	\n\t"\
		"vsubpd	%%ymm2		,%%ymm6,%%ymm6		\n\t		vmovaps	%%ymm10		,%%ymm11			\n\t"\
		"vsubpd	%%ymm3		,%%ymm7,%%ymm7		\n\t		vmovaps	%%ymm12		,%%ymm13			\n\t"\
		"vmovaps	%%ymm0		,    (%%rax)	\n\t		vaddpd	%%ymm8 		,%%ymm10,%%ymm10	\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rax)	\n\t		vsubpd	%%ymm11		,%%ymm8 ,%%ymm8 	\n\t"\
		"vmovaps	%%ymm6		,    (%%rbx)	\n\t		vaddpd	%%ymm9 		,%%ymm12,%%ymm12	\n\t"\
		"vmovaps	%%ymm7		,0x20(%%rbx)	\n\t		vsubpd	%%ymm13		,%%ymm9 ,%%ymm9 	\n\t"\
		"movq		%[add2]	,%%rax				\n\t		vmovaps	%%ymm10		,    (%%r10)	\n\t"\
		"movq		%[add6]	,%%rbx				\n\t		vmovaps	%%ymm12		,0x20(%%r10)	\n\t"\
		"movq		%[c2]	,%%rcx				\n\t		vmovaps	%%ymm8 		,    (%%r11)	\n\t"\
		"movq		%[c6]	,%%rdx				\n\t		vmovaps	%%ymm9 		,0x20(%%r11)	\n\t"\
		"vmovaps	    (%%rax)	,%%ymm0			\n\t		movq		%[add3]	,%%r10				\n\t"\
		"vmovaps	0x20(%%rax)	,%%ymm2			\n\t		movq		%[add7]	,%%r11				\n\t"\
		"vmovaps	    (%%rax)	,%%ymm1			\n\t		movq		%[c3]	,%%r12				\n\t"\
		"vmovaps	0x20(%%rax)	,%%ymm3			\n\t		movq		%[c7]	,%%r13				\n\t"\
		"vmulpd	    (%%rcx)	,%%ymm0,%%ymm0		\n\t		vmovaps	    (%%r10)	,%%ymm8 			\n\t"\
		"vmulpd	0x20(%%rcx)	,%%ymm2,%%ymm2		\n\t		vmovaps	0x20(%%r10)	,%%ymm10			\n\t"\
		"vmulpd	0x20(%%rcx)	,%%ymm1,%%ymm1		\n\t		vmovaps	    (%%r10)	,%%ymm9 			\n\t"\
		"vmulpd	    (%%rcx)	,%%ymm3,%%ymm3		\n\t		vmovaps	0x20(%%r10)	,%%ymm11			\n\t"\
		"vsubpd	%%ymm2		,%%ymm0,%%ymm0		\n\t		vmulpd	    (%%r12)	,%%ymm8 ,%%ymm8 	\n\t"\
		"vaddpd	%%ymm3		,%%ymm1,%%ymm1		\n\t		vmulpd	0x20(%%r12)	,%%ymm10,%%ymm10	\n\t"\
		"vmovaps	    (%%rbx)	,%%ymm2			\n\t		vmulpd	0x20(%%r12)	,%%ymm9 ,%%ymm9 	\n\t"\
		"vmovaps	0x20(%%rbx)	,%%ymm3			\n\t		vmulpd	    (%%r12)	,%%ymm11,%%ymm11	\n\t"\
		"vmovaps	    (%%rbx)	,%%ymm4			\n\t		vsubpd	%%ymm10		,%%ymm8 ,%%ymm8 	\n\t"\
		"vmovaps	0x20(%%rbx)	,%%ymm5			\n\t		vaddpd	%%ymm11		,%%ymm9 ,%%ymm9 	\n\t"\
		"vmulpd	    (%%rdx)	,%%ymm2,%%ymm2		\n\t		vmovaps	    (%%r11)	,%%ymm10			\n\t"\
		"vmulpd	0x20(%%rdx)	,%%ymm3,%%ymm3		\n\t		vmovaps	0x20(%%r11)	,%%ymm11			\n\t"\
		"vmulpd	0x20(%%rdx)	,%%ymm4,%%ymm4		\n\t		vmovaps	    (%%r11)	,%%ymm12			\n\t"\
		"vmulpd	    (%%rdx)	,%%ymm5,%%ymm5		\n\t		vmovaps	0x20(%%r11)	,%%ymm13			\n\t"\
		"vsubpd	%%ymm3		,%%ymm2,%%ymm2		\n\t		vmulpd	    (%%r13)	,%%ymm10,%%ymm10	\n\t"\
		"vaddpd	%%ymm5		,%%ymm4,%%ymm4		\n\t		vmulpd	0x20(%%r13)	,%%ymm11,%%ymm11	\n\t"\
		"vmovaps	%%ymm2		,%%ymm3			\n\t		vmulpd	0x20(%%r13)	,%%ymm12,%%ymm12	\n\t"\
		"vmovaps	%%ymm4		,%%ymm5			\n\t		vmulpd	    (%%r13)	,%%ymm13,%%ymm13	\n\t"\
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
		"movq		%[add0]	,%%rax				\n\t		movq		%[add4]	,%%r10				\n\t"\
		"movq		%[add2]	,%%rbx				\n\t		movq		%[add6]	,%%r11				\n\t"\
		"vmovaps	    (%%rax)	,%%ymm0			\n\t		vmovaps	    (%%r10)	,%%ymm8 			\n\t"\
		"vmovaps	0x20(%%rax)	,%%ymm1			\n\t		vmovaps	0x20(%%r10)	,%%ymm9 			\n\t"\
		"vmovaps	%%ymm0		,%%ymm4			\n\t		vmovaps	%%ymm8 		,%%ymm12			\n\t"\
		"vmovaps	%%ymm1		,%%ymm5			\n\t		vmovaps	%%ymm9 		,%%ymm13			\n\t"\
		"vaddpd	    (%%rbx)	,%%ymm0,%%ymm0		\n\t		vsubpd	0x20(%%r11)	,%%ymm8 ,%%ymm8 	\n\t"\
		"vsubpd	    (%%rbx)	,%%ymm4,%%ymm4		\n\t		vaddpd	0x20(%%r11)	,%%ymm12,%%ymm12	\n\t"\
		"vaddpd	0x20(%%rbx)	,%%ymm1,%%ymm1		\n\t		vaddpd	    (%%r11)	,%%ymm9 ,%%ymm9 	\n\t"\
		"vsubpd	0x20(%%rbx)	,%%ymm5,%%ymm5		\n\t		vsubpd	    (%%r11)	,%%ymm13,%%ymm13	\n\t"\
		"vmovaps	%%ymm0		,    (%%rax)	\n\t		vmovaps	%%ymm8 		,    (%%r10)	\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rax)	\n\t		vmovaps	%%ymm9 		,0x20(%%r10)	\n\t"\
		"vmovaps	%%ymm4		,    (%%rbx)	\n\t		vmovaps	%%ymm12		,    (%%r11)	\n\t"\
		"vmovaps	%%ymm5		,0x20(%%rbx)	\n\t		vmovaps	%%ymm13		,0x20(%%r11)	\n\t"\
		"movq		%[add1]	,%%rcx				\n\t		movq		%[add5]	,%%r12				\n\t"\
		"movq		%[add3]	,%%rdx				\n\t		movq		%[add7]	,%%r13				\n\t"\
		"vmovaps	    (%%rcx)	,%%ymm2			\n\t		vmovaps	    (%%r12)	,%%ymm10			\n\t"\
		"vmovaps	0x20(%%rcx)	,%%ymm3			\n\t		vmovaps	0x20(%%r12)	,%%ymm11			\n\t"\
		"vmovaps	%%ymm2		,%%ymm6			\n\t		vmovaps	%%ymm10		,%%ymm14			\n\t"\
		"vmovaps	%%ymm3		,%%ymm7			\n\t		vmovaps	%%ymm11		,%%ymm15			\n\t"\
		"vaddpd	    (%%rdx)	,%%ymm2,%%ymm2		\n\t		vsubpd	0x20(%%r13)	,%%ymm10,%%ymm10	\n\t"\
		"vsubpd	    (%%rdx)	,%%ymm6,%%ymm6		\n\t		vaddpd	0x20(%%r13)	,%%ymm14,%%ymm14	\n\t"\
		"vaddpd	0x20(%%rdx)	,%%ymm3,%%ymm3		\n\t		vaddpd	    (%%r13)	,%%ymm11,%%ymm11	\n\t"\
		"vsubpd	0x20(%%rdx)	,%%ymm7,%%ymm7		\n\t		vsubpd	    (%%r13)	,%%ymm15,%%ymm15	\n\t"\
		"vsubpd	%%ymm2		,%%ymm0,%%ymm0		\n\t		vmovaps	%%ymm12		,    (%%r13)	\n\t"\
		"vsubpd	%%ymm3		,%%ymm1,%%ymm1		\n\t		vmovaps	%%ymm13		,0x20(%%r13)	\n\t"\
		"vsubpd	%%ymm7		,%%ymm4,%%ymm4		\n\t		movq		%[isrt2]	,%%r10			\n\t"\
		"vsubpd	%%ymm6		,%%ymm5,%%ymm5		\n\t		vmovaps	%%ymm10		,%%ymm13			\n\t"\
		"vaddpd	    (%%rax)	,%%ymm2,%%ymm2		\n\t		vsubpd	%%ymm11		,%%ymm10,%%ymm10	\n\t"\
		"vaddpd	0x20(%%rax)	,%%ymm3,%%ymm3		\n\t		vaddpd	%%ymm11		,%%ymm13,%%ymm13	\n\t"\
		"vaddpd	    (%%rbx)	,%%ymm7,%%ymm7		\n\t		vmulpd	    (%%r10)	,%%ymm10,%%ymm10	\n\t"\
		"vaddpd	0x20(%%rbx)	,%%ymm6,%%ymm6		\n\t		vmulpd	    (%%r10)	,%%ymm13,%%ymm13	\n\t"\
		"vmovaps	%%ymm2		,    (%%rax)	\n\t		vmovaps	0x20(%%r13)	,%%ymm11			\n\t"\
		"vmovaps	%%ymm3		,0x20(%%rax)	\n\t		vmovaps	%%ymm15		,%%ymm12			\n\t"\
		"vmovaps	%%ymm4		,    (%%rbx)	\n\t		vaddpd	%%ymm14		,%%ymm12,%%ymm12	\n\t"\
		"vmovaps	%%ymm6		,0x20(%%rbx)	\n\t		vsubpd	%%ymm14		,%%ymm15,%%ymm15	\n\t"\
		"vmovaps	%%ymm0		,    (%%rcx)	\n\t		vmulpd	(%%r10)		,%%ymm12,%%ymm12	\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rcx)	\n\t		vmulpd	(%%r10)		,%%ymm15,%%ymm15	\n\t"\
		"vmovaps	%%ymm7		,    (%%rdx)	\n\t		vmovaps		(%%r13)	,%%ymm14			\n\t"\
		"vmovaps	%%ymm5		,0x20(%%rdx)	\n\t		movq		%[add4]	,%%r10				\n\t"\
		"													vsubpd	%%ymm10		,%%ymm8 ,%%ymm8 	\n\t"\
		"													vsubpd	%%ymm13		,%%ymm9 ,%%ymm9 	\n\t"\
		"													vsubpd	%%ymm12		,%%ymm14,%%ymm14	\n\t"\
		"													vsubpd	%%ymm15		,%%ymm11,%%ymm11	\n\t"\
		"													vaddpd	    (%%r10)	,%%ymm10,%%ymm10	\n\t"\
		"													vaddpd	0x20(%%r10)	,%%ymm13,%%ymm13	\n\t"\
		"													vaddpd	    (%%r11)	,%%ymm12,%%ymm12	\n\t"\
		"													vaddpd	0x20(%%r11)	,%%ymm15,%%ymm15	\n\t"\
		"													vmovaps	%%ymm10		,    (%%r10)	\n\t"\
		"													vmovaps	%%ymm13		,0x20(%%r10)	\n\t"\
		"													vmovaps	%%ymm14		,    (%%r11)	\n\t"\
		"													vmovaps	%%ymm11		,0x20(%%r11)	\n\t"\
		"													vmovaps	%%ymm8 		,    (%%r12)	\n\t"\
		"													vmovaps	%%ymm9 		,0x20(%%r12)	\n\t"\
		"													vmovaps	%%ymm12		,    (%%r13)	\n\t"\
		"													vmovaps	%%ymm15		,0x20(%%r13)	\n\t"\
		:					/* outputs: none */\
		: [add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[add1] "m" (Xadd1)\
		 ,[add2] "m" (Xadd2)\
		 ,[add3] "m" (Xadd3)\
		 ,[add4] "m" (Xadd4)\
		 ,[add5] "m" (Xadd5)\
		 ,[add6] "m" (Xadd6)\
		 ,[add7] "m" (Xadd7)\
		 ,[isrt2] "m" (Xisrt2)\
		 ,[c1] "m" (Xc1)\
		 ,[c2] "m" (Xc2)\
		 ,[c3] "m" (Xc3)\
		 ,[c4] "m" (Xc4)\
		 ,[c5] "m" (Xc5)\
		 ,[c6] "m" (Xc6)\
		 ,[c7] "m" (Xc7)\
		: "cc","memory","rax","rbx","rcx","rdx","r10","r11","r12","r13","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
	);\
	}


	#define SSE2_RADIX8_DIT_TWIDDLE(Xadd0,Xadd1,Xadd2,Xadd3,Xadd4,Xadd5,Xadd6,Xadd7,Xisrt2,Xc1,Xc2,Xc3,Xc4,Xc5,Xc6,Xc7)\
	{\
	__asm__ volatile (\
		/*** 2nd of 2 length-4 subtransforms gets done first: ***/\
		/* SSE2_RADIX4_DIT_0TWIDDLE_2NDOFTWO(add4-7):	SSE2_RADIX4_DIT_0TWIDDLE(add0-3): */\
		"movq		%[add4]		,%%rax			\n\t		movq		%[add0]		,%%r10		\n\t"\
		"movq		%[add5]		,%%rbx			\n\t		movq		%[add1]		,%%r11		\n\t"\
		"vmovaps		(%%rax)	,%%ymm0			\n\t		vmovaps			(%%r10)	,%%ymm8 	\n\t"\
		"vmovaps	0x20(%%rax)	,%%ymm1			\n\t		vmovaps		0x20(%%r10)	,%%ymm9 	\n\t"\
		"vmovaps	%%ymm0		,%%ymm2			\n\t		vmovaps		%%ymm8 		,%%ymm10	\n\t"\
		"vmovaps	%%ymm1		,%%ymm3			\n\t		vmovaps		%%ymm9 		,%%ymm11	\n\t"\
		"vaddpd		(%%rbx)	,%%ymm2,%%ymm2		\n\t		vaddpd		(%%r11)	,%%ymm10,%%ymm10\n\t"\
		"vaddpd	0x20(%%rbx)	,%%ymm3,%%ymm3		\n\t		vaddpd	0x20(%%r11)	,%%ymm11,%%ymm11\n\t"\
		"vsubpd		(%%rbx)	,%%ymm0,%%ymm0		\n\t		vsubpd		(%%r11)	,%%ymm8 ,%%ymm8 \n\t"\
		"vsubpd	0x20(%%rbx)	,%%ymm1,%%ymm1		\n\t		vsubpd	0x20(%%r11)	,%%ymm9 ,%%ymm9 \n\t"\
		"movq		%[add6]		,%%rcx			\n\t		movq		%[add2]		,%%r12		\n\t"\
		"movq		%[add7]		,%%rdx			\n\t		movq		%[add3]		,%%r13		\n\t"\
		"vmovaps			(%%rcx)	,%%ymm4		\n\t		vmovaps			(%%r12)	,%%ymm12	\n\t"\
		"vmovaps		0x20(%%rcx)	,%%ymm5		\n\t		vmovaps		0x20(%%r12)	,%%ymm13	\n\t"/* t9  in ymm6, needed below in RHS! */\
		"vmovaps		%%ymm4		,%%ymm6		\n\t		vmovaps		%%ymm12		,%%ymm14	\n\t"/* t10 in ymm7, needed below in RHS! */\
		"vmovaps		%%ymm5		,%%ymm7		\n\t		vmovaps		%%ymm13		,%%ymm15	\n\t"\
		"vaddpd		(%%rdx)	,%%ymm6,%%ymm6		\n\t		vaddpd		(%%r13)	,%%ymm14,%%ymm14\n\t"\
		"vaddpd	0x20(%%rdx)	,%%ymm7,%%ymm7		\n\t		vaddpd	0x20(%%r13)	,%%ymm15,%%ymm15\n\t"\
		"vsubpd		(%%rdx)	,%%ymm4,%%ymm4		\n\t		vsubpd		(%%r13)	,%%ymm12,%%ymm12\n\t"\
		"vsubpd	0x20(%%rdx)	,%%ymm5,%%ymm5		\n\t		vsubpd	0x20(%%r13)	,%%ymm13,%%ymm13\n\t"\
		"vmovaps		%%ymm6	,    (%%rcx)	\n\t		vmovaps		%%ymm14		,    (%%r12)\n\t"\
		"vmovaps		%%ymm7	,0x20(%%rcx)	\n\t		vmovaps		%%ymm15		,0x20(%%r12)\n\t"\
		"vmovaps		%%ymm4	,    (%%rdx)	\n\t		vaddpd	%%ymm10		,%%ymm14,%%ymm14\n\t"/* ymm14 <- ~t1 */\
		"vmovaps		%%ymm5	,0x20(%%rdx)	\n\t		vaddpd	%%ymm11		,%%ymm15,%%ymm15\n\t"/* ymm15 <- ~t2 */\
		"vaddpd	%%ymm2		,%%ymm6,%%ymm6		\n\t		vsubpd		(%%r12)	,%%ymm10,%%ymm10\n\t"/* ymm10 <- ~t5 */\
		"vaddpd	%%ymm3		,%%ymm7,%%ymm7		\n\t		vsubpd	0x20(%%r12)	,%%ymm11,%%ymm11\n\t"/* ymm11 <- ~t6 */\
		"vsubpd		(%%rcx)	,%%ymm2,%%ymm2		\n\t		vaddpd	%%ymm6,%%ymm14,%%ymm14		\n\t"/* t1+t9 */\
		"vsubpd	0x20(%%rcx)	,%%ymm3,%%ymm3		\n\t		vaddpd	%%ymm7,%%ymm15,%%ymm15		\n\t"/* t2+t10*/\
		"vmovaps		%%ymm2	,    (%%rcx)	\n\t		vaddpd	%%ymm6,%%ymm6 ,%%ymm6		\n\t"/* 2*t9  */\
		"vmovaps		%%ymm3	,0x20(%%rcx)	\n\t		vaddpd	%%ymm7,%%ymm7 ,%%ymm7		\n\t"/* 2*t10 */\
		"vmovaps		%%ymm4	,%%ymm2			\n\t		vmovaps		%%ymm14		,    (%%r10)\n\t"/* a[j1   ], DONE. */\
		"vmovaps		%%ymm5	,%%ymm3			\n\t		vmovaps		%%ymm15		,0x20(%%r10)\n\t"/* a[j2   ], DONE. */\
		"vaddpd	%%ymm0		,%%ymm5,%%ymm5		\n\t		vsubpd	%%ymm6,%%ymm14,%%ymm14		\n\t"/* t1-t9  = [t1+t9 ] - 2*t9  */\
		"vsubpd	%%ymm3		,%%ymm0,%%ymm0		\n\t		vsubpd	%%ymm7,%%ymm15,%%ymm15		\n\t"/* t2-t10 = [t2+t10] - 2*t10 */\
		"vaddpd	%%ymm1		,%%ymm4,%%ymm4		\n\t		vmovaps		%%ymm12		,%%ymm6		\n\t"/* ymm6<- copy of t7 */\
		"vsubpd	%%ymm2		,%%ymm1,%%ymm1		\n\t		vmovaps		%%ymm13		,%%ymm7		\n\t"/* ymm7<- copy of t8 */\
		"vmovaps		%%ymm5	,%%ymm2			\n\t		vaddpd	%%ymm8 		,%%ymm13,%%ymm13\n\t"/* ymm13<- ~t3 */\
		"vmovaps		%%ymm1	,%%ymm3			\n\t		vsubpd	%%ymm7		,%%ymm8 ,%%ymm8	\n\t"/* ymm8 <- ~t7 */\
		"vaddpd	%%ymm1		,%%ymm5,%%ymm5		\n\t		vaddpd	%%ymm9 		,%%ymm12,%%ymm12\n\t"/* ymm12<- ~t8 */\
		"movq		%[isrt2]	,%%rsi			\n\t		vsubpd	%%ymm6		,%%ymm9 ,%%ymm9 \n\t"/* ymm9 <- ~t4 */\
		"vmovaps			(%%rsi)	,%%ymm1		\n\t"\
		"vsubpd	%%ymm3		,%%ymm2,%%ymm2		\n\t"		/* Combine Outputs 0,4 of two half-transforms, as these are ready: */\
		"vmulpd	%%ymm1		,%%ymm5,%%ymm5		\n\t		movq		%[c4]		,%%rsi		\n\t"\
		"vmulpd	%%ymm1		,%%ymm2,%%ymm2		\n\t		vmovaps		%%ymm14,%%ymm6			\n\t"\
		"vmovaps		%%ymm0	,%%ymm3			\n\t		vmovaps		%%ymm15,%%ymm7			\n\t"\
		"vmovaps		%%ymm5	,    (%%rbx)	\n\t		vmulpd	0x20(%%rsi)	,%%ymm6	,%%ymm6 \n\t"\
		"vmovaps		%%ymm2	,0x20(%%rbx)	\n\t		vmulpd	0x20(%%rsi)	,%%ymm7	,%%ymm7 \n\t"\
		"vmovaps		%%ymm4	,%%ymm5			\n\t		vmulpd		(%%rsi)	,%%ymm14,%%ymm14\n\t"\
		"vaddpd	%%ymm4		,%%ymm0,%%ymm0		\n\t		vmulpd		(%%rsi)	,%%ymm15,%%ymm15\n\t"\
		"vsubpd	%%ymm5		,%%ymm3,%%ymm3		\n\t		vaddpd	%%ymm7 ,%%ymm14,%%ymm14		\n\t"\
		"vmulpd	%%ymm1		,%%ymm0,%%ymm0		\n\t		vsubpd	%%ymm6 ,%%ymm15,%%ymm15		\n\t"\
		"vmulpd	%%ymm1		,%%ymm3,%%ymm3		\n\t		vmovaps		%%ymm14,    (%%rax)		\n\t"\
		"vmovaps		%%ymm0	,    (%%rdx)	\n\t		vmovaps		%%ymm15,0x20(%%rax)		\n\t"\
		"vmovaps		%%ymm3	,0x20(%%rdx)	\n\t"\
	/* Now combine the two half-transforms & store outputs back into original array slots. */\
	/* add0-7 in r10,11,12,13,ax,bx,cx,dx; SIMD Registers 0-7,14-15 FREE: */\
		/* Outputs 1,5: Use ymm 4,5,9,13,14,15			Outputs 2,6: Use ymm 0,1,2,3,6,7,10,11 : */\
		"movq		%[c2],%%rsi		\n\t"/* c6 = c2+2, c1 = c2+4, c5 = c2+6 */\
		"vmovaps		    (%%rbx)	,%%ymm4 	\n\t		vmovaps		    (%%rcx)	,%%ymm2		\n\t"\
		"vmovaps		0x20(%%rbx)	,%%ymm5 	\n\t		vmovaps		0x20(%%rcx)	,%%ymm3		\n\t"\
		"vmovaps		%%ymm13	,%%ymm14		\n\t		vmovaps		%%ymm10		,%%ymm6		\n\t"\
		"vmovaps		%%ymm9 	,%%ymm15		\n\t		vmovaps		%%ymm11		,%%ymm7		\n\t"\
		"vaddpd	%%ymm4 	,%%ymm13,%%ymm13		\n\t		vaddpd	%%ymm3	,%%ymm10,%%ymm10	\n\t"\
		"vsubpd	%%ymm5 	,%%ymm9 ,%%ymm9 		\n\t		vsubpd	%%ymm2	,%%ymm11,%%ymm11	\n\t"\
		"vsubpd	%%ymm4 	,%%ymm14,%%ymm14		\n\t"	/*	vmovaps		    (%%rsi)	,%%ymm0	Need ymm0 to replace ymm9  below */\
		"vaddpd	%%ymm5 	,%%ymm15,%%ymm15		\n\t"	/*	vmovaps		0x20(%%rsi)	,%%ymm1	Need ymm1 to replace ymm13 below */\
		"vmovaps		%%ymm14	,%%ymm4 		\n\t		vsubpd	%%ymm3	,%%ymm6	,%%ymm6		\n\t"\
		"vmovaps		%%ymm15	,%%ymm5 		\n\t		vaddpd	%%ymm2	,%%ymm7	,%%ymm7		\n\t"/* ymm2,3 free */\
		"vmovaps		%%ymm13	,%%ymm14		\n\t		vmovaps		%%ymm10		,%%ymm0 	\n\t"\
		"vmovaps		%%ymm9 	,%%ymm15		\n\t		vmovaps		%%ymm11		,%%ymm1		\n\t"\
		"vmulpd	0x80(%%rsi)	,%%ymm13,%%ymm13	\n\t		vmulpd		(%%rsi)	,%%ymm10,%%ymm10\n\t"\
		"vmulpd	0x80(%%rsi)	,%%ymm9 ,%%ymm9 	\n\t		vmulpd		(%%rsi)	,%%ymm11,%%ymm11\n\t"\
		"vmulpd	0xa0(%%rsi)	,%%ymm14,%%ymm14	\n\t		vmulpd	0x20(%%rsi)	,%%ymm0 ,%%ymm0 \n\t"\
		"vmulpd	0xa0(%%rsi)	,%%ymm15,%%ymm15	\n\t		vmulpd	0x20(%%rsi)	,%%ymm1	,%%ymm1	\n\t"\
		"vsubpd	%%ymm14		,%%ymm9 ,%%ymm9 	\n\t		vsubpd	%%ymm0 		,%%ymm11,%%ymm11\n\t"\
		"vaddpd	%%ymm15		,%%ymm13,%%ymm13	\n\t		vaddpd	%%ymm1		,%%ymm10,%%ymm10\n\t"\
		"vmovaps		%%ymm9 	,0x20(%%r11)	\n\t		vmovaps		0x40(%%rsi)	,%%ymm2		\n\t"\
		"vmovaps		%%ymm13	,    (%%r11)	\n\t		vmovaps		0x60(%%rsi)	,%%ymm3		\n\t"\
		"vmovaps		%%ymm4 	,%%ymm13		\n\t		vmovaps		%%ymm11		,0x20(%%r12)\n\t"\
		"vmovaps		%%ymm5 	,%%ymm9 		\n\t		vmovaps		%%ymm10		,    (%%r12)\n\t"\
		"vmovaps		%%ymm13	,%%ymm14		\n\t		vmovaps		%%ymm6		,%%ymm0 	\n\t"\
		"vmovaps		%%ymm9 	,%%ymm15		\n\t		vmovaps		%%ymm7		,%%ymm1		\n\t"\
		"vmulpd	0xc0(%%rsi)	,%%ymm13,%%ymm13	\n\t		vmulpd	%%ymm2	,%%ymm6	,%%ymm6 	\n\t"\
		"vmulpd	0xc0(%%rsi)	,%%ymm9 ,%%ymm9 	\n\t		vmulpd	%%ymm2	,%%ymm7	,%%ymm7		\n\t"\
		"vmulpd	0xe0(%%rsi)	,%%ymm14,%%ymm14	\n\t		vmulpd	%%ymm3	,%%ymm0 ,%%ymm0		\n\t"\
		"vmulpd	0xe0(%%rsi)	,%%ymm15,%%ymm15	\n\t		vmulpd	%%ymm3	,%%ymm1	,%%ymm1		\n\t"\
		"vsubpd	%%ymm14		,%%ymm9 ,%%ymm9 	\n\t		vsubpd	%%ymm0 	,%%ymm7	,%%ymm7		\n\t"\
		"vaddpd	%%ymm15		,%%ymm13,%%ymm13	\n\t		vaddpd	%%ymm1	,%%ymm6	,%%ymm6		\n\t"\
		"vmovaps		%%ymm9 	,0x20(%%rbx)	\n\t		vmovaps		%%ymm7		,0x20(%%rcx)\n\t"\
		"vmovaps		%%ymm13	,    (%%rbx)	\n\t		vmovaps		%%ymm6		,    (%%rcx)\n\t"\
		/* Outputs 3,7: All SIMD regs except for ymm8,12 (containing one of the 2 complex data being butterflied) free: */\
		"movq		%[c3]		,%%rsi			\n\t"/* c7 = c3+2 */\
		"vmovaps		    (%%rdx)	,%%ymm4		\n\t"\
		"vmovaps		0x20(%%rdx)	,%%ymm5		\n\t"\
		"vmovaps		%%ymm8 	,%%ymm14		\n\t"\
		"vmovaps		%%ymm12	,%%ymm15		\n\t"\
		"vsubpd	%%ymm5	,%%ymm8 ,%%ymm8			\n\t		vsubpd	%%ymm4	,%%ymm12,%%ymm12	\n\t"\
		"vmovaps		    (%%rsi)	,%%ymm0		\n\t"\
		"vmovaps		0x20(%%rsi)	,%%ymm1		\n\t"\
		"vaddpd	%%ymm5	,%%ymm14,%%ymm14		\n\t		vmovaps		0x40(%%rsi)	,%%ymm2		\n\t"\
		"vaddpd	%%ymm4	,%%ymm15,%%ymm15		\n\t		vmovaps		0x60(%%rsi)	,%%ymm3		\n\t"\
		"vmovaps		%%ymm8 	,%%ymm9 		\n\t		vmovaps		%%ymm14		,%%ymm6 	\n\t"\
		"vmovaps		%%ymm12	,%%ymm13		\n\t		vmovaps		%%ymm15		,%%ymm7		\n\t"\
		"vmulpd	%%ymm0	,%%ymm8 ,%%ymm8 		\n\t		vmulpd	%%ymm2	,%%ymm14,%%ymm14	\n\t"\
		"vmulpd	%%ymm0	,%%ymm12,%%ymm12		\n\t		vmulpd	%%ymm2	,%%ymm15,%%ymm15	\n\t"\
		"vmulpd	%%ymm1	,%%ymm9 ,%%ymm9 		\n\t		vmulpd	%%ymm3	,%%ymm6 ,%%ymm6 	\n\t"\
		"vmulpd	%%ymm1	,%%ymm13,%%ymm13		\n\t		vmulpd	%%ymm3	,%%ymm7	,%%ymm7		\n\t"\
		"vsubpd	%%ymm9 		,%%ymm12,%%ymm12	\n\t		vsubpd	%%ymm6 	,%%ymm15,%%ymm15	\n\t"\
		"vaddpd	%%ymm13		,%%ymm8 ,%%ymm8 	\n\t		vaddpd	%%ymm7	,%%ymm14,%%ymm14	\n\t"\
		"vmovaps		%%ymm12	,0x20(%%r13)	\n\t		vmovaps		%%ymm15		,0x20(%%rdx)\n\t"\
		"vmovaps		%%ymm8 	,    (%%r13)	\n\t		vmovaps		%%ymm14		,    (%%rdx)\n\t"\
		:					/* outputs: none */\
		: [add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[add1] "m" (Xadd1)\
		 ,[add2] "m" (Xadd2)\
		 ,[add3] "m" (Xadd3)\
		 ,[add4] "m" (Xadd4)\
		 ,[add5] "m" (Xadd5)\
		 ,[add6] "m" (Xadd6)\
		 ,[add7] "m" (Xadd7)\
		 ,[isrt2] "m" (Xisrt2)\
		 ,[c1] "m" (Xc1)\
		 ,[c2] "m" (Xc2)\
		 ,[c3] "m" (Xc3)\
		 ,[c4] "m" (Xc4)\
		 ,[c5] "m" (Xc5)\
		 ,[c6] "m" (Xc6)\
		 ,[c7] "m" (Xc7)\
		: "cc","memory","rax","rbx","rcx","rdx","rsi","r10","r11","r12","r13","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
	);\
	}

  #endif	// USE_64BIT_ASM_STYLE ?

#elif defined(USE_SSE2)

  #define USE_64BIT_ASM_STYLE	1	// My x86 timings indicate fancier versions below using ymm0-15 for the radix-16 DFT is faster.

  #if !USE_64BIT_ASM_STYLE	// False: Use simple 64-bit-ified version of the analogous 32-bit ASM macros, i.e. using just ymm0-7.

  #if 1	// Improved version, reduce load/stores:

	#define SSE2_RADIX8_DIF_TWIDDLE(Xadd0,Xadd1,Xadd2,Xadd3,Xadd4,Xadd5,Xadd6,Xadd7,Xisrt2,Xc1,Xc2,Xc3,Xc4,Xc5,Xc6,Xc7)\
	{\
	__asm__ volatile (\
		/* Block 0,4: */\
		"movq		%[add0]	,%%rax			\n\t"\
		"movq		%[add4]	,%%rbx			\n\t"\
		"movq		%[c4]	,%%rcx			\n\t"\
		"movaps		    (%%rbx)	,%%xmm2		\n\t"\
		"movaps		0x10(%%rbx)	,%%xmm3		\n\t"\
		"movaps		    (%%rcx)	,%%xmm0		\n\t"\
		"movaps		0x10(%%rcx)	,%%xmm1		\n\t"\
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
		"movq		%[add2]	,%%rax			\n\t"\
		"movq		%[add6]	,%%rbx			\n\t"\
		"movq		%[c2]	,%%rcx			\n\t"\
		"movq		%[c6]	,%%rdx			\n\t"\
		"movaps		    (%%rax)	,%%xmm0		\n\t"\
		"movaps		0x10(%%rax)	,%%xmm2		\n\t"\
		"movaps		    (%%rcx)	,%%xmm6		\n\t"\
		"movaps		0x10(%%rcx)	,%%xmm7		\n\t"\
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
		"movaps		    (%%rdx)	,%%xmm6		\n\t"\
		"movaps		0x10(%%rdx)	,%%xmm7		\n\t"\
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
		"movq		%[add1]	,%%rax			\n\t"\
		"movq		%[add5]	,%%rbx			\n\t"\
		"movq		%[c1]	,%%rcx			\n\t"\
		"movq		%[c5]	,%%rdx			\n\t"\
		"movaps		    (%%rax)	,%%xmm0		\n\t"\
		"movaps		0x10(%%rax)	,%%xmm2		\n\t"\
		"movaps		    (%%rcx)	,%%xmm6		\n\t"\
		"movaps		0x10(%%rcx)	,%%xmm7		\n\t"\
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
		"movaps		    (%%rdx)	,%%xmm6		\n\t"\
		"movaps		0x10(%%rdx)	,%%xmm7		\n\t"\
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
		"movq		%[add3]	,%%rax			\n\t"\
		"movq		%[add7]	,%%rbx			\n\t"\
		"movq		%[c3]	,%%rcx			\n\t"\
		"movq		%[c7]	,%%rdx			\n\t"\
		"movaps		    (%%rax)	,%%xmm0		\n\t"\
		"movaps		0x10(%%rax)	,%%xmm2		\n\t"\
		"movaps		    (%%rcx)	,%%xmm6		\n\t"\
		"movaps		0x10(%%rcx)	,%%xmm7		\n\t"\
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
		"movaps		    (%%rdx)	,%%xmm6		\n\t"\
		"movaps		0x10(%%rdx)	,%%xmm7		\n\t"\
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
		"movq		%[add0]	,%%rax			\n\t"\
		"movq		%[add2]	,%%rbx			\n\t"\
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
		"movq		%[add1]	,%%rcx			\n\t"\
		"movq		%[add3]	,%%rdx			\n\t"\
		"movaps		    (%%rcx)	,%%xmm2		\n\t"\
		"movaps		0x10(%%rcx)	,%%xmm3		\n\t"\
		"movaps		    (%%rdx)	,%%xmm8		\n\t"/* Use xtra 64-bit reg xmm8 here */\
		"movaps		0x10(%%rdx)	,%%xmm9		\n\t"/* Use xtra 64-bit reg xmm9 here */\
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
		"movaps		%%xmm2		,    (%%rax)	\n\t"\
		"movaps		%%xmm3		,0x10(%%rax)	\n\t"\
		"movaps		%%xmm4		,    (%%rbx)	\n\t"\
		"movaps		%%xmm6		,0x10(%%rbx)	\n\t"\
		"movaps		%%xmm0		,    (%%rcx)	\n\t"\
		"movaps		%%xmm1		,0x10(%%rcx)	\n\t"\
		"movaps		%%xmm7		,    (%%rdx)	\n\t"\
		"movaps		%%xmm5		,0x10(%%rdx)	\n\t"\
		/* Combo 2: */\
		"movq		%[add4]	,%%rax			\n\t"\
		"movq		%[add6]	,%%rbx			\n\t"\
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
		"movq		%[add5]	,%%rcx			\n\t"\
		"movq		%[add7]	,%%rdx			\n\t"\
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
		"movq		%[isrt2]	,%%rsi		\n\t"\
		"movaps		    (%%rsi)	,%%xmm8		\n\t"/* Use xtra 64-bit reg xmm8 here */\
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
		"movaps			(%%rdx)	,%%xmm6		\n\t"\
		"subpd		%%xmm2		,%%xmm0		\n\t"\
		"subpd		%%xmm5		,%%xmm1		\n\t"\
		"subpd		%%xmm4		,%%xmm6		\n\t"\
		"subpd		%%xmm7		,%%xmm3		\n\t"\
		"addpd			(%%rax)	,%%xmm2		\n\t"\
		"addpd		0x10(%%rax)	,%%xmm5		\n\t"\
		"addpd			(%%rbx)	,%%xmm4		\n\t"\
		"addpd		0x10(%%rbx)	,%%xmm7		\n\t"\
		"movaps		%%xmm2		,    (%%rax)	\n\t"\
		"movaps		%%xmm5		,0x10(%%rax)	\n\t"\
		"movaps		%%xmm6		,    (%%rbx)	\n\t"\
		"movaps		%%xmm3		,0x10(%%rbx)	\n\t"\
		"movaps		%%xmm0		,    (%%rcx)	\n\t"\
		"movaps		%%xmm1		,0x10(%%rcx)	\n\t"\
		"movaps		%%xmm4		,    (%%rdx)	\n\t"\
		"movaps		%%xmm7		,0x10(%%rdx)	\n\t"\
		:					/* outputs: none */\
		: [add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[add1] "m" (Xadd1)\
		 ,[add2] "m" (Xadd2)\
		 ,[add3] "m" (Xadd3)\
		 ,[add4] "m" (Xadd4)\
		 ,[add5] "m" (Xadd5)\
		 ,[add6] "m" (Xadd6)\
		 ,[add7] "m" (Xadd7)\
		 ,[isrt2] "m" (Xisrt2)\
		 ,[c1] "m" (Xc1)\
		 ,[c2] "m" (Xc2)\
		 ,[c3] "m" (Xc3)\
		 ,[c4] "m" (Xc4)\
		 ,[c5] "m" (Xc5)\
		 ,[c6] "m" (Xc6)\
		 ,[c7] "m" (Xc7)\
		: "cc","memory","rax","rbx","rcx","rdx","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9"		/* Clobbered registers */\
	);\
	}

  #else	// Original version:

	#define SSE2_RADIX8_DIF_TWIDDLE(Xadd0,Xadd1,Xadd2,Xadd3,Xadd4,Xadd5,Xadd6,Xadd7,Xisrt2,Xc1,Xc2,Xc3,Xc4,Xc5,Xc6,Xc7)\
	{\
	__asm__ volatile (\
		"movq		%[add0]	,%%rax			\n\t"\
		"movq		%[add4]	,%%rbx			\n\t"\
		"movq		%[c4]	,%%rcx			\n\t"\
		"movaps		    (%%rax)	,%%xmm0		\n\t"\
		"movaps		0x10(%%rax)	,%%xmm1		\n\t"\
		"movaps		    (%%rax)	,%%xmm6		\n\t"\
		"movaps		0x10(%%rax)	,%%xmm7		\n\t"\
		"movaps		    (%%rbx)	,%%xmm2		\n\t"\
		"movaps		0x10(%%rbx)	,%%xmm3		\n\t"\
		"movaps		    (%%rbx)	,%%xmm4		\n\t"\
		"movaps		0x10(%%rbx)	,%%xmm5		\n\t"\
		"mulpd			(%%rcx)	,%%xmm2		\n\t"\
		"mulpd			(%%rcx)	,%%xmm3		\n\t"\
		"mulpd		0x10(%%rcx)	,%%xmm4		\n\t"\
		"mulpd		0x10(%%rcx)	,%%xmm5		\n\t"\
		"subpd		%%xmm5		,%%xmm2		\n\t"\
		"addpd		%%xmm4		,%%xmm3		\n\t"\
		"addpd		%%xmm2		,%%xmm0		\n\t"\
		"addpd		%%xmm3		,%%xmm1		\n\t"\
		"subpd		%%xmm2		,%%xmm6		\n\t"\
		"subpd		%%xmm3		,%%xmm7		\n\t"\
		"movaps		%%xmm0		,    (%%rax)	\n\t"\
		"movaps		%%xmm1		,0x10(%%rax)	\n\t"\
		"movaps		%%xmm6		,    (%%rbx)	\n\t"\
		"movaps		%%xmm7		,0x10(%%rbx)	\n\t"\
		"movq		%[add2]	,%%rax			\n\t"\
		"movq		%[add6]	,%%rbx			\n\t"\
		"movq		%[c2]	,%%rcx			\n\t"\
		"movq		%[c6]	,%%rdx			\n\t"\
		"movaps		    (%%rax)	,%%xmm0		\n\t"\
		"movaps		0x10(%%rax)	,%%xmm2		\n\t"\
		"movaps		    (%%rax)	,%%xmm1		\n\t"\
		"movaps		0x10(%%rax)	,%%xmm3		\n\t"\
		"mulpd			(%%rcx)	,%%xmm0		\n\t"\
		"mulpd		0x10(%%rcx)	,%%xmm2		\n\t"\
		"mulpd		0x10(%%rcx)	,%%xmm1		\n\t"\
		"mulpd			(%%rcx)	,%%xmm3		\n\t"\
		"subpd		%%xmm2		,%%xmm0		\n\t"\
		"addpd		%%xmm3		,%%xmm1		\n\t"\
		"movaps		    (%%rbx)	,%%xmm2		\n\t"\
		"movaps		0x10(%%rbx)	,%%xmm3		\n\t"\
		"movaps		    (%%rbx)	,%%xmm4		\n\t"\
		"movaps		0x10(%%rbx)	,%%xmm5		\n\t"\
		"mulpd			(%%rdx)	,%%xmm2		\n\t"\
		"mulpd		0x10(%%rdx)	,%%xmm3		\n\t"\
		"mulpd		0x10(%%rdx)	,%%xmm4		\n\t"\
		"mulpd			(%%rdx)	,%%xmm5		\n\t"\
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
		"movq		%[add1]	,%%rax			\n\t"\
		"movq		%[add5]	,%%rbx			\n\t"\
		"movq		%[c1]	,%%rcx			\n\t"\
		"movq		%[c5]	,%%rdx			\n\t"\
		"movaps		    (%%rax)	,%%xmm0		\n\t"\
		"movaps		0x10(%%rax)	,%%xmm2		\n\t"\
		"movaps		    (%%rax)	,%%xmm1		\n\t"\
		"movaps		0x10(%%rax)	,%%xmm3		\n\t"\
		"mulpd			(%%rcx)	,%%xmm0		\n\t"\
		"mulpd		0x10(%%rcx)	,%%xmm2		\n\t"\
		"mulpd		0x10(%%rcx)	,%%xmm1		\n\t"\
		"mulpd			(%%rcx)	,%%xmm3		\n\t"\
		"subpd		%%xmm2		,%%xmm0		\n\t"\
		"addpd		%%xmm3		,%%xmm1		\n\t"\
		"movaps		    (%%rbx)	,%%xmm2		\n\t"\
		"movaps		0x10(%%rbx)	,%%xmm3		\n\t"\
		"movaps		    (%%rbx)	,%%xmm4		\n\t"\
		"movaps		0x10(%%rbx)	,%%xmm5		\n\t"\
		"mulpd			(%%rdx)	,%%xmm2		\n\t"\
		"mulpd		0x10(%%rdx)	,%%xmm3		\n\t"\
		"mulpd		0x10(%%rdx)	,%%xmm4		\n\t"\
		"mulpd			(%%rdx)	,%%xmm5		\n\t"\
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
		"movq		%[add3]	,%%rax			\n\t"\
		"movq		%[add7]	,%%rbx			\n\t"\
		"movq		%[c3]	,%%rcx			\n\t"\
		"movq		%[c7]	,%%rdx			\n\t"\
		"movaps		    (%%rax)	,%%xmm0		\n\t"\
		"movaps		0x10(%%rax)	,%%xmm2		\n\t"\
		"movaps		    (%%rax)	,%%xmm1		\n\t"\
		"movaps		0x10(%%rax)	,%%xmm3		\n\t"\
		"mulpd			(%%rcx)	,%%xmm0		\n\t"\
		"mulpd		0x10(%%rcx)	,%%xmm2		\n\t"\
		"mulpd		0x10(%%rcx)	,%%xmm1		\n\t"\
		"mulpd			(%%rcx)	,%%xmm3		\n\t"\
		"subpd		%%xmm2		,%%xmm0		\n\t"\
		"addpd		%%xmm3		,%%xmm1		\n\t"\
		"movaps		    (%%rbx)	,%%xmm2		\n\t"\
		"movaps		0x10(%%rbx)	,%%xmm3		\n\t"\
		"movaps		    (%%rbx)	,%%xmm4		\n\t"\
		"movaps		0x10(%%rbx)	,%%xmm5		\n\t"\
		"mulpd			(%%rdx)	,%%xmm2		\n\t"\
		"mulpd		0x10(%%rdx)	,%%xmm3		\n\t"\
		"mulpd		0x10(%%rdx)	,%%xmm4		\n\t"\
		"mulpd			(%%rdx)	,%%xmm5		\n\t"\
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
		"movq		%[add0]	,%%rax			\n\t"\
		"movq		%[add2]	,%%rbx			\n\t"\
		"movaps		    (%%rax)	,%%xmm0		\n\t"\
		"movaps		0x10(%%rax)	,%%xmm1		\n\t"\
		"movaps		%%xmm0		,%%xmm4		\n\t"\
		"movaps		%%xmm1		,%%xmm5		\n\t"\
		"addpd			(%%rbx)	,%%xmm0		\n\t"\
		"subpd			(%%rbx)	,%%xmm4		\n\t"\
		"addpd		0x10(%%rbx)	,%%xmm1		\n\t"\
		"subpd		0x10(%%rbx)	,%%xmm5		\n\t"\
		"movaps		%%xmm0		,    (%%rax)	\n\t"\
		"movaps		%%xmm1		,0x10(%%rax)	\n\t"\
		"movaps		%%xmm4		,    (%%rbx)	\n\t"\
		"movaps		%%xmm5		,0x10(%%rbx)	\n\t"\
		"movq		%[add1]	,%%rcx			\n\t"\
		"movq		%[add3]	,%%rdx			\n\t"\
		"movaps		    (%%rcx)	,%%xmm2		\n\t"\
		"movaps		0x10(%%rcx)	,%%xmm3		\n\t"\
		"movaps		%%xmm2		,%%xmm6		\n\t"\
		"movaps		%%xmm3		,%%xmm7		\n\t"\
		"addpd			(%%rdx)	,%%xmm2		\n\t"\
		"subpd			(%%rdx)	,%%xmm6		\n\t"\
		"addpd		0x10(%%rdx)	,%%xmm3		\n\t"\
		"subpd		0x10(%%rdx)	,%%xmm7		\n\t"\
		"subpd		%%xmm2		,%%xmm0		\n\t"\
		"subpd		%%xmm3		,%%xmm1		\n\t"\
		"subpd		%%xmm7		,%%xmm4		\n\t"\
		"subpd		%%xmm6		,%%xmm5		\n\t"\
		"addpd			(%%rax)	,%%xmm2		\n\t"\
		"addpd		0x10(%%rax)	,%%xmm3		\n\t"\
		"addpd			(%%rbx)	,%%xmm7		\n\t"\
		"addpd		0x10(%%rbx)	,%%xmm6		\n\t"\
		"movaps		%%xmm2		,    (%%rax)	\n\t"\
		"movaps		%%xmm3		,0x10(%%rax)	\n\t"\
		"movaps		%%xmm4		,    (%%rbx)	\n\t"\
		"movaps		%%xmm6		,0x10(%%rbx)	\n\t"\
		"movaps		%%xmm0		,    (%%rcx)	\n\t"\
		"movaps		%%xmm1		,0x10(%%rcx)	\n\t"\
		"movaps		%%xmm7		,    (%%rdx)	\n\t"\
		"movaps		%%xmm5		,0x10(%%rdx)	\n\t"\
		"movq		%[add4]	,%%rax			\n\t"\
		"movq		%[add6]	,%%rbx			\n\t"\
		"movaps		    (%%rax)	,%%xmm0		\n\t"\
		"movaps		0x10(%%rax)	,%%xmm1		\n\t"\
		"movaps		%%xmm0		,%%xmm4		\n\t"\
		"movaps		%%xmm1		,%%xmm5		\n\t"\
		"subpd		0x10(%%rbx)	,%%xmm0		\n\t"\
		"addpd		0x10(%%rbx)	,%%xmm4		\n\t"\
		"addpd			(%%rbx)	,%%xmm1		\n\t"\
		"subpd			(%%rbx)	,%%xmm5		\n\t"\
		"movaps		%%xmm0		,    (%%rax)	\n\t"\
		"movaps		%%xmm1		,0x10(%%rax)	\n\t"\
		"movaps		%%xmm4		,    (%%rbx)	\n\t"\
		"movaps		%%xmm5		,0x10(%%rbx)	\n\t"\
		"movq		%[add5]	,%%rcx			\n\t"\
		"movq		%[add7]	,%%rdx			\n\t"\
		"movaps		    (%%rcx)	,%%xmm2		\n\t"\
		"movaps		0x10(%%rcx)	,%%xmm3		\n\t"\
		"movaps		%%xmm2		,%%xmm6		\n\t"\
		"movaps		%%xmm3		,%%xmm7		\n\t"\
		"subpd		0x10(%%rdx)	,%%xmm2		\n\t"\
		"addpd		0x10(%%rdx)	,%%xmm6		\n\t"\
		"addpd			(%%rdx)	,%%xmm3		\n\t"\
		"subpd			(%%rdx)	,%%xmm7		\n\t"\
		"movaps		%%xmm4		,    (%%rdx)	\n\t"\
		"movaps		%%xmm5		,0x10(%%rdx)	\n\t"\
		"movq		%[isrt2]	,%%rax		\n\t"\
		"movaps		%%xmm2		,%%xmm5		\n\t"\
		"subpd		%%xmm3		,%%xmm2		\n\t"\
		"addpd		%%xmm3		,%%xmm5		\n\t"\
		"mulpd			(%%rax)	,%%xmm2		\n\t"\
		"mulpd			(%%rax)	,%%xmm5		\n\t"\
		"movaps		0x10(%%rdx)	,%%xmm3		\n\t"\
		"movaps		%%xmm7		,%%xmm4		\n\t"\
		"addpd		%%xmm6		,%%xmm4		\n\t"\
		"subpd		%%xmm6		,%%xmm7		\n\t"\
		"mulpd			(%%rax)	,%%xmm4		\n\t"\
		"mulpd			(%%rax)	,%%xmm7		\n\t"\
		"movaps			(%%rdx)	,%%xmm6		\n\t"\
		"movq		%[add4]	,%%rax			\n\t"\
		"subpd		%%xmm2		,%%xmm0		\n\t"\
		"subpd		%%xmm5		,%%xmm1		\n\t"\
		"subpd		%%xmm4		,%%xmm6		\n\t"\
		"subpd		%%xmm7		,%%xmm3		\n\t"\
		"addpd			(%%rax)	,%%xmm2		\n\t"\
		"addpd		0x10(%%rax)	,%%xmm5		\n\t"\
		"addpd			(%%rbx)	,%%xmm4		\n\t"\
		"addpd		0x10(%%rbx)	,%%xmm7		\n\t"\
		"movaps		%%xmm2		,    (%%rax)	\n\t"\
		"movaps		%%xmm5		,0x10(%%rax)	\n\t"\
		"movaps		%%xmm6		,    (%%rbx)	\n\t"\
		"movaps		%%xmm3		,0x10(%%rbx)	\n\t"\
		"movaps		%%xmm0		,    (%%rcx)	\n\t"\
		"movaps		%%xmm1		,0x10(%%rcx)	\n\t"\
		"movaps		%%xmm4		,    (%%rdx)	\n\t"\
		"movaps		%%xmm7		,0x10(%%rdx)	\n\t"\
		:					/* outputs: none */\
		: [add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[add1] "m" (Xadd1)\
		 ,[add2] "m" (Xadd2)\
		 ,[add3] "m" (Xadd3)\
		 ,[add4] "m" (Xadd4)\
		 ,[add5] "m" (Xadd5)\
		 ,[add6] "m" (Xadd6)\
		 ,[add7] "m" (Xadd7)\
		 ,[isrt2] "m" (Xisrt2)\
		 ,[c1] "m" (Xc1)\
		 ,[c2] "m" (Xc2)\
		 ,[c3] "m" (Xc3)\
		 ,[c4] "m" (Xc4)\
		 ,[c5] "m" (Xc5)\
		 ,[c6] "m" (Xc6)\
		 ,[c7] "m" (Xc7)\
		: "cc","memory","rax","rbx","rcx","rdx","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
	);\
	}

  #endif	// orig/new DIF

	#define SSE2_RADIX8_DIT_TWIDDLE(Xadd0,Xadd1,Xadd2,Xadd3,Xadd4,Xadd5,Xadd6,Xadd7,Xisrt2,Xc1,Xc2,Xc3,Xc4,Xc5,Xc6,Xc7)\
	{\
	__asm__ volatile (\
		"movq		%[add4]		,%%rax		\n\t"\
		"movq		%[add5]		,%%rbx		\n\t"\
		"movaps			(%%rax)	,%%xmm0		\n\t"\
		"movaps		0x10(%%rax)	,%%xmm1		\n\t"\
		"movaps		%%xmm0		,%%xmm2		\n\t"\
		"movaps		%%xmm1		,%%xmm3		\n\t"\
		"addpd			(%%rbx)	,%%xmm2		\n\t"\
		"addpd		0x10(%%rbx)	,%%xmm3		\n\t"\
		"subpd			(%%rbx)	,%%xmm0		\n\t"\
		"subpd		0x10(%%rbx)	,%%xmm1		\n\t"\
		"movq		%[add6]		,%%rcx		\n\t"\
		"movq		%[add7]		,%%rdx		\n\t"\
		"movaps			(%%rcx)	,%%xmm4		\n\t"\
		"movaps		0x10(%%rcx)	,%%xmm5		\n\t"\
		"movaps		%%xmm4		,%%xmm6		\n\t"\
		"movaps		%%xmm5		,%%xmm7		\n\t"\
		"addpd			(%%rdx)	,%%xmm6		\n\t"\
		"addpd		0x10(%%rdx)	,%%xmm7		\n\t"\
		"subpd			(%%rdx)	,%%xmm4		\n\t"\
		"subpd		0x10(%%rdx)	,%%xmm5		\n\t"\
		"movaps		%%xmm6		,    (%%rcx)\n\t"\
		"movaps		%%xmm7		,0x10(%%rcx)\n\t"\
		"movaps		%%xmm4		,    (%%rdx)\n\t"\
		"movaps		%%xmm5		,0x10(%%rdx)\n\t"\
		"addpd		%%xmm2		,%%xmm6		\n\t"\
		"addpd		%%xmm3		,%%xmm7		\n\t"\
		"subpd			(%%rcx)	,%%xmm2		\n\t"\
		"subpd		0x10(%%rcx)	,%%xmm3		\n\t"\
		"movaps		%%xmm2		,    (%%rcx)\n\t"\
		"movaps		%%xmm3		,0x10(%%rcx)\n\t"\
		"movaps		%%xmm4		,%%xmm2		\n\t"\
		"movaps		%%xmm5		,%%xmm3		\n\t"\
		"addpd		%%xmm0		,%%xmm5		\n\t"\
		"subpd		%%xmm3		,%%xmm0		\n\t"\
		"addpd		%%xmm1		,%%xmm4		\n\t"\
		"subpd		%%xmm2		,%%xmm1		\n\t"\
		"movaps		%%xmm5		,%%xmm2		\n\t"\
		"movaps		%%xmm1		,%%xmm3		\n\t"\
		"addpd		%%xmm1		,%%xmm5		\n\t"\
		"movq		%[isrt2]	,%%rcx		\n\t"\
		"movaps			(%%rcx)	,%%xmm1		\n\t"\
		"subpd		%%xmm3		,%%xmm2		\n\t"\
		"mulpd		%%xmm1		,%%xmm5		\n\t"\
		"mulpd		%%xmm1		,%%xmm2		\n\t"\
		"movaps		%%xmm0		,%%xmm3		\n\t"\
		"movaps		%%xmm5		,    (%%rbx)\n\t"\
		"movaps		%%xmm4		,%%xmm5		\n\t"\
		"addpd		%%xmm4		,%%xmm0		\n\t"\
		"movaps		%%xmm2		,0x10(%%rbx)\n\t"\
		"subpd		%%xmm5		,%%xmm3		\n\t"\
		"mulpd		%%xmm1		,%%xmm0		\n\t"\
		"mulpd		%%xmm1		,%%xmm3		\n\t"\
		"movaps		%%xmm0		,    (%%rdx)\n\t"\
		"movaps		%%xmm3		,0x10(%%rdx)\n\t"\
		"movq		%[add0]		,%%rax		\n\t"\
		"movq		%[add1]		,%%rbx		\n\t"\
		"movaps			(%%rax)	,%%xmm0		\n\t"\
		"movaps		0x10(%%rax)	,%%xmm1		\n\t"\
		"movaps		%%xmm0		,%%xmm2		\n\t"\
		"movaps		%%xmm1		,%%xmm3		\n\t"\
		"addpd			(%%rbx)	,%%xmm2		\n\t"\
		"addpd		0x10(%%rbx)	,%%xmm3		\n\t"\
		"subpd			(%%rbx)	,%%xmm0		\n\t"\
		"subpd		0x10(%%rbx)	,%%xmm1		\n\t"\
		"movaps		%%xmm6		,    (%%rax)\n\t"\
		"movaps		%%xmm7		,0x10(%%rax)\n\t"\
		"addpd		%%xmm6		,%%xmm6		\n\t"\
		"addpd		%%xmm7		,%%xmm7		\n\t"\
		"movq		%[add4]		,%%rcx		\n\t"\
		"movaps		%%xmm6		,    (%%rcx)\n\t"\
		"movaps		%%xmm7		,0x10(%%rcx)\n\t"\
		"movq		%[add2]		,%%rcx		\n\t"\
		"movq		%[add3]		,%%rdx		\n\t"\
		"movaps			(%%rcx)	,%%xmm4		\n\t"\
		"movaps		0x10(%%rcx)	,%%xmm5		\n\t"\
		"movaps		%%xmm4		,%%xmm6		\n\t"\
		"movaps		%%xmm5		,%%xmm7		\n\t"\
		"addpd			(%%rdx)	,%%xmm6		\n\t"\
		"addpd		0x10(%%rdx)	,%%xmm7		\n\t"\
		"subpd			(%%rdx)	,%%xmm4		\n\t"\
		"subpd		0x10(%%rdx)	,%%xmm5		\n\t"\
		"movaps		%%xmm6		,    (%%rcx)\n\t"\
		"movaps		%%xmm7		,0x10(%%rcx)\n\t"\
		"addpd		%%xmm2		,%%xmm6		\n\t"\
		"addpd		%%xmm3		,%%xmm7		\n\t"\
		"subpd			(%%rcx)	,%%xmm2		\n\t"\
		"subpd		0x10(%%rcx)	,%%xmm3		\n\t"\
		"addpd			(%%rax)	,%%xmm6		\n\t"\
		"addpd		0x10(%%rax)	,%%xmm7		\n\t"\
		"movaps		%%xmm6		,    (%%rax)\n\t"\
		"movaps		%%xmm7		,0x10(%%rax)\n\t"\
		"movq		%[add4]		,%%rbx		\n\t"\
		"subpd			(%%rbx)	,%%xmm6		\n\t"\
		"subpd		0x10(%%rbx)	,%%xmm7		\n\t"\
		"movaps		%%xmm6		,    (%%rbx)\n\t"\
		"movaps		%%xmm7		,0x10(%%rbx)\n\t"\
		"movaps		%%xmm6		,    (%%rdx)\n\t"\
		"movaps		%%xmm7		,0x10(%%rdx)\n\t"\
		"movaps		%%xmm4		,%%xmm6		\n\t"\
		"movaps		%%xmm5		,%%xmm7		\n\t"\
		"addpd		%%xmm0		,%%xmm5		\n\t"\
		"subpd		%%xmm7		,%%xmm0		\n\t"\
		"addpd		%%xmm1		,%%xmm4		\n\t"\
		"subpd		%%xmm6		,%%xmm1		\n\t"\
	/* Now combine the two half-transforms & store outputs back into original array slots: */\
		"movq		%[c4]		,%%rcx		\n\t"\
		"movaps		    (%%rbx)	,%%xmm6		\n\t"\
		"movaps		0x10(%%rbx)	,%%xmm7		\n\t"\
		"mulpd		0x10(%%rcx)	,%%xmm6		\n\t"\
		"mulpd		0x10(%%rcx)	,%%xmm7		\n\t"\
		"movaps		%%xmm6		,    (%%rbx)\n\t"\
		"movaps		%%xmm7		,0x10(%%rbx)\n\t"\
		"movaps		    (%%rdx)	,%%xmm6		\n\t"\
		"movaps		0x10(%%rdx)	,%%xmm7		\n\t"\
		"mulpd			(%%rcx)	,%%xmm6		\n\t"\
		"mulpd			(%%rcx)	,%%xmm7		\n\t"\
		"addpd		0x10(%%rbx)	,%%xmm6		\n\t"\
		"subpd			(%%rbx)	,%%xmm7		\n\t"\
		"movaps		%%xmm6		,    (%%rbx)\n\t"\
		"movaps		%%xmm7		,0x10(%%rbx)\n\t"\
		"movq		%[add1]		,%%rax		\n\t"\
		"movq		%[add5]		,%%rbx		\n\t"\
		"movq		%[c1]		,%%rcx		\n\t"\
		"movq		%[c5]		,%%rdx		\n\t"\
		"movaps		%%xmm5		,%%xmm6		\n\t"\
		"movaps		%%xmm1		,%%xmm7		\n\t"\
		"addpd			(%%rbx)	,%%xmm5		\n\t"\
		"subpd		0x10(%%rbx)	,%%xmm1		\n\t"\
		"subpd			(%%rbx)	,%%xmm6		\n\t"\
		"addpd		0x10(%%rbx)	,%%xmm7		\n\t"\
		"movaps		%%xmm6		,    (%%rbx)\n\t"\
		"movaps		%%xmm7		,0x10(%%rbx)\n\t"\
		"movaps		%%xmm5		,%%xmm6		\n\t"\
		"movaps		%%xmm1		,%%xmm7		\n\t"\
		"mulpd			(%%rcx)	,%%xmm5		\n\t"\
		"mulpd			(%%rcx)	,%%xmm1		\n\t"\
		"mulpd		0x10(%%rcx)	,%%xmm6		\n\t"\
		"mulpd		0x10(%%rcx)	,%%xmm7		\n\t"\
		"subpd		%%xmm6		,%%xmm1		\n\t"\
		"addpd		%%xmm7		,%%xmm5		\n\t"\
		"movaps		%%xmm1		,0x10(%%rax)\n\t"\
		"movaps		%%xmm5		,    (%%rax)\n\t"\
		"movaps		    (%%rbx)	,%%xmm5		\n\t"\
		"movaps		0x10(%%rbx)	,%%xmm1		\n\t"\
		"movaps		%%xmm5		,%%xmm6		\n\t"\
		"movaps		%%xmm1		,%%xmm7		\n\t"\
		"mulpd			(%%rdx)	,%%xmm5		\n\t"\
		"mulpd			(%%rdx)	,%%xmm1		\n\t"\
		"mulpd		0x10(%%rdx)	,%%xmm6		\n\t"\
		"mulpd		0x10(%%rdx)	,%%xmm7		\n\t"\
		"subpd		%%xmm6		,%%xmm1		\n\t"\
		"addpd		%%xmm7		,%%xmm5		\n\t"\
		"movaps		%%xmm1		,0x10(%%rbx)\n\t"\
		"movaps		%%xmm5		,    (%%rbx)\n\t"\
		"movq		%[add2]		,%%rax		\n\t"\
		"movq		%[add6]		,%%rbx		\n\t"\
		"movq		%[c2]		,%%rcx		\n\t"\
		"movq		%[c6]		,%%rdx		\n\t"\
		"movaps		%%xmm2		,%%xmm6		\n\t"\
		"movaps		%%xmm3		,%%xmm7		\n\t"\
		"addpd		0x10(%%rbx)	,%%xmm2		\n\t"\
		"subpd			(%%rbx)	,%%xmm3		\n\t"\
		"subpd		0x10(%%rbx)	,%%xmm6		\n\t"\
		"addpd			(%%rbx)	,%%xmm7		\n\t"\
		"movaps		%%xmm2		,%%xmm1		\n\t"\
		"movaps		%%xmm3		,%%xmm5		\n\t"\
		"mulpd			(%%rcx)	,%%xmm2		\n\t"\
		"mulpd			(%%rcx)	,%%xmm3		\n\t"\
		"mulpd		0x10(%%rcx)	,%%xmm1		\n\t"\
		"mulpd		0x10(%%rcx)	,%%xmm5		\n\t"\
		"subpd		%%xmm1		,%%xmm3		\n\t"\
		"addpd		%%xmm5		,%%xmm2		\n\t"\
		"movaps		%%xmm3		,0x10(%%rax)\n\t"\
		"movaps		%%xmm2		,    (%%rax)\n\t"\
		"movaps		%%xmm6		,%%xmm1		\n\t"\
		"movaps		%%xmm7		,%%xmm5		\n\t"\
		"mulpd			(%%rdx)	,%%xmm6		\n\t"\
		"mulpd			(%%rdx)	,%%xmm7		\n\t"\
		"mulpd		0x10(%%rdx)	,%%xmm1		\n\t"\
		"mulpd		0x10(%%rdx)	,%%xmm5		\n\t"\
		"subpd		%%xmm1		,%%xmm7		\n\t"\
		"addpd		%%xmm5		,%%xmm6		\n\t"\
		"movaps		%%xmm7		,0x10(%%rbx)\n\t"\
		"movaps		%%xmm6		,    (%%rbx)\n\t"\
		"movq		%[add3]		,%%rax		\n\t"\
		"movq		%[add7]		,%%rbx		\n\t"\
		"movq		%[c3]		,%%rcx		\n\t"\
		"movq		%[c7]		,%%rdx		\n\t"\
		"movaps		%%xmm0		,%%xmm6		\n\t"\
		"movaps		%%xmm4		,%%xmm7		\n\t"\
		"subpd		0x10(%%rbx)	,%%xmm0		\n\t"\
		"subpd			(%%rbx)	,%%xmm4		\n\t"\
		"addpd		0x10(%%rbx)	,%%xmm6		\n\t"\
		"addpd			(%%rbx)	,%%xmm7		\n\t"\
		"movaps		%%xmm0		,%%xmm1		\n\t"\
		"movaps		%%xmm4		,%%xmm5		\n\t"\
		"mulpd			(%%rcx)	,%%xmm0		\n\t"\
		"mulpd			(%%rcx)	,%%xmm4		\n\t"\
		"mulpd		0x10(%%rcx)	,%%xmm1		\n\t"\
		"mulpd		0x10(%%rcx)	,%%xmm5		\n\t"\
		"subpd		%%xmm1		,%%xmm4		\n\t"\
		"addpd		%%xmm5		,%%xmm0		\n\t"\
		"movaps		%%xmm4		,0x10(%%rax)\n\t"\
		"movaps		%%xmm0		,    (%%rax)\n\t"\
		"movaps		%%xmm6		,%%xmm1		\n\t"\
		"movaps		%%xmm7		,%%xmm5		\n\t"\
		"mulpd			(%%rdx)	,%%xmm6		\n\t"\
		"mulpd			(%%rdx)	,%%xmm7		\n\t"\
		"mulpd		0x10(%%rdx)	,%%xmm1		\n\t"\
		"mulpd		0x10(%%rdx)	,%%xmm5		\n\t"\
		"subpd		%%xmm1		,%%xmm7		\n\t"\
		"addpd		%%xmm5		,%%xmm6		\n\t"\
		"movaps		%%xmm7		,0x10(%%rbx)\n\t"\
		"movaps		%%xmm6		,    (%%rbx)\n\t"\
		:					/* outputs: none */\
		: [add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[add1] "m" (Xadd1)\
		 ,[add2] "m" (Xadd2)\
		 ,[add3] "m" (Xadd3)\
		 ,[add4] "m" (Xadd4)\
		 ,[add5] "m" (Xadd5)\
		 ,[add6] "m" (Xadd6)\
		 ,[add7] "m" (Xadd7)\
		 ,[isrt2] "m" (Xisrt2)\
		 ,[c1] "m" (Xc1)\
		 ,[c2] "m" (Xc2)\
		 ,[c3] "m" (Xc3)\
		 ,[c4] "m" (Xc4)\
		 ,[c5] "m" (Xc5)\
		 ,[c6] "m" (Xc6)\
		 ,[c7] "m" (Xc7)\
		: "cc","memory","rax","rbx","rcx","rdx","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
	);\
	}

  #else // USE_64BIT_ASM_STYLE = True: Deeper 64-bit-ified version of the above 32-bit ASM macros, using all of ymm0-15

	// To-do: propagate load-saving tweaks from 8-reg version to this one:
	#define SSE2_RADIX8_DIF_TWIDDLE(Xadd0,Xadd1,Xadd2,Xadd3,Xadd4,Xadd5,Xadd6,Xadd7,Xisrt2,Xc1,Xc2,Xc3,Xc4,Xc5,Xc6,Xc7)\
	{\
	__asm__ volatile (\
		"												movq		%[add1]	,%%r10			\n\t"\
		"												movq		%[add5]	,%%r11			\n\t"\
		"												movq		%[c1]	,%%r12			\n\t"\
		"												movq		%[c5]	,%%r13			\n\t"\
		"												movaps		    (%%r10)	,%%xmm8 	\n\t"\
		"movq		%[add0]	,%%rax			\n\t		movaps		0x10(%%r10)	,%%xmm10	\n\t"\
		"movq		%[add4]	,%%rbx			\n\t		movaps		    (%%r10)	,%%xmm9 	\n\t"\
		"movq		%[c4]	,%%rcx			\n\t		movaps		0x10(%%r10)	,%%xmm11	\n\t"\
		"movaps		    (%%rax)	,%%xmm0		\n\t		mulpd			(%%r12)	,%%xmm8 	\n\t"\
		"movaps		0x10(%%rax)	,%%xmm1		\n\t		mulpd		0x10(%%r12)	,%%xmm10	\n\t"\
		"movaps		    (%%rax)	,%%xmm6		\n\t		mulpd		0x10(%%r12)	,%%xmm9 	\n\t"\
		"movaps		0x10(%%rax)	,%%xmm7		\n\t		mulpd			(%%r12)	,%%xmm11	\n\t"\
		"movaps		    (%%rbx)	,%%xmm2		\n\t		subpd		%%xmm10		,%%xmm8 	\n\t"\
		"movaps		0x10(%%rbx)	,%%xmm3		\n\t		addpd		%%xmm11		,%%xmm9 	\n\t"\
		"movaps		    (%%rbx)	,%%xmm4		\n\t		movaps		    (%%r11)	,%%xmm10	\n\t"\
		"movaps		0x10(%%rbx)	,%%xmm5		\n\t		movaps		0x10(%%r11)	,%%xmm11	\n\t"\
		"mulpd			(%%rcx)	,%%xmm2		\n\t		movaps		    (%%r11)	,%%xmm12	\n\t"\
		"mulpd			(%%rcx)	,%%xmm3		\n\t		movaps		0x10(%%r11)	,%%xmm13	\n\t"\
		"mulpd		0x10(%%rcx)	,%%xmm4		\n\t		mulpd			(%%r13)	,%%xmm10	\n\t"\
		"mulpd		0x10(%%rcx)	,%%xmm5		\n\t		mulpd		0x10(%%r13)	,%%xmm11	\n\t"\
		"subpd		%%xmm5		,%%xmm2		\n\t		mulpd		0x10(%%r13)	,%%xmm12	\n\t"\
		"addpd		%%xmm4		,%%xmm3		\n\t		mulpd			(%%r13)	,%%xmm13	\n\t"\
		"addpd		%%xmm2		,%%xmm0		\n\t		subpd		%%xmm11		,%%xmm10	\n\t"\
		"addpd		%%xmm3		,%%xmm1		\n\t		addpd		%%xmm13		,%%xmm12	\n\t"\
		"subpd		%%xmm2		,%%xmm6		\n\t		movaps		%%xmm10		,%%xmm11	\n\t"\
		"subpd		%%xmm3		,%%xmm7		\n\t		movaps		%%xmm12		,%%xmm13	\n\t"\
		"movaps		%%xmm0	,    (%%rax)	\n\t		addpd		%%xmm8 		,%%xmm10	\n\t"\
		"movaps		%%xmm1	,0x10(%%rax)	\n\t		subpd		%%xmm11		,%%xmm8 	\n\t"\
		"movaps		%%xmm6	,    (%%rbx)	\n\t		addpd		%%xmm9 		,%%xmm12	\n\t"\
		"movaps		%%xmm7	,0x10(%%rbx)	\n\t		subpd		%%xmm13		,%%xmm9 	\n\t"\
		"movq		%[add2]	,%%rax			\n\t		movaps		%%xmm10	,    (%%r10)	\n\t"\
		"movq		%[add6]	,%%rbx			\n\t		movaps		%%xmm12	,0x10(%%r10)	\n\t"\
		"movq		%[c2]	,%%rcx			\n\t		movaps		%%xmm8 	,    (%%r11)	\n\t"\
		"movq		%[c6]	,%%rdx			\n\t		movaps		%%xmm9 	,0x10(%%r11)	\n\t"\
		"movaps		    (%%rax)	,%%xmm0		\n\t		movq		%[add3]	,%%r10			\n\t"\
		"movaps		0x10(%%rax)	,%%xmm2		\n\t		movq		%[add7]	,%%r11			\n\t"\
		"movaps		    (%%rax)	,%%xmm1		\n\t		movq		%[c3]	,%%r12			\n\t"\
		"movaps		0x10(%%rax)	,%%xmm3		\n\t		movq		%[c7]	,%%r13			\n\t"\
		"mulpd			(%%rcx)	,%%xmm0		\n\t		movaps		    (%%r10)	,%%xmm8 	\n\t"\
		"mulpd		0x10(%%rcx)	,%%xmm2		\n\t		movaps		0x10(%%r10)	,%%xmm10	\n\t"\
		"mulpd		0x10(%%rcx)	,%%xmm1		\n\t		movaps		    (%%r10)	,%%xmm9 	\n\t"\
		"mulpd			(%%rcx)	,%%xmm3		\n\t		movaps		0x10(%%r10)	,%%xmm11	\n\t"\
		"subpd		%%xmm2		,%%xmm0		\n\t		mulpd			(%%r12)	,%%xmm8 	\n\t"\
		"addpd		%%xmm3		,%%xmm1		\n\t		mulpd		0x10(%%r12)	,%%xmm10	\n\t"\
		"movaps		    (%%rbx)	,%%xmm2		\n\t		mulpd		0x10(%%r12)	,%%xmm9 	\n\t"\
		"movaps		0x10(%%rbx)	,%%xmm3		\n\t		mulpd			(%%r12)	,%%xmm11	\n\t"\
		"movaps		    (%%rbx)	,%%xmm4		\n\t		subpd		%%xmm10		,%%xmm8 	\n\t"\
		"movaps		0x10(%%rbx)	,%%xmm5		\n\t		addpd		%%xmm11		,%%xmm9 	\n\t"\
		"mulpd			(%%rdx)	,%%xmm2		\n\t		movaps		    (%%r11)	,%%xmm10	\n\t"\
		"mulpd		0x10(%%rdx)	,%%xmm3		\n\t		movaps		0x10(%%r11)	,%%xmm11	\n\t"\
		"mulpd		0x10(%%rdx)	,%%xmm4		\n\t		movaps		    (%%r11)	,%%xmm12	\n\t"\
		"mulpd			(%%rdx)	,%%xmm5		\n\t		movaps		0x10(%%r11)	,%%xmm13	\n\t"\
		"subpd		%%xmm3		,%%xmm2		\n\t		mulpd			(%%r13)	,%%xmm10	\n\t"\
		"addpd		%%xmm5		,%%xmm4		\n\t		mulpd		0x10(%%r13)	,%%xmm11	\n\t"\
		"movaps		%%xmm2		,%%xmm3		\n\t		mulpd		0x10(%%r13)	,%%xmm12	\n\t"\
		"movaps		%%xmm4		,%%xmm5		\n\t		mulpd			(%%r13)	,%%xmm13	\n\t"\
		"addpd		%%xmm0		,%%xmm2		\n\t		subpd		%%xmm11		,%%xmm10	\n\t"\
		"subpd		%%xmm3		,%%xmm0		\n\t		addpd		%%xmm13		,%%xmm12	\n\t"\
		"addpd		%%xmm1		,%%xmm4		\n\t		movaps		%%xmm10		,%%xmm11	\n\t"\
		"subpd		%%xmm5		,%%xmm1		\n\t		movaps		%%xmm12		,%%xmm13	\n\t"\
		"movaps		%%xmm2	,    (%%rax)	\n\t		addpd		%%xmm8 		,%%xmm10	\n\t"\
		"movaps		%%xmm4	,0x10(%%rax)	\n\t		subpd		%%xmm11		,%%xmm8 	\n\t"\
		"movaps		%%xmm0	,    (%%rbx)	\n\t		addpd		%%xmm9 		,%%xmm12	\n\t"\
		"movaps		%%xmm1	,0x10(%%rbx)	\n\t		subpd		%%xmm13		,%%xmm9 	\n\t"\
		"												movaps		%%xmm10	,    (%%r10)	\n\t"\
		"												movaps		%%xmm12	,0x10(%%r10)	\n\t"\
		"												movaps		%%xmm8 	,    (%%r11)	\n\t"\
		"												movaps		%%xmm9 	,0x10(%%r11)	\n\t"\
/* combine to get 2 length-4 output subtransforms... */\
		"movq		%[add0]	,%%rax			\n\t		movq		%[add4]	,%%r10			\n\t"\
		"movq		%[add2]	,%%rbx			\n\t		movq		%[add6]	,%%r11			\n\t"\
		"movaps		    (%%rax)	,%%xmm0		\n\t		movaps		    (%%r10)	,%%xmm8 	\n\t"\
		"movaps		0x10(%%rax)	,%%xmm1		\n\t		movaps		0x10(%%r10)	,%%xmm9 	\n\t"\
		"movaps		%%xmm0		,%%xmm4		\n\t		movaps		%%xmm8 		,%%xmm12	\n\t"\
		"movaps		%%xmm1		,%%xmm5		\n\t		movaps		%%xmm9 		,%%xmm13	\n\t"\
		"addpd			(%%rbx)	,%%xmm0		\n\t		subpd		0x10(%%r11)	,%%xmm8 	\n\t"\
		"subpd			(%%rbx)	,%%xmm4		\n\t		addpd		0x10(%%r11)	,%%xmm12	\n\t"\
		"addpd		0x10(%%rbx)	,%%xmm1		\n\t		addpd			(%%r11)	,%%xmm9 	\n\t"\
		"subpd		0x10(%%rbx)	,%%xmm5		\n\t		subpd			(%%r11)	,%%xmm13	\n\t"\
		"movaps		%%xmm0	,    (%%rax)	\n\t		movaps		%%xmm8 	,    (%%r10)	\n\t"\
		"movaps		%%xmm1	,0x10(%%rax)	\n\t		movaps		%%xmm9 	,0x10(%%r10)	\n\t"\
		"movaps		%%xmm4	,    (%%rbx)	\n\t		movaps		%%xmm12	,    (%%r11)	\n\t"\
		"movaps		%%xmm5	,0x10(%%rbx)	\n\t		movaps		%%xmm13	,0x10(%%r11)	\n\t"\
		"movq		%[add1]	,%%rcx			\n\t		movq		%[add5]	,%%r12			\n\t"\
		"movq		%[add3]	,%%rdx			\n\t		movq		%[add7]	,%%r13			\n\t"\
		"movaps		    (%%rcx)	,%%xmm2		\n\t		movaps		    (%%r12)	,%%xmm10	\n\t"\
		"movaps		0x10(%%rcx)	,%%xmm3		\n\t		movaps		0x10(%%r12)	,%%xmm11	\n\t"\
		"movaps		%%xmm2		,%%xmm6		\n\t		movaps		%%xmm10		,%%xmm14	\n\t"\
		"movaps		%%xmm3		,%%xmm7		\n\t		movaps		%%xmm11		,%%xmm15	\n\t"\
		"addpd			(%%rdx)	,%%xmm2		\n\t		subpd		0x10(%%r13)	,%%xmm10	\n\t"\
		"subpd			(%%rdx)	,%%xmm6		\n\t		addpd		0x10(%%r13)	,%%xmm14	\n\t"\
		"addpd		0x10(%%rdx)	,%%xmm3		\n\t		addpd			(%%r13)	,%%xmm11	\n\t"\
		"subpd		0x10(%%rdx)	,%%xmm7		\n\t		subpd			(%%r13)	,%%xmm15	\n\t"\
		"subpd		%%xmm2		,%%xmm0		\n\t		movaps		%%xmm12	,    (%%r13)	\n\t"\
		"subpd		%%xmm3		,%%xmm1		\n\t		movaps		%%xmm13	,0x10(%%r13)	\n\t"\
		"subpd		%%xmm7		,%%xmm4		\n\t		movq		%[isrt2]	,%%r10		\n\t"\
		"subpd		%%xmm6		,%%xmm5		\n\t		movaps		%%xmm10		,%%xmm13	\n\t"\
		"addpd			(%%rax)	,%%xmm2		\n\t		subpd		%%xmm11		,%%xmm10	\n\t"\
		"addpd		0x10(%%rax)	,%%xmm3		\n\t		addpd		%%xmm11		,%%xmm13	\n\t"\
		"addpd			(%%rbx)	,%%xmm7		\n\t		mulpd			(%%r10)	,%%xmm10	\n\t"\
		"addpd		0x10(%%rbx)	,%%xmm6		\n\t		mulpd			(%%r10)	,%%xmm13	\n\t"\
		"movaps		%%xmm2	,    (%%rax)	\n\t		movaps		0x10(%%r13)	,%%xmm11	\n\t"\
		"movaps		%%xmm3	,0x10(%%rax)	\n\t		movaps		%%xmm15		,%%xmm12	\n\t"\
		"movaps		%%xmm4	,    (%%rbx)	\n\t		addpd		%%xmm14		,%%xmm12	\n\t"\
		"movaps		%%xmm6	,0x10(%%rbx)	\n\t		subpd		%%xmm14		,%%xmm15	\n\t"\
		"movaps		%%xmm0	,    (%%rcx)	\n\t		mulpd			(%%r10)	,%%xmm12	\n\t"\
		"movaps		%%xmm1	,0x10(%%rcx)	\n\t		mulpd			(%%r10)	,%%xmm15	\n\t"\
		"movaps		%%xmm7	,    (%%rdx)	\n\t		movaps			(%%r13)	,%%xmm14	\n\t"\
		"movaps		%%xmm5	,0x10(%%rdx)	\n\t		movq		%[add4]	,%%r10			\n\t"\
		"												subpd		%%xmm10		,%%xmm8 	\n\t"\
		"												subpd		%%xmm13		,%%xmm9 	\n\t"\
		"												subpd		%%xmm12		,%%xmm14	\n\t"\
		"												subpd		%%xmm15		,%%xmm11	\n\t"\
		"												addpd			(%%r10)	,%%xmm10	\n\t"\
		"												addpd		0x10(%%r10)	,%%xmm13	\n\t"\
		"												addpd			(%%r11)	,%%xmm12	\n\t"\
		"												addpd		0x10(%%r11)	,%%xmm15	\n\t"\
		"												movaps		%%xmm10	,    (%%r10)	\n\t"\
		"												movaps		%%xmm13	,0x10(%%r10)	\n\t"\
		"												movaps		%%xmm14	,    (%%r11)	\n\t"\
		"												movaps		%%xmm11	,0x10(%%r11)	\n\t"\
		"												movaps		%%xmm8 	,    (%%r12)	\n\t"\
		"												movaps		%%xmm9 	,0x10(%%r12)	\n\t"\
		"												movaps		%%xmm12	,    (%%r13)	\n\t"\
		"												movaps		%%xmm15	,0x10(%%r13)	\n\t"\
		:					/* outputs: none */\
		: [add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[add1] "m" (Xadd1)\
		 ,[add2] "m" (Xadd2)\
		 ,[add3] "m" (Xadd3)\
		 ,[add4] "m" (Xadd4)\
		 ,[add5] "m" (Xadd5)\
		 ,[add6] "m" (Xadd6)\
		 ,[add7] "m" (Xadd7)\
		 ,[isrt2] "m" (Xisrt2)\
		 ,[c1] "m" (Xc1)\
		 ,[c2] "m" (Xc2)\
		 ,[c3] "m" (Xc3)\
		 ,[c4] "m" (Xc4)\
		 ,[c5] "m" (Xc5)\
		 ,[c6] "m" (Xc6)\
		 ,[c7] "m" (Xc7)\
		: "cc","memory","rax","rbx","rcx","rdx","r10","r11","r12","r13","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
	);\
	}

	// Note: The SIMD radix-8 DIF is a very basic one, no opts re. memory movement, etc. Thus DIF has 65% more load/stores than DIT
	// DIF SIMD opcount: 140 load/store [56 implicit], 66 add/sub, 32 mul
	// DIT SIMD opcount:  85 load/store [36 implicit], 68 add/sub, 32 mul	<*** To-Do: propagate same optimizations used in DIT to DIF!

	#define SSE2_RADIX8_DIT_TWIDDLE(Xadd0,Xadd1,Xadd2,Xadd3,Xadd4,Xadd5,Xadd6,Xadd7,Xisrt2,Xc1,Xc2,Xc3,Xc4,Xc5,Xc6,Xc7)\
	{\
	__asm__ volatile (\
		/*** 2nd of 2 length-4 subtransforms gets done first: ***/\
		/* SSE2_RADIX4_DIT_0TWIDDLE_2NDOFTWO(add4-7):	SSE2_RADIX4_DIT_0TWIDDLE(add0-3): */\
		"movq		%[add4]		,%%rax		\n\t		movq		%[add0]		,%%r10		\n\t"\
		"movq		%[add5]		,%%rbx		\n\t		movq		%[add1]		,%%r11		\n\t"\
		"movaps		(%%rax)	,%%xmm0			\n\t		movaps			(%%r10)	,%%xmm8 		\n\t"\
		"movaps	0x10(%%rax)	,%%xmm1			\n\t		movaps		0x10(%%r10)	,%%xmm9 		\n\t"\
		"movaps	%%xmm0		,%%xmm2			\n\t		movaps		%%xmm8 		,%%xmm10		\n\t"\
		"movaps	%%xmm1		,%%xmm3			\n\t		movaps		%%xmm9 		,%%xmm11		\n\t"\
		"addpd			(%%rbx)	,%%xmm2		\n\t		addpd			(%%r11)	,%%xmm10		\n\t"\
		"addpd		0x10(%%rbx)	,%%xmm3		\n\t		addpd		0x10(%%r11)	,%%xmm11		\n\t"\
		"subpd			(%%rbx)	,%%xmm0		\n\t		subpd			(%%r11)	,%%xmm8 		\n\t"\
		"subpd		0x10(%%rbx)	,%%xmm1		\n\t		subpd		0x10(%%r11)	,%%xmm9 		\n\t"\
		"movq		%[add6]		,%%rcx		\n\t		movq		%[add2]		,%%r12			\n\t"\
		"movq		%[add7]		,%%rdx		\n\t		movq		%[add3]		,%%r13			\n\t"\
		"movaps			(%%rcx)	,%%xmm4		\n\t		movaps			(%%r12)	,%%xmm12		\n\t"\
		"movaps		0x10(%%rcx)	,%%xmm5		\n\t		movaps		0x10(%%r12)	,%%xmm13		\n\t"/* t9  in xmm6, needed below in RHS! */\
		"movaps		%%xmm4		,%%xmm6		\n\t		movaps		%%xmm12		,%%xmm14		\n\t"/* t10 in xmm7, needed below in RHS! */\
		"movaps		%%xmm5		,%%xmm7		\n\t		movaps		%%xmm13		,%%xmm15		\n\t"\
		"addpd			(%%rdx)	,%%xmm6		\n\t		addpd			(%%r13)	,%%xmm14		\n\t"\
		"addpd		0x10(%%rdx)	,%%xmm7		\n\t		addpd		0x10(%%r13)	,%%xmm15		\n\t"\
		"subpd			(%%rdx)	,%%xmm4		\n\t		subpd			(%%r13)	,%%xmm12		\n\t"\
		"subpd		0x10(%%rdx)	,%%xmm5		\n\t		subpd		0x10(%%r13)	,%%xmm13		\n\t"\
		"movaps		%%xmm6		,    (%%rcx)\n\t		movaps		%%xmm14		,    (%%r12)	\n\t"\
		"movaps		%%xmm7		,0x10(%%rcx)\n\t		movaps		%%xmm15		,0x10(%%r12)	\n\t"\
		"movaps		%%xmm4		,    (%%rdx)\n\t		addpd		%%xmm10		,%%xmm14		\n\t"/* xmm14 <- ~t1 */\
		"movaps		%%xmm5		,0x10(%%rdx)\n\t		addpd		%%xmm11		,%%xmm15		\n\t"/* xmm15 <- ~t2 */\
		"addpd		%%xmm2		,%%xmm6		\n\t		subpd			(%%r12)	,%%xmm10		\n\t"/* xmm10 <- ~t5 */\
		"addpd		%%xmm3		,%%xmm7		\n\t		subpd		0x10(%%r12)	,%%xmm11		\n\t"/* xmm11 <- ~t6 */\
		"subpd			(%%rcx)	,%%xmm2		\n\t		addpd		%%xmm6,%%xmm14		\n\t"/* t1+t9 */\
		"subpd		0x10(%%rcx)	,%%xmm3		\n\t		addpd		%%xmm7,%%xmm15		\n\t"/* t2+t10*/\
		"movaps		%%xmm2		,    (%%rcx)\n\t		addpd		%%xmm6,%%xmm6		\n\t"/* 2*t9  */\
		"movaps		%%xmm3		,0x10(%%rcx)\n\t		addpd		%%xmm7,%%xmm7		\n\t"/* 2*t10 */\
		"movaps		%%xmm4		,%%xmm2		\n\t		movaps		%%xmm14		,    (%%r10)	\n\t"/* a[j1   ], DONE. */\
		"movaps		%%xmm5		,%%xmm3		\n\t		movaps		%%xmm15		,0x10(%%r10)	\n\t"/* a[j2   ], DONE. */\
		"addpd		%%xmm0		,%%xmm5		\n\t		subpd		%%xmm6,%%xmm14		\n\t"/* t1-t9  = [t1+t9 ] - 2*t9  */\
		"subpd		%%xmm3		,%%xmm0		\n\t		subpd		%%xmm7,%%xmm15		\n\t"/* t2-t10 = [t2+t10] - 2*t10 */\
		"addpd		%%xmm1		,%%xmm4		\n\t		movaps		%%xmm12		,%%xmm6		\n\t"/* xmm6<- copy of t7 */\
		"subpd		%%xmm2		,%%xmm1		\n\t		movaps		%%xmm13		,%%xmm7		\n\t"/* xmm7<- copy of t8 */\
		"movaps		%%xmm5		,%%xmm2		\n\t		addpd		%%xmm8 		,%%xmm13		\n\t"/* xmm13<- ~t3 */\
		"movaps		%%xmm1		,%%xmm3		\n\t		subpd		%%xmm7		,%%xmm8 		\n\t"/* xmm8 <- ~t7 */\
		"addpd		%%xmm1		,%%xmm5		\n\t		addpd		%%xmm9 		,%%xmm12		\n\t"/* xmm12<- ~t8 */\
		"movq		%[isrt2]	,%%rsi		\n\t		subpd		%%xmm6		,%%xmm9 		\n\t"/* xmm9 <- ~t4 */\
		"movaps			(%%rsi)	,%%xmm1		\n\t"\
		"subpd		%%xmm3		,%%xmm2		\n\t"		/* Combine Outputs 0,4 of two half-transforms, as these are ready: */\
		"mulpd		%%xmm1		,%%xmm5		\n\t		movq		%[c4]		,%%rsi		\n\t"\
		"mulpd		%%xmm1		,%%xmm2		\n\t		movaps		%%xmm14,%%xmm6			\n\t"\
		"movaps		%%xmm0		,%%xmm3		\n\t		movaps		%%xmm15,%%xmm7			\n\t"\
		"movaps		%%xmm5		,    (%%rbx)\n\t		mulpd		0x10(%%rsi)	,%%xmm6		\n\t"\
		"movaps		%%xmm2		,0x10(%%rbx)\n\t		mulpd		0x10(%%rsi)	,%%xmm7		\n\t"\
		"movaps		%%xmm4		,%%xmm5		\n\t		mulpd			(%%rsi)	,%%xmm14	\n\t"\
		"addpd		%%xmm4		,%%xmm0		\n\t		mulpd			(%%rsi)	,%%xmm15	\n\t"\
		"subpd		%%xmm5		,%%xmm3		\n\t		addpd		%%xmm7 ,%%xmm14			\n\t"\
		"mulpd		%%xmm1		,%%xmm0		\n\t		subpd		%%xmm6 ,%%xmm15			\n\t"\
		"mulpd		%%xmm1		,%%xmm3		\n\t		movaps		%%xmm14,    (%%rax)		\n\t"\
		"movaps		%%xmm0		,    (%%rdx)\n\t		movaps		%%xmm15,0x10(%%rax)		\n\t"\
		"movaps		%%xmm3		,0x10(%%rdx)\n\t"\
	/* Now combine the two half-transforms & store outputs back into original array slots. */\
	/* add0-7 in r10,11,12,13,ax,bx,cx,dx; SIMD Registers 0-7,14-15 FREE: */\
		/* Outputs 1,5: Use xmm 4,5,9,13,14,15			Outputs 2,6: Use xmm 0,1,2,3,6,7,10,11 : */\
		"movq		%[c2],%%rsi		\n\t"/* c6 = c2+2, c1 = c2+4, c5 = c2+6 */\
		"movaps		    (%%rbx)	,%%xmm4 	\n\t		movaps		    (%%rcx)	,%%xmm2		\n\t"\
		"movaps		0x10(%%rbx)	,%%xmm5 	\n\t		movaps		0x10(%%rcx)	,%%xmm3		\n\t"\
		"movaps		%%xmm13		,%%xmm14	\n\t		movaps		%%xmm10		,%%xmm6		\n\t"\
		"movaps		%%xmm9 		,%%xmm15	\n\t		movaps		%%xmm11		,%%xmm7		\n\t"\
		"addpd		%%xmm4 	,%%xmm13		\n\t		addpd		%%xmm3	,%%xmm10		\n\t"\
		"subpd		%%xmm5 	,%%xmm9 		\n\t		subpd		%%xmm2	,%%xmm11		\n\t"\
		"subpd		%%xmm4 	,%%xmm14		\n\t"	/*	movaps		    (%%rsi)	,%%xmm0	Need xmm0 to replace xmm9  below */\
		"addpd		%%xmm5 	,%%xmm15		\n\t"	/*	movaps		0x10(%%rsi)	,%%xmm1	Need xmm1 to replace xmm13 below */\
		"movaps		%%xmm14		,%%xmm4 	\n\t		subpd		%%xmm3	,%%xmm6		\n\t"\
		"movaps		%%xmm15		,%%xmm5 	\n\t		addpd		%%xmm2	,%%xmm7		\n\t"/* xmm2,3 free */\
		"movaps		%%xmm13		,%%xmm14	\n\t		movaps		%%xmm10		,%%xmm0 		\n\t"\
		"movaps		%%xmm9 		,%%xmm15	\n\t		movaps		%%xmm11		,%%xmm1		\n\t"\
		"mulpd		0x40(%%rsi)	,%%xmm13	\n\t		mulpd			(%%rsi)	,%%xmm10		\n\t"\
		"mulpd		0x40(%%rsi)	,%%xmm9 	\n\t		mulpd			(%%rsi)	,%%xmm11		\n\t"\
		"mulpd		0x50(%%rsi)	,%%xmm14	\n\t		mulpd		0x10(%%rsi)	,%%xmm0 		\n\t"\
		"mulpd		0x50(%%rsi)	,%%xmm15	\n\t		mulpd		0x10(%%rsi)	,%%xmm1		\n\t"\
		"subpd		%%xmm14		,%%xmm9 	\n\t		subpd		%%xmm0 		,%%xmm11		\n\t"\
		"addpd		%%xmm15		,%%xmm13	\n\t		addpd		%%xmm1		,%%xmm10		\n\t"\
		"movaps		%%xmm9 		,0x10(%%r11)\n\t		movaps		0x20(%%rsi)	,%%xmm2		\n\t"\
		"movaps		%%xmm13		,    (%%r11)\n\t		movaps		0x30(%%rsi)	,%%xmm3		\n\t"\
		"movaps		%%xmm4 	,%%xmm13		\n\t		movaps		%%xmm11		,0x10(%%r12)	\n\t"\
		"movaps		%%xmm5 	,%%xmm9 		\n\t		movaps		%%xmm10		,    (%%r12)	\n\t"\
		"movaps		%%xmm13		,%%xmm14	\n\t		movaps		%%xmm6		,%%xmm0 		\n\t"\
		"movaps		%%xmm9 		,%%xmm15	\n\t		movaps		%%xmm7		,%%xmm1		\n\t"\
		"mulpd		0x60(%%rsi)	,%%xmm13	\n\t		mulpd		%%xmm2	,%%xmm6		\n\t"\
		"mulpd		0x60(%%rsi)	,%%xmm9 	\n\t		mulpd		%%xmm2	,%%xmm7		\n\t"\
		"mulpd		0x70(%%rsi)	,%%xmm14	\n\t		mulpd		%%xmm3	,%%xmm0 		\n\t"\
		"mulpd		0x70(%%rsi)	,%%xmm15	\n\t		mulpd		%%xmm3	,%%xmm1		\n\t"\
		"subpd		%%xmm14		,%%xmm9 	\n\t		subpd		%%xmm0 		,%%xmm7		\n\t"\
		"addpd		%%xmm15		,%%xmm13	\n\t		addpd		%%xmm1		,%%xmm6		\n\t"\
		"movaps		%%xmm9 		,0x10(%%rbx)\n\t		movaps		%%xmm7		,0x10(%%rcx)	\n\t"\
		"movaps		%%xmm13		,    (%%rbx)\n\t		movaps		%%xmm6		,    (%%rcx)	\n\t"\
		/* Outputs 3,7: All SIMD regs except for xmm8,12 (containing one of the 2 complex data being butterflied) free: */\
		"movq		%[c3]		,%%rsi		\n\t"/* c7 = c3+2 */\
		"movaps		    (%%rdx)	,%%xmm4		\n\t"\
		"movaps		0x10(%%rdx)	,%%xmm5		\n\t"\
		"movaps		%%xmm8 		,%%xmm14	\n\t"\
		"movaps		%%xmm12		,%%xmm15	\n\t"\
		"subpd		%%xmm5	,%%xmm8 		\n\t		subpd		%%xmm4	,%%xmm12		\n\t"\
		"movaps		    (%%rsi)	,%%xmm0		\n\t"\
		"movaps		0x10(%%rsi)	,%%xmm1		\n\t"\
		"addpd		%%xmm5	,%%xmm14		\n\t		movaps		0x20(%%rsi)	,%%xmm2		\n\t"\
		"addpd		%%xmm4	,%%xmm15		\n\t		movaps		0x30(%%rsi)	,%%xmm3		\n\t"\
		"movaps		%%xmm8 		,%%xmm9 	\n\t		movaps		%%xmm14		,%%xmm6 	\n\t"\
		"movaps		%%xmm12		,%%xmm13	\n\t		movaps		%%xmm15		,%%xmm7		\n\t"\
		"mulpd		%%xmm0	,%%xmm8 		\n\t		mulpd		%%xmm2	,%%xmm14		\n\t"\
		"mulpd		%%xmm0	,%%xmm12		\n\t		mulpd		%%xmm2	,%%xmm15		\n\t"\
		"mulpd		%%xmm1	,%%xmm9 		\n\t		mulpd		%%xmm3	,%%xmm6 		\n\t"\
		"mulpd		%%xmm1	,%%xmm13		\n\t		mulpd		%%xmm3	,%%xmm7			\n\t"\
		"subpd		%%xmm9 		,%%xmm12	\n\t		subpd		%%xmm6 		,%%xmm15	\n\t"\
		"addpd		%%xmm13		,%%xmm8 	\n\t		addpd		%%xmm7		,%%xmm14	\n\t"\
		"movaps		%%xmm12		,0x10(%%r13)\n\t		movaps		%%xmm15		,0x10(%%rdx)\n\t"\
		"movaps		%%xmm8 		,    (%%r13)\n\t		movaps		%%xmm14		,    (%%rdx)\n\t"\
		:					/* outputs: none */\
		: [add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[add1] "m" (Xadd1)\
		 ,[add2] "m" (Xadd2)\
		 ,[add3] "m" (Xadd3)\
		 ,[add4] "m" (Xadd4)\
		 ,[add5] "m" (Xadd5)\
		 ,[add6] "m" (Xadd6)\
		 ,[add7] "m" (Xadd7)\
		 ,[isrt2] "m" (Xisrt2)\
		 ,[c1] "m" (Xc1)\
		 ,[c2] "m" (Xc2)\
		 ,[c3] "m" (Xc3)\
		 ,[c4] "m" (Xc4)\
		 ,[c5] "m" (Xc5)\
		 ,[c6] "m" (Xc6)\
		 ,[c7] "m" (Xc7)\
		: "cc","memory","rax","rbx","rcx","rdx","rsi","r10","r11","r12","r13","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
	);\
	}

  #endif	// USE_64BIT_ASM_STYLE ?

#endif	// AVX / SSE2 toggle

#endif	/* radix8_dif_dit_pass_gcc_h_included */

