/*******************************************************************************
*                                                                              *
*   (C) 1997-2016 by Ernst W. Mayer.                                           *
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
#ifndef radix16_ditN_cy_dif1_gcc_h_included
#define radix16_ditN_cy_dif1_gcc_h_included

#ifdef USE_AVX512

	#define SSE2_RADIX16_DIT_NOTWIDDLE(Xadd0,Xp1,Xp2,Xp3,Xp4,Xr1,Xr3,Xr5,Xr7,Xr9,Xr11,Xr17,Xr25,Xisrt2,Xcc0)\
	{\
	__asm__ volatile (\
		"movq	%[__add0],%%rax		\n\t"\
		"movslq	%[__p1],%%rbx		\n\t"\
		"movslq	%[__p2],%%rcx		\n\t"\
		"movslq	%[__p3],%%rdx		\n\t"\
		"movslq	%[__p4],%%rdi		\n\t"\
		"shlq	$3,%%rbx			\n\t"\
		"shlq	$3,%%rcx			\n\t"\
		"shlq	$3,%%rdx			\n\t"\
		"shlq	$3,%%rdi			\n\t"\
		"addq	%%rax,%%rbx			\n\t"\
		"addq	%%rax,%%rcx			\n\t"\
		"addq	%%rax,%%rdx			\n\t"\
		/* SSE2_RADIX4_DIT_0TWIDDLE_B(r1 ): */\
		"vmovaps	     (%%rax),%%zmm2	\n\t"\
		"vmovaps	     (%%rcx),%%zmm6	\n\t"\
		"vmovaps	0x040(%%rax),%%zmm3	\n\t"\
		"vmovaps	0x040(%%rcx),%%zmm7	\n\t"\
		"vmovaps	     (%%rbx),%%zmm0	\n\t"\
		"vmovaps	     (%%rdx),%%zmm4	\n\t"\
		"vmovaps	0x040(%%rbx),%%zmm1	\n\t"\
		"vmovaps	0x040(%%rdx),%%zmm5	\n\t"\
		"vsubpd	%%zmm0,%%zmm2,%%zmm2		\n\t"\
		"vsubpd	%%zmm4,%%zmm6,%%zmm6		\n\t"\
		"vsubpd	%%zmm1,%%zmm3,%%zmm3		\n\t"\
		"vsubpd	%%zmm5,%%zmm7,%%zmm7		\n\t"\
		"vaddpd	%%zmm0,%%zmm0,%%zmm0		\n\t"\
		"vaddpd	%%zmm4,%%zmm4,%%zmm4		\n\t"\
		"vaddpd	%%zmm1,%%zmm1,%%zmm1		\n\t"\
		"vaddpd	%%zmm5,%%zmm5,%%zmm5		\n\t"\
		"vaddpd	%%zmm2,%%zmm0,%%zmm0		\n\t"\
		"vaddpd	%%zmm6,%%zmm4,%%zmm4		\n\t"\
		"vaddpd	%%zmm3,%%zmm1,%%zmm1		\n\t"\
		"vaddpd	%%zmm7,%%zmm5,%%zmm5		\n\t"\
		"vsubpd	%%zmm4,%%zmm0,%%zmm0		\n\t"\
		"vsubpd	%%zmm7,%%zmm2,%%zmm2		\n\t"\
		"vsubpd	%%zmm5,%%zmm1,%%zmm1		\n\t"\
		"vsubpd	%%zmm6,%%zmm3,%%zmm3		\n\t"\
		"movq	%[__r1], %%rsi		\n\t"\
		"vmovaps	%%zmm0,0x100(%%rsi)	\n\t"\
		"vmovaps	%%zmm2,0x180(%%rsi)	\n\t"\
		"vmovaps	%%zmm1,0x140(%%rsi)	\n\t"\
		"vmovaps	%%zmm3,0x0c0(%%rsi)	\n\t"\
		"vaddpd	%%zmm4,%%zmm4,%%zmm4		\n\t"\
		"vaddpd	%%zmm7,%%zmm7,%%zmm7		\n\t"\
		"vaddpd	%%zmm5,%%zmm5,%%zmm5		\n\t"\
		"vaddpd	%%zmm6,%%zmm6,%%zmm6		\n\t"\
		"vaddpd	%%zmm0,%%zmm4,%%zmm4		\n\t"\
		"vaddpd	%%zmm2,%%zmm7,%%zmm7		\n\t"\
		"vaddpd	%%zmm1,%%zmm5,%%zmm5		\n\t"\
		"vaddpd	%%zmm3,%%zmm6,%%zmm6		\n\t"\
		"vmovaps	%%zmm4,     (%%rsi)	\n\t"\
		"vmovaps	%%zmm7,0x080(%%rsi)	\n\t"\
		"vmovaps	%%zmm5,0x040(%%rsi)	\n\t"\
		"vmovaps	%%zmm6,0x1c0(%%rsi)	\n\t"\
		"addq	%%rdi,%%rax			\n\t"\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	%%rdi,%%rcx			\n\t"\
		"addq	%%rdi,%%rdx			\n\t"\
		/* SSE2_RADIX4_DIT_0TWIDDLE_B(r9 ): */\
		"vmovaps	     (%%rax),%%zmm2	\n\t"\
		"vmovaps	     (%%rcx),%%zmm6	\n\t"\
		"vmovaps	0x040(%%rax),%%zmm3	\n\t"\
		"vmovaps	0x040(%%rcx),%%zmm7	\n\t"\
		"vmovaps	     (%%rbx),%%zmm0	\n\t"\
		"vmovaps	     (%%rdx),%%zmm4	\n\t"\
		"vmovaps	0x040(%%rbx),%%zmm1	\n\t"\
		"vmovaps	0x040(%%rdx),%%zmm5	\n\t"\
		"vsubpd	%%zmm0,%%zmm2,%%zmm2		\n\t"\
		"vsubpd	%%zmm4,%%zmm6,%%zmm6		\n\t"\
		"vsubpd	%%zmm1,%%zmm3,%%zmm3		\n\t"\
		"vsubpd	%%zmm5,%%zmm7,%%zmm7		\n\t"\
		"vaddpd	%%zmm0,%%zmm0,%%zmm0		\n\t"\
		"vaddpd	%%zmm4,%%zmm4,%%zmm4		\n\t"\
		"vaddpd	%%zmm1,%%zmm1,%%zmm1		\n\t"\
		"vaddpd	%%zmm5,%%zmm5,%%zmm5		\n\t"\
		"vaddpd	%%zmm2,%%zmm0,%%zmm0		\n\t"\
		"vaddpd	%%zmm6,%%zmm4,%%zmm4		\n\t"\
		"vaddpd	%%zmm3,%%zmm1,%%zmm1		\n\t"\
		"vaddpd	%%zmm7,%%zmm5,%%zmm5		\n\t"\
		"vsubpd	%%zmm4,%%zmm0,%%zmm0		\n\t"\
		"vsubpd	%%zmm7,%%zmm2,%%zmm2		\n\t"\
		"vsubpd	%%zmm5,%%zmm1,%%zmm1		\n\t"\
		"vsubpd	%%zmm6,%%zmm3,%%zmm3		\n\t"\
		"movq	%[__r9], %%rsi		\n\t"\
		"vmovaps	%%zmm0,0x100(%%rsi)	\n\t"\
		"vmovaps	%%zmm2,0x180(%%rsi)	\n\t"\
		"vmovaps	%%zmm1,0x140(%%rsi)	\n\t"\
		"vmovaps	%%zmm3,0x0c0(%%rsi)	\n\t"\
		"vaddpd	%%zmm4,%%zmm4,%%zmm4		\n\t"\
		"vaddpd	%%zmm7,%%zmm7,%%zmm7		\n\t"\
		"vaddpd	%%zmm5,%%zmm5,%%zmm5		\n\t"\
		"vaddpd	%%zmm6,%%zmm6,%%zmm6		\n\t"\
		"vaddpd	%%zmm0,%%zmm4,%%zmm4		\n\t"\
		"vaddpd	%%zmm2,%%zmm7,%%zmm7		\n\t"\
		"vaddpd	%%zmm1,%%zmm5,%%zmm5		\n\t"\
		"vaddpd	%%zmm3,%%zmm6,%%zmm6		\n\t"\
		"vmovaps	%%zmm4,     (%%rsi)	\n\t"\
		"vmovaps	%%zmm7,0x080(%%rsi)	\n\t"\
		"vmovaps	%%zmm5,0x040(%%rsi)	\n\t"\
		"vmovaps	%%zmm6,0x1c0(%%rsi)	\n\t"\
		"addq	%%rdi,%%rax			\n\t"\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	%%rdi,%%rcx			\n\t"\
		"addq	%%rdi,%%rdx			\n\t"\
		/* SSE2_RADIX4_DIT_0TWIDDLE_B(r17): */\
		"vmovaps	     (%%rax),%%zmm2	\n\t"\
		"vmovaps	     (%%rcx),%%zmm6	\n\t"\
		"vmovaps	0x040(%%rax),%%zmm3	\n\t"\
		"vmovaps	0x040(%%rcx),%%zmm7	\n\t"\
		"vmovaps	     (%%rbx),%%zmm0	\n\t"\
		"vmovaps	     (%%rdx),%%zmm4	\n\t"\
		"vmovaps	0x040(%%rbx),%%zmm1	\n\t"\
		"vmovaps	0x040(%%rdx),%%zmm5	\n\t"\
		"vsubpd	%%zmm0,%%zmm2,%%zmm2		\n\t"\
		"vsubpd	%%zmm4,%%zmm6,%%zmm6		\n\t"\
		"vsubpd	%%zmm1,%%zmm3,%%zmm3		\n\t"\
		"vsubpd	%%zmm5,%%zmm7,%%zmm7		\n\t"\
		"vaddpd	%%zmm0,%%zmm0,%%zmm0		\n\t"\
		"vaddpd	%%zmm4,%%zmm4,%%zmm4		\n\t"\
		"vaddpd	%%zmm1,%%zmm1,%%zmm1		\n\t"\
		"vaddpd	%%zmm5,%%zmm5,%%zmm5		\n\t"\
		"vaddpd	%%zmm2,%%zmm0,%%zmm0		\n\t"\
		"vaddpd	%%zmm6,%%zmm4,%%zmm4		\n\t"\
		"vaddpd	%%zmm3,%%zmm1,%%zmm1		\n\t"\
		"vaddpd	%%zmm7,%%zmm5,%%zmm5		\n\t"\
		"vsubpd	%%zmm4,%%zmm0,%%zmm0		\n\t"\
		"vsubpd	%%zmm7,%%zmm2,%%zmm2		\n\t"\
		"vsubpd	%%zmm5,%%zmm1,%%zmm1		\n\t"\
		"vsubpd	%%zmm6,%%zmm3,%%zmm3		\n\t"\
		"movq	%[__r17], %%rsi		\n\t"\
		"vmovaps	%%zmm0,0x100(%%rsi)	\n\t"\
		"vmovaps	%%zmm2,0x180(%%rsi)	\n\t"\
		"vmovaps	%%zmm1,0x140(%%rsi)	\n\t"\
		"vmovaps	%%zmm3,0x0c0(%%rsi)	\n\t"\
		"vaddpd	%%zmm4,%%zmm4,%%zmm4		\n\t"\
		"vaddpd	%%zmm7,%%zmm7,%%zmm7		\n\t"\
		"vaddpd	%%zmm5,%%zmm5,%%zmm5		\n\t"\
		"vaddpd	%%zmm6,%%zmm6,%%zmm6		\n\t"\
		"vaddpd	%%zmm0,%%zmm4,%%zmm4		\n\t"\
		"vaddpd	%%zmm2,%%zmm7,%%zmm7		\n\t"\
		"vaddpd	%%zmm1,%%zmm5,%%zmm5		\n\t"\
		"vaddpd	%%zmm3,%%zmm6,%%zmm6		\n\t"\
		"vmovaps	%%zmm4,     (%%rsi)	\n\t"\
		"vmovaps	%%zmm7,0x080(%%rsi)	\n\t"\
		"vmovaps	%%zmm5,0x040(%%rsi)	\n\t"\
		"vmovaps	%%zmm6,0x1c0(%%rsi)	\n\t"\
		"addq	%%rdi,%%rax			\n\t"\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	%%rdi,%%rcx			\n\t"\
		"addq	%%rdi,%%rdx			\n\t"\
		/* SSE2_RADIX4_DIT_0TWIDDLE_B(r25): */\
		"vmovaps	     (%%rax),%%zmm2	\n\t"\
		"vmovaps	     (%%rcx),%%zmm6	\n\t"\
		"vmovaps	0x040(%%rax),%%zmm3	\n\t"\
		"vmovaps	0x040(%%rcx),%%zmm7	\n\t"\
		"vmovaps	     (%%rbx),%%zmm0	\n\t"\
		"vmovaps	     (%%rdx),%%zmm4	\n\t"\
		"vmovaps	0x040(%%rbx),%%zmm1	\n\t"\
		"vmovaps	0x040(%%rdx),%%zmm5	\n\t"\
		"vsubpd	%%zmm0,%%zmm2,%%zmm2		\n\t"\
		"vsubpd	%%zmm4,%%zmm6,%%zmm6		\n\t"\
		"vsubpd	%%zmm1,%%zmm3,%%zmm3		\n\t"\
		"vsubpd	%%zmm5,%%zmm7,%%zmm7		\n\t"\
		"vaddpd	%%zmm0,%%zmm0,%%zmm0		\n\t"\
		"vaddpd	%%zmm4,%%zmm4,%%zmm4		\n\t"\
		"vaddpd	%%zmm1,%%zmm1,%%zmm1		\n\t"\
		"vaddpd	%%zmm5,%%zmm5,%%zmm5		\n\t"\
		"vaddpd	%%zmm2,%%zmm0,%%zmm0		\n\t"\
		"vaddpd	%%zmm6,%%zmm4,%%zmm4		\n\t"\
		"vaddpd	%%zmm3,%%zmm1,%%zmm1		\n\t"\
		"vaddpd	%%zmm7,%%zmm5,%%zmm5		\n\t"\
		"vsubpd	%%zmm4,%%zmm0,%%zmm0		\n\t"\
		"vsubpd	%%zmm7,%%zmm2,%%zmm2		\n\t"\
		"vsubpd	%%zmm5,%%zmm1,%%zmm1		\n\t"\
		"vsubpd	%%zmm6,%%zmm3,%%zmm3		\n\t"\
		"movq	%[__r25], %%rsi		\n\t"\
		"vmovaps	%%zmm0,0x100(%%rsi)	\n\t"\
		"vmovaps	%%zmm2,0x180(%%rsi)	\n\t"\
		"vmovaps	%%zmm1,0x140(%%rsi)	\n\t"\
		"vmovaps	%%zmm3,0x0c0(%%rsi)	\n\t"\
		"vaddpd	%%zmm4,%%zmm4,%%zmm4		\n\t"\
		"vaddpd	%%zmm7,%%zmm7,%%zmm7		\n\t"\
		"vaddpd	%%zmm5,%%zmm5,%%zmm5		\n\t"\
		"vaddpd	%%zmm6,%%zmm6,%%zmm6		\n\t"\
		"vaddpd	%%zmm0,%%zmm4,%%zmm4		\n\t"\
		"vaddpd	%%zmm2,%%zmm7,%%zmm7		\n\t"\
		"vaddpd	%%zmm1,%%zmm5,%%zmm5		\n\t"\
		"vaddpd	%%zmm3,%%zmm6,%%zmm6		\n\t"\
		"vmovaps	%%zmm4,     (%%rsi)	\n\t"\
		"vmovaps	%%zmm7,0x080(%%rsi)	\n\t"\
		"vmovaps	%%zmm5,0x040(%%rsi)	\n\t"\
		"vmovaps	%%zmm6,0x1c0(%%rsi)	\n\t"\
	/*** Now do four pass-2 4-DFTs, inputs from local stor, outputs back to main array: ***/\
		"movq	%[__r1],%%rax		\n\t"\
		"vmovaps	     (%%rax),%%zmm0	\n\t"\
		"vmovaps	0x040(%%rax),%%zmm1	\n\t"\
		"vmovaps	0x200(%%rax),%%zmm2	\n\t"\
		"vmovaps	0x240(%%rax),%%zmm3	\n\t"\
		"vsubpd	0x200(%%rax),%%zmm0,%%zmm0	\n\t"\
		"vsubpd	0x240(%%rax),%%zmm1,%%zmm1	\n\t"\
		"vaddpd	     (%%rax),%%zmm2,%%zmm2	\n\t"\
		"vaddpd	0x040(%%rax),%%zmm3,%%zmm3	\n\t"\
		"vmovaps	0x400(%%rax),%%zmm4	\n\t"\
		"vmovaps	0x440(%%rax),%%zmm5	\n\t"\
		"vmovaps	0x600(%%rax),%%zmm6	\n\t"\
		"vmovaps	0x640(%%rax),%%zmm7	\n\t"\
		"vsubpd	0x600(%%rax),%%zmm4,%%zmm4	\n\t"\
		"vsubpd	0x640(%%rax),%%zmm5,%%zmm5	\n\t"\
		"vaddpd	0x400(%%rax),%%zmm6,%%zmm6	\n\t"\
		"vaddpd	0x440(%%rax),%%zmm7,%%zmm7	\n\t"\
		"vsubpd	      %%zmm6,%%zmm2,%%zmm2		\n\t"\
		"vsubpd	      %%zmm7,%%zmm3,%%zmm3		\n\t"\
		"vmovaps	%%zmm2,0x400(%%rax)	\n\t"\
		"vmovaps	%%zmm3,0x440(%%rax)	\n\t"\
		"vaddpd	      %%zmm6,%%zmm6,%%zmm6		\n\t"\
		"vaddpd	      %%zmm7,%%zmm7,%%zmm7		\n\t"\
		"vaddpd	      %%zmm2,%%zmm6,%%zmm6		\n\t"\
		"vaddpd	      %%zmm3,%%zmm7,%%zmm7		\n\t"\
		"vmovaps	%%zmm6,     (%%rax)	\n\t"\
		"vmovaps	%%zmm7,0x040(%%rax)	\n\t"\
		"vsubpd	      %%zmm5,%%zmm0,%%zmm0		\n\t"\
		"vsubpd	      %%zmm4,%%zmm1,%%zmm1		\n\t"\
		"vmovaps	%%zmm0,0x600(%%rax)	\n\t"\
		"vmovaps	%%zmm1,0x240(%%rax)	\n\t"\
		"vaddpd	      %%zmm5,%%zmm5,%%zmm5		\n\t"\
		"vaddpd	      %%zmm4,%%zmm4,%%zmm4		\n\t"\
		"vaddpd	      %%zmm0,%%zmm5,%%zmm5		\n\t"\
		"vaddpd	      %%zmm1,%%zmm4,%%zmm4		\n\t"\
		"vmovaps	%%zmm5,0x200(%%rax)	\n\t"\
		"vmovaps	%%zmm4,0x640(%%rax)	\n\t"\
		"movq	%[__r5],%%rax		\n\t"\
		"movq	%[__isrt2],%%rbx	\n\t"\
		"vmovaps	(%%rbx),%%zmm2		\n\t"\
		"vmovaps	0x400(%%rax),%%zmm4	\n\t"\
		"vmovaps	0x440(%%rax),%%zmm5	\n\t"\
		"vmovaps	0x600(%%rax),%%zmm0	\n\t"\
		"vmovaps	0x640(%%rax),%%zmm1	\n\t"\
		"vaddpd	0x440(%%rax),%%zmm4,%%zmm4	\n\t"\
		"vsubpd	0x400(%%rax),%%zmm5,%%zmm5	\n\t"\
		"vsubpd	0x640(%%rax),%%zmm0,%%zmm0	\n\t"\
		"vaddpd	0x600(%%rax),%%zmm1,%%zmm1	\n\t"\
		"vmulpd	      %%zmm2,%%zmm4,%%zmm4		\n\t"\
		"vmulpd	      %%zmm2,%%zmm5,%%zmm5		\n\t"\
		"vmulpd	      %%zmm2,%%zmm0,%%zmm0		\n\t"\
		"vmulpd	      %%zmm2,%%zmm1,%%zmm1		\n\t"\
		"vmovaps	%%zmm4,%%zmm6		\n\t"\
		"vmovaps	%%zmm5,%%zmm7		\n\t"\
		"vsubpd	      %%zmm0,%%zmm4,%%zmm4		\n\t"\
		"vsubpd	      %%zmm1,%%zmm5,%%zmm5		\n\t"\
		"vaddpd	      %%zmm0,%%zmm6,%%zmm6		\n\t"\
		"vaddpd	      %%zmm1,%%zmm7,%%zmm7		\n\t"\
		"vmovaps	     (%%rax),%%zmm0	\n\t"\
		"vmovaps	0x040(%%rax),%%zmm1	\n\t"\
		"vmovaps	0x200(%%rax),%%zmm2	\n\t"\
		"vmovaps	0x240(%%rax),%%zmm3	\n\t"\
		"vsubpd	0x240(%%rax),%%zmm0,%%zmm0	\n\t"\
		"vsubpd	0x200(%%rax),%%zmm1,%%zmm1	\n\t"\
		"vaddpd	     (%%rax),%%zmm3,%%zmm3	\n\t"\
		"vaddpd	0x040(%%rax),%%zmm2,%%zmm2	\n\t"\
		"vsubpd	      %%zmm4,%%zmm3,%%zmm3		\n\t"\
		"vsubpd	      %%zmm5,%%zmm1,%%zmm1		\n\t"\
		"vmovaps	%%zmm3,0x400(%%rax)	\n\t"\
		"vmovaps	%%zmm1,0x440(%%rax)	\n\t"\
		"vaddpd	      %%zmm4,%%zmm4,%%zmm4		\n\t"\
		"vaddpd	      %%zmm5,%%zmm5,%%zmm5		\n\t"\
		"vaddpd	      %%zmm3,%%zmm4,%%zmm4		\n\t"\
		"vaddpd	      %%zmm1,%%zmm5,%%zmm5		\n\t"\
		"vmovaps	%%zmm4,     (%%rax)	\n\t"\
		"vmovaps	%%zmm5,0x040(%%rax)	\n\t"\
		"vsubpd	      %%zmm7,%%zmm0,%%zmm0		\n\t"\
		"vsubpd	      %%zmm6,%%zmm2,%%zmm2		\n\t"\
		"vmovaps	%%zmm0,0x600(%%rax)	\n\t"\
		"vmovaps	%%zmm2,0x240(%%rax)	\n\t"\
		"vaddpd	      %%zmm7,%%zmm7,%%zmm7		\n\t"\
		"vaddpd	      %%zmm6,%%zmm6,%%zmm6		\n\t"\
		"vaddpd	      %%zmm0,%%zmm7,%%zmm7		\n\t"\
		"vaddpd	      %%zmm2,%%zmm6,%%zmm6		\n\t"\
		"vmovaps	%%zmm7,0x200(%%rax)	\n\t"\
		"vmovaps	%%zmm6,0x640(%%rax)	\n\t"\
		"movq	%[__r3],%%rax		\n\t"\
		"movq	%[__isrt2],%%rbx	\n\t"\
		"movq	%[__cc0],%%rcx		\n\t"\
		"vmovaps	0x600(%%rax),%%zmm0	\n\t"\
		"vmovaps	0x640(%%rax),%%zmm1	\n\t"\
		"vmovaps	0x600(%%rax),%%zmm2	\n\t"\
		"vmovaps	0x640(%%rax),%%zmm3	\n\t"\
		"vmulpd	0x040(%%rcx),%%zmm0,%%zmm0	\n\t"\
		"vmulpd	0x040(%%rcx),%%zmm1,%%zmm1	\n\t"\
		"vmulpd	     (%%rcx),%%zmm2,%%zmm2	\n\t"\
		"vmulpd	     (%%rcx),%%zmm3,%%zmm3	\n\t"\
		"vsubpd	      %%zmm2,%%zmm1,%%zmm1		\n\t"\
		"vaddpd	      %%zmm3,%%zmm0,%%zmm0		\n\t"\
		"vmovaps	0x400(%%rax),%%zmm4	\n\t"\
		"vmovaps	0x440(%%rax),%%zmm5	\n\t"\
		"vmovaps	0x400(%%rax),%%zmm6	\n\t"\
		"vmovaps	0x440(%%rax),%%zmm7	\n\t"\
		"vmulpd	     (%%rcx),%%zmm4,%%zmm4	\n\t"\
		"vmulpd	     (%%rcx),%%zmm5,%%zmm5	\n\t"\
		"vmulpd	0x040(%%rcx),%%zmm6,%%zmm6	\n\t"\
		"vmulpd	0x040(%%rcx),%%zmm7,%%zmm7	\n\t"\
		"vsubpd	      %%zmm6,%%zmm5,%%zmm5		\n\t"\
		"vaddpd	      %%zmm7,%%zmm4,%%zmm4		\n\t"\
		"vmovaps	%%zmm5,%%zmm7		\n\t"\
		"vmovaps	%%zmm4,%%zmm6		\n\t"\
		"vaddpd	      %%zmm0,%%zmm4,%%zmm4		\n\t"\
		"vaddpd	      %%zmm1,%%zmm5,%%zmm5		\n\t"\
		"vsubpd	      %%zmm0,%%zmm6,%%zmm6		\n\t"\
		"vsubpd	      %%zmm1,%%zmm7,%%zmm7		\n\t"\
		"vmovaps	0x200(%%rax),%%zmm2	\n\t"\
		"vmovaps	0x240(%%rax),%%zmm3	\n\t"\
		"vmovaps	     (%%rax),%%zmm0	\n\t"\
		"vmovaps	0x040(%%rax),%%zmm1	\n\t"\
		"vaddpd	0x240(%%rax),%%zmm2,%%zmm2	\n\t"\
		"vsubpd	0x200(%%rax),%%zmm3,%%zmm3	\n\t"\
		"vmulpd	     (%%rbx),%%zmm2,%%zmm2	\n\t"\
		"vmulpd	     (%%rbx),%%zmm3,%%zmm3	\n\t"\
		"vsubpd	      %%zmm2,%%zmm0,%%zmm0		\n\t"\
		"vsubpd	      %%zmm3,%%zmm1,%%zmm1		\n\t"\
		"vaddpd	      %%zmm2,%%zmm2,%%zmm2		\n\t"\
		"vaddpd	      %%zmm3,%%zmm3,%%zmm3		\n\t"\
		"vaddpd	      %%zmm0,%%zmm2,%%zmm2		\n\t"\
		"vaddpd	      %%zmm1,%%zmm3,%%zmm3		\n\t"\
		"vsubpd	      %%zmm4,%%zmm2,%%zmm2		\n\t"\
		"vsubpd	      %%zmm5,%%zmm3,%%zmm3		\n\t"\
		"vmovaps	%%zmm2,0x400(%%rax)	\n\t"\
		"vmovaps	%%zmm3,0x440(%%rax)	\n\t"\
		"vaddpd	      %%zmm4,%%zmm4,%%zmm4		\n\t"\
		"vaddpd	      %%zmm5,%%zmm5,%%zmm5		\n\t"\
		"vaddpd	      %%zmm2,%%zmm4,%%zmm4		\n\t"\
		"vaddpd	      %%zmm3,%%zmm5,%%zmm5		\n\t"\
		"vmovaps	%%zmm4,     (%%rax)	\n\t"\
		"vmovaps	%%zmm5,0x040(%%rax)	\n\t"\
		"vsubpd	      %%zmm7,%%zmm0,%%zmm0		\n\t"\
		"vsubpd	      %%zmm6,%%zmm1,%%zmm1		\n\t"\
		"vmovaps	%%zmm0,0x600(%%rax)	\n\t"\
		"vmovaps	%%zmm1,0x240(%%rax)	\n\t"\
		"vaddpd	      %%zmm7,%%zmm7,%%zmm7		\n\t"\
		"vaddpd	      %%zmm6,%%zmm6,%%zmm6		\n\t"\
		"vaddpd	      %%zmm0,%%zmm7,%%zmm7		\n\t"\
		"vaddpd	      %%zmm1,%%zmm6,%%zmm6		\n\t"\
		"vmovaps	%%zmm7,0x200(%%rax)	\n\t"\
		"vmovaps	%%zmm6,0x640(%%rax)	\n\t"\
		"movq	%[__r7],%%rax		\n\t"\
		"movq	%[__isrt2],%%rbx	\n\t"\
		"movq	%[__cc0],%%rcx		\n\t"\
		"vmovaps	0x600(%%rax),%%zmm0	\n\t"\
		"vmovaps	0x640(%%rax),%%zmm1	\n\t"\
		"vmovaps	0x600(%%rax),%%zmm2	\n\t"\
		"vmovaps	0x640(%%rax),%%zmm3	\n\t"\
		"vmulpd	     (%%rcx),%%zmm0,%%zmm0	\n\t"\
		"vmulpd	     (%%rcx),%%zmm1,%%zmm1	\n\t"\
		"vmulpd	0x040(%%rcx),%%zmm2,%%zmm2	\n\t"\
		"vmulpd	0x040(%%rcx),%%zmm3,%%zmm3	\n\t"\
		"vsubpd	      %%zmm2,%%zmm1,%%zmm1		\n\t"\
		"vaddpd	      %%zmm3,%%zmm0,%%zmm0		\n\t"\
		"vmovaps	0x400(%%rax),%%zmm4	\n\t"\
		"vmovaps	0x440(%%rax),%%zmm5	\n\t"\
		"vmovaps	0x400(%%rax),%%zmm6	\n\t"\
		"vmovaps	0x440(%%rax),%%zmm7	\n\t"\
		"vmulpd	0x040(%%rcx),%%zmm4,%%zmm4	\n\t"\
		"vmulpd	0x040(%%rcx),%%zmm5,%%zmm5	\n\t"\
		"vmulpd	     (%%rcx),%%zmm6,%%zmm6	\n\t"\
		"vmulpd	     (%%rcx),%%zmm7,%%zmm7	\n\t"\
		"vsubpd	      %%zmm6,%%zmm5,%%zmm5		\n\t"\
		"vaddpd	      %%zmm7,%%zmm4,%%zmm4		\n\t"\
		"vmovaps	%%zmm5,%%zmm7		\n\t"\
		"vmovaps	%%zmm4,%%zmm6		\n\t"\
		"vaddpd	      %%zmm0,%%zmm4,%%zmm4		\n\t"\
		"vaddpd	      %%zmm1,%%zmm5,%%zmm5		\n\t"\
		"vsubpd	      %%zmm0,%%zmm6,%%zmm6		\n\t"\
		"vsubpd	      %%zmm1,%%zmm7,%%zmm7		\n\t"\
		"vmovaps	0x200(%%rax),%%zmm2	\n\t"\
		"vmovaps	0x240(%%rax),%%zmm3	\n\t"\
		"vmovaps	     (%%rax),%%zmm0	\n\t"\
		"vmovaps	0x040(%%rax),%%zmm1	\n\t"\
		"vsubpd	0x240(%%rax),%%zmm2,%%zmm2	\n\t"\
		"vaddpd	0x200(%%rax),%%zmm3,%%zmm3	\n\t"\
		"vmulpd	     (%%rbx),%%zmm2,%%zmm2	\n\t"\
		"vmulpd	     (%%rbx),%%zmm3,%%zmm3	\n\t"\
		"vsubpd	      %%zmm2,%%zmm0,%%zmm0		\n\t"\
		"vsubpd	      %%zmm3,%%zmm1,%%zmm1		\n\t"\
		"vaddpd	      %%zmm2,%%zmm2,%%zmm2		\n\t"\
		"vaddpd	      %%zmm3,%%zmm3,%%zmm3		\n\t"\
		"vaddpd	      %%zmm0,%%zmm2,%%zmm2		\n\t"\
		"vaddpd	      %%zmm1,%%zmm3,%%zmm3		\n\t"\
		"vsubpd	      %%zmm6,%%zmm0,%%zmm0		\n\t"\
		"vsubpd	      %%zmm7,%%zmm1,%%zmm1		\n\t"\
		"vmovaps	%%zmm0,0x400(%%rax)	\n\t"\
		"vmovaps	%%zmm1,0x440(%%rax)	\n\t"\
		"vaddpd	      %%zmm6,%%zmm6,%%zmm6		\n\t"\
		"vaddpd	      %%zmm7,%%zmm7,%%zmm7		\n\t"\
		"vaddpd	      %%zmm0,%%zmm6,%%zmm6		\n\t"\
		"vaddpd	      %%zmm1,%%zmm7,%%zmm7		\n\t"\
		"vmovaps	%%zmm6,     (%%rax)	\n\t"\
		"vmovaps	%%zmm7,0x040(%%rax)	\n\t"\
		"vsubpd	      %%zmm5,%%zmm2,%%zmm2		\n\t"\
		"vsubpd	      %%zmm4,%%zmm3,%%zmm3		\n\t"\
		"vmovaps	%%zmm2,0x600(%%rax)	\n\t"\
		"vmovaps	%%zmm3,0x240(%%rax)	\n\t"\
		"vaddpd	      %%zmm5,%%zmm5,%%zmm5		\n\t"\
		"vaddpd	      %%zmm4,%%zmm4,%%zmm4		\n\t"\
		"vaddpd	      %%zmm2,%%zmm5,%%zmm5		\n\t"\
		"vaddpd	      %%zmm3,%%zmm4,%%zmm4		\n\t"\
		"vmovaps	%%zmm5,0x200(%%rax)	\n\t"\
		"vmovaps	%%zmm4,0x640(%%rax)	\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__p1] "m" (Xp1)\
		 ,[__p2] "m" (Xp2)\
		 ,[__p3] "m" (Xp3)\
		 ,[__p4] "m" (Xp4)\
		 ,[__r1] "m" (Xr1)\
		 ,[__r3] "m" (Xr3)\
		 ,[__r5] "m" (Xr5)\
		 ,[__r7] "m" (Xr7)\
		 ,[__r9] "m" (Xr9)\
		 ,[__r17] "m" (Xr17)\
		 ,[__r25] "m" (Xr25)\
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__cc0] "m" (Xcc0)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
	);\
	}


	#define SSE2_RADIX16_DIF_NOTWIDDLE(Xadd0,Xp1,Xp2,Xp3,Xp4,Xr1,Xr3,Xr5,Xr7,Xr9,Xr11,Xr13,Xr15,Xr17,Xr19,Xr21,Xr23,Xr25,Xr27,Xr29,Xr31,Xisrt2,Xcc0)\
	{\
	__asm__ volatile (\
		/* SSE2_RADIX4_DIF_IN_PLACE(r1 , r17, r9 , r25): */\
		"movq	%[__r1] ,%%rax		\n\t"\
		"movq	%[__r17],%%rbx		\n\t"\
		"movq	%[__r9] ,%%rcx		\n\t"\
		"movq	%[__r25],%%rdx		\n\t"\
		"vmovaps	     (%%rax),%%zmm0	\n\t"\
		"vmovaps	0x040(%%rax),%%zmm1	\n\t"\
		"vmovaps	     (%%rax),%%zmm2	\n\t"\
		"vmovaps	0x040(%%rax),%%zmm3	\n\t"\
		"vaddpd	     (%%rbx),%%zmm0,%%zmm0	\n\t"\
		"vaddpd	0x040(%%rbx),%%zmm1,%%zmm1	\n\t"\
		"vsubpd	     (%%rbx),%%zmm2,%%zmm2	\n\t"\
		"vsubpd	0x040(%%rbx),%%zmm3,%%zmm3	\n\t"\
		"vmovaps	     (%%rcx),%%zmm4	\n\t"\
		"vmovaps	0x040(%%rcx),%%zmm5	\n\t"\
		"vmovaps	     (%%rcx),%%zmm6	\n\t"\
		"vmovaps	0x040(%%rcx),%%zmm7	\n\t"\
		"vaddpd	     (%%rdx),%%zmm4,%%zmm4	\n\t"\
		"vaddpd	0x040(%%rdx),%%zmm5,%%zmm5	\n\t"\
		"vsubpd	     (%%rdx),%%zmm6,%%zmm6	\n\t"\
		"vsubpd	0x040(%%rdx),%%zmm7,%%zmm7	\n\t"\
		"vsubpd	      %%zmm4,%%zmm0,%%zmm0		\n\t"\
		"vsubpd	      %%zmm5,%%zmm1,%%zmm1		\n\t"\
		"vmovaps	%%zmm0,     (%%rbx)	\n\t"\
		"vmovaps	%%zmm1,0x040(%%rbx)	\n\t"\
		"vaddpd	      %%zmm4,%%zmm4,%%zmm4		\n\t"\
		"vaddpd	      %%zmm5,%%zmm5,%%zmm5		\n\t"\
		"vaddpd	      %%zmm0,%%zmm4,%%zmm4		\n\t"\
		"vaddpd	      %%zmm1,%%zmm5,%%zmm5		\n\t"\
		"vmovaps	%%zmm4,     (%%rax)	\n\t"\
		"vmovaps	%%zmm5,0x040(%%rax)	\n\t"\
		"vsubpd	      %%zmm7,%%zmm2,%%zmm2		\n\t"\
		"vsubpd	      %%zmm6,%%zmm3,%%zmm3		\n\t"\
		"vmovaps	%%zmm2,     (%%rcx)	\n\t"\
		"vmovaps	%%zmm3,0x040(%%rdx)	\n\t"\
		"vaddpd	      %%zmm7,%%zmm7,%%zmm7		\n\t"\
		"vaddpd	      %%zmm6,%%zmm6,%%zmm6		\n\t"\
		"vaddpd	      %%zmm2,%%zmm7,%%zmm7		\n\t"\
		"vaddpd	      %%zmm3,%%zmm6,%%zmm6		\n\t"\
		"vmovaps	%%zmm7,     (%%rdx)	\n\t"\
		"vmovaps	%%zmm6,0x040(%%rcx)	\n\t"\
		/* SSE2_RADIX4_DIF_IN_PLACE(r5 , r21, r13, r29): */\
		"movq	%[__r5] ,%%rax		\n\t"\
		"movq	%[__r21],%%rbx		\n\t"\
		"movq	%[__r13],%%rcx		\n\t"\
		"movq	%[__r29],%%rdx		\n\t"\
		"vmovaps	     (%%rax),%%zmm0	\n\t"\
		"vmovaps	0x040(%%rax),%%zmm1	\n\t"\
		"vmovaps	     (%%rax),%%zmm2	\n\t"\
		"vmovaps	0x040(%%rax),%%zmm3	\n\t"\
		"vaddpd	     (%%rbx),%%zmm0,%%zmm0	\n\t"\
		"vaddpd	0x040(%%rbx),%%zmm1,%%zmm1	\n\t"\
		"vsubpd	     (%%rbx),%%zmm2,%%zmm2	\n\t"\
		"vsubpd	0x040(%%rbx),%%zmm3,%%zmm3	\n\t"\
		"vmovaps	     (%%rcx),%%zmm4	\n\t"\
		"vmovaps	0x040(%%rcx),%%zmm5	\n\t"\
		"vmovaps	     (%%rcx),%%zmm6	\n\t"\
		"vmovaps	0x040(%%rcx),%%zmm7	\n\t"\
		"vaddpd	     (%%rdx),%%zmm4,%%zmm4	\n\t"\
		"vaddpd	0x040(%%rdx),%%zmm5,%%zmm5	\n\t"\
		"vsubpd	     (%%rdx),%%zmm6,%%zmm6	\n\t"\
		"vsubpd	0x040(%%rdx),%%zmm7,%%zmm7	\n\t"\
		"vsubpd	      %%zmm4,%%zmm0,%%zmm0		\n\t"\
		"vsubpd	      %%zmm5,%%zmm1,%%zmm1		\n\t"\
		"vmovaps	%%zmm0,     (%%rbx)	\n\t"\
		"vmovaps	%%zmm1,0x040(%%rbx)	\n\t"\
		"vaddpd	      %%zmm4,%%zmm4,%%zmm4		\n\t"\
		"vaddpd	      %%zmm5,%%zmm5,%%zmm5		\n\t"\
		"vaddpd	      %%zmm0,%%zmm4,%%zmm4		\n\t"\
		"vaddpd	      %%zmm1,%%zmm5,%%zmm5		\n\t"\
		"vmovaps	%%zmm4,     (%%rax)	\n\t"\
		"vmovaps	%%zmm5,0x040(%%rax)	\n\t"\
		"vsubpd	      %%zmm7,%%zmm2,%%zmm2		\n\t"\
		"vsubpd	      %%zmm6,%%zmm3,%%zmm3		\n\t"\
		"vmovaps	%%zmm2,     (%%rcx)	\n\t"\
		"vmovaps	%%zmm3,0x040(%%rdx)	\n\t"\
		"vaddpd	      %%zmm7,%%zmm7,%%zmm7		\n\t"\
		"vaddpd	      %%zmm6,%%zmm6,%%zmm6		\n\t"\
		"vaddpd	      %%zmm2,%%zmm7,%%zmm7		\n\t"\
		"vaddpd	      %%zmm3,%%zmm6,%%zmm6		\n\t"\
		"vmovaps	%%zmm7,     (%%rdx)	\n\t"\
		"vmovaps	%%zmm6,0x040(%%rcx)	\n\t"\
		/* SSE2_RADIX4_DIF_IN_PLACE(r3 , r19, r11, r27): */\
		"movq	%[__r3] ,%%rax		\n\t"\
		"movq	%[__r19],%%rbx		\n\t"\
		"movq	%[__r11],%%rcx		\n\t"\
		"movq	%[__r27],%%rdx		\n\t"\
		"vmovaps	     (%%rax),%%zmm0	\n\t"\
		"vmovaps	0x040(%%rax),%%zmm1	\n\t"\
		"vmovaps	     (%%rax),%%zmm2	\n\t"\
		"vmovaps	0x040(%%rax),%%zmm3	\n\t"\
		"vaddpd	     (%%rbx),%%zmm0,%%zmm0	\n\t"\
		"vaddpd	0x040(%%rbx),%%zmm1,%%zmm1	\n\t"\
		"vsubpd	     (%%rbx),%%zmm2,%%zmm2	\n\t"\
		"vsubpd	0x040(%%rbx),%%zmm3,%%zmm3	\n\t"\
		"vmovaps	     (%%rcx),%%zmm4	\n\t"\
		"vmovaps	0x040(%%rcx),%%zmm5	\n\t"\
		"vmovaps	     (%%rcx),%%zmm6	\n\t"\
		"vmovaps	0x040(%%rcx),%%zmm7	\n\t"\
		"vaddpd	     (%%rdx),%%zmm4,%%zmm4	\n\t"\
		"vaddpd	0x040(%%rdx),%%zmm5,%%zmm5	\n\t"\
		"vsubpd	     (%%rdx),%%zmm6,%%zmm6	\n\t"\
		"vsubpd	0x040(%%rdx),%%zmm7,%%zmm7	\n\t"\
		"vsubpd	      %%zmm4,%%zmm0,%%zmm0		\n\t"\
		"vsubpd	      %%zmm5,%%zmm1,%%zmm1		\n\t"\
		"vmovaps	%%zmm0,     (%%rbx)	\n\t"\
		"vmovaps	%%zmm1,0x040(%%rbx)	\n\t"\
		"vaddpd	      %%zmm4,%%zmm4,%%zmm4		\n\t"\
		"vaddpd	      %%zmm5,%%zmm5,%%zmm5		\n\t"\
		"vaddpd	      %%zmm0,%%zmm4,%%zmm4		\n\t"\
		"vaddpd	      %%zmm1,%%zmm5,%%zmm5		\n\t"\
		"vmovaps	%%zmm4,     (%%rax)	\n\t"\
		"vmovaps	%%zmm5,0x040(%%rax)	\n\t"\
		"vsubpd	      %%zmm7,%%zmm2,%%zmm2		\n\t"\
		"vsubpd	      %%zmm6,%%zmm3,%%zmm3		\n\t"\
		"vmovaps	%%zmm2,     (%%rcx)	\n\t"\
		"vmovaps	%%zmm3,0x040(%%rdx)	\n\t"\
		"vaddpd	      %%zmm7,%%zmm7,%%zmm7		\n\t"\
		"vaddpd	      %%zmm6,%%zmm6,%%zmm6		\n\t"\
		"vaddpd	      %%zmm2,%%zmm7,%%zmm7		\n\t"\
		"vaddpd	      %%zmm3,%%zmm6,%%zmm6		\n\t"\
		"vmovaps	%%zmm7,     (%%rdx)	\n\t"\
		"vmovaps	%%zmm6,0x040(%%rcx)	\n\t"\
		/* SSE2_RADIX4_DIF_IN_PLACE(r7 , r23, r15, r31): */\
		"movq	%[__r7] ,%%rax		\n\t"\
		"movq	%[__r23],%%rbx		\n\t"\
		"movq	%[__r15],%%rcx		\n\t"\
		"movq	%[__r31],%%rdx		\n\t"\
		"vmovaps	     (%%rax),%%zmm0	\n\t"\
		"vmovaps	0x040(%%rax),%%zmm1	\n\t"\
		"vmovaps	     (%%rax),%%zmm2	\n\t"\
		"vmovaps	0x040(%%rax),%%zmm3	\n\t"\
		"vaddpd	     (%%rbx),%%zmm0,%%zmm0	\n\t"\
		"vaddpd	0x040(%%rbx),%%zmm1,%%zmm1	\n\t"\
		"vsubpd	     (%%rbx),%%zmm2,%%zmm2	\n\t"\
		"vsubpd	0x040(%%rbx),%%zmm3,%%zmm3	\n\t"\
		"vmovaps	     (%%rcx),%%zmm4	\n\t"\
		"vmovaps	0x040(%%rcx),%%zmm5	\n\t"\
		"vmovaps	     (%%rcx),%%zmm6	\n\t"\
		"vmovaps	0x040(%%rcx),%%zmm7	\n\t"\
		"vaddpd	     (%%rdx),%%zmm4,%%zmm4	\n\t"\
		"vaddpd	0x040(%%rdx),%%zmm5,%%zmm5	\n\t"\
		"vsubpd	     (%%rdx),%%zmm6,%%zmm6	\n\t"\
		"vsubpd	0x040(%%rdx),%%zmm7,%%zmm7	\n\t"\
		"vsubpd	      %%zmm4,%%zmm0,%%zmm0		\n\t"\
		"vsubpd	      %%zmm5,%%zmm1,%%zmm1		\n\t"\
		"vmovaps	%%zmm0,     (%%rbx)	\n\t"\
		"vmovaps	%%zmm1,0x040(%%rbx)	\n\t"\
		"vaddpd	      %%zmm4,%%zmm4,%%zmm4		\n\t"\
		"vaddpd	      %%zmm5,%%zmm5,%%zmm5		\n\t"\
		"vaddpd	      %%zmm0,%%zmm4,%%zmm4		\n\t"\
		"vaddpd	      %%zmm1,%%zmm5,%%zmm5		\n\t"\
		"vmovaps	%%zmm4,     (%%rax)	\n\t"\
		"vmovaps	%%zmm5,0x040(%%rax)	\n\t"\
		"vsubpd	      %%zmm7,%%zmm2,%%zmm2		\n\t"\
		"vsubpd	      %%zmm6,%%zmm3,%%zmm3		\n\t"\
		"vmovaps	%%zmm2,     (%%rcx)	\n\t"\
		"vmovaps	%%zmm3,0x040(%%rdx)	\n\t"\
		"vaddpd	      %%zmm7,%%zmm7,%%zmm7		\n\t"\
		"vaddpd	      %%zmm6,%%zmm6,%%zmm6		\n\t"\
		"vaddpd	      %%zmm2,%%zmm7,%%zmm7		\n\t"\
		"vaddpd	      %%zmm3,%%zmm6,%%zmm6		\n\t"\
		"vmovaps	%%zmm7,     (%%rdx)	\n\t"\
		"vmovaps	%%zmm6,0x040(%%rcx)	\n\t"\
	/*** Now do four pass-2 4-DFTs, inputs from local store, outputs back to main array: ***/\
		"movq	%[__add0],%%rax		\n\t"\
		"movslq	%[__p1],%%rbx		\n\t"\
		"movslq	%[__p2],%%rcx		\n\t"\
		"movslq	%[__p3],%%rdx		\n\t"\
		"movslq	%[__p4],%%rdi		\n\t"\
		"shlq	$3,%%rbx			\n\t"\
		"shlq	$3,%%rcx			\n\t"\
		"shlq	$3,%%rdx			\n\t"\
		"shlq	$3,%%rdi			\n\t"\
		"addq	%%rax,%%rbx			\n\t"\
		"addq	%%rax,%%rcx			\n\t"\
		"addq	%%rax,%%rdx			\n\t"\
		"movq	%[__r1],%%rsi		\n\t"\
		"vmovaps	     (%%rsi),%%zmm0	\n\t"\
		"vmovaps	0x040(%%rsi),%%zmm1	\n\t"\
		"vmovaps	0x100(%%rsi),%%zmm2	\n\t"\
		"vmovaps	0x140(%%rsi),%%zmm3	\n\t"\
		"vsubpd	0x100(%%rsi),%%zmm0,%%zmm0	\n\t"\
		"vsubpd	0x140(%%rsi),%%zmm1,%%zmm1	\n\t"\
		"vaddpd	     (%%rsi),%%zmm2,%%zmm2	\n\t"\
		"vaddpd	0x040(%%rsi),%%zmm3,%%zmm3	\n\t"\
		"vmovaps	0x080(%%rsi),%%zmm4	\n\t"\
		"vmovaps	0x0c0(%%rsi),%%zmm5	\n\t"\
		"vmovaps	0x180(%%rsi),%%zmm6	\n\t"\
		"vmovaps	0x1c0(%%rsi),%%zmm7	\n\t"\
		"vsubpd	0x180(%%rsi),%%zmm4,%%zmm4	\n\t"\
		"vsubpd	0x1c0(%%rsi),%%zmm5,%%zmm5	\n\t"\
		"vaddpd	0x080(%%rsi),%%zmm6,%%zmm6	\n\t"\
		"vaddpd	0x0c0(%%rsi),%%zmm7,%%zmm7	\n\t"\
		"vsubpd	      %%zmm6,%%zmm2,%%zmm2		\n\t"\
		"vsubpd	      %%zmm7,%%zmm3,%%zmm3		\n\t"\
		"vaddpd	      %%zmm6,%%zmm6,%%zmm6		\n\t"\
		"vaddpd	      %%zmm7,%%zmm7,%%zmm7		\n\t"\
		"vmovaps	%%zmm2,    (%%rbx)	\n\t"\
		"vmovaps	%%zmm3,0x040(%%rbx)	\n\t"\
		"vaddpd	      %%zmm2,%%zmm6,%%zmm6	\n\t"\
		"vaddpd	      %%zmm3,%%zmm7,%%zmm7	\n\t"\
		"vmovaps	%%zmm6,     (%%rax)	\n\t"\
		"vmovaps	%%zmm7,0x040(%%rax)	\n\t"\
		"vsubpd	      %%zmm5,%%zmm0,%%zmm0	\n\t"\
		"vsubpd	      %%zmm4,%%zmm1,%%zmm1	\n\t"\
		"vaddpd	      %%zmm5,%%zmm5,%%zmm5	\n\t"\
		"vaddpd	      %%zmm4,%%zmm4,%%zmm4	\n\t"\
		"vmovaps	%%zmm0,     (%%rcx)	\n\t"\
		"vmovaps	%%zmm1,0x040(%%rdx)	\n\t"\
		"vaddpd	      %%zmm0,%%zmm5,%%zmm5	\n\t"\
		"vaddpd	      %%zmm1,%%zmm4,%%zmm4	\n\t"\
		"vmovaps	%%zmm5,     (%%rdx)	\n\t"\
		"vmovaps	%%zmm4,0x040(%%rcx)	\n\t"\
		"addq	%%rdi,%%rax			\n\t"\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	%%rdi,%%rcx			\n\t"\
		"addq	%%rdi,%%rdx			\n\t"\
		"movq	%[__r17],%%rsi		\n\t"\
		"vmovaps	     (%%rsi),%%zmm0	\n\t"\
		"vmovaps	0x040(%%rsi),%%zmm1	\n\t"\
		"vmovaps	0x100(%%rsi),%%zmm2	\n\t"\
		"vmovaps	0x140(%%rsi),%%zmm3	\n\t"\
		"vsubpd	0x140(%%rsi),%%zmm0,%%zmm0	\n\t"\
		"vsubpd	0x100(%%rsi),%%zmm1,%%zmm1	\n\t"\
		"vaddpd	0x040(%%rsi),%%zmm2,%%zmm2	\n\t"\
		"vaddpd	     (%%rsi),%%zmm3,%%zmm3	\n\t"\
		"vmovaps	0x080(%%rsi),%%zmm4	\n\t"\
		"vmovaps	0x0c0(%%rsi),%%zmm5	\n\t"\
		"vmovaps	0x180(%%rsi),%%zmm6	\n\t"\
		"vmovaps	0x1c0(%%rsi),%%zmm7	\n\t"\
		"vsubpd	0x0c0(%%rsi),%%zmm4,%%zmm4	\n\t"\
		"vaddpd	0x080(%%rsi),%%zmm5,%%zmm5	\n\t"\
		"vaddpd	0x1c0(%%rsi),%%zmm6,%%zmm6	\n\t"\
		"vsubpd	0x180(%%rsi),%%zmm7,%%zmm7	\n\t"\
		"movq	%[__isrt2],%%rsi	\n\t"\
		"vmulpd	     (%%rsi),%%zmm4,%%zmm4		\n\t"\
		"vmulpd	     (%%rsi),%%zmm5,%%zmm5		\n\t"\
		"vmulpd	     (%%rsi),%%zmm6,%%zmm6		\n\t"\
		"vmulpd	     (%%rsi),%%zmm7,%%zmm7		\n\t"\
		"vsubpd	      %%zmm6,%%zmm4,%%zmm4		\n\t"\
		"vsubpd	      %%zmm7,%%zmm5,%%zmm5		\n\t"\
		"vaddpd	      %%zmm6,%%zmm6,%%zmm6		\n\t"\
		"vaddpd	      %%zmm7,%%zmm7,%%zmm7		\n\t"\
		"vaddpd	      %%zmm4,%%zmm6,%%zmm6		\n\t"\
		"vaddpd	      %%zmm5,%%zmm7,%%zmm7		\n\t"\
		"vsubpd	      %%zmm4,%%zmm0,%%zmm0		\n\t"\
		"vsubpd	      %%zmm5,%%zmm2,%%zmm2		\n\t"\
		"vaddpd	      %%zmm4,%%zmm4,%%zmm4		\n\t"\
		"vaddpd	      %%zmm5,%%zmm5,%%zmm5		\n\t"\
		"vmovaps	%%zmm0,    (%%rbx)	\n\t"\
		"vmovaps	%%zmm2,0x040(%%rbx)	\n\t"\
		"vaddpd	      %%zmm0,%%zmm4,%%zmm4	\n\t"\
		"vaddpd	      %%zmm2,%%zmm5,%%zmm5	\n\t"\
		"vmovaps	%%zmm4,     (%%rax)	\n\t"\
		"vmovaps	%%zmm5,0x040(%%rax)	\n\t"\
		"vsubpd	      %%zmm7,%%zmm3,%%zmm3	\n\t"\
		"vsubpd	      %%zmm6,%%zmm1,%%zmm1	\n\t"\
		"vaddpd	      %%zmm7,%%zmm7,%%zmm7	\n\t"\
		"vaddpd	      %%zmm6,%%zmm6,%%zmm6	\n\t"\
		"vmovaps	%%zmm3,     (%%rcx)	\n\t"\
		"vmovaps	%%zmm1,0x040(%%rdx)	\n\t"\
		"vaddpd	      %%zmm3,%%zmm7,%%zmm7	\n\t"\
		"vaddpd	      %%zmm1,%%zmm6,%%zmm6	\n\t"\
		"vmovaps	%%zmm7,     (%%rdx)	\n\t"\
		"vmovaps	%%zmm6,0x040(%%rcx)	\n\t"\
		"addq	%%rdi,%%rax			\n\t"\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	%%rdi,%%rcx			\n\t"\
		"addq	%%rdi,%%rdx			\n\t"\
		"movq	%[__r9] ,%%rsi		\n\t"\
		"movq	%[__cc0],%%rdi		\n\t"\
		"vmovaps	0x080(%%rsi),%%zmm4	\n\t"\
		"vmovaps	0x180(%%rsi),%%zmm6	\n\t"\
		"vmovaps	0x0c0(%%rsi),%%zmm5	\n\t"\
		"vmovaps	0x1c0(%%rsi),%%zmm7	\n\t"\
		"vmovaps	%%zmm4,%%zmm0		\n\t"\
		"vmovaps	%%zmm6,%%zmm2		\n\t"\
		"vmovaps	%%zmm5,%%zmm1		\n\t"\
		"vmovaps	%%zmm7,%%zmm3		\n\t"\
		"vmulpd	     (%%rdi),%%zmm4,%%zmm4	\n\t"\
		"vmulpd	0x040(%%rdi),%%zmm6,%%zmm6	\n\t"\
		"vmulpd	0x040(%%rdi),%%zmm1,%%zmm1	\n\t"\
		"vmulpd	     (%%rdi),%%zmm3,%%zmm3	\n\t"\
		"vmulpd	     (%%rdi),%%zmm5,%%zmm5	\n\t"\
		"vmulpd	0x040(%%rdi),%%zmm7,%%zmm7	\n\t"\
		"vmulpd	0x040(%%rdi),%%zmm0,%%zmm0	\n\t"\
		"vmulpd	     (%%rdi),%%zmm2,%%zmm2	\n\t"\
		"vsubpd	      %%zmm1,%%zmm4,%%zmm4		\n\t"\
		"vsubpd	      %%zmm3,%%zmm6,%%zmm6		\n\t"\
		"vaddpd	      %%zmm0,%%zmm5,%%zmm5		\n\t"\
		"vaddpd	      %%zmm2,%%zmm7,%%zmm7		\n\t"\
		"vsubpd	      %%zmm6,%%zmm4,%%zmm4		\n\t"\
		"vsubpd	      %%zmm7,%%zmm5,%%zmm5		\n\t"\
		"vaddpd	      %%zmm6,%%zmm6,%%zmm6		\n\t"\
		"vaddpd	      %%zmm7,%%zmm7,%%zmm7		\n\t"\
		"vaddpd	      %%zmm4,%%zmm6,%%zmm6		\n\t"\
		"vaddpd	      %%zmm5,%%zmm7,%%zmm7		\n\t"\
		"movq	%[__isrt2],%%rdi	\n\t"\
		"vmovaps	0x100(%%rsi),%%zmm2	\n\t"\
		"vmovaps	0x140(%%rsi),%%zmm3	\n\t"\
		"vsubpd	0x140(%%rsi),%%zmm2,%%zmm2	\n\t"\
		"vaddpd	0x100(%%rsi),%%zmm3,%%zmm3	\n\t"\
		"vmulpd	     (%%rdi),%%zmm2,%%zmm2		\n\t"\
		"vmulpd	     (%%rdi),%%zmm3,%%zmm3		\n\t"\
		"vmovaps	     (%%rsi),%%zmm0	\n\t"\
		"vmovaps	0x040(%%rsi),%%zmm1	\n\t"\
		"movslq	%[__p4],%%rdi		\n\t"\
		"shlq	$3,%%rdi			\n\t"\
		"vsubpd	      %%zmm2,%%zmm0,%%zmm0		\n\t"\
		"vsubpd	      %%zmm3,%%zmm1,%%zmm1		\n\t"\
		"vaddpd	     (%%rsi),%%zmm2,%%zmm2	\n\t"\
		"vaddpd	0x040(%%rsi),%%zmm3,%%zmm3	\n\t"\
		"vsubpd	      %%zmm6,%%zmm2,%%zmm2		\n\t"\
		"vsubpd	      %%zmm7,%%zmm3,%%zmm3		\n\t"\
		"vaddpd	      %%zmm6,%%zmm6,%%zmm6		\n\t"\
		"vaddpd	      %%zmm7,%%zmm7,%%zmm7		\n\t"\
		"vmovaps	%%zmm2,    (%%rbx)	\n\t"\
		"vmovaps	%%zmm3,0x040(%%rbx)	\n\t"\
		"vaddpd	      %%zmm2,%%zmm6,%%zmm6	\n\t"\
		"vaddpd	      %%zmm3,%%zmm7,%%zmm7	\n\t"\
		"vmovaps	%%zmm6,     (%%rax)	\n\t"\
		"vmovaps	%%zmm7,0x040(%%rax)	\n\t"\
		"vsubpd	      %%zmm5,%%zmm0,%%zmm0	\n\t"\
		"vsubpd	      %%zmm4,%%zmm1,%%zmm1	\n\t"\
		"vaddpd	      %%zmm5,%%zmm5,%%zmm5	\n\t"\
		"vaddpd	      %%zmm4,%%zmm4,%%zmm4	\n\t"\
		"vmovaps	%%zmm0,     (%%rcx)	\n\t"\
		"vmovaps	%%zmm1,0x040(%%rdx)	\n\t"\
		"vaddpd	      %%zmm0,%%zmm5,%%zmm5	\n\t"\
		"vaddpd	      %%zmm1,%%zmm4,%%zmm4	\n\t"\
		"vmovaps	%%zmm5,     (%%rdx)	\n\t"\
		"vmovaps	%%zmm4,0x040(%%rcx)	\n\t"\
		"addq	%%rdi,%%rax			\n\t"\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	%%rdi,%%rcx			\n\t"\
		"addq	%%rdi,%%rdx			\n\t"\
		"movq	%[__r25],%%rsi		\n\t"\
		"movq	%[__cc0],%%rdi		\n\t"\
		"vmovaps	0x080(%%rsi),%%zmm4	\n\t"\
		"vmovaps	0x180(%%rsi),%%zmm6	\n\t"\
		"vmovaps	0x0c0(%%rsi),%%zmm5	\n\t"\
		"vmovaps	0x1c0(%%rsi),%%zmm7	\n\t"\
		"vmovaps	%%zmm4,%%zmm0		\n\t"\
		"vmovaps	%%zmm6,%%zmm2		\n\t"\
		"vmovaps	%%zmm5,%%zmm1		\n\t"\
		"vmovaps	%%zmm7,%%zmm3		\n\t"\
		"vmulpd	0x040(%%rdi),%%zmm4,%%zmm4	\n\t"\
		"vmulpd	     (%%rdi),%%zmm6,%%zmm6	\n\t"\
		"vmulpd	     (%%rdi),%%zmm1,%%zmm1	\n\t"\
		"vmulpd	0x040(%%rdi),%%zmm3,%%zmm3	\n\t"\
		"vmulpd	0x040(%%rdi),%%zmm5,%%zmm5	\n\t"\
		"vmulpd	     (%%rdi),%%zmm7,%%zmm7	\n\t"\
		"vmulpd	     (%%rdi),%%zmm0,%%zmm0	\n\t"\
		"vmulpd	0x040(%%rdi),%%zmm2,%%zmm2	\n\t"\
		"vsubpd	      %%zmm1,%%zmm4,%%zmm4		\n\t"\
		"vsubpd	      %%zmm3,%%zmm6,%%zmm6		\n\t"\
		"vaddpd	      %%zmm0,%%zmm5,%%zmm5		\n\t"\
		"vaddpd	      %%zmm2,%%zmm7,%%zmm7		\n\t"\
		"vsubpd	      %%zmm6,%%zmm4,%%zmm4		\n\t"\
		"vsubpd	      %%zmm7,%%zmm5,%%zmm5		\n\t"\
		"vaddpd	      %%zmm6,%%zmm6,%%zmm6		\n\t"\
		"vaddpd	      %%zmm7,%%zmm7,%%zmm7		\n\t"\
		"vaddpd	      %%zmm4,%%zmm6,%%zmm6		\n\t"\
		"vaddpd	      %%zmm5,%%zmm7,%%zmm7		\n\t"\
		"movq	%[__isrt2],%%rdi	\n\t"\
		"vmovaps	0x100(%%rsi),%%zmm2	\n\t"\
		"vmovaps	0x140(%%rsi),%%zmm3	\n\t"\
		"vaddpd	0x140(%%rsi),%%zmm2,%%zmm2	\n\t"\
		"vsubpd	0x100(%%rsi),%%zmm3,%%zmm3	\n\t"\
		"vmulpd	     (%%rdi),%%zmm2,%%zmm2		\n\t"\
		"vmulpd	     (%%rdi),%%zmm3,%%zmm3		\n\t"\
		"vmovaps	     (%%rsi),%%zmm0	\n\t"\
		"vmovaps	0x040(%%rsi),%%zmm1	\n\t"\
		"vsubpd	      %%zmm2,%%zmm0,%%zmm0		\n\t"\
		"vsubpd	      %%zmm3,%%zmm1,%%zmm1		\n\t"\
		"vaddpd	     (%%rsi),%%zmm2,%%zmm2	\n\t"\
		"vaddpd	0x040(%%rsi),%%zmm3,%%zmm3	\n\t"\
		"vsubpd	      %%zmm4,%%zmm0,%%zmm0		\n\t"\
		"vsubpd	      %%zmm5,%%zmm1,%%zmm1		\n\t"\
		"vaddpd	      %%zmm4,%%zmm4,%%zmm4		\n\t"\
		"vaddpd	      %%zmm5,%%zmm5,%%zmm5		\n\t"\
		"vmovaps	%%zmm0,     (%%rbx)	\n\t"\
		"vmovaps	%%zmm1,0x040(%%rbx)	\n\t"\
		"vaddpd	      %%zmm0,%%zmm4,%%zmm4	\n\t"\
		"vaddpd	      %%zmm1,%%zmm5,%%zmm5	\n\t"\
		"vmovaps	%%zmm4,     (%%rax)	\n\t"\
		"vmovaps	%%zmm5,0x040(%%rax)	\n\t"\
		"vsubpd	      %%zmm7,%%zmm2,%%zmm2	\n\t"\
		"vsubpd	      %%zmm6,%%zmm3,%%zmm3	\n\t"\
		"vaddpd	      %%zmm7,%%zmm7,%%zmm7	\n\t"\
		"vaddpd	      %%zmm6,%%zmm6,%%zmm6	\n\t"\
		"vmovaps	%%zmm2,     (%%rcx)	\n\t"\
		"vmovaps	%%zmm3,0x040(%%rdx)	\n\t"\
		"vaddpd	      %%zmm2,%%zmm7,%%zmm7	\n\t"\
		"vaddpd	      %%zmm3,%%zmm6,%%zmm6	\n\t"\
		"vmovaps	%%zmm7,     (%%rdx)	\n\t"\
		"vmovaps	%%zmm6,0x040(%%rcx)	\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__p1] "m" (Xp1)\
		 ,[__p2] "m" (Xp2)\
		 ,[__p3] "m" (Xp3)\
		 ,[__p4] "m" (Xp4)\
		 ,[__r1] "m" (Xr1)\
		 ,[__r3] "m" (Xr3)\
		 ,[__r5] "m" (Xr5)\
		 ,[__r7] "m" (Xr7)\
		 ,[__r9] "m" (Xr9)\
		 ,[__r11] "m" (Xr11)\
		 ,[__r13] "m" (Xr13)\
		 ,[__r15] "m" (Xr15)\
		 ,[__r17] "m" (Xr17)\
		 ,[__r19] "m" (Xr19)\
		 ,[__r21] "m" (Xr21)\
		 ,[__r23] "m" (Xr23)\
		 ,[__r25] "m" (Xr25)\
		 ,[__r27] "m" (Xr27)\
		 ,[__r29] "m" (Xr29)\
		 ,[__r31] "m" (Xr31)\
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__cc0] "m" (Xcc0)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
	);\
	}

#elif defined(USE_AVX)	// AVX and AVX2 both use 256-bit registers

	#define SSE2_RADIX16_DIT_NOTWIDDLE(Xadd0,Xp1,Xp2,Xp3,Xp4,Xr1,Xr3,Xr5,Xr7,Xr9,Xr11,Xr17,Xr25,Xisrt2,Xcc0)\
	{\
	__asm__ volatile (\
		"movq	%[__add0],%%rax		\n\t"\
		"movslq	%[__p1],%%rbx		\n\t"\
		"movslq	%[__p2],%%rcx		\n\t"\
		"movslq	%[__p3],%%rdx		\n\t"\
		"movslq	%[__p4],%%rdi		\n\t"\
		"shlq	$3,%%rbx			\n\t"\
		"shlq	$3,%%rcx			\n\t"\
		"shlq	$3,%%rdx			\n\t"\
		"shlq	$3,%%rdi			\n\t"\
		"addq	%%rax,%%rbx			\n\t"\
		"addq	%%rax,%%rcx			\n\t"\
		"addq	%%rax,%%rdx			\n\t"\
		"/* SSE2_RADIX4_DIT_0TWIDDLE_B(r1 ): */\n\t"\
		"vmovaps	     (%%rax),%%ymm2	\n\t"\
		"vmovaps	     (%%rcx),%%ymm6	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm3	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm7	\n\t"\
		"vmovaps	     (%%rbx),%%ymm0	\n\t"\
		"vmovaps	     (%%rdx),%%ymm4	\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm1	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm5	\n\t"\
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
		"movq	%[__r1], %%rsi		\n\t"\
		"vmovaps	%%ymm0,0x080(%%rsi)	\n\t"\
		"vmovaps	%%ymm2,0x0c0(%%rsi)	\n\t"\
		"vmovaps	%%ymm1,0x0a0(%%rsi)	\n\t"\
		"vmovaps	%%ymm3,0x060(%%rsi)	\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	%%ymm0,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	%%ymm2,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	%%ymm1,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	%%ymm3,%%ymm6,%%ymm6		\n\t"\
		"vmovaps	%%ymm4,     (%%rsi)	\n\t"\
		"vmovaps	%%ymm7,0x040(%%rsi)	\n\t"\
		"vmovaps	%%ymm5,0x020(%%rsi)	\n\t"\
		"vmovaps	%%ymm6,0x0e0(%%rsi)	\n\t"\
		"addq	%%rdi,%%rax			\n\t"\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	%%rdi,%%rcx			\n\t"\
		"addq	%%rdi,%%rdx			\n\t"\
		"/* SSE2_RADIX4_DIT_0TWIDDLE_B(r9 ): */\n\t"\
		"vmovaps	     (%%rax),%%ymm2	\n\t"\
		"vmovaps	     (%%rcx),%%ymm6	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm3	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm7	\n\t"\
		"vmovaps	     (%%rbx),%%ymm0	\n\t"\
		"vmovaps	     (%%rdx),%%ymm4	\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm1	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm5	\n\t"\
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
		"movq	%[__r9], %%rsi		\n\t"\
		"vmovaps	%%ymm0,0x080(%%rsi)	\n\t"\
		"vmovaps	%%ymm2,0x0c0(%%rsi)	\n\t"\
		"vmovaps	%%ymm1,0x0a0(%%rsi)	\n\t"\
		"vmovaps	%%ymm3,0x060(%%rsi)	\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	%%ymm0,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	%%ymm2,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	%%ymm1,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	%%ymm3,%%ymm6,%%ymm6		\n\t"\
		"vmovaps	%%ymm4,     (%%rsi)	\n\t"\
		"vmovaps	%%ymm7,0x040(%%rsi)	\n\t"\
		"vmovaps	%%ymm5,0x020(%%rsi)	\n\t"\
		"vmovaps	%%ymm6,0x0e0(%%rsi)	\n\t"\
		"addq	%%rdi,%%rax			\n\t"\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	%%rdi,%%rcx			\n\t"\
		"addq	%%rdi,%%rdx			\n\t"\
		"/* SSE2_RADIX4_DIT_0TWIDDLE_B(r17): */\n\t"\
		"vmovaps	     (%%rax),%%ymm2	\n\t"\
		"vmovaps	     (%%rcx),%%ymm6	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm3	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm7	\n\t"\
		"vmovaps	     (%%rbx),%%ymm0	\n\t"\
		"vmovaps	     (%%rdx),%%ymm4	\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm1	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm5	\n\t"\
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
		"movq	%[__r17], %%rsi		\n\t"\
		"vmovaps	%%ymm0,0x080(%%rsi)	\n\t"\
		"vmovaps	%%ymm2,0x0c0(%%rsi)	\n\t"\
		"vmovaps	%%ymm1,0x0a0(%%rsi)	\n\t"\
		"vmovaps	%%ymm3,0x060(%%rsi)	\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	%%ymm0,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	%%ymm2,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	%%ymm1,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	%%ymm3,%%ymm6,%%ymm6		\n\t"\
		"vmovaps	%%ymm4,     (%%rsi)	\n\t"\
		"vmovaps	%%ymm7,0x040(%%rsi)	\n\t"\
		"vmovaps	%%ymm5,0x020(%%rsi)	\n\t"\
		"vmovaps	%%ymm6,0x0e0(%%rsi)	\n\t"\
		"addq	%%rdi,%%rax			\n\t"\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	%%rdi,%%rcx			\n\t"\
		"addq	%%rdi,%%rdx			\n\t"\
		"/* SSE2_RADIX4_DIT_0TWIDDLE_B(r25): */\n\t"\
		"vmovaps	     (%%rax),%%ymm2	\n\t"\
		"vmovaps	     (%%rcx),%%ymm6	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm3	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm7	\n\t"\
		"vmovaps	     (%%rbx),%%ymm0	\n\t"\
		"vmovaps	     (%%rdx),%%ymm4	\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm1	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm5	\n\t"\
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
		"movq	%[__r25], %%rsi		\n\t"\
		"vmovaps	%%ymm0,0x080(%%rsi)	\n\t"\
		"vmovaps	%%ymm2,0x0c0(%%rsi)	\n\t"\
		"vmovaps	%%ymm1,0x0a0(%%rsi)	\n\t"\
		"vmovaps	%%ymm3,0x060(%%rsi)	\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	%%ymm0,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	%%ymm2,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	%%ymm1,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	%%ymm3,%%ymm6,%%ymm6		\n\t"\
		"vmovaps	%%ymm4,     (%%rsi)	\n\t"\
		"vmovaps	%%ymm7,0x040(%%rsi)	\n\t"\
		"vmovaps	%%ymm5,0x020(%%rsi)	\n\t"\
		"vmovaps	%%ymm6,0x0e0(%%rsi)	\n\t"\
		"movq	%[__r1],%%rax		\n\t"\
		"vmovaps	     (%%rax),%%ymm0	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1	\n\t"\
		"vmovaps	0x100(%%rax),%%ymm2	\n\t"\
		"vmovaps	0x120(%%rax),%%ymm3	\n\t"\
		"vsubpd	0x100(%%rax),%%ymm0,%%ymm0	\n\t"\
		"vsubpd	0x120(%%rax),%%ymm1,%%ymm1	\n\t"\
		"vaddpd	     (%%rax),%%ymm2,%%ymm2	\n\t"\
		"vaddpd	0x020(%%rax),%%ymm3,%%ymm3	\n\t"\
		"vmovaps	0x200(%%rax),%%ymm4	\n\t"\
		"vmovaps	0x220(%%rax),%%ymm5	\n\t"\
		"vmovaps	0x300(%%rax),%%ymm6	\n\t"\
		"vmovaps	0x320(%%rax),%%ymm7	\n\t"\
		"vsubpd	0x300(%%rax),%%ymm4,%%ymm4	\n\t"\
		"vsubpd	0x320(%%rax),%%ymm5,%%ymm5	\n\t"\
		"vaddpd	0x200(%%rax),%%ymm6,%%ymm6	\n\t"\
		"vaddpd	0x220(%%rax),%%ymm7,%%ymm7	\n\t"\
		"vsubpd	      %%ymm6,%%ymm2,%%ymm2		\n\t"\
		"vsubpd	      %%ymm7,%%ymm3,%%ymm3		\n\t"\
		"vmovaps	%%ymm2,0x200(%%rax)	\n\t"\
		"vmovaps	%%ymm3,0x220(%%rax)	\n\t"\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm2,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm3,%%ymm7,%%ymm7		\n\t"\
		"vmovaps	%%ymm6,     (%%rax)	\n\t"\
		"vmovaps	%%ymm7,0x020(%%rax)	\n\t"\
		"vsubpd	      %%ymm5,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	      %%ymm4,%%ymm1,%%ymm1		\n\t"\
		"vmovaps	%%ymm0,0x300(%%rax)	\n\t"\
		"vmovaps	%%ymm1,0x120(%%rax)	\n\t"\
		"vaddpd	      %%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	      %%ymm0,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm1,%%ymm4,%%ymm4		\n\t"\
		"vmovaps	%%ymm5,0x100(%%rax)	\n\t"\
		"vmovaps	%%ymm4,0x320(%%rax)	\n\t"\
		"movq	%[__r5],%%rax		\n\t"\
		"movq	%[__isrt2],%%rbx	\n\t"\
		"vmovaps	(%%rbx),%%ymm2		\n\t"\
		"vmovaps	0x200(%%rax),%%ymm4	\n\t"\
		"vmovaps	0x220(%%rax),%%ymm5	\n\t"\
		"vmovaps	0x300(%%rax),%%ymm0	\n\t"\
		"vmovaps	0x320(%%rax),%%ymm1	\n\t"\
		"vaddpd	0x220(%%rax),%%ymm4,%%ymm4	\n\t"\
		"vsubpd	0x200(%%rax),%%ymm5,%%ymm5	\n\t"\
		"vsubpd	0x320(%%rax),%%ymm0,%%ymm0	\n\t"\
		"vaddpd	0x300(%%rax),%%ymm1,%%ymm1	\n\t"\
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
		"vmovaps	0x100(%%rax),%%ymm2	\n\t"\
		"vmovaps	0x120(%%rax),%%ymm3	\n\t"\
		"vsubpd	0x120(%%rax),%%ymm0,%%ymm0	\n\t"\
		"vsubpd	0x100(%%rax),%%ymm1,%%ymm1	\n\t"\
		"vaddpd	     (%%rax),%%ymm3,%%ymm3	\n\t"\
		"vaddpd	0x020(%%rax),%%ymm2,%%ymm2	\n\t"\
		"vsubpd	      %%ymm4,%%ymm3,%%ymm3		\n\t"\
		"vsubpd	      %%ymm5,%%ymm1,%%ymm1		\n\t"\
		"vmovaps	%%ymm3,0x200(%%rax)	\n\t"\
		"vmovaps	%%ymm1,0x220(%%rax)	\n\t"\
		"vaddpd	      %%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	      %%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm3,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	      %%ymm1,%%ymm5,%%ymm5		\n\t"\
		"vmovaps	%%ymm4,     (%%rax)	\n\t"\
		"vmovaps	%%ymm5,0x020(%%rax)	\n\t"\
		"vsubpd	      %%ymm7,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	      %%ymm6,%%ymm2,%%ymm2		\n\t"\
		"vmovaps	%%ymm0,0x300(%%rax)	\n\t"\
		"vmovaps	%%ymm2,0x120(%%rax)	\n\t"\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm0,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm2,%%ymm6,%%ymm6		\n\t"\
		"vmovaps	%%ymm7,0x100(%%rax)	\n\t"\
		"vmovaps	%%ymm6,0x320(%%rax)	\n\t"\
		"movq	%[__r3],%%rax		\n\t"\
		"movq	%[__isrt2],%%rbx	\n\t"\
		"movq	%[__cc0],%%rcx		\n\t"\
		"vmovaps	0x300(%%rax),%%ymm0	\n\t"\
		"vmovaps	0x320(%%rax),%%ymm1	\n\t"\
		"vmovaps	0x300(%%rax),%%ymm2	\n\t"\
		"vmovaps	0x320(%%rax),%%ymm3	\n\t"\
		"vmulpd	0x020(%%rcx),%%ymm0,%%ymm0	\n\t"\
		"vmulpd	0x020(%%rcx),%%ymm1,%%ymm1	\n\t"\
		"vmulpd	     (%%rcx),%%ymm2,%%ymm2	\n\t"\
		"vmulpd	     (%%rcx),%%ymm3,%%ymm3	\n\t"\
		"vsubpd	      %%ymm2,%%ymm1,%%ymm1		\n\t"\
		"vaddpd	      %%ymm3,%%ymm0,%%ymm0		\n\t"\
		"vmovaps	0x200(%%rax),%%ymm4	\n\t"\
		"vmovaps	0x220(%%rax),%%ymm5	\n\t"\
		"vmovaps	0x200(%%rax),%%ymm6	\n\t"\
		"vmovaps	0x220(%%rax),%%ymm7	\n\t"\
		"vmulpd	     (%%rcx),%%ymm4,%%ymm4	\n\t"\
		"vmulpd	     (%%rcx),%%ymm5,%%ymm5	\n\t"\
		"vmulpd	0x020(%%rcx),%%ymm6,%%ymm6	\n\t"\
		"vmulpd	0x020(%%rcx),%%ymm7,%%ymm7	\n\t"\
		"vsubpd	      %%ymm6,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm7,%%ymm4,%%ymm4		\n\t"\
		"vmovaps	%%ymm5,%%ymm7		\n\t"\
		"vmovaps	%%ymm4,%%ymm6		\n\t"\
		"vaddpd	      %%ymm0,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	      %%ymm1,%%ymm5,%%ymm5		\n\t"\
		"vsubpd	      %%ymm0,%%ymm6,%%ymm6		\n\t"\
		"vsubpd	      %%ymm1,%%ymm7,%%ymm7		\n\t"\
		"vmovaps	0x100(%%rax),%%ymm2	\n\t"\
		"vmovaps	0x120(%%rax),%%ymm3	\n\t"\
		"vmovaps	     (%%rax),%%ymm0	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1	\n\t"\
		"vaddpd	0x120(%%rax),%%ymm2,%%ymm2	\n\t"\
		"vsubpd	0x100(%%rax),%%ymm3,%%ymm3	\n\t"\
		"vmulpd	     (%%rbx),%%ymm2,%%ymm2	\n\t"\
		"vmulpd	     (%%rbx),%%ymm3,%%ymm3	\n\t"\
		"vsubpd	      %%ymm2,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	      %%ymm3,%%ymm1,%%ymm1		\n\t"\
		"vaddpd	      %%ymm2,%%ymm2,%%ymm2		\n\t"\
		"vaddpd	      %%ymm3,%%ymm3,%%ymm3		\n\t"\
		"vaddpd	      %%ymm0,%%ymm2,%%ymm2		\n\t"\
		"vaddpd	      %%ymm1,%%ymm3,%%ymm3		\n\t"\
		"vsubpd	      %%ymm4,%%ymm2,%%ymm2		\n\t"\
		"vsubpd	      %%ymm5,%%ymm3,%%ymm3		\n\t"\
		"vmovaps	%%ymm2,0x200(%%rax)	\n\t"\
		"vmovaps	%%ymm3,0x220(%%rax)	\n\t"\
		"vaddpd	      %%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	      %%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm2,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	      %%ymm3,%%ymm5,%%ymm5		\n\t"\
		"vmovaps	%%ymm4,     (%%rax)	\n\t"\
		"vmovaps	%%ymm5,0x020(%%rax)	\n\t"\
		"vsubpd	      %%ymm7,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	      %%ymm6,%%ymm1,%%ymm1		\n\t"\
		"vmovaps	%%ymm0,0x300(%%rax)	\n\t"\
		"vmovaps	%%ymm1,0x120(%%rax)	\n\t"\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm0,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm1,%%ymm6,%%ymm6		\n\t"\
		"vmovaps	%%ymm7,0x100(%%rax)	\n\t"\
		"vmovaps	%%ymm6,0x320(%%rax)	\n\t"\
		"movq	%[__r7],%%rax		\n\t"\
		"movq	%[__isrt2],%%rbx	\n\t"\
		"movq	%[__cc0],%%rcx		\n\t"\
		"vmovaps	0x300(%%rax),%%ymm0	\n\t"\
		"vmovaps	0x320(%%rax),%%ymm1	\n\t"\
		"vmovaps	0x300(%%rax),%%ymm2	\n\t"\
		"vmovaps	0x320(%%rax),%%ymm3	\n\t"\
		"vmulpd	     (%%rcx),%%ymm0,%%ymm0	\n\t"\
		"vmulpd	     (%%rcx),%%ymm1,%%ymm1	\n\t"\
		"vmulpd	0x020(%%rcx),%%ymm2,%%ymm2	\n\t"\
		"vmulpd	0x020(%%rcx),%%ymm3,%%ymm3	\n\t"\
		"vsubpd	      %%ymm2,%%ymm1,%%ymm1		\n\t"\
		"vaddpd	      %%ymm3,%%ymm0,%%ymm0		\n\t"\
		"vmovaps	0x200(%%rax),%%ymm4	\n\t"\
		"vmovaps	0x220(%%rax),%%ymm5	\n\t"\
		"vmovaps	0x200(%%rax),%%ymm6	\n\t"\
		"vmovaps	0x220(%%rax),%%ymm7	\n\t"\
		"vmulpd	0x020(%%rcx),%%ymm4,%%ymm4	\n\t"\
		"vmulpd	0x020(%%rcx),%%ymm5,%%ymm5	\n\t"\
		"vmulpd	     (%%rcx),%%ymm6,%%ymm6	\n\t"\
		"vmulpd	     (%%rcx),%%ymm7,%%ymm7	\n\t"\
		"vsubpd	      %%ymm6,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm7,%%ymm4,%%ymm4		\n\t"\
		"vmovaps	%%ymm5,%%ymm7		\n\t"\
		"vmovaps	%%ymm4,%%ymm6		\n\t"\
		"vaddpd	      %%ymm0,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	      %%ymm1,%%ymm5,%%ymm5		\n\t"\
		"vsubpd	      %%ymm0,%%ymm6,%%ymm6		\n\t"\
		"vsubpd	      %%ymm1,%%ymm7,%%ymm7		\n\t"\
		"vmovaps	0x100(%%rax),%%ymm2	\n\t"\
		"vmovaps	0x120(%%rax),%%ymm3	\n\t"\
		"vmovaps	     (%%rax),%%ymm0	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1	\n\t"\
		"vsubpd	0x120(%%rax),%%ymm2,%%ymm2	\n\t"\
		"vaddpd	0x100(%%rax),%%ymm3,%%ymm3	\n\t"\
		"vmulpd	     (%%rbx),%%ymm2,%%ymm2	\n\t"\
		"vmulpd	     (%%rbx),%%ymm3,%%ymm3	\n\t"\
		"vsubpd	      %%ymm2,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	      %%ymm3,%%ymm1,%%ymm1		\n\t"\
		"vaddpd	      %%ymm2,%%ymm2,%%ymm2		\n\t"\
		"vaddpd	      %%ymm3,%%ymm3,%%ymm3		\n\t"\
		"vaddpd	      %%ymm0,%%ymm2,%%ymm2		\n\t"\
		"vaddpd	      %%ymm1,%%ymm3,%%ymm3		\n\t"\
		"vsubpd	      %%ymm6,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	      %%ymm7,%%ymm1,%%ymm1		\n\t"\
		"vmovaps	%%ymm0,0x200(%%rax)	\n\t"\
		"vmovaps	%%ymm1,0x220(%%rax)	\n\t"\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm0,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm1,%%ymm7,%%ymm7		\n\t"\
		"vmovaps	%%ymm6,     (%%rax)	\n\t"\
		"vmovaps	%%ymm7,0x020(%%rax)	\n\t"\
		"vsubpd	      %%ymm5,%%ymm2,%%ymm2		\n\t"\
		"vsubpd	      %%ymm4,%%ymm3,%%ymm3		\n\t"\
		"vmovaps	%%ymm2,0x300(%%rax)	\n\t"\
		"vmovaps	%%ymm3,0x120(%%rax)	\n\t"\
		"vaddpd	      %%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	      %%ymm2,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm3,%%ymm4,%%ymm4		\n\t"\
		"vmovaps	%%ymm5,0x100(%%rax)	\n\t"\
		"vmovaps	%%ymm4,0x320(%%rax)	\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__p1] "m" (Xp1)\
		 ,[__p2] "m" (Xp2)\
		 ,[__p3] "m" (Xp3)\
		 ,[__p4] "m" (Xp4)\
		 ,[__r1] "m" (Xr1)\
		 ,[__r3] "m" (Xr3)\
		 ,[__r5] "m" (Xr5)\
		 ,[__r7] "m" (Xr7)\
		 ,[__r9] "m" (Xr9)\
		 ,[__r17] "m" (Xr17)\
		 ,[__r25] "m" (Xr25)\
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__cc0] "m" (Xcc0)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
	);\
	}


	#define SSE2_RADIX16_DIF_NOTWIDDLE(Xadd0,Xp1,Xp2,Xp3,Xp4,Xr1,Xr3,Xr5,Xr7,Xr9,Xr11,Xr13,Xr15,Xr17,Xr19,Xr21,Xr23,Xr25,Xr27,Xr29,Xr31,Xisrt2,Xcc0)\
	{\
	__asm__ volatile (\
		"/* SSE2_RADIX4_DIF_IN_PLACE(r1 , r17, r9 , r25): */\n\t"\
		"movq	%[__r1] ,%%rax		\n\t"\
		"movq	%[__r17],%%rbx		\n\t"\
		"movq	%[__r9] ,%%rcx		\n\t"\
		"movq	%[__r25],%%rdx		\n\t"\
		"vmovaps	     (%%rax),%%ymm0	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1	\n\t"\
		"vmovaps	     (%%rax),%%ymm2	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm3	\n\t"\
		"vaddpd	     (%%rbx),%%ymm0,%%ymm0	\n\t"\
		"vaddpd	0x020(%%rbx),%%ymm1,%%ymm1	\n\t"\
		"vsubpd	     (%%rbx),%%ymm2,%%ymm2	\n\t"\
		"vsubpd	0x020(%%rbx),%%ymm3,%%ymm3	\n\t"\
		"vmovaps	     (%%rcx),%%ymm4	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5	\n\t"\
		"vmovaps	     (%%rcx),%%ymm6	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm7	\n\t"\
		"vaddpd	     (%%rdx),%%ymm4,%%ymm4	\n\t"\
		"vaddpd	0x020(%%rdx),%%ymm5,%%ymm5	\n\t"\
		"vsubpd	     (%%rdx),%%ymm6,%%ymm6	\n\t"\
		"vsubpd	0x020(%%rdx),%%ymm7,%%ymm7	\n\t"\
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
		"/* SSE2_RADIX4_DIF_IN_PLACE(r5 , r21, r13, r29): */\n\t"\
		"movq	%[__r5] ,%%rax		\n\t"\
		"movq	%[__r21],%%rbx		\n\t"\
		"movq	%[__r13],%%rcx		\n\t"\
		"movq	%[__r29],%%rdx		\n\t"\
		"vmovaps	     (%%rax),%%ymm0	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1	\n\t"\
		"vmovaps	     (%%rax),%%ymm2	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm3	\n\t"\
		"vaddpd	     (%%rbx),%%ymm0,%%ymm0	\n\t"\
		"vaddpd	0x020(%%rbx),%%ymm1,%%ymm1	\n\t"\
		"vsubpd	     (%%rbx),%%ymm2,%%ymm2	\n\t"\
		"vsubpd	0x020(%%rbx),%%ymm3,%%ymm3	\n\t"\
		"vmovaps	     (%%rcx),%%ymm4	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5	\n\t"\
		"vmovaps	     (%%rcx),%%ymm6	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm7	\n\t"\
		"vaddpd	     (%%rdx),%%ymm4,%%ymm4	\n\t"\
		"vaddpd	0x020(%%rdx),%%ymm5,%%ymm5	\n\t"\
		"vsubpd	     (%%rdx),%%ymm6,%%ymm6	\n\t"\
		"vsubpd	0x020(%%rdx),%%ymm7,%%ymm7	\n\t"\
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
		"/* SSE2_RADIX4_DIF_IN_PLACE(r3 , r19, r11, r27): */\n\t"\
		"movq	%[__r3] ,%%rax		\n\t"\
		"movq	%[__r19],%%rbx		\n\t"\
		"movq	%[__r11],%%rcx		\n\t"\
		"movq	%[__r27],%%rdx		\n\t"\
		"vmovaps	     (%%rax),%%ymm0	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1	\n\t"\
		"vmovaps	     (%%rax),%%ymm2	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm3	\n\t"\
		"vaddpd	     (%%rbx),%%ymm0,%%ymm0	\n\t"\
		"vaddpd	0x020(%%rbx),%%ymm1,%%ymm1	\n\t"\
		"vsubpd	     (%%rbx),%%ymm2,%%ymm2	\n\t"\
		"vsubpd	0x020(%%rbx),%%ymm3,%%ymm3	\n\t"\
		"vmovaps	     (%%rcx),%%ymm4	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5	\n\t"\
		"vmovaps	     (%%rcx),%%ymm6	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm7	\n\t"\
		"vaddpd	     (%%rdx),%%ymm4,%%ymm4	\n\t"\
		"vaddpd	0x020(%%rdx),%%ymm5,%%ymm5	\n\t"\
		"vsubpd	     (%%rdx),%%ymm6,%%ymm6	\n\t"\
		"vsubpd	0x020(%%rdx),%%ymm7,%%ymm7	\n\t"\
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
		"/* SSE2_RADIX4_DIF_IN_PLACE(r7 , r23, r15, r31): */\n\t"\
		"movq	%[__r7] ,%%rax		\n\t"\
		"movq	%[__r23],%%rbx		\n\t"\
		"movq	%[__r15],%%rcx		\n\t"\
		"movq	%[__r31],%%rdx		\n\t"\
		"vmovaps	     (%%rax),%%ymm0	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1	\n\t"\
		"vmovaps	     (%%rax),%%ymm2	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm3	\n\t"\
		"vaddpd	     (%%rbx),%%ymm0,%%ymm0	\n\t"\
		"vaddpd	0x020(%%rbx),%%ymm1,%%ymm1	\n\t"\
		"vsubpd	     (%%rbx),%%ymm2,%%ymm2	\n\t"\
		"vsubpd	0x020(%%rbx),%%ymm3,%%ymm3	\n\t"\
		"vmovaps	     (%%rcx),%%ymm4	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5	\n\t"\
		"vmovaps	     (%%rcx),%%ymm6	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm7	\n\t"\
		"vaddpd	     (%%rdx),%%ymm4,%%ymm4	\n\t"\
		"vaddpd	0x020(%%rdx),%%ymm5,%%ymm5	\n\t"\
		"vsubpd	     (%%rdx),%%ymm6,%%ymm6	\n\t"\
		"vsubpd	0x020(%%rdx),%%ymm7,%%ymm7	\n\t"\
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
		"movq	%[__add0],%%rax		\n\t"\
		"movslq	%[__p1],%%rbx		\n\t"\
		"movslq	%[__p2],%%rcx		\n\t"\
		"movslq	%[__p3],%%rdx		\n\t"\
		"movslq	%[__p4],%%rdi		\n\t"\
		"shlq	$3,%%rbx			\n\t"\
		"shlq	$3,%%rcx			\n\t"\
		"shlq	$3,%%rdx			\n\t"\
		"shlq	$3,%%rdi			\n\t"\
		"addq	%%rax,%%rbx			\n\t"\
		"addq	%%rax,%%rcx			\n\t"\
		"addq	%%rax,%%rdx			\n\t"\
		"movq	%[__r1],%%rsi		\n\t"\
		"vmovaps	     (%%rsi),%%ymm0	\n\t"\
		"vmovaps	0x020(%%rsi),%%ymm1	\n\t"\
		"vmovaps	0x080(%%rsi),%%ymm2	\n\t"\
		"vmovaps	0x0a0(%%rsi),%%ymm3	\n\t"\
		"vsubpd	0x080(%%rsi),%%ymm0,%%ymm0	\n\t"\
		"vsubpd	0x0a0(%%rsi),%%ymm1,%%ymm1	\n\t"\
		"vaddpd	     (%%rsi),%%ymm2,%%ymm2	\n\t"\
		"vaddpd	0x020(%%rsi),%%ymm3,%%ymm3	\n\t"\
		"vmovaps	0x040(%%rsi),%%ymm4	\n\t"\
		"vmovaps	0x060(%%rsi),%%ymm5	\n\t"\
		"vmovaps	0x0c0(%%rsi),%%ymm6	\n\t"\
		"vmovaps	0x0e0(%%rsi),%%ymm7	\n\t"\
		"vsubpd	0x0c0(%%rsi),%%ymm4,%%ymm4	\n\t"\
		"vsubpd	0x0e0(%%rsi),%%ymm5,%%ymm5	\n\t"\
		"vaddpd	0x040(%%rsi),%%ymm6,%%ymm6	\n\t"\
		"vaddpd	0x060(%%rsi),%%ymm7,%%ymm7	\n\t"\
		"vsubpd	      %%ymm6,%%ymm2,%%ymm2		\n\t"\
		"vsubpd	      %%ymm7,%%ymm3,%%ymm3		\n\t"\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vmovaps	%%ymm2,    (%%rbx)	\n\t"\
		"vmovaps	%%ymm3,0x020(%%rbx)	\n\t"\
		"vaddpd	      %%ymm2,%%ymm6,%%ymm6	\n\t"\
		"vaddpd	      %%ymm3,%%ymm7,%%ymm7	\n\t"\
		"vmovaps	%%ymm6,     (%%rax)	\n\t"\
		"vmovaps	%%ymm7,0x020(%%rax)	\n\t"\
		"vsubpd	      %%ymm5,%%ymm0,%%ymm0	\n\t"\
		"vsubpd	      %%ymm4,%%ymm1,%%ymm1	\n\t"\
		"vaddpd	      %%ymm5,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	      %%ymm4,%%ymm4,%%ymm4	\n\t"\
		"vmovaps	%%ymm0,     (%%rcx)	\n\t"\
		"vmovaps	%%ymm1,0x020(%%rdx)	\n\t"\
		"vaddpd	      %%ymm0,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	      %%ymm1,%%ymm4,%%ymm4	\n\t"\
		"vmovaps	%%ymm5,     (%%rdx)	\n\t"\
		"vmovaps	%%ymm4,0x020(%%rcx)	\n\t"\
		"addq	%%rdi,%%rax			\n\t"\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	%%rdi,%%rcx			\n\t"\
		"addq	%%rdi,%%rdx			\n\t"\
		"movq	%[__r17],%%rsi		\n\t"\
		"vmovaps	     (%%rsi),%%ymm0	\n\t"\
		"vmovaps	0x020(%%rsi),%%ymm1	\n\t"\
		"vmovaps	0x080(%%rsi),%%ymm2	\n\t"\
		"vmovaps	0x0a0(%%rsi),%%ymm3	\n\t"\
		"vsubpd	0x0a0(%%rsi),%%ymm0,%%ymm0	\n\t"\
		"vsubpd	0x080(%%rsi),%%ymm1,%%ymm1	\n\t"\
		"vaddpd	0x020(%%rsi),%%ymm2,%%ymm2	\n\t"\
		"vaddpd	     (%%rsi),%%ymm3,%%ymm3	\n\t"\
		"vmovaps	0x040(%%rsi),%%ymm4	\n\t"\
		"vmovaps	0x060(%%rsi),%%ymm5	\n\t"\
		"vmovaps	0x0c0(%%rsi),%%ymm6	\n\t"\
		"vmovaps	0x0e0(%%rsi),%%ymm7	\n\t"\
		"vsubpd	0x060(%%rsi),%%ymm4,%%ymm4	\n\t"\
		"vaddpd	0x040(%%rsi),%%ymm5,%%ymm5	\n\t"\
		"vaddpd	0x0e0(%%rsi),%%ymm6,%%ymm6	\n\t"\
		"vsubpd	0x0c0(%%rsi),%%ymm7,%%ymm7	\n\t"\
		"movq	%[__isrt2],%%rsi	\n\t"\
		"vmulpd	     (%%rsi),%%ymm4,%%ymm4		\n\t"\
		"vmulpd	     (%%rsi),%%ymm5,%%ymm5		\n\t"\
		"vmulpd	     (%%rsi),%%ymm6,%%ymm6		\n\t"\
		"vmulpd	     (%%rsi),%%ymm7,%%ymm7		\n\t"\
		"vsubpd	      %%ymm6,%%ymm4,%%ymm4		\n\t"\
		"vsubpd	      %%ymm7,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm4,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm5,%%ymm7,%%ymm7		\n\t"\
		"vsubpd	      %%ymm4,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	      %%ymm5,%%ymm2,%%ymm2		\n\t"\
		"vaddpd	      %%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	      %%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vmovaps	%%ymm0,    (%%rbx)	\n\t"\
		"vmovaps	%%ymm2,0x020(%%rbx)	\n\t"\
		"vaddpd	      %%ymm0,%%ymm4,%%ymm4	\n\t"\
		"vaddpd	      %%ymm2,%%ymm5,%%ymm5	\n\t"\
		"vmovaps	%%ymm4,     (%%rax)	\n\t"\
		"vmovaps	%%ymm5,0x020(%%rax)	\n\t"\
		"vsubpd	      %%ymm7,%%ymm3,%%ymm3	\n\t"\
		"vsubpd	      %%ymm6,%%ymm1,%%ymm1	\n\t"\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7	\n\t"\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6	\n\t"\
		"vmovaps	%%ymm3,     (%%rcx)	\n\t"\
		"vmovaps	%%ymm1,0x020(%%rdx)	\n\t"\
		"vaddpd	      %%ymm3,%%ymm7,%%ymm7	\n\t"\
		"vaddpd	      %%ymm1,%%ymm6,%%ymm6	\n\t"\
		"vmovaps	%%ymm7,     (%%rdx)	\n\t"\
		"vmovaps	%%ymm6,0x020(%%rcx)	\n\t"\
		"addq	%%rdi,%%rax			\n\t"\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	%%rdi,%%rcx			\n\t"\
		"addq	%%rdi,%%rdx			\n\t"\
		"movq	%[__r9] ,%%rsi		\n\t"\
		"movq	%[__cc0],%%rdi		\n\t"\
		"vmovaps	0x040(%%rsi),%%ymm4	\n\t"\
		"vmovaps	0x0c0(%%rsi),%%ymm6	\n\t"\
		"vmovaps	0x060(%%rsi),%%ymm5	\n\t"\
		"vmovaps	0x0e0(%%rsi),%%ymm7	\n\t"\
		"vmovaps	%%ymm4,%%ymm0		\n\t"\
		"vmovaps	%%ymm6,%%ymm2		\n\t"\
		"vmovaps	%%ymm5,%%ymm1		\n\t"\
		"vmovaps	%%ymm7,%%ymm3		\n\t"\
		"vmulpd	     (%%rdi),%%ymm4,%%ymm4	\n\t"\
		"vmulpd	0x020(%%rdi),%%ymm6,%%ymm6	\n\t"\
		"vmulpd	0x020(%%rdi),%%ymm1,%%ymm1	\n\t"\
		"vmulpd	     (%%rdi),%%ymm3,%%ymm3	\n\t"\
		"vmulpd	     (%%rdi),%%ymm5,%%ymm5	\n\t"\
		"vmulpd	0x020(%%rdi),%%ymm7,%%ymm7	\n\t"\
		"vmulpd	0x020(%%rdi),%%ymm0,%%ymm0	\n\t"\
		"vmulpd	     (%%rdi),%%ymm2,%%ymm2	\n\t"\
		"vsubpd	      %%ymm1,%%ymm4,%%ymm4		\n\t"\
		"vsubpd	      %%ymm3,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm0,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm2,%%ymm7,%%ymm7		\n\t"\
		"vsubpd	      %%ymm6,%%ymm4,%%ymm4		\n\t"\
		"vsubpd	      %%ymm7,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm4,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm5,%%ymm7,%%ymm7		\n\t"\
		"movq	%[__isrt2],%%rdi	\n\t"\
		"vmovaps	0x080(%%rsi),%%ymm2	\n\t"\
		"vmovaps	0x0a0(%%rsi),%%ymm3	\n\t"\
		"vsubpd	0x0a0(%%rsi),%%ymm2,%%ymm2	\n\t"\
		"vaddpd	0x080(%%rsi),%%ymm3,%%ymm3	\n\t"\
		"vmulpd	     (%%rdi),%%ymm2,%%ymm2		\n\t"\
		"vmulpd	     (%%rdi),%%ymm3,%%ymm3		\n\t"\
		"vmovaps	     (%%rsi),%%ymm0	\n\t"\
		"vmovaps	0x020(%%rsi),%%ymm1	\n\t"\
		"movslq	%[__p4],%%rdi		\n\t"\
		"shlq	$3,%%rdi			\n\t"\
		"vsubpd	      %%ymm2,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	      %%ymm3,%%ymm1,%%ymm1		\n\t"\
		"vaddpd	     (%%rsi),%%ymm2,%%ymm2	\n\t"\
		"vaddpd	0x020(%%rsi),%%ymm3,%%ymm3	\n\t"\
		"vsubpd	      %%ymm6,%%ymm2,%%ymm2		\n\t"\
		"vsubpd	      %%ymm7,%%ymm3,%%ymm3		\n\t"\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vmovaps	%%ymm2,    (%%rbx)	\n\t"\
		"vmovaps	%%ymm3,0x020(%%rbx)	\n\t"\
		"vaddpd	      %%ymm2,%%ymm6,%%ymm6	\n\t"\
		"vaddpd	      %%ymm3,%%ymm7,%%ymm7	\n\t"\
		"vmovaps	%%ymm6,     (%%rax)	\n\t"\
		"vmovaps	%%ymm7,0x020(%%rax)	\n\t"\
		"vsubpd	      %%ymm5,%%ymm0,%%ymm0	\n\t"\
		"vsubpd	      %%ymm4,%%ymm1,%%ymm1	\n\t"\
		"vaddpd	      %%ymm5,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	      %%ymm4,%%ymm4,%%ymm4	\n\t"\
		"vmovaps	%%ymm0,     (%%rcx)	\n\t"\
		"vmovaps	%%ymm1,0x020(%%rdx)	\n\t"\
		"vaddpd	      %%ymm0,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	      %%ymm1,%%ymm4,%%ymm4	\n\t"\
		"vmovaps	%%ymm5,     (%%rdx)	\n\t"\
		"vmovaps	%%ymm4,0x020(%%rcx)	\n\t"\
		"addq	%%rdi,%%rax			\n\t"\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	%%rdi,%%rcx			\n\t"\
		"addq	%%rdi,%%rdx			\n\t"\
		"movq	%[__r25],%%rsi		\n\t"\
		"movq	%[__cc0],%%rdi		\n\t"\
		"vmovaps	0x040(%%rsi),%%ymm4	\n\t"\
		"vmovaps	0x0c0(%%rsi),%%ymm6	\n\t"\
		"vmovaps	0x060(%%rsi),%%ymm5	\n\t"\
		"vmovaps	0x0e0(%%rsi),%%ymm7	\n\t"\
		"vmovaps	%%ymm4,%%ymm0		\n\t"\
		"vmovaps	%%ymm6,%%ymm2		\n\t"\
		"vmovaps	%%ymm5,%%ymm1		\n\t"\
		"vmovaps	%%ymm7,%%ymm3		\n\t"\
		"vmulpd	0x020(%%rdi),%%ymm4,%%ymm4	\n\t"\
		"vmulpd	     (%%rdi),%%ymm6,%%ymm6	\n\t"\
		"vmulpd	     (%%rdi),%%ymm1,%%ymm1	\n\t"\
		"vmulpd	0x020(%%rdi),%%ymm3,%%ymm3	\n\t"\
		"vmulpd	0x020(%%rdi),%%ymm5,%%ymm5	\n\t"\
		"vmulpd	     (%%rdi),%%ymm7,%%ymm7	\n\t"\
		"vmulpd	     (%%rdi),%%ymm0,%%ymm0	\n\t"\
		"vmulpd	0x020(%%rdi),%%ymm2,%%ymm2	\n\t"\
		"vsubpd	      %%ymm1,%%ymm4,%%ymm4		\n\t"\
		"vsubpd	      %%ymm3,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm0,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm2,%%ymm7,%%ymm7		\n\t"\
		"vsubpd	      %%ymm6,%%ymm4,%%ymm4		\n\t"\
		"vsubpd	      %%ymm7,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm4,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm5,%%ymm7,%%ymm7		\n\t"\
		"movq	%[__isrt2],%%rdi	\n\t"\
		"vmovaps	0x080(%%rsi),%%ymm2	\n\t"\
		"vmovaps	0x0a0(%%rsi),%%ymm3	\n\t"\
		"vaddpd	0x0a0(%%rsi),%%ymm2,%%ymm2	\n\t"\
		"vsubpd	0x080(%%rsi),%%ymm3,%%ymm3	\n\t"\
		"vmulpd	     (%%rdi),%%ymm2,%%ymm2		\n\t"\
		"vmulpd	     (%%rdi),%%ymm3,%%ymm3		\n\t"\
		"vmovaps	     (%%rsi),%%ymm0	\n\t"\
		"vmovaps	0x020(%%rsi),%%ymm1	\n\t"\
		"vsubpd	      %%ymm2,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	      %%ymm3,%%ymm1,%%ymm1		\n\t"\
		"vaddpd	     (%%rsi),%%ymm2,%%ymm2	\n\t"\
		"vaddpd	0x020(%%rsi),%%ymm3,%%ymm3	\n\t"\
		"vsubpd	      %%ymm4,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	      %%ymm5,%%ymm1,%%ymm1		\n\t"\
		"vaddpd	      %%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	      %%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vmovaps	%%ymm0,     (%%rbx)	\n\t"\
		"vmovaps	%%ymm1,0x020(%%rbx)	\n\t"\
		"vaddpd	      %%ymm0,%%ymm4,%%ymm4	\n\t"\
		"vaddpd	      %%ymm1,%%ymm5,%%ymm5	\n\t"\
		"vmovaps	%%ymm4,     (%%rax)	\n\t"\
		"vmovaps	%%ymm5,0x020(%%rax)	\n\t"\
		"vsubpd	      %%ymm7,%%ymm2,%%ymm2	\n\t"\
		"vsubpd	      %%ymm6,%%ymm3,%%ymm3	\n\t"\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7	\n\t"\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6	\n\t"\
		"vmovaps	%%ymm2,     (%%rcx)	\n\t"\
		"vmovaps	%%ymm3,0x020(%%rdx)	\n\t"\
		"vaddpd	      %%ymm2,%%ymm7,%%ymm7	\n\t"\
		"vaddpd	      %%ymm3,%%ymm6,%%ymm6	\n\t"\
		"vmovaps	%%ymm7,     (%%rdx)	\n\t"\
		"vmovaps	%%ymm6,0x020(%%rcx)	\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__p1] "m" (Xp1)\
		 ,[__p2] "m" (Xp2)\
		 ,[__p3] "m" (Xp3)\
		 ,[__p4] "m" (Xp4)\
		 ,[__r1] "m" (Xr1)\
		 ,[__r3] "m" (Xr3)\
		 ,[__r5] "m" (Xr5)\
		 ,[__r7] "m" (Xr7)\
		 ,[__r9] "m" (Xr9)\
		 ,[__r11] "m" (Xr11)\
		 ,[__r13] "m" (Xr13)\
		 ,[__r15] "m" (Xr15)\
		 ,[__r17] "m" (Xr17)\
		 ,[__r19] "m" (Xr19)\
		 ,[__r21] "m" (Xr21)\
		 ,[__r23] "m" (Xr23)\
		 ,[__r25] "m" (Xr25)\
		 ,[__r27] "m" (Xr27)\
		 ,[__r29] "m" (Xr29)\
		 ,[__r31] "m" (Xr31)\
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__cc0] "m" (Xcc0)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
	);\
	}

#elif OS_BITS == 64	// 64-bit SSE2

	#define SSE2_RADIX16_DIT_NOTWIDDLE(Xadd0,Xp1,Xp2,Xp3,Xp4,Xr1,Xr3,Xr5,Xr7,Xr9,Xr11,Xr17,Xr25,Xisrt2,Xcc0)\
	{\
	__asm__ volatile (\
		"movq	%[__add0],%%rax		\n\t"\
		"movslq	%[__p1],%%rbx		\n\t"\
		"movslq	%[__p2],%%rcx		\n\t"\
		"movslq	%[__p3],%%rdx		\n\t"\
		"movslq	%[__p4],%%rdi		\n\t"\
		"shlq	$3,%%rbx			\n\t"\
		"shlq	$3,%%rcx			\n\t"\
		"shlq	$3,%%rdx			\n\t"\
		"shlq	$3,%%rdi			\n\t"\
		"addq	%%rax,%%rbx			\n\t"\
		"addq	%%rax,%%rcx			\n\t"\
		"addq	%%rax,%%rdx			\n\t"\
		"/* SSE2_RADIX4_DIT_0TWIDDLE_B(r1 ): */\n\t"\
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
		"movq	%[__r1], %%rsi		\n\t"\
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
		"addq	%%rdi,%%rax			\n\t"\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	%%rdi,%%rcx			\n\t"\
		"addq	%%rdi,%%rdx			\n\t"\
		"/* SSE2_RADIX4_DIT_0TWIDDLE_B(r9 ): */\n\t"\
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
		"movq	%[__r9], %%rsi		\n\t"\
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
		"addq	%%rdi,%%rax			\n\t"\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	%%rdi,%%rcx			\n\t"\
		"addq	%%rdi,%%rdx			\n\t"\
		"/* SSE2_RADIX4_DIT_0TWIDDLE_B(r17): */\n\t"\
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
		"movq	%[__r17], %%rsi		\n\t"\
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
		"addq	%%rdi,%%rax			\n\t"\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	%%rdi,%%rcx			\n\t"\
		"addq	%%rdi,%%rdx			\n\t"\
		"/* SSE2_RADIX4_DIT_0TWIDDLE_B(r25): */\n\t"\
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
		"movq	%[__r25], %%rsi		\n\t"\
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
		"movq	%[__r1],%%rax		\n\t"\
		"movaps	     (%%rax),%%xmm0	\n\t"\
		"movaps	0x010(%%rax),%%xmm1	\n\t"\
		"movaps	0x080(%%rax),%%xmm2	\n\t"\
		"movaps	0x090(%%rax),%%xmm3	\n\t"\
		"subpd	0x080(%%rax),%%xmm0	\n\t"\
		"subpd	0x090(%%rax),%%xmm1	\n\t"\
		"addpd	     (%%rax),%%xmm2	\n\t"\
		"addpd	0x010(%%rax),%%xmm3	\n\t"\
		"movaps	0x100(%%rax),%%xmm4	\n\t"\
		"movaps	0x110(%%rax),%%xmm5	\n\t"\
		"movaps	0x180(%%rax),%%xmm6	\n\t"\
		"movaps	0x190(%%rax),%%xmm7	\n\t"\
		"subpd	0x180(%%rax),%%xmm4	\n\t"\
		"subpd	0x190(%%rax),%%xmm5	\n\t"\
		"addpd	0x100(%%rax),%%xmm6	\n\t"\
		"addpd	0x110(%%rax),%%xmm7	\n\t"\
		"subpd	%%xmm6,%%xmm2		\n\t"\
		"subpd	%%xmm7,%%xmm3		\n\t"\
		"movaps	%%xmm2,0x100(%%rax)	\n\t"\
		"movaps	%%xmm3,0x110(%%rax)	\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm2,%%xmm6		\n\t"\
		"addpd	%%xmm3,%%xmm7		\n\t"\
		"movaps	%%xmm6,     (%%rax)	\n\t"\
		"movaps	%%xmm7,0x010(%%rax)	\n\t"\
		"subpd	%%xmm5,%%xmm0		\n\t"\
		"subpd	%%xmm4,%%xmm1		\n\t"\
		"movaps	%%xmm0,0x180(%%rax)	\n\t"\
		"movaps	%%xmm1,0x090(%%rax)	\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm0,%%xmm5		\n\t"\
		"addpd	%%xmm1,%%xmm4		\n\t"\
		"movaps	%%xmm5,0x080(%%rax)	\n\t"\
		"movaps	%%xmm4,0x190(%%rax)	\n\t"\
		"movq	%[__r5],%%rax		\n\t"\
		"movq	%[__isrt2],%%rbx	\n\t"\
		"movaps	(%%rbx),%%xmm2		\n\t"\
		"movaps	0x100(%%rax),%%xmm4	\n\t"\
		"movaps	0x110(%%rax),%%xmm5	\n\t"\
		"movaps	0x180(%%rax),%%xmm0	\n\t"\
		"movaps	0x190(%%rax),%%xmm1	\n\t"\
		"addpd	0x110(%%rax),%%xmm4	\n\t"\
		"subpd	0x100(%%rax),%%xmm5	\n\t"\
		"subpd	0x190(%%rax),%%xmm0	\n\t"\
		"addpd	0x180(%%rax),%%xmm1	\n\t"\
		"mulpd	%%xmm2,%%xmm4		\n\t"\
		"mulpd	%%xmm2,%%xmm5		\n\t"\
		"mulpd	%%xmm2,%%xmm0		\n\t"\
		"mulpd	%%xmm2,%%xmm1		\n\t"\
		"movaps	%%xmm4,%%xmm6		\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t"\
		"subpd	%%xmm0,%%xmm4		\n\t"\
		"subpd	%%xmm1,%%xmm5		\n\t"\
		"addpd	%%xmm0,%%xmm6		\n\t"\
		"addpd	%%xmm1,%%xmm7		\n\t"\
		"movaps	     (%%rax),%%xmm0	\n\t"\
		"movaps	0x010(%%rax),%%xmm1	\n\t"\
		"movaps	0x080(%%rax),%%xmm2	\n\t"\
		"movaps	0x090(%%rax),%%xmm3	\n\t"\
		"subpd	0x090(%%rax),%%xmm0	\n\t"\
		"subpd	0x080(%%rax),%%xmm1	\n\t"\
		"addpd	     (%%rax),%%xmm3	\n\t"\
		"addpd	0x010(%%rax),%%xmm2	\n\t"\
		"subpd	%%xmm4,%%xmm3		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"movaps	%%xmm3,0x100(%%rax)	\n\t"\
		"movaps	%%xmm1,0x110(%%rax)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm3,%%xmm4		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"movaps	%%xmm4,     (%%rax)	\n\t"\
		"movaps	%%xmm5,0x010(%%rax)	\n\t"\
		"subpd	%%xmm7,%%xmm0		\n\t"\
		"subpd	%%xmm6,%%xmm2		\n\t"\
		"movaps	%%xmm0,0x180(%%rax)	\n\t"\
		"movaps	%%xmm2,0x090(%%rax)	\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm7		\n\t"\
		"addpd	%%xmm2,%%xmm6		\n\t"\
		"movaps	%%xmm7,0x080(%%rax)	\n\t"\
		"movaps	%%xmm6,0x190(%%rax)	\n\t"\
		"movq	%[__r3],%%rax		\n\t"\
		"movq	%[__isrt2],%%rbx	\n\t"\
		"movq	%[__cc0],%%rcx		\n\t"\
		"movaps	0x180(%%rax),%%xmm0	\n\t"\
		"movaps	0x190(%%rax),%%xmm1	\n\t"\
		"movaps	0x180(%%rax),%%xmm2	\n\t"\
		"movaps	0x190(%%rax),%%xmm3	\n\t"\
		"mulpd	0x010(%%rcx),%%xmm0	\n\t"\
		"mulpd	0x010(%%rcx),%%xmm1	\n\t"\
		"mulpd	     (%%rcx),%%xmm2	\n\t"\
		"mulpd	     (%%rcx),%%xmm3	\n\t"\
		"subpd	%%xmm2,%%xmm1		\n\t"\
		"addpd	%%xmm3,%%xmm0		\n\t"\
		"movaps	0x100(%%rax),%%xmm4	\n\t"\
		"movaps	0x110(%%rax),%%xmm5	\n\t"\
		"movaps	0x100(%%rax),%%xmm6	\n\t"\
		"movaps	0x110(%%rax),%%xmm7	\n\t"\
		"mulpd	     (%%rcx),%%xmm4	\n\t"\
		"mulpd	     (%%rcx),%%xmm5	\n\t"\
		"mulpd	0x010(%%rcx),%%xmm6	\n\t"\
		"mulpd	0x010(%%rcx),%%xmm7	\n\t"\
		"subpd	%%xmm6,%%xmm5		\n\t"\
		"addpd	%%xmm7,%%xmm4		\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t"\
		"movaps	%%xmm4,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"subpd	%%xmm0,%%xmm6		\n\t"\
		"subpd	%%xmm1,%%xmm7		\n\t"\
		"movaps	0x080(%%rax),%%xmm2	\n\t"\
		"movaps	0x090(%%rax),%%xmm3	\n\t"\
		"movaps	     (%%rax),%%xmm0	\n\t"\
		"movaps	0x010(%%rax),%%xmm1	\n\t"\
		"addpd	0x090(%%rax),%%xmm2	\n\t"\
		"subpd	0x080(%%rax),%%xmm3	\n\t"\
		"mulpd	     (%%rbx),%%xmm2	\n\t"\
		"mulpd	     (%%rbx),%%xmm3	\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"subpd	%%xmm4,%%xmm2		\n\t"\
		"subpd	%%xmm5,%%xmm3		\n\t"\
		"movaps	%%xmm2,0x100(%%rax)	\n\t"\
		"movaps	%%xmm3,0x110(%%rax)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm2,%%xmm4		\n\t"\
		"addpd	%%xmm3,%%xmm5		\n\t"\
		"movaps	%%xmm4,     (%%rax)	\n\t"\
		"movaps	%%xmm5,0x010(%%rax)	\n\t"\
		"subpd	%%xmm7,%%xmm0		\n\t"\
		"subpd	%%xmm6,%%xmm1		\n\t"\
		"movaps	%%xmm0,0x180(%%rax)	\n\t"\
		"movaps	%%xmm1,0x090(%%rax)	\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm7		\n\t"\
		"addpd	%%xmm1,%%xmm6		\n\t"\
		"movaps	%%xmm7,0x080(%%rax)	\n\t"\
		"movaps	%%xmm6,0x190(%%rax)	\n\t"\
		"movq	%[__r7],%%rax		\n\t"\
		"movq	%[__isrt2],%%rbx	\n\t"\
		"movq	%[__cc0],%%rcx		\n\t"\
		"movaps	0x180(%%rax),%%xmm0	\n\t"\
		"movaps	0x190(%%rax),%%xmm1	\n\t"\
		"movaps	0x180(%%rax),%%xmm2	\n\t"\
		"movaps	0x190(%%rax),%%xmm3	\n\t"\
		"mulpd	     (%%rcx),%%xmm0	\n\t"\
		"mulpd	     (%%rcx),%%xmm1	\n\t"\
		"mulpd	0x010(%%rcx),%%xmm2	\n\t"\
		"mulpd	0x010(%%rcx),%%xmm3	\n\t"\
		"subpd	%%xmm2,%%xmm1		\n\t"\
		"addpd	%%xmm3,%%xmm0		\n\t"\
		"movaps	0x100(%%rax),%%xmm4	\n\t"\
		"movaps	0x110(%%rax),%%xmm5	\n\t"\
		"movaps	0x100(%%rax),%%xmm6	\n\t"\
		"movaps	0x110(%%rax),%%xmm7	\n\t"\
		"mulpd	0x010(%%rcx),%%xmm4	\n\t"\
		"mulpd	0x010(%%rcx),%%xmm5	\n\t"\
		"mulpd	     (%%rcx),%%xmm6	\n\t"\
		"mulpd	     (%%rcx),%%xmm7	\n\t"\
		"subpd	%%xmm6,%%xmm5		\n\t"\
		"addpd	%%xmm7,%%xmm4		\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t"\
		"movaps	%%xmm4,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"subpd	%%xmm0,%%xmm6		\n\t"\
		"subpd	%%xmm1,%%xmm7		\n\t"\
		"movaps	0x080(%%rax),%%xmm2	\n\t"\
		"movaps	0x090(%%rax),%%xmm3	\n\t"\
		"movaps	     (%%rax),%%xmm0	\n\t"\
		"movaps	0x010(%%rax),%%xmm1	\n\t"\
		"subpd	0x090(%%rax),%%xmm2	\n\t"\
		"addpd	0x080(%%rax),%%xmm3	\n\t"\
		"mulpd	     (%%rbx),%%xmm2	\n\t"\
		"mulpd	     (%%rbx),%%xmm3	\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"subpd	%%xmm6,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm1		\n\t"\
		"movaps	%%xmm0,0x100(%%rax)	\n\t"\
		"movaps	%%xmm1,0x110(%%rax)	\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm0,%%xmm6		\n\t"\
		"addpd	%%xmm1,%%xmm7		\n\t"\
		"movaps	%%xmm6,     (%%rax)	\n\t"\
		"movaps	%%xmm7,0x010(%%rax)	\n\t"\
		"subpd	%%xmm5,%%xmm2		\n\t"\
		"subpd	%%xmm4,%%xmm3		\n\t"\
		"movaps	%%xmm2,0x180(%%rax)	\n\t"\
		"movaps	%%xmm3,0x090(%%rax)	\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm2,%%xmm5		\n\t"\
		"addpd	%%xmm3,%%xmm4		\n\t"\
		"movaps	%%xmm5,0x080(%%rax)	\n\t"\
		"movaps	%%xmm4,0x190(%%rax)	\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__p1] "m" (Xp1)\
		 ,[__p2] "m" (Xp2)\
		 ,[__p3] "m" (Xp3)\
		 ,[__p4] "m" (Xp4)\
		 ,[__r1] "m" (Xr1)\
		 ,[__r3] "m" (Xr3)\
		 ,[__r5] "m" (Xr5)\
		 ,[__r7] "m" (Xr7)\
		 ,[__r9] "m" (Xr9)\
		 ,[__r17] "m" (Xr17)\
		 ,[__r25] "m" (Xr25)\
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__cc0] "m" (Xcc0)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
	);\
	}

	#define SSE2_RADIX16_DIF_NOTWIDDLE(Xadd0,Xp1,Xp2,Xp3,Xp4,Xr1,Xr3,Xr5,Xr7,Xr9,Xr11,Xr13,Xr15,Xr17,Xr19,Xr21,Xr23,Xr25,Xr27,Xr29,Xr31,Xisrt2,Xcc0)\
	{\
	__asm__ volatile (\
		"/* SSE2_RADIX4_DIF_IN_PLACE(r1 , r17, r9 , r25): */\n\t"\
		"movq	%[__r1] ,%%rax		\n\t"\
		"movq	%[__r17],%%rbx		\n\t"\
		"movq	%[__r9] ,%%rcx		\n\t"\
		"movq	%[__r25],%%rdx		\n\t"\
		"movaps	    (%%rax),%%xmm0	\n\t"\
		"movaps	0x10(%%rax),%%xmm1	\n\t"\
		"movaps	    (%%rax),%%xmm2	\n\t"\
		"movaps	0x10(%%rax),%%xmm3	\n\t"\
		"addpd	    (%%rbx),%%xmm0	\n\t"\
		"addpd	0x10(%%rbx),%%xmm1	\n\t"\
		"subpd	    (%%rbx),%%xmm2	\n\t"\
		"subpd	0x10(%%rbx),%%xmm3	\n\t"\
		"movaps	    (%%rcx),%%xmm4	\n\t"\
		"movaps	0x10(%%rcx),%%xmm5	\n\t"\
		"movaps	    (%%rcx),%%xmm6	\n\t"\
		"movaps	0x10(%%rcx),%%xmm7	\n\t"\
		"addpd	    (%%rdx),%%xmm4	\n\t"\
		"addpd	0x10(%%rdx),%%xmm5	\n\t"\
		"subpd	    (%%rdx),%%xmm6	\n\t"\
		"subpd	0x10(%%rdx),%%xmm7	\n\t"\
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
		"/* SSE2_RADIX4_DIF_IN_PLACE(r5 , r21, r13, r29): */\n\t"\
		"movq	%[__r5] ,%%rax		\n\t"\
		"movq	%[__r21],%%rbx		\n\t"\
		"movq	%[__r13],%%rcx		\n\t"\
		"movq	%[__r29],%%rdx		\n\t"\
		"movaps	    (%%rax),%%xmm0	\n\t"\
		"movaps	0x10(%%rax),%%xmm1	\n\t"\
		"movaps	    (%%rax),%%xmm2	\n\t"\
		"movaps	0x10(%%rax),%%xmm3	\n\t"\
		"addpd	    (%%rbx),%%xmm0	\n\t"\
		"addpd	0x10(%%rbx),%%xmm1	\n\t"\
		"subpd	    (%%rbx),%%xmm2	\n\t"\
		"subpd	0x10(%%rbx),%%xmm3	\n\t"\
		"movaps	    (%%rcx),%%xmm4	\n\t"\
		"movaps	0x10(%%rcx),%%xmm5	\n\t"\
		"movaps	    (%%rcx),%%xmm6	\n\t"\
		"movaps	0x10(%%rcx),%%xmm7	\n\t"\
		"addpd	    (%%rdx),%%xmm4	\n\t"\
		"addpd	0x10(%%rdx),%%xmm5	\n\t"\
		"subpd	    (%%rdx),%%xmm6	\n\t"\
		"subpd	0x10(%%rdx),%%xmm7	\n\t"\
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
		"/* SSE2_RADIX4_DIF_IN_PLACE(r3 , r19, r11, r27): */\n\t"\
		"movq	%[__r3] ,%%rax		\n\t"\
		"movq	%[__r19],%%rbx		\n\t"\
		"movq	%[__r11],%%rcx		\n\t"\
		"movq	%[__r27],%%rdx		\n\t"\
		"movaps	    (%%rax),%%xmm0	\n\t"\
		"movaps	0x10(%%rax),%%xmm1	\n\t"\
		"movaps	    (%%rax),%%xmm2	\n\t"\
		"movaps	0x10(%%rax),%%xmm3	\n\t"\
		"addpd	    (%%rbx),%%xmm0	\n\t"\
		"addpd	0x10(%%rbx),%%xmm1	\n\t"\
		"subpd	    (%%rbx),%%xmm2	\n\t"\
		"subpd	0x10(%%rbx),%%xmm3	\n\t"\
		"movaps	    (%%rcx),%%xmm4	\n\t"\
		"movaps	0x10(%%rcx),%%xmm5	\n\t"\
		"movaps	    (%%rcx),%%xmm6	\n\t"\
		"movaps	0x10(%%rcx),%%xmm7	\n\t"\
		"addpd	    (%%rdx),%%xmm4	\n\t"\
		"addpd	0x10(%%rdx),%%xmm5	\n\t"\
		"subpd	    (%%rdx),%%xmm6	\n\t"\
		"subpd	0x10(%%rdx),%%xmm7	\n\t"\
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
		"/* SSE2_RADIX4_DIF_IN_PLACE(r7 , r23, r15, r31): */\n\t"\
		"movq	%[__r7] ,%%rax		\n\t"\
		"movq	%[__r23],%%rbx		\n\t"\
		"movq	%[__r15],%%rcx		\n\t"\
		"movq	%[__r31],%%rdx		\n\t"\
		"movaps	    (%%rax),%%xmm0	\n\t"\
		"movaps	0x10(%%rax),%%xmm1	\n\t"\
		"movaps	    (%%rax),%%xmm2	\n\t"\
		"movaps	0x10(%%rax),%%xmm3	\n\t"\
		"addpd	    (%%rbx),%%xmm0	\n\t"\
		"addpd	0x10(%%rbx),%%xmm1	\n\t"\
		"subpd	    (%%rbx),%%xmm2	\n\t"\
		"subpd	0x10(%%rbx),%%xmm3	\n\t"\
		"movaps	    (%%rcx),%%xmm4	\n\t"\
		"movaps	0x10(%%rcx),%%xmm5	\n\t"\
		"movaps	    (%%rcx),%%xmm6	\n\t"\
		"movaps	0x10(%%rcx),%%xmm7	\n\t"\
		"addpd	    (%%rdx),%%xmm4	\n\t"\
		"addpd	0x10(%%rdx),%%xmm5	\n\t"\
		"subpd	    (%%rdx),%%xmm6	\n\t"\
		"subpd	0x10(%%rdx),%%xmm7	\n\t"\
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
		"movq	%[__add0],%%rax		\n\t"\
		"movslq	%[__p1],%%rbx		\n\t"\
		"movslq	%[__p2],%%rcx		\n\t"\
		"movslq	%[__p3],%%rdx		\n\t"\
		"movslq	%[__p4],%%rdi		\n\t"\
		"shlq	$3,%%rbx			\n\t"\
		"shlq	$3,%%rcx			\n\t"\
		"shlq	$3,%%rdx			\n\t"\
		"shlq	$3,%%rdi			\n\t"\
		"addq	%%rax,%%rbx			\n\t"\
		"addq	%%rax,%%rcx			\n\t"\
		"addq	%%rax,%%rdx			\n\t"\
		"movq	%[__r1],%%rsi		\n\t"\
		"movaps	     (%%rsi),%%xmm0	\n\t"\
		"movaps	0x010(%%rsi),%%xmm1	\n\t"\
		"movaps	0x040(%%rsi),%%xmm2	\n\t"\
		"movaps	0x050(%%rsi),%%xmm3	\n\t"\
		"subpd	0x040(%%rsi),%%xmm0	\n\t"\
		"subpd	0x050(%%rsi),%%xmm1	\n\t"\
		"addpd	     (%%rsi),%%xmm2	\n\t"\
		"addpd	0x010(%%rsi),%%xmm3	\n\t"\
		"movaps	0x020(%%rsi),%%xmm4	\n\t"\
		"movaps	0x030(%%rsi),%%xmm5	\n\t"\
		"movaps	0x060(%%rsi),%%xmm6	\n\t"\
		"movaps	0x070(%%rsi),%%xmm7	\n\t"\
		"subpd	0x060(%%rsi),%%xmm4	\n\t"\
		"subpd	0x070(%%rsi),%%xmm5	\n\t"\
		"addpd	0x020(%%rsi),%%xmm6	\n\t"\
		"addpd	0x030(%%rsi),%%xmm7	\n\t"\
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
		"addq	%%rdi,%%rax			\n\t"\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	%%rdi,%%rcx			\n\t"\
		"addq	%%rdi,%%rdx			\n\t"\
		"movq	%[__r17],%%rsi		\n\t"\
		"movaps	     (%%rsi),%%xmm0	\n\t"\
		"movaps	0x010(%%rsi),%%xmm1	\n\t"\
		"movaps	0x040(%%rsi),%%xmm2	\n\t"\
		"movaps	0x050(%%rsi),%%xmm3	\n\t"\
		"subpd	0x050(%%rsi),%%xmm0	\n\t"\
		"subpd	0x040(%%rsi),%%xmm1	\n\t"\
		"addpd	0x010(%%rsi),%%xmm2	\n\t"\
		"addpd	     (%%rsi),%%xmm3	\n\t"\
		"movaps	0x020(%%rsi),%%xmm4	\n\t"\
		"movaps	0x030(%%rsi),%%xmm5	\n\t"\
		"movaps	0x060(%%rsi),%%xmm6	\n\t"\
		"movaps	0x070(%%rsi),%%xmm7	\n\t"\
		"subpd	0x030(%%rsi),%%xmm4	\n\t"\
		"addpd	0x020(%%rsi),%%xmm5	\n\t"\
		"addpd	0x070(%%rsi),%%xmm6	\n\t"\
		"subpd	0x060(%%rsi),%%xmm7	\n\t"\
		"movq	%[__isrt2],%%rsi	\n\t"\
		"mulpd	(%%rsi),%%xmm4		\n\t"\
		"mulpd	(%%rsi),%%xmm5		\n\t"\
		"mulpd	(%%rsi),%%xmm6		\n\t"\
		"mulpd	(%%rsi),%%xmm7		\n\t"\
		"subpd	%%xmm6,%%xmm4		\n\t"\
		"subpd	%%xmm7,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm6		\n\t"\
		"addpd	%%xmm5,%%xmm7		\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm2		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"movaps	%%xmm0,    (%%rbx)	\n\t"\
		"movaps	%%xmm2,0x10(%%rbx)	\n\t"\
		"addpd	%%xmm0,		%%xmm4	\n\t"\
		"addpd	%%xmm2,		%%xmm5	\n\t"\
		"movaps	%%xmm4,    (%%rax)	\n\t"\
		"movaps	%%xmm5,0x10(%%rax)	\n\t"\
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
		"addq	%%rdi,%%rax			\n\t"\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	%%rdi,%%rcx			\n\t"\
		"addq	%%rdi,%%rdx			\n\t"\
		"movq	%[__r9] ,%%rsi		\n\t"\
		"movq	%[__cc0],%%rdi		\n\t"\
		"movaps	0x020(%%rsi),%%xmm4	\n\t"\
		"movaps	0x060(%%rsi),%%xmm6	\n\t"\
		"movaps	0x030(%%rsi),%%xmm5	\n\t"\
		"movaps	0x070(%%rsi),%%xmm7	\n\t"\
		"movaps	%%xmm4,%%xmm0		\n\t"\
		"movaps	%%xmm6,%%xmm2		\n\t"\
		"movaps	%%xmm5,%%xmm1		\n\t"\
		"movaps	%%xmm7,%%xmm3		\n\t"\
		"mulpd	     (%%rdi),%%xmm4	\n\t"\
		"mulpd	0x010(%%rdi),%%xmm6	\n\t"\
		"mulpd	0x010(%%rdi),%%xmm1	\n\t"\
		"mulpd	     (%%rdi),%%xmm3	\n\t"\
		"mulpd	     (%%rdi),%%xmm5	\n\t"\
		"mulpd	0x010(%%rdi),%%xmm7	\n\t"\
		"mulpd	0x010(%%rdi),%%xmm0	\n\t"\
		"mulpd	     (%%rdi),%%xmm2	\n\t"\
		"subpd	%%xmm1,%%xmm4		\n\t"\
		"subpd	%%xmm3,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm5		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"subpd	%%xmm6,%%xmm4		\n\t"\
		"subpd	%%xmm7,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm6		\n\t"\
		"addpd	%%xmm5,%%xmm7		\n\t"\
		"movq	%[__isrt2],%%rdi	\n\t"\
		"movaps	0x040(%%rsi),%%xmm2	\n\t"\
		"movaps	0x050(%%rsi),%%xmm3	\n\t"\
		"subpd	0x050(%%rsi),%%xmm2	\n\t"\
		"addpd	0x040(%%rsi),%%xmm3	\n\t"\
		"mulpd	(%%rdi),%%xmm2		\n\t"\
		"mulpd	(%%rdi),%%xmm3		\n\t"\
		"movaps	     (%%rsi),%%xmm0	\n\t"\
		"movaps	0x010(%%rsi),%%xmm1	\n\t"\
		"movslq	%[__p4],%%rdi		\n\t"\
		"shlq	$3,%%rdi			\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"addpd	     (%%rsi),%%xmm2	\n\t"\
		"addpd	0x010(%%rsi),%%xmm3	\n\t"\
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
		"addq	%%rdi,%%rax			\n\t"\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	%%rdi,%%rcx			\n\t"\
		"addq	%%rdi,%%rdx			\n\t"\
		"movq	%[__r25],%%rsi		\n\t"\
		"movq	%[__cc0],%%rdi		\n\t"\
		"movaps	0x020(%%rsi),%%xmm4	\n\t"\
		"movaps	0x060(%%rsi),%%xmm6	\n\t"\
		"movaps	0x030(%%rsi),%%xmm5	\n\t"\
		"movaps	0x070(%%rsi),%%xmm7	\n\t"\
		"movaps	%%xmm4,%%xmm0		\n\t"\
		"movaps	%%xmm6,%%xmm2		\n\t"\
		"movaps	%%xmm5,%%xmm1		\n\t"\
		"movaps	%%xmm7,%%xmm3		\n\t"\
		"mulpd	0x010(%%rdi),%%xmm4	\n\t"\
		"mulpd	     (%%rdi),%%xmm6	\n\t"\
		"mulpd	     (%%rdi),%%xmm1	\n\t"\
		"mulpd	0x010(%%rdi),%%xmm3	\n\t"\
		"mulpd	0x010(%%rdi),%%xmm5	\n\t"\
		"mulpd	     (%%rdi),%%xmm7	\n\t"\
		"mulpd	     (%%rdi),%%xmm0	\n\t"\
		"mulpd	0x010(%%rdi),%%xmm2	\n\t"\
		"subpd	%%xmm1,%%xmm4		\n\t"\
		"subpd	%%xmm3,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm5		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"subpd	%%xmm6,%%xmm4		\n\t"\
		"subpd	%%xmm7,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm6		\n\t"\
		"addpd	%%xmm5,%%xmm7		\n\t"\
		"movq	%[__isrt2],%%rdi	\n\t"\
		"movaps	0x040(%%rsi),%%xmm2	\n\t"\
		"movaps	0x050(%%rsi),%%xmm3	\n\t"\
		"addpd	0x050(%%rsi),%%xmm2	\n\t"\
		"subpd	0x040(%%rsi),%%xmm3	\n\t"\
		"mulpd	(%%rdi),%%xmm2		\n\t"\
		"mulpd	(%%rdi),%%xmm3		\n\t"\
		"movaps	     (%%rsi),%%xmm0	\n\t"\
		"movaps	0x010(%%rsi),%%xmm1	\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"addpd	     (%%rsi),%%xmm2	\n\t"\
		"addpd	0x010(%%rsi),%%xmm3	\n\t"\
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
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__p1] "m" (Xp1)\
		 ,[__p2] "m" (Xp2)\
		 ,[__p3] "m" (Xp3)\
		 ,[__p4] "m" (Xp4)\
		 ,[__r1] "m" (Xr1)\
		 ,[__r3] "m" (Xr3)\
		 ,[__r5] "m" (Xr5)\
		 ,[__r7] "m" (Xr7)\
		 ,[__r9] "m" (Xr9)\
		 ,[__r11] "m" (Xr11)\
		 ,[__r13] "m" (Xr13)\
		 ,[__r15] "m" (Xr15)\
		 ,[__r17] "m" (Xr17)\
		 ,[__r19] "m" (Xr19)\
		 ,[__r21] "m" (Xr21)\
		 ,[__r23] "m" (Xr23)\
		 ,[__r25] "m" (Xr25)\
		 ,[__r27] "m" (Xr27)\
		 ,[__r29] "m" (Xr29)\
		 ,[__r31] "m" (Xr31)\
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__cc0] "m" (Xcc0)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
	);\
	}

#elif OS_BITS == 32	// 32-bit SSE2

	#define SSE2_RADIX16_DIT_NOTWIDDLE(Xadd0,Xp1,Xp2,Xp3,Xp4,Xr1,Xr3,Xr5,Xr7,Xr9,Xr11,Xr17,Xr25,Xisrt2,Xcc0)\
	{\
	__asm__ volatile (\
	"pushl %%ebx	\n\t"/* Explicit save/restore of PIC register */\
		"movl	%[__add0],%%eax		\n\t"\
		"movl	%[__p1],%%ebx		\n\t"\
		"movl	%[__p2],%%ecx		\n\t"\
		"movl	%[__p3],%%edx		\n\t"\
		"movl	%[__p4],%%edi		\n\t"\
		"shll	$3,%%ebx			\n\t"\
		"shll	$3,%%ecx			\n\t"\
		"shll	$3,%%edx			\n\t"\
		"shll	$3,%%edi			\n\t"\
		"addl	%%eax,%%ebx			\n\t"\
		"addl	%%eax,%%ecx			\n\t"\
		"addl	%%eax,%%edx			\n\t"\
		"/* SSE2_RADIX4_DIT_0TWIDDLE_B(r1 ): */\n\t"\
		"movaps	    (%%eax),%%xmm2	\n\t"\
		"movaps	    (%%ecx),%%xmm6	\n\t"\
		"movaps	0x10(%%eax),%%xmm3	\n\t"\
		"movaps	0x10(%%ecx),%%xmm7	\n\t"\
		"movaps	    (%%ebx),%%xmm0	\n\t"\
		"movaps	    (%%edx),%%xmm4	\n\t"\
		"movaps	0x10(%%ebx),%%xmm1	\n\t"\
		"movaps	0x10(%%edx),%%xmm5	\n\t"\
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
		"movl	%[__r1], %%esi		\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm2		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"subpd	%%xmm6,%%xmm3		\n\t"\
		"movaps	%%xmm0,0x040(%%esi)	\n\t"\
		"movaps	%%xmm2,0x060(%%esi)	\n\t"\
		"movaps	%%xmm1,0x050(%%esi)	\n\t"\
		"movaps	%%xmm3,0x030(%%esi)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"addpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	%%xmm4,     (%%esi)	\n\t"\
		"movaps	%%xmm7,0x020(%%esi)	\n\t"\
		"movaps	%%xmm5,0x010(%%esi)	\n\t"\
		"movaps	%%xmm6,0x070(%%esi)	\n\t"\
		"addl	%%edi,%%eax			\n\t"\
		"addl	%%edi,%%ebx			\n\t"\
		"addl	%%edi,%%ecx			\n\t"\
		"addl	%%edi,%%edx			\n\t"\
		"/* SSE2_RADIX4_DIT_0TWIDDLE_B(r9 ): */\n\t"\
		"movaps	    (%%eax),%%xmm2	\n\t"\
		"movaps	    (%%ecx),%%xmm6	\n\t"\
		"movaps	0x10(%%eax),%%xmm3	\n\t"\
		"movaps	0x10(%%ecx),%%xmm7	\n\t"\
		"movaps	    (%%ebx),%%xmm0	\n\t"\
		"movaps	    (%%edx),%%xmm4	\n\t"\
		"movaps	0x10(%%ebx),%%xmm1	\n\t"\
		"movaps	0x10(%%edx),%%xmm5	\n\t"\
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
		"movl	%[__r9], %%esi		\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm2		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"subpd	%%xmm6,%%xmm3		\n\t"\
		"movaps	%%xmm0,0x040(%%esi)	\n\t"\
		"movaps	%%xmm2,0x060(%%esi)	\n\t"\
		"movaps	%%xmm1,0x050(%%esi)	\n\t"\
		"movaps	%%xmm3,0x030(%%esi)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"addpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	%%xmm4,     (%%esi)	\n\t"\
		"movaps	%%xmm7,0x020(%%esi)	\n\t"\
		"movaps	%%xmm5,0x010(%%esi)	\n\t"\
		"movaps	%%xmm6,0x070(%%esi)	\n\t"\
		"addl	%%edi,%%eax			\n\t"\
		"addl	%%edi,%%ebx			\n\t"\
		"addl	%%edi,%%ecx			\n\t"\
		"addl	%%edi,%%edx			\n\t"\
		"/* SSE2_RADIX4_DIT_0TWIDDLE_B(r17): */\n\t"\
		"movaps	    (%%eax),%%xmm2	\n\t"\
		"movaps	    (%%ecx),%%xmm6	\n\t"\
		"movaps	0x10(%%eax),%%xmm3	\n\t"\
		"movaps	0x10(%%ecx),%%xmm7	\n\t"\
		"movaps	    (%%ebx),%%xmm0	\n\t"\
		"movaps	    (%%edx),%%xmm4	\n\t"\
		"movaps	0x10(%%ebx),%%xmm1	\n\t"\
		"movaps	0x10(%%edx),%%xmm5	\n\t"\
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
		"movl	%[__r17], %%esi		\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm2		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"subpd	%%xmm6,%%xmm3		\n\t"\
		"movaps	%%xmm0,0x040(%%esi)	\n\t"\
		"movaps	%%xmm2,0x060(%%esi)	\n\t"\
		"movaps	%%xmm1,0x050(%%esi)	\n\t"\
		"movaps	%%xmm3,0x030(%%esi)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"addpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	%%xmm4,     (%%esi)	\n\t"\
		"movaps	%%xmm7,0x020(%%esi)	\n\t"\
		"movaps	%%xmm5,0x010(%%esi)	\n\t"\
		"movaps	%%xmm6,0x070(%%esi)	\n\t"\
		"addl	%%edi,%%eax			\n\t"\
		"addl	%%edi,%%ebx			\n\t"\
		"addl	%%edi,%%ecx			\n\t"\
		"addl	%%edi,%%edx			\n\t"\
		"/* SSE2_RADIX4_DIT_0TWIDDLE_B(r25): */\n\t"\
		"movaps	    (%%eax),%%xmm2	\n\t"\
		"movaps	    (%%ecx),%%xmm6	\n\t"\
		"movaps	0x10(%%eax),%%xmm3	\n\t"\
		"movaps	0x10(%%ecx),%%xmm7	\n\t"\
		"movaps	    (%%ebx),%%xmm0	\n\t"\
		"movaps	    (%%edx),%%xmm4	\n\t"\
		"movaps	0x10(%%ebx),%%xmm1	\n\t"\
		"movaps	0x10(%%edx),%%xmm5	\n\t"\
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
		"movl	%[__r25], %%esi		\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm2		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"subpd	%%xmm6,%%xmm3		\n\t"\
		"movaps	%%xmm0,0x040(%%esi)	\n\t"\
		"movaps	%%xmm2,0x060(%%esi)	\n\t"\
		"movaps	%%xmm1,0x050(%%esi)	\n\t"\
		"movaps	%%xmm3,0x030(%%esi)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"addpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	%%xmm4,     (%%esi)	\n\t"\
		"movaps	%%xmm7,0x020(%%esi)	\n\t"\
		"movaps	%%xmm5,0x010(%%esi)	\n\t"\
		"movaps	%%xmm6,0x070(%%esi)	\n\t"\
		"movl	%[__r1],%%eax		\n\t"\
		"movaps	     (%%eax),%%xmm0	\n\t"\
		"movaps	0x010(%%eax),%%xmm1	\n\t"\
		"movaps	0x080(%%eax),%%xmm2	\n\t"\
		"movaps	0x090(%%eax),%%xmm3	\n\t"\
		"subpd	0x080(%%eax),%%xmm0	\n\t"\
		"subpd	0x090(%%eax),%%xmm1	\n\t"\
		"addpd	     (%%eax),%%xmm2	\n\t"\
		"addpd	0x010(%%eax),%%xmm3	\n\t"\
		"movaps	0x100(%%eax),%%xmm4	\n\t"\
		"movaps	0x110(%%eax),%%xmm5	\n\t"\
		"movaps	0x180(%%eax),%%xmm6	\n\t"\
		"movaps	0x190(%%eax),%%xmm7	\n\t"\
		"subpd	0x180(%%eax),%%xmm4	\n\t"\
		"subpd	0x190(%%eax),%%xmm5	\n\t"\
		"addpd	0x100(%%eax),%%xmm6	\n\t"\
		"addpd	0x110(%%eax),%%xmm7	\n\t"\
		"subpd	%%xmm6,%%xmm2		\n\t"\
		"subpd	%%xmm7,%%xmm3		\n\t"\
		"movaps	%%xmm2,0x100(%%eax)	\n\t"\
		"movaps	%%xmm3,0x110(%%eax)	\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm2,%%xmm6		\n\t"\
		"addpd	%%xmm3,%%xmm7		\n\t"\
		"movaps	%%xmm6,     (%%eax)	\n\t"\
		"movaps	%%xmm7,0x010(%%eax)	\n\t"\
		"subpd	%%xmm5,%%xmm0		\n\t"\
		"subpd	%%xmm4,%%xmm1		\n\t"\
		"movaps	%%xmm0,0x180(%%eax)	\n\t"\
		"movaps	%%xmm1,0x090(%%eax)	\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm0,%%xmm5		\n\t"\
		"addpd	%%xmm1,%%xmm4		\n\t"\
		"movaps	%%xmm5,0x080(%%eax)	\n\t"\
		"movaps	%%xmm4,0x190(%%eax)	\n\t"\
		"movl	%[__r5],%%eax		\n\t"\
		"movl	%[__isrt2],%%ebx	\n\t"\
		"movaps	(%%ebx),%%xmm2		\n\t"\
		"movaps	0x100(%%eax),%%xmm4	\n\t"\
		"movaps	0x110(%%eax),%%xmm5	\n\t"\
		"movaps	0x180(%%eax),%%xmm0	\n\t"\
		"movaps	0x190(%%eax),%%xmm1	\n\t"\
		"addpd	0x110(%%eax),%%xmm4	\n\t"\
		"subpd	0x100(%%eax),%%xmm5	\n\t"\
		"subpd	0x190(%%eax),%%xmm0	\n\t"\
		"addpd	0x180(%%eax),%%xmm1	\n\t"\
		"mulpd	%%xmm2,%%xmm4		\n\t"\
		"mulpd	%%xmm2,%%xmm5		\n\t"\
		"mulpd	%%xmm2,%%xmm0		\n\t"\
		"mulpd	%%xmm2,%%xmm1		\n\t"\
		"movaps	%%xmm4,%%xmm6		\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t"\
		"subpd	%%xmm0,%%xmm4		\n\t"\
		"subpd	%%xmm1,%%xmm5		\n\t"\
		"addpd	%%xmm0,%%xmm6		\n\t"\
		"addpd	%%xmm1,%%xmm7		\n\t"\
		"movaps	     (%%eax),%%xmm0	\n\t"\
		"movaps	0x010(%%eax),%%xmm1	\n\t"\
		"movaps	0x080(%%eax),%%xmm2	\n\t"\
		"movaps	0x090(%%eax),%%xmm3	\n\t"\
		"subpd	0x090(%%eax),%%xmm0	\n\t"\
		"subpd	0x080(%%eax),%%xmm1	\n\t"\
		"addpd	     (%%eax),%%xmm3	\n\t"\
		"addpd	0x010(%%eax),%%xmm2	\n\t"\
		"subpd	%%xmm4,%%xmm3		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"movaps	%%xmm3,0x100(%%eax)	\n\t"\
		"movaps	%%xmm1,0x110(%%eax)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm3,%%xmm4		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"movaps	%%xmm4,     (%%eax)	\n\t"\
		"movaps	%%xmm5,0x010(%%eax)	\n\t"\
		"subpd	%%xmm7,%%xmm0		\n\t"\
		"subpd	%%xmm6,%%xmm2		\n\t"\
		"movaps	%%xmm0,0x180(%%eax)	\n\t"\
		"movaps	%%xmm2,0x090(%%eax)	\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm7		\n\t"\
		"addpd	%%xmm2,%%xmm6		\n\t"\
		"movaps	%%xmm7,0x080(%%eax)	\n\t"\
		"movaps	%%xmm6,0x190(%%eax)	\n\t"\
		"movl	%[__r3],%%eax		\n\t"\
		"movl	%[__isrt2],%%ebx	\n\t"\
		"movl	%[__cc0],%%ecx		\n\t"\
		"movaps	0x180(%%eax),%%xmm0	\n\t"\
		"movaps	0x190(%%eax),%%xmm1	\n\t"\
		"movaps	0x180(%%eax),%%xmm2	\n\t"\
		"movaps	0x190(%%eax),%%xmm3	\n\t"\
		"mulpd	0x010(%%ecx),%%xmm0	\n\t"\
		"mulpd	0x010(%%ecx),%%xmm1	\n\t"\
		"mulpd	     (%%ecx),%%xmm2	\n\t"\
		"mulpd	     (%%ecx),%%xmm3	\n\t"\
		"subpd	%%xmm2,%%xmm1		\n\t"\
		"addpd	%%xmm3,%%xmm0		\n\t"\
		"movaps	0x100(%%eax),%%xmm4	\n\t"\
		"movaps	0x110(%%eax),%%xmm5	\n\t"\
		"movaps	0x100(%%eax),%%xmm6	\n\t"\
		"movaps	0x110(%%eax),%%xmm7	\n\t"\
		"mulpd	     (%%ecx),%%xmm4	\n\t"\
		"mulpd	     (%%ecx),%%xmm5	\n\t"\
		"mulpd	0x010(%%ecx),%%xmm6	\n\t"\
		"mulpd	0x010(%%ecx),%%xmm7	\n\t"\
		"subpd	%%xmm6,%%xmm5		\n\t"\
		"addpd	%%xmm7,%%xmm4		\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t"\
		"movaps	%%xmm4,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"subpd	%%xmm0,%%xmm6		\n\t"\
		"subpd	%%xmm1,%%xmm7		\n\t"\
		"movaps	0x080(%%eax),%%xmm2	\n\t"\
		"movaps	0x090(%%eax),%%xmm3	\n\t"\
		"movaps	     (%%eax),%%xmm0	\n\t"\
		"movaps	0x010(%%eax),%%xmm1	\n\t"\
		"addpd	0x090(%%eax),%%xmm2	\n\t"\
		"subpd	0x080(%%eax),%%xmm3	\n\t"\
		"mulpd	     (%%ebx),%%xmm2	\n\t"\
		"mulpd	     (%%ebx),%%xmm3	\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"subpd	%%xmm4,%%xmm2		\n\t"\
		"subpd	%%xmm5,%%xmm3		\n\t"\
		"movaps	%%xmm2,0x100(%%eax)	\n\t"\
		"movaps	%%xmm3,0x110(%%eax)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm2,%%xmm4		\n\t"\
		"addpd	%%xmm3,%%xmm5		\n\t"\
		"movaps	%%xmm4,     (%%eax)	\n\t"\
		"movaps	%%xmm5,0x010(%%eax)	\n\t"\
		"subpd	%%xmm7,%%xmm0		\n\t"\
		"subpd	%%xmm6,%%xmm1		\n\t"\
		"movaps	%%xmm0,0x180(%%eax)	\n\t"\
		"movaps	%%xmm1,0x090(%%eax)	\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm7		\n\t"\
		"addpd	%%xmm1,%%xmm6		\n\t"\
		"movaps	%%xmm7,0x080(%%eax)	\n\t"\
		"movaps	%%xmm6,0x190(%%eax)	\n\t"\
		"movl	%[__r7],%%eax		\n\t"\
		"movl	%[__isrt2],%%ebx	\n\t"\
		"movl	%[__cc0],%%ecx		\n\t"\
		"movaps	0x180(%%eax),%%xmm0	\n\t"\
		"movaps	0x190(%%eax),%%xmm1	\n\t"\
		"movaps	0x180(%%eax),%%xmm2	\n\t"\
		"movaps	0x190(%%eax),%%xmm3	\n\t"\
		"mulpd	     (%%ecx),%%xmm0	\n\t"\
		"mulpd	     (%%ecx),%%xmm1	\n\t"\
		"mulpd	0x010(%%ecx),%%xmm2	\n\t"\
		"mulpd	0x010(%%ecx),%%xmm3	\n\t"\
		"subpd	%%xmm2,%%xmm1		\n\t"\
		"addpd	%%xmm3,%%xmm0		\n\t"\
		"movaps	0x100(%%eax),%%xmm4	\n\t"\
		"movaps	0x110(%%eax),%%xmm5	\n\t"\
		"movaps	0x100(%%eax),%%xmm6	\n\t"\
		"movaps	0x110(%%eax),%%xmm7	\n\t"\
		"mulpd	0x010(%%ecx),%%xmm4	\n\t"\
		"mulpd	0x010(%%ecx),%%xmm5	\n\t"\
		"mulpd	     (%%ecx),%%xmm6	\n\t"\
		"mulpd	     (%%ecx),%%xmm7	\n\t"\
		"subpd	%%xmm6,%%xmm5		\n\t"\
		"addpd	%%xmm7,%%xmm4		\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t"\
		"movaps	%%xmm4,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"subpd	%%xmm0,%%xmm6		\n\t"\
		"subpd	%%xmm1,%%xmm7		\n\t"\
		"movaps	0x080(%%eax),%%xmm2	\n\t"\
		"movaps	0x090(%%eax),%%xmm3	\n\t"\
		"movaps	     (%%eax),%%xmm0	\n\t"\
		"movaps	0x010(%%eax),%%xmm1	\n\t"\
		"subpd	0x090(%%eax),%%xmm2	\n\t"\
		"addpd	0x080(%%eax),%%xmm3	\n\t"\
		"mulpd	     (%%ebx),%%xmm2	\n\t"\
		"mulpd	     (%%ebx),%%xmm3	\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"subpd	%%xmm6,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm1		\n\t"\
		"movaps	%%xmm0,0x100(%%eax)	\n\t"\
		"movaps	%%xmm1,0x110(%%eax)	\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm0,%%xmm6		\n\t"\
		"addpd	%%xmm1,%%xmm7		\n\t"\
		"movaps	%%xmm6,     (%%eax)	\n\t"\
		"movaps	%%xmm7,0x010(%%eax)	\n\t"\
		"subpd	%%xmm5,%%xmm2		\n\t"\
		"subpd	%%xmm4,%%xmm3		\n\t"\
		"movaps	%%xmm2,0x180(%%eax)	\n\t"\
		"movaps	%%xmm3,0x090(%%eax)	\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm2,%%xmm5		\n\t"\
		"addpd	%%xmm3,%%xmm4		\n\t"\
		"movaps	%%xmm5,0x080(%%eax)	\n\t"\
		"movaps	%%xmm4,0x190(%%eax)	\n\t"\
	"popl %%ebx	\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__p1] "m" (Xp1)\
		 ,[__p2] "m" (Xp2)\
		 ,[__p3] "m" (Xp3)\
		 ,[__p4] "m" (Xp4)\
		 ,[__r1] "m" (Xr1)\
		 ,[__r3] "m" (Xr3)\
		 ,[__r5] "m" (Xr5)\
		 ,[__r7] "m" (Xr7)\
		 ,[__r9] "m" (Xr9)\
		 ,[__r17] "m" (Xr17)\
		 ,[__r25] "m" (Xr25)\
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__cc0] "m" (Xcc0)\
		: "cc","memory","eax",/*"ebx",*/"ecx","edx","edi","esi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"	/* Clobbered registers */\
	);\
	}


	#define SSE2_RADIX16_DIF_NOTWIDDLE(Xadd0,Xp1,Xp2,Xp3,Xp4,Xr1,Xr3,Xr5,Xr7,Xr9,Xr11,Xr13,Xr15,Xr17,Xr19,Xr21,Xr23,Xr25,Xr27,Xr29,Xr31,Xisrt2,Xcc0)\
	{\
	__asm__ volatile (\
	"pushl %%ebx	\n\t"/* Explicit save/restore of PIC register */\
		"/* SSE2_RADIX4_DIF_IN_PLACE(r1 , r17, r9 , r25): */\n\t"\
		"movl	%[__r1] ,%%eax		\n\t"\
		"movl	%[__r17],%%ebx		\n\t"\
		"movl	%[__r9] ,%%ecx		\n\t"\
		"movl	%[__r25],%%edx		\n\t"\
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"movaps	    (%%eax),%%xmm2	\n\t"\
		"movaps	0x10(%%eax),%%xmm3	\n\t"\
		"addpd	    (%%ebx),%%xmm0	\n\t"\
		"addpd	0x10(%%ebx),%%xmm1	\n\t"\
		"subpd	    (%%ebx),%%xmm2	\n\t"\
		"subpd	0x10(%%ebx),%%xmm3	\n\t"\
		"movaps	    (%%ecx),%%xmm4	\n\t"\
		"movaps	0x10(%%ecx),%%xmm5	\n\t"\
		"movaps	    (%%ecx),%%xmm6	\n\t"\
		"movaps	0x10(%%ecx),%%xmm7	\n\t"\
		"addpd	    (%%edx),%%xmm4	\n\t"\
		"addpd	0x10(%%edx),%%xmm5	\n\t"\
		"subpd	    (%%edx),%%xmm6	\n\t"\
		"subpd	0x10(%%edx),%%xmm7	\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"movaps	%%xmm0,     (%%ebx)	\n\t"\
		"movaps	%%xmm1,0x010(%%ebx)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"movaps	%%xmm4,     (%%eax)	\n\t"\
		"movaps	%%xmm5,0x010(%%eax)	\n\t"\
		"subpd	%%xmm7,%%xmm2		\n\t"\
		"subpd	%%xmm6,%%xmm3		\n\t"\
		"movaps	%%xmm2,     (%%ecx)	\n\t"\
		"movaps	%%xmm3,0x010(%%edx)	\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	%%xmm7,     (%%edx)	\n\t"\
		"movaps	%%xmm6,0x010(%%ecx)	\n\t"\
		"/* SSE2_RADIX4_DIF_IN_PLACE(r5 , r21, r13, r29): */\n\t"\
		"movl	%[__r5] ,%%eax		\n\t"\
		"movl	%[__r21],%%ebx		\n\t"\
		"movl	%[__r13],%%ecx		\n\t"\
		"movl	%[__r29],%%edx		\n\t"\
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"movaps	    (%%eax),%%xmm2	\n\t"\
		"movaps	0x10(%%eax),%%xmm3	\n\t"\
		"addpd	    (%%ebx),%%xmm0	\n\t"\
		"addpd	0x10(%%ebx),%%xmm1	\n\t"\
		"subpd	    (%%ebx),%%xmm2	\n\t"\
		"subpd	0x10(%%ebx),%%xmm3	\n\t"\
		"movaps	    (%%ecx),%%xmm4	\n\t"\
		"movaps	0x10(%%ecx),%%xmm5	\n\t"\
		"movaps	    (%%ecx),%%xmm6	\n\t"\
		"movaps	0x10(%%ecx),%%xmm7	\n\t"\
		"addpd	    (%%edx),%%xmm4	\n\t"\
		"addpd	0x10(%%edx),%%xmm5	\n\t"\
		"subpd	    (%%edx),%%xmm6	\n\t"\
		"subpd	0x10(%%edx),%%xmm7	\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"movaps	%%xmm0,     (%%ebx)	\n\t"\
		"movaps	%%xmm1,0x010(%%ebx)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"movaps	%%xmm4,     (%%eax)	\n\t"\
		"movaps	%%xmm5,0x010(%%eax)	\n\t"\
		"subpd	%%xmm7,%%xmm2		\n\t"\
		"subpd	%%xmm6,%%xmm3		\n\t"\
		"movaps	%%xmm2,     (%%ecx)	\n\t"\
		"movaps	%%xmm3,0x010(%%edx)	\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	%%xmm7,     (%%edx)	\n\t"\
		"movaps	%%xmm6,0x010(%%ecx)	\n\t"\
		"/* SSE2_RADIX4_DIF_IN_PLACE(r3 , r19, r11, r27): */\n\t"\
		"movl	%[__r3] ,%%eax		\n\t"\
		"movl	%[__r19],%%ebx		\n\t"\
		"movl	%[__r11],%%ecx		\n\t"\
		"movl	%[__r27],%%edx		\n\t"\
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"movaps	    (%%eax),%%xmm2	\n\t"\
		"movaps	0x10(%%eax),%%xmm3	\n\t"\
		"addpd	    (%%ebx),%%xmm0	\n\t"\
		"addpd	0x10(%%ebx),%%xmm1	\n\t"\
		"subpd	    (%%ebx),%%xmm2	\n\t"\
		"subpd	0x10(%%ebx),%%xmm3	\n\t"\
		"movaps	    (%%ecx),%%xmm4	\n\t"\
		"movaps	0x10(%%ecx),%%xmm5	\n\t"\
		"movaps	    (%%ecx),%%xmm6	\n\t"\
		"movaps	0x10(%%ecx),%%xmm7	\n\t"\
		"addpd	    (%%edx),%%xmm4	\n\t"\
		"addpd	0x10(%%edx),%%xmm5	\n\t"\
		"subpd	    (%%edx),%%xmm6	\n\t"\
		"subpd	0x10(%%edx),%%xmm7	\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"movaps	%%xmm0,     (%%ebx)	\n\t"\
		"movaps	%%xmm1,0x010(%%ebx)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"movaps	%%xmm4,     (%%eax)	\n\t"\
		"movaps	%%xmm5,0x010(%%eax)	\n\t"\
		"subpd	%%xmm7,%%xmm2		\n\t"\
		"subpd	%%xmm6,%%xmm3		\n\t"\
		"movaps	%%xmm2,     (%%ecx)	\n\t"\
		"movaps	%%xmm3,0x010(%%edx)	\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	%%xmm7,     (%%edx)	\n\t"\
		"movaps	%%xmm6,0x010(%%ecx)	\n\t"\
		"/* SSE2_RADIX4_DIF_IN_PLACE(r7 , r23, r15, r31): */\n\t"\
		"movl	%[__r7] ,%%eax		\n\t"\
		"movl	%[__r23],%%ebx		\n\t"\
		"movl	%[__r15],%%ecx		\n\t"\
		"movl	%[__r31],%%edx		\n\t"\
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"movaps	    (%%eax),%%xmm2	\n\t"\
		"movaps	0x10(%%eax),%%xmm3	\n\t"\
		"addpd	    (%%ebx),%%xmm0	\n\t"\
		"addpd	0x10(%%ebx),%%xmm1	\n\t"\
		"subpd	    (%%ebx),%%xmm2	\n\t"\
		"subpd	0x10(%%ebx),%%xmm3	\n\t"\
		"movaps	    (%%ecx),%%xmm4	\n\t"\
		"movaps	0x10(%%ecx),%%xmm5	\n\t"\
		"movaps	    (%%ecx),%%xmm6	\n\t"\
		"movaps	0x10(%%ecx),%%xmm7	\n\t"\
		"addpd	    (%%edx),%%xmm4	\n\t"\
		"addpd	0x10(%%edx),%%xmm5	\n\t"\
		"subpd	    (%%edx),%%xmm6	\n\t"\
		"subpd	0x10(%%edx),%%xmm7	\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"movaps	%%xmm0,     (%%ebx)	\n\t"\
		"movaps	%%xmm1,0x010(%%ebx)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"movaps	%%xmm4,     (%%eax)	\n\t"\
		"movaps	%%xmm5,0x010(%%eax)	\n\t"\
		"subpd	%%xmm7,%%xmm2		\n\t"\
		"subpd	%%xmm6,%%xmm3		\n\t"\
		"movaps	%%xmm2,     (%%ecx)	\n\t"\
		"movaps	%%xmm3,0x010(%%edx)	\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	%%xmm7,     (%%edx)	\n\t"\
		"movaps	%%xmm6,0x010(%%ecx)	\n\t"\
		"movl	%[__add0],%%eax		\n\t"\
		"movl	%[__p1],%%ebx		\n\t"\
		"movl	%[__p2],%%ecx		\n\t"\
		"movl	%[__p3],%%edx		\n\t"\
		"movl	%[__p4],%%edi		\n\t"\
		"shll	$3,%%ebx			\n\t"\
		"shll	$3,%%ecx			\n\t"\
		"shll	$3,%%edx			\n\t"\
		"shll	$3,%%edi			\n\t"\
		"addl	%%eax,%%ebx			\n\t"\
		"addl	%%eax,%%ecx			\n\t"\
		"addl	%%eax,%%edx			\n\t"\
		"movl	%[__r1],%%esi		\n\t"\
		"movaps	     (%%esi),%%xmm0	\n\t"\
		"movaps	0x010(%%esi),%%xmm1	\n\t"\
		"movaps	0x040(%%esi),%%xmm2	\n\t"\
		"movaps	0x050(%%esi),%%xmm3	\n\t"\
		"subpd	0x040(%%esi),%%xmm0	\n\t"\
		"subpd	0x050(%%esi),%%xmm1	\n\t"\
		"addpd	     (%%esi),%%xmm2	\n\t"\
		"addpd	0x010(%%esi),%%xmm3	\n\t"\
		"movaps	0x020(%%esi),%%xmm4	\n\t"\
		"movaps	0x030(%%esi),%%xmm5	\n\t"\
		"movaps	0x060(%%esi),%%xmm6	\n\t"\
		"movaps	0x070(%%esi),%%xmm7	\n\t"\
		"subpd	0x060(%%esi),%%xmm4	\n\t"\
		"subpd	0x070(%%esi),%%xmm5	\n\t"\
		"addpd	0x020(%%esi),%%xmm6	\n\t"\
		"addpd	0x030(%%esi),%%xmm7	\n\t"\
		"subpd	%%xmm6,%%xmm2		\n\t"\
		"subpd	%%xmm7,%%xmm3		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"movaps	%%xmm2,    (%%ebx)	\n\t"\
		"movaps	%%xmm3,0x10(%%ebx)	\n\t"\
		"addpd	%%xmm2,		%%xmm6	\n\t"\
		"addpd	%%xmm3,		%%xmm7	\n\t"\
		"movaps	%%xmm6,    (%%eax)	\n\t"\
		"movaps	%%xmm7,0x10(%%eax)	\n\t"\
		"subpd	%%xmm5,		%%xmm0	\n\t"\
		"subpd	%%xmm4,		%%xmm1	\n\t"\
		"addpd	%%xmm5,		%%xmm5	\n\t"\
		"addpd	%%xmm4,		%%xmm4	\n\t"\
		"movaps	%%xmm0,    (%%ecx)	\n\t"\
		"movaps	%%xmm1,0x10(%%edx)	\n\t"\
		"addpd	%%xmm0,		%%xmm5	\n\t"\
		"addpd	%%xmm1,		%%xmm4	\n\t"\
		"movaps	%%xmm5,    (%%edx)	\n\t"\
		"movaps	%%xmm4,0x10(%%ecx)	\n\t"\
		"addl	%%edi,%%eax			\n\t"\
		"addl	%%edi,%%ebx			\n\t"\
		"addl	%%edi,%%ecx			\n\t"\
		"addl	%%edi,%%edx			\n\t"\
		"movl	%[__r17],%%esi		\n\t"\
		"movaps	     (%%esi),%%xmm0	\n\t"\
		"movaps	0x010(%%esi),%%xmm1	\n\t"\
		"movaps	0x040(%%esi),%%xmm2	\n\t"\
		"movaps	0x050(%%esi),%%xmm3	\n\t"\
		"subpd	0x050(%%esi),%%xmm0	\n\t"\
		"subpd	0x040(%%esi),%%xmm1	\n\t"\
		"addpd	0x010(%%esi),%%xmm2	\n\t"\
		"addpd	     (%%esi),%%xmm3	\n\t"\
		"movaps	0x020(%%esi),%%xmm4	\n\t"\
		"movaps	0x030(%%esi),%%xmm5	\n\t"\
		"movaps	0x060(%%esi),%%xmm6	\n\t"\
		"movaps	0x070(%%esi),%%xmm7	\n\t"\
		"subpd	0x030(%%esi),%%xmm4	\n\t"\
		"addpd	0x020(%%esi),%%xmm5	\n\t"\
		"addpd	0x070(%%esi),%%xmm6	\n\t"\
		"subpd	0x060(%%esi),%%xmm7	\n\t"\
		"movl	%[__isrt2],%%esi	\n\t"\
		"mulpd	(%%esi),%%xmm4		\n\t"\
		"mulpd	(%%esi),%%xmm5		\n\t"\
		"mulpd	(%%esi),%%xmm6		\n\t"\
		"mulpd	(%%esi),%%xmm7		\n\t"\
		"subpd	%%xmm6,%%xmm4		\n\t"\
		"subpd	%%xmm7,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm6		\n\t"\
		"addpd	%%xmm5,%%xmm7		\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm2		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"movaps	%%xmm0,    (%%ebx)	\n\t"\
		"movaps	%%xmm2,0x10(%%ebx)	\n\t"\
		"addpd	%%xmm0,		%%xmm4	\n\t"\
		"addpd	%%xmm2,		%%xmm5	\n\t"\
		"movaps	%%xmm4,    (%%eax)	\n\t"\
		"movaps	%%xmm5,0x10(%%eax)	\n\t"\
		"subpd	%%xmm7,		%%xmm3	\n\t"\
		"subpd	%%xmm6,		%%xmm1	\n\t"\
		"addpd	%%xmm7,		%%xmm7	\n\t"\
		"addpd	%%xmm6,		%%xmm6	\n\t"\
		"movaps	%%xmm3,    (%%ecx)	\n\t"\
		"movaps	%%xmm1,0x10(%%edx)	\n\t"\
		"addpd	%%xmm3,		%%xmm7	\n\t"\
		"addpd	%%xmm1,		%%xmm6	\n\t"\
		"movaps	%%xmm7,    (%%edx)	\n\t"\
		"movaps	%%xmm6,0x10(%%ecx)	\n\t"\
		"addl	%%edi,%%eax			\n\t"\
		"addl	%%edi,%%ebx			\n\t"\
		"addl	%%edi,%%ecx			\n\t"\
		"addl	%%edi,%%edx			\n\t"\
		"movl	%[__r9] ,%%esi		\n\t"\
		"movl	%[__cc0],%%edi		\n\t"\
		"movaps	0x020(%%esi),%%xmm4	\n\t"\
		"movaps	0x060(%%esi),%%xmm6	\n\t"\
		"movaps	0x030(%%esi),%%xmm5	\n\t"\
		"movaps	0x070(%%esi),%%xmm7	\n\t"\
		"movaps	%%xmm4,%%xmm0		\n\t"\
		"movaps	%%xmm6,%%xmm2		\n\t"\
		"movaps	%%xmm5,%%xmm1		\n\t"\
		"movaps	%%xmm7,%%xmm3		\n\t"\
		"mulpd	     (%%edi),%%xmm4	\n\t"\
		"mulpd	0x010(%%edi),%%xmm6	\n\t"\
		"mulpd	0x010(%%edi),%%xmm1	\n\t"\
		"mulpd	     (%%edi),%%xmm3	\n\t"\
		"mulpd	     (%%edi),%%xmm5	\n\t"\
		"mulpd	0x010(%%edi),%%xmm7	\n\t"\
		"mulpd	0x010(%%edi),%%xmm0	\n\t"\
		"mulpd	     (%%edi),%%xmm2	\n\t"\
		"subpd	%%xmm1,%%xmm4		\n\t"\
		"subpd	%%xmm3,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm5		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"subpd	%%xmm6,%%xmm4		\n\t"\
		"subpd	%%xmm7,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm6		\n\t"\
		"addpd	%%xmm5,%%xmm7		\n\t"\
		"movl	%[__isrt2],%%edi	\n\t"\
		"movaps	0x040(%%esi),%%xmm2	\n\t"\
		"movaps	0x050(%%esi),%%xmm3	\n\t"\
		"subpd	0x050(%%esi),%%xmm2	\n\t"\
		"addpd	0x040(%%esi),%%xmm3	\n\t"\
		"mulpd	(%%edi),%%xmm2		\n\t"\
		"mulpd	(%%edi),%%xmm3		\n\t"\
		"movaps	     (%%esi),%%xmm0	\n\t"\
		"movaps	0x010(%%esi),%%xmm1	\n\t"\
		"movl	%[__p4],%%edi		\n\t"\
		"shll	$3,%%edi			\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"addpd	     (%%esi),%%xmm2	\n\t"\
		"addpd	0x010(%%esi),%%xmm3	\n\t"\
		"subpd	%%xmm6,%%xmm2		\n\t"\
		"subpd	%%xmm7,%%xmm3		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"movaps	%%xmm2,    (%%ebx)	\n\t"\
		"movaps	%%xmm3,0x10(%%ebx)	\n\t"\
		"addpd	%%xmm2,		%%xmm6	\n\t"\
		"addpd	%%xmm3,		%%xmm7	\n\t"\
		"movaps	%%xmm6,    (%%eax)	\n\t"\
		"movaps	%%xmm7,0x10(%%eax)	\n\t"\
		"subpd	%%xmm5,		%%xmm0	\n\t"\
		"subpd	%%xmm4,		%%xmm1	\n\t"\
		"addpd	%%xmm5,		%%xmm5	\n\t"\
		"addpd	%%xmm4,		%%xmm4	\n\t"\
		"movaps	%%xmm0,    (%%ecx)	\n\t"\
		"movaps	%%xmm1,0x10(%%edx)	\n\t"\
		"addpd	%%xmm0,		%%xmm5	\n\t"\
		"addpd	%%xmm1,		%%xmm4	\n\t"\
		"movaps	%%xmm5,    (%%edx)	\n\t"\
		"movaps	%%xmm4,0x10(%%ecx)	\n\t"\
		"addl	%%edi,%%eax			\n\t"\
		"addl	%%edi,%%ebx			\n\t"\
		"addl	%%edi,%%ecx			\n\t"\
		"addl	%%edi,%%edx			\n\t"\
		"movl	%[__r25],%%esi		\n\t"\
		"movl	%[__cc0],%%edi		\n\t"\
		"movaps	0x020(%%esi),%%xmm4	\n\t"\
		"movaps	0x060(%%esi),%%xmm6	\n\t"\
		"movaps	0x030(%%esi),%%xmm5	\n\t"\
		"movaps	0x070(%%esi),%%xmm7	\n\t"\
		"movaps	%%xmm4,%%xmm0		\n\t"\
		"movaps	%%xmm6,%%xmm2		\n\t"\
		"movaps	%%xmm5,%%xmm1		\n\t"\
		"movaps	%%xmm7,%%xmm3		\n\t"\
		"mulpd	0x010(%%edi),%%xmm4	\n\t"\
		"mulpd	     (%%edi),%%xmm6	\n\t"\
		"mulpd	     (%%edi),%%xmm1	\n\t"\
		"mulpd	0x010(%%edi),%%xmm3	\n\t"\
		"mulpd	0x010(%%edi),%%xmm5	\n\t"\
		"mulpd	     (%%edi),%%xmm7	\n\t"\
		"mulpd	     (%%edi),%%xmm0	\n\t"\
		"mulpd	0x010(%%edi),%%xmm2	\n\t"\
		"subpd	%%xmm1,%%xmm4		\n\t"\
		"subpd	%%xmm3,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm5		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"subpd	%%xmm6,%%xmm4		\n\t"\
		"subpd	%%xmm7,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm6		\n\t"\
		"addpd	%%xmm5,%%xmm7		\n\t"\
		"movl	%[__isrt2],%%edi	\n\t"\
		"movaps	0x040(%%esi),%%xmm2	\n\t"\
		"movaps	0x050(%%esi),%%xmm3	\n\t"\
		"addpd	0x050(%%esi),%%xmm2	\n\t"\
		"subpd	0x040(%%esi),%%xmm3	\n\t"\
		"mulpd	(%%edi),%%xmm2		\n\t"\
		"mulpd	(%%edi),%%xmm3		\n\t"\
		"movaps	     (%%esi),%%xmm0	\n\t"\
		"movaps	0x010(%%esi),%%xmm1	\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"addpd	     (%%esi),%%xmm2	\n\t"\
		"addpd	0x010(%%esi),%%xmm3	\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"movaps	%%xmm0,    (%%ebx)	\n\t"\
		"movaps	%%xmm1,0x10(%%ebx)	\n\t"\
		"addpd	%%xmm0,		%%xmm4	\n\t"\
		"addpd	%%xmm1,		%%xmm5	\n\t"\
		"movaps	%%xmm4,    (%%eax)	\n\t"\
		"movaps	%%xmm5,0x10(%%eax)	\n\t"\
		"subpd	%%xmm7,		%%xmm2	\n\t"\
		"subpd	%%xmm6,		%%xmm3	\n\t"\
		"addpd	%%xmm7,		%%xmm7	\n\t"\
		"addpd	%%xmm6,		%%xmm6	\n\t"\
		"movaps	%%xmm2,    (%%ecx)	\n\t"\
		"movaps	%%xmm3,0x10(%%edx)	\n\t"\
		"addpd	%%xmm2,		%%xmm7	\n\t"\
		"addpd	%%xmm3,		%%xmm6	\n\t"\
		"movaps	%%xmm7,    (%%edx)	\n\t"\
		"movaps	%%xmm6,0x10(%%ecx)	\n\t"\
	"popl %%ebx	\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__p1] "m" (Xp1)\
		 ,[__p2] "m" (Xp2)\
		 ,[__p3] "m" (Xp3)\
		 ,[__p4] "m" (Xp4)\
		 ,[__r1] "m" (Xr1)\
		 ,[__r3] "m" (Xr3)\
		 ,[__r5] "m" (Xr5)\
		 ,[__r7] "m" (Xr7)\
		 ,[__r9] "m" (Xr9)\
		 ,[__r11] "m" (Xr11)\
		 ,[__r13] "m" (Xr13)\
		 ,[__r15] "m" (Xr15)\
		 ,[__r17] "m" (Xr17)\
		 ,[__r19] "m" (Xr19)\
		 ,[__r21] "m" (Xr21)\
		 ,[__r23] "m" (Xr23)\
		 ,[__r25] "m" (Xr25)\
		 ,[__r27] "m" (Xr27)\
		 ,[__r29] "m" (Xr29)\
		 ,[__r31] "m" (Xr31)\
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__cc0] "m" (Xcc0)\
		: "cc","memory","eax",/*"ebx",*/"ecx","edx","edi","esi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"	/* Clobbered registers */\
	);\
	}

#else

	#error Unhandled combination of #defs!

#endif	// SSE2 or AVX?

#endif	/* radix16_ditN_cy_dif1_gcc_h_included */

