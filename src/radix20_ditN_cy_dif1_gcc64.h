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
#ifndef radix20_ditN_cy_dif1_gcc_h_included
#define radix20_ditN_cy_dif1_gcc_h_included

  #ifdef USE_AVX

	#define	SSE2_RADIX20_DIT_NOTWIDDLE(Xadd0,Xp01,Xp04,Xr00,Xr10,Xr20,Xr30,Xcc1,Xa0,Xa1,Xa2,Xa3,Xa4,Xb0,Xb1,Xb2,Xb3,Xb4,Xc0,Xc1,Xc2,Xc3,Xc4,Xd0,Xd1,Xd2,Xd3,Xd4)\
	{\
	__asm__ volatile (\
		"/* Outputs in SSE2 modes are temps 2*5*16 = 10*16 = 0x140 bytes apart: */\n\t"\
		"/*	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_C(abcd,r00)	*/		/*	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_C(dcba,r02)	*/\n\t"\
		"movq	%[__add0],%%rax	/* Use eax as base address */	\n\t	movq	$0x140,%%rdi		\n\t"\
		"/*	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_C(abcd,r00)	*/		\n\t	movslq	%[__p04],%%r8 		\n\t"\
		"movq	%%rax	,%%rdx									\n\t	shlq	$3,%%r8 			\n\t"\
		"movslq	%[__p01],%%rsi									\n\t	movq	$0xa00,%%r9 		\n\t"\
		"shlq	$3,%%rsi										\n\t	movq	%%r8 ,%%r10			\n\t"\
		"addq	%%rsi,%%rdx										\n\t	movq	%%r8 ,%%r11			\n\t"\
		"movq	%%rdx,%%rbx		/* add1 = add0+p01 */			\n\t	movq	%%r8 ,%%r12			\n\t"\
		"addq	%%rsi,%%rdx										\n\t	movq	%%r8 ,%%r13			\n\t"\
		"movq	%%rdx,%%rcx		/* add3 = add0+p02 */			\n\t	addq	%%rax,%%r10		/* &a[j1+p04] */\n\t"\
		"addq	%%rsi,%%rdx		/* add2 = add0+p03 */			\n\t	addq	%%rbx,%%r11			\n\t"\
		"vmovaps	0x20(%%rax),%%ymm3							\n\t	movq	%[__r00],%%rsi		\n\t"\
		"vmovaps	    (%%rax),%%ymm2							\n\t	addq	%%rcx,%%r12			\n\t"\
		"vmovaps	    (%%rbx),%%ymm0							\n\t	addq	%%rdx,%%r13			\n\t"\
		"vmovaps	0x20(%%rbx),%%ymm1							\n\t	addq	%%rsi,%%r9	/* two */	\n\t"\
		"vmovaps	    (%%rcx),%%ymm6							\n\t	vmovaps	0x20(%%r13),%%ymm11	\n\t"\
		"vmovaps	0x20(%%rcx),%%ymm7							\n\t	vmovaps	    (%%r13),%%ymm10	\n\t"\
		"vmovaps	0x20(%%rdx),%%ymm5							\n\t	vmovaps	    (%%r12),%%ymm8 	\n\t"\
		"vmovaps	    (%%rdx),%%ymm4							\n\t	vmovaps	0x20(%%r12),%%ymm9 	\n\t"\
		"vsubpd	%%ymm0,%%ymm2,%%ymm2							\n\t	vmovaps	    (%%r11),%%ymm14	\n\t"\
		"vsubpd	%%ymm4,%%ymm6,%%ymm6							\n\t	vmovaps	0x20(%%r11),%%ymm15	\n\t"\
		"vsubpd	%%ymm1,%%ymm3,%%ymm3							\n\t	vmovaps	0x20(%%r10),%%ymm13	\n\t"\
		"vsubpd	%%ymm5,%%ymm7,%%ymm7							\n\t	vmovaps	    (%%r10),%%ymm12	\n\t"\
		"vmulpd	(%%r9),%%ymm0,%%ymm0							\n\t	vsubpd	%%ymm8 ,%%ymm10,%%ymm10		\n\t"\
		"vmulpd	(%%r9),%%ymm4,%%ymm4							\n\t	vsubpd	%%ymm12,%%ymm14,%%ymm14		\n\t"\
		"vmulpd	(%%r9),%%ymm1,%%ymm1							\n\t	vsubpd	%%ymm9 ,%%ymm11,%%ymm11		\n\t"\
		"vmulpd	(%%r9),%%ymm5,%%ymm5							\n\t	vsubpd	%%ymm13,%%ymm15,%%ymm15		\n\t"\
		"vaddpd	%%ymm2,%%ymm0,%%ymm0							\n\t	vmulpd	(%%r9) ,%%ymm8 ,%%ymm8 		\n\t"\
		"vaddpd	%%ymm6,%%ymm4,%%ymm4							\n\t	vmulpd	(%%r9) ,%%ymm12,%%ymm12		\n\t"\
		"vaddpd	%%ymm3,%%ymm1,%%ymm1							\n\t	vmulpd	(%%r9) ,%%ymm9 ,%%ymm9 		\n\t"\
		"vaddpd	%%ymm7,%%ymm5,%%ymm5							\n\t	vmulpd	(%%r9) ,%%ymm13,%%ymm13		\n\t"\
		"/* Finish radix-4 butterfly: */\n\t	addq	%%rsi,%%rdi	\n\t"\
		"vsubpd	%%ymm4,%%ymm0,%%ymm0							\n\t	vaddpd	%%ymm10,%%ymm8 ,%%ymm8 		\n\t"\
		"vsubpd	%%ymm7,%%ymm2,%%ymm2							\n\t	vaddpd	%%ymm14,%%ymm12,%%ymm12		\n\t"\
		"vsubpd	%%ymm5,%%ymm1,%%ymm1							\n\t	vaddpd	%%ymm11,%%ymm9 ,%%ymm9 		\n\t"\
		"vsubpd	%%ymm6,%%ymm3,%%ymm3							\n\t	vaddpd	%%ymm15,%%ymm13,%%ymm13		\n\t"\
		"vmovaps	%%ymm0,0x280(%%rsi)							\n\t	vsubpd	%%ymm12,%%ymm8 ,%%ymm8 		\n\t"\
		"vmovaps	%%ymm2,     (%%rdi)							\n\t	vsubpd	%%ymm15,%%ymm10,%%ymm10		\n\t"\
		"vmovaps	%%ymm1,0x2a0(%%rsi)							\n\t	vsubpd	%%ymm13,%%ymm9 ,%%ymm9 		\n\t"\
		"vmovaps	%%ymm3,0x2a0(%%rdi)							\n\t	vsubpd	%%ymm14,%%ymm11,%%ymm11		\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4							\n\t	vmovaps	%%ymm8 ,0x2c0(%%rsi)	\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7							\n\t	vmovaps	%%ymm10,0x2c0(%%rdi)	\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5							\n\t	vmovaps	%%ymm9 ,0x2e0(%%rsi)	\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6							\n\t	vmovaps	%%ymm11,0x060(%%rdi)	\n\t"\
		"vaddpd	%%ymm0,%%ymm4,%%ymm4							\n\t	vmulpd	(%%r9) ,%%ymm12,%%ymm12		\n\t"\
		"vaddpd	%%ymm2,%%ymm7,%%ymm7							\n\t	vmulpd	(%%r9) ,%%ymm15,%%ymm15		\n\t"\
		"vaddpd	%%ymm1,%%ymm5,%%ymm5							\n\t	vmulpd	(%%r9) ,%%ymm13,%%ymm13		\n\t"\
		"vaddpd	%%ymm3,%%ymm6,%%ymm6							\n\t	vmulpd	(%%r9) ,%%ymm14,%%ymm14		\n\t"\
		"vmovaps	%%ymm4,     (%%rsi)							\n\t	vaddpd	%%ymm8 ,%%ymm12,%%ymm12		\n\t"\
		"vmovaps	%%ymm7,0x280(%%rdi)							\n\t	vaddpd	%%ymm10,%%ymm15,%%ymm15		\n\t"\
		"vmovaps	%%ymm5,0x020(%%rsi)							\n\t	vaddpd	%%ymm9 ,%%ymm13,%%ymm13		\n\t"\
		"vmovaps	%%ymm6,0x020(%%rdi)							\n\t	vaddpd	%%ymm11,%%ymm14,%%ymm14		\n\t"\
		"addq	%%r8 ,%%r8 		/* p08 */	\n\t"\
		"addq	%%r8 ,%%rax		/* &a[j1+p08] */				\n\t	vmovaps	%%ymm12,0x040(%%rsi)	\n\t"\
		"addq	%%r8 ,%%rbx										\n\t	vmovaps	%%ymm15,0x040(%%rdi)	\n\t"\
		"addq	%%r8 ,%%rcx										\n\t	vmovaps	%%ymm13,0x060(%%rsi)	\n\t"\
		"addq	%%r8 ,%%rdx										\n\t	vmovaps	%%ymm14,0x2e0(%%rdi)	\n\t"\
		"/*	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_C(bacd,r04)	*/		\n\t	/*	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_C(cdab,r06)	*/\n\t"\
		"vmovaps	0x20(%%rbx),%%ymm3							\n\t	addq	%%r8 ,%%r10		/* &a[j1+p12] */\n\t"\
		"vmovaps	    (%%rbx),%%ymm2							\n\t	addq	%%r8 ,%%r11			\n\t"\
		"vmovaps	    (%%rax),%%ymm0							\n\t	addq	%%r8 ,%%r12			\n\t"\
		"vmovaps	0x20(%%rax),%%ymm1							\n\t	addq	%%r8 ,%%r13			\n\t"\
		"vmovaps	    (%%rcx),%%ymm6							\n\t	vmovaps	0x20(%%r12),%%ymm11	\n\t"\
		"vmovaps	0x20(%%rcx),%%ymm7							\n\t	vmovaps	    (%%r12),%%ymm10	\n\t"\
		"vmovaps	0x20(%%rdx),%%ymm5							\n\t	vmovaps	    (%%r13),%%ymm8 	\n\t"\
		"vmovaps	    (%%rdx),%%ymm4							\n\t	vmovaps	0x20(%%r13),%%ymm9 	\n\t"\
		"vsubpd	%%ymm0,%%ymm2,%%ymm2							\n\t	vmovaps	    (%%r10),%%ymm14	\n\t"\
		"vsubpd	%%ymm4,%%ymm6,%%ymm6							\n\t	vmovaps	0x20(%%r10),%%ymm15	\n\t"\
		"vsubpd	%%ymm1,%%ymm3,%%ymm3							\n\t	vmovaps	0x20(%%r11),%%ymm13	\n\t"\
		"vsubpd	%%ymm5,%%ymm7,%%ymm7							\n\t	vmovaps	    (%%r11),%%ymm12	\n\t"\
		"vmulpd	(%%r9),%%ymm0,%%ymm0							\n\t	vsubpd	%%ymm8 ,%%ymm10,%%ymm10		\n\t"\
		"vmulpd	(%%r9),%%ymm4,%%ymm4							\n\t	vsubpd	%%ymm12,%%ymm14,%%ymm14		\n\t"\
		"vmulpd	(%%r9),%%ymm1,%%ymm1							\n\t	vsubpd	%%ymm9 ,%%ymm11,%%ymm11		\n\t"\
		"vmulpd	(%%r9),%%ymm5,%%ymm5							\n\t	vsubpd	%%ymm13,%%ymm15,%%ymm15		\n\t"\
		"vaddpd	%%ymm2,%%ymm0,%%ymm0							\n\t	vmulpd	(%%r9) ,%%ymm8 ,%%ymm8 		\n\t"\
		"vaddpd	%%ymm6,%%ymm4,%%ymm4							\n\t	vmulpd	(%%r9) ,%%ymm12,%%ymm12		\n\t"\
		"vaddpd	%%ymm3,%%ymm1,%%ymm1							\n\t	vmulpd	(%%r9) ,%%ymm9 ,%%ymm9 		\n\t"\
		"vaddpd	%%ymm7,%%ymm5,%%ymm5							\n\t	vmulpd	(%%r9) ,%%ymm13,%%ymm13		\n\t"\
		"/* Finish radix-4 butterfly: */						\n\t	vaddpd	%%ymm10,%%ymm8 ,%%ymm8 		\n\t"\
		"addq	$0x080,%%rsi	/* r04 */						\n\t	vaddpd	%%ymm14,%%ymm12,%%ymm12		\n\t"\
		"vsubpd	%%ymm4,%%ymm0,%%ymm0							\n\t	vaddpd	%%ymm11,%%ymm9 ,%%ymm9 		\n\t"\
		"vsubpd	%%ymm7,%%ymm2,%%ymm2							\n\t	vaddpd	%%ymm15,%%ymm13,%%ymm13		\n\t"\
		"vsubpd	%%ymm5,%%ymm1,%%ymm1							\n\t	/* Finish radix-4 butterfly: */\n\t"\
		"vsubpd	%%ymm6,%%ymm3,%%ymm3							\n\t	addq	$0x080,%%rdi		\n\t"\
		"vmovaps	%%ymm0,0x280(%%rsi)							\n\t	vsubpd	%%ymm12,%%ymm8 ,%%ymm8 		\n\t"\
		"vmovaps	%%ymm2,0x280(%%rdi)							\n\t	vsubpd	%%ymm15,%%ymm10,%%ymm10		\n\t"\
		"vmovaps	%%ymm1,0x2a0(%%rsi)							\n\t	vsubpd	%%ymm13,%%ymm9 ,%%ymm9 		\n\t"\
		"vmovaps	%%ymm3,0x020(%%rdi)							\n\t	vsubpd	%%ymm14,%%ymm11,%%ymm11		\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4							\n\t	vmovaps	%%ymm8 ,0x2c0(%%rsi)	\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7							\n\t	vmovaps	%%ymm10,0x2c0(%%rdi)	\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5							\n\t	vmovaps	%%ymm9 ,0x2e0(%%rsi)	\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6							\n\t	vmovaps	%%ymm11,0x060(%%rdi)	\n\t"\
		"vaddpd	%%ymm0,%%ymm4,%%ymm4							\n\t	vmulpd	(%%r9) ,%%ymm12,%%ymm12		\n\t"\
		"vaddpd	%%ymm2,%%ymm7,%%ymm7							\n\t	vmulpd	(%%r9) ,%%ymm15,%%ymm15		\n\t"\
		"vaddpd	%%ymm1,%%ymm5,%%ymm5							\n\t	vmulpd	(%%r9) ,%%ymm13,%%ymm13		\n\t"\
		"vaddpd	%%ymm3,%%ymm6,%%ymm6							\n\t	vmulpd	(%%r9) ,%%ymm14,%%ymm14		\n\t"\
		"vmovaps	%%ymm4,     (%%rsi)							\n\t	vaddpd	%%ymm8 ,%%ymm12,%%ymm12		\n\t"\
		"vmovaps	%%ymm7,     (%%rdi)							\n\t	vaddpd	%%ymm10,%%ymm15,%%ymm15		\n\t"\
		"vmovaps	%%ymm5,0x020(%%rsi)							\n\t	vaddpd	%%ymm9 ,%%ymm13,%%ymm13		\n\t"\
		"vmovaps	%%ymm6,0x2a0(%%rdi)							\n\t	vaddpd	%%ymm11,%%ymm14,%%ymm14		\n\t"\
		"														\n\t	vmovaps	%%ymm12,0x040(%%rsi)	\n\t"\
		"														\n\t	vmovaps	%%ymm15,0x040(%%rdi)	\n\t"\
		"														\n\t	vmovaps	%%ymm13,0x060(%%rsi)	\n\t"\
		"														\n\t	vmovaps	%%ymm14,0x2e0(%%rdi)	\n\t"\
		"/*	add0,1,2,3 = &a[j1+p16]+p0,1,3,2 */\n\t"\
		"addq	%%r8 ,%%rax		/* &a[j1+p16] */\n\t"\
		"addq	%%r8 ,%%rbx			\n\t"\
		"addq	%%r8 ,%%rcx			\n\t"\
		"addq	%%r8 ,%%rdx			\n\t"\
		"/*	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_C(abdc,r08)	*/\n\t"\
		"vmovaps	0x20(%%rax),%%ymm3	\n\t"\
		"vmovaps	    (%%rax),%%ymm2	\n\t"\
		"vmovaps	    (%%rbx),%%ymm0	\n\t"\
		"vmovaps	0x20(%%rbx),%%ymm1	\n\t"\
		"vmovaps	    (%%rdx),%%ymm6	\n\t"\
		"vmovaps	0x20(%%rdx),%%ymm7	\n\t"\
		"vmovaps	0x20(%%rcx),%%ymm5	\n\t"\
		"vmovaps	    (%%rcx),%%ymm4	\n\t"\
		"addq	$0x080,%%rsi	/* r08 */\n\t"\
		"vsubpd	%%ymm0,%%ymm2,%%ymm2		\n\t"\
		"vsubpd	%%ymm4,%%ymm6,%%ymm6		\n\t"\
		"vsubpd	%%ymm1,%%ymm3,%%ymm3		\n\t"\
		"vsubpd	%%ymm5,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	%%ymm0,%%ymm0,%%ymm0		\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	%%ymm1,%%ymm1,%%ymm1		\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5							\n\t	movq	%%r9,%%r14	/* two */\n\t"\
		"vaddpd	%%ymm2,%%ymm0,%%ymm0		\n\t"\
		"vaddpd	%%ymm6,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	%%ymm3,%%ymm1,%%ymm1							\n\t/*** This needs to read r18 *after* the last radix-4 writes it! ***/\n\t"\
		"vaddpd	%%ymm7,%%ymm5,%%ymm5							\n\t/* Radix-5 DFT uses adjacent temps, i.e. stride = 2*16 bytes: */\n\t"\
		"/* Finish radix-4 butterfly: */						\n\t/*	SSE2_RADIX_05_DFT_0TWIDDLE(r10,r12,r14,r16,r18,cc1,Xb0,Xb1,Xb2,Xb3,Xb4)	*/\n\t"\
		"addq	$0x080,%%rdi	/* r18 */						\n\t	movq	%[__r10],%%r8 		\n\t"\
		"vsubpd	%%ymm4,%%ymm0,%%ymm0							\n\t	movq	%%r8 ,%%r10			\n\t"\
		"vsubpd	%%ymm7,%%ymm2,%%ymm2							\n\t	movq	%%r8 ,%%r11			\n\t"\
		"vsubpd	%%ymm5,%%ymm1,%%ymm1							\n\t	movq	%%r8 ,%%r12			\n\t"\
		"vsubpd	%%ymm6,%%ymm3,%%ymm3							\n\t/*	movq	%%r8 ,%%r13		*/	\n\t"\
		"vmovaps	%%ymm0,0x280(%%rsi)	/* r28.r */				\n\t	addq	$0x040,%%r10		\n\t"\
		"vmovaps	%%ymm2,0x280(%%rdi)	/* r38.r */				\n\t	addq	$0x080,%%r11		\n\t"\
		"vmovaps	%%ymm1,0x2a0(%%rsi)	/* r28.i */				\n\t	addq	$0x0c0,%%r12		\n\t"\
		"/*vaps	%%ymm3,0x020(%%rdi)	** r18.i */					\n\t/*	addq	$0x100,%%r13	*/	\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4							\n\t	movq	%[__b0],%%r9 		\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7							\n\t	vmovaps	    (%%r10),%%ymm8 	\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5							\n\t	vmovaps	0x20(%%r10),%%ymm9 	\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6							\n\t	vmovaps	    (%%r11),%%ymm10	\n\t"\
		"vaddpd	%%ymm0,%%ymm4,%%ymm4							\n\t	vmovaps	0x20(%%r11),%%ymm11	\n\t"\
		"vaddpd	%%ymm2,%%ymm7,%%ymm7							\n\t	vmovaps	    (%%r12),%%ymm12	\n\t"\
		"vaddpd	%%ymm1,%%ymm5,%%ymm5							\n\t	vmovaps	0x20(%%r12),%%ymm13	\n\t"\
		"vaddpd	%%ymm3,%%ymm6,%%ymm6							\n\t/*	vmovaps	    (%%r13),%%ymm14	// r18.r */\n\t"\
		"/*vaps	%%ymm4,     (%%rsi)	// r08.r */					\n\t/*	vmovaps	0x20(%%r13),%%ymm15	// r18.i */\n\t"\
		"/*vaps	%%ymm7,     (%%rdi)	// r18.r */					\n\t	vsubpd	%%ymm7 ,%%ymm8 ,%%ymm8 		\n\t"\
		"/*vaps	%%ymm5,0x020(%%rsi)	// r08.i */					\n\t	vsubpd	%%ymm3 ,%%ymm9 ,%%ymm9 		\n\t"\
		"vmovaps	%%ymm6,0x2a0(%%rdi)	/* r38.i */				\n\t	vaddpd	%%ymm7 ,%%ymm7 ,%%ymm7 		\n\t"\
		"														\n\t	vaddpd	%%ymm3 ,%%ymm3 ,%%ymm3 		\n\t"\
		"/*	SSE2_RADIX_05_DFT_0TWIDDLE(r00...)	*/				\n\t	vaddpd	%%ymm8 ,%%ymm7 ,%%ymm7 		\n\t"\
		"movq	%[__r00],%%rsi									\n\t	vaddpd	%%ymm9 ,%%ymm3 ,%%ymm3 		\n\t"\
		"movq	%%rsi,%%rax										\n\t	vsubpd	%%ymm12,%%ymm10,%%ymm10		\n\t	vmovaps	%%ymm7 ,%%ymm14		\n\t"\
		"movq	%%rsi,%%rbx										\n\t	vsubpd	%%ymm13,%%ymm11,%%ymm11		\n\t	vmovaps	%%ymm3 ,%%ymm15		\n\t"\
		"movq	%%rsi,%%rcx		\n\t	vmovaps	%%ymm4 ,%%ymm6	\n\t	vaddpd	%%ymm12,%%ymm12,%%ymm12		\n\t"\
		"/*vq	%%rsi,%%rdx	*/	\n\t	vmovaps	%%ymm5 ,%%ymm7	\n\t	vaddpd	%%ymm13,%%ymm13,%%ymm13		\n\t"\
		"addq	$0x040,%%rax									\n\t	vaddpd	%%ymm10,%%ymm12,%%ymm12		\n\t"\
		"addq	$0x080,%%rbx									\n\t	vaddpd	%%ymm11,%%ymm13,%%ymm13		\n\t"\
		"addq	$0x0c0,%%rcx									\n\t	movq	%[__cc1],%%r10		\n\t"\
		"/*dq	$0x100,%%rdx	*/								\n\t	vsubpd	%%ymm12,%%ymm14,%%ymm14		\n\t"\
		"movq	%[__a0],%%rdi									\n\t	vsubpd	%%ymm13,%%ymm15,%%ymm15		\n\t"\
		"vmovaps	    (%%rax),%%ymm0							\n\t	vaddpd	%%ymm12,%%ymm12,%%ymm12		\n\t"\
		"vmovaps	0x20(%%rax),%%ymm1							\n\t	vaddpd	%%ymm13,%%ymm13,%%ymm13		\n\t"\
		"vmovaps	    (%%rbx),%%ymm2							\n\t	vaddpd	%%ymm14,%%ymm12,%%ymm12		\n\t"\
		"vmovaps	0x20(%%rbx),%%ymm3							\n\t	vaddpd	%%ymm15,%%ymm13,%%ymm13		\n\t"\
		"vmovaps	    (%%rcx),%%ymm4							\n\t	vaddpd	    (%%r8 ),%%ymm12,%%ymm12	\n\t"\
		"vmovaps	0x20(%%rcx),%%ymm5							\n\t	vaddpd	0x20(%%r8 ),%%ymm13,%%ymm13	\n\t"\
		"/*vaps	    (%%rdx),%%ymm6	// r08.r */					\n\t	vmovaps	%%ymm12,    (%%r9 )	\n\t"\
		"/*vaps	0x20(%%rdx),%%ymm7	// r08.i */					\n\t	vmovaps	%%ymm13,0x20(%%r9 )	\n\t"\
		"vsubpd	%%ymm6,%%ymm0,%%ymm0							\n\t	vmulpd	0x20(%%r10),%%ymm14,%%ymm14	\n\t"\
		"vsubpd	%%ymm7,%%ymm1,%%ymm1							\n\t	vmulpd	0x20(%%r10),%%ymm15,%%ymm15	\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6							\n\t	vsubpd	    (%%r8 ),%%ymm12,%%ymm12	\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7							\n\t	vsubpd	0x20(%%r8 ),%%ymm13,%%ymm13	\n\t"\
		"vaddpd	%%ymm0,%%ymm6,%%ymm6							\n\t	vmulpd	    (%%r10),%%ymm12,%%ymm12	\n\t"\
		"vaddpd	%%ymm1,%%ymm7,%%ymm7							\n\t	vmulpd	    (%%r10),%%ymm13,%%ymm13	\n\t"\
		"vsubpd	%%ymm4,%%ymm2,%%ymm2							\n\t	vaddpd	    (%%r9 ),%%ymm12,%%ymm12	\n\t"\
		"vsubpd	%%ymm5,%%ymm3,%%ymm3							\n\t	vaddpd	0x20(%%r9 ),%%ymm13,%%ymm13	\n\t"\
		"vmulpd	(%%r14),%%ymm4,%%ymm4							\n\t	vsubpd	%%ymm14,%%ymm12,%%ymm12		\n\t"\
		"vmulpd	(%%r14),%%ymm5,%%ymm5							\n\t	vsubpd	%%ymm15,%%ymm13,%%ymm13		\n\t"\
		"vaddpd	%%ymm2,%%ymm4,%%ymm4							\n\t	vmulpd	(%%r14),%%ymm14,%%ymm14		\n\t"\
		"vaddpd	%%ymm3,%%ymm5,%%ymm5							\n\t	vmulpd	(%%r14),%%ymm15,%%ymm15		\n\t"\
		"movq	%[__cc1],%%rax									\n\t	vaddpd	%%ymm12,%%ymm14,%%ymm14		\n\t"\
		"vsubpd	%%ymm4,%%ymm6,%%ymm6							\n\t	vaddpd	%%ymm13,%%ymm15,%%ymm15		\n\t"\
		"vsubpd	%%ymm5,%%ymm7,%%ymm7							\n\t	vmovaps	%%ymm12,    (%%r8 )	\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4							\n\t	vmovaps	%%ymm13,0x20(%%r8 )	\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5							\n\t	vmovaps	%%ymm8 ,%%ymm12		\n\t"\
		"vaddpd	%%ymm6,%%ymm4,%%ymm4							\n\t	vmovaps	%%ymm9 ,%%ymm13		\n\t"\
		"vaddpd	%%ymm7,%%ymm5,%%ymm5							\n\t	vsubpd	%%ymm10,%%ymm8 ,%%ymm8 		\n\t"\
		"vaddpd	    (%%rsi),%%ymm4,%%ymm4						\n\t	vsubpd	%%ymm11,%%ymm9 ,%%ymm9 		\n\t"\
		"vaddpd	0x20(%%rsi),%%ymm5,%%ymm5						\n\t	vmulpd	0x40(%%r10),%%ymm8 ,%%ymm8 	\n\t"\
		"vmovaps	%%ymm4,    (%%rdi)							\n\t	vmulpd	0x40(%%r10),%%ymm9 ,%%ymm9 	\n\t"\
		"vmovaps	%%ymm5,0x20(%%rdi)							\n\t	vmulpd	0x60(%%r10),%%ymm10,%%ymm10	\n\t"\
		"vmulpd	0x20(%%rax),%%ymm6,%%ymm6						\n\t	vmulpd	0x60(%%r10),%%ymm11,%%ymm11	\n\t"\
		"vmulpd	0x20(%%rax),%%ymm7,%%ymm7						\n\t	vmulpd	0x80(%%r10),%%ymm12,%%ymm12	\n\t"\
		"vsubpd	    (%%rsi),%%ymm4,%%ymm4						\n\t	vmulpd	0x80(%%r10),%%ymm13,%%ymm13	\n\t"\
		"vsubpd	0x20(%%rsi),%%ymm5,%%ymm5						\n\t	vaddpd	%%ymm8 ,%%ymm10,%%ymm10		\n\t"\
		"vmulpd	    (%%rax),%%ymm4,%%ymm4						\n\t	vaddpd	%%ymm9 ,%%ymm11,%%ymm11		\n\t"\
		"vmulpd	    (%%rax),%%ymm5,%%ymm5						\n\t	vsubpd	%%ymm12,%%ymm8 ,%%ymm8 		\n\t"\
		"vaddpd	    (%%rdi),%%ymm4,%%ymm4						\n\t	vsubpd	%%ymm13,%%ymm9 ,%%ymm9 		\n\t"\
		"vaddpd	0x20(%%rdi),%%ymm5,%%ymm5						\n\t	vmovaps	    (%%r8 ),%%ymm12	\n\t"\
		"vsubpd	%%ymm6,%%ymm4,%%ymm4							\n\t	vmovaps	0x20(%%r8 ),%%ymm13	\n\t"\
		"vsubpd	%%ymm7,%%ymm5,%%ymm5							\n\t	movq	%[__b1],%%r10		\n\t"\
		"vmulpd	(%%r14),%%ymm6,%%ymm6							\n\t	vsubpd	%%ymm11,%%ymm14,%%ymm14		\n\t"\
		"vmulpd	(%%r14),%%ymm7,%%ymm7							\n\t	vsubpd	%%ymm10,%%ymm15,%%ymm15		\n\t"\
		"vaddpd	%%ymm4,%%ymm6,%%ymm6							\n\t	vmulpd	(%%r14),%%ymm11,%%ymm11		\n\t"\
		"vaddpd	%%ymm5,%%ymm7,%%ymm7							\n\t	vmulpd	(%%r14),%%ymm10,%%ymm10		\n\t"\
		"vmovaps	%%ymm4,    (%%rsi)							\n\t	movq	%[__b4],%%r13		\n\t"\
		"vmovaps	%%ymm5,0x20(%%rsi)							\n\t	vmovaps	%%ymm14,    (%%r10)	\n\t"\
		"vmovaps	%%ymm0,%%ymm4								\n\t	vmovaps	%%ymm15,0x20(%%r13)	\n\t"\
		"vmovaps	%%ymm1,%%ymm5								\n\t	vaddpd	%%ymm14,%%ymm11,%%ymm11		\n\t"\
		"vsubpd	%%ymm2,%%ymm0,%%ymm0							\n\t	vaddpd	%%ymm15,%%ymm10,%%ymm10		\n\t"\
		"vsubpd	%%ymm3,%%ymm1,%%ymm1							\n\t	vmovaps	%%ymm11,    (%%r13)	\n\t"\
		"vmulpd	0x40(%%rax),%%ymm0,%%ymm0						\n\t	vmovaps	%%ymm10,0x20(%%r10)	\n\t"\
		"vmulpd	0x40(%%rax),%%ymm1,%%ymm1						\n\t	movq	%[__b2],%%r11		\n\t"\
		"vmulpd	0x60(%%rax),%%ymm2,%%ymm2						\n\t	movq	%[__b3],%%r12		\n\t"\
		"vmulpd	0x60(%%rax),%%ymm3,%%ymm3						\n\t	\n\t"\
		"vmulpd	0x80(%%rax),%%ymm4,%%ymm4						\n\t	vsubpd	%%ymm9 ,%%ymm12,%%ymm12		\n\t"\
		"vmulpd	0x80(%%rax),%%ymm5,%%ymm5						\n\t	vsubpd	%%ymm8 ,%%ymm13,%%ymm13		\n\t"\
		"vaddpd	%%ymm0,%%ymm2,%%ymm2							\n\t	vmulpd	(%%r14),%%ymm9 ,%%ymm9 		\n\t"\
		"vaddpd	%%ymm1,%%ymm3,%%ymm3							\n\t	vmulpd	(%%r14),%%ymm8 ,%%ymm8 		\n\t"\
		"vsubpd	%%ymm4,%%ymm0,%%ymm0							\n\t	vmovaps	%%ymm12,    (%%r11)	\n\t"\
		"vsubpd	%%ymm5,%%ymm1,%%ymm1							\n\t	vmovaps	%%ymm13,0x20(%%r12)	\n\t"\
		"vmovaps	    (%%rsi),%%ymm4							\n\t	vaddpd	%%ymm12,%%ymm9 ,%%ymm9 		\n\t"\
		"vmovaps	0x20(%%rsi),%%ymm5							\n\t	vaddpd	%%ymm13,%%ymm8 ,%%ymm8 		\n\t"\
		"movq	%[__a1],%%rax									\n\t	vmovaps	%%ymm9 ,    (%%r12)	\n\t"\
		"movq	%[__a4],%%rdx									\n\t	vmovaps	%%ymm8 ,0x20(%%r11)	\n\t"\
		"vsubpd	%%ymm3,%%ymm6,%%ymm6							\n\t	/*	SSE2_RADIX_05_DFT_0TWIDDLE(r30)	*/\n\t"\
		"vsubpd	%%ymm2,%%ymm7,%%ymm7							\n\t	movq	%[__r30],%%r8 		\n\t"\
		"vaddpd	%%ymm3,%%ymm3,%%ymm3							\n\t	movq	%%r8 ,%%r10			\n\t"\
		"vaddpd	%%ymm2,%%ymm2,%%ymm2							\n\t	movq	%%r8 ,%%r11			\n\t"\
		"vmovaps	%%ymm6,    (%%rax)							\n\t	movq	%%r8 ,%%r12			\n\t"\
		"vmovaps	%%ymm7,0x20(%%rdx)							\n\t	movq	%%r8 ,%%r13			\n\t"\
		"vaddpd	%%ymm6,%%ymm3,%%ymm3							\n\t	addq	$0x040,%%r10		\n\t"\
		"vaddpd	%%ymm7,%%ymm2,%%ymm2							\n\t	addq	$0x080,%%r11		\n\t"\
		"vmovaps	%%ymm3,    (%%rdx)							\n\t	addq	$0x0c0,%%r12		\n\t"\
		"vmovaps	%%ymm2,0x20(%%rax)							\n\t	addq	$0x100,%%r13		\n\t"\
		"movq	%[__a2],%%rbx									\n\t	movq	%[__d0],%%r9 		\n\t"\
		"movq	%[__a3],%%rcx									\n\t	vmovaps	    (%%r10),%%ymm8 	\n\t"\
		"vsubpd	%%ymm1,%%ymm4,%%ymm4							\n\t	vmovaps	0x20(%%r10),%%ymm9 	\n\t"\
		"vsubpd	%%ymm0,%%ymm5,%%ymm5							\n\t	vmovaps	    (%%r11),%%ymm10	\n\t"\
		"vaddpd	%%ymm1,%%ymm1,%%ymm1							\n\t	vmovaps	0x20(%%r11),%%ymm11	\n\t"\
		"vaddpd	%%ymm0,%%ymm0,%%ymm0							\n\t	vmovaps	    (%%r12),%%ymm12	\n\t"\
		"vmovaps	%%ymm4,    (%%rbx)							\n\t	vmovaps	0x20(%%r12),%%ymm13	\n\t"\
		"vmovaps	%%ymm5,0x20(%%rcx)							\n\t	vmovaps	    (%%r13),%%ymm14	\n\t"\
		"vaddpd	%%ymm4,%%ymm1,%%ymm1							\n\t	vmovaps	0x20(%%r13),%%ymm15	\n\t"\
		"vaddpd	%%ymm5,%%ymm0,%%ymm0							\n\t	vsubpd	%%ymm14,%%ymm8 ,%%ymm8 		\n\t"\
		"vmovaps	%%ymm1,    (%%rcx)							\n\t	vsubpd	%%ymm15,%%ymm9 ,%%ymm9 		\n\t"\
		"vmovaps	%%ymm0,0x20(%%rbx)							\n\t	vaddpd	%%ymm14,%%ymm14,%%ymm14		\n\t"\
		"														\n\t	vaddpd	%%ymm15,%%ymm15,%%ymm15		\n\t"\
		"/*	SSE2_RADIX_05_DFT_0TWIDDLE(r20)	*/					\n\t	vaddpd	%%ymm8 ,%%ymm14,%%ymm14		\n\t"\
		"movq	%[__r20],%%rsi									\n\t	vaddpd	%%ymm9 ,%%ymm15,%%ymm15		\n\t"\
		"movq	%%rsi,%%rax										\n\t	vsubpd	%%ymm12,%%ymm10,%%ymm10		\n\t"\
		"movq	%%rsi,%%rbx										\n\t	vsubpd	%%ymm13,%%ymm11,%%ymm11		\n\t"\
		"movq	%%rsi,%%rcx										\n\t	vaddpd	%%ymm12,%%ymm12,%%ymm12		\n\t"\
		"movq	%%rsi,%%rdx										\n\t	vaddpd	%%ymm13,%%ymm13,%%ymm13		\n\t"\
		"addq	$0x040,%%rax									\n\t	vaddpd	%%ymm10,%%ymm12,%%ymm12		\n\t"\
		"addq	$0x080,%%rbx									\n\t	vaddpd	%%ymm11,%%ymm13,%%ymm13		\n\t"\
		"addq	$0x0c0,%%rcx									\n\t	movq	%[__cc1],%%r10		\n\t"\
		"addq	$0x100,%%rdx									\n\t	vsubpd	%%ymm12,%%ymm14,%%ymm14		\n\t"\
		"movq	%[__c0],%%rdi									\n\t	vsubpd	%%ymm13,%%ymm15,%%ymm15		\n\t"\
		"vmovaps	    (%%rax),%%ymm0							\n\t	vaddpd	%%ymm12,%%ymm12,%%ymm12		\n\t"\
		"vmovaps	0x20(%%rax),%%ymm1							\n\t	vaddpd	%%ymm13,%%ymm13,%%ymm13		\n\t"\
		"vmovaps	    (%%rbx),%%ymm2							\n\t	vaddpd	%%ymm14,%%ymm12,%%ymm12		\n\t"\
		"vmovaps	0x20(%%rbx),%%ymm3							\n\t	vaddpd	%%ymm15,%%ymm13,%%ymm13		\n\t"\
		"vmovaps	    (%%rcx),%%ymm4							\n\t	vaddpd	    (%%r8 ),%%ymm12,%%ymm12	\n\t"\
		"vmovaps	0x20(%%rcx),%%ymm5							\n\t	vaddpd	0x20(%%r8 ),%%ymm13,%%ymm13	\n\t"\
		"vmovaps	    (%%rdx),%%ymm6							\n\t	vmovaps	%%ymm12,    (%%r9 )	\n\t"\
		"vmovaps	0x20(%%rdx),%%ymm7							\n\t	vmovaps	%%ymm13,0x20(%%r9 )	\n\t"\
		"vsubpd	%%ymm6,%%ymm0,%%ymm0							\n\t	vmulpd	0x20(%%r10),%%ymm14,%%ymm14	\n\t"\
		"vsubpd	%%ymm7,%%ymm1,%%ymm1							\n\t	vmulpd	0x20(%%r10),%%ymm15,%%ymm15	\n\t"\
		"vmulpd	(%%r14),%%ymm6,%%ymm6							\n\t	vsubpd	    (%%r8 ),%%ymm12,%%ymm12	\n\t"\
		"vmulpd	(%%r14),%%ymm7,%%ymm7							\n\t	vsubpd	0x20(%%r8 ),%%ymm13,%%ymm13	\n\t"\
		"vaddpd	%%ymm0,%%ymm6,%%ymm6							\n\t	vmulpd	    (%%r10),%%ymm12,%%ymm12	\n\t"\
		"vaddpd	%%ymm1,%%ymm7,%%ymm7							\n\t	vmulpd	    (%%r10),%%ymm13,%%ymm13	\n\t"\
		"vsubpd	%%ymm4,%%ymm2,%%ymm2							\n\t	vaddpd	    (%%r9 ),%%ymm12,%%ymm12	\n\t"\
		"vsubpd	%%ymm5,%%ymm3,%%ymm3							\n\t	vaddpd	0x20(%%r9 ),%%ymm13,%%ymm13	\n\t"\
		"vmulpd	(%%r14),%%ymm4,%%ymm4							\n\t	vsubpd	%%ymm14,%%ymm12,%%ymm12		\n\t"\
		"vmulpd	(%%r14),%%ymm5,%%ymm5							\n\t	vsubpd	%%ymm15,%%ymm13,%%ymm13		\n\t"\
		"vaddpd	%%ymm2,%%ymm4,%%ymm4							\n\t	vmulpd	(%%r14),%%ymm14,%%ymm14		\n\t"\
		"vaddpd	%%ymm3,%%ymm5,%%ymm5							\n\t	vmulpd	(%%r14),%%ymm15,%%ymm15		\n\t"\
		"movq	%[__cc1],%%rax									\n\t	vaddpd	%%ymm12,%%ymm14,%%ymm14		\n\t"\
		"vsubpd	%%ymm4,%%ymm6,%%ymm6							\n\t	vaddpd	%%ymm13,%%ymm15,%%ymm15		\n\t"\
		"vsubpd	%%ymm5,%%ymm7,%%ymm7							\n\t	vmovaps	%%ymm12,    (%%r8 )	\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4							\n\t	vmovaps	%%ymm13,0x20(%%r8 )	\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5							\n\t	vmovaps	%%ymm8 ,%%ymm12		\n\t"\
		"vaddpd	%%ymm6,%%ymm4,%%ymm4							\n\t	vmovaps	%%ymm9 ,%%ymm13		\n\t"\
		"vaddpd	%%ymm7,%%ymm5,%%ymm5							\n\t	vsubpd	%%ymm10,%%ymm8 ,%%ymm8 		\n\t"\
		"vaddpd		(%%rsi),%%ymm4,%%ymm4						\n\t	vsubpd	%%ymm11,%%ymm9 ,%%ymm9 		\n\t"\
		"vaddpd	0x20(%%rsi),%%ymm5,%%ymm5						\n\t	vmulpd	0x40(%%r10),%%ymm8 ,%%ymm8 	\n\t"\
		"vmovaps	%%ymm4,    (%%rdi)							\n\t	vmulpd	0x40(%%r10),%%ymm9 ,%%ymm9 	\n\t"\
		"vmovaps	%%ymm5,0x20(%%rdi)							\n\t	vmulpd	0x60(%%r10),%%ymm10,%%ymm10	\n\t"\
		"vmulpd	0x20(%%rax),%%ymm6,%%ymm6						\n\t	vmulpd	0x60(%%r10),%%ymm11,%%ymm11	\n\t"\
		"vmulpd	0x20(%%rax),%%ymm7,%%ymm7						\n\t	vmulpd	0x80(%%r10),%%ymm12,%%ymm12	\n\t"\
		"vsubpd		(%%rsi),%%ymm4,%%ymm4						\n\t	vmulpd	0x80(%%r10),%%ymm13,%%ymm13	\n\t"\
		"vsubpd	0x20(%%rsi),%%ymm5,%%ymm5						\n\t	vaddpd	%%ymm8 ,%%ymm10,%%ymm10		\n\t"\
		"vmulpd		(%%rax),%%ymm4,%%ymm4						\n\t	vaddpd	%%ymm9 ,%%ymm11,%%ymm11		\n\t"\
		"vmulpd		(%%rax),%%ymm5,%%ymm5						\n\t	vsubpd	%%ymm12,%%ymm8 ,%%ymm8 		\n\t"\
		"vaddpd		(%%rdi),%%ymm4,%%ymm4						\n\t	vsubpd	%%ymm13,%%ymm9 ,%%ymm9 		\n\t"\
		"vaddpd	0x20(%%rdi),%%ymm5,%%ymm5						\n\t	vmovaps	    (%%r8 ),%%ymm12	\n\t"\
		"vsubpd	%%ymm6,%%ymm4,%%ymm4							\n\t	vmovaps	0x20(%%r8 ),%%ymm13	\n\t"\
		"vsubpd	%%ymm7,%%ymm5,%%ymm5							\n\t	movq	%[__d1],%%r10		\n\t"\
		"vmulpd	(%%r14),%%ymm6,%%ymm6							\n\t	vsubpd	%%ymm11,%%ymm14,%%ymm14		\n\t"\
		"vmulpd	(%%r14),%%ymm7,%%ymm7							\n\t	vsubpd	%%ymm10,%%ymm15,%%ymm15		\n\t"\
		"vaddpd	%%ymm4,%%ymm6,%%ymm6							\n\t	vmulpd	(%%r14),%%ymm11,%%ymm11		\n\t"\
		"vaddpd	%%ymm5,%%ymm7,%%ymm7							\n\t	vmulpd	(%%r14),%%ymm10,%%ymm10		\n\t"\
		"vmovaps	%%ymm4,    (%%rsi)							\n\t	movq	%[__d4],%%r13		\n\t"\
		"vmovaps	%%ymm5,0x20(%%rsi)							\n\t	vmovaps	%%ymm14,    (%%r10)	\n\t"\
		"vmovaps	%%ymm0,%%ymm4								\n\t	vmovaps	%%ymm15,0x20(%%r13)	\n\t"\
		"vmovaps	%%ymm1,%%ymm5								\n\t	vaddpd	%%ymm14,%%ymm11,%%ymm11		\n\t"\
		"vsubpd	%%ymm2,%%ymm0,%%ymm0							\n\t	vaddpd	%%ymm15,%%ymm10,%%ymm10		\n\t"\
		"vsubpd	%%ymm3,%%ymm1,%%ymm1							\n\t	vmovaps	%%ymm11,    (%%r13)	\n\t"\
		"vmulpd	0x40(%%rax),%%ymm0,%%ymm0						\n\t	vmovaps	%%ymm10,0x20(%%r10)	\n\t"\
		"vmulpd	0x40(%%rax),%%ymm1,%%ymm1						\n\t	movq	%[__d2],%%r11		\n\t"\
		"vmulpd	0x60(%%rax),%%ymm2,%%ymm2						\n\t	movq	%[__d3],%%r12		\n\t"\
		"vmulpd	0x60(%%rax),%%ymm3,%%ymm3						\n\t	vsubpd	%%ymm9 ,%%ymm12,%%ymm12		\n\t"\
		"vmulpd	0x80(%%rax),%%ymm4,%%ymm4						\n\t	vsubpd	%%ymm8 ,%%ymm13,%%ymm13		\n\t"\
		"vmulpd	0x80(%%rax),%%ymm5,%%ymm5						\n\t	\n\t"\
		"vaddpd	%%ymm0,%%ymm2,%%ymm2							\n\t	vmulpd	(%%r14),%%ymm9 ,%%ymm9 		\n\t"\
		"vaddpd	%%ymm1,%%ymm3,%%ymm3							\n\t	vmulpd	(%%r14),%%ymm8 ,%%ymm8 		\n\t"\
		"vsubpd	%%ymm4,%%ymm0,%%ymm0							\n\t	vmovaps	%%ymm12,    (%%r11)	\n\t"\
		"vsubpd	%%ymm5,%%ymm1,%%ymm1							\n\t	vmovaps	%%ymm13,0x20(%%r12)	\n\t"\
		"vmovaps	    (%%rsi),%%ymm4							\n\t	vaddpd	%%ymm12,%%ymm9 ,%%ymm9 		\n\t"\
		"vmovaps	0x20(%%rsi),%%ymm5							\n\t	vaddpd	%%ymm13,%%ymm8 ,%%ymm8 		\n\t"\
		"movq	%[__c1],%%rax									\n\t	vmovaps	%%ymm9 ,    (%%r12)	\n\t"\
		"movq	%[__c4],%%rdx									\n\t	vmovaps	%%ymm8 ,0x20(%%r11)	\n\t"\
		"vsubpd	%%ymm3,%%ymm6,%%ymm6							\n\t"\
		"vsubpd	%%ymm2,%%ymm7,%%ymm7							\n\t"\
		"vaddpd	%%ymm3,%%ymm3,%%ymm3							\n\t"\
		"vaddpd	%%ymm2,%%ymm2,%%ymm2							\n\t"\
		"vmovaps	%%ymm6,    (%%rax)							\n\t"\
		"vmovaps	%%ymm7,0x20(%%rdx)							\n\t"\
		"vaddpd	%%ymm6,%%ymm3,%%ymm3							\n\t"\
		"vaddpd	%%ymm7,%%ymm2,%%ymm2							\n\t"\
		"vmovaps	%%ymm3,    (%%rdx)							\n\t"\
		"vmovaps	%%ymm2,0x20(%%rax)							\n\t"\
		"movq	%[__c2],%%rbx									\n\t"\
		"movq	%[__c3],%%rcx									\n\t"\
		"vsubpd	%%ymm1,%%ymm4,%%ymm4							\n\t"\
		"vsubpd	%%ymm0,%%ymm5,%%ymm5							\n\t"\
		"vaddpd	%%ymm1,%%ymm1,%%ymm1							\n\t"\
		"vaddpd	%%ymm0,%%ymm0,%%ymm0							\n\t"\
		"vmovaps	%%ymm4,    (%%rbx)							\n\t"\
		"vmovaps	%%ymm5,0x20(%%rcx)							\n\t"\
		"vaddpd	%%ymm4,%%ymm1,%%ymm1							\n\t"\
		"vaddpd	%%ymm5,%%ymm0,%%ymm0							\n\t"\
		"vmovaps	%%ymm1,    (%%rcx)							\n\t"\
		"vmovaps	%%ymm0,0x20(%%rbx)							\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		,[__p01] "m" (Xp01)\
		,[__p04] "m" (Xp04)\
		,[__r00] "m" (Xr00)\
		,[__r10] "m" (Xr10)\
		,[__r20] "m" (Xr20)\
		,[__r30] "m" (Xr30)\
		,[__cc1] "m" (Xcc1)\
		,[__a0] "m" (Xa0)\
		,[__a1] "m" (Xa1)\
		,[__a2] "m" (Xa2)\
		,[__a3] "m" (Xa3)\
		,[__a4] "m" (Xa4)\
		,[__b0] "m" (Xb0)\
		,[__b1] "m" (Xb1)\
		,[__b2] "m" (Xb2)\
		,[__b3] "m" (Xb3)\
		,[__b4] "m" (Xb4)\
		,[__c0] "m" (Xc0)\
		,[__c1] "m" (Xc1)\
		,[__c2] "m" (Xc2)\
		,[__c3] "m" (Xc3)\
		,[__c4] "m" (Xc4)\
		,[__d0] "m" (Xd0)\
		,[__d1] "m" (Xd1)\
		,[__d2] "m" (Xd2)\
		,[__d3] "m" (Xd3)\
		,[__d4] "m" (Xd4)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r8","r9","r10","r11","r12","r13","r14","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
	);\
	}


	#define	SSE2_RADIX20_DIF_NOTWIDDLE(Xadd0,Xp01,Xp04,Xp08,Xp16,Xa0,Xa1,Xa2,Xa3,Xa4,Xb0,Xb1,Xb2,Xb3,Xb4,Xc0,Xc1,Xc2,Xc3,Xc4,Xd0,Xd1,Xd2,Xd3,Xd4,Xcc1,Xr00,Xr10,Xr20,Xr30)\
	{\
	__asm__ volatile (\
		"/*	SSE2_RADIX_05_DFT_0TWIDDLE(Xa0,r00*/\n\t"\
		"movq	%[__r00],%%rdi					\n\t"\
		"movq	%[__a0],%%rsi					\n\t"\
		"movq	%[__a1],%%rax					\n\t	movq	$0xa00,%%r14		\n\t"\
		"movq	%[__a2],%%rbx					\n\t"\
		"movq	%[__a3],%%rcx					\n\t"\
		"movq	%[__a4],%%rdx					\n\t	addq	%%rdi,%%r14	/* two */\n\t"\
		"vmovaps	    (%%rax),%%ymm0			\n\t"\
		"vmovaps	0x20(%%rax),%%ymm1			\n\t	/*	SSE2_RADIX_05_DFT_0TWIDDLE(Xc0,r20)	*/\n\t"\
		"vmovaps	    (%%rbx),%%ymm2			\n\t	movq	%[__c0],%%r9 		\n\t"\
		"vmovaps	0x20(%%rbx),%%ymm3			\n\t	movq	%[__c1],%%r10		\n\t"\
		"vmovaps	    (%%rcx),%%ymm4			\n\t	movq	%[__c2],%%r11		\n\t"\
		"vmovaps	0x20(%%rcx),%%ymm5			\n\t	movq	%[__c3],%%r12		\n\t"\
		"vmovaps	    (%%rdx),%%ymm6			\n\t	movq	%[__c4],%%r13		\n\t"\
		"vmovaps	0x20(%%rdx),%%ymm7			\n\t	movq	%[__r20],%%r8 		\n\t"\
		"vsubpd	%%ymm6,%%ymm0,%%ymm0			\n\t	vmovaps	    (%%r10),%%ymm8 	\n\t"\
		"vsubpd	%%ymm7,%%ymm1,%%ymm1			\n\t	vmovaps	0x20(%%r10),%%ymm9 	\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6			\n\t	vmovaps	    (%%r11),%%ymm10	\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7			\n\t	vmovaps	0x20(%%r11),%%ymm11	\n\t"\
		"vaddpd	%%ymm0,%%ymm6,%%ymm6			\n\t	vmovaps	    (%%r12),%%ymm12	\n\t"\
		"vaddpd	%%ymm1,%%ymm7,%%ymm7			\n\t	vmovaps	0x20(%%r12),%%ymm13	\n\t"\
		"vsubpd	%%ymm4,%%ymm2,%%ymm2			\n\t	vmovaps	    (%%r13),%%ymm14	\n\t"\
		"vsubpd	%%ymm5,%%ymm3,%%ymm3			\n\t	vmovaps	0x20(%%r13),%%ymm15	\n\t"\
		"vmulpd	(%%r14),%%ymm4,%%ymm4			\n\t	vsubpd	%%ymm14,%%ymm8 ,%%ymm8 		\n\t"\
		"vmulpd	(%%r14),%%ymm5,%%ymm5			\n\t	vsubpd	%%ymm15,%%ymm9 ,%%ymm9 		\n\t"\
		"vaddpd	%%ymm2,%%ymm4,%%ymm4			\n\t	vmulpd	(%%r14),%%ymm14,%%ymm14		\n\t"\
		"vaddpd	%%ymm3,%%ymm5,%%ymm5			\n\t	vmulpd	(%%r14),%%ymm15,%%ymm15		\n\t"\
		"movq	%[__cc1],%%rax					\n\t	vaddpd	%%ymm8 ,%%ymm14,%%ymm14		\n\t"\
		"vsubpd	%%ymm4,%%ymm6,%%ymm6			\n\t	vaddpd	%%ymm9 ,%%ymm15,%%ymm15		\n\t"\
		"vsubpd	%%ymm5,%%ymm7,%%ymm7			\n\t	vsubpd	%%ymm12,%%ymm10,%%ymm10		\n\t"\
		"vmulpd	(%%r14),%%ymm4,%%ymm4			\n\t	vsubpd	%%ymm13,%%ymm11,%%ymm11		\n\t"\
		"vmulpd	(%%r14),%%ymm5,%%ymm5			\n\t	movq	%[__cc1],%%r10		\n\t"\
		"vaddpd	%%ymm6,%%ymm4,%%ymm4			\n\t	vmulpd	(%%r14),%%ymm12,%%ymm12		\n\t"\
		"vaddpd	%%ymm7,%%ymm5,%%ymm5			\n\t	vmulpd	(%%r14),%%ymm13,%%ymm13		\n\t"\
		"vaddpd	    (%%rsi),%%ymm4,%%ymm4		\n\t	vaddpd	%%ymm10,%%ymm12,%%ymm12		\n\t"\
		"vaddpd	0x20(%%rsi),%%ymm5,%%ymm5		\n\t	vaddpd	%%ymm11,%%ymm13,%%ymm13		\n\t"\
		"vmovaps	%%ymm4,    (%%rdi)			\n\t	vsubpd	%%ymm12,%%ymm14,%%ymm14		\n\t"\
		"vmovaps	%%ymm5,0x20(%%rdi)			\n\t	vsubpd	%%ymm13,%%ymm15,%%ymm15		\n\t"\
		"vmulpd	0x20(%%rax),%%ymm6,%%ymm6		\n\t	vaddpd	%%ymm12,%%ymm12,%%ymm12		\n\t"\
		"vmulpd	0x20(%%rax),%%ymm7,%%ymm7		\n\t	vaddpd	%%ymm13,%%ymm13,%%ymm13		\n\t"\
		"vsubpd		(%%rsi),%%ymm4,%%ymm4		\n\t	vaddpd	%%ymm14,%%ymm12,%%ymm12		\n\t"\
		"vsubpd	0x20(%%rsi),%%ymm5,%%ymm5		\n\t	vaddpd	%%ymm15,%%ymm13,%%ymm13		\n\t"\
		"vmulpd		(%%rax),%%ymm4,%%ymm4		\n\t	vaddpd	    (%%r9 ),%%ymm12,%%ymm12	\n\t"\
		"vmulpd		(%%rax),%%ymm5,%%ymm5		\n\t	vaddpd	0x20(%%r9 ),%%ymm13,%%ymm13	\n\t"\
		"vaddpd		(%%rdi),%%ymm4,%%ymm4		\n\t	vmovaps	%%ymm12,    (%%r8 )	\n\t"\
		"vaddpd	0x20(%%rdi),%%ymm5,%%ymm5		\n\t	vmovaps	%%ymm13,0x20(%%r8 )	\n\t"\
		"vsubpd	%%ymm6,%%ymm4,%%ymm4			\n\t	vmulpd	0x20(%%r10),%%ymm14,%%ymm14	\n\t"\
		"vsubpd	%%ymm7,%%ymm5,%%ymm5			\n\t	vmulpd	0x20(%%r10),%%ymm15,%%ymm15	\n\t"\
		"vmulpd	(%%r14),%%ymm6,%%ymm6			\n\t	vsubpd	    (%%r9 ),%%ymm12,%%ymm12	\n\t"\
		"vmulpd	(%%r14),%%ymm7,%%ymm7			\n\t	vsubpd	0x20(%%r9 ),%%ymm13,%%ymm13	\n\t"\
		"vaddpd	%%ymm4,%%ymm6,%%ymm6			\n\t	vmulpd	    (%%r10),%%ymm12,%%ymm12	\n\t"\
		"vaddpd	%%ymm5,%%ymm7,%%ymm7			\n\t	vmulpd	    (%%r10),%%ymm13,%%ymm13	\n\t"\
		"vmovaps	%%ymm4,    (%%rsi)			\n\t	vaddpd	    (%%r8 ),%%ymm12,%%ymm12	\n\t"\
		"vmovaps	%%ymm5,0x20(%%rsi)			\n\t	vaddpd	0x20(%%r8 ),%%ymm13,%%ymm13	\n\t"\
		"vmovaps	%%ymm0,%%ymm4				\n\t	vsubpd	%%ymm14,%%ymm12,%%ymm12		\n\t"\
		"vmovaps	%%ymm1,%%ymm5				\n\t	vsubpd	%%ymm15,%%ymm13,%%ymm13		\n\t"\
		"vsubpd	%%ymm2,%%ymm0,%%ymm0			\n\t	vmulpd	(%%r14),%%ymm14,%%ymm14		\n\t"\
		"vsubpd	%%ymm3,%%ymm1,%%ymm1			\n\t	vmulpd	(%%r14),%%ymm15,%%ymm15		\n\t"\
		"vmulpd	0x40(%%rax),%%ymm0,%%ymm0		\n\t	vaddpd	%%ymm12,%%ymm14,%%ymm14		\n\t"\
		"vmulpd	0x40(%%rax),%%ymm1,%%ymm1		\n\t	vaddpd	%%ymm13,%%ymm15,%%ymm15		\n\t"\
		"vmulpd	0x60(%%rax),%%ymm2,%%ymm2		\n\t	vmovaps	%%ymm12,    (%%r9 )	\n\t"\
		"vmulpd	0x60(%%rax),%%ymm3,%%ymm3		\n\t	vmovaps	%%ymm13,0x20(%%r9 )	\n\t"\
		"vmulpd	0x80(%%rax),%%ymm4,%%ymm4		\n\t	vmovaps	%%ymm8 ,%%ymm12		\n\t"\
		"vmulpd	0x80(%%rax),%%ymm5,%%ymm5		\n\t	vmovaps	%%ymm9 ,%%ymm13		\n\t"\
		"vaddpd	%%ymm0,%%ymm2,%%ymm2			\n\t	vsubpd	%%ymm10,%%ymm8 ,%%ymm8 		\n\t"\
		"vaddpd	%%ymm1,%%ymm3,%%ymm3			\n\t	vsubpd	%%ymm11,%%ymm9 ,%%ymm9 		\n\t"\
		"vsubpd	%%ymm4,%%ymm0,%%ymm0			\n\t	vmulpd	0x40(%%r10),%%ymm8 ,%%ymm8 	\n\t"\
		"vsubpd	%%ymm5,%%ymm1,%%ymm1			\n\t	vmulpd	0x40(%%r10),%%ymm9 ,%%ymm9 	\n\t"\
		"vmovaps	    (%%rsi),%%ymm4			\n\t	vmulpd	0x60(%%r10),%%ymm10,%%ymm10	\n\t"\
		"vmovaps	0x20(%%rsi),%%ymm5			\n\t	vmulpd	0x60(%%r10),%%ymm11,%%ymm11	\n\t"\
		"vsubpd	%%ymm3,%%ymm6,%%ymm6			\n\t	vmulpd	0x80(%%r10),%%ymm12,%%ymm12	\n\t"\
		"vsubpd	%%ymm2,%%ymm7,%%ymm7			\n\t	vmulpd	0x80(%%r10),%%ymm13,%%ymm13	\n\t"\
		"vaddpd	%%ymm3,%%ymm3,%%ymm3			\n\t	vaddpd	%%ymm8 ,%%ymm10,%%ymm10		\n\t"\
		"vaddpd	%%ymm2,%%ymm2,%%ymm2			\n\t	vaddpd	%%ymm9 ,%%ymm11,%%ymm11		\n\t"\
		"vmovaps	%%ymm6,0x040(%%rdi)			\n\t	vsubpd	%%ymm12,%%ymm8 ,%%ymm8 		\n\t"\
		"vmovaps	%%ymm7,0x0e0(%%rdi)			\n\t	vsubpd	%%ymm13,%%ymm9 ,%%ymm9 		\n\t"\
		"vaddpd	%%ymm6,%%ymm3,%%ymm3			\n\t	vmovaps	    (%%r9 ),%%ymm12	\n\t"\
		"vaddpd	%%ymm7,%%ymm2,%%ymm2			\n\t	vmovaps	0x20(%%r9 ),%%ymm13	\n\t"\
		"vmovaps	%%ymm3,0x0c0(%%rdi)			\n\t	vsubpd	%%ymm11,%%ymm14,%%ymm14		\n\t"\
		"vmovaps	%%ymm2,0x060(%%rdi)			\n\t	vsubpd	%%ymm10,%%ymm15,%%ymm15		\n\t"\
		"vsubpd	%%ymm1,%%ymm4,%%ymm4			\n\t	vaddpd	%%ymm11,%%ymm11,%%ymm11		\n\t"\
		"vsubpd	%%ymm0,%%ymm5,%%ymm5			\n\t	vaddpd	%%ymm10,%%ymm10,%%ymm10		\n\t"\
		"vaddpd	%%ymm1,%%ymm1,%%ymm1			\n\t	vmovaps	%%ymm14,0x040(%%r8 )\n\t"\
		"vaddpd	%%ymm0,%%ymm0,%%ymm0			\n\t	vmovaps	%%ymm15,0x0e0(%%r8 )\n\t"\
		"vmovaps	%%ymm4,0x080(%%rdi)			\n\t	vaddpd	%%ymm14,%%ymm11,%%ymm11		\n\t"\
		"vmovaps	%%ymm5,0x120(%%rdi)			\n\t	vaddpd	%%ymm15,%%ymm10,%%ymm10		\n\t"\
		"vaddpd	%%ymm4,%%ymm1,%%ymm1			\n\t	vmovaps	%%ymm11,0x0c0(%%r8 )\n\t"\
		"vaddpd	%%ymm5,%%ymm0,%%ymm0			\n\t	vmovaps	%%ymm10,0x060(%%r8 )\n\t"\
		"vmovaps	%%ymm1,0x100(%%rdi)			\n\t	vsubpd	%%ymm9 ,%%ymm12,%%ymm12		\n\t"\
		"vmovaps	%%ymm0,0x0a0(%%rdi)			\n\t	vsubpd	%%ymm8 ,%%ymm13,%%ymm13		\n\t"\
		"										\n\t	vaddpd	%%ymm9 ,%%ymm9 ,%%ymm9 		\n\t"\
		"/*	SSE2_RADIX_05_DFT_0TWIDDLE(Xb0,r10*/\n\t	vaddpd	%%ymm8 ,%%ymm8 ,%%ymm8 		\n\t"\
		"movq	%[__b0],%%rsi					\n\t	vmovaps	%%ymm12,0x080(%%r8 )\n\t"\
		"movq	%[__b1],%%rax					\n\t	vmovaps	%%ymm13,0x120(%%r8 )\n\t"\
		"movq	%[__b2],%%rbx					\n\t	vaddpd	%%ymm12,%%ymm9 ,%%ymm9 		\n\t"\
		"movq	%[__b3],%%rcx					\n\t	vaddpd	%%ymm13,%%ymm8 ,%%ymm8 		\n\t"\
		"movq	%[__b4],%%rdx					\n\t	vmovaps	%%ymm9 ,0x100(%%r8 )\n\t"\
		"movq	%[__r10],%%rdi					\n\t	vmovaps	%%ymm8 ,0x0a0(%%r8 )\n\t"\
		"vmovaps	    (%%rax),%%ymm0			\n\t	\n\t"\
		"vmovaps	0x20(%%rax),%%ymm1			\n\t	/*	SSE2_RADIX_05_DFT_0TWIDDLE(Xd0,r30)	*/\n\t"\
		"vmovaps	    (%%rbx),%%ymm2			\n\t	movq	%[__d0],%%r9 		\n\t"\
		"vmovaps	0x20(%%rbx),%%ymm3			\n\t	movq	%[__d1],%%r10		\n\t"\
		"vmovaps	    (%%rcx),%%ymm4			\n\t	movq	%[__d2],%%r11		\n\t"\
		"vmovaps	0x20(%%rcx),%%ymm5			\n\t	movq	%[__d3],%%r12		\n\t"\
		"vmovaps	    (%%rdx),%%ymm6			\n\t	movq	%[__d4],%%r13		\n\t"\
		"vmovaps	0x20(%%rdx),%%ymm7			\n\t	movq	%[__r30],%%r8 		\n\t"\
		"vsubpd	%%ymm6,%%ymm0,%%ymm0			\n\t	vmovaps	    (%%r10),%%ymm8 	\n\t"\
		"vsubpd	%%ymm7,%%ymm1,%%ymm1			\n\t	vmovaps	0x20(%%r10),%%ymm9 	\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6			\n\t	vmovaps	    (%%r11),%%ymm10	\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7			\n\t	vmovaps	0x20(%%r11),%%ymm11	\n\t"\
		"vaddpd	%%ymm0,%%ymm6,%%ymm6			\n\t	vmovaps	    (%%r12),%%ymm12	\n\t"\
		"vaddpd	%%ymm1,%%ymm7,%%ymm7			\n\t	vmovaps	0x20(%%r12),%%ymm13	\n\t"\
		"vsubpd	%%ymm4,%%ymm2,%%ymm2			\n\t	vmovaps	    (%%r13),%%ymm14	\n\t"\
		"vsubpd	%%ymm5,%%ymm3,%%ymm3			\n\t	vmovaps	0x20(%%r13),%%ymm15	\n\t"\
		"vmulpd	(%%r14),%%ymm4,%%ymm4			\n\t	vsubpd	%%ymm14,%%ymm8 ,%%ymm8 		\n\t"\
		"vmulpd	(%%r14),%%ymm5,%%ymm5			\n\t	vsubpd	%%ymm15,%%ymm9 ,%%ymm9 		\n\t"\
		"vaddpd	%%ymm2,%%ymm4,%%ymm4			\n\t	vmulpd	(%%r14),%%ymm14,%%ymm14		\n\t"\
		"vaddpd	%%ymm3,%%ymm5,%%ymm5			\n\t	vmulpd	(%%r14),%%ymm15,%%ymm15		\n\t"\
		"movq	%[__cc1],%%rax					\n\t	vaddpd	%%ymm8 ,%%ymm14,%%ymm14		\n\t"\
		"vsubpd	%%ymm4,%%ymm6,%%ymm6			\n\t	vaddpd	%%ymm9 ,%%ymm15,%%ymm15		\n\t"\
		"vsubpd	%%ymm5,%%ymm7,%%ymm7			\n\t	vsubpd	%%ymm12,%%ymm10,%%ymm10		\n\t"\
		"vmulpd	(%%r14),%%ymm4,%%ymm4			\n\t	vsubpd	%%ymm13,%%ymm11,%%ymm11		\n\t"\
		"vmulpd	(%%r14),%%ymm5,%%ymm5			\n\t	movq	%[__cc1],%%r10		\n\t"\
		"vaddpd	%%ymm6,%%ymm4,%%ymm4			\n\t	vmulpd	(%%r14),%%ymm12,%%ymm12		\n\t"\
		"vaddpd	%%ymm7,%%ymm5,%%ymm5			\n\t	vmulpd	(%%r14),%%ymm13,%%ymm13		\n\t"\
		"vaddpd	    (%%rsi),%%ymm4,%%ymm4		\n\t	vaddpd	%%ymm10,%%ymm12,%%ymm12		\n\t"\
		"vaddpd	0x20(%%rsi),%%ymm5,%%ymm5		\n\t	vaddpd	%%ymm11,%%ymm13,%%ymm13		\n\t"\
		"vmovaps	%%ymm4,    (%%rdi)			\n\t	vsubpd	%%ymm12,%%ymm14,%%ymm14		\n\t"\
		"vmovaps	%%ymm5,0x20(%%rdi)			\n\t	vsubpd	%%ymm13,%%ymm15,%%ymm15		\n\t"\
		"vmulpd	0x20(%%rax),%%ymm6,%%ymm6		\n\t	vaddpd	%%ymm12,%%ymm12,%%ymm12		\n\t"\
		"vmulpd	0x20(%%rax),%%ymm7,%%ymm7		\n\t	vaddpd	%%ymm13,%%ymm13,%%ymm13		\n\t"\
		"vsubpd		(%%rsi),%%ymm4,%%ymm4		\n\t	vaddpd	%%ymm14,%%ymm12,%%ymm12		\n\t"\
		"vsubpd	0x20(%%rsi),%%ymm5,%%ymm5		\n\t	vaddpd	%%ymm15,%%ymm13,%%ymm13		\n\t"\
		"vmulpd		(%%rax),%%ymm4,%%ymm4		\n\t	vaddpd	    (%%r9 ),%%ymm12,%%ymm12	\n\t"\
		"vmulpd		(%%rax),%%ymm5,%%ymm5		\n\t	vaddpd	0x20(%%r9 ),%%ymm13,%%ymm13	\n\t"\
		"vaddpd		(%%rdi),%%ymm4,%%ymm4		\n\t	vmovaps	%%ymm12,    (%%r8 )	\n\t"\
		"vaddpd	0x20(%%rdi),%%ymm5,%%ymm5		\n\t	vmovaps	%%ymm13,0x20(%%r8 )	\n\t"\
		"vsubpd	%%ymm6,%%ymm4,%%ymm4			\n\t	vmulpd	0x20(%%r10),%%ymm14,%%ymm14	\n\t"\
		"vsubpd	%%ymm7,%%ymm5,%%ymm5			\n\t	vmulpd	0x20(%%r10),%%ymm15,%%ymm15	\n\t"\
		"vmulpd	(%%r14),%%ymm6,%%ymm6			\n\t	vsubpd	    (%%r9 ),%%ymm12,%%ymm12	\n\t"\
		"vmulpd	(%%r14),%%ymm7,%%ymm7			\n\t	vsubpd	0x20(%%r9 ),%%ymm13,%%ymm13	\n\t"\
		"vaddpd	%%ymm4,%%ymm6,%%ymm6			\n\t	vmulpd	    (%%r10),%%ymm12,%%ymm12	\n\t"\
		"vaddpd	%%ymm5,%%ymm7,%%ymm7			\n\t	vmulpd	    (%%r10),%%ymm13,%%ymm13	\n\t"\
		"vmovaps	%%ymm4,    (%%rsi)			\n\t	vaddpd	    (%%r8 ),%%ymm12,%%ymm12	\n\t"\
		"vmovaps	%%ymm5,0x20(%%rsi)			\n\t	vaddpd	0x20(%%r8 ),%%ymm13,%%ymm13	\n\t"\
		"vmovaps	%%ymm0,%%ymm4				\n\t	vsubpd	%%ymm14,%%ymm12,%%ymm12		\n\t"\
		"vmovaps	%%ymm1,%%ymm5				\n\t	vsubpd	%%ymm15,%%ymm13,%%ymm13		\n\t"\
		"vsubpd	%%ymm2,%%ymm0,%%ymm0			\n\t	vmulpd	(%%r14),%%ymm14,%%ymm14		\n\t"\
		"vsubpd	%%ymm3,%%ymm1,%%ymm1			\n\t	vmulpd	(%%r14),%%ymm15,%%ymm15		\n\t"\
		"vmulpd	0x40(%%rax),%%ymm0,%%ymm0		\n\t	vaddpd	%%ymm12,%%ymm14,%%ymm14		\n\t"\
		"vmulpd	0x40(%%rax),%%ymm1,%%ymm1		\n\t	vaddpd	%%ymm13,%%ymm15,%%ymm15		\n\t"\
		"vmulpd	0x60(%%rax),%%ymm2,%%ymm2		\n\t	vmovaps	%%ymm12,    (%%r9 )	\n\t"\
		"vmulpd	0x60(%%rax),%%ymm3,%%ymm3		\n\t	vmovaps	%%ymm13,0x20(%%r9 )	\n\t"\
		"vmulpd	0x80(%%rax),%%ymm4,%%ymm4		\n\t	vmovaps	%%ymm8 ,%%ymm12		\n\t"\
		"vmulpd	0x80(%%rax),%%ymm5,%%ymm5		\n\t	vmovaps	%%ymm9 ,%%ymm13		\n\t"\
		"vaddpd	%%ymm0,%%ymm2,%%ymm2			\n\t	vsubpd	%%ymm10,%%ymm8 ,%%ymm8 		\n\t"\
		"vaddpd	%%ymm1,%%ymm3,%%ymm3			\n\t	vsubpd	%%ymm11,%%ymm9 ,%%ymm9 		\n\t"\
		"vsubpd	%%ymm4,%%ymm0,%%ymm0			\n\t	vmulpd	0x40(%%r10),%%ymm8 ,%%ymm8 	\n\t"\
		"vsubpd	%%ymm5,%%ymm1,%%ymm1			\n\t	vmulpd	0x40(%%r10),%%ymm9 ,%%ymm9 	\n\t"\
		"vmovaps	    (%%rsi),%%ymm4			\n\t	vmulpd	0x60(%%r10),%%ymm10,%%ymm10	\n\t"\
		"vmovaps	0x20(%%rsi),%%ymm5			\n\t	vmulpd	0x60(%%r10),%%ymm11,%%ymm11	\n\t"\
		"vsubpd	%%ymm3,%%ymm6,%%ymm6			\n\t	vmulpd	0x80(%%r10),%%ymm12,%%ymm12	\n\t"\
		"vsubpd	%%ymm2,%%ymm7,%%ymm7			\n\t	vmulpd	0x80(%%r10),%%ymm13,%%ymm13	\n\t"\
		"vaddpd	%%ymm3,%%ymm3,%%ymm3			\n\t	vaddpd	%%ymm8 ,%%ymm10,%%ymm10		\n\t"\
		"vaddpd	%%ymm2,%%ymm2,%%ymm2			\n\t	vaddpd	%%ymm9 ,%%ymm11,%%ymm11		\n\t"\
		"vmovaps	%%ymm6,0x040(%%rdi)			\n\t	vsubpd	%%ymm12,%%ymm8 ,%%ymm8 		\n\t"\
		"vmovaps	%%ymm7,0x0e0(%%rdi)			\n\t	vsubpd	%%ymm13,%%ymm9 ,%%ymm9 		\n\t"\
		"vaddpd	%%ymm6,%%ymm3,%%ymm3			\n\t	vmovaps	    (%%r9 ),%%ymm12	\n\t"\
		"vaddpd	%%ymm7,%%ymm2,%%ymm2			\n\t	vmovaps	0x20(%%r9 ),%%ymm13	\n\t"\
		"vmovaps	%%ymm3,0x0c0(%%rdi)			\n\t	vsubpd	%%ymm11,%%ymm14,%%ymm14		\n\t"\
		"vmovaps	%%ymm2,0x060(%%rdi)			\n\t	vsubpd	%%ymm10,%%ymm15,%%ymm15		\n\t"\
		"vsubpd	%%ymm1,%%ymm4,%%ymm4			\n\t	vaddpd	%%ymm11,%%ymm11,%%ymm11		\n\t"\
		"vsubpd	%%ymm0,%%ymm5,%%ymm5			\n\t	vaddpd	%%ymm10,%%ymm10,%%ymm10		\n\t"\
		"vaddpd	%%ymm1,%%ymm1,%%ymm1			\n\t	vmovaps	%%ymm14,0x040(%%r8 )\n\t"\
		"vaddpd	%%ymm0,%%ymm0,%%ymm0			\n\t	vmovaps	%%ymm15,0x0e0(%%r8 )\n\t"\
		"vmovaps	%%ymm4,0x080(%%rdi)			\n\t	vaddpd	%%ymm14,%%ymm11,%%ymm11		\n\t"\
		"vmovaps	%%ymm5,0x120(%%rdi)			\n\t	vaddpd	%%ymm15,%%ymm10,%%ymm10		\n\t"\
		"vaddpd	%%ymm4,%%ymm1,%%ymm1			\n\t	vmovaps	%%ymm11,0x0c0(%%r8 )\n\t"\
		"vaddpd	%%ymm5,%%ymm0,%%ymm0			\n\t	vmovaps	%%ymm10,0x060(%%r8 )\n\t"\
		"vmovaps	%%ymm1,0x100(%%rdi)			\n\t	vsubpd	%%ymm9 ,%%ymm12,%%ymm12		\n\t"\
		"vmovaps	%%ymm0,0x0a0(%%rdi)			\n\t	vsubpd	%%ymm8 ,%%ymm13,%%ymm13		\n\t"\
		"										\n\t	vaddpd	%%ymm9 ,%%ymm9 ,%%ymm9 		\n\t"\
		"movq	%[__add0],%%rax					\n\t	vaddpd	%%ymm8 ,%%ymm8 ,%%ymm8 		\n\t"\
		"movq	%[__add0],%%rdx					\n\t	vmovaps	%%ymm12,0x080(%%r8 )\n\t"\
		"movslq	%[__p01],%%rsi					\n\t	vmovaps	%%ymm13,0x120(%%r8 )\n\t"\
		"shlq	$3,%%rsi						\n\t	vaddpd	%%ymm12,%%ymm9 ,%%ymm9 		\n\t"\
		"addq	%%rsi,%%rdx						\n\t	vaddpd	%%ymm13,%%ymm8 ,%%ymm8 		\n\t"\
		"movq	%%rdx,%%rbx						\n\t	vmovaps	%%ymm9 ,0x100(%%r8 )\n\t"\
		"addq	%%rsi,%%rdx						\n\t	vmovaps	%%ymm8 ,0x0a0(%%r8 )\n\t"\
		"movq	%%rdx,%%rcx						\n\t	\n\t"\
		"addq	%%rsi,%%rdx						\n\t	movq	%[__r00],%%rdi		\n\t"\
		"movslq	%[__p04],%%r8 					\n\t	movslq	%[__p16],%%r10		\n\t"\
		"shlq	$3,%%r8 						\n\t	shlq	$3,%%r10			\n\t"\
		"/*	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_C(r00,abdc) */\n\t"\
		"movq	$0x140,%%rsi					\n\t	/*	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_C(r02,abdc) */\n\t"\
		"addq	%%rdi,%%rsi						\n\t	movq	%%r10,%%r11				\n\t"\
		"vmovaps	     (%%rdi),%%ymm4			\n\t	movq	%%r10,%%r12				\n\t"\
		"vmovaps	     (%%rsi),%%ymm6			\n\t	movq	%%r10,%%r13				\n\t"\
		"vmovaps	0x020(%%rdi),%%ymm5			\n\t	addq	%%rax,%%r10		/* &a[j1+p16] */\n\t"\
		"vmovaps	0x020(%%rsi),%%ymm7			\n\t	addq	%%rbx,%%r11				\n\t"\
		"vmovaps	0x280(%%rdi),%%ymm0			\n\t	addq	%%rcx,%%r12				\n\t"\
		"vmovaps	0x280(%%rsi),%%ymm2			\n\t	addq	%%rdx,%%r13				\n\t"\
		"vmovaps	0x2a0(%%rdi),%%ymm1			\n\t	addq	$0x40,%%rdi		/* r02 */	\n\t"\
		"vmovaps	0x2a0(%%rsi),%%ymm3			\n\t	addq	$0x40,%%rsi				\n\t"\
		"vsubpd	%%ymm0,%%ymm4,%%ymm4			\n\t	vmovaps	     (%%rdi),%%ymm12	\n\t"\
		"vsubpd	%%ymm2,%%ymm6,%%ymm6			\n\t	vmovaps	     (%%rsi),%%ymm14	\n\t"\
		"vsubpd	%%ymm1,%%ymm5,%%ymm5			\n\t	vmovaps	0x020(%%rdi),%%ymm13	\n\t"\
		"vsubpd	%%ymm3,%%ymm7,%%ymm7			\n\t	vmovaps	0x020(%%rsi),%%ymm15	\n\t"\
		"vaddpd	%%ymm0,%%ymm0,%%ymm0			\n\t	vmovaps	0x280(%%rdi),%%ymm8 	\n\t"\
		"vaddpd	%%ymm2,%%ymm2,%%ymm2			\n\t	vmovaps	0x280(%%rsi),%%ymm10	\n\t"\
		"vaddpd	%%ymm1,%%ymm1,%%ymm1			\n\t	vmovaps	0x2a0(%%rdi),%%ymm9 	\n\t"\
		"vaddpd	%%ymm3,%%ymm3,%%ymm3			\n\t	vmovaps	0x2a0(%%rsi),%%ymm11	\n\t"\
		"vaddpd	%%ymm4,%%ymm0,%%ymm0			\n\t	vsubpd	%%ymm8 ,%%ymm12,%%ymm12			\n\t"\
		"vaddpd	%%ymm6,%%ymm2,%%ymm2			\n\t	vsubpd	%%ymm10,%%ymm14,%%ymm14			\n\t"\
		"vaddpd	%%ymm5,%%ymm1,%%ymm1			\n\t	vsubpd	%%ymm9 ,%%ymm13,%%ymm13			\n\t"\
		"vaddpd	%%ymm7,%%ymm3,%%ymm3			\n\t	vsubpd	%%ymm11,%%ymm15,%%ymm15			\n\t"\
		"/* Finish radix-4 butterfly: */		\n\t	vmulpd	(%%r14),%%ymm8 ,%%ymm8 			\n\t"\
		"vsubpd	%%ymm2,%%ymm0,%%ymm0			\n\t	vmulpd	(%%r14),%%ymm10,%%ymm10			\n\t"\
		"vsubpd	%%ymm7,%%ymm4,%%ymm4			\n\t	vmulpd	(%%r14),%%ymm9 ,%%ymm9 			\n\t"\
		"vsubpd	%%ymm3,%%ymm1,%%ymm1			\n\t	vmulpd	(%%r14),%%ymm11,%%ymm11			\n\t"\
		"vsubpd	%%ymm6,%%ymm5,%%ymm5			\n\t	vaddpd	%%ymm12,%%ymm8 ,%%ymm8 			\n\t"\
		"vmovaps	%%ymm0,     (%%rbx)			\n\t	vaddpd	%%ymm14,%%ymm10,%%ymm10			\n\t"\
		"vmovaps	%%ymm4,     (%%rdx)			\n\t	vaddpd	%%ymm13,%%ymm9 ,%%ymm9 			\n\t"\
		"vmovaps	%%ymm1,0x020(%%rbx)			\n\t	vaddpd	%%ymm15,%%ymm11,%%ymm11			\n\t"\
		"vmovaps	%%ymm5,0x020(%%rcx)			\n\t	/* Finish radix-4 butterfly and store results into main-array slots: */\n\t"\
		"vmulpd	(%%r14),%%ymm2,%%ymm2			\n\t	vsubpd	%%ymm10,%%ymm8 ,%%ymm8 			\n\t"\
		"vmulpd	(%%r14),%%ymm7,%%ymm7			\n\t	vsubpd	%%ymm15,%%ymm12,%%ymm12			\n\t"\
		"vmulpd	(%%r14),%%ymm3,%%ymm3			\n\t	vsubpd	%%ymm11,%%ymm9 ,%%ymm9 			\n\t"\
		"vmulpd	(%%r14),%%ymm6,%%ymm6			\n\t	vsubpd	%%ymm14,%%ymm13,%%ymm13			\n\t"\
		"vaddpd	%%ymm0,%%ymm2,%%ymm2			\n\t	vmovaps	%%ymm8 ,     (%%r11)	\n\t"\
		"vaddpd	%%ymm4,%%ymm7,%%ymm7			\n\t	vmovaps	%%ymm12,     (%%r13)	\n\t"\
		"vaddpd	%%ymm1,%%ymm3,%%ymm3			\n\t	vmovaps	%%ymm9 ,0x020(%%r11)	\n\t"\
		"vaddpd	%%ymm5,%%ymm6,%%ymm6			\n\t	vmovaps	%%ymm13,0x020(%%r12)	\n\t"\
		"vmovaps	%%ymm2,     (%%rax)			\n\t	vaddpd	%%ymm10,%%ymm10,%%ymm10			\n\t"\
		"vmovaps	%%ymm7,     (%%rcx)			\n\t	vaddpd	%%ymm15,%%ymm15,%%ymm15			\n\t"\
		"vmovaps	%%ymm3,0x020(%%rax)			\n\t	vaddpd	%%ymm11,%%ymm11,%%ymm11			\n\t"\
		"vmovaps	%%ymm6,0x020(%%rdx)			\n\t	vaddpd	%%ymm14,%%ymm14,%%ymm14			\n\t"\
		"										\n\t	vaddpd	%%ymm8 ,%%ymm10,%%ymm10			\n\t"\
		"/*	add0,1,2,3 = &a[j1+p04]+p3,2,1,0 */	\n\t	vaddpd	%%ymm12,%%ymm15,%%ymm15			\n\t"\
		"addq	%%r8 ,%%rax		/* &a[j1+p04] */\n\t	vaddpd	%%ymm9 ,%%ymm11,%%ymm11			\n\t"\
		"addq	%%r8 ,%%rbx						\n\t	vaddpd	%%ymm13,%%ymm14,%%ymm14			\n\t"\
		"addq	%%r8 ,%%rcx						\n\t	vmovaps	%%ymm10,     (%%r10)	\n\t"\
		"addq	%%r8 ,%%rdx						\n\t	vmovaps	%%ymm15,     (%%r12)	\n\t"\
		"/*	RADIX4_DIF_0TWIDDLE_STRIDE_C(r06) */\n\t	vmovaps	%%ymm11,0x020(%%r10)	\n\t"\
		"addq	$0x80,%%rdi	/* r06 */			\n\t	vmovaps	%%ymm14,0x020(%%r13)	\n\t"\
		"addq	$0x80,%%rsi						\n\t	\n\t"\
		"vmovaps	     (%%rdi),%%ymm4			\n\t	/* add0,1,2,3 = &a[j1+p12]+p2,3,0,1 */\n\t"\
		"vmovaps	     (%%rsi),%%ymm6			\n\t	subq	%%r8 ,%%r10		/* &a[j1+p12] */\n\t"\
		"vmovaps	0x020(%%rdi),%%ymm5			\n\t	subq	%%r8 ,%%r11				\n\t"\
		"vmovaps	0x020(%%rsi),%%ymm7			\n\t	subq	%%r8 ,%%r12				\n\t"\
		"vmovaps	0x280(%%rdi),%%ymm0			\n\t	subq	%%r8 ,%%r13				\n\t"\
		"vmovaps	0x280(%%rsi),%%ymm2			\n\t	/*	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_C(r04,cdab) */\n\t"\
		"vmovaps	0x2a0(%%rdi),%%ymm1			\n\t	subq	$0x40,%%rdi		/* r04 */	\n\t"\
		"vmovaps	0x2a0(%%rsi),%%ymm3			\n\t	subq	$0x40,%%rsi				\n\t"\
		"vsubpd	%%ymm0,%%ymm4,%%ymm4			\n\t	vmovaps	     (%%rdi),%%ymm12	\n\t"\
		"vsubpd	%%ymm2,%%ymm6,%%ymm6			\n\t	vmovaps	     (%%rsi),%%ymm14	\n\t"\
		"vsubpd	%%ymm1,%%ymm5,%%ymm5			\n\t	vmovaps	0x020(%%rdi),%%ymm13	\n\t"\
		"vsubpd	%%ymm3,%%ymm7,%%ymm7			\n\t	vmovaps	0x020(%%rsi),%%ymm15	\n\t"\
		"vaddpd	%%ymm0,%%ymm0,%%ymm0			\n\t	vmovaps	0x280(%%rdi),%%ymm8 	\n\t"\
		"vaddpd	%%ymm2,%%ymm2,%%ymm2			\n\t	vmovaps	0x280(%%rsi),%%ymm10	\n\t"\
		"vaddpd	%%ymm1,%%ymm1,%%ymm1			\n\t	vmovaps	0x2a0(%%rdi),%%ymm9 	\n\t"\
		"vaddpd	%%ymm3,%%ymm3,%%ymm3			\n\t	vmovaps	0x2a0(%%rsi),%%ymm11	\n\t"\
		"vaddpd	%%ymm4,%%ymm0,%%ymm0			\n\t	vsubpd	%%ymm8 ,%%ymm12,%%ymm12			\n\t"\
		"vaddpd	%%ymm6,%%ymm2,%%ymm2			\n\t	vsubpd	%%ymm10,%%ymm14,%%ymm14			\n\t"\
		"vaddpd	%%ymm5,%%ymm1,%%ymm1			\n\t	vsubpd	%%ymm9 ,%%ymm13,%%ymm13			\n\t"\
		"vaddpd	%%ymm7,%%ymm3,%%ymm3			\n\t	vsubpd	%%ymm11,%%ymm15,%%ymm15			\n\t"\
		"/* Finish radix-4 butterfly: */		\n\t	vmulpd	(%%r14),%%ymm8 ,%%ymm8 			\n\t"\
		"vsubpd	%%ymm2,%%ymm0,%%ymm0			\n\t	vmulpd	(%%r14),%%ymm10,%%ymm10			\n\t"\
		"vsubpd	%%ymm7,%%ymm4,%%ymm4			\n\t	vmulpd	(%%r14),%%ymm9 ,%%ymm9 			\n\t"\
		"vsubpd	%%ymm3,%%ymm1,%%ymm1			\n\t	vmulpd	(%%r14),%%ymm11,%%ymm11			\n\t"\
		"vsubpd	%%ymm6,%%ymm5,%%ymm5			\n\t	vaddpd	%%ymm12,%%ymm8 ,%%ymm8 			\n\t"\
		"vmovaps	%%ymm0,     (%%rcx)			\n\t	vaddpd	%%ymm14,%%ymm10,%%ymm10			\n\t"\
		"vmovaps	%%ymm4,     (%%rbx)			\n\t	vaddpd	%%ymm13,%%ymm9 ,%%ymm9 			\n\t"\
		"vmovaps	%%ymm1,0x020(%%rcx)			\n\t	vaddpd	%%ymm15,%%ymm11,%%ymm11			\n\t"\
		"vmovaps	%%ymm5,0x020(%%rax)			\n\t	/* Finish radix-4 butterfly and store results into main-array slots: */\n\t"\
		"vmulpd	(%%r14),%%ymm2,%%ymm2			\n\t	vsubpd	%%ymm10,%%ymm8 ,%%ymm8 			\n\t"\
		"vmulpd	(%%r14),%%ymm7,%%ymm7			\n\t	vsubpd	%%ymm15,%%ymm12,%%ymm12			\n\t"\
		"vmulpd	(%%r14),%%ymm3,%%ymm3			\n\t	vsubpd	%%ymm11,%%ymm9 ,%%ymm9 			\n\t"\
		"vmulpd	(%%r14),%%ymm6,%%ymm6			\n\t	vsubpd	%%ymm14,%%ymm13,%%ymm13			\n\t"\
		"vaddpd	%%ymm0,%%ymm2,%%ymm2			\n\t	vmovaps	%%ymm8 ,     (%%r13)	\n\t"\
		"vaddpd	%%ymm4,%%ymm7,%%ymm7			\n\t	vmovaps	%%ymm12,     (%%r10)	\n\t"\
		"vaddpd	%%ymm1,%%ymm3,%%ymm3			\n\t	vmovaps	%%ymm9 ,0x020(%%r13)	\n\t"\
		"vaddpd	%%ymm5,%%ymm6,%%ymm6			\n\t	vmovaps	%%ymm13,0x020(%%r11)	\n\t"\
		"vmovaps	%%ymm2,     (%%rdx)			\n\t	vaddpd	%%ymm10,%%ymm10,%%ymm10			\n\t"\
		"vmovaps	%%ymm7,     (%%rax)			\n\t	vaddpd	%%ymm15,%%ymm15,%%ymm15			\n\t"\
		"vmovaps	%%ymm3,0x020(%%rdx)			\n\t	vaddpd	%%ymm11,%%ymm11,%%ymm11			\n\t"\
		"vmovaps	%%ymm6,0x020(%%rbx)			\n\t	vaddpd	%%ymm14,%%ymm14,%%ymm14			\n\t"\
		"										\n\t	vaddpd	%%ymm8 ,%%ymm10,%%ymm10			\n\t"\
		"/*	add0,1,2,3 = &a[j1+p08]+p1,0,2,3 */	\n\t	vaddpd	%%ymm12,%%ymm15,%%ymm15			\n\t"\
		"addq	%%r8 ,%%rax		/* &a[j1+p08] */\n\t	vaddpd	%%ymm9 ,%%ymm11,%%ymm11			\n\t"\
		"addq	%%r8 ,%%rbx						\n\t	vaddpd	%%ymm13,%%ymm14,%%ymm14			\n\t"\
		"addq	%%r8 ,%%rcx						\n\t	vmovaps	%%ymm10,     (%%r12)	\n\t"\
		"addq	%%r8 ,%%rdx						\n\t	vmovaps	%%ymm15,     (%%r11)	\n\t"\
		"/*	RADIX4_DIF_0TWIDDLE_STRIDE_C(r08) */\n\t	vmovaps	%%ymm11,0x020(%%r12)	\n\t"\
		"addq	$0x80,%%rdi	/* r08 */			\n\t	vmovaps	%%ymm14,0x020(%%r10)	\n\t"\
		"addq	$0x80,%%rsi						\n\t"\
		"vmovaps	     (%%rdi),%%ymm4			\n\t"\
		"vmovaps	     (%%rsi),%%ymm6			\n\t"\
		"vmovaps	0x020(%%rdi),%%ymm5			\n\t"\
		"vmovaps	0x020(%%rsi),%%ymm7			\n\t"\
		"vmovaps	0x280(%%rdi),%%ymm0			\n\t"\
		"vmovaps	0x280(%%rsi),%%ymm2			\n\t"\
		"vmovaps	0x2a0(%%rdi),%%ymm1			\n\t"\
		"vmovaps	0x2a0(%%rsi),%%ymm3			\n\t"\
		"vsubpd	%%ymm0,%%ymm4,%%ymm4			\n\t"\
		"vsubpd	%%ymm2,%%ymm6,%%ymm6			\n\t"\
		"vsubpd	%%ymm1,%%ymm5,%%ymm5			\n\t"\
		"vsubpd	%%ymm3,%%ymm7,%%ymm7			\n\t"\
		"vaddpd	%%ymm0,%%ymm0,%%ymm0			\n\t"\
		"vaddpd	%%ymm2,%%ymm2,%%ymm2			\n\t"\
		"vaddpd	%%ymm1,%%ymm1,%%ymm1			\n\t"\
		"vaddpd	%%ymm3,%%ymm3,%%ymm3			\n\t"\
		"vaddpd	%%ymm4,%%ymm0,%%ymm0			\n\t"\
		"vaddpd	%%ymm6,%%ymm2,%%ymm2			\n\t"\
		"vaddpd	%%ymm5,%%ymm1,%%ymm1			\n\t"\
		"vaddpd	%%ymm7,%%ymm3,%%ymm3			\n\t"\
		"/* Finish radix-4 butterfly: */		\n\t"\
		"vsubpd	%%ymm2,%%ymm0,%%ymm0			\n\t"\
		"vsubpd	%%ymm7,%%ymm4,%%ymm4			\n\t"\
		"vsubpd	%%ymm3,%%ymm1,%%ymm1			\n\t"\
		"vsubpd	%%ymm6,%%ymm5,%%ymm5			\n\t"\
		"vmovaps	%%ymm0,     (%%rax)			\n\t"\
		"vmovaps	%%ymm4,     (%%rcx)			\n\t"\
		"vmovaps	%%ymm1,0x020(%%rax)			\n\t"\
		"vmovaps	%%ymm5,0x020(%%rdx)			\n\t"\
		"vaddpd	%%ymm2,%%ymm2,%%ymm2			\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7			\n\t"\
		"vaddpd	%%ymm3,%%ymm3,%%ymm3			\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6			\n\t"\
		"vaddpd	%%ymm0,%%ymm2,%%ymm2			\n\t"\
		"vaddpd	%%ymm4,%%ymm7,%%ymm7			\n\t"\
		"vaddpd	%%ymm1,%%ymm3,%%ymm3			\n\t"\
		"vaddpd	%%ymm5,%%ymm6,%%ymm6			\n\t"\
		"vmovaps	%%ymm2,     (%%rbx)			\n\t"\
		"vmovaps	%%ymm7,     (%%rdx)			\n\t"\
		"vmovaps	%%ymm3,0x020(%%rbx)			\n\t"\
		"vmovaps	%%ymm6,0x020(%%rcx)			\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__p01] "m" (Xp01)\
		 ,[__p04] "m" (Xp04)\
		 ,[__p08] "m" (Xp08)\
		 ,[__p16] "m" (Xp16)\
		 ,[__a0] "m" (Xa0)\
		 ,[__a1] "m" (Xa1)\
		 ,[__a2] "m" (Xa2)\
		 ,[__a3] "m" (Xa3)\
		 ,[__a4] "m" (Xa4)\
		 ,[__b0] "m" (Xb0)\
		 ,[__b1] "m" (Xb1)\
		 ,[__b2] "m" (Xb2)\
		 ,[__b3] "m" (Xb3)\
		 ,[__b4] "m" (Xb4)\
		 ,[__c0] "m" (Xc0)\
		 ,[__c1] "m" (Xc1)\
		 ,[__c2] "m" (Xc2)\
		 ,[__c3] "m" (Xc3)\
		 ,[__c4] "m" (Xc4)\
		 ,[__d0] "m" (Xd0)\
		 ,[__d1] "m" (Xd1)\
		 ,[__d2] "m" (Xd2)\
		 ,[__d3] "m" (Xd3)\
		 ,[__d4] "m" (Xd4)\
		 ,[__cc1] "m" (Xcc1)\
		 ,[__r00] "m" (Xr00)\
		 ,[__r10] "m" (Xr10)\
		 ,[__r20] "m" (Xr20)\
		 ,[__r30] "m" (Xr30)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r8","r9","r10","r11","r12","r13","r14","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
	);\
	}

  #else

	#define	SSE2_RADIX20_DIT_NOTWIDDLE(Xadd0,Xp01,Xp04,Xr00,Xr10,Xr20,Xr30,Xcc1,Xa0,Xa1,Xa2,Xa3,Xa4,Xb0,Xb1,Xb2,Xb3,Xb4,Xc0,Xc1,Xc2,Xc3,Xc4,Xd0,Xd1,Xd2,Xd3,Xd4)\
	{\
	__asm__ volatile (\
		"/* Outputs in SSE2 modes are temps 2*5*16 = 10*16 = 0x0a0 bytes apart: */\n\t"\
		"/*	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_C(abcd,r00)	*/		/*	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_C(dcba,r02)	*/\n\t"\
		"movq	%[__add0],%%rax	/* Use eax as base address */	\n\t	movq	$0x0a0,%%rdi		\n\t"\
		"/*	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_C(abcd,r00)	*/		\n\t	movslq	%[__p04],%%r8 		\n\t"\
		"movq	%%rax	,%%rdx									\n\t	shlq	$3,%%r8 			\n\t"\
		"movslq	%[__p01],%%rsi									\n\t	movq	$0x500,%%r9 		\n\t"\
		"shlq	$3,%%rsi										\n\t	movq	%%r8 ,%%r10			\n\t"\
		"addq	%%rsi,%%rdx										\n\t	movq	%%r8 ,%%r11			\n\t"\
		"movq	%%rdx,%%rbx		/* add1 = add0+p01 */			\n\t	movq	%%r8 ,%%r12			\n\t"\
		"addq	%%rsi,%%rdx										\n\t	movq	%%r8 ,%%r13			\n\t"\
		"movq	%%rdx,%%rcx		/* add3 = add0+p02 */			\n\t	addq	%%rax,%%r10		/* &a[j1+p04] */\n\t"\
		"addq	%%rsi,%%rdx		/* add2 = add0+p03 */			\n\t	addq	%%rbx,%%r11			\n\t"\
		"movaps	0x10(%%rax),%%xmm3								\n\t	movq	%[__r00],%%rsi		\n\t"\
		"movaps	    (%%rax),%%xmm2								\n\t	addq	%%rcx,%%r12			\n\t"\
		"movaps	    (%%rbx),%%xmm0								\n\t	addq	%%rdx,%%r13			\n\t"\
		"movaps	0x10(%%rbx),%%xmm1								\n\t	addq	%%rsi,%%r9	/* two */	\n\t"\
		"movaps	    (%%rcx),%%xmm6								\n\t	movaps	0x10(%%r13),%%xmm11	\n\t"\
		"movaps	0x10(%%rcx),%%xmm7								\n\t	movaps	    (%%r13),%%xmm10	\n\t"\
		"movaps	0x10(%%rdx),%%xmm5								\n\t	movaps	    (%%r12),%%xmm8 	\n\t"\
		"movaps	    (%%rdx),%%xmm4								\n\t	movaps	0x10(%%r12),%%xmm9 	\n\t"\
		"subpd	%%xmm0,%%xmm2									\n\t	movaps	    (%%r11),%%xmm14	\n\t"\
		"subpd	%%xmm4,%%xmm6									\n\t	movaps	0x10(%%r11),%%xmm15	\n\t"\
		"subpd	%%xmm1,%%xmm3									\n\t	movaps	0x10(%%r10),%%xmm13	\n\t"\
		"subpd	%%xmm5,%%xmm7									\n\t	movaps	    (%%r10),%%xmm12	\n\t"\
		"mulpd	(%%r9),%%xmm0									\n\t	subpd	%%xmm8 ,%%xmm10		\n\t"\
		"mulpd	(%%r9),%%xmm4									\n\t	subpd	%%xmm12,%%xmm14		\n\t"\
		"mulpd	(%%r9),%%xmm1									\n\t	subpd	%%xmm9 ,%%xmm11		\n\t"\
		"mulpd	(%%r9),%%xmm5									\n\t	subpd	%%xmm13,%%xmm15		\n\t"\
		"addpd	%%xmm2,%%xmm0									\n\t	mulpd	(%%r9),%%xmm8 		\n\t"\
		"addpd	%%xmm6,%%xmm4									\n\t	mulpd	(%%r9),%%xmm12		\n\t"\
		"addpd	%%xmm3,%%xmm1									\n\t	mulpd	(%%r9),%%xmm9 		\n\t"\
		"addpd	%%xmm7,%%xmm5									\n\t	mulpd	(%%r9),%%xmm13		\n\t"\
		"/* Finish radix-4 butterfly: */						\n\t	/* Finish radix-4 butterfly: */	addq	%%rsi,%%rdi	\n\t"\
		"subpd	%%xmm4,%%xmm0									\n\t	addpd	%%xmm10,%%xmm8 		\n\t"\
		"subpd	%%xmm7,%%xmm2									\n\t	addpd	%%xmm14,%%xmm12		\n\t"\
		"subpd	%%xmm5,%%xmm1									\n\t	addpd	%%xmm11,%%xmm9 		\n\t"\
		"subpd	%%xmm6,%%xmm3									\n\t	addpd	%%xmm15,%%xmm13		\n\t"\
		"movaps	%%xmm0,0x140(%%rsi)								\n\t	subpd	%%xmm12,%%xmm8 		\n\t"\
		"movaps	%%xmm2,     (%%rdi)								\n\t	subpd	%%xmm15,%%xmm10		\n\t"\
		"movaps	%%xmm1,0x150(%%rsi)								\n\t	subpd	%%xmm13,%%xmm9 		\n\t"\
		"movaps	%%xmm3,0x150(%%rdi)								\n\t	subpd	%%xmm14,%%xmm11		\n\t"\
		"addpd	%%xmm4,%%xmm4									\n\t	movaps	%%xmm8 ,0x160(%%rsi)	\n\t"\
		"addpd	%%xmm7,%%xmm7									\n\t	movaps	%%xmm10,0x160(%%rdi)	\n\t"\
		"addpd	%%xmm5,%%xmm5									\n\t	movaps	%%xmm9 ,0x170(%%rsi)	\n\t"\
		"addpd	%%xmm6,%%xmm6									\n\t	movaps	%%xmm11,0x030(%%rdi)	\n\t"\
		"addpd	%%xmm0,%%xmm4									\n\t	mulpd	(%%r9),%%xmm12		\n\t"\
		"addpd	%%xmm2,%%xmm7									\n\t	mulpd	(%%r9),%%xmm15		\n\t"\
		"addpd	%%xmm1,%%xmm5									\n\t	mulpd	(%%r9),%%xmm13		\n\t"\
		"addpd	%%xmm3,%%xmm6									\n\t	mulpd	(%%r9),%%xmm14		\n\t"\
		"movaps	%%xmm4,     (%%rsi)								\n\t	addpd	%%xmm8 ,%%xmm12		\n\t"\
		"movaps	%%xmm7,0x140(%%rdi)								\n\t	addpd	%%xmm10,%%xmm15		\n\t"\
		"movaps	%%xmm5,0x010(%%rsi)								\n\t	addpd	%%xmm9 ,%%xmm13		\n\t"\
		"movaps	%%xmm6,0x010(%%rdi)								\n\t	addpd	%%xmm11,%%xmm14		\n\t"\
		"addq	%%r8 ,%%r8 		/* p08 */	\n\t"\
		"addq	%%r8 ,%%rax		/* &a[j1+p08] */				\n\t	movaps	%%xmm12,0x020(%%rsi)	\n\t"\
		"addq	%%r8 ,%%rbx										\n\t	movaps	%%xmm15,0x020(%%rdi)	\n\t"\
		"addq	%%r8 ,%%rcx										\n\t	movaps	%%xmm13,0x030(%%rsi)	\n\t"\
		"addq	%%r8 ,%%rdx										\n\t	movaps	%%xmm14,0x170(%%rdi)	\n\t"\
		"/*	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_C(bacd,r04)	*/		\n\t	/*	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_C(cdab,r06)	*/\n\t"\
		"movaps	0x10(%%rbx),%%xmm3								\n\t	addq	%%r8 ,%%r10		/* &a[j1+p12] */\n\t"\
		"movaps	    (%%rbx),%%xmm2								\n\t	addq	%%r8 ,%%r11			\n\t"\
		"movaps	    (%%rax),%%xmm0								\n\t	addq	%%r8 ,%%r12			\n\t"\
		"movaps	0x10(%%rax),%%xmm1								\n\t	addq	%%r8 ,%%r13			\n\t"\
		"movaps	    (%%rcx),%%xmm6								\n\t	movaps	0x10(%%r12),%%xmm11	\n\t"\
		"movaps	0x10(%%rcx),%%xmm7								\n\t	movaps	    (%%r12),%%xmm10	\n\t"\
		"movaps	0x10(%%rdx),%%xmm5								\n\t	movaps	    (%%r13),%%xmm8 	\n\t"\
		"movaps	    (%%rdx),%%xmm4								\n\t	movaps	0x10(%%r13),%%xmm9 	\n\t"\
		"subpd	%%xmm0,%%xmm2									\n\t	movaps	    (%%r10),%%xmm14	\n\t"\
		"subpd	%%xmm4,%%xmm6									\n\t	movaps	0x10(%%r10),%%xmm15	\n\t"\
		"subpd	%%xmm1,%%xmm3									\n\t	movaps	0x10(%%r11),%%xmm13	\n\t"\
		"subpd	%%xmm5,%%xmm7									\n\t	movaps	    (%%r11),%%xmm12	\n\t"\
		"mulpd	(%%r9),%%xmm0									\n\t	subpd	%%xmm8 ,%%xmm10		\n\t"\
		"mulpd	(%%r9),%%xmm4									\n\t	subpd	%%xmm12,%%xmm14		\n\t"\
		"mulpd	(%%r9),%%xmm1									\n\t	subpd	%%xmm9 ,%%xmm11		\n\t"\
		"mulpd	(%%r9),%%xmm5									\n\t	subpd	%%xmm13,%%xmm15		\n\t"\
		"addpd	%%xmm2,%%xmm0									\n\t	mulpd	(%%r9),%%xmm8 		\n\t"\
		"addpd	%%xmm6,%%xmm4									\n\t	mulpd	(%%r9),%%xmm12		\n\t"\
		"addpd	%%xmm3,%%xmm1									\n\t	mulpd	(%%r9),%%xmm9 		\n\t"\
		"addpd	%%xmm7,%%xmm5									\n\t	mulpd	(%%r9),%%xmm13		\n\t"\
		"/* Finish radix-4 butterfly: */						\n\t	addpd	%%xmm10,%%xmm8 		\n\t"\
		"addq	$0x040,%%rsi	/* r04 */						\n\t	addpd	%%xmm14,%%xmm12		\n\t"\
		"subpd	%%xmm4,%%xmm0									\n\t	addpd	%%xmm11,%%xmm9 		\n\t"\
		"subpd	%%xmm7,%%xmm2									\n\t	addpd	%%xmm15,%%xmm13		\n\t"\
		"subpd	%%xmm5,%%xmm1									\n\t	/* Finish radix-4 butterfly: */\n\t"\
		"subpd	%%xmm6,%%xmm3									\n\t	addq	$0x040,%%rdi		\n\t"\
		"movaps	%%xmm0,0x140(%%rsi)								\n\t	subpd	%%xmm12,%%xmm8 		\n\t"\
		"movaps	%%xmm2,0x140(%%rdi)								\n\t	subpd	%%xmm15,%%xmm10		\n\t"\
		"movaps	%%xmm1,0x150(%%rsi)								\n\t	subpd	%%xmm13,%%xmm9 		\n\t"\
		"movaps	%%xmm3,0x010(%%rdi)								\n\t	subpd	%%xmm14,%%xmm11		\n\t"\
		"addpd	%%xmm4,%%xmm4									\n\t	movaps	%%xmm8 ,0x160(%%rsi)	\n\t"\
		"addpd	%%xmm7,%%xmm7									\n\t	movaps	%%xmm10,0x160(%%rdi)	\n\t"\
		"addpd	%%xmm5,%%xmm5									\n\t	movaps	%%xmm9 ,0x170(%%rsi)	\n\t"\
		"addpd	%%xmm6,%%xmm6									\n\t	movaps	%%xmm11,0x030(%%rdi)	\n\t"\
		"addpd	%%xmm0,%%xmm4									\n\t	mulpd	(%%r9),%%xmm12		\n\t"\
		"addpd	%%xmm2,%%xmm7									\n\t	mulpd	(%%r9),%%xmm15		\n\t"\
		"addpd	%%xmm1,%%xmm5									\n\t	mulpd	(%%r9),%%xmm13		\n\t"\
		"addpd	%%xmm3,%%xmm6									\n\t	mulpd	(%%r9),%%xmm14		\n\t"\
		"movaps	%%xmm4,     (%%rsi)								\n\t	addpd	%%xmm8 ,%%xmm12		\n\t"\
		"movaps	%%xmm7,     (%%rdi)								\n\t	addpd	%%xmm10,%%xmm15		\n\t"\
		"movaps	%%xmm5,0x010(%%rsi)								\n\t	addpd	%%xmm9 ,%%xmm13		\n\t"\
		"movaps	%%xmm6,0x150(%%rdi)								\n\t	addpd	%%xmm11,%%xmm14		\n\t"\
		"														\n\t	movaps	%%xmm12,0x020(%%rsi)	\n\t"\
		"														\n\t	movaps	%%xmm15,0x020(%%rdi)	\n\t"\
		"														\n\t	movaps	%%xmm13,0x030(%%rsi)	\n\t"\
		"														\n\t	movaps	%%xmm14,0x170(%%rdi)	\n\t"\
		"/*	add0,1,2,3 = &a[j1+p16]+p0,1,3,2 */\n\t"\
		"addq	%%r8 ,%%rax		/* &a[j1+p16] */\n\t"\
		"addq	%%r8 ,%%rbx			\n\t"\
		"addq	%%r8 ,%%rcx			\n\t"\
		"addq	%%r8 ,%%rdx			\n\t"\
		"/*	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_C(abdc,r08)	*/\n\t"\
		"movaps	0x10(%%rax),%%xmm3	\n\t"\
		"movaps	    (%%rax),%%xmm2	\n\t"\
		"movaps	    (%%rbx),%%xmm0	\n\t"\
		"movaps	0x10(%%rbx),%%xmm1	\n\t"\
		"movaps	    (%%rdx),%%xmm6	\n\t"\
		"movaps	0x10(%%rdx),%%xmm7	\n\t"\
		"movaps	0x10(%%rcx),%%xmm5	\n\t"\
		"movaps	    (%%rcx),%%xmm4	\n\t"\
		"addq	$0x040,%%rsi	/* r08 */\n\t"\
		"subpd	%%xmm0,%%xmm2		\n\t"\
		"subpd	%%xmm4,%%xmm6		\n\t"\
		"subpd	%%xmm1,%%xmm3		\n\t"\
		"subpd	%%xmm5,%%xmm7		\n\t"\
		"addpd	%%xmm0,%%xmm0		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm1,%%xmm1		\n\t"\
		"addpd	%%xmm5,%%xmm5									\n\t	movq	%%r9,%%r14	/* two */\n\t"\
		"addpd	%%xmm2,%%xmm0		\n\t"\
		"addpd	%%xmm6,%%xmm4		\n\t"\
		"addpd	%%xmm3,%%xmm1									\n\t/*** This needs to read r18 *after* the last radix-4 writes it! ***/\n\t"\
		"addpd	%%xmm7,%%xmm5									\n\t/* Radix-5 DFT uses adjacent temps, i.e. stride = 2*16 bytes: */\n\t"\
		"/* Finish radix-4 butterfly: */						\n\t/*	SSE2_RADIX_05_DFT_0TWIDDLE(r10,r12,r14,r16,r18,cc1,Xb0,Xb1,Xb2,Xb3,Xb4)	*/\n\t"\
		"addq	$0x040,%%rdi	/* r18 */						\n\t	movq	%[__r10],%%r8 		\n\t"\
		"subpd	%%xmm4,%%xmm0									\n\t	movq	%%r8 ,%%r10			\n\t"\
		"subpd	%%xmm7,%%xmm2									\n\t	movq	%%r8 ,%%r11			\n\t"\
		"subpd	%%xmm5,%%xmm1									\n\t	movq	%%r8 ,%%r12			\n\t"\
		"subpd	%%xmm6,%%xmm3									\n\t/*	movq	%%r8 ,%%r13		*/	\n\t"\
		"movaps	%%xmm0,0x140(%%rsi)	/* r28.r */					\n\t	addq	$0x020,%%r10		\n\t"\
		"movaps	%%xmm2,0x140(%%rdi)	/* r38.r */					\n\t	addq	$0x040,%%r11		\n\t"\
		"movaps	%%xmm1,0x150(%%rsi)	/* r28.i */					\n\t	addq	$0x060,%%r12		\n\t"\
		"/*vaps	%%xmm3,0x010(%%rdi)	** r18.i */					\n\t/*	addq	$0x080,%%r13	*/	\n\t"\
		"addpd	%%xmm4,%%xmm4									\n\t	movq	%[__b0],%%r9 		\n\t"\
		"addpd	%%xmm7,%%xmm7									\n\t	movaps	    (%%r10),%%xmm8 	\n\t"\
		"addpd	%%xmm5,%%xmm5									\n\t	movaps	0x10(%%r10),%%xmm9 	\n\t"\
		"addpd	%%xmm6,%%xmm6									\n\t	movaps	    (%%r11),%%xmm10	\n\t"\
		"addpd	%%xmm0,%%xmm4									\n\t	movaps	0x10(%%r11),%%xmm11	\n\t"\
		"addpd	%%xmm2,%%xmm7									\n\t	movaps	    (%%r12),%%xmm12	\n\t"\
		"addpd	%%xmm1,%%xmm5									\n\t	movaps	0x10(%%r12),%%xmm13	\n\t"\
		"addpd	%%xmm3,%%xmm6									\n\t/*	movaps	    (%%r13),%%xmm14	// r18.r */\n\t"\
		"/*vaps	%%xmm4,     (%%rsi)	// r08.r */					\n\t/*	movaps	0x10(%%r13),%%xmm15	// r18.i */\n\t"\
		"/*vaps	%%xmm7,     (%%rdi)	// r18.r */					\n\t	subpd	%%xmm7,%%xmm8 		\n\t"\
		"/*vaps	%%xmm5,0x010(%%rsi)	// r08.i */					\n\t	subpd	%%xmm3,%%xmm9 		\n\t"\
		"movaps	%%xmm6,0x150(%%rdi)	/* r38.i */					\n\t	addpd	%%xmm7,%%xmm7		\n\t"\
		"														\n\t	addpd	%%xmm3,%%xmm3		\n\t"\
		"/*	SSE2_RADIX_05_DFT_0TWIDDLE(r00...)	*/				\n\t	addpd	%%xmm8 ,%%xmm7		\n\t"\
		"movq	%[__r00],%%rsi									\n\t	addpd	%%xmm9 ,%%xmm3		\n\t"\
		"movq	%%rsi,%%rax										\n\t	subpd	%%xmm12,%%xmm10		\n\t	movaps	%%xmm7 ,%%xmm14		\n\t"\
		"movq	%%rsi,%%rbx										\n\t	subpd	%%xmm13,%%xmm11		\n\t	movaps	%%xmm3 ,%%xmm15		\n\t"\
		"movq	%%rsi,%%rcx		\n\t	movaps	%%xmm4 ,%%xmm6	\n\t	addpd	%%xmm12,%%xmm12		\n\t"\
		"/*vq	%%rsi,%%rdx	*/	\n\t	movaps	%%xmm5 ,%%xmm7	\n\t	addpd	%%xmm13,%%xmm13		\n\t"\
		"addq	$0x020,%%rax									\n\t	addpd	%%xmm10,%%xmm12		\n\t"\
		"addq	$0x040,%%rbx									\n\t	addpd	%%xmm11,%%xmm13		\n\t"\
		"addq	$0x060,%%rcx									\n\t	movq	%[__cc1],%%r10		\n\t"\
		"/*dq	$0x080,%%rdx	*/								\n\t	subpd	%%xmm12,%%xmm14		\n\t"\
		"movq	%[__a0],%%rdi									\n\t	subpd	%%xmm13,%%xmm15		\n\t"\
		"movaps	    (%%rax),%%xmm0								\n\t	addpd	%%xmm12,%%xmm12		\n\t"\
		"movaps	0x10(%%rax),%%xmm1								\n\t	addpd	%%xmm13,%%xmm13		\n\t"\
		"movaps	    (%%rbx),%%xmm2								\n\t	addpd	%%xmm14,%%xmm12		\n\t"\
		"movaps	0x10(%%rbx),%%xmm3								\n\t	addpd	%%xmm15,%%xmm13		\n\t"\
		"movaps	    (%%rcx),%%xmm4								\n\t	addpd	    (%%r8 ),%%xmm12	\n\t"\
		"movaps	0x10(%%rcx),%%xmm5								\n\t	addpd	0x10(%%r8 ),%%xmm13	\n\t"\
		"/*vaps	    (%%rdx),%%xmm6	// r08.r */					\n\t	movaps	%%xmm12,    (%%r9 )	\n\t"\
		"/*vaps	0x10(%%rdx),%%xmm7	// r08.i */					\n\t	movaps	%%xmm13,0x10(%%r9 )	\n\t"\
		"subpd	%%xmm6,%%xmm0									\n\t	mulpd	0x10(%%r10),%%xmm14	\n\t"\
		"subpd	%%xmm7,%%xmm1									\n\t	mulpd	0x10(%%r10),%%xmm15	\n\t"\
		"addpd	%%xmm6,%%xmm6									\n\t	subpd	    (%%r8 ),%%xmm12	\n\t"\
		"addpd	%%xmm7,%%xmm7									\n\t	subpd	0x10(%%r8 ),%%xmm13	\n\t"\
		"addpd	%%xmm0,%%xmm6									\n\t	mulpd	    (%%r10),%%xmm12	\n\t"\
		"addpd	%%xmm1,%%xmm7									\n\t	mulpd	    (%%r10),%%xmm13	\n\t"\
		"subpd	%%xmm4,%%xmm2									\n\t	addpd	    (%%r9 ),%%xmm12	\n\t"\
		"subpd	%%xmm5,%%xmm3									\n\t	addpd	0x10(%%r9 ),%%xmm13	\n\t"\
		"mulpd	(%%r14),%%xmm4									\n\t	subpd	%%xmm14,%%xmm12		\n\t"\
		"mulpd	(%%r14),%%xmm5									\n\t	subpd	%%xmm15,%%xmm13		\n\t"\
		"addpd	%%xmm2,%%xmm4									\n\t	mulpd	(%%r14),%%xmm14		\n\t"\
		"addpd	%%xmm3,%%xmm5									\n\t	mulpd	(%%r14),%%xmm15		\n\t"\
		"movq	%[__cc1],%%rax									\n\t	addpd	%%xmm12,%%xmm14		\n\t"\
		"subpd	%%xmm4,%%xmm6									\n\t	addpd	%%xmm13,%%xmm15		\n\t"\
		"subpd	%%xmm5,%%xmm7									\n\t	movaps	%%xmm12,    (%%r8 )	\n\t"\
		"addpd	%%xmm4,%%xmm4									\n\t	movaps	%%xmm13,0x10(%%r8 )	\n\t"\
		"addpd	%%xmm5,%%xmm5									\n\t	movaps	%%xmm8 ,%%xmm12		\n\t"\
		"addpd	%%xmm6,%%xmm4									\n\t	movaps	%%xmm9 ,%%xmm13		\n\t"\
		"addpd	%%xmm7,%%xmm5									\n\t	subpd	%%xmm10,%%xmm8 		\n\t"\
		"addpd	    (%%rsi),%%xmm4								\n\t	subpd	%%xmm11,%%xmm9 		\n\t"\
		"addpd	0x10(%%rsi),%%xmm5								\n\t	mulpd	0x20(%%r10),%%xmm8 	\n\t"\
		"movaps	%%xmm4,    (%%rdi)								\n\t	mulpd	0x20(%%r10),%%xmm9 	\n\t"\
		"movaps	%%xmm5,0x10(%%rdi)								\n\t	mulpd	0x30(%%r10),%%xmm10	\n\t"\
		"mulpd	0x10(%%rax),%%xmm6								\n\t	mulpd	0x30(%%r10),%%xmm11	\n\t"\
		"mulpd	0x10(%%rax),%%xmm7								\n\t	mulpd	0x40(%%r10),%%xmm12	\n\t"\
		"subpd	    (%%rsi),%%xmm4								\n\t	mulpd	0x40(%%r10),%%xmm13	\n\t"\
		"subpd	0x10(%%rsi),%%xmm5								\n\t	addpd	%%xmm8 ,%%xmm10		\n\t"\
		"mulpd	    (%%rax),%%xmm4								\n\t	addpd	%%xmm9 ,%%xmm11		\n\t"\
		"mulpd	    (%%rax),%%xmm5								\n\t	subpd	%%xmm12,%%xmm8 		\n\t"\
		"addpd	    (%%rdi),%%xmm4								\n\t	subpd	%%xmm13,%%xmm9 		\n\t"\
		"addpd	0x10(%%rdi),%%xmm5								\n\t	movaps	    (%%r8 ),%%xmm12	\n\t"\
		"subpd	%%xmm6,%%xmm4									\n\t	movaps	0x10(%%r8 ),%%xmm13	\n\t"\
		"subpd	%%xmm7,%%xmm5									\n\t	movq	%[__b1],%%r10		\n\t"\
		"mulpd	(%%r14),%%xmm6									\n\t	subpd	%%xmm11,%%xmm14		\n\t"\
		"mulpd	(%%r14),%%xmm7									\n\t	subpd	%%xmm10,%%xmm15		\n\t"\
		"addpd	%%xmm4,%%xmm6									\n\t	mulpd	(%%r14),%%xmm11		\n\t"\
		"addpd	%%xmm5,%%xmm7									\n\t	mulpd	(%%r14),%%xmm10		\n\t"\
		"movaps	%%xmm4,    (%%rsi)								\n\t	movq	%[__b4],%%r13		\n\t"\
		"movaps	%%xmm5,0x10(%%rsi)								\n\t	movaps	%%xmm14,    (%%r10)	\n\t"\
		"movaps	%%xmm0,%%xmm4									\n\t	movaps	%%xmm15,0x10(%%r13)	\n\t"\
		"movaps	%%xmm1,%%xmm5									\n\t	addpd	%%xmm14,%%xmm11		\n\t"\
		"subpd	%%xmm2,%%xmm0									\n\t	addpd	%%xmm15,%%xmm10		\n\t"\
		"subpd	%%xmm3,%%xmm1									\n\t	movaps	%%xmm11,    (%%r13)	\n\t"\
		"mulpd	0x20(%%rax),%%xmm0								\n\t	movaps	%%xmm10,0x10(%%r10)	\n\t"\
		"mulpd	0x20(%%rax),%%xmm1								\n\t	movq	%[__b2],%%r11		\n\t"\
		"mulpd	0x30(%%rax),%%xmm2								\n\t	movq	%[__b3],%%r12		\n\t"\
		"mulpd	0x30(%%rax),%%xmm3								\n\t								\n\t"\
		"mulpd	0x40(%%rax),%%xmm4								\n\t	subpd	%%xmm9 ,%%xmm12		\n\t"\
		"mulpd	0x40(%%rax),%%xmm5								\n\t	subpd	%%xmm8 ,%%xmm13		\n\t"\
		"addpd	%%xmm0,%%xmm2									\n\t	mulpd	(%%r14),%%xmm9 		\n\t"\
		"addpd	%%xmm1,%%xmm3									\n\t	mulpd	(%%r14),%%xmm8 		\n\t"\
		"subpd	%%xmm4,%%xmm0									\n\t	movaps	%%xmm12,    (%%r11)	\n\t"\
		"subpd	%%xmm5,%%xmm1									\n\t	movaps	%%xmm13,0x10(%%r12)	\n\t"\
		"movaps	    (%%rsi),%%xmm4								\n\t	addpd	%%xmm12,%%xmm9 		\n\t"\
		"movaps	0x10(%%rsi),%%xmm5								\n\t	addpd	%%xmm13,%%xmm8 		\n\t"\
		"movq	%[__a1],%%rax									\n\t	movaps	%%xmm9 ,    (%%r12)	\n\t"\
		"movq	%[__a4],%%rdx									\n\t	movaps	%%xmm8 ,0x10(%%r11)	\n\t"\
		"subpd	%%xmm3,%%xmm6									\n\t	/*	SSE2_RADIX_05_DFT_0TWIDDLE(r30)	*/\n\t"\
		"subpd	%%xmm2,%%xmm7									\n\t	movq	%[__r30],%%r8 		\n\t"\
		"addpd	%%xmm3,%%xmm3									\n\t	movq	%%r8 ,%%r10			\n\t"\
		"addpd	%%xmm2,%%xmm2									\n\t	movq	%%r8 ,%%r11			\n\t"\
		"movaps	%%xmm6,    (%%rax)								\n\t	movq	%%r8 ,%%r12			\n\t"\
		"movaps	%%xmm7,0x10(%%rdx)								\n\t	movq	%%r8 ,%%r13			\n\t"\
		"addpd	%%xmm6,%%xmm3									\n\t	addq	$0x020,%%r10		\n\t"\
		"addpd	%%xmm7,%%xmm2									\n\t	addq	$0x040,%%r11		\n\t"\
		"movaps	%%xmm3,    (%%rdx)								\n\t	addq	$0x060,%%r12		\n\t"\
		"movaps	%%xmm2,0x10(%%rax)								\n\t	addq	$0x080,%%r13		\n\t"\
		"movq	%[__a2],%%rbx									\n\t	movq	%[__d0],%%r9 		\n\t"\
		"movq	%[__a3],%%rcx									\n\t	movaps	    (%%r10),%%xmm8 	\n\t"\
		"subpd	%%xmm1,%%xmm4									\n\t	movaps	0x10(%%r10),%%xmm9 	\n\t"\
		"subpd	%%xmm0,%%xmm5									\n\t	movaps	    (%%r11),%%xmm10	\n\t"\
		"addpd	%%xmm1,%%xmm1									\n\t	movaps	0x10(%%r11),%%xmm11	\n\t"\
		"addpd	%%xmm0,%%xmm0									\n\t	movaps	    (%%r12),%%xmm12	\n\t"\
		"movaps	%%xmm4,    (%%rbx)								\n\t	movaps	0x10(%%r12),%%xmm13	\n\t"\
		"movaps	%%xmm5,0x10(%%rcx)								\n\t	movaps	    (%%r13),%%xmm14	\n\t"\
		"addpd	%%xmm4,%%xmm1									\n\t	movaps	0x10(%%r13),%%xmm15	\n\t"\
		"addpd	%%xmm5,%%xmm0									\n\t	subpd	%%xmm14,%%xmm8 		\n\t"\
		"movaps	%%xmm1,    (%%rcx)								\n\t	subpd	%%xmm15,%%xmm9 		\n\t"\
		"movaps	%%xmm0,0x10(%%rbx)								\n\t	addpd	%%xmm14,%%xmm14		\n\t"\
		"														\n\t	addpd	%%xmm15,%%xmm15		\n\t"\
		"/*	SSE2_RADIX_05_DFT_0TWIDDLE(r20)	*/					\n\t	addpd	%%xmm8 ,%%xmm14		\n\t"\
		"movq	%[__r20],%%rsi									\n\t	addpd	%%xmm9 ,%%xmm15		\n\t"\
		"movq	%%rsi,%%rax										\n\t	subpd	%%xmm12,%%xmm10		\n\t"\
		"movq	%%rsi,%%rbx										\n\t	subpd	%%xmm13,%%xmm11		\n\t"\
		"movq	%%rsi,%%rcx										\n\t	addpd	%%xmm12,%%xmm12		\n\t"\
		"movq	%%rsi,%%rdx										\n\t	addpd	%%xmm13,%%xmm13		\n\t"\
		"addq	$0x020,%%rax									\n\t	addpd	%%xmm10,%%xmm12		\n\t"\
		"addq	$0x040,%%rbx									\n\t	addpd	%%xmm11,%%xmm13		\n\t"\
		"addq	$0x060,%%rcx									\n\t	movq	%[__cc1],%%r10		\n\t"\
		"addq	$0x080,%%rdx									\n\t	subpd	%%xmm12,%%xmm14		\n\t"\
		"movq	%[__c0],%%rdi									\n\t	subpd	%%xmm13,%%xmm15		\n\t"\
		"movaps	    (%%rax),%%xmm0								\n\t	addpd	%%xmm12,%%xmm12		\n\t"\
		"movaps	0x10(%%rax),%%xmm1								\n\t	addpd	%%xmm13,%%xmm13		\n\t"\
		"movaps	    (%%rbx),%%xmm2								\n\t	addpd	%%xmm14,%%xmm12		\n\t"\
		"movaps	0x10(%%rbx),%%xmm3								\n\t	addpd	%%xmm15,%%xmm13		\n\t"\
		"movaps	    (%%rcx),%%xmm4								\n\t	addpd	    (%%r8 ),%%xmm12	\n\t"\
		"movaps	0x10(%%rcx),%%xmm5								\n\t	addpd	0x10(%%r8 ),%%xmm13	\n\t"\
		"movaps	    (%%rdx),%%xmm6								\n\t	movaps	%%xmm12,    (%%r9 )	\n\t"\
		"movaps	0x10(%%rdx),%%xmm7								\n\t	movaps	%%xmm13,0x10(%%r9 )	\n\t"\
		"subpd	%%xmm6,%%xmm0									\n\t	mulpd	0x10(%%r10),%%xmm14	\n\t"\
		"subpd	%%xmm7,%%xmm1									\n\t	mulpd	0x10(%%r10),%%xmm15	\n\t"\
		"mulpd	(%%r14),%%xmm6									\n\t	subpd	    (%%r8 ),%%xmm12	\n\t"\
		"mulpd	(%%r14),%%xmm7									\n\t	subpd	0x10(%%r8 ),%%xmm13	\n\t"\
		"addpd	%%xmm0,%%xmm6									\n\t	mulpd	    (%%r10),%%xmm12	\n\t"\
		"addpd	%%xmm1,%%xmm7									\n\t	mulpd	    (%%r10),%%xmm13	\n\t"\
		"subpd	%%xmm4,%%xmm2									\n\t	addpd	    (%%r9 ),%%xmm12	\n\t"\
		"subpd	%%xmm5,%%xmm3									\n\t	addpd	0x10(%%r9 ),%%xmm13	\n\t"\
		"mulpd	(%%r14),%%xmm4									\n\t	subpd	%%xmm14,%%xmm12		\n\t"\
		"mulpd	(%%r14),%%xmm5									\n\t	subpd	%%xmm15,%%xmm13		\n\t"\
		"addpd	%%xmm2,%%xmm4									\n\t	mulpd	(%%r14),%%xmm14		\n\t"\
		"addpd	%%xmm3,%%xmm5									\n\t	mulpd	(%%r14),%%xmm15		\n\t"\
		"movq	%[__cc1],%%rax									\n\t	addpd	%%xmm12,%%xmm14		\n\t"\
		"subpd	%%xmm4,%%xmm6									\n\t	addpd	%%xmm13,%%xmm15		\n\t"\
		"subpd	%%xmm5,%%xmm7									\n\t	movaps	%%xmm12,    (%%r8 )	\n\t"\
		"addpd	%%xmm4,%%xmm4									\n\t	movaps	%%xmm13,0x10(%%r8 )	\n\t"\
		"addpd	%%xmm5,%%xmm5									\n\t	movaps	%%xmm8 ,%%xmm12		\n\t"\
		"addpd	%%xmm6,%%xmm4									\n\t	movaps	%%xmm9 ,%%xmm13		\n\t"\
		"addpd	%%xmm7,%%xmm5									\n\t	subpd	%%xmm10,%%xmm8 		\n\t"\
		"addpd	    (%%rsi),%%xmm4								\n\t	subpd	%%xmm11,%%xmm9 		\n\t"\
		"addpd	0x10(%%rsi),%%xmm5								\n\t	mulpd	0x20(%%r10),%%xmm8 	\n\t"\
		"movaps	%%xmm4,    (%%rdi)								\n\t	mulpd	0x20(%%r10),%%xmm9 	\n\t"\
		"movaps	%%xmm5,0x10(%%rdi)								\n\t	mulpd	0x30(%%r10),%%xmm10	\n\t"\
		"mulpd	0x10(%%rax),%%xmm6								\n\t	mulpd	0x30(%%r10),%%xmm11	\n\t"\
		"mulpd	0x10(%%rax),%%xmm7								\n\t	mulpd	0x40(%%r10),%%xmm12	\n\t"\
		"subpd	    (%%rsi),%%xmm4								\n\t	mulpd	0x40(%%r10),%%xmm13	\n\t"\
		"subpd	0x10(%%rsi),%%xmm5								\n\t	addpd	%%xmm8 ,%%xmm10		\n\t"\
		"mulpd	    (%%rax),%%xmm4								\n\t	addpd	%%xmm9 ,%%xmm11		\n\t"\
		"mulpd	    (%%rax),%%xmm5								\n\t	subpd	%%xmm12,%%xmm8 		\n\t"\
		"addpd	    (%%rdi),%%xmm4								\n\t	subpd	%%xmm13,%%xmm9 		\n\t"\
		"addpd	0x10(%%rdi),%%xmm5								\n\t	movaps	    (%%r8 ),%%xmm12	\n\t"\
		"subpd	%%xmm6,%%xmm4									\n\t	movaps	0x10(%%r8 ),%%xmm13	\n\t"\
		"subpd	%%xmm7,%%xmm5									\n\t	movq	%[__d1],%%r10		\n\t"\
		"mulpd	(%%r14),%%xmm6									\n\t	subpd	%%xmm11,%%xmm14		\n\t"\
		"mulpd	(%%r14),%%xmm7									\n\t	subpd	%%xmm10,%%xmm15		\n\t"\
		"addpd	%%xmm4,%%xmm6									\n\t	mulpd	(%%r14),%%xmm11		\n\t"\
		"addpd	%%xmm5,%%xmm7									\n\t	mulpd	(%%r14),%%xmm10		\n\t"\
		"movaps	%%xmm4,    (%%rsi)								\n\t	movq	%[__d4],%%r13		\n\t"\
		"movaps	%%xmm5,0x10(%%rsi)								\n\t	movaps	%%xmm14,    (%%r10)	\n\t"\
		"movaps	%%xmm0,%%xmm4									\n\t	movaps	%%xmm15,0x10(%%r13)	\n\t"\
		"movaps	%%xmm1,%%xmm5									\n\t	addpd	%%xmm14,%%xmm11		\n\t"\
		"subpd	%%xmm2,%%xmm0									\n\t	addpd	%%xmm15,%%xmm10		\n\t"\
		"subpd	%%xmm3,%%xmm1									\n\t	movaps	%%xmm11,    (%%r13)	\n\t"\
		"mulpd	0x20(%%rax),%%xmm0								\n\t	movaps	%%xmm10,0x10(%%r10)	\n\t"\
		"mulpd	0x20(%%rax),%%xmm1								\n\t	movq	%[__d2],%%r11		\n\t"\
		"mulpd	0x30(%%rax),%%xmm2								\n\t	movq	%[__d3],%%r12		\n\t"\
		"mulpd	0x30(%%rax),%%xmm3								\n\t	subpd	%%xmm9 ,%%xmm12		\n\t"\
		"mulpd	0x40(%%rax),%%xmm4								\n\t	subpd	%%xmm8 ,%%xmm13		\n\t"\
		"mulpd	0x40(%%rax),%%xmm5								\n\t								\n\t"\
		"addpd	%%xmm0,%%xmm2									\n\t	mulpd	(%%r14),%%xmm9 		\n\t"\
		"addpd	%%xmm1,%%xmm3									\n\t	mulpd	(%%r14),%%xmm8 		\n\t"\
		"subpd	%%xmm4,%%xmm0									\n\t	movaps	%%xmm12,    (%%r11)	\n\t"\
		"subpd	%%xmm5,%%xmm1									\n\t	movaps	%%xmm13,0x10(%%r12)	\n\t"\
		"movaps	    (%%rsi),%%xmm4								\n\t	addpd	%%xmm12,%%xmm9 		\n\t"\
		"movaps	0x10(%%rsi),%%xmm5								\n\t	addpd	%%xmm13,%%xmm8 		\n\t"\
		"movq	%[__c1],%%rax									\n\t	movaps	%%xmm9 ,    (%%r12)	\n\t"\
		"movq	%[__c4],%%rdx									\n\t	movaps	%%xmm8 ,0x10(%%r11)	\n\t"\
		"subpd	%%xmm3,%%xmm6									\n\t"\
		"subpd	%%xmm2,%%xmm7									\n\t"\
		"addpd	%%xmm3,%%xmm3									\n\t"\
		"addpd	%%xmm2,%%xmm2									\n\t"\
		"movaps	%%xmm6,    (%%rax)								\n\t"\
		"movaps	%%xmm7,0x10(%%rdx)								\n\t"\
		"addpd	%%xmm6,%%xmm3									\n\t"\
		"addpd	%%xmm7,%%xmm2									\n\t"\
		"movaps	%%xmm3,    (%%rdx)								\n\t"\
		"movaps	%%xmm2,0x10(%%rax)								\n\t"\
		"movq	%[__c2],%%rbx									\n\t"\
		"movq	%[__c3],%%rcx									\n\t"\
		"subpd	%%xmm1,%%xmm4									\n\t"\
		"subpd	%%xmm0,%%xmm5									\n\t"\
		"addpd	%%xmm1,%%xmm1									\n\t"\
		"addpd	%%xmm0,%%xmm0									\n\t"\
		"movaps	%%xmm4,    (%%rbx)								\n\t"\
		"movaps	%%xmm5,0x10(%%rcx)								\n\t"\
		"addpd	%%xmm4,%%xmm1									\n\t"\
		"addpd	%%xmm5,%%xmm0									\n\t"\
		"movaps	%%xmm1,    (%%rcx)								\n\t"\
		"movaps	%%xmm0,0x10(%%rbx)								\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__p01] "m" (Xp01)\
		 ,[__p04] "m" (Xp04)\
		 ,[__r00] "m" (Xr00)\
		 ,[__r10] "m" (Xr10)\
		 ,[__r20] "m" (Xr20)\
		 ,[__r30] "m" (Xr30)\
		 ,[__cc1] "m" (Xcc1)\
		 ,[__a0] "m" (Xa0)\
		 ,[__a1] "m" (Xa1)\
		 ,[__a2] "m" (Xa2)\
		 ,[__a3] "m" (Xa3)\
		 ,[__a4] "m" (Xa4)\
		 ,[__b0] "m" (Xb0)\
		 ,[__b1] "m" (Xb1)\
		 ,[__b2] "m" (Xb2)\
		 ,[__b3] "m" (Xb3)\
		 ,[__b4] "m" (Xb4)\
		 ,[__c0] "m" (Xc0)\
		 ,[__c1] "m" (Xc1)\
		 ,[__c2] "m" (Xc2)\
		 ,[__c3] "m" (Xc3)\
		 ,[__c4] "m" (Xc4)\
		 ,[__d0] "m" (Xd0)\
		 ,[__d1] "m" (Xd1)\
		 ,[__d2] "m" (Xd2)\
		 ,[__d3] "m" (Xd3)\
		 ,[__d4] "m" (Xd4)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r8","r9","r10","r11","r12","r13","r14","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
	);\
	}


	#define	SSE2_RADIX20_DIF_NOTWIDDLE(Xadd0,Xp01,Xp04,Xp08,Xp16,Xa0,Xa1,Xa2,Xa3,Xa4,Xb0,Xb1,Xb2,Xb3,Xb4,Xc0,Xc1,Xc2,Xc3,Xc4,Xd0,Xd1,Xd2,Xd3,Xd4,Xcc1,Xr00,Xr10,Xr20,Xr30)\
	{\
	__asm__ volatile (\
		"/*	SSE2_RADIX_05_DFT_0TWIDDLE(Xa0,r00*/\n\t"\
		"movq	%[__r00],%%rdi					\n\t"\
		"movq	%[__a0],%%rsi					\n\t"\
		"movq	%[__a1],%%rax					\n\t	movq	$0x500,%%r14		\n\t"\
		"movq	%[__a2],%%rbx					\n\t"\
		"movq	%[__a3],%%rcx					\n\t"\
		"movq	%[__a4],%%rdx					\n\t	addq	%%rdi,%%r14	/* two */\n\t"\
		"movaps	    (%%rax),%%xmm0				\n\t"\
		"movaps	0x10(%%rax),%%xmm1				\n\t	/*	SSE2_RADIX_05_DFT_0TWIDDLE(Xc0,r20)	*/\n\t"\
		"movaps	    (%%rbx),%%xmm2				\n\t	movq	%[__c0],%%r9 		\n\t"\
		"movaps	0x10(%%rbx),%%xmm3				\n\t	movq	%[__c1],%%r10		\n\t"\
		"movaps	    (%%rcx),%%xmm4				\n\t	movq	%[__c2],%%r11		\n\t"\
		"movaps	0x10(%%rcx),%%xmm5				\n\t	movq	%[__c3],%%r12		\n\t"\
		"movaps	    (%%rdx),%%xmm6				\n\t	movq	%[__c4],%%r13		\n\t"\
		"movaps	0x10(%%rdx),%%xmm7				\n\t	movq	%[__r20],%%r8 		\n\t"\
		"subpd	%%xmm6,%%xmm0					\n\t	movaps	    (%%r10),%%xmm8 	\n\t"\
		"subpd	%%xmm7,%%xmm1					\n\t	movaps	0x10(%%r10),%%xmm9 	\n\t"\
		"addpd	%%xmm6,%%xmm6					\n\t	movaps	    (%%r11),%%xmm10	\n\t"\
		"addpd	%%xmm7,%%xmm7					\n\t	movaps	0x10(%%r11),%%xmm11	\n\t"\
		"addpd	%%xmm0,%%xmm6					\n\t	movaps	    (%%r12),%%xmm12	\n\t"\
		"addpd	%%xmm1,%%xmm7					\n\t	movaps	0x10(%%r12),%%xmm13	\n\t"\
		"subpd	%%xmm4,%%xmm2					\n\t	movaps	    (%%r13),%%xmm14	\n\t"\
		"subpd	%%xmm5,%%xmm3					\n\t	movaps	0x10(%%r13),%%xmm15	\n\t"\
		"mulpd	(%%r14),%%xmm4					\n\t	subpd	%%xmm14,%%xmm8 		\n\t"\
		"mulpd	(%%r14),%%xmm5					\n\t	subpd	%%xmm15,%%xmm9 		\n\t"\
		"addpd	%%xmm2,%%xmm4					\n\t	mulpd	(%%r14),%%xmm14		\n\t"\
		"addpd	%%xmm3,%%xmm5					\n\t	mulpd	(%%r14),%%xmm15		\n\t"\
		"movq	%[__cc1],%%rax					\n\t	addpd	%%xmm8 ,%%xmm14		\n\t"\
		"subpd	%%xmm4,%%xmm6					\n\t	addpd	%%xmm9 ,%%xmm15		\n\t"\
		"subpd	%%xmm5,%%xmm7					\n\t	subpd	%%xmm12,%%xmm10		\n\t"\
		"mulpd	(%%r14),%%xmm4					\n\t	subpd	%%xmm13,%%xmm11		\n\t"\
		"mulpd	(%%r14),%%xmm5					\n\t	movq	%[__cc1],%%r10		\n\t"\
		"addpd	%%xmm6,%%xmm4					\n\t	mulpd	(%%r14),%%xmm12		\n\t"\
		"addpd	%%xmm7,%%xmm5					\n\t	mulpd	(%%r14),%%xmm13		\n\t"\
		"addpd	    (%%rsi),%%xmm4				\n\t	addpd	%%xmm10,%%xmm12		\n\t"\
		"addpd	0x10(%%rsi),%%xmm5				\n\t	addpd	%%xmm11,%%xmm13		\n\t"\
		"movaps	%%xmm4,    (%%rdi)				\n\t	subpd	%%xmm12,%%xmm14		\n\t"\
		"movaps	%%xmm5,0x10(%%rdi)				\n\t	subpd	%%xmm13,%%xmm15		\n\t"\
		"mulpd	0x10(%%rax),%%xmm6				\n\t	addpd	%%xmm12,%%xmm12		\n\t"\
		"mulpd	0x10(%%rax),%%xmm7				\n\t	addpd	%%xmm13,%%xmm13		\n\t"\
		"subpd	    (%%rsi),%%xmm4				\n\t	addpd	%%xmm14,%%xmm12		\n\t"\
		"subpd	0x10(%%rsi),%%xmm5				\n\t	addpd	%%xmm15,%%xmm13		\n\t"\
		"mulpd	    (%%rax),%%xmm4				\n\t	addpd	    (%%r9 ),%%xmm12	\n\t"\
		"mulpd	    (%%rax),%%xmm5				\n\t	addpd	0x10(%%r9 ),%%xmm13	\n\t"\
		"addpd	    (%%rdi),%%xmm4				\n\t	movaps	%%xmm12,    (%%r8 )	\n\t"\
		"addpd	0x10(%%rdi),%%xmm5				\n\t	movaps	%%xmm13,0x10(%%r8 )	\n\t"\
		"subpd	%%xmm6,%%xmm4					\n\t	mulpd	0x10(%%r10),%%xmm14	\n\t"\
		"subpd	%%xmm7,%%xmm5					\n\t	mulpd	0x10(%%r10),%%xmm15	\n\t"\
		"mulpd	(%%r14),%%xmm6					\n\t	subpd	    (%%r9 ),%%xmm12	\n\t"\
		"mulpd	(%%r14),%%xmm7					\n\t	subpd	0x10(%%r9 ),%%xmm13	\n\t"\
		"addpd	%%xmm4,%%xmm6					\n\t	mulpd	    (%%r10),%%xmm12	\n\t"\
		"addpd	%%xmm5,%%xmm7					\n\t	mulpd	    (%%r10),%%xmm13	\n\t"\
		"movaps	%%xmm4,    (%%rsi)				\n\t	addpd	    (%%r8 ),%%xmm12	\n\t"\
		"movaps	%%xmm5,0x10(%%rsi)				\n\t	addpd	0x10(%%r8 ),%%xmm13	\n\t"\
		"movaps	%%xmm0,%%xmm4					\n\t	subpd	%%xmm14,%%xmm12		\n\t"\
		"movaps	%%xmm1,%%xmm5					\n\t	subpd	%%xmm15,%%xmm13		\n\t"\
		"subpd	%%xmm2,%%xmm0					\n\t	mulpd	(%%r14),%%xmm14		\n\t"\
		"subpd	%%xmm3,%%xmm1					\n\t	mulpd	(%%r14),%%xmm15		\n\t"\
		"mulpd	0x20(%%rax),%%xmm0				\n\t	addpd	%%xmm12,%%xmm14		\n\t"\
		"mulpd	0x20(%%rax),%%xmm1				\n\t	addpd	%%xmm13,%%xmm15		\n\t"\
		"mulpd	0x30(%%rax),%%xmm2				\n\t	movaps	%%xmm12,    (%%r9 )	\n\t"\
		"mulpd	0x30(%%rax),%%xmm3				\n\t	movaps	%%xmm13,0x10(%%r9 )	\n\t"\
		"mulpd	0x40(%%rax),%%xmm4				\n\t	movaps	%%xmm8 ,%%xmm12		\n\t"\
		"mulpd	0x40(%%rax),%%xmm5				\n\t	movaps	%%xmm9 ,%%xmm13		\n\t"\
		"addpd	%%xmm0,%%xmm2					\n\t	subpd	%%xmm10,%%xmm8 		\n\t"\
		"addpd	%%xmm1,%%xmm3					\n\t	subpd	%%xmm11,%%xmm9 		\n\t"\
		"subpd	%%xmm4,%%xmm0					\n\t	mulpd	0x20(%%r10),%%xmm8 	\n\t"\
		"subpd	%%xmm5,%%xmm1					\n\t	mulpd	0x20(%%r10),%%xmm9 	\n\t"\
		"movaps	    (%%rsi),%%xmm4				\n\t	mulpd	0x30(%%r10),%%xmm10	\n\t"\
		"movaps	0x10(%%rsi),%%xmm5				\n\t	mulpd	0x30(%%r10),%%xmm11	\n\t"\
		"subpd	%%xmm3,%%xmm6					\n\t	mulpd	0x40(%%r10),%%xmm12	\n\t"\
		"subpd	%%xmm2,%%xmm7					\n\t	mulpd	0x40(%%r10),%%xmm13	\n\t"\
		"addpd	%%xmm3,%%xmm3					\n\t	addpd	%%xmm8 ,%%xmm10		\n\t"\
		"addpd	%%xmm2,%%xmm2					\n\t	addpd	%%xmm9 ,%%xmm11		\n\t"\
		"movaps	%%xmm6,0x020(%%rdi)				\n\t	subpd	%%xmm12,%%xmm8 		\n\t"\
		"movaps	%%xmm7,0x070(%%rdi)				\n\t	subpd	%%xmm13,%%xmm9 		\n\t"\
		"addpd	%%xmm6,%%xmm3					\n\t	movaps	    (%%r9 ),%%xmm12	\n\t"\
		"addpd	%%xmm7,%%xmm2					\n\t	movaps	0x10(%%r9 ),%%xmm13	\n\t"\
		"movaps	%%xmm3,0x060(%%rdi)				\n\t	subpd	%%xmm11,%%xmm14		\n\t"\
		"movaps	%%xmm2,0x030(%%rdi)				\n\t	subpd	%%xmm10,%%xmm15		\n\t"\
		"subpd	%%xmm1,%%xmm4					\n\t	addpd	%%xmm11,%%xmm11		\n\t"\
		"subpd	%%xmm0,%%xmm5					\n\t	addpd	%%xmm10,%%xmm10		\n\t"\
		"addpd	%%xmm1,%%xmm1					\n\t	movaps	%%xmm14,0x020(%%r8 )\n\t"\
		"addpd	%%xmm0,%%xmm0					\n\t	movaps	%%xmm15,0x070(%%r8 )\n\t"\
		"movaps	%%xmm4,0x040(%%rdi)				\n\t	addpd	%%xmm14,%%xmm11		\n\t"\
		"movaps	%%xmm5,0x090(%%rdi)				\n\t	addpd	%%xmm15,%%xmm10		\n\t"\
		"addpd	%%xmm4,%%xmm1					\n\t	movaps	%%xmm11,0x060(%%r8 )\n\t"\
		"addpd	%%xmm5,%%xmm0					\n\t	movaps	%%xmm10,0x030(%%r8 )\n\t"\
		"movaps	%%xmm1,0x080(%%rdi)				\n\t	subpd	%%xmm9 ,%%xmm12		\n\t"\
		"movaps	%%xmm0,0x050(%%rdi)				\n\t	subpd	%%xmm8 ,%%xmm13		\n\t"\
		"										\n\t	addpd	%%xmm9 ,%%xmm9 		\n\t"\
		"/*	SSE2_RADIX_05_DFT_0TWIDDLE(Xb0,r10*/\n\t	addpd	%%xmm8 ,%%xmm8 		\n\t"\
		"movq	%[__b0],%%rsi					\n\t	movaps	%%xmm12,0x040(%%r8 )\n\t"\
		"movq	%[__b1],%%rax					\n\t	movaps	%%xmm13,0x090(%%r8 )\n\t"\
		"movq	%[__b2],%%rbx					\n\t	addpd	%%xmm12,%%xmm9 		\n\t"\
		"movq	%[__b3],%%rcx					\n\t	addpd	%%xmm13,%%xmm8 		\n\t"\
		"movq	%[__b4],%%rdx					\n\t	movaps	%%xmm9 ,0x080(%%r8 )\n\t"\
		"movq	%[__r10],%%rdi					\n\t	movaps	%%xmm8 ,0x050(%%r8 )\n\t"\
		"movaps	    (%%rax),%%xmm0				\n\t								\n\t"\
		"movaps	0x10(%%rax),%%xmm1				\n\t	/*	SSE2_RADIX_05_DFT_0TWIDDLE(Xd0,r30)	*/\n\t"\
		"movaps	    (%%rbx),%%xmm2				\n\t	movq	%[__d0],%%r9 		\n\t"\
		"movaps	0x10(%%rbx),%%xmm3				\n\t	movq	%[__d1],%%r10		\n\t"\
		"movaps	    (%%rcx),%%xmm4				\n\t	movq	%[__d2],%%r11		\n\t"\
		"movaps	0x10(%%rcx),%%xmm5				\n\t	movq	%[__d3],%%r12		\n\t"\
		"movaps	    (%%rdx),%%xmm6				\n\t	movq	%[__d4],%%r13		\n\t"\
		"movaps	0x10(%%rdx),%%xmm7				\n\t	movq	%[__r30],%%r8 		\n\t"\
		"subpd	%%xmm6,%%xmm0					\n\t	movaps	    (%%r10),%%xmm8 	\n\t"\
		"subpd	%%xmm7,%%xmm1					\n\t	movaps	0x10(%%r10),%%xmm9 	\n\t"\
		"addpd	%%xmm6,%%xmm6					\n\t	movaps	    (%%r11),%%xmm10	\n\t"\
		"addpd	%%xmm7,%%xmm7					\n\t	movaps	0x10(%%r11),%%xmm11	\n\t"\
		"addpd	%%xmm0,%%xmm6					\n\t	movaps	    (%%r12),%%xmm12	\n\t"\
		"addpd	%%xmm1,%%xmm7					\n\t	movaps	0x10(%%r12),%%xmm13	\n\t"\
		"subpd	%%xmm4,%%xmm2					\n\t	movaps	    (%%r13),%%xmm14	\n\t"\
		"subpd	%%xmm5,%%xmm3					\n\t	movaps	0x10(%%r13),%%xmm15	\n\t"\
		"mulpd	(%%r14),%%xmm4					\n\t	subpd	%%xmm14,%%xmm8 		\n\t"\
		"mulpd	(%%r14),%%xmm5					\n\t	subpd	%%xmm15,%%xmm9 		\n\t"\
		"addpd	%%xmm2,%%xmm4					\n\t	mulpd	(%%r14),%%xmm14		\n\t"\
		"addpd	%%xmm3,%%xmm5					\n\t	mulpd	(%%r14),%%xmm15		\n\t"\
		"movq	%[__cc1],%%rax					\n\t	addpd	%%xmm8 ,%%xmm14		\n\t"\
		"subpd	%%xmm4,%%xmm6					\n\t	addpd	%%xmm9 ,%%xmm15		\n\t"\
		"subpd	%%xmm5,%%xmm7					\n\t	subpd	%%xmm12,%%xmm10		\n\t"\
		"mulpd	(%%r14),%%xmm4					\n\t	subpd	%%xmm13,%%xmm11		\n\t"\
		"mulpd	(%%r14),%%xmm5					\n\t	movq	%[__cc1],%%r10		\n\t"\
		"addpd	%%xmm6,%%xmm4					\n\t	mulpd	(%%r14),%%xmm12		\n\t"\
		"addpd	%%xmm7,%%xmm5					\n\t	mulpd	(%%r14),%%xmm13		\n\t"\
		"addpd	    (%%rsi),%%xmm4				\n\t	addpd	%%xmm10,%%xmm12		\n\t"\
		"addpd	0x10(%%rsi),%%xmm5				\n\t	addpd	%%xmm11,%%xmm13		\n\t"\
		"movaps	%%xmm4,    (%%rdi)				\n\t	subpd	%%xmm12,%%xmm14		\n\t"\
		"movaps	%%xmm5,0x10(%%rdi)				\n\t	subpd	%%xmm13,%%xmm15		\n\t"\
		"mulpd	0x10(%%rax),%%xmm6				\n\t	addpd	%%xmm12,%%xmm12		\n\t"\
		"mulpd	0x10(%%rax),%%xmm7				\n\t	addpd	%%xmm13,%%xmm13		\n\t"\
		"subpd	    (%%rsi),%%xmm4				\n\t	addpd	%%xmm14,%%xmm12		\n\t"\
		"subpd	0x10(%%rsi),%%xmm5				\n\t	addpd	%%xmm15,%%xmm13		\n\t"\
		"mulpd	    (%%rax),%%xmm4				\n\t	addpd	    (%%r9 ),%%xmm12	\n\t"\
		"mulpd	    (%%rax),%%xmm5				\n\t	addpd	0x10(%%r9 ),%%xmm13	\n\t"\
		"addpd	    (%%rdi),%%xmm4				\n\t	movaps	%%xmm12,    (%%r8 )	\n\t"\
		"addpd	0x10(%%rdi),%%xmm5				\n\t	movaps	%%xmm13,0x10(%%r8 )	\n\t"\
		"subpd	%%xmm6,%%xmm4					\n\t	mulpd	0x10(%%r10),%%xmm14	\n\t"\
		"subpd	%%xmm7,%%xmm5					\n\t	mulpd	0x10(%%r10),%%xmm15	\n\t"\
		"mulpd	(%%r14),%%xmm6					\n\t	subpd	    (%%r9 ),%%xmm12	\n\t"\
		"mulpd	(%%r14),%%xmm7					\n\t	subpd	0x10(%%r9 ),%%xmm13	\n\t"\
		"addpd	%%xmm4,%%xmm6					\n\t	mulpd	    (%%r10),%%xmm12	\n\t"\
		"addpd	%%xmm5,%%xmm7					\n\t	mulpd	    (%%r10),%%xmm13	\n\t"\
		"movaps	%%xmm4,    (%%rsi)				\n\t	addpd	    (%%r8 ),%%xmm12	\n\t"\
		"movaps	%%xmm5,0x10(%%rsi)				\n\t	addpd	0x10(%%r8 ),%%xmm13	\n\t"\
		"movaps	%%xmm0,%%xmm4					\n\t	subpd	%%xmm14,%%xmm12		\n\t"\
		"movaps	%%xmm1,%%xmm5					\n\t	subpd	%%xmm15,%%xmm13		\n\t"\
		"subpd	%%xmm2,%%xmm0					\n\t	mulpd	(%%r14),%%xmm14		\n\t"\
		"subpd	%%xmm3,%%xmm1					\n\t	mulpd	(%%r14),%%xmm15		\n\t"\
		"mulpd	0x20(%%rax),%%xmm0				\n\t	addpd	%%xmm12,%%xmm14		\n\t"\
		"mulpd	0x20(%%rax),%%xmm1				\n\t	addpd	%%xmm13,%%xmm15		\n\t"\
		"mulpd	0x30(%%rax),%%xmm2				\n\t	movaps	%%xmm12,    (%%r9 )	\n\t"\
		"mulpd	0x30(%%rax),%%xmm3				\n\t	movaps	%%xmm13,0x10(%%r9 )	\n\t"\
		"mulpd	0x40(%%rax),%%xmm4				\n\t	movaps	%%xmm8 ,%%xmm12		\n\t"\
		"mulpd	0x40(%%rax),%%xmm5				\n\t	movaps	%%xmm9 ,%%xmm13		\n\t"\
		"addpd	%%xmm0,%%xmm2					\n\t	subpd	%%xmm10,%%xmm8 		\n\t"\
		"addpd	%%xmm1,%%xmm3					\n\t	subpd	%%xmm11,%%xmm9 		\n\t"\
		"subpd	%%xmm4,%%xmm0					\n\t	mulpd	0x20(%%r10),%%xmm8 	\n\t"\
		"subpd	%%xmm5,%%xmm1					\n\t	mulpd	0x20(%%r10),%%xmm9 	\n\t"\
		"movaps	    (%%rsi),%%xmm4				\n\t	mulpd	0x30(%%r10),%%xmm10	\n\t"\
		"movaps	0x10(%%rsi),%%xmm5				\n\t	mulpd	0x30(%%r10),%%xmm11	\n\t"\
		"subpd	%%xmm3,%%xmm6					\n\t	mulpd	0x40(%%r10),%%xmm12	\n\t"\
		"subpd	%%xmm2,%%xmm7					\n\t	mulpd	0x40(%%r10),%%xmm13	\n\t"\
		"addpd	%%xmm3,%%xmm3					\n\t	addpd	%%xmm8 ,%%xmm10		\n\t"\
		"addpd	%%xmm2,%%xmm2					\n\t	addpd	%%xmm9 ,%%xmm11		\n\t"\
		"movaps	%%xmm6,0x020(%%rdi)				\n\t	subpd	%%xmm12,%%xmm8 		\n\t"\
		"movaps	%%xmm7,0x070(%%rdi)				\n\t	subpd	%%xmm13,%%xmm9 		\n\t"\
		"addpd	%%xmm6,%%xmm3					\n\t	movaps	    (%%r9 ),%%xmm12	\n\t"\
		"addpd	%%xmm7,%%xmm2					\n\t	movaps	0x10(%%r9 ),%%xmm13	\n\t"\
		"movaps	%%xmm3,0x060(%%rdi)				\n\t	subpd	%%xmm11,%%xmm14		\n\t"\
		"movaps	%%xmm2,0x030(%%rdi)				\n\t	subpd	%%xmm10,%%xmm15		\n\t"\
		"subpd	%%xmm1,%%xmm4					\n\t	addpd	%%xmm11,%%xmm11		\n\t"\
		"subpd	%%xmm0,%%xmm5					\n\t	addpd	%%xmm10,%%xmm10		\n\t"\
		"addpd	%%xmm1,%%xmm1					\n\t	movaps	%%xmm14,0x020(%%r8 )\n\t"\
		"addpd	%%xmm0,%%xmm0					\n\t	movaps	%%xmm15,0x070(%%r8 )\n\t"\
		"movaps	%%xmm4,0x040(%%rdi)				\n\t	addpd	%%xmm14,%%xmm11		\n\t"\
		"movaps	%%xmm5,0x090(%%rdi)				\n\t	addpd	%%xmm15,%%xmm10		\n\t"\
		"addpd	%%xmm4,%%xmm1					\n\t	movaps	%%xmm11,0x060(%%r8 )\n\t"\
		"addpd	%%xmm5,%%xmm0					\n\t	movaps	%%xmm10,0x030(%%r8 )\n\t"\
		"movaps	%%xmm1,0x080(%%rdi)				\n\t	subpd	%%xmm9 ,%%xmm12		\n\t"\
		"movaps	%%xmm0,0x050(%%rdi)				\n\t	subpd	%%xmm8 ,%%xmm13		\n\t"\
		"										\n\t	addpd	%%xmm9 ,%%xmm9 		\n\t"\
		"movq	%[__add0],%%rax					\n\t	addpd	%%xmm8 ,%%xmm8 		\n\t"\
		"movq	%[__add0],%%rdx					\n\t	movaps	%%xmm12,0x040(%%r8 )\n\t"\
		"movslq	%[__p01],%%rsi					\n\t	movaps	%%xmm13,0x090(%%r8 )\n\t"\
		"shlq	$3,%%rsi						\n\t	addpd	%%xmm12,%%xmm9 		\n\t"\
		"addq	%%rsi,%%rdx						\n\t	addpd	%%xmm13,%%xmm8 		\n\t"\
		"movq	%%rdx,%%rbx						\n\t	movaps	%%xmm9 ,0x080(%%r8 )\n\t"\
		"addq	%%rsi,%%rdx						\n\t	movaps	%%xmm8 ,0x050(%%r8 )\n\t"\
		"movq	%%rdx,%%rcx						\n\t								\n\t"\
		"addq	%%rsi,%%rdx						\n\t	movq	%[__r00],%%rdi		\n\t"\
		"movslq	%[__p04],%%r8 					\n\t	movslq	%[__p16],%%r10		\n\t"\
		"shlq	$3,%%r8 						\n\t	shlq	$3,%%r10			\n\t"\
		"/*	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_C(r00, 0x0a0, 0x140, eax,ebx,edx,ecx) */\n\t"\
		"movq	$0x0a0,%%rsi					\n\t	/*	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_C(r02, 0x0a0, 0x140, eax,ebx,edx,ecx) */\n\t"\
		"addq	%%rdi,%%rsi						\n\t	movq	%%r10,%%r11				\n\t"\
		"movaps	     (%%rdi),%%xmm4				\n\t	movq	%%r10,%%r12				\n\t"\
		"movaps	     (%%rsi),%%xmm6				\n\t	movq	%%r10,%%r13				\n\t"\
		"movaps	0x010(%%rdi),%%xmm5				\n\t	addq	%%rax,%%r10	/* &a[j1+p16] */\n\t"\
		"movaps	0x010(%%rsi),%%xmm7				\n\t	addq	%%rbx,%%r11				\n\t"\
		"movaps	0x140(%%rdi),%%xmm0				\n\t	addq	%%rcx,%%r12				\n\t"\
		"movaps	0x140(%%rsi),%%xmm2				\n\t	addq	%%rdx,%%r13				\n\t"\
		"movaps	0x150(%%rdi),%%xmm1				\n\t	addq	$0x20,%%rdi	/* r02 */	\n\t"\
		"movaps	0x150(%%rsi),%%xmm3				\n\t	addq	$0x20,%%rsi				\n\t"\
		"subpd	%%xmm0,%%xmm4					\n\t	movaps	     (%%rdi),%%xmm12	\n\t"\
		"subpd	%%xmm2,%%xmm6					\n\t	movaps	     (%%rsi),%%xmm14	\n\t"\
		"subpd	%%xmm1,%%xmm5					\n\t	movaps	0x010(%%rdi),%%xmm13	\n\t"\
		"subpd	%%xmm3,%%xmm7					\n\t	movaps	0x010(%%rsi),%%xmm15	\n\t"\
		"addpd	%%xmm0,%%xmm0					\n\t	movaps	0x140(%%rdi),%%xmm8 	\n\t"\
		"addpd	%%xmm2,%%xmm2					\n\t	movaps	0x140(%%rsi),%%xmm10	\n\t"\
		"addpd	%%xmm1,%%xmm1					\n\t	movaps	0x150(%%rdi),%%xmm9 	\n\t"\
		"addpd	%%xmm3,%%xmm3					\n\t	movaps	0x150(%%rsi),%%xmm11	\n\t"\
		"addpd	%%xmm4,%%xmm0					\n\t	subpd	%%xmm8 ,%%xmm12			\n\t"\
		"addpd	%%xmm6,%%xmm2					\n\t	subpd	%%xmm10,%%xmm14			\n\t"\
		"addpd	%%xmm5,%%xmm1					\n\t	subpd	%%xmm9 ,%%xmm13			\n\t"\
		"addpd	%%xmm7,%%xmm3					\n\t	subpd	%%xmm11,%%xmm15			\n\t"\
		"/* Finish radix-4 butterfly: */		\n\t	mulpd	(%%r14),%%xmm8 			\n\t"\
		"subpd	%%xmm2,%%xmm0					\n\t	mulpd	(%%r14),%%xmm10			\n\t"\
		"subpd	%%xmm7,%%xmm4					\n\t	mulpd	(%%r14),%%xmm9 			\n\t"\
		"subpd	%%xmm3,%%xmm1					\n\t	mulpd	(%%r14),%%xmm11			\n\t"\
		"subpd	%%xmm6,%%xmm5					\n\t	addpd	%%xmm12,%%xmm8 			\n\t"\
		"movaps	%%xmm0,     (%%rbx)				\n\t	addpd	%%xmm14,%%xmm10			\n\t"\
		"movaps	%%xmm4,     (%%rdx)				\n\t	addpd	%%xmm13,%%xmm9 			\n\t"\
		"movaps	%%xmm1,0x010(%%rbx)				\n\t	addpd	%%xmm15,%%xmm11			\n\t"\
		"movaps	%%xmm5,0x010(%%rcx)				\n\t	/* Finish radix-4 butterfly and store results into main-array slots: */\n\t"\
		"mulpd	(%%r14),%%xmm2					\n\t	subpd	%%xmm10,%%xmm8 			\n\t"\
		"mulpd	(%%r14),%%xmm7					\n\t	subpd	%%xmm15,%%xmm12			\n\t"\
		"mulpd	(%%r14),%%xmm3					\n\t	subpd	%%xmm11,%%xmm9 			\n\t"\
		"mulpd	(%%r14),%%xmm6					\n\t	subpd	%%xmm14,%%xmm13			\n\t"\
		"addpd	%%xmm0,%%xmm2					\n\t	movaps	%%xmm8 ,     (%%r11)	\n\t"\
		"addpd	%%xmm4,%%xmm7					\n\t	movaps	%%xmm12,     (%%r13)	\n\t"\
		"addpd	%%xmm1,%%xmm3					\n\t	movaps	%%xmm9 ,0x010(%%r11)	\n\t"\
		"addpd	%%xmm5,%%xmm6					\n\t	movaps	%%xmm13,0x010(%%r12)	\n\t"\
		"movaps	%%xmm2,     (%%rax)				\n\t	addpd	%%xmm10,%%xmm10			\n\t"\
		"movaps	%%xmm7,     (%%rcx)				\n\t	addpd	%%xmm15,%%xmm15			\n\t"\
		"movaps	%%xmm3,0x010(%%rax)				\n\t	addpd	%%xmm11,%%xmm11			\n\t"\
		"movaps	%%xmm6,0x010(%%rdx)				\n\t	addpd	%%xmm14,%%xmm14			\n\t"\
		"										\n\t	addpd	%%xmm8 ,%%xmm10			\n\t"\
		"/*	add0,1,2,3 = &a[j1+p04]+p3,2,1,0 */	\n\t	addpd	%%xmm12,%%xmm15			\n\t"\
		"addq	%%r8 ,%%rax		/* &a[j1+p04] */\n\t	addpd	%%xmm9 ,%%xmm11			\n\t"\
		"addq	%%r8 ,%%rbx						\n\t	addpd	%%xmm13,%%xmm14			\n\t"\
		"addq	%%r8 ,%%rcx						\n\t	movaps	%%xmm10,     (%%r10)	\n\t"\
		"addq	%%r8 ,%%rdx						\n\t	movaps	%%xmm15,     (%%r12)	\n\t"\
		"/*	RADIX4_DIF_0TWIDDLE_STRIDE_C(r06) */\n\t	movaps	%%xmm11,0x010(%%r10)	\n\t"\
		"addq	$0x40,%%rdi	/* r06 */			\n\t	movaps	%%xmm14,0x010(%%r13)	\n\t"\
		"addq	$0x40,%%rsi						\n\t									\n\t"\
		"movaps	     (%%rdi),%%xmm4				\n\t	/*	add0,1,2,3 = &a[j1+p12]+p2,3,0,1 */\n\t"\
		"movaps	     (%%rsi),%%xmm6				\n\t	subq	%%r8 ,%%r10		/* &a[j1+p12] */\n\t"\
		"movaps	0x010(%%rdi),%%xmm5				\n\t	subq	%%r8 ,%%r11				\n\t"\
		"movaps	0x010(%%rsi),%%xmm7				\n\t	subq	%%r8 ,%%r12				\n\t"\
		"movaps	0x140(%%rdi),%%xmm0				\n\t	subq	%%r8 ,%%r13				\n\t"\
		"movaps	0x140(%%rsi),%%xmm2				\n\t	/*	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_C(r04, 0x0a0, 0x140, ecx,edx,eax,ebx) */\n\t"\
		"movaps	0x150(%%rdi),%%xmm1				\n\t	subq	$0x20,%%rdi	/* r04 */	\n\t"\
		"movaps	0x150(%%rsi),%%xmm3				\n\t	subq	$0x20,%%rsi				\n\t"\
		"subpd	%%xmm0,%%xmm4					\n\t	movaps	     (%%rdi),%%xmm12	\n\t"\
		"subpd	%%xmm2,%%xmm6					\n\t	movaps	     (%%rsi),%%xmm14	\n\t"\
		"subpd	%%xmm1,%%xmm5					\n\t	movaps	0x010(%%rdi),%%xmm13	\n\t"\
		"subpd	%%xmm3,%%xmm7					\n\t	movaps	0x010(%%rsi),%%xmm15	\n\t"\
		"addpd	%%xmm0,%%xmm0					\n\t	movaps	0x140(%%rdi),%%xmm8 	\n\t"\
		"addpd	%%xmm2,%%xmm2					\n\t	movaps	0x140(%%rsi),%%xmm10	\n\t"\
		"addpd	%%xmm1,%%xmm1					\n\t	movaps	0x150(%%rdi),%%xmm9 	\n\t"\
		"addpd	%%xmm3,%%xmm3					\n\t	movaps	0x150(%%rsi),%%xmm11	\n\t"\
		"addpd	%%xmm4,%%xmm0					\n\t	subpd	%%xmm8 ,%%xmm12			\n\t"\
		"addpd	%%xmm6,%%xmm2					\n\t	subpd	%%xmm10,%%xmm14			\n\t"\
		"addpd	%%xmm5,%%xmm1					\n\t	subpd	%%xmm9 ,%%xmm13			\n\t"\
		"addpd	%%xmm7,%%xmm3					\n\t	subpd	%%xmm11,%%xmm15			\n\t"\
		"/* Finish radix-4 butterfly: */		\n\t	mulpd	(%%r14),%%xmm8 			\n\t"\
		"subpd	%%xmm2,%%xmm0					\n\t	mulpd	(%%r14),%%xmm10			\n\t"\
		"subpd	%%xmm7,%%xmm4					\n\t	mulpd	(%%r14),%%xmm9 			\n\t"\
		"subpd	%%xmm3,%%xmm1					\n\t	mulpd	(%%r14),%%xmm11			\n\t"\
		"subpd	%%xmm6,%%xmm5					\n\t	addpd	%%xmm12,%%xmm8 			\n\t"\
		"movaps	%%xmm0,     (%%rcx)				\n\t	addpd	%%xmm14,%%xmm10			\n\t"\
		"movaps	%%xmm4,     (%%rbx)				\n\t	addpd	%%xmm13,%%xmm9 			\n\t"\
		"movaps	%%xmm1,0x010(%%rcx)				\n\t	addpd	%%xmm15,%%xmm11			\n\t"\
		"movaps	%%xmm5,0x010(%%rax)				\n\t	/* Finish radix-4 butterfly and store results into main-array slots: */\n\t"\
		"mulpd	(%%r14),%%xmm2					\n\t	subpd	%%xmm10,%%xmm8 			\n\t"\
		"mulpd	(%%r14),%%xmm7					\n\t	subpd	%%xmm15,%%xmm12			\n\t"\
		"mulpd	(%%r14),%%xmm3					\n\t	subpd	%%xmm11,%%xmm9 			\n\t"\
		"mulpd	(%%r14),%%xmm6					\n\t	subpd	%%xmm14,%%xmm13			\n\t"\
		"addpd	%%xmm0,%%xmm2					\n\t	movaps	%%xmm8 ,     (%%r13)	\n\t"\
		"addpd	%%xmm4,%%xmm7					\n\t	movaps	%%xmm12,     (%%r10)	\n\t"\
		"addpd	%%xmm1,%%xmm3					\n\t	movaps	%%xmm9 ,0x010(%%r13)	\n\t"\
		"addpd	%%xmm5,%%xmm6					\n\t	movaps	%%xmm13,0x010(%%r11)	\n\t"\
		"movaps	%%xmm2,     (%%rdx)				\n\t	addpd	%%xmm10,%%xmm10			\n\t"\
		"movaps	%%xmm7,     (%%rax)				\n\t	addpd	%%xmm15,%%xmm15			\n\t"\
		"movaps	%%xmm3,0x010(%%rdx)				\n\t	addpd	%%xmm11,%%xmm11			\n\t"\
		"movaps	%%xmm6,0x010(%%rbx)				\n\t	addpd	%%xmm14,%%xmm14			\n\t"\
		"										\n\t	addpd	%%xmm8 ,%%xmm10			\n\t"\
		"/*	add0,1,2,3 = &a[j1+p08]+p1,0,2,3 */	\n\t	addpd	%%xmm12,%%xmm15			\n\t"\
		"addq	%%r8 ,%%rax		/* &a[j1+p08] */\n\t	addpd	%%xmm9 ,%%xmm11			\n\t"\
		"addq	%%r8 ,%%rbx						\n\t	addpd	%%xmm13,%%xmm14			\n\t"\
		"addq	%%r8 ,%%rcx						\n\t	movaps	%%xmm10,     (%%r12)	\n\t"\
		"addq	%%r8 ,%%rdx						\n\t	movaps	%%xmm15,     (%%r11)	\n\t"\
		"/*	RADIX4_DIF_0TWIDDLE_STRIDE_C(r08) */\n\t	movaps	%%xmm11,0x010(%%r12)	\n\t"\
		"addq	$0x40,%%rdi	/* r08 */			\n\t	movaps	%%xmm14,0x010(%%r10)	\n\t"\
		"addq	$0x40,%%rsi						\n\t"\
		"movaps	     (%%rdi),%%xmm4				\n\t"\
		"movaps	     (%%rsi),%%xmm6				\n\t"\
		"movaps	0x010(%%rdi),%%xmm5				\n\t"\
		"movaps	0x010(%%rsi),%%xmm7				\n\t"\
		"movaps	0x140(%%rdi),%%xmm0				\n\t"\
		"movaps	0x140(%%rsi),%%xmm2				\n\t"\
		"movaps	0x150(%%rdi),%%xmm1				\n\t"\
		"movaps	0x150(%%rsi),%%xmm3				\n\t"\
		"subpd	%%xmm0,%%xmm4					\n\t"\
		"subpd	%%xmm2,%%xmm6					\n\t"\
		"subpd	%%xmm1,%%xmm5					\n\t"\
		"subpd	%%xmm3,%%xmm7					\n\t"\
		"addpd	%%xmm0,%%xmm0					\n\t"\
		"addpd	%%xmm2,%%xmm2					\n\t"\
		"addpd	%%xmm1,%%xmm1					\n\t"\
		"addpd	%%xmm3,%%xmm3					\n\t"\
		"addpd	%%xmm4,%%xmm0					\n\t"\
		"addpd	%%xmm6,%%xmm2					\n\t"\
		"addpd	%%xmm5,%%xmm1					\n\t"\
		"addpd	%%xmm7,%%xmm3					\n\t"\
		"/* Finish radix-4 butterfly: */		\n\t"\
		"subpd	%%xmm2,%%xmm0					\n\t"\
		"subpd	%%xmm7,%%xmm4					\n\t"\
		"subpd	%%xmm3,%%xmm1					\n\t"\
		"subpd	%%xmm6,%%xmm5					\n\t"\
		"movaps	%%xmm0,     (%%rax)				\n\t"\
		"movaps	%%xmm4,     (%%rcx)				\n\t"\
		"movaps	%%xmm1,0x010(%%rax)				\n\t"\
		"movaps	%%xmm5,0x010(%%rdx)				\n\t"\
		"addpd	%%xmm2,%%xmm2					\n\t"\
		"addpd	%%xmm7,%%xmm7					\n\t"\
		"addpd	%%xmm3,%%xmm3					\n\t"\
		"addpd	%%xmm6,%%xmm6					\n\t"\
		"addpd	%%xmm0,%%xmm2					\n\t"\
		"addpd	%%xmm4,%%xmm7					\n\t"\
		"addpd	%%xmm1,%%xmm3					\n\t"\
		"addpd	%%xmm5,%%xmm6					\n\t"\
		"movaps	%%xmm2,     (%%rbx)				\n\t"\
		"movaps	%%xmm7,     (%%rdx)				\n\t"\
		"movaps	%%xmm3,0x010(%%rbx)				\n\t"\
		"movaps	%%xmm6,0x010(%%rcx)				\n\t"\
		"										\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__p01] "m" (Xp01)\
		 ,[__p04] "m" (Xp04)\
		 ,[__p08] "m" (Xp08)\
		 ,[__p16] "m" (Xp16)\
		 ,[__a0] "m" (Xa0)\
		 ,[__a1] "m" (Xa1)\
		 ,[__a2] "m" (Xa2)\
		 ,[__a3] "m" (Xa3)\
		 ,[__a4] "m" (Xa4)\
		 ,[__b0] "m" (Xb0)\
		 ,[__b1] "m" (Xb1)\
		 ,[__b2] "m" (Xb2)\
		 ,[__b3] "m" (Xb3)\
		 ,[__b4] "m" (Xb4)\
		 ,[__c0] "m" (Xc0)\
		 ,[__c1] "m" (Xc1)\
		 ,[__c2] "m" (Xc2)\
		 ,[__c3] "m" (Xc3)\
		 ,[__c4] "m" (Xc4)\
		 ,[__d0] "m" (Xd0)\
		 ,[__d1] "m" (Xd1)\
		 ,[__d2] "m" (Xd2)\
		 ,[__d3] "m" (Xd3)\
		 ,[__d4] "m" (Xd4)\
		 ,[__cc1] "m" (Xcc1)\
		 ,[__r00] "m" (Xr00)\
		 ,[__r10] "m" (Xr10)\
		 ,[__r20] "m" (Xr20)\
		 ,[__r30] "m" (Xr30)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r8","r9","r10","r11","r12","r13","r14","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
	);\
	}

  #endif	// AVX/SSE2 toggle

#endif	/* radix20_ditN_cy_dif1_gcc_h_included */

