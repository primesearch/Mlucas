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
#ifndef radix32_wrapper_square_gcc_h_included
#define radix32_wrapper_square_gcc_h_included

#ifdef USE_AVX

	#define SSE2_RADIX32_WRAPPER_DIF(Xadd0,Xadd1,Xadd2,Xadd3,Xr00,Xr10,Xr20,Xr30,Xisrt2,Xcc0,Xc00,Xc01,Xc02,Xc03,Xc05,Xc07)\
	{\
	__asm__ volatile (\
	"/*************************************************************/\n\t"\
	"/* SSE2_RADIX16_WRAPPER_DIF, 1st set of inputs:              */\n\t"\
	"/*************************************************************/\n\t"\
		"movq	%[__r00] ,%%rsi	\n\t"\
		"movq	%[__add0],%%rax\n\t"\
		"movq	%[__add1],%%rbx\n\t"\
		"movq	%[__add2],%%rcx\n\t"\
		"movq	%[__add3],%%rdx\n\t"\
	"/**** Start with 4-way interleaving: ****/\n\t"\
	"/* a[j+p0]: Inputs from [add0,add0+0x200,add1,add1+0x200]+0x0. Outputs into r0 +0/1, 16/17, 32/33, 48/49: */	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps	     (%%rax),%%ymm1						\n\t		vmovaps	0x020(%%rax),%%ymm3								\n\t"\
		"vmovaps	     (%%rbx),%%ymm5						\n\t		vmovaps	0x020(%%rbx),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps	     (%%rcx),%%ymm6						\n\t		vmovaps	0x020(%%rcx),%%ymm14							\n\t"\
		"vmovaps	     (%%rdx),%%ymm5						\n\t		vmovaps	0x020(%%rdx),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm6 ,%%ymm0				\n\t		vshufpd	$15 ,%%ymm13,%%ymm14 ,%%ymm2					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm6 ,%%ymm6				\n\t		vshufpd	$0  ,%%ymm13,%%ymm14 ,%%ymm14					\n\t"\
		"vperm2f128 $32 ,%%ymm0 ,%%ymm7 ,%%ymm5				\n\t		vperm2f128 $32 ,%%ymm2 ,%%ymm15,%%ymm13		/* outB	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm0 ,%%ymm7 ,%%ymm7				\n\t		vperm2f128 $49 ,%%ymm2 ,%%ymm15,%%ymm15		/* outD	*/	\n\t"\
		"vperm2f128 $32 ,%%ymm6 ,%%ymm1 ,%%ymm0				\n\t		vperm2f128 $32 ,%%ymm14,%%ymm3 ,%%ymm2 		/* outA	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm6 ,%%ymm1 ,%%ymm1				\n\t		vperm2f128 $49 ,%%ymm14,%%ymm3 ,%%ymm3 		/* outC	*/	\n\t"\
		"vmovaps 	%%ymm5 ,0x400(%%rsi)					\n\t		vmovaps %%ymm13,0x420(%%rsi)				/* outB	*/	\n\t"\
		"vmovaps 	%%ymm7 ,0x600(%%rsi)					\n\t		vmovaps %%ymm15,0x620(%%rsi)				/* outD	*/	\n\t"\
		"vmovaps 	%%ymm0 ,     (%%rsi)					\n\t		vmovaps %%ymm2 ,0x020(%%rsi)				/* outA	*/	\n\t"\
		"vmovaps 	%%ymm1 ,0x200(%%rsi)					\n\t		vmovaps %%ymm3 ,0x220(%%rsi)				/* outC	*/	\n\t"\
	"/* a[j+p2]: Inputs from [add0,add0+0x200,add1,add1+0x200]+0x40. Outputs into **r8** +0/1, 16/17, 32/33, 48/49: */	\n\t"\
		"addq	$0x40,%%rax	\n\t"\
		"addq	$0x40,%%rbx	\n\t"\
		"addq	$0x40,%%rcx	\n\t"\
		"addq	$0x40,%%rdx	\n\t"\
		"addq	$0x100,%%rsi	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps	     (%%rax),%%ymm1						\n\t		vmovaps	0x020(%%rax),%%ymm3								\n\t"\
		"vmovaps	     (%%rbx),%%ymm5						\n\t		vmovaps	0x020(%%rbx),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps	     (%%rcx),%%ymm6						\n\t		vmovaps	0x020(%%rcx),%%ymm14							\n\t"\
		"vmovaps	     (%%rdx),%%ymm5						\n\t		vmovaps	0x020(%%rdx),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm6 ,%%ymm0				\n\t		vshufpd	$15 ,%%ymm13,%%ymm14 ,%%ymm2					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm6 ,%%ymm6				\n\t		vshufpd	$0  ,%%ymm13,%%ymm14 ,%%ymm14					\n\t"\
		"vperm2f128 $32 ,%%ymm0 ,%%ymm7 ,%%ymm5				\n\t		vperm2f128 $32 ,%%ymm2 ,%%ymm15,%%ymm13		/* outB	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm0 ,%%ymm7 ,%%ymm7				\n\t		vperm2f128 $49 ,%%ymm2 ,%%ymm15,%%ymm15		/* outD	*/	\n\t"\
		"vperm2f128 $32 ,%%ymm6 ,%%ymm1 ,%%ymm0				\n\t		vperm2f128 $32 ,%%ymm14,%%ymm3 ,%%ymm2 		/* outA	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm6 ,%%ymm1 ,%%ymm1				\n\t		vperm2f128 $49 ,%%ymm14,%%ymm3 ,%%ymm3 		/* outC	*/	\n\t"\
		"vmovaps 	%%ymm5 ,0x400(%%rsi)					\n\t		vmovaps %%ymm13,0x420(%%rsi)				/* outB	*/	\n\t"\
		"vmovaps 	%%ymm7 ,0x600(%%rsi)					\n\t		vmovaps %%ymm15,0x620(%%rsi)				/* outD	*/	\n\t"\
		"vmovaps 	%%ymm0 ,     (%%rsi)					\n\t		vmovaps %%ymm2 ,0x020(%%rsi)				/* outA	*/	\n\t"\
		"vmovaps 	%%ymm1 ,0x200(%%rsi)					\n\t		vmovaps %%ymm3 ,0x220(%%rsi)				/* outC	*/	\n\t"\
	"/* a[j+p4]: Inputs from [add0,add0+0x200,add1,add1+0x200]+0x80. Outputs into r4 +0/1, 16/17, 32/33, 48/49: */	\n\t"\
		"addq	$0x40,%%rax	\n\t"\
		"addq	$0x40,%%rbx	\n\t"\
		"addq	$0x40,%%rcx	\n\t"\
		"addq	$0x40,%%rdx	\n\t"\
		"subq	$0x80,%%rsi	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps	     (%%rax),%%ymm1						\n\t		vmovaps	0x020(%%rax),%%ymm3								\n\t"\
		"vmovaps	     (%%rbx),%%ymm5						\n\t		vmovaps	0x020(%%rbx),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps	     (%%rcx),%%ymm6						\n\t		vmovaps	0x020(%%rcx),%%ymm14							\n\t"\
		"vmovaps	     (%%rdx),%%ymm5						\n\t		vmovaps	0x020(%%rdx),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm6 ,%%ymm0				\n\t		vshufpd	$15 ,%%ymm13,%%ymm14 ,%%ymm2					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm6 ,%%ymm6				\n\t		vshufpd	$0  ,%%ymm13,%%ymm14 ,%%ymm14					\n\t"\
		"vperm2f128 $32 ,%%ymm0 ,%%ymm7 ,%%ymm5				\n\t		vperm2f128 $32 ,%%ymm2 ,%%ymm15,%%ymm13		/* outB	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm0 ,%%ymm7 ,%%ymm7				\n\t		vperm2f128 $49 ,%%ymm2 ,%%ymm15,%%ymm15		/* outD	*/	\n\t"\
		"vperm2f128 $32 ,%%ymm6 ,%%ymm1 ,%%ymm0				\n\t		vperm2f128 $32 ,%%ymm14,%%ymm3 ,%%ymm2 		/* outA	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm6 ,%%ymm1 ,%%ymm1				\n\t		vperm2f128 $49 ,%%ymm14,%%ymm3 ,%%ymm3 		/* outC	*/	\n\t"\
		"vmovaps 	%%ymm5 ,0x400(%%rsi)					\n\t		vmovaps %%ymm13,0x420(%%rsi)				/* outB	*/	\n\t"\
		"vmovaps 	%%ymm7 ,0x600(%%rsi)					\n\t		vmovaps %%ymm15,0x620(%%rsi)				/* outD	*/	\n\t"\
		"vmovaps 	%%ymm0 ,     (%%rsi)					\n\t		vmovaps %%ymm2 ,0x020(%%rsi)				/* outA	*/	\n\t"\
		"vmovaps 	%%ymm1 ,0x200(%%rsi)					\n\t		vmovaps %%ymm3 ,0x220(%%rsi)				/* outC	*/	\n\t"\
	"/* a[j+p6]: Inputs from [add0,add0+0x200,add1,add1+0x200]+0xc0. Outputs into **r12** +0/1, 16/17, 32/33, 48/49: */	\n\t"\
		"addq	$0x40,%%rax	\n\t"\
		"addq	$0x40,%%rbx	\n\t"\
		"addq	$0x40,%%rcx	\n\t"\
		"addq	$0x40,%%rdx	\n\t"\
		"addq	$0x100,%%rsi	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps	     (%%rax),%%ymm1						\n\t		vmovaps	0x020(%%rax),%%ymm3								\n\t"\
		"vmovaps	     (%%rbx),%%ymm5						\n\t		vmovaps	0x020(%%rbx),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps	     (%%rcx),%%ymm6						\n\t		vmovaps	0x020(%%rcx),%%ymm14							\n\t"\
		"vmovaps	     (%%rdx),%%ymm5						\n\t		vmovaps	0x020(%%rdx),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm6 ,%%ymm0				\n\t		vshufpd	$15 ,%%ymm13,%%ymm14 ,%%ymm2					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm6 ,%%ymm6				\n\t		vshufpd	$0  ,%%ymm13,%%ymm14 ,%%ymm14					\n\t"\
		"vperm2f128 $32 ,%%ymm0 ,%%ymm7 ,%%ymm5				\n\t		vperm2f128 $32 ,%%ymm2 ,%%ymm15,%%ymm13		/* outB	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm0 ,%%ymm7 ,%%ymm7				\n\t		vperm2f128 $49 ,%%ymm2 ,%%ymm15,%%ymm15		/* outD	*/	\n\t"\
		"vperm2f128 $32 ,%%ymm6 ,%%ymm1 ,%%ymm0				\n\t		vperm2f128 $32 ,%%ymm14,%%ymm3 ,%%ymm2 		/* outA	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm6 ,%%ymm1 ,%%ymm1				\n\t		vperm2f128 $49 ,%%ymm14,%%ymm3 ,%%ymm3 		/* outC	*/	\n\t"\
		"vmovaps 	%%ymm5 ,0x400(%%rsi)					\n\t		vmovaps %%ymm13,0x420(%%rsi)				/* outB	*/	\n\t"\
		"vmovaps 	%%ymm7 ,0x600(%%rsi)					\n\t		vmovaps %%ymm15,0x620(%%rsi)				/* outD	*/	\n\t"\
		"vmovaps 	%%ymm0 ,     (%%rsi)					\n\t		vmovaps %%ymm2 ,0x020(%%rsi)				/* outA	*/	\n\t"\
		"vmovaps 	%%ymm1 ,0x200(%%rsi)					\n\t		vmovaps %%ymm3 ,0x220(%%rsi)				/* outC	*/	\n\t"\
	"/* a[j+p8]: Inputs from [add0,add0+0x200,add1,add1+0x200]+0x100. Outputs into **r2** +0/1, 16/17, 32/33, 48/49: */	\n\t"\
		"addq	$0x40,%%rax	\n\t"\
		"addq	$0x40,%%rbx	\n\t"\
		"addq	$0x40,%%rcx	\n\t"\
		"addq	$0x40,%%rdx	\n\t"\
		"subq	$0x140,%%rsi	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps	     (%%rax),%%ymm1						\n\t		vmovaps	0x020(%%rax),%%ymm3								\n\t"\
		"vmovaps	     (%%rbx),%%ymm5						\n\t		vmovaps	0x020(%%rbx),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps	     (%%rcx),%%ymm6						\n\t		vmovaps	0x020(%%rcx),%%ymm14							\n\t"\
		"vmovaps	     (%%rdx),%%ymm5						\n\t		vmovaps	0x020(%%rdx),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm6 ,%%ymm0				\n\t		vshufpd	$15 ,%%ymm13,%%ymm14 ,%%ymm2					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm6 ,%%ymm6				\n\t		vshufpd	$0  ,%%ymm13,%%ymm14 ,%%ymm14					\n\t"\
		"vperm2f128 $32 ,%%ymm0 ,%%ymm7 ,%%ymm5				\n\t		vperm2f128 $32 ,%%ymm2 ,%%ymm15,%%ymm13		/* outB	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm0 ,%%ymm7 ,%%ymm7				\n\t		vperm2f128 $49 ,%%ymm2 ,%%ymm15,%%ymm15		/* outD	*/	\n\t"\
		"vperm2f128 $32 ,%%ymm6 ,%%ymm1 ,%%ymm0				\n\t		vperm2f128 $32 ,%%ymm14,%%ymm3 ,%%ymm2 		/* outA	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm6 ,%%ymm1 ,%%ymm1				\n\t		vperm2f128 $49 ,%%ymm14,%%ymm3 ,%%ymm3 		/* outC	*/	\n\t"\
		"vmovaps 	%%ymm5 ,0x400(%%rsi)					\n\t		vmovaps %%ymm13,0x420(%%rsi)				/* outB	*/	\n\t"\
		"vmovaps 	%%ymm7 ,0x600(%%rsi)					\n\t		vmovaps %%ymm15,0x620(%%rsi)				/* outD	*/	\n\t"\
		"vmovaps 	%%ymm0 ,     (%%rsi)					\n\t		vmovaps %%ymm2 ,0x020(%%rsi)				/* outA	*/	\n\t"\
		"vmovaps 	%%ymm1 ,0x200(%%rsi)					\n\t		vmovaps %%ymm3 ,0x220(%%rsi)				/* outC	*/	\n\t"\
	"/* a[j+p10]: Inputs from [add0,add0+0x200,add1,add1+0x200]+0x140. Outputs into r10 +0/1, 16/17, 32/33, 48/49: */	\n\t"\
		"addq	$0x40,%%rax	\n\t"\
		"addq	$0x40,%%rbx	\n\t"\
		"addq	$0x40,%%rcx	\n\t"\
		"addq	$0x40,%%rdx	\n\t"\
		"addq	$0x100,%%rsi	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps	     (%%rax),%%ymm1						\n\t		vmovaps	0x020(%%rax),%%ymm3								\n\t"\
		"vmovaps	     (%%rbx),%%ymm5						\n\t		vmovaps	0x020(%%rbx),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps	     (%%rcx),%%ymm6						\n\t		vmovaps	0x020(%%rcx),%%ymm14							\n\t"\
		"vmovaps	     (%%rdx),%%ymm5						\n\t		vmovaps	0x020(%%rdx),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm6 ,%%ymm0				\n\t		vshufpd	$15 ,%%ymm13,%%ymm14 ,%%ymm2					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm6 ,%%ymm6				\n\t		vshufpd	$0  ,%%ymm13,%%ymm14 ,%%ymm14					\n\t"\
		"vperm2f128 $32 ,%%ymm0 ,%%ymm7 ,%%ymm5				\n\t		vperm2f128 $32 ,%%ymm2 ,%%ymm15,%%ymm13		/* outB	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm0 ,%%ymm7 ,%%ymm7				\n\t		vperm2f128 $49 ,%%ymm2 ,%%ymm15,%%ymm15		/* outD	*/	\n\t"\
		"vperm2f128 $32 ,%%ymm6 ,%%ymm1 ,%%ymm0				\n\t		vperm2f128 $32 ,%%ymm14,%%ymm3 ,%%ymm2 		/* outA	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm6 ,%%ymm1 ,%%ymm1				\n\t		vperm2f128 $49 ,%%ymm14,%%ymm3 ,%%ymm3 		/* outC	*/	\n\t"\
		"vmovaps 	%%ymm5 ,0x400(%%rsi)					\n\t		vmovaps %%ymm13,0x420(%%rsi)				/* outB	*/	\n\t"\
		"vmovaps 	%%ymm7 ,0x600(%%rsi)					\n\t		vmovaps %%ymm15,0x620(%%rsi)				/* outD	*/	\n\t"\
		"vmovaps 	%%ymm0 ,     (%%rsi)					\n\t		vmovaps %%ymm2 ,0x020(%%rsi)				/* outA	*/	\n\t"\
		"vmovaps 	%%ymm1 ,0x200(%%rsi)					\n\t		vmovaps %%ymm3 ,0x220(%%rsi)				/* outC	*/	\n\t"\
	"/* a[j+p12]: Inputs from [add0,add0+0x200,add1,add1+0x200]+0x180. Outputs into **r6** +0/1, 16/17, 32/33, 48/49: */	\n\t"\
		"addq	$0x40,%%rax	\n\t"\
		"addq	$0x40,%%rbx	\n\t"\
		"addq	$0x40,%%rcx	\n\t"\
		"addq	$0x40,%%rdx	\n\t"\
		"subq	$0x80,%%rsi	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps	     (%%rax),%%ymm1						\n\t		vmovaps	0x020(%%rax),%%ymm3								\n\t"\
		"vmovaps	     (%%rbx),%%ymm5						\n\t		vmovaps	0x020(%%rbx),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps	     (%%rcx),%%ymm6						\n\t		vmovaps	0x020(%%rcx),%%ymm14							\n\t"\
		"vmovaps	     (%%rdx),%%ymm5						\n\t		vmovaps	0x020(%%rdx),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm6 ,%%ymm0				\n\t		vshufpd	$15 ,%%ymm13,%%ymm14 ,%%ymm2					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm6 ,%%ymm6				\n\t		vshufpd	$0  ,%%ymm13,%%ymm14 ,%%ymm14					\n\t"\
		"vperm2f128 $32 ,%%ymm0 ,%%ymm7 ,%%ymm5				\n\t		vperm2f128 $32 ,%%ymm2 ,%%ymm15,%%ymm13		/* outB	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm0 ,%%ymm7 ,%%ymm7				\n\t		vperm2f128 $49 ,%%ymm2 ,%%ymm15,%%ymm15		/* outD	*/	\n\t"\
		"vperm2f128 $32 ,%%ymm6 ,%%ymm1 ,%%ymm0				\n\t		vperm2f128 $32 ,%%ymm14,%%ymm3 ,%%ymm2 		/* outA	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm6 ,%%ymm1 ,%%ymm1				\n\t		vperm2f128 $49 ,%%ymm14,%%ymm3 ,%%ymm3 		/* outC	*/	\n\t"\
		"vmovaps 	%%ymm5 ,0x400(%%rsi)					\n\t		vmovaps %%ymm13,0x420(%%rsi)				/* outB	*/	\n\t"\
		"vmovaps 	%%ymm7 ,0x600(%%rsi)					\n\t		vmovaps %%ymm15,0x620(%%rsi)				/* outD	*/	\n\t"\
		"vmovaps 	%%ymm0 ,     (%%rsi)					\n\t		vmovaps %%ymm2 ,0x020(%%rsi)				/* outA	*/	\n\t"\
		"vmovaps 	%%ymm1 ,0x200(%%rsi)					\n\t		vmovaps %%ymm3 ,0x220(%%rsi)				/* outC	*/	\n\t"\
	"/* a[j+p14]: Inputs from [add0,add0+0x200,add1,add1+0x200]+0x1c0. Outputs into r14 +0/1, 16/17, 32/33, 48/49: */	\n\t"\
		"addq	$0x40,%%rax	\n\t"\
		"addq	$0x40,%%rbx	\n\t"\
		"addq	$0x40,%%rcx	\n\t"\
		"addq	$0x40,%%rdx	\n\t"\
		"addq	$0x100,%%rsi	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps	     (%%rax),%%ymm1						\n\t		vmovaps	0x020(%%rax),%%ymm3								\n\t"\
		"vmovaps	     (%%rbx),%%ymm5						\n\t		vmovaps	0x020(%%rbx),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps	     (%%rcx),%%ymm6						\n\t		vmovaps	0x020(%%rcx),%%ymm14							\n\t"\
		"vmovaps	     (%%rdx),%%ymm5						\n\t		vmovaps	0x020(%%rdx),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm6 ,%%ymm0				\n\t		vshufpd	$15 ,%%ymm13,%%ymm14 ,%%ymm2					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm6 ,%%ymm6				\n\t		vshufpd	$0  ,%%ymm13,%%ymm14 ,%%ymm14					\n\t"\
		"vperm2f128 $32 ,%%ymm0 ,%%ymm7 ,%%ymm5				\n\t		vperm2f128 $32 ,%%ymm2 ,%%ymm15,%%ymm13		/* outB	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm0 ,%%ymm7 ,%%ymm7				\n\t		vperm2f128 $49 ,%%ymm2 ,%%ymm15,%%ymm15		/* outD	*/	\n\t"\
		"vperm2f128 $32 ,%%ymm6 ,%%ymm1 ,%%ymm0				\n\t		vperm2f128 $32 ,%%ymm14,%%ymm3 ,%%ymm2 		/* outA	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm6 ,%%ymm1 ,%%ymm1				\n\t		vperm2f128 $49 ,%%ymm14,%%ymm3 ,%%ymm3 		/* outC	*/	\n\t"\
		"vmovaps 	%%ymm5 ,0x400(%%rsi)					\n\t		vmovaps %%ymm13,0x420(%%rsi)				/* outB	*/	\n\t"\
		"vmovaps 	%%ymm7 ,0x600(%%rsi)					\n\t		vmovaps %%ymm15,0x620(%%rsi)				/* outD	*/	\n\t"\
		"vmovaps 	%%ymm0 ,     (%%rsi)					\n\t		vmovaps %%ymm2 ,0x020(%%rsi)				/* outA	*/	\n\t"\
		"vmovaps 	%%ymm1 ,0x200(%%rsi)					\n\t		vmovaps %%ymm3 ,0x220(%%rsi)				/* outC	*/	\n\t"\
		"\n\t"\
		"/************************************************************************/\n\t"\
		"/* Forward DIF radix-32 pass on the interleaved block1 and block2 data: */\n\t"\
		"/************************************************************************/\n\t"\
		"/*...Block 0:	*/									\n\t"\
		"movq	%[__isrt2],%%rsi							\n\t		leaq	0x900(%%rsi),%%rdi	/* two */	\n\t"\
		"/*****	SSE2_RADIX4_DIF_4WRAPPER(c00,c08,c10,c18,r00) ****	SSE2_RADIX4_DIF_4WRAPPER_2NDOFTWO(c04,c0C,c14,c1C,r08)	*****/\n\t"\
		"movq		%[__r00]	,%%rcx						\n\t		/*addq		$0x100,%%rcx // __r08 */	\n\t"\
		"movq		%[__c00]	,%%rdx						\n\t		/*addq		$0x100,%%rdx // __c04 */	\n\t"\
		"vmovaps		     (%%rcx),%%ymm0					\n\t		vmovaps		0x100(%%rcx),%%ymm8			\n\t"\
		"vmovaps		0x020(%%rcx),%%ymm1					\n\t		vmovaps		0x120(%%rcx),%%ymm9			\n\t"\
		"vmovaps		%%ymm0		,%%ymm2					\n\t		vmovaps		%%ymm8,%%ymm10				\n\t"\
		"vmovaps		%%ymm1		,%%ymm3					\n\t		vmovaps		%%ymm9,%%ymm11				\n\t"\
		"vmulpd		     (%%rdx),%%ymm0,%%ymm0				\n\t		vmulpd		0x100(%%rdx),%%ymm8	,%%ymm8			\n\t"\
		"vmulpd		     (%%rdx),%%ymm1,%%ymm1				\n\t		vmulpd		0x100(%%rdx),%%ymm9	,%%ymm9			\n\t"\
		"vmulpd		0x020(%%rdx),%%ymm2,%%ymm2				\n\t		vmulpd		0x120(%%rdx),%%ymm10,%%ymm10		\n\t"\
		"vmulpd		0x020(%%rdx),%%ymm3,%%ymm3				\n\t		vmulpd		0x120(%%rdx),%%ymm11,%%ymm11		\n\t"\
		"vaddpd		%%ymm2		,%%ymm1,%%ymm1				\n\t		vaddpd		%%ymm10		,%%ymm9	,%%ymm9				\n\t"\
		"vsubpd		%%ymm3		,%%ymm0,%%ymm0				\n\t		vsubpd		%%ymm11		,%%ymm8	,%%ymm8				\n\t"\
		"vmovaps		%%ymm0		,%%ymm2					\n\t		vmovaps		%%ymm8,%%ymm10						\n\t"\
		"vmovaps		%%ymm1		,%%ymm3					\n\t		vmovaps		%%ymm9,%%ymm11						\n\t"\
		"addq		$0x040		,%%rdx						\n\t"\
		"vmovaps		0x040(%%rcx),%%ymm4					\n\t		vmovaps		0x140(%%rcx),%%ymm12				\n\t"\
		"vmovaps		0x060(%%rcx),%%ymm5					\n\t		vmovaps		0x160(%%rcx),%%ymm13				\n\t"\
		"vmovaps		%%ymm4		,%%ymm6					\n\t		vmovaps		%%ymm12,%%ymm14						\n\t"\
		"vmovaps		%%ymm5		,%%ymm7					\n\t		vmovaps		%%ymm13,%%ymm15						\n\t"\
		"vmulpd		     (%%rdx),%%ymm4,%%ymm4				\n\t		vmulpd		0x100(%%rdx),%%ymm12,%%ymm12		\n\t"\
		"vmulpd		     (%%rdx),%%ymm5,%%ymm5				\n\t		vmulpd		0x100(%%rdx),%%ymm13,%%ymm13		\n\t"\
		"vmulpd		0x020(%%rdx),%%ymm6,%%ymm6				\n\t		vmulpd		0x120(%%rdx),%%ymm14,%%ymm14		\n\t"\
		"vmulpd		0x020(%%rdx),%%ymm7,%%ymm7				\n\t		vmulpd		0x120(%%rdx),%%ymm15,%%ymm15		\n\t"\
		"vaddpd		%%ymm6		,%%ymm5,%%ymm5				\n\t		vaddpd		%%ymm14		,%%ymm13,%%ymm13				\n\t"\
		"vsubpd		%%ymm7		,%%ymm4,%%ymm4				\n\t		vsubpd		%%ymm15		,%%ymm12,%%ymm12				\n\t"\
		"vaddpd		%%ymm4		,%%ymm0,%%ymm0				\n\t		vaddpd		%%ymm12		,%%ymm8	,%%ymm8				\n\t"\
		"vaddpd		%%ymm5		,%%ymm1,%%ymm1				\n\t		vaddpd		%%ymm13		,%%ymm9	,%%ymm9				\n\t"\
		"vsubpd		%%ymm4		,%%ymm2,%%ymm2				\n\t		vsubpd		%%ymm12		,%%ymm10,%%ymm10				\n\t"\
		"vsubpd		%%ymm5		,%%ymm3,%%ymm3				\n\t		vsubpd		%%ymm13		,%%ymm11,%%ymm11				\n\t"\
		"addq		$0x080		,%%rdx						\n\t"\
		"vmovaps		0x0c0(%%rcx),%%ymm4					\n\t		vmovaps		0x1c0(%%rcx),%%ymm12				\n\t"\
		"vmovaps		0x0e0(%%rcx),%%ymm5					\n\t		vmovaps		0x1e0(%%rcx),%%ymm13				\n\t"\
		"vmovaps		%%ymm4		,%%ymm6					\n\t		vmovaps		%%ymm12,%%ymm14						\n\t"\
		"vmovaps		%%ymm5		,%%ymm7					\n\t		vmovaps		%%ymm13,%%ymm15						\n\t"\
		"vmulpd		     (%%rdx),%%ymm4,%%ymm4				\n\t		vmulpd		0x100(%%rdx),%%ymm12,%%ymm12		\n\t"\
		"vmulpd		     (%%rdx),%%ymm5,%%ymm5				\n\t		vmulpd		0x100(%%rdx),%%ymm13,%%ymm13		\n\t"\
		"vmulpd		0x020(%%rdx),%%ymm6,%%ymm6				\n\t		vmulpd		0x120(%%rdx),%%ymm14,%%ymm14		\n\t"\
		"vmulpd		0x020(%%rdx),%%ymm7,%%ymm7				\n\t		vmulpd		0x120(%%rdx),%%ymm15,%%ymm15		\n\t"\
		"vaddpd		%%ymm6		,%%ymm5,%%ymm5				\n\t		vaddpd		%%ymm14		,%%ymm13,%%ymm13				\n\t"\
		"vsubpd		%%ymm7		,%%ymm4,%%ymm4				\n\t		vsubpd		%%ymm15		,%%ymm12,%%ymm12				\n\t"\
		"vmovaps		%%ymm5		,0x020(%%rcx)			\n\t		vmovaps		%%ymm13,0x120(%%rcx)				\n\t"\
		"vmovaps		%%ymm4		,     (%%rcx)			\n\t		vmovaps		%%ymm12,0x100(%%rcx)				\n\t"\
		"subq		$0x040		,%%rdx						\n\t"\
		"vmovaps		0x080(%%rcx),%%ymm4					\n\t		vmovaps		0x180(%%rcx),%%ymm12				\n\t"\
		"vmovaps		0x0a0(%%rcx),%%ymm5					\n\t		vmovaps		0x1a0(%%rcx),%%ymm13				\n\t"\
		"vmovaps		%%ymm4		,%%ymm6					\n\t		vmovaps		%%ymm12,%%ymm14						\n\t"\
		"vmovaps		%%ymm5		,%%ymm7					\n\t		vmovaps		%%ymm13,%%ymm15						\n\t"\
		"vmulpd		     (%%rdx),%%ymm4,%%ymm4				\n\t		vmulpd		0x100(%%rdx),%%ymm12,%%ymm12		\n\t"\
		"vmulpd		     (%%rdx),%%ymm5,%%ymm5				\n\t		vmulpd		0x100(%%rdx),%%ymm13,%%ymm13		\n\t"\
		"vmulpd		0x020(%%rdx),%%ymm6,%%ymm6				\n\t		vmulpd		0x120(%%rdx),%%ymm14,%%ymm14		\n\t"\
		"vmulpd		0x020(%%rdx),%%ymm7,%%ymm7				\n\t		vmulpd		0x120(%%rdx),%%ymm15,%%ymm15		\n\t"\
		"vaddpd		%%ymm6		,%%ymm5,%%ymm5				\n\t		vaddpd		%%ymm14		,%%ymm13,%%ymm13				\n\t"\
		"vsubpd		%%ymm7		,%%ymm4,%%ymm4				\n\t		vsubpd		%%ymm15		,%%ymm12,%%ymm12				\n\t"\
		"vmovaps		%%ymm5		,%%ymm7					\n\t		vmovaps		%%ymm13,%%ymm15						\n\t"\
		"vmovaps		%%ymm4		,%%ymm6					\n\t		vmovaps		%%ymm12,%%ymm14						\n\t"\
		"vsubpd		     (%%rcx),%%ymm4,%%ymm4				\n\t		vsubpd		0x100(%%rcx),%%ymm12,%%ymm12		\n\t"\
		"vsubpd		0x020(%%rcx),%%ymm5,%%ymm5				\n\t		vsubpd		0x120(%%rcx),%%ymm13,%%ymm13		\n\t"\
		"vaddpd		     (%%rcx),%%ymm6,%%ymm6				\n\t		vaddpd		0x100(%%rcx),%%ymm14,%%ymm14		\n\t"\
		"vaddpd		0x020(%%rcx),%%ymm7,%%ymm7				\n\t		vaddpd		0x120(%%rcx),%%ymm15,%%ymm15		\n\t"\
		"vsubpd		%%ymm6		,%%ymm0,%%ymm0				\n\t		vsubpd		%%ymm14		,%%ymm8	,%%ymm8				\n\t"\
		"vsubpd		%%ymm7		,%%ymm1,%%ymm1				\n\t		vsubpd		%%ymm15		,%%ymm9	,%%ymm9				\n\t"\
		"vmulpd		(%%rdi)		,%%ymm6,%%ymm6				\n\t		vmulpd		(%%rdi)		,%%ymm14,%%ymm14				\n\t"\
		"vmulpd		(%%rdi)		,%%ymm7,%%ymm7				\n\t		vmulpd		(%%rdi)		,%%ymm15,%%ymm15				\n\t"\
		"vaddpd		%%ymm0		,%%ymm6,%%ymm6				\n\t		vaddpd		%%ymm8		,%%ymm14,%%ymm14				\n\t"\
		"vaddpd		%%ymm1		,%%ymm7,%%ymm7				\n\t		vaddpd		%%ymm9		,%%ymm15,%%ymm15				\n\t"\
		"/*vmovaps	%%ymm6		,     (%%rcx)	*/			\n\t		vmovaps		%%ymm14,0x100(%%rcx)				\n\t"\
		"/*vmovaps	%%ymm7		,0x020(%%rcx)	*/			\n\t		vmovaps		%%ymm15,0x120(%%rcx)				\n\t"\
		"vsubpd		%%ymm5		,%%ymm2,%%ymm2				\n\t		vsubpd		%%ymm13		,%%ymm10,%%ymm10				\n\t"\
		"vsubpd		%%ymm4		,%%ymm3,%%ymm3				\n\t		vsubpd		%%ymm12		,%%ymm11,%%ymm11				\n\t"\
		"/*vmovaps	%%ymm2		,0x040(%%rcx)	*/			\n\t		vmulpd		(%%rdi)		,%%ymm13,%%ymm13				\n\t"\
		"/*vmovaps	%%ymm3		,0x0e0(%%rcx)	*/			\n\t		vmulpd		(%%rdi)		,%%ymm12,%%ymm12				\n\t"\
		"vmulpd		(%%rdi)		,%%ymm5,%%ymm5				\n\t		vaddpd		%%ymm10		,%%ymm13,%%ymm13				\n\t"\
		"vmulpd		(%%rdi)		,%%ymm4,%%ymm4				\n\t		vaddpd		%%ymm11		,%%ymm12,%%ymm12				\n\t"\
		"vaddpd		%%ymm2		,%%ymm5,%%ymm5				\n\t		vmovaps		%%ymm10,%%ymm14						\n\t"\
		"vaddpd		%%ymm3		,%%ymm4,%%ymm4				\n\t		vmovaps		%%ymm13,%%ymm15						\n\t"\
		"/*vmovaps	%%ymm5		,0x0c0(%%rcx)	*/			\n\t		vsubpd		%%ymm12		,%%ymm10,%%ymm10				\n\t"\
		"/*vmovaps	%%ymm4		,0x060(%%rcx)	*/						vsubpd		%%ymm11		,%%ymm13,%%ymm13				\n\t"\
		"																vaddpd		%%ymm12		,%%ymm14,%%ymm14				\n\t"\
		"																vaddpd		%%ymm11		,%%ymm15,%%ymm15				\n\t"\
		"																vmulpd		(%%rsi)		,%%ymm10,%%ymm10				\n\t"\
		"																vmulpd		(%%rsi)		,%%ymm13,%%ymm13				\n\t"\
		"																vmulpd		(%%rsi)		,%%ymm14,%%ymm14				\n\t"\
		"																vmulpd		(%%rsi)		,%%ymm15,%%ymm15				\n\t"\
		"/*****	SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r00) *****/\n\t"\
		"																vmovaps		0x100(%%rcx),%%ymm11		\n\t"\
		"																vmovaps		0x120(%%rcx),%%ymm12		\n\t"\
		"vsubpd		%%ymm11		,%%ymm6	,%%ymm6				\n\t		vsubpd		%%ymm10		,%%ymm2	,%%ymm2			\n\t"\
		"vsubpd		%%ymm9		,%%ymm0	,%%ymm0				\n\t		vsubpd		%%ymm15		,%%ymm5	,%%ymm5			\n\t"\
		"vsubpd		%%ymm12		,%%ymm7	,%%ymm7				\n\t		vsubpd		%%ymm14		,%%ymm4	,%%ymm4			\n\t"\
		"vsubpd		%%ymm8		,%%ymm1	,%%ymm1				\n\t		vsubpd		%%ymm13		,%%ymm3	,%%ymm3			\n\t"\
		"vmulpd		(%%rdi)		,%%ymm11,%%ymm11			\n\t		vmulpd		(%%rdi)		,%%ymm10,%%ymm10		\n\t"\
		"vmulpd		(%%rdi)		,%%ymm9	,%%ymm9				\n\t		vmulpd		(%%rdi)		,%%ymm15,%%ymm15		\n\t"\
		"vmulpd		(%%rdi)		,%%ymm12,%%ymm12			\n\t		vmulpd		(%%rdi)		,%%ymm14,%%ymm14		\n\t"\
		"vmulpd		(%%rdi)		,%%ymm8	,%%ymm8				\n\t		vmulpd		(%%rdi)		,%%ymm13,%%ymm13		\n\t"\
		"vaddpd		%%ymm6		,%%ymm11,%%ymm11			\n\t		vaddpd		%%ymm2		,%%ymm10,%%ymm10		\n\t"\
		"vaddpd		%%ymm0		,%%ymm9	,%%ymm9				\n\t		vaddpd		%%ymm5		,%%ymm15,%%ymm15		\n\t"\
		"vaddpd		%%ymm7		,%%ymm12,%%ymm12			\n\t		vaddpd		%%ymm4		,%%ymm14,%%ymm14		\n\t"\
		"vaddpd		%%ymm1		,%%ymm8	,%%ymm8				\n\t		vaddpd		%%ymm3		,%%ymm13,%%ymm13		\n\t"\
		"vmovaps		%%ymm6		,0x100(%%rcx)			\n\t		vmovaps		%%ymm2		,0x140(%%rcx)	\n\t"\
		"vmovaps		%%ymm0		,0x080(%%rcx)			\n\t		vmovaps		%%ymm5		,0x0c0(%%rcx)	\n\t"\
		"vmovaps		%%ymm7		,0x120(%%rcx)			\n\t		vmovaps		%%ymm4		,0x160(%%rcx)	\n\t"\
		"vmovaps		%%ymm1		,0x1a0(%%rcx)			\n\t		vmovaps		%%ymm3		,0x1e0(%%rcx)	\n\t"\
		"vmovaps		%%ymm11		,     (%%rcx)			\n\t		vmovaps		%%ymm10		,0x040(%%rcx)	\n\t"\
		"vmovaps		%%ymm9		,0x180(%%rcx)			\n\t		vmovaps		%%ymm15		,0x1c0(%%rcx)	\n\t"\
		"vmovaps		%%ymm12		,0x020(%%rcx)			\n\t		vmovaps		%%ymm14		,0x060(%%rcx)	\n\t"\
		"vmovaps		%%ymm8		,0x0a0(%%rcx)			\n\t		vmovaps		%%ymm13		,0x0e0(%%rcx)	\n\t"\
		"\n\t"\
		"/*...Block 2:	*/\n\t"\
		"/*****	SSE2_RADIX4_DIF_4WRAPPER(c02,c0A,c12,c1A,r10) ****	SSE2_RADIX4_DIF_4WRAPPER_2NDOFTWO(c06,c0E,c16,c1E,r18)	*****/\n\t"\
		"movq		%[__r10]	,%%rcx						\n\t		/*addq		$0x100,%%rcx // __r18 */	\n\t"\
		"movq		%[__c02]	,%%rdx						\n\t		/*addq		$0x100,%%rdx // __c06 */	\n\t"\
		"vmovaps		     (%%rcx),%%ymm0					\n\t		vmovaps		0x100(%%rcx),%%ymm8			\n\t"\
		"vmovaps		0x020(%%rcx),%%ymm1					\n\t		vmovaps		0x120(%%rcx),%%ymm9			\n\t"\
		"vmovaps		%%ymm0		,%%ymm2					\n\t		vmovaps		%%ymm8,%%ymm10				\n\t"\
		"vmovaps		%%ymm1		,%%ymm3					\n\t		vmovaps		%%ymm9,%%ymm11				\n\t"\
		"vmulpd		     (%%rdx),%%ymm0,%%ymm0				\n\t		vmulpd		0x100(%%rdx),%%ymm8	,%%ymm8			\n\t"\
		"vmulpd		     (%%rdx),%%ymm1,%%ymm1				\n\t		vmulpd		0x100(%%rdx),%%ymm9	,%%ymm9			\n\t"\
		"vmulpd		0x020(%%rdx),%%ymm2,%%ymm2				\n\t		vmulpd		0x120(%%rdx),%%ymm10,%%ymm10		\n\t"\
		"vmulpd		0x020(%%rdx),%%ymm3,%%ymm3				\n\t		vmulpd		0x120(%%rdx),%%ymm11,%%ymm11		\n\t"\
		"vaddpd		%%ymm2		,%%ymm1,%%ymm1				\n\t		vaddpd		%%ymm10		,%%ymm9	,%%ymm9				\n\t"\
		"vsubpd		%%ymm3		,%%ymm0,%%ymm0				\n\t		vsubpd		%%ymm11		,%%ymm8	,%%ymm8				\n\t"\
		"vmovaps		%%ymm0		,%%ymm2					\n\t		vmovaps		%%ymm8,%%ymm10						\n\t"\
		"vmovaps		%%ymm1		,%%ymm3					\n\t		vmovaps		%%ymm9,%%ymm11						\n\t"\
		"addq		$0x040		,%%rdx						\n\t"\
		"vmovaps		0x040(%%rcx),%%ymm4					\n\t		vmovaps		0x140(%%rcx),%%ymm12				\n\t"\
		"vmovaps		0x060(%%rcx),%%ymm5					\n\t		vmovaps		0x160(%%rcx),%%ymm13				\n\t"\
		"vmovaps		%%ymm4		,%%ymm6					\n\t		vmovaps		%%ymm12,%%ymm14						\n\t"\
		"vmovaps		%%ymm5		,%%ymm7					\n\t		vmovaps		%%ymm13,%%ymm15						\n\t"\
		"vmulpd		     (%%rdx),%%ymm4,%%ymm4				\n\t		vmulpd		0x100(%%rdx),%%ymm12,%%ymm12		\n\t"\
		"vmulpd		     (%%rdx),%%ymm5,%%ymm5				\n\t		vmulpd		0x100(%%rdx),%%ymm13,%%ymm13		\n\t"\
		"vmulpd		0x020(%%rdx),%%ymm6,%%ymm6				\n\t		vmulpd		0x120(%%rdx),%%ymm14,%%ymm14		\n\t"\
		"vmulpd		0x020(%%rdx),%%ymm7,%%ymm7				\n\t		vmulpd		0x120(%%rdx),%%ymm15,%%ymm15		\n\t"\
		"vaddpd		%%ymm6		,%%ymm5,%%ymm5				\n\t		vaddpd		%%ymm14		,%%ymm13,%%ymm13				\n\t"\
		"vsubpd		%%ymm7		,%%ymm4,%%ymm4				\n\t		vsubpd		%%ymm15		,%%ymm12,%%ymm12				\n\t"\
		"vaddpd		%%ymm4		,%%ymm0,%%ymm0				\n\t		vaddpd		%%ymm12		,%%ymm8	,%%ymm8				\n\t"\
		"vaddpd		%%ymm5		,%%ymm1,%%ymm1				\n\t		vaddpd		%%ymm13		,%%ymm9	,%%ymm9				\n\t"\
		"vsubpd		%%ymm4		,%%ymm2,%%ymm2				\n\t		vsubpd		%%ymm12		,%%ymm10,%%ymm10				\n\t"\
		"vsubpd		%%ymm5		,%%ymm3,%%ymm3				\n\t		vsubpd		%%ymm13		,%%ymm11,%%ymm11				\n\t"\
		"addq		$0x080		,%%rdx						\n\t"\
		"vmovaps		0x0c0(%%rcx),%%ymm4					\n\t		vmovaps		0x1c0(%%rcx),%%ymm12				\n\t"\
		"vmovaps		0x0e0(%%rcx),%%ymm5					\n\t		vmovaps		0x1e0(%%rcx),%%ymm13				\n\t"\
		"vmovaps		%%ymm4		,%%ymm6					\n\t		vmovaps		%%ymm12,%%ymm14						\n\t"\
		"vmovaps		%%ymm5		,%%ymm7					\n\t		vmovaps		%%ymm13,%%ymm15						\n\t"\
		"vmulpd		     (%%rdx),%%ymm4,%%ymm4				\n\t		vmulpd		0x100(%%rdx),%%ymm12,%%ymm12		\n\t"\
		"vmulpd		     (%%rdx),%%ymm5,%%ymm5				\n\t		vmulpd		0x100(%%rdx),%%ymm13,%%ymm13		\n\t"\
		"vmulpd		0x020(%%rdx),%%ymm6,%%ymm6				\n\t		vmulpd		0x120(%%rdx),%%ymm14,%%ymm14		\n\t"\
		"vmulpd		0x020(%%rdx),%%ymm7,%%ymm7				\n\t		vmulpd		0x120(%%rdx),%%ymm15,%%ymm15		\n\t"\
		"vaddpd		%%ymm6		,%%ymm5,%%ymm5				\n\t		vaddpd		%%ymm14		,%%ymm13,%%ymm13				\n\t"\
		"vsubpd		%%ymm7		,%%ymm4,%%ymm4				\n\t		vsubpd		%%ymm15		,%%ymm12,%%ymm12				\n\t"\
		"vmovaps		%%ymm5		,0x020(%%rcx)			\n\t		vmovaps		%%ymm13,0x120(%%rcx)				\n\t"\
		"vmovaps		%%ymm4		,     (%%rcx)			\n\t		vmovaps		%%ymm12,0x100(%%rcx)				\n\t"\
		"subq		$0x040		,%%rdx						\n\t"\
		"vmovaps		0x080(%%rcx),%%ymm4					\n\t		vmovaps		0x180(%%rcx),%%ymm12				\n\t"\
		"vmovaps		0x0a0(%%rcx),%%ymm5					\n\t		vmovaps		0x1a0(%%rcx),%%ymm13				\n\t"\
		"vmovaps		%%ymm4		,%%ymm6					\n\t		vmovaps		%%ymm12,%%ymm14						\n\t"\
		"vmovaps		%%ymm5		,%%ymm7					\n\t		vmovaps		%%ymm13,%%ymm15						\n\t"\
		"vmulpd		     (%%rdx),%%ymm4,%%ymm4				\n\t		vmulpd		0x100(%%rdx),%%ymm12,%%ymm12		\n\t"\
		"vmulpd		     (%%rdx),%%ymm5,%%ymm5				\n\t		vmulpd		0x100(%%rdx),%%ymm13,%%ymm13		\n\t"\
		"vmulpd		0x020(%%rdx),%%ymm6,%%ymm6				\n\t		vmulpd		0x120(%%rdx),%%ymm14,%%ymm14		\n\t"\
		"vmulpd		0x020(%%rdx),%%ymm7,%%ymm7				\n\t		vmulpd		0x120(%%rdx),%%ymm15,%%ymm15		\n\t"\
		"vaddpd		%%ymm6		,%%ymm5,%%ymm5				\n\t		vaddpd		%%ymm14		,%%ymm13,%%ymm13				\n\t"\
		"vsubpd		%%ymm7		,%%ymm4,%%ymm4				\n\t		vsubpd		%%ymm15		,%%ymm12,%%ymm12				\n\t"\
		"vmovaps		%%ymm5		,%%ymm7					\n\t		vmovaps		%%ymm13,%%ymm15						\n\t"\
		"vmovaps		%%ymm4		,%%ymm6					\n\t		vmovaps		%%ymm12,%%ymm14						\n\t"\
		"vsubpd		     (%%rcx),%%ymm4,%%ymm4				\n\t		vsubpd		0x100(%%rcx),%%ymm12,%%ymm12		\n\t"\
		"vsubpd		0x020(%%rcx),%%ymm5,%%ymm5				\n\t		vsubpd		0x120(%%rcx),%%ymm13,%%ymm13		\n\t"\
		"vaddpd		     (%%rcx),%%ymm6,%%ymm6				\n\t		vaddpd		0x100(%%rcx),%%ymm14,%%ymm14		\n\t"\
		"vaddpd		0x020(%%rcx),%%ymm7,%%ymm7				\n\t		vaddpd		0x120(%%rcx),%%ymm15,%%ymm15		\n\t"\
		"vsubpd		%%ymm6		,%%ymm0,%%ymm0				\n\t		vsubpd		%%ymm14		,%%ymm8	,%%ymm8				\n\t"\
		"vsubpd		%%ymm7		,%%ymm1,%%ymm1				\n\t		vsubpd		%%ymm15		,%%ymm9	,%%ymm9				\n\t"\
		"vmulpd		(%%rdi)		,%%ymm6,%%ymm6				\n\t		vmulpd		(%%rdi)		,%%ymm14,%%ymm14				\n\t"\
		"vmulpd		(%%rdi)		,%%ymm7,%%ymm7				\n\t		vmulpd		(%%rdi)		,%%ymm15,%%ymm15				\n\t"\
		"vaddpd		%%ymm0		,%%ymm6,%%ymm6				\n\t		vaddpd		%%ymm8 		,%%ymm14,%%ymm14				\n\t"\
		"vaddpd		%%ymm1		,%%ymm7,%%ymm7				\n\t		vaddpd		%%ymm9 		,%%ymm15,%%ymm15				\n\t"\
		"/*vmovaps	%%ymm6		,     (%%rcx)	*/			\n\t		vmovaps		%%ymm14,0x100(%%rcx)				\n\t"\
		"/*vmovaps	%%ymm7		,0x020(%%rcx)	*/			\n\t		vmovaps		%%ymm15,0x120(%%rcx)				\n\t"\
		"vsubpd		%%ymm5		,%%ymm2,%%ymm2				\n\t		vsubpd		%%ymm13		,%%ymm10,%%ymm10				\n\t"\
		"vsubpd		%%ymm4		,%%ymm3,%%ymm3				\n\t		vsubpd		%%ymm12		,%%ymm11,%%ymm11				\n\t"\
		"/*vmovaps	%%ymm2		,0x040(%%rcx)	*/			\n\t		vmulpd		(%%rdi)		,%%ymm13,%%ymm13				\n\t"\
		"/*vmovaps	%%ymm3		,0x0e0(%%rcx)	*/			\n\t		vmulpd		(%%rdi)		,%%ymm12,%%ymm12				\n\t"\
		"vmulpd		(%%rdi)		,%%ymm5,%%ymm5				\n\t		vaddpd		%%ymm10		,%%ymm13,%%ymm13				\n\t"\
		"vmulpd		(%%rdi)		,%%ymm4,%%ymm4				\n\t		vaddpd		%%ymm11		,%%ymm12,%%ymm12				\n\t"\
		"vaddpd		%%ymm2		,%%ymm5,%%ymm5				\n\t		vmovaps		%%ymm10,%%ymm14						\n\t"\
		"vaddpd		%%ymm3		,%%ymm4,%%ymm4				\n\t		vmovaps		%%ymm13,%%ymm15						\n\t"\
		"/*vmovaps	%%ymm5		,0x0c0(%%rcx)	*/			\n\t		vsubpd		%%ymm12		,%%ymm10,%%ymm10				\n\t"\
		"/*vmovaps	%%ymm4		,0x060(%%rcx)	*/			\n\t		vsubpd		%%ymm11		,%%ymm13,%%ymm13				\n\t"\
		"																vaddpd		%%ymm12		,%%ymm14,%%ymm14				\n\t"\
		"																vaddpd		%%ymm11		,%%ymm15,%%ymm15				\n\t"\
		"																vmulpd		(%%rsi)		,%%ymm10,%%ymm10	/* isrt2 */	\n\t"\
		"																vmulpd		(%%rsi)		,%%ymm13,%%ymm13				\n\t"\
		"																vmulpd		(%%rsi)		,%%ymm14,%%ymm14				\n\t"\
		"																vmulpd		(%%rsi)		,%%ymm15,%%ymm15				\n\t"\
		"/*****	SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r10) *****/\n\t"\
		"																vmovaps		0x100(%%rcx),%%ymm11		\n\t"\
		"																vmovaps		0x120(%%rcx),%%ymm12		\n\t"\
		"vsubpd		%%ymm11		,%%ymm6	,%%ymm6				\n\t		vsubpd		%%ymm10		,%%ymm2	,%%ymm2			\n\t"\
		"vsubpd		%%ymm9		,%%ymm0	,%%ymm0				\n\t		vsubpd		%%ymm15		,%%ymm5	,%%ymm5			\n\t"\
		"vsubpd		%%ymm12		,%%ymm7	,%%ymm7				\n\t		vsubpd		%%ymm14		,%%ymm4	,%%ymm4			\n\t"\
		"vsubpd		%%ymm8		,%%ymm1	,%%ymm1				\n\t		vsubpd		%%ymm13		,%%ymm3	,%%ymm3			\n\t"\
		"vmulpd		(%%rdi)		,%%ymm11,%%ymm11			\n\t		vmulpd		(%%rdi)		,%%ymm10,%%ymm10		\n\t"\
		"vmulpd		(%%rdi)		,%%ymm9	,%%ymm9				\n\t		vmulpd		(%%rdi)		,%%ymm15,%%ymm15		\n\t"\
		"vmulpd		(%%rdi)		,%%ymm12,%%ymm12			\n\t		vmulpd		(%%rdi)		,%%ymm14,%%ymm14		\n\t"\
		"vmulpd		(%%rdi)		,%%ymm8	,%%ymm8				\n\t		vmulpd		(%%rdi)		,%%ymm13,%%ymm13		\n\t"\
		"vaddpd		%%ymm6		,%%ymm11,%%ymm11			\n\t		vaddpd		%%ymm2		,%%ymm10,%%ymm10		\n\t"\
		"vaddpd		%%ymm0		,%%ymm9	,%%ymm9				\n\t		vaddpd		%%ymm5		,%%ymm15,%%ymm15		\n\t"\
		"vaddpd		%%ymm7		,%%ymm12,%%ymm12			\n\t		vaddpd		%%ymm4		,%%ymm14,%%ymm14		\n\t"\
		"vaddpd		%%ymm1		,%%ymm8	,%%ymm8				\n\t		vaddpd		%%ymm3		,%%ymm13,%%ymm13		\n\t"\
		"vmovaps		%%ymm6		,0x100(%%rcx)			\n\t		vmovaps		%%ymm2		,0x140(%%rcx)	\n\t"\
		"vmovaps		%%ymm0		,0x080(%%rcx)			\n\t		vmovaps		%%ymm5		,0x0c0(%%rcx)	\n\t"\
		"vmovaps		%%ymm7		,0x120(%%rcx)			\n\t		vmovaps		%%ymm4		,0x160(%%rcx)	\n\t"\
		"vmovaps		%%ymm1		,0x1a0(%%rcx)			\n\t		vmovaps		%%ymm3		,0x1e0(%%rcx)	\n\t"\
		"vmovaps		%%ymm11		,     (%%rcx)			\n\t		vmovaps		%%ymm10		,0x040(%%rcx)	\n\t"\
		"vmovaps		%%ymm9		,0x180(%%rcx)			\n\t		vmovaps		%%ymm15		,0x1c0(%%rcx)	\n\t"\
		"vmovaps		%%ymm12		,0x020(%%rcx)			\n\t		vmovaps		%%ymm14		,0x060(%%rcx)	\n\t"\
		"vmovaps		%%ymm8		,0x0a0(%%rcx)			\n\t		vmovaps		%%ymm13		,0x0e0(%%rcx)	\n\t"\
		"/********************************************************************************************************\n\t"\
		" Next 2 blocks operate on odd-indexed elements from the unpck*pd commands which we stored to temporaries:\n\t"\
		"********************************************************************************************************/\n\t"\
		"/*...Block 3:	*/\n\t"\
		"/*	SSE2_RADIX4_DIF_4TWIDDLE         (r20,r24,r22,r26,r20,c01) */\n\t"\
		"addq		$0x200		,%%rcx		/*****	SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO(r28,r2C,r2A,r2E,r28,c05)	*****/\n\t"\
		"movq		%[__c01]	,%%rbx						\n\t	/*	movq		%[__c05]	,%%rbx	*/		\n\t"\
		"movq		%%rcx		,%%rax						\n\t	/*	addq		$0x080		,%%rax	*/		\n\t"\
		"addq		$0x040		,%%rcx						\n\t	/*	addq		$0x080		,%%rcx	*/		\n\t"\
		"vmovaps	     (%%rax),%%ymm0	\n\t movq %%rax,%%rdx \n\t	vmovaps		0x100(%%rax),%%ymm8			\n\t"\
		"vmovaps	     (%%rcx),%%ymm4						\n\t		vmovaps		0x100(%%rcx),%%ymm12		\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1						\n\t		vmovaps		0x120(%%rax),%%ymm9			\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5						\n\t		vmovaps		0x120(%%rcx),%%ymm13		\n\t"\
		"vmovaps	     (%%rbx),%%ymm6						\n\t		vmovaps		0x100(%%rbx),%%ymm14		\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm7						\n\t		vmovaps		0x120(%%rbx),%%ymm15		\n\t"\
		"vmovaps	%%ymm0		,%%ymm2						\n\t		vmovaps		%%ymm8		,%%ymm10		\n\t"\
		"vmovaps	%%ymm1		,%%ymm3						\n\t		vmovaps		%%ymm9		,%%ymm11		\n\t"\
		"vmulpd		%%ymm6		,%%ymm0,%%ymm0				\n\t		vmulpd		%%ymm14		,%%ymm8	,%%ymm8			\n\t"\
		"vmulpd		%%ymm6		,%%ymm1,%%ymm1				\n\t		vmulpd		%%ymm14		,%%ymm9	,%%ymm9			\n\t"\
		"vmulpd		%%ymm7		,%%ymm2,%%ymm2				\n\t		vmulpd		%%ymm15		,%%ymm10,%%ymm10		\n\t"\
		"vmulpd		%%ymm7		,%%ymm3,%%ymm3				\n\t		vmulpd		%%ymm15		,%%ymm11,%%ymm11		\n\t"\
		"vmovaps		%%ymm4		,%%ymm6					\n\t		vmovaps		%%ymm12		,%%ymm14				\n\t"\
		"vaddpd		%%ymm2		,%%ymm1,%%ymm1				\n\t		vaddpd		%%ymm10		,%%ymm9	,%%ymm9			\n\t"\
		"vmovaps		%%ymm5		,%%ymm7					\n\t		vmovaps		%%ymm13		,%%ymm15				\n\t"\
		"vmulpd		0x040(%%rbx),%%ymm4,%%ymm4				\n\t		vmulpd		0x140(%%rbx),%%ymm12,%%ymm12		\n\t"\
		"vsubpd		%%ymm3		,%%ymm0,%%ymm0				\n\t		vsubpd		%%ymm11		,%%ymm8	,%%ymm8			\n\t"\
		"vmulpd		0x040(%%rbx),%%ymm5,%%ymm5				\n\t		vmulpd		0x140(%%rbx),%%ymm13,%%ymm13		\n\t"\
		"vmulpd		0x060(%%rbx),%%ymm6,%%ymm6				\n\t		vmulpd		0x160(%%rbx),%%ymm14,%%ymm14		\n\t"\
		"vmovaps		%%ymm0		,%%ymm2					\n\t		vmovaps		%%ymm8		,%%ymm10				\n\t"\
		"vmulpd		0x060(%%rbx),%%ymm7,%%ymm7				\n\t		vmulpd		0x160(%%rbx),%%ymm15,%%ymm15		\n\t"\
		"vaddpd		%%ymm6		,%%ymm5,%%ymm5				\n\t		vaddpd		%%ymm14		,%%ymm13,%%ymm13		\n\t"\
		"vmovaps		%%ymm1		,%%ymm3					\n\t		vmovaps		%%ymm9		,%%ymm11				\n\t"\
		"vsubpd		%%ymm7		,%%ymm4,%%ymm4				\n\t		vsubpd		%%ymm15		,%%ymm12,%%ymm12		\n\t"\
		"addq		$0x080		,%%rcx						\n\t		vaddpd		%%ymm12		,%%ymm8	,%%ymm8			\n\t"\
		"addq		$0x0c0		,%%rbx						\n\t		vaddpd		%%ymm13		,%%ymm9	,%%ymm9			\n\t"\
		"vmovaps		     (%%rcx),%%ymm6					\n\t		vsubpd		%%ymm12		,%%ymm10,%%ymm10		\n\t"\
		"vmovaps		0x020(%%rcx),%%ymm7					\n\t		vsubpd		%%ymm13		,%%ymm11,%%ymm11		\n\t"\
		"vaddpd		%%ymm4		,%%ymm0,%%ymm0				\n\t		vmovaps		0x100(%%rcx),%%ymm12				\n\t"\
		"vaddpd		%%ymm5		,%%ymm1,%%ymm1				\n\t		vmovaps		0x120(%%rcx),%%ymm13				\n\t"\
		"vsubpd		%%ymm4		,%%ymm2,%%ymm2				\n\t		vmovaps		0x100(%%rcx),%%ymm14				\n\t"\
		"vsubpd		%%ymm5		,%%ymm3,%%ymm3				\n\t		vmovaps		0x120(%%rcx),%%ymm15				\n\t"\
		"vmovaps		%%ymm6		,%%ymm4					\n\t		vmulpd		0x100(%%rbx),%%ymm12,%%ymm12		\n\t"\
		"vmovaps		%%ymm7		,%%ymm5					\n\t		vmulpd		0x100(%%rbx),%%ymm13,%%ymm13		\n\t"\
		"vmulpd		     (%%rbx),%%ymm4,%%ymm4				\n\t		vmulpd		0x120(%%rbx),%%ymm14,%%ymm14		\n\t"\
		"vmulpd		     (%%rbx),%%ymm5,%%ymm5				\n\t		vmulpd		0x120(%%rbx),%%ymm15,%%ymm15		\n\t"\
		"vmulpd		0x020(%%rbx),%%ymm6,%%ymm6				\n\t		vaddpd		%%ymm14		,%%ymm13,%%ymm13		\n\t"\
		"vmulpd		0x020(%%rbx),%%ymm7,%%ymm7				\n\t		vsubpd		%%ymm15		,%%ymm12,%%ymm12		\n\t"\
		"/*movq		%%rax		,%%rdx		*/				\n\t		vmovaps		%%ymm13,0x120(%%rdx)				\n\t"\
		"vaddpd		%%ymm6		,%%ymm5,%%ymm5				\n\t		vmovaps		%%ymm12,0x100(%%rdx)				\n\t"\
		"vsubpd		%%ymm7		,%%ymm4,%%ymm4				\n\t"\
		"vmovaps		%%ymm5		,0x020(%%rdx)			\n\t		addq	$0x080,%%rax							\n\t"\
		"vmovaps		%%ymm4		,     (%%rdx)			\n\t		subq	$0x040,%%rbx							\n\t"\
		"vmovaps		     (%%rax),%%ymm4					\n\t		vmovaps		0x100(%%rax),%%ymm12				\n\t"\
		"vmovaps		0x020(%%rax),%%ymm5					\n\t		vmovaps		0x120(%%rax),%%ymm13				\n\t"\
		"vmovaps		%%ymm4		,%%ymm6					\n\t		vmovaps		0x100(%%rax),%%ymm14				\n\t"\
		"vmovaps		%%ymm5		,%%ymm7					\n\t		vmovaps		0x120(%%rax),%%ymm15				\n\t"\
		"vmulpd		     (%%rbx),%%ymm4,%%ymm4				\n\t		vmulpd		0x100(%%rbx),%%ymm12,%%ymm12		\n\t"\
		"vmulpd		     (%%rbx),%%ymm5,%%ymm5				\n\t		vmulpd		0x100(%%rbx),%%ymm13,%%ymm13		\n\t"\
		"vmulpd		0x020(%%rbx),%%ymm6,%%ymm6				\n\t		vmulpd		0x120(%%rbx),%%ymm14,%%ymm14		\n\t"\
		"vmulpd		0x020(%%rbx),%%ymm7,%%ymm7				\n\t		vmulpd		0x120(%%rbx),%%ymm15,%%ymm15		\n\t"\
		"vaddpd		%%ymm6		,%%ymm5,%%ymm5				\n\t		vaddpd		%%ymm14		,%%ymm13,%%ymm13		\n\t"\
		"vsubpd		%%ymm7		,%%ymm4,%%ymm4				\n\t		vsubpd		%%ymm15		,%%ymm12,%%ymm12		\n\t"\
		"vmovaps		%%ymm5		,%%ymm7					\n\t		vmovaps		%%ymm13		,%%ymm15				\n\t"\
		"vmovaps		%%ymm4		,%%ymm6					\n\t		vmovaps		%%ymm12		,%%ymm14				\n\t"\
		"vsubpd		     (%%rdx),%%ymm4,%%ymm4				\n\t		vsubpd		0x100(%%rdx),%%ymm12,%%ymm12		\n\t"\
		"vsubpd		0x020(%%rdx),%%ymm5,%%ymm5				\n\t		vsubpd		0x120(%%rdx),%%ymm13,%%ymm13		\n\t"\
		"vaddpd		     (%%rdx),%%ymm6,%%ymm6				\n\t		vaddpd		0x100(%%rdx),%%ymm14,%%ymm14		\n\t"\
		"vaddpd		0x020(%%rdx),%%ymm7,%%ymm7				\n\t		vaddpd		0x120(%%rdx),%%ymm15,%%ymm15		\n\t"\
		"vsubpd		%%ymm6		,%%ymm0,%%ymm0				\n\t		vsubpd		%%ymm14		,%%ymm8	,%%ymm8			\n\t"\
		"vsubpd		%%ymm5		,%%ymm2,%%ymm2				\n\t		vsubpd		%%ymm15		,%%ymm9	,%%ymm9			\n\t"\
		"vsubpd		%%ymm7		,%%ymm1,%%ymm1				\n\t"\
		"vsubpd		%%ymm4		,%%ymm3,%%ymm3				\n\t		vsubpd		%%ymm13		,%%ymm10,%%ymm10		\n\t"\
		"/*vmovaps		%%ymm0		,0x080(%%rdx)	*/		\n\t"\
		"/*vmovaps		%%ymm2		,0x040(%%rdx)	*/		\n\t		vsubpd		%%ymm12		,%%ymm11,%%ymm11		\n\t"\
		"/*vmovaps		%%ymm1		,0x0a0(%%rdx)	*/		\n\t		vmulpd		(%%rdi)		,%%ymm14,%%ymm14		\n\t"\
		"/*vmovaps		%%ymm3		,0x0e0(%%rdx)	*/		\n\t		vmulpd		(%%rdi)		,%%ymm13,%%ymm13		\n\t"\
		"vmulpd		(%%rdi)		,%%ymm6,%%ymm6				\n\t		vmulpd		(%%rdi)		,%%ymm15,%%ymm15		\n\t"\
		"vmulpd		(%%rdi)		,%%ymm5,%%ymm5				\n\t		vmulpd		(%%rdi)		,%%ymm12,%%ymm12		\n\t"\
		"vmulpd		(%%rdi)		,%%ymm7,%%ymm7				\n\t		vaddpd		%%ymm8		,%%ymm14,%%ymm14		\n\t"\
		"vmulpd		(%%rdi)		,%%ymm4,%%ymm4				\n\t		vaddpd		%%ymm10		,%%ymm13,%%ymm13		\n\t"\
		"vaddpd		%%ymm0		,%%ymm6,%%ymm6				\n\t		vaddpd		%%ymm9		,%%ymm15,%%ymm15		\n\t"\
		"vaddpd		%%ymm2		,%%ymm5,%%ymm5				\n\t		vaddpd		%%ymm11		,%%ymm12,%%ymm12		\n\t"\
		"vaddpd		%%ymm1		,%%ymm7,%%ymm7				\n\t		vmovaps		%%ymm14,0x100(%%rdx)				\n\t"\
		"vaddpd		%%ymm3		,%%ymm4,%%ymm4				\n\t		vmovaps		%%ymm15,0x120(%%rdx)				\n\t"\
		"/*vmovaps		%%ymm6		,     (%%rdx)	*/		\n\t		vmovaps		%%ymm10,%%ymm14						\n\t"\
		"/*vmovaps		%%ymm5		,0x0c0(%%rdx)	*/		\n\t		vmovaps		%%ymm13,%%ymm15						\n\t"\
		"/*vmovaps		%%ymm7		,0x020(%%rdx)	*/		\n\t		vsubpd		%%ymm12		,%%ymm10,%%ymm10		\n\t"\
		"/*vmovaps		%%ymm4		,0x060(%%rdx)	*/		\n\t		vsubpd		%%ymm11		,%%ymm13,%%ymm13		\n\t"\
		"													\n\t		vaddpd		%%ymm12		,%%ymm14,%%ymm14		\n\t"\
		"													\n\t		vaddpd		%%ymm11		,%%ymm15,%%ymm15		\n\t"\
		"													\n\t		vmulpd		(%%rsi)		,%%ymm10,%%ymm10		\n\t"\
		"													\n\t		vmulpd		(%%rsi)		,%%ymm13,%%ymm13		\n\t"\
		"													\n\t		vmulpd		(%%rsi)		,%%ymm14,%%ymm14		\n\t"\
		"													\n\t		vmulpd		(%%rsi)		,%%ymm15,%%ymm15		\n\t"\
		"/*****	SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r20) *****/\n\t"\
		"																vmovaps		0x100(%%rdx),%%ymm11		\n\t"\
		"																vmovaps		0x120(%%rdx),%%ymm12		\n\t"\
		"vsubpd		%%ymm11		,%%ymm6	,%%ymm6				\n\t		vsubpd		%%ymm10		,%%ymm2	,%%ymm2			\n\t"\
		"vsubpd		%%ymm9		,%%ymm0	,%%ymm0				\n\t		vsubpd		%%ymm15		,%%ymm5	,%%ymm5			\n\t"\
		"vsubpd		%%ymm12		,%%ymm7	,%%ymm7				\n\t		vsubpd		%%ymm14		,%%ymm4	,%%ymm4			\n\t"\
		"vsubpd		%%ymm8		,%%ymm1	,%%ymm1				\n\t		vsubpd		%%ymm13		,%%ymm3	,%%ymm3			\n\t"\
		"vmulpd		(%%rdi)		,%%ymm11,%%ymm11			\n\t		vmulpd		(%%rdi)		,%%ymm10,%%ymm10		\n\t"\
		"vmulpd		(%%rdi)		,%%ymm9	,%%ymm9				\n\t		vmulpd		(%%rdi)		,%%ymm15,%%ymm15		\n\t"\
		"vmulpd		(%%rdi)		,%%ymm12,%%ymm12			\n\t		vmulpd		(%%rdi)		,%%ymm14,%%ymm14		\n\t"\
		"vmulpd		(%%rdi)		,%%ymm8	,%%ymm8				\n\t		vmulpd		(%%rdi)		,%%ymm13,%%ymm13		\n\t"\
		"vaddpd		%%ymm6		,%%ymm11,%%ymm11			\n\t		vaddpd		%%ymm2		,%%ymm10,%%ymm10		\n\t"\
		"vaddpd		%%ymm0		,%%ymm9	,%%ymm9				\n\t		vaddpd		%%ymm5		,%%ymm15,%%ymm15		\n\t"\
		"vaddpd		%%ymm7		,%%ymm12,%%ymm12			\n\t		vaddpd		%%ymm4		,%%ymm14,%%ymm14		\n\t"\
		"vaddpd		%%ymm1		,%%ymm8	,%%ymm8				\n\t		vaddpd		%%ymm3		,%%ymm13,%%ymm13		\n\t"\
		"vmovaps		%%ymm6		,0x100(%%rdx)			\n\t		vmovaps		%%ymm2		,0x140(%%rdx)	\n\t"\
		"vmovaps		%%ymm0		,0x080(%%rdx)			\n\t		vmovaps		%%ymm5		,0x0c0(%%rdx)	\n\t"\
		"vmovaps		%%ymm7		,0x120(%%rdx)			\n\t		vmovaps		%%ymm4		,0x160(%%rdx)	\n\t"\
		"vmovaps		%%ymm1		,0x1a0(%%rdx)			\n\t		vmovaps		%%ymm3		,0x1e0(%%rdx)	\n\t"\
		"vmovaps		%%ymm11		,     (%%rdx)			\n\t		vmovaps		%%ymm10		,0x040(%%rdx)	\n\t"\
		"vmovaps		%%ymm9		,0x180(%%rdx)			\n\t		vmovaps		%%ymm15		,0x1c0(%%rdx)	\n\t"\
		"vmovaps		%%ymm12		,0x020(%%rdx)			\n\t		vmovaps		%%ymm14		,0x060(%%rdx)	\n\t"\
		"vmovaps		%%ymm8		,0x0a0(%%rdx)			\n\t		vmovaps		%%ymm13		,0x0e0(%%rdx)	\n\t"\
		"/*...Block 4:	*/\n\t"\
		"/*****	SSE2_RADIX4_DIF_4TWIDDLE         (r30,r34,r32,r36,r30,c03)	*****/\n\t"\
		"movq		%[__c03]	,%%rbx					\n\t"\
		"movq		%[__r30]	,%%rax					\n\t"\
		"movq		%%rax		,%%rcx		/*****	SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO(r38,r3C,r3A,r3E,r38,c07)	*****/\n\t"\
		"addq		$0x040		,%%rcx					\n\t"\
		"vmovaps		 (%%rax),%%ymm0	\n\t	movq %%rax,%%rdx \n\t	vmovaps		0x100(%%rax),%%ymm8			\n\t"\
		"vmovaps		 (%%rcx),%%ymm4						\n\t		vmovaps		0x100(%%rcx),%%ymm12		\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1						\n\t		vmovaps		0x120(%%rax),%%ymm9			\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5						\n\t		vmovaps		0x120(%%rcx),%%ymm13		\n\t"\
		"vmovaps		 (%%rbx),%%ymm6						\n\t		vmovaps		0x100(%%rbx),%%ymm14		\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm7						\n\t		vmovaps		0x120(%%rbx),%%ymm15		\n\t"\
		"vmovaps	%%ymm0		,%%ymm2						\n\t		vmovaps		%%ymm8		,%%ymm10		\n\t"\
		"vmovaps	%%ymm1		,%%ymm3						\n\t		vmovaps		%%ymm9		,%%ymm11		\n\t"\
		"vmulpd		%%ymm6		,%%ymm0,%%ymm0				\n\t		vmulpd		%%ymm14		,%%ymm8	,%%ymm8			\n\t"\
		"vmulpd		%%ymm6		,%%ymm1,%%ymm1				\n\t		vmulpd		%%ymm14		,%%ymm9	,%%ymm9			\n\t"\
		"vmulpd		%%ymm7		,%%ymm2,%%ymm2				\n\t		vmulpd		%%ymm15		,%%ymm10,%%ymm10		\n\t"\
		"vmulpd		%%ymm7		,%%ymm3,%%ymm3				\n\t		vmulpd		%%ymm15		,%%ymm11,%%ymm11		\n\t"\
		"vmovaps		%%ymm4		,%%ymm6					\n\t		vmovaps		%%ymm12		,%%ymm14				\n\t"\
		"vaddpd		%%ymm2		,%%ymm1,%%ymm1				\n\t		vaddpd		%%ymm10		,%%ymm9	,%%ymm9			\n\t"\
		"vmovaps		%%ymm5		,%%ymm7					\n\t		vmovaps		%%ymm13		,%%ymm15				\n\t"\
		"vmulpd		0x040(%%rbx),%%ymm4,%%ymm4				\n\t		vmulpd		0x140(%%rbx),%%ymm12,%%ymm12		\n\t"\
		"vsubpd		%%ymm3		,%%ymm0,%%ymm0				\n\t		vsubpd		%%ymm11		,%%ymm8	,%%ymm8			\n\t"\
		"vmulpd		0x040(%%rbx),%%ymm5,%%ymm5				\n\t		vmulpd		0x140(%%rbx),%%ymm13,%%ymm13		\n\t"\
		"vmulpd		0x060(%%rbx),%%ymm6,%%ymm6				\n\t		vmulpd		0x160(%%rbx),%%ymm14,%%ymm14		\n\t"\
		"vmovaps		%%ymm0		,%%ymm2					\n\t		vmovaps		%%ymm8		,%%ymm10				\n\t"\
		"vmulpd		0x060(%%rbx),%%ymm7,%%ymm7				\n\t		vmulpd		0x160(%%rbx),%%ymm15,%%ymm15		\n\t"\
		"vaddpd		%%ymm6		,%%ymm5,%%ymm5				\n\t		vaddpd		%%ymm14		,%%ymm13,%%ymm13		\n\t"\
		"vmovaps		%%ymm1		,%%ymm3					\n\t		vmovaps		%%ymm9		,%%ymm11				\n\t"\
		"vsubpd		%%ymm7		,%%ymm4,%%ymm4				\n\t		vsubpd		%%ymm15		,%%ymm12,%%ymm12		\n\t"\
		"addq		$0x080		,%%rcx						\n\t		vaddpd		%%ymm12		,%%ymm8	,%%ymm8			\n\t"\
		"addq		$0x0c0		,%%rbx						\n\t		vaddpd		%%ymm13		,%%ymm9	,%%ymm9			\n\t"\
		"vmovaps		     (%%rcx),%%ymm6					\n\t		vsubpd		%%ymm12		,%%ymm10,%%ymm10		\n\t"\
		"vmovaps		0x020(%%rcx),%%ymm7					\n\t		vsubpd		%%ymm13		,%%ymm11,%%ymm11		\n\t"\
		"vaddpd		%%ymm4		,%%ymm0,%%ymm0				\n\t		vmovaps		0x100(%%rcx),%%ymm12				\n\t"\
		"vaddpd		%%ymm5		,%%ymm1,%%ymm1				\n\t		vmovaps		0x120(%%rcx),%%ymm13				\n\t"\
		"vsubpd		%%ymm4		,%%ymm2,%%ymm2				\n\t		vmovaps		0x100(%%rcx),%%ymm14				\n\t"\
		"vsubpd		%%ymm5		,%%ymm3,%%ymm3				\n\t		vmovaps		0x120(%%rcx),%%ymm15				\n\t"\
		"vmovaps		%%ymm6		,%%ymm4					\n\t		vmulpd		0x100(%%rbx),%%ymm12,%%ymm12		\n\t"\
		"vmovaps		%%ymm7		,%%ymm5					\n\t		vmulpd		0x100(%%rbx),%%ymm13,%%ymm13		\n\t"\
		"vmulpd		     (%%rbx),%%ymm4,%%ymm4				\n\t		vmulpd		0x120(%%rbx),%%ymm14,%%ymm14		\n\t"\
		"vmulpd		     (%%rbx),%%ymm5,%%ymm5				\n\t		vmulpd		0x120(%%rbx),%%ymm15,%%ymm15		\n\t"\
		"vmulpd		0x020(%%rbx),%%ymm6,%%ymm6				\n\t		vaddpd		%%ymm14		,%%ymm13,%%ymm13		\n\t"\
		"vmulpd		0x020(%%rbx),%%ymm7,%%ymm7				\n\t		vsubpd		%%ymm15		,%%ymm12,%%ymm12		\n\t"\
		"/*movq		%%rax		,%%rdx,%%rdx		*/		\n\t		vmovaps		%%ymm13,0x120(%%rdx)				\n\t"\
		"vaddpd		%%ymm6		,%%ymm5,%%ymm5				\n\t		vmovaps		%%ymm12,0x100(%%rdx)				\n\t"\
		"vsubpd		%%ymm7		,%%ymm4,%%ymm4				\n\t"\
		"vmovaps		%%ymm5		,0x020(%%rdx)			\n\t		addq	$0x080,%%rax							\n\t"\
		"vmovaps		%%ymm4		,     (%%rdx)			\n\t		subq	$0x040,%%rbx							\n\t"\
		"vmovaps		     (%%rax),%%ymm4					\n\t		vmovaps		0x100(%%rax),%%ymm12				\n\t"\
		"vmovaps		0x020(%%rax),%%ymm5					\n\t		vmovaps		0x120(%%rax),%%ymm13				\n\t"\
		"vmovaps		%%ymm4		,%%ymm6					\n\t		vmovaps		0x100(%%rax),%%ymm14				\n\t"\
		"vmovaps		%%ymm5		,%%ymm7					\n\t		vmovaps		0x120(%%rax),%%ymm15				\n\t"\
		"vmulpd		     (%%rbx),%%ymm4,%%ymm4				\n\t		vmulpd		0x100(%%rbx),%%ymm12,%%ymm12		\n\t"\
		"vmulpd		     (%%rbx),%%ymm5,%%ymm5				\n\t		vmulpd		0x100(%%rbx),%%ymm13,%%ymm13		\n\t"\
		"vmulpd		0x020(%%rbx),%%ymm6,%%ymm6				\n\t		vmulpd		0x120(%%rbx),%%ymm14,%%ymm14		\n\t"\
		"vmulpd		0x020(%%rbx),%%ymm7,%%ymm7				\n\t		vmulpd		0x120(%%rbx),%%ymm15,%%ymm15		\n\t"\
		"vaddpd		%%ymm6		,%%ymm5,%%ymm5				\n\t		vaddpd		%%ymm14		,%%ymm13,%%ymm13		\n\t"\
		"vsubpd		%%ymm7		,%%ymm4,%%ymm4				\n\t		vsubpd		%%ymm15		,%%ymm12,%%ymm12		\n\t"\
		"vmovaps		%%ymm5		,%%ymm7					\n\t		vmovaps		%%ymm13		,%%ymm15				\n\t"\
		"vmovaps		%%ymm4		,%%ymm6					\n\t		vmovaps		%%ymm12		,%%ymm14				\n\t"\
		"vsubpd		     (%%rdx),%%ymm4,%%ymm4				\n\t		vsubpd		0x100(%%rdx),%%ymm12,%%ymm12		\n\t"\
		"vsubpd		0x020(%%rdx),%%ymm5,%%ymm5				\n\t		vsubpd		0x120(%%rdx),%%ymm13,%%ymm13		\n\t"\
		"vaddpd		     (%%rdx),%%ymm6,%%ymm6				\n\t		vaddpd		0x100(%%rdx),%%ymm14,%%ymm14		\n\t"\
		"vaddpd		0x020(%%rdx),%%ymm7,%%ymm7				\n\t		vaddpd		0x120(%%rdx),%%ymm15,%%ymm15		\n\t"\
		"vsubpd		%%ymm6		,%%ymm0,%%ymm0				\n\t		vsubpd		%%ymm14		,%%ymm8	,%%ymm8			\n\t"\
		"vsubpd		%%ymm5		,%%ymm2,%%ymm2				\n\t		vsubpd		%%ymm15		,%%ymm9	,%%ymm9			\n\t"\
		"vsubpd		%%ymm7		,%%ymm1,%%ymm1				\n\t"\
		"vsubpd		%%ymm4		,%%ymm3,%%ymm3				\n\t		vsubpd		%%ymm13		,%%ymm10,%%ymm10		\n\t"\
		"/*vmovaps		%%ymm0		,0x080(%%rdx)	*/		\n\t"\
		"/*vmovaps		%%ymm2		,0x040(%%rdx)	*/		\n\t		vsubpd		%%ymm12		,%%ymm11,%%ymm11		\n\t"\
		"/*vmovaps		%%ymm1		,0x0a0(%%rdx)	*/		\n\t		vmulpd		(%%rdi)		,%%ymm14,%%ymm14		\n\t"\
		"/*vmovaps		%%ymm3		,0x0e0(%%rdx)	*/		\n\t		vmulpd		(%%rdi)		,%%ymm13,%%ymm13		\n\t"\
		"vmulpd		(%%rdi)		,%%ymm6,%%ymm6				\n\t		vmulpd		(%%rdi)		,%%ymm15,%%ymm15		\n\t"\
		"vmulpd		(%%rdi)		,%%ymm5,%%ymm5				\n\t		vmulpd		(%%rdi)		,%%ymm12,%%ymm12		\n\t"\
		"vmulpd		(%%rdi)		,%%ymm7,%%ymm7				\n\t		vaddpd		%%ymm8		,%%ymm14,%%ymm14		\n\t"\
		"vmulpd		(%%rdi)		,%%ymm4,%%ymm4				\n\t		vaddpd		%%ymm10		,%%ymm13,%%ymm13		\n\t"\
		"vaddpd		%%ymm0		,%%ymm6,%%ymm6				\n\t		vaddpd		%%ymm9		,%%ymm15,%%ymm15		\n\t"\
		"vaddpd		%%ymm2		,%%ymm5,%%ymm5				\n\t		vaddpd		%%ymm11		,%%ymm12,%%ymm12		\n\t"\
		"vaddpd		%%ymm1		,%%ymm7,%%ymm7				\n\t		vmovaps		%%ymm14,0x100(%%rdx)				\n\t"\
		"vaddpd		%%ymm3		,%%ymm4,%%ymm4				\n\t		vmovaps		%%ymm15,0x120(%%rdx)				\n\t"\
		"/*vmovaps		%%ymm6		,     (%%rdx)	*/		\n\t		vmovaps		%%ymm10,%%ymm14						\n\t"\
		"/*vmovaps		%%ymm5		,0x0c0(%%rdx)	*/		\n\t		vmovaps		%%ymm13,%%ymm15						\n\t"\
		"/*vmovaps		%%ymm7		,0x020(%%rdx)	*/		\n\t		vsubpd		%%ymm12		,%%ymm10,%%ymm10		\n\t"\
		"/*vmovaps		%%ymm4		,0x060(%%rdx)	*/		\n\t		vsubpd		%%ymm11		,%%ymm13,%%ymm13		\n\t"\
		"													\n\t		vaddpd		%%ymm12		,%%ymm14,%%ymm14		\n\t"\
		"													\n\t		vaddpd		%%ymm11		,%%ymm15,%%ymm15		\n\t"\
		"													\n\t		vmulpd		(%%rsi)		,%%ymm10,%%ymm10		\n\t"\
		"													\n\t		vmulpd		(%%rsi)		,%%ymm13,%%ymm13		\n\t"\
		"													\n\t		vmulpd		(%%rsi)		,%%ymm14,%%ymm14		\n\t"\
		"													\n\t		vmulpd		(%%rsi)		,%%ymm15,%%ymm15		\n\t"\
		"/*****	SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r20) *****/\n\t"\
		"																vmovaps		0x100(%%rdx),%%ymm11		\n\t"\
		"																vmovaps		0x120(%%rdx),%%ymm12		\n\t"\
		"vsubpd		%%ymm11		,%%ymm6	,%%ymm6				\n\t		vsubpd		%%ymm10		,%%ymm2	,%%ymm2			\n\t"\
		"vsubpd		%%ymm9		,%%ymm0	,%%ymm0				\n\t		vsubpd		%%ymm15		,%%ymm5	,%%ymm5			\n\t"\
		"vsubpd		%%ymm12		,%%ymm7	,%%ymm7				\n\t		vsubpd		%%ymm14		,%%ymm4	,%%ymm4			\n\t"\
		"vsubpd		%%ymm8		,%%ymm1	,%%ymm1				\n\t		vsubpd		%%ymm13		,%%ymm3	,%%ymm3			\n\t"\
		"vmulpd		(%%rdi)		,%%ymm11,%%ymm11			\n\t		vmulpd		(%%rdi)		,%%ymm10,%%ymm10		\n\t"\
		"vmulpd		(%%rdi)		,%%ymm9	,%%ymm9				\n\t		vmulpd		(%%rdi)		,%%ymm15,%%ymm15		\n\t"\
		"vmulpd		(%%rdi)		,%%ymm12,%%ymm12			\n\t		vmulpd		(%%rdi)		,%%ymm14,%%ymm14		\n\t"\
		"vmulpd		(%%rdi)		,%%ymm8	,%%ymm8				\n\t		vmulpd		(%%rdi)		,%%ymm13,%%ymm13		\n\t"\
		"vaddpd		%%ymm6		,%%ymm11,%%ymm11			\n\t		vaddpd		%%ymm2		,%%ymm10,%%ymm10		\n\t"\
		"vaddpd		%%ymm0		,%%ymm9	,%%ymm9				\n\t		vaddpd		%%ymm5		,%%ymm15,%%ymm15		\n\t"\
		"vaddpd		%%ymm7		,%%ymm12,%%ymm12			\n\t		vaddpd		%%ymm4		,%%ymm14,%%ymm14		\n\t"\
		"vaddpd		%%ymm1		,%%ymm8	,%%ymm8				\n\t		vaddpd		%%ymm3		,%%ymm13,%%ymm13		\n\t"\
		"vmovaps		%%ymm6		,0x100(%%rdx)			\n\t		vmovaps		%%ymm2		,0x140(%%rdx)	\n\t"\
		"vmovaps		%%ymm0		,0x080(%%rdx)			\n\t		vmovaps		%%ymm5		,0x0c0(%%rdx)	\n\t"\
		"vmovaps		%%ymm7		,0x120(%%rdx)			\n\t		vmovaps		%%ymm4		,0x160(%%rdx)	\n\t"\
		"vmovaps		%%ymm1		,0x1a0(%%rdx)			\n\t		vmovaps		%%ymm3		,0x1e0(%%rdx)	\n\t"\
		"vmovaps		%%ymm11		,     (%%rdx)			\n\t		vmovaps		%%ymm10		,0x040(%%rdx)	\n\t"\
		"vmovaps		%%ymm9		,0x180(%%rdx)			\n\t		vmovaps		%%ymm15		,0x1c0(%%rdx)	\n\t"\
		"vmovaps		%%ymm12		,0x020(%%rdx)			\n\t		vmovaps		%%ymm14		,0x060(%%rdx)	\n\t"\
		"vmovaps		%%ymm8		,0x0a0(%%rdx)			\n\t		vmovaps		%%ymm13		,0x0e0(%%rdx)	\n\t"\
		"/**********************************************************************************/\n\t"\
		"/*...and now do eight radix-4 transforms, including the internal twiddle factors: */\n\t"\
		"/**********************************************************************************/\n\t"\
		"movq		%[__isrt2]		,%%rsi		\n\t"\
		"/*...Block 1: t00,t10,t20,t30	*/					\n\t		/*...Block 5: t08,t18,t28,t38	*/		\n\t"\
		"movq		%[__r00]		,%%rax					\n\t		movq		$0x100,%%r10			\n\t"\
		"movq		%[__r10]		,%%rbx					\n\t		movq		$0x100,%%r11			\n\t"\
		"movq		%[__r20]		,%%rcx					\n\t		movq		$0x100,%%r12			\n\t"\
		"movq		%[__r30]		,%%rdx					\n\t		movq		$0x100,%%r13			\n\t"\
		"vmovaps		      (%%rax)	,%%ymm0				\n\t		addq		%%rax ,%%r10			\n\t"\
		"vmovaps		 0x020(%%rax)	,%%ymm1				\n\t		addq		%%rbx ,%%r11			\n\t"\
		"vmovaps		      (%%rbx)	,%%ymm2				\n\t		addq		%%rcx ,%%r12			\n\t"\
		"vmovaps		 0x020(%%rbx)	,%%ymm3				\n\t		addq		%%rdx ,%%r13			\n\t"\
		"vaddpd		       %%ymm0	,%%ymm2,%%ymm2			\n\t		vmovaps		    (%%r10),%%ymm8			\n\t"\
		"vaddpd		       %%ymm1	,%%ymm3,%%ymm3			\n\t		vmovaps		0x20(%%r10),%%ymm9			\n\t"\
		"vsubpd		      (%%rbx)	,%%ymm0,%%ymm0			\n\t		vmovaps		    (%%r11),%%ymm10			\n\t"\
		"vsubpd		 0x020(%%rbx)	,%%ymm1,%%ymm1			\n\t		vmovaps		0x20(%%r11),%%ymm11			\n\t"\
		"vmovaps		      (%%rcx)	,%%ymm4				\n\t		vaddpd		     %%ymm9 ,%%ymm10,%%ymm10	\n\t"\
		"vmovaps		 0x020(%%rcx)	,%%ymm5				\n\t		vaddpd		     %%ymm8 ,%%ymm11,%%ymm11	\n\t"\
		"vmovaps		      (%%rdx)	,%%ymm6				\n\t		vsubpd		0x020(%%r11),%%ymm8	,%%ymm8		\n\t"\
		"vmovaps		 0x020(%%rdx)	,%%ymm7				\n\t		vsubpd		     (%%r11),%%ymm9	,%%ymm9		\n\t"\
		"vaddpd		       %%ymm4	,%%ymm6,%%ymm6			\n\t		vmovaps		    (%%r12),%%ymm12					\n\t"\
		"vaddpd		       %%ymm5	,%%ymm7,%%ymm7			\n\t		vmovaps		0x20(%%r12),%%ymm13					\n\t"\
		"vsubpd		      (%%rdx)	,%%ymm4,%%ymm4			\n\t		vmovaps		    (%%r13),%%ymm14					\n\t"\
		"vsubpd		 0x020(%%rdx)	,%%ymm5,%%ymm5			\n\t		vmovaps		0x20(%%r13),%%ymm15					\n\t"\
		"vsubpd		%%ymm6			,%%ymm2,%%ymm2			\n\t		vsubpd		     %%ymm13,%%ymm12,%%ymm12	\n\t"\
		"vsubpd		%%ymm7			,%%ymm3,%%ymm3			\n\t		vaddpd		    (%%r12)	,%%ymm13,%%ymm13	\n\t"\
		"vmulpd		(%%rdi)			,%%ymm6,%%ymm6			\n\t		vaddpd		     %%ymm15,%%ymm14,%%ymm14	\n\t"\
		"vmulpd		(%%rdi)			,%%ymm7,%%ymm7			\n\t		vsubpd		    (%%r13)	,%%ymm15,%%ymm15	\n\t"\
		"vmovaps		%%ymm2		,      (%%rcx)			\n\t		vmulpd		(%%rsi)		,%%ymm12,%%ymm12		\n\t"\
		"vmovaps		%%ymm3		, 0x020(%%rcx)			\n\t		vmulpd		(%%rsi)		,%%ymm13,%%ymm13		\n\t"\
		"vaddpd		%%ymm2			,%%ymm6,%%ymm6			\n\t		vmulpd		(%%rsi)		,%%ymm14,%%ymm14		\n\t"\
		"vaddpd		%%ymm3			,%%ymm7,%%ymm7			\n\t		vmulpd		(%%rsi)		,%%ymm15,%%ymm15		\n\t"\
		"vmovaps		%%ymm6		,      (%%rax)			\n\t		vsubpd		%%ymm14		,%%ymm12,%%ymm12		\n\t"\
		"vmovaps		%%ymm7		, 0x020(%%rax)			\n\t		vsubpd		%%ymm15		,%%ymm13,%%ymm13		\n\t"\
		"vsubpd		%%ymm5			,%%ymm0,%%ymm0			\n\t		vmulpd		(%%rdi)		,%%ymm14,%%ymm14		\n\t"\
		"vsubpd		%%ymm4			,%%ymm1,%%ymm1			\n\t		vmulpd		(%%rdi)		,%%ymm15,%%ymm15		\n\t"\
		"vmulpd		(%%rdi)			,%%ymm5,%%ymm5			\n\t		vaddpd		%%ymm12		,%%ymm14,%%ymm14		\n\t"\
		"vmulpd		(%%rdi)			,%%ymm4,%%ymm4			\n\t		vaddpd		%%ymm13		,%%ymm15,%%ymm15		\n\t"\
		"vmovaps		%%ymm0		,      (%%rbx)			\n\t		vsubpd		%%ymm12		,%%ymm8	,%%ymm8			\n\t"\
		"vmovaps		%%ymm1		, 0x020(%%rdx)			\n\t		vsubpd		%%ymm13		,%%ymm10,%%ymm10		\n\t"\
		"vaddpd		%%ymm0			,%%ymm5,%%ymm5			\n\t		vmulpd		(%%rdi)		,%%ymm12,%%ymm12		\n\t"\
		"vaddpd		%%ymm1			,%%ymm4,%%ymm4			\n\t		vmulpd		(%%rdi)		,%%ymm13,%%ymm13		\n\t"\
		"vmovaps		%%ymm5		,      (%%rdx)			\n\t		vmovaps		%%ymm8 ,     (%%r12)				\n\t"\
		"vmovaps		%%ymm4		, 0x020(%%rbx)			\n\t		vmovaps		%%ymm10,0x020(%%r12)				\n\t"\
		"/*...Block 3: t04,t14,t24,t34	*/					\n\t"\
		"addq		$0x080		,%%rax						\n\t		vaddpd		%%ymm8		,%%ymm12,%%ymm12		\n\t"\
		"addq		$0x080		,%%rbx						\n\t		vaddpd		%%ymm10		,%%ymm13,%%ymm13		\n\t"\
		"addq		$0x080		,%%rcx						\n\t		vmovaps		%%ymm12,     (%%r10)				\n\t"\
		"addq		$0x080		,%%rdx						\n\t		vmovaps		%%ymm13,0x020(%%r10)				\n\t"\
		"vmovaps		0x020(%%rsi)	,%%ymm8		/* c */	\n\t		vsubpd		%%ymm15		,%%ymm11,%%ymm11		\n\t"\
		"vmovaps		0x040(%%rsi)	,%%ymm10	/* s */	\n\t		vsubpd		%%ymm14		,%%ymm9	,%%ymm9			\n\t"\
		"vmovaps		     (%%rcx)	,%%ymm4				\n\t		vmulpd		(%%rdi)		,%%ymm15,%%ymm15		\n\t"\
		"vmovaps		     (%%rdx)	,%%ymm6				\n\t		vmulpd		(%%rdi)		,%%ymm14,%%ymm14		\n\t"\
		"vmovaps		0x020(%%rcx)	,%%ymm5				\n\t		vmovaps		%%ymm11,     (%%r11)				\n\t"\
		"vmovaps		0x020(%%rdx)	,%%ymm7				\n\t		vmovaps		%%ymm9 ,0x020(%%r13)				\n\t"\
		"vmovaps		      %%ymm4	,%%ymm0				\n\t		vaddpd		%%ymm11		,%%ymm15,%%ymm15		\n\t"\
		"vmovaps		      %%ymm6	,%%ymm2				\n\t		vaddpd		%%ymm9		,%%ymm14,%%ymm14		\n\t"\
		"vmovaps		      %%ymm5	,%%ymm1				\n\t		vmovaps		%%ymm15,     (%%r13)				\n\t"\
		"vmovaps		      %%ymm7	,%%ymm3				\n\t		vmovaps		%%ymm14,0x020(%%r11)				\n\t"\
		"													/*...Block 7: t0C,t1C,t2C,t3C	*/					\n\t"\
		"vmulpd		%%ymm8		,%%ymm4,%%ymm4				\n\t		addq		$0x080,%%r10					\n\t"\
		"vmulpd		%%ymm10		,%%ymm6,%%ymm6				\n\t		addq		$0x080,%%r11					\n\t"\
		"vmulpd		%%ymm10		,%%ymm1,%%ymm1				\n\t		addq		$0x080,%%r12					\n\t"\
		"vmulpd		%%ymm8		,%%ymm3,%%ymm3				\n\t		addq		$0x080,%%r13					\n\t"\
		"vmulpd		%%ymm8		,%%ymm5,%%ymm5				\n\t		vmovaps		     (%%r12),%%ymm12				\n\t"\
		"vmulpd		%%ymm10		,%%ymm7,%%ymm7				\n\t		vmovaps		     (%%r13),%%ymm14				\n\t"\
		"vmulpd		%%ymm10		,%%ymm0,%%ymm0				\n\t		vmovaps		0x020(%%r12),%%ymm13				\n\t"\
		"vmulpd		%%ymm8		,%%ymm2,%%ymm2				\n\t		vmovaps		0x020(%%r13),%%ymm15				\n\t"\
		"vsubpd		%%ymm1		,%%ymm4,%%ymm4				\n\t		vmovaps		     %%ymm13,%%ymm9					\n\t"\
		"vsubpd		%%ymm3		,%%ymm6,%%ymm6				\n\t		vmovaps		     %%ymm15,%%ymm11				\n\t"\
		"vaddpd		%%ymm0		,%%ymm5,%%ymm5				\n\t		vmulpd		 %%ymm10	,%%ymm12,%%ymm12			\n\t"\
		"vaddpd		%%ymm2		,%%ymm7,%%ymm7				\n\t		vmulpd		 %%ymm8		,%%ymm14,%%ymm14			\n\t"\
		"vsubpd		%%ymm6		,%%ymm4,%%ymm4				\n\t		vmulpd		 %%ymm8		,%%ymm9	,%%ymm9				\n\t"\
		"vsubpd		%%ymm7		,%%ymm5,%%ymm5				\n\t		vmulpd		 %%ymm10	,%%ymm11,%%ymm11			\n\t"\
		"vmulpd		(%%rdi)		,%%ymm6,%%ymm6				\n\t		vmulpd		 %%ymm10	,%%ymm13,%%ymm13			\n\t"\
		"vmulpd		(%%rdi)		,%%ymm7,%%ymm7				\n\t		vmulpd		 %%ymm8		,%%ymm15,%%ymm15			\n\t"\
		"vaddpd		%%ymm4		,%%ymm6,%%ymm6				\n\t		vmulpd		(%%r12)		,%%ymm8	,%%ymm8				\n\t"\
		"vaddpd		%%ymm5		,%%ymm7,%%ymm7				\n\t		vmulpd		(%%r13)		,%%ymm10,%%ymm10			\n\t"\
		"vmovaps		      (%%rbx)	,%%ymm2				\n\t		vsubpd		%%ymm9		,%%ymm12,%%ymm12		\n\t"\
		"vmovaps		 0x020(%%rbx)	,%%ymm3				\n\t		vsubpd		%%ymm11		,%%ymm14,%%ymm14		\n\t"\
		"vsubpd		0x020(%%rbx),%%ymm2,%%ymm2				\n\t		vaddpd		%%ymm8		,%%ymm13,%%ymm13		\n\t"\
		"vaddpd		     (%%rbx),%%ymm3,%%ymm3				\n\t		vaddpd		%%ymm10		,%%ymm15,%%ymm15		\n\t"\
		"vmulpd		     (%%rsi),%%ymm2,%%ymm2				\n\t		vsubpd		%%ymm14		,%%ymm12,%%ymm12		\n\t"\
		"vmulpd		     (%%rsi),%%ymm3,%%ymm3				\n\t		vsubpd		%%ymm15		,%%ymm13,%%ymm13		\n\t"\
		"vmovaps		      (%%rax)	,%%ymm0				\n\t		vmulpd		(%%rdi)		,%%ymm14,%%ymm14		\n\t"\
		"vmovaps		 0x020(%%rax)	,%%ymm1				\n\t		vmulpd		(%%rdi)		,%%ymm15,%%ymm15		\n\t"\
		"vsubpd		%%ymm2		,%%ymm0,%%ymm0				\n\t		vaddpd		%%ymm12		,%%ymm14,%%ymm14		\n\t"\
		"vsubpd		%%ymm3		,%%ymm1,%%ymm1				\n\t		vaddpd		%%ymm13		,%%ymm15,%%ymm15		\n\t"\
		"vaddpd		     (%%rax),%%ymm2,%%ymm2				\n\t		vmovaps		     (%%r11),%%ymm10				\n\t"\
		"vaddpd		0x020(%%rax),%%ymm3,%%ymm3				\n\t		vmovaps		0x020(%%r11),%%ymm11				\n\t"\
		"vsubpd		%%ymm6		,%%ymm2,%%ymm2				\n\t		vaddpd		0x020(%%r11),%%ymm10,%%ymm10	\n\t"\
		"vsubpd		%%ymm7		,%%ymm3,%%ymm3				\n\t		vsubpd		     (%%r11),%%ymm11,%%ymm11	\n\t"\
		"vmulpd		(%%rdi)		,%%ymm6,%%ymm6				\n\t		vmulpd		     (%%rsi),%%ymm10,%%ymm10	\n\t"\
		"vmulpd		(%%rdi)		,%%ymm7,%%ymm7				\n\t		vmulpd		     (%%rsi),%%ymm11,%%ymm11	\n\t"\
		"vmovaps		%%ymm2		,      (%%rcx)			\n\t		vmovaps		      (%%r10),%%ymm8				\n\t"\
		"vmovaps		%%ymm3		, 0x020(%%rcx)			\n\t		vmovaps		 0x020(%%r10),%%ymm9				\n\t"\
		"vaddpd		%%ymm2		,%%ymm6,%%ymm6				\n\t		vsubpd		%%ymm10		,%%ymm8	,%%ymm8			\n\t"\
		"vaddpd		%%ymm3		,%%ymm7,%%ymm7				\n\t		vsubpd		%%ymm11		,%%ymm9	,%%ymm9			\n\t"\
		"vmovaps		%%ymm6		,      (%%rax)			\n\t		vaddpd		     (%%r10),%%ymm10,%%ymm10	\n\t"\
		"vmovaps		%%ymm7		, 0x020(%%rax)			\n\t		vaddpd		0x020(%%r10),%%ymm11,%%ymm11	\n\t"\
		"vsubpd		%%ymm5		,%%ymm0,%%ymm0				\n\t		vsubpd		%%ymm12		,%%ymm8	,%%ymm8			\n\t"\
		"vsubpd		%%ymm4		,%%ymm1,%%ymm1				\n\t		vsubpd		%%ymm13		,%%ymm9	,%%ymm9			\n\t"\
		"vmulpd		(%%rdi)		,%%ymm5,%%ymm5				\n\t		vmulpd		(%%rdi)		,%%ymm12,%%ymm12		\n\t"\
		"vmulpd		(%%rdi)		,%%ymm4,%%ymm4				\n\t		vmulpd		(%%rdi)		,%%ymm13,%%ymm13		\n\t"\
		"vmovaps		%%ymm0		,      (%%rbx)			\n\t		vmovaps		%%ymm8 ,     (%%r12)				\n\t"\
		"vmovaps		%%ymm1		, 0x020(%%rdx)			\n\t		vmovaps		%%ymm9 ,0x020(%%r12)				\n\t"\
		"vaddpd		%%ymm0		,%%ymm5,%%ymm5				\n\t		vaddpd		%%ymm8		,%%ymm12,%%ymm12		\n\t"\
		"vaddpd		%%ymm1		,%%ymm4,%%ymm4				\n\t		vaddpd		%%ymm9		,%%ymm13,%%ymm13		\n\t"\
		"vmovaps	%%ymm5		,      (%%rdx)				\n\t		vmovaps		%%ymm12,     (%%r10)				\n\t"\
		"vmovaps	%%ymm4		, 0x020(%%rbx)				\n\t		vmovaps		%%ymm13,0x020(%%r10)				\n\t"\
		"/*...Block 2: t02,t12,t22,t32	*/\n\t"\
		"subq		$0x040		,%%rax						\n\t		vsubpd		%%ymm15		,%%ymm10,%%ymm10		\n\t"\
		"subq		$0x040		,%%rbx						\n\t		vsubpd		%%ymm14		,%%ymm11,%%ymm11		\n\t"\
		"subq		$0x040		,%%rcx						\n\t		vmulpd		(%%rdi)		,%%ymm15,%%ymm15		\n\t"\
		"subq		$0x040		,%%rdx						\n\t		vmulpd		(%%rdi)		,%%ymm14,%%ymm14		\n\t"\
		"addq		$0x060		,%%rsi	/* cc1 */			\n\t		vmovaps		%%ymm10,     (%%r11)				\n\t"\
		"vmovaps	     (%%rcx)	,%%ymm4					\n\t		vmovaps		%%ymm11,0x020(%%r13)				\n\t"\
		"vmovaps	     (%%rdx)	,%%ymm6					\n\t		vaddpd		%%ymm10		,%%ymm15,%%ymm15		\n\t"\
		"vmovaps	0x020(%%rcx)	,%%ymm5					\n\t		vaddpd		%%ymm11		,%%ymm14,%%ymm14		\n\t"\
		"vmovaps	0x020(%%rdx)	,%%ymm7					\n\t		vmovaps		%%ymm15,     (%%r13)		\n\t"\
		"vmovaps	     (%%rcx)	,%%ymm0					\n\t		vmovaps		%%ymm14,0x020(%%r11)		\n\t"\
		"																/*...Block 6: t0A,t1A,t2A,t3A	*/\n\t"\
		"vmovaps	     (%%rdx)	,%%ymm2					\n\t		subq		$0x040		,%%r10			\n\t"\
		"vmovaps	0x020(%%rcx)	,%%ymm1					\n\t		subq		$0x040		,%%r11			\n\t"\
		"vmovaps	0x020(%%rdx)	,%%ymm3					\n\t		subq		$0x040		,%%r12			\n\t"\
		"vmulpd		     (%%rsi)	,%%ymm4,%%ymm4			\n\t		subq		$0x040		,%%r13			\n\t"\
		"vmulpd		0x040(%%rsi)	,%%ymm6,%%ymm6			\n\t		vmovaps		      (%%r12)	,%%ymm12	\n\t"\
		"vmulpd		0x020(%%rsi)	,%%ymm1,%%ymm1			\n\t		vmovaps		      (%%r13)	,%%ymm14	\n\t"\
		"vmulpd		0x060(%%rsi)	,%%ymm3,%%ymm3			\n\t		vmovaps		 0x020(%%r12)	,%%ymm13	\n\t"\
		"vmulpd		     (%%rsi)	,%%ymm5,%%ymm5			\n\t		vmovaps		 0x020(%%r13)	,%%ymm15	\n\t"\
		"vmulpd		0x040(%%rsi)	,%%ymm7,%%ymm7			\n\t		vmovaps		      (%%r12)	,%%ymm8		\n\t"\
		"vmulpd		0x020(%%rsi)	,%%ymm0,%%ymm0			\n\t		vmovaps		      (%%r13)	,%%ymm10	\n\t"\
		"vmulpd		0x060(%%rsi)	,%%ymm2,%%ymm2			\n\t		vmovaps		 0x020(%%r12)	,%%ymm9		\n\t"\
		"vsubpd		%%ymm1			,%%ymm4,%%ymm4			\n\t		vmovaps		 0x020(%%r13)	,%%ymm11	\n\t"\
		"vsubpd		%%ymm3			,%%ymm6,%%ymm6			\n\t		vmulpd		0x060(%%rsi),%%ymm12,%%ymm12	\n\t"\
		"vaddpd		%%ymm0			,%%ymm5,%%ymm5			\n\t		vmulpd		     (%%rsi),%%ymm14,%%ymm14	\n\t"\
		"vaddpd		%%ymm2			,%%ymm7,%%ymm7			\n\t		vmulpd		0x040(%%rsi),%%ymm9	,%%ymm9		\n\t"\
		"vsubpd		%%ymm6			,%%ymm4,%%ymm4			\n\t		vmulpd		0x020(%%rsi),%%ymm11,%%ymm11	\n\t"\
		"vsubpd		%%ymm7			,%%ymm5,%%ymm5			\n\t		vmulpd		0x060(%%rsi),%%ymm13,%%ymm13	\n\t"\
		"vmulpd		(%%rdi)			,%%ymm6,%%ymm6			\n\t		vmulpd		     (%%rsi),%%ymm15,%%ymm15	\n\t"\
		"vmulpd		(%%rdi)			,%%ymm7,%%ymm7			\n\t		vmulpd		0x040(%%rsi),%%ymm8	,%%ymm8		\n\t"\
		"vaddpd		%%ymm4			,%%ymm6,%%ymm6			\n\t		vmulpd		0x020(%%rsi),%%ymm10,%%ymm10	\n\t"\
		"vaddpd		%%ymm5			,%%ymm7,%%ymm7			\n\t		vsubpd		%%ymm9		,%%ymm12,%%ymm12		\n\t"\
		"vmovaps		      (%%rbx)	,%%ymm2				\n\t		vaddpd		%%ymm11		,%%ymm14,%%ymm14		\n\t"\
		"vmovaps		 0x020(%%rbx)	,%%ymm0				\n\t		vaddpd		%%ymm8		,%%ymm13,%%ymm13		\n\t"\
		"vmovaps		      (%%rbx)	,%%ymm1				\n\t		vsubpd		%%ymm10		,%%ymm15,%%ymm15		\n\t"\
		"vmovaps		 0x020(%%rbx)	,%%ymm3				\n\t		vsubpd		%%ymm14		,%%ymm12,%%ymm12		\n\t"\
		"vmulpd		-0x040(%%rsi)	,%%ymm2,%%ymm2			\n\t		vsubpd		%%ymm15		,%%ymm13,%%ymm13		\n\t"\
		"vmulpd		-0x020(%%rsi)	,%%ymm0,%%ymm0			\n\t		vmulpd		(%%rdi)		,%%ymm14,%%ymm14		\n\t"\
		"vmulpd		-0x040(%%rsi)	,%%ymm3,%%ymm3			\n\t		vmulpd		(%%rdi)		,%%ymm15,%%ymm15		\n\t"\
		"vmulpd		-0x020(%%rsi)	,%%ymm1,%%ymm1			\n\t		vaddpd		%%ymm12		,%%ymm14,%%ymm14		\n\t"\
		"vsubpd		%%ymm0			,%%ymm2,%%ymm2			\n\t		vaddpd		%%ymm13		,%%ymm15,%%ymm15		\n\t"\
		"vaddpd		%%ymm1			,%%ymm3,%%ymm3			\n\t		vmovaps		     (%%r11),%%ymm10				\n\t"\
		"vmovaps		      (%%rax)	,%%ymm0				\n\t		vmovaps		0x020(%%r11),%%ymm8					\n\t"\
		"vmovaps		 0x020(%%rax)	,%%ymm1				\n\t		vmovaps		     (%%r11),%%ymm9					\n\t"\
		"vsubpd		%%ymm2			,%%ymm0,%%ymm0			\n\t		vmovaps		0x020(%%r11),%%ymm11				\n\t"\
		"vsubpd		%%ymm3			,%%ymm1,%%ymm1			\n\t		vmulpd		-0x20(%%rsi),%%ymm10,%%ymm10	\n\t"\
		"vaddpd		      (%%rax)	,%%ymm2,%%ymm2			\n\t		vmulpd		-0x40(%%rsi),%%ymm8	,%%ymm8		\n\t"\
		"vaddpd		 0x020(%%rax)	,%%ymm3,%%ymm3			\n\t		vmulpd		-0x20(%%rsi),%%ymm11,%%ymm11	\n\t"\
		"vsubpd		%%ymm6			,%%ymm2,%%ymm2			\n\t		vmulpd		-0x40(%%rsi),%%ymm9	,%%ymm9		\n\t"\
		"vsubpd		%%ymm7			,%%ymm3,%%ymm3			\n\t		vaddpd		%%ymm8		,%%ymm10,%%ymm10		\n\t"\
		"vmulpd		(%%rdi)			,%%ymm6,%%ymm6			\n\t		vsubpd		%%ymm9		,%%ymm11,%%ymm11		\n\t"\
		"vmulpd		(%%rdi)			,%%ymm7,%%ymm7			\n\t		vmovaps		     (%%r10),%%ymm8					\n\t"\
		"vmovaps		%%ymm2		,      (%%rcx)			\n\t		vmovaps		0x020(%%r10),%%ymm9					\n\t"\
		"vmovaps		%%ymm3		, 0x020(%%rcx)			\n\t		vsubpd		%%ymm10		,%%ymm8	,%%ymm8			\n\t"\
		"vaddpd		%%ymm2			,%%ymm6,%%ymm6			\n\t		vsubpd		%%ymm11		,%%ymm9	,%%ymm9			\n\t"\
		"vaddpd		%%ymm3			,%%ymm7,%%ymm7			\n\t		vaddpd		     (%%r10),%%ymm10,%%ymm10	\n\t"\
		"vmovaps		%%ymm6		,      (%%rax)			\n\t		vaddpd		0x020(%%r10),%%ymm11,%%ymm11	\n\t"\
		"vmovaps		%%ymm7		, 0x020(%%rax)			\n\t		vsubpd		%%ymm12		,%%ymm8	,%%ymm8			\n\t"\
		"vsubpd		%%ymm5			,%%ymm0,%%ymm0			\n\t		vsubpd		%%ymm13		,%%ymm9	,%%ymm9			\n\t"\
		"vsubpd		%%ymm4			,%%ymm1,%%ymm1			\n\t		vmulpd		(%%rdi)		,%%ymm12,%%ymm12		\n\t"\
		"vmulpd		(%%rdi)			,%%ymm5,%%ymm5			\n\t		vmulpd		(%%rdi)		,%%ymm13,%%ymm13		\n\t"\
		"vmulpd		(%%rdi)			,%%ymm4,%%ymm4			\n\t		vmovaps		%%ymm8 ,     (%%r12)				\n\t"\
		"vmovaps		%%ymm0		,      (%%rbx)			\n\t		vmovaps		%%ymm9 ,0x020(%%r12)				\n\t"\
		"vmovaps		%%ymm1		, 0x020(%%rdx)			\n\t		vaddpd		%%ymm8		,%%ymm12,%%ymm12		\n\t"\
		"vaddpd		%%ymm0			,%%ymm5,%%ymm5			\n\t		vaddpd		%%ymm9		,%%ymm13,%%ymm13		\n\t"\
		"vaddpd		%%ymm1			,%%ymm4,%%ymm4			\n\t		vmovaps		%%ymm12,     (%%r10)				\n\t"\
		"vmovaps		%%ymm5		,      (%%rdx)			\n\t		vmovaps		%%ymm13,0x020(%%r10)				\n\t"\
		"vmovaps		%%ymm4		, 0x020(%%rbx)			\n\t		vsubpd		%%ymm15		,%%ymm10,%%ymm10		\n\t"\
		"/*...Block 4: t06,t16,t26,t36	*/\n\t"\
		"addq		$0x080		,%%rax						\n\t		vsubpd		%%ymm14		,%%ymm11,%%ymm11		\n\t"\
		"addq		$0x080		,%%rbx						\n\t		vaddpd		%%ymm15		,%%ymm15,%%ymm15		\n\t"\
		"addq		$0x080		,%%rcx						\n\t		vaddpd		%%ymm14		,%%ymm14,%%ymm14		\n\t"\
		"addq		$0x080		,%%rdx						\n\t		vmovaps		%%ymm10,     (%%r11)				\n\t"\
		"vmovaps		      (%%rcx)	,%%ymm4				\n\t		vmovaps		%%ymm11,0x020(%%r13)				\n\t"\
		"vmovaps		      (%%rdx)	,%%ymm6				\n\t		vaddpd		%%ymm10		,%%ymm15,%%ymm15		\n\t"\
		"vmovaps		 0x020(%%rcx)	,%%ymm5				\n\t		vaddpd		%%ymm11		,%%ymm14,%%ymm14		\n\t"\
		"vmovaps		 0x020(%%rdx)	,%%ymm7				\n\t		vmovaps		%%ymm15,     (%%r13)		\n\t"\
		"vmovaps		      (%%rcx)	,%%ymm0				\n\t		vmovaps		%%ymm14,0x020(%%r11)		\n\t"\
		"																/*...Block 8: t0E,t1E,t2E,t3E	*/		\n\t"\
		"vmovaps		      (%%rdx)	,%%ymm2				\n\t		addq		$0x080		,%%r10			\n\t"\
		"vmovaps		 0x020(%%rcx)	,%%ymm1				\n\t		addq		$0x080		,%%r11			\n\t"\
		"vmovaps		 0x020(%%rdx)	,%%ymm3				\n\t		addq		$0x080		,%%r12			\n\t"\
		"vmulpd		 0x040(%%rsi)	,%%ymm4,%%ymm4			\n\t		addq		$0x080		,%%r13			\n\t"\
		"vmulpd		 0x020(%%rsi)	,%%ymm6,%%ymm6			\n\t		vmovaps		      (%%r12)	,%%ymm12	\n\t"\
		"vmulpd		 0x060(%%rsi)	,%%ymm1,%%ymm1			\n\t		vmovaps		      (%%r13)	,%%ymm14	\n\t"\
		"vmulpd		      (%%rsi)	,%%ymm3,%%ymm3			\n\t		vmovaps		 0x020(%%r12)	,%%ymm13	\n\t"\
		"vmulpd		 0x040(%%rsi)	,%%ymm5,%%ymm5			\n\t		vmovaps		 0x020(%%r13)	,%%ymm15	\n\t"\
		"vmulpd		 0x020(%%rsi)	,%%ymm7,%%ymm7			\n\t		vmovaps		      (%%r12)	,%%ymm8		\n\t"\
		"vmulpd		 0x060(%%rsi)	,%%ymm0,%%ymm0			\n\t		vmovaps		      (%%r13)	,%%ymm10	\n\t"\
		"vmulpd		      (%%rsi)	,%%ymm2,%%ymm2			\n\t		vmovaps		 0x020(%%r12)	,%%ymm9		\n\t"\
		"vsubpd		%%ymm1			,%%ymm4,%%ymm4			\n\t		vmovaps		 0x020(%%r13)	,%%ymm11	\n\t"\
		"vaddpd		%%ymm3			,%%ymm6,%%ymm6			\n\t		vmulpd		0x020(%%rsi),%%ymm12,%%ymm12	\n\t"\
		"vaddpd		%%ymm0			,%%ymm5,%%ymm5			\n\t		vmulpd		0x060(%%rsi),%%ymm14,%%ymm14	\n\t"\
		"vsubpd		%%ymm2			,%%ymm7,%%ymm7			\n\t		vmulpd		     (%%rsi),%%ymm9	,%%ymm9		\n\t"\
		"vsubpd		%%ymm6			,%%ymm4,%%ymm4			\n\t		vmulpd		0x040(%%rsi),%%ymm11,%%ymm11	\n\t"\
		"vsubpd		%%ymm7			,%%ymm5,%%ymm5			\n\t		vmulpd		0x020(%%rsi),%%ymm13,%%ymm13	\n\t"\
		"vmulpd		(%%rdi)			,%%ymm6,%%ymm6			\n\t		vmulpd		0x060(%%rsi),%%ymm15,%%ymm15	\n\t"\
		"vmulpd		(%%rdi)			,%%ymm7,%%ymm7			\n\t		vmulpd		     (%%rsi),%%ymm8	,%%ymm8		\n\t"\
		"vaddpd		%%ymm4			,%%ymm6,%%ymm6			\n\t		vmulpd		0x040(%%rsi),%%ymm10,%%ymm10	\n\t"\
		"vaddpd		%%ymm5			,%%ymm7,%%ymm7			\n\t		vsubpd		%%ymm9		,%%ymm12,%%ymm12		\n\t"\
		"vmovaps		      (%%rbx)	,%%ymm2				\n\t		vsubpd		%%ymm11		,%%ymm14,%%ymm14		\n\t"\
		"vmovaps		 0x020(%%rbx)	,%%ymm0				\n\t		vaddpd		%%ymm8		,%%ymm13,%%ymm13		\n\t"\
		"vmovaps		      (%%rbx)	,%%ymm1				\n\t		vaddpd		%%ymm10		,%%ymm15,%%ymm15		\n\t"\
		"vmovaps		 0x020(%%rbx)	,%%ymm3				\n\t		vsubpd		%%ymm14		,%%ymm12,%%ymm12		\n\t"\
		"vmulpd		-0x020(%%rsi)	,%%ymm2,%%ymm2			\n\t		vsubpd		%%ymm15		,%%ymm13,%%ymm13		\n\t"\
		"vmulpd		-0x040(%%rsi)	,%%ymm0,%%ymm0			\n\t		vmulpd		(%%rdi)		,%%ymm14,%%ymm14		\n\t"\
		"vmulpd		-0x020(%%rsi)	,%%ymm3,%%ymm3			\n\t		vmulpd		(%%rdi)		,%%ymm15,%%ymm15		\n\t"\
		"vmulpd		-0x040(%%rsi)	,%%ymm1,%%ymm1			\n\t		vaddpd		%%ymm12		,%%ymm14,%%ymm14		\n\t"\
		"vsubpd		%%ymm0			,%%ymm2,%%ymm2			\n\t		vaddpd		%%ymm13		,%%ymm15,%%ymm15		\n\t"\
		"vaddpd		%%ymm1			,%%ymm3,%%ymm3			\n\t		vmovaps		     (%%r11),%%ymm10				\n\t"\
		"vmovaps		      (%%rax)	,%%ymm0				\n\t		vmovaps		0x020(%%r11),%%ymm8					\n\t"\
		"vmovaps		 0x020(%%rax)	,%%ymm1				\n\t		vmovaps		     (%%r11),%%ymm9					\n\t"\
		"vsubpd		%%ymm2			,%%ymm0,%%ymm0			\n\t		vmovaps		0x020(%%r11),%%ymm11				\n\t"\
		"vsubpd		%%ymm3			,%%ymm1,%%ymm1			\n\t		vmulpd		-0x40(%%rsi),%%ymm10,%%ymm10	\n\t"\
		"vaddpd		      (%%rax)	,%%ymm2,%%ymm2			\n\t		vmulpd		-0x20(%%rsi),%%ymm8	,%%ymm8		\n\t"\
		"vaddpd		 0x020(%%rax)	,%%ymm3,%%ymm3			\n\t		vmulpd		-0x40(%%rsi),%%ymm11,%%ymm11	\n\t"\
		"vsubpd		%%ymm4			,%%ymm2,%%ymm2			\n\t		vmulpd		-0x20(%%rsi),%%ymm9	,%%ymm9		\n\t"\
		"vsubpd		%%ymm5			,%%ymm3,%%ymm3			\n\t		vaddpd		%%ymm8		,%%ymm10,%%ymm10		\n\t"\
		"vmulpd		(%%rdi)			,%%ymm4,%%ymm4			\n\t		vsubpd		%%ymm9		,%%ymm11,%%ymm11		\n\t"\
		"vmulpd		(%%rdi)			,%%ymm5,%%ymm5			\n\t		vmovaps		     (%%r10),%%ymm8					\n\t"\
		"vmovaps		%%ymm2		,      (%%rcx)			\n\t		vmovaps		0x020(%%r10),%%ymm9					\n\t"\
		"vmovaps		%%ymm3		, 0x020(%%rcx)			\n\t		vsubpd		%%ymm10		,%%ymm8	,%%ymm8			\n\t"\
		"vaddpd		%%ymm2			,%%ymm4,%%ymm4			\n\t		vsubpd		%%ymm11		,%%ymm9	,%%ymm9			\n\t"\
		"vaddpd		%%ymm3			,%%ymm5,%%ymm5			\n\t		vaddpd		     (%%r10),%%ymm10,%%ymm10	\n\t"\
		"vmovaps		%%ymm4		,      (%%rax)			\n\t		vaddpd		0x020(%%r10),%%ymm11,%%ymm11	\n\t"\
		"vmovaps		%%ymm5		, 0x020(%%rax)			\n\t		vsubpd		%%ymm12		,%%ymm8	,%%ymm8			\n\t"\
		"vsubpd		%%ymm7			,%%ymm0,%%ymm0			\n\t		vsubpd		%%ymm13		,%%ymm9	,%%ymm9			\n\t"\
		"vsubpd		%%ymm6			,%%ymm1,%%ymm1			\n\t		vmulpd		(%%rdi)		,%%ymm12,%%ymm12		\n\t"\
		"vmulpd		(%%rdi)			,%%ymm7,%%ymm7			\n\t		vmulpd		(%%rdi)		,%%ymm13,%%ymm13		\n\t"\
		"vmulpd		(%%rdi)			,%%ymm6,%%ymm6			\n\t		vmovaps		%%ymm8 ,     (%%r12)				\n\t"\
		"vmovaps		%%ymm0		,      (%%rbx)			\n\t		vmovaps		%%ymm9 ,0x020(%%r12)				\n\t"\
		"vmovaps		%%ymm1		, 0x020(%%rdx)			\n\t		vaddpd		%%ymm8		,%%ymm12,%%ymm12		\n\t"\
		"vaddpd		%%ymm0			,%%ymm7,%%ymm7			\n\t		vaddpd		%%ymm9		,%%ymm13,%%ymm13		\n\t"\
		"vaddpd		%%ymm1			,%%ymm6,%%ymm6			\n\t		vmovaps		%%ymm12,     (%%r10)				\n\t"\
		"vmovaps		%%ymm7		,      (%%rdx)			\n\t		vmovaps		%%ymm13,0x020(%%r10)				\n\t"\
		"vmovaps		%%ymm6		, 0x020(%%rbx)			\n\t		vsubpd		%%ymm15		,%%ymm10,%%ymm10		\n\t"\
		"																vsubpd		%%ymm14		,%%ymm11,%%ymm11		\n\t"\
		"																vmulpd		(%%rdi)		,%%ymm15,%%ymm15		\n\t"\
		"																vmulpd		(%%rdi)		,%%ymm14,%%ymm14		\n\t"\
		"																vmovaps		%%ymm10,     (%%r11)				\n\t"\
		"																vmovaps		%%ymm11,0x020(%%r13)				\n\t"\
		"																vaddpd		%%ymm10		,%%ymm15,%%ymm15		\n\t"\
		"																vaddpd		%%ymm11		,%%ymm14,%%ymm14		\n\t"\
		"																vmovaps		%%ymm15,     (%%r13)		\n\t"\
		"																vmovaps		%%ymm14,0x020(%%r11)		\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__add1] "m" (Xadd1)\
		 ,[__add2] "m" (Xadd2)\
		 ,[__add3] "m" (Xadd3)\
		 ,[__r00] "m" (Xr00)\
		 ,[__r10] "m" (Xr10)\
		 ,[__r20] "m" (Xr20)\
		 ,[__r30] "m" (Xr30)\
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__cc0] "m" (Xcc0)\
		 ,[__c00] "m" (Xc00)\
		 ,[__c01] "m" (Xc01)\
		 ,[__c02] "m" (Xc02)\
		 ,[__c03] "m" (Xc03)\
		 ,[__c05] "m" (Xc05)\
		 ,[__c07] "m" (Xc07)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r10","r11","r12","r13","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"		/* Clobbered registers */\
	);\
	}

	#define FOO(Xadd0,Xadd1,Xr00,Xr10,Xr20,Xr30,Xisrt2,Xcc0,Xc00,Xc01,Xc02,Xc03,Xc05,Xc07)\
	{\
	__asm__ volatile (\
		"movq	%[__r00] ,%%rcx	\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__add1] "m" (Xadd1)\
		 ,[__r00] "m" (Xr00)\
		 ,[__r10] "m" (Xr10)\
		 ,[__r20] "m" (Xr20)\
		 ,[__r30] "m" (Xr30)\
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__cc0] "m" (Xcc0)\
		 ,[__c00] "m" (Xc00)\
		 ,[__c01] "m" (Xc01)\
		 ,[__c02] "m" (Xc02)\
		 ,[__c03] "m" (Xc03)\
		 ,[__c05] "m" (Xc05)\
		 ,[__c07] "m" (Xc07)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r10","r11","r12","r13","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"		/* Clobbered registers */\
	);\
	}

	#define SSE2_RADIX32_WRAPPER_DIT(Xadd0,Xadd1,Xadd2,Xadd3,Xisrt2,Xr00,Xr08,Xr10,Xr20,Xr28,Xr30,Xc01,Xc02,Xc04,Xc06,Xc08,Xc0A,Xc0C,Xc0E,Xc10,Xc12,Xc14,Xc16,Xc18,Xc1A,Xc1C,Xc1E)\
	{\
	__asm__ volatile (\
		"/************************************************************************/\n\t"\
		"/*...And do an inverse DIT radix-32 pass on the squared-data blocks.	*/\n\t"\
		"/************************************************************************/\n\t"\
		"/*...Block 1: */								\n\t"\
		"movq		%[__isrt2]		,%%rsi				\n\t"\
		"movq		%[__r00]		,%%rax				\n\t"\
		"movq		%%rax		,%%rbx					\n\t"\
		"movq		%%rax		,%%rcx					\n\t"\
		"movq		%%rax		,%%rdx					\n\t"\
		"addq		$0x400		,%%rbx					\n\t"\
		"addq		$0x200		,%%rcx					\n\t"\
		"addq		$0x600		,%%rdx					\n\t"\
		"/*****	SSE2_RADIX4_DIT_IN_PLACE()	*****/		\n\t		/*...Block 2 has tmp-addresses offset +0x80 w.r.to Block 1:	*/\n\t"\
		"vmovaps	      (%%rax)	,%%ymm0				\n\t		vmovaps		 0x080(%%rax)	,%%ymm8 			\n\t"\
		"vmovaps	 0x020(%%rax)	,%%ymm1				\n\t		vmovaps		 0x0a0(%%rax)	,%%ymm9 			\n\t"\
		"vmovaps	      (%%rax)	,%%ymm2				\n\t		vmovaps		 0x080(%%rax)	,%%ymm10			\n\t"\
		"vmovaps	 0x020(%%rax)	,%%ymm3				\n\t		vmovaps		 0x0a0(%%rax)	,%%ymm11			\n\t"\
		"vaddpd		      (%%rbx)	,%%ymm0,%%ymm0		\n\t		vaddpd		 0x080(%%rbx)	,%%ymm8 ,%%ymm8 	\n\t"\
		"vaddpd		 0x020(%%rbx)	,%%ymm1,%%ymm1		\n\t		vaddpd		 0x0a0(%%rbx)	,%%ymm9 ,%%ymm9 	\n\t"\
		"vsubpd		      (%%rbx)	,%%ymm2,%%ymm2		\n\t		vsubpd		 0x080(%%rbx)	,%%ymm10,%%ymm10	\n\t"\
		"vsubpd		 0x020(%%rbx)	,%%ymm3,%%ymm3		\n\t		vsubpd		 0x0a0(%%rbx)	,%%ymm11,%%ymm11	\n\t"\
		"vmovaps	      (%%rcx)	,%%ymm4				\n\t		vmovaps		 0x080(%%rcx)	,%%ymm12			\n\t"\
		"vmovaps	 0x020(%%rcx)	,%%ymm5				\n\t		vmovaps		 0x0a0(%%rcx)	,%%ymm13			\n\t"\
		"vmovaps	      (%%rcx)	,%%ymm6				\n\t		vmovaps		 0x080(%%rcx)	,%%ymm14			\n\t"\
		"vmovaps	 0x020(%%rcx)	,%%ymm7				\n\t		vmovaps		 0x0a0(%%rcx)	,%%ymm15			\n\t"\
		"vaddpd		      (%%rdx)	,%%ymm4,%%ymm4		\n\t		vaddpd		 0x080(%%rdx)	,%%ymm12,%%ymm12	\n\t"\
		"vaddpd		 0x020(%%rdx)	,%%ymm5,%%ymm5		\n\t		vaddpd		 0x0a0(%%rdx)	,%%ymm13,%%ymm13	\n\t"\
		"vsubpd		      (%%rdx)	,%%ymm6,%%ymm6		\n\t		vsubpd		 0x080(%%rdx)	,%%ymm14,%%ymm14	\n\t"\
		"vsubpd		 0x020(%%rdx)	,%%ymm7,%%ymm7		\n\t		vsubpd		 0x0a0(%%rdx)	,%%ymm15,%%ymm15	\n\t"\
		"vsubpd			%%ymm4		,%%ymm0,%%ymm0		\n\t		vsubpd			%%ymm12		,%%ymm8 ,%%ymm8 	\n\t"\
		"vsubpd			%%ymm5		,%%ymm1,%%ymm1		\n\t		vsubpd			%%ymm13		,%%ymm9 ,%%ymm9 	\n\t"\
		"vmovaps		%%ymm0		,      (%%rbx)		\n\t		vmovaps			%%ymm8 		, 0x080(%%rbx)			\n\t"\
		"vmovaps		%%ymm1		, 0x020(%%rbx)		\n\t		vmovaps			%%ymm9 		, 0x0a0(%%rbx)			\n\t"\
		"vaddpd			%%ymm4		,%%ymm4,%%ymm4		\n\t		vaddpd			%%ymm12		,%%ymm12,%%ymm12	\n\t"\
		"vaddpd			%%ymm5		,%%ymm5,%%ymm5		\n\t		vaddpd			%%ymm13		,%%ymm13,%%ymm13	\n\t"\
		"vaddpd			%%ymm0		,%%ymm4,%%ymm4		\n\t		vaddpd			%%ymm8 		,%%ymm12,%%ymm12	\n\t"\
		"vaddpd			%%ymm1		,%%ymm5,%%ymm5		\n\t		vaddpd			%%ymm9 		,%%ymm13,%%ymm13	\n\t"\
		"vmovaps		%%ymm4		,      (%%rax)		\n\t		vmovaps			%%ymm12		, 0x080(%%rax)			\n\t"\
		"vmovaps		%%ymm5		, 0x020(%%rax)		\n\t		vmovaps			%%ymm13		, 0x0a0(%%rax)			\n\t"\
		"vsubpd			%%ymm7		,%%ymm2,%%ymm2		\n\t		vsubpd			%%ymm15		,%%ymm10,%%ymm10	\n\t"\
		"vsubpd			%%ymm6		,%%ymm3,%%ymm3		\n\t		vsubpd			%%ymm14		,%%ymm11,%%ymm11	\n\t"\
		"vmovaps		%%ymm2		,      (%%rdx)		\n\t		vmovaps			%%ymm10		, 0x080(%%rdx)			\n\t"\
		"vmovaps		%%ymm3		, 0x020(%%rcx)		\n\t		vmovaps			%%ymm11		, 0x0a0(%%rcx)			\n\t"\
		"vaddpd			%%ymm7		,%%ymm7,%%ymm7		\n\t		vaddpd			%%ymm15		,%%ymm15,%%ymm15	\n\t"\
		"vaddpd			%%ymm6		,%%ymm6,%%ymm6		\n\t		vaddpd			%%ymm14		,%%ymm14,%%ymm14	\n\t"\
		"vaddpd			%%ymm2		,%%ymm7,%%ymm7		\n\t		vaddpd			%%ymm10		,%%ymm15,%%ymm15	\n\t"\
		"vaddpd			%%ymm3		,%%ymm6,%%ymm6		\n\t		vaddpd			%%ymm11		,%%ymm14,%%ymm14	\n\t"\
		"vmovaps		%%ymm7		,      (%%rcx)		\n\t		vmovaps			%%ymm15		, 0x080(%%rcx)			\n\t"\
		"vmovaps		%%ymm6		, 0x020(%%rdx)		\n\t		vmovaps			%%ymm14		, 0x0a0(%%rdx)			\n\t"\
		"addq		$0x100		,%%rax					\n\t"\
		"addq		$0x100		,%%rbx					\n\t"\
		"addq		$0x100		,%%rcx					\n\t"\
		"addq		$0x100		,%%rdx					\n\t"\
		"vmovaps	      (%%rax)	,%%ymm0				\n\t		vmovaps		 0x080(%%rax)	,%%ymm8 			\n\t"\
		"vmovaps	 0x020(%%rax)	,%%ymm1				\n\t		vmovaps		 0x0a0(%%rax)	,%%ymm9 			\n\t"\
		"vmovaps	      (%%rax)	,%%ymm2				\n\t		vmovaps		 0x080(%%rax)	,%%ymm10			\n\t"\
		"vmovaps	 0x020(%%rax)	,%%ymm3				\n\t		vmovaps		 0x0a0(%%rax)	,%%ymm11			\n\t"\
		"vaddpd		      (%%rbx)	,%%ymm0,%%ymm0		\n\t		vaddpd		 0x080(%%rbx)	,%%ymm8 ,%%ymm8 	\n\t"\
		"vaddpd		 0x020(%%rbx)	,%%ymm1,%%ymm1		\n\t		vaddpd		 0x0a0(%%rbx)	,%%ymm9 ,%%ymm9 	\n\t"\
		"vsubpd		      (%%rbx)	,%%ymm2,%%ymm2		\n\t		vsubpd		 0x080(%%rbx)	,%%ymm10,%%ymm10	\n\t"\
		"vsubpd		 0x020(%%rbx)	,%%ymm3,%%ymm3		\n\t		vsubpd		 0x0a0(%%rbx)	,%%ymm11,%%ymm11	\n\t"\
		"vmovaps	      (%%rcx)	,%%ymm4				\n\t		vmovaps		 0x080(%%rcx)	,%%ymm12			\n\t"\
		"vmovaps	 0x020(%%rcx)	,%%ymm5				\n\t		vmovaps		 0x0a0(%%rcx)	,%%ymm13			\n\t"\
		"vmovaps	      (%%rcx)	,%%ymm6				\n\t		vmovaps		 0x080(%%rcx)	,%%ymm14			\n\t"\
		"vmovaps	 0x020(%%rcx)	,%%ymm7				\n\t		vmovaps		 0x0a0(%%rcx)	,%%ymm15			\n\t"\
		"vaddpd		      (%%rdx)	,%%ymm4,%%ymm4		\n\t		vaddpd		 0x080(%%rdx)	,%%ymm12,%%ymm12	\n\t"\
		"vaddpd		 0x020(%%rdx)	,%%ymm5,%%ymm5		\n\t		vaddpd		 0x0a0(%%rdx)	,%%ymm13,%%ymm13	\n\t"\
		"vsubpd		      (%%rdx)	,%%ymm6,%%ymm6		\n\t		vsubpd		 0x080(%%rdx)	,%%ymm14,%%ymm14	\n\t"\
		"vsubpd		 0x020(%%rdx)	,%%ymm7,%%ymm7		\n\t		vsubpd		 0x0a0(%%rdx)	,%%ymm15,%%ymm15	\n\t"\
		"vsubpd			%%ymm4		,%%ymm0,%%ymm0		\n\t		vsubpd			%%ymm12		,%%ymm8 ,%%ymm8 	\n\t"\
		"vsubpd			%%ymm5		,%%ymm1,%%ymm1		\n\t		vsubpd			%%ymm13		,%%ymm9 ,%%ymm9 	\n\t"\
		"vmovaps		%%ymm0		,      (%%rbx)		\n\t		vmovaps			%%ymm8 		, 0x080(%%rbx)			\n\t"\
		"vmovaps		%%ymm1		, 0x020(%%rbx)		\n\t		vmovaps			%%ymm9 		, 0x0a0(%%rbx)			\n\t"\
		"vaddpd			%%ymm4		,%%ymm4,%%ymm4		\n\t		vaddpd			%%ymm12		,%%ymm12,%%ymm12	\n\t"\
		"vaddpd			%%ymm5		,%%ymm5,%%ymm5		\n\t		vaddpd			%%ymm13		,%%ymm13,%%ymm13	\n\t"\
		"vaddpd			%%ymm0		,%%ymm4,%%ymm4		\n\t		vaddpd			%%ymm8 		,%%ymm12,%%ymm12	\n\t"\
		"vaddpd			%%ymm1		,%%ymm5,%%ymm5		\n\t		vaddpd			%%ymm9 		,%%ymm13,%%ymm13	\n\t"\
		"vmovaps		%%ymm4		,      (%%rax)		\n\t		vmovaps			%%ymm12		, 0x080(%%rax)			\n\t"\
		"vmovaps		%%ymm5		, 0x020(%%rax)		\n\t		vmovaps			%%ymm13		, 0x0a0(%%rax)			\n\t"\
		"vsubpd			%%ymm7		,%%ymm2,%%ymm2		\n\t		vsubpd			%%ymm15		,%%ymm10,%%ymm10	\n\t"\
		"vsubpd			%%ymm6		,%%ymm3,%%ymm3		\n\t		vsubpd			%%ymm14		,%%ymm11,%%ymm11	\n\t"\
		"vaddpd			%%ymm7		,%%ymm7,%%ymm7		\n\t		vaddpd			%%ymm15		,%%ymm15,%%ymm15	\n\t"\
		"vaddpd			%%ymm6		,%%ymm6,%%ymm6		\n\t		vaddpd			%%ymm14		,%%ymm14,%%ymm14	\n\t"\
		"vaddpd			%%ymm2		,%%ymm7,%%ymm7		\n\t		vaddpd			%%ymm10		,%%ymm15,%%ymm15	\n\t"\
		"vaddpd			%%ymm3		,%%ymm6,%%ymm6		\n\t		vaddpd			%%ymm11		,%%ymm14,%%ymm14	\n\t"\
		"vmovaps		%%ymm3		,%%ymm0				\n\t		vmovaps			%%ymm11		,%%ymm8 				\n\t"\
		"vmovaps		%%ymm6		,%%ymm1				\n\t		vmovaps			%%ymm14		,%%ymm9 				\n\t"\
		"vsubpd			%%ymm7		,%%ymm3,%%ymm3		\n\t		vsubpd			%%ymm15		,%%ymm11,%%ymm11	\n\t"\
		"vsubpd			%%ymm2		,%%ymm6,%%ymm6		\n\t		vsubpd			%%ymm10		,%%ymm14,%%ymm14	\n\t"\
		"vaddpd			%%ymm7		,%%ymm0,%%ymm0		\n\t		vaddpd			%%ymm15		,%%ymm8 ,%%ymm8 	\n\t"\
		"vaddpd			%%ymm2		,%%ymm1,%%ymm1		\n\t		vaddpd			%%ymm10		,%%ymm9 ,%%ymm9 	\n\t"\
		"vmulpd		      (%%rsi)	,%%ymm3,%%ymm3		\n\t		vmulpd		      (%%rsi)	,%%ymm11,%%ymm11	\n\t"\
		"vmulpd		      (%%rsi)	,%%ymm6,%%ymm6		\n\t		vmulpd		      (%%rsi)	,%%ymm14,%%ymm14	\n\t"\
		"vmulpd		      (%%rsi)	,%%ymm0,%%ymm0		\n\t		vmulpd		      (%%rsi)	,%%ymm8 ,%%ymm8 	\n\t"\
		"vmulpd		      (%%rsi)	,%%ymm1,%%ymm1		\n\t		vmulpd		      (%%rsi)	,%%ymm9 ,%%ymm9 	\n\t"\
		"vmovaps		%%ymm3		, 0x020(%%rcx)		\n\t		vmovaps			%%ymm11		, 0x0a0(%%rcx)			\n\t"\
		"vmovaps		%%ymm6		, 0x020(%%rdx)		\n\t		vmovaps			%%ymm14		, 0x0a0(%%rdx)			\n\t"\
		"vmovaps		%%ymm0		,      (%%rcx)		\n\t		vmovaps			%%ymm8 		, 0x080(%%rcx)			\n\t"\
		"vmovaps		%%ymm1		,      (%%rdx)		\n\t		vmovaps			%%ymm9 		, 0x080(%%rdx)			\n\t"\
		"/*****	SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS *****/\n\t		/*****	SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS *****/\n\t"\
		"/*        (r00,r10,r20,r30,r08,r18,r28,r38): */\n\t		/*        (r04,r14,r24,r34,r0C,r1C,r2C,r3C):  */\n\t"\
		"vmovaps		-0x100(%%rax)	,%%ymm0			\n\t		vmovaps		-0x080(%%rax)	,%%ymm8 			\n\t"\
		"vmovaps		-0x100(%%rbx)	,%%ymm4			\n\t		vmovaps		-0x080(%%rbx)	,%%ymm12			\n\t"\
		"vmovaps		-0x0e0(%%rax)	,%%ymm1			\n\t		vmovaps		-0x060(%%rax)	,%%ymm9 			\n\t"\
		"vmovaps		-0x0e0(%%rbx)	,%%ymm5			\n\t		vmovaps		-0x060(%%rbx)	,%%ymm13			\n\t"\
		"vmovaps		      (%%rax)	,%%ymm2			\n\t		vmovaps		 0x080(%%rax)	,%%ymm10			\n\t"\
		"vmovaps		 0x020(%%rbx)	,%%ymm7			\n\t		vmovaps		 0x0a0(%%rbx)	,%%ymm15			\n\t"\
		"vmovaps		 0x020(%%rax)	,%%ymm3			\n\t		vmovaps		 0x0a0(%%rax)	,%%ymm11			\n\t"\
		"vmovaps		      (%%rbx)	,%%ymm6			\n\t		vmovaps		 0x080(%%rbx)	,%%ymm14			\n\t"\
		"vsubpd			%%ymm2		,%%ymm0,%%ymm0		\n\t		vsubpd		%%ymm10		,%%ymm8 ,%%ymm8 		\n\t"\
		"vsubpd			%%ymm7		,%%ymm4,%%ymm4		\n\t		vsubpd		%%ymm15		,%%ymm12,%%ymm12		\n\t"\
		"vsubpd			%%ymm3		,%%ymm1,%%ymm1		\n\t		vsubpd		%%ymm11		,%%ymm9 ,%%ymm9 		\n\t"\
		"vsubpd			%%ymm6		,%%ymm5,%%ymm5		\n\t		vsubpd		%%ymm14		,%%ymm13,%%ymm13		\n\t"\
		"vaddpd			%%ymm2		,%%ymm2,%%ymm2		\n\t		vaddpd		%%ymm10		,%%ymm10,%%ymm10		\n\t"\
		"vaddpd			%%ymm7		,%%ymm7,%%ymm7		\n\t		vaddpd		%%ymm15		,%%ymm15,%%ymm15		\n\t"\
		"vaddpd			%%ymm3		,%%ymm3,%%ymm3		\n\t		vaddpd		%%ymm11		,%%ymm11,%%ymm11		\n\t"\
		"vaddpd			%%ymm6		,%%ymm6,%%ymm6		\n\t		vaddpd		%%ymm14		,%%ymm14,%%ymm14		\n\t"\
		"vaddpd			%%ymm0		,%%ymm2,%%ymm2		\n\t		vaddpd		%%ymm8 		,%%ymm10,%%ymm10		\n\t"\
		"vaddpd			%%ymm4		,%%ymm7,%%ymm7		\n\t		vaddpd		%%ymm12		,%%ymm15,%%ymm15		\n\t"\
		"vaddpd			%%ymm1		,%%ymm3,%%ymm3		\n\t		vaddpd		%%ymm9 		,%%ymm11,%%ymm11		\n\t"\
		"vaddpd			%%ymm5		,%%ymm6,%%ymm6		\n\t		vaddpd		%%ymm13		,%%ymm14,%%ymm14		\n\t"\
		"vmovaps		%%ymm0		,      (%%rax)		\n\t		vmovaps		%%ymm8 		, 0x080(%%rax)			\n\t"\
		"vmovaps		%%ymm4		,      (%%rbx)		\n\t		vmovaps		%%ymm12		, 0x080(%%rbx)			\n\t"\
		"vmovaps		%%ymm1		, 0x020(%%rax)		\n\t		vmovaps		%%ymm9 		, 0x0a0(%%rax)			\n\t"\
		"vmovaps		%%ymm5		,-0x0e0(%%rbx)		\n\t		vmovaps		%%ymm13		,-0x060(%%rbx)			\n\t"\
		"vmovaps		%%ymm2		,-0x100(%%rax)		\n\t		vmovaps		%%ymm10		,-0x080(%%rax)			\n\t"\
		"vmovaps		%%ymm7		,-0x100(%%rbx)		\n\t		vmovaps		%%ymm15		,-0x080(%%rbx)			\n\t"\
		"vmovaps		%%ymm3		,-0x0e0(%%rax)		\n\t		vmovaps		%%ymm11		,-0x060(%%rax)			\n\t"\
		"vmovaps		%%ymm6		, 0x020(%%rbx)		\n\t		vmovaps		%%ymm14		, 0x0a0(%%rbx)			\n\t"\
		"vmovaps		-0x100(%%rcx)	,%%ymm0			\n\t		vmovaps		-0x080(%%rcx)	,%%ymm8 			\n\t"\
		"vmovaps		-0x100(%%rdx)	,%%ymm4			\n\t		vmovaps		-0x080(%%rdx)	,%%ymm12			\n\t"\
		"vmovaps		-0x0e0(%%rcx)	,%%ymm1			\n\t		vmovaps		-0x060(%%rcx)	,%%ymm9 			\n\t"\
		"vmovaps		-0x0e0(%%rdx)	,%%ymm5			\n\t		vmovaps		-0x060(%%rdx)	,%%ymm13			\n\t"\
		"vmovaps		      (%%rcx)	,%%ymm2			\n\t		vmovaps		 0x080(%%rcx)	,%%ymm10			\n\t"\
		"vmovaps		 0x020(%%rdx)	,%%ymm7			\n\t		vmovaps		 0x0a0(%%rdx)	,%%ymm15			\n\t"\
		"vmovaps		 0x020(%%rcx)	,%%ymm3			\n\t		vmovaps		 0x0a0(%%rcx)	,%%ymm11			\n\t"\
		"vmovaps		      (%%rdx)	,%%ymm6			\n\t		vmovaps		 0x080(%%rdx)	,%%ymm14			\n\t"\
		"vsubpd			%%ymm2		,%%ymm0,%%ymm0		\n\t		vsubpd		%%ymm10		,%%ymm8 ,%%ymm8 		\n\t"\
		"vsubpd			%%ymm7		,%%ymm4,%%ymm4		\n\t		vsubpd		%%ymm15		,%%ymm12,%%ymm12		\n\t"\
		"vsubpd			%%ymm3		,%%ymm1,%%ymm1		\n\t		vsubpd		%%ymm11		,%%ymm9 ,%%ymm9 		\n\t"\
		"vsubpd			%%ymm6		,%%ymm5,%%ymm5		\n\t		vsubpd		%%ymm14		,%%ymm13,%%ymm13		\n\t"\
		"vaddpd			%%ymm2		,%%ymm2,%%ymm2		\n\t		vaddpd		%%ymm10		,%%ymm10,%%ymm10		\n\t"\
		"vaddpd			%%ymm7		,%%ymm7,%%ymm7		\n\t		vaddpd		%%ymm15		,%%ymm15,%%ymm15		\n\t"\
		"vaddpd			%%ymm3		,%%ymm3,%%ymm3		\n\t		vaddpd		%%ymm11		,%%ymm11,%%ymm11		\n\t"\
		"vaddpd			%%ymm6		,%%ymm6,%%ymm6		\n\t		vaddpd		%%ymm14		,%%ymm14,%%ymm14		\n\t"\
		"vaddpd			%%ymm0		,%%ymm2,%%ymm2		\n\t		vaddpd		%%ymm8 		,%%ymm10,%%ymm10		\n\t"\
		"vaddpd			%%ymm4		,%%ymm7,%%ymm7		\n\t		vaddpd		%%ymm12		,%%ymm15,%%ymm15		\n\t"\
		"vaddpd			%%ymm1		,%%ymm3,%%ymm3		\n\t		vaddpd		%%ymm9 		,%%ymm11,%%ymm11		\n\t"\
		"vaddpd			%%ymm5		,%%ymm6,%%ymm6		\n\t		vaddpd		%%ymm13		,%%ymm14,%%ymm14		\n\t"\
		"vmovaps		%%ymm0		,      (%%rcx)		\n\t		vmovaps		%%ymm8 		, 0x080(%%rcx)			\n\t"\
		"vmovaps		%%ymm4		,      (%%rdx)		\n\t		vmovaps		%%ymm12		, 0x080(%%rdx)			\n\t"\
		"vmovaps		%%ymm1		, 0x020(%%rcx)		\n\t		vmovaps		%%ymm9 		, 0x0a0(%%rcx)			\n\t"\
		"vmovaps		%%ymm5		,-0x0e0(%%rdx)		\n\t		vmovaps		%%ymm13		,-0x060(%%rdx)			\n\t"\
		"vmovaps		%%ymm2		,-0x100(%%rcx)		\n\t		vmovaps		%%ymm10		,-0x080(%%rcx)			\n\t"\
		"vmovaps		%%ymm7		,-0x100(%%rdx)		\n\t		vmovaps		%%ymm15		,-0x080(%%rdx)			\n\t"\
		"vmovaps		%%ymm3		,-0x0e0(%%rcx)		\n\t		vmovaps		%%ymm11		,-0x060(%%rcx)			\n\t"\
		"vmovaps		%%ymm6		, 0x020(%%rdx)		\n\t		vmovaps		%%ymm14		, 0x0a0(%%rdx)			\n\t"\
		"/*...Blocks 3,4 have tmp-addresses offset +0x40 w.r.to Blocks 1,2, respectively (thus +0x100-0x0c0 = +0x040: */\n\t"\
		"subq		$0xc0		,%%rax					\n\t"\
		"subq		$0xc0		,%%rbx					\n\t"\
		"subq		$0xc0		,%%rcx					\n\t"\
		"subq		$0xc0		,%%rdx					\n\t"\
		"vmovaps	      (%%rax)	,%%ymm0				\n\t		vmovaps		 0x080(%%rax)	,%%ymm8 			\n\t"\
		"vmovaps	 0x020(%%rax)	,%%ymm1				\n\t		vmovaps		 0x0a0(%%rax)	,%%ymm9 			\n\t"\
		"vmovaps	      (%%rax)	,%%ymm2				\n\t		vmovaps		 0x080(%%rax)	,%%ymm10			\n\t"\
		"vmovaps	 0x020(%%rax)	,%%ymm3				\n\t		vmovaps		 0x0a0(%%rax)	,%%ymm11			\n\t"\
		"vaddpd		      (%%rbx)	,%%ymm0,%%ymm0		\n\t		vaddpd		 0x080(%%rbx)	,%%ymm8 ,%%ymm8 	\n\t"\
		"vaddpd		 0x020(%%rbx)	,%%ymm1,%%ymm1		\n\t		vaddpd		 0x0a0(%%rbx)	,%%ymm9 ,%%ymm9 	\n\t"\
		"vsubpd		      (%%rbx)	,%%ymm2,%%ymm2		\n\t		vsubpd		 0x080(%%rbx)	,%%ymm10,%%ymm10	\n\t"\
		"vsubpd		 0x020(%%rbx)	,%%ymm3,%%ymm3		\n\t		vsubpd		 0x0a0(%%rbx)	,%%ymm11,%%ymm11	\n\t"\
		"vmovaps	      (%%rcx)	,%%ymm4				\n\t		vmovaps		 0x080(%%rcx)	,%%ymm12			\n\t"\
		"vmovaps	 0x020(%%rcx)	,%%ymm5				\n\t		vmovaps		 0x0a0(%%rcx)	,%%ymm13			\n\t"\
		"vmovaps	      (%%rcx)	,%%ymm6				\n\t		vmovaps		 0x080(%%rcx)	,%%ymm14			\n\t"\
		"vmovaps	 0x020(%%rcx)	,%%ymm7				\n\t		vmovaps		 0x0a0(%%rcx)	,%%ymm15			\n\t"\
		"vaddpd		      (%%rdx)	,%%ymm4,%%ymm4		\n\t		vaddpd		 0x080(%%rdx)	,%%ymm12,%%ymm12	\n\t"\
		"vaddpd		 0x020(%%rdx)	,%%ymm5,%%ymm5		\n\t		vaddpd		 0x0a0(%%rdx)	,%%ymm13,%%ymm13	\n\t"\
		"vsubpd		      (%%rdx)	,%%ymm6,%%ymm6		\n\t		vsubpd		 0x080(%%rdx)	,%%ymm14,%%ymm14	\n\t"\
		"vsubpd		 0x020(%%rdx)	,%%ymm7,%%ymm7		\n\t		vsubpd		 0x0a0(%%rdx)	,%%ymm15,%%ymm15	\n\t"\
		"vsubpd			%%ymm4		,%%ymm0,%%ymm0		\n\t		vsubpd			%%ymm12		,%%ymm8 ,%%ymm8 	\n\t"\
		"vsubpd			%%ymm5		,%%ymm1,%%ymm1		\n\t		vsubpd			%%ymm13		,%%ymm9 ,%%ymm9 	\n\t"\
		"vmovaps		%%ymm0		,      (%%rbx)		\n\t		vmovaps			%%ymm8 		, 0x080(%%rbx)			\n\t"\
		"vmovaps		%%ymm1		, 0x020(%%rbx)		\n\t		vmovaps			%%ymm9 		, 0x0a0(%%rbx)			\n\t"\
		"vaddpd			%%ymm4		,%%ymm4,%%ymm4		\n\t		vaddpd			%%ymm12		,%%ymm12,%%ymm12	\n\t"\
		"vaddpd			%%ymm5		,%%ymm5,%%ymm5		\n\t		vaddpd			%%ymm13		,%%ymm13,%%ymm13	\n\t"\
		"vaddpd			%%ymm0		,%%ymm4,%%ymm4		\n\t		vaddpd			%%ymm8 		,%%ymm12,%%ymm12	\n\t"\
		"vaddpd			%%ymm1		,%%ymm5,%%ymm5		\n\t		vaddpd			%%ymm9 		,%%ymm13,%%ymm13	\n\t"\
		"vmovaps		%%ymm4		,      (%%rax)		\n\t		vmovaps			%%ymm12		, 0x080(%%rax)			\n\t"\
		"vmovaps		%%ymm5		, 0x020(%%rax)		\n\t		vmovaps			%%ymm13		, 0x0a0(%%rax)			\n\t"\
		"vsubpd			%%ymm7		,%%ymm2,%%ymm2		\n\t		vsubpd			%%ymm15		,%%ymm10,%%ymm10	\n\t"\
		"vsubpd			%%ymm6		,%%ymm3,%%ymm3		\n\t		vsubpd			%%ymm14		,%%ymm11,%%ymm11	\n\t"\
		"vmovaps		%%ymm2		,      (%%rdx)		\n\t		vmovaps			%%ymm10		, 0x080(%%rdx)			\n\t"\
		"vmovaps		%%ymm3		, 0x020(%%rcx)		\n\t		vmovaps			%%ymm11		, 0x0a0(%%rcx)			\n\t"\
		"vaddpd			%%ymm7		,%%ymm7,%%ymm7		\n\t		vaddpd			%%ymm15		,%%ymm15,%%ymm15	\n\t"\
		"vaddpd			%%ymm6		,%%ymm6,%%ymm6		\n\t		vaddpd			%%ymm14		,%%ymm14,%%ymm14	\n\t"\
		"vaddpd			%%ymm2		,%%ymm7,%%ymm7		\n\t		vaddpd			%%ymm10		,%%ymm15,%%ymm15	\n\t"\
		"vaddpd			%%ymm3		,%%ymm6,%%ymm6		\n\t		vaddpd			%%ymm11		,%%ymm14,%%ymm14	\n\t"\
		"vmovaps		%%ymm7		,      (%%rcx)		\n\t		vmovaps			%%ymm15		, 0x080(%%rcx)			\n\t"\
		"vmovaps		%%ymm6		, 0x020(%%rdx)		\n\t		vmovaps			%%ymm14		, 0x0a0(%%rdx)			\n\t"\
		"addq		$0x100		,%%rax					\n\t"\
		"addq		$0x100		,%%rbx					\n\t"\
		"addq		$0x100		,%%rcx					\n\t"\
		"addq		$0x100		,%%rdx					\n\t"\
		"vmovaps		  (%%rax)	,%%ymm0				\n\t		vmovaps		 0x080(%%rax)	,%%ymm8 			\n\t"\
		"vmovaps	 0x020(%%rax)	,%%ymm1				\n\t		vmovaps		 0x0a0(%%rax)	,%%ymm9 			\n\t"\
		"vmovaps		  (%%rax)	,%%ymm2				\n\t		vmovaps		 0x080(%%rax)	,%%ymm10			\n\t"\
		"vmovaps	 0x020(%%rax)	,%%ymm3				\n\t		vmovaps		 0x0a0(%%rax)	,%%ymm11			\n\t"\
		"vaddpd		      (%%rbx)	,%%ymm0,%%ymm0		\n\t		vaddpd		 0x080(%%rbx)	,%%ymm8 ,%%ymm8 	\n\t"\
		"vaddpd		 0x020(%%rbx)	,%%ymm1,%%ymm1		\n\t		vaddpd		 0x0a0(%%rbx)	,%%ymm9 ,%%ymm9 	\n\t"\
		"vsubpd		      (%%rbx)	,%%ymm2,%%ymm2		\n\t		vsubpd		 0x080(%%rbx)	,%%ymm10,%%ymm10	\n\t"\
		"vsubpd		 0x020(%%rbx)	,%%ymm3,%%ymm3		\n\t		vsubpd		 0x0a0(%%rbx)	,%%ymm11,%%ymm11	\n\t"\
		"vmovaps	      (%%rcx)	,%%ymm4				\n\t		vmovaps		 0x080(%%rcx)	,%%ymm12			\n\t"\
		"vmovaps	 0x020(%%rcx)	,%%ymm5				\n\t		vmovaps		 0x0a0(%%rcx)	,%%ymm13			\n\t"\
		"vmovaps	      (%%rcx)	,%%ymm6				\n\t		vmovaps		 0x080(%%rcx)	,%%ymm14			\n\t"\
		"vmovaps	 0x020(%%rcx)	,%%ymm7				\n\t		vmovaps		 0x0a0(%%rcx)	,%%ymm15			\n\t"\
		"vaddpd		      (%%rdx)	,%%ymm4,%%ymm4		\n\t		vaddpd		 0x080(%%rdx)	,%%ymm12,%%ymm12	\n\t"\
		"vaddpd		 0x020(%%rdx)	,%%ymm5,%%ymm5		\n\t		vaddpd		 0x0a0(%%rdx)	,%%ymm13,%%ymm13	\n\t"\
		"vsubpd		      (%%rdx)	,%%ymm6,%%ymm6		\n\t		vsubpd		 0x080(%%rdx)	,%%ymm14,%%ymm14	\n\t"\
		"vsubpd		 0x020(%%rdx)	,%%ymm7,%%ymm7		\n\t		vsubpd		 0x0a0(%%rdx)	,%%ymm15,%%ymm15	\n\t"\
		"vsubpd			%%ymm4		,%%ymm0,%%ymm0		\n\t		vsubpd			%%ymm12		,%%ymm8 ,%%ymm8 	\n\t"\
		"vsubpd			%%ymm5		,%%ymm1,%%ymm1		\n\t		vsubpd			%%ymm13		,%%ymm9 ,%%ymm9 	\n\t"\
		"vmovaps		%%ymm0		,      (%%rbx)		\n\t		vmovaps			%%ymm8 		, 0x080(%%rbx)			\n\t"\
		"vmovaps		%%ymm1		, 0x020(%%rbx)		\n\t		vmovaps			%%ymm9 		, 0x0a0(%%rbx)			\n\t"\
		"vaddpd			%%ymm4		,%%ymm4,%%ymm4		\n\t		vaddpd			%%ymm12		,%%ymm12,%%ymm12	\n\t"\
		"vaddpd			%%ymm5		,%%ymm5,%%ymm5		\n\t		vaddpd			%%ymm13		,%%ymm13,%%ymm13	\n\t"\
		"vaddpd			%%ymm0		,%%ymm4,%%ymm4		\n\t		vaddpd			%%ymm8 		,%%ymm12,%%ymm12	\n\t"\
		"vaddpd			%%ymm1		,%%ymm5,%%ymm5		\n\t		vaddpd			%%ymm9 		,%%ymm13,%%ymm13	\n\t"\
		"vmovaps		%%ymm4		,      (%%rax)		\n\t		vmovaps			%%ymm12		, 0x080(%%rax)			\n\t"\
		"vmovaps		%%ymm5		, 0x020(%%rax)		\n\t		vmovaps			%%ymm13		, 0x0a0(%%rax)			\n\t"\
		"vsubpd			%%ymm7		,%%ymm2,%%ymm2		\n\t		vsubpd			%%ymm15		,%%ymm10,%%ymm10	\n\t"\
		"vsubpd			%%ymm6		,%%ymm3,%%ymm3		\n\t		vsubpd			%%ymm14		,%%ymm11,%%ymm11	\n\t"\
		"vaddpd			%%ymm7		,%%ymm7,%%ymm7		\n\t		vaddpd			%%ymm15		,%%ymm15,%%ymm15	\n\t"\
		"vaddpd			%%ymm6		,%%ymm6,%%ymm6		\n\t		vaddpd			%%ymm14		,%%ymm14,%%ymm14	\n\t"\
		"vaddpd			%%ymm2		,%%ymm7,%%ymm7		\n\t		vaddpd			%%ymm10		,%%ymm15,%%ymm15	\n\t"\
		"vaddpd			%%ymm3		,%%ymm6,%%ymm6		\n\t		vaddpd			%%ymm11		,%%ymm14,%%ymm14	\n\t"\
		"vmovaps		%%ymm3		,%%ymm0				\n\t		vmovaps			%%ymm11		,%%ymm8 				\n\t"\
		"vmovaps		%%ymm6		,%%ymm1				\n\t		vmovaps			%%ymm14		,%%ymm9 				\n\t"\
		"vsubpd			%%ymm7		,%%ymm3,%%ymm3		\n\t		vsubpd			%%ymm15		,%%ymm11,%%ymm11	\n\t"\
		"vsubpd			%%ymm2		,%%ymm6,%%ymm6		\n\t		vsubpd			%%ymm10		,%%ymm14,%%ymm14	\n\t"\
		"vaddpd			%%ymm7		,%%ymm0,%%ymm0		\n\t		vaddpd			%%ymm15		,%%ymm8 ,%%ymm8 	\n\t"\
		"vaddpd			%%ymm2		,%%ymm1,%%ymm1		\n\t		vaddpd			%%ymm10		,%%ymm9 ,%%ymm9 	\n\t"\
		"vmulpd		      (%%rsi)	,%%ymm3,%%ymm3		\n\t		vmulpd			  (%%rsi)	,%%ymm11,%%ymm11	\n\t"\
		"vmulpd		      (%%rsi)	,%%ymm6,%%ymm6		\n\t		vmulpd			  (%%rsi)	,%%ymm14,%%ymm14	\n\t"\
		"vmulpd		      (%%rsi)	,%%ymm0,%%ymm0		\n\t		vmulpd			  (%%rsi)	,%%ymm8 ,%%ymm8 	\n\t"\
		"vmulpd		      (%%rsi)	,%%ymm1,%%ymm1		\n\t		vmulpd			  (%%rsi)	,%%ymm9 ,%%ymm9 	\n\t"\
		"vmovaps		%%ymm3		, 0x020(%%rcx)		\n\t		vmovaps			%%ymm11		, 0x0a0(%%rcx)			\n\t"\
		"vmovaps		%%ymm6		, 0x020(%%rdx)		\n\t		vmovaps			%%ymm14		, 0x0a0(%%rdx)			\n\t"\
		"vmovaps		%%ymm0		,      (%%rcx)		\n\t		vmovaps			%%ymm8 		, 0x080(%%rcx)			\n\t"\
		"vmovaps		%%ymm1		,      (%%rdx)		\n\t		vmovaps			%%ymm9 		, 0x080(%%rdx)			\n\t"\
		"/*****	SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS *****/\n\t		/*****	SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS *****/\n\t"\
		"/*        (r02,r12,r22,r32,r0A,r1A,r2A,r3A): */\n\t		/*        (r06,r16,r26,r36,r0E,r1E,r2E,r3E):  */\n\t"\
		"vmovaps		-0x100(%%rax)	,%%ymm0			\n\t		vmovaps		-0x080(%%rax)	,%%ymm8 			\n\t"\
		"vmovaps		-0x100(%%rbx)	,%%ymm4			\n\t		vmovaps		-0x080(%%rbx)	,%%ymm12			\n\t"\
		"vmovaps		-0x0e0(%%rax)	,%%ymm1			\n\t		vmovaps		-0x060(%%rax)	,%%ymm9 			\n\t"\
		"vmovaps		-0x0e0(%%rbx)	,%%ymm5			\n\t		vmovaps		-0x060(%%rbx)	,%%ymm13			\n\t"\
		"vmovaps		      (%%rax)	,%%ymm2			\n\t		vmovaps		 0x080(%%rax)	,%%ymm10			\n\t"\
		"vmovaps		 0x020(%%rbx)	,%%ymm7			\n\t		vmovaps		 0x0a0(%%rbx)	,%%ymm15			\n\t"\
		"vmovaps		 0x020(%%rax)	,%%ymm3			\n\t		vmovaps		 0x0a0(%%rax)	,%%ymm11			\n\t"\
		"vmovaps		      (%%rbx)	,%%ymm6			\n\t		vmovaps		 0x080(%%rbx)	,%%ymm14			\n\t"\
		"vsubpd			%%ymm2		,%%ymm0,%%ymm0		\n\t		vsubpd		%%ymm10		,%%ymm8 ,%%ymm8 		\n\t"\
		"vsubpd			%%ymm7		,%%ymm4,%%ymm4		\n\t		vsubpd		%%ymm15		,%%ymm12,%%ymm12		\n\t"\
		"vsubpd			%%ymm3		,%%ymm1,%%ymm1		\n\t		vsubpd		%%ymm11		,%%ymm9 ,%%ymm9 		\n\t"\
		"vsubpd			%%ymm6		,%%ymm5,%%ymm5		\n\t		vsubpd		%%ymm14		,%%ymm13,%%ymm13		\n\t"\
		"vaddpd			%%ymm2		,%%ymm2,%%ymm2		\n\t		vaddpd		%%ymm10		,%%ymm10,%%ymm10		\n\t"\
		"vaddpd			%%ymm7		,%%ymm7,%%ymm7		\n\t		vaddpd		%%ymm15		,%%ymm15,%%ymm15		\n\t"\
		"vaddpd			%%ymm3		,%%ymm3,%%ymm3		\n\t		vaddpd		%%ymm11		,%%ymm11,%%ymm11		\n\t"\
		"vaddpd			%%ymm6		,%%ymm6,%%ymm6		\n\t		vaddpd		%%ymm14		,%%ymm14,%%ymm14		\n\t"\
		"vaddpd			%%ymm0		,%%ymm2,%%ymm2		\n\t		vaddpd		%%ymm8 		,%%ymm10,%%ymm10		\n\t"\
		"vaddpd			%%ymm4		,%%ymm7,%%ymm7		\n\t		vaddpd		%%ymm12		,%%ymm15,%%ymm15		\n\t"\
		"vaddpd			%%ymm1		,%%ymm3,%%ymm3		\n\t		vaddpd		%%ymm9 		,%%ymm11,%%ymm11		\n\t"\
		"vaddpd			%%ymm5		,%%ymm6,%%ymm6		\n\t		vaddpd		%%ymm13		,%%ymm14,%%ymm14		\n\t"\
		"vmovaps		%%ymm0		,      (%%rax)		\n\t		vmovaps		%%ymm8 		, 0x080(%%rax)			\n\t"\
		"vmovaps		%%ymm4		,      (%%rbx)		\n\t		vmovaps		%%ymm12		, 0x080(%%rbx)			\n\t"\
		"vmovaps		%%ymm1		, 0x020(%%rax)		\n\t		vmovaps		%%ymm9 		, 0x0a0(%%rax)			\n\t"\
		"vmovaps		%%ymm5		,-0x0e0(%%rbx)		\n\t		vmovaps		%%ymm13		,-0x060(%%rbx)			\n\t"\
		"vmovaps		%%ymm2		,-0x100(%%rax)		\n\t		vmovaps		%%ymm10		,-0x080(%%rax)			\n\t"\
		"vmovaps		%%ymm7		,-0x100(%%rbx)		\n\t		vmovaps		%%ymm15		,-0x080(%%rbx)			\n\t"\
		"vmovaps		%%ymm3		,-0x0e0(%%rax)		\n\t		vmovaps		%%ymm11		,-0x060(%%rax)			\n\t"\
		"vmovaps		%%ymm6		, 0x020(%%rbx)		\n\t		vmovaps		%%ymm14		, 0x0a0(%%rbx)			\n\t"\
		"vmovaps		-0x100(%%rcx)	,%%ymm0			\n\t		vmovaps		-0x080(%%rcx)	,%%ymm8 			\n\t"\
		"vmovaps		-0x100(%%rdx)	,%%ymm4			\n\t		vmovaps		-0x080(%%rdx)	,%%ymm12			\n\t"\
		"vmovaps		-0x0e0(%%rcx)	,%%ymm1			\n\t		vmovaps		-0x060(%%rcx)	,%%ymm9 			\n\t"\
		"vmovaps		-0x0e0(%%rdx)	,%%ymm5			\n\t		vmovaps		-0x060(%%rdx)	,%%ymm13			\n\t"\
		"vmovaps		      (%%rcx)	,%%ymm2			\n\t		vmovaps		 0x080(%%rcx)	,%%ymm10			\n\t"\
		"vmovaps		 0x020(%%rdx)	,%%ymm7			\n\t		vmovaps		 0x0a0(%%rdx)	,%%ymm15			\n\t"\
		"vmovaps		 0x020(%%rcx)	,%%ymm3			\n\t		vmovaps		 0x0a0(%%rcx)	,%%ymm11			\n\t"\
		"vmovaps		      (%%rdx)	,%%ymm6			\n\t		vmovaps		 0x080(%%rdx)	,%%ymm14			\n\t"\
		"vsubpd			%%ymm2		,%%ymm0,%%ymm0		\n\t		vsubpd		%%ymm10		,%%ymm8 ,%%ymm8 		\n\t"\
		"vsubpd			%%ymm7		,%%ymm4,%%ymm4		\n\t		vsubpd		%%ymm15		,%%ymm12,%%ymm12		\n\t"\
		"vsubpd			%%ymm3		,%%ymm1,%%ymm1		\n\t		vsubpd		%%ymm11		,%%ymm9 ,%%ymm9 		\n\t"\
		"vsubpd			%%ymm6		,%%ymm5,%%ymm5		\n\t		vsubpd		%%ymm14		,%%ymm13,%%ymm13		\n\t"\
		"vaddpd			%%ymm2		,%%ymm2,%%ymm2		\n\t		vaddpd		%%ymm10		,%%ymm10,%%ymm10		\n\t"\
		"vaddpd			%%ymm7		,%%ymm7,%%ymm7		\n\t		vaddpd		%%ymm15		,%%ymm15,%%ymm15		\n\t"\
		"vaddpd			%%ymm3		,%%ymm3,%%ymm3		\n\t		vaddpd		%%ymm11		,%%ymm11,%%ymm11		\n\t"\
		"vaddpd			%%ymm6		,%%ymm6,%%ymm6		\n\t		vaddpd		%%ymm14		,%%ymm14,%%ymm14		\n\t"\
		"vaddpd			%%ymm0		,%%ymm2,%%ymm2		\n\t		vaddpd		%%ymm8 		,%%ymm10,%%ymm10		\n\t"\
		"vaddpd			%%ymm4		,%%ymm7,%%ymm7		\n\t		vaddpd		%%ymm12		,%%ymm15,%%ymm15		\n\t"\
		"vaddpd			%%ymm1		,%%ymm3,%%ymm3		\n\t		vaddpd		%%ymm9 		,%%ymm11,%%ymm11		\n\t"\
		"vaddpd			%%ymm5		,%%ymm6,%%ymm6		\n\t		vaddpd		%%ymm13		,%%ymm14,%%ymm14		\n\t"\
		"vmovaps		%%ymm0		,      (%%rcx)		\n\t		vmovaps		%%ymm8 		, 0x080(%%rcx)			\n\t"\
		"vmovaps		%%ymm4		,      (%%rdx)		\n\t		vmovaps		%%ymm12		, 0x080(%%rdx)			\n\t"\
		"vmovaps		%%ymm1		, 0x020(%%rcx)		\n\t		vmovaps		%%ymm9 		, 0x0a0(%%rcx)			\n\t"\
		"vmovaps		%%ymm5		,-0x0e0(%%rdx)		\n\t		vmovaps		%%ymm13		,-0x060(%%rdx)			\n\t"\
		"vmovaps		%%ymm2		,-0x100(%%rcx)		\n\t		vmovaps		%%ymm10		,-0x080(%%rcx)			\n\t"\
		"vmovaps		%%ymm7		,-0x100(%%rdx)		\n\t		vmovaps		%%ymm15		,-0x080(%%rdx)			\n\t"\
		"vmovaps		%%ymm3		,-0x0e0(%%rcx)		\n\t		vmovaps		%%ymm11		,-0x060(%%rcx)			\n\t"\
		"vmovaps		%%ymm6		, 0x020(%%rdx)		\n\t		vmovaps		%%ymm14		, 0x0a0(%%rdx)			\n\t"\
		"/***************************************************************************************/\n\t"\
		"/* Now do eight more radix-4 transforms, including the internal and external twiddles: */\n\t"\
		"/***************************************************************************************/\n\t"\
		/* Using the upper block(s) of the main array for temp-storage in the section below led to a nasty AVX bug to track down: */\
		/* In fermat-mod mode the 4 block addresses in ascending order are add0,1,2,3 with no 'gaps' between blocks, whereas for */\
		/* mersenne-mod the addresses in asc. order are add0,2,3,1 with a gap between contiguous-data-block pairs 0,2 and 3,1. Thus */\
		/* for fermat-mod we need [add2] as the base-address of the 'high-half' block for temp-storage; for mersenne-mod we need [add3]. */\
		/* In both cases we have that add2 < add3 so instead use (add2 - add1): > 0 for fermat-mod, < 0 for mersenne - to differentiate: */\
		"movq	%[__add2],%%rsi		\n\t"/* destroyable copy of add2 */\
		"movq	%[__add2],%%rbx		\n\t"\
		"subq	%[__add1],%%rsi		\n\t"/* rsi = (add2 - add1); if this yields a borrow (i.e. sets CF) it's mersenne, else fermat. */\
		"cmovcq %[__add3],%%rbx	\n\t" /* if CF set (i.e. h > l), copy source [add3] into dest (rbx), else leave dest = [add2]. */\
		"\n\t"\
		"/************************************************/	\n\t		/************************************************/	\n\t"\
		"/*...Block 2: t02,t12,t22,t32 -> r10,14,12,16:  */	\n\t		/*...Block 6: t0A,t1A,t2A,t3A -> r18,1C,1A,1E:  */	\n\t"\
		"/************************************************/	\n\t		/************************************************/	\n\t"\
		"movq		%[__isrt2]		,%%rsi					\n\t"\
		"movq		%[__r10]		,%%rax	/* base-addr in rcol = c05/r18, so rax/rdi offset +0x100 vs lcol */\n\t"\
		"movq		%%rsi			,%%rcx					\n\t"\
		"movq		%%rsi			,%%rdx					\n\t"\
		"movq		%[__c01]		,%%rdi					\n\t"\
		"addq		$0x020			,%%rcx	/* cc0 */		\n\t"\
		"addq		$0x060			,%%rdx	/* cc1 */		\n\t"\
		"vmovaps		 0x040(%%rax)	,%%ymm4				\n\t		vmovaps		 0x140(%%rax)	,%%ymm12				\n\t"\
		"vmovaps		 0x0c0(%%rax)	,%%ymm0				\n\t		vmovaps		 0x1c0(%%rax)	,%%ymm8 				\n\t"\
		"vmovaps		 0x060(%%rax)	,%%ymm5				\n\t		vmovaps		 0x160(%%rax)	,%%ymm13				\n\t"\
		"vmovaps		 0x0e0(%%rax)	,%%ymm1				\n\t		vmovaps		 0x1e0(%%rax)	,%%ymm9 				\n\t"\
		"vmovaps		 0x040(%%rax)	,%%ymm6				\n\t		vmovaps		 0x140(%%rax)	,%%ymm14				\n\t"\
		"vmovaps		 0x0c0(%%rax)	,%%ymm2				\n\t		vmovaps		 0x1c0(%%rax)	,%%ymm10				\n\t"\
		"vmovaps		 0x060(%%rax)	,%%ymm7				\n\t		vmovaps		 0x160(%%rax)	,%%ymm15				\n\t"\
		"vmovaps		 0x0e0(%%rax)	,%%ymm3				\n\t		vmovaps		 0x1e0(%%rax)	,%%ymm11				\n\t"\
		"vmulpd		      (%%rdx)	,%%ymm4,%%ymm4			\n\t		vmulpd		 0x060(%%rdx)	,%%ymm12,%%ymm12			\n\t"\
		"vmulpd		 0x040(%%rdx)	,%%ymm0,%%ymm0			\n\t		vmulpd		      (%%rdx)	,%%ymm8 ,%%ymm8 			\n\t"\
		"vmulpd		      (%%rdx)	,%%ymm5,%%ymm5			\n\t		vmulpd		 0x060(%%rdx)	,%%ymm13,%%ymm13			\n\t"\
		"vmulpd		 0x040(%%rdx)	,%%ymm1,%%ymm1			\n\t		vmulpd		      (%%rdx)	,%%ymm9 ,%%ymm9 			\n\t"\
		"vmulpd		 0x020(%%rdx)	,%%ymm6,%%ymm6			\n\t		vmulpd		 0x040(%%rdx)	,%%ymm14,%%ymm14			\n\t"\
		"vmulpd		 0x060(%%rdx)	,%%ymm2,%%ymm2			\n\t		vmulpd		 0x020(%%rdx)	,%%ymm10,%%ymm10			\n\t"\
		"vmulpd		 0x020(%%rdx)	,%%ymm7,%%ymm7			\n\t		vmulpd		 0x040(%%rdx)	,%%ymm15,%%ymm15			\n\t"\
		"vmulpd		 0x060(%%rdx)	,%%ymm3,%%ymm3			\n\t		vmulpd		 0x020(%%rdx)	,%%ymm11,%%ymm11			\n\t"\
		"vsubpd		%%ymm6			,%%ymm5,%%ymm5			\n\t		vsubpd		%%ymm14			,%%ymm13,%%ymm13			\n\t"\
		"vsubpd		%%ymm2			,%%ymm1,%%ymm1			\n\t		vaddpd		%%ymm10			,%%ymm9 ,%%ymm9 			\n\t"\
		"vaddpd		%%ymm7			,%%ymm4,%%ymm4			\n\t		vaddpd		%%ymm15			,%%ymm12,%%ymm12			\n\t"\
		"vaddpd		%%ymm3			,%%ymm0,%%ymm0			\n\t		vsubpd		%%ymm11			,%%ymm8 ,%%ymm8 			\n\t"\
		"vmovaps		%%ymm5		,%%ymm7					\n\t		vmovaps		%%ymm13		,%%ymm15					\n\t"\
		"vmovaps		%%ymm4		,%%ymm6					\n\t		vmovaps		%%ymm12		,%%ymm14					\n\t"\
		"vsubpd		%%ymm0			,%%ymm4,%%ymm4			\n\t		vaddpd		%%ymm8 			,%%ymm12,%%ymm12			\n\t"\
		"vsubpd		%%ymm1			,%%ymm5,%%ymm5			\n\t		vaddpd		%%ymm9 			,%%ymm13,%%ymm13			\n\t"\
		"vaddpd		%%ymm0			,%%ymm6,%%ymm6			\n\t		vsubpd		%%ymm8 			,%%ymm14,%%ymm14			\n\t"\
		"vaddpd		%%ymm1			,%%ymm7,%%ymm7			\n\t		vsubpd		%%ymm9 			,%%ymm15,%%ymm15			\n\t"\
		"vmovaps		 0x080(%%rax)	,%%ymm2				\n\t		vmovaps		 0x180(%%rax)	,%%ymm10				\n\t"\
		"vmovaps		 0x0a0(%%rax)	,%%ymm3				\n\t		vmovaps		 0x1a0(%%rax)	,%%ymm11				\n\t"\
		"vmovaps		 0x080(%%rax)	,%%ymm0				\n\t		vmovaps		 0x180(%%rax)	,%%ymm8 				\n\t"\
		"vmovaps		 0x0a0(%%rax)	,%%ymm1				\n\t		vmovaps		 0x1a0(%%rax)	,%%ymm9 				\n\t"\
		"vmulpd		      (%%rcx)	,%%ymm2,%%ymm2			\n\t		vmulpd		 0x020(%%rcx)	,%%ymm10,%%ymm10			\n\t"\
		"vmulpd		 0x020(%%rcx)	,%%ymm1,%%ymm1			\n\t		vmulpd		      (%%rcx)	,%%ymm9 ,%%ymm9 			\n\t"\
		"vmulpd		      (%%rcx)	,%%ymm3,%%ymm3			\n\t		vmulpd		 0x020(%%rcx)	,%%ymm11,%%ymm11			\n\t"\
		"vmulpd		 0x020(%%rcx)	,%%ymm0,%%ymm0			\n\t		vmulpd		      (%%rcx)	,%%ymm8 ,%%ymm8 			\n\t"\
		"vaddpd		%%ymm1			,%%ymm2,%%ymm2			\n\t		vsubpd		%%ymm9 			,%%ymm10,%%ymm10			\n\t"\
		"vsubpd		%%ymm0			,%%ymm3,%%ymm3			\n\t		vaddpd		%%ymm8 			,%%ymm11,%%ymm11			\n\t"\
		"vmovaps		      (%%rax)	,%%ymm0				\n\t		vmovaps		 0x100(%%rax)	,%%ymm8 				\n\t"\
		"vmovaps		 0x020(%%rax)	,%%ymm1				\n\t		vmovaps		 0x120(%%rax)	,%%ymm9 				\n\t"\
		"vsubpd		%%ymm2			,%%ymm0,%%ymm0			\n\t		vsubpd		%%ymm10			,%%ymm8 ,%%ymm8 			\n\t"\
		"vsubpd		%%ymm3			,%%ymm1,%%ymm1			\n\t		vsubpd		%%ymm11			,%%ymm9 ,%%ymm9 			\n\t"\
		"vaddpd		%%ymm2			,%%ymm2,%%ymm2			\n\t		vaddpd		%%ymm10			,%%ymm10,%%ymm10			\n\t"\
		"vaddpd		%%ymm3			,%%ymm3,%%ymm3			\n\t		vaddpd		%%ymm11			,%%ymm11,%%ymm11			\n\t"\
		"vaddpd		%%ymm0			,%%ymm2,%%ymm2			\n\t		vaddpd		%%ymm8 			,%%ymm10,%%ymm10			\n\t"\
		"vaddpd		%%ymm1			,%%ymm3,%%ymm3			\n\t		vaddpd		%%ymm9 			,%%ymm11,%%ymm11			\n\t"\
		"vsubpd		%%ymm6			,%%ymm2,%%ymm2			\n\t		vsubpd		%%ymm14			,%%ymm8 ,%%ymm8 			\n\t"\
		"vsubpd		%%ymm7			,%%ymm3,%%ymm3			\n\t		vsubpd		%%ymm15			,%%ymm9 ,%%ymm9 			\n\t"\
		"vaddpd		%%ymm6			,%%ymm6,%%ymm6			\n\t		vaddpd		%%ymm14			,%%ymm14,%%ymm14			\n\t"\
		"vaddpd		%%ymm7			,%%ymm7,%%ymm7			\n\t		vaddpd		%%ymm15			,%%ymm15,%%ymm15			\n\t"\
		"vaddpd		%%ymm2			,%%ymm6,%%ymm6			\n\t		vaddpd		%%ymm8 			,%%ymm14,%%ymm14			\n\t"\
		"vaddpd		%%ymm3			,%%ymm7,%%ymm7			\n\t		vaddpd		%%ymm9 			,%%ymm15,%%ymm15			\n\t"\
		"vmovaps		%%ymm2		, 0x040(%%rax)			\n\t		vmovaps		%%ymm8 		, 0x140(%%rax)			\n\t"\
		"vmovaps		%%ymm3		, 0x060(%%rax)			\n\t		vmovaps		%%ymm9 		, 0x160(%%rax)			\n\t"\
		"vmovaps		%%ymm6		,%%ymm2					\n\t		vmovaps		%%ymm14		,%%ymm8 					\n\t"\
		"vmovaps		%%ymm7		,%%ymm3					\n\t		vmovaps		%%ymm15		,%%ymm9 					\n\t"\
		"vmulpd		      (%%rdi)	,%%ymm6,%%ymm6			\n\t		vmulpd		 0x100(%%rdi)	,%%ymm14,%%ymm14			\n\t"\
		"vmulpd		      (%%rdi)	,%%ymm7,%%ymm7			\n\t		vmulpd		 0x100(%%rdi)	,%%ymm15,%%ymm15			\n\t"\
		"vmulpd		 0x020(%%rdi)	,%%ymm2,%%ymm2			\n\t		vmulpd		 0x120(%%rdi)	,%%ymm8 ,%%ymm8 			\n\t"\
		"vmulpd		 0x020(%%rdi)	,%%ymm3,%%ymm3			\n\t		vmulpd		 0x120(%%rdi)	,%%ymm9 ,%%ymm9 			\n\t"\
		"vsubpd		%%ymm2			,%%ymm7,%%ymm7			\n\t		vsubpd		%%ymm8 			,%%ymm15,%%ymm15			\n\t"\
		"vaddpd		%%ymm3			,%%ymm6,%%ymm6			\n\t		vaddpd		%%ymm9 			,%%ymm14,%%ymm14			\n\t"\
		"vmovaps		%%ymm7		, 0x020(%%rbx)			\n\t		vmovaps		%%ymm15		, 0x0a0(%%rbx)			\n\t"\
		"vmovaps		%%ymm6		,      (%%rbx)			\n\t		vmovaps		%%ymm14		, 0x080(%%rbx)			\n\t"\
		"vmovaps		 0x040(%%rax)	,%%ymm6				\n\t		vmovaps		 0x140(%%rax)	,%%ymm14				\n\t"\
		"vmovaps		 0x060(%%rax)	,%%ymm7				\n\t		vmovaps		 0x160(%%rax)	,%%ymm15				\n\t"\
		"vmovaps		%%ymm6		,%%ymm2					\n\t		vmovaps		%%ymm14		,%%ymm8 					\n\t"\
		"vmovaps		%%ymm7		,%%ymm3					\n\t		vmovaps		%%ymm15		,%%ymm9 					\n\t"\
		"vmulpd		 0x040(%%rdi)	,%%ymm6,%%ymm6			\n\t		vmulpd		 0x140(%%rdi)	,%%ymm14,%%ymm14			\n\t"\
		"vmulpd		 0x040(%%rdi)	,%%ymm7,%%ymm7			\n\t		vmulpd		 0x140(%%rdi)	,%%ymm15,%%ymm15			\n\t"\
		"vmulpd		 0x060(%%rdi)	,%%ymm2,%%ymm2			\n\t		vmulpd		 0x160(%%rdi)	,%%ymm8 ,%%ymm8 			\n\t"\
		"vmulpd		 0x060(%%rdi)	,%%ymm3,%%ymm3			\n\t		vmulpd		 0x160(%%rdi)	,%%ymm9 ,%%ymm9 			\n\t"\
		"vsubpd		%%ymm2			,%%ymm7,%%ymm7			\n\t		vsubpd		%%ymm8 			,%%ymm15,%%ymm15			\n\t"\
		"vaddpd		%%ymm3			,%%ymm6,%%ymm6			\n\t		vaddpd		%%ymm9 			,%%ymm14,%%ymm14			\n\t"\
		"vmovaps		%%ymm7		, 0x220(%%rbx)			\n\t		vmovaps		%%ymm15		, 0x2a0(%%rbx)			\n\t"\
		"vmovaps		%%ymm6		, 0x200(%%rbx)			\n\t		vmovaps		%%ymm14		, 0x280(%%rbx)			\n\t"\
		"addq		$0x080			,%%rdi					\n\t"\
		"vsubpd		%%ymm5			,%%ymm0,%%ymm0			\n\t		vsubpd		%%ymm13			,%%ymm10,%%ymm10			\n\t"\
		"vsubpd		%%ymm4			,%%ymm1,%%ymm1			\n\t		vsubpd		%%ymm12			,%%ymm11,%%ymm11			\n\t"\
		"vaddpd		%%ymm5			,%%ymm5,%%ymm5			\n\t		vaddpd		%%ymm13			,%%ymm13,%%ymm13			\n\t"\
		"vaddpd		%%ymm4			,%%ymm4,%%ymm4			\n\t		vaddpd		%%ymm12			,%%ymm12,%%ymm12			\n\t"\
		"vaddpd		%%ymm0			,%%ymm5,%%ymm5			\n\t		vaddpd		%%ymm10			,%%ymm13,%%ymm13			\n\t"\
		"vaddpd		%%ymm1			,%%ymm4,%%ymm4			\n\t		vaddpd		%%ymm11			,%%ymm12,%%ymm12			\n\t"\
		"vmovaps		%%ymm5		,%%ymm2					\n\t		vmovaps		%%ymm13		,%%ymm8 					\n\t"\
		"vmovaps		%%ymm1		,%%ymm3					\n\t		vmovaps		%%ymm11		,%%ymm9 					\n\t"\
		"vmulpd		      (%%rdi)	,%%ymm5,%%ymm5			\n\t		vmulpd		 0x100(%%rdi)	,%%ymm13,%%ymm13			\n\t"\
		"vmulpd		      (%%rdi)	,%%ymm1,%%ymm1			\n\t		vmulpd		 0x100(%%rdi)	,%%ymm11,%%ymm11			\n\t"\
		"vmulpd		 0x020(%%rdi)	,%%ymm2,%%ymm2			\n\t		vmulpd		 0x120(%%rdi)	,%%ymm8 ,%%ymm8 			\n\t"\
		"vmulpd		 0x020(%%rdi)	,%%ymm3,%%ymm3			\n\t		vmulpd		 0x120(%%rdi)	,%%ymm9 ,%%ymm9 			\n\t"\
		"vsubpd		%%ymm2			,%%ymm1,%%ymm1			\n\t		vsubpd		%%ymm8 			,%%ymm11,%%ymm11			\n\t"\
		"vaddpd		%%ymm3			,%%ymm5,%%ymm5			\n\t		vaddpd		%%ymm9 			,%%ymm13,%%ymm13			\n\t"\
		"vmovaps		%%ymm1		, 0x120(%%rbx)			\n\t		vmovaps		%%ymm11		, 0x1a0(%%rbx)			\n\t"\
		"vmovaps		%%ymm5		, 0x100(%%rbx)			\n\t		vmovaps		%%ymm13		, 0x180(%%rbx)			\n\t"\
		"vmovaps		%%ymm0		,%%ymm2					\n\t		vmovaps		%%ymm10		,%%ymm8 					\n\t"\
		"vmovaps		%%ymm4		,%%ymm3					\n\t		vmovaps		%%ymm12		,%%ymm9 					\n\t"\
		"vmulpd		 0x040(%%rdi)	,%%ymm0,%%ymm0			\n\t		vmulpd		 0x140(%%rdi)	,%%ymm10,%%ymm10			\n\t"\
		"vmulpd		 0x040(%%rdi)	,%%ymm4,%%ymm4			\n\t		vmulpd		 0x140(%%rdi)	,%%ymm12,%%ymm12			\n\t"\
		"vmulpd		 0x060(%%rdi)	,%%ymm2,%%ymm2			\n\t		vmulpd		 0x160(%%rdi)	,%%ymm8 ,%%ymm8 			\n\t"\
		"vmulpd		 0x060(%%rdi)	,%%ymm3,%%ymm3			\n\t		vmulpd		 0x160(%%rdi)	,%%ymm9 ,%%ymm9 			\n\t"\
		"vsubpd		%%ymm2			,%%ymm4,%%ymm4			\n\t		vsubpd		%%ymm8 			,%%ymm12,%%ymm12			\n\t"\
		"vaddpd		%%ymm3			,%%ymm0,%%ymm0			\n\t		vaddpd		%%ymm9 			,%%ymm10,%%ymm10			\n\t"\
		"vmovaps		%%ymm4		, 0x320(%%rbx)			\n\t		vmovaps		%%ymm12		, 0x3a0(%%rbx)			\n\t"\
		"vmovaps		%%ymm0		, 0x300(%%rbx)			\n\t		vmovaps		%%ymm10		, 0x380(%%rbx)			\n\t"\
		"\n\t"\
		"/************************************************/	\n\t		/************************************************/	\n\t"\
		"/*...Block 4: t06,t16,t26,t36 -> r30,34,32,36:  */	\n\t		/*...Block 8: t0E,t1E,t2E,t3E -> r38,3C,3A,3E:  */	\n\t"\
		"/************************************************/	\n\t		/************************************************/	\n\t"\
		"addq		$0x400		,%%rax						\n\t		addq		$0x100		,%%rdi						\n\t"\
		"vmovaps		 0x040(%%rax)	,%%ymm4				\n\t		vmovaps		 0x140(%%rax)	,%%ymm12				\n\t"\
		"vmovaps		 0x0c0(%%rax)	,%%ymm0				\n\t		vmovaps		 0x1c0(%%rax)	,%%ymm8 				\n\t"\
		"vmovaps		 0x060(%%rax)	,%%ymm5				\n\t		vmovaps		 0x160(%%rax)	,%%ymm13				\n\t"\
		"vmovaps		 0x0e0(%%rax)	,%%ymm1				\n\t		vmovaps		 0x1e0(%%rax)	,%%ymm9 				\n\t"\
		"vmovaps		 0x040(%%rax)	,%%ymm6				\n\t		vmovaps		 0x140(%%rax)	,%%ymm14				\n\t"\
		"vmovaps		 0x0c0(%%rax)	,%%ymm2				\n\t		vmovaps		 0x1c0(%%rax)	,%%ymm10				\n\t"\
		"vmovaps		 0x060(%%rax)	,%%ymm7				\n\t		vmovaps		 0x160(%%rax)	,%%ymm15				\n\t"\
		"vmovaps		 0x0e0(%%rax)	,%%ymm3				\n\t		vmovaps		 0x1e0(%%rax)	,%%ymm11				\n\t"\
		"vmulpd		 0x040(%%rdx)	,%%ymm4,%%ymm4			\n\t		vmulpd		 0x020(%%rdx)	,%%ymm12,%%ymm12			\n\t"\
		"vmulpd		 0x020(%%rdx)	,%%ymm0,%%ymm0			\n\t		vmulpd		 0x060(%%rdx)	,%%ymm8 ,%%ymm8 			\n\t"\
		"vmulpd		 0x040(%%rdx)	,%%ymm5,%%ymm5			\n\t		vmulpd		 0x020(%%rdx)	,%%ymm13,%%ymm13			\n\t"\
		"vmulpd		 0x020(%%rdx)	,%%ymm1,%%ymm1			\n\t		vmulpd		 0x060(%%rdx)	,%%ymm9 ,%%ymm9 			\n\t"\
		"vmulpd		 0x060(%%rdx)	,%%ymm6,%%ymm6			\n\t		vmulpd		      (%%rdx)	,%%ymm14,%%ymm14			\n\t"\
		"vmulpd		      (%%rdx)	,%%ymm2,%%ymm2			\n\t		vmulpd		 0x040(%%rdx)	,%%ymm10,%%ymm10			\n\t"\
		"vmulpd		 0x060(%%rdx)	,%%ymm7,%%ymm7			\n\t		vmulpd		      (%%rdx)	,%%ymm15,%%ymm15			\n\t"\
		"vmulpd		      (%%rdx)	,%%ymm3,%%ymm3			\n\t		vmulpd		 0x040(%%rdx)	,%%ymm11,%%ymm11			\n\t"\
		"vsubpd		%%ymm6			,%%ymm5,%%ymm5			\n\t		vsubpd		%%ymm14			,%%ymm13,%%ymm13			\n\t"\
		"vaddpd		%%ymm2			,%%ymm1,%%ymm1			\n\t		vsubpd		%%ymm10			,%%ymm9 ,%%ymm9 			\n\t"\
		"vaddpd		%%ymm7			,%%ymm4,%%ymm4			\n\t		vaddpd		%%ymm15			,%%ymm12,%%ymm12			\n\t"\
		"vsubpd		%%ymm3			,%%ymm0,%%ymm0			\n\t		vaddpd		%%ymm11			,%%ymm8 ,%%ymm8 			\n\t"\
		"vmovaps		%%ymm5		,%%ymm7					\n\t		vmovaps		%%ymm13		,%%ymm15					\n\t"\
		"vmovaps		%%ymm4		,%%ymm6					\n\t		vmovaps		%%ymm12		,%%ymm14					\n\t"\
		"vaddpd		%%ymm0			,%%ymm4,%%ymm4			\n\t		vaddpd		%%ymm8 			,%%ymm12,%%ymm12			\n\t"\
		"vaddpd		%%ymm1			,%%ymm5,%%ymm5			\n\t		vaddpd		%%ymm9 			,%%ymm13,%%ymm13			\n\t"\
		"vsubpd		%%ymm0			,%%ymm6,%%ymm6			\n\t		vsubpd		%%ymm8 			,%%ymm14,%%ymm14			\n\t"\
		"vsubpd		%%ymm1			,%%ymm7,%%ymm7			\n\t		vsubpd		%%ymm9 			,%%ymm15,%%ymm15			\n\t"\
		"vmovaps		 0x080(%%rax)	,%%ymm2				\n\t		vmovaps		 0x180(%%rax)	,%%ymm10				\n\t"\
		"vmovaps		 0x0a0(%%rax)	,%%ymm3				\n\t		vmovaps		 0x1a0(%%rax)	,%%ymm11				\n\t"\
		"vmovaps		 0x080(%%rax)	,%%ymm0				\n\t		vmovaps		 0x180(%%rax)	,%%ymm8 				\n\t"\
		"vmovaps		 0x0a0(%%rax)	,%%ymm1				\n\t		vmovaps		 0x1a0(%%rax)	,%%ymm9 				\n\t"\
		"vmulpd		 0x020(%%rcx)	,%%ymm2,%%ymm2			\n\t		vmulpd		      (%%rcx)	,%%ymm10,%%ymm10			\n\t"\
		"vmulpd		      (%%rcx)	,%%ymm1,%%ymm1			\n\t		vmulpd		 0x020(%%rcx)	,%%ymm9 ,%%ymm9 			\n\t"\
		"vmulpd		 0x020(%%rcx)	,%%ymm3,%%ymm3			\n\t		vmulpd		      (%%rcx)	,%%ymm11,%%ymm11			\n\t"\
		"vmulpd		      (%%rcx)	,%%ymm0,%%ymm0			\n\t		vmulpd		 0x020(%%rcx)	,%%ymm8 ,%%ymm8 			\n\t"\
		"vaddpd		%%ymm1			,%%ymm2,%%ymm2			\n\t		vsubpd		%%ymm9 			,%%ymm10,%%ymm10			\n\t"\
		"vsubpd		%%ymm0			,%%ymm3,%%ymm3			\n\t		vaddpd		%%ymm8 			,%%ymm11,%%ymm11			\n\t"\
		"vmovaps		      (%%rax)	,%%ymm0				\n\t		vmovaps		 0x100(%%rax)	,%%ymm8 				\n\t"\
		"vmovaps		 0x020(%%rax)	,%%ymm1				\n\t		vmovaps		 0x120(%%rax)	,%%ymm9 				\n\t"\
		"vsubpd		%%ymm2			,%%ymm0,%%ymm0			\n\t		vsubpd		%%ymm10			,%%ymm8 ,%%ymm8 			\n\t"\
		"vsubpd		%%ymm3			,%%ymm1,%%ymm1			\n\t		vsubpd		%%ymm11			,%%ymm9 ,%%ymm9 			\n\t"\
		"vaddpd		%%ymm2			,%%ymm2,%%ymm2			\n\t		vaddpd		%%ymm10			,%%ymm10,%%ymm10			\n\t"\
		"vaddpd		%%ymm3			,%%ymm3,%%ymm3			\n\t		vaddpd		%%ymm11			,%%ymm11,%%ymm11			\n\t"\
		"vaddpd		%%ymm0			,%%ymm2,%%ymm2			\n\t		vaddpd		%%ymm8 			,%%ymm10,%%ymm10			\n\t"\
		"vaddpd		%%ymm1			,%%ymm3,%%ymm3			\n\t		vaddpd		%%ymm9 			,%%ymm11,%%ymm11			\n\t"\
		"addq		$0x080		,%%rdi						\n\t"\
		"vsubpd		%%ymm6			,%%ymm2,%%ymm2			\n\t		vsubpd		%%ymm14			,%%ymm8 ,%%ymm8 			\n\t"\
		"vsubpd		%%ymm7			,%%ymm3,%%ymm3			\n\t		vsubpd		%%ymm15			,%%ymm9 ,%%ymm9 			\n\t"\
		"vaddpd		%%ymm6			,%%ymm6,%%ymm6			\n\t		vaddpd		%%ymm14			,%%ymm14,%%ymm14			\n\t"\
		"vaddpd		%%ymm7			,%%ymm7,%%ymm7			\n\t		vaddpd		%%ymm15			,%%ymm15,%%ymm15			\n\t"\
		"vaddpd		%%ymm2			,%%ymm6,%%ymm6			\n\t		vaddpd		%%ymm8 			,%%ymm14,%%ymm14			\n\t"\
		"vaddpd		%%ymm3			,%%ymm7,%%ymm7			\n\t		vaddpd		%%ymm9 			,%%ymm15,%%ymm15			\n\t"\
		"vmovaps		%%ymm2		, 0x040(%%rax)			\n\t		vmovaps		%%ymm8 		, 0x140(%%rax)			\n\t"\
		"vmovaps		%%ymm3		, 0x060(%%rax)			\n\t		vmovaps		%%ymm9 		, 0x160(%%rax)			\n\t"\
		"vmovaps		%%ymm6		,%%ymm2					\n\t		vmovaps		%%ymm14		,%%ymm8 					\n\t"\
		"vmovaps		%%ymm7		,%%ymm3					\n\t		vmovaps		%%ymm15		,%%ymm9 					\n\t"\
		"vmulpd		      (%%rdi)	,%%ymm6,%%ymm6			\n\t		vmulpd		 0x100(%%rdi)	,%%ymm14,%%ymm14			\n\t"\
		"vmulpd		      (%%rdi)	,%%ymm7,%%ymm7			\n\t		vmulpd		 0x100(%%rdi)	,%%ymm15,%%ymm15			\n\t"\
		"vmulpd		 0x020(%%rdi)	,%%ymm2,%%ymm2			\n\t		vmulpd		 0x120(%%rdi)	,%%ymm8 ,%%ymm8 			\n\t"\
		"vmulpd		 0x020(%%rdi)	,%%ymm3,%%ymm3			\n\t		vmulpd		 0x120(%%rdi)	,%%ymm9 ,%%ymm9 			\n\t"\
		"vsubpd		%%ymm2			,%%ymm7,%%ymm7			\n\t		vsubpd		%%ymm8 			,%%ymm15,%%ymm15			\n\t"\
		"vaddpd		%%ymm3			,%%ymm6,%%ymm6			\n\t		vaddpd		%%ymm9 			,%%ymm14,%%ymm14			\n\t"\
		"vmovaps		%%ymm7		, 0x060(%%rbx)			\n\t		vmovaps		%%ymm15		, 0x0e0(%%rbx)			\n\t"\
		"vmovaps		%%ymm6		, 0x040(%%rbx)			\n\t		vmovaps		%%ymm14		, 0x0c0(%%rbx)			\n\t"\
		"vmovaps		 0x040(%%rax)	,%%ymm6				\n\t		vmovaps		 0x140(%%rax)	,%%ymm14				\n\t"\
		"vmovaps		 0x060(%%rax)	,%%ymm7				\n\t		vmovaps		 0x160(%%rax)	,%%ymm15				\n\t"\
		"vmovaps		%%ymm6		,%%ymm2					\n\t		vmovaps		%%ymm14		,%%ymm8 					\n\t"\
		"vmovaps		%%ymm7		,%%ymm3					\n\t		vmovaps		%%ymm15		,%%ymm9 					\n\t"\
		"vmulpd		 0x040(%%rdi)	,%%ymm6,%%ymm6			\n\t		vmulpd		 0x140(%%rdi)	,%%ymm14,%%ymm14			\n\t"\
		"vmulpd		 0x040(%%rdi)	,%%ymm7,%%ymm7			\n\t		vmulpd		 0x140(%%rdi)	,%%ymm15,%%ymm15			\n\t"\
		"vmulpd		 0x060(%%rdi)	,%%ymm2,%%ymm2			\n\t		vmulpd		 0x160(%%rdi)	,%%ymm8 ,%%ymm8 			\n\t"\
		"vmulpd		 0x060(%%rdi)	,%%ymm3,%%ymm3			\n\t		vmulpd		 0x160(%%rdi)	,%%ymm9 ,%%ymm9 			\n\t"\
		"vsubpd		%%ymm2			,%%ymm7,%%ymm7			\n\t		vsubpd		%%ymm8 			,%%ymm15,%%ymm15			\n\t"\
		"vaddpd		%%ymm3			,%%ymm6,%%ymm6			\n\t		vaddpd		%%ymm9 			,%%ymm14,%%ymm14			\n\t"\
		"vmovaps		%%ymm7		, 0x260(%%rbx)			\n\t		vmovaps		%%ymm15		, 0x2e0(%%rbx)			\n\t"\
		"vmovaps		%%ymm6		, 0x240(%%rbx)			\n\t		vmovaps		%%ymm14		, 0x2c0(%%rbx)			\n\t"\
		"addq		$0x080		,%%rdi						\n\t"\
		"vsubpd		%%ymm5			,%%ymm0,%%ymm0			\n\t		vsubpd		%%ymm13			,%%ymm10,%%ymm10			\n\t"\
		"vsubpd		%%ymm4			,%%ymm1,%%ymm1			\n\t		vsubpd		%%ymm12			,%%ymm11,%%ymm11			\n\t"\
		"vaddpd		%%ymm5			,%%ymm5,%%ymm5			\n\t		vaddpd		%%ymm13			,%%ymm13,%%ymm13			\n\t"\
		"vaddpd		%%ymm4			,%%ymm4,%%ymm4			\n\t		vaddpd		%%ymm12			,%%ymm12,%%ymm12			\n\t"\
		"vaddpd		%%ymm0			,%%ymm5,%%ymm5			\n\t		vaddpd		%%ymm10			,%%ymm13,%%ymm13			\n\t"\
		"vaddpd		%%ymm1			,%%ymm4,%%ymm4			\n\t		vaddpd		%%ymm11			,%%ymm12,%%ymm12			\n\t"\
		"vmovaps		%%ymm5		,%%ymm2					\n\t		vmovaps		%%ymm13		,%%ymm8 					\n\t"\
		"vmovaps		%%ymm1		,%%ymm3					\n\t		vmovaps		%%ymm11		,%%ymm9 					\n\t"\
		"vmulpd		      (%%rdi)	,%%ymm5,%%ymm5			\n\t		vmulpd		 0x100(%%rdi)	,%%ymm13,%%ymm13			\n\t"\
		"vmulpd		      (%%rdi)	,%%ymm1,%%ymm1			\n\t		vmulpd		 0x100(%%rdi)	,%%ymm11,%%ymm11			\n\t"\
		"vmulpd		 0x020(%%rdi)	,%%ymm2,%%ymm2			\n\t		vmulpd		 0x120(%%rdi)	,%%ymm8 ,%%ymm8 			\n\t"\
		"vmulpd		 0x020(%%rdi)	,%%ymm3,%%ymm3			\n\t		vmulpd		 0x120(%%rdi)	,%%ymm9 ,%%ymm9 			\n\t"\
		"vsubpd		%%ymm2			,%%ymm1,%%ymm1			\n\t		vsubpd		%%ymm8 			,%%ymm11,%%ymm11			\n\t"\
		"vaddpd		%%ymm3			,%%ymm5,%%ymm5			\n\t		vaddpd		%%ymm9 			,%%ymm13,%%ymm13			\n\t"\
		"vmovaps		%%ymm1		, 0x160(%%rbx)			\n\t		vmovaps		%%ymm11		, 0x1e0(%%rbx)			\n\t"\
		"vmovaps		%%ymm5		, 0x140(%%rbx)			\n\t		vmovaps		%%ymm13		, 0x1c0(%%rbx)			\n\t"\
		"vmovaps		%%ymm0		,%%ymm2					\n\t		vmovaps		%%ymm10		,%%ymm8 					\n\t"\
		"vmovaps		%%ymm4		,%%ymm3					\n\t		vmovaps		%%ymm12		,%%ymm9 					\n\t"\
		"vmulpd		 0x040(%%rdi)	,%%ymm0,%%ymm0			\n\t		vmulpd		 0x140(%%rdi)	,%%ymm10,%%ymm10			\n\t"\
		"vmulpd		 0x040(%%rdi)	,%%ymm4,%%ymm4			\n\t		vmulpd		 0x140(%%rdi)	,%%ymm12,%%ymm12			\n\t"\
		"vmulpd		 0x060(%%rdi)	,%%ymm2,%%ymm2			\n\t		vmulpd		 0x160(%%rdi)	,%%ymm8 ,%%ymm8 			\n\t"\
		"vmulpd		 0x060(%%rdi)	,%%ymm3,%%ymm3			\n\t		vmulpd		 0x160(%%rdi)	,%%ymm9 ,%%ymm9 			\n\t"\
		"vsubpd		%%ymm2			,%%ymm4,%%ymm4			\n\t		vsubpd			%%ymm8 		,%%ymm12,%%ymm12			\n\t"\
		"vaddpd		%%ymm3			,%%ymm0,%%ymm0			\n\t		vaddpd			%%ymm9 		,%%ymm10,%%ymm10			\n\t"\
		"vmovaps		%%ymm4		, 0x360(%%rbx)			\n\t		vmovaps		%%ymm12		, 0x3e0(%%rbx)			\n\t"\
		"vmovaps		%%ymm0		, 0x340(%%rbx)			\n\t		vmovaps		%%ymm10		, 0x3c0(%%rbx)			\n\t"\
		"\n\t"\
		"/************************************************/	\n\t		/************************************************/	\n\t"\
		"/*...Block 1: t00,t10,t20,t30 -> r00,04,02,06:  */	\n\t		/*...Block 5: t08,t18,t28,t38 -> r08,0C,0A,0E:  */	\n\t"\
		"/************************************************/	\n\t		/************************************************/	\n\t"\
		"movq	%[__r00],%%rdx	/* base-addr in rcol = r08, so rdx+0x100 in rcol */	\n\t	vmovaps	(%%rsi),%%ymm10	/* isrt2 */	\n\t"\
		"vmovaps		      (%%rdx)	,%%ymm0				\n\t		vmovaps		 0x140(%%rdx)	,%%ymm12				\n\t"\
		"vmovaps		 0x020(%%rdx)	,%%ymm1				\n\t		vmovaps		 0x160(%%rdx)	,%%ymm13				\n\t"\
		"vmovaps		 0x080(%%rdx)	,%%ymm2				\n\t		vmovaps		 0x1c0(%%rdx)	,%%ymm8 				\n\t"\
		"vmovaps		 0x0a0(%%rdx)	,%%ymm3				\n\t		vmovaps		 0x1e0(%%rdx)	,%%ymm9 				\n\t"\
		"vsubpd		 0x080(%%rdx)	,%%ymm0,%%ymm0			\n\t		vaddpd		 0x160(%%rdx)	,%%ymm12,%%ymm12			\n\t"\
		"vsubpd		 0x0a0(%%rdx)	,%%ymm1,%%ymm1			\n\t		vsubpd		 0x140(%%rdx)	,%%ymm13,%%ymm13			\n\t"\
		"vaddpd		      (%%rdx)	,%%ymm2,%%ymm2			\n\t		vsubpd		 0x1e0(%%rdx)	,%%ymm8 ,%%ymm8 			\n\t"\
		"vaddpd		 0x020(%%rdx)	,%%ymm3,%%ymm3			\n\t		vaddpd		 0x1c0(%%rdx)	,%%ymm9 ,%%ymm9 			\n\t"\
		"vmovaps		 0x040(%%rdx)	,%%ymm4				\n\t		vmulpd			%%ymm10		,%%ymm12,%%ymm12			\n\t"\
		"vmovaps		 0x060(%%rdx)	,%%ymm5				\n\t		vmulpd			%%ymm10		,%%ymm13,%%ymm13			\n\t"\
		"vmovaps		 0x0c0(%%rdx)	,%%ymm6				\n\t		vmulpd			%%ymm10		,%%ymm8 ,%%ymm8 			\n\t"\
		"vmovaps		 0x0e0(%%rdx)	,%%ymm7				\n\t		vmulpd			%%ymm10		,%%ymm9 ,%%ymm9 			\n\t"\
		"vsubpd		 0x0c0(%%rdx)	,%%ymm4,%%ymm4			\n\t		vmovaps		%%ymm12		,%%ymm14					\n\t"\
		"vsubpd		 0x0e0(%%rdx)	,%%ymm5,%%ymm5			\n\t		vmovaps		%%ymm13		,%%ymm15					\n\t"\
		"vaddpd		 0x040(%%rdx)	,%%ymm6,%%ymm6			\n\t		vsubpd			%%ymm8 		,%%ymm12,%%ymm12			\n\t"\
		"vaddpd		 0x060(%%rdx)	,%%ymm7,%%ymm7			\n\t		vsubpd			%%ymm9 		,%%ymm13,%%ymm13			\n\t"\
		"/* base-twiddle in l/rcol = c00/c04, so rcx+0x100 in rcol */	vaddpd			%%ymm8 		,%%ymm14,%%ymm14			\n\t"\
		"movq		%[__c10]	,%%rcx						\n\t		vaddpd			%%ymm9 		,%%ymm15,%%ymm15			\n\t"\
		"vaddpd			%%ymm6		,%%ymm2,%%ymm2			\n\t		vmovaps		 0x100(%%rdx)	,%%ymm8 				\n\t"\
		"vaddpd			%%ymm7		,%%ymm3,%%ymm3			\n\t		vmovaps		 0x120(%%rdx)	,%%ymm9 				\n\t"\
		"vmovaps		%%ymm2		,      (%%rdx)			\n\t		vmovaps		 0x180(%%rdx)	,%%ymm10				\n\t"\
		"vmovaps		%%ymm3		, 0x020(%%rdx)			\n\t		vmovaps		 0x1a0(%%rdx)	,%%ymm11				\n\t"\
		"vaddpd			%%ymm6		,%%ymm6,%%ymm6			\n\t		vsubpd		 0x1a0(%%rdx)	,%%ymm8 ,%%ymm8 			\n\t"\
		"vaddpd			%%ymm7		,%%ymm7,%%ymm7			\n\t		vsubpd		 0x180(%%rdx)	,%%ymm9 ,%%ymm9 			\n\t"\
		"vsubpd			%%ymm6		,%%ymm2,%%ymm2			\n\t		vaddpd		 0x100(%%rdx)	,%%ymm11,%%ymm11			\n\t"\
		"vsubpd			%%ymm7		,%%ymm3,%%ymm3			\n\t		vaddpd		 0x120(%%rdx)	,%%ymm10,%%ymm10			\n\t"\
		"vmovaps		%%ymm2		,%%ymm6					\n\t		vsubpd			%%ymm12		,%%ymm11,%%ymm11			\n\t"\
		"vmovaps		%%ymm3		,%%ymm7					\n\t		vsubpd			%%ymm13		,%%ymm9 ,%%ymm9 			\n\t"\
		"vmulpd		      (%%rcx)	,%%ymm2,%%ymm2	/* c10 */	\n\t	vaddpd			%%ymm12		,%%ymm12,%%ymm12			\n\t"\
		"vmulpd		      (%%rcx)	,%%ymm3,%%ymm3			\n\t		vaddpd			%%ymm13		,%%ymm13,%%ymm13			\n\t"\
		"vmulpd		 0x020(%%rcx)	,%%ymm6,%%ymm6			\n\t		vaddpd			%%ymm11		,%%ymm12,%%ymm12			\n\t"\
		"vmulpd		 0x020(%%rcx)	,%%ymm7,%%ymm7			\n\t		vaddpd			%%ymm9 		,%%ymm13,%%ymm13			\n\t"\
		"vsubpd			%%ymm6		,%%ymm3,%%ymm3			\n\t		vmovaps		%%ymm11		, 0x140(%%rdx)			\n\t"\
		"vaddpd			%%ymm7		,%%ymm2,%%ymm2			\n\t		vmovaps		%%ymm9 		, 0x160(%%rdx)			\n\t"\
		"subq $0x40,%%rcx	/* put c00 in rcx to ease bookkeeping*/\n\t	vmovaps		%%ymm12		,%%ymm11					\n\t"\
		"/* add0,1 in rax,rbx; __r00 in rdx: */							vmovaps		%%ymm13		,%%ymm9 					\n\t"\
		"/* For each complex output octet, complex pairs having */		vmulpd		 0x100(%%rcx)	,%%ymm12,%%ymm12	/* c04 */\n\t"\
		"/* reads from offsets 0x0..,0x1..,0x2..,0x3.. go into  */		vmulpd		 0x100(%%rcx)	,%%ymm13,%%ymm13			\n\t"\
		"/* local-mem pairs rXY + 00/10, 04/14, 02/12, 06/16.   */		vmulpd		 0x120(%%rcx)	,%%ymm11,%%ymm11			\n\t"\
		"/* For 1st octet we read from offsets [0x2..,0x0..],   */		vmulpd		 0x120(%%rcx)	,%%ymm9 ,%%ymm9 			\n\t"\
		"/* [0x1..,0x3], other 3 octets use order [0,2],[1,3].  */		vsubpd			%%ymm11		,%%ymm13,%%ymm13			\n\t"\
		"vmovaps	0x220(%%rbx),%%ymm7						\n\t		vaddpd			%%ymm9 		,%%ymm12,%%ymm12			\n\t"\
		"vmovaps	0x200(%%rbx),%%ymm6						\n\t		vmovaps	0x0a0(%%rbx),%%ymm11						\n\t"\
		"vmovaps	%%ymm3,0x060(%%rdx)						\n\t		vmovaps	0x080(%%rbx),%%ymm9 						\n\t"\
		"vmovaps	%%ymm2,0x040(%%rdx)		\n\t/* r02,03 */			vmovaps	%%ymm13,0x120(%%rdx)						\n\t"\
		"vmovaps	%%ymm7,0x260(%%rdx)						\n\t		vmovaps	%%ymm12,0x100(%%rdx)			\n\t"/* r08,09 */\
		"vmovaps	%%ymm6,0x240(%%rdx)		\n\t/* r12,13 */			vmovaps	%%ymm11,0x320(%%rdx)						\n\t"\
		"vmovaps		 0x020(%%rdx)	,%%ymm3				\n\t		vmovaps	%%ymm9 ,0x300(%%rdx)			\n\t"/* r18,19 */\
		"vmovaps		      (%%rdx)	,%%ymm2				\n\t		vmovaps		0x140(%%rdx)	,%%ymm12				\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm7						\n\t		vmovaps		0x160(%%rdx)	,%%ymm13				\n\t"\
		"vmovaps	0x000(%%rbx),%%ymm6						\n\t		vmovaps		%%ymm12		,%%ymm11					\n\t"\
		"vmovaps	%%ymm3,0x020(%%rdx)						\n\t		vmovaps		%%ymm13		,%%ymm9 					\n\t"\
		"vmovaps	%%ymm2,     (%%rdx)		\n\t/* r00,01 */			vmulpd		 0x140(%%rcx)	,%%ymm12,%%ymm12	/* c14 */\n\t"\
		"vmovaps	%%ymm7,0x220(%%rdx)						\n\t		vmulpd		 0x140(%%rcx)	,%%ymm13,%%ymm13			\n\t"\
		"vmovaps	%%ymm6,0x200(%%rdx)		\n\t/* r10,11 */			vmulpd		 0x160(%%rcx)	,%%ymm11,%%ymm11			\n\t"\
		"vaddpd			%%ymm5		,%%ymm0,%%ymm0			\n\t		vmulpd		 0x160(%%rcx)	,%%ymm9 ,%%ymm9 			\n\t"\
		"vsubpd			%%ymm4		,%%ymm1,%%ymm1			\n\t		vsubpd			%%ymm11		,%%ymm13,%%ymm13			\n\t"\
		"vmovaps		%%ymm0		,%%ymm2					\n\t		vaddpd			%%ymm9 		,%%ymm12,%%ymm12			\n\t"\
		"vmovaps		%%ymm1		,%%ymm3					\n\t		vmovaps	0x2a0(%%rbx),%%ymm11						\n\t"\
		"vaddpd			%%ymm5		,%%ymm5,%%ymm5			\n\t		vmovaps	0x280(%%rbx),%%ymm9 						\n\t"\
		"vaddpd			%%ymm4		,%%ymm4,%%ymm4			\n\t		vmovaps	%%ymm13,0x160(%%rdx)						\n\t"\
		"vmovaps		%%ymm0		,%%ymm6					\n\t		vmovaps	%%ymm12,0x140(%%rdx)			\n\t"/* r0a,0b */\
		"vmovaps		%%ymm1		,%%ymm7					\n\t		vmovaps	%%ymm11,0x360(%%rdx)						\n\t"\
		"vmulpd		 0x080(%%rcx)	,%%ymm2,%%ymm2	/* c08 */	\n\t	vmovaps	%%ymm9 ,0x340(%%rdx)			\n\t"/* r1a,1b */\
		"vmulpd		 0x080(%%rcx)	,%%ymm3,%%ymm3			\n\t		vsubpd			%%ymm15		,%%ymm8 ,%%ymm8 			\n\t"\
		"vmulpd		 0x0a0(%%rcx)	,%%ymm6,%%ymm6			\n\t		vsubpd			%%ymm14		,%%ymm10,%%ymm10			\n\t"\
		"vmulpd		 0x0a0(%%rcx)	,%%ymm7,%%ymm7			\n\t		vaddpd			%%ymm15		,%%ymm15,%%ymm15			\n\t"\
		"vsubpd			%%ymm6		,%%ymm3,%%ymm3			\n\t		vaddpd			%%ymm14		,%%ymm14,%%ymm14			\n\t"\
		"vaddpd			%%ymm7		,%%ymm2,%%ymm2			\n\t		vaddpd			%%ymm8 		,%%ymm15,%%ymm15			\n\t"\
		"vmovaps	0x120(%%rbx),%%ymm7						\n\t		vaddpd			%%ymm10		,%%ymm14,%%ymm14			\n\t"\
		"vmovaps	0x100(%%rbx),%%ymm6						\n\t		vmovaps		%%ymm15		,%%ymm12					\n\t"\
		"vmovaps	%%ymm3,0x0a0(%%rdx)						\n\t		vmovaps		%%ymm10		,%%ymm13					\n\t"\
		"vmovaps	%%ymm2,0x080(%%rdx)		\n\t/* r04,05 */			vmulpd		 0x180(%%rcx)	,%%ymm15,%%ymm15	/* c0C */\n\t"\
		"vmovaps	%%ymm7,0x2a0(%%rdx)						\n\t		vmulpd		 0x180(%%rcx)	,%%ymm10,%%ymm10			\n\t"\
		"vmovaps	%%ymm6,0x280(%%rdx)		\n\t/* r14,15 */			vmulpd		 0x1a0(%%rcx)	,%%ymm12,%%ymm12			\n\t"\
		"vsubpd			%%ymm5		,%%ymm0,%%ymm0			\n\t		vmulpd		 0x1a0(%%rcx)	,%%ymm13,%%ymm13			\n\t"\
		"vaddpd			%%ymm4		,%%ymm1,%%ymm1			\n\t		vsubpd			%%ymm12		,%%ymm10,%%ymm10			\n\t"\
		"vmovaps		%%ymm0		,%%ymm6					\n\t		vaddpd			%%ymm13		,%%ymm15,%%ymm15			\n\t"\
		"vmovaps		%%ymm1		,%%ymm7					\n\t		vmovaps	0x1a0(%%rbx),%%ymm13						\n\t"\
		"vmulpd		 0x0c0(%%rcx)	,%%ymm0,%%ymm0	/* c18 */	\n\t	vmovaps	0x180(%%rbx),%%ymm12						\n\t"\
		"vmulpd		 0x0c0(%%rcx)	,%%ymm1,%%ymm1			\n\t		vmovaps	%%ymm10,0x1a0(%%rdx)						\n\t"\
		"vmulpd		 0x0e0(%%rcx)	,%%ymm6,%%ymm6			\n\t		vmovaps	%%ymm15,0x180(%%rdx)			\n\t"/* r0c,0d */\
		"vmulpd		 0x0e0(%%rcx)	,%%ymm7,%%ymm7			\n\t		vmovaps	%%ymm13,0x3a0(%%rdx)						\n\t"\
		"vsubpd			%%ymm6		,%%ymm1,%%ymm1			\n\t		vmovaps	%%ymm12,0x380(%%rdx)			\n\t"/* r1c,1d */\
		"vaddpd			%%ymm7		,%%ymm0,%%ymm0			\n\t		vmovaps		%%ymm8 		,%%ymm12					\n\t"\
		"vmovaps	0x320(%%rbx),%%ymm7						\n\t		vmovaps		%%ymm14		,%%ymm13					\n\t"\
		"vmovaps	0x300(%%rbx),%%ymm6						\n\t		vmulpd		 0x1c0(%%rcx)	,%%ymm8 ,%%ymm8 	/* c1C */\n\t"\
		"vmovaps	%%ymm1,0x0e0(%%rdx)						\n\t		vmulpd		 0x1c0(%%rcx)	,%%ymm14,%%ymm14			\n\t"\
		"vmovaps	%%ymm0,0x0c0(%%rdx)		\n\t/* r06,07 */			vmulpd		 0x1e0(%%rcx)	,%%ymm12,%%ymm12			\n\t"\
		"vmovaps	%%ymm7,0x2e0(%%rdx)						\n\t		vmulpd		 0x1e0(%%rcx)	,%%ymm13,%%ymm13			\n\t"\
		"vmovaps	%%ymm6,0x2c0(%%rdx)		\n\t/* r16,17 */			vsubpd			%%ymm12		,%%ymm14,%%ymm14			\n\t"\
		"																vaddpd			%%ymm13		,%%ymm8 ,%%ymm8 			\n\t"\
		"																vmovaps	0x3a0(%%rbx),%%ymm13						\n\t"\
		"																vmovaps	0x380(%%rbx),%%ymm12						\n\t"\
		"																vmovaps	%%ymm14,0x1e0(%%rdx)						\n\t"\
		"																vmovaps	%%ymm8 ,0x1c0(%%rdx)			\n\t"/* r0e,0f */\
		"																vmovaps	%%ymm13,0x3e0(%%rdx)						\n\t"\
		"																vmovaps	%%ymm12,0x3c0(%%rdx)			\n\t"/* r1e,1f */\
		"\n\t"\
		"/************************************************/	\n\t		/************************************************/	\n\t"\
		"/*...Block 3: t04,t14,t24,t34 -> r20,24,22,26:  */	\n\t		/*...Block 7: t0C,t1C,t2C,t3C -> r28,2C,2A,2E:  */	\n\t"\
		"/************************************************/	\n\t		/************************************************/	\n\t"\
		"movq		%[__r20],%%rdx								/* base-addr in rcol = r28, so rdx offset +0x100 vs lcol */	\n\t"\
		"leaq	0x020(%%rsi),%%rcx	/* cc0; Note cc0/ss0 are shared between lcol/rcol, so no rcx-offset until get to twiddles*/\n\t"\
		"vmovaps		 0x040(%%rdx)	,%%ymm4				\n\t		vmovaps		 0x140(%%rdx)	,%%ymm12				\n\t"\
		"vmovaps		 0x0c0(%%rdx)	,%%ymm0				\n\t		vmovaps		 0x1c0(%%rdx)	,%%ymm8 				\n\t"\
		"vmovaps		 0x060(%%rdx)	,%%ymm5				\n\t		vmovaps		 0x160(%%rdx)	,%%ymm13				\n\t"\
		"vmovaps		 0x0e0(%%rdx)	,%%ymm1				\n\t		vmovaps		 0x1e0(%%rdx)	,%%ymm9 				\n\t"\
		"vmovaps		 0x040(%%rdx)	,%%ymm6				\n\t		vmovaps		 0x140(%%rdx)	,%%ymm14				\n\t"\
		"vmovaps		 0x0c0(%%rdx)	,%%ymm2				\n\t		vmovaps		 0x1c0(%%rdx)	,%%ymm10				\n\t"\
		"vmovaps		 0x060(%%rdx)	,%%ymm7				\n\t		vmovaps		 0x160(%%rdx)	,%%ymm15				\n\t"\
		"vmovaps		 0x0e0(%%rdx)	,%%ymm3				\n\t		vmovaps		 0x1e0(%%rdx)	,%%ymm11				\n\t"\
		"vmulpd		      (%%rcx)	,%%ymm4,%%ymm4			\n\t		vmulpd		 0x020(%%rcx)	,%%ymm12,%%ymm12			\n\t"\
		"vmulpd		 0x020(%%rcx)	,%%ymm0,%%ymm0			\n\t		vmulpd		      (%%rcx)	,%%ymm8 ,%%ymm8 			\n\t"\
		"vmulpd		      (%%rcx)	,%%ymm5,%%ymm5			\n\t		vmulpd		 0x020(%%rcx)	,%%ymm13,%%ymm13			\n\t"\
		"vmulpd		 0x020(%%rcx)	,%%ymm1,%%ymm1			\n\t		vmulpd		      (%%rcx)	,%%ymm9 ,%%ymm9 			\n\t"\
		"vmulpd		 0x020(%%rcx)	,%%ymm6,%%ymm6			\n\t		vmulpd		      (%%rcx)	,%%ymm14,%%ymm14			\n\t"\
		"vmulpd		      (%%rcx)	,%%ymm2,%%ymm2			\n\t		vmulpd		 0x020(%%rcx)	,%%ymm10,%%ymm10			\n\t"\
		"vmulpd		 0x020(%%rcx)	,%%ymm7,%%ymm7			\n\t		vmulpd		      (%%rcx)	,%%ymm15,%%ymm15			\n\t"\
		"vmulpd		      (%%rcx)	,%%ymm3,%%ymm3			\n\t		vmulpd		 0x020(%%rcx)	,%%ymm11,%%ymm11			\n\t"\
		"vsubpd			%%ymm6		,%%ymm5,%%ymm5			\n\t		vsubpd			%%ymm14		,%%ymm13,%%ymm13			\n\t"\
		"vsubpd			%%ymm2		,%%ymm1,%%ymm1			\n\t		vsubpd			%%ymm10		,%%ymm9 ,%%ymm9 			\n\t"\
		"vaddpd			%%ymm7		,%%ymm4,%%ymm4			\n\t		vaddpd			%%ymm15		,%%ymm12,%%ymm12			\n\t"\
		"vaddpd			%%ymm3		,%%ymm0,%%ymm0			\n\t		vaddpd			%%ymm11		,%%ymm8 ,%%ymm8 			\n\t"\
		"vmovaps		%%ymm5		,%%ymm7					\n\t		vmovaps		%%ymm13		,%%ymm15					\n\t"\
		"vmovaps		%%ymm4		,%%ymm6					\n\t		vmovaps		%%ymm12		,%%ymm14					\n\t"\
		"vaddpd			%%ymm0		,%%ymm4,%%ymm4			\n\t		vaddpd			%%ymm8 		,%%ymm12,%%ymm12			\n\t"\
		"vaddpd			%%ymm1		,%%ymm5,%%ymm5			\n\t		vaddpd			%%ymm9 		,%%ymm13,%%ymm13			\n\t"\
		"vsubpd			%%ymm0		,%%ymm6,%%ymm6			\n\t		vsubpd			%%ymm8 		,%%ymm14,%%ymm14			\n\t"\
		"vsubpd			%%ymm1		,%%ymm7,%%ymm7			\n\t		vsubpd			%%ymm9 		,%%ymm15,%%ymm15			\n\t"\
		"vmovaps		 0x080(%%rdx)	,%%ymm2				\n\t		vmovaps		 0x180(%%rdx)	,%%ymm10				\n\t"\
		"vmovaps		 0x0a0(%%rdx)	,%%ymm3				\n\t		vmovaps		 0x1a0(%%rdx)	,%%ymm11				\n\t"\
		"vmovaps		      (%%rdx)	,%%ymm0				\n\t		vmovaps		 0x100(%%rdx)	,%%ymm8 				\n\t"\
		"vmovaps		 0x020(%%rdx)	,%%ymm1				\n\t		vmovaps		 0x120(%%rdx)	,%%ymm9 				\n\t"\
		"vaddpd		 0x0a0(%%rdx)	,%%ymm2,%%ymm2			\n\t		vsubpd		 0x1a0(%%rdx)	,%%ymm10,%%ymm10			\n\t"\
		"vsubpd		 0x080(%%rdx)	,%%ymm3,%%ymm3			\n\t		vaddpd		 0x180(%%rdx)	,%%ymm11,%%ymm11			\n\t"\
		"vmulpd		      (%%rsi)	,%%ymm2,%%ymm2			\n\t		vmulpd		      (%%rsi)	,%%ymm10,%%ymm10			\n\t"\
		"vmulpd		      (%%rsi)	,%%ymm3,%%ymm3			\n\t		vmulpd		      (%%rsi)	,%%ymm11,%%ymm11			\n\t"\
		"vsubpd			%%ymm2		,%%ymm0,%%ymm0			\n\t		vsubpd			%%ymm10		,%%ymm8 ,%%ymm8 			\n\t"\
		"vsubpd			%%ymm3		,%%ymm1,%%ymm1			\n\t		vsubpd			%%ymm11		,%%ymm9 ,%%ymm9 			\n\t"\
		"vaddpd			%%ymm2		,%%ymm2,%%ymm2			\n\t		vaddpd			%%ymm10		,%%ymm10,%%ymm10			\n\t"\
		"vaddpd			%%ymm3		,%%ymm3,%%ymm3			\n\t		vaddpd			%%ymm11		,%%ymm11,%%ymm11			\n\t"\
		"vaddpd			%%ymm0		,%%ymm2,%%ymm2			\n\t		vaddpd			%%ymm8 		,%%ymm10,%%ymm10			\n\t"\
		"vaddpd			%%ymm1		,%%ymm3,%%ymm3			\n\t		vaddpd			%%ymm9 		,%%ymm11,%%ymm11			\n\t"\
		"movq		%[__c02]		,%%rcx				/* base-twiddle addr in rcol = c06, so rcx offset +0x100 vs lcol */	\n\t"\
		"vsubpd			%%ymm4		,%%ymm2,%%ymm2			\n\t		vsubpd			%%ymm14		,%%ymm8 ,%%ymm8 			\n\t"\
		"vsubpd			%%ymm5		,%%ymm3,%%ymm3			\n\t		vsubpd			%%ymm15		,%%ymm9 ,%%ymm9 			\n\t"\
		"vaddpd			%%ymm4		,%%ymm4,%%ymm4			\n\t		vaddpd			%%ymm14		,%%ymm14,%%ymm14			\n\t"\
		"vaddpd			%%ymm5		,%%ymm5,%%ymm5			\n\t		vaddpd			%%ymm15		,%%ymm15,%%ymm15			\n\t"\
		"vaddpd			%%ymm2		,%%ymm4,%%ymm4			\n\t		vaddpd			%%ymm8 		,%%ymm14,%%ymm14			\n\t"\
		"vaddpd			%%ymm3		,%%ymm5,%%ymm5			\n\t		vaddpd			%%ymm9 		,%%ymm15,%%ymm15			\n\t"\
		"vmovaps		%%ymm2		, 0x040(%%rdx)			\n\t		vmovaps		%%ymm8 		, 0x140(%%rdx)			\n\t"\
		"vmovaps		%%ymm3		, 0x060(%%rdx)			\n\t		vmovaps		%%ymm9 		, 0x160(%%rdx)			\n\t"\
		"vmovaps		%%ymm4		,%%ymm2					\n\t		vmovaps		%%ymm14		,%%ymm8 					\n\t"\
		"vmovaps		%%ymm5		,%%ymm3					\n\t		vmovaps		%%ymm15		,%%ymm9 					\n\t"\
		"vmulpd		      (%%rcx)	,%%ymm4,%%ymm4			\n\t		vmulpd		 0x100(%%rcx)	,%%ymm14,%%ymm14			\n\t"\
		"vmulpd		      (%%rcx)	,%%ymm5,%%ymm5			\n\t		vmulpd		 0x100(%%rcx)	,%%ymm15,%%ymm15			\n\t"\
		"vmulpd		 0x020(%%rcx)	,%%ymm2,%%ymm2			\n\t		vmulpd		 0x120(%%rcx)	,%%ymm8 ,%%ymm8 			\n\t"\
		"vmulpd		 0x020(%%rcx)	,%%ymm3,%%ymm3			\n\t		vmulpd		 0x120(%%rcx)	,%%ymm9 ,%%ymm9 			\n\t"\
		"vsubpd			%%ymm2		,%%ymm5,%%ymm5			\n\t		vsubpd			%%ymm8 		,%%ymm15,%%ymm15			\n\t"\
		"vaddpd			%%ymm3		,%%ymm4,%%ymm4			\n\t		vaddpd			%%ymm9 		,%%ymm14,%%ymm14			\n\t"\
		"vmovaps	0x060(%%rbx),%%ymm3						\n\t		vmovaps	0x0e0(%%rbx),%%ymm9 						\n\t"\
		"vmovaps	0x040(%%rbx),%%ymm2						\n\t		vmovaps	0x0c0(%%rbx),%%ymm8 						\n\t"\
		"vmovaps	%%ymm5,0x020(%%rdx)						\n\t		vmovaps	%%ymm15,0x120(%%rdx)						\n\t"\
		"vmovaps	%%ymm4,     (%%rdx)			\n\t/* r20,21 */		vmovaps	%%ymm14,0x100(%%rdx)			\n\t"/* r28,29 */\
		"vmovaps	%%ymm3,0x220(%%rdx)						\n\t		vmovaps	%%ymm9 ,0x320(%%rdx)						\n\t"\
		"vmovaps	%%ymm2,0x200(%%rdx)			\n\t/* r30,31 */		vmovaps	%%ymm8 ,0x300(%%rdx)			\n\t"/* r38,39 */\
		"movq		%[__c12]		,%%rcx					\n\t"		/* rcol uses c16 */\
		"vmovaps		 0x040(%%rdx)	,%%ymm4				\n\t		vmovaps		 0x140(%%rdx)	,%%ymm14				\n\t"\
		"vmovaps		 0x060(%%rdx)	,%%ymm5				\n\t		vmovaps		 0x160(%%rdx)	,%%ymm15				\n\t"\
		"vmovaps		%%ymm4		,%%ymm2					\n\t		vmovaps		%%ymm14		,%%ymm8 					\n\t"\
		"vmovaps		%%ymm5		,%%ymm3					\n\t		vmovaps		%%ymm15		,%%ymm9 					\n\t"\
		"vmulpd		      (%%rcx)	,%%ymm4,%%ymm4			\n\t		vmulpd		 0x100(%%rcx)	,%%ymm14,%%ymm14			\n\t"\
		"vmulpd		      (%%rcx)	,%%ymm5,%%ymm5			\n\t		vmulpd		 0x100(%%rcx)	,%%ymm15,%%ymm15			\n\t"\
		"vmulpd		 0x020(%%rcx)	,%%ymm2,%%ymm2			\n\t		vmulpd		 0x120(%%rcx)	,%%ymm8 ,%%ymm8 			\n\t"\
		"vmulpd		 0x020(%%rcx)	,%%ymm3,%%ymm3			\n\t		vmulpd		 0x120(%%rcx)	,%%ymm9 ,%%ymm9 			\n\t"\
		"vsubpd			%%ymm2		,%%ymm5,%%ymm5			\n\t		vsubpd			%%ymm8 		,%%ymm15,%%ymm15			\n\t"\
		"vaddpd			%%ymm3		,%%ymm4,%%ymm4			\n\t		vaddpd			%%ymm9 		,%%ymm14,%%ymm14			\n\t"\
		"vmovaps	0x260(%%rbx),%%ymm3						\n\t		vmovaps	0x2e0(%%rbx),%%ymm9 						\n\t"\
		"vmovaps	0x240(%%rbx),%%ymm2						\n\t		vmovaps	0x2c0(%%rbx),%%ymm8 						\n\t"\
		"vmovaps	%%ymm5,0x060(%%rdx)						\n\t		vmovaps	%%ymm15,0x160(%%rdx)						\n\t"\
		"vmovaps	%%ymm4,0x040(%%rdx)			\n\t/* r22,23 */		vmovaps	%%ymm14,0x140(%%rdx)			\n\t"/* r2a,2b */\
		"vmovaps	%%ymm3,0x260(%%rdx)						\n\t		vmovaps	%%ymm9 ,0x360(%%rdx)						\n\t"\
		"vmovaps	%%ymm2,0x240(%%rdx)			\n\t/* r32,33 */		vmovaps	%%ymm8 ,0x340(%%rdx)			\n\t"/* r3a,3b */\
		"movq		%[__c0A]		,%%rcx					\n\t"		/* rcol uses c0E */\
		"vsubpd			%%ymm7		,%%ymm0,%%ymm0			\n\t		vsubpd			%%ymm13		,%%ymm10,%%ymm10			\n\t"\
		"vsubpd			%%ymm6		,%%ymm1,%%ymm1			\n\t		vsubpd			%%ymm12		,%%ymm11,%%ymm11			\n\t"\
		"vaddpd			%%ymm7		,%%ymm7,%%ymm7			\n\t		vaddpd			%%ymm13		,%%ymm13,%%ymm13			\n\t"\
		"vaddpd			%%ymm6		,%%ymm6,%%ymm6			\n\t		vaddpd			%%ymm12		,%%ymm12,%%ymm12			\n\t"\
		"vaddpd			%%ymm0		,%%ymm7,%%ymm7			\n\t		vaddpd			%%ymm10		,%%ymm13,%%ymm13			\n\t"\
		"vaddpd			%%ymm1		,%%ymm6,%%ymm6			\n\t		vaddpd			%%ymm11		,%%ymm12,%%ymm12			\n\t"\
		"vmovaps		%%ymm7		,%%ymm4					\n\t		vmovaps		%%ymm13		,%%ymm8 					\n\t"\
		"vmovaps		%%ymm1		,%%ymm5					\n\t		vmovaps		%%ymm11		,%%ymm9 					\n\t"\
		"vmulpd		      (%%rcx)	,%%ymm7,%%ymm7			\n\t		vmulpd		 0x100(%%rcx)	,%%ymm13,%%ymm13			\n\t"\
		"vmulpd		      (%%rcx)	,%%ymm1,%%ymm1			\n\t		vmulpd		 0x100(%%rcx)	,%%ymm11,%%ymm11			\n\t"\
		"vmulpd		 0x020(%%rcx)	,%%ymm4,%%ymm4			\n\t		vmulpd		 0x120(%%rcx)	,%%ymm8 ,%%ymm8 			\n\t"\
		"vmulpd		 0x020(%%rcx)	,%%ymm5,%%ymm5			\n\t		vmulpd		 0x120(%%rcx)	,%%ymm9 ,%%ymm9 			\n\t"\
		"vsubpd			%%ymm4		,%%ymm1,%%ymm1			\n\t		vsubpd			%%ymm8 		,%%ymm11,%%ymm11			\n\t"\
		"vaddpd			%%ymm5		,%%ymm7,%%ymm7			\n\t		vaddpd			%%ymm9 		,%%ymm13,%%ymm13			\n\t"\
		"vmovaps	0x160(%%rbx),%%ymm5						\n\t		vmovaps	0x1e0(%%rbx),%%ymm9 						\n\t"\
		"vmovaps	0x140(%%rbx),%%ymm4						\n\t		vmovaps	0x1c0(%%rbx),%%ymm8 						\n\t"\
		"vmovaps	%%ymm1,0x0a0(%%rdx)						\n\t		vmovaps	%%ymm11,0x1a0(%%rdx)						\n\t"\
		"vmovaps	%%ymm7,0x080(%%rdx)			\n\t/* r24,25 */		vmovaps	%%ymm13,0x180(%%rdx)			\n\t"/* r2c,2d */\
		"vmovaps	%%ymm5,0x2a0(%%rdx)						\n\t		vmovaps	%%ymm9 ,0x3a0(%%rdx)						\n\t"\
		"vmovaps	%%ymm4,0x280(%%rdx)			\n\t/* r34,35 */		vmovaps	%%ymm8 ,0x380(%%rdx)			\n\t"/* r3c,3d */\
		"movq		%[__c1A]		,%%rcx					\n\t"		/* rcol uses c1E */\
		"vmovaps		%%ymm0		,%%ymm4					\n\t		vmovaps		%%ymm10		,%%ymm8 					\n\t"\
		"vmovaps		%%ymm6		,%%ymm5					\n\t		vmovaps		%%ymm12		,%%ymm9 					\n\t"\
		"vmulpd		      (%%rcx)	,%%ymm0,%%ymm0			\n\t		vmulpd		 0x100(%%rcx)	,%%ymm10,%%ymm10			\n\t"\
		"vmulpd		      (%%rcx)	,%%ymm6,%%ymm6			\n\t		vmulpd		 0x100(%%rcx)	,%%ymm12,%%ymm12			\n\t"\
		"vmulpd		 0x020(%%rcx)	,%%ymm4,%%ymm4			\n\t		vmulpd		 0x120(%%rcx)	,%%ymm8 ,%%ymm8 			\n\t"\
		"vmulpd		 0x020(%%rcx)	,%%ymm5,%%ymm5			\n\t		vmulpd		 0x120(%%rcx)	,%%ymm9 ,%%ymm9 			\n\t"\
		"vsubpd			%%ymm4		,%%ymm6,%%ymm6			\n\t		vsubpd			%%ymm8 		,%%ymm12,%%ymm12			\n\t"\
		"vaddpd			%%ymm5		,%%ymm0,%%ymm0			\n\t		vaddpd			%%ymm9 		,%%ymm10,%%ymm10			\n\t"\
		"vmovaps	0x360(%%rbx),%%ymm5						\n\t		vmovaps	0x3e0(%%rbx),%%ymm9 						\n\t"\
		"vmovaps	0x340(%%rbx),%%ymm4						\n\t		vmovaps	0x3c0(%%rbx),%%ymm8 						\n\t"\
		"vmovaps	%%ymm6,0x0e0(%%rdx)						\n\t		vmovaps	%%ymm12,0x1e0(%%rdx)						\n\t"\
		"vmovaps	%%ymm0,0x0c0(%%rdx)			\n\t/* r26,27 */		vmovaps	%%ymm10,0x1c0(%%rdx)			\n\t"/* r2e,2f */\
		"vmovaps	%%ymm5,0x2e0(%%rdx)						\n\t		vmovaps	%%ymm9 ,0x3e0(%%rdx)						\n\t"\
		"vmovaps	%%ymm4,0x2c0(%%rdx)			\n\t/* r36,37 */		vmovaps	%%ymm8 ,0x3c0(%%rdx)			\n\t"/* r3e,3f */\
/*==========================*/"\n\t"\
	"/**** Finish with 4-way 'un'terleaving: ****/\n\t"\
		"movq	%[__r00] ,%%rsi\n\t"\
		"movq	%[__add0],%%rax\n\t"\
		"movq	%[__add1],%%rbx\n\t"\
		"movq	%[__add2],%%rcx\n\t"\
		"movq	%[__add3],%%rdx\n\t"\
	"/* a[j+p0]: Inputs from r00 +0/1, 8/9, 16/17, 24/25. Outputs into [add0,add0+0x100,add1,add1+0x100]+0x0: */	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps	     (%%rsi),%%ymm1						\n\t		vmovaps	0x020(%%rsi),%%ymm3								\n\t"\
		"vmovaps	0x200(%%rsi),%%ymm5						\n\t		vmovaps	0x220(%%rsi),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps	0x400(%%rsi),%%ymm6						\n\t		vmovaps	0x420(%%rsi),%%ymm14							\n\t"\
		"vmovaps	0x600(%%rsi),%%ymm5						\n\t		vmovaps	0x620(%%rsi),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm6 ,%%ymm0				\n\t		vshufpd	$15 ,%%ymm13,%%ymm14 ,%%ymm2					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm6 ,%%ymm6				\n\t		vshufpd	$0  ,%%ymm13,%%ymm14 ,%%ymm14					\n\t"\
		"vperm2f128 $32 ,%%ymm0 ,%%ymm7 ,%%ymm5				\n\t		vperm2f128 $32 ,%%ymm2 ,%%ymm15,%%ymm13		/* outB	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm0 ,%%ymm7 ,%%ymm7				\n\t		vperm2f128 $49 ,%%ymm2 ,%%ymm15,%%ymm15		/* outD	*/	\n\t"\
		"vperm2f128 $32 ,%%ymm6 ,%%ymm1 ,%%ymm0				\n\t		vperm2f128 $32 ,%%ymm14,%%ymm3 ,%%ymm2 		/* outA	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm6 ,%%ymm1 ,%%ymm1				\n\t		vperm2f128 $49 ,%%ymm14,%%ymm3 ,%%ymm3 		/* outC	*/	\n\t"\
		"vmovaps 	%%ymm5 ,(%%rbx)							\n\t		vmovaps %%ymm13,0x020(%%rbx)				/* outB	*/	\n\t"\
		"vmovaps 	%%ymm7 ,(%%rdx)							\n\t		vmovaps %%ymm15,0x020(%%rdx)				/* outD	*/	\n\t"\
		"vmovaps 	%%ymm0 ,(%%rax)							\n\t		vmovaps %%ymm2 ,0x020(%%rax)				/* outA	*/	\n\t"\
		"vmovaps 	%%ymm1 ,(%%rcx)							\n\t		vmovaps %%ymm3 ,0x020(%%rcx)				/* outC	*/	\n\t"\
	"/* a[j+p2]: Inputs from r08 +0/1, 8/9, 16/17, 24/25. Outputs into [add0,add0+0x100,add1,add1+0x100]+0x0: */	\n\t"\
		"addq	$0x100,%%rsi	\n\t"\
		"addq	$0x40,%%rax	\n\t"\
		"addq	$0x40,%%rbx	\n\t"\
		"addq	$0x40,%%rcx	\n\t"\
		"addq	$0x40,%%rdx	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps	     (%%rsi),%%ymm1						\n\t		vmovaps	0x020(%%rsi),%%ymm3								\n\t"\
		"vmovaps	0x200(%%rsi),%%ymm5						\n\t		vmovaps	0x220(%%rsi),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps	0x400(%%rsi),%%ymm6						\n\t		vmovaps	0x420(%%rsi),%%ymm14							\n\t"\
		"vmovaps	0x600(%%rsi),%%ymm5						\n\t		vmovaps	0x620(%%rsi),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm6 ,%%ymm0				\n\t		vshufpd	$15 ,%%ymm13,%%ymm14 ,%%ymm2					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm6 ,%%ymm6				\n\t		vshufpd	$0  ,%%ymm13,%%ymm14 ,%%ymm14					\n\t"\
		"vperm2f128 $32 ,%%ymm0 ,%%ymm7 ,%%ymm5				\n\t		vperm2f128 $32 ,%%ymm2 ,%%ymm15,%%ymm13		/* outB	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm0 ,%%ymm7 ,%%ymm7				\n\t		vperm2f128 $49 ,%%ymm2 ,%%ymm15,%%ymm15		/* outD	*/	\n\t"\
		"vperm2f128 $32 ,%%ymm6 ,%%ymm1 ,%%ymm0				\n\t		vperm2f128 $32 ,%%ymm14,%%ymm3 ,%%ymm2 		/* outA	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm6 ,%%ymm1 ,%%ymm1				\n\t		vperm2f128 $49 ,%%ymm14,%%ymm3 ,%%ymm3 		/* outC	*/	\n\t"\
		"vmovaps 	%%ymm5 ,(%%rbx)							\n\t		vmovaps %%ymm13,0x020(%%rbx)				/* outB	*/	\n\t"\
		"vmovaps 	%%ymm7 ,(%%rdx)							\n\t		vmovaps %%ymm15,0x020(%%rdx)				/* outD	*/	\n\t"\
		"vmovaps 	%%ymm0 ,(%%rax)							\n\t		vmovaps %%ymm2 ,0x020(%%rax)				/* outA	*/	\n\t"\
		"vmovaps 	%%ymm1 ,(%%rcx)							\n\t		vmovaps %%ymm3 ,0x020(%%rcx)				/* outC	*/	\n\t"\
	"/* a[j+p4]: Inputs from r04 +0/1, 8/9, 16/17, 24/25. Outputs into [add0,add0+0x100,add1,add1+0x100]+0x0: */	\n\t"\
		"subq	$0x80,%%rsi	\n\t"\
		"addq	$0x40,%%rax	\n\t"\
		"addq	$0x40,%%rbx	\n\t"\
		"addq	$0x40,%%rcx	\n\t"\
		"addq	$0x40,%%rdx	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps	     (%%rsi),%%ymm1						\n\t		vmovaps	0x020(%%rsi),%%ymm3								\n\t"\
		"vmovaps	0x200(%%rsi),%%ymm5						\n\t		vmovaps	0x220(%%rsi),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps	0x400(%%rsi),%%ymm6						\n\t		vmovaps	0x420(%%rsi),%%ymm14							\n\t"\
		"vmovaps	0x600(%%rsi),%%ymm5						\n\t		vmovaps	0x620(%%rsi),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm6 ,%%ymm0				\n\t		vshufpd	$15 ,%%ymm13,%%ymm14 ,%%ymm2					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm6 ,%%ymm6				\n\t		vshufpd	$0  ,%%ymm13,%%ymm14 ,%%ymm14					\n\t"\
		"vperm2f128 $32 ,%%ymm0 ,%%ymm7 ,%%ymm5				\n\t		vperm2f128 $32 ,%%ymm2 ,%%ymm15,%%ymm13		/* outB	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm0 ,%%ymm7 ,%%ymm7				\n\t		vperm2f128 $49 ,%%ymm2 ,%%ymm15,%%ymm15		/* outD	*/	\n\t"\
		"vperm2f128 $32 ,%%ymm6 ,%%ymm1 ,%%ymm0				\n\t		vperm2f128 $32 ,%%ymm14,%%ymm3 ,%%ymm2 		/* outA	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm6 ,%%ymm1 ,%%ymm1				\n\t		vperm2f128 $49 ,%%ymm14,%%ymm3 ,%%ymm3 		/* outC	*/	\n\t"\
		"vmovaps 	%%ymm5 ,(%%rbx)							\n\t		vmovaps %%ymm13,0x020(%%rbx)				/* outB	*/	\n\t"\
		"vmovaps 	%%ymm7 ,(%%rdx)							\n\t		vmovaps %%ymm15,0x020(%%rdx)				/* outD	*/	\n\t"\
		"vmovaps 	%%ymm0 ,(%%rax)							\n\t		vmovaps %%ymm2 ,0x020(%%rax)				/* outA	*/	\n\t"\
		"vmovaps 	%%ymm1 ,(%%rcx)							\n\t		vmovaps %%ymm3 ,0x020(%%rcx)				/* outC	*/	\n\t"\
	"/* a[j+p6]: Inputs from r0c +0/1, 8/9, 16/17, 24/25. Outputs into [add0,add0+0x100,add1,add1+0x100]+0x0: */	\n\t"\
		"addq	$0x100,%%rsi	\n\t"\
		"addq	$0x40,%%rax	\n\t"\
		"addq	$0x40,%%rbx	\n\t"\
		"addq	$0x40,%%rcx	\n\t"\
		"addq	$0x40,%%rdx	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps	     (%%rsi),%%ymm1						\n\t		vmovaps	0x020(%%rsi),%%ymm3								\n\t"\
		"vmovaps	0x200(%%rsi),%%ymm5						\n\t		vmovaps	0x220(%%rsi),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps	0x400(%%rsi),%%ymm6						\n\t		vmovaps	0x420(%%rsi),%%ymm14							\n\t"\
		"vmovaps	0x600(%%rsi),%%ymm5						\n\t		vmovaps	0x620(%%rsi),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm6 ,%%ymm0				\n\t		vshufpd	$15 ,%%ymm13,%%ymm14 ,%%ymm2					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm6 ,%%ymm6				\n\t		vshufpd	$0  ,%%ymm13,%%ymm14 ,%%ymm14					\n\t"\
		"vperm2f128 $32 ,%%ymm0 ,%%ymm7 ,%%ymm5				\n\t		vperm2f128 $32 ,%%ymm2 ,%%ymm15,%%ymm13		/* outB	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm0 ,%%ymm7 ,%%ymm7				\n\t		vperm2f128 $49 ,%%ymm2 ,%%ymm15,%%ymm15		/* outD	*/	\n\t"\
		"vperm2f128 $32 ,%%ymm6 ,%%ymm1 ,%%ymm0				\n\t		vperm2f128 $32 ,%%ymm14,%%ymm3 ,%%ymm2 		/* outA	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm6 ,%%ymm1 ,%%ymm1				\n\t		vperm2f128 $49 ,%%ymm14,%%ymm3 ,%%ymm3 		/* outC	*/	\n\t"\
		"vmovaps 	%%ymm5 ,(%%rbx)							\n\t		vmovaps %%ymm13,0x020(%%rbx)				/* outB	*/	\n\t"\
		"vmovaps 	%%ymm7 ,(%%rdx)							\n\t		vmovaps %%ymm15,0x020(%%rdx)				/* outD	*/	\n\t"\
		"vmovaps 	%%ymm0 ,(%%rax)							\n\t		vmovaps %%ymm2 ,0x020(%%rax)				/* outA	*/	\n\t"\
		"vmovaps 	%%ymm1 ,(%%rcx)							\n\t		vmovaps %%ymm3 ,0x020(%%rcx)				/* outC	*/	\n\t"\
	"/* a[j+p8]: Inputs from r02 +0/1, 8/9, 16/17, 24/25. Outputs into [add0,add0+0x100,add1,add1+0x100]+0x0: */	\n\t"\
		"subq	$0x140,%%rsi	\n\t"\
		"addq	$0x40,%%rax	\n\t"\
		"addq	$0x40,%%rbx	\n\t"\
		"addq	$0x40,%%rcx	\n\t"\
		"addq	$0x40,%%rdx	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps	     (%%rsi),%%ymm1						\n\t		vmovaps	0x020(%%rsi),%%ymm3								\n\t"\
		"vmovaps	0x200(%%rsi),%%ymm5						\n\t		vmovaps	0x220(%%rsi),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps	0x400(%%rsi),%%ymm6						\n\t		vmovaps	0x420(%%rsi),%%ymm14							\n\t"\
		"vmovaps	0x600(%%rsi),%%ymm5						\n\t		vmovaps	0x620(%%rsi),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm6 ,%%ymm0				\n\t		vshufpd	$15 ,%%ymm13,%%ymm14 ,%%ymm2					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm6 ,%%ymm6				\n\t		vshufpd	$0  ,%%ymm13,%%ymm14 ,%%ymm14					\n\t"\
		"vperm2f128 $32 ,%%ymm0 ,%%ymm7 ,%%ymm5				\n\t		vperm2f128 $32 ,%%ymm2 ,%%ymm15,%%ymm13		/* outB	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm0 ,%%ymm7 ,%%ymm7				\n\t		vperm2f128 $49 ,%%ymm2 ,%%ymm15,%%ymm15		/* outD	*/	\n\t"\
		"vperm2f128 $32 ,%%ymm6 ,%%ymm1 ,%%ymm0				\n\t		vperm2f128 $32 ,%%ymm14,%%ymm3 ,%%ymm2 		/* outA	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm6 ,%%ymm1 ,%%ymm1				\n\t		vperm2f128 $49 ,%%ymm14,%%ymm3 ,%%ymm3 		/* outC	*/	\n\t"\
		"vmovaps 	%%ymm5 ,(%%rbx)							\n\t		vmovaps %%ymm13,0x020(%%rbx)				/* outB	*/	\n\t"\
		"vmovaps 	%%ymm7 ,(%%rdx)							\n\t		vmovaps %%ymm15,0x020(%%rdx)				/* outD	*/	\n\t"\
		"vmovaps 	%%ymm0 ,(%%rax)							\n\t		vmovaps %%ymm2 ,0x020(%%rax)				/* outA	*/	\n\t"\
		"vmovaps 	%%ymm1 ,(%%rcx)							\n\t		vmovaps %%ymm3 ,0x020(%%rcx)				/* outC	*/	\n\t"\
	"/* a[j+p10]: Inputs from r0a +0/1, 8/9, 16/17, 24/25. Outputs into [add0,add0+0x100,add1,add1+0x100]+0x0: */	\n\t"\
		"addq	$0x100,%%rsi	\n\t"\
		"addq	$0x40,%%rax	\n\t"\
		"addq	$0x40,%%rbx	\n\t"\
		"addq	$0x40,%%rcx	\n\t"\
		"addq	$0x40,%%rdx	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps	     (%%rsi),%%ymm1						\n\t		vmovaps	0x020(%%rsi),%%ymm3								\n\t"\
		"vmovaps	0x200(%%rsi),%%ymm5						\n\t		vmovaps	0x220(%%rsi),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps	0x400(%%rsi),%%ymm6						\n\t		vmovaps	0x420(%%rsi),%%ymm14							\n\t"\
		"vmovaps	0x600(%%rsi),%%ymm5						\n\t		vmovaps	0x620(%%rsi),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm6 ,%%ymm0				\n\t		vshufpd	$15 ,%%ymm13,%%ymm14 ,%%ymm2					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm6 ,%%ymm6				\n\t		vshufpd	$0  ,%%ymm13,%%ymm14 ,%%ymm14					\n\t"\
		"vperm2f128 $32 ,%%ymm0 ,%%ymm7 ,%%ymm5				\n\t		vperm2f128 $32 ,%%ymm2 ,%%ymm15,%%ymm13		/* outB	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm0 ,%%ymm7 ,%%ymm7				\n\t		vperm2f128 $49 ,%%ymm2 ,%%ymm15,%%ymm15		/* outD	*/	\n\t"\
		"vperm2f128 $32 ,%%ymm6 ,%%ymm1 ,%%ymm0				\n\t		vperm2f128 $32 ,%%ymm14,%%ymm3 ,%%ymm2 		/* outA	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm6 ,%%ymm1 ,%%ymm1				\n\t		vperm2f128 $49 ,%%ymm14,%%ymm3 ,%%ymm3 		/* outC	*/	\n\t"\
		"vmovaps 	%%ymm5 ,(%%rbx)							\n\t		vmovaps %%ymm13,0x020(%%rbx)				/* outB	*/	\n\t"\
		"vmovaps 	%%ymm7 ,(%%rdx)							\n\t		vmovaps %%ymm15,0x020(%%rdx)				/* outD	*/	\n\t"\
		"vmovaps 	%%ymm0 ,(%%rax)							\n\t		vmovaps %%ymm2 ,0x020(%%rax)				/* outA	*/	\n\t"\
		"vmovaps 	%%ymm1 ,(%%rcx)							\n\t		vmovaps %%ymm3 ,0x020(%%rcx)				/* outC	*/	\n\t"\
	"/* a[j+p12]: Inputs from r06 +0/1, 8/9, 16/17, 24/25. Outputs into [add0,add0+0x100,add1,add1+0x100]+0x0: */	\n\t"\
		"subq	$0x80,%%rsi	\n\t"\
		"addq	$0x40,%%rax	\n\t"\
		"addq	$0x40,%%rbx	\n\t"\
		"addq	$0x40,%%rcx	\n\t"\
		"addq	$0x40,%%rdx	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps	     (%%rsi),%%ymm1						\n\t		vmovaps	0x020(%%rsi),%%ymm3								\n\t"\
		"vmovaps	0x200(%%rsi),%%ymm5						\n\t		vmovaps	0x220(%%rsi),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps	0x400(%%rsi),%%ymm6						\n\t		vmovaps	0x420(%%rsi),%%ymm14							\n\t"\
		"vmovaps	0x600(%%rsi),%%ymm5						\n\t		vmovaps	0x620(%%rsi),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm6 ,%%ymm0				\n\t		vshufpd	$15 ,%%ymm13,%%ymm14 ,%%ymm2					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm6 ,%%ymm6				\n\t		vshufpd	$0  ,%%ymm13,%%ymm14 ,%%ymm14					\n\t"\
		"vperm2f128 $32 ,%%ymm0 ,%%ymm7 ,%%ymm5				\n\t		vperm2f128 $32 ,%%ymm2 ,%%ymm15,%%ymm13		/* outB	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm0 ,%%ymm7 ,%%ymm7				\n\t		vperm2f128 $49 ,%%ymm2 ,%%ymm15,%%ymm15		/* outD	*/	\n\t"\
		"vperm2f128 $32 ,%%ymm6 ,%%ymm1 ,%%ymm0				\n\t		vperm2f128 $32 ,%%ymm14,%%ymm3 ,%%ymm2 		/* outA	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm6 ,%%ymm1 ,%%ymm1				\n\t		vperm2f128 $49 ,%%ymm14,%%ymm3 ,%%ymm3 		/* outC	*/	\n\t"\
		"vmovaps 	%%ymm5 ,(%%rbx)							\n\t		vmovaps %%ymm13,0x020(%%rbx)				/* outB	*/	\n\t"\
		"vmovaps 	%%ymm7 ,(%%rdx)							\n\t		vmovaps %%ymm15,0x020(%%rdx)				/* outD	*/	\n\t"\
		"vmovaps 	%%ymm0 ,(%%rax)							\n\t		vmovaps %%ymm2 ,0x020(%%rax)				/* outA	*/	\n\t"\
		"vmovaps 	%%ymm1 ,(%%rcx)							\n\t		vmovaps %%ymm3 ,0x020(%%rcx)				/* outC	*/	\n\t"\
	"/* a[j+p14]: Inputs from r0e +0/1, 8/9, 16/17, 24/25. Outputs into [add0,add0+0x100,add1,add1+0x100]+0x0: */	\n\t"\
		"addq	$0x100,%%rsi	\n\t"\
		"addq	$0x40,%%rax	\n\t"\
		"addq	$0x40,%%rbx	\n\t"\
		"addq	$0x40,%%rcx	\n\t"\
		"addq	$0x40,%%rdx	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps	     (%%rsi),%%ymm1						\n\t		vmovaps	0x020(%%rsi),%%ymm3								\n\t"\
		"vmovaps	0x200(%%rsi),%%ymm5						\n\t		vmovaps	0x220(%%rsi),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps	0x400(%%rsi),%%ymm6						\n\t		vmovaps	0x420(%%rsi),%%ymm14							\n\t"\
		"vmovaps	0x600(%%rsi),%%ymm5						\n\t		vmovaps	0x620(%%rsi),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm6 ,%%ymm0				\n\t		vshufpd	$15 ,%%ymm13,%%ymm14 ,%%ymm2					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm6 ,%%ymm6				\n\t		vshufpd	$0  ,%%ymm13,%%ymm14 ,%%ymm14					\n\t"\
		"vperm2f128 $32 ,%%ymm0 ,%%ymm7 ,%%ymm5				\n\t		vperm2f128 $32 ,%%ymm2 ,%%ymm15,%%ymm13		/* outB	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm0 ,%%ymm7 ,%%ymm7				\n\t		vperm2f128 $49 ,%%ymm2 ,%%ymm15,%%ymm15		/* outD	*/	\n\t"\
		"vperm2f128 $32 ,%%ymm6 ,%%ymm1 ,%%ymm0				\n\t		vperm2f128 $32 ,%%ymm14,%%ymm3 ,%%ymm2 		/* outA	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm6 ,%%ymm1 ,%%ymm1				\n\t		vperm2f128 $49 ,%%ymm14,%%ymm3 ,%%ymm3 		/* outC	*/	\n\t"\
		"vmovaps 	%%ymm5 ,(%%rbx)							\n\t		vmovaps %%ymm13,0x020(%%rbx)				/* outB	*/	\n\t"\
		"vmovaps 	%%ymm7 ,(%%rdx)							\n\t		vmovaps %%ymm15,0x020(%%rdx)				/* outD	*/	\n\t"\
		"vmovaps 	%%ymm0 ,(%%rax)							\n\t		vmovaps %%ymm2 ,0x020(%%rax)				/* outA	*/	\n\t"\
		"vmovaps 	%%ymm1 ,(%%rcx)							\n\t		vmovaps %%ymm3 ,0x020(%%rcx)				/* outC	*/	\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__add1] "m" (Xadd1)\
		 ,[__add2] "m" (Xadd2)\
		 ,[__add3] "m" (Xadd3)\
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__r00] "m" (Xr00)\
		 ,[__r08] "m" (Xr08)\
		 ,[__r10] "m" (Xr10)\
		 ,[__r20] "m" (Xr20)\
		 ,[__r28] "m" (Xr28)\
		 ,[__r30] "m" (Xr30)\
		 ,[__c01] "m" (Xc01)\
		 ,[__c02] "m" (Xc02)\
		 ,[__c04] "m" (Xc04)\
		 ,[__c06] "m" (Xc06)\
		 ,[__c08] "m" (Xc08)\
		 ,[__c0A] "m" (Xc0A)\
		 ,[__c0C] "m" (Xc0C)\
		 ,[__c0E] "m" (Xc0E)\
		 ,[__c10] "m" (Xc10)\
		 ,[__c12] "m" (Xc12)\
		 ,[__c14] "m" (Xc14)\
		 ,[__c16] "m" (Xc16)\
		 ,[__c18] "m" (Xc18)\
		 ,[__c1A] "m" (Xc1A)\
		 ,[__c1C] "m" (Xc1C)\
		 ,[__c1E] "m" (Xc1E)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"		/* Clobbered registers */\
	);\
	}

	#define BAR(Xadd0,Xadd1,Xisrt2,Xr00,Xr08,Xr10,Xr20,Xr28,Xr30,Xc00,Xc01,Xc02,Xc03,Xc04,Xc05,Xc06,Xc07,Xc08,Xc0A,Xc0C,Xc0E,Xc10,Xc12,Xc14,Xc16,Xc18,Xc1A,Xc1C,Xc1E)\
	{\
	__asm__ volatile (\
		"\n\t"\
		"movq	%[__isrt2],%%rbx	\n\t	addq $0x920,%%rbx	/* DEBUG: Dump register contents into scratch-data slots: */\n\t"\
		"vmovaps	%%ymm0 ,0x000(%%rbx)	\n\t"\
		"vmovaps	%%ymm1 ,0x020(%%rbx)	\n\t"\
		"vmovaps	%%ymm2 ,0x040(%%rbx)	\n\t"\
		"vmovaps	%%ymm3 ,0x060(%%rbx)	\n\t"\
		"vmovaps	%%ymm4 ,0x080(%%rbx)	\n\t"\
		"vmovaps	%%ymm5 ,0x0a0(%%rbx)	\n\t"\
		"vmovaps	%%ymm6 ,0x0c0(%%rbx)	\n\t"\
		"vmovaps	%%ymm7 ,0x0e0(%%rbx)	\n\t"\
		"vmovaps	%%ymm8 ,0x100(%%rbx)	\n\t"\
		"vmovaps	%%ymm9 ,0x120(%%rbx)	\n\t"\
		"vmovaps	%%ymm10,0x140(%%rbx)	\n\t"\
		"vmovaps	%%ymm11,0x160(%%rbx)	\n\t"\
		"vmovaps	%%ymm12,0x180(%%rbx)	\n\t"\
		"vmovaps	%%ymm13,0x1a0(%%rbx)	\n\t"\
		"vmovaps	%%ymm14,0x1c0(%%rbx)	\n\t"\
		"vmovaps	%%ymm15,0x1e0(%%rbx)	\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__add1] "m" (Xadd1)\
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__r00] "m" (Xr00)\
		 ,[__r08] "m" (Xr08)\
		 ,[__r10] "m" (Xr10)\
		 ,[__r20] "m" (Xr20)\
		 ,[__r28] "m" (Xr28)\
		 ,[__r30] "m" (Xr30)\
		 ,[__c00] "m" (Xc00)\
		 ,[__c01] "m" (Xc01)\
		 ,[__c02] "m" (Xc02)\
		 ,[__c03] "m" (Xc03)\
		 ,[__c04] "m" (Xc04)\
		 ,[__c05] "m" (Xc05)\
		 ,[__c06] "m" (Xc06)\
		 ,[__c07] "m" (Xc07)\
		 ,[__c08] "m" (Xc08)\
		 ,[__c0A] "m" (Xc0A)\
		 ,[__c0C] "m" (Xc0C)\
		 ,[__c0E] "m" (Xc0E)\
		 ,[__c10] "m" (Xc10)\
		 ,[__c12] "m" (Xc12)\
		 ,[__c14] "m" (Xc14)\
		 ,[__c16] "m" (Xc16)\
		 ,[__c18] "m" (Xc18)\
		 ,[__c1A] "m" (Xc1A)\
		 ,[__c1C] "m" (Xc1C)\
		 ,[__c1E] "m" (Xc1E)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"		/* Clobbered registers */\
	);\
	}

#elif defined(USE_SSE2)

	#define SSE2_RADIX32_WRAPPER_DIF(Xadd0,Xadd1,Xr00,Xr10,Xr20,Xr30,Xisrt2,Xcc0,Xc00,Xc01,Xc02,Xc03,Xc05,Xc07)\
	{\
	__asm__ volatile (\
		"/************************************************************************/\n\t"\
		"/* Forward DIF radix-32 pass on the interleaved block1 and block2 data: */\n\t"\
		"/************************************************************************/\n\t"\
		"/*...Block 0:	*/									\n\t		movq	%[__isrt2],%%rsi		\n\t"\
		"movq		%[__add0]	,%%rax						\n\t		movq		%%rsi,%%rdi			\n\t"\
		"movq		%[__add1]	,%%rbx						\n\t		addq	$0x480,%%rdi	/* two */	\n\t"\
		"/*****	SSE2_RADIX4_DIF_4WRAPPER(c00,c08,c10,c18,r00) ****	SSE2_RADIX4_DIF_4WRAPPER_2NDOFTWO(c04,c0C,c14,c1C,r08)	*****/\n\t"\
		"movq		%[__r00]	,%%rcx						\n\t		/*addq		$0x080,%%rcx // __r08 */	\n\t"\
		"movq		%[__c00]	,%%rdx						\n\t		/*addq		$0x080,%%rdx // __c04 */	\n\t"\
		"movaps		     (%%rax),%%xmm6						\n\t		movaps		0x40 (%%rax),%%xmm14		\n\t"\
		"movaps		     (%%rax),%%xmm0						\n\t		movaps		0x40 (%%rax),%%xmm8			\n\t"\
		"unpckhpd	     (%%rbx),%%xmm6						\n\t		unpckhpd	0x40 (%%rbx),%%xmm14		\n\t"\
		"unpcklpd	     (%%rbx),%%xmm0						\n\t		unpcklpd	0x40 (%%rbx),%%xmm8			\n\t"\
		"movaps		%%xmm6		,0x200(%%rcx)				\n\t		movaps		%%xmm14,0x280(%%rcx)		\n\t"\
		"movaps		0x010(%%rax),%%xmm7						\n\t		movaps		0x50 (%%rax),%%xmm15		\n\t"\
		"movaps		0x010(%%rax),%%xmm1						\n\t		movaps		0x50 (%%rax),%%xmm9			\n\t"\
		"unpckhpd	0x010(%%rbx),%%xmm7						\n\t		unpckhpd	0x50 (%%rbx),%%xmm15		\n\t"\
		"unpcklpd	0x010(%%rbx),%%xmm1						\n\t		unpcklpd	0x50 (%%rbx),%%xmm9			\n\t"\
		"movaps		%%xmm7		,0x210(%%rcx)				\n\t		movaps		%%xmm15,0x290(%%rcx)		\n\t"\
		"movaps		%%xmm0		,%%xmm2						\n\t		movaps		%%xmm8,%%xmm10				\n\t"\
		"movaps		%%xmm1		,%%xmm3						\n\t		movaps		%%xmm9,%%xmm11				\n\t"\
		"mulpd		     (%%rdx),%%xmm0						\n\t		mulpd		0x080(%%rdx),%%xmm8			\n\t"\
		"mulpd		     (%%rdx),%%xmm1						\n\t		mulpd		0x080(%%rdx),%%xmm9			\n\t"\
		"mulpd		0x010(%%rdx),%%xmm2						\n\t		mulpd		0x090(%%rdx),%%xmm10		\n\t"\
		"mulpd		0x010(%%rdx),%%xmm3						\n\t		mulpd		0x090(%%rdx),%%xmm11		\n\t"\
		"addpd		%%xmm2		,%%xmm1						\n\t		addpd		%%xmm10,%%xmm9				\n\t"\
		"subpd		%%xmm3		,%%xmm0						\n\t		subpd		%%xmm11,%%xmm8				\n\t"\
		"movaps		%%xmm0		,%%xmm2						\n\t		movaps		%%xmm8,%%xmm10				\n\t"\
		"movaps		%%xmm1		,%%xmm3						\n\t		movaps		%%xmm9,%%xmm11				\n\t"\
		"addq		$0x020		,%%rdx																			\n\t"\
		"movaps		0x100(%%rax),%%xmm6						\n\t		movaps		0x140(%%rax),%%xmm14		\n\t"\
		"movaps		0x100(%%rax),%%xmm4						\n\t		movaps		0x140(%%rax),%%xmm12		\n\t"\
		"unpckhpd	0x100(%%rbx),%%xmm6						\n\t		unpckhpd	0x140(%%rbx),%%xmm14		\n\t"\
		"unpcklpd	0x100(%%rbx),%%xmm4						\n\t		unpcklpd	0x140(%%rbx),%%xmm12		\n\t"\
		"movaps		%%xmm6		,0x220(%%rcx)				\n\t		movaps		%%xmm14,0x2a0(%%rcx)		\n\t"\
		"movaps		0x110(%%rax),%%xmm7						\n\t		movaps		0x150(%%rax),%%xmm15		\n\t"\
		"movaps		0x110(%%rax),%%xmm5						\n\t		movaps		0x150(%%rax),%%xmm13		\n\t"\
		"unpckhpd	0x110(%%rbx),%%xmm7						\n\t		unpckhpd	0x150(%%rbx),%%xmm15		\n\t"\
		"unpcklpd	0x110(%%rbx),%%xmm5						\n\t		unpcklpd	0x150(%%rbx),%%xmm13		\n\t"\
		"movaps		%%xmm7		,0x230(%%rcx)				\n\t		movaps		%%xmm15,0x2b0(%%rcx)		\n\t"\
		"movaps		%%xmm4		,%%xmm6						\n\t		movaps		%%xmm12,%%xmm14				\n\t"\
		"movaps		%%xmm5		,%%xmm7						\n\t		movaps		%%xmm13,%%xmm15				\n\t"\
		"mulpd		     (%%rdx),%%xmm4						\n\t		mulpd		0x080(%%rdx),%%xmm12		\n\t"\
		"mulpd		     (%%rdx),%%xmm5						\n\t		mulpd		0x080(%%rdx),%%xmm13		\n\t"\
		"mulpd		0x010(%%rdx),%%xmm6						\n\t		mulpd		0x090(%%rdx),%%xmm14		\n\t"\
		"mulpd		0x010(%%rdx),%%xmm7						\n\t		mulpd		0x090(%%rdx),%%xmm15		\n\t"\
		"addpd		%%xmm6		,%%xmm5						\n\t		addpd		%%xmm14,%%xmm13				\n\t"\
		"subpd		%%xmm7		,%%xmm4						\n\t		subpd		%%xmm15,%%xmm12				\n\t"\
		"addpd		%%xmm4		,%%xmm0						\n\t		addpd		%%xmm12,%%xmm8				\n\t"\
		"addpd		%%xmm5		,%%xmm1						\n\t		addpd		%%xmm13,%%xmm9				\n\t"\
		"subpd		%%xmm4		,%%xmm2						\n\t		subpd		%%xmm12,%%xmm10				\n\t"\
		"subpd		%%xmm5		,%%xmm3						\n\t		subpd		%%xmm13,%%xmm11				\n\t"\
		"addq		$0x040		,%%rdx																			\n\t"\
		"movaps		0x180(%%rax),%%xmm6						\n\t		movaps		0x1c0(%%rax),%%xmm14		\n\t"\
		"movaps		0x180(%%rax),%%xmm4						\n\t		movaps		0x1c0(%%rax),%%xmm12		\n\t"\
		"unpckhpd	0x180(%%rbx),%%xmm6						\n\t		unpckhpd	0x1c0(%%rbx),%%xmm14		\n\t"\
		"unpcklpd	0x180(%%rbx),%%xmm4						\n\t		unpcklpd	0x1c0(%%rbx),%%xmm12		\n\t"\
		"movaps		%%xmm6		,0x260(%%rcx)				\n\t		movaps		%%xmm14,0x2e0(%%rcx)		\n\t"\
		"movaps		0x190(%%rax),%%xmm7						\n\t		movaps		0x1d0(%%rax),%%xmm15		\n\t"\
		"movaps		0x190(%%rax),%%xmm5						\n\t		movaps		0x1d0(%%rax),%%xmm13		\n\t"\
		"unpckhpd	0x190(%%rbx),%%xmm7						\n\t		unpckhpd	0x1d0(%%rbx),%%xmm15		\n\t"\
		"unpcklpd	0x190(%%rbx),%%xmm5						\n\t		unpcklpd	0x1d0(%%rbx),%%xmm13		\n\t"\
		"movaps		%%xmm7		,0x270(%%rcx)				\n\t		movaps		%%xmm15,0x2f0(%%rcx)		\n\t"\
		"movaps		%%xmm4		,%%xmm6						\n\t		movaps		%%xmm12,%%xmm14				\n\t"\
		"movaps		%%xmm5		,%%xmm7						\n\t		movaps		%%xmm13,%%xmm15				\n\t"\
		"mulpd		     (%%rdx),%%xmm4						\n\t		mulpd		0x080(%%rdx),%%xmm12		\n\t"\
		"mulpd		     (%%rdx),%%xmm5						\n\t		mulpd		0x080(%%rdx),%%xmm13		\n\t"\
		"mulpd		0x010(%%rdx),%%xmm6						\n\t		mulpd		0x090(%%rdx),%%xmm14		\n\t"\
		"mulpd		0x010(%%rdx),%%xmm7						\n\t		mulpd		0x090(%%rdx),%%xmm15		\n\t"\
		"addpd		%%xmm6		,%%xmm5						\n\t		addpd		%%xmm14,%%xmm13				\n\t"\
		"subpd		%%xmm7		,%%xmm4						\n\t		subpd		%%xmm15,%%xmm12				\n\t"\
		"movaps		%%xmm5		,0x010(%%rcx)				\n\t		movaps		%%xmm13,0x090(%%rcx)		\n\t"\
		"movaps		%%xmm4		,     (%%rcx)				\n\t		movaps		%%xmm12,0x080(%%rcx)		\n\t"\
		"subq		$0x020		,%%rdx																			\n\t"\
		"movaps		0x080(%%rax),%%xmm6						\n\t		movaps		0x0c0(%%rax),%%xmm14		\n\t"\
		"movaps		0x080(%%rax),%%xmm4						\n\t		movaps		0x0c0(%%rax),%%xmm12		\n\t"\
		"unpckhpd	0x080(%%rbx),%%xmm6						\n\t		unpckhpd	0x0c0(%%rbx),%%xmm14		\n\t"\
		"unpcklpd	0x080(%%rbx),%%xmm4						\n\t		unpcklpd	0x0c0(%%rbx),%%xmm12		\n\t"\
		"movaps		%%xmm6		,0x240(%%rcx)				\n\t		movaps		%%xmm14,0x2c0(%%rcx)		\n\t"\
		"movaps		0x090(%%rax),%%xmm7						\n\t		movaps		0x0d0(%%rax),%%xmm15		\n\t"\
		"movaps		0x090(%%rax),%%xmm5						\n\t		movaps		0x0d0(%%rax),%%xmm13		\n\t"\
		"unpckhpd	0x090(%%rbx),%%xmm7						\n\t		unpckhpd	0x0d0(%%rbx),%%xmm15		\n\t"\
		"unpcklpd	0x090(%%rbx),%%xmm5						\n\t		unpcklpd	0x0d0(%%rbx),%%xmm13		\n\t"\
		"movaps		%%xmm7		,0x250(%%rcx)				\n\t		movaps		%%xmm15,0x2d0(%%rcx)		\n\t"\
		"movaps		%%xmm4		,%%xmm6						\n\t		movaps		%%xmm12,%%xmm14				\n\t"\
		"movaps		%%xmm5		,%%xmm7						\n\t		movaps		%%xmm13,%%xmm15				\n\t"\
		"mulpd		     (%%rdx),%%xmm4						\n\t		mulpd		0x080(%%rdx),%%xmm12		\n\t"\
		"mulpd		     (%%rdx),%%xmm5						\n\t		mulpd		0x080(%%rdx),%%xmm13		\n\t"\
		"mulpd		0x010(%%rdx),%%xmm6						\n\t		mulpd		0x090(%%rdx),%%xmm14		\n\t"\
		"mulpd		0x010(%%rdx),%%xmm7						\n\t		mulpd		0x090(%%rdx),%%xmm15		\n\t"\
		"addpd		%%xmm6		,%%xmm5						\n\t		addpd		%%xmm14,%%xmm13				\n\t"\
		"subpd		%%xmm7		,%%xmm4						\n\t		subpd		%%xmm15,%%xmm12				\n\t"\
		"movaps		%%xmm5		,%%xmm7						\n\t		movaps		%%xmm13,%%xmm15				\n\t"\
		"movaps		%%xmm4		,%%xmm6						\n\t		movaps		%%xmm12,%%xmm14				\n\t"\
		"subpd		     (%%rcx),%%xmm4						\n\t		subpd		0x080(%%rcx),%%xmm12		\n\t"\
		"subpd		0x010(%%rcx),%%xmm5						\n\t		subpd		0x090(%%rcx),%%xmm13		\n\t"\
		"addpd		     (%%rcx),%%xmm6						\n\t		addpd		0x080(%%rcx),%%xmm14		\n\t"\
		"addpd		0x010(%%rcx),%%xmm7						\n\t		addpd		0x090(%%rcx),%%xmm15		\n\t"\
		"subpd		%%xmm6		,%%xmm0						\n\t		subpd		%%xmm14,%%xmm8				\n\t"\
		"subpd		%%xmm7		,%%xmm1						\n\t		subpd		%%xmm15,%%xmm9				\n\t"\
		"/*movaps	%%xmm0		,0x040(%%rcx)	*/			\n\t	/*	movaps		%%xmm8,0x0c0(%%rcx)		*/	\n\t"\
		"/*movaps	%%xmm1		,0x050(%%rcx)	*/			\n\t	/*	movaps		%%xmm9,0x0d0(%%rcx)		*/	\n\t"\
		"mulpd		(%%rdi)		,%%xmm6						\n\t		mulpd		(%%rdi),%%xmm14				\n\t"\
		"mulpd		(%%rdi)		,%%xmm7						\n\t		mulpd		(%%rdi),%%xmm15				\n\t"\
		"addpd		%%xmm0		,%%xmm6						\n\t		addpd		%%xmm8,%%xmm14				\n\t"\
		"addpd		%%xmm1		,%%xmm7						\n\t		addpd		%%xmm9,%%xmm15				\n\t"\
		"/*movaps	%%xmm6		,     (%%rcx)	*/			\n\t		movaps		%%xmm14,0x080(%%rcx)		\n\t"\
		"/*movaps	%%xmm7		,0x010(%%rcx)	*/			\n\t		movaps		%%xmm15,0x090(%%rcx)		\n\t"\
		"subpd		%%xmm5		,%%xmm2						\n\t		subpd		%%xmm13,%%xmm10				\n\t"\
		"subpd		%%xmm4		,%%xmm3						\n\t		subpd		%%xmm12,%%xmm11				\n\t"\
		"/*movaps	%%xmm2		,0x020(%%rcx)	*/			\n\t		mulpd		(%%rdi),%%xmm13				\n\t"\
		"/*movaps	%%xmm3		,0x070(%%rcx)	*/			\n\t		mulpd		(%%rdi),%%xmm12				\n\t"\
		"mulpd		(%%rdi)		,%%xmm5						\n\t		addpd		%%xmm10,%%xmm13				\n\t"\
		"mulpd		(%%rdi)		,%%xmm4						\n\t		addpd		%%xmm11,%%xmm12				\n\t"\
		"addpd		%%xmm2		,%%xmm5						\n\t		movaps		%%xmm10,%%xmm14				\n\t"\
		"addpd		%%xmm3		,%%xmm4						\n\t		movaps		%%xmm13,%%xmm15				\n\t"\
		"/*movaps	%%xmm5		,0x060(%%rcx)	*/			\n\t		subpd		%%xmm12,%%xmm10				\n\t"\
		"/*movaps	%%xmm4		,0x030(%%rcx)	*/						subpd		%%xmm11,%%xmm13				\n\t"\
		"																addpd		%%xmm12,%%xmm14				\n\t"\
		"																addpd		%%xmm11,%%xmm15				\n\t"\
		"																mulpd		(%%rsi),%%xmm10				\n\t"\
		"																mulpd		(%%rsi),%%xmm13				\n\t"\
		"																mulpd		(%%rsi),%%xmm14				\n\t"\
		"																mulpd		(%%rsi),%%xmm15				\n\t"\
		"															/*	movaps		%%xmm10,0x0a0(%%rcx)	*/	\n\t"\
		"															/*	movaps		%%xmm13,0x0e0(%%rcx)	*/	\n\t"\
		"															/*	movaps		%%xmm14,0x0b0(%%rcx)	*/	\n\t"\
		"															/*	movaps		%%xmm15,0x0f0(%%rcx)	*/	\n\t"\
		"/*****	SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r00)	*****/\n\t"\
		"																movaps		0x080(%%rcx),%%xmm11		\n\t"\
		"																movaps		0x090(%%rcx),%%xmm12		\n\t"\
		"subpd		%%xmm11		,%%xmm6						\n\t		subpd		%%xmm10		,%%xmm2			\n\t"\
		"subpd		%%xmm9		,%%xmm0						\n\t		subpd		%%xmm15		,%%xmm5			\n\t"\
		"subpd		%%xmm12		,%%xmm7						\n\t		subpd		%%xmm14		,%%xmm4			\n\t"\
		"subpd		%%xmm8		,%%xmm1						\n\t		subpd		%%xmm13		,%%xmm3			\n\t"\
		"mulpd		(%%rdi)		,%%xmm11					\n\t		mulpd		(%%rdi)		,%%xmm10		\n\t"\
		"mulpd		(%%rdi)		,%%xmm9						\n\t		mulpd		(%%rdi)		,%%xmm15		\n\t"\
		"mulpd		(%%rdi)		,%%xmm12					\n\t		mulpd		(%%rdi)		,%%xmm14		\n\t"\
		"mulpd		(%%rdi)		,%%xmm8						\n\t		mulpd		(%%rdi)		,%%xmm13		\n\t"\
		"addpd		%%xmm6		,%%xmm11					\n\t		addpd		%%xmm2		,%%xmm10		\n\t"\
		"addpd		%%xmm0		,%%xmm9						\n\t		addpd		%%xmm5		,%%xmm15		\n\t"\
		"addpd		%%xmm7		,%%xmm12					\n\t		addpd		%%xmm4		,%%xmm14		\n\t"\
		"addpd		%%xmm1		,%%xmm8						\n\t		addpd		%%xmm3		,%%xmm13		\n\t"\
		"movaps		%%xmm6		,0x080(%%rcx)				\n\t		movaps		%%xmm2		,0x0a0(%%rcx)	\n\t"\
		"movaps		%%xmm0		,0x040(%%rcx)				\n\t		movaps		%%xmm5		,0x060(%%rcx)	\n\t"\
		"movaps		%%xmm7		,0x090(%%rcx)				\n\t		movaps		%%xmm4		,0x0b0(%%rcx)	\n\t"\
		"movaps		%%xmm1		,0x0d0(%%rcx)				\n\t		movaps		%%xmm3		,0x0f0(%%rcx)	\n\t"\
		"movaps		%%xmm11		,     (%%rcx)				\n\t		movaps		%%xmm10		,0x020(%%rcx)	\n\t"\
		"movaps		%%xmm9		,0x0c0(%%rcx)				\n\t		movaps		%%xmm15		,0x0e0(%%rcx)	\n\t"\
		"movaps		%%xmm12		,0x010(%%rcx)				\n\t		movaps		%%xmm14		,0x030(%%rcx)	\n\t"\
		"movaps		%%xmm8		,0x050(%%rcx)				\n\t		movaps		%%xmm13		,0x070(%%rcx)	\n\t"\
		"\n\t"\
		"/*...Block 2:	*/\n\t"\
		"addq		$0x20		,%%rax\n\t"\
		"addq		$0x20		,%%rbx\n\t"\
		"/*****	SSE2_RADIX4_DIF_4WRAPPER(c02,c0A,c12,c1A,r10) ****	SSE2_RADIX4_DIF_4WRAPPER_2NDOFTWO(c06,c0E,c16,c1E,r18)	*****/\n\t"\
		"movq		%[__r10]	,%%rcx						\n\t		/*addq		$0x080,%%rcx // __r18 */	\n\t"\
		"movq		%[__c02]	,%%rdx						\n\t		/*addq		$0x080,%%rdx // __c06 */	\n\t"\
		"movaps		     (%%rax),%%xmm6						\n\t		movaps		0x40 (%%rax),%%xmm14		\n\t"\
		"movaps		     (%%rax),%%xmm0						\n\t		movaps		0x40 (%%rax),%%xmm8			\n\t"\
		"unpckhpd	     (%%rbx),%%xmm6						\n\t		unpckhpd	0x40 (%%rbx),%%xmm14		\n\t"\
		"unpcklpd	     (%%rbx),%%xmm0						\n\t		unpcklpd	0x40 (%%rbx),%%xmm8			\n\t"\
		"movaps		%%xmm6		,0x200(%%rcx)				\n\t		movaps		%%xmm14,0x280(%%rcx)		\n\t"\
		"movaps		0x010(%%rax),%%xmm7						\n\t		movaps		0x50 (%%rax),%%xmm15		\n\t"\
		"movaps		0x010(%%rax),%%xmm1						\n\t		movaps		0x50 (%%rax),%%xmm9			\n\t"\
		"unpckhpd	0x010(%%rbx),%%xmm7						\n\t		unpckhpd	0x50 (%%rbx),%%xmm15		\n\t"\
		"unpcklpd	0x010(%%rbx),%%xmm1						\n\t		unpcklpd	0x50 (%%rbx),%%xmm9			\n\t"\
		"movaps		%%xmm7		,0x210(%%rcx)				\n\t		movaps		%%xmm15,0x290(%%rcx)		\n\t"\
		"movaps		%%xmm0		,%%xmm2						\n\t		movaps		%%xmm8,%%xmm10				\n\t"\
		"movaps		%%xmm1		,%%xmm3						\n\t		movaps		%%xmm9,%%xmm11				\n\t"\
		"mulpd		     (%%rdx),%%xmm0						\n\t		mulpd		0x080(%%rdx),%%xmm8			\n\t"\
		"mulpd		     (%%rdx),%%xmm1						\n\t		mulpd		0x080(%%rdx),%%xmm9			\n\t"\
		"mulpd		0x010(%%rdx),%%xmm2						\n\t		mulpd		0x090(%%rdx),%%xmm10		\n\t"\
		"mulpd		0x010(%%rdx),%%xmm3						\n\t		mulpd		0x090(%%rdx),%%xmm11		\n\t"\
		"addpd		%%xmm2		,%%xmm1						\n\t		addpd		%%xmm10,%%xmm9				\n\t"\
		"subpd		%%xmm3		,%%xmm0						\n\t		subpd		%%xmm11,%%xmm8				\n\t"\
		"movaps		%%xmm0		,%%xmm2						\n\t		movaps		%%xmm8,%%xmm10				\n\t"\
		"movaps		%%xmm1		,%%xmm3						\n\t		movaps		%%xmm9,%%xmm11				\n\t"\
		"addq		$0x020		,%%rdx																			\n\t"\
		"movaps		0x100(%%rax),%%xmm6						\n\t		movaps		0x140(%%rax),%%xmm14		\n\t"\
		"movaps		0x100(%%rax),%%xmm4						\n\t		movaps		0x140(%%rax),%%xmm12		\n\t"\
		"unpckhpd	0x100(%%rbx),%%xmm6						\n\t		unpckhpd	0x140(%%rbx),%%xmm14		\n\t"\
		"unpcklpd	0x100(%%rbx),%%xmm4						\n\t		unpcklpd	0x140(%%rbx),%%xmm12		\n\t"\
		"movaps		%%xmm6		,0x220(%%rcx)				\n\t		movaps		%%xmm14,0x2a0(%%rcx)		\n\t"\
		"movaps		0x110(%%rax),%%xmm7						\n\t		movaps		0x150(%%rax),%%xmm15		\n\t"\
		"movaps		0x110(%%rax),%%xmm5						\n\t		movaps		0x150(%%rax),%%xmm13		\n\t"\
		"unpckhpd	0x110(%%rbx),%%xmm7						\n\t		unpckhpd	0x150(%%rbx),%%xmm15		\n\t"\
		"unpcklpd	0x110(%%rbx),%%xmm5						\n\t		unpcklpd	0x150(%%rbx),%%xmm13		\n\t"\
		"movaps		%%xmm7		,0x230(%%rcx)				\n\t		movaps		%%xmm15,0x2b0(%%rcx)		\n\t"\
		"movaps		%%xmm4		,%%xmm6						\n\t		movaps		%%xmm12,%%xmm14				\n\t"\
		"movaps		%%xmm5		,%%xmm7						\n\t		movaps		%%xmm13,%%xmm15				\n\t"\
		"mulpd		     (%%rdx),%%xmm4						\n\t		mulpd		0x080(%%rdx),%%xmm12		\n\t"\
		"mulpd		     (%%rdx),%%xmm5						\n\t		mulpd		0x080(%%rdx),%%xmm13		\n\t"\
		"mulpd		0x010(%%rdx),%%xmm6						\n\t		mulpd		0x090(%%rdx),%%xmm14		\n\t"\
		"mulpd		0x010(%%rdx),%%xmm7						\n\t		mulpd		0x090(%%rdx),%%xmm15		\n\t"\
		"addpd		%%xmm6		,%%xmm5						\n\t		addpd		%%xmm14,%%xmm13				\n\t"\
		"subpd		%%xmm7		,%%xmm4						\n\t		subpd		%%xmm15,%%xmm12				\n\t"\
		"addpd		%%xmm4		,%%xmm0						\n\t		addpd		%%xmm12,%%xmm8				\n\t"\
		"addpd		%%xmm5		,%%xmm1						\n\t		addpd		%%xmm13,%%xmm9				\n\t"\
		"subpd		%%xmm4		,%%xmm2						\n\t		subpd		%%xmm12,%%xmm10				\n\t"\
		"subpd		%%xmm5		,%%xmm3						\n\t		subpd		%%xmm13,%%xmm11				\n\t"\
		"addq		$0x040		,%%rdx																			\n\t"\
		"movaps		0x180(%%rax),%%xmm6						\n\t		movaps		0x1c0(%%rax),%%xmm14		\n\t"\
		"movaps		0x180(%%rax),%%xmm4						\n\t		movaps		0x1c0(%%rax),%%xmm12		\n\t"\
		"unpckhpd	0x180(%%rbx),%%xmm6						\n\t		unpckhpd	0x1c0(%%rbx),%%xmm14		\n\t"\
		"unpcklpd	0x180(%%rbx),%%xmm4						\n\t		unpcklpd	0x1c0(%%rbx),%%xmm12		\n\t"\
		"movaps		%%xmm6		,0x260(%%rcx)				\n\t		movaps		%%xmm14,0x2e0(%%rcx)		\n\t"\
		"movaps		0x190(%%rax),%%xmm7						\n\t		movaps		0x1d0(%%rax),%%xmm15		\n\t"\
		"movaps		0x190(%%rax),%%xmm5						\n\t		movaps		0x1d0(%%rax),%%xmm13		\n\t"\
		"unpckhpd	0x190(%%rbx),%%xmm7						\n\t		unpckhpd	0x1d0(%%rbx),%%xmm15		\n\t"\
		"unpcklpd	0x190(%%rbx),%%xmm5						\n\t		unpcklpd	0x1d0(%%rbx),%%xmm13		\n\t"\
		"movaps		%%xmm7		,0x270(%%rcx)				\n\t		movaps		%%xmm15,0x2f0(%%rcx)		\n\t"\
		"movaps		%%xmm4		,%%xmm6						\n\t		movaps		%%xmm12,%%xmm14				\n\t"\
		"movaps		%%xmm5		,%%xmm7						\n\t		movaps		%%xmm13,%%xmm15				\n\t"\
		"mulpd		     (%%rdx),%%xmm4						\n\t		mulpd		0x080(%%rdx),%%xmm12		\n\t"\
		"mulpd		     (%%rdx),%%xmm5						\n\t		mulpd		0x080(%%rdx),%%xmm13		\n\t"\
		"mulpd		0x010(%%rdx),%%xmm6						\n\t		mulpd		0x090(%%rdx),%%xmm14		\n\t"\
		"mulpd		0x010(%%rdx),%%xmm7						\n\t		mulpd		0x090(%%rdx),%%xmm15		\n\t"\
		"addpd		%%xmm6		,%%xmm5						\n\t		addpd		%%xmm14,%%xmm13				\n\t"\
		"subpd		%%xmm7		,%%xmm4						\n\t		subpd		%%xmm15,%%xmm12				\n\t"\
		"movaps		%%xmm5		,0x010(%%rcx)				\n\t		movaps		%%xmm13,0x090(%%rcx)		\n\t"\
		"movaps		%%xmm4		,     (%%rcx)				\n\t		movaps		%%xmm12,0x080(%%rcx)		\n\t"\
		"subq		$0x020		,%%rdx																			\n\t"\
		"movaps		0x080(%%rax),%%xmm6						\n\t		movaps		0x0c0(%%rax),%%xmm14		\n\t"\
		"movaps		0x080(%%rax),%%xmm4						\n\t		movaps		0x0c0(%%rax),%%xmm12		\n\t"\
		"unpckhpd	0x080(%%rbx),%%xmm6						\n\t		unpckhpd	0x0c0(%%rbx),%%xmm14		\n\t"\
		"unpcklpd	0x080(%%rbx),%%xmm4						\n\t		unpcklpd	0x0c0(%%rbx),%%xmm12		\n\t"\
		"movaps		%%xmm6		,0x240(%%rcx)				\n\t		movaps		%%xmm14,0x2c0(%%rcx)		\n\t"\
		"movaps		0x090(%%rax),%%xmm7						\n\t		movaps		0x0d0(%%rax),%%xmm15		\n\t"\
		"movaps		0x090(%%rax),%%xmm5						\n\t		movaps		0x0d0(%%rax),%%xmm13		\n\t"\
		"unpckhpd	0x090(%%rbx),%%xmm7						\n\t		unpckhpd	0x0d0(%%rbx),%%xmm15		\n\t"\
		"unpcklpd	0x090(%%rbx),%%xmm5						\n\t		unpcklpd	0x0d0(%%rbx),%%xmm13		\n\t"\
		"movaps		%%xmm7		,0x250(%%rcx)				\n\t		movaps		%%xmm15,0x2d0(%%rcx)		\n\t"\
		"movaps		%%xmm4		,%%xmm6						\n\t		movaps		%%xmm12,%%xmm14				\n\t"\
		"movaps		%%xmm5		,%%xmm7						\n\t		movaps		%%xmm13,%%xmm15				\n\t"\
		"mulpd		     (%%rdx),%%xmm4						\n\t		mulpd		0x080(%%rdx),%%xmm12		\n\t"\
		"mulpd		     (%%rdx),%%xmm5						\n\t		mulpd		0x080(%%rdx),%%xmm13		\n\t"\
		"mulpd		0x010(%%rdx),%%xmm6						\n\t		mulpd		0x090(%%rdx),%%xmm14		\n\t"\
		"mulpd		0x010(%%rdx),%%xmm7						\n\t		mulpd		0x090(%%rdx),%%xmm15		\n\t"\
		"addpd		%%xmm6		,%%xmm5						\n\t		addpd		%%xmm14,%%xmm13				\n\t"\
		"subpd		%%xmm7		,%%xmm4						\n\t		subpd		%%xmm15,%%xmm12				\n\t"\
		"movaps		%%xmm5		,%%xmm7						\n\t		movaps		%%xmm13,%%xmm15				\n\t"\
		"movaps		%%xmm4		,%%xmm6						\n\t		movaps		%%xmm12,%%xmm14				\n\t"\
		"subpd		     (%%rcx),%%xmm4						\n\t		subpd		0x080(%%rcx),%%xmm12		\n\t"\
		"subpd		0x010(%%rcx),%%xmm5						\n\t		subpd		0x090(%%rcx),%%xmm13		\n\t"\
		"addpd		     (%%rcx),%%xmm6						\n\t		addpd		0x080(%%rcx),%%xmm14		\n\t"\
		"addpd		0x010(%%rcx),%%xmm7						\n\t		addpd		0x090(%%rcx),%%xmm15		\n\t"\
		"subpd		%%xmm6		,%%xmm0						\n\t		subpd		%%xmm14,%%xmm8				\n\t"\
		"subpd		%%xmm7		,%%xmm1						\n\t		subpd		%%xmm15,%%xmm9				\n\t"\
		"/*movaps	%%xmm0		,0x040(%%rcx)	*/			\n\t	/*	movaps		%%xmm8,0x0c0(%%rcx)		*/	\n\t"\
		"/*movaps	%%xmm1		,0x050(%%rcx)	*/			\n\t	/*	movaps		%%xmm9,0x0d0(%%rcx)		*/	\n\t"\
		"mulpd		(%%rdi)		,%%xmm6						\n\t		mulpd		(%%rdi),%%xmm14				\n\t"\
		"mulpd		(%%rdi)		,%%xmm7						\n\t		mulpd		(%%rdi),%%xmm15				\n\t"\
		"addpd		%%xmm0		,%%xmm6						\n\t		addpd		%%xmm8,%%xmm14				\n\t"\
		"addpd		%%xmm1		,%%xmm7						\n\t		addpd		%%xmm9,%%xmm15				\n\t"\
		"/*movaps	%%xmm6		,     (%%rcx)	*/			\n\t		movaps		%%xmm14,0x080(%%rcx)		\n\t"\
		"/*movaps	%%xmm7		,0x010(%%rcx)	*/			\n\t		movaps		%%xmm15,0x090(%%rcx)		\n\t"\
		"subpd		%%xmm5		,%%xmm2						\n\t		subpd		%%xmm13,%%xmm10				\n\t"\
		"subpd		%%xmm4		,%%xmm3						\n\t		subpd		%%xmm12,%%xmm11				\n\t"\
		"/*movaps	%%xmm2		,0x020(%%rcx)	*/			\n\t		mulpd		(%%rdi),%%xmm13				\n\t"\
		"/*movaps	%%xmm3		,0x070(%%rcx)	*/			\n\t		mulpd		(%%rdi),%%xmm12				\n\t"\
		"mulpd		(%%rdi)		,%%xmm5						\n\t		addpd		%%xmm10,%%xmm13				\n\t"\
		"mulpd		(%%rdi)		,%%xmm4						\n\t		addpd		%%xmm11,%%xmm12				\n\t"\
		"addpd		%%xmm2		,%%xmm5						\n\t		movaps		%%xmm10,%%xmm14				\n\t"\
		"addpd		%%xmm3		,%%xmm4						\n\t		movaps		%%xmm13,%%xmm15				\n\t"\
		"/*movaps	%%xmm5		,0x060(%%rcx)	*/			\n\t		subpd		%%xmm12,%%xmm10				\n\t"\
		"/*movaps	%%xmm4		,0x030(%%rcx)	*/			\n\t		subpd		%%xmm11,%%xmm13				\n\t"\
		"																addpd		%%xmm12,%%xmm14				\n\t"\
		"																addpd		%%xmm11,%%xmm15				\n\t"\
		"																mulpd		(%%rsi),%%xmm10	/* isrt2 */	\n\t"\
		"																mulpd		(%%rsi),%%xmm13				\n\t"\
		"																mulpd		(%%rsi),%%xmm14				\n\t"\
		"																mulpd		(%%rsi),%%xmm15				\n\t"\
		"															/*	movaps		%%xmm10,0x0a0(%%rcx)	*/	\n\t"\
		"															/*	movaps		%%xmm13,0x0e0(%%rcx)	*/	\n\t"\
		"															/*	movaps		%%xmm14,0x0b0(%%rcx)	*/	\n\t"\
		"															/*	movaps		%%xmm15,0x0f0(%%rcx)	*/	\n\t"\
		"/*****	SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r10)	*****/\n\t"\
		"																movaps		0x080(%%rcx),%%xmm11		\n\t"\
		"																movaps		0x090(%%rcx),%%xmm12		\n\t"\
		"subpd		%%xmm11		,%%xmm6						\n\t		subpd		%%xmm10		,%%xmm2			\n\t"\
		"subpd		%%xmm9		,%%xmm0						\n\t		subpd		%%xmm15		,%%xmm5			\n\t"\
		"subpd		%%xmm12		,%%xmm7						\n\t		subpd		%%xmm14		,%%xmm4			\n\t"\
		"subpd		%%xmm8		,%%xmm1						\n\t		subpd		%%xmm13		,%%xmm3			\n\t"\
		"mulpd		(%%rdi)	,%%xmm11						\n\t		mulpd		(%%rdi)		,%%xmm10		\n\t"\
		"mulpd		(%%rdi)	,%%xmm9							\n\t		mulpd		(%%rdi)		,%%xmm15		\n\t"\
		"mulpd		(%%rdi)	,%%xmm12						\n\t		mulpd		(%%rdi)		,%%xmm14		\n\t"\
		"mulpd		(%%rdi)	,%%xmm8							\n\t		mulpd		(%%rdi)		,%%xmm13		\n\t"\
		"addpd		%%xmm6		,%%xmm11					\n\t		addpd		%%xmm2		,%%xmm10		\n\t"\
		"addpd		%%xmm0		,%%xmm9						\n\t		addpd		%%xmm5		,%%xmm15		\n\t"\
		"addpd		%%xmm7		,%%xmm12					\n\t		addpd		%%xmm4		,%%xmm14		\n\t"\
		"addpd		%%xmm1		,%%xmm8						\n\t		addpd		%%xmm3		,%%xmm13		\n\t"\
		"movaps		%%xmm6		,0x080(%%rcx)				\n\t		movaps		%%xmm2		,0x0a0(%%rcx)	\n\t"\
		"movaps		%%xmm0		,0x040(%%rcx)				\n\t		movaps		%%xmm5		,0x060(%%rcx)	\n\t"\
		"movaps		%%xmm7		,0x090(%%rcx)				\n\t		movaps		%%xmm4		,0x0b0(%%rcx)	\n\t"\
		"movaps		%%xmm1		,0x0d0(%%rcx)				\n\t		movaps		%%xmm3		,0x0f0(%%rcx)	\n\t"\
		"movaps		%%xmm11		,     (%%rcx)				\n\t		movaps		%%xmm10		,0x020(%%rcx)	\n\t"\
		"movaps		%%xmm9		,0x0c0(%%rcx)				\n\t		movaps		%%xmm15		,0x0e0(%%rcx)	\n\t"\
		"movaps		%%xmm12		,0x010(%%rcx)				\n\t		movaps		%%xmm14		,0x030(%%rcx)	\n\t"\
		"movaps		%%xmm8		,0x050(%%rcx)				\n\t		movaps		%%xmm13		,0x070(%%rcx)	\n\t"\
		"/********************************************************************************************************\n\t"\
		" Next 2 blocks operate on odd-indexed elements from the unpck*pd commands which we stored to temporaries:\n\t"\
		"********************************************************************************************************/\n\t"\
		"/*...Block 3:	*/\n\t"\
		"/*	SSE2_RADIX4_DIF_4TWIDDLE         (r20,r24,r22,r26,r20,c01) */\n\t"\
		"addq		$0x100		,%%rcx		/*****	SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO(r28,r2C,r2A,r2E,r28,c05)	*****/\n\t"\
		"movq		%[__c01]	,%%rbx						\n\t	/*	movq		%[__c05]	,%%rbx	*/		\n\t"\
		"movq		%%rcx		,%%rax						\n\t	/*	addq		$0x040		,%%rax	*/		\n\t"\
		"addq		$0x020		,%%rcx						\n\t	/*	addq		$0x040		,%%rcx	*/		\n\t"\
		"movaps		     (%%rax),%%xmm0	\n\t	movq %%rax,%%rdx \n\t	movaps		0x080(%%rax),%%xmm8			\n\t"\
		"movaps		     (%%rcx),%%xmm4						\n\t		movaps		0x080(%%rcx),%%xmm12		\n\t"\
		"movaps		0x010(%%rax),%%xmm1						\n\t		movaps		0x090(%%rax),%%xmm9			\n\t"\
		"movaps		0x010(%%rcx),%%xmm5						\n\t		movaps		0x090(%%rcx),%%xmm13		\n\t"\
		"movaps		     (%%rbx),%%xmm6						\n\t		movaps		0x080(%%rbx),%%xmm14		\n\t"\
		"movaps		0x010(%%rbx),%%xmm7						\n\t		movaps		0x090(%%rbx),%%xmm15		\n\t"\
		"movaps		%%xmm0		,%%xmm2						\n\t		movaps		%%xmm8		,%%xmm10		\n\t"\
		"movaps		%%xmm1		,%%xmm3						\n\t		movaps		%%xmm9		,%%xmm11		\n\t"\
		"mulpd		%%xmm6		,%%xmm0						\n\t		mulpd		%%xmm14		,%%xmm8			\n\t"\
		"mulpd		%%xmm6		,%%xmm1						\n\t		mulpd		%%xmm14		,%%xmm9			\n\t"\
		"mulpd		%%xmm7		,%%xmm2						\n\t		mulpd		%%xmm15		,%%xmm10		\n\t"\
		"mulpd		%%xmm7		,%%xmm3						\n\t		mulpd		%%xmm15		,%%xmm11		\n\t"\
		"movaps		%%xmm4		,%%xmm6						\n\t		movaps		%%xmm12		,%%xmm14		\n\t"\
		"addpd		%%xmm2		,%%xmm1						\n\t		addpd		%%xmm10		,%%xmm9			\n\t"\
		"movaps		%%xmm5		,%%xmm7						\n\t		movaps		%%xmm13		,%%xmm15		\n\t"\
		"mulpd		0x20 (%%rbx),%%xmm4						\n\t		mulpd		0xa0 (%%rbx),%%xmm12		\n\t"\
		"subpd		%%xmm3		,%%xmm0						\n\t		subpd		%%xmm11		,%%xmm8			\n\t"\
		"mulpd		0x20 (%%rbx),%%xmm5						\n\t		mulpd		0xa0 (%%rbx),%%xmm13		\n\t"\
		"mulpd		0x30 (%%rbx),%%xmm6						\n\t		mulpd		0xb0 (%%rbx),%%xmm14		\n\t"\
		"movaps		%%xmm0		,%%xmm2						\n\t		movaps		%%xmm8		,%%xmm10		\n\t"\
		"mulpd		0x30 (%%rbx),%%xmm7						\n\t		mulpd		0xb0 (%%rbx),%%xmm15		\n\t"\
		"addpd		%%xmm6		,%%xmm5						\n\t		addpd		%%xmm14		,%%xmm13		\n\t"\
		"movaps		%%xmm1		,%%xmm3						\n\t		movaps		%%xmm9		,%%xmm11		\n\t"\
		"subpd		%%xmm7		,%%xmm4						\n\t		subpd		%%xmm15		,%%xmm12		\n\t"\
		"addq		$0x040		,%%rcx						\n\t		addpd		%%xmm12		,%%xmm8			\n\t"\
		"addq		$0x60		,%%rbx						\n\t		addpd		%%xmm13		,%%xmm9			\n\t"\
		"movaps		     (%%rcx),%%xmm6						\n\t		subpd		%%xmm12		,%%xmm10		\n\t"\
		"movaps		0x010(%%rcx),%%xmm7						\n\t		subpd		%%xmm13		,%%xmm11		\n\t"\
		"addpd		%%xmm4		,%%xmm0						\n\t		movaps		0x080(%%rcx),%%xmm12		\n\t"\
		"addpd		%%xmm5		,%%xmm1						\n\t		movaps		0x090(%%rcx),%%xmm13		\n\t"\
		"subpd		%%xmm4		,%%xmm2						\n\t		movaps		0x080(%%rcx),%%xmm14		\n\t"\
		"subpd		%%xmm5		,%%xmm3						\n\t		movaps		0x090(%%rcx),%%xmm15		\n\t"\
		"movaps		%%xmm6		,%%xmm4						\n\t		mulpd		0x080(%%rbx),%%xmm12		\n\t"\
		"movaps		%%xmm7		,%%xmm5						\n\t		mulpd		0x080(%%rbx),%%xmm13		\n\t"\
		"mulpd		     (%%rbx),%%xmm4						\n\t		mulpd		0x090(%%rbx),%%xmm14		\n\t"\
		"mulpd		     (%%rbx),%%xmm5						\n\t		mulpd		0x090(%%rbx),%%xmm15		\n\t"\
		"mulpd		0x010(%%rbx),%%xmm6						\n\t		addpd		%%xmm14		,%%xmm13		\n\t"\
		"mulpd		0x010(%%rbx),%%xmm7						\n\t		subpd		%%xmm15		,%%xmm12		\n\t"\
		"/*movq		%%rax		,%%rdx		*/				\n\t		movaps		%%xmm13		,0x090(%%rdx)	\n\t"\
		"addpd		%%xmm6		,%%xmm5						\n\t		movaps		%%xmm12		,0x080(%%rdx)	\n\t"\
		"subpd		%%xmm7		,%%xmm4						\n\t	/*	subq	$0x20	,%%rbx */				\n\t"\
		"movaps		%%xmm5		,0x010(%%rdx)	\n\t	addq	$0x040		,%%rax								\n\t"\
		"movaps		%%xmm4		,     (%%rdx)	\n\t	subq	$0x20		,%%rbx								\n\t"\
		"movaps		     (%%rax),%%xmm4						\n\t		movaps		0x080(%%rax),%%xmm12		\n\t"\
		"movaps		0x010(%%rax),%%xmm5						\n\t		movaps		0x090(%%rax),%%xmm13		\n\t"\
		"movaps		%%xmm4		,%%xmm6						\n\t		movaps		0x080(%%rax),%%xmm14		\n\t"\
		"movaps		%%xmm5		,%%xmm7						\n\t		movaps		0x090(%%rax),%%xmm15		\n\t"\
		"mulpd		     (%%rbx),%%xmm4						\n\t		mulpd		0x080(%%rbx),%%xmm12		\n\t"\
		"mulpd		     (%%rbx),%%xmm5						\n\t		mulpd		0x080(%%rbx),%%xmm13		\n\t"\
		"mulpd		0x010(%%rbx),%%xmm6						\n\t		mulpd		0x090(%%rbx),%%xmm14		\n\t"\
		"mulpd		0x010(%%rbx),%%xmm7						\n\t		mulpd		0x090(%%rbx),%%xmm15		\n\t"\
		"addpd		%%xmm6		,%%xmm5						\n\t		addpd		%%xmm14		,%%xmm13		\n\t"\
		"subpd		%%xmm7		,%%xmm4						\n\t		subpd		%%xmm15		,%%xmm12		\n\t"\
		"movaps		%%xmm5		,%%xmm7						\n\t		movaps		%%xmm13		,%%xmm15		\n\t"\
		"movaps		%%xmm4		,%%xmm6						\n\t		movaps		%%xmm12		,%%xmm14		\n\t"\
		"subpd		     (%%rdx),%%xmm4						\n\t		subpd		0x080(%%rdx),%%xmm12		\n\t"\
		"subpd		0x010(%%rdx),%%xmm5						\n\t		subpd		0x090(%%rdx),%%xmm13		\n\t"\
		"addpd		     (%%rdx),%%xmm6						\n\t		addpd		0x080(%%rdx),%%xmm14		\n\t"\
		"addpd		0x010(%%rdx),%%xmm7						\n\t		addpd		0x090(%%rdx),%%xmm15		\n\t"\
		"subpd		%%xmm6		,%%xmm0						\n\t		subpd		%%xmm14		,%%xmm8			\n\t"\
		"subpd		%%xmm5		,%%xmm2						\n\t		subpd		%%xmm15		,%%xmm9			\n\t"\
		"subpd		%%xmm7		,%%xmm1						\n\t	/*	movaps		%%xmm8		,0x0c0(%%rdx)*/	\n\t"\
		"subpd		%%xmm4		,%%xmm3						\n\t		subpd		%%xmm13		,%%xmm10		\n\t"\
		"/*movaps		%%xmm0		,0x040(%%rdx)	*/		\n\t	/*	movaps		%%xmm9		,0x0d0(%%rdx)*/	\n\t"\
		"/*movaps		%%xmm2		,0x020(%%rdx)	*/		\n\t		subpd		%%xmm12		,%%xmm11		\n\t"\
		"/*movaps		%%xmm1		,0x050(%%rdx)	*/		\n\t		mulpd		(%%rdi)		,%%xmm14		\n\t"\
		"/*movaps		%%xmm3		,0x070(%%rdx)	*/		\n\t		mulpd		(%%rdi)		,%%xmm13		\n\t"\
		"mulpd		(%%rdi)		,%%xmm6						\n\t		mulpd		(%%rdi)		,%%xmm15		\n\t"\
		"mulpd		(%%rdi)		,%%xmm5						\n\t		mulpd		(%%rdi)		,%%xmm12		\n\t"\
		"mulpd		(%%rdi)		,%%xmm7						\n\t		addpd		%%xmm8		,%%xmm14		\n\t"\
		"mulpd		(%%rdi)		,%%xmm4						\n\t		addpd		%%xmm10		,%%xmm13		\n\t"\
		"addpd		%%xmm0		,%%xmm6						\n\t		addpd		%%xmm9		,%%xmm15		\n\t"\
		"addpd		%%xmm2		,%%xmm5						\n\t		addpd		%%xmm11		,%%xmm12		\n\t"\
		"addpd		%%xmm1		,%%xmm7						\n\t		movaps		%%xmm14		,0x080(%%rdx)	\n\t"\
		"addpd		%%xmm3		,%%xmm4						\n\t		movaps		%%xmm15		,0x090(%%rdx)	\n\t"\
		"/*movaps		%%xmm6		,     (%%rdx)	*/		\n\t		movaps		%%xmm10		,%%xmm14		\n\t"\
		"/*movaps		%%xmm5		,0x060(%%rdx)	*/		\n\t		movaps		%%xmm13		,%%xmm15		\n\t"\
		"/*movaps		%%xmm7		,0x010(%%rdx)	*/		\n\t		subpd		%%xmm12		,%%xmm10		\n\t"\
		"/*movaps		%%xmm4		,0x030(%%rdx)	*/		\n\t		subpd		%%xmm11		,%%xmm13		\n\t"\
		"													\n\t		addpd		%%xmm12		,%%xmm14		\n\t"\
		"													\n\t		addpd		%%xmm11		,%%xmm15		\n\t"\
		"													\n\t		mulpd		(%%rsi)		,%%xmm10		\n\t"\
		"													\n\t		mulpd		(%%rsi)		,%%xmm13		\n\t"\
		"													\n\t		mulpd		(%%rsi)		,%%xmm14		\n\t"\
		"													\n\t		mulpd		(%%rsi)		,%%xmm15		\n\t"\
		"													\n\t	/*	movaps		%%xmm10		,0x0a0(%%rdx)*/	\n\t"\
		"													\n\t	/*	movaps		%%xmm13		,0x0e0(%%rdx)*/	\n\t"\
		"													\n\t	/*	movaps		%%xmm14		,0x0b0(%%rdx)*/	\n\t"\
		"													\n\t	/*	movaps		%%xmm15		,0x0f0(%%rdx)*/	\n\t"\
		"/*****	SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r20)	*****/\n\t"\
		"																movaps		0x080(%%rdx),%%xmm11		\n\t"\
		"																movaps		0x090(%%rdx),%%xmm12		\n\t"\
		"subpd		%%xmm11		,%%xmm6						\n\t		subpd		%%xmm10		,%%xmm2			\n\t"\
		"subpd		%%xmm9		,%%xmm0						\n\t		subpd		%%xmm15		,%%xmm5			\n\t"\
		"subpd		%%xmm12		,%%xmm7						\n\t		subpd		%%xmm14		,%%xmm4			\n\t"\
		"subpd		%%xmm8		,%%xmm1						\n\t		subpd		%%xmm13		,%%xmm3			\n\t"\
		"mulpd		(%%rdi)		,%%xmm11					\n\t		mulpd		(%%rdi)		,%%xmm10		\n\t"\
		"mulpd		(%%rdi)		,%%xmm9						\n\t		mulpd		(%%rdi)		,%%xmm15		\n\t"\
		"mulpd		(%%rdi)		,%%xmm12					\n\t		mulpd		(%%rdi)		,%%xmm14		\n\t"\
		"mulpd		(%%rdi)		,%%xmm8						\n\t		mulpd		(%%rdi)		,%%xmm13		\n\t"\
		"addpd		%%xmm6		,%%xmm11					\n\t		addpd		%%xmm2		,%%xmm10		\n\t"\
		"addpd		%%xmm0		,%%xmm9						\n\t		addpd		%%xmm5		,%%xmm15		\n\t"\
		"addpd		%%xmm7		,%%xmm12					\n\t		addpd		%%xmm4		,%%xmm14		\n\t"\
		"addpd		%%xmm1		,%%xmm8						\n\t		addpd		%%xmm3		,%%xmm13		\n\t"\
		"movaps		%%xmm6		,0x080(%%rdx)				\n\t		movaps		%%xmm2		,0x0a0(%%rdx)	\n\t"\
		"movaps		%%xmm0		,0x040(%%rdx)				\n\t		movaps		%%xmm5		,0x060(%%rdx)	\n\t"\
		"movaps		%%xmm7		,0x090(%%rdx)				\n\t		movaps		%%xmm4		,0x0b0(%%rdx)	\n\t"\
		"movaps		%%xmm1		,0x0d0(%%rdx)				\n\t		movaps		%%xmm3		,0x0f0(%%rdx)	\n\t"\
		"movaps		%%xmm11		,     (%%rdx)				\n\t		movaps		%%xmm10		,0x020(%%rdx)	\n\t"\
		"movaps		%%xmm9		,0x0c0(%%rdx)				\n\t		movaps		%%xmm15		,0x0e0(%%rdx)	\n\t"\
		"movaps		%%xmm12		,0x010(%%rdx)				\n\t		movaps		%%xmm14		,0x030(%%rdx)	\n\t"\
		"movaps		%%xmm8		,0x050(%%rdx)				\n\t		movaps		%%xmm13		,0x070(%%rdx)	\n\t"\
		"/*...Block 4:	*/\n\t"\
		"/*****	SSE2_RADIX4_DIF_4TWIDDLE         (r30,r34,r32,r36,r30,c03)	*****/\n\t"\
		"movq		%[__c03]	,%%rbx					\n\t"\
		"movq		%[__r30]	,%%rax					\n\t"\
		"movq		%%rax		,%%rcx		/*****	SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO(r38,r3C,r3A,r3E,r38,c07)	*****/\n\t"\
		"addq		$0x020		,%%rcx					\n\t"\
		"movaps		     (%%rax),%%xmm0	\n\t	movq %%rax,%%rdx \n\t	movaps		0x080(%%rax),%%xmm8			\n\t"\
		"movaps		     (%%rcx),%%xmm4						\n\t		movaps		0x080(%%rcx),%%xmm12		\n\t"\
		"movaps		0x010(%%rax),%%xmm1						\n\t		movaps		0x090(%%rax),%%xmm9			\n\t"\
		"movaps		0x010(%%rcx),%%xmm5						\n\t		movaps		0x090(%%rcx),%%xmm13		\n\t"\
		"movaps		     (%%rbx),%%xmm6						\n\t		movaps		0x080(%%rbx),%%xmm14		\n\t"\
		"movaps		0x010(%%rbx),%%xmm7						\n\t		movaps		0x090(%%rbx),%%xmm15		\n\t"\
		"movaps		%%xmm0		,%%xmm2						\n\t		movaps		%%xmm8		,%%xmm10		\n\t"\
		"movaps		%%xmm1		,%%xmm3						\n\t		movaps		%%xmm9		,%%xmm11		\n\t"\
		"mulpd		%%xmm6		,%%xmm0						\n\t		mulpd		%%xmm14		,%%xmm8			\n\t"\
		"mulpd		%%xmm6		,%%xmm1						\n\t		mulpd		%%xmm14		,%%xmm9			\n\t"\
		"mulpd		%%xmm7		,%%xmm2						\n\t		mulpd		%%xmm15		,%%xmm10		\n\t"\
		"mulpd		%%xmm7		,%%xmm3						\n\t		mulpd		%%xmm15		,%%xmm11		\n\t"\
		"movaps		%%xmm4		,%%xmm6						\n\t		movaps		%%xmm12		,%%xmm14		\n\t"\
		"addpd		%%xmm2		,%%xmm1						\n\t		addpd		%%xmm10		,%%xmm9			\n\t"\
		"movaps		%%xmm5		,%%xmm7						\n\t		movaps		%%xmm13		,%%xmm15		\n\t"\
		"mulpd		0x20 (%%rbx),%%xmm4						\n\t		mulpd		0xa0 (%%rbx),%%xmm12		\n\t"\
		"subpd		%%xmm3		,%%xmm0						\n\t		subpd		%%xmm11		,%%xmm8			\n\t"\
		"mulpd		0x20 (%%rbx),%%xmm5						\n\t		mulpd		0xa0 (%%rbx),%%xmm13		\n\t"\
		"mulpd		0x30 (%%rbx),%%xmm6						\n\t		mulpd		0xb0 (%%rbx),%%xmm14		\n\t"\
		"movaps		%%xmm0		,%%xmm2						\n\t		movaps		%%xmm8		,%%xmm10		\n\t"\
		"mulpd		0x30 (%%rbx),%%xmm7						\n\t		mulpd		0xb0 (%%rbx),%%xmm15		\n\t"\
		"addpd		%%xmm6		,%%xmm5						\n\t		addpd		%%xmm14		,%%xmm13		\n\t"\
		"movaps		%%xmm1		,%%xmm3						\n\t		movaps		%%xmm9		,%%xmm11		\n\t"\
		"subpd		%%xmm7		,%%xmm4						\n\t		subpd		%%xmm15		,%%xmm12		\n\t"\
		"addq		$0x040		,%%rcx						\n\t		addpd		%%xmm12		,%%xmm8			\n\t"\
		"addq		$0x60		,%%rbx						\n\t		addpd		%%xmm13		,%%xmm9			\n\t"\
		"movaps		     (%%rcx),%%xmm6						\n\t		subpd		%%xmm12		,%%xmm10		\n\t"\
		"movaps		0x010(%%rcx),%%xmm7						\n\t		subpd		%%xmm13		,%%xmm11		\n\t"\
		"addpd		%%xmm4		,%%xmm0						\n\t		movaps		0x080(%%rcx),%%xmm12		\n\t"\
		"addpd		%%xmm5		,%%xmm1						\n\t		movaps		0x090(%%rcx),%%xmm13		\n\t"\
		"subpd		%%xmm4		,%%xmm2						\n\t		movaps		0x080(%%rcx),%%xmm14		\n\t"\
		"subpd		%%xmm5		,%%xmm3						\n\t		movaps		0x090(%%rcx),%%xmm15		\n\t"\
		"movaps		%%xmm6		,%%xmm4						\n\t		mulpd		0x080(%%rbx),%%xmm12		\n\t"\
		"movaps		%%xmm7		,%%xmm5						\n\t		mulpd		0x080(%%rbx),%%xmm13		\n\t"\
		"mulpd		     (%%rbx),%%xmm4						\n\t		mulpd		0x090(%%rbx),%%xmm14		\n\t"\
		"mulpd		     (%%rbx),%%xmm5						\n\t		mulpd		0x090(%%rbx),%%xmm15		\n\t"\
		"mulpd		0x010(%%rbx),%%xmm6						\n\t		addpd		%%xmm14		,%%xmm13		\n\t"\
		"mulpd		0x010(%%rbx),%%xmm7						\n\t		subpd		%%xmm15		,%%xmm12		\n\t"\
		"/*movq		%%rax		,%%rdx		*/				\n\t		movaps		%%xmm13		,0x090(%%rdx)	\n\t"\
		"addpd		%%xmm6		,%%xmm5						\n\t		movaps		%%xmm12		,0x080(%%rdx)	\n\t"\
		"subpd		%%xmm7		,%%xmm4						\n\t	/*	subq	$0x20	,%%rbx */				\n\t"\
		"movaps		%%xmm5		,0x010(%%rdx)	\n\t	addq	$0x040		,%%rax								\n\t"\
		"movaps		%%xmm4		,     (%%rdx)	\n\t	subq	$0x20		,%%rbx								\n\t"\
		"movaps		     (%%rax),%%xmm4						\n\t		movaps		0x080(%%rax),%%xmm12		\n\t"\
		"movaps		0x010(%%rax),%%xmm5						\n\t		movaps		0x090(%%rax),%%xmm13		\n\t"\
		"movaps		%%xmm4		,%%xmm6						\n\t		movaps		0x080(%%rax),%%xmm14		\n\t"\
		"movaps		%%xmm5		,%%xmm7						\n\t		movaps		0x090(%%rax),%%xmm15		\n\t"\
		"mulpd		     (%%rbx),%%xmm4						\n\t		mulpd		0x080(%%rbx),%%xmm12		\n\t"\
		"mulpd		     (%%rbx),%%xmm5						\n\t		mulpd		0x080(%%rbx),%%xmm13		\n\t"\
		"mulpd		0x010(%%rbx),%%xmm6						\n\t		mulpd		0x090(%%rbx),%%xmm14		\n\t"\
		"mulpd		0x010(%%rbx),%%xmm7						\n\t		mulpd		0x090(%%rbx),%%xmm15		\n\t"\
		"addpd		%%xmm6		,%%xmm5						\n\t		addpd		%%xmm14		,%%xmm13		\n\t"\
		"subpd		%%xmm7		,%%xmm4						\n\t		subpd		%%xmm15		,%%xmm12		\n\t"\
		"movaps		%%xmm5		,%%xmm7						\n\t		movaps		%%xmm13		,%%xmm15		\n\t"\
		"movaps		%%xmm4		,%%xmm6						\n\t		movaps		%%xmm12		,%%xmm14		\n\t"\
		"subpd		     (%%rdx),%%xmm4						\n\t		subpd		0x080(%%rdx),%%xmm12		\n\t"\
		"subpd		0x010(%%rdx),%%xmm5						\n\t		subpd		0x090(%%rdx),%%xmm13		\n\t"\
		"addpd		     (%%rdx),%%xmm6						\n\t		addpd		0x080(%%rdx),%%xmm14		\n\t"\
		"addpd		0x010(%%rdx),%%xmm7						\n\t		addpd		0x090(%%rdx),%%xmm15		\n\t"\
		"subpd		%%xmm6		,%%xmm0						\n\t		subpd		%%xmm14		,%%xmm8			\n\t"\
		"subpd		%%xmm5		,%%xmm2						\n\t		subpd		%%xmm15		,%%xmm9			\n\t"\
		"subpd		%%xmm7		,%%xmm1						\n\t	/*	movaps		%%xmm8		,0x0c0(%%rdx)*/	\n\t"\
		"subpd		%%xmm4		,%%xmm3						\n\t		subpd		%%xmm13		,%%xmm10		\n\t"\
		"/*movaps		%%xmm0		,0x040(%%rdx)	*/		\n\t	/*	movaps		%%xmm9		,0x0d0(%%rdx)*/	\n\t"\
		"/*movaps		%%xmm2		,0x020(%%rdx)	*/		\n\t		subpd		%%xmm12		,%%xmm11		\n\t"\
		"/*movaps		%%xmm1		,0x050(%%rdx)	*/		\n\t		mulpd		(%%rdi)		,%%xmm14		\n\t"\
		"/*movaps		%%xmm3		,0x070(%%rdx)	*/		\n\t		mulpd		(%%rdi)		,%%xmm13		\n\t"\
		"mulpd		(%%rdi)		,%%xmm6						\n\t		mulpd		(%%rdi)		,%%xmm15		\n\t"\
		"mulpd		(%%rdi)		,%%xmm5						\n\t		mulpd		(%%rdi)		,%%xmm12		\n\t"\
		"mulpd		(%%rdi)		,%%xmm7						\n\t		addpd		%%xmm8		,%%xmm14		\n\t"\
		"mulpd		(%%rdi)		,%%xmm4						\n\t		addpd		%%xmm10		,%%xmm13		\n\t"\
		"addpd		%%xmm0		,%%xmm6						\n\t		addpd		%%xmm9		,%%xmm15		\n\t"\
		"addpd		%%xmm2		,%%xmm5						\n\t		addpd		%%xmm11		,%%xmm12		\n\t"\
		"addpd		%%xmm1		,%%xmm7						\n\t		movaps		%%xmm14		,0x080(%%rdx)	\n\t"\
		"addpd		%%xmm3		,%%xmm4						\n\t		movaps		%%xmm15		,0x090(%%rdx)	\n\t"\
		"/*movaps		%%xmm6		,     (%%rdx)	*/		\n\t		movaps		%%xmm10		,%%xmm14		\n\t"\
		"/*movaps		%%xmm5		,0x060(%%rdx)	*/		\n\t		movaps		%%xmm13		,%%xmm15		\n\t"\
		"/*movaps		%%xmm7		,0x010(%%rdx)	*/		\n\t		subpd		%%xmm12		,%%xmm10		\n\t"\
		"/*movaps		%%xmm4		,0x030(%%rdx)	*/		\n\t		subpd		%%xmm11		,%%xmm13		\n\t"\
		"													\n\t		addpd		%%xmm12		,%%xmm14		\n\t"\
		"													\n\t		addpd		%%xmm11		,%%xmm15		\n\t"\
		"													\n\t		mulpd		(%%rsi)		,%%xmm10		\n\t"\
		"													\n\t		mulpd		(%%rsi)		,%%xmm13		\n\t"\
		"													\n\t		mulpd		(%%rsi)		,%%xmm14		\n\t"\
		"													\n\t		mulpd		(%%rsi)		,%%xmm15		\n\t"\
		"													\n\t	/*	movaps		%%xmm10		,0x0a0(%%rdx)*/	\n\t"\
		"													\n\t	/*	movaps		%%xmm13		,0x0e0(%%rdx)*/	\n\t"\
		"													\n\t	/*	movaps		%%xmm14		,0x0b0(%%rdx)*/	\n\t"\
		"													\n\t	/*	movaps		%%xmm15		,0x0f0(%%rdx)*/	\n\t"\
		"/*****	SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r20)	*****/\n\t"\
		"																movaps		0x080(%%rdx),%%xmm11		\n\t"\
		"																movaps		0x090(%%rdx),%%xmm12		\n\t"\
		"subpd		%%xmm11		,%%xmm6						\n\t		subpd		%%xmm10		,%%xmm2			\n\t"\
		"subpd		%%xmm9		,%%xmm0						\n\t		subpd		%%xmm15		,%%xmm5			\n\t"\
		"subpd		%%xmm12		,%%xmm7						\n\t		subpd		%%xmm14		,%%xmm4			\n\t"\
		"subpd		%%xmm8		,%%xmm1						\n\t		subpd		%%xmm13		,%%xmm3			\n\t"\
		"mulpd		(%%rdi)		,%%xmm11					\n\t		mulpd		(%%rdi)		,%%xmm10		\n\t"\
		"mulpd		(%%rdi)		,%%xmm9						\n\t		mulpd		(%%rdi)		,%%xmm15		\n\t"\
		"mulpd		(%%rdi)		,%%xmm12					\n\t		mulpd		(%%rdi)		,%%xmm14		\n\t"\
		"mulpd		(%%rdi)		,%%xmm8						\n\t		mulpd		(%%rdi)		,%%xmm13		\n\t"\
		"addpd		%%xmm6		,%%xmm11					\n\t		addpd		%%xmm2		,%%xmm10		\n\t"\
		"addpd		%%xmm0		,%%xmm9						\n\t		addpd		%%xmm5		,%%xmm15		\n\t"\
		"addpd		%%xmm7		,%%xmm12					\n\t		addpd		%%xmm4		,%%xmm14		\n\t"\
		"addpd		%%xmm1		,%%xmm8						\n\t		addpd		%%xmm3		,%%xmm13		\n\t"\
		"movaps		%%xmm6		,0x080(%%rdx)				\n\t		movaps		%%xmm2		,0x0a0(%%rdx)	\n\t"\
		"movaps		%%xmm0		,0x040(%%rdx)				\n\t		movaps		%%xmm5		,0x060(%%rdx)	\n\t"\
		"movaps		%%xmm7		,0x090(%%rdx)				\n\t		movaps		%%xmm4		,0x0b0(%%rdx)	\n\t"\
		"movaps		%%xmm1		,0x0d0(%%rdx)				\n\t		movaps		%%xmm3		,0x0f0(%%rdx)	\n\t"\
		"movaps		%%xmm11		,     (%%rdx)				\n\t		movaps		%%xmm10		,0x020(%%rdx)	\n\t"\
		"movaps		%%xmm9		,0x0c0(%%rdx)				\n\t		movaps		%%xmm15		,0x0e0(%%rdx)	\n\t"\
		"movaps		%%xmm12		,0x010(%%rdx)				\n\t		movaps		%%xmm14		,0x030(%%rdx)	\n\t"\
		"movaps		%%xmm8		,0x050(%%rdx)				\n\t		movaps		%%xmm13		,0x070(%%rdx)	\n\t"\
		"/**********************************************************************************/\n\t"\
		"/*...and now do eight radix-4 transforms, including the internal twiddle factots: */\n\t"\
		"/**********************************************************************************/\n\t"\
		"movq		%[__isrt2]		,%%rsi		\n\t"\
		"/*...Block 1: t00,t10,t20,t30	*/			\n\t	/*...Block 5: t08,t18,t28,t38	*/		\n\t"\
		"movq		%[__r00]		,%%rax			\n\t	movq		$0x080		,%%r10			\n\t"\
		"movq		%[__r10]		,%%rbx			\n\t	movq		$0x080		,%%r11			\n\t"\
		"movq		%[__r20]		,%%rcx			\n\t	movq		$0x080		,%%r12			\n\t"\
		"movq		%[__r30]		,%%rdx			\n\t	movq		$0x080		,%%r13			\n\t"\
		"movaps		      (%%rax)	,%%xmm0			\n\t	addq		%%rax		,%%r10			\n\t"\
		"movaps		 0x010(%%rax)	,%%xmm1			\n\t	addq		%%rbx		,%%r11			\n\t"\
		"movaps		      (%%rbx)	,%%xmm2			\n\t	addq		%%rcx		,%%r12			\n\t"\
		"movaps		 0x010(%%rbx)	,%%xmm3			\n\t	addq		%%rdx		,%%r13			\n\t"\
		"addpd		       %%xmm0	,%%xmm2			\n\t	movaps		      (%%r10)	,%%xmm8		\n\t"\
		"addpd		       %%xmm1	,%%xmm3			\n\t	movaps		 0x010(%%r10)	,%%xmm9		\n\t"\
		"subpd		      (%%rbx)	,%%xmm0			\n\t	movaps		      (%%r11)	,%%xmm10	\n\t"\
		"subpd		 0x010(%%rbx)	,%%xmm1			\n\t	movaps		 0x010(%%r11)	,%%xmm11	\n\t"\
		"movaps		      (%%rcx)	,%%xmm4			\n\t	addpd		       %%xmm9	,%%xmm10	\n\t"\
		"movaps		 0x010(%%rcx)	,%%xmm5			\n\t	addpd		       %%xmm8	,%%xmm11	\n\t"\
		"movaps		      (%%rdx)	,%%xmm6			\n\t	subpd		 0x010(%%r11)	,%%xmm8		\n\t"\
		"movaps		 0x010(%%rdx)	,%%xmm7			\n\t	subpd		      (%%r11)	,%%xmm9		\n\t"\
		"addpd		       %%xmm4	,%%xmm6			\n\t	movaps		      (%%r12)	,%%xmm12	\n\t"\
		"addpd		       %%xmm5	,%%xmm7			\n\t	movaps		 0x010(%%r12)	,%%xmm13	\n\t"\
		"subpd		      (%%rdx)	,%%xmm4			\n\t	movaps		      (%%r13)	,%%xmm14	\n\t"\
		"subpd		 0x010(%%rdx)	,%%xmm5			\n\t	movaps		 0x010(%%r13)	,%%xmm15	\n\t"\
		"subpd		%%xmm6		,%%xmm2				\n\t	subpd		       %%xmm13	,%%xmm12	\n\t"\
		"subpd		%%xmm7		,%%xmm3				\n\t	addpd		      (%%r12)	,%%xmm13	\n\t"\
		"mulpd		(%%rdi)		,%%xmm6				\n\t	addpd		       %%xmm15	,%%xmm14	\n\t"\
		"mulpd		(%%rdi)		,%%xmm7				\n\t	subpd		      (%%r13)	,%%xmm15	\n\t"\
		"movaps		%%xmm2		,      (%%rcx)		\n\t	mulpd		(%%rsi)		,%%xmm12		\n\t"\
		"movaps		%%xmm3		, 0x010(%%rcx)		\n\t	mulpd		(%%rsi)		,%%xmm13		\n\t"\
		"addpd		%%xmm2		,%%xmm6				\n\t	mulpd		(%%rsi)		,%%xmm14		\n\t"\
		"addpd		%%xmm3		,%%xmm7				\n\t	mulpd		(%%rsi)		,%%xmm15		\n\t"\
		"movaps		%%xmm6		,      (%%rax)		\n\t	subpd		%%xmm14		,%%xmm12		\n\t"\
		"movaps		%%xmm7		, 0x010(%%rax)		\n\t	subpd		%%xmm15		,%%xmm13		\n\t"\
		"subpd		%%xmm5		,%%xmm0				\n\t	mulpd		(%%rdi)		,%%xmm14		\n\t"\
		"subpd		%%xmm4		,%%xmm1				\n\t	mulpd		(%%rdi)		,%%xmm15		\n\t"\
		"mulpd		(%%rdi)		,%%xmm5				\n\t	addpd		%%xmm12		,%%xmm14		\n\t"\
		"mulpd		(%%rdi)		,%%xmm4				\n\t	addpd		%%xmm13		,%%xmm15		\n\t"\
		"movaps		%%xmm0		,      (%%rbx)		\n\t	subpd		%%xmm12		,%%xmm8			\n\t"\
		"movaps		%%xmm1		, 0x010(%%rdx)		\n\t	subpd		%%xmm13		,%%xmm10		\n\t"\
		"addpd		%%xmm0		,%%xmm5				\n\t	mulpd		(%%rdi)		,%%xmm12		\n\t"\
		"addpd		%%xmm1		,%%xmm4				\n\t	mulpd		(%%rdi)		,%%xmm13		\n\t"\
		"movaps		%%xmm5		,      (%%rdx)		\n\t	movaps		%%xmm8		,      (%%r12)	\n\t"\
		"movaps		%%xmm4		, 0x010(%%rbx)		\n\t	movaps		%%xmm10		, 0x010(%%r12)	\n\t"\
		"/*...Block 3: t04,t14,t24,t34	*/			\n\t"\
		"addq		$0x040		,%%rax				\n\t	addpd		%%xmm8		,%%xmm12		\n\t"\
		"addq		$0x040		,%%rbx				\n\t	addpd		%%xmm10		,%%xmm13		\n\t"\
		"addq		$0x040		,%%rcx				\n\t	movaps		%%xmm12		,      (%%r10)	\n\t"\
		"addq		$0x040		,%%rdx				\n\t	movaps		%%xmm13		, 0x010(%%r10)	\n\t"\
		"movaps		0x10(%%rsi)	,%%xmm8		/* c */	\n\t	subpd		%%xmm15		,%%xmm11		\n\t"\
		"movaps		0x20(%%rsi)	,%%xmm10	/* s */	\n\t	subpd		%%xmm14		,%%xmm9			\n\t"\
		"movaps		      (%%rcx)	,%%xmm4			\n\t	mulpd		(%%rdi)		,%%xmm15		\n\t"\
		"movaps		      (%%rdx)	,%%xmm6			\n\t	mulpd		(%%rdi)		,%%xmm14		\n\t"\
		"movaps		 0x010(%%rcx)	,%%xmm5			\n\t	movaps		%%xmm11		,      (%%r11)	\n\t"\
		"movaps		 0x010(%%rdx)	,%%xmm7			\n\t	movaps		%%xmm9		, 0x010(%%r13)	\n\t"\
		"movaps		       %%xmm4	,%%xmm0			\n\t	addpd		%%xmm11		,%%xmm15		\n\t"\
		"movaps		       %%xmm6	,%%xmm2			\n\t	addpd		%%xmm9		,%%xmm14		\n\t"\
		"movaps		       %%xmm5	,%%xmm1			\n\t	movaps		%%xmm15		,      (%%r13)	\n\t"\
		"movaps		       %%xmm7	,%%xmm3			\n\t	movaps		%%xmm14		, 0x010(%%r11)	\n\t"\
		"													/*...Block 7: t0C,t1C,t2C,t3C	*/		\n\t"\
		"mulpd		 %%xmm8	,%%xmm4					\n\t	addq		$0x040		,%%r10			\n\t"\
		"mulpd		 %%xmm10,%%xmm6					\n\t	addq		$0x040		,%%r11			\n\t"\
		"mulpd		 %%xmm10,%%xmm1					\n\t	addq		$0x040		,%%r12			\n\t"\
		"mulpd		 %%xmm8	,%%xmm3					\n\t	addq		$0x040		,%%r13			\n\t"\
		"mulpd		 %%xmm8	,%%xmm5					\n\t	movaps		      (%%r12)	,%%xmm12	\n\t"\
		"mulpd		 %%xmm10,%%xmm7					\n\t	movaps		      (%%r13)	,%%xmm14	\n\t"\
		"mulpd		 %%xmm10,%%xmm0					\n\t	movaps		 0x010(%%r12)	,%%xmm13	\n\t"\
		"mulpd		 %%xmm8	,%%xmm2					\n\t	movaps		 0x010(%%r13)	,%%xmm15	\n\t"\
		"subpd		%%xmm1		,%%xmm4				\n\t	movaps		       %%xmm13	,%%xmm9		\n\t"\
		"subpd		%%xmm3		,%%xmm6				\n\t	movaps		       %%xmm15	,%%xmm11	\n\t"\
		"addpd		%%xmm0		,%%xmm5				\n\t	mulpd		 %%xmm10,%%xmm12			\n\t"\
		"addpd		%%xmm2		,%%xmm7				\n\t	mulpd		 %%xmm8	,%%xmm14			\n\t"\
		"subpd		%%xmm6		,%%xmm4				\n\t	mulpd		 %%xmm8	,%%xmm9				\n\t"\
		"subpd		%%xmm7		,%%xmm5				\n\t	mulpd		 %%xmm10,%%xmm11			\n\t"\
		"mulpd		(%%rdi)		,%%xmm6				\n\t	mulpd		 %%xmm10,%%xmm13			\n\t"\
		"mulpd		(%%rdi)		,%%xmm7				\n\t	mulpd		 %%xmm8	,%%xmm15			\n\t"\
		"addpd		%%xmm4		,%%xmm6				\n\t	mulpd		(%%r12)	,%%xmm8				\n\t"\
		"addpd		%%xmm5		,%%xmm7				\n\t	mulpd		(%%r13)	,%%xmm10			\n\t"\
		"movaps		      (%%rbx)	,%%xmm2			\n\t	subpd		%%xmm9		,%%xmm12		\n\t"\
		"movaps		 0x010(%%rbx)	,%%xmm3			\n\t	subpd		%%xmm11		,%%xmm14		\n\t"\
		"subpd		 0x010(%%rbx)	,%%xmm2			\n\t	addpd		%%xmm8		,%%xmm13		\n\t"\
		"addpd		      (%%rbx)	,%%xmm3			\n\t	addpd		%%xmm10		,%%xmm15		\n\t"\
		"mulpd		      (%%rsi)	,%%xmm2			\n\t	subpd		%%xmm14		,%%xmm12		\n\t"\
		"mulpd		      (%%rsi)	,%%xmm3			\n\t	subpd		%%xmm15		,%%xmm13		\n\t"\
		"movaps		      (%%rax)	,%%xmm0			\n\t	mulpd		(%%rdi)		,%%xmm14		\n\t"\
		"movaps		 0x010(%%rax)	,%%xmm1			\n\t	mulpd		(%%rdi)		,%%xmm15		\n\t"\
		"subpd		%%xmm2		,%%xmm0				\n\t	addpd		%%xmm12		,%%xmm14		\n\t"\
		"subpd		%%xmm3		,%%xmm1				\n\t	addpd		%%xmm13		,%%xmm15		\n\t"\
		"addpd		      (%%rax)	,%%xmm2			\n\t	movaps		      (%%r11)	,%%xmm10	\n\t"\
		"addpd		 0x010(%%rax)	,%%xmm3			\n\t	movaps		 0x010(%%r11)	,%%xmm11	\n\t"\
		"subpd		%%xmm6		,%%xmm2				\n\t	addpd		 0x010(%%r11)	,%%xmm10	\n\t"\
		"subpd		%%xmm7		,%%xmm3				\n\t	subpd		      (%%r11)	,%%xmm11	\n\t"\
		"mulpd		(%%rdi)		,%%xmm6				\n\t	mulpd		      (%%rsi)	,%%xmm10	\n\t"\
		"mulpd		(%%rdi)		,%%xmm7				\n\t	mulpd		      (%%rsi)	,%%xmm11	\n\t"\
		"movaps		%%xmm2		,      (%%rcx)		\n\t	movaps		      (%%r10)	,%%xmm8		\n\t"\
		"movaps		%%xmm3		, 0x010(%%rcx)		\n\t	movaps		 0x010(%%r10)	,%%xmm9		\n\t"\
		"addpd		%%xmm2		,%%xmm6				\n\t	subpd		%%xmm10		,%%xmm8			\n\t"\
		"addpd		%%xmm3		,%%xmm7				\n\t	subpd		%%xmm11		,%%xmm9			\n\t"\
		"movaps		%%xmm6		,      (%%rax)		\n\t	addpd		      (%%r10)	,%%xmm10	\n\t"\
		"movaps		%%xmm7		, 0x010(%%rax)		\n\t	addpd		 0x010(%%r10)	,%%xmm11	\n\t"\
		"subpd		%%xmm5		,%%xmm0				\n\t	subpd		%%xmm12		,%%xmm8			\n\t"\
		"subpd		%%xmm4		,%%xmm1				\n\t	subpd		%%xmm13		,%%xmm9			\n\t"\
		"mulpd		(%%rdi)		,%%xmm5				\n\t	mulpd		(%%rdi)		,%%xmm12		\n\t"\
		"mulpd		(%%rdi)		,%%xmm4				\n\t	mulpd		(%%rdi)		,%%xmm13		\n\t"\
		"movaps		%%xmm0		,      (%%rbx)		\n\t	movaps		%%xmm8		,      (%%r12)	\n\t"\
		"movaps		%%xmm1		, 0x010(%%rdx)		\n\t	movaps		%%xmm9		, 0x010(%%r12)	\n\t"\
		"addpd		%%xmm0		,%%xmm5				\n\t	addpd		%%xmm8		,%%xmm12		\n\t"\
		"addpd		%%xmm1		,%%xmm4				\n\t	addpd		%%xmm9		,%%xmm13		\n\t"\
		"movaps		%%xmm5		,      (%%rdx)		\n\t	movaps		%%xmm12		,      (%%r10)	\n\t"\
		"movaps		%%xmm4		, 0x010(%%rbx)		\n\t	movaps		%%xmm13		, 0x010(%%r10)	\n\t"\
		"/*...Block 2: t02,t12,t22,t32	*/\n\t"\
		"subq		$0x020		,%%rax				\n\t	subpd		%%xmm15		,%%xmm10		\n\t"\
		"subq		$0x020		,%%rbx				\n\t	subpd		%%xmm14		,%%xmm11		\n\t"\
		"subq		$0x020		,%%rcx				\n\t	mulpd		(%%rdi)		,%%xmm15		\n\t"\
		"subq		$0x020		,%%rdx				\n\t	mulpd		(%%rdi)		,%%xmm14		\n\t"\
		"addq		$0x30		,%%rsi	/* cc1 */	\n\t	movaps		%%xmm10		,      (%%r11)	\n\t"\
		"movaps		      (%%rcx)	,%%xmm4			\n\t	movaps		%%xmm11		, 0x010(%%r13)	\n\t"\
		"movaps		      (%%rdx)	,%%xmm6			\n\t	addpd		%%xmm10		,%%xmm15		\n\t"\
		"movaps		 0x010(%%rcx)	,%%xmm5			\n\t	addpd		%%xmm11		,%%xmm14		\n\t"\
		"movaps		 0x010(%%rdx)	,%%xmm7			\n\t	movaps		%%xmm15		,      (%%r13)	\n\t"\
		"movaps		      (%%rcx)	,%%xmm0			\n\t	movaps		%%xmm14		, 0x010(%%r11)	\n\t"\
		"													/*...Block 6: t0A,t1A,t2A,t3A	*/\n\t"\
		"movaps		      (%%rdx)	,%%xmm2			\n\t	subq		$0x020		,%%r10			\n\t"\
		"movaps		 0x010(%%rcx)	,%%xmm1			\n\t	subq		$0x020		,%%r11			\n\t"\
		"movaps		 0x010(%%rdx)	,%%xmm3			\n\t	subq		$0x020		,%%r12			\n\t"\
		"mulpd		      (%%rsi)	,%%xmm4			\n\t	subq		$0x020		,%%r13			\n\t"\
		"mulpd		 0x20 (%%rsi)	,%%xmm6			\n\t	movaps		      (%%r12)	,%%xmm12	\n\t"\
		"mulpd		 0x010(%%rsi)	,%%xmm1			\n\t	movaps		      (%%r13)	,%%xmm14	\n\t"\
		"mulpd		 0x30 (%%rsi)	,%%xmm3			\n\t	movaps		 0x010(%%r12)	,%%xmm13	\n\t"\
		"mulpd		      (%%rsi)	,%%xmm5			\n\t	movaps		 0x010(%%r13)	,%%xmm15	\n\t"\
		"mulpd		 0x20 (%%rsi)	,%%xmm7			\n\t	movaps		      (%%r12)	,%%xmm8		\n\t"\
		"mulpd		 0x010(%%rsi)	,%%xmm0			\n\t	movaps		      (%%r13)	,%%xmm10	\n\t"\
		"mulpd		 0x30 (%%rsi)	,%%xmm2			\n\t	movaps		 0x010(%%r12)	,%%xmm9		\n\t"\
		"subpd		%%xmm1		,%%xmm4				\n\t	movaps		 0x010(%%r13)	,%%xmm11	\n\t"\
		"subpd		%%xmm3		,%%xmm6				\n\t	mulpd		 0x30 (%%rsi)	,%%xmm12	\n\t"\
		"addpd		%%xmm0		,%%xmm5				\n\t	mulpd		      (%%rsi)	,%%xmm14	\n\t"\
		"addpd		%%xmm2		,%%xmm7				\n\t	mulpd		 0x20 (%%rsi)	,%%xmm9		\n\t"\
		"subpd		%%xmm6		,%%xmm4				\n\t	mulpd		 0x010(%%rsi)	,%%xmm11	\n\t"\
		"subpd		%%xmm7		,%%xmm5				\n\t	mulpd		 0x30 (%%rsi)	,%%xmm13	\n\t"\
		"mulpd		(%%rdi)		,%%xmm6				\n\t	mulpd		      (%%rsi)	,%%xmm15	\n\t"\
		"mulpd		(%%rdi)		,%%xmm7				\n\t	mulpd		 0x20 (%%rsi)	,%%xmm8		\n\t"\
		"addpd		%%xmm4		,%%xmm6				\n\t	mulpd		 0x010(%%rsi)	,%%xmm10	\n\t"\
		"addpd		%%xmm5		,%%xmm7				\n\t	subpd		%%xmm9		,%%xmm12		\n\t"\
		"movaps		      (%%rbx)	,%%xmm2			\n\t	addpd		%%xmm11		,%%xmm14		\n\t"\
		"movaps		 0x010(%%rbx)	,%%xmm0			\n\t	addpd		%%xmm8		,%%xmm13		\n\t"\
		"movaps		      (%%rbx)	,%%xmm1			\n\t	subpd		%%xmm10		,%%xmm15		\n\t"\
		"movaps		 0x010(%%rbx)	,%%xmm3			\n\t	subpd		%%xmm14		,%%xmm12		\n\t"\
		"mulpd		-0x020(%%rsi)	,%%xmm2			\n\t	subpd		%%xmm15		,%%xmm13		\n\t"\
		"mulpd		-0x010(%%rsi)	,%%xmm0			\n\t	mulpd		(%%rdi)		,%%xmm14		\n\t"\
		"mulpd		-0x020(%%rsi)	,%%xmm3			\n\t	mulpd		(%%rdi)		,%%xmm15		\n\t"\
		"mulpd		-0x010(%%rsi)	,%%xmm1			\n\t	addpd		%%xmm12		,%%xmm14		\n\t"\
		"subpd		%%xmm0		,%%xmm2				\n\t	addpd		%%xmm13		,%%xmm15		\n\t"\
		"addpd		%%xmm1		,%%xmm3				\n\t	movaps		      (%%r11)	,%%xmm10	\n\t"\
		"movaps		      (%%rax)	,%%xmm0			\n\t	movaps		 0x010(%%r11)	,%%xmm8		\n\t"\
		"movaps		 0x010(%%rax)	,%%xmm1			\n\t	movaps		      (%%r11)	,%%xmm9		\n\t"\
		"subpd		%%xmm2		,%%xmm0				\n\t	movaps		 0x010(%%r11)	,%%xmm11	\n\t"\
		"subpd		%%xmm3		,%%xmm1				\n\t	mulpd		-0x010(%%rsi)	,%%xmm10	\n\t"\
		"addpd		      (%%rax)	,%%xmm2			\n\t	mulpd		-0x020(%%rsi)	,%%xmm8		\n\t"\
		"addpd		 0x010(%%rax)	,%%xmm3			\n\t	mulpd		-0x010(%%rsi)	,%%xmm11	\n\t"\
		"subpd		%%xmm6		,%%xmm2				\n\t	mulpd		-0x020(%%rsi)	,%%xmm9		\n\t"\
		"subpd		%%xmm7		,%%xmm3				\n\t	addpd		%%xmm8		,%%xmm10		\n\t"\
		"mulpd		(%%rdi)		,%%xmm6				\n\t	subpd		%%xmm9		,%%xmm11		\n\t"\
		"mulpd		(%%rdi)		,%%xmm7				\n\t	movaps		      (%%r10)	,%%xmm8		\n\t"\
		"movaps		%%xmm2		,      (%%rcx)		\n\t	movaps		 0x010(%%r10)	,%%xmm9		\n\t"\
		"movaps		%%xmm3		, 0x010(%%rcx)		\n\t	subpd		%%xmm10		,%%xmm8			\n\t"\
		"addpd		%%xmm2		,%%xmm6				\n\t	subpd		%%xmm11		,%%xmm9			\n\t"\
		"addpd		%%xmm3		,%%xmm7				\n\t	addpd		      (%%r10)	,%%xmm10	\n\t"\
		"movaps		%%xmm6		,      (%%rax)		\n\t	addpd		 0x010(%%r10)	,%%xmm11	\n\t"\
		"movaps		%%xmm7		, 0x010(%%rax)		\n\t	subpd		%%xmm12		,%%xmm8			\n\t"\
		"subpd		%%xmm5		,%%xmm0				\n\t	subpd		%%xmm13		,%%xmm9			\n\t"\
		"subpd		%%xmm4		,%%xmm1				\n\t	mulpd		(%%rdi)		,%%xmm12		\n\t"\
		"mulpd		(%%rdi)		,%%xmm5				\n\t	mulpd		(%%rdi)		,%%xmm13		\n\t"\
		"mulpd		(%%rdi)		,%%xmm4				\n\t	movaps		%%xmm8		,      (%%r12)	\n\t"\
		"movaps		%%xmm0		,      (%%rbx)		\n\t	movaps		%%xmm9		, 0x010(%%r12)	\n\t"\
		"movaps		%%xmm1		, 0x010(%%rdx)		\n\t	addpd		%%xmm8		,%%xmm12		\n\t"\
		"addpd		%%xmm0		,%%xmm5				\n\t	addpd		%%xmm9		,%%xmm13		\n\t"\
		"addpd		%%xmm1		,%%xmm4				\n\t	movaps		%%xmm12		,      (%%r10)	\n\t"\
		"movaps		%%xmm5		,      (%%rdx)		\n\t	movaps		%%xmm13		, 0x010(%%r10)	\n\t"\
		"movaps		%%xmm4		, 0x010(%%rbx)		\n\t	subpd		%%xmm15		,%%xmm10		\n\t"\
		"/*...Block 4: t06,t16,t26,t36	*/\n\t"\
		"addq		$0x040		,%%rax				\n\t	subpd		%%xmm14		,%%xmm11		\n\t"\
		"addq		$0x040		,%%rbx				\n\t	addpd		%%xmm15		,%%xmm15		\n\t"\
		"addq		$0x040		,%%rcx				\n\t	addpd		%%xmm14		,%%xmm14		\n\t"\
		"addq		$0x040		,%%rdx				\n\t	movaps		%%xmm10		,      (%%r11)	\n\t"\
		"movaps		      (%%rcx)	,%%xmm4			\n\t	movaps		%%xmm11		, 0x010(%%r13)	\n\t"\
		"movaps		      (%%rdx)	,%%xmm6			\n\t	addpd		%%xmm10		,%%xmm15		\n\t"\
		"movaps		 0x010(%%rcx)	,%%xmm5			\n\t	addpd		%%xmm11		,%%xmm14		\n\t"\
		"movaps		 0x010(%%rdx)	,%%xmm7			\n\t	movaps		%%xmm15		,      (%%r13)	\n\t"\
		"movaps		      (%%rcx)	,%%xmm0			\n\t	movaps		%%xmm14		, 0x010(%%r11)	\n\t"\
		"													/*...Block 8: t0E,t1E,t2E,t3E	*/		\n\t"\
		"movaps		      (%%rdx)	,%%xmm2			\n\t	addq		$0x040		,%%r10			\n\t"\
		"movaps		 0x010(%%rcx)	,%%xmm1			\n\t	addq		$0x040		,%%r11			\n\t"\
		"movaps		 0x010(%%rdx)	,%%xmm3			\n\t	addq		$0x040		,%%r12			\n\t"\
		"mulpd		 0x20 (%%rsi)	,%%xmm4			\n\t	addq		$0x040		,%%r13			\n\t"\
		"mulpd		 0x010(%%rsi)	,%%xmm6			\n\t	movaps		      (%%r12)	,%%xmm12	\n\t"\
		"mulpd		 0x30 (%%rsi)	,%%xmm1			\n\t	movaps		      (%%r13)	,%%xmm14	\n\t"\
		"mulpd		      (%%rsi)	,%%xmm3			\n\t	movaps		 0x010(%%r12)	,%%xmm13	\n\t"\
		"mulpd		 0x20 (%%rsi)	,%%xmm5			\n\t	movaps		 0x010(%%r13)	,%%xmm15	\n\t"\
		"mulpd		 0x010(%%rsi)	,%%xmm7			\n\t	movaps		      (%%r12)	,%%xmm8		\n\t"\
		"mulpd		 0x30 (%%rsi)	,%%xmm0			\n\t	movaps		      (%%r13)	,%%xmm10	\n\t"\
		"mulpd		      (%%rsi)	,%%xmm2			\n\t	movaps		 0x010(%%r12)	,%%xmm9		\n\t"\
		"subpd		%%xmm1		,%%xmm4				\n\t	movaps		 0x010(%%r13)	,%%xmm11	\n\t"\
		"addpd		%%xmm3		,%%xmm6				\n\t	mulpd		 0x010(%%rsi)	,%%xmm12	\n\t"\
		"addpd		%%xmm0		,%%xmm5				\n\t	mulpd		 0x30 (%%rsi)	,%%xmm14	\n\t"\
		"subpd		%%xmm2		,%%xmm7				\n\t	mulpd		      (%%rsi)	,%%xmm9		\n\t"\
		"subpd		%%xmm6		,%%xmm4				\n\t	mulpd		 0x20 (%%rsi)	,%%xmm11	\n\t"\
		"subpd		%%xmm7		,%%xmm5				\n\t	mulpd		 0x010(%%rsi)	,%%xmm13	\n\t"\
		"mulpd		(%%rdi)		,%%xmm6				\n\t	mulpd		 0x30 (%%rsi)	,%%xmm15	\n\t"\
		"mulpd		(%%rdi)		,%%xmm7				\n\t	mulpd		      (%%rsi)	,%%xmm8		\n\t"\
		"addpd		%%xmm4		,%%xmm6				\n\t	mulpd		 0x20 (%%rsi)	,%%xmm10	\n\t"\
		"addpd		%%xmm5		,%%xmm7				\n\t	subpd		%%xmm9		,%%xmm12		\n\t"\
		"movaps		      (%%rbx)	,%%xmm2			\n\t	subpd		%%xmm11		,%%xmm14		\n\t"\
		"movaps		 0x010(%%rbx)	,%%xmm0			\n\t	addpd		%%xmm8		,%%xmm13		\n\t"\
		"movaps		      (%%rbx)	,%%xmm1			\n\t	addpd		%%xmm10		,%%xmm15		\n\t"\
		"movaps		 0x010(%%rbx)	,%%xmm3			\n\t	subpd		%%xmm14		,%%xmm12		\n\t"\
		"mulpd		-0x010(%%rsi)	,%%xmm2			\n\t	subpd		%%xmm15		,%%xmm13		\n\t"\
		"mulpd		-0x020(%%rsi)	,%%xmm0			\n\t	mulpd		(%%rdi)		,%%xmm14		\n\t"\
		"mulpd		-0x010(%%rsi)	,%%xmm3			\n\t	mulpd		(%%rdi)		,%%xmm15		\n\t"\
		"mulpd		-0x020(%%rsi)	,%%xmm1			\n\t	addpd		%%xmm12		,%%xmm14		\n\t"\
		"subpd		%%xmm0		,%%xmm2				\n\t	addpd		%%xmm13		,%%xmm15		\n\t"\
		"addpd		%%xmm1		,%%xmm3				\n\t	movaps		      (%%r11)	,%%xmm10	\n\t"\
		"movaps		      (%%rax)	,%%xmm0			\n\t	movaps		 0x010(%%r11)	,%%xmm8		\n\t"\
		"movaps		 0x010(%%rax)	,%%xmm1			\n\t	movaps		      (%%r11)	,%%xmm9		\n\t"\
		"subpd		%%xmm2		,%%xmm0				\n\t	movaps		 0x010(%%r11)	,%%xmm11	\n\t"\
		"subpd		%%xmm3		,%%xmm1				\n\t	mulpd		-0x020(%%rsi)	,%%xmm10	\n\t"\
		"addpd		      (%%rax)	,%%xmm2			\n\t	mulpd		-0x010(%%rsi)	,%%xmm8		\n\t"\
		"addpd		 0x010(%%rax)	,%%xmm3			\n\t	mulpd		-0x020(%%rsi)	,%%xmm11	\n\t"\
		"subpd		%%xmm4		,%%xmm2				\n\t	mulpd		-0x010(%%rsi)	,%%xmm9		\n\t"\
		"subpd		%%xmm5		,%%xmm3				\n\t	addpd		%%xmm8		,%%xmm10		\n\t"\
		"mulpd		(%%rdi)		,%%xmm4				\n\t	subpd		%%xmm9		,%%xmm11		\n\t"\
		"mulpd		(%%rdi)		,%%xmm5				\n\t	movaps		      (%%r10)	,%%xmm8		\n\t"\
		"movaps		%%xmm2		,      (%%rcx)		\n\t	movaps		 0x010(%%r10)	,%%xmm9		\n\t"\
		"movaps		%%xmm3		, 0x010(%%rcx)		\n\t	subpd		%%xmm10		,%%xmm8			\n\t"\
		"addpd		%%xmm2		,%%xmm4				\n\t	subpd		%%xmm11		,%%xmm9			\n\t"\
		"addpd		%%xmm3		,%%xmm5				\n\t	addpd		      (%%r10)	,%%xmm10	\n\t"\
		"movaps		%%xmm4		,      (%%rax)		\n\t	addpd		 0x010(%%r10)	,%%xmm11	\n\t"\
		"movaps		%%xmm5		, 0x010(%%rax)		\n\t	subpd		%%xmm12		,%%xmm8			\n\t"\
		"subpd		%%xmm7		,%%xmm0				\n\t	subpd		%%xmm13		,%%xmm9			\n\t"\
		"subpd		%%xmm6		,%%xmm1				\n\t	mulpd		(%%rdi)		,%%xmm12		\n\t"\
		"mulpd		(%%rdi)		,%%xmm7				\n\t	mulpd		(%%rdi)		,%%xmm13		\n\t"\
		"mulpd		(%%rdi)		,%%xmm6				\n\t	movaps		%%xmm8		,      (%%r12)	\n\t"\
		"movaps		%%xmm0		,      (%%rbx)		\n\t	movaps		%%xmm9		, 0x010(%%r12)	\n\t"\
		"movaps		%%xmm1		, 0x010(%%rdx)		\n\t	addpd		%%xmm8		,%%xmm12		\n\t"\
		"addpd		%%xmm0		,%%xmm7				\n\t	addpd		%%xmm9		,%%xmm13		\n\t"\
		"addpd		%%xmm1		,%%xmm6				\n\t	movaps		%%xmm12		,      (%%r10)	\n\t"\
		"movaps		%%xmm7		,      (%%rdx)		\n\t	movaps		%%xmm13		, 0x010(%%r10)	\n\t"\
		"movaps		%%xmm6		, 0x010(%%rbx)		\n\t	subpd		%%xmm15		,%%xmm10		\n\t"\
		"													subpd		%%xmm14		,%%xmm11		\n\t"\
		"													mulpd		(%%rdi)		,%%xmm15		\n\t"\
		"													mulpd		(%%rdi)		,%%xmm14		\n\t"\
		"													movaps		%%xmm10		,      (%%r11)	\n\t"\
		"													movaps		%%xmm11		, 0x010(%%r13)	\n\t"\
		"													addpd		%%xmm10		,%%xmm15		\n\t"\
		"													addpd		%%xmm11		,%%xmm14		\n\t"\
		"													movaps		%%xmm15		,      (%%r13)	\n\t"\
		"													movaps		%%xmm14		, 0x010(%%r11)	\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__add1] "m" (Xadd1)\
		 ,[__r00] "m" (Xr00)\
		 ,[__r10] "m" (Xr10)\
		 ,[__r20] "m" (Xr20)\
		 ,[__r30] "m" (Xr30)\
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__cc0] "m" (Xcc0)\
		 ,[__c00] "m" (Xc00)\
		 ,[__c01] "m" (Xc01)\
		 ,[__c02] "m" (Xc02)\
		 ,[__c03] "m" (Xc03)\
		 ,[__c05] "m" (Xc05)\
		 ,[__c07] "m" (Xc07)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r10","r11","r12","r13","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"		/* Clobbered registers */\
	);\
	}

	#define FOO(Xadd0,Xadd1,Xr00,Xr10,Xr20,Xr30,Xisrt2,Xcc0,Xc00,Xc01,Xc02,Xc03,Xc05,Xc07)\
	{\
	__asm__ volatile (\
		"movq	%[__r00] ,%%rcx	\n\t"\
		"/***************************************************************************************/\n\t"\
		"/* AVX DEBUG: Excerpt just the permute/local-store parts of forward DIF radix-32 pass: */\n\t"\
		"/***************************************************************************************/\n\t"\
		"/*...Block 0:	*/									\n\t		movq	%[__isrt2],%%rsi		\n\t"\
		"movq		%[__add0]	,%%rax						\n\t		movq		%%rsi,%%rdi			\n\t"\
		"movq		%[__add1]	,%%rbx						\n\t		addq	$0x480,%%rdi	/* two */	\n\t"\
		"/*****	SSE2_RADIX4_DIF_4WRAPPER(c00,c08,c10,c18,r00) ****	SSE2_RADIX4_DIF_4WRAPPER_2NDOFTWO(c04,c0C,c14,c1C,r08)	*****/\n\t"\
		"movq		%[__r00]	,%%rcx						\n\t		/*addq		$0x080,%%rcx // __r08 */	\n\t"\
		"movq		%[__c00]	,%%rdx						\n\t		/*addq		$0x080,%%rdx // __c04 */	\n\t"\
		"movaps		     (%%rax),%%xmm6						\n\t		movaps		0x40 (%%rax),%%xmm14		\n\t"\
		"movaps		     (%%rax),%%xmm0						\n\t		movaps		0x40 (%%rax),%%xmm8			\n\t"\
		"unpckhpd	     (%%rbx),%%xmm6						\n\t		unpckhpd	0x40 (%%rbx),%%xmm14		\n\t"\
		"unpcklpd	     (%%rbx),%%xmm0						\n\t		unpcklpd	0x40 (%%rbx),%%xmm8			\n\t"\
		"movaps		%%xmm6		,0x200(%%rcx)				\n\t		movaps		%%xmm14,0x280(%%rcx)		\n\t"\
		"movaps		%%xmm0		,     (%%rcx)				\n\t		movaps		%%xmm8 ,0x080(%%rcx)		\n\t"\
		"movaps		0x010(%%rax),%%xmm7						\n\t		movaps		0x50 (%%rax),%%xmm15		\n\t"\
		"movaps		0x010(%%rax),%%xmm1						\n\t		movaps		0x50 (%%rax),%%xmm9			\n\t"\
		"unpckhpd	0x010(%%rbx),%%xmm7						\n\t		unpckhpd	0x50 (%%rbx),%%xmm15		\n\t"\
		"unpcklpd	0x010(%%rbx),%%xmm1						\n\t		unpcklpd	0x50 (%%rbx),%%xmm9			\n\t"\
		"movaps		%%xmm7		,0x210(%%rcx)				\n\t		movaps		%%xmm15,0x290(%%rcx)		\n\t"\
		"movaps		%%xmm1		,0x010(%%rcx)				\n\t		movaps		%%xmm9 ,0x090(%%rcx)		\n\t"\
		"addq		$0x020		,%%rdx																			\n\t"\
		"movaps		0x100(%%rax),%%xmm6						\n\t		movaps		0x140(%%rax),%%xmm14		\n\t"\
		"movaps		0x100(%%rax),%%xmm4						\n\t		movaps		0x140(%%rax),%%xmm12		\n\t"\
		"unpckhpd	0x100(%%rbx),%%xmm6						\n\t		unpckhpd	0x140(%%rbx),%%xmm14		\n\t"\
		"unpcklpd	0x100(%%rbx),%%xmm4						\n\t		unpcklpd	0x140(%%rbx),%%xmm12		\n\t"\
		"movaps		%%xmm6		,0x220(%%rcx)				\n\t		movaps		%%xmm14,0x2a0(%%rcx)		\n\t"\
		"movaps		%%xmm4		,0x020(%%rcx)				\n\t		movaps		%%xmm12,0x0a0(%%rcx)		\n\t"\
		"movaps		0x110(%%rax),%%xmm7						\n\t		movaps		0x150(%%rax),%%xmm15		\n\t"\
		"movaps		0x110(%%rax),%%xmm5						\n\t		movaps		0x150(%%rax),%%xmm13		\n\t"\
		"unpckhpd	0x110(%%rbx),%%xmm7						\n\t		unpckhpd	0x150(%%rbx),%%xmm15		\n\t"\
		"unpcklpd	0x110(%%rbx),%%xmm5						\n\t		unpcklpd	0x150(%%rbx),%%xmm13		\n\t"\
		"movaps		%%xmm7		,0x230(%%rcx)				\n\t		movaps		%%xmm15,0x2b0(%%rcx)		\n\t"\
		"movaps		%%xmm5		,0x030(%%rcx)				\n\t		movaps		%%xmm13,0x0b0(%%rcx)		\n\t"\
		"addq		$0x040		,%%rdx																			\n\t"\
		"movaps		0x180(%%rax),%%xmm6						\n\t		movaps		0x1c0(%%rax),%%xmm14		\n\t"\
		"movaps		0x180(%%rax),%%xmm4						\n\t		movaps		0x1c0(%%rax),%%xmm12		\n\t"\
		"unpckhpd	0x180(%%rbx),%%xmm6						\n\t		unpckhpd	0x1c0(%%rbx),%%xmm14		\n\t"\
		"unpcklpd	0x180(%%rbx),%%xmm4						\n\t		unpcklpd	0x1c0(%%rbx),%%xmm12		\n\t"\
		"movaps		%%xmm6		,0x260(%%rcx)				\n\t		movaps		%%xmm14,0x2e0(%%rcx)		\n\t"\
		"movaps		%%xmm4		,0x060(%%rcx)				\n\t		movaps		%%xmm12,0x0e0(%%rcx)		\n\t"\
		"movaps		0x190(%%rax),%%xmm7						\n\t		movaps		0x1d0(%%rax),%%xmm15		\n\t"\
		"movaps		0x190(%%rax),%%xmm5						\n\t		movaps		0x1d0(%%rax),%%xmm13		\n\t"\
		"unpckhpd	0x190(%%rbx),%%xmm7						\n\t		unpckhpd	0x1d0(%%rbx),%%xmm15		\n\t"\
		"unpcklpd	0x190(%%rbx),%%xmm5						\n\t		unpcklpd	0x1d0(%%rbx),%%xmm13		\n\t"\
		"movaps		%%xmm7		,0x270(%%rcx)				\n\t		movaps		%%xmm15,0x2f0(%%rcx)		\n\t"\
		"movaps		%%xmm5		,0x070(%%rcx)				\n\t		movaps		%%xmm13,0x0f0(%%rcx)		\n\t"\
		"subq		$0x020		,%%rdx																			\n\t"\
		"movaps		0x080(%%rax),%%xmm6						\n\t		movaps		0x0c0(%%rax),%%xmm14		\n\t"\
		"movaps		0x080(%%rax),%%xmm4						\n\t		movaps		0x0c0(%%rax),%%xmm12		\n\t"\
		"unpckhpd	0x080(%%rbx),%%xmm6						\n\t		unpckhpd	0x0c0(%%rbx),%%xmm14		\n\t"\
		"unpcklpd	0x080(%%rbx),%%xmm4						\n\t		unpcklpd	0x0c0(%%rbx),%%xmm12		\n\t"\
		"movaps		%%xmm6		,0x240(%%rcx)				\n\t		movaps		%%xmm14,0x2c0(%%rcx)		\n\t"\
		"movaps		%%xmm4		,0x040(%%rcx)				\n\t		movaps		%%xmm12,0x0c0(%%rcx)		\n\t"\
		"movaps		0x090(%%rax),%%xmm7						\n\t		movaps		0x0d0(%%rax),%%xmm15		\n\t"\
		"movaps		0x090(%%rax),%%xmm5						\n\t		movaps		0x0d0(%%rax),%%xmm13		\n\t"\
		"unpckhpd	0x090(%%rbx),%%xmm7						\n\t		unpckhpd	0x0d0(%%rbx),%%xmm15		\n\t"\
		"unpcklpd	0x090(%%rbx),%%xmm5						\n\t		unpcklpd	0x0d0(%%rbx),%%xmm13		\n\t"\
		"movaps		%%xmm7		,0x250(%%rcx)				\n\t		movaps		%%xmm15,0x2d0(%%rcx)		\n\t"\
		"movaps		%%xmm5		,0x050(%%rcx)				\n\t		movaps		%%xmm13,0x0d0(%%rcx)		\n\t"\
		"\n\t"\
		"/*...Block 2:	*/\n\t"\
		"addq		$0x20		,%%rax\n\t"\
		"addq		$0x20		,%%rbx\n\t"\
		"/*****	SSE2_RADIX4_DIF_4WRAPPER(c02,c0A,c12,c1A,r10) ****	SSE2_RADIX4_DIF_4WRAPPER_2NDOFTWO(c06,c0E,c16,c1E,r18)	*****/\n\t"\
		"movq		%[__r10]	,%%rcx						\n\t		/*addq		$0x080,%%rcx // __r18 */	\n\t"\
		"movq		%[__c02]	,%%rdx						\n\t		/*addq		$0x080,%%rdx // __c06 */	\n\t"\
		"movaps		     (%%rax),%%xmm6						\n\t		movaps		0x40 (%%rax),%%xmm14		\n\t"\
		"movaps		     (%%rax),%%xmm0						\n\t		movaps		0x40 (%%rax),%%xmm8			\n\t"\
		"unpckhpd	     (%%rbx),%%xmm6						\n\t		unpckhpd	0x40 (%%rbx),%%xmm14		\n\t"\
		"unpcklpd	     (%%rbx),%%xmm0						\n\t		unpcklpd	0x40 (%%rbx),%%xmm8			\n\t"\
		"movaps		%%xmm6		,0x200(%%rcx)				\n\t		movaps		%%xmm14,0x280(%%rcx)		\n\t"\
		"movaps		%%xmm0		,0x000(%%rcx)				\n\t		movaps		%%xmm8 ,0x080(%%rcx)		\n\t"\
		"movaps		0x010(%%rax),%%xmm7						\n\t		movaps		0x50 (%%rax),%%xmm15		\n\t"\
		"movaps		0x010(%%rax),%%xmm1						\n\t		movaps		0x50 (%%rax),%%xmm9			\n\t"\
		"unpckhpd	0x010(%%rbx),%%xmm7						\n\t		unpckhpd	0x50 (%%rbx),%%xmm15		\n\t"\
		"unpcklpd	0x010(%%rbx),%%xmm1						\n\t		unpcklpd	0x50 (%%rbx),%%xmm9			\n\t"\
		"movaps		%%xmm7		,0x210(%%rcx)				\n\t		movaps		%%xmm15,0x290(%%rcx)		\n\t"\
		"movaps		%%xmm1		,0x010(%%rcx)				\n\t		movaps		%%xmm9 ,0x090(%%rcx)		\n\t"\
		"addq		$0x020		,%%rdx																			\n\t"\
		"movaps		0x100(%%rax),%%xmm6						\n\t		movaps		0x140(%%rax),%%xmm14		\n\t"\
		"movaps		0x100(%%rax),%%xmm4						\n\t		movaps		0x140(%%rax),%%xmm12		\n\t"\
		"unpckhpd	0x100(%%rbx),%%xmm6						\n\t		unpckhpd	0x140(%%rbx),%%xmm14		\n\t"\
		"unpcklpd	0x100(%%rbx),%%xmm4						\n\t		unpcklpd	0x140(%%rbx),%%xmm12		\n\t"\
		"movaps		%%xmm6		,0x220(%%rcx)				\n\t		movaps		%%xmm14,0x2a0(%%rcx)		\n\t"\
		"movaps		%%xmm4		,0x020(%%rcx)				\n\t		movaps		%%xmm12,0x0a0(%%rcx)		\n\t"\
		"movaps		0x110(%%rax),%%xmm7						\n\t		movaps		0x150(%%rax),%%xmm15		\n\t"\
		"movaps		0x110(%%rax),%%xmm5						\n\t		movaps		0x150(%%rax),%%xmm13		\n\t"\
		"unpckhpd	0x110(%%rbx),%%xmm7						\n\t		unpckhpd	0x150(%%rbx),%%xmm15		\n\t"\
		"unpcklpd	0x110(%%rbx),%%xmm5						\n\t		unpcklpd	0x150(%%rbx),%%xmm13		\n\t"\
		"movaps		%%xmm7		,0x230(%%rcx)				\n\t		movaps		%%xmm15,0x2b0(%%rcx)		\n\t"\
		"movaps		%%xmm5		,0x030(%%rcx)				\n\t		movaps		%%xmm13,0x0b0(%%rcx)		\n\t"\
		"addq		$0x040		,%%rdx																			\n\t"\
		"movaps		0x180(%%rax),%%xmm6						\n\t		movaps		0x1c0(%%rax),%%xmm14		\n\t"\
		"movaps		0x180(%%rax),%%xmm4						\n\t		movaps		0x1c0(%%rax),%%xmm12		\n\t"\
		"unpckhpd	0x180(%%rbx),%%xmm6						\n\t		unpckhpd	0x1c0(%%rbx),%%xmm14		\n\t"\
		"unpcklpd	0x180(%%rbx),%%xmm4						\n\t		unpcklpd	0x1c0(%%rbx),%%xmm12		\n\t"\
		"movaps		%%xmm6		,0x260(%%rcx)				\n\t		movaps		%%xmm14,0x2e0(%%rcx)		\n\t"\
		"movaps		%%xmm4		,0x060(%%rcx)				\n\t		movaps		%%xmm12,0x0e0(%%rcx)		\n\t"\
		"movaps		0x190(%%rax),%%xmm7						\n\t		movaps		0x1d0(%%rax),%%xmm15		\n\t"\
		"movaps		0x190(%%rax),%%xmm5						\n\t		movaps		0x1d0(%%rax),%%xmm13		\n\t"\
		"unpckhpd	0x190(%%rbx),%%xmm7						\n\t		unpckhpd	0x1d0(%%rbx),%%xmm15		\n\t"\
		"unpcklpd	0x190(%%rbx),%%xmm5						\n\t		unpcklpd	0x1d0(%%rbx),%%xmm13		\n\t"\
		"movaps		%%xmm7		,0x270(%%rcx)				\n\t		movaps		%%xmm15,0x2f0(%%rcx)		\n\t"\
		"movaps		%%xmm5		,0x070(%%rcx)				\n\t		movaps		%%xmm13,0x0f0(%%rcx)		\n\t"\
		"subq		$0x020		,%%rdx																			\n\t"\
		"movaps		0x080(%%rax),%%xmm6						\n\t		movaps		0x0c0(%%rax),%%xmm14		\n\t"\
		"movaps		0x080(%%rax),%%xmm4						\n\t		movaps		0x0c0(%%rax),%%xmm12		\n\t"\
		"unpckhpd	0x080(%%rbx),%%xmm6						\n\t		unpckhpd	0x0c0(%%rbx),%%xmm14		\n\t"\
		"unpcklpd	0x080(%%rbx),%%xmm4						\n\t		unpcklpd	0x0c0(%%rbx),%%xmm12		\n\t"\
		"movaps		%%xmm6		,0x240(%%rcx)				\n\t		movaps		%%xmm14,0x2c0(%%rcx)		\n\t"\
		"movaps		%%xmm4		,0x040(%%rcx)				\n\t		movaps		%%xmm12,0x0c0(%%rcx)		\n\t"\
		"movaps		0x090(%%rax),%%xmm7						\n\t		movaps		0x0d0(%%rax),%%xmm15		\n\t"\
		"movaps		0x090(%%rax),%%xmm5						\n\t		movaps		0x0d0(%%rax),%%xmm13		\n\t"\
		"unpckhpd	0x090(%%rbx),%%xmm7						\n\t		unpckhpd	0x0d0(%%rbx),%%xmm15		\n\t"\
		"unpcklpd	0x090(%%rbx),%%xmm5						\n\t		unpcklpd	0x0d0(%%rbx),%%xmm13		\n\t"\
		"movaps		%%xmm7		,0x250(%%rcx)				\n\t		movaps		%%xmm15,0x2d0(%%rcx)		\n\t"\
		"movaps		%%xmm5		,0x050(%%rcx)				\n\t		movaps		%%xmm13,0x0d0(%%rcx)		\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__add1] "m" (Xadd1)\
		 ,[__r00] "m" (Xr00)\
		 ,[__r10] "m" (Xr10)\
		 ,[__r20] "m" (Xr20)\
		 ,[__r30] "m" (Xr30)\
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__cc0] "m" (Xcc0)\
		 ,[__c00] "m" (Xc00)\
		 ,[__c01] "m" (Xc01)\
		 ,[__c02] "m" (Xc02)\
		 ,[__c03] "m" (Xc03)\
		 ,[__c05] "m" (Xc05)\
		 ,[__c07] "m" (Xc07)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r10","r11","r12","r13","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"		/* Clobbered registers */\
	);\
	}

	#define SSE2_RADIX32_WRAPPER_DIT(Xadd0,Xadd1,Xisrt2,Xr00,Xr08,Xr10,Xr20,Xr28,Xr30,Xc01,Xc02,Xc04,Xc06,Xc08,Xc0A,Xc0C,Xc0E,Xc10,Xc12,Xc14,Xc16,Xc18,Xc1A,Xc1C,Xc1E)\
	{\
	__asm__ volatile (\
		"/************************************************************************/\n\t"\
		"/*...And do an inverse DIT radix-32 pass on the squared-data blocks.	*/\n\t"\
		"/************************************************************************/\n\t"\
		"/*...Block 1: */\n\t"\
		"movq		%[__isrt2]	,%%rsi	\n\t"\
		"movq		%[__r00]	,%%rax	\n\t"\
		"leaq		0x200(%%rax),%%rbx	\n\t"\
		"leaq		0x100(%%rax),%%rcx	\n\t"\
		"leaq		0x300(%%rax),%%rdx	\n\t"\
		"/*****	SSE2_RADIX4_DIT_IN_PLACE()	*****/	\n\t		/*...Block 2 has tmp-addresses offset +0x40 w.r.to Block 1:	*/\n\t"\
		"movaps		      (%%rax)	,%%xmm0			\n\t		movaps		 0x040(%%rax)	,%%xmm8 \n\t"\
		"movaps		 0x010(%%rax)	,%%xmm1			\n\t		movaps		 0x050(%%rax)	,%%xmm9 \n\t"\
		"movaps		      (%%rax)	,%%xmm2			\n\t		movaps		 0x040(%%rax)	,%%xmm10\n\t"\
		"movaps		 0x010(%%rax)	,%%xmm3			\n\t		movaps		 0x050(%%rax)	,%%xmm11\n\t"\
		"addpd		      (%%rbx)	,%%xmm0			\n\t		addpd		 0x040(%%rbx)	,%%xmm8 \n\t"\
		"addpd		 0x010(%%rbx)	,%%xmm1			\n\t		addpd		 0x050(%%rbx)	,%%xmm9 \n\t"\
		"subpd		      (%%rbx)	,%%xmm2			\n\t		subpd		 0x040(%%rbx)	,%%xmm10\n\t"\
		"subpd		 0x010(%%rbx)	,%%xmm3			\n\t		subpd		 0x050(%%rbx)	,%%xmm11\n\t"\
		"movaps		      (%%rcx)	,%%xmm4			\n\t		movaps		 0x040(%%rcx)	,%%xmm12\n\t"\
		"movaps		 0x010(%%rcx)	,%%xmm5			\n\t		movaps		 0x050(%%rcx)	,%%xmm13\n\t"\
		"movaps		      (%%rcx)	,%%xmm6			\n\t		movaps		 0x040(%%rcx)	,%%xmm14\n\t"\
		"movaps		 0x010(%%rcx)	,%%xmm7			\n\t		movaps		 0x050(%%rcx)	,%%xmm15\n\t"\
		"addpd		      (%%rdx)	,%%xmm4			\n\t		addpd		 0x040(%%rdx)	,%%xmm12\n\t"\
		"addpd		 0x010(%%rdx)	,%%xmm5			\n\t		addpd		 0x050(%%rdx)	,%%xmm13\n\t"\
		"subpd		      (%%rdx)	,%%xmm6			\n\t		subpd		 0x040(%%rdx)	,%%xmm14\n\t"\
		"subpd		 0x010(%%rdx)	,%%xmm7			\n\t		subpd		 0x050(%%rdx)	,%%xmm15\n\t"\
		"subpd		%%xmm4		,%%xmm0				\n\t		subpd		%%xmm12		,%%xmm8 \n\t"\
		"subpd		%%xmm5		,%%xmm1				\n\t		subpd		%%xmm13		,%%xmm9 \n\t"\
		"movaps		%%xmm0		,      (%%rbx)		\n\t		movaps		%%xmm8 		, 0x040(%%rbx)\n\t"\
		"movaps		%%xmm1		, 0x010(%%rbx)		\n\t		movaps		%%xmm9 		, 0x050(%%rbx)\n\t"\
		"addpd		%%xmm4		,%%xmm4				\n\t		addpd		%%xmm12		,%%xmm12\n\t"\
		"addpd		%%xmm5		,%%xmm5				\n\t		addpd		%%xmm13		,%%xmm13\n\t"\
		"addpd		%%xmm0		,%%xmm4				\n\t		addpd		%%xmm8 		,%%xmm12\n\t"\
		"addpd		%%xmm1		,%%xmm5				\n\t		addpd		%%xmm9 		,%%xmm13\n\t"\
		"movaps		%%xmm4		,      (%%rax)		\n\t		movaps		%%xmm12		, 0x040(%%rax)\n\t"\
		"movaps		%%xmm5		, 0x010(%%rax)		\n\t		movaps		%%xmm13		, 0x050(%%rax)\n\t"\
		"subpd		%%xmm7		,%%xmm2				\n\t		subpd		%%xmm15		,%%xmm10\n\t"\
		"subpd		%%xmm6		,%%xmm3				\n\t		subpd		%%xmm14		,%%xmm11\n\t"\
		"movaps		%%xmm2		,      (%%rdx)		\n\t		movaps		%%xmm10		, 0x040(%%rdx)\n\t"\
		"movaps		%%xmm3		, 0x010(%%rcx)		\n\t		movaps		%%xmm11		, 0x050(%%rcx)\n\t"\
		"addpd		%%xmm7		,%%xmm7				\n\t		addpd		%%xmm15		,%%xmm15\n\t"\
		"addpd		%%xmm6		,%%xmm6				\n\t		addpd		%%xmm14		,%%xmm14\n\t"\
		"addpd		%%xmm2		,%%xmm7				\n\t		addpd		%%xmm10		,%%xmm15\n\t"\
		"addpd		%%xmm3		,%%xmm6				\n\t		addpd		%%xmm11		,%%xmm14\n\t"\
		"movaps		%%xmm7		,      (%%rcx)		\n\t		movaps		%%xmm15		, 0x040(%%rcx)\n\t"\
		"movaps		%%xmm6		, 0x010(%%rdx)		\n\t		movaps		%%xmm14		, 0x050(%%rdx)\n\t"\
		"addq		$0x080		,%%rax				\n\t"\
		"addq		$0x080		,%%rbx				\n\t"\
		"addq		$0x080		,%%rcx				\n\t"\
		"addq		$0x080		,%%rdx				\n\t"\
		"movaps		      (%%rax)	,%%xmm0			\n\t		movaps		 0x040(%%rax)	,%%xmm8 \n\t"\
		"movaps		 0x010(%%rax)	,%%xmm1			\n\t		movaps		 0x050(%%rax)	,%%xmm9 \n\t"\
		"movaps		      (%%rax)	,%%xmm2			\n\t		movaps		 0x040(%%rax)	,%%xmm10\n\t"\
		"movaps		 0x010(%%rax)	,%%xmm3			\n\t		movaps		 0x050(%%rax)	,%%xmm11\n\t"\
		"addpd		      (%%rbx)	,%%xmm0			\n\t		addpd		 0x040(%%rbx)	,%%xmm8 \n\t"\
		"addpd		 0x010(%%rbx)	,%%xmm1			\n\t		addpd		 0x050(%%rbx)	,%%xmm9 \n\t"\
		"subpd		      (%%rbx)	,%%xmm2			\n\t		subpd		 0x040(%%rbx)	,%%xmm10\n\t"\
		"subpd		 0x010(%%rbx)	,%%xmm3			\n\t		subpd		 0x050(%%rbx)	,%%xmm11\n\t"\
		"movaps		      (%%rcx)	,%%xmm4			\n\t		movaps		 0x040(%%rcx)	,%%xmm12\n\t"\
		"movaps		 0x010(%%rcx)	,%%xmm5			\n\t		movaps		 0x050(%%rcx)	,%%xmm13\n\t"\
		"movaps		      (%%rcx)	,%%xmm6			\n\t		movaps		 0x040(%%rcx)	,%%xmm14\n\t"\
		"movaps		 0x010(%%rcx)	,%%xmm7			\n\t		movaps		 0x050(%%rcx)	,%%xmm15\n\t"\
		"addpd		      (%%rdx)	,%%xmm4			\n\t		addpd		 0x040(%%rdx)	,%%xmm12\n\t"\
		"addpd		 0x010(%%rdx)	,%%xmm5			\n\t		addpd		 0x050(%%rdx)	,%%xmm13\n\t"\
		"subpd		      (%%rdx)	,%%xmm6			\n\t		subpd		 0x040(%%rdx)	,%%xmm14\n\t"\
		"subpd		 0x010(%%rdx)	,%%xmm7			\n\t		subpd		 0x050(%%rdx)	,%%xmm15\n\t"\
		"subpd		%%xmm4		,%%xmm0				\n\t		subpd		%%xmm12		,%%xmm8 \n\t"\
		"subpd		%%xmm5		,%%xmm1				\n\t		subpd		%%xmm13		,%%xmm9 \n\t"\
		"movaps		%%xmm0		,      (%%rbx)		\n\t		movaps		%%xmm8 		, 0x040(%%rbx)\n\t"\
		"movaps		%%xmm1		, 0x010(%%rbx)		\n\t		movaps		%%xmm9 		, 0x050(%%rbx)\n\t"\
		"addpd		%%xmm4		,%%xmm4				\n\t		addpd		%%xmm12		,%%xmm12\n\t"\
		"addpd		%%xmm5		,%%xmm5				\n\t		addpd		%%xmm13		,%%xmm13\n\t"\
		"addpd		%%xmm0		,%%xmm4				\n\t		addpd		%%xmm8 		,%%xmm12\n\t"\
		"addpd		%%xmm1		,%%xmm5				\n\t		addpd		%%xmm9 		,%%xmm13\n\t"\
		"movaps		%%xmm4		,      (%%rax)		\n\t		movaps		%%xmm12		, 0x040(%%rax)\n\t"\
		"movaps		%%xmm5		, 0x010(%%rax)		\n\t		movaps		%%xmm13		, 0x050(%%rax)\n\t"\
		"subpd		%%xmm7		,%%xmm2				\n\t		subpd		%%xmm15		,%%xmm10\n\t"\
		"subpd		%%xmm6		,%%xmm3				\n\t		subpd		%%xmm14		,%%xmm11\n\t"\
		"addpd		%%xmm7		,%%xmm7				\n\t		addpd		%%xmm15		,%%xmm15\n\t"\
		"addpd		%%xmm6		,%%xmm6				\n\t		addpd		%%xmm14		,%%xmm14\n\t"\
		"addpd		%%xmm2		,%%xmm7				\n\t		addpd		%%xmm10		,%%xmm15\n\t"\
		"addpd		%%xmm3		,%%xmm6				\n\t		addpd		%%xmm11		,%%xmm14\n\t"\
		"movaps		%%xmm3		,%%xmm0				\n\t		movaps		%%xmm11		,%%xmm8 \n\t"\
		"movaps		%%xmm6		,%%xmm1				\n\t		movaps		%%xmm14		,%%xmm9 \n\t"\
		"subpd		%%xmm7		,%%xmm3				\n\t		subpd		%%xmm15		,%%xmm11\n\t"\
		"subpd		%%xmm2		,%%xmm6				\n\t		subpd		%%xmm10		,%%xmm14\n\t"\
		"addpd		%%xmm7		,%%xmm0				\n\t		addpd		%%xmm15		,%%xmm8 \n\t"\
		"addpd		%%xmm2		,%%xmm1				\n\t		addpd		%%xmm10		,%%xmm9 \n\t"\
		"mulpd		      (%%rsi)	,%%xmm3			\n\t		mulpd		      (%%rsi)	,%%xmm11\n\t"\
		"mulpd		      (%%rsi)	,%%xmm6			\n\t		mulpd		      (%%rsi)	,%%xmm14\n\t"\
		"mulpd		      (%%rsi)	,%%xmm0			\n\t		mulpd		      (%%rsi)	,%%xmm8 \n\t"\
		"mulpd		      (%%rsi)	,%%xmm1			\n\t		mulpd		      (%%rsi)	,%%xmm9 \n\t"\
		"movaps		%%xmm3		, 0x010(%%rcx)		\n\t		movaps		%%xmm11		, 0x050(%%rcx)\n\t"\
		"movaps		%%xmm6		, 0x010(%%rdx)		\n\t		movaps		%%xmm14		, 0x050(%%rdx)\n\t"\
		"movaps		%%xmm0		,      (%%rcx)		\n\t		movaps		%%xmm8 		, 0x040(%%rcx)\n\t"\
		"movaps		%%xmm1		,      (%%rdx)		\n\t		movaps		%%xmm9 		, 0x040(%%rdx)\n\t"\
		"/***************************	SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS ************************************/\n\t"\
		"movaps		-0x080(%%rax)	,%%xmm0			\n\t		movaps		-0x040(%%rax)	,%%xmm8 \n\t"\
		"movaps		-0x080(%%rbx)	,%%xmm4			\n\t		movaps		-0x040(%%rbx)	,%%xmm12\n\t"\
		"movaps		-0x070(%%rax)	,%%xmm1			\n\t		movaps		-0x030(%%rax)	,%%xmm9 \n\t"\
		"movaps		-0x070(%%rbx)	,%%xmm5			\n\t		movaps		-0x030(%%rbx)	,%%xmm13\n\t"\
		"movaps		      (%%rax)	,%%xmm2			\n\t		movaps		 0x040(%%rax)	,%%xmm10\n\t"\
		"movaps		 0x010(%%rbx)	,%%xmm7			\n\t		movaps		 0x050(%%rbx)	,%%xmm15\n\t"\
		"movaps		 0x010(%%rax)	,%%xmm3			\n\t		movaps		 0x050(%%rax)	,%%xmm11\n\t"\
		"movaps		      (%%rbx)	,%%xmm6			\n\t		movaps		 0x040(%%rbx)	,%%xmm14\n\t"\
		"subpd		%%xmm2		,%%xmm0				\n\t		subpd		%%xmm10		,%%xmm8 \n\t"\
		"subpd		%%xmm7		,%%xmm4				\n\t		subpd		%%xmm15		,%%xmm12\n\t"\
		"subpd		%%xmm3		,%%xmm1				\n\t		subpd		%%xmm11		,%%xmm9 \n\t"\
		"subpd		%%xmm6		,%%xmm5				\n\t		subpd		%%xmm14		,%%xmm13\n\t"\
		"addpd		%%xmm2		,%%xmm2				\n\t		addpd		%%xmm10		,%%xmm10\n\t"\
		"addpd		%%xmm7		,%%xmm7				\n\t		addpd		%%xmm15		,%%xmm15\n\t"\
		"addpd		%%xmm3		,%%xmm3				\n\t		addpd		%%xmm11		,%%xmm11\n\t"\
		"addpd		%%xmm6		,%%xmm6				\n\t		addpd		%%xmm14		,%%xmm14\n\t"\
		"addpd		%%xmm0		,%%xmm2				\n\t		addpd		%%xmm8 		,%%xmm10\n\t"\
		"addpd		%%xmm4		,%%xmm7				\n\t		addpd		%%xmm12		,%%xmm15\n\t"\
		"addpd		%%xmm1		,%%xmm3				\n\t		addpd		%%xmm9 		,%%xmm11\n\t"\
		"addpd		%%xmm5		,%%xmm6				\n\t		addpd		%%xmm13		,%%xmm14\n\t"\
		"movaps		%%xmm0		,      (%%rax)		\n\t		movaps		%%xmm8 		, 0x040(%%rax)\n\t"\
		"movaps		%%xmm4		,      (%%rbx)		\n\t		movaps		%%xmm12		, 0x040(%%rbx)\n\t"\
		"movaps		%%xmm1		, 0x010(%%rax)		\n\t		movaps		%%xmm9 		, 0x050(%%rax)\n\t"\
		"movaps		%%xmm5		,-0x070(%%rbx)		\n\t		movaps		%%xmm13		,-0x030(%%rbx)\n\t"\
		"movaps		%%xmm2		,-0x080(%%rax)		\n\t		movaps		%%xmm10		,-0x040(%%rax)\n\t"\
		"movaps		%%xmm7		,-0x080(%%rbx)		\n\t		movaps		%%xmm15		,-0x040(%%rbx)\n\t"\
		"movaps		%%xmm3		,-0x070(%%rax)		\n\t		movaps		%%xmm11		,-0x030(%%rax)\n\t"\
		"movaps		%%xmm6		, 0x010(%%rbx)		\n\t		movaps		%%xmm14		, 0x050(%%rbx)\n\t"\
		"movaps		-0x080(%%rcx)	,%%xmm0			\n\t		movaps		-0x040(%%rcx)	,%%xmm8 \n\t"\
		"movaps		-0x080(%%rdx)	,%%xmm4			\n\t		movaps		-0x040(%%rdx)	,%%xmm12\n\t"\
		"movaps		-0x070(%%rcx)	,%%xmm1			\n\t		movaps		-0x030(%%rcx)	,%%xmm9 \n\t"\
		"movaps		-0x070(%%rdx)	,%%xmm5			\n\t		movaps		-0x030(%%rdx)	,%%xmm13\n\t"\
		"movaps		      (%%rcx)	,%%xmm2			\n\t		movaps		 0x040(%%rcx)	,%%xmm10\n\t"\
		"movaps		 0x010(%%rdx)	,%%xmm7			\n\t		movaps		 0x050(%%rdx)	,%%xmm15\n\t"\
		"movaps		 0x010(%%rcx)	,%%xmm3			\n\t		movaps		 0x050(%%rcx)	,%%xmm11\n\t"\
		"movaps		      (%%rdx)	,%%xmm6			\n\t		movaps		 0x040(%%rdx)	,%%xmm14\n\t"\
		"subpd		%%xmm2		,%%xmm0				\n\t		subpd		%%xmm10		,%%xmm8 \n\t"\
		"subpd		%%xmm7		,%%xmm4				\n\t		subpd		%%xmm15		,%%xmm12\n\t"\
		"subpd		%%xmm3		,%%xmm1				\n\t		subpd		%%xmm11		,%%xmm9 \n\t"\
		"subpd		%%xmm6		,%%xmm5				\n\t		subpd		%%xmm14		,%%xmm13\n\t"\
		"addpd		%%xmm2		,%%xmm2				\n\t		addpd		%%xmm10		,%%xmm10\n\t"\
		"addpd		%%xmm7		,%%xmm7				\n\t		addpd		%%xmm15		,%%xmm15\n\t"\
		"addpd		%%xmm3		,%%xmm3				\n\t		addpd		%%xmm11		,%%xmm11\n\t"\
		"addpd		%%xmm6		,%%xmm6				\n\t		addpd		%%xmm14		,%%xmm14\n\t"\
		"addpd		%%xmm0		,%%xmm2				\n\t		addpd		%%xmm8 		,%%xmm10\n\t"\
		"addpd		%%xmm4		,%%xmm7				\n\t		addpd		%%xmm12		,%%xmm15\n\t"\
		"addpd		%%xmm1		,%%xmm3				\n\t		addpd		%%xmm9 		,%%xmm11\n\t"\
		"addpd		%%xmm5		,%%xmm6				\n\t		addpd		%%xmm13		,%%xmm14\n\t"\
		"movaps		%%xmm0		,      (%%rcx)		\n\t		movaps		%%xmm8 		, 0x040(%%rcx)\n\t"\
		"movaps		%%xmm4		,      (%%rdx)		\n\t		movaps		%%xmm12		, 0x040(%%rdx)\n\t"\
		"movaps		%%xmm1		, 0x010(%%rcx)		\n\t		movaps		%%xmm9 		, 0x050(%%rcx)\n\t"\
		"movaps		%%xmm5		,-0x070(%%rdx)		\n\t		movaps		%%xmm13		,-0x030(%%rdx)\n\t"\
		"movaps		%%xmm2		,-0x080(%%rcx)		\n\t		movaps		%%xmm10		,-0x040(%%rcx)\n\t"\
		"movaps		%%xmm7		,-0x080(%%rdx)		\n\t		movaps		%%xmm15		,-0x040(%%rdx)\n\t"\
		"movaps		%%xmm3		,-0x070(%%rcx)		\n\t		movaps		%%xmm11		,-0x030(%%rcx)\n\t"\
		"movaps		%%xmm6		, 0x010(%%rdx)		\n\t		movaps		%%xmm14		, 0x050(%%rdx)\n\t"\
		"\n\t"\
		"/*...Blocks 3,4 have tmp-addresses offset +0x40 w.r.to Blocks 1,2, respectively (thus +0x100-0x0c0 = +0x040: */\n\t"\
		"subq		$0x60		,%%rax				\n\t"\
		"subq		$0x60		,%%rbx				\n\t"\
		"subq		$0x60		,%%rcx				\n\t"\
		"subq		$0x60		,%%rdx				\n\t"\
		"movaps		      (%%rax)	,%%xmm0			\n\t		movaps		 0x040(%%rax)	,%%xmm8 \n\t"\
		"movaps		 0x010(%%rax)	,%%xmm1			\n\t		movaps		 0x050(%%rax)	,%%xmm9 \n\t"\
		"movaps		      (%%rax)	,%%xmm2			\n\t		movaps		 0x040(%%rax)	,%%xmm10\n\t"\
		"movaps		 0x010(%%rax)	,%%xmm3			\n\t		movaps		 0x050(%%rax)	,%%xmm11\n\t"\
		"addpd		      (%%rbx)	,%%xmm0			\n\t		addpd		 0x040(%%rbx)	,%%xmm8 \n\t"\
		"addpd		 0x010(%%rbx)	,%%xmm1			\n\t		addpd		 0x050(%%rbx)	,%%xmm9 \n\t"\
		"subpd		      (%%rbx)	,%%xmm2			\n\t		subpd		 0x040(%%rbx)	,%%xmm10\n\t"\
		"subpd		 0x010(%%rbx)	,%%xmm3			\n\t		subpd		 0x050(%%rbx)	,%%xmm11\n\t"\
		"movaps		      (%%rcx)	,%%xmm4			\n\t		movaps		 0x040(%%rcx)	,%%xmm12\n\t"\
		"movaps		 0x010(%%rcx)	,%%xmm5			\n\t		movaps		 0x050(%%rcx)	,%%xmm13\n\t"\
		"movaps		      (%%rcx)	,%%xmm6			\n\t		movaps		 0x040(%%rcx)	,%%xmm14\n\t"\
		"movaps		 0x010(%%rcx)	,%%xmm7			\n\t		movaps		 0x050(%%rcx)	,%%xmm15\n\t"\
		"addpd		      (%%rdx)	,%%xmm4			\n\t		addpd		 0x040(%%rdx)	,%%xmm12\n\t"\
		"addpd		 0x010(%%rdx)	,%%xmm5			\n\t		addpd		 0x050(%%rdx)	,%%xmm13\n\t"\
		"subpd		      (%%rdx)	,%%xmm6			\n\t		subpd		 0x040(%%rdx)	,%%xmm14\n\t"\
		"subpd		 0x010(%%rdx)	,%%xmm7			\n\t		subpd		 0x050(%%rdx)	,%%xmm15\n\t"\
		"subpd		%%xmm4		,%%xmm0				\n\t		subpd		%%xmm12		,%%xmm8 \n\t"\
		"subpd		%%xmm5		,%%xmm1				\n\t		subpd		%%xmm13		,%%xmm9 \n\t"\
		"movaps		%%xmm0		,      (%%rbx)		\n\t		movaps		%%xmm8 		, 0x040(%%rbx)\n\t"\
		"movaps		%%xmm1		, 0x010(%%rbx)		\n\t		movaps		%%xmm9 		, 0x050(%%rbx)\n\t"\
		"addpd		%%xmm4		,%%xmm4				\n\t		addpd		%%xmm12		,%%xmm12\n\t"\
		"addpd		%%xmm5		,%%xmm5				\n\t		addpd		%%xmm13		,%%xmm13\n\t"\
		"addpd		%%xmm0		,%%xmm4				\n\t		addpd		%%xmm8 		,%%xmm12\n\t"\
		"addpd		%%xmm1		,%%xmm5				\n\t		addpd		%%xmm9 		,%%xmm13\n\t"\
		"movaps		%%xmm4		,      (%%rax)		\n\t		movaps		%%xmm12		, 0x040(%%rax)\n\t"\
		"movaps		%%xmm5		, 0x010(%%rax)		\n\t		movaps		%%xmm13		, 0x050(%%rax)\n\t"\
		"subpd		%%xmm7		,%%xmm2				\n\t		subpd		%%xmm15		,%%xmm10\n\t"\
		"subpd		%%xmm6		,%%xmm3				\n\t		subpd		%%xmm14		,%%xmm11\n\t"\
		"movaps		%%xmm2		,      (%%rdx)		\n\t		movaps		%%xmm10		, 0x040(%%rdx)\n\t"\
		"movaps		%%xmm3		, 0x010(%%rcx)		\n\t		movaps		%%xmm11		, 0x050(%%rcx)\n\t"\
		"addpd		%%xmm7		,%%xmm7				\n\t		addpd		%%xmm15		,%%xmm15\n\t"\
		"addpd		%%xmm6		,%%xmm6				\n\t		addpd		%%xmm14		,%%xmm14\n\t"\
		"addpd		%%xmm2		,%%xmm7				\n\t		addpd		%%xmm10		,%%xmm15\n\t"\
		"addpd		%%xmm3		,%%xmm6				\n\t		addpd		%%xmm11		,%%xmm14\n\t"\
		"movaps		%%xmm7		,      (%%rcx)		\n\t		movaps		%%xmm15		, 0x040(%%rcx)\n\t"\
		"movaps		%%xmm6		, 0x010(%%rdx)		\n\t		movaps		%%xmm14		, 0x050(%%rdx)\n\t"\
		"addq		$0x080		,%%rax				\n\t"\
		"addq		$0x080		,%%rbx				\n\t"\
		"addq		$0x080		,%%rcx				\n\t"\
		"addq		$0x080		,%%rdx				\n\t"\
		"movaps		      (%%rax)	,%%xmm0			\n\t		movaps		 0x040(%%rax)	,%%xmm8 \n\t"\
		"movaps		 0x010(%%rax)	,%%xmm1			\n\t		movaps		 0x050(%%rax)	,%%xmm9 \n\t"\
		"movaps		      (%%rax)	,%%xmm2			\n\t		movaps		 0x040(%%rax)	,%%xmm10\n\t"\
		"movaps		 0x010(%%rax)	,%%xmm3			\n\t		movaps		 0x050(%%rax)	,%%xmm11\n\t"\
		"addpd		      (%%rbx)	,%%xmm0			\n\t		addpd		 0x040(%%rbx)	,%%xmm8 \n\t"\
		"addpd		 0x010(%%rbx)	,%%xmm1			\n\t		addpd		 0x050(%%rbx)	,%%xmm9 \n\t"\
		"subpd		      (%%rbx)	,%%xmm2			\n\t		subpd		 0x040(%%rbx)	,%%xmm10\n\t"\
		"subpd		 0x010(%%rbx)	,%%xmm3			\n\t		subpd		 0x050(%%rbx)	,%%xmm11\n\t"\
		"movaps		      (%%rcx)	,%%xmm4			\n\t		movaps		 0x040(%%rcx)	,%%xmm12\n\t"\
		"movaps		 0x010(%%rcx)	,%%xmm5			\n\t		movaps		 0x050(%%rcx)	,%%xmm13\n\t"\
		"movaps		      (%%rcx)	,%%xmm6			\n\t		movaps		 0x040(%%rcx)	,%%xmm14\n\t"\
		"movaps		 0x010(%%rcx)	,%%xmm7			\n\t		movaps		 0x050(%%rcx)	,%%xmm15\n\t"\
		"addpd		      (%%rdx)	,%%xmm4			\n\t		addpd		 0x040(%%rdx)	,%%xmm12\n\t"\
		"addpd		 0x010(%%rdx)	,%%xmm5			\n\t		addpd		 0x050(%%rdx)	,%%xmm13\n\t"\
		"subpd		      (%%rdx)	,%%xmm6			\n\t		subpd		 0x040(%%rdx)	,%%xmm14\n\t"\
		"subpd		 0x010(%%rdx)	,%%xmm7			\n\t		subpd		 0x050(%%rdx)	,%%xmm15\n\t"\
		"subpd		%%xmm4		,%%xmm0				\n\t		subpd		%%xmm12		,%%xmm8 \n\t"\
		"subpd		%%xmm5		,%%xmm1				\n\t		subpd		%%xmm13		,%%xmm9 \n\t"\
		"movaps		%%xmm0		,      (%%rbx)		\n\t		movaps		%%xmm8 		, 0x040(%%rbx)\n\t"\
		"movaps		%%xmm1		, 0x010(%%rbx)		\n\t		movaps		%%xmm9 		, 0x050(%%rbx)\n\t"\
		"addpd		%%xmm4		,%%xmm4				\n\t		addpd		%%xmm12		,%%xmm12\n\t"\
		"addpd		%%xmm5		,%%xmm5				\n\t		addpd		%%xmm13		,%%xmm13\n\t"\
		"addpd		%%xmm0		,%%xmm4				\n\t		addpd		%%xmm8 		,%%xmm12\n\t"\
		"addpd		%%xmm1		,%%xmm5				\n\t		addpd		%%xmm9 		,%%xmm13\n\t"\
		"movaps		%%xmm4		,      (%%rax)		\n\t		movaps		%%xmm12		, 0x040(%%rax)\n\t"\
		"movaps		%%xmm5		, 0x010(%%rax)		\n\t		movaps		%%xmm13		, 0x050(%%rax)\n\t"\
		"subpd		%%xmm7		,%%xmm2				\n\t		subpd		%%xmm15		,%%xmm10\n\t"\
		"subpd		%%xmm6		,%%xmm3				\n\t		subpd		%%xmm14		,%%xmm11\n\t"\
		"addpd		%%xmm7		,%%xmm7				\n\t		addpd		%%xmm15		,%%xmm15\n\t"\
		"addpd		%%xmm6		,%%xmm6				\n\t		addpd		%%xmm14		,%%xmm14\n\t"\
		"addpd		%%xmm2		,%%xmm7				\n\t		addpd		%%xmm10		,%%xmm15\n\t"\
		"addpd		%%xmm3		,%%xmm6				\n\t		addpd		%%xmm11		,%%xmm14\n\t"\
		"movaps		%%xmm3		,%%xmm0				\n\t		movaps		%%xmm11		,%%xmm8 \n\t"\
		"movaps		%%xmm6		,%%xmm1				\n\t		movaps		%%xmm14		,%%xmm9 \n\t"\
		"subpd		%%xmm7		,%%xmm3				\n\t		subpd		%%xmm15		,%%xmm11\n\t"\
		"subpd		%%xmm2		,%%xmm6				\n\t		subpd		%%xmm10		,%%xmm14\n\t"\
		"addpd		%%xmm7		,%%xmm0				\n\t		addpd		%%xmm15		,%%xmm8 \n\t"\
		"addpd		%%xmm2		,%%xmm1				\n\t		addpd		%%xmm10		,%%xmm9 \n\t"\
		"mulpd		      (%%rsi)	,%%xmm3			\n\t		mulpd		      (%%rsi)	,%%xmm11\n\t"\
		"mulpd		      (%%rsi)	,%%xmm6			\n\t		mulpd		      (%%rsi)	,%%xmm14\n\t"\
		"mulpd		      (%%rsi)	,%%xmm0			\n\t		mulpd		      (%%rsi)	,%%xmm8 \n\t"\
		"mulpd		      (%%rsi)	,%%xmm1			\n\t		mulpd		      (%%rsi)	,%%xmm9 \n\t"\
		"movaps		%%xmm3		, 0x010(%%rcx)		\n\t		movaps		%%xmm11		, 0x050(%%rcx)\n\t"\
		"movaps		%%xmm6		, 0x010(%%rdx)		\n\t		movaps		%%xmm14		, 0x050(%%rdx)\n\t"\
		"movaps		%%xmm0		,      (%%rcx)		\n\t		movaps		%%xmm8 		, 0x040(%%rcx)\n\t"\
		"movaps		%%xmm1		,      (%%rdx)		\n\t		movaps		%%xmm9 		, 0x040(%%rdx)\n\t"\
		"/***************************	SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS ************************************/\n\t"\
		"movaps		-0x080(%%rax)	,%%xmm0			\n\t		movaps		-0x040(%%rax)	,%%xmm8 \n\t"\
		"movaps		-0x080(%%rbx)	,%%xmm4			\n\t		movaps		-0x040(%%rbx)	,%%xmm12\n\t"\
		"movaps		-0x070(%%rax)	,%%xmm1			\n\t		movaps		-0x030(%%rax)	,%%xmm9 \n\t"\
		"movaps		-0x070(%%rbx)	,%%xmm5			\n\t		movaps		-0x030(%%rbx)	,%%xmm13\n\t"\
		"movaps		      (%%rax)	,%%xmm2			\n\t		movaps		 0x040(%%rax)	,%%xmm10\n\t"\
		"movaps		 0x010(%%rbx)	,%%xmm7			\n\t		movaps		 0x050(%%rbx)	,%%xmm15\n\t"\
		"movaps		 0x010(%%rax)	,%%xmm3			\n\t		movaps		 0x050(%%rax)	,%%xmm11\n\t"\
		"movaps		      (%%rbx)	,%%xmm6			\n\t		movaps		 0x040(%%rbx)	,%%xmm14\n\t"\
		"subpd		%%xmm2		,%%xmm0				\n\t		subpd		%%xmm10		,%%xmm8 \n\t"\
		"subpd		%%xmm7		,%%xmm4				\n\t		subpd		%%xmm15		,%%xmm12\n\t"\
		"subpd		%%xmm3		,%%xmm1				\n\t		subpd		%%xmm11		,%%xmm9 \n\t"\
		"subpd		%%xmm6		,%%xmm5				\n\t		subpd		%%xmm14		,%%xmm13\n\t"\
		"addpd		%%xmm2		,%%xmm2				\n\t		addpd		%%xmm10		,%%xmm10\n\t"\
		"addpd		%%xmm7		,%%xmm7				\n\t		addpd		%%xmm15		,%%xmm15\n\t"\
		"addpd		%%xmm3		,%%xmm3				\n\t		addpd		%%xmm11		,%%xmm11\n\t"\
		"addpd		%%xmm6		,%%xmm6				\n\t		addpd		%%xmm14		,%%xmm14\n\t"\
		"addpd		%%xmm0		,%%xmm2				\n\t		addpd		%%xmm8 		,%%xmm10\n\t"\
		"addpd		%%xmm4		,%%xmm7				\n\t		addpd		%%xmm12		,%%xmm15\n\t"\
		"addpd		%%xmm1		,%%xmm3				\n\t		addpd		%%xmm9 		,%%xmm11\n\t"\
		"addpd		%%xmm5		,%%xmm6				\n\t		addpd		%%xmm13		,%%xmm14\n\t"\
		"movaps		%%xmm0		,      (%%rax)		\n\t		movaps		%%xmm8 		, 0x040(%%rax)\n\t"\
		"movaps		%%xmm4		,      (%%rbx)		\n\t		movaps		%%xmm12		, 0x040(%%rbx)\n\t"\
		"movaps		%%xmm1		, 0x010(%%rax)		\n\t		movaps		%%xmm9 		, 0x050(%%rax)\n\t"\
		"movaps		%%xmm5		,-0x070(%%rbx)		\n\t		movaps		%%xmm13		,-0x030(%%rbx)\n\t"\
		"movaps		%%xmm2		,-0x080(%%rax)		\n\t		movaps		%%xmm10		,-0x040(%%rax)\n\t"\
		"movaps		%%xmm7		,-0x080(%%rbx)		\n\t		movaps		%%xmm15		,-0x040(%%rbx)\n\t"\
		"movaps		%%xmm3		,-0x070(%%rax)		\n\t		movaps		%%xmm11		,-0x030(%%rax)\n\t"\
		"movaps		%%xmm6		, 0x010(%%rbx)		\n\t		movaps		%%xmm14		, 0x050(%%rbx)\n\t"\
		"movaps		-0x080(%%rcx)	,%%xmm0			\n\t		movaps		-0x040(%%rcx)	,%%xmm8 \n\t"\
		"movaps		-0x080(%%rdx)	,%%xmm4			\n\t		movaps		-0x040(%%rdx)	,%%xmm12\n\t"\
		"movaps		-0x070(%%rcx)	,%%xmm1			\n\t		movaps		-0x030(%%rcx)	,%%xmm9 \n\t"\
		"movaps		-0x070(%%rdx)	,%%xmm5			\n\t		movaps		-0x030(%%rdx)	,%%xmm13\n\t"\
		"movaps		      (%%rcx)	,%%xmm2			\n\t		movaps		 0x040(%%rcx)	,%%xmm10\n\t"\
		"movaps		 0x010(%%rdx)	,%%xmm7			\n\t		movaps		 0x050(%%rdx)	,%%xmm15\n\t"\
		"movaps		 0x010(%%rcx)	,%%xmm3			\n\t		movaps		 0x050(%%rcx)	,%%xmm11\n\t"\
		"movaps		      (%%rdx)	,%%xmm6			\n\t		movaps		 0x040(%%rdx)	,%%xmm14\n\t"\
		"subpd		%%xmm2		,%%xmm0				\n\t		subpd		%%xmm10		,%%xmm8 \n\t"\
		"subpd		%%xmm7		,%%xmm4				\n\t		subpd		%%xmm15		,%%xmm12\n\t"\
		"subpd		%%xmm3		,%%xmm1				\n\t		subpd		%%xmm11		,%%xmm9 \n\t"\
		"subpd		%%xmm6		,%%xmm5				\n\t		subpd		%%xmm14		,%%xmm13\n\t"\
		"addpd		%%xmm2		,%%xmm2				\n\t		addpd		%%xmm10		,%%xmm10\n\t"\
		"addpd		%%xmm7		,%%xmm7				\n\t		addpd		%%xmm15		,%%xmm15\n\t"\
		"addpd		%%xmm3		,%%xmm3				\n\t		addpd		%%xmm11		,%%xmm11\n\t"\
		"addpd		%%xmm6		,%%xmm6				\n\t		addpd		%%xmm14		,%%xmm14\n\t"\
		"addpd		%%xmm0		,%%xmm2				\n\t		addpd		%%xmm8 		,%%xmm10\n\t"\
		"addpd		%%xmm4		,%%xmm7				\n\t		addpd		%%xmm12		,%%xmm15\n\t"\
		"addpd		%%xmm1		,%%xmm3				\n\t		addpd		%%xmm9 		,%%xmm11\n\t"\
		"addpd		%%xmm5		,%%xmm6				\n\t		addpd		%%xmm13		,%%xmm14\n\t"\
		"movaps		%%xmm0		,      (%%rcx)		\n\t		movaps		%%xmm8 		, 0x040(%%rcx)\n\t"\
		"movaps		%%xmm4		,      (%%rdx)		\n\t		movaps		%%xmm12		, 0x040(%%rdx)\n\t"\
		"movaps		%%xmm1		, 0x010(%%rcx)		\n\t		movaps		%%xmm9 		, 0x050(%%rcx)\n\t"\
		"movaps		%%xmm5		,-0x070(%%rdx)		\n\t		movaps		%%xmm13		,-0x030(%%rdx)\n\t"\
		"movaps		%%xmm2		,-0x080(%%rcx)		\n\t		movaps		%%xmm10		,-0x040(%%rcx)\n\t"\
		"movaps		%%xmm7		,-0x080(%%rdx)		\n\t		movaps		%%xmm15		,-0x040(%%rdx)\n\t"\
		"movaps		%%xmm3		,-0x070(%%rcx)		\n\t		movaps		%%xmm11		,-0x030(%%rcx)\n\t"\
		"movaps		%%xmm6		, 0x010(%%rdx)		\n\t		movaps		%%xmm14		, 0x050(%%rdx)\n\t"\
		"/***************************************************************************************/\n\t"\
		"/* Now do eight more radix-4 transforms, including the internal and external twiddles: */\n\t"\
		"/***************************************************************************************/\n\t"\
		"/***********************************************/	\n\t	/************************************************/	\n\t"\
		"/*...Block 2: t02,t12,t22,t32 -> r10,14,12,16: */	\n\t	/*...Block 6: t0A,t1A,t2A,t3A -> r18,1C,1A,1E:  */	\n\t"\
		"/***********************************************/	\n\t	/************************************************/	\n\t"\
		"movq		%[__isrt2]		,%%rsi	\n\t"\
		"movq		%[__r10]		,%%rax	/* base-addr in rcol = c05/r18, so rax/rdi offset +0x80 vs lcol */\n\t"\
		"movq		%[__add1]		,%%rbx	\n\t"\
		"movq		%%rsi			,%%rcx	\n\t"\
		"movq		%%rsi			,%%rdx	\n\t"\
		"movq		%[__c01]		,%%rdi	\n\t"\
		"addq		$0x010			,%%rcx	/* cc0 */	\n\t"\
		"addq		$0x030			,%%rdx	/* cc1 */	\n\t"\
		"movaps		 0x020(%%rax)	,%%xmm4				\n\t		movaps		 0x0a0(%%rax)	,%%xmm12	\n\t"\
		"movaps		 0x060(%%rax)	,%%xmm0				\n\t		movaps		 0x0e0(%%rax)	,%%xmm8 	\n\t"\
		"movaps		 0x030(%%rax)	,%%xmm5				\n\t		movaps		 0x0b0(%%rax)	,%%xmm13	\n\t"\
		"movaps		 0x070(%%rax)	,%%xmm1				\n\t		movaps		 0x0f0(%%rax)	,%%xmm9 	\n\t"\
		"movaps		 0x020(%%rax)	,%%xmm6				\n\t		movaps		 0x0a0(%%rax)	,%%xmm14	\n\t"\
		"movaps		 0x060(%%rax)	,%%xmm2				\n\t		movaps		 0x0e0(%%rax)	,%%xmm10	\n\t"\
		"movaps		 0x030(%%rax)	,%%xmm7				\n\t		movaps		 0x0b0(%%rax)	,%%xmm15	\n\t"\
		"movaps		 0x070(%%rax)	,%%xmm3				\n\t		movaps		 0x0f0(%%rax)	,%%xmm11	\n\t"\
		"mulpd		      (%%rdx)	,%%xmm4				\n\t		mulpd		 0x30 (%%rdx)	,%%xmm12	\n\t"\
		"mulpd		 0x20 (%%rdx)	,%%xmm0				\n\t		mulpd		      (%%rdx)	,%%xmm8 	\n\t"\
		"mulpd		      (%%rdx)	,%%xmm5				\n\t		mulpd		 0x30 (%%rdx)	,%%xmm13	\n\t"\
		"mulpd		 0x20 (%%rdx)	,%%xmm1				\n\t		mulpd		      (%%rdx)	,%%xmm9 	\n\t"\
		"mulpd		 0x010(%%rdx)	,%%xmm6				\n\t		mulpd		 0x20 (%%rdx)	,%%xmm14	\n\t"\
		"mulpd		 0x30 (%%rdx)	,%%xmm2				\n\t		mulpd		 0x010(%%rdx)	,%%xmm10	\n\t"\
		"mulpd		 0x010(%%rdx)	,%%xmm7				\n\t		mulpd		 0x20 (%%rdx)	,%%xmm15	\n\t"\
		"mulpd		 0x30 (%%rdx)	,%%xmm3				\n\t		mulpd		 0x010(%%rdx)	,%%xmm11	\n\t"\
		"subpd		%%xmm6		,%%xmm5					\n\t		subpd		%%xmm14		,%%xmm13		\n\t"\
		"subpd		%%xmm2		,%%xmm1					\n\t		addpd		%%xmm10		,%%xmm9 		\n\t"\
		"addpd		%%xmm7		,%%xmm4					\n\t		addpd		%%xmm15		,%%xmm12		\n\t"\
		"addpd		%%xmm3		,%%xmm0					\n\t		subpd		%%xmm11		,%%xmm8 		\n\t"\
		"movaps		%%xmm5		,%%xmm7					\n\t		movaps		%%xmm13		,%%xmm15		\n\t"\
		"movaps		%%xmm4		,%%xmm6					\n\t		movaps		%%xmm12		,%%xmm14		\n\t"\
		"subpd		%%xmm0		,%%xmm4					\n\t		addpd		%%xmm8 		,%%xmm12		\n\t"\
		"subpd		%%xmm1		,%%xmm5					\n\t		addpd		%%xmm9 		,%%xmm13		\n\t"\
		"addpd		%%xmm0		,%%xmm6					\n\t		subpd		%%xmm8 		,%%xmm14		\n\t"\
		"addpd		%%xmm1		,%%xmm7					\n\t		subpd		%%xmm9 		,%%xmm15		\n\t"\
		"movaps		 0x040(%%rax)	,%%xmm2				\n\t		movaps		 0x0c0(%%rax)	,%%xmm10	\n\t"\
		"movaps		 0x050(%%rax)	,%%xmm3				\n\t		movaps		 0x0d0(%%rax)	,%%xmm11	\n\t"\
		"movaps		 0x040(%%rax)	,%%xmm0				\n\t		movaps		 0x0c0(%%rax)	,%%xmm8 	\n\t"\
		"movaps		 0x050(%%rax)	,%%xmm1				\n\t		movaps		 0x0d0(%%rax)	,%%xmm9 	\n\t"\
		"mulpd		      (%%rcx)	,%%xmm2				\n\t		mulpd		 0x010(%%rcx)	,%%xmm10	\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm1				\n\t		mulpd		      (%%rcx)	,%%xmm9 	\n\t"\
		"mulpd		      (%%rcx)	,%%xmm3				\n\t		mulpd		 0x010(%%rcx)	,%%xmm11	\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm0				\n\t		mulpd		      (%%rcx)	,%%xmm8 	\n\t"\
		"addpd		%%xmm1		,%%xmm2					\n\t		subpd		%%xmm9 		,%%xmm10		\n\t"\
		"subpd		%%xmm0		,%%xmm3					\n\t		addpd		%%xmm8 		,%%xmm11		\n\t"\
		"movaps		      (%%rax)	,%%xmm0				\n\t		movaps		 0x080(%%rax)	,%%xmm8 	\n\t"\
		"movaps		 0x010(%%rax)	,%%xmm1				\n\t		movaps		 0x090(%%rax)	,%%xmm9 	\n\t"\
		"subpd		%%xmm2		,%%xmm0					\n\t		subpd		%%xmm10		,%%xmm8 		\n\t"\
		"subpd		%%xmm3		,%%xmm1					\n\t		subpd		%%xmm11		,%%xmm9 		\n\t"\
		"addpd		%%xmm2		,%%xmm2					\n\t		addpd		%%xmm10		,%%xmm10		\n\t"\
		"addpd		%%xmm3		,%%xmm3					\n\t		addpd		%%xmm11		,%%xmm11		\n\t"\
		"addpd		%%xmm0		,%%xmm2					\n\t		addpd		%%xmm8 		,%%xmm10		\n\t"\
		"addpd		%%xmm1		,%%xmm3					\n\t		addpd		%%xmm9 		,%%xmm11		\n\t"\
		"subpd		%%xmm6		,%%xmm2					\n\t		subpd		%%xmm14		,%%xmm8 		\n\t"\
		"subpd		%%xmm7		,%%xmm3					\n\t		subpd		%%xmm15		,%%xmm9 		\n\t"\
		"addpd		%%xmm6		,%%xmm6					\n\t		addpd		%%xmm14		,%%xmm14		\n\t"\
		"addpd		%%xmm7		,%%xmm7					\n\t		addpd		%%xmm15		,%%xmm15		\n\t"\
		"addpd		%%xmm2		,%%xmm6					\n\t		addpd		%%xmm8 		,%%xmm14		\n\t"\
		"addpd		%%xmm3		,%%xmm7					\n\t		addpd		%%xmm9 		,%%xmm15		\n\t"\
		"movaps		%%xmm2		, 0x020(%%rax)			\n\t		movaps		%%xmm8 		, 0x0a0(%%rax)	\n\t"\
		"movaps		%%xmm3		, 0x030(%%rax)			\n\t		movaps		%%xmm9 		, 0x0b0(%%rax)	\n\t"\
		"movaps		%%xmm6		,%%xmm2					\n\t		movaps		%%xmm14		,%%xmm8 		\n\t"\
		"movaps		%%xmm7		,%%xmm3					\n\t		movaps		%%xmm15		,%%xmm9 		\n\t"\
		"mulpd		      (%%rdi)	,%%xmm6				\n\t		mulpd		 0x080(%%rdi)	,%%xmm14	\n\t"\
		"mulpd		      (%%rdi)	,%%xmm7				\n\t		mulpd		 0x080(%%rdi)	,%%xmm15	\n\t"\
		"mulpd		 0x010(%%rdi)	,%%xmm2				\n\t		mulpd		 0x090(%%rdi)	,%%xmm8 	\n\t"\
		"mulpd		 0x010(%%rdi)	,%%xmm3				\n\t		mulpd		 0x090(%%rdi)	,%%xmm9 	\n\t"\
		"subpd		%%xmm2		,%%xmm7					\n\t		subpd		%%xmm8 		,%%xmm15		\n\t"\
		"addpd		%%xmm3		,%%xmm6					\n\t		addpd		%%xmm9 		,%%xmm14		\n\t"\
		"movaps		%%xmm7		, 0x010(%%rbx)			\n\t		movaps		%%xmm15		, 0x050(%%rbx)	\n\t"\
		"movaps		%%xmm6		,      (%%rbx)			\n\t		movaps		%%xmm14		, 0x040(%%rbx)	\n\t"\
		"movaps		 0x020(%%rax)	,%%xmm6				\n\t		movaps		 0x0a0(%%rax)	,%%xmm14	\n\t"\
		"movaps		 0x030(%%rax)	,%%xmm7				\n\t		movaps		 0x0b0(%%rax)	,%%xmm15	\n\t"\
		"movaps		%%xmm6		,%%xmm2					\n\t		movaps		%%xmm14		,%%xmm8 		\n\t"\
		"movaps		%%xmm7		,%%xmm3					\n\t		movaps		%%xmm15		,%%xmm9 		\n\t"\
		"mulpd		 0x20 (%%rdi)	,%%xmm6				\n\t		mulpd		 0xa0 (%%rdi)	,%%xmm14	\n\t"\
		"mulpd		 0x20 (%%rdi)	,%%xmm7				\n\t		mulpd		 0xa0 (%%rdi)	,%%xmm15	\n\t"\
		"mulpd		 0x30 (%%rdi)	,%%xmm2				\n\t		mulpd		 0xb0 (%%rdi)	,%%xmm8 	\n\t"\
		"mulpd		 0x30 (%%rdi)	,%%xmm3				\n\t		mulpd		 0xb0 (%%rdi)	,%%xmm9 	\n\t"\
		"subpd		%%xmm2		,%%xmm7					\n\t		subpd		%%xmm8 		,%%xmm15		\n\t"\
		"addpd		%%xmm3		,%%xmm6					\n\t		addpd		%%xmm9 		,%%xmm14		\n\t"\
		"movaps		%%xmm7		, 0x110(%%rbx)			\n\t		movaps		%%xmm15		, 0x150(%%rbx)	\n\t"\
		"movaps		%%xmm6		, 0x100(%%rbx)			\n\t		movaps		%%xmm14		, 0x140(%%rbx)	\n\t"\
		"addq		$0x040		,%%rdi					\n\t"\
		"subpd		%%xmm5		,%%xmm0					\n\t		subpd		%%xmm13		,%%xmm10		\n\t"\
		"subpd		%%xmm4		,%%xmm1					\n\t		subpd		%%xmm12		,%%xmm11		\n\t"\
		"addpd		%%xmm5		,%%xmm5					\n\t		addpd		%%xmm13		,%%xmm13		\n\t"\
		"addpd		%%xmm4		,%%xmm4					\n\t		addpd		%%xmm12		,%%xmm12		\n\t"\
		"addpd		%%xmm0		,%%xmm5					\n\t		addpd		%%xmm10		,%%xmm13		\n\t"\
		"addpd		%%xmm1		,%%xmm4					\n\t		addpd		%%xmm11		,%%xmm12		\n\t"\
		"movaps		%%xmm5		,%%xmm2					\n\t		movaps		%%xmm13		,%%xmm8 		\n\t"\
		"movaps		%%xmm1		,%%xmm3					\n\t		movaps		%%xmm11		,%%xmm9 		\n\t"\
		"mulpd		      (%%rdi)	,%%xmm5				\n\t		mulpd		 0x080(%%rdi)	,%%xmm13	\n\t"\
		"mulpd		      (%%rdi)	,%%xmm1				\n\t		mulpd		 0x080(%%rdi)	,%%xmm11	\n\t"\
		"mulpd		 0x010(%%rdi)	,%%xmm2				\n\t		mulpd		 0x090(%%rdi)	,%%xmm8 	\n\t"\
		"mulpd		 0x010(%%rdi)	,%%xmm3				\n\t		mulpd		 0x090(%%rdi)	,%%xmm9 	\n\t"\
		"subpd		%%xmm2		,%%xmm1					\n\t		subpd		%%xmm8 		,%%xmm11		\n\t"\
		"addpd		%%xmm3		,%%xmm5					\n\t		addpd		%%xmm9 		,%%xmm13		\n\t"\
		"movaps		%%xmm1		, 0x090(%%rbx)			\n\t		movaps		%%xmm11		, 0x0d0(%%rbx)	\n\t"\
		"movaps		%%xmm5		, 0x080(%%rbx)			\n\t		movaps		%%xmm13		, 0x0c0(%%rbx)	\n\t"\
		"movaps		%%xmm0		,%%xmm2					\n\t		movaps		%%xmm10		,%%xmm8 		\n\t"\
		"movaps		%%xmm4		,%%xmm3					\n\t		movaps		%%xmm12		,%%xmm9 		\n\t"\
		"mulpd		 0x20 (%%rdi)	,%%xmm0				\n\t		mulpd		 0xa0 (%%rdi)	,%%xmm10	\n\t"\
		"mulpd		 0x20 (%%rdi)	,%%xmm4				\n\t		mulpd		 0xa0 (%%rdi)	,%%xmm12	\n\t"\
		"mulpd		 0x30 (%%rdi)	,%%xmm2				\n\t		mulpd		 0xb0 (%%rdi)	,%%xmm8 	\n\t"\
		"mulpd		 0x30 (%%rdi)	,%%xmm3				\n\t		mulpd		 0xb0 (%%rdi)	,%%xmm9 	\n\t"\
		"subpd		%%xmm2		,%%xmm4					\n\t		subpd		%%xmm8 		,%%xmm12		\n\t"\
		"addpd		%%xmm3		,%%xmm0					\n\t		addpd		%%xmm9 		,%%xmm10		\n\t"\
		"movaps		%%xmm4		, 0x190(%%rbx)			\n\t		movaps		%%xmm12		, 0x1d0(%%rbx)	\n\t"\
		"movaps		%%xmm0		, 0x180(%%rbx)			\n\t		movaps		%%xmm10		, 0x1c0(%%rbx)	\n\t"\
		"\n\t"\
		"/************************************************/	\n\t	/************************************************/	\n\t"\
		"/*...Block 4: t06,t16,t26,t36 -> r30,34,32,36:  */	\n\t	/*...Block 8: t0E,t1E,t2E,t3E -> r38,3C,3A,3E:  */	\n\t"\
		"/************************************************/	\n\t	/************************************************/	\n\t"\
		"addq		$0x200		,%%rax					\n\t		addq		$0x080		,%%rdi			\n\t"\
		"movaps		 0x020(%%rax)	,%%xmm4				\n\t		movaps		 0x0a0(%%rax)	,%%xmm12	\n\t"\
		"movaps		 0x060(%%rax)	,%%xmm0				\n\t		movaps		 0x0e0(%%rax)	,%%xmm8 	\n\t"\
		"movaps		 0x030(%%rax)	,%%xmm5				\n\t		movaps		 0x0b0(%%rax)	,%%xmm13	\n\t"\
		"movaps		 0x070(%%rax)	,%%xmm1				\n\t		movaps		 0x0f0(%%rax)	,%%xmm9 	\n\t"\
		"movaps		 0x020(%%rax)	,%%xmm6				\n\t		movaps		 0x0a0(%%rax)	,%%xmm14	\n\t"\
		"movaps		 0x060(%%rax)	,%%xmm2				\n\t		movaps		 0x0e0(%%rax)	,%%xmm10	\n\t"\
		"movaps		 0x030(%%rax)	,%%xmm7				\n\t		movaps		 0x0b0(%%rax)	,%%xmm15	\n\t"\
		"movaps		 0x070(%%rax)	,%%xmm3				\n\t		movaps		 0x0f0(%%rax)	,%%xmm11	\n\t"\
		"mulpd		 0x20 (%%rdx)	,%%xmm4				\n\t		mulpd		 0x010(%%rdx)	,%%xmm12	\n\t"\
		"mulpd		 0x010(%%rdx)	,%%xmm0				\n\t		mulpd		 0x30 (%%rdx)	,%%xmm8 	\n\t"\
		"mulpd		 0x20 (%%rdx)	,%%xmm5				\n\t		mulpd		 0x010(%%rdx)	,%%xmm13	\n\t"\
		"mulpd		 0x010(%%rdx)	,%%xmm1				\n\t		mulpd		 0x30 (%%rdx)	,%%xmm9 	\n\t"\
		"mulpd		 0x30 (%%rdx)	,%%xmm6				\n\t		mulpd		      (%%rdx)	,%%xmm14	\n\t"\
		"mulpd		      (%%rdx)	,%%xmm2				\n\t		mulpd		 0x20 (%%rdx)	,%%xmm10	\n\t"\
		"mulpd		 0x30 (%%rdx)	,%%xmm7				\n\t		mulpd		      (%%rdx)	,%%xmm15	\n\t"\
		"mulpd		      (%%rdx)	,%%xmm3				\n\t		mulpd		 0x20 (%%rdx)	,%%xmm11	\n\t"\
		"subpd		%%xmm6		,%%xmm5					\n\t		subpd		%%xmm14		,%%xmm13		\n\t"\
		"addpd		%%xmm2		,%%xmm1					\n\t		subpd		%%xmm10		,%%xmm9 		\n\t"\
		"addpd		%%xmm7		,%%xmm4					\n\t		addpd		%%xmm15		,%%xmm12		\n\t"\
		"subpd		%%xmm3		,%%xmm0					\n\t		addpd		%%xmm11		,%%xmm8 		\n\t"\
		"movaps		%%xmm5		,%%xmm7					\n\t		movaps		%%xmm13		,%%xmm15		\n\t"\
		"movaps		%%xmm4		,%%xmm6					\n\t		movaps		%%xmm12		,%%xmm14		\n\t"\
		"addpd		%%xmm0		,%%xmm4					\n\t		addpd		%%xmm8 		,%%xmm12		\n\t"\
		"addpd		%%xmm1		,%%xmm5					\n\t		addpd		%%xmm9 		,%%xmm13		\n\t"\
		"subpd		%%xmm0		,%%xmm6					\n\t		subpd		%%xmm8 		,%%xmm14		\n\t"\
		"subpd		%%xmm1		,%%xmm7					\n\t		subpd		%%xmm9 		,%%xmm15		\n\t"\
		"movaps		 0x040(%%rax)	,%%xmm2				\n\t		movaps		 0x0c0(%%rax)	,%%xmm10	\n\t"\
		"movaps		 0x050(%%rax)	,%%xmm3				\n\t		movaps		 0x0d0(%%rax)	,%%xmm11	\n\t"\
		"movaps		 0x040(%%rax)	,%%xmm0				\n\t		movaps		 0x0c0(%%rax)	,%%xmm8 	\n\t"\
		"movaps		 0x050(%%rax)	,%%xmm1				\n\t		movaps		 0x0d0(%%rax)	,%%xmm9 	\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm2				\n\t		mulpd		      (%%rcx)	,%%xmm10	\n\t"\
		"mulpd		      (%%rcx)	,%%xmm1				\n\t		mulpd		 0x010(%%rcx)	,%%xmm9 	\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm3				\n\t		mulpd		      (%%rcx)	,%%xmm11	\n\t"\
		"mulpd		      (%%rcx)	,%%xmm0				\n\t		mulpd		 0x010(%%rcx)	,%%xmm8 	\n\t"\
		"addpd		%%xmm1		,%%xmm2					\n\t		subpd		%%xmm9 		,%%xmm10		\n\t"\
		"subpd		%%xmm0		,%%xmm3					\n\t		addpd		%%xmm8 		,%%xmm11		\n\t"\
		"movaps		      (%%rax)	,%%xmm0				\n\t		movaps		 0x080(%%rax)	,%%xmm8 	\n\t"\
		"movaps		 0x010(%%rax)	,%%xmm1				\n\t		movaps		 0x090(%%rax)	,%%xmm9 	\n\t"\
		"subpd		%%xmm2		,%%xmm0					\n\t		subpd		%%xmm10		,%%xmm8 		\n\t"\
		"subpd		%%xmm3		,%%xmm1					\n\t		subpd		%%xmm11		,%%xmm9 		\n\t"\
		"addpd		%%xmm2		,%%xmm2					\n\t		addpd		%%xmm10		,%%xmm10		\n\t"\
		"addpd		%%xmm3		,%%xmm3					\n\t		addpd		%%xmm11		,%%xmm11		\n\t"\
		"addpd		%%xmm0		,%%xmm2					\n\t		addpd		%%xmm8 		,%%xmm10		\n\t"\
		"addpd		%%xmm1		,%%xmm3					\n\t		addpd		%%xmm9 		,%%xmm11		\n\t"\
		"addq		$0x040		,%%rdi					\n\t"\
		"subpd		%%xmm6		,%%xmm2					\n\t		subpd		%%xmm14		,%%xmm8 		\n\t"\
		"subpd		%%xmm7		,%%xmm3					\n\t		subpd		%%xmm15		,%%xmm9 		\n\t"\
		"addpd		%%xmm6		,%%xmm6					\n\t		addpd		%%xmm14		,%%xmm14		\n\t"\
		"addpd		%%xmm7		,%%xmm7					\n\t		addpd		%%xmm15		,%%xmm15		\n\t"\
		"addpd		%%xmm2		,%%xmm6					\n\t		addpd		%%xmm8 		,%%xmm14		\n\t"\
		"addpd		%%xmm3		,%%xmm7					\n\t		addpd		%%xmm9 		,%%xmm15		\n\t"\
		"movaps		%%xmm2		, 0x020(%%rax)			\n\t		movaps		%%xmm8 		, 0x0a0(%%rax)	\n\t"\
		"movaps		%%xmm3		, 0x030(%%rax)			\n\t		movaps		%%xmm9 		, 0x0b0(%%rax)	\n\t"\
		"movaps		%%xmm6		,%%xmm2					\n\t		movaps		%%xmm14		,%%xmm8 		\n\t"\
		"movaps		%%xmm7		,%%xmm3					\n\t		movaps		%%xmm15		,%%xmm9 		\n\t"\
		"mulpd		      (%%rdi)	,%%xmm6				\n\t		mulpd		 0x080(%%rdi)	,%%xmm14	\n\t"\
		"mulpd		      (%%rdi)	,%%xmm7				\n\t		mulpd		 0x080(%%rdi)	,%%xmm15	\n\t"\
		"mulpd		 0x010(%%rdi)	,%%xmm2				\n\t		mulpd		 0x090(%%rdi)	,%%xmm8 	\n\t"\
		"mulpd		 0x010(%%rdi)	,%%xmm3				\n\t		mulpd		 0x090(%%rdi)	,%%xmm9 	\n\t"\
		"subpd		%%xmm2		,%%xmm7					\n\t		subpd		%%xmm8 		,%%xmm15		\n\t"\
		"addpd		%%xmm3		,%%xmm6					\n\t		addpd		%%xmm9 		,%%xmm14		\n\t"\
		"movaps		%%xmm7		, 0x030(%%rbx)			\n\t		movaps		%%xmm15		, 0x070(%%rbx)	\n\t"\
		"movaps		%%xmm6		, 0x020(%%rbx)			\n\t		movaps		%%xmm14		, 0x060(%%rbx)	\n\t"\
		"movaps		 0x020(%%rax)	,%%xmm6				\n\t		movaps		 0x0a0(%%rax)	,%%xmm14	\n\t"\
		"movaps		 0x030(%%rax)	,%%xmm7				\n\t		movaps		 0x0b0(%%rax)	,%%xmm15	\n\t"\
		"movaps		%%xmm6		,%%xmm2					\n\t		movaps		%%xmm14		,%%xmm8 		\n\t"\
		"movaps		%%xmm7		,%%xmm3					\n\t		movaps		%%xmm15		,%%xmm9 		\n\t"\
		"mulpd		 0x20 (%%rdi)	,%%xmm6				\n\t		mulpd		 0xa0 (%%rdi)	,%%xmm14	\n\t"\
		"mulpd		 0x20 (%%rdi)	,%%xmm7				\n\t		mulpd		 0xa0 (%%rdi)	,%%xmm15	\n\t"\
		"mulpd		 0x30 (%%rdi)	,%%xmm2				\n\t		mulpd		 0xb0 (%%rdi)	,%%xmm8 	\n\t"\
		"mulpd		 0x30 (%%rdi)	,%%xmm3				\n\t		mulpd		 0xb0 (%%rdi)	,%%xmm9 	\n\t"\
		"subpd		%%xmm2		,%%xmm7					\n\t		subpd		%%xmm8 		,%%xmm15		\n\t"\
		"addpd		%%xmm3		,%%xmm6					\n\t		addpd		%%xmm9 		,%%xmm14		\n\t"\
		"movaps		%%xmm7		, 0x130(%%rbx)			\n\t		movaps		%%xmm15		, 0x170(%%rbx)	\n\t"\
		"movaps		%%xmm6		, 0x120(%%rbx)			\n\t		movaps		%%xmm14		, 0x160(%%rbx)	\n\t"\
		"addq		$0x040		,%%rdi					\n\t"\
		"subpd		%%xmm5		,%%xmm0					\n\t		subpd		%%xmm13		,%%xmm10		\n\t"\
		"subpd		%%xmm4		,%%xmm1					\n\t		subpd		%%xmm12		,%%xmm11		\n\t"\
		"addpd		%%xmm5		,%%xmm5					\n\t		addpd		%%xmm13		,%%xmm13		\n\t"\
		"addpd		%%xmm4		,%%xmm4					\n\t		addpd		%%xmm12		,%%xmm12		\n\t"\
		"addpd		%%xmm0		,%%xmm5					\n\t		addpd		%%xmm10		,%%xmm13		\n\t"\
		"addpd		%%xmm1		,%%xmm4					\n\t		addpd		%%xmm11		,%%xmm12		\n\t"\
		"movaps		%%xmm5		,%%xmm2					\n\t		movaps		%%xmm13		,%%xmm8 		\n\t"\
		"movaps		%%xmm1		,%%xmm3					\n\t		movaps		%%xmm11		,%%xmm9 		\n\t"\
		"mulpd		      (%%rdi)	,%%xmm5				\n\t		mulpd		 0x080(%%rdi)	,%%xmm13	\n\t"\
		"mulpd		      (%%rdi)	,%%xmm1				\n\t		mulpd		 0x080(%%rdi)	,%%xmm11	\n\t"\
		"mulpd		 0x010(%%rdi)	,%%xmm2				\n\t		mulpd		 0x090(%%rdi)	,%%xmm8 	\n\t"\
		"mulpd		 0x010(%%rdi)	,%%xmm3				\n\t		mulpd		 0x090(%%rdi)	,%%xmm9 	\n\t"\
		"subpd		%%xmm2		,%%xmm1					\n\t		subpd		%%xmm8 		,%%xmm11		\n\t"\
		"addpd		%%xmm3		,%%xmm5					\n\t		addpd		%%xmm9 		,%%xmm13		\n\t"\
		"movaps		%%xmm1		, 0x0b0(%%rbx)			\n\t		movaps		%%xmm11		, 0x0f0(%%rbx)	\n\t"\
		"movaps		%%xmm5		, 0x0a0(%%rbx)			\n\t		movaps		%%xmm13		, 0x0e0(%%rbx)	\n\t"\
		"movaps		%%xmm0		,%%xmm2					\n\t		movaps		%%xmm10		,%%xmm8 		\n\t"\
		"movaps		%%xmm4		,%%xmm3					\n\t		movaps		%%xmm12		,%%xmm9 		\n\t"\
		"mulpd		 0x20 (%%rdi)	,%%xmm0				\n\t		mulpd		 0xa0 (%%rdi)	,%%xmm10	\n\t"\
		"mulpd		 0x20 (%%rdi)	,%%xmm4				\n\t		mulpd		 0xa0 (%%rdi)	,%%xmm12	\n\t"\
		"mulpd		 0x30 (%%rdi)	,%%xmm2				\n\t		mulpd		 0xb0 (%%rdi)	,%%xmm8 	\n\t"\
		"mulpd		 0x30 (%%rdi)	,%%xmm3				\n\t		mulpd		 0xb0 (%%rdi)	,%%xmm9 	\n\t"\
		"subpd		%%xmm2		,%%xmm4					\n\t		subpd		%%xmm8 		,%%xmm12		\n\t"\
		"addpd		%%xmm3		,%%xmm0					\n\t		addpd		%%xmm9 		,%%xmm10		\n\t"\
		"movaps		%%xmm4		, 0x1b0(%%rbx)			\n\t		movaps		%%xmm12		, 0x1f0(%%rbx)	\n\t"\
		"movaps		%%xmm0		, 0x1a0(%%rbx)			\n\t		movaps		%%xmm10		, 0x1e0(%%rbx)	\n\t"\
		"/************************************************/	\n\t		/************************************************/	\n\t"\
		"/*...Block 1: t00,t10,t20,t30 -> r00,04,02,06:  */	\n\t		/*...Block 5: t08,t18,t28,t38 -> r08,0C,0A,0E:  */	\n\t"\
		"/************************************************/	\n\t		/************************************************/	\n\t"\
		"movq	%[__r00],%%rdx	/* base-addr in rcol = c04/r08, so rcx,rdx offset +0x80 in rcol */	\n\t	movaps	(%%rsi),%%xmm10	/* isrt2 */	\n\t"\
		"movaps		      (%%rdx)	,%%xmm0					\n\t		movaps		 0x0a0(%%rdx)	,%%xmm12				\n\t"\
		"movaps		 0x010(%%rdx)	,%%xmm1					\n\t		movaps		 0x0b0(%%rdx)	,%%xmm13				\n\t"\
		"movaps		 0x040(%%rdx)	,%%xmm2					\n\t		movaps		 0x0e0(%%rdx)	,%%xmm8 				\n\t"\
		"movaps		 0x050(%%rdx)	,%%xmm3					\n\t		movaps		 0x0f0(%%rdx)	,%%xmm9 				\n\t"\
		"subpd		 0x040(%%rdx)	,%%xmm0					\n\t		addpd		 0x0b0(%%rdx)	,%%xmm12				\n\t"\
		"subpd		 0x050(%%rdx)	,%%xmm1					\n\t		subpd		 0x0a0(%%rdx)	,%%xmm13				\n\t"\
		"addpd		      (%%rdx)	,%%xmm2					\n\t		subpd		 0x0f0(%%rdx)	,%%xmm8 				\n\t"\
		"addpd		 0x010(%%rdx)	,%%xmm3					\n\t		addpd		 0x0e0(%%rdx)	,%%xmm9 				\n\t"\
		"movaps		 0x020(%%rdx)	,%%xmm4					\n\t		mulpd		%%xmm10		,%%xmm12					\n\t"\
		"movaps		 0x030(%%rdx)	,%%xmm5					\n\t		mulpd		%%xmm10		,%%xmm13					\n\t"\
		"movaps		 0x060(%%rdx)	,%%xmm6					\n\t		mulpd		%%xmm10		,%%xmm8 					\n\t"\
		"movaps		 0x070(%%rdx)	,%%xmm7					\n\t		mulpd		%%xmm10		,%%xmm9 					\n\t"\
		"subpd		 0x060(%%rdx)	,%%xmm4					\n\t		movaps		%%xmm12		,%%xmm14					\n\t"\
		"subpd		 0x070(%%rdx)	,%%xmm5					\n\t		movaps		%%xmm13		,%%xmm15					\n\t"\
		"addpd		 0x020(%%rdx)	,%%xmm6					\n\t		subpd		%%xmm8 		,%%xmm12					\n\t"\
		"addpd		 0x030(%%rdx)	,%%xmm7					\n\t		subpd		%%xmm9 		,%%xmm13					\n\t"\
		"movq		%[__add0]		,%%rax					\n\t		addpd		%%xmm8 		,%%xmm14					\n\t"\
		"movq		%[__c10]		,%%rcx					\n\t		addpd		%%xmm9 		,%%xmm15					\n\t"\
		"addpd		%%xmm6		,%%xmm2						\n\t		movaps		 0x080(%%rdx)	,%%xmm8 				\n\t"\
		"addpd		%%xmm7		,%%xmm3						\n\t		movaps		 0x090(%%rdx)	,%%xmm9 				\n\t"\
		"movaps		%%xmm2		,      (%%rdx)				\n\t		movaps		 0x0c0(%%rdx)	,%%xmm10				\n\t"\
		"movaps		%%xmm3		, 0x010(%%rdx)				\n\t		movaps		 0x0d0(%%rdx)	,%%xmm11				\n\t"\
		"addpd		%%xmm6		,%%xmm6						\n\t		subpd		 0x0d0(%%rdx)	,%%xmm8 				\n\t"\
		"addpd		%%xmm7		,%%xmm7						\n\t		subpd		 0x0c0(%%rdx)	,%%xmm9 				\n\t"\
		"subpd		%%xmm6		,%%xmm2						\n\t		addpd		 0x080(%%rdx)	,%%xmm11				\n\t"\
		"subpd		%%xmm7		,%%xmm3						\n\t		addpd		 0x090(%%rdx)	,%%xmm10				\n\t"\
		"movaps		%%xmm2		,%%xmm6						\n\t		subpd		%%xmm12		,%%xmm11					\n\t"\
		"movaps		%%xmm3		,%%xmm7						\n\t		subpd		%%xmm13		,%%xmm9 					\n\t"\
		"mulpd		      (%%rcx)	,%%xmm2					\n\t		addpd		%%xmm12		,%%xmm12					\n\t"\
		"mulpd		      (%%rcx)	,%%xmm3					\n\t		addpd		%%xmm13		,%%xmm13					\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm6					\n\t		addpd		%%xmm11		,%%xmm12					\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm7					\n\t		addpd		%%xmm9 		,%%xmm13					\n\t"\
		"subpd		%%xmm6		,%%xmm3						\n\t		movaps		%%xmm11		, 0x080(%%rdx)				\n\t"\
		"addpd		%%xmm7		,%%xmm2						\n\t		movaps		%%xmm9 		, 0x090(%%rdx)				\n\t"\
		"subq	$0x20,%%rcx	/* put c00 in rcx to ease bookkeeping*/\n\t	movaps		%%xmm12		,%%xmm11					\n\t"\
		"movq		%[__add1]		,%%rbx					\n\t		movaps		%%xmm13		,%%xmm9 					\n\t"\
		"movaps		%%xmm3		,%%xmm7						\n\t		mulpd		 0x080(%%rcx)	,%%xmm12	/* c04 */	\n\t"\
		"movaps		%%xmm2		,%%xmm6						\n\t		mulpd		 0x080(%%rcx)	,%%xmm13				\n\t"\
		"unpckhpd	 0x110(%%rbx)	,%%xmm7					\n\t		mulpd		 0x090(%%rcx)	,%%xmm11				\n\t"\
		"unpcklpd	 0x110(%%rbx)	,%%xmm3					\n\t		mulpd		 0x090(%%rcx)	,%%xmm9 				\n\t"\
		"movaps		%%xmm7		, 0x110(%%rbx)				\n\t		subpd		%%xmm11		,%%xmm13					\n\t"\
		"unpckhpd	 0x100(%%rbx)	,%%xmm6					\n\t		addpd		%%xmm9 		,%%xmm12					\n\t"\
		"unpcklpd	 0x100(%%rbx)	,%%xmm2					\n\t		movaps		%%xmm13		,%%xmm11					\n\t"\
		"movaps		%%xmm6		, 0x100(%%rbx)				\n\t		movaps		%%xmm12		,%%xmm9 					\n\t"\
		"movaps		%%xmm3		, 0x110(%%rax)				\n\t		unpckhpd	 0x050(%%rbx)	,%%xmm11				\n\t"\
		"movaps		%%xmm2		, 0x100(%%rax)				\n\t		unpcklpd	 0x050(%%rbx)	,%%xmm13				\n\t"\
		"movaps		 0x010(%%rdx)	,%%xmm3					\n\t		movaps		%%xmm11		, 0x050(%%rbx)				\n\t"\
		"movaps		      (%%rdx)	,%%xmm2					\n\t		unpckhpd	 0x040(%%rbx)	,%%xmm9 				\n\t"\
		"movaps		%%xmm3		,%%xmm7						\n\t		unpcklpd	 0x040(%%rbx)	,%%xmm12				\n\t"\
		"movaps		%%xmm2		,%%xmm6						\n\t		movaps		%%xmm9 		, 0x040(%%rbx)				\n\t"\
		"unpckhpd	 0x010(%%rbx)	,%%xmm7					\n\t		movaps		%%xmm13		, 0x050(%%rax)				\n\t"\
		"unpcklpd	 0x010(%%rbx)	,%%xmm3					\n\t		movaps		%%xmm12		, 0x040(%%rax)				\n\t"\
		"movaps		%%xmm7		, 0x010(%%rbx)				\n\t		movaps		 0x080(%%rdx)	,%%xmm12				\n\t"\
		"unpckhpd	      (%%rbx)	,%%xmm6					\n\t		movaps		 0x090(%%rdx)	,%%xmm13				\n\t"\
		"unpcklpd	      (%%rbx)	,%%xmm2					\n\t		movaps		%%xmm12		,%%xmm11					\n\t"\
		"movaps		%%xmm6		,      (%%rbx)				\n\t		movaps		%%xmm13		,%%xmm9 					\n\t"\
		"movaps		%%xmm3		, 0x010(%%rax)				\n\t		mulpd		 0x0a0(%%rcx)	,%%xmm12	/* c14 */	\n\t"\
		"movaps		%%xmm2		,      (%%rax)				\n\t		mulpd		 0x0a0(%%rcx)	,%%xmm13				\n\t"\
		"addpd		%%xmm5		,%%xmm0						\n\t		mulpd		 0x0b0(%%rcx)	,%%xmm11				\n\t"\
		"subpd		%%xmm4		,%%xmm1						\n\t		mulpd		 0x0b0(%%rcx)	,%%xmm9 				\n\t"\
		"movaps		%%xmm0		,%%xmm2						\n\t		subpd		%%xmm11		,%%xmm13					\n\t"\
		"movaps		%%xmm1		,%%xmm3						\n\t		addpd		%%xmm9 		,%%xmm12					\n\t"\
		"addpd		%%xmm5		,%%xmm5						\n\t		movaps		%%xmm13		,%%xmm11					\n\t"\
		"addpd		%%xmm4		,%%xmm4						\n\t		movaps		%%xmm12		,%%xmm9 					\n\t"\
		"movaps		%%xmm0		,%%xmm6						\n\t		unpckhpd	 0x150(%%rbx)	,%%xmm11				\n\t"\
		"movaps		%%xmm1		,%%xmm7						\n\t		unpcklpd	 0x150(%%rbx)	,%%xmm13				\n\t"\
		"mulpd		 0x040(%%rcx)	,%%xmm2		/* c08 */	\n\t		movaps		%%xmm11		, 0x150(%%rbx)				\n\t"\
		"mulpd		 0x040(%%rcx)	,%%xmm3					\n\t		unpckhpd	 0x140(%%rbx)	,%%xmm9 				\n\t"\
		"mulpd		 0x050(%%rcx)	,%%xmm6					\n\t		unpcklpd	 0x140(%%rbx)	,%%xmm12				\n\t"\
		"mulpd		 0x050(%%rcx)	,%%xmm7					\n\t		movaps		%%xmm9 		, 0x140(%%rbx)				\n\t"\
		"subpd		%%xmm6		,%%xmm3						\n\t		movaps		%%xmm13		, 0x150(%%rax)				\n\t"\
		"addpd		%%xmm7		,%%xmm2						\n\t		movaps		%%xmm12		, 0x140(%%rax)				\n\t"\
		"movaps		%%xmm3		,%%xmm7						\n\t		subpd		%%xmm15		,%%xmm8 					\n\t"\
		"movaps		%%xmm2		,%%xmm6						\n\t		subpd		%%xmm14		,%%xmm10					\n\t"\
		"unpckhpd	 0x090(%%rbx)	,%%xmm7					\n\t		addpd		%%xmm15		,%%xmm15					\n\t"\
		"unpcklpd	 0x090(%%rbx)	,%%xmm3					\n\t		addpd		%%xmm14		,%%xmm14					\n\t"\
		"movaps		%%xmm7		, 0x090(%%rbx)				\n\t		addpd		%%xmm8 		,%%xmm15					\n\t"\
		"unpckhpd	 0x080(%%rbx)	,%%xmm6					\n\t		addpd		%%xmm10		,%%xmm14					\n\t"\
		"unpcklpd	 0x080(%%rbx)	,%%xmm2					\n\t		movaps		%%xmm15		,%%xmm12					\n\t"\
		"movaps		%%xmm6		, 0x080(%%rbx)				\n\t		movaps		%%xmm10		,%%xmm13					\n\t"\
		"movaps		%%xmm3		, 0x090(%%rax)				\n\t		mulpd		 0x0c0(%%rcx)	,%%xmm15	/* c0C */	\n\t"\
		"movaps		%%xmm2		, 0x080(%%rax)				\n\t		mulpd		 0x0c0(%%rcx)	,%%xmm10				\n\t"\
		"subpd		%%xmm5		,%%xmm0						\n\t		mulpd		 0x0d0(%%rcx)	,%%xmm12				\n\t"\
		"addpd		%%xmm4		,%%xmm1						\n\t		mulpd		 0x0d0(%%rcx)	,%%xmm13				\n\t"\
		"movaps		%%xmm0		,%%xmm6						\n\t		subpd		%%xmm12		,%%xmm10					\n\t"\
		"movaps		%%xmm1		,%%xmm7						\n\t		addpd		%%xmm13		,%%xmm15					\n\t"\
		"mulpd		 0x060(%%rcx)	,%%xmm0		/* c18 */	\n\t		movaps		%%xmm10		,%%xmm13					\n\t"\
		"mulpd		 0x060(%%rcx)	,%%xmm1					\n\t		movaps		%%xmm15		,%%xmm12					\n\t"\
		"mulpd		 0x070(%%rcx)	,%%xmm6					\n\t		unpckhpd	 0x0d0(%%rbx)	,%%xmm13				\n\t"\
		"mulpd		 0x070(%%rcx)	,%%xmm7					\n\t		unpcklpd	 0x0d0(%%rbx)	,%%xmm10				\n\t"\
		"subpd		%%xmm6		,%%xmm1						\n\t		movaps		%%xmm13		, 0x0d0(%%rbx)				\n\t"\
		"addpd		%%xmm7		,%%xmm0						\n\t		unpckhpd	 0x0c0(%%rbx)	,%%xmm12				\n\t"\
		"movaps		%%xmm1		,%%xmm7						\n\t		unpcklpd	 0x0c0(%%rbx)	,%%xmm15				\n\t"\
		"movaps		%%xmm0		,%%xmm6						\n\t		movaps		%%xmm12		, 0x0c0(%%rbx)				\n\t"\
		"unpckhpd	 0x190(%%rbx)	,%%xmm7					\n\t		movaps		%%xmm10		, 0x0d0(%%rax)				\n\t"\
		"unpcklpd	 0x190(%%rbx)	,%%xmm1					\n\t		movaps		%%xmm15		, 0x0c0(%%rax)				\n\t"\
		"movaps		%%xmm7		, 0x190(%%rbx)				\n\t		movaps		%%xmm8 		,%%xmm12					\n\t"\
		"unpckhpd	 0x180(%%rbx)	,%%xmm6					\n\t		movaps		%%xmm14		,%%xmm13					\n\t"\
		"unpcklpd	 0x180(%%rbx)	,%%xmm0					\n\t		mulpd		 0x0e0(%%rcx)	,%%xmm8 	/* c1C */	\n\t"\
		"movaps		%%xmm6		, 0x180(%%rbx)				\n\t		mulpd		 0x0e0(%%rcx)	,%%xmm14				\n\t"\
		"movaps		%%xmm1		, 0x190(%%rax)				\n\t		mulpd		 0x0f0(%%rcx)	,%%xmm12				\n\t"\
		"movaps		%%xmm0		, 0x180(%%rax)				\n\t		mulpd		 0x0f0(%%rcx)	,%%xmm13				\n\t"\
		"																subpd		%%xmm12		,%%xmm14					\n\t"\
		"																addpd		%%xmm13		,%%xmm8 					\n\t"\
		"																movaps		%%xmm14		,%%xmm13					\n\t"\
		"																movaps		%%xmm8 		,%%xmm12					\n\t"\
		"																unpckhpd	 0x1d0(%%rbx)	,%%xmm13				\n\t"\
		"																unpcklpd	 0x1d0(%%rbx)	,%%xmm14				\n\t"\
		"																movaps		%%xmm13		, 0x1d0(%%rbx)				\n\t"\
		"																unpckhpd	 0x1c0(%%rbx)	,%%xmm12				\n\t"\
		"																unpcklpd	 0x1c0(%%rbx)	,%%xmm8 				\n\t"\
		"																movaps		%%xmm12		, 0x1c0(%%rbx)				\n\t"\
		"																movaps		%%xmm14		, 0x1d0(%%rax)				\n\t"\
		"																movaps		%%xmm8 		, 0x1c0(%%rax)				\n\t"\
		"/************************************************/	\n\t		/************************************************/	\n\t"\
		"/*...Block 3: t04,t14,t24,t34 -> r20,24,22,26:  */	\n\t		/*...Block 7: t0C,t1C,t2C,t3C -> r28,2C,2A,2E:  */	\n\t"\
		"/************************************************/	\n\t		/************************************************/	\n\t"\
		"movq		%[__r20],%%rdx								/* base-addr in rcol = r28, so rdx offset +0x80 vs lcol */	\n\t"\
		"leaq	0x010(%%rsi),%%rcx	/* cc0; Note cc0/ss0 are shared between lcol/rcol, so no rcx-offset until get to twiddles*/\n\t"\
		"movaps		 0x020(%%rdx)	,%%xmm4					\n\t		movaps		 0x0a0(%%rdx)	,%%xmm12				\n\t"\
		"movaps		 0x060(%%rdx)	,%%xmm0					\n\t		movaps		 0x0e0(%%rdx)	,%%xmm8 				\n\t"\
		"movaps		 0x030(%%rdx)	,%%xmm5					\n\t		movaps		 0x0b0(%%rdx)	,%%xmm13				\n\t"\
		"movaps		 0x070(%%rdx)	,%%xmm1					\n\t		movaps		 0x0f0(%%rdx)	,%%xmm9 				\n\t"\
		"movaps		 0x020(%%rdx)	,%%xmm6					\n\t		movaps		 0x0a0(%%rdx)	,%%xmm14				\n\t"\
		"movaps		 0x060(%%rdx)	,%%xmm2					\n\t		movaps		 0x0e0(%%rdx)	,%%xmm10				\n\t"\
		"movaps		 0x030(%%rdx)	,%%xmm7					\n\t		movaps		 0x0b0(%%rdx)	,%%xmm15				\n\t"\
		"movaps		 0x070(%%rdx)	,%%xmm3					\n\t		movaps		 0x0f0(%%rdx)	,%%xmm11				\n\t"\
		"mulpd		      (%%rcx)	,%%xmm4					\n\t		mulpd		 0x010(%%rcx)	,%%xmm12				\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm0					\n\t		mulpd		      (%%rcx)	,%%xmm8 				\n\t"\
		"mulpd		      (%%rcx)	,%%xmm5					\n\t		mulpd		 0x010(%%rcx)	,%%xmm13				\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm1					\n\t		mulpd		      (%%rcx)	,%%xmm9 				\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm6					\n\t		mulpd		      (%%rcx)	,%%xmm14				\n\t"\
		"mulpd		      (%%rcx)	,%%xmm2					\n\t		mulpd		 0x010(%%rcx)	,%%xmm10				\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm7					\n\t		mulpd		      (%%rcx)	,%%xmm15				\n\t"\
		"mulpd		      (%%rcx)	,%%xmm3					\n\t		mulpd		 0x010(%%rcx)	,%%xmm11				\n\t"\
		"subpd		%%xmm6		,%%xmm5						\n\t		subpd		%%xmm14		,%%xmm13					\n\t"\
		"subpd		%%xmm2		,%%xmm1						\n\t		subpd		%%xmm10		,%%xmm9 					\n\t"\
		"addpd		%%xmm7		,%%xmm4						\n\t		addpd		%%xmm15		,%%xmm12					\n\t"\
		"addpd		%%xmm3		,%%xmm0						\n\t		addpd		%%xmm11		,%%xmm8 					\n\t"\
		"movaps		%%xmm5		,%%xmm7						\n\t		movaps		%%xmm13		,%%xmm15					\n\t"\
		"movaps		%%xmm4		,%%xmm6						\n\t		movaps		%%xmm12		,%%xmm14					\n\t"\
		"addpd		%%xmm0		,%%xmm4						\n\t		addpd		%%xmm8 		,%%xmm12					\n\t"\
		"addpd		%%xmm1		,%%xmm5						\n\t		addpd		%%xmm9 		,%%xmm13					\n\t"\
		"subpd		%%xmm0		,%%xmm6						\n\t		subpd		%%xmm8 		,%%xmm14					\n\t"\
		"subpd		%%xmm1		,%%xmm7						\n\t		subpd		%%xmm9 		,%%xmm15					\n\t"\
		"movaps		 0x040(%%rdx)	,%%xmm2					\n\t		movaps		 0x0c0(%%rdx)	,%%xmm10				\n\t"\
		"movaps		 0x050(%%rdx)	,%%xmm3					\n\t		movaps		 0x0d0(%%rdx)	,%%xmm11				\n\t"\
		"movaps		      (%%rdx)	,%%xmm0					\n\t		movaps		 0x080(%%rdx)	,%%xmm8 				\n\t"\
		"movaps		 0x010(%%rdx)	,%%xmm1					\n\t		movaps		 0x090(%%rdx)	,%%xmm9 				\n\t"\
		"addpd		 0x050(%%rdx)	,%%xmm2					\n\t		subpd		 0x0d0(%%rdx)	,%%xmm10				\n\t"\
		"subpd		 0x040(%%rdx)	,%%xmm3					\n\t		addpd		 0x0c0(%%rdx)	,%%xmm11				\n\t"\
		"mulpd		      (%%rsi)	,%%xmm2					\n\t		mulpd		      (%%rsi)	,%%xmm10				\n\t"\
		"mulpd		      (%%rsi)	,%%xmm3					\n\t		mulpd		      (%%rsi)	,%%xmm11				\n\t"\
		"subpd		%%xmm2		,%%xmm0						\n\t		subpd		%%xmm10		,%%xmm8 					\n\t"\
		"subpd		%%xmm3		,%%xmm1						\n\t		subpd		%%xmm11		,%%xmm9 					\n\t"\
		"addpd		%%xmm2		,%%xmm2						\n\t		addpd		%%xmm10		,%%xmm10					\n\t"\
		"addpd		%%xmm3		,%%xmm3						\n\t		addpd		%%xmm11		,%%xmm11					\n\t"\
		"addpd		%%xmm0		,%%xmm2						\n\t		addpd		%%xmm8 		,%%xmm10					\n\t"\
		"addpd		%%xmm1		,%%xmm3						\n\t		addpd		%%xmm9 		,%%xmm11					\n\t"\
		"movq		%[__add0]		,%%rax					\n\t"\
		"movq		%[__c02]		,%%rcx					\n\t"\
		"subpd		%%xmm4		,%%xmm2						\n\t		subpd		%%xmm14		,%%xmm8 					\n\t"\
		"subpd		%%xmm5		,%%xmm3						\n\t		subpd		%%xmm15		,%%xmm9 					\n\t"\
		"addpd		%%xmm4		,%%xmm4						\n\t		addpd		%%xmm14		,%%xmm14					\n\t"\
		"addpd		%%xmm5		,%%xmm5						\n\t		addpd		%%xmm15		,%%xmm15					\n\t"\
		"addpd		%%xmm2		,%%xmm4						\n\t		addpd		%%xmm8 		,%%xmm14					\n\t"\
		"addpd		%%xmm3		,%%xmm5						\n\t		addpd		%%xmm9 		,%%xmm15					\n\t"\
		"movaps		%%xmm2		,      (%%rdx)				\n\t		movaps		%%xmm8 		, 0x080(%%rdx)				\n\t"\
		"movaps		%%xmm3		, 0x010(%%rdx)				\n\t		movaps		%%xmm9 		, 0x090(%%rdx)				\n\t"\
		"movaps		%%xmm4		,%%xmm2						\n\t		movaps		%%xmm14		,%%xmm8 					\n\t"\
		"movaps		%%xmm5		,%%xmm3						\n\t		movaps		%%xmm15		,%%xmm9 					\n\t"\
		"mulpd		      (%%rcx)	,%%xmm4					\n\t		mulpd		 0x080(%%rcx)	,%%xmm14	/* c06 */	\n\t"\
		"mulpd		      (%%rcx)	,%%xmm5					\n\t		mulpd		 0x080(%%rcx)	,%%xmm15				\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm2					\n\t		mulpd		 0x090(%%rcx)	,%%xmm8 				\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm3					\n\t		mulpd		 0x090(%%rcx)	,%%xmm9 				\n\t"\
		"subpd		%%xmm2		,%%xmm5						\n\t		subpd		%%xmm8 		,%%xmm15					\n\t"\
		"addpd		%%xmm3		,%%xmm4						\n\t		addpd		%%xmm9 		,%%xmm14					\n\t"\
		"movq		%[__add1]		,%%rbx					\n\t"\
		"movaps		%%xmm5		,%%xmm3						\n\t		movaps		%%xmm15		,%%xmm9 					\n\t"\
		"movaps		%%xmm4		,%%xmm2						\n\t		movaps		%%xmm14		,%%xmm8 					\n\t"\
		"unpckhpd	 0x030(%%rbx)	,%%xmm3					\n\t		unpckhpd	 0x070(%%rbx)	,%%xmm9 				\n\t"\
		"unpcklpd	 0x030(%%rbx)	,%%xmm5					\n\t		unpcklpd	 0x070(%%rbx)	,%%xmm15				\n\t"\
		"movaps		%%xmm3		, 0x030(%%rbx)				\n\t		movaps		%%xmm9 		, 0x070(%%rbx)				\n\t"\
		"unpckhpd	 0x020(%%rbx)	,%%xmm2					\n\t		unpckhpd	 0x060(%%rbx)	,%%xmm8 				\n\t"\
		"unpcklpd	 0x020(%%rbx)	,%%xmm4					\n\t		unpcklpd	 0x060(%%rbx)	,%%xmm14				\n\t"\
		"movaps		%%xmm2		, 0x020(%%rbx)				\n\t		movaps		%%xmm8 		, 0x060(%%rbx)				\n\t"\
		"movaps		%%xmm5		, 0x030(%%rax)				\n\t		movaps		%%xmm15		, 0x070(%%rax)				\n\t"\
		"movaps		%%xmm4		, 0x020(%%rax)				\n\t		movaps		%%xmm14		, 0x060(%%rax)				\n\t"\
		"movq		%[__c12]		,%%rcx					\n\t"\
		"movaps		      (%%rdx)	,%%xmm4					\n\t		movaps		 0x080(%%rdx)	,%%xmm14				\n\t"\
		"movaps		 0x010(%%rdx)	,%%xmm5					\n\t		movaps		 0x090(%%rdx)	,%%xmm15				\n\t"\
		"movaps		%%xmm4		,%%xmm2						\n\t		movaps		%%xmm14		,%%xmm8 					\n\t"\
		"movaps		%%xmm5		,%%xmm3						\n\t		movaps		%%xmm15		,%%xmm9 					\n\t"\
		"mulpd		      (%%rcx)	,%%xmm4					\n\t		mulpd		 0x080(%%rcx)	,%%xmm14	/* c16 */	\n\t"\
		"mulpd		      (%%rcx)	,%%xmm5					\n\t		mulpd		 0x080(%%rcx)	,%%xmm15				\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm2					\n\t		mulpd		 0x090(%%rcx)	,%%xmm8 				\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm3					\n\t		mulpd		 0x090(%%rcx)	,%%xmm9 				\n\t"\
		"subpd		%%xmm2		,%%xmm5						\n\t		subpd		%%xmm8 		,%%xmm15					\n\t"\
		"addpd		%%xmm3		,%%xmm4						\n\t		addpd		%%xmm9 		,%%xmm14					\n\t"\
		"movaps		%%xmm5		,%%xmm3						\n\t		movaps		%%xmm15		,%%xmm9 					\n\t"\
		"movaps		%%xmm4		,%%xmm2						\n\t		movaps		%%xmm14		,%%xmm8 					\n\t"\
		"unpckhpd	 0x130(%%rbx)	,%%xmm3					\n\t		unpckhpd	 0x170(%%rbx)	,%%xmm9 				\n\t"\
		"unpcklpd	 0x130(%%rbx)	,%%xmm5					\n\t		unpcklpd	 0x170(%%rbx)	,%%xmm15				\n\t"\
		"movaps		%%xmm3		, 0x130(%%rbx)				\n\t		movaps		%%xmm9 		, 0x170(%%rbx)				\n\t"\
		"unpckhpd	 0x120(%%rbx)	,%%xmm2					\n\t		unpckhpd	 0x160(%%rbx)	,%%xmm8 				\n\t"\
		"unpcklpd	 0x120(%%rbx)	,%%xmm4					\n\t		unpcklpd	 0x160(%%rbx)	,%%xmm14				\n\t"\
		"movaps		%%xmm2		, 0x120(%%rbx)				\n\t		movaps		%%xmm8 		, 0x160(%%rbx)				\n\t"\
		"movaps		%%xmm5		, 0x130(%%rax)				\n\t		movaps		%%xmm15		, 0x170(%%rax)				\n\t"\
		"movaps		%%xmm4		, 0x120(%%rax)				\n\t		movaps		%%xmm14		, 0x160(%%rax)				\n\t"\
		"movq		%[__c0A]		,%%rcx					\n\t"\
		"subpd		%%xmm7		,%%xmm0						\n\t		subpd		%%xmm13		,%%xmm10					\n\t"\
		"subpd		%%xmm6		,%%xmm1						\n\t		subpd		%%xmm12		,%%xmm11					\n\t"\
		"addpd		%%xmm7		,%%xmm7						\n\t		addpd		%%xmm13		,%%xmm13					\n\t"\
		"addpd		%%xmm6		,%%xmm6						\n\t		addpd		%%xmm12		,%%xmm12					\n\t"\
		"addpd		%%xmm0		,%%xmm7						\n\t		addpd		%%xmm10		,%%xmm13					\n\t"\
		"addpd		%%xmm1		,%%xmm6						\n\t		addpd		%%xmm11		,%%xmm12					\n\t"\
		"movaps		%%xmm7		,%%xmm4						\n\t		movaps		%%xmm13		,%%xmm8 					\n\t"\
		"movaps		%%xmm1		,%%xmm5						\n\t		movaps		%%xmm11		,%%xmm9 					\n\t"\
		"mulpd		      (%%rcx)	,%%xmm7					\n\t		mulpd		 0x080(%%rcx)	,%%xmm13	/* c0E */	\n\t"\
		"mulpd		      (%%rcx)	,%%xmm1					\n\t		mulpd		 0x080(%%rcx)	,%%xmm11				\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm4					\n\t		mulpd		 0x090(%%rcx)	,%%xmm8 				\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm5					\n\t		mulpd		 0x090(%%rcx)	,%%xmm9 				\n\t"\
		"subpd		%%xmm4		,%%xmm1						\n\t		subpd		%%xmm8 		,%%xmm11					\n\t"\
		"addpd		%%xmm5		,%%xmm7						\n\t		addpd		%%xmm9 		,%%xmm13					\n\t"\
		"movaps		%%xmm1		,%%xmm5						\n\t		movaps		%%xmm11		,%%xmm9 					\n\t"\
		"movaps		%%xmm7		,%%xmm4						\n\t		movaps		%%xmm13		,%%xmm8 					\n\t"\
		"unpckhpd	 0x0b0(%%rbx)	,%%xmm5					\n\t		unpckhpd	 0x0f0(%%rbx)	,%%xmm9 				\n\t"\
		"unpcklpd	 0x0b0(%%rbx)	,%%xmm1					\n\t		unpcklpd	 0x0f0(%%rbx)	,%%xmm11				\n\t"\
		"movaps		%%xmm5		, 0x0b0(%%rbx)				\n\t		movaps		%%xmm9 		, 0x0f0(%%rbx)				\n\t"\
		"unpckhpd	 0x0a0(%%rbx)	,%%xmm4					\n\t		unpckhpd	 0x0e0(%%rbx)	,%%xmm8 				\n\t"\
		"unpcklpd	 0x0a0(%%rbx)	,%%xmm7					\n\t		unpcklpd	 0x0e0(%%rbx)	,%%xmm13				\n\t"\
		"movaps		%%xmm4		, 0x0a0(%%rbx)				\n\t		movaps		%%xmm8 		, 0x0e0(%%rbx)				\n\t"\
		"movaps		%%xmm1		, 0x0b0(%%rax)				\n\t		movaps		%%xmm11		, 0x0f0(%%rax)				\n\t"\
		"movaps		%%xmm7		, 0x0a0(%%rax)				\n\t		movaps		%%xmm13		, 0x0e0(%%rax)				\n\t"\
		"movq		%[__c1A]		,%%rcx					\n\t"\
		"movaps		%%xmm0		,%%xmm4						\n\t		movaps		%%xmm10		,%%xmm8 					\n\t"\
		"movaps		%%xmm6		,%%xmm5						\n\t		movaps		%%xmm12		,%%xmm9 					\n\t"\
		"mulpd		      (%%rcx)	,%%xmm0					\n\t		mulpd		 0x080(%%rcx)	,%%xmm10	/* c1E */	\n\t"\
		"mulpd		      (%%rcx)	,%%xmm6					\n\t		mulpd		 0x080(%%rcx)	,%%xmm12				\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm4					\n\t		mulpd		 0x090(%%rcx)	,%%xmm8 				\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm5					\n\t		mulpd		 0x090(%%rcx)	,%%xmm9 				\n\t"\
		"subpd		%%xmm4		,%%xmm6						\n\t		subpd		%%xmm8 		,%%xmm12					\n\t"\
		"addpd		%%xmm5		,%%xmm0						\n\t		addpd		%%xmm9 		,%%xmm10					\n\t"\
		"movaps		%%xmm6		,%%xmm5						\n\t		movaps		%%xmm12		,%%xmm9 					\n\t"\
		"movaps		%%xmm0		,%%xmm4						\n\t		movaps		%%xmm10		,%%xmm8 					\n\t"\
		"unpckhpd	 0x1b0(%%rbx)	,%%xmm5					\n\t		unpckhpd	 0x1f0(%%rbx)	,%%xmm9 				\n\t"\
		"unpcklpd	 0x1b0(%%rbx)	,%%xmm6					\n\t		unpcklpd	 0x1f0(%%rbx)	,%%xmm12				\n\t"\
		"movaps		%%xmm5		, 0x1b0(%%rbx)				\n\t		movaps		%%xmm9 		, 0x1f0(%%rbx)				\n\t"\
		"unpckhpd	 0x1a0(%%rbx)	,%%xmm4					\n\t		unpckhpd	 0x1e0(%%rbx)	,%%xmm8 				\n\t"\
		"unpcklpd	 0x1a0(%%rbx)	,%%xmm0					\n\t		unpcklpd	 0x1e0(%%rbx)	,%%xmm10				\n\t"\
		"movaps		%%xmm4		, 0x1a0(%%rbx)				\n\t		movaps		%%xmm8 		, 0x1e0(%%rbx)				\n\t"\
		"movaps		%%xmm6		, 0x1b0(%%rax)				\n\t		movaps		%%xmm12		, 0x1f0(%%rax)				\n\t"\
		"movaps		%%xmm0		, 0x1a0(%%rax)				\n\t		movaps		%%xmm10		, 0x1e0(%%rax)				\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__add1] "m" (Xadd1)\
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__r00] "m" (Xr00)\
		 ,[__r08] "m" (Xr08)\
		 ,[__r10] "m" (Xr10)\
		 ,[__r20] "m" (Xr20)\
		 ,[__r28] "m" (Xr28)\
		 ,[__r30] "m" (Xr30)\
		 ,[__c01] "m" (Xc01)\
		 ,[__c02] "m" (Xc02)\
		 ,[__c04] "m" (Xc04)\
		 ,[__c06] "m" (Xc06)\
		 ,[__c08] "m" (Xc08)\
		 ,[__c0A] "m" (Xc0A)\
		 ,[__c0C] "m" (Xc0C)\
		 ,[__c0E] "m" (Xc0E)\
		 ,[__c10] "m" (Xc10)\
		 ,[__c12] "m" (Xc12)\
		 ,[__c14] "m" (Xc14)\
		 ,[__c16] "m" (Xc16)\
		 ,[__c18] "m" (Xc18)\
		 ,[__c1A] "m" (Xc1A)\
		 ,[__c1C] "m" (Xc1C)\
		 ,[__c1E] "m" (Xc1E)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"		/* Clobbered registers */\
	);\
	}

	#define BAR(Xadd0,Xadd1,Xisrt2,Xr00,Xr08,Xr10,Xr20,Xr28,Xr30,Xc00,Xc01,Xc02,Xc03,Xc04,Xc05,Xc06,Xc07,Xc08,Xc0A,Xc0C,Xc0E,Xc10,Xc12,Xc14,Xc16,Xc18,Xc1A,Xc1C,Xc1E)\
	{\
	__asm__ volatile (\
		"movq		%[__add1]		,%%rbx\n\t"\
/*==========================*/"\n\t"\
"\n\t"/* AVX debug: everywhere we have an unpckhpd/unpcklpd pair followed by a pair of writes-to-memory, copy the same data into local memory: */\
"\n\t"/* add0,1 in rax,rbx; __r00 in rdx: */\
		"/*movaps		%%xmm3		,%%xmm7*/\n\t"\
		"/*movaps		%%xmm2		,%%xmm6*/\n\t"\
		"/*unpckhpd	 0x110(%%rbx)	,%%xmm7*/\n\t"\
		"/*unpcklpd	 0x110(%%rbx)	,%%xmm3*/\n\t"\
		"/*movaps		%%xmm7		, 0x110(%%rbx)*/\n\t"\
		"/*unpckhpd	 0x100(%%rbx)	,%%xmm6*/\n\t"\
		"/*unpcklpd	 0x100(%%rbx)	,%%xmm2*/\n\t"\
		"/*movaps		%%xmm6		, 0x100(%%rbx)*/\n\t"\
		"/*movaps		%%xmm3		, 0x110(%%rax)*/\n\t"\
		"/*movaps		%%xmm2		, 0x100(%%rax)*/\n\t"\
"\n\t"/* For each complex output octet, the complex pairs involving memory reads from offsets 0x00..,0x08..,0x10..,0x18.. go into local-mem pairs rXY+00/10,04/14,02/12,06/16: */\
"\n\t"/* For the first output octet we do memory reads from offsets [0x10..,0x00..],[0x08..,0x18], the other 3 octets use order [0x00..,0x10..],[0x08..,0x18]: */\
"movaps	0x110(%%rbx),%%xmm7		\n\t"\
"movaps	0x100(%%rbx),%%xmm6		\n\t"\
"movaps	%%xmm3,0x030(%%rdx)		\n\t"\
"movaps	%%xmm2,0x020(%%rdx)		\n\t"/* r02,03 */\
"movaps	%%xmm7,0x130(%%rdx)		\n\t"\
"movaps	%%xmm6,0x120(%%rdx)		\n\t"/* r12,13 */\
"\n\t"\
		"movaps		 0x010(%%rdx)	,%%xmm3\n\t"\
		"movaps		      (%%rdx)	,%%xmm2\n\t"\
		"/*movaps		%%xmm3		,%%xmm7*/\n\t"\
		"/*movaps		%%xmm2		,%%xmm6*/\n\t"\
		"/*unpckhpd	 0x010(%%rbx)	,%%xmm7*/\n\t"\
		"/*unpcklpd	 0x010(%%rbx)	,%%xmm3*/\n\t"\
		"/*movaps		%%xmm7		, 0x010(%%rbx)*/\n\t"\
		"/*unpckhpd	      (%%rbx)	,%%xmm6*/\n\t"\
		"/*unpcklpd	      (%%rbx)	,%%xmm2*/\n\t"\
		"/*movaps		%%xmm6		,      (%%rbx)*/\n\t"\
		"/*movaps		%%xmm3		, 0x010(%%rax)*/\n\t"\
		"/*movaps		%%xmm2		,      (%%rax)*/\n\t"\
"movaps	0x010(%%rbx),%%xmm7		\n\t"\
"movaps	0x000(%%rbx),%%xmm6		\n\t"\
"movaps	%%xmm3,0x010(%%rdx)		\n\t"\
"movaps	%%xmm2,     (%%rdx)		\n\t"/* r00,01 */\
"movaps	%%xmm7,0x110(%%rdx)		\n\t"\
"movaps	%%xmm6,0x100(%%rdx)		\n\t"/* r10,11 */\
"\n\t"\
		"movq		%[__c08]		,%%rcx\n\t"\
		"addpd		%%xmm5		,%%xmm0\n\t"\
		"subpd		%%xmm4		,%%xmm1\n\t"\
		"movaps		%%xmm0		,%%xmm2\n\t"\
		"movaps		%%xmm1		,%%xmm3\n\t"\
		"addpd		%%xmm5		,%%xmm5\n\t"\
		"addpd		%%xmm4		,%%xmm4\n\t"\
		"movaps		%%xmm0		,%%xmm6\n\t"\
		"movaps		%%xmm1		,%%xmm7\n\t"\
		"mulpd		      (%%rcx)	,%%xmm2\n\t"\
		"mulpd		      (%%rcx)	,%%xmm3\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm6\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm7\n\t"\
		"subpd		%%xmm6		,%%xmm3\n\t"\
		"addpd		%%xmm7		,%%xmm2\n\t"\
		"/*movaps		%%xmm3		,%%xmm7*/\n\t"\
		"/*movaps		%%xmm2		,%%xmm6*/\n\t"\
		"/*unpckhpd	 0x090(%%rbx)	,%%xmm7*/\n\t"\
		"/*unpcklpd	 0x090(%%rbx)	,%%xmm3*/\n\t"\
		"/*movaps		%%xmm7		, 0x090(%%rbx)*/\n\t"\
		"/*unpckhpd	 0x080(%%rbx)	,%%xmm6*/\n\t"\
		"/*unpcklpd	 0x080(%%rbx)	,%%xmm2*/\n\t"\
		"/*movaps		%%xmm6		, 0x080(%%rbx)*/\n\t"\
		"/*movaps		%%xmm3		, 0x090(%%rax)*/\n\t"\
		"/*movaps		%%xmm2		, 0x080(%%rax)*/\n\t"\
"movaps	0x090(%%rbx),%%xmm7		\n\t"\
"movaps	0x080(%%rbx),%%xmm6		\n\t"\
"movaps	%%xmm3,0x050(%%rdx)		\n\t"\
"movaps	%%xmm2,0x040(%%rdx)		\n\t"/* r04,05 */\
"movaps	%%xmm7,0x150(%%rdx)		\n\t"\
"movaps	%%xmm6,0x140(%%rdx)		\n\t"/* r14,15 */\
"\n\t"\
		"movq		%[__c18]		,%%rcx\n\t"\
		"subpd		%%xmm5		,%%xmm0\n\t"\
		"addpd		%%xmm4		,%%xmm1\n\t"\
		"movaps		%%xmm0		,%%xmm6\n\t"\
		"movaps		%%xmm1		,%%xmm7\n\t"\
		"mulpd		      (%%rcx)	,%%xmm0\n\t"\
		"mulpd		      (%%rcx)	,%%xmm1\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm6\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm7\n\t"\
		"subpd		%%xmm6		,%%xmm1\n\t"\
		"addpd		%%xmm7		,%%xmm0\n\t"\
		"/*movaps		%%xmm1		,%%xmm7*/\n\t"\
		"/*movaps		%%xmm0		,%%xmm6*/\n\t"\
		"/*unpckhpd	 0x190(%%rbx)	,%%xmm7*/\n\t"\
		"/*unpcklpd	 0x190(%%rbx)	,%%xmm1*/\n\t"\
		"/*movaps		%%xmm7		, 0x190(%%rbx)*/\n\t"\
		"/*unpckhpd	 0x180(%%rbx)	,%%xmm6*/\n\t"\
		"/*unpcklpd	 0x180(%%rbx)	,%%xmm0*/\n\t"\
		"/*movaps		%%xmm6		, 0x180(%%rbx)*/\n\t"\
		"/*movaps		%%xmm1		, 0x190(%%rax)*/\n\t"\
		"/*movaps		%%xmm0		, 0x180(%%rax)*/\n\t"\
"movaps	0x190(%%rbx),%%xmm7		\n\t"\
"movaps	0x180(%%rbx),%%xmm6		\n\t"\
"movaps	%%xmm1,0x070(%%rdx)		\n\t"\
"movaps	%%xmm0,0x060(%%rdx)		\n\t"/* r06,07 */\
"movaps	%%xmm7,0x170(%%rdx)		\n\t"\
"movaps	%%xmm6,0x160(%%rdx)		\n\t"/* r16,17 */\
"\n\t"\
		"/************************************************/\n\t"\
		"/*...Block 5: t08,t18,t28,t38 -> r08,0C,0A,0E	*/\n\t"\
		"/************************************************/\n\t"\
		"movq		%[__r08]		,%%rdx\n\t"\
		"movaps		 (%%rsi)	,%%xmm2\n\t	/* isrt2 */"\
		"movaps		 0x020(%%rdx)	,%%xmm4\n\t"\
		"movaps		 0x030(%%rdx)	,%%xmm5\n\t"\
		"movaps		 0x060(%%rdx)	,%%xmm0\n\t"\
		"movaps		 0x070(%%rdx)	,%%xmm1\n\t"\
		"addpd		 0x030(%%rdx)	,%%xmm4\n\t"\
		"subpd		 0x020(%%rdx)	,%%xmm5\n\t"\
		"subpd		 0x070(%%rdx)	,%%xmm0\n\t"\
		"addpd		 0x060(%%rdx)	,%%xmm1\n\t"\
		"mulpd		%%xmm2		,%%xmm4\n\t"\
		"mulpd		%%xmm2		,%%xmm5\n\t"\
		"mulpd		%%xmm2		,%%xmm0\n\t"\
		"mulpd		%%xmm2		,%%xmm1\n\t"\
		"movaps		%%xmm4		,%%xmm6\n\t"\
		"movaps		%%xmm5		,%%xmm7\n\t"\
		"subpd		%%xmm0		,%%xmm4\n\t"\
		"subpd		%%xmm1		,%%xmm5\n\t"\
		"addpd		%%xmm0		,%%xmm6\n\t"\
		"addpd		%%xmm1		,%%xmm7\n\t"\
		"movaps		      (%%rdx)	,%%xmm0\n\t"\
		"movaps		 0x010(%%rdx)	,%%xmm1\n\t"\
		"movaps		 0x040(%%rdx)	,%%xmm2\n\t"\
		"movaps		 0x050(%%rdx)	,%%xmm3\n\t"\
		"subpd		 0x050(%%rdx)	,%%xmm0\n\t"\
		"subpd		 0x040(%%rdx)	,%%xmm1\n\t"\
		"addpd		      (%%rdx)	,%%xmm3\n\t"\
		"addpd		 0x010(%%rdx)	,%%xmm2\n\t"\
		"movq		%[__add0]		,%%rax\n\t"\
		"movq		%[__c04]		,%%rcx\n\t"\
		"subpd		%%xmm4		,%%xmm3\n\t"\
		"subpd		%%xmm5		,%%xmm1\n\t"\
		"addpd		%%xmm4		,%%xmm4\n\t"\
		"addpd		%%xmm5		,%%xmm5\n\t"\
		"addpd		%%xmm3		,%%xmm4\n\t"\
		"addpd		%%xmm1		,%%xmm5\n\t"\
	"\n\t"/* Switch these spills from rdx + 0,1 to + 2,3 to avoid them being overwritten by the r08,09 tmp-stores in the debug below: */\
		"movaps		%%xmm3		, 0x020(%%rdx)\n\t"\
		"movaps		%%xmm1		, 0x030(%%rdx)\n\t"\
		"movaps		%%xmm4		,%%xmm3\n\t"\
		"movaps		%%xmm5		,%%xmm1\n\t"\
		"mulpd		      (%%rcx)	,%%xmm4\n\t"\
		"mulpd		      (%%rcx)	,%%xmm5\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm3\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm1\n\t"\
		"subpd		%%xmm3		,%%xmm5\n\t"\
		"addpd		%%xmm1		,%%xmm4\n\t"\
		"movq		%[__add1]		,%%rbx\n\t"\
		"/*movaps		%%xmm5		,%%xmm3*/\n\t"\
		"/*movaps		%%xmm4		,%%xmm1*/\n\t"\
		"/*unpckhpd	 0x050(%%rbx)	,%%xmm3*/\n\t"\
		"/*unpcklpd	 0x050(%%rbx)	,%%xmm5*/\n\t"\
		"/*movaps		%%xmm3		, 0x050(%%rbx)*/\n\t"\
		"/*unpckhpd	 0x040(%%rbx)	,%%xmm1*/\n\t"\
		"/*unpcklpd	 0x040(%%rbx)	,%%xmm4*/\n\t"\
		"/*movaps		%%xmm1		, 0x040(%%rbx)*/\n\t"\
		"/*movaps		%%xmm5		, 0x050(%%rax)*/\n\t"\
		"/*movaps		%%xmm4		, 0x040(%%rax)*/\n\t"\
"movaps	0x050(%%rbx),%%xmm3		\n\t"\
"movaps	0x040(%%rbx),%%xmm1		\n\t"\
"movaps	%%xmm5,0x030(%%rdx)		\n\t"\
"movaps	%%xmm4,0x020(%%rdx)		\n\t"/* r0a,0b */\
"movaps	%%xmm3,0x130(%%rdx)		\n\t"\
"movaps	%%xmm1,0x120(%%rdx)		\n\t"/* r1a,1b */\
"\n\t"\
		"movq		%[__c14]		,%%rcx\n\t"\
		"movaps		 0x020(%%rdx)	,%%xmm4\n\t"\
		"movaps		 0x030(%%rdx)	,%%xmm5\n\t"\
		"movaps		%%xmm4		,%%xmm3\n\t"\
		"movaps		%%xmm5		,%%xmm1\n\t"\
		"mulpd		      (%%rcx)	,%%xmm4\n\t"\
		"mulpd		      (%%rcx)	,%%xmm5\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm3\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm1\n\t"\
		"subpd		%%xmm3		,%%xmm5\n\t"\
		"addpd		%%xmm1		,%%xmm4\n\t"\
		"/*movaps		%%xmm5		,%%xmm3*/\n\t"\
		"/*movaps		%%xmm4		,%%xmm1*/\n\t"\
		"/*unpckhpd	 0x150(%%rbx)	,%%xmm3*/\n\t"\
		"/*unpcklpd	 0x150(%%rbx)	,%%xmm5*/\n\t"\
		"/*movaps		%%xmm3		, 0x150(%%rbx)*/\n\t"\
		"/*unpckhpd	 0x140(%%rbx)	,%%xmm1*/\n\t"\
		"/*unpcklpd	 0x140(%%rbx)	,%%xmm4*/\n\t"\
		"/*movaps		%%xmm1		, 0x140(%%rbx)*/\n\t"\
		"/*movaps		%%xmm5		, 0x150(%%rax)*/\n\t"\
		"/*movaps		%%xmm4		, 0x140(%%rax)*/\n\t"\
"movaps	0x150(%%rbx),%%xmm3		\n\t"\
"movaps	0x140(%%rbx),%%xmm1		\n\t"\
"movaps	%%xmm5,0x010(%%rdx)		\n\t"\
"movaps	%%xmm4,     (%%rdx)		\n\t"/* r08,09 */\
"movaps	%%xmm3,0x110(%%rdx)		\n\t"\
"movaps	%%xmm1,0x100(%%rdx)		\n\t"/* r18,19 */\
"\n\t"\
		"movq		%[__c0C]		,%%rcx\n\t"\
		"subpd		%%xmm7		,%%xmm0\n\t"\
		"subpd		%%xmm6		,%%xmm2\n\t"\
		"addpd		%%xmm7		,%%xmm7\n\t"\
		"addpd		%%xmm6		,%%xmm6\n\t"\
		"addpd		%%xmm0		,%%xmm7\n\t"\
		"addpd		%%xmm2		,%%xmm6\n\t"\
		"movaps		%%xmm7		,%%xmm4\n\t"\
		"movaps		%%xmm2		,%%xmm5\n\t"\
		"mulpd		      (%%rcx)	,%%xmm7\n\t"\
		"mulpd		      (%%rcx)	,%%xmm2\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm4\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm5\n\t"\
		"subpd		%%xmm4		,%%xmm2\n\t"\
		"addpd		%%xmm5		,%%xmm7\n\t"\
		"/*movaps		%%xmm2		,%%xmm5*/\n\t"\
		"/*movaps		%%xmm7		,%%xmm4*/\n\t"\
		"/*unpckhpd	 0x0d0(%%rbx)	,%%xmm5*/\n\t"\
		"/*unpcklpd	 0x0d0(%%rbx)	,%%xmm2*/\n\t"\
		"/*movaps		%%xmm5		, 0x0d0(%%rbx)*/\n\t"\
		"/*unpckhpd	 0x0c0(%%rbx)	,%%xmm4*/\n\t"\
		"/*unpcklpd	 0x0c0(%%rbx)	,%%xmm7*/\n\t"\
		"/*movaps		%%xmm4		, 0x0c0(%%rbx)*/\n\t"\
		"/*movaps		%%xmm2		, 0x0d0(%%rax)*/\n\t"\
		"/*movaps		%%xmm7		, 0x0c0(%%rax)*/\n\t"\
"movaps	0x0d0(%%rbx),%%xmm5		\n\t"\
"movaps	0x0c0(%%rbx),%%xmm4		\n\t"\
"movaps	%%xmm2,0x050(%%rdx)		\n\t"\
"movaps	%%xmm7,0x040(%%rdx)		\n\t"/* r0c,0d */\
"movaps	%%xmm5,0x150(%%rdx)		\n\t"\
"movaps	%%xmm4,0x140(%%rdx)		\n\t"/* r1c,1d */\
"\n\t"\
		"movq		%[__c1C]		,%%rcx\n\t"\
		"movaps		%%xmm0		,%%xmm4\n\t"\
		"movaps		%%xmm6		,%%xmm5\n\t"\
		"mulpd		      (%%rcx)	,%%xmm0\n\t"\
		"mulpd		      (%%rcx)	,%%xmm6\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm4\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm5\n\t"\
		"subpd		%%xmm4		,%%xmm6\n\t"\
		"addpd		%%xmm5		,%%xmm0\n\t"\
		"/*movaps		%%xmm6		,%%xmm5*/\n\t"\
		"/*movaps		%%xmm0		,%%xmm4*/\n\t"\
		"/*unpckhpd	 0x1d0(%%rbx)	,%%xmm5*/\n\t"\
		"/*unpcklpd	 0x1d0(%%rbx)	,%%xmm6*/\n\t"\
		"/*movaps		%%xmm5		, 0x1d0(%%rbx)*/\n\t"\
		"/*unpckhpd	 0x1c0(%%rbx)	,%%xmm4*/\n\t"\
		"/*unpcklpd	 0x1c0(%%rbx)	,%%xmm0*/\n\t"\
		"/*movaps		%%xmm4		, 0x1c0(%%rbx)*/\n\t"\
		"/*movaps		%%xmm6		, 0x1d0(%%rax)*/\n\t"\
		"/*movaps		%%xmm0		, 0x1c0(%%rax)*/\n\t"\
"movaps	0x1d0(%%rbx),%%xmm5		\n\t"\
"movaps	0x1c0(%%rbx),%%xmm4		\n\t"\
"movaps	%%xmm6,0x070(%%rdx)		\n\t"\
"movaps	%%xmm0,0x060(%%rdx)		\n\t"/* r0e,0f */\
"movaps	%%xmm5,0x170(%%rdx)		\n\t"\
"movaps	%%xmm4,0x160(%%rdx)		\n\t"/* r1e,1f */\
"\n\t"\
		"/************************************************/\n\t"\
		"/*...Block 3: t04,t14,t24,t34 -> r20,24,22,26	*/\n\t"\
		"/************************************************/\n\t"\
		"movq		%[__r20]		,%%rdx\n\t"\
		"movq		%%rsi			,%%rcx\n\t"\
		"addq		$0x010			,%%rcx	/* cc0 */\n\t"\
		"movaps		 0x020(%%rdx)	,%%xmm4\n\t"\
		"movaps		 0x060(%%rdx)	,%%xmm0\n\t"\
		"movaps		 0x030(%%rdx)	,%%xmm5\n\t"\
		"movaps		 0x070(%%rdx)	,%%xmm1\n\t"\
		"movaps		 0x020(%%rdx)	,%%xmm6\n\t"\
		"movaps		 0x060(%%rdx)	,%%xmm2\n\t"\
		"movaps		 0x030(%%rdx)	,%%xmm7\n\t"\
		"movaps		 0x070(%%rdx)	,%%xmm3\n\t"\
		"mulpd		      (%%rcx)	,%%xmm4\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm0\n\t"\
		"mulpd		      (%%rcx)	,%%xmm5\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm1\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm6\n\t"\
		"mulpd		      (%%rcx)	,%%xmm2\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm7\n\t"\
		"mulpd		      (%%rcx)	,%%xmm3\n\t"\
		"subpd		%%xmm6		,%%xmm5\n\t"\
		"subpd		%%xmm2		,%%xmm1\n\t"\
		"addpd		%%xmm7		,%%xmm4\n\t"\
		"addpd		%%xmm3		,%%xmm0\n\t"\
		"movaps		%%xmm5		,%%xmm7\n\t"\
		"movaps		%%xmm4		,%%xmm6\n\t"\
		"addpd		%%xmm0		,%%xmm4\n\t"\
		"addpd		%%xmm1		,%%xmm5\n\t"\
		"subpd		%%xmm0		,%%xmm6\n\t"\
		"subpd		%%xmm1		,%%xmm7\n\t"\
		"movaps		 0x040(%%rdx)	,%%xmm2\n\t"\
		"movaps		 0x050(%%rdx)	,%%xmm3\n\t"\
		"movaps		      (%%rdx)	,%%xmm0\n\t"\
		"movaps		 0x010(%%rdx)	,%%xmm1\n\t"\
		"addpd		 0x050(%%rdx)	,%%xmm2\n\t"\
		"subpd		 0x040(%%rdx)	,%%xmm3\n\t"\
		"mulpd		      (%%rsi)	,%%xmm2\n\t"\
		"mulpd		      (%%rsi)	,%%xmm3\n\t"\
		"subpd		%%xmm2		,%%xmm0\n\t"\
		"subpd		%%xmm3		,%%xmm1\n\t"\
		"addpd		%%xmm2		,%%xmm2\n\t"\
		"addpd		%%xmm3		,%%xmm3\n\t"\
		"addpd		%%xmm0		,%%xmm2\n\t"\
		"addpd		%%xmm1		,%%xmm3\n\t"\
		"movq		%[__add0]		,%%rax\n\t"\
		"movq		%[__c02]		,%%rcx\n\t"\
		"subpd		%%xmm4		,%%xmm2\n\t"\
		"subpd		%%xmm5		,%%xmm3\n\t"\
		"addpd		%%xmm4		,%%xmm4\n\t"\
		"addpd		%%xmm5		,%%xmm5\n\t"\
		"addpd		%%xmm2		,%%xmm4\n\t"\
		"addpd		%%xmm3		,%%xmm5\n\t"\
	"\n\t"/* Switch these spills from rdx + 0,1 to + 2,3 to avoid them being overwritten by the r20,21 tmp-stores in the debug below: */\
		"movaps		%%xmm2		, 0x020(%%rdx)\n\t"\
		"movaps		%%xmm3		, 0x030(%%rdx)\n\t"\
		"movaps		%%xmm4		,%%xmm2\n\t"\
		"movaps		%%xmm5		,%%xmm3\n\t"\
		"mulpd		      (%%rcx)	,%%xmm4\n\t"\
		"mulpd		      (%%rcx)	,%%xmm5\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm2\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm3\n\t"\
		"subpd		%%xmm2		,%%xmm5\n\t"\
		"addpd		%%xmm3		,%%xmm4\n\t"\
		"movq		%[__add1]		,%%rbx\n\t"\
		"/*movaps		%%xmm5		,%%xmm3*/\n\t"\
		"/*movaps		%%xmm4		,%%xmm2*/\n\t"\
		"/*unpckhpd	 0x030(%%rbx)	,%%xmm3*/\n\t"\
		"/*unpcklpd	 0x030(%%rbx)	,%%xmm5*/\n\t"\
		"/*movaps		%%xmm3		, 0x030(%%rbx)*/\n\t"\
		"/*unpckhpd	 0x020(%%rbx)	,%%xmm2*/\n\t"\
		"/*unpcklpd	 0x020(%%rbx)	,%%xmm4*/\n\t"\
		"/*movaps		%%xmm2		, 0x020(%%rbx)*/\n\t"\
		"/*movaps		%%xmm5		, 0x030(%%rax)*/\n\t"\
		"/*movaps		%%xmm4		, 0x020(%%rax)*/\n\t"\
"movaps	0x030(%%rbx),%%xmm3		\n\t"\
"movaps	0x020(%%rbx),%%xmm2		\n\t"\
"movaps	%%xmm5,0x030(%%rdx)		\n\t"\
"movaps	%%xmm4,0x020(%%rdx)		\n\t"/* r22,23 */\
"movaps	%%xmm3,0x130(%%rdx)		\n\t"\
"movaps	%%xmm2,0x120(%%rdx)		\n\t"/* r32,33 */\
"\n\t"\
		"movq		%[__c12]		,%%rcx\n\t"\
		"movaps		 0x020(%%rdx)	,%%xmm4\n\t"\
		"movaps		 0x030(%%rdx)	,%%xmm5\n\t"\
		"movaps		%%xmm4		,%%xmm2\n\t"\
		"movaps		%%xmm5		,%%xmm3\n\t"\
		"mulpd		      (%%rcx)	,%%xmm4\n\t"\
		"mulpd		      (%%rcx)	,%%xmm5\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm2\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm3\n\t"\
		"subpd		%%xmm2		,%%xmm5\n\t"\
		"addpd		%%xmm3		,%%xmm4\n\t"\
		"/*movaps		%%xmm5		,%%xmm3*/\n\t"\
		"/*movaps		%%xmm4		,%%xmm2*/\n\t"\
		"/*unpckhpd	 0x130(%%rbx)	,%%xmm3*/\n\t"\
		"/*unpcklpd	 0x130(%%rbx)	,%%xmm5*/\n\t"\
		"/*movaps		%%xmm3		, 0x130(%%rbx)*/\n\t"\
		"/*unpckhpd	 0x120(%%rbx)	,%%xmm2*/\n\t"\
		"/*unpcklpd	 0x120(%%rbx)	,%%xmm4*/\n\t"\
		"/*movaps		%%xmm2		, 0x120(%%rbx)*/\n\t"\
		"/*movaps		%%xmm5		, 0x130(%%rax)*/\n\t"\
		"/*movaps		%%xmm4		, 0x120(%%rax)*/\n\t"\
"movaps	0x130(%%rbx),%%xmm3		\n\t"\
"movaps	0x120(%%rbx),%%xmm2		\n\t"\
"movaps	%%xmm5,0x010(%%rdx)		\n\t"\
"movaps	%%xmm4,     (%%rdx)		\n\t"/* r20,21 */\
"movaps	%%xmm3,0x110(%%rdx)		\n\t"\
"movaps	%%xmm2,0x100(%%rdx)		\n\t"/* r30,31 */\
"\n\t"\
		"movq		%[__c0A]		,%%rcx\n\t"\
		"subpd		%%xmm7		,%%xmm0\n\t"\
		"subpd		%%xmm6		,%%xmm1\n\t"\
		"addpd		%%xmm7		,%%xmm7\n\t"\
		"addpd		%%xmm6		,%%xmm6\n\t"\
		"addpd		%%xmm0		,%%xmm7\n\t"\
		"addpd		%%xmm1		,%%xmm6\n\t"\
		"movaps		%%xmm7		,%%xmm4\n\t"\
		"movaps		%%xmm1		,%%xmm5\n\t"\
		"mulpd		      (%%rcx)	,%%xmm7\n\t"\
		"mulpd		      (%%rcx)	,%%xmm1\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm4\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm5\n\t"\
		"subpd		%%xmm4		,%%xmm1\n\t"\
		"addpd		%%xmm5		,%%xmm7\n\t"\
		"/*movaps		%%xmm1		,%%xmm5*/\n\t"\
		"/*movaps		%%xmm7		,%%xmm4*/\n\t"\
		"/*unpckhpd	 0x0b0(%%rbx)	,%%xmm5*/\n\t"\
		"/*unpcklpd	 0x0b0(%%rbx)	,%%xmm1*/\n\t"\
		"/*movaps		%%xmm5		, 0x0b0(%%rbx)*/\n\t"\
		"/*unpckhpd	 0x0a0(%%rbx)	,%%xmm4*/\n\t"\
		"/*unpcklpd	 0x0a0(%%rbx)	,%%xmm7*/\n\t"\
		"/*movaps		%%xmm4		, 0x0a0(%%rbx)*/\n\t"\
		"/*movaps		%%xmm1		, 0x0b0(%%rax)*/\n\t"\
		"/*movaps		%%xmm7		, 0x0a0(%%rax)*/\n\t"\
"movaps	0x0b0(%%rbx),%%xmm5		\n\t"\
"movaps	0x0a0(%%rbx),%%xmm4		\n\t"\
"movaps	%%xmm1,0x050(%%rdx)		\n\t"\
"movaps	%%xmm7,0x040(%%rdx)		\n\t"/* r24,25 */\
"movaps	%%xmm5,0x150(%%rdx)		\n\t"\
"movaps	%%xmm4,0x140(%%rdx)		\n\t"/* r34,35 */\
"\n\t"\
		"movq		%[__c1A]		,%%rcx\n\t"\
		"movaps		%%xmm0		,%%xmm4\n\t"\
		"movaps		%%xmm6		,%%xmm5\n\t"\
		"mulpd		      (%%rcx)	,%%xmm0\n\t"\
		"mulpd		      (%%rcx)	,%%xmm6\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm4\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm5\n\t"\
		"subpd		%%xmm4		,%%xmm6\n\t"\
		"addpd		%%xmm5		,%%xmm0\n\t"\
		"/*movaps		%%xmm6		,%%xmm5*/\n\t"\
		"/*movaps		%%xmm0		,%%xmm4*/\n\t"\
		"/*unpckhpd	 0x1b0(%%rbx)	,%%xmm5*/\n\t"\
		"/*unpcklpd	 0x1b0(%%rbx)	,%%xmm6*/\n\t"\
		"/*movaps		%%xmm5		, 0x1b0(%%rbx)*/\n\t"\
		"/*unpckhpd	 0x1a0(%%rbx)	,%%xmm4*/\n\t"\
		"/*unpcklpd	 0x1a0(%%rbx)	,%%xmm0*/\n\t"\
		"/*movaps		%%xmm4		, 0x1a0(%%rbx)*/\n\t"\
		"/*movaps		%%xmm6		, 0x1b0(%%rax)*/\n\t"\
		"/*movaps		%%xmm0		, 0x1a0(%%rax)*/\n\t"\
"movaps	0x1b0(%%rbx),%%xmm5		\n\t"\
"movaps	0x1a0(%%rbx),%%xmm4		\n\t"\
"movaps	%%xmm6,0x070(%%rdx)		\n\t"\
"movaps	%%xmm0,0x060(%%rdx)		\n\t"/* r26,27 */\
"movaps	%%xmm5,0x170(%%rdx)		\n\t"\
"movaps	%%xmm4,0x160(%%rdx)		\n\t"/* r36,37 */\
"\n\t"\
		"/************************************************/\n\t"\
		"/*...Block 7: t0C,t1C,t2C,t3C -> r28,2C,2A,2E	*/\n\t"\
		"/************************************************/\n\t"\
		"movq		%[__r28]		,%%rdx\n\t"\
		"movq		%%rsi			,%%rcx\n\t"\
		"addq		$0x010			,%%rcx	/* cc0 */\n\t"\
		"movaps		 0x020(%%rdx)	,%%xmm4\n\t"\
		"movaps		 0x060(%%rdx)	,%%xmm0\n\t"\
		"movaps		 0x030(%%rdx)	,%%xmm5\n\t"\
		"movaps		 0x070(%%rdx)	,%%xmm1\n\t"\
		"movaps		 0x020(%%rdx)	,%%xmm6\n\t"\
		"movaps		 0x060(%%rdx)	,%%xmm2\n\t"\
		"movaps		 0x030(%%rdx)	,%%xmm7\n\t"\
		"movaps		 0x070(%%rdx)	,%%xmm3\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm4\n\t"\
		"mulpd		      (%%rcx)	,%%xmm0\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm5\n\t"\
		"mulpd		      (%%rcx)	,%%xmm1\n\t"\
		"mulpd		      (%%rcx)	,%%xmm6\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm2\n\t"\
		"mulpd		      (%%rcx)	,%%xmm7\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm3\n\t"\
		"subpd		%%xmm6		,%%xmm5\n\t"\
		"subpd		%%xmm2		,%%xmm1\n\t"\
		"addpd		%%xmm7		,%%xmm4\n\t"\
		"addpd		%%xmm3		,%%xmm0\n\t"\
		"movaps		%%xmm5		,%%xmm7\n\t"\
		"movaps		%%xmm4		,%%xmm6\n\t"\
		"addpd		%%xmm0		,%%xmm4\n\t"\
		"addpd		%%xmm1		,%%xmm5\n\t"\
		"subpd		%%xmm0		,%%xmm6\n\t"\
		"subpd		%%xmm1		,%%xmm7\n\t"\
		"movaps		 0x040(%%rdx)	,%%xmm2\n\t"\
		"movaps		 0x050(%%rdx)	,%%xmm3\n\t"\
		"movaps		      (%%rdx)	,%%xmm0\n\t"\
		"movaps		 0x010(%%rdx)	,%%xmm1\n\t"\
		"subpd		 0x050(%%rdx)	,%%xmm2\n\t"\
		"addpd		 0x040(%%rdx)	,%%xmm3\n\t"\
		"mulpd		      (%%rsi)	,%%xmm2\n\t"\
		"mulpd		      (%%rsi)	,%%xmm3\n\t"\
		"subpd		%%xmm2		,%%xmm0\n\t"\
		"subpd		%%xmm3		,%%xmm1\n\t"\
		"addpd		%%xmm2		,%%xmm2\n\t"\
		"addpd		%%xmm3		,%%xmm3\n\t"\
		"addpd		%%xmm0		,%%xmm2\n\t"\
		"addpd		%%xmm1		,%%xmm3\n\t"\
		"movq		%[__add0]		,%%rax\n\t"\
		"movq		%[__c06]		,%%rcx\n\t"\
		"subpd		%%xmm6		,%%xmm0\n\t"\
		"subpd		%%xmm7		,%%xmm1\n\t"\
		"addpd		%%xmm6		,%%xmm6\n\t"\
		"addpd		%%xmm7		,%%xmm7\n\t"\
		"addpd		%%xmm0		,%%xmm6\n\t"\
		"addpd		%%xmm1		,%%xmm7\n\t"\
	"\n\t"/* Switch these spills from rdx + 0,1 to + 2,3 to avoid them being overwritten by the r28,29 tmp-stores in the debug below: */\
		"movaps		%%xmm0		, 0x020(%%rdx)\n\t"\
		"movaps		%%xmm1		, 0x030(%%rdx)\n\t"\
		"movaps		%%xmm6		,%%xmm0\n\t"\
		"movaps		%%xmm7		,%%xmm1\n\t"\
		"mulpd		      (%%rcx)	,%%xmm6\n\t"\
		"mulpd		      (%%rcx)	,%%xmm7\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm0\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm1\n\t"\
		"subpd		%%xmm0		,%%xmm7\n\t"\
		"addpd		%%xmm1		,%%xmm6\n\t"\
		"movq		%[__add1]		,%%rbx\n\t"\
		"/*movaps		%%xmm7		,%%xmm1*/\n\t"\
		"/*movaps		%%xmm6		,%%xmm0*/\n\t"\
		"/*unpckhpd	 0x070(%%rbx)	,%%xmm1*/\n\t"\
		"/*unpcklpd	 0x070(%%rbx)	,%%xmm7*/\n\t"\
		"/*movaps		%%xmm1		, 0x070(%%rbx)*/\n\t"\
		"/*unpckhpd	 0x060(%%rbx)	,%%xmm0*/\n\t"\
		"/*unpcklpd	 0x060(%%rbx)	,%%xmm6*/\n\t"\
		"/*movaps		%%xmm0		, 0x060(%%rbx)*/\n\t"\
		"/*movaps		%%xmm7		, 0x070(%%rax)*/\n\t"\
		"/*movaps		%%xmm6		, 0x060(%%rax)*/\n\t"\
"movaps	0x070(%%rbx),%%xmm1		\n\t"\
"movaps	0x060(%%rbx),%%xmm0		\n\t"\
"movaps	%%xmm7,0x030(%%rdx)		\n\t"\
"movaps	%%xmm6,0x020(%%rdx)		\n\t"/* r2a,2b */\
"movaps	%%xmm1,0x130(%%rdx)		\n\t"\
"movaps	%%xmm0,0x120(%%rdx)		\n\t"/* r3a,3b */\
"\n\t"\
		"movq		%[__c16]		,%%rcx\n\t"\
		"movaps		 0x020(%%rdx)	,%%xmm6\n\t"\
		"movaps		 0x030(%%rdx)	,%%xmm7\n\t"\
		"movaps		%%xmm6		,%%xmm0\n\t"\
		"movaps		%%xmm7		,%%xmm1\n\t"\
		"mulpd		      (%%rcx)	,%%xmm6\n\t"\
		"mulpd		      (%%rcx)	,%%xmm7\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm0\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm1\n\t"\
		"subpd		%%xmm0		,%%xmm7\n\t"\
		"addpd		%%xmm1		,%%xmm6\n\t"\
		"/*movaps		%%xmm7		,%%xmm1*/\n\t"\
		"/*movaps		%%xmm6		,%%xmm0*/\n\t"\
		"/*unpckhpd	 0x170(%%rbx)	,%%xmm1*/\n\t"\
		"/*unpcklpd	 0x170(%%rbx)	,%%xmm7*/\n\t"\
		"/*movaps		%%xmm1		, 0x170(%%rbx)*/\n\t"\
		"/*unpckhpd	 0x160(%%rbx)	,%%xmm0*/\n\t"\
		"/*unpcklpd	 0x160(%%rbx)	,%%xmm6*/\n\t"\
		"/*movaps		%%xmm0		, 0x160(%%rbx)*/\n\t"\
		"/*movaps		%%xmm7		, 0x170(%%rax)*/\n\t"\
		"/*movaps		%%xmm6		, 0x160(%%rax)*/\n\t"\
"movaps	0x170(%%rbx),%%xmm1		\n\t"\
"movaps	0x160(%%rbx),%%xmm0		\n\t"\
"movaps	%%xmm7,0x010(%%rdx)		\n\t"\
"movaps	%%xmm6,     (%%rdx)		\n\t"/* r28,29 */\
"movaps	%%xmm1,0x110(%%rdx)		\n\t"\
"movaps	%%xmm0,0x100(%%rdx)		\n\t"/* r38,39 */\
"\n\t"\
		"movq		%[__c0E]		,%%rcx\n\t"\
		"subpd		%%xmm5		,%%xmm2\n\t"\
		"subpd		%%xmm4		,%%xmm3\n\t"\
		"addpd		%%xmm5		,%%xmm5\n\t"\
		"addpd		%%xmm4		,%%xmm4\n\t"\
		"addpd		%%xmm2		,%%xmm5\n\t"\
		"addpd		%%xmm3		,%%xmm4\n\t"\
		"movaps		%%xmm5		,%%xmm0\n\t"\
		"movaps		%%xmm3		,%%xmm1\n\t"\
		"mulpd		      (%%rcx)	,%%xmm5\n\t"\
		"mulpd		      (%%rcx)	,%%xmm3\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm0\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm1\n\t"\
		"subpd		%%xmm0		,%%xmm3\n\t"\
		"addpd		%%xmm1		,%%xmm5\n\t"\
		"/*movaps		%%xmm3		,%%xmm1*/\n\t"\
		"/*movaps		%%xmm5		,%%xmm0*/\n\t"\
		"/*unpckhpd	 0x0f0(%%rbx)	,%%xmm1*/\n\t"\
		"/*unpcklpd	 0x0f0(%%rbx)	,%%xmm3*/\n\t"\
		"/*movaps		%%xmm1		, 0x0f0(%%rbx)*/\n\t"\
		"/*unpckhpd	 0x0e0(%%rbx)	,%%xmm0*/\n\t"\
		"/*unpcklpd	 0x0e0(%%rbx)	,%%xmm5*/\n\t"\
		"/*movaps		%%xmm0		, 0x0e0(%%rbx)*/\n\t"\
		"/*movaps		%%xmm3		, 0x0f0(%%rax)*/\n\t"\
		"/*movaps		%%xmm5		, 0x0e0(%%rax)*/\n\t"\
"movaps	0x0f0(%%rbx),%%xmm1		\n\t"\
"movaps	0x0e0(%%rbx),%%xmm0		\n\t"\
"movaps	%%xmm3,0x050(%%rdx)		\n\t"\
"movaps	%%xmm5,0x040(%%rdx)		\n\t"/* r2c,2d */\
"movaps	%%xmm1,0x150(%%rdx)		\n\t"\
"movaps	%%xmm0,0x140(%%rdx)		\n\t"/* r3c,3d */\
"\n\t"\
		"movq		%[__c1E]		,%%rcx\n\t"\
		"movaps		%%xmm2		,%%xmm0\n\t"\
		"movaps		%%xmm4		,%%xmm1\n\t"\
		"mulpd		      (%%rcx)	,%%xmm2\n\t"\
		"mulpd		      (%%rcx)	,%%xmm4\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm0\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm1\n\t"\
		"subpd		%%xmm0		,%%xmm4\n\t"\
		"addpd		%%xmm1		,%%xmm2\n\t"\
		"/*movaps		%%xmm4		,%%xmm1*/\n\t"\
		"/*movaps		%%xmm2		,%%xmm0*/\n\t"\
		"/*unpckhpd	 0x1f0(%%rbx)	,%%xmm1*/\n\t"\
		"/*unpcklpd	 0x1f0(%%rbx)	,%%xmm4*/\n\t"\
		"/*movaps		%%xmm1		, 0x1f0(%%rbx)*/\n\t"\
		"/*unpckhpd	 0x1e0(%%rbx)	,%%xmm0*/\n\t"\
		"/*unpcklpd	 0x1e0(%%rbx)	,%%xmm2*/\n\t"\
		"/*movaps		%%xmm0		, 0x1e0(%%rbx)*/\n\t"\
		"/*movaps		%%xmm4		, 0x1f0(%%rax)*/\n\t"\
		"/*movaps		%%xmm2		, 0x1e0(%%rax)*/\n\t"\
"movaps	0x1f0(%%rbx),%%xmm1		\n\t"\
"movaps	0x1e0(%%rbx),%%xmm0		\n\t"\
"movaps	%%xmm4,0x070(%%rdx)		\n\t"\
"movaps	%%xmm2,0x060(%%rdx)		\n\t"/* r2e,2f */\
"movaps	%%xmm1,0x170(%%rdx)		\n\t"\
"movaps	%%xmm0,0x160(%%rdx)		\n\t"/* r3e,3f */\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__add1] "m" (Xadd1)\
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__r00] "m" (Xr00)\
		 ,[__r08] "m" (Xr08)\
		 ,[__r10] "m" (Xr10)\
		 ,[__r20] "m" (Xr20)\
		 ,[__r28] "m" (Xr28)\
		 ,[__r30] "m" (Xr30)\
		 ,[__c01] "m" (Xc01)\
		 ,[__c02] "m" (Xc02)\
		 ,[__c04] "m" (Xc04)\
		 ,[__c06] "m" (Xc06)\
		 ,[__c08] "m" (Xc08)\
		 ,[__c0A] "m" (Xc0A)\
		 ,[__c0C] "m" (Xc0C)\
		 ,[__c0E] "m" (Xc0E)\
		 ,[__c10] "m" (Xc10)\
		 ,[__c12] "m" (Xc12)\
		 ,[__c14] "m" (Xc14)\
		 ,[__c16] "m" (Xc16)\
		 ,[__c18] "m" (Xc18)\
		 ,[__c1A] "m" (Xc1A)\
		 ,[__c1C] "m" (Xc1C)\
		 ,[__c1E] "m" (Xc1E)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"		/* Clobbered registers */\
	);\
	}

#endif	// SSE2 or AVX?

#endif	/* radix32_wrapper_square_gcc_h_included */

