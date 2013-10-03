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
#ifndef radix16_wrapper_square_gcc_h_included
#define radix16_wrapper_square_gcc_h_included

#ifdef USE_AVX	// See the corresponding radix16_dyadic_square header file for notes re. the AVX-mode data layout

	#define SSE2_RADIX16_WRAPPER_DIF(Xadd0,Xadd1,Xadd2,Xadd3,Xr1,Xr9,Xr17,Xr25,Xisrt2,Xcc0,Xc1,Xc2,Xc3,Xc4,Xc5,Xc6,Xc7,Xc8,Xc9,Xc10,Xc11,Xc12,Xc13,Xc14,Xc15)\
	{\
	__asm__ volatile (\
	"/*************************************************************/\n\t"\
	"/*                  1st set of inputs:                       */\n\t"\
	"/*************************************************************/\n\t"\
		"movq	%[__r1] ,%%rsi	\n\t"\
		"movq	%[__add0],%%rax\n\t"\
		"movq	%[__add1],%%rbx\n\t"\
		"movq	%[__add2],%%rcx\n\t"\
		"movq	%[__add3],%%rdx\n\t"\
	"/**** Start with 4-way interleaving: ****/\n\t"\
	"/* a[j+p0]: Inputs from [add0,add0+0x100,add1,add1+0x100]+0x0. Outputs into r1 +0/1, 8/9, 16/17, 24/25: */	\n\t"\
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
		"vmovaps 	%%ymm5 ,0x200(%%rsi)					\n\t		vmovaps %%ymm13,0x220(%%rsi)				/* outB	*/	\n\t"\
		"vmovaps 	%%ymm7 ,0x300(%%rsi)					\n\t		vmovaps %%ymm15,0x320(%%rsi)				/* outD	*/	\n\t"\
		"vmovaps 	%%ymm0 ,     (%%rsi)					\n\t		vmovaps %%ymm2 ,0x020(%%rsi)				/* outA	*/	\n\t"\
		"vmovaps 	%%ymm1 ,0x100(%%rsi)					\n\t		vmovaps %%ymm3 ,0x120(%%rsi)				/* outC	*/	\n\t"\
	"/* a[j+p2]: Inputs from [add0,add0+0x100,add1,add1+0x100]+0x40. Outputs into r3 +0/1, 8/9, 16/17, 24/25: */		\n\t"\
		"addq	$0x40,%%rax	\n\t"\
		"addq	$0x40,%%rbx	\n\t"\
		"addq	$0x40,%%rcx	\n\t"\
		"addq	$0x40,%%rdx	\n\t"\
		"addq	$0x40,%%rsi	\n\t"\
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
		"vmovaps 	%%ymm5 ,0x200(%%rsi)					\n\t		vmovaps %%ymm13,0x220(%%rsi)				/* outB	*/	\n\t"\
		"vmovaps 	%%ymm7 ,0x300(%%rsi)					\n\t		vmovaps %%ymm15,0x320(%%rsi)				/* outD	*/	\n\t"\
		"vmovaps 	%%ymm0 ,     (%%rsi)					\n\t		vmovaps %%ymm2 ,0x020(%%rsi)				/* outA	*/	\n\t"\
		"vmovaps 	%%ymm1 ,0x100(%%rsi)					\n\t		vmovaps %%ymm3 ,0x120(%%rsi)				/* outC	*/	\n\t"\
	"/* a[j+p4]: Inputs from [add0,add0+0x100,add1,add1+0x100]+0x80. Outputs into r5 +0/1, 8/9, 16/17, 24/25: */		\n\t"\
		"addq	$0x40,%%rax	\n\t"\
		"addq	$0x40,%%rbx	\n\t"\
		"addq	$0x40,%%rcx	\n\t"\
		"addq	$0x40,%%rdx	\n\t"\
		"addq	$0x40,%%rsi	\n\t"\
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
		"vmovaps 	%%ymm5 ,0x200(%%rsi)					\n\t		vmovaps %%ymm13,0x220(%%rsi)				/* outB	*/	\n\t"\
		"vmovaps 	%%ymm7 ,0x300(%%rsi)					\n\t		vmovaps %%ymm15,0x320(%%rsi)				/* outD	*/	\n\t"\
		"vmovaps 	%%ymm0 ,     (%%rsi)					\n\t		vmovaps %%ymm2 ,0x020(%%rsi)				/* outA	*/	\n\t"\
		"vmovaps 	%%ymm1 ,0x100(%%rsi)					\n\t		vmovaps %%ymm3 ,0x120(%%rsi)				/* outC	*/	\n\t"\
	"/* a[j+p6]: Inputs from [add0,add0+0x100,add1,add1+0x100]+0xc0. Outputs into r7 +0/1, 8/9, 16/17, 24/25: */		\n\t"\
		"addq	$0x40,%%rax	\n\t"\
		"addq	$0x40,%%rbx	\n\t"\
		"addq	$0x40,%%rcx	\n\t"\
		"addq	$0x40,%%rdx	\n\t"\
		"addq	$0x40,%%rsi	\n\t"\
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
		"vmovaps 	%%ymm5 ,0x200(%%rsi)					\n\t		vmovaps %%ymm13,0x220(%%rsi)				/* outB	*/	\n\t"\
		"vmovaps 	%%ymm7 ,0x300(%%rsi)					\n\t		vmovaps %%ymm15,0x320(%%rsi)				/* outD	*/	\n\t"\
		"vmovaps 	%%ymm0 ,     (%%rsi)					\n\t		vmovaps %%ymm2 ,0x020(%%rsi)				/* outA	*/	\n\t"\
		"vmovaps 	%%ymm1 ,0x100(%%rsi)					\n\t		vmovaps %%ymm3 ,0x120(%%rsi)				/* outC	*/	\n\t"\
	"/*****************/\n\t"\
	"/* Radix-16 DIF: */\n\t"\
	"/*****************/\n\t"\
		"movq	%[__r1] ,%%rcx	\n\t"\
	"/*...Block 1: */														\n\t	/*...Block 2: */\n\t"\
		"leaq	0x560(%%rcx),%%rdi	/* c2 */\n\t"\
		"leaq	0x4e0(%%rcx),%%rdx	/* c4 */\n\t"\
		"vmovaps	0x040(%%rcx),%%ymm0	/* ymm0 <-     a[jt+p4] */			\n\t		vmovaps		0x100(%%rcx),%%ymm8	/* ymm10 <-     a[jt+p2] */			\n\t"\
		"vmovaps	0x060(%%rcx),%%ymm1	/* ymm1 <-     a[jp+p4] */			\n\t		vmovaps		0x120(%%rcx),%%ymm9	/* ymm11 <-     a[jp+p2] */			\n\t"\
		"vmovaps	%%ymm0		,%%ymm2	/* ymm2 <- cpy a[jt+p4] */			\n\t		vmovaps		%%ymm8 	,%%ymm10	/* ymm10 <- cpy a[jt+p2] */			\n\t"\
		"vmovaps	%%ymm1		,%%ymm3	/* ymm3 <- cpy a[jp+p4] */			\n\t		vmovaps		%%ymm9 	,%%ymm11	/* ymm11 <- cpy a[jp+p2] */			\n\t"\
	"/***************************************************************************/	\n\t"\
	"/*** From hereon, things are identical to the code in radix16_dif_pass: ****/	\n\t"\
	"/***************************************************************************/	\n\t"\
		"vmulpd		(%%rdx)		,%%ymm0,%%ymm0		/* a[jt+p4]*c4 */		\n\t		vmulpd		     (%%rdi)	,%%ymm8 ,%%ymm8 		/* a[jt+p2]*c2 */	\n\t"\
		"vmulpd		(%%rdx)		,%%ymm1,%%ymm1		/* a[jp+p4]*c4 */		\n\t		vmulpd		     (%%rdi)	,%%ymm9 ,%%ymm9 		/* a[jp+p2]*c2 */	\n\t"\
		"vmulpd		0x020(%%rdx),%%ymm2,%%ymm2		/* a[jt+p4]*s4 */		\n\t		vmulpd		0x020(%%rdi),%%ymm10,%%ymm10		/* a[jt+p2]*s2 */	\n\t"\
		"vmulpd		0x020(%%rdx),%%ymm3,%%ymm3		/* a[jp+p4]*s4 */		\n\t		vmulpd		0x020(%%rdi),%%ymm11,%%ymm11		/* a[jp+p2]*s2 */	\n\t"\
		"vaddpd		%%ymm2		,%%ymm1,%%ymm1	/* ymm1 <- t6 */			\n\t		vaddpd		%%ymm10		,%%ymm9 ,%%ymm9 	/* ymm9  <- t10*/		\n\t"\
		"vsubpd		%%ymm3		,%%ymm0,%%ymm0	/* ymm0 <- t5 */			\n\t		vsubpd		%%ymm11		,%%ymm8 ,%%ymm8 	/* ymm8  <- t9 */		\n\t"\
		"vmovaps	%%ymm1		,%%ymm3		/* ymm3 <- cpy t6 */			\n\t		vmovaps		%%ymm9 		,%%ymm11		/* ymm11 <- cpy t10*/		\n\t"\
		"vmovaps	%%ymm0		,%%ymm2		/* ymm2 <- cpy t5 */			\n\t		vmovaps		%%ymm8 		,%%ymm10		/* ymm10 <- cpy t9 */		\n\t"\
		"addq		$0x040,%%rdx	/* c12 */								\n\t		addq	$0x040,%%rdi	/* c10 */		\n\t"\
		"vmovaps	0x0c0(%%rcx),%%ymm4			/* ymm4 <-     a[jt+p12] */	\n\t		vmovaps		0x180(%%rcx),%%ymm12			/* ymm12 <-     a[jt+p10] */\n\t"\
		"vmovaps	0x0e0(%%rcx),%%ymm5			/* ymm5 <-     a[jp+p12] */	\n\t		vmovaps		0x1a0(%%rcx),%%ymm13			/* ymm13 <-     a[jp+p10] */\n\t"\
		"vmovaps	%%ymm4		,%%ymm6			/* ymm6 <- cpy a[jt+p12] */	\n\t		vmovaps			%%ymm12	,%%ymm14			/* ymm14 <- cpy a[jt+p10] */\n\t"\
		"vmovaps	%%ymm5		,%%ymm7			/* ymm7 <- cpy a[jp+p12] */	\n\t		vmovaps			%%ymm13	,%%ymm15			/* ymm15 <- cpy a[jp+p10] */\n\t"\
		"vmulpd		(%%rdx)		,%%ymm4,%%ymm4		/* a[jt+p12]*c12 */		\n\t		vmulpd		     (%%rdi)	,%%ymm12,%%ymm12		/* a[jt+p10]*c10 */		\n\t"\
		"vmulpd		(%%rdx)		,%%ymm5,%%ymm5		/* a[jp+p12]*c12 */		\n\t		vmulpd		     (%%rdi)	,%%ymm13,%%ymm13		/* a[jp+p10]*c10 */		\n\t"\
		"vmulpd		0x020(%%rdx),%%ymm6,%%ymm6		/* a[jt+p12]*s12 */		\n\t		vmulpd		0x020(%%rdi),%%ymm14,%%ymm14		/* a[jt+p10]*s10 */		\n\t"\
		"vmulpd		0x020(%%rdx),%%ymm7,%%ymm7		/* a[jp+p12]*s12 */		\n\t		vmulpd		0x020(%%rdi),%%ymm15,%%ymm15		/* a[jp+p10]*s10 */		\n\t"\
		"vaddpd		%%ymm6		,%%ymm5,%%ymm5	/* ymm5 <- it */			\n\t		vaddpd		%%ymm14		,%%ymm13,%%ymm13	/* ymm13 <- it */						\n\t"\
		"vsubpd		%%ymm7		,%%ymm4,%%ymm4	/* ymm4 <- rt; ymm6,7 free*/\n\t		vsubpd		%%ymm15		,%%ymm12,%%ymm12	/* ymm12 <- rt		ymm14,7 free */		\n\t"\
		"vaddpd		%%ymm4		,%%ymm0,%%ymm0		/* ~t5 <- t5 +rt */		\n\t		vaddpd		%%ymm12		,%%ymm8 ,%%ymm8  	/* ~t13<- t13+rt */						\n\t"\
		"vaddpd		%%ymm5		,%%ymm1,%%ymm1		/* ~t6 <- t6 +it */		\n\t		vaddpd		%%ymm13		,%%ymm9 ,%%ymm9  	/* ~t14<- t14+it */						\n\t"\
		"vsubpd		%%ymm4		,%%ymm2,%%ymm2		/* ~t7 <- t5 -rt */		\n\t		vsubpd		%%ymm12		,%%ymm10,%%ymm10 	/* ~t15<- t13-rt */						\n\t"\
		"vsubpd		%%ymm5		,%%ymm3,%%ymm3		/* ~t8 <- t6 -it */		\n\t		vsubpd		%%ymm13		,%%ymm11,%%ymm11 	/* ~t16<- t14-it	ymm12,13 free */	\n\t"\
		"\n\t"\
		"/* Now do the p0,8 combo: */\n\t									/* Do the p6,14 combo - do p14 first so registers come out in same order as for p2,10 */\n\t"\
		"leaq	0x4a0(%%rcx),%%rdx	/* c8 */								\n\t		leaq	0x620(%%rcx),%%rdi	/* c14 */\n\t"\
		"vmovaps	0x080(%%rcx)	,%%ymm4		/* a[jt+p8 ] */				\n\t		vmovaps		0x1c0(%%rcx),%%ymm12		/* a[jt+p14] */				\n\t"\
		"vmovaps	0x0a0(%%rcx)	,%%ymm5		/* a[jp+p8 ] */				\n\t		vmovaps		0x1e0(%%rcx),%%ymm13		/* a[jp+p14] */				\n\t"\
		"																				vmovaps			%%ymm12	,%%ymm14		/* ymm14 <- cpy a[jt+p14] */\n\t"\
		"vmovaps	%%ymm4		,%%ymm6	/* ymm6 <- cpy a[jt+p8] */			\n\t		vmovaps			%%ymm13	,%%ymm15		/* ymm15 <- cpy a[jp+p14] */\n\t"\
		"vmovaps	%%ymm5		,%%ymm7	/* ymm7 <- cpy a[jp+p8] */			\n\t		vmulpd		     (%%rdi)	,%%ymm12,%%ymm12		/* a[jt+p14]*c14 */			\n\t"\
		"																				vmulpd		     (%%rdi)	,%%ymm13,%%ymm13		/* a[jp+p14]*c14 */			\n\t"\
		"vmulpd		(%%rdx)		,%%ymm4,%%ymm4		/* a[jt+p8]*c8 */		\n\t		vmulpd		0x020(%%rdi),%%ymm14,%%ymm14		/* a[jt+p14]*s14 */			\n\t"\
		"vmulpd		(%%rdx)		,%%ymm5,%%ymm5		/* a[jp+p8]*c8 */		\n\t		vmulpd		0x020(%%rdi),%%ymm15,%%ymm15		/* a[jp+p14]*s14 */			\n\t"\
		"vmulpd		0x020(%%rdx),%%ymm6,%%ymm6		/* a[jt+p8]*s8 */		\n\t		vaddpd		%%ymm14		,%%ymm13,%%ymm13		/* ymm13 <- it */			\n\t"\
		"vmulpd		0x020(%%rdx),%%ymm7,%%ymm7		/* a[jp+p8]*s8 */		\n\t		vsubpd		%%ymm15		,%%ymm12,%%ymm12		/* ymm12 <- rt ymm14,15 free */	\n\t"\
		"vaddpd		%%ymm6		,%%ymm5,%%ymm5	/* ymm5 <- it */			\n\t		vmovaps		%%ymm13		,0x1e0(%%rcx)	/* Store it in t16*/		\n\t"\
		"vsubpd		%%ymm7		,%%ymm4,%%ymm4	/* ymm4 <- rt; ymm6,7 free*/\n\t		vmovaps		%%ymm12		,0x1c0(%%rcx)	/* Store rt in t15*/		\n\t"\
		"																				subq	$0x040,%%rdi	/* c6  */	\n\t"\
	"/* Real parts:*/\n\t															/* Real parts: */\n\t"\
		"vmovaps		 (%%rcx),%%ymm6		/* a[jt    ] */					\n\t		vmovaps		0x140(%%rcx),%%ymm12		/* a[jt+p6 ] */				\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm7		/* a[jp    ] */					\n\t		vmovaps		0x160(%%rcx),%%ymm13		/* a[jp+p6 ] */				\n\t"\
		"																				vmovaps			%%ymm12	,%%ymm14		/* ymm14 <- cpy a[jt+p6] */		\n\t"\
		"																				vmovaps			%%ymm13	,%%ymm15		/* ymm15 <- cpy a[jp+p6] */		\n\t"\
		"vsubpd		%%ymm4	,%%ymm6,%%ymm6	/* ~t3 <- t1 -rt */				\n\t		vmulpd		     (%%rdi)	,%%ymm12,%%ymm12		/* a[jt+p6]*c6 */			\n\t"\
		"vsubpd		%%ymm5	,%%ymm7,%%ymm7	/* ~t4 <- t2 -it */				\n\t		vmulpd		     (%%rdi)	,%%ymm13,%%ymm13		/* a[jp+p6]*c6 */			\n\t"\
		"vaddpd		%%ymm4	,%%ymm4,%%ymm4	/*          2*rt */				\n\t		vmulpd		0x020(%%rdi),%%ymm14,%%ymm14		/* a[jt+p6]*s6 */			\n\t"\
		"vaddpd		%%ymm5	,%%ymm5,%%ymm5	/*          2*it */				\n\t		vmulpd		0x020(%%rdi),%%ymm15,%%ymm15		/* a[jp+p6]*s6 */			\n\t"\
		"vaddpd		%%ymm6	,%%ymm4,%%ymm4	/* ~t1 <- t1 +rt */				\n\t		vaddpd		%%ymm14		,%%ymm13,%%ymm13	/* ymm13 <- t14*/				\n\t"\
		"vaddpd		%%ymm7	,%%ymm5,%%ymm5	/* ~t2 <- t2 +it; ymm4,5 free */\n\t		vsubpd		%%ymm15		,%%ymm12,%%ymm12	/* ymm12 <- t13*/				\n\t"\
		"																				vmovaps		%%ymm13		,%%ymm15		/* ymm15 <- cpy t14*/			\n\t"\
		"																				vmovaps		%%ymm12		,%%ymm14		/* ymm14 <- cpy t13*/			\n\t"\
		"																				vsubpd		0x1c0(%%rcx),%%ymm12,%%ymm12		/* ~t15<- t13-rt */			\n\t"\
		"																				vsubpd		0x1e0(%%rcx),%%ymm13,%%ymm13		/* ~t16<- t14-it */			\n\t"\
		"																				vaddpd		0x1c0(%%rcx),%%ymm14,%%ymm14		/* ~t13<- t13+rt */			\n\t"\
		"																				vaddpd		0x1e0(%%rcx),%%ymm15,%%ymm15		/* ~t14<- t14+it */			\n\t"\
		"\n\t"\
	"/* Finish radix-4 butterfly and store results into temporary-array slots: */\n\t	/* Finish radix-4 butterfly and store results into temporary-array slots: */\n\t"\
		"vsubpd		%%ymm0		,%%ymm4,%%ymm4	/*~t5 =t1 -t5 */			\n\t		vsubpd		%%ymm14		,%%ymm8,%%ymm8 	/*~t13*/						\n\t"\
		"vsubpd		%%ymm1		,%%ymm5,%%ymm5	/*~t6 =t2 -t6 */			\n\t		vsubpd		%%ymm15		,%%ymm9,%%ymm9 	/*~t14*/						\n\t"\
		"vmovaps	%%ymm4		,0x080(%%rcx)	/* a[jt+p8 ] <- ~t5 */		\n\t		vmovaps		%%ymm8 		,0x180(%%rcx)	/* a[jt+p8 ] <- ~t13*/		\n\t"\
		"vmovaps	%%ymm5		,0x0a0(%%rcx)	/* a[jp+p8 ] <- ~t6 */		\n\t		vmovaps		%%ymm9 		,0x1a0(%%rcx)	/* a[jp+p8 ] <- ~t14*/		\n\t"\
		"vaddpd		%%ymm0		,%%ymm0,%%ymm0		/* 2*t5 */				\n\t		vaddpd		%%ymm14		,%%ymm14,%%ymm14	/* 2*t13*/						\n\t"\
		"vaddpd		%%ymm1		,%%ymm1,%%ymm1		/* 2*t6 */				\n\t		vaddpd		%%ymm15		,%%ymm15,%%ymm15	/* 2*t14*/						\n\t"\
		"vaddpd		%%ymm4		,%%ymm0,%%ymm0		/*~t1 =t1 +t5 */		\n\t		vaddpd		%%ymm8 		,%%ymm14,%%ymm14	/*~t9 */						\n\t"\
		"vaddpd		%%ymm5		,%%ymm1,%%ymm1		/*~t2 =t2 +t6 */		\n\t		vaddpd		%%ymm9 		,%%ymm15,%%ymm15	/*~t10*/						\n\t"\
		"vmovaps	%%ymm0		,     (%%rcx)	/* a[jt    ] <- ~t1 */		\n\t		vmovaps		%%ymm14		,0x100(%%rcx)	/* a[jt    ] <- ~t9 */		\n\t"\
		"vmovaps	%%ymm1		,0x020(%%rcx)	/* a[jp    ] <- ~t2 */		\n\t		vmovaps		%%ymm15		,0x120(%%rcx)	/* a[jp    ] <- ~t10*/		\n\t"\
		"\n\t"\
		"vsubpd		%%ymm3		,%%ymm6,%%ymm6	/*~t3 =t3 -t8 */			\n\t		vsubpd		%%ymm13		,%%ymm10,%%ymm10	/*~t11*/						\n\t"\
		"vsubpd		%%ymm2		,%%ymm7,%%ymm7	/*~t8 =t4 -t7 */			\n\t		vsubpd		%%ymm12		,%%ymm11,%%ymm11	/*~t16*/						\n\t"\
		"vmovaps	%%ymm6		,0x040(%%rcx)	/* a[jt+p4 ] <- ~t3 */		\n\t		vmovaps		%%ymm10		,0x140(%%rcx)	/* a[jt+p4 ] <- ~t11*/		\n\t"\
		"vmovaps	%%ymm7		,0x0e0(%%rcx)	/* a[jp+p12] <- ~t8 */		\n\t		vmovaps		%%ymm11		,0x1e0(%%rcx)	/* a[jp+p12] <- ~t16*/		\n\t"\
		"vaddpd		%%ymm3		,%%ymm3,%%ymm3		/* 2*t8 */				\n\t		vaddpd		%%ymm13		,%%ymm13,%%ymm13	/* 2*t16*/						\n\t"\
		"vaddpd		%%ymm2		,%%ymm2,%%ymm2		/* 2*t7 */				\n\t		vaddpd		%%ymm12		,%%ymm12,%%ymm12	/* 2*t15*/						\n\t"\
		"vaddpd		%%ymm6		,%%ymm3,%%ymm3		/*~t7 =t3 +t8 */		\n\t		vaddpd		%%ymm10		,%%ymm13,%%ymm13	/*~t15*/						\n\t"\
		"vaddpd		%%ymm7		,%%ymm2,%%ymm2		/*~t4 =t4 +t7 */		\n\t		vaddpd		%%ymm11		,%%ymm12,%%ymm12	/*~t12*/						\n\t"\
		"vmovaps	%%ymm3		,0x0c0(%%rcx)	/* a[jt+p12] <- ~t7 */		\n\t		vmovaps		%%ymm13		,0x1c0(%%rcx)	/* a[jt+p12] <- ~t15*/		\n\t"\
		"vmovaps	%%ymm2		,0x060(%%rcx)	/* a[jp+p4 ] <- ~t4 */		\n\t		vmovaps		%%ymm12		,0x160(%%rcx)	/* a[jp+p4 ] <- ~t12*/		\n\t"\
		"\n\t"\
	"/*...Block 3: */														\n\t	/*...Block 4: */\n\t"\
	"/*	SSE2_RADIX4_DIF_4TWIDDLE(r17,r21,r19,r23, r17,c1): */				\n\t	/* SSE2_RADIX4_DIF_4TWIDDLE(r25,r29,r27,r31, r25,c3): */\n\t"\
		"/* Do the p0,p8 combo: */											\n\t	/* Do the p0,p8 combo: */		\n\t"\
		"leaq	0x660(%%rcx),%%rbx	/* c1 */								\n\t/* All __r and __c pointers incr by +0x100 in rcol w.r.to lcol: */\n\t"\
		"addq	$0x200,%%rcx		/* r17 */								\n\t		/* r25 */\n\t"\
		"\n\t																\n\t		/* c3  */\n\t"\
		"vmovaps		 (%%rcx),%%ymm0		/* a[jt   ] */					\n\t		vmovaps		0x100(%%rcx),%%ymm8 		/* a[jt    ] */\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm1		/* a[jp   ] */					\n\t		vmovaps		0x120(%%rcx),%%ymm9 		/* a[jp    ] */\n\t"\
		"vmovaps	0x080(%%rcx),%%ymm4		/* a[jt+p8 ] */					\n\t		vmovaps		0x180(%%rcx),%%ymm12		/* a[jt+p8 ] */\n\t"\
		"vmovaps	0x0a0(%%rcx),%%ymm5		/* a[jp+p8 ] */					\n\t		vmovaps		0x1a0(%%rcx),%%ymm13		/* a[jp+p8 ] */\n\t"\
		"vmovaps		 (%%rbx),%%ymm6		/* c0 */						\n\t		vmovaps		0x100(%%rbx),%%ymm14		/* c0 */\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm7		/* s0 */						\n\t		vmovaps		0x120(%%rbx),%%ymm15		/* s0 */\n\t"\
		"vmovaps	%%ymm0		,%%ymm2		/* ymm2 <- cpy a[jt   ] */		\n\t		vmovaps		%%ymm8 		,%%ymm10		/* ymm10 <- cpy a[jt   ] */\n\t"\
		"vmovaps	%%ymm1		,%%ymm3		/* ymm3 <- cpy a[jp   ] */		\n\t		vmovaps		%%ymm9 		,%%ymm11		/* ymm11 <- cpy a[jp   ] */\n\t"\
		"\n\t																\n\t"\
		"vmulpd   	%%ymm6		,%%ymm0,%%ymm0		/* a[jt   ]*c0 */		\n\t		vmulpd   	%%ymm14		,%%ymm8 ,%%ymm8 		/* a[jt   ]*c0 */\n\t"\
		"vmulpd   	%%ymm6		,%%ymm1,%%ymm1		/* a[jp   ]*c0 */		\n\t		vmulpd   	%%ymm14		,%%ymm9 ,%%ymm9 		/* a[jp   ]*c0 */\n\t"\
		"vmulpd   	%%ymm7		,%%ymm2,%%ymm2		/* a[jt   ]*s0 */		\n\t		vmulpd   	%%ymm15		,%%ymm10,%%ymm10		/* a[jt   ]*s0 */\n\t"\
		"vmulpd   	%%ymm7		,%%ymm3,%%ymm3	/* a[jp   ]*s0 ymm6,7 free*/\n\t		vmulpd   	%%ymm15		,%%ymm11,%%ymm11		/* a[jp   ]*s0	ymm14,7 free */\n\t"\
		"vmovaps	%%ymm4		,%%ymm6		/* ymm6 <- cpy a[jt+p8 ] */		\n\t		vmovaps		%%ymm12		,%%ymm14		/* ymm14 <- cpy a[jt+p8 ] */\n\t"\
		"vmovaps	%%ymm5		,%%ymm7		/* ymm7 <- cpy a[jp+p8 ] */		\n\t		vmovaps		%%ymm13		,%%ymm15		/* ymm15 <- cpy a[jp+p8 ] */\n\t"\
		"vaddpd   	%%ymm2		,%%ymm1,%%ymm1		/* ymm1 <- t2 */		\n\t		vaddpd   	%%ymm10		,%%ymm9 ,%%ymm9 		/* ymm9  <- t2 */			\n\t"\
		"vsubpd   	%%ymm3		,%%ymm0,%%ymm0		/* ymm0 <- t1 */		\n\t		vsubpd   	%%ymm11		,%%ymm8 ,%%ymm8 		/* ymm8  <- t1 */			\n\t"\
		"vmulpd		0x040(%%rbx),%%ymm4,%%ymm4		/* a[jt+p8 ]*c8 */		\n\t		vmulpd		0x140(%%rbx),%%ymm12,%%ymm12		/* a[jt+p8 ]*c8 */\n\t"\
		"vmulpd		0x040(%%rbx),%%ymm5,%%ymm5		/* a[jp+p8 ]*c8 */		\n\t		vmulpd		0x140(%%rbx),%%ymm13,%%ymm13		/* a[jp+p8 ]*c8 */\n\t"\
		"vmulpd		0x060(%%rbx),%%ymm6,%%ymm6		/* a[jt+p8 ]*s8 */		\n\t		vmulpd		0x160(%%rbx),%%ymm14,%%ymm14		/* a[jt+p8 ]*s8 */\n\t"\
		"vmulpd		0x060(%%rbx),%%ymm7,%%ymm7		/* a[jp+p8 ]*s8 */		\n\t		vmulpd		0x160(%%rbx),%%ymm15,%%ymm15		/* a[jp+p8 ]*s8 */\n\t"\
		"vmovaps	%%ymm0		,%%ymm2		/* ymm2 <- cpy t1 */			\n\t		vmovaps		%%ymm8 		,%%ymm10		/* ymm10 <- cpy t1 */		\n\t"\
		"vmovaps	%%ymm1		,%%ymm3		/* ymm3 <- cpy t2 */			\n\t		vmovaps		%%ymm9 		,%%ymm11		/* ymm11 <- cpy t2 */		\n\t"\
		"vaddpd		%%ymm6	    ,%%ymm5,%%ymm5		/* ymm5 <- it */		\n\t		vaddpd		%%ymm14	    ,%%ymm13,%%ymm13		/* ymm13 <- it */\n\t"\
		"vsubpd		%%ymm7		,%%ymm4,%%ymm4	/* ymm4 <- rt ymm6,7 free*/	\n\t		vsubpd		%%ymm15		,%%ymm12,%%ymm12		/* ymm12 <- rt    ymm14,7 free */\n\t"\
		"addq	$0x0c0 ,%%rbx	/* c13 */									\n\t		/* c15 */\n\t"\
		"vmovaps	0x0c0(%%rcx),%%ymm6		/* a[jt+p12] */					\n\t		vmovaps		0x1c0(%%rcx),%%ymm14		/* a[jt+p12] */\n\t"\
		"vmovaps	0x0e0(%%rcx),%%ymm7		/* a[jp+p12] */					\n\t		vmovaps		0x1e0(%%rcx),%%ymm15		/* a[jp+p12] */\n\t"\
		"vaddpd		%%ymm4	,%%ymm0,%%ymm0		/* ~t1 <- t1 +rt */			\n\t		vaddpd		%%ymm12	,%%ymm8 ,%%ymm8 		/* ~t1 <- t1 +rt */\n\t"\
		"vaddpd		%%ymm5	,%%ymm1,%%ymm1		/* ~t2 <- t2 +it */			\n\t		vaddpd		%%ymm13	,%%ymm9 ,%%ymm9 		/* ~t2 <- t2 +it */\n\t"\
		"vsubpd		%%ymm4	,%%ymm2,%%ymm2		/* ~t3 <- t1 -rt */			\n\t		vsubpd		%%ymm12	,%%ymm10,%%ymm10		/* ~t3 <- t1 -rt */\n\t"\
		"vsubpd		%%ymm5	,%%ymm3,%%ymm3	/* ~t4 <- t2 -it ymm4,5 free*/	\n\t		vsubpd		%%ymm13	,%%ymm11,%%ymm11		/* ~t4 <- t2 -it	ymm12,5 free */\n\t"\
		"\n\t"\
	"/* Do the p4,12 combo: */												\n\t	/* Do the p4,12 combo: */\n\t"\
		"vmovaps	%%ymm6		,%%ymm4		/* ymm4 <- cpy a[jt+p12] */		\n\t		vmovaps		%%ymm14		,%%ymm12		/* ymm12 <- cpy a[jt+p12] */\n\t"\
		"vmovaps	%%ymm7		,%%ymm5		/* ymm5 <- cpy a[jp+p12] */		\n\t		vmovaps		%%ymm15		,%%ymm13		/* ymm13 <- cpy a[jp+p12] */\n\t"\
		"vmulpd		(%%rbx)		,%%ymm4,%%ymm4		/* a[jt+p12]*c12 */		\n\t		vmulpd		0x100(%%rbx),%%ymm12,%%ymm12		/* a[jt+p12]*c12 */\n\t"\
		"vmulpd		(%%rbx)		,%%ymm5,%%ymm5		/* a[jp+p12]*c12 */		\n\t		vmulpd		0x100(%%rbx),%%ymm13,%%ymm13		/* a[jp+p12]*c12 */\n\t"\
		"vmulpd		0x020(%%rbx),%%ymm6,%%ymm6		/* a[jt+p12]*s12 */		\n\t		vmulpd		0x120(%%rbx),%%ymm14,%%ymm14		/* a[jt+p12]*s12 */\n\t"\
		"vmulpd		0x020(%%rbx),%%ymm7,%%ymm7		/* a[jp+p12]*s12 */		\n\t		vmulpd		0x120(%%rbx),%%ymm15,%%ymm15		/* a[jp+p12]*s12 */\n\t"\
		"vaddpd		%%ymm6		,%%ymm5,%%ymm5	/* ymm5 <- it */			\n\t		vaddpd		%%ymm14		,%%ymm13,%%ymm13	/* ymm13 <- it */\n\t"\
		"vsubpd		%%ymm7		,%%ymm4,%%ymm4	/* ymm4 <- rt */			\n\t		vsubpd		%%ymm15		,%%ymm12,%%ymm12	/* ymm12 <- rt */\n\t"\
		"vmovaps	%%ymm5		,0x020(%%rcx)	/* store it */				\n\t		vmovaps		%%ymm13		,0x120(%%rcx)	/* store it */\n\t"\
		"vmovaps	%%ymm4		,     (%%rcx)	/* store rt */				\n\t		vmovaps		%%ymm12		,0x100(%%rcx)	/* store rt */\n\t"\
		"\n\t"\
		"subq	$0x040 ,%%rbx	/* c5 */									\n\t		/* c7  */\n\t"\
		"vmovaps	0x040(%%rcx),%%ymm4		/* a[jt+p4] */					\n\t		vmovaps		0x140(%%rcx),%%ymm12		/* a[jt+p4] */\n\t"\
		"vmovaps	0x060(%%rcx),%%ymm5		/* a[jp+p4] */					\n\t		vmovaps		0x160(%%rcx),%%ymm13		/* a[jp+p4] */\n\t"\
		"vmovaps		%%ymm4	,%%ymm6		/* ymm4 <- cpy a[jt+p4] */		\n\t		vmovaps			%%ymm12	,%%ymm14		/* ymm12 <- cpy a[jt+p4] */\n\t"\
		"vmovaps		%%ymm5	,%%ymm7		/* ymm5 <- cpy a[jp+p4] */		\n\t		vmovaps			%%ymm13	,%%ymm15		/* ymm13 <- cpy a[jp+p4] */\n\t"\
		"vmulpd		     (%%rbx),%%ymm4,%%ymm4		/* a[jt+p4]*c4 */		\n\t		vmulpd		0x100(%%rbx),%%ymm12,%%ymm12		/* a[jt+p4]*c4 */\n\t"\
		"vmulpd		     (%%rbx),%%ymm5,%%ymm5		/* a[jp+p4]*c4 */		\n\t		vmulpd		0x100(%%rbx),%%ymm13,%%ymm13		/* a[jp+p4]*c4 */\n\t"\
		"vmulpd		0x020(%%rbx),%%ymm6,%%ymm6		/* a[jt+p4]*s4 */		\n\t		vmulpd		0x120(%%rbx),%%ymm14,%%ymm14		/* a[jt+p4]*s4 */\n\t"\
		"vmulpd		0x020(%%rbx),%%ymm7,%%ymm7		/* a[jp+p4]*s4 */		\n\t		vmulpd		0x120(%%rbx),%%ymm15,%%ymm15		/* a[jp+p4]*s4 */\n\t"\
		"vaddpd		%%ymm6		,%%ymm5,%%ymm5	/* ymm5 <- t6 */			\n\t		vaddpd		%%ymm14		,%%ymm13,%%ymm13	/* ymm13 <- t6 */\n\t"\
		"vsubpd		%%ymm7		,%%ymm4,%%ymm4	/* ymm4 <- t5 ymm6,7 free */\n\t		vsubpd		%%ymm15		,%%ymm12,%%ymm12	/* ymm12 <- t5 	ymm14,7 free */\n\t"\
		"vmovaps	%%ymm5		,%%ymm7		/* ymm7 <- cpy t6 */			\n\t		vmovaps		%%ymm13		,%%ymm15		/* ymm15 <- cpy t6 */\n\t"\
		"vmovaps	%%ymm4		,%%ymm6		/* ymm6 <- cpy t5 */			\n\t		vmovaps		%%ymm12		,%%ymm14		/* ymm14 <- cpy t5 */\n\t"\
		"vsubpd		     (%%rcx),%%ymm4,%%ymm4		/* ~t7 <- t5 -rt */		\n\t		vsubpd		0x100(%%rcx),%%ymm12,%%ymm12		/* ~t7 <- t5 -rt */\n\t"\
		"vsubpd		0x020(%%rcx),%%ymm5,%%ymm5		/* ~t8 <- t6 -it */		\n\t		vsubpd		0x120(%%rcx),%%ymm13,%%ymm13		/* ~t8 <- t6 -it */\n\t"\
		"vaddpd		     (%%rcx),%%ymm6,%%ymm6		/* ~t5 <- t5 +rt */		\n\t		vaddpd		0x100(%%rcx),%%ymm14,%%ymm14		/* ~t5 <- t5 +rt */\n\t"\
		"vaddpd		0x020(%%rcx),%%ymm7,%%ymm7		/* ~t6 <- t6 +it */		\n\t		vaddpd		0x120(%%rcx),%%ymm15,%%ymm15		/* ~t6 <- t6 +it */\n\t"\
		"\n\t"\
	"/* Finish radix-4 butterfly and store results into temp-array slots: */\n\t"\
		"vsubpd		%%ymm6,%%ymm0,%%ymm0	/*~t5 */						\n\t		vsubpd		%%ymm14,%%ymm8 ,%%ymm8 	/*~t5 */						\n\t"\
		"vsubpd		%%ymm7,%%ymm1,%%ymm1	/*~t6 */						\n\t		vsubpd		%%ymm15,%%ymm9 ,%%ymm9 	/*~t6 */						\n\t"\
		"vsubpd		%%ymm5,%%ymm2,%%ymm2	/*~t3 */						\n\t		vsubpd		%%ymm13,%%ymm10,%%ymm10	/*~t3 */\n\t"\
		"vsubpd		%%ymm4,%%ymm3,%%ymm3	/*~t8 */						\n\t		vsubpd		%%ymm12,%%ymm11,%%ymm11	/*~t8 */\n\t"\
		"vmovaps	%%ymm0,0x080(%%rcx)	/* a[jt+p8 ] <- ~t5 */				\n\t		vmovaps		%%ymm8 ,0x180(%%rcx)	/* a[jt+p8 ] <- ~t5 */	\n\t"\
		"vmovaps	%%ymm1,0x0a0(%%rcx)	/* a[jp+p8 ] <- ~t6 */				\n\t		vmovaps		%%ymm9 ,0x1a0(%%rcx)	/* a[jp+p8 ] <- ~t6 */	\n\t"\
		"vmovaps	%%ymm2,0x040(%%rcx)	/* a[jt+p4 ] <- ~t3 */				\n\t		vmovaps		%%ymm10,0x140(%%rcx)	/* a[jt+p4 ] <- ~t3 */\n\t"\
		"vmovaps	%%ymm3,0x0e0(%%rcx)	/* a[jp+p12] <- ~t8 */				\n\t		vmovaps		%%ymm11,0x1e0(%%rcx)	/* a[jp+p12] <- ~t8 */\n\t"\
		"vaddpd		%%ymm6,%%ymm6,%%ymm6	/* 2*t5 */						\n\t		vaddpd		%%ymm14,%%ymm14,%%ymm14	/* 2*t5 */						\n\t"\
		"vaddpd		%%ymm7,%%ymm7,%%ymm7	/* 2*t6 */						\n\t		vaddpd		%%ymm15,%%ymm15,%%ymm15	/* 2*t6 */						\n\t"\
		"vaddpd		%%ymm5,%%ymm5,%%ymm5	/* 2*t8 */						\n\t		vaddpd		%%ymm13,%%ymm13,%%ymm13	/* 2*t8 */\n\t"\
		"vaddpd		%%ymm4,%%ymm4,%%ymm4	/* 2*t7 */						\n\t		vaddpd		%%ymm12,%%ymm12,%%ymm12	/* 2*t7 */\n\t"\
		"vaddpd		%%ymm0,%%ymm6,%%ymm6	/*~t1 */						\n\t		vaddpd		%%ymm8 ,%%ymm14,%%ymm14	/*~t1 */						\n\t"\
		"vaddpd		%%ymm1,%%ymm7,%%ymm7	/*~t2 */						\n\t		vaddpd		%%ymm9 ,%%ymm15,%%ymm15	/*~t2 */						\n\t"\
		"vaddpd		%%ymm2,%%ymm5,%%ymm5	/*~t7 */						\n\t		vaddpd		%%ymm10,%%ymm13,%%ymm13	/*~t7 */\n\t"\
		"vaddpd		%%ymm3,%%ymm4,%%ymm4	/*~t4 */						\n\t		vaddpd		%%ymm11,%%ymm12,%%ymm12	/*~t4 */\n\t"\
		"vmovaps	%%ymm6,     (%%rcx)/* a[jt    ] <- ~t1 */				\n\t		vmovaps		%%ymm14,0x100(%%rcx)	/* a[jt    ] <- ~t1 */	\n\t"\
		"vmovaps	%%ymm7,0x020(%%rcx)	/* a[jp    ] <- ~t2 */				\n\t		vmovaps		%%ymm15,0x120(%%rcx)	/* a[jp    ] <- ~t2 */	\n\t"\
		"vmovaps	%%ymm5,0x0c0(%%rcx)	/* a[jt+p12] <- ~t7 */				\n\t		vmovaps		%%ymm13,0x1c0(%%rcx)	/* a[jt+p12] <- ~t7 */\n\t"\
		"vmovaps	%%ymm4,0x060(%%rcx)	/* a[jp+p4 ] <- ~t4 */				\n\t		vmovaps		%%ymm12,0x160(%%rcx)	/* a[jp+p4 ] <- ~t4 */\n\t"\
		"\n\t"\
"/**************************************************************************************/\n\t"\
"/*...and now do four more radix-4 transforms, including the internal twiddle factors: */\n\t"\
"/**************************************************************************************/\n\t"\
		"\n\t"\
"/*...Block 1: t1,9,17,25 */\n\t											/*...Block 3: t5,13,21,29: All rax-offsets incr +0x40 in rcol w.r.to lcol: */\n\t"\
		"movq	%[__r1],%%rax							\n\t				vmovaps		0x080(%%rax),%%ymm8 		/* t5  */\n\t"\
		"movq	%[__r9],%%rbx							\n\t				vmovaps		0x0a0(%%rax),%%ymm9 		/* t6  */\n\t"\
		"movq	%[__isrt2],%%rsi						\n\t				vmovaps		0x1a0(%%rax),%%ymm11		/* t14 */\n\t"\
		"vmovaps		 (%%rax),%%ymm0		/* t1  */\n\t					vmovaps		0x180(%%rax),%%ymm10		/* t13 */\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1		/* t2  */\n\t					vsubpd		%%ymm11,%%ymm8 ,%%ymm8 		/* t5 =t5 -t14 */\n\t"\
		"vmovaps	0x100(%%rax),%%ymm2		/* t9  */\n\t					vsubpd		%%ymm10,%%ymm9 ,%%ymm9 		/* t14=t6 -t13 */\n\t"\
		"vmovaps	0x120(%%rax),%%ymm3		/* t14 */\n\t					vaddpd		%%ymm11,%%ymm11,%%ymm11		/* 2.t14 */\n\t"\
		"\n\t																vaddpd		%%ymm10,%%ymm10,%%ymm10		/* 2.t13 */\n\t"\
		"vsubpd		%%ymm2,%%ymm0,%%ymm0		/* t9 =t1 -t9  */\n\t		vaddpd		%%ymm8 ,%%ymm11,%%ymm11		/* t13=t5 +t14 */\n\t"\
		"vsubpd		%%ymm3,%%ymm1,%%ymm1		/* t14=t2 -t14 */\n\t		vaddpd		%%ymm9 ,%%ymm10,%%ymm10		/* t6 =t6 +t13 */\n\t"\
		"vaddpd		%%ymm2,%%ymm2,%%ymm2		/* 2.t9  */\n\t				vmovaps		0x280(%%rax),%%ymm12		/* t21 */\n\t"\
		"vaddpd		%%ymm3,%%ymm3,%%ymm3		/* 2.t14 */\n\t				vmovaps		0x2a0(%%rax),%%ymm13		/* t22 */\n\t"\
		"vaddpd		%%ymm0,%%ymm2,%%ymm2		/* t1 =t1 +t9  */\n\t		vmovaps		0x380(%%rax),%%ymm14		/* t29 */\n\t"\
		"vaddpd		%%ymm1,%%ymm3,%%ymm3		/* t2 =t2 +t14 */\n\t		vmovaps		0x3a0(%%rax),%%ymm15		/* t30 */\n\t"\
		"\n\t																vsubpd		0x2a0(%%rax),%%ymm12,%%ymm12		/* t21-t22 */\n\t"\
		"vmovaps	0x200(%%rax),%%ymm4		/* t17 */\n\t					vaddpd		0x280(%%rax),%%ymm13,%%ymm13		/* t22+t21 */\n\t"\
		"vmovaps	0x220(%%rax),%%ymm5		/* t18 */\n\t					vmulpd		(%%rsi)		,%%ymm12,%%ymm12	/* t21 = (t21-t22)*ISRT2 */\n\t"\
		"vmovaps	0x300(%%rax),%%ymm6		/* t25 */\n\t					vmulpd		(%%rsi)		,%%ymm13,%%ymm13	/* t22 = (t22+t21)*ISRT2 */\n\t"\
		"vmovaps	0x320(%%rax),%%ymm7		/* t26 */\n\t					vaddpd		0x3a0(%%rax),%%ymm14,%%ymm14		/* t29+t30 */\n\t"\
		"\n\t																vsubpd		0x380(%%rax),%%ymm15,%%ymm15		/* t30-t29 */\n\t"\
		"vsubpd		%%ymm6,%%ymm4,%%ymm4		/* t25=t17-t25 */\n\t		vmulpd		(%%rsi),%%ymm14,%%ymm14	/*  rt = (t29+t30)*ISRT2 */\n\t"\
		"vsubpd		%%ymm7,%%ymm5,%%ymm5		/* t26=t18-t26 */\n\t		vmulpd		(%%rsi),%%ymm15,%%ymm15	/*  it = (t30-t29)*ISRT2 */\n\t"\
		"vaddpd		%%ymm6,%%ymm6,%%ymm6		/* 2.t25 */\n\t				vsubpd		%%ymm14,%%ymm12,%%ymm12		/* t21=t21-rt */\n\t"\
		"vaddpd		%%ymm7,%%ymm7,%%ymm7		/* 2.t26 */\n\t				vsubpd		%%ymm15,%%ymm13,%%ymm13		/* t22=t22-it */\n\t"\
		"vaddpd		%%ymm4,%%ymm6,%%ymm6		/* t17=t17+t25 */\n\t		vaddpd		%%ymm14,%%ymm14,%%ymm14		/*      2* rt */\n\t"\
		"vaddpd		%%ymm5,%%ymm7,%%ymm7		/* t18=t18+t26 */\n\t		vaddpd		%%ymm15,%%ymm15,%%ymm15		/*      2* it */\n\t"\
		"\n\t																vaddpd		%%ymm12,%%ymm14,%%ymm14		/* t29=t21+rt */\n\t"\
		"\n\t																vaddpd		%%ymm13,%%ymm15,%%ymm15		/* t30=t22+it */\n\t"\
		"movq	%[__r17],%%rcx						\n\t		movq		%[__r25],%%rdx			\n\t"\
		"vsubpd		%%ymm6,%%ymm2,%%ymm2						\n\t		vsubpd		%%ymm12,%%ymm8 ,%%ymm8 			\n\t"\
		"vsubpd		%%ymm7,%%ymm3,%%ymm3						\n\t		vsubpd		%%ymm13,%%ymm10,%%ymm10			\n\t"\
		"vaddpd		%%ymm6,%%ymm6,%%ymm6						\n\t		vaddpd		%%ymm12,%%ymm12,%%ymm12			\n\t"\
		"vaddpd		%%ymm7,%%ymm7,%%ymm7						\n\t		vaddpd		%%ymm13,%%ymm13,%%ymm13			\n\t"\
		"vmovaps		%%ymm2,     (%%rcx)						\n\t		vmovaps		%%ymm8 ,0x080(%%rcx)		\n\t"\
		"vmovaps		%%ymm3,0x020(%%rcx)						\n\t		vmovaps		%%ymm10,0x0a0(%%rcx)		\n\t"\
		"vaddpd		%%ymm2,%%ymm6,%%ymm6						\n\t		vaddpd		%%ymm8 ,%%ymm12,%%ymm12			\n\t"\
		"vaddpd		%%ymm3,%%ymm7,%%ymm7						\n\t		vaddpd		%%ymm10,%%ymm13,%%ymm13			\n\t"\
		"vmovaps		%%ymm6,     (%%rax)						\n\t		vmovaps		%%ymm12,0x080(%%rax)		\n\t"\
		"vmovaps		%%ymm7,0x020(%%rax)						\n\t		vmovaps		%%ymm13,0x0a0(%%rax)		\n\t"\
		"vsubpd		%%ymm5,%%ymm0,%%ymm0						\n\t		vsubpd		%%ymm15,%%ymm11,%%ymm11			\n\t"\
		"vsubpd		%%ymm4,%%ymm1,%%ymm1						\n\t		vsubpd		%%ymm14,%%ymm9 ,%%ymm9 			\n\t"\
		"vaddpd		%%ymm5,%%ymm5,%%ymm5						\n\t		vaddpd		%%ymm15,%%ymm15,%%ymm15			\n\t"\
		"vaddpd		%%ymm4,%%ymm4,%%ymm4						\n\t		vaddpd		%%ymm14,%%ymm14,%%ymm14			\n\t"\
		"vmovaps		%%ymm0,     (%%rbx)						\n\t		vmovaps		%%ymm11,0x080(%%rbx)		\n\t"\
		"vmovaps		%%ymm1,0x020(%%rdx)						\n\t		vmovaps		%%ymm9 ,0x0a0(%%rdx)		\n\t"\
		"vaddpd		%%ymm0,%%ymm5,%%ymm5						\n\t		vaddpd		%%ymm11,%%ymm15,%%ymm15			\n\t"\
		"vaddpd		%%ymm1,%%ymm4,%%ymm4						\n\t		vaddpd		%%ymm9 ,%%ymm14,%%ymm14			\n\t"\
		"vmovaps		%%ymm5,     (%%rdx)						\n\t		vmovaps		%%ymm15,0x080(%%rdx)		\n\t"\
		"vmovaps		%%ymm4,0x020(%%rbx)						\n\t		vmovaps		%%ymm14,0x0a0(%%rbx)		\n\t"\
		"\n\t"\
"/*...Block 2: t3,11,19,27 */\n\t"\
		"addq	$0x040,%%rax	/* r3  */\n\t"\
		"addq	$0x040,%%rbx							\n\t"\
		"addq	$0x040,%%rcx							\n\t"\
		"addq	$0x040,%%rdx							\n\t"\
		"leaq	0x020(%%rsi),%%rdi	/* cc0, from isrt2 */\n\t"\
		"\n\t																/*...Block 4: t7,15,23,31 */\n\t"\
		"vmovaps	0x200(%%rax),%%ymm4		/* t19 */			\n\t		vmovaps		0x280(%%rax),%%ymm12		/* t23 */\n\t"\
		"vmovaps	0x220(%%rax),%%ymm5		/* t20 */			\n\t		vmovaps		0x2a0(%%rax),%%ymm13		/* t24 */\n\t"\
		"vmovaps	0x300(%%rax),%%ymm6	/* t27 */\n\t						vmovaps		0x380(%%rax),%%ymm14		/* t31 */\n\t"\
		"vmovaps	0x320(%%rax),%%ymm7	/* t28 */\n\t						vmovaps		0x3a0(%%rax),%%ymm15		/* t32 */\n\t"\
		"vmovaps	0x200(%%rax),%%ymm0		/* copy t19 */		\n\t		vmovaps		0x280(%%rax),%%ymm8 		/* copy t23 */\n\t"\
		"vmovaps	0x220(%%rax),%%ymm1		/* copy t20 */		\n\t		vmovaps		0x2a0(%%rax),%%ymm9 		/* copy t24 */\n\t"\
		"vmovaps	0x300(%%rax),%%ymm2	/* copy t27 */\n\t					vmovaps		0x380(%%rax),%%ymm10		/* copy t31 */\n\t"\
		"vmovaps	0x320(%%rax),%%ymm3	/* copy t28 */\n\t					vmovaps		0x3a0(%%rax),%%ymm11		/* copy t32 */\n\t"\
		"\n\t"\
		"vmulpd		     (%%rdi),%%ymm4,%%ymm4		/* t19*c */	\n\t		vmulpd		0x020(%%rdi),%%ymm12,%%ymm12		/* t23*s */\n\t"\
		"vmulpd		0x020(%%rdi),%%ymm1,%%ymm1		/* t20*s */	\n\t		vmulpd		     (%%rdi),%%ymm9 ,%%ymm9 		/* t24*c */\n\t"\
		"vmulpd		0x020(%%rdi),%%ymm6,%%ymm6		/* t27*s */	\n\t		vmulpd		     (%%rdi),%%ymm14,%%ymm14		/* t31*c */\n\t"\
		"vmulpd		     (%%rdi),%%ymm3,%%ymm3		/* t28*c */	\n\t		vmulpd		0x020(%%rdi),%%ymm11,%%ymm11		/* t32*s */\n\t"\
		"vmulpd		     (%%rdi),%%ymm5,%%ymm5		/* t20*c */	\n\t		vmulpd		0x020(%%rdi),%%ymm13,%%ymm13		/* t24*s */\n\t"\
		"vmulpd		0x020(%%rdi),%%ymm7,%%ymm7		/* t28*s */	\n\t		vmulpd		     (%%rdi),%%ymm8 ,%%ymm8 		/* t23*c */\n\t"\
		"vmulpd		0x020(%%rdi),%%ymm0,%%ymm0		/* t19*s */	\n\t		vmulpd		     (%%rdi),%%ymm15,%%ymm15		/* t32*c */\n\t"\
		"vmulpd		     (%%rdi),%%ymm2,%%ymm2		/* t27*c */	\n\t		vmulpd		0x020(%%rdi),%%ymm10,%%ymm10		/* t31*s */\n\t"\
		"vsubpd		%%ymm1,%%ymm4,%%ymm4	/* ~t19 */			\n\t		vsubpd		%%ymm9 ,%%ymm12,%%ymm12		/* ~t23 */\n\t"\
		"vsubpd		%%ymm3,%%ymm6,%%ymm6	/* rt */\n\t					vaddpd		%%ymm8 ,%%ymm13,%%ymm13		/* ~t24 */\n\t"\
		"vaddpd		%%ymm0,%%ymm5,%%ymm5	/* ~t20 */			\n\t		vsubpd		%%ymm11,%%ymm14,%%ymm14		/* rt */\n\t"\
		"vaddpd		%%ymm2,%%ymm7,%%ymm7	/* it */\n\t					vaddpd		%%ymm10,%%ymm15,%%ymm15		/* it */\n\t"\
		"\n\t"\
		"vsubpd		%%ymm6,%%ymm4,%%ymm4		/*~t27=t19-rt */\n\t		vsubpd		%%ymm14,%%ymm12,%%ymm12		/*~t23=t23-rt */\n\t"\
		"vsubpd		%%ymm7,%%ymm5,%%ymm5		/*~t28=t20-it */\n\t		vsubpd		%%ymm15,%%ymm13,%%ymm13		/*~t24=t24-it */\n\t"\
		"vaddpd		%%ymm6,%%ymm6,%%ymm6		/*      2* rt */\n\t		vaddpd		%%ymm14,%%ymm14,%%ymm14		/*      2* rt */\n\t"\
		"vaddpd		%%ymm7,%%ymm7,%%ymm7		/*      2* it */\n\t		vaddpd		%%ymm15,%%ymm15,%%ymm15		/*      2* it */\n\t"\
		"vaddpd		%%ymm4,%%ymm6,%%ymm6		/*~t19=t19+rt */\n\t		vaddpd		%%ymm12,%%ymm14,%%ymm14		/*~t31=t23+rt */\n\t"\
		"vaddpd		%%ymm5,%%ymm7,%%ymm7		/*~t20=t20+it */\n\t		vaddpd		%%ymm13,%%ymm15,%%ymm15		/*~t32=t24+it */\n\t"\
		"\n\t"\
		"vmovaps	0x100(%%rax),%%ymm2		/* t11 */\n\t					vmovaps		0x180(%%rax),%%ymm10		/* t15 */\n\t"\
		"vmovaps	0x120(%%rax),%%ymm3		/* t12 */\n\t					vmovaps		0x1a0(%%rax),%%ymm11		/* t16 */\n\t"\
		"vsubpd		0x120(%%rax),%%ymm2,%%ymm2		/* t11-t12 */\n\t		vaddpd		0x1a0(%%rax),%%ymm10,%%ymm10	/* t15+t16 */\n\t"\
		"vaddpd		0x100(%%rax),%%ymm3,%%ymm3		/* t12+t11 */\n\t		vsubpd		0x180(%%rax),%%ymm11,%%ymm11	/* t16-t15 */\n\t"\
		"vmulpd		(%%rsi),%%ymm2,%%ymm2	/* rt = (t11-t12)*ISRT2 */\n\t	vmulpd		(%%rsi)		,%%ymm10,%%ymm10	/* rt = (t15+t16)*ISRT2 */\n\t"\
		"vmulpd		(%%rsi),%%ymm3,%%ymm3	/* it = (t12+t11)*ISRT2 */\n\t	vmulpd		(%%rsi)		,%%ymm11,%%ymm11	/* it = (t16-t15)*ISRT2 */\n\t"\
		"\n\t"\
		"vmovaps	     (%%rax),%%ymm0		/* t3  */\n\t					vmovaps		0x080(%%rax),%%ymm8 		/* t7  */\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1		/* t4  */\n\t					vmovaps		0x0a0(%%rax),%%ymm9 		/* t8  */\n\t"\
		"\n\t"\
		"vsubpd		      %%ymm2,%%ymm0,%%ymm0		/*~t11=t3 -rt */\n\t	vsubpd		     %%ymm10,%%ymm8 ,%%ymm8 	/*~t7 =t7 -rt */\n\t"\
		"vsubpd		      %%ymm3,%%ymm1,%%ymm1		/*~t12=t4 -it */\n\t	vsubpd		     %%ymm11,%%ymm9 ,%%ymm9 	/*~t8 =t8 -it */\n\t"\
		"vaddpd		     (%%rax),%%ymm2,%%ymm2		/*~t3 =rt +t3 */\n\t	vaddpd		0x080(%%rax),%%ymm10,%%ymm10	/*~t15=rt +t7 */\n\t"\
		"vaddpd		0x020(%%rax),%%ymm3,%%ymm3		/*~t4 =it +t4 */\n\t	vaddpd		0x0a0(%%rax),%%ymm11,%%ymm11	/*~t16=it +t8 */\n\t"\
		"\n\t"\
		"vsubpd		%%ymm6,%%ymm2,%%ymm2		/* t3 -t19 */\n\t			vsubpd		%%ymm12		,%%ymm8 ,%%ymm8 	/* t7 -t23 */\n\t"\
		"vsubpd		%%ymm7,%%ymm3,%%ymm3		/* t4 -t20 */\n\t			vsubpd		%%ymm13		,%%ymm9 ,%%ymm9 	/* t8 -t24 */\n\t"\
		"vaddpd		%%ymm6,%%ymm6,%%ymm6		/*   2*t19 */\n\t			vaddpd		%%ymm12		,%%ymm12,%%ymm12	/*   2*t23 */\n\t"\
		"vaddpd		%%ymm7,%%ymm7,%%ymm7		/*   2*t20 */\n\t			vaddpd		%%ymm13		,%%ymm13,%%ymm13	/*   2*t24 */\n\t"\
		"vmovaps		%%ymm2,     (%%rcx)						\n\t		vmovaps		%%ymm8 ,0x080(%%rcx)		\n\t"\
		"vmovaps		%%ymm3,0x020(%%rcx)						\n\t		vmovaps		%%ymm9 ,0x0a0(%%rcx)		\n\t"\
		"vaddpd		%%ymm2,%%ymm6,%%ymm6						\n\t		vaddpd		%%ymm8 ,%%ymm12,%%ymm12			\n\t"\
		"vaddpd		%%ymm3,%%ymm7,%%ymm7						\n\t		vaddpd		%%ymm9 ,%%ymm13,%%ymm13			\n\t"\
		"vmovaps		%%ymm6,     (%%rax)						\n\t		vmovaps		%%ymm12,0x080(%%rax)		\n\t"\
		"vmovaps		%%ymm7,0x020(%%rax)						\n\t		vmovaps		%%ymm13,0x0a0(%%rax)		\n\t"\
		"vsubpd		%%ymm5,%%ymm0,%%ymm0						\n\t		vsubpd		%%ymm15,%%ymm10,%%ymm10			\n\t"\
		"vsubpd		%%ymm4,%%ymm1,%%ymm1						\n\t		vsubpd		%%ymm14,%%ymm11,%%ymm11			\n\t"\
		"vaddpd		%%ymm5,%%ymm5,%%ymm5						\n\t		vaddpd		%%ymm15,%%ymm15,%%ymm15			\n\t"\
		"vaddpd		%%ymm4,%%ymm4,%%ymm4						\n\t		vaddpd		%%ymm14,%%ymm14,%%ymm14			\n\t"\
		"vmovaps		%%ymm0,     (%%rbx)						\n\t		vmovaps		%%ymm10,0x080(%%rbx)		\n\t"\
		"vmovaps		%%ymm1,0x020(%%rdx)						\n\t		vmovaps		%%ymm11,0x0a0(%%rdx)		\n\t"\
		"vaddpd		%%ymm0,%%ymm5,%%ymm5						\n\t		vaddpd		%%ymm10,%%ymm15,%%ymm15			\n\t"\
		"vaddpd		%%ymm1,%%ymm4,%%ymm4						\n\t		vaddpd		%%ymm11,%%ymm14,%%ymm14			\n\t"\
		"vmovaps		%%ymm5,     (%%rdx)						\n\t		vmovaps		%%ymm15,0x080(%%rdx)		\n\t"\
		"vmovaps		%%ymm4,0x020(%%rbx)						\n\t		vmovaps		%%ymm14,0x0a0(%%rbx)		\n\t"\
		 :					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__add1] "m" (Xadd1)\
		 ,[__add2] "m" (Xadd2)\
		 ,[__add3] "m" (Xadd3)\
		 ,[__r1] "m" (Xr1)\
		 ,[__r9] "m" (Xr9)\
		 ,[__r17] "m" (Xr17)\
		 ,[__r25] "m" (Xr25)\
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__cc0] "m" (Xcc0)\
		 ,[__c1] "m" (Xc1)\
		 ,[__c2] "m" (Xc2)\
		 ,[__c3] "m" (Xc3)\
		 ,[__c4] "m" (Xc4)\
		 ,[__c5] "m" (Xc5)\
		 ,[__c6] "m" (Xc6)\
		 ,[__c7] "m" (Xc7)\
		 ,[__c8] "m" (Xc8)\
		 ,[__c9] "m" (Xc9)\
		 ,[__c10] "m" (Xc10)\
		 ,[__c11] "m" (Xc11)\
		 ,[__c12] "m" (Xc12)\
		 ,[__c13] "m" (Xc13)\
		 ,[__c14] "m" (Xc14)\
		 ,[__c15] "m" (Xc15)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"		/* Clobbered registers */\
	);\
	}

	/****** To-Do: Need a 16-register version! ****************/
	#define SSE2_RADIX16_WRAPPER_DIT(Xadd0,Xadd1,Xadd2,Xadd3,Xr1,Xr9,Xr17,Xr25,Xisrt2,Xcc0,Xc1,Xc2,Xc3,Xc4,Xc5,Xc6,Xc7,Xc8,Xc9,Xc10,Xc11,Xc12,Xc13,Xc14,Xc15)\
	{\
	__asm__ volatile (\
		"/*...Block 1: r1,9,17,25 */				\n\t"\
		"movq	%[__r1],%%rax					\n\t"\
		"movq	%%rax,%%rbx						\n\t"\
		"movq	%%rax,%%rcx						\n\t"\
		"movq	%%rax,%%rdx						\n\t"\
		"addq	$0x200,%%rbx					\n\t"\
		"addq	$0x100,%%rcx					\n\t"\
		"addq	$0x300,%%rdx					\n\t"\
		"/* SSE2_RADIX4_DIT_IN_PLACE(): */		\n\t"\
		"vmovaps	     (%%rax),%%ymm0				\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1				\n\t"\
		"vmovaps	     (%%rax),%%ymm2				\n\t"\
		"vmovaps	0x020(%%rax),%%ymm3				\n\t"\
		"vaddpd	     (%%rbx),%%ymm0,%%ymm0				\n\t"\
		"vaddpd	0x020(%%rbx),%%ymm1,%%ymm1				\n\t"\
		"vsubpd	     (%%rbx),%%ymm2,%%ymm2				\n\t"\
		"vsubpd	0x020(%%rbx),%%ymm3,%%ymm3				\n\t"\
		"vmovaps	     (%%rcx),%%ymm4				\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5				\n\t"\
		"vmovaps	     (%%rcx),%%ymm6				\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm7				\n\t"\
		"vaddpd	     (%%rdx),%%ymm4,%%ymm4				\n\t"\
		"vaddpd	0x020(%%rdx),%%ymm5,%%ymm5				\n\t"\
		"vsubpd	     (%%rdx),%%ymm6,%%ymm6				\n\t"\
		"vsubpd	0x020(%%rdx),%%ymm7,%%ymm7				\n\t"\
		"vsubpd	      %%ymm4,%%ymm0,%%ymm0					\n\t"\
		"vsubpd	      %%ymm5,%%ymm1,%%ymm1					\n\t"\
		"vmovaps	%%ymm0,     (%%rbx)				\n\t"\
		"vmovaps	%%ymm1,0x020(%%rbx)				\n\t"\
		"vaddpd	      %%ymm4,%%ymm4,%%ymm4					\n\t"\
		"vaddpd	      %%ymm5,%%ymm5,%%ymm5					\n\t"\
		"vaddpd	      %%ymm0,%%ymm4,%%ymm4					\n\t"\
		"vaddpd	      %%ymm1,%%ymm5,%%ymm5					\n\t"\
		"vmovaps	%%ymm4,     (%%rax)				\n\t"\
		"vmovaps	%%ymm5,0x020(%%rax)				\n\t"\
		"vsubpd	      %%ymm7,%%ymm2,%%ymm2					\n\t"\
		"vsubpd	      %%ymm6,%%ymm3,%%ymm3					\n\t"\
		"vmovaps	%%ymm2,     (%%rdx)				\n\t"\
		"vmovaps	%%ymm3,0x020(%%rcx)				\n\t"\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7					\n\t"\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6					\n\t"\
		"vaddpd	      %%ymm2,%%ymm7,%%ymm7					\n\t"\
		"vaddpd	      %%ymm3,%%ymm6,%%ymm6					\n\t"\
		"vmovaps	%%ymm7,     (%%rcx)				\n\t"\
		"vmovaps	%%ymm6,0x020(%%rdx)				\n\t"\
		"/*...Block 2: r5,13,21,29 */				\n\t"\
		"addq	$0x080,%%rax					\n\t"\
		"addq	$0x080,%%rbx					\n\t"\
		"addq	$0x080,%%rcx					\n\t"\
		"addq	$0x080,%%rdx					\n\t"\
		"/* SSE2_RADIX4_DIT_IN_PLACE(): */		\n\t"\
		"vmovaps	     (%%rax),%%ymm0				\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1				\n\t"\
		"vmovaps	     (%%rax),%%ymm2				\n\t"\
		"vmovaps	0x020(%%rax),%%ymm3				\n\t"\
		"vaddpd	     (%%rbx),%%ymm0,%%ymm0				\n\t"\
		"vaddpd	0x020(%%rbx),%%ymm1,%%ymm1				\n\t"\
		"vsubpd	     (%%rbx),%%ymm2,%%ymm2				\n\t"\
		"vsubpd	0x020(%%rbx),%%ymm3,%%ymm3				\n\t"\
		"vmovaps	     (%%rcx),%%ymm4				\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5				\n\t"\
		"vmovaps	     (%%rcx),%%ymm6				\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm7				\n\t"\
		"vaddpd	     (%%rdx),%%ymm4,%%ymm4				\n\t"\
		"vaddpd	0x020(%%rdx),%%ymm5,%%ymm5				\n\t"\
		"vsubpd	     (%%rdx),%%ymm6,%%ymm6				\n\t"\
		"vsubpd	0x020(%%rdx),%%ymm7,%%ymm7				\n\t"\
		"vsubpd	      %%ymm4,%%ymm0,%%ymm0					\n\t"\
		"vsubpd	      %%ymm5,%%ymm1,%%ymm1					\n\t"\
		"vmovaps	%%ymm0,     (%%rbx)				\n\t"\
		"vmovaps	%%ymm1,0x020(%%rbx)				\n\t"\
		"vaddpd	      %%ymm4,%%ymm4,%%ymm4					\n\t"\
		"vaddpd	      %%ymm5,%%ymm5,%%ymm5					\n\t"\
		"vaddpd	      %%ymm0,%%ymm4,%%ymm4					\n\t"\
		"vaddpd	      %%ymm1,%%ymm5,%%ymm5					\n\t"\
		"vmovaps	%%ymm4,     (%%rax)				\n\t"\
		"vmovaps	%%ymm5,0x020(%%rax)				\n\t"\
		"vsubpd	      %%ymm7,%%ymm2,%%ymm2					\n\t"\
		"vsubpd	      %%ymm6,%%ymm3,%%ymm3					\n\t"\
		"vmovaps	%%ymm2,     (%%rdx)				\n\t"\
		"vmovaps	%%ymm3,0x020(%%rcx)				\n\t"\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7					\n\t"\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6					\n\t"\
		"vaddpd	      %%ymm2,%%ymm7,%%ymm7					\n\t"\
		"vaddpd	      %%ymm3,%%ymm6,%%ymm6					\n\t"\
		"vmovaps	%%ymm7,     (%%rcx)				\n\t"\
		"vmovaps	%%ymm6,0x020(%%rdx)				\n\t"\
		"/*...Block 3: r3,11,19,27 */				\n\t"\
		"subq	$0x040,%%rax					\n\t"\
		"subq	$0x040,%%rbx					\n\t"\
		"subq	$0x040,%%rcx					\n\t"\
		"subq	$0x040,%%rdx					\n\t"\
		"/* SSE2_RADIX4_DIT_IN_PLACE(): */		\n\t"\
		"vmovaps	     (%%rax),%%ymm0				\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1				\n\t"\
		"vmovaps	     (%%rax),%%ymm2				\n\t"\
		"vmovaps	0x020(%%rax),%%ymm3				\n\t"\
		"vaddpd	     (%%rbx),%%ymm0,%%ymm0				\n\t"\
		"vaddpd	0x020(%%rbx),%%ymm1,%%ymm1				\n\t"\
		"vsubpd	     (%%rbx),%%ymm2,%%ymm2				\n\t"\
		"vsubpd	0x020(%%rbx),%%ymm3,%%ymm3				\n\t"\
		"vmovaps	     (%%rcx),%%ymm4				\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5				\n\t"\
		"vmovaps	     (%%rcx),%%ymm6				\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm7				\n\t"\
		"vaddpd	     (%%rdx),%%ymm4,%%ymm4				\n\t"\
		"vaddpd	0x020(%%rdx),%%ymm5,%%ymm5				\n\t"\
		"vsubpd	     (%%rdx),%%ymm6,%%ymm6				\n\t"\
		"vsubpd	0x020(%%rdx),%%ymm7,%%ymm7				\n\t"\
		"vsubpd	      %%ymm4,%%ymm0,%%ymm0					\n\t"\
		"vsubpd	      %%ymm5,%%ymm1,%%ymm1					\n\t"\
		"vmovaps	%%ymm0,     (%%rbx)				\n\t"\
		"vmovaps	%%ymm1,0x020(%%rbx)				\n\t"\
		"vaddpd	      %%ymm4,%%ymm4,%%ymm4					\n\t"\
		"vaddpd	      %%ymm5,%%ymm5,%%ymm5					\n\t"\
		"vaddpd	      %%ymm0,%%ymm4,%%ymm4					\n\t"\
		"vaddpd	      %%ymm1,%%ymm5,%%ymm5					\n\t"\
		"vmovaps	%%ymm4,     (%%rax)				\n\t"\
		"vmovaps	%%ymm5,0x020(%%rax)				\n\t"\
		"vsubpd	      %%ymm7,%%ymm2,%%ymm2					\n\t"\
		"vsubpd	      %%ymm6,%%ymm3,%%ymm3					\n\t"\
		"vmovaps	%%ymm2,     (%%rdx)				\n\t"\
		"vmovaps	%%ymm3,0x020(%%rcx)				\n\t"\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7					\n\t"\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6					\n\t"\
		"vaddpd	      %%ymm2,%%ymm7,%%ymm7					\n\t"\
		"vaddpd	      %%ymm3,%%ymm6,%%ymm6					\n\t"\
		"vmovaps	%%ymm7,     (%%rcx)				\n\t"\
		"vmovaps	%%ymm6,0x020(%%rdx)				\n\t"\
		"/*...Block 4: r7,15,23,31 */				\n\t"\
		"addq	$0x080,%%rax					\n\t"\
		"addq	$0x080,%%rbx					\n\t"\
		"addq	$0x080,%%rcx					\n\t"\
		"addq	$0x080,%%rdx					\n\t"\
		"/* SSE2_RADIX4_DIT_IN_PLACE(): */		\n\t"\
		"vmovaps	     (%%rax),%%ymm0				\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1				\n\t"\
		"vmovaps	     (%%rax),%%ymm2				\n\t"\
		"vmovaps	0x020(%%rax),%%ymm3				\n\t"\
		"vaddpd	     (%%rbx),%%ymm0,%%ymm0				\n\t"\
		"vaddpd	0x020(%%rbx),%%ymm1,%%ymm1				\n\t"\
		"vsubpd	     (%%rbx),%%ymm2,%%ymm2				\n\t"\
		"vsubpd	0x020(%%rbx),%%ymm3,%%ymm3				\n\t"\
		"vmovaps	     (%%rcx),%%ymm4				\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5				\n\t"\
		"vmovaps	     (%%rcx),%%ymm6				\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm7				\n\t"\
		"vaddpd	     (%%rdx),%%ymm4,%%ymm4				\n\t"\
		"vaddpd	0x020(%%rdx),%%ymm5,%%ymm5				\n\t"\
		"vsubpd	     (%%rdx),%%ymm6,%%ymm6				\n\t"\
		"vsubpd	0x020(%%rdx),%%ymm7,%%ymm7				\n\t"\
		"vsubpd	      %%ymm4,%%ymm0,%%ymm0					\n\t"\
		"vsubpd	      %%ymm5,%%ymm1,%%ymm1					\n\t"\
		"vmovaps	%%ymm0,     (%%rbx)				\n\t"\
		"vmovaps	%%ymm1,0x020(%%rbx)				\n\t"\
		"vaddpd	      %%ymm4,%%ymm4,%%ymm4					\n\t"\
		"vaddpd	      %%ymm5,%%ymm5,%%ymm5					\n\t"\
		"vaddpd	      %%ymm0,%%ymm4,%%ymm4					\n\t"\
		"vaddpd	      %%ymm1,%%ymm5,%%ymm5					\n\t"\
		"vmovaps	%%ymm4,     (%%rax)				\n\t"\
		"vmovaps	%%ymm5,0x020(%%rax)				\n\t"\
		"vsubpd	      %%ymm7,%%ymm2,%%ymm2					\n\t"\
		"vsubpd	      %%ymm6,%%ymm3,%%ymm3					\n\t"\
		"vmovaps	%%ymm2,     (%%rdx)				\n\t"\
		"vmovaps	%%ymm3,0x020(%%rcx)				\n\t"\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7					\n\t"\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6					\n\t"\
		"vaddpd	      %%ymm2,%%ymm7,%%ymm7					\n\t"\
		"vaddpd	      %%ymm3,%%ymm6,%%ymm6					\n\t"\
		"vmovaps	%%ymm7,     (%%rcx)				\n\t"\
		"vmovaps	%%ymm6,0x020(%%rdx)				\n\t"\
		"\n\t"\
	"/***************************************************************************************************/\n\t"\
	"/*..and now do four more radix-4 transforms, including the internal and external twiddle factors.  */\n\t"\
	"/*  Write even-index 16-byte output pairs to a(j1), odd-index to a(j2), unpack same as on inputs.  */\n\t"\
	"/*  We do the last 2 radix-4 blocks first, to make the unpack-interleaving of outputs smoother.    */\n\t"\
	"/***************************************************************************************************/\n\t"\
	"/* Main-array addresses still in add0,1, no need to re-init. 		All rax-offsets incr +0x200 in rcol w.r.to lcol: */\n\t"\
		"/*...Block 3: t3,11,19,27 -> r9,13,11,15: */	\n\t				/*...Block 4: t7,15,23,31 -> r25,29,27,31: */\n\t"\
		"movq		%[__r9],%%rax				\n\t"\
		"movq		%[__isrt2],%%rbx			\n\t"\
		"movq		%[__cc0],%%rcx				\n\t"\
		"vmovaps	0x040(%%rax),%%ymm4			\n\t						vmovaps		0x240(%%rax),%%ymm12		\n\t"\
		"vmovaps	0x0c0(%%rax),%%ymm0			\n\t						vmovaps		0x2c0(%%rax),%%ymm8 		\n\t"\
		"vmovaps	0x060(%%rax),%%ymm5			\n\t						vmovaps		0x260(%%rax),%%ymm13		\n\t"\
		"vmovaps	0x0e0(%%rax),%%ymm1			\n\t						vmovaps		0x2e0(%%rax),%%ymm9 		\n\t"\
		"vmulpd		     (%%rcx),%%ymm4,%%ymm4	\n\t						vmulpd		0x020(%%rcx),%%ymm12,%%ymm12		\n\t"\
		"vmulpd		0x020(%%rcx),%%ymm0,%%ymm0	\n\t						vmulpd		     (%%rcx),%%ymm8 ,%%ymm8 		\n\t"\
		"vmulpd		     (%%rcx),%%ymm5,%%ymm5	\n\t						vmulpd		0x020(%%rcx),%%ymm13,%%ymm13		\n\t"\
		"vmulpd		0x020(%%rcx),%%ymm1,%%ymm1	\n\t						vmulpd		     (%%rcx),%%ymm9 ,%%ymm9 		\n\t"\
		"vmovaps	0x040(%%rax),%%ymm6			\n\t						vmovaps		0x240(%%rax),%%ymm14		\n\t"\
		"vmovaps	0x0c0(%%rax),%%ymm2			\n\t						vmovaps		0x2c0(%%rax),%%ymm10		\n\t"\
		"vmovaps	0x060(%%rax),%%ymm7			\n\t						vmovaps		0x260(%%rax),%%ymm15		\n\t"\
		"vmovaps	0x0e0(%%rax),%%ymm3			\n\t						vmovaps		0x2e0(%%rax),%%ymm11		\n\t"\
		"vmulpd		0x020(%%rcx),%%ymm6,%%ymm6	\n\t						vmulpd		     (%%rcx),%%ymm14,%%ymm14		\n\t"\
		"vmulpd		     (%%rcx),%%ymm2,%%ymm2	\n\t						vmulpd		0x020(%%rcx),%%ymm10,%%ymm10		\n\t"\
		"vmulpd		0x020(%%rcx),%%ymm7,%%ymm7	\n\t						vmulpd		     (%%rcx),%%ymm15,%%ymm15		\n\t"\
		"vmulpd		     (%%rcx),%%ymm3,%%ymm3	\n\t						vmulpd		0x020(%%rcx),%%ymm11,%%ymm11		\n\t"\
		"vsubpd		%%ymm6,%%ymm5,%%ymm5		\n\t						vsubpd		%%ymm14,%%ymm13,%%ymm13		\n\t"\
		"vsubpd		%%ymm2,%%ymm1,%%ymm1		\n\t						vsubpd		%%ymm10,%%ymm9 ,%%ymm9 		\n\t"\
		"vaddpd		%%ymm7,%%ymm4,%%ymm4		\n\t						vaddpd		%%ymm15,%%ymm12,%%ymm12		\n\t"\
		"vaddpd		%%ymm3,%%ymm0,%%ymm0		\n\t						vaddpd		%%ymm11,%%ymm8 ,%%ymm8 		\n\t"\
		"vmovaps	%%ymm5,%%ymm7				\n\t						vmovaps		%%ymm13,%%ymm15				\n\t"\
		"vmovaps	%%ymm4,%%ymm6				\n\t						vmovaps		%%ymm12,%%ymm14				\n\t"\
		"vaddpd		%%ymm0,%%ymm4,%%ymm4		\n\t						vaddpd		%%ymm8 ,%%ymm12,%%ymm12		\n\t"\
		"vaddpd		%%ymm1,%%ymm5,%%ymm5		\n\t						vaddpd		%%ymm9 ,%%ymm13,%%ymm13		\n\t"\
		"vsubpd		%%ymm0,%%ymm6,%%ymm6		\n\t						vsubpd		%%ymm8 ,%%ymm14,%%ymm14		\n\t"\
		"vsubpd		%%ymm1,%%ymm7,%%ymm7		\n\t						vsubpd		%%ymm9 ,%%ymm15,%%ymm15		\n\t"\
		"vmovaps	0x080(%%rax),%%ymm2			\n\t						vmovaps		0x280(%%rax),%%ymm10		\n\t"\
		"vmovaps	0x0a0(%%rax),%%ymm3			\n\t						vmovaps		0x2a0(%%rax),%%ymm11		\n\t"\
		"vmovaps	     (%%rax),%%ymm0			\n\t						vmovaps		0x200(%%rax),%%ymm8 		\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1			\n\t						vmovaps		0x220(%%rax),%%ymm9 		\n\t"\
		"vaddpd		0x0a0(%%rax),%%ymm2,%%ymm2	\n\t						vsubpd		0x2a0(%%rax),%%ymm10,%%ymm10		\n\t"\
		"vsubpd		0x080(%%rax),%%ymm3,%%ymm3	\n\t						vaddpd		0x280(%%rax),%%ymm11,%%ymm11		\n\t"\
		"vmulpd		     (%%rbx),%%ymm2,%%ymm2	\n\t						vmulpd		(%%rbx),%%ymm10,%%ymm10		\n\t"\
		"vmulpd		     (%%rbx),%%ymm3,%%ymm3	\n\t						vmulpd		(%%rbx),%%ymm11,%%ymm11		\n\t"\
		"vsubpd		%%ymm2,%%ymm0,%%ymm0		\n\t						vsubpd		%%ymm10,%%ymm8 ,%%ymm8 		\n\t"\
		"vsubpd		%%ymm3,%%ymm1,%%ymm1		\n\t						vsubpd		%%ymm11,%%ymm9 ,%%ymm9 		\n\t"\
		"vaddpd		%%ymm2,%%ymm2,%%ymm2		\n\t						vaddpd		%%ymm10,%%ymm10,%%ymm10		\n\t"\
		"vaddpd		%%ymm3,%%ymm3,%%ymm3		\n\t						vaddpd		%%ymm11,%%ymm11,%%ymm11		\n\t"\
		"vaddpd		%%ymm0,%%ymm2,%%ymm2		\n\t						vaddpd		%%ymm8 ,%%ymm10,%%ymm10		\n\t"\
		"vaddpd		%%ymm1,%%ymm3,%%ymm3		\n\t						vaddpd		%%ymm9 ,%%ymm11,%%ymm11		\n\t"\
		"addq	$0x240,%%rcx		/* c1 from cc0  ; c3  in rcol */\n\t	vsubpd		%%ymm14,%%ymm8 ,%%ymm8 		\n\t"\
		"leaq	0x2a0(%%rbx),%%rdx	/* c9 from isrt2; c11 in rcol */\n\t	vsubpd		%%ymm15,%%ymm9 ,%%ymm9 		\n\t"\
		/* Using the upper block(s) of the main array for temp-storage in the section below led to a nasty AVX bug to track down: */\
		/* In fermat-mod mode the 4 block addresses in ascending order are add0,1,2,3 with no 'gaps' between blocks, whereas for */\
		/* mersenne-mod the addresses in asc. order are add0,2,3,1 with a gap between contiguous-data-block pairs 0,2 and 3,1. Thus */\
		/* for fermat-mod we need [add2] as the base-address of the 'high-half' block for temp-storage; for mersenne-mod we need [add3]. */\
		/* In both cases we have that add2 < add3 so instead use (add2 - add1): > 0 for fermat-mod, < 0 for mersenne - to differentiate: */\
		"movq	%[__add2],%%rsi		\n\t"/* destroyable copy of add2 */\
		"movq	%[__add2],%%rbx		\n\t"\
		"subq	%[__add1],%%rsi		\n\t"/* rsi = (add2 - add1); if this yields a borrow (i.e. sets CF) it's mersenne, else fermat. */\
		"cmovcq %[__add3],%%rbx	\n\t" /* if CF set (i.e. h > l), copy source [add3] into dest (rbx), else leave dest = [add2]. */\
		/* rbx shared between rcol/lcol; rcx/rdx-offsets incr +0x80 in rcol for rest of block: */\
		"vsubpd		%%ymm4,%%ymm2,%%ymm2		\n\t						vaddpd		%%ymm14,%%ymm14,%%ymm14		\n\t"\
		"vsubpd		%%ymm5,%%ymm3,%%ymm3		\n\t						vaddpd		%%ymm15,%%ymm15,%%ymm15		\n\t"\
		"vaddpd		%%ymm4,%%ymm4,%%ymm4		\n\t						vaddpd		%%ymm8 ,%%ymm14,%%ymm14		\n\t"\
		"vaddpd		%%ymm5,%%ymm5,%%ymm5		\n\t						vaddpd		%%ymm9 ,%%ymm15,%%ymm15		\n\t"\
		"vaddpd		%%ymm2,%%ymm4,%%ymm4		\n\t						vmovaps		%%ymm8 ,0x200(%%rax)		\n\t"\
		"vaddpd		%%ymm3,%%ymm5,%%ymm5		\n\t						vmovaps		%%ymm9 ,0x220(%%rax)		\n\t"\
		"vmovaps	%%ymm2,     (%%rax)			\n\t						vmovaps		%%ymm14,%%ymm8 		\n\t"\
		"vmovaps	%%ymm3,0x020(%%rax)			\n\t						vmovaps		%%ymm15,%%ymm9 		\n\t"\
		"vmovaps	%%ymm4,%%ymm2				\n\t						vmulpd		0x100(%%rcx),%%ymm14,%%ymm14		\n\t"\
		"vmovaps	%%ymm5,%%ymm3				\n\t						vmulpd		0x100(%%rcx),%%ymm15,%%ymm15		\n\t"\
		"vmulpd		     (%%rcx),%%ymm4,%%ymm4	\n\t						vmulpd		0x120(%%rcx),%%ymm8 ,%%ymm8 		\n\t"\
		"vmulpd		     (%%rcx),%%ymm5,%%ymm5	\n\t						vmulpd		0x120(%%rcx),%%ymm9 ,%%ymm9 		\n\t"\
		"vmulpd		0x020(%%rcx),%%ymm2,%%ymm2	\n\t						vsubpd		%%ymm8      ,%%ymm15,%%ymm15		\n\t"\
		"vmulpd		0x020(%%rcx),%%ymm3,%%ymm3	\n\t						vaddpd		%%ymm9      ,%%ymm14,%%ymm14		\n\t"\
		"vsubpd		%%ymm2,%%ymm5,%%ymm5		\n\t						vmovaps		%%ymm15,0x060(%%rbx)		\n\t"\
		"vaddpd		%%ymm3,%%ymm4,%%ymm4		\n\t						vmovaps		%%ymm14,0x040(%%rbx)		\n\t"\
		"vmovaps	%%ymm5,0x020(%%rbx)			\n\t						vmovaps		0x200(%%rax),%%ymm14		\n\t"\
		"vmovaps	%%ymm4,     (%%rbx)			\n\t						vmovaps		0x220(%%rax),%%ymm15		\n\t"\
		"vmovaps		 (%%rax),%%ymm4			\n\t						vmovaps		%%ymm14,%%ymm8 		\n\t"\
		"vmovaps	0x020(%%rax),%%ymm5			\n\t						vmovaps		%%ymm15,%%ymm9 		\n\t"\
		"vmovaps	%%ymm4,%%ymm2				\n\t						vmulpd		0x100(%%rdx),%%ymm14,%%ymm14		\n\t"\
		"vmovaps	%%ymm5,%%ymm3				\n\t						vmulpd		0x100(%%rdx),%%ymm15,%%ymm15		\n\t"\
		"vmulpd		     (%%rdx),%%ymm4,%%ymm4	\n\t						vmulpd		0x120(%%rdx),%%ymm8 ,%%ymm8 		\n\t"\
		"vmulpd		     (%%rdx),%%ymm5,%%ymm5	\n\t						vmulpd		0x120(%%rdx),%%ymm9 ,%%ymm9 		\n\t"\
		"vmulpd		0x020(%%rdx),%%ymm2,%%ymm2	\n\t						vsubpd		%%ymm8      ,%%ymm15,%%ymm15		\n\t"\
		"vmulpd		0x020(%%rdx),%%ymm3,%%ymm3	\n\t						vaddpd		%%ymm9      ,%%ymm14,%%ymm14		\n\t"\
		"vsubpd		%%ymm2,%%ymm5,%%ymm5		\n\t						vmovaps		%%ymm15,0x160(%%rbx)		\n\t"\
		"vaddpd		%%ymm3,%%ymm4,%%ymm4		\n\t						vmovaps		%%ymm14,0x140(%%rbx)		\n\t"\
		"vmovaps	%%ymm5,0x120(%%rbx)			\n\t						vsubpd		%%ymm13,%%ymm10,%%ymm10		\n\t"\
		"vmovaps	%%ymm4,0x100(%%rbx)			\n\t						vsubpd		%%ymm12,%%ymm11,%%ymm11		\n\t"\
		"addq	$0x080,%%rcx	/* c5  from c1; c7  in rcol */\n\t			vaddpd		%%ymm13,%%ymm13,%%ymm13		\n\t"\
		"addq	$0x080,%%rdx	/* c13 from c9; c15 in rcol */\n\t			vaddpd		%%ymm12,%%ymm12,%%ymm12		\n\t"\
		"vsubpd		%%ymm7,%%ymm0,%%ymm0		\n\t						vaddpd		%%ymm10,%%ymm13,%%ymm13		\n\t"\
		"vsubpd		%%ymm6,%%ymm1,%%ymm1		\n\t						vaddpd		%%ymm11,%%ymm12,%%ymm12		\n\t"\
		"vaddpd		%%ymm7,%%ymm7,%%ymm7		\n\t						vmovaps		%%ymm13,%%ymm8 		\n\t"\
		"vaddpd		%%ymm6,%%ymm6,%%ymm6		\n\t						vmovaps		%%ymm11,%%ymm9 		\n\t"\
		"vaddpd		%%ymm0,%%ymm7,%%ymm7		\n\t						vmulpd		0x100(%%rcx),%%ymm13,%%ymm13		\n\t"\
		"vaddpd		%%ymm1,%%ymm6,%%ymm6		\n\t						vmulpd		0x100(%%rcx),%%ymm11,%%ymm11		\n\t"\
		"vmovaps	%%ymm7,%%ymm4				\n\t						vmulpd		0x120(%%rcx),%%ymm8 ,%%ymm8 		\n\t"\
		"vmovaps	%%ymm1,%%ymm5				\n\t						vmulpd		0x120(%%rcx),%%ymm9 ,%%ymm9 		\n\t"\
		"vmulpd		     (%%rcx),%%ymm7,%%ymm7	\n\t						vsubpd		%%ymm8      ,%%ymm11,%%ymm11		\n\t"\
		"vmulpd		     (%%rcx),%%ymm1,%%ymm1	\n\t						vaddpd		%%ymm9      ,%%ymm13,%%ymm13		\n\t"\
		"vmulpd		0x020(%%rcx),%%ymm4,%%ymm4	\n\t						vmovaps		%%ymm11,0x0e0(%%rbx)		\n\t"\
		"vmulpd		0x020(%%rcx),%%ymm5,%%ymm5	\n\t						vmovaps		%%ymm13,0x0c0(%%rbx)		\n\t"\
		"vsubpd		%%ymm4,%%ymm1,%%ymm1		\n\t						vmovaps		%%ymm10,%%ymm8 		\n\t"\
		"vaddpd		%%ymm5,%%ymm7,%%ymm7		\n\t						vmovaps		%%ymm12,%%ymm9 		\n\t"\
		"vmovaps	%%ymm1,0x0a0(%%rbx)			\n\t						vmulpd		0x100(%%rdx),%%ymm10,%%ymm10		\n\t"\
		"vmovaps	%%ymm7,0x080(%%rbx)			\n\t						vmulpd		0x100(%%rdx),%%ymm12,%%ymm12		\n\t"\
		"vmovaps	%%ymm0,%%ymm4				\n\t						vmulpd		0x120(%%rdx),%%ymm8 ,%%ymm8 		\n\t"\
		"vmovaps	%%ymm6,%%ymm5				\n\t						vmulpd		0x120(%%rdx),%%ymm9 ,%%ymm9 		\n\t"\
		"vmulpd		     (%%rdx),%%ymm0,%%ymm0	\n\t						vsubpd		%%ymm8      ,%%ymm12,%%ymm12		\n\t"\
		"vmulpd		     (%%rdx),%%ymm6,%%ymm6	\n\t						vaddpd		%%ymm9      ,%%ymm10,%%ymm10		\n\t"\
		"vmulpd		0x020(%%rdx),%%ymm4,%%ymm4	\n\t						vmovaps		%%ymm12,0x1e0(%%rbx)		\n\t"\
		"vmulpd		0x020(%%rdx),%%ymm5,%%ymm5	\n\t						vmovaps		%%ymm10,0x1c0(%%rbx)		\n\t"\
		"vsubpd		%%ymm4,%%ymm6,%%ymm6		\n\t						/*...Block 3: t5,13,21,29 -> r17,21,19,23: */		\n\t"\
		"vaddpd		%%ymm5,%%ymm0,%%ymm0		\n\t						movq	%[__r1],%%rax		/* r17 in rcol */		\n\t"\
		"vmovaps	%%ymm6,0x1a0(%%rbx)			\n\t						movq		%[__isrt2],%%rdi		\n\t"\
		"vmovaps	%%ymm0,0x180(%%rbx)			\n\t						vmovaps		(%%rdi),%%ymm10		\n\t"\
		"																	vmovaps		0x240(%%rax),%%ymm12		\n\t"\
		"/*...Block 1: t1,9,17,25 -> r1,5,3,7: */		\n\t				vmovaps		0x260(%%rax),%%ymm13		\n\t"\
		"vmovaps	     (%%rax),%%ymm0			\n\t						vmovaps		0x2c0(%%rax),%%ymm8 		\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1			\n\t						vmovaps		0x2e0(%%rax),%%ymm9 		\n\t"\
		"vmovaps	0x080(%%rax),%%ymm2			\n\t						vaddpd		0x260(%%rax),%%ymm12,%%ymm12		\n\t"\
		"vmovaps	0x0a0(%%rax),%%ymm3			\n\t						vsubpd		0x240(%%rax),%%ymm13,%%ymm13		\n\t"\
		"vsubpd		0x080(%%rax),%%ymm0,%%ymm0	\n\t						vsubpd		0x2e0(%%rax),%%ymm8 ,%%ymm8 		\n\t"\
		"vsubpd		0x0a0(%%rax),%%ymm1,%%ymm1	\n\t						vaddpd		0x2c0(%%rax),%%ymm9 ,%%ymm9 		\n\t"\
		"vaddpd		     (%%rax),%%ymm2,%%ymm2	\n\t						vmulpd		%%ymm10,%%ymm12,%%ymm12		\n\t"\
		"vaddpd		0x020(%%rax),%%ymm3,%%ymm3	\n\t						vmulpd		%%ymm10,%%ymm13,%%ymm13		\n\t"\
		"vmovaps	0x040(%%rax),%%ymm4			\n\t						vmulpd		%%ymm10,%%ymm8 ,%%ymm8 		\n\t"\
		"vmovaps	0x060(%%rax),%%ymm5			\n\t						vmulpd		%%ymm10,%%ymm9 ,%%ymm9 		\n\t"\
		"vmovaps	0x0c0(%%rax),%%ymm6			\n\t						vmovaps		%%ymm12,%%ymm14		\n\t"\
		"vmovaps	0x0e0(%%rax),%%ymm7			\n\t						vmovaps		%%ymm13,%%ymm15		\n\t"\
		"vsubpd		0x0c0(%%rax),%%ymm4,%%ymm4	\n\t						vsubpd		%%ymm8 ,%%ymm12,%%ymm12		\n\t"\
		"vsubpd		0x0e0(%%rax),%%ymm5,%%ymm5	\n\t						vsubpd		%%ymm9 ,%%ymm13,%%ymm13		\n\t"\
		"vaddpd		0x040(%%rax),%%ymm6,%%ymm6	\n\t						vaddpd		%%ymm8 ,%%ymm14,%%ymm14		\n\t"\
		"vaddpd		0x060(%%rax),%%ymm7,%%ymm7	\n\t						vaddpd		%%ymm9 ,%%ymm15,%%ymm15		\n\t"\
		"movq	%%rbx,%%rcx	\n\t"/* rbx saves 'high-half' block address; move that to rcx and put 'low-half' block address into rbx */\
		"movq	%[__add0],%%rbx					\n\t						vmovaps		0x200(%%rax),%%ymm8 		\n\t"\
		"subq	$0x280,%%rdx	/* c8 from c13 */\n\t						vmovaps		0x220(%%rax),%%ymm9 		\n\t"\
		"vaddpd		%%ymm6,%%ymm2,%%ymm2		\n\t						vmovaps		0x280(%%rax),%%ymm10		\n\t"\
		"vaddpd		%%ymm7,%%ymm3,%%ymm3		\n\t						vmovaps		0x2a0(%%rax),%%ymm11		\n\t"\
		"vmovaps	%%ymm2,     (%%rbx)			\n\t						vsubpd		0x2a0(%%rax),%%ymm8 ,%%ymm8 		\n\t"\
		"vmovaps	%%ymm3,0x020(%%rbx)			\n\t						vsubpd		0x280(%%rax),%%ymm9 ,%%ymm9 		\n\t"\
		"vaddpd		%%ymm6,%%ymm6,%%ymm6		\n\t						vaddpd		0x200(%%rax),%%ymm11,%%ymm11		\n\t"\
		"vaddpd		%%ymm7,%%ymm7,%%ymm7		\n\t						vaddpd		0x220(%%rax),%%ymm10,%%ymm10		\n\t"\
		"vsubpd		%%ymm6,%%ymm2,%%ymm2		\n\t						vsubpd		%%ymm12,%%ymm11,%%ymm11		\n\t"\
		"vsubpd		%%ymm7,%%ymm3,%%ymm3		\n\t						vsubpd		%%ymm13,%%ymm9 ,%%ymm9 		\n\t"\
		"vmovaps	%%ymm2,%%ymm6				\n\t						vaddpd		%%ymm12,%%ymm12,%%ymm12		\n\t"\
		"vmovaps	%%ymm3,%%ymm7				\n\t						vaddpd		%%ymm13,%%ymm13,%%ymm13		\n\t"\
		"vmulpd		     (%%rdx),%%ymm2,%%ymm2	\n\t						vaddpd		%%ymm11,%%ymm12,%%ymm12		\n\t"\
		"vmulpd		     (%%rdx),%%ymm3,%%ymm3	\n\t						vaddpd		%%ymm9 ,%%ymm13,%%ymm13		\n\t"\
		"vmulpd		0x020(%%rdx),%%ymm6,%%ymm6	\n\t						vmovaps		%%ymm11,     (%%rax)		\n\t"\
		"vmulpd		0x020(%%rdx),%%ymm7,%%ymm7	\n\t						vmovaps		%%ymm9 ,0x020(%%rax)		\n\t"\
		"vsubpd		%%ymm6,%%ymm3,%%ymm3		\n\t						vmovaps		%%ymm12,%%ymm11		\n\t"\
		"vaddpd		%%ymm7,%%ymm2,%%ymm2		\n\t						vmovaps		%%ymm13,%%ymm9 		\n\t"\
/*==========================*/"\n\t"\
		"																	vmulpd		0x0c0(%%rdx),%%ymm12,%%ymm12	/* c2 */		\n\t"\
		"																	vmulpd		0x0c0(%%rdx),%%ymm13,%%ymm13		\n\t"\
		"																	vmulpd		0x0e0(%%rdx),%%ymm11,%%ymm11		\n\t"\
		"																	vmulpd		0x0e0(%%rdx),%%ymm9 ,%%ymm9 		\n\t"\
		"																	vsubpd		%%ymm11     ,%%ymm13,%%ymm13		\n\t"\
		"																	vaddpd		%%ymm9      ,%%ymm12,%%ymm12		\n\t"\
"vmovaps	0x120(%%rcx),%%ymm7		\n\t"\
"vmovaps	0x100(%%rcx),%%ymm6		\n\t"\
"vmovaps	%%ymm3,0x120(%%rax)	/* r9 */\n\t"\
"vmovaps	%%ymm2,0x100(%%rax)	/* r8 */\n\t"\
"vmovaps	%%ymm7,0x320(%%rax)	/* r25 */\n\t"\
"vmovaps	%%ymm6,0x300(%%rax)	/* r24 */\n\t"\
																	"vmovaps	0x060(%%rcx),%%ymm11	\n\t"\
																	"vmovaps	0x040(%%rcx),%%ymm9		\n\t"\
																	"vmovaps	%%ymm13,0x060(%%rax)	/* r3 */\n\t"\
																	"vmovaps	%%ymm12,0x040(%%rax)	/* r2 */\n\t"\
																	"vmovaps	%%ymm11,0x260(%%rax)	/* r19 */\n\t"\
																	"vmovaps	%%ymm9 ,0x240(%%rax)	/* r18 */\n\t"\
		"																	addq	$0x040,%%rdx		/* c4 in lcol; c10 in rcol*/		\n\t"\
		"																	vmovaps		     (%%rax),%%ymm12		\n\t"\
		"																	vmovaps		0x020(%%rax),%%ymm13		\n\t"\
		"																	vmovaps		%%ymm12,%%ymm11		\n\t"\
		"																	vmovaps		%%ymm13,%%ymm9 		\n\t"\
		"																	vmulpd		0x0c0(%%rdx),%%ymm12,%%ymm12		\n\t"\
		"																	vmulpd		0x0c0(%%rdx),%%ymm13,%%ymm13		\n\t"\
		"																	vmulpd		0x0e0(%%rdx),%%ymm11,%%ymm11		\n\t"\
		"																	vmulpd		0x0e0(%%rdx),%%ymm9 ,%%ymm9 		\n\t"\
"vmovaps	0x020(%%rbx),%%ymm3		\n\t"\
"vmovaps	0x000(%%rbx),%%ymm2		\n\t"\
"vmovaps	0x020(%%rcx),%%ymm7		\n\t"\
"vmovaps	0x000(%%rcx),%%ymm6		\n\t"\
"vmovaps	%%ymm3,0x020(%%rax)	/* r1 */\n\t"\
"vmovaps	%%ymm2,0x000(%%rax)	/* r0 */\n\t"\
"vmovaps	%%ymm7,0x220(%%rax)	/* r17 */\n\t"\
"vmovaps	%%ymm6,0x200(%%rax)	/* r16 */\n\t"\
		"vaddpd		%%ymm5,%%ymm0,%%ymm0		\n\t						vsubpd		%%ymm11     ,%%ymm13,%%ymm13		\n\t"\
		"vsubpd		%%ymm4,%%ymm1,%%ymm1		\n\t						vaddpd		%%ymm9      ,%%ymm12,%%ymm12		\n\t"\
		"vmovaps	%%ymm0,%%ymm2				\n\t"\
		"vmovaps	%%ymm1,%%ymm3				\n\t"\
		"vaddpd		%%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vaddpd		%%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vmovaps	%%ymm0,%%ymm6				\n\t"\
		"vmovaps	%%ymm1,%%ymm7				\n\t"\
		"vmulpd		     (%%rdx),%%ymm2,%%ymm2	\n\t"\
		"vmulpd		     (%%rdx),%%ymm3,%%ymm3	\n\t"\
		"vmulpd		0x020(%%rdx),%%ymm6,%%ymm6	\n\t"\
		"vmulpd		0x020(%%rdx),%%ymm7,%%ymm7	\n\t"\
																	"vmovaps	0x160(%%rcx),%%ymm11	\n\t"\
																	"vmovaps	0x140(%%rcx),%%ymm9		\n\t"\
																	"vmovaps	%%ymm13,0x160(%%rax)	/* r11 */\n\t"\
																	"vmovaps	%%ymm12,0x140(%%rax)	/* r10 */\n\t"\
																	"vmovaps	%%ymm11,0x360(%%rax)	/* r27 */\n\t"\
																	"vmovaps	%%ymm9 ,0x340(%%rax)	/* r26 */\n\t"\
		"addq	$0x040,%%rdx		/* c12 in lcol; c6 in rcol*/\n\t		vsubpd		%%ymm15,%%ymm8 ,%%ymm8 		\n\t"\
		"vsubpd		%%ymm6,%%ymm3,%%ymm3		\n\t						vsubpd		%%ymm14,%%ymm10,%%ymm10		\n\t"\
		"vaddpd		%%ymm7,%%ymm2,%%ymm2		\n\t						vaddpd		%%ymm15,%%ymm15,%%ymm15		\n\t"\
		"																	vaddpd		%%ymm14,%%ymm14,%%ymm14		\n\t"\
		"																	vaddpd		%%ymm8 ,%%ymm15,%%ymm15		\n\t"\
		"																	vaddpd		%%ymm10,%%ymm14,%%ymm14		\n\t"\
		"																	vmovaps		%%ymm15,%%ymm12		\n\t"\
		"																	vmovaps		%%ymm10,%%ymm13		\n\t"\
		"																	vmulpd		0x0c0(%%rdx),%%ymm15,%%ymm15		\n\t"\
		"																	vmulpd		0x0c0(%%rdx),%%ymm10,%%ymm10		\n\t"\
		"																	vmulpd		0x0e0(%%rdx),%%ymm12,%%ymm12		\n\t"\
		"																	vmulpd		0x0e0(%%rdx),%%ymm13,%%ymm13		\n\t"\
		"																	vsubpd		%%ymm12     ,%%ymm10,%%ymm10		\n\t"\
"vmovaps	0x0a0(%%rcx),%%ymm7		\n\t"\
"vmovaps	0x080(%%rcx),%%ymm6		\n\t"\
"vmovaps	%%ymm3,0x0a0(%%rax)	/* r5 */\n\t"\
"vmovaps	%%ymm2,0x080(%%rax)	/* r4 */\n\t"\
"vmovaps	%%ymm7,0x2a0(%%rax)	/* r21 */\n\t"\
"vmovaps	%%ymm6,0x280(%%rax)	/* r20 */\n\t"\
		"vsubpd		%%ymm5,%%ymm0,%%ymm0		\n\t						vaddpd		%%ymm13     ,%%ymm15,%%ymm15		\n\t"\
		"vaddpd		%%ymm4,%%ymm1,%%ymm1		\n\t"\
		"vmovaps	%%ymm0,%%ymm6				\n\t"\
		"vmovaps	%%ymm1,%%ymm7				\n\t"\
		"vmulpd		     (%%rdx),%%ymm0,%%ymm0	\n\t"\
		"vmulpd		     (%%rdx),%%ymm1,%%ymm1	\n\t"\
		"vmulpd		0x020(%%rdx),%%ymm6,%%ymm6	\n\t"\
		"vmulpd		0x020(%%rdx),%%ymm7,%%ymm7	\n\t"\
		"vsubpd		%%ymm6,%%ymm1,%%ymm1		\n\t"\
		"vaddpd		%%ymm7,%%ymm0,%%ymm0		\n\t"\
																	"vmovaps	0x0e0(%%rcx),%%ymm13	\n\t"\
																	"vmovaps	0x0c0(%%rcx),%%ymm12	\n\t"\
																	"vmovaps	%%ymm10,0x0e0(%%rax)	/* r7 */\n\t"\
																	"vmovaps	%%ymm15,0x0c0(%%rax)	/* r6 */\n\t"\
																	"vmovaps	%%ymm13,0x2e0(%%rax)	/* r23 */\n\t"\
																	"vmovaps	%%ymm12,0x2c0(%%rax)	/* r22 */\n\t"\
		"																	vmovaps		%%ymm8 ,%%ymm12		\n\t"\
		"																	vmovaps		%%ymm14,%%ymm13		\n\t"\
		"																	vmulpd		0x100(%%rdx),%%ymm8 ,%%ymm8 	/* c14 */		\n\t"\
		"																	vmulpd		0x100(%%rdx),%%ymm14,%%ymm14		\n\t"\
		"																	vmulpd		0x120(%%rdx),%%ymm12,%%ymm12		\n\t"\
		"																	vmulpd		0x120(%%rdx),%%ymm13,%%ymm13		\n\t"\
		"																	vsubpd		%%ymm12     ,%%ymm14,%%ymm14		\n\t"\
		"																	vaddpd		%%ymm13     ,%%ymm8 ,%%ymm8 		\n\t"\
"vmovaps	0x1a0(%%rcx),%%ymm7				\n\t"\
"vmovaps	0x180(%%rcx),%%ymm6				\n\t"\
"vmovaps	%%ymm1,0x1a0(%%rax)	/* r13 */	\n\t"\
"vmovaps	%%ymm0,0x180(%%rax)	/* r12 */	\n\t"\
"vmovaps	%%ymm7,0x3a0(%%rax)	/* r29 */	\n\t"\
"vmovaps	%%ymm6,0x380(%%rax)	/* r28 */	\n\t"\
																	"vmovaps	0x1e0(%%rcx),%%ymm13	\n\t"\
																	"vmovaps	0x1c0(%%rcx),%%ymm12	\n\t"\
																	"vmovaps	%%ymm14,0x1e0(%%rax)	/* r15 */\n\t"\
																	"vmovaps	%%ymm8 ,0x1c0(%%rax)	/* r14 */\n\t"\
																	"vmovaps	%%ymm13,0x3e0(%%rax)	/* r31 */\n\t"\
																	"vmovaps	%%ymm12,0x3c0(%%rax)	/* r30 */\n\t"\
/*==========================*/"\n\t"\
	"/**** Finish with 4-way 'un'terleaving: ****/\n\t"\
		"movq	%[__r1] ,%%rsi	\n\t"\
		"movq	%[__add0],%%rax\n\t"\
		"movq	%[__add1],%%rbx\n\t"\
		"movq	%[__add2],%%rcx\n\t"\
		"movq	%[__add3],%%rdx\n\t"\
	"/* a[j+p0]: Inputs from r1 +0/1, 16/17, 2/3, 18/19. Outputs into [add0,add0+0x100,add1,add1+0x100]+0x0: */	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps	     (%%rsi),%%ymm1						\n\t		vmovaps	0x020(%%rsi),%%ymm3								\n\t"\
		"vmovaps	0x200(%%rsi),%%ymm5						\n\t		vmovaps	0x220(%%rsi),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps	0x040(%%rsi),%%ymm6						\n\t		vmovaps	0x060(%%rsi),%%ymm14							\n\t"\
		"vmovaps	0x240(%%rsi),%%ymm5						\n\t		vmovaps	0x260(%%rsi),%%ymm13							\n\t"\
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
/* DEBUG: Dump register contents back into above local-store slots: */\
"vmovaps	%%ymm0 ,0x000(%%rsi)	\n\t	vmovaps	%%ymm2 ,0x020(%%rsi)	\n\t"\
"vmovaps	%%ymm5 ,0x040(%%rsi)	\n\t	vmovaps	%%ymm13,0x060(%%rsi)	\n\t"\
"vmovaps	%%ymm1 ,0x200(%%rsi)	\n\t	vmovaps	%%ymm3 ,0x220(%%rsi)	\n\t"\
"vmovaps	%%ymm7 ,0x240(%%rsi)	\n\t	vmovaps	%%ymm15,0x260(%%rsi)	\n\t"\
	"/* a[j+p2]: Inputs from r3 +0/1, 16/17, 2/3, 18/19. Outputs into [add0,add0+0x100,add1,add1+0x100]+0x40: */		\n\t"\
		"addq	$0x80,%%rsi	\n\t"\
		"addq	$0x40,%%rax	\n\t"\
		"addq	$0x40,%%rbx	\n\t"\
		"addq	$0x40,%%rcx	\n\t"\
		"addq	$0x40,%%rdx	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps	     (%%rsi),%%ymm1						\n\t		vmovaps	0x020(%%rsi),%%ymm3								\n\t"\
		"vmovaps	0x200(%%rsi),%%ymm5						\n\t		vmovaps	0x220(%%rsi),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps	0x040(%%rsi),%%ymm6						\n\t		vmovaps	0x060(%%rsi),%%ymm14							\n\t"\
		"vmovaps	0x240(%%rsi),%%ymm5						\n\t		vmovaps	0x260(%%rsi),%%ymm13							\n\t"\
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
/* DEBUG: Dump register contents back into above local-store slots: */\
"vmovaps	%%ymm0 ,0x000(%%rsi)	\n\t	vmovaps	%%ymm2 ,0x020(%%rsi)	\n\t"\
"vmovaps	%%ymm5 ,0x040(%%rsi)	\n\t	vmovaps	%%ymm13,0x060(%%rsi)	\n\t"\
"vmovaps	%%ymm1 ,0x200(%%rsi)	\n\t	vmovaps	%%ymm3 ,0x220(%%rsi)	\n\t"\
"vmovaps	%%ymm7 ,0x240(%%rsi)	\n\t	vmovaps	%%ymm15,0x260(%%rsi)	\n\t"\
	"/* a[j+p4]: Inputs from r5 +0/1, 16/17, 2/3, 18/19. Outputs into [add0,add0+0x100,add1,add1+0x100]+0x80: */		\n\t"\
		"addq	$0x80,%%rsi	\n\t"\
		"addq	$0x40,%%rax	\n\t"\
		"addq	$0x40,%%rbx	\n\t"\
		"addq	$0x40,%%rcx	\n\t"\
		"addq	$0x40,%%rdx	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps	     (%%rsi),%%ymm1						\n\t		vmovaps	0x020(%%rsi),%%ymm3								\n\t"\
		"vmovaps	0x200(%%rsi),%%ymm5						\n\t		vmovaps	0x220(%%rsi),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps	0x040(%%rsi),%%ymm6						\n\t		vmovaps	0x060(%%rsi),%%ymm14							\n\t"\
		"vmovaps	0x240(%%rsi),%%ymm5						\n\t		vmovaps	0x260(%%rsi),%%ymm13							\n\t"\
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
/* DEBUG: Dump register contents back into above local-store slots: */\
"vmovaps	%%ymm0 ,0x000(%%rsi)	\n\t	vmovaps	%%ymm2 ,0x020(%%rsi)	\n\t"\
"vmovaps	%%ymm5 ,0x040(%%rsi)	\n\t	vmovaps	%%ymm13,0x060(%%rsi)	\n\t"\
"vmovaps	%%ymm1 ,0x200(%%rsi)	\n\t	vmovaps	%%ymm3 ,0x220(%%rsi)	\n\t"\
"vmovaps	%%ymm7 ,0x240(%%rsi)	\n\t	vmovaps	%%ymm15,0x260(%%rsi)	\n\t"\
	"/* a[j+p6]: Inputs from r7 +0/1, 16/17, 2/3, 18/19. Outputs into [add0,add0+0x100,add1,add1+0x100]+0xc0: */		\n\t"\
		"addq	$0x80,%%rsi	\n\t"\
		"addq	$0x40,%%rax	\n\t"\
		"addq	$0x40,%%rbx	\n\t"\
		"addq	$0x40,%%rcx	\n\t"\
		"addq	$0x40,%%rdx	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps	     (%%rsi),%%ymm1						\n\t		vmovaps	0x020(%%rsi),%%ymm3								\n\t"\
		"vmovaps	0x200(%%rsi),%%ymm5						\n\t		vmovaps	0x220(%%rsi),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps	0x040(%%rsi),%%ymm6						\n\t		vmovaps	0x060(%%rsi),%%ymm14							\n\t"\
		"vmovaps	0x240(%%rsi),%%ymm5						\n\t		vmovaps	0x260(%%rsi),%%ymm13							\n\t"\
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
/* DEBUG: Dump register contents back into above local-store slots: */\
"vmovaps	%%ymm0 ,0x000(%%rsi)	\n\t	vmovaps	%%ymm2 ,0x020(%%rsi)	\n\t"\
"vmovaps	%%ymm5 ,0x040(%%rsi)	\n\t	vmovaps	%%ymm13,0x060(%%rsi)	\n\t"\
"vmovaps	%%ymm1 ,0x200(%%rsi)	\n\t	vmovaps	%%ymm3 ,0x220(%%rsi)	\n\t"\
"vmovaps	%%ymm7 ,0x240(%%rsi)	\n\t	vmovaps	%%ymm15,0x260(%%rsi)	\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__add1] "m" (Xadd1)\
		 ,[__add2] "m" (Xadd2)\
		 ,[__add3] "m" (Xadd3)\
		 ,[__r1] "m" (Xr1)\
		 ,[__r9] "m" (Xr9)\
		 ,[__r17] "m" (Xr17)\
		 ,[__r25] "m" (Xr25)\
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__cc0] "m" (Xcc0)\
		 ,[__c1] "m" (Xc1)\
		 ,[__c2] "m" (Xc2)\
		 ,[__c3] "m" (Xc3)\
		 ,[__c4] "m" (Xc4)\
		 ,[__c5] "m" (Xc5)\
		 ,[__c6] "m" (Xc6)\
		 ,[__c7] "m" (Xc7)\
		 ,[__c8] "m" (Xc8)\
		 ,[__c9] "m" (Xc9)\
		 ,[__c10] "m" (Xc10)\
		 ,[__c11] "m" (Xc11)\
		 ,[__c12] "m" (Xc12)\
		 ,[__c13] "m" (Xc13)\
		 ,[__c14] "m" (Xc14)\
		 ,[__c15] "m" (Xc15)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
	);\
	}

	// Debug dump-registers macro: Insert body into needed location above
	#define FOO(Xaddr)\
	{\
	__asm__ volatile (\
		"\n\t"\
		"movq	%[__addr],%%rax	/* DEBUG: Dump register contents into trig-data slots: */\n\t"\
		"vmovaps	%%ymm0 ,0x000(%%rax)	\n\t"\
		"vmovaps	%%ymm1 ,0x020(%%rax)	\n\t"\
		"vmovaps	%%ymm2 ,0x040(%%rax)	\n\t"\
		"vmovaps	%%ymm3 ,0x060(%%rax)	\n\t"\
		"vmovaps	%%ymm4 ,0x080(%%rax)	\n\t"\
		"vmovaps	%%ymm5 ,0x0a0(%%rax)	\n\t"\
		"vmovaps	%%ymm6 ,0x0c0(%%rax)	\n\t"\
		"vmovaps	%%ymm7 ,0x0e0(%%rax)	\n\t"\
		"vmovaps	%%ymm8 ,0x100(%%rax)	\n\t"\
		"vmovaps	%%ymm9 ,0x120(%%rax)	\n\t"\
		"vmovaps	%%ymm10,0x140(%%rax)	\n\t"\
		"vmovaps	%%ymm11,0x160(%%rax)	\n\t"\
		"vmovaps	%%ymm12,0x180(%%rax)	\n\t"\
		"vmovaps	%%ymm13,0x1a0(%%rax)	\n\t"\
		"vmovaps	%%ymm14,0x1c0(%%rax)	\n\t"\
		"vmovaps	%%ymm15,0x1e0(%%rax)	\n\t"\
		:					/* outputs: none */\
		: [__addr] "m" (Xaddr)	/* All inputs from memory addresses here */\
		: "cc","memory","rax","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
	);\
	}

#elif defined(USE_SSE2)

  #define USE_64BIT_ASM_STYLE	1	// My x86 timings indicate fancier versions below using xmm0-15 for the radix-32 DFT is faster.

  #if !USE_64BIT_ASM_STYLE	// Default is to use simple 64-bit-ified version of the analogous 32-bit ASM macros, i.e. using just xmm0-7.

	#define SSE2_RADIX16_WRAPPER_DIF(Xadd0,Xadd1,Xr1,Xr9,Xr17,Xr25,Xisrt2,Xcc0,Xc1,Xc2,Xc3,Xc4,Xc5,Xc6,Xc7,Xc8,Xc9,Xc10,Xc11,Xc12,Xc13,Xc14,Xc15)\
	{\
	__asm__ volatile (\
	"/*************************************************************/\n\t"\
	"/*                  1st set of inputs:                       */\n\t"\
	"/*************************************************************/\n\t"\
	"/*...Block 1: */\n\t"\
		"movq		%[__add0],%%rax\n\t"\
		"movq		%[__add1],%%rbx\n\t"\
		"movq		%[__r1] ,%%rcx\n\t"\
		"movq		%[__c4] ,%%rdx\n\t"\
		"\n\t"\
	"/* Real parts: */\n\t"\
		"movaps		0x40(%%rax)	,%%xmm6		/* a[j1+p4], this is the scratch xmm register  */\n\t"\
		"movaps		0x40(%%rbx)	,%%xmm2		/* a[j2+p4] */\n\t"\
		"movaps		%%xmm6		,%%xmm0		/* a[j1+p4] copy, his is the active  xmm register */\n\t"\
		"movaps		%%xmm2		,%%xmm3		/* a[j2+p4] copy */\n\t"\
		"unpckhpd	%%xmm2,%%xmm6\n\t"\
		"unpcklpd	%%xmm3,%%xmm0\n\t"\
		"movaps		%%xmm6		,0x140(%%rcx)	/* Store hi real in t21 */\n\t"\
	"/* Imag parts: */\n\t"\
		"movaps		0x50(%%rax)	,%%xmm7\n\t"\
		"movaps		0x50(%%rbx)	,%%xmm4\n\t"\
		"movaps		%%xmm7,%%xmm1\n\t"\
		"movaps		%%xmm4,%%xmm5\n\t"\
		"unpckhpd	%%xmm4,%%xmm7\n\t"\
		"unpcklpd	%%xmm5,%%xmm1\n\t"\
		"movaps		%%xmm7		,0x150(%%rcx)	/* Store hi imag in t22 */\n\t"\
		"\n\t"\
		"movaps		%%xmm0		,%%xmm2	/* xmm2 <- cpy a[jt+p4] */\n\t"\
		"movaps		%%xmm1		,%%xmm3	/* xmm3 <- cpy a[jp+p4] */\n\t"\
	"/*************************************************************************************/\n\t"\
	"/******** From here on, things are identical to the code in radix16_dif_pass: ********/\n\t"\
	"/*************************************************************************************/\n\t"\
		"mulpd		    (%%rdx)	,%%xmm0		/* a[jt+p4]*c4 */\n\t"\
		"mulpd		    (%%rdx)	,%%xmm1		/* a[jp+p4]*c4 */\n\t"\
		"mulpd		0x10(%%rdx)	,%%xmm2		/* a[jt+p4]*s4 */\n\t"\
		"mulpd		0x10(%%rdx)	,%%xmm3		/* a[jp+p4]*s4 */\n\t"\
		"addpd		%%xmm2		,%%xmm1	/* xmm1 <- t6 */\n\t"\
		"subpd		%%xmm3		,%%xmm0	/* xmm0 <- t5 */\n\t"\
		"movaps		%%xmm1		,%%xmm3	/* xmm3 <- cpy t6 */\n\t"\
		"movaps		%%xmm0		,%%xmm2	/* xmm2 <- cpy t5 */\n\t"\
		"\n\t"\
		"movq		%[__c12],%%rdx\n\t"\
	"/* Real parts: */\n\t"\
		"movaps		0xc0(%%rax)	,%%xmm6		/* a[j1+p12], this is the scratch xmm register  */\n\t"\
		"movaps		0xc0(%%rax)	,%%xmm4		/* a[j1+p12], this is the active  xmm register */\n\t"\
		"unpckhpd	0xc0(%%rbx)	,%%xmm6		/* a[j2+p12] gets read twice */\n\t"\
		"unpcklpd	0xc0(%%rbx)	,%%xmm4		/* a[jt+p12] */\n\t"\
		"movaps		%%xmm6		,0x160(%%rcx)	/* Store hi real in t23 */\n\t"\
	"/* Imag parts: */\n\t"\
		"movaps		0xd0(%%rax)	,%%xmm7\n\t"\
		"movaps		0xd0(%%rax)	,%%xmm5\n\t"\
		"unpckhpd	0xd0(%%rbx)	,%%xmm7\n\t"\
		"unpcklpd	0xd0(%%rbx)	,%%xmm5		/* a[jp+p12] */\n\t"\
		"movaps		%%xmm7		,0x170(%%rcx)	/* Store hi imag in t24 */\n\t"\
		"\n\t"\
		"movaps		%%xmm4		,%%xmm6	/* xmm6 <- cpy a[jt+p12] */\n\t"\
		"movaps		%%xmm5		,%%xmm7	/* xmm7 <- cpy a[jp+p12] */\n\t"\
		"\n\t"\
		"mulpd		    (%%rdx)	,%%xmm4		/* a[jt+p12]*c12 */\n\t"\
		"mulpd		    (%%rdx)	,%%xmm5		/* a[jp+p12]*c12 */\n\t"\
		"mulpd		0x10(%%rdx)	,%%xmm6		/* a[jt+p12]*s12 */\n\t"\
		"mulpd		0x10(%%rdx)	,%%xmm7		/* a[jp+p12]*s12 */\n\t"\
		"addpd		%%xmm6		,%%xmm5	/* xmm5 <- it */\n\t"\
		"subpd		%%xmm7		,%%xmm4	/* xmm4 <- rt		xmm6,7 free */\n\t"\
		"\n\t"\
		"addpd		%%xmm4		,%%xmm0	/* ~t5 <- t5 +rt */\n\t"\
		"addpd		%%xmm5		,%%xmm1	/* ~t6 <- t6 +it */\n\t"\
		"subpd		%%xmm4		,%%xmm2	/* ~t7 <- t5 -rt */\n\t"\
		"subpd		%%xmm5		,%%xmm3	/* ~t8 <- t6 -it	xmm4,5 free */\n\t"\
		"\n\t"\
		"/* Now do the p0,8 combo: */\n\t"\
		"movq		%[__c8] ,%%rdx\n\t"\
	"/* Real parts: */\n\t"\
		"movaps		0x80(%%rax)	,%%xmm6		/* a[j1+p8 ], this is the scratch xmm register  */\n\t"\
		"movaps		0x80(%%rax)	,%%xmm4		/* a[j1+p8 ], this is the active  xmm register */\n\t"\
		"unpckhpd	0x80(%%rbx)	,%%xmm6		/* a[j2+p8 ] gets read twice */\n\t"\
		"unpcklpd	0x80(%%rbx)	,%%xmm4		/* a[jt+p8 ] */\n\t"\
		"movaps			%%xmm6	,0x120(%%rcx)	/* Store hi real in t19 */\n\t"\
	"/* Imag parts: */\n\t"\
		"movaps		0x90(%%rax)	,%%xmm7\n\t"\
		"movaps		0x90(%%rax)	,%%xmm5\n\t"\
		"unpckhpd	0x90(%%rbx)	,%%xmm7\n\t"\
		"unpcklpd	0x90(%%rbx)	,%%xmm5		/* a[jp+p8 ] */\n\t"\
		"movaps		%%xmm7		,0x130(%%rcx)	/* Store hi imag in t20 */\n\t"\
		"\n\t"\
		"movaps		%%xmm4		,%%xmm6	/* xmm6 <- cpy a[jt+p8] */\n\t"\
		"movaps		%%xmm5		,%%xmm7	/* xmm7 <- cpy a[jp+p8] */\n\t"\
		"\n\t"\
		"mulpd		    (%%rdx)	,%%xmm4		/* a[jt+p8]*c8 */\n\t"\
		"mulpd		    (%%rdx)	,%%xmm5		/* a[jp+p8]*c8 */\n\t"\
		"mulpd		0x10(%%rdx)	,%%xmm6		/* a[jt+p8]*s8 */\n\t"\
		"mulpd		0x10(%%rdx)	,%%xmm7		/* a[jp+p8]*s8 */\n\t"\
		"addpd		%%xmm6		,%%xmm5	/* xmm5 <- it */\n\t"\
		"subpd		%%xmm7		,%%xmm4	/* xmm4 <- rt 	xmm6,7 free - stick t1,t2 in those */\n\t"\
		"\n\t"\
	"/* Real parts: */\n\t"\
		"movaps		    (%%rax)	,%%xmm6		/* a[j1    ], this is the scratch xmm register  */\n\t"\
		"movaps		    (%%rax)	,%%xmm7		/* a[j1    ], this is the active  xmm register */\n\t"\
		"unpckhpd	    (%%rbx)	,%%xmm6		/* a[j2    ] gets read twice */\n\t"\
		"unpcklpd	    (%%rbx)	,%%xmm7		/* a[jt] = t1*/\n\t"\
		"movaps		%%xmm6		,0x100(%%rcx)	/* Store hi real in t17 */\n\t"\
		"movaps		%%xmm7		,     (%%rcx)	/* Store active  in t1  */\n\t"\
	"/* Imag parts: */\n\t"\
		"movaps		0x10(%%rax)	,%%xmm6\n\t"\
		"movaps		0x10(%%rax)	,%%xmm7\n\t"\
		"unpckhpd	0x10(%%rbx)	,%%xmm6\n\t"\
		"unpcklpd	0x10(%%rbx)	,%%xmm7		/* a[jp] = t2*/\n\t"\
		"movaps		%%xmm6		,0x110(%%rcx)	/* Store hi imag in t18... */\n\t"\
		"movaps		    (%%rcx)	,%%xmm6		/* ...and reload t1. */\n\t"\
		"\n\t"\
		"subpd		%%xmm4		,%%xmm6	/* ~t3 <- t1 -rt */\n\t"\
		"subpd		%%xmm5		,%%xmm7	/* ~t4 <- t2 -it */\n\t"\
		"addpd		%%xmm4		,%%xmm4	/*          2*rt */\n\t"\
		"addpd		%%xmm5		,%%xmm5	/*          2*it */\n\t"\
		"addpd		%%xmm6		,%%xmm4	/* ~t1 <- t1 +rt */\n\t"\
		"addpd		%%xmm7		,%%xmm5	/* ~t2 <- t2 +it	xmm4,5 free */\n\t"\
		"\n\t"\
	"/* Finish radix-4 butterfly and store results into temporary-array slots: */\n\t"\
		"movq		%[__r1] ,%%rax\n\t"\
	"/*\n\t"\
	"~t5 =t1 -t5;		~t1 =t1 +t5;\n\t"\
	"~t6 =t2 -t6;		~t2 =t2 +t6;\n\t"\
	"*/\n\t"\
		"subpd		%%xmm0		,%%xmm4	/*~t5 =t1 -t5 */\n\t"\
		"subpd		%%xmm1		,%%xmm5	/*~t6 =t2 -t6 */\n\t"\
		"movaps		%%xmm4		,0x040(%%rax)	/* a[jt+p8 ] <- ~t5 */\n\t"\
		"movaps		%%xmm5		,0x050(%%rax)	/* a[jp+p8 ] <- ~t6 */\n\t"\
		"addpd		%%xmm0		,%%xmm0			/* 2*t5 */\n\t"\
		"addpd		%%xmm1		,%%xmm1			/* 2*t6 */\n\t"\
		"addpd		%%xmm4		,%%xmm0			/*~t1 =t1 +t5 */\n\t"\
		"addpd		%%xmm5		,%%xmm1			/*~t2 =t2 +t6 */\n\t"\
		"movaps		%%xmm0		,     (%%rax)	/* a[jt    ] <- ~t1 */\n\t"\
		"movaps		%%xmm1		,0x010(%%rax)	/* a[jp    ] <- ~t2 */\n\t"\
		"\n\t"\
	"/*\n\t"\
	"~t7 =t3 +t8;		~t3 =t3 -t8;\n\t"\
	"~t8 =t4 -t7;		~t4 =t4 +t7;\n\t"\
	"*/\n\t"\
		"subpd		%%xmm3		,%%xmm6	/*~t3 =t3 -t8 */\n\t"\
		"subpd		%%xmm2		,%%xmm7	/*~t8 =t4 -t7 */\n\t"\
		"movaps		%%xmm6		,0x020(%%rax)	/* a[jt+p4 ] <- ~t3 */\n\t"\
		"movaps		%%xmm7		,0x070(%%rax)	/* a[jp+p12] <- ~t8 */\n\t"\
		"addpd		%%xmm3		,%%xmm3			/* 2*t8 */\n\t"\
		"addpd		%%xmm2		,%%xmm2			/* 2*t7 */\n\t"\
		"addpd		%%xmm6		,%%xmm3			/*~t7 =t3 +t8 */\n\t"\
		"addpd		%%xmm7		,%%xmm2			/*~t4 =t4 +t7 */\n\t"\
		"movaps		%%xmm3		,0x060(%%rax)	/* a[jt+p12] <- ~t7 */\n\t"\
		"movaps		%%xmm2		,0x030(%%rax)	/* a[jp+p4 ] <- ~t4 */\n\t"\
		"\n\t"\
"/*...Block 2:		Cost: 46 MOVapd, 16 UNPCKHPD, 28 ADD/SUBpd, 16 MULpd */\n\t"\
		"\n\t"\
		"movq		%[__add0],%%rax\n\t"\
		"movq		%[__add1],%%rbx\n\t"\
		"movq		%[__r9] ,%%rcx\n\t"\
		"\n\t"\
"/* Do the p2,10 combo: */\n\t"\
		"movq		%[__c2] ,%%rdx\n\t"\
	"/* Real parts: */\n\t"\
		"movaps		0x20(%%rax)	,%%xmm6		/* a[j1+p2 ], this is the scratch xmm register */\n\t"\
		"movaps		0x20(%%rax)	,%%xmm0		/* a[j1+p2 ], this is the active  xmm register */\n\t"\
		"unpckhpd	0x20(%%rbx)	,%%xmm6		/* a[j2+p2 ] gets read twice */\n\t"\
		"unpcklpd	0x20(%%rbx)	,%%xmm0		/* a[jt+p2 ] */\n\t"\
		"movaps			%%xmm6	,0x100(%%rcx)	/* Store hi real in t9 +16 */\n\t"\
	"/* Imag parts: */\n\t"\
		"movaps		0x30(%%rax)	,%%xmm7\n\t"\
		"movaps		0x30(%%rax)	,%%xmm1\n\t"\
		"unpckhpd	0x30(%%rbx)	,%%xmm7\n\t"\
		"unpcklpd	0x30(%%rbx)	,%%xmm1		/* a[jp+p2 ] */\n\t"\
		"movaps		%%xmm7		,0x110(%%rcx)	/* Store hi imag in t10+16 */\n\t"\
		"\n\t"\
		"movaps		%%xmm0		,%%xmm2	/* xmm2 <- cpy a[jt+p2] */\n\t"\
		"movaps		%%xmm1		,%%xmm3	/* xmm3 <- cpy a[jp+p2] */\n\t"\
		"\n\t"\
		"mulpd		    (%%rdx)	,%%xmm0		/* a[jt+p2]*c2 */\n\t"\
		"mulpd		    (%%rdx)	,%%xmm1		/* a[jp+p2]*c2 */\n\t"\
		"mulpd		0x10(%%rdx)	,%%xmm2		/* a[jt+p2]*s2 */\n\t"\
		"mulpd		0x10(%%rdx)	,%%xmm3		/* a[jp+p2]*s2 */\n\t"\
		"addpd		%%xmm2		,%%xmm1	/* xmm1 <- t10*/\n\t"\
		"subpd		%%xmm3		,%%xmm0	/* xmm0 <- t9 */\n\t"\
		"movaps		%%xmm1		,%%xmm3	/* xmm3 <- cpy t10*/\n\t"\
		"movaps		%%xmm0		,%%xmm2	/* xmm2 <- cpy t9 */\n\t"\
		"\n\t"\
		"movq		%[__c10],%%rdx\n\t"\
	"/* Real parts: */\n\t"\
		"movaps		0xa0(%%rax)	,%%xmm6		/* a[j1+p10], this is the scratch xmm register  */\n\t"\
		"movaps		0xa0(%%rax)	,%%xmm4		/* a[j1+p10], this is the active  xmm register */\n\t"\
		"unpckhpd	0xa0(%%rbx)	,%%xmm6		/* a[j2+p10] gets read twice */\n\t"\
		"unpcklpd	0xa0(%%rbx)	,%%xmm4		/* a[jt+p10] */\n\t"\
		"movaps			%%xmm6	,0x120(%%rcx)	/* Store hi real in t11+16 */\n\t"\
	"/* Imag parts: */\n\t"\
		"movaps		0xb0(%%rax)	,%%xmm7\n\t"\
		"movaps		0xb0(%%rax)	,%%xmm5\n\t"\
		"unpckhpd	0xb0(%%rbx)	,%%xmm7\n\t"\
		"unpcklpd	0xb0(%%rbx)	,%%xmm5		/* a[jp+p10] */\n\t"\
		"movaps			%%xmm7	,0x130(%%rcx)	/* Store hi imag in t12+16 */\n\t"\
		"\n\t"\
		"movaps			%%xmm4	,%%xmm6	/* xmm6 <- cpy a[jt+p10] */\n\t"\
		"movaps			%%xmm5	,%%xmm7	/* xmm7 <- cpy a[jp+p10] */\n\t"\
		"\n\t"\
		"mulpd		    (%%rdx)	,%%xmm4		/* a[jt+p10]*c10 */\n\t"\
		"mulpd		    (%%rdx)	,%%xmm5		/* a[jp+p10]*c10 */\n\t"\
		"mulpd		0x10(%%rdx)	,%%xmm6		/* a[jt+p10]*s10 */\n\t"\
		"mulpd		0x10(%%rdx)	,%%xmm7		/* a[jp+p10]*s10 */\n\t"\
		"addpd		%%xmm6		,%%xmm5	/* xmm5 <- it */\n\t"\
		"subpd		%%xmm7		,%%xmm4	/* xmm4 <- rt		xmm6,7 free */\n\t"\
		"\n\t"\
		"addpd		%%xmm4		,%%xmm0	/* ~t13<- t13+rt */\n\t"\
		"addpd		%%xmm5		,%%xmm1	/* ~t14<- t14+it */\n\t"\
		"subpd		%%xmm4		,%%xmm2	/* ~t15<- t13-rt */\n\t"\
		"subpd		%%xmm5		,%%xmm3	/* ~t16<- t14-it	xmm4,5 free */\n\t"\
		"\n\t"\
"/* Do the p6,14 combo - do p14 first so register assignments come out in same relative order as for p2,10 */\n\t"\
		"movq		%[__c14],%%rdx\n\t"\
	"/* Real parts: */\n\t"\
		"movaps		0xe0(%%rax)	,%%xmm6		/* a[j1+p14], this is the scratch xmm register  */\n\t"\
		"movaps		0xe0(%%rax)	,%%xmm4		/* a[j1+p14], this is the active  xmm register */\n\t"\
		"unpckhpd	0xe0(%%rbx)	,%%xmm6		/* a[j2+p14] gets read twice */\n\t"\
		"unpcklpd	0xe0(%%rbx)	,%%xmm4		/* a[jt+p14] */\n\t"\
		"movaps			%%xmm6	,0x160(%%rcx)	/* Store hi real in t15+16 */\n\t"\
	"/* Imag parts: */\n\t"\
		"movaps		0xf0(%%rax)	,%%xmm7\n\t"\
		"movaps		0xf0(%%rax)	,%%xmm5\n\t"\
		"unpckhpd	0xf0(%%rbx)	,%%xmm7\n\t"\
		"unpcklpd	0xf0(%%rbx)	,%%xmm5		/* a[jp+p14] */\n\t"\
		"movaps			%%xmm7	,0x170(%%rcx)	/* Store hi imag in t16+16 */\n\t"\
		"\n\t"\
		"movaps			%%xmm4	,%%xmm6		/* xmm6 <- cpy a[jt+p14] */\n\t"\
		"movaps			%%xmm5	,%%xmm7		/* xmm7 <- cpy a[jp+p14] */\n\t"\
		"\n\t"\
		"mulpd		    (%%rdx)	,%%xmm4		/* a[jt+p14]*c14 */\n\t"\
		"mulpd		    (%%rdx)	,%%xmm5		/* a[jp+p14]*c14 */\n\t"\
		"mulpd		0x10(%%rdx)	,%%xmm6		/* a[jt+p14]*s14 */\n\t"\
		"mulpd		0x10(%%rdx)	,%%xmm7		/* a[jp+p14]*s14 */\n\t"\
		"addpd		%%xmm6		,%%xmm5		/* xmm5 <- it */\n\t"\
		"subpd		%%xmm7		,%%xmm4		/* xmm4 <- rt		xmm6,7 free */\n\t"\
		"movaps		%%xmm5		,0x070(%%rcx)	/* Store it in t16*/\n\t"\
		"movaps		%%xmm4		,0x060(%%rcx)	/* Store rt in t15*/\n\t"\
		"\n\t"\
		"movq		%[__c6] ,%%rdx\n\t"\
	"/* Real parts: */\n\t"\
		"movaps		0x60(%%rax)	,%%xmm6		/* a[j1+p6 ], this is the scratch xmm register  */\n\t"\
		"movaps		0x60(%%rax)	,%%xmm4		/* a[j1+p6 ], this is the active  xmm register */\n\t"\
		"unpckhpd	0x60(%%rbx)	,%%xmm6		/* a[j2+p6 ] gets read twice */\n\t"\
		"unpcklpd	0x60(%%rbx)	,%%xmm4		/* a[jt+p6 ] */\n\t"\
		"movaps			%%xmm6	,0x140(%%rcx)	/* Store hi real in t13+16 */\n\t"\
	"/* Imag parts: */\n\t"\
		"movaps		0x70(%%rax)	,%%xmm7\n\t"\
		"movaps		0x70(%%rax)	,%%xmm5\n\t"\
		"unpckhpd	0x70(%%rbx)	,%%xmm7\n\t"\
		"unpcklpd	0x70(%%rbx)	,%%xmm5		/* a[jp+p6 ] */\n\t"\
		"movaps			%%xmm7	,0x150(%%rcx)	/* Store hi imag in t14+16 */\n\t"\
		"\n\t"\
		"movaps			%%xmm4	,%%xmm6	/* xmm6 <- cpy a[jt+p6] */\n\t"\
		"movaps			%%xmm5	,%%xmm7	/* xmm7 <- cpy a[jp+p6] */\n\t"\
		"\n\t"\
		"mulpd		    (%%rdx)	,%%xmm4		/* a[jt+p6]*c6 */\n\t"\
		"mulpd		    (%%rdx)	,%%xmm5		/* a[jp+p6]*c6 */\n\t"\
		"mulpd		0x10(%%rdx)	,%%xmm6		/* a[jt+p6]*s6 */\n\t"\
		"mulpd		0x10(%%rdx)	,%%xmm7		/* a[jp+p6]*s6 */\n\t"\
		"addpd		%%xmm6		,%%xmm5	/* xmm5 <- t14*/\n\t"\
		"subpd		%%xmm7		,%%xmm4	/* xmm4 <- t13*/\n\t"\
		"movaps		%%xmm5		,%%xmm7	/* xmm7 <- cpy t14*/\n\t"\
		"movaps		%%xmm4		,%%xmm6	/* xmm6 <- cpy t13*/\n\t"\
		"\n\t"\
		"subpd		0x060(%%rcx)	,%%xmm4		/* ~t15<- t13-rt */\n\t"\
		"subpd		0x070(%%rcx)	,%%xmm5		/* ~t16<- t14-it */\n\t"\
		"addpd		0x060(%%rcx)	,%%xmm6		/* ~t13<- t13+rt */\n\t"\
		"addpd		0x070(%%rcx)	,%%xmm7		/* ~t14<- t14+it */\n\t"\
		"\n\t"\
	"/* Finish radix-4 butterfly and store results into temporary-array slots: */\n\t"\
	"/*\n\t"\
	"~t13=t9 -t5;		~t9 =t9 +t5;\n\t"\
	"~t14=t10-t6;		~t10=t10+t6;\n\t"\
	"*/\n\t"\
		"subpd		%%xmm6		,%%xmm0	/*~t13*/\n\t"\
		"subpd		%%xmm7		,%%xmm1	/*~t14*/\n\t"\
		"movaps		%%xmm0		,0x040(%%rcx)	/* a[jt+p8 ] <- ~t13*/\n\t"\
		"movaps		%%xmm1		,0x050(%%rcx)	/* a[jp+p8 ] <- ~t14*/\n\t"\
		"addpd		%%xmm6		,%%xmm6	/* 2*t13*/\n\t"\
		"addpd		%%xmm7		,%%xmm7	/* 2*t14*/\n\t"\
		"addpd		%%xmm0		,%%xmm6	/*~t9 */\n\t"\
		"addpd		%%xmm1		,%%xmm7	/*~t10*/\n\t"\
		"movaps		%%xmm6		,     (%%rcx)	/* a[jt    ] <- ~t9 */\n\t"\
		"movaps		%%xmm7		,0x010(%%rcx)	/* a[jp    ] <- ~t10*/\n\t"\
		"\n\t"\
	"/*\n\t"\
	"~t15=t11+t8;		~t11=t11-t8;\n\t"\
	"~t16=t12-t7;		~t12=t12+t7;\n\t"\
	"*/\n\t"\
		"subpd		%%xmm5		,%%xmm2	/*~t11*/\n\t"\
		"subpd		%%xmm4		,%%xmm3	/*~t16*/\n\t"\
		"movaps		%%xmm2		,0x020(%%rcx)	/* a[jt+p4 ] <- ~t11*/\n\t"\
		"movaps		%%xmm3		,0x070(%%rcx)	/* a[jp+p12] <- ~t16*/\n\t"\
		"addpd		%%xmm5		,%%xmm5	/* 2*t16*/\n\t"\
		"addpd		%%xmm4		,%%xmm4	/* 2*t15*/\n\t"\
		"addpd		%%xmm2		,%%xmm5	/*~t15*/\n\t"\
		"addpd		%%xmm3		,%%xmm4	/*~t12*/\n\t"\
		"movaps		%%xmm5		,0x060(%%rcx)	/* a[jt+p12] <- ~t15*/\n\t"\
		"movaps		%%xmm4		,0x030(%%rcx)	/* a[jp+p4 ] <- ~t12*/\n\t"\
		"\n\t"\
"/******************************************************************************************************************************/\n\t"\
"/* Each of the above 2 radix-4 blocks is unique, but the ensuing 2 blocks                                                     */\n\t"\
"/* [operating on the odd-indexed elements from the unpck*pd commands which were stored to temporaries can use a common macro: */\n\t"\
"/******************************************************************************************************************************/\n\t"\
		"\n\t"\
"/*...Block 3: */\n\t"\
		"\n\t"\
"/*	SSE2_RADIX4_DIF_4TWIDDLE(r17,r21,r19,r23, r17,c1): */\n\t"\
	"/* Do the p0,p8 combo: */\n\t"\
		"movq		%[__r17],%%rax\n\t"\
		"movq		%[__c1] ,%%rbx\n\t"\
		"movq		%%rax   ,%%rcx\n\t"\
		"addq		$0x20   ,%%rcx	/* r19 */\n\t"\
		"\n\t"\
		"movaps		    (%%rax)	,%%xmm0		/* a[jt   ] */				\n\t	movaps	    (%%rcx),%%xmm4		/* a[jt+p8 ] */\n\t"\
		"movaps		0x10(%%rax)	,%%xmm1		/* a[jp   ] */				\n\t	movaps	0x10(%%rcx),%%xmm5		/* a[jp+p8 ] */\n\t"\
		"movaps		    (%%rbx)	,%%xmm6		/* c0 */\n\t"\
		"movaps		0x10(%%rbx)	,%%xmm7		/* s0 */\n\t"\
		"movaps		%%xmm0		,%%xmm2		/* xmm2 <- cpy a[jt   ] */\n\t"\
		"movaps		%%xmm1		,%%xmm3		/* xmm3 <- cpy a[jp   ] */\n\t"\
		"\n\t"\
		"mulpd   	%%xmm6		,%%xmm0		/* a[jt   ]*c0 */\n\t"\
		"mulpd   	%%xmm6		,%%xmm1		/* a[jp   ]*c0 */\n\t"\
		"mulpd   	%%xmm7		,%%xmm2		/* a[jt   ]*s0 */\n\t"\
		"mulpd   	%%xmm7		,%%xmm3		/* a[jp   ]*s0	xmm6,7 free */\n\t"\
		"																\n\t	movaps	%%xmm4		,%%xmm6		/* xmm6 <- cpy a[jt+p8 ] */\n\t"\
		"addpd   	%%xmm2		,%%xmm1		/* xmm1 <- t2 */			\n\t	movaps	%%xmm5		,%%xmm7		/* xmm7 <- cpy a[jp+p8 ] */\n\t"\
		"																\n\t	mulpd   0x20(%%rbx)	,%%xmm4		/* a[jt+p8 ]*c8 */\n\t"\
		"subpd   	%%xmm3		,%%xmm0		/* xmm0 <- t1 */			\n\t	mulpd   0x20(%%rbx)	,%%xmm5		/* a[jp+p8 ]*c8 */\n\t"\
		"																\n\t	mulpd   0x30(%%rbx)	,%%xmm6		/* a[jt+p8 ]*s8 */\n\t"\
		"movaps		%%xmm0		,%%xmm2		/* xmm2 <- cpy t1 */		\n\t	mulpd	0x30(%%rbx)	,%%xmm7		/* a[jp+p8 ]*s8 */\n\t"\
		"																\n\t	addpd   %%xmm6	    ,%%xmm5		/* xmm5 <- it */\n\t"\
		"movaps		%%xmm1		,%%xmm3		/* xmm3 <- cpy t2 */		\n\t	subpd	%%xmm7		,%%xmm4		/* xmm4 <- rt    xmm6,7 free */\n\t"\
		"\n\t"\
			"															\n\t	addq	$0x40 ,%%rcx	/* r23 */\n\t"\
		"																\n\t	addq	$0x60 ,%%rbx\n\t"\
		"																\n\t	movaps	    (%%rcx),%%xmm6		/* a[jt+p12] */\n\t"\
		"																\n\t	movaps	0x10(%%rcx),%%xmm7		/* a[jp+p12] */\n\t"\
		"\n\t"\
		"addpd		%%xmm4		,%%xmm0		/* ~t1 <- t1 +rt */\n\t"\
		"addpd		%%xmm5		,%%xmm1		/* ~t2 <- t2 +it */\n\t"\
		"subpd		%%xmm4		,%%xmm2		/* ~t3 <- t1 -rt */\n\t"\
		"subpd		%%xmm5		,%%xmm3		/* ~t4 <- t2 -it	xmm4,5 free */\n\t"\
		"\n\t"\
	"/* Do the p4,12 combo: */\n\t"\
		"movaps		%%xmm6		,%%xmm4		/* xmm4 <- cpy a[jt+p12] */\n\t"\
		"movaps		%%xmm7		,%%xmm5		/* xmm5 <- cpy a[jp+p12] */\n\t"\
		"\n\t"\
		"mulpd		    (%%rbx)	,%%xmm4		/* a[jt+p12]*c12 */\n\t"\
		"mulpd		    (%%rbx)	,%%xmm5		/* a[jp+p12]*c12 */\n\t"\
		"mulpd		0x10(%%rbx)	,%%xmm6		/* a[jt+p12]*s12 */\n\t"\
		"mulpd		0x10(%%rbx)	,%%xmm7		/* a[jp+p12]*s12 */\n\t"\
		"movq		%%rax,%%rdx	/* r17 */\n\t"\
		"addpd		%%xmm6		,%%xmm5	/* xmm5 <- it */\n\t"\
		"subpd		%%xmm7		,%%xmm4	/* xmm4 <- rt */\n\t"\
		"movaps		%%xmm5		,0x010(%%rdx)	/* store it */\n\t"\
		"movaps		%%xmm4		,     (%%rdx)	/* store rt */\n\t"\
		"\n\t"\
		"addq		$0x40 ,%%rax	/* r21 */\n\t"\
		"subq		$0x20 ,%%rbx\n\t"\
		"movaps		    (%%rax)	,%%xmm4		/* a[jt+p4] */\n\t"\
		"movaps		0x10(%%rax)	,%%xmm5		/* a[jp+p4] */\n\t"\
		"movaps			%%xmm4	,%%xmm6		/* xmm4 <- cpy a[jt+p4] */\n\t"\
		"movaps			%%xmm5	,%%xmm7		/* xmm5 <- cpy a[jp+p4] */\n\t"\
		"\n\t"\
		"mulpd		    (%%rbx)	,%%xmm4		/* a[jt+p4]*c4 */\n\t"\
		"mulpd		    (%%rbx)	,%%xmm5		/* a[jp+p4]*c4 */\n\t"\
		"mulpd		0x10(%%rbx)	,%%xmm6		/* a[jt+p4]*s4 */\n\t"\
		"mulpd		0x10(%%rbx)	,%%xmm7		/* a[jp+p4]*s4 */\n\t"\
		"addpd		%%xmm6		,%%xmm5	/* xmm5 <- t6 */\n\t"\
		"subpd		%%xmm7		,%%xmm4	/* xmm4 <- t5 	xmm6,7 free */\n\t"\
		"movaps		%%xmm5		,%%xmm7	/* xmm7 <- cpy t6 */\n\t"\
		"movaps		%%xmm4		,%%xmm6	/* xmm6 <- cpy t5 */\n\t"\
		"\n\t"\
		"subpd		     (%%rdx),%%xmm4		/* ~t7 <- t5 -rt */\n\t"\
		"subpd		0x010(%%rdx),%%xmm5		/* ~t8 <- t6 -it */\n\t"\
		"addpd		     (%%rdx),%%xmm6		/* ~t5 <- t5 +rt */\n\t"\
		"addpd		0x010(%%rdx),%%xmm7		/* ~t6 <- t6 +it */\n\t"\
		"\n\t"\
	"/* Finish radix-4 butterfly and store results into temporary-array slots: */\n\t"\
		"subpd		%%xmm6		,%%xmm0	/*~t5 */						\n\t	subpd	%%xmm5,%%xmm2			/*~t3 */\n\t"\
		"subpd		%%xmm7		,%%xmm1	/*~t6 */						\n\t	subpd	%%xmm4,%%xmm3			/*~t8 */\n\t"\
		"movaps		%%xmm0		,0x040(%%rdx)	/* a[jt+p8 ] <- ~t5 */	\n\t	movaps	%%xmm2,0x020(%%rdx)		/* a[jt+p4 ] <- ~t3 */\n\t"\
		"movaps		%%xmm1		,0x050(%%rdx)	/* a[jp+p8 ] <- ~t6 */	\n\t	movaps	%%xmm3,0x070(%%rdx)		/* a[jp+p12] <- ~t8 */\n\t"\
		"addpd		%%xmm6		,%%xmm6	/* 2*t5 */						\n\t	addpd	%%xmm5,%%xmm5			/* 2*t8 */\n\t"\
		"addpd		%%xmm7		,%%xmm7	/* 2*t6 */						\n\t	addpd	%%xmm4,%%xmm4			/* 2*t7 */\n\t"\
		"addpd		%%xmm0		,%%xmm6	/*~t1 */						\n\t	addpd	%%xmm2,%%xmm5			/*~t7 */\n\t"\
		"addpd		%%xmm1		,%%xmm7	/*~t2 */						\n\t	addpd	%%xmm3,%%xmm4			/*~t4 */\n\t"\
		"movaps		%%xmm6		,     (%%rdx)	/* a[jt    ] <- ~t1 */	\n\t	movaps	%%xmm5,0x060(%%rdx)		/* a[jt+p12] <- ~t7 */\n\t"\
		"movaps		%%xmm7		,0x010(%%rdx)	/* a[jp    ] <- ~t2 */	\n\t	movaps	%%xmm4,0x030(%%rdx)		/* a[jp+p4 ] <- ~t4 */\n\t"\
		"\n\t"\
"/*...Block 4: */\n\t"\
		"\n\t"\
"/* SSE2_RADIX4_DIF_4TWIDDLE(r25,r29,r27,r31, r25,c3): */\n\t"\
	"/* Do the p0,p8 combo: */\n\t"\
		"movq		%[__r25],%%rax\n\t"\
		"movq		%[__c3] ,%%rbx\n\t"\
		"movq		%%rax   ,%%rcx\n\t"\
		"addq		$0x20   ,%%rcx	/* r27 */\n\t"\
		"\n\t"\
		"movaps		    (%%rax)	,%%xmm0		/* a[jt   ] */				\n\t	movaps	    (%%rcx),%%xmm4		/* a[jt+p8 ] */\n\t"\
		"movaps		0x10(%%rax)	,%%xmm1		/* a[jp   ] */				\n\t	movaps	0x10(%%rcx),%%xmm5		/* a[jp+p8 ] */\n\t"\
		"movaps		    (%%rbx)	,%%xmm6		/* c0 */\n\t"\
		"movaps		0x10(%%rbx)	,%%xmm7		/* s0 */\n\t"\
		"movaps		%%xmm0		,%%xmm2		/* xmm2 <- cpy a[jt   ] */\n\t"\
		"movaps		%%xmm1		,%%xmm3		/* xmm3 <- cpy a[jp   ] */\n\t"\
		"\n\t"\
		"mulpd   	%%xmm6		,%%xmm0		/* a[jt   ]*c0 */\n\t"\
		"mulpd   	%%xmm6		,%%xmm1		/* a[jp   ]*c0 */\n\t"\
		"mulpd   	%%xmm7		,%%xmm2		/* a[jt   ]*s0 */\n\t"\
		"mulpd   	%%xmm7		,%%xmm3		/* a[jp   ]*s0	xmm6,7 free */\n\t"\
		"																\n\t	movaps	%%xmm4		,%%xmm6		/* xmm6 <- cpy a[jt+p8 ] */\n\t"\
		"addpd   	%%xmm2		,%%xmm1		/* xmm1 <- t2 */			\n\t	movaps	%%xmm5		,%%xmm7		/* xmm7 <- cpy a[jp+p8 ] */\n\t"\
		"																\n\t	mulpd   0x20(%%rbx)	,%%xmm4		/* a[jt+p8 ]*c8 */\n\t"\
		"subpd   	%%xmm3		,%%xmm0		/* xmm0 <- t1 */			\n\t	mulpd   0x20(%%rbx)	,%%xmm5		/* a[jp+p8 ]*c8 */\n\t"\
		"																\n\t	mulpd   0x30(%%rbx)	,%%xmm6		/* a[jt+p8 ]*s8 */\n\t"\
		"movaps		%%xmm0		,%%xmm2		/* xmm2 <- cpy t1 */		\n\t	mulpd	0x30(%%rbx)	,%%xmm7		/* a[jp+p8 ]*s8 */\n\t"\
		"																\n\t	addpd   %%xmm6	    ,%%xmm5		/* xmm5 <- it */\n\t"\
		"movaps		%%xmm1		,%%xmm3		/* xmm3 <- cpy t2 */		\n\t	subpd	%%xmm7		,%%xmm4		/* xmm4 <- rt    xmm6,7 free */\n\t"\
		"\n\t"\
			"															\n\t	addq	$0x40 ,%%rcx	/* r31 */\n\t"\
		"																\n\t	addq	$0x60 ,%%rbx\n\t"\
		"																\n\t	movaps	    (%%rcx),%%xmm6		/* a[jt+p12] */\n\t"\
		"																\n\t	movaps	0x10(%%rcx),%%xmm7		/* a[jp+p12] */\n\t"\
		"\n\t"\
		"addpd		%%xmm4		,%%xmm0		/* ~t1 <- t1 +rt */\n\t"\
		"addpd		%%xmm5		,%%xmm1		/* ~t2 <- t2 +it */\n\t"\
		"subpd		%%xmm4		,%%xmm2		/* ~t3 <- t1 -rt */\n\t"\
		"subpd		%%xmm5		,%%xmm3		/* ~t4 <- t2 -it	xmm4,5 free */\n\t"\
		"\n\t"\
	"/* Do the p4,12 combo: */\n\t"\
		"movaps		%%xmm6		,%%xmm4		/* xmm4 <- cpy a[jt+p12] */\n\t"\
		"movaps		%%xmm7		,%%xmm5		/* xmm5 <- cpy a[jp+p12] */\n\t"\
		"\n\t"\
		"mulpd		    (%%rbx)	,%%xmm4		/* a[jt+p12]*c12 */\n\t"\
		"mulpd		    (%%rbx)	,%%xmm5		/* a[jp+p12]*c12 */\n\t"\
		"mulpd		0x10(%%rbx)	,%%xmm6		/* a[jt+p12]*s12 */\n\t"\
		"mulpd		0x10(%%rbx)	,%%xmm7		/* a[jp+p12]*s12 */\n\t"\
		"movq		%%rax,%%rdx	/* r25 */\n\t"\
		"addpd		%%xmm6		,%%xmm5	/* xmm5 <- it */\n\t"\
		"subpd		%%xmm7		,%%xmm4	/* xmm4 <- rt */\n\t"\
		"movaps		%%xmm5		,0x010(%%rdx)	/* store it */\n\t"\
		"movaps		%%xmm4		,     (%%rdx)	/* store rt */\n\t"\
		"\n\t"\
		"addq		$0x40 ,%%rax	/* r29 */\n\t"\
		"subq		$0x20 ,%%rbx\n\t"\
		"movaps		    (%%rax)	,%%xmm4		/* a[jt+p4] */\n\t"\
		"movaps		0x10(%%rax)	,%%xmm5		/* a[jp+p4] */\n\t"\
		"movaps			%%xmm4	,%%xmm6		/* xmm4 <- cpy a[jt+p4] */\n\t"\
		"movaps			%%xmm5	,%%xmm7		/* xmm5 <- cpy a[jp+p4] */\n\t"\
		"\n\t"\
		"mulpd		    (%%rbx)	,%%xmm4		/* a[jt+p4]*c4 */\n\t"\
		"mulpd		    (%%rbx)	,%%xmm5		/* a[jp+p4]*c4 */\n\t"\
		"mulpd		0x10(%%rbx)	,%%xmm6		/* a[jt+p4]*s4 */\n\t"\
		"mulpd		0x10(%%rbx)	,%%xmm7		/* a[jp+p4]*s4 */\n\t"\
		"addpd		%%xmm6		,%%xmm5	/* xmm5 <- t6 */\n\t"\
		"subpd		%%xmm7		,%%xmm4	/* xmm4 <- t5 	xmm6,7 free */\n\t"\
		"movaps		%%xmm5		,%%xmm7	/* xmm7 <- cpy t6 */\n\t"\
		"movaps		%%xmm4		,%%xmm6	/* xmm6 <- cpy t5 */\n\t"\
		"\n\t"\
		"subpd		     (%%rdx),%%xmm4		/* ~t7 <- t5 -rt */\n\t"\
		"subpd		0x010(%%rdx),%%xmm5		/* ~t8 <- t6 -it */\n\t"\
		"addpd		     (%%rdx),%%xmm6		/* ~t5 <- t5 +rt */\n\t"\
		"addpd		0x010(%%rdx),%%xmm7		/* ~t6 <- t6 +it */\n\t"\
		"\n\t"\
	"/* Finish radix-4 butterfly and store results into temporary-array slots: */\n\t"\
		"subpd		%%xmm6		,%%xmm0	/*~t5 */						\n\t	subpd	%%xmm5,%%xmm2			/*~t3 */\n\t"\
		"subpd		%%xmm7		,%%xmm1	/*~t6 */						\n\t	subpd	%%xmm4,%%xmm3			/*~t8 */\n\t"\
		"movaps		%%xmm0		,0x040(%%rdx)	/* a[jt+p8 ] <- ~t5 */	\n\t	movaps	%%xmm2,0x020(%%rdx)		/* a[jt+p4 ] <- ~t3 */\n\t"\
		"movaps		%%xmm1		,0x050(%%rdx)	/* a[jp+p8 ] <- ~t6 */	\n\t	movaps	%%xmm3,0x070(%%rdx)		/* a[jp+p12] <- ~t8 */\n\t"\
		"addpd		%%xmm6		,%%xmm6	/* 2*t5 */						\n\t	addpd	%%xmm5,%%xmm5			/* 2*t8 */\n\t"\
		"addpd		%%xmm7		,%%xmm7	/* 2*t6 */						\n\t	addpd	%%xmm4,%%xmm4			/* 2*t7 */\n\t"\
		"addpd		%%xmm0		,%%xmm6	/*~t1 */						\n\t	addpd	%%xmm2,%%xmm5			/*~t7 */\n\t"\
		"addpd		%%xmm1		,%%xmm7	/*~t2 */						\n\t	addpd	%%xmm3,%%xmm4			/*~t4 */\n\t"\
		"movaps		%%xmm6		,     (%%rdx)	/* a[jt    ] <- ~t1 */	\n\t	movaps	%%xmm5,0x060(%%rdx)		/* a[jt+p12] <- ~t7 */\n\t"\
		"movaps		%%xmm7		,0x010(%%rdx)	/* a[jp    ] <- ~t2 */	\n\t	movaps	%%xmm4,0x030(%%rdx)		/* a[jp+p4 ] <- ~t4 */\n\t"\
		"\n\t"\
"/**************************************************************************************/\n\t"\
"/*...and now do four more radix-4 transforms, including the internal twiddle factors: */\n\t"\
"/**************************************************************************************/\n\t"\
		"\n\t"\
"/*...Block 1: t1,9,17,25 */\n\t"\
		"movq		%[__r1] ,%%rax\n\t"\
		"movq		%[__r9] ,%%rbx\n\t"\
		"movq		%[__r17],%%rcx\n\t"\
		"movq		%[__r25],%%rdx\n\t"\
		"\n\t"\
		"movaps		     (%%rax)	,%%xmm0		/* t1  */\n\t"\
		"movaps		0x010(%%rax)	,%%xmm1		/* t2  */\n\t"\
		"movaps		     (%%rbx)	,%%xmm2		/* t9  */\n\t"\
		"movaps		0x010(%%rbx)	,%%xmm3		/* t10 */\n\t"\
		"\n\t"\
		"subpd		     (%%rbx)	,%%xmm0		/* t9 =t1 -rt */\n\t"\
		"subpd		0x010(%%rbx)	,%%xmm1		/* t10=t2 -it */\n\t"\
		"addpd		     (%%rax)	,%%xmm2		/* t1 =t1 +rt */\n\t"\
		"addpd		0x010(%%rax)	,%%xmm3		/* t2 =t2 +it */\n\t"\
		"\n\t"\
		"movaps		     (%%rcx)	,%%xmm4		/* t17 */\n\t"\
		"movaps		0x010(%%rcx)	,%%xmm5		/* t18 */\n\t"\
		"movaps		     (%%rdx)	,%%xmm6		/* t25 */\n\t"\
		"movaps		0x010(%%rdx)	,%%xmm7		/* t26 */\n\t"\
		"\n\t"\
		"subpd		     (%%rdx)	,%%xmm4		/* t25=t17-rt */\n\t"\
		"subpd		0x010(%%rdx)	,%%xmm5		/* t26=t18-it */\n\t"\
		"addpd		     (%%rcx)	,%%xmm6		/* t17=t17+rt */\n\t"\
		"addpd		0x010(%%rcx)	,%%xmm7		/* t18=t18+it */\n\t"\
		"\n\t"\
		"subpd		%%xmm6		,%%xmm2		/* t1  <- t1 -t17 */\n\t"\
		"subpd		%%xmm7		,%%xmm3		/* t2  <- t2 -t18 */\n\t"\
		"addpd		%%xmm6		,%%xmm6		/*          2*t17 */\n\t"\
		"addpd		%%xmm7		,%%xmm7		/*          2*t18 */\n\t"\
		"movaps		%%xmm2		,    (%%rcx)/* a[jt+p1 ], store in t17 */\n\t"\
		"movaps		%%xmm3		,0x10(%%rcx)/* a[jp+p1 ], store in t18 */\n\t"\
		"addpd		%%xmm2		,%%xmm6		/* t17 <- t1 +t17 */\n\t"\
		"addpd		%%xmm3		,%%xmm7		/* t18 <- t2 +t18 */\n\t"\
		"movaps		%%xmm6		,    (%%rax)/* a[jt+p0 ], store in t0  */\n\t"\
		"movaps		%%xmm7		,0x10(%%rax)/* a[jp+p0 ], store in t1  */\n\t"\
		"\n\t"\
		"subpd		%%xmm5		,%%xmm0		/* t9  <- t9 -t26 */\n\t"\
		"subpd		%%xmm4		,%%xmm1		/* t10 <- t10-t25 */\n\t"\
		"addpd		%%xmm5		,%%xmm5		/*          2*t26 */\n\t"\
		"addpd		%%xmm4		,%%xmm4		/*          2*t25 */\n\t"\
		"movaps		%%xmm0		,    (%%rbx)/* a[jt+p2 ], store in t9  */\n\t"\
		"movaps		%%xmm1		,0x10(%%rdx)/* a[jp+p3 ], store in t26 */\n\t"\
		"addpd		%%xmm0		,%%xmm5		/* t26 <- t9 +t26 */\n\t"\
		"addpd		%%xmm1		,%%xmm4		/* t25 <- t10+t25 */\n\t"\
		"movaps		%%xmm5		,    (%%rdx)/* a[jt+p3 ], store in t25 */\n\t"\
		"movaps		%%xmm4		,0x10(%%rbx)/* a[jp+p2 ], store in t10 */\n\t"\
		"\n\t"\
"/*...Block 3: t5,13,21,29 */\n\t"\
		"addq		$0x40,%%rax	/* r5  */	\n\t"\
		"addq		$0x40,%%rbx	/* r13 */	\n\t"\
		"addq		$0x40,%%rcx	/* r21 */	\n\t"\
		"addq		$0x40,%%rdx	/* r29 */	\n\t"\
		"\n\t"\
		"movaps		     (%%rax)	,%%xmm0		/* t5  */\n\t"\
		"movaps		0x010(%%rax)	,%%xmm1		/* t6  */\n\t"\
		"movaps		0x080(%%rax)	,%%xmm2		/* t13 */\n\t"\
		"movaps		0x090(%%rax)	,%%xmm3		/* t14 */\n\t"\
		"\n\t"\
		"subpd		0x090(%%rax)	,%%xmm0		/* t5 =t5 -t14*/\n\t"\
		"subpd		0x080(%%rax)	,%%xmm1		/* t14=t6 -t13*/\n\t"\
		"addpd		0x010(%%rax)	,%%xmm2		/* t6 =t13+t6 */\n\t"\
		"addpd		     (%%rax)	,%%xmm3		/* t13=t14+t5 */\n\t"\
		"movq		%[__isrt2],%%rsi\n\t"\
		"\n\t"\
		"movaps		0x100(%%rax)	,%%xmm4		/* t21 */\n\t"\
		"movaps		0x110(%%rax)	,%%xmm5		/* t22 */\n\t"\
		"movaps		0x180(%%rax)	,%%xmm6		/* t29 */\n\t"\
		"movaps		0x190(%%rax)	,%%xmm7		/* t30 */\n\t"\
		"\n\t"\
		"subpd		0x110(%%rax)	,%%xmm4		/* t21-t22 */\n\t"\
		"addpd		0x100(%%rax)	,%%xmm5		/* t22+t21 */\n\t"\
		"mulpd		(%%rsi)		,%%xmm4	/* t21 = (t21-t22)*ISRT2 */\n\t"\
		"mulpd		(%%rsi)		,%%xmm5	/* t22 = (t22+t21)*ISRT2 */\n\t"\
		"\n\t"\
		"addpd		0x190(%%rax)	,%%xmm6		/* t29+t30 */\n\t"\
		"subpd		0x180(%%rax)	,%%xmm7		/* t30-t29 */\n\t"\
		"mulpd		(%%rsi)		,%%xmm6	/*  rt = (t29+t30)*ISRT2 */\n\t"\
		"mulpd		(%%rsi)		,%%xmm7	/*  it = (t30-t29)*ISRT2 */\n\t"\
		"\n\t"\
		"subpd		%%xmm6		,%%xmm4		/* t21=t21-rt */\n\t"\
		"subpd		%%xmm7		,%%xmm5		/* t22=t22-it */\n\t"\
		"addpd		%%xmm6		,%%xmm6		/*      2* rt */\n\t"\
		"addpd		%%xmm7		,%%xmm7		/*      2* it */\n\t"\
		"addpd		%%xmm4		,%%xmm6		/* t29=t21+rt */\n\t"\
		"addpd		%%xmm5		,%%xmm7		/* t30=t22+it */\n\t"\
		"\n\t"\
		"subpd		%%xmm4		,%%xmm0		/* t5 -t21 */\n\t"\
		"subpd		%%xmm5		,%%xmm2		/* t6 -t22 */\n\t"\
		"addpd		%%xmm4		,%%xmm4		/*   2*t21 */\n\t"\
		"addpd		%%xmm5		,%%xmm5		/*   2*t22 */\n\t"\
		"\n\t"\
		"movaps		%%xmm0		,    (%%rcx)/* a[jt+p1 ] */\n\t"\
		"movaps		%%xmm2		,0x10(%%rcx)/* a[jp+p1 ] */\n\t"\
		"addpd		%%xmm0		,%%xmm4		/* t5 +t21 */\n\t"\
		"addpd		%%xmm2		,%%xmm5		/* t6 +t22 */\n\t"\
		"movaps		%%xmm4		,    (%%rax)/* a[jt+p0 ] */\n\t"\
		"movaps		%%xmm5		,0x10(%%rax)/* a[jp+p0 ] */\n\t"\
		"\n\t"\
		"subpd		%%xmm7		,%%xmm3		/* t13-t30 */\n\t"\
		"subpd		%%xmm6		,%%xmm1		/* t14-t29 */\n\t"\
		"addpd		%%xmm7		,%%xmm7		/*   2*t30 */\n\t"\
		"addpd		%%xmm6		,%%xmm6		/*   2*t29 */\n\t"\
		"movaps		%%xmm3		,    (%%rbx)/* a[jt+p2 ] */\n\t"\
		"movaps		%%xmm1		,0x10(%%rdx)/* a[jp+p3 ] */\n\t"\
		"addpd		%%xmm3		,%%xmm7		/* t13+t30 */\n\t"\
		"addpd		%%xmm1		,%%xmm6		/* t14+t29 */\n\t"\
		"movaps		%%xmm7		,    (%%rdx)/* a[jt+p3 ] */\n\t"\
		"movaps		%%xmm6		,0x10(%%rbx)/* a[jp+p2 ] */\n\t"\
		"\n\t"\
"/*...Block 2: t3,11,19,27 */\n\t"\
		"subq		$0x20,%%rax	/* r3  */	\n\t"\
		"subq		$0x20,%%rbx	/* r11 */	\n\t"\
		"subq		$0x20,%%rcx	/* r19 */	\n\t"\
		"subq		$0x20,%%rdx	/* r27 */	\n\t"\
		"movq		%[__cc0],%%rdi\n\t"\
		"\n\t"\
		"movaps		0x100(%%rax),%%xmm4		/* t19 */			\n\t	movaps	0x180(%%rax),%%xmm6	/* t27 */\n\t"\
		"movaps		0x110(%%rax),%%xmm5		/* t20 */			\n\t	movaps	0x190(%%rax),%%xmm7	/* t28 */\n\t"\
		"movaps		0x100(%%rax),%%xmm0		/* copy t19 */		\n\t	movaps	0x180(%%rax),%%xmm2	/* copy t27 */\n\t"\
		"movaps		0x110(%%rax),%%xmm1		/* copy t20 */		\n\t	movaps	0x190(%%rax),%%xmm3	/* copy t28 */\n\t"\
		"\n\t"\
		"mulpd		    (%%rdi)	,%%xmm4		/* t19*c */			\n\t	mulpd	0x10(%%rdi)	,%%xmm6	/* t27*s */\n\t"\
		"mulpd		0x10(%%rdi)	,%%xmm1		/* t20*s */			\n\t	mulpd	    (%%rdi)	,%%xmm3	/* t28*c */\n\t"\
		"mulpd		    (%%rdi)	,%%xmm5		/* t20*c */			\n\t	mulpd	0x10(%%rdi)	,%%xmm7	/* t28*s */\n\t"\
		"mulpd		0x10(%%rdi)	,%%xmm0		/* t19*s */			\n\t	mulpd	    (%%rdi)	,%%xmm2	/* t27*c */\n\t"\
		"subpd		%%xmm1		,%%xmm4	/* ~t19 */				\n\t	subpd	%%xmm3,%%xmm6	/* rt */\n\t"\
		"addpd		%%xmm0		,%%xmm5	/* ~t20 */				\n\t	addpd	%%xmm2,%%xmm7	/* it */\n\t"\
		"\n\t"\
		"subpd		%%xmm6		,%%xmm4		/*~t27=t19-rt */\n\t"\
		"subpd		%%xmm7		,%%xmm5		/*~t28=t20-it */\n\t"\
		"addpd		%%xmm6		,%%xmm6		/*      2* rt */\n\t"\
		"addpd		%%xmm7		,%%xmm7		/*      2* it */\n\t"\
		"addpd		%%xmm4		,%%xmm6		/*~t19=t19+rt */\n\t"\
		"addpd		%%xmm5		,%%xmm7		/*~t20=t20+it */\n\t"\
		"\n\t"\
		"movaps		0x080(%%rax)	,%%xmm2		/* t11 */\n\t"\
		"movaps		0x090(%%rax)	,%%xmm3		/* t12 */\n\t"\
		"subpd		0x090(%%rax)	,%%xmm2		/* t11-t12 */\n\t"\
		"addpd		0x080(%%rax)	,%%xmm3		/* t12+t11 */\n\t"\
		"mulpd		(%%rsi)		,%%xmm2	/* rt = (t11-t12)*ISRT2 */\n\t"\
		"mulpd		(%%rsi)		,%%xmm3	/* it = (t12+t11)*ISRT2 */\n\t"\
		"\n\t"\
		"movaps		     (%%rax)	,%%xmm0		/* t3  */\n\t"\
		"movaps		0x010(%%rax)	,%%xmm1		/* t4  */\n\t"\
		"\n\t"\
		"subpd		%%xmm2,%%xmm0					/*~t11=t3 -rt */\n\t"\
		"subpd		%%xmm3,%%xmm1					/*~t12=t4 -it */\n\t"\
		"addpd		     (%%rax)	,%%xmm2		/*~t3 =rt +t3 */\n\t"\
		"addpd		0x010(%%rax)	,%%xmm3		/*~t4 =it +t4 */\n\t"\
		"\n\t"\
		"\n\t"\
		"subpd		%%xmm6		,%%xmm2		/* t3 -t19 */\n\t"\
		"subpd		%%xmm7		,%%xmm3		/* t4 -t20 */\n\t"\
		"addpd		%%xmm6		,%%xmm6		/*   2*t19 */\n\t"\
		"addpd		%%xmm7		,%%xmm7		/*   2*t20 */\n\t"\
		"movaps		%%xmm2		,    (%%rcx)/* a[jt+p1 ] */\n\t"\
		"movaps		%%xmm3		,0x10(%%rcx)/* a[jp+p1 ] */\n\t"\
		"addpd		%%xmm2		,%%xmm6		/* t3 +t19 */\n\t"\
		"addpd		%%xmm3		,%%xmm7		/* t4 +t20 */\n\t"\
		"movaps		%%xmm6		,    (%%rax)/* a[jt+p0 ] */\n\t"\
		"movaps		%%xmm7		,0x10(%%rax)/* a[jp+p0 ] */\n\t"\
		"\n\t"\
		"subpd		%%xmm5		,%%xmm0		/* t11-t28 */\n\t"\
		"subpd		%%xmm4		,%%xmm1		/* t12-t27 */\n\t"\
		"addpd		%%xmm5		,%%xmm5		/*          2*t28 */\n\t"\
		"addpd		%%xmm4		,%%xmm4		/*          2*t27 */\n\t"\
		"movaps		%%xmm0		,    (%%rbx)/* a[jt+p2 ] */\n\t"\
		"movaps		%%xmm1		,0x10(%%rdx)/* a[jp+p3 ] */\n\t"\
		"addpd		%%xmm0		,		%%xmm5/* t11+t28 */\n\t"\
		"addpd		%%xmm1		,		%%xmm4/* t12+t27 */\n\t"\
		"movaps		%%xmm5		,    (%%rdx)/* a[jt+p3 ] */\n\t"\
		"movaps		%%xmm4		,0x10(%%rbx)/* a[jp+p2 ] */\n\t"\
		"\n\t"\
"/*...Block 4: t7,15,23,31 */\n\t"\
		"addq		$0x40,%%rax	/* r7  */	\n\t"\
		"addq		$0x40,%%rbx	/* r15 */	\n\t"\
		"addq		$0x40,%%rcx	/* r23 */	\n\t"\
		"addq		$0x40,%%rdx	/* r31 */	\n\t"\
		"\n\t"\
		"movaps		0x100(%%rax),%%xmm4		/* t23 */			\n\t	movaps	0x180(%%rax),%%xmm6		/* t31 */\n\t"\
		"movaps		0x110(%%rax),%%xmm5		/* t24 */			\n\t	movaps	0x190(%%rax),%%xmm7		/* t32 */\n\t"\
		"movaps		0x100(%%rax),%%xmm0		/* copy t23 */		\n\t	movaps	0x180(%%rax),%%xmm2		/* copy t31 */\n\t"\
		"movaps		0x110(%%rax),%%xmm1		/* copy t24 */		\n\t	movaps	0x190(%%rax),%%xmm3		/* copy t32 */\n\t"\
		"\n\t"\
		"mulpd		0x10(%%rdi)	,%%xmm4		/* t23*s */			\n\t	mulpd	    (%%rdi)	,%%xmm6		/* t31*c */\n\t"\
		"mulpd		    (%%rdi)	,%%xmm1		/* t24*c */			\n\t	mulpd	0x10(%%rdi)	,%%xmm3		/* t32*s */\n\t"\
		"mulpd		0x10(%%rdi)	,%%xmm5		/* t24*s */			\n\t	mulpd	    (%%rdi)	,%%xmm7		/* t32*c */\n\t"\
		"mulpd		    (%%rdi)	,%%xmm0		/* t23*c */			\n\t	mulpd	0x10(%%rdi)	,%%xmm2		/* t31*s */\n\t"\
		"subpd		%%xmm1		,%%xmm4	/* ~t23 */				\n\t	subpd	%%xmm3		,%%xmm6		/* rt */\n\t"\
		"addpd		%%xmm0		,%%xmm5	/* ~t24 */				\n\t	addpd	%%xmm2		,%%xmm7		/* it */\n\t"\
		"\n\t"\
		"subpd		%%xmm6		,%%xmm4		/*~t23=t23-rt */\n\t"\
		"subpd		%%xmm7		,%%xmm5		/*~t24=t24-it */\n\t"\
		"addpd		%%xmm6		,%%xmm6		/*      2* rt */\n\t"\
		"addpd		%%xmm7		,%%xmm7		/*      2* it */\n\t"\
		"addpd		%%xmm4		,%%xmm6		/*~t31=t23+rt */\n\t"\
		"addpd		%%xmm5		,%%xmm7		/*~t32=t24+it */\n\t"\
		"\n\t"\
		"movaps		0x080(%%rax)	,%%xmm2		/* t15 */\n\t"\
		"movaps		0x090(%%rax)	,%%xmm3		/* t16 */\n\t"\
		"addpd		0x090(%%rax)	,%%xmm2		/* t15+t16 */\n\t"\
		"subpd		0x080(%%rax)	,%%xmm3		/* t16-t15 */\n\t"\
		"mulpd		(%%rsi)		,%%xmm2	/* rt = (t15+t16)*ISRT2 */\n\t"\
		"mulpd		(%%rsi)		,%%xmm3	/* it = (t16-t15)*ISRT2 */\n\t"\
		"\n\t"\
		"movaps		     (%%rax)	,%%xmm0		/* t7  */\n\t"\
		"movaps		0x010(%%rax)	,%%xmm1		/* t8  */\n\t"\
		"\n\t"\
		"subpd		%%xmm2,%%xmm0					/*~t7 =t7 -rt */\n\t"\
		"subpd		%%xmm3,%%xmm1					/*~t8 =t8 -it */\n\t"\
		"addpd		     (%%rax)	,%%xmm2		/*~t15=rt +t7 */\n\t"\
		"addpd		0x010(%%rax)	,%%xmm3		/*~t16=it +t8 */\n\t"\
		"\n\t"\
		"subpd		%%xmm4		,%%xmm0		/* t7 -t23 */\n\t"\
		"subpd		%%xmm5		,%%xmm1		/* t8 -t24 */\n\t"\
		"addpd		%%xmm4		,%%xmm4		/*   2*t23 */\n\t"\
		"addpd		%%xmm5		,%%xmm5		/*   2*t24 */\n\t"\
		"movaps		%%xmm0		,    (%%rcx)/* a[jt+p1 ] */\n\t"\
		"movaps		%%xmm1		,0x10(%%rcx)/* a[jp+p1 ] */\n\t"\
		"addpd		%%xmm0		,%%xmm4		/* t7 +t23 */\n\t"\
		"addpd		%%xmm1		,%%xmm5		/* t8 +t24 */\n\t"\
		"movaps		%%xmm4		,    (%%rax)/* a[jt+p0 ] */\n\t"\
		"movaps		%%xmm5		,0x10(%%rax)/* a[jp+p0 ] */\n\t"\
		"\n\t"\
		"subpd		%%xmm7		,%%xmm2		/* t15-t32 */\n\t"\
		"subpd		%%xmm6		,%%xmm3		/* t16-t31 */\n\t"\
		"addpd		%%xmm7		,%%xmm7		/*   2*t32 */\n\t"\
		"addpd		%%xmm6		,%%xmm6		/*   2*t31 */\n\t"\
		"movaps		%%xmm2		,    (%%rbx)/* a[jt+p2 ] */\n\t"\
		"movaps		%%xmm3		,0x10(%%rdx)/* a[jp+p3 ] */\n\t"\
		"addpd		%%xmm2		,%%xmm7		/* t15+t32 */\n\t"\
		"addpd		%%xmm3		,%%xmm6		/* t16+t31 */\n\t"\
		"movaps		%%xmm7		,    (%%rdx)/* a[jt+p3 ] */\n\t"\
		"movaps		%%xmm6		,0x10(%%rbx)/* a[jp+p2 ] */\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__add1] "m" (Xadd1)\
		 ,[__r1] "m" (Xr1)\
		 ,[__r9] "m" (Xr9)\
		 ,[__r17] "m" (Xr17)\
		 ,[__r25] "m" (Xr25)\
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__cc0] "m" (Xcc0)\
		 ,[__c1] "m" (Xc1)\
		 ,[__c2] "m" (Xc2)\
		 ,[__c3] "m" (Xc3)\
		 ,[__c4] "m" (Xc4)\
		 ,[__c5] "m" (Xc5)\
		 ,[__c6] "m" (Xc6)\
		 ,[__c7] "m" (Xc7)\
		 ,[__c8] "m" (Xc8)\
		 ,[__c9] "m" (Xc9)\
		 ,[__c10] "m" (Xc10)\
		 ,[__c11] "m" (Xc11)\
		 ,[__c12] "m" (Xc12)\
		 ,[__c13] "m" (Xc13)\
		 ,[__c14] "m" (Xc14)\
		 ,[__c15] "m" (Xc15)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
	);\
	}

	#define SSE2_RADIX16_WRAPPER_DIT(Xadd0,Xadd1,Xr1,Xr9,Xr17,Xr25,Xisrt2,Xcc0,Xc1,Xc2,Xc3,Xc4,Xc5,Xc6,Xc7,Xc8,Xc9,Xc10,Xc11,Xc12,Xc13,Xc14,Xc15)\
	{\
	__asm__ volatile (\
		"/*...Block 1: */						\n\t"\
		"movq	%[__r1],%%rax					\n\t"\
		"movq	%%rax,%%rbx						\n\t"\
		"movq	%%rax,%%rcx						\n\t"\
		"movq	%%rax,%%rdx						\n\t"\
		"addq	$0x100,%%rbx					\n\t"\
		"addq	$0x080,%%rcx					\n\t"\
		"addq	$0x180,%%rdx					\n\t"\
		"/* SSE2_RADIX4_DIT_IN_PLACE(): */		\n\t"\
		"movaps	    (%%rax),%%xmm0				\n\t"\
		"movaps	0x10(%%rax),%%xmm1				\n\t"\
		"movaps	    (%%rax),%%xmm2				\n\t"\
		"movaps	0x10(%%rax),%%xmm3				\n\t"\
		"addpd	    (%%rbx),%%xmm0				\n\t"\
		"addpd	0x10(%%rbx),%%xmm1				\n\t"\
		"subpd	    (%%rbx),%%xmm2				\n\t"\
		"subpd	0x10(%%rbx),%%xmm3				\n\t"\
		"movaps	    (%%rcx),%%xmm4				\n\t"\
		"movaps	0x10(%%rcx),%%xmm5				\n\t"\
		"movaps	    (%%rcx),%%xmm6				\n\t"\
		"movaps	0x10(%%rcx),%%xmm7				\n\t"\
		"addpd	    (%%rdx),%%xmm4				\n\t"\
		"addpd	0x10(%%rdx),%%xmm5				\n\t"\
		"subpd	    (%%rdx),%%xmm6				\n\t"\
		"subpd	0x10(%%rdx),%%xmm7				\n\t"\
		"subpd	%%xmm4,%%xmm0					\n\t"\
		"subpd	%%xmm5,%%xmm1					\n\t"\
		"movaps	%%xmm0,     (%%rbx)				\n\t"\
		"movaps	%%xmm1,0x010(%%rbx)				\n\t"\
		"addpd	%%xmm4,%%xmm4					\n\t"\
		"addpd	%%xmm5,%%xmm5					\n\t"\
		"addpd	%%xmm0,%%xmm4					\n\t"\
		"addpd	%%xmm1,%%xmm5					\n\t"\
		"movaps	%%xmm4,     (%%rax)				\n\t"\
		"movaps	%%xmm5,0x010(%%rax)				\n\t"\
		"subpd	%%xmm7,%%xmm2					\n\t"\
		"subpd	%%xmm6,%%xmm3					\n\t"\
		"movaps	%%xmm2,     (%%rdx)				\n\t"\
		"movaps	%%xmm3,0x010(%%rcx)				\n\t"\
		"addpd	%%xmm7,%%xmm7					\n\t"\
		"addpd	%%xmm6,%%xmm6					\n\t"\
		"addpd	%%xmm2,%%xmm7					\n\t"\
		"addpd	%%xmm3,%%xmm6					\n\t"\
		"movaps	%%xmm7,     (%%rcx)				\n\t"\
		"movaps	%%xmm6,0x010(%%rdx)				\n\t"\
		"/*...Block 2: */						\n\t"\
		"addq	$0x040,%%rax					\n\t"\
		"addq	$0x040,%%rbx					\n\t"\
		"addq	$0x040,%%rcx					\n\t"\
		"addq	$0x040,%%rdx					\n\t"\
		"/* SSE2_RADIX4_DIT_IN_PLACE(): */		\n\t"\
		"movaps	    (%%rax),%%xmm0				\n\t"\
		"movaps	0x10(%%rax),%%xmm1				\n\t"\
		"movaps	    (%%rax),%%xmm2				\n\t"\
		"movaps	0x10(%%rax),%%xmm3				\n\t"\
		"addpd	    (%%rbx),%%xmm0				\n\t"\
		"addpd	0x10(%%rbx),%%xmm1				\n\t"\
		"subpd	    (%%rbx),%%xmm2				\n\t"\
		"subpd	0x10(%%rbx),%%xmm3				\n\t"\
		"movaps	    (%%rcx),%%xmm4				\n\t"\
		"movaps	0x10(%%rcx),%%xmm5				\n\t"\
		"movaps	    (%%rcx),%%xmm6				\n\t"\
		"movaps	0x10(%%rcx),%%xmm7				\n\t"\
		"addpd	    (%%rdx),%%xmm4				\n\t"\
		"addpd	0x10(%%rdx),%%xmm5				\n\t"\
		"subpd	    (%%rdx),%%xmm6				\n\t"\
		"subpd	0x10(%%rdx),%%xmm7				\n\t"\
		"subpd	%%xmm4,%%xmm0					\n\t"\
		"subpd	%%xmm5,%%xmm1					\n\t"\
		"movaps	%%xmm0,     (%%rbx)				\n\t"\
		"movaps	%%xmm1,0x010(%%rbx)				\n\t"\
		"addpd	%%xmm4,%%xmm4					\n\t"\
		"addpd	%%xmm5,%%xmm5					\n\t"\
		"addpd	%%xmm0,%%xmm4					\n\t"\
		"addpd	%%xmm1,%%xmm5					\n\t"\
		"movaps	%%xmm4,     (%%rax)				\n\t"\
		"movaps	%%xmm5,0x010(%%rax)				\n\t"\
		"subpd	%%xmm7,%%xmm2					\n\t"\
		"subpd	%%xmm6,%%xmm3					\n\t"\
		"movaps	%%xmm2,     (%%rdx)				\n\t"\
		"movaps	%%xmm3,0x010(%%rcx)				\n\t"\
		"addpd	%%xmm7,%%xmm7					\n\t"\
		"addpd	%%xmm6,%%xmm6					\n\t"\
		"addpd	%%xmm2,%%xmm7					\n\t"\
		"addpd	%%xmm3,%%xmm6					\n\t"\
		"movaps	%%xmm7,     (%%rcx)				\n\t"\
		"movaps	%%xmm6,0x010(%%rdx)				\n\t"\
		"/*...Block 3: */						\n\t"\
		"subq	$0x020,%%rax					\n\t"\
		"subq	$0x020,%%rbx					\n\t"\
		"subq	$0x020,%%rcx					\n\t"\
		"subq	$0x020,%%rdx					\n\t"\
		"/* SSE2_RADIX4_DIT_IN_PLACE(): */		\n\t"\
		"movaps	    (%%rax),%%xmm0				\n\t"\
		"movaps	0x10(%%rax),%%xmm1				\n\t"\
		"movaps	    (%%rax),%%xmm2				\n\t"\
		"movaps	0x10(%%rax),%%xmm3				\n\t"\
		"addpd	    (%%rbx),%%xmm0				\n\t"\
		"addpd	0x10(%%rbx),%%xmm1				\n\t"\
		"subpd	    (%%rbx),%%xmm2				\n\t"\
		"subpd	0x10(%%rbx),%%xmm3				\n\t"\
		"movaps	    (%%rcx),%%xmm4				\n\t"\
		"movaps	0x10(%%rcx),%%xmm5				\n\t"\
		"movaps	    (%%rcx),%%xmm6				\n\t"\
		"movaps	0x10(%%rcx),%%xmm7				\n\t"\
		"addpd	    (%%rdx),%%xmm4				\n\t"\
		"addpd	0x10(%%rdx),%%xmm5				\n\t"\
		"subpd	    (%%rdx),%%xmm6				\n\t"\
		"subpd	0x10(%%rdx),%%xmm7				\n\t"\
		"subpd	%%xmm4,%%xmm0					\n\t"\
		"subpd	%%xmm5,%%xmm1					\n\t"\
		"movaps	%%xmm0,     (%%rbx)				\n\t"\
		"movaps	%%xmm1,0x010(%%rbx)				\n\t"\
		"addpd	%%xmm4,%%xmm4					\n\t"\
		"addpd	%%xmm5,%%xmm5					\n\t"\
		"addpd	%%xmm0,%%xmm4					\n\t"\
		"addpd	%%xmm1,%%xmm5					\n\t"\
		"movaps	%%xmm4,     (%%rax)				\n\t"\
		"movaps	%%xmm5,0x010(%%rax)				\n\t"\
		"subpd	%%xmm7,%%xmm2					\n\t"\
		"subpd	%%xmm6,%%xmm3					\n\t"\
		"movaps	%%xmm2,     (%%rdx)				\n\t"\
		"movaps	%%xmm3,0x010(%%rcx)				\n\t"\
		"addpd	%%xmm7,%%xmm7					\n\t"\
		"addpd	%%xmm6,%%xmm6					\n\t"\
		"addpd	%%xmm2,%%xmm7					\n\t"\
		"addpd	%%xmm3,%%xmm6					\n\t"\
		"movaps	%%xmm7,     (%%rcx)				\n\t"\
		"movaps	%%xmm6,0x010(%%rdx)				\n\t"\
		"/*...Block 4: */						\n\t"\
		"addq	$0x040,%%rax					\n\t"\
		"addq	$0x040,%%rbx					\n\t"\
		"addq	$0x040,%%rcx					\n\t"\
		"addq	$0x040,%%rdx					\n\t"\
		"/* SSE2_RADIX4_DIT_IN_PLACE(): */		\n\t"\
		"movaps	    (%%rax),%%xmm0				\n\t"\
		"movaps	0x10(%%rax),%%xmm1				\n\t"\
		"movaps	    (%%rax),%%xmm2				\n\t"\
		"movaps	0x10(%%rax),%%xmm3				\n\t"\
		"addpd	    (%%rbx),%%xmm0				\n\t"\
		"addpd	0x10(%%rbx),%%xmm1				\n\t"\
		"subpd	    (%%rbx),%%xmm2				\n\t"\
		"subpd	0x10(%%rbx),%%xmm3				\n\t"\
		"movaps	    (%%rcx),%%xmm4				\n\t"\
		"movaps	0x10(%%rcx),%%xmm5				\n\t"\
		"movaps	    (%%rcx),%%xmm6				\n\t"\
		"movaps	0x10(%%rcx),%%xmm7				\n\t"\
		"addpd	    (%%rdx),%%xmm4				\n\t"\
		"addpd	0x10(%%rdx),%%xmm5				\n\t"\
		"subpd	    (%%rdx),%%xmm6				\n\t"\
		"subpd	0x10(%%rdx),%%xmm7				\n\t"\
		"subpd	%%xmm4,%%xmm0					\n\t"\
		"subpd	%%xmm5,%%xmm1					\n\t"\
		"movaps	%%xmm0,     (%%rbx)				\n\t"\
		"movaps	%%xmm1,0x010(%%rbx)				\n\t"\
		"addpd	%%xmm4,%%xmm4					\n\t"\
		"addpd	%%xmm5,%%xmm5					\n\t"\
		"addpd	%%xmm0,%%xmm4					\n\t"\
		"addpd	%%xmm1,%%xmm5					\n\t"\
		"movaps	%%xmm4,     (%%rax)				\n\t"\
		"movaps	%%xmm5,0x010(%%rax)				\n\t"\
		"subpd	%%xmm7,%%xmm2					\n\t"\
		"subpd	%%xmm6,%%xmm3					\n\t"\
		"movaps	%%xmm2,     (%%rdx)				\n\t"\
		"movaps	%%xmm3,0x010(%%rcx)				\n\t"\
		"addpd	%%xmm7,%%xmm7					\n\t"\
		"addpd	%%xmm6,%%xmm6					\n\t"\
		"addpd	%%xmm2,%%xmm7					\n\t"\
		"addpd	%%xmm3,%%xmm6					\n\t"\
		"movaps	%%xmm7,     (%%rcx)				\n\t"\
		"movaps	%%xmm6,0x010(%%rdx)				\n\t"\
		"/***************************************************************************************************/\n\t"\
		"/*..and now do four more radix-4 transforms, including the internal and external twiddle factors.  */\n\t"\
		"/*  Write even-index 16-byte output pairs to a(j1), odd-index to a(j2), unpack same as on inputs.  */\n\t"\
		"/*  We do the last 2 radix-4 blocks first, to make the unpack-interleaving of outputs smoother.    */\n\t"\
		"/***************************************************************************************************/\n\t"\
		"/* Main-array addresses still in add0,1, no need to re-init: */									  \n\t"\
		"/*...Block 3: t3,11,19,27 -> r9,13,11,15: */	\n\t"\
		"movq		%[__r9],%%rax				\n\t"\
		"movq		%[__isrt2],%%rbx			\n\t"\
		"movq		%[__cc0],%%rcx				\n\t"\
		"movaps		0x020(%%rax),%%xmm4			\n\t"\
		"movaps		0x060(%%rax),%%xmm0			\n\t"\
		"movaps		0x030(%%rax),%%xmm5			\n\t"\
		"movaps		0x070(%%rax),%%xmm1			\n\t"\
		"movaps		0x020(%%rax),%%xmm6			\n\t"\
		"movaps		0x060(%%rax),%%xmm2			\n\t"\
		"movaps		0x030(%%rax),%%xmm7			\n\t"\
		"movaps		0x070(%%rax),%%xmm3			\n\t"\
		"mulpd		    (%%rcx),%%xmm4			\n\t"\
		"mulpd		0x10(%%rcx),%%xmm0			\n\t"\
		"mulpd		    (%%rcx),%%xmm5			\n\t"\
		"mulpd		0x10(%%rcx),%%xmm1			\n\t"\
		"mulpd		0x10(%%rcx),%%xmm6			\n\t"\
		"mulpd		    (%%rcx),%%xmm2			\n\t"\
		"mulpd		0x10(%%rcx),%%xmm7			\n\t"\
		"mulpd		    (%%rcx),%%xmm3			\n\t"\
		"subpd		%%xmm6,%%xmm5				\n\t"\
		"subpd		%%xmm2,%%xmm1				\n\t"\
		"addpd		%%xmm7,%%xmm4				\n\t"\
		"addpd		%%xmm3,%%xmm0				\n\t"\
		"movaps		%%xmm5,%%xmm7				\n\t"\
		"movaps		%%xmm4,%%xmm6				\n\t"\
		"addpd		%%xmm0,%%xmm4				\n\t"\
		"addpd		%%xmm1,%%xmm5				\n\t"\
		"subpd		%%xmm0,%%xmm6				\n\t"\
		"subpd		%%xmm1,%%xmm7				\n\t"\
		"movaps		0x040(%%rax),%%xmm2			\n\t"\
		"movaps		0x050(%%rax),%%xmm3			\n\t"\
		"movaps		     (%%rax),%%xmm0			\n\t"\
		"movaps		0x010(%%rax),%%xmm1			\n\t"\
		"addpd		0x050(%%rax),%%xmm2			\n\t"\
		"subpd		0x040(%%rax),%%xmm3			\n\t"\
		"mulpd		(%%rbx),%%xmm2				\n\t"\
		"mulpd		(%%rbx),%%xmm3				\n\t"\
		"subpd		%%xmm2,%%xmm0				\n\t"\
		"subpd		%%xmm3,%%xmm1				\n\t"\
		"addpd		%%xmm2,%%xmm2				\n\t"\
		"addpd		%%xmm3,%%xmm3				\n\t"\
		"addpd		%%xmm0,%%xmm2				\n\t"\
		"addpd		%%xmm1,%%xmm3				\n\t"\
		"movq		%[__add1],%%rbx				\n\t"\
		"movq		%[__c1],%%rcx				\n\t"\
		"movq		%[__c9],%%rdx				\n\t"\
		"subpd		%%xmm4,%%xmm2				\n\t"\
		"subpd		%%xmm5,%%xmm3				\n\t"\
		"addpd		%%xmm4,%%xmm4				\n\t"\
		"addpd		%%xmm5,%%xmm5				\n\t"\
		"addpd		%%xmm2,%%xmm4				\n\t"\
		"addpd		%%xmm3,%%xmm5				\n\t"\
		"movaps		%%xmm2,     (%%rax)			\n\t"\
		"movaps		%%xmm3,0x010(%%rax)			\n\t"\
		"movaps		%%xmm4,%%xmm2				\n\t"\
		"movaps		%%xmm5,%%xmm3				\n\t"\
		"mulpd		    (%%rcx),%%xmm4			\n\t"\
		"mulpd		    (%%rcx),%%xmm5			\n\t"\
		"mulpd		0x10(%%rcx),%%xmm2			\n\t"\
		"mulpd		0x10(%%rcx),%%xmm3			\n\t"\
		"subpd		%%xmm2,%%xmm5				\n\t"\
		"addpd		%%xmm3,%%xmm4				\n\t"\
		"movaps		%%xmm5,0x10(%%rbx)			\n\t"\
		"movaps		%%xmm4,    (%%rbx)			\n\t"\
		"movaps		     (%%rax),%%xmm4			\n\t"\
		"movaps		0x010(%%rax),%%xmm5			\n\t"\
		"movaps		%%xmm4,%%xmm2				\n\t"\
		"movaps		%%xmm5,%%xmm3				\n\t"\
		"mulpd		    (%%rdx),%%xmm4			\n\t"\
		"mulpd		    (%%rdx),%%xmm5			\n\t"\
		"mulpd		0x10(%%rdx),%%xmm2			\n\t"\
		"mulpd		0x10(%%rdx),%%xmm3			\n\t"\
		"subpd		%%xmm2,%%xmm5				\n\t"\
		"addpd		%%xmm3,%%xmm4				\n\t"\
		"movaps		%%xmm5,0x90(%%rbx)			\n\t"\
		"movaps		%%xmm4,0x80(%%rbx)			\n\t"\
		"movq		%[__c5],%%rcx				\n\t"\
		"movq		%[__c13],%%rdx				\n\t"\
		"subpd		%%xmm7,%%xmm0				\n\t"\
		"subpd		%%xmm6,%%xmm1				\n\t"\
		"addpd		%%xmm7,%%xmm7				\n\t"\
		"addpd		%%xmm6,%%xmm6				\n\t"\
		"addpd		%%xmm0,%%xmm7				\n\t"\
		"addpd		%%xmm1,%%xmm6				\n\t"\
		"movaps		%%xmm7,%%xmm4				\n\t"\
		"movaps		%%xmm1,%%xmm5				\n\t"\
		"mulpd		    (%%rcx),%%xmm7			\n\t"\
		"mulpd		    (%%rcx),%%xmm1			\n\t"\
		"mulpd		0x10(%%rcx),%%xmm4			\n\t"\
		"mulpd		0x10(%%rcx),%%xmm5			\n\t"\
		"subpd		%%xmm4,%%xmm1				\n\t"\
		"addpd		%%xmm5,%%xmm7				\n\t"\
		"movaps		%%xmm1,0x50(%%rbx)			\n\t"\
		"movaps		%%xmm7,0x40(%%rbx)			\n\t"\
		"movaps		%%xmm0,%%xmm4				\n\t"\
		"movaps		%%xmm6,%%xmm5				\n\t"\
		"mulpd		    (%%rdx),%%xmm0			\n\t"\
		"mulpd		    (%%rdx),%%xmm6			\n\t"\
		"mulpd		0x10(%%rdx),%%xmm4			\n\t"\
		"mulpd		0x10(%%rdx),%%xmm5			\n\t"\
		"subpd		%%xmm4,%%xmm6				\n\t"\
		"addpd		%%xmm5,%%xmm0				\n\t"\
		"movaps		%%xmm6,0xd0(%%rbx)			\n\t"\
		"movaps		%%xmm0,0xc0(%%rbx)			\n\t"\
		"/*...Block 4: t7,15,23,31 -> r25,29,27,31: */	\n\t"\
		"movq		%[__r25],%%rax				\n\t"\
		"movq		%[__isrt2],%%rbx			\n\t"\
		"movq		%[__cc0],%%rcx				\n\t"\
		"movaps		0x020(%%rax),%%xmm4			\n\t"\
		"movaps		0x060(%%rax),%%xmm0			\n\t"\
		"movaps		0x030(%%rax),%%xmm5			\n\t"\
		"movaps		0x070(%%rax),%%xmm1			\n\t"\
		"movaps		0x020(%%rax),%%xmm6			\n\t"\
		"movaps		0x060(%%rax),%%xmm2			\n\t"\
		"movaps		0x030(%%rax),%%xmm7			\n\t"\
		"movaps		0x070(%%rax),%%xmm3			\n\t"\
		"mulpd		0x10(%%rcx),%%xmm4			\n\t"\
		"mulpd		    (%%rcx),%%xmm0			\n\t"\
		"mulpd		0x10(%%rcx),%%xmm5			\n\t"\
		"mulpd		    (%%rcx),%%xmm1			\n\t"\
		"mulpd		    (%%rcx),%%xmm6			\n\t"\
		"mulpd		0x10(%%rcx),%%xmm2			\n\t"\
		"mulpd		    (%%rcx),%%xmm7			\n\t"\
		"mulpd		0x10(%%rcx),%%xmm3			\n\t"\
		"subpd		%%xmm6,%%xmm5				\n\t"\
		"subpd		%%xmm2,%%xmm1				\n\t"\
		"addpd		%%xmm7,%%xmm4				\n\t"\
		"addpd		%%xmm3,%%xmm0				\n\t"\
		"movaps		%%xmm5,%%xmm7				\n\t"\
		"movaps		%%xmm4,%%xmm6				\n\t"\
		"addpd		%%xmm0,%%xmm4				\n\t"\
		"addpd		%%xmm1,%%xmm5				\n\t"\
		"subpd		%%xmm0,%%xmm6				\n\t"\
		"subpd		%%xmm1,%%xmm7				\n\t"\
		"movaps		0x040(%%rax),%%xmm2			\n\t"\
		"movaps		0x050(%%rax),%%xmm3			\n\t"\
		"movaps		     (%%rax),%%xmm0			\n\t"\
		"movaps		0x010(%%rax),%%xmm1			\n\t"\
		"subpd		0x050(%%rax),%%xmm2			\n\t"\
		"addpd		0x040(%%rax),%%xmm3			\n\t"\
		"mulpd		(%%rbx),%%xmm2				\n\t"\
		"mulpd		(%%rbx),%%xmm3				\n\t"\
		"subpd		%%xmm2,%%xmm0				\n\t"\
		"subpd		%%xmm3,%%xmm1				\n\t"\
		"addpd		%%xmm2,%%xmm2				\n\t"\
		"addpd		%%xmm3,%%xmm3				\n\t"\
		"addpd		%%xmm0,%%xmm2				\n\t"\
		"addpd		%%xmm1,%%xmm3				\n\t"\
		"movq		%[__add1],%%rbx				\n\t"\
		"movq		%[__c3],%%rcx				\n\t"\
		"movq		%[__c11],%%rdx				\n\t"\
		"subpd		%%xmm6,%%xmm0				\n\t"\
		"subpd		%%xmm7,%%xmm1				\n\t"\
		"addpd		%%xmm6,%%xmm6				\n\t"\
		"addpd		%%xmm7,%%xmm7				\n\t"\
		"addpd		%%xmm0,%%xmm6				\n\t"\
		"addpd		%%xmm1,%%xmm7				\n\t"\
		"movaps		%%xmm0,     (%%rax)			\n\t"\
		"movaps		%%xmm1,0x010(%%rax)			\n\t"\
		"movaps		%%xmm6,%%xmm0				\n\t"\
		"movaps		%%xmm7,%%xmm1				\n\t"\
		"mulpd		    (%%rcx),%%xmm6			\n\t"\
		"mulpd		    (%%rcx),%%xmm7			\n\t"\
		"mulpd		0x10(%%rcx),%%xmm0			\n\t"\
		"mulpd		0x10(%%rcx),%%xmm1			\n\t"\
		"subpd		%%xmm0,%%xmm7				\n\t"\
		"addpd		%%xmm1,%%xmm6				\n\t"\
		"movaps		%%xmm7,0x30(%%rbx)			\n\t"\
		"movaps		%%xmm6,0x20(%%rbx)			\n\t"\
		"movaps		     (%%rax),%%xmm6			\n\t"\
		"movaps		0x010(%%rax),%%xmm7			\n\t"\
		"movaps		%%xmm6,%%xmm0				\n\t"\
		"movaps		%%xmm7,%%xmm1				\n\t"\
		"mulpd		    (%%rdx),%%xmm6			\n\t"\
		"mulpd		    (%%rdx),%%xmm7			\n\t"\
		"mulpd		0x10(%%rdx),%%xmm0			\n\t"\
		"mulpd		0x10(%%rdx),%%xmm1			\n\t"\
		"subpd		%%xmm0,%%xmm7				\n\t"\
		"addpd		%%xmm1,%%xmm6				\n\t"\
		"movaps		%%xmm7,0xb0(%%rbx)			\n\t"\
		"movaps		%%xmm6,0xa0(%%rbx)			\n\t"\
		"movq		%[__c7],%%rcx				\n\t"\
		"movq		%[__c15],%%rdx				\n\t"\
		"subpd		%%xmm5,%%xmm2				\n\t"\
		"subpd		%%xmm4,%%xmm3				\n\t"\
		"addpd		%%xmm5,%%xmm5				\n\t"\
		"addpd		%%xmm4,%%xmm4				\n\t"\
		"addpd		%%xmm2,%%xmm5				\n\t"\
		"addpd		%%xmm3,%%xmm4				\n\t"\
		"movaps		%%xmm5,%%xmm0				\n\t"\
		"movaps		%%xmm3,%%xmm1				\n\t"\
		"mulpd		    (%%rcx),%%xmm5			\n\t"\
		"mulpd		    (%%rcx),%%xmm3			\n\t"\
		"mulpd		0x10(%%rcx),%%xmm0			\n\t"\
		"mulpd		0x10(%%rcx),%%xmm1			\n\t"\
		"subpd		%%xmm0,%%xmm3				\n\t"\
		"addpd		%%xmm1,%%xmm5				\n\t"\
		"movaps		%%xmm3,0x70(%%rbx)			\n\t"\
		"movaps		%%xmm5,0x60(%%rbx)			\n\t"\
		"movaps		%%xmm2,%%xmm0				\n\t"\
		"movaps		%%xmm4,%%xmm1				\n\t"\
		"mulpd		    (%%rdx),%%xmm2			\n\t"\
		"mulpd		    (%%rdx),%%xmm4			\n\t"\
		"mulpd		0x10(%%rdx),%%xmm0			\n\t"\
		"mulpd		0x10(%%rdx),%%xmm1			\n\t"\
		"subpd		%%xmm0,%%xmm4				\n\t"\
		"addpd		%%xmm1,%%xmm2				\n\t"\
		"movaps		%%xmm4,0xf0(%%rbx)			\n\t"\
		"movaps		%%xmm2,0xe0(%%rbx)			\n\t"\
		"/*...Block 1: t1,9,17,25 -> r1,5,3,7: */		\n\t"\
		"movq		%[__r1],%%rax				\n\t"\
		"movaps		     (%%rax),%%xmm0			\n\t"\
		"movaps		0x010(%%rax),%%xmm1			\n\t"\
		"movaps		0x040(%%rax),%%xmm2			\n\t"\
		"movaps		0x050(%%rax),%%xmm3			\n\t"\
		"subpd		0x040(%%rax),%%xmm0			\n\t"\
		"subpd		0x050(%%rax),%%xmm1			\n\t"\
		"addpd		     (%%rax),%%xmm2			\n\t"\
		"addpd		0x010(%%rax),%%xmm3			\n\t"\
		"movaps		0x020(%%rax),%%xmm4			\n\t"\
		"movaps		0x030(%%rax),%%xmm5			\n\t"\
		"movaps		0x060(%%rax),%%xmm6			\n\t"\
		"movaps		0x070(%%rax),%%xmm7			\n\t"\
		"subpd		0x060(%%rax),%%xmm4			\n\t"\
		"subpd		0x070(%%rax),%%xmm5			\n\t"\
		"addpd		0x020(%%rax),%%xmm6			\n\t"\
		"addpd		0x030(%%rax),%%xmm7			\n\t"\
		"movq		%[__add0],%%rax				\n\t"\
		"movq		%[__c8],%%rdx				\n\t"\
		"addpd		%%xmm6,%%xmm2				\n\t"\
		"addpd		%%xmm7,%%xmm3				\n\t"\
		"movaps		%%xmm2,    (%%rax)			\n\t"\
		"movaps		%%xmm3,0x10(%%rax)			\n\t"\
		"addpd		%%xmm6,%%xmm6				\n\t"\
		"addpd		%%xmm7,%%xmm7				\n\t"\
		"subpd		%%xmm6,%%xmm2				\n\t"\
		"subpd		%%xmm7,%%xmm3				\n\t"\
		"movaps		%%xmm2,%%xmm6				\n\t"\
		"movaps		%%xmm3,%%xmm7				\n\t"\
		"mulpd		    (%%rdx),%%xmm2			\n\t"\
		"mulpd		    (%%rdx),%%xmm3			\n\t"\
		"mulpd		0x10(%%rdx),%%xmm6			\n\t"\
		"mulpd		0x10(%%rdx),%%xmm7			\n\t"\
		"subpd		%%xmm6,%%xmm3				\n\t"\
		"addpd		%%xmm7,%%xmm2				\n\t"\
		"movq		%[__add1],%%rcx				\n\t"\
		"movaps		%%xmm3,%%xmm7				\n\t"\
		"movaps		%%xmm2,%%xmm6				\n\t"\
		"unpckhpd	0x90(%%rcx),%%xmm7			\n\t"\
		"unpcklpd	0x90(%%rcx),%%xmm3			\n\t"\
		"movaps		%%xmm7,0x90(%%rcx)			\n\t"\
		"unpckhpd	0x80(%%rcx),%%xmm6			\n\t"\
		"unpcklpd	0x80(%%rcx),%%xmm2			\n\t"\
		"movaps		%%xmm6,0x80(%%rcx)			\n\t"\
		"movaps		%%xmm3,0x90(%%rax)			\n\t"\
		"movaps		%%xmm2,0x80(%%rax)			\n\t"\
		"movaps		0x10(%%rax),%%xmm3			\n\t"\
		"movaps		(%%rax),%%xmm2				\n\t"\
		"movaps		%%xmm3,%%xmm7				\n\t"\
		"movaps		%%xmm2,%%xmm6				\n\t"\
		"unpckhpd	0x10(%%rcx),%%xmm7			\n\t"\
		"unpcklpd	0x10(%%rcx),%%xmm3			\n\t"\
		"movaps		%%xmm7,0x10(%%rcx)			\n\t"\
		"unpckhpd	(%%rcx),%%xmm6				\n\t"\
		"unpcklpd	(%%rcx),%%xmm2				\n\t"\
		"movaps		%%xmm6,    (%%rcx)			\n\t"\
		"movaps		%%xmm3,0x10(%%rax)			\n\t"\
		"movaps		%%xmm2,    (%%rax)			\n\t"\
		"movq		%[__c4],%%rcx				\n\t"\
		"movq		%[__c12],%%rdx				\n\t"\
		"addpd		%%xmm5,%%xmm0				\n\t"\
		"subpd		%%xmm4,%%xmm1				\n\t"\
		"movaps		%%xmm0,%%xmm2				\n\t"\
		"movaps		%%xmm1,%%xmm3				\n\t"\
		"addpd		%%xmm5,%%xmm5				\n\t"\
		"addpd		%%xmm4,%%xmm4				\n\t"\
		"movaps		%%xmm0,%%xmm6				\n\t"\
		"movaps		%%xmm1,%%xmm7				\n\t"\
		"mulpd		    (%%rcx),%%xmm2			\n\t"\
		"mulpd		    (%%rcx),%%xmm3			\n\t"\
		"mulpd		0x10(%%rcx),%%xmm6			\n\t"\
		"mulpd		0x10(%%rcx),%%xmm7			\n\t"\
		"subpd		%%xmm6,%%xmm3				\n\t"\
		"addpd		%%xmm7,%%xmm2				\n\t"\
		"movq		%[__add1],%%rcx				\n\t"\
		"movaps		%%xmm3,%%xmm7				\n\t"\
		"movaps		%%xmm2,%%xmm6				\n\t"\
		"unpckhpd	0x50(%%rcx),%%xmm7			\n\t"\
		"unpcklpd	0x50(%%rcx),%%xmm3			\n\t"\
		"movaps		%%xmm7,0x50(%%rcx)			\n\t"\
		"unpckhpd	0x40(%%rcx),%%xmm6			\n\t"\
		"unpcklpd	0x40(%%rcx),%%xmm2			\n\t"\
		"movaps		%%xmm6,0x40(%%rcx)			\n\t"\
		"movaps		%%xmm3,0x50(%%rax)			\n\t"\
		"movaps		%%xmm2,0x40(%%rax)			\n\t"\
		"subpd		%%xmm5,%%xmm0				\n\t"\
		"addpd		%%xmm4,%%xmm1				\n\t"\
		"movaps		%%xmm0,%%xmm6				\n\t"\
		"movaps		%%xmm1,%%xmm7				\n\t"\
		"mulpd		    (%%rdx),%%xmm0			\n\t"\
		"mulpd		    (%%rdx),%%xmm1			\n\t"\
		"mulpd		0x10(%%rdx),%%xmm6			\n\t"\
		"mulpd		0x10(%%rdx),%%xmm7			\n\t"\
		"subpd		%%xmm6,%%xmm1				\n\t"\
		"addpd		%%xmm7,%%xmm0				\n\t"\
		"movq		%[__add1],%%rcx				\n\t"\
		"movaps		%%xmm1,%%xmm7				\n\t"\
		"movaps		%%xmm0,%%xmm6				\n\t"\
		"unpckhpd	0xd0(%%rcx),%%xmm7			\n\t"\
		"unpcklpd	0xd0(%%rcx),%%xmm1			\n\t"\
		"movaps		%%xmm7,0xd0(%%rcx)			\n\t"\
		"unpckhpd	0xc0(%%rcx),%%xmm6			\n\t"\
		"unpcklpd	0xc0(%%rcx),%%xmm0			\n\t"\
		"movaps		%%xmm6,0xc0(%%rcx)			\n\t"\
		"movaps		%%xmm1,0xd0(%%rax)			\n\t"\
		"movaps		%%xmm0,0xc0(%%rax)			\n\t"\
		"/*...Block 2: t5,13,21,29 -> r17,21,19,23: */	\n\t"\
		"movq		%[__r17],%%rax				\n\t"\
		"movq		%[__isrt2],%%rbx			\n\t"\
		"movaps		(%%rbx),%%xmm2				\n\t"\
		"movaps		0x020(%%rax),%%xmm4			\n\t"\
		"movaps		0x030(%%rax),%%xmm5			\n\t"\
		"movaps		0x060(%%rax),%%xmm0			\n\t"\
		"movaps		0x070(%%rax),%%xmm1			\n\t"\
		"addpd		0x030(%%rax),%%xmm4			\n\t"\
		"subpd		0x020(%%rax),%%xmm5			\n\t"\
		"subpd		0x070(%%rax),%%xmm0			\n\t"\
		"addpd		0x060(%%rax),%%xmm1			\n\t"\
		"mulpd		%%xmm2,%%xmm4				\n\t"\
		"mulpd		%%xmm2,%%xmm5				\n\t"\
		"mulpd		%%xmm2,%%xmm0				\n\t"\
		"mulpd		%%xmm2,%%xmm1				\n\t"\
		"movaps		%%xmm4,%%xmm6				\n\t"\
		"movaps		%%xmm5,%%xmm7				\n\t"\
		"subpd		%%xmm0,%%xmm4				\n\t"\
		"subpd		%%xmm1,%%xmm5				\n\t"\
		"addpd		%%xmm0,%%xmm6				\n\t"\
		"addpd		%%xmm1,%%xmm7				\n\t"\
		"movaps		     (%%rax),%%xmm0			\n\t"\
		"movaps		0x010(%%rax),%%xmm1			\n\t"\
		"movaps		0x040(%%rax),%%xmm2			\n\t"\
		"movaps		0x050(%%rax),%%xmm3			\n\t"\
		"subpd		0x050(%%rax),%%xmm0			\n\t"\
		"subpd		0x040(%%rax),%%xmm1			\n\t"\
		"addpd		     (%%rax),%%xmm3			\n\t"\
		"addpd		0x010(%%rax),%%xmm2			\n\t"\
		"movq		%[__add0],%%rbx				\n\t"\
		"movq		%[__c2],%%rcx				\n\t"\
		"movq		%[__c10],%%rdx				\n\t"\
		"subpd		%%xmm4,%%xmm3				\n\t"\
		"subpd		%%xmm5,%%xmm1				\n\t"\
		"addpd		%%xmm4,%%xmm4				\n\t"\
		"addpd		%%xmm5,%%xmm5				\n\t"\
		"addpd		%%xmm3,%%xmm4				\n\t"\
		"addpd		%%xmm1,%%xmm5				\n\t"\
		"movaps		%%xmm3,     (%%rax)			\n\t"\
		"movaps		%%xmm1,0x010(%%rax)			\n\t"\
		"movaps		%%xmm4,%%xmm3				\n\t"\
		"movaps		%%xmm5,%%xmm1				\n\t"\
		"mulpd		    (%%rcx),%%xmm4			\n\t"\
		"mulpd		    (%%rcx),%%xmm5			\n\t"\
		"mulpd		0x10(%%rcx),%%xmm3			\n\t"\
		"mulpd		0x10(%%rcx),%%xmm1			\n\t"\
		"subpd		%%xmm3,%%xmm5				\n\t"\
		"addpd		%%xmm1,%%xmm4				\n\t"\
		"movq		%[__add1],%%rcx				\n\t"\
		"movaps		%%xmm5,%%xmm3				\n\t"\
		"movaps		%%xmm4,%%xmm1				\n\t"\
		"unpckhpd	0x30(%%rcx),%%xmm3			\n\t"\
		"unpcklpd	0x30(%%rcx),%%xmm5			\n\t"\
		"movaps		%%xmm3,0x30(%%rcx)			\n\t"\
		"unpckhpd	0x20(%%rcx),%%xmm1			\n\t"\
		"unpcklpd	0x20(%%rcx),%%xmm4			\n\t"\
		"movaps		%%xmm1,0x20(%%rcx)			\n\t"\
		"movaps		%%xmm5,0x30(%%rbx)			\n\t"\
		"movaps		%%xmm4,0x20(%%rbx)			\n\t"\
		"movaps		     (%%rax),%%xmm4			\n\t"\
		"movaps		0x010(%%rax),%%xmm5			\n\t"\
		"movaps		%%xmm4,%%xmm3				\n\t"\
		"movaps		%%xmm5,%%xmm1				\n\t"\
		"mulpd		    (%%rdx),%%xmm4			\n\t"\
		"mulpd		    (%%rdx),%%xmm5			\n\t"\
		"mulpd		0x10(%%rdx),%%xmm3			\n\t"\
		"mulpd		0x10(%%rdx),%%xmm1			\n\t"\
		"subpd		%%xmm3,%%xmm5				\n\t"\
		"addpd		%%xmm1,%%xmm4				\n\t"\
		"movaps		%%xmm5,%%xmm3				\n\t"\
		"movaps		%%xmm4,%%xmm1				\n\t"\
		"unpckhpd	0xb0(%%rcx),%%xmm3			\n\t"\
		"unpcklpd	0xb0(%%rcx),%%xmm5			\n\t"\
		"movaps		%%xmm3,0xb0(%%rcx)			\n\t"\
		"unpckhpd	0xa0(%%rcx),%%xmm1			\n\t"\
		"unpcklpd	0xa0(%%rcx),%%xmm4			\n\t"\
		"movaps		%%xmm1,0xa0(%%rcx)			\n\t"\
		"movaps		%%xmm5,0xb0(%%rbx)			\n\t"\
		"movaps		%%xmm4,0xa0(%%rbx)			\n\t"\
		"movq		%[__c6],%%rcx				\n\t"\
		"movq		%[__c14],%%rdx				\n\t"\
		"subpd		%%xmm7,%%xmm0				\n\t"\
		"subpd		%%xmm6,%%xmm2				\n\t"\
		"addpd		%%xmm7,%%xmm7				\n\t"\
		"addpd		%%xmm6,%%xmm6				\n\t"\
		"addpd		%%xmm0,%%xmm7				\n\t"\
		"addpd		%%xmm2,%%xmm6				\n\t"\
		"movaps		%%xmm7,%%xmm4				\n\t"\
		"movaps		%%xmm2,%%xmm5				\n\t"\
		"mulpd		    (%%rcx),%%xmm7			\n\t"\
		"mulpd		    (%%rcx),%%xmm2			\n\t"\
		"mulpd		0x10(%%rcx),%%xmm4			\n\t"\
		"mulpd		0x10(%%rcx),%%xmm5			\n\t"\
		"subpd		%%xmm4,%%xmm2				\n\t"\
		"addpd		%%xmm5,%%xmm7				\n\t"\
		"movq		%[__add1],%%rcx				\n\t"\
		"movaps		%%xmm2,%%xmm5				\n\t"\
		"movaps		%%xmm7,%%xmm4				\n\t"\
		"unpckhpd	0x70(%%rcx),%%xmm5			\n\t"\
		"unpcklpd	0x70(%%rcx),%%xmm2			\n\t"\
		"movaps		%%xmm5,0x70(%%rcx)			\n\t"\
		"unpckhpd	0x60(%%rcx),%%xmm4			\n\t"\
		"unpcklpd	0x60(%%rcx),%%xmm7			\n\t"\
		"movaps		%%xmm4,0x60(%%rcx)			\n\t"\
		"movaps		%%xmm2,0x70(%%rbx)			\n\t"\
		"movaps		%%xmm7,0x60(%%rbx)			\n\t"\
		"movaps		%%xmm0,%%xmm4				\n\t"\
		"movaps		%%xmm6,%%xmm5				\n\t"\
		"mulpd		    (%%rdx),%%xmm0			\n\t"\
		"mulpd		    (%%rdx),%%xmm6			\n\t"\
		"mulpd		0x10(%%rdx),%%xmm4			\n\t"\
		"mulpd		0x10(%%rdx),%%xmm5			\n\t"\
		"subpd		%%xmm4,%%xmm6				\n\t"\
		"addpd		%%xmm5,%%xmm0				\n\t"\
		"movaps		%%xmm6,%%xmm5				\n\t"\
		"movaps		%%xmm0,%%xmm4				\n\t"\
		"unpckhpd	0xf0(%%rcx),%%xmm5			\n\t"\
		"unpcklpd	0xf0(%%rcx),%%xmm6			\n\t"\
		"movaps		%%xmm5,0xf0(%%rcx)			\n\t"\
		"unpckhpd	0xe0(%%rcx),%%xmm4			\n\t"\
		"unpcklpd	0xe0(%%rcx),%%xmm0			\n\t"\
		"movaps		%%xmm4,0xe0(%%rcx)			\n\t"\
		"movaps		%%xmm6,0xf0(%%rbx)			\n\t"\
		"movaps		%%xmm0,0xe0(%%rbx)			\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__add1] "m" (Xadd1)\
		 ,[__r1] "m" (Xr1)\
		 ,[__r9] "m" (Xr9)\
		 ,[__r17] "m" (Xr17)\
		 ,[__r25] "m" (Xr25)\
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__cc0] "m" (Xcc0)\
		 ,[__c1] "m" (Xc1)\
		 ,[__c2] "m" (Xc2)\
		 ,[__c3] "m" (Xc3)\
		 ,[__c4] "m" (Xc4)\
		 ,[__c5] "m" (Xc5)\
		 ,[__c6] "m" (Xc6)\
		 ,[__c7] "m" (Xc7)\
		 ,[__c8] "m" (Xc8)\
		 ,[__c9] "m" (Xc9)\
		 ,[__c10] "m" (Xc10)\
		 ,[__c11] "m" (Xc11)\
		 ,[__c12] "m" (Xc12)\
		 ,[__c13] "m" (Xc13)\
		 ,[__c14] "m" (Xc14)\
		 ,[__c15] "m" (Xc15)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
	);\
	}

  #else // USE_64BIT_ASM_STYLE - Deeper 64-bit-ified version of the above 32-bit ASM macros, using all of xmm0-15

	#define SSE2_RADIX16_WRAPPER_DIF(Xadd0,Xadd1,Xr1,Xr9,Xr17,Xr25,Xisrt2,Xcc0,Xc1,Xc2,Xc3,Xc4,Xc5,Xc6,Xc7,Xc8,Xc9,Xc10,Xc11,Xc12,Xc13,Xc14,Xc15)\
	{\
	__asm__ volatile (\
	"/*************************************************************/\n\t"\
	"/*                  1st set of inputs:                       */\n\t"\
	"/*************************************************************/\n\t"\
		"/*...Block 1: */										\n\t	/*...Block 2: */					\n\t"\
		"movq		%[__add0],%%rax								\n\t	movq		%[__r9] ,%%r10			\n\t"\
		"movq		%[__add1],%%rbx								\n\t	movq		%[__c2] ,%%r11			\n\t"\
		"movq		%[__r1] ,%%rcx								\n\t	movaps		0x20(%%rax),%%xmm14		\n\t"\
		"movq		%[__c4] ,%%rdx								\n\t	movaps		0x20(%%rax),%%xmm8		\n\t"\
		"movaps		0x40(%%rax),%%xmm6							\n\t	unpckhpd	0x20(%%rbx),%%xmm14		\n\t"\
		"movaps		0x40(%%rbx),%%xmm2							\n\t	unpcklpd	0x20(%%rbx),%%xmm8		\n\t"\
		"movaps		%%xmm6,%%xmm0								\n\t	movaps			%%xmm14,0x100(%%r10)\n\t"\
		"movaps		%%xmm2,%%xmm3								\n\t	movaps		0x30(%%rax),%%xmm15		\n\t"\
		"unpckhpd	%%xmm2,%%xmm6								\n\t	movaps		0x30(%%rax),%%xmm9		\n\t"\
		"unpcklpd	%%xmm3,%%xmm0								\n\t	unpckhpd	0x30(%%rbx),%%xmm15		\n\t"\
		"movaps		%%xmm6,0x140(%%rcx)							\n\t	unpcklpd	0x30(%%rbx),%%xmm9		\n\t"\
		"movaps		0x50(%%rax),%%xmm7							\n\t	movaps		%%xmm15,0x110(%%r10)	\n\t"\
		"movaps		0x50(%%rbx),%%xmm4							\n\t	movaps		%%xmm8,%%xmm10			\n\t"\
		"movaps		%%xmm7,%%xmm1								\n\t	movaps		%%xmm9,%%xmm11			\n\t"\
		"movaps		%%xmm4,%%xmm5								\n\t	mulpd		    (%%r11),%%xmm8		\n\t"\
		"unpckhpd	%%xmm4,%%xmm7								\n\t	mulpd		    (%%r11),%%xmm9		\n\t"\
		"unpcklpd	%%xmm5,%%xmm1								\n\t	mulpd		0x10(%%r11),%%xmm10		\n\t"\
		"movaps		%%xmm7,0x150(%%rcx)							\n\t	mulpd		0x10(%%r11),%%xmm11		\n\t"\
		"movaps		%%xmm0,%%xmm2								\n\t	addpd		%%xmm10,%%xmm9			\n\t"\
		"movaps		%%xmm1,%%xmm3								\n\t	subpd		%%xmm11,%%xmm8			\n\t"\
		"/* Rest identical to code in radix16_dif_pass: */		\n\t	movaps		%%xmm9,%%xmm11			\n\t"\
		"mulpd		    (%%rdx),%%xmm0							\n\t	movaps		%%xmm8,%%xmm10			\n\t"\
		"mulpd		    (%%rdx),%%xmm1							\n\t	movq		%[__c10],%%r11			\n\t"\
		"mulpd		0x10(%%rdx),%%xmm2							\n\t	movaps		0xa0(%%rax),%%xmm14		\n\t"\
		"mulpd		0x10(%%rdx),%%xmm3							\n\t	movaps		0xa0(%%rax),%%xmm12		\n\t"\
		"addpd		%%xmm2,%%xmm1								\n\t	unpckhpd	0xa0(%%rbx),%%xmm14		\n\t"\
		"subpd		%%xmm3,%%xmm0								\n\t	unpcklpd	0xa0(%%rbx),%%xmm12		\n\t"\
		"movaps		%%xmm1,%%xmm3								\n\t	movaps			%%xmm14,0x120(%%r10)\n\t"\
		"movaps		%%xmm0,%%xmm2								\n\t	movaps		0xb0(%%rax),%%xmm15		\n\t"\
		"movq		%[__c12],%%rdx								\n\t	movaps		0xb0(%%rax),%%xmm13		\n\t"\
		"movaps		0xc0(%%rax),%%xmm6							\n\t	unpckhpd	0xb0(%%rbx),%%xmm15		\n\t"\
		"movaps		0xc0(%%rax),%%xmm4							\n\t	unpcklpd	0xb0(%%rbx),%%xmm13		\n\t"\
		"unpckhpd	0xc0(%%rbx),%%xmm6							\n\t	movaps			%%xmm15,0x130(%%r10)\n\t"\
		"unpcklpd	0xc0(%%rbx),%%xmm4							\n\t	movaps			%%xmm12,%%xmm14		\n\t"\
		"movaps		%%xmm6,0x160(%%rcx)							\n\t	movaps			%%xmm13,%%xmm15		\n\t"\
		"movaps		0xd0(%%rax),%%xmm7							\n\t	mulpd		    (%%r11),%%xmm12		\n\t"\
		"movaps		0xd0(%%rax),%%xmm5							\n\t	mulpd		    (%%r11),%%xmm13		\n\t"\
		"unpckhpd	0xd0(%%rbx),%%xmm7							\n\t	mulpd		0x10(%%r11),%%xmm14		\n\t"\
		"unpcklpd	0xd0(%%rbx),%%xmm5							\n\t	mulpd		0x10(%%r11),%%xmm15		\n\t"\
		"movaps		%%xmm7,0x170(%%rcx)							\n\t	addpd		%%xmm14,%%xmm13			\n\t"\
		"movaps		%%xmm4,%%xmm6								\n\t	subpd		%%xmm15,%%xmm12			\n\t"\
		"movaps		%%xmm5,%%xmm7								\n\t	addpd		%%xmm12,%%xmm8			\n\t"\
		"mulpd		    (%%rdx),%%xmm4							\n\t	addpd		%%xmm13,%%xmm9			\n\t"\
		"mulpd		    (%%rdx),%%xmm5							\n\t	subpd		%%xmm12,%%xmm10			\n\t"\
		"mulpd		0x10(%%rdx),%%xmm6							\n\t	subpd		%%xmm13,%%xmm11			\n\t"\
		"mulpd		0x10(%%rdx),%%xmm7							\n\t	movq		%[__c14],%%r11			\n\t"\
		"addpd		%%xmm6,%%xmm5								\n\t	movaps		0xe0(%%rax),%%xmm14		\n\t"\
		"subpd		%%xmm7,%%xmm4								\n\t	movaps		0xe0(%%rax),%%xmm12		\n\t"\
		"addpd		%%xmm4,%%xmm0								\n\t	unpckhpd	0xe0(%%rbx),%%xmm14		\n\t"\
		"addpd		%%xmm5,%%xmm1								\n\t	unpcklpd	0xe0(%%rbx),%%xmm12		\n\t"\
		"subpd		%%xmm4,%%xmm2								\n\t	movaps			%%xmm14,0x160(%%r10)\n\t"\
		"subpd		%%xmm5,%%xmm3								\n\t	movaps		0xf0(%%rax),%%xmm15		\n\t"\
		"movq		%[__c8] ,%%rdx								\n\t	movaps		0xf0(%%rax),%%xmm13		\n\t"\
		"movaps		0x80(%%rax),%%xmm6							\n\t	unpckhpd	0xf0(%%rbx),%%xmm15		\n\t"\
		"movaps		0x80(%%rax),%%xmm4							\n\t	unpcklpd	0xf0(%%rbx),%%xmm13		\n\t"\
		"unpckhpd	0x80(%%rbx),%%xmm6							\n\t	movaps			%%xmm15,0x170(%%r10)\n\t"\
		"unpcklpd	0x80(%%rbx),%%xmm4							\n\t	movaps			%%xmm12,%%xmm14		\n\t"\
		"movaps			%%xmm6,0x120(%%rcx)						\n\t	movaps			%%xmm13,%%xmm15		\n\t"\
		"movaps		0x90(%%rax),%%xmm7							\n\t	mulpd		    (%%r11),%%xmm12		\n\t"\
		"movaps		0x90(%%rax),%%xmm5							\n\t	mulpd		    (%%r11),%%xmm13		\n\t"\
		"unpckhpd	0x90(%%rbx),%%xmm7							\n\t	mulpd		0x10(%%r11),%%xmm14		\n\t"\
		"unpcklpd	0x90(%%rbx),%%xmm5							\n\t	mulpd		0x10(%%r11),%%xmm15		\n\t"\
		"movaps		%%xmm7,0x130(%%rcx)							\n\t	addpd		%%xmm14,%%xmm13			\n\t"\
		"movaps		%%xmm4,%%xmm6								\n\t	subpd		%%xmm15,%%xmm12			\n\t"\
		"movaps		%%xmm5,%%xmm7								\n\t	movaps		%%xmm13,0x070(%%r10)	\n\t"\
		"mulpd		    (%%rdx),%%xmm4							\n\t	movaps		%%xmm12,0x060(%%r10)	\n\t"\
		"mulpd		    (%%rdx),%%xmm5							\n\t	movq		%[__c6] ,%%r11			\n\t"\
		"mulpd		0x10(%%rdx),%%xmm6							\n\t	movaps		0x60(%%rax),%%xmm14		\n\t"\
		"mulpd		0x10(%%rdx),%%xmm7							\n\t	movaps		0x60(%%rax),%%xmm12		\n\t"\
		"addpd		%%xmm6,%%xmm5								\n\t	unpckhpd	0x60(%%rbx),%%xmm14		\n\t"\
		"subpd		%%xmm7,%%xmm4								\n\t	unpcklpd	0x60(%%rbx),%%xmm12		\n\t"\
		"movaps		    (%%rax),%%xmm6							\n\t	movaps			%%xmm14,0x140(%%r10)\n\t"\
		"movaps		    (%%rax),%%xmm7							\n\t	movaps		0x70(%%rax),%%xmm15		\n\t"\
		"unpckhpd	    (%%rbx),%%xmm6							\n\t	movaps		0x70(%%rax),%%xmm13		\n\t"\
		"unpcklpd	    (%%rbx),%%xmm7							\n\t	unpckhpd	0x70(%%rbx),%%xmm15		\n\t"\
		"movaps		%%xmm6,0x100(%%rcx)							\n\t	unpcklpd	0x70(%%rbx),%%xmm13		\n\t"\
		"movaps		%%xmm7,     (%%rcx)							\n\t	movaps			%%xmm15,0x150(%%r10)\n\t"\
		"movaps		0x10(%%rax),%%xmm6							\n\t	movaps			%%xmm12,%%xmm14		\n\t"\
		"movaps		0x10(%%rax),%%xmm7							\n\t	movaps			%%xmm13,%%xmm15		\n\t"\
		"unpckhpd	0x10(%%rbx),%%xmm6							\n\t	mulpd		    (%%r11),%%xmm12		\n\t"\
		"unpcklpd	0x10(%%rbx),%%xmm7							\n\t	mulpd		    (%%r11),%%xmm13		\n\t"\
		"movaps		%%xmm6,0x110(%%rcx)							\n\t	mulpd		0x10(%%r11),%%xmm14		\n\t"\
		"movaps		    (%%rcx),%%xmm6							\n\t	mulpd		0x10(%%r11),%%xmm15		\n\t"\
		"subpd		%%xmm4,%%xmm6								\n\t	addpd		%%xmm14,%%xmm13			\n\t"\
		"subpd		%%xmm5,%%xmm7								\n\t	subpd		%%xmm15,%%xmm12			\n\t"\
		"addpd		%%xmm4,%%xmm4								\n\t	movaps		%%xmm13,%%xmm15			\n\t"\
		"addpd		%%xmm5,%%xmm5								\n\t	movaps		%%xmm12,%%xmm14			\n\t"\
		"addpd		%%xmm6,%%xmm4								\n\t	subpd		0x060(%%r10),%%xmm12	\n\t"\
		"addpd		%%xmm7,%%xmm5								\n\t	subpd		0x070(%%r10),%%xmm13	\n\t"\
		"/* Finish radix-4 butterfly: */						\n\t	addpd		0x060(%%r10),%%xmm14	\n\t"\
		"subpd		%%xmm0,%%xmm4								\n\t	addpd		0x070(%%r10),%%xmm15	\n\t"\
		"subpd		%%xmm1,%%xmm5								\n\t	/* Finish radix-4 butterfly: */		\n\t"\
		"movaps		%%xmm4,0x040(%%rcx)							\n\t	subpd		%%xmm14,%%xmm8			\n\t"\
		"movaps		%%xmm5,0x050(%%rcx)							\n\t	subpd		%%xmm15,%%xmm9			\n\t"\
		"addpd		%%xmm0,%%xmm0								\n\t	movaps		%%xmm8,0x040(%%r10)		\n\t"\
		"addpd		%%xmm1,%%xmm1								\n\t	movaps		%%xmm9,0x050(%%r10)		\n\t"\
		"addpd		%%xmm4,%%xmm0								\n\t	addpd		%%xmm14,%%xmm14			\n\t"\
		"addpd		%%xmm5,%%xmm1								\n\t	addpd		%%xmm15,%%xmm15			\n\t"\
		"movaps		%%xmm0,     (%%rcx)							\n\t	addpd		%%xmm8,%%xmm14			\n\t"\
		"movaps		%%xmm1,0x010(%%rcx)							\n\t	addpd		%%xmm9,%%xmm15			\n\t"\
		"subpd		%%xmm3,%%xmm6								\n\t	movaps		%%xmm14,     (%%r10)	\n\t"\
		"subpd		%%xmm2,%%xmm7								\n\t	movaps		%%xmm15,0x010(%%r10)	\n\t"\
		"movaps		%%xmm6,0x020(%%rcx)							\n\t	subpd		%%xmm13,%%xmm10			\n\t"\
		"movaps		%%xmm7,0x070(%%rcx)							\n\t	subpd		%%xmm12,%%xmm11			\n\t"\
		"addpd		%%xmm3,%%xmm3								\n\t	movaps		%%xmm10,0x020(%%r10)	\n\t"\
		"addpd		%%xmm2,%%xmm2								\n\t	movaps		%%xmm11,0x070(%%r10)	\n\t"\
		"addpd		%%xmm6,%%xmm3								\n\t	addpd		%%xmm13,%%xmm13			\n\t"\
		"addpd		%%xmm7,%%xmm2								\n\t	addpd		%%xmm12,%%xmm12			\n\t"\
		"movaps		%%xmm3,0x060(%%rcx)							\n\t	addpd		%%xmm10,%%xmm13			\n\t"\
		"movaps		%%xmm2,0x030(%%rcx)							\n\t	addpd		%%xmm11,%%xmm12			\n\t"\
		"														\n\t	movaps		%%xmm13,0x060(%%r10)	\n\t"\
		"														\n\t	movaps		%%xmm12,0x030(%%r10)	\n\t"\
	"/****************************************************************************************************/	\n\t"\
	"/* Each of the above 2 radix-4 blocks is unique, but the ensuing 2 blocks [operating on odd-indexed */	\n\t"\
	"/* elements from the unpck*pd commands which were stored to temporaries] can use a common macro:    */	\n\t"\
	"/****************************************************************************************************/	\n\t"\
		"/*...Block 3: */										\n\t	/*...Block 4: */					\n\t"\
		"/*	SSE2_RADIX4_DIF_4TWIDDLE(r17,r21,r19,r23, c1): */	\n\t	/* SSE2_RADIX4_DIF_4TWIDDLE(r25,r29,r27,r31, c3): */\n\t"\
		"/* Do the p0,p8 combo: */								\n\t	/* Do the p0,p8 combo: */			\n\t"\
		"movq		%[__r17],%%rax								\n\t	movq		%[__r25],%%r10			\n\t"\
		"movq		%[__c1] ,%%rbx								\n\t	movq		%[__c3] ,%%r11			\n\t"\
		"movq		%%rax   ,%%rcx								\n\t	movq		%%r10   ,%%r12			\n\t"\
		"addq		$0x20   ,%%rcx								\n\t	addq		$0x20   ,%%r12			\n\t"\
		"movaps		    (%%rax),%%xmm0							\n\t	movaps		    (%%r10),%%xmm8 		\n\t"\
		"movaps	        (%%rcx),%%xmm4							\n\t	movaps	        (%%r12),%%xmm12		\n\t"\
		"movaps		0x10(%%rax),%%xmm1							\n\t	movaps		0x10(%%r10),%%xmm9 		\n\t"\
		"movaps		0x10(%%rcx),%%xmm5							\n\t	movaps		0x10(%%r12),%%xmm13		\n\t"\
		"movaps		    (%%rbx),%%xmm6							\n\t	movaps		    (%%r11),%%xmm14		\n\t"\
		"movaps		0x10(%%rbx),%%xmm7							\n\t	movaps		0x10(%%r11),%%xmm15		\n\t"\
		"movaps		%%xmm0,%%xmm2								\n\t	movaps		%%xmm8 ,%%xmm10			\n\t"\
		"movaps		%%xmm1,%%xmm3								\n\t	movaps		%%xmm9 ,%%xmm11			\n\t"\
		"mulpd		%%xmm6,%%xmm0								\n\t	mulpd   	%%xmm14,%%xmm8 			\n\t"\
		"mulpd		%%xmm6,%%xmm1								\n\t	mulpd   	%%xmm14,%%xmm9 			\n\t"\
		"mulpd		%%xmm7,%%xmm2								\n\t	mulpd   	%%xmm15,%%xmm10			\n\t"\
		"mulpd		%%xmm7,%%xmm3								\n\t	mulpd   	%%xmm15,%%xmm11			\n\t"\
		"movaps		%%xmm4,%%xmm6								\n\t	movaps		%%xmm12,%%xmm14			\n\t"\
		"addpd		%%xmm2,%%xmm1								\n\t	addpd   	%%xmm10,%%xmm9 			\n\t"\
		"movaps		%%xmm5,%%xmm7								\n\t	movaps		%%xmm13,%%xmm15			\n\t"\
		"mulpd		0x20(%%rbx),%%xmm4							\n\t	mulpd		0x20(%%r11),%%xmm12		\n\t"\
		"subpd		%%xmm3,%%xmm0								\n\t	subpd		%%xmm11,%%xmm8 			\n\t"\
		"mulpd		0x20(%%rbx),%%xmm5							\n\t	mulpd		0x20(%%r11),%%xmm13		\n\t"\
		"mulpd		0x30(%%rbx),%%xmm6							\n\t	mulpd		0x30(%%r11),%%xmm14		\n\t"\
		"movaps		%%xmm0,%%xmm2								\n\t	movaps		%%xmm8 ,%%xmm10			\n\t"\
		"mulpd		0x30(%%rbx),%%xmm7							\n\t	mulpd		0x30(%%r11),%%xmm15		\n\t"\
		"addpd		%%xmm6,%%xmm5								\n\t	addpd		%%xmm14,%%xmm13			\n\t"\
		"movaps		%%xmm1,%%xmm3								\n\t	movaps		%%xmm9 ,%%xmm11			\n\t"\
		"subpd		%%xmm7,%%xmm4								\n\t	subpd		%%xmm15,%%xmm12			\n\t"\
		"addq		$0x40,%%rcx									\n\t	addq		$0x40	,%%r12			\n\t"\
		"addq		$0x60,%%rbx									\n\t	addq		$0x60	,%%r11			\n\t"\
		"movaps			(%%rcx),%%xmm6							\n\t	movaps		    (%%r12),%%xmm14		\n\t"\
		"movaps		0x10(%%rcx),%%xmm7							\n\t	movaps		0x10(%%r12),%%xmm15		\n\t"\
		"addpd		%%xmm4,%%xmm0								\n\t	addpd		%%xmm12,%%xmm8 			\n\t"\
		"addpd		%%xmm5,%%xmm1								\n\t	addpd		%%xmm13,%%xmm9 			\n\t"\
		"subpd		%%xmm4,%%xmm2								\n\t	subpd		%%xmm12,%%xmm10			\n\t"\
		"subpd		%%xmm5,%%xmm3								\n\t	subpd		%%xmm13,%%xmm11			\n\t"\
		"/* Do the p4,12 combo: */								\n\t	/* Do the p4,12 combo: */			\n\t"\
		"movaps		%%xmm6,%%xmm4								\n\t	movaps		%%xmm14,%%xmm12			\n\t"\
		"movaps		%%xmm7,%%xmm5								\n\t	movaps		%%xmm15,%%xmm13			\n\t"\
		"mulpd		    (%%rbx),%%xmm4							\n\t	mulpd		    (%%r11),%%xmm12		\n\t"\
		"mulpd		    (%%rbx),%%xmm5							\n\t	mulpd		    (%%r11),%%xmm13		\n\t"\
		"mulpd		0x10(%%rbx),%%xmm6							\n\t	mulpd		0x10(%%r11),%%xmm14		\n\t"\
		"mulpd		0x10(%%rbx),%%xmm7							\n\t	mulpd		0x10(%%r11),%%xmm15		\n\t"\
		"movq		%%rax,%%rdx									\n\t	movq		%%r10,%%r13				\n\t"\
		"addpd		%%xmm6,%%xmm5								\n\t	addpd		%%xmm14,%%xmm13			\n\t"\
		"subpd		%%xmm7,%%xmm4								\n\t	subpd		%%xmm15,%%xmm12			\n\t"\
		"movaps		%%xmm5,0x010(%%rdx)							\n\t	movaps		%%xmm13,0x010(%%r13)	\n\t"\
		"movaps		%%xmm4,     (%%rdx)							\n\t	movaps		%%xmm12,     (%%r13)	\n\t"\
		"addq		$0x40 ,%%rax								\n\t	addq		$0x40 ,%%r10			\n\t"\
		"subq		$0x20 ,%%rbx								\n\t	subq		$0x20 ,%%r11			\n\t"\
		"movaps		    (%%rax),%%xmm4							\n\t	movaps		    (%%r10),%%xmm12		\n\t"\
		"movaps		0x10(%%rax),%%xmm5							\n\t	movaps		0x10(%%r10),%%xmm13		\n\t"\
		"movaps			%%xmm4,%%xmm6							\n\t	movaps			%%xmm12,%%xmm14		\n\t"\
		"movaps			%%xmm5,%%xmm7							\n\t	movaps			%%xmm13,%%xmm15		\n\t"\
		"mulpd		    (%%rbx),%%xmm4							\n\t	mulpd		    (%%r11),%%xmm12		\n\t"\
		"mulpd		    (%%rbx),%%xmm5							\n\t	mulpd		    (%%r11),%%xmm13		\n\t"\
		"mulpd		0x10(%%rbx),%%xmm6							\n\t	mulpd		0x10(%%r11),%%xmm14		\n\t"\
		"mulpd		0x10(%%rbx),%%xmm7							\n\t	mulpd		0x10(%%r11),%%xmm15		\n\t"\
		"addpd		%%xmm6,%%xmm5								\n\t	addpd		%%xmm14,%%xmm13			\n\t"\
		"subpd		%%xmm7,%%xmm4								\n\t	subpd		%%xmm15,%%xmm12			\n\t"\
		"movaps		%%xmm5,%%xmm7								\n\t	movaps		%%xmm13,%%xmm15			\n\t"\
		"movaps		%%xmm4,%%xmm6								\n\t	movaps		%%xmm12,%%xmm14			\n\t"\
		"subpd		     (%%rdx),%%xmm4							\n\t	subpd		     (%%r13),%%xmm12	\n\t"\
		"subpd		0x010(%%rdx),%%xmm5							\n\t	subpd		0x010(%%r13),%%xmm13	\n\t"\
		"addpd		     (%%rdx),%%xmm6							\n\t	addpd		     (%%r13),%%xmm14	\n\t"\
		"addpd		0x010(%%rdx),%%xmm7							\n\t	addpd		0x010(%%r13),%%xmm15	\n\t"\
		"/* Finish radix-4 butterfly: */						\n\t	/* Finish radix-4 butterfly: */		\n\t"\
		"subpd		%%xmm6,%%xmm0								\n\t	subpd		%%xmm14,%%xmm8 			\n\t"\
		"subpd		%%xmm5,%%xmm2								\n\t	subpd		%%xmm13,%%xmm10			\n\t"\
		"subpd		%%xmm7,%%xmm1								\n\t	subpd		%%xmm15,%%xmm9 			\n\t"\
		"subpd		%%xmm4,%%xmm3								\n\t	subpd		%%xmm12,%%xmm11			\n\t"\
		"movaps		%%xmm0,0x040(%%rdx)							\n\t	movaps		%%xmm8 ,0x040(%%r13)	\n\t"\
		"movaps		%%xmm2,0x020(%%rdx)							\n\t	movaps		%%xmm10,0x020(%%r13)	\n\t"\
		"movaps		%%xmm1,0x050(%%rdx)							\n\t	movaps		%%xmm9 ,0x050(%%r13)	\n\t"\
		"movaps		%%xmm3,0x070(%%rdx)							\n\t	movaps		%%xmm11,0x070(%%r13)	\n\t"\
		"addpd		%%xmm6,%%xmm6								\n\t	addpd		%%xmm14,%%xmm14			\n\t"\
		"addpd		%%xmm5,%%xmm5								\n\t	addpd		%%xmm13,%%xmm13			\n\t"\
		"addpd		%%xmm7,%%xmm7								\n\t	addpd		%%xmm15,%%xmm15			\n\t"\
		"addpd		%%xmm4,%%xmm4								\n\t	addpd		%%xmm12,%%xmm12			\n\t"\
		"addpd		%%xmm0,%%xmm6								\n\t	addpd		%%xmm8 ,%%xmm14			\n\t"\
		"addpd		%%xmm2,%%xmm5								\n\t	addpd		%%xmm10,%%xmm13			\n\t"\
		"addpd		%%xmm1,%%xmm7								\n\t	addpd		%%xmm9 ,%%xmm15			\n\t"\
		"addpd		%%xmm3,%%xmm4								\n\t	addpd		%%xmm11,%%xmm12			\n\t"\
		"movaps		%%xmm6,     (%%rdx)							\n\t	movaps		%%xmm14,     (%%r13)	\n\t"\
		"movaps		%%xmm5,0x060(%%rdx)							\n\t	movaps		%%xmm13,0x060(%%r13)	\n\t"\
		"movaps		%%xmm7,0x010(%%rdx)							\n\t	movaps		%%xmm15,0x010(%%r13)	\n\t"\
		"movaps		%%xmm4,0x030(%%rdx)							\n\t	movaps		%%xmm12,0x030(%%r13)	\n\t"\
	"/**************************************************************************************/\n\t"\
	"/*...and now do four more radix-4 transforms, including the internal twiddle factors: */\n\t"\
	"/**************************************************************************************/\n\t"\
		"/*...Block 1: t1,9,17,25 */		\n\t		/*...Block 3: t5,13,21,29 */		\n\t"\
		"movq		%[__r1] ,%%rax						\n\t"\
		"movq		%[__r9] ,%%rbx						\n\t"\
		"movq		%[__r17],%%rcx						\n\t"\
		"movq		%[__r25],%%rdx			\n\t		movq		%[__isrt2],%%rsi		\n\t"\
		"movaps		     (%%rax),%%xmm0		\n\t		movaps		0x040(%%rax),%%xmm8		\n\t"\
		"movaps		0x010(%%rax),%%xmm1		\n\t		movaps		0x050(%%rax),%%xmm9		\n\t"\
		"movaps		     (%%rbx),%%xmm2		\n\t		movaps		0x0c0(%%rax),%%xmm10	\n\t"\
		"movaps		0x010(%%rbx),%%xmm3		\n\t		movaps		0x0d0(%%rax),%%xmm11	\n\t"\
		"subpd		     (%%rbx),%%xmm0		\n\t		subpd		0x0d0(%%rax),%%xmm8		\n\t"\
		"subpd		0x010(%%rbx),%%xmm1		\n\t		subpd		0x0c0(%%rax),%%xmm9		\n\t"\
		"addpd		     (%%rax),%%xmm2		\n\t		addpd		0x050(%%rax),%%xmm10	\n\t"\
		"addpd		0x010(%%rax),%%xmm3		\n\t		addpd		0x040(%%rax),%%xmm11	\n\t"\
		"movaps		     (%%rcx),%%xmm4		\n\t		movaps		0x140(%%rax),%%xmm12	\n\t"\
		"movaps		0x010(%%rcx),%%xmm5		\n\t		movaps		0x150(%%rax),%%xmm13	\n\t"\
		"movaps		     (%%rdx),%%xmm6		\n\t		movaps		0x1c0(%%rax),%%xmm14	\n\t"\
		"movaps		0x010(%%rdx),%%xmm7		\n\t		movaps		0x1d0(%%rax),%%xmm15	\n\t"\
		"subpd		     (%%rdx),%%xmm4		\n\t		subpd		0x150(%%rax),%%xmm12	\n\t"\
		"subpd		0x010(%%rdx),%%xmm5		\n\t		addpd		0x140(%%rax),%%xmm13	\n\t"\
		"addpd		     (%%rcx),%%xmm6		\n\t		mulpd		(%%rsi),%%xmm12			\n\t"\
		"addpd		0x010(%%rcx),%%xmm7		\n\t		mulpd		(%%rsi),%%xmm13			\n\t"\
		"subpd		%%xmm6,%%xmm2			\n\t		addpd		0x1d0(%%rax),%%xmm14	\n\t"\
		"subpd		%%xmm7,%%xmm3			\n\t		subpd		0x1c0(%%rax),%%xmm15	\n\t"\
		"addpd		%%xmm6,%%xmm6			\n\t		mulpd		(%%rsi),%%xmm14			\n\t"\
		"addpd		%%xmm7,%%xmm7			\n\t		mulpd		(%%rsi),%%xmm15			\n\t"\
		"movaps		%%xmm2,    (%%rcx)		\n\t		subpd		%%xmm14,%%xmm12			\n\t"\
		"movaps		%%xmm3,0x10(%%rcx)		\n\t		subpd		%%xmm15,%%xmm13			\n\t"\
		"addpd		%%xmm2,%%xmm6			\n\t		addpd		%%xmm14,%%xmm14			\n\t"\
		"addpd		%%xmm3,%%xmm7			\n\t		addpd		%%xmm15,%%xmm15			\n\t"\
		"movaps		%%xmm6,    (%%rax)		\n\t		addpd		%%xmm12,%%xmm14			\n\t"\
		"movaps		%%xmm7,0x10(%%rax)		\n\t		addpd		%%xmm13,%%xmm15			\n\t"\
		"subpd		%%xmm5,%%xmm0			\n\t		subpd		%%xmm12,%%xmm8			\n\t"\
		"subpd		%%xmm4,%%xmm1			\n\t		subpd		%%xmm13,%%xmm10			\n\t"\
		"addpd		%%xmm5,%%xmm5			\n\t		addpd		%%xmm12,%%xmm12			\n\t"\
		"addpd		%%xmm4,%%xmm4			\n\t		addpd		%%xmm13,%%xmm13			\n\t"\
		"movaps		%%xmm0,    (%%rbx)		\n\t		movaps		%%xmm8,0x40(%%rcx)		\n\t"\
		"movaps		%%xmm1,0x10(%%rdx)		\n\t		movaps		%%xmm10,0x50(%%rcx)		\n\t"\
		"addpd		%%xmm0,%%xmm5			\n\t		addpd		%%xmm8,%%xmm12			\n\t"\
		"addpd		%%xmm1,%%xmm4			\n\t		addpd		%%xmm10,%%xmm13			\n\t"\
		"movaps		%%xmm5,    (%%rdx)		\n\t		movaps		%%xmm12,0x40(%%rax)		\n\t"\
		"movaps		%%xmm4,0x10(%%rbx)		\n\t		movaps		%%xmm13,0x50(%%rax)		\n\t"\
		"												subpd		%%xmm15,%%xmm11			\n\t"\
		"												subpd		%%xmm14,%%xmm9			\n\t"\
		"												addpd		%%xmm15,%%xmm15			\n\t"\
		"												addpd		%%xmm14,%%xmm14			\n\t"\
		"												movaps		%%xmm11,0x40(%%rbx)		\n\t"\
		"												movaps		%%xmm9,0x50(%%rdx)		\n\t"\
		"												addpd		%%xmm11,%%xmm15			\n\t"\
		"												addpd		%%xmm9,%%xmm14			\n\t"\
		"												movaps		%%xmm15,0x40(%%rdx)		\n\t"\
		"												movaps		%%xmm14,0x50(%%rbx)		\n\t"\
		"/*...Block 2: t3,11,19,27 */		\n\t		/*...Block 4: t7,15,23,31 */		\n\t"\
		"addq		$0x20,%%rax							\n\t"\
		"addq		$0x20,%%rbx							\n\t"\
		"addq		$0x20,%%rcx							\n\t"\
		"addq		$0x20,%%rdx							\n\t"\
		"movq		%[__cc0],%%rdi						\n\t"\
		"movaps		0x100(%%rax),%%xmm4		\n\t		movaps		0x140(%%rax),%%xmm12		\n\t"\
		"movaps		0x180(%%rax),%%xmm6		\n\t		movaps		0x1c0(%%rax),%%xmm14		\n\t"\
		"movaps		0x110(%%rax),%%xmm5		\n\t		movaps		0x150(%%rax),%%xmm13		\n\t"\
		"movaps		0x190(%%rax),%%xmm7		\n\t		movaps		0x1d0(%%rax),%%xmm15		\n\t"\
		"movaps		0x100(%%rax),%%xmm0		\n\t		movaps		0x140(%%rax),%%xmm8		\n\t"\
		"movaps		0x180(%%rax),%%xmm2		\n\t		movaps		0x1c0(%%rax),%%xmm10		\n\t"\
		"movaps		0x110(%%rax),%%xmm1		\n\t		movaps		0x150(%%rax),%%xmm9		\n\t"\
		"movaps		0x190(%%rax),%%xmm3		\n\t		movaps		0x1d0(%%rax),%%xmm11		\n\t"\
		"mulpd		    (%%rdi),%%xmm4		\n\t		mulpd		0x10(%%rdi),%%xmm12		\n\t"\
		"mulpd		0x10(%%rdi),%%xmm6		\n\t		mulpd		    (%%rdi),%%xmm14		\n\t"\
		"mulpd		0x10(%%rdi),%%xmm1		\n\t		mulpd		    (%%rdi),%%xmm9		\n\t"\
		"mulpd		    (%%rdi),%%xmm3		\n\t		mulpd		0x10(%%rdi),%%xmm11		\n\t"\
		"mulpd		    (%%rdi),%%xmm5		\n\t		mulpd		0x10(%%rdi),%%xmm13		\n\t"\
		"mulpd		0x10(%%rdi),%%xmm7		\n\t		mulpd		    (%%rdi),%%xmm15		\n\t"\
		"mulpd		0x10(%%rdi),%%xmm0		\n\t		mulpd		    (%%rdi),%%xmm8		\n\t"\
		"mulpd		(%%rdi),%%xmm2			\n\t		mulpd		0x10(%%rdi),%%xmm10		\n\t"\
		"subpd		%%xmm1,%%xmm4			\n\t		subpd		%%xmm9,%%xmm12			\n\t"\
		"subpd		%%xmm3,%%xmm6			\n\t		subpd		%%xmm11,%%xmm14			\n\t"\
		"addpd		%%xmm0,%%xmm5			\n\t		addpd		%%xmm8,%%xmm13			\n\t"\
		"addpd		%%xmm2,%%xmm7			\n\t		addpd		%%xmm10,%%xmm15			\n\t"\
		"subpd		%%xmm6,%%xmm4			\n\t		subpd		%%xmm14,%%xmm12			\n\t"\
		"subpd		%%xmm7,%%xmm5			\n\t		subpd		%%xmm15,%%xmm13			\n\t"\
		"addpd		%%xmm6,%%xmm6			\n\t		addpd		%%xmm14,%%xmm14			\n\t"\
		"addpd		%%xmm7,%%xmm7			\n\t		addpd		%%xmm15,%%xmm15			\n\t"\
		"addpd		%%xmm4,%%xmm6			\n\t		addpd		%%xmm12,%%xmm14			\n\t"\
		"addpd		%%xmm5,%%xmm7			\n\t		addpd		%%xmm13,%%xmm15			\n\t"\
		"movaps		0x080(%%rax),%%xmm2		\n\t		movaps		0x0c0(%%rax),%%xmm10		\n\t"\
		"movaps		0x090(%%rax),%%xmm3		\n\t		movaps		0x0d0(%%rax),%%xmm11		\n\t"\
		"subpd		0x090(%%rax),%%xmm2		\n\t		addpd		0x0d0(%%rax),%%xmm10		\n\t"\
		"addpd		0x080(%%rax),%%xmm3		\n\t		subpd		0x0c0(%%rax),%%xmm11		\n\t"\
		"mulpd		(%%rsi),%%xmm2			\n\t		mulpd		(%%rsi),%%xmm10			\n\t"\
		"mulpd		(%%rsi),%%xmm3			\n\t		mulpd		(%%rsi),%%xmm11			\n\t"\
		"movaps		     (%%rax),%%xmm0		\n\t		movaps		0x040(%%rax),%%xmm8		\n\t"\
		"movaps		0x010(%%rax),%%xmm1		\n\t		movaps		0x050(%%rax),%%xmm9		\n\t"\
		"subpd		%%xmm2,%%xmm0			\n\t		subpd		%%xmm10,%%xmm8			\n\t"\
		"subpd		%%xmm3,%%xmm1			\n\t		subpd		%%xmm11,%%xmm9			\n\t"\
		"addpd		     (%%rax),%%xmm2		\n\t		addpd		0x040(%%rax),%%xmm10		\n\t"\
		"addpd		0x010(%%rax),%%xmm3		\n\t		addpd		0x050(%%rax),%%xmm11		\n\t"\
		"subpd		%%xmm6,%%xmm2			\n\t		subpd		%%xmm12,%%xmm8			\n\t"\
		"subpd		%%xmm7,%%xmm3			\n\t		subpd		%%xmm13,%%xmm9			\n\t"\
		"addpd		%%xmm6,%%xmm6			\n\t		addpd		%%xmm12,%%xmm12			\n\t"\
		"addpd		%%xmm7,%%xmm7			\n\t		addpd		%%xmm13,%%xmm13			\n\t"\
		"movaps		%%xmm2,    (%%rcx)		\n\t		movaps		%%xmm8,0x40(%%rcx)		\n\t"\
		"movaps		%%xmm3,0x10(%%rcx)		\n\t		movaps		%%xmm9,0x50(%%rcx)		\n\t"\
		"addpd		%%xmm2,%%xmm6			\n\t		addpd		%%xmm8,%%xmm12			\n\t"\
		"addpd		%%xmm3,%%xmm7			\n\t		addpd		%%xmm9,%%xmm13			\n\t"\
		"movaps		%%xmm6,    (%%rax)		\n\t		movaps		%%xmm12,0x40(%%rax)		\n\t"\
		"movaps		%%xmm7,0x10(%%rax)		\n\t		movaps		%%xmm13,0x50(%%rax)		\n\t"\
		"subpd		%%xmm5,%%xmm0			\n\t		subpd		%%xmm15,%%xmm10			\n\t"\
		"subpd		%%xmm4,%%xmm1			\n\t		subpd		%%xmm14,%%xmm11			\n\t"\
		"addpd		%%xmm5,%%xmm5			\n\t		addpd		%%xmm15,%%xmm15			\n\t"\
		"addpd		%%xmm4,%%xmm4			\n\t		addpd		%%xmm14,%%xmm14			\n\t"\
		"movaps		%%xmm0,    (%%rbx)		\n\t		movaps		%%xmm10,0x40(%%rbx)		\n\t"\
		"movaps		%%xmm1,0x10(%%rdx)		\n\t		movaps		%%xmm11,0x50(%%rdx)		\n\t"\
		"addpd		%%xmm0,		%%xmm5		\n\t		addpd		%%xmm10,%%xmm15			\n\t"\
		"addpd		%%xmm1,		%%xmm4		\n\t		addpd		%%xmm11,%%xmm14			\n\t"\
		"movaps		%%xmm5,    (%%rdx)		\n\t		movaps		%%xmm15,0x40(%%rdx)		\n\t"\
		"movaps		%%xmm4,0x10(%%rbx)		\n\t		movaps		%%xmm14,0x50(%%rbx)		\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__add1] "m" (Xadd1)\
		 ,[__r1] "m" (Xr1)\
		 ,[__r9] "m" (Xr9)\
		 ,[__r17] "m" (Xr17)\
		 ,[__r25] "m" (Xr25)\
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__cc0] "m" (Xcc0)\
		 ,[__c1] "m" (Xc1)\
		 ,[__c2] "m" (Xc2)\
		 ,[__c3] "m" (Xc3)\
		 ,[__c4] "m" (Xc4)\
		 ,[__c5] "m" (Xc5)\
		 ,[__c6] "m" (Xc6)\
		 ,[__c7] "m" (Xc7)\
		 ,[__c8] "m" (Xc8)\
		 ,[__c9] "m" (Xc9)\
		 ,[__c10] "m" (Xc10)\
		 ,[__c11] "m" (Xc11)\
		 ,[__c12] "m" (Xc12)\
		 ,[__c13] "m" (Xc13)\
		 ,[__c14] "m" (Xc14)\
		 ,[__c15] "m" (Xc15)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r10","r11","r12","r13","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"		/* Clobbered registers */\
	);\
	}

	#define SSE2_RADIX16_WRAPPER_DIT(Xadd0,Xadd1,Xr1,Xr9,Xr17,Xr25,Xisrt2,Xcc0,Xc1,Xc2,Xc3,Xc4,Xc5,Xc6,Xc7,Xc8,Xc9,Xc10,Xc11,Xc12,Xc13,Xc14,Xc15)\
	{\
	__asm__ volatile (\
		"/*...Block 1: */						\n\t"\
		"movq	%[__r1],%%rax					\n\t"\
		"movq	%%rax,%%rbx						\n\t"\
		"movq	%%rax,%%rcx						\n\t"\
		"movq	%%rax,%%rdx						\n\t"\
		"addq	$0x100,%%rbx					\n\t"\
		"addq	$0x080,%%rcx					\n\t"\
		"addq	$0x180,%%rdx					\n\t		/*...Block 3: */\n\t"\
		"/* SSE2_RADIX4_DIT_IN_PLACE(): */		\n\t		/* SSE2_RADIX4_DIT_IN_PLACE(): */		\n\t"\
		"movaps	    (%%rax),%%xmm0				\n\t		movaps	0x20(%%rax),%%xmm8 				\n\t"\
		"movaps	0x10(%%rax),%%xmm1				\n\t		movaps	0x30(%%rax),%%xmm9 				\n\t"\
		"movaps	    (%%rax),%%xmm2				\n\t		movaps	0x20(%%rax),%%xmm10				\n\t"\
		"movaps	0x10(%%rax),%%xmm3				\n\t		movaps	0x30(%%rax),%%xmm11				\n\t"\
		"addpd	    (%%rbx),%%xmm0				\n\t		addpd	0x20(%%rbx),%%xmm8 				\n\t"\
		"addpd	0x10(%%rbx),%%xmm1				\n\t		addpd	0x30(%%rbx),%%xmm9 				\n\t"\
		"subpd	    (%%rbx),%%xmm2				\n\t		subpd	0x20(%%rbx),%%xmm10				\n\t"\
		"subpd	0x10(%%rbx),%%xmm3				\n\t		subpd	0x30(%%rbx),%%xmm11				\n\t"\
		"movaps	    (%%rcx),%%xmm4				\n\t		movaps	0x20(%%rcx),%%xmm12				\n\t"\
		"movaps	0x10(%%rcx),%%xmm5				\n\t		movaps	0x30(%%rcx),%%xmm13				\n\t"\
		"movaps	    (%%rcx),%%xmm6				\n\t		movaps	0x20(%%rcx),%%xmm14				\n\t"\
		"movaps	0x10(%%rcx),%%xmm7				\n\t		movaps	0x30(%%rcx),%%xmm15				\n\t"\
		"addpd	    (%%rdx),%%xmm4				\n\t		addpd	0x20(%%rdx),%%xmm12				\n\t"\
		"addpd	0x10(%%rdx),%%xmm5				\n\t		addpd	0x30(%%rdx),%%xmm13				\n\t"\
		"subpd	    (%%rdx),%%xmm6				\n\t		subpd	0x20(%%rdx),%%xmm14				\n\t"\
		"subpd	0x10(%%rdx),%%xmm7				\n\t		subpd	0x30(%%rdx),%%xmm15				\n\t"\
		"subpd	%%xmm4,%%xmm0					\n\t		subpd	%%xmm12,%%xmm8 					\n\t"\
		"subpd	%%xmm5,%%xmm1					\n\t		subpd	%%xmm13,%%xmm9 					\n\t"\
		"movaps	%%xmm0,    (%%rbx)				\n\t		movaps	%%xmm8 ,0x20(%%rbx)				\n\t"\
		"movaps	%%xmm1,0x10(%%rbx)				\n\t		movaps	%%xmm9 ,0x30(%%rbx)				\n\t"\
		"addpd	%%xmm4,%%xmm4					\n\t		addpd	%%xmm12,%%xmm12					\n\t"\
		"addpd	%%xmm5,%%xmm5					\n\t		addpd	%%xmm13,%%xmm13					\n\t"\
		"addpd	%%xmm0,%%xmm4					\n\t		addpd	%%xmm8 ,%%xmm12					\n\t"\
		"addpd	%%xmm1,%%xmm5					\n\t		addpd	%%xmm9 ,%%xmm13					\n\t"\
		"movaps	%%xmm4,    (%%rax)				\n\t		movaps	%%xmm12,0x20(%%rax)				\n\t"\
		"movaps	%%xmm5,0x10(%%rax)				\n\t		movaps	%%xmm13,0x30(%%rax)				\n\t"\
		"subpd	%%xmm7,%%xmm2					\n\t		subpd	%%xmm15,%%xmm10					\n\t"\
		"subpd	%%xmm6,%%xmm3					\n\t		subpd	%%xmm14,%%xmm11					\n\t"\
		"movaps	%%xmm2,    (%%rdx)				\n\t		movaps	%%xmm10,0x20(%%rdx)				\n\t"\
		"movaps	%%xmm3,0x10(%%rcx)				\n\t		movaps	%%xmm11,0x30(%%rcx)				\n\t"\
		"addpd	%%xmm7,%%xmm7					\n\t		addpd	%%xmm15,%%xmm15					\n\t"\
		"addpd	%%xmm6,%%xmm6					\n\t		addpd	%%xmm14,%%xmm14					\n\t"\
		"addpd	%%xmm2,%%xmm7					\n\t		addpd	%%xmm10,%%xmm15					\n\t"\
		"addpd	%%xmm3,%%xmm6					\n\t		addpd	%%xmm11,%%xmm14					\n\t"\
		"movaps	%%xmm7,    (%%rcx)				\n\t		movaps	%%xmm15,0x20(%%rcx)				\n\t"\
		"movaps	%%xmm6,0x10(%%rdx)				\n\t		movaps	%%xmm14,0x30(%%rdx)				\n\t"\
		"/*...Block 2: */						\n\t		/*...Block 4: */						\n\t"\
		"addq	$0x040,%%rax					\n\t"\
		"addq	$0x040,%%rbx					\n\t"\
		"addq	$0x040,%%rcx					\n\t"\
		"addq	$0x040,%%rdx					\n\t"\
		"/* SSE2_RADIX4_DIT_IN_PLACE(): */		\n\t		/* SSE2_RADIX4_DIT_IN_PLACE(): */		\n\t"\
		"movaps	    (%%rax),%%xmm0				\n\t		movaps	0x20(%%rax),%%xmm8 				\n\t"\
		"movaps	0x10(%%rax),%%xmm1				\n\t		movaps	0x30(%%rax),%%xmm9 				\n\t"\
		"movaps	    (%%rax),%%xmm2				\n\t		movaps	0x20(%%rax),%%xmm10				\n\t"\
		"movaps	0x10(%%rax),%%xmm3				\n\t		movaps	0x30(%%rax),%%xmm11				\n\t"\
		"addpd	    (%%rbx),%%xmm0				\n\t		addpd	0x20(%%rbx),%%xmm8 				\n\t"\
		"addpd	0x10(%%rbx),%%xmm1				\n\t		addpd	0x30(%%rbx),%%xmm9 				\n\t"\
		"subpd	    (%%rbx),%%xmm2				\n\t		subpd	0x20(%%rbx),%%xmm10				\n\t"\
		"subpd	0x10(%%rbx),%%xmm3				\n\t		subpd	0x30(%%rbx),%%xmm11				\n\t"\
		"movaps	    (%%rcx),%%xmm4				\n\t		movaps	0x20(%%rcx),%%xmm12				\n\t"\
		"movaps	0x10(%%rcx),%%xmm5				\n\t		movaps	0x30(%%rcx),%%xmm13				\n\t"\
		"movaps	    (%%rcx),%%xmm6				\n\t		movaps	0x20(%%rcx),%%xmm14				\n\t"\
		"movaps	0x10(%%rcx),%%xmm7				\n\t		movaps	0x30(%%rcx),%%xmm15				\n\t"\
		"addpd	    (%%rdx),%%xmm4				\n\t		addpd	0x20(%%rdx),%%xmm12				\n\t"\
		"addpd	0x10(%%rdx),%%xmm5				\n\t		addpd	0x30(%%rdx),%%xmm13				\n\t"\
		"subpd	    (%%rdx),%%xmm6				\n\t		subpd	0x20(%%rdx),%%xmm14				\n\t"\
		"subpd	0x10(%%rdx),%%xmm7				\n\t		subpd	0x30(%%rdx),%%xmm15				\n\t"\
		"subpd	%%xmm4,%%xmm0					\n\t		subpd	%%xmm12,%%xmm8 					\n\t"\
		"subpd	%%xmm5,%%xmm1					\n\t		subpd	%%xmm13,%%xmm9 					\n\t"\
		"movaps	%%xmm0,    (%%rbx)				\n\t		movaps	%%xmm8 ,0x20(%%rbx)				\n\t"\
		"movaps	%%xmm1,0x10(%%rbx)				\n\t		movaps	%%xmm9 ,0x30(%%rbx)				\n\t"\
		"addpd	%%xmm4,%%xmm4					\n\t		addpd	%%xmm12,%%xmm12					\n\t"\
		"addpd	%%xmm5,%%xmm5					\n\t		addpd	%%xmm13,%%xmm13					\n\t"\
		"addpd	%%xmm0,%%xmm4					\n\t		addpd	%%xmm8 ,%%xmm12					\n\t"\
		"addpd	%%xmm1,%%xmm5					\n\t		addpd	%%xmm9 ,%%xmm13					\n\t"\
		"movaps	%%xmm4,    (%%rax)				\n\t		movaps	%%xmm12,0x20(%%rax)				\n\t"\
		"movaps	%%xmm5,0x10(%%rax)				\n\t		movaps	%%xmm13,0x30(%%rax)				\n\t"\
		"subpd	%%xmm7,%%xmm2					\n\t		subpd	%%xmm15,%%xmm10					\n\t"\
		"subpd	%%xmm6,%%xmm3					\n\t		subpd	%%xmm14,%%xmm11					\n\t"\
		"movaps	%%xmm2,    (%%rdx)				\n\t		movaps	%%xmm10,0x20(%%rdx)				\n\t"\
		"movaps	%%xmm3,0x10(%%rcx)				\n\t		movaps	%%xmm11,0x30(%%rcx)				\n\t"\
		"addpd	%%xmm7,%%xmm7					\n\t		addpd	%%xmm15,%%xmm15					\n\t"\
		"addpd	%%xmm6,%%xmm6					\n\t		addpd	%%xmm14,%%xmm14					\n\t"\
		"addpd	%%xmm2,%%xmm7					\n\t		addpd	%%xmm10,%%xmm15					\n\t"\
		"addpd	%%xmm3,%%xmm6					\n\t		addpd	%%xmm11,%%xmm14					\n\t"\
		"movaps	%%xmm7,    (%%rcx)				\n\t		movaps	%%xmm15,0x20(%%rcx)				\n\t"\
		"movaps	%%xmm6,0x10(%%rdx)				\n\t		movaps	%%xmm14,0x30(%%rdx)				\n\t"\
		"/***************************************************************************************************/\n\t"\
		"/*..and now do four more radix-4 transforms, including the internal and external twiddle factors.  */\n\t"\
		"/*  Write even-index 16-byte output pairs to a(j1), odd-index to a(j2), unpack same as on inputs.  */\n\t"\
		"/*  We do the last 2 radix-4 blocks first, to make the unpack-interleaving of outputs smoother.    */\n\t"\
		"/***************************************************************************************************/\n\t"\
		"/* Main-array addresses still in add0,1, no need to re-init: */									  \n\t"\
		"/*...Block 3: t3,11,19,27 -> r9,13,11,15: */	/*...Block 4: t7,15,23,31 -> r25,29,27,31: */	\n\t"\
		"movq		%[__r9],%%rax				\n\t	"/*movq		%[__r25],%%rax	Instead incr rcol rax offsets by +0x100 */\
		"movq		%[__isrt2],%%rbx			\n\t"\
		"movq		%[__cc0],%%rcx				\n\t"\
		"movaps		0x020(%%rax),%%xmm4			\n\t		movaps		0x120(%%rax),%%xmm12		\n\t"\
		"movaps		0x060(%%rax),%%xmm0			\n\t		movaps		0x160(%%rax),%%xmm8 		\n\t"\
		"movaps		0x030(%%rax),%%xmm5			\n\t		movaps		0x130(%%rax),%%xmm13		\n\t"\
		"movaps		0x070(%%rax),%%xmm1			\n\t		movaps		0x170(%%rax),%%xmm9 		\n\t"\
		"movaps		0x020(%%rax),%%xmm6			\n\t		movaps		0x120(%%rax),%%xmm14		\n\t"\
		"movaps		0x060(%%rax),%%xmm2			\n\t		movaps		0x160(%%rax),%%xmm10		\n\t"\
		"movaps		0x030(%%rax),%%xmm7			\n\t		movaps		0x130(%%rax),%%xmm15		\n\t"\
		"movaps		0x070(%%rax),%%xmm3			\n\t		movaps		0x170(%%rax),%%xmm11		\n\t"\
		"mulpd		    (%%rcx),%%xmm4			\n\t		mulpd		0x10(%%rcx),%%xmm12			\n\t"\
		"mulpd		0x10(%%rcx),%%xmm0			\n\t		mulpd		    (%%rcx),%%xmm8 			\n\t"\
		"mulpd		    (%%rcx),%%xmm5			\n\t		mulpd		0x10(%%rcx),%%xmm13			\n\t"\
		"mulpd		0x10(%%rcx),%%xmm1			\n\t		mulpd		    (%%rcx),%%xmm9 			\n\t"\
		"mulpd		0x10(%%rcx),%%xmm6			\n\t		mulpd		    (%%rcx),%%xmm14			\n\t"\
		"mulpd		    (%%rcx),%%xmm2			\n\t		mulpd		0x10(%%rcx),%%xmm10			\n\t"\
		"mulpd		0x10(%%rcx),%%xmm7			\n\t		mulpd		    (%%rcx),%%xmm15			\n\t"\
		"mulpd		    (%%rcx),%%xmm3			\n\t		mulpd		0x10(%%rcx),%%xmm11			\n\t"\
		"subpd		%%xmm6,%%xmm5				\n\t		subpd		%%xmm14,%%xmm13				\n\t"\
		"subpd		%%xmm2,%%xmm1				\n\t		subpd		%%xmm10,%%xmm9 				\n\t"\
		"addpd		%%xmm7,%%xmm4				\n\t		addpd		%%xmm15,%%xmm12				\n\t"\
		"addpd		%%xmm3,%%xmm0				\n\t		addpd		%%xmm11,%%xmm8 				\n\t"\
		"movaps		%%xmm5,%%xmm7				\n\t		movaps		%%xmm13,%%xmm15				\n\t"\
		"movaps		%%xmm4,%%xmm6				\n\t		movaps		%%xmm12,%%xmm14				\n\t"\
		"addpd		%%xmm0,%%xmm4				\n\t		addpd		%%xmm8 ,%%xmm12				\n\t"\
		"addpd		%%xmm1,%%xmm5				\n\t		addpd		%%xmm9 ,%%xmm13				\n\t"\
		"subpd		%%xmm0,%%xmm6				\n\t		subpd		%%xmm8 ,%%xmm14				\n\t"\
		"subpd		%%xmm1,%%xmm7				\n\t		subpd		%%xmm9 ,%%xmm15				\n\t"\
		"movaps		0x040(%%rax),%%xmm2			\n\t		movaps		0x140(%%rax),%%xmm10		\n\t"\
		"movaps		0x050(%%rax),%%xmm3			\n\t		movaps		0x150(%%rax),%%xmm11		\n\t"\
		"movaps		     (%%rax),%%xmm0			\n\t		movaps		0x100(%%rax),%%xmm8 		\n\t"\
		"movaps		0x010(%%rax),%%xmm1			\n\t		movaps		0x110(%%rax),%%xmm9 		\n\t"\
		"addpd		0x050(%%rax),%%xmm2			\n\t		subpd		0x150(%%rax),%%xmm10		\n\t"\
		"subpd		0x040(%%rax),%%xmm3			\n\t		addpd		0x140(%%rax),%%xmm11		\n\t"\
		"mulpd		(%%rbx),%%xmm2				\n\t		mulpd		(%%rbx),%%xmm10				\n\t"\
		"mulpd		(%%rbx),%%xmm3				\n\t		mulpd		(%%rbx),%%xmm11				\n\t"\
		"subpd		%%xmm2,%%xmm0				\n\t		subpd		%%xmm10,%%xmm8 				\n\t"\
		"subpd		%%xmm3,%%xmm1				\n\t		subpd		%%xmm11,%%xmm9 				\n\t"\
		"addpd		%%xmm2,%%xmm2				\n\t		addpd		%%xmm10,%%xmm10				\n\t"\
		"addpd		%%xmm3,%%xmm3				\n\t		addpd		%%xmm11,%%xmm11				\n\t"\
		"addpd		%%xmm0,%%xmm2				\n\t		addpd		%%xmm8 ,%%xmm10				\n\t"\
		"addpd		%%xmm1,%%xmm3				\n\t		addpd		%%xmm9 ,%%xmm11				\n\t"\
		"movq		%[__add1],%%rbx				\n\t"\
		"movq		%[__c1],%%rcx				\n\t		"/*movq		%[__c3],%%rcx	c3  = c1 + 8 */\
		"movq		%[__c9],%%rdx				\n\t		"/*movq		%[__c11],%%rdx	c11 = c9 + 8 */\
		"subpd		%%xmm4,%%xmm2				\n\t		subpd		%%xmm14,%%xmm8 				\n\t"\
		"subpd		%%xmm5,%%xmm3				\n\t		subpd		%%xmm15,%%xmm9 				\n\t"\
		"addpd		%%xmm4,%%xmm4				\n\t		addpd		%%xmm14,%%xmm14				\n\t"\
		"addpd		%%xmm5,%%xmm5				\n\t		addpd		%%xmm15,%%xmm15				\n\t"\
		"addpd		%%xmm2,%%xmm4				\n\t		addpd		%%xmm8 ,%%xmm14				\n\t"\
		"addpd		%%xmm3,%%xmm5				\n\t		addpd		%%xmm9 ,%%xmm15				\n\t"\
		"movaps		%%xmm2,     (%%rax)			\n\t		movaps		%%xmm8 ,0x100(%%rax)		\n\t"\
		"movaps		%%xmm3,0x010(%%rax)			\n\t		movaps		%%xmm9 ,0x110(%%rax)		\n\t"\
		"movaps		%%xmm4,%%xmm2				\n\t		movaps		%%xmm14,%%xmm8 				\n\t"\
		"movaps		%%xmm5,%%xmm3				\n\t		movaps		%%xmm15,%%xmm9 				\n\t"\
		"mulpd		    (%%rcx),%%xmm4			\n\t		mulpd		0x80(%%rcx),%%xmm14			\n\t"\
		"mulpd		    (%%rcx),%%xmm5			\n\t		mulpd		0x80(%%rcx),%%xmm15			\n\t"\
		"mulpd		0x10(%%rcx),%%xmm2			\n\t		mulpd		0x90(%%rcx),%%xmm8 			\n\t"\
		"mulpd		0x10(%%rcx),%%xmm3			\n\t		mulpd		0x90(%%rcx),%%xmm9 			\n\t"\
		"subpd		%%xmm2,%%xmm5				\n\t		subpd		%%xmm8 ,%%xmm15				\n\t"\
		"addpd		%%xmm3,%%xmm4				\n\t		addpd		%%xmm9 ,%%xmm14				\n\t"\
/* Tmps for later unpack into add1+0,1,4,5,8,9,c,d: */	/* Tmps for later unpack into add1+2,3,6,7,a,b,e,f: */\
		"movaps		%%xmm5,0x10(%%rbx)			\n\t		movaps		%%xmm15,0x30(%%rbx)			\n\t"\
		"movaps		%%xmm4,    (%%rbx)			\n\t		movaps		%%xmm14,0x20(%%rbx)			\n\t"\
		"movaps		     (%%rax),%%xmm4			\n\t		movaps		0x100(%%rax),%%xmm14		\n\t"\
		"movaps		0x010(%%rax),%%xmm5			\n\t		movaps		0x110(%%rax),%%xmm15		\n\t"\
		"movaps		%%xmm4,%%xmm2				\n\t		movaps		%%xmm14,%%xmm8 				\n\t"\
		"movaps		%%xmm5,%%xmm3				\n\t		movaps		%%xmm15,%%xmm9 				\n\t"\
		"mulpd		    (%%rdx),%%xmm4			\n\t		mulpd		0x80(%%rdx),%%xmm14			\n\t"\
		"mulpd		    (%%rdx),%%xmm5			\n\t		mulpd		0x80(%%rdx),%%xmm15			\n\t"\
		"mulpd		0x10(%%rdx),%%xmm2			\n\t		mulpd		0x90(%%rdx),%%xmm8 			\n\t"\
		"mulpd		0x10(%%rdx),%%xmm3			\n\t		mulpd		0x90(%%rdx),%%xmm9 			\n\t"\
		"subpd		%%xmm2,%%xmm5				\n\t		subpd		%%xmm8 ,%%xmm15				\n\t"\
		"addpd		%%xmm3,%%xmm4				\n\t		addpd		%%xmm9 ,%%xmm14				\n\t"\
		"movaps		%%xmm5,0x90(%%rbx)			\n\t		movaps		%%xmm15,0xb0(%%rbx)			\n\t"\
		"movaps		%%xmm4,0x80(%%rbx)			\n\t		movaps		%%xmm14,0xa0(%%rbx)			\n\t"\
		"movq		%[__c5],%%rcx				\n\t		"/*movq		%[__c7],%%rcx	c7  = c5  + 8 */\
		"movq		%[__c13],%%rdx				\n\t		"/*movq		%[__c15],%%rdx	c15 = c13 + 8 */\
		"subpd		%%xmm7,%%xmm0				\n\t		subpd		%%xmm13,%%xmm10				\n\t"\
		"subpd		%%xmm6,%%xmm1				\n\t		subpd		%%xmm12,%%xmm11				\n\t"\
		"addpd		%%xmm7,%%xmm7				\n\t		addpd		%%xmm13,%%xmm13				\n\t"\
		"addpd		%%xmm6,%%xmm6				\n\t		addpd		%%xmm12,%%xmm12				\n\t"\
		"addpd		%%xmm0,%%xmm7				\n\t		addpd		%%xmm10,%%xmm13				\n\t"\
		"addpd		%%xmm1,%%xmm6				\n\t		addpd		%%xmm11,%%xmm12				\n\t"\
		"movaps		%%xmm7,%%xmm4				\n\t		movaps		%%xmm13,%%xmm8 				\n\t"\
		"movaps		%%xmm1,%%xmm5				\n\t		movaps		%%xmm11,%%xmm9 				\n\t"\
		"mulpd		    (%%rcx),%%xmm7			\n\t		mulpd		0x80(%%rcx),%%xmm13			\n\t"\
		"mulpd		    (%%rcx),%%xmm1			\n\t		mulpd		0x80(%%rcx),%%xmm11			\n\t"\
		"mulpd		0x10(%%rcx),%%xmm4			\n\t		mulpd		0x90(%%rcx),%%xmm8 			\n\t"\
		"mulpd		0x10(%%rcx),%%xmm5			\n\t		mulpd		0x90(%%rcx),%%xmm9 			\n\t"\
		"subpd		%%xmm4,%%xmm1				\n\t		subpd		%%xmm8 ,%%xmm11				\n\t"\
		"addpd		%%xmm5,%%xmm7				\n\t		addpd		%%xmm9 ,%%xmm13				\n\t"\
		"movaps		%%xmm1,0x50(%%rbx)			\n\t		movaps		%%xmm11,0x70(%%rbx)			\n\t"\
		"movaps		%%xmm7,0x40(%%rbx)			\n\t		movaps		%%xmm13,0x60(%%rbx)			\n\t"\
		"movaps		%%xmm0,%%xmm4				\n\t		movaps		%%xmm10,%%xmm8 				\n\t"\
		"movaps		%%xmm6,%%xmm5				\n\t		movaps		%%xmm12,%%xmm9 				\n\t"\
		"mulpd		    (%%rdx),%%xmm0			\n\t		mulpd		0x80(%%rdx),%%xmm10			\n\t"\
		"mulpd		    (%%rdx),%%xmm6			\n\t		mulpd		0x80(%%rdx),%%xmm12			\n\t"\
		"mulpd		0x10(%%rdx),%%xmm4			\n\t		mulpd		0x90(%%rdx),%%xmm8 			\n\t"\
		"mulpd		0x10(%%rdx),%%xmm5			\n\t		mulpd		0x90(%%rdx),%%xmm9 			\n\t"\
		"subpd		%%xmm4,%%xmm6				\n\t		subpd		%%xmm8 ,%%xmm12				\n\t"\
		"addpd		%%xmm5,%%xmm0				\n\t		addpd		%%xmm9 ,%%xmm10				\n\t"\
		"movaps		%%xmm6,0xd0(%%rbx)			\n\t		movaps		%%xmm12,0xf0(%%rbx)			\n\t"\
		"movaps		%%xmm0,0xc0(%%rbx)			\n\t		movaps		%%xmm10,0xe0(%%rbx)			\n\t"\
		"/*...Block 1: t1,9,17,25 -> r1,5,3,7 */\n\t		/*...Block 3: t5,13,21,29 -> r17,21,19,23: */	\n\t"\
		"movq		%[__isrt2],%%rbx			\n\t"\
		"movq		%[__r1],%%rax				\n\t		movaps		0x120(%%rax),%%xmm12		\n\t"\
		"movaps		     (%%rax),%%xmm0			\n\t		movaps		0x130(%%rax),%%xmm13		\n\t"\
		"movaps		0x010(%%rax),%%xmm1			\n\t		movaps		0x160(%%rax),%%xmm8 		\n\t"\
		"movaps		0x040(%%rax),%%xmm2			\n\t		movaps		0x170(%%rax),%%xmm9 		\n\t"\
		"movaps		0x050(%%rax),%%xmm3			\n\t		movaps		(%%rbx),%%xmm10	/* isrt2 */	\n\t"\
		"subpd		0x040(%%rax),%%xmm0			\n\t		addpd		0x130(%%rax),%%xmm12		\n\t"\
		"subpd		0x050(%%rax),%%xmm1			\n\t		subpd		0x120(%%rax),%%xmm13		\n\t"\
		"addpd		     (%%rax),%%xmm2			\n\t		subpd		0x170(%%rax),%%xmm8 		\n\t"\
		"addpd		0x010(%%rax),%%xmm3			\n\t		addpd		0x160(%%rax),%%xmm9 		\n\t"\
		"movaps		0x020(%%rax),%%xmm4			\n\t		mulpd		%%xmm10,%%xmm12				\n\t"\
		"movaps		0x030(%%rax),%%xmm5			\n\t		mulpd		%%xmm10,%%xmm13				\n\t"\
		"movaps		0x060(%%rax),%%xmm6			\n\t		mulpd		%%xmm10,%%xmm8 				\n\t"\
		"movaps		0x070(%%rax),%%xmm7			\n\t		mulpd		%%xmm10,%%xmm9 				\n\t"\
		"subpd		0x060(%%rax),%%xmm4			\n\t		movaps		%%xmm12,%%xmm14				\n\t"\
		"subpd		0x070(%%rax),%%xmm5			\n\t		movaps		%%xmm13,%%xmm15				\n\t"\
		"addpd		0x020(%%rax),%%xmm6			\n\t		subpd		%%xmm8 ,%%xmm12				\n\t"\
		"addpd		0x030(%%rax),%%xmm7			\n\t		subpd		%%xmm9 ,%%xmm13				\n\t"\
		"movq		%[__add0],%%rbx				\n\t		addpd		%%xmm8 ,%%xmm14				\n\t"\
		"movq		%[__c8],%%rdx				\n\t		addpd		%%xmm9 ,%%xmm15				\n\t"\
		"addpd		%%xmm6,%%xmm2				\n\t		movaps		0x100(%%rax),%%xmm8 		\n\t"\
		"addpd		%%xmm7,%%xmm3				\n\t		movaps		0x110(%%rax),%%xmm9 		\n\t"\
		"movaps		%%xmm2,    (%%rbx)			\n\t		movaps		0x140(%%rax),%%xmm10		\n\t"\
		"movaps		%%xmm3,0x10(%%rbx)			\n\t		movaps		0x150(%%rax),%%xmm11		\n\t"\
		"addpd		%%xmm6,%%xmm6				\n\t		subpd		0x150(%%rax),%%xmm8 		\n\t"\
		"addpd		%%xmm7,%%xmm7				\n\t		subpd		0x140(%%rax),%%xmm9 		\n\t"\
		"subpd		%%xmm6,%%xmm2				\n\t		addpd		0x100(%%rax),%%xmm11		\n\t"\
		"subpd		%%xmm7,%%xmm3				\n\t		addpd		0x110(%%rax),%%xmm10		\n\t"\
		"movaps		%%xmm2,%%xmm6				\n\t		subpd		%%xmm12,%%xmm11				\n\t"\
		"movaps		%%xmm3,%%xmm7				\n\t		subpd		%%xmm13,%%xmm9 				\n\t"\
		"mulpd		    (%%rdx),%%xmm2			\n\t		addpd		%%xmm12,%%xmm12				\n\t"\
		"mulpd		    (%%rdx),%%xmm3			\n\t		addpd		%%xmm13,%%xmm13				\n\t"\
		"mulpd		0x10(%%rdx),%%xmm6			\n\t		addpd		%%xmm11,%%xmm12				\n\t"\
		"mulpd		0x10(%%rdx),%%xmm7			\n\t		addpd		%%xmm9 ,%%xmm13				\n\t"\
		"subpd		%%xmm6,%%xmm3				\n\t		movq %[__c10],%%rdx	\n\t"/*movq %[__c2],%%rcx	Use c2 = c10-2 since lcol needs rcx */\
		"addpd		%%xmm7,%%xmm2				\n\t		movaps		%%xmm11,0x100(%%rax)		\n\t"\
		"movq		%[__add1],%%rcx				\n\t		movaps		%%xmm9 ,0x110(%%rax)		\n\t"\
/*********************************************************************************************/\
/************ ALL STORES FROM HERE ON ARE FINAL **********************************************/\
/*********************************************************************************************/\
		"movaps		%%xmm3,%%xmm7				\n\t		movaps		%%xmm12,%%xmm11				\n\t"\
		"movaps		%%xmm2,%%xmm6				\n\t		movaps		%%xmm13,%%xmm9 				\n\t"\
		"unpckhpd	0x90(%%rcx),%%xmm7			\n\t		mulpd		-0x20(%%rdx),%%xmm12		\n\t"\
		"unpcklpd	0x90(%%rcx),%%xmm3			\n\t		mulpd		-0x20(%%rdx),%%xmm13		\n\t"\
		"movaps		%%xmm7,0x90(%%rcx)			\n\t		mulpd		-0x10(%%rdx),%%xmm11		\n\t"\
		"unpckhpd	0x80(%%rcx),%%xmm6			\n\t		mulpd		-0x10(%%rdx),%%xmm9 		\n\t"\
		"unpcklpd	0x80(%%rcx),%%xmm2			\n\t		subpd		%%xmm11,%%xmm13				\n\t"\
		"movaps		%%xmm6,0x80(%%rcx)			\n\t		addpd		%%xmm9 ,%%xmm12				\n\t"\
		"movaps		%%xmm3,0x90(%%rbx)			\n\t		movaps		%%xmm13,%%xmm11				\n\t"\
		"movaps		%%xmm2,0x80(%%rbx)			\n\t		movaps		%%xmm12,%%xmm9 				\n\t"\
/* DEBUG: Dump output-register contents back into corresponding local-store slots: */\
"movaps	%%xmm7,0x090(%%rax)	\n\t"\
"movaps	%%xmm6,0x080(%%rax)	\n\t"\
"movaps	%%xmm3,0x190(%%rax)	\n\t"\
"movaps	%%xmm2,0x180(%%rax)	\n\t"\
		"movaps		0x10(%%rbx),%%xmm3			\n\t		unpckhpd	0x30(%%rcx),%%xmm11			\n\t"\
		"movaps		    (%%rbx),%%xmm2			\n\t		unpcklpd	0x30(%%rcx),%%xmm13			\n\t"\
		"movaps		%%xmm3,%%xmm7				\n\t		movaps		%%xmm11,0x30(%%rcx)			\n\t"\
		"movaps		%%xmm2,%%xmm6				\n\t		unpckhpd	0x20(%%rcx),%%xmm9 			\n\t"\
		"unpckhpd	0x10(%%rcx),%%xmm7			\n\t		unpcklpd	0x20(%%rcx),%%xmm12			\n\t"\
		"unpcklpd	0x10(%%rcx),%%xmm3			\n\t		movaps		%%xmm9 ,0x20(%%rcx)			\n\t"\
		"movaps		%%xmm7,0x10(%%rcx)			\n\t		movaps		%%xmm13,0x30(%%rbx)			\n\t"\
		"unpckhpd	    (%%rcx),%%xmm6			\n\t		movaps		%%xmm12,0x20(%%rbx)			\n\t"\
													/* DEBUG: Dump output-register contents back into corresponding local-store slots: */\
													"movaps	%%xmm11,0x030(%%rax)	\n\t"\
													"movaps	%%xmm9 ,0x020(%%rax)	\n\t"\
													"movaps	%%xmm13,0x130(%%rax)	\n\t"\
													"movaps	%%xmm12,0x120(%%rax)	\n\t"\
		"unpcklpd	    (%%rcx),%%xmm2			\n\t		movaps		0x100(%%rax),%%xmm12		\n\t"\
		"movaps		%%xmm6,    (%%rcx)			\n\t		movaps		0x110(%%rax),%%xmm13		\n\t"\
		"movaps		%%xmm3,0x10(%%rbx)			\n\t		movaps		%%xmm12,%%xmm11				\n\t"\
		"movaps		%%xmm2,    (%%rbx)			\n\t		movaps		%%xmm13,%%xmm9 				\n\t"\
/* DEBUG: Dump output-register contents back into corresponding local-store slots: */\
"movaps	%%xmm7,0x010(%%rax)	\n\t"\
"movaps	%%xmm6,0x000(%%rax)	\n\t"\
"movaps	%%xmm3,0x110(%%rax)	\n\t"\
"movaps	%%xmm2,0x100(%%rax)	\n\t"\
		"addpd		%%xmm5,%%xmm0				\n\t		mulpd		    (%%rdx),%%xmm12			\n\t"\
		"subpd		%%xmm4,%%xmm1				\n\t		mulpd		    (%%rdx),%%xmm13			\n\t"\
		"movaps		%%xmm0,%%xmm2				\n\t		mulpd		0x10(%%rdx),%%xmm11			\n\t"\
		"movaps		%%xmm1,%%xmm3				\n\t		mulpd		0x10(%%rdx),%%xmm9 			\n\t"\
		"movq %[__c12],%%rdx	/* c4 = c12-2 */\n\t		subpd		%%xmm11,%%xmm13				\n\t"\
		"addpd		%%xmm5,%%xmm5				\n\t		addpd		%%xmm9 ,%%xmm12				\n\t"\
		"addpd		%%xmm4,%%xmm4				\n\t		movaps		%%xmm13,%%xmm11				\n\t"\
		"movaps		%%xmm0,%%xmm6				\n\t		movaps		%%xmm12,%%xmm9 				\n\t"\
		"movaps		%%xmm1,%%xmm7				\n\t		unpckhpd	0xb0(%%rcx),%%xmm11			\n\t"\
		"mulpd		-0x20(%%rdx),%%xmm2			\n\t		unpcklpd	0xb0(%%rcx),%%xmm13			\n\t"\
		"mulpd		-0x20(%%rdx),%%xmm3			\n\t		movaps		%%xmm11,0xb0(%%rcx)			\n\t"\
		"mulpd		-0x10(%%rdx),%%xmm6			\n\t		unpckhpd	0xa0(%%rcx),%%xmm9 			\n\t"\
		"mulpd		-0x10(%%rdx),%%xmm7			\n\t		unpcklpd	0xa0(%%rcx),%%xmm12			\n\t"\
		"subpd		%%xmm6,%%xmm3				\n\t		movaps		%%xmm9 ,0xa0(%%rcx)			\n\t"\
		"addpd		%%xmm7,%%xmm2				\n\t		movaps		%%xmm13,0xb0(%%rbx)			\n\t"\
		"movaps		%%xmm3,%%xmm7				\n\t		movaps		%%xmm12,0xa0(%%rbx)			\n\t"\
													/* DEBUG: Dump output-register contents back into corresponding local-store slots: */\
													"movaps	%%xmm11,0x0b0(%%rax)	\n\t"\
													"movaps	%%xmm9 ,0x0a0(%%rax)	\n\t"\
													"movaps	%%xmm13,0x1b0(%%rax)	\n\t"\
													"movaps	%%xmm12,0x1a0(%%rax)	\n\t"\
		"movaps		%%xmm2,%%xmm6				\n\t		subpd		%%xmm15,%%xmm8 				\n\t"\
		"unpckhpd	0x50(%%rcx),%%xmm7			\n\t		subpd		%%xmm14,%%xmm10				\n\t"\
		"unpcklpd	0x50(%%rcx),%%xmm3			\n\t		addpd		%%xmm15,%%xmm15				\n\t"\
		"movaps		%%xmm7,0x50(%%rcx)			\n\t		addpd		%%xmm14,%%xmm14				\n\t"\
		"unpckhpd	0x40(%%rcx),%%xmm6			\n\t		addpd		%%xmm8 ,%%xmm15				\n\t"\
		"unpcklpd	0x40(%%rcx),%%xmm2			\n\t		addpd		%%xmm10,%%xmm14				\n\t"\
		"movaps		%%xmm6,0x40(%%rcx)			\n\t		movaps		%%xmm15,%%xmm12				\n\t"\
		"movaps		%%xmm3,0x50(%%rbx)			\n\t		movaps		%%xmm10,%%xmm13				\n\t"\
		"movaps		%%xmm2,0x40(%%rbx)			\n\t		mulpd		0x60(%%rdx),%%xmm15	/* c6  = c12 + 6 */\n\t"\
/* DEBUG: Dump output-register contents back into corresponding local-store slots: */\
"movaps	%%xmm7,0x050(%%rax)	\n\t"\
"movaps	%%xmm6,0x040(%%rax)	\n\t"\
"movaps	%%xmm3,0x150(%%rax)	\n\t"\
"movaps	%%xmm2,0x140(%%rax)	\n\t"\
		"subpd		%%xmm5,%%xmm0				\n\t		mulpd		0x60(%%rdx),%%xmm10			\n\t"\
		"addpd		%%xmm4,%%xmm1				\n\t		mulpd		0x70(%%rdx),%%xmm12			\n\t"\
		"movaps		%%xmm0,%%xmm6				\n\t		mulpd		0x70(%%rdx),%%xmm13			\n\t"\
		"movaps		%%xmm1,%%xmm7				\n\t		subpd		%%xmm12,%%xmm10				\n\t"\
		"mulpd		    (%%rdx),%%xmm0			\n\t		addpd		%%xmm13,%%xmm15				\n\t"\
		"mulpd		    (%%rdx),%%xmm1			\n\t		movaps		%%xmm10,%%xmm13				\n\t"\
		"mulpd		0x10(%%rdx),%%xmm6			\n\t		movaps		%%xmm15,%%xmm12				\n\t"\
		"mulpd		0x10(%%rdx),%%xmm7			\n\t		unpckhpd	0x70(%%rcx),%%xmm13			\n\t"\
		"subpd		%%xmm6,%%xmm1				\n\t		unpcklpd	0x70(%%rcx),%%xmm10			\n\t"\
		"addpd		%%xmm7,%%xmm0				\n\t		movaps		%%xmm13,0x70(%%rcx)			\n\t"\
		"movaps		%%xmm1,%%xmm7				\n\t		unpckhpd	0x60(%%rcx),%%xmm12			\n\t"\
		"movaps		%%xmm0,%%xmm6				\n\t		unpcklpd	0x60(%%rcx),%%xmm15			\n\t"\
		"unpckhpd	0xd0(%%rcx),%%xmm7			\n\t		movaps		%%xmm12,0x60(%%rcx)			\n\t"\
		"unpcklpd	0xd0(%%rcx),%%xmm1			\n\t		movaps		%%xmm10,0x70(%%rbx)			\n\t"\
		"movaps		%%xmm7,0xd0(%%rcx)			\n\t		movaps		%%xmm15,0x60(%%rbx)			\n\t"\
													/* DEBUG: Dump output-register contents back into corresponding local-store slots: */\
													"movaps	%%xmm13,0x070(%%rax)	\n\t"\
													"movaps	%%xmm12,0x060(%%rax)	\n\t"\
													"movaps	%%xmm10,0x170(%%rax)	\n\t"\
													"movaps	%%xmm15,0x160(%%rax)	\n\t"\
		"unpckhpd	0xc0(%%rcx),%%xmm6			\n\t		movaps		%%xmm8 ,%%xmm12				\n\t"\
		"unpcklpd	0xc0(%%rcx),%%xmm0			\n\t		movaps		%%xmm14,%%xmm13				\n\t"\
		"movaps		%%xmm6,0xc0(%%rcx)			\n\t		mulpd		0x80(%%rdx),%%xmm8 	/* c14 = c12 + 8 */\n\t"\
		"movaps		%%xmm1,0xd0(%%rbx)			\n\t		mulpd		0x80(%%rdx),%%xmm14			\n\t"\
		"movaps		%%xmm0,0xc0(%%rbx)			\n\t		mulpd		0x90(%%rdx),%%xmm12			\n\t"\
/* DEBUG: Dump output-register contents back into corresponding local-store slots: */\
"movaps	%%xmm7,0x0d0(%%rax)	\n\t"\
"movaps	%%xmm6,0x0c0(%%rax)	\n\t"\
"movaps	%%xmm1,0x1d0(%%rax)	\n\t"\
"movaps	%%xmm0,0x1c0(%%rax)	\n\t"\
		"													mulpd		0x90(%%rdx),%%xmm13			\n\t"\
		"													subpd		%%xmm12,%%xmm14				\n\t"\
		"													addpd		%%xmm13,%%xmm8 				\n\t"\
		"													movaps		%%xmm14,%%xmm13				\n\t"\
		"													movaps		%%xmm8 ,%%xmm12				\n\t"\
		"													unpckhpd	0xf0(%%rcx),%%xmm13			\n\t"\
		"													unpcklpd	0xf0(%%rcx),%%xmm14			\n\t"\
		"													movaps		%%xmm13,0xf0(%%rcx)			\n\t"\
		"													unpckhpd	0xe0(%%rcx),%%xmm12			\n\t"\
		"													unpcklpd	0xe0(%%rcx),%%xmm8 			\n\t"\
		"													movaps		%%xmm12,0xe0(%%rcx)			\n\t"\
		"													movaps		%%xmm14,0xf0(%%rbx)			\n\t"\
		"													movaps		%%xmm8 ,0xe0(%%rbx)			\n\t"\
													/* DEBUG: Dump output-register contents back into corresponding local-store slots: */\
													"movaps	%%xmm13,0x0f0(%%rax)	\n\t"\
													"movaps	%%xmm12,0x0e0(%%rax)	\n\t"\
													"movaps	%%xmm14,0x1f0(%%rax)	\n\t"\
													"movaps	%%xmm8 ,0x1e0(%%rax)	\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__add1] "m" (Xadd1)\
		 ,[__r1] "m" (Xr1)\
		 ,[__r9] "m" (Xr9)\
		 ,[__r17] "m" (Xr17)\
		 ,[__r25] "m" (Xr25)\
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__cc0] "m" (Xcc0)\
		 ,[__c1] "m" (Xc1)\
		 ,[__c2] "m" (Xc2)\
		 ,[__c3] "m" (Xc3)\
		 ,[__c4] "m" (Xc4)\
		 ,[__c5] "m" (Xc5)\
		 ,[__c6] "m" (Xc6)\
		 ,[__c7] "m" (Xc7)\
		 ,[__c8] "m" (Xc8)\
		 ,[__c9] "m" (Xc9)\
		 ,[__c10] "m" (Xc10)\
		 ,[__c11] "m" (Xc11)\
		 ,[__c12] "m" (Xc12)\
		 ,[__c13] "m" (Xc13)\
		 ,[__c14] "m" (Xc14)\
		 ,[__c15] "m" (Xc15)\
		: "cc","memory","rax","rbx","rcx","rdx","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"		/* Clobbered registers */\
	);\
	}

  #endif // USE_64BIT_ASM_STYLE

	// Debug dump-registers macro: Insert body into needed location above
	#define FOO(Xaddr)\
	{\
	__asm__ volatile (\
		"\n\t"\
		"movq	%[__addr],%%rax	/* DEBUG: Dump register contents into trig-data slots: */\n\t"\
		"movaps	%%xmm0 ,0x00(%%rax)	\n\t"\
		"movaps	%%xmm1 ,0x10(%%rax)	\n\t"\
		"movaps	%%xmm2 ,0x20(%%rax)	\n\t"\
		"movaps	%%xmm3 ,0x30(%%rax)	\n\t"\
		"movaps	%%xmm4 ,0x40(%%rax)	\n\t"\
		"movaps	%%xmm5 ,0x50(%%rax)	\n\t"\
		"movaps	%%xmm6 ,0x60(%%rax)	\n\t"\
		"movaps	%%xmm7 ,0x70(%%rax)	\n\t"\
		"movaps	%%xmm8 ,0x80(%%rax)	\n\t"\
		"movaps	%%xmm9 ,0x90(%%rax)	\n\t"\
		"movaps	%%xmm10,0xa0(%%rax)	\n\t"\
		"movaps	%%xmm11,0xb0(%%rax)	\n\t"\
		"movaps	%%xmm12,0xc0(%%rax)	\n\t"\
		"movaps	%%xmm13,0xd0(%%rax)	\n\t"\
		"movaps	%%xmm14,0xe0(%%rax)	\n\t"\
		"movaps	%%xmm15,0xf0(%%rax)	\n\t"\
		:					/* outputs: none */\
		: [__addr] "m" (Xaddr)	/* All inputs from memory addresses here */\
		: "cc","memory","rax","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
	);\
	}

#endif	// SSE2 or AVX?

#endif	/* radix16_wrapper_square_gcc_h_included */

