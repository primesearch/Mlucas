/*******************************************************************************
*                                                                             *
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
#ifndef radix16_dyadic_square_gcc_h_included
#define radix16_dyadic_square_gcc_h_included

#ifdef USE_AVX

	/* In AVX mode, the input data are arranged in memory like so, where we view things in 32-byte chunks:

		&a[j1]	a0.re,a1.re,a2.re,a3.re	&a[j1+32]	b0.re,b1.re,b2.re,b3.re	 &a[j1+64]	c0.re,c1.re,c2.re,c3.re	 &a[j1+96]	d0.re,d1.re,d2.re,d3.re
		+0x020	a0.im,a1.im,a2.im,a3.im		+0x020	b0.im,b1.im,b2.im,b3.im		+0x020	c0.im,c1.im,c2.im,c3.im		+0x020	d0.im,d1.im,d2.im,d3.im
		+0x040	a4.re,a5.re,a6.re,a7.re		+0x040	b4.re,b5.re,b6.re,b7.re		+0x040	c4.re,c5.re,c6.re,c7.re		+0x040	d4.re,d5.re,d6.re,d7.re
		+0x060	a4.im,a5.im,a6.im,a7.im		+0x060	b4.im,b5.im,b6.im,b7.im		+0x060	c4.im,c5.im,c6.im,c7.im		+0x060	d4.im,d5.im,d6.im,d7.im
		+0x080	a8.re,a9.re,aA.re,aB.re		+0x080	b8.re,b9.re,bA.re,bB.re		+0x080	c8.re,c9.re,cA.re,cB.re		+0x080	d8.re,d9.re,dA.re,dB.re
		+0x0a0	a8.im,a9.im,aA.im,aB.im		+0x0a0	b8.im,b9.im,bA.im,bB.im		+0x0a0	c8.im,c9.im,cA.im,cB.im		+0x0a0	d8.im,d9.im,dA.im,dB.im
		+0x0c0	aC.re,aD.re,aE.re,aF.re		+0x0c0	bC.re,bD.re,bE.re,bF.re		+0x0c0	cC.re,cD.re,cE.re,cF.re		+0x0c0	dC.re,dD.re,dE.re,dF.re
		+0x0e0	aC.im,aD.im,aE.im,aF.im		+0x0e0	bC.im,bD.im,bE.im,bF.im		+0x0e0	cC.im,cD.im,cE.im,cF.im		+0x0e0	dC.im,dD.im,dE.im,dF.im

	Thus, we have 16 input complex data quartets, whose subelements need to be properly interleaved prior to the radix-16 complex DFT.
	This is analogous to the SSE2 case, but there our 16 vector-complex inputs had Re/Im parts 128 bits wide, so we interleaved
	2 such at a time, using one for the ensuing radix-4 sub-DFT computation, and dumping the other to the part of local memory
	where it will end up after also being radix-4-DFT-butterflied with 3 other such data.

	In the AVX case, we have 256-bit 4-double-containing Re/Im data, and need to read in 4 such at a time, do the 4-way
	interleaving described below, and write 2 to local memory for later use, using the other 2 for immediate computation.
	In other words, we collect 2 pairs at a time of unpckhpd/unpckhld events from the SSE2 code (both from the same radix-4 DFT)
	and replace them with a single 4-way interleaving of the type detailed below.

	We can break this into a bunch of smaller interleave steps, in each of which have 4 register quartets-of-doubles
	a0-3, b0-3, c0-3, d0-3, which we need to interleave like so (= a 4x4 matrix transposition, if one considers each
	input register to hold a row of the matrix):

		[a0,a1,a2,a3]         |\     [a0,b0,c0,d0]
		[b0,b1,b2,b3]    -----  \    [a1,b1,c1,d1]
		[c0,c1,c2,c3]    -----  /    [a2,b2,c2,d2]
		[d0,d1,d2,d3]         |/     [a3,b3,c3,d3]

	George says he uses this sequence for this transposition operation:

		vmovaps	ymm1, [srcreg]				;; R1											a0,a1,a2,a3
		vmovaps	ymm7, [srcreg+d1]			;; R2											b0,b1,b2,b3
		vshufpd	ymm0, ymm1, ymm7, 15		;; Shuffle R1 and R2 to create R1/R2 hi			ymm0:a1,b1,a3,b3
		vshufpd	ymm1, ymm1, ymm7, 0			;; Shuffle R1 and R2 to create R1/R2 low		ymm1:a0,b0,a2,b2

		vmovaps	ymm2, [srcreg+d2]			;; R3											c0,c1,c2,c3
		vmovaps	ymm7, [srcreg+d2+d1]		;; R4											d0,d1,d2,d3
		vshufpd	ymm3, ymm2, ymm7, 15		;; Shuffle R3 and R4 to create R3/R4 hi			ymm3:c1,d1,c3,d3
		vshufpd	ymm2, ymm2, ymm7, 0			;; Shuffle R3 and R4 to create R3/R4 low		ymm2:c0,d0,c2,d2

		vperm2f128 ymm4, ymm0, ymm3, 32		;; Shuffle R1/R2 hi and R3/R4 hi (new R2)		ymm4:[a1,b1][c1,d1]	imm0:1=0,4:5=2 -> [in1.lo,in2.lo] = [ymm0.lo,ymm3.lo]
		vperm2f128 ymm0, ymm0, ymm3, 49		;; Shuffle R1/R2 hi and R3/R4 hi (new R4)		ymm0:[a3,b3][c3,d3]	imm0:1=1,4:5=3 -> [in1.hi,in2.hi] = [ymm0.hi,ymm3.hi]

		vperm2f128 ymm3, ymm1, ymm2, 32		;; Shuffle R1/R2 low and R3/R4 low (new R1)		ymm3:[a0,b0][c0,d0]	imm0:1=0,4:5=2 -> [in1.lo,in2.lo] = [ymm1.lo,ymm2.lo]
		vperm2f128 ymm1, ymm1, ymm2, 49		;; Shuffle R1/R2 low and R3/R4 low (new R3)		ymm1:[a2,b2][c2,d2]	imm0:1=1,4:5=3 -> [in1.hi,in2.hi] = [ymm1.hi,ymm2.hi]

	I want output-reg index to match subscripts of final a-d set stored in the register, so swap ymm indices 0134 -> 3201:

											// Data exit in ymm0-3; ymmX,Y are any 2 other registers
		vmovaps	ymm2, [srcreg]				// a0,a1,a2,a3
		vmovaps	ymmX, [srcreg+d1]			// b0,b1,b2,b3
		vshufpd	ymm3, ymm2, ymmX, 15		// ymm3:a1,b1,a3,b3
		vshufpd	ymm2, ymm2, ymmX, 0			// ymm2:a0,b0,a2,b2

		vmovaps	ymmY, [srcreg+d2]			// c0,c1,c2,c3
		vmovaps	ymmX, [srcreg+d2+d1]		// d0,d1,d2,d3
		vshufpd	ymm0, ymmY, ymmX, 15		// ymm0:c1,d1,c3,d3
		vshufpd	ymmY, ymmY, ymmX, 0			// ymmY:c0,d0,c2,d2

		vperm2f128 ymm1, ymm3, ymm0, 32		// ymm1:[a1,b1][c1,d1]
		vperm2f128 ymm3, ymm3, ymm0, 49		// ymm3:[a3,b3][c3,d3]

		vperm2f128 ymm0, ymm2, ymmY, 32		// ymm0:[a0,b0][c0,d0]
		vperm2f128 ymm2, ymm2, ymmY, 49		// ymm2:[a2,b2][c2,d2]
	*/
	#define SSE2_RADIX16_DIF_DYADIC_DIT(Xadd0,Xadd1,Xr1,Xisrt2,Xpfetch_dist)\
	{\
	__asm__ volatile (\
	"/*************************************************************/\n\t"\
	"/* SSE2_RADIX16_WRAPPER_DIF, 1st set of inputs:              */\n\t"\
	"/*************************************************************/\n\t"\
		"movq	%[__r1] ,%%rcx	\n\t"\
		"movq	%[__add0],%%rax\n\t"\
		"movq	%[__add1],%%rbx\n\t"\
	"movslq	%[__pfetch_dist],%%r15	\n\t"\
	"leaq	(%%rax,%%r15,8),%%r14	\n\t"	/* Block 1 [base-address + data-fetch-ahead index] */\
	"leaq	(%%rbx,%%r15,8),%%r15	\n\t"	/* Block 2 [base-address + data-fetch-ahead index] */\
	"prefetcht1	(%%r14)\n\t"\
	"/**** Start with 4-way interleaving: ****/\n\t"\
	"/* a[j+p0]: Inputs from [add0,add0+0x100,add1,add1+0x100]+0x0. Outputs into r1 +0/1, 8/9, 16/17, 24/25: */	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps	     (%%rax),%%ymm1						\n\t		vmovaps	0x020(%%rax),%%ymm3								\n\t"\
		"vmovaps	0x100(%%rax),%%ymm5						\n\t		vmovaps	0x120(%%rax),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps	     (%%rbx),%%ymm6						\n\t		vmovaps	0x020(%%rbx),%%ymm14							\n\t"\
		"vmovaps	0x100(%%rbx),%%ymm5						\n\t		vmovaps	0x120(%%rbx),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm6 ,%%ymm0				\n\t		vshufpd	$15 ,%%ymm13,%%ymm14 ,%%ymm2					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm6 ,%%ymm6				\n\t		vshufpd	$0  ,%%ymm13,%%ymm14 ,%%ymm14					\n\t"\
		"vperm2f128 $32 ,%%ymm0 ,%%ymm7 ,%%ymm5				\n\t		vperm2f128 $32 ,%%ymm2 ,%%ymm15,%%ymm13		/* outB	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm0 ,%%ymm7 ,%%ymm7				\n\t		vperm2f128 $49 ,%%ymm2 ,%%ymm15,%%ymm15		/* outD	*/	\n\t"\
		"vperm2f128 $32 ,%%ymm6 ,%%ymm1 ,%%ymm0				\n\t		vperm2f128 $32 ,%%ymm14,%%ymm3 ,%%ymm2 		/* outA	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm6 ,%%ymm1 ,%%ymm1				\n\t		vperm2f128 $49 ,%%ymm14,%%ymm3 ,%%ymm3 		/* outC	*/	\n\t"\
		"vmovaps 	%%ymm5 ,0x200(%%rcx)					\n\t		vmovaps %%ymm13,0x220(%%rcx)				/* outB	*/	\n\t"\
		"vmovaps 	%%ymm7 ,0x300(%%rcx)					\n\t		vmovaps %%ymm15,0x320(%%rcx)				/* outD	*/	\n\t"\
		"vmovaps 	%%ymm0 ,     (%%rcx)					\n\t		vmovaps %%ymm2 ,0x020(%%rcx)				/* outA	*/	\n\t"\
		"vmovaps 	%%ymm1 ,0x100(%%rcx)					\n\t		vmovaps %%ymm3 ,0x120(%%rcx)				/* outC	*/	\n\t"\
	"/* a[j+p2]: Inputs from [add0,add0+0x100,add1,add1+0x100]+0x40. Outputs into r3 +0/1, 8/9, 16/17, 24/25: */		\n\t"\
		"addq	$0x40,%%rax	\n\t"\
		"addq	$0x40,%%rbx	\n\t"\
		"addq	$0x40,%%rcx	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps	     (%%rax),%%ymm1						\n\t		vmovaps	0x020(%%rax),%%ymm3								\n\t"\
		"vmovaps	0x100(%%rax),%%ymm5						\n\t		vmovaps	0x120(%%rax),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps	     (%%rbx),%%ymm6						\n\t		vmovaps	0x020(%%rbx),%%ymm14							\n\t"\
		"vmovaps	0x100(%%rbx),%%ymm5						\n\t		vmovaps	0x120(%%rbx),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm6 ,%%ymm0				\n\t		vshufpd	$15 ,%%ymm13,%%ymm14 ,%%ymm2					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm6 ,%%ymm6				\n\t		vshufpd	$0  ,%%ymm13,%%ymm14 ,%%ymm14					\n\t"\
		"vperm2f128 $32 ,%%ymm0 ,%%ymm7 ,%%ymm5				\n\t		vperm2f128 $32 ,%%ymm2 ,%%ymm15,%%ymm13		/* outB	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm0 ,%%ymm7 ,%%ymm7				\n\t		vperm2f128 $49 ,%%ymm2 ,%%ymm15,%%ymm15		/* outD	*/	\n\t"\
		"vperm2f128 $32 ,%%ymm6 ,%%ymm1 ,%%ymm0				\n\t		vperm2f128 $32 ,%%ymm14,%%ymm3 ,%%ymm2 		/* outA	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm6 ,%%ymm1 ,%%ymm1				\n\t		vperm2f128 $49 ,%%ymm14,%%ymm3 ,%%ymm3 		/* outC	*/	\n\t"\
		"vmovaps 	%%ymm5 ,0x200(%%rcx)					\n\t		vmovaps %%ymm13,0x220(%%rcx)				/* outB	*/	\n\t"\
		"vmovaps 	%%ymm7 ,0x300(%%rcx)					\n\t		vmovaps %%ymm15,0x320(%%rcx)				/* outD	*/	\n\t"\
		"vmovaps 	%%ymm0 ,     (%%rcx)					\n\t		vmovaps %%ymm2 ,0x020(%%rcx)				/* outA	*/	\n\t"\
		"vmovaps 	%%ymm1 ,0x100(%%rcx)					\n\t		vmovaps %%ymm3 ,0x120(%%rcx)				/* outC	*/	\n\t"\
	"/* a[j+p4]: Inputs from [add0,add0+0x100,add1,add1+0x100]+0x80. Outputs into r5 +0/1, 8/9, 16/17, 24/25: */		\n\t"\
		"addq	$0x40,%%rax	\n\t"\
		"addq	$0x40,%%rbx	\n\t"\
		"addq	$0x40,%%rcx	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps	     (%%rax),%%ymm1						\n\t		vmovaps	0x020(%%rax),%%ymm3								\n\t"\
		"vmovaps	0x100(%%rax),%%ymm5						\n\t		vmovaps	0x120(%%rax),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps	     (%%rbx),%%ymm6						\n\t		vmovaps	0x020(%%rbx),%%ymm14							\n\t"\
		"vmovaps	0x100(%%rbx),%%ymm5						\n\t		vmovaps	0x120(%%rbx),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm6 ,%%ymm0				\n\t		vshufpd	$15 ,%%ymm13,%%ymm14 ,%%ymm2					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm6 ,%%ymm6				\n\t		vshufpd	$0  ,%%ymm13,%%ymm14 ,%%ymm14					\n\t"\
		"vperm2f128 $32 ,%%ymm0 ,%%ymm7 ,%%ymm5				\n\t		vperm2f128 $32 ,%%ymm2 ,%%ymm15,%%ymm13		/* outB	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm0 ,%%ymm7 ,%%ymm7				\n\t		vperm2f128 $49 ,%%ymm2 ,%%ymm15,%%ymm15		/* outD	*/	\n\t"\
		"vperm2f128 $32 ,%%ymm6 ,%%ymm1 ,%%ymm0				\n\t		vperm2f128 $32 ,%%ymm14,%%ymm3 ,%%ymm2 		/* outA	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm6 ,%%ymm1 ,%%ymm1				\n\t		vperm2f128 $49 ,%%ymm14,%%ymm3 ,%%ymm3 		/* outC	*/	\n\t"\
		"vmovaps 	%%ymm5 ,0x200(%%rcx)					\n\t		vmovaps %%ymm13,0x220(%%rcx)				/* outB	*/	\n\t"\
		"vmovaps 	%%ymm7 ,0x300(%%rcx)					\n\t		vmovaps %%ymm15,0x320(%%rcx)				/* outD	*/	\n\t"\
		"vmovaps 	%%ymm0 ,     (%%rcx)					\n\t		vmovaps %%ymm2 ,0x020(%%rcx)				/* outA	*/	\n\t"\
		"vmovaps 	%%ymm1 ,0x100(%%rcx)					\n\t		vmovaps %%ymm3 ,0x120(%%rcx)				/* outC	*/	\n\t"\
	"/* a[j+p6]: Inputs from [add0,add0+0x100,add1,add1+0x100]+0xc0. Outputs into r7 +0/1, 8/9, 16/17, 24/25: */		\n\t"\
		"addq	$0x40,%%rax	\n\t"\
		"addq	$0x40,%%rbx	\n\t"\
		"addq	$0x40,%%rcx	\n\t"\
	"prefetcht1	(%%r15)\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps	     (%%rax),%%ymm1						\n\t		vmovaps	0x020(%%rax),%%ymm3								\n\t"\
		"vmovaps	0x100(%%rax),%%ymm5						\n\t		vmovaps	0x120(%%rax),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps	     (%%rbx),%%ymm6						\n\t		vmovaps	0x020(%%rbx),%%ymm14							\n\t"\
		"vmovaps	0x100(%%rbx),%%ymm5						\n\t		vmovaps	0x120(%%rbx),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm6 ,%%ymm0				\n\t		vshufpd	$15 ,%%ymm13,%%ymm14 ,%%ymm2					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm6 ,%%ymm6				\n\t		vshufpd	$0  ,%%ymm13,%%ymm14 ,%%ymm14					\n\t"\
		"vperm2f128 $32 ,%%ymm0 ,%%ymm7 ,%%ymm5				\n\t		vperm2f128 $32 ,%%ymm2 ,%%ymm15,%%ymm13		/* outB	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm0 ,%%ymm7 ,%%ymm7				\n\t		vperm2f128 $49 ,%%ymm2 ,%%ymm15,%%ymm15		/* outD	*/	\n\t"\
		"vperm2f128 $32 ,%%ymm6 ,%%ymm1 ,%%ymm0				\n\t		vperm2f128 $32 ,%%ymm14,%%ymm3 ,%%ymm2 		/* outA	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm6 ,%%ymm1 ,%%ymm1				\n\t		vperm2f128 $49 ,%%ymm14,%%ymm3 ,%%ymm3 		/* outC	*/	\n\t"\
		"vmovaps 	%%ymm5 ,0x200(%%rcx)					\n\t		vmovaps %%ymm13,0x220(%%rcx)				/* outB	*/	\n\t"\
		"vmovaps 	%%ymm7 ,0x300(%%rcx)					\n\t		vmovaps %%ymm15,0x320(%%rcx)				/* outD	*/	\n\t"\
		"vmovaps 	%%ymm0 ,     (%%rcx)					\n\t		vmovaps %%ymm2 ,0x020(%%rcx)				/* outA	*/	\n\t"\
		"vmovaps 	%%ymm1 ,0x100(%%rcx)					\n\t		vmovaps %%ymm3 ,0x120(%%rcx)				/* outC	*/	\n\t"\
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
		"vmulpd		(%%rdx)		,%%ymm0,%%ymm0		/* a[jt+p4]*c4 */		\n\t		vmulpd		    (%%rdi)	,%%ymm8 ,%%ymm8 		/* a[jt+p2]*c2 */	\n\t"\
		"vmulpd		(%%rdx)		,%%ymm1,%%ymm1		/* a[jp+p4]*c4 */		\n\t		vmulpd		    (%%rdi)	,%%ymm9 ,%%ymm9 		/* a[jp+p2]*c2 */	\n\t"\
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
		"vmulpd		(%%rdx)		,%%ymm4,%%ymm4		/* a[jt+p12]*c12 */		\n\t		vmulpd		    (%%rdi)	,%%ymm12,%%ymm12		/* a[jt+p10]*c10 */		\n\t"\
		"vmulpd		(%%rdx)		,%%ymm5,%%ymm5		/* a[jp+p12]*c12 */		\n\t		vmulpd		    (%%rdi)	,%%ymm13,%%ymm13		/* a[jp+p10]*c10 */		\n\t"\
		"vmulpd		0x020(%%rdx),%%ymm6,%%ymm6		/* a[jt+p12]*s12 */		\n\t		vmulpd		0x020(%%rdi),%%ymm14,%%ymm14		/* a[jt+p10]*s10 */		\n\t"\
		"vmulpd		0x020(%%rdx),%%ymm7,%%ymm7		/* a[jp+p12]*s12 */		\n\t		vmulpd		0x020(%%rdi),%%ymm15,%%ymm15		/* a[jp+p10]*s10 */		\n\t"\
		"vaddpd		%%ymm6		,%%ymm5,%%ymm5	/* ymm5 <- it */			\n\t		vaddpd		%%ymm14		,%%ymm13,%%ymm13	/* ymm13 <- it */						\n\t"\
		"vsubpd		%%ymm7		,%%ymm4,%%ymm4	/* ymm4 <- rt; ymm6,7 free*/\n\t		vsubpd		%%ymm15		,%%ymm12,%%ymm12	/* ymm12 <- rt		ymm14,7 free */		\n\t"\
		"vaddpd		%%ymm4		,%%ymm0,%%ymm0		/* ~t5 <- t5 +rt */		\n\t		vaddpd		%%ymm12		,%%ymm8 ,%%ymm8  	/* ~t13<- t13+rt */						\n\t"\
		"vaddpd		%%ymm5		,%%ymm1,%%ymm1		/* ~t6 <- t6 +it */		\n\t		vaddpd		%%ymm13		,%%ymm9 ,%%ymm9  	/* ~t14<- t14+it */						\n\t"\
		"vsubpd		%%ymm4		,%%ymm2,%%ymm2		/* ~t7 <- t5 -rt */		\n\t		vsubpd		%%ymm12		,%%ymm10,%%ymm10 	/* ~t15<- t13-rt */						\n\t"\
		"vsubpd		%%ymm5		,%%ymm3,%%ymm3		/* ~t8 <- t6 -it */		\n\t		vsubpd		%%ymm13		,%%ymm11,%%ymm11 	/* ~t16<- t14-it	ymm12,13 free */	\n\t"\
	"prefetcht1	0x80(%%r14)\n\t"\
		"\n\t"\
		"/* Now do the p0,8 combo: */\n\t									/* Do the p6,14 combo - do p14 first so registers come out in same order as for p2,10 */\n\t"\
		"leaq	0x4a0(%%rcx),%%rdx	/* c8 */								\n\t		leaq	0x620(%%rcx),%%rdi	/* c14 */\n\t"\
		"vmovaps	0x080(%%rcx)	,%%ymm4		/* a[jt+p8 ] */				\n\t		vmovaps		0x1c0(%%rcx),%%ymm12		/* a[jt+p14] */				\n\t"\
		"vmovaps	0x0a0(%%rcx)	,%%ymm5		/* a[jp+p8 ] */				\n\t		vmovaps		0x1e0(%%rcx),%%ymm13		/* a[jp+p14] */				\n\t"\
		"																				vmovaps			%%ymm12	,%%ymm14		/* ymm14 <- cpy a[jt+p14] */\n\t"\
		"vmovaps	%%ymm4		,%%ymm6	/* ymm6 <- cpy a[jt+p8] */			\n\t		vmovaps			%%ymm13	,%%ymm15		/* ymm15 <- cpy a[jp+p14] */\n\t"\
		"vmovaps	%%ymm5		,%%ymm7	/* ymm7 <- cpy a[jp+p8] */			\n\t		vmulpd		    (%%rdi)	,%%ymm12,%%ymm12		/* a[jt+p14]*c14 */			\n\t"\
		"																				vmulpd		    (%%rdi)	,%%ymm13,%%ymm13		/* a[jp+p14]*c14 */			\n\t"\
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
		"vsubpd		%%ymm4	,%%ymm6,%%ymm6	/* ~t3 <- t1 -rt */				\n\t		vmulpd		    (%%rdi)	,%%ymm12,%%ymm12		/* a[jt+p6]*c6 */			\n\t"\
		"vsubpd		%%ymm5	,%%ymm7,%%ymm7	/* ~t4 <- t2 -it */				\n\t		vmulpd		    (%%rdi)	,%%ymm13,%%ymm13		/* a[jp+p6]*c6 */			\n\t"\
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
	"prefetcht1	0x80(%%r15)\n\t"\
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
		"vmulpd		    (%%rbx)	,%%ymm4,%%ymm4		/* a[jt+p4]*c4 */		\n\t		vmulpd		0x100(%%rbx),%%ymm12,%%ymm12		/* a[jt+p4]*c4 */\n\t"\
		"vmulpd		    (%%rbx)	,%%ymm5,%%ymm5		/* a[jp+p4]*c4 */		\n\t		vmulpd		0x100(%%rbx),%%ymm13,%%ymm13		/* a[jp+p4]*c4 */\n\t"\
		"vmulpd		0x020(%%rbx),%%ymm6,%%ymm6		/* a[jt+p4]*s4 */		\n\t		vmulpd		0x120(%%rbx),%%ymm14,%%ymm14		/* a[jt+p4]*s4 */\n\t"\
		"vmulpd		0x020(%%rbx),%%ymm7,%%ymm7		/* a[jp+p4]*s4 */		\n\t		vmulpd		0x120(%%rbx),%%ymm15,%%ymm15		/* a[jp+p4]*s4 */\n\t"\
	"prefetcht1	0x100(%%r14)\n\t"\
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
		"movq	%[__r1] ,%%rax\n\t											vmovaps		0x080(%%rax),%%ymm8 		/* t5  */\n\t"\
		"movq	%[__isrt2],%%rsi\n\t										vmovaps		0x0a0(%%rax),%%ymm9 		/* t6  */\n\t"\
		"\n\t																vmovaps		0x1a0(%%rax),%%ymm11		/* t14 */\n\t"\
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
	"prefetcht1	0x100(%%r15)\n\t"\
		"vsubpd		%%ymm6,%%ymm4,%%ymm4		/* t25=t17-t25 */\n\t		vmulpd		(%%rsi),%%ymm14,%%ymm14	/*  rt = (t29+t30)*ISRT2 */\n\t"\
		"vsubpd		%%ymm7,%%ymm5,%%ymm5		/* t26=t18-t26 */\n\t		vmulpd		(%%rsi),%%ymm15,%%ymm15	/*  it = (t30-t29)*ISRT2 */\n\t"\
		"vaddpd		%%ymm6,%%ymm6,%%ymm6		/* 2.t25 */\n\t				vsubpd		%%ymm14,%%ymm12,%%ymm12		/* t21=t21-rt */\n\t"\
		"vaddpd		%%ymm7,%%ymm7,%%ymm7		/* 2.t26 */\n\t				vsubpd		%%ymm15,%%ymm13,%%ymm13		/* t22=t22-it */\n\t"\
		"vaddpd		%%ymm4,%%ymm6,%%ymm6		/* t17=t17+t25 */\n\t		vaddpd		%%ymm14,%%ymm14,%%ymm14		/*      2* rt */\n\t"\
		"vaddpd		%%ymm5,%%ymm7,%%ymm7		/* t18=t18+t26 */\n\t		vaddpd		%%ymm15,%%ymm15,%%ymm15		/*      2* it */\n\t"\
		"\n\t																vaddpd		%%ymm12,%%ymm14,%%ymm14		/* t29=t21+rt */\n\t"\
		"vsubpd		%%ymm6,%%ymm2,%%ymm2		/* t1  <- t1 -t17 */\n\t	vaddpd		%%ymm13,%%ymm15,%%ymm15		/* t30=t22+it */\n\t"\
		"vsubpd		%%ymm7,%%ymm3,%%ymm3		/* t2  <- t2 -t18 */\n\t	vsubpd		%%ymm12,%%ymm8 ,%%ymm8 		/* t5 -t21 */\n\t"\
		"vaddpd		%%ymm6,%%ymm6,%%ymm6		/*          2*t17 */\n\t	vsubpd		%%ymm13,%%ymm10,%%ymm10		/* t6 -t22 */\n\t"\
		"vaddpd		%%ymm7,%%ymm7,%%ymm7		/*          2*t18 */\n\t	vaddpd		%%ymm12,%%ymm12,%%ymm12		/*   2*t21 */\n\t"\
		"vaddpd		%%ymm2,%%ymm6,%%ymm6		/* t17 <- t1 +t17 */\n\t	vaddpd		%%ymm13,%%ymm13,%%ymm13		/*   2*t22 */\n\t"\
		"vaddpd		%%ymm3,%%ymm7,%%ymm7		/* t18 <- t2 +t18 */\n\t	vaddpd		%%ymm8 ,%%ymm12,%%ymm12		/* t5 +t21 */\n\t"\
		"/* x in 2, y in 3, spill 6 and use as tmp: */\n\t					vaddpd		%%ymm10,%%ymm13,%%ymm13		/* t6 +t22 */\n\t"\
		"vmovaps	%%ymm6,(%%rax)	/* tmp-store t17 in t0 */\n\t			/* x in 0, y in 2, spill 4 and use as tmp: */\n\t"\
		"vmovaps	%%ymm2,%%ymm6	/* cpy x */\n\t							vmovaps		%%ymm12,0x080(%%rax)	/* tmp-store t17 in t0 */\n\t"\
		"vsubpd		%%ymm3,%%ymm2,%%ymm2	/* x-y */\n\t					vmovaps		%%ymm8 ,%%ymm12	/* cpy x */\n\t"\
		"vaddpd		%%ymm3,%%ymm3,%%ymm3	/* 2*y */\n\t					vsubpd		%%ymm10,%%ymm8 ,%%ymm8 	/* x-y */\n\t"\
		"vmulpd		%%ymm3,%%ymm6,%%ymm6	/* 2xy */\n\t					vaddpd		%%ymm10,%%ymm10,%%ymm10	/* 2*y */\n\t"\
		"vaddpd		%%ymm2,%%ymm3,%%ymm3	/* x+y */\n\t					vmulpd		%%ymm10,%%ymm12,%%ymm12	/* 2xy */\n\t"\
		"vmulpd		%%ymm3,%%ymm2,%%ymm2	/* x^2-y^2 */\n\t				vaddpd		%%ymm8 ,%%ymm10,%%ymm10	/* x+y */\n\t"\
		"vmovaps	%%ymm6,0x020(%%rax)	/* a[jp+p1 ], store in t18 */\n\t	vmulpd		%%ymm10,%%ymm8 ,%%ymm8 	/* x^2-y^2 */\n\t"\
		"vmovaps	     (%%rax),%%ymm6	/* a[jt+p0 ], reload */\n\t			vmovaps		%%ymm12,0x0a0(%%rax)	/* a[jp+p1 ], store in t18 */\n\t"\
		"vmovaps	%%ymm2,     (%%rax)	/* a[jt+p1 ], store in t17 */\n\t	vmovaps	    0x080(%%rax),%%ymm12	/* a[jt+p0 ], reload */\n\t"\
		"/* Have 2 free regs for remaining 3 squarings: */\n\t				vmovaps		%%ymm8 ,0x080(%%rax)	/* a[jt+p1 ], store in t17 */\n\t"\
		"vmovaps	%%ymm6,%%ymm2	/* cpy x */\n\t							/* Have 2 free regs for remaining 3 squarings: */\n\t"\
		"vmovaps	%%ymm6,%%ymm3	/* cpy x */\n\t							vmovaps		%%ymm12,%%ymm8 	\n\t"\
		"vaddpd		%%ymm7,%%ymm6,%%ymm6	/* x+y */\n\t					vmovaps		%%ymm12,%%ymm10	\n\t"\
		"vsubpd		%%ymm7,%%ymm2,%%ymm2	/* x-y */\n\t					vaddpd		%%ymm13,%%ymm12,%%ymm12	\n\t"\
		"vaddpd		%%ymm7,%%ymm7,%%ymm7	/* 2*y */\n\t					vsubpd		%%ymm13,%%ymm8 ,%%ymm8 	\n\t"\
		"vmulpd		%%ymm2,%%ymm6,%%ymm6	/* x^2-y^2 */\n\t				vaddpd		%%ymm13,%%ymm13,%%ymm13	\n\t"\
		"vmulpd		%%ymm3,%%ymm7,%%ymm7	/* 2xy */\n\t					vmulpd		%%ymm8 ,%%ymm12,%%ymm12	\n\t"\
		"\n\t																vmulpd		%%ymm10,%%ymm13,%%ymm13	\n\t"\
		"vsubpd		%%ymm5,%%ymm0,%%ymm0		/* t9  <- t9 -t26 */\n\t	vsubpd		%%ymm15,%%ymm11,%%ymm11	\n\t"\
		"vsubpd		%%ymm4,%%ymm1,%%ymm1		/* t10 <- t10-t25 */\n\t	vsubpd		%%ymm14,%%ymm9 ,%%ymm9 	\n\t"\
		"vaddpd		%%ymm5,%%ymm5,%%ymm5		/*          2*t26 */\n\t	vaddpd		%%ymm15,%%ymm15,%%ymm15	\n\t"\
		"vaddpd		%%ymm4,%%ymm4,%%ymm4		/*          2*t25 */\n\t	vaddpd		%%ymm14,%%ymm14,%%ymm14	\n\t"\
		"vaddpd		%%ymm0,%%ymm5,%%ymm5		/* t26 <- t9 +t26 */\n\t	vaddpd		%%ymm11,%%ymm15,%%ymm15	\n\t"\
		"vaddpd		%%ymm1,%%ymm4,%%ymm4		/* t25 <- t10+t25 */\n\t	vaddpd		%%ymm9 ,%%ymm14,%%ymm14	\n\t"\
		"vmovaps	%%ymm5,%%ymm2	/* cpy x */\n\t							vmovaps		%%ymm15,%%ymm8 	\n\t"\
		"vmovaps	%%ymm5,%%ymm3	/* cpy x */\n\t							vmovaps		%%ymm15,%%ymm10	\n\t"\
		"vaddpd		%%ymm1,%%ymm5,%%ymm5	/* x+y */\n\t					vaddpd		%%ymm9 ,%%ymm15,%%ymm15	\n\t"\
		"vsubpd		%%ymm1,%%ymm2,%%ymm2	/* x-y */\n\t					vsubpd		%%ymm9 ,%%ymm8 ,%%ymm8 	\n\t"\
		"vaddpd		%%ymm1,%%ymm1,%%ymm1	/* 2*y */\n\t					vaddpd		%%ymm9 ,%%ymm9 ,%%ymm9 	\n\t"\
		"vmulpd		%%ymm2,%%ymm5,%%ymm5	/* x^2-y^2 */\n\t				vmulpd		%%ymm8 ,%%ymm15,%%ymm15	\n\t"\
		"vmulpd		%%ymm3,%%ymm1,%%ymm1	/* 2xy */\n\t					vmulpd		%%ymm10,%%ymm9 ,%%ymm9 	\n\t"\
		"\n\t"\
		"vmovaps	%%ymm0,%%ymm2	/* cpy x */\n\t							vmovaps		%%ymm11,%%ymm8 	\n\t"\
		"vmovaps	%%ymm0,%%ymm3	/* cpy x */\n\t							vmovaps		%%ymm11,%%ymm10	\n\t"\
		"vaddpd		%%ymm4,%%ymm0,%%ymm0	/* x+y */\n\t					vaddpd		%%ymm14,%%ymm11,%%ymm11	\n\t"\
		"vsubpd		%%ymm4,%%ymm2,%%ymm2	/* x-y */\n\t					vsubpd		%%ymm14,%%ymm8 ,%%ymm8 	\n\t"\
		"vaddpd		%%ymm4,%%ymm4,%%ymm4	/* 2*y */\n\t					vaddpd		%%ymm14,%%ymm14,%%ymm14	\n\t"\
		"vmulpd		%%ymm2,%%ymm0,%%ymm0	/* x^2-y^2 */\n\t				vmulpd		%%ymm8 ,%%ymm11,%%ymm11	\n\t"\
		"vmulpd		%%ymm3,%%ymm4,%%ymm4	/* 2xy */\n\t					vmulpd		%%ymm10,%%ymm14,%%ymm14	\n\t"\
		"\n\t"\
		"vmovaps	     (%%rax),%%ymm2	/* a[jt+p1 ], reload */\n\t			vmovaps	0x080(%%rax),%%ymm8 	/* a[jt+p1 ], reload */\n\t"\
		"vmovaps	0x020(%%rax),%%ymm3	/* a[jp+p1 ], reload */\n\t			vmovaps	0x0a0(%%rax),%%ymm10	/* a[jp+p1 ], reload */\n\t"\
	"prefetcht1	0x180(%%r14)\n\t"\
	"/* SSE2_RADIX4_DIT_IN_PLACE_C(ymm6,7,2,3,0,4,5,1) - stores all 8 memlocs*/\n\t	/* SSE2_RADIX4_DIT_IN_PLACE_C(ymm12,5,0,2,3,6,7,1): */\n\t"\
		"vsubpd		%%ymm2,%%ymm6,%%ymm6	/* t3 */		\n\t			vsubpd	%%ymm8 ,%%ymm12,%%ymm12		\n\t"\
		"vsubpd		%%ymm3,%%ymm7,%%ymm7	/* t4 */		\n\t			vsubpd	%%ymm10,%%ymm13,%%ymm13		\n\t"\
		"vaddpd		%%ymm2,%%ymm2,%%ymm2	/* 2*y */		\n\t			vaddpd	%%ymm8 ,%%ymm8 ,%%ymm8 		\n\t"\
		"vaddpd		%%ymm3,%%ymm3,%%ymm3	/* 2*y */		\n\t			vaddpd	%%ymm10,%%ymm10,%%ymm10		\n\t"\
		"vaddpd		%%ymm6,%%ymm2,%%ymm2	/* t1 */		\n\t			vaddpd	%%ymm12,%%ymm8 ,%%ymm8 		\n\t"\
		"vaddpd		%%ymm7,%%ymm3,%%ymm3	/* t2 */		\n\t			vaddpd	%%ymm13,%%ymm10,%%ymm10		\n\t"\
		"vsubpd		%%ymm5,%%ymm0,%%ymm0	/* t7 */		\n\t			vsubpd	%%ymm15,%%ymm11,%%ymm11		\n\t"\
		"vsubpd		%%ymm1,%%ymm4,%%ymm4	/* t8 */		\n\t			vsubpd	%%ymm9 ,%%ymm14,%%ymm14		\n\t"\
		"vsubpd		%%ymm0,%%ymm7,%%ymm7	/* ~t4 <- t4 -t7 */	\n\t		vsubpd	%%ymm11,%%ymm13,%%ymm13		\n\t"\
		"vaddpd		%%ymm5,%%ymm5,%%ymm5	/* 2*y */		\n\t			vaddpd	%%ymm15,%%ymm15,%%ymm15		\n\t"\
		"vsubpd		%%ymm4,%%ymm6,%%ymm6	/* ~t7 <- t3 -t8 */	\n\t		vsubpd	%%ymm14,%%ymm12,%%ymm12		\n\t"\
		"vaddpd		%%ymm1,%%ymm1,%%ymm1	/* 2*y */		\n\t			vaddpd	%%ymm9 ,%%ymm9 ,%%ymm9 		\n\t"\
		"vmovaps	%%ymm7,0x120(%%rax)	/* <- ~t4 */\n\t					vmovaps	%%ymm13,0x1a0(%%rax)		\n\t"\
		"vaddpd		%%ymm0,%%ymm5,%%ymm5	/* t4 */		\n\t			vaddpd	%%ymm11,%%ymm15,%%ymm15		\n\t"\
		"vmovaps	%%ymm6,0x300(%%rax)	/* <- ~t7 */\n\t					vmovaps	%%ymm12,0x380(%%rax)		\n\t"\
		"vaddpd		%%ymm4,%%ymm1,%%ymm1	/* t5 */		\n\t			vaddpd	%%ymm14,%%ymm9 ,%%ymm9 		\n\t"\
		"vaddpd		%%ymm4,%%ymm4,%%ymm4	/*          2*t8 */		\n\t	vaddpd	%%ymm14,%%ymm14,%%ymm14		\n\t"\
		"vsubpd		%%ymm5,%%ymm2,%%ymm2	/* ~t5 <- t1 -t5 */		\n\t	vsubpd	%%ymm15,%%ymm8 ,%%ymm8 		\n\t"\
		"vaddpd		%%ymm0,%%ymm0,%%ymm0	/*          2*t7 */		\n\t	vaddpd	%%ymm11,%%ymm11,%%ymm11		\n\t"\
		"vsubpd		%%ymm1,%%ymm3,%%ymm3	/* ~t6 <- t2 -t6 */		\n\t	vsubpd	%%ymm9 ,%%ymm10,%%ymm10		\n\t"\
		"vaddpd		%%ymm6,%%ymm4,%%ymm4	/* ~t3 <- t3 +t8 */		\n\t	vaddpd	%%ymm12,%%ymm14,%%ymm14		\n\t"\
		"vmovaps	%%ymm2,0x200(%%rax)	/* <- ~t5 */	\n\t				vmovaps	%%ymm8 ,0x280(%%rax)		\n\t"\
		"vaddpd		%%ymm7,%%ymm0,%%ymm0	/* ~t8 <- t4 +t7 */		\n\t	vaddpd	%%ymm13,%%ymm11,%%ymm11		\n\t"\
		"vmovaps	%%ymm3,0x220(%%rax)	/* <- ~t6 */	\n\t				vmovaps	%%ymm10,0x2a0(%%rax)		\n\t"\
		"vmovaps	%%ymm4,0x100(%%rax)	/* <- ~t3 */	\n\t				vmovaps	%%ymm14,0x180(%%rax)		\n\t"\
		"vaddpd		%%ymm5,%%ymm5,%%ymm5	/*          2*t5 */		\n\t	vaddpd	%%ymm15,%%ymm15,%%ymm15		\n\t"\
		"vmovaps	%%ymm0,0x320(%%rax)	/* <- ~t8 */	\n\t				vmovaps	%%ymm11,0x3a0(%%rax)		\n\t"\
		"vaddpd		%%ymm1,%%ymm1,%%ymm1	/*          2*t6 */		\n\t	vaddpd	%%ymm9 ,%%ymm9 ,%%ymm9 		\n\t"\
		"vaddpd		%%ymm2,%%ymm5,%%ymm5	/* ~t1 <- t1 +t5 */		\n\t	vaddpd	%%ymm8 ,%%ymm15,%%ymm15		\n\t"\
		"vaddpd		%%ymm3,%%ymm1,%%ymm1	/* ~t2 <- t2 +t6 */		\n\t	vaddpd	%%ymm10,%%ymm9 ,%%ymm9 		\n\t"\
		"vmovaps	%%ymm5,     (%%rax)	/* <- ~t1 */	\n\t				vmovaps	%%ymm15,0x080(%%rax)		\n\t"\
		"vmovaps	%%ymm1,0x020(%%rax)	/* <- ~t2 */	\n\t				vmovaps	%%ymm9 ,0x0a0(%%rax)		\n\t"\
		"\n\t"\
"/*...Block 2: t3,11,19,27 */\n\t"\
		"addq	$0x040,%%rax	/* r3  */\n\t"\
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
	"prefetcht1	0x180(%%r15)\n\t"\
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
		"vaddpd		%%ymm2,%%ymm6,%%ymm6		/* t3 +t19 */\n\t			vaddpd		%%ymm8 		,%%ymm12,%%ymm12	/* t7 +t23 */\n\t"\
		"vaddpd		%%ymm3,%%ymm7,%%ymm7		/* t4 +t20 */\n\t			vaddpd		%%ymm9 		,%%ymm13,%%ymm13	/* t8 +t24 */\n\t"\
		"/* x in 2, y in 3, spill 6 and use as tmp: */\n\t					/* x in 0, y in 1, spill 4 and use as tmp: */\n\t"\
		"vmovaps	%%ymm6,(%%rax)	/* tmp-store t17 in t0 */\n\t			vmovaps	%%ymm12,0x080(%%rax)	/* tmp-store t17 in t0 */\n\t"\
		"vmovaps	%%ymm2,%%ymm6	/* cpy x */\n\t							vmovaps	%%ymm8 ,%%ymm12	/* cpy x */\n\t"\
		"vsubpd		%%ymm3,%%ymm2,%%ymm2	/* x-y */\n\t					vsubpd	%%ymm9 ,%%ymm8 ,%%ymm8 	/* x-y */\n\t"\
		"vaddpd		%%ymm3,%%ymm3,%%ymm3	/* 2*y */\n\t					vaddpd	%%ymm9 ,%%ymm9 ,%%ymm9 	/* 2*y */\n\t"\
		"vmulpd		%%ymm3,%%ymm6,%%ymm6	/* 2xy */\n\t					vmulpd	%%ymm9 ,%%ymm12,%%ymm12	/* 2xy */\n\t"\
		"vaddpd		%%ymm2,%%ymm3,%%ymm3	/* x+y */\n\t					vaddpd	%%ymm8 ,%%ymm9 ,%%ymm9 	/* x+y */\n\t"\
		"vmulpd		%%ymm3,%%ymm2,%%ymm2	/* x^2-y^2 */\n\t				vmulpd	%%ymm9 ,%%ymm8 ,%%ymm8 	/* x^2-y^2 */\n\t"\
		"vmovaps	%%ymm6,0x020(%%rax)	/* a[jp+p1 ], store in t18 */\n\t	vmovaps		%%ymm12,0x0a0(%%rax)	/* a[jp+p1 ], store in t18 */\n\t"\
		"vmovaps	     (%%rax),%%ymm6	/* a[jt+p0 ], reload */\n\t			vmovaps	    0x080(%%rax),%%ymm12	/* a[jt+p0 ], reload */\n\t"\
		"vmovaps	%%ymm2,     (%%rax)	/* a[jt+p1 ], store in t17 */\n\t	vmovaps		%%ymm8 ,0x080(%%rax)	/* a[jt+p1 ], store in t17 */\n\t"\
		"/* Have 2 free regs for remaining 3 squarings: */\n\t				/* Have 2 free regs for remaining 3 squarings: */\n\t"\
		"vmovaps	%%ymm6,%%ymm2	\n\t									vmovaps	%%ymm12,%%ymm8 		\n\t"\
		"vmovaps	%%ymm6,%%ymm3	\n\t									vmovaps	%%ymm12,%%ymm9 		\n\t"\
		"vaddpd		%%ymm7,%%ymm6,%%ymm6	\n\t							vaddpd	%%ymm13,%%ymm12,%%ymm12		\n\t"\
		"vsubpd		%%ymm7,%%ymm2,%%ymm2	\n\t							vsubpd	%%ymm13,%%ymm8 ,%%ymm8 		\n\t"\
		"vaddpd		%%ymm7,%%ymm7,%%ymm7	\n\t							vaddpd	%%ymm13,%%ymm13,%%ymm13		\n\t"\
		"vmulpd		%%ymm2,%%ymm6,%%ymm6	\n\t							vmulpd	%%ymm8 ,%%ymm12,%%ymm12		\n\t"\
		"vmulpd		%%ymm3,%%ymm7,%%ymm7	\n\t							vmulpd	%%ymm9 ,%%ymm13,%%ymm13		\n\t"\
		"\n\t"\
		"vsubpd		%%ymm5,%%ymm0,%%ymm0	\n\t							vsubpd	%%ymm15,%%ymm10,%%ymm10		\n\t"\
		"vsubpd		%%ymm4,%%ymm1,%%ymm1	\n\t							vsubpd	%%ymm14,%%ymm11,%%ymm11		\n\t"\
		"vaddpd		%%ymm5,%%ymm5,%%ymm5	\n\t							vaddpd	%%ymm15,%%ymm15,%%ymm15		\n\t"\
		"vaddpd		%%ymm4,%%ymm4,%%ymm4	\n\t							vaddpd	%%ymm14,%%ymm14,%%ymm14		\n\t"\
		"vaddpd		%%ymm0,%%ymm5,%%ymm5	\n\t							vaddpd	%%ymm10,%%ymm15,%%ymm15		\n\t"\
		"vaddpd		%%ymm1,%%ymm4,%%ymm4	\n\t							vaddpd	%%ymm11,%%ymm14,%%ymm14		\n\t"\
		"vmovaps	%%ymm5,%%ymm2	\n\t									vmovaps	%%ymm15,%%ymm8 		\n\t"\
		"vmovaps	%%ymm5,%%ymm3	\n\t									vmovaps	%%ymm15,%%ymm9 		\n\t"\
		"vaddpd		%%ymm1,%%ymm5,%%ymm5	\n\t							vaddpd	%%ymm11,%%ymm15,%%ymm15		\n\t"\
		"vsubpd		%%ymm1,%%ymm2,%%ymm2	\n\t							vsubpd	%%ymm11,%%ymm8 ,%%ymm8 		\n\t"\
		"vaddpd		%%ymm1,%%ymm1,%%ymm1	\n\t							vaddpd	%%ymm11,%%ymm11,%%ymm11		\n\t"\
		"vmulpd		%%ymm2,%%ymm5,%%ymm5	\n\t							vmulpd	%%ymm8 ,%%ymm15,%%ymm15		\n\t"\
		"vmulpd		%%ymm3,%%ymm1,%%ymm1	\n\t							vmulpd	%%ymm9 ,%%ymm11,%%ymm11		\n\t"\
	"prefetcht1	0x200(%%r14)\n\t"\
		"\n\t"\
		"vmovaps	%%ymm0,%%ymm2	\n\t									vmovaps	%%ymm10,%%ymm8 		\n\t"\
		"vmovaps	%%ymm0,%%ymm3	\n\t									vmovaps	%%ymm10,%%ymm9 		\n\t"\
		"vaddpd		%%ymm4,%%ymm0,%%ymm0	\n\t							vaddpd	%%ymm14,%%ymm10,%%ymm10		\n\t"\
		"vsubpd		%%ymm4,%%ymm2,%%ymm2	\n\t							vsubpd	%%ymm14,%%ymm8 ,%%ymm8 		\n\t"\
		"vaddpd		%%ymm4,%%ymm4,%%ymm4	\n\t							vaddpd	%%ymm14,%%ymm14,%%ymm14		\n\t"\
		"vmulpd		%%ymm2,%%ymm0,%%ymm0	\n\t							vmulpd	%%ymm8 ,%%ymm10,%%ymm10		\n\t"\
		"vmulpd		%%ymm3,%%ymm4,%%ymm4	\n\t							vmulpd	%%ymm9 ,%%ymm14,%%ymm14		\n\t"\
		"\n\t"\
		"vmovaps	      (%%rax),%%ymm2	/* a[jt+p1 ], reload */\n\t		vmovaps	 0x080(%%rax),%%ymm8 	/* a[jt+p1 ], reload */		\n\t"\
		"vmovaps	 0x020(%%rax),%%ymm3	/* a[jp+p1 ], reload */\n\t		vmovaps	 0x0a0(%%rax),%%ymm9 	/* a[jp+p1 ], reload */		\n\t"\
	"/* SSE2_RADIX4_DIT_IN_PLACE_C(ymm6,ymm7,ymm2,ymm3,ymm0,ymm4,ymm5,ymm1) SSE2_RADIX4_DIT_IN_PLACE_C(ymm12,ymm13,ymm8 ,ymm9 ,ymm10,ymm14,ymm15,ymm11) - This stores all 8 memlocs */		\n\t"\
		"vsubpd		%%ymm2,%%ymm6,%%ymm6	\n\t							vsubpd	%%ymm8 ,%%ymm12,%%ymm12		\n\t"\
		"vsubpd		%%ymm3,%%ymm7,%%ymm7	\n\t							vsubpd	%%ymm9 ,%%ymm13,%%ymm13		\n\t"\
		"vaddpd		%%ymm2,%%ymm2,%%ymm2	\n\t							vaddpd	%%ymm8 ,%%ymm8 ,%%ymm8 		\n\t"\
		"vaddpd		%%ymm3,%%ymm3,%%ymm3	\n\t							vaddpd	%%ymm9 ,%%ymm9 ,%%ymm9 		\n\t"\
		"vaddpd		%%ymm6,%%ymm2,%%ymm2	\n\t							vaddpd	%%ymm12,%%ymm8 ,%%ymm8 		\n\t"\
		"vaddpd		%%ymm7,%%ymm3,%%ymm3	\n\t							vaddpd	%%ymm13,%%ymm9 ,%%ymm9 		\n\t"\
		"vsubpd		%%ymm5,%%ymm0,%%ymm0	\n\t							vsubpd	%%ymm15,%%ymm10,%%ymm10		\n\t"\
		"vsubpd		%%ymm1,%%ymm4,%%ymm4	\n\t							vsubpd	%%ymm11,%%ymm14,%%ymm14		\n\t"\
		"vsubpd		%%ymm0,%%ymm7,%%ymm7	\n\t							vsubpd	%%ymm10,%%ymm13,%%ymm13		\n\t"\
		"vaddpd		%%ymm5,%%ymm5,%%ymm5	\n\t							vaddpd	%%ymm15,%%ymm15,%%ymm15		\n\t"\
		"vsubpd		%%ymm4,%%ymm6,%%ymm6	\n\t							vsubpd	%%ymm14,%%ymm12,%%ymm12		\n\t"\
		"vaddpd		%%ymm1,%%ymm1,%%ymm1	\n\t							vaddpd	%%ymm11,%%ymm11,%%ymm11		\n\t"\
		"vmovaps	%%ymm7,0x120(%%rax)	\n\t								vmovaps	%%ymm13,0x1a0(%%rax)		\n\t"\
		"vaddpd		%%ymm0,%%ymm5,%%ymm5	\n\t							vaddpd	%%ymm10,%%ymm15,%%ymm15		\n\t"\
		"vmovaps	%%ymm6,0x300(%%rax)	\n\t								vmovaps	%%ymm12,0x380(%%rax)		\n\t"\
		"vaddpd		%%ymm4,%%ymm1,%%ymm1	\n\t							vaddpd	%%ymm14,%%ymm11,%%ymm11		\n\t"\
		"vaddpd		%%ymm4,%%ymm4,%%ymm4	\n\t							vaddpd	%%ymm14,%%ymm14,%%ymm14		\n\t"\
		"vsubpd		%%ymm5,%%ymm2,%%ymm2	\n\t							vsubpd	%%ymm15,%%ymm8 ,%%ymm8 		\n\t"\
		"vaddpd		%%ymm0,%%ymm0,%%ymm0	\n\t							vaddpd	%%ymm10,%%ymm10,%%ymm10		\n\t"\
		"vsubpd		%%ymm1,%%ymm3,%%ymm3	\n\t							vsubpd	%%ymm11,%%ymm9 ,%%ymm9 		\n\t"\
		"vaddpd		%%ymm6,%%ymm4,%%ymm4	\n\t							vaddpd	%%ymm12,%%ymm14,%%ymm14		\n\t"\
		"vmovaps	%%ymm2,0x200(%%rax)	\n\t								vmovaps	%%ymm8 ,0x280(%%rax)		\n\t"\
		"vaddpd		%%ymm7,%%ymm0,%%ymm0	\n\t							vaddpd	%%ymm13,%%ymm10,%%ymm10		\n\t"\
		"vmovaps	%%ymm3,0x220(%%rax)	\n\t								vmovaps	%%ymm9 ,0x2a0(%%rax)		\n\t"\
		"vmovaps	%%ymm4,0x100(%%rax)	\n\t								vmovaps	%%ymm14,0x180(%%rax)		\n\t"\
		"vaddpd		%%ymm5,%%ymm5,%%ymm5	\n\t							vaddpd	%%ymm15,%%ymm15,%%ymm15		\n\t"\
		"vmovaps	%%ymm0,0x320(%%rax)	\n\t								vmovaps	%%ymm10,0x3a0(%%rax)		\n\t"\
		"vaddpd		%%ymm1,%%ymm1,%%ymm1	\n\t							vaddpd	%%ymm11,%%ymm11,%%ymm11		\n\t"\
		"vaddpd		%%ymm2,%%ymm5,%%ymm5	\n\t							vaddpd	%%ymm8 ,%%ymm15,%%ymm15		\n\t"\
		"vaddpd		%%ymm3,%%ymm1,%%ymm1	\n\t							vaddpd	%%ymm9 ,%%ymm11,%%ymm11		\n\t"\
		"vmovaps	%%ymm5,     (%%rax)	\n\t								vmovaps	%%ymm15,0x080(%%rax)		\n\t"\
		"vmovaps	%%ymm1,0x020(%%rax)	\n\t								vmovaps	%%ymm11,0x0a0(%%rax)		\n\t"\
		"\n\t"\
	"/***************************************************************************************************/\n\t"\
	"/*..and now do four more radix-4 transforms, including the internal and external twiddle factors.  */\n\t"\
	"/*  Write even-index 16-byte output pairs to a(j1), odd-index to a(j2), unpack same as on inputs.  */\n\t"\
	"/*  We do the last 2 radix-4 blocks first, to make the unpack-interleaving of outputs smoother.    */\n\t"\
	"/***************************************************************************************************/\n\t"\
	"/* Main-array addresses still in add0,1, no need to re-init: */\n\t"\
		"/*...Block 3: t3,11,19,27 -> r9,13,11,15: */	\n\t				/*...Block 4: t7,15,23,31 -> r25,29,27,31: */\n\t"\
		"addq	$0x0c0,%%rax		/* r9 */	\n\t						/* r25 */\n\t"\
		"movq	%[__isrt2],%%rbx				\n\t						/* All rax-offsets incr +0x200 in rcol w.r.to lcol: */\n\t"\
		"movq	%%rdi,%%rcx		/* cc0 */		\n\t"\
	"prefetcht1	0x200(%%r15)\n\t"\
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
		"movq	%[__add1],%%rbx	/* rbx shared between rcol/lcol; rcx/rdx-offsets incr +0x80 in rcol for rest of block: */\n\t"\
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
	"prefetcht1	0x280(%%r14)\n\t"\
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
	"prefetcht1	0x280(%%r15)\n\t"\
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
		"movq	%[__add1],%%rcx				\n\t							vmulpd		0x0c0(%%rdx),%%ymm12,%%ymm12	/* c2 */		\n\t"\
		"																	vmulpd		0x0c0(%%rdx),%%ymm13,%%ymm13		\n\t"\
		"																	vmulpd		0x0e0(%%rdx),%%ymm11,%%ymm11		\n\t"\
		"																	vmulpd		0x0e0(%%rdx),%%ymm9 ,%%ymm9 		\n\t"\
		"																	vsubpd		%%ymm11     ,%%ymm13,%%ymm13		\n\t"\
		"																	vaddpd		%%ymm9      ,%%ymm12,%%ymm12		\n\t"\
"vmovapd	0x120(%%rcx),%%ymm7		\n\t"\
"vmovapd	0x100(%%rcx),%%ymm6		\n\t"\
"vmovapd	%%ymm3,0x120(%%rax)	/* r9 */\n\t"\
"vmovapd	%%ymm2,0x100(%%rax)	/* r8 */\n\t"\
"vmovapd	%%ymm7,0x320(%%rax)	/* r25 */\n\t"\
"vmovapd	%%ymm6,0x300(%%rax)	/* r24 */\n\t"\
																	"vmovapd	0x060(%%rcx),%%ymm11	\n\t"\
																	"vmovapd	0x040(%%rcx),%%ymm9		\n\t"\
																	"vmovapd	%%ymm13,0x060(%%rax)	/* r3 */\n\t"\
																	"vmovapd	%%ymm12,0x040(%%rax)	/* r2 */\n\t"\
																	"vmovapd	%%ymm11,0x260(%%rax)	/* r19 */\n\t"\
																	"vmovapd	%%ymm9 ,0x240(%%rax)	/* r18 */\n\t"\
		"																	addq	$0x040,%%rdx		/* c4 in lcol; c10 in rcol*/		\n\t"\
		"																	vmovapd		     (%%rax),%%ymm12		\n\t"\
		"																	vmovapd		0x020(%%rax),%%ymm13		\n\t"\
		"																	vmovapd		%%ymm12,%%ymm11		\n\t"\
		"																	vmovapd		%%ymm13,%%ymm9 		\n\t"\
		"																	vmulpd		0x0c0(%%rdx),%%ymm12,%%ymm12		\n\t"\
		"																	vmulpd		0x0c0(%%rdx),%%ymm13,%%ymm13		\n\t"\
		"																	vmulpd		0x0e0(%%rdx),%%ymm11,%%ymm11		\n\t"\
		"																	vmulpd		0x0e0(%%rdx),%%ymm9 ,%%ymm9 		\n\t"\
"vmovapd	0x020(%%rbx),%%ymm3		\n\t"\
"vmovapd	0x000(%%rbx),%%ymm2		\n\t"\
"vmovapd	0x020(%%rcx),%%ymm7		\n\t"\
"vmovapd	0x000(%%rcx),%%ymm6		\n\t"\
"vmovapd	%%ymm3,0x020(%%rax)	/* r1 */\n\t"\
"vmovapd	%%ymm2,0x000(%%rax)	/* r0 */\n\t"\
"vmovapd	%%ymm7,0x220(%%rax)	/* r17 */\n\t"\
"vmovapd	%%ymm6,0x200(%%rax)	/* r16 */\n\t"\
		"vaddpd		%%ymm5,%%ymm0,%%ymm0		\n\t						vsubpd		%%ymm11     ,%%ymm13,%%ymm13		\n\t"\
		"vsubpd		%%ymm4,%%ymm1,%%ymm1		\n\t						vaddpd		%%ymm9      ,%%ymm12,%%ymm12		\n\t"\
		"vmovapd	%%ymm0,%%ymm2				\n\t"\
		"vmovapd	%%ymm1,%%ymm3				\n\t"\
		"vaddpd		%%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vaddpd		%%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vmovapd	%%ymm0,%%ymm6				\n\t"\
		"vmovapd	%%ymm1,%%ymm7				\n\t"\
		"vmulpd		     (%%rdx),%%ymm2,%%ymm2	\n\t"\
		"vmulpd		     (%%rdx),%%ymm3,%%ymm3	\n\t"\
		"vmulpd		0x020(%%rdx),%%ymm6,%%ymm6	\n\t"\
		"vmulpd		0x020(%%rdx),%%ymm7,%%ymm7	\n\t"\
	"prefetcht1	0x300(%%r14)\n\t"\
																	"vmovapd	0x160(%%rcx),%%ymm11	\n\t"\
																	"vmovapd	0x140(%%rcx),%%ymm9		\n\t"\
																	"vmovapd	%%ymm13,0x160(%%rax)	/* r11 */\n\t"\
																	"vmovapd	%%ymm12,0x140(%%rax)	/* r10 */\n\t"\
																	"vmovapd	%%ymm11,0x360(%%rax)	/* r27 */\n\t"\
																	"vmovapd	%%ymm9 ,0x340(%%rax)	/* r26 */\n\t"\
		"addq	$0x040,%%rdx		/* c12 in lcol; c6 in rcol*/\n\t		vsubpd		%%ymm15,%%ymm8 ,%%ymm8 		\n\t"\
		"vsubpd		%%ymm6,%%ymm3,%%ymm3		\n\t						vsubpd		%%ymm14,%%ymm10,%%ymm10		\n\t"\
		"vaddpd		%%ymm7,%%ymm2,%%ymm2		\n\t						vaddpd		%%ymm15,%%ymm15,%%ymm15		\n\t"\
		"																	vaddpd		%%ymm14,%%ymm14,%%ymm14		\n\t"\
		"																	vaddpd		%%ymm8 ,%%ymm15,%%ymm15		\n\t"\
		"																	vaddpd		%%ymm10,%%ymm14,%%ymm14		\n\t"\
		"																	vmovapd		%%ymm15,%%ymm12		\n\t"\
		"																	vmovapd		%%ymm10,%%ymm13		\n\t"\
		"																	vmulpd		0x0c0(%%rdx),%%ymm15,%%ymm15		\n\t"\
		"																	vmulpd		0x0c0(%%rdx),%%ymm10,%%ymm10		\n\t"\
		"																	vmulpd		0x0e0(%%rdx),%%ymm12,%%ymm12		\n\t"\
		"																	vmulpd		0x0e0(%%rdx),%%ymm13,%%ymm13		\n\t"\
		"																	vsubpd		%%ymm12     ,%%ymm10,%%ymm10		\n\t"\
"vmovapd	0x0a0(%%rcx),%%ymm7		\n\t"\
"vmovapd	0x080(%%rcx),%%ymm6		\n\t"\
"vmovapd	%%ymm3,0x0a0(%%rax)	/* r5 */\n\t"\
"vmovapd	%%ymm2,0x080(%%rax)	/* r4 */\n\t"\
"vmovapd	%%ymm7,0x2a0(%%rax)	/* r21 */\n\t"\
"vmovapd	%%ymm6,0x280(%%rax)	/* r20 */\n\t"\
		"vsubpd		%%ymm5,%%ymm0,%%ymm0		\n\t						vaddpd		%%ymm13     ,%%ymm15,%%ymm15		\n\t"\
		"vaddpd		%%ymm4,%%ymm1,%%ymm1		\n\t"\
		"vmovapd	%%ymm0,%%ymm6				\n\t"\
		"vmovapd	%%ymm1,%%ymm7				\n\t"\
		"vmulpd		     (%%rdx),%%ymm0,%%ymm0	\n\t"\
		"vmulpd		     (%%rdx),%%ymm1,%%ymm1	\n\t"\
		"vmulpd		0x020(%%rdx),%%ymm6,%%ymm6	\n\t"\
		"vmulpd		0x020(%%rdx),%%ymm7,%%ymm7	\n\t"\
		"vsubpd		%%ymm6,%%ymm1,%%ymm1		\n\t"\
		"vaddpd		%%ymm7,%%ymm0,%%ymm0		\n\t"\
																	"vmovapd	0x0e0(%%rcx),%%ymm13	\n\t"\
																	"vmovapd	0x0c0(%%rcx),%%ymm12	\n\t"\
																	"vmovapd	%%ymm10,0x0e0(%%rax)	/* r7 */\n\t"\
																	"vmovapd	%%ymm15,0x0c0(%%rax)	/* r6 */\n\t"\
																	"vmovapd	%%ymm13,0x2e0(%%rax)	/* r23 */\n\t"\
																	"vmovapd	%%ymm12,0x2c0(%%rax)	/* r22 */\n\t"\
		"																	vmovapd		%%ymm8 ,%%ymm12		\n\t"\
		"																	vmovapd		%%ymm14,%%ymm13		\n\t"\
		"																	vmulpd		0x100(%%rdx),%%ymm8 ,%%ymm8 	/* c14 */		\n\t"\
		"																	vmulpd		0x100(%%rdx),%%ymm14,%%ymm14		\n\t"\
		"																	vmulpd		0x120(%%rdx),%%ymm12,%%ymm12		\n\t"\
		"																	vmulpd		0x120(%%rdx),%%ymm13,%%ymm13		\n\t"\
		"																	vsubpd		%%ymm12     ,%%ymm14,%%ymm14		\n\t"\
		"																	vaddpd		%%ymm13     ,%%ymm8 ,%%ymm8 		\n\t"\
"vmovapd	0x1a0(%%rcx),%%ymm7				\n\t"\
"vmovapd	0x180(%%rcx),%%ymm6				\n\t"\
"vmovapd	%%ymm1,0x1a0(%%rax)	/* r13 */	\n\t"\
"vmovapd	%%ymm0,0x180(%%rax)	/* r12 */	\n\t"\
"vmovapd	%%ymm7,0x3a0(%%rax)	/* r29 */	\n\t"\
"vmovapd	%%ymm6,0x380(%%rax)	/* r28 */	\n\t"\
	"prefetcht1	0x300(%%r15)\n\t"\
																	"vmovapd	0x1e0(%%rcx),%%ymm13	\n\t"\
																	"vmovapd	0x1c0(%%rcx),%%ymm12	\n\t"\
																	"vmovapd	%%ymm14,0x1e0(%%rax)	/* r15 */\n\t"\
																	"vmovapd	%%ymm8 ,0x1c0(%%rax)	/* r14 */\n\t"\
																	"vmovapd	%%ymm13,0x3e0(%%rax)	/* r31 */\n\t"\
																	"vmovapd	%%ymm12,0x3c0(%%rax)	/* r30 */\n\t"\
/*==========================*/"\n\t"\
	"/**** Finish with 4-way 'un'terleaving: ****/\n\t"\
	"/* a[j+p0]: Inputs from r1 +0/1, 16/17, 2/3, 18/19. Outputs into [add0,add0+0x100,add1,add1+0x100]+0x0: */	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovapd	     (%%rax),%%ymm1						\n\t		vmovapd	0x020(%%rax),%%ymm3								\n\t"\
		"vmovapd	0x200(%%rax),%%ymm5						\n\t		vmovapd	0x220(%%rax),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovapd	0x040(%%rax),%%ymm6						\n\t		vmovapd	0x060(%%rax),%%ymm14							\n\t"\
		"vmovapd	0x240(%%rax),%%ymm5						\n\t		vmovapd	0x260(%%rax),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm6 ,%%ymm0				\n\t		vshufpd	$15 ,%%ymm13,%%ymm14 ,%%ymm2					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm6 ,%%ymm6				\n\t		vshufpd	$0  ,%%ymm13,%%ymm14 ,%%ymm14					\n\t"\
		"vperm2f128 $32 ,%%ymm0 ,%%ymm7 ,%%ymm5				\n\t		vperm2f128 $32 ,%%ymm2 ,%%ymm15,%%ymm13		/* outB	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm0 ,%%ymm7 ,%%ymm7				\n\t		vperm2f128 $49 ,%%ymm2 ,%%ymm15,%%ymm15		/* outD	*/	\n\t"\
		"vperm2f128 $32 ,%%ymm6 ,%%ymm1 ,%%ymm0				\n\t		vperm2f128 $32 ,%%ymm14,%%ymm3 ,%%ymm2 		/* outA	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm6 ,%%ymm1 ,%%ymm1				\n\t		vperm2f128 $49 ,%%ymm14,%%ymm3 ,%%ymm3 		/* outC	*/	\n\t"\
		"vmovapd 	%%ymm5 ,0x100(%%rbx)					\n\t		vmovapd %%ymm13,0x120(%%rbx)				/* outB	*/	\n\t"\
		"vmovapd 	%%ymm7 ,0x300(%%rbx)					\n\t		vmovapd %%ymm15,0x320(%%rbx)				/* outD	*/	\n\t"\
		"vmovapd 	%%ymm0 ,     (%%rbx)					\n\t		vmovapd %%ymm2 ,0x020(%%rbx)				/* outA	*/	\n\t"\
		"vmovapd 	%%ymm1 ,0x200(%%rbx)					\n\t		vmovapd %%ymm3 ,0x220(%%rbx)				/* outC	*/	\n\t"\
	"/* a[j+p2]: Inputs from r3 +0/1, 16/17, 2/3, 18/19. Outputs into [add0,add0+0x100,add1,add1+0x100]+0x40: */		\n\t"\
		"addq	$0x80,%%rax	\n\t"\
		"addq	$0x40,%%rbx	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovapd	     (%%rax),%%ymm1						\n\t		vmovapd	0x020(%%rax),%%ymm3								\n\t"\
		"vmovapd	0x200(%%rax),%%ymm5						\n\t		vmovapd	0x220(%%rax),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovapd	0x040(%%rax),%%ymm6						\n\t		vmovapd	0x060(%%rax),%%ymm14							\n\t"\
		"vmovapd	0x240(%%rax),%%ymm5						\n\t		vmovapd	0x260(%%rax),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm6 ,%%ymm0				\n\t		vshufpd	$15 ,%%ymm13,%%ymm14 ,%%ymm2					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm6 ,%%ymm6				\n\t		vshufpd	$0  ,%%ymm13,%%ymm14 ,%%ymm14					\n\t"\
		"vperm2f128 $32 ,%%ymm0 ,%%ymm7 ,%%ymm5				\n\t		vperm2f128 $32 ,%%ymm2 ,%%ymm15,%%ymm13		/* outB	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm0 ,%%ymm7 ,%%ymm7				\n\t		vperm2f128 $49 ,%%ymm2 ,%%ymm15,%%ymm15		/* outD	*/	\n\t"\
		"vperm2f128 $32 ,%%ymm6 ,%%ymm1 ,%%ymm0				\n\t		vperm2f128 $32 ,%%ymm14,%%ymm3 ,%%ymm2 		/* outA	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm6 ,%%ymm1 ,%%ymm1				\n\t		vperm2f128 $49 ,%%ymm14,%%ymm3 ,%%ymm3 		/* outC	*/	\n\t"\
		"vmovapd 	%%ymm5 ,0x100(%%rbx)					\n\t		vmovapd %%ymm13,0x120(%%rbx)				/* outB	*/	\n\t"\
		"vmovapd 	%%ymm7 ,0x300(%%rbx)					\n\t		vmovapd %%ymm15,0x320(%%rbx)				/* outD	*/	\n\t"\
		"vmovapd 	%%ymm0 ,     (%%rbx)					\n\t		vmovapd %%ymm2 ,0x020(%%rbx)				/* outA	*/	\n\t"\
		"vmovapd 	%%ymm1 ,0x200(%%rbx)					\n\t		vmovapd %%ymm3 ,0x220(%%rbx)				/* outC	*/	\n\t"\
	"prefetcht1	0x380(%%r14)\n\t"\
	"/* a[j+p4]: Inputs from r5 +0/1, 16/17, 2/3, 18/19. Outputs into [add0,add0+0x100,add1,add1+0x100]+0x80: */		\n\t"\
		"addq	$0x80,%%rax	\n\t"\
		"addq	$0x40,%%rbx	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovapd	     (%%rax),%%ymm1						\n\t		vmovapd	0x020(%%rax),%%ymm3								\n\t"\
		"vmovapd	0x200(%%rax),%%ymm5						\n\t		vmovapd	0x220(%%rax),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovapd	0x040(%%rax),%%ymm6						\n\t		vmovapd	0x060(%%rax),%%ymm14							\n\t"\
		"vmovapd	0x240(%%rax),%%ymm5						\n\t		vmovapd	0x260(%%rax),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm6 ,%%ymm0				\n\t		vshufpd	$15 ,%%ymm13,%%ymm14 ,%%ymm2					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm6 ,%%ymm6				\n\t		vshufpd	$0  ,%%ymm13,%%ymm14 ,%%ymm14					\n\t"\
		"vperm2f128 $32 ,%%ymm0 ,%%ymm7 ,%%ymm5				\n\t		vperm2f128 $32 ,%%ymm2 ,%%ymm15,%%ymm13		/* outB	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm0 ,%%ymm7 ,%%ymm7				\n\t		vperm2f128 $49 ,%%ymm2 ,%%ymm15,%%ymm15		/* outD	*/	\n\t"\
		"vperm2f128 $32 ,%%ymm6 ,%%ymm1 ,%%ymm0				\n\t		vperm2f128 $32 ,%%ymm14,%%ymm3 ,%%ymm2 		/* outA	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm6 ,%%ymm1 ,%%ymm1				\n\t		vperm2f128 $49 ,%%ymm14,%%ymm3 ,%%ymm3 		/* outC	*/	\n\t"\
		"vmovapd 	%%ymm5 ,0x100(%%rbx)					\n\t		vmovapd %%ymm13,0x120(%%rbx)				/* outB	*/	\n\t"\
		"vmovapd 	%%ymm7 ,0x300(%%rbx)					\n\t		vmovapd %%ymm15,0x320(%%rbx)				/* outD	*/	\n\t"\
		"vmovapd 	%%ymm0 ,     (%%rbx)					\n\t		vmovapd %%ymm2 ,0x020(%%rbx)				/* outA	*/	\n\t"\
		"vmovapd 	%%ymm1 ,0x200(%%rbx)					\n\t		vmovapd %%ymm3 ,0x220(%%rbx)				/* outC	*/	\n\t"\
	"/* a[j+p6]: Inputs from r7 +0/1, 16/17, 2/3, 18/19. Outputs into [add0,add0+0x100,add1,add1+0x100]+0xc0: */		\n\t"\
		"addq	$0x80,%%rax	\n\t"\
		"addq	$0x40,%%rbx	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovapd	     (%%rax),%%ymm1						\n\t		vmovapd	0x020(%%rax),%%ymm3								\n\t"\
		"vmovapd	0x200(%%rax),%%ymm5						\n\t		vmovapd	0x220(%%rax),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovapd	0x040(%%rax),%%ymm6						\n\t		vmovapd	0x060(%%rax),%%ymm14							\n\t"\
		"vmovapd	0x240(%%rax),%%ymm5						\n\t		vmovapd	0x260(%%rax),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm6 ,%%ymm0				\n\t		vshufpd	$15 ,%%ymm13,%%ymm14 ,%%ymm2					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm6 ,%%ymm6				\n\t		vshufpd	$0  ,%%ymm13,%%ymm14 ,%%ymm14					\n\t"\
		"vperm2f128 $32 ,%%ymm0 ,%%ymm7 ,%%ymm5				\n\t		vperm2f128 $32 ,%%ymm2 ,%%ymm15,%%ymm13		/* outB	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm0 ,%%ymm7 ,%%ymm7				\n\t		vperm2f128 $49 ,%%ymm2 ,%%ymm15,%%ymm15		/* outD	*/	\n\t"\
		"vperm2f128 $32 ,%%ymm6 ,%%ymm1 ,%%ymm0				\n\t		vperm2f128 $32 ,%%ymm14,%%ymm3 ,%%ymm2 		/* outA	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm6 ,%%ymm1 ,%%ymm1				\n\t		vperm2f128 $49 ,%%ymm14,%%ymm3 ,%%ymm3 		/* outC	*/	\n\t"\
		"vmovapd 	%%ymm5 ,0x100(%%rbx)					\n\t		vmovapd %%ymm13,0x120(%%rbx)				/* outB	*/	\n\t"\
		"vmovapd 	%%ymm7 ,0x300(%%rbx)					\n\t		vmovapd %%ymm15,0x320(%%rbx)				/* outD	*/	\n\t"\
		"vmovapd 	%%ymm0 ,     (%%rbx)					\n\t		vmovapd %%ymm2 ,0x020(%%rbx)				/* outA	*/	\n\t"\
		"vmovapd 	%%ymm1 ,0x200(%%rbx)					\n\t		vmovapd %%ymm3 ,0x220(%%rbx)				/* outC	*/	\n\t"\
	"prefetcht1	0x380(%%r15)\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__add1] "m" (Xadd1)\
		 ,[__r1] "m" (Xr1)\
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__pfetch_dist] "m" (Xpfetch_dist)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r14","r15","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
	);\
	}

	#define FOO(Xadd0,Xadd1,Xr1,Xisrt2)\
	{\
	__asm__ volatile (\
		"movq	%[__add1],%%rcx				\n\t							vmulpd		0x0c0(%%rdx),%%ymm12,%%ymm12	/* c2 */		\n\t"\
		"/*vmovaps	%%ymm3,%%ymm7				*/\n\t						vmulpd		0x0c0(%%rdx),%%ymm13,%%ymm13		\n\t"\
		"/*vmovaps	%%ymm2,%%ymm6				*/\n\t						vmulpd		0x0e0(%%rdx),%%ymm11,%%ymm11		\n\t"\
		"/*unpckhpd	0x120(%%rcx),%%ymm7			*/\n\t						vmulpd		0x0e0(%%rdx),%%ymm9 ,%%ymm9 		\n\t"\
		"/*unpcklpd	0x120(%%rcx),%%ymm3			*/\n\t						vsubpd		%%ymm11     ,%%ymm13,%%ymm13		\n\t"\
		"/*vmovaps	%%ymm7,0x120(%%rcx)			*/\n\t						vaddpd		%%ymm9      ,%%ymm12,%%ymm12		\n\t"\
		"/*addq	$0x040,%%rdx	 		c4 in lcol; c10 in rcol */\n\t		/*vmovaps		%%ymm13,%%ymm11				*/\n\t"\
		"/*unpckhpd	0x100(%%rcx),%%ymm6			*/\n\t						/*vmovaps		%%ymm12,%%ymm9 				*/\n\t"\
		"/*unpcklpd	0x100(%%rcx),%%ymm2			*/\n\t						/*unpckhpd	0x060(%%rcx),%%ymm11		*/\n\t"\
		"/*vmovaps	%%ymm6,0x100(%%rcx)			*/\n\t						/*unpcklpd	0x060(%%rcx),%%ymm13		*/\n\t"\
		"/*vmovaps	%%ymm3,0x120(%%rbx)			*/\n\t						/*vmovaps		%%ymm11,0x060(%%rcx)		*/\n\t"\
		"/*vmovaps	%%ymm2,0x100(%%rbx)			*/\n\t						/*unpckhpd	0x040(%%rcx),%%ymm9 		*/\n\t"\
"\n\t"/* AVX debug: everywhere we have an unpckhpd/unpcklpd pair followed by a pair of writes-to-memory, copy the same data into local memory: */\
"vmovaps	0x120(%%rcx),%%ymm7		\n\t"\
"vmovaps	0x100(%%rcx),%%ymm6		\n\t"\
"vmovaps	%%ymm3,0x120(%%rax)	/* r9 */\n\t"\
"vmovaps	%%ymm2,0x100(%%rax)	/* r8 */\n\t"\
"vmovaps	%%ymm7,0x320(%%rax)	/* r25 */\n\t"\
"vmovaps	%%ymm6,0x300(%%rax)	/* r24 */\n\t"\
		"/*vmovaps	0x020(%%rbx),%%ymm3			*/\n\t						/*unpcklpd	0x040(%%rcx),%%ymm12		*/\n\t"\
		"/*vmovaps	     (%%rbx),%%ymm2			*/\n\t						/*vmovaps		%%ymm9 ,0x040(%%rcx)		*/\n\t"\
		"/*vmovaps	%%ymm3,%%ymm7				*/\n\t						/*vmovaps		%%ymm13,0x060(%%rbx)		*/\n\t"\
		"/*vmovaps	%%ymm2,%%ymm6				*/\n\t						/*vmovaps		%%ymm12,0x040(%%rbx)		*/\n\t"\
																	"vmovaps	0x060(%%rcx),%%ymm11	\n\t"\
																	"vmovaps	0x040(%%rcx),%%ymm9		\n\t"\
																	"vmovaps	%%ymm13,0x060(%%rax)	/* r3 */\n\t"\
																	"vmovaps	%%ymm12,0x040(%%rax)	/* r2 */\n\t"\
																	"vmovaps	%%ymm11,0x260(%%rax)	/* r19 */\n\t"\
																	"vmovaps	%%ymm9 ,0x240(%%rax)	/* r18 */\n\t"\
		"																	addq	$0x040,%%rdx		/* c4 in lcol; c10 in rcol*/		\n\t"\
		"/*unpckhpd	0x020(%%rcx),%%ymm7			*/\n\t						vmovaps		     (%%rax),%%ymm12		\n\t"\
		"/*unpcklpd	0x020(%%rcx),%%ymm3			*/\n\t						vmovaps		0x020(%%rax),%%ymm13		\n\t"\
		"/*vmovaps	%%ymm7,0x020(%%rcx)			*/\n\t						vmovaps		%%ymm12,%%ymm11		\n\t"\
		"/*unpckhpd	(%%rcx),%%ymm6				*/\n\t						vmovaps		%%ymm13,%%ymm9 		\n\t"\
		"/*unpcklpd	(%%rcx),%%ymm2				*/\n\t						vmulpd		0x0c0(%%rdx),%%ymm12,%%ymm12		\n\t"\
		"/*vmovaps	%%ymm6,    (%%rcx)			*/\n\t						vmulpd		0x0c0(%%rdx),%%ymm13,%%ymm13		\n\t"\
		"/*vmovaps	%%ymm3,0x020(%%rbx)			*/\n\t						vmulpd		0x0e0(%%rdx),%%ymm11,%%ymm11		\n\t"\
		"/*vmovaps	%%ymm2,    (%%rbx)			*/\n\t						vmulpd		0x0e0(%%rdx),%%ymm9 ,%%ymm9 		\n\t"\
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
		"vmovaps	%%ymm0,%%ymm2				\n\t						/*vmovaps		%%ymm13,%%ymm11				*/\n\t"\
		"vmovaps	%%ymm1,%%ymm3				\n\t						/*vmovaps		%%ymm12,%%ymm9 				*/\n\t"\
		"vaddpd		%%ymm5,%%ymm5,%%ymm5		\n\t						/*unpckhpd	0x160(%%rcx),%%ymm11		*/\n\t"\
		"vaddpd		%%ymm4,%%ymm4,%%ymm4		\n\t						/*unpcklpd	0x160(%%rcx),%%ymm13		*/\n\t"\
		"vmovaps	%%ymm0,%%ymm6				\n\t						/*vmovaps		%%ymm11,0x160(%%rcx)		*/\n\t"\
		"vmovaps	%%ymm1,%%ymm7				\n\t						/*unpckhpd	0x140(%%rcx),%%ymm9 		*/\n\t"\
		"vmulpd		     (%%rdx),%%ymm2,%%ymm2	\n\t						/*unpcklpd	0x140(%%rcx),%%ymm12		*/\n\t"\
		"vmulpd		     (%%rdx),%%ymm3,%%ymm3	\n\t						/*vmovaps		%%ymm9 ,0x140(%%rcx)		*/\n\t"\
		"vmulpd		0x020(%%rdx),%%ymm6,%%ymm6	\n\t						/*vmovaps		%%ymm13,0x160(%%rbx)		*/\n\t"\
		"vmulpd		0x020(%%rdx),%%ymm7,%%ymm7	\n\t						/*vmovaps		%%ymm12,0x140(%%rbx)		*/\n\t"\
																	"vmovaps	0x160(%%rcx),%%ymm11	\n\t"\
																	"vmovaps	0x140(%%rcx),%%ymm9		\n\t"\
																	"vmovaps	%%ymm13,0x160(%%rax)	/* r11 */\n\t"\
																	"vmovaps	%%ymm12,0x140(%%rax)	/* r10 */\n\t"\
																	"vmovaps	%%ymm11,0x360(%%rax)	/* r27 */\n\t"\
																	"vmovaps	%%ymm9 ,0x340(%%rax)	/* r26 */\n\t"\
		"addq	$0x040,%%rdx		/* c12 in lcol; c6 in rcol*/\n\t		vsubpd		%%ymm15,%%ymm8 ,%%ymm8 		\n\t"\
		"vsubpd		%%ymm6,%%ymm3,%%ymm3		\n\t						vsubpd		%%ymm14,%%ymm10,%%ymm10		\n\t"\
		"vaddpd		%%ymm7,%%ymm2,%%ymm2		\n\t						vaddpd		%%ymm15,%%ymm15,%%ymm15		\n\t"\
		"/*vmovaps	%%ymm3,%%ymm7				*/\n\t						vaddpd		%%ymm14,%%ymm14,%%ymm14		\n\t"\
		"/*vmovaps	%%ymm2,%%ymm6				*/\n\t						vaddpd		%%ymm8 ,%%ymm15,%%ymm15		\n\t"\
		"/*unpckhpd	0x0a0(%%rcx),%%ymm7			*/\n\t						vaddpd		%%ymm10,%%ymm14,%%ymm14		\n\t"\
		"/*unpcklpd	0x0a0(%%rcx),%%ymm3			*/\n\t						vmovaps		%%ymm15,%%ymm12		\n\t"\
		"/*vmovaps	%%ymm7,0x0a0(%%rcx)			*/\n\t						vmovaps		%%ymm10,%%ymm13		\n\t"\
		"/*unpckhpd	0x080(%%rcx),%%ymm6			*/\n\t						vmulpd		0x0c0(%%rdx),%%ymm15,%%ymm15		\n\t"\
		"/*unpcklpd	0x080(%%rcx),%%ymm2			*/\n\t						vmulpd		0x0c0(%%rdx),%%ymm10,%%ymm10		\n\t"\
		"/*vmovaps	%%ymm6,0x080(%%rcx)			*/\n\t						vmulpd		0x0e0(%%rdx),%%ymm12,%%ymm12		\n\t"\
		"/*vmovaps	%%ymm3,0x0a0(%%rbx)			*/\n\t						vmulpd		0x0e0(%%rdx),%%ymm13,%%ymm13		\n\t"\
		"/*vmovaps	%%ymm2,0x080(%%rbx)			*/\n\t						vsubpd		%%ymm12     ,%%ymm10,%%ymm10		\n\t"\
"vmovaps	0x0a0(%%rcx),%%ymm7		\n\t"\
"vmovaps	0x080(%%rcx),%%ymm6		\n\t"\
"vmovaps	%%ymm3,0x0a0(%%rax)	/* r5 */\n\t"\
"vmovaps	%%ymm2,0x080(%%rax)	/* r4 */\n\t"\
"vmovaps	%%ymm7,0x2a0(%%rax)	/* r21 */\n\t"\
"vmovaps	%%ymm6,0x280(%%rax)	/* r20 */\n\t"\
		"vsubpd		%%ymm5,%%ymm0,%%ymm0		\n\t						vaddpd		%%ymm13     ,%%ymm15,%%ymm15		\n\t"\
		"vaddpd		%%ymm4,%%ymm1,%%ymm1		\n\t						/*vmovaps		%%ymm10,%%ymm13			*/\n\t"\
		"vmovaps	%%ymm0,%%ymm6				\n\t						/*vmovaps		%%ymm15,%%ymm12			*/\n\t"\
		"vmovaps	%%ymm1,%%ymm7				\n\t						/*unpckhpd	0x0e0(%%rcx),%%ymm13		*/\n\t"\
		"vmulpd		     (%%rdx),%%ymm0,%%ymm0	\n\t						/*unpcklpd	0x0e0(%%rcx),%%ymm10		*/\n\t"\
		"vmulpd		     (%%rdx),%%ymm1,%%ymm1	\n\t						/*vmovaps		%%ymm13,0x0e0(%%rcx)	*/\n\t"\
		"vmulpd		0x020(%%rdx),%%ymm6,%%ymm6	\n\t						/*unpckhpd	0x0c0(%%rcx),%%ymm12		*/\n\t"\
		"vmulpd		0x020(%%rdx),%%ymm7,%%ymm7	\n\t						/*unpcklpd	0x0c0(%%rcx),%%ymm15		*/\n\t"\
		"vsubpd		%%ymm6,%%ymm1,%%ymm1		\n\t						/*vmovaps		%%ymm12,0x0c0(%%rcx)	*/\n\t"\
		"vaddpd		%%ymm7,%%ymm0,%%ymm0		\n\t						/*vmovaps		%%ymm10,0x0e0(%%rbx)	*/\n\t"\
		"/*vmovaps	%%ymm1,%%ymm7				*/\n\t						/*vmovaps		%%ymm15,0x0c0(%%rbx)	*/\n\t"\
																	"vmovaps	0x0e0(%%rcx),%%ymm13	\n\t"\
																	"vmovaps	0x0c0(%%rcx),%%ymm12	\n\t"\
																	"vmovaps	%%ymm10,0x0e0(%%rax)	/* r7 */\n\t"\
																	"vmovaps	%%ymm15,0x0c0(%%rax)	/* r6 */\n\t"\
																	"vmovaps	%%ymm13,0x2e0(%%rax)	/* r23 */\n\t"\
																	"vmovaps	%%ymm12,0x2c0(%%rax)	/* r22 */\n\t"\
		"/*vmovaps	%%ymm0,%%ymm6				*/\n\t						vmovaps		%%ymm8 ,%%ymm12		\n\t"\
		"/*unpckhpd	0x1a0(%%rcx),%%ymm7			*/\n\t						vmovaps		%%ymm14,%%ymm13		\n\t"\
		"/*unpcklpd	0x1a0(%%rcx),%%ymm1			*/\n\t						vmulpd		0x100(%%rdx),%%ymm8 ,%%ymm8 	/* c14 */		\n\t"\
		"/*vmovaps	%%ymm7,0x1a0(%%rcx)			*/\n\t						vmulpd		0x100(%%rdx),%%ymm14,%%ymm14		\n\t"\
		"/*unpckhpd	0x180(%%rcx),%%ymm6			*/\n\t						vmulpd		0x120(%%rdx),%%ymm12,%%ymm12		\n\t"\
		"/*unpcklpd	0x180(%%rcx),%%ymm0			*/\n\t						vmulpd		0x120(%%rdx),%%ymm13,%%ymm13		\n\t"\
		"/*vmovaps	%%ymm6,0x180(%%rcx)			*/\n\t						vsubpd		%%ymm12     ,%%ymm14,%%ymm14		\n\t"\
		"/*vmovaps	%%ymm1,0x1a0(%%rbx)			*/\n\t						vaddpd		%%ymm13     ,%%ymm8 ,%%ymm8 		\n\t"\
		"/*vmovaps	%%ymm0,0x180(%%rbx)			*/\n\t						/*vmovaps		%%ymm14,%%ymm13				*/\n\t"\
"vmovaps	0x1a0(%%rcx),%%ymm7		\n\t"\
"vmovaps	0x180(%%rcx),%%ymm6		\n\t"\
"vmovaps	%%ymm1,0x1a0(%%rax)	/* r13 */\n\t"\
"vmovaps	%%ymm0,0x180(%%rax)	/* r12 */\n\t"\
"vmovaps	%%ymm7,0x3a0(%%rax)	/* r29 */\n\t"\
"vmovaps	%%ymm6,0x380(%%rax)	/* r28 */\n\t"\
		"																	/*vmovaps		%%ymm8 ,%%ymm12				*/\n\t"\
		"																	/*unpckhpd	0x1e0(%%rcx),%%ymm13		*/\n\t"\
		"																	/*unpcklpd	0x1e0(%%rcx),%%ymm14		*/\n\t"\
		"																	/*vmovaps		%%ymm13,0x1e0(%%rcx)		*/\n\t"\
		"																	/*unpckhpd	0x1c0(%%rcx),%%ymm12		*/\n\t"\
		"																	/*unpcklpd	0x1c0(%%rcx),%%ymm8 		*/\n\t"\
		"																	/*vmovaps		%%ymm12,0x1c0(%%rcx)		*/\n\t"\
		"																	/*vmovaps		%%ymm14,0x1e0(%%rbx)		*/\n\t"\
		"																	/*vmovaps		%%ymm8 ,0x1c0(%%rbx)		*/\n\t"\
																	"vmovaps	0x1e0(%%rcx),%%ymm13	\n\t"\
																	"vmovaps	0x1c0(%%rcx),%%ymm12	\n\t"\
																	"vmovaps	%%ymm14,0x1e0(%%rax)	/* r15 */\n\t"\
																	"vmovaps	%%ymm8 ,0x1c0(%%rax)	/* r14 */\n\t"\
																	"vmovaps	%%ymm13,0x3e0(%%rax)	/* r31 */\n\t"\
																	"vmovaps	%%ymm12,0x3c0(%%rax)	/* r30 */\n\t"\
		"\n\t"\
		"movq	%[__isrt2],%%rbx	/* DEBUG: Dump register contents into trig-data slots: */\n\t"\
		"vmovaps	%%ymm0 ,0x600(%%rbx)	\n\t"\
		"vmovaps	%%ymm1 ,0x620(%%rbx)	\n\t"\
		"vmovaps	%%ymm2 ,0x640(%%rbx)	\n\t"\
		"vmovaps	%%ymm3 ,0x660(%%rbx)	\n\t"\
		"vmovaps	%%ymm4 ,0x680(%%rbx)	\n\t"\
		"vmovaps	%%ymm5 ,0x6a0(%%rbx)	\n\t"\
		"vmovaps	%%ymm6 ,0x6c0(%%rbx)	\n\t"\
		"vmovaps	%%ymm7 ,0x6e0(%%rbx)	\n\t"\
		"vmovaps	%%ymm8 ,0x700(%%rbx)	\n\t"\
		"vmovaps	%%ymm9 ,0x720(%%rbx)	\n\t"\
		"vmovaps	%%ymm10,0x740(%%rbx)	\n\t"\
		"vmovaps	%%ymm11,0x760(%%rbx)	\n\t"\
		"vmovaps	%%ymm12,0x780(%%rbx)	\n\t"\
		"vmovaps	%%ymm13,0x7a0(%%rbx)	\n\t"\
		"vmovaps	%%ymm14,0x7c0(%%rbx)	\n\t"\
		"vmovaps	%%ymm15,0x7e0(%%rbx)	\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__add1] "m" (Xadd1)\
		 ,[__r1] "m" (Xr1)\
		 ,[__isrt2] "m" (Xisrt2)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
	);\
	}

#elif defined(USE_SSE2)

	// Fused radix-16-DIF/dyadic-square/radix-16-DIT pass enabled by the FULLY_FUSED flag in radix16_dyadic_square.c:
	#define SSE2_RADIX16_DIF_DYADIC_DIT(Xadd0,Xadd1,Xr1,Xisrt2,Xpfetch_dist)\
	{\
	__asm__ volatile (\
	"/*************************************************************/\n\t"\
	"/* SSE2_RADIX16_WRAPPER_DIF, 1st set of inputs:              */\n\t"\
	"/*************************************************************/\n\t"\
		"movq	%[__r1] ,%%rcx	\n\t"\
	"/*...Block 1: Do the p4,12 combo: */\n\t										/*...Block 2:		Cost: 46 MOVapd, 16 UNPCKHPD, 28 ADD/SUBpd, 16 MULpd */\n\t"\
		"movq	%[__add0],%%rax\n\t														/* Do the p2,10 combo: */\n\t"\
		"movq	%[__add1],%%rbx\n\t"\
	"movslq	%[__pfetch_dist],%%r15	\n\t"\
	"leaq	(%%rax,%%r15,8),%%r14	\n\t"	/* Block 1 [base-address + data-fetch-ahead index] */\
	"leaq	(%%rbx,%%r15,8),%%r15	\n\t"	/* Block 2 [base-address + data-fetch-ahead index] */\
	"prefetcht1	(%%r14)\n\t"\
		"leaq	0x2b0(%%rcx),%%rdi	/* c2 */\n\t										movaps		0x20(%%rax)	,%%xmm14		/* a[j1+p2 ] */\n\t"\
		"leaq	0x270(%%rcx),%%rdx	/* c4 */\n\t										movaps		0x20(%%rbx)	,%%xmm13		/* a[j2+p2 ], xmm13 used for scratch */\n\t"\
		"\n\t																			movaps		%%xmm14		,%%xmm8 		\n\t"\
	"/* Real parts: */\n\t																unpckhpd	%%xmm13		,%%xmm14		\n\t"\
		"movaps		0x40(%%rax)	,%%xmm6		/* a[j1+p4] */\n\t							unpcklpd	%%xmm13		,%%xmm8 		\n\t"\
		"movaps		0x40(%%rbx)	,%%xmm5		/* a[j2+p4] */\n\t							movaps			%%xmm14	,0x180(%%rcx)	/* Store hi real in t9 +16 */\n\t"\
		"movaps		%%xmm6		,%%xmm0		/* a[j1+p4] copy */\n\t						movaps		0x30(%%rax)	,%%xmm15		\n\t"\
		"unpckhpd	%%xmm5,%%xmm6		\n\t											movaps		0x30(%%rbx)	,%%xmm14		\n\t"\
		"unpcklpd	%%xmm5,%%xmm0		\n\t											movaps		%%xmm15		,%%xmm9 		\n\t"\
		"movaps		%%xmm6		,0x140(%%rcx)	/* Store hi real in t21 */\n\t			unpckhpd	%%xmm14		,%%xmm15		\n\t"\
	"/* Imag parts: */\n\t																unpcklpd	%%xmm14		,%%xmm9 		\n\t"\
		"movaps		0x50(%%rax)	,%%xmm7	\n\t											movaps		%%xmm15		,0x190(%%rcx)	/* Store hi imag in t10+16 */\n\t"\
		"movaps		0x50(%%rbx)	,%%xmm5	\n\t											movaps		%%xmm8 		,%%xmm10	/* xmm10 <- cpy a[jt+p2] */\n\t"\
		"movaps		%%xmm7,%%xmm1		\n\t											movaps		%%xmm9 		,%%xmm11	/* xmm11 <- cpy a[jp+p2] */\n\t"\
		"unpckhpd	%%xmm5,%%xmm7		\n\t											mulpd		    (%%rdi)	,%%xmm8 		/* a[jt+p2]*c2 */\n\t"\
		"unpcklpd	%%xmm5,%%xmm1		\n\t											mulpd		    (%%rdi)	,%%xmm9 		/* a[jp+p2]*c2 */\n\t"\
		"movaps		%%xmm7		,0x150(%%rcx)	/* Store hi imag in t22 */\n\t			mulpd		0x10(%%rdi)	,%%xmm10		/* a[jt+p2]*s2 */\n\t"\
		"\n\t																			mulpd		0x10(%%rdi)	,%%xmm11		/* a[jp+p2]*s2 */\n\t"\
		"movaps		%%xmm0		,%%xmm2	/* xmm2 <- cpy a[jt+p4] */\n\t					addpd		%%xmm10		,%%xmm9 	/* xmm9  <- t10*/\n\t"\
		"movaps		%%xmm1		,%%xmm3	/* xmm3 <- cpy a[jp+p4] */\n\t					subpd		%%xmm11		,%%xmm8 	/* xmm8  <- t9 */\n\t"\
	"/***************************************************************************/\n\t	movaps		%%xmm9 		,%%xmm11	/* xmm11 <- cpy t10*/\n\t"\
	"/*** From here on, things are identical to the code in radix16_dif_pass: ***/\n\t	movaps		%%xmm8 		,%%xmm10	/* xmm10 <- cpy t9 */\n\t"\
	"/***************************************************************************/\n\t	addq	$0x20,%%rdi	/* c10 */	\n\t"\
		"mulpd		    (%%rdx)	,%%xmm0		/* a[jt+p4]*c4 */\n\t						movaps		0xa0(%%rax)	,%%xmm14		/* a[j1+p10] */\n\t"\
		"mulpd		    (%%rdx)	,%%xmm1		/* a[jp+p4]*c4 */\n\t						movaps		0xa0(%%rbx)	,%%xmm13		/* a[j2+p10], xmm13 used for scratch space */\n\t"\
		"mulpd		0x10(%%rdx)	,%%xmm2		/* a[jt+p4]*s4 */\n\t						movaps		%%xmm14		,%%xmm12		\n\t"\
		"mulpd		0x10(%%rdx)	,%%xmm3		/* a[jp+p4]*s4 */\n\t						unpckhpd	%%xmm13		,%%xmm14		\n\t"\
		"addpd		%%xmm2		,%%xmm1	/* xmm1 <- t6 */\n\t							unpcklpd	%%xmm13		,%%xmm12		\n\t"\
		"subpd		%%xmm3		,%%xmm0	/* xmm0 <- t5 */\n\t							movaps		%%xmm14		,0x1a0(%%rcx)	/* Store hi real in t11+16 */\n\t"\
		"movaps		%%xmm1		,%%xmm3	/* xmm3 <- cpy t6 */\n\t						\n\t"\
		"movaps		%%xmm0		,%%xmm2	/* xmm2 <- cpy t5 */\n\t						movaps		0xb0(%%rax)	,%%xmm15		\n\t"\
		"\n\t																			movaps		0xb0(%%rbx)	,%%xmm14		/* Here xmm14 is the tmp-register */\n\t"\
		"addq		$0x20,%%rdx	/* c12 */	\n\t										movaps		%%xmm15		,%%xmm13		\n\t"\
	"/* Real parts: */\n\t																unpckhpd	%%xmm14		,%%xmm15		\n\t"\
		"movaps		0xc0(%%rax)	,%%xmm6		/* a[j1+p12] */\n\t							unpcklpd	%%xmm14		,%%xmm13		\n\t"\
		"movaps		0xc0(%%rbx)	,%%xmm5		/* a[j2+p12], xmm5 used for scratch */\n\t	movaps			%%xmm15	,0x1b0(%%rcx)	/* Store hi imag in t12+16 */\n\t"\
		"movaps		%%xmm6		,%%xmm4		\n\t										movaps			%%xmm12	,%%xmm14	/* xmm14 <- cpy a[jt+p10] */\n\t"\
		"unpckhpd	%%xmm5		,%%xmm6		\n\t										movaps			%%xmm13	,%%xmm15	/* xmm15 <- cpy a[jp+p10] */\n\t"\
		"unpcklpd	%%xmm5		,%%xmm4		\n\t										mulpd		    (%%rdi)	,%%xmm12		/* a[jt+p10]*c10 */\n\t"\
		"movaps		%%xmm6		,0x160(%%rcx)	/* Store hi real in t23 */\n\t			mulpd		    (%%rdi)	,%%xmm13		/* a[jp+p10]*c10 */\n\t"\
	"prefetcht1	(%%r15)\n\t"\
	"/* Imag parts: */\n\t																mulpd		0x10(%%rdi)	,%%xmm14		/* a[jt+p10]*s10 */\n\t"\
		"movaps		0xd0(%%rax)	,%%xmm7		\n\t										mulpd		0x10(%%rdi)	,%%xmm15		/* a[jp+p10]*s10 */\n\t"\
		"movaps		0xd0(%%rbx)	,%%xmm6		/* Here xmm6 is the tmp-register */\n\t		addpd		%%xmm14		,%%xmm13	/* xmm13 <- it */\n\t"\
		"movaps		%%xmm7		,%%xmm5		\n\t										subpd		%%xmm15		,%%xmm12	/* xmm12 <- rt		xmm14,7 free */\n\t"\
		"unpckhpd	%%xmm6		,%%xmm7		\n\t										addpd		%%xmm12		,%%xmm8 	/* ~t13<- t13+rt */\n\t"\
		"unpcklpd	%%xmm6		,%%xmm5		\n\t										addpd		%%xmm13		,%%xmm9 	/* ~t14<- t14+it */\n\t"\
		"movaps		%%xmm7		,0x170(%%rcx)	/* Store hi imag in t24 */\n\t			subpd		%%xmm12		,%%xmm10	/* ~t15<- t13-rt */\n\t"\
		"movaps		%%xmm4		,%%xmm6	/* xmm6 <- cpy a[jt+p12] */\n\t					subpd		%%xmm13		,%%xmm11	/* ~t16<- t14-it	xmm12,5 free */\n\t"\
		"movaps		%%xmm5		,%%xmm7	/* xmm7 <- cpy a[jp+p12] */\n\t					\n\t"\
		"mulpd		    (%%rdx)	,%%xmm4		/* a[jt+p12]*c12 */\n\t						/* Do the p6,14 combo - do p14 first so registers come out in same order as for p2,10 */\n\t"\
		"mulpd		    (%%rdx)	,%%xmm5		/* a[jp+p12]*c12 */\n\t						leaq	0x310(%%rcx),%%rdi	/* c14 */\n\t"\
		"mulpd		0x10(%%rdx)	,%%xmm6		/* a[jt+p12]*s12 */\n\t						\n\t"\
		"mulpd		0x10(%%rdx)	,%%xmm7		/* a[jp+p12]*s12 */\n\t						movaps		0xe0(%%rax)	,%%xmm14		/* a[j1+p14] */\n\t"\
		"addpd		%%xmm6		,%%xmm5	/* xmm5 <- it */\n\t							movaps		0xe0(%%rbx)	,%%xmm13		/* a[j2+p14], xmm13 used for scratch space */\n\t"\
		"subpd		%%xmm7		,%%xmm4	/* xmm4 <- rt		xmm6,7 free */\n\t			movaps		%%xmm14		,%%xmm12		\n\t"\
		"addpd		%%xmm4		,%%xmm0	/* ~t5 <- t5 +rt */\n\t							unpckhpd	%%xmm13		,%%xmm14		\n\t"\
		"addpd		%%xmm5		,%%xmm1	/* ~t6 <- t6 +it */\n\t							unpcklpd	%%xmm13		,%%xmm12		\n\t"\
		"subpd		%%xmm4		,%%xmm2	/* ~t7 <- t5 -rt */\n\t							movaps		%%xmm14		,0x1e0(%%rcx)	/* Store hi real in t15+16 */\n\t"\
		"subpd		%%xmm5		,%%xmm3	/* ~t8 <- t6 -it	xmm4,5 free */\n\t			\n\t"\
		"\n\t																			movaps		0xf0(%%rax)	,%%xmm15		\n\t"\
		"/* Now do the p0,8 combo: */\n\t												movaps		0xf0(%%rbx)	,%%xmm14		/* Here xmm14 is the tmp-register */\n\t"\
		"leaq	0x250(%%rcx),%%rdx	/* c8 */\n\t										movaps		%%xmm15		,%%xmm13		\n\t"\
	"/* Real parts: */\n\t																unpckhpd	%%xmm14		,%%xmm15		\n\t"\
		"movaps		0x80(%%rax)	,%%xmm6		/* a[j1+p8 ] */\n\t							unpcklpd	%%xmm14		,%%xmm13		\n\t"\
		"movaps		0x80(%%rbx)	,%%xmm5		/* a[j2+p8 ], xmm5 used for scratch */\n\t	movaps			%%xmm15	,0x1f0(%%rcx)	/* Store hi imag in t16+16 */\n\t"\
		"movaps		%%xmm6		,%%xmm4		\n\t										\n\t"\
		"unpckhpd	%%xmm5		,%%xmm6		\n\t										movaps			%%xmm12	,%%xmm14		/* xmm14 <- cpy a[jt+p14] */\n\t"\
		"unpcklpd	%%xmm5		,%%xmm4		\n\t										movaps			%%xmm13	,%%xmm15		/* xmm15 <- cpy a[jp+p14] */\n\t"\
		"movaps		%%xmm6		,0x120(%%rcx)	/* Store hi real in t19 */\n\t			\n\t"\
	"/* Imag parts: */\n\t																mulpd		    (%%rdi)	,%%xmm12		/* a[jt+p14]*c14 */\n\t"\
		"movaps		0x90(%%rax)	,%%xmm7		\n\t										mulpd		    (%%rdi)	,%%xmm13		/* a[jp+p14]*c14 */\n\t"\
		"movaps		0x90(%%rbx)	,%%xmm6		/* Here xmm6 is the tmp-register */\n\t		mulpd		0x10(%%rdi)	,%%xmm14		/* a[jt+p14]*s14 */\n\t"\
		"movaps		%%xmm7		,%%xmm5		\n\t										mulpd		0x10(%%rdi)	,%%xmm15		/* a[jp+p14]*s14 */\n\t"\
		"unpckhpd	%%xmm6		,%%xmm7		\n\t										addpd		%%xmm14		,%%xmm13		/* xmm13 <- it */\n\t"\
		"unpcklpd	%%xmm6		,%%xmm5		\n\t										subpd		%%xmm15		,%%xmm12		/* xmm12 <- rt		xmm14,7 free */\n\t"\
		"movaps		%%xmm7		,0x130(%%rcx)	/* Store hi imag in t20 */\n\t			movaps		%%xmm13		,0x0f0(%%rcx)	/* Store it in t16*/\n\t"\
		"\n\t																			movaps		%%xmm12		,0x0e0(%%rcx)	/* Store rt in t15*/\n\t"\
	"prefetcht1	0x40(%%r14)\n\t"\
		"movaps		%%xmm4		,%%xmm6	/* xmm6 <- cpy a[jt+p8] */\n\t					\n\t"\
		"movaps		%%xmm5		,%%xmm7	/* xmm7 <- cpy a[jp+p8] */\n\t					subq	$0x20,%%rdi	/* c6  */	\n\t"\
		"\n\t																			/* Real parts: */\n\t"\
		"mulpd		    (%%rdx)	,%%xmm4		/* a[jt+p8]*c8 */\n\t						movaps		0x60(%%rax)	,%%xmm14		/* a[j1+p6 ], this is the scratch register  */\n\t"\
		"mulpd		    (%%rdx)	,%%xmm5		/* a[jp+p8]*c8 */\n\t						movaps		0x60(%%rax)	,%%xmm12		/* a[j1+p6 ], this is the active  register */\n\t"\
		"mulpd		0x10(%%rdx)	,%%xmm6		/* a[jt+p8]*s8 */\n\t						unpckhpd	0x60(%%rbx)	,%%xmm14		/* a[j2+p6 ] gets read twice */\n\t"\
		"mulpd		0x10(%%rdx)	,%%xmm7		/* a[jp+p8]*s8 */\n\t						unpcklpd	0x60(%%rbx)	,%%xmm12		/* a[jt+p6 ] */\n\t"\
		"addpd		%%xmm6	,%%xmm5	/* xmm5 <- it */\n\t								movaps			%%xmm14	,0x1c0(%%rcx)	/* Store hi real in t13+16 */\n\t"\
		"subpd		%%xmm7	,%%xmm4	/* xmm4 <- rt;xmm6,7 free - store t1,t2 */\n\t	/* Imag parts: */\n\t"\
		"\n\t																			movaps		0x70(%%rax)	,%%xmm15\n\t"\
	"/* Real parts: */\n\t																movaps		0x70(%%rax)	,%%xmm13\n\t"\
		"movaps		    (%%rax)	,%%xmm6		/* a[j1    ], scratch register  */\n\t		unpckhpd	0x70(%%rbx)	,%%xmm15\n\t"\
		"movaps		    (%%rax)	,%%xmm7		/* a[j1    ], active  register */\n\t		unpcklpd	0x70(%%rbx)	,%%xmm13		/* a[jp+p6 ] */\n\t"\
		"unpckhpd	    (%%rbx)	,%%xmm6		/* a[j2    ] gets read twice */\n\t			movaps			%%xmm15	,0x1d0(%%rcx)	/* Store hi imag in t14+16 */\n\t"\
		"unpcklpd	    (%%rbx)	,%%xmm7		/* a[jt] = t1*/\n\t							\n\t"\
		"movaps		%%xmm6		,0x100(%%rcx)	/* Store hi real in t17 */\n\t			movaps			%%xmm12	,%%xmm14	/* xmm14 <- cpy a[jt+p6] */\n\t"\
		"movaps		%%xmm7		,     (%%rcx)	/* Store active  in t1  */\n\t			movaps			%%xmm13	,%%xmm15	/* xmm15 <- cpy a[jp+p6] */\n\t"\
	"/* Imag parts: */\n\t																\n\t"\
		"movaps		0x10(%%rax)	,%%xmm6\n\t												mulpd		    (%%rdi)	,%%xmm12		/* a[jt+p6]*c6 */\n\t"\
		"movaps		0x10(%%rax)	,%%xmm7\n\t												mulpd		    (%%rdi)	,%%xmm13		/* a[jp+p6]*c6 */\n\t"\
		"unpckhpd	0x10(%%rbx)	,%%xmm6\n\t												mulpd		0x10(%%rdi)	,%%xmm14		/* a[jt+p6]*s6 */\n\t"\
		"unpcklpd	0x10(%%rbx)	,%%xmm7		/* a[jp] = t2*/\n\t							mulpd		0x10(%%rdi)	,%%xmm15		/* a[jp+p6]*s6 */\n\t"\
		"movaps		%%xmm6		,0x110(%%rcx)	/* Store hi imag in t18... */\n\t		addpd		%%xmm14		,%%xmm13	/* xmm13 <- t14*/\n\t"\
		"movaps		    (%%rcx)	,%%xmm6		/* ...and reload t1. */\n\t					subpd		%%xmm15		,%%xmm12	/* xmm12 <- t13*/\n\t"\
		"\n\t																			movaps		%%xmm13		,%%xmm15	/* xmm15 <- cpy t14*/\n\t"\
		"subpd		%%xmm4		,%%xmm6	/* ~t3 <- t1 -rt */\n\t							movaps		%%xmm12		,%%xmm14	/* xmm14 <- cpy t13*/\n\t"\
		"subpd		%%xmm5		,%%xmm7	/* ~t4 <- t2 -it */\n\t							\n\t"\
		"addpd		%%xmm4		,%%xmm4	/*          2*rt */\n\t							subpd		0x0e0(%%rcx)	,%%xmm12		/* ~t15<- t13-rt */\n\t"\
		"addpd		%%xmm5		,%%xmm5	/*          2*it */\n\t							subpd		0x0f0(%%rcx)	,%%xmm13		/* ~t16<- t14-it */\n\t"\
		"addpd		%%xmm6		,%%xmm4	/* ~t1 <- t1 +rt */\n\t							addpd		0x0e0(%%rcx)	,%%xmm14		/* ~t13<- t13+rt */\n\t"\
		"addpd		%%xmm7		,%%xmm5	/* ~t2 <- t2 +it	xmm4,5 free */\n\t			addpd		0x0f0(%%rcx)	,%%xmm15		/* ~t14<- t14+it */\n\t"\
		"\n\t																			\n\t"\
	"/* Finish radix-4 butterfly and store results into temporary-array slots: */\n\t	/* Finish radix-4 butterfly and store results into temporary-array slots: */\n\t"\
		"subpd		%%xmm0		,%%xmm4	/*~t5 =t1 -t5 */\n\t							subpd		%%xmm14		,%%xmm8 	/*~t13*/\n\t"\
		"subpd		%%xmm1		,%%xmm5	/*~t6 =t2 -t6 */\n\t							subpd		%%xmm15		,%%xmm9 	/*~t14*/\n\t"\
		"movaps		%%xmm4		,0x040(%%rcx)	/* a[jt+p8 ] <- ~t5 */\n\t				movaps		%%xmm8 		,0x0c0(%%rcx)	/* a[jt+p8 ] <- ~t13*/\n\t"\
		"movaps		%%xmm5		,0x050(%%rcx)	/* a[jp+p8 ] <- ~t6 */\n\t				movaps		%%xmm9 		,0x0d0(%%rcx)	/* a[jp+p8 ] <- ~t14*/\n\t"\
		"addpd		%%xmm0		,%%xmm0			/* 2*t5 */\n\t							addpd		%%xmm14		,%%xmm14	/* 2*t13*/\n\t"\
		"addpd		%%xmm1		,%%xmm1			/* 2*t6 */\n\t							addpd		%%xmm15		,%%xmm15	/* 2*t14*/\n\t"\
		"addpd		%%xmm4		,%%xmm0			/*~t1 =t1 +t5 */\n\t					addpd		%%xmm8 		,%%xmm14	/*~t9 */\n\t"\
		"addpd		%%xmm5		,%%xmm1			/*~t2 =t2 +t6 */\n\t					addpd		%%xmm9 		,%%xmm15	/*~t10*/\n\t"\
		"movaps		%%xmm0		,     (%%rcx)	/* a[jt    ] <- ~t1 */\n\t				movaps		%%xmm14		,0x080(%%rcx)	/* a[jt    ] <- ~t9 */\n\t"\
		"movaps		%%xmm1		,0x010(%%rcx)	/* a[jp    ] <- ~t2 */\n\t				movaps		%%xmm15		,0x090(%%rcx)	/* a[jp    ] <- ~t10*/\n\t"\
	"prefetcht1	0x40(%%r15)\n\t"\
		"\n\t																			\n\t"\
		"subpd		%%xmm3		,%%xmm6	/*~t3 =t3 -t8 */\n\t							subpd		%%xmm13		,%%xmm10	/*~t11*/\n\t"\
		"subpd		%%xmm2		,%%xmm7	/*~t8 =t4 -t7 */\n\t							subpd		%%xmm12		,%%xmm11	/*~t16*/\n\t"\
		"movaps		%%xmm6		,0x020(%%rcx)	/* a[jt+p4 ] <- ~t3 */\n\t				movaps		%%xmm10		,0x0a0(%%rcx)	/* a[jt+p4 ] <- ~t11*/\n\t"\
		"movaps		%%xmm7		,0x070(%%rcx)	/* a[jp+p12] <- ~t8 */\n\t				movaps		%%xmm11		,0x0f0(%%rcx)	/* a[jp+p12] <- ~t16*/\n\t"\
		"addpd		%%xmm3		,%%xmm3			/* 2*t8 */\n\t							addpd		%%xmm13		,%%xmm13	/* 2*t16*/\n\t"\
		"addpd		%%xmm2		,%%xmm2			/* 2*t7 */\n\t							addpd		%%xmm12		,%%xmm12	/* 2*t15*/\n\t"\
		"addpd		%%xmm6		,%%xmm3			/*~t7 =t3 +t8 */\n\t					addpd		%%xmm10		,%%xmm13	/*~t15*/\n\t"\
		"addpd		%%xmm7		,%%xmm2			/*~t4 =t4 +t7 */\n\t					addpd		%%xmm11		,%%xmm12	/*~t12*/\n\t"\
		"movaps		%%xmm3		,0x060(%%rcx)	/* a[jt+p12] <- ~t7 */\n\t				movaps		%%xmm13		,0x0e0(%%rcx)	/* a[jt+p12] <- ~t15*/\n\t"\
		"movaps		%%xmm2		,0x030(%%rcx)	/* a[jp+p4 ] <- ~t4 */\n\t				movaps		%%xmm12		,0x0b0(%%rcx)	/* a[jp+p4 ] <- ~t12*/\n\t"\
		"\n\t"\
"/******************************************************************************************************************************/\n\t"\
"/* Each of the above 2 radix-4 blocks is unique, but the ensuing 2 blocks                                                     */\n\t"\
"/* [operating on the odd-indexed elements from the unpck*pd commands which were stored to temporaries can use a common macro: */\n\t"\
"/******************************************************************************************************************************/\n\t"\
		"\n\t"\
"/*...Block 3: */															\n\t/*...Block 4: */\n\t"\
		"\n\t																\n\t"\
"/*	SSE2_RADIX4_DIF_4TWIDDLE(r17,r21,r19,r23, r17,c1): */					\n\t/* SSE2_RADIX4_DIF_4TWIDDLE(r25,r29,r27,r31, r25,c3): */\n\t"\
	"/* Do the p0,p8 combo: */												\n\t	/* Do the p0,p8 combo: */		\n\t"\
		"leaq	0x100(%%rcx),%%rax	/* r17 */								\n\t		/* All __r and __c pointers incr by +0x80 in rcol w.r.to lcol: */\n\t"\
		"leaq	0x330(%%rcx),%%rbx	/* c1 */								\n\t		/* r25 */\n\t"\
		"addq	$0x120,%%rcx		/* r19 */								\n\t		/* c3  */\n\t"\
		"\n\t																\n\t		/* r27 */\n\t"\
		"movaps		    (%%rax)	,%%xmm0		/* a[jt   ] */					\n\t		movaps		0x80(%%rax)	,%%xmm8 		/* a[jt   ] */				\n\t"\
		"movaps		0x10(%%rax)	,%%xmm1		/* a[jp   ] */					\n\t		movaps		0x90(%%rax)	,%%xmm9 		/* a[jp   ] */				\n\t"\
		"movaps		    (%%rcx)	,%%xmm4		/* a[jt+p8 ] */					\n\t		movaps		0x80(%%rcx)	,%%xmm12		/* a[jt+p8 ] */\n\t"\
		"movaps		0x10(%%rcx)	,%%xmm5		/* a[jp+p8 ] */					\n\t		movaps		0x90(%%rcx)	,%%xmm13		/* a[jp+p8 ] */\n\t"\
		"movaps		    (%%rbx)	,%%xmm6		/* c0 */						\n\t		movaps		0x80(%%rbx)	,%%xmm14		/* c0 */\n\t"\
		"movaps		0x10(%%rbx)	,%%xmm7		/* s0 */						\n\t		movaps		0x90(%%rbx)	,%%xmm15		/* s0 */\n\t"\
		"movaps		%%xmm0		,%%xmm2		/* xmm2 <- cpy a[jt   ] */		\n\t		movaps		%%xmm8 		,%%xmm10		/* xmm10 <- cpy a[jt   ] */\n\t"\
		"movaps		%%xmm1		,%%xmm3		/* xmm3 <- cpy a[jp   ] */		\n\t		movaps		%%xmm9 		,%%xmm11		/* xmm11 <- cpy a[jp   ] */\n\t"\
		"\n\t																\n\t"\
		"mulpd   	%%xmm6		,%%xmm0		/* a[jt   ]*c0 */				\n\t		mulpd   	%%xmm14		,%%xmm8 		/* a[jt   ]*c0 */\n\t"\
		"mulpd   	%%xmm6		,%%xmm1		/* a[jp   ]*c0 */				\n\t		mulpd   	%%xmm14		,%%xmm9 		/* a[jp   ]*c0 */\n\t"\
		"mulpd   	%%xmm7		,%%xmm2		/* a[jt   ]*s0 */				\n\t		mulpd   	%%xmm15		,%%xmm10		/* a[jt   ]*s0 */\n\t"\
		"mulpd   	%%xmm7		,%%xmm3		/* a[jp   ]*s0	xmm6,7 free */	\n\t		mulpd   	%%xmm15		,%%xmm11		/* a[jp   ]*s0	xmm14,7 free */\n\t"\
		"movaps		%%xmm4		,%%xmm6		/* xmm6 <- cpy a[jt+p8 ] */		\n\t		movaps		%%xmm12		,%%xmm14		/* xmm14 <- cpy a[jt+p8 ] */\n\t"\
		"movaps		%%xmm5		,%%xmm7		/* xmm7 <- cpy a[jp+p8 ] */		\n\t		movaps		%%xmm13		,%%xmm15		/* xmm15 <- cpy a[jp+p8 ] */\n\t"\
		"addpd   	%%xmm2		,%%xmm1		/* xmm1 <- t2 */				\n\t		addpd   	%%xmm10		,%%xmm9 		/* xmm9  <- t2 */			\n\t"\
		"subpd   	%%xmm3		,%%xmm0		/* xmm0 <- t1 */				\n\t		subpd   	%%xmm11		,%%xmm8 		/* xmm8  <- t1 */			\n\t"\
	"prefetcht1	0x80(%%r14)\n\t"\
		"mulpd		0x20(%%rbx)	,%%xmm4		/* a[jt+p8 ]*c8 */				\n\t		mulpd		0xa0(%%rbx)	,%%xmm12		/* a[jt+p8 ]*c8 */\n\t"\
		"mulpd		0x20(%%rbx)	,%%xmm5		/* a[jp+p8 ]*c8 */				\n\t		mulpd		0xa0(%%rbx)	,%%xmm13		/* a[jp+p8 ]*c8 */\n\t"\
		"mulpd		0x30(%%rbx)	,%%xmm6		/* a[jt+p8 ]*s8 */				\n\t		mulpd		0xb0(%%rbx)	,%%xmm14		/* a[jt+p8 ]*s8 */\n\t"\
		"mulpd		0x30(%%rbx)	,%%xmm7		/* a[jp+p8 ]*s8 */				\n\t		mulpd		0xb0(%%rbx)	,%%xmm15		/* a[jp+p8 ]*s8 */\n\t"\
		"movaps		%%xmm0		,%%xmm2		/* xmm2 <- cpy t1 */			\n\t		movaps		%%xmm8 		,%%xmm10		/* xmm10 <- cpy t1 */		\n\t"\
		"movaps		%%xmm1		,%%xmm3		/* xmm3 <- cpy t2 */			\n\t		movaps		%%xmm9 		,%%xmm11		/* xmm11 <- cpy t2 */		\n\t"\
		"addpd		%%xmm6	    ,%%xmm5		/* xmm5 <- it */				\n\t		addpd		%%xmm14	    ,%%xmm13		/* xmm13 <- it */\n\t"\
		"subpd		%%xmm7		,%%xmm4		/* xmm4 <- rt    xmm6,7 free */	\n\t		subpd		%%xmm15		,%%xmm12		/* xmm12 <- rt    xmm14,7 free */\n\t"\
		"addq	$0x40 ,%%rcx	/* r23 */									\n\t		/* r31 */\n\t"\
		"addq	$0x60 ,%%rbx	/* c13 */									\n\t		/* c15 */\n\t"\
		"movaps		    (%%rcx)	,%%xmm6		/* a[jt+p12] */					\n\t		movaps	0x80(%%rcx),%%xmm14		/* a[jt+p12] */\n\t"\
		"movaps		0x10(%%rcx)	,%%xmm7		/* a[jp+p12] */					\n\t		movaps	0x90(%%rcx),%%xmm15		/* a[jp+p12] */\n\t"\
		"addpd		%%xmm4		,%%xmm0		/* ~t1 <- t1 +rt */				\n\t		addpd		%%xmm12		,%%xmm8 		/* ~t1 <- t1 +rt */\n\t"\
		"addpd		%%xmm5		,%%xmm1		/* ~t2 <- t2 +it */				\n\t		addpd		%%xmm13		,%%xmm9 		/* ~t2 <- t2 +it */\n\t"\
		"subpd		%%xmm4		,%%xmm2		/* ~t3 <- t1 -rt */				\n\t		subpd		%%xmm12		,%%xmm10		/* ~t3 <- t1 -rt */\n\t"\
		"subpd		%%xmm5		,%%xmm3		/* ~t4 <- t2 -it xmm4,5 free */	\n\t		subpd		%%xmm13		,%%xmm11		/* ~t4 <- t2 -it	xmm12,5 free */\n\t"\
		"\n\t																\n\t"\
	"/* Do the p4,12 combo: */												\n\t	/* Do the p4,12 combo: */\n\t"\
		"movaps		%%xmm6		,%%xmm4		/* xmm4 <- cpy a[jt+p12] */		\n\t		movaps		%%xmm14		,%%xmm12		/* xmm12 <- cpy a[jt+p12] */\n\t"\
		"movaps		%%xmm7		,%%xmm5		/* xmm5 <- cpy a[jp+p12] */		\n\t		movaps		%%xmm15		,%%xmm13		/* xmm13 <- cpy a[jp+p12] */\n\t"\
		"mulpd		    (%%rbx)	,%%xmm4		/* a[jt+p12]*c12 */				\n\t		mulpd		0x80(%%rbx)	,%%xmm12		/* a[jt+p12]*c12 */\n\t"\
		"mulpd		    (%%rbx)	,%%xmm5		/* a[jp+p12]*c12 */				\n\t		mulpd		0x80(%%rbx)	,%%xmm13		/* a[jp+p12]*c12 */\n\t"\
		"mulpd		0x10(%%rbx)	,%%xmm6		/* a[jt+p12]*s12 */				\n\t		mulpd		0x90(%%rbx)	,%%xmm14		/* a[jt+p12]*s12 */\n\t"\
		"mulpd		0x10(%%rbx)	,%%xmm7		/* a[jp+p12]*s12 */				\n\t		mulpd		0x90(%%rbx)	,%%xmm15		/* a[jp+p12]*s12 */\n\t"\
		"movq		%%rax,%%rdx	/* r17 */									\n\t		/* r25 */\n\t"\
		"addpd		%%xmm6		,%%xmm5	/* xmm5 <- it */					\n\t		addpd		%%xmm14		,%%xmm13	/* xmm13 <- it */\n\t"\
		"subpd		%%xmm7		,%%xmm4	/* xmm4 <- rt */					\n\t		subpd		%%xmm15		,%%xmm12	/* xmm12 <- rt */\n\t"\
		"movaps		%%xmm5		,0x10(%%rdx)	/* store it */				\n\t		movaps		%%xmm13		,0x90(%%rdx)	/* store it */\n\t"\
		"movaps		%%xmm4		,    (%%rdx)	/* store rt */				\n\t		movaps		%%xmm12		,0x80(%%rdx)	/* store rt */\n\t"\
		"\n\t																\n\t"\
		"addq	$0x40 ,%%rax	/* r21 */									\n\t		/* r29 */\n\t"\
		"subq	$0x20 ,%%rbx	/* c5 */									\n\t		/* c7  */\n\t"\
		"movaps		    (%%rax)	,%%xmm4		/* a[jt+p4] */					\n\t		movaps		0x80(%%rax)	,%%xmm12		/* a[jt+p4] */\n\t"\
		"movaps		0x10(%%rax)	,%%xmm5		/* a[jp+p4] */					\n\t		movaps		0x90(%%rax)	,%%xmm13		/* a[jp+p4] */\n\t"\
		"movaps			%%xmm4	,%%xmm6		/* xmm4 <- cpy a[jt+p4] */		\n\t		movaps			%%xmm12	,%%xmm14		/* xmm12 <- cpy a[jt+p4] */\n\t"\
		"movaps			%%xmm5	,%%xmm7		/* xmm5 <- cpy a[jp+p4] */		\n\t		movaps			%%xmm13	,%%xmm15		/* xmm13 <- cpy a[jp+p4] */\n\t"\
		"mulpd		    (%%rbx)	,%%xmm4		/* a[jt+p4]*c4 */				\n\t		mulpd		0x80(%%rbx)	,%%xmm12		/* a[jt+p4]*c4 */\n\t"\
		"mulpd		    (%%rbx)	,%%xmm5		/* a[jp+p4]*c4 */				\n\t		mulpd		0x80(%%rbx)	,%%xmm13		/* a[jp+p4]*c4 */\n\t"\
		"mulpd		0x10(%%rbx)	,%%xmm6		/* a[jt+p4]*s4 */				\n\t		mulpd		0x90(%%rbx)	,%%xmm14		/* a[jt+p4]*s4 */\n\t"\
		"mulpd		0x10(%%rbx)	,%%xmm7		/* a[jp+p4]*s4 */				\n\t		mulpd		0x90(%%rbx)	,%%xmm15		/* a[jp+p4]*s4 */\n\t"\
		"addpd		%%xmm6		,%%xmm5	/* xmm5 <- t6 */					\n\t		addpd		%%xmm14		,%%xmm13	/* xmm13 <- t6 */\n\t"\
		"subpd		%%xmm7		,%%xmm4	/* xmm4 <- t5 	xmm6,7 free */		\n\t		subpd		%%xmm15		,%%xmm12	/* xmm12 <- t5 	xmm14,7 free */\n\t"\
		"movaps		%%xmm5		,%%xmm7	/* xmm7 <- cpy t6 */				\n\t		movaps		%%xmm13		,%%xmm15	/* xmm15 <- cpy t6 */\n\t"\
		"movaps		%%xmm4		,%%xmm6	/* xmm6 <- cpy t5 */				\n\t		movaps		%%xmm12		,%%xmm14	/* xmm14 <- cpy t5 */\n\t"\
	"prefetcht1	0x80(%%r15)\n\t"\
		"subpd		    (%%rdx),%%xmm4		/* ~t7 <- t5 -rt */				\n\t		subpd		0x80(%%rdx)	,%%xmm12		/* ~t7 <- t5 -rt */\n\t"\
		"subpd		0x10(%%rdx),%%xmm5		/* ~t8 <- t6 -it */				\n\t		subpd		0x90(%%rdx)	,%%xmm13		/* ~t8 <- t6 -it */\n\t"\
		"addpd		    (%%rdx),%%xmm6		/* ~t5 <- t5 +rt */				\n\t		addpd		0x80(%%rdx)	,%%xmm14		/* ~t5 <- t5 +rt */\n\t"\
		"addpd		0x10(%%rdx),%%xmm7		/* ~t6 <- t6 +it */				\n\t		addpd		0x90(%%rdx)	,%%xmm15		/* ~t6 <- t6 +it */\n\t"\
		"\n\t																\n\t"\
	"/* Finish radix-4 butterfly and store results into temp-array slots: */\n\t	/* Finish radix-4 butterfly and store results into temporary-array slots: */\n\t"\
		"subpd		%%xmm6,%%xmm0	/*~t5 */								\n\t		subpd		%%xmm14,%%xmm8 	/*~t5 */						\n\t"\
		"subpd		%%xmm7,%%xmm1	/*~t6 */								\n\t		subpd		%%xmm15,%%xmm9 	/*~t6 */						\n\t"\
		"subpd		%%xmm5,%%xmm2	/*~t3 */								\n\t		subpd		%%xmm13,%%xmm10	/*~t3 */\n\t"\
		"subpd		%%xmm4,%%xmm3	/*~t8 */								\n\t		subpd		%%xmm12,%%xmm11	/*~t8 */\n\t"\
		"movaps		%%xmm0,0x40(%%rdx)	/* a[jt+p8 ] <- ~t5 */				\n\t		movaps		%%xmm8 ,0xc0(%%rdx)	/* a[jt+p8 ] <- ~t5 */	\n\t"\
		"movaps		%%xmm1,0x50(%%rdx)	/* a[jp+p8 ] <- ~t6 */				\n\t		movaps		%%xmm9 ,0xd0(%%rdx)	/* a[jp+p8 ] <- ~t6 */	\n\t"\
		"movaps		%%xmm2,0x20(%%rdx)	/* a[jt+p4 ] <- ~t3 */				\n\t		movaps		%%xmm10,0xa0(%%rdx)	/* a[jt+p4 ] <- ~t3 */\n\t"\
		"movaps		%%xmm3,0x70(%%rdx)	/* a[jp+p12] <- ~t8 */				\n\t		movaps		%%xmm11,0xf0(%%rdx)	/* a[jp+p12] <- ~t8 */\n\t"\
		"addpd		%%xmm6,%%xmm6	/* 2*t5 */								\n\t		addpd		%%xmm14,%%xmm14	/* 2*t5 */						\n\t"\
		"addpd		%%xmm7,%%xmm7	/* 2*t6 */								\n\t		addpd		%%xmm15,%%xmm15	/* 2*t6 */						\n\t"\
		"addpd		%%xmm5,%%xmm5	/* 2*t8 */								\n\t		addpd		%%xmm13,%%xmm13	/* 2*t8 */\n\t"\
		"addpd		%%xmm4,%%xmm4	/* 2*t7 */								\n\t		addpd		%%xmm12,%%xmm12	/* 2*t7 */\n\t"\
		"addpd		%%xmm0,%%xmm6	/*~t1 */								\n\t		addpd		%%xmm8 ,%%xmm14	/*~t1 */						\n\t"\
		"addpd		%%xmm1,%%xmm7	/*~t2 */								\n\t		addpd		%%xmm9 ,%%xmm15	/*~t2 */						\n\t"\
		"addpd		%%xmm2,%%xmm5	/*~t7 */								\n\t		addpd		%%xmm10,%%xmm13	/*~t7 */\n\t"\
		"addpd		%%xmm3,%%xmm4	/*~t4 */								\n\t		addpd		%%xmm11,%%xmm12	/*~t4 */\n\t"\
		"movaps		%%xmm6,    (%%rdx)	/* a[jt    ] <- ~t1 */				\n\t		movaps		%%xmm14,0x80(%%rdx)	/* a[jt    ] <- ~t1 */	\n\t"\
		"movaps		%%xmm7,0x10(%%rdx)	/* a[jp    ] <- ~t2 */				\n\t		movaps		%%xmm15,0x90(%%rdx)	/* a[jp    ] <- ~t2 */	\n\t"\
		"movaps		%%xmm5,0x60(%%rdx)	/* a[jt+p12] <- ~t7 */				\n\t		movaps		%%xmm13,0xe0(%%rdx)	/* a[jt+p12] <- ~t7 */\n\t"\
		"movaps		%%xmm4,0x30(%%rdx)	/* a[jp+p4 ] <- ~t4 */				\n\t		movaps		%%xmm12,0xb0(%%rdx)	/* a[jp+p4 ] <- ~t4 */\n\t"\
		"\n\t"\
"/**************************************************************************************/\n\t"\
"/*...and now do four more radix-4 transforms, including the internal twiddle factors: */\n\t"\
"/**************************************************************************************/\n\t"\
		"\n\t"\
"/*...Block 1: t1,9,17,25 */\n\t												/*...Block 3: t5,13,21,29: All rax-offsets incr +0x40 in rcol w.r.to lcol: */\n\t"\
		"movq	%[__r1] ,%%rax\n\t														movaps		0x040(%%rax),%%xmm8 		/* t5  */\n\t"\
		"movq	%[__isrt2],%%rsi\n\t													movaps		0x050(%%rax),%%xmm9 		/* t6  */\n\t"\
		"\n\t																			movaps		0x0d0(%%rax),%%xmm11		/* t14 */\n\t"\
		"movaps		     (%%rax),%%xmm0		/* t1  */\n\t								movaps		0x0c0(%%rax),%%xmm10		/* t13 */\n\t"\
		"movaps		0x010(%%rax),%%xmm1		/* t2  */\n\t								subpd		%%xmm11,%%xmm8 		/* t5 =t5 -t14 */\n\t"\
		"movaps		0x080(%%rax),%%xmm2		/* t9  */\n\t								subpd		%%xmm10,%%xmm9 		/* t14=t6 -t13 */\n\t"\
		"movaps		0x090(%%rax),%%xmm3		/* t14 */\n\t								addpd		%%xmm11,%%xmm11		/* 2.t14 */\n\t"\
	"prefetcht1	0xc0(%%r14)\n\t"\
		"\n\t																			addpd		%%xmm10,%%xmm10		/* 2.t13 */\n\t"\
		"subpd		%%xmm2,%%xmm0		/* t9 =t1 -t9  */\n\t							addpd		%%xmm8 ,%%xmm11		/* t13=t5 +t14 */\n\t"\
		"subpd		%%xmm3,%%xmm1		/* t14=t2 -t14 */\n\t							addpd		%%xmm9 ,%%xmm10		/* t6 =t6 +t13 */\n\t"\
		"addpd		%%xmm2,%%xmm2		/* 2.t9  */\n\t									movaps		0x140(%%rax)	,%%xmm12		/* t21 */\n\t"\
		"addpd		%%xmm3,%%xmm3		/* 2.t14 */\n\t									movaps		0x150(%%rax)	,%%xmm13		/* t22 */\n\t"\
		"addpd		%%xmm0,%%xmm2		/* t1 =t1 +t9  */\n\t							movaps		0x1c0(%%rax)	,%%xmm14		/* t29 */\n\t"\
		"addpd		%%xmm1,%%xmm3		/* t2 =t2 +t14 */\n\t							movaps		0x1d0(%%rax)	,%%xmm15		/* t30 */\n\t"\
		"\n\t																			subpd		0x150(%%rax)	,%%xmm12		/* t21-t22 */\n\t"\
		"movaps		0x100(%%rax),%%xmm4		/* t17 */\n\t								addpd		0x140(%%rax)	,%%xmm13		/* t22+t21 */\n\t"\
		"movaps		0x110(%%rax),%%xmm5		/* t18 */\n\t								mulpd		(%%rsi)		,%%xmm12	/* t21 = (t21-t22)*ISRT2 */\n\t"\
		"movaps		0x180(%%rax),%%xmm6		/* t25 */\n\t								mulpd		(%%rsi)		,%%xmm13	/* t22 = (t22+t21)*ISRT2 */\n\t"\
		"movaps		0x190(%%rax),%%xmm7		/* t26 */\n\t								addpd		0x1d0(%%rax)	,%%xmm14		/* t29+t30 */\n\t"\
		"\n\t																			subpd		0x1c0(%%rax)	,%%xmm15		/* t30-t29 */\n\t"\
		"subpd		%%xmm6,%%xmm4		/* t25=t17-t25 */\n\t							mulpd		(%%rsi)		,%%xmm14	/*  rt = (t29+t30)*ISRT2 */\n\t"\
		"subpd		%%xmm7,%%xmm5		/* t26=t18-t26 */\n\t							mulpd		(%%rsi)		,%%xmm15	/*  it = (t30-t29)*ISRT2 */\n\t"\
		"addpd		%%xmm6,%%xmm6		/* 2.t25 */\n\t									subpd		%%xmm14		,%%xmm12		/* t21=t21-rt */\n\t"\
		"addpd		%%xmm7,%%xmm7		/* 2.t26 */\n\t									subpd		%%xmm15		,%%xmm13		/* t22=t22-it */\n\t"\
		"addpd		%%xmm4,%%xmm6		/* t17=t17+t25 */\n\t							addpd		%%xmm14		,%%xmm14		/*      2* rt */\n\t"\
		"addpd		%%xmm5,%%xmm7		/* t18=t18+t26 */\n\t							addpd		%%xmm15		,%%xmm15		/*      2* it */\n\t"\
		"\n\t																			addpd		%%xmm12		,%%xmm14		/* t29=t21+rt */\n\t"\
		"subpd		%%xmm6		,%%xmm2		/* t1  <- t1 -t17 */\n\t					addpd		%%xmm13		,%%xmm15		/* t30=t22+it */\n\t"\
		"subpd		%%xmm7		,%%xmm3		/* t2  <- t2 -t18 */\n\t					subpd		%%xmm12		,%%xmm8 		/* t5 -t21 */\n\t"\
		"addpd		%%xmm6		,%%xmm6		/*          2*t17 */\n\t					subpd		%%xmm13		,%%xmm10		/* t6 -t22 */\n\t"\
		"addpd		%%xmm7		,%%xmm7		/*          2*t18 */\n\t					addpd		%%xmm12		,%%xmm12		/*   2*t21 */\n\t"\
		"addpd		%%xmm2		,%%xmm6		/* t17 <- t1 +t17 */\n\t					addpd		%%xmm13		,%%xmm13		/*   2*t22 */\n\t"\
		"addpd		%%xmm3		,%%xmm7		/* t18 <- t2 +t18 */\n\t					addpd		%%xmm8 		,%%xmm12		/* t5 +t21 */\n\t"\
		"/* x in 2, y in 3, spill 6 and use as tmp: */\n\t								addpd		%%xmm10		,%%xmm13		/* t6 +t22 */\n\t"\
		"movaps	%%xmm6,(%%rax)	/* tmp-store t17 in t0 */\n\t							/* x in 0, y in 2, spill 4 and use as tmp: */\n\t"\
		"movaps	%%xmm2,%%xmm6	/* cpy x */\n\t											movaps	%%xmm12,0x040(%%rax)	/* tmp-store t17 in t0 */\n\t"\
		"subpd	%%xmm3,%%xmm2	/* x-y */\n\t											movaps	%%xmm8 ,%%xmm12	/* cpy x */\n\t"\
		"addpd	%%xmm3,%%xmm3	/* 2*y */\n\t											subpd	%%xmm10,%%xmm8 	/* x-y */\n\t"\
		"mulpd	%%xmm3,%%xmm6	/* 2xy */\n\t											addpd	%%xmm10,%%xmm10	/* 2*y */\n\t"\
		"addpd	%%xmm2,%%xmm3	/* x+y */\n\t											mulpd	%%xmm10,%%xmm12	/* 2xy */\n\t"\
		"mulpd	%%xmm3,%%xmm2	/* x^2-y^2 */\n\t										addpd	%%xmm8 ,%%xmm10	/* x+y */\n\t"\
		"movaps		%%xmm6,0x010(%%rax)	/* a[jp+p1 ], store in t18 */\n\t				mulpd	%%xmm10,%%xmm8 	/* x^2-y^2 */\n\t"\
		"movaps	         (%%rax),%%xmm6	/* a[jt+p0 ], reload */\n\t						movaps		%%xmm12,0x050(%%rax)	/* a[jp+p1 ], store in t18 */\n\t"\
		"movaps		%%xmm2,     (%%rax)	/* a[jt+p1 ], store in t17 */\n\t				movaps	    0x040(%%rax),%%xmm12	/* a[jt+p0 ], reload */\n\t"\
		"/* Have 2 free regs for remaining 3 squarings: */\n\t							movaps		%%xmm8 ,0x040(%%rax)	/* a[jt+p1 ], store in t17 */\n\t"\
		"movaps	%%xmm6,%%xmm2	/* cpy x */\n\t											/* Have 2 free regs for remaining 3 squarings: */\n\t"\
		"movaps	%%xmm6,%%xmm3	/* cpy x */\n\t											movaps	%%xmm12,%%xmm8 	\n\t"\
		"addpd	%%xmm7,%%xmm6	/* x+y */\n\t											movaps	%%xmm12,%%xmm10	\n\t"\
		"subpd	%%xmm7,%%xmm2	/* x-y */\n\t											addpd	%%xmm13,%%xmm12	\n\t"\
		"addpd	%%xmm7,%%xmm7	/* 2*y */\n\t											subpd	%%xmm13,%%xmm8 	\n\t"\
		"mulpd	%%xmm2,%%xmm6	/* x^2-y^2 */\n\t										addpd	%%xmm13,%%xmm13	\n\t"\
		"mulpd	%%xmm3,%%xmm7	/* 2xy */\n\t											mulpd	%%xmm8 ,%%xmm12	\n\t"\
	"prefetcht1	0xc0(%%r15)\n\t"\
		"\n\t																			mulpd	%%xmm10,%%xmm13	\n\t"\
		"subpd		%%xmm5		,%%xmm0		/* t9  <- t9 -t26 */\n\t					subpd		%%xmm15,%%xmm11	\n\t"\
		"subpd		%%xmm4		,%%xmm1		/* t10 <- t10-t25 */\n\t					subpd		%%xmm14,%%xmm9 	\n\t"\
		"addpd		%%xmm5		,%%xmm5		/*          2*t26 */\n\t					addpd		%%xmm15,%%xmm15	\n\t"\
		"addpd		%%xmm4		,%%xmm4		/*          2*t25 */\n\t					addpd		%%xmm14,%%xmm14	\n\t"\
		"addpd		%%xmm0		,%%xmm5		/* t26 <- t9 +t26 */\n\t					addpd		%%xmm11,%%xmm15	\n\t"\
		"addpd		%%xmm1		,%%xmm4		/* t25 <- t10+t25 */\n\t					addpd		%%xmm9 ,%%xmm14	\n\t"\
		"movaps	%%xmm5,%%xmm2	/* cpy x */\n\t											movaps	%%xmm15,%%xmm8 	\n\t"\
		"movaps	%%xmm5,%%xmm3	/* cpy x */\n\t											movaps	%%xmm15,%%xmm10	\n\t"\
		"addpd	%%xmm1,%%xmm5	/* x+y */\n\t											addpd	%%xmm9 ,%%xmm15	\n\t"\
		"subpd	%%xmm1,%%xmm2	/* x-y */\n\t											subpd	%%xmm9 ,%%xmm8 	\n\t"\
		"addpd	%%xmm1,%%xmm1	/* 2*y */\n\t											addpd	%%xmm9 ,%%xmm9 	\n\t"\
		"mulpd	%%xmm2,%%xmm5	/* x^2-y^2 */\n\t										mulpd	%%xmm8 ,%%xmm15	\n\t"\
		"mulpd	%%xmm3,%%xmm1	/* 2xy */\n\t											mulpd	%%xmm10,%%xmm9 	\n\t"\
		"\n\t																			\n\t"\
		"movaps	%%xmm0,%%xmm2	/* cpy x */\n\t											movaps	%%xmm11,%%xmm8 	\n\t"\
		"movaps	%%xmm0,%%xmm3	/* cpy x */\n\t											movaps	%%xmm11,%%xmm10	\n\t"\
		"addpd	%%xmm4,%%xmm0	/* x+y */\n\t											addpd	%%xmm14,%%xmm11	\n\t"\
		"subpd	%%xmm4,%%xmm2	/* x-y */\n\t											subpd	%%xmm14,%%xmm8 	\n\t"\
		"addpd	%%xmm4,%%xmm4	/* 2*y */\n\t											addpd	%%xmm14,%%xmm14	\n\t"\
		"mulpd	%%xmm2,%%xmm0	/* x^2-y^2 */\n\t										mulpd	%%xmm8 ,%%xmm11	\n\t"\
		"mulpd	%%xmm3,%%xmm4	/* 2xy */\n\t											mulpd	%%xmm10,%%xmm14	\n\t"\
		"\n\t																			\n\t"\
		"movaps	     (%%rax),%%xmm2	/* a[jt+p1 ], reload */\n\t							movaps	0x040(%%rax),%%xmm8 	/* a[jt+p1 ], reload */\n\t"\
		"movaps	0x010(%%rax),%%xmm3	/* a[jp+p1 ], reload */\n\t							movaps	0x050(%%rax),%%xmm10	/* a[jp+p1 ], reload */\n\t"\
	"/* SSE2_RADIX4_DIT_IN_PLACE_C(xmm6,7,2,3,0,4,5,1) - stores all 8 memlocs */\n\t/* SSE2_RADIX4_DIT_IN_PLACE_C(xmm12,5,0,2,3,6,7,1): */\n\t"\
		"subpd	%%xmm2,%%xmm6	/* t3 */		\n\t									subpd	%%xmm8 ,%%xmm12	\n\t"\
		"subpd	%%xmm3,%%xmm7	/* t4 */		\n\t									subpd	%%xmm10,%%xmm13	\n\t"\
		"addpd	%%xmm2,%%xmm2	/* 2*y */		\n\t									addpd	%%xmm8 ,%%xmm8 	\n\t"\
		"addpd	%%xmm3,%%xmm3	/* 2*y */		\n\t									addpd	%%xmm10,%%xmm10	\n\t"\
		"addpd	%%xmm6,%%xmm2	/* t1 */		\n\t									addpd	%%xmm12,%%xmm8 	\n\t"\
		"addpd	%%xmm7,%%xmm3	/* t2 */		\n\t									addpd	%%xmm13,%%xmm10	\n\t"\
		"subpd	%%xmm5,%%xmm0	/* t7 */		\n\t									subpd	%%xmm15,%%xmm11	\n\t"\
		"subpd	%%xmm1,%%xmm4	/* t8 */		\n\t									subpd	%%xmm9 ,%%xmm14	\n\t"\
		"subpd	%%xmm0,%%xmm7	/* ~t4 <- t4 -t7 */	\n\t								subpd	%%xmm11,%%xmm13	\n\t"\
		"addpd	%%xmm5,%%xmm5	/* 2*y */		\n\t									addpd	%%xmm15,%%xmm15	\n\t"\
		"subpd	%%xmm4,%%xmm6	/* ~t7 <- t3 -t8 */	\n\t								subpd	%%xmm14,%%xmm12	\n\t"\
		"addpd	%%xmm1,%%xmm1	/* 2*y */		\n\t									addpd	%%xmm9 ,%%xmm9 	\n\t"\
		"movaps	%%xmm7,0x090(%%rax)	/* <- ~t4 */\n\t									movaps	%%xmm13,0x0d0(%%rax)	\n\t"\
		"addpd	%%xmm0,%%xmm5	/* t4 */		\n\t									addpd	%%xmm11,%%xmm15	\n\t"\
	"prefetcht1	0x100(%%r14)\n\t"\
		"movaps	%%xmm6,0x180(%%rax)	/* <- ~t7 */\n\t									movaps	%%xmm12,0x1c0(%%rax)	\n\t"\
		"addpd	%%xmm4,%%xmm1	/* t5 */		\n\t									addpd	%%xmm14,%%xmm9 	\n\t"\
		"addpd	%%xmm4,%%xmm4	/*          2*t8 */		\n\t							addpd	%%xmm14,%%xmm14	\n\t"\
		"subpd	%%xmm5,%%xmm2	/* ~t5 <- t1 -t5 */		\n\t							subpd	%%xmm15,%%xmm8 	\n\t"\
		"addpd	%%xmm0,%%xmm0	/*          2*t7 */		\n\t							addpd	%%xmm11,%%xmm11	\n\t"\
		"subpd	%%xmm1,%%xmm3	/* ~t6 <- t2 -t6 */		\n\t							subpd	%%xmm9 ,%%xmm10	\n\t"\
		"addpd	%%xmm6,%%xmm4	/* ~t3 <- t3 +t8 */		\n\t							addpd	%%xmm12,%%xmm14	\n\t"\
		"movaps	%%xmm2,0x100(%%rax)	/* <- ~t5 */	\n\t								movaps	%%xmm8 ,0x140(%%rax)	\n\t"\
		"addpd	%%xmm7,%%xmm0	/* ~t8 <- t4 +t7 */		\n\t							addpd	%%xmm13,%%xmm11	\n\t"\
		"movaps	%%xmm3,0x110(%%rax)	/* <- ~t6 */	\n\t								movaps	%%xmm10,0x150(%%rax)	\n\t"\
		"movaps	%%xmm4,0x080(%%rax)	/* <- ~t3 */	\n\t								movaps	%%xmm14,0x0c0(%%rax)	\n\t"\
		"addpd	%%xmm5,%%xmm5	/*          2*t5 */		\n\t							addpd	%%xmm15,%%xmm15	\n\t"\
		"movaps	%%xmm0,0x190(%%rax)	/* <- ~t8 */	\n\t								movaps	%%xmm11,0x1d0(%%rax)	\n\t"\
		"addpd	%%xmm1,%%xmm1	/*          2*t6 */		\n\t							addpd	%%xmm9 ,%%xmm9 	\n\t"\
		"addpd	%%xmm2,%%xmm5	/* ~t1 <- t1 +t5 */		\n\t							addpd	%%xmm8 ,%%xmm15	\n\t"\
		"addpd	%%xmm3,%%xmm1	/* ~t2 <- t2 +t6 */		\n\t							addpd	%%xmm10,%%xmm9 	\n\t"\
		"movaps	%%xmm5,     (%%rax)	/* <- ~t1 */	\n\t								movaps	%%xmm15,0x040(%%rax)	\n\t"\
		"movaps	%%xmm1,0x010(%%rax)	/* <- ~t2 */	\n\t								movaps	%%xmm9 ,0x050(%%rax)	\n\t"\
		"\n\t"\
"/*...Block 2: t3,11,19,27 */\n\t"\
		"addq	$0x20,%%rax	/* r3  */	\n\t"\
		"leaq	0x10(%%rsi),%%rdi	/* cc0, from isrt2 */\n\t"\
		"\n\t																	/*...Block 4: t7,15,23,31 */\n\t"\
		"movaps		0x100(%%rax),%%xmm4		/* t19 */			\n\t					movaps		0x140(%%rax),%%xmm12		/* t23 */			\n\t"\
		"movaps		0x110(%%rax),%%xmm5		/* t20 */			\n\t					movaps		0x150(%%rax),%%xmm13		/* t24 */			\n\t"\
		"movaps		0x180(%%rax),%%xmm6	/* t27 */\n\t									movaps		0x1c0(%%rax),%%xmm14		/* t31 */\n\t"\
		"movaps		0x190(%%rax),%%xmm7	/* t28 */\n\t									movaps		0x1d0(%%rax),%%xmm15		/* t32 */\n\t"\
		"movaps		0x100(%%rax),%%xmm0		/* copy t19 */		\n\t					movaps		0x140(%%rax),%%xmm8 		/* copy t23 */		\n\t"\
		"movaps		0x110(%%rax),%%xmm1		/* copy t20 */		\n\t					movaps		0x150(%%rax),%%xmm9 		/* copy t24 */		\n\t"\
		"movaps		0x180(%%rax),%%xmm2	/* copy t27 */\n\t								movaps		0x1c0(%%rax),%%xmm10		/* copy t31 */\n\t"\
		"movaps		0x190(%%rax),%%xmm3	/* copy t28 */\n\t								movaps		0x1d0(%%rax),%%xmm11		/* copy t32 */\n\t"\
		"\n\t																			\n\t"\
		"mulpd		    (%%rdi)	,%%xmm4		/* t19*c */			\n\t					mulpd		0x10(%%rdi)	,%%xmm12		/* t23*s */			\n\t"\
		"mulpd		0x10(%%rdi)	,%%xmm1		/* t20*s */			\n\t					mulpd		    (%%rdi)	,%%xmm9 		/* t24*c */			\n\t"\
		"mulpd		0x10(%%rdi)	,%%xmm6		/* t27*s */\n\t								mulpd		    (%%rdi)	,%%xmm14		/* t31*c */\n\t"\
		"mulpd			(%%rdi)	,%%xmm3		/* t28*c */\n\t								mulpd		0x10(%%rdi)	,%%xmm11		/* t32*s */\n\t"\
		"mulpd		    (%%rdi)	,%%xmm5		/* t20*c */			\n\t					mulpd		0x10(%%rdi)	,%%xmm13		/* t24*s */			\n\t"\
		"mulpd		0x10(%%rdi)	,%%xmm7		/* t28*s */\n\t								mulpd		    (%%rdi)	,%%xmm8 		/* t23*c */			\n\t"\
		"mulpd		0x10(%%rdi)	,%%xmm0		/* t19*s */			\n\t					mulpd		    (%%rdi)	,%%xmm15		/* t32*c */\n\t"\
		"mulpd			(%%rdi)	,%%xmm2		/* t27*c */\n\t								mulpd		0x10(%%rdi)	,%%xmm10		/* t31*s */\n\t"\
		"subpd		%%xmm1,%%xmm4	/* ~t19 */				\n\t						subpd		%%xmm9 		,%%xmm12	/* ~t23 */				\n\t"\
		"subpd		%%xmm3,%%xmm6	/* rt */\n\t										addpd		%%xmm8 		,%%xmm13	/* ~t24 */				\n\t"\
		"addpd		%%xmm0,%%xmm5	/* ~t20 */				\n\t						subpd		%%xmm11		,%%xmm14		/* rt */\n\t"\
		"addpd		%%xmm2,%%xmm7	/* it */\n\t										addpd		%%xmm10		,%%xmm15		/* it */\n\t"\
	"prefetcht1	0x100(%%r15)\n\t"\
		"\n\t																			\n\t"\
		"subpd		%%xmm6,%%xmm4		/*~t27=t19-rt */\n\t							subpd		%%xmm14		,%%xmm12		/*~t23=t23-rt */\n\t"\
		"subpd		%%xmm7,%%xmm5		/*~t28=t20-it */\n\t							subpd		%%xmm15		,%%xmm13		/*~t24=t24-it */\n\t"\
		"addpd		%%xmm6,%%xmm6		/*      2* rt */\n\t							addpd		%%xmm14		,%%xmm14		/*      2* rt */\n\t"\
		"addpd		%%xmm7,%%xmm7		/*      2* it */\n\t							addpd		%%xmm15		,%%xmm15		/*      2* it */\n\t"\
		"addpd		%%xmm4,%%xmm6		/*~t19=t19+rt */\n\t							addpd		%%xmm12		,%%xmm14		/*~t31=t23+rt */\n\t"\
		"addpd		%%xmm5,%%xmm7		/*~t20=t20+it */\n\t							addpd		%%xmm13		,%%xmm15		/*~t32=t24+it */\n\t"\
		"\n\t																			\n\t"\
		"movaps		0x080(%%rax),%%xmm2		/* t11 */\n\t								movaps		0x0c0(%%rax)	,%%xmm10		/* t15 */\n\t"\
		"movaps		0x090(%%rax),%%xmm3		/* t12 */\n\t								movaps		0x0d0(%%rax)	,%%xmm11		/* t16 */\n\t"\
		"subpd		0x090(%%rax),%%xmm2		/* t11-t12 */\n\t							addpd		0x0d0(%%rax)	,%%xmm10		/* t15+t16 */\n\t"\
		"addpd		0x080(%%rax),%%xmm3		/* t12+t11 */\n\t							subpd		0x0c0(%%rax)	,%%xmm11		/* t16-t15 */\n\t"\
		"mulpd		(%%rsi),%%xmm2	/* rt = (t11-t12)*ISRT2 */\n\t						mulpd		(%%rsi)		,%%xmm10	/* rt = (t15+t16)*ISRT2 */\n\t"\
		"mulpd		(%%rsi),%%xmm3	/* it = (t12+t11)*ISRT2 */\n\t						mulpd		(%%rsi)		,%%xmm11	/* it = (t16-t15)*ISRT2 */\n\t"\
		"\n\t																			\n\t"\
		"movaps		     (%%rax),%%xmm0		/* t3  */\n\t								movaps		0x040(%%rax)	,%%xmm8 		/* t7  */\n\t"\
		"movaps		0x010(%%rax),%%xmm1		/* t4  */\n\t								movaps		0x050(%%rax)	,%%xmm9 		/* t8  */\n\t"\
		"\n\t																			\n\t"\
		"subpd		%%xmm2,%%xmm0					/*~t11=t3 -rt */\n\t				subpd		%%xmm10,%%xmm8 					/*~t7 =t7 -rt */\n\t"\
		"subpd		%%xmm3,%%xmm1					/*~t12=t4 -it */\n\t				subpd		%%xmm11,%%xmm9 					/*~t8 =t8 -it */\n\t"\
		"addpd		     (%%rax),%%xmm2		/*~t3 =rt +t3 */\n\t						addpd		0x040(%%rax)	,%%xmm10		/*~t15=rt +t7 */\n\t"\
		"addpd		0x010(%%rax),%%xmm3		/*~t4 =it +t4 */\n\t						addpd		0x050(%%rax)	,%%xmm11		/*~t16=it +t8 */\n\t"\
		"\n\t																			\n\t"\
		"subpd		%%xmm6,%%xmm2		/* t3 -t19 */\n\t								subpd		%%xmm12		,%%xmm8 		/* t7 -t23 */\n\t"\
		"subpd		%%xmm7,%%xmm3		/* t4 -t20 */\n\t								subpd		%%xmm13		,%%xmm9 		/* t8 -t24 */\n\t"\
		"addpd		%%xmm6,%%xmm6		/*   2*t19 */\n\t								addpd		%%xmm12		,%%xmm12		/*   2*t23 */\n\t"\
		"addpd		%%xmm7,%%xmm7		/*   2*t20 */\n\t								addpd		%%xmm13		,%%xmm13		/*   2*t24 */\n\t"\
		"addpd		%%xmm2,%%xmm6		/* t3 +t19 */\n\t								addpd		%%xmm8 		,%%xmm12		/* t7 +t23 */\n\t"\
		"addpd		%%xmm3,%%xmm7		/* t4 +t20 */\n\t								addpd		%%xmm9 		,%%xmm13		/* t8 +t24 */\n\t"\
		"/* x in 2, y in 3, spill 6 and use as tmp: */\n\t								/* x in 0, y in 1, spill 4 and use as tmp: */\n\t"\
		"movaps	%%xmm6,(%%rax)	/* tmp-store t17 in t0 */\n\t							movaps	%%xmm12,0x040(%%rax)	/* tmp-store t17 in t0 */\n\t"\
		"movaps	%%xmm2,%%xmm6	/* cpy x */\n\t											movaps	%%xmm8 ,%%xmm12	/* cpy x */\n\t"\
		"subpd	%%xmm3,%%xmm2	/* x-y */\n\t											subpd	%%xmm9 ,%%xmm8 	/* x-y */\n\t"\
		"addpd	%%xmm3,%%xmm3	/* 2*y */\n\t											addpd	%%xmm9 ,%%xmm9 	/* 2*y */\n\t"\
		"mulpd	%%xmm3,%%xmm6	/* 2xy */\n\t											mulpd	%%xmm9 ,%%xmm12	/* 2xy */\n\t"\
		"addpd	%%xmm2,%%xmm3	/* x+y */\n\t											addpd	%%xmm8 ,%%xmm9 	/* x+y */\n\t"\
		"mulpd	%%xmm3,%%xmm2	/* x^2-y^2 */\n\t										mulpd	%%xmm9 ,%%xmm8 	/* x^2-y^2 */\n\t"\
		"movaps		%%xmm6,0x010(%%rax)	/* a[jp+p1 ], store in t18 */\n\t				movaps		%%xmm12,0x050(%%rax)	/* a[jp+p1 ], store in t18 */\n\t"\
		"movaps	         (%%rax),%%xmm6	/* a[jt+p0 ], reload */\n\t						movaps	    0x040(%%rax),%%xmm12	/* a[jt+p0 ], reload */\n\t"\
		"movaps		%%xmm2,     (%%rax)	/* a[jt+p1 ], store in t17 */\n\t				movaps		%%xmm8 ,0x040(%%rax)	/* a[jt+p1 ], store in t17 */\n\t"\
	"prefetcht1	0x140(%%r14)\n\t"\
		"/* Have 2 free regs for remaining 3 squarings: */\n\t							/* Have 2 free regs for remaining 3 squarings: */\n\t"\
		"movaps	%%xmm6,%%xmm2	\n\t													movaps	%%xmm12,%%xmm8 	\n\t"\
		"movaps	%%xmm6,%%xmm3	\n\t													movaps	%%xmm12,%%xmm9 	\n\t"\
		"addpd	%%xmm7,%%xmm6	\n\t													addpd	%%xmm13,%%xmm12	\n\t"\
		"subpd	%%xmm7,%%xmm2	\n\t													subpd	%%xmm13,%%xmm8 	\n\t"\
		"addpd	%%xmm7,%%xmm7	\n\t													addpd	%%xmm13,%%xmm13	\n\t"\
		"mulpd	%%xmm2,%%xmm6	\n\t													mulpd	%%xmm8 ,%%xmm12	\n\t"\
		"mulpd	%%xmm3,%%xmm7	\n\t													mulpd	%%xmm9 ,%%xmm13	\n\t"\
		"\n\t																			\n\t"\
		"subpd		%%xmm5,%%xmm0	\n\t												subpd		%%xmm15,%%xmm10	\n\t"\
		"subpd		%%xmm4,%%xmm1	\n\t												subpd		%%xmm14,%%xmm11	\n\t"\
		"addpd		%%xmm5,%%xmm5	\n\t												addpd		%%xmm15,%%xmm15	\n\t"\
		"addpd		%%xmm4,%%xmm4	\n\t												addpd		%%xmm14,%%xmm14	\n\t"\
		"addpd		%%xmm0,%%xmm5	\n\t												addpd		%%xmm10,%%xmm15	\n\t"\
		"addpd		%%xmm1,%%xmm4	\n\t												addpd		%%xmm11,%%xmm14	\n\t"\
		"movaps	%%xmm5,%%xmm2	\n\t													movaps	%%xmm15,%%xmm8 	\n\t"\
		"movaps	%%xmm5,%%xmm3	\n\t													movaps	%%xmm15,%%xmm9 	\n\t"\
		"addpd	%%xmm1,%%xmm5	\n\t													addpd	%%xmm11,%%xmm15	\n\t"\
		"subpd	%%xmm1,%%xmm2	\n\t													subpd	%%xmm11,%%xmm8 	\n\t"\
		"addpd	%%xmm1,%%xmm1	\n\t													addpd	%%xmm11,%%xmm11	\n\t"\
		"mulpd	%%xmm2,%%xmm5	\n\t													mulpd	%%xmm8 ,%%xmm15	\n\t"\
		"mulpd	%%xmm3,%%xmm1	\n\t													mulpd	%%xmm9 ,%%xmm11	\n\t"\
		"\n\t																			\n\t"\
		"movaps	%%xmm0,%%xmm2	\n\t													movaps	%%xmm10,%%xmm8 	\n\t"\
		"movaps	%%xmm0,%%xmm3	\n\t													movaps	%%xmm10,%%xmm9 	\n\t"\
		"addpd	%%xmm4,%%xmm0	\n\t													addpd	%%xmm14,%%xmm10	\n\t"\
		"subpd	%%xmm4,%%xmm2	\n\t													subpd	%%xmm14,%%xmm8 	\n\t"\
		"addpd	%%xmm4,%%xmm4	\n\t													addpd	%%xmm14,%%xmm14	\n\t"\
		"mulpd	%%xmm2,%%xmm0	\n\t													mulpd	%%xmm8 ,%%xmm10	\n\t"\
		"mulpd	%%xmm3,%%xmm4	\n\t													mulpd	%%xmm9 ,%%xmm14	\n\t"\
		"\n\t																			\n\t"\
		"movaps	      (%%rax),%%xmm2	/* a[jt+p1 ], reload */\n\t						movaps	 0x040(%%rax),%%xmm8 	/* a[jt+p1 ], reload */\n\t"\
		"movaps	 0x010(%%rax),%%xmm3	/* a[jp+p1 ], reload */\n\t						movaps	 0x050(%%rax),%%xmm9 	/* a[jp+p1 ], reload */\n\t"\
	"/* SSE2_RADIX4_DIT_IN_PLACE_C(xmm6,xmm7,xmm2,xmm3,xmm0,xmm4,xmm5,xmm1) - Th	/* SSE2_RADIX4_DIT_IN_PLACE_C(xmm12,xmm13,xmm8 ,xmm9 ,xmm10,xmm14,xmm15,xmm11) - This stores all 8 memlocs */\n\t"\
		"subpd	%%xmm2,%%xmm6	\n\t													subpd	%%xmm8 ,%%xmm12	\n\t"\
		"subpd	%%xmm3,%%xmm7	\n\t													subpd	%%xmm9 ,%%xmm13	\n\t"\
		"addpd	%%xmm2,%%xmm2	\n\t													addpd	%%xmm8 ,%%xmm8 	\n\t"\
		"addpd	%%xmm3,%%xmm3	\n\t													addpd	%%xmm9 ,%%xmm9 	\n\t"\
		"addpd	%%xmm6,%%xmm2	\n\t													addpd	%%xmm12,%%xmm8 	\n\t"\
		"addpd	%%xmm7,%%xmm3	\n\t													addpd	%%xmm13,%%xmm9 	\n\t"\
		"subpd	%%xmm5,%%xmm0	\n\t													subpd	%%xmm15,%%xmm10	\n\t"\
		"subpd	%%xmm1,%%xmm4	\n\t													subpd	%%xmm11,%%xmm14	\n\t"\
		"subpd	%%xmm0,%%xmm7	\n\t													subpd	%%xmm10,%%xmm13	\n\t"\
		"addpd	%%xmm5,%%xmm5	\n\t													addpd	%%xmm15,%%xmm15	\n\t"\
		"subpd	%%xmm4,%%xmm6	\n\t													subpd	%%xmm14,%%xmm12	\n\t"\
		"addpd	%%xmm1,%%xmm1	\n\t													addpd	%%xmm11,%%xmm11	\n\t"\
	"prefetcht1	0x140(%%r15)\n\t"\
		"movaps	%%xmm7,0x090(%%rax)	\n\t												movaps	%%xmm13,0x0d0(%%rax)	\n\t"\
		"addpd	%%xmm0,%%xmm5	\n\t													addpd	%%xmm10,%%xmm15	\n\t"\
		"movaps	%%xmm6,0x180(%%rax)	\n\t												movaps	%%xmm12,0x1c0(%%rax)	\n\t"\
		"addpd	%%xmm4,%%xmm1	\n\t													addpd	%%xmm14,%%xmm11	\n\t"\
		"addpd	%%xmm4,%%xmm4	\n\t													addpd	%%xmm14,%%xmm14	\n\t"\
		"subpd	%%xmm5,%%xmm2	\n\t													subpd	%%xmm15,%%xmm8 	\n\t"\
		"addpd	%%xmm0,%%xmm0	\n\t													addpd	%%xmm10,%%xmm10	\n\t"\
		"subpd	%%xmm1,%%xmm3	\n\t													subpd	%%xmm11,%%xmm9 	\n\t"\
		"addpd	%%xmm6,%%xmm4	\n\t													addpd	%%xmm12,%%xmm14	\n\t"\
		"movaps	%%xmm2,0x100(%%rax)	\n\t												movaps	%%xmm8 ,0x140(%%rax)	\n\t"\
		"addpd	%%xmm7,%%xmm0	\n\t													addpd	%%xmm13,%%xmm10	\n\t"\
		"movaps	%%xmm3,0x110(%%rax)	\n\t												movaps	%%xmm9 ,0x150(%%rax)	\n\t"\
		"movaps	%%xmm4,0x080(%%rax)	\n\t												movaps	%%xmm14,0x0c0(%%rax)	\n\t"\
		"addpd	%%xmm5,%%xmm5	\n\t													addpd	%%xmm15,%%xmm15	\n\t"\
		"movaps	%%xmm0,0x190(%%rax)	\n\t												movaps	%%xmm10,0x1d0(%%rax)	\n\t"\
		"addpd	%%xmm1,%%xmm1	\n\t													addpd	%%xmm11,%%xmm11	\n\t"\
		"addpd	%%xmm2,%%xmm5	\n\t													addpd	%%xmm8 ,%%xmm15	\n\t"\
		"addpd	%%xmm3,%%xmm1	\n\t													addpd	%%xmm9 ,%%xmm11	\n\t"\
		"movaps	%%xmm5,     (%%rax)	\n\t												movaps	%%xmm15,0x040(%%rax)	\n\t"\
		"movaps	%%xmm1,0x010(%%rax)	\n\t												movaps	%%xmm11,0x050(%%rax)	\n\t"\
		"\n\t"\
	"/***************************************************************************************************/\n\t"\
	"/*..and now do four more radix-4 transforms, including the internal and external twiddle factors.  */\n\t"\
	"/*  Write even-index 16-byte output pairs to a(j1), odd-index to a(j2), unpack same as on inputs.  */\n\t"\
	"/*  We do the last 2 radix-4 blocks first, to make the unpack-interleaving of outputs smoother.    */\n\t"\
	"/***************************************************************************************************/\n\t"\
	"/* Main-array addresses still in add0,1, no need to re-init: */\n\t"\
		"/*...Block 3: t3,11,19,27 -> r9,13,11,15: */	\n\t		/*...Block 4: t7,15,23,31 -> r25,29,27,31: */	\n\t"\
		"addq	$0x60,%%rax		/* r9 */	\n\t							/* r25 */					\n\t"\
		"movq	%[__isrt2],%%rbx			\n\t					/* All rax-offsets incr +0x100 in rcol w.r.to lcol: */\n\t"\
		"movq	%%rdi,%%rcx		/* cc0 */	\n\t					\n\t"\
		"movaps		0x020(%%rax),%%xmm4			\n\t				movaps		0x120(%%rax),%%xmm12			\n\t"\
		"movaps		0x060(%%rax),%%xmm0			\n\t				movaps		0x160(%%rax),%%xmm8 			\n\t"\
		"movaps		0x030(%%rax),%%xmm5			\n\t				movaps		0x130(%%rax),%%xmm13			\n\t"\
		"movaps		0x070(%%rax),%%xmm1			\n\t				movaps		0x170(%%rax),%%xmm9 			\n\t"\
		"movaps		0x020(%%rax),%%xmm6			\n\t				movaps		0x120(%%rax),%%xmm14			\n\t"\
		"movaps		0x060(%%rax),%%xmm2			\n\t				movaps		0x160(%%rax),%%xmm10			\n\t"\
		"movaps		0x030(%%rax),%%xmm7			\n\t				movaps		0x130(%%rax),%%xmm15			\n\t"\
		"movaps		0x070(%%rax),%%xmm3			\n\t				movaps		0x170(%%rax),%%xmm11			\n\t"\
		"mulpd		    (%%rcx),%%xmm4			\n\t				mulpd		0x10(%%rcx),%%xmm12			\n\t"\
		"mulpd		0x10(%%rcx),%%xmm0			\n\t				mulpd		    (%%rcx),%%xmm8 			\n\t"\
		"mulpd		    (%%rcx),%%xmm5			\n\t				mulpd		0x10(%%rcx),%%xmm13			\n\t"\
		"mulpd		0x10(%%rcx),%%xmm1			\n\t				mulpd		    (%%rcx),%%xmm9 			\n\t"\
		"mulpd		0x10(%%rcx),%%xmm6			\n\t				mulpd		    (%%rcx),%%xmm14			\n\t"\
		"mulpd		    (%%rcx),%%xmm2			\n\t				mulpd		0x10(%%rcx),%%xmm10			\n\t"\
		"mulpd		0x10(%%rcx),%%xmm7			\n\t				mulpd		    (%%rcx),%%xmm15			\n\t"\
		"mulpd		    (%%rcx),%%xmm3			\n\t				mulpd		0x10(%%rcx),%%xmm11			\n\t"\
	"prefetcht1	0x180(%%r14)\n\t"\
		"subpd		%%xmm6,%%xmm5				\n\t				subpd		%%xmm14,%%xmm13				\n\t"\
		"subpd		%%xmm2,%%xmm1				\n\t				subpd		%%xmm10,%%xmm9 				\n\t"\
		"addpd		%%xmm7,%%xmm4				\n\t				addpd		%%xmm15,%%xmm12				\n\t"\
		"addpd		%%xmm3,%%xmm0				\n\t				addpd		%%xmm11,%%xmm8 				\n\t"\
		"movaps		%%xmm5,%%xmm7				\n\t				movaps		%%xmm13,%%xmm15				\n\t"\
		"movaps		%%xmm4,%%xmm6				\n\t				movaps		%%xmm12,%%xmm14				\n\t"\
		"addpd		%%xmm0,%%xmm4				\n\t				addpd		%%xmm8 ,%%xmm12				\n\t"\
		"addpd		%%xmm1,%%xmm5				\n\t				addpd		%%xmm9 ,%%xmm13				\n\t"\
		"subpd		%%xmm0,%%xmm6				\n\t				subpd		%%xmm8 ,%%xmm14				\n\t"\
		"subpd		%%xmm1,%%xmm7				\n\t				subpd		%%xmm9 ,%%xmm15				\n\t"\
		"movaps		0x040(%%rax),%%xmm2			\n\t				movaps		0x140(%%rax),%%xmm10			\n\t"\
		"movaps		0x050(%%rax),%%xmm3			\n\t				movaps		0x150(%%rax),%%xmm11			\n\t"\
		"movaps		     (%%rax),%%xmm0			\n\t				movaps		0x100(%%rax),%%xmm8 			\n\t"\
		"movaps		0x010(%%rax),%%xmm1			\n\t				movaps		0x110(%%rax),%%xmm9 			\n\t"\
		"addpd		0x050(%%rax),%%xmm2			\n\t				subpd		0x150(%%rax),%%xmm10			\n\t"\
		"subpd		0x040(%%rax),%%xmm3			\n\t				addpd		0x140(%%rax),%%xmm11			\n\t"\
		"mulpd		(%%rbx),%%xmm2				\n\t				mulpd		(%%rbx),%%xmm10				\n\t"\
		"mulpd		(%%rbx),%%xmm3				\n\t				mulpd		(%%rbx),%%xmm11				\n\t"\
		"subpd		%%xmm2,%%xmm0				\n\t				subpd		%%xmm10,%%xmm8 				\n\t"\
		"subpd		%%xmm3,%%xmm1				\n\t				subpd		%%xmm11,%%xmm9 				\n\t"\
		"addpd		%%xmm2,%%xmm2				\n\t				addpd		%%xmm10,%%xmm10				\n\t"\
		"addpd		%%xmm3,%%xmm3				\n\t				addpd		%%xmm11,%%xmm11				\n\t"\
		"addpd		%%xmm0,%%xmm2				\n\t				addpd		%%xmm8 ,%%xmm10				\n\t"\
		"addpd		%%xmm1,%%xmm3				\n\t				addpd		%%xmm9 ,%%xmm11				\n\t"\
		"addq	$0x120,%%rcx		/* c1 from cc0  ; c3  in rcol */\n\t	subpd		%%xmm14,%%xmm8 				\n\t"\
		"leaq	0x150(%%rbx),%%rdx	/* c9 from isrt2; c11 in rcol */\n\t	subpd		%%xmm15,%%xmm9 				\n\t"\
		"movq	%[__add1],%%rbx	/* rbx shared between rcol/lcol; rcx/rdx-offsets incr +0x80 in rcol for rest of block: */\n\t"\
		"subpd		%%xmm4,%%xmm2				\n\t				addpd		%%xmm14,%%xmm14				\n\t"\
		"subpd		%%xmm5,%%xmm3				\n\t				addpd		%%xmm15,%%xmm15				\n\t"\
		"addpd		%%xmm4,%%xmm4				\n\t				addpd		%%xmm8 ,%%xmm14				\n\t"\
		"addpd		%%xmm5,%%xmm5				\n\t				addpd		%%xmm9 ,%%xmm15				\n\t"\
		"addpd		%%xmm2,%%xmm4				\n\t				movaps		%%xmm8 ,0x100(%%rax)			\n\t"\
		"addpd		%%xmm3,%%xmm5				\n\t				movaps		%%xmm9 ,0x110(%%rax)			\n\t"\
		"movaps		%%xmm2,     (%%rax)			\n\t				movaps		%%xmm14,%%xmm8 				\n\t"\
		"movaps		%%xmm3,0x010(%%rax)			\n\t				movaps		%%xmm15,%%xmm9 				\n\t"\
		"movaps		%%xmm4,%%xmm2				\n\t				mulpd		0x80(%%rcx),%%xmm14			\n\t"\
		"movaps		%%xmm5,%%xmm3				\n\t				mulpd		0x80(%%rcx),%%xmm15			\n\t"\
		"mulpd		    (%%rcx),%%xmm4			\n\t				mulpd		0x90(%%rcx),%%xmm8 			\n\t"\
		"mulpd		    (%%rcx),%%xmm5			\n\t				mulpd		0x90(%%rcx),%%xmm9 			\n\t"\
		"mulpd		0x10(%%rcx),%%xmm2			\n\t				subpd		%%xmm8 ,%%xmm15				\n\t"\
		"mulpd		0x10(%%rcx),%%xmm3			\n\t				addpd		%%xmm9 ,%%xmm14				\n\t"\
	"prefetcht1	0x180(%%r15)\n\t"\
		"subpd		%%xmm2,%%xmm5				\n\t				movaps		%%xmm15,0x30(%%rbx)			\n\t"\
		"addpd		%%xmm3,%%xmm4				\n\t				movaps		%%xmm14,0x20(%%rbx)			\n\t"\
		"movaps		%%xmm5,0x10(%%rbx)			\n\t				movaps		0x100(%%rax),%%xmm14			\n\t"\
		"movaps		%%xmm4,    (%%rbx)			\n\t				movaps		0x110(%%rax),%%xmm15			\n\t"\
		"movaps		     (%%rax),%%xmm4			\n\t				movaps		%%xmm14,%%xmm8 				\n\t"\
		"movaps		0x010(%%rax),%%xmm5			\n\t				movaps		%%xmm15,%%xmm9 				\n\t"\
		"movaps		%%xmm4,%%xmm2				\n\t				mulpd		0x80(%%rdx),%%xmm14			\n\t"\
		"movaps		%%xmm5,%%xmm3				\n\t				mulpd		0x80(%%rdx),%%xmm15			\n\t"\
		"mulpd		    (%%rdx),%%xmm4			\n\t				mulpd		0x90(%%rdx),%%xmm8 			\n\t"\
		"mulpd		    (%%rdx),%%xmm5			\n\t				mulpd		0x90(%%rdx),%%xmm9 			\n\t"\
		"mulpd		0x10(%%rdx),%%xmm2			\n\t				subpd		%%xmm8 ,%%xmm15				\n\t"\
		"mulpd		0x10(%%rdx),%%xmm3			\n\t				addpd		%%xmm9 ,%%xmm14				\n\t"\
		"subpd		%%xmm2,%%xmm5				\n\t				movaps		%%xmm15,0xb0(%%rbx)			\n\t"\
		"addpd		%%xmm3,%%xmm4				\n\t				movaps		%%xmm14,0xa0(%%rbx)			\n\t"\
		"movaps		%%xmm5,0x90(%%rbx)			\n\t				subpd		%%xmm13,%%xmm10				\n\t"\
		"movaps		%%xmm4,0x80(%%rbx)			\n\t				subpd		%%xmm12,%%xmm11				\n\t"\
		"addq	$0x40,%%rcx	/* c5  from c1; c7  in rcol */\n\t		addpd		%%xmm13,%%xmm13				\n\t"\
		"addq	$0x40,%%rdx	/* c13 from c9; c15 in rcol */\n\t		addpd		%%xmm12,%%xmm12				\n\t"\
		"subpd		%%xmm7,%%xmm0				\n\t				addpd		%%xmm10,%%xmm13				\n\t"\
		"subpd		%%xmm6,%%xmm1				\n\t				addpd		%%xmm11,%%xmm12				\n\t"\
		"addpd		%%xmm7,%%xmm7				\n\t				movaps		%%xmm13,%%xmm8 				\n\t"\
		"addpd		%%xmm6,%%xmm6				\n\t				movaps		%%xmm11,%%xmm9 				\n\t"\
		"addpd		%%xmm0,%%xmm7				\n\t				mulpd		0x80(%%rcx),%%xmm13			\n\t"\
		"addpd		%%xmm1,%%xmm6				\n\t				mulpd		0x80(%%rcx),%%xmm11			\n\t"\
		"movaps		%%xmm7,%%xmm4				\n\t				mulpd		0x90(%%rcx),%%xmm8 			\n\t"\
		"movaps		%%xmm1,%%xmm5				\n\t				mulpd		0x90(%%rcx),%%xmm9 			\n\t"\
		"mulpd		    (%%rcx),%%xmm7			\n\t				subpd		%%xmm8 ,%%xmm11				\n\t"\
		"mulpd		    (%%rcx),%%xmm1			\n\t				addpd		%%xmm9 ,%%xmm13				\n\t"\
		"mulpd		0x10(%%rcx),%%xmm4			\n\t				movaps		%%xmm11,0x70(%%rbx)			\n\t"\
		"mulpd		0x10(%%rcx),%%xmm5			\n\t				movaps		%%xmm13,0x60(%%rbx)			\n\t"\
		"subpd		%%xmm4,%%xmm1				\n\t				movaps		%%xmm10,%%xmm8 				\n\t"\
		"addpd		%%xmm5,%%xmm7				\n\t				movaps		%%xmm12,%%xmm9 				\n\t"\
		"movaps		%%xmm1,0x50(%%rbx)			\n\t				mulpd		0x80(%%rdx),%%xmm10			\n\t"\
		"movaps		%%xmm7,0x40(%%rbx)			\n\t				mulpd		0x80(%%rdx),%%xmm12			\n\t"\
		"movaps		%%xmm0,%%xmm4				\n\t				mulpd		0x90(%%rdx),%%xmm8 			\n\t"\
		"movaps		%%xmm6,%%xmm5				\n\t				mulpd		0x90(%%rdx),%%xmm9 			\n\t"\
		"mulpd		    (%%rdx),%%xmm0			\n\t				subpd		%%xmm8 ,%%xmm12				\n\t"\
		"mulpd		    (%%rdx),%%xmm6			\n\t				addpd		%%xmm9 ,%%xmm10				\n\t"\
		"mulpd		0x10(%%rdx),%%xmm4			\n\t				movaps		%%xmm12,0xf0(%%rbx)			\n\t"\
		"mulpd		0x10(%%rdx),%%xmm5			\n\t				movaps		%%xmm10,0xe0(%%rbx)			\n\t"\
	"prefetcht1	0x1c0(%%r14)\n\t"\
		"subpd		%%xmm4,%%xmm6				\n\t				/*...Block 3: t5,13,21,29 -> r17,21,19,23: */	\n\t"\
		"addpd		%%xmm5,%%xmm0				\n\t				movq	%[__r1],%%rax		/* r17 in rcol */\n\t"\
		"movaps		%%xmm6,0xd0(%%rbx)			\n\t				movq		%[__isrt2],%%rdi			\n\t"\
		"movaps		%%xmm0,0xc0(%%rbx)			\n\t				movaps		(%%rdi),%%xmm10				\n\t"\
		"\n\t														movaps		0x120(%%rax),%%xmm12			\n\t"\
		"/*...Block 1: t1,9,17,25 -> r1,5,3,7: */		\n\t		movaps		0x130(%%rax),%%xmm13			\n\t"\
		"movaps		     (%%rax),%%xmm0			\n\t				movaps		0x160(%%rax),%%xmm8 			\n\t"\
		"movaps		0x010(%%rax),%%xmm1			\n\t				movaps		0x170(%%rax),%%xmm9 			\n\t"\
		"movaps		0x040(%%rax),%%xmm2			\n\t				addpd		0x130(%%rax),%%xmm12			\n\t"\
		"movaps		0x050(%%rax),%%xmm3			\n\t				subpd		0x120(%%rax),%%xmm13			\n\t"\
		"subpd		0x040(%%rax),%%xmm0			\n\t				subpd		0x170(%%rax),%%xmm8 			\n\t"\
		"subpd		0x050(%%rax),%%xmm1			\n\t				addpd		0x160(%%rax),%%xmm9 			\n\t"\
		"addpd		     (%%rax),%%xmm2			\n\t				mulpd		%%xmm10,%%xmm12				\n\t"\
		"addpd		0x010(%%rax),%%xmm3			\n\t				mulpd		%%xmm10,%%xmm13				\n\t"\
		"movaps		0x020(%%rax),%%xmm4			\n\t				mulpd		%%xmm10,%%xmm8 				\n\t"\
		"movaps		0x030(%%rax),%%xmm5			\n\t				mulpd		%%xmm10,%%xmm9 				\n\t"\
		"movaps		0x060(%%rax),%%xmm6			\n\t				movaps		%%xmm12,%%xmm14				\n\t"\
		"movaps		0x070(%%rax),%%xmm7			\n\t				movaps		%%xmm13,%%xmm15				\n\t"\
		"subpd		0x060(%%rax),%%xmm4			\n\t				subpd		%%xmm8 ,%%xmm12				\n\t"\
		"subpd		0x070(%%rax),%%xmm5			\n\t				subpd		%%xmm9 ,%%xmm13				\n\t"\
		"addpd		0x020(%%rax),%%xmm6			\n\t				addpd		%%xmm8 ,%%xmm14				\n\t"\
		"addpd		0x030(%%rax),%%xmm7			\n\t				addpd		%%xmm9 ,%%xmm15				\n\t"\
		"movq	%[__add0],%%rbx				\n\t					movaps		0x100(%%rax),%%xmm8 			\n\t"\
		"subq	$0x140,%%rdx	/* c8 from c13 */\n\t				movaps		0x110(%%rax),%%xmm9 			\n\t"\
		"addpd		%%xmm6,%%xmm2				\n\t				movaps		0x140(%%rax),%%xmm10			\n\t"\
		"addpd		%%xmm7,%%xmm3				\n\t				movaps		0x150(%%rax),%%xmm11			\n\t"\
		"movaps		%%xmm2,    (%%rbx)			\n\t				subpd		0x150(%%rax),%%xmm8 			\n\t"\
		"movaps		%%xmm3,0x10(%%rbx)			\n\t				subpd		0x140(%%rax),%%xmm9 			\n\t"\
		"addpd		%%xmm6,%%xmm6				\n\t				addpd		0x100(%%rax),%%xmm11			\n\t"\
		"addpd		%%xmm7,%%xmm7				\n\t				addpd		0x110(%%rax),%%xmm10			\n\t"\
		"subpd		%%xmm6,%%xmm2				\n\t				subpd		%%xmm12,%%xmm11				\n\t"\
		"subpd		%%xmm7,%%xmm3				\n\t				subpd		%%xmm13,%%xmm9 				\n\t"\
		"movaps		%%xmm2,%%xmm6				\n\t				addpd		%%xmm12,%%xmm12				\n\t"\
		"movaps		%%xmm3,%%xmm7				\n\t				addpd		%%xmm13,%%xmm13				\n\t"\
		"mulpd		    (%%rdx),%%xmm2			\n\t				addpd		%%xmm11,%%xmm12				\n\t"\
		"mulpd		    (%%rdx),%%xmm3			\n\t				addpd		%%xmm9 ,%%xmm13				\n\t"\
		"mulpd		0x10(%%rdx),%%xmm6			\n\t				movaps		%%xmm11,     (%%rax)			\n\t"\
		"mulpd		0x10(%%rdx),%%xmm7			\n\t				movaps		%%xmm9 ,0x010(%%rax)			\n\t"\
		"subpd		%%xmm6,%%xmm3				\n\t				movaps		%%xmm12,%%xmm11				\n\t"\
		"addpd		%%xmm7,%%xmm2				\n\t				movaps		%%xmm13,%%xmm9 				\n\t"\
		"movq	%[__add1],%%rcx				\n\t					mulpd		0x60(%%rdx),%%xmm12	/* c2 */\n\t"\
		"movaps		%%xmm3,%%xmm7				\n\t				mulpd		0x60(%%rdx),%%xmm13			\n\t"\
		"movaps		%%xmm2,%%xmm6				\n\t				mulpd		0x70(%%rdx),%%xmm11			\n\t"\
		"unpckhpd	0x90(%%rcx),%%xmm7			\n\t				mulpd		0x70(%%rdx),%%xmm9 			\n\t"\
		"unpcklpd	0x90(%%rcx),%%xmm3			\n\t				subpd		%%xmm11,%%xmm13				\n\t"\
		"movaps		%%xmm7,0x90(%%rcx)			\n\t				addpd		%%xmm9 ,%%xmm12				\n\t"\
	"prefetcht1	0x1c0(%%r15)\n\t"\
		"addq	$0x20,%%rdx		/* c4 in lcol; c10 in rcol*/\n\t	movaps		%%xmm13,%%xmm11				\n\t"\
		"unpckhpd	0x80(%%rcx),%%xmm6			\n\t				movaps		%%xmm12,%%xmm9 				\n\t"\
		"unpcklpd	0x80(%%rcx),%%xmm2			\n\t				unpckhpd	0x30(%%rcx),%%xmm11			\n\t"\
		"movaps		%%xmm6,0x80(%%rcx)			\n\t				unpcklpd	0x30(%%rcx),%%xmm13			\n\t"\
		"movaps		%%xmm3,0x90(%%rbx)			\n\t				movaps		%%xmm11,0x30(%%rcx)			\n\t"\
		"movaps		%%xmm2,0x80(%%rbx)			\n\t				unpckhpd	0x20(%%rcx),%%xmm9 			\n\t"\
		"movaps		0x10(%%rbx),%%xmm3			\n\t				unpcklpd	0x20(%%rcx),%%xmm12			\n\t"\
		"movaps		(%%rbx),%%xmm2				\n\t				movaps		%%xmm9 ,0x20(%%rcx)			\n\t"\
		"movaps		%%xmm3,%%xmm7				\n\t				movaps		%%xmm13,0x30(%%rbx)			\n\t"\
		"movaps		%%xmm2,%%xmm6				\n\t				movaps		%%xmm12,0x20(%%rbx)			\n\t"\
		"unpckhpd	0x10(%%rcx),%%xmm7			\n\t				movaps		     (%%rax),%%xmm12			\n\t"\
		"unpcklpd	0x10(%%rcx),%%xmm3			\n\t				movaps		0x010(%%rax),%%xmm13			\n\t"\
		"movaps		%%xmm7,0x10(%%rcx)			\n\t				movaps		%%xmm12,%%xmm11				\n\t"\
		"unpckhpd	(%%rcx),%%xmm6				\n\t				movaps		%%xmm13,%%xmm9 				\n\t"\
		"unpcklpd	(%%rcx),%%xmm2				\n\t				mulpd		0x60(%%rdx),%%xmm12			\n\t"\
		"movaps		%%xmm6,    (%%rcx)			\n\t				mulpd		0x60(%%rdx),%%xmm13			\n\t"\
		"movaps		%%xmm3,0x10(%%rbx)			\n\t				mulpd		0x70(%%rdx),%%xmm11			\n\t"\
		"movaps		%%xmm2,    (%%rbx)			\n\t				mulpd		0x70(%%rdx),%%xmm9 			\n\t"\
		"addpd		%%xmm5,%%xmm0				\n\t				subpd		%%xmm11,%%xmm13				\n\t"\
		"subpd		%%xmm4,%%xmm1				\n\t				addpd		%%xmm9 ,%%xmm12				\n\t"\
		"movaps		%%xmm0,%%xmm2				\n\t				movaps		%%xmm13,%%xmm11				\n\t"\
		"movaps		%%xmm1,%%xmm3				\n\t				movaps		%%xmm12,%%xmm9 				\n\t"\
		"addpd		%%xmm5,%%xmm5				\n\t				unpckhpd	0xb0(%%rcx),%%xmm11			\n\t"\
		"addpd		%%xmm4,%%xmm4				\n\t				unpcklpd	0xb0(%%rcx),%%xmm13			\n\t"\
		"movaps		%%xmm0,%%xmm6				\n\t				movaps		%%xmm11,0xb0(%%rcx)			\n\t"\
		"movaps		%%xmm1,%%xmm7				\n\t				unpckhpd	0xa0(%%rcx),%%xmm9 			\n\t"\
		"mulpd		    (%%rdx),%%xmm2			\n\t				unpcklpd	0xa0(%%rcx),%%xmm12			\n\t"\
		"mulpd		    (%%rdx),%%xmm3			\n\t				movaps		%%xmm9 ,0xa0(%%rcx)			\n\t"\
		"mulpd		0x10(%%rdx),%%xmm6			\n\t				movaps		%%xmm13,0xb0(%%rbx)			\n\t"\
		"mulpd		0x10(%%rdx),%%xmm7			\n\t				movaps		%%xmm12,0xa0(%%rbx)			\n\t"\
		"addq	$0x20,%%rdx		/* c12 in lcol; c6 in rcol*/\n\t	subpd		%%xmm15,%%xmm8 				\n\t"\
		"subpd		%%xmm6,%%xmm3				\n\t				subpd		%%xmm14,%%xmm10				\n\t"\
		"addpd		%%xmm7,%%xmm2				\n\t				addpd		%%xmm15,%%xmm15				\n\t"\
		"movaps		%%xmm3,%%xmm7				\n\t				addpd		%%xmm14,%%xmm14				\n\t"\
		"movaps		%%xmm2,%%xmm6				\n\t				addpd		%%xmm8 ,%%xmm15				\n\t"\
		"unpckhpd	0x50(%%rcx),%%xmm7			\n\t				addpd		%%xmm10,%%xmm14				\n\t"\
		"unpcklpd	0x50(%%rcx),%%xmm3			\n\t				movaps		%%xmm15,%%xmm12				\n\t"\
		"movaps		%%xmm7,0x50(%%rcx)			\n\t				movaps		%%xmm10,%%xmm13				\n\t"\
		"unpckhpd	0x40(%%rcx),%%xmm6			\n\t				mulpd		0x60(%%rdx),%%xmm15			\n\t"\
		"unpcklpd	0x40(%%rcx),%%xmm2			\n\t				mulpd		0x60(%%rdx),%%xmm10			\n\t"\
		"movaps		%%xmm6,0x40(%%rcx)			\n\t				mulpd		0x70(%%rdx),%%xmm12			\n\t"\
		"movaps		%%xmm3,0x50(%%rbx)			\n\t				mulpd		0x70(%%rdx),%%xmm13			\n\t"\
		"movaps		%%xmm2,0x40(%%rbx)			\n\t				subpd		%%xmm12,%%xmm10				\n\t"\
		"subpd		%%xmm5,%%xmm0				\n\t				addpd		%%xmm13,%%xmm15				\n\t"\
		"addpd		%%xmm4,%%xmm1				\n\t				movaps		%%xmm10,%%xmm13				\n\t"\
		"movaps		%%xmm0,%%xmm6				\n\t				movaps		%%xmm15,%%xmm12				\n\t"\
		"movaps		%%xmm1,%%xmm7				\n\t				unpckhpd	0x70(%%rcx),%%xmm13			\n\t"\
		"mulpd		    (%%rdx),%%xmm0			\n\t				unpcklpd	0x70(%%rcx),%%xmm10			\n\t"\
		"mulpd		    (%%rdx),%%xmm1			\n\t				movaps		%%xmm13,0x70(%%rcx)			\n\t"\
		"mulpd		0x10(%%rdx),%%xmm6			\n\t				unpckhpd	0x60(%%rcx),%%xmm12			\n\t"\
		"mulpd		0x10(%%rdx),%%xmm7			\n\t				unpcklpd	0x60(%%rcx),%%xmm15			\n\t"\
		"subpd		%%xmm6,%%xmm1				\n\t				movaps		%%xmm12,0x60(%%rcx)			\n\t"\
		"addpd		%%xmm7,%%xmm0				\n\t				movaps		%%xmm10,0x70(%%rbx)			\n\t"\
		"movaps		%%xmm1,%%xmm7				\n\t				movaps		%%xmm15,0x60(%%rbx)			\n\t"\
		"movaps		%%xmm0,%%xmm6				\n\t				movaps		%%xmm8 ,%%xmm12				\n\t"\
		"unpckhpd	0xd0(%%rcx),%%xmm7			\n\t				movaps		%%xmm14,%%xmm13				\n\t"\
		"unpcklpd	0xd0(%%rcx),%%xmm1			\n\t				mulpd		0x80(%%rdx),%%xmm8 	/* c14 */\n\t"\
		"movaps		%%xmm7,0xd0(%%rcx)			\n\t				mulpd		0x80(%%rdx),%%xmm14			\n\t"\
		"unpckhpd	0xc0(%%rcx),%%xmm6			\n\t				mulpd		0x90(%%rdx),%%xmm12			\n\t"\
		"unpcklpd	0xc0(%%rcx),%%xmm0			\n\t				mulpd		0x90(%%rdx),%%xmm13			\n\t"\
		"movaps		%%xmm6,0xc0(%%rcx)			\n\t				subpd		%%xmm12,%%xmm14				\n\t"\
		"movaps		%%xmm1,0xd0(%%rbx)			\n\t				addpd		%%xmm13,%%xmm8 				\n\t"\
		"movaps		%%xmm0,0xc0(%%rbx)			\n\t				movaps		%%xmm14,%%xmm13				\n\t"\
		"															movaps		%%xmm8 ,%%xmm12				\n\t"\
		"															unpckhpd	0xf0(%%rcx),%%xmm13			\n\t"\
		"															unpcklpd	0xf0(%%rcx),%%xmm14			\n\t"\
		"															movaps		%%xmm13,0xf0(%%rcx)			\n\t"\
		"															unpckhpd	0xe0(%%rcx),%%xmm12			\n\t"\
		"															unpcklpd	0xe0(%%rcx),%%xmm8 			\n\t"\
		"															movaps		%%xmm12,0xe0(%%rcx)			\n\t"\
		"															movaps		%%xmm14,0xf0(%%rbx)			\n\t"\
		"															movaps		%%xmm8 ,0xe0(%%rbx)			\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__add1] "m" (Xadd1)\
		 ,[__r1] "m" (Xr1)\
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__pfetch_dist] "m" (Xpfetch_dist)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r14","r15","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
	);\
	}

	#define FOO(Xadd0,Xadd1,Xr1,Xisrt2)\
	{\
	__asm__ volatile (\
		"movq	%[__add1],%%rcx				\n\t							mulpd		0x060(%%rdx),%%xmm12	/* c2 */		\n\t"\
		"/*movaps	%%xmm3,%%xmm7				*/\n\t						mulpd		0x060(%%rdx),%%xmm13		\n\t"\
		"/*movaps	%%xmm2,%%xmm6				*/\n\t						mulpd		0x070(%%rdx),%%xmm11		\n\t"\
		"/*unpckhpd	0x090(%%rcx),%%xmm7			*/\n\t						mulpd		0x070(%%rdx),%%xmm9 		\n\t"\
		"/*unpcklpd	0x090(%%rcx),%%xmm3			*/\n\t						subpd		%%xmm11     ,%%xmm13		\n\t"\
		"/*movaps	%%xmm7,0x090(%%rcx)			*/\n\t						addpd		%%xmm9      ,%%xmm12		\n\t"\
		"/*addq	$0x020,%%rdx	 		c4 in lcol; c10 in rcol */\n\t		/*movaps		%%xmm13,%%xmm11				*/\n\t"\
		"/*unpckhpd	0x080(%%rcx),%%xmm6			*/\n\t						/*movaps		%%xmm12,%%xmm9 				*/\n\t"\
		"/*unpcklpd	0x080(%%rcx),%%xmm2			*/\n\t						/*unpckhpd	0x030(%%rcx),%%xmm11		*/\n\t"\
		"/*movaps	%%xmm6,0x080(%%rcx)			*/\n\t						/*unpcklpd	0x030(%%rcx),%%xmm13		*/\n\t"\
		"/*movaps	%%xmm3,0x090(%%rbx)			*/\n\t						/*movaps		%%xmm11,0x030(%%rcx)		*/\n\t"\
		"/*movaps	%%xmm2,0x080(%%rbx)			*/\n\t						/*unpckhpd	0x020(%%rcx),%%xmm9 		*/\n\t"\
"\n\t"/* AVX debug: everywhere we have an unpckhpd/unpcklpd pair followed by a pair of writes-to-memory, copy the same data into local memory: */\
"movaps	0x090(%%rcx),%%xmm7		\n\t"\
"movaps	0x080(%%rcx),%%xmm6		\n\t"\
"movaps	%%xmm3,0x090(%%rax)		\n\t"\
"movaps	%%xmm2,0x080(%%rax)		\n\t"\
"movaps	%%xmm7,0x190(%%rax)		\n\t"\
"movaps	%%xmm6,0x180(%%rax)		\n\t"\
		"/*movaps	0x010(%%rbx),%%xmm3			*/\n\t						/*unpcklpd	0x020(%%rcx),%%xmm12		*/\n\t"\
		"/*movaps	     (%%rbx),%%xmm2			*/\n\t						/*movaps		%%xmm9 ,0x020(%%rcx)		*/\n\t"\
		"/*movaps	%%xmm3,%%xmm7				*/\n\t						/*movaps		%%xmm13,0x030(%%rbx)		*/\n\t"\
		"/*movaps	%%xmm2,%%xmm6				*/\n\t						/*movaps		%%xmm12,0x020(%%rbx)		*/\n\t"\
																	"movaps	0x030(%%rcx),%%xmm11	\n\t"\
																	"movaps	0x020(%%rcx),%%xmm9		\n\t"\
																	"movaps	%%xmm13,0x030(%%rax)	\n\t"\
																	"movaps	%%xmm12,0x020(%%rax)	\n\t"\
																	"movaps	%%xmm11,0x130(%%rax)	\n\t"\
																	"movaps	%%xmm9 ,0x120(%%rax)	\n\t"\
		"																	addq	$0x020,%%rdx		/* c4 in lcol; c10 in rcol*/		\n\t"\
		"/*unpckhpd	0x010(%%rcx),%%xmm7			*/\n\t						movaps		     (%%rax),%%xmm12		\n\t"\
		"/*unpcklpd	0x010(%%rcx),%%xmm3			*/\n\t						movaps		0x010(%%rax),%%xmm13		\n\t"\
		"/*movaps	%%xmm7,0x010(%%rcx)			*/\n\t						movaps		%%xmm12,%%xmm11		\n\t"\
		"/*unpckhpd	(%%rcx),%%xmm6				*/\n\t						movaps		%%xmm13,%%xmm9 		\n\t"\
		"/*unpcklpd	(%%rcx),%%xmm2				*/\n\t						mulpd		0x060(%%rdx),%%xmm12		\n\t"\
		"/*movaps	%%xmm6,    (%%rcx)			*/\n\t						mulpd		0x060(%%rdx),%%xmm13		\n\t"\
		"/*movaps	%%xmm3,0x010(%%rbx)			*/\n\t						mulpd		0x070(%%rdx),%%xmm11		\n\t"\
		"/*movaps	%%xmm2,    (%%rbx)			*/\n\t						mulpd		0x070(%%rdx),%%xmm9 		\n\t"\
"movaps	0x010(%%rbx),%%xmm3		\n\t"\
"movaps	0x000(%%rbx),%%xmm2		\n\t"\
"movaps	0x010(%%rcx),%%xmm7		\n\t"\
"movaps	0x000(%%rcx),%%xmm6		\n\t"\
"movaps	%%xmm3,0x010(%%rax)		\n\t"\
"movaps	%%xmm2,0x000(%%rax)		\n\t"\
"movaps	%%xmm7,0x110(%%rax)		\n\t"\
"movaps	%%xmm6,0x100(%%rax)		\n\t"\
		"addpd		%%xmm5,%%xmm0		\n\t						subpd		%%xmm11     ,%%xmm13		\n\t"\
		"subpd		%%xmm4,%%xmm1		\n\t						addpd		%%xmm9      ,%%xmm12		\n\t"\
		"movaps	%%xmm0,%%xmm2		\n\t						/*movaps		%%xmm13,%%xmm11				*/\n\t"\
		"movaps	%%xmm1,%%xmm3		\n\t						/*movaps		%%xmm12,%%xmm9 				*/\n\t"\
		"addpd		%%xmm5,%%xmm5		\n\t						/*unpckhpd	0x0b0(%%rcx),%%xmm11		*/\n\t"\
		"addpd		%%xmm4,%%xmm4		\n\t						/*unpcklpd	0x0b0(%%rcx),%%xmm13		*/\n\t"\
		"movaps	%%xmm0,%%xmm6				\n\t						/*movaps		%%xmm11,0x0b0(%%rcx)		*/\n\t"\
		"movaps	%%xmm1,%%xmm7				\n\t						/*unpckhpd	0x0a0(%%rcx),%%xmm9 		*/\n\t"\
		"mulpd		     (%%rdx),%%xmm2	\n\t						/*unpcklpd	0x0a0(%%rcx),%%xmm12		*/\n\t"\
		"mulpd		     (%%rdx),%%xmm3	\n\t						/*movaps		%%xmm9 ,0x0a0(%%rcx)		*/\n\t"\
		"mulpd		0x010(%%rdx),%%xmm6	\n\t						/*movaps		%%xmm13,0x0b0(%%rbx)		*/\n\t"\
		"mulpd		0x010(%%rdx),%%xmm7	\n\t						/*movaps		%%xmm12,0x0a0(%%rbx)		*/\n\t"\
																	"movaps	0x0b0(%%rcx),%%xmm11	\n\t"\
																	"movaps	0x0a0(%%rcx),%%xmm9		\n\t"\
																	"movaps	%%xmm13,0x0b0(%%rax)	\n\t"\
																	"movaps	%%xmm12,0x0a0(%%rax)	\n\t"\
																	"movaps	%%xmm11,0x1b0(%%rax)	\n\t"\
																	"movaps	%%xmm9 ,0x1a0(%%rax)	\n\t"\
		"addq	$0x020,%%rdx	/* c12 in lcol; c6 in rcol*/\n\t	subpd		%%xmm15,%%xmm8 		\n\t"\
		"subpd		%%xmm6,%%xmm3		\n\t						subpd		%%xmm14,%%xmm10		\n\t"\
		"addpd		%%xmm7,%%xmm2		\n\t						addpd		%%xmm15,%%xmm15		\n\t"\
		"/*movaps	%%xmm3,%%xmm7				*/\n\t				addpd		%%xmm14,%%xmm14		\n\t"\
		"/*movaps	%%xmm2,%%xmm6				*/\n\t				addpd		%%xmm8 ,%%xmm15		\n\t"\
		"/*unpckhpd	0x050(%%rcx),%%xmm7			*/\n\t				addpd		%%xmm10,%%xmm14		\n\t"\
		"/*unpcklpd	0x050(%%rcx),%%xmm3			*/\n\t				movaps		%%xmm15,%%xmm12		\n\t"\
		"/*movaps	%%xmm7,0x050(%%rcx)			*/\n\t				movaps		%%xmm10,%%xmm13		\n\t"\
		"/*unpckhpd	0x040(%%rcx),%%xmm6			*/\n\t				mulpd		0x060(%%rdx),%%xmm15		\n\t"\
		"/*unpcklpd	0x040(%%rcx),%%xmm2			*/\n\t				mulpd		0x060(%%rdx),%%xmm10		\n\t"\
		"/*movaps	%%xmm6,0x040(%%rcx)			*/\n\t				mulpd		0x070(%%rdx),%%xmm12		\n\t"\
		"/*movaps	%%xmm3,0x050(%%rbx)			*/\n\t				mulpd		0x070(%%rdx),%%xmm13		\n\t"\
		"/*movaps	%%xmm2,0x040(%%rbx)			*/\n\t				subpd		%%xmm12     ,%%xmm10		\n\t"\
"movaps	0x050(%%rcx),%%xmm7		\n\t"\
"movaps	0x040(%%rcx),%%xmm6		\n\t"\
"movaps	%%xmm3,0x050(%%rax)		\n\t"\
"movaps	%%xmm2,0x040(%%rax)		\n\t"\
"movaps	%%xmm7,0x150(%%rax)		\n\t"\
"movaps	%%xmm6,0x140(%%rax)		\n\t"\
		"subpd		%%xmm5,%%xmm0		\n\t						addpd		%%xmm13     ,%%xmm15		\n\t"\
		"addpd		%%xmm4,%%xmm1		\n\t						/*movaps		%%xmm10,%%xmm13			*/\n\t"\
		"movaps	%%xmm0,%%xmm6				\n\t					/*movaps		%%xmm15,%%xmm12			*/\n\t"\
		"movaps	%%xmm1,%%xmm7				\n\t					/*unpckhpd	0x070(%%rcx),%%xmm13		*/\n\t"\
		"mulpd		     (%%rdx),%%xmm0	\n\t						/*unpcklpd	0x070(%%rcx),%%xmm10		*/\n\t"\
		"mulpd		     (%%rdx),%%xmm1	\n\t						/*movaps		%%xmm13,0x070(%%rcx)	*/\n\t"\
		"mulpd		0x010(%%rdx),%%xmm6	\n\t						/*unpckhpd	0x060(%%rcx),%%xmm12		*/\n\t"\
		"mulpd		0x010(%%rdx),%%xmm7	\n\t						/*unpcklpd	0x060(%%rcx),%%xmm15		*/\n\t"\
		"subpd		%%xmm6,%%xmm1		\n\t						/*movaps		%%xmm12,0x060(%%rcx)	*/\n\t"\
		"addpd		%%xmm7,%%xmm0		\n\t						/*movaps		%%xmm10,0x070(%%rbx)	*/\n\t"\
		"/*movaps	%%xmm1,%%xmm7		*/\n\t						/*movaps		%%xmm15,0x060(%%rbx)	*/\n\t"\
																	"movaps	0x070(%%rcx),%%xmm13	\n\t"\
																	"movaps	0x060(%%rcx),%%xmm12	\n\t"\
																	"movaps	%%xmm10,0x070(%%rax)	\n\t"\
																	"movaps	%%xmm15,0x060(%%rax)	\n\t"\
																	"movaps	%%xmm13,0x170(%%rax)	\n\t"\
																	"movaps	%%xmm12,0x160(%%rax)	\n\t"\
		"/*movaps	%%xmm0,%%xmm6				*/\n\t						movaps		%%xmm8 ,%%xmm12		\n\t"\
		"/*unpckhpd	0x0d0(%%rcx),%%xmm7			*/\n\t						movaps		%%xmm14,%%xmm13		\n\t"\
		"/*unpcklpd	0x0d0(%%rcx),%%xmm1			*/\n\t						mulpd		0x080(%%rdx),%%xmm8 	/* c14 */		\n\t"\
		"/*movaps	%%xmm7,0x0d0(%%rcx)			*/\n\t						mulpd		0x080(%%rdx),%%xmm14		\n\t"\
		"/*unpckhpd	0x0c0(%%rcx),%%xmm6			*/\n\t						mulpd		0x090(%%rdx),%%xmm12		\n\t"\
		"/*unpcklpd	0x0c0(%%rcx),%%xmm0			*/\n\t						mulpd		0x090(%%rdx),%%xmm13		\n\t"\
		"/*movaps	%%xmm6,0x0c0(%%rcx)			*/\n\t						subpd		%%xmm12     ,%%xmm14		\n\t"\
		"/*movaps	%%xmm1,0x0d0(%%rbx)			*/\n\t						addpd		%%xmm13     ,%%xmm8 		\n\t"\
		"/*movaps	%%xmm0,0x0c0(%%rbx)			*/\n\t						/*movaps		%%xmm14,%%xmm13				*/\n\t"\
"movaps	0x0d0(%%rcx),%%xmm7		\n\t"\
"movaps	0x0c0(%%rcx),%%xmm6		\n\t"\
"movaps	%%xmm1,0x0d0(%%rax)		\n\t"\
"movaps	%%xmm0,0x0c0(%%rax)		\n\t"\
"movaps	%%xmm7,0x1d0(%%rax)		\n\t"\
"movaps	%%xmm6,0x1c0(%%rax)		\n\t"\
		"																	/*movaps		%%xmm8 ,%%xmm12				*/\n\t"\
		"																	/*unpckhpd	0x0f0(%%rcx),%%xmm13		*/\n\t"\
		"																	/*unpcklpd	0x0f0(%%rcx),%%xmm14		*/\n\t"\
		"																	/*movaps		%%xmm13,0x0f0(%%rcx)		*/\n\t"\
		"																	/*unpckhpd	0x0e0(%%rcx),%%xmm12		*/\n\t"\
		"																	/*unpcklpd	0x0e0(%%rcx),%%xmm8 		*/\n\t"\
		"																	/*movaps		%%xmm12,0x0e0(%%rcx)		*/\n\t"\
		"																	/*movaps		%%xmm14,0x0f0(%%rbx)		*/\n\t"\
		"																	/*movaps		%%xmm8 ,0x0e0(%%rbx)		*/\n\t"\
																	"movaps	0x0f0(%%rcx),%%xmm13	\n\t"\
																	"movaps	0x0e0(%%rcx),%%xmm12	\n\t"\
																	"movaps	%%xmm14,0x0f0(%%rax)	\n\t"\
																	"movaps	%%xmm8 ,0x0e0(%%rax)	\n\t"\
																	"movaps	%%xmm13,0x1f0(%%rax)	\n\t"\
																	"movaps	%%xmm12,0x1e0(%%rax)	\n\t"\
		"\n\t"\
		"movq	%[__isrt2],%%rbx	/* DEBUG: Dump register contents into trig-data slots: */\n\t"\
		"movaps	%%xmm0 ,0x300(%%rbx)	\n\t"\
		"movaps	%%xmm1 ,0x310(%%rbx)	\n\t"\
		"movaps	%%xmm2 ,0x320(%%rbx)	\n\t"\
		"movaps	%%xmm3 ,0x330(%%rbx)	\n\t"\
		"movaps	%%xmm4 ,0x340(%%rbx)	\n\t"\
		"movaps	%%xmm5 ,0x350(%%rbx)	\n\t"\
		"movaps	%%xmm6 ,0x360(%%rbx)	\n\t"\
		"movaps	%%xmm7 ,0x370(%%rbx)	\n\t"\
		"movaps	%%xmm8 ,0x380(%%rbx)	\n\t"\
		"movaps	%%xmm9 ,0x390(%%rbx)	\n\t"\
		"movaps	%%xmm10,0x3a0(%%rbx)	\n\t"\
		"movaps	%%xmm11,0x3b0(%%rbx)	\n\t"\
		"movaps	%%xmm12,0x3c0(%%rbx)	\n\t"\
		"movaps	%%xmm13,0x3d0(%%rbx)	\n\t"\
		"movaps	%%xmm14,0x3e0(%%rbx)	\n\t"\
		"movaps	%%xmm15,0x3f0(%%rbx)	\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__add1] "m" (Xadd1)\
		 ,[__r1] "m" (Xr1)\
		 ,[__isrt2] "m" (Xisrt2)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
	);\
	}

	#define BAR(Xadd0,Xadd1,Xr1,Xisrt2)\
	{\
	__asm__ volatile (\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__add1] "m" (Xadd1)\
		 ,[__r1] "m" (Xr1)\
		 ,[__isrt2] "m" (Xisrt2)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
	);\
	}

#endif	// SSE2 or AVX?

#endif	/* radix16_dyadic_square_gcc_h_included */

