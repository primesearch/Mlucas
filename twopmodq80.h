/*******************************************************************************
*                                                                              *
*   (C) 1997-2012 by Ernst W. Mayer.                                           *
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
#ifndef twopmodq80_h_included
#define twopmodq80_h_included

#include "masterdefs.h"
#include "types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* 4-operand SSE2 and 8-operand AVX modmul require 64-bit GCC build: */

#if defined(USE_SSE2) && defined(COMPILER_TYPE_GCC) && (OS_BITS == 64)

  #ifdef USE_AVX

	#define SSE2_twopmodq78_modmul_q16(Xaq0,Xaqinv0,Xax0,Xtwo26i,Xpshift,Xj)\
	{\
	__asm__ volatile (\
	/* SQR_LOHI78_3WORD_DOUBLE_q16(). Inputs [a|c|e|g]lo0-2 enter in ymm[0-2,4-6,8-10,12-14], resp.: */\
		"movq	%[__ax0],%%rax		\n\t"\
		"movq	%[__two26i],%%rbx	\n\t"\
		"vmovaps	-0x40(%%rbx),%%ymm3	\n\t"/* two13i */\
			/*** STREAM 0 ***/							/*** STREAM 1 ***/							/*** STREAM 2 ***/							/*** STREAM 3 ***/					\
		"/* fx0,1,2 assumed in ymm0,1,2 on loop entry */\n\t"\
			"vmulpd	%%ymm3,%%ymm0,%%ymm0		\n\t	vmulpd	%%ymm3,%%ymm4,%%ymm4		\n\t	vmulpd	%%ymm3,%%ymm8,%%ymm8		\n\t	vmulpd	%%ymm3,%%ymm12,%%ymm12		\n\t"/* scale [fx0 . 2^-13] */\
			"vmulpd	%%ymm3,%%ymm1,%%ymm1		\n\t	vmulpd	%%ymm3,%%ymm5,%%ymm5		\n\t	vmulpd	%%ymm3,%%ymm9,%%ymm9		\n\t	vmulpd	%%ymm3,%%ymm13,%%ymm13		\n\t"/* scale [fx1 . 2^-13] */\
			"vmulpd	%%ymm3,%%ymm2,%%ymm2		\n\t	vmulpd	%%ymm3,%%ymm6,%%ymm6		\n\t	vmulpd	%%ymm3,%%ymm10,%%ymm10		\n\t	vmulpd	%%ymm3,%%ymm14,%%ymm14		\n\t"/* scale [fx2 . 2^-13] */\
			"vmovaps	%%ymm0,%%ymm3			\n\t	vmovaps	%%ymm4,%%ymm7				\n\t	vmovaps	%%ymm8,%%ymm11				\n\t	vmovaps	%%ymm12,%%ymm15				\n\t"/* cpy of fx0 */\
			"vaddpd	%%ymm3,%%ymm3,%%ymm3		\n\t	vaddpd	%%ymm7,%%ymm7,%%ymm7		\n\t	vaddpd	%%ymm11,%%ymm11,%%ymm11		\n\t	vaddpd	%%ymm15,%%ymm15,%%ymm15		\n\t"/* 2.fx0 */\
			"vmovaps	%%ymm2,0x100(%%rax)		\n\t	vmovaps	%%ymm6,0x120(%%rax)			\n\t	vmovaps	%%ymm10,0x140(%%rax)		\n\t	vmovaps	%%ymm14,0x160(%%rax)		\n\t"/* Store fx2 to free up a register */\
			"/* Digit 0: */\n\t"\
			"vmulpd	%%ymm0,%%ymm0,%%ymm0		\n\t	vmulpd	%%ymm4,%%ymm4,%%ymm4		\n\t	vmulpd	%%ymm8,%%ymm8,%%ymm8		\n\t	vmulpd	%%ymm12,%%ymm12,%%ymm12		\n\t"/* [  fx0 *= fx0] / 2^26 */\
			"vroundpd $0,%%ymm0,%%ymm2			\n\t	vroundpd $0,%%ymm4,%%ymm6			\n\t	vroundpd $0,%%ymm8,%%ymm10			\n\t	vroundpd $0,%%ymm12,%%ymm14			\n\t"/* Init: FCY = dnint(fx0*fx0) */\
			"vsubpd	%%ymm2,%%ymm0,%%ymm0		\n\t	vsubpd	%%ymm6,%%ymm4,%%ymm4		\n\t	vsubpd	%%ymm10,%%ymm8,%%ymm8		\n\t	vsubpd	%%ymm14,%%ymm12,%%ymm12		\n\t"/* fprod0 -= FCY */\
			"vmulpd	-0x20(%%rbx),%%ymm0,%%ymm0	\n\t	vmulpd	-0x20(%%rbx),%%ymm4,%%ymm4	\n\t	vmulpd	-0x20(%%rbx),%%ymm8,%%ymm8	\n\t	vmulpd	-0x20(%%rbx),%%ymm12,%%ymm12\n\t"/* fprod0 *= two26f */\
			"vmulpd	     (%%rbx),%%ymm2,%%ymm2	\n\t	vmulpd	     (%%rbx),%%ymm6,%%ymm6	\n\t	vmulpd	     (%%rbx),%%ymm10,%%ymm10\n\t	vmulpd	     (%%rbx),%%ymm14,%%ymm14\n\t"/* FCY    *= two26i */\
			"vmovaps	%%ymm0,    (%%rax)		\n\t	vmovaps	%%ymm4,0x20(%%rax)			\n\t	vmovaps	%%ymm8,0x40(%%rax)			\n\t	vmovaps	%%ymm12,0x60(%%rax)			\n\t"/* Store fprod0 to free up a register */\
			"/* Digit 1: */\n\t"\
			"vmovaps	%%ymm1,%%ymm0			\n\t	vmovaps	%%ymm5,%%ymm4				\n\t	vmovaps	%%ymm9,%%ymm8				\n\t	vmovaps	%%ymm13,%%ymm12				\n\t"/* cpy of fx1 */\
			"vmulpd	%%ymm3,%%ymm0,%%ymm0		\n\t	vmulpd	%%ymm7,%%ymm4,%%ymm4		\n\t	vmulpd	%%ymm11,%%ymm8,%%ymm8		\n\t	vmulpd	%%ymm15,%%ymm12,%%ymm12		\n\t"/* [fx1 *= 2.fx0] / 2^26 */\
			"vaddpd	%%ymm2,%%ymm0,%%ymm0		\n\t	vaddpd	%%ymm6,%%ymm4,%%ymm4		\n\t	vaddpd	%%ymm10,%%ymm8,%%ymm8		\n\t	vaddpd	%%ymm14,%%ymm12,%%ymm12		\n\t"/* fprod1 = 2.fx0*fx1 + FCY */\
			"vroundpd $0,%%ymm0,%%ymm2			\n\t	vroundpd $0,%%ymm4,%%ymm6			\n\t	vroundpd $0,%%ymm8,%%ymm10			\n\t	vroundpd $0,%%ymm12,%%ymm14			\n\t"/* FCY */\
			"vsubpd	%%ymm2,%%ymm0,%%ymm0		\n\t	vsubpd	%%ymm6,%%ymm4,%%ymm4		\n\t	vsubpd	%%ymm10,%%ymm8,%%ymm8		\n\t	vsubpd	%%ymm14,%%ymm12,%%ymm12		\n\t"/* fprod1 -= FCY */\
			"vmulpd	-0x20(%%rbx),%%ymm0,%%ymm0	\n\t	vmulpd	-0x20(%%rbx),%%ymm4,%%ymm4	\n\t	vmulpd	-0x20(%%rbx),%%ymm8,%%ymm8	\n\t	vmulpd	-0x20(%%rbx),%%ymm12,%%ymm12\n\t"/* fprod1 *= two26f */\
			"vmulpd	     (%%rbx),%%ymm2,%%ymm2	\n\t	vmulpd	     (%%rbx),%%ymm6,%%ymm6	\n\t	vmulpd	     (%%rbx),%%ymm10,%%ymm10\n\t	vmulpd	     (%%rbx),%%ymm14,%%ymm14\n\t"/* FCY    *= two26i */\
			"vmovaps	%%ymm0,0x80(%%rax)		\n\t	vmovaps	%%ymm4,0xa0(%%rax)			\n\t	vmovaps	%%ymm8,0xc0(%%rax)			\n\t	vmovaps	%%ymm12,0xe0(%%rax)			\n\t"/* Store fprod1 to free up a register */\
			"/* Digit 2: */\n\t"\
			"vmovaps	%%ymm1,%%ymm0			\n\t	vmovaps	%%ymm5,%%ymm4				\n\t	vmovaps	%%ymm9,%%ymm8				\n\t	vmovaps	%%ymm13,%%ymm12				\n\t"/* cpy of fx1 */\
			"vaddpd	%%ymm0,%%ymm0,%%ymm0		\n\t	vaddpd	%%ymm4,%%ymm4,%%ymm4		\n\t	vaddpd	%%ymm8,%%ymm8,%%ymm8		\n\t	vaddpd	%%ymm12,%%ymm12,%%ymm12		\n\t"/* 2.fx1 */\
			"vmulpd	%%ymm1,%%ymm1,%%ymm1		\n\t	vmulpd	%%ymm5,%%ymm5,%%ymm5		\n\t	vmulpd	%%ymm9,%%ymm9,%%ymm9		\n\t	vmulpd	%%ymm13,%%ymm13,%%ymm13		\n\t"/* [  fx1 *= fx1] / 2^26 */\
			"vmulpd	0x100(%%rax),%%ymm3,%%ymm3	\n\t	vmulpd	0x120(%%rax),%%ymm7,%%ymm7	\n\t	vmulpd	0x140(%%rax),%%ymm11,%%ymm11\n\t	vmulpd	0x160(%%rax),%%ymm15,%%ymm15\n\t"/* [2.fx0 *= fx2] / 2^26 */\
			"vaddpd	%%ymm2,%%ymm1,%%ymm1		\n\t	vaddpd	%%ymm6,%%ymm5,%%ymm5		\n\t	vaddpd	%%ymm10,%%ymm9,%%ymm9		\n\t	vaddpd	%%ymm14,%%ymm13,%%ymm13		\n\t"/* fx1*fx1 += FCY */\
			"vaddpd	%%ymm3,%%ymm1,%%ymm1		\n\t	vaddpd	%%ymm7,%%ymm5,%%ymm5		\n\t	vaddpd	%%ymm11,%%ymm9,%%ymm9		\n\t	vaddpd	%%ymm15,%%ymm13,%%ymm13		\n\t"/* fprod2 = 2.fx0*fx2 + fx1*fx1 + FCY */\
			"vroundpd $1,%%ymm1,%%ymm2			\n\t	vroundpd $1,%%ymm5,%%ymm6			\n\t	vroundpd $1,%%ymm9,%%ymm10			\n\t	vroundpd $1,%%ymm13,%%ymm14			\n\t"/* FCY = floor(fprod2) */\
			"vsubpd	%%ymm2,%%ymm1,%%ymm1		\n\t	vsubpd	%%ymm6,%%ymm5,%%ymm5		\n\t	vsubpd	%%ymm10,%%ymm9,%%ymm9		\n\t	vsubpd	%%ymm14,%%ymm13,%%ymm13		\n\t"/* fprod2 -= FCY */\
			"vmulpd	-0x20(%%rbx),%%ymm1,%%ymm1	\n\t	vmulpd	-0x20(%%rbx),%%ymm5,%%ymm5	\n\t	vmulpd	-0x20(%%rbx),%%ymm9,%%ymm9	\n\t	vmulpd	-0x20(%%rbx),%%ymm13,%%ymm13\n\t"/* fprod2 *= two26f */\
			"vmulpd	     (%%rbx),%%ymm2,%%ymm2	\n\t	vmulpd	     (%%rbx),%%ymm6,%%ymm6	\n\t	vmulpd	     (%%rbx),%%ymm10,%%ymm10\n\t	vmulpd	     (%%rbx),%%ymm14,%%ymm14\n\t"/* FCY    *= two26i */\
			"/* Digit 3: */\n\t"\
			"vmovaps	0x100(%%rax),%%ymm3		\n\t	vmovaps	0x120(%%rax),%%ymm7			\n\t	vmovaps	0x140(%%rax),%%ymm11		\n\t	vmovaps	0x160(%%rax),%%ymm15		\n\t"/* Reload fx2 */\
			"vmulpd	%%ymm3,%%ymm0,%%ymm0		\n\t	vmulpd	%%ymm7,%%ymm4,%%ymm4		\n\t	vmulpd	%%ymm11,%%ymm8,%%ymm8		\n\t	vmulpd	%%ymm15,%%ymm12,%%ymm12		\n\t"/* [2.fx1 *= fx2] / 2^26 */\
			"vaddpd	%%ymm2,%%ymm0,%%ymm0		\n\t	vaddpd	%%ymm6,%%ymm4,%%ymm4		\n\t	vaddpd	%%ymm10,%%ymm8,%%ymm8		\n\t	vaddpd	%%ymm14,%%ymm12,%%ymm12		\n\t"/* fprod3 = 2.fx1*fx2 + FCY */\
			"vroundpd $0,%%ymm0,%%ymm2			\n\t	vroundpd $0,%%ymm4,%%ymm6			\n\t	vroundpd $0,%%ymm8,%%ymm10			\n\t	vroundpd $0,%%ymm12,%%ymm14			\n\t"\
			"vsubpd	%%ymm2,%%ymm0,%%ymm0		\n\t	vsubpd	%%ymm6,%%ymm4,%%ymm4		\n\t	vsubpd	%%ymm10,%%ymm8,%%ymm8		\n\t	vsubpd	%%ymm14,%%ymm12,%%ymm12		\n\t"/* fprod3 -= FCY */\
			"vmulpd	-0x20(%%rbx),%%ymm0,%%ymm0	\n\t	vmulpd	-0x20(%%rbx),%%ymm4,%%ymm4	\n\t	vmulpd	-0x20(%%rbx),%%ymm8,%%ymm8	\n\t	vmulpd	-0x20(%%rbx),%%ymm12,%%ymm12\n\t"/* fprod3 *= two26f */\
			"vmulpd	     (%%rbx),%%ymm2,%%ymm2	\n\t	vmulpd	     (%%rbx),%%ymm6,%%ymm6	\n\t	vmulpd	     (%%rbx),%%ymm10,%%ymm10\n\t	vmulpd	     (%%rbx),%%ymm14,%%ymm14\n\t"/* FCY    *= two26i */\
			"/* Digit 4: */\n\t"\
			"vmulpd	%%ymm3,%%ymm3,%%ymm3		\n\t	vmulpd	%%ymm7,%%ymm7,%%ymm7		\n\t	vmulpd	%%ymm11,%%ymm11,%%ymm11		\n\t	vmulpd	%%ymm15,%%ymm15,%%ymm15		\n\t"/* [  fx2 *= fx2] / 2^26 */\
			"vaddpd	%%ymm2,%%ymm3,%%ymm3		\n\t	vaddpd	%%ymm6,%%ymm7,%%ymm7		\n\t	vaddpd	%%ymm10,%%ymm11,%%ymm11		\n\t	vaddpd	%%ymm14,%%ymm15,%%ymm15		\n\t"/* fprod4 = fx2*fx2 + fcy */\
			"vroundpd $0,%%ymm3,%%ymm2			\n\t	vroundpd $0,%%ymm7,%%ymm6			\n\t	vroundpd $0,%%ymm11,%%ymm10			\n\t	vroundpd $0,%%ymm15,%%ymm14			\n\t"\
			"vsubpd	%%ymm2,%%ymm3,%%ymm3		\n\t	vsubpd	%%ymm6,%%ymm7,%%ymm7		\n\t	vsubpd	%%ymm10,%%ymm11,%%ymm11		\n\t	vsubpd	%%ymm14,%%ymm15,%%ymm15		\n\t"/* fprod4 -= fcy */\
			"vmulpd	-0x20(%%rbx),%%ymm3,%%ymm3	\n\t	vmulpd	-0x20(%%rbx),%%ymm7,%%ymm7	\n\t	vmulpd	-0x20(%%rbx),%%ymm11,%%ymm11\n\t	vmulpd	-0x20(%%rbx),%%ymm15,%%ymm15\n\t"/* fprod4 *= two26f */\
			"/* Digit 5 = the carry. CY in ymm2; fhi0,1 in ymm0,ymm3 */\n\t"\
			"vmulpd	-0x20(%%rbx),%%ymm2,%%ymm2	\n\t	vmulpd	-0x20(%%rbx),%%ymm6,%%ymm6	\n\t	vmulpd	-0x20(%%rbx),%%ymm10,%%ymm10\n\t	vmulpd	-0x20(%%rbx),%%ymm14,%%ymm14\n\t"/* fhi2 * two26f */\
			"vaddpd	%%ymm2,%%ymm3,%%ymm3		\n\t	vaddpd	%%ymm6,%%ymm7,%%ymm7		\n\t	vaddpd	%%ymm10,%%ymm11,%%ymm11		\n\t	vaddpd	%%ymm14,%%ymm15,%%ymm15		\n\t"/* fhi, top 52 bits; ymm2 FREE */\
			"vmovaps	%%ymm0,0x180(%%rax)		\n\t	vmovaps	%%ymm4,0x1a0(%%rax)			\n\t	vmovaps	%%ymm8,0x1c0(%%rax)			\n\t	vmovaps	%%ymm12,0x1e0(%%rax)		\n\t"/* Store fhi0,1 = fprod3,4 to free up 2 more registers */\
			"vmovaps	%%ymm3,0x200(%%rax)		\n\t	vmovaps	%%ymm7,0x220(%%rax)			\n\t	vmovaps	%%ymm11,0x240(%%rax)		\n\t	vmovaps	%%ymm15,0x260(%%rax)		\n\t"/* Recall: fx0,1 in eax+0x[0,2]0, fx2 in ymm1. */\
		/* MULL78 section below needs flo0,1,2 in ymm0,1,2: */\
		/* flo = MULL78_3WORD_DOUBLE_q2(flo, fqinv): */\
			"movq	%[__aqinv0],%%rdx	\n\t"\
			"/* Digit 0: */\n\t"\
			"vmovaps	    (%%rax),%%ymm0		\n\t	vmovaps	0x20(%%rax),%%ymm4			\n\t	vmovaps	0x40(%%rax),%%ymm8			\n\t	vmovaps	0x60(%%rax),%%ymm12			\n\t"/* Reload flo0 */\
			"vmovaps	%%ymm1,%%ymm2			\n\t	vmovaps	%%ymm5,%%ymm6				\n\t	vmovaps	%%ymm9,%%ymm10				\n\t	vmovaps	%%ymm13,%%ymm14				\n\t"/* Move flo2 into ymm2 */\
			"vmovaps	%%ymm0,%%ymm3			\n\t	vmovaps	%%ymm4,%%ymm7				\n\t	vmovaps	%%ymm8,%%ymm11				\n\t	vmovaps	%%ymm12,%%ymm15				\n\t"/* ymm3 = cpy of x0 */\
			"vmulpd	    (%%rdx),%%ymm0,%%ymm0	\n\t	vmulpd	0x20(%%rdx),%%ymm4,%%ymm4	\n\t	vmulpd	0x40(%%rdx),%%ymm8,%%ymm8	\n\t	vmulpd	0x60(%%rdx),%%ymm12,%%ymm12	\n\t"/* fprod0 = x0*y0 */\
			"vroundpd $0,%%ymm0,%%ymm1			\n\t	vroundpd $0,%%ymm4,%%ymm5			\n\t	vroundpd $0,%%ymm8,%%ymm9			\n\t	vroundpd $0,%%ymm12,%%ymm13			\n\t"\
			"vsubpd	%%ymm1,%%ymm0,%%ymm0		\n\t	vsubpd	%%ymm5,%%ymm4,%%ymm4		\n\t	vsubpd	%%ymm9,%%ymm8,%%ymm8		\n\t	vsubpd	%%ymm13,%%ymm12,%%ymm12		\n\t"/* fprod0 -= fcy */\
			"vmulpd	-0x20(%%rbx),%%ymm0,%%ymm0	\n\t	vmulpd	-0x20(%%rbx),%%ymm4,%%ymm4	\n\t	vmulpd	-0x20(%%rbx),%%ymm8,%%ymm8	\n\t	vmulpd	-0x20(%%rbx),%%ymm12,%%ymm12\n\t"/* fprod0 *= two26f */\
			"vmulpd	     (%%rbx),%%ymm1,%%ymm1	\n\t	vmulpd	     (%%rbx),%%ymm5,%%ymm5	\n\t	vmulpd	     (%%rbx),%%ymm9,%%ymm9	\n\t	vmulpd	     (%%rbx),%%ymm13,%%ymm13\n\t"/* fcy    *= two26i */\
			"vmovaps	%%ymm0,    (%%rax)		\n\t	vmovaps	%%ymm4,0x20(%%rax)			\n\t	vmovaps	%%ymm8,0x40(%%rax)			\n\t	vmovaps	%%ymm12,0x60(%%rax)			\n\t"/* Store fprod0 to free up a register */\
			"vmovaps	%%ymm3,%%ymm0			\n\t	vmovaps	%%ymm7,%%ymm4				\n\t	vmovaps	%%ymm11,%%ymm8				\n\t	vmovaps	%%ymm15,%%ymm12				\n\t"/* ymm0 = cpy of x0 */\
			"/* Digit 1: */\n\t"\
			"vmulpd	0x80(%%rdx),%%ymm3,%%ymm3	\n\t	vmulpd	0xa0(%%rdx),%%ymm7,%%ymm7	\n\t	vmulpd	0xc0(%%rdx),%%ymm11,%%ymm11	\n\t	vmulpd	0xe0(%%rdx),%%ymm15,%%ymm15	\n\t"/* x0*y1 */\
			"vaddpd	%%ymm1,%%ymm3,%%ymm3		\n\t	vaddpd	%%ymm5,%%ymm7,%%ymm7		\n\t	vaddpd	%%ymm9,%%ymm11,%%ymm11		\n\t	vaddpd	%%ymm13,%%ymm15,%%ymm15		\n\t"/* x0*y1 + fcy; ymm1 FREE */\
			"vmovaps	0x80(%%rax),%%ymm1		\n\t	vmovaps	0xa0(%%rax),%%ymm5			\n\t	vmovaps	0xc0(%%rax),%%ymm9			\n\t	vmovaps	0xe0(%%rax),%%ymm13			\n\t"/* Reload flo1 */\
			"vmulpd	    (%%rdx),%%ymm1,%%ymm1	\n\t	vmulpd	0x20(%%rdx),%%ymm5,%%ymm5	\n\t	vmulpd	0x40(%%rdx),%%ymm9,%%ymm9	\n\t	vmulpd	0x60(%%rdx),%%ymm13,%%ymm13	\n\t"/* x1*y0 */\
			"vaddpd	%%ymm3,%%ymm1,%%ymm1		\n\t	vaddpd	%%ymm7,%%ymm5,%%ymm5		\n\t	vaddpd	%%ymm11,%%ymm9,%%ymm9		\n\t	vaddpd	%%ymm15,%%ymm13,%%ymm13		\n\t"/* fprod1 = x0*y1 + x1*y0; ymm3 FREE */\
			"vroundpd $0,%%ymm1,%%ymm3			\n\t	vroundpd $0,%%ymm5,%%ymm7			\n\t	vroundpd $0,%%ymm9,%%ymm11			\n\t	vroundpd $0,%%ymm13,%%ymm15			\n\t"\
			"vsubpd	%%ymm3,%%ymm1,%%ymm1		\n\t	vsubpd	%%ymm7,%%ymm5,%%ymm5		\n\t	vsubpd	%%ymm11,%%ymm9,%%ymm9		\n\t	vsubpd	%%ymm15,%%ymm13,%%ymm13		\n\t"/* fprod1 -= fcy */\
			"vmulpd	-0x20(%%rbx),%%ymm1,%%ymm1	\n\t	vmulpd	-0x20(%%rbx),%%ymm5,%%ymm5	\n\t	vmulpd	-0x20(%%rbx),%%ymm9,%%ymm9	\n\t	vmulpd	-0x20(%%rbx),%%ymm13,%%ymm13\n\t"/* fprod0 *= two26f */\
			"vmulpd	     (%%rbx),%%ymm3,%%ymm3	\n\t	vmulpd	     (%%rbx),%%ymm7,%%ymm7	\n\t	vmulpd	     (%%rbx),%%ymm11,%%ymm11\n\t	vmulpd	     (%%rbx),%%ymm15,%%ymm15\n\t"/* fcy    *= two26i */\
			"vmovaps	%%ymm1,0x100(%%rax)		\n\t	vmovaps	%%ymm5,0x120(%%rax)			\n\t	vmovaps	%%ymm9,0x140(%%rax)			\n\t	vmovaps	%%ymm13,0x160(%%rax)		\n\t"/* Store fprod1 [in unused flo2 slot] to free up a register */\
			"vmovaps	0x80(%%rax),%%ymm1		\n\t	vmovaps	0xa0(%%rax),%%ymm5			\n\t	vmovaps	0xc0(%%rax),%%ymm9			\n\t	vmovaps	0xe0(%%rax),%%ymm13			\n\t"/* Reload flo1 */\
			"/* Digit 2: Require output to be nonnegative, so leave MSW unbalanced: */\n\t"\
			"vmulpd	    (%%rdx),%%ymm2,%%ymm2	\n\t	vmulpd	0x20(%%rdx),%%ymm6,%%ymm6	\n\t	vmulpd	0x40(%%rdx),%%ymm10,%%ymm10	\n\t	vmulpd	0x60(%%rdx),%%ymm14,%%ymm14	\n\t"/* x2*y0 */\
			"vmulpd	0x80(%%rdx),%%ymm1,%%ymm1	\n\t	vmulpd	0xa0(%%rdx),%%ymm5,%%ymm5	\n\t	vmulpd	0xc0(%%rdx),%%ymm9,%%ymm9	\n\t	vmulpd	0xe0(%%rdx),%%ymm13,%%ymm13	\n\t"/* x1*y1 */\
			"vmulpd	0x100(%%rdx),%%ymm0,%%ymm0	\n\t	vmulpd	0x120(%%rdx),%%ymm4,%%ymm4	\n\t	vmulpd	0x140(%%rdx),%%ymm8,%%ymm8	\n\t	vmulpd	0x160(%%rdx),%%ymm12,%%ymm12\n\t"/* x0*y2 */\
			"vaddpd	%%ymm1,%%ymm2,%%ymm2		\n\t	vaddpd	%%ymm5,%%ymm6,%%ymm6		\n\t	vaddpd	%%ymm9,%%ymm10,%%ymm10		\n\t	vaddpd	%%ymm13,%%ymm14,%%ymm14		\n\t"/* x1*y1 + x2*y0; ymm1 FREE */\
			"vaddpd	%%ymm0,%%ymm3,%%ymm3		\n\t	vaddpd	%%ymm4,%%ymm7,%%ymm7		\n\t	vaddpd	%%ymm8,%%ymm11,%%ymm11		\n\t	vaddpd	%%ymm12,%%ymm15,%%ymm15		\n\t"/* x0*y2 + fcy; ymm0 FREE */\
			"vaddpd	%%ymm3,%%ymm2,%%ymm2		\n\t	vaddpd	%%ymm7,%%ymm6,%%ymm6		\n\t	vaddpd	%%ymm11,%%ymm10,%%ymm10		\n\t	vaddpd	%%ymm15,%%ymm14,%%ymm14		\n\t"/* fprod2; ymm3 FREE */\
			"vroundpd $1,%%ymm2,%%ymm3			\n\t	vroundpd $1,%%ymm6,%%ymm7			\n\t	vroundpd $1,%%ymm10,%%ymm11			\n\t	vroundpd $1,%%ymm14,%%ymm15			\n\t"/* fcy = floor(fprod2) */\
			"vsubpd	%%ymm3,%%ymm2,%%ymm2		\n\t	vsubpd	%%ymm7,%%ymm6,%%ymm6		\n\t	vsubpd	%%ymm11,%%ymm10,%%ymm10		\n\t	vsubpd	%%ymm15,%%ymm14,%%ymm14		\n\t"/* fprod2 -= fcy */\
			"vmulpd	-0x20(%%rbx),%%ymm2,%%ymm2	\n\t	vmulpd	-0x20(%%rbx),%%ymm6,%%ymm6	\n\t	vmulpd	-0x20(%%rbx),%%ymm10,%%ymm10\n\t	vmulpd	-0x20(%%rbx),%%ymm14,%%ymm14\n\t"/* fprod2 *= two26f */\
			"vmovaps	    (%%rax),%%ymm0		\n\t	vmovaps	0x20(%%rax),%%ymm4			\n\t	vmovaps	0x40(%%rax),%%ymm8			\n\t	vmovaps	0x60(%%rax),%%ymm12			\n\t"/* Reload fprod0 */\
			"vmovaps	0x100(%%rax),%%ymm1		\n\t	vmovaps	0x120(%%rax),%%ymm5			\n\t	vmovaps	0x140(%%rax),%%ymm9			\n\t	vmovaps	0x160(%%rax),%%ymm13		\n\t"/* Reload fprod1 [in flo2 slot] */\
		/* MULH96(q,lo,lo) --> lo =(q*lo)/2^78: flo0,1,2 enter in ymm0,1,2. In the comments x = q, y = lo: */\
			"movq	%[__aq0],%%rdx		\n\t"\
			"/* Digit 0: */\n\t"\
			"vmovaps	%%ymm0,%%ymm3			\n\t	vmovaps	%%ymm4,%%ymm7				\n\t	vmovaps	%%ymm8,%%ymm11				\n\t	vmovaps	%%ymm12,%%ymm15				\n\t"/* ymm3 = cpy of y0 */\
			"vmulpd	    (%%rdx),%%ymm0,%%ymm0	\n\t	vmulpd	0x20(%%rdx),%%ymm4,%%ymm4	\n\t	vmulpd	0x40(%%rdx),%%ymm8,%%ymm8	\n\t	vmulpd	0x60(%%rdx),%%ymm12,%%ymm12	\n\t"/* fprod0 = y0*x0 */\
			"vmulpd	    (%%rbx),%%ymm0,%%ymm0	\n\t	vmulpd	    (%%rbx),%%ymm4,%%ymm4	\n\t	vmulpd	    (%%rbx),%%ymm8,%%ymm8	\n\t	vmulpd	    (%%rbx),%%ymm12,%%ymm12	\n\t"/* CY    *= two26i */\
			"/* Digit 1: */\n\t"\
			"vmulpd	0x80(%%rdx),%%ymm3,%%ymm3	\n\t	vmulpd	0xa0(%%rdx),%%ymm7,%%ymm7	\n\t	vmulpd	0xc0(%%rdx),%%ymm11,%%ymm11	\n\t	vmulpd	0xe0(%%rdx),%%ymm15,%%ymm15	\n\t"/* y0 *= x1 */\
			"vaddpd	%%ymm0,%%ymm3,%%ymm3		\n\t	vaddpd	%%ymm4,%%ymm7,%%ymm7		\n\t	vaddpd	%%ymm8,%%ymm11,%%ymm11		\n\t	vaddpd	%%ymm12,%%ymm15,%%ymm15		\n\t"/* y0*x1 + CY; ymm0 FREE */\
			"vmulpd	    (%%rdx),%%ymm1,%%ymm1	\n\t	vmulpd	0x20(%%rdx),%%ymm5,%%ymm5	\n\t	vmulpd	0x40(%%rdx),%%ymm9,%%ymm9	\n\t	vmulpd	0x60(%%rdx),%%ymm13,%%ymm13	\n\t"/* y1 *= x0 */\
			"vaddpd	%%ymm3,%%ymm1,%%ymm1		\n\t	vaddpd	%%ymm7,%%ymm5,%%ymm5		\n\t	vaddpd	%%ymm11,%%ymm9,%%ymm9		\n\t	vaddpd	%%ymm15,%%ymm13,%%ymm13		\n\t"/* fprod1 = x0*y1 + x1*y0; ymm3 FREE */\
			"vmulpd	    (%%rbx),%%ymm1,%%ymm1	\n\t	vmulpd	    (%%rbx),%%ymm5,%%ymm5	\n\t	vmulpd	    (%%rbx),%%ymm9,%%ymm9	\n\t	vmulpd	    (%%rbx),%%ymm13,%%ymm13	\n\t"/* CY    *= two26i */\
			"/* Digit 2: Require output to be nonnegative, so leave MSW unbalanced: */\n\t"\
			"vmovaps	%%ymm2,%%ymm3			\n\t	vmovaps	%%ymm6,%%ymm7				\n\t	vmovaps	%%ymm10,%%ymm11				\n\t	vmovaps	%%ymm14,%%ymm15				\n\t"/* ymm3 = cpy of y2 */\
			"vmovaps	0x100(%%rax),%%ymm0		\n\t	vmovaps	0x120(%%rax),%%ymm4			\n\t	vmovaps	0x140(%%rax),%%ymm8			\n\t	vmovaps	0x160(%%rax),%%ymm12		\n\t"/* Reload y1 [in flo2 slot] */\
			"vmulpd	    (%%rdx),%%ymm2,%%ymm2	\n\t	vmulpd	0x20(%%rdx),%%ymm6,%%ymm6	\n\t	vmulpd	0x40(%%rdx),%%ymm10,%%ymm10	\n\t	vmulpd	0x60(%%rdx),%%ymm14,%%ymm14	\n\t"/* y2 *= x0 */\
			"vmulpd	0x80(%%rdx),%%ymm0,%%ymm0	\n\t	vmulpd	0xa0(%%rdx),%%ymm4,%%ymm4	\n\t	vmulpd	0xc0(%%rdx),%%ymm8,%%ymm8	\n\t	vmulpd	0xe0(%%rdx),%%ymm12,%%ymm12	\n\t"/* y1 *= x1 */\
			"vaddpd	%%ymm1,%%ymm2,%%ymm2		\n\t	vaddpd	%%ymm5,%%ymm6,%%ymm6		\n\t	vaddpd	%%ymm9,%%ymm10,%%ymm10		\n\t	vaddpd	%%ymm13,%%ymm14,%%ymm14		\n\t"/* y2*x0 + CY; ymm1 FREE */\
			"vmovaps	    (%%rax),%%ymm1		\n\t	vmovaps	0x20(%%rax),%%ymm5			\n\t	vmovaps	0x40(%%rax),%%ymm9			\n\t	vmovaps	0x60(%%rax),%%ymm13			\n\t"/* Reload y0 */\
			"vmulpd	0x100(%%rdx),%%ymm1,%%ymm1	\n\t	vmulpd	0x120(%%rdx),%%ymm5,%%ymm5	\n\t	vmulpd	0x140(%%rdx),%%ymm9,%%ymm9	\n\t	vmulpd	0x160(%%rdx),%%ymm13,%%ymm13\n\t"/* y0 *= x2 */\
			"vaddpd	%%ymm0,%%ymm1,%%ymm1		\n\t	vaddpd	%%ymm4,%%ymm5,%%ymm5		\n\t	vaddpd	%%ymm8,%%ymm9,%%ymm9		\n\t	vaddpd	%%ymm12,%%ymm13,%%ymm13		\n\t"/* y1*x1 + y0*x2; ymm0 FREE */\
			"vaddpd	%%ymm1,%%ymm2,%%ymm2		\n\t	vaddpd	%%ymm5,%%ymm6,%%ymm6		\n\t	vaddpd	%%ymm9,%%ymm10,%%ymm10		\n\t	vaddpd	%%ymm13,%%ymm14,%%ymm14		\n\t"/* fprod2; ymm1 FREE */\
			"vroundpd $1,%%ymm2,%%ymm2			\n\t	vroundpd $1,%%ymm6,%%ymm6			\n\t	vroundpd $1,%%ymm10,%%ymm10			\n\t	vroundpd $1,%%ymm14,%%ymm14			\n\t"/* CY = floor(fprod2) */\
			"vmulpd	    (%%rbx),%%ymm2,%%ymm2	\n\t	vmulpd	    (%%rbx),%%ymm6,%%ymm6	\n\t	vmulpd	    (%%rbx),%%ymm10,%%ymm10	\n\t	vmulpd	    (%%rbx),%%ymm14,%%ymm14	\n\t"/* CY    *= two26i */\
		/* Precompute partial products needed for upper half: */\
			"vmovaps	%%ymm3,%%ymm1			\n\t	vmovaps	%%ymm7,%%ymm5				\n\t	vmovaps	%%ymm11,%%ymm9				\n\t	vmovaps	%%ymm15,%%ymm13				\n\t"/* ymm1 = cpy of y2 */\
			"vmovaps	0x100(%%rax),%%ymm0		\n\t	vmovaps	0x120(%%rax),%%ymm4			\n\t	vmovaps	0x140(%%rax),%%ymm8			\n\t	vmovaps	0x160(%%rax),%%ymm12		\n\t"/* Reload y1 [in flo2 slot] */\
			"vmulpd	0x100(%%rdx),%%ymm0,%%ymm0	\n\t	vmulpd	0x120(%%rdx),%%ymm4,%%ymm4	\n\t	vmulpd	0x140(%%rdx),%%ymm8,%%ymm8	\n\t	vmulpd	0x160(%%rdx),%%ymm12,%%ymm12\n\t"/* y1 *= x2 */\
			"vmulpd	0x80(%%rdx),%%ymm1,%%ymm1	\n\t	vmulpd	0xa0(%%rdx),%%ymm5,%%ymm5	\n\t	vmulpd	0xc0(%%rdx),%%ymm9,%%ymm9	\n\t	vmulpd	0xe0(%%rdx),%%ymm13,%%ymm13	\n\t"/* y2 *= x1 */\
			"vmulpd	0x100(%%rdx),%%ymm3,%%ymm3	\n\t	vmulpd	0x120(%%rdx),%%ymm7,%%ymm7	\n\t	vmulpd	0x140(%%rdx),%%ymm11,%%ymm11\n\t	vmulpd	0x160(%%rdx),%%ymm15,%%ymm15\n\t"/* y2 *= x2 */\
			"/* Digit 3: */\n\t"\
			"vaddpd	%%ymm2,%%ymm1,%%ymm1		\n\t	vaddpd	%%ymm6,%%ymm5,%%ymm5		\n\t	vaddpd	%%ymm10,%%ymm9,%%ymm9		\n\t	vaddpd	%%ymm14,%%ymm13,%%ymm13		\n\t"/* y2*x1 += CY */\
			"vaddpd	%%ymm1,%%ymm0,%%ymm0		\n\t	vaddpd	%%ymm5,%%ymm4,%%ymm4		\n\t	vaddpd	%%ymm9,%%ymm8,%%ymm8		\n\t	vaddpd	%%ymm13,%%ymm12,%%ymm12		\n\t"/* fprod3 = y1*x2 + y2*x1 + CY */\
			"vroundpd $0,%%ymm0,%%ymm2			\n\t	vroundpd $0,%%ymm4,%%ymm6			\n\t	vroundpd $0,%%ymm8,%%ymm10			\n\t	vroundpd $0,%%ymm12,%%ymm14			\n\t"\
			"vsubpd	%%ymm2,%%ymm0,%%ymm0		\n\t	vsubpd	%%ymm6,%%ymm4,%%ymm4		\n\t	vsubpd	%%ymm10,%%ymm8,%%ymm8		\n\t	vsubpd	%%ymm14,%%ymm12,%%ymm12		\n\t"/* fprod3 -= CY */\
			"vmulpd	-0x20(%%rbx),%%ymm0,%%ymm0	\n\t	vmulpd	-0x20(%%rbx),%%ymm4,%%ymm4	\n\t	vmulpd	-0x20(%%rbx),%%ymm8,%%ymm8	\n\t	vmulpd	-0x20(%%rbx),%%ymm12,%%ymm12\n\t"/* fprod3 *= two26f */\
			"vmulpd	     (%%rbx),%%ymm2,%%ymm2	\n\t	vmulpd	     (%%rbx),%%ymm6,%%ymm6	\n\t	vmulpd	     (%%rbx),%%ymm10,%%ymm10\n\t	vmulpd	     (%%rbx),%%ymm14,%%ymm14\n\t"/* CY    *= two26i */\
			"/* Digit 4: */\n\t"\
			"vaddpd	%%ymm2,%%ymm3,%%ymm3		\n\t	vaddpd	%%ymm6,%%ymm7,%%ymm7		\n\t	vaddpd	%%ymm10,%%ymm11,%%ymm11		\n\t	vaddpd	%%ymm14,%%ymm15,%%ymm15		\n\t"/* CY = y2*x2 + CY */\
			"vmovaps	%%ymm3,%%ymm1			\n\t	vmovaps	%%ymm7,%%ymm5				\n\t	vmovaps	%%ymm11,%%ymm9				\n\t	vmovaps	%%ymm15,%%ymm13				\n\t"/* fprod4 = cpy of CY */\
			"vroundpd $0,%%ymm3,%%ymm3			\n\t	vroundpd $0,%%ymm7,%%ymm7			\n\t	vroundpd $0,%%ymm11,%%ymm11			\n\t	vroundpd $0,%%ymm15,%%ymm15			\n\t"\
			"vsubpd	%%ymm3,%%ymm1,%%ymm1		\n\t	vsubpd	%%ymm7,%%ymm5,%%ymm5		\n\t	vsubpd	%%ymm11,%%ymm9,%%ymm9		\n\t	vsubpd	%%ymm15,%%ymm13,%%ymm13		\n\t"/* fprod4 -= CY. Note we use the carry in the next section! */\
			"vmulpd	-0x20(%%rbx),%%ymm1,%%ymm1	\n\t	vmulpd	-0x20(%%rbx),%%ymm5,%%ymm5	\n\t	vmulpd	-0x20(%%rbx),%%ymm9,%%ymm9	\n\t	vmulpd	-0x20(%%rbx),%%ymm13,%%ymm13\n\t"/* fprod4 *= two26f */\
		/* If h < l, calculate h-l+q; otherwise h-l. Use leading 52 bits to approximate the full 78-bit compare. Result is in [0, q)	. */\
			"vmovaps	0x200(%%rax),%%ymm2		\n\t	vmovaps	0x220(%%rax),%%ymm6			\n\t	vmovaps	0x240(%%rax),%%ymm10		\n\t	vmovaps	0x260(%%rax),%%ymm14		\n\t"/* fhi, top 52 bits */\
			"vmulpd	-0x20(%%rbx),%%ymm3,%%ymm3	\n\t	vmulpd	-0x20(%%rbx),%%ymm7,%%ymm7	\n\t	vmulpd	-0x20(%%rbx),%%ymm11,%%ymm11\n\t	vmulpd	-0x20(%%rbx),%%ymm15,%%ymm15\n\t"/* flo2 *= two26f */\
			"vaddpd	%%ymm1,%%ymm3,%%ymm3		\n\t	vaddpd	%%ymm5,%%ymm7,%%ymm7		\n\t	vaddpd	%%ymm9,%%ymm11,%%ymm11		\n\t	vaddpd	%%ymm13,%%ymm15,%%ymm15		\n\t"/* flo, top 52 bits */\
			"vmovaps	0x180(%%rax),%%ymm1		\n\t	vmovaps	0x1a0(%%rax),%%ymm5			\n\t	vmovaps	0x1c0(%%rax),%%ymm9			\n\t	vmovaps	0x1e0(%%rax),%%ymm13		\n\t"/* fhi, low 26 bits */\
			"vcmppd	$0x1,%%ymm3,%%ymm2,%%ymm2	\n\t	vcmppd	$0x1,%%ymm7,%%ymm6,%%ymm6	\n\t	vcmppd	$0x1,%%ymm11,%%ymm10,%%ymm10\n\t	vcmppd	$0x1,%%ymm15,%%ymm14,%%ymm14\n\t"/* bitmask = (fhi < flo)? */\
			"vmovaps	%%ymm2,    (%%rax)		\n\t	vmovaps	%%ymm6,0x20(%%rax)			\n\t	vmovaps	%%ymm10,0x40(%%rax)			\n\t	vmovaps	%%ymm14,0x60(%%rax)			\n\t"/* Store bitmask to free up a register */\
			"vmovaps	0x200(%%rax),%%ymm2		\n\t	vmovaps	0x220(%%rax),%%ymm6			\n\t	vmovaps	0x240(%%rax),%%ymm10		\n\t	vmovaps	0x260(%%rax),%%ymm14		\n\t"/* fhi, top 52 bits */\
			"vsubpd	%%ymm0,%%ymm1,%%ymm1		\n\t	vsubpd	%%ymm4,%%ymm5,%%ymm5		\n\t	vsubpd	%%ymm8,%%ymm9,%%ymm9		\n\t	vsubpd	%%ymm12,%%ymm13,%%ymm13		\n\t"/* (fhi - flo), low 26 bits, ymm0 FREE */\
			"vsubpd	%%ymm3,%%ymm2,%%ymm2		\n\t	vsubpd	%%ymm7,%%ymm6,%%ymm6		\n\t	vsubpd	%%ymm11,%%ymm10,%%ymm10		\n\t	vsubpd	%%ymm15,%%ymm14,%%ymm14		\n\t"/* (fhi - flo), top 52 bits, ymm3 FREE */\
			"vmovaps	-0x20(%%rbx),%%ymm3		\n\t	vmovaps	-0x20(%%rbx),%%ymm7			\n\t	vmovaps	-0x20(%%rbx),%%ymm11		\n\t	vmovaps	-0x20(%%rbx),%%ymm15		\n\t"/* 2^26 */\
			"vmulpd	    (%%rdx),%%ymm3,%%ymm3	\n\t	vmulpd	0x20(%%rdx),%%ymm7,%%ymm7	\n\t	vmulpd	0x40(%%rdx),%%ymm11,%%ymm11	\n\t	vmulpd	0x60(%%rdx),%%ymm15,%%ymm15	\n\t"/* fq, low 26 bits */\
			"vmovaps	0x180(%%rdx),%%ymm0		\n\t	vmovaps	0x1a0(%%rdx),%%ymm4			\n\t	vmovaps	0x1c0(%%rdx),%%ymm8			\n\t	vmovaps	0x1e0(%%rdx),%%ymm12		\n\t"/* qhi52 :+ fq, top 52 bits */\
			"vandpd	    (%%rax),%%ymm0,%%ymm0	\n\t	vandpd	0x20(%%rax),%%ymm4,%%ymm4	\n\t	vandpd	0x40(%%rax),%%ymm8,%%ymm8	\n\t	vandpd	0x60(%%rax),%%ymm12,%%ymm12	\n\t"/* qhi52 & bitmask */\
			"vaddpd	%%ymm0,%%ymm2,%%ymm2		\n\t	vaddpd	%%ymm4,%%ymm6,%%ymm6		\n\t	vaddpd	%%ymm8,%%ymm10,%%ymm10		\n\t	vaddpd	%%ymm12,%%ymm14,%%ymm14		\n\t"/* xhi = (h-l)hi + (qhi52 & bitmask) */\
			"vmovaps	%%ymm3,0x80(%%rax)		\n\t	vmovaps	%%ymm7,0xa0(%%rax)			\n\t	vmovaps	%%ymm11,0xc0(%%rax)			\n\t	vmovaps	%%ymm15,0xe0(%%rax)			\n\t"/* Store qlo26 to free up a register */\
			"vandpd	    (%%rax),%%ymm3,%%ymm3	\n\t	vandpd	0x20(%%rax),%%ymm7,%%ymm7	\n\t	vandpd	0x40(%%rax),%%ymm11,%%ymm11	\n\t	vandpd	0x60(%%rax),%%ymm15,%%ymm15	\n\t"/* qlo26 &= bitmask */\
			"vaddpd	%%ymm3,%%ymm1,%%ymm1		\n\t	vaddpd	%%ymm7,%%ymm5,%%ymm5		\n\t	vaddpd	%%ymm11,%%ymm9,%%ymm9		\n\t	vaddpd	%%ymm15,%%ymm13,%%ymm13		\n\t"/* xlo = (h-l)lo + (qlo26 & bitmask) */\
			"/* xhi52,xlo26 in ymm2,ymm1; qhi52,qlo26 in 0x180(%%rdx),0x80(%%rax): */\n\t"\
		/* if((pshift >> j) & (uint64)1) { */\
			"movslq	%[__pshift],%%rax		\n\t"\
			"movslq	%[__j],%%rcx			\n\t"\
			"shrq	%%cl,%%rax				\n\t"\
			"andq	$0x1,%%rax				\n\t"\
		"je twopmodq78_3wdq16_gcc64	\n\t"\
			"movq	%[__ax0],%%rax			\n\t"\
			"vaddpd	%%ymm2,%%ymm2,%%ymm2		\n\t	vaddpd	%%ymm6,%%ymm6,%%ymm6		\n\t	vaddpd	%%ymm10,%%ymm10,%%ymm10		\n\t	vaddpd	%%ymm14,%%ymm14,%%ymm14		\n\t"/* top 52 bits */\
			"vaddpd	%%ymm1,%%ymm1,%%ymm1		\n\t	vaddpd	%%ymm5,%%ymm5,%%ymm5		\n\t	vaddpd	%%ymm9,%%ymm9,%%ymm9		\n\t	vaddpd	%%ymm13,%%ymm13,%%ymm13		\n\t"/* low 26 bits */\
			"/* If x > q, subtract q: */\n\t"\
			"vmovaps	0x180(%%rdx),%%ymm0		\n\t	vmovaps	0x1a0(%%rdx),%%ymm4			\n\t	vmovaps	0x1c0(%%rdx),%%ymm8			\n\t	vmovaps	0x1e0(%%rdx),%%ymm12		\n\t"/* qhi52 */\
			"vmovaps	%%ymm0,%%ymm3			\n\t	vmovaps	%%ymm4,%%ymm7				\n\t	vmovaps	%%ymm8,%%ymm11				\n\t	vmovaps	%%ymm12,%%ymm15				\n\t"/* cpy of qhi */\
			"vcmppd	$0x2,%%ymm2,%%ymm3,%%ymm3	\n\t	vcmppd	$0x2,%%ymm6,%%ymm7,%%ymm7	\n\t	vcmppd	$0x2,%%ymm10,%%ymm11,%%ymm11\n\t	vcmppd	$0x2,%%ymm14,%%ymm15,%%ymm15\n\t"/* bitmask = (qhi <= xhi)? */\
			"vandpd	%%ymm3,%%ymm0,%%ymm0		\n\t	vandpd	%%ymm7,%%ymm4,%%ymm4		\n\t	vandpd	%%ymm11,%%ymm8,%%ymm8		\n\t	vandpd	%%ymm15,%%ymm12,%%ymm12		\n\t"/* qhi52 & bitmask */\
			"vandpd	0x80(%%rax),%%ymm3,%%ymm3	\n\t	vandpd	0xa0(%%rax),%%ymm7,%%ymm7	\n\t	vandpd	0xc0(%%rax),%%ymm11,%%ymm11	\n\t	vandpd	0xe0(%%rax),%%ymm15,%%ymm15	\n\t"/* qlo26 & bitmask */\
			"vsubpd	%%ymm0,%%ymm2,%%ymm2		\n\t	vsubpd	%%ymm4,%%ymm6,%%ymm6		\n\t	vsubpd	%%ymm8,%%ymm10,%%ymm10		\n\t	vsubpd	%%ymm12,%%ymm14,%%ymm14		\n\t"/* x % q, top 52 bits */\
			"vsubpd	%%ymm3,%%ymm1,%%ymm1		\n\t	vsubpd	%%ymm7,%%ymm5,%%ymm5		\n\t	vsubpd	%%ymm11,%%ymm9,%%ymm9		\n\t	vsubpd	%%ymm15,%%ymm13,%%ymm13		\n\t"/* x % q, low 26 bits */\
		"twopmodq78_3wdq16_gcc64:	\n\t"\
		/* } */\
		/* Normalize the result: */\
			"vmovaps	     (%%rbx),%%ymm3		\n\t	vmovaps	-0x20(%%rbx),%%ymm7		\n\t"/* two26i,f */\
			"vmulpd	%%ymm3,%%ymm1,%%ymm1		\n\t	vmulpd	%%ymm3,%%ymm5,%%ymm5		\n\t	vmulpd	%%ymm3,%%ymm9,%%ymm9		\n\t	vmulpd	%%ymm3,%%ymm13,%%ymm13		\n\t"/* xlo *= two26i */\
			"vmulpd	%%ymm3,%%ymm2,%%ymm2		\n\t	vmulpd	%%ymm3,%%ymm6,%%ymm6		\n\t	vmulpd	%%ymm3,%%ymm10,%%ymm10		\n\t	vmulpd	%%ymm3,%%ymm14,%%ymm14		\n\t"/* xhi *= two26i */\
			"vmovaps	%%ymm1,%%ymm0			\n\t	vmovaps	%%ymm5,%%ymm4				\n\t	vmovaps	%%ymm9,%%ymm8				\n\t	vmovaps	%%ymm13,%%ymm12				\n\t"/* Init: fcy = cpy of ~xlo */\
			"vroundpd $0,%%ymm1,%%ymm1			\n\t	vroundpd $0,%%ymm5,%%ymm5			\n\t	vroundpd $0,%%ymm9,%%ymm9			\n\t	vroundpd $0,%%ymm13,%%ymm13			\n\t"/* fcy */\
			"vsubpd	%%ymm1,%%ymm0,%%ymm0		\n\t	vsubpd	%%ymm5,%%ymm4,%%ymm4		\n\t	vsubpd	%%ymm9,%%ymm8,%%ymm8		\n\t	vsubpd	%%ymm13,%%ymm12,%%ymm12		\n\t"/* fx0 -= fcy */\
			"vmulpd	%%ymm7,%%ymm0,%%ymm0		\n\t	vmulpd	%%ymm7,%%ymm4,%%ymm4		\n\t	vmulpd	%%ymm7,%%ymm8,%%ymm8		\n\t	vmulpd	%%ymm7,%%ymm12,%%ymm12		\n\t"/* fx0 *= two26f */\
			"vmulpd	%%ymm3,%%ymm1,%%ymm1		\n\t	vmulpd	%%ymm3,%%ymm5,%%ymm5		\n\t	vmulpd	%%ymm3,%%ymm9,%%ymm9		\n\t	vmulpd	%%ymm3,%%ymm13,%%ymm13		\n\t"/* fcy *= two26i */\
			"vaddpd	%%ymm1,%%ymm2,%%ymm2		\n\t	vaddpd	%%ymm5,%%ymm6,%%ymm6		\n\t	vaddpd	%%ymm9,%%ymm10,%%ymm10		\n\t	vaddpd	%%ymm13,%%ymm14,%%ymm14		\n\t"/* Add carry into xhi */\
			"vmovaps	%%ymm2,%%ymm1			\n\t	vmovaps	%%ymm6,%%ymm5				\n\t	vmovaps	%%ymm10,%%ymm9				\n\t	vmovaps	%%ymm14,%%ymm13				\n\t"/* Init: fcy = cpy of ~xhi */\
			"vroundpd $0,%%ymm2,%%ymm2			\n\t	vroundpd $0,%%ymm6,%%ymm6			\n\t	vroundpd $0,%%ymm10,%%ymm10			\n\t	vroundpd $0,%%ymm14,%%ymm14			\n\t"/* fx2 = fcy (no divide by 2^26) */\
			"vsubpd	%%ymm2,%%ymm1,%%ymm1		\n\t	vsubpd	%%ymm6,%%ymm5,%%ymm5		\n\t	vsubpd	%%ymm10,%%ymm9,%%ymm9		\n\t	vsubpd	%%ymm14,%%ymm13,%%ymm13		\n\t"/* fx1 -= fcy */\
			"vmulpd	%%ymm7,%%ymm1,%%ymm1		\n\t	vmulpd	%%ymm7,%%ymm5,%%ymm5		\n\t	vmulpd	%%ymm7,%%ymm9,%%ymm9		\n\t	vmulpd	%%ymm7,%%ymm13,%%ymm13		\n\t"/* fx1 *= two26f */\
			:					/* outputs: none */\
			: [__aq0]	 "m" (Xaq0)	/* All inputs from memory addresses here */\
			 ,[__aqinv0] "m" (Xaqinv0)	\
			 ,[__ax0]	 "m" (Xax0)		\
			 ,[__two26i] "m" (Xtwo26i)	\
			 ,[__pshift] "m" (Xpshift)	\
			 ,[__j]		 "m" (Xj)		\
			: "cc","memory","cl","rax","rbx","rcx","rdx","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
		);\
		}


		/* Address-offset map: _q4 [imagined AVX impl] compared to _q8:
			OLD:		NEW:
	SSE	AVX	fq0			aq0
	20	40	fq1   		aq1
	40	80	fq2   		aq2
	60	c0	fqhi52		aqhi52
	80	100	fqinv0		aqinv0
	a0	140	fqinv1		aqinv1
	c0	180	fqinv2		aqinv2
	e0	1c0	fx0   		ax0
	100	200	fx1   		ax1
	120	240	fx2   		ax2
	140	280	flo0  		[allocated but unnamed]	<*** These higher address offsets need modification:
	160	2c0	flo1  		[allocated but unnamed]
	180	300	flo2  		n/a
	1a0	340	fhi0  	Use ax0[1c0]
	1c0	380	fhi1  	Use ax1[200]
	1e0	3c0	fhi2  	Use ax2[240]
	200	400	two13i		300
	220	440	two26f		340
	240	480	two26i		380
	260	4c0	sse2_rnd	n/a
	280	500	half		n/a

		OLD:
			"vmovaps	%%ymm3,0x340(%%rdx)				\n\t	vmovaps	%%ymm11,0x360(%%rdx)		\n\t"// Spill hi0 into ax0 //\
			"vmovaps	%%ymm4,0x380(%%rdx)				\n\t	vmovaps	%%ymm12,0x3a0(%%rdx)		\n\t"// Spill hi1 into ax1 //\
		NEW:
			"vmovaps	%%ymm3,0x1c0(%%rdx)				\n\t	vmovaps	%%ymm11,0x1e0(%%rdx)		\n\t"// Spill hi0 into ax0 //\
			"vmovaps	%%ymm4,0x200(%%rdx)				\n\t	vmovaps	%%ymm12,0x220(%%rdx)		\n\t"// Spill hi1 into ax1 //\
		*/

		#define SSE2_twopmodq78_modmul_q8(Xfq0,Xpshift,Xj)\
		{\
		__asm__ volatile (\
			"movq	%[__fq0],%%rdx		\n\t"\
			"/* SQR_LOHI78_3WORD_DOUBLE_q4(): */\n\t"\
			"vmovaps	0x300(%%rdx),%%ymm6	/* two13i */\n\t"\
			"/* fx0,1,2 assumed in ymm0,2,4 on loop entry */\n\t"\
			"vmulpd		%%ymm6,%%ymm0,%%ymm0			\n\t	vmulpd	%%ymm6 ,%%ymm8	,%%ymm8		\n\t"\
			"vmulpd		%%ymm6,%%ymm2,%%ymm2			\n\t	vmulpd	%%ymm6 ,%%ymm10,%%ymm10		\n\t"\
			"vmulpd		%%ymm6,%%ymm4,%%ymm4			\n\t	vmulpd	%%ymm6 ,%%ymm12,%%ymm12		\n\t"\
			"vmovaps	%%ymm0,%%ymm1					\n\t	vmovaps	%%ymm8 ,%%ymm9 				\n\t"\
			"vmovaps	%%ymm2,%%ymm3					\n\t	vmovaps	%%ymm10,%%ymm11				\n\t"\
			"vaddpd		%%ymm1,%%ymm1,%%ymm1			\n\t	vaddpd	%%ymm9 ,%%ymm9 ,%%ymm9 		\n\t"\
			"vaddpd		%%ymm3,%%ymm3,%%ymm3			\n\t	vaddpd	%%ymm11,%%ymm11,%%ymm11		\n\t"\
			"vmovaps	%%ymm1,%%ymm5					\n\t	vmovaps	%%ymm9 ,%%ymm13				\n\t"\
			"vmulpd		%%ymm0,%%ymm0,%%ymm0			\n\t	vmulpd	%%ymm8 ,%%ymm8 ,%%ymm8 		\n\t"\
			"vmulpd		%%ymm2,%%ymm1,%%ymm1			\n\t	vmulpd	%%ymm10,%%ymm9 ,%%ymm9 		\n\t"\
			"vmulpd		%%ymm2,%%ymm2,%%ymm2			\n\t	vmulpd	%%ymm10,%%ymm10,%%ymm10		\n\t"\
			"vmulpd		%%ymm4,%%ymm5,%%ymm5			\n\t	vmulpd	%%ymm12,%%ymm13,%%ymm13		\n\t"\
			"vmulpd		%%ymm4,%%ymm3,%%ymm3			\n\t	vmulpd	%%ymm12,%%ymm11,%%ymm11		\n\t"\
			"vmulpd		%%ymm4,%%ymm4,%%ymm4			\n\t	vmulpd	%%ymm12,%%ymm12,%%ymm12		\n\t"\
			"vaddpd		%%ymm5,%%ymm2,%%ymm2			\n\t	vaddpd	%%ymm13,%%ymm10,%%ymm10		/* Move this part of Digit 2 computation here to free up ymm5,13: */\n\t"\
			"/* Digit 0: */\n\t"\
			"vmovaps	0x340(%%rdx),%%ymm15			\n\t	vmovaps	0x320(%%rdx),%%ymm5			/* ymm15,5 = two26i,f */\n\t"\
			"vroundpd	$0,%%ymm0,%%ymm6				\n\t	vroundpd	$0,%%ymm8 ,%%ymm14		\n\t"\
			"vsubpd		%%ymm6 ,%%ymm0,%%ymm0			\n\t	vsubpd	%%ymm14,%%ymm8 ,%%ymm8 		\n\t"\
			"vmulpd		%%ymm5 ,%%ymm0,%%ymm0			\n\t	vmulpd	%%ymm5 ,%%ymm8 ,%%ymm8 		\n\t"\
			"vmulpd		%%ymm15,%%ymm6,%%ymm6			\n\t	vmulpd	%%ymm15,%%ymm14,%%ymm14		\n\t"\
			"/* Digit 1: */\n\t"\
			"vaddpd		%%ymm6,%%ymm1,%%ymm1			\n\t	vaddpd	%%ymm14,%%ymm9,%%ymm9		\n\t"\
			"vroundpd	$0,%%ymm1,%%ymm6				\n\t	vroundpd $0,%%ymm9 ,%%ymm14			\n\t"\
			"vsubpd		%%ymm6,%%ymm1,%%ymm1			\n\t	vsubpd	%%ymm14,%%ymm9 ,%%ymm9 		\n\t"\
			"vmulpd		%%ymm5,%%ymm1,%%ymm1			\n\t	vmulpd	%%ymm5 ,%%ymm9 ,%%ymm9 		\n\t"\
			"vmulpd		%%ymm15,%%ymm6,%%ymm6			\n\t	vmulpd	%%ymm15,%%ymm14,%%ymm14		\n\t"\
			"/* Digit 2: Require both hi and lo half of output to be nonnegative, so leave unbalanced: */\n\t"\
			"vaddpd		%%ymm6,%%ymm2,%%ymm2			\n\t	vaddpd	%%ymm14,%%ymm10,%%ymm10		\n\t"\
			"vroundpd	$1,%%ymm2,%%ymm6				\n\t	vroundpd	$1,%%ymm10,%%ymm14		\n\t"\
			"vsubpd		%%ymm6,%%ymm2,%%ymm2			\n\t	vsubpd	%%ymm14,%%ymm10,%%ymm10		\n\t"\
			"vmulpd		%%ymm5,%%ymm2,%%ymm2			\n\t	vmulpd	%%ymm5 ,%%ymm10,%%ymm10		\n\t"\
			"vmulpd		%%ymm15,%%ymm6,%%ymm6			\n\t	vmulpd	%%ymm15,%%ymm14,%%ymm14		\n\t"\
			"/* Digit 3: */\n\t"\
			"vaddpd		%%ymm6,%%ymm3,%%ymm3			\n\t	vaddpd	%%ymm14,%%ymm11,%%ymm11		\n\t"\
			"vroundpd	$0,%%ymm3,%%ymm6				\n\t	vroundpd	$0,%%ymm11,%%ymm14		\n\t"\
			"vsubpd		%%ymm6,%%ymm3,%%ymm3			\n\t	vsubpd	%%ymm14,%%ymm11,%%ymm11		\n\t"\
			"vmulpd		%%ymm5,%%ymm3,%%ymm3			\n\t	vmulpd	%%ymm5 ,%%ymm11,%%ymm11		\n\t"\
			"vmulpd		%%ymm15,%%ymm6,%%ymm6			\n\t	vmulpd	%%ymm15,%%ymm14,%%ymm14		\n\t"\
			"/* Digit 4: */\n\t"\
			"vaddpd		%%ymm6,%%ymm4,%%ymm4			\n\t	vaddpd	%%ymm14,%%ymm12,%%ymm12		\n\t"\
			"vroundpd	$0,%%ymm4,%%ymm6				\n\t	vroundpd	$0,%%ymm12,%%ymm14		\n\t"\
			"vsubpd		%%ymm6,%%ymm4,%%ymm4			\n\t	vsubpd	%%ymm14,%%ymm12,%%ymm12		\n\t"\
			"vmulpd		%%ymm5,%%ymm4,%%ymm4			\n\t	vmulpd	%%ymm5 ,%%ymm12,%%ymm12		\n\t"\
			"/* Digit 5 = the carry. flo0,1,2 in ymm0,1,2; fhi0,1,2 in ymm3,4,6 */\n\t"\
			"vmulpd		%%ymm5,%%ymm6,%%ymm6			\n\t	vmulpd	%%ymm5 ,%%ymm14,%%ymm14		\n\t"\
			"vaddpd		%%ymm6,%%ymm4,%%ymm4			\n\t	vaddpd	%%ymm14,%%ymm12,%%ymm12		\n\t"\
			"vmovaps	%%ymm3,0x1c0(%%rdx)				\n\t	vmovaps	%%ymm11,0x1e0(%%rdx)		\n\t"/* Spill hi0 into ax0 */\
			"vmovaps	%%ymm4,0x200(%%rdx)				\n\t	vmovaps	%%ymm12,0x220(%%rdx)		\n\t"/* Spill hi1 into ax1 */\
			"/* flo = MULL78_3WORD_DOUBLE_q2(flo, fqinv): */\n\t"\
			"/* Digit 0: */\n\t"\
			"vmovaps	%%ymm0,%%ymm3					\n\t	vmovaps	%%ymm8 ,%%ymm11				\n\t"\
			"vmovaps	%%ymm0,%%ymm4					\n\t	vmovaps	%%ymm8 ,%%ymm12				\n\t"\
			"vmulpd	  0x100(%%rdx),%%ymm0,%%ymm0		\n\t	vmulpd 0x120(%%rdx),%%ymm8 ,%%ymm8 	\n\t"\
			"vroundpd	$0,%%ymm0,%%ymm6				\n\t	vroundpd	$0,%%ymm8 ,%%ymm14		\n\t"\
			"vsubpd		%%ymm6,%%ymm0,%%ymm0			\n\t	vsubpd	%%ymm14,%%ymm8 ,%%ymm8 		\n\t"\
			"vmulpd		%%ymm5,%%ymm0,%%ymm0			\n\t	vmulpd	%%ymm5 ,%%ymm8 ,%%ymm8 		\n\t"\
			"vmulpd		%%ymm15,%%ymm6,%%ymm6			\n\t	vmulpd	%%ymm15,%%ymm14,%%ymm14		\n\t"\
			"/* Digit 1: */\n\t"\
			"vmulpd	  0x140(%%rdx),%%ymm3,%%ymm3		\n\t	vmulpd 0x160(%%rdx),%%ymm11,%%ymm11	\n\t"\
			"vaddpd		%%ymm6,%%ymm3,%%ymm3			\n\t	vaddpd	%%ymm14,%%ymm11,%%ymm11		\n\t"\
			"vmovaps	%%ymm1,%%ymm6					\n\t	vmovaps	%%ymm9 ,%%ymm14				\n\t"\
			"vmulpd	  0x100(%%rdx),%%ymm1,%%ymm1		\n\t	vmulpd 0x120(%%rdx),%%ymm9 ,%%ymm9 	\n\t"\
			"vaddpd		%%ymm3,%%ymm1,%%ymm1			\n\t	vaddpd	%%ymm11,%%ymm9 ,%%ymm9 		\n\t"\
			"vroundpd	$0,%%ymm1,%%ymm3				\n\t	vroundpd	$0,%%ymm9 ,%%ymm11		\n\t"\
			"vsubpd		%%ymm3,%%ymm1,%%ymm1			\n\t	vsubpd	%%ymm11,%%ymm9 ,%%ymm9 		\n\t"\
			"vmulpd		%%ymm5,%%ymm1,%%ymm1			\n\t	vmulpd	%%ymm5 ,%%ymm9 ,%%ymm9 		\n\t"\
			"vmulpd		%%ymm15,%%ymm3,%%ymm3			\n\t	vmulpd	%%ymm15,%%ymm11,%%ymm11		\n\t"\
			"/* Digit 2: Require output to be nonnegative, so leave MSW unbalanced: */\n\t"\
			"vmulpd	  0x100(%%rdx),%%ymm2,%%ymm2		\n\t	vmulpd 0x120(%%rdx),%%ymm10,%%ymm10	\n\t"\
			"vmulpd	  0x140(%%rdx),%%ymm6,%%ymm6		\n\t	vmulpd 0x160(%%rdx),%%ymm14,%%ymm14	\n\t"\
			"vmulpd	  0x180(%%rdx),%%ymm4,%%ymm4		\n\t	vmulpd 0x1a0(%%rdx),%%ymm12,%%ymm12	\n\t"\
			"vaddpd		%%ymm6,%%ymm2,%%ymm2			\n\t	vaddpd	%%ymm14,%%ymm10,%%ymm10		\n\t"\
			"vaddpd		%%ymm4,%%ymm3,%%ymm3			\n\t	vaddpd	%%ymm12,%%ymm11,%%ymm11		\n\t"\
			"vaddpd		%%ymm3,%%ymm2,%%ymm2			\n\t	vaddpd	%%ymm11,%%ymm10,%%ymm10		\n\t"\
			"vroundpd	$1,%%ymm2,%%ymm3				\n\t	vroundpd	$1,%%ymm10,%%ymm11		\n\t"\
			"vsubpd		%%ymm3,%%ymm2,%%ymm2			\n\t	vsubpd	%%ymm11,%%ymm10,%%ymm10		\n\t"\
			"vmulpd		%%ymm5,%%ymm2,%%ymm2			\n\t	vmulpd	%%ymm5 ,%%ymm10,%%ymm10		\n\t"\
			"/* MULH96(q,lo,lo) --> lo = (q*lo)/2^78 */\n\t"\
			"/*** movq	%[__fq0],%%rdx		Moved to start of asm-block ***/\n\t"\
			"/* Digit 0: */\n\t"\
			"vmovaps	%%ymm0,%%ymm3					\n\t	vmovaps	%%ymm8 ,%%ymm11				\n\t"\
			"vmovaps	%%ymm0,%%ymm4					\n\t	vmovaps	%%ymm8 ,%%ymm12				\n\t"\
			"vmulpd	   (%%rdx),%%ymm0,%%ymm0			\n\t	vmulpd	0x20(%%rdx),%%ymm8 ,%%ymm8 	\n\t"\
			"vroundpd	$0,%%ymm0,%%ymm0				\n\t	vroundpd	$0,%%ymm8 ,%%ymm8		\n\t"\
			"vmulpd		%%ymm15,%%ymm0,%%ymm0			\n\t	vmulpd	%%ymm15,%%ymm8 ,%%ymm8 		\n\t"\
			"/* Digit 1: */\n\t"\
			"vmovaps	%%ymm1,%%ymm6					\n\t	vmovaps	%%ymm9 ,%%ymm14				\n\t"\
			"vmulpd	   0x40(%%rdx),%%ymm3,%%ymm3		\n\t	vmulpd	0x60(%%rdx),%%ymm11,%%ymm11	\n\t"\
			"vaddpd		%%ymm0,%%ymm3,%%ymm3			\n\t	vaddpd	%%ymm8 ,%%ymm11,%%ymm11		\n\t"\
			"vmulpd	   (%%rdx),%%ymm1,%%ymm1			\n\t	vmulpd	0x20(%%rdx),%%ymm9 ,%%ymm9 	\n\t"\
			"vaddpd		%%ymm3,%%ymm1,%%ymm1			\n\t	vaddpd	%%ymm11,%%ymm9 ,%%ymm9 		\n\t"\
			"vroundpd	$0,%%ymm1,%%ymm1				\n\t	vroundpd	$0,%%ymm9 ,%%ymm9		\n\t"\
			"vmulpd		%%ymm15,%%ymm1,%%ymm1			\n\t	vmulpd	%%ymm15,%%ymm9 ,%%ymm9 		\n\t"\
			"/* Digit 2: Require output to be nonnegative, so leave MSW unbalanced: */\n\t"\
			"vmovaps	%%ymm2,%%ymm3					\n\t	vmovaps	%%ymm10,%%ymm11				\n\t"\
			"vmovaps	%%ymm6,%%ymm0					\n\t	vmovaps	%%ymm14,%%ymm8				\n\t"\
			"vmulpd	   (%%rdx),%%ymm2,%%ymm2			\n\t	vmulpd	0x20(%%rdx),%%ymm10,%%ymm10	\n\t"\
			"vmulpd	   0x40(%%rdx),%%ymm6,%%ymm6		\n\t	vmulpd	0x60(%%rdx),%%ymm14,%%ymm14	\n\t"\
			"vmulpd	   0x80(%%rdx),%%ymm4,%%ymm4		\n\t	vmulpd	0xa0(%%rdx),%%ymm12,%%ymm12	\n\t"\
			"vaddpd		%%ymm6,%%ymm2,%%ymm2			\n\t	vaddpd	%%ymm14,%%ymm10,%%ymm10		\n\t"\
			"vaddpd		%%ymm4,%%ymm1,%%ymm1			\n\t	vaddpd	%%ymm12,%%ymm9 ,%%ymm9 		\n\t"\
			"vaddpd		%%ymm1,%%ymm2,%%ymm2			\n\t	vaddpd	%%ymm9 ,%%ymm10,%%ymm10		\n\t"\
			"vroundpd	$1,%%ymm2,%%ymm2				\n\t	vroundpd	$1,%%ymm10,%%ymm10		\n\t"\
			"vmulpd		%%ymm15,%%ymm2,%%ymm2			\n\t	vmulpd	%%ymm15,%%ymm10,%%ymm10		\n\t"\
			"/* Precompute all the needed partial products: */\n\t"\
			"vmovaps	%%ymm3,%%ymm1					\n\t	vmovaps	%%ymm11,%%ymm9				\n\t"\
			"vmulpd	   0x80(%%rdx),%%ymm0,%%ymm0		\n\t	vmulpd	0xa0(%%rdx),%%ymm8 ,%%ymm8 	\n\t"\
			"vmulpd	   0x40(%%rdx),%%ymm1,%%ymm1		\n\t	vmulpd	0x60(%%rdx),%%ymm9 ,%%ymm9 	\n\t"\
			"vmulpd	   0x80(%%rdx),%%ymm3,%%ymm3		\n\t	vmulpd	0xa0(%%rdx),%%ymm11,%%ymm11	\n\t"\
			"/* Digit 3: */\n\t"\
			"vaddpd		%%ymm2,%%ymm1,%%ymm1			\n\t	vaddpd	%%ymm10,%%ymm9 ,%%ymm9 		\n\t"\
			"vaddpd		%%ymm1,%%ymm0,%%ymm0			\n\t	vaddpd	%%ymm9 ,%%ymm8 ,%%ymm8 		\n\t"\
			"vroundpd	$0,%%ymm0,%%ymm6				\n\t	vroundpd	$0,%%ymm8 ,%%ymm14		\n\t"\
			"vsubpd		%%ymm6,%%ymm0,%%ymm0			\n\t	vsubpd	%%ymm14,%%ymm8 ,%%ymm8 		\n\t"\
			"vmulpd		%%ymm5,%%ymm0,%%ymm0			\n\t	vmulpd	%%ymm5 ,%%ymm8 ,%%ymm8 		\n\t"\
			"vmulpd		%%ymm15,%%ymm6,%%ymm6			\n\t	vmulpd	%%ymm15,%%ymm14,%%ymm14		\n\t"\
			"/* Digit 4: */\n\t"\
			"vaddpd		%%ymm6,%%ymm3,%%ymm3			\n\t	vaddpd	%%ymm14,%%ymm11,%%ymm11		\n\t"\
			"vmovaps	%%ymm3,%%ymm1					\n\t	vmovaps	%%ymm11,%%ymm9				\n\t"\
			"vroundpd	$0,%%ymm3,%%ymm3				\n\t	vroundpd	$0,%%ymm11,%%ymm11		\n\t"\
			"vsubpd		%%ymm3,%%ymm1,%%ymm1			\n\t	vsubpd	%%ymm11,%%ymm9 ,%%ymm9 		\n\t"\
			"vmulpd		%%ymm5,%%ymm1,%%ymm1			\n\t	vmulpd	%%ymm5 ,%%ymm9 ,%%ymm9 		\n\t"\
			"/* If h < l, calculate h-l+q; otherwise h-l. Use leading 52 bits to approximate the full 78-bit compare. Result is in [0, q). */\n\t"\
			"vmovaps	%%ymm5,%%ymm13	/* Need a copy of two26f, stored in %%ymm5 */\n\t"\
			"vmovaps	0x200(%%rdx),%%ymm2				\n\t	vmovaps	0x220(%%rdx),%%ymm10		\n\t"/* Restore hi1 from ax1 */\
			"vmovaps	%%ymm2,%%ymm6					\n\t	vmovaps	%%ymm10,%%ymm14				\n\t"\
			"vmovaps	0xc0(%%rdx),%%ymm4				\n\t	vmovaps	0xe0(%%rdx),%%ymm12			\n\t"\
			"vmulpd		%%ymm5,%%ymm3,%%ymm3			\n\t	vmulpd	%%ymm5 ,%%ymm11,%%ymm11		\n\t"\
			"vmulpd	   (%%rdx),%%ymm5,%%ymm5			\n\t	vmulpd	0x20(%%rdx),%%ymm13,%%ymm13	\n\t"\
			"vaddpd		%%ymm1,%%ymm3,%%ymm3			\n\t	vaddpd	%%ymm9 ,%%ymm11,%%ymm11		\n\t"\
			"vmovaps	0x1c0(%%rdx),%%ymm1				\n\t	vmovaps	0x1e0(%%rdx),%%ymm9			\n\t"/* Restore hi0 from ax0 */\
			"vcmppd	   $0x1,%%ymm3,%%ymm6,%%ymm6		\n\t	vcmppd $0x1,%%ymm11,%%ymm14,%%ymm14	\n\t"\
			"vsubpd		%%ymm0,%%ymm1,%%ymm1			\n\t	vsubpd	%%ymm8 ,%%ymm9 ,%%ymm9 		\n\t"\
			"vsubpd		%%ymm3,%%ymm2,%%ymm2			\n\t	vsubpd	%%ymm11,%%ymm10,%%ymm10		\n\t"\
			"vmovaps	%%ymm4,%%ymm0					\n\t	vmovaps	%%ymm12,%%ymm8				\n\t"\
			"vandpd		%%ymm6,%%ymm4,%%ymm4			\n\t	vandpd	%%ymm14,%%ymm12,%%ymm12		\n\t"\
			"vaddpd		%%ymm4,%%ymm2,%%ymm2			\n\t	vaddpd	%%ymm12,%%ymm10,%%ymm10		\n\t"\
			"vandpd		%%ymm5,%%ymm6,%%ymm6			\n\t	vandpd	%%ymm13,%%ymm14,%%ymm14		\n\t"\
			"vaddpd		%%ymm6,%%ymm1,%%ymm1			\n\t	vaddpd	%%ymm14,%%ymm9 ,%%ymm9 		\n\t"\
			"/* qlo26 in ymm5, qhi52 in ymm0*/\n\t"\
		"/* if((pshift >> j) & (uint64)1) { */\n\t"\
			"movslq	%[__pshift],%%rax		\n\t"\
			"movslq	%[__j],%%rcx			\n\t"\
			"shrq	%%cl,%%rax				\n\t"\
			"andq	$0x1,%%rax				\n\t"\
		"je	twopmodq78_3wdq8	\n\t"\
			"/* Double: Use that high part < 2^52 (strictly >= 0): xhi52,xlo26 in ymm2,ymm1; qhi52,qlo26 in ymm0,ymm5 */\n\t"\
			"	vmovaps		%%ymm0,%%ymm6				\n\t	vmovaps		%%ymm8 ,%%ymm14			/* cpy of qhi */\n\t"\
			"	vaddpd		%%ymm2,%%ymm2,%%ymm2		\n\t	vaddpd		%%ymm10,%%ymm10,%%ymm10	/* top 52 bits */\n\t"\
			"	vaddpd		%%ymm1,%%ymm1,%%ymm1		\n\t	vaddpd		%%ymm9 ,%%ymm9 ,%%ymm9 	/* low 26 bits */\n\t"\
			"/* If x > q, subtract q: */\n\t"\
			"	vcmppd $0x2,%%ymm2,%%ymm6,%%ymm6		\n\t	vcmppd $0x2,%%ymm10,%%ymm14,%%ymm14	/* bitmask = (qhi <= xhi) */\n\t"\
			"	vandpd		%%ymm6,%%ymm0,%%ymm0		\n\t	vandpd		%%ymm14,%%ymm8 ,%%ymm8 	/* qhi52 & bitmask */\n\t"\
			"	vandpd		%%ymm6,%%ymm5,%%ymm5		\n\t	vandpd		%%ymm14,%%ymm13,%%ymm13	/* qlo26 & bitmask */\n\t"\
			"	vsubpd		%%ymm0,%%ymm2,%%ymm2		\n\t	vsubpd		%%ymm8 ,%%ymm10,%%ymm10	/* x mod q, top 52 bits */\n\t"\
			"	vsubpd		%%ymm5,%%ymm1,%%ymm1		\n\t	vsubpd		%%ymm13,%%ymm9 ,%%ymm9 	/* x mod q, low 26 bits */\n\t"\
		"twopmodq78_3wdq8:	\n\t"\
		"/* } */\n\t"\
			"/* Normalize the result: */\n\t"\
			"vmovaps	0x320(%%rdx),%%ymm4				\n\t	vmovaps	0x320(%%rdx),%%ymm12		\n\t"\
			"vmovaps		%%ymm15,%%ymm3				\n\t	vmovaps		%%ymm15,%%ymm11			\n\t"\
			"vmulpd		%%ymm3,%%ymm1,%%ymm1			\n\t	vmulpd	%%ymm11,%%ymm9 ,%%ymm9 		\n\t"\
			"vmulpd		%%ymm3,%%ymm2,%%ymm2			\n\t	vmulpd	%%ymm11,%%ymm10,%%ymm10		\n\t"\
			"vmovaps	%%ymm1,%%ymm0					\n\t	vmovaps	%%ymm9 ,%%ymm8 				\n\t"\
			"vroundpd	$0,%%ymm1,%%ymm1				\n\t	vroundpd	$0,%%ymm9 ,%%ymm9		\n\t"\
			"vsubpd		%%ymm1,%%ymm0,%%ymm0			\n\t	vsubpd	%%ymm9 ,%%ymm8 ,%%ymm8 		\n\t"\
			"vmulpd		%%ymm4,%%ymm0,%%ymm0			\n\t	vmulpd	%%ymm12,%%ymm8 ,%%ymm8 		\n\t"\
			"vmulpd		%%ymm3,%%ymm1,%%ymm1			\n\t	vmulpd	%%ymm11,%%ymm9 ,%%ymm9 		\n\t"\
			"vaddpd		%%ymm1,%%ymm2,%%ymm2			\n\t	vaddpd	%%ymm9 ,%%ymm10,%%ymm10		\n\t"\
			"vmovaps	%%ymm2,%%ymm1					\n\t	vmovaps	%%ymm10,%%ymm9				\n\t"\
			"vroundpd	$0,%%ymm2,%%ymm2				\n\t	vroundpd	$0,%%ymm10,%%ymm10		\n\t"\
			"vsubpd		%%ymm2,%%ymm1,%%ymm1			\n\t	vsubpd	%%ymm10,%%ymm9 ,%%ymm9 		\n\t"\
			"vmulpd		%%ymm4,%%ymm1,%%ymm1			\n\t	vmulpd	%%ymm12,%%ymm9 ,%%ymm9 		\n\t"\
			"/* Move high 2 words of result into input registers expected by start of loop body: */\n\t"\
			"vmovaps	%%ymm2,%%ymm4 /* fx2 */			\n\t	vmovaps	%%ymm10,%%ymm12	/* hx2 */	\n\t"\
			"vmovaps	%%ymm1,%%ymm2 /* fx1 */			\n\t	vmovaps	%%ymm9 ,%%ymm10	/* hx1 */	\n\t"\
			:					/* outputs: none */\
			: [__fq0] "m" (Xfq0)	/* All inputs from memory addresses here */\
			 ,[__pshift] "m" (Xpshift)	\
			 ,[__j]		"m" (Xj)		\
			: "cl","rax","rcx","rdx","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
		);\
		}

  #else	// SSE2:

		#define SSE2_twopmodq78_modmul_q8(Xaq0,Xaqinv0,Xax0,Xtwo26i,Xpshift,Xj)\
		{\
		__asm__ volatile (\
		/* SQR_LOHI78_3WORD_DOUBLE_q8(). Inputs [a|c|e|g]lo0-2 enter in xmm[0-2,4-6,8-10,12-14], resp.: */\
			"movq	%[__ax0],%%rax		\n\t"\
			"movq	%[__two26i],%%rbx	\n\t"\
			"movaps	-0x20(%%rbx),%%xmm3	\n\t"/* two13i */\
			/*** STREAM 0 ***/				/*** STREAM 1 ***/						/*** STREAM 2 ***/						/*** STREAM 3 ***/					\
			"/* fx0,1,2 assumed in xmm0,1,2 on loop entry */\n\t"\
			"mulpd	%%xmm3,%%xmm0			\n\t	mulpd	%%xmm3,%%xmm4			\n\t	mulpd	%%xmm3,%%xmm8			\n\t	mulpd	%%xmm3,%%xmm12		\n\t"/* scale [fx0 . 2^-13] */\
			"mulpd	%%xmm3,%%xmm1			\n\t	mulpd	%%xmm3,%%xmm5			\n\t	mulpd	%%xmm3,%%xmm9			\n\t	mulpd	%%xmm3,%%xmm13		\n\t"/* scale [fx1 . 2^-13] */\
			"mulpd	%%xmm3,%%xmm2			\n\t	mulpd	%%xmm3,%%xmm6			\n\t	mulpd	%%xmm3,%%xmm10			\n\t	mulpd	%%xmm3,%%xmm14		\n\t"/* scale [fx2 . 2^-13] */\
			"movaps	%%xmm0,%%xmm3			\n\t	movaps	%%xmm4,%%xmm7			\n\t	movaps	%%xmm8,%%xmm11			\n\t	movaps	%%xmm12,%%xmm15		\n\t"/* cpy of fx0 */\
			"addpd	%%xmm3,%%xmm3			\n\t	addpd	%%xmm7,%%xmm7			\n\t	addpd	%%xmm11,%%xmm11			\n\t	addpd	%%xmm15,%%xmm15		\n\t"/* 2.fx0 */\
			"movaps	%%xmm2,0x80(%%rax)		\n\t	movaps	%%xmm6,0x90(%%rax)		\n\t	movaps	%%xmm10,0xa0(%%rax)		\n\t	movaps	%%xmm14,0xb0(%%rax)	\n\t"/* Store fx2 to free up a register */\
			"/* Digit 0: */\n\t"\
			"mulpd	%%xmm0,%%xmm0			\n\t	mulpd	%%xmm4,%%xmm4			\n\t	mulpd	%%xmm8,%%xmm8			\n\t	mulpd	%%xmm12,%%xmm12		\n\t"/* [  fx0 *= fx0] / 2^26 */\
			"roundpd $0,%%xmm0,%%xmm2		\n\t	roundpd $0,%%xmm4,%%xmm6		\n\t	roundpd $0,%%xmm8,%%xmm10		\n\t	roundpd $0,%%xmm12,%%xmm14	\n\t"/* Init: FCY = dnint(fx0*fx0) */\
			"subpd	%%xmm2,%%xmm0			\n\t	subpd	%%xmm6,%%xmm4			\n\t	subpd	%%xmm10,%%xmm8			\n\t	subpd	%%xmm14,%%xmm12		\n\t"/* fprod0 -= FCY */\
			"mulpd	-0x10(%%rbx),%%xmm0		\n\t	mulpd	-0x10(%%rbx),%%xmm4		\n\t	mulpd	-0x10(%%rbx),%%xmm8		\n\t	mulpd	-0x10(%%rbx),%%xmm12\n\t"/* fprod0 *= two26f */\
			"mulpd	     (%%rbx),%%xmm2		\n\t	mulpd	     (%%rbx),%%xmm6		\n\t	mulpd	     (%%rbx),%%xmm10	\n\t	mulpd	     (%%rbx),%%xmm14\n\t"/* FCY    *= two26i */\
			"movaps	%%xmm0,    (%%rax)		\n\t	movaps	%%xmm4,0x10(%%rax)		\n\t	movaps	%%xmm8,0x20(%%rax)		\n\t	movaps	%%xmm12,0x30(%%rax)	\n\t"/* Store fprod0 to free up a register */\
			"/* Digit 1: */\n\t"\
			"movaps	%%xmm1,%%xmm0			\n\t	movaps	%%xmm5,%%xmm4			\n\t	movaps	%%xmm9,%%xmm8			\n\t	movaps	%%xmm13,%%xmm12		\n\t"/* cpy of fx1 */\
			"mulpd	%%xmm3,%%xmm0			\n\t	mulpd	%%xmm7,%%xmm4			\n\t	mulpd	%%xmm11,%%xmm8			\n\t	mulpd	%%xmm15,%%xmm12		\n\t"/* [fx1 *= 2.fx0] / 2^26 */\
			"addpd	%%xmm2,%%xmm0			\n\t	addpd	%%xmm6,%%xmm4			\n\t	addpd	%%xmm10,%%xmm8			\n\t	addpd	%%xmm14,%%xmm12		\n\t"/* fprod1 = 2.fx0*fx1 + FCY */\
			"roundpd $0,%%xmm0,%%xmm2		\n\t	roundpd $0,%%xmm4,%%xmm6		\n\t	roundpd $0,%%xmm8,%%xmm10		\n\t	roundpd $0,%%xmm12,%%xmm14	\n\t"/* FCY */\
			"subpd	%%xmm2,%%xmm0			\n\t	subpd	%%xmm6,%%xmm4			\n\t	subpd	%%xmm10,%%xmm8			\n\t	subpd	%%xmm14,%%xmm12		\n\t"/* fprod1 -= FCY */\
			"mulpd	-0x10(%%rbx),%%xmm0		\n\t	mulpd	-0x10(%%rbx),%%xmm4		\n\t	mulpd	-0x10(%%rbx),%%xmm8		\n\t	mulpd	-0x10(%%rbx),%%xmm12\n\t"/* fprod1 *= two26f */\
			"mulpd	     (%%rbx),%%xmm2		\n\t	mulpd	     (%%rbx),%%xmm6		\n\t	mulpd	     (%%rbx),%%xmm10	\n\t	mulpd	     (%%rbx),%%xmm14\n\t"/* FCY    *= two26i */\
			"movaps	%%xmm0,0x40(%%rax)		\n\t	movaps	%%xmm4,0x50(%%rax)		\n\t	movaps	%%xmm8,0x60(%%rax)		\n\t	movaps	%%xmm12,0x70(%%rax)	\n\t"/* Store fprod1 to free up a register */\
			"/* Digit 2: */\n\t"\
			"movaps	%%xmm1,%%xmm0			\n\t	movaps	%%xmm5,%%xmm4			\n\t	movaps	%%xmm9,%%xmm8			\n\t	movaps	%%xmm13,%%xmm12		\n\t"/* cpy of fx1 */\
			"addpd	%%xmm0,%%xmm0			\n\t	addpd	%%xmm4,%%xmm4			\n\t	addpd	%%xmm8,%%xmm8			\n\t	addpd	%%xmm12,%%xmm12		\n\t"/* 2.fx1 */\
			"mulpd	%%xmm1,%%xmm1			\n\t	mulpd	%%xmm5,%%xmm5			\n\t	mulpd	%%xmm9,%%xmm9			\n\t	mulpd	%%xmm13,%%xmm13		\n\t"/* [  fx1 *= fx1] / 2^26 */\
			"mulpd	0x80(%%rax),%%xmm3		\n\t	mulpd	0x90(%%rax),%%xmm7		\n\t	mulpd	0xa0(%%rax),%%xmm11		\n\t	mulpd	0xb0(%%rax),%%xmm15	\n\t"/* [2.fx0 *= fx2] / 2^26 */\
			"addpd	%%xmm2,%%xmm1			\n\t	addpd	%%xmm6,%%xmm5			\n\t	addpd	%%xmm10,%%xmm9			\n\t	addpd	%%xmm14,%%xmm13		\n\t"/* fx1*fx1 += FCY */\
			"addpd	%%xmm3,%%xmm1			\n\t	addpd	%%xmm7,%%xmm5			\n\t	addpd	%%xmm11,%%xmm9			\n\t	addpd	%%xmm15,%%xmm13		\n\t"/* fprod2 = 2.fx0*fx2 + fx1*fx1 + FCY */\
			"roundpd $1,%%xmm1,%%xmm2		\n\t	roundpd $1,%%xmm5,%%xmm6		\n\t	roundpd $1,%%xmm9,%%xmm10		\n\t	roundpd $1,%%xmm13,%%xmm14	\n\t"/* FCY = floor(fprod2) */\
			"subpd	%%xmm2,%%xmm1			\n\t	subpd	%%xmm6,%%xmm5			\n\t	subpd	%%xmm10,%%xmm9			\n\t	subpd	%%xmm14,%%xmm13		\n\t"/* fprod2 -= FCY */\
			"mulpd	-0x10(%%rbx),%%xmm1		\n\t	mulpd	-0x10(%%rbx),%%xmm5		\n\t	mulpd	-0x10(%%rbx),%%xmm9		\n\t	mulpd	-0x10(%%rbx),%%xmm13\n\t"/* fprod2 *= two26f */\
			"mulpd	     (%%rbx),%%xmm2		\n\t	mulpd	     (%%rbx),%%xmm6		\n\t	mulpd	     (%%rbx),%%xmm10	\n\t	mulpd	     (%%rbx),%%xmm14\n\t"/* FCY    *= two26i */\
			"/* Digit 3: */\n\t"\
			"movaps	0x80(%%rax),%%xmm3		\n\t	movaps	0x90(%%rax),%%xmm7		\n\t	movaps	0xa0(%%rax),%%xmm11		\n\t	movaps	0xb0(%%rax),%%xmm15	\n\t"/* Reload fx2 */\
			"mulpd	%%xmm3,%%xmm0			\n\t	mulpd	%%xmm7,%%xmm4			\n\t	mulpd	%%xmm11,%%xmm8			\n\t	mulpd	%%xmm15,%%xmm12		\n\t"/* [2.fx1 *= fx2] / 2^26 */\
			"addpd	%%xmm2,%%xmm0			\n\t	addpd	%%xmm6,%%xmm4			\n\t	addpd	%%xmm10,%%xmm8			\n\t	addpd	%%xmm14,%%xmm12		\n\t"/* fprod3 = 2.fx1*fx2 + FCY */\
			"roundpd $0,%%xmm0,%%xmm2		\n\t	roundpd $0,%%xmm4,%%xmm6		\n\t	roundpd $0,%%xmm8,%%xmm10		\n\t	roundpd $0,%%xmm12,%%xmm14	\n\t"\
			"subpd	%%xmm2,%%xmm0			\n\t	subpd	%%xmm6,%%xmm4			\n\t	subpd	%%xmm10,%%xmm8			\n\t	subpd	%%xmm14,%%xmm12		\n\t"/* fprod3 -= FCY */\
			"mulpd	-0x10(%%rbx),%%xmm0		\n\t	mulpd	-0x10(%%rbx),%%xmm4		\n\t	mulpd	-0x10(%%rbx),%%xmm8		\n\t	mulpd	-0x10(%%rbx),%%xmm12\n\t"/* fprod3 *= two26f */\
			"mulpd	     (%%rbx),%%xmm2		\n\t	mulpd	     (%%rbx),%%xmm6		\n\t	mulpd	     (%%rbx),%%xmm10	\n\t	mulpd	     (%%rbx),%%xmm14\n\t"/* FCY    *= two26i */\
			"/* Digit 4: */\n\t"\
			"mulpd	%%xmm3,%%xmm3			\n\t	mulpd	%%xmm7,%%xmm7			\n\t	mulpd	%%xmm11,%%xmm11			\n\t	mulpd	%%xmm15,%%xmm15		\n\t"/* [  fx2 *= fx2] / 2^26 */\
			"addpd	%%xmm2,%%xmm3			\n\t	addpd	%%xmm6,%%xmm7			\n\t	addpd	%%xmm10,%%xmm11			\n\t	addpd	%%xmm14,%%xmm15		\n\t"/* fprod4 = fx2*fx2 + fcy */\
			"roundpd $0,%%xmm3,%%xmm2		\n\t	roundpd $0,%%xmm7,%%xmm6		\n\t	roundpd $0,%%xmm11,%%xmm10		\n\t	roundpd $0,%%xmm15,%%xmm14	\n\t"\
			"subpd	%%xmm2,%%xmm3			\n\t	subpd	%%xmm6,%%xmm7			\n\t	subpd	%%xmm10,%%xmm11			\n\t	subpd	%%xmm14,%%xmm15		\n\t"/* fprod4 -= fcy */\
			"mulpd	-0x10(%%rbx),%%xmm3		\n\t	mulpd	-0x10(%%rbx),%%xmm7		\n\t	mulpd	-0x10(%%rbx),%%xmm11	\n\t	mulpd	-0x10(%%rbx),%%xmm15\n\t"/* fprod4 *= two26f */\
			"/* Digit 5 = the carry. CY in xmm2; fhi0,1 in xmm0,xmm3 */\n\t"\
			"mulpd	-0x10(%%rbx),%%xmm2		\n\t	mulpd	-0x10(%%rbx),%%xmm6		\n\t	mulpd	-0x10(%%rbx),%%xmm10	\n\t	mulpd	-0x10(%%rbx),%%xmm14\n\t"/* fhi2 * two26f */\
			"addpd	%%xmm2,%%xmm3			\n\t	addpd	%%xmm6,%%xmm7			\n\t	addpd	%%xmm10,%%xmm11			\n\t	addpd	%%xmm14,%%xmm15		\n\t"/* fhi, top 52 bits; xmm2 FREE */\
			"movaps	%%xmm0,0xc0(%%rax)		\n\t	movaps	%%xmm4,0xd0(%%rax)		\n\t	movaps	%%xmm8,0xe0(%%rax)		\n\t	movaps	%%xmm12,0xf0(%%rax)	\n\t"/* Store fhi0,1 = fprod3,4 to free up 2 more registers */\
			"movaps	%%xmm3,0x100(%%rax)		\n\t	movaps	%%xmm7,0x110(%%rax)		\n\t	movaps	%%xmm11,0x120(%%rax)	\n\t	movaps	%%xmm15,0x130(%%rax)\n\t"/* Recall: fx0,1 in eax+0x[0,2]0, fx2 in xmm1. */\
		/* MULL78 section below needs flo0,1,2 in xmm0,1,2: */\
		/* flo = MULL78_3WORD_DOUBLE_q2(flo, fqinv): */\
			"movq	%[__aqinv0],%%rdx	\n\t"\
			"/* Digit 0: */\n\t"\
			"movaps	    (%%rax),%%xmm0		\n\t	movaps	0x10(%%rax),%%xmm4		\n\t	movaps	0x20(%%rax),%%xmm8		\n\t	movaps	0x30(%%rax),%%xmm12	\n\t"/* Reload flo0 */\
			"movaps	%%xmm1,%%xmm2			\n\t	movaps	%%xmm5,%%xmm6			\n\t	movaps	%%xmm9,%%xmm10			\n\t	movaps	%%xmm13,%%xmm14		\n\t"/* Move flo2 into xmm2 */\
			"movaps	%%xmm0,%%xmm3			\n\t	movaps	%%xmm4,%%xmm7			\n\t	movaps	%%xmm8,%%xmm11			\n\t	movaps	%%xmm12,%%xmm15		\n\t"/* xmm3 = cpy of x0 */\
			"mulpd	    (%%rdx),%%xmm0		\n\t	mulpd	0x10(%%rdx),%%xmm4		\n\t	mulpd	0x20(%%rdx),%%xmm8		\n\t	mulpd	0x30(%%rdx),%%xmm12	\n\t"/* fprod0 = x0*y0 */\
			"roundpd $0,%%xmm0,%%xmm1		\n\t	roundpd $0,%%xmm4,%%xmm5		\n\t	roundpd $0,%%xmm8,%%xmm9		\n\t	roundpd $0,%%xmm12,%%xmm13	\n\t"\
			"subpd	%%xmm1,%%xmm0			\n\t	subpd	%%xmm5,%%xmm4			\n\t	subpd	%%xmm9,%%xmm8			\n\t	subpd	%%xmm13,%%xmm12		\n\t"/* fprod0 -= fcy */\
			"mulpd	-0x10(%%rbx),%%xmm0		\n\t	mulpd	-0x10(%%rbx),%%xmm4		\n\t	mulpd	-0x10(%%rbx),%%xmm8		\n\t	mulpd	-0x10(%%rbx),%%xmm12\n\t"/* fprod0 *= two26f */\
			"mulpd	     (%%rbx),%%xmm1		\n\t	mulpd	     (%%rbx),%%xmm5		\n\t	mulpd	     (%%rbx),%%xmm9		\n\t	mulpd	     (%%rbx),%%xmm13\n\t"/* fcy    *= two26i */\
			"movaps	%%xmm0,    (%%rax)		\n\t	movaps	%%xmm4,0x10(%%rax)		\n\t	movaps	%%xmm8,0x20(%%rax)		\n\t	movaps	%%xmm12,0x30(%%rax)	\n\t"/* Store fprod0 to free up a register */\
			"movaps	%%xmm3,%%xmm0			\n\t	movaps	%%xmm7,%%xmm4			\n\t	movaps	%%xmm11,%%xmm8			\n\t	movaps	%%xmm15,%%xmm12		\n\t"/* xmm0 = cpy of x0 */\
			"/* Digit 1: */\n\t"\
			"mulpd	0x40(%%rdx),%%xmm3		\n\t	mulpd	0x50(%%rdx),%%xmm7		\n\t	mulpd	0x60(%%rdx),%%xmm11		\n\t	mulpd	0x70(%%rdx),%%xmm15	\n\t"/* x0*y1 */\
			"addpd	%%xmm1,%%xmm3			\n\t	addpd	%%xmm5,%%xmm7			\n\t	addpd	%%xmm9,%%xmm11			\n\t	addpd	%%xmm13,%%xmm15		\n\t"/* x0*y1 + fcy; xmm1 FREE */\
			"movaps	0x40(%%rax),%%xmm1		\n\t	movaps	0x50(%%rax),%%xmm5		\n\t	movaps	0x60(%%rax),%%xmm9		\n\t	movaps	0x70(%%rax),%%xmm13	\n\t"/* Reload flo1 */\
			"mulpd	    (%%rdx),%%xmm1		\n\t	mulpd	0x10(%%rdx),%%xmm5		\n\t	mulpd	0x20(%%rdx),%%xmm9		\n\t	mulpd	0x30(%%rdx),%%xmm13	\n\t"/* x1*y0 */\
			"addpd	%%xmm3,%%xmm1			\n\t	addpd	%%xmm7,%%xmm5			\n\t	addpd	%%xmm11,%%xmm9			\n\t	addpd	%%xmm15,%%xmm13		\n\t"/* fprod1 = x0*y1 + x1*y0; xmm3 FREE */\
			"roundpd $0,%%xmm1,%%xmm3		\n\t	roundpd $0,%%xmm5,%%xmm7		\n\t	roundpd $0,%%xmm9,%%xmm11		\n\t	roundpd $0,%%xmm13,%%xmm15	\n\t"\
			"subpd	%%xmm3,%%xmm1			\n\t	subpd	%%xmm7,%%xmm5			\n\t	subpd	%%xmm11,%%xmm9			\n\t	subpd	%%xmm15,%%xmm13		\n\t"/* fprod1 -= fcy */\
			"mulpd	-0x10(%%rbx),%%xmm1		\n\t	mulpd	-0x10(%%rbx),%%xmm5		\n\t	mulpd	-0x10(%%rbx),%%xmm9		\n\t	mulpd	-0x10(%%rbx),%%xmm13\n\t"/* fprod0 *= two26f */\
			"mulpd	     (%%rbx),%%xmm3		\n\t	mulpd	     (%%rbx),%%xmm7		\n\t	mulpd	     (%%rbx),%%xmm11	\n\t	mulpd	     (%%rbx),%%xmm15\n\t"/* fcy    *= two26i */\
			"movaps	%%xmm1,0x80(%%rax)		\n\t	movaps	%%xmm5,0x90(%%rax)		\n\t	movaps	%%xmm9,0xa0(%%rax)		\n\t	movaps	%%xmm13,0xb0(%%rax)	\n\t"/* Store fprod1 [in unused flo2 slot] to free up a register */\
			"movaps	0x40(%%rax),%%xmm1		\n\t	movaps	0x50(%%rax),%%xmm5		\n\t	movaps	0x60(%%rax),%%xmm9		\n\t	movaps	0x70(%%rax),%%xmm13	\n\t"/* Reload flo1 */\
			"/* Digit 2: Require output to be nonnegative, so leave MSW unbalanced: */\n\t"\
			"mulpd	    (%%rdx),%%xmm2		\n\t	mulpd	0x10(%%rdx),%%xmm6		\n\t	mulpd	0x20(%%rdx),%%xmm10		\n\t	mulpd	0x30(%%rdx),%%xmm14	\n\t"/* x2*y0 */\
			"mulpd	0x40(%%rdx),%%xmm1		\n\t	mulpd	0x50(%%rdx),%%xmm5		\n\t	mulpd	0x60(%%rdx),%%xmm9		\n\t	mulpd	0x70(%%rdx),%%xmm13	\n\t"/* x1*y1 */\
			"mulpd	0x80(%%rdx),%%xmm0		\n\t	mulpd	0x90(%%rdx),%%xmm4		\n\t	mulpd	0xa0(%%rdx),%%xmm8		\n\t	mulpd	0xb0(%%rdx),%%xmm12	\n\t"/* x0*y2 */\
			"addpd	%%xmm1,%%xmm2			\n\t	addpd	%%xmm5,%%xmm6			\n\t	addpd	%%xmm9,%%xmm10			\n\t	addpd	%%xmm13,%%xmm14		\n\t"/* x1*y1 + x2*y0; xmm1 FREE */\
			"addpd	%%xmm0,%%xmm3			\n\t	addpd	%%xmm4,%%xmm7			\n\t	addpd	%%xmm8,%%xmm11			\n\t	addpd	%%xmm12,%%xmm15		\n\t"/* x0*y2 + fcy; xmm0 FREE */\
			"addpd	%%xmm3,%%xmm2			\n\t	addpd	%%xmm7,%%xmm6			\n\t	addpd	%%xmm11,%%xmm10			\n\t	addpd	%%xmm15,%%xmm14		\n\t"/* fprod2; xmm3 FREE */\
			"roundpd $1,%%xmm2,%%xmm3		\n\t	roundpd $1,%%xmm6,%%xmm7		\n\t	roundpd $1,%%xmm10,%%xmm11		\n\t	roundpd $1,%%xmm14,%%xmm15	\n\t"/* fcy = floor(fprod2) */\
			"subpd	%%xmm3,%%xmm2			\n\t	subpd	%%xmm7,%%xmm6			\n\t	subpd	%%xmm11,%%xmm10			\n\t	subpd	%%xmm15,%%xmm14		\n\t"/* fprod2 -= fcy */\
			"mulpd	-0x10(%%rbx),%%xmm2		\n\t	mulpd	-0x10(%%rbx),%%xmm6		\n\t	mulpd	-0x10(%%rbx),%%xmm10	\n\t	mulpd	-0x10(%%rbx),%%xmm14\n\t"/* fprod2 *= two26f */\
			"movaps	    (%%rax),%%xmm0		\n\t	movaps	0x10(%%rax),%%xmm4		\n\t	movaps	0x20(%%rax),%%xmm8		\n\t	movaps	0x30(%%rax),%%xmm12	\n\t"/* Reload fprod0 */\
			"movaps	0x80(%%rax),%%xmm1		\n\t	movaps	0x90(%%rax),%%xmm5		\n\t	movaps	0xa0(%%rax),%%xmm9		\n\t	movaps	0xb0(%%rax),%%xmm13	\n\t"/* Reload fprod1 [in flo2 slot] */\
		/* MULH96(q,lo,lo) --> lo =(q*lo)/2^78: flo0,1,2 enter in xmm0,1,2. In the comments x = q, y = lo: */\
			"movq	%[__aq0],%%rdx		\n\t"\
			"/* Digit 0: */\n\t"\
			"movaps	%%xmm0,%%xmm3			\n\t	movaps	%%xmm4,%%xmm7			\n\t	movaps	%%xmm8,%%xmm11			\n\t	movaps	%%xmm12,%%xmm15		\n\t"/* xmm3 = cpy of y0 */\
			"mulpd	    (%%rdx),%%xmm0		\n\t	mulpd	0x10(%%rdx),%%xmm4		\n\t	mulpd	0x20(%%rdx),%%xmm8		\n\t	mulpd	0x30(%%rdx),%%xmm12	\n\t"/* fprod0 = y0*x0 */\
			"mulpd	    (%%rbx),%%xmm0		\n\t	mulpd	    (%%rbx),%%xmm4		\n\t	mulpd	    (%%rbx),%%xmm8		\n\t	mulpd	    (%%rbx),%%xmm12	\n\t"/* CY    *= two26i */\
			"/* Digit 1: */\n\t"\
			"mulpd	0x40(%%rdx),%%xmm3		\n\t	mulpd	0x50(%%rdx),%%xmm7		\n\t	mulpd	0x60(%%rdx),%%xmm11		\n\t	mulpd	0x70(%%rdx),%%xmm15	\n\t"/* y0 *= x1 */\
			"addpd	%%xmm0,%%xmm3			\n\t	addpd	%%xmm4,%%xmm7			\n\t	addpd	%%xmm8,%%xmm11			\n\t	addpd	%%xmm12,%%xmm15		\n\t"/* y0*x1 + CY; xmm0 FREE */\
			"mulpd	    (%%rdx),%%xmm1		\n\t	mulpd	0x10(%%rdx),%%xmm5		\n\t	mulpd	0x20(%%rdx),%%xmm9		\n\t	mulpd	0x30(%%rdx),%%xmm13	\n\t"/* y1 *= x0 */\
			"addpd	%%xmm3,%%xmm1			\n\t	addpd	%%xmm7,%%xmm5			\n\t	addpd	%%xmm11,%%xmm9			\n\t	addpd	%%xmm15,%%xmm13		\n\t"/* fprod1 = x0*y1 + x1*y0; xmm3 FREE */\
			"mulpd	    (%%rbx),%%xmm1		\n\t	mulpd	    (%%rbx),%%xmm5		\n\t	mulpd	    (%%rbx),%%xmm9		\n\t	mulpd	    (%%rbx),%%xmm13	\n\t"/* CY    *= two26i */\
			"/* Digit 2: Require output to be nonnegative, so leave MSW unbalanced: */\n\t"\
			"movaps	%%xmm2,%%xmm3			\n\t	movaps	%%xmm6,%%xmm7			\n\t	movaps	%%xmm10,%%xmm11			\n\t	movaps	%%xmm14,%%xmm15		\n\t"/* xmm3 = cpy of y2 */\
			"movaps	0x80(%%rax),%%xmm0		\n\t	movaps	0x90(%%rax),%%xmm4		\n\t	movaps	0xa0(%%rax),%%xmm8		\n\t	movaps	0xb0(%%rax),%%xmm12	\n\t"/* Reload y1 [in flo2 slot] */\
			"mulpd	    (%%rdx),%%xmm2		\n\t	mulpd	0x10(%%rdx),%%xmm6		\n\t	mulpd	0x20(%%rdx),%%xmm10		\n\t	mulpd	0x30(%%rdx),%%xmm14	\n\t"/* y2 *= x0 */\
			"mulpd	0x40(%%rdx),%%xmm0		\n\t	mulpd	0x50(%%rdx),%%xmm4		\n\t	mulpd	0x60(%%rdx),%%xmm8		\n\t	mulpd	0x70(%%rdx),%%xmm12	\n\t"/* y1 *= x1 */\
			"addpd	%%xmm1,%%xmm2			\n\t	addpd	%%xmm5,%%xmm6			\n\t	addpd	%%xmm9,%%xmm10			\n\t	addpd	%%xmm13,%%xmm14		\n\t"/* y2*x0 + CY; xmm1 FREE */\
			"movaps	    (%%rax),%%xmm1		\n\t	movaps	0x10(%%rax),%%xmm5		\n\t	movaps	0x20(%%rax),%%xmm9		\n\t	movaps	0x30(%%rax),%%xmm13	\n\t"/* Reload y0 */\
			"mulpd	0x80(%%rdx),%%xmm1		\n\t	mulpd	0x90(%%rdx),%%xmm5		\n\t	mulpd	0xa0(%%rdx),%%xmm9		\n\t	mulpd	0xb0(%%rdx),%%xmm13	\n\t"/* y0 *= x2 */\
			"addpd	%%xmm0,%%xmm1			\n\t	addpd	%%xmm4,%%xmm5			\n\t	addpd	%%xmm8,%%xmm9			\n\t	addpd	%%xmm12,%%xmm13		\n\t"/* y1*x1 + y0*x2; xmm0 FREE */\
			"addpd	%%xmm1,%%xmm2			\n\t	addpd	%%xmm5,%%xmm6			\n\t	addpd	%%xmm9,%%xmm10			\n\t	addpd	%%xmm13,%%xmm14		\n\t"/* fprod2; xmm1 FREE */\
			"roundpd $1,%%xmm2,%%xmm2		\n\t	roundpd $1,%%xmm6,%%xmm6		\n\t	roundpd $1,%%xmm10,%%xmm10		\n\t	roundpd $1,%%xmm14,%%xmm14	\n\t"/* CY = floor(fprod2) */\
			"mulpd	    (%%rbx),%%xmm2		\n\t	mulpd	    (%%rbx),%%xmm6		\n\t	mulpd	    (%%rbx),%%xmm10		\n\t	mulpd	    (%%rbx),%%xmm14	\n\t"/* CY    *= two26i */\
		/* Precompute partial products needed for upper half: */\
			"movaps	%%xmm3,%%xmm1			\n\t	movaps	%%xmm7,%%xmm5			\n\t	movaps	%%xmm11,%%xmm9			\n\t	movaps	%%xmm15,%%xmm13		\n\t"/* xmm1 = cpy of y2 */\
			"movaps	0x80(%%rax),%%xmm0		\n\t	movaps	0x90(%%rax),%%xmm4		\n\t	movaps	0xa0(%%rax),%%xmm8		\n\t	movaps	0xb0(%%rax),%%xmm12	\n\t"/* Reload y1 [in flo2 slot] */\
			"mulpd	0x80(%%rdx),%%xmm0		\n\t	mulpd	0x90(%%rdx),%%xmm4		\n\t	mulpd	0xa0(%%rdx),%%xmm8		\n\t	mulpd	0xb0(%%rdx),%%xmm12	\n\t"/* y1 *= x2 */\
			"mulpd	0x40(%%rdx),%%xmm1		\n\t	mulpd	0x50(%%rdx),%%xmm5		\n\t	mulpd	0x60(%%rdx),%%xmm9		\n\t	mulpd	0x70(%%rdx),%%xmm13	\n\t"/* y2 *= x1 */\
			"mulpd	0x80(%%rdx),%%xmm3		\n\t	mulpd	0x90(%%rdx),%%xmm7		\n\t	mulpd	0xa0(%%rdx),%%xmm11		\n\t	mulpd	0xb0(%%rdx),%%xmm15	\n\t"/* y2 *= x2 */\
			"/* Digit 3: */\n\t"\
			"addpd	%%xmm2,%%xmm1			\n\t	addpd	%%xmm6,%%xmm5			\n\t	addpd	%%xmm10,%%xmm9			\n\t	addpd	%%xmm14,%%xmm13		\n\t"/* y2*x1 += CY */\
			"addpd	%%xmm1,%%xmm0			\n\t	addpd	%%xmm5,%%xmm4			\n\t	addpd	%%xmm9,%%xmm8			\n\t	addpd	%%xmm13,%%xmm12		\n\t"/* fprod3 = y1*x2 + y2*x1 + CY */\
			"roundpd $0,%%xmm0,%%xmm2		\n\t	roundpd $0,%%xmm4,%%xmm6		\n\t	roundpd $0,%%xmm8,%%xmm10		\n\t	roundpd $0,%%xmm12,%%xmm14	\n\t"\
			"subpd	%%xmm2,%%xmm0			\n\t	subpd	%%xmm6,%%xmm4			\n\t	subpd	%%xmm10,%%xmm8			\n\t	subpd	%%xmm14,%%xmm12		\n\t"/* fprod3 -= CY */\
			"mulpd	-0x10(%%rbx),%%xmm0		\n\t	mulpd	-0x10(%%rbx),%%xmm4		\n\t	mulpd	-0x10(%%rbx),%%xmm8		\n\t	mulpd	-0x10(%%rbx),%%xmm12\n\t"/* fprod3 *= two26f */\
			"mulpd	     (%%rbx),%%xmm2		\n\t	mulpd	     (%%rbx),%%xmm6		\n\t	mulpd	     (%%rbx),%%xmm10	\n\t	mulpd	     (%%rbx),%%xmm14\n\t"/* CY    *= two26i */\
			"/* Digit 4: */\n\t"\
			"addpd	%%xmm2,%%xmm3			\n\t	addpd	%%xmm6,%%xmm7			\n\t	addpd	%%xmm10,%%xmm11			\n\t	addpd	%%xmm14,%%xmm15		\n\t"/* CY = y2*x2 + CY */\
			"movaps	%%xmm3,%%xmm1			\n\t	movaps	%%xmm7,%%xmm5			\n\t	movaps	%%xmm11,%%xmm9			\n\t	movaps	%%xmm15,%%xmm13		\n\t"/* fprod4 = cpy of CY */\
			"roundpd $0,%%xmm3,%%xmm3		\n\t	roundpd $0,%%xmm7,%%xmm7		\n\t	roundpd $0,%%xmm11,%%xmm11		\n\t	roundpd $0,%%xmm15,%%xmm15	\n\t"\
			"subpd	%%xmm3,%%xmm1			\n\t	subpd	%%xmm7,%%xmm5			\n\t	subpd	%%xmm11,%%xmm9			\n\t	subpd	%%xmm15,%%xmm13		\n\t"/* fprod4 -= CY. Note we use the carry in the next section! */\
			"mulpd	-0x10(%%rbx),%%xmm1		\n\t	mulpd	-0x10(%%rbx),%%xmm5		\n\t	mulpd	-0x10(%%rbx),%%xmm9		\n\t	mulpd	-0x10(%%rbx),%%xmm13\n\t"/* fprod4 *= two26f */\
		/* If h < l, calculate h-l+q; otherwise h-l. Use leading 52 bits to approximate the full 78-bit compare. Result is in [0, q). */\
			"movaps	0x100(%%rax),%%xmm2		\n\t	movaps	0x110(%%rax),%%xmm6		\n\t	movaps	0x120(%%rax),%%xmm10	\n\t	movaps	0x130(%%rax),%%xmm14\n\t"/* fhi, top 52 bits */\
			"mulpd	-0x10(%%rbx),%%xmm3		\n\t	mulpd	-0x10(%%rbx),%%xmm7		\n\t	mulpd	-0x10(%%rbx),%%xmm11	\n\t	mulpd	-0x10(%%rbx),%%xmm15\n\t"/* flo2 *= two26f */\
			"addpd	%%xmm1,%%xmm3			\n\t	addpd	%%xmm5,%%xmm7			\n\t	addpd	%%xmm9,%%xmm11			\n\t	addpd	%%xmm13,%%xmm15		\n\t"/* flo, top 52 bits */\
			"movaps	0xc0(%%rax),%%xmm1		\n\t	movaps	0xd0(%%rax),%%xmm5		\n\t	movaps	0xe0(%%rax),%%xmm9		\n\t	movaps	0xf0(%%rax),%%xmm13	\n\t"/* fhi, low 26 bits */\
			"cmppd	$0x1,%%xmm3,%%xmm2		\n\t	cmppd	$0x1,%%xmm7,%%xmm6		\n\t	cmppd	$0x1,%%xmm11,%%xmm10	\n\t	cmppd	$0x1,%%xmm15,%%xmm14\n\t"/* bitmask = (fhi < flo)? */\
			"movaps	%%xmm2,    (%%rax)		\n\t	movaps	%%xmm6,0x10(%%rax)		\n\t	movaps	%%xmm10,0x20(%%rax)		\n\t	movaps	%%xmm14,0x30(%%rax)		\n\t"/* Store bitmask to free up a register */\
			"movaps	0x100(%%rax),%%xmm2		\n\t	movaps	0x110(%%rax),%%xmm6		\n\t	movaps	0x120(%%rax),%%xmm10	\n\t	movaps	0x130(%%rax),%%xmm14\n\t"/* fhi, top 52 bits */\
			"subpd	%%xmm0,%%xmm1			\n\t	subpd	%%xmm4,%%xmm5			\n\t	subpd	%%xmm8,%%xmm9			\n\t	subpd	%%xmm12,%%xmm13		\n\t"/* (fhi - flo), low 26 bits, xmm0 FREE */\
			"subpd	%%xmm3,%%xmm2			\n\t	subpd	%%xmm7,%%xmm6			\n\t	subpd	%%xmm11,%%xmm10			\n\t	subpd	%%xmm15,%%xmm14		\n\t"/* (fhi - flo), top 52 bits, xmm3 FREE */\
			"movaps	-0x10(%%rbx),%%xmm3		\n\t	movaps	-0x10(%%rbx),%%xmm7		\n\t	movaps	-0x10(%%rbx),%%xmm11	\n\t	movaps	-0x10(%%rbx),%%xmm15\n\t"/* 2^26 */\
			"mulpd	    (%%rdx),%%xmm3		\n\t	mulpd	0x10(%%rdx),%%xmm7		\n\t	mulpd	0x20(%%rdx),%%xmm11		\n\t	mulpd	0x30(%%rdx),%%xmm15	\n\t"/* fq, low 26 bits */\
			"movaps	0xc0(%%rdx),%%xmm0		\n\t	movaps	0xd0(%%rdx),%%xmm4		\n\t	movaps	0xe0(%%rdx),%%xmm8		\n\t	movaps	0xf0(%%rdx),%%xmm12	\n\t"/* qhi52 :+ fq, top 52 bits */\
			"andpd	(%%rax),%%xmm0			\n\t	andpd	0x10(%%rax),%%xmm4		\n\t	andpd	0x20(%%rax),%%xmm8		\n\t	andpd	0x30(%%rax),%%xmm12	\n\t"/* qhi52 & bitmask */\
			"addpd	%%xmm0,%%xmm2			\n\t	addpd	%%xmm4,%%xmm6			\n\t	addpd	%%xmm8,%%xmm10			\n\t	addpd	%%xmm12,%%xmm14		\n\t"/* xhi = (h-l)hi + (qhi52 & bitmask) */\
			"movaps	%%xmm3,0x40(%%rax)		\n\t	movaps	%%xmm7,0x50(%%rax)		\n\t	movaps	%%xmm11,0x60(%%rax)		\n\t	movaps	%%xmm15,0x70(%%rax)	\n\t"/* Store qlo26 to free up a register */\
			"andpd	(%%rax),%%xmm3			\n\t	andpd	0x10(%%rax),%%xmm7		\n\t	andpd	0x20(%%rax),%%xmm11		\n\t	andpd	0x30(%%rax),%%xmm15	\n\t"/* qlo26 &= bitmask */\
			"addpd	%%xmm3,%%xmm1			\n\t	addpd	%%xmm7,%%xmm5			\n\t	addpd	%%xmm11,%%xmm9			\n\t	addpd	%%xmm15,%%xmm13		\n\t"/* xlo = (h-l)lo + (qlo26 & bitmask) */\
			"/* xhi52,xlo26 in xmm2,xmm1; qhi52,qlo26 in 0xc0(%%rdx),0x40(%%rax): */\n\t"\
		/* if((pshift >> j) & (uint64)1) { */\
			"movslq	%[__pshift],%%rax	\n\t"\
			"movslq	%[__j],%%rcx		\n\t"\
			"shrq	%%cl,%%rax			\n\t"\
			"andq	$0x1,%%rax			\n\t"\
		/* Double: Use that high part < 2^52 (strictly >= 0): xhi52,xlo26 in xmm2,xmm1; qhi52,qlo26 in 0xc0(%%rdx),0x40(%%rax): */\
		"je twopmodq78_3wdq8_gcc64	\n\t"\
			"movq	%[__ax0],%%rax		\n\t"\
			"addpd	%%xmm2,%%xmm2		\n\t	addpd	%%xmm6,%%xmm6		\n\t	addpd	%%xmm10,%%xmm10		\n\t	addpd	%%xmm14,%%xmm14		\n\t"/* top 52 bits */\
			"addpd	%%xmm1,%%xmm1		\n\t	addpd	%%xmm5,%%xmm5		\n\t	addpd	%%xmm9,%%xmm9		\n\t	addpd	%%xmm13,%%xmm13		\n\t"/* low 26 bits */\
			"/* If x > q, subtract q: */\n\t"\
			"movaps	0xc0(%%rdx),%%xmm0	\n\t	movaps	0xd0(%%rdx),%%xmm4	\n\t	movaps	0xe0(%%rdx),%%xmm8	\n\t	movaps	0xf0(%%rdx),%%xmm12	\n\t"/* qhi52 */\
			"movaps	%%xmm0,%%xmm3		\n\t	movaps	%%xmm4,%%xmm7		\n\t	movaps	%%xmm8,%%xmm11		\n\t	movaps	%%xmm12,%%xmm15		\n\t"/* cpy of qhi */\
			"cmppd	$0x2,%%xmm2,%%xmm3	\n\t	cmppd	$0x2,%%xmm6,%%xmm7	\n\t	cmppd	$0x2,%%xmm10,%%xmm11\n\t	cmppd	$0x2,%%xmm14,%%xmm15\n\t"/* bitmask = (qhi <= xhi)? */\
			"andpd	%%xmm3,%%xmm0		\n\t	andpd	%%xmm7,%%xmm4		\n\t	andpd	%%xmm11,%%xmm8		\n\t	andpd	%%xmm15,%%xmm12		\n\t"/* qhi52 & bitmask */\
			"andpd	0x40(%%rax),%%xmm3	\n\t	andpd	0x50(%%rax),%%xmm7	\n\t	andpd	0x60(%%rax),%%xmm11	\n\t	andpd	0x70(%%rax),%%xmm15	\n\t"/* qlo26 & bitmask */\
			"subpd	%%xmm0,%%xmm2		\n\t	subpd	%%xmm4,%%xmm6		\n\t	subpd	%%xmm8,%%xmm10		\n\t	subpd	%%xmm12,%%xmm14		\n\t"/* x % q, top 52 bits */\
			"subpd	%%xmm3,%%xmm1		\n\t	subpd	%%xmm7,%%xmm5		\n\t	subpd	%%xmm11,%%xmm9		\n\t	subpd	%%xmm15,%%xmm13		\n\t"/* x % q, low 26 bits */\
		"twopmodq78_3wdq8_gcc64:	\n\t"\
		/* } */\
		/* Normalize the result: */\
			"movaps	     (%%rbx),%%xmm3		\n\t	movaps	-0x10(%%rbx),%%xmm7		\n\t"/* two26i,f */\
			"mulpd	%%xmm3,%%xmm1			\n\t	mulpd	%%xmm3,%%xmm5			\n\t	mulpd	%%xmm3,%%xmm9			\n\t	mulpd	%%xmm3,%%xmm13		\n\t"/* xlo *= two26i */\
			"mulpd	%%xmm3,%%xmm2			\n\t	mulpd	%%xmm3,%%xmm6			\n\t	mulpd	%%xmm3,%%xmm10			\n\t	mulpd	%%xmm3,%%xmm14		\n\t"/* xhi *= two26i */\
			"movaps	%%xmm1,%%xmm0			\n\t	movaps	%%xmm5,%%xmm4			\n\t	movaps	%%xmm9,%%xmm8			\n\t	movaps	%%xmm13,%%xmm12		\n\t"/* Init: fcy = cpy of ~xlo */\
			"roundpd $0,%%xmm1,%%xmm1		\n\t	roundpd $0,%%xmm5,%%xmm5		\n\t	roundpd $0,%%xmm9,%%xmm9		\n\t	roundpd $0,%%xmm13,%%xmm13	\n\t"/* fcy */\
			"subpd	%%xmm1,%%xmm0			\n\t	subpd	%%xmm5,%%xmm4			\n\t	subpd	%%xmm9,%%xmm8			\n\t	subpd	%%xmm13,%%xmm12		\n\t"/* fx0 -= fcy */\
			"mulpd	%%xmm7,%%xmm0			\n\t	mulpd	%%xmm7,%%xmm4			\n\t	mulpd	%%xmm7,%%xmm8			\n\t	mulpd	%%xmm7,%%xmm12		\n\t"/* fx0 *= two26f */\
			"mulpd	%%xmm3,%%xmm1			\n\t	mulpd	%%xmm3,%%xmm5			\n\t	mulpd	%%xmm3,%%xmm9			\n\t	mulpd	%%xmm3,%%xmm13		\n\t"/* fcy *= two26i */\
			"addpd	%%xmm1,%%xmm2			\n\t	addpd	%%xmm5,%%xmm6			\n\t	addpd	%%xmm9,%%xmm10			\n\t	addpd	%%xmm13,%%xmm14		\n\t"/* Add carry into xhi */\
			"movaps	%%xmm2,%%xmm1			\n\t	movaps	%%xmm6,%%xmm5			\n\t	movaps	%%xmm10,%%xmm9			\n\t	movaps	%%xmm14,%%xmm13		\n\t"/* Init: fcy = cpy of ~xhi */\
			"roundpd $0,%%xmm2,%%xmm2		\n\t	roundpd $0,%%xmm6,%%xmm6		\n\t	roundpd $0,%%xmm10,%%xmm10		\n\t	roundpd $0,%%xmm14,%%xmm14	\n\t"/* fx2 = fcy (no divide by 2^26) */\
			"subpd	%%xmm2,%%xmm1			\n\t	subpd	%%xmm6,%%xmm5			\n\t	subpd	%%xmm10,%%xmm9			\n\t	subpd	%%xmm14,%%xmm13		\n\t"/* fx1 -= fcy */\
			"mulpd	%%xmm7,%%xmm1			\n\t	mulpd	%%xmm7,%%xmm5			\n\t	mulpd	%%xmm7,%%xmm9			\n\t	mulpd	%%xmm7,%%xmm13		\n\t"/* fx1 *= two26f */\
			:					/* outputs: none */\
			: [__aq0]	 "m" (Xaq0)	/* All inputs from memory addresses here */\
			 ,[__aqinv0] "m" (Xaqinv0)	\
			 ,[__ax0]	 "m" (Xax0)		\
			 ,[__two26i] "m" (Xtwo26i)	\
			 ,[__pshift] "m" (Xpshift)	\
			 ,[__j]		 "m" (Xj)		\
			: "cc","memory","cl","rax","rbx","rcx","rdx","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
		);\
		}

		//=======================
		/* Branchless version of the conditional doubling in the above macro:
			"xorq	$0x1,%%rax			\n\t"\
			"movd	%%rax,%%ymm0		\n\t"\
			"pshufd	$0,%%ymm0,%%ymm0	\n\t"\
			"cvtdq2pd	%%ymm0,%%ymm0	\n\t"\
			"xorpd	%%ymm3,%%ymm3		\n\t"\
			"cmppd	$0x0,%%ymm3,%%ymm0	\n\t"\
			"movq	%[__ax0],%%rax			\n\t	movaps	%%ymm0,%%ymm4			\n\t	movaps	%%ymm0,%%ymm8			\n\t	movaps	%%ymm0,%%ymm12	\n\t"\
		// Double: Use that high part < 2^52 (strictly >= 0): xhi52,xlo26 in ymm2,ymm1; qhi52,qlo26 in 0x180(%%rdx),0x80(%%rax): //\
			"movaps	%%ymm1,%%ymm3			\n\t	movaps	%%ymm5,%%ymm7			\n\t	movaps	%%ymm9,%%ymm11			\n\t	movaps	%%ymm13,%%ymm15		\n\t"// cpy of xlo //\
			"andpd	%%ymm0,%%ymm3			\n\t	andpd	%%ymm4,%%ymm7			\n\t	andpd	%%ymm8,%%ymm11			\n\t	andpd	%%ymm12,%%ymm15		\n\t"// xlo_ & bitmask //\
			"andpd	%%ymm2,%%ymm0			\n\t	andpd	%%ymm6,%%ymm4			\n\t	andpd	%%ymm10,%%ymm8			\n\t	andpd	%%ymm14,%%ymm12		\n\t"// xhi_ & bitmask, overwrite bitmask with result //\
			"addpd	%%ymm3,%%ymm1			\n\t	addpd	%%ymm7,%%ymm5			\n\t	addpd	%%ymm11,%%ymm9			\n\t	addpd	%%ymm15,%%ymm13		\n\t"// low 26 bits //\
			"addpd	%%ymm0,%%ymm2			\n\t	addpd	%%ymm4,%%ymm6			\n\t	addpd	%%ymm8,%%ymm10			\n\t	addpd	%%ymm12,%%ymm14		\n\t"// top 52 bits //\
			"// If x > q, subtract q: //\n\t"\
			"movaps	0x180(%%rdx),%%ymm0		\n\t	movaps	0x1a0(%%rdx),%%ymm4		\n\t	movaps	0x1c0(%%rdx),%%ymm8		\n\t	movaps	0x1e0(%%rdx),%%ymm12	\n\t"// qhi52 //\
			"movaps	%%ymm0,%%ymm3			\n\t	movaps	%%ymm4,%%ymm7			\n\t	movaps	%%ymm8,%%ymm11			\n\t	movaps	%%ymm12,%%ymm15		\n\t"// cpy of qhi //\
			"cmppd	$0x2,%%ymm2,%%ymm3		\n\t	cmppd	$0x2,%%ymm6,%%ymm7		\n\t	cmppd	$0x2,%%ymm10,%%ymm11	\n\t	cmppd	$0x2,%%ymm14,%%ymm15\n\t"// bitmask = (qhi <= xhi)? //\
			"andpd	%%ymm3,%%ymm0			\n\t	andpd	%%ymm7,%%ymm4			\n\t	andpd	%%ymm11,%%ymm8			\n\t	andpd	%%ymm15,%%ymm12		\n\t"// qhi52 & bitmask //\
			"andpd	0x80(%%rax),%%ymm3		\n\t	andpd	0xa0(%%rax),%%ymm7		\n\t	andpd	0xc0(%%rax),%%ymm11		\n\t	andpd	0xe0(%%rax),%%ymm15	\n\t"// qlo26 & bitmask //\
			"subpd	%%ymm0,%%ymm2			\n\t	subpd	%%ymm4,%%ymm6			\n\t	subpd	%%ymm8,%%ymm10			\n\t	subpd	%%ymm12,%%ymm14		\n\t"// x % q, top 52 bits //\
			"subpd	%%ymm3,%%ymm1			\n\t	subpd	%%ymm7,%%ymm5			\n\t	subpd	%%ymm11,%%ymm9			\n\t	subpd	%%ymm15,%%ymm13		\n\t"// x % q, low 26 bits //\
		*/
		//========================


  #endif // AVX / SSE2

	// Select between add/sub-round-const for DNINT (works for all SSE2 versions) and vroundpd (SSE4+ only):
	#ifdef USE_ROUNDPD

		#define SSE2_twopmodq78_modmul_q4(Xfq0,Xpshift,Xj)\
		{\
		__asm__ volatile (\
			"movq	%[__fq0],%%rdx		\n\t"\
		"/* SQR_LOHI78_3WORD_DOUBLE_LEAVE_HIGH52_UNNORMALIZED_q4(): */\n\t"\
			"movaps	0x200(%%rdx),%%xmm6	/* two13i */	\n\t"\
			"/* fx0,1,2 assumed in xmm0,2,4 on loop entry */\n\t"\
			"mulpd	%%xmm6,%%xmm0			\n\t	mulpd	%%xmm6 ,%%xmm8			\n\t"\
			"mulpd	%%xmm6,%%xmm2			\n\t	mulpd	%%xmm6 ,%%xmm10			\n\t"\
			"mulpd	%%xmm6,%%xmm4			\n\t	mulpd	%%xmm6 ,%%xmm12			\n\t"\
			"movaps	%%xmm0,%%xmm1			\n\t	movaps	%%xmm8 ,%%xmm9			\n\t"\
			"movaps	%%xmm2,%%xmm3			\n\t	movaps	%%xmm10,%%xmm11			\n\t"\
			"addpd	%%xmm1,%%xmm1			\n\t	addpd	%%xmm9 ,%%xmm9			\n\t"\
			"addpd	%%xmm3,%%xmm3			\n\t	addpd	%%xmm11,%%xmm11			\n\t"\
			"movaps	%%xmm1,%%xmm5			\n\t	movaps	%%xmm9 ,%%xmm13			\n\t"\
			"mulpd	%%xmm0,%%xmm0			\n\t	mulpd	%%xmm8 ,%%xmm8			\n\t"\
			"mulpd	%%xmm2,%%xmm1			\n\t	mulpd	%%xmm10,%%xmm9			\n\t"\
			"mulpd	%%xmm2,%%xmm2			\n\t	mulpd	%%xmm10,%%xmm10			\n\t"\
			"mulpd	%%xmm4,%%xmm5			\n\t	mulpd	%%xmm12,%%xmm13			\n\t"\
			"mulpd	%%xmm4,%%xmm3			\n\t	mulpd	%%xmm12,%%xmm11			\n\t"\
			"mulpd	%%xmm4,%%xmm4			\n\t	mulpd	%%xmm12,%%xmm12			\n\t"\
			"/* Move this part of Digit 2 computation here to free up xmm5,13: */	\n\t"\
			"addpd	%%xmm5,%%xmm2			\n\t	addpd	%%xmm13,%%xmm10			\n\t"\
			"/* Digit 0: */					\n\t"\
			"movaps 0x220(%%rdx),%%xmm15	\n\t	movaps	0x210(%%rdx),%%xmm5		/* xmm15,5 = two26i,f */\n\t"\
			"roundpd	$0,%%xmm0,%%xmm6	\n\t	roundpd	$0,%%xmm8,%%xmm14			\n\t"\
			"subpd	%%xmm6,%%xmm0			\n\t	subpd	%%xmm14,%%xmm8			\n\t"\
			"mulpd	%%xmm15,%%xmm6			\n\t	mulpd	%%xmm15,%%xmm14	\n\t"\
			"mulpd	%%xmm5 ,%%xmm0			\n\t	mulpd	%%xmm5 ,%%xmm8		\n\t"\
			"/* Digit 1: */					\n\t"\
			"addpd	%%xmm6,%%xmm1			\n\t	addpd	%%xmm14,%%xmm9			\n\t"\
			"roundpd	$0,%%xmm1,%%xmm6	\n\t	roundpd	$0,%%xmm9,%%xmm14			\n\t"\
			"subpd	%%xmm6,%%xmm1			\n\t	subpd	%%xmm14,%%xmm9			\n\t"\
			"mulpd	%%xmm15,%%xmm6			\n\t	mulpd	%%xmm15,%%xmm14	\n\t"\
			"mulpd	%%xmm5 ,%%xmm1			\n\t	mulpd	%%xmm5 ,%%xmm9		\n\t"\
			"/* Digit 2: Require both hi and lo half of output to be nonnegative, so leave unbalanced: */	\n\t"\
			"addpd	%%xmm6,%%xmm2			\n\t	addpd	%%xmm14,%%xmm10			\n\t"\
			"movaps	%%xmm2,%%xmm6			\n\t	movaps	%%xmm10,%%xmm14			\n\t"\
			"subpd	0x240(%%rdx),%%xmm6		\n\t	subpd	0x240(%%rdx),%%xmm14	\n\t"\
			"roundpd	$0,%%xmm6,%%xmm6	\n\t	roundpd	$0,%%xmm14,%%xmm14			\n\t"\
			"subpd	%%xmm6,%%xmm2			\n\t	subpd	%%xmm14,%%xmm10			\n\t"\
			"mulpd	%%xmm15,%%xmm6			\n\t	mulpd	%%xmm15,%%xmm14			\n\t"\
			"mulpd	%%xmm5 ,%%xmm2			\n\t	mulpd	%%xmm5 ,%%xmm10			\n\t"\
			"/* Digit 3: */					\n\t"\
			"addpd	%%xmm6,%%xmm3			\n\t	addpd	%%xmm14,%%xmm11			\n\t"\
			"roundpd	$0,%%xmm3,%%xmm6	\n\t	roundpd	$0,%%xmm11,%%xmm14			\n\t"\
			"subpd	%%xmm6,%%xmm3			\n\t	subpd	%%xmm14,%%xmm11			\n\t"\
			"mulpd	%%xmm15,%%xmm6			\n\t	mulpd	%%xmm15,%%xmm14			\n\t"\
			"mulpd	%%xmm5 ,%%xmm3			\n\t	mulpd	%%xmm5 ,%%xmm11			\n\t"\
			"/* Digits 4+5: */				\n\t"\
			"addpd	%%xmm6,%%xmm4			\n\t	addpd	%%xmm14,%%xmm12			\n\t"\
			"mulpd	%%xmm5,%%xmm4			\n\t	mulpd	%%xmm5 ,%%xmm12			\n\t"\
			"/* No carry in this variant. flo0,1,2 in xmm0,1,2; fhi0,[1,2] in xmm3,4 */	\n\t"\
			"movaps	%%xmm3,0x0e0(%%rdx)		\n\t	movaps	%%xmm11,0x0f0(%%rdx)	\n\t"\
			"movaps	%%xmm4,0x100(%%rdx)		\n\t	movaps	%%xmm12,0x110(%%rdx)	\n\t"\
		"/* flo = MULL78_3WORD_DOUBLE_q2(flo, fqinv): */						\n\t"\
			"/* Digit 0: */					\n\t"\
			"movaps	%%xmm0,%%xmm3			\n\t	movaps	%%xmm8 ,%%xmm11			\n\t"\
			"movaps	%%xmm0,%%xmm4			\n\t	movaps	%%xmm8 ,%%xmm12			\n\t"\
			"mulpd	0x80(%%rdx),%%xmm0		\n\t	mulpd	0x90(%%rdx),%%xmm8		\n\t"\
			"roundpd	$0,%%xmm0,%%xmm6	\n\t	roundpd	$0,%%xmm8,%%xmm14			\n\t"\
			"subpd	%%xmm6,%%xmm0			\n\t	subpd	%%xmm14,%%xmm8			\n\t"\
			"mulpd	%%xmm15,%%xmm6			\n\t	mulpd	%%xmm15,%%xmm14			\n\t"\
			"mulpd	%%xmm5 ,%%xmm0			\n\t	mulpd	%%xmm5 ,%%xmm8			\n\t"\
			"/* Digit 1: */					\n\t"\
			"mulpd	0xa0(%%rdx),%%xmm3		\n\t	mulpd	0xb0(%%rdx),%%xmm11		\n\t"\
			"addpd	%%xmm6,%%xmm3			\n\t	addpd	%%xmm14,%%xmm11			\n\t"\
			"movaps	%%xmm1,%%xmm6			\n\t	movaps	%%xmm9 ,%%xmm14			\n\t"\
			"mulpd	0x80(%%rdx),%%xmm1		\n\t	mulpd	0x90(%%rdx),%%xmm9		\n\t"\
			"addpd	%%xmm3,%%xmm1			\n\t	addpd	%%xmm11,%%xmm9			\n\t"\
			"roundpd	$0,%%xmm1,%%xmm3	\n\t	roundpd	$0,%%xmm9,%%xmm11			\n\t"\
			"subpd	%%xmm3,%%xmm1			\n\t	subpd	%%xmm11,%%xmm9			\n\t"\
			"mulpd	%%xmm5 ,%%xmm1			\n\t	mulpd	%%xmm5 ,%%xmm9			\n\t"\
			"mulpd	%%xmm15,%%xmm3			\n\t	mulpd	%%xmm15,%%xmm11			\n\t"\
			"/* Digit 2: Require output to be nonnegative, so leave MSW unbalanced: */	\n\t"\
			"mulpd	0x80(%%rdx),%%xmm2		\n\t	mulpd	0x90(%%rdx),%%xmm10		\n\t"\
			"mulpd	0xa0(%%rdx),%%xmm6		\n\t	mulpd	0xb0(%%rdx),%%xmm14		\n\t"\
			"mulpd	0xc0(%%rdx),%%xmm4		\n\t	mulpd	0xd0(%%rdx),%%xmm12		\n\t"\
			"addpd	%%xmm6,%%xmm2			\n\t	addpd	%%xmm14,%%xmm10			\n\t"\
			"addpd	%%xmm4,%%xmm3			\n\t	addpd	%%xmm12,%%xmm11			\n\t"\
			"addpd	%%xmm3,%%xmm2			\n\t	addpd	%%xmm11,%%xmm10			\n\t"\
			"movaps	%%xmm2,%%xmm3			\n\t	movaps	%%xmm10,%%xmm11			\n\t"\
			"subpd	0x240(%%rdx),%%xmm3		\n\t	subpd	0x240(%%rdx),%%xmm11	\n\t"\
			"roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd	$0,%%xmm11,%%xmm11			\n\t"\
			"subpd	%%xmm3,%%xmm2			\n\t	subpd	%%xmm11,%%xmm10			\n\t"\
			"mulpd	%%xmm5,%%xmm2			\n\t	mulpd	%%xmm5 ,%%xmm10			\n\t"\
			"/* MULH96(q,lo,lo) --> lo = (q*lo)/2^78 */								\n\t"\
			"/*** movq	%[__fq0],%%rdx		Moved to start of asm-block ***/\n\t"\
			"/* Digit 0: */					\n\t"\
			"movaps	%%xmm0,%%xmm3			\n\t	movaps	%%xmm8 ,%%xmm11			\n\t"/* xmm3 = cpy of y0 */\
			"movaps	%%xmm0,%%xmm4			\n\t	movaps	%%xmm8 ,%%xmm12			\n\t"/* xmm4 = cpy of y0 */\
			"mulpd	   (%%rdx),%%xmm0		\n\t	mulpd	0x10(%%rdx),%%xmm8		\n\t"/* fprod0 = y0*x0; fq0 in rdx */\
			"mulpd	   %%xmm15,%%xmm0		\n\t	mulpd	   %%xmm15,%%xmm8		\n\t"/* fcy *= two26i */\
			"/* Digit 1: */					\n\t"\
			"movaps	%%xmm1,%%xmm6			\n\t	movaps	%%xmm9 ,%%xmm14			\n\t"/* xmm6 = cpy of y1 */\
			"mulpd	0x20(%%rdx),%%xmm3		\n\t	mulpd	0x30(%%rdx),%%xmm11		\n\t"/* y0*x1 */\
			"addpd	%%xmm0,%%xmm3			\n\t	addpd	%%xmm8 ,%%xmm11			\n\t"/* y0*x1 + fcy; xmm0 FREE */\
			"mulpd	   (%%rdx),%%xmm1		\n\t	mulpd	0x10(%%rdx),%%xmm9		\n\t"/* y1*x0 */\
			"addpd	%%xmm3,%%xmm1			\n\t	addpd	%%xmm11,%%xmm9			\n\t"/* fprod1 = x0*y1 + x1*y0; xmm3 FREE */\
			"mulpd	   %%xmm15,%%xmm1		\n\t	mulpd	   %%xmm15,%%xmm9		\n\t"/* fcy *= two26i */\
			"/* Digit 2: Require output to be nonnegative, so leave MSW unbalanced: */	\n\t"\
			"movaps	%%xmm2,%%xmm3			\n\t	movaps	%%xmm10,%%xmm11			\n\t"/* xmm3 = cpy of y2 */\
			"movaps	%%xmm6,%%xmm0			\n\t	movaps	%%xmm14,%%xmm8			\n\t"/* xmm0 = cpy of y1 */\
			"mulpd	   (%%rdx),%%xmm2		\n\t	mulpd	0x10(%%rdx),%%xmm10		\n\t"/* y2*x0 */\
			"mulpd	0x20(%%rdx),%%xmm6		\n\t	mulpd	0x30(%%rdx),%%xmm14		\n\t"/* y1*x1 */\
			"mulpd	0x40(%%rdx),%%xmm4		\n\t	mulpd	0x50(%%rdx),%%xmm12		\n\t"/* y0*x2 */\
			"addpd	%%xmm6,%%xmm2			\n\t	addpd	%%xmm14,%%xmm10			\n\t"/* y1*x1 + y2*x0; xmm6 FREE */\
			"addpd	%%xmm4,%%xmm1			\n\t	addpd	%%xmm12,%%xmm9			\n\t"/* y0*x2 + fcy; xmm4 FREE */\
			"addpd	%%xmm1,%%xmm2			\n\t	addpd	%%xmm9 ,%%xmm10			\n\t"/* fprod2; xmm1 FREE */\
			"subpd	0x240(%%rdx),%%xmm2		\n\t	subpd	0x240(%%rdx),%%xmm10	\n\t"/* fprod2 - 0.5 */\
			"roundpd	$0,%%xmm2,%%xmm2	\n\t	roundpd	$0,%%xmm10,%%xmm10		\n\t"/* fcy */\
			"mulpd	   %%xmm15,%%xmm2		\n\t	mulpd	   %%xmm15,%%xmm10		\n\t"/* fcy *= two26i */\
			"/* Precompute all the needed partial products: */						\n\t"\
			"movaps	%%xmm3,%%xmm1			\n\t	movaps	%%xmm11,%%xmm9			\n\t"/* xmm1 = cpy of y2 */\
			"mulpd	0x40(%%rdx),%%xmm0		\n\t	mulpd	0x50(%%rdx),%%xmm8		\n\t"/* y1*x2 */\
			"mulpd	0x20(%%rdx),%%xmm1		\n\t	mulpd	0x30(%%rdx),%%xmm9		\n\t"/* y2*x1 */\
			"mulpd	0x40(%%rdx),%%xmm3		\n\t	mulpd	0x50(%%rdx),%%xmm11		\n\t"/* y2*x2 */\
			"/* Digit 3: */					\n\t"\
			"addpd	%%xmm2,%%xmm1			\n\t	addpd	%%xmm10,%%xmm9			\n\t"/* fy2*fx1 + fcy */\
			"addpd	%%xmm1,%%xmm0			\n\t	addpd	%%xmm9 ,%%xmm8			\n\t"/* fprod3 = fy1*fx2 + fy2*fx1 + fcy */\
			"roundpd	$0,%%xmm0,%%xmm6	\n\t	roundpd	$0,%%xmm8,%%xmm14		\n\t"/* fcy */\
			"subpd	%%xmm6,%%xmm0			\n\t	subpd	%%xmm14,%%xmm8			\n\t"/* fprod3 -= fcy */\
			"mulpd	   %%xmm15,%%xmm6		\n\t	mulpd	   %%xmm15,%%xmm14		\n\t"/* fcy    *= two26i */\
			"mulpd	%%xmm5,%%xmm0			\n\t	mulpd	%%xmm5 ,%%xmm8			\n\t"/* fprod3 *= two26f */\
			"/* Digit 4: */					\n\t"\
			"addpd	%%xmm6,%%xmm3			\n\t	addpd	%%xmm14,%%xmm11			\n\t"/* fcy = fy2*fx2 + fcy */\
			"movaps	%%xmm3,%%xmm1			\n\t	movaps	%%xmm11,%%xmm9			\n\t"/* fprod4 = cpy of fcy */\
			"roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd	$0,%%xmm11,%%xmm11		\n\t"/* fcy */\
			"subpd	%%xmm3,%%xmm1			\n\t	subpd	%%xmm11,%%xmm9			\n\t"/* fprod4 -= fcy */\
			"mulpd	%%xmm5,%%xmm1			\n\t	mulpd	%%xmm5 ,%%xmm9			\n\t"/* fprod4 *= two26f */\
			"/* If h < l, calculate h-l+q; otherwise h-l. Use leading 52 bits to approximate the full 78-bit compare. Result is in [0, q). */\n\t"\
			"movaps %%xmm5,%%xmm13	/* Need a copy of two26f, stored in %%xmm5 */\n\t"\
			"movaps	0x100(%%rdx),%%xmm2		\n\t	movaps	0x110(%%rdx),%%xmm10	\n\t"\
			"movaps	%%xmm2,%%xmm6			\n\t	movaps	%%xmm10,%%xmm14			\n\t"\
			"movaps	0x60(%%rdx),%%xmm4		\n\t	movaps	0x70(%%rdx),%%xmm12		\n\t"\
			"mulpd	%%xmm5,%%xmm3			\n\t	mulpd	%%xmm5 ,%%xmm11			\n\t"\
			"mulpd	   (%%rdx),%%xmm5		\n\t	mulpd	0x10(%%rdx),%%xmm13		\n\t"\
			"addpd	%%xmm1,%%xmm3			\n\t	addpd	%%xmm9 ,%%xmm11			\n\t"\
			"movaps	0x0e0(%%rdx),%%xmm1		\n\t	movaps	0x0f0(%%rdx),%%xmm9		\n\t"\
			"cmppd	$0x1,%%xmm3,%%xmm6		\n\t	cmppd	$0x1,%%xmm11,%%xmm14	\n\t"\
			"subpd	%%xmm0,%%xmm1			\n\t	subpd	%%xmm8 ,%%xmm9			\n\t"\
			"subpd	%%xmm3,%%xmm2			\n\t	subpd	%%xmm11,%%xmm10			\n\t"\
			"movaps	%%xmm4,%%xmm0			\n\t	movaps	%%xmm12,%%xmm8			\n\t"\
			"andpd	%%xmm6,%%xmm4			\n\t	andpd	%%xmm14,%%xmm12			\n\t"\
			"addpd	%%xmm4,%%xmm2			\n\t	addpd	%%xmm12,%%xmm10			\n\t"\
			"andpd	%%xmm5,%%xmm6			\n\t	andpd	%%xmm13,%%xmm14			\n\t"\
			"addpd	%%xmm6,%%xmm1			\n\t	addpd	%%xmm14,%%xmm9			\n\t"\
			"/* qlo26 in xmm5, qhi52 in xmm0 */		\n\t"\
			"/* if((pshift >> j) & (uint64)1) { */	\n\t"\
			"movslq	%[__pshift],%%rax		\n\t"\
			"movslq	%[__j],%%rcx			\n\t"\
			"shrq	%%cl,%%rax				\n\t"\
			"andq	$0x1,%%rax				\n\t"\
		"je	twopmodq78_3wdq4			\n\t"\
		"/* Double: Use that high part < 2^52 (strictly >= 0): xhi52,xlo26 in xmm2,xmm1; qhi52,qlo26 in xmm0,xmm5 */\n\t"\
		"	movaps		%%xmm0,%%xmm6		\n\t	movaps		%%xmm8 ,%%xmm14	/* cpy of qhi */			\n\t"\
		"	addpd		%%xmm2,%%xmm2		\n\t	addpd		%%xmm10,%%xmm10	/* top 52 bits */			\n\t"\
		"	addpd		%%xmm1,%%xmm1		\n\t	addpd		%%xmm9 ,%%xmm9	/* low 26 bits */			\n\t"\
		"/* If x > q, subtract q: */		\n\t"\
		"	cmppd	$0x2,%%xmm2,%%xmm6		\n\t	cmppd	$0x2,%%xmm10,%%xmm14	/* bitmask = (qhi <= xhi) */\n\t"\
		"	andpd		%%xmm6,%%xmm0		\n\t	andpd		%%xmm14,%%xmm8	/* qhi52 & bitmask */		\n\t"\
		"	andpd		%%xmm6,%%xmm5		\n\t	andpd		%%xmm14,%%xmm13	/* qlo26 & bitmask */		\n\t"\
		"	subpd		%%xmm0,%%xmm2		\n\t	subpd		%%xmm8 ,%%xmm10	/* x mod q, top 52 bits */	\n\t"\
		"	subpd		%%xmm5,%%xmm1		\n\t	subpd		%%xmm13,%%xmm9	/* x mod q, low 26 bits */	\n\t"\
		"twopmodq78_3wdq4:					\n\t"\
			"/* } */						\n\t"\
			"/* Normalize the result: */	\n\t"\
			"movaps	0x210(%%rdx),%%xmm4		\n\t"\
			"movaps		%%xmm15,%%xmm3		\n\t"\
			"mulpd	%%xmm3,%%xmm1			\n\t	mulpd	%%xmm3 ,%%xmm9			\n\t"\
			"mulpd	%%xmm3,%%xmm2			\n\t	mulpd	%%xmm3 ,%%xmm10			\n\t"\
			"movaps	%%xmm1,%%xmm0			\n\t	movaps	%%xmm9 ,%%xmm8			\n\t"\
			"roundpd	$0,%%xmm1,%%xmm1	\n\t	roundpd	$0,%%xmm9 ,%%xmm9			\n\t"\
			"subpd	%%xmm1,%%xmm0			\n\t	subpd	%%xmm9 ,%%xmm8			\n\t"\
			"mulpd	%%xmm4,%%xmm0			\n\t	mulpd	%%xmm4 ,%%xmm8			\n\t"\
			"mulpd	%%xmm3,%%xmm1			\n\t	mulpd	%%xmm3 ,%%xmm9			\n\t"\
			"addpd	%%xmm1,%%xmm2			\n\t	addpd	%%xmm9 ,%%xmm10			\n\t"\
			"movaps	%%xmm2,%%xmm1			\n\t	movaps	%%xmm10,%%xmm9			\n\t"\
			"roundpd	$0,%%xmm2,%%xmm2	\n\t	roundpd	$0,%%xmm10,%%xmm10			\n\t"\
			"subpd	%%xmm2,%%xmm1			\n\t	subpd	%%xmm10,%%xmm9			\n\t"\
			"mulpd	%%xmm4,%%xmm1			\n\t	mulpd	%%xmm4 ,%%xmm9			\n\t"\
			"/* Move high 2 words of result into input registers expected by start of loop body: */\n\t"\
			"movaps	%%xmm2,%%xmm4 /* fx2 */	\n\t	movaps	%%xmm10,%%xmm12	/* hx2 */	\n\t"\
			"movaps	%%xmm1,%%xmm2 /* fx1 */	\n\t	movaps	%%xmm9 ,%%xmm10	/* hx1 */	\n\t"\
			:					/* outputs: none */\
			: [__fq0] "m" (Xfq0)	/* All inputs from memory addresses here */\
			,[__pshift] "m" (Xpshift)	\
			,[__j]		"m" (Xj)		\
			: "cl","rax","rcx","rdx","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
		);\
		}

	#else	// ifndef(USE_ROUNDPD)

		#define SSE2_twopmodq78_modmul_q4(Xfq0,Xpshift,Xj)\
		{\
		__asm__ volatile (\
			"movq	%[__fq0],%%rdx		\n\t"\
			"/* SQR_LOHI78_3WORD_DOUBLE_q4(): */\n\t"\
			"movaps	0x200(%%rdx),%%xmm6	/* two13i */	\n\t"\
			"/* fx0,1,2 assumed in xmm0,2,4 on loop entry */\n\t"\
			"mulpd	%%xmm6,%%xmm0			\n\t	mulpd	%%xmm6 ,%%xmm8			\n\t"\
			"mulpd	%%xmm6,%%xmm2			\n\t	mulpd	%%xmm6 ,%%xmm10			\n\t"\
			"mulpd	%%xmm6,%%xmm4			\n\t	mulpd	%%xmm6 ,%%xmm12			\n\t"\
			"movaps	0x230(%%rdx),%%xmm7		/* xmm7 = rnd_const shared between both columns*/\n\t"\
			"movaps	%%xmm0,%%xmm1			\n\t	movaps	%%xmm8 ,%%xmm9			\n\t"\
			"movaps	%%xmm2,%%xmm3			\n\t	movaps	%%xmm10,%%xmm11			\n\t"\
			"addpd	%%xmm1,%%xmm1			\n\t	addpd	%%xmm9 ,%%xmm9			\n\t"\
			"addpd	%%xmm3,%%xmm3			\n\t	addpd	%%xmm11,%%xmm11			\n\t"\
			"movaps	%%xmm1,%%xmm5			\n\t	movaps	%%xmm9 ,%%xmm13			\n\t"\
			"mulpd	%%xmm0,%%xmm0			\n\t	mulpd	%%xmm8 ,%%xmm8			\n\t"\
			"mulpd	%%xmm2,%%xmm1			\n\t	mulpd	%%xmm10,%%xmm9			\n\t"\
			"mulpd	%%xmm2,%%xmm2			\n\t	mulpd	%%xmm10,%%xmm10			\n\t"\
			"mulpd	%%xmm4,%%xmm5			\n\t	mulpd	%%xmm12,%%xmm13			\n\t"\
			"mulpd	%%xmm4,%%xmm3			\n\t	mulpd	%%xmm12,%%xmm11			\n\t"\
			"mulpd	%%xmm4,%%xmm4			\n\t	mulpd	%%xmm12,%%xmm12			\n\t"\
			"/* Move this part of Digit 2 computation here to free up xmm5,13: */	\n\t"\
			"addpd	%%xmm5,%%xmm2			\n\t	addpd	%%xmm13,%%xmm10			\n\t"\
			"/* Digit 0: */					\n\t"\
			"movaps 0x220(%%rdx),%%xmm15	\n\t	movaps	0x210(%%rdx),%%xmm5		/* xmm15,5 = two26i,f */\n\t"\
			"movaps	%%xmm0,%%xmm6			\n\t	movaps	%%xmm8 ,%%xmm14			\n\t"\
			"addpd	%%xmm7,%%xmm6			\n\t	addpd	%%xmm7 ,%%xmm14			\n\t"\
			"subpd	%%xmm7,%%xmm6			\n\t	subpd	%%xmm7 ,%%xmm14			\n\t"\
			"subpd	%%xmm6,%%xmm0			\n\t	subpd	%%xmm14,%%xmm8			\n\t"\
			"mulpd	%%xmm5,%%xmm0		\n\t	mulpd	%%xmm5,%%xmm8		\n\t"\
			"mulpd		%%xmm15,%%xmm6		\n\t	mulpd		%%xmm15,%%xmm14	\n\t"\
			"/* Digit 1: */					\n\t"\
			"addpd	%%xmm6,%%xmm1			\n\t	addpd	%%xmm14,%%xmm9			\n\t"\
			"movaps	%%xmm1,%%xmm6			\n\t	movaps	%%xmm9 ,%%xmm14			\n\t"\
			"addpd	%%xmm7,%%xmm6			\n\t	addpd	%%xmm7 ,%%xmm14			\n\t"\
			"subpd	%%xmm7,%%xmm6			\n\t	subpd	%%xmm7 ,%%xmm14			\n\t"\
			"subpd	%%xmm6,%%xmm1			\n\t	subpd	%%xmm14,%%xmm9			\n\t"\
			"mulpd	%%xmm5,%%xmm1		\n\t	mulpd	%%xmm5,%%xmm9		\n\t"\
			"mulpd		%%xmm15,%%xmm6		\n\t	mulpd		%%xmm15,%%xmm14	\n\t"\
			"/* Digit 2: Require both hi and lo half of output to be nonnegative, so leave unbalanced: */	\n\t"\
			"addpd	%%xmm6,%%xmm2			\n\t	addpd	%%xmm14,%%xmm10			\n\t"\
			"movaps	%%xmm2,%%xmm6			\n\t	movaps	%%xmm10,%%xmm14			\n\t"\
			"subpd	0x240(%%rdx),%%xmm6		\n\t	subpd	0x240(%%rdx),%%xmm14	\n\t"\
			"addpd	%%xmm7,%%xmm6			\n\t	addpd	%%xmm7 ,%%xmm14			\n\t"\
			"subpd	%%xmm7,%%xmm6			\n\t	subpd	%%xmm7 ,%%xmm14			\n\t"\
			"subpd	%%xmm6,%%xmm2			\n\t	subpd	%%xmm14,%%xmm10			\n\t"\
			"mulpd	%%xmm5,%%xmm2			\n\t	mulpd	%%xmm5 ,%%xmm10			\n\t"\
			"mulpd	%%xmm15,%%xmm6			\n\t	mulpd	%%xmm15,%%xmm14			\n\t"\
			"/* Digit 3: */					\n\t"\
			"addpd	%%xmm6,%%xmm3			\n\t	addpd	%%xmm14,%%xmm11			\n\t"\
			"movaps	%%xmm3,%%xmm6			\n\t	movaps	%%xmm11,%%xmm14			\n\t"\
			"addpd	%%xmm7,%%xmm6			\n\t	addpd	%%xmm7 ,%%xmm14			\n\t"\
			"subpd	%%xmm7,%%xmm6			\n\t	subpd	%%xmm7 ,%%xmm14			\n\t"\
			"subpd	%%xmm6,%%xmm3			\n\t	subpd	%%xmm14,%%xmm11			\n\t"\
			"mulpd	%%xmm5,%%xmm3			\n\t	mulpd	%%xmm5 ,%%xmm11			\n\t"\
			"mulpd	%%xmm15,%%xmm6			\n\t	mulpd	%%xmm15,%%xmm14			\n\t"\
			"/* Digit 4: */					\n\t"\
			"addpd	%%xmm6,%%xmm4			\n\t	addpd	%%xmm14,%%xmm12			\n\t"\
			"movaps	%%xmm4,%%xmm6			\n\t	movaps	%%xmm12,%%xmm14			\n\t"\
			"addpd	%%xmm7,%%xmm6			\n\t	addpd	%%xmm7 ,%%xmm14			\n\t"\
			"subpd	%%xmm7,%%xmm6			\n\t	subpd	%%xmm7 ,%%xmm14			\n\t"\
			"subpd	%%xmm6,%%xmm4			\n\t	subpd	%%xmm14,%%xmm12			\n\t"\
			"mulpd	%%xmm5,%%xmm4			\n\t	mulpd	%%xmm5 ,%%xmm12			\n\t"\
			"/* Digit 5 = the carry. flo0,1,2 in xmm0,1,2; fhi0,1,2 in xmm3,4,6 */	\n\t"\
			"mulpd	%%xmm5,%%xmm6			\n\t	mulpd	%%xmm5 ,%%xmm14			\n\t"\
			"addpd	%%xmm6,%%xmm4			\n\t	addpd	%%xmm14,%%xmm12			\n\t"\
			"movaps	%%xmm3,0x1a0(%%rdx)		\n\t	movaps	%%xmm11,0x1b0(%%rdx)	\n\t"\
			"movaps	%%xmm4,0x1c0(%%rdx)		\n\t	movaps	%%xmm12,0x1d0(%%rdx)	\n\t"\
			"/* flo = MULL78_3WORD_DOUBLE_q2(flo, fqinv): */						\n\t"\
			"/* Digit 0: */					\n\t"\
			"movaps	%%xmm0,%%xmm3			\n\t	movaps	%%xmm8 ,%%xmm11			\n\t"\
			"movaps	%%xmm0,%%xmm4			\n\t	movaps	%%xmm8 ,%%xmm12			\n\t"\
			"mulpd	0x80(%%rdx),%%xmm0		\n\t	mulpd	0x90(%%rdx),%%xmm8		\n\t"\
			"movaps	%%xmm0,%%xmm6			\n\t	movaps	%%xmm8 ,%%xmm14			\n\t"\
			"addpd	%%xmm7,%%xmm6			\n\t	addpd	%%xmm7 ,%%xmm14			\n\t"\
			"subpd	%%xmm7,%%xmm6			\n\t	subpd	%%xmm7 ,%%xmm14			\n\t"\
			"subpd	%%xmm6,%%xmm0			\n\t	subpd	%%xmm14,%%xmm8			\n\t"\
			"mulpd	%%xmm5,%%xmm0			\n\t	mulpd	%%xmm5 ,%%xmm8			\n\t"\
			"mulpd	%%xmm15,%%xmm6			\n\t	mulpd	%%xmm15,%%xmm14			\n\t"\
			"/* Digit 1: */					\n\t"\
			"mulpd	0xa0(%%rdx),%%xmm3		\n\t	mulpd	0xb0(%%rdx),%%xmm11		\n\t"\
			"addpd	%%xmm6,%%xmm3			\n\t	addpd	%%xmm14,%%xmm11			\n\t"\
			"movaps	%%xmm1,%%xmm6			\n\t	movaps	%%xmm9 ,%%xmm14			\n\t"\
			"mulpd	0x80(%%rdx),%%xmm1		\n\t	mulpd	0x90(%%rdx),%%xmm9		\n\t"\
			"addpd	%%xmm3,%%xmm1			\n\t	addpd	%%xmm11,%%xmm9			\n\t"\
			"movaps	%%xmm1,%%xmm3			\n\t	movaps	%%xmm9 ,%%xmm11			\n\t"\
			"addpd	%%xmm7,%%xmm3			\n\t	addpd	%%xmm7 ,%%xmm11			\n\t"\
			"subpd	%%xmm7,%%xmm3			\n\t	subpd	%%xmm7 ,%%xmm11			\n\t"\
			"subpd	%%xmm3,%%xmm1			\n\t	subpd	%%xmm11,%%xmm9			\n\t"\
			"mulpd	%%xmm5,%%xmm1			\n\t	mulpd	%%xmm5 ,%%xmm9			\n\t"\
			"mulpd	%%xmm15,%%xmm3			\n\t	mulpd	%%xmm15,%%xmm11			\n\t"\
			"/* Digit 2: Require output to be nonnegative, so leave MSW unbalanced: */	\n\t"\
			"mulpd	0x80(%%rdx),%%xmm2		\n\t	mulpd	0x90(%%rdx),%%xmm10		\n\t"\
			"mulpd	0xa0(%%rdx),%%xmm6		\n\t	mulpd	0xb0(%%rdx),%%xmm14		\n\t"\
			"mulpd	0xc0(%%rdx),%%xmm4		\n\t	mulpd	0xd0(%%rdx),%%xmm12		\n\t"\
			"addpd	%%xmm6,%%xmm2			\n\t	addpd	%%xmm14,%%xmm10			\n\t"\
			"addpd	%%xmm4,%%xmm3			\n\t	addpd	%%xmm12,%%xmm11			\n\t"\
			"addpd	%%xmm3,%%xmm2			\n\t	addpd	%%xmm11,%%xmm10			\n\t"\
			"movaps	%%xmm2,%%xmm3			\n\t	movaps	%%xmm10,%%xmm11			\n\t"\
			"subpd	0x240(%%rdx),%%xmm3		\n\t	subpd	0x240(%%rdx),%%xmm11	\n\t"\
			"addpd	%%xmm7,%%xmm3			\n\t	addpd	%%xmm7 ,%%xmm11			\n\t"\
			"subpd	%%xmm7,%%xmm3			\n\t	subpd	%%xmm7 ,%%xmm11			\n\t"\
			"subpd	%%xmm3,%%xmm2			\n\t	subpd	%%xmm11,%%xmm10			\n\t"\
			"mulpd	%%xmm5,%%xmm2			\n\t	mulpd	%%xmm5 ,%%xmm10			\n\t"\
			"/* MULH96(q,lo,lo) --> lo = (q*lo)/2^78 */								\n\t"\
			"/*** movq	%[__fq0],%%rdx		Moved to start of asm-block ***/\n\t"\
			"/* Digit 0: */					\n\t"\
			"movaps	%%xmm0,%%xmm3			\n\t	movaps	%%xmm8 ,%%xmm11			\n\t"\
			"movaps	%%xmm0,%%xmm4			\n\t	movaps	%%xmm8 ,%%xmm12			\n\t"\
			"mulpd	   (%%rdx),%%xmm0		\n\t	mulpd	0x10(%%rdx),%%xmm8		\n\t"\
			"addpd	%%xmm7,%%xmm0			\n\t	addpd	%%xmm7 ,%%xmm8			\n\t"\
			"subpd	%%xmm7,%%xmm0			\n\t	subpd	%%xmm7 ,%%xmm8			\n\t"\
			"mulpd	   %%xmm15,%%xmm0		\n\t	mulpd	   %%xmm15,%%xmm8		\n\t"\
			"/* Digit 1: */					\n\t"\
			"movaps	%%xmm1,%%xmm6			\n\t	movaps	%%xmm9 ,%%xmm14			\n\t"\
			"mulpd	0x20(%%rdx),%%xmm3		\n\t	mulpd	0x30(%%rdx),%%xmm11		\n\t"\
			"addpd	%%xmm0,%%xmm3			\n\t	addpd	%%xmm8 ,%%xmm11			\n\t"\
			"mulpd	   (%%rdx),%%xmm1		\n\t	mulpd	0x10(%%rdx),%%xmm9		\n\t"\
			"addpd	%%xmm3,%%xmm1			\n\t	addpd	%%xmm11,%%xmm9			\n\t"\
			"addpd	%%xmm7,%%xmm1			\n\t	addpd	%%xmm7 ,%%xmm9			\n\t"\
			"subpd	%%xmm7,%%xmm1			\n\t	subpd	%%xmm7 ,%%xmm9			\n\t"\
			"mulpd	   %%xmm15,%%xmm1		\n\t	mulpd	   %%xmm15,%%xmm9		\n\t"\
			"/* Digit 2: Require output to be nonnegative, so leave MSW unbalanced: */	\n\t"\
			"movaps	%%xmm2,%%xmm3			\n\t	movaps	%%xmm10,%%xmm11			\n\t"\
			"movaps	%%xmm6,%%xmm0			\n\t	movaps	%%xmm14,%%xmm8			\n\t"\
			"mulpd	   (%%rdx),%%xmm2		\n\t	mulpd	0x10(%%rdx),%%xmm10		\n\t"\
			"mulpd	0x20(%%rdx),%%xmm6		\n\t	mulpd	0x30(%%rdx),%%xmm14		\n\t"\
			"mulpd	0x40(%%rdx),%%xmm4		\n\t	mulpd	0x50(%%rdx),%%xmm12		\n\t"\
			"addpd	%%xmm6,%%xmm2			\n\t	addpd	%%xmm14,%%xmm10			\n\t"\
			"addpd	%%xmm4,%%xmm1			\n\t	addpd	%%xmm12,%%xmm9			\n\t"\
			"addpd	%%xmm1,%%xmm2			\n\t	addpd	%%xmm9 ,%%xmm10			\n\t"\
			"subpd	0x240(%%rdx),%%xmm2		\n\t	subpd	0x240(%%rdx),%%xmm10		\n\t"\
			"addpd	%%xmm7,%%xmm2			\n\t	addpd	%%xmm7 ,%%xmm10			\n\t"\
			"subpd	%%xmm7,%%xmm2			\n\t	subpd	%%xmm7 ,%%xmm10			\n\t"\
			"mulpd	   %%xmm15,%%xmm2		\n\t	mulpd	   %%xmm15,%%xmm10		\n\t"\
			"/* Precompute all the needed partial products: */						\n\t"\
			"movaps	%%xmm3,%%xmm1			\n\t	movaps	%%xmm11,%%xmm9			\n\t"\
			"mulpd	0x40(%%rdx),%%xmm0		\n\t	mulpd	0x50(%%rdx),%%xmm8		\n\t"\
			"mulpd	0x20(%%rdx),%%xmm1		\n\t	mulpd	0x30(%%rdx),%%xmm9		\n\t"\
			"mulpd	0x40(%%rdx),%%xmm3		\n\t	mulpd	0x50(%%rdx),%%xmm11		\n\t"\
			"/* Digit 3: */					\n\t"\
			"addpd	%%xmm2,%%xmm1			\n\t	addpd	%%xmm10,%%xmm9			\n\t"\
			"addpd	%%xmm1,%%xmm0			\n\t	addpd	%%xmm9 ,%%xmm8			\n\t"\
			"movaps	%%xmm0,%%xmm6			\n\t	movaps	%%xmm8 ,%%xmm14			\n\t"\
			"addpd	%%xmm7,%%xmm6			\n\t	addpd	%%xmm7 ,%%xmm14			\n\t"\
			"subpd	%%xmm7,%%xmm6			\n\t	subpd	%%xmm7 ,%%xmm14			\n\t"\
			"subpd	%%xmm6,%%xmm0			\n\t	subpd	%%xmm14,%%xmm8			\n\t"\
			"mulpd	%%xmm5,%%xmm0			\n\t	mulpd	%%xmm5 ,%%xmm8			\n\t"\
			"mulpd	   %%xmm15,%%xmm6		\n\t	mulpd	   %%xmm15,%%xmm14		\n\t"\
			"/* Digit 4: */					\n\t"\
			"addpd	%%xmm6,%%xmm3			\n\t	addpd	%%xmm14,%%xmm11			\n\t"\
			"movaps	%%xmm3,%%xmm1			\n\t	movaps	%%xmm11,%%xmm9			\n\t"\
			"addpd	%%xmm7,%%xmm3			\n\t	addpd	%%xmm7 ,%%xmm11			\n\t"\
			"subpd	%%xmm7,%%xmm3			\n\t	subpd	%%xmm7 ,%%xmm11			\n\t"\
			"subpd	%%xmm3,%%xmm1			\n\t	subpd	%%xmm11,%%xmm9			\n\t"\
			"mulpd	%%xmm5,%%xmm1			\n\t	mulpd	%%xmm5 ,%%xmm9			\n\t"\
			"/* If h < l, calculate h-l+q; otherwise h-l. Use leading 52 bits to approximate the full 78-bit compare. Result is in [0, q). */\n\t"\
			"movaps %%xmm5,%%xmm13	/* Need a copy of two26f, stored in %%xmm5 */\n\t"\
			"movaps	0x1c0(%%rdx),%%xmm2		\n\t	movaps	0x1d0(%%rdx),%%xmm10	\n\t"\
			"movaps	%%xmm2,%%xmm6			\n\t	movaps	%%xmm10,%%xmm14			\n\t"\
			"movaps	0x60(%%rdx),%%xmm4		\n\t	movaps	0x70(%%rdx),%%xmm12		\n\t"\
			"mulpd	%%xmm5,%%xmm3			\n\t	mulpd	%%xmm5 ,%%xmm11			\n\t"\
			"mulpd	   (%%rdx),%%xmm5		\n\t	mulpd	0x10(%%rdx),%%xmm13		\n\t"\
			"addpd	%%xmm1,%%xmm3			\n\t	addpd	%%xmm9 ,%%xmm11			\n\t"\
			"movaps	0x1a0(%%rdx),%%xmm1		\n\t	movaps	0x1b0(%%rdx),%%xmm9		\n\t"\
			"cmppd	$0x1,%%xmm3,%%xmm6		\n\t	cmppd	$0x1,%%xmm11,%%xmm14	\n\t"\
			"subpd	%%xmm0,%%xmm1			\n\t	subpd	%%xmm8 ,%%xmm9			\n\t"\
			"subpd	%%xmm3,%%xmm2			\n\t	subpd	%%xmm11,%%xmm10			\n\t"\
			"movaps	%%xmm4,%%xmm0			\n\t	movaps	%%xmm12,%%xmm8			\n\t"\
			"andpd	%%xmm6,%%xmm4			\n\t	andpd	%%xmm14,%%xmm12			\n\t"\
			"addpd	%%xmm4,%%xmm2			\n\t	addpd	%%xmm12,%%xmm10			\n\t"\
			"andpd	%%xmm5,%%xmm6			\n\t	andpd	%%xmm13,%%xmm14			\n\t"\
			"addpd	%%xmm6,%%xmm1			\n\t	addpd	%%xmm14,%%xmm9			\n\t"\
			"/* qlo26 in xmm5, qhi52 in xmm0 */		\n\t"\
			"/* if((pshift >> j) & (uint64)1) { */	\n\t"\
			"movslq	%[__pshift],%%rax		\n\t"\
			"movslq	%[__j],%%rcx			\n\t"\
			"shrq	%%cl,%%rax				\n\t"\
			"andq	$0x1,%%rax				\n\t"\
		"je	twopmodq78_3wdq4			\n\t"\
		"/* Double: Use that high part < 2^52 (strictly >= 0): xhi52,xlo26 in xmm2,xmm1; qhi52,qlo26 in xmm0,xmm5 */\n\t"\
		"	movaps		%%xmm0,%%xmm6		\n\t	movaps		%%xmm8 ,%%xmm14	/* cpy of qhi */			\n\t"\
		"	addpd		%%xmm2,%%xmm2		\n\t	addpd		%%xmm10,%%xmm10	/* top 52 bits */			\n\t"\
		"	addpd		%%xmm1,%%xmm1		\n\t	addpd		%%xmm9 ,%%xmm9	/* low 26 bits */			\n\t"\
		"/* If x > q, subtract q: */		\n\t"\
		"	cmppd	$0x2,%%xmm2,%%xmm6		\n\t	cmppd	$0x2,%%xmm10,%%xmm14	/* bitmask = (qhi <= xhi) */\n\t"\
		"	andpd		%%xmm6,%%xmm0		\n\t	andpd		%%xmm14,%%xmm8	/* qhi52 & bitmask */		\n\t"\
		"	andpd		%%xmm6,%%xmm5		\n\t	andpd		%%xmm14,%%xmm13	/* qlo26 & bitmask */		\n\t"\
		"	subpd		%%xmm0,%%xmm2		\n\t	subpd		%%xmm8 ,%%xmm10	/* x mod q, top 52 bits */	\n\t"\
		"	subpd		%%xmm5,%%xmm1		\n\t	subpd		%%xmm13,%%xmm9	/* x mod q, low 26 bits */	\n\t"\
		"twopmodq78_3wdq4:					\n\t"\
			"/* } */						\n\t"\
			"/* Normalize the result: */	\n\t"\
			"movaps	0x210(%%rdx),%%xmm4		\n\t	movaps	0x210(%%rdx),%%xmm12	\n\t"\
			"movaps		%%xmm15,%%xmm3		\n\t	movaps		%%xmm15,%%xmm11	\n\t"\
			"mulpd	%%xmm3,%%xmm1			\n\t	mulpd	%%xmm11,%%xmm9			\n\t"\
			"mulpd	%%xmm3,%%xmm2			\n\t	mulpd	%%xmm11,%%xmm10			\n\t"\
			"movaps	%%xmm1,%%xmm0			\n\t	movaps	%%xmm9 ,%%xmm8			\n\t"\
			"addpd	%%xmm7,%%xmm1			\n\t	addpd	%%xmm7 ,%%xmm9			\n\t"\
			"subpd	%%xmm7,%%xmm1			\n\t	subpd	%%xmm7 ,%%xmm9			\n\t"\
			"subpd	%%xmm1,%%xmm0			\n\t	subpd	%%xmm9 ,%%xmm8			\n\t"\
			"mulpd	%%xmm4,%%xmm0			\n\t	mulpd	%%xmm12,%%xmm8			\n\t"\
			"mulpd	%%xmm3,%%xmm1			\n\t	mulpd	%%xmm11,%%xmm9			\n\t"\
			"addpd	%%xmm1,%%xmm2			\n\t	addpd	%%xmm9 ,%%xmm10			\n\t"\
			"movaps	%%xmm2,%%xmm1			\n\t	movaps	%%xmm10,%%xmm9			\n\t"\
			"addpd	%%xmm7,%%xmm2			\n\t	addpd	%%xmm7 ,%%xmm10			\n\t"\
			"subpd	%%xmm7,%%xmm2			\n\t	subpd	%%xmm7 ,%%xmm10			\n\t"\
			"subpd	%%xmm2,%%xmm1			\n\t	subpd	%%xmm10,%%xmm9			\n\t"\
			"mulpd	%%xmm4,%%xmm1			\n\t	mulpd	%%xmm12,%%xmm9			\n\t"\
			"/* Move high 2 words of result into input registers expected by start of loop body: */\n\t"\
			"movaps	%%xmm2,%%xmm4 /* fx2 */	\n\t	movaps	%%xmm10,%%xmm12	/* hx2 */	\n\t"\
			"movaps	%%xmm1,%%xmm2 /* fx1 */	\n\t	movaps	%%xmm9 ,%%xmm10	/* hx1 */	\n\t"\
			:					/* outputs: none */\
			: [__fq0] "m" (Xfq0)	/* All inputs from memory addresses here */\
			,[__pshift] "m" (Xpshift)	\
			,[__j]		"m" (Xj)		\
			: "cl","rax","rcx","rdx","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
		);\
		}

	#endif	// ifdef(USE_ROUNDPD)

#endif	// 64-bit GCC/SIMD

/****************************************************/
/* Useful 78-bit-int-via-3-floating-doubles Macros: */
/****************************************************/

/* Converts a 78-bit unsigned input __x (stored in a uint96)
to balanced-digit floating-point form. Outputs have the following size ranges:

	fword0,1 in [-2^25, +2^25]
	fword2   in [   -1, +2^26]
*/
#define CVT_UINT78_3WORD_DOUBLE(__x, __fword0, __fword1, __fword2)\
{\
	uint64 __tmp64;\
	int64 __itmp, __cy;\
	\
	DBG_ASSERT(HERE, (__x.d1 >> 14) == 0, "Input > 78-bit limit!");\
	\
	/* Digit 0: */\
	__tmp64 = __x.d0;\
	__itmp = __tmp64 & 0x0000000003ffffff;\
	/* Is the current digit >= (base/2)? */\
	__cy = (int64)((uint64)__itmp>>25);	/* Cast to unsigned to ensure logical right-shift */\
	/* If yes, balance it by subtracting the base: */\
	/* RHS terms must be signed to prevent integer underflow-on-subtract: */\
	__fword0 = (double)(__itmp - (__cy<<26));\
	\
	/* Digit 1: */\
	__tmp64 = (__tmp64 >> 26) + __cy;\
	__itmp = __tmp64 & 0x0000000003ffffff;\
	/* Is the current digit >= (base/2)? */\
	__cy = (int64)((uint64)__itmp>>25);	/* Cast to unsigned to ensure logical right-shift */\
	/* If yes, balance it by subtracting the base: */\
	/* RHS terms must be signed to prevent integer underflow-on-subtract: */\
	__fword1 = (double)(__itmp - (__cy<<26));\
	\
	/* Digit 2: */\
	__tmp64 = (__tmp64 >> 26) + __cy;\
	__tmp64 += (__x.d1 << 12);	/* 12 = (64 - 2*26) */\
	/* No balanced-digit normalization of MSW: */\
	__fword2 = (double)__tmp64;\
	\
	DBG_ASSERT(HERE, __fword2 <= TWO26FLOAT, "MSW > TWO26FLOAT");\
}

/* Converts a 78-bit unsigned input __x (stored in balanced-digit
floating-point form) to a uint96. Assumes the FP input is properly normalized.
*/
#define CVT78_3WORD_DOUBLE_UINT96(__fword0, __fword1, __fword2, __x)\
{\
	int64 __itmp, __cy;\
	\
	/* Cast current digit to int64 form, subtracting any borrow from previous digit: */\
	__itmp = (int64)__fword0;\
	if(__itmp < 0)	/* If current digit < 0, add the base and set carry = -1	*/\
	{\
		__itmp += TWO26FLOAT;\
		DBG_ASSERT(HERE, __itmp >= 0, "Normalized digit still < 0!");\
		__cy = -1;\
	}\
	else\
	{\
		__cy = 0;\
	}\
	__x.d0 = (uint64)__itmp;\
\
	/* Digit 1: */\
	__itmp = (int64)__fword1 +  __cy;\
	if(__itmp < 0)\
	{\
		__itmp += TWO26FLOAT;\
		DBG_ASSERT(HERE, __itmp >= 0, "Normalized digit still < 0!");\
		__cy = -1;\
	}\
	else\
	{\
		__cy = 0;\
	}\
	__x.d0 += ((uint64)__itmp << 26);\
\
	/* Digit 2: */\
	__itmp = (int64)__fword2 +  __cy;\
	if(__itmp < 0)\
	{\
		__itmp += TWO26FLOAT;\
		DBG_ASSERT(HERE, __itmp >= 0, "Normalized digit still < 0!");\
		__cy = -1;\
	}\
	else\
	{\
		__cy = 0;\
	}\
	__x.d0 += ((uint64)__itmp << 52);\
	__x.d1  = ((uint64)__itmp >> 12) & 0x0000000000003fff;	/* Only case where we really need the (uint64) cast */\
	\
	DBG_ASSERT(HERE, (__x.d1 >> 14) == 0, "Output > 78-bit limit!");\
	DBG_ASSERT(HERE,  __cy          == 0, "Nonzero exit carry!");\
}

/* Takes a 78-bit unsigned input __x stored in balanced-digit floating-point form
and renormalizes with respect to the balanced-digit base.
*/
#define NORMALIZE78_3WORD_DOUBLE(__x0, __x1, __x2)\
{\
	double __fcy;\
	\
	/* Digit 0: */\
	__fcy = DNINT(__x0*TWO26FLINV);\
	__x0 -= __fcy*TWO26FLOAT;\
	\
	/* Digit 1: */\
	__x1 += __fcy;\
	__fcy = DNINT(__x1*TWO26FLINV);\
	__x1 -= __fcy*TWO26FLOAT;\
	\
	/* Digit 2: */\
	__x2 += __fcy;\
	\
	DBG_ASSERT(HERE, __x2 <= TWO26FLOAT, "MSW > TWO26FLOAT");\
	DBG_ASSERT(HERE, __x2 >= 0         , "MSW < 0!");\
}

/* Takes a 156-bit unsigned input __x stored in balanced-digit floating-point form
and renormalizes with respect to an 78-bit = (26,26,26)-bit balanced-digit base.
Because we expect that we may wind up using the upper and lower halves of the result
separately, we require the MSW of each to be nonnegative, i.e. we don't balance __x2.
*/
#define NORMALIZE156_6WORD_DOUBLE(__x0, __x1, __x2, __x3, __x4, __x5)\
{\
	double __fcy;\
	\
	/* Digit 0: */\
	__fcy = DNINT(__x0*TWO26FLINV);\
	__x0 -= __fcy*TWO26FLOAT;\
	\
	/* Digit 1: */\
	__x1 += __fcy;\
	__fcy = DNINT(__x1*TWO26FLINV);\
	__x1 -= __fcy*TWO26FLOAT;\
	\
	/* Digit 2: */\
	__x2 += __fcy;\
	__fcy = DNINT(__x2*TWO26FLINV);\
	__x2 -= __fcy*TWO26FLOAT;\
	if(__x2 < 0)\
	{\
		__x2 += TWO26FLOAT;\
		__fcy--;\
	}\
	\
	/* Digit 3: */\
	__x3 += __fcy;\
	__fcy = DNINT(__x3*TWO26FLINV);\
	__x3 -= __fcy*TWO26FLOAT;\
	\
	/* Digit 4: */\
	__x4 += __fcy;\
	__fcy = DNINT(__x4*TWO26FLINV);\
	__x4 -= __fcy*TWO26FLOAT;\
	\
	/* Digit 5: */\
	__x5 += __fcy;\
	\
	DBG_ASSERT(HERE, __x2 >= 0         , "_x2 < 0!");\
	DBG_ASSERT(HERE, __x2 <= TWO26FLOAT, "_x2 > TWO26FLOAT");\
	DBG_ASSERT(HERE, __x5 >= 0         , "MSW < 0!");\
	DBG_ASSERT(HERE, __x5 <= TWO26FLOAT, "MSW > TWO26FLOAT");\
}

/**********************************************************************************/
/* Balanced-digit FP arithmetic with base 2^26 allows MULTIPLY inputs up to 2^78. */
/**********************************************************************************/

/* Product of 78-bit x and y, stored in semi-balanced-digit 3-word-double-form with respect to base B := 2^b.
The 3-word-double inout operands are "semi-balanced", which in general means the MSW is nonnegative (and is
xero iff the entire operand is 0), but in this branchless optimized form means components in the following size ranges:

		fword0,1 in [-B/2, +B/2]
		fword2   in [  -1, +B/2]

	Let x = x0 + x1*B + x2*B^2,
	and y = y0 + y1*B + y2*B^2.

In terms of the 6 output coefficients (of which we need only the upper 3), here are the resulting size bounds
prior to the carry-step renormalization. We neglect the incoming carry (i.e. the (...).hi + ... term) in w1-5
to a first approximation, because we shall soon see that it cannot possibly affect the result:

	Output coefficient                                      Range
	----------------------------------------------------    ------------------------------
	w0 = (x0*y0).lo                                         [-B^2/4, +B^2/4]
	w1 = (x0*y0).hi + (x0*y1 + x1*y0).lo                    [-B^2/4, +B^2/4]*2
	w2 = (x0*y1 + x1*y0).hi + (x0*y2 + x1*y1 + x2*y0).lo    [-B^2/2, +B^2/2]*2 + [-B^2/4, +B^2/4]
	w3 = (x0*y2 + x1*y1 + x2*y0).hi + (x2*y1 + x1*y2).lo    [-B^2/4, +B^2/4]*2
	w4 = (x2*y1 + x1*y2).hi + (x2*y2).lo                    [-B^2/4, +B^2/4]
	w5 = (x2*y2).hi ,

where x.(lo,hi) denote the lower and upper 26 bits of the 52-bit product term x.

It is clear the w2 term is the one which limits B - This is in [-B^2, +B^2]*(3/4) < 2^53, giving b_max = 26.
*/

/*...Square of a 78-bit input __x .
Lower and upper halves of __x^2 are returned in __lo and __hi, respectively.
Current version needs 16 FMUL, 12 FADD, several cast-between int-and-double, several ALU.

Because we expect that we may wind up using the upper and lower halves of the result
separately, we require the MSW of each to be nonnegative, i.e. we don't balance __x2.

ASSUMES:
	- None of the input and output addresses coincide;
*/
#define SQR_LOHI78_3WORD_DOUBLE(\
  __fx0, __fx1, __fx2, __fprod0, __fprod1, __fprod2, __fprod3, __fprod4, __fprod5\
)\
{\
	double __f2x0 = __fx0 + __fx0, __f2x1 = __fx1 + __fx1;\
	double __fcy;\
	uint32 __itmp;\
	\
	DBG_ASSERT(HERE, __fx0 < TWO26FLOAT, "x0 > TWO26FLOAT");\
	DBG_ASSERT(HERE, __fx1 < TWO26FLOAT, "x1 > TWO26FLOAT");\
	DBG_ASSERT(HERE, __fx2 < TWO26FLOAT, "x2 > TWO26FLOAT");\
	\
	/* Digit 0: */\
	__fprod0  =  __fx0*__fx0;\
	__fcy    = DNINT(__fprod0*TWO26FLINV);\
	__fprod0 -= __fcy*TWO26FLOAT;\
	\
	/* Digit 1: */\
	__fprod1  = __f2x0*__fx1 + __fcy;\
	__fcy    = DNINT(__fprod1*TWO26FLINV);\
	__fprod1 -= __fcy*TWO26FLOAT;\
	\
	/* Digit 2: */\
	__fprod2  = __f2x0*__fx2 + __fx1*__fx1 + __fcy;\
	__fcy    = DNINT(__fprod2*TWO26FLINV);\
	__fprod2 -= __fcy*TWO26FLOAT;\
	/* Branchless sequence to unbalance the __fprod2 term: */\
	__itmp    = ((__fprod2*TWO26FLOAT + __fprod1) < 0);\
	__fcy   -= (double)__itmp;\
	__fprod2 += (double)(__itmp << 26);\
	\
	/* Digit 3: */\
	__fprod3  = __f2x1*__fx2 + __fcy;\
	__fcy    = DNINT(__fprod3*TWO26FLINV);\
	__fprod3 -= __fcy*TWO26FLOAT;\
	\
	/* Digit 4: */\
	__fprod4  =  __fx2*__fx2 + __fcy;\
	__fcy    = DNINT(__fprod4*TWO26FLINV);\
	__fprod4 -= __fcy*TWO26FLOAT;\
	\
	/* Digit 5: */\
	__fprod5  = __fcy;\
	\
	DBG_ASSERT(HERE, __fprod2 >= 0         , "_x2 < 0!");\
	DBG_ASSERT(HERE, __fprod2 <= TWO26FLOAT, "_x2 > TWO26FLOAT");\
	DBG_ASSERT(HERE, __fprod5 >= 0         , "MSW < 0!");\
	DBG_ASSERT(HERE, __fprod5 <= TWO26FLOAT, "MSW > TWO26FLOAT");\
}

#define SQR_LOHI78_3WORD_DOUBLE_q2(\
  __fx0,__fx1,__fx2, __fprod0,__fprod1,__fprod2,__fprod3,__fprod4,__fprod5\
, __gx0,__gx1,__gx2, __gprod0,__gprod1,__gprod2,__gprod3,__gprod4,__gprod5\
)\
{\
	double __f2x0 = __fx0 + __fx0, __f2x1 = __fx1 + __fx1;\
	double __g2x0 = __gx0 + __gx0, __g2x1 = __gx1 + __gx1;\
	double __fcy, __gcy;\
	uint32 __itmp, __jtmp;\
	\
	/* Digit 0: */\
	__fprod0  =  __fx0*__fx0;							__gprod0  =  __gx0*__gx0;\
	__fcy     = DNINT(__fprod0*TWO26FLINV);				__gcy     = DNINT(__gprod0*TWO26FLINV);\
	__fprod0 -= __fcy*TWO26FLOAT;						__gprod0 -= __gcy*TWO26FLOAT;\
	/* Digit 1: */\
	__fprod1  = __f2x0*__fx1 + __fcy;					__gprod1  = __g2x0*__gx1 + __gcy;\
	__fcy     = DNINT(__fprod1*TWO26FLINV);				__gcy     = DNINT(__gprod1*TWO26FLINV);\
	__fprod1 -= __fcy*TWO26FLOAT;						__gprod1 -= __gcy*TWO26FLOAT;\
	/* Digit 2: */\
	__fprod2  = __f2x0*__fx2 + __fx1*__fx1 + __fcy;		__gprod2  = __g2x0*__gx2 + __gx1*__gx1 + __gcy;\
	__fcy     = DNINT(__fprod2*TWO26FLINV);				__gcy     = DNINT(__gprod2*TWO26FLINV);\
	__fprod2 -= __fcy*TWO26FLOAT;						__gprod2 -= __gcy*TWO26FLOAT;\
	__itmp    = ((__fprod2*TWO26FLOAT + __fprod1) < 0);	__jtmp    = ((__gprod2*TWO26FLOAT + __gprod1) < 0);\
	__fcy    -= (double)__itmp;							__gcy    -= (double)__jtmp;\
	__fprod2 += (double)(__itmp << 26);					__gprod2 += (double)(__jtmp << 26);\
	/* Digit 3: */\
	__fprod3  = __f2x1*__fx2 + __fcy;					__gprod3  = __g2x1*__gx2 + __gcy;\
	__fcy     = DNINT(__fprod3*TWO26FLINV);				__gcy     = DNINT(__gprod3*TWO26FLINV);\
	__fprod3 -= __fcy*TWO26FLOAT;						__gprod3 -= __gcy*TWO26FLOAT;\
	/* Digit 4: */\
	__fprod4  =  __fx2*__fx2 + __fcy;					__gprod4  =  __gx2*__gx2 + __gcy;\
	__fcy     = DNINT(__fprod4*TWO26FLINV);				__gcy     = DNINT(__gprod4*TWO26FLINV);\
	__fprod4 -= __fcy*TWO26FLOAT;						__gprod4 -= __gcy*TWO26FLOAT;\
	/* Digit 5: */\
	__fprod5  = __fcy;									__gprod5  = __gcy;\
}

#define SQR_LOHI78_3WORD_DOUBLE_q4(\
  __fx0,__fx1,__fx2, __fprod0,__fprod1,__fprod2,__fprod3,__fprod4,__fcy\
, __gx0,__gx1,__gx2, __gprod0,__gprod1,__gprod2,__gprod3,__gprod4,__gcy\
, __hx0,__hx1,__hx2, __hprod0,__hprod1,__hprod2,__hprod3,__hprod4,__hcy\
, __ix0,__ix1,__ix2, __iprod0,__iprod1,__iprod2,__iprod3,__iprod4,__icy\
)\
{\
	double __f2x0 = __fx0 + __fx0, __f2x1 = __fx1 + __fx1;\
	double __g2x0 = __gx0 + __gx0, __g2x1 = __gx1 + __gx1;\
	double __h2x0 = __hx0 + __hx0, __h2x1 = __hx1 + __hx1;\
	double __i2x0 = __ix0 + __ix0, __i2x1 = __ix1 + __ix1;\
	uint32 __itmp, __jtmp, __ktmp, __ltmp;\
	\
	/* Digit 0: */\
	__fprod0  =  __fx0*__fx0;							__gprod0  =  __gx0*__gx0;							__hprod0  =  __hx0*__hx0;							__iprod0  =  __ix0*__ix0;\
	__fcy     = DNINT(__fprod0*TWO26FLINV);				__gcy     = DNINT(__gprod0*TWO26FLINV);				__hcy     = DNINT(__hprod0*TWO26FLINV);				__icy     = DNINT(__iprod0*TWO26FLINV);\
	__fprod0 -= __fcy*TWO26FLOAT;						__gprod0 -= __gcy*TWO26FLOAT;						__hprod0 -= __hcy*TWO26FLOAT;						__iprod0 -= __icy*TWO26FLOAT;\
	/* Digit 1: */\
	__fprod1  = __f2x0*__fx1 + __fcy;					__gprod1  = __g2x0*__gx1 + __gcy;					__hprod1  = __h2x0*__hx1 + __hcy;					__iprod1  = __i2x0*__ix1 + __icy;\
	__fcy     = DNINT(__fprod1*TWO26FLINV);				__gcy     = DNINT(__gprod1*TWO26FLINV);				__hcy     = DNINT(__hprod1*TWO26FLINV);				__icy     = DNINT(__iprod1*TWO26FLINV);\
	__fprod1 -= __fcy*TWO26FLOAT;						__gprod1 -= __gcy*TWO26FLOAT;						__hprod1 -= __hcy*TWO26FLOAT;						__iprod1 -= __icy*TWO26FLOAT;\
	/* Digit 2: */\
	__fprod2  = __f2x0*__fx2 + __fx1*__fx1 + __fcy;		__gprod2  = __g2x0*__gx2 + __gx1*__gx1 + __gcy;		__hprod2  = __h2x0*__hx2 + __hx1*__hx1 + __hcy;		__iprod2  = __i2x0*__ix2 + __ix1*__ix1 + __icy;\
	__fcy     = DNINT(__fprod2*TWO26FLINV);				__gcy     = DNINT(__gprod2*TWO26FLINV);				__hcy     = DNINT(__hprod2*TWO26FLINV);				__icy     = DNINT(__iprod2*TWO26FLINV);\
	__fprod2 -= __fcy*TWO26FLOAT;						__gprod2 -= __gcy*TWO26FLOAT;						__hprod2 -= __hcy*TWO26FLOAT;						__iprod2 -= __icy*TWO26FLOAT;\
__itmp = (__fprod2 < 0);							__jtmp = (__gprod2 < 0);							__ktmp = (__hprod2 < 0);							__ltmp = (__iprod2 < 0);\
/*	__itmp    = ((__fprod2*TWO26FLOAT + __fprod1) < 0);	__jtmp    = ((__gprod2*TWO26FLOAT + __gprod1) < 0);	__ktmp    = ((__hprod2*TWO26FLOAT + __hprod1) < 0);	__ltmp    = ((__iprod2*TWO26FLOAT + __iprod1) < 0);*/\
	__fcy    -= (double)__itmp;							__gcy    -= (double)__jtmp;							__hcy    -= (double)__ktmp;							__icy    -= (double)__ltmp;\
	__fprod2 += (double)(__itmp << 26);					__gprod2 += (double)(__jtmp << 26);					__hprod2 += (double)(__ktmp << 26);					__iprod2 += (double)(__ltmp << 26);\
	/* Digit 3: */\
	__fprod3  = __f2x1*__fx2 + __fcy;					__gprod3  = __g2x1*__gx2 + __gcy;					__hprod3  = __h2x1*__hx2 + __hcy;					__iprod3  = __i2x1*__ix2 + __icy;\
	__fcy     = DNINT(__fprod3*TWO26FLINV);				__gcy     = DNINT(__gprod3*TWO26FLINV);				__hcy     = DNINT(__hprod3*TWO26FLINV);				__icy     = DNINT(__iprod3*TWO26FLINV);\
	__fprod3 -= __fcy*TWO26FLOAT;						__gprod3 -= __gcy*TWO26FLOAT;						__hprod3 -= __hcy*TWO26FLOAT;						__iprod3 -= __icy*TWO26FLOAT;\
	/* Digit 4,5: */\
	__fprod4  =  __fx2*__fx2 + __fcy;					__gprod4  =  __gx2*__gx2 + __gcy;					__hprod4  =  __hx2*__hx2 + __hcy;					__iprod4  =  __ix2*__ix2 + __icy;\
	__fcy     = DNINT(__fprod4*TWO26FLINV);				__gcy     = DNINT(__gprod4*TWO26FLINV);				__hcy     = DNINT(__hprod4*TWO26FLINV);				__icy     = DNINT(__iprod4*TWO26FLINV);\
	__fprod4 -= __fcy*TWO26FLOAT;						__gprod4 -= __gcy*TWO26FLOAT;						__hprod4 -= __hcy*TWO26FLOAT;						__iprod4 -= __icy*TWO26FLOAT;\
}

#define SQR_LOHI78_3WORD_DOUBLE_LEAVE_HIGH52_UNNORMALIZED_q4(\
  __fx0,__fx1,__fx2, __fprod0,__fprod1,__fprod2,__fprod3,__fprod4\
, __gx0,__gx1,__gx2, __gprod0,__gprod1,__gprod2,__gprod3,__gprod4\
, __hx0,__hx1,__hx2, __hprod0,__hprod1,__hprod2,__hprod3,__hprod4\
, __ix0,__ix1,__ix2, __iprod0,__iprod1,__iprod2,__iprod3,__iprod4\
)\
{\
	double __f2x0 = __fx0 + __fx0, __f2x1 = __fx1 + __fx1;\
	double __g2x0 = __gx0 + __gx0, __g2x1 = __gx1 + __gx1;\
	double __h2x0 = __hx0 + __hx0, __h2x1 = __hx1 + __hx1;\
	double __i2x0 = __ix0 + __ix0, __i2x1 = __ix1 + __ix1;\
	double __fcy, __gcy, __hcy, __icy;\
	uint32 __itmp, __jtmp, __ktmp, __ltmp;\
	\
	/* Digit 0: */\
	__fprod0  =  __fx0*__fx0;						__gprod0  =  __gx0*__gx0;						__hprod0  =  __hx0*__hx0;						__iprod0  =  __ix0*__ix0;\
	__fcy     = DNINT(__fprod0*TWO26FLINV);			__gcy     = DNINT(__gprod0*TWO26FLINV);			__hcy     = DNINT(__hprod0*TWO26FLINV);			__icy     = DNINT(__iprod0*TWO26FLINV);\
	__fprod0 -= __fcy*TWO26FLOAT;					__gprod0 -= __gcy*TWO26FLOAT;					__hprod0 -= __hcy*TWO26FLOAT;					__iprod0 -= __icy*TWO26FLOAT;\
	/* Digit 1: */\
	__fprod1  = __f2x0*__fx1 + __fcy;				__gprod1  = __g2x0*__gx1 + __gcy;				__hprod1  = __h2x0*__hx1 + __hcy;				__iprod1  = __i2x0*__ix1 + __icy;\
	__fcy     = DNINT(__fprod1*TWO26FLINV);			__gcy     = DNINT(__gprod1*TWO26FLINV);			__hcy     = DNINT(__hprod1*TWO26FLINV);			__icy     = DNINT(__iprod1*TWO26FLINV);\
	__fprod1 -= __fcy*TWO26FLOAT;					__gprod1 -= __gcy*TWO26FLOAT;					__hprod1 -= __hcy*TWO26FLOAT;					__iprod1 -= __icy*TWO26FLOAT;\
	/* Digit 2: */\
	__fprod2  = __f2x0*__fx2 + __fx1*__fx1 + __fcy;	__gprod2  = __g2x0*__gx2 + __gx1*__gx1 + __gcy;	__hprod2  = __h2x0*__hx2 + __hx1*__hx1 + __hcy;	__iprod2  = __i2x0*__ix2 + __ix1*__ix1 + __icy;\
	__fcy     = DNINT(__fprod2*TWO26FLINV);			__gcy     = DNINT(__gprod2*TWO26FLINV);			__hcy     = DNINT(__hprod2*TWO26FLINV);			__icy     = DNINT(__iprod2*TWO26FLINV);\
	__fprod2 -= __fcy*TWO26FLOAT;					__gprod2 -= __gcy*TWO26FLOAT;					__hprod2 -= __hcy*TWO26FLOAT;					__iprod2 -= __icy*TWO26FLOAT;\
	__itmp = ((__fprod2*TWO26FLOAT + __fprod1) < 0);__jtmp = ((__gprod2*TWO26FLOAT + __gprod1) < 0);__ktmp = ((__hprod2*TWO26FLOAT + __hprod1) < 0);__ltmp = ((__iprod2*TWO26FLOAT + __iprod1) < 0);\
	__fcy    -= (double)__itmp;						__gcy    -= (double)__jtmp;						__hcy    -= (double)__ktmp;						__icy    -= (double)__ltmp;\
	__fprod2 += (double)(__itmp << 26);				__gprod2 += (double)(__jtmp << 26);				__hprod2 += (double)(__ktmp << 26);				__iprod2 += (double)(__ltmp << 26);\
	/* Digit 3: */\
	__fprod3  = __f2x1*__fx2 + __fcy;				__gprod3  = __g2x1*__gx2 + __gcy;				__hprod3  = __h2x1*__hx2 + __hcy;				__iprod3  = __i2x1*__ix2 + __icy;\
	__fcy     = DNINT(__fprod3*TWO26FLINV);			__gcy     = DNINT(__gprod3*TWO26FLINV);			__hcy     = DNINT(__hprod3*TWO26FLINV);			__icy     = DNINT(__iprod3*TWO26FLINV);\
	__fprod3 -= __fcy*TWO26FLOAT;					__gprod3 -= __gcy*TWO26FLOAT;					__hprod3 -= __hcy*TWO26FLOAT;					__iprod3 -= __icy*TWO26FLOAT;\
	/* Digits 4,5 remain in a 52-bit double: */\
	__fprod4  =  __fx2*__fx2 + __fcy;				__gprod4  =  __gx2*__gx2 + __gcy;				__hprod4  =  __hx2*__hx2 + __hcy;				__iprod4  =  __ix2*__ix2 + __icy;\
}

/* Lower half of __x * __y .
Current version needs 12 FMUL, 14 FADD.

Because we may desire to overwrite one of the two sets of inputs with the outputs,
we code so that any or all of __X, __Y and __LO may have the same addresses.
*/
#define MULL78_3WORD_DOUBLE(__x0,__x1,__x2, __y0,__y1,__y2, __lo0,__lo1,__lo2)\
{\
	double __fcy, __prod0, __prod1, __prod2;\
	uint32 __itmp;\
	\
	DBG_ASSERT(HERE, __x0 < TWO26FLOAT, "x0 > TWO26FLOAT");\
	DBG_ASSERT(HERE, __x1 < TWO26FLOAT, "x1 > TWO26FLOAT");\
	DBG_ASSERT(HERE, __x2 < TWO26FLOAT, "x2 > TWO26FLOAT");\
	\
	DBG_ASSERT(HERE, __y0 < TWO26FLOAT, "y0 > TWO26FLOAT");\
	DBG_ASSERT(HERE, __y1 < TWO26FLOAT, "y1 > TWO26FLOAT");\
	DBG_ASSERT(HERE, __y2 < TWO26FLOAT, "y2 > TWO26FLOAT");\
	\
	/* Precompute all the needed partial products: */\
	__prod0  = __x0*__y0;\
	__prod1  = __x0*__y1 + __x1*__y0;\
	__prod2  = __x0*__y2 + __x1*__y1 + __x2*__y0;\
	\
	/* Digit 0: */\
	__fcy    = DNINT(__prod0*TWO26FLINV);\
	__lo0    = __prod0 - __fcy*TWO26FLOAT;\
	\
	/* Digit 1: */\
	__prod1 += __fcy;\
	__fcy    = DNINT(__prod1*TWO26FLINV);\
	__lo1    = __prod1 - __fcy*TWO26FLOAT;\
	\
	/* Digit 2: */\
	__prod2 += __fcy;\
	__fcy    = DNINT(__prod2*TWO26FLINV);\
	__lo2    = __prod2 - __fcy*TWO26FLOAT;\
	/* Branchless sequence to unbalance the __lo2 term: */\
	__itmp   = (__lo2 < 0);\
	__lo2   += (double)(__itmp << 26);\
	/* Require output to be nonnegative, so leave MSW unbalanced: */\
	DBG_ASSERT(HERE, __lo2 >= 0, "MSW < 0!");\
}

#define MULL78_3WORD_DOUBLE_q2(\
  __fx0,__fx1,__fx2, __fy0,__fy1,__fy2, __flo0,__flo1,__flo2\
, __gx0,__gx1,__gx2, __gy0,__gy1,__gy2, __glo0,__glo1,__glo2\
)\
{\
	double __fcy, __fprod0, __fprod1, __fprod2;\
	double __gcy, __gprod0, __gprod1, __gprod2;\
	uint32 __itmp, __jtmp;\
	\
	/* Precompute all the needed partial products: */\
	__fprod0  = __fx0*__fy0;								__gprod0  = __gx0*__gy0;\
	__fprod1  = __fx0*__fy1 + __fx1*__fy0;					__gprod1  = __gx0*__gy1 + __gx1*__gy0;\
	__fprod2  = __fx0*__fy2 + __fx1*__fy1 + __fx2*__fy0;	__gprod2  = __gx0*__gy2 + __gx1*__gy1 + __gx2*__gy0;\
	/* Digit 0: */\
	__fcy     = DNINT(__fprod0*TWO26FLINV);					__gcy     = DNINT(__gprod0*TWO26FLINV);\
	__flo0    = __fprod0 - __fcy*TWO26FLOAT;				__glo0    = __gprod0 - __gcy*TWO26FLOAT;\
	/* Digit 1: */\
	__fprod1 += __fcy;										__gprod1 += __gcy;\
	__fcy     = DNINT(__fprod1*TWO26FLINV);					__gcy     = DNINT(__gprod1*TWO26FLINV);\
	__flo1    = __fprod1 - __fcy*TWO26FLOAT;				__glo1    = __gprod1 - __gcy*TWO26FLOAT;\
	/* Digit 2: */\
	__fprod2 += __fcy;										__gprod2 += __gcy;\
	__fcy     = DNINT(__fprod2*TWO26FLINV);					__gcy     = DNINT(__gprod2*TWO26FLINV);\
	__flo2    = __fprod2 - __fcy*TWO26FLOAT;				__glo2    = __gprod2 - __gcy*TWO26FLOAT;\
	/* Branchless sequence to unbalance the __fprod2 term: */\
	__itmp    = (__flo2 < 0);								__jtmp    = (__glo2 < 0);\
	__flo2 += (double)(__itmp << 26);						__glo2 += (double)(__jtmp << 26);\
}

#define MULL78_3WORD_DOUBLE_q4(\
  __fx0,__fx1,__fx2, __fy0,__fy1,__fy2, __flo0,__flo1,__flo2\
, __gx0,__gx1,__gx2, __gy0,__gy1,__gy2, __glo0,__glo1,__glo2\
, __hx0,__hx1,__hx2, __hy0,__hy1,__hy2, __hlo0,__hlo1,__hlo2\
, __ix0,__ix1,__ix2, __iy0,__iy1,__iy2, __ilo0,__ilo1,__ilo2\
)\
{\
	double __fprod0, __fprod1, __fprod2;\
	double __gprod0, __gprod1, __gprod2;\
	double __hprod0, __hprod1, __hprod2;\
	double __iprod0, __iprod1, __iprod2;\
	double __fcy, __gcy, __hcy, __icy;\
	uint32 __itmp, __jtmp, __ktmp, __ltmp;\
	\
	/* Precompute all the needed partial products: */\
	__fprod0  = __fx0*__fy0;								__gprod0  = __gx0*__gy0;								__hprod0  = __hx0*__hy0;								__iprod0  = __ix0*__iy0;							\
	__fprod1  = __fx0*__fy1 + __fx1*__fy0;					__gprod1  = __gx0*__gy1 + __gx1*__gy0;					__hprod1  = __hx0*__hy1 + __hx1*__hy0;					__iprod1  = __ix0*__iy1 + __ix1*__iy0;				\
	__fprod2  = __fx0*__fy2 + __fx1*__fy1 + __fx2*__fy0;	__gprod2  = __gx0*__gy2 + __gx1*__gy1 + __gx2*__gy0;	__hprod2  = __hx0*__hy2 + __hx1*__hy1 + __hx2*__hy0;	__iprod2  = __ix0*__iy2 + __ix1*__iy1 + __ix2*__iy0;\
	/* Digit 0: */\
	__fcy     = DNINT(__fprod0*TWO26FLINV);					__gcy     = DNINT(__gprod0*TWO26FLINV);					__hcy     = DNINT(__hprod0*TWO26FLINV);					__icy     = DNINT(__iprod0*TWO26FLINV);\
	__flo0    = __fprod0 - __fcy*TWO26FLOAT;				__glo0    = __gprod0 - __gcy*TWO26FLOAT;				__hlo0    = __hprod0 - __hcy*TWO26FLOAT;				__ilo0    = __iprod0 - __icy*TWO26FLOAT;\
	/* Digit 1: */\
	__fprod1 += __fcy;										__gprod1 += __gcy;										__hprod1 += __hcy;										__iprod1 += __icy;\
	__fcy     = DNINT(__fprod1*TWO26FLINV);					__gcy     = DNINT(__gprod1*TWO26FLINV);					__hcy     = DNINT(__hprod1*TWO26FLINV);					__icy     = DNINT(__iprod1*TWO26FLINV);\
	__flo1    = __fprod1 - __fcy*TWO26FLOAT;				__glo1    = __gprod1 - __gcy*TWO26FLOAT;				__hlo1    = __hprod1 - __hcy*TWO26FLOAT;				__ilo1    = __iprod1 - __icy*TWO26FLOAT;\
	/* Digit 2: */\
	__fprod2 += __fcy;										__gprod2 += __gcy;										__hprod2 += __hcy;										__iprod2 += __icy;\
	__fcy     = DNINT(__fprod2*TWO26FLINV);					__gcy     = DNINT(__gprod2*TWO26FLINV);					__hcy     = DNINT(__hprod2*TWO26FLINV);					__icy     = DNINT(__iprod2*TWO26FLINV);\
	__flo2    = __fprod2 - __fcy*TWO26FLOAT;				__glo2    = __gprod2 - __gcy*TWO26FLOAT;				__hlo2    = __hprod2 - __hcy*TWO26FLOAT;				__ilo2    = __iprod2 - __icy*TWO26FLOAT;\
	/* Branchless sequence to unbalance the __fprod2 term: */\
	__itmp    = (__flo2 < 0);								__jtmp    = (__glo2 < 0);								__ktmp    = (__hlo2 < 0);								__ltmp    = (__ilo2 < 0);\
	__flo2 += (double)(__itmp << 26);						__glo2 += (double)(__jtmp << 26);						__hlo2 += (double)(__ktmp << 26);						__ilo2 += (double)(__ltmp << 26);\
}

/* Upper half of __x * __y .
Current version needs 16 FMUL, 22 FADD.

NOTE: 52x26-bit version below (in the _q4 versdion of this macro) needs just 11 FMUL, 11 FADD, 1 compare, 2 casts.

Because we may desire to overwrite one of the two sets of inputs with the outputs,
we code so that any or all of __X, __Y and __LO may have the same addresses.
*/
#define MULH78_3WORD_DOUBLE(__x0,__x1,__x2, __y0,__y1,__y2, __hi0,__hi1,__hi2)\
{\
	double __fcy, __tmp, __prod3, __prod4;\
	uint32 __itmp;\
	\
	DBG_ASSERT(HERE, __x0 < TWO26FLOAT, "x0 > TWO26FLOAT");\
	DBG_ASSERT(HERE, __x1 < TWO26FLOAT, "x1 > TWO26FLOAT");\
	DBG_ASSERT(HERE, __x2 < TWO26FLOAT, "x2 > TWO26FLOAT");\
	\
	DBG_ASSERT(HERE, __y0 < TWO26FLOAT, "y0 > TWO26FLOAT");\
	DBG_ASSERT(HERE, __y1 < TWO26FLOAT, "y1 > TWO26FLOAT");\
	DBG_ASSERT(HERE, __y2 < TWO26FLOAT, "y2 > TWO26FLOAT");\
	\
	/* Digit 0: */\
	__tmp  =  __x0*__y0;\
	__fcy  = DNINT(__tmp*TWO26FLINV);\
	\
	/* Digit 1: */\
	__tmp  = __x0*__y1 + __x1*__y0 + __fcy;\
	__fcy  = DNINT(__tmp*TWO26FLINV);\
	\
	/* Digit 2: */\
	__tmp  = __x0*__y2 + __x1*__y1 + __x2*__y0 + __fcy;\
	__fcy  = DNINT(__tmp*TWO26FLINV);\
	__tmp -= __fcy*TWO26FLOAT;\
	/* Branchless sequence to unbalance the __prod2 term: */\
	__itmp = (__tmp < 0);\
	__fcy -= (double)__itmp;\
	/* Require low half to be nonnegative, so leave this term unbalanced: */\
	/*if(__tmp < 0)	*/\
	/*{				*/\
	/*	__fcy--;	*/\
	/*}				*/\
	\
	/* At this point the possibility of same-address in-and-outputs comes into play: */\
	/* Precompute all the needed partial products: */\
	__prod3 = __x1*__y2 + __x2*__y1 + __fcy;\
	__prod4 = __x2*__y2;\
	\
	/* Digit 3: */\
	__fcy    = DNINT(__prod3*TWO26FLINV);\
	__hi0    = __prod3 - __fcy*TWO26FLOAT;\
	\
	/* Digit 4: */\
	__prod4 += __fcy;\
	__fcy    = DNINT(__prod4*TWO26FLINV);\
	__hi1    = __prod4 - __fcy*TWO26FLOAT;\
	\
	/* Digit 5: */\
	__hi2    = __fcy;\
	\
	DBG_ASSERT(HERE, __hi2 >= 0, "MSW < 0!");\
}

#define MULH78_3WORD_DOUBLE_q2(\
  __fx0,__fx1,__fx2, __fy0,__fy1,__fy2, __fhi0,__fhi1,__fhi2\
, __gx0,__gx1,__gx2, __gy0,__gy1,__gy2, __ghi0,__ghi1,__ghi2\
  )\
{\
	double __fcy, __ftmp, __fprod3, __fprod4;\
	double __gcy, __gtmp, __gprod3, __gprod4;\
	uint32 __itmp,__jtmp;\
	\
	/* Digit 0: */\
	__ftmp =  __fx0*__fy0;										__gtmp =  __gx0*__gy0;\
	__fcy  = DNINT(__ftmp*TWO26FLINV);							__gcy  = DNINT(__gtmp*TWO26FLINV);\
	\
	/* Digit 1: */\
	__ftmp = __fx0*__fy1 + __fx1*__fy0 + __fcy;					__gtmp = __gx0*__gy1 + __gx1*__gy0 + __gcy;\
	__fcy  = DNINT(__ftmp*TWO26FLINV);							__gcy  = DNINT(__gtmp*TWO26FLINV);\
	\
	/* Digit 2: */\
	__ftmp = __fx0*__fy2 + __fx1*__fy1 + __fx2*__fy0 + __fcy;	__gtmp = __gx0*__gy2 + __gx1*__gy1 + __gx2*__gy0 + __gcy;\
	__fcy  = DNINT(__ftmp*TWO26FLINV);							__gcy  = DNINT(__gtmp*TWO26FLINV);\
	__ftmp-= __fcy*TWO26FLOAT;									__gtmp-= __gcy*TWO26FLOAT;\
	__itmp = (__ftmp < 0);										__jtmp = (__gtmp < 0);\
	__fcy -= (double)__itmp;									__gcy -= (double)__jtmp;\
	\
	/* At this point the possibility of same-address in-and-outputs comes into play: */\
	/* Precompute all the needed partial products: */\
	__fprod3 = __fx1*__fy2 + __fx2*__fy1 + __fcy;				__gprod3 = __gx1*__gy2 + __gx2*__gy1 + __gcy;\
	__fprod4 = __fx2*__fy2;										__gprod4 = __gx2*__gy2;\
	\
	/* Digit 3: */\
	__fcy  = DNINT(__fprod3*TWO26FLINV);						__gcy  = DNINT(__gprod3*TWO26FLINV);\
	__fhi0 = __fprod3 - __fcy*TWO26FLOAT;						__ghi0 = __gprod3 - __gcy*TWO26FLOAT;\
	\
	/* Digit 4: */\
	__fprod4 += __fcy;											__gprod4 += __gcy;\
	__fcy  = DNINT(__fprod4*TWO26FLINV);						__gcy  = DNINT(__gprod4*TWO26FLINV);\
	__fhi1 = __fprod4 - __fcy*TWO26FLOAT;						__ghi1 = __gprod4 - __gcy*TWO26FLOAT;\
	\
	/* Digit 5: */\
	__fhi2 = __fcy;												__ghi2 = __gcy;\
}

/* 3 Versions of MULH78: First, the more-expensive exact version: Cost = 19 FMUL, 14 FADD, 5 DNINT: */
	#define MULH78_3WORD_DOUBLE_q4(\
	 __fx0,__fx1,__fx2, __fy0,__fy1,__fy2, __fhi0,__fhi1,__fhi2\
	, __gx0,__gx1,__gx2, __gy0,__gy1,__gy2, __ghi0,__ghi1,__ghi2\
	, __hx0,__hx1,__hx2, __hy0,__hy1,__hy2, __hhi0,__hhi1,__hhi2\
	, __ix0,__ix1,__ix2, __iy0,__iy1,__iy2, __ihi0,__ihi1,__ihi2\
	 )\
	{\
		double __fprod0, __fprod1, __fprod2, __fprod3, __fprod4,__fcy;\
		double __gprod0, __gprod1, __gprod2, __gprod3, __gprod4,__gcy;\
		double __hprod0, __hprod1, __hprod2, __hprod3, __hprod4,__hcy;\
		double __iprod0, __iprod1, __iprod2, __iprod3, __iprod4,__icy;\
		uint32 __itmp, __jtmp, __ktmp, __ltmp;\
		/* Digit 0: */\
		__fprod0 =  __fx0*__fy0;								__gprod0 =  __gx0*__gy0;								__hprod0 =  __hx0*__hy0;								__iprod0 =  __ix0*__iy0;\
		__fcy     = DNINT(__fprod0*TWO26FLINV);					__gcy     = DNINT(__gprod0*TWO26FLINV);					__hcy     = DNINT(__hprod0*TWO26FLINV);					__icy     = DNINT(__iprod0*TWO26FLINV);\
		/* Digit 1: */\
		__fprod1 = __fx0*__fy1 + __fx1*__fy0 + __fcy;			__gprod1 = __gx0*__gy1 + __gx1*__gy0 + __gcy;			__hprod1 = __hx0*__hy1 + __hx1*__hy0 + __hcy;			__iprod1 = __ix0*__iy1 + __ix1*__iy0 + __icy;\
		__fcy     = DNINT(__fprod1*TWO26FLINV);					__gcy     = DNINT(__gprod1*TWO26FLINV);					__hcy     = DNINT(__hprod1*TWO26FLINV);					__icy     = DNINT(__iprod1*TWO26FLINV);\
		__fprod1 -= __fcy*TWO26FLOAT;							__gprod1 -= __gcy*TWO26FLOAT;							__hprod1 -= __hcy*TWO26FLOAT;							__iprod1 -= __icy*TWO26FLOAT;\
		/* Digit 2: */\
		__fprod2 = __fx0*__fy2+__fx1*__fy1+__fx2*__fy0+__fcy;	__gprod2 = __gx0*__gy2+__gx1*__gy1+__gx2*__gy0+__gcy;	__hprod2 = __hx0*__hy2+__hx1*__hy1+__hx2*__hy0+__hcy;	__iprod2 = __ix0*__iy2+__ix1*__iy1+__ix2*__iy0+__icy;\
		__fcy     = DNINT(__fprod2*TWO26FLINV);					__gcy     = DNINT(__gprod2*TWO26FLINV);					__hcy     = DNINT(__hprod2*TWO26FLINV);					__icy     = DNINT(__iprod2*TWO26FLINV);\
		__fprod2-= __fcy*TWO26FLOAT;							__gprod2-= __gcy*TWO26FLOAT;							__hprod2-= __hcy*TWO26FLOAT;							__iprod2-= __icy*TWO26FLOAT;\
		__itmp    = ((__fprod2*TWO26FLOAT + __fprod1) < 0);		__jtmp    = ((__gprod2*TWO26FLOAT + __gprod1) < 0);		__ktmp    = ((__hprod2*TWO26FLOAT + __hprod1) < 0);		__ltmp    = ((__iprod2*TWO26FLOAT + __iprod1) < 0);\
		__fcy    -= (double)__itmp;								__gcy    -= (double)__jtmp;								__hcy    -= (double)__ktmp;								__icy    -= (double)__ltmp;\
		/* At this point the possibility of same-address in-and-outputs comes into play: */\
		/* Precompute all the needed partial products: */\
		__fprod3 = __fx1*__fy2 + __fx2*__fy1 + __fcy;			__gprod3 = __gx1*__gy2 + __gx2*__gy1 + __gcy;			__hprod3 = 	__hx1*__hy2 + __hx2*__hy1 + __hcy;			__iprod3 = __ix1*__iy2 + __ix2*__iy1 + __icy;\
		__fprod4 = __fx2*__fy2;									__gprod4 = __gx2*__gy2;									__hprod4 = 	__hx2*__hy2;								__iprod4 = __ix2*__iy2;\
		/* Digit 3: */\
		__fcy     = DNINT(__fprod3*TWO26FLINV);					__gcy     = DNINT(__gprod3*TWO26FLINV);					__hcy     = DNINT(__hprod3*TWO26FLINV);					__icy     = DNINT(__iprod3*TWO26FLINV);\
		__fhi0 = __fprod3 - __fcy*TWO26FLOAT;					__ghi0 = __gprod3 - __gcy*TWO26FLOAT;					__hhi0 = 	__hprod3 - __hcy*TWO26FLOAT;				__ihi0 = __iprod3 - __icy*TWO26FLOAT;\
		/* Digit 4: */\
		__fprod4 += __fcy;										__gprod4 += __gcy;										__hprod4 += 	__hcy;									__iprod4 += __icy;\
		__fcy     = DNINT(__fprod4*TWO26FLINV);					__gcy     = DNINT(__gprod4*TWO26FLINV);					__hcy     = DNINT(__hprod4*TWO26FLINV);					__icy     = DNINT(__iprod4*TWO26FLINV);\
		__fhi1 = __fprod4 - __fcy*TWO26FLOAT;					__ghi1 = __gprod4 - __gcy*TWO26FLOAT;					__hhi1 = 	__hprod4 - __hcy*TWO26FLOAT;				__ihi1 = __iprod4 - __icy*TWO26FLOAT;\
		__fhi2    = __fcy;                      				__ghi2    = __gcy;   	                				__hhi2    = __hcy;    	               				__ihi2    = __icy;\
	}

	/* Cheaper version which uses approximate carry into upper half. Cost = 15 FMUL, 11 FADD, 3 DNINT: */
	#define MULH78_3WORD_DOUBLE_q4_v2(\
	 __fx0,__fx1,__fx2, __fy0,__fy1,__fy2, __fhi0,__fhi1,__fhi2\
	, __gx0,__gx1,__gx2, __gy0,__gy1,__gy2, __ghi0,__ghi1,__ghi2\
	, __hx0,__hx1,__hx2, __hy0,__hy1,__hy2, __hhi0,__hhi1,__hhi2\
	, __ix0,__ix1,__ix2, __iy0,__iy1,__iy2, __ihi0,__ihi1,__ihi2\
	 )\
	{\
		double __fcy, __ftmp,__fprod3, __fprod4;\
		double __gcy, __gtmp,__gprod3, __gprod4;\
		double __hcy, __htmp,__hprod3, __hprod4;\
		double __icy, __itmp,__iprod3, __iprod4;\
		\
		/* Digit 1: */\
		__ftmp = __fx0*__fy1 + __fx1*__fy0        ;				__gtmp = __gx0*__gy1 + __gx1*__gy0        ;				__htmp = __hx0*__hy1 + __hx1*__hy0        ;				__itmp = __ix0*__iy1 + __ix1*__iy0        ;\
		__fcy     =      (__ftmp*TWO26FLINV);					__gcy     =      (__gtmp*TWO26FLINV);					__hcy     =      (__htmp*TWO26FLINV);					__icy     =      (__itmp*TWO26FLINV);\
		/* Digit 2: */\
		__ftmp = __fx0*__fy2+__fx1*__fy1+__fx2*__fy0+__fcy;		__gtmp = __gx0*__gy2+__gx1*__gy1+__gx2*__gy0+__gcy;		__htmp = __hx0*__hy2+__hx1*__hy1+__hx2*__hy0+__hcy;		__itmp = __ix0*__iy2+__ix1*__iy1+__ix2*__iy0+__icy;\
		__fcy     = DNINT(__ftmp*TWO26FLINV);					__gcy     = DNINT(__gtmp*TWO26FLINV);					__hcy     = DNINT(__htmp*TWO26FLINV);					__icy     = DNINT(__itmp*TWO26FLINV);\
		__ftmp-= __fcy*TWO26FLOAT;								__gtmp-= __gcy*TWO26FLOAT;								__htmp-= __hcy*TWO26FLOAT;								__itmp-= __icy*TWO26FLOAT;\
		__fcy -= (double)(__ftmp < 0);							__gcy -= (double)(__gtmp < 0);							__hcy -= (double)(__htmp < 0);							__icy -= (double)(__itmp < 0);\
		/* At this point the possibility of same-address in-and-outputs comes into play: */\
		/* Precompute all the needed partial products: */\
		__fhi0 = __fx1*__fy2 + __fx2*__fy1 + __fcy;				__ghi0 = __gx1*__gy2 + __gx2*__gy1 + __gcy;				__hhi0 = __hx1*__hy2 + __hx2*__hy1 + __hcy;				__ihi0 = __ix1*__iy2 + __ix2*__iy1 + __icy;\
		__fhi1 = __fx2*__fy2;									__ghi1 = __gx2*__gy2;									__hhi1 = __hx2*__hy2;									__ihi1 = __ix2*__iy2;\
		/* Digit 3: */\
		__fcy     = DNINT(__fhi0*TWO26FLINV);					__gcy     = DNINT(__ghi0*TWO26FLINV);					__hcy     = DNINT(__hhi0*TWO26FLINV);					__icy     = DNINT(__ihi0*TWO26FLINV);\
		__fhi0 -= __fcy*TWO26FLOAT;								__ghi0 -= __gcy*TWO26FLOAT;								__hhi0 -= __hcy*TWO26FLOAT;								__ihi0 -= __icy*TWO26FLOAT;\
		/* Digit 4,5: */\
		__fhi1 += __fcy;										__ghi1 += __gcy;										__hhi1 += __hcy;										__ihi1 += __icy;\
		__fhi2     = DNINT(__fhi1*TWO26FLINV);					__ghi2     = DNINT(__ghi1*TWO26FLINV);					__hhi2     = DNINT(__hhi1*TWO26FLINV);					__ihi2     = DNINT(__ihi1*TWO26FLINV);\
		__fhi1 -= __fhi2*TWO26FLOAT;							__ghi1 -= __ghi2*TWO26FLOAT;							__hhi1 -= __hhi2*TWO26FLOAT;							__ihi1 -= __ihi2*TWO26FLOAT;\
	}

	/* Even-Cheaper version which leaves high 52 output bits as-is. Cost = 13 FMUL, 10 FADD, 2 DNINT: */
	#define MULH78_3WORD_DOUBLE_LEAVE_HIGH52_UNNORMALIZED_q4(\
	 __fx0,__fx1,__fx2, __fy0,__fy1,__fy2, __fhi0,__fhi1\
	, __gx0,__gx1,__gx2, __gy0,__gy1,__gy2, __ghi0,__ghi1\
	, __hx0,__hx1,__hx2, __hy0,__hy1,__hy2, __hhi0,__hhi1\
	, __ix0,__ix1,__ix2, __iy0,__iy1,__iy2, __ihi0,__ihi1\
	 )\
	{\
		double __ftmp, __gtmp, __htmp, __itmp;\
		double __fcy, __gcy, __hcy, __icy;\
		\
		/* Digit 1: */\
		__ftmp = __fx0*__fy1;				__gtmp = __gx0*__gy1;				__htmp = __hx0*__hy1;				__itmp = __ix0*__iy1;				\
		__ftmp+= __fx1*__fy0;				__gtmp+= __gx1*__gy0;				__htmp+= __hx1*__hy0;				__itmp+= __ix1*__iy0;				\
		__ftmp*= TWO26FLINV;				__gtmp*= TWO26FLINV;				__htmp*= TWO26FLINV;				__itmp*= TWO26FLINV;				\
		/* Digit 2: */\
		__ftmp+= __fx0*__fy2;				__gtmp+= __gx0*__gy2;				__htmp+= __hx0*__hy2;				__itmp+= __ix0*__iy2;				\
		__ftmp+= __fx1*__fy1;				__gtmp+= __gx1*__gy1;				__htmp+= __hx1*__hy1;				__itmp+= __ix1*__iy1;				\
		__ftmp+= __fx2*__fy0;				__gtmp+= __gx2*__gy0;				__htmp+= __hx2*__hy0;				__itmp+= __ix2*__iy0;				\
		__fcy  = DNINT(__ftmp*TWO26FLINV);	__gcy  = DNINT(__gtmp*TWO26FLINV);	__hcy  = DNINT(__htmp*TWO26FLINV);	__icy  = DNINT(__itmp*TWO26FLINV);	\
		__ftmp-= __fcy*TWO26FLOAT;			__gtmp-= __gcy*TWO26FLOAT;			__htmp-= __hcy*TWO26FLOAT;			__itmp-= __icy*TWO26FLOAT;			\
		__fcy -= (double)(__ftmp < 0);		__gcy -= (double)(__gtmp < 0);		__hcy -= (double)(__htmp < 0);		__icy -= (double)(__itmp < 0);		\
		/* At this point the possibility of same-address in-and-outputs comes into play: */\
		/* Precompute all the needed partial products: */\
		__ftmp = __fx1*__fy2 + __fcy;		__gtmp = __gx1*__gy2 + __gcy;		__htmp = __hx1*__hy2 + __hcy;		__itmp = __ix1*__iy2 + __icy;		\
		__ftmp+= __fx2*__fy1;				__gtmp+= __gx2*__gy1;				__htmp+= __hx2*__hy1;				__itmp+= __ix2*__iy1;				\
		__fhi1 = __fx2*__fy2;				__ghi1 = __gx2*__gy2;				__hhi1 = __hx2*__hy2;				__ihi1 = __ix2*__iy2;				\
		/* Digit 3: */\
		__fcy  = DNINT(__ftmp*TWO26FLINV);	__gcy  = DNINT(__gtmp*TWO26FLINV);	__hcy  = DNINT(__htmp*TWO26FLINV);	__icy  = DNINT(__itmp*TWO26FLINV);	\
		__fhi0 = __ftmp - __fcy*TWO26FLOAT;	__ghi0 = __gtmp - __gcy*TWO26FLOAT;	__hhi0 = __htmp - __hcy*TWO26FLOAT;	__ihi0 = __itmp - __icy*TWO26FLOAT;	\
		/* Digits 4,5 remain in a 52-bit double: */\
		__fhi1 += __fcy;					__ghi1 += __gcy;					__hhi1 += __hcy;					__ihi1 += __icy;					\
	}

#if 0
/*
void foobar(double __fx52, double __fy2, double __fx2, double __fy52)
{
	uint64  a,b,c,d,lo26;
	uint128 i128,j128,out128;
a = (uint64)__fx52; b = (uint64)__fy2; c = (uint64)__fx2; d = (uint64)__fy52;\
MUL_LOHI64(a,b,i128.d0,i128.d1);\
MUL_LOHI64(c,d,j128.d0,j128.d1);\
ADD128(i128,j128,j128);\
out128.d0 = (uint64)__fcy; out128.d1 = (uint64)0ull;\
ADD128(j128,out128,out128);\
lo26 = out128.d0 & 0x0000000003FFFFFFull;\
RSHIFT128(out128, 26, out128);\
fprintf(stderr,"exact<52:77>, <78:129> = %20llu, %20llu\n",out128.d0,out128.d1);\
}
*/
	#define MULH78_3WORD_DOUBLE_q4(\
	 __fx0,__fx1,__fx2, __fy0,__fy1,__fy2, __fhi0,__fhi1,__fhi2\
	, __gx0,__gx1,__gx2, __gy0,__gy1,__gy2, __ghi0,__ghi1,__ghi2\
	, __hx0,__hx1,__hx2, __hy0,__hy1,__hy2, __hhi0,__hhi1,__hhi2\
	, __ix0,__ix1,__ix2, __iy0,__iy1,__iy2, __ihi0,__ihi1,__ihi2\
	 )\
	{\
		uint64  a,b,c,d,lo26;\
		uint128 i128,j128,out128;\
		double __ftmp, __fprod4;\
		double __gtmp, __gprod4;\
		double __htmp, __hprod4;\
		double __itmp, __iprod4;\
		double __fcy, __gcy, __hcy, __icy;\
		double __fx52 = __fx0 + __fx1*TWO26FLOAT;\
		double __gx52 = __gx0 + __gx1*TWO26FLOAT;\
		double __hx52 = __hx0 + __hx1*TWO26FLOAT;\
		double __ix52 = __ix0 + __ix1*TWO26FLOAT;\
		double __fy52 = __fy0 + __fy1*TWO26FLOAT;\
		double __gy52 = __gy0 + __gy1*TWO26FLOAT;\
		double __hy52 = __hy0 + __hy1*TWO26FLOAT;\
		double __iy52 = __iy0 + __iy1*TWO26FLOAT;\
		/* Bottom 104 bits, with upper 52 bits, i.e. <52:103> correct: */\
		__ftmp = __fx52*__fy52;\
		__gtmp = __gx52*__gy52;\
		__htmp = __hx52*__hy52;\
		__itmp = __ix52*__iy52;\
		__fcy  = __ftmp*TWO52FLINV;\
		__gcy  = __gtmp*TWO52FLINV;\
		__hcy  = __htmp*TWO52FLINV;\
		__icy  = __itmp*TWO52FLINV;\
		/* 78-bit subproducts starting at bit 52, with upper 52 bits, i.e. <78:129> correct: */\
a = (uint64)__fx52; b = (uint64)__fy2; c = (uint64)__fx2; d = (uint64)__fy52;\
MUL_LOHI64(a,b,i128.d0,i128.d1);\
MUL_LOHI64(c,d,j128.d0,j128.d1);\
ADD128(i128,j128,j128);\
out128.d0 = (uint64)__fcy; out128.d1 = (uint64)0ull;\
ADD128(j128,out128,out128);\
lo26 = out128.d0 & 0x0000000003FFFFFFull;\
RSHIFT128(out128, 26, out128);\
fprintf(stderr,"exact<52:77>, <78:129> = %20llu, %20llu\n",out128.d0,out128.d1);\
		__ftmp = __fx52*__fy2+__fx2*__fy52+__fcy;\
		__gtmp = __gx52*__gy2+__gx2*__gy52+__gcy;\
		__htmp = __hx52*__hy2+__hx2*__hy52+__hcy;\
		__itmp = __ix52*__iy2+__ix2*__iy52+__icy;\
		__fhi0 = DNINT(__ftmp*TWO26FLINV);	/* bits <78:129> */\
/*fprintf(stderr,"__fx1*__fy2 + __fx2*__fy1, __fcy = %25.2f, %25.2f\n", (__fx52*__fy2+__fx2*__fy52)*TWO26FLINV, __fcy*TWO26FLINV);*/\
		__ghi0 = DNINT(__gtmp*TWO26FLINV);\
		__hhi0 = DNINT(__htmp*TWO26FLINV);\
		__ihi0 = DNINT(__itmp*TWO26FLINV);\
		__ftmp-= __fhi0*TWO26FLOAT;\
		__gtmp-= __ghi0*TWO26FLOAT;\
		__htmp-= __hhi0*TWO26FLOAT;\
		__itmp-= __ihi0*TWO26FLOAT;\
		__fhi0 -= (double)(__ftmp < 0);\
/*fprintf(stderr,"__fhi0*TWO26FLOAT, ?<0 = %20.5f, %20.5f\n", __ftmp, __fhi0);*/\
fprintf(stderr,"bits <52:77>, <78:129> = %25.2f, %25.2f\n", __ftmp, __fhi0);\
		__ghi0 -= (double)(__gtmp < 0);\
		__hhi0 -= (double)(__htmp < 0);\
		__ihi0 -= (double)(__itmp < 0);\
		/* At this point the possibility of same-address in-and-outputs comes into play: */\
		/* Precompute all the needed partial products: */\
		__fprod4 = __fx2*__fy2;\
		__gprod4 = __gx2*__gy2;\
		__hprod4 = __hx2*__hy2;\
		__iprod4 = __ix2*__iy2;\
		/* Digit 3: */\
		__fcy  = DNINT(__fhi0*TWO26FLINV);\
/*fprintf(stderr,"__fhi0*TWO26FLINV, rnd = %20.5f, %20.5f\n", __fhi0*TWO26FLINV, __fcy);*/\
		__gcy  = DNINT(__ghi0*TWO26FLINV);\
		__hcy  = DNINT(__hhi0*TWO26FLINV);\
		__icy  = DNINT(__ihi0*TWO26FLINV);\
		__fhi0 -= __fcy*TWO26FLOAT;\
		__ghi0 -= __gcy*TWO26FLOAT;\
		__hhi0 -= __hcy*TWO26FLOAT;\
		__ihi0 -= __icy*TWO26FLOAT;\
		/* Digit 4: */\
		__fprod4 += __fcy;\
		__gprod4 += __gcy;\
		__hprod4 += __hcy;\
		__iprod4 += __icy;\
		/* Digit 5: */\
		__fhi2 = DNINT(__fprod4*TWO26FLINV);\
		__ghi2 = DNINT(__gprod4*TWO26FLINV);\
		__hhi2 = DNINT(__hprod4*TWO26FLINV);\
		__ihi2 = DNINT(__iprod4*TWO26FLINV);\
		__fhi1 = __fprod4 - __fhi2*TWO26FLOAT;\
		__ghi1 = __gprod4 - __ghi2*TWO26FLOAT;\
		__hhi1 = __hprod4 - __hhi2*TWO26FLOAT;\
		__ihi1 = __iprod4 - __ihi2*TWO26FLOAT;\
	}

/*
Testing 63-bit factors...
bits <52:77>, <78:129> =               -8388608.00,        766800360225472.00
outs =      -15353152.00000,        4956677.00000,            771.00000
bits <52:77>, <78:129> =               12582912.00,        230767036792465.00
outs =      -12517743.00000,        6610254.00000,            232.00000
bits <52:77>, <78:129> =               -8388608.00,        475070094820834.00
outs =        4113890.00000,      -15984694.00000,            478.00000
bits <52:77>, <78:129> =              -25165824.00,        963535107048456.00
outs =       -3619832.00000,       -5519355.00000,            969.00000
bits <52:77>, <78:129> =               16777216.00,       1565116592251788.00
outs =      -24943732.00000,      -12839927.00000,           1574.00000
bits <52:77>, <78:129> =              -25690112.00,         54122087152284.00
outs =       -3704164.00000,       28471982.00000,             54.00000
bits <52:77>, <78:129> =                      0.00,        266929927896509.00
outs =       -7848515.00000,       31434267.00000,            268.00000
bits <52:77>, <78:129> =                      0.00,        870051964696220.00
outs =      -28658020.00000,       -9711616.00000,            875.00000
bits <52:77>, <78:129> =              -16777216.00,       1621177384534797.00
outs =      -31489267.00000,       16559967.00000,           1630.00000
bits <52:77>, <78:129> =               25165824.00,        550079871260311.00
outs =      -11331945.00000,        9939198.00000,            553.00000
bits <52:77>, <78:129> =              -33554432.00,       1411719962872699.00
outs =      -18297989.00000,      -24486390.00000,           1420.00000
bits <52:77>, <78:129> =                8388608.00,        756830204706573.00
outs =       -8274163.00000,        3013905.00000,            761.00000
bits <52:77>, <78:129> =               -4194304.00,        320046132618562.00
outs =         779586.00000,       -9755764.00000,            322.00000
bits <52:77>, <78:129> =              -30408704.00,         65659214510417.00
outs =      -30918319.00000,         732337.00000,             66.00000
bits <52:77>, <78:129> =               20971520.00,        482301445212355.00
outs =      -28243773.00000,       -2449501.00000,            485.00000
bits <52:77>, <78:129> =               33554432.00,        705144653961898.00
outs =       10312362.00000,        7043103.00000,            709.00000
bits <52:77>, <78:129> =              -25165824.00,        829375789605950.00
outs =       22225982.00000,       -3938278.00000,            834.00000
*/
#endif

#define	CMPLT78_3WORD_DOUBLE(__x0,__x1,__x2, __y0,__y1,__y2)	(__x2 < __y2 || (__x2 == __y2 && __x1 < __y1) || (__x2 == __y2 && __x1 == __y1 && __x0 < __y0))

#define	CMPGT78_3WORD_DOUBLE(__x0,__x1,__x2, __y0,__y1,__y2)	(__x2 > __y2 || (__x2 == __y2 && __x1 > __y1) || (__x2 == __y2 && __x1 == __y1 && __x0 > __y0))

#define ADD78_3WORD_DOUBLE(__x0,__x1,__x2, __y0,__y1,__y2, __z0,__z1,__z2)\
{\
	__z0 = __x0 + __y0;\
	__z1 = __x1 + __y1;\
	__z2 = __x2 + __y2;\
}

#define SUB78_3WORD_DOUBLE(__x0,__x1,__x2, __y0,__y1,__y2, __z0,__z1,__z2)\
{\
	__z0 = __x0 - __y0;\
	__z1 = __x1 - __y1;\
	__z2 = __x2 - __y2;\
}

#ifdef __cplusplus
}
#endif

#endif	/* twopmodq80_h_included */

