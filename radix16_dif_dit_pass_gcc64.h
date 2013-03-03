/*******************************************************************************
*                                                                              *
*   (C) 1997-2009 by Ernst W. Mayer.                                           *
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
#ifndef radix16_dif_dit_pass_gcc_h_included
#define radix16_dif_dit_pass_gcc_h_included

  #define USE_64BIT_ASM_STYLE	1	// My x86 timings indicate fancier versions below using xmm0-15 for the radix-32 DFT is faster.

  #if !USE_64BIT_ASM_STYLE	// False: Use simple 64-bit-ified version of the analogous 32-bit ASM macros, i.e. using just xmm0-7.

	#define SSE2_RADIX16_DIF_TWIDDLE(Xadd0,Xp1,Xp2,Xp3,Xp4,Xp8,Xp12,Xr1,Xisrt2)\
	{\
	__asm__ volatile (\
		"movq	%[__add0],%%rax		\n\t"\
		"movslq	%[__p4],%%rbx		\n\t"\
		"movslq	%[__p8],%%rcx		\n\t"\
		"movslq	%[__p12],%%rdx		\n\t"\
		"shlq	$3,%%rbx			\n\t"\
		"shlq	$3,%%rcx			\n\t"\
		"shlq	$3,%%rdx			\n\t"\
		"movq	%[__isrt2],%%rsi 		\n\t"\
		"addq	%%rax,%%rbx			\n\t"\
		"addq	%%rax,%%rcx			\n\t"\
		"addq	%%rax,%%rdx			\n\t"\
		"/* SSE2_RADIX4_DIF_4TWIDDLE_B(r1 ,c0) */\n\t"\
		"/* Do	the p0,p8 combo: */	\n\t"\
		"addq	$0x30,%%rsi 	/* c0 */\n\t"\
		"movaps	    (%%rax),%%xmm0	\n\t"\
		"movaps	    (%%rcx),%%xmm4	\n\t"\
		"movaps	0x10(%%rax),%%xmm1	\n\t"\
		"movaps	0x10(%%rcx),%%xmm5	\n\t"\
		"movaps	    (%%rsi),%%xmm6	\n\t"\
		"movaps	0x10(%%rsi),%%xmm7	\n\t"\
		"movaps	%%xmm0,%%xmm2		\n\t"\
		"movaps	%%xmm1,%%xmm3		\n\t"\
		"mulpd	%%xmm6,%%xmm0		\n\t"\
		"mulpd	%%xmm6,%%xmm1		\n\t"\
		"mulpd	%%xmm7,%%xmm2		\n\t"\
		"mulpd	%%xmm7,%%xmm3		\n\t"\
		"movaps	%%xmm4,%%xmm6		\n\t"\
		"addpd	%%xmm2,%%xmm1		\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t"\
		"mulpd	0x20(%%rsi),%%xmm4	\n\t"\
		"subpd	%%xmm3,%%xmm0		\n\t"\
		"mulpd	0x20(%%rsi),%%xmm5	\n\t"\
		"mulpd	0x30(%%rsi),%%xmm6	\n\t"\
		"movaps	%%xmm0,%%xmm2		\n\t"\
		"mulpd	0x30(%%rsi),%%xmm7	\n\t"\
		"addpd	%%xmm6,%%xmm5		\n\t"\
		"movaps	%%xmm1,%%xmm3		\n\t"\
		"subpd	%%xmm7,%%xmm4		\n\t"\
		"addpd	%%xmm4,%%xmm0		\n\t"\
		"addpd	%%xmm5,%%xmm1		\n\t"\
		"subpd	%%xmm4,%%xmm2		\n\t"\
		"subpd	%%xmm5,%%xmm3		\n\t"\
		"/* Do	the p4,12 combo: */	\n\t"\
		"movaps	    (%%rdx),%%xmm4	\n\t"\
		"movaps	0x10(%%rdx),%%xmm5	\n\t"\
		"movaps	    (%%rdx),%%xmm6	\n\t"\
		"movaps	0x10(%%rdx),%%xmm7	\n\t"\
		"mulpd	0x60(%%rsi),%%xmm4	/* c12 */\n\t"\
		"mulpd	0x60(%%rsi),%%xmm5	\n\t"\
		"mulpd	0x70(%%rsi),%%xmm6	\n\t"\
		"mulpd	0x70(%%rsi),%%xmm7	\n\t"\
		"movq	%[__r1],%%rdi 		\n\t"\
		"addpd	%%xmm6,%%xmm5		\n\t"\
		"subpd	%%xmm7,%%xmm4		\n\t"\
		"movaps	%%xmm5,0x010(%%rdi)	\n\t"\
		"movaps	%%xmm4,     (%%rdi)	\n\t"\
		"addq	$0x40,%%rsi 		/* c4 */\n\t"\
		"movaps	    (%%rbx),%%xmm4	\n\t"\
		"movaps	0x10(%%rbx),%%xmm5	\n\t"\
		"movaps	    (%%rbx),%%xmm6	\n\t"\
		"movaps	0x10(%%rbx),%%xmm7	\n\t"\
		"mulpd	    (%%rsi),%%xmm4	\n\t"\
		"mulpd	    (%%rsi),%%xmm5	\n\t"\
		"mulpd	0x10(%%rsi),%%xmm6	\n\t"\
		"mulpd	0x10(%%rsi),%%xmm7	\n\t"\
		"addpd	%%xmm6,%%xmm5		\n\t"\
		"subpd	%%xmm7,%%xmm4		\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t"\
		"movaps	%%xmm4,%%xmm6		\n\t"\
		"subpd	     (%%rdi),%%xmm4	\n\t"\
		"subpd	0x010(%%rdi),%%xmm5	\n\t"\
		"addpd	     (%%rdi),%%xmm6	\n\t"\
		"addpd	0x010(%%rdi),%%xmm7	\n\t"\
		"/* Finish radix-4 butterfly and store results into temporary-array slots: */\n\t"\
		"subpd	%%xmm6,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm2		\n\t"\
		"subpd	%%xmm7,%%xmm1		\n\t"\
		"subpd	%%xmm4,%%xmm3		\n\t"\
		"movaps	%%xmm0,0x040(%%rdi)	\n\t"\
		"movaps	%%xmm2,0x020(%%rdi)	\n\t"\
		"movaps	%%xmm1,0x050(%%rdi)	\n\t"\
		"movaps	%%xmm3,0x070(%%rdi)	\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm0,%%xmm6		\n\t"\
		"addpd	%%xmm2,%%xmm5		\n\t"\
		"addpd	%%xmm1,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm4		\n\t"\
		"movaps	%%xmm6,     (%%rdi)	\n\t"\
		"movaps	%%xmm5,0x060(%%rdi)	\n\t"\
		"movaps	%%xmm7,0x010(%%rdi)	\n\t"\
		"movaps	%%xmm4,0x030(%%rdi)	\n\t"\
		"movslq	%[__p2],%%rdi		\n\t"\
		"shlq	$3,%%rdi			\n\t"\
		"addq	%%rdi,%%rax			\n\t"\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	%%rdi,%%rcx			\n\t"\
		"addq	%%rdi,%%rdx			\n\t"\
		"/* SSE2_RADIX4_DIF_4TWIDDLE_B(r9 ,c2) */\n\t"\
		"/* Do	the p0,p8 combo: */	\n\t"\
		"addq	$0x40,%%rsi 		/* c2 */\n\t"\
		"movaps	    (%%rax),%%xmm0	\n\t"\
		"movaps	    (%%rcx),%%xmm4	\n\t"\
		"movaps	0x10(%%rax),%%xmm1	\n\t"\
		"movaps	0x10(%%rcx),%%xmm5	\n\t"\
		"movaps	    (%%rsi),%%xmm6	\n\t"\
		"movaps	0x10(%%rsi),%%xmm7	\n\t"\
		"movaps	%%xmm0,%%xmm2		\n\t"\
		"movaps	%%xmm1,%%xmm3		\n\t"\
		"mulpd	%%xmm6,%%xmm0		\n\t"\
		"mulpd	%%xmm6,%%xmm1		\n\t"\
		"mulpd	%%xmm7,%%xmm2		\n\t"\
		"mulpd	%%xmm7,%%xmm3		\n\t"\
		"movaps	%%xmm4,%%xmm6		\n\t"\
		"addpd	%%xmm2,%%xmm1		\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t"\
		"mulpd	0x20(%%rsi),%%xmm4	\n\t"\
		"subpd	%%xmm3,%%xmm0		\n\t"\
		"mulpd	0x20(%%rsi),%%xmm5	\n\t"\
		"mulpd	0x30(%%rsi),%%xmm6	\n\t"\
		"movaps	%%xmm0,%%xmm2		\n\t"\
		"mulpd	0x30(%%rsi),%%xmm7	\n\t"\
		"addpd	%%xmm6,%%xmm5		\n\t"\
		"movaps	%%xmm1,%%xmm3		\n\t"\
		"subpd	%%xmm7,%%xmm4		\n\t"\
		"addpd	%%xmm4,%%xmm0		\n\t"\
		"addpd	%%xmm5,%%xmm1		\n\t"\
		"subpd	%%xmm4,%%xmm2		\n\t"\
		"subpd	%%xmm5,%%xmm3		\n\t"\
		"/* Do	the p4,12 combo: */	\n\t"\
		"addq	$0x60,%%rsi 	/* c14 */\n\t"\
		"movaps	    (%%rdx),%%xmm4	\n\t"\
		"movaps	0x10(%%rdx),%%xmm5	\n\t"\
		"movaps	    (%%rdx),%%xmm6	\n\t"\
		"movaps	0x10(%%rdx),%%xmm7	\n\t"\
		"mulpd	    (%%rsi),%%xmm4	\n\t"\
		"mulpd	    (%%rsi),%%xmm5	\n\t"\
		"mulpd	0x10(%%rsi),%%xmm6	\n\t"\
		"mulpd	0x10(%%rsi),%%xmm7	\n\t"\
		"subq	$0x290,%%rsi 	/* r9 */\n\t"\
		"addpd	%%xmm6,%%xmm5		\n\t"\
		"subpd	%%xmm7,%%xmm4		\n\t"\
		"movaps	%%xmm5,0x10(%%rsi)	/* r9.im */\n\t"\
		"movaps	%%xmm4,    (%%rsi)	/* r9.re */\n\t"\
		"movaps	    (%%rbx),%%xmm4	\n\t"\
		"movaps	0x10(%%rbx),%%xmm5	\n\t"\
		"movaps	    (%%rbx),%%xmm6	\n\t"\
		"movaps	0x10(%%rbx),%%xmm7	\n\t"\
		"mulpd	0x270(%%rsi),%%xmm4	/* c6, using r9 as base ptr */\n\t"\
		"mulpd	0x270(%%rsi),%%xmm5	\n\t"\
		"mulpd	0x280(%%rsi),%%xmm6	\n\t"\
		"mulpd	0x280(%%rsi),%%xmm7	\n\t"\
		"addpd	%%xmm6,%%xmm5		\n\t"\
		"subpd	%%xmm7,%%xmm4		\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t"\
		"movaps	%%xmm4,%%xmm6		\n\t"\
		"subpd	     (%%rsi),%%xmm4	\n\t"\
		"subpd	0x010(%%rsi),%%xmm5	\n\t"\
		"addpd	     (%%rsi),%%xmm6	\n\t"\
		"addpd	0x010(%%rsi),%%xmm7	\n\t"\
		"/* Finish radix-4 butterfly and store results into temporary-array slots: */\n\t"\
		"subpd	%%xmm6,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm2		\n\t"\
		"subpd	%%xmm7,%%xmm1		\n\t"\
		"subpd	%%xmm4,%%xmm3		\n\t"\
		"movaps	%%xmm0,0x040(%%rsi)	\n\t"\
		"movaps	%%xmm2,0x020(%%rsi)	\n\t"\
		"movaps	%%xmm1,0x050(%%rsi)	\n\t"\
		"movaps	%%xmm3,0x070(%%rsi)	\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm0,%%xmm6		\n\t"\
		"addpd	%%xmm2,%%xmm5		\n\t"\
		"addpd	%%xmm1,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm4		\n\t"\
		"movaps	%%xmm6,     (%%rsi)	\n\t"\
		"movaps	%%xmm5,0x060(%%rsi)	\n\t"\
		"movaps	%%xmm7,0x010(%%rsi)	\n\t"\
		"movaps	%%xmm4,0x030(%%rsi)	\n\t"\
		"subq	%%rdi,%%rax			\n\t"\
		"subq	%%rdi,%%rbx			\n\t"\
		"subq	%%rdi,%%rcx			\n\t"\
		"subq	%%rdi,%%rdx			\n\t"\
		"movslq	%[__p1],%%rdi		\n\t"\
		"shlq	$3,%%rdi			\n\t"\
		"addq	%%rdi,%%rax			\n\t"\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	%%rdi,%%rcx			\n\t"\
		"addq	%%rdi,%%rdx			\n\t"\
		"/* SSE2_RADIX4_DIF_4TWIDDLE_B(r17,c1) */\n\t"\
		"/* Do	the p0,p8 combo: */	\n\t"\
		"addq	$0x2b0,%%rsi 	/* c1, from r9 */\n\t"\
		"movaps	    (%%rax),%%xmm0	\n\t"\
		"movaps	    (%%rcx),%%xmm4	\n\t"\
		"movaps	0x10(%%rax),%%xmm1	\n\t"\
		"movaps	0x10(%%rcx),%%xmm5	\n\t"\
		"movaps	    (%%rsi),%%xmm6	\n\t"\
		"movaps	0x10(%%rsi),%%xmm7	\n\t"\
		"movaps	%%xmm0,%%xmm2		\n\t"\
		"movaps	%%xmm1,%%xmm3		\n\t"\
		"mulpd	%%xmm6,%%xmm0		\n\t"\
		"mulpd	%%xmm6,%%xmm1		\n\t"\
		"mulpd	%%xmm7,%%xmm2		\n\t"\
		"mulpd	%%xmm7,%%xmm3		\n\t"\
		"movaps	%%xmm4,%%xmm6		\n\t"\
		"addpd	%%xmm2,%%xmm1		\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t"\
		"mulpd	0x20(%%rsi) ,%%xmm4	\n\t"\
		"subpd	%%xmm3,%%xmm0		\n\t"\
		"mulpd	0x20(%%rsi),%%xmm5	\n\t"\
		"mulpd	0x30(%%rsi),%%xmm6	\n\t"\
		"movaps	%%xmm0,%%xmm2		\n\t"\
		"mulpd	0x30(%%rsi),%%xmm7	\n\t"\
		"addpd	%%xmm6,%%xmm5		\n\t"\
		"movaps	%%xmm1,%%xmm3		\n\t"\
		"subpd	%%xmm7,%%xmm4		\n\t"\
		"addpd	%%xmm4,%%xmm0		\n\t"\
		"addpd	%%xmm5,%%xmm1		\n\t"\
		"subpd	%%xmm4,%%xmm2		\n\t"\
		"subpd	%%xmm5,%%xmm3		\n\t"\
		"/* Do	the p4,12 combo: */	\n\t"\
		"addq	$0x60,%%rsi 	/* c13 */\n\t"\
		"movaps	    (%%rdx),%%xmm4	\n\t"\
		"movaps	0x10(%%rdx),%%xmm5	\n\t"\
		"movaps	    (%%rdx),%%xmm6	\n\t"\
		"movaps	0x10(%%rdx),%%xmm7	\n\t"\
		"mulpd	    (%%rsi),%%xmm4	\n\t"\
		"mulpd	    (%%rsi),%%xmm5	\n\t"\
		"mulpd	0x10(%%rsi),%%xmm6	\n\t"\
		"mulpd	0x10(%%rsi),%%xmm7	\n\t"\
		"subq	$0x290,%%rsi 	/* r17 */\n\t"\
		"addpd	%%xmm6,%%xmm5		\n\t"\
		"subpd	%%xmm7,%%xmm4		\n\t"\
		"movaps	%%xmm5,0x010(%%rsi)	\n\t"\
		"movaps	%%xmm4,     (%%rsi)	\n\t"\
		"movaps	    (%%rbx),%%xmm4	\n\t"\
		"movaps	0x10(%%rbx),%%xmm5	\n\t"\
		"movaps	    (%%rbx),%%xmm6	\n\t"\
		"movaps	0x10(%%rbx),%%xmm7	\n\t"\
		"mulpd	0x270(%%rsi),%%xmm4	/* c5, using r17 as base ptr */\n\t"\
		"mulpd	0x270(%%rsi),%%xmm5	\n\t"\
		"mulpd	0x280(%%rsi),%%xmm6	\n\t"\
		"mulpd	0x280(%%rsi),%%xmm7	\n\t"\
		"addpd	%%xmm6,%%xmm5		\n\t"\
		"subpd	%%xmm7,%%xmm4		\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t"\
		"movaps	%%xmm4,%%xmm6		\n\t"\
		"subpd	     (%%rsi),%%xmm4	\n\t"\
		"subpd	0x010(%%rsi),%%xmm5	\n\t"\
		"addpd	     (%%rsi),%%xmm6	\n\t"\
		"addpd	0x010(%%rsi),%%xmm7	\n\t"\
		"/* Finish radix-4 butterfly and store results into temporary-array slots: */\n\t"\
		"subpd	%%xmm6,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm2		\n\t"\
		"subpd	%%xmm7,%%xmm1		\n\t"\
		"subpd	%%xmm4,%%xmm3		\n\t"\
		"movaps	%%xmm0,0x040(%%rsi)	\n\t"\
		"movaps	%%xmm2,0x020(%%rsi)	\n\t"\
		"movaps	%%xmm1,0x050(%%rsi)	\n\t"\
		"movaps	%%xmm3,0x070(%%rsi)	\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm0,%%xmm6		\n\t"\
		"addpd	%%xmm2,%%xmm5		\n\t"\
		"addpd	%%xmm1,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm4		\n\t"\
		"movaps	%%xmm6,     (%%rsi)	\n\t"\
		"movaps	%%xmm5,0x060(%%rsi)	\n\t"\
		"movaps	%%xmm7,0x010(%%rsi)	\n\t"\
		"movaps	%%xmm4,0x030(%%rsi)	\n\t"\
		"subq	%%rdi,%%rax			\n\t"\
		"subq	%%rdi,%%rbx			\n\t"\
		"subq	%%rdi,%%rcx			\n\t"\
		"subq	%%rdi,%%rdx			\n\t"\
		"movslq	%[__p3],%%rdi		\n\t"\
		"shlq	$3,%%rdi			\n\t"\
		"addq	%%rdi,%%rax			\n\t"\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	%%rdi,%%rcx			\n\t"\
		"addq	%%rdi,%%rdx			\n\t"\
		"/* SSE2_RADIX4_DIF_4TWIDDLE_B(r25,c3) */\n\t"\
		"/* Do	the p0,p8 combo: */	\n\t"\
		"addq	$0x2b0,%%rsi 	/* c3, from r17 */\n\t"\
		"movaps	    (%%rax),%%xmm0	\n\t"\
		"movaps	    (%%rcx),%%xmm4	\n\t"\
		"movaps	0x10(%%rax),%%xmm1	\n\t"\
		"movaps	0x10(%%rcx),%%xmm5	\n\t"\
		"movaps	    (%%rsi),%%xmm6	\n\t"\
		"movaps	0x10(%%rsi),%%xmm7	\n\t"\
		"movaps	%%xmm0,%%xmm2		\n\t"\
		"movaps	%%xmm1,%%xmm3		\n\t"\
		"mulpd	%%xmm6,%%xmm0		\n\t"\
		"mulpd	%%xmm6,%%xmm1		\n\t"\
		"mulpd	%%xmm7,%%xmm2		\n\t"\
		"mulpd	%%xmm7,%%xmm3		\n\t"\
		"movaps	%%xmm4,%%xmm6		\n\t"\
		"addpd	%%xmm2,%%xmm1		\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t"\
		"mulpd	0x20(%%rsi),%%xmm4	\n\t"\
		"subpd	%%xmm3,%%xmm0		\n\t"\
		"mulpd	0x20(%%rsi),%%xmm5	\n\t"\
		"mulpd	0x30(%%rsi),%%xmm6	\n\t"\
		"movaps	%%xmm0,%%xmm2		\n\t"\
		"mulpd	0x30(%%rsi),%%xmm7	\n\t"\
		"addpd	%%xmm6,%%xmm5		\n\t"\
		"movaps	%%xmm1,%%xmm3		\n\t"\
		"subpd	%%xmm7,%%xmm4		\n\t"\
		"addpd	%%xmm4,%%xmm0		\n\t"\
		"addpd	%%xmm5,%%xmm1		\n\t"\
		"subpd	%%xmm4,%%xmm2		\n\t"\
		"subpd	%%xmm5,%%xmm3		\n\t"\
		"/* Do	the p4,12 combo: */	\n\t"\
		"addq	$0x60,%%rsi 	/* c15 */\n\t"\
		"movaps	    (%%rdx),%%xmm4	\n\t"\
		"movaps	0x10(%%rdx),%%xmm5	\n\t"\
		"movaps	    (%%rdx),%%xmm6	\n\t"\
		"movaps	0x10(%%rdx),%%xmm7	\n\t"\
		"mulpd	    (%%rsi),%%xmm4	\n\t"\
		"mulpd	    (%%rsi),%%xmm5	\n\t"\
		"mulpd	0x10(%%rsi),%%xmm6	\n\t"\
		"mulpd	0x10(%%rsi),%%xmm7	\n\t"\
		"subq	$0x290,%%rsi 	/* r25 */\n\t"\
		"addpd	%%xmm6,%%xmm5		\n\t"\
		"subpd	%%xmm7,%%xmm4		\n\t"\
		"movaps	%%xmm5,0x010(%%rsi)	\n\t"\
		"movaps	%%xmm4,     (%%rsi)	\n\t"\
		"movaps	    (%%rbx),%%xmm4	\n\t"\
		"movaps	0x10(%%rbx),%%xmm5	\n\t"\
		"movaps	    (%%rbx),%%xmm6	\n\t"\
		"movaps	0x10(%%rbx),%%xmm7	\n\t"\
		"mulpd	0x270(%%rsi),%%xmm4	/* c7, using r25 as base ptr */\n\t"\
		"mulpd	0x270(%%rsi),%%xmm5	\n\t"\
		"mulpd	0x280(%%rsi),%%xmm6	\n\t"\
		"mulpd	0x280(%%rsi),%%xmm7	\n\t"\
		"addpd	%%xmm6,%%xmm5		\n\t"\
		"subpd	%%xmm7,%%xmm4		\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t"\
		"movaps	%%xmm4,%%xmm6		\n\t"\
		"subpd	     (%%rsi),%%xmm4	\n\t"\
		"subpd	0x010(%%rsi),%%xmm5	\n\t"\
		"addpd	     (%%rsi),%%xmm6	\n\t"\
		"addpd	0x010(%%rsi),%%xmm7	\n\t"\
		"/* Finish radix-4 butterfly and store results into temporary-array slots: */\n\t"\
		"subpd	%%xmm6,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm2		\n\t"\
		"subpd	%%xmm7,%%xmm1		\n\t"\
		"subpd	%%xmm4,%%xmm3		\n\t"\
		"movaps	%%xmm0,0x040(%%rsi)	\n\t"\
		"movaps	%%xmm2,0x020(%%rsi)	\n\t"\
		"movaps	%%xmm1,0x050(%%rsi)	\n\t"\
		"movaps	%%xmm3,0x070(%%rsi)	\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm0,%%xmm6		\n\t"\
		"addpd	%%xmm2,%%xmm5		\n\t"\
		"addpd	%%xmm1,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm4		\n\t"\
		"movaps	%%xmm6,     (%%rsi)	\n\t"\
		"movaps	%%xmm5,0x060(%%rsi)	\n\t"\
		"movaps	%%xmm7,0x010(%%rsi)	\n\t"\
		"movaps	%%xmm4,0x030(%%rsi)	\n\t"\
		"/*************************************************************************************/\n\t"\
		"/*  And now do four more radix-4 transforms, including the internal twiddle factors: */\n\t"\
		"/*************************************************************************************/\n\t"\
		"/* Block 1: */				\n\t"\
		"movq	%[__add0],%%rax		\n\t"\
		"movslq	%[__p1],%%rbx		\n\t"\
		"movslq	%[__p2],%%rcx		\n\t"\
		"shlq	$3,%%rbx			\n\t"\
		"shlq	$3,%%rcx			\n\t"\
		"movq	%%rdi,%%rdx	/* p3-ptr already in rdi */\n\t"\
		"addq	%%rax,%%rbx			\n\t"\
		"addq	%%rax,%%rcx			\n\t"\
		"addq	%%rax,%%rdx			\n\t"\
		"movq	%[__r1],%%rsi		\n\t"\
		"movaps	     (%%rsi),%%xmm0	\n\t"\
		"movaps	0x100(%%rsi),%%xmm4	\n\t"\
		"movaps	0x010(%%rsi),%%xmm1	\n\t"\
		"movaps	0x110(%%rsi),%%xmm5	\n\t"\
		"movaps	0x080(%%rsi),%%xmm2	\n\t"\
		"movaps	0x180(%%rsi),%%xmm6	\n\t"\
		"movaps	0x090(%%rsi),%%xmm3	\n\t"\
		"movaps	0x190(%%rsi),%%xmm7	\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm6,%%xmm4		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"subpd	%%xmm7,%%xmm5		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm4,%%xmm6		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"addpd	%%xmm5,%%xmm7		\n\t"\
		"subpd	%%xmm6,%%xmm2		\n\t"\
		"subpd	%%xmm7,%%xmm3		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"movaps	%%xmm2,    (%%rbx)	\n\t"\
		"movaps	%%xmm3,0x10(%%rbx)	\n\t"\
		"addpd	%%xmm2,	%%xmm6		\n\t"\
		"addpd	%%xmm3,	%%xmm7		\n\t"\
		"movaps	%%xmm6,    (%%rax)	\n\t"\
		"movaps	%%xmm7,0x10(%%rax)	\n\t"\
		"subpd	%%xmm5,	%%xmm0		\n\t"\
		"subpd	%%xmm4,	%%xmm1		\n\t"\
		"addpd	%%xmm5,	%%xmm5		\n\t"\
		"addpd	%%xmm4,	%%xmm4		\n\t"\
		"movaps	%%xmm0,    (%%rcx)	\n\t"\
		"movaps	%%xmm1,0x10(%%rdx)	\n\t"\
		"addpd	%%xmm0,	%%xmm5		\n\t"\
		"addpd	%%xmm1,	%%xmm4		\n\t"\
		"movaps	%%xmm5,    (%%rdx)	\n\t"\
		"movaps	%%xmm4,0x10(%%rcx)	\n\t"\
		"/* Block 3: */				\n\t"\
		"movslq	%[__p4],%%rdi		\n\t"\
		"shlq	$3,%%rdi			\n\t"\
		"movaps	0x200(%%rsi),%%xmm3	/* isrt2, using r1 as base-ptr */\n\t"\
		"addq	%%rdi,%%rax			\n\t"\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	%%rdi,%%rcx			\n\t"\
		"addq	%%rdi,%%rdx			\n\t"\
		"addq	$0x40,%%rsi		/* r5 */\n\t"\
		"movaps	0x100(%%rsi),%%xmm4	/* r5 */\n\t"\
		"movaps	0x110(%%rsi),%%xmm5	\n\t"\
		"movaps	0x180(%%rsi),%%xmm6	\n\t"\
		"movaps	0x190(%%rsi),%%xmm7	\n\t"\
		"mulpd	%%xmm3,%%xmm4		\n\t"\
		"movaps	     (%%rsi),%%xmm0	\n\t"\
		"mulpd	%%xmm3,%%xmm5		\n\t"\
		"movaps	0x010(%%rsi),%%xmm1	\n\t"\
		"mulpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	0x080(%%rsi),%%xmm2	\n\t"\
		"mulpd	%%xmm3,%%xmm7		\n\t"\
		"movaps	0x090(%%rsi),%%xmm3	\n\t"\
		"subpd	%%xmm3,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm4		\n\t"\
		"subpd	%%xmm2,%%xmm1		\n\t"\
		"subpd	%%xmm6,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm3		\n\t"\
		"addpd	%%xmm4,%%xmm5		\n\t"\
		"addpd	%%xmm1,%%xmm2		\n\t"\
		"addpd	%%xmm7,%%xmm6		\n\t"\
		"subpd	%%xmm6,%%xmm4		\n\t"\
		"subpd	%%xmm7,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm6		\n\t"\
		"addpd	%%xmm5,%%xmm7		\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm3		\n\t"\
		"subpd	%%xmm5,%%xmm2		\n\t"\
		"subpd	%%xmm6,%%xmm1		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"movaps	%%xmm0,    (%%rbx)	\n\t"\
		"movaps	%%xmm3,    (%%rcx)	\n\t"\
		"movaps	%%xmm2,0x10(%%rbx)	\n\t"\
		"movaps	%%xmm1,0x10(%%rdx)	\n\t"\
		"addpd	%%xmm0,	%%xmm4		\n\t"\
		"addpd	%%xmm3,	%%xmm7		\n\t"\
		"addpd	%%xmm2,	%%xmm5		\n\t"\
		"addpd	%%xmm1,	%%xmm6		\n\t"\
		"movaps	%%xmm4,    (%%rax)	\n\t"\
		"movaps	%%xmm7,    (%%rdx)	\n\t"\
		"movaps	%%xmm5,0x10(%%rax)	\n\t"\
		"movaps	%%xmm6,0x10(%%rcx)	\n\t"\
		"/* Block 2: */				\n\t"\
		"addq	%%rdi,%%rax			\n\t"\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	%%rdi,%%rcx			\n\t"\
		"addq	%%rdi,%%rdx			\n\t"\
		"subq	$0x20,%%rsi		/* r3 */\n\t"\
		"movaps	0x100(%%rsi),%%xmm4	\n\t"\
		"movaps	0x110(%%rsi),%%xmm5	\n\t"\
		"movaps	0x1f0(%%rsi),%%xmm3	/* cc0, using r3 as base-ptr */\n\t"\
		"movaps	0x200(%%rsi),%%xmm2	/* ss0, using r3 as base-ptr */\n\t"\
		"movaps	%%xmm4,%%xmm6		\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t"\
		"mulpd	%%xmm3,%%xmm4		\n\t"\
		"mulpd	%%xmm3,%%xmm5		\n\t"\
		"mulpd	%%xmm2,%%xmm6		\n\t"\
		"movaps	0x180(%%rsi),%%xmm0	\n\t"\
		"mulpd	%%xmm2,%%xmm7		\n\t"\
		"movaps	0x190(%%rsi),%%xmm1	\n\t"\
		"addpd	%%xmm6,%%xmm5		\n\t"\
		"movaps	%%xmm0,%%xmm6		\n\t"\
		"subpd	%%xmm7,%%xmm4		\n\t"\
		"movaps	%%xmm1,%%xmm7		\n\t"\
		"mulpd	%%xmm2,%%xmm6		\n\t"\
		"mulpd	%%xmm2,%%xmm7		\n\t"\
		"mulpd	%%xmm3,%%xmm0		\n\t"\
		"mulpd	%%xmm3,%%xmm1		\n\t"\
		"addpd	%%xmm0,%%xmm7		\n\t"\
		"subpd	%%xmm1,%%xmm6		\n\t"\
		"movaps	%%xmm4,%%xmm2		\n\t"\
		"movaps	%%xmm5,%%xmm3		\n\t"\
		"subpd	%%xmm6,%%xmm4		\n\t"\
		"subpd	%%xmm7,%%xmm5		\n\t"\
		"addpd	%%xmm2,%%xmm6		\n\t"\
		"addpd	%%xmm3,%%xmm7		\n\t"\
		"movaps	0x080(%%rsi),%%xmm2	\n\t"\
		"movaps	0x090(%%rsi),%%xmm3	\n\t"\
		"movaps	0x1e0(%%rsi),%%xmm1	/* isrt2, using r3 as base-ptr */\n\t"\
		"movaps	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm2		\n\t"\
		"addpd	%%xmm0,%%xmm3		\n\t"\
		"mulpd	%%xmm1,%%xmm2		\n\t"\
		"mulpd	%%xmm1,%%xmm3		\n\t"\
		"movaps	     (%%rsi),%%xmm0	\n\t"\
		"movaps	0x010(%%rsi),%%xmm1	\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"subpd	%%xmm6,%%xmm2		\n\t"\
		"subpd	%%xmm5,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm3		\n\t"\
		"subpd	%%xmm4,%%xmm1		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"movaps	%%xmm2,    (%%rbx)	\n\t"\
		"movaps	%%xmm0,    (%%rcx)	\n\t"\
		"movaps	%%xmm3,0x10(%%rbx)	\n\t"\
		"movaps	%%xmm1,0x10(%%rdx)	\n\t"\
		"addpd	%%xmm2,	%%xmm6		\n\t"\
		"addpd	%%xmm0,	%%xmm5		\n\t"\
		"addpd	%%xmm3,	%%xmm7		\n\t"\
		"addpd	%%xmm1,	%%xmm4		\n\t"\
		"movaps	%%xmm6,    (%%rax)	\n\t"\
		"movaps	%%xmm5,    (%%rdx)	\n\t"\
		"movaps	%%xmm7,0x10(%%rax)	\n\t"\
		"movaps	%%xmm4,0x10(%%rcx)	\n\t"\
		"/* Block 4: */				\n\t"\
		"addq	%%rdi,%%rax			\n\t"\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	%%rdi,%%rcx			\n\t"\
		"addq	%%rdi,%%rdx			\n\t"\
		"addq	$0x40,%%rsi		/* r7 */\n\t"\
		"movaps	0x100(%%rsi),%%xmm4	\n\t"\
		"movaps	0x110(%%rsi),%%xmm5	\n\t"\
		"movaps	0x1b0(%%rsi),%%xmm2	/* cc0, using r7 as base-ptr */\n\t"\
		"movaps	0x1c0(%%rsi),%%xmm3	/* ss0, using r7 as base-ptr */\n\t"\
		"movaps	%%xmm4,%%xmm6		\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t"\
		"mulpd	%%xmm3,%%xmm4		\n\t"\
		"mulpd	%%xmm3,%%xmm5		\n\t"\
		"mulpd	%%xmm2,%%xmm6		\n\t"\
		"movaps	0x180(%%rsi),%%xmm0	\n\t"\
		"mulpd	%%xmm2,%%xmm7		\n\t"\
		"movaps	0x190(%%rsi),%%xmm1	\n\t"\
		"addpd	%%xmm6,%%xmm5		\n\t"\
		"movaps	%%xmm0,%%xmm6		\n\t"\
		"subpd	%%xmm7,%%xmm4		\n\t"\
		"movaps	%%xmm1,%%xmm7		\n\t"\
		"mulpd	%%xmm2,%%xmm6		\n\t"\
		"mulpd	%%xmm2,%%xmm7		\n\t"\
		"mulpd	%%xmm3,%%xmm0		\n\t"\
		"mulpd	%%xmm3,%%xmm1		\n\t"\
		"addpd	%%xmm0,%%xmm7		\n\t"\
		"subpd	%%xmm1,%%xmm6		\n\t"\
		"movaps	%%xmm4,%%xmm2		\n\t"\
		"movaps	%%xmm5,%%xmm3		\n\t"\
		"subpd	%%xmm6,%%xmm4		\n\t"\
		"subpd	%%xmm7,%%xmm5		\n\t"\
		"addpd	%%xmm2,%%xmm6		\n\t"\
		"addpd	%%xmm3,%%xmm7		\n\t"\
		"movaps	0x080(%%rsi),%%xmm2	\n\t"\
		"movaps	0x090(%%rsi),%%xmm3	\n\t"\
		"movaps	0x1a0(%%rsi),%%xmm1		/* isrt2, using r7 as base-ptr */\n\t"\
		"movaps	%%xmm2,%%xmm0		\n\t"\
		"addpd	%%xmm3,%%xmm2		\n\t"\
		"subpd	%%xmm0,%%xmm3		\n\t"\
		"mulpd	%%xmm1,%%xmm2		\n\t"\
		"mulpd	%%xmm1,%%xmm3		\n\t"\
		"movaps	     (%%rsi),%%xmm0	\n\t"\
		"movaps	0x010(%%rsi),%%xmm1	\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm2		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"subpd	%%xmm6,%%xmm3		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"movaps	%%xmm0,    (%%rbx)	\n\t"\
		"movaps	%%xmm2,    (%%rcx)	\n\t"\
		"movaps	%%xmm1,0x10(%%rbx)	\n\t"\
		"movaps	%%xmm3,0x10(%%rdx)	\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"addpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	%%xmm4,    (%%rax)	\n\t"\
		"movaps	%%xmm7,    (%%rdx)	\n\t"\
		"movaps	%%xmm5,0x10(%%rax)	\n\t"\
		"movaps	%%xmm6,0x10(%%rcx)	\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__p1] "m" (Xp1)\
		 ,[__p2] "m" (Xp2)\
		 ,[__p3] "m" (Xp3)\
		 ,[__p4] "m" (Xp4)\
		 ,[__p8] "m" (Xp8)\
		 ,[__p12] "m" (Xp12)\
		 ,[__r1] "m" (Xr1)\
		 ,[__isrt2] "m" (Xisrt2)\
		: "rax","rbx","rcx","rdx","rdi","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
	);\
	}

	#define SSE2_RADIX16_DIT_TWIDDLE(Xadd0,Xp1,Xp2,Xp3,Xp4,Xp8,Xr1,Xisrt2)\
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
		"/* SSE2_RADIX4_DIT_0TWIDDLE_B(r1) */\n\t"\
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
		"movq	%[__r1],%%rsi		\n\t"\
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
		"/* SSE2_RADIX4_DIT_0TWIDDLE_B(r9) */\n\t"\
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
		"addq	$0x80,%%rsi	/* r9  */\n\t"\
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
		"/* SSE2_RADIX4_DIT_0TWIDDLE_B(r17) */\n\t"\
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
		"addq	$0x80,%%rsi	/* r17 */\n\t"\
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
		"/* SSE2_RADIX4_DIT_0TWIDDLE_B(r25) */\n\t"\
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
		"addq	$0x80,%%rsi	/* r25 */\n\t"\
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
		"/****************************************************************************************************\n\t"\
		"/*..and now do four more radix-4 transforms, including the internal and external twiddle factors:  */\n\t"\
		"/***************************************************************************************************/\n\t"\
	"/*...Block 1: 24 MOVapd, 26 ADD/SUBpd, 12 MULpd */\n\t"\
		"movq	%[__add0],%%rax		\n\t"\
		"movslq	%[__p4],%%rbx		\n\t"\
		"movq	%[__r1],%%rcx		\n\t"\
		"leaq	0x080(%%rcx),%%rdx	/* r9 */\n\t"\
		"movslq	%[__p8],%%rdi		\n\t"\
		"shlq	$3,%%rbx			\n\t"\
		"shlq	$3,%%rdi			\n\t"\
		"addq	%%rax,%%rbx			\n\t"\
		"movaps	     (%%rdx),%%xmm2	\n\t"\
		"movaps	0x100(%%rdx),%%xmm4	\n\t"\
		"movaps	0x010(%%rdx),%%xmm3	\n\t"\
		"movaps	0x110(%%rdx),%%xmm5	\n\t"\
		"movaps	     (%%rcx),%%xmm0	\n\t"\
		"movaps	0x100(%%rcx),%%xmm6	\n\t"\
		"movaps	0x010(%%rcx),%%xmm1	\n\t"\
		"movaps	0x110(%%rcx),%%xmm7	\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm4,%%xmm6		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"subpd	%%xmm5,%%xmm7		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm6,%%xmm4		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"addpd	%%xmm7,%%xmm5		\n\t"\
		"/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c0) */\n\t"\
		"addq	$0xb0,%%rsi /* c0, incr from r25 */\n\t"\
		"subpd	%%xmm4,%%xmm2		\n\t"\
		"subpd	%%xmm7,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm3		\n\t"\
		"subpd	%%xmm6,%%xmm1		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm2,%%xmm4		\n\t"\
		"addpd	%%xmm0,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm5		\n\t"\
		"addpd	%%xmm1,%%xmm6		\n\t"\
		"movaps	%%xmm2,     (%%rcx)	\n\t"\
		"movaps	%%xmm0,0x010(%%rdx)	\n\t"\
		"movaps	%%xmm3,0x010(%%rcx)	\n\t"\
		"movaps	%%xmm6,     (%%rdx)	\n\t"\
		"movaps	%%xmm4,%%xmm2		\n\t"\
		"movaps	%%xmm7,%%xmm0		\n\t"\
		"movaps	%%xmm5,%%xmm3		\n\t"\
		"movaps	%%xmm1,%%xmm6		\n\t"\
		"mulpd	    (%%rsi),%%xmm4	\n\t"\
		"mulpd	0x20(%%rsi),%%xmm7	\n\t"\
		"mulpd	    (%%rsi),%%xmm5	\n\t"\
		"mulpd	0x20(%%rsi),%%xmm1	\n\t"\
		"mulpd	0x10(%%rsi),%%xmm2	\n\t"\
		"mulpd	0x30(%%rsi),%%xmm0	\n\t"\
		"mulpd	0x10(%%rsi),%%xmm3	\n\t"\
		"mulpd	0x30(%%rsi),%%xmm6	\n\t"\
		"subpd	%%xmm2,%%xmm5		\n\t"\
		"subpd	%%xmm0,%%xmm1		\n\t"\
		"addpd	%%xmm3,%%xmm4		\n\t"\
		"addpd	%%xmm6,%%xmm7		\n\t"\
		"movaps	%%xmm5,0x10(%%rax)	\n\t"\
		"movaps	%%xmm1,0x10(%%rbx)	\n\t"\
		"movaps	%%xmm4,    (%%rax)	\n\t"\
		"movaps	%%xmm7,    (%%rbx)	\n\t"\
		"addq	%%rdi,%%rax			\n\t"\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	$0x40,%%rsi	/* c8 */\n\t"\
		"movaps	     (%%rcx),%%xmm4	\n\t"\
		"movaps	0x010(%%rdx),%%xmm0	\n\t"\
		"movaps	0x010(%%rcx),%%xmm5	\n\t"\
		"movaps	     (%%rdx),%%xmm6	\n\t"\
		"movaps	%%xmm4,%%xmm2		\n\t"\
		"movaps	%%xmm0,%%xmm1		\n\t"\
		"movaps	%%xmm5,%%xmm3		\n\t"\
		"movaps	%%xmm6,%%xmm7		\n\t"\
		"mulpd	    (%%rsi),%%xmm4	\n\t"\
		"mulpd	0x20(%%rsi),%%xmm0	\n\t"\
		"mulpd	    (%%rsi),%%xmm5	\n\t"\
		"mulpd	0x20(%%rsi),%%xmm6	\n\t"\
		"mulpd	0x10(%%rsi),%%xmm2	\n\t"\
		"mulpd	0x30(%%rsi),%%xmm1	\n\t"\
		"mulpd	0x10(%%rsi),%%xmm3	\n\t"\
		"mulpd	0x30(%%rsi),%%xmm7	\n\t"\
		"subpd	%%xmm2,%%xmm5		\n\t"\
		"subpd	%%xmm1,%%xmm6		\n\t"\
		"addpd	%%xmm3,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm0		\n\t"\
		"movaps	%%xmm5,0x10(%%rax)	\n\t"\
		"movaps	%%xmm6,0x10(%%rbx)	\n\t"\
		"movaps	%%xmm4,    (%%rax)	\n\t"\
		"movaps	%%xmm0,    (%%rbx)	\n\t"\
	"/*...Block 2: 34 MOVapd, 36 ADD/SUBpd, 26 MULpd */\n\t"\
		"movq	%[__add0],%%rax		\n\t"\
		"movslq	%[__p1],%%rsi		\n\t"\
		"movslq	%[__p4],%%rbx		\n\t"\
		"shlq	$3,%%rsi			\n\t"\
		"shlq	$3,%%rbx			\n\t"\
		"addq	%%rsi,%%rax			\n\t"\
		"addq	%%rax,%%rbx			\n\t"\
		"addq	$0x120,%%rcx	/* r3 +0x100 */\n\t"\
		"addq	$0x120,%%rdx	/* r11+0x100 */\n\t"\
		"movq	%[__isrt2],%%rsi		\n\t"\
		"movaps	     (%%rcx),%%xmm4		\n\t"\
		"movaps	0x010(%%rcx),%%xmm5	\n\t"\
		"movaps	0x10(%%rsi),%%xmm2	/* cc0 */\n\t"\
		"movaps	0x20(%%rsi),%%xmm3	/* ss0 */\n\t"\
		"movaps	%%xmm4,%%xmm6		\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t"\
		"mulpd	%%xmm2,%%xmm4		\n\t"\
		"mulpd	%%xmm2,%%xmm5		\n\t"\
		"mulpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	     (%%rdx),%%xmm0	\n\t"\
		"mulpd	%%xmm3,%%xmm7		\n\t"\
		"movaps	0x010(%%rdx),%%xmm1	\n\t"\
		"subpd	%%xmm6,%%xmm5		\n\t"\
		"movaps	%%xmm0,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm4		\n\t"\
		"movaps	%%xmm1,%%xmm7		\n\t"\
		"mulpd	%%xmm3,%%xmm0		\n\t"\
		"mulpd	%%xmm3,%%xmm1		\n\t"\
		"mulpd	%%xmm2,%%xmm6		\n\t"\
		"mulpd	%%xmm2,%%xmm7		\n\t"\
		"subpd	%%xmm6,%%xmm1		\n\t"\
		"addpd	%%xmm7,%%xmm0		\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t"\
		"movaps	%%xmm4,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"subpd	%%xmm0,%%xmm6		\n\t"\
		"subpd	%%xmm1,%%xmm7		\n\t"\
		"subq	$0x100,%%rcx		\n\t"\
		"subq	$0x100,%%rdx		\n\t"\
		"movaps	     (%%rdx),%%xmm2	\n\t"\
		"movaps	0x010(%%rdx),%%xmm3	\n\t"\
		"movaps	(%%rsi),%%xmm1		\n\t"\
		"movaps	%%xmm3,%%xmm0		\n\t"\
		"subpd	%%xmm2,%%xmm3		\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"mulpd	%%xmm1,%%xmm2		\n\t"\
		"mulpd	%%xmm1,%%xmm3		\n\t"\
		"movaps	     (%%rcx),%%xmm0	\n\t"\
		"movaps	0x010(%%rcx),%%xmm1	\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c1) */\n\t"\
		"addq	$0x130,%%rsi	/* c1 */\n\t"\
		"subpd	%%xmm4,%%xmm2		\n\t"\
		"subpd	%%xmm7,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm3		\n\t"\
		"subpd	%%xmm6,%%xmm1		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm2,%%xmm4		\n\t"\
		"addpd	%%xmm0,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm5		\n\t"\
		"addpd	%%xmm1,%%xmm6		\n\t"\
		"movaps	%%xmm2,     (%%rcx)	\n\t"\
		"movaps	%%xmm0,0x010(%%rdx)	\n\t"\
		"movaps	%%xmm3,0x010(%%rcx)	\n\t"\
		"movaps	%%xmm6,     (%%rdx)	\n\t"\
		"movaps	%%xmm4,%%xmm2		\n\t"\
		"movaps	%%xmm7,%%xmm0		\n\t"\
		"movaps	%%xmm5,%%xmm3		\n\t"\
		"movaps	%%xmm1,%%xmm6		\n\t"\
		"mulpd	    (%%rsi),%%xmm4	\n\t"\
		"mulpd	0x20(%%rsi),%%xmm7	\n\t"\
		"mulpd	    (%%rsi),%%xmm5	\n\t"\
		"mulpd	0x20(%%rsi),%%xmm1	\n\t"\
		"mulpd	0x10(%%rsi),%%xmm2	\n\t"\
		"mulpd	0x30(%%rsi),%%xmm0	\n\t"\
		"mulpd	0x10(%%rsi),%%xmm3	\n\t"\
		"mulpd	0x30(%%rsi),%%xmm6	\n\t"\
		"subpd	%%xmm2,%%xmm5		\n\t"\
		"subpd	%%xmm0,%%xmm1		\n\t"\
		"addpd	%%xmm3,%%xmm4		\n\t"\
		"addpd	%%xmm6,%%xmm7		\n\t"\
		"movaps	%%xmm5,0x10(%%rax)	\n\t"\
		"movaps	%%xmm1,0x10(%%rbx)	\n\t"\
		"movaps	%%xmm4,    (%%rax)	\n\t"\
		"movaps	%%xmm7,    (%%rbx)	\n\t"\
		"addq	%%rdi,%%rax			\n\t"\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	$0x40,%%rsi	/* c9 */\n\t"\
		"movaps	     (%%rcx),%%xmm4	\n\t"\
		"movaps	0x010(%%rdx),%%xmm0	\n\t"\
		"movaps	0x010(%%rcx),%%xmm5	\n\t"\
		"movaps	     (%%rdx),%%xmm6	\n\t"\
		"movaps	%%xmm4,%%xmm2		\n\t"\
		"movaps	%%xmm0,%%xmm1		\n\t"\
		"movaps	%%xmm5,%%xmm3		\n\t"\
		"movaps	%%xmm6,%%xmm7		\n\t"\
		"mulpd	    (%%rsi),%%xmm4	\n\t"\
		"mulpd	0x20(%%rsi),%%xmm0	\n\t"\
		"mulpd	    (%%rsi),%%xmm5	\n\t"\
		"mulpd	0x20(%%rsi),%%xmm6	\n\t"\
		"mulpd	0x10(%%rsi),%%xmm2	\n\t"\
		"mulpd	0x30(%%rsi),%%xmm1	\n\t"\
		"mulpd	0x10(%%rsi),%%xmm3	\n\t"\
		"mulpd	0x30(%%rsi),%%xmm7	\n\t"\
		"subpd	%%xmm2,%%xmm5		\n\t"\
		"subpd	%%xmm1,%%xmm6		\n\t"\
		"addpd	%%xmm3,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm0		\n\t"\
		"movaps	%%xmm5,0x10(%%rax)	\n\t"\
		"movaps	%%xmm6,0x10(%%rbx)	\n\t"\
		"movaps	%%xmm4,    (%%rax)	\n\t"\
		"movaps	%%xmm0,    (%%rbx)	\n\t"\
	"/*...Block 3: 31 MOVapd, 32 ADD/SUBpd, 20 MULpd */\n\t"\
		"movq	%[__add0],%%rax		\n\t"\
		"movslq	%[__p2],%%rsi		\n\t"\
		"movslq	%[__p4],%%rbx		\n\t"\
		"shlq	$3,%%rsi			\n\t"\
		"shlq	$3,%%rbx			\n\t"\
		"addq	%%rsi,%%rax			\n\t"\
		"addq	%%rax,%%rbx			\n\t"\
		"addq	$0x120,%%rcx	/* r5 +0x100 */\n\t"\
		"addq	$0x120,%%rdx	/* r13+0x100 */\n\t"\
		"movq	%[__isrt2],%%rsi	\n\t"\
		"movaps	(%%rsi),%%xmm2	/* isrt2 */\n\t"\
		"movaps	     (%%rcx),%%xmm4	\n\t"\
		"movaps	0x010(%%rcx),%%xmm5	\n\t"\
		"movaps	     (%%rdx),%%xmm6	\n\t"\
		"movaps	0x010(%%rdx),%%xmm7	\n\t"\
		"subq	$0x100,%%rcx		\n\t"\
		"subq	$0x100,%%rdx		\n\t"\
		"mulpd	%%xmm2,%%xmm4		\n\t"\
		"mulpd	%%xmm2,%%xmm5		\n\t"\
		"mulpd	%%xmm2,%%xmm6		\n\t"\
		"mulpd	%%xmm2,%%xmm7		\n\t"\
		"subpd	%%xmm4,%%xmm5		\n\t"\
		"movaps	     (%%rcx),%%xmm0	\n\t"\
		"subpd	%%xmm7,%%xmm6		\n\t"\
		"movaps	0x010(%%rdx),%%xmm2	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"movaps	0x010(%%rcx),%%xmm3	\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"movaps	     (%%rdx),%%xmm1	\n\t"\
		"addpd	%%xmm5,%%xmm4		\n\t"\
		"addpd	%%xmm6,%%xmm7		\n\t"\
		"subpd	%%xmm6,%%xmm4		\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm5		\n\t"\
		"subpd	%%xmm1,%%xmm3		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm1,%%xmm1		\n\t"\
		"addpd	%%xmm4,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm5,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm1		\n\t"\
		"/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c2) */\n\t"\
		"addq	$0x0b0,%%rsi	/* c2 */\n\t"\
		"subpd	%%xmm4,%%xmm2		\n\t"\
		"subpd	%%xmm7,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm3		\n\t"\
		"subpd	%%xmm6,%%xmm1		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm2,%%xmm4		\n\t"\
		"addpd	%%xmm0,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm5		\n\t"\
		"addpd	%%xmm1,%%xmm6		\n\t"\
		"movaps	%%xmm2,     (%%rcx)	\n\t"\
		"movaps	%%xmm0,0x010(%%rdx)	\n\t"\
		"movaps	%%xmm3,0x010(%%rcx)	\n\t"\
		"movaps	%%xmm6,     (%%rdx)	\n\t"\
		"movaps	%%xmm4,%%xmm2		\n\t"\
		"movaps	%%xmm7,%%xmm0		\n\t"\
		"movaps	%%xmm5,%%xmm3		\n\t"\
		"movaps	%%xmm1,%%xmm6		\n\t"\
		"mulpd	    (%%rsi),%%xmm4	\n\t"\
		"mulpd	0x20(%%rsi),%%xmm7	\n\t"\
		"mulpd	    (%%rsi),%%xmm5	\n\t"\
		"mulpd	0x20(%%rsi),%%xmm1	\n\t"\
		"mulpd	0x10(%%rsi),%%xmm2	\n\t"\
		"mulpd	0x30(%%rsi),%%xmm0	\n\t"\
		"mulpd	0x10(%%rsi),%%xmm3	\n\t"\
		"mulpd	0x30(%%rsi),%%xmm6	\n\t"\
		"subpd	%%xmm2,%%xmm5		\n\t"\
		"subpd	%%xmm0,%%xmm1		\n\t"\
		"addpd	%%xmm3,%%xmm4		\n\t"\
		"addpd	%%xmm6,%%xmm7		\n\t"\
		"movaps	%%xmm5,0x10(%%rax)	\n\t"\
		"movaps	%%xmm1,0x10(%%rbx)	\n\t"\
		"movaps	%%xmm4,    (%%rax)	\n\t"\
		"movaps	%%xmm7,    (%%rbx)	\n\t"\
		"addq	%%rdi,%%rax			\n\t"\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	$0x40,%%rsi	/* c10 */\n\t"\
		"movaps	     (%%rcx),%%xmm4	\n\t"\
		"movaps	0x010(%%rdx),%%xmm0	\n\t"\
		"movaps	0x010(%%rcx),%%xmm5	\n\t"\
		"movaps	     (%%rdx),%%xmm6	\n\t"\
		"movaps	%%xmm4,%%xmm2		\n\t"\
		"movaps	%%xmm0,%%xmm1		\n\t"\
		"movaps	%%xmm5,%%xmm3		\n\t"\
		"movaps	%%xmm6,%%xmm7		\n\t"\
		"mulpd	    (%%rsi),%%xmm4	\n\t"\
		"mulpd	0x20(%%rsi),%%xmm0	\n\t"\
		"mulpd	    (%%rsi),%%xmm5	\n\t"\
		"mulpd	0x20(%%rsi),%%xmm6	\n\t"\
		"mulpd	0x10(%%rsi),%%xmm2	\n\t"\
		"mulpd	0x30(%%rsi),%%xmm1	\n\t"\
		"mulpd	0x10(%%rsi),%%xmm3	\n\t"\
		"mulpd	0x30(%%rsi),%%xmm7	\n\t"\
		"subpd	%%xmm2,%%xmm5		\n\t"\
		"subpd	%%xmm1,%%xmm6		\n\t"\
		"addpd	%%xmm3,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm0		\n\t"\
		"movaps	%%xmm5,0x10(%%rax)	\n\t"\
		"movaps	%%xmm6,0x10(%%rbx)	\n\t"\
		"movaps	%%xmm4,    (%%rax)	\n\t"\
		"movaps	%%xmm0,    (%%rbx)	\n\t"\
	"/*...Block 4: 34 MOVapd, 36 ADD/SUBpd, 26 MULpd */\n\t"\
		"movq	%[__add0],%%rax		\n\t"\
		"movslq	%[__p3],%%rsi		\n\t"\
		"movslq	%[__p4],%%rbx		\n\t"\
		"shlq	$3,%%rsi			\n\t"\
		"shlq	$3,%%rbx			\n\t"\
		"addq	%%rsi,%%rax			\n\t"\
		"addq	%%rax,%%rbx			\n\t"\
		"addq	$0x120,%%rcx	/* r7 +0x100 */\n\t"\
		"addq	$0x120,%%rdx	/* r15+0x100 */\n\t"\
		"movq	%[__isrt2],%%rsi		\n\t"\
		"movaps	     (%%rcx),%%xmm4	\n\t"\
		"movaps	0x010(%%rcx),%%xmm5	\n\t"\
		"movaps	0x10(%%rsi),%%xmm3	/* cc0 */\n\t"\
		"movaps	0x20(%%rsi),%%xmm2	/* ss0 */\n\t"\
		"movaps	%%xmm4,%%xmm6		\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t"\
		"mulpd	%%xmm2,%%xmm4		\n\t"\
		"mulpd	%%xmm2,%%xmm5		\n\t"\
		"mulpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	     (%%rdx),%%xmm0	\n\t"\
		"mulpd	%%xmm3,%%xmm7		\n\t"\
		"movaps	0x010(%%rdx),%%xmm1	\n\t"\
		"subpd	%%xmm6,%%xmm5		\n\t"\
		"movaps	%%xmm0,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm4		\n\t"\
		"movaps	%%xmm1,%%xmm7		\n\t"\
		"mulpd	%%xmm3,%%xmm0		\n\t"\
		"mulpd	%%xmm3,%%xmm1		\n\t"\
		"mulpd	%%xmm2,%%xmm6		\n\t"\
		"mulpd	%%xmm2,%%xmm7		\n\t"\
		"subpd	%%xmm6,%%xmm1		\n\t"\
		"addpd	%%xmm7,%%xmm0		\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t"\
		"movaps	%%xmm4,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm6		\n\t"\
		"addpd	%%xmm1,%%xmm7		\n\t"\
		"subpd	%%xmm0,%%xmm4		\n\t"\
		"subpd	%%xmm1,%%xmm5		\n\t"\
		"subq	$0x100,%%rcx		\n\t"\
		"subq	$0x100,%%rdx		\n\t"\
		"movaps	     (%%rdx),%%xmm0	\n\t"\
		"movaps	0x010(%%rdx),%%xmm1	\n\t"\
		"movaps	(%%rsi),%%xmm3		\n\t"\
		"movaps	%%xmm0,%%xmm2		\n\t"\
		"subpd	%%xmm1,%%xmm0		\n\t"\
		"addpd	%%xmm2,%%xmm1		\n\t"\
		"mulpd	%%xmm3,%%xmm0		\n\t"\
		"mulpd	%%xmm3,%%xmm1		\n\t"\
		"movaps	     (%%rcx),%%xmm2	\n\t"\
		"movaps	0x010(%%rcx),%%xmm3	\n\t"\
		"subpd	%%xmm0,%%xmm2		\n\t"\
		"subpd	%%xmm1,%%xmm3		\n\t"\
		"addpd	%%xmm0,%%xmm0		\n\t"\
		"addpd	%%xmm1,%%xmm1		\n\t"\
		"addpd	%%xmm2,%%xmm0		\n\t"\
		"addpd	%%xmm3,%%xmm1		\n\t"\
		"/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c3) */\n\t"\
		"addq	$0x1b0,%%rsi	/* c3 */\n\t"\
		"subpd	%%xmm4,%%xmm2		\n\t"\
		"subpd	%%xmm7,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm3		\n\t"\
		"subpd	%%xmm6,%%xmm1		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm2,%%xmm4		\n\t"\
		"addpd	%%xmm0,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm5		\n\t"\
		"addpd	%%xmm1,%%xmm6		\n\t"\
		"movaps	%%xmm2,     (%%rcx)	\n\t"\
		"movaps	%%xmm0,0x010(%%rdx)	\n\t"\
		"movaps	%%xmm3,0x010(%%rcx)	\n\t"\
		"movaps	%%xmm6,     (%%rdx)	\n\t"\
		"movaps	%%xmm4,%%xmm2		\n\t"\
		"movaps	%%xmm7,%%xmm0		\n\t"\
		"movaps	%%xmm5,%%xmm3		\n\t"\
		"movaps	%%xmm1,%%xmm6		\n\t"\
		"mulpd	    (%%rsi),%%xmm4	\n\t"\
		"mulpd	0x20(%%rsi),%%xmm7	\n\t"\
		"mulpd	    (%%rsi),%%xmm5	\n\t"\
		"mulpd	0x20(%%rsi),%%xmm1	\n\t"\
		"mulpd	0x10(%%rsi),%%xmm2	\n\t"\
		"mulpd	0x30(%%rsi),%%xmm0	\n\t"\
		"mulpd	0x10(%%rsi),%%xmm3	\n\t"\
		"mulpd	0x30(%%rsi),%%xmm6	\n\t"\
		"subpd	%%xmm2,%%xmm5		\n\t"\
		"subpd	%%xmm0,%%xmm1		\n\t"\
		"addpd	%%xmm3,%%xmm4		\n\t"\
		"addpd	%%xmm6,%%xmm7		\n\t"\
		"movaps	%%xmm5,0x10(%%rax)	\n\t"\
		"movaps	%%xmm1,0x10(%%rbx)	\n\t"\
		"movaps	%%xmm4,    (%%rax)	\n\t"\
		"movaps	%%xmm7,    (%%rbx)	\n\t"\
		"addq	%%rdi,%%rax			\n\t"\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	$0x40,%%rsi	/* c11 */\n\t"\
		"movaps	     (%%rcx),%%xmm4	\n\t"\
		"movaps	0x010(%%rdx),%%xmm0	\n\t"\
		"movaps	0x010(%%rcx),%%xmm5	\n\t"\
		"movaps	     (%%rdx),%%xmm6	\n\t"\
		"movaps	%%xmm4,%%xmm2		\n\t"\
		"movaps	%%xmm0,%%xmm1		\n\t"\
		"movaps	%%xmm5,%%xmm3		\n\t"\
		"movaps	%%xmm6,%%xmm7		\n\t"\
		"mulpd	    (%%rsi),%%xmm4	\n\t"\
		"mulpd	0x20(%%rsi),%%xmm0	\n\t"\
		"mulpd	    (%%rsi),%%xmm5	\n\t"\
		"mulpd	0x20(%%rsi),%%xmm6	\n\t"\
		"mulpd	0x10(%%rsi),%%xmm2	\n\t"\
		"mulpd	0x30(%%rsi),%%xmm1	\n\t"\
		"mulpd	0x10(%%rsi),%%xmm3	\n\t"\
		"mulpd	0x30(%%rsi),%%xmm7	\n\t"\
		"subpd	%%xmm2,%%xmm5		\n\t"\
		"subpd	%%xmm1,%%xmm6		\n\t"\
		"addpd	%%xmm3,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm0		\n\t"\
		"movaps	%%xmm5,0x10(%%rax)	\n\t"\
		"movaps	%%xmm6,0x10(%%rbx)	\n\t"\
		"movaps	%%xmm4,    (%%rax)	\n\t"\
		"movaps	%%xmm0,    (%%rbx)	\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__p1] "m" (Xp1)\
		 ,[__p2] "m" (Xp2)\
		 ,[__p3] "m" (Xp3)\
		 ,[__p4] "m" (Xp4)\
		 ,[__p8] "m" (Xp8)\
		 ,[__r1] "m" (Xr1)\
		 ,[__isrt2] "m" (Xisrt2)\
		: "rax","rbx","rcx","rdx","rdi","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
	);\
	}

  #else // USE_64BIT_ASM_STYLE = True: Deeper 64-bit-ified version of the above 32-bit ASM macros, using all of xmm0-15

	#define SSE2_RADIX16_DIF_TWIDDLE(Xadd0,Xp1,Xp2,Xp3,Xp4,Xp8,Xp12,Xr1,Xisrt2)\
	{\
	__asm__ volatile (\
		"movq	%[__r1],%%r9						\n\t"\
		"movq	%[__add0],%%rax						\n\t"\
		"movslq	%[__p4],%%rbx						\n\t"\
		"movslq	%[__p8],%%rcx						\n\t	movq	%[__isrt2],%%r8			\n\t"\
		"movslq	%[__p12],%%rdx						\n\t	movslq	%[__p1],%%r10			\n\t"\
		"shlq	$3,%%rbx							\n\t	movslq	%[__p2],%%rdi			\n\t"\
		"shlq	$3,%%rcx							\n\t	shlq	$3,%%r10				\n\t"\
		"shlq	$3,%%rdx							\n\t	shlq	$3,%%rdi				\n\t"\
		"addq	%%rax,%%rbx							\n\t	addq	$0x230,%%r8	/* two */	\n\t"\
		"addq	%%rax,%%rcx							\n\t	movq	%%r10,%%r11				\n\t"\
		"addq	%%rax,%%rdx							\n\t	movq	%%r10,%%r12				\n\t"\
		"/* SSE2_RADIX4_DIF_4TWIDDLE_B(r1 ,c0) */	\n\t	movq	%%r10,%%r13				\n\t"\
		"/* Do	the p0,p8 combo: */					\n\t	addq	%%rax,%%r10				\n\t"\
		"leaq	0x230(%%r9),%%rsi	/* c0, from r1 */\n\t	addq	%%rbx,%%r11				\n\t"\
		"movaps	    (%%rax),%%xmm0					\n\t	addq	%%rcx,%%r12				\n\t"\
		"movaps	    (%%rcx),%%xmm4					\n\t	addq	%%rdx,%%r13				\n\t"\
		"movaps	0x10(%%rax),%%xmm1					\n\t	/* SSE2_RADIX4_DIF_4TWIDDLE_B(r17,c1) */\n\t"\
		"movaps	0x10(%%rcx),%%xmm5					\n\t	/* Do	the p0,p8 combo: */		\n\t"\
		"movaps	    (%%rsi),%%xmm6		/* c0 */	\n\t	movaps	    (%%r10),%%xmm8 		\n\t"\
		"movaps	0x10(%%rsi),%%xmm7					\n\t	movaps	    (%%r12),%%xmm12		\n\t"\
		"movaps	%%xmm0,%%xmm2						\n\t	movaps	0x10(%%r10),%%xmm9 		\n\t"\
		"movaps	%%xmm1,%%xmm3						\n\t	movaps	0x10(%%r12),%%xmm13		\n\t"\
		"mulpd	%%xmm6,%%xmm0						\n\t	movaps	0x100(%%rsi),%%xmm14	/* c1 */\n\t"\
		"mulpd	%%xmm6,%%xmm1						\n\t	movaps	0x110(%%rsi),%%xmm15	\n\t"\
		"mulpd	%%xmm7,%%xmm2						\n\t	movaps	%%xmm8 ,%%xmm10			\n\t"\
		"mulpd	%%xmm7,%%xmm3						\n\t	movaps	%%xmm9 ,%%xmm11			\n\t"\
		"movaps	%%xmm4,%%xmm6						\n\t	mulpd	%%xmm14,%%xmm8 			\n\t"\
		"addpd	%%xmm2,%%xmm1						\n\t	mulpd	%%xmm14,%%xmm9 			\n\t"\
		"movaps	%%xmm5,%%xmm7						\n\t	mulpd	%%xmm15,%%xmm10			\n\t"\
		"mulpd	0x20(%%rsi),%%xmm4		/* c8 */	\n\t	mulpd	%%xmm15,%%xmm11			\n\t"\
		"subpd	%%xmm3,%%xmm0						\n\t	movaps	%%xmm12,%%xmm14			\n\t"\
		"mulpd	0x20(%%rsi),%%xmm5					\n\t	addpd	%%xmm10,%%xmm9 			\n\t"\
		"mulpd	0x30(%%rsi),%%xmm6					\n\t	movaps	%%xmm13,%%xmm15			\n\t"\
		"movaps	%%xmm0,%%xmm2						\n\t	mulpd	0x120(%%rsi),%%xmm12	/* c9 */\n\t"\
		"mulpd	0x30(%%rsi),%%xmm7					\n\t	subpd	%%xmm11,%%xmm8 			\n\t"\
		"addpd	%%xmm6,%%xmm5						\n\t	mulpd	0x120(%%rsi),%%xmm13	\n\t"\
		"movaps	%%xmm1,%%xmm3						\n\t	mulpd	0x130(%%rsi),%%xmm14	\n\t"\
		"subpd	%%xmm7,%%xmm4						\n\t	movaps	%%xmm8 ,%%xmm10			\n\t"\
		"addpd	%%xmm4,%%xmm0						\n\t	mulpd	0x130(%%rsi),%%xmm15	\n\t"\
		"addpd	%%xmm5,%%xmm1						\n\t	addpd	%%xmm14,%%xmm13			\n\t"\
		"subpd	%%xmm4,%%xmm2						\n\t	movaps	%%xmm9 ,%%xmm11			\n\t"\
		"subpd	%%xmm5,%%xmm3						\n\t	subpd	%%xmm15,%%xmm12			\n\t"\
		"/* Do	the p4,12 combo: */					\n\t	addpd	%%xmm12,%%xmm8 			\n\t"\
		"movaps	    (%%rdx),%%xmm4					\n\t	addpd	%%xmm13,%%xmm9 			\n\t"\
		"movaps	0x10(%%rdx),%%xmm5					\n\t	subpd	%%xmm12,%%xmm10			\n\t"\
		"movaps	    (%%rdx),%%xmm6					\n\t	subpd	%%xmm13,%%xmm11			\n\t"\
		"movaps	0x10(%%rdx),%%xmm7					\n\t	/* Do	the p4,12 combo: */		\n\t"\
		"mulpd	0x60(%%rsi),%%xmm4		/* c12*/	\n\t	movaps	    (%%r13),%%xmm12		\n\t"\
		"mulpd	0x60(%%rsi),%%xmm5					\n\t	movaps	0x10(%%r13),%%xmm13		\n\t"\
		"mulpd	0x70(%%rsi),%%xmm6					\n\t	movaps	    (%%r13),%%xmm14		\n\t"\
		"mulpd	0x70(%%rsi),%%xmm7					\n\t	movaps	0x10(%%r13),%%xmm15		\n\t"\
		"addpd	%%xmm6,%%xmm5						\n\t	mulpd	0x160(%%rsi),%%xmm12	/* c13*/\n\t"\
		"subpd	%%xmm7,%%xmm4						\n\t	mulpd	0x160(%%rsi),%%xmm13	\n\t"\
		"movaps	%%xmm5,0x010(%%r9)					\n\t	mulpd	0x170(%%rsi),%%xmm14	\n\t"\
		"movaps	%%xmm4,     (%%r9)					\n\t	mulpd	0x170(%%rsi),%%xmm15	\n\t"\
		"movaps	    (%%rbx),%%xmm4					\n\t	addpd	%%xmm14,%%xmm13			\n\t"\
		"movaps	0x10(%%rbx),%%xmm5					\n\t	subpd	%%xmm15,%%xmm12			\n\t"\
		"movaps	    (%%rbx),%%xmm6					\n\t	movaps	%%xmm13,0x110(%%r9)		/* r17 */\n\t"\
		"movaps	0x10(%%rbx),%%xmm7					\n\t	movaps	%%xmm12,0x100(%%r9)		\n\t"\
		"mulpd	0x40(%%rsi),%%xmm4		/* c4 */	\n\t	movaps	    (%%r11),%%xmm12		\n\t"\
		"mulpd	0x40(%%rsi),%%xmm5					\n\t	movaps	0x10(%%r11),%%xmm13		\n\t"\
		"mulpd	0x50(%%rsi),%%xmm6					\n\t	movaps	    (%%r11),%%xmm14		\n\t"\
		"mulpd	0x50(%%rsi),%%xmm7					\n\t	movaps	0x10(%%r11),%%xmm15		\n\t"\
		"addpd	%%xmm6,%%xmm5						\n\t	mulpd	0x140(%%rsi),%%xmm12	/* c5 */\n\t"\
		"subpd	%%xmm7,%%xmm4						\n\t	mulpd	0x140(%%rsi),%%xmm13	\n\t"\
		"movaps	%%xmm5,%%xmm7						\n\t	mulpd	0x150(%%rsi),%%xmm14	\n\t"\
		"movaps	%%xmm4,%%xmm6						\n\t	mulpd	0x150(%%rsi),%%xmm15	\n\t"\
		"subpd	     (%%r9),%%xmm4					\n\t	addpd	%%xmm14,%%xmm13			\n\t"\
		"subpd	0x010(%%r9),%%xmm5					\n\t	subpd	%%xmm15,%%xmm12			\n\t"\
		"addpd	     (%%r9),%%xmm6					\n\t	movaps	%%xmm13,%%xmm15			\n\t"\
		"addpd	0x010(%%r9),%%xmm7					\n\t	movaps	%%xmm12,%%xmm14			\n\t"\
		"/* Finish radix-4 butterfly, store: */		\n\t	subpd	0x100(%%r9),%%xmm12		\n\t"\
		"subpd	%%xmm6,%%xmm0						\n\t	subpd	0x110(%%r9),%%xmm13		\n\t"\
		"subpd	%%xmm5,%%xmm2						\n\t	addpd	0x100(%%r9),%%xmm14		\n\t"\
		"subpd	%%xmm7,%%xmm1						\n\t	addpd	0x110(%%r9),%%xmm15		\n\t"\
		"subpd	%%xmm4,%%xmm3						\n\t	/* Finish radix-4 butterfly and store results into temporary-array slots: */\n\t"\
		"movaps	%%xmm0,0x040(%%r9)					\n\t	subpd	%%xmm14,%%xmm8 			\n\t"\
		"movaps	%%xmm2,0x020(%%r9)					\n\t	subpd	%%xmm13,%%xmm10			\n\t"\
		"movaps	%%xmm1,0x050(%%r9)					\n\t	subpd	%%xmm15,%%xmm9 			\n\t"\
		"movaps	%%xmm3,0x070(%%r9)					\n\t	subpd	%%xmm12,%%xmm11			\n\t"\
		"addpd	%%xmm6,%%xmm6						\n\t	movaps	%%xmm8 ,0x140(%%r9)		\n\t"\
		"addpd	%%xmm5,%%xmm5						\n\t	movaps	%%xmm10,0x120(%%r9)		\n\t"\
		"addpd	%%xmm7,%%xmm7						\n\t	movaps	%%xmm9 ,0x150(%%r9)		\n\t"\
		"addpd	%%xmm4,%%xmm4						\n\t	movaps	%%xmm11,0x170(%%r9)		\n\t"\
		"addpd	%%xmm0,%%xmm6						\n\t	mulpd	(%%r8),%%xmm14			\n\t"\
		"addpd	%%xmm2,%%xmm5						\n\t	mulpd	(%%r8),%%xmm13			\n\t"\
		"addpd	%%xmm1,%%xmm7						\n\t	mulpd	(%%r8),%%xmm15			\n\t"\
		"addpd	%%xmm3,%%xmm4						\n\t	mulpd	(%%r8),%%xmm12			\n\t"\
		"movaps	%%xmm6,     (%%r9)					\n\t	addpd	%%xmm8 ,%%xmm14			\n\t"\
		"movaps	%%xmm5,0x060(%%r9)					\n\t	addpd	%%xmm10,%%xmm13			\n\t"\
		"movaps	%%xmm7,0x010(%%r9)					\n\t	addpd	%%xmm9 ,%%xmm15			\n\t"\
		"movaps	%%xmm4,0x030(%%r9)					\n\t	addpd	%%xmm11,%%xmm12			\n\t"\
		"addq	%%rdi,%%rax							\n\t	movaps	%%xmm14,0x100(%%r9)		\n\t"\
		"addq	%%rdi,%%rbx							\n\t	movaps	%%xmm13,0x160(%%r9)		\n\t"\
		"addq	%%rdi,%%rcx							\n\t	movaps	%%xmm15,0x110(%%r9)		\n\t"\
		"addq	%%rdi,%%rdx							\n\t	movaps	%%xmm12,0x130(%%r9)		\n\t"\
		"/* SSE2_RADIX4_DIF_4TWIDDLE_B(r9 ,c2) */	\n\t	addq	%%rdi,%%r10				\n\t"\
		"/* Do	the p0,p8 combo: */					\n\t	addq	%%rdi,%%r11				\n\t"\
		"addq	$0x80,%%rsi 		/* c2 */		\n\t	addq	%%rdi,%%r12				\n\t"\
		"movaps	    (%%rax),%%xmm0					\n\t	addq	%%rdi,%%r13				\n\t"\
		"movaps	    (%%rcx),%%xmm4					\n\t	/* SSE2_RADIX4_DIF_4TWIDDLE_B(r25,c3) */\n\t"\
		"movaps	0x10(%%rax),%%xmm1					\n\t	/* Do	the p0,p8 combo: */		\n\t"\
		"movaps	0x10(%%rcx),%%xmm5					\n\t	movaps	    (%%r10),%%xmm8 		\n\t"\
		"movaps	    (%%rsi),%%xmm6					\n\t	movaps	    (%%r12),%%xmm12		\n\t"\
		"movaps	0x10(%%rsi),%%xmm7					\n\t	movaps	0x10(%%r10),%%xmm9 		\n\t"\
		"movaps	%%xmm0,%%xmm2						\n\t	movaps	0x10(%%r12),%%xmm13		\n\t"\
		"movaps	%%xmm1,%%xmm3						\n\t	movaps	0x100(%%rsi),%%xmm14	/* c3 */\n\t"\
		"mulpd	%%xmm6,%%xmm0						\n\t	movaps	0x110(%%rsi),%%xmm15	\n\t"\
		"mulpd	%%xmm6,%%xmm1						\n\t	movaps	%%xmm8 ,%%xmm10			\n\t"\
		"mulpd	%%xmm7,%%xmm2						\n\t	movaps	%%xmm9 ,%%xmm11			\n\t"\
		"mulpd	%%xmm7,%%xmm3						\n\t	mulpd	 %%xmm14,%%xmm8 		\n\t"\
		"movaps	%%xmm4,%%xmm6						\n\t	mulpd	 %%xmm14,%%xmm9 		\n\t"\
		"addpd	%%xmm2,%%xmm1						\n\t	mulpd	 %%xmm15,%%xmm10		\n\t"\
		"movaps	%%xmm5,%%xmm7						\n\t	mulpd	 %%xmm15,%%xmm11		\n\t"\
		"mulpd	0x20(%%rsi),%%xmm4	/* c10*/		\n\t	movaps	%%xmm12,%%xmm14			\n\t"\
		"subpd	%%xmm3,%%xmm0						\n\t	addpd	 %%xmm10,%%xmm9 		\n\t"\
		"mulpd	0x20(%%rsi),%%xmm5					\n\t	movaps	%%xmm13,%%xmm15			\n\t"\
		"mulpd	0x30(%%rsi),%%xmm6					\n\t	mulpd	0x120(%%rsi),%%xmm12	/* c11*/\n\t"\
		"movaps	%%xmm0,%%xmm2						\n\t	subpd	 %%xmm11,%%xmm8 		\n\t"\
		"mulpd	0x30(%%rsi),%%xmm7					\n\t	mulpd	0x120(%%rsi),%%xmm13	\n\t"\
		"addpd	%%xmm6,%%xmm5						\n\t	mulpd	0x130(%%rsi),%%xmm14	\n\t"\
		"movaps	%%xmm1,%%xmm3						\n\t	movaps	%%xmm8 ,%%xmm10			\n\t"\
		"subpd	%%xmm7,%%xmm4						\n\t	mulpd	0x130(%%rsi),%%xmm15	\n\t"\
		"addpd	%%xmm4,%%xmm0						\n\t	addpd	 %%xmm14,%%xmm13		\n\t"\
		"addpd	%%xmm5,%%xmm1						\n\t	movaps	%%xmm9 ,%%xmm11			\n\t"\
		"subpd	%%xmm4,%%xmm2						\n\t	subpd	%%xmm15,%%xmm12			\n\t"\
		"subpd	%%xmm5,%%xmm3						\n\t	addpd	%%xmm12,%%xmm8 			\n\t"\
		"/* Do	the p4,12 combo: */					\n\t	addpd	%%xmm13,%%xmm9 			\n\t"\
		"addq	$0x80,%%r9 			/* r9 */		\n\t	subpd	%%xmm12,%%xmm10			\n\t"\
		"movaps	    (%%rdx),%%xmm4					\n\t	subpd	%%xmm13,%%xmm11			\n\t"\
		"movaps	0x10(%%rdx),%%xmm5					\n\t	/* Do	the p4,12 combo: */		\n\t"\
		"movaps	    (%%rdx),%%xmm6					\n\t	movaps	    (%%r13),%%xmm12		\n\t"\
		"movaps	0x10(%%rdx),%%xmm7					\n\t	movaps	0x10(%%r13),%%xmm13		\n\t"\
		"mulpd	0x60(%%rsi),%%xmm4	/* c14*/		\n\t	movaps	    (%%r13),%%xmm14		\n\t"\
		"mulpd	0x60(%%rsi),%%xmm5					\n\t	movaps	0x10(%%r13),%%xmm15		\n\t"\
		"mulpd	0x70(%%rsi),%%xmm6					\n\t	mulpd	0x160(%%rsi),%%xmm12	/* c15*/\n\t"\
		"mulpd	0x70(%%rsi),%%xmm7					\n\t	mulpd	0x160(%%rsi),%%xmm13	\n\t"\
		"addpd	%%xmm6,%%xmm5						\n\t	mulpd	0x170(%%rsi),%%xmm14	\n\t"\
		"subpd	%%xmm7,%%xmm4						\n\t	mulpd	0x170(%%rsi),%%xmm15	\n\t"\
		"movaps	%%xmm5,0x010(%%r9)					\n\t	addpd	%%xmm14,%%xmm13			\n\t"\
		"movaps	%%xmm4,     (%%r9)					\n\t	subpd	%%xmm15,%%xmm12			\n\t"\
		"movaps	    (%%rbx),%%xmm4					\n\t	movaps	%%xmm13,0x110(%%r9)		/* r25*/\n\t"\
		"movaps	0x10(%%rbx),%%xmm5					\n\t	movaps	%%xmm12,0x100(%%r9)		\n\t"\
		"movaps	    (%%rbx),%%xmm6					\n\t	movaps	    (%%r11),%%xmm12		\n\t"\
		"movaps	0x10(%%rbx),%%xmm7					\n\t	movaps	0x10(%%r11),%%xmm13		\n\t"\
		"mulpd	0x40(%%rsi),%%xmm4	/* c6 */		\n\t	movaps	    (%%r11),%%xmm14		\n\t"\
		"mulpd	0x40(%%rsi),%%xmm5					\n\t	movaps	0x10(%%r11),%%xmm15		\n\t"\
		"mulpd	0x50(%%rsi),%%xmm6					\n\t	mulpd	0x140(%%rsi),%%xmm12	/* c7 */\n\t"\
		"mulpd	0x50(%%rsi),%%xmm7					\n\t	mulpd	0x140(%%rsi),%%xmm13	\n\t"\
		"addpd	%%xmm6,%%xmm5						\n\t	mulpd	0x150(%%rsi),%%xmm14	\n\t"\
		"subpd	%%xmm7,%%xmm4						\n\t	mulpd	0x150(%%rsi),%%xmm15	\n\t"\
		"movaps	%%xmm5,%%xmm7						\n\t	addpd	%%xmm14,%%xmm13			\n\t"\
		"movaps	%%xmm4,%%xmm6						\n\t	subpd	%%xmm15,%%xmm12			\n\t"\
		"subpd	     (%%r9),%%xmm4					\n\t	movaps	%%xmm13,%%xmm15			\n\t"\
		"subpd	0x010(%%r9),%%xmm5					\n\t	movaps	%%xmm12,%%xmm14			\n\t"\
		"addpd	     (%%r9),%%xmm6					\n\t	subpd	0x100(%%r9),%%xmm12		\n\t"\
		"addpd	0x010(%%r9),%%xmm7					\n\t	subpd	0x110(%%r9),%%xmm13		\n\t"\
		"/* Finish radix-4 butterfly, store: */		\n\t	addpd	0x100(%%r9),%%xmm14		\n\t"\
		"subpd	%%xmm6,%%xmm0						\n\t	addpd	0x110(%%r9),%%xmm15		\n\t"\
		"subpd	%%xmm5,%%xmm2						\n\t	/* Finish radix-4 butterfly and store results into temporary-array slots: */\n\t"\
		"subpd	%%xmm7,%%xmm1						\n\t	subpd	%%xmm14,%%xmm8 			\n\t"\
		"subpd	%%xmm4,%%xmm3						\n\t	subpd	%%xmm13,%%xmm10			\n\t"\
		"movaps	%%xmm0,0x040(%%r9)					\n\t	subpd	%%xmm15,%%xmm9 			\n\t"\
		"movaps	%%xmm2,0x020(%%r9)					\n\t	subpd	%%xmm12,%%xmm11			\n\t"\
		"movaps	%%xmm1,0x050(%%r9)					\n\t	movaps	%%xmm8 ,0x140(%%r9)		\n\t"\
		"movaps	%%xmm3,0x070(%%r9)					\n\t	movaps	%%xmm10,0x120(%%r9)		\n\t"\
		"addpd	%%xmm6,%%xmm6						\n\t	movaps	%%xmm9 ,0x150(%%r9)		\n\t"\
		"addpd	%%xmm5,%%xmm5						\n\t	movaps	%%xmm11,0x170(%%r9)		\n\t"\
		"addpd	%%xmm7,%%xmm7						\n\t	mulpd	(%%r8),%%xmm14			\n\t"\
		"addpd	%%xmm4,%%xmm4						\n\t	mulpd	(%%r8),%%xmm13			\n\t"\
		"addpd	%%xmm0,%%xmm6						\n\t	mulpd	(%%r8),%%xmm15			\n\t"\
		"addpd	%%xmm2,%%xmm5						\n\t	mulpd	(%%r8),%%xmm12			\n\t"\
		"addpd	%%xmm1,%%xmm7						\n\t	addpd	%%xmm8 ,%%xmm14			\n\t"\
		"addpd	%%xmm3,%%xmm4						\n\t	addpd	%%xmm10,%%xmm13			\n\t"\
		"movaps	%%xmm6,     (%%r9)					\n\t	addpd	%%xmm9 ,%%xmm15			\n\t"\
		"movaps	%%xmm5,0x060(%%r9)					\n\t	addpd	%%xmm11,%%xmm12			\n\t"\
		"movaps	%%xmm7,0x010(%%r9)					\n\t	movaps	%%xmm14,0x100(%%r9)		\n\t"\
		"movaps	%%xmm4,0x030(%%r9)					\n\t	movaps	%%xmm13,0x160(%%r9)		\n\t"\
		"											\n\t	movaps	%%xmm15,0x110(%%r9)		\n\t"\
		"											\n\t	movaps	%%xmm12,0x130(%%r9)		\n\t"\
		"/*************************************************************************************/\n\t"\
		"/*  And now do four more radix-4 transforms, including the internal twiddle factors: */\n\t"\
		"/*************************************************************************************/\n\t"\
		"/* Block 1 - p2 in rdi, add0+p2,3 already in rax,r10: */	\n\t"\
		"movq	%%r10,%%rbx		/* cpy add0+p3*/	\n\t"\
		"movq	%%rax,%%rcx		/* add0+p2 */		\n\t"\
		"movq	%%r10,%%rdx		/* add0+p3*/		/* Block 3: */				\n\t"\
		"subq	%%rdi,%%rax		/* add0    */		\n\t	movslq	%[__p4],%%r10			\n\t"\
		"subq	%%rdi,%%rbx		/* add0+p1 */		\n\t	movslq	%[__p8],%%rdi			\n\t"\
		"movq	%[__isrt2],%%rsi					\n\t	shlq	$3,%%r10				\n\t"\
		"subq	$0x80,%%r9		/* r1 */			\n\t	shlq	$3,%%rdi				\n\t"\
		"movaps	     (%%r9),%%xmm0					\n\t	movq	%%r10,%%r11				\n\t"\
		"movaps	0x100(%%r9),%%xmm4					\n\t	movq	%%r10,%%r12				\n\t"\
		"movaps	0x010(%%r9),%%xmm1					\n\t	movq	%%r10,%%r13				\n\t"\
		"movaps	0x110(%%r9),%%xmm5					\n\t	addq	%%rax,%%r10				\n\t"\
		"movaps	0x080(%%r9),%%xmm2					\n\t	addq	%%rbx,%%r11				\n\t"\
		"movaps	0x180(%%r9),%%xmm6					\n\t	addq	%%rcx,%%r12				\n\t"\
		"movaps	0x090(%%r9),%%xmm3					\n\t	addq	%%rdx,%%r13				\n\t"\
		"movaps	0x190(%%r9),%%xmm7					\n\t	movaps	(%%rsi),%%xmm11	/* isrt2 */	\n\t"\
		"subpd	%%xmm2,%%xmm0						\n\t	movaps	0x140(%%r9),%%xmm12	/* r5 */\n\t"\
		"subpd	%%xmm6,%%xmm4						\n\t	movaps	0x150(%%r9),%%xmm13	\n\t"\
		"subpd	%%xmm3,%%xmm1						\n\t	movaps	0x1c0(%%r9),%%xmm14	\n\t"\
		"subpd	%%xmm7,%%xmm5						\n\t	movaps	0x1d0(%%r9),%%xmm15	\n\t"\
		"addpd	%%xmm2,%%xmm2						\n\t	mulpd	%%xmm11,%%xmm12		\n\t"\
		"addpd	%%xmm6,%%xmm6						\n\t	movaps	0x040(%%r9),%%xmm8 	\n\t"\
		"addpd	%%xmm3,%%xmm3						\n\t	mulpd	%%xmm11,%%xmm13		\n\t"\
		"addpd	%%xmm7,%%xmm7						\n\t	movaps	0x050(%%r9),%%xmm9 	\n\t"\
		"addpd	%%xmm0,%%xmm2						\n\t	mulpd	%%xmm11,%%xmm14		\n\t"\
		"addpd	%%xmm4,%%xmm6						\n\t	movaps	0x0c0(%%r9),%%xmm10	\n\t"\
		"addpd	%%xmm1,%%xmm3						\n\t	mulpd	%%xmm11,%%xmm15		\n\t"\
		"addpd	%%xmm5,%%xmm7						\n\t	movaps	0x0d0(%%r9),%%xmm11	\n\t"\
		"subpd	%%xmm6,%%xmm2						\n\t	subpd	%%xmm11,%%xmm8 		\n\t"\
		"subpd	%%xmm7,%%xmm3						\n\t	subpd	%%xmm13,%%xmm12		\n\t"\
		"mulpd	(%%r8),%%xmm6						\n\t	subpd	%%xmm10,%%xmm9 		\n\t"\
		"mulpd	(%%r8),%%xmm7						\n\t	subpd	%%xmm14,%%xmm15		\n\t"\
		"movaps	%%xmm2,    (%%rbx)					\n\t	addpd	%%xmm11,%%xmm11		\n\t"\
		"movaps	%%xmm3,0x10(%%rbx)					\n\t	addpd	%%xmm13,%%xmm13		\n\t"\
		"addpd	%%xmm2,	%%xmm6						\n\t	mulpd	(%%r8),%%xmm10		\n\t"\
		"addpd	%%xmm3,	%%xmm7						\n\t	mulpd	(%%r8),%%xmm14		\n\t"\
		"movaps	%%xmm6,    (%%rax)					\n\t	addpd	%%xmm8 ,%%xmm11		\n\t"\
		"movaps	%%xmm7,0x10(%%rax)					\n\t	addpd	%%xmm12,%%xmm13		\n\t"\
		"subpd	%%xmm5,	%%xmm0						\n\t	addpd	%%xmm9 ,%%xmm10		\n\t"\
		"subpd	%%xmm4,	%%xmm1						\n\t	addpd	%%xmm15,%%xmm14		\n\t"\
		"mulpd	(%%r8),	%%xmm5						\n\t	subpd	%%xmm14,%%xmm12		\n\t"\
		"mulpd	(%%r8),	%%xmm4						\n\t	subpd	%%xmm15,%%xmm13		\n\t"\
		"movaps	%%xmm0,    (%%rcx)					\n\t	addpd	%%xmm14,%%xmm14		\n\t"\
		"movaps	%%xmm1,0x10(%%rdx)					\n\t	addpd	%%xmm15,%%xmm15		\n\t"\
		"addpd	%%xmm0,	%%xmm5						\n\t	addpd	%%xmm12,%%xmm14		\n\t"\
		"addpd	%%xmm1,	%%xmm4						\n\t	addpd	%%xmm13,%%xmm15		\n\t"\
		"movaps	%%xmm5,    (%%rdx)					\n\t	subpd	%%xmm12,%%xmm8 		\n\t"\
		"movaps	%%xmm4,0x10(%%rcx)					\n\t	subpd	%%xmm15,%%xmm11		\n\t"\
		"											\n\t	subpd	%%xmm13,%%xmm10		\n\t"\
		"											\n\t	subpd	%%xmm14,%%xmm9 		\n\t"\
		"											\n\t	mulpd	(%%r8),%%xmm12		\n\t"\
		"											\n\t	mulpd	(%%r8),%%xmm15		\n\t"\
		"											\n\t	mulpd	(%%r8),%%xmm13		\n\t"\
		"											\n\t	mulpd	(%%r8),%%xmm14		\n\t"\
		"											\n\t	movaps	%%xmm8 ,    (%%r11)	\n\t"\
		"											\n\t	movaps	%%xmm11,    (%%r12)	\n\t"\
		"											\n\t	movaps	%%xmm10,0x10(%%r11)	\n\t"\
		"											\n\t	movaps	%%xmm9 ,0x10(%%r13)	\n\t"\
		"/* Block 2: */								\n\t	addpd	%%xmm8 ,	%%xmm12	\n\t"\
		"addq	%%rdi,%%rax							\n\t	addpd	%%xmm11,	%%xmm15	\n\t"\
		"addq	%%rdi,%%rbx							\n\t	addpd	%%xmm10,	%%xmm13	\n\t"\
		"addq	%%rdi,%%rcx							\n\t	addpd	%%xmm9 ,	%%xmm14	\n\t"\
		"addq	%%rdi,%%rdx							\n\t	movaps	%%xmm12,    (%%r10)	\n\t"\
		"addq	$0x20,%%r9	/* r3 */				\n\t	movaps	%%xmm15,    (%%r13)	\n\t"\
		"movaps	0x100(%%r9),%%xmm4					\n\t	movaps	%%xmm13,0x10(%%r10)	\n\t"\
		"movaps	0x110(%%r9),%%xmm5					\n\t	movaps	%%xmm14,0x10(%%r12)	\n\t"\
		"/* Share cc0/ss0 between 2 halves: */		\n\t	/* Block 4: */				\n\t"\
		"movaps	0x10(%%rsi),%%xmm11	/* cc0 */		\n\t	addq	%%rdi,%%r10			\n\t"\
		"movaps	0x20(%%rsi),%%xmm10	/* ss0 */		\n\t	addq	%%rdi,%%r11			\n\t"\
		"movaps	%%xmm4,%%xmm6						\n\t	addq	%%rdi,%%r12			\n\t"\
		"movaps	%%xmm5,%%xmm7						\n\t	addq	%%rdi,%%r13			\n\t"\
		"mulpd	%%xmm11,%%xmm4						\n\t	movaps	0x140(%%r9),%%xmm12	/* r7 */\n\t"\
		"mulpd	%%xmm11,%%xmm5						\n\t	movaps	0x150(%%r9),%%xmm13	\n\t"\
		"mulpd	%%xmm10,%%xmm6						\n\t	movaps	%%xmm12,%%xmm14		\n\t"\
		"movaps	0x180(%%r9),%%xmm0					\n\t	movaps	%%xmm13,%%xmm15		\n\t"\
		"mulpd	%%xmm10,%%xmm7						\n\t	mulpd	%%xmm10,%%xmm12		\n\t"\
		"movaps	0x190(%%r9),%%xmm1					\n\t	mulpd	%%xmm10,%%xmm13		\n\t"\
		"addpd	%%xmm6,%%xmm5						\n\t	mulpd	%%xmm11,%%xmm14		\n\t"\
		"movaps	%%xmm0,%%xmm6						\n\t	movaps	0x1c0(%%r9),%%xmm8 	\n\t"\
		"subpd	%%xmm7,%%xmm4						\n\t	mulpd	%%xmm11,%%xmm15		\n\t"\
		"movaps	%%xmm1,%%xmm7						\n\t	movaps	0x1d0(%%r9),%%xmm9 	\n\t"\
		"mulpd	%%xmm10,%%xmm6						\n\t	addpd	%%xmm14,%%xmm13		\n\t"\
		"mulpd	%%xmm10,%%xmm7						\n\t	movaps	%%xmm8 ,%%xmm14		\n\t"\
		"mulpd	%%xmm11,%%xmm0						\n\t	subpd	%%xmm15,%%xmm12		\n\t"\
		"mulpd	%%xmm11,%%xmm1						\n\t	movaps	%%xmm9 ,%%xmm15		\n\t"\
		"addpd	%%xmm0,%%xmm7						\n\t	mulpd	%%xmm11,%%xmm14		\n\t"\
		"subpd	%%xmm1,%%xmm6						\n\t	mulpd	%%xmm11,%%xmm15		\n\t"\
		"movaps	%%xmm4,%%xmm2						\n\t	mulpd	%%xmm10,%%xmm8 		\n\t"\
		"movaps	%%xmm5,%%xmm3						\n\t	mulpd	%%xmm10,%%xmm9 		\n\t"\
		"subpd	%%xmm6,%%xmm4						\n\t	addpd	%%xmm8 ,%%xmm15		\n\t"\
		"subpd	%%xmm7,%%xmm5						\n\t	subpd	%%xmm9 ,%%xmm14		\n\t"\
		"addpd	%%xmm2,%%xmm6						\n\t	movaps	%%xmm12,%%xmm10		\n\t"\
		"addpd	%%xmm3,%%xmm7						\n\t	movaps	%%xmm13,%%xmm11		\n\t"\
		"movaps	0x080(%%r9),%%xmm2					\n\t	subpd	%%xmm14,%%xmm12		\n\t"\
		"movaps	0x090(%%r9),%%xmm3					\n\t	subpd	%%xmm15,%%xmm13		\n\t"\
		"movaps	(%%rsi),%%xmm1	/* isrt2 */			\n\t	addpd	%%xmm10,%%xmm14		\n\t"\
		"movaps	%%xmm2,%%xmm0						\n\t	addpd	%%xmm11,%%xmm15		\n\t"\
		"subpd	%%xmm3,%%xmm2						\n\t	movaps	0x0c0(%%r9),%%xmm10	\n\t"\
		"addpd	%%xmm0,%%xmm3						\n\t	movaps	0x0d0(%%r9),%%xmm11	\n\t"\
		"mulpd	%%xmm1,%%xmm2						\n\t	movaps	(%%rsi),%%xmm9 	/* isrt2 */\n\t"\
		"mulpd	%%xmm1,%%xmm3						\n\t	movaps	%%xmm10,%%xmm8 		\n\t"\
		"movaps	     (%%r9),%%xmm0					\n\t	addpd	%%xmm11,%%xmm10		\n\t"\
		"movaps	0x010(%%r9),%%xmm1					\n\t	subpd	%%xmm8 ,%%xmm11		\n\t"\
		"subpd	%%xmm2,%%xmm0						\n\t	mulpd	%%xmm9 ,%%xmm10		\n\t"\
		"subpd	%%xmm3,%%xmm1						\n\t	mulpd	%%xmm9 ,%%xmm11		\n\t"\
		"addpd	%%xmm2,%%xmm2						\n\t	movaps	0x040(%%r9),%%xmm8 	\n\t"\
		"addpd	%%xmm3,%%xmm3						\n\t	movaps	0x050(%%r9),%%xmm9 	\n\t"\
		"addpd	%%xmm0,%%xmm2						\n\t	subpd	%%xmm10,%%xmm8 		\n\t"\
		"addpd	%%xmm1,%%xmm3						\n\t	subpd	%%xmm11,%%xmm9 		\n\t"\
		"subpd	%%xmm6,%%xmm2						\n\t	mulpd	(%%r8),%%xmm10		\n\t"\
		"subpd	%%xmm5,%%xmm0						\n\t	mulpd	(%%r8),%%xmm11		\n\t"\
		"subpd	%%xmm7,%%xmm3						\n\t	addpd	%%xmm8 ,%%xmm10		\n\t"\
		"subpd	%%xmm4,%%xmm1						\n\t	addpd	%%xmm9 ,%%xmm11		\n\t"\
		"mulpd	(%%r8),%%xmm6						\n\t	subpd	%%xmm12,%%xmm8 		\n\t"\
		"mulpd	(%%r8),%%xmm5						\n\t	subpd	%%xmm15,%%xmm10		\n\t"\
		"mulpd	(%%r8),%%xmm7						\n\t	subpd	%%xmm13,%%xmm9 		\n\t"\
		"mulpd	(%%r8),%%xmm4						\n\t	subpd	%%xmm14,%%xmm11		\n\t"\
		"movaps	%%xmm2,    (%%rbx)					\n\t	addpd	%%xmm12,%%xmm12		\n\t"\
		"movaps	%%xmm0,    (%%rcx)					\n\t	addpd	%%xmm15,%%xmm15		\n\t"\
		"movaps	%%xmm3,0x10(%%rbx)					\n\t	addpd	%%xmm13,%%xmm13		\n\t"\
		"movaps	%%xmm1,0x10(%%rdx)					\n\t	addpd	%%xmm14,%%xmm14		\n\t"\
		"addpd	%%xmm2,	%%xmm6						\n\t	movaps	%%xmm8 ,    (%%r11)	\n\t"\
		"addpd	%%xmm0,	%%xmm5						\n\t	movaps	%%xmm10,    (%%r12)	\n\t"\
		"addpd	%%xmm3,	%%xmm7						\n\t	movaps	%%xmm9 ,0x10(%%r11)	\n\t"\
		"addpd	%%xmm1,	%%xmm4						\n\t	movaps	%%xmm11,0x10(%%r13)	\n\t"\
		"movaps	%%xmm6,    (%%rax)					\n\t	addpd	%%xmm8 ,%%xmm12		\n\t"\
		"movaps	%%xmm5,    (%%rdx)					\n\t	addpd	%%xmm10,%%xmm15		\n\t"\
		"movaps	%%xmm7,0x10(%%rax)					\n\t	addpd	%%xmm9 ,%%xmm13		\n\t"\
		"movaps	%%xmm4,0x10(%%rcx)					\n\t	addpd	%%xmm11,%%xmm14		\n\t"\
		"													movaps	%%xmm12,    (%%r10)	\n\t"\
		"													movaps	%%xmm15,    (%%r13)	\n\t"\
		"													movaps	%%xmm13,0x10(%%r10)	\n\t"\
		"													movaps	%%xmm14,0x10(%%r12)	\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__p1] "m" (Xp1)\
		 ,[__p2] "m" (Xp2)\
		 ,[__p3] "m" (Xp3)\
		 ,[__p4] "m" (Xp4)\
		 ,[__p8] "m" (Xp8)\
		 ,[__p12] "m" (Xp12)\
		 ,[__r1] "m" (Xr1)\
		 ,[__isrt2] "m" (Xisrt2)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r8","r9","r10","r11","r12","r13","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
	);\
	}

	#define SSE2_RADIX16_DIT_TWIDDLE(Xadd0,Xp1,Xp2,Xp3,Xp4,Xp8,Xr1,Xisrt2)\
	{\
	__asm__ volatile (\
		"movq	%[__isrt2],%%r8		\n\t"\
		"movq	%[__add0],%%rax		\n\t"\
		"movslq	%[__p1],%%rbx		\n\t"\
		"movslq	%[__p2],%%rcx		\n\t"\
		"movslq	%[__p3],%%rdx		\n\t"\
		"movslq	%[__p4],%%rdi		\n\t"\
		"shlq	$3,%%rbx			\n\t"\
		"shlq	$3,%%rcx			\n\t"\
		"shlq	$3,%%rdx			\n\t"\
		"shlq	$3,%%rdi			\n\t"\
		"addq	$0x230,%%r8	/* two */\n\t"\
		"addq	%%rax,%%rbx			\n\t"\
		"addq	%%rax,%%rcx			\n\t"\
		"addq	%%rax,%%rdx			\n\t			movslq	%[__p8],%%r10		\n\t"\
		"/* SSE2_RADIX4_DIT_0TWIDDLE_B(r1) */\n\t	shlq	$3,%%r10			\n\t"\
		"movq	%[__r1],%%rsi			\n\t		movq	%%r10,%%r11			\n\t"\
		"movaps	    (%%rax),%%xmm2		\n\t		movq	%%r10,%%r12			\n\t"\
		"movaps	    (%%rcx),%%xmm6		\n\t		movq	%%r10,%%r13			\n\t"\
		"movaps	0x10(%%rax),%%xmm3		\n\t		addq	%%rax,%%r10			\n\t"\
		"movaps	0x10(%%rcx),%%xmm7		\n\t		addq	%%rbx,%%r11			\n\t"\
		"movaps	    (%%rbx),%%xmm0		\n\t		addq	%%rcx,%%r12			\n\t"\
		"movaps	    (%%rdx),%%xmm4		\n\t		addq	%%rdx,%%r13			\n\t"\
		"movaps	0x10(%%rbx),%%xmm1		\n\t		/* SSE2_RADIX4_DIT_0TWIDDLE_B(r17) */\n\t"\
		"movaps	0x10(%%rdx),%%xmm5		\n\t		leaq	0x100(%%rsi),%%r14	/* r17 */\n\t"\
		"subpd	%%xmm0,%%xmm2			\n\t		movaps	    (%%r10),%%xmm10	\n\t"\
		"subpd	%%xmm4,%%xmm6			\n\t		movaps	    (%%r12),%%xmm14	\n\t"\
		"subpd	%%xmm1,%%xmm3			\n\t		movaps	0x10(%%r10),%%xmm11	\n\t"\
		"subpd	%%xmm5,%%xmm7			\n\t		movaps	0x10(%%r12),%%xmm15	\n\t"\
		"addpd	%%xmm0,%%xmm0			\n\t		movaps	    (%%r11),%%xmm8	\n\t"\
		"addpd	%%xmm4,%%xmm4			\n\t		movaps	    (%%r13),%%xmm12	\n\t"\
		"addpd	%%xmm1,%%xmm1			\n\t		movaps	0x10(%%r11),%%xmm9	\n\t"\
		"addpd	%%xmm5,%%xmm5			\n\t		movaps	0x10(%%r13),%%xmm13	\n\t"\
		"addpd	%%xmm2,%%xmm0			\n\t		subpd	%%xmm8 ,%%xmm10		\n\t"\
		"addpd	%%xmm6,%%xmm4			\n\t		subpd	%%xmm12,%%xmm14		\n\t"\
		"addpd	%%xmm3,%%xmm1			\n\t		subpd	%%xmm9 ,%%xmm11		\n\t"\
		"addpd	%%xmm7,%%xmm5			\n\t		subpd	%%xmm13,%%xmm15		\n\t"\
		"subpd	%%xmm4,%%xmm0			\n\t		mulpd	(%%r8),%%xmm8		\n\t"\
		"subpd	%%xmm7,%%xmm2			\n\t		mulpd	(%%r8),%%xmm12		\n\t"\
		"subpd	%%xmm5,%%xmm1			\n\t		mulpd	(%%r8),%%xmm9		\n\t"\
		"subpd	%%xmm6,%%xmm3			\n\t		mulpd	(%%r8),%%xmm13		\n\t"\
		"movaps	%%xmm0,0x040(%%rsi)		\n\t		addpd	%%xmm10,%%xmm8		\n\t"\
		"movaps	%%xmm2,0x060(%%rsi)		\n\t		addpd	%%xmm14,%%xmm12		\n\t"\
		"movaps	%%xmm1,0x050(%%rsi)		\n\t		addpd	%%xmm11,%%xmm9		\n\t"\
		"movaps	%%xmm3,0x030(%%rsi)		\n\t		addpd	%%xmm15,%%xmm13		\n\t"\
		"mulpd	(%%r8),%%xmm4			\n\t		subpd	%%xmm12,%%xmm8		\n\t"\
		"mulpd	(%%r8),%%xmm7			\n\t		subpd	%%xmm15,%%xmm10		\n\t"\
		"mulpd	(%%r8),%%xmm5			\n\t		subpd	%%xmm13,%%xmm9		\n\t"\
		"mulpd	(%%r8),%%xmm6			\n\t		subpd	%%xmm14,%%xmm11		\n\t"\
		"addpd	%%xmm0,%%xmm4			\n\t		movaps	%%xmm8 ,0x040(%%r14)\n\t"\
		"addpd	%%xmm2,%%xmm7			\n\t		movaps	%%xmm10,0x060(%%r14)\n\t"\
		"addpd	%%xmm1,%%xmm5			\n\t		movaps	%%xmm9 ,0x050(%%r14)\n\t"\
		"addpd	%%xmm3,%%xmm6			\n\t		movaps	%%xmm11,0x030(%%r14)\n\t"\
		"movaps	%%xmm4,     (%%rsi)		\n\t		addpd	%%xmm12,%%xmm12		\n\t"\
		"movaps	%%xmm7,0x020(%%rsi)		\n\t		addpd	%%xmm15,%%xmm15		\n\t"\
		"movaps	%%xmm5,0x010(%%rsi)		\n\t		addpd	%%xmm13,%%xmm13		\n\t"\
		"movaps	%%xmm6,0x070(%%rsi)		\n\t		addpd	%%xmm14,%%xmm14		\n\t"\
		"addq	%%rdi,%%rax				\n\t		addpd	%%xmm8,%%xmm12		\n\t"\
		"addq	%%rdi,%%rbx				\n\t		addpd	%%xmm10,%%xmm15		\n\t"\
		"addq	%%rdi,%%rcx				\n\t		addpd	%%xmm9,%%xmm13		\n\t"\
		"addq	%%rdi,%%rdx				\n\t		addpd	%%xmm11,%%xmm14		\n\t"\
		"/* SSE2_RADIX4_DIT_0TWIDDLE_B(r9) */\n\t	movaps	%%xmm12,     (%%r14)\n\t"\
		"addq	$0x80,%%rsi		/* r9 */\n\t		movaps	%%xmm15,0x020(%%r14)\n\t"\
		"movaps	    (%%rax),%%xmm2		\n\t		movaps	%%xmm13,0x010(%%r14)\n\t"\
		"movaps	    (%%rcx),%%xmm6		\n\t		movaps	%%xmm14,0x070(%%r14)\n\t"\
		"movaps	0x10(%%rax),%%xmm3		\n\t		addq	%%rdi,%%r10			\n\t"\
		"movaps	0x10(%%rcx),%%xmm7		\n\t		addq	%%rdi,%%r11			\n\t"\
		"movaps	    (%%rbx),%%xmm0		\n\t		addq	%%rdi,%%r12			\n\t"\
		"movaps	    (%%rdx),%%xmm4		\n\t		addq	%%rdi,%%r13			\n\t"\
		"movaps	0x10(%%rbx),%%xmm1		\n\t		/* SSE2_RADIX4_DIT_0TWIDDLE_B(r25) */\n\t"\
		"movaps	0x10(%%rdx),%%xmm5		\n\t		leaq	0x100(%%rsi),%%r14	/* r25 */\n\t"\
		"subpd	%%xmm0,%%xmm2			\n\t		movaps	    (%%r10),%%xmm10	\n\t"\
		"subpd	%%xmm4,%%xmm6			\n\t		movaps	    (%%r12),%%xmm14	\n\t"\
		"subpd	%%xmm1,%%xmm3			\n\t		movaps	0x10(%%r10),%%xmm11	\n\t"\
		"subpd	%%xmm5,%%xmm7			\n\t		movaps	0x10(%%r12),%%xmm15	\n\t"\
		"addpd	%%xmm0,%%xmm0			\n\t		movaps	    (%%r11),%%xmm8	\n\t"\
		"addpd	%%xmm4,%%xmm4			\n\t		movaps	    (%%r13),%%xmm12	\n\t"\
		"addpd	%%xmm1,%%xmm1			\n\t		movaps	0x10(%%r11),%%xmm9	\n\t"\
		"addpd	%%xmm5,%%xmm5			\n\t		movaps	0x10(%%r13),%%xmm13	\n\t"\
		"addpd	%%xmm2,%%xmm0			\n\t		subpd	%%xmm8 ,%%xmm10		\n\t"\
		"addpd	%%xmm6,%%xmm4			\n\t		subpd	%%xmm12,%%xmm14		\n\t"\
		"addpd	%%xmm3,%%xmm1			\n\t		subpd	%%xmm9 ,%%xmm11		\n\t"\
		"addpd	%%xmm7,%%xmm5			\n\t		subpd	%%xmm13,%%xmm15		\n\t"\
		"subpd	%%xmm4,%%xmm0			\n\t		mulpd	(%%r8),%%xmm8		\n\t"\
		"subpd	%%xmm7,%%xmm2			\n\t		mulpd	(%%r8),%%xmm12		\n\t"\
		"subpd	%%xmm5,%%xmm1			\n\t		mulpd	(%%r8),%%xmm9		\n\t"\
		"subpd	%%xmm6,%%xmm3			\n\t		mulpd	(%%r8),%%xmm13		\n\t"\
		"movaps	%%xmm0,0x040(%%rsi)		\n\t		addpd	%%xmm10,%%xmm8		\n\t"\
		"movaps	%%xmm2,0x060(%%rsi)		\n\t		addpd	%%xmm14,%%xmm12		\n\t"\
		"movaps	%%xmm1,0x050(%%rsi)		\n\t		addpd	%%xmm11,%%xmm9		\n\t"\
		"movaps	%%xmm3,0x030(%%rsi)		\n\t		addpd	%%xmm15,%%xmm13		\n\t"\
		"mulpd	(%%r8),%%xmm4			\n\t		subpd	%%xmm12,%%xmm8		\n\t"\
		"mulpd	(%%r8),%%xmm7			\n\t		subpd	%%xmm15,%%xmm10		\n\t"\
		"mulpd	(%%r8),%%xmm5			\n\t		subpd	%%xmm13,%%xmm9		\n\t"\
		"mulpd	(%%r8),%%xmm6			\n\t		subpd	%%xmm14,%%xmm11		\n\t"\
		"addpd	%%xmm0,%%xmm4			\n\t		movaps	%%xmm8 ,0x040(%%r14)\n\t"\
		"addpd	%%xmm2,%%xmm7			\n\t		movaps	%%xmm10,0x060(%%r14)\n\t"\
		"addpd	%%xmm1,%%xmm5			\n\t		movaps	%%xmm9 ,0x050(%%r14)\n\t"\
		"addpd	%%xmm3,%%xmm6			\n\t		movaps	%%xmm11,0x030(%%r14)\n\t"\
		"movaps	%%xmm4,     (%%rsi)		\n\t		addpd	%%xmm12,%%xmm12		\n\t"\
		"movaps	%%xmm7,0x020(%%rsi)		\n\t		addpd	%%xmm15,%%xmm15		\n\t"\
		"movaps	%%xmm5,0x010(%%rsi)		\n\t		addpd	%%xmm13,%%xmm13		\n\t"\
		"movaps	%%xmm6,0x070(%%rsi)		\n\t		addpd	%%xmm14,%%xmm14		\n\t"\
		"											addpd	%%xmm8 ,%%xmm12		\n\t"\
		"											addpd	%%xmm10,%%xmm15		\n\t"\
		"											addpd	%%xmm9 ,%%xmm13		\n\t"\
		"											addpd	%%xmm11,%%xmm14		\n\t"\
		"											movaps	%%xmm12,     (%%r14)\n\t"\
		"											movaps	%%xmm15,0x020(%%r14)\n\t"\
		"											movaps	%%xmm13,0x010(%%r14)\n\t"\
		"											movaps	%%xmm14,0x070(%%r14)\n\t"\
		"/***************************************************************************************************/\n\t"\
		"/*..and now do four more radix-4 transforms, including the internal and external twiddle factors:  */\n\t"\
		"/***************************************************************************************************/\n\t"\
	"/*...Block 1: 24 MOVapd, 26 ADD/SUBpd, 12 MULpd */	/*...Block 3: 31 MOVapd, 32 ADD/SUBpd, 20 MULpd */\n\t"\
		"movq	%[__add0],%%rax							\n\t	movslq	%[__p2],%%r10			\n\t"\
		"movslq	%[__p4],%%rbx							\n\t	leaq	-0x40(%%rsi),%%r12	/* r5  */\n\t"\
		"movslq	%[__p8],%%rdi							\n\t	leaq	 0x40(%%rsi),%%r13	/* r13 */\n\t"\
		"shlq	$3,%%rbx								\n\t	movq	%%rbx,%%r11		/* p4 */\n\t"\
		"shlq	$3,%%rdi								\n\t	shlq	$3,%%r10		/* p2 */\n\t"\
		"movq	%[__r1],%%rcx							\n\t	movaps	-0x230(%%r8),%%xmm10	/* isrt2 */\n\t"\
		"movq	%%rsi,%%rdx		/* r9 */				\n\t	movaps	0x100(%%r12),%%xmm12	\n\t"\
		"			addq	%%rax,%%r10	/* a[j+p2 ] */	\n\t	movaps	0x110(%%r12),%%xmm13	\n\t"\
		"addq	%%rax,%%rbx	/* a[j+p4 ] */				\n\t	movaps	0x100(%%r13),%%xmm14	\n\t"\
		"			addq	%%r10,%%r11	/* a[j+p6 ] */	\n\t	movaps	0x110(%%r13),%%xmm15	\n\t"\
		"movaps	     (%%rdx),%%xmm2						\n\t	mulpd	%%xmm10,%%xmm12			\n\t"\
		"movaps	0x100(%%rdx),%%xmm4						\n\t	mulpd	%%xmm10,%%xmm13			\n\t"\
		"movaps	0x010(%%rdx),%%xmm3						\n\t	mulpd	%%xmm10,%%xmm14			\n\t"\
		"movaps	0x110(%%rdx),%%xmm5						\n\t	mulpd	%%xmm10,%%xmm15			\n\t"\
		"movaps	     (%%rcx),%%xmm0						\n\t	movaps	     (%%r12),%%xmm8		\n\t"\
		"movaps	0x100(%%rcx),%%xmm6						\n\t	movaps	0x010(%%r13),%%xmm10	\n\t"\
		"movaps	0x010(%%rcx),%%xmm1						\n\t	movaps	0x010(%%r12),%%xmm11	\n\t"\
		"movaps	0x110(%%rcx),%%xmm7						\n\t	movaps	     (%%r13),%%xmm9		\n\t"\
		"subpd	%%xmm2,%%xmm0							\n\t	subpd	%%xmm12,%%xmm13			\n\t"\
		"subpd	%%xmm4,%%xmm6							\n\t	subpd	%%xmm15,%%xmm14			\n\t"\
		"subpd	%%xmm3,%%xmm1							\n\t	mulpd	(%%r8),%%xmm12			\n\t"\
		"subpd	%%xmm5,%%xmm7							\n\t	mulpd	(%%r8),%%xmm15			\n\t"\
		"mulpd	(%%r8),%%xmm2							\n\t	addpd	%%xmm13,%%xmm12			\n\t"\
		"mulpd	(%%r8),%%xmm4							\n\t	addpd	%%xmm14,%%xmm15			\n\t"\
		"mulpd	(%%r8),%%xmm3							\n\t	subpd	%%xmm14,%%xmm12			\n\t"\
		"mulpd	(%%r8),%%xmm5							\n\t	subpd	%%xmm10,%%xmm8			\n\t"\
		"addpd	%%xmm0,%%xmm2							\n\t	subpd	%%xmm15,%%xmm13			\n\t"\
		"addpd	%%xmm6,%%xmm4							\n\t	subpd	%%xmm9 ,%%xmm11			\n\t"\
		"addpd	%%xmm1,%%xmm3							\n\t	mulpd	(%%r8),%%xmm14			\n\t"\
		"addpd	%%xmm7,%%xmm5							\n\t	mulpd	(%%r8),%%xmm10			\n\t"\
		"/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c0) */	\n\t	mulpd	(%%r8),%%xmm15			\n\t"\
		"addq	$0x1b0,%%rsi	/* c0, from r9 */		\n\t	mulpd	(%%r8),%%xmm9			\n\t"\
		"subpd	%%xmm4,%%xmm2							\n\t	addpd	%%xmm12,%%xmm14			\n\t"\
		"subpd	%%xmm7,%%xmm0							\n\t	addpd	%%xmm8 ,%%xmm10			\n\t"\
		"subpd	%%xmm5,%%xmm3							\n\t	addpd	%%xmm13,%%xmm15			\n\t"\
		"subpd	%%xmm6,%%xmm1							\n\t	addpd	%%xmm11,%%xmm9			\n\t"\
		"mulpd	(%%r8),%%xmm4							\n\t	/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c2) */\n\t"\
		"mulpd	(%%r8),%%xmm7							\n\t	addq	$0x130,%%r14	/* c2, from r25 */\n\t"\
		"mulpd	(%%r8),%%xmm5							\n\t	subpd	%%xmm12,%%xmm10			\n\t"\
		"mulpd	(%%r8),%%xmm6							\n\t	subpd	%%xmm15,%%xmm8			\n\t"\
		"addpd	%%xmm2,%%xmm4							\n\t	subpd	%%xmm13,%%xmm11			\n\t"\
		"addpd	%%xmm0,%%xmm7							\n\t	subpd	%%xmm14,%%xmm9			\n\t"\
		"addpd	%%xmm3,%%xmm5							\n\t	mulpd	(%%r8),%%xmm12			\n\t"\
		"addpd	%%xmm1,%%xmm6							\n\t	mulpd	(%%r8),%%xmm15			\n\t"\
		"movaps	%%xmm2,     (%%rcx)						\n\t	addpd	%%xmm13,%%xmm13			\n\t"\
		"movaps	%%xmm0,0x010(%%rdx)						\n\t	addpd	%%xmm14,%%xmm14			\n\t"\
		"movaps	%%xmm3,0x010(%%rcx)						\n\t	addpd	%%xmm10,%%xmm12			\n\t"\
		"movaps	%%xmm6,     (%%rdx)						\n\t	addpd	%%xmm8 ,%%xmm15			\n\t"\
		"movaps	%%xmm4,%%xmm2							\n\t	addpd	%%xmm11,%%xmm13			\n\t"\
		"movaps	%%xmm7,%%xmm0							\n\t	addpd	%%xmm9 ,%%xmm14			\n\t"\
		"movaps	%%xmm5,%%xmm3							\n\t	movaps	%%xmm10,     (%%r12)	\n\t"\
		"movaps	%%xmm1,%%xmm6							\n\t	movaps	%%xmm8 ,0x010(%%r13)	\n\t"\
		"mulpd	    (%%rsi),%%xmm4						\n\t	movaps	%%xmm11,0x010(%%r12)	\n\t"\
		"mulpd	0x20(%%rsi),%%xmm7						\n\t	movaps	%%xmm14,     (%%r13)	\n\t"\
		"mulpd	    (%%rsi),%%xmm5						\n\t	movaps	%%xmm12,%%xmm10			\n\t"\
		"mulpd	0x20(%%rsi),%%xmm1						\n\t	movaps	%%xmm15,%%xmm8			\n\t"\
		"mulpd	0x10(%%rsi),%%xmm2						\n\t	movaps	%%xmm13,%%xmm11			\n\t"\
		"mulpd	0x30(%%rsi),%%xmm0						\n\t	movaps	%%xmm9 ,%%xmm14			\n\t"\
		"mulpd	0x10(%%rsi),%%xmm3						\n\t	mulpd	    (%%r14),%%xmm12		\n\t"\
		"mulpd	0x30(%%rsi),%%xmm6						\n\t	mulpd	0x20(%%r14),%%xmm15		\n\t"\
		"subpd	%%xmm2,%%xmm5							\n\t	mulpd	    (%%r14),%%xmm13		\n\t"\
		"subpd	%%xmm0,%%xmm1							\n\t	mulpd	0x20(%%r14),%%xmm9		\n\t"\
		"addpd	%%xmm3,%%xmm4							\n\t	mulpd	0x10(%%r14),%%xmm10		\n\t"\
		"addpd	%%xmm6,%%xmm7							\n\t	mulpd	0x30(%%r14),%%xmm8		\n\t"\
		"movaps	%%xmm5,0x10(%%rax)						\n\t	mulpd	0x10(%%r14),%%xmm11		\n\t"\
		"movaps	%%xmm1,0x10(%%rbx)						\n\t	mulpd	0x30(%%r14),%%xmm14		\n\t"\
		"movaps	%%xmm4,    (%%rax)						\n\t	subpd	%%xmm10,%%xmm13			\n\t"\
		"movaps	%%xmm7,    (%%rbx)						\n\t	subpd	%%xmm8 ,%%xmm9			\n\t"\
		"addq	%%rdi,%%rax	/* a[j+p8 ] */				\n\t	addpd	%%xmm11,%%xmm12			\n\t"\
		"addq	%%rdi,%%rbx	/* a[j+p12] */				\n\t	addpd	%%xmm14,%%xmm15			\n\t"\
		"addq	$0x40,%%rsi	/* c4 */					\n\t	movaps	%%xmm13,0x10(%%r10)		\n\t"\
		"movaps	     (%%rcx),%%xmm4						\n\t	movaps	%%xmm9 ,0x10(%%r11)		\n\t"\
		"movaps	0x010(%%rdx),%%xmm0						\n\t	movaps	%%xmm12,    (%%r10)		\n\t"\
		"movaps	0x010(%%rcx),%%xmm5						\n\t	movaps	%%xmm15,    (%%r11)		\n\t"\
		"movaps	     (%%rdx),%%xmm6						\n\t	addq	%%rdi,%%r10	/* a[j+p10] */	\n\t"\
		"movaps	%%xmm4,%%xmm2							\n\t	addq	%%rdi,%%r11	/* a[j+p14] */	\n\t"\
		"movaps	%%xmm0,%%xmm1							\n\t	addq	$0x40,%%r14	/* c2 */	\n\t"\
		"movaps	%%xmm5,%%xmm3							\n\t	movaps	     (%%r12),%%xmm12	\n\t"\
		"movaps	%%xmm6,%%xmm7							\n\t	movaps	0x010(%%r13),%%xmm8		\n\t"\
		"mulpd	    (%%rsi),%%xmm4						\n\t	movaps	0x010(%%r12),%%xmm13	\n\t"\
		"mulpd	0x20(%%rsi),%%xmm0						\n\t	movaps	     (%%r13),%%xmm14	\n\t"\
		"mulpd	    (%%rsi),%%xmm5						\n\t	movaps	%%xmm12,%%xmm10			\n\t"\
		"mulpd	0x20(%%rsi),%%xmm6						\n\t	movaps	%%xmm8 ,%%xmm9			\n\t"\
		"mulpd	0x10(%%rsi),%%xmm2						\n\t	movaps	%%xmm13,%%xmm11			\n\t"\
		"mulpd	0x30(%%rsi),%%xmm1						\n\t	movaps	%%xmm14,%%xmm15			\n\t"\
		"mulpd	0x10(%%rsi),%%xmm3						\n\t	mulpd	    (%%r14),%%xmm12		\n\t"\
		"mulpd	0x30(%%rsi),%%xmm7						\n\t	mulpd	0x20(%%r14),%%xmm8		\n\t"\
		"subpd	%%xmm2,%%xmm5							\n\t	mulpd	    (%%r14),%%xmm13		\n\t"\
		"subpd	%%xmm1,%%xmm6							\n\t	mulpd	0x20(%%r14),%%xmm14		\n\t"\
		"addpd	%%xmm3,%%xmm4							\n\t	mulpd	0x10(%%r14),%%xmm10		\n\t"\
		"addpd	%%xmm7,%%xmm0							\n\t	mulpd	0x30(%%r14),%%xmm9		\n\t"\
		"movaps	%%xmm5,0x10(%%rax)						\n\t	mulpd	0x10(%%r14),%%xmm11		\n\t"\
		"movaps	%%xmm6,0x10(%%rbx)						\n\t	mulpd	0x30(%%r14),%%xmm15		\n\t"\
		"movaps	%%xmm4,    (%%rax)						\n\t	subpd	%%xmm10,%%xmm13			\n\t"\
		"movaps	%%xmm0,    (%%rbx)						\n\t	subpd	%%xmm9 ,%%xmm14			\n\t"\
		"/*...Block 2: 34 MOVapd, 36 ADD/SUBpd, 26 MULpd */\n\t	addpd	%%xmm11,%%xmm12			\n\t"\
		"movslq	%[__p1],%%r14							\n\t	addpd	%%xmm15,%%xmm8			\n\t"\
		"subq	%%rdi,%%rax	/* a[j+p0 ] */				\n\t	movaps	%%xmm13,0x10(%%r10)		\n\t"\
		"subq	%%rdi,%%rbx	/* a[j+p4 ] */				\n\t	movaps	%%xmm14,0x10(%%r11)		\n\t"\
		"shlq	$3,%%r14								\n\t	movaps	%%xmm12,    (%%r10)		\n\t"\
		"addq	%%r14,%%rax	/* a[j+p1 ] */				\n\t	movaps	%%xmm8 ,    (%%r11)		\n\t"\
		"addq	%%r14,%%rbx	/* a[j+p5 ] */				\n\t/*...Block 4: 34 MOVapd, 36 ADD/SUBpd, 26 MULpd */\n\t"\
		"addq	$0x20,%%rcx			/* r3 , from r1 */	\n\t	subq	%%rdi,%%r10	/* a[j+p2 ] */	\n\t"\
		"leaq	0x80(%%rcx),%%rdx	/* r11, from r3 */	\n\t	subq	%%rdi,%%r11	/* a[j+p6 ] */	\n\t"\
		"subq	$0x60,%%rsi		/* cc0, from c4 */		\n\t	addq	%%r14,%%r10	/* a[j+p3 ] */	\n\t"\
		"movaps	0x100(%%rcx),%%xmm4						\n\t	addq	%%r14,%%r11	/* a[j+p7 ] */	\n\t"\
		"movaps	0x110(%%rcx),%%xmm5						\n\t	addq	$0x20,%%r12	/* r7 , from r5  */\n\t"\
		"movaps	     (%%rsi),%%xmm2						\n\t	addq	$0x20,%%r13	/* r15, from r13 */\n\t"\
		"movaps	0x010(%%rsi),%%xmm3						\n\t	addq	$0x100,%%r12			\n\t"\
		"movaps	%%xmm4,%%xmm6							\n\t	addq	$0x100,%%r13			\n\t"\
		"movaps	%%xmm5,%%xmm7							\n\t	movaps	     (%%r12),%%xmm12	\n\t"\
		"mulpd	%%xmm2,%%xmm4							\n\t	movaps	0x010(%%r12),%%xmm13	\n\t"\
		"mulpd	%%xmm2,%%xmm5							\n\t	movaps	     (%%rsi),%%xmm11	\n\t"\
		"mulpd	%%xmm3,%%xmm6							\n\t	movaps	0x010(%%rsi),%%xmm10	\n\t"\
		"movaps	0x100(%%rdx),%%xmm0						\n\t	movaps	%%xmm12,%%xmm14			\n\t"\
		"mulpd	%%xmm3,%%xmm7							\n\t	movaps	%%xmm13,%%xmm15			\n\t"\
		"movaps	0x110(%%rdx),%%xmm1						\n\t	mulpd	%%xmm10,%%xmm12			\n\t"\
		"subpd	%%xmm6,%%xmm5							\n\t	mulpd	%%xmm10,%%xmm13			\n\t"\
		"movaps	%%xmm0,%%xmm6							\n\t	mulpd	%%xmm11,%%xmm14			\n\t"\
		"addpd	%%xmm7,%%xmm4							\n\t	movaps	     (%%r13),%%xmm8		\n\t"\
		"movaps	%%xmm1,%%xmm7							\n\t	mulpd	%%xmm11,%%xmm15			\n\t"\
		"mulpd	%%xmm3,%%xmm0							\n\t	movaps	0x010(%%r13),%%xmm9		\n\t"\
		"mulpd	%%xmm3,%%xmm1							\n\t	subpd	%%xmm14,%%xmm13			\n\t"\
		"mulpd	%%xmm2,%%xmm6							\n\t	movaps	%%xmm8 ,%%xmm14			\n\t"\
		"mulpd	%%xmm2,%%xmm7							\n\t	addpd	%%xmm15,%%xmm12			\n\t"\
		"subpd	%%xmm6,%%xmm1							\n\t	movaps	%%xmm9 ,%%xmm15			\n\t"\
		"addpd	%%xmm7,%%xmm0							\n\t	mulpd	%%xmm11,%%xmm8			\n\t"\
		"movaps	%%xmm5,%%xmm7							\n\t	mulpd	%%xmm11,%%xmm9			\n\t"\
		"movaps	%%xmm4,%%xmm6							\n\t	mulpd	%%xmm10,%%xmm14			\n\t"\
		"addpd	%%xmm0,%%xmm4							\n\t	mulpd	%%xmm10,%%xmm15			\n\t"\
		"addpd	%%xmm1,%%xmm5							\n\t	subpd	%%xmm14,%%xmm9			\n\t"\
		"subpd	%%xmm0,%%xmm6							\n\t	addpd	%%xmm15,%%xmm8			\n\t"\
		"subpd	%%xmm1,%%xmm7							\n\t	movaps	%%xmm13,%%xmm15			\n\t"\
		"movaps	     (%%rdx),%%xmm2						\n\t	movaps	%%xmm12,%%xmm14			\n\t"\
		"movaps	0x010(%%rdx),%%xmm3						\n\t	addpd	%%xmm8 ,%%xmm14			\n\t"\
		"movaps	-0x230(%%r8),%%xmm1						\n\t	addpd	%%xmm9 ,%%xmm15			\n\t"\
		"movaps	%%xmm3,%%xmm0							\n\t	subpd	%%xmm8 ,%%xmm12			\n\t"\
		"subpd	%%xmm2,%%xmm3							\n\t	subpd	%%xmm9 ,%%xmm13			\n\t"\
		"addpd	%%xmm0,%%xmm2							\n\t	subq	$0x100,%%r12			\n\t"\
		"mulpd	%%xmm1,%%xmm2							\n\t	subq	$0x100,%%r13			\n\t"\
		"mulpd	%%xmm1,%%xmm3							\n\t	movaps	     (%%r13),%%xmm8		\n\t"\
		"movaps	     (%%rcx),%%xmm0						\n\t	movaps	0x010(%%r13),%%xmm9		\n\t"\
		"movaps	0x010(%%rcx),%%xmm1						\n\t	movaps	-0x230(%%r8),%%xmm11	\n\t"\
		"subpd	%%xmm2,%%xmm0							\n\t	movaps	%%xmm8 ,%%xmm10			\n\t"\
		"subpd	%%xmm3,%%xmm1							\n\t	subpd	%%xmm9 ,%%xmm8			\n\t"\
		"mulpd	(%%r8),%%xmm2							\n\t	addpd	%%xmm10,%%xmm9			\n\t"\
		"mulpd	(%%r8),%%xmm3							\n\t	mulpd	%%xmm11,%%xmm8			\n\t"\
		"addpd	%%xmm0,%%xmm2							\n\t	mulpd	%%xmm11,%%xmm9			\n\t"\
		"addpd	%%xmm1,%%xmm3							\n\t	movaps	     (%%r12),%%xmm10	\n\t"\
		"/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c1) */	\n\t	movaps	0x010(%%r12),%%xmm11	\n\t"\
		"addq	$0x120,%%rsi	/* c1, from cc0 */		\n\t	subpd	%%xmm8 ,%%xmm10			\n\t"\
		"subpd	%%xmm4,%%xmm2							\n\t	subpd	%%xmm9 ,%%xmm11			\n\t"\
		"subpd	%%xmm7,%%xmm0							\n\t	mulpd	(%%r8),%%xmm8			\n\t"\
		"subpd	%%xmm5,%%xmm3							\n\t	mulpd	(%%r8),%%xmm9			\n\t"\
		"subpd	%%xmm6,%%xmm1							\n\t	addpd	%%xmm10,%%xmm8			\n\t"\
		"mulpd	(%%r8),%%xmm4							\n\t	addpd	%%xmm11,%%xmm9			\n\t"\
		"mulpd	(%%r8),%%xmm7							\n\t/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c3) */\n\t"\
		"mulpd	(%%r8),%%xmm5							\n\t	leaq	0x80(%%rsi),%%r14	/* c3, from c1 */\n\t"\
		"mulpd	(%%r8),%%xmm6							\n\t	subpd	%%xmm12,%%xmm10			\n\t"\
		"addpd	%%xmm2,%%xmm4							\n\t	subpd	%%xmm15,%%xmm8			\n\t"\
		"addpd	%%xmm0,%%xmm7							\n\t	subpd	%%xmm13,%%xmm11			\n\t"\
		"addpd	%%xmm3,%%xmm5							\n\t	subpd	%%xmm14,%%xmm9			\n\t"\
		"addpd	%%xmm1,%%xmm6							\n\t	addpd	%%xmm12,%%xmm12			\n\t"\
		"movaps	%%xmm2,     (%%rcx)						\n\t	addpd	%%xmm15,%%xmm15			\n\t"\
		"movaps	%%xmm0,0x010(%%rdx)						\n\t	addpd	%%xmm13,%%xmm13			\n\t"\
		"movaps	%%xmm3,0x010(%%rcx)						\n\t	addpd	%%xmm14,%%xmm14			\n\t"\
		"movaps	%%xmm6,     (%%rdx)						\n\t	addpd	%%xmm10,%%xmm12			\n\t"\
		"movaps	%%xmm4,%%xmm2							\n\t	addpd	%%xmm8 ,%%xmm15			\n\t"\
		"movaps	%%xmm7,%%xmm0							\n\t	addpd	%%xmm11,%%xmm13			\n\t"\
		"movaps	%%xmm5,%%xmm3							\n\t	addpd	%%xmm9 ,%%xmm14			\n\t"\
		"movaps	%%xmm1,%%xmm6							\n\t	movaps	%%xmm10,     (%%r12)	\n\t"\
		"mulpd	    (%%rsi),%%xmm4						\n\t	movaps	%%xmm8 ,0x010(%%r13)	\n\t"\
		"mulpd	0x20(%%rsi),%%xmm7						\n\t	movaps	%%xmm11,0x010(%%r12)	\n\t"\
		"mulpd	    (%%rsi),%%xmm5						\n\t	movaps	%%xmm14,     (%%r13)	\n\t"\
		"mulpd	0x20(%%rsi),%%xmm1						\n\t	movaps	%%xmm12,%%xmm10			\n\t"\
		"mulpd	0x10(%%rsi),%%xmm2						\n\t	movaps	%%xmm15,%%xmm8			\n\t"\
		"mulpd	0x30(%%rsi),%%xmm0						\n\t	movaps	%%xmm13,%%xmm11			\n\t"\
		"mulpd	0x10(%%rsi),%%xmm3						\n\t	movaps	%%xmm9 ,%%xmm14			\n\t"\
		"mulpd	0x30(%%rsi),%%xmm6						\n\t	mulpd	    (%%r14),%%xmm12		\n\t"\
		"subpd	%%xmm2,%%xmm5							\n\t	mulpd	0x20(%%r14),%%xmm15		\n\t"\
		"subpd	%%xmm0,%%xmm1							\n\t	mulpd	    (%%r14),%%xmm13		\n\t"\
		"addpd	%%xmm3,%%xmm4							\n\t	mulpd	0x20(%%r14),%%xmm9		\n\t"\
		"addpd	%%xmm6,%%xmm7							\n\t	mulpd	0x10(%%r14),%%xmm10		\n\t"\
		"movaps	%%xmm5,0x10(%%rax)						\n\t	mulpd	0x30(%%r14),%%xmm8		\n\t"\
		"movaps	%%xmm1,0x10(%%rbx)						\n\t	mulpd	0x10(%%r14),%%xmm11		\n\t"\
		"movaps	%%xmm4,    (%%rax)						\n\t	mulpd	0x30(%%r14),%%xmm14		\n\t"\
		"movaps	%%xmm7,    (%%rbx)						\n\t	subpd	%%xmm10,%%xmm13			\n\t"\
		"addq	%%rdi,%%rax								\n\t	subpd	%%xmm8 ,%%xmm9			\n\t"\
		"addq	%%rdi,%%rbx								\n\t	addpd	%%xmm11,%%xmm12			\n\t"\
		"addq	$0x40,%%rsi	/* c2 */					\n\t	addpd	%%xmm14,%%xmm15			\n\t"\
		"movaps	     (%%rcx),%%xmm4						\n\t	movaps	%%xmm13,0x10(%%r10)		\n\t"\
		"movaps	0x010(%%rdx),%%xmm0						\n\t	movaps	%%xmm9 ,0x10(%%r11)		\n\t"\
		"movaps	0x010(%%rcx),%%xmm5						\n\t	movaps	%%xmm12,    (%%r10)		\n\t"\
		"movaps	     (%%rdx),%%xmm6						\n\t	movaps	%%xmm15,    (%%r11)		\n\t"\
		"movaps	%%xmm4,%%xmm2							\n\t	addq	%%rdi,%%r10				\n\t"\
		"movaps	%%xmm0,%%xmm1							\n\t	addq	%%rdi,%%r11				\n\t"\
		"movaps	%%xmm5,%%xmm3							\n\t	addq	$0x40,%%r14	/* c2 */	\n\t"\
		"movaps	%%xmm6,%%xmm7							\n\t	movaps	     (%%r12),%%xmm12	\n\t"\
		"mulpd	    (%%rsi),%%xmm4						\n\t	movaps	0x010(%%r13),%%xmm8		\n\t"\
		"mulpd	0x20(%%rsi),%%xmm0						\n\t	movaps	0x010(%%r12),%%xmm13	\n\t"\
		"mulpd	    (%%rsi),%%xmm5						\n\t	movaps	     (%%r13),%%xmm14	\n\t"\
		"mulpd	0x20(%%rsi),%%xmm6						\n\t	movaps	%%xmm12,%%xmm10			\n\t"\
		"mulpd	0x10(%%rsi),%%xmm2						\n\t	movaps	%%xmm8 ,%%xmm9			\n\t"\
		"mulpd	0x30(%%rsi),%%xmm1						\n\t	movaps	%%xmm13,%%xmm11			\n\t"\
		"mulpd	0x10(%%rsi),%%xmm3						\n\t	movaps	%%xmm14,%%xmm15			\n\t"\
		"mulpd	0x30(%%rsi),%%xmm7						\n\t	mulpd	    (%%r14),%%xmm12		\n\t"\
		"subpd	%%xmm2,%%xmm5							\n\t	mulpd	0x20(%%r14),%%xmm8		\n\t"\
		"subpd	%%xmm1,%%xmm6							\n\t	mulpd	    (%%r14),%%xmm13		\n\t"\
		"addpd	%%xmm3,%%xmm4							\n\t	mulpd	0x20(%%r14),%%xmm14		\n\t"\
		"addpd	%%xmm7,%%xmm0							\n\t	mulpd	0x10(%%r14),%%xmm10		\n\t"\
		"movaps	%%xmm5,0x10(%%rax)						\n\t	mulpd	0x30(%%r14),%%xmm9		\n\t"\
		"movaps	%%xmm6,0x10(%%rbx)						\n\t	mulpd	0x10(%%r14),%%xmm11		\n\t"\
		"movaps	%%xmm4,    (%%rax)						\n\t	mulpd	0x30(%%r14),%%xmm15		\n\t"\
		"movaps	%%xmm0,    (%%rbx)						\n\t	subpd	%%xmm10,%%xmm13			\n\t"\
		"												\n\t	subpd	%%xmm9 ,%%xmm14			\n\t"\
		"												\n\t	addpd	%%xmm11,%%xmm12			\n\t"\
		"												\n\t	addpd	%%xmm15,%%xmm8			\n\t"\
		"												\n\t	movaps	%%xmm13,0x10(%%r10)		\n\t"\
		"												\n\t	movaps	%%xmm14,0x10(%%r11)		\n\t"\
		"												\n\t	movaps	%%xmm12,    (%%r10)		\n\t"\
		"												\n\t	movaps	%%xmm8 ,    (%%r11)		\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__p1] "m" (Xp1)\
		 ,[__p2] "m" (Xp2)\
		 ,[__p3] "m" (Xp3)\
		 ,[__p4] "m" (Xp4)\
		 ,[__p8] "m" (Xp8)\
		 ,[__r1] "m" (Xr1)\
		 ,[__isrt2] "m" (Xisrt2)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r8","r10","r11","r12","r13","r14","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
	);\
	}

  #endif // USE_64BIT_ASM_STYLE

#endif	/* radix16_dif_dit_pass_gcc_h_included */

