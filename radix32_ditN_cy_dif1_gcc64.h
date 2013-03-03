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
#ifndef radix32_ditN_cy_dif1_gcc_h_included
#define radix32_ditN_cy_dif1_gcc_h_included

	#define SSE2_RADIX32_DIT_NOTWIDDLE(Xadd,Xp01,Xp02,Xp03,Xp04,Xp08,Xp10,Xp18,Xr00,Xisrt2,Xcc0)\
	{\
	__asm__ volatile (\
		"/*...Block 1: */\n\t"\
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
		"/* SSE2_RADIX4_DIT_0TWIDDLE(add0, add1, add2, add3, r00): */\n\t"\
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
		"/* SSE2_RADIX4_DIT_0TWIDDLE_2NDOFTWO(add4, add5, add6, add7, r08): */\n\t"\
		"addq	$0x80,%%rsi			\n\t"\
		"movslq	%[__p04],%%rdi		\n\t"\
		"shlq	$3,%%rdi			\n\t"\
		"addq	%%rdi,%%rax			\n\t"\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	%%rdi,%%rcx			\n\t"\
		"addq	%%rdi,%%rdx			\n\t"\
		"\n\t"\
		"movaps	    (%%rax),%%xmm2	\n\t"\
		"movaps	    (%%rcx),%%xmm6	\n\t"\
		"movaps	0x10(%%rax),%%xmm3	\n\t"\
		"movaps	0x10(%%rcx),%%xmm7	\n\t"\
		"movaps	    (%%rbx),%%xmm0	\n\t"\
		"movaps	    (%%rdx),%%xmm4	\n\t"\
		"movaps	0x10(%%rbx),%%xmm1	\n\t"\
		"movaps	0x10(%%rdx),%%xmm5	\n\t"\
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
		"movaps	(%%rdi),%%xmm5	/* ISRT2 */\n\t"\
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
		"/* SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS(r00,r02,r04,r06,r08,r0A,r0C,r0E): */\n\t"\
		"movaps		-0x080(%%rsi),%%xmm0	/* r00 */\n\t"\
		"movaps		-0x040(%%rsi),%%xmm4	/* r04 */\n\t"\
		"movaps		-0x070(%%rsi),%%xmm1	/* r01 */\n\t"\
		"movaps		-0x030(%%rsi),%%xmm5	/* r05 */\n\t"\
		"movaps		      (%%rsi),%%xmm2	/* r08 */\n\t"\
		"movaps		 0x050(%%rsi),%%xmm7	/* r0D */\n\t"\
		"movaps		 0x010(%%rsi),%%xmm3	/* r09 */\n\t"\
		"movaps		 0x040(%%rsi),%%xmm6	/* r0C */\n\t"\
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
		"movaps		-0x060(%%rsi),%%xmm0	/* r02 */\n\t"\
		"movaps		-0x020(%%rsi),%%xmm4	/* r06 */\n\t"\
		"movaps		-0x050(%%rsi),%%xmm1	/* r03 */\n\t"\
		"movaps		-0x010(%%rsi),%%xmm5	/* r07 */\n\t"\
		"movaps		 0x020(%%rsi),%%xmm2	/* r0A */\n\t"\
		"movaps		 0x070(%%rsi),%%xmm7	/* r0F */\n\t"\
		"movaps		 0x030(%%rsi),%%xmm3	/* r0B */\n\t"\
		"movaps		 0x060(%%rsi),%%xmm6	/* r0E */\n\t"\
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
		"/*...Block 2: */\n\t"\
		"addq	$0x80,%%rsi			\n\t"\
		"movslq	%[__p04],%%rdi		\n\t"\
		"shlq	$3,%%rdi			\n\t"\
		"addq	%%rdi,%%rax			\n\t"\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	%%rdi,%%rcx			\n\t"\
		"addq	%%rdi,%%rdx			\n\t"\
		"\n\t"\
		"/* SSE2_RADIX4_DIT_0TWIDDLE(add0, add1, add2, add3, r10): */\n\t"\
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
		"/* SSE2_RADIX4_DIT_0TWIDDLE_2NDOFTWO(add4, add5, add6, add7, r18): */\n\t"\
		"addq	$0x80,%%rsi			\n\t"\
		"movslq	%[__p04],%%rdi		\n\t"\
		"shlq	$3,%%rdi			\n\t"\
		"addq	%%rdi,%%rax			\n\t"\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	%%rdi,%%rcx			\n\t"\
		"addq	%%rdi,%%rdx			\n\t"\
		"\n\t"\
		"movaps	    (%%rax),%%xmm2	\n\t"\
		"movaps	    (%%rcx),%%xmm6	\n\t"\
		"movaps	0x10(%%rax),%%xmm3	\n\t"\
		"movaps	0x10(%%rcx),%%xmm7	\n\t"\
		"movaps	    (%%rbx),%%xmm0	\n\t"\
		"movaps	    (%%rdx),%%xmm4	\n\t"\
		"movaps	0x10(%%rbx),%%xmm1	\n\t"\
		"movaps	0x10(%%rdx),%%xmm5	\n\t"\
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
		"movaps	(%%rdi),%%xmm5	/* ISRT2 */\n\t"\
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
		"/* SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS(r10,r12,r14,r16,r18,r1A,r1C,r1E): */\n\t"\
		"movaps		-0x080(%%rsi),%%xmm0	/* r10 */\n\t"\
		"movaps		-0x040(%%rsi),%%xmm4	/* r14 */\n\t"\
		"movaps		-0x070(%%rsi),%%xmm1	/* r11 */\n\t"\
		"movaps		-0x030(%%rsi),%%xmm5	/* r15 */\n\t"\
		"movaps		      (%%rsi),%%xmm2	/* r18 */\n\t"\
		"movaps		 0x050(%%rsi),%%xmm7	/* r1D */\n\t"\
		"movaps		 0x010(%%rsi),%%xmm3	/* r19 */\n\t"\
		"movaps		 0x040(%%rsi),%%xmm6	/* r1C */\n\t"\
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
		"movaps		-0x060(%%rsi),%%xmm0	/* r12 */\n\t"\
		"movaps		-0x020(%%rsi),%%xmm4	/* r16 */\n\t"\
		"movaps		-0x050(%%rsi),%%xmm1	/* r13 */\n\t"\
		"movaps		-0x010(%%rsi),%%xmm5	/* r17 */\n\t"\
		"movaps		 0x020(%%rsi),%%xmm2	/* r1A */\n\t"\
		"movaps		 0x070(%%rsi),%%xmm7	/* r1F */\n\t"\
		"movaps		 0x030(%%rsi),%%xmm3	/* r1B */\n\t"\
		"movaps		 0x060(%%rsi),%%xmm6	/* r1E */\n\t"\
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
		"/*...Block 3: */\n\t"\
		"addq	$0x80,%%rsi			\n\t"\
		"movslq	%[__p04],%%rdi		\n\t"\
		"shlq	$3,%%rdi			\n\t"\
		"addq	%%rdi,%%rax			\n\t"\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	%%rdi,%%rcx			\n\t"\
		"addq	%%rdi,%%rdx			\n\t"\
		"\n\t"\
		"/* SSE2_RADIX4_DIT_0TWIDDLE(add0, add1, add2, add3, r20): */\n\t"\
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
		"/* SSE2_RADIX4_DIT_0TWIDDLE_2NDOFTWO(add4, add5, add6, add7, r28): */\n\t"\
		"addq	$0x80,%%rsi			\n\t"\
		"movslq	%[__p04],%%rdi		\n\t"\
		"shlq	$3,%%rdi			\n\t"\
		"addq	%%rdi,%%rax			\n\t"\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	%%rdi,%%rcx			\n\t"\
		"addq	%%rdi,%%rdx			\n\t"\
		"\n\t"\
		"movaps	    (%%rax),%%xmm2	\n\t"\
		"movaps	    (%%rcx),%%xmm6	\n\t"\
		"movaps	0x10(%%rax),%%xmm3	\n\t"\
		"movaps	0x10(%%rcx),%%xmm7	\n\t"\
		"movaps	    (%%rbx),%%xmm0	\n\t"\
		"movaps	    (%%rdx),%%xmm4	\n\t"\
		"movaps	0x10(%%rbx),%%xmm1	\n\t"\
		"movaps	0x10(%%rdx),%%xmm5	\n\t"\
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
		"movaps	(%%rdi),%%xmm5	/* ISRT2 */\n\t"\
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
		"/* SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS(r20,r22,r24,r26,r28,r2A,r2C,r2E): */\n\t"\
		"movaps		-0x080(%%rsi),%%xmm0	/* r20 */\n\t"\
		"movaps		-0x040(%%rsi),%%xmm4	/* r24 */\n\t"\
		"movaps		-0x070(%%rsi),%%xmm1	/* r21 */\n\t"\
		"movaps		-0x030(%%rsi),%%xmm5	/* r25 */\n\t"\
		"movaps		      (%%rsi),%%xmm2	/* r28 */\n\t"\
		"movaps		 0x050(%%rsi),%%xmm7	/* r2D */\n\t"\
		"movaps		 0x010(%%rsi),%%xmm3	/* r29 */\n\t"\
		"movaps		 0x040(%%rsi),%%xmm6	/* r2C */\n\t"\
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
		"movaps		-0x060(%%rsi),%%xmm0	/* r22 */\n\t"\
		"movaps		-0x020(%%rsi),%%xmm4	/* r26 */\n\t"\
		"movaps		-0x050(%%rsi),%%xmm1	/* r23 */\n\t"\
		"movaps		-0x010(%%rsi),%%xmm5	/* r27 */\n\t"\
		"movaps		 0x020(%%rsi),%%xmm2	/* r2A */\n\t"\
		"movaps		 0x070(%%rsi),%%xmm7	/* r2F */\n\t"\
		"movaps		 0x030(%%rsi),%%xmm3	/* r2B */\n\t"\
		"movaps		 0x060(%%rsi),%%xmm6	/* r2E */\n\t"\
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
		"/*...Block 4: */\n\t"\
		"addq	$0x80,%%rsi			\n\t"\
		"movslq	%[__p04],%%rdi		\n\t"\
		"shlq	$3,%%rdi			\n\t"\
		"addq	%%rdi,%%rax			\n\t"\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	%%rdi,%%rcx			\n\t"\
		"addq	%%rdi,%%rdx			\n\t"\
		"\n\t"\
		"/* SSE2_RADIX4_DIT_0TWIDDLE(add0, add1, add2, add3, r30): */\n\t"\
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
		"/* SSE2_RADIX4_DIT_0TWIDDLE_2NDOFTWO(add4, add5, add6, add7, r38): */\n\t"\
		"addq	$0x80,%%rsi			\n\t"\
		"movslq	%[__p04],%%rdi		\n\t"\
		"shlq	$3,%%rdi			\n\t"\
		"addq	%%rdi,%%rax			\n\t"\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	%%rdi,%%rcx			\n\t"\
		"addq	%%rdi,%%rdx			\n\t"\
		"\n\t"\
		"movaps	    (%%rax),%%xmm2	\n\t"\
		"movaps	    (%%rcx),%%xmm6	\n\t"\
		"movaps	0x10(%%rax),%%xmm3	\n\t"\
		"movaps	0x10(%%rcx),%%xmm7	\n\t"\
		"movaps	    (%%rbx),%%xmm0	\n\t"\
		"movaps	    (%%rdx),%%xmm4	\n\t"\
		"movaps	0x10(%%rbx),%%xmm1	\n\t"\
		"movaps	0x10(%%rdx),%%xmm5	\n\t"\
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
		"movaps	(%%rdi),%%xmm5	/* ISRT2 */\n\t"\
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
		"/* SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS(r30,r32,r34,r36,r38,r3A,r3C,r3E): */\n\t"\
		"movaps		-0x080(%%rsi),%%xmm0	/* r30 */\n\t"\
		"movaps		-0x040(%%rsi),%%xmm4	/* r34 */\n\t"\
		"movaps		-0x070(%%rsi),%%xmm1	/* r31 */\n\t"\
		"movaps		-0x030(%%rsi),%%xmm5	/* r35 */\n\t"\
		"movaps		      (%%rsi),%%xmm2	/* r38 */\n\t"\
		"movaps		 0x050(%%rsi),%%xmm7	/* r3D */\n\t"\
		"movaps		 0x010(%%rsi),%%xmm3	/* r39 */\n\t"\
		"movaps		 0x040(%%rsi),%%xmm6	/* r3C */\n\t"\
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
		"movaps		-0x060(%%rsi),%%xmm0	/* r32 */\n\t"\
		"movaps		-0x020(%%rsi),%%xmm4	/* r36 */\n\t"\
		"movaps		-0x050(%%rsi),%%xmm1	/* r33 */\n\t"\
		"movaps		-0x010(%%rsi),%%xmm5	/* r37 */\n\t"\
		"movaps		 0x020(%%rsi),%%xmm2	/* r3A */\n\t"\
		"movaps		 0x070(%%rsi),%%xmm7	/* r3F */\n\t"\
		"movaps		 0x030(%%rsi),%%xmm3	/* r3B */\n\t"\
		"movaps		 0x060(%%rsi),%%xmm6	/* r3E */\n\t"\
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
	"/*...and now do eight radix-4 transforms, including the internal twiddle factors:	*/\n\t"\
		"/*...Block 1: r00,r10,r20,r30	*/\n\t"\
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
		"movaps	    (%%rax),%%xmm0		\n\t"\
		"movaps	0x10(%%rax),%%xmm1		\n\t"\
		"movaps	    (%%rbx),%%xmm2		\n\t"\
		"movaps	0x10(%%rbx),%%xmm3		\n\t"\
		"\n\t"\
		"subpd	    (%%rbx),%%xmm0		\n\t"\
		"subpd	0x10(%%rbx),%%xmm1		\n\t"\
		"addpd	    (%%rax),%%xmm2		\n\t"\
		"addpd	0x10(%%rax),%%xmm3		\n\t"\
		"\n\t"\
		"movaps	    (%%rcx),%%xmm4		\n\t"\
		"movaps	0x10(%%rcx),%%xmm5		\n\t"\
		"movaps	    (%%rdx),%%xmm6		\n\t"\
		"movaps	0x10(%%rdx),%%xmm7		\n\t"\
		"\n\t"\
		"subpd	    (%%rdx),%%xmm4		\n\t"\
		"subpd	0x10(%%rdx),%%xmm5		\n\t"\
		"addpd	    (%%rcx),%%xmm6		\n\t"\
		"addpd	0x10(%%rcx),%%xmm7		\n\t"\
		"\n\t"\
		"/* SSE2_RADIX4_DIT_LASTPAIR_IN_PLACE_STRIDE(xmm2,xmm3,xmm0,xmm1,xmm6,xmm7,xmm4,xmm5): swap xmm[01]<->[23], xmm[45]<->[67] */\n\t"\
		"addpd	%%xmm6,%%xmm2		\n\t"\
		"addpd	%%xmm7,%%xmm3		\n\t"\
		"movaps	%%xmm2,    (%%rax)	\n\t"\
		"movaps	%%xmm3,0x10(%%rax)	\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"subpd	%%xmm6,%%xmm2		\n\t"\
		"subpd	%%xmm7,%%xmm3		\n\t"\
		"movaps	%%xmm2,    (%%rcx)	\n\t"\
		"movaps	%%xmm3,0x10(%%rcx)	\n\t"\
		"\n\t"\
		"addpd	%%xmm5,%%xmm0		\n\t"\
		"subpd	%%xmm4,%%xmm1		\n\t"\
		"movaps	%%xmm0,    (%%rbx)	\n\t"\
		"movaps	%%xmm1,0x10(%%rbx)	\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"subpd	%%xmm5,%%xmm0		\n\t"\
		"addpd	%%xmm4,%%xmm1		\n\t"\
		"movaps	%%xmm0,    (%%rdx)	\n\t"\
		"movaps	%%xmm1,0x10(%%rdx)	\n\t"\
		"\n\t"\
		"/*...Block 5: r08,r18,r28,r38	*/\n\t"\
		"\n\t"\
		"addq	$0x80,%%rax		/* r08 */\n\t"\
		"addq	$0x80,%%rbx		\n\t"\
		"addq	$0x80,%%rcx		\n\t"\
		"addq	$0x80,%%rdx		\n\t"\
		"movaps	(%%rdi),%%xmm2	/* isrt2 */\n\t"\
		"\n\t"\
		"movaps	    (%%rcx),%%xmm4	\n\t"\
		"movaps	0x10(%%rcx),%%xmm5	\n\t"\
		"movaps	    (%%rdx),%%xmm0	\n\t"\
		"movaps	0x10(%%rdx),%%xmm1	\n\t"\
		"\n\t"\
		"addpd	0x10(%%rcx),%%xmm4	\n\t"\
		"subpd	    (%%rcx),%%xmm5	\n\t"\
		"subpd	0x10(%%rdx),%%xmm0	\n\t"\
		"addpd	    (%%rdx),%%xmm1	\n\t"\
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
		"movaps	    (%%rax),%%xmm0	\n\t"\
		"movaps	0x10(%%rax),%%xmm1	\n\t"\
		"movaps	    (%%rbx),%%xmm2	\n\t"\
		"movaps	0x10(%%rbx),%%xmm3	\n\t"\
		"\n\t"\
		"subpd	0x10(%%rbx),%%xmm0	\n\t"\
		"subpd	    (%%rbx),%%xmm1	\n\t"\
		"addpd	    (%%rax),%%xmm3	\n\t"\
		"addpd	0x10(%%rax),%%xmm2	\n\t"\
		"\n\t"\
		"/* SSE2_RADIX4_DIT_LASTPAIR_IN_PLACE_STRIDE(xmm3,xmm1,xmm0,xmm2,xmm4,xmm5,xmm6,xmm7): swap xmm0123<->3102 */\n\t"\
		"addpd	%%xmm4,%%xmm3		\n\t"\
		"addpd	%%xmm5,%%xmm1		\n\t"\
		"movaps	%%xmm3,    (%%rax)	\n\t"\
		"movaps	%%xmm1,0x10(%%rax)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"subpd	%%xmm4,%%xmm3		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"movaps	%%xmm3,    (%%rcx)	\n\t"\
		"movaps	%%xmm1,0x10(%%rcx)	\n\t"\
		"\n\t"\
		"addpd	%%xmm7,%%xmm0		\n\t"\
		"subpd	%%xmm6,%%xmm2		\n\t"\
		"movaps	%%xmm0,    (%%rbx)	\n\t"\
		"movaps	%%xmm2,0x10(%%rbx)	\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"subpd	%%xmm7,%%xmm0		\n\t"\
		"addpd	%%xmm6,%%xmm2		\n\t"\
		"movaps	%%xmm0,    (%%rdx)	\n\t"\
		"movaps	%%xmm2,0x10(%%rdx)	\n\t"\
		"\n\t"\
		"/*...Block 3: r04,r14,r24,r34	*/\n\t"\
		"\n\t"\
		"subq	$0x40,%%rax 	/* r04 */\n\t"\
		"subq	$0x40,%%rbx		\n\t"\
		"subq	$0x40,%%rcx		\n\t"\
		"subq	$0x40,%%rdx		\n\t"\
		"\n\t"\
		"movaps	    (%%rcx),%%xmm4	\n\t"\
		"movaps	    (%%rdx),%%xmm0	\n\t"\
		"movaps	0x10(%%rcx),%%xmm5	\n\t"\
		"movaps	0x10(%%rdx),%%xmm1	\n\t"\
		"movaps	    (%%rcx),%%xmm6	\n\t"\
		"movaps	    (%%rdx),%%xmm2	\n\t"\
		"movaps	0x10(%%rcx),%%xmm7	\n\t"\
		"movaps	0x10(%%rdx),%%xmm3	\n\t"\
		"\n\t"\
		"mulpd	    (%%rsi),%%xmm4	\n\t"\
		"mulpd	0x10(%%rsi),%%xmm0	\n\t"\
		"mulpd	    (%%rsi),%%xmm5	\n\t"\
		"mulpd	0x10(%%rsi),%%xmm1	\n\t"\
		"mulpd	0x10(%%rsi),%%xmm6	\n\t"\
		"mulpd	    (%%rsi),%%xmm2	\n\t"\
		"mulpd	0x10(%%rsi),%%xmm7	\n\t"\
		"mulpd	    (%%rsi),%%xmm3	\n\t"\
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
		"movaps	    (%%rbx),%%xmm2	\n\t"\
		"movaps	0x10(%%rbx),%%xmm3	\n\t"\
		"movaps	    (%%rax),%%xmm0	\n\t"\
		"movaps	0x10(%%rax),%%xmm1	\n\t"\
		"addpd	0x10(%%rbx),%%xmm2	\n\t"\
		"subpd	    (%%rbx),%%xmm3	\n\t"\
		"mulpd	(%%rdi),%%xmm2		\n\t"\
		"mulpd	(%%rdi),%%xmm3		\n\t"\
		"\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"\n\t"\
		"/* SSE2_RADIX4_DIT_LASTPAIR_IN_PLACE_STRIDE(xmm2,xmm3,xmm0,xmm1,xmm4,xmm5,xmm6,xmm7): swap xmm[01]<->[23] */\n\t"\
		"addpd	%%xmm4,%%xmm2		\n\t"\
		"addpd	%%xmm5,%%xmm3		\n\t"\
		"movaps	%%xmm2,    (%%rax)	\n\t"\
		"movaps	%%xmm3,0x10(%%rax)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"subpd	%%xmm4,%%xmm2		\n\t"\
		"subpd	%%xmm5,%%xmm3		\n\t"\
		"movaps	%%xmm2,    (%%rcx)	\n\t"\
		"movaps	%%xmm3,0x10(%%rcx)	\n\t"\
		"\n\t"\
		"addpd	%%xmm7,%%xmm0		\n\t"\
		"subpd	%%xmm6,%%xmm1		\n\t"\
		"movaps	%%xmm0,    (%%rbx)	\n\t"\
		"movaps	%%xmm1,0x10(%%rbx)	\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"subpd	%%xmm7,%%xmm0		\n\t"\
		"addpd	%%xmm6,%%xmm1		\n\t"\
		"movaps	%%xmm0,    (%%rdx)	\n\t"\
		"movaps	%%xmm1,0x10(%%rdx)	\n\t"\
		"\n\t"\
		"/*...Block 7: r0C,r1C,r2C,r3C	*/\n\t"\
		"\n\t"\
		"addq	$0x80,%%rax 	/* r0C */\n\t"\
		"addq	$0x80,%%rbx		\n\t"\
		"addq	$0x80,%%rcx		\n\t"\
		"addq	$0x80,%%rdx		\n\t"\
		"\n\t"\
		"movaps	    (%%rcx),%%xmm4	\n\t"\
		"movaps	    (%%rdx),%%xmm0	\n\t"\
		"movaps	0x10(%%rcx),%%xmm5	\n\t"\
		"movaps	0x10(%%rdx),%%xmm1	\n\t"\
		"movaps	    (%%rcx),%%xmm6	\n\t"\
		"movaps	    (%%rdx),%%xmm2	\n\t"\
		"movaps	0x10(%%rcx),%%xmm7	\n\t"\
		"movaps	0x10(%%rdx),%%xmm3	\n\t"\
		"\n\t"\
		"mulpd	0x10(%%rsi),%%xmm4	\n\t"\
		"mulpd	    (%%rsi),%%xmm0	\n\t"\
		"mulpd	0x10(%%rsi),%%xmm5	\n\t"\
		"mulpd	    (%%rsi),%%xmm1	\n\t"\
		"mulpd	    (%%rsi),%%xmm6	\n\t"\
		"mulpd	0x10(%%rsi),%%xmm2	\n\t"\
		"mulpd	    (%%rsi),%%xmm7	\n\t"\
		"mulpd	0x10(%%rsi),%%xmm3	\n\t"\
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
		"movaps	    (%%rbx),%%xmm2	\n\t"\
		"movaps	0x10(%%rbx),%%xmm3	\n\t"\
		"movaps	    (%%rax),%%xmm0	\n\t"\
		"movaps	0x10(%%rax),%%xmm1	\n\t"\
		"subpd	0x10(%%rbx),%%xmm2	\n\t"\
		"addpd	    (%%rbx),%%xmm3	\n\t"\
		"mulpd	(%%rdi),%%xmm2		\n\t"\
		"mulpd	(%%rdi),%%xmm3		\n\t"\
		"\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"\n\t"\
		"/* SSE2_RADIX4_DIT_LASTPAIR_IN_PLACE_STRIDE(xmm0,xmm1,xmm2,xmm3,xmm6,xmm7,xmm4,xmm5): swap xmm[45]<->[67] */\n\t"\
		"addpd	%%xmm6,%%xmm0		\n\t"\
		"addpd	%%xmm7,%%xmm1		\n\t"\
		"movaps	%%xmm0,    (%%rax)	\n\t"\
		"movaps	%%xmm1,0x10(%%rax)	\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"subpd	%%xmm6,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm1		\n\t"\
		"movaps	%%xmm0,    (%%rcx)	\n\t"\
		"movaps	%%xmm1,0x10(%%rcx)	\n\t"\
		"\n\t"\
		"addpd	%%xmm5,%%xmm2		\n\t"\
		"subpd	%%xmm4,%%xmm3		\n\t"\
		"movaps	%%xmm2,    (%%rbx)	\n\t"\
		"movaps	%%xmm3,0x10(%%rbx)	\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"subpd	%%xmm5,%%xmm2		\n\t"\
		"addpd	%%xmm4,%%xmm3		\n\t"\
		"movaps	%%xmm2,    (%%rdx)	\n\t"\
		"movaps	%%xmm3,0x10(%%rdx)	\n\t"\
		"\n\t"\
		"/*...Block 2: r02,r12,r22,r32	*/\n\t"\
		"\n\t"\
		"subq	$0xa0,%%rax 	/* r02 */\n\t"\
		"subq	$0xa0,%%rbx			\n\t"\
		"subq	$0xa0,%%rcx			\n\t"\
		"subq	$0xa0,%%rdx			\n\t"\
		"addq	$0x30,%%rdi /* cc1 */\n\t"\
		"addq	$0x40,%%rsi /* cc3 */\n\t"\
		"\n\t"\
		"movaps	    (%%rcx),%%xmm4	\n\t"\
		"movaps	    (%%rdx),%%xmm0	\n\t"\
		"movaps	0x10(%%rcx),%%xmm5	\n\t"\
		"movaps	0x10(%%rdx),%%xmm1	\n\t"\
		"movaps	    (%%rcx),%%xmm6	\n\t"\
		"movaps	    (%%rdx),%%xmm2	\n\t"\
		"movaps	0x10(%%rcx),%%xmm7	\n\t"\
		"movaps	0x10(%%rdx),%%xmm3	\n\t"\
		"\n\t"\
		"mulpd	    (%%rdi),%%xmm4	\n\t"\
		"mulpd	    (%%rsi),%%xmm0	\n\t"\
		"mulpd	    (%%rdi),%%xmm5	\n\t"\
		"mulpd	    (%%rsi),%%xmm1	\n\t"\
		"mulpd	0x10(%%rdi),%%xmm6	\n\t"\
		"mulpd	0x10(%%rsi),%%xmm2	\n\t"\
		"mulpd	0x10(%%rdi),%%xmm7	\n\t"\
		"mulpd	0x10(%%rsi),%%xmm3	\n\t"\
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
		"subq	$0x40,%%rsi /* cc0 */\n\t"\
		"movaps	    (%%rbx),%%xmm2	\n\t"\
		"movaps	0x10(%%rbx),%%xmm3	\n\t"\
		"movaps	    (%%rbx),%%xmm0	\n\t"\
		"movaps	0x10(%%rbx),%%xmm1	\n\t"\
		"\n\t"\
		"mulpd	    (%%rsi),%%xmm2	\n\t"\
		"mulpd	0x10(%%rsi),%%xmm1	\n\t"\
		"mulpd	    (%%rsi),%%xmm3	\n\t"\
		"mulpd	0x10(%%rsi),%%xmm0	\n\t"\
		"addpd	%%xmm1,%%xmm2		\n\t"\
		"subpd	%%xmm0,%%xmm3		\n\t"\
		"\n\t"\
		"movaps	    (%%rax),%%xmm0	\n\t"\
		"movaps	0x10(%%rax),%%xmm1	\n\t"\
		"\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"\n\t"\
		"/* SSE2_RADIX4_DIT_LASTPAIR_IN_PLACE_STRIDE(xmm2,xmm3,xmm0,xmm1,xmm6,xmm7,xmm4,xmm5): swap xmm[01]<->[23], xmm[45]<->[67] */\n\t"\
		"addpd	%%xmm6,%%xmm2		\n\t"\
		"addpd	%%xmm7,%%xmm3		\n\t"\
		"movaps	%%xmm2,    (%%rax)	\n\t"\
		"movaps	%%xmm3,0x10(%%rax)	\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"subpd	%%xmm6,%%xmm2		\n\t"\
		"subpd	%%xmm7,%%xmm3		\n\t"\
		"movaps	%%xmm2,    (%%rcx)	\n\t"\
		"movaps	%%xmm3,0x10(%%rcx)	\n\t"\
		"\n\t"\
		"addpd	%%xmm5,%%xmm0		\n\t"\
		"subpd	%%xmm4,%%xmm1		\n\t"\
		"movaps	%%xmm0,    (%%rbx)	\n\t"\
		"movaps	%%xmm1,0x10(%%rbx)	\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"subpd	%%xmm5,%%xmm0		\n\t"\
		"addpd	%%xmm4,%%xmm1		\n\t"\
		"movaps	%%xmm0,    (%%rdx)	\n\t"\
		"movaps	%%xmm1,0x10(%%rdx)	\n\t"\
		"\n\t"\
		"/*...Block 6: r0A,r1A,r2A,r3A	*/\n\t"\
		"\n\t"\
		"addq	$0x80,%%rax 	/* r0A */\n\t"\
		"addq	$0x80,%%rbx			\n\t"\
		"addq	$0x80,%%rcx			\n\t"\
		"addq	$0x80,%%rdx			\n\t"\
		"addq	$0x40,%%rsi /* cc3 */\n\t"\
		"\n\t"\
		"movaps	    (%%rcx),%%xmm4	\n\t"\
		"movaps	    (%%rdx),%%xmm0	\n\t"\
		"movaps	0x10(%%rcx),%%xmm5	\n\t"\
		"movaps	0x10(%%rdx),%%xmm1	\n\t"\
		"movaps	    (%%rcx),%%xmm6	\n\t"\
		"movaps	    (%%rdx),%%xmm2	\n\t"\
		"movaps	0x10(%%rcx),%%xmm7	\n\t"\
		"movaps	0x10(%%rdx),%%xmm3	\n\t"\
		"\n\t"\
		"mulpd	0x10(%%rsi),%%xmm4	\n\t"\
		"mulpd	    (%%rdi),%%xmm0	\n\t"\
		"mulpd	0x10(%%rsi),%%xmm5	\n\t"\
		"mulpd	    (%%rdi),%%xmm1	\n\t"\
		"mulpd	    (%%rsi),%%xmm6	\n\t"\
		"mulpd	0x10(%%rdi),%%xmm2	\n\t"\
		"mulpd	    (%%rsi),%%xmm7	\n\t"\
		"mulpd	0x10(%%rdi),%%xmm3	\n\t"\
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
		"subq	$0x40,%%rsi /* cc0 */\n\t"\
		"movaps	    (%%rbx),%%xmm2	\n\t"\
		"movaps	0x10(%%rbx),%%xmm3	\n\t"\
		"movaps	    (%%rbx),%%xmm0	\n\t"\
		"movaps	0x10(%%rbx),%%xmm1	\n\t"\
		"\n\t"\
		"mulpd	0x10(%%rsi),%%xmm2	\n\t"\
		"mulpd	    (%%rsi),%%xmm1	\n\t"\
		"mulpd	0x10(%%rsi),%%xmm3	\n\t"\
		"mulpd	    (%%rsi),%%xmm0	\n\t"\
		"subpd	%%xmm1,%%xmm2		\n\t"\
		"addpd	%%xmm0,%%xmm3		\n\t"\
		"\n\t"\
		"movaps	    (%%rax),%%xmm0	\n\t"\
		"movaps	0x10(%%rax),%%xmm1	\n\t"\
		"\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"\n\t"\
		"/* SSE2_RADIX4_DIT_LASTPAIR_IN_PLACE_STRIDE(xmm0,xmm1,xmm2,xmm3,xmm6,xmm7,xmm4,xmm5): swap xmm[45]<->[67] */\n\t"\
		"addpd	%%xmm6,%%xmm0		\n\t"\
		"addpd	%%xmm7,%%xmm1		\n\t"\
		"movaps	%%xmm0,    (%%rax)	\n\t"\
		"movaps	%%xmm1,0x10(%%rax)	\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"subpd	%%xmm6,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm1		\n\t"\
		"movaps	%%xmm0,    (%%rcx)	\n\t"\
		"movaps	%%xmm1,0x10(%%rcx)	\n\t"\
		"\n\t"\
		"addpd	%%xmm5,%%xmm2		\n\t"\
		"subpd	%%xmm4,%%xmm3		\n\t"\
		"movaps	%%xmm2,    (%%rbx)	\n\t"\
		"movaps	%%xmm3,0x10(%%rbx)	\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"subpd	%%xmm5,%%xmm2		\n\t"\
		"addpd	%%xmm4,%%xmm3		\n\t"\
		"movaps	%%xmm2,    (%%rdx)	\n\t"\
		"movaps	%%xmm3,0x10(%%rdx)	\n\t"\
		"\n\t"\
		"/*...Block 4: r06,r16,r26,r36	*/\n\t"\
		"\n\t"\
		"subq	$0x40,%%rax 	/* r06 */\n\t"\
		"subq	$0x40,%%rbx			\n\t"\
		"subq	$0x40,%%rcx			\n\t"\
		"subq	$0x40,%%rdx			\n\t"\
		"addq	$0x40,%%rsi /* cc3 */\n\t"\
		"\n\t"\
		"movaps	    (%%rcx),%%xmm4	\n\t"\
		"movaps	    (%%rdx),%%xmm0	\n\t"\
		"movaps	0x10(%%rcx),%%xmm5	\n\t"\
		"movaps	0x10(%%rdx),%%xmm1	\n\t"\
		"movaps	    (%%rcx),%%xmm6	\n\t"\
		"movaps	    (%%rdx),%%xmm2	\n\t"\
		"movaps	0x10(%%rcx),%%xmm7	\n\t"\
		"movaps	0x10(%%rdx),%%xmm3	\n\t"\
		"\n\t"\
		"mulpd	    (%%rsi),%%xmm4	\n\t"\
		"mulpd	0x10(%%rdi),%%xmm0	\n\t"\
		"mulpd	    (%%rsi),%%xmm5	\n\t"\
		"mulpd	0x10(%%rdi),%%xmm1	\n\t"\
		"mulpd	0x10(%%rsi),%%xmm6	\n\t"\
		"mulpd	    (%%rdi),%%xmm2	\n\t"\
		"mulpd	0x10(%%rsi),%%xmm7	\n\t"\
		"mulpd	    (%%rdi),%%xmm3	\n\t"\
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
		"subq	$0x40,%%rsi /* cc0 */\n\t"\
		"movaps	    (%%rbx),%%xmm2	\n\t"\
		"movaps	0x10(%%rbx),%%xmm3	\n\t"\
		"movaps	    (%%rbx),%%xmm0	\n\t"\
		"movaps	0x10(%%rbx),%%xmm1	\n\t"\
		"\n\t"\
		"mulpd	0x10(%%rsi),%%xmm2	\n\t"\
		"mulpd	    (%%rsi),%%xmm1	\n\t"\
		"mulpd	0x10(%%rsi),%%xmm3	\n\t"\
		"mulpd	    (%%rsi),%%xmm0	\n\t"\
		"addpd	%%xmm1,%%xmm2		\n\t"\
		"subpd	%%xmm0,%%xmm3		\n\t"\
		"\n\t"\
		"movaps	    (%%rax),%%xmm0	\n\t"\
		"movaps	0x10(%%rax),%%xmm1	\n\t"\
		"\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"\n\t"\
		"/* SSE2_RADIX4_DIT_LASTPAIR_IN_PLACE_STRIDE(xmm2,xmm3,xmm0,xmm1,xmm6,xmm7,xmm4,xmm5): swap xmm[01]<->[23], xmm[45]<->[67] */\n\t"\
		"addpd	%%xmm6,%%xmm2		\n\t"\
		"addpd	%%xmm7,%%xmm3		\n\t"\
		"movaps	%%xmm2,    (%%rax)	\n\t"\
		"movaps	%%xmm3,0x10(%%rax)	\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"subpd	%%xmm6,%%xmm2		\n\t"\
		"subpd	%%xmm7,%%xmm3		\n\t"\
		"movaps	%%xmm2,    (%%rcx)	\n\t"\
		"movaps	%%xmm3,0x10(%%rcx)	\n\t"\
		"\n\t"\
		"addpd	%%xmm5,%%xmm0		\n\t"\
		"subpd	%%xmm4,%%xmm1		\n\t"\
		"movaps	%%xmm0,    (%%rbx)	\n\t"\
		"movaps	%%xmm1,0x10(%%rbx)	\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"subpd	%%xmm5,%%xmm0		\n\t"\
		"addpd	%%xmm4,%%xmm1		\n\t"\
		"movaps	%%xmm0,    (%%rdx)	\n\t"\
		"movaps	%%xmm1,0x10(%%rdx)	\n\t"\
		"\n\t"\
		"/*...Block 8: r0E,r1E,r2E,r3E	*/\n\t"\
		"\n\t"\
		"addq	$0x80,%%rax 	/* r0E */\n\t"\
		"addq	$0x80,%%rbx		\n\t"\
		"addq	$0x80,%%rcx		\n\t"\
		"addq	$0x80,%%rdx		\n\t"\
		"addq	$0x40,%%rsi /* cc3 */\n\t"\
		"\n\t"\
		"movaps	    (%%rcx),%%xmm4	\n\t"\
		"movaps	    (%%rdx),%%xmm0	\n\t"\
		"movaps	0x10(%%rcx),%%xmm5	\n\t"\
		"movaps	0x10(%%rdx),%%xmm1	\n\t"\
		"movaps	    (%%rcx),%%xmm6	\n\t"\
		"movaps	    (%%rdx),%%xmm2	\n\t"\
		"movaps	0x10(%%rcx),%%xmm7	\n\t"\
		"movaps	0x10(%%rdx),%%xmm3	\n\t"\
		"\n\t"\
		"mulpd	0x10(%%rdi),%%xmm4	\n\t"\
		"mulpd	0x10(%%rsi),%%xmm0	\n\t"\
		"mulpd	0x10(%%rdi),%%xmm5	\n\t"\
		"mulpd	0x10(%%rsi),%%xmm1	\n\t"\
		"mulpd	    (%%rdi),%%xmm6	\n\t"\
		"mulpd	    (%%rsi),%%xmm2	\n\t"\
		"mulpd	    (%%rdi),%%xmm7	\n\t"\
		"mulpd	    (%%rsi),%%xmm3	\n\t"\
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
		"subq	$0x40,%%rsi /* cc0 */\n\t"\
		"movaps	    (%%rbx),%%xmm2	\n\t"\
		"movaps	0x10(%%rbx),%%xmm3	\n\t"\
		"movaps	    (%%rbx),%%xmm0	\n\t"\
		"movaps	0x10(%%rbx),%%xmm1	\n\t"\
		"\n\t"\
		"mulpd	    (%%rsi),%%xmm2	\n\t"\
		"mulpd	0x10(%%rsi),%%xmm1	\n\t"\
		"mulpd	    (%%rsi),%%xmm3	\n\t"\
		"mulpd	0x10(%%rsi),%%xmm0	\n\t"\
		"subpd	%%xmm1,%%xmm2		\n\t"\
		"addpd	%%xmm0,%%xmm3		\n\t"\
		"\n\t"\
		"movaps	    (%%rax),%%xmm0	\n\t"\
		"movaps	0x10(%%rax),%%xmm1	\n\t"\
		"\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"\n\t"\
		"/* SSE2_RADIX4_DIT_LASTPAIR_IN_PLACE_STRIDE(xmm0,xmm1,xmm2,xmm3,xmm6,xmm7,xmm4,xmm5): swap xmm[45]<->[67] */\n\t"\
		"addpd	%%xmm6,%%xmm0		\n\t"\
		"addpd	%%xmm7,%%xmm1		\n\t"\
		"movaps	%%xmm0,    (%%rax)	\n\t"\
		"movaps	%%xmm1,0x10(%%rax)	\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"subpd	%%xmm6,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm1		\n\t"\
		"movaps	%%xmm0,    (%%rcx)	\n\t"\
		"movaps	%%xmm1,0x10(%%rcx)	\n\t"\
		"\n\t"\
		"addpd	%%xmm5,%%xmm2		\n\t"\
		"subpd	%%xmm4,%%xmm3		\n\t"\
		"movaps	%%xmm2,    (%%rbx)	\n\t"\
		"movaps	%%xmm3,0x10(%%rbx)	\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"subpd	%%xmm5,%%xmm2		\n\t"\
		"addpd	%%xmm4,%%xmm3		\n\t"\
		"movaps	%%xmm2,    (%%rdx)	\n\t"\
		"movaps	%%xmm3,0x10(%%rdx)	\n\t"\
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
		"/* SSE2_RADIX4_DIF_IN_PLACE(r00,r20,r10,r30): */\n\t"\
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
		"movaps	    (%%rax),%%xmm0	\n\t"\
		"movaps	0x10(%%rax),%%xmm1	\n\t"\
		"movaps	    (%%rax),%%xmm2	\n\t"\
		"movaps	0x10(%%rax),%%xmm3	\n\t"\
		"\n\t"\
		"addpd	    (%%rbx),%%xmm0	\n\t"\
		"addpd	0x10(%%rbx),%%xmm1	\n\t"\
		"subpd	    (%%rbx),%%xmm2	\n\t"\
		"subpd	0x10(%%rbx),%%xmm3	\n\t"\
		"\n\t"\
		"movaps	    (%%rcx),%%xmm4	\n\t"\
		"movaps	0x10(%%rcx),%%xmm5	\n\t"\
		"movaps	    (%%rcx),%%xmm6	\n\t"\
		"movaps	0x10(%%rcx),%%xmm7	\n\t"\
		"\n\t"\
		"addpd	    (%%rdx),%%xmm4	\n\t"\
		"addpd	0x10(%%rdx),%%xmm5	\n\t"\
		"subpd	    (%%rdx),%%xmm6	\n\t"\
		"subpd	0x10(%%rdx),%%xmm7	\n\t"\
		"/* Finish radix-4 butterfly and store results into temporary-array slots: */\n\t"\
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
		"/* SSE2_RADIX4_DIF_IN_PLACE_2NDOFTWO(r08,r28,r18,r38): */\n\t"\
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
		"/* Finish radix-4 butterfly and store results into temporary-array slots: */\n\t"\
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
		"movaps	(%%rdi),%%xmm0	/* isrt2 */\n\t"\
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
		"/* SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS(r00,r10,r20,r30,r08,r18,r28,r38): */\n\t"\
		"subq	$0x080,%%rax	/* r00 */\n\t"\
		"subq	$0x200,%%rbx	/* r08 */\n\t"\
		"addq	$0x080,%%rcx	/* r20 */\n\t"\
		"subq	$0x100,%%rdx	/* r28 */\n\t"\
		"movaps	    (%%rax),%%xmm0	\n\t"\
		"movaps	    (%%rcx),%%xmm4	\n\t"\
		"movaps	0x10(%%rax),%%xmm1	\n\t"\
		"movaps	0x10(%%rcx),%%xmm5	\n\t"\
		"movaps	    (%%rbx),%%xmm2	\n\t"\
		"movaps	0x10(%%rdx),%%xmm7	\n\t"\
		"movaps	0x10(%%rbx),%%xmm3	\n\t"\
		"movaps	    (%%rdx),%%xmm6	\n\t"\
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
		"movaps	%%xmm0,    (%%rbx)	\n\t"\
		"movaps	%%xmm4,    (%%rcx)	\n\t"\
		"movaps	%%xmm1,0x10(%%rbx)	\n\t"\
		"movaps	%%xmm5,0x10(%%rdx)	\n\t"\
		"movaps	%%xmm2,    (%%rax)	\n\t"\
		"movaps	%%xmm7,    (%%rdx)	\n\t"\
		"movaps	%%xmm3,0x10(%%rax)	\n\t"\
		"movaps	%%xmm6,0x10(%%rcx)	\n\t"\
		"\n\t"\
		"addq	$0x100,%%rax	/* r10 */	\n\t"\
		"addq	$0x100,%%rbx	/* r18 */	\n\t"\
		"addq	$0x100,%%rcx	/* r30 */	\n\t"\
		"addq	$0x100,%%rdx	/* r38 */	\n\t"\
		"\n\t"\
		"movaps	    (%%rax),%%xmm0	\n\t"\
		"movaps	    (%%rcx),%%xmm4	\n\t"\
		"movaps	0x10(%%rax),%%xmm1	\n\t"\
		"movaps	0x10(%%rcx),%%xmm5	\n\t"\
		"movaps	    (%%rbx),%%xmm2	\n\t"\
		"movaps	0x10(%%rdx),%%xmm7	\n\t"\
		"movaps	0x10(%%rbx),%%xmm3	\n\t"\
		"movaps	    (%%rdx),%%xmm6	\n\t"\
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
		"movaps	%%xmm0,    (%%rbx)	\n\t"\
		"movaps	%%xmm4,    (%%rcx)	\n\t"\
		"movaps	%%xmm1,0x10(%%rbx)	\n\t"\
		"movaps	%%xmm5,0x10(%%rdx)	\n\t"\
		"movaps	%%xmm2,    (%%rax)	\n\t"\
		"movaps	%%xmm7,    (%%rdx)	\n\t"\
		"movaps	%%xmm3,0x10(%%rax)	\n\t"\
		"movaps	%%xmm6,0x10(%%rcx)	\n\t"\
		"\n\t"\
		"/* SSE2_RADIX4_DIF_IN_PLACE(r04,r24,r14,r34): */\n\t"\
		"\n\t"\
		"subq	$0x0C0,%%rax	/* r04 */	\n\t"\
		"addq	$0x0C0,%%rbx	/* r24 */	\n\t"\
		"subq	$0x1C0,%%rcx	/* r14 */	\n\t"\
		"subq	$0x040,%%rdx	/* r34 */	\n\t"\
		"\n\t"\
		"movaps	    (%%rax),%%xmm0	\n\t"\
		"movaps	0x10(%%rax),%%xmm1	\n\t"\
		"movaps	    (%%rax),%%xmm2	\n\t"\
		"movaps	0x10(%%rax),%%xmm3	\n\t"\
		"\n\t"\
		"addpd	    (%%rbx),%%xmm0	\n\t"\
		"addpd	0x10(%%rbx),%%xmm1	\n\t"\
		"subpd	    (%%rbx),%%xmm2	\n\t"\
		"subpd	0x10(%%rbx),%%xmm3	\n\t"\
		"\n\t"\
		"movaps	    (%%rcx),%%xmm4	\n\t"\
		"movaps	0x10(%%rcx),%%xmm5	\n\t"\
		"movaps	    (%%rcx),%%xmm6	\n\t"\
		"movaps	0x10(%%rcx),%%xmm7	\n\t"\
		"\n\t"\
		"addpd	    (%%rdx),%%xmm4	\n\t"\
		"addpd	0x10(%%rdx),%%xmm5	\n\t"\
		"subpd	    (%%rdx),%%xmm6	\n\t"\
		"subpd	0x10(%%rdx),%%xmm7	\n\t"\
		"/* Finish radix-4 butterfly and store results into temporary-array slots: */\n\t"\
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
		"/* SSE2_RADIX4_DIF_IN_PLACE_2NDOFTWO(r0C,r2C,r1C,r3C): */\n\t"\
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
		"/* Finish radix-4 butterfly and store results into temporary-array slots: */\n\t"\
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
		"movaps	(%%rdi),%%xmm0	/* isrt2 */\n\t"\
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
		"/* SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS(r04,r14,r24,r34,r0C,r1C,r2C,r3C): */\n\t"\
		"subq	$0x080,%%rax	/* r04 */\n\t"\
		"subq	$0x200,%%rbx	/* r0C */\n\t"\
		"addq	$0x080,%%rcx	/* r24 */\n\t"\
		"subq	$0x100,%%rdx	/* r2C */\n\t"\
		"movaps	    (%%rax),%%xmm0	\n\t"\
		"movaps	    (%%rcx),%%xmm4	\n\t"\
		"movaps	0x10(%%rax),%%xmm1	\n\t"\
		"movaps	0x10(%%rcx),%%xmm5	\n\t"\
		"movaps	    (%%rbx),%%xmm2	\n\t"\
		"movaps	0x10(%%rdx),%%xmm7	\n\t"\
		"movaps	0x10(%%rbx),%%xmm3	\n\t"\
		"movaps	    (%%rdx),%%xmm6	\n\t"\
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
		"movaps	%%xmm0,    (%%rbx)	\n\t"\
		"movaps	%%xmm4,    (%%rcx)	\n\t"\
		"movaps	%%xmm1,0x10(%%rbx)	\n\t"\
		"movaps	%%xmm5,0x10(%%rdx)	\n\t"\
		"movaps	%%xmm2,    (%%rax)	\n\t"\
		"movaps	%%xmm7,    (%%rdx)	\n\t"\
		"movaps	%%xmm3,0x10(%%rax)	\n\t"\
		"movaps	%%xmm6,0x10(%%rcx)	\n\t"\
		"\n\t"\
		"addq	$0x100,%%rax	/* r14 */	\n\t"\
		"addq	$0x100,%%rbx	/* r1C */	\n\t"\
		"addq	$0x100,%%rcx	/* r34 */	\n\t"\
		"addq	$0x100,%%rdx	/* r3C */	\n\t"\
		"\n\t"\
		"movaps	    (%%rax),%%xmm0	\n\t"\
		"movaps	    (%%rcx),%%xmm4	\n\t"\
		"movaps	0x10(%%rax),%%xmm1	\n\t"\
		"movaps	0x10(%%rcx),%%xmm5	\n\t"\
		"movaps	    (%%rbx),%%xmm2	\n\t"\
		"movaps	0x10(%%rdx),%%xmm7	\n\t"\
		"movaps	0x10(%%rbx),%%xmm3	\n\t"\
		"movaps	    (%%rdx),%%xmm6	\n\t"\
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
		"movaps	%%xmm0,    (%%rbx)	\n\t"\
		"movaps	%%xmm4,    (%%rcx)	\n\t"\
		"movaps	%%xmm1,0x10(%%rbx)	\n\t"\
		"movaps	%%xmm5,0x10(%%rdx)	\n\t"\
		"movaps	%%xmm2,    (%%rax)	\n\t"\
		"movaps	%%xmm7,    (%%rdx)	\n\t"\
		"movaps	%%xmm3,0x10(%%rax)	\n\t"\
		"movaps	%%xmm6,0x10(%%rcx)	\n\t"\
		"\n\t"\
		"/* SSE2_RADIX4_DIF_IN_PLACE(r02,r22,r12,r32): */\n\t"\
		"\n\t"\
		"subq	$0x120,%%rax	/* r02 */	\n\t"\
		"addq	$0x060,%%rbx	/* r22 */	\n\t"\
		"subq	$0x220,%%rcx	/* r12 */	\n\t"\
		"subq	$0x0a0,%%rdx	/* r32 */	\n\t"\
		"\n\t"\
		"movaps	    (%%rax),%%xmm0	\n\t"\
		"movaps	0x10(%%rax),%%xmm1	\n\t"\
		"movaps	    (%%rax),%%xmm2	\n\t"\
		"movaps	0x10(%%rax),%%xmm3	\n\t"\
		"\n\t"\
		"addpd	    (%%rbx),%%xmm0	\n\t"\
		"addpd	0x10(%%rbx),%%xmm1	\n\t"\
		"subpd	    (%%rbx),%%xmm2	\n\t"\
		"subpd	0x10(%%rbx),%%xmm3	\n\t"\
		"\n\t"\
		"movaps	    (%%rcx),%%xmm4	\n\t"\
		"movaps	0x10(%%rcx),%%xmm5	\n\t"\
		"movaps	    (%%rcx),%%xmm6	\n\t"\
		"movaps	0x10(%%rcx),%%xmm7	\n\t"\
		"\n\t"\
		"addpd	    (%%rdx),%%xmm4	\n\t"\
		"addpd	0x10(%%rdx),%%xmm5	\n\t"\
		"subpd	    (%%rdx),%%xmm6	\n\t"\
		"subpd	0x10(%%rdx),%%xmm7	\n\t"\
		"/* Finish radix-4 butterfly and store results into temporary-array slots: */\n\t"\
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
		"/* SSE2_RADIX4_DIF_IN_PLACE_2NDOFTWO(r0A,r2A,r1A,r3A): */\n\t"\
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
		"/* Finish radix-4 butterfly and store results into temporary-array slots: */\n\t"\
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
		"movaps	(%%rdi),%%xmm0	/* isrt2 */\n\t"\
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
		"/* SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS(r02,r12,r22,r32,r0A,r1A,r2A,r3A): */\n\t"\
		"subq	$0x080,%%rax	/* r02 */\n\t"\
		"subq	$0x200,%%rbx	/* r0A */\n\t"\
		"addq	$0x080,%%rcx	/* r22 */\n\t"\
		"subq	$0x100,%%rdx	/* r2A */\n\t"\
		"movaps	    (%%rax),%%xmm0	\n\t"\
		"movaps	    (%%rcx),%%xmm4	\n\t"\
		"movaps	0x10(%%rax),%%xmm1	\n\t"\
		"movaps	0x10(%%rcx),%%xmm5	\n\t"\
		"movaps	    (%%rbx),%%xmm2	\n\t"\
		"movaps	0x10(%%rdx),%%xmm7	\n\t"\
		"movaps	0x10(%%rbx),%%xmm3	\n\t"\
		"movaps	    (%%rdx),%%xmm6	\n\t"\
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
		"movaps	%%xmm0,    (%%rbx)	\n\t"\
		"movaps	%%xmm4,    (%%rcx)	\n\t"\
		"movaps	%%xmm1,0x10(%%rbx)	\n\t"\
		"movaps	%%xmm5,0x10(%%rdx)	\n\t"\
		"movaps	%%xmm2,    (%%rax)	\n\t"\
		"movaps	%%xmm7,    (%%rdx)	\n\t"\
		"movaps	%%xmm3,0x10(%%rax)	\n\t"\
		"movaps	%%xmm6,0x10(%%rcx)	\n\t"\
		"\n\t"\
		"addq	$0x100,%%rax	/* r12 */	\n\t"\
		"addq	$0x100,%%rbx	/* r1A */	\n\t"\
		"addq	$0x100,%%rcx	/* r32 */	\n\t"\
		"addq	$0x100,%%rdx	/* r3A */	\n\t"\
		"\n\t"\
		"movaps	    (%%rax),%%xmm0	\n\t"\
		"movaps	    (%%rcx),%%xmm4	\n\t"\
		"movaps	0x10(%%rax),%%xmm1	\n\t"\
		"movaps	0x10(%%rcx),%%xmm5	\n\t"\
		"movaps	    (%%rbx),%%xmm2	\n\t"\
		"movaps	0x10(%%rdx),%%xmm7	\n\t"\
		"movaps	0x10(%%rbx),%%xmm3	\n\t"\
		"movaps	    (%%rdx),%%xmm6	\n\t"\
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
		"movaps	%%xmm0,    (%%rbx)	\n\t"\
		"movaps	%%xmm4,    (%%rcx)	\n\t"\
		"movaps	%%xmm1,0x10(%%rbx)	\n\t"\
		"movaps	%%xmm5,0x10(%%rdx)	\n\t"\
		"movaps	%%xmm2,    (%%rax)	\n\t"\
		"movaps	%%xmm7,    (%%rdx)	\n\t"\
		"movaps	%%xmm3,0x10(%%rax)	\n\t"\
		"movaps	%%xmm6,0x10(%%rcx)	\n\t"\
		"\n\t"\
		"/* SSE2_RADIX4_DIF_IN_PLACE(r06,r26,r16,r36): */\n\t"\
		"\n\t"\
		"subq	$0x0C0,%%rax	/* r06 */	\n\t"\
		"addq	$0x0C0,%%rbx	/* r26 */	\n\t"\
		"subq	$0x1C0,%%rcx	/* r16 */	\n\t"\
		"subq	$0x040,%%rdx	/* r36 */	\n\t"\
		"\n\t"\
		"movaps	    (%%rax),%%xmm0	\n\t"\
		"movaps	0x10(%%rax),%%xmm1	\n\t"\
		"movaps	    (%%rax),%%xmm2	\n\t"\
		"movaps	0x10(%%rax),%%xmm3	\n\t"\
		"\n\t"\
		"addpd	    (%%rbx),%%xmm0	\n\t"\
		"addpd	0x10(%%rbx),%%xmm1	\n\t"\
		"subpd	    (%%rbx),%%xmm2	\n\t"\
		"subpd	0x10(%%rbx),%%xmm3	\n\t"\
		"\n\t"\
		"movaps	    (%%rcx),%%xmm4	\n\t"\
		"movaps	0x10(%%rcx),%%xmm5	\n\t"\
		"movaps	    (%%rcx),%%xmm6	\n\t"\
		"movaps	0x10(%%rcx),%%xmm7	\n\t"\
		"\n\t"\
		"addpd	    (%%rdx),%%xmm4	\n\t"\
		"addpd	0x10(%%rdx),%%xmm5	\n\t"\
		"subpd	    (%%rdx),%%xmm6	\n\t"\
		"subpd	0x10(%%rdx),%%xmm7	\n\t"\
		"/* Finish radix-4 butterfly and store results into temporary-array slots: */\n\t"\
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
		"/* SSE2_RADIX4_DIF_IN_PLACE_2NDOFTWO(r0E,r2E,r1E,r3E): */\n\t"\
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
		"/* Finish radix-4 butterfly and store results into temporary-array slots: */\n\t"\
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
		"movaps	(%%rdi),%%xmm0	/* isrt2 */\n\t"\
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
		"/* SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS(r06,r16,r26,r36,r0E,r1E,r2E,r3E): */\n\t"\
		"subq	$0x080,%%rax	/* r02 */\n\t"\
		"subq	$0x200,%%rbx	/* r0A */\n\t"\
		"addq	$0x080,%%rcx	/* r22 */\n\t"\
		"subq	$0x100,%%rdx	/* r2A */\n\t"\
		"movaps	    (%%rax),%%xmm0	\n\t"\
		"movaps	    (%%rcx),%%xmm4	\n\t"\
		"movaps	0x10(%%rax),%%xmm1	\n\t"\
		"movaps	0x10(%%rcx),%%xmm5	\n\t"\
		"movaps	    (%%rbx),%%xmm2	\n\t"\
		"movaps	0x10(%%rdx),%%xmm7	\n\t"\
		"movaps	0x10(%%rbx),%%xmm3	\n\t"\
		"movaps	    (%%rdx),%%xmm6	\n\t"\
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
		"movaps	%%xmm0,    (%%rbx)	\n\t"\
		"movaps	%%xmm4,    (%%rcx)	\n\t"\
		"movaps	%%xmm1,0x10(%%rbx)	\n\t"\
		"movaps	%%xmm5,0x10(%%rdx)	\n\t"\
		"movaps	%%xmm2,    (%%rax)	\n\t"\
		"movaps	%%xmm7,    (%%rdx)	\n\t"\
		"movaps	%%xmm3,0x10(%%rax)	\n\t"\
		"movaps	%%xmm6,0x10(%%rcx)	\n\t"\
		"\n\t"\
		"addq	$0x100,%%rax	/* r12 */	\n\t"\
		"addq	$0x100,%%rbx	/* r1A */	\n\t"\
		"addq	$0x100,%%rcx	/* r32 */	\n\t"\
		"addq	$0x100,%%rdx	/* r3A */	\n\t"\
		"\n\t"\
		"movaps	    (%%rax),%%xmm0	\n\t"\
		"movaps	    (%%rcx),%%xmm4	\n\t"\
		"movaps	0x10(%%rax),%%xmm1	\n\t"\
		"movaps	0x10(%%rcx),%%xmm5	\n\t"\
		"movaps	    (%%rbx),%%xmm2	\n\t"\
		"movaps	0x10(%%rdx),%%xmm7	\n\t"\
		"movaps	0x10(%%rbx),%%xmm3	\n\t"\
		"movaps	    (%%rdx),%%xmm6	\n\t"\
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
		"movaps	%%xmm0,    (%%rbx)	\n\t"\
		"movaps	%%xmm4,    (%%rcx)	\n\t"\
		"movaps	%%xmm1,0x10(%%rbx)	\n\t"\
		"movaps	%%xmm5,0x10(%%rdx)	\n\t"\
		"movaps	%%xmm2,    (%%rax)	\n\t"\
		"movaps	%%xmm7,    (%%rdx)	\n\t"\
		"movaps	%%xmm3,0x10(%%rax)	\n\t"\
		"movaps	%%xmm6,0x10(%%rcx)	\n\t"\
		"\n\t"\
		"/*...and now do eight radix-4 transforms, including the internal twiddle factors: */\n\t"\
		"/*...Block 1: t00,t10,t20,t30 in r00,04,02,06 - note swapped middle 2 indices! */\n\t"\
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
		"movaps	    (%%rsi),%%xmm0	/* t00 */\n\t"\
		"movaps	0x20(%%rsi),%%xmm4	/* t20 */\n\t"\
		"movaps	0x10(%%rsi),%%xmm1	/* t01 */\n\t"\
		"movaps	0x30(%%rsi),%%xmm5	/* t21 */\n\t"\
		"movaps	0x40(%%rsi),%%xmm2	/* t10 */\n\t"\
		"movaps	0x60(%%rsi),%%xmm6	/* t30 */\n\t"\
		"movaps	0x50(%%rsi),%%xmm3	/* t11 */\n\t"\
		"movaps	0x70(%%rsi),%%xmm7	/* t31 */\n\t"\
		"subpd	0x40(%%rsi),%%xmm0	/* t10=t00-rt */\n\t"\
		"subpd	0x60(%%rsi),%%xmm4	/* t30=t20-rt */\n\t"\
		"subpd	0x50(%%rsi),%%xmm1	/* t11=t01-it */\n\t"\
		"subpd	0x70(%%rsi),%%xmm5	/* t31=t21-it */\n\t"\
		"addpd	    (%%rsi),%%xmm2	/* t00=t00+rt */\n\t"\
		"addpd	0x20(%%rsi),%%xmm6	/* t20=t20+rt */\n\t"\
		"addpd	0x10(%%rsi),%%xmm3	/* t01=t01+it */\n\t"\
		"addpd	0x30(%%rsi),%%xmm7	/* t21=t21+it */\n\t"\
		"subpd	%%xmm6,%%xmm2		/* t00 <- t00-t20 */\n\t"\
		"subpd	%%xmm5,%%xmm0		/* t10 <- t10-t31 */\n\t"\
		"subpd	%%xmm7,%%xmm3		/* t01 <- t01-t21 */\n\t"\
		"subpd	%%xmm4,%%xmm1		/* t11 <- t11-t30 */\n\t"\
		"addpd	%%xmm6,%%xmm6		/*          2*t20 */\n\t"\
		"addpd	%%xmm5,%%xmm5		/*          2*t31 */\n\t"\
		"addpd	%%xmm7,%%xmm7		/*          2*t21 */\n\t"\
		"addpd	%%xmm4,%%xmm4		/*          2*t30 */\n\t"\
		"movaps	%%xmm2,    (%%rbx)	/* a(jt+p1 ) */\n\t"\
		"movaps	%%xmm0,    (%%rcx)	/* a(jt+p2 ) */\n\t"\
		"movaps	%%xmm3,0x10(%%rbx)	/* a(jp+p1 ) */\n\t"\
		"movaps	%%xmm1,0x10(%%rdx)	/* a(jp+p3 ) */\n\t"\
		"addpd	%%xmm2,%%xmm6		/* t20 <- t00+t20 */\n\t"\
		"addpd	%%xmm0,%%xmm5		/* t31 <- t10+t31 */\n\t"\
		"addpd	%%xmm3,%%xmm7		/* t21 <- t01+t21 */\n\t"\
		"addpd	%%xmm1,%%xmm4		/* t30 <- t11+t30 */\n\t"\
		"movaps	%%xmm6,    (%%rax)	/* a(jt+p0 ) */\n\t"\
		"movaps	%%xmm5,    (%%rdx)	/* a(jt+p3 ) */\n\t"\
		"movaps	%%xmm7,0x10(%%rax)	/* a(jp+p0 ) */\n\t"\
		"movaps	%%xmm4,0x10(%%rcx)	/* a(jp+p2 ) */\n\t"\
		"/*...Block 5: t08,t18,t28,t38	*/\n\t"\
		"addq	$0x80,%%rsi			/* r08 */\n\t"\
		"movslq	%[__p04],%%rdi		\n\t"\
		"shlq	$3,%%rdi			\n\t"\
		"addq	%%rdi,%%rax			/* ap04 */\n\t"\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	%%rdi,%%rcx			\n\t"\
		"addq	%%rdi,%%rdx			\n\t"\
		"movq	%[__isrt2],%%rdi	\n\t"\
		"\n\t"\
		"movaps	0x20(%%rsi),%%xmm4	/* t28 */\n\t"\
		"movaps	0x30(%%rsi),%%xmm5	/* t29 */\n\t"\
		"movaps	0x60(%%rsi),%%xmm6	/* t38 */\n\t"\
		"movaps	0x70(%%rsi),%%xmm7	/* t39 */\n\t"\
		"movaps	    (%%rsi),%%xmm0	/* t08 */\n\t"\
		"subpd	0x30(%%rsi),%%xmm4	/* t28-t29 */\n\t"\
		"movaps	0x10(%%rsi),%%xmm1	/* t09 */\n\t"\
		"addpd	0x20(%%rsi),%%xmm5	/* t29+t28 */\n\t"\
		"movaps	0x40(%%rsi),%%xmm2	/* t18 */\n\t"\
		"mulpd	(%%rdi),%%xmm4		/* t28 = (t28-t29)*ISRT2 */\n\t"\
		"movaps	0x50(%%rsi),%%xmm3	/* t19 */\n\t"\
		"mulpd	(%%rdi),%%xmm5		/* t29 = (t29+t28)*ISRT2 */\n\t"\
		"subpd	0x50(%%rsi),%%xmm0	/* t08=t08-t19*/\n\t"\
		"addpd	0x70(%%rsi),%%xmm6	/* t38+t39 */\n\t"\
		"subpd	0x40(%%rsi),%%xmm1	/* t19=t09-t18*/\n\t"\
		"subpd	0x60(%%rsi),%%xmm7	/* t39-t38 */\n\t"\
		"addpd	0x10(%%rsi),%%xmm2	/* t09=t18+t09*/\n\t"\
		"mulpd	(%%rdi),%%xmm6		/*  rt = (t38+t39)*ISRT2 */\n\t"\
		"addpd	    (%%rsi),%%xmm3	/* t18=t19+t08*/\n\t"\
		"mulpd	(%%rdi),%%xmm7		/*  it = (t39-t38)*ISRT2 */\n\t"\
		"\n\t"\
		"subpd	%%xmm6,%%xmm4		/* t28=t28-rt */\n\t"\
		"subpd	%%xmm7,%%xmm5		/* t29=t29-it */\n\t"\
		"addpd	%%xmm6,%%xmm6		/*      2* rt */\n\t"\
		"addpd	%%xmm7,%%xmm7		/*      2* it */\n\t"\
		"addpd	%%xmm4,%%xmm6		/* t38=t28+rt */\n\t"\
		"addpd	%%xmm5,%%xmm7		/* t39=t29+it */\n\t"\
		"subpd	%%xmm4,%%xmm0		/* t08-t28 */\n\t"\
		"subpd	%%xmm5,%%xmm2		/* t09-t29 */\n\t"\
		"addpd	%%xmm4,%%xmm4		/*   2*t28 */\n\t"\
		"addpd	%%xmm5,%%xmm5		/*   2*t29 */\n\t"\
		"\n\t"\
		"subpd	%%xmm7,%%xmm3		/* t18-t39 */\n\t"\
		"subpd	%%xmm6,%%xmm1		/* t19-t38 */\n\t"\
		"addpd	%%xmm7,%%xmm7		/*   2*t39 */\n\t"\
		"addpd	%%xmm6,%%xmm6		/*   2*t38 */\n\t"\
		"movaps	%%xmm0,    (%%rbx)	/* a(jt+p1 ) */\n\t"\
		"movaps	%%xmm2,0x10(%%rbx)	/* a(jp+p1 ) */\n\t"\
		"movaps	%%xmm3,    (%%rcx)	/* a(jt+p2 ) */\n\t"\
		"movaps	%%xmm1,0x10(%%rdx)	/* a(jp+p3 ) */\n\t"\
		"addpd	%%xmm0,%%xmm4		/* t08+t28 */\n\t"\
		"addpd	%%xmm2,%%xmm5		/* t09+t29 */\n\t"\
		"addpd	%%xmm3,%%xmm7		/* t18+t39 */\n\t"\
		"addpd	%%xmm1,%%xmm6		/* t19+t38 */\n\t"\
		"movaps	%%xmm4,    (%%rax)	/* a(jt+p0 ) */\n\t"\
		"movaps	%%xmm5,0x10(%%rax)	/* a(jp+p0 ) */\n\t"\
		"movaps	%%xmm7,    (%%rdx)	/* a(jt+p3 ) */\n\t"\
		"movaps	%%xmm6,0x10(%%rcx)	/* a(jp+p2 ) */\n\t"\
		"/*...Block 3: t04,t14,t24,t34	*/\n\t"\
		"addq	$0x180,%%rsi		/* r20 */\n\t"\
		"movslq	%[__p04],%%rdi		\n\t"\
		"shlq	$3,%%rdi			\n\t"\
		"addq	%%rdi,%%rax			/* ap08 */\n\t"\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	%%rdi,%%rcx			\n\t"\
		"addq	%%rdi,%%rdx			\n\t"\
		"movq	%[__cc0],%%rdi	\n\t"\
		"\n\t"\
		"movaps	0x20(%%rsi),%%xmm4	/* t24 */\n\t"\
		"movaps	0x60(%%rsi),%%xmm6	/* t34 */\n\t"\
		"movaps	0x30(%%rsi),%%xmm5	/* t25 */\n\t"\
		"movaps	0x70(%%rsi),%%xmm7	/* t35 */\n\t"\
		"movaps	0x20(%%rsi),%%xmm0	/* copy t24 */\n\t"\
		"movaps	0x60(%%rsi),%%xmm2	/* copy t34 */\n\t"\
		"movaps	0x30(%%rsi),%%xmm1	/* copy t25 */\n\t"\
		"movaps	0x70(%%rsi),%%xmm3	/* copy t35 */\n\t"\
		"mulpd	    (%%rdi),%%xmm4	/* t24*c */\n\t"\
		"mulpd	0x10(%%rdi),%%xmm6	/* t34*s */\n\t"\
		"mulpd	0x10(%%rdi),%%xmm1	/* t25*s */\n\t"\
		"mulpd	    (%%rdi),%%xmm3	/* t35*c */\n\t"\
		"mulpd	    (%%rdi),%%xmm5	/* t25*c */\n\t"\
		"mulpd	0x10(%%rdi),%%xmm7	/* t35*s */\n\t"\
		"mulpd	0x10(%%rdi),%%xmm0	/* t24*s */\n\t"\
		"mulpd	    (%%rdi),%%xmm2	/* t34*c */\n\t"\
		"subpd	%%xmm1,%%xmm4		/* ~t24 */\n\t"\
		"subpd	%%xmm3,%%xmm6		/* rt */\n\t"\
		"addpd	%%xmm0,%%xmm5		/* ~t25 */\n\t"\
		"addpd	%%xmm2,%%xmm7		/* it */\n\t"\
		"subq	$0x10,%%rdi			/* isrt2 */\n\t"\
		"\n\t"\
		"movaps	0x40(%%rsi),%%xmm2	/* t14 */\n\t"\
		"subpd	%%xmm6,%%xmm4		/*~t34=t24-rt */\n\t"\
		"movaps	0x50(%%rsi),%%xmm3	/* t15 */\n\t"\
		"subpd	%%xmm7,%%xmm5		/*~t35=t25-it */\n\t"\
		"subpd	0x50(%%rsi),%%xmm2	/* t14-t15 */\n\t"\
		"addpd	%%xmm6,%%xmm6		/*      2* rt */\n\t"\
		"addpd	0x40(%%rsi),%%xmm3	/* t15+t14 */\n\t"\
		"addpd	%%xmm7,%%xmm7		/*      2* it */\n\t"\
		"mulpd	(%%rdi),%%xmm2		/* rt = (t14-t15)*ISRT2 */\n\t"\
		"addpd	%%xmm4,%%xmm6		/*~t24=t24+rt */\n\t"\
		"mulpd	(%%rdi),%%xmm3		/* it = (t15+t14)*ISRT2 */\n\t"\
		"addpd	%%xmm5,%%xmm7		/*~t25=t25+it */\n\t"\
		"movaps	    (%%rsi),%%xmm0	/* t04 */\n\t"\
		"movaps	0x10(%%rsi),%%xmm1	/* t05 */\n\t"\
		"subpd	%%xmm2,%%xmm0		/*~t14=t04-rt */\n\t"\
		"subpd	%%xmm3,%%xmm1		/*~t15=t05-it */\n\t"\
		"addpd	    (%%rsi),%%xmm2	/*~t04=rt +t04*/\n\t"\
		"addpd	0x10(%%rsi),%%xmm3	/*~t05=it +t05*/\n\t"\
		"subpd	%%xmm6,%%xmm2		/* t04-t24 */\n\t"\
		"subpd	%%xmm5,%%xmm0		/* t14-t35 */\n\t"\
		"subpd	%%xmm7,%%xmm3		/* t05-t25 */\n\t"\
		"subpd	%%xmm4,%%xmm1		/* t15-t34 */\n\t"\
		"addpd	%%xmm6,%%xmm6		/*   2*t24 */\n\t"\
		"addpd	%%xmm5,%%xmm5		/*          2*t35 */\n\t"\
		"addpd	%%xmm7,%%xmm7		/*   2*t25 */\n\t"\
		"addpd	%%xmm4,%%xmm4		/*          2*t34 */\n\t"\
		"movaps	%%xmm2,    (%%rbx)	/* a(jt+p1 ) */\n\t"\
		"movaps	%%xmm0,    (%%rcx)	/* a(jt+p2 ) */\n\t"\
		"movaps	%%xmm3,0x10(%%rbx)	/* a(jp+p1 ) */\n\t"\
		"movaps	%%xmm1,0x10(%%rdx)	/* a(jp+p3 ) */\n\t"\
		"addpd	%%xmm2,%%xmm6		/* t04+t24 */\n\t"\
		"addpd	%%xmm0,%%xmm5		/* t14+t35 */\n\t"\
		"addpd	%%xmm3,%%xmm7		/* t05+t25 */\n\t"\
		"addpd	%%xmm1,%%xmm4		/* t15+t34 */\n\t"\
		"movaps	%%xmm6,    (%%rax)	/* a(jt+p0 ) */\n\t"\
		"movaps	%%xmm5,    (%%rdx)	/* a(jt+p3 ) */\n\t"\
		"movaps	%%xmm7,0x10(%%rax)	/* a(jp+p0 ) */\n\t"\
		"movaps	%%xmm4,0x10(%%rcx)	/* a(jp+p2 ) */\n\t"\
		"/*...Block 7: t0C,t1C,t2C,t3C	*/\n\t"\
		"addq	$0x80,%%rsi			/* r28 */\n\t"\
		"movslq	%[__p04],%%rdi		\n\t"\
		"shlq	$3,%%rdi			\n\t"\
		"addq	%%rdi,%%rax			/* ap0C */\n\t"\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	%%rdi,%%rcx			\n\t"\
		"addq	%%rdi,%%rdx			\n\t"\
		"movq	%[__cc0],%%rdi	\n\t"\
		"\n\t"\
		"movaps	0x20(%%rsi),%%xmm4	/* t2C */\n\t"\
		"movaps	0x60(%%rsi),%%xmm6	/* t3C */\n\t"\
		"movaps	0x30(%%rsi),%%xmm5	/* t2D */\n\t"\
		"movaps	0x70(%%rsi),%%xmm7	/* t3D */\n\t"\
		"movaps	0x20(%%rsi),%%xmm0	/* copy t2C */\n\t"\
		"movaps	0x60(%%rsi),%%xmm2	/* copy t3C */\n\t"\
		"movaps	0x30(%%rsi),%%xmm1	/* copy t2D */\n\t"\
		"movaps	0x70(%%rsi),%%xmm3	/* copy t3D */\n\t"\
		"mulpd	0x10(%%rdi),%%xmm4	/* t2C*s */\n\t"\
		"mulpd	    (%%rdi),%%xmm6	/* t3C*c */\n\t"\
		"mulpd	    (%%rdi),%%xmm1	/* t2D*c */\n\t"\
		"mulpd	0x10(%%rdi),%%xmm3	/* t3D*s */\n\t"\
		"mulpd	0x10(%%rdi),%%xmm5	/* t2D*s */\n\t"\
		"mulpd	    (%%rdi),%%xmm7	/* t3D*c */\n\t"\
		"mulpd	    (%%rdi),%%xmm0	/* t2C*c */\n\t"\
		"mulpd	0x10(%%rdi),%%xmm2	/* t3C*s */\n\t"\
		"subpd	%%xmm1,%%xmm4		/* ~t24 */\n\t"\
		"subpd	%%xmm3,%%xmm6		/* rt */\n\t"\
		"addpd	%%xmm0,%%xmm5		/* ~t25 */\n\t"\
		"addpd	%%xmm2,%%xmm7		/* it */\n\t"\
		"subq	$0x10,%%rdi			/* isrt2 */\n\t"\
		"\n\t"\
		"movaps	0x40(%%rsi),%%xmm2	/* t14 */\n\t"\
		"subpd	%%xmm6,%%xmm4		/*~t2C=t2C-rt */\n\t"\
		"movaps	0x50(%%rsi),%%xmm3	/* t1D */\n\t"\
		"subpd	%%xmm7,%%xmm5		/*~t2D=t2D-it */\n\t"\
		"addpd	0x50(%%rsi),%%xmm2	/* t1C+t1D */\n\t"\
		"addpd	%%xmm6,%%xmm6		/*      2* rt */\n\t"\
		"subpd	0x40(%%rsi),%%xmm3	/* t1D-t1C */\n\t"\
		"addpd	%%xmm7,%%xmm7		/*      2* it */\n\t"\
		"mulpd	(%%rdi),%%xmm2		/* rt = (t1C+t1D)*ISRT2 */\n\t"\
		"addpd	%%xmm4,%%xmm6		/*~t3C=t2C+rt */\n\t"\
		"mulpd	(%%rdi),%%xmm3		/* it = (t1D-t1C)*ISRT2 */\n\t"\
		"addpd	%%xmm5,%%xmm7		/*~t3D=t2D+it */\n\t"\
		"movaps	    (%%rsi),%%xmm0	/* t0C */\n\t"\
		"movaps	0x10(%%rsi),%%xmm1	/* t0D */\n\t"\
		"subpd	%%xmm2,%%xmm0		/*~t0C=t0C-rt */\n\t"\
		"subpd	%%xmm3,%%xmm1		/*~t0D=t0D-it */\n\t"\
		"addpd	    (%%rsi),%%xmm2	/*~t1C=rt +t0C*/\n\t"\
		"addpd	0x10(%%rsi),%%xmm3	/*~t1D=it +t0D*/\n\t"\
		"subpd	%%xmm4,%%xmm0		/* t0C-t2C */\n\t"\
		"subpd	%%xmm7,%%xmm2		/* t1C-t3D */\n\t"\
		"subpd	%%xmm5,%%xmm1		/* t0D-t2D */\n\t"\
		"subpd	%%xmm6,%%xmm3		/* t1D-t3C */\n\t"\
		"addpd	%%xmm4,%%xmm4		/*   2*t2C */\n\t"\
		"addpd	%%xmm7,%%xmm7		/*   2*t3D */\n\t"\
		"addpd	%%xmm5,%%xmm5		/*   2*t2D */\n\t"\
		"addpd	%%xmm6,%%xmm6		/*   2*t3C */\n\t"\
		"movaps	%%xmm0,    (%%rbx)	/* a(jt+p1 ) */\n\t"\
		"movaps	%%xmm2,    (%%rcx)	/* a(jt+p2 ) */\n\t"\
		"movaps	%%xmm1,0x10(%%rbx)	/* a(jp+p1 ) */\n\t"\
		"movaps	%%xmm3,0x10(%%rdx)	/* a(jp+p3 ) */\n\t"\
		"addpd	%%xmm0,%%xmm4		/* t0C+t2C */\n\t"\
		"addpd	%%xmm2,%%xmm7		/* t1C+t3D */\n\t"\
		"addpd	%%xmm1,%%xmm5		/* t0D+t2D */\n\t"\
		"addpd	%%xmm3,%%xmm6		/* t1D+t3C */\n\t"\
		"movaps	%%xmm4,    (%%rax)	/* a(jt+p0 ) */\n\t"\
		"movaps	%%xmm7,    (%%rdx)	/* a(jt+p3 ) */\n\t"\
		"movaps	%%xmm5,0x10(%%rax)	/* a(jp+p0 ) */\n\t"\
		"movaps	%%xmm6,0x10(%%rcx)	/* a(jp+p2 ) */\n\t"\
		"/*...Block 2: t02,t12,t22,t32	*/\n\t"\
		"subq	$0x180,%%rsi		/* r10 */\n\t"\
		"movslq	%[__p04],%%rdi		\n\t"\
		"shlq	$3,%%rdi			\n\t"\
		"addq	%%rdi,%%rax			/* ap10 */\n\t"\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	%%rdi,%%rcx			\n\t"\
		"addq	%%rdi,%%rdx			\n\t"\
		"movq	%[__cc0],%%rdi	\n\t"\
		"\n\t"\
		"movaps	0x20(%%rsi),%%xmm4	/* t22 */\n\t"\
		"movaps	0x60(%%rsi),%%xmm6	/* t32 */\n\t"\
		"movaps	0x30(%%rsi),%%xmm5	/* t23 */\n\t"\
		"movaps	0x70(%%rsi),%%xmm7	/* t33 */\n\t"\
		"movaps	0x20(%%rsi),%%xmm0	/* copy t22 */\n\t"\
		"movaps	0x60(%%rsi),%%xmm2	/* copy t32 */\n\t"\
		"movaps	0x30(%%rsi),%%xmm1	/* copy t23 */\n\t"\
		"movaps	0x70(%%rsi),%%xmm3	/* copy t33 */\n\t"\
		"mulpd	0x20(%%rdi),%%xmm4	/* t22*c32_1 */\n\t"\
		"mulpd	0x40(%%rdi),%%xmm6	/* t32*c32_3 */\n\t"\
		"mulpd	0x30(%%rdi),%%xmm1	/* t23*s32_1 */\n\t"\
		"mulpd	0x50(%%rdi),%%xmm3	/* t33*s32_3 */\n\t"\
		"mulpd	0x20(%%rdi),%%xmm5	/* t23*c32_1 */\n\t"\
		"mulpd	0x40(%%rdi),%%xmm7	/* t33*c32_3 */\n\t"\
		"mulpd	0x30(%%rdi),%%xmm0	/* t22*s32_1 */\n\t"\
		"mulpd	0x50(%%rdi),%%xmm2	/* t32*s32_3 */\n\t"\
		"subpd	%%xmm1,%%xmm4		/* ~t22 */\n\t"\
		"subpd	%%xmm3,%%xmm6		/* rt */\n\t"\
		"addpd	%%xmm0,%%xmm5		/* ~t23 */\n\t"\
		"addpd	%%xmm2,%%xmm7		/* it */\n\t"\
		"\n\t"\
		"movaps	0x40(%%rsi),%%xmm2	/* t12 */\n\t"\
		"movaps	0x50(%%rsi),%%xmm0	/* t13 */\n\t"\
		"movaps	0x40(%%rsi),%%xmm1	/* copy t12 */\n\t"\
		"movaps	0x50(%%rsi),%%xmm3	/* copy t13 */\n\t"\
		"subpd	%%xmm6,%%xmm4		/*~t32=t22-rt */\n\t"\
		"mulpd	    (%%rdi),%%xmm2	/* t12*c */\n\t"\
		"subpd	%%xmm7,%%xmm5		/*~t33=t23-it */\n\t"\
		"mulpd	0x10(%%rdi),%%xmm0	/* t13*s */\n\t"\
		"addpd	%%xmm6,%%xmm6		/*      2* rt */\n\t"\
		"mulpd	    (%%rdi),%%xmm3	/* t13*c */\n\t"\
		"addpd	%%xmm7,%%xmm7		/*      2* it */\n\t"\
		"mulpd	0x10(%%rdi),%%xmm1	/* t12*s */\n\t"\
		"addpd	%%xmm4,%%xmm6		/*~t22=t22+rt */\n\t"\
		"subpd	%%xmm0,%%xmm2		/* rt */\n\t"\
		"addpd	%%xmm5,%%xmm7		/*~t23=t23+it */\n\t"\
		"addpd	%%xmm1,%%xmm3		/* it */\n\t"\
		"\n\t"\
		"movaps	    (%%rsi),%%xmm0	/* t02 */\n\t"\
		"movaps	0x10(%%rsi),%%xmm1	/* t03 */\n\t"\
		"subpd	%%xmm2,%%xmm0		/*~t12=t02-rt */\n\t"\
		"subpd	%%xmm3,%%xmm1		/*~t13=t03-it */\n\t"\
		"addpd	    (%%rsi),%%xmm2	/*~t02=rt+t02 */\n\t"\
		"addpd	0x10(%%rsi),%%xmm3	/*~t03=it+t03 */\n\t"\
		"subpd	%%xmm6,%%xmm2		/* t02-t22 */\n\t"\
		"subpd	%%xmm5,%%xmm0		/* t12-t33 */\n\t"\
		"subpd	%%xmm7,%%xmm3		/* t03-t23 */\n\t"\
		"subpd	%%xmm4,%%xmm1		/* t13-t32 */\n\t"\
		"addpd	%%xmm6,%%xmm6		/*   2*t22 */\n\t"\
		"addpd	%%xmm5,%%xmm5		/*   2*t33 */\n\t"\
		"addpd	%%xmm7,%%xmm7		/*   2*t23 */\n\t"\
		"addpd	%%xmm4,%%xmm4		/*   2*t32 */\n\t"\
		"movaps	%%xmm2,    (%%rbx)	/* a(jt+p1 ) */\n\t"\
		"movaps	%%xmm0,    (%%rcx)	/* a(jt+p2 ) */\n\t"\
		"movaps	%%xmm3,0x10(%%rbx)	/* a(jp+p1 ) */\n\t"\
		"movaps	%%xmm1,0x10(%%rdx)	/* a(jp+p3 ) */\n\t"\
		"addpd	%%xmm2,%%xmm6		/* t02+t22 */\n\t"\
		"addpd	%%xmm0,%%xmm5		/* t12+t33 */\n\t"\
		"addpd	%%xmm3,%%xmm7		/* t03+t23 */\n\t"\
		"addpd	%%xmm1,%%xmm4		/* t13+t32 */\n\t"\
		"movaps	%%xmm6,    (%%rax)	/* a(jt+p0 ) */\n\t"\
		"movaps	%%xmm5,    (%%rdx)	/* a(jt+p3 ) */\n\t"\
		"movaps	%%xmm7,0x10(%%rax)	/* a(jp+p0 ) */\n\t"\
		"movaps	%%xmm4,0x10(%%rcx)	/* a(jp+p2 ) */\n\t"\
		"/*...Block 6: t0A,t1A,t2A,t3A	*/\n\t"\
		"addq	$0x80,%%rsi			/* r18 */\n\t"\
		"movslq	%[__p04],%%rdi		\n\t"\
		"shlq	$3,%%rdi			\n\t"\
		"addq	%%rdi,%%rax			/* ap14 */\n\t"\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	%%rdi,%%rcx			\n\t"\
		"addq	%%rdi,%%rdx			\n\t"\
		"movq	%[__cc0],%%rdi	\n\t"\
		"\n\t"\
		"movaps	0x20(%%rsi),%%xmm4	/* t2A */\n\t"\
		"movaps	0x60(%%rsi),%%xmm6	/* t3A */\n\t"\
		"movaps	0x30(%%rsi),%%xmm5	/* t2B */\n\t"\
		"movaps	0x70(%%rsi),%%xmm7	/* t3B */\n\t"\
		"movaps	0x20(%%rsi),%%xmm0	/* copy t2A */\n\t"\
		"movaps	0x60(%%rsi),%%xmm2	/* copy t3A */\n\t"\
		"movaps	0x30(%%rsi),%%xmm1	/* copy t2B */\n\t"\
		"movaps	0x70(%%rsi),%%xmm3	/* copy t3B */\n\t"\
		"mulpd	0x50(%%rdi),%%xmm4	/* t2A*s32_3 */\n\t"\
		"mulpd	0x20(%%rdi),%%xmm6	/* t3A*c32_1 */\n\t"\
		"mulpd	0x40(%%rdi),%%xmm1	/* t2B*c32_3 */\n\t"\
		"mulpd	0x30(%%rdi),%%xmm3	/* t3B*s32_1 */\n\t"\
		"mulpd	0x50(%%rdi),%%xmm5	/* t2B*s32_3 */\n\t"\
		"mulpd	0x20(%%rdi),%%xmm7	/* t3B*c32_1 */\n\t"\
		"mulpd	0x40(%%rdi),%%xmm0	/* t2A*c32_3 */\n\t"\
		"mulpd	0x30(%%rdi),%%xmm2	/* t3A*s32_1 */\n\t"\
		"subpd	%%xmm1,%%xmm4		/* ~t2A */\n\t"\
		"addpd	%%xmm3,%%xmm6		/* rt */\n\t"\
		"addpd	%%xmm0,%%xmm5		/* ~t2B */\n\t"\
		"subpd	%%xmm2,%%xmm7		/* it */\n\t"\
		"\n\t"\
		"movaps	0x40(%%rsi),%%xmm2	/* t1A */\n\t"\
		"movaps	0x50(%%rsi),%%xmm0	/* t1B */\n\t"\
		"movaps	0x40(%%rsi),%%xmm1	/* copy t1A */\n\t"\
		"movaps	0x50(%%rsi),%%xmm3	/* copy t1B */\n\t"\
		"subpd	%%xmm6,%%xmm4		/*~t2A=t2A-rt */\n\t"\
		"mulpd	0x10(%%rdi),%%xmm2	/* t1A*s */\n\t"\
		"subpd	%%xmm7,%%xmm5		/*~t2B=t2B-it */\n\t"\
		"mulpd	    (%%rdi),%%xmm0	/* t1B*c */\n\t"\
		"addpd	%%xmm6,%%xmm6		/*      2* rt */\n\t"\
		"mulpd	0x10(%%rdi),%%xmm3	/* t1B*s */\n\t"\
		"addpd	%%xmm7,%%xmm7		/*      2* it */\n\t"\
		"mulpd	    (%%rdi),%%xmm1	/* t1A*c */\n\t"\
		"addpd	%%xmm4,%%xmm6		/*~t3A=t2A+rt */\n\t"\
		"addpd	%%xmm0,%%xmm2		/* rt */\n\t"\
		"addpd	%%xmm5,%%xmm7		/*~t3B=t2B+it */\n\t"\
		"subpd	%%xmm1,%%xmm3		/* it */\n\t"\
		"\n\t"\
		"movaps	    (%%rsi),%%xmm0	/* t0A */\n\t"\
		"movaps	0x10(%%rsi),%%xmm1	/* t0B */\n\t"\
		"subpd	%%xmm2,%%xmm0		/*~t0A=t0A-rt */\n\t"\
		"subpd	%%xmm3,%%xmm1		/*~t0B=t0B-it */\n\t"\
		"addpd	    (%%rsi),%%xmm2	/*~t1A=rt+t0A */\n\t"\
		"addpd	0x10(%%rsi),%%xmm3	/*~t1B=it+t0B */\n\t"\
		"subpd	%%xmm4,%%xmm0		/* t0A-t2A */\n\t"\
		"subpd	%%xmm7,%%xmm2		/* t1A-t3B */\n\t"\
		"subpd	%%xmm5,%%xmm1		/* t0B-t2B */\n\t"\
		"subpd	%%xmm6,%%xmm3		/* t1B-t3A */\n\t"\
		"addpd	%%xmm4,%%xmm4		/*   2*t2A */\n\t"\
		"addpd	%%xmm7,%%xmm7		/*   2*t3B */\n\t"\
		"addpd	%%xmm5,%%xmm5		/*   2*t2B */\n\t"\
		"addpd	%%xmm6,%%xmm6		/*   2*t3A */\n\t"\
		"movaps	%%xmm0,    (%%rbx)	/* a(jt+p1 ) */\n\t"\
		"movaps	%%xmm2,    (%%rcx)	/* a(jt+p2 ) */\n\t"\
		"movaps	%%xmm1,0x10(%%rbx)	/* a(jp+p1 ) */\n\t"\
		"movaps	%%xmm3,0x10(%%rdx)	/* a(jp+p3 ) */\n\t"\
		"addpd	%%xmm0,%%xmm4		/* t0A+t2A */\n\t"\
		"addpd	%%xmm2,%%xmm7		/* t1A+t3B */\n\t"\
		"addpd	%%xmm1,%%xmm5		/* t0B+t2B */\n\t"\
		"addpd	%%xmm3,%%xmm6		/* t1B+t3A */\n\t"\
		"movaps	%%xmm4,    (%%rax)	/* a(jt+p0 ) */\n\t"\
		"movaps	%%xmm7,    (%%rdx)	/* a(jt+p3 ) */\n\t"\
		"movaps	%%xmm5,0x10(%%rax)	/* a(jp+p0 ) */\n\t"\
		"movaps	%%xmm6,0x10(%%rcx)	/* a(jp+p2 ) */\n\t"\
		"/*...Block 4: t06,t16,t26,t36	*/\n\t"\
		"addq	$0x180,%%rsi		/* r30 */\n\t"\
		"movslq	%[__p04],%%rdi		\n\t"\
		"shlq	$3,%%rdi			\n\t"\
		"addq	%%rdi,%%rax			/* ap18 */\n\t"\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	%%rdi,%%rcx			\n\t"\
		"addq	%%rdi,%%rdx			\n\t"\
		"movq	%[__cc0],%%rdi	\n\t"\
		"\n\t"\
		"movaps	0x20(%%rsi),%%xmm4	/* t26 */\n\t"\
		"movaps	0x60(%%rsi),%%xmm6	/* t36 */\n\t"\
		"movaps	0x30(%%rsi),%%xmm5	/* t27 */\n\t"\
		"movaps	0x70(%%rsi),%%xmm7	/* t37 */\n\t"\
		"movaps	0x20(%%rsi),%%xmm0	/* copy t26 */\n\t"\
		"movaps	0x60(%%rsi),%%xmm2	/* copy t36 */\n\t"\
		"movaps	0x30(%%rsi),%%xmm1	/* copy t27 */\n\t"\
		"movaps	0x70(%%rsi),%%xmm3	/* copy t37 */\n\t"\
		"mulpd	0x40(%%rdi),%%xmm4	/* t26*s32_3 */\n\t"\
		"mulpd	0x30(%%rdi),%%xmm6	/* t36*s32_1 */\n\t"\
		"mulpd	0x50(%%rdi),%%xmm1	/* t27*s32_3 */\n\t"\
		"mulpd	0x20(%%rdi),%%xmm3	/* t37*c32_1 */\n\t"\
		"mulpd	0x40(%%rdi),%%xmm5	/* t27*c32_3 */\n\t"\
		"mulpd	0x30(%%rdi),%%xmm7	/* t37*s32_1 */\n\t"\
		"mulpd	0x50(%%rdi),%%xmm0	/* t26*s32_3 */\n\t"\
		"mulpd	0x20(%%rdi),%%xmm2	/* t36*c32_1 */\n\t"\
		"subpd	%%xmm1,%%xmm4		/* ~t26 */\n\t"\
		"addpd	%%xmm3,%%xmm6		/* rt */\n\t"\
		"addpd	%%xmm0,%%xmm5		/* ~t27 */\n\t"\
		"subpd	%%xmm2,%%xmm7		/* it */\n\t"\
		"\n\t"\
		"movaps	0x40(%%rsi),%%xmm2	/* t16 */\n\t"\
		"movaps	0x50(%%rsi),%%xmm0	/* t17 */\n\t"\
		"movaps	0x40(%%rsi),%%xmm1	/* copy t16 */\n\t"\
		"movaps	0x50(%%rsi),%%xmm3	/* copy t17 */\n\t"\
		"subpd	%%xmm6,%%xmm4		/*~t26=t26-rt */\n\t"\
		"mulpd	0x10(%%rdi),%%xmm2	/* t16*s */\n\t"\
		"subpd	%%xmm7,%%xmm5		/*~t27=t27-it */\n\t"\
		"mulpd	    (%%rdi),%%xmm0	/* t17*c */\n\t"\
		"addpd	%%xmm6,%%xmm6		/*      2* rt */\n\t"\
		"mulpd	0x10(%%rdi),%%xmm3	/* t17*s */\n\t"\
		"addpd	%%xmm7,%%xmm7		/*      2* it */\n\t"\
		"mulpd	    (%%rdi),%%xmm1	/* t16*c */\n\t"\
		"addpd	%%xmm4,%%xmm6		/*~t36=t26+rt */\n\t"\
		"subpd	%%xmm0,%%xmm2		/* rt */\n\t"\
		"addpd	%%xmm5,%%xmm7		/*~t37=t27+it */\n\t"\
		"addpd	%%xmm1,%%xmm3		/* it */\n\t"\
		"\n\t"\
		"movaps	    (%%rsi),%%xmm0	/* t06 */\n\t"\
		"movaps	0x10(%%rsi),%%xmm1	/* t07 */\n\t"\
		"subpd	%%xmm2,%%xmm0		/*~t16=t06-rt */\n\t"\
		"subpd	%%xmm3,%%xmm1		/*~t17=t07-it */\n\t"\
		"addpd	    (%%rsi),%%xmm2	/*~t06=rt+t06 */\n\t"\
		"addpd	0x10(%%rsi),%%xmm3	/*~t07=it+t07 */\n\t"\
		"subpd	%%xmm4,%%xmm2		/* t06-t26 */\n\t"\
		"subpd	%%xmm7,%%xmm0		/* t16-t37 */\n\t"\
		"subpd	%%xmm5,%%xmm3		/* t07-t27 */\n\t"\
		"subpd	%%xmm6,%%xmm1		/* t17-t36 */\n\t"\
		"addpd	%%xmm4,%%xmm4		/*   2*t26 */\n\t"\
		"addpd	%%xmm7,%%xmm7		/*   2*t37 */\n\t"\
		"addpd	%%xmm5,%%xmm5		/*   2*t27 */\n\t"\
		"addpd	%%xmm6,%%xmm6		/*   2*t36 */\n\t"\
		"movaps	%%xmm2,    (%%rbx)	/* a(jt+p1 ) */\n\t"\
		"movaps	%%xmm0,    (%%rcx)	/* a(jt+p2 ) */\n\t"\
		"movaps	%%xmm3,0x10(%%rbx)	/* a(jp+p1 ) */\n\t"\
		"movaps	%%xmm1,0x10(%%rdx)	/* a(jp+p3 ) */\n\t"\
		"addpd	%%xmm2,%%xmm4		/* t06+t26 */\n\t"\
		"addpd	%%xmm0,%%xmm7		/* t16+t37 */\n\t"\
		"addpd	%%xmm3,%%xmm5		/* t07+t27 */\n\t"\
		"addpd	%%xmm1,%%xmm6		/* t17+t36 */\n\t"\
		"movaps	%%xmm4,    (%%rax)	/* a(jt+p0 ) */\n\t"\
		"movaps	%%xmm7,    (%%rdx)	/* a(jt+p3 ) */\n\t"\
		"movaps	%%xmm5,0x10(%%rax)	/* a(jp+p0 ) */\n\t"\
		"movaps	%%xmm6,0x10(%%rcx)	/* a(jp+p2 ) */\n\t"\
		"/*...Block 8: t0E,t1E,t2E,t3E	*/\n\t"\
		"addq	$0x80,%%rsi			/* r38 */\n\t"\
		"movslq	%[__p04],%%rdi		\n\t"\
		"shlq	$3,%%rdi			\n\t"\
		"addq	%%rdi,%%rax			/* ap1C */\n\t"\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	%%rdi,%%rcx			\n\t"\
		"addq	%%rdi,%%rdx			\n\t"\
		"movq	%[__cc0],%%rdi	\n\t"\
		"\n\t"\
		"movaps	0x20(%%rsi),%%xmm4	/* t2E */\n\t"\
		"movaps	0x60(%%rsi),%%xmm6	/* t3E */\n\t"\
		"movaps	0x30(%%rsi),%%xmm5	/* t2F */\n\t"\
		"movaps	0x70(%%rsi),%%xmm7	/* t3F */\n\t"\
		"movaps	0x20(%%rsi),%%xmm0	/* copy t2E */\n\t"\
		"movaps	0x60(%%rsi),%%xmm2	/* copy t3E */\n\t"\
		"movaps	0x30(%%rsi),%%xmm1	/* copy t2F */\n\t"\
		"movaps	0x70(%%rsi),%%xmm3	/* copy t3F */\n\t"\
		"mulpd	0x30(%%rdi),%%xmm4	/* t2E*s32_1 */\n\t"\
		"mulpd	0x50(%%rdi),%%xmm6	/* t3E*c32_3 */\n\t"\
		"mulpd	0x20(%%rdi),%%xmm1	/* t2F*c32_1 */\n\t"\
		"mulpd	0x40(%%rdi),%%xmm3	/* t3F*s32_3 */\n\t"\
		"mulpd	0x30(%%rdi),%%xmm5	/* t2F*s32_1 */\n\t"\
		"mulpd	0x50(%%rdi),%%xmm7	/* t3F*c32_3 */\n\t"\
		"mulpd	0x20(%%rdi),%%xmm0	/* t2E*c32_1 */\n\t"\
		"mulpd	0x40(%%rdi),%%xmm2	/* t3E*s32_3 */\n\t"\
		"subpd	%%xmm1,%%xmm4		/* ~t2E */\n\t"\
		"subpd	%%xmm3,%%xmm6		/* rt */\n\t"\
		"addpd	%%xmm0,%%xmm5		/* ~t2F */\n\t"\
		"addpd	%%xmm2,%%xmm7		/* it */\n\t"\
		"\n\t"\
		"movaps	0x40(%%rsi),%%xmm2	/* t1E */\n\t"\
		"movaps	0x50(%%rsi),%%xmm0	/* t1F */\n\t"\
		"movaps	0x40(%%rsi),%%xmm1	/* copy t1E */\n\t"\
		"movaps	0x50(%%rsi),%%xmm3	/* copy t1F */\n\t"\
		"subpd	%%xmm6,%%xmm4		/*~t2E=t2E-rt */\n\t"\
		"mulpd	    (%%rdi),%%xmm2	/* t1E*c */\n\t"\
		"subpd	%%xmm7,%%xmm5		/*~t2F=t2F-it */\n\t"\
		"mulpd	0x10(%%rdi),%%xmm0	/* t1F*s */\n\t"\
		"addpd	%%xmm6,%%xmm6		/*      2* rt */\n\t"\
		"mulpd	    (%%rdi),%%xmm3	/* t1F*c */\n\t"\
		"addpd	%%xmm7,%%xmm7		/*      2* it */\n\t"\
		"mulpd	0x10(%%rdi),%%xmm1	/* t1E*s */\n\t"\
		"addpd	%%xmm4,%%xmm6		/*~t3E=t2E+rt */\n\t"\
		"addpd	%%xmm0,%%xmm2		/* rt */\n\t"\
		"addpd	%%xmm5,%%xmm7		/*~t3F=t2F+it */\n\t"\
		"subpd	%%xmm1,%%xmm3		/* it */\n\t"\
		"\n\t"\
		"movaps	    (%%rsi),%%xmm0	/* t0E */\n\t"\
		"movaps	0x10(%%rsi),%%xmm1	/* t0F */\n\t"\
		"subpd	%%xmm2,%%xmm0		/*~t0E=t0E-rt */\n\t"\
		"subpd	%%xmm3,%%xmm1		/*~t0F=t0F-it */\n\t"\
		"addpd	    (%%rsi),%%xmm2	/*~t1E=rt+t0E */\n\t"\
		"addpd	0x10(%%rsi),%%xmm3	/*~t1F=it+t0F */\n\t"\
		"subpd	%%xmm4,%%xmm0		/* t0E-t2E */\n\t"\
		"subpd	%%xmm7,%%xmm2		/* t1E-t3F */\n\t"\
		"subpd	%%xmm5,%%xmm1		/* t0F-t2F */\n\t"\
		"subpd	%%xmm6,%%xmm3		/* t1F-t3E */\n\t"\
		"addpd	%%xmm4,%%xmm4		/*   2*t2E */\n\t"\
		"addpd	%%xmm7,%%xmm7		/*   2*t3F */\n\t"\
		"addpd	%%xmm5,%%xmm5		/*   2*t2F */\n\t"\
		"addpd	%%xmm6,%%xmm6		/*   2*t3E */\n\t"\
		"movaps	%%xmm0,    (%%rbx)	/* a(jt+p1 ) */\n\t"\
		"movaps	%%xmm2,    (%%rcx)	/* a(jt+p2 ) */\n\t"\
		"movaps	%%xmm1,0x10(%%rbx)	/* a(jp+p1 ) */\n\t"\
		"movaps	%%xmm3,0x10(%%rdx)	/* a(jp+p3 ) */\n\t"\
		"addpd	%%xmm0,%%xmm4		/* t0E+t2E */\n\t"\
		"addpd	%%xmm2,%%xmm7		/* t1E+t3F */\n\t"\
		"addpd	%%xmm1,%%xmm5		/* t0F+t2F */\n\t"\
		"addpd	%%xmm3,%%xmm6		/* t1F+t3E */\n\t"\
		"movaps	%%xmm4,    (%%rax)	/* a(jt+p0 ) */\n\t"\
		"movaps	%%xmm7,    (%%rdx)	/* a(jt+p3 ) */\n\t"\
		"movaps	%%xmm5,0x10(%%rax)	/* a(jp+p0 ) */\n\t"\
		"movaps	%%xmm6,0x10(%%rcx)	/* a(jp+p2 ) */\n\t"\
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

#endif	/* radix32_ditN_cy_dif1_gcc_h_included */

