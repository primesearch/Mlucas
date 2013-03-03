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
#ifndef radix24_ditN_cy_dif1_gcc_h_included
#define radix24_ditN_cy_dif1_gcc_h_included

	/* NOTE: Must use out-array for temp storage here, otherwise input-index permutations hose outputs */
	#define	SSE2_RADIX24_DIT_NOTWIDDLE(Xadd,Xp01,Xp02,Xp03,Xp04,Xp08,Xp16,Xout,Xisrt2,Xcc3)\
	{\
	__asm__ volatile (\
		"/* SSE2_RADIX8_DIT_0TWIDDLE(add0+p[0,1,3,2,7,6,5,4], s1p00r) */\n\t"\
			"movslq	%[__p04],%%rdi	\n\t"\
			"movq	%[__add],%%rax	/* Use eax as base address throughout */\n\t"\
			"shlq	$3,%%rdi		\n\t"\
			"movslq	%[__p01],%%rbx	\n\t"\
			"movslq	%[__p02],%%rcx	\n\t"\
			"movslq	%[__p03],%%rdx	\n\t"\
			"addq	%%rdi,%%rax		/* eax <- add0+p04 */\n\t"\
			"shlq	$3,%%rbx		/* Pointer offset for floating doubles */\n\t"\
			"shlq	$3,%%rcx		\n\t"\
			"shlq	$3,%%rdx		\n\t"\
			"addq	%%rax,%%rbx		/* eax <- add0+p05 */\n\t"\
			"addq	%%rax,%%rcx		/* ebx <- add0+p06 */\n\t"\
			"addq	%%rax,%%rdx		/* ecx <- add0+p07 */\n\t"\
			"movq	%[__out],%%rsi	/* s1p00r */\n\t"\
			"/* MSVC macro assumes add8+p[7,6,5,4] in eax,ebx,ecx,edx, but here get add0+p[4,5,6,7], so replace eax <-> edx and ebx <-> ecx: */\n\t"\
			"/* Do the p0,p4 combo: */\n\t"\
			"movaps	    (%%rdx),%%xmm0			\n\t"\
			"movaps	0x10(%%rdx),%%xmm1			\n\t"\
			"movaps	%%xmm0,%%xmm2				\n\t"\
			"movaps	%%xmm1,%%xmm3				\n\t"\
			"addpd	    (%%rcx),%%xmm2			\n\t"\
			"addpd	0x10(%%rcx),%%xmm3			\n\t"\
			"subpd	    (%%rcx),%%xmm0			\n\t"\
			"subpd	0x10(%%rcx),%%xmm1			\n\t"\
			"\n\t"\
			"movaps	    (%%rbx),%%xmm4			\n\t"\
			"movaps	0x10(%%rbx),%%xmm5			\n\t"\
			"movaps	%%xmm4,%%xmm6				\n\t"\
			"movaps	%%xmm5,%%xmm7				\n\t"\
			"addpd	    (%%rax),%%xmm6			\n\t"\
			"addpd	0x10(%%rax),%%xmm7			\n\t"\
			"subpd	    (%%rax),%%xmm4			\n\t"\
			"subpd	0x10(%%rax),%%xmm5			\n\t"\
			"/* Copy t6r,i into main-array slot add6 */\n\t"\
			"movaps	%%xmm6,0xc0(%%rsi)			\n\t"\
			"movaps	%%xmm7,0xd0(%%rsi)			\n\t"\
			"/* Copy t7r,i into main-array slot add7 */\n\t"\
			"movaps	%%xmm4,0xe0(%%rsi)			\n\t"\
			"movaps	%%xmm5,0xf0(%%rsi)			\n\t"\
			"\n\t"\
			"/* Move outputs t5r,i into a(j1,j2+p5), first doing the addsub and mul by ISRT2: */\n\t"\
			"/* Move outputs t7r,i into a(j1,j2+p7), first doing the addsub and mul by ISRT2: */\n\t"\
			"\n\t"\
			"addpd	%%xmm2,%%xmm6				\n\t"\
			"addpd	%%xmm3,%%xmm7				\n\t"\
			"subpd	0xc0(%%rsi),%%xmm2			\n\t"\
			"subpd	0xd0(%%rsi),%%xmm3			\n\t"\
			"/* Move t6r,i into a(j1,j2+p6) */	\n\t"\
			"movaps	%%xmm2,0xc0(%%rsi)			\n\t"\
			"movaps	%%xmm3,0xd0(%%rsi)			\n\t"\
			"\n\t"\
			"movaps	%%xmm4,%%xmm2				\n\t"\
			"movaps	%%xmm5,%%xmm3				\n\t"\
			"addpd	%%xmm0,%%xmm5				\n\t"\
			"subpd	%%xmm3,%%xmm0				\n\t"\
			"addpd	%%xmm1,%%xmm4				\n\t"\
			"subpd	%%xmm2,%%xmm1				\n\t"\
			"\n\t"\
			"movaps	%%xmm5,%%xmm2				\n\t"\
			"movaps	%%xmm1,%%xmm3				\n\t"\
			"addpd	%%xmm1,%%xmm5				\n\t"\
			"movaps	0x300(%%rsi),%%xmm1	/* isrt2 */	\n\t"\
			"subpd	%%xmm3,%%xmm2				\n\t"\
			"mulpd	%%xmm1,%%xmm5				\n\t"\
			"mulpd	%%xmm1,%%xmm2				\n\t"\
			"movaps	%%xmm0,%%xmm3				\n\t"\
			"\n\t"\
			"movaps	%%xmm5,0xa0(%%rsi)			\n\t"\
			"movaps	%%xmm4,%%xmm5				\n\t"\
			"addpd	%%xmm4,%%xmm0				\n\t"\
			"movaps	%%xmm2,0xb0(%%rsi)			\n\t"\
			"subpd	%%xmm5,%%xmm3				\n\t"\
			"mulpd	%%xmm1,%%xmm0				\n\t"\
			"mulpd	%%xmm1,%%xmm3				\n\t"\
			"movaps	%%xmm0,0xe0(%%rsi)			\n\t"\
			"movaps	%%xmm3,0xf0(%%rsi)			\n\t"\
			"\n\t"\
			"/* MSVC macro assumes add0+p[0,1,3,2] in eax,ebx,ecx,edx, but here get add0+p[0,1,2,3], so replace ecx <-> edx: */\n\t"\
			"subq	%%rdi,%%rax		/* add2 = add0     */\n\t"\
			"subq	%%rdi,%%rbx		/* add3 = add0+p01 */\n\t"\
			"subq	%%rdi,%%rcx		/* add1 = add0+p02 */\n\t"\
			"subq	%%rdi,%%rdx		/* add0 = add0+p03 */\n\t"\
			"\n\t"\
			"movaps	    (%%rax),%%xmm0			\n\t"\
			"movaps	0x10(%%rax),%%xmm1			\n\t"\
			"movaps	%%xmm0,%%xmm2				\n\t"\
			"movaps	%%xmm1,%%xmm3				\n\t"\
			"addpd	    (%%rbx),%%xmm2			\n\t"\
			"addpd	0x10(%%rbx),%%xmm3			\n\t"\
			"subpd	    (%%rbx),%%xmm0			\n\t"\
			"subpd	0x10(%%rbx),%%xmm1			\n\t"\
			"\n\t"\
			"/* Move   t4r,i into temp(j1+p0) in anticipation of final outputs (t0+t4)r,i which will go there: */\n\t"\
			"movaps	%%xmm6,    (%%rsi)			\n\t"\
			"movaps	%%xmm7,0x10(%%rsi)			\n\t"\
			"addpd	%%xmm6,%%xmm6				\n\t"\
			"addpd	%%xmm7,%%xmm7				\n\t"\
			"/* Move 2*t4r,i into temp(j1+p4) in anticipation of final outputs (t0-t4)r,i which will go there: */\n\t"\
			"movaps	%%xmm6,0x80(%%rsi)			\n\t"\
			"movaps	%%xmm7,0x90(%%rsi)			\n\t"\
			"\n\t"\
			"movaps	    (%%rdx),%%xmm4			\n\t"\
			"movaps	0x10(%%rdx),%%xmm5			\n\t"\
			"movaps	%%xmm4,%%xmm6				\n\t"\
			"movaps	%%xmm5,%%xmm7				\n\t"\
			"addpd	    (%%rcx),%%xmm6			\n\t"\
			"addpd	0x10(%%rcx),%%xmm7			\n\t"\
			"subpd	    (%%rcx),%%xmm4			\n\t"\
			"subpd	0x10(%%rcx),%%xmm5			\n\t"\
			"/* Copy t5,6 into temp(j1+p2) */	\n\t"\
			"movaps	%%xmm6,0x40(%%rsi)			\n\t"\
			"movaps	%%xmm7,0x50(%%rsi)			\n\t"\
			"addpd	%%xmm2,%%xmm6				\n\t"\
			"addpd	%%xmm3,%%xmm7				\n\t"\
			"subpd	0x40(%%rsi),%%xmm2			\n\t"\
			"subpd	0x50(%%rsi),%%xmm3			\n\t"\
			"\n\t"\
			"/* Compute and dump first 2 outputs now, in order to free up 2 registers: */\n\t"\
			"addpd	    (%%rsi),%%xmm6			\n\t"\
			"addpd	0x10(%%rsi),%%xmm7			\n\t"\
			"movaps	%%xmm6,    (%%rsi)			\n\t"\
			"movaps	%%xmm7,0x10(%%rsi)			\n\t"\
			"\n\t"\
			"subpd	0x80(%%rsi),%%xmm6			\n\t"\
			"subpd	0x90(%%rsi),%%xmm7			\n\t"\
			"movaps	%%xmm6,0x80(%%rsi)			\n\t"\
			"movaps	%%xmm7,0x90(%%rsi)			\n\t"\
			"\n\t"\
			"movaps	%%xmm4,%%xmm6				\n\t"\
			"movaps	%%xmm5,%%xmm7				\n\t"\
			"addpd	%%xmm0,%%xmm5				\n\t"\
			"subpd	%%xmm7,%%xmm0				\n\t"\
			"addpd	%%xmm1,%%xmm4				\n\t"\
			"subpd	%%xmm6,%%xmm1				\n\t"\
			"\n\t"\
			"movaps	%%xmm2,%%xmm6				\n\t"\
			"movaps	%%xmm3,%%xmm7				\n\t"\
			"addpd	0xd0(%%rsi),%%xmm2			\n\t"\
			"subpd	0xc0(%%rsi),%%xmm3			\n\t"\
			"subpd	0xd0(%%rsi),%%xmm6			\n\t"\
			"addpd	0xc0(%%rsi),%%xmm7			\n\t"\
			"movaps	%%xmm2,0xc0(%%rsi)			\n\t"\
			"movaps	%%xmm3,0xd0(%%rsi)			\n\t"\
			"movaps	%%xmm6,0x40(%%rsi)			\n\t"\
			"movaps	%%xmm7,0x50(%%rsi)			\n\t"\
			"\n\t"\
			"movaps	%%xmm5,%%xmm2				\n\t"\
			"movaps	%%xmm0,%%xmm6				\n\t"\
			"movaps	%%xmm1,%%xmm3				\n\t"\
			"movaps	%%xmm4,%%xmm7				\n\t"\
			"addpd	0xa0(%%rsi),%%xmm5			\n\t"\
			"subpd	0xf0(%%rsi),%%xmm0			\n\t"\
			"subpd	0xb0(%%rsi),%%xmm1			\n\t"\
			"subpd	0xe0(%%rsi),%%xmm4			\n\t"\
			"subpd	0xa0(%%rsi),%%xmm2			\n\t"\
			"addpd	0xf0(%%rsi),%%xmm6			\n\t"\
			"addpd	0xb0(%%rsi),%%xmm3			\n\t"\
			"addpd	0xe0(%%rsi),%%xmm7			\n\t"\
			"movaps	%%xmm5,0xe0(%%rsi)			\n\t"\
			"movaps	%%xmm0,0xa0(%%rsi)			\n\t"\
			"movaps	%%xmm1,0xf0(%%rsi)			\n\t"\
			"movaps	%%xmm4,0xb0(%%rsi)			\n\t"\
			"movaps	%%xmm2,0x60(%%rsi)			\n\t"\
			"movaps	%%xmm6,0x20(%%rsi)			\n\t"\
			"movaps	%%xmm3,0x70(%%rsi)			\n\t"\
			"movaps	%%xmm7,0x30(%%rsi)			\n\t"\
			"\n\t"\
		"/* SSE2_RADIX8_DIT_0TWIDDLE([add0+p08]+p[5,4,6,7,1,0,2,3], s1p08r) */\n\t"\
			"movslq	%[__p08],%%rdi	\n\t"\
			"shlq	$3,%%rdi		/* Pointer offset for floating doubles */\n\t"\
			"addq	%%rdi,%%rax		/* edx <- add8     */\n\t"\
			"addq	%%rdi,%%rbx		/* ecx <- add8+p01 */\n\t"\
			"addq	%%rdi,%%rcx		/* ebx <- add8+p02 */\n\t"\
			"addq	%%rdi,%%rdx		/* eax <- add8+p03 */\n\t"\
			"addq	$0x100,%%rsi	/* s1p08r */\n\t"\
			"/* MSVC macro assumes add8+p[0,1,2,3] in eax,ebx,ecx,edx, but here get add+p[1,0,2,3], so replace eax <-> ebx: */\n\t"\
			"/* Do the p0,p4 combo: */\n\t"\
			"movaps	    (%%rbx),%%xmm0			\n\t"\
			"movaps	0x10(%%rbx),%%xmm1			\n\t"\
			"movaps	%%xmm0,%%xmm2				\n\t"\
			"movaps	%%xmm1,%%xmm3				\n\t"\
			"addpd	    (%%rax),%%xmm2			\n\t"\
			"addpd	0x10(%%rax),%%xmm3			\n\t"\
			"subpd	    (%%rax),%%xmm0			\n\t"\
			"subpd	0x10(%%rax),%%xmm1			\n\t"\
			"\n\t"\
			"movaps	    (%%rcx),%%xmm4			\n\t"\
			"movaps	0x10(%%rcx),%%xmm5			\n\t"\
			"movaps	%%xmm4,%%xmm6				\n\t"\
			"movaps	%%xmm5,%%xmm7				\n\t"\
			"addpd	    (%%rdx),%%xmm6			\n\t"\
			"addpd	0x10(%%rdx),%%xmm7			\n\t"\
			"subpd	    (%%rdx),%%xmm4			\n\t"\
			"subpd	0x10(%%rdx),%%xmm5			\n\t"\
			"/* Copy t6r,i into main-array slot add6 */\n\t"\
			"movaps	%%xmm6,0xc0(%%rsi)			\n\t"\
			"movaps	%%xmm7,0xd0(%%rsi)			\n\t"\
			"/* Copy t7r,i into main-array slot add7 */\n\t"\
			"movaps	%%xmm4,0xe0(%%rsi)			\n\t"\
			"movaps	%%xmm5,0xf0(%%rsi)			\n\t"\
			"\n\t"\
			"/* Move outputs t5r,i into a(j1,j2+p5), first doing the addsub and mul by ISRT2: */\n\t"\
			"/* Move outputs t7r,i into a(j1,j2+p7), first doing the addsub and mul by ISRT2: */\n\t"\
			"\n\t"\
			"addpd	%%xmm2,%%xmm6				\n\t"\
			"addpd	%%xmm3,%%xmm7				\n\t"\
			"subpd	0xc0(%%rsi),%%xmm2			\n\t"\
			"subpd	0xd0(%%rsi),%%xmm3			\n\t"\
			"/* Move t6r,i into a(j1,j2+p6) */	\n\t"\
			"movaps	%%xmm2,0xc0(%%rsi)			\n\t"\
			"movaps	%%xmm3,0xd0(%%rsi)			\n\t"\
			"\n\t"\
			"movaps	%%xmm4,%%xmm2				\n\t"\
			"movaps	%%xmm5,%%xmm3				\n\t"\
			"addpd	%%xmm0,%%xmm5				\n\t"\
			"subpd	%%xmm3,%%xmm0				\n\t"\
			"addpd	%%xmm1,%%xmm4				\n\t"\
			"subpd	%%xmm2,%%xmm1				\n\t"\
			"\n\t"\
			"movaps	%%xmm5,%%xmm2				\n\t"\
			"movaps	%%xmm1,%%xmm3				\n\t"\
			"addpd	%%xmm1,%%xmm5				\n\t"\
			"movaps	0x200(%%rsi),%%xmm1	/* isrt2 */	\n\t"\
			"subpd	%%xmm3,%%xmm2				\n\t"\
			"mulpd	%%xmm1,%%xmm5				\n\t"\
			"mulpd	%%xmm1,%%xmm2				\n\t"\
			"movaps	%%xmm0,%%xmm3				\n\t"\
			"\n\t"\
			"movaps	%%xmm5,0xa0(%%rsi)			\n\t"\
			"movaps	%%xmm4,%%xmm5				\n\t"\
			"addpd	%%xmm4,%%xmm0				\n\t"\
			"movaps	%%xmm2,0xb0(%%rsi)			\n\t"\
			"subpd	%%xmm5,%%xmm3				\n\t"\
			"mulpd	%%xmm1,%%xmm0				\n\t"\
			"mulpd	%%xmm1,%%xmm3				\n\t"\
			"movaps	%%xmm0,0xe0(%%rsi)			\n\t"\
			"movaps	%%xmm3,0xf0(%%rsi)			\n\t"\
			"\n\t"\
			"/* MSVC macro assumes add8+p[5,4,6,7] in eax,ebx,ecx,edx, but here get add8+p[4,5,6,7], so replace eax <-> ebx: */\n\t"\
			"movslq	%[__p04],%%rdi	\n\t"\
			"shlq	$3,%%rdi		/* Pointer offset for floating doubles */\n\t"\
			"addq	%%rdi,%%rax		/* add0 = add8+p04 */\n\t"\
			"addq	%%rdi,%%rbx		/* add1 = add8+p05 */\n\t"\
			"addq	%%rdi,%%rcx		/* add3 = add8+p06 */\n\t"\
			"addq	%%rdi,%%rdx		/* add2 = add8+p07 */\n\t"\
			"\n\t"\
			"movaps	    (%%rbx),%%xmm0			\n\t"\
			"movaps	0x10(%%rbx),%%xmm1			\n\t"\
			"movaps	%%xmm0,%%xmm2				\n\t"\
			"movaps	%%xmm1,%%xmm3				\n\t"\
			"addpd	    (%%rax),%%xmm2			\n\t"\
			"addpd	0x10(%%rax),%%xmm3			\n\t"\
			"subpd	    (%%rax),%%xmm0			\n\t"\
			"subpd	0x10(%%rax),%%xmm1			\n\t"\
			"\n\t"\
			"/* Move   t4r,i into temp(j1+p0) in anticipation of final outputs (t0+t4)r,i which will go there: */\n\t"\
			"movaps	%%xmm6,    (%%rsi)			\n\t"\
			"movaps	%%xmm7,0x10(%%rsi)			\n\t"\
			"addpd	%%xmm6,%%xmm6				\n\t"\
			"addpd	%%xmm7,%%xmm7				\n\t"\
			"/* Move 2*t4r,i into temp(j1+p4) in anticipation of final outputs (t0-t4)r,i which will go there: */\n\t"\
			"movaps	%%xmm6,0x80(%%rsi)			\n\t"\
			"movaps	%%xmm7,0x90(%%rsi)			\n\t"\
			"\n\t"\
			"movaps	    (%%rcx),%%xmm4			\n\t"\
			"movaps	0x10(%%rcx),%%xmm5			\n\t"\
			"movaps	%%xmm4,%%xmm6				\n\t"\
			"movaps	%%xmm5,%%xmm7				\n\t"\
			"addpd	    (%%rdx),%%xmm6			\n\t"\
			"addpd	0x10(%%rdx),%%xmm7			\n\t"\
			"subpd	    (%%rdx),%%xmm4			\n\t"\
			"subpd	0x10(%%rdx),%%xmm5			\n\t"\
			"/* Copy t5,6 into temp(j1+p2) */	\n\t"\
			"movaps	%%xmm6,0x40(%%rsi)			\n\t"\
			"movaps	%%xmm7,0x50(%%rsi)			\n\t"\
			"addpd	%%xmm2,%%xmm6				\n\t"\
			"addpd	%%xmm3,%%xmm7				\n\t"\
			"subpd	0x40(%%rsi),%%xmm2			\n\t"\
			"subpd	0x50(%%rsi),%%xmm3			\n\t"\
			"\n\t"\
			"/* Compute and dump first 2 outputs now, in order to free up 2 registers: */\n\t"\
			"addpd	    (%%rsi),%%xmm6			\n\t"\
			"addpd	0x10(%%rsi),%%xmm7			\n\t"\
			"movaps	%%xmm6,    (%%rsi)			\n\t"\
			"movaps	%%xmm7,0x10(%%rsi)			\n\t"\
			"\n\t"\
			"subpd	0x80(%%rsi),%%xmm6			\n\t"\
			"subpd	0x90(%%rsi),%%xmm7			\n\t"\
			"movaps	%%xmm6,0x80(%%rsi)			\n\t"\
			"movaps	%%xmm7,0x90(%%rsi)			\n\t"\
			"\n\t"\
			"movaps	%%xmm4,%%xmm6				\n\t"\
			"movaps	%%xmm5,%%xmm7				\n\t"\
			"addpd	%%xmm0,%%xmm5				\n\t"\
			"subpd	%%xmm7,%%xmm0				\n\t"\
			"addpd	%%xmm1,%%xmm4				\n\t"\
			"subpd	%%xmm6,%%xmm1				\n\t"\
			"\n\t"\
			"movaps	%%xmm2,%%xmm6				\n\t"\
			"movaps	%%xmm3,%%xmm7				\n\t"\
			"addpd	0xd0(%%rsi),%%xmm2			\n\t"\
			"subpd	0xc0(%%rsi),%%xmm3			\n\t"\
			"subpd	0xd0(%%rsi),%%xmm6			\n\t"\
			"addpd	0xc0(%%rsi),%%xmm7			\n\t"\
			"movaps	%%xmm2,0xc0(%%rsi)			\n\t"\
			"movaps	%%xmm3,0xd0(%%rsi)			\n\t"\
			"movaps	%%xmm6,0x40(%%rsi)			\n\t"\
			"movaps	%%xmm7,0x50(%%rsi)			\n\t"\
			"\n\t"\
			"movaps	%%xmm5,%%xmm2				\n\t"\
			"movaps	%%xmm0,%%xmm6				\n\t"\
			"movaps	%%xmm1,%%xmm3				\n\t"\
			"movaps	%%xmm4,%%xmm7				\n\t"\
			"addpd	0xa0(%%rsi),%%xmm5			\n\t"\
			"subpd	0xf0(%%rsi),%%xmm0			\n\t"\
			"subpd	0xb0(%%rsi),%%xmm1			\n\t"\
			"subpd	0xe0(%%rsi),%%xmm4			\n\t"\
			"subpd	0xa0(%%rsi),%%xmm2			\n\t"\
			"addpd	0xf0(%%rsi),%%xmm6			\n\t"\
			"addpd	0xb0(%%rsi),%%xmm3			\n\t"\
			"addpd	0xe0(%%rsi),%%xmm7			\n\t"\
			"movaps	%%xmm5,0xe0(%%rsi)			\n\t"\
			"movaps	%%xmm0,0xa0(%%rsi)			\n\t"\
			"movaps	%%xmm1,0xf0(%%rsi)			\n\t"\
			"movaps	%%xmm4,0xb0(%%rsi)			\n\t"\
			"movaps	%%xmm2,0x60(%%rsi)			\n\t"\
			"movaps	%%xmm6,0x20(%%rsi)			\n\t"\
			"movaps	%%xmm3,0x70(%%rsi)			\n\t"\
			"movaps	%%xmm7,0x30(%%rsi)			\n\t"\
			"\n\t"\
		"/* SSE2_RADIX8_DIT_0TWIDDLE([add0+p16]+p[2,3,0,1,4,5,7,6], s1p16r) */\n\t"\
			"movslq	%[__p08],%%rdi	\n\t"\
			"shlq	$3,%%rdi		/* Pointer offset for floating doubles */\n\t"\
			"subq	%%rdi,%%rax		/* edx <- add0 +p04 */\n\t"\
			"subq	%%rdi,%%rbx		/* ecx <- add0 +p05 */\n\t"\
			"subq	%%rdi,%%rcx		/* ebx <- add0 +p06 */\n\t"\
			"subq	%%rdi,%%rdx		/* eax <- add0 +p07 */\n\t"\
			"movslq	%[__p16],%%rdi	\n\t"\
			"shlq	$3,%%rdi		/* Pointer offset for floating doubles */\n\t"\
			"addq	%%rdi,%%rax		/* edx <- add16+p04 */\n\t"\
			"addq	%%rdi,%%rbx		/* ecx <- add16+p05 */\n\t"\
			"addq	%%rdi,%%rcx		/* ebx <- add16+p06 */\n\t"\
			"addq	%%rdi,%%rdx		/* eax <- add16+p07 */\n\t"\
			"addq	$0x100,%%rsi	/* s1p16r */\n\t"\
			"/* MSVC macro assumes add16+p[4,5,7,6] in eax,ebx,ecx,edx, but here get add16+p[4,5,6,7], so replace ecx <-> edx: */\n\t"\
			"/* Do the p0,p4 combo: */\n\t"\
			"movaps	    (%%rax),%%xmm0			\n\t"\
			"movaps	0x10(%%rax),%%xmm1			\n\t"\
			"movaps	%%xmm0,%%xmm2				\n\t"\
			"movaps	%%xmm1,%%xmm3				\n\t"\
			"addpd	    (%%rbx),%%xmm2			\n\t"\
			"addpd	0x10(%%rbx),%%xmm3			\n\t"\
			"subpd	    (%%rbx),%%xmm0			\n\t"\
			"subpd	0x10(%%rbx),%%xmm1			\n\t"\
			"\n\t"\
			"movaps	    (%%rdx),%%xmm4			\n\t"\
			"movaps	0x10(%%rdx),%%xmm5			\n\t"\
			"movaps	%%xmm4,%%xmm6				\n\t"\
			"movaps	%%xmm5,%%xmm7				\n\t"\
			"addpd	    (%%rcx),%%xmm6			\n\t"\
			"addpd	0x10(%%rcx),%%xmm7			\n\t"\
			"subpd	    (%%rcx),%%xmm4			\n\t"\
			"subpd	0x10(%%rcx),%%xmm5			\n\t"\
			"/* Copy t6r,i into main-array slot add6 */\n\t"\
			"movaps	%%xmm6,0xc0(%%rsi)			\n\t"\
			"movaps	%%xmm7,0xd0(%%rsi)			\n\t"\
			"/* Copy t7r,i into main-array slot add7 */\n\t"\
			"movaps	%%xmm4,0xe0(%%rsi)			\n\t"\
			"movaps	%%xmm5,0xf0(%%rsi)			\n\t"\
			"\n\t"\
			"/* Move outputs t5r,i into a(j1,j2+p5), first doing the addsub and mul by ISRT2: */\n\t"\
			"/* Move outputs t7r,i into a(j1,j2+p7), first doing the addsub and mul by ISRT2: */\n\t"\
			"\n\t"\
			"addpd	%%xmm2,%%xmm6				\n\t"\
			"addpd	%%xmm3,%%xmm7				\n\t"\
			"subpd	0xc0(%%rsi),%%xmm2			\n\t"\
			"subpd	0xd0(%%rsi),%%xmm3			\n\t"\
			"/* Move t6r,i into a(j1,j2+p6) */	\n\t"\
			"movaps	%%xmm2,0xc0(%%rsi)			\n\t"\
			"movaps	%%xmm3,0xd0(%%rsi)			\n\t"\
			"\n\t"\
			"movaps	%%xmm4,%%xmm2				\n\t"\
			"movaps	%%xmm5,%%xmm3				\n\t"\
			"addpd	%%xmm0,%%xmm5				\n\t"\
			"subpd	%%xmm3,%%xmm0				\n\t"\
			"addpd	%%xmm1,%%xmm4				\n\t"\
			"subpd	%%xmm2,%%xmm1				\n\t"\
			"\n\t"\
			"movaps	%%xmm5,%%xmm2				\n\t"\
			"movaps	%%xmm1,%%xmm3				\n\t"\
			"addpd	%%xmm1,%%xmm5				\n\t"\
			"movaps	0x100(%%rsi),%%xmm1	/* isrt2 */	\n\t"\
			"subpd	%%xmm3,%%xmm2				\n\t"\
			"mulpd	%%xmm1,%%xmm5				\n\t"\
			"mulpd	%%xmm1,%%xmm2				\n\t"\
			"movaps	%%xmm0,%%xmm3				\n\t"\
			"\n\t"\
			"movaps	%%xmm5,0xa0(%%rsi)			\n\t"\
			"movaps	%%xmm4,%%xmm5				\n\t"\
			"addpd	%%xmm4,%%xmm0				\n\t"\
			"movaps	%%xmm2,0xb0(%%rsi)			\n\t"\
			"subpd	%%xmm5,%%xmm3				\n\t"\
			"mulpd	%%xmm1,%%xmm0				\n\t"\
			"mulpd	%%xmm1,%%xmm3				\n\t"\
			"movaps	%%xmm0,0xe0(%%rsi)			\n\t"\
			"movaps	%%xmm3,0xf0(%%rsi)			\n\t"\
			"\n\t"\
			"/* MSVC macro assumes add16+p[2,3,0,1] in eax,ebx,ecx,edx, but here get add16+p[0,1,2,3], so replace eax <-> ecx, ebx <-> edx: */\n\t"\
			"movslq	%[__p04],%%rdi	\n\t"\
			"shlq	$3,%%rdi		/* Pointer offset for floating doubles */\n\t"\
			"subq	%%rdi,%%rax		/* add0 = add8+p04 */\n\t"\
			"subq	%%rdi,%%rbx		/* add1 = add8+p05 */\n\t"\
			"subq	%%rdi,%%rcx		/* add3 = add8+p06 */\n\t"\
			"subq	%%rdi,%%rdx		/* add2 = add8+p07 */\n\t"\
			"\n\t"\
			"movaps	    (%%rcx),%%xmm0			\n\t"\
			"movaps	0x10(%%rcx),%%xmm1			\n\t"\
			"movaps	%%xmm0,%%xmm2				\n\t"\
			"movaps	%%xmm1,%%xmm3				\n\t"\
			"addpd	    (%%rdx),%%xmm2			\n\t"\
			"addpd	0x10(%%rdx),%%xmm3			\n\t"\
			"subpd	    (%%rdx),%%xmm0			\n\t"\
			"subpd	0x10(%%rdx),%%xmm1			\n\t"\
			"\n\t"\
			"/* Move   t4r,i into temp(j1+p0) in anticipation of final outputs (t0+t4)r,i which will go there: */\n\t"\
			"movaps	%%xmm6,    (%%rsi)			\n\t"\
			"movaps	%%xmm7,0x10(%%rsi)			\n\t"\
			"addpd	%%xmm6,%%xmm6				\n\t"\
			"addpd	%%xmm7,%%xmm7				\n\t"\
			"/* Move 2*t4r,i into temp(j1+p4) in anticipation of final outputs (t0-t4)r,i which will go there: */\n\t"\
			"movaps	%%xmm6,0x80(%%rsi)			\n\t"\
			"movaps	%%xmm7,0x90(%%rsi)			\n\t"\
			"\n\t"\
			"movaps	    (%%rax),%%xmm4			\n\t"\
			"movaps	0x10(%%rax),%%xmm5			\n\t"\
			"movaps	%%xmm4,%%xmm6				\n\t"\
			"movaps	%%xmm5,%%xmm7				\n\t"\
			"addpd	    (%%rbx),%%xmm6			\n\t"\
			"addpd	0x10(%%rbx),%%xmm7			\n\t"\
			"subpd	    (%%rbx),%%xmm4			\n\t"\
			"subpd	0x10(%%rbx),%%xmm5			\n\t"\
			"/* Copy t5,6 into temp(j1+p2) */	\n\t"\
			"movaps	%%xmm6,0x40(%%rsi)			\n\t"\
			"movaps	%%xmm7,0x50(%%rsi)			\n\t"\
			"addpd	%%xmm2,%%xmm6				\n\t"\
			"addpd	%%xmm3,%%xmm7				\n\t"\
			"subpd	0x40(%%rsi),%%xmm2			\n\t"\
			"subpd	0x50(%%rsi),%%xmm3			\n\t"\
			"\n\t"\
			"/* Compute and dump first 2 outputs now, in order to free up 2 registers: */\n\t"\
			"addpd	    (%%rsi),%%xmm6			\n\t"\
			"addpd	0x10(%%rsi),%%xmm7			\n\t"\
			"movaps	%%xmm6,    (%%rsi)			\n\t"\
			"movaps	%%xmm7,0x10(%%rsi)			\n\t"\
			"\n\t"\
			"subpd	0x80(%%rsi),%%xmm6			\n\t"\
			"subpd	0x90(%%rsi),%%xmm7			\n\t"\
			"movaps	%%xmm6,0x80(%%rsi)			\n\t"\
			"movaps	%%xmm7,0x90(%%rsi)			\n\t"\
			"\n\t"\
			"movaps	%%xmm4,%%xmm6				\n\t"\
			"movaps	%%xmm5,%%xmm7				\n\t"\
			"addpd	%%xmm0,%%xmm5				\n\t"\
			"subpd	%%xmm7,%%xmm0				\n\t"\
			"addpd	%%xmm1,%%xmm4				\n\t"\
			"subpd	%%xmm6,%%xmm1				\n\t"\
			"\n\t"\
			"movaps	%%xmm2,%%xmm6				\n\t"\
			"movaps	%%xmm3,%%xmm7				\n\t"\
			"addpd	0xd0(%%rsi),%%xmm2			\n\t"\
			"subpd	0xc0(%%rsi),%%xmm3			\n\t"\
			"subpd	0xd0(%%rsi),%%xmm6			\n\t"\
			"addpd	0xc0(%%rsi),%%xmm7			\n\t"\
			"movaps	%%xmm2,0xc0(%%rsi)			\n\t"\
			"movaps	%%xmm3,0xd0(%%rsi)			\n\t"\
			"movaps	%%xmm6,0x40(%%rsi)			\n\t"\
			"movaps	%%xmm7,0x50(%%rsi)			\n\t"\
			"\n\t"\
			"movaps	%%xmm5,%%xmm2				\n\t"\
			"movaps	%%xmm0,%%xmm6				\n\t"\
			"movaps	%%xmm1,%%xmm3				\n\t"\
			"movaps	%%xmm4,%%xmm7				\n\t"\
			"addpd	0xa0(%%rsi),%%xmm5			\n\t"\
			"subpd	0xf0(%%rsi),%%xmm0			\n\t"\
			"subpd	0xb0(%%rsi),%%xmm1			\n\t"\
			"subpd	0xe0(%%rsi),%%xmm4			\n\t"\
			"subpd	0xa0(%%rsi),%%xmm2			\n\t"\
			"addpd	0xf0(%%rsi),%%xmm6			\n\t"\
			"addpd	0xb0(%%rsi),%%xmm3			\n\t"\
			"addpd	0xe0(%%rsi),%%xmm7			\n\t"\
			"movaps	%%xmm5,0xe0(%%rsi)			\n\t"\
			"movaps	%%xmm0,0xa0(%%rsi)			\n\t"\
			"movaps	%%xmm1,0xf0(%%rsi)			\n\t"\
			"movaps	%%xmm4,0xb0(%%rsi)			\n\t"\
			"movaps	%%xmm2,0x60(%%rsi)			\n\t"\
			"movaps	%%xmm6,0x20(%%rsi)			\n\t"\
			"movaps	%%xmm3,0x70(%%rsi)			\n\t"\
			"movaps	%%xmm7,0x30(%%rsi)			\n\t"\
			"\n\t"\
			"movq	%[__cc3],%%rdx				\n\t"\
			"movq	 %%rsi,%%rax	/* cpy s1p16r */\n\t"\
			"movq	 %%rsi,%%rbx	/* cpy s1p16r */\n\t"\
			"movq	 %%rsi,%%rcx	/* cpy s1p16r */\n\t"\
			"\n\t"\
		"/* SSE2_RADIX_03_DFT(s1p00r,s1p08r,s1p16r,cc3,s1p00r,s1p16r,s1p08r) */\n\t"\
			"subq	$0x200,%%rax	/* s1p00r */\n\t"\
			"subq	$0x100,%%rbx	/* s1p08r */\n\t"\
			"\n\t"\
			"movaps	    (%%rbx),%%xmm2			\n\t"\
			"movaps	0x10(%%rbx),%%xmm3			\n\t"\
			"movaps	    (%%rax),%%xmm0			\n\t"\
			"movaps	0x10(%%rax),%%xmm1			\n\t"\
			"movaps	    (%%rcx),%%xmm6			\n\t"\
			"movaps	0x10(%%rcx),%%xmm7			\n\t"\
			"movaps	%%xmm2,%%xmm4				\n\t"\
			"movaps	%%xmm3,%%xmm5				\n\t"\
			"\n\t"\
			"addpd	%%xmm6,%%xmm2				\n\t"\
			"addpd	%%xmm7,%%xmm3				\n\t"\
			"subpd	%%xmm6,%%xmm4				\n\t"\
			"subpd	%%xmm7,%%xmm5				\n\t"\
			"addpd	%%xmm2,%%xmm0				\n\t"\
			"addpd	%%xmm3,%%xmm1				\n\t"\
			"movaps	    (%%rdx),%%xmm6			\n\t"\
			"movaps	0x10(%%rdx),%%xmm7			\n\t"\
			"movaps	%%xmm0,     (%%rax)			\n\t"\
			"movaps	%%xmm1,0x010(%%rax)			\n\t"\
			"\n\t"\
			"mulpd	%%xmm6,%%xmm2				\n\t"\
			"mulpd	%%xmm6,%%xmm3				\n\t"\
			"mulpd	%%xmm7,%%xmm4				\n\t"\
			"mulpd	%%xmm7,%%xmm5				\n\t"\
			"addpd	%%xmm0,%%xmm2				\n\t"\
			"addpd	%%xmm1,%%xmm3				\n\t"\
			"\n\t"\
			"movaps	%%xmm2,%%xmm0				\n\t"\
			"movaps	%%xmm3,%%xmm1				\n\t"\
			"\n\t"\
			"subpd	%%xmm5,%%xmm2				\n\t"\
			"addpd	%%xmm4,%%xmm3				\n\t"\
			"addpd	%%xmm5,%%xmm0				\n\t"\
			"subpd	%%xmm4,%%xmm1				\n\t"\
			"\n\t"\
			"movaps	%%xmm2,     (%%rcx)	/*b<>c*/\n\t"\
			"movaps	%%xmm3,0x010(%%rcx)			\n\t"\
			"movaps	%%xmm0,     (%%rbx)			\n\t"\
			"movaps	%%xmm1,0x010(%%rbx)			\n\t"\
			"\n\t"\
		"/* SSE2_RADIX_03_DFT(s1p01r,s1p09r,s1p17r,cc3,s1p09r,s1p01r,s1p17r) */\n\t"\
			"addq	$0x020,%%rax	/* s1p01r */\n\t"\
			"addq	$0x020,%%rbx	/* s1p09r */\n\t"\
			"addq	$0x020,%%rcx	/* s1p17r */\n\t"\
			"\n\t"\
			"movaps	    (%%rbx),%%xmm2			\n\t"\
			"movaps	0x10(%%rbx),%%xmm3			\n\t"\
			"movaps	    (%%rax),%%xmm0			\n\t"\
			"movaps	0x10(%%rax),%%xmm1			\n\t"\
			"movaps	    (%%rcx),%%xmm6			\n\t"\
			"movaps	0x10(%%rcx),%%xmm7			\n\t"\
			"movaps	%%xmm2,%%xmm4				\n\t"\
			"movaps	%%xmm3,%%xmm5				\n\t"\
			"\n\t"\
			"addpd	%%xmm6,%%xmm2				\n\t"\
			"addpd	%%xmm7,%%xmm3				\n\t"\
			"subpd	%%xmm6,%%xmm4				\n\t"\
			"subpd	%%xmm7,%%xmm5				\n\t"\
			"addpd	%%xmm2,%%xmm0				\n\t"\
			"addpd	%%xmm3,%%xmm1				\n\t"\
			"movaps	    (%%rdx),%%xmm6			\n\t"\
			"movaps	0x10(%%rdx),%%xmm7			\n\t"\
			"movaps	%%xmm0,     (%%rbx)	/*a<>b*/\n\t"\
			"movaps	%%xmm1,0x010(%%rbx)			\n\t"\
			"\n\t"\
			"mulpd	%%xmm6,%%xmm2				\n\t"\
			"mulpd	%%xmm6,%%xmm3				\n\t"\
			"mulpd	%%xmm7,%%xmm4				\n\t"\
			"mulpd	%%xmm7,%%xmm5				\n\t"\
			"addpd	%%xmm0,%%xmm2				\n\t"\
			"addpd	%%xmm1,%%xmm3				\n\t"\
			"\n\t"\
			"movaps	%%xmm2,%%xmm0				\n\t"\
			"movaps	%%xmm3,%%xmm1				\n\t"\
			"\n\t"\
			"subpd	%%xmm5,%%xmm2				\n\t"\
			"addpd	%%xmm4,%%xmm3				\n\t"\
			"addpd	%%xmm5,%%xmm0				\n\t"\
			"subpd	%%xmm4,%%xmm1				\n\t"\
			"\n\t"\
			"movaps	%%xmm2,     (%%rax)			\n\t"\
			"movaps	%%xmm3,0x010(%%rax)			\n\t"\
			"movaps	%%xmm0,     (%%rcx)			\n\t"\
			"movaps	%%xmm1,0x010(%%rcx)			\n\t"\
			"\n\t"\
		"/* SSE2_RADIX_03_DFT(s1p02r,s1p10r,s1p18r,cc3,s1p18r,s1p10r,s1p02r) */\n\t"\
			"addq	$0x020,%%rax	\n\t"\
			"addq	$0x020,%%rbx	\n\t"\
			"addq	$0x020,%%rcx	\n\t"\
			"\n\t"\
			"movaps	    (%%rbx),%%xmm2			\n\t"\
			"movaps	0x10(%%rbx),%%xmm3			\n\t"\
			"movaps	    (%%rax),%%xmm0			\n\t"\
			"movaps	0x10(%%rax),%%xmm1			\n\t"\
			"movaps	    (%%rcx),%%xmm6			\n\t"\
			"movaps	0x10(%%rcx),%%xmm7			\n\t"\
			"movaps	%%xmm2,%%xmm4				\n\t"\
			"movaps	%%xmm3,%%xmm5				\n\t"\
			"\n\t"\
			"addpd	%%xmm6,%%xmm2				\n\t"\
			"addpd	%%xmm7,%%xmm3				\n\t"\
			"subpd	%%xmm6,%%xmm4				\n\t"\
			"subpd	%%xmm7,%%xmm5				\n\t"\
			"addpd	%%xmm2,%%xmm0				\n\t"\
			"addpd	%%xmm3,%%xmm1				\n\t"\
			"movaps	    (%%rdx),%%xmm6			\n\t"\
			"movaps	0x10(%%rdx),%%xmm7			\n\t"\
			"movaps	%%xmm0,     (%%rcx)	/*a<>c*/\n\t"\
			"movaps	%%xmm1,0x010(%%rcx)			\n\t"\
			"\n\t"\
			"mulpd	%%xmm6,%%xmm2				\n\t"\
			"mulpd	%%xmm6,%%xmm3				\n\t"\
			"mulpd	%%xmm7,%%xmm4				\n\t"\
			"mulpd	%%xmm7,%%xmm5				\n\t"\
			"addpd	%%xmm0,%%xmm2				\n\t"\
			"addpd	%%xmm1,%%xmm3				\n\t"\
			"\n\t"\
			"movaps	%%xmm2,%%xmm0				\n\t"\
			"movaps	%%xmm3,%%xmm1				\n\t"\
			"\n\t"\
			"subpd	%%xmm5,%%xmm2				\n\t"\
			"addpd	%%xmm4,%%xmm3				\n\t"\
			"addpd	%%xmm5,%%xmm0				\n\t"\
			"subpd	%%xmm4,%%xmm1				\n\t"\
			"\n\t"\
			"movaps	%%xmm2,     (%%rbx)			\n\t"\
			"movaps	%%xmm3,0x010(%%rbx)			\n\t"\
			"movaps	%%xmm0,     (%%rax)			\n\t"\
			"movaps	%%xmm1,0x010(%%rax)			\n\t"\
			"\n\t"\
		"/* SSE2_RADIX_03_DFT(s1p03r,s1p11r,s1p19r,cc3,s1p03r,s1p19r,s1p11r) */\n\t"\
			"addq	$0x020,%%rax	\n\t"\
			"addq	$0x020,%%rbx	\n\t"\
			"addq	$0x020,%%rcx	\n\t"\
			"\n\t"\
			"movaps	    (%%rbx),%%xmm2			\n\t"\
			"movaps	0x10(%%rbx),%%xmm3			\n\t"\
			"movaps	    (%%rax),%%xmm0			\n\t"\
			"movaps	0x10(%%rax),%%xmm1			\n\t"\
			"movaps	    (%%rcx),%%xmm6			\n\t"\
			"movaps	0x10(%%rcx),%%xmm7			\n\t"\
			"movaps	%%xmm2,%%xmm4				\n\t"\
			"movaps	%%xmm3,%%xmm5				\n\t"\
			"\n\t"\
			"addpd	%%xmm6,%%xmm2				\n\t"\
			"addpd	%%xmm7,%%xmm3				\n\t"\
			"subpd	%%xmm6,%%xmm4				\n\t"\
			"subpd	%%xmm7,%%xmm5				\n\t"\
			"addpd	%%xmm2,%%xmm0				\n\t"\
			"addpd	%%xmm3,%%xmm1				\n\t"\
			"movaps	    (%%rdx),%%xmm6			\n\t"\
			"movaps	0x10(%%rdx),%%xmm7			\n\t"\
			"movaps	%%xmm0,     (%%rax)	/*b<>c*/\n\t"\
			"movaps	%%xmm1,0x010(%%rax)			\n\t"\
			"\n\t"\
			"mulpd	%%xmm6,%%xmm2				\n\t"\
			"mulpd	%%xmm6,%%xmm3				\n\t"\
			"mulpd	%%xmm7,%%xmm4				\n\t"\
			"mulpd	%%xmm7,%%xmm5				\n\t"\
			"addpd	%%xmm0,%%xmm2				\n\t"\
			"addpd	%%xmm1,%%xmm3				\n\t"\
			"\n\t"\
			"movaps	%%xmm2,%%xmm0				\n\t"\
			"movaps	%%xmm3,%%xmm1				\n\t"\
			"\n\t"\
			"subpd	%%xmm5,%%xmm2				\n\t"\
			"addpd	%%xmm4,%%xmm3				\n\t"\
			"addpd	%%xmm5,%%xmm0				\n\t"\
			"subpd	%%xmm4,%%xmm1				\n\t"\
			"\n\t"\
			"movaps	%%xmm2,     (%%rcx)			\n\t"\
			"movaps	%%xmm3,0x010(%%rcx)			\n\t"\
			"movaps	%%xmm0,     (%%rbx)			\n\t"\
			"movaps	%%xmm1,0x010(%%rbx)			\n\t"\
			"\n\t"\
		"/* SSE2_RADIX_03_DFT(s1p04r,s1p12r,s1p20r,cc3,s1p12r,s1p04r,s1p20r) */\n\t"\
			"addq	$0x020,%%rax	\n\t"\
			"addq	$0x020,%%rbx	\n\t"\
			"addq	$0x020,%%rcx	\n\t"\
			"\n\t"\
			"movaps	    (%%rbx),%%xmm2			\n\t"\
			"movaps	0x10(%%rbx),%%xmm3			\n\t"\
			"movaps	    (%%rax),%%xmm0			\n\t"\
			"movaps	0x10(%%rax),%%xmm1			\n\t"\
			"movaps	    (%%rcx),%%xmm6			\n\t"\
			"movaps	0x10(%%rcx),%%xmm7			\n\t"\
			"movaps	%%xmm2,%%xmm4				\n\t"\
			"movaps	%%xmm3,%%xmm5				\n\t"\
			"\n\t"\
			"addpd	%%xmm6,%%xmm2				\n\t"\
			"addpd	%%xmm7,%%xmm3				\n\t"\
			"subpd	%%xmm6,%%xmm4				\n\t"\
			"subpd	%%xmm7,%%xmm5				\n\t"\
			"addpd	%%xmm2,%%xmm0				\n\t"\
			"addpd	%%xmm3,%%xmm1				\n\t"\
			"movaps	    (%%rdx),%%xmm6			\n\t"\
			"movaps	0x10(%%rdx),%%xmm7			\n\t"\
			"movaps	%%xmm0,     (%%rbx)	/*a<>b*/\n\t"\
			"movaps	%%xmm1,0x010(%%rbx)			\n\t"\
			"\n\t"\
			"mulpd	%%xmm6,%%xmm2				\n\t"\
			"mulpd	%%xmm6,%%xmm3				\n\t"\
			"mulpd	%%xmm7,%%xmm4				\n\t"\
			"mulpd	%%xmm7,%%xmm5				\n\t"\
			"addpd	%%xmm0,%%xmm2				\n\t"\
			"addpd	%%xmm1,%%xmm3				\n\t"\
			"\n\t"\
			"movaps	%%xmm2,%%xmm0				\n\t"\
			"movaps	%%xmm3,%%xmm1				\n\t"\
			"\n\t"\
			"subpd	%%xmm5,%%xmm2				\n\t"\
			"addpd	%%xmm4,%%xmm3				\n\t"\
			"addpd	%%xmm5,%%xmm0				\n\t"\
			"subpd	%%xmm4,%%xmm1				\n\t"\
			"\n\t"\
			"movaps	%%xmm2,     (%%rax)			\n\t"\
			"movaps	%%xmm3,0x010(%%rax)			\n\t"\
			"movaps	%%xmm0,     (%%rcx)			\n\t"\
			"movaps	%%xmm1,0x010(%%rcx)			\n\t"\
			"\n\t"\
		"/* SSE2_RADIX_03_DFT(s1p05r,s1p13r,s1p21r,cc3,s1p21r,s1p13r,s1p05r) */\n\t"\
			"addq	$0x020,%%rax	\n\t"\
			"addq	$0x020,%%rbx	\n\t"\
			"addq	$0x020,%%rcx	\n\t"\
			"\n\t"\
			"movaps	    (%%rbx),%%xmm2			\n\t"\
			"movaps	0x10(%%rbx),%%xmm3			\n\t"\
			"movaps	    (%%rax),%%xmm0			\n\t"\
			"movaps	0x10(%%rax),%%xmm1			\n\t"\
			"movaps	    (%%rcx),%%xmm6			\n\t"\
			"movaps	0x10(%%rcx),%%xmm7			\n\t"\
			"movaps	%%xmm2,%%xmm4				\n\t"\
			"movaps	%%xmm3,%%xmm5				\n\t"\
			"\n\t"\
			"addpd	%%xmm6,%%xmm2				\n\t"\
			"addpd	%%xmm7,%%xmm3				\n\t"\
			"subpd	%%xmm6,%%xmm4				\n\t"\
			"subpd	%%xmm7,%%xmm5				\n\t"\
			"addpd	%%xmm2,%%xmm0				\n\t"\
			"addpd	%%xmm3,%%xmm1				\n\t"\
			"movaps	    (%%rdx),%%xmm6			\n\t"\
			"movaps	0x10(%%rdx),%%xmm7			\n\t"\
			"movaps	%%xmm0,     (%%rcx)	/*a<>c*/\n\t"\
			"movaps	%%xmm1,0x010(%%rcx)			\n\t"\
			"\n\t"\
			"mulpd	%%xmm6,%%xmm2				\n\t"\
			"mulpd	%%xmm6,%%xmm3				\n\t"\
			"mulpd	%%xmm7,%%xmm4				\n\t"\
			"mulpd	%%xmm7,%%xmm5				\n\t"\
			"addpd	%%xmm0,%%xmm2				\n\t"\
			"addpd	%%xmm1,%%xmm3				\n\t"\
			"\n\t"\
			"movaps	%%xmm2,%%xmm0				\n\t"\
			"movaps	%%xmm3,%%xmm1				\n\t"\
			"\n\t"\
			"subpd	%%xmm5,%%xmm2				\n\t"\
			"addpd	%%xmm4,%%xmm3				\n\t"\
			"addpd	%%xmm5,%%xmm0				\n\t"\
			"subpd	%%xmm4,%%xmm1				\n\t"\
			"\n\t"\
			"movaps	%%xmm2,     (%%rbx)			\n\t"\
			"movaps	%%xmm3,0x010(%%rbx)			\n\t"\
			"movaps	%%xmm0,     (%%rax)			\n\t"\
			"movaps	%%xmm1,0x010(%%rax)			\n\t"\
			"\n\t"\
		"/* SSE2_RADIX_03_DFT(s1p06r,s1p14r,s1p22r,cc3,s1p06r,s1p22r,s1p14r) */\n\t"\
			"addq	$0x020,%%rax	\n\t"\
			"addq	$0x020,%%rbx	\n\t"\
			"addq	$0x020,%%rcx	\n\t"\
			"\n\t"\
			"movaps	    (%%rbx),%%xmm2			\n\t"\
			"movaps	0x10(%%rbx),%%xmm3			\n\t"\
			"movaps	    (%%rax),%%xmm0			\n\t"\
			"movaps	0x10(%%rax),%%xmm1			\n\t"\
			"movaps	    (%%rcx),%%xmm6			\n\t"\
			"movaps	0x10(%%rcx),%%xmm7			\n\t"\
			"movaps	%%xmm2,%%xmm4				\n\t"\
			"movaps	%%xmm3,%%xmm5				\n\t"\
			"\n\t"\
			"addpd	%%xmm6,%%xmm2				\n\t"\
			"addpd	%%xmm7,%%xmm3				\n\t"\
			"subpd	%%xmm6,%%xmm4				\n\t"\
			"subpd	%%xmm7,%%xmm5				\n\t"\
			"addpd	%%xmm2,%%xmm0				\n\t"\
			"addpd	%%xmm3,%%xmm1				\n\t"\
			"movaps	    (%%rdx),%%xmm6			\n\t"\
			"movaps	0x10(%%rdx),%%xmm7			\n\t"\
			"movaps	%%xmm0,     (%%rax)	/*b<>c*/\n\t"\
			"movaps	%%xmm1,0x010(%%rax)			\n\t"\
			"\n\t"\
			"mulpd	%%xmm6,%%xmm2				\n\t"\
			"mulpd	%%xmm6,%%xmm3				\n\t"\
			"mulpd	%%xmm7,%%xmm4				\n\t"\
			"mulpd	%%xmm7,%%xmm5				\n\t"\
			"addpd	%%xmm0,%%xmm2				\n\t"\
			"addpd	%%xmm1,%%xmm3				\n\t"\
			"\n\t"\
			"movaps	%%xmm2,%%xmm0				\n\t"\
			"movaps	%%xmm3,%%xmm1				\n\t"\
			"\n\t"\
			"subpd	%%xmm5,%%xmm2				\n\t"\
			"addpd	%%xmm4,%%xmm3				\n\t"\
			"addpd	%%xmm5,%%xmm0				\n\t"\
			"subpd	%%xmm4,%%xmm1				\n\t"\
			"\n\t"\
			"movaps	%%xmm2,     (%%rcx)			\n\t"\
			"movaps	%%xmm3,0x010(%%rcx)			\n\t"\
			"movaps	%%xmm0,     (%%rbx)			\n\t"\
			"movaps	%%xmm1,0x010(%%rbx)			\n\t"\
			"\n\t"\
		"/* SSE2_RADIX_03_DFT(s1p07r,s1p15r,s1p23r,cc3,s1p15r,s1p07r,s1p23r) */\n\t"\
			"addq	$0x020,%%rax	\n\t"\
			"addq	$0x020,%%rbx	\n\t"\
			"addq	$0x020,%%rcx	\n\t"\
			"\n\t"\
			"movaps	    (%%rbx),%%xmm2			\n\t"\
			"movaps	0x10(%%rbx),%%xmm3			\n\t"\
			"movaps	    (%%rax),%%xmm0			\n\t"\
			"movaps	0x10(%%rax),%%xmm1			\n\t"\
			"movaps	    (%%rcx),%%xmm6			\n\t"\
			"movaps	0x10(%%rcx),%%xmm7			\n\t"\
			"movaps	%%xmm2,%%xmm4				\n\t"\
			"movaps	%%xmm3,%%xmm5				\n\t"\
			"\n\t"\
			"addpd	%%xmm6,%%xmm2				\n\t"\
			"addpd	%%xmm7,%%xmm3				\n\t"\
			"subpd	%%xmm6,%%xmm4				\n\t"\
			"subpd	%%xmm7,%%xmm5				\n\t"\
			"addpd	%%xmm2,%%xmm0				\n\t"\
			"addpd	%%xmm3,%%xmm1				\n\t"\
			"movaps	    (%%rdx),%%xmm6			\n\t"\
			"movaps	0x10(%%rdx),%%xmm7			\n\t"\
			"movaps	%%xmm0,     (%%rbx)	/*a<>b*/\n\t"\
			"movaps	%%xmm1,0x010(%%rbx)			\n\t"\
			"\n\t"\
			"mulpd	%%xmm6,%%xmm2				\n\t"\
			"mulpd	%%xmm6,%%xmm3				\n\t"\
			"mulpd	%%xmm7,%%xmm4				\n\t"\
			"mulpd	%%xmm7,%%xmm5				\n\t"\
			"addpd	%%xmm0,%%xmm2				\n\t"\
			"addpd	%%xmm1,%%xmm3				\n\t"\
			"\n\t"\
			"movaps	%%xmm2,%%xmm0				\n\t"\
			"movaps	%%xmm3,%%xmm1				\n\t"\
			"\n\t"\
			"subpd	%%xmm5,%%xmm2				\n\t"\
			"addpd	%%xmm4,%%xmm3				\n\t"\
			"addpd	%%xmm5,%%xmm0				\n\t"\
			"subpd	%%xmm4,%%xmm1				\n\t"\
			"\n\t"\
			"movaps	%%xmm2,     (%%rax)			\n\t"\
			"movaps	%%xmm3,0x010(%%rax)			\n\t"\
			"movaps	%%xmm0,     (%%rcx)			\n\t"\
			"movaps	%%xmm1,0x010(%%rcx)			\n\t"\
			"\n\t"\
			:					/* outputs: none */\
			: [__add] "m" (Xadd)	/* All inputs from memory addresses here */\
			 ,[__p01] "m" (Xp01)\
			 ,[__p02] "m" (Xp02)\
			 ,[__p03] "m" (Xp03)\
			 ,[__p04] "m" (Xp04)\
			 ,[__p08] "m" (Xp08)\
			 ,[__p16] "m" (Xp16)\
			 ,[__out] "m" (Xout)\
			 ,[__isrt2] "m" (Xisrt2)\
			 ,[__cc3] "m" (Xcc3)\
			: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
		);\
		}

	#define	SSE2_RADIX24_DIF_NOTWIDDLE(Xadd,Xp01,Xp02,Xp03,Xp04,Xp05,Xp06,Xp07,Xp08,Xp16,Xout,Xisrt2,Xcc3)\
	{\
	__asm__ volatile (\
			"movq	%[__out],%%rax	/* s1p00r */\n\t"\
			"movq	 %%rax,%%rbx	/* s1p00r */\n\t"\
			"movq	 %%rax,%%rcx	/* s1p08r */\n\t"\
			"movq	%[__cc3],%%rdx				\n\t"\
			"addq	$0x100,%%rbx	/* s1p08r */\n\t"\
			"addq	$0x200,%%rcx	/* s1p16r */\n\t"\
			"\n\t"\
		"/* On the DIF side it's always an *input* pair to the radix-3 DFT that gets swapped: */\n\t"\
		"/* SSE2_RADIX_03_DFT(s1p00r,s1p16r,s1p08r,cc3,s1p00r,s1p08r,s1p16r) */\n\t"\
			"movaps	    (%%rcx),%%xmm2	/*b<>c*/\n\t"\
			"movaps	0x10(%%rcx),%%xmm3			\n\t"\
			"movaps	    (%%rax),%%xmm0			\n\t"\
			"movaps	0x10(%%rax),%%xmm1			\n\t"\
			"movaps	    (%%rbx),%%xmm6			\n\t"\
			"movaps	0x10(%%rbx),%%xmm7			\n\t"\
			"movaps	%%xmm2,%%xmm4				\n\t"\
			"movaps	%%xmm3,%%xmm5				\n\t"\
			"\n\t"\
			"addpd	%%xmm6,%%xmm2				\n\t"\
			"addpd	%%xmm7,%%xmm3				\n\t"\
			"subpd	%%xmm6,%%xmm4				\n\t"\
			"subpd	%%xmm7,%%xmm5				\n\t"\
			"addpd	%%xmm2,%%xmm0				\n\t"\
			"addpd	%%xmm3,%%xmm1				\n\t"\
			"movaps	    (%%rdx),%%xmm6			\n\t"\
			"movaps	0x10(%%rdx),%%xmm7			\n\t"\
			"movaps	%%xmm0,     (%%rax)			\n\t"\
			"movaps	%%xmm1,0x010(%%rax)			\n\t"\
			"\n\t"\
			"mulpd	%%xmm6,%%xmm2				\n\t"\
			"mulpd	%%xmm6,%%xmm3				\n\t"\
			"mulpd	%%xmm7,%%xmm4				\n\t"\
			"mulpd	%%xmm7,%%xmm5				\n\t"\
			"addpd	%%xmm0,%%xmm2				\n\t"\
			"addpd	%%xmm1,%%xmm3				\n\t"\
			"\n\t"\
			"movaps	%%xmm2,%%xmm0				\n\t"\
			"movaps	%%xmm3,%%xmm1				\n\t"\
			"\n\t"\
			"subpd	%%xmm5,%%xmm2				\n\t"\
			"addpd	%%xmm4,%%xmm3				\n\t"\
			"addpd	%%xmm5,%%xmm0				\n\t"\
			"subpd	%%xmm4,%%xmm1				\n\t"\
			"\n\t"\
			"movaps	%%xmm2,     (%%rbx)			\n\t"\
			"movaps	%%xmm3,0x010(%%rbx)			\n\t"\
			"movaps	%%xmm0,     (%%rcx)			\n\t"\
			"movaps	%%xmm1,0x010(%%rcx)			\n\t"\
			"\n\t"\
		"/* SSE2_RADIX_03_DFT(s1p09r,s1p01r,s1p17r,cc3,s1p01r,s1p09r,s1p17r) */\n\t"\
			"addq	$0x020,%%rax	/* s1p01r */\n\t"\
			"addq	$0x020,%%rbx	/* s1p09r */\n\t"\
			"addq	$0x020,%%rcx	/* s1p17r */\n\t"\
			"\n\t"\
			"movaps	    (%%rax),%%xmm2	/*a<>b*/\n\t"\
			"movaps	0x10(%%rax),%%xmm3			\n\t"\
			"movaps	    (%%rbx),%%xmm0			\n\t"\
			"movaps	0x10(%%rbx),%%xmm1			\n\t"\
			"movaps	    (%%rcx),%%xmm6			\n\t"\
			"movaps	0x10(%%rcx),%%xmm7			\n\t"\
			"movaps	%%xmm2,%%xmm4				\n\t"\
			"movaps	%%xmm3,%%xmm5				\n\t"\
			"\n\t"\
			"addpd	%%xmm6,%%xmm2				\n\t"\
			"addpd	%%xmm7,%%xmm3				\n\t"\
			"subpd	%%xmm6,%%xmm4				\n\t"\
			"subpd	%%xmm7,%%xmm5				\n\t"\
			"addpd	%%xmm2,%%xmm0				\n\t"\
			"addpd	%%xmm3,%%xmm1				\n\t"\
			"movaps	    (%%rdx),%%xmm6			\n\t"\
			"movaps	0x10(%%rdx),%%xmm7			\n\t"\
			"movaps	%%xmm0,     (%%rax)			\n\t"\
			"movaps	%%xmm1,0x010(%%rax)			\n\t"\
			"\n\t"\
			"mulpd	%%xmm6,%%xmm2				\n\t"\
			"mulpd	%%xmm6,%%xmm3				\n\t"\
			"mulpd	%%xmm7,%%xmm4				\n\t"\
			"mulpd	%%xmm7,%%xmm5				\n\t"\
			"addpd	%%xmm0,%%xmm2				\n\t"\
			"addpd	%%xmm1,%%xmm3				\n\t"\
			"\n\t"\
			"movaps	%%xmm2,%%xmm0				\n\t"\
			"movaps	%%xmm3,%%xmm1				\n\t"\
			"\n\t"\
			"subpd	%%xmm5,%%xmm2				\n\t"\
			"addpd	%%xmm4,%%xmm3				\n\t"\
			"addpd	%%xmm5,%%xmm0				\n\t"\
			"subpd	%%xmm4,%%xmm1				\n\t"\
			"\n\t"\
			"movaps	%%xmm2,     (%%rbx)			\n\t"\
			"movaps	%%xmm3,0x010(%%rbx)			\n\t"\
			"movaps	%%xmm0,     (%%rcx)			\n\t"\
			"movaps	%%xmm1,0x010(%%rcx)			\n\t"\
			"\n\t"\
		"/* SSE2_RADIX_03_DFT(s1p18r,s1p10r,s1p02r,cc3,s1p02r,s1p10r,s1p18r) */\n\t"\
			"addq	$0x020,%%rax	/* s1p02r */\n\t"\
			"addq	$0x020,%%rbx	/* s1p10r */\n\t"\
			"addq	$0x020,%%rcx	/* s1p18r */\n\t"\
			"\n\t"\
			"movaps	    (%%rbx),%%xmm2	/*a<>c*/\n\t"\
			"movaps	0x10(%%rbx),%%xmm3			\n\t"\
			"movaps	    (%%rcx),%%xmm0			\n\t"\
			"movaps	0x10(%%rcx),%%xmm1			\n\t"\
			"movaps	    (%%rax),%%xmm6			\n\t"\
			"movaps	0x10(%%rax),%%xmm7			\n\t"\
			"movaps	%%xmm2,%%xmm4				\n\t"\
			"movaps	%%xmm3,%%xmm5				\n\t"\
			"\n\t"\
			"addpd	%%xmm6,%%xmm2				\n\t"\
			"addpd	%%xmm7,%%xmm3				\n\t"\
			"subpd	%%xmm6,%%xmm4				\n\t"\
			"subpd	%%xmm7,%%xmm5				\n\t"\
			"addpd	%%xmm2,%%xmm0				\n\t"\
			"addpd	%%xmm3,%%xmm1				\n\t"\
			"movaps	    (%%rdx),%%xmm6			\n\t"\
			"movaps	0x10(%%rdx),%%xmm7			\n\t"\
			"movaps	%%xmm0,     (%%rax)			\n\t"\
			"movaps	%%xmm1,0x010(%%rax)			\n\t"\
			"\n\t"\
			"mulpd	%%xmm6,%%xmm2				\n\t"\
			"mulpd	%%xmm6,%%xmm3				\n\t"\
			"mulpd	%%xmm7,%%xmm4				\n\t"\
			"mulpd	%%xmm7,%%xmm5				\n\t"\
			"addpd	%%xmm0,%%xmm2				\n\t"\
			"addpd	%%xmm1,%%xmm3				\n\t"\
			"\n\t"\
			"movaps	%%xmm2,%%xmm0				\n\t"\
			"movaps	%%xmm3,%%xmm1				\n\t"\
			"\n\t"\
			"subpd	%%xmm5,%%xmm2				\n\t"\
			"addpd	%%xmm4,%%xmm3				\n\t"\
			"addpd	%%xmm5,%%xmm0				\n\t"\
			"subpd	%%xmm4,%%xmm1				\n\t"\
			"\n\t"\
			"movaps	%%xmm2,     (%%rbx)			\n\t"\
			"movaps	%%xmm3,0x010(%%rbx)			\n\t"\
			"movaps	%%xmm0,     (%%rcx)			\n\t"\
			"movaps	%%xmm1,0x010(%%rcx)			\n\t"\
			"\n\t"\
		"/* SSE2_RADIX_03_DFT(s1p03r,s1p19r,s1p11r,cc3,s1p03r,s1p11r,s1p19r) */\n\t"\
			"addq	$0x020,%%rax	/* s1p03r */\n\t"\
			"addq	$0x020,%%rbx	/* s1p11r */\n\t"\
			"addq	$0x020,%%rcx	/* s1p19r */\n\t"\
			"\n\t"\
			"movaps	    (%%rcx),%%xmm2	/*b<>c*/\n\t"\
			"movaps	0x10(%%rcx),%%xmm3			\n\t"\
			"movaps	    (%%rax),%%xmm0			\n\t"\
			"movaps	0x10(%%rax),%%xmm1			\n\t"\
			"movaps	    (%%rbx),%%xmm6			\n\t"\
			"movaps	0x10(%%rbx),%%xmm7			\n\t"\
			"movaps	%%xmm2,%%xmm4				\n\t"\
			"movaps	%%xmm3,%%xmm5				\n\t"\
			"\n\t"\
			"addpd	%%xmm6,%%xmm2				\n\t"\
			"addpd	%%xmm7,%%xmm3				\n\t"\
			"subpd	%%xmm6,%%xmm4				\n\t"\
			"subpd	%%xmm7,%%xmm5				\n\t"\
			"addpd	%%xmm2,%%xmm0				\n\t"\
			"addpd	%%xmm3,%%xmm1				\n\t"\
			"movaps	    (%%rdx),%%xmm6			\n\t"\
			"movaps	0x10(%%rdx),%%xmm7			\n\t"\
			"movaps	%%xmm0,     (%%rax)			\n\t"\
			"movaps	%%xmm1,0x010(%%rax)			\n\t"\
			"\n\t"\
			"mulpd	%%xmm6,%%xmm2				\n\t"\
			"mulpd	%%xmm6,%%xmm3				\n\t"\
			"mulpd	%%xmm7,%%xmm4				\n\t"\
			"mulpd	%%xmm7,%%xmm5				\n\t"\
			"addpd	%%xmm0,%%xmm2				\n\t"\
			"addpd	%%xmm1,%%xmm3				\n\t"\
			"\n\t"\
			"movaps	%%xmm2,%%xmm0				\n\t"\
			"movaps	%%xmm3,%%xmm1				\n\t"\
			"\n\t"\
			"subpd	%%xmm5,%%xmm2				\n\t"\
			"addpd	%%xmm4,%%xmm3				\n\t"\
			"addpd	%%xmm5,%%xmm0				\n\t"\
			"subpd	%%xmm4,%%xmm1				\n\t"\
			"\n\t"\
			"movaps	%%xmm2,     (%%rbx)			\n\t"\
			"movaps	%%xmm3,0x010(%%rbx)			\n\t"\
			"movaps	%%xmm0,     (%%rcx)			\n\t"\
			"movaps	%%xmm1,0x010(%%rcx)			\n\t"\
			"\n\t"\
		"/* SSE2_RADIX_03_DFT(s1p12r,s1p04r,s1p20r,cc3,s1p04r,s1p12r,s1p20r) */\n\t"\
			"addq	$0x020,%%rax	/* s1p04r */\n\t"\
			"addq	$0x020,%%rbx	/* s1p12r */\n\t"\
			"addq	$0x020,%%rcx	/* s1p20r */\n\t"\
			"\n\t"\
			"movaps	    (%%rax),%%xmm2	/*a<>b*/\n\t"\
			"movaps	0x10(%%rax),%%xmm3			\n\t"\
			"movaps	    (%%rbx),%%xmm0			\n\t"\
			"movaps	0x10(%%rbx),%%xmm1			\n\t"\
			"movaps	    (%%rcx),%%xmm6			\n\t"\
			"movaps	0x10(%%rcx),%%xmm7			\n\t"\
			"movaps	%%xmm2,%%xmm4				\n\t"\
			"movaps	%%xmm3,%%xmm5				\n\t"\
			"\n\t"\
			"addpd	%%xmm6,%%xmm2				\n\t"\
			"addpd	%%xmm7,%%xmm3				\n\t"\
			"subpd	%%xmm6,%%xmm4				\n\t"\
			"subpd	%%xmm7,%%xmm5				\n\t"\
			"addpd	%%xmm2,%%xmm0				\n\t"\
			"addpd	%%xmm3,%%xmm1				\n\t"\
			"movaps	    (%%rdx),%%xmm6			\n\t"\
			"movaps	0x10(%%rdx),%%xmm7			\n\t"\
			"movaps	%%xmm0,     (%%rax)			\n\t"\
			"movaps	%%xmm1,0x010(%%rax)			\n\t"\
			"\n\t"\
			"mulpd	%%xmm6,%%xmm2				\n\t"\
			"mulpd	%%xmm6,%%xmm3				\n\t"\
			"mulpd	%%xmm7,%%xmm4				\n\t"\
			"mulpd	%%xmm7,%%xmm5				\n\t"\
			"addpd	%%xmm0,%%xmm2				\n\t"\
			"addpd	%%xmm1,%%xmm3				\n\t"\
			"\n\t"\
			"movaps	%%xmm2,%%xmm0				\n\t"\
			"movaps	%%xmm3,%%xmm1				\n\t"\
			"\n\t"\
			"subpd	%%xmm5,%%xmm2				\n\t"\
			"addpd	%%xmm4,%%xmm3				\n\t"\
			"addpd	%%xmm5,%%xmm0				\n\t"\
			"subpd	%%xmm4,%%xmm1				\n\t"\
			"\n\t"\
			"movaps	%%xmm2,     (%%rbx)			\n\t"\
			"movaps	%%xmm3,0x010(%%rbx)			\n\t"\
			"movaps	%%xmm0,     (%%rcx)			\n\t"\
			"movaps	%%xmm1,0x010(%%rcx)			\n\t"\
			"\n\t"\
		"/* SSE2_RADIX_03_DFT(s1p21r,s1p13r,s1p05r,cc3,s1p05r,s1p13r,s1p21r) */\n\t"\
			"addq	$0x020,%%rax	/* s1p05r */\n\t"\
			"addq	$0x020,%%rbx	/* s1p13r */\n\t"\
			"addq	$0x020,%%rcx	/* s1p21r */\n\t"\
			"\n\t"\
			"movaps	    (%%rbx),%%xmm2	/*a<>c*/\n\t"\
			"movaps	0x10(%%rbx),%%xmm3			\n\t"\
			"movaps	    (%%rcx),%%xmm0			\n\t"\
			"movaps	0x10(%%rcx),%%xmm1			\n\t"\
			"movaps	    (%%rax),%%xmm6			\n\t"\
			"movaps	0x10(%%rax),%%xmm7			\n\t"\
			"movaps	%%xmm2,%%xmm4				\n\t"\
			"movaps	%%xmm3,%%xmm5				\n\t"\
			"\n\t"\
			"addpd	%%xmm6,%%xmm2				\n\t"\
			"addpd	%%xmm7,%%xmm3				\n\t"\
			"subpd	%%xmm6,%%xmm4				\n\t"\
			"subpd	%%xmm7,%%xmm5				\n\t"\
			"addpd	%%xmm2,%%xmm0				\n\t"\
			"addpd	%%xmm3,%%xmm1				\n\t"\
			"movaps	    (%%rdx),%%xmm6			\n\t"\
			"movaps	0x10(%%rdx),%%xmm7			\n\t"\
			"movaps	%%xmm0,     (%%rax)			\n\t"\
			"movaps	%%xmm1,0x010(%%rax)			\n\t"\
			"\n\t"\
			"mulpd	%%xmm6,%%xmm2				\n\t"\
			"mulpd	%%xmm6,%%xmm3				\n\t"\
			"mulpd	%%xmm7,%%xmm4				\n\t"\
			"mulpd	%%xmm7,%%xmm5				\n\t"\
			"addpd	%%xmm0,%%xmm2				\n\t"\
			"addpd	%%xmm1,%%xmm3				\n\t"\
			"\n\t"\
			"movaps	%%xmm2,%%xmm0				\n\t"\
			"movaps	%%xmm3,%%xmm1				\n\t"\
			"\n\t"\
			"subpd	%%xmm5,%%xmm2				\n\t"\
			"addpd	%%xmm4,%%xmm3				\n\t"\
			"addpd	%%xmm5,%%xmm0				\n\t"\
			"subpd	%%xmm4,%%xmm1				\n\t"\
			"\n\t"\
			"movaps	%%xmm2,     (%%rbx)			\n\t"\
			"movaps	%%xmm3,0x010(%%rbx)			\n\t"\
			"movaps	%%xmm0,     (%%rcx)			\n\t"\
			"movaps	%%xmm1,0x010(%%rcx)			\n\t"\
			"\n\t"\
		"/* SSE2_RADIX_03_DFT(s1p06r,s1p22r,s1p14r,cc3,s1p06r,s1p14r,s1p22r) */\n\t"\
			"addq	$0x020,%%rax	/* s1p06r */\n\t"\
			"addq	$0x020,%%rbx	/* s1p14r */\n\t"\
			"addq	$0x020,%%rcx	/* s1p22r */\n\t"\
			"\n\t"\
			"movaps	    (%%rcx),%%xmm2	/*b<>c*/\n\t"\
			"movaps	0x10(%%rcx),%%xmm3			\n\t"\
			"movaps	    (%%rax),%%xmm0			\n\t"\
			"movaps	0x10(%%rax),%%xmm1			\n\t"\
			"movaps	    (%%rbx),%%xmm6			\n\t"\
			"movaps	0x10(%%rbx),%%xmm7			\n\t"\
			"movaps	%%xmm2,%%xmm4				\n\t"\
			"movaps	%%xmm3,%%xmm5				\n\t"\
			"\n\t"\
			"addpd	%%xmm6,%%xmm2				\n\t"\
			"addpd	%%xmm7,%%xmm3				\n\t"\
			"subpd	%%xmm6,%%xmm4				\n\t"\
			"subpd	%%xmm7,%%xmm5				\n\t"\
			"addpd	%%xmm2,%%xmm0				\n\t"\
			"addpd	%%xmm3,%%xmm1				\n\t"\
			"movaps	    (%%rdx),%%xmm6			\n\t"\
			"movaps	0x10(%%rdx),%%xmm7			\n\t"\
			"movaps	%%xmm0,     (%%rax)			\n\t"\
			"movaps	%%xmm1,0x010(%%rax)			\n\t"\
			"\n\t"\
			"mulpd	%%xmm6,%%xmm2				\n\t"\
			"mulpd	%%xmm6,%%xmm3				\n\t"\
			"mulpd	%%xmm7,%%xmm4				\n\t"\
			"mulpd	%%xmm7,%%xmm5				\n\t"\
			"addpd	%%xmm0,%%xmm2				\n\t"\
			"addpd	%%xmm1,%%xmm3				\n\t"\
			"\n\t"\
			"movaps	%%xmm2,%%xmm0				\n\t"\
			"movaps	%%xmm3,%%xmm1				\n\t"\
			"\n\t"\
			"subpd	%%xmm5,%%xmm2				\n\t"\
			"addpd	%%xmm4,%%xmm3				\n\t"\
			"addpd	%%xmm5,%%xmm0				\n\t"\
			"subpd	%%xmm4,%%xmm1				\n\t"\
			"\n\t"\
			"movaps	%%xmm2,     (%%rbx)			\n\t"\
			"movaps	%%xmm3,0x010(%%rbx)			\n\t"\
			"movaps	%%xmm0,     (%%rcx)			\n\t"\
			"movaps	%%xmm1,0x010(%%rcx)			\n\t"\
			"\n\t"\
		"/* SSE2_RADIX_03_DFT(s1p15r,s1p07r,s1p23r,cc3,s1p07r,s1p15r,s1p23r) */\n\t"\
			"addq	$0x020,%%rax	/* s1p07r */\n\t"\
			"addq	$0x020,%%rbx	/* s1p15r */\n\t"\
			"addq	$0x020,%%rcx	/* s1p23r */\n\t"\
			"\n\t"\
			"movaps	    (%%rax),%%xmm2	/*a<>b*/\n\t"\
			"movaps	0x10(%%rax),%%xmm3			\n\t"\
			"movaps	    (%%rbx),%%xmm0			\n\t"\
			"movaps	0x10(%%rbx),%%xmm1			\n\t"\
			"movaps	    (%%rcx),%%xmm6			\n\t"\
			"movaps	0x10(%%rcx),%%xmm7			\n\t"\
			"movaps	%%xmm2,%%xmm4				\n\t"\
			"movaps	%%xmm3,%%xmm5				\n\t"\
			"\n\t"\
			"addpd	%%xmm6,%%xmm2				\n\t"\
			"addpd	%%xmm7,%%xmm3				\n\t"\
			"subpd	%%xmm6,%%xmm4				\n\t"\
			"subpd	%%xmm7,%%xmm5				\n\t"\
			"addpd	%%xmm2,%%xmm0				\n\t"\
			"addpd	%%xmm3,%%xmm1				\n\t"\
			"movaps	    (%%rdx),%%xmm6			\n\t"\
			"movaps	0x10(%%rdx),%%xmm7			\n\t"\
			"movaps	%%xmm0,     (%%rax)			\n\t"\
			"movaps	%%xmm1,0x010(%%rax)			\n\t"\
			"\n\t"\
			"mulpd	%%xmm6,%%xmm2				\n\t"\
			"mulpd	%%xmm6,%%xmm3				\n\t"\
			"mulpd	%%xmm7,%%xmm4				\n\t"\
			"mulpd	%%xmm7,%%xmm5				\n\t"\
			"addpd	%%xmm0,%%xmm2				\n\t"\
			"addpd	%%xmm1,%%xmm3				\n\t"\
			"\n\t"\
			"movaps	%%xmm2,%%xmm0				\n\t"\
			"movaps	%%xmm3,%%xmm1				\n\t"\
			"\n\t"\
			"subpd	%%xmm5,%%xmm2				\n\t"\
			"addpd	%%xmm4,%%xmm3				\n\t"\
			"addpd	%%xmm5,%%xmm0				\n\t"\
			"subpd	%%xmm4,%%xmm1				\n\t"\
			"\n\t"\
			"movaps	%%xmm2,     (%%rbx)			\n\t"\
			"movaps	%%xmm3,0x010(%%rbx)			\n\t"\
			"movaps	%%xmm0,     (%%rcx)			\n\t"\
			"movaps	%%xmm1,0x010(%%rcx)			\n\t"\
			"\n\t"\
			"movq	%[__isrt2],%%rsi			\n\t"\
			"\n\t"\
		"/* For the radix-8 DIF DFTs, the input offsets always have the same pattern; outputs are permuted */\n\t"\
		"/* SSE2_RADIX8_DIF_0TWIDDLE( i[0-7] = s1p00r + 0x[0a4e82c6]0, o[0-7] = add0 + p[01235476]) */\n\t"\
			"movq	%[__out],%%rax	/* i0 = s1p00r */\n\t"\
			"movq	$0x40  ,%%rbx	/* i2 */	\n\t"\
			"movq	$0x80  ,%%rcx	/* i4 */	\n\t"\
			"movq	$0xc0  ,%%rdx	/* i6 */	\n\t"\
			"addq	%%rax,%%rbx					\n\t"\
			"addq	%%rax,%%rcx					\n\t"\
			"addq	%%rax,%%rdx					\n\t"\
			"\n\t"\
			"/* Do the p0,p4 combo: */			\n\t"\
			"\n\t"\
			"movaps	    (%%rax),%%xmm0			\n\t"\
			"movaps	    (%%rcx),%%xmm4			\n\t"\
			"movaps	0x10(%%rax),%%xmm1			\n\t"\
			"movaps	0x10(%%rcx),%%xmm5			\n\t"\
			"movaps	%%xmm0,%%xmm2				\n\t"\
			"movaps	%%xmm1,%%xmm3				\n\t"\
			"\n\t"\
			"addpd	%%xmm4,%%xmm0				\n\t"\
			"addpd	%%xmm5,%%xmm1				\n\t"\
			"subpd	%%xmm4,%%xmm2				\n\t"\
			"subpd	%%xmm5,%%xmm3				\n\t"\
			"\n\t"\
			"/* Do the p2,6 combo: */			\n\t"\
			"\n\t"\
			"movaps	    (%%rbx),%%xmm4			\n\t"\
			"movaps	0x10(%%rbx),%%xmm5			\n\t"\
			"movaps	%%xmm4,%%xmm6				\n\t"\
			"movaps	%%xmm5,%%xmm7				\n\t"\
			"\n\t"\
			"addpd	     (%%rdx),%%xmm4			\n\t"\
			"addpd	0x010(%%rdx),%%xmm5			\n\t"\
			"subpd	     (%%rdx),%%xmm6			\n\t"\
			"subpd	0x010(%%rdx),%%xmm7			\n\t"\
			"\n\t"\
			"subpd	%%xmm4,%%xmm0				\n\t"\
			"subpd	%%xmm7,%%xmm2				\n\t"\
			"subpd	%%xmm5,%%xmm1				\n\t"\
			"subpd	%%xmm6,%%xmm3				\n\t"\
			"movaps	%%xmm0,     (%%rcx)			\n\t"\
			"movaps	%%xmm2,     (%%rbx)			\n\t"\
			"movaps	%%xmm1,0x010(%%rcx)			\n\t"\
			"movaps	%%xmm3,0x010(%%rdx)			\n\t"\
			"addpd	%%xmm4,%%xmm4				\n\t"\
			"addpd	%%xmm7,%%xmm7				\n\t"\
			"addpd	%%xmm5,%%xmm5				\n\t"\
			"addpd	%%xmm6,%%xmm6				\n\t"\
			"addpd	%%xmm0,%%xmm4				\n\t"\
			"addpd	%%xmm2,%%xmm7				\n\t"\
			"addpd	%%xmm1,%%xmm5				\n\t"\
			"addpd	%%xmm3,%%xmm6				\n\t"\
			"movaps	%%xmm4,     (%%rax)			\n\t"\
			"movaps	%%xmm7,     (%%rdx)			\n\t"\
			"movaps	%%xmm5,0x010(%%rax)			\n\t"\
			"movaps	%%xmm6,0x010(%%rbx)			\n\t"\
			"\n\t"\
			"/* Inline of SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO_B(r8,c1): Do 2nd of 2 radix-4 subtransforms, but now keep outputs in registers: */\n\t"\
			"\n\t"\
			"movq	$0xe0,%%rbx		/* i3 */	\n\t"\
			"movq	$0x20,%%rcx		/* i5 */	\n\t"\
			"movq	$0x60,%%rdx		/* i7 */	\n\t"\
			"addq	%%rax,%%rbx					\n\t"\
			"addq	%%rax,%%rcx					\n\t"\
			"addq	%%rax,%%rdx					\n\t"\
			"addq	$0xa0,%%rax		/* i1 */	\n\t"\
			"\n\t"\
			"/* Do the p1,p5 combo: */			\n\t"\
			"\n\t"\
			"movaps	    (%%rax),%%xmm0			\n\t"\
			"movaps	    (%%rcx),%%xmm4			\n\t"\
			"movaps	0x10(%%rax),%%xmm1			\n\t"\
			"movaps	0x10(%%rcx),%%xmm5			\n\t"\
			"movaps	%%xmm0,%%xmm2				\n\t"\
			"movaps	%%xmm1,%%xmm3				\n\t"\
			"\n\t"\
			"addpd	%%xmm4,%%xmm0				\n\t"\
			"addpd	%%xmm5,%%xmm1				\n\t"\
			"subpd	%%xmm4,%%xmm2				\n\t"\
			"subpd	%%xmm5,%%xmm3				\n\t"\
			"\n\t"\
			"/* Do the p3,7 combo: */			\n\t"\
			"\n\t"\
			"movaps	    (%%rbx),%%xmm4			\n\t"\
			"movaps	0x10(%%rbx),%%xmm5			\n\t"\
			"movaps	%%xmm4,%%xmm6				\n\t"\
			"movaps	%%xmm5,%%xmm7				\n\t"\
			"\n\t"\
			"subpd	    (%%rdx),%%xmm4			\n\t"\
			"subpd	0x10(%%rdx),%%xmm5			\n\t"\
			"addpd	    (%%rdx),%%xmm6			\n\t"\
			"addpd	0x10(%%rdx),%%xmm7			\n\t"\
			"\n\t"\
			"/* Finish radix-4 butterfly and store just the 1st of the 8 outputs into output-array slots: */\n\t"\
			"\n\t"\
			"subpd	%%xmm6,%%xmm0				\n\t"\
			"subpd	%%xmm7,%%xmm1				\n\t"\
			"subpd	%%xmm5,%%xmm2				\n\t"\
			"subpd	%%xmm4,%%xmm3				\n\t"\
			"addpd	%%xmm6,%%xmm6				\n\t"\
			"addpd	%%xmm5,%%xmm5				\n\t"\
			"addpd	%%xmm7,%%xmm7				\n\t"\
			"addpd	%%xmm4,%%xmm4				\n\t"\
			"addpd	%%xmm0,%%xmm6				\n\t"\
			"addpd	%%xmm2,%%xmm5				\n\t"\
			"addpd	%%xmm1,%%xmm7				\n\t"\
			"addpd	%%xmm3,%%xmm4				\n\t"\
			"movaps	%%xmm6,    (%%rax)			\n\t"\
			"movaps	%%xmm7,0x10(%%rax)			\n\t"\
			"\n\t"\
			"movaps	%%xmm2,%%xmm6				\n\t"\
			"movaps	%%xmm5,%%xmm7				\n\t"\
			"subpd	%%xmm4,%%xmm2				\n\t"\
			"subpd	%%xmm3,%%xmm5				\n\t"\
			"addpd	%%xmm6,%%xmm4				\n\t"\
			"addpd	%%xmm7,%%xmm3				\n\t"\
			"\n\t"\
			"movaps	(%%rsi),%%xmm6				\n\t"\
			"mulpd	%%xmm6,%%xmm2				\n\t"\
			"mulpd	%%xmm6,%%xmm5				\n\t"\
			"mulpd	%%xmm6,%%xmm4				\n\t"\
			"mulpd	%%xmm6,%%xmm3				\n\t"\
			"\n\t"\
			"/* Inline the contents of SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_A(r0): Combine the 2 radix-4 subtransforms and write outputs to temp array: */\n\t"\
			"/***** t0,1,2,3,4,5,6,7 in __i0,2,4,6 (eax,ebx,ecx,edx) *****/\n\t"\
			"/***** t8,9,a,b,c,d,e,f in __i1,xmm2,4,0,1,5,3; xmm6,7 FREE */\n\t"\
			"subq	$0xa0,%%rax		/* i0 */	\n\t"\
			"movq	%%rax , %%rdi				\n\t"\
			"\n\t"\
			"movslq	%[__p05],%%rax	\n\t"\
			"movslq	%[__p04],%%rbx	\n\t"\
			"movslq	%[__p07],%%rcx	\n\t"\
			"movslq	%[__p06],%%rdx	\n\t"\
			"shlq	$3,%%rax		/* Pointer offset for floating doubles */\n\t"\
			"shlq	$3,%%rbx		\n\t"\
			"shlq	$3,%%rcx		\n\t"\
			"shlq	$3,%%rdx		\n\t"\
			"addq	%[__add],%%rax	/* o4 */\n\t"\
			"addq	%[__add],%%rbx	/* o5 */\n\t"\
			"addq	%[__add],%%rcx	/* o6 */\n\t"\
			"addq	%[__add],%%rdx	/* o7 */\n\t"\
			"\n\t"\
			"movaps	0x40(%%rdi),%%xmm6		\n\t"\
			"movaps	0x50(%%rdi),%%xmm7		\n\t"\
			"subpd   %%xmm2,%%xmm6				\n\t"\
			"subpd   %%xmm4,%%xmm7				\n\t"\
			"addpd   %%xmm2,%%xmm2				\n\t"\
			"addpd   %%xmm4,%%xmm4				\n\t"\
			"addpd   %%xmm6,%%xmm2				\n\t"\
			"addpd   %%xmm7,%%xmm4				\n\t"\
			"\n\t"\
			"movaps	%%xmm6,    (%%rbx)			\n\t"\
			"movaps	%%xmm7,0x10(%%rbx)			\n\t"\
			"movaps	%%xmm2,    (%%rax)			\n\t"\
			"movaps	%%xmm4,0x10(%%rax)			\n\t"\
			"\n\t"\
			"/* Can't do 2nd set side-by-side with 1st, since only 2 registers (xmm6,7) available for temp storage: */\n\t"\
			"movaps	0xc0(%%rdi),%%xmm6		\n\t"\
			"movaps	0xd0(%%rdi),%%xmm7		\n\t"\
			"subpd   %%xmm3,%%xmm6				\n\t"\
			"subpd   %%xmm5,%%xmm7				\n\t"\
			"addpd   %%xmm3,%%xmm3				\n\t"\
			"addpd   %%xmm5,%%xmm5				\n\t"\
			"addpd   %%xmm6,%%xmm3				\n\t"\
			"addpd   %%xmm7,%%xmm5				\n\t"\
			"\n\t"\
			"movaps	%%xmm6,    (%%rcx)			\n\t"\
			"movaps	%%xmm7,0x10(%%rdx)			\n\t"\
			"movaps	%%xmm3,    (%%rdx)			\n\t"\
			"movaps	%%xmm5,0x10(%%rcx)			\n\t"\
			"\n\t"\
			"subq	%[__add],%%rbx	/* p4 */\n\t"\
			"subq	%%rbx,%%rax	/* o1 = add+p1 */\n\t"\
			"subq	%%rbx,%%rcx	/* o2 = add+p3 */\n\t"\
			"subq	%%rbx,%%rdx	/* o3 = add+p2 */\n\t"\
			"movq	%[__add],%%rbx	/* o0 */\n\t"\
			"\n\t"\
			"movaps	    (%%rdi),%%xmm6			\n\t"\
			"movaps	0x80(%%rdi),%%xmm4		\n\t"\
			"movaps	0x10(%%rdi),%%xmm7			\n\t"\
			"movaps	0x90(%%rdi),%%xmm5		\n\t"\
			"movaps	0xa0(%%rdi),%%xmm2		\n\t"\
			"movaps	0xb0(%%rdi),%%xmm3		\n\t"\
			"subpd   %%xmm2,%%xmm6				\n\t"\
			"subpd   %%xmm1,%%xmm4				\n\t"\
			"subpd   %%xmm3,%%xmm7				\n\t"\
			"subpd   %%xmm0,%%xmm5				\n\t"\
			"addpd   %%xmm2,%%xmm2				\n\t"\
			"addpd   %%xmm1,%%xmm1				\n\t"\
			"addpd   %%xmm3,%%xmm3				\n\t"\
			"addpd   %%xmm0,%%xmm0				\n\t"\
			"addpd   %%xmm6,%%xmm2				\n\t"\
			"addpd   %%xmm4,%%xmm1				\n\t"\
			"addpd   %%xmm7,%%xmm3				\n\t"\
			"addpd   %%xmm5,%%xmm0				\n\t"\
			"\n\t"\
			"movaps	%%xmm6,    (%%rax)	/* a<>b */\n\t"\
			"movaps	%%xmm4,    (%%rdx)	/* c<>d */\n\t"\
			"movaps	%%xmm7,0x10(%%rax)			\n\t"\
			"movaps	%%xmm5,0x10(%%rcx)			\n\t"\
			"movaps	%%xmm2,    (%%rbx)			\n\t"\
			"movaps	%%xmm1,    (%%rcx)			\n\t"\
			"movaps	%%xmm3,0x10(%%rbx)			\n\t"\
			"movaps	%%xmm0,0x10(%%rdx)			\n\t"\
			"\n\t"\
		"/* SSE2_RADIX8_DIF_0TWIDDLE( i[0-7] = s1p16r + 0x[0a4e82c6]0, o[0-7] = add8 + p[54762310]) */\n\t"\
			"movq	%%rdi,%%rax	/* s1p00r */\n\t"\
			"addq	$0x200 ,%%rax	/* i0 = s1p16r */	\n\t"\
			"movq	$0x40  ,%%rbx	/* i2 */\n\t"\
			"movq	$0x80  ,%%rcx	/* i4 */\n\t"\
			"movq	$0xc0  ,%%rdx	/* i6 */\n\t"\
			"addq	%%rax,%%rbx				\n\t"\
			"addq	%%rax,%%rcx				\n\t"\
			"addq	%%rax,%%rdx				\n\t"\
			"\n\t"\
			"/* Do the p0,p4 combo: */		\n\t"\
			"\n\t"\
			"movaps	    (%%rax),%%xmm0		\n\t"\
			"movaps	    (%%rcx),%%xmm4		\n\t"\
			"movaps	0x10(%%rax),%%xmm1		\n\t"\
			"movaps	0x10(%%rcx),%%xmm5		\n\t"\
			"movaps	%%xmm0,%%xmm2			\n\t"\
			"movaps	%%xmm1,%%xmm3			\n\t"\
			"\n\t"\
			"addpd	%%xmm4,%%xmm0			\n\t"\
			"addpd	%%xmm5,%%xmm1			\n\t"\
			"subpd	%%xmm4,%%xmm2			\n\t"\
			"subpd	%%xmm5,%%xmm3			\n\t"\
			"\n\t"\
			"/* Do the p2,6 combo: */		\n\t"\
			"\n\t"\
			"movaps	    (%%rbx),%%xmm4		\n\t"\
			"movaps	0x10(%%rbx),%%xmm5		\n\t"\
			"movaps	%%xmm4,%%xmm6			\n\t"\
			"movaps	%%xmm5,%%xmm7			\n\t"\
			"\n\t"\
			"addpd	     (%%rdx),%%xmm4		\n\t"\
			"addpd	0x010(%%rdx),%%xmm5		\n\t"\
			"subpd	     (%%rdx),%%xmm6		\n\t"\
			"subpd	0x010(%%rdx),%%xmm7		\n\t"\
			"\n\t"\
			"subpd	%%xmm4,%%xmm0			\n\t"\
			"subpd	%%xmm7,%%xmm2			\n\t"\
			"subpd	%%xmm5,%%xmm1			\n\t"\
			"subpd	%%xmm6,%%xmm3			\n\t"\
			"movaps	%%xmm0,     (%%rcx)		\n\t"\
			"movaps	%%xmm2,     (%%rbx)		\n\t"\
			"movaps	%%xmm1,0x010(%%rcx)		\n\t"\
			"movaps	%%xmm3,0x010(%%rdx)		\n\t"\
			"addpd	%%xmm4,%%xmm4			\n\t"\
			"addpd	%%xmm7,%%xmm7			\n\t"\
			"addpd	%%xmm5,%%xmm5			\n\t"\
			"addpd	%%xmm6,%%xmm6			\n\t"\
			"addpd	%%xmm0,%%xmm4			\n\t"\
			"addpd	%%xmm2,%%xmm7			\n\t"\
			"addpd	%%xmm1,%%xmm5			\n\t"\
			"addpd	%%xmm3,%%xmm6			\n\t"\
			"movaps	%%xmm4,     (%%rax)		\n\t"\
			"movaps	%%xmm7,     (%%rdx)		\n\t"\
			"movaps	%%xmm5,0x010(%%rax)		\n\t"\
			"movaps	%%xmm6,0x010(%%rbx)		\n\t"\
			"\n\t"\
			"/* Inline of SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO_B(r8,c1): Do 2nd of 2 radix-4 subtransforms, but now keep outputs in registers: */\n\t"\
			"\n\t"\
			"movq	$0xe0,%%rbx		/* i3 */\n\t"\
			"movq	$0x20,%%rcx		/* i5 */\n\t"\
			"movq	$0x60,%%rdx		/* i7 */\n\t"\
			"addq	%%rax,%%rbx				\n\t"\
			"addq	%%rax,%%rcx				\n\t"\
			"addq	%%rax,%%rdx				\n\t"\
			"addq	$0xa0,%%rax		/* i1 */\n\t"\
			"\n\t"\
			"/* Do the p1,p5 combo: */		\n\t"\
			"\n\t"\
			"movaps	    (%%rax),%%xmm0		\n\t"\
			"movaps	    (%%rcx),%%xmm4		\n\t"\
			"movaps	0x10(%%rax),%%xmm1		\n\t"\
			"movaps	0x10(%%rcx),%%xmm5		\n\t"\
			"movaps	%%xmm0,%%xmm2			\n\t"\
			"movaps	%%xmm1,%%xmm3			\n\t"\
			"\n\t"\
			"addpd	%%xmm4,%%xmm0			\n\t"\
			"addpd	%%xmm5,%%xmm1			\n\t"\
			"subpd	%%xmm4,%%xmm2			\n\t"\
			"subpd	%%xmm5,%%xmm3			\n\t"\
			"\n\t"\
			"/* Do the p3,7 combo: */		\n\t"\
			"\n\t"\
			"movaps	    (%%rbx),%%xmm4		\n\t"\
			"movaps	0x10(%%rbx),%%xmm5		\n\t"\
			"movaps	%%xmm4,%%xmm6			\n\t"\
			"movaps	%%xmm5,%%xmm7			\n\t"\
			"\n\t"\
			"subpd	    (%%rdx),%%xmm4		\n\t"\
			"subpd	0x10(%%rdx),%%xmm5		\n\t"\
			"addpd	    (%%rdx),%%xmm6		\n\t"\
			"addpd	0x10(%%rdx),%%xmm7		\n\t"\
			"\n\t"\
			"/* Finish radix-4 butterfly andstore just the 1st of the 8 outputs into output-array slots: */\n\t"\
			"\n\t"\
			"subpd	%%xmm6,%%xmm0			\n\t"\
			"subpd	%%xmm7,%%xmm1			\n\t"\
			"subpd	%%xmm5,%%xmm2			\n\t"\
			"subpd	%%xmm4,%%xmm3			\n\t"\
			"addpd	%%xmm6,%%xmm6			\n\t"\
			"addpd	%%xmm5,%%xmm5			\n\t"\
			"addpd	%%xmm7,%%xmm7			\n\t"\
			"addpd	%%xmm4,%%xmm4			\n\t"\
			"addpd	%%xmm0,%%xmm6			\n\t"\
			"addpd	%%xmm2,%%xmm5			\n\t"\
			"addpd	%%xmm1,%%xmm7			\n\t"\
			"addpd	%%xmm3,%%xmm4			\n\t"\
			"movaps	%%xmm6,    (%%rax)		\n\t"\
			"movaps	%%xmm7,0x10(%%rax)		\n\t"\
			"\n\t"\
			"movaps	%%xmm2,%%xmm6			\n\t"\
			"movaps	%%xmm5,%%xmm7			\n\t"\
			"subpd	%%xmm4,%%xmm2			\n\t"\
			"subpd	%%xmm3,%%xmm5			\n\t"\
			"addpd	%%xmm6,%%xmm4			\n\t"\
			"addpd	%%xmm7,%%xmm3			\n\t"\
			"\n\t"\
			"movaps	(%%rsi),%%xmm6			\n\t"\
			"mulpd	%%xmm6,%%xmm2			\n\t"\
			"mulpd	%%xmm6,%%xmm5			\n\t"\
			"mulpd	%%xmm6,%%xmm4			\n\t"\
			"mulpd	%%xmm6,%%xmm3			\n\t"\
			"\n\t"\
			"/* Inline the contents of SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_A(r0): Combine the 2 radix-4 subtransforms and write outputs to temp array: */\n\t"\
			"/***** t0,1,2,3,4,5,6,7 in __i0,2,4,6 (eax,ebx,ecx,edx) *****/\n\t"\
			"/***** t8,9,a,b,c,d,e,f in __i1,xmm2,4,0,1,5,3; xmm6,7 FREE */\n\t"\
			"subq	$0xa0,%%rax		/* i0 */	\n\t"\
			"movq	%%rax , %%rdi			\n\t"\
			"\n\t"\
			"movslq	%[__p08],%%rdx			\n\t"\
			"shlq	$3,%%rdx		/* Pointer offset for floating doubles */\n\t"\
			"movslq	%[__p02],%%rax			\n\t"\
			"movslq	%[__p03],%%rbx			\n\t"\
			"movslq	%[__p01],%%rcx			\n\t"\
			"addq	%[__add],%%rdx	/* o7 */\n\t"\
			"shlq	$3,%%rax		/* Pointer offset for floating doubles */\n\t"\
			"shlq	$3,%%rbx				\n\t"\
			"shlq	$3,%%rcx				\n\t"\
			"addq	%%rdx,%%rax	/* o4 */	\n\t"\
			"addq	%%rdx,%%rbx	/* o5 */	\n\t"\
			"addq	%%rdx,%%rcx	/* o6 */	\n\t"\
			"\n\t"\
			"movaps	0x40(%%rdi),%%xmm6		\n\t"\
			"movaps	0x50(%%rdi),%%xmm7		\n\t"\
			"subpd   %%xmm2,%%xmm6			\n\t"\
			"subpd   %%xmm4,%%xmm7			\n\t"\
			"addpd   %%xmm2,%%xmm2			\n\t"\
			"addpd   %%xmm4,%%xmm4			\n\t"\
			"addpd   %%xmm6,%%xmm2			\n\t"\
			"addpd   %%xmm7,%%xmm4			\n\t"\
			"\n\t"\
			"movaps	%%xmm6,    (%%rbx)		\n\t"\
			"movaps	%%xmm7,0x10(%%rbx)		\n\t"\
			"movaps	%%xmm2,    (%%rax)		\n\t"\
			"movaps	%%xmm4,0x10(%%rax)		\n\t"\
			"\n\t"\
			"/* Can't do 2nd set side-by-side with 1st, since only 2 registers (xmm6,7) available for temp storage: */\n\t"\
			"movaps	0xc0(%%rdi),%%xmm6		\n\t"\
			"movaps	0xd0(%%rdi),%%xmm7		\n\t"\
			"subpd   %%xmm3,%%xmm6			\n\t"\
			"subpd   %%xmm5,%%xmm7			\n\t"\
			"addpd   %%xmm3,%%xmm3			\n\t"\
			"addpd   %%xmm5,%%xmm5			\n\t"\
			"addpd   %%xmm6,%%xmm3			\n\t"\
			"addpd   %%xmm7,%%xmm5			\n\t"\
			"\n\t"\
			"movaps	%%xmm6,    (%%rcx)		\n\t"\
			"movaps	%%xmm7,0x10(%%rdx)		\n\t"\
			"movaps	%%xmm3,    (%%rdx)		\n\t"\
			"movaps	%%xmm5,0x10(%%rcx)		\n\t"\
			"\n\t"\
			"movslq	%[__p04],%%rsi			\n\t"\
			"shlq	$3,%%rsi		/* Pointer offset for floating doubles */\n\t"\
			"addq	%%rsi,%%rax	/* o3 = add8+p6 */\n\t"\
			"addq	%%rsi,%%rbx	/* o2 = add8+p7 */\n\t"\
			"addq	%%rsi,%%rcx	/* o0 = add8+p5 */\n\t"\
			"addq	%%rsi,%%rdx	/* o1 = add8+p4 */\n\t"\
			"\n\t"\
			"movq	%[__isrt2],%%rsi	/* Restore isrt2 ptr */\n\t"\
			"\n\t"\
			"movaps	    (%%rdi),%%xmm6		\n\t"\
			"movaps	0x80(%%rdi),%%xmm4		\n\t"\
			"movaps	0x10(%%rdi),%%xmm7		\n\t"\
			"movaps	0x90(%%rdi),%%xmm5		\n\t"\
			"movaps	0xa0(%%rdi),%%xmm2		\n\t"\
			"movaps	0xb0(%%rdi),%%xmm3		\n\t"\
			"subpd   %%xmm2,%%xmm6			\n\t"\
			"subpd   %%xmm1,%%xmm4			\n\t"\
			"subpd   %%xmm3,%%xmm7			\n\t"\
			"subpd   %%xmm0,%%xmm5			\n\t"\
			"addpd   %%xmm2,%%xmm2			\n\t"\
			"addpd   %%xmm1,%%xmm1			\n\t"\
			"addpd   %%xmm3,%%xmm3			\n\t"\
			"addpd   %%xmm0,%%xmm0			\n\t"\
			"addpd   %%xmm6,%%xmm2			\n\t"\
			"addpd   %%xmm4,%%xmm1			\n\t"\
			"addpd   %%xmm7,%%xmm3			\n\t"\
			"addpd   %%xmm5,%%xmm0			\n\t"\
			"\n\t"\
			"movaps	%%xmm6,    (%%rdx)	/* abcd -> cdba */\n\t"\
			"movaps	%%xmm4,    (%%rbx)		\n\t"\
			"movaps	%%xmm7,0x10(%%rdx)		\n\t"\
			"movaps	%%xmm5,0x10(%%rax)		\n\t"\
			"movaps	%%xmm2,    (%%rcx)		\n\t"\
			"movaps	%%xmm1,    (%%rax)		\n\t"\
			"movaps	%%xmm3,0x10(%%rcx)		\n\t"\
			"movaps	%%xmm0,0x10(%%rbx)		\n\t"\
			"\n\t"\
		"/* SSE2_RADIX8_DIF_0TWIDDLE( i[0-7] = s1p08r + 0x[0a4e82c6]0, o[0-7] = add16+ p[23107645]) */\n\t"\
			"movq	%%rdi,%%rax	/* s1p00r */\n\t"\
			"subq	$0x100 ,%%rax	/* i0 = s1p08r */	\n\t"\
			"movq	$0x40  ,%%rbx	/* i2 */\n\t"\
			"movq	$0x80  ,%%rcx	/* i4 */\n\t"\
			"movq	$0xc0  ,%%rdx	/* i6 */\n\t"\
			"addq	%%rax,%%rbx				\n\t"\
			"addq	%%rax,%%rcx				\n\t"\
			"addq	%%rax,%%rdx				\n\t"\
			"\n\t"\
			"/* Do the p0,p4 combo: */		\n\t"\
			"\n\t"\
			"movaps	    (%%rax),%%xmm0		\n\t"\
			"movaps	    (%%rcx),%%xmm4		\n\t"\
			"movaps	0x10(%%rax),%%xmm1		\n\t"\
			"movaps	0x10(%%rcx),%%xmm5		\n\t"\
			"movaps	%%xmm0,%%xmm2			\n\t"\
			"movaps	%%xmm1,%%xmm3			\n\t"\
			"\n\t"\
			"addpd	%%xmm4,%%xmm0			\n\t"\
			"addpd	%%xmm5,%%xmm1			\n\t"\
			"subpd	%%xmm4,%%xmm2			\n\t"\
			"subpd	%%xmm5,%%xmm3			\n\t"\
			"\n\t"\
			"/* Do the p2,6 combo: */		\n\t"\
			"\n\t"\
			"movaps	    (%%rbx),%%xmm4		\n\t"\
			"movaps	0x10(%%rbx),%%xmm5		\n\t"\
			"movaps	%%xmm4,%%xmm6			\n\t"\
			"movaps	%%xmm5,%%xmm7			\n\t"\
			"\n\t"\
			"addpd	     (%%rdx),%%xmm4		\n\t"\
			"addpd	0x010(%%rdx),%%xmm5		\n\t"\
			"subpd	     (%%rdx),%%xmm6		\n\t"\
			"subpd	0x010(%%rdx),%%xmm7		\n\t"\
			"\n\t"\
			"subpd	%%xmm4,%%xmm0			\n\t"\
			"subpd	%%xmm7,%%xmm2			\n\t"\
			"subpd	%%xmm5,%%xmm1			\n\t"\
			"subpd	%%xmm6,%%xmm3			\n\t"\
			"movaps	%%xmm0,     (%%rcx)		\n\t"\
			"movaps	%%xmm2,     (%%rbx)		\n\t"\
			"movaps	%%xmm1,0x010(%%rcx)		\n\t"\
			"movaps	%%xmm3,0x010(%%rdx)		\n\t"\
			"addpd	%%xmm4,%%xmm4			\n\t"\
			"addpd	%%xmm7,%%xmm7			\n\t"\
			"addpd	%%xmm5,%%xmm5			\n\t"\
			"addpd	%%xmm6,%%xmm6			\n\t"\
			"addpd	%%xmm0,%%xmm4			\n\t"\
			"addpd	%%xmm2,%%xmm7			\n\t"\
			"addpd	%%xmm1,%%xmm5			\n\t"\
			"addpd	%%xmm3,%%xmm6			\n\t"\
			"movaps	%%xmm4,     (%%rax)		\n\t"\
			"movaps	%%xmm7,     (%%rdx)		\n\t"\
			"movaps	%%xmm5,0x010(%%rax)		\n\t"\
			"movaps	%%xmm6,0x010(%%rbx)		\n\t"\
			"\n\t"\
			"/* Inline of SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO_B(r8,c1): Do 2nd of 2 radix-4 subtransforms, but now keep outputs in registers: */\n\t"\
			"\n\t"\
			"movq	$0xe0,%%rbx		/* i3 */\n\t"\
			"movq	$0x20,%%rcx		/* i5 */\n\t"\
			"movq	$0x60,%%rdx		/* i7 */\n\t"\
			"addq	%%rax,%%rbx				\n\t"\
			"addq	%%rax,%%rcx				\n\t"\
			"addq	%%rax,%%rdx				\n\t"\
			"addq	$0xa0,%%rax		/* i1 */\n\t"\
			"\n\t"\
			"/* Do the p1,p5 combo: */		\n\t"\
			"\n\t"\
			"movaps	    (%%rax),%%xmm0		\n\t"\
			"movaps	    (%%rcx),%%xmm4		\n\t"\
			"movaps	0x10(%%rax),%%xmm1		\n\t"\
			"movaps	0x10(%%rcx),%%xmm5		\n\t"\
			"movaps	%%xmm0,%%xmm2			\n\t"\
			"movaps	%%xmm1,%%xmm3			\n\t"\
			"\n\t"\
			"addpd	%%xmm4,%%xmm0			\n\t"\
			"addpd	%%xmm5,%%xmm1			\n\t"\
			"subpd	%%xmm4,%%xmm2			\n\t"\
			"subpd	%%xmm5,%%xmm3			\n\t"\
			"\n\t"\
			"/* Do the p3,7 combo: */		\n\t"\
			"\n\t"\
			"movaps	    (%%rbx),%%xmm4		\n\t"\
			"movaps	0x10(%%rbx),%%xmm5		\n\t"\
			"movaps	%%xmm4,%%xmm6			\n\t"\
			"movaps	%%xmm5,%%xmm7			\n\t"\
			"\n\t"\
			"subpd	    (%%rdx),%%xmm4		\n\t"\
			"subpd	0x10(%%rdx),%%xmm5		\n\t"\
			"addpd	    (%%rdx),%%xmm6		\n\t"\
			"addpd	0x10(%%rdx),%%xmm7		\n\t"\
			"\n\t"\
			"/* Finish radix-4 butterfly andstore just the 1st of the 8 outputs into output-array slots: */\n\t"\
			"\n\t"\
			"subpd	%%xmm6,%%xmm0			\n\t"\
			"subpd	%%xmm7,%%xmm1			\n\t"\
			"subpd	%%xmm5,%%xmm2			\n\t"\
			"subpd	%%xmm4,%%xmm3			\n\t"\
			"addpd	%%xmm6,%%xmm6			\n\t"\
			"addpd	%%xmm5,%%xmm5			\n\t"\
			"addpd	%%xmm7,%%xmm7			\n\t"\
			"addpd	%%xmm4,%%xmm4			\n\t"\
			"addpd	%%xmm0,%%xmm6			\n\t"\
			"addpd	%%xmm2,%%xmm5			\n\t"\
			"addpd	%%xmm1,%%xmm7			\n\t"\
			"addpd	%%xmm3,%%xmm4			\n\t"\
			"movaps	%%xmm6,    (%%rax)		\n\t"\
			"movaps	%%xmm7,0x10(%%rax)		\n\t"\
			"\n\t"\
			"movaps	%%xmm2,%%xmm6			\n\t"\
			"movaps	%%xmm5,%%xmm7			\n\t"\
			"subpd	%%xmm4,%%xmm2			\n\t"\
			"subpd	%%xmm3,%%xmm5			\n\t"\
			"addpd	%%xmm6,%%xmm4			\n\t"\
			"addpd	%%xmm7,%%xmm3			\n\t"\
			"\n\t"\
			"movaps	(%%rsi),%%xmm6			\n\t"\
			"mulpd	%%xmm6,%%xmm2			\n\t"\
			"mulpd	%%xmm6,%%xmm5			\n\t"\
			"mulpd	%%xmm6,%%xmm4			\n\t"\
			"mulpd	%%xmm6,%%xmm3			\n\t"\
			"\n\t"\
			"/* Inline the contents of SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_A(r0): Combine the 2 radix-4 subtransforms and write outputs to temp array: */\n\t"\
			"/***** t0,1,2,3,4,5,6,7 in __i0,2,4,6 (eax,ebx,ecx,edx) *****/\n\t"\
			"/***** t8,9,a,b,c,d,e,f in __i1,xmm2,4,0,1,5,3; xmm6,7 FREE */\n\t"\
			"subq	$0xa0,%%rax		/* i0 */	\n\t"\
			"movq	%%rax , %%rdi			\n\t"\
			"\n\t"\
			"movslq	%[__p04],%%rcx			\n\t"\
			"movslq	%[__p16],%%rdx			\n\t"\
			"addq	%%rdx,%%rcx		/* p20 */\n\t"\
			"movslq	%[__p03],%%rax			\n\t"\
			"movslq	%[__p02],%%rbx			\n\t"\
			"movslq	%[__p01],%%rdx			\n\t"\
			"shlq	$3,%%rcx				\n\t"\
			"shlq	$3,%%rax		/* Pointer offset for floating doubles */\n\t"\
			"shlq	$3,%%rbx				\n\t"\
			"shlq	$3,%%rdx				\n\t"\
			"addq	%[__add],%%rcx	/* o6 = add16+p4 */\n\t"\
			"addq	%%rcx,%%rax		/* o4 = add16+p7 */\n\t"\
			"addq	%%rcx,%%rbx		/* o5 = add16+p6 */\n\t"\
			"addq	%%rcx,%%rdx		/* o7 = add16+p5 */\n\t"\
			"\n\t"\
			"movaps	0x40(%%rdi),%%xmm6		\n\t"\
			"movaps	0x50(%%rdi),%%xmm7		\n\t"\
			"subpd   %%xmm2,%%xmm6			\n\t"\
			"subpd   %%xmm4,%%xmm7			\n\t"\
			"addpd   %%xmm2,%%xmm2			\n\t"\
			"addpd   %%xmm4,%%xmm4			\n\t"\
			"addpd   %%xmm6,%%xmm2			\n\t"\
			"addpd   %%xmm7,%%xmm4			\n\t"\
			"\n\t"\
			"movaps	%%xmm6,    (%%rbx)		\n\t"\
			"movaps	%%xmm7,0x10(%%rbx)		\n\t"\
			"movaps	%%xmm2,    (%%rax)		\n\t"\
			"movaps	%%xmm4,0x10(%%rax)		\n\t"\
			"\n\t"\
			"/* Can't do 2nd set side-by-side with 1st, since only 2 registers (xmm6,7) available for temp storage: */\n\t"\
			"movaps	0xc0(%%rdi),%%xmm6		\n\t"\
			"movaps	0xd0(%%rdi),%%xmm7		\n\t"\
			"subpd   %%xmm3,%%xmm6			\n\t"\
			"subpd   %%xmm5,%%xmm7			\n\t"\
			"addpd   %%xmm3,%%xmm3			\n\t"\
			"addpd   %%xmm5,%%xmm5			\n\t"\
			"addpd   %%xmm6,%%xmm3			\n\t"\
			"addpd   %%xmm7,%%xmm5			\n\t"\
			"\n\t"\
			"movaps	%%xmm6,    (%%rcx)		\n\t"\
			"movaps	%%xmm7,0x10(%%rdx)		\n\t"\
			"movaps	%%xmm3,    (%%rdx)		\n\t"\
			"movaps	%%xmm5,0x10(%%rcx)		\n\t"\
			"\n\t"\
			"movslq	%[__p04],%%rsi		/* overwrite __isrt2 ptr */\n\t"\
			"shlq	$3,%%rsi		/* Pointer offset for floating doubles */\n\t"\
			"subq	%%rsi,%%rax	/* o1 = add16+p3 */\n\t"\
			"subq	%%rsi,%%rbx	/* o0 = add16+p2 */\n\t"\
			"subq	%%rsi,%%rcx	/* o3 = add16+p0 */\n\t"\
			"subq	%%rsi,%%rdx	/* o2 = add16+p1 */\n\t"\
			"\n\t"\
			"movaps	    (%%rdi),%%xmm6		\n\t"\
			"movaps	0x80(%%rdi),%%xmm4		\n\t"\
			"movaps	0x10(%%rdi),%%xmm7		\n\t"\
			"movaps	0x90(%%rdi),%%xmm5		\n\t"\
			"movaps	0xa0(%%rdi),%%xmm2		\n\t"\
			"movaps	0xb0(%%rdi),%%xmm3		\n\t"\
			"subpd   %%xmm2,%%xmm6			\n\t"\
			"subpd   %%xmm1,%%xmm4			\n\t"\
			"subpd   %%xmm3,%%xmm7			\n\t"\
			"subpd   %%xmm0,%%xmm5			\n\t"\
			"addpd   %%xmm2,%%xmm2			\n\t"\
			"addpd   %%xmm1,%%xmm1			\n\t"\
			"addpd   %%xmm3,%%xmm3			\n\t"\
			"addpd   %%xmm0,%%xmm0			\n\t"\
			"addpd   %%xmm6,%%xmm2			\n\t"\
			"addpd   %%xmm4,%%xmm1			\n\t"\
			"addpd   %%xmm7,%%xmm3			\n\t"\
			"addpd   %%xmm5,%%xmm0			\n\t"\
			"\n\t"\
			"movaps	%%xmm6,    (%%rax)	/* abcd -> badc */\n\t"\
			"movaps	%%xmm4,    (%%rdx)		\n\t"\
			"movaps	%%xmm7,0x10(%%rax)		\n\t"\
			"movaps	%%xmm5,0x10(%%rcx)		\n\t"\
			"movaps	%%xmm2,    (%%rbx)		\n\t"\
			"movaps	%%xmm1,    (%%rcx)		\n\t"\
			"movaps	%%xmm3,0x10(%%rbx)		\n\t"\
			"movaps	%%xmm0,0x10(%%rdx)		\n\t"\
			"\n\t"\
			:					/* outputs: none */\
			: [__add] "m" (Xadd)	/* All inputs from memory addresses here */\
			 ,[__p01] "m" (Xp01)\
			 ,[__p02] "m" (Xp02)\
			 ,[__p03] "m" (Xp03)\
			 ,[__p04] "m" (Xp04)\
			 ,[__p05] "m" (Xp05)\
			 ,[__p06] "m" (Xp06)\
			 ,[__p07] "m" (Xp07)\
			 ,[__p08] "m" (Xp08)\
			 ,[__p16] "m" (Xp16)\
			 ,[__out] "m" (Xout)\
			 ,[__isrt2] "m" (Xisrt2)\
			 ,[__cc3] "m" (Xcc3)\
			: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
		);\
		}

#endif	/* radix24_ditN_cy_dif1_gcc_h_included */

