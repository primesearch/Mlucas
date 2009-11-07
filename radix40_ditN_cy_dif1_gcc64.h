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
#ifndef radix40_ditN_cy_dif1_gcc_h_included
#define radix40_ditN_cy_dif1_gcc_h_included

	#define	SSE2_RADIX40_DIT_NOTWIDDLE(Xadd,Xp01,Xp02,Xp03,Xp04,Xp08,Xp16,Xp24,Xp32,Xr00,Xcc1,Xa0,Xa1,Xa2,Xa3,Xa4)\
	{\
	__asm__ volatile (\
	"/* SSE2_RADIX8_DIT_0TWIDDLE(add0+p[0,1,3,2,7,6,5,4], r00) */\n\t"\
		"movslq	%[__p04],%%rdi	\n\t"\
		"movq	%[__add],%%rax	/* Use eax as base address throughout */\n\t"\
		"shlq	$3,%%rdi		\n\t"\
		"movslq	%[__p01],%%r15	\n\t"\
		"movslq	%[__p02],%%rcx	\n\t"\
		"movslq	%[__p03],%%rdx	\n\t"\
		"addq	%%rdi,%%rax		/* eax <- add0+p04 */\n\t"\
		"shlq	$3,%%r15		/* Pointer offset for floating doubles */\n\t"\
		"shlq	$3,%%rcx		\n\t"\
		"shlq	$3,%%rdx		\n\t"\
		"addq	%%rax,%%r15		/* ebx <- add0+p05 */\n\t"\
		"addq	%%rax,%%rcx		/* ecx <- add0+p06 */\n\t"\
		"addq	%%rax,%%rdx		/* edx <- add0+p07 */\n\t"\
		"movq	%[__r00],%%rsi	\n\t"\
		"/* MSVC macro assumes add8+p[7,6,5,4] in eax,ebx,ecx,edx, but here get add0+p[4,5,6,7] in those registers, so replace eax <-> edx and ebx <-> ecx: */\n\t"\
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
		"movaps	    (%%r15),%%xmm4			\n\t"\
		"movaps	0x10(%%r15),%%xmm5			\n\t"\
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
		"movaps	0xa00(%%rsi),%%xmm1	/* isrt2 */	\n\t"\
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
		"subq	%%rdi,%%r15		/* add3 = add0+p01 */\n\t"\
		"subq	%%rdi,%%rcx		/* add1 = add0+p02 */\n\t"\
		"subq	%%rdi,%%rdx		/* add0 = add0+p03 */\n\t"\
		"\n\t"\
		"movaps	    (%%rax),%%xmm0			\n\t"\
		"movaps	0x10(%%rax),%%xmm1			\n\t"\
		"movaps	%%xmm0,%%xmm2				\n\t"\
		"movaps	%%xmm1,%%xmm3				\n\t"\
		"addpd	    (%%r15),%%xmm2			\n\t"\
		"addpd	0x10(%%r15),%%xmm3			\n\t"\
		"subpd	    (%%r15),%%xmm0			\n\t"\
		"subpd	0x10(%%r15),%%xmm1			\n\t"\
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
	"/* SSE2_RADIX8_DIT_0TWIDDLE([add0+p16]+p[3,2,1,0,5,4,6,7], r16) */\n\t"\
		"movslq	%[__p16],%%r10	\n\t"\
		"shlq	$3,%%r10		/* Pointer offset for floating doubles */\n\t"\
		"addq	%%rdi,%%r10		/*         p04+p16 */\n\t"\
		"addq	%%r10,%%rax		/* edx <- add16+p04 */\n\t"\
		"addq	%%r10,%%r15		/* ecx <- add16+p05 */\n\t"\
		"addq	%%r10,%%rcx		/* ebx <- add16+p06 */\n\t"\
		"addq	%%r10,%%rdx		/* eax <- add16+p07 */\n\t"\
		"addq	$0x100,%%rsi	/* r16 */\n\t"\
		"/* MSVC macro assumes add8+p[5,4,6,7] in eax,ebx,ecx,edx, but here get add+p[4,5,6,7], so replace eax <-> ebx: */\n\t"\
		"/* Do the p0,p4 combo: */\n\t"\
		"movaps	    (%%r15),%%xmm0			\n\t"\
		"movaps	0x10(%%r15),%%xmm1			\n\t"\
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
		"movaps	0x900(%%rsi),%%xmm1	/* isrt2 */	\n\t"\
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
		"/* MSVC macro assumes add8+p[3,2,1,0] in eax,ebx,ecx,edx, but here get add8+p[0,1,2,3], so replace eax <-> edx and ebx <-> ecx: */\n\t"\
		"subq	%%rdi,%%rax		/* add0 = add     */\n\t"\
		"subq	%%rdi,%%r15		/* add1 = add+p01 */\n\t"\
		"subq	%%rdi,%%rcx		/* add3 = add+p02 */\n\t"\
		"subq	%%rdi,%%rdx		/* add2 = add+p03 */\n\t"\
		"\n\t"\
		"movaps	    (%%rdx),%%xmm0			\n\t"\
		"movaps	0x10(%%rdx),%%xmm1			\n\t"\
		"movaps	%%xmm0,%%xmm2				\n\t"\
		"movaps	%%xmm1,%%xmm3				\n\t"\
		"addpd	    (%%rcx),%%xmm2			\n\t"\
		"addpd	0x10(%%rcx),%%xmm3			\n\t"\
		"subpd	    (%%rcx),%%xmm0			\n\t"\
		"subpd	0x10(%%rcx),%%xmm1			\n\t"\
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
		"movaps	    (%%r15),%%xmm4			\n\t"\
		"movaps	0x10(%%r15),%%xmm5			\n\t"\
		"movaps	%%xmm4,%%xmm6				\n\t"\
		"movaps	%%xmm5,%%xmm7				\n\t"\
		"addpd	    (%%rax),%%xmm6			\n\t"\
		"addpd	0x10(%%rax),%%xmm7			\n\t"\
		"subpd	    (%%rax),%%xmm4			\n\t"\
		"subpd	0x10(%%rax),%%xmm5			\n\t"\
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
	"/* SSE2_RADIX8_DIT_0TWIDDLE([add0+p32]+p[1,0,2,3,6,7,4,5], r32) */\n\t"\
		"subq	%%rdi,%%r10		/* p20 - p04 = p16 */\n\t"\
		"subq	%%r10,%%rax		/* edx <- add0 +p00 */\n\t"\
		"subq	%%r10,%%r15		/* ecx <- add0 +p01 */\n\t"\
		"subq	%%r10,%%rcx		/* ebx <- add0 +p02 */\n\t"\
		"subq	%%r10,%%rdx		/* eax <- add0 +p03 */\n\t"\
		"movslq	%[__p32],%%r10	\n\t"\
		"shlq	$3,%%r10		/* Pointer offset for floating doubles */\n\t"\
		"addq	%%rdi,%%r10		/* p36 = p04 + p32 */\n\t"\
		"addq	%%r10,%%rax		/* edx <- add32+p04 */\n\t"\
		"addq	%%r10,%%r15		/* ecx <- add32+p05 */\n\t"\
		"addq	%%r10,%%rcx		/* ebx <- add32+p06 */\n\t"\
		"addq	%%r10,%%rdx		/* eax <- add32+p07 */\n\t"\
		"addq	$0x100,%%rsi	/* r32 */\n\t"\
		"/* MSVC macro assumes add8+p[6,7,4,5] in eax,ebx,ecx,edx, but here get add+p[4,5,6,7], so replace [eax,ebx] <-> [ecx,edx]: */\n\t"\
		"/* Do the p0,p4 combo: */\n\t"\
		"movaps	    (%%rcx),%%xmm0			\n\t"\
		"movaps	0x10(%%rcx),%%xmm1			\n\t"\
		"movaps	%%xmm0,%%xmm2				\n\t"\
		"movaps	%%xmm1,%%xmm3				\n\t"\
		"addpd	    (%%rdx),%%xmm2			\n\t"\
		"addpd	0x10(%%rdx),%%xmm3			\n\t"\
		"subpd	    (%%rdx),%%xmm0			\n\t"\
		"subpd	0x10(%%rdx),%%xmm1			\n\t"\
		"\n\t"\
		"movaps	    (%%rax),%%xmm4			\n\t"\
		"movaps	0x10(%%rax),%%xmm5			\n\t"\
		"movaps	%%xmm4,%%xmm6				\n\t"\
		"movaps	%%xmm5,%%xmm7				\n\t"\
		"addpd	    (%%r15),%%xmm6			\n\t"\
		"addpd	0x10(%%r15),%%xmm7			\n\t"\
		"subpd	    (%%r15),%%xmm4			\n\t"\
		"subpd	0x10(%%r15),%%xmm5			\n\t"\
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
		"movaps	0x800(%%rsi),%%xmm1	/* isrt2 */	\n\t"\
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
		"/* MSVC macro assumes add8+p[1,0,2,3] in eax,ebx,ecx,edx, but here get add8+p[0,1,2,3], so replace eax <-> ebx: */\n\t"\
		"subq	%%rdi,%%rax		/* add0 = add     */\n\t"\
		"subq	%%rdi,%%r15		/* add1 = add+p01 */\n\t"\
		"subq	%%rdi,%%rcx		/* add3 = add+p02 */\n\t"\
		"subq	%%rdi,%%rdx		/* add2 = add+p03 */\n\t"\
		"\n\t"\
		"movaps	    (%%r15),%%xmm0			\n\t"\
		"movaps	0x10(%%r15),%%xmm1			\n\t"\
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
	"/* SSE2_RADIX8_DIT_0TWIDDLE([add0+p08]+p[6,7,4,5,2,3,0,1], r48) */\n\t"\
		"subq	%%rdi,%%r10		/* p36 - p04 = p32 */\n\t"\
		"subq	%%r10,%%rax		/* edx <- add0 +p00 */\n\t"\
		"subq	%%r10,%%r15		/* ecx <- add0 +p01 */\n\t"\
		"subq	%%r10,%%rcx		/* ebx <- add0 +p02 */\n\t"\
		"subq	%%r10,%%rdx		/* eax <- add0 +p03 */\n\t"\
		"movslq	%[__p08],%%r10	\n\t"\
		"shlq	$3,%%r10		/* Pointer offset for floating doubles */\n\t"\
		"addq	%%r10,%%rax		/* edx <- add08+p00 */\n\t"\
		"addq	%%r10,%%r15		/* ecx <- add08+p01 */\n\t"\
		"addq	%%r10,%%rcx		/* ebx <- add08+p02 */\n\t"\
		"addq	%%r10,%%rdx		/* eax <- add08+p03 */\n\t"\
		"addq	$0x100,%%rsi	/* r48 */\n\t"\
		"/* MSVC macro assumes add8+p[2,3,0,1] in eax,ebx,ecx,edx, but here get add+p[0,1,2,3], so replace [eax,ebx] <-> [ecx,edx]: */\n\t"\
		"/* Do the p0,p4 combo: */\n\t"\
		"movaps	    (%%rcx),%%xmm0			\n\t"\
		"movaps	0x10(%%rcx),%%xmm1			\n\t"\
		"movaps	%%xmm0,%%xmm2				\n\t"\
		"movaps	%%xmm1,%%xmm3				\n\t"\
		"addpd	    (%%rdx),%%xmm2			\n\t"\
		"addpd	0x10(%%rdx),%%xmm3			\n\t"\
		"subpd	    (%%rdx),%%xmm0			\n\t"\
		"subpd	0x10(%%rdx),%%xmm1			\n\t"\
		"\n\t"\
		"movaps	    (%%rax),%%xmm4			\n\t"\
		"movaps	0x10(%%rax),%%xmm5			\n\t"\
		"movaps	%%xmm4,%%xmm6				\n\t"\
		"movaps	%%xmm5,%%xmm7				\n\t"\
		"addpd	    (%%r15),%%xmm6			\n\t"\
		"addpd	0x10(%%r15),%%xmm7			\n\t"\
		"subpd	    (%%r15),%%xmm4			\n\t"\
		"subpd	0x10(%%r15),%%xmm5			\n\t"\
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
		"movaps	0x700(%%rsi),%%xmm1	/* isrt2 */	\n\t"\
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
		"/* MSVC macro assumes add8+p[6,7,4,5] in eax,ebx,ecx,edx, but here get add+p[4,5,6,7], so replace [eax,ebx] <-> [ecx,edx]: */\n\t"\
		"addq	%%rdi,%%rax		/* add0 = add+p04 */\n\t"\
		"addq	%%rdi,%%r15		/* add1 = add+p05 */\n\t"\
		"addq	%%rdi,%%rcx		/* add3 = add+p06 */\n\t"\
		"addq	%%rdi,%%rdx		/* add2 = add+p07 */\n\t"\
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
		"addpd	    (%%r15),%%xmm6			\n\t"\
		"addpd	0x10(%%r15),%%xmm7			\n\t"\
		"subpd	    (%%r15),%%xmm4			\n\t"\
		"subpd	0x10(%%r15),%%xmm5			\n\t"\
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
	"/* SSE2_RADIX8_DIT_0TWIDDLE([add0+p24]+p[4,5,7,6,0,1,3,2], r64) */\n\t"\
		"movslq	%[__p16],%%r10	\n\t"\
		"shlq	$3,%%r10		/* Pointer offset for floating doubles */\n\t"\
		"subq	%%rdi,%%r10		/* p16 - p04 = p12 */\n\t"\
		"addq	%%r10,%%rax		/* edx <- add24+p00 */\n\t"\
		"addq	%%r10,%%r15		/* ecx <- add24+p01 */\n\t"\
		"addq	%%r10,%%rcx		/* ebx <- add24+p02 */\n\t"\
		"addq	%%r10,%%rdx		/* eax <- add24+p03 */\n\t"\
		"addq	$0x100,%%rsi	/* r64 */\n\t"\
		"/* MSVC macro assumes add8+p[0,1,3,2] in eax,ebx,ecx,edx, but here get add+p[0,1,2,3], so replace ecx <-> edx: */\n\t"\
		"/* Do the p0,p4 combo: */\n\t"\
		"movaps	    (%%rax),%%xmm0			\n\t"\
		"movaps	0x10(%%rax),%%xmm1			\n\t"\
		"movaps	%%xmm0,%%xmm2				\n\t"\
		"movaps	%%xmm1,%%xmm3				\n\t"\
		"addpd	    (%%r15),%%xmm2			\n\t"\
		"addpd	0x10(%%r15),%%xmm3			\n\t"\
		"subpd	    (%%r15),%%xmm0			\n\t"\
		"subpd	0x10(%%r15),%%xmm1			\n\t"\
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
		"movaps	0x600(%%rsi),%%xmm1	/* isrt2 */	\n\t"\
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
		"/* MSVC macro assumes add8+p[4,5,7,6] in eax,ebx,ecx,edx, but here get add+p[4,5,6,7], so replace ecx <-> edx: */\n\t"\
		"addq	%%rdi,%%rax		/* add0 = add+p04 */\n\t"\
		"addq	%%rdi,%%r15		/* add1 = add+p05 */\n\t"\
		"addq	%%rdi,%%rcx		/* add3 = add+p06 */\n\t"\
		"addq	%%rdi,%%rdx		/* add2 = add+p07 */\n\t"\
		"\n\t"\
		"movaps	    (%%rax),%%xmm0			\n\t"\
		"movaps	0x10(%%rax),%%xmm1			\n\t"\
		"movaps	%%xmm0,%%xmm2				\n\t"\
		"movaps	%%xmm1,%%xmm3				\n\t"\
		"addpd	    (%%r15),%%xmm2			\n\t"\
		"addpd	0x10(%%r15),%%xmm3			\n\t"\
		"subpd	    (%%r15),%%xmm0			\n\t"\
		"subpd	0x10(%%r15),%%xmm1			\n\t"\
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
	"/*********************************************************************/\n\t"\
	"/******************     Now do 8 radix-5 DFTs:    ********************/\n\t"\
	"/*********************************************************************/\n\t"\
		"/* SSE2_RADIX_05_DFT_0TWIDDLE(r00,r16,r32,r48,r64,cc1,s1p00r,s1p16r,s1p32r,s1p08r,s1p24r) */\n\t"\
		"movq	%[__r00],%%rsi		\n\t"\
		"movq	%%rsi,%%rax			\n\t"\
		"movq	%%rsi,%%r15			\n\t"\
		"movq	%%rsi,%%rcx			\n\t"\
		"movq	%%rsi,%%rdx			\n\t"\
		"addq	$0x100,%%rax		\n\t"\
		"addq	$0x200,%%r15		\n\t"\
		"addq	$0x300,%%rcx		\n\t"\
		"addq	$0x400,%%rdx		\n\t"\
		"movq	%[__a0],%%rdi	/* Out0 <- s1p00r */\n\t"\
		"movaps	    (%%rax),%%xmm0	\n\t"\
		"movaps	0x10(%%rax),%%xmm1	\n\t"\
		"movaps	    (%%r15),%%xmm2	\n\t"\
		"movaps	0x10(%%r15),%%xmm3	\n\t"\
		"movaps	    (%%rcx),%%xmm4	\n\t"\
		"movaps	0x10(%%rcx),%%xmm5	\n\t"\
		"movaps	    (%%rdx),%%xmm6	\n\t"\
		"movaps	0x10(%%rdx),%%xmm7	\n\t"\
		"subpd	%%xmm6,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm1		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm0,%%xmm6		\n\t"\
		"addpd	%%xmm1,%%xmm7		\n\t"\
		"subpd	%%xmm4,%%xmm2		\n\t"\
		"subpd	%%xmm5,%%xmm3		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm2,%%xmm4		\n\t"\
		"addpd	%%xmm3,%%xmm5		\n\t"\
		"/*movq	%[__cc1],%%rax		Need rax for a1, instead use that cc1 = r00 + 0xa10 */\n\t"\
		"subpd	%%xmm4,%%xmm6		\n\t"\
		"subpd	%%xmm5,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm5		\n\t"\
		"addpd	    (%%rsi),%%xmm4	\n\t"\
		"addpd	0x10(%%rsi),%%xmm5	\n\t"\
		"movaps	%%xmm4,    (%%rdi)	\n\t"\
		"movaps	%%xmm5,0x10(%%rdi)	\n\t"\
		"mulpd	0xa20(%%rsi),%%xmm6	\n\t"\
		"mulpd	0xa20(%%rsi),%%xmm7	\n\t"\
		"subpd	     (%%rsi),%%xmm4	\n\t"\
		"subpd	0x010(%%rsi),%%xmm5	\n\t"\
		"mulpd	0xa10(%%rsi),%%xmm4	\n\t"\
		"mulpd	0xa10(%%rsi),%%xmm5	\n\t"\
		"addpd	     (%%rdi),%%xmm4	\n\t"\
		"addpd	0x010(%%rdi),%%xmm5	\n\t"\
		"subpd	%%xmm6,%%xmm4		\n\t"\
		"subpd	%%xmm7,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm6		\n\t"\
		"addpd	%%xmm5,%%xmm7		\n\t"\
		"movaps	%%xmm4,    (%%rsi)	\n\t"\
		"movaps	%%xmm5,0x10(%%rsi)	\n\t"\
		"movaps	%%xmm0,%%xmm4		\n\t"\
		"movaps	%%xmm1,%%xmm5		\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"mulpd	0xa30(%%rsi),%%xmm0	\n\t"\
		"mulpd	0xa30(%%rsi),%%xmm1	\n\t"\
		"mulpd	0xa40(%%rsi),%%xmm2	\n\t"\
		"mulpd	0xa40(%%rsi),%%xmm3	\n\t"\
		"mulpd	0xa50(%%rsi),%%xmm4	\n\t"\
		"mulpd	0xa50(%%rsi),%%xmm5	\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"movaps	    (%%rsi),%%xmm4	\n\t"\
		"movaps	0x10(%%rsi),%%xmm5	\n\t"\
		"movq	%[__a2],%%rax	/* Out1 <- s1p16r */\n\t"\
		"movq	%[__a3],%%rdx	/* Out4 <- s1p24r */\n\t"\
		"subpd	%%xmm3,%%xmm6		\n\t"\
		"subpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"movaps	%%xmm6,    (%%rax)	\n\t"\
		"movaps	%%xmm7,0x10(%%rdx)	\n\t"\
		"addpd	%%xmm6,%%xmm3		\n\t"\
		"addpd	%%xmm7,%%xmm2		\n\t"\
		"movaps	%%xmm3,    (%%rdx)	\n\t"\
		"movaps	%%xmm2,0x10(%%rax)	\n\t"\
		"movq	%[__a4],%%r15	/* Out2 <- s1p32r */\n\t"\
		"movq	%[__a1],%%rcx	/* Out3 <- s1p08r */\n\t"\
		"subpd	%%xmm1,%%xmm4		\n\t"\
		"subpd	%%xmm0,%%xmm5		\n\t"\
		"addpd	%%xmm1,%%xmm1		\n\t"\
		"addpd	%%xmm0,%%xmm0		\n\t"\
		"movaps	%%xmm4,    (%%r15)	\n\t"\
		"movaps	%%xmm5,0x10(%%rcx)	\n\t"\
		"addpd	%%xmm4,%%xmm1		\n\t"\
		"addpd	%%xmm5,%%xmm0		\n\t"\
		"movaps	%%xmm1,    (%%rcx)	\n\t"\
		"movaps	%%xmm0,0x10(%%r15)	\n\t"\
		"							\n\t"\
		"/* SSE2_RADIX_05_DFT_0TWIDDLE(r02,r18,r34,r50,r66,cc1,s1p25r,s1p01r,s1p17r,s1p33r,s1p09r) */\n\t"\
		"addq	$0x20,%%rsi	/* r02 */\n\t"\
		"movq	%%rsi,%%rax			\n\t"\
		"movq	%%rsi,%%r15			\n\t"\
		"movq	%%rsi,%%rcx			\n\t"\
		"movq	%%rsi,%%rdx			\n\t"\
		"addq	$0x100,%%rax		\n\t"\
		"addq	$0x200,%%r15		\n\t"\
		"addq	$0x300,%%rcx		\n\t"\
		"addq	$0x400,%%rdx		\n\t"\
		"movq	%[__a3],%%rdi	/* s1p24r */\n\t"\
		"addq	$0x20,%%rdi	/* Out0 <- s1p25r */\n\t"\
		"movaps	    (%%rax),%%xmm0	\n\t"\
		"movaps	0x10(%%rax),%%xmm1	\n\t"\
		"movaps	    (%%r15),%%xmm2	\n\t"\
		"movaps	0x10(%%r15),%%xmm3	\n\t"\
		"movaps	    (%%rcx),%%xmm4	\n\t"\
		"movaps	0x10(%%rcx),%%xmm5	\n\t"\
		"movaps	    (%%rdx),%%xmm6	\n\t"\
		"movaps	0x10(%%rdx),%%xmm7	\n\t"\
		"subpd	%%xmm6,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm1		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm0,%%xmm6		\n\t"\
		"addpd	%%xmm1,%%xmm7		\n\t"\
		"subpd	%%xmm4,%%xmm2		\n\t"\
		"subpd	%%xmm5,%%xmm3		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm2,%%xmm4		\n\t"\
		"addpd	%%xmm3,%%xmm5		\n\t"\
		"/*movq	%[__cc1],%%rax		Need rax for a1, instead use that cc1 = r02 + 0x9f0 */\n\t"\
		"subpd	%%xmm4,%%xmm6		\n\t"\
		"subpd	%%xmm5,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm5		\n\t"\
		"addpd	    (%%rsi),%%xmm4	\n\t"\
		"addpd	0x10(%%rsi),%%xmm5	\n\t"\
		"movaps	%%xmm4,    (%%rdi)	\n\t"\
		"movaps	%%xmm5,0x10(%%rdi)	\n\t"\
		"mulpd	0xa00(%%rsi),%%xmm6	\n\t"\
		"mulpd	0xa00(%%rsi),%%xmm7	\n\t"\
		"subpd	     (%%rsi),%%xmm4	\n\t"\
		"subpd	0x010(%%rsi),%%xmm5	\n\t"\
		"mulpd	0x9f0(%%rsi),%%xmm4	\n\t"\
		"mulpd	0x9f0(%%rsi),%%xmm5	\n\t"\
		"addpd	     (%%rdi),%%xmm4	\n\t"\
		"addpd	0x010(%%rdi),%%xmm5	\n\t"\
		"subpd	%%xmm6,%%xmm4		\n\t"\
		"subpd	%%xmm7,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm6		\n\t"\
		"addpd	%%xmm5,%%xmm7		\n\t"\
		"movaps	%%xmm4,    (%%rsi)	\n\t"\
		"movaps	%%xmm5,0x10(%%rsi)	\n\t"\
		"movaps	%%xmm0,%%xmm4		\n\t"\
		"movaps	%%xmm1,%%xmm5		\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"mulpd	0xa10(%%rsi),%%xmm0	\n\t"\
		"mulpd	0xa10(%%rsi),%%xmm1	\n\t"\
		"mulpd	0xa20(%%rsi),%%xmm2	\n\t"\
		"mulpd	0xa20(%%rsi),%%xmm3	\n\t"\
		"mulpd	0xa30(%%rsi),%%xmm4	\n\t"\
		"mulpd	0xa30(%%rsi),%%xmm5	\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"movaps	    (%%rsi),%%xmm4	\n\t"\
		"movaps	0x10(%%rsi),%%xmm5	\n\t"\
		"movq	%[__a0],%%rax	/* Out1 <- s1p01r */\n\t"\
		"movq	%[__a1],%%rdx	/* Out4 <- s1p09r */\n\t"\
		"addq	$0x20,%%rax			\n\t"\
		"addq	$0x20,%%rdx			\n\t"\
		"subpd	%%xmm3,%%xmm6		\n\t"\
		"subpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"movaps	%%xmm6,    (%%rax)	\n\t"\
		"movaps	%%xmm7,0x10(%%rdx)	\n\t"\
		"addpd	%%xmm6,%%xmm3		\n\t"\
		"addpd	%%xmm7,%%xmm2		\n\t"\
		"movaps	%%xmm3,    (%%rdx)	\n\t"\
		"movaps	%%xmm2,0x10(%%rax)	\n\t"\
		"movq	%[__a2],%%r15	/* Out2 <- s1p17r */\n\t"\
		"movq	%[__a4],%%rcx	/* Out3 <- s1p33r */\n\t"\
		"addq	$0x20,%%r15			\n\t"\
		"addq	$0x20,%%rcx			\n\t"\
		"subpd	%%xmm1,%%xmm4		\n\t"\
		"subpd	%%xmm0,%%xmm5		\n\t"\
		"addpd	%%xmm1,%%xmm1		\n\t"\
		"addpd	%%xmm0,%%xmm0		\n\t"\
		"movaps	%%xmm4,    (%%r15)	\n\t"\
		"movaps	%%xmm5,0x10(%%rcx)	\n\t"\
		"addpd	%%xmm4,%%xmm1		\n\t"\
		"addpd	%%xmm5,%%xmm0		\n\t"\
		"movaps	%%xmm1,    (%%rcx)	\n\t"\
		"movaps	%%xmm0,0x10(%%r15)	\n\t"\
		"							\n\t"\
		"/* SSE2_RADIX_05_DFT_0TWIDDLE(r04,r20,r36,r52,r68,cc1,s1p10r,s1p26r,s1p02r,s1p18r,s1p34r) */\n\t"\
		"addq	$0x20,%%rsi	/* r04 */\n\t"\
		"movq	%%rsi,%%rax			\n\t"\
		"movq	%%rsi,%%r15			\n\t"\
		"movq	%%rsi,%%rcx			\n\t"\
		"movq	%%rsi,%%rdx			\n\t"\
		"addq	$0x100,%%rax		\n\t"\
		"addq	$0x200,%%r15		\n\t"\
		"addq	$0x300,%%rcx		\n\t"\
		"addq	$0x400,%%rdx		\n\t"\
		"movq	%[__a1],%%rdi	/* Out0 <- s1p10r */\n\t"\
		"addq	$0x40,%%rdi			\n\t"\
		"movaps	    (%%rax),%%xmm0	\n\t"\
		"movaps	0x10(%%rax),%%xmm1	\n\t"\
		"movaps	    (%%r15),%%xmm2	\n\t"\
		"movaps	0x10(%%r15),%%xmm3	\n\t"\
		"movaps	    (%%rcx),%%xmm4	\n\t"\
		"movaps	0x10(%%rcx),%%xmm5	\n\t"\
		"movaps	    (%%rdx),%%xmm6	\n\t"\
		"movaps	0x10(%%rdx),%%xmm7	\n\t"\
		"subpd	%%xmm6,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm1		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm0,%%xmm6		\n\t"\
		"addpd	%%xmm1,%%xmm7		\n\t"\
		"subpd	%%xmm4,%%xmm2		\n\t"\
		"subpd	%%xmm5,%%xmm3		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm2,%%xmm4		\n\t"\
		"addpd	%%xmm3,%%xmm5		\n\t"\
		"/*movq	%[__cc1],%%rax		Need rax for a1, instead use that cc1 = r04 + 0x9d0 */\n\t"\
		"subpd	%%xmm4,%%xmm6		\n\t"\
		"subpd	%%xmm5,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm5		\n\t"\
		"addpd	    (%%rsi),%%xmm4	\n\t"\
		"addpd	0x10(%%rsi),%%xmm5	\n\t"\
		"movaps	%%xmm4,    (%%rdi)	\n\t"\
		"movaps	%%xmm5,0x10(%%rdi)	\n\t"\
		"mulpd	0x9e0(%%rsi),%%xmm6	\n\t"\
		"mulpd	0x9e0(%%rsi),%%xmm7	\n\t"\
		"subpd	     (%%rsi),%%xmm4	\n\t"\
		"subpd	0x010(%%rsi),%%xmm5	\n\t"\
		"mulpd	0x9d0(%%rsi),%%xmm4	\n\t"\
		"mulpd	0x9d0(%%rsi),%%xmm5	\n\t"\
		"addpd	     (%%rdi),%%xmm4	\n\t"\
		"addpd	0x010(%%rdi),%%xmm5	\n\t"\
		"subpd	%%xmm6,%%xmm4		\n\t"\
		"subpd	%%xmm7,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm6		\n\t"\
		"addpd	%%xmm5,%%xmm7		\n\t"\
		"movaps	%%xmm4,    (%%rsi)	\n\t"\
		"movaps	%%xmm5,0x10(%%rsi)	\n\t"\
		"movaps	%%xmm0,%%xmm4		\n\t"\
		"movaps	%%xmm1,%%xmm5		\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"mulpd	0x9f0(%%rsi),%%xmm0	\n\t"\
		"mulpd	0x9f0(%%rsi),%%xmm1	\n\t"\
		"mulpd	0xa00(%%rsi),%%xmm2	\n\t"\
		"mulpd	0xa00(%%rsi),%%xmm3	\n\t"\
		"mulpd	0xa10(%%rsi),%%xmm4	\n\t"\
		"mulpd	0xa10(%%rsi),%%xmm5	\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"movaps	    (%%rsi),%%xmm4	\n\t"\
		"movaps	0x10(%%rsi),%%xmm5	\n\t"\
		"movq	%[__a3],%%rax	/* Out1 <- s1p26r */\n\t"\
		"movq	%[__a4],%%rdx	/* Out4 <- s1p34r */\n\t"\
		"addq	$0x40,%%rax			\n\t"\
		"addq	$0x40,%%rdx			\n\t"\
		"subpd	%%xmm3,%%xmm6		\n\t"\
		"subpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"movaps	%%xmm6,    (%%rax)	\n\t"\
		"movaps	%%xmm7,0x10(%%rdx)	\n\t"\
		"addpd	%%xmm6,%%xmm3		\n\t"\
		"addpd	%%xmm7,%%xmm2		\n\t"\
		"movaps	%%xmm3,    (%%rdx)	\n\t"\
		"movaps	%%xmm2,0x10(%%rax)	\n\t"\
		"movq	%[__a0],%%r15	/* Out2 <- s1p02r */\n\t"\
		"movq	%[__a2],%%rcx	/* Out3 <- s1p18r */\n\t"\
		"addq	$0x40,%%r15			\n\t"\
		"addq	$0x40,%%rcx			\n\t"\
		"subpd	%%xmm1,%%xmm4		\n\t"\
		"subpd	%%xmm0,%%xmm5		\n\t"\
		"addpd	%%xmm1,%%xmm1		\n\t"\
		"addpd	%%xmm0,%%xmm0		\n\t"\
		"movaps	%%xmm4,    (%%r15)	\n\t"\
		"movaps	%%xmm5,0x10(%%rcx)	\n\t"\
		"addpd	%%xmm4,%%xmm1		\n\t"\
		"addpd	%%xmm5,%%xmm0		\n\t"\
		"movaps	%%xmm1,    (%%rcx)	\n\t"\
		"movaps	%%xmm0,0x10(%%r15)	\n\t"\
		"							\n\t"\
		"/* SSE2_RADIX_05_DFT_0TWIDDLE(r06,r22,r38,r54,r70,cc1,s1p35r,s1p11r,s1p27r,s1p03r,s1p19r) */\n\t"\
		"addq	$0x20,%%rsi	/* r06 */\n\t"\
		"movq	%%rsi,%%rax			\n\t"\
		"movq	%%rsi,%%r15			\n\t"\
		"movq	%%rsi,%%rcx			\n\t"\
		"movq	%%rsi,%%rdx			\n\t"\
		"addq	$0x100,%%rax		\n\t"\
		"addq	$0x200,%%r15		\n\t"\
		"addq	$0x300,%%rcx		\n\t"\
		"addq	$0x400,%%rdx		\n\t"\
		"movq	%[__a4],%%rdi	/* Out0 <- s1p34r */\n\t"\
		"addq	$0x60,%%rdi			\n\t"\
		"movaps	    (%%rax),%%xmm0	\n\t"\
		"movaps	0x10(%%rax),%%xmm1	\n\t"\
		"movaps	    (%%r15),%%xmm2	\n\t"\
		"movaps	0x10(%%r15),%%xmm3	\n\t"\
		"movaps	    (%%rcx),%%xmm4	\n\t"\
		"movaps	0x10(%%rcx),%%xmm5	\n\t"\
		"movaps	    (%%rdx),%%xmm6	\n\t"\
		"movaps	0x10(%%rdx),%%xmm7	\n\t"\
		"subpd	%%xmm6,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm1		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm0,%%xmm6		\n\t"\
		"addpd	%%xmm1,%%xmm7		\n\t"\
		"subpd	%%xmm4,%%xmm2		\n\t"\
		"subpd	%%xmm5,%%xmm3		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm2,%%xmm4		\n\t"\
		"addpd	%%xmm3,%%xmm5		\n\t"\
		"/*movq	%[__cc1],%%rax		Need rax for a1, instead use that cc1 = r06 + 0x9b0 */\n\t"\
		"subpd	%%xmm4,%%xmm6		\n\t"\
		"subpd	%%xmm5,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm5		\n\t"\
		"addpd	    (%%rsi),%%xmm4	\n\t"\
		"addpd	0x10(%%rsi),%%xmm5	\n\t"\
		"movaps	%%xmm4,    (%%rdi)	\n\t"\
		"movaps	%%xmm5,0x10(%%rdi)	\n\t"\
		"mulpd	0x9c0(%%rsi),%%xmm6	\n\t"\
		"mulpd	0x9c0(%%rsi),%%xmm7	\n\t"\
		"subpd	     (%%rsi),%%xmm4	\n\t"\
		"subpd	0x010(%%rsi),%%xmm5	\n\t"\
		"mulpd	0x9b0(%%rsi),%%xmm4	\n\t"\
		"mulpd	0x9b0(%%rsi),%%xmm5	\n\t"\
		"addpd	     (%%rdi),%%xmm4	\n\t"\
		"addpd	0x010(%%rdi),%%xmm5	\n\t"\
		"subpd	%%xmm6,%%xmm4		\n\t"\
		"subpd	%%xmm7,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm6		\n\t"\
		"addpd	%%xmm5,%%xmm7		\n\t"\
		"movaps	%%xmm4,    (%%rsi)	\n\t"\
		"movaps	%%xmm5,0x10(%%rsi)	\n\t"\
		"movaps	%%xmm0,%%xmm4		\n\t"\
		"movaps	%%xmm1,%%xmm5		\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"mulpd	0x9d0(%%rsi),%%xmm0	\n\t"\
		"mulpd	0x9d0(%%rsi),%%xmm1	\n\t"\
		"mulpd	0x9e0(%%rsi),%%xmm2	\n\t"\
		"mulpd	0x9e0(%%rsi),%%xmm3	\n\t"\
		"mulpd	0x9f0(%%rsi),%%xmm4	\n\t"\
		"mulpd	0x9f0(%%rsi),%%xmm5	\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"movaps	    (%%rsi),%%xmm4	\n\t"\
		"movaps	0x10(%%rsi),%%xmm5	\n\t"\
		"movq	%[__a1],%%rax	/* Out1 <- s1p11r */\n\t"\
		"movq	%[__a2],%%rdx	/* Out4 <- s1p19r */\n\t"\
		"addq	$0x60,%%rax			\n\t"\
		"addq	$0x60,%%rdx			\n\t"\
		"subpd	%%xmm3,%%xmm6		\n\t"\
		"subpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"movaps	%%xmm6,    (%%rax)	\n\t"\
		"movaps	%%xmm7,0x10(%%rdx)	\n\t"\
		"addpd	%%xmm6,%%xmm3		\n\t"\
		"addpd	%%xmm7,%%xmm2		\n\t"\
		"movaps	%%xmm3,    (%%rdx)	\n\t"\
		"movaps	%%xmm2,0x10(%%rax)	\n\t"\
		"movq	%[__a3],%%r15	/* Out2 <- s1p27r */\n\t"\
		"movq	%[__a0],%%rcx	/* Out3 <- s1p03r */\n\t"\
		"addq	$0x60,%%r15			\n\t"\
		"addq	$0x60,%%rcx			\n\t"\
		"subpd	%%xmm1,%%xmm4		\n\t"\
		"subpd	%%xmm0,%%xmm5		\n\t"\
		"addpd	%%xmm1,%%xmm1		\n\t"\
		"addpd	%%xmm0,%%xmm0		\n\t"\
		"movaps	%%xmm4,    (%%r15)	\n\t"\
		"movaps	%%xmm5,0x10(%%rcx)	\n\t"\
		"addpd	%%xmm4,%%xmm1		\n\t"\
		"addpd	%%xmm5,%%xmm0		\n\t"\
		"movaps	%%xmm1,    (%%rcx)	\n\t"\
		"movaps	%%xmm0,0x10(%%r15)	\n\t"\
		"							\n\t"\
		"/* SSE2_RADIX_05_DFT_0TWIDDLE(r08,r24,r40,r56,r72,cc1,s1p20r,s1p36r,s1p12r,s1p28r,s1p04r) */\n\t"\
		"addq	$0x20,%%rsi	/* r08 */\n\t"\
		"movq	%%rsi,%%rax			\n\t"\
		"movq	%%rsi,%%r15			\n\t"\
		"movq	%%rsi,%%rcx			\n\t"\
		"movq	%%rsi,%%rdx			\n\t"\
		"addq	$0x100,%%rax		\n\t"\
		"addq	$0x200,%%r15		\n\t"\
		"addq	$0x300,%%rcx		\n\t"\
		"addq	$0x400,%%rdx		\n\t"\
		"movq	%[__a2],%%rdi	/* Out0 <- s1p20r */\n\t"\
		"addq	$0x80,%%rdi			\n\t"\
		"movaps	    (%%rax),%%xmm0	\n\t"\
		"movaps	0x10(%%rax),%%xmm1	\n\t"\
		"movaps	    (%%r15),%%xmm2	\n\t"\
		"movaps	0x10(%%r15),%%xmm3	\n\t"\
		"movaps	    (%%rcx),%%xmm4	\n\t"\
		"movaps	0x10(%%rcx),%%xmm5	\n\t"\
		"movaps	    (%%rdx),%%xmm6	\n\t"\
		"movaps	0x10(%%rdx),%%xmm7	\n\t"\
		"subpd	%%xmm6,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm1		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm0,%%xmm6		\n\t"\
		"addpd	%%xmm1,%%xmm7		\n\t"\
		"subpd	%%xmm4,%%xmm2		\n\t"\
		"subpd	%%xmm5,%%xmm3		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm2,%%xmm4		\n\t"\
		"addpd	%%xmm3,%%xmm5		\n\t"\
		"/*movq	%[__cc1],%%rax		Need rax for a1, instead use that cc1 = r08 + 0x990 */\n\t"\
		"subpd	%%xmm4,%%xmm6		\n\t"\
		"subpd	%%xmm5,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm5		\n\t"\
		"addpd	    (%%rsi),%%xmm4	\n\t"\
		"addpd	0x10(%%rsi),%%xmm5	\n\t"\
		"movaps	%%xmm4,    (%%rdi)	\n\t"\
		"movaps	%%xmm5,0x10(%%rdi)	\n\t"\
		"mulpd	0x9a0(%%rsi),%%xmm6	\n\t"\
		"mulpd	0x9a0(%%rsi),%%xmm7	\n\t"\
		"subpd	     (%%rsi),%%xmm4	\n\t"\
		"subpd	0x010(%%rsi),%%xmm5	\n\t"\
		"mulpd	0x990(%%rsi),%%xmm4	\n\t"\
		"mulpd	0x990(%%rsi),%%xmm5	\n\t"\
		"addpd	     (%%rdi),%%xmm4	\n\t"\
		"addpd	0x010(%%rdi),%%xmm5	\n\t"\
		"subpd	%%xmm6,%%xmm4		\n\t"\
		"subpd	%%xmm7,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm6		\n\t"\
		"addpd	%%xmm5,%%xmm7		\n\t"\
		"movaps	%%xmm4,    (%%rsi)	\n\t"\
		"movaps	%%xmm5,0x10(%%rsi)	\n\t"\
		"movaps	%%xmm0,%%xmm4		\n\t"\
		"movaps	%%xmm1,%%xmm5		\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"mulpd	0x9b0(%%rsi),%%xmm0	\n\t"\
		"mulpd	0x9b0(%%rsi),%%xmm1	\n\t"\
		"mulpd	0x9c0(%%rsi),%%xmm2	\n\t"\
		"mulpd	0x9c0(%%rsi),%%xmm3	\n\t"\
		"mulpd	0x9d0(%%rsi),%%xmm4	\n\t"\
		"mulpd	0x9d0(%%rsi),%%xmm5	\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"movaps	    (%%rsi),%%xmm4	\n\t"\
		"movaps	0x10(%%rsi),%%xmm5	\n\t"\
		"movq	%[__a4],%%rax	/* Out1 <- s1p36r */\n\t"\
		"movq	%[__a0],%%rdx	/* Out4 <- s1p04r */\n\t"\
		"addq	$0x80,%%rax			\n\t"\
		"addq	$0x80,%%rdx			\n\t"\
		"subpd	%%xmm3,%%xmm6		\n\t"\
		"subpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"movaps	%%xmm6,    (%%rax)	\n\t"\
		"movaps	%%xmm7,0x10(%%rdx)	\n\t"\
		"addpd	%%xmm6,%%xmm3		\n\t"\
		"addpd	%%xmm7,%%xmm2		\n\t"\
		"movaps	%%xmm3,    (%%rdx)	\n\t"\
		"movaps	%%xmm2,0x10(%%rax)	\n\t"\
		"movq	%[__a1],%%r15	/* Out2 <- s1p12r */\n\t"\
		"movq	%[__a3],%%rcx	/* Out3 <- s1p28r */\n\t"\
		"addq	$0x80,%%r15			\n\t"\
		"addq	$0x80,%%rcx			\n\t"\
		"subpd	%%xmm1,%%xmm4		\n\t"\
		"subpd	%%xmm0,%%xmm5		\n\t"\
		"addpd	%%xmm1,%%xmm1		\n\t"\
		"addpd	%%xmm0,%%xmm0		\n\t"\
		"movaps	%%xmm4,    (%%r15)	\n\t"\
		"movaps	%%xmm5,0x10(%%rcx)	\n\t"\
		"addpd	%%xmm4,%%xmm1		\n\t"\
		"addpd	%%xmm5,%%xmm0		\n\t"\
		"movaps	%%xmm1,    (%%rcx)	\n\t"\
		"movaps	%%xmm0,0x10(%%r15)	\n\t"\
		"							\n\t"\
		"/* SSE2_RADIX_05_DFT_0TWIDDLE(r10,r26,r42,r58,r74,cc1,s1p05r,s1p21r,s1p37r,s1p13r,s1p29r) */\n\t"\
		"addq	$0x20,%%rsi	/* r10 */\n\t"\
		"movq	%%rsi,%%rax			\n\t"\
		"movq	%%rsi,%%r15			\n\t"\
		"movq	%%rsi,%%rcx			\n\t"\
		"movq	%%rsi,%%rdx			\n\t"\
		"addq	$0x100,%%rax		\n\t"\
		"addq	$0x200,%%r15		\n\t"\
		"addq	$0x300,%%rcx		\n\t"\
		"addq	$0x400,%%rdx		\n\t"\
		"movq	%[__a0],%%rdi	/* Out0 <- s1p05r */\n\t"\
		"addq	$0xa0,%%rdi			\n\t"\
		"movaps	    (%%rax),%%xmm0	\n\t"\
		"movaps	0x10(%%rax),%%xmm1	\n\t"\
		"movaps	    (%%r15),%%xmm2	\n\t"\
		"movaps	0x10(%%r15),%%xmm3	\n\t"\
		"movaps	    (%%rcx),%%xmm4	\n\t"\
		"movaps	0x10(%%rcx),%%xmm5	\n\t"\
		"movaps	    (%%rdx),%%xmm6	\n\t"\
		"movaps	0x10(%%rdx),%%xmm7	\n\t"\
		"subpd	%%xmm6,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm1		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm0,%%xmm6		\n\t"\
		"addpd	%%xmm1,%%xmm7		\n\t"\
		"subpd	%%xmm4,%%xmm2		\n\t"\
		"subpd	%%xmm5,%%xmm3		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm2,%%xmm4		\n\t"\
		"addpd	%%xmm3,%%xmm5		\n\t"\
		"/*movq	%[__cc1],%%rax		Need rax for a1, instead use that cc1 = r08 + 0x970 */\n\t"\
		"subpd	%%xmm4,%%xmm6		\n\t"\
		"subpd	%%xmm5,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm5		\n\t"\
		"addpd	    (%%rsi),%%xmm4	\n\t"\
		"addpd	0x10(%%rsi),%%xmm5	\n\t"\
		"movaps	%%xmm4,    (%%rdi)	\n\t"\
		"movaps	%%xmm5,0x10(%%rdi)	\n\t"\
		"mulpd	0x980(%%rsi),%%xmm6	\n\t"\
		"mulpd	0x980(%%rsi),%%xmm7	\n\t"\
		"subpd	     (%%rsi),%%xmm4	\n\t"\
		"subpd	0x010(%%rsi),%%xmm5	\n\t"\
		"mulpd	0x970(%%rsi),%%xmm4	\n\t"\
		"mulpd	0x970(%%rsi),%%xmm5	\n\t"\
		"addpd	     (%%rdi),%%xmm4	\n\t"\
		"addpd	0x010(%%rdi),%%xmm5	\n\t"\
		"subpd	%%xmm6,%%xmm4		\n\t"\
		"subpd	%%xmm7,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm6		\n\t"\
		"addpd	%%xmm5,%%xmm7		\n\t"\
		"movaps	%%xmm4,    (%%rsi)	\n\t"\
		"movaps	%%xmm5,0x10(%%rsi)	\n\t"\
		"movaps	%%xmm0,%%xmm4		\n\t"\
		"movaps	%%xmm1,%%xmm5		\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"mulpd	0x990(%%rsi),%%xmm0	\n\t"\
		"mulpd	0x990(%%rsi),%%xmm1	\n\t"\
		"mulpd	0x9a0(%%rsi),%%xmm2	\n\t"\
		"mulpd	0x9a0(%%rsi),%%xmm3	\n\t"\
		"mulpd	0x9b0(%%rsi),%%xmm4	\n\t"\
		"mulpd	0x9b0(%%rsi),%%xmm5	\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"movaps	    (%%rsi),%%xmm4	\n\t"\
		"movaps	0x10(%%rsi),%%xmm5	\n\t"\
		"movq	%[__a2],%%rax	/* Out1 <- s1p21r */\n\t"\
		"movq	%[__a3],%%rdx	/* Out4 <- s1p29r */\n\t"\
		"addq	$0xa0,%%rax			\n\t"\
		"addq	$0xa0,%%rdx			\n\t"\
		"subpd	%%xmm3,%%xmm6		\n\t"\
		"subpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"movaps	%%xmm6,    (%%rax)	\n\t"\
		"movaps	%%xmm7,0x10(%%rdx)	\n\t"\
		"addpd	%%xmm6,%%xmm3		\n\t"\
		"addpd	%%xmm7,%%xmm2		\n\t"\
		"movaps	%%xmm3,    (%%rdx)	\n\t"\
		"movaps	%%xmm2,0x10(%%rax)	\n\t"\
		"movq	%[__a4],%%r15	/* Out2 <- s1p37r */\n\t"\
		"movq	%[__a1],%%rcx	/* Out3 <- s1p13r */\n\t"\
		"addq	$0xa0,%%r15			\n\t"\
		"addq	$0xa0,%%rcx			\n\t"\
		"subpd	%%xmm1,%%xmm4		\n\t"\
		"subpd	%%xmm0,%%xmm5		\n\t"\
		"addpd	%%xmm1,%%xmm1		\n\t"\
		"addpd	%%xmm0,%%xmm0		\n\t"\
		"movaps	%%xmm4,    (%%r15)	\n\t"\
		"movaps	%%xmm5,0x10(%%rcx)	\n\t"\
		"addpd	%%xmm4,%%xmm1		\n\t"\
		"addpd	%%xmm5,%%xmm0		\n\t"\
		"movaps	%%xmm1,    (%%rcx)	\n\t"\
		"movaps	%%xmm0,0x10(%%r15)	\n\t"\
		"							\n\t"\
		"/* SSE2_RADIX_05_DFT_0TWIDDLE(r12,r28,r44,r60,r76,cc1,s1p30r,s1p06r,s1p22r,s1p38r,s1p14r) */\n\t"\
		"addq	$0x20,%%rsi	/* r12 */\n\t"\
		"movq	%%rsi,%%rax			\n\t"\
		"movq	%%rsi,%%r15			\n\t"\
		"movq	%%rsi,%%rcx			\n\t"\
		"movq	%%rsi,%%rdx			\n\t"\
		"addq	$0x100,%%rax		\n\t"\
		"addq	$0x200,%%r15		\n\t"\
		"addq	$0x300,%%rcx		\n\t"\
		"addq	$0x400,%%rdx		\n\t"\
		"movq	%[__a3],%%rdi	/* Out0 <- s1p30r */\n\t"\
		"addq	$0xc0,%%rdi			\n\t"\
		"movaps	    (%%rax),%%xmm0	\n\t"\
		"movaps	0x10(%%rax),%%xmm1	\n\t"\
		"movaps	    (%%r15),%%xmm2	\n\t"\
		"movaps	0x10(%%r15),%%xmm3	\n\t"\
		"movaps	    (%%rcx),%%xmm4	\n\t"\
		"movaps	0x10(%%rcx),%%xmm5	\n\t"\
		"movaps	    (%%rdx),%%xmm6	\n\t"\
		"movaps	0x10(%%rdx),%%xmm7	\n\t"\
		"subpd	%%xmm6,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm1		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm0,%%xmm6		\n\t"\
		"addpd	%%xmm1,%%xmm7		\n\t"\
		"subpd	%%xmm4,%%xmm2		\n\t"\
		"subpd	%%xmm5,%%xmm3		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm2,%%xmm4		\n\t"\
		"addpd	%%xmm3,%%xmm5		\n\t"\
		"/*movq	%[__cc1],%%rax		Need rax for a1, instead use that cc1 = r08 + 0x950 */\n\t"\
		"subpd	%%xmm4,%%xmm6		\n\t"\
		"subpd	%%xmm5,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm5		\n\t"\
		"addpd	    (%%rsi),%%xmm4	\n\t"\
		"addpd	0x10(%%rsi),%%xmm5	\n\t"\
		"movaps	%%xmm4,    (%%rdi)	\n\t"\
		"movaps	%%xmm5,0x10(%%rdi)	\n\t"\
		"mulpd	0x960(%%rsi),%%xmm6	\n\t"\
		"mulpd	0x960(%%rsi),%%xmm7	\n\t"\
		"subpd	     (%%rsi),%%xmm4	\n\t"\
		"subpd	0x010(%%rsi),%%xmm5	\n\t"\
		"mulpd	0x950(%%rsi),%%xmm4	\n\t"\
		"mulpd	0x950(%%rsi),%%xmm5	\n\t"\
		"addpd	     (%%rdi),%%xmm4	\n\t"\
		"addpd	0x010(%%rdi),%%xmm5	\n\t"\
		"subpd	%%xmm6,%%xmm4		\n\t"\
		"subpd	%%xmm7,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm6		\n\t"\
		"addpd	%%xmm5,%%xmm7		\n\t"\
		"movaps	%%xmm4,    (%%rsi)	\n\t"\
		"movaps	%%xmm5,0x10(%%rsi)	\n\t"\
		"movaps	%%xmm0,%%xmm4		\n\t"\
		"movaps	%%xmm1,%%xmm5		\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"mulpd	0x970(%%rsi),%%xmm0	\n\t"\
		"mulpd	0x970(%%rsi),%%xmm1	\n\t"\
		"mulpd	0x980(%%rsi),%%xmm2	\n\t"\
		"mulpd	0x980(%%rsi),%%xmm3	\n\t"\
		"mulpd	0x990(%%rsi),%%xmm4	\n\t"\
		"mulpd	0x990(%%rsi),%%xmm5	\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"movaps	    (%%rsi),%%xmm4	\n\t"\
		"movaps	0x10(%%rsi),%%xmm5	\n\t"\
		"movq	%[__a0],%%rax	/* Out1 <- s1p06r */\n\t"\
		"movq	%[__a1],%%rdx	/* Out4 <- s1p14r */\n\t"\
		"addq	$0xc0,%%rax			\n\t"\
		"addq	$0xc0,%%rdx			\n\t"\
		"subpd	%%xmm3,%%xmm6		\n\t"\
		"subpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"movaps	%%xmm6,    (%%rax)	\n\t"\
		"movaps	%%xmm7,0x10(%%rdx)	\n\t"\
		"addpd	%%xmm6,%%xmm3		\n\t"\
		"addpd	%%xmm7,%%xmm2		\n\t"\
		"movaps	%%xmm3,    (%%rdx)	\n\t"\
		"movaps	%%xmm2,0x10(%%rax)	\n\t"\
		"movq	%[__a2],%%r15	/* Out2 <- s1p22r */\n\t"\
		"movq	%[__a4],%%rcx	/* Out3 <- s1p38r */\n\t"\
		"addq	$0xc0,%%r15			\n\t"\
		"addq	$0xc0,%%rcx			\n\t"\
		"subpd	%%xmm1,%%xmm4		\n\t"\
		"subpd	%%xmm0,%%xmm5		\n\t"\
		"addpd	%%xmm1,%%xmm1		\n\t"\
		"addpd	%%xmm0,%%xmm0		\n\t"\
		"movaps	%%xmm4,    (%%r15)	\n\t"\
		"movaps	%%xmm5,0x10(%%rcx)	\n\t"\
		"addpd	%%xmm4,%%xmm1		\n\t"\
		"addpd	%%xmm5,%%xmm0		\n\t"\
		"movaps	%%xmm1,    (%%rcx)	\n\t"\
		"movaps	%%xmm0,0x10(%%r15)	\n\t"\
		"							\n\t"\
		"/* SSE2_RADIX_05_DFT_0TWIDDLE(r14,r30,r46,r62,r78,cc1,s1p15r,s1p31r,s1p07r,s1p23r,s1p39r) */\n\t"\
		"addq	$0x20,%%rsi	/* r14 */\n\t"\
		"movq	%%rsi,%%rax			\n\t"\
		"movq	%%rsi,%%r15			\n\t"\
		"movq	%%rsi,%%rcx			\n\t"\
		"movq	%%rsi,%%rdx			\n\t"\
		"addq	$0x100,%%rax		\n\t"\
		"addq	$0x200,%%r15		\n\t"\
		"addq	$0x300,%%rcx		\n\t"\
		"addq	$0x400,%%rdx		\n\t"\
		"movq	%[__a1],%%rdi	/* Out0 <- s1p15r */\n\t"\
		"addq	$0xe0,%%rdi			\n\t"\
		"movaps	    (%%rax),%%xmm0	\n\t"\
		"movaps	0x10(%%rax),%%xmm1	\n\t"\
		"movaps	    (%%r15),%%xmm2	\n\t"\
		"movaps	0x10(%%r15),%%xmm3	\n\t"\
		"movaps	    (%%rcx),%%xmm4	\n\t"\
		"movaps	0x10(%%rcx),%%xmm5	\n\t"\
		"movaps	    (%%rdx),%%xmm6	\n\t"\
		"movaps	0x10(%%rdx),%%xmm7	\n\t"\
		"subpd	%%xmm6,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm1		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm0,%%xmm6		\n\t"\
		"addpd	%%xmm1,%%xmm7		\n\t"\
		"subpd	%%xmm4,%%xmm2		\n\t"\
		"subpd	%%xmm5,%%xmm3		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm2,%%xmm4		\n\t"\
		"addpd	%%xmm3,%%xmm5		\n\t"\
		"/*movq	%[__cc1],%%rax		Need rax for a1, instead use that cc1 = r08 + 0x990 */\n\t"\
		"subpd	%%xmm4,%%xmm6		\n\t"\
		"subpd	%%xmm5,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm5		\n\t"\
		"addpd	    (%%rsi),%%xmm4	\n\t"\
		"addpd	0x10(%%rsi),%%xmm5	\n\t"\
		"movaps	%%xmm4,    (%%rdi)	\n\t"\
		"movaps	%%xmm5,0x10(%%rdi)	\n\t"\
		"mulpd	0x940(%%rsi),%%xmm6	\n\t"\
		"mulpd	0x940(%%rsi),%%xmm7	\n\t"\
		"subpd	     (%%rsi),%%xmm4	\n\t"\
		"subpd	0x010(%%rsi),%%xmm5	\n\t"\
		"mulpd	0x930(%%rsi),%%xmm4	\n\t"\
		"mulpd	0x930(%%rsi),%%xmm5	\n\t"\
		"addpd	     (%%rdi),%%xmm4	\n\t"\
		"addpd	0x010(%%rdi),%%xmm5	\n\t"\
		"subpd	%%xmm6,%%xmm4		\n\t"\
		"subpd	%%xmm7,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm6		\n\t"\
		"addpd	%%xmm5,%%xmm7		\n\t"\
		"movaps	%%xmm4,    (%%rsi)	\n\t"\
		"movaps	%%xmm5,0x10(%%rsi)	\n\t"\
		"movaps	%%xmm0,%%xmm4		\n\t"\
		"movaps	%%xmm1,%%xmm5		\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"mulpd	0x950(%%rsi),%%xmm0	\n\t"\
		"mulpd	0x950(%%rsi),%%xmm1	\n\t"\
		"mulpd	0x960(%%rsi),%%xmm2	\n\t"\
		"mulpd	0x960(%%rsi),%%xmm3	\n\t"\
		"mulpd	0x970(%%rsi),%%xmm4	\n\t"\
		"mulpd	0x970(%%rsi),%%xmm5	\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"movaps	    (%%rsi),%%xmm4	\n\t"\
		"movaps	0x10(%%rsi),%%xmm5	\n\t"\
		"movq	%[__a3],%%rax	/* Out1 <- s1p31r */\n\t"\
		"movq	%[__a4],%%rdx	/* Out4 <- s1p39r */\n\t"\
		"addq	$0xe0,%%rax			\n\t"\
		"addq	$0xe0,%%rdx			\n\t"\
		"subpd	%%xmm3,%%xmm6		\n\t"\
		"subpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"movaps	%%xmm6,    (%%rax)	\n\t"\
		"movaps	%%xmm7,0x10(%%rdx)	\n\t"\
		"addpd	%%xmm6,%%xmm3		\n\t"\
		"addpd	%%xmm7,%%xmm2		\n\t"\
		"movaps	%%xmm3,    (%%rdx)	\n\t"\
		"movaps	%%xmm2,0x10(%%rax)	\n\t"\
		"movq	%[__a0],%%r15	/* Out2 <- s1p07r */\n\t"\
		"movq	%[__a2],%%rcx	/* Out3 <- s1p23r */\n\t"\
		"addq	$0xe0,%%r15			\n\t"\
		"addq	$0xe0,%%rcx			\n\t"\
		"subpd	%%xmm1,%%xmm4		\n\t"\
		"subpd	%%xmm0,%%xmm5		\n\t"\
		"addpd	%%xmm1,%%xmm1		\n\t"\
		"addpd	%%xmm0,%%xmm0		\n\t"\
		"movaps	%%xmm4,    (%%r15)	\n\t"\
		"movaps	%%xmm5,0x10(%%rcx)	\n\t"\
		"addpd	%%xmm4,%%xmm1		\n\t"\
		"addpd	%%xmm5,%%xmm0		\n\t"\
		"movaps	%%xmm1,    (%%rcx)	\n\t"\
		"movaps	%%xmm0,0x10(%%r15)	\n\t"\
		"							\n\t"\
		:					/* outputs: none */\
		: [__add] "m" (Xadd)	/* All inputs from memory addresses here */\
		 ,[__p01] "m" (Xp01)\
		 ,[__p02] "m" (Xp02)\
		 ,[__p03] "m" (Xp03)\
		 ,[__p04] "m" (Xp04)\
		 ,[__p08] "m" (Xp08)\
		 ,[__p16] "m" (Xp16)\
		 ,[__p24] "m" (Xp24)\
		 ,[__p32] "m" (Xp32)\
		 ,[__r00] "m" (Xr00)\
		 ,[__cc1] "m" (Xcc1)\
		 ,[__a0] "m" (Xa0)\
		 ,[__a1] "m" (Xa1)\
		 ,[__a2] "m" (Xa2)\
		 ,[__a3] "m" (Xa3)\
		 ,[__a4] "m" (Xa4)\
		: "rax","r15","rcx","rdx","rdi","rsi","r10","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
	);\
	}

	#define	SSE2_RADIX40_DIF_NOTWIDDLE(Xadd,Xp01,Xp02,Xp03,Xp04,Xp05,Xp06,Xp07,Xp08,Xp16,Xp24,Xp32,Xr00,Xcc1,Xa0,Xa1,Xa2,Xa3,Xa4)\
	{\
	__asm__ volatile (\
		"/* SSE2_RADIX_05_DFT_0TWIDDLE(s1p00r,s1p32r,s1p24r,s1p16r,s1p08r,cc1,r00,r16,r32,r48,r64) */\n\t"\
		"movq	%[__a0],%%rsi	/* In0 <- s1p00r */	\n\t"\
		"movq	%[__r00],%%rdi	/* Out0 <- r00 */	\n\t"\
		"movaps	0x400(%%rsi),%%xmm0	/* In1 <- s1p32r */\n\t"\
		"movaps	0x410(%%rsi),%%xmm1	\n\t"\
		"movq	%%rdi,%%rax			\n\t"\
		"movaps	0x300(%%rsi),%%xmm2	/* In2 <- s1p24r */\n\t"\
		"movaps	0x310(%%rsi),%%xmm3	\n\t"\
		"movq	%%rdi,%%r15			\n\t"\
		"movaps	0x200(%%rsi),%%xmm4	/* In3 <- s1p16r */\n\t"\
		"movaps	0x210(%%rsi),%%xmm5	\n\t"\
		"movq	%%rdi,%%rcx			\n\t"\
		"movaps	0x100(%%rsi),%%xmm6	/* In4 <- s1p08r */\n\t"\
		"movaps	0x110(%%rsi),%%xmm7	\n\t"\
		"movq	%%rdi,%%rdx			\n\t"\
		"subpd	%%xmm6,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm1		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm0,%%xmm6		\n\t"\
		"addpd	%%xmm1,%%xmm7		\n\t"\
		"subpd	%%xmm4,%%xmm2		\n\t"\
		"subpd	%%xmm5,%%xmm3		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm2,%%xmm4		\n\t"\
		"addpd	%%xmm3,%%xmm5		\n\t"\
		"/*movq	%[__cc1],%%rax		Need rax for a1, instead use that cc1 = r00 + 0xa10 */\n\t"\
		"subpd	%%xmm4,%%xmm6		\n\t"\
		"subpd	%%xmm5,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm5		\n\t"\
		"addpd	    (%%rsi),%%xmm4	\n\t"\
		"addpd	0x10(%%rsi),%%xmm5	\n\t"\
		"movaps	%%xmm4,    (%%rdi)	\n\t"\
		"movaps	%%xmm5,0x10(%%rdi)	\n\t"\
		"mulpd	0xa20(%%rdi),%%xmm6	\n\t"\
		"mulpd	0xa20(%%rdi),%%xmm7	\n\t"\
		"subpd	     (%%rsi),%%xmm4	\n\t"\
		"subpd	0x010(%%rsi),%%xmm5	\n\t"\
		"mulpd	0xa10(%%rdi),%%xmm4	\n\t"\
		"mulpd	0xa10(%%rdi),%%xmm5	\n\t"\
		"addpd	     (%%rdi),%%xmm4	\n\t"\
		"addpd	0x010(%%rdi),%%xmm5	\n\t"\
		"subpd	%%xmm6,%%xmm4		\n\t"\
		"subpd	%%xmm7,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm6		\n\t"\
		"addpd	%%xmm5,%%xmm7		\n\t"\
		"movaps	%%xmm4,    (%%rsi)	\n\t"\
		"movaps	%%xmm5,0x10(%%rsi)	\n\t"\
		"movaps	%%xmm0,%%xmm4		\n\t"\
		"movaps	%%xmm1,%%xmm5		\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"mulpd	0xa30(%%rdi),%%xmm0	\n\t"\
		"mulpd	0xa30(%%rdi),%%xmm1	\n\t"\
		"mulpd	0xa40(%%rdi),%%xmm2	\n\t"\
		"mulpd	0xa40(%%rdi),%%xmm3	\n\t"\
		"mulpd	0xa50(%%rdi),%%xmm4	\n\t"\
		"mulpd	0xa50(%%rdi),%%xmm5	\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"movaps	    (%%rsi),%%xmm4	\n\t"\
		"movaps	0x10(%%rsi),%%xmm5	\n\t"\
		"addq	$0x100,%%rax	/* Out1 <- r16 */\n\t"\
		"addq	$0x400,%%rdx	/* Out4 <- r64 */\n\t"\
		"subpd	%%xmm3,%%xmm6		\n\t"\
		"subpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"movaps	%%xmm6,    (%%rax)	\n\t"\
		"movaps	%%xmm7,0x10(%%rdx)	\n\t"\
		"addpd	%%xmm6,%%xmm3		\n\t"\
		"addpd	%%xmm7,%%xmm2		\n\t"\
		"movaps	%%xmm3,    (%%rdx)	\n\t"\
		"movaps	%%xmm2,0x10(%%rax)	\n\t"\
		"addq	$0x200,%%r15	/* Out2 <- r32 */\n\t"\
		"addq	$0x300,%%rcx	/* Out3 <- r48 */\n\t"\
		"subpd	%%xmm1,%%xmm4		\n\t"\
		"subpd	%%xmm0,%%xmm5		\n\t"\
		"addpd	%%xmm1,%%xmm1		\n\t"\
		"addpd	%%xmm0,%%xmm0		\n\t"\
		"movaps	%%xmm4,    (%%r15)	\n\t"\
		"movaps	%%xmm5,0x10(%%rcx)	\n\t"\
		"addpd	%%xmm4,%%xmm1		\n\t"\
		"addpd	%%xmm5,%%xmm0		\n\t"\
		"movaps	%%xmm1,    (%%rcx)	\n\t"\
		"movaps	%%xmm0,0x10(%%r15)	\n\t"\
		"							\n\t"\
		"/* SSE2_RADIX_05_DFT_0TWIDDLE(s1p35r,s1p27r,s1p19r,s1p11r,s1p03r,cc1,r02,r18,r34,r50,r66) */\n\t"\
		"addq	$0x460,%%rsi	/* In0 <- s1p35r */	\n\t"\
		"addq	$0x20,%%rdi		/* Out0 <- r02 */	\n\t"\
		"movaps	-0x100(%%rsi),%%xmm0	/* In1 <- s1p27r */\n\t"\
		"movaps	-0x0f0(%%rsi),%%xmm1	\n\t"\
		"addq	$0x20,%%rax			\n\t"\
		"movaps	-0x200(%%rsi),%%xmm2	/* In2 <- s1p19r */\n\t"\
		"movaps	-0x1f0(%%rsi),%%xmm3	\n\t"\
		"addq	$0x20,%%r15			\n\t"\
		"movaps	-0x300(%%rsi),%%xmm4	/* In3 <- s1p11r */\n\t"\
		"movaps	-0x2f0(%%rsi),%%xmm5	\n\t"\
		"addq	$0x20,%%rcx			\n\t"\
		"movaps	-0x400(%%rsi),%%xmm6	/* In4 <- s1p03r */\n\t"\
		"movaps	-0x3f0(%%rsi),%%xmm7	\n\t"\
		"addq	$0x20,%%rdx			\n\t"\
		"subpd	%%xmm6,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm1		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm0,%%xmm6		\n\t"\
		"addpd	%%xmm1,%%xmm7		\n\t"\
		"subpd	%%xmm4,%%xmm2		\n\t"\
		"subpd	%%xmm5,%%xmm3		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm2,%%xmm4		\n\t"\
		"addpd	%%xmm3,%%xmm5		\n\t"\
		"/*movq	%[__cc1],%%rax		Need rax for a1, instead use that cc1 = r02 + 0x9f0 */\n\t"\
		"subpd	%%xmm4,%%xmm6		\n\t"\
		"subpd	%%xmm5,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm5		\n\t"\
		"addpd	    (%%rsi),%%xmm4	\n\t"\
		"addpd	0x10(%%rsi),%%xmm5	\n\t"\
		"movaps	%%xmm4,    (%%rdi)	\n\t"\
		"movaps	%%xmm5,0x10(%%rdi)	\n\t"\
		"mulpd	0xa00(%%rdi),%%xmm6	\n\t"\
		"mulpd	0xa00(%%rdi),%%xmm7	\n\t"\
		"subpd	     (%%rsi),%%xmm4	\n\t"\
		"subpd	0x010(%%rsi),%%xmm5	\n\t"\
		"mulpd	0x9f0(%%rdi),%%xmm4	\n\t"\
		"mulpd	0x9f0(%%rdi),%%xmm5	\n\t"\
		"addpd	     (%%rdi),%%xmm4	\n\t"\
		"addpd	0x010(%%rdi),%%xmm5	\n\t"\
		"subpd	%%xmm6,%%xmm4		\n\t"\
		"subpd	%%xmm7,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm6		\n\t"\
		"addpd	%%xmm5,%%xmm7		\n\t"\
		"movaps	%%xmm4,    (%%rsi)	\n\t"\
		"movaps	%%xmm5,0x10(%%rsi)	\n\t"\
		"movaps	%%xmm0,%%xmm4		\n\t"\
		"movaps	%%xmm1,%%xmm5		\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"mulpd	0xa10(%%rdi),%%xmm0	\n\t"\
		"mulpd	0xa10(%%rdi),%%xmm1	\n\t"\
		"mulpd	0xa20(%%rdi),%%xmm2	\n\t"\
		"mulpd	0xa20(%%rdi),%%xmm3	\n\t"\
		"mulpd	0xa30(%%rdi),%%xmm4	\n\t"\
		"mulpd	0xa30(%%rdi),%%xmm5	\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"movaps	    (%%rsi),%%xmm4	\n\t"\
		"movaps	0x10(%%rsi),%%xmm5	\n\t"\
		"subpd	%%xmm3,%%xmm6		\n\t"\
		"subpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"movaps	%%xmm6,    (%%rax)	\n\t"\
		"movaps	%%xmm7,0x10(%%rdx)	\n\t"\
		"addpd	%%xmm6,%%xmm3		\n\t"\
		"addpd	%%xmm7,%%xmm2		\n\t"\
		"movaps	%%xmm3,    (%%rdx)	\n\t"\
		"movaps	%%xmm2,0x10(%%rax)	\n\t"\
		"subpd	%%xmm1,%%xmm4		\n\t"\
		"subpd	%%xmm0,%%xmm5		\n\t"\
		"addpd	%%xmm1,%%xmm1		\n\t"\
		"addpd	%%xmm0,%%xmm0		\n\t"\
		"movaps	%%xmm4,    (%%r15)	\n\t"\
		"movaps	%%xmm5,0x10(%%rcx)	\n\t"\
		"addpd	%%xmm4,%%xmm1		\n\t"\
		"addpd	%%xmm5,%%xmm0		\n\t"\
		"movaps	%%xmm1,    (%%rcx)	\n\t"\
		"movaps	%%xmm0,0x10(%%r15)	\n\t"\
		"							\n\t"\
		"/* SSE2_RADIX_05_DFT_0TWIDDLE(s1p30r,s1p22r,s1p14r,s1p06r,s1p38r,cc1,r04,r20,r36,r52,r68) */\n\t"\
		"subq	$0x0a0,%%rsi	/* In0 <- s1p30r */	\n\t"\
		"addq	$0x20,%%rdi		/* Out0 <- r04 */	\n\t"\
		"movaps	-0x100(%%rsi),%%xmm0	/* In1 <- s1p22r */\n\t"\
		"movaps	-0x0f0(%%rsi),%%xmm1	\n\t"\
		"addq	$0x20,%%rax			\n\t"\
		"movaps	-0x200(%%rsi),%%xmm2	/* In2 <- s1p14r */\n\t"\
		"movaps	-0x1f0(%%rsi),%%xmm3	\n\t"\
		"addq	$0x20,%%r15			\n\t"\
		"movaps	-0x300(%%rsi),%%xmm4	/* In3 <- s1p06r */\n\t"\
		"movaps	-0x2f0(%%rsi),%%xmm5	\n\t"\
		"addq	$0x20,%%rcx			\n\t"\
		"movaps	 0x100(%%rsi),%%xmm6	/* In4 <- s1p38r */\n\t"\
		"movaps	 0x110(%%rsi),%%xmm7	\n\t"\
		"addq	$0x20,%%rdx			\n\t"\
		"subpd	%%xmm6,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm1		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm0,%%xmm6		\n\t"\
		"addpd	%%xmm1,%%xmm7		\n\t"\
		"subpd	%%xmm4,%%xmm2		\n\t"\
		"subpd	%%xmm5,%%xmm3		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm2,%%xmm4		\n\t"\
		"addpd	%%xmm3,%%xmm5		\n\t"\
		"/*movq	%[__cc1],%%rax		Need rax for a1, instead use that cc1 = r04 + 0x9d0 */\n\t"\
		"subpd	%%xmm4,%%xmm6		\n\t"\
		"subpd	%%xmm5,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm5		\n\t"\
		"addpd	    (%%rsi),%%xmm4	\n\t"\
		"addpd	0x10(%%rsi),%%xmm5	\n\t"\
		"movaps	%%xmm4,    (%%rdi)	\n\t"\
		"movaps	%%xmm5,0x10(%%rdi)	\n\t"\
		"mulpd	0x9e0(%%rdi),%%xmm6	\n\t"\
		"mulpd	0x9e0(%%rdi),%%xmm7	\n\t"\
		"subpd	     (%%rsi),%%xmm4	\n\t"\
		"subpd	0x010(%%rsi),%%xmm5	\n\t"\
		"mulpd	0x9d0(%%rdi),%%xmm4	\n\t"\
		"mulpd	0x9d0(%%rdi),%%xmm5	\n\t"\
		"addpd	     (%%rdi),%%xmm4	\n\t"\
		"addpd	0x010(%%rdi),%%xmm5	\n\t"\
		"subpd	%%xmm6,%%xmm4		\n\t"\
		"subpd	%%xmm7,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm6		\n\t"\
		"addpd	%%xmm5,%%xmm7		\n\t"\
		"movaps	%%xmm4,    (%%rsi)	\n\t"\
		"movaps	%%xmm5,0x10(%%rsi)	\n\t"\
		"movaps	%%xmm0,%%xmm4		\n\t"\
		"movaps	%%xmm1,%%xmm5		\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"mulpd	0x9f0(%%rdi),%%xmm0	\n\t"\
		"mulpd	0x9f0(%%rdi),%%xmm1	\n\t"\
		"mulpd	0xa00(%%rdi),%%xmm2	\n\t"\
		"mulpd	0xa00(%%rdi),%%xmm3	\n\t"\
		"mulpd	0xa10(%%rdi),%%xmm4	\n\t"\
		"mulpd	0xa10(%%rdi),%%xmm5	\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"movaps	    (%%rsi),%%xmm4	\n\t"\
		"movaps	0x10(%%rsi),%%xmm5	\n\t"\
		"subpd	%%xmm3,%%xmm6		\n\t"\
		"subpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"movaps	%%xmm6,    (%%rax)	\n\t"\
		"movaps	%%xmm7,0x10(%%rdx)	\n\t"\
		"addpd	%%xmm6,%%xmm3		\n\t"\
		"addpd	%%xmm7,%%xmm2		\n\t"\
		"movaps	%%xmm3,    (%%rdx)	\n\t"\
		"movaps	%%xmm2,0x10(%%rax)	\n\t"\
		"subpd	%%xmm1,%%xmm4		\n\t"\
		"subpd	%%xmm0,%%xmm5		\n\t"\
		"addpd	%%xmm1,%%xmm1		\n\t"\
		"addpd	%%xmm0,%%xmm0		\n\t"\
		"movaps	%%xmm4,    (%%r15)	\n\t"\
		"movaps	%%xmm5,0x10(%%rcx)	\n\t"\
		"addpd	%%xmm4,%%xmm1		\n\t"\
		"addpd	%%xmm5,%%xmm0		\n\t"\
		"movaps	%%xmm1,    (%%rcx)	\n\t"\
		"movaps	%%xmm0,0x10(%%r15)	\n\t"\
		"							\n\t"\
		"/* SSE2_RADIX_05_DFT_0TWIDDLE(s1p25r,s1p17r,s1p09r,s1p01r,s1p33r,cc1,r06,r22,r38,r54,r70) */\n\t"\
		"subq	$0x0a0,%%rsi	/* In0 <- s1p25r */	\n\t"\
		"addq	$0x20,%%rdi		/* Out0 <- r06 */	\n\t"\
		"movaps	-0x100(%%rsi),%%xmm0	/* In1 <- s1p17r */\n\t"\
		"movaps	-0x0f0(%%rsi),%%xmm1	\n\t"\
		"addq	$0x20,%%rax			\n\t"\
		"movaps	-0x200(%%rsi),%%xmm2	/* In2 <- s1p09r */\n\t"\
		"movaps	-0x1f0(%%rsi),%%xmm3	\n\t"\
		"addq	$0x20,%%r15			\n\t"\
		"movaps	-0x300(%%rsi),%%xmm4	/* In3 <- s1p01r */\n\t"\
		"movaps	-0x2f0(%%rsi),%%xmm5	\n\t"\
		"addq	$0x20,%%rcx			\n\t"\
		"movaps	 0x100(%%rsi),%%xmm6	/* In4 <- s1p33r */\n\t"\
		"movaps	 0x110(%%rsi),%%xmm7	\n\t"\
		"addq	$0x20,%%rdx			\n\t"\
		"subpd	%%xmm6,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm1		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm0,%%xmm6		\n\t"\
		"addpd	%%xmm1,%%xmm7		\n\t"\
		"subpd	%%xmm4,%%xmm2		\n\t"\
		"subpd	%%xmm5,%%xmm3		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm2,%%xmm4		\n\t"\
		"addpd	%%xmm3,%%xmm5		\n\t"\
		"/*movq	%[__cc1],%%rax		Need rax for a1, instead use that cc1 = r06 + 0x9b0 */\n\t"\
		"subpd	%%xmm4,%%xmm6		\n\t"\
		"subpd	%%xmm5,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm5		\n\t"\
		"addpd	    (%%rsi),%%xmm4	\n\t"\
		"addpd	0x10(%%rsi),%%xmm5	\n\t"\
		"movaps	%%xmm4,    (%%rdi)	\n\t"\
		"movaps	%%xmm5,0x10(%%rdi)	\n\t"\
		"mulpd	0x9c0(%%rdi),%%xmm6	\n\t"\
		"mulpd	0x9c0(%%rdi),%%xmm7	\n\t"\
		"subpd	     (%%rsi),%%xmm4	\n\t"\
		"subpd	0x010(%%rsi),%%xmm5	\n\t"\
		"mulpd	0x9b0(%%rdi),%%xmm4	\n\t"\
		"mulpd	0x9b0(%%rdi),%%xmm5	\n\t"\
		"addpd	     (%%rdi),%%xmm4	\n\t"\
		"addpd	0x010(%%rdi),%%xmm5	\n\t"\
		"subpd	%%xmm6,%%xmm4		\n\t"\
		"subpd	%%xmm7,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm6		\n\t"\
		"addpd	%%xmm5,%%xmm7		\n\t"\
		"movaps	%%xmm4,    (%%rsi)	\n\t"\
		"movaps	%%xmm5,0x10(%%rsi)	\n\t"\
		"movaps	%%xmm0,%%xmm4		\n\t"\
		"movaps	%%xmm1,%%xmm5		\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"mulpd	0x9d0(%%rdi),%%xmm0	\n\t"\
		"mulpd	0x9d0(%%rdi),%%xmm1	\n\t"\
		"mulpd	0x9e0(%%rdi),%%xmm2	\n\t"\
		"mulpd	0x9e0(%%rdi),%%xmm3	\n\t"\
		"mulpd	0x9f0(%%rdi),%%xmm4	\n\t"\
		"mulpd	0x9f0(%%rdi),%%xmm5	\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"movaps	    (%%rsi),%%xmm4	\n\t"\
		"movaps	0x10(%%rsi),%%xmm5	\n\t"\
		"subpd	%%xmm3,%%xmm6		\n\t"\
		"subpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"movaps	%%xmm6,    (%%rax)	\n\t"\
		"movaps	%%xmm7,0x10(%%rdx)	\n\t"\
		"addpd	%%xmm6,%%xmm3		\n\t"\
		"addpd	%%xmm7,%%xmm2		\n\t"\
		"movaps	%%xmm3,    (%%rdx)	\n\t"\
		"movaps	%%xmm2,0x10(%%rax)	\n\t"\
		"subpd	%%xmm1,%%xmm4		\n\t"\
		"subpd	%%xmm0,%%xmm5		\n\t"\
		"addpd	%%xmm1,%%xmm1		\n\t"\
		"addpd	%%xmm0,%%xmm0		\n\t"\
		"movaps	%%xmm4,    (%%r15)	\n\t"\
		"movaps	%%xmm5,0x10(%%rcx)	\n\t"\
		"addpd	%%xmm4,%%xmm1		\n\t"\
		"addpd	%%xmm5,%%xmm0		\n\t"\
		"movaps	%%xmm1,    (%%rcx)	\n\t"\
		"movaps	%%xmm0,0x10(%%r15)	\n\t"\
		"							\n\t"\
		"/* SSE2_RADIX_05_DFT_0TWIDDLE(s1p20r,s1p12r,s1p04r,s1p36r,s1p28r,cc1,r08,r24,r40,r56,r72) */\n\t"\
		"subq	$0x0a0,%%rsi	/* In0 <- s1p20r */	\n\t"\
		"addq	$0x20,%%rdi		/* Out0 <- r08 */	\n\t"\
		"movaps	-0x100(%%rsi),%%xmm0	/* In1 <- s1p12r */\n\t"\
		"movaps	-0x0f0(%%rsi),%%xmm1	\n\t"\
		"addq	$0x20,%%rax			\n\t"\
		"movaps	-0x200(%%rsi),%%xmm2	/* In2 <- s1p04r */\n\t"\
		"movaps	-0x1f0(%%rsi),%%xmm3	\n\t"\
		"addq	$0x20,%%r15			\n\t"\
		"movaps	 0x200(%%rsi),%%xmm4	/* In3 <- s1p36r */\n\t"\
		"movaps	 0x210(%%rsi),%%xmm5	\n\t"\
		"addq	$0x20,%%rcx			\n\t"\
		"movaps	 0x100(%%rsi),%%xmm6	/* In4 <- s1p28r */\n\t"\
		"movaps	 0x110(%%rsi),%%xmm7	\n\t"\
		"addq	$0x20,%%rdx			\n\t"\
		"subpd	%%xmm6,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm1		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm0,%%xmm6		\n\t"\
		"addpd	%%xmm1,%%xmm7		\n\t"\
		"subpd	%%xmm4,%%xmm2		\n\t"\
		"subpd	%%xmm5,%%xmm3		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm2,%%xmm4		\n\t"\
		"addpd	%%xmm3,%%xmm5		\n\t"\
		"/*movq	%[__cc1],%%rax		Need rax for a1, instead use that cc1 = r08 + 0x990 */\n\t"\
		"subpd	%%xmm4,%%xmm6		\n\t"\
		"subpd	%%xmm5,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm5		\n\t"\
		"addpd	    (%%rsi),%%xmm4	\n\t"\
		"addpd	0x10(%%rsi),%%xmm5	\n\t"\
		"movaps	%%xmm4,    (%%rdi)	\n\t"\
		"movaps	%%xmm5,0x10(%%rdi)	\n\t"\
		"mulpd	0x9a0(%%rdi),%%xmm6	\n\t"\
		"mulpd	0x9a0(%%rdi),%%xmm7	\n\t"\
		"subpd	     (%%rsi),%%xmm4	\n\t"\
		"subpd	0x010(%%rsi),%%xmm5	\n\t"\
		"mulpd	0x990(%%rdi),%%xmm4	\n\t"\
		"mulpd	0x990(%%rdi),%%xmm5	\n\t"\
		"addpd	     (%%rdi),%%xmm4	\n\t"\
		"addpd	0x010(%%rdi),%%xmm5	\n\t"\
		"subpd	%%xmm6,%%xmm4		\n\t"\
		"subpd	%%xmm7,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm6		\n\t"\
		"addpd	%%xmm5,%%xmm7		\n\t"\
		"movaps	%%xmm4,    (%%rsi)	\n\t"\
		"movaps	%%xmm5,0x10(%%rsi)	\n\t"\
		"movaps	%%xmm0,%%xmm4		\n\t"\
		"movaps	%%xmm1,%%xmm5		\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"mulpd	0x9b0(%%rdi),%%xmm0	\n\t"\
		"mulpd	0x9b0(%%rdi),%%xmm1	\n\t"\
		"mulpd	0x9c0(%%rdi),%%xmm2	\n\t"\
		"mulpd	0x9c0(%%rdi),%%xmm3	\n\t"\
		"mulpd	0x9d0(%%rdi),%%xmm4	\n\t"\
		"mulpd	0x9d0(%%rdi),%%xmm5	\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"movaps	    (%%rsi),%%xmm4	\n\t"\
		"movaps	0x10(%%rsi),%%xmm5	\n\t"\
		"subpd	%%xmm3,%%xmm6		\n\t"\
		"subpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"movaps	%%xmm6,    (%%rax)	\n\t"\
		"movaps	%%xmm7,0x10(%%rdx)	\n\t"\
		"addpd	%%xmm6,%%xmm3		\n\t"\
		"addpd	%%xmm7,%%xmm2		\n\t"\
		"movaps	%%xmm3,    (%%rdx)	\n\t"\
		"movaps	%%xmm2,0x10(%%rax)	\n\t"\
		"subpd	%%xmm1,%%xmm4		\n\t"\
		"subpd	%%xmm0,%%xmm5		\n\t"\
		"addpd	%%xmm1,%%xmm1		\n\t"\
		"addpd	%%xmm0,%%xmm0		\n\t"\
		"movaps	%%xmm4,    (%%r15)	\n\t"\
		"movaps	%%xmm5,0x10(%%rcx)	\n\t"\
		"addpd	%%xmm4,%%xmm1		\n\t"\
		"addpd	%%xmm5,%%xmm0		\n\t"\
		"movaps	%%xmm1,    (%%rcx)	\n\t"\
		"movaps	%%xmm0,0x10(%%r15)	\n\t"\
		"							\n\t"\
		"/* SSE2_RADIX_05_DFT_0TWIDDLE(s1p15r,s1p07r,s1p39r,s1p31r,s1p23r,cc1,r10,r26,r42,r58,r74) */\n\t"\
		"subq	$0x0a0,%%rsi	/* In0 <- s1p15r */	\n\t"\
		"addq	$0x20,%%rdi		/* Out0 <- r10 */	\n\t"\
		"movaps	-0x100(%%rsi),%%xmm0	/* In1 <- s1p07r */\n\t"\
		"movaps	-0x0f0(%%rsi),%%xmm1	\n\t"\
		"addq	$0x20,%%rax			\n\t"\
		"movaps	 0x300(%%rsi),%%xmm2	/* In2 <- s1p39r */\n\t"\
		"movaps	 0x310(%%rsi),%%xmm3	\n\t"\
		"addq	$0x20,%%r15			\n\t"\
		"movaps	 0x200(%%rsi),%%xmm4	/* In3 <- s1p31r */\n\t"\
		"movaps	 0x210(%%rsi),%%xmm5	\n\t"\
		"addq	$0x20,%%rcx			\n\t"\
		"movaps	 0x100(%%rsi),%%xmm6	/* In4 <- s1p23r */\n\t"\
		"movaps	 0x110(%%rsi),%%xmm7	\n\t"\
		"addq	$0x20,%%rdx			\n\t"\
		"subpd	%%xmm6,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm1		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm0,%%xmm6		\n\t"\
		"addpd	%%xmm1,%%xmm7		\n\t"\
		"subpd	%%xmm4,%%xmm2		\n\t"\
		"subpd	%%xmm5,%%xmm3		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm2,%%xmm4		\n\t"\
		"addpd	%%xmm3,%%xmm5		\n\t"\
		"/*movq	%[__cc1],%%rax		Need rax for a1, instead use that cc1 = r10 + 0x970 */\n\t"\
		"subpd	%%xmm4,%%xmm6		\n\t"\
		"subpd	%%xmm5,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm5		\n\t"\
		"addpd	    (%%rsi),%%xmm4	\n\t"\
		"addpd	0x10(%%rsi),%%xmm5	\n\t"\
		"movaps	%%xmm4,    (%%rdi)	\n\t"\
		"movaps	%%xmm5,0x10(%%rdi)	\n\t"\
		"mulpd	0x980(%%rdi),%%xmm6	\n\t"\
		"mulpd	0x980(%%rdi),%%xmm7	\n\t"\
		"subpd	     (%%rsi),%%xmm4	\n\t"\
		"subpd	0x010(%%rsi),%%xmm5	\n\t"\
		"mulpd	0x970(%%rdi),%%xmm4	\n\t"\
		"mulpd	0x970(%%rdi),%%xmm5	\n\t"\
		"addpd	     (%%rdi),%%xmm4	\n\t"\
		"addpd	0x010(%%rdi),%%xmm5	\n\t"\
		"subpd	%%xmm6,%%xmm4		\n\t"\
		"subpd	%%xmm7,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm6		\n\t"\
		"addpd	%%xmm5,%%xmm7		\n\t"\
		"movaps	%%xmm4,    (%%rsi)	\n\t"\
		"movaps	%%xmm5,0x10(%%rsi)	\n\t"\
		"movaps	%%xmm0,%%xmm4		\n\t"\
		"movaps	%%xmm1,%%xmm5		\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"mulpd	0x990(%%rdi),%%xmm0	\n\t"\
		"mulpd	0x990(%%rdi),%%xmm1	\n\t"\
		"mulpd	0x9a0(%%rdi),%%xmm2	\n\t"\
		"mulpd	0x9a0(%%rdi),%%xmm3	\n\t"\
		"mulpd	0x9b0(%%rdi),%%xmm4	\n\t"\
		"mulpd	0x9b0(%%rdi),%%xmm5	\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"movaps	    (%%rsi),%%xmm4	\n\t"\
		"movaps	0x10(%%rsi),%%xmm5	\n\t"\
		"subpd	%%xmm3,%%xmm6		\n\t"\
		"subpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"movaps	%%xmm6,    (%%rax)	\n\t"\
		"movaps	%%xmm7,0x10(%%rdx)	\n\t"\
		"addpd	%%xmm6,%%xmm3		\n\t"\
		"addpd	%%xmm7,%%xmm2		\n\t"\
		"movaps	%%xmm3,    (%%rdx)	\n\t"\
		"movaps	%%xmm2,0x10(%%rax)	\n\t"\
		"subpd	%%xmm1,%%xmm4		\n\t"\
		"subpd	%%xmm0,%%xmm5		\n\t"\
		"addpd	%%xmm1,%%xmm1		\n\t"\
		"addpd	%%xmm0,%%xmm0		\n\t"\
		"movaps	%%xmm4,    (%%r15)	\n\t"\
		"movaps	%%xmm5,0x10(%%rcx)	\n\t"\
		"addpd	%%xmm4,%%xmm1		\n\t"\
		"addpd	%%xmm5,%%xmm0		\n\t"\
		"movaps	%%xmm1,    (%%rcx)	\n\t"\
		"movaps	%%xmm0,0x10(%%r15)	\n\t"\
		"							\n\t"\
		"/* SSE2_RADIX_05_DFT_0TWIDDLE(s1p10r,s1p02r,s1p34r,s1p26r,s1p18r,cc1,r12,r28,r44,r60,r76) */\n\t"\
		"subq	$0x0a0,%%rsi	/* In0 <- s1p10r */	\n\t"\
		"addq	$0x20,%%rdi		/* Out0 <- r12 */	\n\t"\
		"movaps	-0x100(%%rsi),%%xmm0	/* In1 <- s1p02r */\n\t"\
		"movaps	-0x0f0(%%rsi),%%xmm1	\n\t"\
		"addq	$0x20,%%rax			\n\t"\
		"movaps	 0x300(%%rsi),%%xmm2	/* In2 <- s1p34r */\n\t"\
		"movaps	 0x310(%%rsi),%%xmm3	\n\t"\
		"addq	$0x20,%%r15			\n\t"\
		"movaps	 0x200(%%rsi),%%xmm4	/* In3 <- s1p26r */\n\t"\
		"movaps	 0x210(%%rsi),%%xmm5	\n\t"\
		"addq	$0x20,%%rcx			\n\t"\
		"movaps	 0x100(%%rsi),%%xmm6	/* In4 <- s1p18r */\n\t"\
		"movaps	 0x110(%%rsi),%%xmm7	\n\t"\
		"addq	$0x20,%%rdx			\n\t"\
		"subpd	%%xmm6,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm1		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm0,%%xmm6		\n\t"\
		"addpd	%%xmm1,%%xmm7		\n\t"\
		"subpd	%%xmm4,%%xmm2		\n\t"\
		"subpd	%%xmm5,%%xmm3		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm2,%%xmm4		\n\t"\
		"addpd	%%xmm3,%%xmm5		\n\t"\
		"/*movq	%[__cc1],%%rax		Need rax for a1, instead use that cc1 = r12 + 0x950 */\n\t"\
		"subpd	%%xmm4,%%xmm6		\n\t"\
		"subpd	%%xmm5,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm5		\n\t"\
		"addpd	    (%%rsi),%%xmm4	\n\t"\
		"addpd	0x10(%%rsi),%%xmm5	\n\t"\
		"movaps	%%xmm4,    (%%rdi)	\n\t"\
		"movaps	%%xmm5,0x10(%%rdi)	\n\t"\
		"mulpd	0x960(%%rdi),%%xmm6	\n\t"\
		"mulpd	0x960(%%rdi),%%xmm7	\n\t"\
		"subpd	     (%%rsi),%%xmm4	\n\t"\
		"subpd	0x010(%%rsi),%%xmm5	\n\t"\
		"mulpd	0x950(%%rdi),%%xmm4	\n\t"\
		"mulpd	0x950(%%rdi),%%xmm5	\n\t"\
		"addpd	     (%%rdi),%%xmm4	\n\t"\
		"addpd	0x010(%%rdi),%%xmm5	\n\t"\
		"subpd	%%xmm6,%%xmm4		\n\t"\
		"subpd	%%xmm7,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm6		\n\t"\
		"addpd	%%xmm5,%%xmm7		\n\t"\
		"movaps	%%xmm4,    (%%rsi)	\n\t"\
		"movaps	%%xmm5,0x10(%%rsi)	\n\t"\
		"movaps	%%xmm0,%%xmm4		\n\t"\
		"movaps	%%xmm1,%%xmm5		\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"mulpd	0x970(%%rdi),%%xmm0	\n\t"\
		"mulpd	0x970(%%rdi),%%xmm1	\n\t"\
		"mulpd	0x980(%%rdi),%%xmm2	\n\t"\
		"mulpd	0x980(%%rdi),%%xmm3	\n\t"\
		"mulpd	0x990(%%rdi),%%xmm4	\n\t"\
		"mulpd	0x990(%%rdi),%%xmm5	\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"movaps	    (%%rsi),%%xmm4	\n\t"\
		"movaps	0x10(%%rsi),%%xmm5	\n\t"\
		"subpd	%%xmm3,%%xmm6		\n\t"\
		"subpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"movaps	%%xmm6,    (%%rax)	\n\t"\
		"movaps	%%xmm7,0x10(%%rdx)	\n\t"\
		"addpd	%%xmm6,%%xmm3		\n\t"\
		"addpd	%%xmm7,%%xmm2		\n\t"\
		"movaps	%%xmm3,    (%%rdx)	\n\t"\
		"movaps	%%xmm2,0x10(%%rax)	\n\t"\
		"subpd	%%xmm1,%%xmm4		\n\t"\
		"subpd	%%xmm0,%%xmm5		\n\t"\
		"addpd	%%xmm1,%%xmm1		\n\t"\
		"addpd	%%xmm0,%%xmm0		\n\t"\
		"movaps	%%xmm4,    (%%r15)	\n\t"\
		"movaps	%%xmm5,0x10(%%rcx)	\n\t"\
		"addpd	%%xmm4,%%xmm1		\n\t"\
		"addpd	%%xmm5,%%xmm0		\n\t"\
		"movaps	%%xmm1,    (%%rcx)	\n\t"\
		"movaps	%%xmm0,0x10(%%r15)	\n\t"\
		"							\n\t"\
		"/* SSE2_RADIX_05_DFT_0TWIDDLE(s1p05r,s1p37r,s1p29r,s1p21r,s1p13r,cc1,r14,r30,r46,r62,r78) */\n\t"\
		"subq	$0x0a0,%%rsi	/* In0 <- s1p05r */	\n\t"\
		"addq	$0x20,%%rdi		/* Out0 <- r14 */	\n\t"\
		"movaps	 0x400(%%rsi),%%xmm0	/* In1 <- s1p37r */\n\t"\
		"movaps	 0x410(%%rsi),%%xmm1	\n\t"\
		"addq	$0x20,%%rax			\n\t"\
		"movaps	 0x300(%%rsi),%%xmm2	/* In2 <- s1p29r */\n\t"\
		"movaps	 0x310(%%rsi),%%xmm3	\n\t"\
		"addq	$0x20,%%r15			\n\t"\
		"movaps	 0x200(%%rsi),%%xmm4	/* In3 <- s1p21r */\n\t"\
		"movaps	 0x210(%%rsi),%%xmm5	\n\t"\
		"addq	$0x20,%%rcx			\n\t"\
		"movaps	 0x100(%%rsi),%%xmm6	/* In4 <- s1p13r */\n\t"\
		"movaps	 0x110(%%rsi),%%xmm7	\n\t"\
		"addq	$0x20,%%rdx			\n\t"\
		"subpd	%%xmm6,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm1		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm0,%%xmm6		\n\t"\
		"addpd	%%xmm1,%%xmm7		\n\t"\
		"subpd	%%xmm4,%%xmm2		\n\t"\
		"subpd	%%xmm5,%%xmm3		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm2,%%xmm4		\n\t"\
		"addpd	%%xmm3,%%xmm5		\n\t"\
		"/*movq	%[__cc1],%%rax		Need rax for a1, instead use that cc1 = r14 + 0x930 */\n\t"\
		"subpd	%%xmm4,%%xmm6		\n\t"\
		"subpd	%%xmm5,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm5		\n\t"\
		"addpd	    (%%rsi),%%xmm4	\n\t"\
		"addpd	0x10(%%rsi),%%xmm5	\n\t"\
		"movaps	%%xmm4,    (%%rdi)	\n\t"\
		"movaps	%%xmm5,0x10(%%rdi)	\n\t"\
		"mulpd	0x940(%%rdi),%%xmm6	\n\t"\
		"mulpd	0x940(%%rdi),%%xmm7	\n\t"\
		"subpd	     (%%rsi),%%xmm4	\n\t"\
		"subpd	0x010(%%rsi),%%xmm5	\n\t"\
		"mulpd	0x930(%%rdi),%%xmm4	\n\t"\
		"mulpd	0x930(%%rdi),%%xmm5	\n\t"\
		"addpd	     (%%rdi),%%xmm4	\n\t"\
		"addpd	0x010(%%rdi),%%xmm5	\n\t"\
		"subpd	%%xmm6,%%xmm4		\n\t"\
		"subpd	%%xmm7,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm6		\n\t"\
		"addpd	%%xmm5,%%xmm7		\n\t"\
		"movaps	%%xmm4,    (%%rsi)	\n\t"\
		"movaps	%%xmm5,0x10(%%rsi)	\n\t"\
		"movaps	%%xmm0,%%xmm4		\n\t"\
		"movaps	%%xmm1,%%xmm5		\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"mulpd	0x950(%%rdi),%%xmm0	\n\t"\
		"mulpd	0x950(%%rdi),%%xmm1	\n\t"\
		"mulpd	0x960(%%rdi),%%xmm2	\n\t"\
		"mulpd	0x960(%%rdi),%%xmm3	\n\t"\
		"mulpd	0x970(%%rdi),%%xmm4	\n\t"\
		"mulpd	0x970(%%rdi),%%xmm5	\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"movaps	    (%%rsi),%%xmm4	\n\t"\
		"movaps	0x10(%%rsi),%%xmm5	\n\t"\
		"subpd	%%xmm3,%%xmm6		\n\t"\
		"subpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"movaps	%%xmm6,    (%%rax)	\n\t"\
		"movaps	%%xmm7,0x10(%%rdx)	\n\t"\
		"addpd	%%xmm6,%%xmm3		\n\t"\
		"addpd	%%xmm7,%%xmm2		\n\t"\
		"movaps	%%xmm3,    (%%rdx)	\n\t"\
		"movaps	%%xmm2,0x10(%%rax)	\n\t"\
		"subpd	%%xmm1,%%xmm4		\n\t"\
		"subpd	%%xmm0,%%xmm5		\n\t"\
		"addpd	%%xmm1,%%xmm1		\n\t"\
		"addpd	%%xmm0,%%xmm0		\n\t"\
		"movaps	%%xmm4,    (%%r15)	\n\t"\
		"movaps	%%xmm5,0x10(%%rcx)	\n\t"\
		"addpd	%%xmm4,%%xmm1		\n\t"\
		"addpd	%%xmm5,%%xmm0		\n\t"\
		"movaps	%%xmm1,    (%%rcx)	\n\t"\
		"movaps	%%xmm0,0x10(%%r15)	\n\t"\
		"							\n\t"\
	"/**************************************************************************************************/\n\t"\
	"/* For the radix-8 DIF DFTs, the input offsets always have the same pattern; outputs are permuted */\n\t"\
	"/**************************************************************************************************/\n\t"\
		"\n\t"\
		"/* SSE2_RADIX8_DIF_0TWIDDLE(r00 + 0x[02468ace]0, add0 + p[01326745]) */\n\t"\
		"movq	%[__r00],%%rax	/* i0 = r00 */\n\t"\
	"movq	%%rax,%%rsi	/* isrt2 */\n\t"\
	"addq	$0xa00,%%rsi\n\t"\
		"movq	$0x40  ,%%r15	/* i2 */	\n\t"\
		"movq	$0x80  ,%%rcx	/* i4 */	\n\t"\
		"movq	$0xc0  ,%%rdx	/* i6 */	\n\t"\
		"addq	%%rax,%%r15					\n\t"\
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
		"movaps	    (%%r15),%%xmm4			\n\t"\
		"movaps	0x10(%%r15),%%xmm5			\n\t"\
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
		"movaps	%%xmm2,     (%%r15)			\n\t"\
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
		"movaps	%%xmm6,0x010(%%r15)			\n\t"\
		"\n\t"\
		"/* Inline of SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO_B(r8,c1): Do 2nd of 2 radix-4 subtransforms, but now keep outputs in registers: */\n\t"\
		"\n\t"\
		"movq	$0x60,%%r15		/* i3 */	\n\t"\
		"movq	$0xa0,%%rcx		/* i5 */	\n\t"\
		"movq	$0xe0,%%rdx		/* i7 */	\n\t"\
		"addq	%%rax,%%r15					\n\t"\
		"addq	%%rax,%%rcx					\n\t"\
		"addq	%%rax,%%rdx					\n\t"\
		"addq	$0x20,%%rax		/* i1 */	\n\t"\
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
		"movaps	    (%%r15),%%xmm4			\n\t"\
		"movaps	0x10(%%r15),%%xmm5			\n\t"\
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
		"movaps	(%%rsi),%%xmm6	/* isrt2 */	\n\t"\
		"mulpd	%%xmm6,%%xmm2				\n\t"\
		"mulpd	%%xmm6,%%xmm5				\n\t"\
		"mulpd	%%xmm6,%%xmm4				\n\t"\
		"mulpd	%%xmm6,%%xmm3				\n\t"\
		"\n\t"\
		"/* Inline the contents of SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_A(r0): Combine the 2 radix-4 subtransforms and write outputs to temp array: */\n\t"\
		"/***** t0,1,2,3,4,5,6,7 in __i0,2,4,6 (eax,ebx,ecx,edx) *****/\n\t"\
		"/***** t8,9,a,b,c,d,e,f in __i1,xmm2,4,0,1,5,3; xmm6,7 FREE */\n\t"\
		"subq	$0x20,%%rax		/* i0 */	\n\t"\
		"movq	%%rax , %%rdi				\n\t"\
		"\n\t"\
		"movslq	%[__p06],%%rax	\n\t"\
		"movslq	%[__p07],%%r15	\n\t"\
		"movslq	%[__p04],%%rcx	\n\t"\
		"movslq	%[__p05],%%rdx	\n\t"\
		"shlq	$3,%%rax		/* Pointer offset for floating doubles */\n\t"\
		"shlq	$3,%%r15		\n\t"\
		"shlq	$3,%%rcx		\n\t"\
		"shlq	$3,%%rdx		\n\t"\
		"addq	%[__add],%%rax	/* o4 */\n\t"\
		"addq	%[__add],%%r15	/* o5 */\n\t"\
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
		"movaps	%%xmm6,    (%%r15)	/* o5r */\n\t"\
		"movaps	%%xmm7,0x10(%%r15)	/* o5i */\n\t"\
		"movaps	%%xmm2,    (%%rax)	/* o4r */\n\t"\
		"movaps	%%xmm4,0x10(%%rax)	/* o4i */\n\t"\
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
		"movaps	%%xmm6,    (%%rcx)	/* o6r */\n\t"\
		"movaps	%%xmm7,0x10(%%rdx)	/* o7i */\n\t"\
		"movaps	%%xmm3,    (%%rdx)	/* o7r */\n\t"\
		"movaps	%%xmm5,0x10(%%rcx)	/* o6i */\n\t"\
		"\n\t"\
		"subq	%[__add],%%rcx	/* p4 */\n\t"\
		"subq	%%rcx,%%rax		/* o3 = add+p2 */\n\t"\
		"subq	%%rcx,%%r15		/* o2 = add+p3 */\n\t"\
		"subq	%%rcx,%%rdx		/* o1 = add+p1 */\n\t"\
		"movq	%[__add],%%rcx	/* o0 = add+p0 */\n\t"\
		"\n\t"\
		"movaps	    (%%rdi),%%xmm6			\n\t"\
		"movaps	0x80(%%rdi),%%xmm4		\n\t"\
		"movaps	0x10(%%rdi),%%xmm7			\n\t"\
		"movaps	0x90(%%rdi),%%xmm5		\n\t"\
		"movaps	0x20(%%rdi),%%xmm2		\n\t"\
		"movaps	0x30(%%rdi),%%xmm3		\n\t"\
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
		"movaps	%%xmm6,    (%%rdx)	/* o1r */\n\t"\
		"movaps	%%xmm4,    (%%r15)	/* o2r */\n\t"\
		"movaps	%%xmm7,0x10(%%rdx)	/* o1i */\n\t"\
		"movaps	%%xmm5,0x10(%%rax)	/* o3i */\n\t"\
		"movaps	%%xmm2,    (%%rcx)	/* o0r */\n\t"\
		"movaps	%%xmm1,    (%%rax)	/* o3r */\n\t"\
		"movaps	%%xmm3,0x10(%%rcx)	/* o0r */\n\t"\
		"movaps	%%xmm0,0x10(%%r15)	/* o2i */\n\t"\
		"\n\t"\
	"/* SSE2_RADIX8_DIF_0TWIDDLE(r16 + 0x[02468ace]0, add32 + p[10237654]) */\n\t"\
		"movq	%[__r00],%%rax	\n\t"\
		"movslq	%[__p32],%%r10	\n\t"\
		"addq	$0x100,%%rax	/* i0 = r16 */\n\t"\
		"shlq	$3,%%r10		/* Pointer offset for floating doubles */\n\t"\
		"addq	%[__add],%%r10	/* add + p32 */\n\t"\
		"movq	$0x40  ,%%r15	/* i2 */	\n\t"\
		"movq	$0x80  ,%%rcx	/* i4 */	\n\t"\
		"movq	$0xc0  ,%%rdx	/* i6 */	\n\t"\
		"addq	%%rax,%%r15					\n\t"\
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
		"movaps	    (%%r15),%%xmm4			\n\t"\
		"movaps	0x10(%%r15),%%xmm5			\n\t"\
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
		"movaps	%%xmm2,     (%%r15)			\n\t"\
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
		"movaps	%%xmm6,0x010(%%r15)			\n\t"\
		"\n\t"\
		"/* Inline of SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO_B(r8,c1): Do 2nd of 2 radix-4 subtransforms, but now keep outputs in registers: */\n\t"\
		"\n\t"\
		"movq	$0x60,%%r15		/* i3 */	\n\t"\
		"movq	$0xa0,%%rcx		/* i5 */	\n\t"\
		"movq	$0xe0,%%rdx		/* i7 */	\n\t"\
		"addq	%%rax,%%r15					\n\t"\
		"addq	%%rax,%%rcx					\n\t"\
		"addq	%%rax,%%rdx					\n\t"\
		"addq	$0x20,%%rax		/* i1 */	\n\t"\
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
		"movaps	    (%%r15),%%xmm4			\n\t"\
		"movaps	0x10(%%r15),%%xmm5			\n\t"\
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
		"movaps	(%%rsi),%%xmm6	/* isrt2 */	\n\t"\
		"mulpd	%%xmm6,%%xmm2				\n\t"\
		"mulpd	%%xmm6,%%xmm5				\n\t"\
		"mulpd	%%xmm6,%%xmm4				\n\t"\
		"mulpd	%%xmm6,%%xmm3				\n\t"\
		"\n\t"\
		"/* Inline the contents of SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_A(r0): Combine the 2 radix-4 subtransforms and write outputs to temp array: */\n\t"\
		"/***** t0,1,2,3,4,5,6,7 in __i0,2,4,6 (eax,ebx,ecx,edx) *****/\n\t"\
		"/***** t8,9,a,b,c,d,e,f in __i1,xmm2,4,0,1,5,3; xmm6,7 FREE */\n\t"\
		"subq	$0x20,%%rax		/* i0 */	\n\t"\
		"movq	%%rax , %%rdi				\n\t"\
		"\n\t"\
		"movslq	%[__p06],%%rax	/* Keep same index pattern as pvs radix-8 block, but swap registers a<->b, c<->d in outputs */\n\t"\
		"movslq	%[__p07],%%r15	\n\t"\
		"movslq	%[__p04],%%rcx	\n\t"\
		"movslq	%[__p05],%%rdx	\n\t"\
		"shlq	$3,%%rax		/* Pointer offset for floating doubles */\n\t"\
		"shlq	$3,%%r15		\n\t"\
		"shlq	$3,%%rcx		\n\t"\
		"shlq	$3,%%rdx		\n\t"\
		"addq	%%r10,%%rax	/* o5 */\n\t"\
		"addq	%%r10,%%r15	/* o4 */\n\t"\
		"addq	%%r10,%%rcx	/* o7 */\n\t"\
		"addq	%%r10,%%rdx	/* o6 */\n\t"\
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
		"movaps	%%xmm6,    (%%rax)	/* o5r */\n\t"\
		"movaps	%%xmm7,0x10(%%rax)	/* o5i */\n\t"\
		"movaps	%%xmm2,    (%%r15)	/* o4r */\n\t"\
		"movaps	%%xmm4,0x10(%%r15)	/* o4i */\n\t"\
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
		"movaps	%%xmm6,    (%%rdx)	/* o6r */\n\t"\
		"movaps	%%xmm7,0x10(%%rcx)	/* o7i */\n\t"\
		"movaps	%%xmm3,    (%%rcx)	/* o7r */\n\t"\
		"movaps	%%xmm5,0x10(%%rdx)	/* o6i */\n\t"\
		"\n\t"\
		"subq	%%r10,%%rcx	/* p4 */\n\t"\
		"subq	%%rcx,%%rax		/* o2 = add+p2 */\n\t"\
		"subq	%%rcx,%%r15		/* o3 = add+p3 */\n\t"\
		"subq	%%rcx,%%rdx		/* o0 = add+p1 */\n\t"\
		"movq	%%r10,%%rcx	/* o1 = add+p0 */\n\t"\
		"\n\t"\
		"movaps	    (%%rdi),%%xmm6			\n\t"\
		"movaps	0x80(%%rdi),%%xmm4		\n\t"\
		"movaps	0x10(%%rdi),%%xmm7			\n\t"\
		"movaps	0x90(%%rdi),%%xmm5		\n\t"\
		"movaps	0x20(%%rdi),%%xmm2		\n\t"\
		"movaps	0x30(%%rdi),%%xmm3		\n\t"\
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
		"movaps	%%xmm6,    (%%rcx)	/* o1r */\n\t"\
		"movaps	%%xmm4,    (%%rax)	/* o2r */\n\t"\
		"movaps	%%xmm7,0x10(%%rcx)	/* o1i */\n\t"\
		"movaps	%%xmm5,0x10(%%r15)	/* o3i */\n\t"\
		"movaps	%%xmm2,    (%%rdx)	/* o0r */\n\t"\
		"movaps	%%xmm1,    (%%r15)	/* o3r */\n\t"\
		"movaps	%%xmm3,0x10(%%rdx)	/* o0r */\n\t"\
		"movaps	%%xmm0,0x10(%%rax)	/* o2i */\n\t"\
		"\n\t"\
	"/* SSE2_RADIX8_DIF_0TWIDDLE(r32 + 0x[02468ace]0, add24 + p[45761023]) */\n\t"\
		"movq	%[__r00],%%rax	\n\t"\
		"movslq	%[__p24],%%r10	\n\t"\
		"addq	$0x200,%%rax	/* i0 = r32 */\n\t"\
		"shlq	$3,%%r10		/* Pointer offset for floating doubles */\n\t"\
		"addq	%[__add],%%r10	/* add + p24 */\n\t"\
		"movq	$0x40  ,%%r15	/* i2 */	\n\t"\
		"movq	$0x80  ,%%rcx	/* i4 */	\n\t"\
		"movq	$0xc0  ,%%rdx	/* i6 */	\n\t"\
		"addq	%%rax,%%r15					\n\t"\
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
		"movaps	    (%%r15),%%xmm4			\n\t"\
		"movaps	0x10(%%r15),%%xmm5			\n\t"\
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
		"movaps	%%xmm2,     (%%r15)			\n\t"\
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
		"movaps	%%xmm6,0x010(%%r15)			\n\t"\
		"\n\t"\
		"/* Inline of SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO_B(r8,c1): Do 2nd of 2 radix-4 subtransforms, but now keep outputs in registers: */\n\t"\
		"\n\t"\
		"movq	$0x60,%%r15		/* i3 */	\n\t"\
		"movq	$0xa0,%%rcx		/* i5 */	\n\t"\
		"movq	$0xe0,%%rdx		/* i7 */	\n\t"\
		"addq	%%rax,%%r15					\n\t"\
		"addq	%%rax,%%rcx					\n\t"\
		"addq	%%rax,%%rdx					\n\t"\
		"addq	$0x20,%%rax		/* i1 */	\n\t"\
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
		"movaps	    (%%r15),%%xmm4			\n\t"\
		"movaps	0x10(%%r15),%%xmm5			\n\t"\
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
		"movaps	(%%rsi),%%xmm6	/* isrt2 */	\n\t"\
		"mulpd	%%xmm6,%%xmm2				\n\t"\
		"mulpd	%%xmm6,%%xmm5				\n\t"\
		"mulpd	%%xmm6,%%xmm4				\n\t"\
		"mulpd	%%xmm6,%%xmm3				\n\t"\
		"\n\t"\
		"/* Inline the contents of SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_A(r0): Combine the 2 radix-4 subtransforms and write outputs to temp array: */\n\t"\
		"/***** t0,1,2,3,4,5,6,7 in __i0,2,4,6 (eax,ebx,ecx,edx) *****/\n\t"\
		"/***** t8,9,a,b,c,d,e,f in __i1,xmm2,4,0,1,5,3; xmm6,7 FREE */\n\t"\
		"subq	$0x20,%%rax		/* i0 */	\n\t"\
		"movq	%%rax , %%rdi				\n\t"\
		"\n\t"\
		"movq	$0x00000,%%rax	\n\t"\
		"movslq	%[__p01],%%r15	\n\t"\
		"movslq	%[__p02],%%rcx	\n\t"\
		"movslq	%[__p03],%%rdx	\n\t"\
		"shlq	$3,%%rax		/* Pointer offset for floating doubles */\n\t"\
		"shlq	$3,%%r15		\n\t"\
		"shlq	$3,%%rcx		\n\t"\
		"shlq	$3,%%rdx		\n\t"\
		"addq	%%r10,%%rax	/* o5 */\n\t"\
		"addq	%%r10,%%r15	/* o4 */\n\t"\
		"addq	%%r10,%%rcx	/* o6 */\n\t"\
		"addq	%%r10,%%rdx	/* o7 */\n\t"\
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
		"movaps	%%xmm6,    (%%rax)	/* o5r */\n\t"\
		"movaps	%%xmm7,0x10(%%rax)	/* o5i */\n\t"\
		"movaps	%%xmm2,    (%%r15)	/* o4r */\n\t"\
		"movaps	%%xmm4,0x10(%%r15)	/* o4i */\n\t"\
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
		"movaps	%%xmm6,    (%%rcx)	/* o6r */\n\t"\
		"movaps	%%xmm7,0x10(%%rdx)	/* o7i */\n\t"\
		"movaps	%%xmm3,    (%%rdx)	/* o7r */\n\t"\
		"movaps	%%xmm5,0x10(%%rcx)	/* o6i */\n\t"\
		"\n\t"\
		"movslq	%[__p04],%%rax	/* p4 */\n\t"\
		"shlq	$3,%%rax		/* Pointer offset for floating doubles */\n\t"\
		"addq	%%rax,%%r15		/* o1 = add+p5 */\n\t"\
		"addq	%%rax,%%rcx		/* o3 = add+p6 */\n\t"\
		"addq	%%rax,%%rdx		/* o2 = add+p7 */\n\t"\
		"addq	%%r10,%%rax	/* o0 = add+p4 */\n\t"\
		"\n\t"\
		"movaps	    (%%rdi),%%xmm6			\n\t"\
		"movaps	0x80(%%rdi),%%xmm4		\n\t"\
		"movaps	0x10(%%rdi),%%xmm7			\n\t"\
		"movaps	0x90(%%rdi),%%xmm5		\n\t"\
		"movaps	0x20(%%rdi),%%xmm2		\n\t"\
		"movaps	0x30(%%rdi),%%xmm3		\n\t"\
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
		"movaps	%%xmm6,    (%%r15)	/* o1r */\n\t"\
		"movaps	%%xmm4,    (%%rdx)	/* o2r */\n\t"\
		"movaps	%%xmm7,0x10(%%r15)	/* o1i */\n\t"\
		"movaps	%%xmm5,0x10(%%rcx)	/* o3i */\n\t"\
		"movaps	%%xmm2,    (%%rax)	/* o0r */\n\t"\
		"movaps	%%xmm1,    (%%rcx)	/* o3r */\n\t"\
		"movaps	%%xmm3,0x10(%%rax)	/* o0r */\n\t"\
		"movaps	%%xmm0,0x10(%%rdx)	/* o2i */\n\t"\
		"\n\t"\
	"/* SSE2_RADIX8_DIF_0TWIDDLE(r48 + 0x[02468ace]0, add16 + p[32104576]) */\n\t"\
		"movq	%[__r00],%%rax	\n\t"\
		"movslq	%[__p16],%%r10	\n\t"\
		"addq	$0x300,%%rax	/* i0 = r48 */\n\t"\
		"shlq	$3,%%r10		/* Pointer offset for floating doubles */\n\t"\
		"addq	%[__add],%%r10	/* add + p16 */\n\t"\
		"movq	$0x40  ,%%r15	/* i2 */	\n\t"\
		"movq	$0x80  ,%%rcx	/* i4 */	\n\t"\
		"movq	$0xc0  ,%%rdx	/* i6 */	\n\t"\
		"addq	%%rax,%%r15					\n\t"\
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
		"movaps	    (%%r15),%%xmm4			\n\t"\
		"movaps	0x10(%%r15),%%xmm5			\n\t"\
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
		"movaps	%%xmm2,     (%%r15)			\n\t"\
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
		"movaps	%%xmm6,0x010(%%r15)			\n\t"\
		"\n\t"\
		"/* Inline of SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO_B(r8,c1): Do 2nd of 2 radix-4 subtransforms, but now keep outputs in registers: */\n\t"\
		"\n\t"\
		"movq	$0x60,%%r15		/* i3 */	\n\t"\
		"movq	$0xa0,%%rcx		/* i5 */	\n\t"\
		"movq	$0xe0,%%rdx		/* i7 */	\n\t"\
		"addq	%%rax,%%r15					\n\t"\
		"addq	%%rax,%%rcx					\n\t"\
		"addq	%%rax,%%rdx					\n\t"\
		"addq	$0x20,%%rax		/* i1 */	\n\t"\
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
		"movaps	    (%%r15),%%xmm4			\n\t"\
		"movaps	0x10(%%r15),%%xmm5			\n\t"\
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
		"movaps	(%%rsi),%%xmm6	/* isrt2 */	\n\t"\
		"mulpd	%%xmm6,%%xmm2				\n\t"\
		"mulpd	%%xmm6,%%xmm5				\n\t"\
		"mulpd	%%xmm6,%%xmm4				\n\t"\
		"mulpd	%%xmm6,%%xmm3				\n\t"\
		"\n\t"\
		"/* Inline the contents of SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_A(r0): Combine the 2 radix-4 subtransforms and write outputs to temp array: */\n\t"\
		"/***** t0,1,2,3,4,5,6,7 in __i0,2,4,6 (eax,ebx,ecx,edx) *****/\n\t"\
		"/***** t8,9,a,b,c,d,e,f in __i1,xmm2,4,0,1,5,3; xmm6,7 FREE */\n\t"\
		"subq	$0x20,%%rax		/* i0 */	\n\t"\
		"movq	%%rax , %%rdi				\n\t"\
		"\n\t"\
	"/* SSE2_RADIX8_DIF_0TWIDDLE(r48 + 0x[02468ace]0, add16 + p[32104576]) */\n\t"\
		"movslq	%[__p04],%%rax	\n\t"\
		"movslq	%[__p05],%%r15	\n\t"\
		"movslq	%[__p06],%%rcx	\n\t"\
		"movslq	%[__p07],%%rdx	\n\t"\
		"shlq	$3,%%rax		/* Pointer offset for floating doubles */\n\t"\
		"shlq	$3,%%r15		\n\t"\
		"shlq	$3,%%rcx		\n\t"\
		"shlq	$3,%%rdx		\n\t"\
		"addq	%%r10,%%rax	/* o4 */\n\t"\
		"addq	%%r10,%%r15	/* o5 */\n\t"\
		"addq	%%r10,%%rcx	/* o7 */\n\t"\
		"addq	%%r10,%%rdx	/* o6 */\n\t"\
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
		"movaps	%%xmm6,    (%%r15)	/* o5r */\n\t"\
		"movaps	%%xmm7,0x10(%%r15)	/* o5i */\n\t"\
		"movaps	%%xmm2,    (%%rax)	/* o4r */\n\t"\
		"movaps	%%xmm4,0x10(%%rax)	/* o4i */\n\t"\
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
		"movaps	%%xmm6,    (%%rdx)	/* o6r */\n\t"\
		"movaps	%%xmm7,0x10(%%rcx)	/* o7i */\n\t"\
		"movaps	%%xmm3,    (%%rcx)	/* o7r */\n\t"\
		"movaps	%%xmm5,0x10(%%rdx)	/* o6i */\n\t"\
		"\n\t"\
		"subq	%%r10,%%rax	/* p4 */\n\t"\
		"subq	%%rax,%%r15		/* o2 = add+p1 */\n\t"\
		"subq	%%rax,%%rcx		/* o1 = add+p2 */\n\t"\
		"subq	%%rax,%%rdx		/* o0 = add+p3 */\n\t"\
		"movq	%%r10,%%rax	/* o3 = add+p0 */\n\t"\
		"\n\t"\
		"movaps	    (%%rdi),%%xmm6			\n\t"\
		"movaps	0x80(%%rdi),%%xmm4		\n\t"\
		"movaps	0x10(%%rdi),%%xmm7			\n\t"\
		"movaps	0x90(%%rdi),%%xmm5		\n\t"\
		"movaps	0x20(%%rdi),%%xmm2		\n\t"\
		"movaps	0x30(%%rdi),%%xmm3		\n\t"\
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
		"movaps	%%xmm6,    (%%rcx)	/* o1r */\n\t"\
		"movaps	%%xmm4,    (%%r15)	/* o2r */\n\t"\
		"movaps	%%xmm7,0x10(%%rcx)	/* o1i */\n\t"\
		"movaps	%%xmm5,0x10(%%rax)	/* o3i */\n\t"\
		"movaps	%%xmm2,    (%%rdx)	/* o0r */\n\t"\
		"movaps	%%xmm1,    (%%rax)	/* o3r */\n\t"\
		"movaps	%%xmm3,0x10(%%rdx)	/* o0r */\n\t"\
		"movaps	%%xmm0,0x10(%%r15)	/* o2i */\n\t"\
		"\n\t"\
	"/* SSE2_RADIX8_DIF_0TWIDDLE(r64 + 0x[02468ace]0, add08 + p[67453210]) */\n\t"\
		"movq	%[__r00],%%rax	\n\t"\
		"movslq	%[__p08],%%r10	\n\t"\
		"addq	$0x400,%%rax	/* i0 = r64 */\n\t"\
		"shlq	$3,%%r10		/* Pointer offset for floating doubles */\n\t"\
		"addq	%[__add],%%r10	/* add + p08 */\n\t"\
		"movq	$0x40  ,%%r15	/* i2 */	\n\t"\
		"movq	$0x80  ,%%rcx	/* i4 */	\n\t"\
		"movq	$0xc0  ,%%rdx	/* i6 */	\n\t"\
		"addq	%%rax,%%r15					\n\t"\
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
		"movaps	    (%%r15),%%xmm4			\n\t"\
		"movaps	0x10(%%r15),%%xmm5			\n\t"\
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
		"movaps	%%xmm2,     (%%r15)			\n\t"\
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
		"movaps	%%xmm6,0x010(%%r15)			\n\t"\
		"\n\t"\
		"/* Inline of SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO_B(r8,c1): Do 2nd of 2 radix-4 subtransforms, but now keep outputs in registers: */\n\t"\
		"\n\t"\
		"movq	$0x60,%%r15		/* i3 */	\n\t"\
		"movq	$0xa0,%%rcx		/* i5 */	\n\t"\
		"movq	$0xe0,%%rdx		/* i7 */	\n\t"\
		"addq	%%rax,%%r15					\n\t"\
		"addq	%%rax,%%rcx					\n\t"\
		"addq	%%rax,%%rdx					\n\t"\
		"addq	$0x20,%%rax		/* i1 */	\n\t"\
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
		"movaps	    (%%r15),%%xmm4			\n\t"\
		"movaps	0x10(%%r15),%%xmm5			\n\t"\
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
		"movaps	(%%rsi),%%xmm6	/* isrt2 */	\n\t"\
		"mulpd	%%xmm6,%%xmm2				\n\t"\
		"mulpd	%%xmm6,%%xmm5				\n\t"\
		"mulpd	%%xmm6,%%xmm4				\n\t"\
		"mulpd	%%xmm6,%%xmm3				\n\t"\
		"\n\t"\
		"/* Inline the contents of SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_A(r0): Combine the 2 radix-4 subtransforms and write outputs to temp array: */\n\t"\
		"/***** t0,1,2,3,4,5,6,7 in __i0,2,4,6 (eax,ebx,ecx,edx) *****/\n\t"\
		"/***** t8,9,a,b,c,d,e,f in __i1,xmm2,4,0,1,5,3; xmm6,7 FREE */\n\t"\
		"subq	$0x20,%%rax		/* i0 */	\n\t"\
		"movq	%%rax , %%rdi				\n\t"\
		"\n\t"\
		"movq	$0x00000,%%rax	\n\t"\
		"movslq	%[__p01],%%r15	\n\t"\
		"movslq	%[__p02],%%rcx	\n\t"\
		"movslq	%[__p03],%%rdx	\n\t"\
		"shlq	$3,%%rax		/* Pointer offset for floating doubles */\n\t"\
		"shlq	$3,%%r15		\n\t"\
		"shlq	$3,%%rcx		\n\t"\
		"shlq	$3,%%rdx		\n\t"\
		"addq	%%r10,%%rax	/* o7 */\n\t"\
		"addq	%%r10,%%r15	/* o6 */\n\t"\
		"addq	%%r10,%%rcx	/* o5 */\n\t"\
		"addq	%%r10,%%rdx	/* o4 */\n\t"\
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
		"movaps	%%xmm6,    (%%rcx)	/* o5r */\n\t"\
		"movaps	%%xmm7,0x10(%%rcx)	/* o5i */\n\t"\
		"movaps	%%xmm2,    (%%rdx)	/* o4r */\n\t"\
		"movaps	%%xmm4,0x10(%%rdx)	/* o4i */\n\t"\
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
		"movaps	%%xmm6,    (%%r15)	/* o6r */\n\t"\
		"movaps	%%xmm7,0x10(%%rax)	/* o7i */\n\t"\
		"movaps	%%xmm3,    (%%rax)	/* o7r */\n\t"\
		"movaps	%%xmm5,0x10(%%r15)	/* o6i */\n\t"\
		"\n\t"\
		"movslq	%[__p04],%%rax	/* p4 */\n\t"\
		"shlq	$3,%%rax		/* Pointer offset for floating doubles */\n\t"\
		"addq	%%rax,%%r15		/* o3 = add+p5 */\n\t"\
		"addq	%%rax,%%rcx		/* o0 = add+p6 */\n\t"\
		"addq	%%rax,%%rdx		/* o1 = add+p7 */\n\t"\
		"addq	%%r10,%%rax	/* o2 = add+p4 */\n\t"\
		"\n\t"\
		"movaps	    (%%rdi),%%xmm6			\n\t"\
		"movaps	0x80(%%rdi),%%xmm4		\n\t"\
		"movaps	0x10(%%rdi),%%xmm7			\n\t"\
		"movaps	0x90(%%rdi),%%xmm5		\n\t"\
		"movaps	0x20(%%rdi),%%xmm2		\n\t"\
		"movaps	0x30(%%rdi),%%xmm3		\n\t"\
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
		"movaps	%%xmm6,    (%%rdx)	/* o1r */\n\t"\
		"movaps	%%xmm4,    (%%rax)	/* o2r */\n\t"\
		"movaps	%%xmm7,0x10(%%rdx)	/* o1i */\n\t"\
		"movaps	%%xmm5,0x10(%%r15)	/* o3i */\n\t"\
		"movaps	%%xmm2,    (%%rcx)	/* o0r */\n\t"\
		"movaps	%%xmm1,    (%%r15)	/* o3r */\n\t"\
		"movaps	%%xmm3,0x10(%%rcx)	/* o0r */\n\t"\
		"movaps	%%xmm0,0x10(%%rax)	/* o2i */\n\t"\
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
		 ,[__p24] "m" (Xp24)\
		 ,[__p32] "m" (Xp32)\
		 ,[__r00] "m" (Xr00)\
		 ,[__cc1] "m" (Xcc1)\
		 ,[__a0] "m" (Xa0)\
		 ,[__a1] "m" (Xa1)\
		 ,[__a2] "m" (Xa2)\
		 ,[__a3] "m" (Xa3)\
		 ,[__a4] "m" (Xa4)\
		: "rax","r15","rcx","rdx","rdi","rsi","r10","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
	);\
	}

#endif	/* radix40_ditN_cy_dif1_gcc_h_included */

