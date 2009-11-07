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
#ifndef radix36_ditN_cy_dif1_gcc_h_included
#define radix36_ditN_cy_dif1_gcc_h_included

	#define	SSE2_RADIX36_DIT_NOTWIDDLE(Xadd,Xp01,Xp02,Xp03,Xp04,Xp08,Xp12,Xp16,Xp20,Xp24,Xp28,Xp32,Xr00,Xout,Xcc1)\
	{\
	__asm__ volatile (\
		"/*	add0,1,3,2 = &a[j1+p00]+p0,1,2,3: SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_C(eax,ebx,edx,ecx, 0x120, 0x240, r00) */\n\t"\
			"movq	%[__r00],%%rsi\n\t"\
			"movq	%[__add],%%rax\n\t"\
			"movslq	%[__p01],%%r15\n\t"\
			"movslq	%[__p02],%%rcx\n\t"\
			"movslq	%[__p03],%%rdx\n\t"\
			"shlq	$3,%%r15		/* Pointer offset for floating doubles */\n\t"\
			"shlq	$3,%%rcx\n\t"\
			"shlq	$3,%%rdx\n\t"\
			"addq	%%rax,%%r15\n\t"\
			"addq	%%rax,%%rcx\n\t"\
			"addq	%%rax,%%rdx\n\t"\
			"/* ecx <-> edx */\n\t"\
			"movaps	     (%%rax),%%xmm2	\n\t"\
			"movaps	     (%%rdx),%%xmm6	\n\t"\
			"movaps	0x010(%%rax),%%xmm3	\n\t"\
			"movaps	0x010(%%rdx),%%xmm7	\n\t"\
			"movaps	     (%%r15),%%xmm0	\n\t"\
			"movaps	     (%%rcx),%%xmm4	\n\t"\
			"movaps	0x010(%%r15),%%xmm1	\n\t"\
			"movaps	0x010(%%rcx),%%xmm5	\n\t"\
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
			"subpd	%%xmm4,%%xmm0		\n\t"\
			"subpd	%%xmm7,%%xmm2		\n\t"\
			"subpd	%%xmm5,%%xmm1		\n\t"\
			"subpd	%%xmm6,%%xmm3		\n\t"\
			"movaps	%%xmm0,0x240(%%rsi)	\n\t"\
			"movaps	%%xmm2,0x360(%%rsi)	\n\t"\
			"movaps	%%xmm1,0x250(%%rsi)	\n\t"\
			"movaps	%%xmm3,0x130(%%rsi)	\n\t"\
			"addpd	%%xmm4,%%xmm4		\n\t"\
			"addpd	%%xmm7,%%xmm7		\n\t"\
			"addpd	%%xmm5,%%xmm5		\n\t"\
			"addpd	%%xmm6,%%xmm6		\n\t"\
			"addpd	%%xmm0,%%xmm4		\n\t"\
			"addpd	%%xmm2,%%xmm7		\n\t"\
			"addpd	%%xmm1,%%xmm5		\n\t"\
			"addpd	%%xmm3,%%xmm6		\n\t"\
			"movaps	%%xmm4,     (%%rsi)	\n\t"\
			"movaps	%%xmm7,0x120(%%rsi)	\n\t"\
			"movaps	%%xmm5,0x010(%%rsi)	\n\t"\
			"movaps	%%xmm6,0x370(%%rsi)	\n\t"\
		"/*	add2,3,0,1 = &a[j1+p04]+p0,1,2,3: SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_C(ecx,edx,eax,ebx, 0x120, 0x240, r04) */\n\t"\
			"addq	$0x040,%%rsi	/* r04 */\n\t"\
			"movslq	%[__p04],%%rdi\n\t"\
			"shlq	$3,%%rdi \n\t"\
			"addq	%%rdi,%%rax	/* add04 */\n\t"\
			"addq	%%rdi,%%r15\n\t"\
			"addq	%%rdi,%%rcx\n\t"\
			"addq	%%rdi,%%rdx\n\t"\
			"/* e[ab]x <-> e[cd]x */\n\t"\
			"movaps	     (%%rcx),%%xmm2	\n\t"\
			"movaps	     (%%rax),%%xmm6	\n\t"\
			"movaps	0x010(%%rcx),%%xmm3	\n\t"\
			"movaps	0x010(%%rax),%%xmm7	\n\t"\
			"movaps	     (%%rdx),%%xmm0	\n\t"\
			"movaps	     (%%r15),%%xmm4	\n\t"\
			"movaps	0x010(%%rdx),%%xmm1	\n\t"\
			"movaps	0x010(%%r15),%%xmm5	\n\t"\
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
			"subpd	%%xmm4,%%xmm0		\n\t"\
			"subpd	%%xmm7,%%xmm2		\n\t"\
			"subpd	%%xmm5,%%xmm1		\n\t"\
			"subpd	%%xmm6,%%xmm3		\n\t"\
			"movaps	%%xmm0,0x240(%%rsi)	\n\t"\
			"movaps	%%xmm2,0x360(%%rsi)	\n\t"\
			"movaps	%%xmm1,0x250(%%rsi)	\n\t"\
			"movaps	%%xmm3,0x130(%%rsi)	\n\t"\
			"addpd	%%xmm4,%%xmm4		\n\t"\
			"addpd	%%xmm7,%%xmm7		\n\t"\
			"addpd	%%xmm5,%%xmm5		\n\t"\
			"addpd	%%xmm6,%%xmm6		\n\t"\
			"addpd	%%xmm0,%%xmm4		\n\t"\
			"addpd	%%xmm2,%%xmm7		\n\t"\
			"addpd	%%xmm1,%%xmm5		\n\t"\
			"addpd	%%xmm3,%%xmm6		\n\t"\
			"movaps	%%xmm4,     (%%rsi)	\n\t"\
			"movaps	%%xmm7,0x120(%%rsi)	\n\t"\
			"movaps	%%xmm5,0x010(%%rsi)	\n\t"\
			"movaps	%%xmm6,0x370(%%rsi)	\n\t"\
		"/*	add1,0,2,3 = &a[j1+p08]+p0,1,2,3: SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_C(ebx,eax,ecx,edx, 0x120, 0x240, r02) */\n\t"\
			"subq	$0x020,%%rsi	/* r02 */\n\t"\
			"addq	%%rdi,%%rax	/* add08 */\n\t"\
			"addq	%%rdi,%%r15\n\t"\
			"addq	%%rdi,%%rcx\n\t"\
			"addq	%%rdi,%%rdx\n\t"\
			"/* eax <-> ebx */\n\t"\
			"movaps	     (%%r15),%%xmm2	\n\t"\
			"movaps	     (%%rcx),%%xmm6	\n\t"\
			"movaps	0x010(%%r15),%%xmm3	\n\t"\
			"movaps	0x010(%%rcx),%%xmm7	\n\t"\
			"movaps	     (%%rax),%%xmm0	\n\t"\
			"movaps	     (%%rdx),%%xmm4	\n\t"\
			"movaps	0x010(%%rax),%%xmm1	\n\t"\
			"movaps	0x010(%%rdx),%%xmm5	\n\t"\
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
			"subpd	%%xmm4,%%xmm0		\n\t"\
			"subpd	%%xmm7,%%xmm2		\n\t"\
			"subpd	%%xmm5,%%xmm1		\n\t"\
			"subpd	%%xmm6,%%xmm3		\n\t"\
			"movaps	%%xmm0,0x240(%%rsi)	\n\t"\
			"movaps	%%xmm2,0x360(%%rsi)	\n\t"\
			"movaps	%%xmm1,0x250(%%rsi)	\n\t"\
			"movaps	%%xmm3,0x130(%%rsi)	\n\t"\
			"addpd	%%xmm4,%%xmm4		\n\t"\
			"addpd	%%xmm7,%%xmm7		\n\t"\
			"addpd	%%xmm5,%%xmm5		\n\t"\
			"addpd	%%xmm6,%%xmm6		\n\t"\
			"addpd	%%xmm0,%%xmm4		\n\t"\
			"addpd	%%xmm2,%%xmm7		\n\t"\
			"addpd	%%xmm1,%%xmm5		\n\t"\
			"addpd	%%xmm3,%%xmm6		\n\t"\
			"movaps	%%xmm4,     (%%rsi)	\n\t"\
			"movaps	%%xmm7,0x120(%%rsi)	\n\t"\
			"movaps	%%xmm5,0x010(%%rsi)	\n\t"\
			"movaps	%%xmm6,0x370(%%rsi)	\n\t"\
		"/*	add3,2,1,0 = &a[j1+p12]+p0,1,2,3: SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_C(edx,ecx,ebx,eax, 0x120, 0x240, r0c) */\n\t"\
			"addq	$0x0a0,%%rsi	/* r0c */\n\t"\
			"addq	%%rdi,%%rax	/* add12 */\n\t"\
			"addq	%%rdi,%%r15\n\t"\
			"addq	%%rdi,%%rcx\n\t"\
			"addq	%%rdi,%%rdx\n\t"\
			"/* eax <-> edx, ebx <-> ecx */\n\t"\
			"movaps	     (%%rdx),%%xmm2	\n\t"\
			"movaps	     (%%r15),%%xmm6	\n\t"\
			"movaps	0x010(%%rdx),%%xmm3	\n\t"\
			"movaps	0x010(%%r15),%%xmm7	\n\t"\
			"movaps	     (%%rcx),%%xmm0	\n\t"\
			"movaps	     (%%rax),%%xmm4	\n\t"\
			"movaps	0x010(%%rcx),%%xmm1	\n\t"\
			"movaps	0x010(%%rax),%%xmm5	\n\t"\
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
			"subpd	%%xmm4,%%xmm0		\n\t"\
			"subpd	%%xmm7,%%xmm2		\n\t"\
			"subpd	%%xmm5,%%xmm1		\n\t"\
			"subpd	%%xmm6,%%xmm3		\n\t"\
			"movaps	%%xmm0,0x240(%%rsi)	\n\t"\
			"movaps	%%xmm2,0x360(%%rsi)	\n\t"\
			"movaps	%%xmm1,0x250(%%rsi)	\n\t"\
			"movaps	%%xmm3,0x130(%%rsi)	\n\t"\
			"addpd	%%xmm4,%%xmm4		\n\t"\
			"addpd	%%xmm7,%%xmm7		\n\t"\
			"addpd	%%xmm5,%%xmm5		\n\t"\
			"addpd	%%xmm6,%%xmm6		\n\t"\
			"addpd	%%xmm0,%%xmm4		\n\t"\
			"addpd	%%xmm2,%%xmm7		\n\t"\
			"addpd	%%xmm1,%%xmm5		\n\t"\
			"addpd	%%xmm3,%%xmm6		\n\t"\
			"movaps	%%xmm4,     (%%rsi)	\n\t"\
			"movaps	%%xmm7,0x120(%%rsi)	\n\t"\
			"movaps	%%xmm5,0x010(%%rsi)	\n\t"\
			"movaps	%%xmm6,0x370(%%rsi)	\n\t"\
		"/*	add0,1,3,2 = &a[j1+p16]+p0,1,2,3: SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_C(eax,ebx,edx,ecx, 0x120, 0x240, r0g) */\n\t"\
			"addq	$0x040,%%rsi	/* r0g */\n\t"\
			"addq	%%rdi,%%rax	/* add16 */\n\t"\
			"addq	%%rdi,%%r15\n\t"\
			"addq	%%rdi,%%rcx\n\t"\
			"addq	%%rdi,%%rdx\n\t"\
			"/* ecx <-> edx */\n\t"\
			"movaps	     (%%rax),%%xmm2	\n\t"\
			"movaps	     (%%rdx),%%xmm6	\n\t"\
			"movaps	0x010(%%rax),%%xmm3	\n\t"\
			"movaps	0x010(%%rdx),%%xmm7	\n\t"\
			"movaps	     (%%r15),%%xmm0	\n\t"\
			"movaps	     (%%rcx),%%xmm4	\n\t"\
			"movaps	0x010(%%r15),%%xmm1	\n\t"\
			"movaps	0x010(%%rcx),%%xmm5	\n\t"\
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
			"subpd	%%xmm4,%%xmm0		\n\t"\
			"subpd	%%xmm7,%%xmm2		\n\t"\
			"subpd	%%xmm5,%%xmm1		\n\t"\
			"subpd	%%xmm6,%%xmm3		\n\t"\
			"movaps	%%xmm0,0x240(%%rsi)	\n\t"\
			"movaps	%%xmm2,0x360(%%rsi)	\n\t"\
			"movaps	%%xmm1,0x250(%%rsi)	\n\t"\
			"movaps	%%xmm3,0x130(%%rsi)	\n\t"\
			"addpd	%%xmm4,%%xmm4		\n\t"\
			"addpd	%%xmm7,%%xmm7		\n\t"\
			"addpd	%%xmm5,%%xmm5		\n\t"\
			"addpd	%%xmm6,%%xmm6		\n\t"\
			"addpd	%%xmm0,%%xmm4		\n\t"\
			"addpd	%%xmm2,%%xmm7		\n\t"\
			"addpd	%%xmm1,%%xmm5		\n\t"\
			"addpd	%%xmm3,%%xmm6		\n\t"\
			"movaps	%%xmm4,     (%%rsi)	\n\t"\
			"movaps	%%xmm7,0x120(%%rsi)	\n\t"\
			"movaps	%%xmm5,0x010(%%rsi)	\n\t"\
			"movaps	%%xmm6,0x370(%%rsi)	\n\t"\
		"/*	add2,3,0,1 = &a[j1+p20]+p0,1,2,3: SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_C(ecx,edx,eax,ebx, 0x120, 0x240, r0e) */\n\t"\
			"subq	$0x020,%%rsi	/* r0e */\n\t"\
			"addq	%%rdi,%%rax	/* add20 */\n\t"\
			"addq	%%rdi,%%r15\n\t"\
			"addq	%%rdi,%%rcx\n\t"\
			"addq	%%rdi,%%rdx\n\t"\
			"/* e[ab]x <-> e[cd]x */\n\t"\
			"movaps	     (%%rcx),%%xmm2	\n\t"\
			"movaps	     (%%rax),%%xmm6	\n\t"\
			"movaps	0x010(%%rcx),%%xmm3	\n\t"\
			"movaps	0x010(%%rax),%%xmm7	\n\t"\
			"movaps	     (%%rdx),%%xmm0	\n\t"\
			"movaps	     (%%r15),%%xmm4	\n\t"\
			"movaps	0x010(%%rdx),%%xmm1	\n\t"\
			"movaps	0x010(%%r15),%%xmm5	\n\t"\
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
			"subpd	%%xmm4,%%xmm0		\n\t"\
			"subpd	%%xmm7,%%xmm2		\n\t"\
			"subpd	%%xmm5,%%xmm1		\n\t"\
			"subpd	%%xmm6,%%xmm3		\n\t"\
			"movaps	%%xmm0,0x240(%%rsi)	\n\t"\
			"movaps	%%xmm2,0x360(%%rsi)	\n\t"\
			"movaps	%%xmm1,0x250(%%rsi)	\n\t"\
			"movaps	%%xmm3,0x130(%%rsi)	\n\t"\
			"addpd	%%xmm4,%%xmm4		\n\t"\
			"addpd	%%xmm7,%%xmm7		\n\t"\
			"addpd	%%xmm5,%%xmm5		\n\t"\
			"addpd	%%xmm6,%%xmm6		\n\t"\
			"addpd	%%xmm0,%%xmm4		\n\t"\
			"addpd	%%xmm2,%%xmm7		\n\t"\
			"addpd	%%xmm1,%%xmm5		\n\t"\
			"addpd	%%xmm3,%%xmm6		\n\t"\
			"movaps	%%xmm4,     (%%rsi)	\n\t"\
			"movaps	%%xmm7,0x120(%%rsi)	\n\t"\
			"movaps	%%xmm5,0x010(%%rsi)	\n\t"\
			"movaps	%%xmm6,0x370(%%rsi)	\n\t"\
		"/*	add1,0,2,3 = &a[j1+p24]+p0,1,2,3: SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_C(ebx,eax,ecx,edx, 0x120, 0x240, r08) */\n\t"\
			"subq	$0x060,%%rsi	/* r08 */\n\t"\
			"addq	%%rdi,%%rax	/* add24 */\n\t"\
			"addq	%%rdi,%%r15\n\t"\
			"addq	%%rdi,%%rcx\n\t"\
			"addq	%%rdi,%%rdx\n\t"\
			"/* eax <-> ebx */\n\t"\
			"movaps	     (%%r15),%%xmm2	\n\t"\
			"movaps	     (%%rcx),%%xmm6	\n\t"\
			"movaps	0x010(%%r15),%%xmm3	\n\t"\
			"movaps	0x010(%%rcx),%%xmm7	\n\t"\
			"movaps	     (%%rax),%%xmm0	\n\t"\
			"movaps	     (%%rdx),%%xmm4	\n\t"\
			"movaps	0x010(%%rax),%%xmm1	\n\t"\
			"movaps	0x010(%%rdx),%%xmm5	\n\t"\
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
			"subpd	%%xmm4,%%xmm0		\n\t"\
			"subpd	%%xmm7,%%xmm2		\n\t"\
			"subpd	%%xmm5,%%xmm1		\n\t"\
			"subpd	%%xmm6,%%xmm3		\n\t"\
			"movaps	%%xmm0,0x240(%%rsi)	\n\t"\
			"movaps	%%xmm2,0x360(%%rsi)	\n\t"\
			"movaps	%%xmm1,0x250(%%rsi)	\n\t"\
			"movaps	%%xmm3,0x130(%%rsi)	\n\t"\
			"addpd	%%xmm4,%%xmm4		\n\t"\
			"addpd	%%xmm7,%%xmm7		\n\t"\
			"addpd	%%xmm5,%%xmm5		\n\t"\
			"addpd	%%xmm6,%%xmm6		\n\t"\
			"addpd	%%xmm0,%%xmm4		\n\t"\
			"addpd	%%xmm2,%%xmm7		\n\t"\
			"addpd	%%xmm1,%%xmm5		\n\t"\
			"addpd	%%xmm3,%%xmm6		\n\t"\
			"movaps	%%xmm4,     (%%rsi)	\n\t"\
			"movaps	%%xmm7,0x120(%%rsi)	\n\t"\
			"movaps	%%xmm5,0x010(%%rsi)	\n\t"\
			"movaps	%%xmm6,0x370(%%rsi)	\n\t"\
		"/*	add3,2,1,0 = &a[j1+p28]+p0,1,2,3: SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_C(edx,ecx,ebx,eax, 0x120, 0x240, r06) */\n\t"\
			"subq	$0x020,%%rsi	/* r04 */\n\t"\
			"addq	%%rdi,%%rax	/* add28 */\n\t"\
			"addq	%%rdi,%%r15\n\t"\
			"addq	%%rdi,%%rcx\n\t"\
			"addq	%%rdi,%%rdx\n\t"\
			"/* eax <-> edx, ebx <-> ecx */\n\t"\
			"movaps	     (%%rdx),%%xmm2	\n\t"\
			"movaps	     (%%r15),%%xmm6	\n\t"\
			"movaps	0x010(%%rdx),%%xmm3	\n\t"\
			"movaps	0x010(%%r15),%%xmm7	\n\t"\
			"movaps	     (%%rcx),%%xmm0	\n\t"\
			"movaps	     (%%rax),%%xmm4	\n\t"\
			"movaps	0x010(%%rcx),%%xmm1	\n\t"\
			"movaps	0x010(%%rax),%%xmm5	\n\t"\
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
			"subpd	%%xmm4,%%xmm0		\n\t"\
			"subpd	%%xmm7,%%xmm2		\n\t"\
			"subpd	%%xmm5,%%xmm1		\n\t"\
			"subpd	%%xmm6,%%xmm3		\n\t"\
			"movaps	%%xmm0,0x240(%%rsi)	\n\t"\
			"movaps	%%xmm2,0x360(%%rsi)	\n\t"\
			"movaps	%%xmm1,0x250(%%rsi)	\n\t"\
			"movaps	%%xmm3,0x130(%%rsi)	\n\t"\
			"addpd	%%xmm4,%%xmm4		\n\t"\
			"addpd	%%xmm7,%%xmm7		\n\t"\
			"addpd	%%xmm5,%%xmm5		\n\t"\
			"addpd	%%xmm6,%%xmm6		\n\t"\
			"addpd	%%xmm0,%%xmm4		\n\t"\
			"addpd	%%xmm2,%%xmm7		\n\t"\
			"addpd	%%xmm1,%%xmm5		\n\t"\
			"addpd	%%xmm3,%%xmm6		\n\t"\
			"movaps	%%xmm4,     (%%rsi)	\n\t"\
			"movaps	%%xmm7,0x120(%%rsi)	\n\t"\
			"movaps	%%xmm5,0x010(%%rsi)	\n\t"\
			"movaps	%%xmm6,0x370(%%rsi)	\n\t"\
		"/*	add0,1,3,2 = &a[j1+p32]+p0,1,2,3: SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_C(eax,ebx,edx,ecx, 0x120, 0x240, r0a) */\n\t"\
			"addq	$0x040,%%rsi	/* r04 */\n\t"\
			"addq	%%rdi,%%rax	/* add32 */\n\t"\
			"addq	%%rdi,%%r15\n\t"\
			"addq	%%rdi,%%rcx\n\t"\
			"addq	%%rdi,%%rdx\n\t"\
			"/* ecx <-> edx */\n\t"\
			"movaps	     (%%rax),%%xmm2	\n\t"\
			"movaps	     (%%rdx),%%xmm6	\n\t"\
			"movaps	0x010(%%rax),%%xmm3	\n\t"\
			"movaps	0x010(%%rdx),%%xmm7	\n\t"\
			"movaps	     (%%r15),%%xmm0	\n\t"\
			"movaps	     (%%rcx),%%xmm4	\n\t"\
			"movaps	0x010(%%r15),%%xmm1	\n\t"\
			"movaps	0x010(%%rcx),%%xmm5	\n\t"\
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
			"subpd	%%xmm4,%%xmm0		\n\t"\
			"subpd	%%xmm7,%%xmm2		\n\t"\
			"subpd	%%xmm5,%%xmm1		\n\t"\
			"subpd	%%xmm6,%%xmm3		\n\t"\
			"movaps	%%xmm0,0x240(%%rsi)	\n\t"\
			"movaps	%%xmm2,0x360(%%rsi)	\n\t"\
			"movaps	%%xmm1,0x250(%%rsi)	\n\t"\
			"movaps	%%xmm3,0x130(%%rsi)	\n\t"\
			"addpd	%%xmm4,%%xmm4		\n\t"\
			"addpd	%%xmm7,%%xmm7		\n\t"\
			"addpd	%%xmm5,%%xmm5		\n\t"\
			"addpd	%%xmm6,%%xmm6		\n\t"\
			"addpd	%%xmm0,%%xmm4		\n\t"\
			"addpd	%%xmm2,%%xmm7		\n\t"\
			"addpd	%%xmm1,%%xmm5		\n\t"\
			"addpd	%%xmm3,%%xmm6		\n\t"\
			"movaps	%%xmm4,     (%%rsi)	\n\t"\
			"movaps	%%xmm7,0x120(%%rsi)	\n\t"\
			"movaps	%%xmm5,0x010(%%rsi)	\n\t"\
			"movaps	%%xmm6,0x370(%%rsi)	\n\t"\
			"\n\t"\
			"/*...and now do 4 radix-9 transforms...*/\n\t"\
			"\n\t"\
			"/* SSE2_RADIX_09_DIT_0TWIDDLE(r00,s1p00r,s1p32r,s1p28r,s1p24r,s1p20r,s1p16r,s1p12r,s1p08r,s1p04r): */\n\t"\
			"movq	%[__out]	,%%rsi 		/* __o0-8 = s1p[00,32,28,24,20,16,12,08,04]r; esi,edi store output addresses throughout */\n\t"\
			"movq	%%rsi		,%%rdi\n\t"\
			"movq	%[__r00]	,%%rax\n\t"\
		"movq	%[__cc1]	,%%rdx 		/* edx stores trig addresses throughout */\n\t"\
			"movq	%%rax		,%%r15\n\t"\
			"movq	%%rax		,%%rcx\n\t"\
			"addq	$0x20		,%%r15\n\t"\
			"addq	$0x40		,%%rcx\n\t"\
		"addq	$0x40		,%%rdx 		/* c3m1 */\n\t"\
			"movaps	    (%%r15)	,%%xmm2\n\t"\
			"movaps	0x10(%%r15)	,%%xmm3\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t"\
			"movaps	    (%%rcx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm7\n\t"\
			"movaps	%%xmm2		,%%xmm4\n\t"\
			"movaps	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm6		,%%xmm2\n\t"\
			"addpd	%%xmm7		,%%xmm3\n\t"\
			"subpd	%%xmm6		,%%xmm4\n\t"\
			"subpd	%%xmm7		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm0\n\t"\
			"addpd	%%xmm3		,%%xmm1\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
			"movaps	%%xmm0		,    (%%rax)\n\t"\
			"movaps	%%xmm1		,0x10(%%rax)\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t"\
			"mulpd	%%xmm7		,%%xmm4\n\t"\
			"mulpd	%%xmm7		,%%xmm5\n\t"\
			"addpd	%%xmm0		,%%xmm2\n\t"\
			"addpd	%%xmm1		,%%xmm3\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t"\
			"addpd	%%xmm5		,%%xmm2\n\t"\
			"subpd	%%xmm4		,%%xmm3\n\t"\
			"subpd	%%xmm5		,%%xmm0\n\t"\
			"addpd	%%xmm4		,%%xmm1\n\t"\
			"movaps	%%xmm2		,    (%%r15)\n\t"\
			"movaps	%%xmm3		,0x10(%%r15)\n\t"\
			"movaps	%%xmm0		,    (%%rcx)\n\t"\
			"movaps	%%xmm1		,0x10(%%rcx)\n\t"\
			"addq	$0x60		,%%rax\n\t"\
			"addq	$0x60		,%%r15\n\t"\
			"addq	$0x60		,%%rcx\n\t"\
			"movaps	    (%%r15)	,%%xmm4\n\t"\
			"movaps	0x10(%%r15)	,%%xmm5\n\t"\
			"movaps	    (%%rcx)	,%%xmm2\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm3\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t"\
			"subpd	%%xmm2		,%%xmm4\n\t"\
			"subpd	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm2\n\t"\
			"addpd	%%xmm3		,%%xmm3\n\t"\
			"addpd	%%xmm4		,%%xmm2\n\t"\
			"addpd	%%xmm5		,%%xmm3\n\t"\
			"addpd	%%xmm2		,%%xmm0\n\t"\
			"addpd	%%xmm3		,%%xmm1\n\t"\
			"movaps	%%xmm0		,    (%%rax)\n\t"\
			"movaps	%%xmm1		,0x10(%%rax)\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t"\
			"mulpd	%%xmm7		,%%xmm4\n\t"\
			"mulpd	%%xmm7		,%%xmm5\n\t"\
			"addpd	%%xmm0		,%%xmm2\n\t"\
			"addpd	%%xmm1		,%%xmm3\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t"\
			"addpd	%%xmm5		,%%xmm2\n\t"\
			"subpd	%%xmm4		,%%xmm3\n\t"\
			"subpd	%%xmm5		,%%xmm0\n\t"\
			"addpd	%%xmm4		,%%xmm1\n\t"\
			"movaps	%%xmm2		,    (%%r15)\n\t"\
			"movaps	%%xmm3		,0x10(%%r15)\n\t"\
			"movaps	%%xmm0		,    (%%rcx)\n\t"\
			"movaps	%%xmm1		,0x10(%%rcx)\n\t"\
			"addq	$0x60		,%%rax\n\t"\
			"addq	$0x60		,%%r15\n\t"\
			"addq	$0x60		,%%rcx\n\t"\
			"movaps	    (%%r15)	,%%xmm4\n\t"\
			"movaps	0x10(%%r15)	,%%xmm5\n\t"\
			"movaps	    (%%rcx)	,%%xmm2\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm3\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t"\
			"subpd	%%xmm2		,%%xmm4\n\t"\
			"subpd	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm2\n\t"\
			"addpd	%%xmm3		,%%xmm3\n\t"\
			"addpd	%%xmm4		,%%xmm2\n\t"\
			"addpd	%%xmm5		,%%xmm3\n\t"\
			"addpd	%%xmm2		,%%xmm0\n\t"\
			"addpd	%%xmm3		,%%xmm1\n\t"\
			"movaps	%%xmm0		,    (%%rax)\n\t"\
			"movaps	%%xmm1		,0x10(%%rax)\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t"\
			"mulpd	%%xmm7		,%%xmm4\n\t"\
			"mulpd	%%xmm7		,%%xmm5\n\t"\
			"addpd	%%xmm0		,%%xmm2\n\t"\
			"addpd	%%xmm1		,%%xmm3\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t"\
			"addpd	%%xmm5		,%%xmm2\n\t"\
			"subpd	%%xmm4		,%%xmm3\n\t"\
			"subpd	%%xmm5		,%%xmm0\n\t"\
			"addpd	%%xmm4		,%%xmm1\n\t"\
			"movaps	%%xmm2		,    (%%r15)\n\t"\
			"movaps	%%xmm3		,0x10(%%r15)\n\t"\
			"movaps	%%xmm0		,    (%%rcx)\n\t"\
			"movaps	%%xmm1		,0x10(%%rcx)\n\t"\
			"/******************************************************************************/\n\t"\
			"/*...and now do three more radix-3 transforms, including the twiddle factors: */\n\t"\
			"/******************************************************************************/\n\t"\
			"subq	$0xc0		,%%rax 		/* __r00 */\n\t"\
			"subq	$0x80		,%%r15 		/* ebx stores r02+0xc0 = r00+0xe0; need r06 = r00+0x60 = r02+0x40 = ebx-0x80 */\n\t"\
			"subq	$0x40		,%%rcx 		/* ecx stores r04+0xc0           ; need r0c = r00+0xc0 = r04+0x80 = ebx-0x40 */\n\t"\
			"movaps	    (%%r15)	,%%xmm4\n\t"\
			"movaps	0x10(%%r15)	,%%xmm5\n\t"\
			"movaps	    (%%rcx)	,%%xmm2\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm3\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t"\
			"subpd	%%xmm2		,%%xmm4\n\t"\
			"subpd	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm2\n\t"\
			"addpd	%%xmm3		,%%xmm3\n\t"\
			"addpd	%%xmm4		,%%xmm2\n\t"\
			"addpd	%%xmm5		,%%xmm3\n\t"\
			"addpd	%%xmm2		,%%xmm0\n\t"\
			"addpd	%%xmm3		,%%xmm1\n\t"\
			"movaps	%%xmm0		,    (%%rsi)\n\t"\
			"movaps	%%xmm1		,0x10(%%rsi)\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t"\
			"mulpd	%%xmm7		,%%xmm4\n\t"\
			"mulpd	%%xmm7		,%%xmm5\n\t"\
			"addpd	%%xmm0		,%%xmm2\n\t"\
			"addpd	%%xmm1		,%%xmm3\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t"\
			"addq	$0x300		,%%rdi 		/* __o3 = s1p24r */\n\t"\
			"addq	$0x180		,%%rsi 		/* __o6 = s1p12r */\n\t"\
			"addpd	%%xmm5		,%%xmm2\n\t"\
			"subpd	%%xmm4		,%%xmm3\n\t"\
			"subpd	%%xmm5		,%%xmm0\n\t"\
			"addpd	%%xmm4		,%%xmm1\n\t"\
			"movaps	%%xmm2		,    (%%rdi)\n\t"\
			"movaps	%%xmm3		,0x10(%%rdi)\n\t"\
			"movaps	%%xmm0		,    (%%rsi)\n\t"\
			"movaps	%%xmm1		,0x10(%%rsi)\n\t"\
			"addq	$0x20		,%%rax\n\t"\
			"addq	$0x20		,%%r15\n\t"\
			"addq	$0x20		,%%rcx\n\t"\
		"subq	$0x40		,%%rdx 		/* c1 */\n\t"\
			"movaps	    (%%r15)	,%%xmm2\n\t"\
			"movaps	0x10(%%r15)	,%%xmm3\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"addpd	%%xmm1		,%%xmm2\n\t"\
			"subpd	%%xmm0		,%%xmm3\n\t"\
		"addq	$0x20		,%%rdx 		/* c2 */\n\t"\
			"movaps	    (%%rcx)	,%%xmm4\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm5\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"addpd	%%xmm1		,%%xmm4\n\t"\
			"subpd	%%xmm0		,%%xmm5\n\t"\
		"addq	$0x20		,%%rdx 		/* c3m1 */\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t"\
			"subq	$0x080		,%%rsi 		/* __o7 = s1p08r */\n\t"\
			"subpd	%%xmm4		,%%xmm2\n\t"\
			"subpd	%%xmm5		,%%xmm3\n\t"\
			"addpd	%%xmm4		,%%xmm4\n\t"\
			"addpd	%%xmm5		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm4\n\t"\
			"addpd	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm4		,%%xmm0\n\t"\
			"addpd	%%xmm5		,%%xmm1\n\t"\
			"movaps	%%xmm0		,    (%%rsi)\n\t"\
			"movaps	%%xmm1		,0x10(%%rsi)\n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t"\
			"mulpd	%%xmm7		,%%xmm2\n\t"\
			"mulpd	%%xmm7		,%%xmm3\n\t"\
			"addpd	%%xmm0		,%%xmm4\n\t"\
			"addpd	%%xmm1		,%%xmm5\n\t"\
			"addq	$0x100		,%%rdi 		/* __o1 = s1p32r */\n\t"\
			"addq	$0x180		,%%rsi 		/* __o4 = s1p20r */\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t"\
			"addpd	%%xmm3		,%%xmm4\n\t"\
			"subpd	%%xmm2		,%%xmm5\n\t"\
			"subpd	%%xmm3		,%%xmm0\n\t"\
			"addpd	%%xmm2		,%%xmm1\n\t"\
			"movaps	%%xmm4		,    (%%rdi)\n\t"\
			"movaps	%%xmm5		,0x10(%%rdi)\n\t"\
			"movaps	%%xmm0		,    (%%rsi)\n\t"\
			"movaps	%%xmm1		,0x10(%%rsi)\n\t"\
			"addq	$0x20		,%%rax\n\t"\
			"addq	$0x20		,%%r15\n\t"\
			"addq	$0x20		,%%rcx\n\t"\
		"subq	$0x20		,%%rdx 		/* c2 */\n\t"\
			"movaps	    (%%r15)	,%%xmm2\n\t"\
			"movaps	0x10(%%r15)	,%%xmm3\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"addpd	%%xmm1		,%%xmm2\n\t"\
			"subpd	%%xmm0		,%%xmm3\n\t"\
		"addq	$0x40		,%%rdx 		/* c4 */\n\t"\
			"movaps	    (%%rcx)	,%%xmm4\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm5\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
		"subq	$0x20		,%%rdx 		/* c3m1 */\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"addpd	%%xmm1		,%%xmm4\n\t"\
			"subpd	%%xmm0		,%%xmm5\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t"\
			"subq	$0x080		,%%rsi 		/* __o5 = s1p16r */\n\t"\
			"subpd	%%xmm4		,%%xmm2\n\t"\
			"subpd	%%xmm5		,%%xmm3\n\t"\
			"addpd	%%xmm4		,%%xmm4\n\t"\
			"addpd	%%xmm5		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm4\n\t"\
			"addpd	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm4		,%%xmm0\n\t"\
			"addpd	%%xmm5		,%%xmm1\n\t"\
			"movaps	%%xmm0		,    (%%rsi)\n\t"\
			"movaps	%%xmm1		,0x10(%%rsi)\n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t"\
			"mulpd	%%xmm7		,%%xmm2\n\t"\
			"mulpd	%%xmm7		,%%xmm3\n\t"\
			"addpd	%%xmm0		,%%xmm4\n\t"\
			"addpd	%%xmm1		,%%xmm5\n\t"\
			"subq	$0x380		,%%rdi 		/* __o8 = s1p04r */\n\t"\
			"addq	$0x180		,%%rsi 		/* __o2 = s1p28r */\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t"\
			"addpd	%%xmm3		,%%xmm4\n\t"\
			"subpd	%%xmm2		,%%xmm5\n\t"\
			"subpd	%%xmm3		,%%xmm0\n\t"\
			"addpd	%%xmm2		,%%xmm1\n\t"\
			"movaps	%%xmm4		,    (%%rdi)\n\t"\
			"movaps	%%xmm5		,0x10(%%rdi)\n\t"\
			"movaps	%%xmm0		,    (%%rsi)\n\t"\
			"movaps	%%xmm1		,0x10(%%rsi)\n\t"\
			"/* SSE2_RADIX_09_DIT_0TWIDDLE(r10,s1p27r,s1p23r,s1p19r,s1p15r,s1p11r,s1p07r,s1p03r,s1p35r,s1p31r): */\n\t"\
			"movq	%[__out]	,%%rsi 		/* __o0-8 = s1p[27,23,19,15,11,07,03,35,31]r; esi,edi store output addresses throughout */\n\t"\
			"addq	$0x360		,%%rsi 		/* s1p27 */\n\t"\
			"movq	%%rsi		,%%rdi\n\t"\
			"movq	%[__r00]	,%%rax\n\t"\
			"addq	$0x120		,%%rax 		/* __r10 */\n\t"\
		"movq	%[__cc1]	,%%rdx 		/* edx stores trig addresses throughout */\n\t"\
			"movq	%%rax		,%%r15\n\t"\
			"movq	%%rax		,%%rcx\n\t"\
			"addq	$0x20		,%%r15\n\t"\
			"addq	$0x40		,%%rcx\n\t"\
		"addq	$0x40		,%%rdx 		/* c3m1 */\n\t"\
			"movaps	    (%%r15)	,%%xmm2\n\t"\
			"movaps	0x10(%%r15)	,%%xmm3\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t"\
			"movaps	    (%%rcx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm7\n\t"\
			"movaps	%%xmm2		,%%xmm4\n\t"\
			"movaps	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm6		,%%xmm2\n\t"\
			"addpd	%%xmm7		,%%xmm3\n\t"\
			"subpd	%%xmm6		,%%xmm4\n\t"\
			"subpd	%%xmm7		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm0\n\t"\
			"addpd	%%xmm3		,%%xmm1\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
			"movaps	%%xmm0		,    (%%rax)\n\t"\
			"movaps	%%xmm1		,0x10(%%rax)\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t"\
			"mulpd	%%xmm7		,%%xmm4\n\t"\
			"mulpd	%%xmm7		,%%xmm5\n\t"\
			"addpd	%%xmm0		,%%xmm2\n\t"\
			"addpd	%%xmm1		,%%xmm3\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t"\
			"addpd	%%xmm5		,%%xmm2\n\t"\
			"subpd	%%xmm4		,%%xmm3\n\t"\
			"subpd	%%xmm5		,%%xmm0\n\t"\
			"addpd	%%xmm4		,%%xmm1\n\t"\
			"movaps	%%xmm2		,    (%%r15)\n\t"\
			"movaps	%%xmm3		,0x10(%%r15)\n\t"\
			"movaps	%%xmm0		,    (%%rcx)\n\t"\
			"movaps	%%xmm1		,0x10(%%rcx)\n\t"\
			"addq	$0x60		,%%rax\n\t"\
			"addq	$0x60		,%%r15\n\t"\
			"addq	$0x60		,%%rcx\n\t"\
			"movaps	    (%%r15)	,%%xmm4\n\t"\
			"movaps	0x10(%%r15)	,%%xmm5\n\t"\
			"movaps	    (%%rcx)	,%%xmm2\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm3\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t"\
			"subpd	%%xmm2		,%%xmm4\n\t"\
			"subpd	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm2\n\t"\
			"addpd	%%xmm3		,%%xmm3\n\t"\
			"addpd	%%xmm4		,%%xmm2\n\t"\
			"addpd	%%xmm5		,%%xmm3\n\t"\
			"addpd	%%xmm2		,%%xmm0\n\t"\
			"addpd	%%xmm3		,%%xmm1\n\t"\
			"movaps	%%xmm0		,    (%%rax)\n\t"\
			"movaps	%%xmm1		,0x10(%%rax)\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t"\
			"mulpd	%%xmm7		,%%xmm4\n\t"\
			"mulpd	%%xmm7		,%%xmm5\n\t"\
			"addpd	%%xmm0		,%%xmm2\n\t"\
			"addpd	%%xmm1		,%%xmm3\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t"\
			"addpd	%%xmm5		,%%xmm2\n\t"\
			"subpd	%%xmm4		,%%xmm3\n\t"\
			"subpd	%%xmm5		,%%xmm0\n\t"\
			"addpd	%%xmm4		,%%xmm1\n\t"\
			"movaps	%%xmm2		,    (%%r15)\n\t"\
			"movaps	%%xmm3		,0x10(%%r15)\n\t"\
			"movaps	%%xmm0		,    (%%rcx)\n\t"\
			"movaps	%%xmm1		,0x10(%%rcx)\n\t"\
			"addq	$0x60		,%%rax\n\t"\
			"addq	$0x60		,%%r15\n\t"\
			"addq	$0x60		,%%rcx\n\t"\
			"movaps	    (%%r15)	,%%xmm4\n\t"\
			"movaps	0x10(%%r15)	,%%xmm5\n\t"\
			"movaps	    (%%rcx)	,%%xmm2\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm3\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t"\
			"subpd	%%xmm2		,%%xmm4\n\t"\
			"subpd	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm2\n\t"\
			"addpd	%%xmm3		,%%xmm3\n\t"\
			"addpd	%%xmm4		,%%xmm2\n\t"\
			"addpd	%%xmm5		,%%xmm3\n\t"\
			"addpd	%%xmm2		,%%xmm0\n\t"\
			"addpd	%%xmm3		,%%xmm1\n\t"\
			"movaps	%%xmm0		,    (%%rax)\n\t"\
			"movaps	%%xmm1		,0x10(%%rax)\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t"\
			"mulpd	%%xmm7		,%%xmm4\n\t"\
			"mulpd	%%xmm7		,%%xmm5\n\t"\
			"addpd	%%xmm0		,%%xmm2\n\t"\
			"addpd	%%xmm1		,%%xmm3\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t"\
			"addpd	%%xmm5		,%%xmm2\n\t"\
			"subpd	%%xmm4		,%%xmm3\n\t"\
			"subpd	%%xmm5		,%%xmm0\n\t"\
			"addpd	%%xmm4		,%%xmm1\n\t"\
			"movaps	%%xmm2		,    (%%r15)\n\t"\
			"movaps	%%xmm3		,0x10(%%r15)\n\t"\
			"movaps	%%xmm0		,    (%%rcx)\n\t"\
			"movaps	%%xmm1		,0x10(%%rcx)\n\t"\
			"/******************************************************************************/\n\t"\
			"/*...and now do three more radix-3 transforms, including the twiddle factors: */\n\t"\
			"/******************************************************************************/\n\t"\
			"subq	$0xc0		,%%rax\n\t"\
			"subq	$0x80		,%%r15\n\t"\
			"subq	$0x40		,%%rcx\n\t"\
			"movaps	    (%%r15)	,%%xmm4\n\t"\
			"movaps	0x10(%%r15)	,%%xmm5\n\t"\
			"movaps	    (%%rcx)	,%%xmm2\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm3\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t"\
			"subpd	%%xmm2		,%%xmm4\n\t"\
			"subpd	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm2\n\t"\
			"addpd	%%xmm3		,%%xmm3\n\t"\
			"addpd	%%xmm4		,%%xmm2\n\t"\
			"addpd	%%xmm5		,%%xmm3\n\t"\
			"addpd	%%xmm2		,%%xmm0\n\t"\
			"addpd	%%xmm3		,%%xmm1\n\t"\
			"movaps	%%xmm0		,    (%%rsi)\n\t"\
			"movaps	%%xmm1		,0x10(%%rsi)\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t"\
			"mulpd	%%xmm7		,%%xmm4\n\t"\
			"mulpd	%%xmm7		,%%xmm5\n\t"\
			"addpd	%%xmm0		,%%xmm2\n\t"\
			"addpd	%%xmm1		,%%xmm3\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t"\
			"subq	$0x180		,%%rdi 		/* __o3 = s1p15r */\n\t"\
			"subq	$0x300		,%%rsi 		/* __o6 = s1p03r */\n\t"\
			"addpd	%%xmm5		,%%xmm2\n\t"\
			"subpd	%%xmm4		,%%xmm3\n\t"\
			"subpd	%%xmm5		,%%xmm0\n\t"\
			"addpd	%%xmm4		,%%xmm1\n\t"\
			"movaps	%%xmm2		,    (%%rdi)\n\t"\
			"movaps	%%xmm3		,0x10(%%rdi)\n\t"\
			"movaps	%%xmm0		,    (%%rsi)\n\t"\
			"movaps	%%xmm1		,0x10(%%rsi)\n\t"\
			"addq	$0x20		,%%rax\n\t"\
			"addq	$0x20		,%%r15\n\t"\
			"addq	$0x20		,%%rcx\n\t"\
		"subq	$0x40		,%%rdx 		/* c1 */\n\t"\
			"movaps	    (%%r15)	,%%xmm2\n\t"\
			"movaps	0x10(%%r15)	,%%xmm3\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"addpd	%%xmm1		,%%xmm2\n\t"\
			"subpd	%%xmm0		,%%xmm3\n\t"\
		"addq	$0x20		,%%rdx 		/* c2 */\n\t"\
			"movaps	    (%%rcx)	,%%xmm4\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm5\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"addpd	%%xmm1		,%%xmm4\n\t"\
			"subpd	%%xmm0		,%%xmm5\n\t"\
		"addq	$0x20		,%%rdx 		/* c3m1 */\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t"\
			"addq	$0x400		,%%rsi 		/* __o7 = s1p35r */\n\t"\
			"subpd	%%xmm4		,%%xmm2\n\t"\
			"subpd	%%xmm5		,%%xmm3\n\t"\
			"addpd	%%xmm4		,%%xmm4\n\t"\
			"addpd	%%xmm5		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm4\n\t"\
			"addpd	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm4		,%%xmm0\n\t"\
			"addpd	%%xmm5		,%%xmm1\n\t"\
			"movaps	%%xmm0		,    (%%rsi)\n\t"\
			"movaps	%%xmm1		,0x10(%%rsi)\n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t"\
			"mulpd	%%xmm7		,%%xmm2\n\t"\
			"mulpd	%%xmm7		,%%xmm3\n\t"\
			"addpd	%%xmm0		,%%xmm4\n\t"\
			"addpd	%%xmm1		,%%xmm5\n\t"\
			"addq	$0x100		,%%rdi 		/* __o1 = s1p23r */\n\t"\
			"subq	$0x300		,%%rsi 		/* __o4 = s1p11r */\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t"\
			"addpd	%%xmm3		,%%xmm4\n\t"\
			"subpd	%%xmm2		,%%xmm5\n\t"\
			"subpd	%%xmm3		,%%xmm0\n\t"\
			"addpd	%%xmm2		,%%xmm1\n\t"\
			"movaps	%%xmm4		,    (%%rdi)\n\t"\
			"movaps	%%xmm5		,0x10(%%rdi)\n\t"\
			"movaps	%%xmm0		,    (%%rsi)\n\t"\
			"movaps	%%xmm1		,0x10(%%rsi)\n\t"\
			"addq	$0x20		,%%rax\n\t"\
			"addq	$0x20		,%%r15\n\t"\
			"addq	$0x20		,%%rcx\n\t"\
		"subq	$0x20		,%%rdx 		/* c2 */\n\t"\
			"movaps	    (%%r15)	,%%xmm2\n\t"\
			"movaps	0x10(%%r15)	,%%xmm3\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"addpd	%%xmm1		,%%xmm2\n\t"\
			"subpd	%%xmm0		,%%xmm3\n\t"\
		"addq	$0x40		,%%rdx 		/* c4 */\n\t"\
			"movaps	    (%%rcx)	,%%xmm4\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm5\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
		"subq	$0x20		,%%rdx 		/* c3m1 */\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"addpd	%%xmm1		,%%xmm4\n\t"\
			"subpd	%%xmm0		,%%xmm5\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t"\
			"subq	$0x080		,%%rsi 		/* __o5 = s1p07r */\n\t"\
			"subpd	%%xmm4		,%%xmm2\n\t"\
			"subpd	%%xmm5		,%%xmm3\n\t"\
			"addpd	%%xmm4		,%%xmm4\n\t"\
			"addpd	%%xmm5		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm4\n\t"\
			"addpd	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm4		,%%xmm0\n\t"\
			"addpd	%%xmm5		,%%xmm1\n\t"\
			"movaps	%%xmm0		,    (%%rsi)\n\t"\
			"movaps	%%xmm1		,0x10(%%rsi)\n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t"\
			"mulpd	%%xmm7		,%%xmm2\n\t"\
			"mulpd	%%xmm7		,%%xmm3\n\t"\
			"addpd	%%xmm0		,%%xmm4\n\t"\
			"addpd	%%xmm1		,%%xmm5\n\t"\
			"addq	$0x100		,%%rdi 		/* __o8 = s1p31r */\n\t"\
			"addq	$0x180		,%%rsi 		/* __o2 = s1p19r */\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t"\
			"addpd	%%xmm3		,%%xmm4\n\t"\
			"subpd	%%xmm2		,%%xmm5\n\t"\
			"subpd	%%xmm3		,%%xmm0\n\t"\
			"addpd	%%xmm2		,%%xmm1\n\t"\
			"movaps	%%xmm4		,    (%%rdi)\n\t"\
			"movaps	%%xmm5		,0x10(%%rdi)\n\t"\
			"movaps	%%xmm0		,    (%%rsi)\n\t"\
			"movaps	%%xmm1		,0x10(%%rsi)\n\t"\
			"/* SSE2_RADIX_09_DIT_0TWIDDLE(r20,s1p18r,s1p14r,s1p10r,s1p06r,s1p02r,s1p34r,s1p30r,s1p26r,s1p22r): */\n\t"\
			"movq	%[__out]	,%%rsi 		/* __o0-8 = s1p[18,14,10,06,02,34,30,26,22]r; esi,edi store output addresses throughout */\n\t"\
			"addq	$0x240		,%%rsi 		/* s1p18 */\n\t"\
			"movq	%%rsi		,%%rdi\n\t"\
			"movq	%[__r00]	,%%rax\n\t"\
			"addq	$0x240		,%%rax 		/* __r20 */\n\t"\
		"movq	%[__cc1]	,%%rdx 		/* edx stores trig addresses throughout */\n\t"\
			"movq	%%rax		,%%r15\n\t"\
			"movq	%%rax		,%%rcx\n\t"\
			"addq	$0x20		,%%r15\n\t"\
			"addq	$0x40		,%%rcx\n\t"\
		"addq	$0x40		,%%rdx 		/* c3m1 */\n\t"\
			"movaps	    (%%r15)	,%%xmm2\n\t"\
			"movaps	0x10(%%r15)	,%%xmm3\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t"\
			"movaps	    (%%rcx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm7\n\t"\
			"movaps	%%xmm2		,%%xmm4\n\t"\
			"movaps	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm6		,%%xmm2\n\t"\
			"addpd	%%xmm7		,%%xmm3\n\t"\
			"subpd	%%xmm6		,%%xmm4\n\t"\
			"subpd	%%xmm7		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm0\n\t"\
			"addpd	%%xmm3		,%%xmm1\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
			"movaps	%%xmm0		,    (%%rax)\n\t"\
			"movaps	%%xmm1		,0x10(%%rax)\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t"\
			"mulpd	%%xmm7		,%%xmm4\n\t"\
			"mulpd	%%xmm7		,%%xmm5\n\t"\
			"addpd	%%xmm0		,%%xmm2\n\t"\
			"addpd	%%xmm1		,%%xmm3\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t"\
			"addpd	%%xmm5		,%%xmm2\n\t"\
			"subpd	%%xmm4		,%%xmm3\n\t"\
			"subpd	%%xmm5		,%%xmm0\n\t"\
			"addpd	%%xmm4		,%%xmm1\n\t"\
			"movaps	%%xmm2		,    (%%r15)\n\t"\
			"movaps	%%xmm3		,0x10(%%r15)\n\t"\
			"movaps	%%xmm0		,    (%%rcx)\n\t"\
			"movaps	%%xmm1		,0x10(%%rcx)\n\t"\
			"addq	$0x60		,%%rax\n\t"\
			"addq	$0x60		,%%r15\n\t"\
			"addq	$0x60		,%%rcx\n\t"\
			"movaps	    (%%r15)	,%%xmm4\n\t"\
			"movaps	0x10(%%r15)	,%%xmm5\n\t"\
			"movaps	    (%%rcx)	,%%xmm2\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm3\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t"\
			"subpd	%%xmm2		,%%xmm4\n\t"\
			"subpd	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm2\n\t"\
			"addpd	%%xmm3		,%%xmm3\n\t"\
			"addpd	%%xmm4		,%%xmm2\n\t"\
			"addpd	%%xmm5		,%%xmm3\n\t"\
			"addpd	%%xmm2		,%%xmm0\n\t"\
			"addpd	%%xmm3		,%%xmm1\n\t"\
			"movaps	%%xmm0		,    (%%rax)\n\t"\
			"movaps	%%xmm1		,0x10(%%rax)\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t"\
			"mulpd	%%xmm7		,%%xmm4\n\t"\
			"mulpd	%%xmm7		,%%xmm5\n\t"\
			"addpd	%%xmm0		,%%xmm2\n\t"\
			"addpd	%%xmm1		,%%xmm3\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t"\
			"addpd	%%xmm5		,%%xmm2\n\t"\
			"subpd	%%xmm4		,%%xmm3\n\t"\
			"subpd	%%xmm5		,%%xmm0\n\t"\
			"addpd	%%xmm4		,%%xmm1\n\t"\
			"movaps	%%xmm2		,    (%%r15)\n\t"\
			"movaps	%%xmm3		,0x10(%%r15)\n\t"\
			"movaps	%%xmm0		,    (%%rcx)\n\t"\
			"movaps	%%xmm1		,0x10(%%rcx)\n\t"\
			"addq	$0x60		,%%rax\n\t"\
			"addq	$0x60		,%%r15\n\t"\
			"addq	$0x60		,%%rcx\n\t"\
			"movaps	    (%%r15)	,%%xmm4\n\t"\
			"movaps	0x10(%%r15)	,%%xmm5\n\t"\
			"movaps	    (%%rcx)	,%%xmm2\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm3\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t"\
			"subpd	%%xmm2		,%%xmm4\n\t"\
			"subpd	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm2\n\t"\
			"addpd	%%xmm3		,%%xmm3\n\t"\
			"addpd	%%xmm4		,%%xmm2\n\t"\
			"addpd	%%xmm5		,%%xmm3\n\t"\
			"addpd	%%xmm2		,%%xmm0\n\t"\
			"addpd	%%xmm3		,%%xmm1\n\t"\
			"movaps	%%xmm0		,    (%%rax)\n\t"\
			"movaps	%%xmm1		,0x10(%%rax)\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t"\
			"mulpd	%%xmm7		,%%xmm4\n\t"\
			"mulpd	%%xmm7		,%%xmm5\n\t"\
			"addpd	%%xmm0		,%%xmm2\n\t"\
			"addpd	%%xmm1		,%%xmm3\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t"\
			"addpd	%%xmm5		,%%xmm2\n\t"\
			"subpd	%%xmm4		,%%xmm3\n\t"\
			"subpd	%%xmm5		,%%xmm0\n\t"\
			"addpd	%%xmm4		,%%xmm1\n\t"\
			"movaps	%%xmm2		,    (%%r15)\n\t"\
			"movaps	%%xmm3		,0x10(%%r15)\n\t"\
			"movaps	%%xmm0		,    (%%rcx)\n\t"\
			"movaps	%%xmm1		,0x10(%%rcx)\n\t"\
			"/******************************************************************************/\n\t"\
			"/*...and now do three more radix-3 transforms, including the twiddle factors: */\n\t"\
			"/******************************************************************************/\n\t"\
			"subq	$0xc0		,%%rax\n\t"\
			"subq	$0x80		,%%r15\n\t"\
			"subq	$0x40		,%%rcx\n\t"\
			"movaps	    (%%r15)	,%%xmm4\n\t"\
			"movaps	0x10(%%r15)	,%%xmm5\n\t"\
			"movaps	    (%%rcx)	,%%xmm2\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm3\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t"\
			"subpd	%%xmm2		,%%xmm4\n\t"\
			"subpd	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm2\n\t"\
			"addpd	%%xmm3		,%%xmm3\n\t"\
			"addpd	%%xmm4		,%%xmm2\n\t"\
			"addpd	%%xmm5		,%%xmm3\n\t"\
			"addpd	%%xmm2		,%%xmm0\n\t"\
			"addpd	%%xmm3		,%%xmm1\n\t"\
			"movaps	%%xmm0		,    (%%rsi)\n\t"\
			"movaps	%%xmm1		,0x10(%%rsi)\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t"\
			"mulpd	%%xmm7		,%%xmm4\n\t"\
			"mulpd	%%xmm7		,%%xmm5\n\t"\
			"addpd	%%xmm0		,%%xmm2\n\t"\
			"addpd	%%xmm1		,%%xmm3\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t"\
			"subq	$0x180		,%%rdi 		/* __o3 = s1p06r */\n\t"\
			"addq	$0x180		,%%rsi 		/* __o6 = s1p30r */\n\t"\
			"addpd	%%xmm5		,%%xmm2\n\t"\
			"subpd	%%xmm4		,%%xmm3\n\t"\
			"subpd	%%xmm5		,%%xmm0\n\t"\
			"addpd	%%xmm4		,%%xmm1\n\t"\
			"movaps	%%xmm2		,    (%%rdi)\n\t"\
			"movaps	%%xmm3		,0x10(%%rdi)\n\t"\
			"movaps	%%xmm0		,    (%%rsi)\n\t"\
			"movaps	%%xmm1		,0x10(%%rsi)\n\t"\
			"addq	$0x20		,%%rax\n\t"\
			"addq	$0x20		,%%r15\n\t"\
			"addq	$0x20		,%%rcx\n\t"\
		"subq	$0x40		,%%rdx 		/* c1 */\n\t"\
			"movaps	    (%%r15)	,%%xmm2\n\t"\
			"movaps	0x10(%%r15)	,%%xmm3\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"addpd	%%xmm1		,%%xmm2\n\t"\
			"subpd	%%xmm0		,%%xmm3\n\t"\
		"addq	$0x20		,%%rdx 		/* c2 */\n\t"\
			"movaps	    (%%rcx)	,%%xmm4\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm5\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"addpd	%%xmm1		,%%xmm4\n\t"\
			"subpd	%%xmm0		,%%xmm5\n\t"\
		"addq	$0x20		,%%rdx 		/* c3m1 */\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t"\
			"subq	$0x080		,%%rsi 		/* __o7 = s1p26r */\n\t"\
			"subpd	%%xmm4		,%%xmm2\n\t"\
			"subpd	%%xmm5		,%%xmm3\n\t"\
			"addpd	%%xmm4		,%%xmm4\n\t"\
			"addpd	%%xmm5		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm4\n\t"\
			"addpd	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm4		,%%xmm0\n\t"\
			"addpd	%%xmm5		,%%xmm1\n\t"\
			"movaps	%%xmm0		,    (%%rsi)\n\t"\
			"movaps	%%xmm1		,0x10(%%rsi)\n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t"\
			"mulpd	%%xmm7		,%%xmm2\n\t"\
			"mulpd	%%xmm7		,%%xmm3\n\t"\
			"addpd	%%xmm0		,%%xmm4\n\t"\
			"addpd	%%xmm1		,%%xmm5\n\t"\
			"addq	$0x100		,%%rdi 		/* __o1 = s1p14r */\n\t"\
			"subq	$0x300		,%%rsi 		/* __o4 = s1p02r */\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t"\
			"addpd	%%xmm3		,%%xmm4\n\t"\
			"subpd	%%xmm2		,%%xmm5\n\t"\
			"subpd	%%xmm3		,%%xmm0\n\t"\
			"addpd	%%xmm2		,%%xmm1\n\t"\
			"movaps	%%xmm4		,    (%%rdi)\n\t"\
			"movaps	%%xmm5		,0x10(%%rdi)\n\t"\
			"movaps	%%xmm0		,    (%%rsi)\n\t"\
			"movaps	%%xmm1		,0x10(%%rsi)\n\t"\
			"addq	$0x20		,%%rax\n\t"\
			"addq	$0x20		,%%r15\n\t"\
			"addq	$0x20		,%%rcx\n\t"\
		"subq	$0x20		,%%rdx 		/* c2 */\n\t"\
			"movaps	    (%%r15)	,%%xmm2\n\t"\
			"movaps	0x10(%%r15)	,%%xmm3\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"addpd	%%xmm1		,%%xmm2\n\t"\
			"subpd	%%xmm0		,%%xmm3\n\t"\
		"addq	$0x40		,%%rdx 		/* c4 */\n\t"\
			"movaps	    (%%rcx)	,%%xmm4\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm5\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
		"subq	$0x20		,%%rdx 		/* c3m1 */\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"addpd	%%xmm1		,%%xmm4\n\t"\
			"subpd	%%xmm0		,%%xmm5\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t"\
			"addq	$0x400		,%%rsi 		/* __o5 = s1p34r */\n\t"\
			"subpd	%%xmm4		,%%xmm2\n\t"\
			"subpd	%%xmm5		,%%xmm3\n\t"\
			"addpd	%%xmm4		,%%xmm4\n\t"\
			"addpd	%%xmm5		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm4\n\t"\
			"addpd	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm4		,%%xmm0\n\t"\
			"addpd	%%xmm5		,%%xmm1\n\t"\
			"movaps	%%xmm0		,    (%%rsi)\n\t"\
			"movaps	%%xmm1		,0x10(%%rsi)\n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t"\
			"mulpd	%%xmm7		,%%xmm2\n\t"\
			"mulpd	%%xmm7		,%%xmm3\n\t"\
			"addpd	%%xmm0		,%%xmm4\n\t"\
			"addpd	%%xmm1		,%%xmm5\n\t"\
			"addq	$0x100		,%%rdi 		/* __o8 = s1p22r */\n\t"\
			"subq	$0x300		,%%rsi 		/* __o2 = s1p10r */\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t"\
			"addpd	%%xmm3		,%%xmm4\n\t"\
			"subpd	%%xmm2		,%%xmm5\n\t"\
			"subpd	%%xmm3		,%%xmm0\n\t"\
			"addpd	%%xmm2		,%%xmm1\n\t"\
			"movaps	%%xmm4		,    (%%rdi)\n\t"\
			"movaps	%%xmm5		,0x10(%%rdi)\n\t"\
			"movaps	%%xmm0		,    (%%rsi)\n\t"\
			"movaps	%%xmm1		,0x10(%%rsi)\n\t"\
			"/* SSE2_RADIX_09_DIT_0TWIDDLE(r30,s1p09r,s1p05r,s1p01r,s1p33r,s1p29r,s1p25r,s1p21r,s1p17r,s1p13r): */\n\t"\
			"movq	%[__out]	,%%rsi 		/* __o0-8 = s1p[09,05,01,33,29,25,21,17,13]r; esi,edi store output addresses throughout */\n\t"\
			"addq	$0x120		,%%rsi 		/* s1p09 */\n\t"\
			"movq	%%rsi		,%%rdi\n\t"\
			"movq	%[__r00]	,%%rax\n\t"\
			"addq	$0x360		,%%rax 		/* __r20 */\n\t"\
		"movq	%[__cc1]	,%%rdx 		/* edx stores trig addresses throughout */\n\t"\
			"movq	%%rax		,%%r15\n\t"\
			"movq	%%rax		,%%rcx\n\t"\
			"addq	$0x20		,%%r15\n\t"\
			"addq	$0x40		,%%rcx\n\t"\
		"addq	$0x40		,%%rdx 		/* c3m1 */\n\t"\
			"movaps	    (%%r15)	,%%xmm2\n\t"\
			"movaps	0x10(%%r15)	,%%xmm3\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t"\
			"movaps	    (%%rcx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm7\n\t"\
			"movaps	%%xmm2		,%%xmm4\n\t"\
			"movaps	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm6		,%%xmm2\n\t"\
			"addpd	%%xmm7		,%%xmm3\n\t"\
			"subpd	%%xmm6		,%%xmm4\n\t"\
			"subpd	%%xmm7		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm0\n\t"\
			"addpd	%%xmm3		,%%xmm1\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
			"movaps	%%xmm0		,    (%%rax)\n\t"\
			"movaps	%%xmm1		,0x10(%%rax)\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t"\
			"mulpd	%%xmm7		,%%xmm4\n\t"\
			"mulpd	%%xmm7		,%%xmm5\n\t"\
			"addpd	%%xmm0		,%%xmm2\n\t"\
			"addpd	%%xmm1		,%%xmm3\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t"\
			"addpd	%%xmm5		,%%xmm2\n\t"\
			"subpd	%%xmm4		,%%xmm3\n\t"\
			"subpd	%%xmm5		,%%xmm0\n\t"\
			"addpd	%%xmm4		,%%xmm1\n\t"\
			"movaps	%%xmm2		,    (%%r15)\n\t"\
			"movaps	%%xmm3		,0x10(%%r15)\n\t"\
			"movaps	%%xmm0		,    (%%rcx)\n\t"\
			"movaps	%%xmm1		,0x10(%%rcx)\n\t"\
			"addq	$0x60		,%%rax\n\t"\
			"addq	$0x60		,%%r15\n\t"\
			"addq	$0x60		,%%rcx\n\t"\
			"movaps	    (%%r15)	,%%xmm4\n\t"\
			"movaps	0x10(%%r15)	,%%xmm5\n\t"\
			"movaps	    (%%rcx)	,%%xmm2\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm3\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t"\
			"subpd	%%xmm2		,%%xmm4\n\t"\
			"subpd	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm2\n\t"\
			"addpd	%%xmm3		,%%xmm3\n\t"\
			"addpd	%%xmm4		,%%xmm2\n\t"\
			"addpd	%%xmm5		,%%xmm3\n\t"\
			"addpd	%%xmm2		,%%xmm0\n\t"\
			"addpd	%%xmm3		,%%xmm1\n\t"\
			"movaps	%%xmm0		,    (%%rax)\n\t"\
			"movaps	%%xmm1		,0x10(%%rax)\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t"\
			"mulpd	%%xmm7		,%%xmm4\n\t"\
			"mulpd	%%xmm7		,%%xmm5\n\t"\
			"addpd	%%xmm0		,%%xmm2\n\t"\
			"addpd	%%xmm1		,%%xmm3\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t"\
			"addpd	%%xmm5		,%%xmm2\n\t"\
			"subpd	%%xmm4		,%%xmm3\n\t"\
			"subpd	%%xmm5		,%%xmm0\n\t"\
			"addpd	%%xmm4		,%%xmm1\n\t"\
			"movaps	%%xmm2		,    (%%r15)\n\t"\
			"movaps	%%xmm3		,0x10(%%r15)\n\t"\
			"movaps	%%xmm0		,    (%%rcx)\n\t"\
			"movaps	%%xmm1		,0x10(%%rcx)\n\t"\
			"addq	$0x60		,%%rax\n\t"\
			"addq	$0x60		,%%r15\n\t"\
			"addq	$0x60		,%%rcx\n\t"\
			"movaps	    (%%r15)	,%%xmm4\n\t"\
			"movaps	0x10(%%r15)	,%%xmm5\n\t"\
			"movaps	    (%%rcx)	,%%xmm2\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm3\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t"\
			"subpd	%%xmm2		,%%xmm4\n\t"\
			"subpd	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm2\n\t"\
			"addpd	%%xmm3		,%%xmm3\n\t"\
			"addpd	%%xmm4		,%%xmm2\n\t"\
			"addpd	%%xmm5		,%%xmm3\n\t"\
			"addpd	%%xmm2		,%%xmm0\n\t"\
			"addpd	%%xmm3		,%%xmm1\n\t"\
			"movaps	%%xmm0		,    (%%rax)\n\t"\
			"movaps	%%xmm1		,0x10(%%rax)\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t"\
			"mulpd	%%xmm7		,%%xmm4\n\t"\
			"mulpd	%%xmm7		,%%xmm5\n\t"\
			"addpd	%%xmm0		,%%xmm2\n\t"\
			"addpd	%%xmm1		,%%xmm3\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t"\
			"addpd	%%xmm5		,%%xmm2\n\t"\
			"subpd	%%xmm4		,%%xmm3\n\t"\
			"subpd	%%xmm5		,%%xmm0\n\t"\
			"addpd	%%xmm4		,%%xmm1\n\t"\
			"movaps	%%xmm2		,    (%%r15)\n\t"\
			"movaps	%%xmm3		,0x10(%%r15)\n\t"\
			"movaps	%%xmm0		,    (%%rcx)\n\t"\
			"movaps	%%xmm1		,0x10(%%rcx)\n\t"\
			"/******************************************************************************/\n\t"\
			"/*...and now do three more radix-3 transforms, including the twiddle factors: */\n\t"\
			"/******************************************************************************/\n\t"\
			"subq	$0xc0		,%%rax\n\t"\
			"subq	$0x80		,%%r15\n\t"\
			"subq	$0x40		,%%rcx\n\t"\
			"movaps	    (%%r15)	,%%xmm4\n\t"\
			"movaps	0x10(%%r15)	,%%xmm5\n\t"\
			"movaps	    (%%rcx)	,%%xmm2\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm3\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t"\
			"subpd	%%xmm2		,%%xmm4\n\t"\
			"subpd	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm2\n\t"\
			"addpd	%%xmm3		,%%xmm3\n\t"\
			"addpd	%%xmm4		,%%xmm2\n\t"\
			"addpd	%%xmm5		,%%xmm3\n\t"\
			"addpd	%%xmm2		,%%xmm0\n\t"\
			"addpd	%%xmm3		,%%xmm1\n\t"\
			"movaps	%%xmm0		,    (%%rsi)\n\t"\
			"movaps	%%xmm1		,0x10(%%rsi)\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t"\
			"mulpd	%%xmm7		,%%xmm4\n\t"\
			"mulpd	%%xmm7		,%%xmm5\n\t"\
			"addpd	%%xmm0		,%%xmm2\n\t"\
			"addpd	%%xmm1		,%%xmm3\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t"\
			"addq	$0x300		,%%rdi 		/* __o3 = s1p33r */\n\t"\
			"addq	$0x180		,%%rsi 		/* __o6 = s1p21r */\n\t"\
			"addpd	%%xmm5		,%%xmm2\n\t"\
			"subpd	%%xmm4		,%%xmm3\n\t"\
			"subpd	%%xmm5		,%%xmm0\n\t"\
			"addpd	%%xmm4		,%%xmm1\n\t"\
			"movaps	%%xmm2		,    (%%rdi)\n\t"\
			"movaps	%%xmm3		,0x10(%%rdi)\n\t"\
			"movaps	%%xmm0		,    (%%rsi)\n\t"\
			"movaps	%%xmm1		,0x10(%%rsi)\n\t"\
			"addq	$0x20		,%%rax\n\t"\
			"addq	$0x20		,%%r15\n\t"\
			"addq	$0x20		,%%rcx\n\t"\
		"subq	$0x40		,%%rdx 		/* c1 */\n\t"\
			"movaps	    (%%r15)	,%%xmm2\n\t"\
			"movaps	0x10(%%r15)	,%%xmm3\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"addpd	%%xmm1		,%%xmm2\n\t"\
			"subpd	%%xmm0		,%%xmm3\n\t"\
		"addq	$0x20		,%%rdx 		/* c2 */\n\t"\
			"movaps	    (%%rcx)	,%%xmm4\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm5\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"addpd	%%xmm1		,%%xmm4\n\t"\
			"subpd	%%xmm0		,%%xmm5\n\t"\
		"addq	$0x20		,%%rdx 		/* c3m1 */\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t"\
			"subq	$0x080		,%%rsi 		/* __o7 = s1p17r */\n\t"\
			"subpd	%%xmm4		,%%xmm2\n\t"\
			"subpd	%%xmm5		,%%xmm3\n\t"\
			"addpd	%%xmm4		,%%xmm4\n\t"\
			"addpd	%%xmm5		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm4\n\t"\
			"addpd	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm4		,%%xmm0\n\t"\
			"addpd	%%xmm5		,%%xmm1\n\t"\
			"movaps	%%xmm0		,    (%%rsi)\n\t"\
			"movaps	%%xmm1		,0x10(%%rsi)\n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t"\
			"mulpd	%%xmm7		,%%xmm2\n\t"\
			"mulpd	%%xmm7		,%%xmm3\n\t"\
			"addpd	%%xmm0		,%%xmm4\n\t"\
			"addpd	%%xmm1		,%%xmm5\n\t"\
			"subq	$0x380		,%%rdi 		/* __o1 = s1p05r */\n\t"\
			"addq	$0x180		,%%rsi 		/* __o4 = s1p29r */\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t"\
			"addpd	%%xmm3		,%%xmm4\n\t"\
			"subpd	%%xmm2		,%%xmm5\n\t"\
			"subpd	%%xmm3		,%%xmm0\n\t"\
			"addpd	%%xmm2		,%%xmm1\n\t"\
			"movaps	%%xmm4		,    (%%rdi)\n\t"\
			"movaps	%%xmm5		,0x10(%%rdi)\n\t"\
			"movaps	%%xmm0		,    (%%rsi)\n\t"\
			"movaps	%%xmm1		,0x10(%%rsi)\n\t"\
			"addq	$0x20		,%%rax\n\t"\
			"addq	$0x20		,%%r15\n\t"\
			"addq	$0x20		,%%rcx\n\t"\
		"subq	$0x20		,%%rdx 		/* c2 */\n\t"\
			"movaps	    (%%r15)	,%%xmm2\n\t"\
			"movaps	0x10(%%r15)	,%%xmm3\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"addpd	%%xmm1		,%%xmm2\n\t"\
			"subpd	%%xmm0		,%%xmm3\n\t"\
		"addq	$0x40		,%%rdx 		/* c4 */\n\t"\
			"movaps	    (%%rcx)	,%%xmm4\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm5\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
		"subq	$0x20		,%%rdx 		/* c3m1 */\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"addpd	%%xmm1		,%%xmm4\n\t"\
			"subpd	%%xmm0		,%%xmm5\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t"\
			"subq	$0x080		,%%rsi 		/* __o5 = s1p25r */\n\t"\
			"subpd	%%xmm4		,%%xmm2\n\t"\
			"subpd	%%xmm5		,%%xmm3\n\t"\
			"addpd	%%xmm4		,%%xmm4\n\t"\
			"addpd	%%xmm5		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm4\n\t"\
			"addpd	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm4		,%%xmm0\n\t"\
			"addpd	%%xmm5		,%%xmm1\n\t"\
			"movaps	%%xmm0		,    (%%rsi)\n\t"\
			"movaps	%%xmm1		,0x10(%%rsi)\n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t"\
			"mulpd	%%xmm7		,%%xmm2\n\t"\
			"mulpd	%%xmm7		,%%xmm3\n\t"\
			"addpd	%%xmm0		,%%xmm4\n\t"\
			"addpd	%%xmm1		,%%xmm5\n\t"\
			"addq	$0x100		,%%rdi 		/* __o8 = s1p13r */\n\t"\
			"subq	$0x300		,%%rsi 		/* __o2 = s1p01r */\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t"\
			"addpd	%%xmm3		,%%xmm4\n\t"\
			"subpd	%%xmm2		,%%xmm5\n\t"\
			"subpd	%%xmm3		,%%xmm0\n\t"\
			"addpd	%%xmm2		,%%xmm1\n\t"\
			"movaps	%%xmm4		,    (%%rdi)\n\t"\
			"movaps	%%xmm5		,0x10(%%rdi)\n\t"\
			"movaps	%%xmm0		,    (%%rsi)\n\t"\
			"movaps	%%xmm1		,0x10(%%rsi)\n\t"\
			"\n\t"\
			:					/* outputs: none */\
			: [__add] "m" (Xadd)	/* All inputs from memory addresses here */\
			 ,[__p01] "m" (Xp01)\
			 ,[__p02] "m" (Xp02)\
			 ,[__p03] "m" (Xp03)\
			 ,[__p04] "m" (Xp04)\
			 ,[__p08] "m" (Xp08)\
			 ,[__p12] "m" (Xp12)\
			 ,[__p16] "m" (Xp16)\
			 ,[__p20] "m" (Xp20)\
			 ,[__p24] "m" (Xp24)\
			 ,[__p28] "m" (Xp28)\
			 ,[__p32] "m" (Xp32)\
			 ,[__r00] "m" (Xr00)\
			 ,[__out] "m" (Xout)\
			 ,[__cc1] "m" (Xcc1)\
			: "rax","r15","rcx","rdx","rdi","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
		);\
	}


	#define	SSE2_RADIX36_DIF_NOTWIDDLE(Xadd,Xp01,Xp02,Xp03,Xp04,Xp08,Xp12,Xp16,Xp20,Xp24,Xp28,Xp32,Xr00,Xin0,Xcc1)\
	{\
	__asm__ volatile (\
			"/* SSE2_RADIX_09_DIF_0TWIDDLE(r00,s1p00r,s1p32r,s1p28r,s1p24r,s1p20r,s1p16r,s1p12r,s1p08r,s1p04r) */\n\t"\
			"movq	%[__in0]	,%%rax 		/* __i0-8 = s1p[00,32,28,24,20,16,12,08,04]r; e[abc]x store output addresses throughout */\n\t"\
			"movq	%[__r00]	,%%rsi\n\t"\
		"movq	%[__cc1]	,%%rdx 		/* edx stores trig addresses throughout */\n\t"\
			"movq	%%rax		,%%r15\n\t"\
			"movq	%%rax		,%%rcx\n\t"\
			"addq	$0x300		,%%r15 		/* __i3 = s1p24r */\n\t"\
			"addq	$0x180		,%%rcx 		/* __i6 = s1p12r */\n\t"\
		"addq	$0x40		,%%rdx 		/* c3m1 */\n\t"\
			"movaps	    (%%r15)	,%%xmm2\n\t"\
			"movaps	0x10(%%r15)	,%%xmm3\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t"\
			"movaps	    (%%rcx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm7\n\t"\
			"movaps	%%xmm2		,%%xmm4\n\t"\
			"movaps	%%xmm3		,%%xmm5\n\t"\
		"movq	%%rsi		,%%rdi\n\t"\
		"addq	$0x20		,%%rdi\n\t"\
			"addpd	%%xmm6		,%%xmm2\n\t"\
			"addpd	%%xmm7		,%%xmm3\n\t"\
			"subpd	%%xmm6		,%%xmm4\n\t"\
			"subpd	%%xmm7		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm0\n\t"\
			"addpd	%%xmm3		,%%xmm1\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
			"movaps	%%xmm0		,    (%%rsi)	/* <- t00 */\n\t"\
			"movaps	%%xmm1		,0x10(%%rsi)	/* <- t01 */\n\t"\
		"addq	$0x40		,%%rsi\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t"\
			"mulpd	%%xmm7		,%%xmm4\n\t"\
			"mulpd	%%xmm7		,%%xmm5\n\t"\
			"addpd	%%xmm0		,%%xmm2\n\t"\
			"addpd	%%xmm1		,%%xmm3\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t"\
			"subpd	%%xmm5		,%%xmm2\n\t"\
			"addpd	%%xmm4		,%%xmm3\n\t"\
			"addpd	%%xmm5		,%%xmm0\n\t"\
			"subpd	%%xmm4		,%%xmm1\n\t"\
			"movaps	%%xmm2		,    (%%rdi)	/* <- t02 */\n\t"\
			"movaps	%%xmm3		,0x10(%%rdi)	/* <- t03 */\n\t"\
			"movaps	%%xmm0		,    (%%rsi)	/* <- t04 */\n\t"\
			"movaps	%%xmm1		,0x10(%%rsi)	/* <- t05 */\n\t"\
			"addq	$0x400		,%%rax 		/* __i1 = s1p32r */\n\t"\
			"subq	$0x080		,%%r15 		/* __i4 = s1p20r */\n\t"\
			"subq	$0x080		,%%rcx 		/* __i7 = s1p08r */\n\t"\
			"movaps	    (%%r15)	,%%xmm4\n\t"\
			"movaps	0x10(%%r15)	,%%xmm5\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t"\
			"movaps	    (%%rcx)	,%%xmm2\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm3\n\t"\
		"addq	$0x40		,%%rdi\n\t"\
			"subpd	%%xmm2		,%%xmm4\n\t"\
			"subpd	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm2\n\t"\
			"addpd	%%xmm3		,%%xmm3\n\t"\
			"addpd	%%xmm4		,%%xmm2\n\t"\
			"addpd	%%xmm5		,%%xmm3\n\t"\
			"addpd	%%xmm2		,%%xmm0\n\t"\
			"addpd	%%xmm3		,%%xmm1\n\t"\
			"movaps	%%xmm0		,    (%%rdi)	/* <- t06 */\n\t"\
			"movaps	%%xmm1		,0x10(%%rdi)	/* <- t07 */\n\t"\
		"addq	$0x40		,%%rsi\n\t"\
		"addq	$0x40		,%%rdi\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t"\
			"mulpd	%%xmm7		,%%xmm4\n\t"\
			"mulpd	%%xmm7		,%%xmm5\n\t"\
			"addpd	%%xmm0		,%%xmm2\n\t"\
			"addpd	%%xmm1		,%%xmm3\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t"\
			"subpd	%%xmm5		,%%xmm2\n\t"\
			"addpd	%%xmm4		,%%xmm3\n\t"\
			"addpd	%%xmm5		,%%xmm0\n\t"\
			"subpd	%%xmm4		,%%xmm1\n\t"\
			"movaps	%%xmm2		,    (%%rsi)	/* <- t08 */\n\t"\
			"movaps	%%xmm3		,0x10(%%rsi)	/* <- t09 */\n\t"\
			"movaps	%%xmm0		,    (%%rdi)	/* <- t0a */\n\t"\
			"movaps	%%xmm1		,0x10(%%rdi)	/* <- t0b */\n\t"\
		"addq	$0x40		,%%rsi\n\t"\
			"subq	$0x080		,%%rax 		/* __i2 = s1p28r */\n\t"\
			"subq	$0x080		,%%r15 		/* __i5 = s1p16r */\n\t"\
			"subq	$0x080		,%%rcx 		/* __i8 = s1p04r */\n\t"\
			"movaps	    (%%r15)	,%%xmm4\n\t"\
			"movaps	0x10(%%r15)	,%%xmm5\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t"\
			"movaps	    (%%rcx)	,%%xmm2\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm3\n\t"\
			"subpd	%%xmm2		,%%xmm4\n\t"\
			"subpd	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm2\n\t"\
			"addpd	%%xmm3		,%%xmm3\n\t"\
			"addpd	%%xmm4		,%%xmm2\n\t"\
			"addpd	%%xmm5		,%%xmm3\n\t"\
			"addpd	%%xmm2		,%%xmm0\n\t"\
			"addpd	%%xmm3		,%%xmm1\n\t"\
			"movaps	%%xmm0		,    (%%rsi)	/* <- t0c */\n\t"\
			"movaps	%%xmm1		,0x10(%%rsi)	/* <- t0d */\n\t"\
		"addq	$0x40		,%%rdi\n\t"\
		"addq	$0x40		,%%rsi\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t"\
			"mulpd	%%xmm7		,%%xmm4\n\t"\
			"mulpd	%%xmm7		,%%xmm5\n\t"\
			"addpd	%%xmm0		,%%xmm2\n\t"\
			"addpd	%%xmm1		,%%xmm3\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t"\
			"subpd	%%xmm5		,%%xmm2\n\t"\
			"addpd	%%xmm4		,%%xmm3\n\t"\
			"addpd	%%xmm5		,%%xmm0\n\t"\
			"subpd	%%xmm4		,%%xmm1\n\t"\
			"movaps	%%xmm2		,    (%%rdi)	/* <- t0e */\n\t"\
			"movaps	%%xmm3		,0x10(%%rdi)	/* <- t0f */\n\t"\
			"movaps	%%xmm0		,    (%%rsi)	/* <- t0g */\n\t"\
			"movaps	%%xmm1		,0x10(%%rsi)	/* <- t0h */\n\t"\
			"/******************************************************************************/\n\t"\
			"/*...and now do three more radix-3 transforms, including the twiddle factors: */\n\t"\
			"/******************************************************************************/\n\t"\
			"movq	%[__r00]	,%%rax 		/* __r00 */\n\t"\
			"movq	%%rax		,%%r15\n\t"\
			"movq	%%rax		,%%rcx\n\t"\
			"addq	$0x60		,%%r15 		/* __r06 */\n\t"\
			"addq	$0xc0		,%%rcx 		/* __r0c */\n\t"\
			"movaps	    (%%r15)	,%%xmm4\n\t"\
			"movaps	0x10(%%r15)	,%%xmm5\n\t"\
			"movaps	    (%%rcx)	,%%xmm2\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm3\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t"\
			"subpd	%%xmm2		,%%xmm4\n\t"\
			"subpd	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm2\n\t"\
			"addpd	%%xmm3		,%%xmm3\n\t"\
			"addpd	%%xmm4		,%%xmm2\n\t"\
			"addpd	%%xmm5		,%%xmm3\n\t"\
			"addpd	%%xmm2		,%%xmm0\n\t"\
			"addpd	%%xmm3		,%%xmm1\n\t"\
			"movaps	%%xmm0		,    (%%rax)\n\t"\
			"movaps	%%xmm1		,0x10(%%rax)\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t"\
			"mulpd	%%xmm7		,%%xmm4\n\t"\
			"mulpd	%%xmm7		,%%xmm5\n\t"\
			"addpd	%%xmm0		,%%xmm2\n\t"\
			"addpd	%%xmm1		,%%xmm3\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t"\
			"subpd	%%xmm5		,%%xmm2\n\t"\
			"addpd	%%xmm4		,%%xmm3\n\t"\
			"addpd	%%xmm5		,%%xmm0\n\t"\
			"subpd	%%xmm4		,%%xmm1\n\t"\
			"movaps	%%xmm2		,    (%%r15)\n\t"\
			"movaps	%%xmm3		,0x10(%%r15)\n\t"\
			"movaps	%%xmm0		,    (%%rcx)\n\t"\
			"movaps	%%xmm1		,0x10(%%rcx)\n\t"\
		"subq	$0x40		,%%rdx 		/* c1 */\n\t"\
			"addq	$0x20		,%%rax 		/* __r02 */\n\t"\
			"addq	$0x20		,%%r15 		/* __r08 */\n\t"\
			"addq	$0x20		,%%rcx 		/* __r0e */\n\t"\
			"movaps	    (%%r15)	,%%xmm2\n\t"\
			"movaps	0x10(%%r15)	,%%xmm3\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
		"addq	$0x20		,%%rdx 		/* c2 */\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"subpd	%%xmm1		,%%xmm2\n\t"\
			"addpd	%%xmm0		,%%xmm3\n\t"\
			"movaps	    (%%rcx)	,%%xmm4\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm5\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
		"addq	$0x20		,%%rdx 		/* c3m1 */\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"subpd	%%xmm1		,%%xmm4\n\t"\
			"addpd	%%xmm0		,%%xmm5\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t"\
			"subpd	%%xmm4		,%%xmm2\n\t"\
			"subpd	%%xmm5		,%%xmm3\n\t"\
			"addpd	%%xmm4		,%%xmm4\n\t"\
			"addpd	%%xmm5		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm4\n\t"\
			"addpd	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm4		,%%xmm0\n\t"\
			"addpd	%%xmm5		,%%xmm1\n\t"\
			"movaps	%%xmm0		,    (%%rax)\n\t"\
			"movaps	%%xmm1		,0x10(%%rax)\n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t"\
			"mulpd	%%xmm7		,%%xmm2\n\t"\
			"mulpd	%%xmm7		,%%xmm3\n\t"\
			"addpd	%%xmm0		,%%xmm4\n\t"\
			"addpd	%%xmm1		,%%xmm5\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t"\
			"subpd	%%xmm3		,%%xmm4\n\t"\
			"addpd	%%xmm2		,%%xmm5\n\t"\
			"addpd	%%xmm3		,%%xmm0\n\t"\
			"subpd	%%xmm2		,%%xmm1\n\t"\
			"movaps	%%xmm4		,    (%%r15)\n\t"\
			"movaps	%%xmm5		,0x10(%%r15)\n\t"\
			"movaps	%%xmm0		,    (%%rcx)\n\t"\
			"movaps	%%xmm1		,0x10(%%rcx)\n\t"\
		"subq	$0x20		,%%rdx 		/* c2 */\n\t"\
			"addq	$0x20		,%%rax 		/* __r04 */\n\t"\
			"addq	$0x20		,%%r15 		/* __r0a */\n\t"\
			"addq	$0x20		,%%rcx 		/* __r0g */\n\t"\
			"movaps	    (%%r15)	,%%xmm2\n\t"\
			"movaps	0x10(%%r15)	,%%xmm3\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
		"addq	$0x40		,%%rdx 		/* c4 */\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"subpd	%%xmm1		,%%xmm2\n\t"\
			"addpd	%%xmm0		,%%xmm3\n\t"\
			"movaps	    (%%rcx)	,%%xmm4\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm5\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
		"subq	$0x20		,%%rdx 		/* c3m1 */\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"subpd	%%xmm1		,%%xmm4\n\t"\
			"addpd	%%xmm0		,%%xmm5\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t"\
			"\n\t"\
			"subpd	%%xmm4		,%%xmm2\n\t"\
			"subpd	%%xmm5		,%%xmm3\n\t"\
			"addpd	%%xmm4		,%%xmm4\n\t"\
			"addpd	%%xmm5		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm4\n\t"\
			"addpd	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm4		,%%xmm0\n\t"\
			"addpd	%%xmm5		,%%xmm1\n\t"\
			"movaps	%%xmm0		,    (%%rax)\n\t"\
			"movaps	%%xmm1		,0x10(%%rax)\n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t"\
			"mulpd	%%xmm7		,%%xmm2\n\t"\
			"mulpd	%%xmm7		,%%xmm3\n\t"\
			"addpd	%%xmm0		,%%xmm4\n\t"\
			"addpd	%%xmm1		,%%xmm5\n\t"\
			"\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t"\
			"subpd	%%xmm3		,%%xmm4\n\t"\
			"addpd	%%xmm2		,%%xmm5\n\t"\
			"addpd	%%xmm3		,%%xmm0\n\t"\
			"subpd	%%xmm2		,%%xmm1\n\t"\
			"movaps	%%xmm4		,    (%%r15)\n\t"\
			"movaps	%%xmm5		,0x10(%%r15)\n\t"\
			"movaps	%%xmm0		,    (%%rcx)\n\t"\
			"movaps	%%xmm1		,0x10(%%rcx)\n\t"\
			"/* SSE2_RADIX_09_DIF_0TWIDDLE(r10,s1p27r,s1p23r,s1p19r,s1p15r,s1p11r,s1p07r,s1p03r,s1p35r,s1p31r) */\n\t"\
			"movq	%[__in0]	,%%rax 		/* __i0-8 = s1p[27,23,19,15,11,07,03,35,31]r; e[abc]x store output addresses throughout */\n\t"\
			"movq	%[__r00]	,%%rsi\n\t"\
			"addq	$0x360		,%%rax 		/* __i0 = s1p27r */\n\t"\
			"addq	$0x120		,%%rsi 		/* r10 */\n\t"\
		"movq	%[__cc1]	,%%rdx 		/* edx stores trig addresses throughout */\n\t"\
			"movq	%%rax		,%%r15\n\t"\
			"movq	%%rax		,%%rcx\n\t"\
			"subq	$0x180		,%%r15 		/* __i3 = s1p15r */\n\t"\
			"subq	$0x300		,%%rcx 		/* __i6 = s1p03r */\n\t"\
		"addq	$0x40		,%%rdx 		/* c3m1 */\n\t"\
			"movaps	    (%%r15)	,%%xmm2\n\t"\
			"movaps	0x10(%%r15)	,%%xmm3\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t"\
			"movaps	    (%%rcx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm7\n\t"\
			"movaps	%%xmm2		,%%xmm4\n\t"\
			"movaps	%%xmm3		,%%xmm5\n\t"\
		"movq	%%rsi		,%%rdi\n\t"\
		"addq	$0x20		,%%rdi\n\t"\
			"addpd	%%xmm6		,%%xmm2\n\t"\
			"addpd	%%xmm7		,%%xmm3\n\t"\
			"subpd	%%xmm6		,%%xmm4\n\t"\
			"subpd	%%xmm7		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm0\n\t"\
			"addpd	%%xmm3		,%%xmm1\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
			"movaps	%%xmm0		,    (%%rsi)	/* <- t00 */\n\t"\
			"movaps	%%xmm1		,0x10(%%rsi)	/* <- t01 */\n\t"\
		"addq	$0x40		,%%rsi\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t"\
			"mulpd	%%xmm7		,%%xmm4\n\t"\
			"mulpd	%%xmm7		,%%xmm5\n\t"\
			"addpd	%%xmm0		,%%xmm2\n\t"\
			"addpd	%%xmm1		,%%xmm3\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t"\
			"subpd	%%xmm5		,%%xmm2\n\t"\
			"addpd	%%xmm4		,%%xmm3\n\t"\
			"addpd	%%xmm5		,%%xmm0\n\t"\
			"subpd	%%xmm4		,%%xmm1\n\t"\
			"movaps	%%xmm2		,    (%%rdi)	/* <- t02 */\n\t"\
			"movaps	%%xmm3		,0x10(%%rdi)	/* <- t03 */\n\t"\
			"movaps	%%xmm0		,    (%%rsi)	/* <- t04 */\n\t"\
			"movaps	%%xmm1		,0x10(%%rsi)	/* <- t05 */\n\t"\
			"subq	$0x080		,%%rax 		/* __i1 = s1p23r */\n\t"\
			"subq	$0x080		,%%r15 		/* __i4 = s1p11r */\n\t"\
			"addq	$0x400		,%%rcx 		/* __i7 = s1p35r */\n\t"\
			"movaps	    (%%r15)	,%%xmm4\n\t"\
			"movaps	0x10(%%r15)	,%%xmm5\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t"\
			"movaps	    (%%rcx)	,%%xmm2\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm3\n\t"\
		"addq	$0x40		,%%rdi\n\t"\
			"subpd	%%xmm2		,%%xmm4\n\t"\
			"subpd	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm2\n\t"\
			"addpd	%%xmm3		,%%xmm3\n\t"\
			"addpd	%%xmm4		,%%xmm2\n\t"\
			"addpd	%%xmm5		,%%xmm3\n\t"\
			"addpd	%%xmm2		,%%xmm0\n\t"\
			"addpd	%%xmm3		,%%xmm1\n\t"\
			"movaps	%%xmm0		,    (%%rdi)	/* <- t06 */\n\t"\
			"movaps	%%xmm1		,0x10(%%rdi)	/* <- t07 */\n\t"\
		"addq	$0x40		,%%rsi\n\t"\
		"addq	$0x40		,%%rdi\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t"\
			"mulpd	%%xmm7		,%%xmm4\n\t"\
			"mulpd	%%xmm7		,%%xmm5\n\t"\
			"addpd	%%xmm0		,%%xmm2\n\t"\
			"addpd	%%xmm1		,%%xmm3\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t"\
			"subpd	%%xmm5		,%%xmm2\n\t"\
			"addpd	%%xmm4		,%%xmm3\n\t"\
			"addpd	%%xmm5		,%%xmm0\n\t"\
			"subpd	%%xmm4		,%%xmm1\n\t"\
			"movaps	%%xmm2		,    (%%rsi)	/* <- t08 */\n\t"\
			"movaps	%%xmm3		,0x10(%%rsi)	/* <- t09 */\n\t"\
			"movaps	%%xmm0		,    (%%rdi)	/* <- t0a */\n\t"\
			"movaps	%%xmm1		,0x10(%%rdi)	/* <- t0b */\n\t"\
		"addq	$0x40		,%%rsi\n\t"\
			"subq	$0x080		,%%rax 		/* __i2 = s1p19r */\n\t"\
			"subq	$0x080		,%%r15 		/* __i5 = s1p07r */\n\t"\
			"subq	$0x080		,%%rcx 		/* __i8 = s1p31r */\n\t"\
			"movaps	    (%%r15)	,%%xmm4\n\t"\
			"movaps	0x10(%%r15)	,%%xmm5\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t"\
			"movaps	    (%%rcx)	,%%xmm2\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm3\n\t"\
			"subpd	%%xmm2		,%%xmm4\n\t"\
			"subpd	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm2\n\t"\
			"addpd	%%xmm3		,%%xmm3\n\t"\
			"addpd	%%xmm4		,%%xmm2\n\t"\
			"addpd	%%xmm5		,%%xmm3\n\t"\
			"addpd	%%xmm2		,%%xmm0\n\t"\
			"addpd	%%xmm3		,%%xmm1\n\t"\
			"movaps	%%xmm0		,    (%%rsi)	/* <- t0c */\n\t"\
			"movaps	%%xmm1		,0x10(%%rsi)	/* <- t0d */\n\t"\
		"addq	$0x40		,%%rdi\n\t"\
		"addq	$0x40		,%%rsi\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t"\
			"mulpd	%%xmm7		,%%xmm4\n\t"\
			"mulpd	%%xmm7		,%%xmm5\n\t"\
			"addpd	%%xmm0		,%%xmm2\n\t"\
			"addpd	%%xmm1		,%%xmm3\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t"\
			"subpd	%%xmm5		,%%xmm2\n\t"\
			"addpd	%%xmm4		,%%xmm3\n\t"\
			"addpd	%%xmm5		,%%xmm0\n\t"\
			"subpd	%%xmm4		,%%xmm1\n\t"\
			"movaps	%%xmm2		,    (%%rdi)	/* <- t0e */\n\t"\
			"movaps	%%xmm3		,0x10(%%rdi)	/* <- t0f */\n\t"\
			"movaps	%%xmm0		,    (%%rsi)	/* <- t0g */\n\t"\
			"movaps	%%xmm1		,0x10(%%rsi)	/* <- t0h */\n\t"\
			"/******************************************************************************/\n\t"\
			"/*...and now do three more radix-3 transforms, including the twiddle factors: */\n\t"\
			"/******************************************************************************/\n\t"\
			"movq	%[__r00]	,%%rax 		/* __r00 */\n\t"\
			"addq	$0x120		,%%rax 		/* r10 */\n\t"\
			"movq	%%rax		,%%r15\n\t"\
			"movq	%%rax		,%%rcx\n\t"\
			"addq	$0x60		,%%r15 		/* __r16 */\n\t"\
			"addq	$0xc0		,%%rcx 		/* __r1c */\n\t"\
			"movaps	    (%%r15)	,%%xmm4\n\t"\
			"movaps	0x10(%%r15)	,%%xmm5\n\t"\
			"movaps	    (%%rcx)	,%%xmm2\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm3\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t"\
			"subpd	%%xmm2		,%%xmm4\n\t"\
			"subpd	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm2\n\t"\
			"addpd	%%xmm3		,%%xmm3\n\t"\
			"addpd	%%xmm4		,%%xmm2\n\t"\
			"addpd	%%xmm5		,%%xmm3\n\t"\
			"addpd	%%xmm2		,%%xmm0\n\t"\
			"addpd	%%xmm3		,%%xmm1\n\t"\
			"movaps	%%xmm0		,    (%%rax)\n\t"\
			"movaps	%%xmm1		,0x10(%%rax)\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t"\
			"mulpd	%%xmm7		,%%xmm4\n\t"\
			"mulpd	%%xmm7		,%%xmm5\n\t"\
			"addpd	%%xmm0		,%%xmm2\n\t"\
			"addpd	%%xmm1		,%%xmm3\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t"\
			"subpd	%%xmm5		,%%xmm2\n\t"\
			"addpd	%%xmm4		,%%xmm3\n\t"\
			"addpd	%%xmm5		,%%xmm0\n\t"\
			"subpd	%%xmm4		,%%xmm1\n\t"\
			"movaps	%%xmm2		,    (%%r15)\n\t"\
			"movaps	%%xmm3		,0x10(%%r15)\n\t"\
			"movaps	%%xmm0		,    (%%rcx)\n\t"\
			"movaps	%%xmm1		,0x10(%%rcx)\n\t"\
		"subq	$0x40		,%%rdx 		/* c1 */\n\t"\
			"addq	$0x20		,%%rax 		/* __r12 */\n\t"\
			"addq	$0x20		,%%r15 		/* __r18 */\n\t"\
			"addq	$0x20		,%%rcx 		/* __r1e */\n\t"\
			"movaps	    (%%r15)	,%%xmm2\n\t"\
			"movaps	0x10(%%r15)	,%%xmm3\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
		"addq	$0x20		,%%rdx 		/* c2 */\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"subpd	%%xmm1		,%%xmm2\n\t"\
			"addpd	%%xmm0		,%%xmm3\n\t"\
			"movaps	    (%%rcx)	,%%xmm4\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm5\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
		"addq	$0x20		,%%rdx 		/* c3m1 */\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"subpd	%%xmm1		,%%xmm4\n\t"\
			"addpd	%%xmm0		,%%xmm5\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t"\
			"subpd	%%xmm4		,%%xmm2\n\t"\
			"subpd	%%xmm5		,%%xmm3\n\t"\
			"addpd	%%xmm4		,%%xmm4\n\t"\
			"addpd	%%xmm5		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm4\n\t"\
			"addpd	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm4		,%%xmm0\n\t"\
			"addpd	%%xmm5		,%%xmm1\n\t"\
			"movaps	%%xmm0		,    (%%rax)\n\t"\
			"movaps	%%xmm1		,0x10(%%rax)\n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t"\
			"mulpd	%%xmm7		,%%xmm2\n\t"\
			"mulpd	%%xmm7		,%%xmm3\n\t"\
			"addpd	%%xmm0		,%%xmm4\n\t"\
			"addpd	%%xmm1		,%%xmm5\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t"\
			"subpd	%%xmm3		,%%xmm4\n\t"\
			"addpd	%%xmm2		,%%xmm5\n\t"\
			"addpd	%%xmm3		,%%xmm0\n\t"\
			"subpd	%%xmm2		,%%xmm1\n\t"\
			"movaps	%%xmm4		,    (%%r15)\n\t"\
			"movaps	%%xmm5		,0x10(%%r15)\n\t"\
			"movaps	%%xmm0		,    (%%rcx)\n\t"\
			"movaps	%%xmm1		,0x10(%%rcx)\n\t"\
		"subq	$0x20		,%%rdx 		/* c2 */\n\t"\
			"addq	$0x20		,%%rax 		/* __r14 */\n\t"\
			"addq	$0x20		,%%r15 		/* __r1a */\n\t"\
			"addq	$0x20		,%%rcx 		/* __r1g */\n\t"\
			"movaps	    (%%r15)	,%%xmm2\n\t"\
			"movaps	0x10(%%r15)	,%%xmm3\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
		"addq	$0x40		,%%rdx 		/* c4 */\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"subpd	%%xmm1		,%%xmm2\n\t"\
			"addpd	%%xmm0		,%%xmm3\n\t"\
			"movaps	    (%%rcx)	,%%xmm4\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm5\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
		"subq	$0x20		,%%rdx 		/* c3m1 */\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"subpd	%%xmm1		,%%xmm4\n\t"\
			"addpd	%%xmm0		,%%xmm5\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t"\
			"\n\t"\
			"subpd	%%xmm4		,%%xmm2\n\t"\
			"subpd	%%xmm5		,%%xmm3\n\t"\
			"addpd	%%xmm4		,%%xmm4\n\t"\
			"addpd	%%xmm5		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm4\n\t"\
			"addpd	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm4		,%%xmm0\n\t"\
			"addpd	%%xmm5		,%%xmm1\n\t"\
			"movaps	%%xmm0		,    (%%rax)\n\t"\
			"movaps	%%xmm1		,0x10(%%rax)\n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t"\
			"mulpd	%%xmm7		,%%xmm2\n\t"\
			"mulpd	%%xmm7		,%%xmm3\n\t"\
			"addpd	%%xmm0		,%%xmm4\n\t"\
			"addpd	%%xmm1		,%%xmm5\n\t"\
			"\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t"\
			"subpd	%%xmm3		,%%xmm4\n\t"\
			"addpd	%%xmm2		,%%xmm5\n\t"\
			"addpd	%%xmm3		,%%xmm0\n\t"\
			"subpd	%%xmm2		,%%xmm1\n\t"\
			"movaps	%%xmm4		,    (%%r15)\n\t"\
			"movaps	%%xmm5		,0x10(%%r15)\n\t"\
			"movaps	%%xmm0		,    (%%rcx)\n\t"\
			"movaps	%%xmm1		,0x10(%%rcx)\n\t"\
			"/* SSE2_RADIX_09_DIF_0TWIDDLE(r20,s1p18r,s1p14r,s1p10r,s1p06r,s1p02r,s1p34r,s1p30r,s1p26r,s1p22r) */\n\t"\
			"movq	%[__in0]	,%%rax 		/* __i0-8 = s1p[18,14,10,06,02,34,30,26,22]r; e[abc]x store output addresses throughout */\n\t"\
			"movq	%[__r00]	,%%rsi\n\t"\
			"addq	$0x240		,%%rax 		/* __i0 = s1p18r */\n\t"\
			"addq	$0x240		,%%rsi 		/* r20 */\n\t"\
		"movq	%[__cc1]	,%%rdx 		/* edx stores trig addresses throughout */\n\t"\
			"movq	%%rax		,%%r15\n\t"\
			"movq	%%rax		,%%rcx\n\t"\
			"subq	$0x180		,%%r15 		/* __i3 = s1p06r */\n\t"\
			"addq	$0x180		,%%rcx 		/* __i6 = s1p30r */\n\t"\
		"addq	$0x40		,%%rdx 		/* c3m1 */\n\t"\
			"movaps	    (%%r15)	,%%xmm2\n\t"\
			"movaps	0x10(%%r15)	,%%xmm3\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t"\
			"movaps	    (%%rcx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm7\n\t"\
			"movaps	%%xmm2		,%%xmm4\n\t"\
			"movaps	%%xmm3		,%%xmm5\n\t"\
		"movq	%%rsi		,%%rdi\n\t"\
		"addq	$0x20		,%%rdi\n\t"\
			"addpd	%%xmm6		,%%xmm2\n\t"\
			"addpd	%%xmm7		,%%xmm3\n\t"\
			"subpd	%%xmm6		,%%xmm4\n\t"\
			"subpd	%%xmm7		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm0\n\t"\
			"addpd	%%xmm3		,%%xmm1\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
			"movaps	%%xmm0		,    (%%rsi)	/* <- t00 */\n\t"\
			"movaps	%%xmm1		,0x10(%%rsi)	/* <- t01 */\n\t"\
		"addq	$0x40		,%%rsi\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t"\
			"mulpd	%%xmm7		,%%xmm4\n\t"\
			"mulpd	%%xmm7		,%%xmm5\n\t"\
			"addpd	%%xmm0		,%%xmm2\n\t"\
			"addpd	%%xmm1		,%%xmm3\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t"\
			"subpd	%%xmm5		,%%xmm2\n\t"\
			"addpd	%%xmm4		,%%xmm3\n\t"\
			"addpd	%%xmm5		,%%xmm0\n\t"\
			"subpd	%%xmm4		,%%xmm1\n\t"\
			"movaps	%%xmm2		,    (%%rdi)	/* <- t02 */\n\t"\
			"movaps	%%xmm3		,0x10(%%rdi)	/* <- t03 */\n\t"\
			"movaps	%%xmm0		,    (%%rsi)	/* <- t04 */\n\t"\
			"movaps	%%xmm1		,0x10(%%rsi)	/* <- t05 */\n\t"\
			"subq	$0x080		,%%rax 		/* __i1 = s1p14r */\n\t"\
			"subq	$0x080		,%%r15 		/* __i4 = s1p02r */\n\t"\
			"subq	$0x080		,%%rcx 		/* __i7 = s1p26r */\n\t"\
			"movaps	    (%%r15)	,%%xmm4\n\t"\
			"movaps	0x10(%%r15)	,%%xmm5\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t"\
			"movaps	    (%%rcx)	,%%xmm2\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm3\n\t"\
		"addq	$0x40		,%%rdi\n\t"\
			"subpd	%%xmm2		,%%xmm4\n\t"\
			"subpd	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm2\n\t"\
			"addpd	%%xmm3		,%%xmm3\n\t"\
			"addpd	%%xmm4		,%%xmm2\n\t"\
			"addpd	%%xmm5		,%%xmm3\n\t"\
			"addpd	%%xmm2		,%%xmm0\n\t"\
			"addpd	%%xmm3		,%%xmm1\n\t"\
			"movaps	%%xmm0		,    (%%rdi)	/* <- t06 */\n\t"\
			"movaps	%%xmm1		,0x10(%%rdi)	/* <- t07 */\n\t"\
		"addq	$0x40		,%%rsi\n\t"\
		"addq	$0x40		,%%rdi\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t"\
			"mulpd	%%xmm7		,%%xmm4\n\t"\
			"mulpd	%%xmm7		,%%xmm5\n\t"\
			"addpd	%%xmm0		,%%xmm2\n\t"\
			"addpd	%%xmm1		,%%xmm3\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t"\
			"subpd	%%xmm5		,%%xmm2\n\t"\
			"addpd	%%xmm4		,%%xmm3\n\t"\
			"addpd	%%xmm5		,%%xmm0\n\t"\
			"subpd	%%xmm4		,%%xmm1\n\t"\
			"movaps	%%xmm2		,    (%%rsi)	/* <- t08 */\n\t"\
			"movaps	%%xmm3		,0x10(%%rsi)	/* <- t09 */\n\t"\
			"movaps	%%xmm0		,    (%%rdi)	/* <- t0a */\n\t"\
			"movaps	%%xmm1		,0x10(%%rdi)	/* <- t0b */\n\t"\
		"addq	$0x40		,%%rsi\n\t"\
			"subq	$0x080		,%%rax 		/* __i2 = s1p10r */\n\t"\
			"addq	$0x400		,%%r15 		/* __i5 = s1p34r */\n\t"\
			"subq	$0x080		,%%rcx 		/* __i8 = s1p22r */\n\t"\
			"movaps	    (%%r15)	,%%xmm4\n\t"\
			"movaps	0x10(%%r15)	,%%xmm5\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t"\
			"movaps	    (%%rcx)	,%%xmm2\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm3\n\t"\
			"subpd	%%xmm2		,%%xmm4\n\t"\
			"subpd	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm2\n\t"\
			"addpd	%%xmm3		,%%xmm3\n\t"\
			"addpd	%%xmm4		,%%xmm2\n\t"\
			"addpd	%%xmm5		,%%xmm3\n\t"\
			"addpd	%%xmm2		,%%xmm0\n\t"\
			"addpd	%%xmm3		,%%xmm1\n\t"\
			"movaps	%%xmm0		,    (%%rsi)	/* <- t0c */\n\t"\
			"movaps	%%xmm1		,0x10(%%rsi)	/* <- t0d */\n\t"\
		"addq	$0x40		,%%rdi\n\t"\
		"addq	$0x40		,%%rsi\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t"\
			"mulpd	%%xmm7		,%%xmm4\n\t"\
			"mulpd	%%xmm7		,%%xmm5\n\t"\
			"addpd	%%xmm0		,%%xmm2\n\t"\
			"addpd	%%xmm1		,%%xmm3\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t"\
			"subpd	%%xmm5		,%%xmm2\n\t"\
			"addpd	%%xmm4		,%%xmm3\n\t"\
			"addpd	%%xmm5		,%%xmm0\n\t"\
			"subpd	%%xmm4		,%%xmm1\n\t"\
			"movaps	%%xmm2		,    (%%rdi)	/* <- t0e */\n\t"\
			"movaps	%%xmm3		,0x10(%%rdi)	/* <- t0f */\n\t"\
			"movaps	%%xmm0		,    (%%rsi)	/* <- t0g */\n\t"\
			"movaps	%%xmm1		,0x10(%%rsi)	/* <- t0h */\n\t"\
			"/******************************************************************************/\n\t"\
			"/*...and now do three more radix-3 transforms, including the twiddle factors: */\n\t"\
			"/******************************************************************************/\n\t"\
			"movq	%[__r00]	,%%rax 		/* __r00 */\n\t"\
			"addq	$0x240		,%%rax 		/* r20 */\n\t"\
			"movq	%%rax		,%%r15\n\t"\
			"movq	%%rax		,%%rcx\n\t"\
			"addq	$0x60		,%%r15 		/* __r26 */\n\t"\
			"addq	$0xc0		,%%rcx 		/* __r2c */\n\t"\
			"movaps	    (%%r15)	,%%xmm4\n\t"\
			"movaps	0x10(%%r15)	,%%xmm5\n\t"\
			"movaps	    (%%rcx)	,%%xmm2\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm3\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t"\
			"subpd	%%xmm2		,%%xmm4\n\t"\
			"subpd	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm2\n\t"\
			"addpd	%%xmm3		,%%xmm3\n\t"\
			"addpd	%%xmm4		,%%xmm2\n\t"\
			"addpd	%%xmm5		,%%xmm3\n\t"\
			"addpd	%%xmm2		,%%xmm0\n\t"\
			"addpd	%%xmm3		,%%xmm1\n\t"\
			"movaps	%%xmm0		,    (%%rax)\n\t"\
			"movaps	%%xmm1		,0x10(%%rax)\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t"\
			"mulpd	%%xmm7		,%%xmm4\n\t"\
			"mulpd	%%xmm7		,%%xmm5\n\t"\
			"addpd	%%xmm0		,%%xmm2\n\t"\
			"addpd	%%xmm1		,%%xmm3\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t"\
			"subpd	%%xmm5		,%%xmm2\n\t"\
			"addpd	%%xmm4		,%%xmm3\n\t"\
			"addpd	%%xmm5		,%%xmm0\n\t"\
			"subpd	%%xmm4		,%%xmm1\n\t"\
			"movaps	%%xmm2		,    (%%r15)\n\t"\
			"movaps	%%xmm3		,0x10(%%r15)\n\t"\
			"movaps	%%xmm0		,    (%%rcx)\n\t"\
			"movaps	%%xmm1		,0x10(%%rcx)\n\t"\
		"subq	$0x40		,%%rdx 		/* c1 */\n\t"\
			"addq	$0x20		,%%rax 		/* __r22 */\n\t"\
			"addq	$0x20		,%%r15 		/* __r28 */\n\t"\
			"addq	$0x20		,%%rcx 		/* __r2e */\n\t"\
			"movaps	    (%%r15)	,%%xmm2\n\t"\
			"movaps	0x10(%%r15)	,%%xmm3\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
		"addq	$0x20		,%%rdx 		/* c2 */\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"subpd	%%xmm1		,%%xmm2\n\t"\
			"addpd	%%xmm0		,%%xmm3\n\t"\
			"movaps	    (%%rcx)	,%%xmm4\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm5\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
		"addq	$0x20		,%%rdx 		/* c3m1 */\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"subpd	%%xmm1		,%%xmm4\n\t"\
			"addpd	%%xmm0		,%%xmm5\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t"\
			"subpd	%%xmm4		,%%xmm2\n\t"\
			"subpd	%%xmm5		,%%xmm3\n\t"\
			"addpd	%%xmm4		,%%xmm4\n\t"\
			"addpd	%%xmm5		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm4\n\t"\
			"addpd	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm4		,%%xmm0\n\t"\
			"addpd	%%xmm5		,%%xmm1\n\t"\
			"movaps	%%xmm0		,    (%%rax)\n\t"\
			"movaps	%%xmm1		,0x10(%%rax)\n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t"\
			"mulpd	%%xmm7		,%%xmm2\n\t"\
			"mulpd	%%xmm7		,%%xmm3\n\t"\
			"addpd	%%xmm0		,%%xmm4\n\t"\
			"addpd	%%xmm1		,%%xmm5\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t"\
			"subpd	%%xmm3		,%%xmm4\n\t"\
			"addpd	%%xmm2		,%%xmm5\n\t"\
			"addpd	%%xmm3		,%%xmm0\n\t"\
			"subpd	%%xmm2		,%%xmm1\n\t"\
			"movaps	%%xmm4		,    (%%r15)\n\t"\
			"movaps	%%xmm5		,0x10(%%r15)\n\t"\
			"movaps	%%xmm0		,    (%%rcx)\n\t"\
			"movaps	%%xmm1		,0x10(%%rcx)\n\t"\
		"subq	$0x20		,%%rdx 		/* c2 */\n\t"\
			"addq	$0x20		,%%rax 		/* __r24 */\n\t"\
			"addq	$0x20		,%%r15 		/* __r2a */\n\t"\
			"addq	$0x20		,%%rcx 		/* __r2g */\n\t"\
			"movaps	    (%%r15)	,%%xmm2\n\t"\
			"movaps	0x10(%%r15)	,%%xmm3\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
		"addq	$0x40		,%%rdx 		/* c4 */\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"subpd	%%xmm1		,%%xmm2\n\t"\
			"addpd	%%xmm0		,%%xmm3\n\t"\
			"movaps	    (%%rcx)	,%%xmm4\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm5\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
		"subq	$0x20		,%%rdx 		/* c3m1 */\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"subpd	%%xmm1		,%%xmm4\n\t"\
			"addpd	%%xmm0		,%%xmm5\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t"\
			"\n\t"\
			"subpd	%%xmm4		,%%xmm2\n\t"\
			"subpd	%%xmm5		,%%xmm3\n\t"\
			"addpd	%%xmm4		,%%xmm4\n\t"\
			"addpd	%%xmm5		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm4\n\t"\
			"addpd	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm4		,%%xmm0\n\t"\
			"addpd	%%xmm5		,%%xmm1\n\t"\
			"movaps	%%xmm0		,    (%%rax)\n\t"\
			"movaps	%%xmm1		,0x10(%%rax)\n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t"\
			"mulpd	%%xmm7		,%%xmm2\n\t"\
			"mulpd	%%xmm7		,%%xmm3\n\t"\
			"addpd	%%xmm0		,%%xmm4\n\t"\
			"addpd	%%xmm1		,%%xmm5\n\t"\
			"\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t"\
			"subpd	%%xmm3		,%%xmm4\n\t"\
			"addpd	%%xmm2		,%%xmm5\n\t"\
			"addpd	%%xmm3		,%%xmm0\n\t"\
			"subpd	%%xmm2		,%%xmm1\n\t"\
			"movaps	%%xmm4		,    (%%r15)\n\t"\
			"movaps	%%xmm5		,0x10(%%r15)\n\t"\
			"movaps	%%xmm0		,    (%%rcx)\n\t"\
			"movaps	%%xmm1		,0x10(%%rcx)\n\t"\
			"/* SSE2_RADIX_09_DIF_0TWIDDLE(r30,s1p09r,s1p05r,s1p01r,s1p33r,s1p29r,s1p25r,s1p21r,s1p17r,s1p13r) */\n\t"\
			"movq	%[__in0]	,%%rax 		/* __i0-8 = s1p[09,05,01,33,29,25,21,17,13]r; e[abc]x store output addresses throughout */\n\t"\
			"movq	%[__r00]	,%%rsi\n\t"\
			"addq	$0x120		,%%rax 		/* __i0 = s1p09r */\n\t"\
			"addq	$0x360		,%%rsi 		/* r30 */\n\t"\
		"movq	%[__cc1]	,%%rdx 		/* edx stores trig addresses throughout */\n\t"\
			"movq	%%rax		,%%r15\n\t"\
			"movq	%%rax		,%%rcx\n\t"\
			"addq	$0x300		,%%r15 		/* __i3 = s1p33r */\n\t"\
			"addq	$0x180		,%%rcx 		/* __i6 = s1p21r */\n\t"\
		"addq	$0x40		,%%rdx 		/* c3m1 */\n\t"\
			"movaps	    (%%r15)	,%%xmm2\n\t"\
			"movaps	0x10(%%r15)	,%%xmm3\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t"\
			"movaps	    (%%rcx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm7\n\t"\
			"movaps	%%xmm2		,%%xmm4\n\t"\
			"movaps	%%xmm3		,%%xmm5\n\t"\
		"movq	%%rsi		,%%rdi\n\t"\
		"addq	$0x20		,%%rdi\n\t"\
			"addpd	%%xmm6		,%%xmm2\n\t"\
			"addpd	%%xmm7		,%%xmm3\n\t"\
			"subpd	%%xmm6		,%%xmm4\n\t"\
			"subpd	%%xmm7		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm0\n\t"\
			"addpd	%%xmm3		,%%xmm1\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
			"movaps	%%xmm0		,    (%%rsi)	/* <- t00 */\n\t"\
			"movaps	%%xmm1		,0x10(%%rsi)	/* <- t01 */\n\t"\
		"addq	$0x40		,%%rsi\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t"\
			"mulpd	%%xmm7		,%%xmm4\n\t"\
			"mulpd	%%xmm7		,%%xmm5\n\t"\
			"addpd	%%xmm0		,%%xmm2\n\t"\
			"addpd	%%xmm1		,%%xmm3\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t"\
			"subpd	%%xmm5		,%%xmm2\n\t"\
			"addpd	%%xmm4		,%%xmm3\n\t"\
			"addpd	%%xmm5		,%%xmm0\n\t"\
			"subpd	%%xmm4		,%%xmm1\n\t"\
			"movaps	%%xmm2		,    (%%rdi)	/* <- t02 */\n\t"\
			"movaps	%%xmm3		,0x10(%%rdi)	/* <- t03 */\n\t"\
			"movaps	%%xmm0		,    (%%rsi)	/* <- t04 */\n\t"\
			"movaps	%%xmm1		,0x10(%%rsi)	/* <- t05 */\n\t"\
			"subq	$0x080		,%%rax 		/* __i1 = s1p05r */\n\t"\
			"subq	$0x080		,%%r15 		/* __i4 = s1p29r */\n\t"\
			"subq	$0x080		,%%rcx 		/* __i7 = s1p17r */\n\t"\
			"movaps	    (%%r15)	,%%xmm4\n\t"\
			"movaps	0x10(%%r15)	,%%xmm5\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t"\
			"movaps	    (%%rcx)	,%%xmm2\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm3\n\t"\
		"addq	$0x40		,%%rdi\n\t"\
			"subpd	%%xmm2		,%%xmm4\n\t"\
			"subpd	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm2\n\t"\
			"addpd	%%xmm3		,%%xmm3\n\t"\
			"addpd	%%xmm4		,%%xmm2\n\t"\
			"addpd	%%xmm5		,%%xmm3\n\t"\
			"addpd	%%xmm2		,%%xmm0\n\t"\
			"addpd	%%xmm3		,%%xmm1\n\t"\
			"movaps	%%xmm0		,    (%%rdi)	/* <- t06 */\n\t"\
			"movaps	%%xmm1		,0x10(%%rdi)	/* <- t07 */\n\t"\
		"addq	$0x40		,%%rsi\n\t"\
		"addq	$0x40		,%%rdi\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t"\
			"mulpd	%%xmm7		,%%xmm4\n\t"\
			"mulpd	%%xmm7		,%%xmm5\n\t"\
			"addpd	%%xmm0		,%%xmm2\n\t"\
			"addpd	%%xmm1		,%%xmm3\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t"\
			"subpd	%%xmm5		,%%xmm2\n\t"\
			"addpd	%%xmm4		,%%xmm3\n\t"\
			"addpd	%%xmm5		,%%xmm0\n\t"\
			"subpd	%%xmm4		,%%xmm1\n\t"\
			"movaps	%%xmm2		,    (%%rsi)	/* <- t08 */\n\t"\
			"movaps	%%xmm3		,0x10(%%rsi)	/* <- t09 */\n\t"\
			"movaps	%%xmm0		,    (%%rdi)	/* <- t0a */\n\t"\
			"movaps	%%xmm1		,0x10(%%rdi)	/* <- t0b */\n\t"\
		"addq	$0x40		,%%rsi\n\t"\
			"subq	$0x080		,%%rax 		/* __i2 = s1p01r */\n\t"\
			"subq	$0x080		,%%r15 		/* __i5 = s1p25r */\n\t"\
			"subq	$0x080		,%%rcx 		/* __i8 = s1p13r */\n\t"\
			"movaps	    (%%r15)	,%%xmm4\n\t"\
			"movaps	0x10(%%r15)	,%%xmm5\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t"\
			"movaps	    (%%rcx)	,%%xmm2\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm3\n\t"\
			"subpd	%%xmm2		,%%xmm4\n\t"\
			"subpd	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm2\n\t"\
			"addpd	%%xmm3		,%%xmm3\n\t"\
			"addpd	%%xmm4		,%%xmm2\n\t"\
			"addpd	%%xmm5		,%%xmm3\n\t"\
			"addpd	%%xmm2		,%%xmm0\n\t"\
			"addpd	%%xmm3		,%%xmm1\n\t"\
			"movaps	%%xmm0		,    (%%rsi)	/* <- t0c */\n\t"\
			"movaps	%%xmm1		,0x10(%%rsi)	/* <- t0d */\n\t"\
		"addq	$0x40		,%%rdi\n\t"\
		"addq	$0x40		,%%rsi\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t"\
			"mulpd	%%xmm7		,%%xmm4\n\t"\
			"mulpd	%%xmm7		,%%xmm5\n\t"\
			"addpd	%%xmm0		,%%xmm2\n\t"\
			"addpd	%%xmm1		,%%xmm3\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t"\
			"subpd	%%xmm5		,%%xmm2\n\t"\
			"addpd	%%xmm4		,%%xmm3\n\t"\
			"addpd	%%xmm5		,%%xmm0\n\t"\
			"subpd	%%xmm4		,%%xmm1\n\t"\
			"movaps	%%xmm2		,    (%%rdi)	/* <- t0e */\n\t"\
			"movaps	%%xmm3		,0x10(%%rdi)	/* <- t0f */\n\t"\
			"movaps	%%xmm0		,    (%%rsi)	/* <- t0g */\n\t"\
			"movaps	%%xmm1		,0x10(%%rsi)	/* <- t0h */\n\t"\
			"/******************************************************************************/\n\t"\
			"/*...and now do three more radix-3 transforms, including the twiddle factors: */\n\t"\
			"/******************************************************************************/\n\t"\
			"movq	%[__r00]	,%%rax 		/* __r00 */\n\t"\
			"addq	$0x360		,%%rax 		/* r30 */\n\t"\
			"movq	%%rax		,%%r15\n\t"\
			"movq	%%rax		,%%rcx\n\t"\
			"addq	$0x60		,%%r15 		/* __r36 */\n\t"\
			"addq	$0xc0		,%%rcx 		/* __r3c */\n\t"\
			"movaps	    (%%r15)	,%%xmm4\n\t"\
			"movaps	0x10(%%r15)	,%%xmm5\n\t"\
			"movaps	    (%%rcx)	,%%xmm2\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm3\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t"\
			"subpd	%%xmm2		,%%xmm4\n\t"\
			"subpd	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm2\n\t"\
			"addpd	%%xmm3		,%%xmm3\n\t"\
			"addpd	%%xmm4		,%%xmm2\n\t"\
			"addpd	%%xmm5		,%%xmm3\n\t"\
			"addpd	%%xmm2		,%%xmm0\n\t"\
			"addpd	%%xmm3		,%%xmm1\n\t"\
			"movaps	%%xmm0		,    (%%rax)\n\t"\
			"movaps	%%xmm1		,0x10(%%rax)\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t"\
			"mulpd	%%xmm7		,%%xmm4\n\t"\
			"mulpd	%%xmm7		,%%xmm5\n\t"\
			"addpd	%%xmm0		,%%xmm2\n\t"\
			"addpd	%%xmm1		,%%xmm3\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t"\
			"subpd	%%xmm5		,%%xmm2\n\t"\
			"addpd	%%xmm4		,%%xmm3\n\t"\
			"addpd	%%xmm5		,%%xmm0\n\t"\
			"subpd	%%xmm4		,%%xmm1\n\t"\
			"movaps	%%xmm2		,    (%%r15)\n\t"\
			"movaps	%%xmm3		,0x10(%%r15)\n\t"\
			"movaps	%%xmm0		,    (%%rcx)\n\t"\
			"movaps	%%xmm1		,0x10(%%rcx)\n\t"\
		"subq	$0x40		,%%rdx 		/* c1 */\n\t"\
			"addq	$0x20		,%%rax 		/* __r32 */\n\t"\
			"addq	$0x20		,%%r15 		/* __r38 */\n\t"\
			"addq	$0x20		,%%rcx 		/* __r3e */\n\t"\
			"movaps	    (%%r15)	,%%xmm2\n\t"\
			"movaps	0x10(%%r15)	,%%xmm3\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
		"addq	$0x20		,%%rdx 		/* c2 */\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"subpd	%%xmm1		,%%xmm2\n\t"\
			"addpd	%%xmm0		,%%xmm3\n\t"\
			"movaps	    (%%rcx)	,%%xmm4\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm5\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
		"addq	$0x20		,%%rdx 		/* c3m1 */\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"subpd	%%xmm1		,%%xmm4\n\t"\
			"addpd	%%xmm0		,%%xmm5\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t"\
			"subpd	%%xmm4		,%%xmm2\n\t"\
			"subpd	%%xmm5		,%%xmm3\n\t"\
			"addpd	%%xmm4		,%%xmm4\n\t"\
			"addpd	%%xmm5		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm4\n\t"\
			"addpd	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm4		,%%xmm0\n\t"\
			"addpd	%%xmm5		,%%xmm1\n\t"\
			"movaps	%%xmm0		,    (%%rax)\n\t"\
			"movaps	%%xmm1		,0x10(%%rax)\n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t"\
			"mulpd	%%xmm7		,%%xmm2\n\t"\
			"mulpd	%%xmm7		,%%xmm3\n\t"\
			"addpd	%%xmm0		,%%xmm4\n\t"\
			"addpd	%%xmm1		,%%xmm5\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t"\
			"subpd	%%xmm3		,%%xmm4\n\t"\
			"addpd	%%xmm2		,%%xmm5\n\t"\
			"addpd	%%xmm3		,%%xmm0\n\t"\
			"subpd	%%xmm2		,%%xmm1\n\t"\
			"movaps	%%xmm4		,    (%%r15)\n\t"\
			"movaps	%%xmm5		,0x10(%%r15)\n\t"\
			"movaps	%%xmm0		,    (%%rcx)\n\t"\
			"movaps	%%xmm1		,0x10(%%rcx)\n\t"\
		"subq	$0x20		,%%rdx 		/* c2 */\n\t"\
			"addq	$0x20		,%%rax 		/* __r34 */\n\t"\
			"addq	$0x20		,%%r15 		/* __r3a */\n\t"\
			"addq	$0x20		,%%rcx 		/* __r3g */\n\t"\
			"movaps	    (%%r15)	,%%xmm2\n\t"\
			"movaps	0x10(%%r15)	,%%xmm3\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
		"addq	$0x40		,%%rdx 		/* c4 */\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"subpd	%%xmm1		,%%xmm2\n\t"\
			"addpd	%%xmm0		,%%xmm3\n\t"\
			"movaps	    (%%rcx)	,%%xmm4\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm5\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
		"subq	$0x20		,%%rdx 		/* c3m1 */\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"subpd	%%xmm1		,%%xmm4\n\t"\
			"addpd	%%xmm0		,%%xmm5\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t"\
			"\n\t"\
			"subpd	%%xmm4		,%%xmm2\n\t"\
			"subpd	%%xmm5		,%%xmm3\n\t"\
			"addpd	%%xmm4		,%%xmm4\n\t"\
			"addpd	%%xmm5		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm4\n\t"\
			"addpd	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm4		,%%xmm0\n\t"\
			"addpd	%%xmm5		,%%xmm1\n\t"\
			"movaps	%%xmm0		,    (%%rax)\n\t"\
			"movaps	%%xmm1		,0x10(%%rax)\n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t"\
			"mulpd	%%xmm7		,%%xmm2\n\t"\
			"mulpd	%%xmm7		,%%xmm3\n\t"\
			"addpd	%%xmm0		,%%xmm4\n\t"\
			"addpd	%%xmm1		,%%xmm5\n\t"\
			"\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t"\
			"subpd	%%xmm3		,%%xmm4\n\t"\
			"addpd	%%xmm2		,%%xmm5\n\t"\
			"addpd	%%xmm3		,%%xmm0\n\t"\
			"subpd	%%xmm2		,%%xmm1\n\t"\
			"movaps	%%xmm4		,    (%%r15)\n\t"\
			"movaps	%%xmm5		,0x10(%%r15)\n\t"\
			"movaps	%%xmm0		,    (%%rcx)\n\t"\
			"movaps	%%xmm1		,0x10(%%rcx)\n\t"\
			"\n\t"\
		"/**********************************/"\
		"/*** And now do 9 radix-4 DFTs: ***/"\
		"/**********************************/"\
			"\n\t"\
			"/* add0,1,3,2 = &a[j1+p00]+p0,1,2,3 */		\n\t"\
			"movq	%[__add],%%rax	/* Must use eax as base address throughout, since that is preserved in SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_B */\n\t"\
			"movq	%[__add],%%rdx						\n\t"\
			"movslq	%[__p01],%%rsi	/* esi will store power-of-2 multiples of p01 throughout */\n\t"\
			"shlq	$3,%%rsi		/* Pointer offset for floating doubles */\n\t"\
			"addq	%%rsi,%%rdx							\n\t"\
			"movq	%%rdx,%%r15		/* add0+p01 */		\n\t"\
			"addq	%%rsi,%%rdx							\n\t"\
			"movq	%%rdx,%%rcx		/* add0+p02 */		\n\t"\
			"addq	%%rsi,%%rdx		/* add0+p03 */		\n\t"\
			"/* SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_C(r00, 0x120, 0x240, eax,ebx,edx,ecx) */\n\t"\
			"movq	%[__r00],%%rdi 	/* edi = __r00 */	\n\t"\
			"movq	$0x120	,%%rsi 						\n\t"\
			"addq	%%rdi	,%%rsi 	/* esi = __r00 + 0x120 */\n\t"\
			"movaps	    (%%rdi),%%xmm4					\n\t"\
			"movaps	    (%%rsi),%%xmm6					\n\t"\
			"movaps	0x10(%%rdi),%%xmm5					\n\t"\
			"movaps	0x10(%%rsi),%%xmm7					\n\t"\
			"addq	$0x240,%%rdi 						\n\t"\
			"addq	$0x240,%%rsi 						\n\t"\
			"movaps	    (%%rdi),%%xmm0					\n\t"\
			"movaps	    (%%rsi),%%xmm2					\n\t"\
			"movaps	0x10(%%rdi),%%xmm1					\n\t"\
			"movaps	0x10(%%rsi),%%xmm3					\n\t"\
			"subpd	%%xmm0,%%xmm4						\n\t"\
			"subpd	%%xmm2,%%xmm6						\n\t"\
			"subpd	%%xmm1,%%xmm5						\n\t"\
			"subpd	%%xmm3,%%xmm7						\n\t"\
			"addpd	%%xmm0,%%xmm0						\n\t"\
			"addpd	%%xmm2,%%xmm2						\n\t"\
			"addpd	%%xmm1,%%xmm1						\n\t"\
			"addpd	%%xmm3,%%xmm3						\n\t"\
			"addpd	%%xmm4,%%xmm0						\n\t"\
			"addpd	%%xmm6,%%xmm2						\n\t"\
			"addpd	%%xmm5,%%xmm1						\n\t"\
			"addpd	%%xmm7,%%xmm3						\n\t"\
			"/* Swap output addresses d <-> c: */		\n\t"\
			"subpd	%%xmm2,%%xmm0						\n\t"\
			"subpd	%%xmm7,%%xmm4						\n\t"\
			"subpd	%%xmm3,%%xmm1						\n\t"\
			"subpd	%%xmm6,%%xmm5						\n\t"\
			"movaps	%%xmm0,    (%%r15)					\n\t"\
			"movaps	%%xmm4,    (%%rdx)					\n\t"\
			"movaps	%%xmm1,0x10(%%r15)					\n\t"\
			"movaps	%%xmm5,0x10(%%rcx)					\n\t"\
			"addpd	%%xmm2,%%xmm2						\n\t"\
			"addpd	%%xmm7,%%xmm7						\n\t"\
			"addpd	%%xmm3,%%xmm3						\n\t"\
			"addpd	%%xmm6,%%xmm6						\n\t"\
			"addpd	%%xmm0,%%xmm2						\n\t"\
			"addpd	%%xmm4,%%xmm7						\n\t"\
			"addpd	%%xmm1,%%xmm3						\n\t"\
			"addpd	%%xmm5,%%xmm6						\n\t"\
			"movaps	%%xmm2,    (%%rax)					\n\t"\
			"movaps	%%xmm7,    (%%rcx)					\n\t"\
			"movaps	%%xmm3,0x10(%%rax)					\n\t"\
			"movaps	%%xmm6,0x10(%%rdx)					\n\t"\
			"\n\t"\
			"/* add0,1,3,2 = &a[j1+p32]+p0,1,2,3 */		\n\t"\
			"movslq	%[__p32],%%rsi						\n\t"\
			"shlq	$3,%%rsi							\n\t"\
			"addq	%%rsi,%%rax		/* &a[j1+p32] */	\n\t"\
			"addq	%%rsi,%%r15							\n\t"\
			"addq	%%rsi,%%rcx							\n\t"\
			"addq	%%rsi,%%rdx							\n\t"\
			"/* SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_C(r02, 0x120, 0x240, eax,ebx,edx,ecx) */\n\t"\
			"subq	$0x220	,%%rdi 	/* edi = __r02 */	\n\t"\
			"movq	$0x120	,%%rsi 						\n\t"\
			"addq	%%rdi	,%%rsi 	/* esi = __r02 + 0x120 */\n\t"\
			"movaps	    (%%rdi),%%xmm4					\n\t"\
			"movaps	    (%%rsi),%%xmm6					\n\t"\
			"movaps	0x10(%%rdi),%%xmm5					\n\t"\
			"movaps	0x10(%%rsi),%%xmm7					\n\t"\
			"addq	$0x240,%%rdi 						\n\t"\
			"addq	$0x240,%%rsi 						\n\t"\
			"movaps	    (%%rdi),%%xmm0					\n\t"\
			"movaps	    (%%rsi),%%xmm2					\n\t"\
			"movaps	0x10(%%rdi),%%xmm1					\n\t"\
			"movaps	0x10(%%rsi),%%xmm3					\n\t"\
			"subpd	%%xmm0,%%xmm4						\n\t"\
			"subpd	%%xmm2,%%xmm6						\n\t"\
			"subpd	%%xmm1,%%xmm5						\n\t"\
			"subpd	%%xmm3,%%xmm7						\n\t"\
			"addpd	%%xmm0,%%xmm0						\n\t"\
			"addpd	%%xmm2,%%xmm2						\n\t"\
			"addpd	%%xmm1,%%xmm1						\n\t"\
			"addpd	%%xmm3,%%xmm3						\n\t"\
			"addpd	%%xmm4,%%xmm0						\n\t"\
			"addpd	%%xmm6,%%xmm2						\n\t"\
			"addpd	%%xmm5,%%xmm1						\n\t"\
			"addpd	%%xmm7,%%xmm3						\n\t"\
			"/* Swap output addresses d <-> c: */		\n\t"\
			"subpd	%%xmm2,%%xmm0						\n\t"\
			"subpd	%%xmm7,%%xmm4						\n\t"\
			"subpd	%%xmm3,%%xmm1						\n\t"\
			"subpd	%%xmm6,%%xmm5						\n\t"\
			"movaps	%%xmm0,    (%%r15)					\n\t"\
			"movaps	%%xmm4,    (%%rdx)					\n\t"\
			"movaps	%%xmm1,0x10(%%r15)					\n\t"\
			"movaps	%%xmm5,0x10(%%rcx)					\n\t"\
			"addpd	%%xmm2,%%xmm2						\n\t"\
			"addpd	%%xmm7,%%xmm7						\n\t"\
			"addpd	%%xmm3,%%xmm3						\n\t"\
			"addpd	%%xmm6,%%xmm6						\n\t"\
			"addpd	%%xmm0,%%xmm2						\n\t"\
			"addpd	%%xmm4,%%xmm7						\n\t"\
			"addpd	%%xmm1,%%xmm3						\n\t"\
			"addpd	%%xmm5,%%xmm6						\n\t"\
			"movaps	%%xmm2,    (%%rax)					\n\t"\
			"movaps	%%xmm7,    (%%rcx)					\n\t"\
			"movaps	%%xmm3,0x10(%%rax)					\n\t"\
			"movaps	%%xmm6,0x10(%%rdx)					\n\t"\
			"\n\t"\
			"/* add2,3,0,1 = &a[j1+p20]+p0,1,2,3 */		\n\t"\
			"movslq	%[__p12],%%rsi						\n\t"\
			"shlq	$3,%%rsi							\n\t"\
			"subq	%%rsi,%%rax		/* &a[j1+p20] */	\n\t"\
			"subq	%%rsi,%%r15							\n\t"\
			"subq	%%rsi,%%rcx							\n\t"\
			"subq	%%rsi,%%rdx							\n\t"\
			"/* SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_C(r04, 0x120, 0x240, ecx,edx,eax,ebx) */\n\t"\
			"subq	$0x220	,%%rdi 	/* edi = __r04 */	\n\t"\
			"movq	$0x120	,%%rsi 						\n\t"\
			"addq	%%rdi	,%%rsi 	/* esi = __r04 + 0x120 */\n\t"\
			"movaps	    (%%rdi),%%xmm4					\n\t"\
			"movaps	    (%%rsi),%%xmm6					\n\t"\
			"movaps	0x10(%%rdi),%%xmm5					\n\t"\
			"movaps	0x10(%%rsi),%%xmm7					\n\t"\
			"addq	$0x240,%%rdi 						\n\t"\
			"addq	$0x240,%%rsi 						\n\t"\
			"movaps	    (%%rdi),%%xmm0					\n\t"\
			"movaps	    (%%rsi),%%xmm2					\n\t"\
			"movaps	0x10(%%rdi),%%xmm1					\n\t"\
			"movaps	0x10(%%rsi),%%xmm3					\n\t"\
			"subpd	%%xmm0,%%xmm4						\n\t"\
			"subpd	%%xmm2,%%xmm6						\n\t"\
			"subpd	%%xmm1,%%xmm5						\n\t"\
			"subpd	%%xmm3,%%xmm7						\n\t"\
			"addpd	%%xmm0,%%xmm0						\n\t"\
			"addpd	%%xmm2,%%xmm2						\n\t"\
			"addpd	%%xmm1,%%xmm1						\n\t"\
			"addpd	%%xmm3,%%xmm3						\n\t"\
			"addpd	%%xmm4,%%xmm0						\n\t"\
			"addpd	%%xmm6,%%xmm2						\n\t"\
			"addpd	%%xmm5,%%xmm1						\n\t"\
			"addpd	%%xmm7,%%xmm3						\n\t"\
			"/* Swap output addresses a <-> c ,b <-> d: */\n\t"\
			"subpd	%%xmm2,%%xmm0						\n\t"\
			"subpd	%%xmm7,%%xmm4						\n\t"\
			"subpd	%%xmm3,%%xmm1						\n\t"\
			"subpd	%%xmm6,%%xmm5						\n\t"\
			"movaps	%%xmm0,    (%%rdx)					\n\t"\
			"movaps	%%xmm4,    (%%rax)					\n\t"\
			"movaps	%%xmm1,0x10(%%rdx)					\n\t"\
			"movaps	%%xmm5,0x10(%%r15)					\n\t"\
			"addpd	%%xmm2,%%xmm2						\n\t"\
			"addpd	%%xmm7,%%xmm7						\n\t"\
			"addpd	%%xmm3,%%xmm3						\n\t"\
			"addpd	%%xmm6,%%xmm6						\n\t"\
			"addpd	%%xmm0,%%xmm2						\n\t"\
			"addpd	%%xmm4,%%xmm7						\n\t"\
			"addpd	%%xmm1,%%xmm3						\n\t"\
			"addpd	%%xmm5,%%xmm6						\n\t"\
			"movaps	%%xmm2,    (%%rcx)					\n\t"\
			"movaps	%%xmm7,    (%%r15)					\n\t"\
			"movaps	%%xmm3,0x10(%%rcx)					\n\t"\
			"movaps	%%xmm6,0x10(%%rax)					\n\t"\
			"\n\t"\
			"/* add1,0,2,3 = &a[j1+p08]+p0,1,2,3 */		\n\t"\
			"movslq	%[__p12],%%rsi						\n\t"\
			"shlq	$3,%%rsi							\n\t"\
			"subq	%%rsi,%%rax		/* &a[j1+p08] */	\n\t"\
			"subq	%%rsi,%%r15							\n\t"\
			"subq	%%rsi,%%rcx							\n\t"\
			"subq	%%rsi,%%rdx							\n\t"\
			"/* SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_C(r06, 0x120, 0x240, ebx,eax,ecx,edx) */\n\t"\
			"subq	$0x220	,%%rdi 	/* edi = __r06 */	\n\t"\
			"movq	$0x120	,%%rsi 						\n\t"\
			"addq	%%rdi	,%%rsi 	/* esi = __r06 + 0x120 */\n\t"\
			"movaps	    (%%rdi),%%xmm4					\n\t"\
			"movaps	    (%%rsi),%%xmm6					\n\t"\
			"movaps	0x10(%%rdi),%%xmm5					\n\t"\
			"movaps	0x10(%%rsi),%%xmm7					\n\t"\
			"addq	$0x240,%%rdi 						\n\t"\
			"addq	$0x240,%%rsi 						\n\t"\
			"movaps	    (%%rdi),%%xmm0					\n\t"\
			"movaps	    (%%rsi),%%xmm2					\n\t"\
			"movaps	0x10(%%rdi),%%xmm1					\n\t"\
			"movaps	0x10(%%rsi),%%xmm3					\n\t"\
			"subpd	%%xmm0,%%xmm4						\n\t"\
			"subpd	%%xmm2,%%xmm6						\n\t"\
			"subpd	%%xmm1,%%xmm5						\n\t"\
			"subpd	%%xmm3,%%xmm7						\n\t"\
			"addpd	%%xmm0,%%xmm0						\n\t"\
			"addpd	%%xmm2,%%xmm2						\n\t"\
			"addpd	%%xmm1,%%xmm1						\n\t"\
			"addpd	%%xmm3,%%xmm3						\n\t"\
			"addpd	%%xmm4,%%xmm0						\n\t"\
			"addpd	%%xmm6,%%xmm2						\n\t"\
			"addpd	%%xmm5,%%xmm1						\n\t"\
			"addpd	%%xmm7,%%xmm3						\n\t"\
			"/* Swap output addresses a <-> b: */		\n\t"\
			"subpd	%%xmm2,%%xmm0						\n\t"\
			"subpd	%%xmm7,%%xmm4						\n\t"\
			"subpd	%%xmm3,%%xmm1						\n\t"\
			"subpd	%%xmm6,%%xmm5						\n\t"\
			"movaps	%%xmm0,    (%%rax)					\n\t"\
			"movaps	%%xmm4,    (%%rcx)					\n\t"\
			"movaps	%%xmm1,0x10(%%rax)					\n\t"\
			"movaps	%%xmm5,0x10(%%rdx)					\n\t"\
			"addpd	%%xmm2,%%xmm2						\n\t"\
			"addpd	%%xmm7,%%xmm7						\n\t"\
			"addpd	%%xmm3,%%xmm3						\n\t"\
			"addpd	%%xmm6,%%xmm6						\n\t"\
			"addpd	%%xmm0,%%xmm2						\n\t"\
			"addpd	%%xmm4,%%xmm7						\n\t"\
			"addpd	%%xmm1,%%xmm3						\n\t"\
			"addpd	%%xmm5,%%xmm6						\n\t"\
			"movaps	%%xmm2,    (%%r15)					\n\t"\
			"movaps	%%xmm7,    (%%rdx)					\n\t"\
			"movaps	%%xmm3,0x10(%%r15)					\n\t"\
			"movaps	%%xmm6,0x10(%%rcx)					\n\t"\
			"\n\t"\
			"/* add3,2,1,0 = &a[j1+p28]+p0,1,2,3 */		\n\t"\
			"movslq	%[__p20],%%rsi						\n\t"\
			"shlq	$3,%%rsi							\n\t"\
			"addq	%%rsi,%%rax		/* &a[j1+p28] */	\n\t"\
			"addq	%%rsi,%%r15							\n\t"\
			"addq	%%rsi,%%rcx							\n\t"\
			"addq	%%rsi,%%rdx							\n\t"\
			"/* SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_C(r08, 0x120, 0x240, edx,ecx,ebx,eax) */\n\t"\
			"subq	$0x220	,%%rdi 	/* edi = __r08 */	\n\t"\
			"movq	$0x120	,%%rsi 						\n\t"\
			"addq	%%rdi	,%%rsi 	/* esi = __r08 + 0x120 */\n\t"\
			"movaps	    (%%rdi),%%xmm4					\n\t"\
			"movaps	    (%%rsi),%%xmm6					\n\t"\
			"movaps	0x10(%%rdi),%%xmm5					\n\t"\
			"movaps	0x10(%%rsi),%%xmm7					\n\t"\
			"addq	$0x240,%%rdi 						\n\t"\
			"addq	$0x240,%%rsi 						\n\t"\
			"movaps	    (%%rdi),%%xmm0					\n\t"\
			"movaps	    (%%rsi),%%xmm2					\n\t"\
			"movaps	0x10(%%rdi),%%xmm1					\n\t"\
			"movaps	0x10(%%rsi),%%xmm3					\n\t"\
			"subpd	%%xmm0,%%xmm4						\n\t"\
			"subpd	%%xmm2,%%xmm6						\n\t"\
			"subpd	%%xmm1,%%xmm5						\n\t"\
			"subpd	%%xmm3,%%xmm7						\n\t"\
			"addpd	%%xmm0,%%xmm0						\n\t"\
			"addpd	%%xmm2,%%xmm2						\n\t"\
			"addpd	%%xmm1,%%xmm1						\n\t"\
			"addpd	%%xmm3,%%xmm3						\n\t"\
			"addpd	%%xmm4,%%xmm0						\n\t"\
			"addpd	%%xmm6,%%xmm2						\n\t"\
			"addpd	%%xmm5,%%xmm1						\n\t"\
			"addpd	%%xmm7,%%xmm3						\n\t"\
			"/* Swap output addresses a <-> d ,b <-> c: */\n\t"\
			"subpd	%%xmm2,%%xmm0						\n\t"\
			"subpd	%%xmm7,%%xmm4						\n\t"\
			"subpd	%%xmm3,%%xmm1						\n\t"\
			"subpd	%%xmm6,%%xmm5						\n\t"\
			"movaps	%%xmm0,    (%%rcx)					\n\t"\
			"movaps	%%xmm4,    (%%r15)					\n\t"\
			"movaps	%%xmm1,0x10(%%rcx)					\n\t"\
			"movaps	%%xmm5,0x10(%%rax)					\n\t"\
			"addpd	%%xmm2,%%xmm2						\n\t"\
			"addpd	%%xmm7,%%xmm7						\n\t"\
			"addpd	%%xmm3,%%xmm3						\n\t"\
			"addpd	%%xmm6,%%xmm6						\n\t"\
			"addpd	%%xmm0,%%xmm2						\n\t"\
			"addpd	%%xmm4,%%xmm7						\n\t"\
			"addpd	%%xmm1,%%xmm3						\n\t"\
			"addpd	%%xmm5,%%xmm6						\n\t"\
			"movaps	%%xmm2,    (%%rdx)					\n\t"\
			"movaps	%%xmm7,    (%%rax)					\n\t"\
			"movaps	%%xmm3,0x10(%%rdx)					\n\t"\
			"movaps	%%xmm6,0x10(%%r15)					\n\t"\
			"\n\t"\
			"/* add0,1,3,2 = &a[j1+p16]+p0,1,2,3 */		\n\t"\
			"movslq	%[__p12],%%rsi						\n\t"\
			"shlq	$3,%%rsi							\n\t"\
			"subq	%%rsi,%%rax		/* &a[j1+p16] */	\n\t"\
			"subq	%%rsi,%%r15							\n\t"\
			"subq	%%rsi,%%rcx							\n\t"\
			"subq	%%rsi,%%rdx							\n\t"\
			"/* SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_C(r0a, 0x120, 0x240, eax,ebx,edx,ecx) */\n\t"\
			"subq	$0x220	,%%rdi 	/* edi = __r0a */	\n\t"\
			"movq	$0x120	,%%rsi 						\n\t"\
			"addq	%%rdi	,%%rsi 	/* esi = __r0a + 0x120 */\n\t"\
			"movaps	    (%%rdi),%%xmm4					\n\t"\
			"movaps	    (%%rsi),%%xmm6					\n\t"\
			"movaps	0x10(%%rdi),%%xmm5					\n\t"\
			"movaps	0x10(%%rsi),%%xmm7					\n\t"\
			"addq	$0x240,%%rdi 						\n\t"\
			"addq	$0x240,%%rsi 						\n\t"\
			"movaps	    (%%rdi),%%xmm0					\n\t"\
			"movaps	    (%%rsi),%%xmm2					\n\t"\
			"movaps	0x10(%%rdi),%%xmm1					\n\t"\
			"movaps	0x10(%%rsi),%%xmm3					\n\t"\
			"subpd	%%xmm0,%%xmm4						\n\t"\
			"subpd	%%xmm2,%%xmm6						\n\t"\
			"subpd	%%xmm1,%%xmm5						\n\t"\
			"subpd	%%xmm3,%%xmm7						\n\t"\
			"addpd	%%xmm0,%%xmm0						\n\t"\
			"addpd	%%xmm2,%%xmm2						\n\t"\
			"addpd	%%xmm1,%%xmm1						\n\t"\
			"addpd	%%xmm3,%%xmm3						\n\t"\
			"addpd	%%xmm4,%%xmm0						\n\t"\
			"addpd	%%xmm6,%%xmm2						\n\t"\
			"addpd	%%xmm5,%%xmm1						\n\t"\
			"addpd	%%xmm7,%%xmm3						\n\t"\
			"/* Swap output addresses d <-> c: */		\n\t"\
			"subpd	%%xmm2,%%xmm0						\n\t"\
			"subpd	%%xmm7,%%xmm4						\n\t"\
			"subpd	%%xmm3,%%xmm1						\n\t"\
			"subpd	%%xmm6,%%xmm5						\n\t"\
			"movaps	%%xmm0,    (%%r15)					\n\t"\
			"movaps	%%xmm4,    (%%rdx)					\n\t"\
			"movaps	%%xmm1,0x10(%%r15)					\n\t"\
			"movaps	%%xmm5,0x10(%%rcx)					\n\t"\
			"addpd	%%xmm2,%%xmm2						\n\t"\
			"addpd	%%xmm7,%%xmm7						\n\t"\
			"addpd	%%xmm3,%%xmm3						\n\t"\
			"addpd	%%xmm6,%%xmm6						\n\t"\
			"addpd	%%xmm0,%%xmm2						\n\t"\
			"addpd	%%xmm4,%%xmm7						\n\t"\
			"addpd	%%xmm1,%%xmm3						\n\t"\
			"addpd	%%xmm5,%%xmm6						\n\t"\
			"movaps	%%xmm2,    (%%rax)					\n\t"\
			"movaps	%%xmm7,    (%%rcx)					\n\t"\
			"movaps	%%xmm3,0x10(%%rax)					\n\t"\
			"movaps	%%xmm6,0x10(%%rdx)					\n\t"\
			"\n\t"\
			"/* add2,3,0,1 = &a[j1+p04]+p0,1,2,3 */		\n\t"\
			"movslq	%[__p12],%%rsi						\n\t"\
			"shlq	$3,%%rsi							\n\t"\
			"subq	%%rsi,%%rax		/* &a[j1+p04] */	\n\t"\
			"subq	%%rsi,%%r15							\n\t"\
			"subq	%%rsi,%%rcx							\n\t"\
			"subq	%%rsi,%%rdx							\n\t"\
			"/* SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_C(r0c, 0x120, 0x240, ecx,edx,eax,ebx) */\n\t"\
			"subq	$0x220	,%%rdi 	/* edi = __r0c */	\n\t"\
			"movq	$0x120	,%%rsi 						\n\t"\
			"addq	%%rdi	,%%rsi 	/* esi = __r0c + 0x120 */\n\t"\
			"movaps	    (%%rdi),%%xmm4					\n\t"\
			"movaps	    (%%rsi),%%xmm6					\n\t"\
			"movaps	0x10(%%rdi),%%xmm5					\n\t"\
			"movaps	0x10(%%rsi),%%xmm7					\n\t"\
			"addq	$0x240,%%rdi 						\n\t"\
			"addq	$0x240,%%rsi 						\n\t"\
			"movaps	    (%%rdi),%%xmm0					\n\t"\
			"movaps	    (%%rsi),%%xmm2					\n\t"\
			"movaps	0x10(%%rdi),%%xmm1					\n\t"\
			"movaps	0x10(%%rsi),%%xmm3					\n\t"\
			"subpd	%%xmm0,%%xmm4						\n\t"\
			"subpd	%%xmm2,%%xmm6						\n\t"\
			"subpd	%%xmm1,%%xmm5						\n\t"\
			"subpd	%%xmm3,%%xmm7						\n\t"\
			"addpd	%%xmm0,%%xmm0						\n\t"\
			"addpd	%%xmm2,%%xmm2						\n\t"\
			"addpd	%%xmm1,%%xmm1						\n\t"\
			"addpd	%%xmm3,%%xmm3						\n\t"\
			"addpd	%%xmm4,%%xmm0						\n\t"\
			"addpd	%%xmm6,%%xmm2						\n\t"\
			"addpd	%%xmm5,%%xmm1						\n\t"\
			"addpd	%%xmm7,%%xmm3						\n\t"\
			"/* Swap output addresses a <-> c ,b <-> d: */\n\t"\
			"subpd	%%xmm2,%%xmm0						\n\t"\
			"subpd	%%xmm7,%%xmm4						\n\t"\
			"subpd	%%xmm3,%%xmm1						\n\t"\
			"subpd	%%xmm6,%%xmm5						\n\t"\
			"movaps	%%xmm0,    (%%rdx)					\n\t"\
			"movaps	%%xmm4,    (%%rax)					\n\t"\
			"movaps	%%xmm1,0x10(%%rdx)					\n\t"\
			"movaps	%%xmm5,0x10(%%r15)					\n\t"\
			"addpd	%%xmm2,%%xmm2						\n\t"\
			"addpd	%%xmm7,%%xmm7						\n\t"\
			"addpd	%%xmm3,%%xmm3						\n\t"\
			"addpd	%%xmm6,%%xmm6						\n\t"\
			"addpd	%%xmm0,%%xmm2						\n\t"\
			"addpd	%%xmm4,%%xmm7						\n\t"\
			"addpd	%%xmm1,%%xmm3						\n\t"\
			"addpd	%%xmm5,%%xmm6						\n\t"\
			"movaps	%%xmm2,    (%%rcx)					\n\t"\
			"movaps	%%xmm7,    (%%r15)					\n\t"\
			"movaps	%%xmm3,0x10(%%rcx)					\n\t"\
			"movaps	%%xmm6,0x10(%%rax)					\n\t"\
			"\n\t"\
			"/* add1,0,2,3 = &a[j1+p24]+p0,1,2,3 */		\n\t"\
			"movslq	%[__p20],%%rsi						\n\t"\
			"shlq	$3,%%rsi							\n\t"\
			"addq	%%rsi,%%rax		/* &a[j1+p24] */	\n\t"\
			"addq	%%rsi,%%r15							\n\t"\
			"addq	%%rsi,%%rcx							\n\t"\
			"addq	%%rsi,%%rdx							\n\t"\
			"/* SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_C(r0e, 0x120, 0x240, ebx,eax,ecx,edx) */\n\t"\
			"subq	$0x220	,%%rdi 	/* edi = __r0e */	\n\t"\
			"movq	$0x120	,%%rsi 						\n\t"\
			"addq	%%rdi	,%%rsi 	/* esi = __r0e + 0x120 */\n\t"\
			"movaps	    (%%rdi),%%xmm4					\n\t"\
			"movaps	    (%%rsi),%%xmm6					\n\t"\
			"movaps	0x10(%%rdi),%%xmm5					\n\t"\
			"movaps	0x10(%%rsi),%%xmm7					\n\t"\
			"addq	$0x240,%%rdi 						\n\t"\
			"addq	$0x240,%%rsi 						\n\t"\
			"movaps	    (%%rdi),%%xmm0					\n\t"\
			"movaps	    (%%rsi),%%xmm2					\n\t"\
			"movaps	0x10(%%rdi),%%xmm1					\n\t"\
			"movaps	0x10(%%rsi),%%xmm3					\n\t"\
			"subpd	%%xmm0,%%xmm4						\n\t"\
			"subpd	%%xmm2,%%xmm6						\n\t"\
			"subpd	%%xmm1,%%xmm5						\n\t"\
			"subpd	%%xmm3,%%xmm7						\n\t"\
			"addpd	%%xmm0,%%xmm0						\n\t"\
			"addpd	%%xmm2,%%xmm2						\n\t"\
			"addpd	%%xmm1,%%xmm1						\n\t"\
			"addpd	%%xmm3,%%xmm3						\n\t"\
			"addpd	%%xmm4,%%xmm0						\n\t"\
			"addpd	%%xmm6,%%xmm2						\n\t"\
			"addpd	%%xmm5,%%xmm1						\n\t"\
			"addpd	%%xmm7,%%xmm3						\n\t"\
			"/* Swap output addresses a <-> b: */		\n\t"\
			"subpd	%%xmm2,%%xmm0						\n\t"\
			"subpd	%%xmm7,%%xmm4						\n\t"\
			"subpd	%%xmm3,%%xmm1						\n\t"\
			"subpd	%%xmm6,%%xmm5						\n\t"\
			"movaps	%%xmm0,    (%%rax)					\n\t"\
			"movaps	%%xmm4,    (%%rcx)					\n\t"\
			"movaps	%%xmm1,0x10(%%rax)					\n\t"\
			"movaps	%%xmm5,0x10(%%rdx)					\n\t"\
			"addpd	%%xmm2,%%xmm2						\n\t"\
			"addpd	%%xmm7,%%xmm7						\n\t"\
			"addpd	%%xmm3,%%xmm3						\n\t"\
			"addpd	%%xmm6,%%xmm6						\n\t"\
			"addpd	%%xmm0,%%xmm2						\n\t"\
			"addpd	%%xmm4,%%xmm7						\n\t"\
			"addpd	%%xmm1,%%xmm3						\n\t"\
			"addpd	%%xmm5,%%xmm6						\n\t"\
			"movaps	%%xmm2,    (%%r15)					\n\t"\
			"movaps	%%xmm7,    (%%rdx)					\n\t"\
			"movaps	%%xmm3,0x10(%%r15)					\n\t"\
			"movaps	%%xmm6,0x10(%%rcx)					\n\t"\
			"\n\t"\
			"/* add,2,1,03 = &a[j1+p12]+p0,1,2,3 */		\n\t"\
			"movslq	%[__p12],%%rsi						\n\t"\
			"shlq	$3,%%rsi							\n\t"\
			"subq	%%rsi,%%rax		/* &a[j1+p12] */	\n\t"\
			"subq	%%rsi,%%r15							\n\t"\
			"subq	%%rsi,%%rcx							\n\t"\
			"subq	%%rsi,%%rdx							\n\t"\
			"/* SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_C(r0g, 0x120, 0x240, edx,ecx,ebx,eax) */\n\t"\
			"subq	$0x220	,%%rdi 	/* edi = __r0g */	\n\t"\
			"movq	$0x120	,%%rsi 						\n\t"\
			"addq	%%rdi	,%%rsi 	/* esi = __r0g + 0x120 */\n\t"\
			"movaps	    (%%rdi),%%xmm4					\n\t"\
			"movaps	    (%%rsi),%%xmm6					\n\t"\
			"movaps	0x10(%%rdi),%%xmm5					\n\t"\
			"movaps	0x10(%%rsi),%%xmm7					\n\t"\
			"addq	$0x240,%%rdi 						\n\t"\
			"addq	$0x240,%%rsi 						\n\t"\
			"movaps	    (%%rdi),%%xmm0					\n\t"\
			"movaps	    (%%rsi),%%xmm2					\n\t"\
			"movaps	0x10(%%rdi),%%xmm1					\n\t"\
			"movaps	0x10(%%rsi),%%xmm3					\n\t"\
			"subpd	%%xmm0,%%xmm4						\n\t"\
			"subpd	%%xmm2,%%xmm6						\n\t"\
			"subpd	%%xmm1,%%xmm5						\n\t"\
			"subpd	%%xmm3,%%xmm7						\n\t"\
			"addpd	%%xmm0,%%xmm0						\n\t"\
			"addpd	%%xmm2,%%xmm2						\n\t"\
			"addpd	%%xmm1,%%xmm1						\n\t"\
			"addpd	%%xmm3,%%xmm3						\n\t"\
			"addpd	%%xmm4,%%xmm0						\n\t"\
			"addpd	%%xmm6,%%xmm2						\n\t"\
			"addpd	%%xmm5,%%xmm1						\n\t"\
			"addpd	%%xmm7,%%xmm3						\n\t"\
			"/* Swap output addresses a <-> d ,b <-> c: */\n\t"\
			"subpd	%%xmm2,%%xmm0						\n\t"\
			"subpd	%%xmm7,%%xmm4						\n\t"\
			"subpd	%%xmm3,%%xmm1						\n\t"\
			"subpd	%%xmm6,%%xmm5						\n\t"\
			"movaps	%%xmm0,    (%%rcx)					\n\t"\
			"movaps	%%xmm4,    (%%r15)					\n\t"\
			"movaps	%%xmm1,0x10(%%rcx)					\n\t"\
			"movaps	%%xmm5,0x10(%%rax)					\n\t"\
			"addpd	%%xmm2,%%xmm2						\n\t"\
			"addpd	%%xmm7,%%xmm7						\n\t"\
			"addpd	%%xmm3,%%xmm3						\n\t"\
			"addpd	%%xmm6,%%xmm6						\n\t"\
			"addpd	%%xmm0,%%xmm2						\n\t"\
			"addpd	%%xmm4,%%xmm7						\n\t"\
			"addpd	%%xmm1,%%xmm3						\n\t"\
			"addpd	%%xmm5,%%xmm6						\n\t"\
			"movaps	%%xmm2,    (%%rdx)					\n\t"\
			"movaps	%%xmm7,    (%%rax)					\n\t"\
			"movaps	%%xmm3,0x10(%%rdx)					\n\t"\
			"movaps	%%xmm6,0x10(%%r15)					\n\t"\
			"\n\t"\
			:					/* outputs: none */\
			: [__add] "m" (Xadd)	/* All inputs from memory addresses here */\
			 ,[__p01] "m" (Xp01)\
			 ,[__p02] "m" (Xp02)\
			 ,[__p03] "m" (Xp03)\
			 ,[__p04] "m" (Xp04)\
			 ,[__p08] "m" (Xp08)\
			 ,[__p12] "m" (Xp12)\
			 ,[__p16] "m" (Xp16)\
			 ,[__p20] "m" (Xp20)\
			 ,[__p24] "m" (Xp24)\
			 ,[__p28] "m" (Xp28)\
			 ,[__p32] "m" (Xp32)\
			 ,[__r00] "m" (Xr00)\
			 ,[__in0] "m" (Xin0)\
			 ,[__cc1] "m" (Xcc1)\
			: "rax","r15","rcx","rdx","rdi","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
		);\
	}

#endif	/* radix36_ditN_cy_dif1_gcc_h_included */

