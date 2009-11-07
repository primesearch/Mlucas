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
#ifndef radix28_ditN_cy_dif1_gcc_h_included
#define radix28_ditN_cy_dif1_gcc_h_included

	#define	SSE2_RADIX28_DIT_NOTWIDDLE(Xadd,Xp01,Xp02,Xp03,Xp04,Xp08,Xp12,Xp16,Xp20,Xp24,Xout,Xcc0)\
	{\
	__asm__ volatile (\
		"/* SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_E(add00+p[0,1,3,2], s1p00r,s1p03r,s1p02r,s1p01r): */\n\t"\
			"movq	%[__out],%%rsi	/* s1p00r */\n\t"\
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
			"movaps	    (%%rax),%%xmm0		/* a(jt   ) */\n\t"\
			"movaps	    (%%rdx),%%xmm4		/* a(jt+p2) */\n\t"\
			"movaps	0x10(%%rax),%%xmm1		/* a(jp   ) */\n\t"\
			"movaps	0x10(%%rdx),%%xmm5		/* a(jp+p2) */\n\t"\
			"movaps	%%xmm0,%%xmm2			/* xmm2 <- cpy a(jt   ) */\n\t"\
			"movaps	%%xmm4,%%xmm6			/* xmm4 <- cpy a(jt+p2) */\n\t"\
			"movaps	%%xmm1,%%xmm3			/* xmm3 <- cpy a(jp   ) */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* xmm5 <- cpy a(jp+p2) */\n\t"\
			"addpd	    (%%r15),%%xmm0		/* t1 */\n\t"\
			"addpd	    (%%rcx),%%xmm4		/* t5 */\n\t"\
			"addpd	0x10(%%r15),%%xmm1		/* t2 */\n\t"\
			"addpd	0x10(%%rcx),%%xmm5		/* t6 */\n\t"\
			"subpd	    (%%r15),%%xmm2		/* t3 */\n\t"\
			"subpd	    (%%rcx),%%xmm6		/* t7 */\n\t"\
			"subpd	0x10(%%r15),%%xmm3		/* t4 */\n\t"\
			"subpd	0x10(%%rcx),%%xmm7		/* t8 */\n\t"\
			"/* Finish radix-4 butterfly and store results into output-array slots: */\n\t"\
			"/* Swap _o1 <-> _o3: 0x20,0x30 <-> 0x60,0x70: */\n\t"\
			"subpd	%%xmm4,%%xmm0			/* ~t5 <- t1 -t5 */\n\t"\
			"subpd	%%xmm7,%%xmm2			/* ~t7 <- t3 -t8 */\n\t"\
			"subpd	%%xmm5,%%xmm1			/* ~t6 <- t2 -t6 */\n\t"\
			"subpd	%%xmm6,%%xmm3			/* ~t4 <- t4 -t7 */\n\t"\
			"movaps	%%xmm0,0x040(%%rsi)		/* <- ~t5 */\n\t"\
			"movaps	%%xmm2,0x020(%%rsi)		/* <- ~t7 */\n\t"\
			"movaps	%%xmm1,0x050(%%rsi)		/* <- ~t6 */\n\t"\
			"movaps	%%xmm3,0x070(%%rsi)		/* <- ~t4 */\n\t"\
			"addpd	%%xmm4,%%xmm4			/*          2*t5 */\n\t"\
			"addpd	%%xmm7,%%xmm7			/*          2*t8 */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*          2*t6 */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*          2*t7 */\n\t"\
			"addpd	%%xmm0,%%xmm4			/* ~t1 <- t1 +t5 */\n\t"\
			"addpd	%%xmm2,%%xmm7			/* ~t3 <- t3 +t8 */\n\t"\
			"addpd	%%xmm1,%%xmm5			/* ~t2 <- t2 +t6 */\n\t"\
			"addpd	%%xmm3,%%xmm6			/* ~t8 <- t4 +t7 */\n\t"\
			"movaps	%%xmm4,     (%%rsi)		/* <- ~t1 */\n\t"\
			"movaps	%%xmm7,0x060(%%rsi)		/* <- ~t3 */\n\t"\
			"movaps	%%xmm5,0x010(%%rsi)		/* <- ~t2 */\n\t"\
			"movaps	%%xmm6,0x030(%%rsi)		/* <- ~t8 */\n\t"\
		"/* SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_E(add12+p[3,2,1,0], s1p04r,s1p07r,s1p06r,s1p05r): */\n\t"\
			"addq	$0x80,%%rsi	/* s1p04r */\n\t"\
			"subq	%%rax,%%r15		/* p01 */\n\t"\
			"subq	%%rax,%%rcx		/* p02 */\n\t"\
			"subq	%%rax,%%rdx		/* p03 */\n\t"\
			"movslq	%[__p12],%%rdi\n\t"\
			"shlq	$3,%%rdi		/* Pointer offset for floating doubles */\n\t"\
			"addq	%%rdi,%%rax		/* add12 */\n\t"\
			"addq	%%rax,%%r15\n\t"\
			"addq	%%rax,%%rcx\n\t"\
			"addq	%%rax,%%rdx\n\t"\
			"/* eax <-> edx, ebx <-> ecx */\n\t"\
			"movaps	    (%%rdx),%%xmm0		/* a(jt   ) */\n\t"\
			"movaps	    (%%r15),%%xmm4		/* a(jt+p2) */\n\t"\
			"movaps	0x10(%%rdx),%%xmm1		/* a(jp   ) */\n\t"\
			"movaps	0x10(%%r15),%%xmm5		/* a(jp+p2) */\n\t"\
			"movaps	%%xmm0,%%xmm2			/* xmm2 <- cpy a(jt   ) */\n\t"\
			"movaps	%%xmm4,%%xmm6			/* xmm4 <- cpy a(jt+p2) */\n\t"\
			"movaps	%%xmm1,%%xmm3			/* xmm3 <- cpy a(jp   ) */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* xmm5 <- cpy a(jp+p2) */\n\t"\
			"addpd	    (%%rcx),%%xmm0		/* t1 */\n\t"\
			"addpd	    (%%rax),%%xmm4		/* t5 */\n\t"\
			"addpd	0x10(%%rcx),%%xmm1		/* t2 */\n\t"\
			"addpd	0x10(%%rax),%%xmm5		/* t6 */\n\t"\
			"subpd	    (%%rcx),%%xmm2		/* t3 */\n\t"\
			"subpd	    (%%rax),%%xmm6		/* t7 */\n\t"\
			"subpd	0x10(%%rcx),%%xmm3		/* t4 */\n\t"\
			"subpd	0x10(%%rax),%%xmm7		/* t8 */\n\t"\
			"/* Finish radix-4 butterfly and store results into output-array slots: */\n\t"\
			"/* Swap _o1 <-> _o3: 0x20,0x30 <-> 0x60,0x70: */\n\t"\
			"subpd	%%xmm4,%%xmm0			/* ~t5 <- t1 -t5 */\n\t"\
			"subpd	%%xmm7,%%xmm2			/* ~t7 <- t3 -t8 */\n\t"\
			"subpd	%%xmm5,%%xmm1			/* ~t6 <- t2 -t6 */\n\t"\
			"subpd	%%xmm6,%%xmm3			/* ~t4 <- t4 -t7 */\n\t"\
			"movaps	%%xmm0,0x040(%%rsi)		/* <- ~t5 */\n\t"\
			"movaps	%%xmm2,0x020(%%rsi)		/* <- ~t7 */\n\t"\
			"movaps	%%xmm1,0x050(%%rsi)		/* <- ~t6 */\n\t"\
			"movaps	%%xmm3,0x070(%%rsi)		/* <- ~t4 */\n\t"\
			"addpd	%%xmm4,%%xmm4			/*          2*t5 */\n\t"\
			"addpd	%%xmm7,%%xmm7			/*          2*t8 */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*          2*t6 */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*          2*t7 */\n\t"\
			"addpd	%%xmm0,%%xmm4			/* ~t1 <- t1 +t5 */\n\t"\
			"addpd	%%xmm2,%%xmm7			/* ~t3 <- t3 +t8 */\n\t"\
			"addpd	%%xmm1,%%xmm5			/* ~t2 <- t2 +t6 */\n\t"\
			"addpd	%%xmm3,%%xmm6			/* ~t8 <- t4 +t7 */\n\t"\
			"movaps	%%xmm4,     (%%rsi)		/* <- ~t1 */\n\t"\
			"movaps	%%xmm7,0x060(%%rsi)		/* <- ~t3 */\n\t"\
			"movaps	%%xmm5,0x010(%%rsi)		/* <- ~t2 */\n\t"\
			"movaps	%%xmm6,0x030(%%rsi)		/* <- ~t8 */\n\t"\
		"/* SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_E(add24+p[1,0,2,3], s1p08r,s1p11r,s1p10r,s1p09r): */\n\t"\
			"addq	$0x80,%%rsi	/* s1p08r */\n\t"\
			"subq	%%rax,%%r15		/* p01 */\n\t"\
			"subq	%%rax,%%rcx		/* p02 */\n\t"\
			"subq	%%rax,%%rdx		/* p03 */\n\t"\
			"subq	%%rdi,%%rax		/* add0 */\n\t"\
			"movslq	%[__p24],%%rdi\n\t"\
			"shlq	$3,%%rdi		/* Pointer offset for floating doubles */\n\t"\
			"addq	%%rdi,%%rax		/* add24 */\n\t"\
			"addq	%%rax,%%r15\n\t"\
			"addq	%%rax,%%rcx\n\t"\
			"addq	%%rax,%%rdx\n\t"\
			"/* eax <-> ebx */\n\t"\
			"movaps	    (%%r15),%%xmm0		/* a(jt   ) */\n\t"\
			"movaps	    (%%rcx),%%xmm4		/* a(jt+p2) */\n\t"\
			"movaps	0x10(%%r15),%%xmm1		/* a(jp   ) */\n\t"\
			"movaps	0x10(%%rcx),%%xmm5		/* a(jp+p2) */\n\t"\
			"movaps	%%xmm0,%%xmm2			/* xmm2 <- cpy a(jt   ) */\n\t"\
			"movaps	%%xmm4,%%xmm6			/* xmm4 <- cpy a(jt+p2) */\n\t"\
			"movaps	%%xmm1,%%xmm3			/* xmm3 <- cpy a(jp   ) */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* xmm5 <- cpy a(jp+p2) */\n\t"\
			"addpd	    (%%rax),%%xmm0		/* t1 */\n\t"\
			"addpd	    (%%rdx),%%xmm4		/* t5 */\n\t"\
			"addpd	0x10(%%rax),%%xmm1		/* t2 */\n\t"\
			"addpd	0x10(%%rdx),%%xmm5		/* t6 */\n\t"\
			"subpd	    (%%rax),%%xmm2		/* t3 */\n\t"\
			"subpd	    (%%rdx),%%xmm6		/* t7 */\n\t"\
			"subpd	0x10(%%rax),%%xmm3		/* t4 */\n\t"\
			"subpd	0x10(%%rdx),%%xmm7		/* t8 */\n\t"\
			"/* Finish radix-4 butterfly and store results into output-array slots: */\n\t"\
			"/* Swap _o1 <-> _o3: 0x20,0x30 <-> 0x60,0x70: */\n\t"\
			"subpd	%%xmm4,%%xmm0			/* ~t5 <- t1 -t5 */\n\t"\
			"subpd	%%xmm7,%%xmm2			/* ~t7 <- t3 -t8 */\n\t"\
			"subpd	%%xmm5,%%xmm1			/* ~t6 <- t2 -t6 */\n\t"\
			"subpd	%%xmm6,%%xmm3			/* ~t4 <- t4 -t7 */\n\t"\
			"movaps	%%xmm0,0x040(%%rsi)		/* <- ~t5 */\n\t"\
			"movaps	%%xmm2,0x020(%%rsi)		/* <- ~t7 */\n\t"\
			"movaps	%%xmm1,0x050(%%rsi)		/* <- ~t6 */\n\t"\
			"movaps	%%xmm3,0x070(%%rsi)		/* <- ~t4 */\n\t"\
			"addpd	%%xmm4,%%xmm4			/*          2*t5 */\n\t"\
			"addpd	%%xmm7,%%xmm7			/*          2*t8 */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*          2*t6 */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*          2*t7 */\n\t"\
			"addpd	%%xmm0,%%xmm4			/* ~t1 <- t1 +t5 */\n\t"\
			"addpd	%%xmm2,%%xmm7			/* ~t3 <- t3 +t8 */\n\t"\
			"addpd	%%xmm1,%%xmm5			/* ~t2 <- t2 +t6 */\n\t"\
			"addpd	%%xmm3,%%xmm6			/* ~t8 <- t4 +t7 */\n\t"\
			"movaps	%%xmm4,     (%%rsi)		/* <- ~t1 */\n\t"\
			"movaps	%%xmm7,0x060(%%rsi)		/* <- ~t3 */\n\t"\
			"movaps	%%xmm5,0x010(%%rsi)		/* <- ~t2 */\n\t"\
			"movaps	%%xmm6,0x030(%%rsi)		/* <- ~t8 */\n\t"\
		"/* SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_E(add08+p[1,0,2,3], s1p12r,s1p15r,s1p14r,s1p13r): */\n\t"\
			"addq	$0x80,%%rsi	/* s1p12r */\n\t"\
			"subq	%%rax,%%r15		/* p01 */\n\t"\
			"subq	%%rax,%%rcx		/* p02 */\n\t"\
			"subq	%%rax,%%rdx		/* p03 */\n\t"\
			"subq	%%rdi,%%rax		/* add0 */\n\t"\
			"movslq	%[__p08],%%rdi\n\t"\
			"shlq	$3,%%rdi		/* Pointer offset for floating doubles */\n\t"\
			"addq	%%rdi,%%rax		/* add08 */\n\t"\
			"addq	%%rax,%%r15\n\t"\
			"addq	%%rax,%%rcx\n\t"\
			"addq	%%rax,%%rdx\n\t"\
			"/* eax <-> %%r15 */\n\t"\
			"movaps	    (%%r15),%%xmm0		/* a(jt   ) */\n\t"\
			"movaps	    (%%rcx),%%xmm4		/* a(jt+p2) */\n\t"\
			"movaps	0x10(%%r15),%%xmm1		/* a(jp   ) */\n\t"\
			"movaps	0x10(%%rcx),%%xmm5		/* a(jp+p2) */\n\t"\
			"movaps	%%xmm0,%%xmm2			/* xmm2 <- cpy a(jt   ) */\n\t"\
			"movaps	%%xmm4,%%xmm6			/* xmm4 <- cpy a(jt+p2) */\n\t"\
			"movaps	%%xmm1,%%xmm3			/* xmm3 <- cpy a(jp   ) */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* xmm5 <- cpy a(jp+p2) */\n\t"\
			"addpd	    (%%rax),%%xmm0		/* t1 */\n\t"\
			"addpd	    (%%rdx),%%xmm4		/* t5 */\n\t"\
			"addpd	0x10(%%rax),%%xmm1		/* t2 */\n\t"\
			"addpd	0x10(%%rdx),%%xmm5		/* t6 */\n\t"\
			"subpd	    (%%rax),%%xmm2		/* t3 */\n\t"\
			"subpd	    (%%rdx),%%xmm6		/* t7 */\n\t"\
			"subpd	0x10(%%rax),%%xmm3		/* t4 */\n\t"\
			"subpd	0x10(%%rdx),%%xmm7		/* t8 */\n\t"\
			"/* Finish radix-4 butterfly and store results into output-array slots: */\n\t"\
			"/* Swap _o1 <-> _o3: 0x20,0x30 <-> 0x60,0x70: */\n\t"\
			"subpd	%%xmm4,%%xmm0			/* ~t5 <- t1 -t5 */\n\t"\
			"subpd	%%xmm7,%%xmm2			/* ~t7 <- t3 -t8 */\n\t"\
			"subpd	%%xmm5,%%xmm1			/* ~t6 <- t2 -t6 */\n\t"\
			"subpd	%%xmm6,%%xmm3			/* ~t4 <- t4 -t7 */\n\t"\
			"movaps	%%xmm0,0x040(%%rsi)		/* <- ~t5 */\n\t"\
			"movaps	%%xmm2,0x020(%%rsi)		/* <- ~t7 */\n\t"\
			"movaps	%%xmm1,0x050(%%rsi)		/* <- ~t6 */\n\t"\
			"movaps	%%xmm3,0x070(%%rsi)		/* <- ~t4 */\n\t"\
			"addpd	%%xmm4,%%xmm4			/*          2*t5 */\n\t"\
			"addpd	%%xmm7,%%xmm7			/*          2*t8 */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*          2*t6 */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*          2*t7 */\n\t"\
			"addpd	%%xmm0,%%xmm4			/* ~t1 <- t1 +t5 */\n\t"\
			"addpd	%%xmm2,%%xmm7			/* ~t3 <- t3 +t8 */\n\t"\
			"addpd	%%xmm1,%%xmm5			/* ~t2 <- t2 +t6 */\n\t"\
			"addpd	%%xmm3,%%xmm6			/* ~t8 <- t4 +t7 */\n\t"\
			"movaps	%%xmm4,     (%%rsi)		/* <- ~t1 */\n\t"\
			"movaps	%%xmm7,0x060(%%rsi)		/* <- ~t3 */\n\t"\
			"movaps	%%xmm5,0x010(%%rsi)		/* <- ~t2 */\n\t"\
			"movaps	%%xmm6,0x030(%%rsi)		/* <- ~t8 */\n\t"\
		"/* SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_E(add20+p[2,3,0,1], s1p16r,s1p19r,s1p18r,s1p17r): */\n\t"\
			"addq	$0x80,%%rsi	/* s1p16r */\n\t"\
			"subq	%%rax,%%r15		/* p01 */\n\t"\
			"subq	%%rax,%%rcx		/* p02 */\n\t"\
			"subq	%%rax,%%rdx		/* p03 */\n\t"\
			"subq	%%rdi,%%rax		/* add0 */\n\t"\
			"movslq	%[__p20],%%rdi\n\t"\
			"shlq	$3,%%rdi		/* Pointer offset for floating doubles */\n\t"\
			"addq	%%rdi,%%rax		/* add20 */\n\t"\
			"addq	%%rax,%%r15\n\t"\
			"addq	%%rax,%%rcx\n\t"\
			"addq	%%rax,%%rdx\n\t"\
			"/* eax <-> %%rcx, %%r15 <-> edx */\n\t"\
			"movaps	    (%%rcx),%%xmm0		/* a(jt   ) */\n\t"\
			"movaps	    (%%rax),%%xmm4		/* a(jt+p2) */\n\t"\
			"movaps	0x10(%%rcx),%%xmm1		/* a(jp   ) */\n\t"\
			"movaps	0x10(%%rax),%%xmm5		/* a(jp+p2) */\n\t"\
			"movaps	%%xmm0,%%xmm2			/* xmm2 <- cpy a(jt   ) */\n\t"\
			"movaps	%%xmm4,%%xmm6			/* xmm4 <- cpy a(jt+p2) */\n\t"\
			"movaps	%%xmm1,%%xmm3			/* xmm3 <- cpy a(jp   ) */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* xmm5 <- cpy a(jp+p2) */\n\t"\
			"addpd	    (%%rdx),%%xmm0		/* t1 */\n\t"\
			"addpd	    (%%r15),%%xmm4		/* t5 */\n\t"\
			"addpd	0x10(%%rdx),%%xmm1		/* t2 */\n\t"\
			"addpd	0x10(%%r15),%%xmm5		/* t6 */\n\t"\
			"subpd	    (%%rdx),%%xmm2		/* t3 */\n\t"\
			"subpd	    (%%r15),%%xmm6		/* t7 */\n\t"\
			"subpd	0x10(%%rdx),%%xmm3		/* t4 */\n\t"\
			"subpd	0x10(%%r15),%%xmm7		/* t8 */\n\t"\
			"/* Finish radix-4 butterfly and store results into output-array slots: */\n\t"\
			"/* Swap _o1 <-> _o3: 0x20,0x30 <-> 0x60,0x70: */\n\t"\
			"subpd	%%xmm4,%%xmm0			/* ~t5 <- t1 -t5 */\n\t"\
			"subpd	%%xmm7,%%xmm2			/* ~t7 <- t3 -t8 */\n\t"\
			"subpd	%%xmm5,%%xmm1			/* ~t6 <- t2 -t6 */\n\t"\
			"subpd	%%xmm6,%%xmm3			/* ~t4 <- t4 -t7 */\n\t"\
			"movaps	%%xmm0,0x040(%%rsi)		/* <- ~t5 */\n\t"\
			"movaps	%%xmm2,0x020(%%rsi)		/* <- ~t7 */\n\t"\
			"movaps	%%xmm1,0x050(%%rsi)		/* <- ~t6 */\n\t"\
			"movaps	%%xmm3,0x070(%%rsi)		/* <- ~t4 */\n\t"\
			"addpd	%%xmm4,%%xmm4			/*          2*t5 */\n\t"\
			"addpd	%%xmm7,%%xmm7			/*          2*t8 */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*          2*t6 */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*          2*t7 */\n\t"\
			"addpd	%%xmm0,%%xmm4			/* ~t1 <- t1 +t5 */\n\t"\
			"addpd	%%xmm2,%%xmm7			/* ~t3 <- t3 +t8 */\n\t"\
			"addpd	%%xmm1,%%xmm5			/* ~t2 <- t2 +t6 */\n\t"\
			"addpd	%%xmm3,%%xmm6			/* ~t8 <- t4 +t7 */\n\t"\
			"movaps	%%xmm4,     (%%rsi)		/* <- ~t1 */\n\t"\
			"movaps	%%xmm7,0x060(%%rsi)		/* <- ~t3 */\n\t"\
			"movaps	%%xmm5,0x010(%%rsi)		/* <- ~t2 */\n\t"\
			"movaps	%%xmm6,0x030(%%rsi)		/* <- ~t8 */\n\t"\
		"/* SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_E(add04+p[2,3,0,1], s1p20r,s1p23r,s1p22r,s1p21r): */\n\t"\
			"addq	$0x80,%%rsi	/* s1p20r */\n\t"\
			"subq	%%rax,%%r15		/* p01 */\n\t"\
			"subq	%%rax,%%rcx		/* p02 */\n\t"\
			"subq	%%rax,%%rdx		/* p03 */\n\t"\
			"subq	%%rdi,%%rax		/* add0 */\n\t"\
			"movslq	%[__p04],%%rdi\n\t"\
			"shlq	$3,%%rdi		/* Pointer offset for floating doubles */\n\t"\
			"addq	%%rdi,%%rax		/* add04 */\n\t"\
			"addq	%%rax,%%r15\n\t"\
			"addq	%%rax,%%rcx\n\t"\
			"addq	%%rax,%%rdx\n\t"\
			"/* eax <-> %%rcx, %%r15 <-> edx */\n\t"\
			"movaps	    (%%rcx),%%xmm0		/* a(jt   ) */\n\t"\
			"movaps	    (%%rax),%%xmm4		/* a(jt+p2) */\n\t"\
			"movaps	0x10(%%rcx),%%xmm1		/* a(jp   ) */\n\t"\
			"movaps	0x10(%%rax),%%xmm5		/* a(jp+p2) */\n\t"\
			"movaps	%%xmm0,%%xmm2			/* xmm2 <- cpy a(jt   ) */\n\t"\
			"movaps	%%xmm4,%%xmm6			/* xmm4 <- cpy a(jt+p2) */\n\t"\
			"movaps	%%xmm1,%%xmm3			/* xmm3 <- cpy a(jp   ) */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* xmm5 <- cpy a(jp+p2) */\n\t"\
			"addpd	    (%%rdx),%%xmm0		/* t1 */\n\t"\
			"addpd	    (%%r15),%%xmm4		/* t5 */\n\t"\
			"addpd	0x10(%%rdx),%%xmm1		/* t2 */\n\t"\
			"addpd	0x10(%%r15),%%xmm5		/* t6 */\n\t"\
			"subpd	    (%%rdx),%%xmm2		/* t3 */\n\t"\
			"subpd	    (%%r15),%%xmm6		/* t7 */\n\t"\
			"subpd	0x10(%%rdx),%%xmm3		/* t4 */\n\t"\
			"subpd	0x10(%%r15),%%xmm7		/* t8 */\n\t"\
			"/* Finish radix-4 butterfly and store results into output-array slots: */\n\t"\
			"/* Swap _o1 <-> _o3: 0x20,0x30 <-> 0x60,0x70: */\n\t"\
			"subpd	%%xmm4,%%xmm0			/* ~t5 <- t1 -t5 */\n\t"\
			"subpd	%%xmm7,%%xmm2			/* ~t7 <- t3 -t8 */\n\t"\
			"subpd	%%xmm5,%%xmm1			/* ~t6 <- t2 -t6 */\n\t"\
			"subpd	%%xmm6,%%xmm3			/* ~t4 <- t4 -t7 */\n\t"\
			"movaps	%%xmm0,0x040(%%rsi)		/* <- ~t5 */\n\t"\
			"movaps	%%xmm2,0x020(%%rsi)		/* <- ~t7 */\n\t"\
			"movaps	%%xmm1,0x050(%%rsi)		/* <- ~t6 */\n\t"\
			"movaps	%%xmm3,0x070(%%rsi)		/* <- ~t4 */\n\t"\
			"addpd	%%xmm4,%%xmm4			/*          2*t5 */\n\t"\
			"addpd	%%xmm7,%%xmm7			/*          2*t8 */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*          2*t6 */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*          2*t7 */\n\t"\
			"addpd	%%xmm0,%%xmm4			/* ~t1 <- t1 +t5 */\n\t"\
			"addpd	%%xmm2,%%xmm7			/* ~t3 <- t3 +t8 */\n\t"\
			"addpd	%%xmm1,%%xmm5			/* ~t2 <- t2 +t6 */\n\t"\
			"addpd	%%xmm3,%%xmm6			/* ~t8 <- t4 +t7 */\n\t"\
			"movaps	%%xmm4,     (%%rsi)		/* <- ~t1 */\n\t"\
			"movaps	%%xmm7,0x060(%%rsi)		/* <- ~t3 */\n\t"\
			"movaps	%%xmm5,0x010(%%rsi)		/* <- ~t2 */\n\t"\
			"movaps	%%xmm6,0x030(%%rsi)		/* <- ~t8 */\n\t"\
		"/* SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_E(add16+p[0,1,3,2], s1p24r,s1p27r,s1p26r,s1p25r): */\n\t"\
			"addq	$0x80,%%rsi	/* s1p24r */\n\t"\
			"subq	%%rax,%%r15		/* p01 */\n\t"\
			"subq	%%rax,%%rcx		/* p02 */\n\t"\
			"subq	%%rax,%%rdx		/* p03 */\n\t"\
			"subq	%%rdi,%%rax		/* add0 */\n\t"\
			"movslq	%[__p16],%%rdi\n\t"\
			"shlq	$3,%%rdi		/* Pointer offset for floating doubles */\n\t"\
			"addq	%%rdi,%%rax		/* add16 */\n\t"\
			"addq	%%rax,%%r15\n\t"\
			"addq	%%rax,%%rcx\n\t"\
			"addq	%%rax,%%rdx\n\t"\
			"/* ecx <-> %%rdx */\n\t"\
			"movaps	    (%%rax),%%xmm0		/* a(jt   ) */\n\t"\
			"movaps	    (%%rdx),%%xmm4		/* a(jt+p2) */\n\t"\
			"movaps	0x10(%%rax),%%xmm1		/* a(jp   ) */\n\t"\
			"movaps	0x10(%%rdx),%%xmm5		/* a(jp+p2) */\n\t"\
			"movaps	%%xmm0,%%xmm2			/* xmm2 <- cpy a(jt   ) */\n\t"\
			"movaps	%%xmm4,%%xmm6			/* xmm4 <- cpy a(jt+p2) */\n\t"\
			"movaps	%%xmm1,%%xmm3			/* xmm3 <- cpy a(jp   ) */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* xmm5 <- cpy a(jp+p2) */\n\t"\
			"addpd	    (%%r15),%%xmm0		/* t1 */\n\t"\
			"addpd	    (%%rcx),%%xmm4		/* t5 */\n\t"\
			"addpd	0x10(%%r15),%%xmm1		/* t2 */\n\t"\
			"addpd	0x10(%%rcx),%%xmm5		/* t6 */\n\t"\
			"subpd	    (%%r15),%%xmm2		/* t3 */\n\t"\
			"subpd	    (%%rcx),%%xmm6		/* t7 */\n\t"\
			"subpd	0x10(%%r15),%%xmm3		/* t4 */\n\t"\
			"subpd	0x10(%%rcx),%%xmm7		/* t8 */\n\t"\
			"/* Finish radix-4 butterfly and store results into output-array slots: */\n\t"\
			"/* Swap _o1 <-> _o3: 0x20,0x30 <-> 0x60,0x70: */\n\t"\
			"subpd	%%xmm4,%%xmm0			/* ~t5 <- t1 -t5 */\n\t"\
			"subpd	%%xmm7,%%xmm2			/* ~t7 <- t3 -t8 */\n\t"\
			"subpd	%%xmm5,%%xmm1			/* ~t6 <- t2 -t6 */\n\t"\
			"subpd	%%xmm6,%%xmm3			/* ~t4 <- t4 -t7 */\n\t"\
			"movaps	%%xmm0,0x040(%%rsi)		/* <- ~t5 */\n\t"\
			"movaps	%%xmm2,0x020(%%rsi)		/* <- ~t7 */\n\t"\
			"movaps	%%xmm1,0x050(%%rsi)		/* <- ~t6 */\n\t"\
			"movaps	%%xmm3,0x070(%%rsi)		/* <- ~t4 */\n\t"\
			"addpd	%%xmm4,%%xmm4			/*          2*t5 */\n\t"\
			"addpd	%%xmm7,%%xmm7			/*          2*t8 */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*          2*t6 */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*          2*t7 */\n\t"\
			"addpd	%%xmm0,%%xmm4			/* ~t1 <- t1 +t5 */\n\t"\
			"addpd	%%xmm2,%%xmm7			/* ~t3 <- t3 +t8 */\n\t"\
			"addpd	%%xmm1,%%xmm5			/* ~t2 <- t2 +t6 */\n\t"\
			"addpd	%%xmm3,%%xmm6			/* ~t8 <- t4 +t7 */\n\t"\
			"movaps	%%xmm4,     (%%rsi)		/* <- ~t1 */\n\t"\
			"movaps	%%xmm7,0x060(%%rsi)		/* <- ~t3 */\n\t"\
			"movaps	%%xmm5,0x010(%%rsi)		/* <- ~t2 */\n\t"\
			"movaps	%%xmm6,0x030(%%rsi)		/* <- ~t8 */\n\t"\
			"\n\t"\
			"/*...and now do 4 radix-7 transforms...*/\n\t"\
			"\n\t"\
		"/* SSE2_RADIX_07_DFT(s1p00r,s1p04r,s1p08r,s1p12r,s1p16r,s1p20r,s1p24r,cc0,s1p00r,s1p08r,s1p16r,s1p24r,s1p04r,s1p12r,s1p20r): */\n\t"\
			"/*\n\t"\
			"t1r=A1r+A6r;	t2r=A2r+A5r;	t3r=A3r+A4r;\n\t"\
			"t6r=A1r-A6r;	t5r=A2r-A5r;	t4r=A3r-A4r;\n\t"\
			"*/\n\t"\
			"movq	%[__out],%%rdi	/* s1p00r */\n\t"\
			"movq	%[__cc0],%%rax	/* cc0 */\n\t"\
			"movq	$0x20,%%r15\n\t"\
			"movq	$0x40,%%rcx\n\t"\
			"movq	$0x60,%%rdx\n\t"\
			"movq	$0x80,%%rsi\n\t"\
			"addq	%%rax,%%r15		/* cc1 */\n\t"\
			"addq	%%rax,%%rcx		/* cc2 */\n\t"\
			"addq	%%rax,%%rdx		/* cc3 */\n\t"\
			"addq	%%rax,%%rsi		/* scratch storage */\n\t"\
			"movaps	0x080(%%rdi),%%xmm6			/* A1r = s1p04r */\n\t"\
			"movaps	0x300(%%rdi),%%xmm1			/* A6r = s1p24r */\n\t"\
			"movaps	0x100(%%rdi),%%xmm5			/* A2r = s1p08r */\n\t"\
			"movaps	0x280(%%rdi),%%xmm2			/* A5r = s1p20r */\n\t"\
			"movaps	0x180(%%rdi),%%xmm4			/* A3r = s1p12r */\n\t"\
			"movaps	0x200(%%rdi),%%xmm3			/* A4r = s1p16r */\n\t"\
			"\n\t"\
			"subpd	%%xmm1,%%xmm6			/* t6r = A1r-A6r */\n\t"\
			"addpd	%%xmm1,%%xmm1			/*         2*A6r */\n\t"\
			"addpd	%%xmm6,%%xmm1			/* t1r = A1r+A6r */\n\t"\
			"\n\t"\
			"subpd	%%xmm2,%%xmm5			/* t5r = A2r-A5r */\n\t"\
			"addpd	%%xmm2,%%xmm2			/*         2*A5r */\n\t"\
			"addpd	%%xmm5,%%xmm2			/* t2r = A2r+A5r */\n\t"\
			"\n\t"\
			"movaps	    (%%rdi),%%xmm0		/* Ar0 */\n\t"\
			"subpd	%%xmm3,%%xmm4			/* t4r = A3r-A4r */\n\t"\
			"addpd	%%xmm3,%%xmm3			/*         2*A4r */\n\t"\
			"addpd	%%xmm4,%%xmm3			/* t3r = A3r+A4r */\n\t"\
			"/*\n\t"\
			"	rt  = t1r+t2r+t3r;			\n\t"\
			"	B0r = rt + A0r;				\n\t"\
			"	t0r = rt*cx0 + A0r;			t3r=(t6r-t4r+t5r)*sx0;\n\t"\
			"	t1r = t1r-t2r;				t6r= t6r-t5r;\n\t"\
			"	t2r = t3r-t2r;				t5r= t4r+t5r;\n\t"\
			"	t3r =(t1r+t2r)*cx3;			t4r=(t5r-t6r)*sx3;\n\t"\
			"	t1r = t1r*cx1;				t6r= t6r*sx1;\n\t"\
			"	t2r = t2r*cx2;				t5r= t5r*sx2;\n\t"\
			"	tt  = t1r-t3r;				t6r= t4r+t6r;\n\t"\
			"	t2r = t2r-t3r;				t5r= t4r-t5r;\n\t"\
			"	\n\t"\
			"	t1r= t0r- tt-t2r;			t4r= t3r-t6r-t5r;\n\t"\
			"	t2r= t0r+t2r;				t5r= t3r+t5r;\n\t"\
			"	t0r= t0r+ tt;				t3r= t3r+t6r;\n\t"\
			"*/\n\t"\
			"movaps	%%xmm0,    (%%rsi)		/* cpy t0 into scratch sincos slot */\n\t"\
			"movaps	%%xmm6,0x10(%%rsi)		/* cpy t6 into scratch sincos slot */\n\t"\
			"addpd	%%xmm1,%%xmm0			/*~A0 = A0+t1 */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* cpy t5 */\n\t"\
			"addpd	%%xmm2,%%xmm3			/*~t3 = t3+t2 */\n\t"\
			"subpd	%%xmm4,%%xmm5			/*~t5 = t5-t4 */\n\t"\
			"subpd	%%xmm2,%%xmm1			/*~t1 = t1-t2 */\n\t"\
			"subpd	%%xmm7,%%xmm6			/*~t6 = t6-t5 */\n\t"\
			"addpd	%%xmm2,%%xmm2			/* 2*t2 */\n\t"\
			"addpd	%%xmm7,%%xmm4			/*~t5 = t4+t5 */\n\t"\
			"addpd	%%xmm3,%%xmm0			/* B0 */\n\t"\
			"addpd	0x10(%%rsi),%%xmm5		/* t3 = (t5-t4)+t6 */\n\t"\
			"subpd	%%xmm2,%%xmm3			/*~t2 =  (t2+t3) - 2*t2 = t3-t2 */\n\t"\
			"movaps	%%xmm4,%%xmm7			/* cpy t5 */\n\t"\
			"movaps	%%xmm0,     (%%rdi)		/* <-B0 = s1p00r, xmm0 FREE */\n\t"\
			"subpd	%%xmm6,%%xmm4			/* t4 = ~t5-~t6 */\n\t"\
			"movaps	%%xmm1,%%xmm2			/* cpy ~t1 */\n\t"\
			"subpd	    (%%rsi),%%xmm0		/* r = B0 - t0 */\n\t"\
			"mulpd	0x10(%%rax),%%xmm5		/*~t3 = t3*sx0 */\n\t"\
			"addpd	%%xmm3,%%xmm2			/* ~t1+~t2 */\n\t"\
			"mulpd	    (%%rcx),%%xmm3		/* t2 = t2*cx2 */\n\t"\
			"mulpd	0x10(%%rdx),%%xmm4		/*~t4 = t4*sx3 */\n\t"\
			"mulpd	    (%%r15),%%xmm1		/* t1 = t1*cx1 */\n\t"\
			"mulpd	0x10(%%r15),%%xmm6		/*~t6 = t6*sx1 */\n\t"\
			"mulpd	    (%%rax),%%xmm0		/* ~r = r*(cx0-1) */\n\t"\
			"mulpd	0x10(%%rcx),%%xmm7		/*~t5 = t5*sx2 */\n\t"\
			"mulpd	    (%%rdx),%%xmm2		/* t3 =(t1+t2)*cx3 */\n\t"\
			"addpd	    (%%rdi),%%xmm0		/* t0 =~r + B0 */\n\t"\
			"addpd	%%xmm4,%%xmm6			/*~t6 = t4+t6 */\n\t"\
			"subpd	%%xmm2,%%xmm1			/* tt = t1-t3 */\n\t"\
			"subpd	%%xmm7,%%xmm4			/*~t5 = t4-t5, xmm7 FREE */\n\t"\
			"subpd	%%xmm2,%%xmm3			/* t2 = t2-t3, xmm2 FREE */\n\t"\
			"\n\t"\
			"movaps	%%xmm0,%%xmm2			/* cpy t0 */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* cpy t3 */\n\t"\
			"addpd	%%xmm1,%%xmm0			/*~t0 = t0+tt */\n\t"\
			"addpd	%%xmm6,%%xmm5			/*~t3 = t3+t6 */\n\t"\
			"addpd	%%xmm3,%%xmm1			/*~tt = tt+t2 */\n\t"\
			"addpd	%%xmm4,%%xmm6			/*      t6+t5 */\n\t"\
			"addpd	%%xmm2,%%xmm3			/*~t2 = t2+t0 */\n\t"\
			"addpd	%%xmm7,%%xmm4			/*~t5 = t5+t3 */\n\t"\
			"subpd	%%xmm1,%%xmm2			/*~t1 = t0-tt-t2 */\n\t"\
			"subpd	%%xmm6,%%xmm7			/*~t4 = t3-t6-t5 */\n\t"\
			"movaps	%%xmm0,0x100(%%rdi)			/* B1 <- t0 = s1p08r */\n\t"\
			"movaps	%%xmm5,0x280(%%rdi)			/* B6 <- t3 = s1p20r */\n\t"\
			"movaps	%%xmm2,0x200(%%rdi)			/* B2 <- t1 = s1p16r */\n\t"\
			"movaps	%%xmm7,0x180(%%rdi)			/* B5 <- t4 = s1p12r */\n\t"\
			"movaps	%%xmm3,0x300(%%rdi)			/* B3 <- t2 = s1p24r */\n\t"\
			"movaps	%%xmm4,0x080(%%rdi)			/* B4 <- t5 = s1p04r */\n\t"\
			"/************************** Imaginary Parts: ******************************************/\n\t"\
			"movaps	0x090(%%rdi),%%xmm6		/* A1i = s1p04i */\n\t"\
			"movaps	0x310(%%rdi),%%xmm1		/* A6i = s1p24i */\n\t"\
			"movaps	0x110(%%rdi),%%xmm5		/* A2i = s1p08i */\n\t"\
			"movaps	0x290(%%rdi),%%xmm2		/* A5i = s1p20i */\n\t"\
			"movaps	0x190(%%rdi),%%xmm4		/* A3i = s1p12i */\n\t"\
			"movaps	0x210(%%rdi),%%xmm3		/* A4i = s1p16i */\n\t"\
			"\n\t"\
			"subpd	%%xmm1,%%xmm6			/* t6i = A1i-A6i */\n\t"\
			"addpd	%%xmm1,%%xmm1			/*         2*A6i */\n\t"\
			"addpd	%%xmm6,%%xmm1			/* t1i = A1i+A6i */\n\t"\
			"\n\t"\
			"subpd	%%xmm2,%%xmm5			/* t5i = A2i-A5i */\n\t"\
			"addpd	%%xmm2,%%xmm2			/*         2*A5i */\n\t"\
			"addpd	%%xmm5,%%xmm2			/* t2i = A2i+A5i */\n\t"\
			"\n\t"\
			"movaps	0x10(%%rdi),%%xmm0		/* Ai0 */\n\t"\
			"subpd	%%xmm3,%%xmm4			/* t4i = A3i-A4i */\n\t"\
			"addpd	%%xmm3,%%xmm3			/*         2*A4i */\n\t"\
			"addpd	%%xmm4,%%xmm3			/* t3i = A3i+A4i */\n\t"\
			"/*\n\t"\
			"	it  = t1i+t2i+t3i;			\n\t"\
			"	B0i = it + A0i;				\n\t"\
			"	t0i = it*cx0 + A0i;			t3i=(t6i-t4i+t5i)*sx0;\n\t"\
			"	t1i = t1i-t2i;				t6i= t6i-t5i;\n\t"\
			"	t2i = t2i-t3i;				t5i= t4i+t5i;\n\t"\
			"	t3i =(t1i-t2i)*cx3;			t4i=(t5i-t6i)*sx3;\n\t"\
			"	t1i = t1i*cx1;				t6i= t6i*sx1;\n\t"\
			"	t2i = t2i*cx2;				t5i= t5i*sx2;\n\t"\
			"	it  = t1i-t3i;				t6i= t4i+t6i;\n\t"\
			"	t2i = t2i-t3i;				t5i= t4i-t5i;\n\t"\
			"	\n\t"\
			"	t1i= t0i- it-t2i;			t4i= t3i-t6i-t5i;\n\t"\
			"	t2i= t0i+t2i;				t5i= t3i+t5i;\n\t"\
			"	t0i= t0i+ it;				t3i= t3i+t6i;\n\t"\
			"*/\n\t"\
			"movaps	%%xmm0,    (%%rsi)		/* cpy t0 into scratch sincos slot */\n\t"\
			"movaps	%%xmm6,0x10(%%rsi)		/* cpy t6 into scratch sincos slot */\n\t"\
			"addpd	%%xmm1,%%xmm0			/*~A0 = A0+t1 */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* cpy t5 */\n\t"\
			"addpd	%%xmm2,%%xmm3			/*~t3 = t3+t2 */\n\t"\
			"subpd	%%xmm4,%%xmm5			/*~t5 = t5-t4 */\n\t"\
			"subpd	%%xmm2,%%xmm1			/*~t1 = t1-t2 */\n\t"\
			"subpd	%%xmm7,%%xmm6			/*~t6 = t6-t5 */\n\t"\
			"addpd	%%xmm2,%%xmm2			/* 2*t2 */\n\t"\
			"addpd	%%xmm7,%%xmm4			/*~t5 = t4+t5 */\n\t"\
			"addpd	%%xmm3,%%xmm0			/* B0 */\n\t"\
			"addpd	0x10(%%rsi),%%xmm5		/* t3 = (t5-t4)+t6 */\n\t"\
			"subpd	%%xmm2,%%xmm3			/*~t2 =  (t2+t3) - 2*t2 = t3-t2 */\n\t"\
			"movaps	%%xmm4,%%xmm7			/* cpy t5 */\n\t"\
			"movaps	%%xmm0,0x010(%%rdi)		/* <-B0 = s1p00i, xmm0 FREE */\n\t"\
			"subpd	%%xmm6,%%xmm4			/* t4 = ~t5-~t6 */\n\t"\
			"movaps	%%xmm1,%%xmm2			/* cpy ~t1 */\n\t"\
			"subpd	    (%%rsi),%%xmm0		/* r = B0 - t0 */\n\t"\
			"mulpd	0x10(%%rax),%%xmm5		/*~t3 = t3*sx0 */\n\t"\
			"addpd	%%xmm3,%%xmm2			/* ~t1+~t2 */\n\t"\
			"mulpd	    (%%rcx),%%xmm3		/* t2 = t2*cx2 */\n\t"\
			"mulpd	0x10(%%rdx),%%xmm4		/*~t4 = t4*sx3 */\n\t"\
			"mulpd	    (%%r15),%%xmm1		/* t1 = t1*cx1 */\n\t"\
			"mulpd	0x10(%%r15),%%xmm6		/*~t6 = t6*sx1 */\n\t"\
			"mulpd	    (%%rax),%%xmm0		/* ~r = r*(cx0-1) */\n\t"\
			"mulpd	0x10(%%rcx),%%xmm7		/*~t5 = t5*sx2 */\n\t"\
			"mulpd	    (%%rdx),%%xmm2		/* t3 =(t1+t2)*cx3 */\n\t"\
			"addpd	0x10(%%rdi),%%xmm0		/* t0 =~r + B0 */\n\t"\
			"addpd	%%xmm4,%%xmm6			/*~t6 = t4+t6 */\n\t"\
			"subpd	%%xmm2,%%xmm1			/* tt = t1-t3 */\n\t"\
			"subpd	%%xmm7,%%xmm4			/*~t5 = t4-t5, xmm7 FREE */\n\t"\
			"subpd	%%xmm2,%%xmm3			/* t2 = t2-t3, xmm2 FREE */\n\t"\
			"movaps	%%xmm0,%%xmm2			/* cpy t0 */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* cpy t3 */\n\t"\
			"addpd	%%xmm1,%%xmm0			/*~t0 = t0+tt */\n\t"\
			"addpd	%%xmm6,%%xmm5			/*~t3 = t3+t6 */\n\t"\
			"addpd	%%xmm3,%%xmm1			/*~tt = tt+t2 */\n\t"\
			"addpd	%%xmm4,%%xmm6			/*      t6+t5 */\n\t"\
			"addpd	%%xmm2,%%xmm3			/*~t2 = t2+t0 */\n\t"\
			"addpd	%%xmm7,%%xmm4			/*~t5 = t5+t3 */\n\t"\
			"subpd	%%xmm1,%%xmm2			/*~t1 = t0-tt-t2, xmm1 FREE */\n\t"\
			"subpd	%%xmm6,%%xmm7			/*~t4 = t3-t6-t5, xmm6 FREE */\n\t"\
			"/*\n\t"\
			"	B1r =t0r-t3i;				B1i*=t0i+t3r;\n\t"\
			"	B2r =t1r-t4i;				B2i*=t1i+t4r;\n\t"\
			"	B3r*=t2r+t5i;				B3i =t2i-t5r;\n\t"\
			"	B4r*=t2r-t5i;				B4i =t2i+t5r;\n\t"\
			"	B5r =t1r+t4i;				B5i*=t1i-t4r;\n\t"\
			"	B6r =t0r+t3i;				B6i*=t0i-t3r;\n\t"\
			"*/\n\t"\
			"/* xmm1,6 FREE */\n\t"\
			"movaps	0x100(%%rdi),%%xmm1		/* t0r */\n\t"\
			"movaps	0x280(%%rdi),%%xmm6		/* t3r */\n\t"\
			"subpd	%%xmm5,%%xmm1			/* B1r =t0r-t3i */\n\t"\
			"subpd	%%xmm6,%%xmm0			/* B6i =t0i-t3r */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*        2*t3i */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*        2*t3r */\n\t"\
			"addpd	%%xmm1,%%xmm5			/* B6r =t0r+t3i */\n\t"\
			"addpd	%%xmm0,%%xmm6			/* B1i =t0i+t3r */\n\t"\
			"movaps	%%xmm1,0x100(%%rdi)		/* <-B1r = s1p08r */\n\t"\
			"movaps	%%xmm0,0x290(%%rdi)		/* <-B6i = s1p20r */\n\t"\
			"movaps	%%xmm5,0x280(%%rdi)		/* <-B6r = s1p20r */\n\t"\
			"movaps	%%xmm6,0x110(%%rdi)		/* <-B1i = s1p08i */\n\t"\
			"\n\t"\
			"movaps	0x200(%%rdi),%%xmm1		/* t1r */\n\t"\
			"movaps	0x180(%%rdi),%%xmm6		/* t4r */\n\t"\
			"subpd	%%xmm7,%%xmm1			/* B2r =t1r-t4i */\n\t"\
			"subpd	%%xmm6,%%xmm2			/* B5i =t1i-t4r */\n\t"\
			"addpd	%%xmm7,%%xmm7			/*        2*t4i */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*        2*t4r */\n\t"\
			"addpd	%%xmm1,%%xmm7			/* B5r =t1r+t4i */\n\t"\
			"addpd	%%xmm2,%%xmm6			/* B2i =t1i+t4r */\n\t"\
			"movaps	%%xmm1,0x200(%%rdi)		/* <-B2r = s1p16r */\n\t"\
			"movaps	%%xmm2,0x190(%%rdi)		/* <-B5i = s1p12r */\n\t"\
			"movaps	%%xmm7,0x180(%%rdi)		/* <-B5r = s1p12r */\n\t"\
			"movaps	%%xmm6,0x210(%%rdi)		/* <-B2i = s1p16i */\n\t"\
			"\n\t"\
			"/* Note the order reversal on this pair of outputs: */\n\t"\
			"movaps	0x300(%%rdi),%%xmm0		/* t2r */\n\t"\
			"movaps	0x080(%%rdi),%%xmm5		/* t5r */\n\t"\
			"subpd	%%xmm4,%%xmm0			/* B4r =t2r-t5i */\n\t"\
			"subpd	%%xmm5,%%xmm3			/* B3i =t2i-t5r */\n\t"\
			"addpd	%%xmm4,%%xmm4			/*        2*t5i */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*        2*t5r */\n\t"\
			"addpd	%%xmm0,%%xmm4			/* B3r =t2r+t5i */\n\t"\
			"addpd	%%xmm3,%%xmm5			/* B4i =t2i+t5r */\n\t"\
			"movaps	%%xmm0,0x080(%%rdi)		/* <-B4r = s1p04r */\n\t"\
			"movaps	%%xmm3,0x310(%%rdi)		/* <-B3i = s1p24r */\n\t"\
			"movaps	%%xmm4,0x300(%%rdi)		/* <-B3r = s1p24r */\n\t"\
			"movaps	%%xmm5,0x090(%%rdi)		/* <-B4i = s1p04i */\n\t"\
			"\n\t"\
		"/* SSE2_RADIX_07_DFT(s1p01r,s1p05r,s1p09r,s1p13r,s1p17r,s1p21r,s1p25r,cc0,s1p21r,s1p01r,s1p09r,s1p17r,s1p25r,s1p05r,s1p13r): */\n\t"\
			"addq	$0x20,%%rdi	/* s1p01r */\n\t"\
			"movaps	0x080(%%rdi),%%xmm6			/* A1r = s1p05r */\n\t"\
			"movaps	0x300(%%rdi),%%xmm1			/* A6r = s1p25r */\n\t"\
			"movaps	0x100(%%rdi),%%xmm5			/* A2r = s1p09r */\n\t"\
			"movaps	0x280(%%rdi),%%xmm2			/* A5r = s1p21r */\n\t"\
			"movaps	0x180(%%rdi),%%xmm4			/* A3r = s1p13r */\n\t"\
			"movaps	0x200(%%rdi),%%xmm3			/* A4r = s1p17r */\n\t"\
			"\n\t"\
			"subpd	%%xmm1,%%xmm6			/* t6r = A1r-A6r */\n\t"\
			"addpd	%%xmm1,%%xmm1			/*         2*A6r */\n\t"\
			"addpd	%%xmm6,%%xmm1			/* t1r = A1r+A6r */\n\t"\
			"\n\t"\
			"subpd	%%xmm2,%%xmm5			/* t5r = A2r-A5r */\n\t"\
			"addpd	%%xmm2,%%xmm2			/*         2*A5r */\n\t"\
			"addpd	%%xmm5,%%xmm2			/* t2r = A2r+A5r */\n\t"\
			"\n\t"\
			"movaps	    (%%rdi),%%xmm0		/* Ar0 */\n\t"\
			"subpd	%%xmm3,%%xmm4			/* t4r = A3r-A4r */\n\t"\
			"addpd	%%xmm3,%%xmm3			/*         2*A4r */\n\t"\
			"addpd	%%xmm4,%%xmm3			/* t3r = A3r+A4r */\n\t"\
			"\n\t"\
			"movaps	%%xmm0,    (%%rsi)		/* cpy t0 into scratch sincos slot */\n\t"\
			"movaps	%%xmm6,0x10(%%rsi)		/* cpy t6 into scratch sincos slot */\n\t"\
			"addpd	%%xmm1,%%xmm0			/*~A0 = A0+t1 */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* cpy t5 */\n\t"\
			"addpd	%%xmm2,%%xmm3			/*~t3 = t3+t2 */\n\t"\
			"subpd	%%xmm4,%%xmm5			/*~t5 = t5-t4 */\n\t"\
			"subpd	%%xmm2,%%xmm1			/*~t1 = t1-t2 */\n\t"\
			"subpd	%%xmm7,%%xmm6			/*~t6 = t6-t5 */\n\t"\
			"addpd	%%xmm2,%%xmm2			/* 2*t2 */\n\t"\
			"addpd	%%xmm7,%%xmm4			/*~t5 = t4+t5 */\n\t"\
			"addpd	%%xmm3,%%xmm0			/* B0 */\n\t"\
			"addpd	0x10(%%rsi),%%xmm5		/* t3 = (t5-t4)+t6 */\n\t"\
			"subpd	%%xmm2,%%xmm3			/*~t2 =  (t2+t3) - 2*t2 = t3-t2 */\n\t"\
			"movaps	%%xmm4,%%xmm7			/* cpy t5 */\n\t"\
			"movaps	%%xmm0,0x280(%%rdi)		/* <-B0 = s1p21r, xmm0 FREE */\n\t"\
			"subpd	%%xmm6,%%xmm4			/* t4 = ~t5-~t6 */\n\t"\
			"movaps	%%xmm1,%%xmm2			/* cpy ~t1 */\n\t"\
			"subpd	    (%%rsi),%%xmm0		/* r = B0 - t0 */\n\t"\
			"mulpd	0x10(%%rax),%%xmm5		/*~t3 = t3*sx0 */\n\t"\
			"addpd	%%xmm3,%%xmm2			/* ~t1+~t2 */\n\t"\
			"mulpd	    (%%rcx),%%xmm3		/* t2 = t2*cx2 */\n\t"\
			"mulpd	0x10(%%rdx),%%xmm4		/*~t4 = t4*sx3 */\n\t"\
			"mulpd	    (%%r15),%%xmm1		/* t1 = t1*cx1 */\n\t"\
			"mulpd	0x10(%%r15),%%xmm6		/*~t6 = t6*sx1 */\n\t"\
			"mulpd	    (%%rax),%%xmm0		/* ~r = r*(cx0-1) */\n\t"\
			"mulpd	0x10(%%rcx),%%xmm7		/*~t5 = t5*sx2 */\n\t"\
			"mulpd	    (%%rdx),%%xmm2		/* t3 =(t1+t2)*cx3 */\n\t"\
			"addpd	0x280(%%rdi),%%xmm0		/* t0 =~r + B0 */\n\t"\
			"addpd	%%xmm4,%%xmm6			/*~t6 = t4+t6 */\n\t"\
			"subpd	%%xmm2,%%xmm1			/* tt = t1-t3 */\n\t"\
			"subpd	%%xmm7,%%xmm4			/*~t5 = t4-t5, xmm7 FREE */\n\t"\
			"subpd	%%xmm2,%%xmm3			/* t2 = t2-t3, xmm2 FREE */\n\t"\
			"\n\t"\
			"movaps	%%xmm0,%%xmm2			/* cpy t0 */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* cpy t3 */\n\t"\
			"addpd	%%xmm1,%%xmm0			/*~t0 = t0+tt */\n\t"\
			"addpd	%%xmm6,%%xmm5			/*~t3 = t3+t6 */\n\t"\
			"addpd	%%xmm3,%%xmm1			/*~tt = tt+t2 */\n\t"\
			"addpd	%%xmm4,%%xmm6			/*      t6+t5 */\n\t"\
			"addpd	%%xmm2,%%xmm3			/*~t2 = t2+t0 */\n\t"\
			"addpd	%%xmm7,%%xmm4			/*~t5 = t5+t3 */\n\t"\
			"subpd	%%xmm1,%%xmm2			/*~t1 = t0-tt-t2 */\n\t"\
			"subpd	%%xmm6,%%xmm7			/*~t4 = t3-t6-t5 */\n\t"\
			"movaps	%%xmm0,     (%%rdi)			/* B1 <- t0 = s1p01r */\n\t"\
			"movaps	%%xmm5,0x180(%%rdi)			/* B6 <- t3 = s1p13r */\n\t"\
			"movaps	%%xmm2,0x100(%%rdi)			/* B2 <- t1 = s1p09r */\n\t"\
			"movaps	%%xmm7,0x080(%%rdi)			/* B5 <- t4 = s1p05r */\n\t"\
			"movaps	%%xmm3,0x200(%%rdi)			/* B3 <- t2 = s1p17r */\n\t"\
			"movaps	%%xmm4,0x300(%%rdi)			/* B4 <- t5 = s1p25r */\n\t"\
			"/************************** Imaginary Parts: ******************************************/\n\t"\
			"movaps	0x090(%%rdi),%%xmm6		/* A1i = s1p05i */\n\t"\
			"movaps	0x310(%%rdi),%%xmm1		/* A6i = s1p25i */\n\t"\
			"movaps	0x110(%%rdi),%%xmm5		/* A2i = s1p09i */\n\t"\
			"movaps	0x290(%%rdi),%%xmm2		/* A5i = s1p21i */\n\t"\
			"movaps	0x190(%%rdi),%%xmm4		/* A3i = s1p13i */\n\t"\
			"movaps	0x210(%%rdi),%%xmm3		/* A4i = s1p17i */\n\t"\
			"\n\t"\
			"subpd	%%xmm1,%%xmm6			/* t6i = A1i-A6i */\n\t"\
			"addpd	%%xmm1,%%xmm1			/*         2*A6i */\n\t"\
			"addpd	%%xmm6,%%xmm1			/* t1i = A1i+A6i */\n\t"\
			"\n\t"\
			"subpd	%%xmm2,%%xmm5			/* t5i = A2i-A5i */\n\t"\
			"addpd	%%xmm2,%%xmm2			/*         2*A5i */\n\t"\
			"addpd	%%xmm5,%%xmm2			/* t2i = A2i+A5i */\n\t"\
			"\n\t"\
			"movaps	0x10(%%rdi),%%xmm0		/* Ai0 */\n\t"\
			"subpd	%%xmm3,%%xmm4			/* t4i = A3i-A4i */\n\t"\
			"addpd	%%xmm3,%%xmm3			/*         2*A4i */\n\t"\
			"addpd	%%xmm4,%%xmm3			/* t3i = A3i+A4i */\n\t"\
			"\n\t"\
			"movaps	%%xmm0,    (%%rsi)		/* cpy t0 into scratch sincos slot */\n\t"\
			"movaps	%%xmm6,0x10(%%rsi)		/* cpy t6 into scratch sincos slot */\n\t"\
			"addpd	%%xmm1,%%xmm0			/*~A0 = A0+t1 */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* cpy t5 */\n\t"\
			"addpd	%%xmm2,%%xmm3			/*~t3 = t3+t2 */\n\t"\
			"subpd	%%xmm4,%%xmm5			/*~t5 = t5-t4 */\n\t"\
			"subpd	%%xmm2,%%xmm1			/*~t1 = t1-t2 */\n\t"\
			"subpd	%%xmm7,%%xmm6			/*~t6 = t6-t5 */\n\t"\
			"addpd	%%xmm2,%%xmm2			/* 2*t2 */\n\t"\
			"addpd	%%xmm7,%%xmm4			/*~t5 = t4+t5 */\n\t"\
			"addpd	%%xmm3,%%xmm0			/* B0 */\n\t"\
			"addpd	0x10(%%rsi),%%xmm5		/* t3 = (t5-t4)+t6 */\n\t"\
			"subpd	%%xmm2,%%xmm3			/*~t2 =  (t2+t3) - 2*t2 = t3-t2 */\n\t"\
			"movaps	%%xmm4,%%xmm7			/* cpy t5 */\n\t"\
			"movaps	%%xmm0,0x290(%%rdi)		/* <-B0 = s1p21i, xmm0 FREE */\n\t"\
			"subpd	%%xmm6,%%xmm4			/* t4 = ~t5-~t6 */\n\t"\
			"movaps	%%xmm1,%%xmm2			/* cpy ~t1 */\n\t"\
			"subpd	    (%%rsi),%%xmm0		/* r = B0 - t0 */\n\t"\
			"mulpd	0x10(%%rax),%%xmm5		/*~t3 = t3*sx0 */\n\t"\
			"addpd	%%xmm3,%%xmm2			/* ~t1+~t2 */\n\t"\
			"mulpd	    (%%rcx),%%xmm3		/* t2 = t2*cx2 */\n\t"\
			"mulpd	0x10(%%rdx),%%xmm4		/*~t4 = t4*sx3 */\n\t"\
			"mulpd	    (%%r15),%%xmm1		/* t1 = t1*cx1 */\n\t"\
			"mulpd	0x10(%%r15),%%xmm6		/*~t6 = t6*sx1 */\n\t"\
			"mulpd	    (%%rax),%%xmm0		/* ~r = r*(cx0-1) */\n\t"\
			"mulpd	0x10(%%rcx),%%xmm7		/*~t5 = t5*sx2 */\n\t"\
			"mulpd	    (%%rdx),%%xmm2		/* t3 =(t1+t2)*cx3 */\n\t"\
			"addpd	0x290(%%rdi),%%xmm0		/* t0 =~r + B0 */\n\t"\
			"addpd	%%xmm4,%%xmm6			/*~t6 = t4+t6 */\n\t"\
			"subpd	%%xmm2,%%xmm1			/* tt = t1-t3 */\n\t"\
			"subpd	%%xmm7,%%xmm4			/*~t5 = t4-t5, xmm7 FREE */\n\t"\
			"subpd	%%xmm2,%%xmm3			/* t2 = t2-t3, xmm2 FREE */\n\t"\
			"movaps	%%xmm0,%%xmm2			/* cpy t0 */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* cpy t3 */\n\t"\
			"addpd	%%xmm1,%%xmm0			/*~t0 = t0+tt */\n\t"\
			"addpd	%%xmm6,%%xmm5			/*~t3 = t3+t6 */\n\t"\
			"addpd	%%xmm3,%%xmm1			/*~tt = tt+t2 */\n\t"\
			"addpd	%%xmm4,%%xmm6			/*      t6+t5 */\n\t"\
			"addpd	%%xmm2,%%xmm3			/*~t2 = t2+t0 */\n\t"\
			"addpd	%%xmm7,%%xmm4			/*~t5 = t5+t3 */\n\t"\
			"subpd	%%xmm1,%%xmm2			/*~t1 = t0-tt-t2, xmm1 FREE */\n\t"\
			"subpd	%%xmm6,%%xmm7			/*~t4 = t3-t6-t5, xmm6 FREE */\n\t"\
			"\n\t"\
			"/* xmm1,6 FREE */\n\t"\
			"movaps	     (%%rdi),%%xmm1		/* t0r */\n\t"\
			"movaps	0x180(%%rdi),%%xmm6		/* t3r */\n\t"\
			"subpd	%%xmm5,%%xmm1			/* B1r =t0r-t3i */\n\t"\
			"subpd	%%xmm6,%%xmm0			/* B6i =t0i-t3r */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*        2*t3i */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*        2*t3r */\n\t"\
			"addpd	%%xmm1,%%xmm5			/* B6r =t0r+t3i */\n\t"\
			"addpd	%%xmm0,%%xmm6			/* B1i =t0i+t3r */\n\t"\
			"movaps	%%xmm1,     (%%rdi)		/* <-B1r = s1p08r */\n\t"\
			"movaps	%%xmm0,0x190(%%rdi)		/* <-B6i = s1p20r */\n\t"\
			"movaps	%%xmm5,0x180(%%rdi)		/* <-B6r = s1p20r */\n\t"\
			"movaps	%%xmm6,0x010(%%rdi)		/* <-B1i = s1p08i */\n\t"\
			"\n\t"\
			"movaps	0x100(%%rdi),%%xmm1		/* t1r */\n\t"\
			"movaps	0x080(%%rdi),%%xmm6		/* t4r */\n\t"\
			"subpd	%%xmm7,%%xmm1			/* B2r =t1r-t4i */\n\t"\
			"subpd	%%xmm6,%%xmm2			/* B5i =t1i-t4r */\n\t"\
			"addpd	%%xmm7,%%xmm7			/*        2*t4i */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*        2*t4r */\n\t"\
			"addpd	%%xmm1,%%xmm7			/* B5r =t1r+t4i */\n\t"\
			"addpd	%%xmm2,%%xmm6			/* B2i =t1i+t4r */\n\t"\
			"movaps	%%xmm1,0x100(%%rdi)		/* <-B2r = s1p16r */\n\t"\
			"movaps	%%xmm2,0x090(%%rdi)		/* <-B5i = s1p12r */\n\t"\
			"movaps	%%xmm7,0x080(%%rdi)		/* <-B5r = s1p12r */\n\t"\
			"movaps	%%xmm6,0x110(%%rdi)		/* <-B2i = s1p16i */\n\t"\
			"\n\t"\
			"/* Note the order reversal on this pair of outputs: */\n\t"\
			"movaps	0x200(%%rdi),%%xmm0		/* t2r */\n\t"\
			"movaps	0x300(%%rdi),%%xmm5		/* t5r */\n\t"\
			"subpd	%%xmm4,%%xmm0			/* B4r =t2r-t5i */\n\t"\
			"subpd	%%xmm5,%%xmm3			/* B3i =t2i-t5r */\n\t"\
			"addpd	%%xmm4,%%xmm4			/*        2*t5i */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*        2*t5r */\n\t"\
			"addpd	%%xmm0,%%xmm4			/* B3r =t2r+t5i */\n\t"\
			"addpd	%%xmm3,%%xmm5			/* B4i =t2i+t5r */\n\t"\
			"movaps	%%xmm0,0x300(%%rdi)		/* <-B4r = s1p04r */\n\t"\
			"movaps	%%xmm3,0x210(%%rdi)		/* <-B3i = s1p24r */\n\t"\
			"movaps	%%xmm4,0x200(%%rdi)		/* <-B3r = s1p24r */\n\t"\
			"movaps	%%xmm5,0x310(%%rdi)		/* <-B4i = s1p04i */\n\t"\
			"\n\t"\
		"/* SSE2_RADIX_07_DFT(s1p02r,s1p06r,s1p10r,s1p14r,s1p18r,s1p22r,s1p26r,cc0,s1p14r,s1p22r,s1p02r,s1p10r,s1p18r,s1p26r,s1p06r): */\n\t"\
			"addq	$0x20,%%rdi	/* s1p02r */\n\t"\
			"movaps	0x080(%%rdi),%%xmm6			/* A1r = s1p06r */\n\t"\
			"movaps	0x300(%%rdi),%%xmm1			/* A6r = s1p26r */\n\t"\
			"movaps	0x100(%%rdi),%%xmm5			/* A2r = s1p10r */\n\t"\
			"movaps	0x280(%%rdi),%%xmm2			/* A5r = s1p22r */\n\t"\
			"movaps	0x180(%%rdi),%%xmm4			/* A3r = s1p14r */\n\t"\
			"movaps	0x200(%%rdi),%%xmm3			/* A4r = s1p18r */\n\t"\
			"\n\t"\
			"subpd	%%xmm1,%%xmm6			/* t6r = A1r-A6r */\n\t"\
			"addpd	%%xmm1,%%xmm1			/*         2*A6r */\n\t"\
			"addpd	%%xmm6,%%xmm1			/* t1r = A1r+A6r */\n\t"\
			"\n\t"\
			"subpd	%%xmm2,%%xmm5			/* t5r = A2r-A5r */\n\t"\
			"addpd	%%xmm2,%%xmm2			/*         2*A5r */\n\t"\
			"addpd	%%xmm5,%%xmm2			/* t2r = A2r+A5r */\n\t"\
			"\n\t"\
			"movaps	    (%%rdi),%%xmm0		/* Ar0 */\n\t"\
			"subpd	%%xmm3,%%xmm4			/* t4r = A3r-A4r */\n\t"\
			"addpd	%%xmm3,%%xmm3			/*         2*A4r */\n\t"\
			"addpd	%%xmm4,%%xmm3			/* t3r = A3r+A4r */\n\t"\
			"\n\t"\
			"movaps	%%xmm0,    (%%rsi)		/* cpy t0 into scratch sincos slot */\n\t"\
			"movaps	%%xmm6,0x10(%%rsi)		/* cpy t6 into scratch sincos slot */\n\t"\
			"addpd	%%xmm1,%%xmm0			/*~A0 = A0+t1 */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* cpy t5 */\n\t"\
			"addpd	%%xmm2,%%xmm3			/*~t3 = t3+t2 */\n\t"\
			"subpd	%%xmm4,%%xmm5			/*~t5 = t5-t4 */\n\t"\
			"subpd	%%xmm2,%%xmm1			/*~t1 = t1-t2 */\n\t"\
			"subpd	%%xmm7,%%xmm6			/*~t6 = t6-t5 */\n\t"\
			"addpd	%%xmm2,%%xmm2			/* 2*t2 */\n\t"\
			"addpd	%%xmm7,%%xmm4			/*~t5 = t4+t5 */\n\t"\
			"addpd	%%xmm3,%%xmm0			/* B0 */\n\t"\
			"addpd	0x10(%%rsi),%%xmm5		/* t3 = (t5-t4)+t6 */\n\t"\
			"subpd	%%xmm2,%%xmm3			/*~t2 =  (t2+t3) - 2*t2 = t3-t2 */\n\t"\
			"movaps	%%xmm4,%%xmm7			/* cpy t5 */\n\t"\
			"movaps	%%xmm0,0x180(%%rdi)		/* <-B0 = s1p14r, xmm0 FREE */\n\t"\
			"subpd	%%xmm6,%%xmm4			/* t4 = ~t5-~t6 */\n\t"\
			"movaps	%%xmm1,%%xmm2			/* cpy ~t1 */\n\t"\
			"subpd	    (%%rsi),%%xmm0		/* r = B0 - t0 */\n\t"\
			"mulpd	0x10(%%rax),%%xmm5		/*~t3 = t3*sx0 */\n\t"\
			"addpd	%%xmm3,%%xmm2			/* ~t1+~t2 */\n\t"\
			"mulpd	    (%%rcx),%%xmm3		/* t2 = t2*cx2 */\n\t"\
			"mulpd	0x10(%%rdx),%%xmm4		/*~t4 = t4*sx3 */\n\t"\
			"mulpd	    (%%r15),%%xmm1		/* t1 = t1*cx1 */\n\t"\
			"mulpd	0x10(%%r15),%%xmm6		/*~t6 = t6*sx1 */\n\t"\
			"mulpd	    (%%rax),%%xmm0		/* ~r = r*(cx0-1) */\n\t"\
			"mulpd	0x10(%%rcx),%%xmm7		/*~t5 = t5*sx2 */\n\t"\
			"mulpd	    (%%rdx),%%xmm2		/* t3 =(t1+t2)*cx3 */\n\t"\
			"addpd	0x180(%%rdi),%%xmm0		/* t0 =~r + B0 */\n\t"\
			"addpd	%%xmm4,%%xmm6			/*~t6 = t4+t6 */\n\t"\
			"subpd	%%xmm2,%%xmm1			/* tt = t1-t3 */\n\t"\
			"subpd	%%xmm7,%%xmm4			/*~t5 = t4-t5, xmm7 FREE */\n\t"\
			"subpd	%%xmm2,%%xmm3			/* t2 = t2-t3, xmm2 FREE */\n\t"\
			"\n\t"\
			"movaps	%%xmm0,%%xmm2			/* cpy t0 */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* cpy t3 */\n\t"\
			"addpd	%%xmm1,%%xmm0			/*~t0 = t0+tt */\n\t"\
			"addpd	%%xmm6,%%xmm5			/*~t3 = t3+t6 */\n\t"\
			"addpd	%%xmm3,%%xmm1			/*~tt = tt+t2 */\n\t"\
			"addpd	%%xmm4,%%xmm6			/*      t6+t5 */\n\t"\
			"addpd	%%xmm2,%%xmm3			/*~t2 = t2+t0 */\n\t"\
			"addpd	%%xmm7,%%xmm4			/*~t5 = t5+t3 */\n\t"\
			"subpd	%%xmm1,%%xmm2			/*~t1 = t0-tt-t2 */\n\t"\
			"subpd	%%xmm6,%%xmm7			/*~t4 = t3-t6-t5 */\n\t"\
			"movaps	%%xmm0,0x280(%%rdi)			/* B1 <- t0 = s1p22r */\n\t"\
			"movaps	%%xmm5,0x080(%%rdi)			/* B6 <- t3 = s1p06r */\n\t"\
			"movaps	%%xmm2,     (%%rdi)			/* B2 <- t1 = s1p02r */\n\t"\
			"movaps	%%xmm7,0x300(%%rdi)			/* B5 <- t4 = s1p26r */\n\t"\
			"movaps	%%xmm3,0x100(%%rdi)			/* B3 <- t2 = s1p10r */\n\t"\
			"movaps	%%xmm4,0x200(%%rdi)			/* B4 <- t5 = s1p18r */\n\t"\
			"/************************** Imaginary Parts: ******************************************/\n\t"\
			"movaps	0x090(%%rdi),%%xmm6		/* A1i = s1p06i */\n\t"\
			"movaps	0x310(%%rdi),%%xmm1		/* A6i = s1p26i */\n\t"\
			"movaps	0x110(%%rdi),%%xmm5		/* A2i = s1p10i */\n\t"\
			"movaps	0x290(%%rdi),%%xmm2		/* A5i = s1p22i */\n\t"\
			"movaps	0x190(%%rdi),%%xmm4		/* A3i = s1p14i */\n\t"\
			"movaps	0x210(%%rdi),%%xmm3		/* A4i = s1p18i */\n\t"\
			"\n\t"\
			"subpd	%%xmm1,%%xmm6			/* t6i = A1i-A6i */\n\t"\
			"addpd	%%xmm1,%%xmm1			/*         2*A6i */\n\t"\
			"addpd	%%xmm6,%%xmm1			/* t1i = A1i+A6i */\n\t"\
			"\n\t"\
			"subpd	%%xmm2,%%xmm5			/* t5i = A2i-A5i */\n\t"\
			"addpd	%%xmm2,%%xmm2			/*         2*A5i */\n\t"\
			"addpd	%%xmm5,%%xmm2			/* t2i = A2i+A5i */\n\t"\
			"\n\t"\
			"movaps	0x10(%%rdi),%%xmm0		/* Ai0 */\n\t"\
			"subpd	%%xmm3,%%xmm4			/* t4i = A3i-A4i */\n\t"\
			"addpd	%%xmm3,%%xmm3			/*         2*A4i */\n\t"\
			"addpd	%%xmm4,%%xmm3			/* t3i = A3i+A4i */\n\t"\
			"\n\t"\
			"movaps	%%xmm0,    (%%rsi)		/* cpy t0 into scratch sincos slot */\n\t"\
			"movaps	%%xmm6,0x10(%%rsi)		/* cpy t6 into scratch sincos slot */\n\t"\
			"addpd	%%xmm1,%%xmm0			/*~A0 = A0+t1 */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* cpy t5 */\n\t"\
			"addpd	%%xmm2,%%xmm3			/*~t3 = t3+t2 */\n\t"\
			"subpd	%%xmm4,%%xmm5			/*~t5 = t5-t4 */\n\t"\
			"subpd	%%xmm2,%%xmm1			/*~t1 = t1-t2 */\n\t"\
			"subpd	%%xmm7,%%xmm6			/*~t6 = t6-t5 */\n\t"\
			"addpd	%%xmm2,%%xmm2			/* 2*t2 */\n\t"\
			"addpd	%%xmm7,%%xmm4			/*~t5 = t4+t5 */\n\t"\
			"addpd	%%xmm3,%%xmm0			/* B0 */\n\t"\
			"addpd	0x10(%%rsi),%%xmm5		/* t3 = (t5-t4)+t6 */\n\t"\
			"subpd	%%xmm2,%%xmm3			/*~t2 =  (t2+t3) - 2*t2 = t3-t2 */\n\t"\
			"movaps	%%xmm4,%%xmm7			/* cpy t5 */\n\t"\
			"movaps	%%xmm0,0x190(%%rdi)		/* <-B0 = s1p14i, xmm0 FREE */\n\t"\
			"subpd	%%xmm6,%%xmm4			/* t4 = ~t5-~t6 */\n\t"\
			"movaps	%%xmm1,%%xmm2			/* cpy ~t1 */\n\t"\
			"subpd	    (%%rsi),%%xmm0		/* r = B0 - t0 */\n\t"\
			"mulpd	0x10(%%rax),%%xmm5		/*~t3 = t3*sx0 */\n\t"\
			"addpd	%%xmm3,%%xmm2			/* ~t1+~t2 */\n\t"\
			"mulpd	    (%%rcx),%%xmm3		/* t2 = t2*cx2 */\n\t"\
			"mulpd	0x10(%%rdx),%%xmm4		/*~t4 = t4*sx3 */\n\t"\
			"mulpd	    (%%r15),%%xmm1		/* t1 = t1*cx1 */\n\t"\
			"mulpd	0x10(%%r15),%%xmm6		/*~t6 = t6*sx1 */\n\t"\
			"mulpd	    (%%rax),%%xmm0		/* ~r = r*(cx0-1) */\n\t"\
			"mulpd	0x10(%%rcx),%%xmm7		/*~t5 = t5*sx2 */\n\t"\
			"mulpd	    (%%rdx),%%xmm2		/* t3 =(t1+t2)*cx3 */\n\t"\
			"addpd	0x190(%%rdi),%%xmm0		/* t0 =~r + B0 */\n\t"\
			"addpd	%%xmm4,%%xmm6			/*~t6 = t4+t6 */\n\t"\
			"subpd	%%xmm2,%%xmm1			/* tt = t1-t3 */\n\t"\
			"subpd	%%xmm7,%%xmm4			/*~t5 = t4-t5, xmm7 FREE */\n\t"\
			"subpd	%%xmm2,%%xmm3			/* t2 = t2-t3, xmm2 FREE */\n\t"\
			"movaps	%%xmm0,%%xmm2			/* cpy t0 */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* cpy t3 */\n\t"\
			"addpd	%%xmm1,%%xmm0			/*~t0 = t0+tt */\n\t"\
			"addpd	%%xmm6,%%xmm5			/*~t3 = t3+t6 */\n\t"\
			"addpd	%%xmm3,%%xmm1			/*~tt = tt+t2 */\n\t"\
			"addpd	%%xmm4,%%xmm6			/*      t6+t5 */\n\t"\
			"addpd	%%xmm2,%%xmm3			/*~t2 = t2+t0 */\n\t"\
			"addpd	%%xmm7,%%xmm4			/*~t5 = t5+t3 */\n\t"\
			"subpd	%%xmm1,%%xmm2			/*~t1 = t0-tt-t2, xmm1 FREE */\n\t"\
			"subpd	%%xmm6,%%xmm7			/*~t4 = t3-t6-t5, xmm6 FREE */\n\t"\
			"\n\t"\
			"/* xmm1,6 FREE */\n\t"\
			"movaps	0x280(%%rdi),%%xmm1		/* t0r */\n\t"\
			"movaps	0x080(%%rdi),%%xmm6		/* t3r */\n\t"\
			"subpd	%%xmm5,%%xmm1			/* B1r =t0r-t3i */\n\t"\
			"subpd	%%xmm6,%%xmm0			/* B6i =t0i-t3r */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*        2*t3i */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*        2*t3r */\n\t"\
			"addpd	%%xmm1,%%xmm5			/* B6r =t0r+t3i */\n\t"\
			"addpd	%%xmm0,%%xmm6			/* B1i =t0i+t3r */\n\t"\
			"movaps	%%xmm1,0x280(%%rdi)		/* <-B1r = s1p08r */\n\t"\
			"movaps	%%xmm0,0x090(%%rdi)		/* <-B6i = s1p20r */\n\t"\
			"movaps	%%xmm5,0x080(%%rdi)		/* <-B6r = s1p20r */\n\t"\
			"movaps	%%xmm6,0x290(%%rdi)		/* <-B1i = s1p08i */\n\t"\
			"\n\t"\
			"movaps	     (%%rdi),%%xmm1		/* t1r */\n\t"\
			"movaps	0x300(%%rdi),%%xmm6		/* t4r */\n\t"\
			"subpd	%%xmm7,%%xmm1			/* B2r =t1r-t4i */\n\t"\
			"subpd	%%xmm6,%%xmm2			/* B5i =t1i-t4r */\n\t"\
			"addpd	%%xmm7,%%xmm7			/*        2*t4i */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*        2*t4r */\n\t"\
			"addpd	%%xmm1,%%xmm7			/* B5r =t1r+t4i */\n\t"\
			"addpd	%%xmm2,%%xmm6			/* B2i =t1i+t4r */\n\t"\
			"movaps	%%xmm1,     (%%rdi)		/* <-B2r = s1p16r */\n\t"\
			"movaps	%%xmm2,0x310(%%rdi)		/* <-B5i = s1p12r */\n\t"\
			"movaps	%%xmm7,0x300(%%rdi)		/* <-B5r = s1p12r */\n\t"\
			"movaps	%%xmm6,0x010(%%rdi)		/* <-B2i = s1p16i */\n\t"\
			"\n\t"\
			"/* Note the order reversal on this pair of outputs: */\n\t"\
			"movaps	0x100(%%rdi),%%xmm0		/* t2r */\n\t"\
			"movaps	0x200(%%rdi),%%xmm5		/* t5r */\n\t"\
			"subpd	%%xmm4,%%xmm0			/* B4r =t2r-t5i */\n\t"\
			"subpd	%%xmm5,%%xmm3			/* B3i =t2i-t5r */\n\t"\
			"addpd	%%xmm4,%%xmm4			/*        2*t5i */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*        2*t5r */\n\t"\
			"addpd	%%xmm0,%%xmm4			/* B3r =t2r+t5i */\n\t"\
			"addpd	%%xmm3,%%xmm5			/* B4i =t2i+t5r */\n\t"\
			"movaps	%%xmm0,0x200(%%rdi)		/* <-B4r = s1p04r */\n\t"\
			"movaps	%%xmm3,0x110(%%rdi)		/* <-B3i = s1p24r */\n\t"\
			"movaps	%%xmm4,0x100(%%rdi)		/* <-B3r = s1p24r */\n\t"\
			"movaps	%%xmm5,0x210(%%rdi)		/* <-B4i = s1p04i */\n\t"\
			"\n\t"\
		"/* SSE2_RADIX_07_DFT(s1p03r,s1p07r,s1p11r,s1p15r,s1p19r,s1p23r,s1p27r,cc0,s1p07r,s1p15r,s1p23r,s1p03r,s1p11r,s1p19r,s1p27r): */\n\t"\
			"addq	$0x20,%%rdi	/* s1p03r */\n\t"\
			"movaps	0x080(%%rdi),%%xmm6			/* A1r = s1p07r */\n\t"\
			"movaps	0x300(%%rdi),%%xmm1			/* A6r = s1p27r */\n\t"\
			"movaps	0x100(%%rdi),%%xmm5			/* A2r = s1p11r */\n\t"\
			"movaps	0x280(%%rdi),%%xmm2			/* A5r = s1p23r */\n\t"\
			"movaps	0x180(%%rdi),%%xmm4			/* A3r = s1p15r */\n\t"\
			"movaps	0x200(%%rdi),%%xmm3			/* A4r = s1p19r */\n\t"\
			"\n\t"\
			"subpd	%%xmm1,%%xmm6			/* t6r = A1r-A6r */\n\t"\
			"addpd	%%xmm1,%%xmm1			/*         2*A6r */\n\t"\
			"addpd	%%xmm6,%%xmm1			/* t1r = A1r+A6r */\n\t"\
			"\n\t"\
			"subpd	%%xmm2,%%xmm5			/* t5r = A2r-A5r */\n\t"\
			"addpd	%%xmm2,%%xmm2			/*         2*A5r */\n\t"\
			"addpd	%%xmm5,%%xmm2			/* t2r = A2r+A5r */\n\t"\
			"\n\t"\
			"movaps	    (%%rdi),%%xmm0		/* Ar0 */\n\t"\
			"subpd	%%xmm3,%%xmm4			/* t4r = A3r-A4r */\n\t"\
			"addpd	%%xmm3,%%xmm3			/*         2*A4r */\n\t"\
			"addpd	%%xmm4,%%xmm3			/* t3r = A3r+A4r */\n\t"\
			"\n\t"\
			"movaps	%%xmm0,    (%%rsi)		/* cpy t0 into scratch sincos slot */\n\t"\
			"movaps	%%xmm6,0x10(%%rsi)		/* cpy t6 into scratch sincos slot */\n\t"\
			"addpd	%%xmm1,%%xmm0			/*~A0 = A0+t1 */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* cpy t5 */\n\t"\
			"addpd	%%xmm2,%%xmm3			/*~t3 = t3+t2 */\n\t"\
			"subpd	%%xmm4,%%xmm5			/*~t5 = t5-t4 */\n\t"\
			"subpd	%%xmm2,%%xmm1			/*~t1 = t1-t2 */\n\t"\
			"subpd	%%xmm7,%%xmm6			/*~t6 = t6-t5 */\n\t"\
			"addpd	%%xmm2,%%xmm2			/* 2*t2 */\n\t"\
			"addpd	%%xmm7,%%xmm4			/*~t5 = t4+t5 */\n\t"\
			"addpd	%%xmm3,%%xmm0			/* B0 */\n\t"\
			"addpd	0x10(%%rsi),%%xmm5		/* t3 = (t5-t4)+t6 */\n\t"\
			"subpd	%%xmm2,%%xmm3			/*~t2 =  (t2+t3) - 2*t2 = t3-t2 */\n\t"\
			"movaps	%%xmm4,%%xmm7			/* cpy t5 */\n\t"\
			"movaps	%%xmm0,0x080(%%rdi)		/* <-B0 = s1p07r, xmm0 FREE */\n\t"\
			"subpd	%%xmm6,%%xmm4			/* t4 = ~t5-~t6 */\n\t"\
			"movaps	%%xmm1,%%xmm2			/* cpy ~t1 */\n\t"\
			"subpd	    (%%rsi),%%xmm0		/* r = B0 - t0 */\n\t"\
			"mulpd	0x10(%%rax),%%xmm5		/*~t3 = t3*sx0 */\n\t"\
			"addpd	%%xmm3,%%xmm2			/* ~t1+~t2 */\n\t"\
			"mulpd	    (%%rcx),%%xmm3		/* t2 = t2*cx2 */\n\t"\
			"mulpd	0x10(%%rdx),%%xmm4		/*~t4 = t4*sx3 */\n\t"\
			"mulpd	    (%%r15),%%xmm1		/* t1 = t1*cx1 */\n\t"\
			"mulpd	0x10(%%r15),%%xmm6		/*~t6 = t6*sx1 */\n\t"\
			"mulpd	    (%%rax),%%xmm0		/* ~r = r*(cx0-1) */\n\t"\
			"mulpd	0x10(%%rcx),%%xmm7		/*~t5 = t5*sx2 */\n\t"\
			"mulpd	    (%%rdx),%%xmm2		/* t3 =(t1+t2)*cx3 */\n\t"\
			"addpd	0x080(%%rdi),%%xmm0		/* t0 =~r + B0 */\n\t"\
			"addpd	%%xmm4,%%xmm6			/*~t6 = t4+t6 */\n\t"\
			"subpd	%%xmm2,%%xmm1			/* tt = t1-t3 */\n\t"\
			"subpd	%%xmm7,%%xmm4			/*~t5 = t4-t5, xmm7 FREE */\n\t"\
			"subpd	%%xmm2,%%xmm3			/* t2 = t2-t3, xmm2 FREE */\n\t"\
			"\n\t"\
			"movaps	%%xmm0,%%xmm2			/* cpy t0 */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* cpy t3 */\n\t"\
			"addpd	%%xmm1,%%xmm0			/*~t0 = t0+tt */\n\t"\
			"addpd	%%xmm6,%%xmm5			/*~t3 = t3+t6 */\n\t"\
			"addpd	%%xmm3,%%xmm1			/*~tt = tt+t2 */\n\t"\
			"addpd	%%xmm4,%%xmm6			/*      t6+t5 */\n\t"\
			"addpd	%%xmm2,%%xmm3			/*~t2 = t2+t0 */\n\t"\
			"addpd	%%xmm7,%%xmm4			/*~t5 = t5+t3 */\n\t"\
			"subpd	%%xmm1,%%xmm2			/*~t1 = t0-tt-t2 */\n\t"\
			"subpd	%%xmm6,%%xmm7			/*~t4 = t3-t6-t5 */\n\t"\
			"movaps	%%xmm0,0x180(%%rdi)			/* B1 <- t0 = s1p15r */\n\t"\
			"movaps	%%xmm5,0x300(%%rdi)			/* B6 <- t3 = s1p27r */\n\t"\
			"movaps	%%xmm2,0x280(%%rdi)			/* B2 <- t1 = s1p23r */\n\t"\
			"movaps	%%xmm7,0x200(%%rdi)			/* B5 <- t4 = s1p19r */\n\t"\
			"movaps	%%xmm3,     (%%rdi)			/* B3 <- t2 = s1p03r */\n\t"\
			"movaps	%%xmm4,0x100(%%rdi)			/* B4 <- t5 = s1p11r */\n\t"\
			"/************************** Imaginary Parts: ******************************************/\n\t"\
			"movaps	0x090(%%rdi),%%xmm6		/* A1i = s1p07i */\n\t"\
			"movaps	0x310(%%rdi),%%xmm1		/* A6i = s1p27i */\n\t"\
			"movaps	0x110(%%rdi),%%xmm5		/* A2i = s1p11i */\n\t"\
			"movaps	0x290(%%rdi),%%xmm2		/* A5i = s1p23i */\n\t"\
			"movaps	0x190(%%rdi),%%xmm4		/* A3i = s1p15i */\n\t"\
			"movaps	0x210(%%rdi),%%xmm3		/* A4i = s1p19i */\n\t"\
			"\n\t"\
			"subpd	%%xmm1,%%xmm6			/* t6i = A1i-A6i */\n\t"\
			"addpd	%%xmm1,%%xmm1			/*         2*A6i */\n\t"\
			"addpd	%%xmm6,%%xmm1			/* t1i = A1i+A6i */\n\t"\
			"\n\t"\
			"subpd	%%xmm2,%%xmm5			/* t5i = A2i-A5i */\n\t"\
			"addpd	%%xmm2,%%xmm2			/*         2*A5i */\n\t"\
			"addpd	%%xmm5,%%xmm2			/* t2i = A2i+A5i */\n\t"\
			"\n\t"\
			"movaps	0x10(%%rdi),%%xmm0		/* Ai0 */\n\t"\
			"subpd	%%xmm3,%%xmm4			/* t4i = A3i-A4i */\n\t"\
			"addpd	%%xmm3,%%xmm3			/*         2*A4i */\n\t"\
			"addpd	%%xmm4,%%xmm3			/* t3i = A3i+A4i */\n\t"\
			"\n\t"\
			"movaps	%%xmm0,    (%%rsi)		/* cpy t0 into scratch sincos slot */\n\t"\
			"movaps	%%xmm6,0x10(%%rsi)		/* cpy t6 into scratch sincos slot */\n\t"\
			"addpd	%%xmm1,%%xmm0			/*~A0 = A0+t1 */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* cpy t5 */\n\t"\
			"addpd	%%xmm2,%%xmm3			/*~t3 = t3+t2 */\n\t"\
			"subpd	%%xmm4,%%xmm5			/*~t5 = t5-t4 */\n\t"\
			"subpd	%%xmm2,%%xmm1			/*~t1 = t1-t2 */\n\t"\
			"subpd	%%xmm7,%%xmm6			/*~t6 = t6-t5 */\n\t"\
			"addpd	%%xmm2,%%xmm2			/* 2*t2 */\n\t"\
			"addpd	%%xmm7,%%xmm4			/*~t5 = t4+t5 */\n\t"\
			"addpd	%%xmm3,%%xmm0			/* B0 */\n\t"\
			"addpd	0x10(%%rsi),%%xmm5		/* t3 = (t5-t4)+t6 */\n\t"\
			"subpd	%%xmm2,%%xmm3			/*~t2 =  (t2+t3) - 2*t2 = t3-t2 */\n\t"\
			"movaps	%%xmm4,%%xmm7			/* cpy t5 */\n\t"\
			"movaps	%%xmm0,0x090(%%rdi)		/* <-B0 = s1p07i, xmm0 FREE */\n\t"\
			"subpd	%%xmm6,%%xmm4			/* t4 = ~t5-~t6 */\n\t"\
			"movaps	%%xmm1,%%xmm2			/* cpy ~t1 */\n\t"\
			"subpd	    (%%rsi),%%xmm0		/* r = B0 - t0 */\n\t"\
			"mulpd	0x10(%%rax),%%xmm5		/*~t3 = t3*sx0 */\n\t"\
			"addpd	%%xmm3,%%xmm2			/* ~t1+~t2 */\n\t"\
			"mulpd	    (%%rcx),%%xmm3		/* t2 = t2*cx2 */\n\t"\
			"mulpd	0x10(%%rdx),%%xmm4		/*~t4 = t4*sx3 */\n\t"\
			"mulpd	    (%%r15),%%xmm1		/* t1 = t1*cx1 */\n\t"\
			"mulpd	0x10(%%r15),%%xmm6		/*~t6 = t6*sx1 */\n\t"\
			"mulpd	    (%%rax),%%xmm0		/* ~r = r*(cx0-1) */\n\t"\
			"mulpd	0x10(%%rcx),%%xmm7		/*~t5 = t5*sx2 */\n\t"\
			"mulpd	    (%%rdx),%%xmm2		/* t3 =(t1+t2)*cx3 */\n\t"\
			"addpd	0x090(%%rdi),%%xmm0		/* t0 =~r + B0 */\n\t"\
			"addpd	%%xmm4,%%xmm6			/*~t6 = t4+t6 */\n\t"\
			"subpd	%%xmm2,%%xmm1			/* tt = t1-t3 */\n\t"\
			"subpd	%%xmm7,%%xmm4			/*~t5 = t4-t5, xmm7 FREE */\n\t"\
			"subpd	%%xmm2,%%xmm3			/* t2 = t2-t3, xmm2 FREE */\n\t"\
			"movaps	%%xmm0,%%xmm2			/* cpy t0 */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* cpy t3 */\n\t"\
			"addpd	%%xmm1,%%xmm0			/*~t0 = t0+tt */\n\t"\
			"addpd	%%xmm6,%%xmm5			/*~t3 = t3+t6 */\n\t"\
			"addpd	%%xmm3,%%xmm1			/*~tt = tt+t2 */\n\t"\
			"addpd	%%xmm4,%%xmm6			/*      t6+t5 */\n\t"\
			"addpd	%%xmm2,%%xmm3			/*~t2 = t2+t0 */\n\t"\
			"addpd	%%xmm7,%%xmm4			/*~t5 = t5+t3 */\n\t"\
			"subpd	%%xmm1,%%xmm2			/*~t1 = t0-tt-t2, xmm1 FREE */\n\t"\
			"subpd	%%xmm6,%%xmm7			/*~t4 = t3-t6-t5, xmm6 FREE */\n\t"\
			"\n\t"\
			"/* xmm1,6 FREE */\n\t"\
			"movaps	0x180(%%rdi),%%xmm1		/* t0r */\n\t"\
			"movaps	0x300(%%rdi),%%xmm6		/* t3r */\n\t"\
			"subpd	%%xmm5,%%xmm1			/* B1r =t0r-t3i */\n\t"\
			"subpd	%%xmm6,%%xmm0			/* B6i =t0i-t3r */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*        2*t3i */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*        2*t3r */\n\t"\
			"addpd	%%xmm1,%%xmm5			/* B6r =t0r+t3i */\n\t"\
			"addpd	%%xmm0,%%xmm6			/* B1i =t0i+t3r */\n\t"\
			"movaps	%%xmm1,0x180(%%rdi)		/* <-B1r = s1p08r */\n\t"\
			"movaps	%%xmm0,0x310(%%rdi)		/* <-B6i = s1p20r */\n\t"\
			"movaps	%%xmm5,0x300(%%rdi)		/* <-B6r = s1p20r */\n\t"\
			"movaps	%%xmm6,0x190(%%rdi)		/* <-B1i = s1p08i */\n\t"\
			"\n\t"\
			"movaps	0x280(%%rdi),%%xmm1		/* t1r */\n\t"\
			"movaps	0x200(%%rdi),%%xmm6		/* t4r */\n\t"\
			"subpd	%%xmm7,%%xmm1			/* B2r =t1r-t4i */\n\t"\
			"subpd	%%xmm6,%%xmm2			/* B5i =t1i-t4r */\n\t"\
			"addpd	%%xmm7,%%xmm7			/*        2*t4i */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*        2*t4r */\n\t"\
			"addpd	%%xmm1,%%xmm7			/* B5r =t1r+t4i */\n\t"\
			"addpd	%%xmm2,%%xmm6			/* B2i =t1i+t4r */\n\t"\
			"movaps	%%xmm1,0x280(%%rdi)		/* <-B2r = s1p16r */\n\t"\
			"movaps	%%xmm2,0x210(%%rdi)		/* <-B5i = s1p12r */\n\t"\
			"movaps	%%xmm7,0x200(%%rdi)		/* <-B5r = s1p12r */\n\t"\
			"movaps	%%xmm6,0x290(%%rdi)		/* <-B2i = s1p16i */\n\t"\
			"\n\t"\
			"/* Note the order reversal on this pair of outputs: */\n\t"\
			"movaps	     (%%rdi),%%xmm0		/* t2r */\n\t"\
			"movaps	0x100(%%rdi),%%xmm5		/* t5r */\n\t"\
			"subpd	%%xmm4,%%xmm0			/* B4r =t2r-t5i */\n\t"\
			"subpd	%%xmm5,%%xmm3			/* B3i =t2i-t5r */\n\t"\
			"addpd	%%xmm4,%%xmm4			/*        2*t5i */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*        2*t5r */\n\t"\
			"addpd	%%xmm0,%%xmm4			/* B3r =t2r+t5i */\n\t"\
			"addpd	%%xmm3,%%xmm5			/* B4i =t2i+t5r */\n\t"\
			"movaps	%%xmm0,0x100(%%rdi)		/* <-B4r = s1p04r */\n\t"\
			"movaps	%%xmm3,0x010(%%rdi)		/* <-B3i = s1p24r */\n\t"\
			"movaps	%%xmm4,     (%%rdi)		/* <-B3r = s1p24r */\n\t"\
			"movaps	%%xmm5,0x110(%%rdi)		/* <-B4i = s1p04i */\n\t"\
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
			 ,[__out] "m" (Xout)\
			 ,[__cc0] "m" (Xcc0)\
			: "rax","r15","rcx","rdx","rdi","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
		);\
	}


	#define	SSE2_RADIX28_DIF_NOTWIDDLE(Xadd,Xp01,Xp02,Xp03,Xp04,Xp08,Xp12,Xp16,Xp20,Xp24,Xout,Xcc0)\
	{\
	__asm__ volatile (\
		"/* SSE2_RADIX_07_DFT(s1p00r,s1p24r,s1p20r,s1p16r,s1p12r,s1p08r,s1p04r,cc0,s1p00r,s1p04r,s1p08r,s1p12r,s1p16r,s1p20r,s1p24r): */\n\t"\
			"movq	%[__out],%%rdi	/* s1p00r */\n\t"\
			"movq	%[__cc0],%%rax	/* cc0 */\n\t"\
			"movq	$0x20,%%r15\n\t"\
			"movq	$0x40,%%rcx\n\t"\
			"movq	$0x60,%%rdx\n\t"\
			"movq	$0x80,%%rsi\n\t"\
			"addq	%%rax,%%r15		/* cc1 */\n\t"\
			"addq	%%rax,%%rcx		/* cc2 */\n\t"\
			"addq	%%rax,%%rdx		/* cc3 */\n\t"\
			"addq	%%rax,%%rsi		/* scratch storage */\n\t"\
			"movaps	0x300(%%rdi),%%xmm6			/* A1r = s1p24r */\n\t"\
			"movaps	0x080(%%rdi),%%xmm1			/* A6r = s1p04r */\n\t"\
			"movaps	0x280(%%rdi),%%xmm5			/* A2r = s1p20r */\n\t"\
			"movaps	0x100(%%rdi),%%xmm2			/* A5r = s1p08r */\n\t"\
			"movaps	0x200(%%rdi),%%xmm4			/* A3r = s1p16r */\n\t"\
			"movaps	0x180(%%rdi),%%xmm3			/* A4r = s1p12r */\n\t"\
			"\n\t"\
			"subpd	%%xmm1,%%xmm6			/* t6r = A1r-A6r */\n\t"\
			"addpd	%%xmm1,%%xmm1			/*         2*A6r */\n\t"\
			"addpd	%%xmm6,%%xmm1			/* t1r = A1r+A6r */\n\t"\
			"\n\t"\
			"subpd	%%xmm2,%%xmm5			/* t5r = A2r-A5r */\n\t"\
			"addpd	%%xmm2,%%xmm2			/*         2*A5r */\n\t"\
			"addpd	%%xmm5,%%xmm2			/* t2r = A2r+A5r */\n\t"\
			"\n\t"\
			"movaps	    (%%rdi),%%xmm0		/* Ar0 */\n\t"\
			"subpd	%%xmm3,%%xmm4			/* t4r = A3r-A4r */\n\t"\
			"addpd	%%xmm3,%%xmm3			/*         2*A4r */\n\t"\
			"addpd	%%xmm4,%%xmm3			/* t3r = A3r+A4r */\n\t"\
			"\n\t"\
			"movaps	%%xmm0,    (%%rsi)		/* cpy t0 into scratch sincos slot */\n\t"\
			"movaps	%%xmm6,0x10(%%rsi)		/* cpy t6 into scratch sincos slot */\n\t"\
			"addpd	%%xmm1,%%xmm0			/*~A0 = A0+t1 */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* cpy t5 */\n\t"\
			"addpd	%%xmm2,%%xmm3			/*~t3 = t3+t2 */\n\t"\
			"subpd	%%xmm4,%%xmm5			/*~t5 = t5-t4 */\n\t"\
			"subpd	%%xmm2,%%xmm1			/*~t1 = t1-t2 */\n\t"\
			"subpd	%%xmm7,%%xmm6			/*~t6 = t6-t5 */\n\t"\
			"addpd	%%xmm2,%%xmm2			/* 2*t2 */\n\t"\
			"addpd	%%xmm7,%%xmm4			/*~t5 = t4+t5 */\n\t"\
			"addpd	%%xmm3,%%xmm0			/* B0 */\n\t"\
			"addpd	0x10(%%rsi),%%xmm5		/* t3 = (t5-t4)+t6 */\n\t"\
			"subpd	%%xmm2,%%xmm3			/*~t2 =  (t2+t3) - 2*t2 = t3-t2 */\n\t"\
			"movaps	%%xmm4,%%xmm7			/* cpy t5 */\n\t"\
			"movaps	%%xmm0,     (%%rdi)		/* <-B0 = s1p00r, xmm0 FREE */\n\t"\
			"subpd	%%xmm6,%%xmm4			/* t4 = ~t5-~t6 */\n\t"\
			"movaps	%%xmm1,%%xmm2			/* cpy ~t1 */\n\t"\
			"subpd	    (%%rsi),%%xmm0		/* r = B0 - t0 */\n\t"\
			"mulpd	0x10(%%rax),%%xmm5		/*~t3 = t3*sx0 */\n\t"\
			"addpd	%%xmm3,%%xmm2			/* ~t1+~t2 */\n\t"\
			"mulpd	    (%%rcx),%%xmm3		/* t2 = t2*cx2 */\n\t"\
			"mulpd	0x10(%%rdx),%%xmm4		/*~t4 = t4*sx3 */\n\t"\
			"mulpd	    (%%r15),%%xmm1		/* t1 = t1*cx1 */\n\t"\
			"mulpd	0x10(%%r15),%%xmm6		/*~t6 = t6*sx1 */\n\t"\
			"mulpd	    (%%rax),%%xmm0		/* ~r = r*(cx0-1) */\n\t"\
			"mulpd	0x10(%%rcx),%%xmm7		/*~t5 = t5*sx2 */\n\t"\
			"mulpd	    (%%rdx),%%xmm2		/* t3 =(t1+t2)*cx3 */\n\t"\
			"addpd	    (%%rdi),%%xmm0		/* t0 =~r + B0 */\n\t"\
			"addpd	%%xmm4,%%xmm6			/*~t6 = t4+t6 */\n\t"\
			"subpd	%%xmm2,%%xmm1			/* tt = t1-t3 */\n\t"\
			"subpd	%%xmm7,%%xmm4			/*~t5 = t4-t5, xmm7 FREE */\n\t"\
			"subpd	%%xmm2,%%xmm3			/* t2 = t2-t3, xmm2 FREE */\n\t"\
			"\n\t"\
			"movaps	%%xmm0,%%xmm2			/* cpy t0 */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* cpy t3 */\n\t"\
			"addpd	%%xmm1,%%xmm0			/*~t0 = t0+tt */\n\t"\
			"addpd	%%xmm6,%%xmm5			/*~t3 = t3+t6 */\n\t"\
			"addpd	%%xmm3,%%xmm1			/*~tt = tt+t2 */\n\t"\
			"addpd	%%xmm4,%%xmm6			/*      t6+t5 */\n\t"\
			"addpd	%%xmm2,%%xmm3			/*~t2 = t2+t0 */\n\t"\
			"addpd	%%xmm7,%%xmm4			/*~t5 = t5+t3 */\n\t"\
			"subpd	%%xmm1,%%xmm2			/*~t1 = t0-tt-t2 */\n\t"\
			"subpd	%%xmm6,%%xmm7			/*~t4 = t3-t6-t5 */\n\t"\
			"movaps	%%xmm0,0x080(%%rdi)			/* B1 <- t0 = s1p04r */\n\t"\
			"movaps	%%xmm5,0x300(%%rdi)			/* B6 <- t3 = s1p24r */\n\t"\
			"movaps	%%xmm2,0x100(%%rdi)			/* B2 <- t1 = s1p08r */\n\t"\
			"movaps	%%xmm7,0x280(%%rdi)			/* B5 <- t4 = s1p20r */\n\t"\
			"movaps	%%xmm3,0x180(%%rdi)			/* B3 <- t2 = s1p12r */\n\t"\
			"movaps	%%xmm4,0x200(%%rdi)			/* B4 <- t5 = s1p16r */\n\t"\
			"/************************** Imaginary Parts: ******************************************/\n\t"\
			"movaps	0x310(%%rdi),%%xmm6		/* A1i = s1p24i */\n\t"\
			"movaps	0x090(%%rdi),%%xmm1		/* A6i = s1p04i */\n\t"\
			"movaps	0x290(%%rdi),%%xmm5		/* A2i = s1p20i */\n\t"\
			"movaps	0x110(%%rdi),%%xmm2		/* A5i = s1p08i */\n\t"\
			"movaps	0x210(%%rdi),%%xmm4		/* A3i = s1p16i */\n\t"\
			"movaps	0x190(%%rdi),%%xmm3		/* A4i = s1p12i */\n\t"\
			"\n\t"\
			"subpd	%%xmm1,%%xmm6			/* t6i = A1i-A6i */\n\t"\
			"addpd	%%xmm1,%%xmm1			/*         2*A6i */\n\t"\
			"addpd	%%xmm6,%%xmm1			/* t1i = A1i+A6i */\n\t"\
			"\n\t"\
			"subpd	%%xmm2,%%xmm5			/* t5i = A2i-A5i */\n\t"\
			"addpd	%%xmm2,%%xmm2			/*         2*A5i */\n\t"\
			"addpd	%%xmm5,%%xmm2			/* t2i = A2i+A5i */\n\t"\
			"\n\t"\
			"movaps	0x10(%%rdi),%%xmm0		/* Ai0 */\n\t"\
			"subpd	%%xmm3,%%xmm4			/* t4i = A3i-A4i */\n\t"\
			"addpd	%%xmm3,%%xmm3			/*         2*A4i */\n\t"\
			"addpd	%%xmm4,%%xmm3			/* t3i = A3i+A4i */\n\t"\
			"\n\t"\
			"movaps	%%xmm0,    (%%rsi)		/* cpy t0 into scratch sincos slot */\n\t"\
			"movaps	%%xmm6,0x10(%%rsi)		/* cpy t6 into scratch sincos slot */\n\t"\
			"addpd	%%xmm1,%%xmm0			/*~A0 = A0+t1 */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* cpy t5 */\n\t"\
			"addpd	%%xmm2,%%xmm3			/*~t3 = t3+t2 */\n\t"\
			"subpd	%%xmm4,%%xmm5			/*~t5 = t5-t4 */\n\t"\
			"subpd	%%xmm2,%%xmm1			/*~t1 = t1-t2 */\n\t"\
			"subpd	%%xmm7,%%xmm6			/*~t6 = t6-t5 */\n\t"\
			"addpd	%%xmm2,%%xmm2			/* 2*t2 */\n\t"\
			"addpd	%%xmm7,%%xmm4			/*~t5 = t4+t5 */\n\t"\
			"addpd	%%xmm3,%%xmm0			/* B0 */\n\t"\
			"addpd	0x10(%%rsi),%%xmm5		/* t3 = (t5-t4)+t6 */\n\t"\
			"subpd	%%xmm2,%%xmm3			/*~t2 =  (t2+t3) - 2*t2 = t3-t2 */\n\t"\
			"movaps	%%xmm4,%%xmm7			/* cpy t5 */\n\t"\
			"movaps	%%xmm0,0x010(%%rdi)		/* <-B0 = s1p00i, xmm0 FREE */\n\t"\
			"subpd	%%xmm6,%%xmm4			/* t4 = ~t5-~t6 */\n\t"\
			"movaps	%%xmm1,%%xmm2			/* cpy ~t1 */\n\t"\
			"subpd	    (%%rsi),%%xmm0		/* r = B0 - t0 */\n\t"\
			"mulpd	0x10(%%rax),%%xmm5		/*~t3 = t3*sx0 */\n\t"\
			"addpd	%%xmm3,%%xmm2			/* ~t1+~t2 */\n\t"\
			"mulpd	    (%%rcx),%%xmm3		/* t2 = t2*cx2 */\n\t"\
			"mulpd	0x10(%%rdx),%%xmm4		/*~t4 = t4*sx3 */\n\t"\
			"mulpd	    (%%r15),%%xmm1		/* t1 = t1*cx1 */\n\t"\
			"mulpd	0x10(%%r15),%%xmm6		/*~t6 = t6*sx1 */\n\t"\
			"mulpd	    (%%rax),%%xmm0		/* ~r = r*(cx0-1) */\n\t"\
			"mulpd	0x10(%%rcx),%%xmm7		/*~t5 = t5*sx2 */\n\t"\
			"mulpd	    (%%rdx),%%xmm2		/* t3 =(t1+t2)*cx3 */\n\t"\
			"addpd	0x10(%%rdi),%%xmm0		/* t0 =~r + B0 */\n\t"\
			"addpd	%%xmm4,%%xmm6			/*~t6 = t4+t6 */\n\t"\
			"subpd	%%xmm2,%%xmm1			/* tt = t1-t3 */\n\t"\
			"subpd	%%xmm7,%%xmm4			/*~t5 = t4-t5, xmm7 FREE */\n\t"\
			"subpd	%%xmm2,%%xmm3			/* t2 = t2-t3, xmm2 FREE */\n\t"\
			"movaps	%%xmm0,%%xmm2			/* cpy t0 */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* cpy t3 */\n\t"\
			"addpd	%%xmm1,%%xmm0			/*~t0 = t0+tt */\n\t"\
			"addpd	%%xmm6,%%xmm5			/*~t3 = t3+t6 */\n\t"\
			"addpd	%%xmm3,%%xmm1			/*~tt = tt+t2 */\n\t"\
			"addpd	%%xmm4,%%xmm6			/*      t6+t5 */\n\t"\
			"addpd	%%xmm2,%%xmm3			/*~t2 = t2+t0 */\n\t"\
			"addpd	%%xmm7,%%xmm4			/*~t5 = t5+t3 */\n\t"\
			"subpd	%%xmm1,%%xmm2			/*~t1 = t0-tt-t2, xmm1 FREE */\n\t"\
			"subpd	%%xmm6,%%xmm7			/*~t4 = t3-t6-t5, xmm6 FREE */\n\t"\
			"\n\t"\
			"/* xmm1,6 FREE */\n\t"\
			"movaps	0x080(%%rdi),%%xmm1		/* t0r */\n\t"\
			"movaps	0x300(%%rdi),%%xmm6		/* t3r */\n\t"\
			"subpd	%%xmm5,%%xmm1			/* B1r =t0r-t3i */\n\t"\
			"subpd	%%xmm6,%%xmm0			/* B6i =t0i-t3r */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*        2*t3i */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*        2*t3r */\n\t"\
			"addpd	%%xmm1,%%xmm5			/* B6r =t0r+t3i */\n\t"\
			"addpd	%%xmm0,%%xmm6			/* B1i =t0i+t3r */\n\t"\
			"movaps	%%xmm1,0x080(%%rdi)		/* <-B1r = s1p08r */\n\t"\
			"movaps	%%xmm0,0x310(%%rdi)		/* <-B6i = s1p20r */\n\t"\
			"movaps	%%xmm5,0x300(%%rdi)		/* <-B6r = s1p20r */\n\t"\
			"movaps	%%xmm6,0x090(%%rdi)		/* <-B1i = s1p08i */\n\t"\
			"\n\t"\
			"movaps	0x100(%%rdi),%%xmm1		/* t1r */\n\t"\
			"movaps	0x280(%%rdi),%%xmm6		/* t4r */\n\t"\
			"subpd	%%xmm7,%%xmm1			/* B2r =t1r-t4i */\n\t"\
			"subpd	%%xmm6,%%xmm2			/* B5i =t1i-t4r */\n\t"\
			"addpd	%%xmm7,%%xmm7			/*        2*t4i */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*        2*t4r */\n\t"\
			"addpd	%%xmm1,%%xmm7			/* B5r =t1r+t4i */\n\t"\
			"addpd	%%xmm2,%%xmm6			/* B2i =t1i+t4r */\n\t"\
			"movaps	%%xmm1,0x100(%%rdi)		/* <-B2r = s1p16r */\n\t"\
			"movaps	%%xmm2,0x290(%%rdi)		/* <-B5i = s1p12r */\n\t"\
			"movaps	%%xmm7,0x280(%%rdi)		/* <-B5r = s1p12r */\n\t"\
			"movaps	%%xmm6,0x110(%%rdi)		/* <-B2i = s1p16i */\n\t"\
			"\n\t"\
			"/* Note the order reversal on this pair of outputs: */\n\t"\
			"movaps	0x180(%%rdi),%%xmm0		/* t2r */\n\t"\
			"movaps	0x200(%%rdi),%%xmm5		/* t5r */\n\t"\
			"subpd	%%xmm4,%%xmm0			/* B4r =t2r-t5i */\n\t"\
			"subpd	%%xmm5,%%xmm3			/* B3i =t2i-t5r */\n\t"\
			"addpd	%%xmm4,%%xmm4			/*        2*t5i */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*        2*t5r */\n\t"\
			"addpd	%%xmm0,%%xmm4			/* B3r =t2r+t5i */\n\t"\
			"addpd	%%xmm3,%%xmm5			/* B4i =t2i+t5r */\n\t"\
			"movaps	%%xmm0,0x200(%%rdi)		/* <-B4r = s1p04r */\n\t"\
			"movaps	%%xmm3,0x190(%%rdi)		/* <-B3i = s1p24r */\n\t"\
			"movaps	%%xmm4,0x180(%%rdi)		/* <-B3r = s1p24r */\n\t"\
			"movaps	%%xmm5,0x210(%%rdi)		/* <-B4i = s1p04i */\n\t"\
			"\n\t"\
		"/* SSE2_RADIX_07_DFT(s1p21r,s1p17r,s1p13r,s1p09r,s1p05r,s1p01r,s1p25r,cc0,s1p01r,s1p05r,s1p09r,s1p13r,s1p17r,s1p21r,s1p25r): */\n\t"\
			"addq	$0x20,%%rdi	/* s1p01r */\n\t"\
			"movaps	0x200(%%rdi),%%xmm6			/* A1r = s1p17r */\n\t"\
			"movaps	0x300(%%rdi),%%xmm1			/* A6r = s1p25r */\n\t"\
			"movaps	0x180(%%rdi),%%xmm5			/* A2r = s1p13r */\n\t"\
			"movaps	     (%%rdi),%%xmm2			/* A5r = s1p01r */\n\t"\
			"movaps	0x100(%%rdi),%%xmm4			/* A3r = s1p09r */\n\t"\
			"movaps	0x080(%%rdi),%%xmm3			/* A4r = s1p05r */\n\t"\
			"\n\t"\
			"subpd	%%xmm1,%%xmm6			/* t6r = A1r-A6r */\n\t"\
			"addpd	%%xmm1,%%xmm1			/*         2*A6r */\n\t"\
			"addpd	%%xmm6,%%xmm1			/* t1r = A1r+A6r */\n\t"\
			"\n\t"\
			"subpd	%%xmm2,%%xmm5			/* t5r = A2r-A5r */\n\t"\
			"addpd	%%xmm2,%%xmm2			/*         2*A5r */\n\t"\
			"addpd	%%xmm5,%%xmm2			/* t2r = A2r+A5r */\n\t"\
			"\n\t"\
			"movaps	0x280(%%rdi),%%xmm0		/* Ar0 = s1p21r */\n\t"\
			"subpd	%%xmm3,%%xmm4			/* t4r = A3r-A4r */\n\t"\
			"addpd	%%xmm3,%%xmm3			/*         2*A4r */\n\t"\
			"addpd	%%xmm4,%%xmm3			/* t3r = A3r+A4r */\n\t"\
			"\n\t"\
			"movaps	%%xmm0,    (%%rsi)		/* cpy t0 into scratch sincos slot */\n\t"\
			"movaps	%%xmm6,0x10(%%rsi)		/* cpy t6 into scratch sincos slot */\n\t"\
			"addpd	%%xmm1,%%xmm0			/*~A0 = A0+t1 */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* cpy t5 */\n\t"\
			"addpd	%%xmm2,%%xmm3			/*~t3 = t3+t2 */\n\t"\
			"subpd	%%xmm4,%%xmm5			/*~t5 = t5-t4 */\n\t"\
			"subpd	%%xmm2,%%xmm1			/*~t1 = t1-t2 */\n\t"\
			"subpd	%%xmm7,%%xmm6			/*~t6 = t6-t5 */\n\t"\
			"addpd	%%xmm2,%%xmm2			/* 2*t2 */\n\t"\
			"addpd	%%xmm7,%%xmm4			/*~t5 = t4+t5 */\n\t"\
			"addpd	%%xmm3,%%xmm0			/* B0 */\n\t"\
			"addpd	0x10(%%rsi),%%xmm5		/* t3 = (t5-t4)+t6 */\n\t"\
			"subpd	%%xmm2,%%xmm3			/*~t2 =  (t2+t3) - 2*t2 = t3-t2 */\n\t"\
			"movaps	%%xmm4,%%xmm7			/* cpy t5 */\n\t"\
			"movaps	%%xmm0,     (%%rdi)		/* <-B0 = s1p01r, xmm0 FREE */\n\t"\
			"subpd	%%xmm6,%%xmm4			/* t4 = ~t5-~t6 */\n\t"\
			"movaps	%%xmm1,%%xmm2			/* cpy ~t1 */\n\t"\
			"subpd	    (%%rsi),%%xmm0		/* r = B0 - t0 */\n\t"\
			"mulpd	0x10(%%rax),%%xmm5		/*~t3 = t3*sx0 */\n\t"\
			"addpd	%%xmm3,%%xmm2			/* ~t1+~t2 */\n\t"\
			"mulpd	    (%%rcx),%%xmm3		/* t2 = t2*cx2 */\n\t"\
			"mulpd	0x10(%%rdx),%%xmm4		/*~t4 = t4*sx3 */\n\t"\
			"mulpd	    (%%r15),%%xmm1		/* t1 = t1*cx1 */\n\t"\
			"mulpd	0x10(%%r15),%%xmm6		/*~t6 = t6*sx1 */\n\t"\
			"mulpd	    (%%rax),%%xmm0		/* ~r = r*(cx0-1) */\n\t"\
			"mulpd	0x10(%%rcx),%%xmm7		/*~t5 = t5*sx2 */\n\t"\
			"mulpd	    (%%rdx),%%xmm2		/* t3 =(t1+t2)*cx3 */\n\t"\
			"addpd	    (%%rdi),%%xmm0		/* t0 =~r + B0 */\n\t"\
			"addpd	%%xmm4,%%xmm6			/*~t6 = t4+t6 */\n\t"\
			"subpd	%%xmm2,%%xmm1			/* tt = t1-t3 */\n\t"\
			"subpd	%%xmm7,%%xmm4			/*~t5 = t4-t5, xmm7 FREE */\n\t"\
			"subpd	%%xmm2,%%xmm3			/* t2 = t2-t3, xmm2 FREE */\n\t"\
			"\n\t"\
			"movaps	%%xmm0,%%xmm2			/* cpy t0 */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* cpy t3 */\n\t"\
			"addpd	%%xmm1,%%xmm0			/*~t0 = t0+tt */\n\t"\
			"addpd	%%xmm6,%%xmm5			/*~t3 = t3+t6 */\n\t"\
			"addpd	%%xmm3,%%xmm1			/*~tt = tt+t2 */\n\t"\
			"addpd	%%xmm4,%%xmm6			/*      t6+t5 */\n\t"\
			"addpd	%%xmm2,%%xmm3			/*~t2 = t2+t0 */\n\t"\
			"addpd	%%xmm7,%%xmm4			/*~t5 = t5+t3 */\n\t"\
			"subpd	%%xmm1,%%xmm2			/*~t1 = t0-tt-t2 */\n\t"\
			"subpd	%%xmm6,%%xmm7			/*~t4 = t3-t6-t5 */\n\t"\
			"movaps	%%xmm0,0x080(%%rdi)			/* B1 <- t0 = s1p05r */\n\t"\
			"movaps	%%xmm5,0x300(%%rdi)			/* B6 <- t3 = s1p25r */\n\t"\
			"movaps	%%xmm2,0x100(%%rdi)			/* B2 <- t1 = s1p09r */\n\t"\
			"movaps	%%xmm7,0x280(%%rdi)			/* B5 <- t4 = s1p21r */\n\t"\
			"movaps	%%xmm3,0x180(%%rdi)			/* B3 <- t2 = s1p13r */\n\t"\
			"movaps	%%xmm4,0x200(%%rdi)			/* B4 <- t5 = s1p17r */\n\t"\
			"/************************** Imaginary Parts: ******************************************/\n\t"\
			"movaps	0x210(%%rdi),%%xmm6		/* A1i = s1p17i */\n\t"\
			"movaps	0x310(%%rdi),%%xmm1		/* A6i = s1p25i */\n\t"\
			"movaps	0x190(%%rdi),%%xmm5		/* A2i = s1p13i */\n\t"\
			"movaps	0x010(%%rdi),%%xmm2		/* A5i = s1p01i */\n\t"\
			"movaps	0x110(%%rdi),%%xmm4		/* A3i = s1p09i */\n\t"\
			"movaps	0x090(%%rdi),%%xmm3		/* A4i = s1p05i */\n\t"\
			"\n\t"\
			"subpd	%%xmm1,%%xmm6			/* t6i = A1i-A6i */\n\t"\
			"addpd	%%xmm1,%%xmm1			/*         2*A6i */\n\t"\
			"addpd	%%xmm6,%%xmm1			/* t1i = A1i+A6i */\n\t"\
			"\n\t"\
			"subpd	%%xmm2,%%xmm5			/* t5i = A2i-A5i */\n\t"\
			"addpd	%%xmm2,%%xmm2			/*         2*A5i */\n\t"\
			"addpd	%%xmm5,%%xmm2			/* t2i = A2i+A5i */\n\t"\
			"\n\t"\
			"movaps	0x290(%%rdi),%%xmm0		/* Ai0 = s1p21i */\n\t"\
			"subpd	%%xmm3,%%xmm4			/* t4i = A3i-A4i */\n\t"\
			"addpd	%%xmm3,%%xmm3			/*         2*A4i */\n\t"\
			"addpd	%%xmm4,%%xmm3			/* t3i = A3i+A4i */\n\t"\
			"\n\t"\
			"movaps	%%xmm0,    (%%rsi)		/* cpy t0 into scratch sincos slot */\n\t"\
			"movaps	%%xmm6,0x10(%%rsi)		/* cpy t6 into scratch sincos slot */\n\t"\
			"addpd	%%xmm1,%%xmm0			/*~A0 = A0+t1 */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* cpy t5 */\n\t"\
			"addpd	%%xmm2,%%xmm3			/*~t3 = t3+t2 */\n\t"\
			"subpd	%%xmm4,%%xmm5			/*~t5 = t5-t4 */\n\t"\
			"subpd	%%xmm2,%%xmm1			/*~t1 = t1-t2 */\n\t"\
			"subpd	%%xmm7,%%xmm6			/*~t6 = t6-t5 */\n\t"\
			"addpd	%%xmm2,%%xmm2			/* 2*t2 */\n\t"\
			"addpd	%%xmm7,%%xmm4			/*~t5 = t4+t5 */\n\t"\
			"addpd	%%xmm3,%%xmm0			/* B0 */\n\t"\
			"addpd	0x10(%%rsi),%%xmm5		/* t3 = (t5-t4)+t6 */\n\t"\
			"subpd	%%xmm2,%%xmm3			/*~t2 =  (t2+t3) - 2*t2 = t3-t2 */\n\t"\
			"movaps	%%xmm4,%%xmm7			/* cpy t5 */\n\t"\
			"movaps	%%xmm0,0x010(%%rdi)		/* <-B0 = s1p00i, xmm0 FREE */\n\t"\
			"subpd	%%xmm6,%%xmm4			/* t4 = ~t5-~t6 */\n\t"\
			"movaps	%%xmm1,%%xmm2			/* cpy ~t1 */\n\t"\
			"subpd	    (%%rsi),%%xmm0		/* r = B0 - t0 */\n\t"\
			"mulpd	0x10(%%rax),%%xmm5		/*~t3 = t3*sx0 */\n\t"\
			"addpd	%%xmm3,%%xmm2			/* ~t1+~t2 */\n\t"\
			"mulpd	    (%%rcx),%%xmm3		/* t2 = t2*cx2 */\n\t"\
			"mulpd	0x10(%%rdx),%%xmm4		/*~t4 = t4*sx3 */\n\t"\
			"mulpd	    (%%r15),%%xmm1		/* t1 = t1*cx1 */\n\t"\
			"mulpd	0x10(%%r15),%%xmm6		/*~t6 = t6*sx1 */\n\t"\
			"mulpd	    (%%rax),%%xmm0		/* ~r = r*(cx0-1) */\n\t"\
			"mulpd	0x10(%%rcx),%%xmm7		/*~t5 = t5*sx2 */\n\t"\
			"mulpd	    (%%rdx),%%xmm2		/* t3 =(t1+t2)*cx3 */\n\t"\
			"addpd	0x10(%%rdi),%%xmm0		/* t0 =~r + B0 */\n\t"\
			"addpd	%%xmm4,%%xmm6			/*~t6 = t4+t6 */\n\t"\
			"subpd	%%xmm2,%%xmm1			/* tt = t1-t3 */\n\t"\
			"subpd	%%xmm7,%%xmm4			/*~t5 = t4-t5, xmm7 FREE */\n\t"\
			"subpd	%%xmm2,%%xmm3			/* t2 = t2-t3, xmm2 FREE */\n\t"\
			"movaps	%%xmm0,%%xmm2			/* cpy t0 */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* cpy t3 */\n\t"\
			"addpd	%%xmm1,%%xmm0			/*~t0 = t0+tt */\n\t"\
			"addpd	%%xmm6,%%xmm5			/*~t3 = t3+t6 */\n\t"\
			"addpd	%%xmm3,%%xmm1			/*~tt = tt+t2 */\n\t"\
			"addpd	%%xmm4,%%xmm6			/*      t6+t5 */\n\t"\
			"addpd	%%xmm2,%%xmm3			/*~t2 = t2+t0 */\n\t"\
			"addpd	%%xmm7,%%xmm4			/*~t5 = t5+t3 */\n\t"\
			"subpd	%%xmm1,%%xmm2			/*~t1 = t0-tt-t2, xmm1 FREE */\n\t"\
			"subpd	%%xmm6,%%xmm7			/*~t4 = t3-t6-t5, xmm6 FREE */\n\t"\
			"\n\t"\
			"/* xmm1,6 FREE */\n\t"\
			"movaps	0x080(%%rdi),%%xmm1		/* t0r */\n\t"\
			"movaps	0x300(%%rdi),%%xmm6		/* t3r */\n\t"\
			"subpd	%%xmm5,%%xmm1			/* B1r =t0r-t3i */\n\t"\
			"subpd	%%xmm6,%%xmm0			/* B6i =t0i-t3r */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*        2*t3i */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*        2*t3r */\n\t"\
			"addpd	%%xmm1,%%xmm5			/* B6r =t0r+t3i */\n\t"\
			"addpd	%%xmm0,%%xmm6			/* B1i =t0i+t3r */\n\t"\
			"movaps	%%xmm1,0x080(%%rdi)		/* <-B1r = s1p08r */\n\t"\
			"movaps	%%xmm0,0x310(%%rdi)		/* <-B6i = s1p20r */\n\t"\
			"movaps	%%xmm5,0x300(%%rdi)		/* <-B6r = s1p20r */\n\t"\
			"movaps	%%xmm6,0x090(%%rdi)		/* <-B1i = s1p08i */\n\t"\
			"\n\t"\
			"movaps	0x100(%%rdi),%%xmm1		/* t1r */\n\t"\
			"movaps	0x280(%%rdi),%%xmm6		/* t4r */\n\t"\
			"subpd	%%xmm7,%%xmm1			/* B2r =t1r-t4i */\n\t"\
			"subpd	%%xmm6,%%xmm2			/* B5i =t1i-t4r */\n\t"\
			"addpd	%%xmm7,%%xmm7			/*        2*t4i */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*        2*t4r */\n\t"\
			"addpd	%%xmm1,%%xmm7			/* B5r =t1r+t4i */\n\t"\
			"addpd	%%xmm2,%%xmm6			/* B2i =t1i+t4r */\n\t"\
			"movaps	%%xmm1,0x100(%%rdi)		/* <-B2r = s1p16r */\n\t"\
			"movaps	%%xmm2,0x290(%%rdi)		/* <-B5i = s1p12r */\n\t"\
			"movaps	%%xmm7,0x280(%%rdi)		/* <-B5r = s1p12r */\n\t"\
			"movaps	%%xmm6,0x110(%%rdi)		/* <-B2i = s1p16i */\n\t"\
			"\n\t"\
			"/* Note the order reversal on this pair of outputs: */\n\t"\
			"movaps	0x180(%%rdi),%%xmm0		/* t2r */\n\t"\
			"movaps	0x200(%%rdi),%%xmm5		/* t5r */\n\t"\
			"subpd	%%xmm4,%%xmm0			/* B4r =t2r-t5i */\n\t"\
			"subpd	%%xmm5,%%xmm3			/* B3i =t2i-t5r */\n\t"\
			"addpd	%%xmm4,%%xmm4			/*        2*t5i */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*        2*t5r */\n\t"\
			"addpd	%%xmm0,%%xmm4			/* B3r =t2r+t5i */\n\t"\
			"addpd	%%xmm3,%%xmm5			/* B4i =t2i+t5r */\n\t"\
			"movaps	%%xmm0,0x200(%%rdi)		/* <-B4r = s1p04r */\n\t"\
			"movaps	%%xmm3,0x190(%%rdi)		/* <-B3i = s1p24r */\n\t"\
			"movaps	%%xmm4,0x180(%%rdi)		/* <-B3r = s1p24r */\n\t"\
			"movaps	%%xmm5,0x210(%%rdi)		/* <-B4i = s1p04i */\n\t"\
			"\n\t"\
		"/* SSE2_RADIX_07_DFT(s1p14r,s1p10r,s1p06r,s1p02r,s1p26r,s1p22r,s1p18r,cc0,s1p02r,s1p06r,s1p10r,s1p14r,s1p18r,s1p22r,s1p26r): */\n\t"\
			"addq	$0x20,%%rdi	/* s1p02r */\n\t"\
			"movaps	0x100(%%rdi),%%xmm6			/* A1r = s1p10r */\n\t"\
			"movaps	0x200(%%rdi),%%xmm1			/* A6r = s1p18r */\n\t"\
			"movaps	0x080(%%rdi),%%xmm5			/* A2r = s1p06r */\n\t"\
			"movaps	0x280(%%rdi),%%xmm2			/* A5r = s1p22r */\n\t"\
			"movaps	     (%%rdi),%%xmm4			/* A3r = s1p02r */\n\t"\
			"movaps	0x300(%%rdi),%%xmm3			/* A4r = s1p26r */\n\t"\
			"\n\t"\
			"subpd	%%xmm1,%%xmm6			/* t6r = A1r-A6r */\n\t"\
			"addpd	%%xmm1,%%xmm1			/*         2*A6r */\n\t"\
			"addpd	%%xmm6,%%xmm1			/* t1r = A1r+A6r */\n\t"\
			"\n\t"\
			"subpd	%%xmm2,%%xmm5			/* t5r = A2r-A5r */\n\t"\
			"addpd	%%xmm2,%%xmm2			/*         2*A5r */\n\t"\
			"addpd	%%xmm5,%%xmm2			/* t2r = A2r+A5r */\n\t"\
			"\n\t"\
			"movaps	0x180(%%rdi),%%xmm0		/* Ar0 = s1p14r */\n\t"\
			"subpd	%%xmm3,%%xmm4			/* t4r = A3r-A4r */\n\t"\
			"addpd	%%xmm3,%%xmm3			/*         2*A4r */\n\t"\
			"addpd	%%xmm4,%%xmm3			/* t3r = A3r+A4r */\n\t"\
			"\n\t"\
			"movaps	%%xmm0,    (%%rsi)		/* cpy t0 into scratch sincos slot */\n\t"\
			"movaps	%%xmm6,0x10(%%rsi)		/* cpy t6 into scratch sincos slot */\n\t"\
			"addpd	%%xmm1,%%xmm0			/*~A0 = A0+t1 */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* cpy t5 */\n\t"\
			"addpd	%%xmm2,%%xmm3			/*~t3 = t3+t2 */\n\t"\
			"subpd	%%xmm4,%%xmm5			/*~t5 = t5-t4 */\n\t"\
			"subpd	%%xmm2,%%xmm1			/*~t1 = t1-t2 */\n\t"\
			"subpd	%%xmm7,%%xmm6			/*~t6 = t6-t5 */\n\t"\
			"addpd	%%xmm2,%%xmm2			/* 2*t2 */\n\t"\
			"addpd	%%xmm7,%%xmm4			/*~t5 = t4+t5 */\n\t"\
			"addpd	%%xmm3,%%xmm0			/* B0 */\n\t"\
			"addpd	0x10(%%rsi),%%xmm5		/* t3 = (t5-t4)+t6 */\n\t"\
			"subpd	%%xmm2,%%xmm3			/*~t2 =  (t2+t3) - 2*t2 = t3-t2 */\n\t"\
			"movaps	%%xmm4,%%xmm7			/* cpy t5 */\n\t"\
			"movaps	%%xmm0,     (%%rdi)		/* <-B0 = s1p01r, xmm0 FREE */\n\t"\
			"subpd	%%xmm6,%%xmm4			/* t4 = ~t5-~t6 */\n\t"\
			"movaps	%%xmm1,%%xmm2			/* cpy ~t1 */\n\t"\
			"subpd	    (%%rsi),%%xmm0		/* r = B0 - t0 */\n\t"\
			"mulpd	0x10(%%rax),%%xmm5		/*~t3 = t3*sx0 */\n\t"\
			"addpd	%%xmm3,%%xmm2			/* ~t1+~t2 */\n\t"\
			"mulpd	    (%%rcx),%%xmm3		/* t2 = t2*cx2 */\n\t"\
			"mulpd	0x10(%%rdx),%%xmm4		/*~t4 = t4*sx3 */\n\t"\
			"mulpd	    (%%r15),%%xmm1		/* t1 = t1*cx1 */\n\t"\
			"mulpd	0x10(%%r15),%%xmm6		/*~t6 = t6*sx1 */\n\t"\
			"mulpd	    (%%rax),%%xmm0		/* ~r = r*(cx0-1) */\n\t"\
			"mulpd	0x10(%%rcx),%%xmm7		/*~t5 = t5*sx2 */\n\t"\
			"mulpd	    (%%rdx),%%xmm2		/* t3 =(t1+t2)*cx3 */\n\t"\
			"addpd	    (%%rdi),%%xmm0		/* t0 =~r + B0 */\n\t"\
			"addpd	%%xmm4,%%xmm6			/*~t6 = t4+t6 */\n\t"\
			"subpd	%%xmm2,%%xmm1			/* tt = t1-t3 */\n\t"\
			"subpd	%%xmm7,%%xmm4			/*~t5 = t4-t5, xmm7 FREE */\n\t"\
			"subpd	%%xmm2,%%xmm3			/* t2 = t2-t3, xmm2 FREE */\n\t"\
			"\n\t"\
			"movaps	%%xmm0,%%xmm2			/* cpy t0 */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* cpy t3 */\n\t"\
			"addpd	%%xmm1,%%xmm0			/*~t0 = t0+tt */\n\t"\
			"addpd	%%xmm6,%%xmm5			/*~t3 = t3+t6 */\n\t"\
			"addpd	%%xmm3,%%xmm1			/*~tt = tt+t2 */\n\t"\
			"addpd	%%xmm4,%%xmm6			/*      t6+t5 */\n\t"\
			"addpd	%%xmm2,%%xmm3			/*~t2 = t2+t0 */\n\t"\
			"addpd	%%xmm7,%%xmm4			/*~t5 = t5+t3 */\n\t"\
			"subpd	%%xmm1,%%xmm2			/*~t1 = t0-tt-t2 */\n\t"\
			"subpd	%%xmm6,%%xmm7			/*~t4 = t3-t6-t5 */\n\t"\
			"movaps	%%xmm0,0x080(%%rdi)			/* B1 <- t0 = s1p05r */\n\t"\
			"movaps	%%xmm5,0x300(%%rdi)			/* B6 <- t3 = s1p25r */\n\t"\
			"movaps	%%xmm2,0x100(%%rdi)			/* B2 <- t1 = s1p09r */\n\t"\
			"movaps	%%xmm7,0x280(%%rdi)			/* B5 <- t4 = s1p21r */\n\t"\
			"movaps	%%xmm3,0x180(%%rdi)			/* B3 <- t2 = s1p13r */\n\t"\
			"movaps	%%xmm4,0x200(%%rdi)			/* B4 <- t5 = s1p17r */\n\t"\
			"/************************** Imaginary Parts: ******************************************/\n\t"\
			"movaps	0x110(%%rdi),%%xmm6		/* A1i = s1p10i */\n\t"\
			"movaps	0x210(%%rdi),%%xmm1		/* A6i = s1p18i */\n\t"\
			"movaps	0x090(%%rdi),%%xmm5		/* A2i = s1p06i */\n\t"\
			"movaps	0x290(%%rdi),%%xmm2		/* A5i = s1p22i */\n\t"\
			"movaps	0x010(%%rdi),%%xmm4		/* A3i = s1p02i */\n\t"\
			"movaps	0x310(%%rdi),%%xmm3		/* A4i = s1p26i */\n\t"\
			"\n\t"\
			"subpd	%%xmm1,%%xmm6			/* t6i = A1i-A6i */\n\t"\
			"addpd	%%xmm1,%%xmm1			/*         2*A6i */\n\t"\
			"addpd	%%xmm6,%%xmm1			/* t1i = A1i+A6i */\n\t"\
			"\n\t"\
			"subpd	%%xmm2,%%xmm5			/* t5i = A2i-A5i */\n\t"\
			"addpd	%%xmm2,%%xmm2			/*         2*A5i */\n\t"\
			"addpd	%%xmm5,%%xmm2			/* t2i = A2i+A5i */\n\t"\
			"\n\t"\
			"movaps	0x190(%%rdi),%%xmm0		/* Ai0 = s1p14i */\n\t"\
			"subpd	%%xmm3,%%xmm4			/* t4i = A3i-A4i */\n\t"\
			"addpd	%%xmm3,%%xmm3			/*         2*A4i */\n\t"\
			"addpd	%%xmm4,%%xmm3			/* t3i = A3i+A4i */\n\t"\
			"\n\t"\
			"movaps	%%xmm0,    (%%rsi)		/* cpy t0 into scratch sincos slot */\n\t"\
			"movaps	%%xmm6,0x10(%%rsi)		/* cpy t6 into scratch sincos slot */\n\t"\
			"addpd	%%xmm1,%%xmm0			/*~A0 = A0+t1 */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* cpy t5 */\n\t"\
			"addpd	%%xmm2,%%xmm3			/*~t3 = t3+t2 */\n\t"\
			"subpd	%%xmm4,%%xmm5			/*~t5 = t5-t4 */\n\t"\
			"subpd	%%xmm2,%%xmm1			/*~t1 = t1-t2 */\n\t"\
			"subpd	%%xmm7,%%xmm6			/*~t6 = t6-t5 */\n\t"\
			"addpd	%%xmm2,%%xmm2			/* 2*t2 */\n\t"\
			"addpd	%%xmm7,%%xmm4			/*~t5 = t4+t5 */\n\t"\
			"addpd	%%xmm3,%%xmm0			/* B0 */\n\t"\
			"addpd	0x10(%%rsi),%%xmm5		/* t3 = (t5-t4)+t6 */\n\t"\
			"subpd	%%xmm2,%%xmm3			/*~t2 =  (t2+t3) - 2*t2 = t3-t2 */\n\t"\
			"movaps	%%xmm4,%%xmm7			/* cpy t5 */\n\t"\
			"movaps	%%xmm0,0x010(%%rdi)		/* <-B0 = s1p00i, xmm0 FREE */\n\t"\
			"subpd	%%xmm6,%%xmm4			/* t4 = ~t5-~t6 */\n\t"\
			"movaps	%%xmm1,%%xmm2			/* cpy ~t1 */\n\t"\
			"subpd	    (%%rsi),%%xmm0		/* r = B0 - t0 */\n\t"\
			"mulpd	0x10(%%rax),%%xmm5		/*~t3 = t3*sx0 */\n\t"\
			"addpd	%%xmm3,%%xmm2			/* ~t1+~t2 */\n\t"\
			"mulpd	    (%%rcx),%%xmm3		/* t2 = t2*cx2 */\n\t"\
			"mulpd	0x10(%%rdx),%%xmm4		/*~t4 = t4*sx3 */\n\t"\
			"mulpd	    (%%r15),%%xmm1		/* t1 = t1*cx1 */\n\t"\
			"mulpd	0x10(%%r15),%%xmm6		/*~t6 = t6*sx1 */\n\t"\
			"mulpd	    (%%rax),%%xmm0		/* ~r = r*(cx0-1) */\n\t"\
			"mulpd	0x10(%%rcx),%%xmm7		/*~t5 = t5*sx2 */\n\t"\
			"mulpd	    (%%rdx),%%xmm2		/* t3 =(t1+t2)*cx3 */\n\t"\
			"addpd	0x10(%%rdi),%%xmm0		/* t0 =~r + B0 */\n\t"\
			"addpd	%%xmm4,%%xmm6			/*~t6 = t4+t6 */\n\t"\
			"subpd	%%xmm2,%%xmm1			/* tt = t1-t3 */\n\t"\
			"subpd	%%xmm7,%%xmm4			/*~t5 = t4-t5, xmm7 FREE */\n\t"\
			"subpd	%%xmm2,%%xmm3			/* t2 = t2-t3, xmm2 FREE */\n\t"\
			"movaps	%%xmm0,%%xmm2			/* cpy t0 */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* cpy t3 */\n\t"\
			"addpd	%%xmm1,%%xmm0			/*~t0 = t0+tt */\n\t"\
			"addpd	%%xmm6,%%xmm5			/*~t3 = t3+t6 */\n\t"\
			"addpd	%%xmm3,%%xmm1			/*~tt = tt+t2 */\n\t"\
			"addpd	%%xmm4,%%xmm6			/*      t6+t5 */\n\t"\
			"addpd	%%xmm2,%%xmm3			/*~t2 = t2+t0 */\n\t"\
			"addpd	%%xmm7,%%xmm4			/*~t5 = t5+t3 */\n\t"\
			"subpd	%%xmm1,%%xmm2			/*~t1 = t0-tt-t2, xmm1 FREE */\n\t"\
			"subpd	%%xmm6,%%xmm7			/*~t4 = t3-t6-t5, xmm6 FREE */\n\t"\
			"\n\t"\
			"/* xmm1,6 FREE */\n\t"\
			"movaps	0x080(%%rdi),%%xmm1		/* t0r */\n\t"\
			"movaps	0x300(%%rdi),%%xmm6		/* t3r */\n\t"\
			"subpd	%%xmm5,%%xmm1			/* B1r =t0r-t3i */\n\t"\
			"subpd	%%xmm6,%%xmm0			/* B6i =t0i-t3r */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*        2*t3i */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*        2*t3r */\n\t"\
			"addpd	%%xmm1,%%xmm5			/* B6r =t0r+t3i */\n\t"\
			"addpd	%%xmm0,%%xmm6			/* B1i =t0i+t3r */\n\t"\
			"movaps	%%xmm1,0x080(%%rdi)		/* <-B1r = s1p08r */\n\t"\
			"movaps	%%xmm0,0x310(%%rdi)		/* <-B6i = s1p20r */\n\t"\
			"movaps	%%xmm5,0x300(%%rdi)		/* <-B6r = s1p20r */\n\t"\
			"movaps	%%xmm6,0x090(%%rdi)		/* <-B1i = s1p08i */\n\t"\
			"\n\t"\
			"movaps	0x100(%%rdi),%%xmm1		/* t1r */\n\t"\
			"movaps	0x280(%%rdi),%%xmm6		/* t4r */\n\t"\
			"subpd	%%xmm7,%%xmm1			/* B2r =t1r-t4i */\n\t"\
			"subpd	%%xmm6,%%xmm2			/* B5i =t1i-t4r */\n\t"\
			"addpd	%%xmm7,%%xmm7			/*        2*t4i */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*        2*t4r */\n\t"\
			"addpd	%%xmm1,%%xmm7			/* B5r =t1r+t4i */\n\t"\
			"addpd	%%xmm2,%%xmm6			/* B2i =t1i+t4r */\n\t"\
			"movaps	%%xmm1,0x100(%%rdi)		/* <-B2r = s1p16r */\n\t"\
			"movaps	%%xmm2,0x290(%%rdi)		/* <-B5i = s1p12r */\n\t"\
			"movaps	%%xmm7,0x280(%%rdi)		/* <-B5r = s1p12r */\n\t"\
			"movaps	%%xmm6,0x110(%%rdi)		/* <-B2i = s1p16i */\n\t"\
			"\n\t"\
			"/* Note the order reversal on this pair of outputs: */\n\t"\
			"movaps	0x180(%%rdi),%%xmm0		/* t2r */\n\t"\
			"movaps	0x200(%%rdi),%%xmm5		/* t5r */\n\t"\
			"subpd	%%xmm4,%%xmm0			/* B4r =t2r-t5i */\n\t"\
			"subpd	%%xmm5,%%xmm3			/* B3i =t2i-t5r */\n\t"\
			"addpd	%%xmm4,%%xmm4			/*        2*t5i */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*        2*t5r */\n\t"\
			"addpd	%%xmm0,%%xmm4			/* B3r =t2r+t5i */\n\t"\
			"addpd	%%xmm3,%%xmm5			/* B4i =t2i+t5r */\n\t"\
			"movaps	%%xmm0,0x200(%%rdi)		/* <-B4r = s1p04r */\n\t"\
			"movaps	%%xmm3,0x190(%%rdi)		/* <-B3i = s1p24r */\n\t"\
			"movaps	%%xmm4,0x180(%%rdi)		/* <-B3r = s1p24r */\n\t"\
			"movaps	%%xmm5,0x210(%%rdi)		/* <-B4i = s1p04i */\n\t"\
			"\n\t"\
		"/* SSE2_RADIX_07_DFT(s1p07r,s1p03r,s1p27r,s1p23r,s1p19r,s1p15r,s1p11r,cc0,s1p03r,s1p07r,s1p11r,s1p15r,s1p19r,s1p23r,s1p27r): */\n\t"\
			"addq	$0x20,%%rdi	/* s1p03r */\n\t"\
			"movaps	     (%%rdi),%%xmm6			/* A1r = s1p03r */\n\t"\
			"movaps	0x100(%%rdi),%%xmm1			/* A6r = s1p11r */\n\t"\
			"movaps	0x300(%%rdi),%%xmm5			/* A2r = s1p27r */\n\t"\
			"movaps	0x180(%%rdi),%%xmm2			/* A5r = s1p15r */\n\t"\
			"movaps	0x280(%%rdi),%%xmm4			/* A3r = s1p23r */\n\t"\
			"movaps	0x200(%%rdi),%%xmm3			/* A4r = s1p19r */\n\t"\
			"\n\t"\
			"subpd	%%xmm1,%%xmm6			/* t6r = A1r-A6r */\n\t"\
			"addpd	%%xmm1,%%xmm1			/*         2*A6r */\n\t"\
			"addpd	%%xmm6,%%xmm1			/* t1r = A1r+A6r */\n\t"\
			"\n\t"\
			"subpd	%%xmm2,%%xmm5			/* t5r = A2r-A5r */\n\t"\
			"addpd	%%xmm2,%%xmm2			/*         2*A5r */\n\t"\
			"addpd	%%xmm5,%%xmm2			/* t2r = A2r+A5r */\n\t"\
			"\n\t"\
			"movaps	0x080(%%rdi),%%xmm0		/* Ar0 = s1p07r */\n\t"\
			"subpd	%%xmm3,%%xmm4			/* t4r = A3r-A4r */\n\t"\
			"addpd	%%xmm3,%%xmm3			/*         2*A4r */\n\t"\
			"addpd	%%xmm4,%%xmm3			/* t3r = A3r+A4r */\n\t"\
			"\n\t"\
			"movaps	%%xmm0,    (%%rsi)		/* cpy t0 into scratch sincos slot */\n\t"\
			"movaps	%%xmm6,0x10(%%rsi)		/* cpy t6 into scratch sincos slot */\n\t"\
			"addpd	%%xmm1,%%xmm0			/*~A0 = A0+t1 */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* cpy t5 */\n\t"\
			"addpd	%%xmm2,%%xmm3			/*~t3 = t3+t2 */\n\t"\
			"subpd	%%xmm4,%%xmm5			/*~t5 = t5-t4 */\n\t"\
			"subpd	%%xmm2,%%xmm1			/*~t1 = t1-t2 */\n\t"\
			"subpd	%%xmm7,%%xmm6			/*~t6 = t6-t5 */\n\t"\
			"addpd	%%xmm2,%%xmm2			/* 2*t2 */\n\t"\
			"addpd	%%xmm7,%%xmm4			/*~t5 = t4+t5 */\n\t"\
			"addpd	%%xmm3,%%xmm0			/* B0 */\n\t"\
			"addpd	0x10(%%rsi),%%xmm5		/* t3 = (t5-t4)+t6 */\n\t"\
			"subpd	%%xmm2,%%xmm3			/*~t2 =  (t2+t3) - 2*t2 = t3-t2 */\n\t"\
			"movaps	%%xmm4,%%xmm7			/* cpy t5 */\n\t"\
			"movaps	%%xmm0,     (%%rdi)		/* <-B0 = s1p01r, xmm0 FREE */\n\t"\
			"subpd	%%xmm6,%%xmm4			/* t4 = ~t5-~t6 */\n\t"\
			"movaps	%%xmm1,%%xmm2			/* cpy ~t1 */\n\t"\
			"subpd	    (%%rsi),%%xmm0		/* r = B0 - t0 */\n\t"\
			"mulpd	0x10(%%rax),%%xmm5		/*~t3 = t3*sx0 */\n\t"\
			"addpd	%%xmm3,%%xmm2			/* ~t1+~t2 */\n\t"\
			"mulpd	    (%%rcx),%%xmm3		/* t2 = t2*cx2 */\n\t"\
			"mulpd	0x10(%%rdx),%%xmm4		/*~t4 = t4*sx3 */\n\t"\
			"mulpd	    (%%r15),%%xmm1		/* t1 = t1*cx1 */\n\t"\
			"mulpd	0x10(%%r15),%%xmm6		/*~t6 = t6*sx1 */\n\t"\
			"mulpd	    (%%rax),%%xmm0		/* ~r = r*(cx0-1) */\n\t"\
			"mulpd	0x10(%%rcx),%%xmm7		/*~t5 = t5*sx2 */\n\t"\
			"mulpd	    (%%rdx),%%xmm2		/* t3 =(t1+t2)*cx3 */\n\t"\
			"addpd	    (%%rdi),%%xmm0		/* t0 =~r + B0 */\n\t"\
			"addpd	%%xmm4,%%xmm6			/*~t6 = t4+t6 */\n\t"\
			"subpd	%%xmm2,%%xmm1			/* tt = t1-t3 */\n\t"\
			"subpd	%%xmm7,%%xmm4			/*~t5 = t4-t5, xmm7 FREE */\n\t"\
			"subpd	%%xmm2,%%xmm3			/* t2 = t2-t3, xmm2 FREE */\n\t"\
			"\n\t"\
			"movaps	%%xmm0,%%xmm2			/* cpy t0 */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* cpy t3 */\n\t"\
			"addpd	%%xmm1,%%xmm0			/*~t0 = t0+tt */\n\t"\
			"addpd	%%xmm6,%%xmm5			/*~t3 = t3+t6 */\n\t"\
			"addpd	%%xmm3,%%xmm1			/*~tt = tt+t2 */\n\t"\
			"addpd	%%xmm4,%%xmm6			/*      t6+t5 */\n\t"\
			"addpd	%%xmm2,%%xmm3			/*~t2 = t2+t0 */\n\t"\
			"addpd	%%xmm7,%%xmm4			/*~t5 = t5+t3 */\n\t"\
			"subpd	%%xmm1,%%xmm2			/*~t1 = t0-tt-t2 */\n\t"\
			"subpd	%%xmm6,%%xmm7			/*~t4 = t3-t6-t5 */\n\t"\
			"movaps	%%xmm0,0x080(%%rdi)			/* B1 <- t0 = s1p05r */\n\t"\
			"movaps	%%xmm5,0x300(%%rdi)			/* B6 <- t3 = s1p25r */\n\t"\
			"movaps	%%xmm2,0x100(%%rdi)			/* B2 <- t1 = s1p09r */\n\t"\
			"movaps	%%xmm7,0x280(%%rdi)			/* B5 <- t4 = s1p21r */\n\t"\
			"movaps	%%xmm3,0x180(%%rdi)			/* B3 <- t2 = s1p13r */\n\t"\
			"movaps	%%xmm4,0x200(%%rdi)			/* B4 <- t5 = s1p17r */\n\t"\
			"/************************** Imaginary Parts: ******************************************/\n\t"\
			"movaps	0x010(%%rdi),%%xmm6			/* A1i = s1p03i */\n\t"\
			"movaps	0x110(%%rdi),%%xmm1			/* A6i = s1p11i */\n\t"\
			"movaps	0x310(%%rdi),%%xmm5			/* A2i = s1p27i */\n\t"\
			"movaps	0x190(%%rdi),%%xmm2			/* A5i = s1p15i */\n\t"\
			"movaps	0x290(%%rdi),%%xmm4			/* A3i = s1p23i */\n\t"\
			"movaps	0x210(%%rdi),%%xmm3			/* A4i = s1p19i */\n\t"\
			"\n\t"\
			"subpd	%%xmm1,%%xmm6			/* t6i = A1i-A6i */\n\t"\
			"addpd	%%xmm1,%%xmm1			/*         2*A6i */\n\t"\
			"addpd	%%xmm6,%%xmm1			/* t1i = A1i+A6i */\n\t"\
			"\n\t"\
			"subpd	%%xmm2,%%xmm5			/* t5i = A2i-A5i */\n\t"\
			"addpd	%%xmm2,%%xmm2			/*         2*A5i */\n\t"\
			"addpd	%%xmm5,%%xmm2			/* t2i = A2i+A5i */\n\t"\
			"\n\t"\
			"movaps	0x090(%%rdi),%%xmm0		/* Ai0 = s1p07i */\n\t"\
			"subpd	%%xmm3,%%xmm4			/* t4i = A3i-A4i */\n\t"\
			"addpd	%%xmm3,%%xmm3			/*         2*A4i */\n\t"\
			"addpd	%%xmm4,%%xmm3			/* t3i = A3i+A4i */\n\t"\
			"\n\t"\
			"movaps	%%xmm0,    (%%rsi)		/* cpy t0 into scratch sincos slot */\n\t"\
			"movaps	%%xmm6,0x10(%%rsi)		/* cpy t6 into scratch sincos slot */\n\t"\
			"addpd	%%xmm1,%%xmm0			/*~A0 = A0+t1 */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* cpy t5 */\n\t"\
			"addpd	%%xmm2,%%xmm3			/*~t3 = t3+t2 */\n\t"\
			"subpd	%%xmm4,%%xmm5			/*~t5 = t5-t4 */\n\t"\
			"subpd	%%xmm2,%%xmm1			/*~t1 = t1-t2 */\n\t"\
			"subpd	%%xmm7,%%xmm6			/*~t6 = t6-t5 */\n\t"\
			"addpd	%%xmm2,%%xmm2			/* 2*t2 */\n\t"\
			"addpd	%%xmm7,%%xmm4			/*~t5 = t4+t5 */\n\t"\
			"addpd	%%xmm3,%%xmm0			/* B0 */\n\t"\
			"addpd	0x10(%%rsi),%%xmm5		/* t3 = (t5-t4)+t6 */\n\t"\
			"subpd	%%xmm2,%%xmm3			/*~t2 =  (t2+t3) - 2*t2 = t3-t2 */\n\t"\
			"movaps	%%xmm4,%%xmm7			/* cpy t5 */\n\t"\
			"movaps	%%xmm0,0x010(%%rdi)		/* <-B0 = s1p00i, xmm0 FREE */\n\t"\
			"subpd	%%xmm6,%%xmm4			/* t4 = ~t5-~t6 */\n\t"\
			"movaps	%%xmm1,%%xmm2			/* cpy ~t1 */\n\t"\
			"subpd	    (%%rsi),%%xmm0		/* r = B0 - t0 */\n\t"\
			"mulpd	0x10(%%rax),%%xmm5		/*~t3 = t3*sx0 */\n\t"\
			"addpd	%%xmm3,%%xmm2			/* ~t1+~t2 */\n\t"\
			"mulpd	    (%%rcx),%%xmm3		/* t2 = t2*cx2 */\n\t"\
			"mulpd	0x10(%%rdx),%%xmm4		/*~t4 = t4*sx3 */\n\t"\
			"mulpd	    (%%r15),%%xmm1		/* t1 = t1*cx1 */\n\t"\
			"mulpd	0x10(%%r15),%%xmm6		/*~t6 = t6*sx1 */\n\t"\
			"mulpd	    (%%rax),%%xmm0		/* ~r = r*(cx0-1) */\n\t"\
			"mulpd	0x10(%%rcx),%%xmm7		/*~t5 = t5*sx2 */\n\t"\
			"mulpd	    (%%rdx),%%xmm2		/* t3 =(t1+t2)*cx3 */\n\t"\
			"addpd	0x10(%%rdi),%%xmm0		/* t0 =~r + B0 */\n\t"\
			"addpd	%%xmm4,%%xmm6			/*~t6 = t4+t6 */\n\t"\
			"subpd	%%xmm2,%%xmm1			/* tt = t1-t3 */\n\t"\
			"subpd	%%xmm7,%%xmm4			/*~t5 = t4-t5, xmm7 FREE */\n\t"\
			"subpd	%%xmm2,%%xmm3			/* t2 = t2-t3, xmm2 FREE */\n\t"\
			"movaps	%%xmm0,%%xmm2			/* cpy t0 */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* cpy t3 */\n\t"\
			"addpd	%%xmm1,%%xmm0			/*~t0 = t0+tt */\n\t"\
			"addpd	%%xmm6,%%xmm5			/*~t3 = t3+t6 */\n\t"\
			"addpd	%%xmm3,%%xmm1			/*~tt = tt+t2 */\n\t"\
			"addpd	%%xmm4,%%xmm6			/*      t6+t5 */\n\t"\
			"addpd	%%xmm2,%%xmm3			/*~t2 = t2+t0 */\n\t"\
			"addpd	%%xmm7,%%xmm4			/*~t5 = t5+t3 */\n\t"\
			"subpd	%%xmm1,%%xmm2			/*~t1 = t0-tt-t2, xmm1 FREE */\n\t"\
			"subpd	%%xmm6,%%xmm7			/*~t4 = t3-t6-t5, xmm6 FREE */\n\t"\
			"\n\t"\
			"/* xmm1,6 FREE */\n\t"\
			"movaps	0x080(%%rdi),%%xmm1		/* t0r */\n\t"\
			"movaps	0x300(%%rdi),%%xmm6		/* t3r */\n\t"\
			"subpd	%%xmm5,%%xmm1			/* B1r =t0r-t3i */\n\t"\
			"subpd	%%xmm6,%%xmm0			/* B6i =t0i-t3r */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*        2*t3i */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*        2*t3r */\n\t"\
			"addpd	%%xmm1,%%xmm5			/* B6r =t0r+t3i */\n\t"\
			"addpd	%%xmm0,%%xmm6			/* B1i =t0i+t3r */\n\t"\
			"movaps	%%xmm1,0x080(%%rdi)		/* <-B1r = s1p08r */\n\t"\
			"movaps	%%xmm0,0x310(%%rdi)		/* <-B6i = s1p20r */\n\t"\
			"movaps	%%xmm5,0x300(%%rdi)		/* <-B6r = s1p20r */\n\t"\
			"movaps	%%xmm6,0x090(%%rdi)		/* <-B1i = s1p08i */\n\t"\
			"\n\t"\
			"movaps	0x100(%%rdi),%%xmm1		/* t1r */\n\t"\
			"movaps	0x280(%%rdi),%%xmm6		/* t4r */\n\t"\
			"subpd	%%xmm7,%%xmm1			/* B2r =t1r-t4i */\n\t"\
			"subpd	%%xmm6,%%xmm2			/* B5i =t1i-t4r */\n\t"\
			"addpd	%%xmm7,%%xmm7			/*        2*t4i */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*        2*t4r */\n\t"\
			"addpd	%%xmm1,%%xmm7			/* B5r =t1r+t4i */\n\t"\
			"addpd	%%xmm2,%%xmm6			/* B2i =t1i+t4r */\n\t"\
			"movaps	%%xmm1,0x100(%%rdi)		/* <-B2r = s1p16r */\n\t"\
			"movaps	%%xmm2,0x290(%%rdi)		/* <-B5i = s1p12r */\n\t"\
			"movaps	%%xmm7,0x280(%%rdi)		/* <-B5r = s1p12r */\n\t"\
			"movaps	%%xmm6,0x110(%%rdi)		/* <-B2i = s1p16i */\n\t"\
			"\n\t"\
			"/* Note the order reversal on this pair of outputs: */\n\t"\
			"movaps	0x180(%%rdi),%%xmm0		/* t2r */\n\t"\
			"movaps	0x200(%%rdi),%%xmm5		/* t5r */\n\t"\
			"subpd	%%xmm4,%%xmm0			/* B4r =t2r-t5i */\n\t"\
			"subpd	%%xmm5,%%xmm3			/* B3i =t2i-t5r */\n\t"\
			"addpd	%%xmm4,%%xmm4			/*        2*t5i */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*        2*t5r */\n\t"\
			"addpd	%%xmm0,%%xmm4			/* B3r =t2r+t5i */\n\t"\
			"addpd	%%xmm3,%%xmm5			/* B4i =t2i+t5r */\n\t"\
			"movaps	%%xmm0,0x200(%%rdi)		/* <-B4r = s1p04r */\n\t"\
			"movaps	%%xmm3,0x190(%%rdi)		/* <-B3i = s1p24r */\n\t"\
			"movaps	%%xmm4,0x180(%%rdi)		/* <-B3r = s1p24r */\n\t"\
			"movaps	%%xmm5,0x210(%%rdi)		/* <-B4i = s1p04i */\n\t"\
			"\n\t"\
		"/* Since doing radix-7 in-place here, inputs of radix-4 are in consecutive memory locs, i.e. 0x20 bytes apart, e.g. the distance between s1p00r and s1p01r: */\n\t"\
			"\n\t"\
		"/* SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_E(s1p00r,s1p01r,s1p02r,s1p03r, add00+p[0,1,2,3]): */\n\t"\
			"movq	%[__out],%%rsi	/* s1p00r */\n\t"\
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
			"\n\t"\
			"movaps	    (%%rsi),%%xmm0		/* in0r = s1p00r */\n\t"\
			"movaps	0x20(%%rsi),%%xmm4		/* in1r = s1p01r */\n\t"\
			"movaps	0x10(%%rsi),%%xmm1		/* in0i = s1p00i */\n\t"\
			"movaps	0x30(%%rsi),%%xmm5		/* in1i = s1p01i */\n\t"\
			"movaps	%%xmm0,%%xmm2			/* xmm2 <- in0r */\n\t"\
			"movaps	%%xmm4,%%xmm6			/* xmm4 <- in1r */\n\t"\
			"movaps	%%xmm1,%%xmm3			/* xmm3 <- in0i */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* xmm5 <- in1i */\n\t"\
			"\n\t"\
			"addpd	0x40(%%rsi),%%xmm0		/* t1 */\n\t"\
			"addpd	0x60(%%rsi),%%xmm4		/* t5 */\n\t"\
			"addpd	0x50(%%rsi),%%xmm1		/* t2 */\n\t"\
			"addpd	0x70(%%rsi),%%xmm5		/* t6 */\n\t"\
			"subpd	0x40(%%rsi),%%xmm2		/* t3 */\n\t"\
			"subpd	0x60(%%rsi),%%xmm6		/* t7 */\n\t"\
			"subpd	0x50(%%rsi),%%xmm3		/* t4 */\n\t"\
			"subpd	0x70(%%rsi),%%xmm7		/* t8 */\n\t"\
			"/* Finish radix-4 butterfly and store results into output-array slots: */\n\t"\
			"subpd	%%xmm4,%%xmm0			/* ~t5 <- t1 -t5 */\n\t"\
			"subpd	%%xmm7,%%xmm2			/* ~t7 <- t3 -t8 */\n\t"\
			"subpd	%%xmm5,%%xmm1			/* ~t6 <- t2 -t6 */\n\t"\
			"subpd	%%xmm6,%%xmm3			/* ~t4 <- t4 -t7 */\n\t"\
			"movaps	%%xmm0,     (%%r15)		/* <- ~t5 */\n\t"\
			"movaps	%%xmm2,     (%%rcx)		/* <- ~t7 */\n\t"\
			"movaps	%%xmm1,0x010(%%r15)		/* <- ~t6 */\n\t"\
			"movaps	%%xmm3,0x010(%%rdx)		/* <- ~t4 */\n\t"\
			"addpd	%%xmm4,%%xmm4			/*          2*t5 */\n\t"\
			"addpd	%%xmm7,%%xmm7			/*          2*t8 */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*          2*t6 */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*          2*t7 */\n\t"\
			"addpd	%%xmm0,%%xmm4			/* ~t1 <- t1 +t5 */\n\t"\
			"addpd	%%xmm2,%%xmm7			/* ~t3 <- t3 +t8 */\n\t"\
			"addpd	%%xmm1,%%xmm5			/* ~t2 <- t2 +t6 */\n\t"\
			"addpd	%%xmm3,%%xmm6			/* ~t8 <- t4 +t7 */\n\t"\
			"movaps	%%xmm4,     (%%rax)		/* <- ~t1 */\n\t"\
			"movaps	%%xmm7,     (%%rdx)		/* <- ~t3 */\n\t"\
			"movaps	%%xmm5,0x010(%%rax)		/* <- ~t2 */\n\t"\
			"movaps	%%xmm6,0x010(%%rcx)		/* <- ~t8 */\n\t"\
		"/* SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_E(s1p04r,s1p05r,s1p06r,s1p07r, add24+p[1,0,3,2]): */\n\t"\
			"addq	$0x80,%%rsi		/* s1p04r */\n\t"\
			"movslq	%[__p24],%%rdi\n\t"\
			"shlq	$3,%%rdi		/* Pointer offset for floating doubles */\n\t"\
			"subq	%%rax,%%r15\n\t"\
			"subq	%%rax,%%rcx\n\t"\
			"subq	%%rax,%%rdx\n\t"\
			"addq	%%rdi,%%rax		/* add24 */\n\t"\
			"addq	%%rax,%%r15\n\t"\
			"addq	%%rax,%%rcx\n\t"\
			"addq	%%rax,%%rdx\n\t"\
			"\n\t"\
			"movaps	    (%%rsi),%%xmm0		/* in0r = s1p00r */\n\t"\
			"movaps	0x20(%%rsi),%%xmm4		/* in1r = s1p01r */\n\t"\
			"movaps	0x10(%%rsi),%%xmm1		/* in0i = s1p00i */\n\t"\
			"movaps	0x30(%%rsi),%%xmm5		/* in1i = s1p01i */\n\t"\
			"movaps	%%xmm0,%%xmm2			/* xmm2 <- in0r */\n\t"\
			"movaps	%%xmm4,%%xmm6			/* xmm4 <- in1r */\n\t"\
			"movaps	%%xmm1,%%xmm3			/* xmm3 <- in0i */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* xmm5 <- in1i */\n\t"\
			"\n\t"\
			"addpd	0x40(%%rsi),%%xmm0		/* t1 */\n\t"\
			"addpd	0x60(%%rsi),%%xmm4		/* t5 */\n\t"\
			"addpd	0x50(%%rsi),%%xmm1		/* t2 */\n\t"\
			"addpd	0x70(%%rsi),%%xmm5		/* t6 */\n\t"\
			"subpd	0x40(%%rsi),%%xmm2		/* t3 */\n\t"\
			"subpd	0x60(%%rsi),%%xmm6		/* t7 */\n\t"\
			"subpd	0x50(%%rsi),%%xmm3		/* t4 */\n\t"\
			"subpd	0x70(%%rsi),%%xmm7		/* t8 */\n\t"\
			"/* Finish radix-4 butterfly and store results into output-array slots: */\n\t"\
			"/* Swap _o1 <-> _o2, _o3 <-> _o4: eax <-> ebx, ecx <-> edx */\n\t"\
			"subpd	%%xmm4,%%xmm0			/* ~t5 <- t1 -t5 */\n\t"\
			"subpd	%%xmm7,%%xmm2			/* ~t7 <- t3 -t8 */\n\t"\
			"subpd	%%xmm5,%%xmm1			/* ~t6 <- t2 -t6 */\n\t"\
			"subpd	%%xmm6,%%xmm3			/* ~t4 <- t4 -t7 */\n\t"\
			"movaps	%%xmm0,     (%%rax)		/* <- ~t5 */\n\t"\
			"movaps	%%xmm2,     (%%rdx)		/* <- ~t7 */\n\t"\
			"movaps	%%xmm1,0x010(%%rax)		/* <- ~t6 */\n\t"\
			"movaps	%%xmm3,0x010(%%rcx)		/* <- ~t4 */\n\t"\
			"addpd	%%xmm4,%%xmm4			/*          2*t5 */\n\t"\
			"addpd	%%xmm7,%%xmm7			/*          2*t8 */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*          2*t6 */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*          2*t7 */\n\t"\
			"addpd	%%xmm0,%%xmm4			/* ~t1 <- t1 +t5 */\n\t"\
			"addpd	%%xmm2,%%xmm7			/* ~t3 <- t3 +t8 */\n\t"\
			"addpd	%%xmm1,%%xmm5			/* ~t2 <- t2 +t6 */\n\t"\
			"addpd	%%xmm3,%%xmm6			/* ~t8 <- t4 +t7 */\n\t"\
			"movaps	%%xmm4,     (%%r15)		/* <- ~t1 */\n\t"\
			"movaps	%%xmm7,     (%%rcx)		/* <- ~t3 */\n\t"\
			"movaps	%%xmm5,0x010(%%r15)		/* <- ~t2 */\n\t"\
			"movaps	%%xmm6,0x010(%%rdx)		/* <- ~t8 */\n\t"\
		"/* SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_E(s1p08r,s1p09r,s1p10r,s1p11r, add20+p[2,3,1,0]): */\n\t"\
			"addq	$0x80,%%rsi		/* s1p08r */\n\t"\
			"subq	%%rax,%%r15\n\t"\
			"subq	%%rax,%%rcx\n\t"\
			"subq	%%rax,%%rdx\n\t"\
			"subq	%%rdi,%%rax		/* add0  */\n\t"\
			"movslq	%[__p20],%%rdi\n\t"\
			"shlq	$3,%%rdi		/* Pointer offset for floating doubles */\n\t"\
			"addq	%%rdi,%%rax		/* add20 */\n\t"\
			"addq	%%rax,%%r15\n\t"\
			"addq	%%rax,%%rcx\n\t"\
			"addq	%%rax,%%rdx\n\t"\
			"\n\t"\
			"movaps	    (%%rsi),%%xmm0		/* in0r = s1p00r */\n\t"\
			"movaps	0x20(%%rsi),%%xmm4		/* in1r = s1p01r */\n\t"\
			"movaps	0x10(%%rsi),%%xmm1		/* in0i = s1p00i */\n\t"\
			"movaps	0x30(%%rsi),%%xmm5		/* in1i = s1p01i */\n\t"\
			"movaps	%%xmm0,%%xmm2			/* xmm2 <- in0r */\n\t"\
			"movaps	%%xmm4,%%xmm6			/* xmm4 <- in1r */\n\t"\
			"movaps	%%xmm1,%%xmm3			/* xmm3 <- in0i */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* xmm5 <- in1i */\n\t"\
			"\n\t"\
			"addpd	0x40(%%rsi),%%xmm0		/* t1 */\n\t"\
			"addpd	0x60(%%rsi),%%xmm4		/* t5 */\n\t"\
			"addpd	0x50(%%rsi),%%xmm1		/* t2 */\n\t"\
			"addpd	0x70(%%rsi),%%xmm5		/* t6 */\n\t"\
			"subpd	0x40(%%rsi),%%xmm2		/* t3 */\n\t"\
			"subpd	0x60(%%rsi),%%xmm6		/* t7 */\n\t"\
			"subpd	0x50(%%rsi),%%xmm3		/* t4 */\n\t"\
			"subpd	0x70(%%rsi),%%xmm7		/* t8 */\n\t"\
			"/* Finish radix-4 butterfly and store results into output-array slots: */\n\t"\
			"/* Swap _o0123 <-> _o2310: e[abcd]x <-> e[cdba]x */\n\t"\
			"subpd	%%xmm4,%%xmm0			/* ~t5 <- t1 -t5 */\n\t"\
			"subpd	%%xmm7,%%xmm2			/* ~t7 <- t3 -t8 */\n\t"\
			"subpd	%%xmm5,%%xmm1			/* ~t6 <- t2 -t6 */\n\t"\
			"subpd	%%xmm6,%%xmm3			/* ~t4 <- t4 -t7 */\n\t"\
			"movaps	%%xmm0,     (%%rdx)		/* <- ~t5 */\n\t"\
			"movaps	%%xmm2,     (%%r15)		/* <- ~t7 */\n\t"\
			"movaps	%%xmm1,0x010(%%rdx)		/* <- ~t6 */\n\t"\
			"movaps	%%xmm3,0x010(%%rax)		/* <- ~t4 */\n\t"\
			"addpd	%%xmm4,%%xmm4			/*          2*t5 */\n\t"\
			"addpd	%%xmm7,%%xmm7			/*          2*t8 */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*          2*t6 */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*          2*t7 */\n\t"\
			"addpd	%%xmm0,%%xmm4			/* ~t1 <- t1 +t5 */\n\t"\
			"addpd	%%xmm2,%%xmm7			/* ~t3 <- t3 +t8 */\n\t"\
			"addpd	%%xmm1,%%xmm5			/* ~t2 <- t2 +t6 */\n\t"\
			"addpd	%%xmm3,%%xmm6			/* ~t8 <- t4 +t7 */\n\t"\
			"movaps	%%xmm4,     (%%rcx)		/* <- ~t1 */\n\t"\
			"movaps	%%xmm7,     (%%rax)		/* <- ~t3 */\n\t"\
			"movaps	%%xmm5,0x010(%%rcx)		/* <- ~t2 */\n\t"\
			"movaps	%%xmm6,0x010(%%r15)		/* <- ~t8 */\n\t"\
		"/* SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_E(s1p12r,s1p13r,s1p14r,s1p15r, add16+p[0,1,2,3]): */\n\t"\
			"addq	$0x80,%%rsi		/* s1p12r */\n\t"\
			"subq	%%rax,%%r15\n\t"\
			"subq	%%rax,%%rcx\n\t"\
			"subq	%%rax,%%rdx\n\t"\
			"subq	%%rdi,%%rax		/* add0  */\n\t"\
			"movslq	%[__p16],%%rdi\n\t"\
			"shlq	$3,%%rdi		/* Pointer offset for floating doubles */\n\t"\
			"addq	%%rdi,%%rax		/* add16 */\n\t"\
			"addq	%%rax,%%r15\n\t"\
			"addq	%%rax,%%rcx\n\t"\
			"addq	%%rax,%%rdx\n\t"\
			"\n\t"\
			"movaps	    (%%rsi),%%xmm0		/* in0r = s1p00r */\n\t"\
			"movaps	0x20(%%rsi),%%xmm4		/* in1r = s1p01r */\n\t"\
			"movaps	0x10(%%rsi),%%xmm1		/* in0i = s1p00i */\n\t"\
			"movaps	0x30(%%rsi),%%xmm5		/* in1i = s1p01i */\n\t"\
			"movaps	%%xmm0,%%xmm2			/* xmm2 <- in0r */\n\t"\
			"movaps	%%xmm4,%%xmm6			/* xmm4 <- in1r */\n\t"\
			"movaps	%%xmm1,%%xmm3			/* xmm3 <- in0i */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* xmm5 <- in1i */\n\t"\
			"\n\t"\
			"addpd	0x40(%%rsi),%%xmm0		/* t1 */\n\t"\
			"addpd	0x60(%%rsi),%%xmm4		/* t5 */\n\t"\
			"addpd	0x50(%%rsi),%%xmm1		/* t2 */\n\t"\
			"addpd	0x70(%%rsi),%%xmm5		/* t6 */\n\t"\
			"subpd	0x40(%%rsi),%%xmm2		/* t3 */\n\t"\
			"subpd	0x60(%%rsi),%%xmm6		/* t7 */\n\t"\
			"subpd	0x50(%%rsi),%%xmm3		/* t4 */\n\t"\
			"subpd	0x70(%%rsi),%%xmm7		/* t8 */\n\t"\
			"/* Finish radix-4 butterfly and store results into output-array slots: */\n\t"\
			"subpd	%%xmm4,%%xmm0			/* ~t5 <- t1 -t5 */\n\t"\
			"subpd	%%xmm7,%%xmm2			/* ~t7 <- t3 -t8 */\n\t"\
			"subpd	%%xmm5,%%xmm1			/* ~t6 <- t2 -t6 */\n\t"\
			"subpd	%%xmm6,%%xmm3			/* ~t4 <- t4 -t7 */\n\t"\
			"movaps	%%xmm0,     (%%r15)		/* <- ~t5 */\n\t"\
			"movaps	%%xmm2,     (%%rcx)		/* <- ~t7 */\n\t"\
			"movaps	%%xmm1,0x010(%%r15)		/* <- ~t6 */\n\t"\
			"movaps	%%xmm3,0x010(%%rdx)		/* <- ~t4 */\n\t"\
			"addpd	%%xmm4,%%xmm4			/*          2*t5 */\n\t"\
			"addpd	%%xmm7,%%xmm7			/*          2*t8 */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*          2*t6 */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*          2*t7 */\n\t"\
			"addpd	%%xmm0,%%xmm4			/* ~t1 <- t1 +t5 */\n\t"\
			"addpd	%%xmm2,%%xmm7			/* ~t3 <- t3 +t8 */\n\t"\
			"addpd	%%xmm1,%%xmm5			/* ~t2 <- t2 +t6 */\n\t"\
			"addpd	%%xmm3,%%xmm6			/* ~t8 <- t4 +t7 */\n\t"\
			"movaps	%%xmm4,     (%%rax)		/* <- ~t1 */\n\t"\
			"movaps	%%xmm7,     (%%rdx)		/* <- ~t3 */\n\t"\
			"movaps	%%xmm5,0x010(%%rax)		/* <- ~t2 */\n\t"\
			"movaps	%%xmm6,0x010(%%rcx)		/* <- ~t8 */\n\t"\
		"/* SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_E(s1p16r,s1p17r,s1p18r,s1p19r, add12+p[3,2,0,1]): */\n\t"\
			"addq	$0x80,%%rsi		/* s1p16r */\n\t"\
			"subq	%%rax,%%r15\n\t"\
			"subq	%%rax,%%rcx\n\t"\
			"subq	%%rax,%%rdx\n\t"\
			"subq	%%rdi,%%rax		/* add0  */\n\t"\
			"movslq	%[__p12],%%rdi\n\t"\
			"shlq	$3,%%rdi		/* Pointer offset for floating doubles */\n\t"\
			"addq	%%rdi,%%rax		/* add12 */\n\t"\
			"addq	%%rax,%%r15\n\t"\
			"addq	%%rax,%%rcx\n\t"\
			"addq	%%rax,%%rdx\n\t"\
			"\n\t"\
			"movaps	    (%%rsi),%%xmm0		/* in0r = s1p00r */\n\t"\
			"movaps	0x20(%%rsi),%%xmm4		/* in1r = s1p01r */\n\t"\
			"movaps	0x10(%%rsi),%%xmm1		/* in0i = s1p00i */\n\t"\
			"movaps	0x30(%%rsi),%%xmm5		/* in1i = s1p01i */\n\t"\
			"movaps	%%xmm0,%%xmm2			/* xmm2 <- in0r */\n\t"\
			"movaps	%%xmm4,%%xmm6			/* xmm4 <- in1r */\n\t"\
			"movaps	%%xmm1,%%xmm3			/* xmm3 <- in0i */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* xmm5 <- in1i */\n\t"\
			"\n\t"\
			"addpd	0x40(%%rsi),%%xmm0		/* t1 */\n\t"\
			"addpd	0x60(%%rsi),%%xmm4		/* t5 */\n\t"\
			"addpd	0x50(%%rsi),%%xmm1		/* t2 */\n\t"\
			"addpd	0x70(%%rsi),%%xmm5		/* t6 */\n\t"\
			"subpd	0x40(%%rsi),%%xmm2		/* t3 */\n\t"\
			"subpd	0x60(%%rsi),%%xmm6		/* t7 */\n\t"\
			"subpd	0x50(%%rsi),%%xmm3		/* t4 */\n\t"\
			"subpd	0x70(%%rsi),%%xmm7		/* t8 */\n\t"\
			"/* Finish radix-4 butterfly and store results into output-array slots: */\n\t"\
			"/* Swap _o0123 <-> _o3201: e[abcd]x <-> e[dcab]x */\n\t"\
			"subpd	%%xmm4,%%xmm0			/* ~t5 <- t1 -t5 */\n\t"\
			"subpd	%%xmm7,%%xmm2			/* ~t7 <- t3 -t8 */\n\t"\
			"subpd	%%xmm5,%%xmm1			/* ~t6 <- t2 -t6 */\n\t"\
			"subpd	%%xmm6,%%xmm3			/* ~t4 <- t4 -t7 */\n\t"\
			"movaps	%%xmm0,     (%%rcx)		/* <- ~t5 */\n\t"\
			"movaps	%%xmm2,     (%%rax)		/* <- ~t7 */\n\t"\
			"movaps	%%xmm1,0x010(%%rcx)		/* <- ~t6 */\n\t"\
			"movaps	%%xmm3,0x010(%%r15)		/* <- ~t4 */\n\t"\
			"addpd	%%xmm4,%%xmm4			/*          2*t5 */\n\t"\
			"addpd	%%xmm7,%%xmm7			/*          2*t8 */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*          2*t6 */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*          2*t7 */\n\t"\
			"addpd	%%xmm0,%%xmm4			/* ~t1 <- t1 +t5 */\n\t"\
			"addpd	%%xmm2,%%xmm7			/* ~t3 <- t3 +t8 */\n\t"\
			"addpd	%%xmm1,%%xmm5			/* ~t2 <- t2 +t6 */\n\t"\
			"addpd	%%xmm3,%%xmm6			/* ~t8 <- t4 +t7 */\n\t"\
			"movaps	%%xmm4,     (%%rdx)		/* <- ~t1 */\n\t"\
			"movaps	%%xmm7,     (%%r15)		/* <- ~t3 */\n\t"\
			"movaps	%%xmm5,0x010(%%rdx)		/* <- ~t2 */\n\t"\
			"movaps	%%xmm6,0x010(%%rax)		/* <- ~t8 */\n\t"\
		"/* SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_E(s1p20r,s1p21r,s1p22r,s1p23r, add08+p[1,0,3,2]): */\n\t"\
			"addq	$0x80,%%rsi		/* s1p20r */\n\t"\
			"subq	%%rax,%%r15\n\t"\
			"subq	%%rax,%%rcx\n\t"\
			"subq	%%rax,%%rdx\n\t"\
			"subq	%%rdi,%%rax		/* add0  */\n\t"\
			"movslq	%[__p08],%%rdi\n\t"\
			"shlq	$3,%%rdi		/* Pointer offset for floating doubles */\n\t"\
			"addq	%%rdi,%%rax		/* add08 */\n\t"\
			"addq	%%rax,%%r15\n\t"\
			"addq	%%rax,%%rcx\n\t"\
			"addq	%%rax,%%rdx\n\t"\
			"\n\t"\
			"movaps	    (%%rsi),%%xmm0		/* in0r = s1p00r */\n\t"\
			"movaps	0x20(%%rsi),%%xmm4		/* in1r = s1p01r */\n\t"\
			"movaps	0x10(%%rsi),%%xmm1		/* in0i = s1p00i */\n\t"\
			"movaps	0x30(%%rsi),%%xmm5		/* in1i = s1p01i */\n\t"\
			"movaps	%%xmm0,%%xmm2			/* xmm2 <- in0r */\n\t"\
			"movaps	%%xmm4,%%xmm6			/* xmm4 <- in1r */\n\t"\
			"movaps	%%xmm1,%%xmm3			/* xmm3 <- in0i */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* xmm5 <- in1i */\n\t"\
			"\n\t"\
			"addpd	0x40(%%rsi),%%xmm0		/* t1 */\n\t"\
			"addpd	0x60(%%rsi),%%xmm4		/* t5 */\n\t"\
			"addpd	0x50(%%rsi),%%xmm1		/* t2 */\n\t"\
			"addpd	0x70(%%rsi),%%xmm5		/* t6 */\n\t"\
			"subpd	0x40(%%rsi),%%xmm2		/* t3 */\n\t"\
			"subpd	0x60(%%rsi),%%xmm6		/* t7 */\n\t"\
			"subpd	0x50(%%rsi),%%xmm3		/* t4 */\n\t"\
			"subpd	0x70(%%rsi),%%xmm7		/* t8 */\n\t"\
			"/* Finish radix-4 butterfly and store results into output-array slots: */\n\t"\
			"/* Swap _o0 <-> _o1, _o2 <-> _o3: eax <-> ebx, ecx <-> edx */\n\t"\
			"subpd	%%xmm4,%%xmm0			/* ~t5 <- t1 -t5 */\n\t"\
			"subpd	%%xmm7,%%xmm2			/* ~t7 <- t3 -t8 */\n\t"\
			"subpd	%%xmm5,%%xmm1			/* ~t6 <- t2 -t6 */\n\t"\
			"subpd	%%xmm6,%%xmm3			/* ~t4 <- t4 -t7 */\n\t"\
			"movaps	%%xmm0,     (%%rax)		/* <- ~t5 */\n\t"\
			"movaps	%%xmm2,     (%%rdx)		/* <- ~t7 */\n\t"\
			"movaps	%%xmm1,0x010(%%rax)		/* <- ~t6 */\n\t"\
			"movaps	%%xmm3,0x010(%%rcx)		/* <- ~t4 */\n\t"\
			"addpd	%%xmm4,%%xmm4			/*          2*t5 */\n\t"\
			"addpd	%%xmm7,%%xmm7			/*          2*t8 */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*          2*t6 */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*          2*t7 */\n\t"\
			"addpd	%%xmm0,%%xmm4			/* ~t1 <- t1 +t5 */\n\t"\
			"addpd	%%xmm2,%%xmm7			/* ~t3 <- t3 +t8 */\n\t"\
			"addpd	%%xmm1,%%xmm5			/* ~t2 <- t2 +t6 */\n\t"\
			"addpd	%%xmm3,%%xmm6			/* ~t8 <- t4 +t7 */\n\t"\
			"movaps	%%xmm4,     (%%r15)		/* <- ~t1 */\n\t"\
			"movaps	%%xmm7,     (%%rcx)		/* <- ~t3 */\n\t"\
			"movaps	%%xmm5,0x010(%%r15)		/* <- ~t2 */\n\t"\
			"movaps	%%xmm6,0x010(%%rdx)		/* <- ~t8 */\n\t"\
		"/* SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_E(s1p24r,s1p25r,s1p26r,s1p27r, add04+p[2,3,1,0]): */\n\t"\
			"addq	$0x80,%%rsi		/* s1p24r */\n\t"\
			"subq	%%rax,%%r15\n\t"\
			"subq	%%rax,%%rcx\n\t"\
			"subq	%%rax,%%rdx\n\t"\
			"subq	%%rdi,%%rax		/* add0  */\n\t"\
			"movslq	%[__p04],%%rdi\n\t"\
			"shlq	$3,%%rdi		/* Pointer offset for floating doubles */\n\t"\
			"addq	%%rdi,%%rax		/* add04 */\n\t"\
			"addq	%%rax,%%r15\n\t"\
			"addq	%%rax,%%rcx\n\t"\
			"addq	%%rax,%%rdx\n\t"\
			"\n\t"\
			"movaps	    (%%rsi),%%xmm0		/* in0r = s1p00r */\n\t"\
			"movaps	0x20(%%rsi),%%xmm4		/* in1r = s1p01r */\n\t"\
			"movaps	0x10(%%rsi),%%xmm1		/* in0i = s1p00i */\n\t"\
			"movaps	0x30(%%rsi),%%xmm5		/* in1i = s1p01i */\n\t"\
			"movaps	%%xmm0,%%xmm2			/* xmm2 <- in0r */\n\t"\
			"movaps	%%xmm4,%%xmm6			/* xmm4 <- in1r */\n\t"\
			"movaps	%%xmm1,%%xmm3			/* xmm3 <- in0i */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* xmm5 <- in1i */\n\t"\
			"\n\t"\
			"addpd	0x40(%%rsi),%%xmm0		/* t1 */\n\t"\
			"addpd	0x60(%%rsi),%%xmm4		/* t5 */\n\t"\
			"addpd	0x50(%%rsi),%%xmm1		/* t2 */\n\t"\
			"addpd	0x70(%%rsi),%%xmm5		/* t6 */\n\t"\
			"subpd	0x40(%%rsi),%%xmm2		/* t3 */\n\t"\
			"subpd	0x60(%%rsi),%%xmm6		/* t7 */\n\t"\
			"subpd	0x50(%%rsi),%%xmm3		/* t4 */\n\t"\
			"subpd	0x70(%%rsi),%%xmm7		/* t8 */\n\t"\
			"/* Finish radix-4 butterfly and store results into output-array slots: */\n\t"\
			"/* Swap _o0123 <-> _o2310: e[abcd]x <-> e[cdba]x */\n\t"\
			"subpd	%%xmm4,%%xmm0			/* ~t5 <- t1 -t5 */\n\t"\
			"subpd	%%xmm7,%%xmm2			/* ~t7 <- t3 -t8 */\n\t"\
			"subpd	%%xmm5,%%xmm1			/* ~t6 <- t2 -t6 */\n\t"\
			"subpd	%%xmm6,%%xmm3			/* ~t4 <- t4 -t7 */\n\t"\
			"movaps	%%xmm0,     (%%rdx)		/* <- ~t5 */\n\t"\
			"movaps	%%xmm2,     (%%r15)		/* <- ~t7 */\n\t"\
			"movaps	%%xmm1,0x010(%%rdx)		/* <- ~t6 */\n\t"\
			"movaps	%%xmm3,0x010(%%rax)		/* <- ~t4 */\n\t"\
			"addpd	%%xmm4,%%xmm4			/*          2*t5 */\n\t"\
			"addpd	%%xmm7,%%xmm7			/*          2*t8 */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*          2*t6 */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*          2*t7 */\n\t"\
			"addpd	%%xmm0,%%xmm4			/* ~t1 <- t1 +t5 */\n\t"\
			"addpd	%%xmm2,%%xmm7			/* ~t3 <- t3 +t8 */\n\t"\
			"addpd	%%xmm1,%%xmm5			/* ~t2 <- t2 +t6 */\n\t"\
			"addpd	%%xmm3,%%xmm6			/* ~t8 <- t4 +t7 */\n\t"\
			"movaps	%%xmm4,     (%%rcx)		/* <- ~t1 */\n\t"\
			"movaps	%%xmm7,     (%%rax)		/* <- ~t3 */\n\t"\
			"movaps	%%xmm5,0x010(%%rcx)		/* <- ~t2 */\n\t"\
			"movaps	%%xmm6,0x010(%%r15)		/* <- ~t8 */\n\t"\
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
			 ,[__out] "m" (Xout)\
			 ,[__cc0] "m" (Xcc0)\
			: "rax","r15","rcx","rdx","rdi","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
		);\
	}

#endif	/* radix28_ditN_cy_dif1_gcc_h_included */
