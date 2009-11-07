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
			"movl	%[__out],%%esi	/* s1p00r */\n\t"\
			"movl	%[__add],%%eax\n\t"\
			"movl	%[__p01],%%ebx\n\t"\
			"movl	%[__p02],%%ecx\n\t"\
			"movl	%[__p03],%%edx\n\t"\
			"shll	$3,%%ebx		/* Pointer offset for floating doubles */\n\t"\
			"shll	$3,%%ecx\n\t"\
			"shll	$3,%%edx\n\t"\
			"addl	%%eax,%%ebx\n\t"\
			"addl	%%eax,%%ecx\n\t"\
			"addl	%%eax,%%edx\n\t"\
			"/* ecx <-> edx */\n\t"\
			"movaps	    (%%eax),%%xmm0		/* a(jt   ) */\n\t"\
			"movaps	    (%%edx),%%xmm4		/* a(jt+p2) */\n\t"\
			"movaps	0x10(%%eax),%%xmm1		/* a(jp   ) */\n\t"\
			"movaps	0x10(%%edx),%%xmm5		/* a(jp+p2) */\n\t"\
			"movaps	%%xmm0,%%xmm2			/* xmm2 <- cpy a(jt   ) */\n\t"\
			"movaps	%%xmm4,%%xmm6			/* xmm4 <- cpy a(jt+p2) */\n\t"\
			"movaps	%%xmm1,%%xmm3			/* xmm3 <- cpy a(jp   ) */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* xmm5 <- cpy a(jp+p2) */\n\t"\
			"addpd	    (%%ebx),%%xmm0		/* t1 */\n\t"\
			"addpd	    (%%ecx),%%xmm4		/* t5 */\n\t"\
			"addpd	0x10(%%ebx),%%xmm1		/* t2 */\n\t"\
			"addpd	0x10(%%ecx),%%xmm5		/* t6 */\n\t"\
			"subpd	    (%%ebx),%%xmm2		/* t3 */\n\t"\
			"subpd	    (%%ecx),%%xmm6		/* t7 */\n\t"\
			"subpd	0x10(%%ebx),%%xmm3		/* t4 */\n\t"\
			"subpd	0x10(%%ecx),%%xmm7		/* t8 */\n\t"\
			"/* Finish radix-4 butterfly and store results into output-array slots: */\n\t"\
			"/* Swap _o1 <-> _o3: 0x20,0x30 <-> 0x60,0x70: */\n\t"\
			"subpd	%%xmm4,%%xmm0			/* ~t5 <- t1 -t5 */\n\t"\
			"subpd	%%xmm7,%%xmm2			/* ~t7 <- t3 -t8 */\n\t"\
			"subpd	%%xmm5,%%xmm1			/* ~t6 <- t2 -t6 */\n\t"\
			"subpd	%%xmm6,%%xmm3			/* ~t4 <- t4 -t7 */\n\t"\
			"movaps	%%xmm0,0x040(%%esi)		/* <- ~t5 */\n\t"\
			"movaps	%%xmm2,0x020(%%esi)		/* <- ~t7 */\n\t"\
			"movaps	%%xmm1,0x050(%%esi)		/* <- ~t6 */\n\t"\
			"movaps	%%xmm3,0x070(%%esi)		/* <- ~t4 */\n\t"\
			"addpd	%%xmm4,%%xmm4			/*          2*t5 */\n\t"\
			"addpd	%%xmm7,%%xmm7			/*          2*t8 */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*          2*t6 */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*          2*t7 */\n\t"\
			"addpd	%%xmm0,%%xmm4			/* ~t1 <- t1 +t5 */\n\t"\
			"addpd	%%xmm2,%%xmm7			/* ~t3 <- t3 +t8 */\n\t"\
			"addpd	%%xmm1,%%xmm5			/* ~t2 <- t2 +t6 */\n\t"\
			"addpd	%%xmm3,%%xmm6			/* ~t8 <- t4 +t7 */\n\t"\
			"movaps	%%xmm4,     (%%esi)		/* <- ~t1 */\n\t"\
			"movaps	%%xmm7,0x060(%%esi)		/* <- ~t3 */\n\t"\
			"movaps	%%xmm5,0x010(%%esi)		/* <- ~t2 */\n\t"\
			"movaps	%%xmm6,0x030(%%esi)		/* <- ~t8 */\n\t"\
		"/* SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_E(add12+p[3,2,1,0], s1p04r,s1p07r,s1p06r,s1p05r): */\n\t"\
			"addl	$0x80,%%esi	/* s1p04r */\n\t"\
			"subl	%%eax,%%ebx		/* p01 */\n\t"\
			"subl	%%eax,%%ecx		/* p02 */\n\t"\
			"subl	%%eax,%%edx		/* p03 */\n\t"\
			"movl	%[__p12],%%edi\n\t"\
			"shll	$3,%%edi		/* Pointer offset for floating doubles */\n\t"\
			"addl	%%edi,%%eax		/* add12 */\n\t"\
			"addl	%%eax,%%ebx\n\t"\
			"addl	%%eax,%%ecx\n\t"\
			"addl	%%eax,%%edx\n\t"\
			"/* eax <-> edx, ebx <-> ecx */\n\t"\
			"movaps	    (%%edx),%%xmm0		/* a(jt   ) */\n\t"\
			"movaps	    (%%ebx),%%xmm4		/* a(jt+p2) */\n\t"\
			"movaps	0x10(%%edx),%%xmm1		/* a(jp   ) */\n\t"\
			"movaps	0x10(%%ebx),%%xmm5		/* a(jp+p2) */\n\t"\
			"movaps	%%xmm0,%%xmm2			/* xmm2 <- cpy a(jt   ) */\n\t"\
			"movaps	%%xmm4,%%xmm6			/* xmm4 <- cpy a(jt+p2) */\n\t"\
			"movaps	%%xmm1,%%xmm3			/* xmm3 <- cpy a(jp   ) */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* xmm5 <- cpy a(jp+p2) */\n\t"\
			"addpd	    (%%ecx),%%xmm0		/* t1 */\n\t"\
			"addpd	    (%%eax),%%xmm4		/* t5 */\n\t"\
			"addpd	0x10(%%ecx),%%xmm1		/* t2 */\n\t"\
			"addpd	0x10(%%eax),%%xmm5		/* t6 */\n\t"\
			"subpd	    (%%ecx),%%xmm2		/* t3 */\n\t"\
			"subpd	    (%%eax),%%xmm6		/* t7 */\n\t"\
			"subpd	0x10(%%ecx),%%xmm3		/* t4 */\n\t"\
			"subpd	0x10(%%eax),%%xmm7		/* t8 */\n\t"\
			"/* Finish radix-4 butterfly and store results into output-array slots: */\n\t"\
			"/* Swap _o1 <-> _o3: 0x20,0x30 <-> 0x60,0x70: */\n\t"\
			"subpd	%%xmm4,%%xmm0			/* ~t5 <- t1 -t5 */\n\t"\
			"subpd	%%xmm7,%%xmm2			/* ~t7 <- t3 -t8 */\n\t"\
			"subpd	%%xmm5,%%xmm1			/* ~t6 <- t2 -t6 */\n\t"\
			"subpd	%%xmm6,%%xmm3			/* ~t4 <- t4 -t7 */\n\t"\
			"movaps	%%xmm0,0x040(%%esi)		/* <- ~t5 */\n\t"\
			"movaps	%%xmm2,0x020(%%esi)		/* <- ~t7 */\n\t"\
			"movaps	%%xmm1,0x050(%%esi)		/* <- ~t6 */\n\t"\
			"movaps	%%xmm3,0x070(%%esi)		/* <- ~t4 */\n\t"\
			"addpd	%%xmm4,%%xmm4			/*          2*t5 */\n\t"\
			"addpd	%%xmm7,%%xmm7			/*          2*t8 */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*          2*t6 */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*          2*t7 */\n\t"\
			"addpd	%%xmm0,%%xmm4			/* ~t1 <- t1 +t5 */\n\t"\
			"addpd	%%xmm2,%%xmm7			/* ~t3 <- t3 +t8 */\n\t"\
			"addpd	%%xmm1,%%xmm5			/* ~t2 <- t2 +t6 */\n\t"\
			"addpd	%%xmm3,%%xmm6			/* ~t8 <- t4 +t7 */\n\t"\
			"movaps	%%xmm4,     (%%esi)		/* <- ~t1 */\n\t"\
			"movaps	%%xmm7,0x060(%%esi)		/* <- ~t3 */\n\t"\
			"movaps	%%xmm5,0x010(%%esi)		/* <- ~t2 */\n\t"\
			"movaps	%%xmm6,0x030(%%esi)		/* <- ~t8 */\n\t"\
		"/* SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_E(add24+p[1,0,2,3], s1p08r,s1p11r,s1p10r,s1p09r): */\n\t"\
			"addl	$0x80,%%esi	/* s1p08r */\n\t"\
			"subl	%%eax,%%ebx		/* p01 */\n\t"\
			"subl	%%eax,%%ecx		/* p02 */\n\t"\
			"subl	%%eax,%%edx		/* p03 */\n\t"\
			"subl	%%edi,%%eax		/* add0 */\n\t"\
			"movl	%[__p24],%%edi\n\t"\
			"shll	$3,%%edi		/* Pointer offset for floating doubles */\n\t"\
			"addl	%%edi,%%eax		/* add24 */\n\t"\
			"addl	%%eax,%%ebx\n\t"\
			"addl	%%eax,%%ecx\n\t"\
			"addl	%%eax,%%edx\n\t"\
			"/* eax <-> ebx */\n\t"\
			"movaps	    (%%ebx),%%xmm0		/* a(jt   ) */\n\t"\
			"movaps	    (%%ecx),%%xmm4		/* a(jt+p2) */\n\t"\
			"movaps	0x10(%%ebx),%%xmm1		/* a(jp   ) */\n\t"\
			"movaps	0x10(%%ecx),%%xmm5		/* a(jp+p2) */\n\t"\
			"movaps	%%xmm0,%%xmm2			/* xmm2 <- cpy a(jt   ) */\n\t"\
			"movaps	%%xmm4,%%xmm6			/* xmm4 <- cpy a(jt+p2) */\n\t"\
			"movaps	%%xmm1,%%xmm3			/* xmm3 <- cpy a(jp   ) */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* xmm5 <- cpy a(jp+p2) */\n\t"\
			"addpd	    (%%eax),%%xmm0		/* t1 */\n\t"\
			"addpd	    (%%edx),%%xmm4		/* t5 */\n\t"\
			"addpd	0x10(%%eax),%%xmm1		/* t2 */\n\t"\
			"addpd	0x10(%%edx),%%xmm5		/* t6 */\n\t"\
			"subpd	    (%%eax),%%xmm2		/* t3 */\n\t"\
			"subpd	    (%%edx),%%xmm6		/* t7 */\n\t"\
			"subpd	0x10(%%eax),%%xmm3		/* t4 */\n\t"\
			"subpd	0x10(%%edx),%%xmm7		/* t8 */\n\t"\
			"/* Finish radix-4 butterfly and store results into output-array slots: */\n\t"\
			"/* Swap _o1 <-> _o3: 0x20,0x30 <-> 0x60,0x70: */\n\t"\
			"subpd	%%xmm4,%%xmm0			/* ~t5 <- t1 -t5 */\n\t"\
			"subpd	%%xmm7,%%xmm2			/* ~t7 <- t3 -t8 */\n\t"\
			"subpd	%%xmm5,%%xmm1			/* ~t6 <- t2 -t6 */\n\t"\
			"subpd	%%xmm6,%%xmm3			/* ~t4 <- t4 -t7 */\n\t"\
			"movaps	%%xmm0,0x040(%%esi)		/* <- ~t5 */\n\t"\
			"movaps	%%xmm2,0x020(%%esi)		/* <- ~t7 */\n\t"\
			"movaps	%%xmm1,0x050(%%esi)		/* <- ~t6 */\n\t"\
			"movaps	%%xmm3,0x070(%%esi)		/* <- ~t4 */\n\t"\
			"addpd	%%xmm4,%%xmm4			/*          2*t5 */\n\t"\
			"addpd	%%xmm7,%%xmm7			/*          2*t8 */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*          2*t6 */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*          2*t7 */\n\t"\
			"addpd	%%xmm0,%%xmm4			/* ~t1 <- t1 +t5 */\n\t"\
			"addpd	%%xmm2,%%xmm7			/* ~t3 <- t3 +t8 */\n\t"\
			"addpd	%%xmm1,%%xmm5			/* ~t2 <- t2 +t6 */\n\t"\
			"addpd	%%xmm3,%%xmm6			/* ~t8 <- t4 +t7 */\n\t"\
			"movaps	%%xmm4,     (%%esi)		/* <- ~t1 */\n\t"\
			"movaps	%%xmm7,0x060(%%esi)		/* <- ~t3 */\n\t"\
			"movaps	%%xmm5,0x010(%%esi)		/* <- ~t2 */\n\t"\
			"movaps	%%xmm6,0x030(%%esi)		/* <- ~t8 */\n\t"\
		"/* SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_E(add08+p[1,0,2,3], s1p12r,s1p15r,s1p14r,s1p13r): */\n\t"\
			"addl	$0x80,%%esi	/* s1p12r */\n\t"\
			"subl	%%eax,%%ebx		/* p01 */\n\t"\
			"subl	%%eax,%%ecx		/* p02 */\n\t"\
			"subl	%%eax,%%edx		/* p03 */\n\t"\
			"subl	%%edi,%%eax		/* add0 */\n\t"\
			"movl	%[__p08],%%edi\n\t"\
			"shll	$3,%%edi		/* Pointer offset for floating doubles */\n\t"\
			"addl	%%edi,%%eax		/* add08 */\n\t"\
			"addl	%%eax,%%ebx\n\t"\
			"addl	%%eax,%%ecx\n\t"\
			"addl	%%eax,%%edx\n\t"\
			"/* eax <-> %%ebx */\n\t"\
			"movaps	    (%%ebx),%%xmm0		/* a(jt   ) */\n\t"\
			"movaps	    (%%ecx),%%xmm4		/* a(jt+p2) */\n\t"\
			"movaps	0x10(%%ebx),%%xmm1		/* a(jp   ) */\n\t"\
			"movaps	0x10(%%ecx),%%xmm5		/* a(jp+p2) */\n\t"\
			"movaps	%%xmm0,%%xmm2			/* xmm2 <- cpy a(jt   ) */\n\t"\
			"movaps	%%xmm4,%%xmm6			/* xmm4 <- cpy a(jt+p2) */\n\t"\
			"movaps	%%xmm1,%%xmm3			/* xmm3 <- cpy a(jp   ) */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* xmm5 <- cpy a(jp+p2) */\n\t"\
			"addpd	    (%%eax),%%xmm0		/* t1 */\n\t"\
			"addpd	    (%%edx),%%xmm4		/* t5 */\n\t"\
			"addpd	0x10(%%eax),%%xmm1		/* t2 */\n\t"\
			"addpd	0x10(%%edx),%%xmm5		/* t6 */\n\t"\
			"subpd	    (%%eax),%%xmm2		/* t3 */\n\t"\
			"subpd	    (%%edx),%%xmm6		/* t7 */\n\t"\
			"subpd	0x10(%%eax),%%xmm3		/* t4 */\n\t"\
			"subpd	0x10(%%edx),%%xmm7		/* t8 */\n\t"\
			"/* Finish radix-4 butterfly and store results into output-array slots: */\n\t"\
			"/* Swap _o1 <-> _o3: 0x20,0x30 <-> 0x60,0x70: */\n\t"\
			"subpd	%%xmm4,%%xmm0			/* ~t5 <- t1 -t5 */\n\t"\
			"subpd	%%xmm7,%%xmm2			/* ~t7 <- t3 -t8 */\n\t"\
			"subpd	%%xmm5,%%xmm1			/* ~t6 <- t2 -t6 */\n\t"\
			"subpd	%%xmm6,%%xmm3			/* ~t4 <- t4 -t7 */\n\t"\
			"movaps	%%xmm0,0x040(%%esi)		/* <- ~t5 */\n\t"\
			"movaps	%%xmm2,0x020(%%esi)		/* <- ~t7 */\n\t"\
			"movaps	%%xmm1,0x050(%%esi)		/* <- ~t6 */\n\t"\
			"movaps	%%xmm3,0x070(%%esi)		/* <- ~t4 */\n\t"\
			"addpd	%%xmm4,%%xmm4			/*          2*t5 */\n\t"\
			"addpd	%%xmm7,%%xmm7			/*          2*t8 */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*          2*t6 */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*          2*t7 */\n\t"\
			"addpd	%%xmm0,%%xmm4			/* ~t1 <- t1 +t5 */\n\t"\
			"addpd	%%xmm2,%%xmm7			/* ~t3 <- t3 +t8 */\n\t"\
			"addpd	%%xmm1,%%xmm5			/* ~t2 <- t2 +t6 */\n\t"\
			"addpd	%%xmm3,%%xmm6			/* ~t8 <- t4 +t7 */\n\t"\
			"movaps	%%xmm4,     (%%esi)		/* <- ~t1 */\n\t"\
			"movaps	%%xmm7,0x060(%%esi)		/* <- ~t3 */\n\t"\
			"movaps	%%xmm5,0x010(%%esi)		/* <- ~t2 */\n\t"\
			"movaps	%%xmm6,0x030(%%esi)		/* <- ~t8 */\n\t"\
		"/* SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_E(add20+p[2,3,0,1], s1p16r,s1p19r,s1p18r,s1p17r): */\n\t"\
			"addl	$0x80,%%esi	/* s1p16r */\n\t"\
			"subl	%%eax,%%ebx		/* p01 */\n\t"\
			"subl	%%eax,%%ecx		/* p02 */\n\t"\
			"subl	%%eax,%%edx		/* p03 */\n\t"\
			"subl	%%edi,%%eax		/* add0 */\n\t"\
			"movl	%[__p20],%%edi\n\t"\
			"shll	$3,%%edi		/* Pointer offset for floating doubles */\n\t"\
			"addl	%%edi,%%eax		/* add20 */\n\t"\
			"addl	%%eax,%%ebx\n\t"\
			"addl	%%eax,%%ecx\n\t"\
			"addl	%%eax,%%edx\n\t"\
			"/* eax <-> %%ecx, %%ebx <-> edx */\n\t"\
			"movaps	    (%%ecx),%%xmm0		/* a(jt   ) */\n\t"\
			"movaps	    (%%eax),%%xmm4		/* a(jt+p2) */\n\t"\
			"movaps	0x10(%%ecx),%%xmm1		/* a(jp   ) */\n\t"\
			"movaps	0x10(%%eax),%%xmm5		/* a(jp+p2) */\n\t"\
			"movaps	%%xmm0,%%xmm2			/* xmm2 <- cpy a(jt   ) */\n\t"\
			"movaps	%%xmm4,%%xmm6			/* xmm4 <- cpy a(jt+p2) */\n\t"\
			"movaps	%%xmm1,%%xmm3			/* xmm3 <- cpy a(jp   ) */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* xmm5 <- cpy a(jp+p2) */\n\t"\
			"addpd	    (%%edx),%%xmm0		/* t1 */\n\t"\
			"addpd	    (%%ebx),%%xmm4		/* t5 */\n\t"\
			"addpd	0x10(%%edx),%%xmm1		/* t2 */\n\t"\
			"addpd	0x10(%%ebx),%%xmm5		/* t6 */\n\t"\
			"subpd	    (%%edx),%%xmm2		/* t3 */\n\t"\
			"subpd	    (%%ebx),%%xmm6		/* t7 */\n\t"\
			"subpd	0x10(%%edx),%%xmm3		/* t4 */\n\t"\
			"subpd	0x10(%%ebx),%%xmm7		/* t8 */\n\t"\
			"/* Finish radix-4 butterfly and store results into output-array slots: */\n\t"\
			"/* Swap _o1 <-> _o3: 0x20,0x30 <-> 0x60,0x70: */\n\t"\
			"subpd	%%xmm4,%%xmm0			/* ~t5 <- t1 -t5 */\n\t"\
			"subpd	%%xmm7,%%xmm2			/* ~t7 <- t3 -t8 */\n\t"\
			"subpd	%%xmm5,%%xmm1			/* ~t6 <- t2 -t6 */\n\t"\
			"subpd	%%xmm6,%%xmm3			/* ~t4 <- t4 -t7 */\n\t"\
			"movaps	%%xmm0,0x040(%%esi)		/* <- ~t5 */\n\t"\
			"movaps	%%xmm2,0x020(%%esi)		/* <- ~t7 */\n\t"\
			"movaps	%%xmm1,0x050(%%esi)		/* <- ~t6 */\n\t"\
			"movaps	%%xmm3,0x070(%%esi)		/* <- ~t4 */\n\t"\
			"addpd	%%xmm4,%%xmm4			/*          2*t5 */\n\t"\
			"addpd	%%xmm7,%%xmm7			/*          2*t8 */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*          2*t6 */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*          2*t7 */\n\t"\
			"addpd	%%xmm0,%%xmm4			/* ~t1 <- t1 +t5 */\n\t"\
			"addpd	%%xmm2,%%xmm7			/* ~t3 <- t3 +t8 */\n\t"\
			"addpd	%%xmm1,%%xmm5			/* ~t2 <- t2 +t6 */\n\t"\
			"addpd	%%xmm3,%%xmm6			/* ~t8 <- t4 +t7 */\n\t"\
			"movaps	%%xmm4,     (%%esi)		/* <- ~t1 */\n\t"\
			"movaps	%%xmm7,0x060(%%esi)		/* <- ~t3 */\n\t"\
			"movaps	%%xmm5,0x010(%%esi)		/* <- ~t2 */\n\t"\
			"movaps	%%xmm6,0x030(%%esi)		/* <- ~t8 */\n\t"\
		"/* SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_E(add04+p[2,3,0,1], s1p20r,s1p23r,s1p22r,s1p21r): */\n\t"\
			"addl	$0x80,%%esi	/* s1p20r */\n\t"\
			"subl	%%eax,%%ebx		/* p01 */\n\t"\
			"subl	%%eax,%%ecx		/* p02 */\n\t"\
			"subl	%%eax,%%edx		/* p03 */\n\t"\
			"subl	%%edi,%%eax		/* add0 */\n\t"\
			"movl	%[__p04],%%edi\n\t"\
			"shll	$3,%%edi		/* Pointer offset for floating doubles */\n\t"\
			"addl	%%edi,%%eax		/* add04 */\n\t"\
			"addl	%%eax,%%ebx\n\t"\
			"addl	%%eax,%%ecx\n\t"\
			"addl	%%eax,%%edx\n\t"\
			"/* eax <-> %%ecx, %%ebx <-> edx */\n\t"\
			"movaps	    (%%ecx),%%xmm0		/* a(jt   ) */\n\t"\
			"movaps	    (%%eax),%%xmm4		/* a(jt+p2) */\n\t"\
			"movaps	0x10(%%ecx),%%xmm1		/* a(jp   ) */\n\t"\
			"movaps	0x10(%%eax),%%xmm5		/* a(jp+p2) */\n\t"\
			"movaps	%%xmm0,%%xmm2			/* xmm2 <- cpy a(jt   ) */\n\t"\
			"movaps	%%xmm4,%%xmm6			/* xmm4 <- cpy a(jt+p2) */\n\t"\
			"movaps	%%xmm1,%%xmm3			/* xmm3 <- cpy a(jp   ) */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* xmm5 <- cpy a(jp+p2) */\n\t"\
			"addpd	    (%%edx),%%xmm0		/* t1 */\n\t"\
			"addpd	    (%%ebx),%%xmm4		/* t5 */\n\t"\
			"addpd	0x10(%%edx),%%xmm1		/* t2 */\n\t"\
			"addpd	0x10(%%ebx),%%xmm5		/* t6 */\n\t"\
			"subpd	    (%%edx),%%xmm2		/* t3 */\n\t"\
			"subpd	    (%%ebx),%%xmm6		/* t7 */\n\t"\
			"subpd	0x10(%%edx),%%xmm3		/* t4 */\n\t"\
			"subpd	0x10(%%ebx),%%xmm7		/* t8 */\n\t"\
			"/* Finish radix-4 butterfly and store results into output-array slots: */\n\t"\
			"/* Swap _o1 <-> _o3: 0x20,0x30 <-> 0x60,0x70: */\n\t"\
			"subpd	%%xmm4,%%xmm0			/* ~t5 <- t1 -t5 */\n\t"\
			"subpd	%%xmm7,%%xmm2			/* ~t7 <- t3 -t8 */\n\t"\
			"subpd	%%xmm5,%%xmm1			/* ~t6 <- t2 -t6 */\n\t"\
			"subpd	%%xmm6,%%xmm3			/* ~t4 <- t4 -t7 */\n\t"\
			"movaps	%%xmm0,0x040(%%esi)		/* <- ~t5 */\n\t"\
			"movaps	%%xmm2,0x020(%%esi)		/* <- ~t7 */\n\t"\
			"movaps	%%xmm1,0x050(%%esi)		/* <- ~t6 */\n\t"\
			"movaps	%%xmm3,0x070(%%esi)		/* <- ~t4 */\n\t"\
			"addpd	%%xmm4,%%xmm4			/*          2*t5 */\n\t"\
			"addpd	%%xmm7,%%xmm7			/*          2*t8 */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*          2*t6 */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*          2*t7 */\n\t"\
			"addpd	%%xmm0,%%xmm4			/* ~t1 <- t1 +t5 */\n\t"\
			"addpd	%%xmm2,%%xmm7			/* ~t3 <- t3 +t8 */\n\t"\
			"addpd	%%xmm1,%%xmm5			/* ~t2 <- t2 +t6 */\n\t"\
			"addpd	%%xmm3,%%xmm6			/* ~t8 <- t4 +t7 */\n\t"\
			"movaps	%%xmm4,     (%%esi)		/* <- ~t1 */\n\t"\
			"movaps	%%xmm7,0x060(%%esi)		/* <- ~t3 */\n\t"\
			"movaps	%%xmm5,0x010(%%esi)		/* <- ~t2 */\n\t"\
			"movaps	%%xmm6,0x030(%%esi)		/* <- ~t8 */\n\t"\
		"/* SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_E(add16+p[0,1,3,2], s1p24r,s1p27r,s1p26r,s1p25r): */\n\t"\
			"addl	$0x80,%%esi	/* s1p24r */\n\t"\
			"subl	%%eax,%%ebx		/* p01 */\n\t"\
			"subl	%%eax,%%ecx		/* p02 */\n\t"\
			"subl	%%eax,%%edx		/* p03 */\n\t"\
			"subl	%%edi,%%eax		/* add0 */\n\t"\
			"movl	%[__p16],%%edi\n\t"\
			"shll	$3,%%edi		/* Pointer offset for floating doubles */\n\t"\
			"addl	%%edi,%%eax		/* add16 */\n\t"\
			"addl	%%eax,%%ebx\n\t"\
			"addl	%%eax,%%ecx\n\t"\
			"addl	%%eax,%%edx\n\t"\
			"/* ecx <-> %%edx */\n\t"\
			"movaps	    (%%eax),%%xmm0		/* a(jt   ) */\n\t"\
			"movaps	    (%%edx),%%xmm4		/* a(jt+p2) */\n\t"\
			"movaps	0x10(%%eax),%%xmm1		/* a(jp   ) */\n\t"\
			"movaps	0x10(%%edx),%%xmm5		/* a(jp+p2) */\n\t"\
			"movaps	%%xmm0,%%xmm2			/* xmm2 <- cpy a(jt   ) */\n\t"\
			"movaps	%%xmm4,%%xmm6			/* xmm4 <- cpy a(jt+p2) */\n\t"\
			"movaps	%%xmm1,%%xmm3			/* xmm3 <- cpy a(jp   ) */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* xmm5 <- cpy a(jp+p2) */\n\t"\
			"addpd	    (%%ebx),%%xmm0		/* t1 */\n\t"\
			"addpd	    (%%ecx),%%xmm4		/* t5 */\n\t"\
			"addpd	0x10(%%ebx),%%xmm1		/* t2 */\n\t"\
			"addpd	0x10(%%ecx),%%xmm5		/* t6 */\n\t"\
			"subpd	    (%%ebx),%%xmm2		/* t3 */\n\t"\
			"subpd	    (%%ecx),%%xmm6		/* t7 */\n\t"\
			"subpd	0x10(%%ebx),%%xmm3		/* t4 */\n\t"\
			"subpd	0x10(%%ecx),%%xmm7		/* t8 */\n\t"\
			"/* Finish radix-4 butterfly and store results into output-array slots: */\n\t"\
			"/* Swap _o1 <-> _o3: 0x20,0x30 <-> 0x60,0x70: */\n\t"\
			"subpd	%%xmm4,%%xmm0			/* ~t5 <- t1 -t5 */\n\t"\
			"subpd	%%xmm7,%%xmm2			/* ~t7 <- t3 -t8 */\n\t"\
			"subpd	%%xmm5,%%xmm1			/* ~t6 <- t2 -t6 */\n\t"\
			"subpd	%%xmm6,%%xmm3			/* ~t4 <- t4 -t7 */\n\t"\
			"movaps	%%xmm0,0x040(%%esi)		/* <- ~t5 */\n\t"\
			"movaps	%%xmm2,0x020(%%esi)		/* <- ~t7 */\n\t"\
			"movaps	%%xmm1,0x050(%%esi)		/* <- ~t6 */\n\t"\
			"movaps	%%xmm3,0x070(%%esi)		/* <- ~t4 */\n\t"\
			"addpd	%%xmm4,%%xmm4			/*          2*t5 */\n\t"\
			"addpd	%%xmm7,%%xmm7			/*          2*t8 */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*          2*t6 */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*          2*t7 */\n\t"\
			"addpd	%%xmm0,%%xmm4			/* ~t1 <- t1 +t5 */\n\t"\
			"addpd	%%xmm2,%%xmm7			/* ~t3 <- t3 +t8 */\n\t"\
			"addpd	%%xmm1,%%xmm5			/* ~t2 <- t2 +t6 */\n\t"\
			"addpd	%%xmm3,%%xmm6			/* ~t8 <- t4 +t7 */\n\t"\
			"movaps	%%xmm4,     (%%esi)		/* <- ~t1 */\n\t"\
			"movaps	%%xmm7,0x060(%%esi)		/* <- ~t3 */\n\t"\
			"movaps	%%xmm5,0x010(%%esi)		/* <- ~t2 */\n\t"\
			"movaps	%%xmm6,0x030(%%esi)		/* <- ~t8 */\n\t"\
			"\n\t"\
			"/*...and now do 4 radix-7 transforms...*/\n\t"\
			"\n\t"\
		"/* SSE2_RADIX_07_DFT(s1p00r,s1p04r,s1p08r,s1p12r,s1p16r,s1p20r,s1p24r,cc0,s1p00r,s1p08r,s1p16r,s1p24r,s1p04r,s1p12r,s1p20r): */\n\t"\
			"/*\n\t"\
			"t1r=A1r+A6r;	t2r=A2r+A5r;	t3r=A3r+A4r;\n\t"\
			"t6r=A1r-A6r;	t5r=A2r-A5r;	t4r=A3r-A4r;\n\t"\
			"*/\n\t"\
			"movl	%[__out],%%edi	/* s1p00r */\n\t"\
			"movl	%[__cc0],%%eax	/* cc0 */\n\t"\
			"movl	$0x20,%%ebx\n\t"\
			"movl	$0x40,%%ecx\n\t"\
			"movl	$0x60,%%edx\n\t"\
			"movl	$0x80,%%esi\n\t"\
			"addl	%%eax,%%ebx		/* cc1 */\n\t"\
			"addl	%%eax,%%ecx		/* cc2 */\n\t"\
			"addl	%%eax,%%edx		/* cc3 */\n\t"\
			"addl	%%eax,%%esi		/* scratch storage */\n\t"\
			"movaps	0x080(%%edi),%%xmm6			/* A1r = s1p04r */\n\t"\
			"movaps	0x300(%%edi),%%xmm1			/* A6r = s1p24r */\n\t"\
			"movaps	0x100(%%edi),%%xmm5			/* A2r = s1p08r */\n\t"\
			"movaps	0x280(%%edi),%%xmm2			/* A5r = s1p20r */\n\t"\
			"movaps	0x180(%%edi),%%xmm4			/* A3r = s1p12r */\n\t"\
			"movaps	0x200(%%edi),%%xmm3			/* A4r = s1p16r */\n\t"\
			"\n\t"\
			"subpd	%%xmm1,%%xmm6			/* t6r = A1r-A6r */\n\t"\
			"addpd	%%xmm1,%%xmm1			/*         2*A6r */\n\t"\
			"addpd	%%xmm6,%%xmm1			/* t1r = A1r+A6r */\n\t"\
			"\n\t"\
			"subpd	%%xmm2,%%xmm5			/* t5r = A2r-A5r */\n\t"\
			"addpd	%%xmm2,%%xmm2			/*         2*A5r */\n\t"\
			"addpd	%%xmm5,%%xmm2			/* t2r = A2r+A5r */\n\t"\
			"\n\t"\
			"movaps	    (%%edi),%%xmm0		/* Ar0 */\n\t"\
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
			"movaps	%%xmm0,    (%%esi)		/* cpy t0 into scratch sincos slot */\n\t"\
			"movaps	%%xmm6,0x10(%%esi)		/* cpy t6 into scratch sincos slot */\n\t"\
			"addpd	%%xmm1,%%xmm0			/*~A0 = A0+t1 */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* cpy t5 */\n\t"\
			"addpd	%%xmm2,%%xmm3			/*~t3 = t3+t2 */\n\t"\
			"subpd	%%xmm4,%%xmm5			/*~t5 = t5-t4 */\n\t"\
			"subpd	%%xmm2,%%xmm1			/*~t1 = t1-t2 */\n\t"\
			"subpd	%%xmm7,%%xmm6			/*~t6 = t6-t5 */\n\t"\
			"addpd	%%xmm2,%%xmm2			/* 2*t2 */\n\t"\
			"addpd	%%xmm7,%%xmm4			/*~t5 = t4+t5 */\n\t"\
			"addpd	%%xmm3,%%xmm0			/* B0 */\n\t"\
			"addpd	0x10(%%esi),%%xmm5		/* t3 = (t5-t4)+t6 */\n\t"\
			"subpd	%%xmm2,%%xmm3			/*~t2 =  (t2+t3) - 2*t2 = t3-t2 */\n\t"\
			"movaps	%%xmm4,%%xmm7			/* cpy t5 */\n\t"\
			"movaps	%%xmm0,     (%%edi)		/* <-B0 = s1p00r, xmm0 FREE */\n\t"\
			"subpd	%%xmm6,%%xmm4			/* t4 = ~t5-~t6 */\n\t"\
			"movaps	%%xmm1,%%xmm2			/* cpy ~t1 */\n\t"\
			"subpd	    (%%esi),%%xmm0		/* r = B0 - t0 */\n\t"\
			"mulpd	0x10(%%eax),%%xmm5		/*~t3 = t3*sx0 */\n\t"\
			"addpd	%%xmm3,%%xmm2			/* ~t1+~t2 */\n\t"\
			"mulpd	    (%%ecx),%%xmm3		/* t2 = t2*cx2 */\n\t"\
			"mulpd	0x10(%%edx),%%xmm4		/*~t4 = t4*sx3 */\n\t"\
			"mulpd	    (%%ebx),%%xmm1		/* t1 = t1*cx1 */\n\t"\
			"mulpd	0x10(%%ebx),%%xmm6		/*~t6 = t6*sx1 */\n\t"\
			"mulpd	    (%%eax),%%xmm0		/* ~r = r*(cx0-1) */\n\t"\
			"mulpd	0x10(%%ecx),%%xmm7		/*~t5 = t5*sx2 */\n\t"\
			"mulpd	    (%%edx),%%xmm2		/* t3 =(t1+t2)*cx3 */\n\t"\
			"addpd	    (%%edi),%%xmm0		/* t0 =~r + B0 */\n\t"\
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
			"movaps	%%xmm0,0x100(%%edi)			/* B1 <- t0 = s1p08r */\n\t"\
			"movaps	%%xmm5,0x280(%%edi)			/* B6 <- t3 = s1p20r */\n\t"\
			"movaps	%%xmm2,0x200(%%edi)			/* B2 <- t1 = s1p16r */\n\t"\
			"movaps	%%xmm7,0x180(%%edi)			/* B5 <- t4 = s1p12r */\n\t"\
			"movaps	%%xmm3,0x300(%%edi)			/* B3 <- t2 = s1p24r */\n\t"\
			"movaps	%%xmm4,0x080(%%edi)			/* B4 <- t5 = s1p04r */\n\t"\
			"/************************** Imaginary Parts: ******************************************/\n\t"\
			"movaps	0x090(%%edi),%%xmm6		/* A1i = s1p04i */\n\t"\
			"movaps	0x310(%%edi),%%xmm1		/* A6i = s1p24i */\n\t"\
			"movaps	0x110(%%edi),%%xmm5		/* A2i = s1p08i */\n\t"\
			"movaps	0x290(%%edi),%%xmm2		/* A5i = s1p20i */\n\t"\
			"movaps	0x190(%%edi),%%xmm4		/* A3i = s1p12i */\n\t"\
			"movaps	0x210(%%edi),%%xmm3		/* A4i = s1p16i */\n\t"\
			"\n\t"\
			"subpd	%%xmm1,%%xmm6			/* t6i = A1i-A6i */\n\t"\
			"addpd	%%xmm1,%%xmm1			/*         2*A6i */\n\t"\
			"addpd	%%xmm6,%%xmm1			/* t1i = A1i+A6i */\n\t"\
			"\n\t"\
			"subpd	%%xmm2,%%xmm5			/* t5i = A2i-A5i */\n\t"\
			"addpd	%%xmm2,%%xmm2			/*         2*A5i */\n\t"\
			"addpd	%%xmm5,%%xmm2			/* t2i = A2i+A5i */\n\t"\
			"\n\t"\
			"movaps	0x10(%%edi),%%xmm0		/* Ai0 */\n\t"\
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
			"movaps	%%xmm0,    (%%esi)		/* cpy t0 into scratch sincos slot */\n\t"\
			"movaps	%%xmm6,0x10(%%esi)		/* cpy t6 into scratch sincos slot */\n\t"\
			"addpd	%%xmm1,%%xmm0			/*~A0 = A0+t1 */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* cpy t5 */\n\t"\
			"addpd	%%xmm2,%%xmm3			/*~t3 = t3+t2 */\n\t"\
			"subpd	%%xmm4,%%xmm5			/*~t5 = t5-t4 */\n\t"\
			"subpd	%%xmm2,%%xmm1			/*~t1 = t1-t2 */\n\t"\
			"subpd	%%xmm7,%%xmm6			/*~t6 = t6-t5 */\n\t"\
			"addpd	%%xmm2,%%xmm2			/* 2*t2 */\n\t"\
			"addpd	%%xmm7,%%xmm4			/*~t5 = t4+t5 */\n\t"\
			"addpd	%%xmm3,%%xmm0			/* B0 */\n\t"\
			"addpd	0x10(%%esi),%%xmm5		/* t3 = (t5-t4)+t6 */\n\t"\
			"subpd	%%xmm2,%%xmm3			/*~t2 =  (t2+t3) - 2*t2 = t3-t2 */\n\t"\
			"movaps	%%xmm4,%%xmm7			/* cpy t5 */\n\t"\
			"movaps	%%xmm0,0x010(%%edi)		/* <-B0 = s1p00i, xmm0 FREE */\n\t"\
			"subpd	%%xmm6,%%xmm4			/* t4 = ~t5-~t6 */\n\t"\
			"movaps	%%xmm1,%%xmm2			/* cpy ~t1 */\n\t"\
			"subpd	    (%%esi),%%xmm0		/* r = B0 - t0 */\n\t"\
			"mulpd	0x10(%%eax),%%xmm5		/*~t3 = t3*sx0 */\n\t"\
			"addpd	%%xmm3,%%xmm2			/* ~t1+~t2 */\n\t"\
			"mulpd	    (%%ecx),%%xmm3		/* t2 = t2*cx2 */\n\t"\
			"mulpd	0x10(%%edx),%%xmm4		/*~t4 = t4*sx3 */\n\t"\
			"mulpd	    (%%ebx),%%xmm1		/* t1 = t1*cx1 */\n\t"\
			"mulpd	0x10(%%ebx),%%xmm6		/*~t6 = t6*sx1 */\n\t"\
			"mulpd	    (%%eax),%%xmm0		/* ~r = r*(cx0-1) */\n\t"\
			"mulpd	0x10(%%ecx),%%xmm7		/*~t5 = t5*sx2 */\n\t"\
			"mulpd	    (%%edx),%%xmm2		/* t3 =(t1+t2)*cx3 */\n\t"\
			"addpd	0x10(%%edi),%%xmm0		/* t0 =~r + B0 */\n\t"\
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
			"movaps	0x100(%%edi),%%xmm1		/* t0r */\n\t"\
			"movaps	0x280(%%edi),%%xmm6		/* t3r */\n\t"\
			"subpd	%%xmm5,%%xmm1			/* B1r =t0r-t3i */\n\t"\
			"subpd	%%xmm6,%%xmm0			/* B6i =t0i-t3r */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*        2*t3i */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*        2*t3r */\n\t"\
			"addpd	%%xmm1,%%xmm5			/* B6r =t0r+t3i */\n\t"\
			"addpd	%%xmm0,%%xmm6			/* B1i =t0i+t3r */\n\t"\
			"movaps	%%xmm1,0x100(%%edi)		/* <-B1r = s1p08r */\n\t"\
			"movaps	%%xmm0,0x290(%%edi)		/* <-B6i = s1p20r */\n\t"\
			"movaps	%%xmm5,0x280(%%edi)		/* <-B6r = s1p20r */\n\t"\
			"movaps	%%xmm6,0x110(%%edi)		/* <-B1i = s1p08i */\n\t"\
			"\n\t"\
			"movaps	0x200(%%edi),%%xmm1		/* t1r */\n\t"\
			"movaps	0x180(%%edi),%%xmm6		/* t4r */\n\t"\
			"subpd	%%xmm7,%%xmm1			/* B2r =t1r-t4i */\n\t"\
			"subpd	%%xmm6,%%xmm2			/* B5i =t1i-t4r */\n\t"\
			"addpd	%%xmm7,%%xmm7			/*        2*t4i */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*        2*t4r */\n\t"\
			"addpd	%%xmm1,%%xmm7			/* B5r =t1r+t4i */\n\t"\
			"addpd	%%xmm2,%%xmm6			/* B2i =t1i+t4r */\n\t"\
			"movaps	%%xmm1,0x200(%%edi)		/* <-B2r = s1p16r */\n\t"\
			"movaps	%%xmm2,0x190(%%edi)		/* <-B5i = s1p12r */\n\t"\
			"movaps	%%xmm7,0x180(%%edi)		/* <-B5r = s1p12r */\n\t"\
			"movaps	%%xmm6,0x210(%%edi)		/* <-B2i = s1p16i */\n\t"\
			"\n\t"\
			"/* Note the order reversal on this pair of outputs: */\n\t"\
			"movaps	0x300(%%edi),%%xmm0		/* t2r */\n\t"\
			"movaps	0x080(%%edi),%%xmm5		/* t5r */\n\t"\
			"subpd	%%xmm4,%%xmm0			/* B4r =t2r-t5i */\n\t"\
			"subpd	%%xmm5,%%xmm3			/* B3i =t2i-t5r */\n\t"\
			"addpd	%%xmm4,%%xmm4			/*        2*t5i */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*        2*t5r */\n\t"\
			"addpd	%%xmm0,%%xmm4			/* B3r =t2r+t5i */\n\t"\
			"addpd	%%xmm3,%%xmm5			/* B4i =t2i+t5r */\n\t"\
			"movaps	%%xmm0,0x080(%%edi)		/* <-B4r = s1p04r */\n\t"\
			"movaps	%%xmm3,0x310(%%edi)		/* <-B3i = s1p24r */\n\t"\
			"movaps	%%xmm4,0x300(%%edi)		/* <-B3r = s1p24r */\n\t"\
			"movaps	%%xmm5,0x090(%%edi)		/* <-B4i = s1p04i */\n\t"\
			"\n\t"\
		"/* SSE2_RADIX_07_DFT(s1p01r,s1p05r,s1p09r,s1p13r,s1p17r,s1p21r,s1p25r,cc0,s1p21r,s1p01r,s1p09r,s1p17r,s1p25r,s1p05r,s1p13r): */\n\t"\
			"addl	$0x20,%%edi	/* s1p01r */\n\t"\
			"movaps	0x080(%%edi),%%xmm6			/* A1r = s1p05r */\n\t"\
			"movaps	0x300(%%edi),%%xmm1			/* A6r = s1p25r */\n\t"\
			"movaps	0x100(%%edi),%%xmm5			/* A2r = s1p09r */\n\t"\
			"movaps	0x280(%%edi),%%xmm2			/* A5r = s1p21r */\n\t"\
			"movaps	0x180(%%edi),%%xmm4			/* A3r = s1p13r */\n\t"\
			"movaps	0x200(%%edi),%%xmm3			/* A4r = s1p17r */\n\t"\
			"\n\t"\
			"subpd	%%xmm1,%%xmm6			/* t6r = A1r-A6r */\n\t"\
			"addpd	%%xmm1,%%xmm1			/*         2*A6r */\n\t"\
			"addpd	%%xmm6,%%xmm1			/* t1r = A1r+A6r */\n\t"\
			"\n\t"\
			"subpd	%%xmm2,%%xmm5			/* t5r = A2r-A5r */\n\t"\
			"addpd	%%xmm2,%%xmm2			/*         2*A5r */\n\t"\
			"addpd	%%xmm5,%%xmm2			/* t2r = A2r+A5r */\n\t"\
			"\n\t"\
			"movaps	    (%%edi),%%xmm0		/* Ar0 */\n\t"\
			"subpd	%%xmm3,%%xmm4			/* t4r = A3r-A4r */\n\t"\
			"addpd	%%xmm3,%%xmm3			/*         2*A4r */\n\t"\
			"addpd	%%xmm4,%%xmm3			/* t3r = A3r+A4r */\n\t"\
			"\n\t"\
			"movaps	%%xmm0,    (%%esi)		/* cpy t0 into scratch sincos slot */\n\t"\
			"movaps	%%xmm6,0x10(%%esi)		/* cpy t6 into scratch sincos slot */\n\t"\
			"addpd	%%xmm1,%%xmm0			/*~A0 = A0+t1 */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* cpy t5 */\n\t"\
			"addpd	%%xmm2,%%xmm3			/*~t3 = t3+t2 */\n\t"\
			"subpd	%%xmm4,%%xmm5			/*~t5 = t5-t4 */\n\t"\
			"subpd	%%xmm2,%%xmm1			/*~t1 = t1-t2 */\n\t"\
			"subpd	%%xmm7,%%xmm6			/*~t6 = t6-t5 */\n\t"\
			"addpd	%%xmm2,%%xmm2			/* 2*t2 */\n\t"\
			"addpd	%%xmm7,%%xmm4			/*~t5 = t4+t5 */\n\t"\
			"addpd	%%xmm3,%%xmm0			/* B0 */\n\t"\
			"addpd	0x10(%%esi),%%xmm5		/* t3 = (t5-t4)+t6 */\n\t"\
			"subpd	%%xmm2,%%xmm3			/*~t2 =  (t2+t3) - 2*t2 = t3-t2 */\n\t"\
			"movaps	%%xmm4,%%xmm7			/* cpy t5 */\n\t"\
			"movaps	%%xmm0,0x280(%%edi)		/* <-B0 = s1p21r, xmm0 FREE */\n\t"\
			"subpd	%%xmm6,%%xmm4			/* t4 = ~t5-~t6 */\n\t"\
			"movaps	%%xmm1,%%xmm2			/* cpy ~t1 */\n\t"\
			"subpd	    (%%esi),%%xmm0		/* r = B0 - t0 */\n\t"\
			"mulpd	0x10(%%eax),%%xmm5		/*~t3 = t3*sx0 */\n\t"\
			"addpd	%%xmm3,%%xmm2			/* ~t1+~t2 */\n\t"\
			"mulpd	    (%%ecx),%%xmm3		/* t2 = t2*cx2 */\n\t"\
			"mulpd	0x10(%%edx),%%xmm4		/*~t4 = t4*sx3 */\n\t"\
			"mulpd	    (%%ebx),%%xmm1		/* t1 = t1*cx1 */\n\t"\
			"mulpd	0x10(%%ebx),%%xmm6		/*~t6 = t6*sx1 */\n\t"\
			"mulpd	    (%%eax),%%xmm0		/* ~r = r*(cx0-1) */\n\t"\
			"mulpd	0x10(%%ecx),%%xmm7		/*~t5 = t5*sx2 */\n\t"\
			"mulpd	    (%%edx),%%xmm2		/* t3 =(t1+t2)*cx3 */\n\t"\
			"addpd	0x280(%%edi),%%xmm0		/* t0 =~r + B0 */\n\t"\
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
			"movaps	%%xmm0,     (%%edi)			/* B1 <- t0 = s1p01r */\n\t"\
			"movaps	%%xmm5,0x180(%%edi)			/* B6 <- t3 = s1p13r */\n\t"\
			"movaps	%%xmm2,0x100(%%edi)			/* B2 <- t1 = s1p09r */\n\t"\
			"movaps	%%xmm7,0x080(%%edi)			/* B5 <- t4 = s1p05r */\n\t"\
			"movaps	%%xmm3,0x200(%%edi)			/* B3 <- t2 = s1p17r */\n\t"\
			"movaps	%%xmm4,0x300(%%edi)			/* B4 <- t5 = s1p25r */\n\t"\
			"/************************** Imaginary Parts: ******************************************/\n\t"\
			"movaps	0x090(%%edi),%%xmm6		/* A1i = s1p05i */\n\t"\
			"movaps	0x310(%%edi),%%xmm1		/* A6i = s1p25i */\n\t"\
			"movaps	0x110(%%edi),%%xmm5		/* A2i = s1p09i */\n\t"\
			"movaps	0x290(%%edi),%%xmm2		/* A5i = s1p21i */\n\t"\
			"movaps	0x190(%%edi),%%xmm4		/* A3i = s1p13i */\n\t"\
			"movaps	0x210(%%edi),%%xmm3		/* A4i = s1p17i */\n\t"\
			"\n\t"\
			"subpd	%%xmm1,%%xmm6			/* t6i = A1i-A6i */\n\t"\
			"addpd	%%xmm1,%%xmm1			/*         2*A6i */\n\t"\
			"addpd	%%xmm6,%%xmm1			/* t1i = A1i+A6i */\n\t"\
			"\n\t"\
			"subpd	%%xmm2,%%xmm5			/* t5i = A2i-A5i */\n\t"\
			"addpd	%%xmm2,%%xmm2			/*         2*A5i */\n\t"\
			"addpd	%%xmm5,%%xmm2			/* t2i = A2i+A5i */\n\t"\
			"\n\t"\
			"movaps	0x10(%%edi),%%xmm0		/* Ai0 */\n\t"\
			"subpd	%%xmm3,%%xmm4			/* t4i = A3i-A4i */\n\t"\
			"addpd	%%xmm3,%%xmm3			/*         2*A4i */\n\t"\
			"addpd	%%xmm4,%%xmm3			/* t3i = A3i+A4i */\n\t"\
			"\n\t"\
			"movaps	%%xmm0,    (%%esi)		/* cpy t0 into scratch sincos slot */\n\t"\
			"movaps	%%xmm6,0x10(%%esi)		/* cpy t6 into scratch sincos slot */\n\t"\
			"addpd	%%xmm1,%%xmm0			/*~A0 = A0+t1 */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* cpy t5 */\n\t"\
			"addpd	%%xmm2,%%xmm3			/*~t3 = t3+t2 */\n\t"\
			"subpd	%%xmm4,%%xmm5			/*~t5 = t5-t4 */\n\t"\
			"subpd	%%xmm2,%%xmm1			/*~t1 = t1-t2 */\n\t"\
			"subpd	%%xmm7,%%xmm6			/*~t6 = t6-t5 */\n\t"\
			"addpd	%%xmm2,%%xmm2			/* 2*t2 */\n\t"\
			"addpd	%%xmm7,%%xmm4			/*~t5 = t4+t5 */\n\t"\
			"addpd	%%xmm3,%%xmm0			/* B0 */\n\t"\
			"addpd	0x10(%%esi),%%xmm5		/* t3 = (t5-t4)+t6 */\n\t"\
			"subpd	%%xmm2,%%xmm3			/*~t2 =  (t2+t3) - 2*t2 = t3-t2 */\n\t"\
			"movaps	%%xmm4,%%xmm7			/* cpy t5 */\n\t"\
			"movaps	%%xmm0,0x290(%%edi)		/* <-B0 = s1p21i, xmm0 FREE */\n\t"\
			"subpd	%%xmm6,%%xmm4			/* t4 = ~t5-~t6 */\n\t"\
			"movaps	%%xmm1,%%xmm2			/* cpy ~t1 */\n\t"\
			"subpd	    (%%esi),%%xmm0		/* r = B0 - t0 */\n\t"\
			"mulpd	0x10(%%eax),%%xmm5		/*~t3 = t3*sx0 */\n\t"\
			"addpd	%%xmm3,%%xmm2			/* ~t1+~t2 */\n\t"\
			"mulpd	    (%%ecx),%%xmm3		/* t2 = t2*cx2 */\n\t"\
			"mulpd	0x10(%%edx),%%xmm4		/*~t4 = t4*sx3 */\n\t"\
			"mulpd	    (%%ebx),%%xmm1		/* t1 = t1*cx1 */\n\t"\
			"mulpd	0x10(%%ebx),%%xmm6		/*~t6 = t6*sx1 */\n\t"\
			"mulpd	    (%%eax),%%xmm0		/* ~r = r*(cx0-1) */\n\t"\
			"mulpd	0x10(%%ecx),%%xmm7		/*~t5 = t5*sx2 */\n\t"\
			"mulpd	    (%%edx),%%xmm2		/* t3 =(t1+t2)*cx3 */\n\t"\
			"addpd	0x290(%%edi),%%xmm0		/* t0 =~r + B0 */\n\t"\
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
			"movaps	     (%%edi),%%xmm1		/* t0r */\n\t"\
			"movaps	0x180(%%edi),%%xmm6		/* t3r */\n\t"\
			"subpd	%%xmm5,%%xmm1			/* B1r =t0r-t3i */\n\t"\
			"subpd	%%xmm6,%%xmm0			/* B6i =t0i-t3r */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*        2*t3i */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*        2*t3r */\n\t"\
			"addpd	%%xmm1,%%xmm5			/* B6r =t0r+t3i */\n\t"\
			"addpd	%%xmm0,%%xmm6			/* B1i =t0i+t3r */\n\t"\
			"movaps	%%xmm1,     (%%edi)		/* <-B1r = s1p08r */\n\t"\
			"movaps	%%xmm0,0x190(%%edi)		/* <-B6i = s1p20r */\n\t"\
			"movaps	%%xmm5,0x180(%%edi)		/* <-B6r = s1p20r */\n\t"\
			"movaps	%%xmm6,0x010(%%edi)		/* <-B1i = s1p08i */\n\t"\
			"\n\t"\
			"movaps	0x100(%%edi),%%xmm1		/* t1r */\n\t"\
			"movaps	0x080(%%edi),%%xmm6		/* t4r */\n\t"\
			"subpd	%%xmm7,%%xmm1			/* B2r =t1r-t4i */\n\t"\
			"subpd	%%xmm6,%%xmm2			/* B5i =t1i-t4r */\n\t"\
			"addpd	%%xmm7,%%xmm7			/*        2*t4i */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*        2*t4r */\n\t"\
			"addpd	%%xmm1,%%xmm7			/* B5r =t1r+t4i */\n\t"\
			"addpd	%%xmm2,%%xmm6			/* B2i =t1i+t4r */\n\t"\
			"movaps	%%xmm1,0x100(%%edi)		/* <-B2r = s1p16r */\n\t"\
			"movaps	%%xmm2,0x090(%%edi)		/* <-B5i = s1p12r */\n\t"\
			"movaps	%%xmm7,0x080(%%edi)		/* <-B5r = s1p12r */\n\t"\
			"movaps	%%xmm6,0x110(%%edi)		/* <-B2i = s1p16i */\n\t"\
			"\n\t"\
			"/* Note the order reversal on this pair of outputs: */\n\t"\
			"movaps	0x200(%%edi),%%xmm0		/* t2r */\n\t"\
			"movaps	0x300(%%edi),%%xmm5		/* t5r */\n\t"\
			"subpd	%%xmm4,%%xmm0			/* B4r =t2r-t5i */\n\t"\
			"subpd	%%xmm5,%%xmm3			/* B3i =t2i-t5r */\n\t"\
			"addpd	%%xmm4,%%xmm4			/*        2*t5i */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*        2*t5r */\n\t"\
			"addpd	%%xmm0,%%xmm4			/* B3r =t2r+t5i */\n\t"\
			"addpd	%%xmm3,%%xmm5			/* B4i =t2i+t5r */\n\t"\
			"movaps	%%xmm0,0x300(%%edi)		/* <-B4r = s1p04r */\n\t"\
			"movaps	%%xmm3,0x210(%%edi)		/* <-B3i = s1p24r */\n\t"\
			"movaps	%%xmm4,0x200(%%edi)		/* <-B3r = s1p24r */\n\t"\
			"movaps	%%xmm5,0x310(%%edi)		/* <-B4i = s1p04i */\n\t"\
			"\n\t"\
		"/* SSE2_RADIX_07_DFT(s1p02r,s1p06r,s1p10r,s1p14r,s1p18r,s1p22r,s1p26r,cc0,s1p14r,s1p22r,s1p02r,s1p10r,s1p18r,s1p26r,s1p06r): */\n\t"\
			"addl	$0x20,%%edi	/* s1p02r */\n\t"\
			"movaps	0x080(%%edi),%%xmm6			/* A1r = s1p06r */\n\t"\
			"movaps	0x300(%%edi),%%xmm1			/* A6r = s1p26r */\n\t"\
			"movaps	0x100(%%edi),%%xmm5			/* A2r = s1p10r */\n\t"\
			"movaps	0x280(%%edi),%%xmm2			/* A5r = s1p22r */\n\t"\
			"movaps	0x180(%%edi),%%xmm4			/* A3r = s1p14r */\n\t"\
			"movaps	0x200(%%edi),%%xmm3			/* A4r = s1p18r */\n\t"\
			"\n\t"\
			"subpd	%%xmm1,%%xmm6			/* t6r = A1r-A6r */\n\t"\
			"addpd	%%xmm1,%%xmm1			/*         2*A6r */\n\t"\
			"addpd	%%xmm6,%%xmm1			/* t1r = A1r+A6r */\n\t"\
			"\n\t"\
			"subpd	%%xmm2,%%xmm5			/* t5r = A2r-A5r */\n\t"\
			"addpd	%%xmm2,%%xmm2			/*         2*A5r */\n\t"\
			"addpd	%%xmm5,%%xmm2			/* t2r = A2r+A5r */\n\t"\
			"\n\t"\
			"movaps	    (%%edi),%%xmm0		/* Ar0 */\n\t"\
			"subpd	%%xmm3,%%xmm4			/* t4r = A3r-A4r */\n\t"\
			"addpd	%%xmm3,%%xmm3			/*         2*A4r */\n\t"\
			"addpd	%%xmm4,%%xmm3			/* t3r = A3r+A4r */\n\t"\
			"\n\t"\
			"movaps	%%xmm0,    (%%esi)		/* cpy t0 into scratch sincos slot */\n\t"\
			"movaps	%%xmm6,0x10(%%esi)		/* cpy t6 into scratch sincos slot */\n\t"\
			"addpd	%%xmm1,%%xmm0			/*~A0 = A0+t1 */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* cpy t5 */\n\t"\
			"addpd	%%xmm2,%%xmm3			/*~t3 = t3+t2 */\n\t"\
			"subpd	%%xmm4,%%xmm5			/*~t5 = t5-t4 */\n\t"\
			"subpd	%%xmm2,%%xmm1			/*~t1 = t1-t2 */\n\t"\
			"subpd	%%xmm7,%%xmm6			/*~t6 = t6-t5 */\n\t"\
			"addpd	%%xmm2,%%xmm2			/* 2*t2 */\n\t"\
			"addpd	%%xmm7,%%xmm4			/*~t5 = t4+t5 */\n\t"\
			"addpd	%%xmm3,%%xmm0			/* B0 */\n\t"\
			"addpd	0x10(%%esi),%%xmm5		/* t3 = (t5-t4)+t6 */\n\t"\
			"subpd	%%xmm2,%%xmm3			/*~t2 =  (t2+t3) - 2*t2 = t3-t2 */\n\t"\
			"movaps	%%xmm4,%%xmm7			/* cpy t5 */\n\t"\
			"movaps	%%xmm0,0x180(%%edi)		/* <-B0 = s1p14r, xmm0 FREE */\n\t"\
			"subpd	%%xmm6,%%xmm4			/* t4 = ~t5-~t6 */\n\t"\
			"movaps	%%xmm1,%%xmm2			/* cpy ~t1 */\n\t"\
			"subpd	    (%%esi),%%xmm0		/* r = B0 - t0 */\n\t"\
			"mulpd	0x10(%%eax),%%xmm5		/*~t3 = t3*sx0 */\n\t"\
			"addpd	%%xmm3,%%xmm2			/* ~t1+~t2 */\n\t"\
			"mulpd	    (%%ecx),%%xmm3		/* t2 = t2*cx2 */\n\t"\
			"mulpd	0x10(%%edx),%%xmm4		/*~t4 = t4*sx3 */\n\t"\
			"mulpd	    (%%ebx),%%xmm1		/* t1 = t1*cx1 */\n\t"\
			"mulpd	0x10(%%ebx),%%xmm6		/*~t6 = t6*sx1 */\n\t"\
			"mulpd	    (%%eax),%%xmm0		/* ~r = r*(cx0-1) */\n\t"\
			"mulpd	0x10(%%ecx),%%xmm7		/*~t5 = t5*sx2 */\n\t"\
			"mulpd	    (%%edx),%%xmm2		/* t3 =(t1+t2)*cx3 */\n\t"\
			"addpd	0x180(%%edi),%%xmm0		/* t0 =~r + B0 */\n\t"\
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
			"movaps	%%xmm0,0x280(%%edi)			/* B1 <- t0 = s1p22r */\n\t"\
			"movaps	%%xmm5,0x080(%%edi)			/* B6 <- t3 = s1p06r */\n\t"\
			"movaps	%%xmm2,     (%%edi)			/* B2 <- t1 = s1p02r */\n\t"\
			"movaps	%%xmm7,0x300(%%edi)			/* B5 <- t4 = s1p26r */\n\t"\
			"movaps	%%xmm3,0x100(%%edi)			/* B3 <- t2 = s1p10r */\n\t"\
			"movaps	%%xmm4,0x200(%%edi)			/* B4 <- t5 = s1p18r */\n\t"\
			"/************************** Imaginary Parts: ******************************************/\n\t"\
			"movaps	0x090(%%edi),%%xmm6		/* A1i = s1p06i */\n\t"\
			"movaps	0x310(%%edi),%%xmm1		/* A6i = s1p26i */\n\t"\
			"movaps	0x110(%%edi),%%xmm5		/* A2i = s1p10i */\n\t"\
			"movaps	0x290(%%edi),%%xmm2		/* A5i = s1p22i */\n\t"\
			"movaps	0x190(%%edi),%%xmm4		/* A3i = s1p14i */\n\t"\
			"movaps	0x210(%%edi),%%xmm3		/* A4i = s1p18i */\n\t"\
			"\n\t"\
			"subpd	%%xmm1,%%xmm6			/* t6i = A1i-A6i */\n\t"\
			"addpd	%%xmm1,%%xmm1			/*         2*A6i */\n\t"\
			"addpd	%%xmm6,%%xmm1			/* t1i = A1i+A6i */\n\t"\
			"\n\t"\
			"subpd	%%xmm2,%%xmm5			/* t5i = A2i-A5i */\n\t"\
			"addpd	%%xmm2,%%xmm2			/*         2*A5i */\n\t"\
			"addpd	%%xmm5,%%xmm2			/* t2i = A2i+A5i */\n\t"\
			"\n\t"\
			"movaps	0x10(%%edi),%%xmm0		/* Ai0 */\n\t"\
			"subpd	%%xmm3,%%xmm4			/* t4i = A3i-A4i */\n\t"\
			"addpd	%%xmm3,%%xmm3			/*         2*A4i */\n\t"\
			"addpd	%%xmm4,%%xmm3			/* t3i = A3i+A4i */\n\t"\
			"\n\t"\
			"movaps	%%xmm0,    (%%esi)		/* cpy t0 into scratch sincos slot */\n\t"\
			"movaps	%%xmm6,0x10(%%esi)		/* cpy t6 into scratch sincos slot */\n\t"\
			"addpd	%%xmm1,%%xmm0			/*~A0 = A0+t1 */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* cpy t5 */\n\t"\
			"addpd	%%xmm2,%%xmm3			/*~t3 = t3+t2 */\n\t"\
			"subpd	%%xmm4,%%xmm5			/*~t5 = t5-t4 */\n\t"\
			"subpd	%%xmm2,%%xmm1			/*~t1 = t1-t2 */\n\t"\
			"subpd	%%xmm7,%%xmm6			/*~t6 = t6-t5 */\n\t"\
			"addpd	%%xmm2,%%xmm2			/* 2*t2 */\n\t"\
			"addpd	%%xmm7,%%xmm4			/*~t5 = t4+t5 */\n\t"\
			"addpd	%%xmm3,%%xmm0			/* B0 */\n\t"\
			"addpd	0x10(%%esi),%%xmm5		/* t3 = (t5-t4)+t6 */\n\t"\
			"subpd	%%xmm2,%%xmm3			/*~t2 =  (t2+t3) - 2*t2 = t3-t2 */\n\t"\
			"movaps	%%xmm4,%%xmm7			/* cpy t5 */\n\t"\
			"movaps	%%xmm0,0x190(%%edi)		/* <-B0 = s1p14i, xmm0 FREE */\n\t"\
			"subpd	%%xmm6,%%xmm4			/* t4 = ~t5-~t6 */\n\t"\
			"movaps	%%xmm1,%%xmm2			/* cpy ~t1 */\n\t"\
			"subpd	    (%%esi),%%xmm0		/* r = B0 - t0 */\n\t"\
			"mulpd	0x10(%%eax),%%xmm5		/*~t3 = t3*sx0 */\n\t"\
			"addpd	%%xmm3,%%xmm2			/* ~t1+~t2 */\n\t"\
			"mulpd	    (%%ecx),%%xmm3		/* t2 = t2*cx2 */\n\t"\
			"mulpd	0x10(%%edx),%%xmm4		/*~t4 = t4*sx3 */\n\t"\
			"mulpd	    (%%ebx),%%xmm1		/* t1 = t1*cx1 */\n\t"\
			"mulpd	0x10(%%ebx),%%xmm6		/*~t6 = t6*sx1 */\n\t"\
			"mulpd	    (%%eax),%%xmm0		/* ~r = r*(cx0-1) */\n\t"\
			"mulpd	0x10(%%ecx),%%xmm7		/*~t5 = t5*sx2 */\n\t"\
			"mulpd	    (%%edx),%%xmm2		/* t3 =(t1+t2)*cx3 */\n\t"\
			"addpd	0x190(%%edi),%%xmm0		/* t0 =~r + B0 */\n\t"\
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
			"movaps	0x280(%%edi),%%xmm1		/* t0r */\n\t"\
			"movaps	0x080(%%edi),%%xmm6		/* t3r */\n\t"\
			"subpd	%%xmm5,%%xmm1			/* B1r =t0r-t3i */\n\t"\
			"subpd	%%xmm6,%%xmm0			/* B6i =t0i-t3r */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*        2*t3i */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*        2*t3r */\n\t"\
			"addpd	%%xmm1,%%xmm5			/* B6r =t0r+t3i */\n\t"\
			"addpd	%%xmm0,%%xmm6			/* B1i =t0i+t3r */\n\t"\
			"movaps	%%xmm1,0x280(%%edi)		/* <-B1r = s1p08r */\n\t"\
			"movaps	%%xmm0,0x090(%%edi)		/* <-B6i = s1p20r */\n\t"\
			"movaps	%%xmm5,0x080(%%edi)		/* <-B6r = s1p20r */\n\t"\
			"movaps	%%xmm6,0x290(%%edi)		/* <-B1i = s1p08i */\n\t"\
			"\n\t"\
			"movaps	     (%%edi),%%xmm1		/* t1r */\n\t"\
			"movaps	0x300(%%edi),%%xmm6		/* t4r */\n\t"\
			"subpd	%%xmm7,%%xmm1			/* B2r =t1r-t4i */\n\t"\
			"subpd	%%xmm6,%%xmm2			/* B5i =t1i-t4r */\n\t"\
			"addpd	%%xmm7,%%xmm7			/*        2*t4i */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*        2*t4r */\n\t"\
			"addpd	%%xmm1,%%xmm7			/* B5r =t1r+t4i */\n\t"\
			"addpd	%%xmm2,%%xmm6			/* B2i =t1i+t4r */\n\t"\
			"movaps	%%xmm1,     (%%edi)		/* <-B2r = s1p16r */\n\t"\
			"movaps	%%xmm2,0x310(%%edi)		/* <-B5i = s1p12r */\n\t"\
			"movaps	%%xmm7,0x300(%%edi)		/* <-B5r = s1p12r */\n\t"\
			"movaps	%%xmm6,0x010(%%edi)		/* <-B2i = s1p16i */\n\t"\
			"\n\t"\
			"/* Note the order reversal on this pair of outputs: */\n\t"\
			"movaps	0x100(%%edi),%%xmm0		/* t2r */\n\t"\
			"movaps	0x200(%%edi),%%xmm5		/* t5r */\n\t"\
			"subpd	%%xmm4,%%xmm0			/* B4r =t2r-t5i */\n\t"\
			"subpd	%%xmm5,%%xmm3			/* B3i =t2i-t5r */\n\t"\
			"addpd	%%xmm4,%%xmm4			/*        2*t5i */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*        2*t5r */\n\t"\
			"addpd	%%xmm0,%%xmm4			/* B3r =t2r+t5i */\n\t"\
			"addpd	%%xmm3,%%xmm5			/* B4i =t2i+t5r */\n\t"\
			"movaps	%%xmm0,0x200(%%edi)		/* <-B4r = s1p04r */\n\t"\
			"movaps	%%xmm3,0x110(%%edi)		/* <-B3i = s1p24r */\n\t"\
			"movaps	%%xmm4,0x100(%%edi)		/* <-B3r = s1p24r */\n\t"\
			"movaps	%%xmm5,0x210(%%edi)		/* <-B4i = s1p04i */\n\t"\
			"\n\t"\
		"/* SSE2_RADIX_07_DFT(s1p03r,s1p07r,s1p11r,s1p15r,s1p19r,s1p23r,s1p27r,cc0,s1p07r,s1p15r,s1p23r,s1p03r,s1p11r,s1p19r,s1p27r): */\n\t"\
			"addl	$0x20,%%edi	/* s1p03r */\n\t"\
			"movaps	0x080(%%edi),%%xmm6			/* A1r = s1p07r */\n\t"\
			"movaps	0x300(%%edi),%%xmm1			/* A6r = s1p27r */\n\t"\
			"movaps	0x100(%%edi),%%xmm5			/* A2r = s1p11r */\n\t"\
			"movaps	0x280(%%edi),%%xmm2			/* A5r = s1p23r */\n\t"\
			"movaps	0x180(%%edi),%%xmm4			/* A3r = s1p15r */\n\t"\
			"movaps	0x200(%%edi),%%xmm3			/* A4r = s1p19r */\n\t"\
			"\n\t"\
			"subpd	%%xmm1,%%xmm6			/* t6r = A1r-A6r */\n\t"\
			"addpd	%%xmm1,%%xmm1			/*         2*A6r */\n\t"\
			"addpd	%%xmm6,%%xmm1			/* t1r = A1r+A6r */\n\t"\
			"\n\t"\
			"subpd	%%xmm2,%%xmm5			/* t5r = A2r-A5r */\n\t"\
			"addpd	%%xmm2,%%xmm2			/*         2*A5r */\n\t"\
			"addpd	%%xmm5,%%xmm2			/* t2r = A2r+A5r */\n\t"\
			"\n\t"\
			"movaps	    (%%edi),%%xmm0		/* Ar0 */\n\t"\
			"subpd	%%xmm3,%%xmm4			/* t4r = A3r-A4r */\n\t"\
			"addpd	%%xmm3,%%xmm3			/*         2*A4r */\n\t"\
			"addpd	%%xmm4,%%xmm3			/* t3r = A3r+A4r */\n\t"\
			"\n\t"\
			"movaps	%%xmm0,    (%%esi)		/* cpy t0 into scratch sincos slot */\n\t"\
			"movaps	%%xmm6,0x10(%%esi)		/* cpy t6 into scratch sincos slot */\n\t"\
			"addpd	%%xmm1,%%xmm0			/*~A0 = A0+t1 */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* cpy t5 */\n\t"\
			"addpd	%%xmm2,%%xmm3			/*~t3 = t3+t2 */\n\t"\
			"subpd	%%xmm4,%%xmm5			/*~t5 = t5-t4 */\n\t"\
			"subpd	%%xmm2,%%xmm1			/*~t1 = t1-t2 */\n\t"\
			"subpd	%%xmm7,%%xmm6			/*~t6 = t6-t5 */\n\t"\
			"addpd	%%xmm2,%%xmm2			/* 2*t2 */\n\t"\
			"addpd	%%xmm7,%%xmm4			/*~t5 = t4+t5 */\n\t"\
			"addpd	%%xmm3,%%xmm0			/* B0 */\n\t"\
			"addpd	0x10(%%esi),%%xmm5		/* t3 = (t5-t4)+t6 */\n\t"\
			"subpd	%%xmm2,%%xmm3			/*~t2 =  (t2+t3) - 2*t2 = t3-t2 */\n\t"\
			"movaps	%%xmm4,%%xmm7			/* cpy t5 */\n\t"\
			"movaps	%%xmm0,0x080(%%edi)		/* <-B0 = s1p07r, xmm0 FREE */\n\t"\
			"subpd	%%xmm6,%%xmm4			/* t4 = ~t5-~t6 */\n\t"\
			"movaps	%%xmm1,%%xmm2			/* cpy ~t1 */\n\t"\
			"subpd	    (%%esi),%%xmm0		/* r = B0 - t0 */\n\t"\
			"mulpd	0x10(%%eax),%%xmm5		/*~t3 = t3*sx0 */\n\t"\
			"addpd	%%xmm3,%%xmm2			/* ~t1+~t2 */\n\t"\
			"mulpd	    (%%ecx),%%xmm3		/* t2 = t2*cx2 */\n\t"\
			"mulpd	0x10(%%edx),%%xmm4		/*~t4 = t4*sx3 */\n\t"\
			"mulpd	    (%%ebx),%%xmm1		/* t1 = t1*cx1 */\n\t"\
			"mulpd	0x10(%%ebx),%%xmm6		/*~t6 = t6*sx1 */\n\t"\
			"mulpd	    (%%eax),%%xmm0		/* ~r = r*(cx0-1) */\n\t"\
			"mulpd	0x10(%%ecx),%%xmm7		/*~t5 = t5*sx2 */\n\t"\
			"mulpd	    (%%edx),%%xmm2		/* t3 =(t1+t2)*cx3 */\n\t"\
			"addpd	0x080(%%edi),%%xmm0		/* t0 =~r + B0 */\n\t"\
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
			"movaps	%%xmm0,0x180(%%edi)			/* B1 <- t0 = s1p15r */\n\t"\
			"movaps	%%xmm5,0x300(%%edi)			/* B6 <- t3 = s1p27r */\n\t"\
			"movaps	%%xmm2,0x280(%%edi)			/* B2 <- t1 = s1p23r */\n\t"\
			"movaps	%%xmm7,0x200(%%edi)			/* B5 <- t4 = s1p19r */\n\t"\
			"movaps	%%xmm3,     (%%edi)			/* B3 <- t2 = s1p03r */\n\t"\
			"movaps	%%xmm4,0x100(%%edi)			/* B4 <- t5 = s1p11r */\n\t"\
			"/************************** Imaginary Parts: ******************************************/\n\t"\
			"movaps	0x090(%%edi),%%xmm6		/* A1i = s1p07i */\n\t"\
			"movaps	0x310(%%edi),%%xmm1		/* A6i = s1p27i */\n\t"\
			"movaps	0x110(%%edi),%%xmm5		/* A2i = s1p11i */\n\t"\
			"movaps	0x290(%%edi),%%xmm2		/* A5i = s1p23i */\n\t"\
			"movaps	0x190(%%edi),%%xmm4		/* A3i = s1p15i */\n\t"\
			"movaps	0x210(%%edi),%%xmm3		/* A4i = s1p19i */\n\t"\
			"\n\t"\
			"subpd	%%xmm1,%%xmm6			/* t6i = A1i-A6i */\n\t"\
			"addpd	%%xmm1,%%xmm1			/*         2*A6i */\n\t"\
			"addpd	%%xmm6,%%xmm1			/* t1i = A1i+A6i */\n\t"\
			"\n\t"\
			"subpd	%%xmm2,%%xmm5			/* t5i = A2i-A5i */\n\t"\
			"addpd	%%xmm2,%%xmm2			/*         2*A5i */\n\t"\
			"addpd	%%xmm5,%%xmm2			/* t2i = A2i+A5i */\n\t"\
			"\n\t"\
			"movaps	0x10(%%edi),%%xmm0		/* Ai0 */\n\t"\
			"subpd	%%xmm3,%%xmm4			/* t4i = A3i-A4i */\n\t"\
			"addpd	%%xmm3,%%xmm3			/*         2*A4i */\n\t"\
			"addpd	%%xmm4,%%xmm3			/* t3i = A3i+A4i */\n\t"\
			"\n\t"\
			"movaps	%%xmm0,    (%%esi)		/* cpy t0 into scratch sincos slot */\n\t"\
			"movaps	%%xmm6,0x10(%%esi)		/* cpy t6 into scratch sincos slot */\n\t"\
			"addpd	%%xmm1,%%xmm0			/*~A0 = A0+t1 */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* cpy t5 */\n\t"\
			"addpd	%%xmm2,%%xmm3			/*~t3 = t3+t2 */\n\t"\
			"subpd	%%xmm4,%%xmm5			/*~t5 = t5-t4 */\n\t"\
			"subpd	%%xmm2,%%xmm1			/*~t1 = t1-t2 */\n\t"\
			"subpd	%%xmm7,%%xmm6			/*~t6 = t6-t5 */\n\t"\
			"addpd	%%xmm2,%%xmm2			/* 2*t2 */\n\t"\
			"addpd	%%xmm7,%%xmm4			/*~t5 = t4+t5 */\n\t"\
			"addpd	%%xmm3,%%xmm0			/* B0 */\n\t"\
			"addpd	0x10(%%esi),%%xmm5		/* t3 = (t5-t4)+t6 */\n\t"\
			"subpd	%%xmm2,%%xmm3			/*~t2 =  (t2+t3) - 2*t2 = t3-t2 */\n\t"\
			"movaps	%%xmm4,%%xmm7			/* cpy t5 */\n\t"\
			"movaps	%%xmm0,0x090(%%edi)		/* <-B0 = s1p07i, xmm0 FREE */\n\t"\
			"subpd	%%xmm6,%%xmm4			/* t4 = ~t5-~t6 */\n\t"\
			"movaps	%%xmm1,%%xmm2			/* cpy ~t1 */\n\t"\
			"subpd	    (%%esi),%%xmm0		/* r = B0 - t0 */\n\t"\
			"mulpd	0x10(%%eax),%%xmm5		/*~t3 = t3*sx0 */\n\t"\
			"addpd	%%xmm3,%%xmm2			/* ~t1+~t2 */\n\t"\
			"mulpd	    (%%ecx),%%xmm3		/* t2 = t2*cx2 */\n\t"\
			"mulpd	0x10(%%edx),%%xmm4		/*~t4 = t4*sx3 */\n\t"\
			"mulpd	    (%%ebx),%%xmm1		/* t1 = t1*cx1 */\n\t"\
			"mulpd	0x10(%%ebx),%%xmm6		/*~t6 = t6*sx1 */\n\t"\
			"mulpd	    (%%eax),%%xmm0		/* ~r = r*(cx0-1) */\n\t"\
			"mulpd	0x10(%%ecx),%%xmm7		/*~t5 = t5*sx2 */\n\t"\
			"mulpd	    (%%edx),%%xmm2		/* t3 =(t1+t2)*cx3 */\n\t"\
			"addpd	0x090(%%edi),%%xmm0		/* t0 =~r + B0 */\n\t"\
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
			"movaps	0x180(%%edi),%%xmm1		/* t0r */\n\t"\
			"movaps	0x300(%%edi),%%xmm6		/* t3r */\n\t"\
			"subpd	%%xmm5,%%xmm1			/* B1r =t0r-t3i */\n\t"\
			"subpd	%%xmm6,%%xmm0			/* B6i =t0i-t3r */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*        2*t3i */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*        2*t3r */\n\t"\
			"addpd	%%xmm1,%%xmm5			/* B6r =t0r+t3i */\n\t"\
			"addpd	%%xmm0,%%xmm6			/* B1i =t0i+t3r */\n\t"\
			"movaps	%%xmm1,0x180(%%edi)		/* <-B1r = s1p08r */\n\t"\
			"movaps	%%xmm0,0x310(%%edi)		/* <-B6i = s1p20r */\n\t"\
			"movaps	%%xmm5,0x300(%%edi)		/* <-B6r = s1p20r */\n\t"\
			"movaps	%%xmm6,0x190(%%edi)		/* <-B1i = s1p08i */\n\t"\
			"\n\t"\
			"movaps	0x280(%%edi),%%xmm1		/* t1r */\n\t"\
			"movaps	0x200(%%edi),%%xmm6		/* t4r */\n\t"\
			"subpd	%%xmm7,%%xmm1			/* B2r =t1r-t4i */\n\t"\
			"subpd	%%xmm6,%%xmm2			/* B5i =t1i-t4r */\n\t"\
			"addpd	%%xmm7,%%xmm7			/*        2*t4i */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*        2*t4r */\n\t"\
			"addpd	%%xmm1,%%xmm7			/* B5r =t1r+t4i */\n\t"\
			"addpd	%%xmm2,%%xmm6			/* B2i =t1i+t4r */\n\t"\
			"movaps	%%xmm1,0x280(%%edi)		/* <-B2r = s1p16r */\n\t"\
			"movaps	%%xmm2,0x210(%%edi)		/* <-B5i = s1p12r */\n\t"\
			"movaps	%%xmm7,0x200(%%edi)		/* <-B5r = s1p12r */\n\t"\
			"movaps	%%xmm6,0x290(%%edi)		/* <-B2i = s1p16i */\n\t"\
			"\n\t"\
			"/* Note the order reversal on this pair of outputs: */\n\t"\
			"movaps	     (%%edi),%%xmm0		/* t2r */\n\t"\
			"movaps	0x100(%%edi),%%xmm5		/* t5r */\n\t"\
			"subpd	%%xmm4,%%xmm0			/* B4r =t2r-t5i */\n\t"\
			"subpd	%%xmm5,%%xmm3			/* B3i =t2i-t5r */\n\t"\
			"addpd	%%xmm4,%%xmm4			/*        2*t5i */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*        2*t5r */\n\t"\
			"addpd	%%xmm0,%%xmm4			/* B3r =t2r+t5i */\n\t"\
			"addpd	%%xmm3,%%xmm5			/* B4i =t2i+t5r */\n\t"\
			"movaps	%%xmm0,0x100(%%edi)		/* <-B4r = s1p04r */\n\t"\
			"movaps	%%xmm3,0x010(%%edi)		/* <-B3i = s1p24r */\n\t"\
			"movaps	%%xmm4,     (%%edi)		/* <-B3r = s1p24r */\n\t"\
			"movaps	%%xmm5,0x110(%%edi)		/* <-B4i = s1p04i */\n\t"\
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
			: "eax","ebx","ecx","edx","edi","esi"		/* Clobbered registers */\
		);\
	}


	#define	SSE2_RADIX28_DIF_NOTWIDDLE(Xadd,Xp01,Xp02,Xp03,Xp04,Xp08,Xp12,Xp16,Xp20,Xp24,Xout,Xcc0)\
	{\
	__asm__ volatile (\
		"/* SSE2_RADIX_07_DFT(s1p00r,s1p24r,s1p20r,s1p16r,s1p12r,s1p08r,s1p04r,cc0,s1p00r,s1p04r,s1p08r,s1p12r,s1p16r,s1p20r,s1p24r): */\n\t"\
			"movl	%[__out],%%edi	/* s1p00r */\n\t"\
			"movl	%[__cc0],%%eax	/* cc0 */\n\t"\
			"movl	$0x20,%%ebx\n\t"\
			"movl	$0x40,%%ecx\n\t"\
			"movl	$0x60,%%edx\n\t"\
			"movl	$0x80,%%esi\n\t"\
			"addl	%%eax,%%ebx		/* cc1 */\n\t"\
			"addl	%%eax,%%ecx		/* cc2 */\n\t"\
			"addl	%%eax,%%edx		/* cc3 */\n\t"\
			"addl	%%eax,%%esi		/* scratch storage */\n\t"\
			"movaps	0x300(%%edi),%%xmm6			/* A1r = s1p24r */\n\t"\
			"movaps	0x080(%%edi),%%xmm1			/* A6r = s1p04r */\n\t"\
			"movaps	0x280(%%edi),%%xmm5			/* A2r = s1p20r */\n\t"\
			"movaps	0x100(%%edi),%%xmm2			/* A5r = s1p08r */\n\t"\
			"movaps	0x200(%%edi),%%xmm4			/* A3r = s1p16r */\n\t"\
			"movaps	0x180(%%edi),%%xmm3			/* A4r = s1p12r */\n\t"\
			"\n\t"\
			"subpd	%%xmm1,%%xmm6			/* t6r = A1r-A6r */\n\t"\
			"addpd	%%xmm1,%%xmm1			/*         2*A6r */\n\t"\
			"addpd	%%xmm6,%%xmm1			/* t1r = A1r+A6r */\n\t"\
			"\n\t"\
			"subpd	%%xmm2,%%xmm5			/* t5r = A2r-A5r */\n\t"\
			"addpd	%%xmm2,%%xmm2			/*         2*A5r */\n\t"\
			"addpd	%%xmm5,%%xmm2			/* t2r = A2r+A5r */\n\t"\
			"\n\t"\
			"movaps	    (%%edi),%%xmm0		/* Ar0 */\n\t"\
			"subpd	%%xmm3,%%xmm4			/* t4r = A3r-A4r */\n\t"\
			"addpd	%%xmm3,%%xmm3			/*         2*A4r */\n\t"\
			"addpd	%%xmm4,%%xmm3			/* t3r = A3r+A4r */\n\t"\
			"\n\t"\
			"movaps	%%xmm0,    (%%esi)		/* cpy t0 into scratch sincos slot */\n\t"\
			"movaps	%%xmm6,0x10(%%esi)		/* cpy t6 into scratch sincos slot */\n\t"\
			"addpd	%%xmm1,%%xmm0			/*~A0 = A0+t1 */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* cpy t5 */\n\t"\
			"addpd	%%xmm2,%%xmm3			/*~t3 = t3+t2 */\n\t"\
			"subpd	%%xmm4,%%xmm5			/*~t5 = t5-t4 */\n\t"\
			"subpd	%%xmm2,%%xmm1			/*~t1 = t1-t2 */\n\t"\
			"subpd	%%xmm7,%%xmm6			/*~t6 = t6-t5 */\n\t"\
			"addpd	%%xmm2,%%xmm2			/* 2*t2 */\n\t"\
			"addpd	%%xmm7,%%xmm4			/*~t5 = t4+t5 */\n\t"\
			"addpd	%%xmm3,%%xmm0			/* B0 */\n\t"\
			"addpd	0x10(%%esi),%%xmm5		/* t3 = (t5-t4)+t6 */\n\t"\
			"subpd	%%xmm2,%%xmm3			/*~t2 =  (t2+t3) - 2*t2 = t3-t2 */\n\t"\
			"movaps	%%xmm4,%%xmm7			/* cpy t5 */\n\t"\
			"movaps	%%xmm0,     (%%edi)		/* <-B0 = s1p00r, xmm0 FREE */\n\t"\
			"subpd	%%xmm6,%%xmm4			/* t4 = ~t5-~t6 */\n\t"\
			"movaps	%%xmm1,%%xmm2			/* cpy ~t1 */\n\t"\
			"subpd	    (%%esi),%%xmm0		/* r = B0 - t0 */\n\t"\
			"mulpd	0x10(%%eax),%%xmm5		/*~t3 = t3*sx0 */\n\t"\
			"addpd	%%xmm3,%%xmm2			/* ~t1+~t2 */\n\t"\
			"mulpd	    (%%ecx),%%xmm3		/* t2 = t2*cx2 */\n\t"\
			"mulpd	0x10(%%edx),%%xmm4		/*~t4 = t4*sx3 */\n\t"\
			"mulpd	    (%%ebx),%%xmm1		/* t1 = t1*cx1 */\n\t"\
			"mulpd	0x10(%%ebx),%%xmm6		/*~t6 = t6*sx1 */\n\t"\
			"mulpd	    (%%eax),%%xmm0		/* ~r = r*(cx0-1) */\n\t"\
			"mulpd	0x10(%%ecx),%%xmm7		/*~t5 = t5*sx2 */\n\t"\
			"mulpd	    (%%edx),%%xmm2		/* t3 =(t1+t2)*cx3 */\n\t"\
			"addpd	    (%%edi),%%xmm0		/* t0 =~r + B0 */\n\t"\
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
			"movaps	%%xmm0,0x080(%%edi)			/* B1 <- t0 = s1p04r */\n\t"\
			"movaps	%%xmm5,0x300(%%edi)			/* B6 <- t3 = s1p24r */\n\t"\
			"movaps	%%xmm2,0x100(%%edi)			/* B2 <- t1 = s1p08r */\n\t"\
			"movaps	%%xmm7,0x280(%%edi)			/* B5 <- t4 = s1p20r */\n\t"\
			"movaps	%%xmm3,0x180(%%edi)			/* B3 <- t2 = s1p12r */\n\t"\
			"movaps	%%xmm4,0x200(%%edi)			/* B4 <- t5 = s1p16r */\n\t"\
			"/************************** Imaginary Parts: ******************************************/\n\t"\
			"movaps	0x310(%%edi),%%xmm6		/* A1i = s1p24i */\n\t"\
			"movaps	0x090(%%edi),%%xmm1		/* A6i = s1p04i */\n\t"\
			"movaps	0x290(%%edi),%%xmm5		/* A2i = s1p20i */\n\t"\
			"movaps	0x110(%%edi),%%xmm2		/* A5i = s1p08i */\n\t"\
			"movaps	0x210(%%edi),%%xmm4		/* A3i = s1p16i */\n\t"\
			"movaps	0x190(%%edi),%%xmm3		/* A4i = s1p12i */\n\t"\
			"\n\t"\
			"subpd	%%xmm1,%%xmm6			/* t6i = A1i-A6i */\n\t"\
			"addpd	%%xmm1,%%xmm1			/*         2*A6i */\n\t"\
			"addpd	%%xmm6,%%xmm1			/* t1i = A1i+A6i */\n\t"\
			"\n\t"\
			"subpd	%%xmm2,%%xmm5			/* t5i = A2i-A5i */\n\t"\
			"addpd	%%xmm2,%%xmm2			/*         2*A5i */\n\t"\
			"addpd	%%xmm5,%%xmm2			/* t2i = A2i+A5i */\n\t"\
			"\n\t"\
			"movaps	0x10(%%edi),%%xmm0		/* Ai0 */\n\t"\
			"subpd	%%xmm3,%%xmm4			/* t4i = A3i-A4i */\n\t"\
			"addpd	%%xmm3,%%xmm3			/*         2*A4i */\n\t"\
			"addpd	%%xmm4,%%xmm3			/* t3i = A3i+A4i */\n\t"\
			"\n\t"\
			"movaps	%%xmm0,    (%%esi)		/* cpy t0 into scratch sincos slot */\n\t"\
			"movaps	%%xmm6,0x10(%%esi)		/* cpy t6 into scratch sincos slot */\n\t"\
			"addpd	%%xmm1,%%xmm0			/*~A0 = A0+t1 */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* cpy t5 */\n\t"\
			"addpd	%%xmm2,%%xmm3			/*~t3 = t3+t2 */\n\t"\
			"subpd	%%xmm4,%%xmm5			/*~t5 = t5-t4 */\n\t"\
			"subpd	%%xmm2,%%xmm1			/*~t1 = t1-t2 */\n\t"\
			"subpd	%%xmm7,%%xmm6			/*~t6 = t6-t5 */\n\t"\
			"addpd	%%xmm2,%%xmm2			/* 2*t2 */\n\t"\
			"addpd	%%xmm7,%%xmm4			/*~t5 = t4+t5 */\n\t"\
			"addpd	%%xmm3,%%xmm0			/* B0 */\n\t"\
			"addpd	0x10(%%esi),%%xmm5		/* t3 = (t5-t4)+t6 */\n\t"\
			"subpd	%%xmm2,%%xmm3			/*~t2 =  (t2+t3) - 2*t2 = t3-t2 */\n\t"\
			"movaps	%%xmm4,%%xmm7			/* cpy t5 */\n\t"\
			"movaps	%%xmm0,0x010(%%edi)		/* <-B0 = s1p00i, xmm0 FREE */\n\t"\
			"subpd	%%xmm6,%%xmm4			/* t4 = ~t5-~t6 */\n\t"\
			"movaps	%%xmm1,%%xmm2			/* cpy ~t1 */\n\t"\
			"subpd	    (%%esi),%%xmm0		/* r = B0 - t0 */\n\t"\
			"mulpd	0x10(%%eax),%%xmm5		/*~t3 = t3*sx0 */\n\t"\
			"addpd	%%xmm3,%%xmm2			/* ~t1+~t2 */\n\t"\
			"mulpd	    (%%ecx),%%xmm3		/* t2 = t2*cx2 */\n\t"\
			"mulpd	0x10(%%edx),%%xmm4		/*~t4 = t4*sx3 */\n\t"\
			"mulpd	    (%%ebx),%%xmm1		/* t1 = t1*cx1 */\n\t"\
			"mulpd	0x10(%%ebx),%%xmm6		/*~t6 = t6*sx1 */\n\t"\
			"mulpd	    (%%eax),%%xmm0		/* ~r = r*(cx0-1) */\n\t"\
			"mulpd	0x10(%%ecx),%%xmm7		/*~t5 = t5*sx2 */\n\t"\
			"mulpd	    (%%edx),%%xmm2		/* t3 =(t1+t2)*cx3 */\n\t"\
			"addpd	0x10(%%edi),%%xmm0		/* t0 =~r + B0 */\n\t"\
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
			"movaps	0x080(%%edi),%%xmm1		/* t0r */\n\t"\
			"movaps	0x300(%%edi),%%xmm6		/* t3r */\n\t"\
			"subpd	%%xmm5,%%xmm1			/* B1r =t0r-t3i */\n\t"\
			"subpd	%%xmm6,%%xmm0			/* B6i =t0i-t3r */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*        2*t3i */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*        2*t3r */\n\t"\
			"addpd	%%xmm1,%%xmm5			/* B6r =t0r+t3i */\n\t"\
			"addpd	%%xmm0,%%xmm6			/* B1i =t0i+t3r */\n\t"\
			"movaps	%%xmm1,0x080(%%edi)		/* <-B1r = s1p08r */\n\t"\
			"movaps	%%xmm0,0x310(%%edi)		/* <-B6i = s1p20r */\n\t"\
			"movaps	%%xmm5,0x300(%%edi)		/* <-B6r = s1p20r */\n\t"\
			"movaps	%%xmm6,0x090(%%edi)		/* <-B1i = s1p08i */\n\t"\
			"\n\t"\
			"movaps	0x100(%%edi),%%xmm1		/* t1r */\n\t"\
			"movaps	0x280(%%edi),%%xmm6		/* t4r */\n\t"\
			"subpd	%%xmm7,%%xmm1			/* B2r =t1r-t4i */\n\t"\
			"subpd	%%xmm6,%%xmm2			/* B5i =t1i-t4r */\n\t"\
			"addpd	%%xmm7,%%xmm7			/*        2*t4i */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*        2*t4r */\n\t"\
			"addpd	%%xmm1,%%xmm7			/* B5r =t1r+t4i */\n\t"\
			"addpd	%%xmm2,%%xmm6			/* B2i =t1i+t4r */\n\t"\
			"movaps	%%xmm1,0x100(%%edi)		/* <-B2r = s1p16r */\n\t"\
			"movaps	%%xmm2,0x290(%%edi)		/* <-B5i = s1p12r */\n\t"\
			"movaps	%%xmm7,0x280(%%edi)		/* <-B5r = s1p12r */\n\t"\
			"movaps	%%xmm6,0x110(%%edi)		/* <-B2i = s1p16i */\n\t"\
			"\n\t"\
			"/* Note the order reversal on this pair of outputs: */\n\t"\
			"movaps	0x180(%%edi),%%xmm0		/* t2r */\n\t"\
			"movaps	0x200(%%edi),%%xmm5		/* t5r */\n\t"\
			"subpd	%%xmm4,%%xmm0			/* B4r =t2r-t5i */\n\t"\
			"subpd	%%xmm5,%%xmm3			/* B3i =t2i-t5r */\n\t"\
			"addpd	%%xmm4,%%xmm4			/*        2*t5i */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*        2*t5r */\n\t"\
			"addpd	%%xmm0,%%xmm4			/* B3r =t2r+t5i */\n\t"\
			"addpd	%%xmm3,%%xmm5			/* B4i =t2i+t5r */\n\t"\
			"movaps	%%xmm0,0x200(%%edi)		/* <-B4r = s1p04r */\n\t"\
			"movaps	%%xmm3,0x190(%%edi)		/* <-B3i = s1p24r */\n\t"\
			"movaps	%%xmm4,0x180(%%edi)		/* <-B3r = s1p24r */\n\t"\
			"movaps	%%xmm5,0x210(%%edi)		/* <-B4i = s1p04i */\n\t"\
			"\n\t"\
		"/* SSE2_RADIX_07_DFT(s1p21r,s1p17r,s1p13r,s1p09r,s1p05r,s1p01r,s1p25r,cc0,s1p01r,s1p05r,s1p09r,s1p13r,s1p17r,s1p21r,s1p25r): */\n\t"\
			"addl	$0x20,%%edi	/* s1p01r */\n\t"\
			"movaps	0x200(%%edi),%%xmm6			/* A1r = s1p17r */\n\t"\
			"movaps	0x300(%%edi),%%xmm1			/* A6r = s1p25r */\n\t"\
			"movaps	0x180(%%edi),%%xmm5			/* A2r = s1p13r */\n\t"\
			"movaps	     (%%edi),%%xmm2			/* A5r = s1p01r */\n\t"\
			"movaps	0x100(%%edi),%%xmm4			/* A3r = s1p09r */\n\t"\
			"movaps	0x080(%%edi),%%xmm3			/* A4r = s1p05r */\n\t"\
			"\n\t"\
			"subpd	%%xmm1,%%xmm6			/* t6r = A1r-A6r */\n\t"\
			"addpd	%%xmm1,%%xmm1			/*         2*A6r */\n\t"\
			"addpd	%%xmm6,%%xmm1			/* t1r = A1r+A6r */\n\t"\
			"\n\t"\
			"subpd	%%xmm2,%%xmm5			/* t5r = A2r-A5r */\n\t"\
			"addpd	%%xmm2,%%xmm2			/*         2*A5r */\n\t"\
			"addpd	%%xmm5,%%xmm2			/* t2r = A2r+A5r */\n\t"\
			"\n\t"\
			"movaps	0x280(%%edi),%%xmm0		/* Ar0 = s1p21r */\n\t"\
			"subpd	%%xmm3,%%xmm4			/* t4r = A3r-A4r */\n\t"\
			"addpd	%%xmm3,%%xmm3			/*         2*A4r */\n\t"\
			"addpd	%%xmm4,%%xmm3			/* t3r = A3r+A4r */\n\t"\
			"\n\t"\
			"movaps	%%xmm0,    (%%esi)		/* cpy t0 into scratch sincos slot */\n\t"\
			"movaps	%%xmm6,0x10(%%esi)		/* cpy t6 into scratch sincos slot */\n\t"\
			"addpd	%%xmm1,%%xmm0			/*~A0 = A0+t1 */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* cpy t5 */\n\t"\
			"addpd	%%xmm2,%%xmm3			/*~t3 = t3+t2 */\n\t"\
			"subpd	%%xmm4,%%xmm5			/*~t5 = t5-t4 */\n\t"\
			"subpd	%%xmm2,%%xmm1			/*~t1 = t1-t2 */\n\t"\
			"subpd	%%xmm7,%%xmm6			/*~t6 = t6-t5 */\n\t"\
			"addpd	%%xmm2,%%xmm2			/* 2*t2 */\n\t"\
			"addpd	%%xmm7,%%xmm4			/*~t5 = t4+t5 */\n\t"\
			"addpd	%%xmm3,%%xmm0			/* B0 */\n\t"\
			"addpd	0x10(%%esi),%%xmm5		/* t3 = (t5-t4)+t6 */\n\t"\
			"subpd	%%xmm2,%%xmm3			/*~t2 =  (t2+t3) - 2*t2 = t3-t2 */\n\t"\
			"movaps	%%xmm4,%%xmm7			/* cpy t5 */\n\t"\
			"movaps	%%xmm0,     (%%edi)		/* <-B0 = s1p01r, xmm0 FREE */\n\t"\
			"subpd	%%xmm6,%%xmm4			/* t4 = ~t5-~t6 */\n\t"\
			"movaps	%%xmm1,%%xmm2			/* cpy ~t1 */\n\t"\
			"subpd	    (%%esi),%%xmm0		/* r = B0 - t0 */\n\t"\
			"mulpd	0x10(%%eax),%%xmm5		/*~t3 = t3*sx0 */\n\t"\
			"addpd	%%xmm3,%%xmm2			/* ~t1+~t2 */\n\t"\
			"mulpd	    (%%ecx),%%xmm3		/* t2 = t2*cx2 */\n\t"\
			"mulpd	0x10(%%edx),%%xmm4		/*~t4 = t4*sx3 */\n\t"\
			"mulpd	    (%%ebx),%%xmm1		/* t1 = t1*cx1 */\n\t"\
			"mulpd	0x10(%%ebx),%%xmm6		/*~t6 = t6*sx1 */\n\t"\
			"mulpd	    (%%eax),%%xmm0		/* ~r = r*(cx0-1) */\n\t"\
			"mulpd	0x10(%%ecx),%%xmm7		/*~t5 = t5*sx2 */\n\t"\
			"mulpd	    (%%edx),%%xmm2		/* t3 =(t1+t2)*cx3 */\n\t"\
			"addpd	    (%%edi),%%xmm0		/* t0 =~r + B0 */\n\t"\
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
			"movaps	%%xmm0,0x080(%%edi)			/* B1 <- t0 = s1p05r */\n\t"\
			"movaps	%%xmm5,0x300(%%edi)			/* B6 <- t3 = s1p25r */\n\t"\
			"movaps	%%xmm2,0x100(%%edi)			/* B2 <- t1 = s1p09r */\n\t"\
			"movaps	%%xmm7,0x280(%%edi)			/* B5 <- t4 = s1p21r */\n\t"\
			"movaps	%%xmm3,0x180(%%edi)			/* B3 <- t2 = s1p13r */\n\t"\
			"movaps	%%xmm4,0x200(%%edi)			/* B4 <- t5 = s1p17r */\n\t"\
			"/************************** Imaginary Parts: ******************************************/\n\t"\
			"movaps	0x210(%%edi),%%xmm6		/* A1i = s1p17i */\n\t"\
			"movaps	0x310(%%edi),%%xmm1		/* A6i = s1p25i */\n\t"\
			"movaps	0x190(%%edi),%%xmm5		/* A2i = s1p13i */\n\t"\
			"movaps	0x010(%%edi),%%xmm2		/* A5i = s1p01i */\n\t"\
			"movaps	0x110(%%edi),%%xmm4		/* A3i = s1p09i */\n\t"\
			"movaps	0x090(%%edi),%%xmm3		/* A4i = s1p05i */\n\t"\
			"\n\t"\
			"subpd	%%xmm1,%%xmm6			/* t6i = A1i-A6i */\n\t"\
			"addpd	%%xmm1,%%xmm1			/*         2*A6i */\n\t"\
			"addpd	%%xmm6,%%xmm1			/* t1i = A1i+A6i */\n\t"\
			"\n\t"\
			"subpd	%%xmm2,%%xmm5			/* t5i = A2i-A5i */\n\t"\
			"addpd	%%xmm2,%%xmm2			/*         2*A5i */\n\t"\
			"addpd	%%xmm5,%%xmm2			/* t2i = A2i+A5i */\n\t"\
			"\n\t"\
			"movaps	0x290(%%edi),%%xmm0		/* Ai0 = s1p21i */\n\t"\
			"subpd	%%xmm3,%%xmm4			/* t4i = A3i-A4i */\n\t"\
			"addpd	%%xmm3,%%xmm3			/*         2*A4i */\n\t"\
			"addpd	%%xmm4,%%xmm3			/* t3i = A3i+A4i */\n\t"\
			"\n\t"\
			"movaps	%%xmm0,    (%%esi)		/* cpy t0 into scratch sincos slot */\n\t"\
			"movaps	%%xmm6,0x10(%%esi)		/* cpy t6 into scratch sincos slot */\n\t"\
			"addpd	%%xmm1,%%xmm0			/*~A0 = A0+t1 */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* cpy t5 */\n\t"\
			"addpd	%%xmm2,%%xmm3			/*~t3 = t3+t2 */\n\t"\
			"subpd	%%xmm4,%%xmm5			/*~t5 = t5-t4 */\n\t"\
			"subpd	%%xmm2,%%xmm1			/*~t1 = t1-t2 */\n\t"\
			"subpd	%%xmm7,%%xmm6			/*~t6 = t6-t5 */\n\t"\
			"addpd	%%xmm2,%%xmm2			/* 2*t2 */\n\t"\
			"addpd	%%xmm7,%%xmm4			/*~t5 = t4+t5 */\n\t"\
			"addpd	%%xmm3,%%xmm0			/* B0 */\n\t"\
			"addpd	0x10(%%esi),%%xmm5		/* t3 = (t5-t4)+t6 */\n\t"\
			"subpd	%%xmm2,%%xmm3			/*~t2 =  (t2+t3) - 2*t2 = t3-t2 */\n\t"\
			"movaps	%%xmm4,%%xmm7			/* cpy t5 */\n\t"\
			"movaps	%%xmm0,0x010(%%edi)		/* <-B0 = s1p00i, xmm0 FREE */\n\t"\
			"subpd	%%xmm6,%%xmm4			/* t4 = ~t5-~t6 */\n\t"\
			"movaps	%%xmm1,%%xmm2			/* cpy ~t1 */\n\t"\
			"subpd	    (%%esi),%%xmm0		/* r = B0 - t0 */\n\t"\
			"mulpd	0x10(%%eax),%%xmm5		/*~t3 = t3*sx0 */\n\t"\
			"addpd	%%xmm3,%%xmm2			/* ~t1+~t2 */\n\t"\
			"mulpd	    (%%ecx),%%xmm3		/* t2 = t2*cx2 */\n\t"\
			"mulpd	0x10(%%edx),%%xmm4		/*~t4 = t4*sx3 */\n\t"\
			"mulpd	    (%%ebx),%%xmm1		/* t1 = t1*cx1 */\n\t"\
			"mulpd	0x10(%%ebx),%%xmm6		/*~t6 = t6*sx1 */\n\t"\
			"mulpd	    (%%eax),%%xmm0		/* ~r = r*(cx0-1) */\n\t"\
			"mulpd	0x10(%%ecx),%%xmm7		/*~t5 = t5*sx2 */\n\t"\
			"mulpd	    (%%edx),%%xmm2		/* t3 =(t1+t2)*cx3 */\n\t"\
			"addpd	0x10(%%edi),%%xmm0		/* t0 =~r + B0 */\n\t"\
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
			"movaps	0x080(%%edi),%%xmm1		/* t0r */\n\t"\
			"movaps	0x300(%%edi),%%xmm6		/* t3r */\n\t"\
			"subpd	%%xmm5,%%xmm1			/* B1r =t0r-t3i */\n\t"\
			"subpd	%%xmm6,%%xmm0			/* B6i =t0i-t3r */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*        2*t3i */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*        2*t3r */\n\t"\
			"addpd	%%xmm1,%%xmm5			/* B6r =t0r+t3i */\n\t"\
			"addpd	%%xmm0,%%xmm6			/* B1i =t0i+t3r */\n\t"\
			"movaps	%%xmm1,0x080(%%edi)		/* <-B1r = s1p08r */\n\t"\
			"movaps	%%xmm0,0x310(%%edi)		/* <-B6i = s1p20r */\n\t"\
			"movaps	%%xmm5,0x300(%%edi)		/* <-B6r = s1p20r */\n\t"\
			"movaps	%%xmm6,0x090(%%edi)		/* <-B1i = s1p08i */\n\t"\
			"\n\t"\
			"movaps	0x100(%%edi),%%xmm1		/* t1r */\n\t"\
			"movaps	0x280(%%edi),%%xmm6		/* t4r */\n\t"\
			"subpd	%%xmm7,%%xmm1			/* B2r =t1r-t4i */\n\t"\
			"subpd	%%xmm6,%%xmm2			/* B5i =t1i-t4r */\n\t"\
			"addpd	%%xmm7,%%xmm7			/*        2*t4i */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*        2*t4r */\n\t"\
			"addpd	%%xmm1,%%xmm7			/* B5r =t1r+t4i */\n\t"\
			"addpd	%%xmm2,%%xmm6			/* B2i =t1i+t4r */\n\t"\
			"movaps	%%xmm1,0x100(%%edi)		/* <-B2r = s1p16r */\n\t"\
			"movaps	%%xmm2,0x290(%%edi)		/* <-B5i = s1p12r */\n\t"\
			"movaps	%%xmm7,0x280(%%edi)		/* <-B5r = s1p12r */\n\t"\
			"movaps	%%xmm6,0x110(%%edi)		/* <-B2i = s1p16i */\n\t"\
			"\n\t"\
			"/* Note the order reversal on this pair of outputs: */\n\t"\
			"movaps	0x180(%%edi),%%xmm0		/* t2r */\n\t"\
			"movaps	0x200(%%edi),%%xmm5		/* t5r */\n\t"\
			"subpd	%%xmm4,%%xmm0			/* B4r =t2r-t5i */\n\t"\
			"subpd	%%xmm5,%%xmm3			/* B3i =t2i-t5r */\n\t"\
			"addpd	%%xmm4,%%xmm4			/*        2*t5i */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*        2*t5r */\n\t"\
			"addpd	%%xmm0,%%xmm4			/* B3r =t2r+t5i */\n\t"\
			"addpd	%%xmm3,%%xmm5			/* B4i =t2i+t5r */\n\t"\
			"movaps	%%xmm0,0x200(%%edi)		/* <-B4r = s1p04r */\n\t"\
			"movaps	%%xmm3,0x190(%%edi)		/* <-B3i = s1p24r */\n\t"\
			"movaps	%%xmm4,0x180(%%edi)		/* <-B3r = s1p24r */\n\t"\
			"movaps	%%xmm5,0x210(%%edi)		/* <-B4i = s1p04i */\n\t"\
			"\n\t"\
		"/* SSE2_RADIX_07_DFT(s1p14r,s1p10r,s1p06r,s1p02r,s1p26r,s1p22r,s1p18r,cc0,s1p02r,s1p06r,s1p10r,s1p14r,s1p18r,s1p22r,s1p26r): */\n\t"\
			"addl	$0x20,%%edi	/* s1p02r */\n\t"\
			"movaps	0x100(%%edi),%%xmm6			/* A1r = s1p10r */\n\t"\
			"movaps	0x200(%%edi),%%xmm1			/* A6r = s1p18r */\n\t"\
			"movaps	0x080(%%edi),%%xmm5			/* A2r = s1p06r */\n\t"\
			"movaps	0x280(%%edi),%%xmm2			/* A5r = s1p22r */\n\t"\
			"movaps	     (%%edi),%%xmm4			/* A3r = s1p02r */\n\t"\
			"movaps	0x300(%%edi),%%xmm3			/* A4r = s1p26r */\n\t"\
			"\n\t"\
			"subpd	%%xmm1,%%xmm6			/* t6r = A1r-A6r */\n\t"\
			"addpd	%%xmm1,%%xmm1			/*         2*A6r */\n\t"\
			"addpd	%%xmm6,%%xmm1			/* t1r = A1r+A6r */\n\t"\
			"\n\t"\
			"subpd	%%xmm2,%%xmm5			/* t5r = A2r-A5r */\n\t"\
			"addpd	%%xmm2,%%xmm2			/*         2*A5r */\n\t"\
			"addpd	%%xmm5,%%xmm2			/* t2r = A2r+A5r */\n\t"\
			"\n\t"\
			"movaps	0x180(%%edi),%%xmm0		/* Ar0 = s1p14r */\n\t"\
			"subpd	%%xmm3,%%xmm4			/* t4r = A3r-A4r */\n\t"\
			"addpd	%%xmm3,%%xmm3			/*         2*A4r */\n\t"\
			"addpd	%%xmm4,%%xmm3			/* t3r = A3r+A4r */\n\t"\
			"\n\t"\
			"movaps	%%xmm0,    (%%esi)		/* cpy t0 into scratch sincos slot */\n\t"\
			"movaps	%%xmm6,0x10(%%esi)		/* cpy t6 into scratch sincos slot */\n\t"\
			"addpd	%%xmm1,%%xmm0			/*~A0 = A0+t1 */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* cpy t5 */\n\t"\
			"addpd	%%xmm2,%%xmm3			/*~t3 = t3+t2 */\n\t"\
			"subpd	%%xmm4,%%xmm5			/*~t5 = t5-t4 */\n\t"\
			"subpd	%%xmm2,%%xmm1			/*~t1 = t1-t2 */\n\t"\
			"subpd	%%xmm7,%%xmm6			/*~t6 = t6-t5 */\n\t"\
			"addpd	%%xmm2,%%xmm2			/* 2*t2 */\n\t"\
			"addpd	%%xmm7,%%xmm4			/*~t5 = t4+t5 */\n\t"\
			"addpd	%%xmm3,%%xmm0			/* B0 */\n\t"\
			"addpd	0x10(%%esi),%%xmm5		/* t3 = (t5-t4)+t6 */\n\t"\
			"subpd	%%xmm2,%%xmm3			/*~t2 =  (t2+t3) - 2*t2 = t3-t2 */\n\t"\
			"movaps	%%xmm4,%%xmm7			/* cpy t5 */\n\t"\
			"movaps	%%xmm0,     (%%edi)		/* <-B0 = s1p01r, xmm0 FREE */\n\t"\
			"subpd	%%xmm6,%%xmm4			/* t4 = ~t5-~t6 */\n\t"\
			"movaps	%%xmm1,%%xmm2			/* cpy ~t1 */\n\t"\
			"subpd	    (%%esi),%%xmm0		/* r = B0 - t0 */\n\t"\
			"mulpd	0x10(%%eax),%%xmm5		/*~t3 = t3*sx0 */\n\t"\
			"addpd	%%xmm3,%%xmm2			/* ~t1+~t2 */\n\t"\
			"mulpd	    (%%ecx),%%xmm3		/* t2 = t2*cx2 */\n\t"\
			"mulpd	0x10(%%edx),%%xmm4		/*~t4 = t4*sx3 */\n\t"\
			"mulpd	    (%%ebx),%%xmm1		/* t1 = t1*cx1 */\n\t"\
			"mulpd	0x10(%%ebx),%%xmm6		/*~t6 = t6*sx1 */\n\t"\
			"mulpd	    (%%eax),%%xmm0		/* ~r = r*(cx0-1) */\n\t"\
			"mulpd	0x10(%%ecx),%%xmm7		/*~t5 = t5*sx2 */\n\t"\
			"mulpd	    (%%edx),%%xmm2		/* t3 =(t1+t2)*cx3 */\n\t"\
			"addpd	    (%%edi),%%xmm0		/* t0 =~r + B0 */\n\t"\
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
			"movaps	%%xmm0,0x080(%%edi)			/* B1 <- t0 = s1p05r */\n\t"\
			"movaps	%%xmm5,0x300(%%edi)			/* B6 <- t3 = s1p25r */\n\t"\
			"movaps	%%xmm2,0x100(%%edi)			/* B2 <- t1 = s1p09r */\n\t"\
			"movaps	%%xmm7,0x280(%%edi)			/* B5 <- t4 = s1p21r */\n\t"\
			"movaps	%%xmm3,0x180(%%edi)			/* B3 <- t2 = s1p13r */\n\t"\
			"movaps	%%xmm4,0x200(%%edi)			/* B4 <- t5 = s1p17r */\n\t"\
			"/************************** Imaginary Parts: ******************************************/\n\t"\
			"movaps	0x110(%%edi),%%xmm6		/* A1i = s1p10i */\n\t"\
			"movaps	0x210(%%edi),%%xmm1		/* A6i = s1p18i */\n\t"\
			"movaps	0x090(%%edi),%%xmm5		/* A2i = s1p06i */\n\t"\
			"movaps	0x290(%%edi),%%xmm2		/* A5i = s1p22i */\n\t"\
			"movaps	0x010(%%edi),%%xmm4		/* A3i = s1p02i */\n\t"\
			"movaps	0x310(%%edi),%%xmm3		/* A4i = s1p26i */\n\t"\
			"\n\t"\
			"subpd	%%xmm1,%%xmm6			/* t6i = A1i-A6i */\n\t"\
			"addpd	%%xmm1,%%xmm1			/*         2*A6i */\n\t"\
			"addpd	%%xmm6,%%xmm1			/* t1i = A1i+A6i */\n\t"\
			"\n\t"\
			"subpd	%%xmm2,%%xmm5			/* t5i = A2i-A5i */\n\t"\
			"addpd	%%xmm2,%%xmm2			/*         2*A5i */\n\t"\
			"addpd	%%xmm5,%%xmm2			/* t2i = A2i+A5i */\n\t"\
			"\n\t"\
			"movaps	0x190(%%edi),%%xmm0		/* Ai0 = s1p14i */\n\t"\
			"subpd	%%xmm3,%%xmm4			/* t4i = A3i-A4i */\n\t"\
			"addpd	%%xmm3,%%xmm3			/*         2*A4i */\n\t"\
			"addpd	%%xmm4,%%xmm3			/* t3i = A3i+A4i */\n\t"\
			"\n\t"\
			"movaps	%%xmm0,    (%%esi)		/* cpy t0 into scratch sincos slot */\n\t"\
			"movaps	%%xmm6,0x10(%%esi)		/* cpy t6 into scratch sincos slot */\n\t"\
			"addpd	%%xmm1,%%xmm0			/*~A0 = A0+t1 */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* cpy t5 */\n\t"\
			"addpd	%%xmm2,%%xmm3			/*~t3 = t3+t2 */\n\t"\
			"subpd	%%xmm4,%%xmm5			/*~t5 = t5-t4 */\n\t"\
			"subpd	%%xmm2,%%xmm1			/*~t1 = t1-t2 */\n\t"\
			"subpd	%%xmm7,%%xmm6			/*~t6 = t6-t5 */\n\t"\
			"addpd	%%xmm2,%%xmm2			/* 2*t2 */\n\t"\
			"addpd	%%xmm7,%%xmm4			/*~t5 = t4+t5 */\n\t"\
			"addpd	%%xmm3,%%xmm0			/* B0 */\n\t"\
			"addpd	0x10(%%esi),%%xmm5		/* t3 = (t5-t4)+t6 */\n\t"\
			"subpd	%%xmm2,%%xmm3			/*~t2 =  (t2+t3) - 2*t2 = t3-t2 */\n\t"\
			"movaps	%%xmm4,%%xmm7			/* cpy t5 */\n\t"\
			"movaps	%%xmm0,0x010(%%edi)		/* <-B0 = s1p00i, xmm0 FREE */\n\t"\
			"subpd	%%xmm6,%%xmm4			/* t4 = ~t5-~t6 */\n\t"\
			"movaps	%%xmm1,%%xmm2			/* cpy ~t1 */\n\t"\
			"subpd	    (%%esi),%%xmm0		/* r = B0 - t0 */\n\t"\
			"mulpd	0x10(%%eax),%%xmm5		/*~t3 = t3*sx0 */\n\t"\
			"addpd	%%xmm3,%%xmm2			/* ~t1+~t2 */\n\t"\
			"mulpd	    (%%ecx),%%xmm3		/* t2 = t2*cx2 */\n\t"\
			"mulpd	0x10(%%edx),%%xmm4		/*~t4 = t4*sx3 */\n\t"\
			"mulpd	    (%%ebx),%%xmm1		/* t1 = t1*cx1 */\n\t"\
			"mulpd	0x10(%%ebx),%%xmm6		/*~t6 = t6*sx1 */\n\t"\
			"mulpd	    (%%eax),%%xmm0		/* ~r = r*(cx0-1) */\n\t"\
			"mulpd	0x10(%%ecx),%%xmm7		/*~t5 = t5*sx2 */\n\t"\
			"mulpd	    (%%edx),%%xmm2		/* t3 =(t1+t2)*cx3 */\n\t"\
			"addpd	0x10(%%edi),%%xmm0		/* t0 =~r + B0 */\n\t"\
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
			"movaps	0x080(%%edi),%%xmm1		/* t0r */\n\t"\
			"movaps	0x300(%%edi),%%xmm6		/* t3r */\n\t"\
			"subpd	%%xmm5,%%xmm1			/* B1r =t0r-t3i */\n\t"\
			"subpd	%%xmm6,%%xmm0			/* B6i =t0i-t3r */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*        2*t3i */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*        2*t3r */\n\t"\
			"addpd	%%xmm1,%%xmm5			/* B6r =t0r+t3i */\n\t"\
			"addpd	%%xmm0,%%xmm6			/* B1i =t0i+t3r */\n\t"\
			"movaps	%%xmm1,0x080(%%edi)		/* <-B1r = s1p08r */\n\t"\
			"movaps	%%xmm0,0x310(%%edi)		/* <-B6i = s1p20r */\n\t"\
			"movaps	%%xmm5,0x300(%%edi)		/* <-B6r = s1p20r */\n\t"\
			"movaps	%%xmm6,0x090(%%edi)		/* <-B1i = s1p08i */\n\t"\
			"\n\t"\
			"movaps	0x100(%%edi),%%xmm1		/* t1r */\n\t"\
			"movaps	0x280(%%edi),%%xmm6		/* t4r */\n\t"\
			"subpd	%%xmm7,%%xmm1			/* B2r =t1r-t4i */\n\t"\
			"subpd	%%xmm6,%%xmm2			/* B5i =t1i-t4r */\n\t"\
			"addpd	%%xmm7,%%xmm7			/*        2*t4i */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*        2*t4r */\n\t"\
			"addpd	%%xmm1,%%xmm7			/* B5r =t1r+t4i */\n\t"\
			"addpd	%%xmm2,%%xmm6			/* B2i =t1i+t4r */\n\t"\
			"movaps	%%xmm1,0x100(%%edi)		/* <-B2r = s1p16r */\n\t"\
			"movaps	%%xmm2,0x290(%%edi)		/* <-B5i = s1p12r */\n\t"\
			"movaps	%%xmm7,0x280(%%edi)		/* <-B5r = s1p12r */\n\t"\
			"movaps	%%xmm6,0x110(%%edi)		/* <-B2i = s1p16i */\n\t"\
			"\n\t"\
			"/* Note the order reversal on this pair of outputs: */\n\t"\
			"movaps	0x180(%%edi),%%xmm0		/* t2r */\n\t"\
			"movaps	0x200(%%edi),%%xmm5		/* t5r */\n\t"\
			"subpd	%%xmm4,%%xmm0			/* B4r =t2r-t5i */\n\t"\
			"subpd	%%xmm5,%%xmm3			/* B3i =t2i-t5r */\n\t"\
			"addpd	%%xmm4,%%xmm4			/*        2*t5i */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*        2*t5r */\n\t"\
			"addpd	%%xmm0,%%xmm4			/* B3r =t2r+t5i */\n\t"\
			"addpd	%%xmm3,%%xmm5			/* B4i =t2i+t5r */\n\t"\
			"movaps	%%xmm0,0x200(%%edi)		/* <-B4r = s1p04r */\n\t"\
			"movaps	%%xmm3,0x190(%%edi)		/* <-B3i = s1p24r */\n\t"\
			"movaps	%%xmm4,0x180(%%edi)		/* <-B3r = s1p24r */\n\t"\
			"movaps	%%xmm5,0x210(%%edi)		/* <-B4i = s1p04i */\n\t"\
			"\n\t"\
		"/* SSE2_RADIX_07_DFT(s1p07r,s1p03r,s1p27r,s1p23r,s1p19r,s1p15r,s1p11r,cc0,s1p03r,s1p07r,s1p11r,s1p15r,s1p19r,s1p23r,s1p27r): */\n\t"\
			"addl	$0x20,%%edi	/* s1p03r */\n\t"\
			"movaps	     (%%edi),%%xmm6			/* A1r = s1p03r */\n\t"\
			"movaps	0x100(%%edi),%%xmm1			/* A6r = s1p11r */\n\t"\
			"movaps	0x300(%%edi),%%xmm5			/* A2r = s1p27r */\n\t"\
			"movaps	0x180(%%edi),%%xmm2			/* A5r = s1p15r */\n\t"\
			"movaps	0x280(%%edi),%%xmm4			/* A3r = s1p23r */\n\t"\
			"movaps	0x200(%%edi),%%xmm3			/* A4r = s1p19r */\n\t"\
			"\n\t"\
			"subpd	%%xmm1,%%xmm6			/* t6r = A1r-A6r */\n\t"\
			"addpd	%%xmm1,%%xmm1			/*         2*A6r */\n\t"\
			"addpd	%%xmm6,%%xmm1			/* t1r = A1r+A6r */\n\t"\
			"\n\t"\
			"subpd	%%xmm2,%%xmm5			/* t5r = A2r-A5r */\n\t"\
			"addpd	%%xmm2,%%xmm2			/*         2*A5r */\n\t"\
			"addpd	%%xmm5,%%xmm2			/* t2r = A2r+A5r */\n\t"\
			"\n\t"\
			"movaps	0x080(%%edi),%%xmm0		/* Ar0 = s1p07r */\n\t"\
			"subpd	%%xmm3,%%xmm4			/* t4r = A3r-A4r */\n\t"\
			"addpd	%%xmm3,%%xmm3			/*         2*A4r */\n\t"\
			"addpd	%%xmm4,%%xmm3			/* t3r = A3r+A4r */\n\t"\
			"\n\t"\
			"movaps	%%xmm0,    (%%esi)		/* cpy t0 into scratch sincos slot */\n\t"\
			"movaps	%%xmm6,0x10(%%esi)		/* cpy t6 into scratch sincos slot */\n\t"\
			"addpd	%%xmm1,%%xmm0			/*~A0 = A0+t1 */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* cpy t5 */\n\t"\
			"addpd	%%xmm2,%%xmm3			/*~t3 = t3+t2 */\n\t"\
			"subpd	%%xmm4,%%xmm5			/*~t5 = t5-t4 */\n\t"\
			"subpd	%%xmm2,%%xmm1			/*~t1 = t1-t2 */\n\t"\
			"subpd	%%xmm7,%%xmm6			/*~t6 = t6-t5 */\n\t"\
			"addpd	%%xmm2,%%xmm2			/* 2*t2 */\n\t"\
			"addpd	%%xmm7,%%xmm4			/*~t5 = t4+t5 */\n\t"\
			"addpd	%%xmm3,%%xmm0			/* B0 */\n\t"\
			"addpd	0x10(%%esi),%%xmm5		/* t3 = (t5-t4)+t6 */\n\t"\
			"subpd	%%xmm2,%%xmm3			/*~t2 =  (t2+t3) - 2*t2 = t3-t2 */\n\t"\
			"movaps	%%xmm4,%%xmm7			/* cpy t5 */\n\t"\
			"movaps	%%xmm0,     (%%edi)		/* <-B0 = s1p01r, xmm0 FREE */\n\t"\
			"subpd	%%xmm6,%%xmm4			/* t4 = ~t5-~t6 */\n\t"\
			"movaps	%%xmm1,%%xmm2			/* cpy ~t1 */\n\t"\
			"subpd	    (%%esi),%%xmm0		/* r = B0 - t0 */\n\t"\
			"mulpd	0x10(%%eax),%%xmm5		/*~t3 = t3*sx0 */\n\t"\
			"addpd	%%xmm3,%%xmm2			/* ~t1+~t2 */\n\t"\
			"mulpd	    (%%ecx),%%xmm3		/* t2 = t2*cx2 */\n\t"\
			"mulpd	0x10(%%edx),%%xmm4		/*~t4 = t4*sx3 */\n\t"\
			"mulpd	    (%%ebx),%%xmm1		/* t1 = t1*cx1 */\n\t"\
			"mulpd	0x10(%%ebx),%%xmm6		/*~t6 = t6*sx1 */\n\t"\
			"mulpd	    (%%eax),%%xmm0		/* ~r = r*(cx0-1) */\n\t"\
			"mulpd	0x10(%%ecx),%%xmm7		/*~t5 = t5*sx2 */\n\t"\
			"mulpd	    (%%edx),%%xmm2		/* t3 =(t1+t2)*cx3 */\n\t"\
			"addpd	    (%%edi),%%xmm0		/* t0 =~r + B0 */\n\t"\
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
			"movaps	%%xmm0,0x080(%%edi)			/* B1 <- t0 = s1p05r */\n\t"\
			"movaps	%%xmm5,0x300(%%edi)			/* B6 <- t3 = s1p25r */\n\t"\
			"movaps	%%xmm2,0x100(%%edi)			/* B2 <- t1 = s1p09r */\n\t"\
			"movaps	%%xmm7,0x280(%%edi)			/* B5 <- t4 = s1p21r */\n\t"\
			"movaps	%%xmm3,0x180(%%edi)			/* B3 <- t2 = s1p13r */\n\t"\
			"movaps	%%xmm4,0x200(%%edi)			/* B4 <- t5 = s1p17r */\n\t"\
			"/************************** Imaginary Parts: ******************************************/\n\t"\
			"movaps	0x010(%%edi),%%xmm6			/* A1i = s1p03i */\n\t"\
			"movaps	0x110(%%edi),%%xmm1			/* A6i = s1p11i */\n\t"\
			"movaps	0x310(%%edi),%%xmm5			/* A2i = s1p27i */\n\t"\
			"movaps	0x190(%%edi),%%xmm2			/* A5i = s1p15i */\n\t"\
			"movaps	0x290(%%edi),%%xmm4			/* A3i = s1p23i */\n\t"\
			"movaps	0x210(%%edi),%%xmm3			/* A4i = s1p19i */\n\t"\
			"\n\t"\
			"subpd	%%xmm1,%%xmm6			/* t6i = A1i-A6i */\n\t"\
			"addpd	%%xmm1,%%xmm1			/*         2*A6i */\n\t"\
			"addpd	%%xmm6,%%xmm1			/* t1i = A1i+A6i */\n\t"\
			"\n\t"\
			"subpd	%%xmm2,%%xmm5			/* t5i = A2i-A5i */\n\t"\
			"addpd	%%xmm2,%%xmm2			/*         2*A5i */\n\t"\
			"addpd	%%xmm5,%%xmm2			/* t2i = A2i+A5i */\n\t"\
			"\n\t"\
			"movaps	0x090(%%edi),%%xmm0		/* Ai0 = s1p07i */\n\t"\
			"subpd	%%xmm3,%%xmm4			/* t4i = A3i-A4i */\n\t"\
			"addpd	%%xmm3,%%xmm3			/*         2*A4i */\n\t"\
			"addpd	%%xmm4,%%xmm3			/* t3i = A3i+A4i */\n\t"\
			"\n\t"\
			"movaps	%%xmm0,    (%%esi)		/* cpy t0 into scratch sincos slot */\n\t"\
			"movaps	%%xmm6,0x10(%%esi)		/* cpy t6 into scratch sincos slot */\n\t"\
			"addpd	%%xmm1,%%xmm0			/*~A0 = A0+t1 */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* cpy t5 */\n\t"\
			"addpd	%%xmm2,%%xmm3			/*~t3 = t3+t2 */\n\t"\
			"subpd	%%xmm4,%%xmm5			/*~t5 = t5-t4 */\n\t"\
			"subpd	%%xmm2,%%xmm1			/*~t1 = t1-t2 */\n\t"\
			"subpd	%%xmm7,%%xmm6			/*~t6 = t6-t5 */\n\t"\
			"addpd	%%xmm2,%%xmm2			/* 2*t2 */\n\t"\
			"addpd	%%xmm7,%%xmm4			/*~t5 = t4+t5 */\n\t"\
			"addpd	%%xmm3,%%xmm0			/* B0 */\n\t"\
			"addpd	0x10(%%esi),%%xmm5		/* t3 = (t5-t4)+t6 */\n\t"\
			"subpd	%%xmm2,%%xmm3			/*~t2 =  (t2+t3) - 2*t2 = t3-t2 */\n\t"\
			"movaps	%%xmm4,%%xmm7			/* cpy t5 */\n\t"\
			"movaps	%%xmm0,0x010(%%edi)		/* <-B0 = s1p00i, xmm0 FREE */\n\t"\
			"subpd	%%xmm6,%%xmm4			/* t4 = ~t5-~t6 */\n\t"\
			"movaps	%%xmm1,%%xmm2			/* cpy ~t1 */\n\t"\
			"subpd	    (%%esi),%%xmm0		/* r = B0 - t0 */\n\t"\
			"mulpd	0x10(%%eax),%%xmm5		/*~t3 = t3*sx0 */\n\t"\
			"addpd	%%xmm3,%%xmm2			/* ~t1+~t2 */\n\t"\
			"mulpd	    (%%ecx),%%xmm3		/* t2 = t2*cx2 */\n\t"\
			"mulpd	0x10(%%edx),%%xmm4		/*~t4 = t4*sx3 */\n\t"\
			"mulpd	    (%%ebx),%%xmm1		/* t1 = t1*cx1 */\n\t"\
			"mulpd	0x10(%%ebx),%%xmm6		/*~t6 = t6*sx1 */\n\t"\
			"mulpd	    (%%eax),%%xmm0		/* ~r = r*(cx0-1) */\n\t"\
			"mulpd	0x10(%%ecx),%%xmm7		/*~t5 = t5*sx2 */\n\t"\
			"mulpd	    (%%edx),%%xmm2		/* t3 =(t1+t2)*cx3 */\n\t"\
			"addpd	0x10(%%edi),%%xmm0		/* t0 =~r + B0 */\n\t"\
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
			"movaps	0x080(%%edi),%%xmm1		/* t0r */\n\t"\
			"movaps	0x300(%%edi),%%xmm6		/* t3r */\n\t"\
			"subpd	%%xmm5,%%xmm1			/* B1r =t0r-t3i */\n\t"\
			"subpd	%%xmm6,%%xmm0			/* B6i =t0i-t3r */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*        2*t3i */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*        2*t3r */\n\t"\
			"addpd	%%xmm1,%%xmm5			/* B6r =t0r+t3i */\n\t"\
			"addpd	%%xmm0,%%xmm6			/* B1i =t0i+t3r */\n\t"\
			"movaps	%%xmm1,0x080(%%edi)		/* <-B1r = s1p08r */\n\t"\
			"movaps	%%xmm0,0x310(%%edi)		/* <-B6i = s1p20r */\n\t"\
			"movaps	%%xmm5,0x300(%%edi)		/* <-B6r = s1p20r */\n\t"\
			"movaps	%%xmm6,0x090(%%edi)		/* <-B1i = s1p08i */\n\t"\
			"\n\t"\
			"movaps	0x100(%%edi),%%xmm1		/* t1r */\n\t"\
			"movaps	0x280(%%edi),%%xmm6		/* t4r */\n\t"\
			"subpd	%%xmm7,%%xmm1			/* B2r =t1r-t4i */\n\t"\
			"subpd	%%xmm6,%%xmm2			/* B5i =t1i-t4r */\n\t"\
			"addpd	%%xmm7,%%xmm7			/*        2*t4i */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*        2*t4r */\n\t"\
			"addpd	%%xmm1,%%xmm7			/* B5r =t1r+t4i */\n\t"\
			"addpd	%%xmm2,%%xmm6			/* B2i =t1i+t4r */\n\t"\
			"movaps	%%xmm1,0x100(%%edi)		/* <-B2r = s1p16r */\n\t"\
			"movaps	%%xmm2,0x290(%%edi)		/* <-B5i = s1p12r */\n\t"\
			"movaps	%%xmm7,0x280(%%edi)		/* <-B5r = s1p12r */\n\t"\
			"movaps	%%xmm6,0x110(%%edi)		/* <-B2i = s1p16i */\n\t"\
			"\n\t"\
			"/* Note the order reversal on this pair of outputs: */\n\t"\
			"movaps	0x180(%%edi),%%xmm0		/* t2r */\n\t"\
			"movaps	0x200(%%edi),%%xmm5		/* t5r */\n\t"\
			"subpd	%%xmm4,%%xmm0			/* B4r =t2r-t5i */\n\t"\
			"subpd	%%xmm5,%%xmm3			/* B3i =t2i-t5r */\n\t"\
			"addpd	%%xmm4,%%xmm4			/*        2*t5i */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*        2*t5r */\n\t"\
			"addpd	%%xmm0,%%xmm4			/* B3r =t2r+t5i */\n\t"\
			"addpd	%%xmm3,%%xmm5			/* B4i =t2i+t5r */\n\t"\
			"movaps	%%xmm0,0x200(%%edi)		/* <-B4r = s1p04r */\n\t"\
			"movaps	%%xmm3,0x190(%%edi)		/* <-B3i = s1p24r */\n\t"\
			"movaps	%%xmm4,0x180(%%edi)		/* <-B3r = s1p24r */\n\t"\
			"movaps	%%xmm5,0x210(%%edi)		/* <-B4i = s1p04i */\n\t"\
			"\n\t"\
		"/* Since doing radix-7 in-place here, inputs of radix-4 are in consecutive memory locs, i.e. 0x20 bytes apart, e.g. the distance between s1p00r and s1p01r: */\n\t"\
			"\n\t"\
		"/* SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_E(s1p00r,s1p01r,s1p02r,s1p03r, add00+p[0,1,2,3]): */\n\t"\
			"movl	%[__out],%%esi	/* s1p00r */\n\t"\
			"movl	%[__add],%%eax\n\t"\
			"movl	%[__p01],%%ebx\n\t"\
			"movl	%[__p02],%%ecx\n\t"\
			"movl	%[__p03],%%edx\n\t"\
			"shll	$3,%%ebx		/* Pointer offset for floating doubles */\n\t"\
			"shll	$3,%%ecx\n\t"\
			"shll	$3,%%edx\n\t"\
			"addl	%%eax,%%ebx\n\t"\
			"addl	%%eax,%%ecx\n\t"\
			"addl	%%eax,%%edx\n\t"\
			"\n\t"\
			"movaps	    (%%esi),%%xmm0		/* in0r = s1p00r */\n\t"\
			"movaps	0x20(%%esi),%%xmm4		/* in1r = s1p01r */\n\t"\
			"movaps	0x10(%%esi),%%xmm1		/* in0i = s1p00i */\n\t"\
			"movaps	0x30(%%esi),%%xmm5		/* in1i = s1p01i */\n\t"\
			"movaps	%%xmm0,%%xmm2			/* xmm2 <- in0r */\n\t"\
			"movaps	%%xmm4,%%xmm6			/* xmm4 <- in1r */\n\t"\
			"movaps	%%xmm1,%%xmm3			/* xmm3 <- in0i */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* xmm5 <- in1i */\n\t"\
			"\n\t"\
			"addpd	0x40(%%esi),%%xmm0		/* t1 */\n\t"\
			"addpd	0x60(%%esi),%%xmm4		/* t5 */\n\t"\
			"addpd	0x50(%%esi),%%xmm1		/* t2 */\n\t"\
			"addpd	0x70(%%esi),%%xmm5		/* t6 */\n\t"\
			"subpd	0x40(%%esi),%%xmm2		/* t3 */\n\t"\
			"subpd	0x60(%%esi),%%xmm6		/* t7 */\n\t"\
			"subpd	0x50(%%esi),%%xmm3		/* t4 */\n\t"\
			"subpd	0x70(%%esi),%%xmm7		/* t8 */\n\t"\
			"/* Finish radix-4 butterfly and store results into output-array slots: */\n\t"\
			"subpd	%%xmm4,%%xmm0			/* ~t5 <- t1 -t5 */\n\t"\
			"subpd	%%xmm7,%%xmm2			/* ~t7 <- t3 -t8 */\n\t"\
			"subpd	%%xmm5,%%xmm1			/* ~t6 <- t2 -t6 */\n\t"\
			"subpd	%%xmm6,%%xmm3			/* ~t4 <- t4 -t7 */\n\t"\
			"movaps	%%xmm0,     (%%ebx)		/* <- ~t5 */\n\t"\
			"movaps	%%xmm2,     (%%ecx)		/* <- ~t7 */\n\t"\
			"movaps	%%xmm1,0x010(%%ebx)		/* <- ~t6 */\n\t"\
			"movaps	%%xmm3,0x010(%%edx)		/* <- ~t4 */\n\t"\
			"addpd	%%xmm4,%%xmm4			/*          2*t5 */\n\t"\
			"addpd	%%xmm7,%%xmm7			/*          2*t8 */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*          2*t6 */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*          2*t7 */\n\t"\
			"addpd	%%xmm0,%%xmm4			/* ~t1 <- t1 +t5 */\n\t"\
			"addpd	%%xmm2,%%xmm7			/* ~t3 <- t3 +t8 */\n\t"\
			"addpd	%%xmm1,%%xmm5			/* ~t2 <- t2 +t6 */\n\t"\
			"addpd	%%xmm3,%%xmm6			/* ~t8 <- t4 +t7 */\n\t"\
			"movaps	%%xmm4,     (%%eax)		/* <- ~t1 */\n\t"\
			"movaps	%%xmm7,     (%%edx)		/* <- ~t3 */\n\t"\
			"movaps	%%xmm5,0x010(%%eax)		/* <- ~t2 */\n\t"\
			"movaps	%%xmm6,0x010(%%ecx)		/* <- ~t8 */\n\t"\
		"/* SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_E(s1p04r,s1p05r,s1p06r,s1p07r, add24+p[1,0,3,2]): */\n\t"\
			"addl	$0x80,%%esi		/* s1p04r */\n\t"\
			"movl	%[__p24],%%edi\n\t"\
			"shll	$3,%%edi		/* Pointer offset for floating doubles */\n\t"\
			"subl	%%eax,%%ebx\n\t"\
			"subl	%%eax,%%ecx\n\t"\
			"subl	%%eax,%%edx\n\t"\
			"addl	%%edi,%%eax		/* add24 */\n\t"\
			"addl	%%eax,%%ebx\n\t"\
			"addl	%%eax,%%ecx\n\t"\
			"addl	%%eax,%%edx\n\t"\
			"\n\t"\
			"movaps	    (%%esi),%%xmm0		/* in0r = s1p00r */\n\t"\
			"movaps	0x20(%%esi),%%xmm4		/* in1r = s1p01r */\n\t"\
			"movaps	0x10(%%esi),%%xmm1		/* in0i = s1p00i */\n\t"\
			"movaps	0x30(%%esi),%%xmm5		/* in1i = s1p01i */\n\t"\
			"movaps	%%xmm0,%%xmm2			/* xmm2 <- in0r */\n\t"\
			"movaps	%%xmm4,%%xmm6			/* xmm4 <- in1r */\n\t"\
			"movaps	%%xmm1,%%xmm3			/* xmm3 <- in0i */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* xmm5 <- in1i */\n\t"\
			"\n\t"\
			"addpd	0x40(%%esi),%%xmm0		/* t1 */\n\t"\
			"addpd	0x60(%%esi),%%xmm4		/* t5 */\n\t"\
			"addpd	0x50(%%esi),%%xmm1		/* t2 */\n\t"\
			"addpd	0x70(%%esi),%%xmm5		/* t6 */\n\t"\
			"subpd	0x40(%%esi),%%xmm2		/* t3 */\n\t"\
			"subpd	0x60(%%esi),%%xmm6		/* t7 */\n\t"\
			"subpd	0x50(%%esi),%%xmm3		/* t4 */\n\t"\
			"subpd	0x70(%%esi),%%xmm7		/* t8 */\n\t"\
			"/* Finish radix-4 butterfly and store results into output-array slots: */\n\t"\
			"/* Swap _o1 <-> _o2, _o3 <-> _o4: eax <-> ebx, ecx <-> edx */\n\t"\
			"subpd	%%xmm4,%%xmm0			/* ~t5 <- t1 -t5 */\n\t"\
			"subpd	%%xmm7,%%xmm2			/* ~t7 <- t3 -t8 */\n\t"\
			"subpd	%%xmm5,%%xmm1			/* ~t6 <- t2 -t6 */\n\t"\
			"subpd	%%xmm6,%%xmm3			/* ~t4 <- t4 -t7 */\n\t"\
			"movaps	%%xmm0,     (%%eax)		/* <- ~t5 */\n\t"\
			"movaps	%%xmm2,     (%%edx)		/* <- ~t7 */\n\t"\
			"movaps	%%xmm1,0x010(%%eax)		/* <- ~t6 */\n\t"\
			"movaps	%%xmm3,0x010(%%ecx)		/* <- ~t4 */\n\t"\
			"addpd	%%xmm4,%%xmm4			/*          2*t5 */\n\t"\
			"addpd	%%xmm7,%%xmm7			/*          2*t8 */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*          2*t6 */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*          2*t7 */\n\t"\
			"addpd	%%xmm0,%%xmm4			/* ~t1 <- t1 +t5 */\n\t"\
			"addpd	%%xmm2,%%xmm7			/* ~t3 <- t3 +t8 */\n\t"\
			"addpd	%%xmm1,%%xmm5			/* ~t2 <- t2 +t6 */\n\t"\
			"addpd	%%xmm3,%%xmm6			/* ~t8 <- t4 +t7 */\n\t"\
			"movaps	%%xmm4,     (%%ebx)		/* <- ~t1 */\n\t"\
			"movaps	%%xmm7,     (%%ecx)		/* <- ~t3 */\n\t"\
			"movaps	%%xmm5,0x010(%%ebx)		/* <- ~t2 */\n\t"\
			"movaps	%%xmm6,0x010(%%edx)		/* <- ~t8 */\n\t"\
		"/* SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_E(s1p08r,s1p09r,s1p10r,s1p11r, add20+p[2,3,1,0]): */\n\t"\
			"addl	$0x80,%%esi		/* s1p08r */\n\t"\
			"subl	%%eax,%%ebx\n\t"\
			"subl	%%eax,%%ecx\n\t"\
			"subl	%%eax,%%edx\n\t"\
			"subl	%%edi,%%eax		/* add0  */\n\t"\
			"movl	%[__p20],%%edi\n\t"\
			"shll	$3,%%edi		/* Pointer offset for floating doubles */\n\t"\
			"addl	%%edi,%%eax		/* add20 */\n\t"\
			"addl	%%eax,%%ebx\n\t"\
			"addl	%%eax,%%ecx\n\t"\
			"addl	%%eax,%%edx\n\t"\
			"\n\t"\
			"movaps	    (%%esi),%%xmm0		/* in0r = s1p00r */\n\t"\
			"movaps	0x20(%%esi),%%xmm4		/* in1r = s1p01r */\n\t"\
			"movaps	0x10(%%esi),%%xmm1		/* in0i = s1p00i */\n\t"\
			"movaps	0x30(%%esi),%%xmm5		/* in1i = s1p01i */\n\t"\
			"movaps	%%xmm0,%%xmm2			/* xmm2 <- in0r */\n\t"\
			"movaps	%%xmm4,%%xmm6			/* xmm4 <- in1r */\n\t"\
			"movaps	%%xmm1,%%xmm3			/* xmm3 <- in0i */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* xmm5 <- in1i */\n\t"\
			"\n\t"\
			"addpd	0x40(%%esi),%%xmm0		/* t1 */\n\t"\
			"addpd	0x60(%%esi),%%xmm4		/* t5 */\n\t"\
			"addpd	0x50(%%esi),%%xmm1		/* t2 */\n\t"\
			"addpd	0x70(%%esi),%%xmm5		/* t6 */\n\t"\
			"subpd	0x40(%%esi),%%xmm2		/* t3 */\n\t"\
			"subpd	0x60(%%esi),%%xmm6		/* t7 */\n\t"\
			"subpd	0x50(%%esi),%%xmm3		/* t4 */\n\t"\
			"subpd	0x70(%%esi),%%xmm7		/* t8 */\n\t"\
			"/* Finish radix-4 butterfly and store results into output-array slots: */\n\t"\
			"/* Swap _o0123 <-> _o2310: e[abcd]x <-> e[cdba]x */\n\t"\
			"subpd	%%xmm4,%%xmm0			/* ~t5 <- t1 -t5 */\n\t"\
			"subpd	%%xmm7,%%xmm2			/* ~t7 <- t3 -t8 */\n\t"\
			"subpd	%%xmm5,%%xmm1			/* ~t6 <- t2 -t6 */\n\t"\
			"subpd	%%xmm6,%%xmm3			/* ~t4 <- t4 -t7 */\n\t"\
			"movaps	%%xmm0,     (%%edx)		/* <- ~t5 */\n\t"\
			"movaps	%%xmm2,     (%%ebx)		/* <- ~t7 */\n\t"\
			"movaps	%%xmm1,0x010(%%edx)		/* <- ~t6 */\n\t"\
			"movaps	%%xmm3,0x010(%%eax)		/* <- ~t4 */\n\t"\
			"addpd	%%xmm4,%%xmm4			/*          2*t5 */\n\t"\
			"addpd	%%xmm7,%%xmm7			/*          2*t8 */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*          2*t6 */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*          2*t7 */\n\t"\
			"addpd	%%xmm0,%%xmm4			/* ~t1 <- t1 +t5 */\n\t"\
			"addpd	%%xmm2,%%xmm7			/* ~t3 <- t3 +t8 */\n\t"\
			"addpd	%%xmm1,%%xmm5			/* ~t2 <- t2 +t6 */\n\t"\
			"addpd	%%xmm3,%%xmm6			/* ~t8 <- t4 +t7 */\n\t"\
			"movaps	%%xmm4,     (%%ecx)		/* <- ~t1 */\n\t"\
			"movaps	%%xmm7,     (%%eax)		/* <- ~t3 */\n\t"\
			"movaps	%%xmm5,0x010(%%ecx)		/* <- ~t2 */\n\t"\
			"movaps	%%xmm6,0x010(%%ebx)		/* <- ~t8 */\n\t"\
		"/* SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_E(s1p12r,s1p13r,s1p14r,s1p15r, add16+p[0,1,2,3]): */\n\t"\
			"addl	$0x80,%%esi		/* s1p12r */\n\t"\
			"subl	%%eax,%%ebx\n\t"\
			"subl	%%eax,%%ecx\n\t"\
			"subl	%%eax,%%edx\n\t"\
			"subl	%%edi,%%eax		/* add0  */\n\t"\
			"movl	%[__p16],%%edi\n\t"\
			"shll	$3,%%edi		/* Pointer offset for floating doubles */\n\t"\
			"addl	%%edi,%%eax		/* add16 */\n\t"\
			"addl	%%eax,%%ebx\n\t"\
			"addl	%%eax,%%ecx\n\t"\
			"addl	%%eax,%%edx\n\t"\
			"\n\t"\
			"movaps	    (%%esi),%%xmm0		/* in0r = s1p00r */\n\t"\
			"movaps	0x20(%%esi),%%xmm4		/* in1r = s1p01r */\n\t"\
			"movaps	0x10(%%esi),%%xmm1		/* in0i = s1p00i */\n\t"\
			"movaps	0x30(%%esi),%%xmm5		/* in1i = s1p01i */\n\t"\
			"movaps	%%xmm0,%%xmm2			/* xmm2 <- in0r */\n\t"\
			"movaps	%%xmm4,%%xmm6			/* xmm4 <- in1r */\n\t"\
			"movaps	%%xmm1,%%xmm3			/* xmm3 <- in0i */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* xmm5 <- in1i */\n\t"\
			"\n\t"\
			"addpd	0x40(%%esi),%%xmm0		/* t1 */\n\t"\
			"addpd	0x60(%%esi),%%xmm4		/* t5 */\n\t"\
			"addpd	0x50(%%esi),%%xmm1		/* t2 */\n\t"\
			"addpd	0x70(%%esi),%%xmm5		/* t6 */\n\t"\
			"subpd	0x40(%%esi),%%xmm2		/* t3 */\n\t"\
			"subpd	0x60(%%esi),%%xmm6		/* t7 */\n\t"\
			"subpd	0x50(%%esi),%%xmm3		/* t4 */\n\t"\
			"subpd	0x70(%%esi),%%xmm7		/* t8 */\n\t"\
			"/* Finish radix-4 butterfly and store results into output-array slots: */\n\t"\
			"subpd	%%xmm4,%%xmm0			/* ~t5 <- t1 -t5 */\n\t"\
			"subpd	%%xmm7,%%xmm2			/* ~t7 <- t3 -t8 */\n\t"\
			"subpd	%%xmm5,%%xmm1			/* ~t6 <- t2 -t6 */\n\t"\
			"subpd	%%xmm6,%%xmm3			/* ~t4 <- t4 -t7 */\n\t"\
			"movaps	%%xmm0,     (%%ebx)		/* <- ~t5 */\n\t"\
			"movaps	%%xmm2,     (%%ecx)		/* <- ~t7 */\n\t"\
			"movaps	%%xmm1,0x010(%%ebx)		/* <- ~t6 */\n\t"\
			"movaps	%%xmm3,0x010(%%edx)		/* <- ~t4 */\n\t"\
			"addpd	%%xmm4,%%xmm4			/*          2*t5 */\n\t"\
			"addpd	%%xmm7,%%xmm7			/*          2*t8 */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*          2*t6 */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*          2*t7 */\n\t"\
			"addpd	%%xmm0,%%xmm4			/* ~t1 <- t1 +t5 */\n\t"\
			"addpd	%%xmm2,%%xmm7			/* ~t3 <- t3 +t8 */\n\t"\
			"addpd	%%xmm1,%%xmm5			/* ~t2 <- t2 +t6 */\n\t"\
			"addpd	%%xmm3,%%xmm6			/* ~t8 <- t4 +t7 */\n\t"\
			"movaps	%%xmm4,     (%%eax)		/* <- ~t1 */\n\t"\
			"movaps	%%xmm7,     (%%edx)		/* <- ~t3 */\n\t"\
			"movaps	%%xmm5,0x010(%%eax)		/* <- ~t2 */\n\t"\
			"movaps	%%xmm6,0x010(%%ecx)		/* <- ~t8 */\n\t"\
		"/* SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_E(s1p16r,s1p17r,s1p18r,s1p19r, add12+p[3,2,0,1]): */\n\t"\
			"addl	$0x80,%%esi		/* s1p16r */\n\t"\
			"subl	%%eax,%%ebx\n\t"\
			"subl	%%eax,%%ecx\n\t"\
			"subl	%%eax,%%edx\n\t"\
			"subl	%%edi,%%eax		/* add0  */\n\t"\
			"movl	%[__p12],%%edi\n\t"\
			"shll	$3,%%edi		/* Pointer offset for floating doubles */\n\t"\
			"addl	%%edi,%%eax		/* add12 */\n\t"\
			"addl	%%eax,%%ebx\n\t"\
			"addl	%%eax,%%ecx\n\t"\
			"addl	%%eax,%%edx\n\t"\
			"\n\t"\
			"movaps	    (%%esi),%%xmm0		/* in0r = s1p00r */\n\t"\
			"movaps	0x20(%%esi),%%xmm4		/* in1r = s1p01r */\n\t"\
			"movaps	0x10(%%esi),%%xmm1		/* in0i = s1p00i */\n\t"\
			"movaps	0x30(%%esi),%%xmm5		/* in1i = s1p01i */\n\t"\
			"movaps	%%xmm0,%%xmm2			/* xmm2 <- in0r */\n\t"\
			"movaps	%%xmm4,%%xmm6			/* xmm4 <- in1r */\n\t"\
			"movaps	%%xmm1,%%xmm3			/* xmm3 <- in0i */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* xmm5 <- in1i */\n\t"\
			"\n\t"\
			"addpd	0x40(%%esi),%%xmm0		/* t1 */\n\t"\
			"addpd	0x60(%%esi),%%xmm4		/* t5 */\n\t"\
			"addpd	0x50(%%esi),%%xmm1		/* t2 */\n\t"\
			"addpd	0x70(%%esi),%%xmm5		/* t6 */\n\t"\
			"subpd	0x40(%%esi),%%xmm2		/* t3 */\n\t"\
			"subpd	0x60(%%esi),%%xmm6		/* t7 */\n\t"\
			"subpd	0x50(%%esi),%%xmm3		/* t4 */\n\t"\
			"subpd	0x70(%%esi),%%xmm7		/* t8 */\n\t"\
			"/* Finish radix-4 butterfly and store results into output-array slots: */\n\t"\
			"/* Swap _o0123 <-> _o3201: e[abcd]x <-> e[dcab]x */\n\t"\
			"subpd	%%xmm4,%%xmm0			/* ~t5 <- t1 -t5 */\n\t"\
			"subpd	%%xmm7,%%xmm2			/* ~t7 <- t3 -t8 */\n\t"\
			"subpd	%%xmm5,%%xmm1			/* ~t6 <- t2 -t6 */\n\t"\
			"subpd	%%xmm6,%%xmm3			/* ~t4 <- t4 -t7 */\n\t"\
			"movaps	%%xmm0,     (%%ecx)		/* <- ~t5 */\n\t"\
			"movaps	%%xmm2,     (%%eax)		/* <- ~t7 */\n\t"\
			"movaps	%%xmm1,0x010(%%ecx)		/* <- ~t6 */\n\t"\
			"movaps	%%xmm3,0x010(%%ebx)		/* <- ~t4 */\n\t"\
			"addpd	%%xmm4,%%xmm4			/*          2*t5 */\n\t"\
			"addpd	%%xmm7,%%xmm7			/*          2*t8 */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*          2*t6 */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*          2*t7 */\n\t"\
			"addpd	%%xmm0,%%xmm4			/* ~t1 <- t1 +t5 */\n\t"\
			"addpd	%%xmm2,%%xmm7			/* ~t3 <- t3 +t8 */\n\t"\
			"addpd	%%xmm1,%%xmm5			/* ~t2 <- t2 +t6 */\n\t"\
			"addpd	%%xmm3,%%xmm6			/* ~t8 <- t4 +t7 */\n\t"\
			"movaps	%%xmm4,     (%%edx)		/* <- ~t1 */\n\t"\
			"movaps	%%xmm7,     (%%ebx)		/* <- ~t3 */\n\t"\
			"movaps	%%xmm5,0x010(%%edx)		/* <- ~t2 */\n\t"\
			"movaps	%%xmm6,0x010(%%eax)		/* <- ~t8 */\n\t"\
		"/* SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_E(s1p20r,s1p21r,s1p22r,s1p23r, add08+p[1,0,3,2]): */\n\t"\
			"addl	$0x80,%%esi		/* s1p20r */\n\t"\
			"subl	%%eax,%%ebx\n\t"\
			"subl	%%eax,%%ecx\n\t"\
			"subl	%%eax,%%edx\n\t"\
			"subl	%%edi,%%eax		/* add0  */\n\t"\
			"movl	%[__p08],%%edi\n\t"\
			"shll	$3,%%edi		/* Pointer offset for floating doubles */\n\t"\
			"addl	%%edi,%%eax		/* add08 */\n\t"\
			"addl	%%eax,%%ebx\n\t"\
			"addl	%%eax,%%ecx\n\t"\
			"addl	%%eax,%%edx\n\t"\
			"\n\t"\
			"movaps	    (%%esi),%%xmm0		/* in0r = s1p00r */\n\t"\
			"movaps	0x20(%%esi),%%xmm4		/* in1r = s1p01r */\n\t"\
			"movaps	0x10(%%esi),%%xmm1		/* in0i = s1p00i */\n\t"\
			"movaps	0x30(%%esi),%%xmm5		/* in1i = s1p01i */\n\t"\
			"movaps	%%xmm0,%%xmm2			/* xmm2 <- in0r */\n\t"\
			"movaps	%%xmm4,%%xmm6			/* xmm4 <- in1r */\n\t"\
			"movaps	%%xmm1,%%xmm3			/* xmm3 <- in0i */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* xmm5 <- in1i */\n\t"\
			"\n\t"\
			"addpd	0x40(%%esi),%%xmm0		/* t1 */\n\t"\
			"addpd	0x60(%%esi),%%xmm4		/* t5 */\n\t"\
			"addpd	0x50(%%esi),%%xmm1		/* t2 */\n\t"\
			"addpd	0x70(%%esi),%%xmm5		/* t6 */\n\t"\
			"subpd	0x40(%%esi),%%xmm2		/* t3 */\n\t"\
			"subpd	0x60(%%esi),%%xmm6		/* t7 */\n\t"\
			"subpd	0x50(%%esi),%%xmm3		/* t4 */\n\t"\
			"subpd	0x70(%%esi),%%xmm7		/* t8 */\n\t"\
			"/* Finish radix-4 butterfly and store results into output-array slots: */\n\t"\
			"/* Swap _o0 <-> _o1, _o2 <-> _o3: eax <-> ebx, ecx <-> edx */\n\t"\
			"subpd	%%xmm4,%%xmm0			/* ~t5 <- t1 -t5 */\n\t"\
			"subpd	%%xmm7,%%xmm2			/* ~t7 <- t3 -t8 */\n\t"\
			"subpd	%%xmm5,%%xmm1			/* ~t6 <- t2 -t6 */\n\t"\
			"subpd	%%xmm6,%%xmm3			/* ~t4 <- t4 -t7 */\n\t"\
			"movaps	%%xmm0,     (%%eax)		/* <- ~t5 */\n\t"\
			"movaps	%%xmm2,     (%%edx)		/* <- ~t7 */\n\t"\
			"movaps	%%xmm1,0x010(%%eax)		/* <- ~t6 */\n\t"\
			"movaps	%%xmm3,0x010(%%ecx)		/* <- ~t4 */\n\t"\
			"addpd	%%xmm4,%%xmm4			/*          2*t5 */\n\t"\
			"addpd	%%xmm7,%%xmm7			/*          2*t8 */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*          2*t6 */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*          2*t7 */\n\t"\
			"addpd	%%xmm0,%%xmm4			/* ~t1 <- t1 +t5 */\n\t"\
			"addpd	%%xmm2,%%xmm7			/* ~t3 <- t3 +t8 */\n\t"\
			"addpd	%%xmm1,%%xmm5			/* ~t2 <- t2 +t6 */\n\t"\
			"addpd	%%xmm3,%%xmm6			/* ~t8 <- t4 +t7 */\n\t"\
			"movaps	%%xmm4,     (%%ebx)		/* <- ~t1 */\n\t"\
			"movaps	%%xmm7,     (%%ecx)		/* <- ~t3 */\n\t"\
			"movaps	%%xmm5,0x010(%%ebx)		/* <- ~t2 */\n\t"\
			"movaps	%%xmm6,0x010(%%edx)		/* <- ~t8 */\n\t"\
		"/* SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_E(s1p24r,s1p25r,s1p26r,s1p27r, add04+p[2,3,1,0]): */\n\t"\
			"addl	$0x80,%%esi		/* s1p24r */\n\t"\
			"subl	%%eax,%%ebx\n\t"\
			"subl	%%eax,%%ecx\n\t"\
			"subl	%%eax,%%edx\n\t"\
			"subl	%%edi,%%eax		/* add0  */\n\t"\
			"movl	%[__p04],%%edi\n\t"\
			"shll	$3,%%edi		/* Pointer offset for floating doubles */\n\t"\
			"addl	%%edi,%%eax		/* add04 */\n\t"\
			"addl	%%eax,%%ebx\n\t"\
			"addl	%%eax,%%ecx\n\t"\
			"addl	%%eax,%%edx\n\t"\
			"\n\t"\
			"movaps	    (%%esi),%%xmm0		/* in0r = s1p00r */\n\t"\
			"movaps	0x20(%%esi),%%xmm4		/* in1r = s1p01r */\n\t"\
			"movaps	0x10(%%esi),%%xmm1		/* in0i = s1p00i */\n\t"\
			"movaps	0x30(%%esi),%%xmm5		/* in1i = s1p01i */\n\t"\
			"movaps	%%xmm0,%%xmm2			/* xmm2 <- in0r */\n\t"\
			"movaps	%%xmm4,%%xmm6			/* xmm4 <- in1r */\n\t"\
			"movaps	%%xmm1,%%xmm3			/* xmm3 <- in0i */\n\t"\
			"movaps	%%xmm5,%%xmm7			/* xmm5 <- in1i */\n\t"\
			"\n\t"\
			"addpd	0x40(%%esi),%%xmm0		/* t1 */\n\t"\
			"addpd	0x60(%%esi),%%xmm4		/* t5 */\n\t"\
			"addpd	0x50(%%esi),%%xmm1		/* t2 */\n\t"\
			"addpd	0x70(%%esi),%%xmm5		/* t6 */\n\t"\
			"subpd	0x40(%%esi),%%xmm2		/* t3 */\n\t"\
			"subpd	0x60(%%esi),%%xmm6		/* t7 */\n\t"\
			"subpd	0x50(%%esi),%%xmm3		/* t4 */\n\t"\
			"subpd	0x70(%%esi),%%xmm7		/* t8 */\n\t"\
			"/* Finish radix-4 butterfly and store results into output-array slots: */\n\t"\
			"/* Swap _o0123 <-> _o2310: e[abcd]x <-> e[cdba]x */\n\t"\
			"subpd	%%xmm4,%%xmm0			/* ~t5 <- t1 -t5 */\n\t"\
			"subpd	%%xmm7,%%xmm2			/* ~t7 <- t3 -t8 */\n\t"\
			"subpd	%%xmm5,%%xmm1			/* ~t6 <- t2 -t6 */\n\t"\
			"subpd	%%xmm6,%%xmm3			/* ~t4 <- t4 -t7 */\n\t"\
			"movaps	%%xmm0,     (%%edx)		/* <- ~t5 */\n\t"\
			"movaps	%%xmm2,     (%%ebx)		/* <- ~t7 */\n\t"\
			"movaps	%%xmm1,0x010(%%edx)		/* <- ~t6 */\n\t"\
			"movaps	%%xmm3,0x010(%%eax)		/* <- ~t4 */\n\t"\
			"addpd	%%xmm4,%%xmm4			/*          2*t5 */\n\t"\
			"addpd	%%xmm7,%%xmm7			/*          2*t8 */\n\t"\
			"addpd	%%xmm5,%%xmm5			/*          2*t6 */\n\t"\
			"addpd	%%xmm6,%%xmm6			/*          2*t7 */\n\t"\
			"addpd	%%xmm0,%%xmm4			/* ~t1 <- t1 +t5 */\n\t"\
			"addpd	%%xmm2,%%xmm7			/* ~t3 <- t3 +t8 */\n\t"\
			"addpd	%%xmm1,%%xmm5			/* ~t2 <- t2 +t6 */\n\t"\
			"addpd	%%xmm3,%%xmm6			/* ~t8 <- t4 +t7 */\n\t"\
			"movaps	%%xmm4,     (%%ecx)		/* <- ~t1 */\n\t"\
			"movaps	%%xmm7,     (%%eax)		/* <- ~t3 */\n\t"\
			"movaps	%%xmm5,0x010(%%ecx)		/* <- ~t2 */\n\t"\
			"movaps	%%xmm6,0x010(%%ebx)		/* <- ~t8 */\n\t"\
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
			: "eax","ebx","ecx","edx","edi","esi"		/* Clobbered registers */\
		);\
	}

#endif	/* radix28_ditN_cy_dif1_gcc_h_included */
