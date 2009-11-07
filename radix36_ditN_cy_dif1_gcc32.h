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
			"movl	%[__r00],%%esi\n\t"\
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
			"movaps	     (%%eax),%%xmm2	\n\t"\
			"movaps	     (%%edx),%%xmm6	\n\t"\
			"movaps	0x010(%%eax),%%xmm3	\n\t"\
			"movaps	0x010(%%edx),%%xmm7	\n\t"\
			"movaps	     (%%ebx),%%xmm0	\n\t"\
			"movaps	     (%%ecx),%%xmm4	\n\t"\
			"movaps	0x010(%%ebx),%%xmm1	\n\t"\
			"movaps	0x010(%%ecx),%%xmm5	\n\t"\
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
			"movaps	%%xmm0,0x240(%%esi)	\n\t"\
			"movaps	%%xmm2,0x360(%%esi)	\n\t"\
			"movaps	%%xmm1,0x250(%%esi)	\n\t"\
			"movaps	%%xmm3,0x130(%%esi)	\n\t"\
			"addpd	%%xmm4,%%xmm4		\n\t"\
			"addpd	%%xmm7,%%xmm7		\n\t"\
			"addpd	%%xmm5,%%xmm5		\n\t"\
			"addpd	%%xmm6,%%xmm6		\n\t"\
			"addpd	%%xmm0,%%xmm4		\n\t"\
			"addpd	%%xmm2,%%xmm7		\n\t"\
			"addpd	%%xmm1,%%xmm5		\n\t"\
			"addpd	%%xmm3,%%xmm6		\n\t"\
			"movaps	%%xmm4,     (%%esi)	\n\t"\
			"movaps	%%xmm7,0x120(%%esi)	\n\t"\
			"movaps	%%xmm5,0x010(%%esi)	\n\t"\
			"movaps	%%xmm6,0x370(%%esi)	\n\t"\
		"/*	add2,3,0,1 = &a[j1+p04]+p0,1,2,3: SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_C(ecx,edx,eax,ebx, 0x120, 0x240, r04) */\n\t"\
			"addl	$0x040,%%esi	/* r04 */\n\t"\
			"movl	%[__p04],%%edi\n\t"\
			"shll	$3,%%edi \n\t"\
			"addl	%%edi,%%eax	/* add04 */\n\t"\
			"addl	%%edi,%%ebx\n\t"\
			"addl	%%edi,%%ecx\n\t"\
			"addl	%%edi,%%edx\n\t"\
			"/* e[ab]x <-> e[cd]x */\n\t"\
			"movaps	     (%%ecx),%%xmm2	\n\t"\
			"movaps	     (%%eax),%%xmm6	\n\t"\
			"movaps	0x010(%%ecx),%%xmm3	\n\t"\
			"movaps	0x010(%%eax),%%xmm7	\n\t"\
			"movaps	     (%%edx),%%xmm0	\n\t"\
			"movaps	     (%%ebx),%%xmm4	\n\t"\
			"movaps	0x010(%%edx),%%xmm1	\n\t"\
			"movaps	0x010(%%ebx),%%xmm5	\n\t"\
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
			"movaps	%%xmm0,0x240(%%esi)	\n\t"\
			"movaps	%%xmm2,0x360(%%esi)	\n\t"\
			"movaps	%%xmm1,0x250(%%esi)	\n\t"\
			"movaps	%%xmm3,0x130(%%esi)	\n\t"\
			"addpd	%%xmm4,%%xmm4		\n\t"\
			"addpd	%%xmm7,%%xmm7		\n\t"\
			"addpd	%%xmm5,%%xmm5		\n\t"\
			"addpd	%%xmm6,%%xmm6		\n\t"\
			"addpd	%%xmm0,%%xmm4		\n\t"\
			"addpd	%%xmm2,%%xmm7		\n\t"\
			"addpd	%%xmm1,%%xmm5		\n\t"\
			"addpd	%%xmm3,%%xmm6		\n\t"\
			"movaps	%%xmm4,     (%%esi)	\n\t"\
			"movaps	%%xmm7,0x120(%%esi)	\n\t"\
			"movaps	%%xmm5,0x010(%%esi)	\n\t"\
			"movaps	%%xmm6,0x370(%%esi)	\n\t"\
		"/*	add1,0,2,3 = &a[j1+p08]+p0,1,2,3: SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_C(ebx,eax,ecx,edx, 0x120, 0x240, r02) */\n\t"\
			"subl	$0x020,%%esi	/* r02 */\n\t"\
			"addl	%%edi,%%eax	/* add08 */\n\t"\
			"addl	%%edi,%%ebx\n\t"\
			"addl	%%edi,%%ecx\n\t"\
			"addl	%%edi,%%edx\n\t"\
			"/* eax <-> ebx */\n\t"\
			"movaps	     (%%ebx),%%xmm2	\n\t"\
			"movaps	     (%%ecx),%%xmm6	\n\t"\
			"movaps	0x010(%%ebx),%%xmm3	\n\t"\
			"movaps	0x010(%%ecx),%%xmm7	\n\t"\
			"movaps	     (%%eax),%%xmm0	\n\t"\
			"movaps	     (%%edx),%%xmm4	\n\t"\
			"movaps	0x010(%%eax),%%xmm1	\n\t"\
			"movaps	0x010(%%edx),%%xmm5	\n\t"\
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
			"movaps	%%xmm0,0x240(%%esi)	\n\t"\
			"movaps	%%xmm2,0x360(%%esi)	\n\t"\
			"movaps	%%xmm1,0x250(%%esi)	\n\t"\
			"movaps	%%xmm3,0x130(%%esi)	\n\t"\
			"addpd	%%xmm4,%%xmm4		\n\t"\
			"addpd	%%xmm7,%%xmm7		\n\t"\
			"addpd	%%xmm5,%%xmm5		\n\t"\
			"addpd	%%xmm6,%%xmm6		\n\t"\
			"addpd	%%xmm0,%%xmm4		\n\t"\
			"addpd	%%xmm2,%%xmm7		\n\t"\
			"addpd	%%xmm1,%%xmm5		\n\t"\
			"addpd	%%xmm3,%%xmm6		\n\t"\
			"movaps	%%xmm4,     (%%esi)	\n\t"\
			"movaps	%%xmm7,0x120(%%esi)	\n\t"\
			"movaps	%%xmm5,0x010(%%esi)	\n\t"\
			"movaps	%%xmm6,0x370(%%esi)	\n\t"\
		"/*	add3,2,1,0 = &a[j1+p12]+p0,1,2,3: SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_C(edx,ecx,ebx,eax, 0x120, 0x240, r0c) */\n\t"\
			"addl	$0x0a0,%%esi	/* r0c */\n\t"\
			"addl	%%edi,%%eax	/* add12 */\n\t"\
			"addl	%%edi,%%ebx\n\t"\
			"addl	%%edi,%%ecx\n\t"\
			"addl	%%edi,%%edx\n\t"\
			"/* eax <-> edx, ebx <-> ecx */\n\t"\
			"movaps	     (%%edx),%%xmm2	\n\t"\
			"movaps	     (%%ebx),%%xmm6	\n\t"\
			"movaps	0x010(%%edx),%%xmm3	\n\t"\
			"movaps	0x010(%%ebx),%%xmm7	\n\t"\
			"movaps	     (%%ecx),%%xmm0	\n\t"\
			"movaps	     (%%eax),%%xmm4	\n\t"\
			"movaps	0x010(%%ecx),%%xmm1	\n\t"\
			"movaps	0x010(%%eax),%%xmm5	\n\t"\
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
			"movaps	%%xmm0,0x240(%%esi)	\n\t"\
			"movaps	%%xmm2,0x360(%%esi)	\n\t"\
			"movaps	%%xmm1,0x250(%%esi)	\n\t"\
			"movaps	%%xmm3,0x130(%%esi)	\n\t"\
			"addpd	%%xmm4,%%xmm4		\n\t"\
			"addpd	%%xmm7,%%xmm7		\n\t"\
			"addpd	%%xmm5,%%xmm5		\n\t"\
			"addpd	%%xmm6,%%xmm6		\n\t"\
			"addpd	%%xmm0,%%xmm4		\n\t"\
			"addpd	%%xmm2,%%xmm7		\n\t"\
			"addpd	%%xmm1,%%xmm5		\n\t"\
			"addpd	%%xmm3,%%xmm6		\n\t"\
			"movaps	%%xmm4,     (%%esi)	\n\t"\
			"movaps	%%xmm7,0x120(%%esi)	\n\t"\
			"movaps	%%xmm5,0x010(%%esi)	\n\t"\
			"movaps	%%xmm6,0x370(%%esi)	\n\t"\
		"/*	add0,1,3,2 = &a[j1+p16]+p0,1,2,3: SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_C(eax,ebx,edx,ecx, 0x120, 0x240, r0g) */\n\t"\
			"addl	$0x040,%%esi	/* r0g */\n\t"\
			"addl	%%edi,%%eax	/* add16 */\n\t"\
			"addl	%%edi,%%ebx\n\t"\
			"addl	%%edi,%%ecx\n\t"\
			"addl	%%edi,%%edx\n\t"\
			"/* ecx <-> edx */\n\t"\
			"movaps	     (%%eax),%%xmm2	\n\t"\
			"movaps	     (%%edx),%%xmm6	\n\t"\
			"movaps	0x010(%%eax),%%xmm3	\n\t"\
			"movaps	0x010(%%edx),%%xmm7	\n\t"\
			"movaps	     (%%ebx),%%xmm0	\n\t"\
			"movaps	     (%%ecx),%%xmm4	\n\t"\
			"movaps	0x010(%%ebx),%%xmm1	\n\t"\
			"movaps	0x010(%%ecx),%%xmm5	\n\t"\
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
			"movaps	%%xmm0,0x240(%%esi)	\n\t"\
			"movaps	%%xmm2,0x360(%%esi)	\n\t"\
			"movaps	%%xmm1,0x250(%%esi)	\n\t"\
			"movaps	%%xmm3,0x130(%%esi)	\n\t"\
			"addpd	%%xmm4,%%xmm4		\n\t"\
			"addpd	%%xmm7,%%xmm7		\n\t"\
			"addpd	%%xmm5,%%xmm5		\n\t"\
			"addpd	%%xmm6,%%xmm6		\n\t"\
			"addpd	%%xmm0,%%xmm4		\n\t"\
			"addpd	%%xmm2,%%xmm7		\n\t"\
			"addpd	%%xmm1,%%xmm5		\n\t"\
			"addpd	%%xmm3,%%xmm6		\n\t"\
			"movaps	%%xmm4,     (%%esi)	\n\t"\
			"movaps	%%xmm7,0x120(%%esi)	\n\t"\
			"movaps	%%xmm5,0x010(%%esi)	\n\t"\
			"movaps	%%xmm6,0x370(%%esi)	\n\t"\
		"/*	add2,3,0,1 = &a[j1+p20]+p0,1,2,3: SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_C(ecx,edx,eax,ebx, 0x120, 0x240, r0e) */\n\t"\
			"subl	$0x020,%%esi	/* r0e */\n\t"\
			"addl	%%edi,%%eax	/* add20 */\n\t"\
			"addl	%%edi,%%ebx\n\t"\
			"addl	%%edi,%%ecx\n\t"\
			"addl	%%edi,%%edx\n\t"\
			"/* e[ab]x <-> e[cd]x */\n\t"\
			"movaps	     (%%ecx),%%xmm2	\n\t"\
			"movaps	     (%%eax),%%xmm6	\n\t"\
			"movaps	0x010(%%ecx),%%xmm3	\n\t"\
			"movaps	0x010(%%eax),%%xmm7	\n\t"\
			"movaps	     (%%edx),%%xmm0	\n\t"\
			"movaps	     (%%ebx),%%xmm4	\n\t"\
			"movaps	0x010(%%edx),%%xmm1	\n\t"\
			"movaps	0x010(%%ebx),%%xmm5	\n\t"\
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
			"movaps	%%xmm0,0x240(%%esi)	\n\t"\
			"movaps	%%xmm2,0x360(%%esi)	\n\t"\
			"movaps	%%xmm1,0x250(%%esi)	\n\t"\
			"movaps	%%xmm3,0x130(%%esi)	\n\t"\
			"addpd	%%xmm4,%%xmm4		\n\t"\
			"addpd	%%xmm7,%%xmm7		\n\t"\
			"addpd	%%xmm5,%%xmm5		\n\t"\
			"addpd	%%xmm6,%%xmm6		\n\t"\
			"addpd	%%xmm0,%%xmm4		\n\t"\
			"addpd	%%xmm2,%%xmm7		\n\t"\
			"addpd	%%xmm1,%%xmm5		\n\t"\
			"addpd	%%xmm3,%%xmm6		\n\t"\
			"movaps	%%xmm4,     (%%esi)	\n\t"\
			"movaps	%%xmm7,0x120(%%esi)	\n\t"\
			"movaps	%%xmm5,0x010(%%esi)	\n\t"\
			"movaps	%%xmm6,0x370(%%esi)	\n\t"\
		"/*	add1,0,2,3 = &a[j1+p24]+p0,1,2,3: SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_C(ebx,eax,ecx,edx, 0x120, 0x240, r08) */\n\t"\
			"subl	$0x060,%%esi	/* r08 */\n\t"\
			"addl	%%edi,%%eax	/* add24 */\n\t"\
			"addl	%%edi,%%ebx\n\t"\
			"addl	%%edi,%%ecx\n\t"\
			"addl	%%edi,%%edx\n\t"\
			"/* eax <-> ebx */\n\t"\
			"movaps	     (%%ebx),%%xmm2	\n\t"\
			"movaps	     (%%ecx),%%xmm6	\n\t"\
			"movaps	0x010(%%ebx),%%xmm3	\n\t"\
			"movaps	0x010(%%ecx),%%xmm7	\n\t"\
			"movaps	     (%%eax),%%xmm0	\n\t"\
			"movaps	     (%%edx),%%xmm4	\n\t"\
			"movaps	0x010(%%eax),%%xmm1	\n\t"\
			"movaps	0x010(%%edx),%%xmm5	\n\t"\
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
			"movaps	%%xmm0,0x240(%%esi)	\n\t"\
			"movaps	%%xmm2,0x360(%%esi)	\n\t"\
			"movaps	%%xmm1,0x250(%%esi)	\n\t"\
			"movaps	%%xmm3,0x130(%%esi)	\n\t"\
			"addpd	%%xmm4,%%xmm4		\n\t"\
			"addpd	%%xmm7,%%xmm7		\n\t"\
			"addpd	%%xmm5,%%xmm5		\n\t"\
			"addpd	%%xmm6,%%xmm6		\n\t"\
			"addpd	%%xmm0,%%xmm4		\n\t"\
			"addpd	%%xmm2,%%xmm7		\n\t"\
			"addpd	%%xmm1,%%xmm5		\n\t"\
			"addpd	%%xmm3,%%xmm6		\n\t"\
			"movaps	%%xmm4,     (%%esi)	\n\t"\
			"movaps	%%xmm7,0x120(%%esi)	\n\t"\
			"movaps	%%xmm5,0x010(%%esi)	\n\t"\
			"movaps	%%xmm6,0x370(%%esi)	\n\t"\
		"/*	add3,2,1,0 = &a[j1+p28]+p0,1,2,3: SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_C(edx,ecx,ebx,eax, 0x120, 0x240, r06) */\n\t"\
			"subl	$0x020,%%esi	/* r04 */\n\t"\
			"addl	%%edi,%%eax	/* add28 */\n\t"\
			"addl	%%edi,%%ebx\n\t"\
			"addl	%%edi,%%ecx\n\t"\
			"addl	%%edi,%%edx\n\t"\
			"/* eax <-> edx, ebx <-> ecx */\n\t"\
			"movaps	     (%%edx),%%xmm2	\n\t"\
			"movaps	     (%%ebx),%%xmm6	\n\t"\
			"movaps	0x010(%%edx),%%xmm3	\n\t"\
			"movaps	0x010(%%ebx),%%xmm7	\n\t"\
			"movaps	     (%%ecx),%%xmm0	\n\t"\
			"movaps	     (%%eax),%%xmm4	\n\t"\
			"movaps	0x010(%%ecx),%%xmm1	\n\t"\
			"movaps	0x010(%%eax),%%xmm5	\n\t"\
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
			"movaps	%%xmm0,0x240(%%esi)	\n\t"\
			"movaps	%%xmm2,0x360(%%esi)	\n\t"\
			"movaps	%%xmm1,0x250(%%esi)	\n\t"\
			"movaps	%%xmm3,0x130(%%esi)	\n\t"\
			"addpd	%%xmm4,%%xmm4		\n\t"\
			"addpd	%%xmm7,%%xmm7		\n\t"\
			"addpd	%%xmm5,%%xmm5		\n\t"\
			"addpd	%%xmm6,%%xmm6		\n\t"\
			"addpd	%%xmm0,%%xmm4		\n\t"\
			"addpd	%%xmm2,%%xmm7		\n\t"\
			"addpd	%%xmm1,%%xmm5		\n\t"\
			"addpd	%%xmm3,%%xmm6		\n\t"\
			"movaps	%%xmm4,     (%%esi)	\n\t"\
			"movaps	%%xmm7,0x120(%%esi)	\n\t"\
			"movaps	%%xmm5,0x010(%%esi)	\n\t"\
			"movaps	%%xmm6,0x370(%%esi)	\n\t"\
		"/*	add0,1,3,2 = &a[j1+p32]+p0,1,2,3: SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_C(eax,ebx,edx,ecx, 0x120, 0x240, r0a) */\n\t"\
			"addl	$0x040,%%esi	/* r04 */\n\t"\
			"addl	%%edi,%%eax	/* add32 */\n\t"\
			"addl	%%edi,%%ebx\n\t"\
			"addl	%%edi,%%ecx\n\t"\
			"addl	%%edi,%%edx\n\t"\
			"/* ecx <-> edx */\n\t"\
			"movaps	     (%%eax),%%xmm2	\n\t"\
			"movaps	     (%%edx),%%xmm6	\n\t"\
			"movaps	0x010(%%eax),%%xmm3	\n\t"\
			"movaps	0x010(%%edx),%%xmm7	\n\t"\
			"movaps	     (%%ebx),%%xmm0	\n\t"\
			"movaps	     (%%ecx),%%xmm4	\n\t"\
			"movaps	0x010(%%ebx),%%xmm1	\n\t"\
			"movaps	0x010(%%ecx),%%xmm5	\n\t"\
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
			"movaps	%%xmm0,0x240(%%esi)	\n\t"\
			"movaps	%%xmm2,0x360(%%esi)	\n\t"\
			"movaps	%%xmm1,0x250(%%esi)	\n\t"\
			"movaps	%%xmm3,0x130(%%esi)	\n\t"\
			"addpd	%%xmm4,%%xmm4		\n\t"\
			"addpd	%%xmm7,%%xmm7		\n\t"\
			"addpd	%%xmm5,%%xmm5		\n\t"\
			"addpd	%%xmm6,%%xmm6		\n\t"\
			"addpd	%%xmm0,%%xmm4		\n\t"\
			"addpd	%%xmm2,%%xmm7		\n\t"\
			"addpd	%%xmm1,%%xmm5		\n\t"\
			"addpd	%%xmm3,%%xmm6		\n\t"\
			"movaps	%%xmm4,     (%%esi)	\n\t"\
			"movaps	%%xmm7,0x120(%%esi)	\n\t"\
			"movaps	%%xmm5,0x010(%%esi)	\n\t"\
			"movaps	%%xmm6,0x370(%%esi)	\n\t"\
			"\n\t"\
			"/*...and now do 4 radix-9 transforms...*/\n\t"\
			"\n\t"\
			"/* SSE2_RADIX_09_DIT_0TWIDDLE(r00,s1p00r,s1p32r,s1p28r,s1p24r,s1p20r,s1p16r,s1p12r,s1p08r,s1p04r): */\n\t"\
			"movl	%[__out]	,%%esi 		/* __o0-8 = s1p[00,32,28,24,20,16,12,08,04]r; esi,edi store output addresses throughout */\n\t"\
			"movl	%%esi		,%%edi\n\t"\
			"movl	%[__r00]	,%%eax\n\t"\
		"movl	%[__cc1]	,%%edx 		/* edx stores trig addresses throughout */\n\t"\
			"movl	%%eax		,%%ebx\n\t"\
			"movl	%%eax		,%%ecx\n\t"\
			"addl	$0x20		,%%ebx\n\t"\
			"addl	$0x40		,%%ecx\n\t"\
		"addl	$0x40		,%%edx 		/* c3m1 */\n\t"\
			"movaps	    (%%ebx)	,%%xmm2\n\t"\
			"movaps	0x10(%%ebx)	,%%xmm3\n\t"\
			"movaps	    (%%eax)	,%%xmm0\n\t"\
			"movaps	0x10(%%eax)	,%%xmm1\n\t"\
			"movaps	    (%%ecx)	,%%xmm6\n\t"\
			"movaps	0x10(%%ecx)	,%%xmm7\n\t"\
			"movaps	%%xmm2		,%%xmm4\n\t"\
			"movaps	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm6		,%%xmm2\n\t"\
			"addpd	%%xmm7		,%%xmm3\n\t"\
			"subpd	%%xmm6		,%%xmm4\n\t"\
			"subpd	%%xmm7		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm0\n\t"\
			"addpd	%%xmm3		,%%xmm1\n\t"\
			"movaps	    (%%edx)	,%%xmm6\n\t"\
			"movaps	0x10(%%edx)	,%%xmm7\n\t"\
			"movaps	%%xmm0		,    (%%eax)\n\t"\
			"movaps	%%xmm1		,0x10(%%eax)\n\t"\
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
			"movaps	%%xmm2		,    (%%ebx)\n\t"\
			"movaps	%%xmm3		,0x10(%%ebx)\n\t"\
			"movaps	%%xmm0		,    (%%ecx)\n\t"\
			"movaps	%%xmm1		,0x10(%%ecx)\n\t"\
			"addl	$0x60		,%%eax\n\t"\
			"addl	$0x60		,%%ebx\n\t"\
			"addl	$0x60		,%%ecx\n\t"\
			"movaps	    (%%ebx)	,%%xmm4\n\t"\
			"movaps	0x10(%%ebx)	,%%xmm5\n\t"\
			"movaps	    (%%ecx)	,%%xmm2\n\t"\
			"movaps	0x10(%%ecx)	,%%xmm3\n\t"\
			"movaps	    (%%eax)	,%%xmm0\n\t"\
			"movaps	0x10(%%eax)	,%%xmm1\n\t"\
			"subpd	%%xmm2		,%%xmm4\n\t"\
			"subpd	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm2\n\t"\
			"addpd	%%xmm3		,%%xmm3\n\t"\
			"addpd	%%xmm4		,%%xmm2\n\t"\
			"addpd	%%xmm5		,%%xmm3\n\t"\
			"addpd	%%xmm2		,%%xmm0\n\t"\
			"addpd	%%xmm3		,%%xmm1\n\t"\
			"movaps	%%xmm0		,    (%%eax)\n\t"\
			"movaps	%%xmm1		,0x10(%%eax)\n\t"\
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
			"movaps	%%xmm2		,    (%%ebx)\n\t"\
			"movaps	%%xmm3		,0x10(%%ebx)\n\t"\
			"movaps	%%xmm0		,    (%%ecx)\n\t"\
			"movaps	%%xmm1		,0x10(%%ecx)\n\t"\
			"addl	$0x60		,%%eax\n\t"\
			"addl	$0x60		,%%ebx\n\t"\
			"addl	$0x60		,%%ecx\n\t"\
			"movaps	    (%%ebx)	,%%xmm4\n\t"\
			"movaps	0x10(%%ebx)	,%%xmm5\n\t"\
			"movaps	    (%%ecx)	,%%xmm2\n\t"\
			"movaps	0x10(%%ecx)	,%%xmm3\n\t"\
			"movaps	    (%%eax)	,%%xmm0\n\t"\
			"movaps	0x10(%%eax)	,%%xmm1\n\t"\
			"subpd	%%xmm2		,%%xmm4\n\t"\
			"subpd	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm2\n\t"\
			"addpd	%%xmm3		,%%xmm3\n\t"\
			"addpd	%%xmm4		,%%xmm2\n\t"\
			"addpd	%%xmm5		,%%xmm3\n\t"\
			"addpd	%%xmm2		,%%xmm0\n\t"\
			"addpd	%%xmm3		,%%xmm1\n\t"\
			"movaps	%%xmm0		,    (%%eax)\n\t"\
			"movaps	%%xmm1		,0x10(%%eax)\n\t"\
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
			"movaps	%%xmm2		,    (%%ebx)\n\t"\
			"movaps	%%xmm3		,0x10(%%ebx)\n\t"\
			"movaps	%%xmm0		,    (%%ecx)\n\t"\
			"movaps	%%xmm1		,0x10(%%ecx)\n\t"\
			"/******************************************************************************/\n\t"\
			"/*...and now do three more radix-3 transforms, including the twiddle factors: */\n\t"\
			"/******************************************************************************/\n\t"\
			"subl	$0xc0		,%%eax 		/* __r00 */\n\t"\
			"subl	$0x80		,%%ebx 		/* ebx stores r02+0xc0 = r00+0xe0; need r06 = r00+0x60 = r02+0x40 = ebx-0x80 */\n\t"\
			"subl	$0x40		,%%ecx 		/* ecx stores r04+0xc0           ; need r0c = r00+0xc0 = r04+0x80 = ebx-0x40 */\n\t"\
			"movaps	    (%%ebx)	,%%xmm4\n\t"\
			"movaps	0x10(%%ebx)	,%%xmm5\n\t"\
			"movaps	    (%%ecx)	,%%xmm2\n\t"\
			"movaps	0x10(%%ecx)	,%%xmm3\n\t"\
			"movaps	    (%%eax)	,%%xmm0\n\t"\
			"movaps	0x10(%%eax)	,%%xmm1\n\t"\
			"subpd	%%xmm2		,%%xmm4\n\t"\
			"subpd	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm2\n\t"\
			"addpd	%%xmm3		,%%xmm3\n\t"\
			"addpd	%%xmm4		,%%xmm2\n\t"\
			"addpd	%%xmm5		,%%xmm3\n\t"\
			"addpd	%%xmm2		,%%xmm0\n\t"\
			"addpd	%%xmm3		,%%xmm1\n\t"\
			"movaps	%%xmm0		,    (%%esi)\n\t"\
			"movaps	%%xmm1		,0x10(%%esi)\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t"\
			"mulpd	%%xmm7		,%%xmm4\n\t"\
			"mulpd	%%xmm7		,%%xmm5\n\t"\
			"addpd	%%xmm0		,%%xmm2\n\t"\
			"addpd	%%xmm1		,%%xmm3\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t"\
			"addl	$0x300		,%%edi 		/* __o3 = s1p24r */\n\t"\
			"addl	$0x180		,%%esi 		/* __o6 = s1p12r */\n\t"\
			"addpd	%%xmm5		,%%xmm2\n\t"\
			"subpd	%%xmm4		,%%xmm3\n\t"\
			"subpd	%%xmm5		,%%xmm0\n\t"\
			"addpd	%%xmm4		,%%xmm1\n\t"\
			"movaps	%%xmm2		,    (%%edi)\n\t"\
			"movaps	%%xmm3		,0x10(%%edi)\n\t"\
			"movaps	%%xmm0		,    (%%esi)\n\t"\
			"movaps	%%xmm1		,0x10(%%esi)\n\t"\
			"addl	$0x20		,%%eax\n\t"\
			"addl	$0x20		,%%ebx\n\t"\
			"addl	$0x20		,%%ecx\n\t"\
		"subl	$0x40		,%%edx 		/* c1 */\n\t"\
			"movaps	    (%%ebx)	,%%xmm2\n\t"\
			"movaps	0x10(%%ebx)	,%%xmm3\n\t"\
			"movaps	    (%%edx)	,%%xmm6\n\t"\
			"movaps	0x10(%%edx)	,%%xmm7\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"addpd	%%xmm1		,%%xmm2\n\t"\
			"subpd	%%xmm0		,%%xmm3\n\t"\
		"addl	$0x20		,%%edx 		/* c2 */\n\t"\
			"movaps	    (%%ecx)	,%%xmm4\n\t"\
			"movaps	0x10(%%ecx)	,%%xmm5\n\t"\
			"movaps	    (%%edx)	,%%xmm6\n\t"\
			"movaps	0x10(%%edx)	,%%xmm7\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"addpd	%%xmm1		,%%xmm4\n\t"\
			"subpd	%%xmm0		,%%xmm5\n\t"\
		"addl	$0x20		,%%edx 		/* c3m1 */\n\t"\
			"movaps	    (%%edx)	,%%xmm6\n\t"\
			"movaps	0x10(%%edx)	,%%xmm7\n\t"\
			"movaps	    (%%eax)	,%%xmm0\n\t"\
			"movaps	0x10(%%eax)	,%%xmm1\n\t"\
			"subl	$0x080		,%%esi 		/* __o7 = s1p08r */\n\t"\
			"subpd	%%xmm4		,%%xmm2\n\t"\
			"subpd	%%xmm5		,%%xmm3\n\t"\
			"addpd	%%xmm4		,%%xmm4\n\t"\
			"addpd	%%xmm5		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm4\n\t"\
			"addpd	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm4		,%%xmm0\n\t"\
			"addpd	%%xmm5		,%%xmm1\n\t"\
			"movaps	%%xmm0		,    (%%esi)\n\t"\
			"movaps	%%xmm1		,0x10(%%esi)\n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t"\
			"mulpd	%%xmm7		,%%xmm2\n\t"\
			"mulpd	%%xmm7		,%%xmm3\n\t"\
			"addpd	%%xmm0		,%%xmm4\n\t"\
			"addpd	%%xmm1		,%%xmm5\n\t"\
			"addl	$0x100		,%%edi 		/* __o1 = s1p32r */\n\t"\
			"addl	$0x180		,%%esi 		/* __o4 = s1p20r */\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t"\
			"addpd	%%xmm3		,%%xmm4\n\t"\
			"subpd	%%xmm2		,%%xmm5\n\t"\
			"subpd	%%xmm3		,%%xmm0\n\t"\
			"addpd	%%xmm2		,%%xmm1\n\t"\
			"movaps	%%xmm4		,    (%%edi)\n\t"\
			"movaps	%%xmm5		,0x10(%%edi)\n\t"\
			"movaps	%%xmm0		,    (%%esi)\n\t"\
			"movaps	%%xmm1		,0x10(%%esi)\n\t"\
			"addl	$0x20		,%%eax\n\t"\
			"addl	$0x20		,%%ebx\n\t"\
			"addl	$0x20		,%%ecx\n\t"\
		"subl	$0x20		,%%edx 		/* c2 */\n\t"\
			"movaps	    (%%ebx)	,%%xmm2\n\t"\
			"movaps	0x10(%%ebx)	,%%xmm3\n\t"\
			"movaps	    (%%edx)	,%%xmm6\n\t"\
			"movaps	0x10(%%edx)	,%%xmm7\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"addpd	%%xmm1		,%%xmm2\n\t"\
			"subpd	%%xmm0		,%%xmm3\n\t"\
		"addl	$0x40		,%%edx 		/* c4 */\n\t"\
			"movaps	    (%%ecx)	,%%xmm4\n\t"\
			"movaps	0x10(%%ecx)	,%%xmm5\n\t"\
			"movaps	    (%%edx)	,%%xmm6\n\t"\
			"movaps	0x10(%%edx)	,%%xmm7\n\t"\
		"subl	$0x20		,%%edx 		/* c3m1 */\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"addpd	%%xmm1		,%%xmm4\n\t"\
			"subpd	%%xmm0		,%%xmm5\n\t"\
			"movaps	    (%%edx)	,%%xmm6\n\t"\
			"movaps	0x10(%%edx)	,%%xmm7\n\t"\
			"movaps	    (%%eax)	,%%xmm0\n\t"\
			"movaps	0x10(%%eax)	,%%xmm1\n\t"\
			"subl	$0x080		,%%esi 		/* __o5 = s1p16r */\n\t"\
			"subpd	%%xmm4		,%%xmm2\n\t"\
			"subpd	%%xmm5		,%%xmm3\n\t"\
			"addpd	%%xmm4		,%%xmm4\n\t"\
			"addpd	%%xmm5		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm4\n\t"\
			"addpd	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm4		,%%xmm0\n\t"\
			"addpd	%%xmm5		,%%xmm1\n\t"\
			"movaps	%%xmm0		,    (%%esi)\n\t"\
			"movaps	%%xmm1		,0x10(%%esi)\n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t"\
			"mulpd	%%xmm7		,%%xmm2\n\t"\
			"mulpd	%%xmm7		,%%xmm3\n\t"\
			"addpd	%%xmm0		,%%xmm4\n\t"\
			"addpd	%%xmm1		,%%xmm5\n\t"\
			"subl	$0x380		,%%edi 		/* __o8 = s1p04r */\n\t"\
			"addl	$0x180		,%%esi 		/* __o2 = s1p28r */\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t"\
			"addpd	%%xmm3		,%%xmm4\n\t"\
			"subpd	%%xmm2		,%%xmm5\n\t"\
			"subpd	%%xmm3		,%%xmm0\n\t"\
			"addpd	%%xmm2		,%%xmm1\n\t"\
			"movaps	%%xmm4		,    (%%edi)\n\t"\
			"movaps	%%xmm5		,0x10(%%edi)\n\t"\
			"movaps	%%xmm0		,    (%%esi)\n\t"\
			"movaps	%%xmm1		,0x10(%%esi)\n\t"\
			"/* SSE2_RADIX_09_DIT_0TWIDDLE(r10,s1p27r,s1p23r,s1p19r,s1p15r,s1p11r,s1p07r,s1p03r,s1p35r,s1p31r): */\n\t"\
			"movl	%[__out]	,%%esi 		/* __o0-8 = s1p[27,23,19,15,11,07,03,35,31]r; esi,edi store output addresses throughout */\n\t"\
			"addl	$0x360		,%%esi 		/* s1p27 */\n\t"\
			"movl	%%esi		,%%edi\n\t"\
			"movl	%[__r00]	,%%eax\n\t"\
			"addl	$0x120		,%%eax 		/* __r10 */\n\t"\
		"movl	%[__cc1]	,%%edx 		/* edx stores trig addresses throughout */\n\t"\
			"movl	%%eax		,%%ebx\n\t"\
			"movl	%%eax		,%%ecx\n\t"\
			"addl	$0x20		,%%ebx\n\t"\
			"addl	$0x40		,%%ecx\n\t"\
		"addl	$0x40		,%%edx 		/* c3m1 */\n\t"\
			"movaps	    (%%ebx)	,%%xmm2\n\t"\
			"movaps	0x10(%%ebx)	,%%xmm3\n\t"\
			"movaps	    (%%eax)	,%%xmm0\n\t"\
			"movaps	0x10(%%eax)	,%%xmm1\n\t"\
			"movaps	    (%%ecx)	,%%xmm6\n\t"\
			"movaps	0x10(%%ecx)	,%%xmm7\n\t"\
			"movaps	%%xmm2		,%%xmm4\n\t"\
			"movaps	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm6		,%%xmm2\n\t"\
			"addpd	%%xmm7		,%%xmm3\n\t"\
			"subpd	%%xmm6		,%%xmm4\n\t"\
			"subpd	%%xmm7		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm0\n\t"\
			"addpd	%%xmm3		,%%xmm1\n\t"\
			"movaps	    (%%edx)	,%%xmm6\n\t"\
			"movaps	0x10(%%edx)	,%%xmm7\n\t"\
			"movaps	%%xmm0		,    (%%eax)\n\t"\
			"movaps	%%xmm1		,0x10(%%eax)\n\t"\
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
			"movaps	%%xmm2		,    (%%ebx)\n\t"\
			"movaps	%%xmm3		,0x10(%%ebx)\n\t"\
			"movaps	%%xmm0		,    (%%ecx)\n\t"\
			"movaps	%%xmm1		,0x10(%%ecx)\n\t"\
			"addl	$0x60		,%%eax\n\t"\
			"addl	$0x60		,%%ebx\n\t"\
			"addl	$0x60		,%%ecx\n\t"\
			"movaps	    (%%ebx)	,%%xmm4\n\t"\
			"movaps	0x10(%%ebx)	,%%xmm5\n\t"\
			"movaps	    (%%ecx)	,%%xmm2\n\t"\
			"movaps	0x10(%%ecx)	,%%xmm3\n\t"\
			"movaps	    (%%eax)	,%%xmm0\n\t"\
			"movaps	0x10(%%eax)	,%%xmm1\n\t"\
			"subpd	%%xmm2		,%%xmm4\n\t"\
			"subpd	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm2\n\t"\
			"addpd	%%xmm3		,%%xmm3\n\t"\
			"addpd	%%xmm4		,%%xmm2\n\t"\
			"addpd	%%xmm5		,%%xmm3\n\t"\
			"addpd	%%xmm2		,%%xmm0\n\t"\
			"addpd	%%xmm3		,%%xmm1\n\t"\
			"movaps	%%xmm0		,    (%%eax)\n\t"\
			"movaps	%%xmm1		,0x10(%%eax)\n\t"\
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
			"movaps	%%xmm2		,    (%%ebx)\n\t"\
			"movaps	%%xmm3		,0x10(%%ebx)\n\t"\
			"movaps	%%xmm0		,    (%%ecx)\n\t"\
			"movaps	%%xmm1		,0x10(%%ecx)\n\t"\
			"addl	$0x60		,%%eax\n\t"\
			"addl	$0x60		,%%ebx\n\t"\
			"addl	$0x60		,%%ecx\n\t"\
			"movaps	    (%%ebx)	,%%xmm4\n\t"\
			"movaps	0x10(%%ebx)	,%%xmm5\n\t"\
			"movaps	    (%%ecx)	,%%xmm2\n\t"\
			"movaps	0x10(%%ecx)	,%%xmm3\n\t"\
			"movaps	    (%%eax)	,%%xmm0\n\t"\
			"movaps	0x10(%%eax)	,%%xmm1\n\t"\
			"subpd	%%xmm2		,%%xmm4\n\t"\
			"subpd	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm2\n\t"\
			"addpd	%%xmm3		,%%xmm3\n\t"\
			"addpd	%%xmm4		,%%xmm2\n\t"\
			"addpd	%%xmm5		,%%xmm3\n\t"\
			"addpd	%%xmm2		,%%xmm0\n\t"\
			"addpd	%%xmm3		,%%xmm1\n\t"\
			"movaps	%%xmm0		,    (%%eax)\n\t"\
			"movaps	%%xmm1		,0x10(%%eax)\n\t"\
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
			"movaps	%%xmm2		,    (%%ebx)\n\t"\
			"movaps	%%xmm3		,0x10(%%ebx)\n\t"\
			"movaps	%%xmm0		,    (%%ecx)\n\t"\
			"movaps	%%xmm1		,0x10(%%ecx)\n\t"\
			"/******************************************************************************/\n\t"\
			"/*...and now do three more radix-3 transforms, including the twiddle factors: */\n\t"\
			"/******************************************************************************/\n\t"\
			"subl	$0xc0		,%%eax\n\t"\
			"subl	$0x80		,%%ebx\n\t"\
			"subl	$0x40		,%%ecx\n\t"\
			"movaps	    (%%ebx)	,%%xmm4\n\t"\
			"movaps	0x10(%%ebx)	,%%xmm5\n\t"\
			"movaps	    (%%ecx)	,%%xmm2\n\t"\
			"movaps	0x10(%%ecx)	,%%xmm3\n\t"\
			"movaps	    (%%eax)	,%%xmm0\n\t"\
			"movaps	0x10(%%eax)	,%%xmm1\n\t"\
			"subpd	%%xmm2		,%%xmm4\n\t"\
			"subpd	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm2\n\t"\
			"addpd	%%xmm3		,%%xmm3\n\t"\
			"addpd	%%xmm4		,%%xmm2\n\t"\
			"addpd	%%xmm5		,%%xmm3\n\t"\
			"addpd	%%xmm2		,%%xmm0\n\t"\
			"addpd	%%xmm3		,%%xmm1\n\t"\
			"movaps	%%xmm0		,    (%%esi)\n\t"\
			"movaps	%%xmm1		,0x10(%%esi)\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t"\
			"mulpd	%%xmm7		,%%xmm4\n\t"\
			"mulpd	%%xmm7		,%%xmm5\n\t"\
			"addpd	%%xmm0		,%%xmm2\n\t"\
			"addpd	%%xmm1		,%%xmm3\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t"\
			"subl	$0x180		,%%edi 		/* __o3 = s1p15r */\n\t"\
			"subl	$0x300		,%%esi 		/* __o6 = s1p03r */\n\t"\
			"addpd	%%xmm5		,%%xmm2\n\t"\
			"subpd	%%xmm4		,%%xmm3\n\t"\
			"subpd	%%xmm5		,%%xmm0\n\t"\
			"addpd	%%xmm4		,%%xmm1\n\t"\
			"movaps	%%xmm2		,    (%%edi)\n\t"\
			"movaps	%%xmm3		,0x10(%%edi)\n\t"\
			"movaps	%%xmm0		,    (%%esi)\n\t"\
			"movaps	%%xmm1		,0x10(%%esi)\n\t"\
			"addl	$0x20		,%%eax\n\t"\
			"addl	$0x20		,%%ebx\n\t"\
			"addl	$0x20		,%%ecx\n\t"\
		"subl	$0x40		,%%edx 		/* c1 */\n\t"\
			"movaps	    (%%ebx)	,%%xmm2\n\t"\
			"movaps	0x10(%%ebx)	,%%xmm3\n\t"\
			"movaps	    (%%edx)	,%%xmm6\n\t"\
			"movaps	0x10(%%edx)	,%%xmm7\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"addpd	%%xmm1		,%%xmm2\n\t"\
			"subpd	%%xmm0		,%%xmm3\n\t"\
		"addl	$0x20		,%%edx 		/* c2 */\n\t"\
			"movaps	    (%%ecx)	,%%xmm4\n\t"\
			"movaps	0x10(%%ecx)	,%%xmm5\n\t"\
			"movaps	    (%%edx)	,%%xmm6\n\t"\
			"movaps	0x10(%%edx)	,%%xmm7\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"addpd	%%xmm1		,%%xmm4\n\t"\
			"subpd	%%xmm0		,%%xmm5\n\t"\
		"addl	$0x20		,%%edx 		/* c3m1 */\n\t"\
			"movaps	    (%%edx)	,%%xmm6\n\t"\
			"movaps	0x10(%%edx)	,%%xmm7\n\t"\
			"movaps	    (%%eax)	,%%xmm0\n\t"\
			"movaps	0x10(%%eax)	,%%xmm1\n\t"\
			"addl	$0x400		,%%esi 		/* __o7 = s1p35r */\n\t"\
			"subpd	%%xmm4		,%%xmm2\n\t"\
			"subpd	%%xmm5		,%%xmm3\n\t"\
			"addpd	%%xmm4		,%%xmm4\n\t"\
			"addpd	%%xmm5		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm4\n\t"\
			"addpd	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm4		,%%xmm0\n\t"\
			"addpd	%%xmm5		,%%xmm1\n\t"\
			"movaps	%%xmm0		,    (%%esi)\n\t"\
			"movaps	%%xmm1		,0x10(%%esi)\n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t"\
			"mulpd	%%xmm7		,%%xmm2\n\t"\
			"mulpd	%%xmm7		,%%xmm3\n\t"\
			"addpd	%%xmm0		,%%xmm4\n\t"\
			"addpd	%%xmm1		,%%xmm5\n\t"\
			"addl	$0x100		,%%edi 		/* __o1 = s1p23r */\n\t"\
			"subl	$0x300		,%%esi 		/* __o4 = s1p11r */\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t"\
			"addpd	%%xmm3		,%%xmm4\n\t"\
			"subpd	%%xmm2		,%%xmm5\n\t"\
			"subpd	%%xmm3		,%%xmm0\n\t"\
			"addpd	%%xmm2		,%%xmm1\n\t"\
			"movaps	%%xmm4		,    (%%edi)\n\t"\
			"movaps	%%xmm5		,0x10(%%edi)\n\t"\
			"movaps	%%xmm0		,    (%%esi)\n\t"\
			"movaps	%%xmm1		,0x10(%%esi)\n\t"\
			"addl	$0x20		,%%eax\n\t"\
			"addl	$0x20		,%%ebx\n\t"\
			"addl	$0x20		,%%ecx\n\t"\
		"subl	$0x20		,%%edx 		/* c2 */\n\t"\
			"movaps	    (%%ebx)	,%%xmm2\n\t"\
			"movaps	0x10(%%ebx)	,%%xmm3\n\t"\
			"movaps	    (%%edx)	,%%xmm6\n\t"\
			"movaps	0x10(%%edx)	,%%xmm7\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"addpd	%%xmm1		,%%xmm2\n\t"\
			"subpd	%%xmm0		,%%xmm3\n\t"\
		"addl	$0x40		,%%edx 		/* c4 */\n\t"\
			"movaps	    (%%ecx)	,%%xmm4\n\t"\
			"movaps	0x10(%%ecx)	,%%xmm5\n\t"\
			"movaps	    (%%edx)	,%%xmm6\n\t"\
			"movaps	0x10(%%edx)	,%%xmm7\n\t"\
		"subl	$0x20		,%%edx 		/* c3m1 */\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"addpd	%%xmm1		,%%xmm4\n\t"\
			"subpd	%%xmm0		,%%xmm5\n\t"\
			"movaps	    (%%edx)	,%%xmm6\n\t"\
			"movaps	0x10(%%edx)	,%%xmm7\n\t"\
			"movaps	    (%%eax)	,%%xmm0\n\t"\
			"movaps	0x10(%%eax)	,%%xmm1\n\t"\
			"subl	$0x080		,%%esi 		/* __o5 = s1p07r */\n\t"\
			"subpd	%%xmm4		,%%xmm2\n\t"\
			"subpd	%%xmm5		,%%xmm3\n\t"\
			"addpd	%%xmm4		,%%xmm4\n\t"\
			"addpd	%%xmm5		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm4\n\t"\
			"addpd	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm4		,%%xmm0\n\t"\
			"addpd	%%xmm5		,%%xmm1\n\t"\
			"movaps	%%xmm0		,    (%%esi)\n\t"\
			"movaps	%%xmm1		,0x10(%%esi)\n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t"\
			"mulpd	%%xmm7		,%%xmm2\n\t"\
			"mulpd	%%xmm7		,%%xmm3\n\t"\
			"addpd	%%xmm0		,%%xmm4\n\t"\
			"addpd	%%xmm1		,%%xmm5\n\t"\
			"addl	$0x100		,%%edi 		/* __o8 = s1p31r */\n\t"\
			"addl	$0x180		,%%esi 		/* __o2 = s1p19r */\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t"\
			"addpd	%%xmm3		,%%xmm4\n\t"\
			"subpd	%%xmm2		,%%xmm5\n\t"\
			"subpd	%%xmm3		,%%xmm0\n\t"\
			"addpd	%%xmm2		,%%xmm1\n\t"\
			"movaps	%%xmm4		,    (%%edi)\n\t"\
			"movaps	%%xmm5		,0x10(%%edi)\n\t"\
			"movaps	%%xmm0		,    (%%esi)\n\t"\
			"movaps	%%xmm1		,0x10(%%esi)\n\t"\
			"/* SSE2_RADIX_09_DIT_0TWIDDLE(r20,s1p18r,s1p14r,s1p10r,s1p06r,s1p02r,s1p34r,s1p30r,s1p26r,s1p22r): */\n\t"\
			"movl	%[__out]	,%%esi 		/* __o0-8 = s1p[18,14,10,06,02,34,30,26,22]r; esi,edi store output addresses throughout */\n\t"\
			"addl	$0x240		,%%esi 		/* s1p18 */\n\t"\
			"movl	%%esi		,%%edi\n\t"\
			"movl	%[__r00]	,%%eax\n\t"\
			"addl	$0x240		,%%eax 		/* __r20 */\n\t"\
		"movl	%[__cc1]	,%%edx 		/* edx stores trig addresses throughout */\n\t"\
			"movl	%%eax		,%%ebx\n\t"\
			"movl	%%eax		,%%ecx\n\t"\
			"addl	$0x20		,%%ebx\n\t"\
			"addl	$0x40		,%%ecx\n\t"\
		"addl	$0x40		,%%edx 		/* c3m1 */\n\t"\
			"movaps	    (%%ebx)	,%%xmm2\n\t"\
			"movaps	0x10(%%ebx)	,%%xmm3\n\t"\
			"movaps	    (%%eax)	,%%xmm0\n\t"\
			"movaps	0x10(%%eax)	,%%xmm1\n\t"\
			"movaps	    (%%ecx)	,%%xmm6\n\t"\
			"movaps	0x10(%%ecx)	,%%xmm7\n\t"\
			"movaps	%%xmm2		,%%xmm4\n\t"\
			"movaps	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm6		,%%xmm2\n\t"\
			"addpd	%%xmm7		,%%xmm3\n\t"\
			"subpd	%%xmm6		,%%xmm4\n\t"\
			"subpd	%%xmm7		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm0\n\t"\
			"addpd	%%xmm3		,%%xmm1\n\t"\
			"movaps	    (%%edx)	,%%xmm6\n\t"\
			"movaps	0x10(%%edx)	,%%xmm7\n\t"\
			"movaps	%%xmm0		,    (%%eax)\n\t"\
			"movaps	%%xmm1		,0x10(%%eax)\n\t"\
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
			"movaps	%%xmm2		,    (%%ebx)\n\t"\
			"movaps	%%xmm3		,0x10(%%ebx)\n\t"\
			"movaps	%%xmm0		,    (%%ecx)\n\t"\
			"movaps	%%xmm1		,0x10(%%ecx)\n\t"\
			"addl	$0x60		,%%eax\n\t"\
			"addl	$0x60		,%%ebx\n\t"\
			"addl	$0x60		,%%ecx\n\t"\
			"movaps	    (%%ebx)	,%%xmm4\n\t"\
			"movaps	0x10(%%ebx)	,%%xmm5\n\t"\
			"movaps	    (%%ecx)	,%%xmm2\n\t"\
			"movaps	0x10(%%ecx)	,%%xmm3\n\t"\
			"movaps	    (%%eax)	,%%xmm0\n\t"\
			"movaps	0x10(%%eax)	,%%xmm1\n\t"\
			"subpd	%%xmm2		,%%xmm4\n\t"\
			"subpd	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm2\n\t"\
			"addpd	%%xmm3		,%%xmm3\n\t"\
			"addpd	%%xmm4		,%%xmm2\n\t"\
			"addpd	%%xmm5		,%%xmm3\n\t"\
			"addpd	%%xmm2		,%%xmm0\n\t"\
			"addpd	%%xmm3		,%%xmm1\n\t"\
			"movaps	%%xmm0		,    (%%eax)\n\t"\
			"movaps	%%xmm1		,0x10(%%eax)\n\t"\
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
			"movaps	%%xmm2		,    (%%ebx)\n\t"\
			"movaps	%%xmm3		,0x10(%%ebx)\n\t"\
			"movaps	%%xmm0		,    (%%ecx)\n\t"\
			"movaps	%%xmm1		,0x10(%%ecx)\n\t"\
			"addl	$0x60		,%%eax\n\t"\
			"addl	$0x60		,%%ebx\n\t"\
			"addl	$0x60		,%%ecx\n\t"\
			"movaps	    (%%ebx)	,%%xmm4\n\t"\
			"movaps	0x10(%%ebx)	,%%xmm5\n\t"\
			"movaps	    (%%ecx)	,%%xmm2\n\t"\
			"movaps	0x10(%%ecx)	,%%xmm3\n\t"\
			"movaps	    (%%eax)	,%%xmm0\n\t"\
			"movaps	0x10(%%eax)	,%%xmm1\n\t"\
			"subpd	%%xmm2		,%%xmm4\n\t"\
			"subpd	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm2\n\t"\
			"addpd	%%xmm3		,%%xmm3\n\t"\
			"addpd	%%xmm4		,%%xmm2\n\t"\
			"addpd	%%xmm5		,%%xmm3\n\t"\
			"addpd	%%xmm2		,%%xmm0\n\t"\
			"addpd	%%xmm3		,%%xmm1\n\t"\
			"movaps	%%xmm0		,    (%%eax)\n\t"\
			"movaps	%%xmm1		,0x10(%%eax)\n\t"\
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
			"movaps	%%xmm2		,    (%%ebx)\n\t"\
			"movaps	%%xmm3		,0x10(%%ebx)\n\t"\
			"movaps	%%xmm0		,    (%%ecx)\n\t"\
			"movaps	%%xmm1		,0x10(%%ecx)\n\t"\
			"/******************************************************************************/\n\t"\
			"/*...and now do three more radix-3 transforms, including the twiddle factors: */\n\t"\
			"/******************************************************************************/\n\t"\
			"subl	$0xc0		,%%eax\n\t"\
			"subl	$0x80		,%%ebx\n\t"\
			"subl	$0x40		,%%ecx\n\t"\
			"movaps	    (%%ebx)	,%%xmm4\n\t"\
			"movaps	0x10(%%ebx)	,%%xmm5\n\t"\
			"movaps	    (%%ecx)	,%%xmm2\n\t"\
			"movaps	0x10(%%ecx)	,%%xmm3\n\t"\
			"movaps	    (%%eax)	,%%xmm0\n\t"\
			"movaps	0x10(%%eax)	,%%xmm1\n\t"\
			"subpd	%%xmm2		,%%xmm4\n\t"\
			"subpd	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm2\n\t"\
			"addpd	%%xmm3		,%%xmm3\n\t"\
			"addpd	%%xmm4		,%%xmm2\n\t"\
			"addpd	%%xmm5		,%%xmm3\n\t"\
			"addpd	%%xmm2		,%%xmm0\n\t"\
			"addpd	%%xmm3		,%%xmm1\n\t"\
			"movaps	%%xmm0		,    (%%esi)\n\t"\
			"movaps	%%xmm1		,0x10(%%esi)\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t"\
			"mulpd	%%xmm7		,%%xmm4\n\t"\
			"mulpd	%%xmm7		,%%xmm5\n\t"\
			"addpd	%%xmm0		,%%xmm2\n\t"\
			"addpd	%%xmm1		,%%xmm3\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t"\
			"subl	$0x180		,%%edi 		/* __o3 = s1p06r */\n\t"\
			"addl	$0x180		,%%esi 		/* __o6 = s1p30r */\n\t"\
			"addpd	%%xmm5		,%%xmm2\n\t"\
			"subpd	%%xmm4		,%%xmm3\n\t"\
			"subpd	%%xmm5		,%%xmm0\n\t"\
			"addpd	%%xmm4		,%%xmm1\n\t"\
			"movaps	%%xmm2		,    (%%edi)\n\t"\
			"movaps	%%xmm3		,0x10(%%edi)\n\t"\
			"movaps	%%xmm0		,    (%%esi)\n\t"\
			"movaps	%%xmm1		,0x10(%%esi)\n\t"\
			"addl	$0x20		,%%eax\n\t"\
			"addl	$0x20		,%%ebx\n\t"\
			"addl	$0x20		,%%ecx\n\t"\
		"subl	$0x40		,%%edx 		/* c1 */\n\t"\
			"movaps	    (%%ebx)	,%%xmm2\n\t"\
			"movaps	0x10(%%ebx)	,%%xmm3\n\t"\
			"movaps	    (%%edx)	,%%xmm6\n\t"\
			"movaps	0x10(%%edx)	,%%xmm7\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"addpd	%%xmm1		,%%xmm2\n\t"\
			"subpd	%%xmm0		,%%xmm3\n\t"\
		"addl	$0x20		,%%edx 		/* c2 */\n\t"\
			"movaps	    (%%ecx)	,%%xmm4\n\t"\
			"movaps	0x10(%%ecx)	,%%xmm5\n\t"\
			"movaps	    (%%edx)	,%%xmm6\n\t"\
			"movaps	0x10(%%edx)	,%%xmm7\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"addpd	%%xmm1		,%%xmm4\n\t"\
			"subpd	%%xmm0		,%%xmm5\n\t"\
		"addl	$0x20		,%%edx 		/* c3m1 */\n\t"\
			"movaps	    (%%edx)	,%%xmm6\n\t"\
			"movaps	0x10(%%edx)	,%%xmm7\n\t"\
			"movaps	    (%%eax)	,%%xmm0\n\t"\
			"movaps	0x10(%%eax)	,%%xmm1\n\t"\
			"subl	$0x080		,%%esi 		/* __o7 = s1p26r */\n\t"\
			"subpd	%%xmm4		,%%xmm2\n\t"\
			"subpd	%%xmm5		,%%xmm3\n\t"\
			"addpd	%%xmm4		,%%xmm4\n\t"\
			"addpd	%%xmm5		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm4\n\t"\
			"addpd	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm4		,%%xmm0\n\t"\
			"addpd	%%xmm5		,%%xmm1\n\t"\
			"movaps	%%xmm0		,    (%%esi)\n\t"\
			"movaps	%%xmm1		,0x10(%%esi)\n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t"\
			"mulpd	%%xmm7		,%%xmm2\n\t"\
			"mulpd	%%xmm7		,%%xmm3\n\t"\
			"addpd	%%xmm0		,%%xmm4\n\t"\
			"addpd	%%xmm1		,%%xmm5\n\t"\
			"addl	$0x100		,%%edi 		/* __o1 = s1p14r */\n\t"\
			"subl	$0x300		,%%esi 		/* __o4 = s1p02r */\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t"\
			"addpd	%%xmm3		,%%xmm4\n\t"\
			"subpd	%%xmm2		,%%xmm5\n\t"\
			"subpd	%%xmm3		,%%xmm0\n\t"\
			"addpd	%%xmm2		,%%xmm1\n\t"\
			"movaps	%%xmm4		,    (%%edi)\n\t"\
			"movaps	%%xmm5		,0x10(%%edi)\n\t"\
			"movaps	%%xmm0		,    (%%esi)\n\t"\
			"movaps	%%xmm1		,0x10(%%esi)\n\t"\
			"addl	$0x20		,%%eax\n\t"\
			"addl	$0x20		,%%ebx\n\t"\
			"addl	$0x20		,%%ecx\n\t"\
		"subl	$0x20		,%%edx 		/* c2 */\n\t"\
			"movaps	    (%%ebx)	,%%xmm2\n\t"\
			"movaps	0x10(%%ebx)	,%%xmm3\n\t"\
			"movaps	    (%%edx)	,%%xmm6\n\t"\
			"movaps	0x10(%%edx)	,%%xmm7\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"addpd	%%xmm1		,%%xmm2\n\t"\
			"subpd	%%xmm0		,%%xmm3\n\t"\
		"addl	$0x40		,%%edx 		/* c4 */\n\t"\
			"movaps	    (%%ecx)	,%%xmm4\n\t"\
			"movaps	0x10(%%ecx)	,%%xmm5\n\t"\
			"movaps	    (%%edx)	,%%xmm6\n\t"\
			"movaps	0x10(%%edx)	,%%xmm7\n\t"\
		"subl	$0x20		,%%edx 		/* c3m1 */\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"addpd	%%xmm1		,%%xmm4\n\t"\
			"subpd	%%xmm0		,%%xmm5\n\t"\
			"movaps	    (%%edx)	,%%xmm6\n\t"\
			"movaps	0x10(%%edx)	,%%xmm7\n\t"\
			"movaps	    (%%eax)	,%%xmm0\n\t"\
			"movaps	0x10(%%eax)	,%%xmm1\n\t"\
			"addl	$0x400		,%%esi 		/* __o5 = s1p34r */\n\t"\
			"subpd	%%xmm4		,%%xmm2\n\t"\
			"subpd	%%xmm5		,%%xmm3\n\t"\
			"addpd	%%xmm4		,%%xmm4\n\t"\
			"addpd	%%xmm5		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm4\n\t"\
			"addpd	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm4		,%%xmm0\n\t"\
			"addpd	%%xmm5		,%%xmm1\n\t"\
			"movaps	%%xmm0		,    (%%esi)\n\t"\
			"movaps	%%xmm1		,0x10(%%esi)\n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t"\
			"mulpd	%%xmm7		,%%xmm2\n\t"\
			"mulpd	%%xmm7		,%%xmm3\n\t"\
			"addpd	%%xmm0		,%%xmm4\n\t"\
			"addpd	%%xmm1		,%%xmm5\n\t"\
			"addl	$0x100		,%%edi 		/* __o8 = s1p22r */\n\t"\
			"subl	$0x300		,%%esi 		/* __o2 = s1p10r */\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t"\
			"addpd	%%xmm3		,%%xmm4\n\t"\
			"subpd	%%xmm2		,%%xmm5\n\t"\
			"subpd	%%xmm3		,%%xmm0\n\t"\
			"addpd	%%xmm2		,%%xmm1\n\t"\
			"movaps	%%xmm4		,    (%%edi)\n\t"\
			"movaps	%%xmm5		,0x10(%%edi)\n\t"\
			"movaps	%%xmm0		,    (%%esi)\n\t"\
			"movaps	%%xmm1		,0x10(%%esi)\n\t"\
			"/* SSE2_RADIX_09_DIT_0TWIDDLE(r30,s1p09r,s1p05r,s1p01r,s1p33r,s1p29r,s1p25r,s1p21r,s1p17r,s1p13r): */\n\t"\
			"movl	%[__out]	,%%esi 		/* __o0-8 = s1p[09,05,01,33,29,25,21,17,13]r; esi,edi store output addresses throughout */\n\t"\
			"addl	$0x120		,%%esi 		/* s1p09 */\n\t"\
			"movl	%%esi		,%%edi\n\t"\
			"movl	%[__r00]	,%%eax\n\t"\
			"addl	$0x360		,%%eax 		/* __r20 */\n\t"\
		"movl	%[__cc1]	,%%edx 		/* edx stores trig addresses throughout */\n\t"\
			"movl	%%eax		,%%ebx\n\t"\
			"movl	%%eax		,%%ecx\n\t"\
			"addl	$0x20		,%%ebx\n\t"\
			"addl	$0x40		,%%ecx\n\t"\
		"addl	$0x40		,%%edx 		/* c3m1 */\n\t"\
			"movaps	    (%%ebx)	,%%xmm2\n\t"\
			"movaps	0x10(%%ebx)	,%%xmm3\n\t"\
			"movaps	    (%%eax)	,%%xmm0\n\t"\
			"movaps	0x10(%%eax)	,%%xmm1\n\t"\
			"movaps	    (%%ecx)	,%%xmm6\n\t"\
			"movaps	0x10(%%ecx)	,%%xmm7\n\t"\
			"movaps	%%xmm2		,%%xmm4\n\t"\
			"movaps	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm6		,%%xmm2\n\t"\
			"addpd	%%xmm7		,%%xmm3\n\t"\
			"subpd	%%xmm6		,%%xmm4\n\t"\
			"subpd	%%xmm7		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm0\n\t"\
			"addpd	%%xmm3		,%%xmm1\n\t"\
			"movaps	    (%%edx)	,%%xmm6\n\t"\
			"movaps	0x10(%%edx)	,%%xmm7\n\t"\
			"movaps	%%xmm0		,    (%%eax)\n\t"\
			"movaps	%%xmm1		,0x10(%%eax)\n\t"\
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
			"movaps	%%xmm2		,    (%%ebx)\n\t"\
			"movaps	%%xmm3		,0x10(%%ebx)\n\t"\
			"movaps	%%xmm0		,    (%%ecx)\n\t"\
			"movaps	%%xmm1		,0x10(%%ecx)\n\t"\
			"addl	$0x60		,%%eax\n\t"\
			"addl	$0x60		,%%ebx\n\t"\
			"addl	$0x60		,%%ecx\n\t"\
			"movaps	    (%%ebx)	,%%xmm4\n\t"\
			"movaps	0x10(%%ebx)	,%%xmm5\n\t"\
			"movaps	    (%%ecx)	,%%xmm2\n\t"\
			"movaps	0x10(%%ecx)	,%%xmm3\n\t"\
			"movaps	    (%%eax)	,%%xmm0\n\t"\
			"movaps	0x10(%%eax)	,%%xmm1\n\t"\
			"subpd	%%xmm2		,%%xmm4\n\t"\
			"subpd	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm2\n\t"\
			"addpd	%%xmm3		,%%xmm3\n\t"\
			"addpd	%%xmm4		,%%xmm2\n\t"\
			"addpd	%%xmm5		,%%xmm3\n\t"\
			"addpd	%%xmm2		,%%xmm0\n\t"\
			"addpd	%%xmm3		,%%xmm1\n\t"\
			"movaps	%%xmm0		,    (%%eax)\n\t"\
			"movaps	%%xmm1		,0x10(%%eax)\n\t"\
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
			"movaps	%%xmm2		,    (%%ebx)\n\t"\
			"movaps	%%xmm3		,0x10(%%ebx)\n\t"\
			"movaps	%%xmm0		,    (%%ecx)\n\t"\
			"movaps	%%xmm1		,0x10(%%ecx)\n\t"\
			"addl	$0x60		,%%eax\n\t"\
			"addl	$0x60		,%%ebx\n\t"\
			"addl	$0x60		,%%ecx\n\t"\
			"movaps	    (%%ebx)	,%%xmm4\n\t"\
			"movaps	0x10(%%ebx)	,%%xmm5\n\t"\
			"movaps	    (%%ecx)	,%%xmm2\n\t"\
			"movaps	0x10(%%ecx)	,%%xmm3\n\t"\
			"movaps	    (%%eax)	,%%xmm0\n\t"\
			"movaps	0x10(%%eax)	,%%xmm1\n\t"\
			"subpd	%%xmm2		,%%xmm4\n\t"\
			"subpd	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm2\n\t"\
			"addpd	%%xmm3		,%%xmm3\n\t"\
			"addpd	%%xmm4		,%%xmm2\n\t"\
			"addpd	%%xmm5		,%%xmm3\n\t"\
			"addpd	%%xmm2		,%%xmm0\n\t"\
			"addpd	%%xmm3		,%%xmm1\n\t"\
			"movaps	%%xmm0		,    (%%eax)\n\t"\
			"movaps	%%xmm1		,0x10(%%eax)\n\t"\
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
			"movaps	%%xmm2		,    (%%ebx)\n\t"\
			"movaps	%%xmm3		,0x10(%%ebx)\n\t"\
			"movaps	%%xmm0		,    (%%ecx)\n\t"\
			"movaps	%%xmm1		,0x10(%%ecx)\n\t"\
			"/******************************************************************************/\n\t"\
			"/*...and now do three more radix-3 transforms, including the twiddle factors: */\n\t"\
			"/******************************************************************************/\n\t"\
			"subl	$0xc0		,%%eax\n\t"\
			"subl	$0x80		,%%ebx\n\t"\
			"subl	$0x40		,%%ecx\n\t"\
			"movaps	    (%%ebx)	,%%xmm4\n\t"\
			"movaps	0x10(%%ebx)	,%%xmm5\n\t"\
			"movaps	    (%%ecx)	,%%xmm2\n\t"\
			"movaps	0x10(%%ecx)	,%%xmm3\n\t"\
			"movaps	    (%%eax)	,%%xmm0\n\t"\
			"movaps	0x10(%%eax)	,%%xmm1\n\t"\
			"subpd	%%xmm2		,%%xmm4\n\t"\
			"subpd	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm2\n\t"\
			"addpd	%%xmm3		,%%xmm3\n\t"\
			"addpd	%%xmm4		,%%xmm2\n\t"\
			"addpd	%%xmm5		,%%xmm3\n\t"\
			"addpd	%%xmm2		,%%xmm0\n\t"\
			"addpd	%%xmm3		,%%xmm1\n\t"\
			"movaps	%%xmm0		,    (%%esi)\n\t"\
			"movaps	%%xmm1		,0x10(%%esi)\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t"\
			"mulpd	%%xmm7		,%%xmm4\n\t"\
			"mulpd	%%xmm7		,%%xmm5\n\t"\
			"addpd	%%xmm0		,%%xmm2\n\t"\
			"addpd	%%xmm1		,%%xmm3\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t"\
			"addl	$0x300		,%%edi 		/* __o3 = s1p33r */\n\t"\
			"addl	$0x180		,%%esi 		/* __o6 = s1p21r */\n\t"\
			"addpd	%%xmm5		,%%xmm2\n\t"\
			"subpd	%%xmm4		,%%xmm3\n\t"\
			"subpd	%%xmm5		,%%xmm0\n\t"\
			"addpd	%%xmm4		,%%xmm1\n\t"\
			"movaps	%%xmm2		,    (%%edi)\n\t"\
			"movaps	%%xmm3		,0x10(%%edi)\n\t"\
			"movaps	%%xmm0		,    (%%esi)\n\t"\
			"movaps	%%xmm1		,0x10(%%esi)\n\t"\
			"addl	$0x20		,%%eax\n\t"\
			"addl	$0x20		,%%ebx\n\t"\
			"addl	$0x20		,%%ecx\n\t"\
		"subl	$0x40		,%%edx 		/* c1 */\n\t"\
			"movaps	    (%%ebx)	,%%xmm2\n\t"\
			"movaps	0x10(%%ebx)	,%%xmm3\n\t"\
			"movaps	    (%%edx)	,%%xmm6\n\t"\
			"movaps	0x10(%%edx)	,%%xmm7\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"addpd	%%xmm1		,%%xmm2\n\t"\
			"subpd	%%xmm0		,%%xmm3\n\t"\
		"addl	$0x20		,%%edx 		/* c2 */\n\t"\
			"movaps	    (%%ecx)	,%%xmm4\n\t"\
			"movaps	0x10(%%ecx)	,%%xmm5\n\t"\
			"movaps	    (%%edx)	,%%xmm6\n\t"\
			"movaps	0x10(%%edx)	,%%xmm7\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"addpd	%%xmm1		,%%xmm4\n\t"\
			"subpd	%%xmm0		,%%xmm5\n\t"\
		"addl	$0x20		,%%edx 		/* c3m1 */\n\t"\
			"movaps	    (%%edx)	,%%xmm6\n\t"\
			"movaps	0x10(%%edx)	,%%xmm7\n\t"\
			"movaps	    (%%eax)	,%%xmm0\n\t"\
			"movaps	0x10(%%eax)	,%%xmm1\n\t"\
			"subl	$0x080		,%%esi 		/* __o7 = s1p17r */\n\t"\
			"subpd	%%xmm4		,%%xmm2\n\t"\
			"subpd	%%xmm5		,%%xmm3\n\t"\
			"addpd	%%xmm4		,%%xmm4\n\t"\
			"addpd	%%xmm5		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm4\n\t"\
			"addpd	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm4		,%%xmm0\n\t"\
			"addpd	%%xmm5		,%%xmm1\n\t"\
			"movaps	%%xmm0		,    (%%esi)\n\t"\
			"movaps	%%xmm1		,0x10(%%esi)\n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t"\
			"mulpd	%%xmm7		,%%xmm2\n\t"\
			"mulpd	%%xmm7		,%%xmm3\n\t"\
			"addpd	%%xmm0		,%%xmm4\n\t"\
			"addpd	%%xmm1		,%%xmm5\n\t"\
			"subl	$0x380		,%%edi 		/* __o1 = s1p05r */\n\t"\
			"addl	$0x180		,%%esi 		/* __o4 = s1p29r */\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t"\
			"addpd	%%xmm3		,%%xmm4\n\t"\
			"subpd	%%xmm2		,%%xmm5\n\t"\
			"subpd	%%xmm3		,%%xmm0\n\t"\
			"addpd	%%xmm2		,%%xmm1\n\t"\
			"movaps	%%xmm4		,    (%%edi)\n\t"\
			"movaps	%%xmm5		,0x10(%%edi)\n\t"\
			"movaps	%%xmm0		,    (%%esi)\n\t"\
			"movaps	%%xmm1		,0x10(%%esi)\n\t"\
			"addl	$0x20		,%%eax\n\t"\
			"addl	$0x20		,%%ebx\n\t"\
			"addl	$0x20		,%%ecx\n\t"\
		"subl	$0x20		,%%edx 		/* c2 */\n\t"\
			"movaps	    (%%ebx)	,%%xmm2\n\t"\
			"movaps	0x10(%%ebx)	,%%xmm3\n\t"\
			"movaps	    (%%edx)	,%%xmm6\n\t"\
			"movaps	0x10(%%edx)	,%%xmm7\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"addpd	%%xmm1		,%%xmm2\n\t"\
			"subpd	%%xmm0		,%%xmm3\n\t"\
		"addl	$0x40		,%%edx 		/* c4 */\n\t"\
			"movaps	    (%%ecx)	,%%xmm4\n\t"\
			"movaps	0x10(%%ecx)	,%%xmm5\n\t"\
			"movaps	    (%%edx)	,%%xmm6\n\t"\
			"movaps	0x10(%%edx)	,%%xmm7\n\t"\
		"subl	$0x20		,%%edx 		/* c3m1 */\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"addpd	%%xmm1		,%%xmm4\n\t"\
			"subpd	%%xmm0		,%%xmm5\n\t"\
			"movaps	    (%%edx)	,%%xmm6\n\t"\
			"movaps	0x10(%%edx)	,%%xmm7\n\t"\
			"movaps	    (%%eax)	,%%xmm0\n\t"\
			"movaps	0x10(%%eax)	,%%xmm1\n\t"\
			"subl	$0x080		,%%esi 		/* __o5 = s1p25r */\n\t"\
			"subpd	%%xmm4		,%%xmm2\n\t"\
			"subpd	%%xmm5		,%%xmm3\n\t"\
			"addpd	%%xmm4		,%%xmm4\n\t"\
			"addpd	%%xmm5		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm4\n\t"\
			"addpd	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm4		,%%xmm0\n\t"\
			"addpd	%%xmm5		,%%xmm1\n\t"\
			"movaps	%%xmm0		,    (%%esi)\n\t"\
			"movaps	%%xmm1		,0x10(%%esi)\n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t"\
			"mulpd	%%xmm7		,%%xmm2\n\t"\
			"mulpd	%%xmm7		,%%xmm3\n\t"\
			"addpd	%%xmm0		,%%xmm4\n\t"\
			"addpd	%%xmm1		,%%xmm5\n\t"\
			"addl	$0x100		,%%edi 		/* __o8 = s1p13r */\n\t"\
			"subl	$0x300		,%%esi 		/* __o2 = s1p01r */\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t"\
			"addpd	%%xmm3		,%%xmm4\n\t"\
			"subpd	%%xmm2		,%%xmm5\n\t"\
			"subpd	%%xmm3		,%%xmm0\n\t"\
			"addpd	%%xmm2		,%%xmm1\n\t"\
			"movaps	%%xmm4		,    (%%edi)\n\t"\
			"movaps	%%xmm5		,0x10(%%edi)\n\t"\
			"movaps	%%xmm0		,    (%%esi)\n\t"\
			"movaps	%%xmm1		,0x10(%%esi)\n\t"\
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
			: "eax","ebx","ecx","edx","edi","esi"		/* Clobbered registers */\
		);\
	}


	#define	SSE2_RADIX36_DIF_NOTWIDDLE(Xadd,Xp01,Xp02,Xp03,Xp04,Xp08,Xp12,Xp16,Xp20,Xp24,Xp28,Xp32,Xr00,Xin0,Xcc1)\
	{\
	__asm__ volatile (\
			"/* SSE2_RADIX_09_DIF_0TWIDDLE(r00,s1p00r,s1p32r,s1p28r,s1p24r,s1p20r,s1p16r,s1p12r,s1p08r,s1p04r) */\n\t"\
			"movl	%[__in0]	,%%eax 		/* __i0-8 = s1p[00,32,28,24,20,16,12,08,04]r; e[abc]x store output addresses throughout */\n\t"\
			"movl	%[__r00]	,%%esi\n\t"\
		"movl	%[__cc1]	,%%edx 		/* edx stores trig addresses throughout */\n\t"\
			"movl	%%eax		,%%ebx\n\t"\
			"movl	%%eax		,%%ecx\n\t"\
			"addl	$0x300		,%%ebx 		/* __i3 = s1p24r */\n\t"\
			"addl	$0x180		,%%ecx 		/* __i6 = s1p12r */\n\t"\
		"addl	$0x40		,%%edx 		/* c3m1 */\n\t"\
			"movaps	    (%%ebx)	,%%xmm2\n\t"\
			"movaps	0x10(%%ebx)	,%%xmm3\n\t"\
			"movaps	    (%%eax)	,%%xmm0\n\t"\
			"movaps	0x10(%%eax)	,%%xmm1\n\t"\
			"movaps	    (%%ecx)	,%%xmm6\n\t"\
			"movaps	0x10(%%ecx)	,%%xmm7\n\t"\
			"movaps	%%xmm2		,%%xmm4\n\t"\
			"movaps	%%xmm3		,%%xmm5\n\t"\
		"movl	%%esi		,%%edi\n\t"\
		"addl	$0x20		,%%edi\n\t"\
			"addpd	%%xmm6		,%%xmm2\n\t"\
			"addpd	%%xmm7		,%%xmm3\n\t"\
			"subpd	%%xmm6		,%%xmm4\n\t"\
			"subpd	%%xmm7		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm0\n\t"\
			"addpd	%%xmm3		,%%xmm1\n\t"\
			"movaps	    (%%edx)	,%%xmm6\n\t"\
			"movaps	0x10(%%edx)	,%%xmm7\n\t"\
			"movaps	%%xmm0		,    (%%esi)	/* <- t00 */\n\t"\
			"movaps	%%xmm1		,0x10(%%esi)	/* <- t01 */\n\t"\
		"addl	$0x40		,%%esi\n\t"\
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
			"movaps	%%xmm2		,    (%%edi)	/* <- t02 */\n\t"\
			"movaps	%%xmm3		,0x10(%%edi)	/* <- t03 */\n\t"\
			"movaps	%%xmm0		,    (%%esi)	/* <- t04 */\n\t"\
			"movaps	%%xmm1		,0x10(%%esi)	/* <- t05 */\n\t"\
			"addl	$0x400		,%%eax 		/* __i1 = s1p32r */\n\t"\
			"subl	$0x080		,%%ebx 		/* __i4 = s1p20r */\n\t"\
			"subl	$0x080		,%%ecx 		/* __i7 = s1p08r */\n\t"\
			"movaps	    (%%ebx)	,%%xmm4\n\t"\
			"movaps	0x10(%%ebx)	,%%xmm5\n\t"\
			"movaps	    (%%eax)	,%%xmm0\n\t"\
			"movaps	0x10(%%eax)	,%%xmm1\n\t"\
			"movaps	    (%%ecx)	,%%xmm2\n\t"\
			"movaps	0x10(%%ecx)	,%%xmm3\n\t"\
		"addl	$0x40		,%%edi\n\t"\
			"subpd	%%xmm2		,%%xmm4\n\t"\
			"subpd	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm2\n\t"\
			"addpd	%%xmm3		,%%xmm3\n\t"\
			"addpd	%%xmm4		,%%xmm2\n\t"\
			"addpd	%%xmm5		,%%xmm3\n\t"\
			"addpd	%%xmm2		,%%xmm0\n\t"\
			"addpd	%%xmm3		,%%xmm1\n\t"\
			"movaps	%%xmm0		,    (%%edi)	/* <- t06 */\n\t"\
			"movaps	%%xmm1		,0x10(%%edi)	/* <- t07 */\n\t"\
		"addl	$0x40		,%%esi\n\t"\
		"addl	$0x40		,%%edi\n\t"\
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
			"movaps	%%xmm2		,    (%%esi)	/* <- t08 */\n\t"\
			"movaps	%%xmm3		,0x10(%%esi)	/* <- t09 */\n\t"\
			"movaps	%%xmm0		,    (%%edi)	/* <- t0a */\n\t"\
			"movaps	%%xmm1		,0x10(%%edi)	/* <- t0b */\n\t"\
		"addl	$0x40		,%%esi\n\t"\
			"subl	$0x080		,%%eax 		/* __i2 = s1p28r */\n\t"\
			"subl	$0x080		,%%ebx 		/* __i5 = s1p16r */\n\t"\
			"subl	$0x080		,%%ecx 		/* __i8 = s1p04r */\n\t"\
			"movaps	    (%%ebx)	,%%xmm4\n\t"\
			"movaps	0x10(%%ebx)	,%%xmm5\n\t"\
			"movaps	    (%%eax)	,%%xmm0\n\t"\
			"movaps	0x10(%%eax)	,%%xmm1\n\t"\
			"movaps	    (%%ecx)	,%%xmm2\n\t"\
			"movaps	0x10(%%ecx)	,%%xmm3\n\t"\
			"subpd	%%xmm2		,%%xmm4\n\t"\
			"subpd	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm2\n\t"\
			"addpd	%%xmm3		,%%xmm3\n\t"\
			"addpd	%%xmm4		,%%xmm2\n\t"\
			"addpd	%%xmm5		,%%xmm3\n\t"\
			"addpd	%%xmm2		,%%xmm0\n\t"\
			"addpd	%%xmm3		,%%xmm1\n\t"\
			"movaps	%%xmm0		,    (%%esi)	/* <- t0c */\n\t"\
			"movaps	%%xmm1		,0x10(%%esi)	/* <- t0d */\n\t"\
		"addl	$0x40		,%%edi\n\t"\
		"addl	$0x40		,%%esi\n\t"\
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
			"movaps	%%xmm2		,    (%%edi)	/* <- t0e */\n\t"\
			"movaps	%%xmm3		,0x10(%%edi)	/* <- t0f */\n\t"\
			"movaps	%%xmm0		,    (%%esi)	/* <- t0g */\n\t"\
			"movaps	%%xmm1		,0x10(%%esi)	/* <- t0h */\n\t"\
			"/******************************************************************************/\n\t"\
			"/*...and now do three more radix-3 transforms, including the twiddle factors: */\n\t"\
			"/******************************************************************************/\n\t"\
			"movl	%[__r00]	,%%eax 		/* __r00 */\n\t"\
			"movl	%%eax		,%%ebx\n\t"\
			"movl	%%eax		,%%ecx\n\t"\
			"addl	$0x60		,%%ebx 		/* __r06 */\n\t"\
			"addl	$0xc0		,%%ecx 		/* __r0c */\n\t"\
			"movaps	    (%%ebx)	,%%xmm4\n\t"\
			"movaps	0x10(%%ebx)	,%%xmm5\n\t"\
			"movaps	    (%%ecx)	,%%xmm2\n\t"\
			"movaps	0x10(%%ecx)	,%%xmm3\n\t"\
			"movaps	    (%%eax)	,%%xmm0\n\t"\
			"movaps	0x10(%%eax)	,%%xmm1\n\t"\
			"subpd	%%xmm2		,%%xmm4\n\t"\
			"subpd	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm2\n\t"\
			"addpd	%%xmm3		,%%xmm3\n\t"\
			"addpd	%%xmm4		,%%xmm2\n\t"\
			"addpd	%%xmm5		,%%xmm3\n\t"\
			"addpd	%%xmm2		,%%xmm0\n\t"\
			"addpd	%%xmm3		,%%xmm1\n\t"\
			"movaps	%%xmm0		,    (%%eax)\n\t"\
			"movaps	%%xmm1		,0x10(%%eax)\n\t"\
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
			"movaps	%%xmm2		,    (%%ebx)\n\t"\
			"movaps	%%xmm3		,0x10(%%ebx)\n\t"\
			"movaps	%%xmm0		,    (%%ecx)\n\t"\
			"movaps	%%xmm1		,0x10(%%ecx)\n\t"\
		"subl	$0x40		,%%edx 		/* c1 */\n\t"\
			"addl	$0x20		,%%eax 		/* __r02 */\n\t"\
			"addl	$0x20		,%%ebx 		/* __r08 */\n\t"\
			"addl	$0x20		,%%ecx 		/* __r0e */\n\t"\
			"movaps	    (%%ebx)	,%%xmm2\n\t"\
			"movaps	0x10(%%ebx)	,%%xmm3\n\t"\
			"movaps	    (%%edx)	,%%xmm6\n\t"\
			"movaps	0x10(%%edx)	,%%xmm7\n\t"\
		"addl	$0x20		,%%edx 		/* c2 */\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"subpd	%%xmm1		,%%xmm2\n\t"\
			"addpd	%%xmm0		,%%xmm3\n\t"\
			"movaps	    (%%ecx)	,%%xmm4\n\t"\
			"movaps	0x10(%%ecx)	,%%xmm5\n\t"\
			"movaps	    (%%edx)	,%%xmm6\n\t"\
			"movaps	0x10(%%edx)	,%%xmm7\n\t"\
		"addl	$0x20		,%%edx 		/* c3m1 */\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"subpd	%%xmm1		,%%xmm4\n\t"\
			"addpd	%%xmm0		,%%xmm5\n\t"\
			"movaps	    (%%edx)	,%%xmm6\n\t"\
			"movaps	0x10(%%edx)	,%%xmm7\n\t"\
			"movaps	    (%%eax)	,%%xmm0\n\t"\
			"movaps	0x10(%%eax)	,%%xmm1\n\t"\
			"subpd	%%xmm4		,%%xmm2\n\t"\
			"subpd	%%xmm5		,%%xmm3\n\t"\
			"addpd	%%xmm4		,%%xmm4\n\t"\
			"addpd	%%xmm5		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm4\n\t"\
			"addpd	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm4		,%%xmm0\n\t"\
			"addpd	%%xmm5		,%%xmm1\n\t"\
			"movaps	%%xmm0		,    (%%eax)\n\t"\
			"movaps	%%xmm1		,0x10(%%eax)\n\t"\
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
			"movaps	%%xmm4		,    (%%ebx)\n\t"\
			"movaps	%%xmm5		,0x10(%%ebx)\n\t"\
			"movaps	%%xmm0		,    (%%ecx)\n\t"\
			"movaps	%%xmm1		,0x10(%%ecx)\n\t"\
		"subl	$0x20		,%%edx 		/* c2 */\n\t"\
			"addl	$0x20		,%%eax 		/* __r04 */\n\t"\
			"addl	$0x20		,%%ebx 		/* __r0a */\n\t"\
			"addl	$0x20		,%%ecx 		/* __r0g */\n\t"\
			"movaps	    (%%ebx)	,%%xmm2\n\t"\
			"movaps	0x10(%%ebx)	,%%xmm3\n\t"\
			"movaps	    (%%edx)	,%%xmm6\n\t"\
			"movaps	0x10(%%edx)	,%%xmm7\n\t"\
		"addl	$0x40		,%%edx 		/* c4 */\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"subpd	%%xmm1		,%%xmm2\n\t"\
			"addpd	%%xmm0		,%%xmm3\n\t"\
			"movaps	    (%%ecx)	,%%xmm4\n\t"\
			"movaps	0x10(%%ecx)	,%%xmm5\n\t"\
			"movaps	    (%%edx)	,%%xmm6\n\t"\
			"movaps	0x10(%%edx)	,%%xmm7\n\t"\
		"subl	$0x20		,%%edx 		/* c3m1 */\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"subpd	%%xmm1		,%%xmm4\n\t"\
			"addpd	%%xmm0		,%%xmm5\n\t"\
			"movaps	    (%%edx)	,%%xmm6\n\t"\
			"movaps	0x10(%%edx)	,%%xmm7\n\t"\
			"movaps	    (%%eax)	,%%xmm0\n\t"\
			"movaps	0x10(%%eax)	,%%xmm1\n\t"\
			"\n\t"\
			"subpd	%%xmm4		,%%xmm2\n\t"\
			"subpd	%%xmm5		,%%xmm3\n\t"\
			"addpd	%%xmm4		,%%xmm4\n\t"\
			"addpd	%%xmm5		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm4\n\t"\
			"addpd	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm4		,%%xmm0\n\t"\
			"addpd	%%xmm5		,%%xmm1\n\t"\
			"movaps	%%xmm0		,    (%%eax)\n\t"\
			"movaps	%%xmm1		,0x10(%%eax)\n\t"\
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
			"movaps	%%xmm4		,    (%%ebx)\n\t"\
			"movaps	%%xmm5		,0x10(%%ebx)\n\t"\
			"movaps	%%xmm0		,    (%%ecx)\n\t"\
			"movaps	%%xmm1		,0x10(%%ecx)\n\t"\
			"/* SSE2_RADIX_09_DIF_0TWIDDLE(r10,s1p27r,s1p23r,s1p19r,s1p15r,s1p11r,s1p07r,s1p03r,s1p35r,s1p31r) */\n\t"\
			"movl	%[__in0]	,%%eax 		/* __i0-8 = s1p[27,23,19,15,11,07,03,35,31]r; e[abc]x store output addresses throughout */\n\t"\
			"movl	%[__r00]	,%%esi\n\t"\
			"addl	$0x360		,%%eax 		/* __i0 = s1p27r */\n\t"\
			"addl	$0x120		,%%esi 		/* r10 */\n\t"\
		"movl	%[__cc1]	,%%edx 		/* edx stores trig addresses throughout */\n\t"\
			"movl	%%eax		,%%ebx\n\t"\
			"movl	%%eax		,%%ecx\n\t"\
			"subl	$0x180		,%%ebx 		/* __i3 = s1p15r */\n\t"\
			"subl	$0x300		,%%ecx 		/* __i6 = s1p03r */\n\t"\
		"addl	$0x40		,%%edx 		/* c3m1 */\n\t"\
			"movaps	    (%%ebx)	,%%xmm2\n\t"\
			"movaps	0x10(%%ebx)	,%%xmm3\n\t"\
			"movaps	    (%%eax)	,%%xmm0\n\t"\
			"movaps	0x10(%%eax)	,%%xmm1\n\t"\
			"movaps	    (%%ecx)	,%%xmm6\n\t"\
			"movaps	0x10(%%ecx)	,%%xmm7\n\t"\
			"movaps	%%xmm2		,%%xmm4\n\t"\
			"movaps	%%xmm3		,%%xmm5\n\t"\
		"movl	%%esi		,%%edi\n\t"\
		"addl	$0x20		,%%edi\n\t"\
			"addpd	%%xmm6		,%%xmm2\n\t"\
			"addpd	%%xmm7		,%%xmm3\n\t"\
			"subpd	%%xmm6		,%%xmm4\n\t"\
			"subpd	%%xmm7		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm0\n\t"\
			"addpd	%%xmm3		,%%xmm1\n\t"\
			"movaps	    (%%edx)	,%%xmm6\n\t"\
			"movaps	0x10(%%edx)	,%%xmm7\n\t"\
			"movaps	%%xmm0		,    (%%esi)	/* <- t00 */\n\t"\
			"movaps	%%xmm1		,0x10(%%esi)	/* <- t01 */\n\t"\
		"addl	$0x40		,%%esi\n\t"\
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
			"movaps	%%xmm2		,    (%%edi)	/* <- t02 */\n\t"\
			"movaps	%%xmm3		,0x10(%%edi)	/* <- t03 */\n\t"\
			"movaps	%%xmm0		,    (%%esi)	/* <- t04 */\n\t"\
			"movaps	%%xmm1		,0x10(%%esi)	/* <- t05 */\n\t"\
			"subl	$0x080		,%%eax 		/* __i1 = s1p23r */\n\t"\
			"subl	$0x080		,%%ebx 		/* __i4 = s1p11r */\n\t"\
			"addl	$0x400		,%%ecx 		/* __i7 = s1p35r */\n\t"\
			"movaps	    (%%ebx)	,%%xmm4\n\t"\
			"movaps	0x10(%%ebx)	,%%xmm5\n\t"\
			"movaps	    (%%eax)	,%%xmm0\n\t"\
			"movaps	0x10(%%eax)	,%%xmm1\n\t"\
			"movaps	    (%%ecx)	,%%xmm2\n\t"\
			"movaps	0x10(%%ecx)	,%%xmm3\n\t"\
		"addl	$0x40		,%%edi\n\t"\
			"subpd	%%xmm2		,%%xmm4\n\t"\
			"subpd	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm2\n\t"\
			"addpd	%%xmm3		,%%xmm3\n\t"\
			"addpd	%%xmm4		,%%xmm2\n\t"\
			"addpd	%%xmm5		,%%xmm3\n\t"\
			"addpd	%%xmm2		,%%xmm0\n\t"\
			"addpd	%%xmm3		,%%xmm1\n\t"\
			"movaps	%%xmm0		,    (%%edi)	/* <- t06 */\n\t"\
			"movaps	%%xmm1		,0x10(%%edi)	/* <- t07 */\n\t"\
		"addl	$0x40		,%%esi\n\t"\
		"addl	$0x40		,%%edi\n\t"\
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
			"movaps	%%xmm2		,    (%%esi)	/* <- t08 */\n\t"\
			"movaps	%%xmm3		,0x10(%%esi)	/* <- t09 */\n\t"\
			"movaps	%%xmm0		,    (%%edi)	/* <- t0a */\n\t"\
			"movaps	%%xmm1		,0x10(%%edi)	/* <- t0b */\n\t"\
		"addl	$0x40		,%%esi\n\t"\
			"subl	$0x080		,%%eax 		/* __i2 = s1p19r */\n\t"\
			"subl	$0x080		,%%ebx 		/* __i5 = s1p07r */\n\t"\
			"subl	$0x080		,%%ecx 		/* __i8 = s1p31r */\n\t"\
			"movaps	    (%%ebx)	,%%xmm4\n\t"\
			"movaps	0x10(%%ebx)	,%%xmm5\n\t"\
			"movaps	    (%%eax)	,%%xmm0\n\t"\
			"movaps	0x10(%%eax)	,%%xmm1\n\t"\
			"movaps	    (%%ecx)	,%%xmm2\n\t"\
			"movaps	0x10(%%ecx)	,%%xmm3\n\t"\
			"subpd	%%xmm2		,%%xmm4\n\t"\
			"subpd	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm2\n\t"\
			"addpd	%%xmm3		,%%xmm3\n\t"\
			"addpd	%%xmm4		,%%xmm2\n\t"\
			"addpd	%%xmm5		,%%xmm3\n\t"\
			"addpd	%%xmm2		,%%xmm0\n\t"\
			"addpd	%%xmm3		,%%xmm1\n\t"\
			"movaps	%%xmm0		,    (%%esi)	/* <- t0c */\n\t"\
			"movaps	%%xmm1		,0x10(%%esi)	/* <- t0d */\n\t"\
		"addl	$0x40		,%%edi\n\t"\
		"addl	$0x40		,%%esi\n\t"\
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
			"movaps	%%xmm2		,    (%%edi)	/* <- t0e */\n\t"\
			"movaps	%%xmm3		,0x10(%%edi)	/* <- t0f */\n\t"\
			"movaps	%%xmm0		,    (%%esi)	/* <- t0g */\n\t"\
			"movaps	%%xmm1		,0x10(%%esi)	/* <- t0h */\n\t"\
			"/******************************************************************************/\n\t"\
			"/*...and now do three more radix-3 transforms, including the twiddle factors: */\n\t"\
			"/******************************************************************************/\n\t"\
			"movl	%[__r00]	,%%eax 		/* __r00 */\n\t"\
			"addl	$0x120		,%%eax 		/* r10 */\n\t"\
			"movl	%%eax		,%%ebx\n\t"\
			"movl	%%eax		,%%ecx\n\t"\
			"addl	$0x60		,%%ebx 		/* __r16 */\n\t"\
			"addl	$0xc0		,%%ecx 		/* __r1c */\n\t"\
			"movaps	    (%%ebx)	,%%xmm4\n\t"\
			"movaps	0x10(%%ebx)	,%%xmm5\n\t"\
			"movaps	    (%%ecx)	,%%xmm2\n\t"\
			"movaps	0x10(%%ecx)	,%%xmm3\n\t"\
			"movaps	    (%%eax)	,%%xmm0\n\t"\
			"movaps	0x10(%%eax)	,%%xmm1\n\t"\
			"subpd	%%xmm2		,%%xmm4\n\t"\
			"subpd	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm2\n\t"\
			"addpd	%%xmm3		,%%xmm3\n\t"\
			"addpd	%%xmm4		,%%xmm2\n\t"\
			"addpd	%%xmm5		,%%xmm3\n\t"\
			"addpd	%%xmm2		,%%xmm0\n\t"\
			"addpd	%%xmm3		,%%xmm1\n\t"\
			"movaps	%%xmm0		,    (%%eax)\n\t"\
			"movaps	%%xmm1		,0x10(%%eax)\n\t"\
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
			"movaps	%%xmm2		,    (%%ebx)\n\t"\
			"movaps	%%xmm3		,0x10(%%ebx)\n\t"\
			"movaps	%%xmm0		,    (%%ecx)\n\t"\
			"movaps	%%xmm1		,0x10(%%ecx)\n\t"\
		"subl	$0x40		,%%edx 		/* c1 */\n\t"\
			"addl	$0x20		,%%eax 		/* __r12 */\n\t"\
			"addl	$0x20		,%%ebx 		/* __r18 */\n\t"\
			"addl	$0x20		,%%ecx 		/* __r1e */\n\t"\
			"movaps	    (%%ebx)	,%%xmm2\n\t"\
			"movaps	0x10(%%ebx)	,%%xmm3\n\t"\
			"movaps	    (%%edx)	,%%xmm6\n\t"\
			"movaps	0x10(%%edx)	,%%xmm7\n\t"\
		"addl	$0x20		,%%edx 		/* c2 */\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"subpd	%%xmm1		,%%xmm2\n\t"\
			"addpd	%%xmm0		,%%xmm3\n\t"\
			"movaps	    (%%ecx)	,%%xmm4\n\t"\
			"movaps	0x10(%%ecx)	,%%xmm5\n\t"\
			"movaps	    (%%edx)	,%%xmm6\n\t"\
			"movaps	0x10(%%edx)	,%%xmm7\n\t"\
		"addl	$0x20		,%%edx 		/* c3m1 */\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"subpd	%%xmm1		,%%xmm4\n\t"\
			"addpd	%%xmm0		,%%xmm5\n\t"\
			"movaps	    (%%edx)	,%%xmm6\n\t"\
			"movaps	0x10(%%edx)	,%%xmm7\n\t"\
			"movaps	    (%%eax)	,%%xmm0\n\t"\
			"movaps	0x10(%%eax)	,%%xmm1\n\t"\
			"subpd	%%xmm4		,%%xmm2\n\t"\
			"subpd	%%xmm5		,%%xmm3\n\t"\
			"addpd	%%xmm4		,%%xmm4\n\t"\
			"addpd	%%xmm5		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm4\n\t"\
			"addpd	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm4		,%%xmm0\n\t"\
			"addpd	%%xmm5		,%%xmm1\n\t"\
			"movaps	%%xmm0		,    (%%eax)\n\t"\
			"movaps	%%xmm1		,0x10(%%eax)\n\t"\
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
			"movaps	%%xmm4		,    (%%ebx)\n\t"\
			"movaps	%%xmm5		,0x10(%%ebx)\n\t"\
			"movaps	%%xmm0		,    (%%ecx)\n\t"\
			"movaps	%%xmm1		,0x10(%%ecx)\n\t"\
		"subl	$0x20		,%%edx 		/* c2 */\n\t"\
			"addl	$0x20		,%%eax 		/* __r14 */\n\t"\
			"addl	$0x20		,%%ebx 		/* __r1a */\n\t"\
			"addl	$0x20		,%%ecx 		/* __r1g */\n\t"\
			"movaps	    (%%ebx)	,%%xmm2\n\t"\
			"movaps	0x10(%%ebx)	,%%xmm3\n\t"\
			"movaps	    (%%edx)	,%%xmm6\n\t"\
			"movaps	0x10(%%edx)	,%%xmm7\n\t"\
		"addl	$0x40		,%%edx 		/* c4 */\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"subpd	%%xmm1		,%%xmm2\n\t"\
			"addpd	%%xmm0		,%%xmm3\n\t"\
			"movaps	    (%%ecx)	,%%xmm4\n\t"\
			"movaps	0x10(%%ecx)	,%%xmm5\n\t"\
			"movaps	    (%%edx)	,%%xmm6\n\t"\
			"movaps	0x10(%%edx)	,%%xmm7\n\t"\
		"subl	$0x20		,%%edx 		/* c3m1 */\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"subpd	%%xmm1		,%%xmm4\n\t"\
			"addpd	%%xmm0		,%%xmm5\n\t"\
			"movaps	    (%%edx)	,%%xmm6\n\t"\
			"movaps	0x10(%%edx)	,%%xmm7\n\t"\
			"movaps	    (%%eax)	,%%xmm0\n\t"\
			"movaps	0x10(%%eax)	,%%xmm1\n\t"\
			"\n\t"\
			"subpd	%%xmm4		,%%xmm2\n\t"\
			"subpd	%%xmm5		,%%xmm3\n\t"\
			"addpd	%%xmm4		,%%xmm4\n\t"\
			"addpd	%%xmm5		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm4\n\t"\
			"addpd	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm4		,%%xmm0\n\t"\
			"addpd	%%xmm5		,%%xmm1\n\t"\
			"movaps	%%xmm0		,    (%%eax)\n\t"\
			"movaps	%%xmm1		,0x10(%%eax)\n\t"\
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
			"movaps	%%xmm4		,    (%%ebx)\n\t"\
			"movaps	%%xmm5		,0x10(%%ebx)\n\t"\
			"movaps	%%xmm0		,    (%%ecx)\n\t"\
			"movaps	%%xmm1		,0x10(%%ecx)\n\t"\
			"/* SSE2_RADIX_09_DIF_0TWIDDLE(r20,s1p18r,s1p14r,s1p10r,s1p06r,s1p02r,s1p34r,s1p30r,s1p26r,s1p22r) */\n\t"\
			"movl	%[__in0]	,%%eax 		/* __i0-8 = s1p[18,14,10,06,02,34,30,26,22]r; e[abc]x store output addresses throughout */\n\t"\
			"movl	%[__r00]	,%%esi\n\t"\
			"addl	$0x240		,%%eax 		/* __i0 = s1p18r */\n\t"\
			"addl	$0x240		,%%esi 		/* r20 */\n\t"\
		"movl	%[__cc1]	,%%edx 		/* edx stores trig addresses throughout */\n\t"\
			"movl	%%eax		,%%ebx\n\t"\
			"movl	%%eax		,%%ecx\n\t"\
			"subl	$0x180		,%%ebx 		/* __i3 = s1p06r */\n\t"\
			"addl	$0x180		,%%ecx 		/* __i6 = s1p30r */\n\t"\
		"addl	$0x40		,%%edx 		/* c3m1 */\n\t"\
			"movaps	    (%%ebx)	,%%xmm2\n\t"\
			"movaps	0x10(%%ebx)	,%%xmm3\n\t"\
			"movaps	    (%%eax)	,%%xmm0\n\t"\
			"movaps	0x10(%%eax)	,%%xmm1\n\t"\
			"movaps	    (%%ecx)	,%%xmm6\n\t"\
			"movaps	0x10(%%ecx)	,%%xmm7\n\t"\
			"movaps	%%xmm2		,%%xmm4\n\t"\
			"movaps	%%xmm3		,%%xmm5\n\t"\
		"movl	%%esi		,%%edi\n\t"\
		"addl	$0x20		,%%edi\n\t"\
			"addpd	%%xmm6		,%%xmm2\n\t"\
			"addpd	%%xmm7		,%%xmm3\n\t"\
			"subpd	%%xmm6		,%%xmm4\n\t"\
			"subpd	%%xmm7		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm0\n\t"\
			"addpd	%%xmm3		,%%xmm1\n\t"\
			"movaps	    (%%edx)	,%%xmm6\n\t"\
			"movaps	0x10(%%edx)	,%%xmm7\n\t"\
			"movaps	%%xmm0		,    (%%esi)	/* <- t00 */\n\t"\
			"movaps	%%xmm1		,0x10(%%esi)	/* <- t01 */\n\t"\
		"addl	$0x40		,%%esi\n\t"\
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
			"movaps	%%xmm2		,    (%%edi)	/* <- t02 */\n\t"\
			"movaps	%%xmm3		,0x10(%%edi)	/* <- t03 */\n\t"\
			"movaps	%%xmm0		,    (%%esi)	/* <- t04 */\n\t"\
			"movaps	%%xmm1		,0x10(%%esi)	/* <- t05 */\n\t"\
			"subl	$0x080		,%%eax 		/* __i1 = s1p14r */\n\t"\
			"subl	$0x080		,%%ebx 		/* __i4 = s1p02r */\n\t"\
			"subl	$0x080		,%%ecx 		/* __i7 = s1p26r */\n\t"\
			"movaps	    (%%ebx)	,%%xmm4\n\t"\
			"movaps	0x10(%%ebx)	,%%xmm5\n\t"\
			"movaps	    (%%eax)	,%%xmm0\n\t"\
			"movaps	0x10(%%eax)	,%%xmm1\n\t"\
			"movaps	    (%%ecx)	,%%xmm2\n\t"\
			"movaps	0x10(%%ecx)	,%%xmm3\n\t"\
		"addl	$0x40		,%%edi\n\t"\
			"subpd	%%xmm2		,%%xmm4\n\t"\
			"subpd	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm2\n\t"\
			"addpd	%%xmm3		,%%xmm3\n\t"\
			"addpd	%%xmm4		,%%xmm2\n\t"\
			"addpd	%%xmm5		,%%xmm3\n\t"\
			"addpd	%%xmm2		,%%xmm0\n\t"\
			"addpd	%%xmm3		,%%xmm1\n\t"\
			"movaps	%%xmm0		,    (%%edi)	/* <- t06 */\n\t"\
			"movaps	%%xmm1		,0x10(%%edi)	/* <- t07 */\n\t"\
		"addl	$0x40		,%%esi\n\t"\
		"addl	$0x40		,%%edi\n\t"\
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
			"movaps	%%xmm2		,    (%%esi)	/* <- t08 */\n\t"\
			"movaps	%%xmm3		,0x10(%%esi)	/* <- t09 */\n\t"\
			"movaps	%%xmm0		,    (%%edi)	/* <- t0a */\n\t"\
			"movaps	%%xmm1		,0x10(%%edi)	/* <- t0b */\n\t"\
		"addl	$0x40		,%%esi\n\t"\
			"subl	$0x080		,%%eax 		/* __i2 = s1p10r */\n\t"\
			"addl	$0x400		,%%ebx 		/* __i5 = s1p34r */\n\t"\
			"subl	$0x080		,%%ecx 		/* __i8 = s1p22r */\n\t"\
			"movaps	    (%%ebx)	,%%xmm4\n\t"\
			"movaps	0x10(%%ebx)	,%%xmm5\n\t"\
			"movaps	    (%%eax)	,%%xmm0\n\t"\
			"movaps	0x10(%%eax)	,%%xmm1\n\t"\
			"movaps	    (%%ecx)	,%%xmm2\n\t"\
			"movaps	0x10(%%ecx)	,%%xmm3\n\t"\
			"subpd	%%xmm2		,%%xmm4\n\t"\
			"subpd	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm2\n\t"\
			"addpd	%%xmm3		,%%xmm3\n\t"\
			"addpd	%%xmm4		,%%xmm2\n\t"\
			"addpd	%%xmm5		,%%xmm3\n\t"\
			"addpd	%%xmm2		,%%xmm0\n\t"\
			"addpd	%%xmm3		,%%xmm1\n\t"\
			"movaps	%%xmm0		,    (%%esi)	/* <- t0c */\n\t"\
			"movaps	%%xmm1		,0x10(%%esi)	/* <- t0d */\n\t"\
		"addl	$0x40		,%%edi\n\t"\
		"addl	$0x40		,%%esi\n\t"\
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
			"movaps	%%xmm2		,    (%%edi)	/* <- t0e */\n\t"\
			"movaps	%%xmm3		,0x10(%%edi)	/* <- t0f */\n\t"\
			"movaps	%%xmm0		,    (%%esi)	/* <- t0g */\n\t"\
			"movaps	%%xmm1		,0x10(%%esi)	/* <- t0h */\n\t"\
			"/******************************************************************************/\n\t"\
			"/*...and now do three more radix-3 transforms, including the twiddle factors: */\n\t"\
			"/******************************************************************************/\n\t"\
			"movl	%[__r00]	,%%eax 		/* __r00 */\n\t"\
			"addl	$0x240		,%%eax 		/* r20 */\n\t"\
			"movl	%%eax		,%%ebx\n\t"\
			"movl	%%eax		,%%ecx\n\t"\
			"addl	$0x60		,%%ebx 		/* __r26 */\n\t"\
			"addl	$0xc0		,%%ecx 		/* __r2c */\n\t"\
			"movaps	    (%%ebx)	,%%xmm4\n\t"\
			"movaps	0x10(%%ebx)	,%%xmm5\n\t"\
			"movaps	    (%%ecx)	,%%xmm2\n\t"\
			"movaps	0x10(%%ecx)	,%%xmm3\n\t"\
			"movaps	    (%%eax)	,%%xmm0\n\t"\
			"movaps	0x10(%%eax)	,%%xmm1\n\t"\
			"subpd	%%xmm2		,%%xmm4\n\t"\
			"subpd	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm2\n\t"\
			"addpd	%%xmm3		,%%xmm3\n\t"\
			"addpd	%%xmm4		,%%xmm2\n\t"\
			"addpd	%%xmm5		,%%xmm3\n\t"\
			"addpd	%%xmm2		,%%xmm0\n\t"\
			"addpd	%%xmm3		,%%xmm1\n\t"\
			"movaps	%%xmm0		,    (%%eax)\n\t"\
			"movaps	%%xmm1		,0x10(%%eax)\n\t"\
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
			"movaps	%%xmm2		,    (%%ebx)\n\t"\
			"movaps	%%xmm3		,0x10(%%ebx)\n\t"\
			"movaps	%%xmm0		,    (%%ecx)\n\t"\
			"movaps	%%xmm1		,0x10(%%ecx)\n\t"\
		"subl	$0x40		,%%edx 		/* c1 */\n\t"\
			"addl	$0x20		,%%eax 		/* __r22 */\n\t"\
			"addl	$0x20		,%%ebx 		/* __r28 */\n\t"\
			"addl	$0x20		,%%ecx 		/* __r2e */\n\t"\
			"movaps	    (%%ebx)	,%%xmm2\n\t"\
			"movaps	0x10(%%ebx)	,%%xmm3\n\t"\
			"movaps	    (%%edx)	,%%xmm6\n\t"\
			"movaps	0x10(%%edx)	,%%xmm7\n\t"\
		"addl	$0x20		,%%edx 		/* c2 */\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"subpd	%%xmm1		,%%xmm2\n\t"\
			"addpd	%%xmm0		,%%xmm3\n\t"\
			"movaps	    (%%ecx)	,%%xmm4\n\t"\
			"movaps	0x10(%%ecx)	,%%xmm5\n\t"\
			"movaps	    (%%edx)	,%%xmm6\n\t"\
			"movaps	0x10(%%edx)	,%%xmm7\n\t"\
		"addl	$0x20		,%%edx 		/* c3m1 */\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"subpd	%%xmm1		,%%xmm4\n\t"\
			"addpd	%%xmm0		,%%xmm5\n\t"\
			"movaps	    (%%edx)	,%%xmm6\n\t"\
			"movaps	0x10(%%edx)	,%%xmm7\n\t"\
			"movaps	    (%%eax)	,%%xmm0\n\t"\
			"movaps	0x10(%%eax)	,%%xmm1\n\t"\
			"subpd	%%xmm4		,%%xmm2\n\t"\
			"subpd	%%xmm5		,%%xmm3\n\t"\
			"addpd	%%xmm4		,%%xmm4\n\t"\
			"addpd	%%xmm5		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm4\n\t"\
			"addpd	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm4		,%%xmm0\n\t"\
			"addpd	%%xmm5		,%%xmm1\n\t"\
			"movaps	%%xmm0		,    (%%eax)\n\t"\
			"movaps	%%xmm1		,0x10(%%eax)\n\t"\
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
			"movaps	%%xmm4		,    (%%ebx)\n\t"\
			"movaps	%%xmm5		,0x10(%%ebx)\n\t"\
			"movaps	%%xmm0		,    (%%ecx)\n\t"\
			"movaps	%%xmm1		,0x10(%%ecx)\n\t"\
		"subl	$0x20		,%%edx 		/* c2 */\n\t"\
			"addl	$0x20		,%%eax 		/* __r24 */\n\t"\
			"addl	$0x20		,%%ebx 		/* __r2a */\n\t"\
			"addl	$0x20		,%%ecx 		/* __r2g */\n\t"\
			"movaps	    (%%ebx)	,%%xmm2\n\t"\
			"movaps	0x10(%%ebx)	,%%xmm3\n\t"\
			"movaps	    (%%edx)	,%%xmm6\n\t"\
			"movaps	0x10(%%edx)	,%%xmm7\n\t"\
		"addl	$0x40		,%%edx 		/* c4 */\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"subpd	%%xmm1		,%%xmm2\n\t"\
			"addpd	%%xmm0		,%%xmm3\n\t"\
			"movaps	    (%%ecx)	,%%xmm4\n\t"\
			"movaps	0x10(%%ecx)	,%%xmm5\n\t"\
			"movaps	    (%%edx)	,%%xmm6\n\t"\
			"movaps	0x10(%%edx)	,%%xmm7\n\t"\
		"subl	$0x20		,%%edx 		/* c3m1 */\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"subpd	%%xmm1		,%%xmm4\n\t"\
			"addpd	%%xmm0		,%%xmm5\n\t"\
			"movaps	    (%%edx)	,%%xmm6\n\t"\
			"movaps	0x10(%%edx)	,%%xmm7\n\t"\
			"movaps	    (%%eax)	,%%xmm0\n\t"\
			"movaps	0x10(%%eax)	,%%xmm1\n\t"\
			"\n\t"\
			"subpd	%%xmm4		,%%xmm2\n\t"\
			"subpd	%%xmm5		,%%xmm3\n\t"\
			"addpd	%%xmm4		,%%xmm4\n\t"\
			"addpd	%%xmm5		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm4\n\t"\
			"addpd	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm4		,%%xmm0\n\t"\
			"addpd	%%xmm5		,%%xmm1\n\t"\
			"movaps	%%xmm0		,    (%%eax)\n\t"\
			"movaps	%%xmm1		,0x10(%%eax)\n\t"\
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
			"movaps	%%xmm4		,    (%%ebx)\n\t"\
			"movaps	%%xmm5		,0x10(%%ebx)\n\t"\
			"movaps	%%xmm0		,    (%%ecx)\n\t"\
			"movaps	%%xmm1		,0x10(%%ecx)\n\t"\
			"/* SSE2_RADIX_09_DIF_0TWIDDLE(r30,s1p09r,s1p05r,s1p01r,s1p33r,s1p29r,s1p25r,s1p21r,s1p17r,s1p13r) */\n\t"\
			"movl	%[__in0]	,%%eax 		/* __i0-8 = s1p[09,05,01,33,29,25,21,17,13]r; e[abc]x store output addresses throughout */\n\t"\
			"movl	%[__r00]	,%%esi\n\t"\
			"addl	$0x120		,%%eax 		/* __i0 = s1p09r */\n\t"\
			"addl	$0x360		,%%esi 		/* r30 */\n\t"\
		"movl	%[__cc1]	,%%edx 		/* edx stores trig addresses throughout */\n\t"\
			"movl	%%eax		,%%ebx\n\t"\
			"movl	%%eax		,%%ecx\n\t"\
			"addl	$0x300		,%%ebx 		/* __i3 = s1p33r */\n\t"\
			"addl	$0x180		,%%ecx 		/* __i6 = s1p21r */\n\t"\
		"addl	$0x40		,%%edx 		/* c3m1 */\n\t"\
			"movaps	    (%%ebx)	,%%xmm2\n\t"\
			"movaps	0x10(%%ebx)	,%%xmm3\n\t"\
			"movaps	    (%%eax)	,%%xmm0\n\t"\
			"movaps	0x10(%%eax)	,%%xmm1\n\t"\
			"movaps	    (%%ecx)	,%%xmm6\n\t"\
			"movaps	0x10(%%ecx)	,%%xmm7\n\t"\
			"movaps	%%xmm2		,%%xmm4\n\t"\
			"movaps	%%xmm3		,%%xmm5\n\t"\
		"movl	%%esi		,%%edi\n\t"\
		"addl	$0x20		,%%edi\n\t"\
			"addpd	%%xmm6		,%%xmm2\n\t"\
			"addpd	%%xmm7		,%%xmm3\n\t"\
			"subpd	%%xmm6		,%%xmm4\n\t"\
			"subpd	%%xmm7		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm0\n\t"\
			"addpd	%%xmm3		,%%xmm1\n\t"\
			"movaps	    (%%edx)	,%%xmm6\n\t"\
			"movaps	0x10(%%edx)	,%%xmm7\n\t"\
			"movaps	%%xmm0		,    (%%esi)	/* <- t00 */\n\t"\
			"movaps	%%xmm1		,0x10(%%esi)	/* <- t01 */\n\t"\
		"addl	$0x40		,%%esi\n\t"\
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
			"movaps	%%xmm2		,    (%%edi)	/* <- t02 */\n\t"\
			"movaps	%%xmm3		,0x10(%%edi)	/* <- t03 */\n\t"\
			"movaps	%%xmm0		,    (%%esi)	/* <- t04 */\n\t"\
			"movaps	%%xmm1		,0x10(%%esi)	/* <- t05 */\n\t"\
			"subl	$0x080		,%%eax 		/* __i1 = s1p05r */\n\t"\
			"subl	$0x080		,%%ebx 		/* __i4 = s1p29r */\n\t"\
			"subl	$0x080		,%%ecx 		/* __i7 = s1p17r */\n\t"\
			"movaps	    (%%ebx)	,%%xmm4\n\t"\
			"movaps	0x10(%%ebx)	,%%xmm5\n\t"\
			"movaps	    (%%eax)	,%%xmm0\n\t"\
			"movaps	0x10(%%eax)	,%%xmm1\n\t"\
			"movaps	    (%%ecx)	,%%xmm2\n\t"\
			"movaps	0x10(%%ecx)	,%%xmm3\n\t"\
		"addl	$0x40		,%%edi\n\t"\
			"subpd	%%xmm2		,%%xmm4\n\t"\
			"subpd	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm2\n\t"\
			"addpd	%%xmm3		,%%xmm3\n\t"\
			"addpd	%%xmm4		,%%xmm2\n\t"\
			"addpd	%%xmm5		,%%xmm3\n\t"\
			"addpd	%%xmm2		,%%xmm0\n\t"\
			"addpd	%%xmm3		,%%xmm1\n\t"\
			"movaps	%%xmm0		,    (%%edi)	/* <- t06 */\n\t"\
			"movaps	%%xmm1		,0x10(%%edi)	/* <- t07 */\n\t"\
		"addl	$0x40		,%%esi\n\t"\
		"addl	$0x40		,%%edi\n\t"\
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
			"movaps	%%xmm2		,    (%%esi)	/* <- t08 */\n\t"\
			"movaps	%%xmm3		,0x10(%%esi)	/* <- t09 */\n\t"\
			"movaps	%%xmm0		,    (%%edi)	/* <- t0a */\n\t"\
			"movaps	%%xmm1		,0x10(%%edi)	/* <- t0b */\n\t"\
		"addl	$0x40		,%%esi\n\t"\
			"subl	$0x080		,%%eax 		/* __i2 = s1p01r */\n\t"\
			"subl	$0x080		,%%ebx 		/* __i5 = s1p25r */\n\t"\
			"subl	$0x080		,%%ecx 		/* __i8 = s1p13r */\n\t"\
			"movaps	    (%%ebx)	,%%xmm4\n\t"\
			"movaps	0x10(%%ebx)	,%%xmm5\n\t"\
			"movaps	    (%%eax)	,%%xmm0\n\t"\
			"movaps	0x10(%%eax)	,%%xmm1\n\t"\
			"movaps	    (%%ecx)	,%%xmm2\n\t"\
			"movaps	0x10(%%ecx)	,%%xmm3\n\t"\
			"subpd	%%xmm2		,%%xmm4\n\t"\
			"subpd	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm2\n\t"\
			"addpd	%%xmm3		,%%xmm3\n\t"\
			"addpd	%%xmm4		,%%xmm2\n\t"\
			"addpd	%%xmm5		,%%xmm3\n\t"\
			"addpd	%%xmm2		,%%xmm0\n\t"\
			"addpd	%%xmm3		,%%xmm1\n\t"\
			"movaps	%%xmm0		,    (%%esi)	/* <- t0c */\n\t"\
			"movaps	%%xmm1		,0x10(%%esi)	/* <- t0d */\n\t"\
		"addl	$0x40		,%%edi\n\t"\
		"addl	$0x40		,%%esi\n\t"\
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
			"movaps	%%xmm2		,    (%%edi)	/* <- t0e */\n\t"\
			"movaps	%%xmm3		,0x10(%%edi)	/* <- t0f */\n\t"\
			"movaps	%%xmm0		,    (%%esi)	/* <- t0g */\n\t"\
			"movaps	%%xmm1		,0x10(%%esi)	/* <- t0h */\n\t"\
			"/******************************************************************************/\n\t"\
			"/*...and now do three more radix-3 transforms, including the twiddle factors: */\n\t"\
			"/******************************************************************************/\n\t"\
			"movl	%[__r00]	,%%eax 		/* __r00 */\n\t"\
			"addl	$0x360		,%%eax 		/* r30 */\n\t"\
			"movl	%%eax		,%%ebx\n\t"\
			"movl	%%eax		,%%ecx\n\t"\
			"addl	$0x60		,%%ebx 		/* __r36 */\n\t"\
			"addl	$0xc0		,%%ecx 		/* __r3c */\n\t"\
			"movaps	    (%%ebx)	,%%xmm4\n\t"\
			"movaps	0x10(%%ebx)	,%%xmm5\n\t"\
			"movaps	    (%%ecx)	,%%xmm2\n\t"\
			"movaps	0x10(%%ecx)	,%%xmm3\n\t"\
			"movaps	    (%%eax)	,%%xmm0\n\t"\
			"movaps	0x10(%%eax)	,%%xmm1\n\t"\
			"subpd	%%xmm2		,%%xmm4\n\t"\
			"subpd	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm2\n\t"\
			"addpd	%%xmm3		,%%xmm3\n\t"\
			"addpd	%%xmm4		,%%xmm2\n\t"\
			"addpd	%%xmm5		,%%xmm3\n\t"\
			"addpd	%%xmm2		,%%xmm0\n\t"\
			"addpd	%%xmm3		,%%xmm1\n\t"\
			"movaps	%%xmm0		,    (%%eax)\n\t"\
			"movaps	%%xmm1		,0x10(%%eax)\n\t"\
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
			"movaps	%%xmm2		,    (%%ebx)\n\t"\
			"movaps	%%xmm3		,0x10(%%ebx)\n\t"\
			"movaps	%%xmm0		,    (%%ecx)\n\t"\
			"movaps	%%xmm1		,0x10(%%ecx)\n\t"\
		"subl	$0x40		,%%edx 		/* c1 */\n\t"\
			"addl	$0x20		,%%eax 		/* __r32 */\n\t"\
			"addl	$0x20		,%%ebx 		/* __r38 */\n\t"\
			"addl	$0x20		,%%ecx 		/* __r3e */\n\t"\
			"movaps	    (%%ebx)	,%%xmm2\n\t"\
			"movaps	0x10(%%ebx)	,%%xmm3\n\t"\
			"movaps	    (%%edx)	,%%xmm6\n\t"\
			"movaps	0x10(%%edx)	,%%xmm7\n\t"\
		"addl	$0x20		,%%edx 		/* c2 */\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"subpd	%%xmm1		,%%xmm2\n\t"\
			"addpd	%%xmm0		,%%xmm3\n\t"\
			"movaps	    (%%ecx)	,%%xmm4\n\t"\
			"movaps	0x10(%%ecx)	,%%xmm5\n\t"\
			"movaps	    (%%edx)	,%%xmm6\n\t"\
			"movaps	0x10(%%edx)	,%%xmm7\n\t"\
		"addl	$0x20		,%%edx 		/* c3m1 */\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"subpd	%%xmm1		,%%xmm4\n\t"\
			"addpd	%%xmm0		,%%xmm5\n\t"\
			"movaps	    (%%edx)	,%%xmm6\n\t"\
			"movaps	0x10(%%edx)	,%%xmm7\n\t"\
			"movaps	    (%%eax)	,%%xmm0\n\t"\
			"movaps	0x10(%%eax)	,%%xmm1\n\t"\
			"subpd	%%xmm4		,%%xmm2\n\t"\
			"subpd	%%xmm5		,%%xmm3\n\t"\
			"addpd	%%xmm4		,%%xmm4\n\t"\
			"addpd	%%xmm5		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm4\n\t"\
			"addpd	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm4		,%%xmm0\n\t"\
			"addpd	%%xmm5		,%%xmm1\n\t"\
			"movaps	%%xmm0		,    (%%eax)\n\t"\
			"movaps	%%xmm1		,0x10(%%eax)\n\t"\
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
			"movaps	%%xmm4		,    (%%ebx)\n\t"\
			"movaps	%%xmm5		,0x10(%%ebx)\n\t"\
			"movaps	%%xmm0		,    (%%ecx)\n\t"\
			"movaps	%%xmm1		,0x10(%%ecx)\n\t"\
		"subl	$0x20		,%%edx 		/* c2 */\n\t"\
			"addl	$0x20		,%%eax 		/* __r34 */\n\t"\
			"addl	$0x20		,%%ebx 		/* __r3a */\n\t"\
			"addl	$0x20		,%%ecx 		/* __r3g */\n\t"\
			"movaps	    (%%ebx)	,%%xmm2\n\t"\
			"movaps	0x10(%%ebx)	,%%xmm3\n\t"\
			"movaps	    (%%edx)	,%%xmm6\n\t"\
			"movaps	0x10(%%edx)	,%%xmm7\n\t"\
		"addl	$0x40		,%%edx 		/* c4 */\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"subpd	%%xmm1		,%%xmm2\n\t"\
			"addpd	%%xmm0		,%%xmm3\n\t"\
			"movaps	    (%%ecx)	,%%xmm4\n\t"\
			"movaps	0x10(%%ecx)	,%%xmm5\n\t"\
			"movaps	    (%%edx)	,%%xmm6\n\t"\
			"movaps	0x10(%%edx)	,%%xmm7\n\t"\
		"subl	$0x20		,%%edx 		/* c3m1 */\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"subpd	%%xmm1		,%%xmm4\n\t"\
			"addpd	%%xmm0		,%%xmm5\n\t"\
			"movaps	    (%%edx)	,%%xmm6\n\t"\
			"movaps	0x10(%%edx)	,%%xmm7\n\t"\
			"movaps	    (%%eax)	,%%xmm0\n\t"\
			"movaps	0x10(%%eax)	,%%xmm1\n\t"\
			"\n\t"\
			"subpd	%%xmm4		,%%xmm2\n\t"\
			"subpd	%%xmm5		,%%xmm3\n\t"\
			"addpd	%%xmm4		,%%xmm4\n\t"\
			"addpd	%%xmm5		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm4\n\t"\
			"addpd	%%xmm3		,%%xmm5\n\t"\
			"addpd	%%xmm4		,%%xmm0\n\t"\
			"addpd	%%xmm5		,%%xmm1\n\t"\
			"movaps	%%xmm0		,    (%%eax)\n\t"\
			"movaps	%%xmm1		,0x10(%%eax)\n\t"\
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
			"movaps	%%xmm4		,    (%%ebx)\n\t"\
			"movaps	%%xmm5		,0x10(%%ebx)\n\t"\
			"movaps	%%xmm0		,    (%%ecx)\n\t"\
			"movaps	%%xmm1		,0x10(%%ecx)\n\t"\
			"\n\t"\
		"/**********************************/"\
		"/*** And now do 9 radix-4 DFTs: ***/"\
		"/**********************************/"\
			"\n\t"\
			"/* add0,1,3,2 = &a[j1+p00]+p0,1,2,3 */		\n\t"\
			"movl	%[__add],%%eax	/* Must use eax as base address throughout, since that is preserved in SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_B */\n\t"\
			"movl	%[__add],%%edx						\n\t"\
			"movl	%[__p01],%%esi	/* esi will store power-of-2 multiples of p01 throughout */\n\t"\
			"shll	$3,%%esi		/* Pointer offset for floating doubles */\n\t"\
			"addl	%%esi,%%edx							\n\t"\
			"movl	%%edx,%%ebx		/* add0+p01 */		\n\t"\
			"addl	%%esi,%%edx							\n\t"\
			"movl	%%edx,%%ecx		/* add0+p02 */		\n\t"\
			"addl	%%esi,%%edx		/* add0+p03 */		\n\t"\
			"/* SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_C(r00, 0x120, 0x240, eax,ebx,edx,ecx) */\n\t"\
			"movl	%[__r00],%%edi 	/* edi = __r00 */	\n\t"\
			"movl	$0x120	,%%esi 						\n\t"\
			"addl	%%edi	,%%esi 	/* esi = __r00 + 0x120 */\n\t"\
			"movaps	    (%%edi),%%xmm4					\n\t"\
			"movaps	    (%%esi),%%xmm6					\n\t"\
			"movaps	0x10(%%edi),%%xmm5					\n\t"\
			"movaps	0x10(%%esi),%%xmm7					\n\t"\
			"addl	$0x240,%%edi 						\n\t"\
			"addl	$0x240,%%esi 						\n\t"\
			"movaps	    (%%edi),%%xmm0					\n\t"\
			"movaps	    (%%esi),%%xmm2					\n\t"\
			"movaps	0x10(%%edi),%%xmm1					\n\t"\
			"movaps	0x10(%%esi),%%xmm3					\n\t"\
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
			"movaps	%%xmm0,    (%%ebx)					\n\t"\
			"movaps	%%xmm4,    (%%edx)					\n\t"\
			"movaps	%%xmm1,0x10(%%ebx)					\n\t"\
			"movaps	%%xmm5,0x10(%%ecx)					\n\t"\
			"addpd	%%xmm2,%%xmm2						\n\t"\
			"addpd	%%xmm7,%%xmm7						\n\t"\
			"addpd	%%xmm3,%%xmm3						\n\t"\
			"addpd	%%xmm6,%%xmm6						\n\t"\
			"addpd	%%xmm0,%%xmm2						\n\t"\
			"addpd	%%xmm4,%%xmm7						\n\t"\
			"addpd	%%xmm1,%%xmm3						\n\t"\
			"addpd	%%xmm5,%%xmm6						\n\t"\
			"movaps	%%xmm2,    (%%eax)					\n\t"\
			"movaps	%%xmm7,    (%%ecx)					\n\t"\
			"movaps	%%xmm3,0x10(%%eax)					\n\t"\
			"movaps	%%xmm6,0x10(%%edx)					\n\t"\
			"\n\t"\
			"/* add0,1,3,2 = &a[j1+p32]+p0,1,2,3 */		\n\t"\
			"movl	%[__p32],%%esi						\n\t"\
			"shll	$3,%%esi							\n\t"\
			"addl	%%esi,%%eax		/* &a[j1+p32] */	\n\t"\
			"addl	%%esi,%%ebx							\n\t"\
			"addl	%%esi,%%ecx							\n\t"\
			"addl	%%esi,%%edx							\n\t"\
			"/* SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_C(r02, 0x120, 0x240, eax,ebx,edx,ecx) */\n\t"\
			"subl	$0x220	,%%edi 	/* edi = __r02 */	\n\t"\
			"movl	$0x120	,%%esi 						\n\t"\
			"addl	%%edi	,%%esi 	/* esi = __r02 + 0x120 */\n\t"\
			"movaps	    (%%edi),%%xmm4					\n\t"\
			"movaps	    (%%esi),%%xmm6					\n\t"\
			"movaps	0x10(%%edi),%%xmm5					\n\t"\
			"movaps	0x10(%%esi),%%xmm7					\n\t"\
			"addl	$0x240,%%edi 						\n\t"\
			"addl	$0x240,%%esi 						\n\t"\
			"movaps	    (%%edi),%%xmm0					\n\t"\
			"movaps	    (%%esi),%%xmm2					\n\t"\
			"movaps	0x10(%%edi),%%xmm1					\n\t"\
			"movaps	0x10(%%esi),%%xmm3					\n\t"\
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
			"movaps	%%xmm0,    (%%ebx)					\n\t"\
			"movaps	%%xmm4,    (%%edx)					\n\t"\
			"movaps	%%xmm1,0x10(%%ebx)					\n\t"\
			"movaps	%%xmm5,0x10(%%ecx)					\n\t"\
			"addpd	%%xmm2,%%xmm2						\n\t"\
			"addpd	%%xmm7,%%xmm7						\n\t"\
			"addpd	%%xmm3,%%xmm3						\n\t"\
			"addpd	%%xmm6,%%xmm6						\n\t"\
			"addpd	%%xmm0,%%xmm2						\n\t"\
			"addpd	%%xmm4,%%xmm7						\n\t"\
			"addpd	%%xmm1,%%xmm3						\n\t"\
			"addpd	%%xmm5,%%xmm6						\n\t"\
			"movaps	%%xmm2,    (%%eax)					\n\t"\
			"movaps	%%xmm7,    (%%ecx)					\n\t"\
			"movaps	%%xmm3,0x10(%%eax)					\n\t"\
			"movaps	%%xmm6,0x10(%%edx)					\n\t"\
			"\n\t"\
			"/* add2,3,0,1 = &a[j1+p20]+p0,1,2,3 */		\n\t"\
			"movl	%[__p12],%%esi						\n\t"\
			"shll	$3,%%esi							\n\t"\
			"subl	%%esi,%%eax		/* &a[j1+p20] */	\n\t"\
			"subl	%%esi,%%ebx							\n\t"\
			"subl	%%esi,%%ecx							\n\t"\
			"subl	%%esi,%%edx							\n\t"\
			"/* SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_C(r04, 0x120, 0x240, ecx,edx,eax,ebx) */\n\t"\
			"subl	$0x220	,%%edi 	/* edi = __r04 */	\n\t"\
			"movl	$0x120	,%%esi 						\n\t"\
			"addl	%%edi	,%%esi 	/* esi = __r04 + 0x120 */\n\t"\
			"movaps	    (%%edi),%%xmm4					\n\t"\
			"movaps	    (%%esi),%%xmm6					\n\t"\
			"movaps	0x10(%%edi),%%xmm5					\n\t"\
			"movaps	0x10(%%esi),%%xmm7					\n\t"\
			"addl	$0x240,%%edi 						\n\t"\
			"addl	$0x240,%%esi 						\n\t"\
			"movaps	    (%%edi),%%xmm0					\n\t"\
			"movaps	    (%%esi),%%xmm2					\n\t"\
			"movaps	0x10(%%edi),%%xmm1					\n\t"\
			"movaps	0x10(%%esi),%%xmm3					\n\t"\
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
			"movaps	%%xmm0,    (%%edx)					\n\t"\
			"movaps	%%xmm4,    (%%eax)					\n\t"\
			"movaps	%%xmm1,0x10(%%edx)					\n\t"\
			"movaps	%%xmm5,0x10(%%ebx)					\n\t"\
			"addpd	%%xmm2,%%xmm2						\n\t"\
			"addpd	%%xmm7,%%xmm7						\n\t"\
			"addpd	%%xmm3,%%xmm3						\n\t"\
			"addpd	%%xmm6,%%xmm6						\n\t"\
			"addpd	%%xmm0,%%xmm2						\n\t"\
			"addpd	%%xmm4,%%xmm7						\n\t"\
			"addpd	%%xmm1,%%xmm3						\n\t"\
			"addpd	%%xmm5,%%xmm6						\n\t"\
			"movaps	%%xmm2,    (%%ecx)					\n\t"\
			"movaps	%%xmm7,    (%%ebx)					\n\t"\
			"movaps	%%xmm3,0x10(%%ecx)					\n\t"\
			"movaps	%%xmm6,0x10(%%eax)					\n\t"\
			"\n\t"\
			"/* add1,0,2,3 = &a[j1+p08]+p0,1,2,3 */		\n\t"\
			"movl	%[__p12],%%esi						\n\t"\
			"shll	$3,%%esi							\n\t"\
			"subl	%%esi,%%eax		/* &a[j1+p08] */	\n\t"\
			"subl	%%esi,%%ebx							\n\t"\
			"subl	%%esi,%%ecx							\n\t"\
			"subl	%%esi,%%edx							\n\t"\
			"/* SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_C(r06, 0x120, 0x240, ebx,eax,ecx,edx) */\n\t"\
			"subl	$0x220	,%%edi 	/* edi = __r06 */	\n\t"\
			"movl	$0x120	,%%esi 						\n\t"\
			"addl	%%edi	,%%esi 	/* esi = __r06 + 0x120 */\n\t"\
			"movaps	    (%%edi),%%xmm4					\n\t"\
			"movaps	    (%%esi),%%xmm6					\n\t"\
			"movaps	0x10(%%edi),%%xmm5					\n\t"\
			"movaps	0x10(%%esi),%%xmm7					\n\t"\
			"addl	$0x240,%%edi 						\n\t"\
			"addl	$0x240,%%esi 						\n\t"\
			"movaps	    (%%edi),%%xmm0					\n\t"\
			"movaps	    (%%esi),%%xmm2					\n\t"\
			"movaps	0x10(%%edi),%%xmm1					\n\t"\
			"movaps	0x10(%%esi),%%xmm3					\n\t"\
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
			"movaps	%%xmm0,    (%%eax)					\n\t"\
			"movaps	%%xmm4,    (%%ecx)					\n\t"\
			"movaps	%%xmm1,0x10(%%eax)					\n\t"\
			"movaps	%%xmm5,0x10(%%edx)					\n\t"\
			"addpd	%%xmm2,%%xmm2						\n\t"\
			"addpd	%%xmm7,%%xmm7						\n\t"\
			"addpd	%%xmm3,%%xmm3						\n\t"\
			"addpd	%%xmm6,%%xmm6						\n\t"\
			"addpd	%%xmm0,%%xmm2						\n\t"\
			"addpd	%%xmm4,%%xmm7						\n\t"\
			"addpd	%%xmm1,%%xmm3						\n\t"\
			"addpd	%%xmm5,%%xmm6						\n\t"\
			"movaps	%%xmm2,    (%%ebx)					\n\t"\
			"movaps	%%xmm7,    (%%edx)					\n\t"\
			"movaps	%%xmm3,0x10(%%ebx)					\n\t"\
			"movaps	%%xmm6,0x10(%%ecx)					\n\t"\
			"\n\t"\
			"/* add3,2,1,0 = &a[j1+p28]+p0,1,2,3 */		\n\t"\
			"movl	%[__p20],%%esi						\n\t"\
			"shll	$3,%%esi							\n\t"\
			"addl	%%esi,%%eax		/* &a[j1+p28] */	\n\t"\
			"addl	%%esi,%%ebx							\n\t"\
			"addl	%%esi,%%ecx							\n\t"\
			"addl	%%esi,%%edx							\n\t"\
			"/* SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_C(r08, 0x120, 0x240, edx,ecx,ebx,eax) */\n\t"\
			"subl	$0x220	,%%edi 	/* edi = __r08 */	\n\t"\
			"movl	$0x120	,%%esi 						\n\t"\
			"addl	%%edi	,%%esi 	/* esi = __r08 + 0x120 */\n\t"\
			"movaps	    (%%edi),%%xmm4					\n\t"\
			"movaps	    (%%esi),%%xmm6					\n\t"\
			"movaps	0x10(%%edi),%%xmm5					\n\t"\
			"movaps	0x10(%%esi),%%xmm7					\n\t"\
			"addl	$0x240,%%edi 						\n\t"\
			"addl	$0x240,%%esi 						\n\t"\
			"movaps	    (%%edi),%%xmm0					\n\t"\
			"movaps	    (%%esi),%%xmm2					\n\t"\
			"movaps	0x10(%%edi),%%xmm1					\n\t"\
			"movaps	0x10(%%esi),%%xmm3					\n\t"\
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
			"movaps	%%xmm0,    (%%ecx)					\n\t"\
			"movaps	%%xmm4,    (%%ebx)					\n\t"\
			"movaps	%%xmm1,0x10(%%ecx)					\n\t"\
			"movaps	%%xmm5,0x10(%%eax)					\n\t"\
			"addpd	%%xmm2,%%xmm2						\n\t"\
			"addpd	%%xmm7,%%xmm7						\n\t"\
			"addpd	%%xmm3,%%xmm3						\n\t"\
			"addpd	%%xmm6,%%xmm6						\n\t"\
			"addpd	%%xmm0,%%xmm2						\n\t"\
			"addpd	%%xmm4,%%xmm7						\n\t"\
			"addpd	%%xmm1,%%xmm3						\n\t"\
			"addpd	%%xmm5,%%xmm6						\n\t"\
			"movaps	%%xmm2,    (%%edx)					\n\t"\
			"movaps	%%xmm7,    (%%eax)					\n\t"\
			"movaps	%%xmm3,0x10(%%edx)					\n\t"\
			"movaps	%%xmm6,0x10(%%ebx)					\n\t"\
			"\n\t"\
			"/* add0,1,3,2 = &a[j1+p16]+p0,1,2,3 */		\n\t"\
			"movl	%[__p12],%%esi						\n\t"\
			"shll	$3,%%esi							\n\t"\
			"subl	%%esi,%%eax		/* &a[j1+p16] */	\n\t"\
			"subl	%%esi,%%ebx							\n\t"\
			"subl	%%esi,%%ecx							\n\t"\
			"subl	%%esi,%%edx							\n\t"\
			"/* SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_C(r0a, 0x120, 0x240, eax,ebx,edx,ecx) */\n\t"\
			"subl	$0x220	,%%edi 	/* edi = __r0a */	\n\t"\
			"movl	$0x120	,%%esi 						\n\t"\
			"addl	%%edi	,%%esi 	/* esi = __r0a + 0x120 */\n\t"\
			"movaps	    (%%edi),%%xmm4					\n\t"\
			"movaps	    (%%esi),%%xmm6					\n\t"\
			"movaps	0x10(%%edi),%%xmm5					\n\t"\
			"movaps	0x10(%%esi),%%xmm7					\n\t"\
			"addl	$0x240,%%edi 						\n\t"\
			"addl	$0x240,%%esi 						\n\t"\
			"movaps	    (%%edi),%%xmm0					\n\t"\
			"movaps	    (%%esi),%%xmm2					\n\t"\
			"movaps	0x10(%%edi),%%xmm1					\n\t"\
			"movaps	0x10(%%esi),%%xmm3					\n\t"\
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
			"movaps	%%xmm0,    (%%ebx)					\n\t"\
			"movaps	%%xmm4,    (%%edx)					\n\t"\
			"movaps	%%xmm1,0x10(%%ebx)					\n\t"\
			"movaps	%%xmm5,0x10(%%ecx)					\n\t"\
			"addpd	%%xmm2,%%xmm2						\n\t"\
			"addpd	%%xmm7,%%xmm7						\n\t"\
			"addpd	%%xmm3,%%xmm3						\n\t"\
			"addpd	%%xmm6,%%xmm6						\n\t"\
			"addpd	%%xmm0,%%xmm2						\n\t"\
			"addpd	%%xmm4,%%xmm7						\n\t"\
			"addpd	%%xmm1,%%xmm3						\n\t"\
			"addpd	%%xmm5,%%xmm6						\n\t"\
			"movaps	%%xmm2,    (%%eax)					\n\t"\
			"movaps	%%xmm7,    (%%ecx)					\n\t"\
			"movaps	%%xmm3,0x10(%%eax)					\n\t"\
			"movaps	%%xmm6,0x10(%%edx)					\n\t"\
			"\n\t"\
			"/* add2,3,0,1 = &a[j1+p04]+p0,1,2,3 */		\n\t"\
			"movl	%[__p12],%%esi						\n\t"\
			"shll	$3,%%esi							\n\t"\
			"subl	%%esi,%%eax		/* &a[j1+p04] */	\n\t"\
			"subl	%%esi,%%ebx							\n\t"\
			"subl	%%esi,%%ecx							\n\t"\
			"subl	%%esi,%%edx							\n\t"\
			"/* SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_C(r0c, 0x120, 0x240, ecx,edx,eax,ebx) */\n\t"\
			"subl	$0x220	,%%edi 	/* edi = __r0c */	\n\t"\
			"movl	$0x120	,%%esi 						\n\t"\
			"addl	%%edi	,%%esi 	/* esi = __r0c + 0x120 */\n\t"\
			"movaps	    (%%edi),%%xmm4					\n\t"\
			"movaps	    (%%esi),%%xmm6					\n\t"\
			"movaps	0x10(%%edi),%%xmm5					\n\t"\
			"movaps	0x10(%%esi),%%xmm7					\n\t"\
			"addl	$0x240,%%edi 						\n\t"\
			"addl	$0x240,%%esi 						\n\t"\
			"movaps	    (%%edi),%%xmm0					\n\t"\
			"movaps	    (%%esi),%%xmm2					\n\t"\
			"movaps	0x10(%%edi),%%xmm1					\n\t"\
			"movaps	0x10(%%esi),%%xmm3					\n\t"\
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
			"movaps	%%xmm0,    (%%edx)					\n\t"\
			"movaps	%%xmm4,    (%%eax)					\n\t"\
			"movaps	%%xmm1,0x10(%%edx)					\n\t"\
			"movaps	%%xmm5,0x10(%%ebx)					\n\t"\
			"addpd	%%xmm2,%%xmm2						\n\t"\
			"addpd	%%xmm7,%%xmm7						\n\t"\
			"addpd	%%xmm3,%%xmm3						\n\t"\
			"addpd	%%xmm6,%%xmm6						\n\t"\
			"addpd	%%xmm0,%%xmm2						\n\t"\
			"addpd	%%xmm4,%%xmm7						\n\t"\
			"addpd	%%xmm1,%%xmm3						\n\t"\
			"addpd	%%xmm5,%%xmm6						\n\t"\
			"movaps	%%xmm2,    (%%ecx)					\n\t"\
			"movaps	%%xmm7,    (%%ebx)					\n\t"\
			"movaps	%%xmm3,0x10(%%ecx)					\n\t"\
			"movaps	%%xmm6,0x10(%%eax)					\n\t"\
			"\n\t"\
			"/* add1,0,2,3 = &a[j1+p24]+p0,1,2,3 */		\n\t"\
			"movl	%[__p20],%%esi						\n\t"\
			"shll	$3,%%esi							\n\t"\
			"addl	%%esi,%%eax		/* &a[j1+p24] */	\n\t"\
			"addl	%%esi,%%ebx							\n\t"\
			"addl	%%esi,%%ecx							\n\t"\
			"addl	%%esi,%%edx							\n\t"\
			"/* SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_C(r0e, 0x120, 0x240, ebx,eax,ecx,edx) */\n\t"\
			"subl	$0x220	,%%edi 	/* edi = __r0e */	\n\t"\
			"movl	$0x120	,%%esi 						\n\t"\
			"addl	%%edi	,%%esi 	/* esi = __r0e + 0x120 */\n\t"\
			"movaps	    (%%edi),%%xmm4					\n\t"\
			"movaps	    (%%esi),%%xmm6					\n\t"\
			"movaps	0x10(%%edi),%%xmm5					\n\t"\
			"movaps	0x10(%%esi),%%xmm7					\n\t"\
			"addl	$0x240,%%edi 						\n\t"\
			"addl	$0x240,%%esi 						\n\t"\
			"movaps	    (%%edi),%%xmm0					\n\t"\
			"movaps	    (%%esi),%%xmm2					\n\t"\
			"movaps	0x10(%%edi),%%xmm1					\n\t"\
			"movaps	0x10(%%esi),%%xmm3					\n\t"\
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
			"movaps	%%xmm0,    (%%eax)					\n\t"\
			"movaps	%%xmm4,    (%%ecx)					\n\t"\
			"movaps	%%xmm1,0x10(%%eax)					\n\t"\
			"movaps	%%xmm5,0x10(%%edx)					\n\t"\
			"addpd	%%xmm2,%%xmm2						\n\t"\
			"addpd	%%xmm7,%%xmm7						\n\t"\
			"addpd	%%xmm3,%%xmm3						\n\t"\
			"addpd	%%xmm6,%%xmm6						\n\t"\
			"addpd	%%xmm0,%%xmm2						\n\t"\
			"addpd	%%xmm4,%%xmm7						\n\t"\
			"addpd	%%xmm1,%%xmm3						\n\t"\
			"addpd	%%xmm5,%%xmm6						\n\t"\
			"movaps	%%xmm2,    (%%ebx)					\n\t"\
			"movaps	%%xmm7,    (%%edx)					\n\t"\
			"movaps	%%xmm3,0x10(%%ebx)					\n\t"\
			"movaps	%%xmm6,0x10(%%ecx)					\n\t"\
			"\n\t"\
			"/* add,2,1,03 = &a[j1+p12]+p0,1,2,3 */		\n\t"\
			"movl	%[__p12],%%esi						\n\t"\
			"shll	$3,%%esi							\n\t"\
			"subl	%%esi,%%eax		/* &a[j1+p12] */	\n\t"\
			"subl	%%esi,%%ebx							\n\t"\
			"subl	%%esi,%%ecx							\n\t"\
			"subl	%%esi,%%edx							\n\t"\
			"/* SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_C(r0g, 0x120, 0x240, edx,ecx,ebx,eax) */\n\t"\
			"subl	$0x220	,%%edi 	/* edi = __r0g */	\n\t"\
			"movl	$0x120	,%%esi 						\n\t"\
			"addl	%%edi	,%%esi 	/* esi = __r0g + 0x120 */\n\t"\
			"movaps	    (%%edi),%%xmm4					\n\t"\
			"movaps	    (%%esi),%%xmm6					\n\t"\
			"movaps	0x10(%%edi),%%xmm5					\n\t"\
			"movaps	0x10(%%esi),%%xmm7					\n\t"\
			"addl	$0x240,%%edi 						\n\t"\
			"addl	$0x240,%%esi 						\n\t"\
			"movaps	    (%%edi),%%xmm0					\n\t"\
			"movaps	    (%%esi),%%xmm2					\n\t"\
			"movaps	0x10(%%edi),%%xmm1					\n\t"\
			"movaps	0x10(%%esi),%%xmm3					\n\t"\
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
			"movaps	%%xmm0,    (%%ecx)					\n\t"\
			"movaps	%%xmm4,    (%%ebx)					\n\t"\
			"movaps	%%xmm1,0x10(%%ecx)					\n\t"\
			"movaps	%%xmm5,0x10(%%eax)					\n\t"\
			"addpd	%%xmm2,%%xmm2						\n\t"\
			"addpd	%%xmm7,%%xmm7						\n\t"\
			"addpd	%%xmm3,%%xmm3						\n\t"\
			"addpd	%%xmm6,%%xmm6						\n\t"\
			"addpd	%%xmm0,%%xmm2						\n\t"\
			"addpd	%%xmm4,%%xmm7						\n\t"\
			"addpd	%%xmm1,%%xmm3						\n\t"\
			"addpd	%%xmm5,%%xmm6						\n\t"\
			"movaps	%%xmm2,    (%%edx)					\n\t"\
			"movaps	%%xmm7,    (%%eax)					\n\t"\
			"movaps	%%xmm3,0x10(%%edx)					\n\t"\
			"movaps	%%xmm6,0x10(%%ebx)					\n\t"\
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
			: "eax","ebx","ecx","edx","edi","esi"		/* Clobbered registers */\
		);\
	}

#endif	/* radix36_ditN_cy_dif1_gcc_h_included */

