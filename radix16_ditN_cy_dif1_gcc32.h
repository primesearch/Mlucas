/*******************************************************************************
*                                                                              *
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
#ifndef radix16_ditN_cy_dif1_gcc_h_included
#define radix16_ditN_cy_dif1_gcc_h_included

	#define SSE2_RADIX16_DIT_NOTWIDDLE(Xadd0,Xp1,Xp2,Xp3,Xp4,Xr1,Xr3,Xr5,Xr7,Xr9,Xr11,Xr17,Xr25,Xisrt2,Xcc0)\
	{\
	__asm__ volatile (\
	"pushl %%ebx	\n\t"/* Explicit save/restore of PIC register */\
		"movl	%[__add0],%%eax		\n\t"\
		"movl	%[__p1],%%ebx		\n\t"\
		"movl	%[__p2],%%ecx		\n\t"\
		"movl	%[__p3],%%edx		\n\t"\
		"movl	%[__p4],%%edi		\n\t"\
		"shll	$3,%%ebx			\n\t"\
		"shll	$3,%%ecx			\n\t"\
		"shll	$3,%%edx			\n\t"\
		"shll	$3,%%edi			\n\t"\
		"addl	%%eax,%%ebx			\n\t"\
		"addl	%%eax,%%ecx			\n\t"\
		"addl	%%eax,%%edx			\n\t"\
		"/* SSE2_RADIX4_DIT_0TWIDDLE_B(r1 ): */\n\t"\
		"movaps	    (%%eax),%%xmm2	\n\t"\
		"movaps	    (%%ecx),%%xmm6	\n\t"\
		"movaps	0x10(%%eax),%%xmm3	\n\t"\
		"movaps	0x10(%%ecx),%%xmm7	\n\t"\
		"movaps	    (%%ebx),%%xmm0	\n\t"\
		"movaps	    (%%edx),%%xmm4	\n\t"\
		"movaps	0x10(%%ebx),%%xmm1	\n\t"\
		"movaps	0x10(%%edx),%%xmm5	\n\t"\
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
		"movl	%[__r1], %%esi		\n\t"\
		"							\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm2		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"subpd	%%xmm6,%%xmm3		\n\t"\
		"movaps	%%xmm0,0x040(%%esi)	\n\t"\
		"movaps	%%xmm2,0x060(%%esi)	\n\t"\
		"movaps	%%xmm1,0x050(%%esi)	\n\t"\
		"movaps	%%xmm3,0x030(%%esi)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"addpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	%%xmm4,     (%%esi)	\n\t"\
		"movaps	%%xmm7,0x020(%%esi)	\n\t"\
		"movaps	%%xmm5,0x010(%%esi)	\n\t"\
		"movaps	%%xmm6,0x070(%%esi)	\n\t"\
		"addl	%%edi,%%eax			\n\t"\
		"addl	%%edi,%%ebx			\n\t"\
		"addl	%%edi,%%ecx			\n\t"\
		"addl	%%edi,%%edx			\n\t"\
		"/* SSE2_RADIX4_DIT_0TWIDDLE_B(r9 ): */\n\t"\
		"movaps	    (%%eax),%%xmm2	\n\t"\
		"movaps	    (%%ecx),%%xmm6	\n\t"\
		"movaps	0x10(%%eax),%%xmm3	\n\t"\
		"movaps	0x10(%%ecx),%%xmm7	\n\t"\
		"movaps	    (%%ebx),%%xmm0	\n\t"\
		"movaps	    (%%edx),%%xmm4	\n\t"\
		"movaps	0x10(%%ebx),%%xmm1	\n\t"\
		"movaps	0x10(%%edx),%%xmm5	\n\t"\
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
		"movl	%[__r9], %%esi		\n\t"\
		"							\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm2		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"subpd	%%xmm6,%%xmm3		\n\t"\
		"movaps	%%xmm0,0x040(%%esi)	\n\t"\
		"movaps	%%xmm2,0x060(%%esi)	\n\t"\
		"movaps	%%xmm1,0x050(%%esi)	\n\t"\
		"movaps	%%xmm3,0x030(%%esi)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"addpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	%%xmm4,     (%%esi)	\n\t"\
		"movaps	%%xmm7,0x020(%%esi)	\n\t"\
		"movaps	%%xmm5,0x010(%%esi)	\n\t"\
		"movaps	%%xmm6,0x070(%%esi)	\n\t"\
		"addl	%%edi,%%eax			\n\t"\
		"addl	%%edi,%%ebx			\n\t"\
		"addl	%%edi,%%ecx			\n\t"\
		"addl	%%edi,%%edx			\n\t"\
		"/* SSE2_RADIX4_DIT_0TWIDDLE_B(r17): */\n\t"\
		"movaps	    (%%eax),%%xmm2	\n\t"\
		"movaps	    (%%ecx),%%xmm6	\n\t"\
		"movaps	0x10(%%eax),%%xmm3	\n\t"\
		"movaps	0x10(%%ecx),%%xmm7	\n\t"\
		"movaps	    (%%ebx),%%xmm0	\n\t"\
		"movaps	    (%%edx),%%xmm4	\n\t"\
		"movaps	0x10(%%ebx),%%xmm1	\n\t"\
		"movaps	0x10(%%edx),%%xmm5	\n\t"\
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
		"movl	%[__r17], %%esi		\n\t"\
		"							\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm2		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"subpd	%%xmm6,%%xmm3		\n\t"\
		"movaps	%%xmm0,0x040(%%esi)	\n\t"\
		"movaps	%%xmm2,0x060(%%esi)	\n\t"\
		"movaps	%%xmm1,0x050(%%esi)	\n\t"\
		"movaps	%%xmm3,0x030(%%esi)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"addpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	%%xmm4,     (%%esi)	\n\t"\
		"movaps	%%xmm7,0x020(%%esi)	\n\t"\
		"movaps	%%xmm5,0x010(%%esi)	\n\t"\
		"movaps	%%xmm6,0x070(%%esi)	\n\t"\
		"addl	%%edi,%%eax			\n\t"\
		"addl	%%edi,%%ebx			\n\t"\
		"addl	%%edi,%%ecx			\n\t"\
		"addl	%%edi,%%edx			\n\t"\
		"/* SSE2_RADIX4_DIT_0TWIDDLE_B(r25): */\n\t"\
		"movaps	    (%%eax),%%xmm2	\n\t"\
		"movaps	    (%%ecx),%%xmm6	\n\t"\
		"movaps	0x10(%%eax),%%xmm3	\n\t"\
		"movaps	0x10(%%ecx),%%xmm7	\n\t"\
		"movaps	    (%%ebx),%%xmm0	\n\t"\
		"movaps	    (%%edx),%%xmm4	\n\t"\
		"movaps	0x10(%%ebx),%%xmm1	\n\t"\
		"movaps	0x10(%%edx),%%xmm5	\n\t"\
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
		"movl	%[__r25], %%esi		\n\t"\
		"							\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm2		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"subpd	%%xmm6,%%xmm3		\n\t"\
		"movaps	%%xmm0,0x040(%%esi)	\n\t"\
		"movaps	%%xmm2,0x060(%%esi)	\n\t"\
		"movaps	%%xmm1,0x050(%%esi)	\n\t"\
		"movaps	%%xmm3,0x030(%%esi)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"addpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	%%xmm4,     (%%esi)	\n\t"\
		"movaps	%%xmm7,0x020(%%esi)	\n\t"\
		"movaps	%%xmm5,0x010(%%esi)	\n\t"\
		"movaps	%%xmm6,0x070(%%esi)	\n\t"\
		"movl	%[__r1],%%eax		\n\t"\
		"							\n\t"\
		"movaps	     (%%eax),%%xmm0	\n\t"\
		"movaps	0x010(%%eax),%%xmm1	\n\t"\
		"movaps	0x080(%%eax),%%xmm2	\n\t"\
		"movaps	0x090(%%eax),%%xmm3	\n\t"\
		"							\n\t"\
		"subpd	0x080(%%eax),%%xmm0	\n\t"\
		"subpd	0x090(%%eax),%%xmm1	\n\t"\
		"addpd	     (%%eax),%%xmm2	\n\t"\
		"addpd	0x010(%%eax),%%xmm3	\n\t"\
		"							\n\t"\
		"movaps	0x100(%%eax),%%xmm4	\n\t"\
		"movaps	0x110(%%eax),%%xmm5	\n\t"\
		"movaps	0x180(%%eax),%%xmm6	\n\t"\
		"movaps	0x190(%%eax),%%xmm7	\n\t"\
		"							\n\t"\
		"subpd	0x180(%%eax),%%xmm4	\n\t"\
		"subpd	0x190(%%eax),%%xmm5	\n\t"\
		"addpd	0x100(%%eax),%%xmm6	\n\t"\
		"addpd	0x110(%%eax),%%xmm7	\n\t"\
		"							\n\t"\
		"subpd	%%xmm6,%%xmm2		\n\t"\
		"subpd	%%xmm7,%%xmm3		\n\t"\
		"movaps	%%xmm2,0x100(%%eax)	\n\t"\
		"movaps	%%xmm3,0x110(%%eax)	\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm2,%%xmm6		\n\t"\
		"addpd	%%xmm3,%%xmm7		\n\t"\
		"movaps	%%xmm6,     (%%eax)	\n\t"\
		"movaps	%%xmm7,0x010(%%eax)	\n\t"\
		"							\n\t"\
		"subpd	%%xmm5,%%xmm0		\n\t"\
		"subpd	%%xmm4,%%xmm1		\n\t"\
		"movaps	%%xmm0,0x180(%%eax)	\n\t"\
		"movaps	%%xmm1,0x090(%%eax)	\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm0,%%xmm5		\n\t"\
		"addpd	%%xmm1,%%xmm4		\n\t"\
		"movaps	%%xmm5,0x080(%%eax)	\n\t"\
		"movaps	%%xmm4,0x190(%%eax)	\n\t"\
		"							\n\t"\
		"movl	%[__r5],%%eax		\n\t"\
		"movl	%[__isrt2],%%ebx	\n\t"\
		"movaps	(%%ebx),%%xmm2		\n\t"\
		"							\n\t"\
		"movaps	0x100(%%eax),%%xmm4	\n\t"\
		"movaps	0x110(%%eax),%%xmm5	\n\t"\
		"movaps	0x180(%%eax),%%xmm0	\n\t"\
		"movaps	0x190(%%eax),%%xmm1	\n\t"\
		"							\n\t"\
		"addpd	0x110(%%eax),%%xmm4	\n\t"\
		"subpd	0x100(%%eax),%%xmm5	\n\t"\
		"subpd	0x190(%%eax),%%xmm0	\n\t"\
		"addpd	0x180(%%eax),%%xmm1	\n\t"\
		"mulpd	%%xmm2,%%xmm4		\n\t"\
		"mulpd	%%xmm2,%%xmm5		\n\t"\
		"mulpd	%%xmm2,%%xmm0		\n\t"\
		"mulpd	%%xmm2,%%xmm1		\n\t"\
		"movaps	%%xmm4,%%xmm6		\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t"\
		"							\n\t"\
		"subpd	%%xmm0,%%xmm4		\n\t"\
		"subpd	%%xmm1,%%xmm5		\n\t"\
		"addpd	%%xmm0,%%xmm6		\n\t"\
		"addpd	%%xmm1,%%xmm7		\n\t"\
		"							\n\t"\
		"movaps	     (%%eax),%%xmm0	\n\t"\
		"movaps	0x010(%%eax),%%xmm1	\n\t"\
		"movaps	0x080(%%eax),%%xmm2	\n\t"\
		"movaps	0x090(%%eax),%%xmm3	\n\t"\
		"							\n\t"\
		"subpd	0x090(%%eax),%%xmm0	\n\t"\
		"subpd	0x080(%%eax),%%xmm1	\n\t"\
		"addpd	     (%%eax),%%xmm3	\n\t"\
		"addpd	0x010(%%eax),%%xmm2	\n\t"\
		"							\n\t"\
		"subpd	%%xmm4,%%xmm3		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"movaps	%%xmm3,0x100(%%eax)	\n\t"\
		"movaps	%%xmm1,0x110(%%eax)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm3,%%xmm4		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"movaps	%%xmm4,     (%%eax)	\n\t"\
		"movaps	%%xmm5,0x010(%%eax)	\n\t"\
		"							\n\t"\
		"subpd	%%xmm7,%%xmm0		\n\t"\
		"subpd	%%xmm6,%%xmm2		\n\t"\
		"movaps	%%xmm0,0x180(%%eax)	\n\t"\
		"movaps	%%xmm2,0x090(%%eax)	\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm7		\n\t"\
		"addpd	%%xmm2,%%xmm6		\n\t"\
		"movaps	%%xmm7,0x080(%%eax)	\n\t"\
		"movaps	%%xmm6,0x190(%%eax)	\n\t"\
		"							\n\t"\
		"movl	%[__r3],%%eax		\n\t"\
		"movl	%[__isrt2],%%ebx	\n\t"\
		"movl	%[__cc0],%%ecx		\n\t"\
		"							\n\t"\
		"movaps	0x180(%%eax),%%xmm0	\n\t"\
		"movaps	0x190(%%eax),%%xmm1	\n\t"\
		"movaps	0x180(%%eax),%%xmm2	\n\t"\
		"movaps	0x190(%%eax),%%xmm3	\n\t"\
		"							\n\t"\
		"mulpd	0x010(%%ecx),%%xmm0	\n\t"\
		"mulpd	0x010(%%ecx),%%xmm1	\n\t"\
		"mulpd	     (%%ecx),%%xmm2	\n\t"\
		"mulpd	     (%%ecx),%%xmm3	\n\t"\
		"subpd	%%xmm2,%%xmm1		\n\t"\
		"addpd	%%xmm3,%%xmm0		\n\t"\
		"							\n\t"\
		"movaps	0x100(%%eax),%%xmm4	\n\t"\
		"movaps	0x110(%%eax),%%xmm5	\n\t"\
		"movaps	0x100(%%eax),%%xmm6	\n\t"\
		"movaps	0x110(%%eax),%%xmm7	\n\t"\
		"							\n\t"\
		"mulpd	     (%%ecx),%%xmm4	\n\t"\
		"mulpd	     (%%ecx),%%xmm5	\n\t"\
		"mulpd	0x010(%%ecx),%%xmm6	\n\t"\
		"mulpd	0x010(%%ecx),%%xmm7	\n\t"\
		"subpd	%%xmm6,%%xmm5		\n\t"\
		"addpd	%%xmm7,%%xmm4		\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t"\
		"movaps	%%xmm4,%%xmm6		\n\t"\
		"							\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"subpd	%%xmm0,%%xmm6		\n\t"\
		"subpd	%%xmm1,%%xmm7		\n\t"\
		"							\n\t"\
		"movaps	0x080(%%eax),%%xmm2	\n\t"\
		"movaps	0x090(%%eax),%%xmm3	\n\t"\
		"movaps	     (%%eax),%%xmm0	\n\t"\
		"movaps	0x010(%%eax),%%xmm1	\n\t"\
		"addpd	0x090(%%eax),%%xmm2	\n\t"\
		"subpd	0x080(%%eax),%%xmm3	\n\t"\
		"mulpd	     (%%ebx),%%xmm2	\n\t"\
		"mulpd	     (%%ebx),%%xmm3	\n\t"\
		"							\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"							\n\t"\
		"subpd	%%xmm4,%%xmm2		\n\t"\
		"subpd	%%xmm5,%%xmm3		\n\t"\
		"movaps	%%xmm2,0x100(%%eax)	\n\t"\
		"movaps	%%xmm3,0x110(%%eax)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm2,%%xmm4		\n\t"\
		"addpd	%%xmm3,%%xmm5		\n\t"\
		"movaps	%%xmm4,     (%%eax)	\n\t"\
		"movaps	%%xmm5,0x010(%%eax)	\n\t"\
		"							\n\t"\
		"subpd	%%xmm7,%%xmm0		\n\t"\
		"subpd	%%xmm6,%%xmm1		\n\t"\
		"movaps	%%xmm0,0x180(%%eax)	\n\t"\
		"movaps	%%xmm1,0x090(%%eax)	\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm7		\n\t"\
		"addpd	%%xmm1,%%xmm6		\n\t"\
		"movaps	%%xmm7,0x080(%%eax)	\n\t"\
		"movaps	%%xmm6,0x190(%%eax)	\n\t"\
		"							\n\t"\
		"movl	%[__r7],%%eax		\n\t"\
		"movl	%[__isrt2],%%ebx	\n\t"\
		"movl	%[__cc0],%%ecx		\n\t"\
		"							\n\t"\
		"movaps	0x180(%%eax),%%xmm0	\n\t"\
		"movaps	0x190(%%eax),%%xmm1	\n\t"\
		"movaps	0x180(%%eax),%%xmm2	\n\t"\
		"movaps	0x190(%%eax),%%xmm3	\n\t"\
		"							\n\t"\
		"mulpd	     (%%ecx),%%xmm0	\n\t"\
		"mulpd	     (%%ecx),%%xmm1	\n\t"\
		"mulpd	0x010(%%ecx),%%xmm2	\n\t"\
		"mulpd	0x010(%%ecx),%%xmm3	\n\t"\
		"subpd	%%xmm2,%%xmm1		\n\t"\
		"addpd	%%xmm3,%%xmm0		\n\t"\
		"							\n\t"\
		"movaps	0x100(%%eax),%%xmm4	\n\t"\
		"movaps	0x110(%%eax),%%xmm5	\n\t"\
		"movaps	0x100(%%eax),%%xmm6	\n\t"\
		"movaps	0x110(%%eax),%%xmm7	\n\t"\
		"							\n\t"\
		"mulpd	0x010(%%ecx),%%xmm4	\n\t"\
		"mulpd	0x010(%%ecx),%%xmm5	\n\t"\
		"mulpd	     (%%ecx),%%xmm6	\n\t"\
		"mulpd	     (%%ecx),%%xmm7	\n\t"\
		"subpd	%%xmm6,%%xmm5		\n\t"\
		"addpd	%%xmm7,%%xmm4		\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t"\
		"movaps	%%xmm4,%%xmm6		\n\t"\
		"							\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"subpd	%%xmm0,%%xmm6		\n\t"\
		"subpd	%%xmm1,%%xmm7		\n\t"\
		"							\n\t"\
		"movaps	0x080(%%eax),%%xmm2	\n\t"\
		"movaps	0x090(%%eax),%%xmm3	\n\t"\
		"movaps	     (%%eax),%%xmm0	\n\t"\
		"movaps	0x010(%%eax),%%xmm1	\n\t"\
		"subpd	0x090(%%eax),%%xmm2	\n\t"\
		"addpd	0x080(%%eax),%%xmm3	\n\t"\
		"mulpd	     (%%ebx),%%xmm2	\n\t"\
		"mulpd	     (%%ebx),%%xmm3	\n\t"\
		"							\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"							\n\t"\
		"subpd	%%xmm6,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm1		\n\t"\
		"movaps	%%xmm0,0x100(%%eax)	\n\t"\
		"movaps	%%xmm1,0x110(%%eax)	\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm0,%%xmm6		\n\t"\
		"addpd	%%xmm1,%%xmm7		\n\t"\
		"movaps	%%xmm6,     (%%eax)	\n\t"\
		"movaps	%%xmm7,0x010(%%eax)	\n\t"\
		"							\n\t"\
		"subpd	%%xmm5,%%xmm2		\n\t"\
		"subpd	%%xmm4,%%xmm3		\n\t"\
		"movaps	%%xmm2,0x180(%%eax)	\n\t"\
		"movaps	%%xmm3,0x090(%%eax)	\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm2,%%xmm5		\n\t"\
		"addpd	%%xmm3,%%xmm4		\n\t"\
		"movaps	%%xmm5,0x080(%%eax)	\n\t"\
		"movaps	%%xmm4,0x190(%%eax)	\n\t"\
	"popl %%ebx	\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__p1] "m" (Xp1)\
		 ,[__p2] "m" (Xp2)\
		 ,[__p3] "m" (Xp3)\
		 ,[__p4] "m" (Xp4)\
		 ,[__r1] "m" (Xr1)\
		 ,[__r3] "m" (Xr3)\
		 ,[__r5] "m" (Xr5)\
		 ,[__r7] "m" (Xr7)\
		 ,[__r9] "m" (Xr9)\
		 ,[__r17] "m" (Xr17)\
		 ,[__r25] "m" (Xr25)\
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__cc0] "m" (Xcc0)\
		: "cc","memory","eax",/*"ebx",*/"ecx","edx","edi","esi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"	/* Clobbered registers */\
	);\
	}


	#define SSE2_RADIX16_DIF_NOTWIDDLE(Xadd0,Xp1,Xp2,Xp3,Xp4,Xr1,Xr3,Xr5,Xr7,Xr9,Xr11,Xr13,Xr15,Xr17,Xr19,Xr21,Xr23,Xr25,Xr27,Xr29,Xr31,Xisrt2,Xcc0)\
	{\
	__asm__ volatile (\
	"pushl %%ebx	\n\t"/* Explicit save/restore of PIC register */\
		"/* SSE2_RADIX4_DIF_IN_PLACE(r1 , r17, r9 , r25): */\n\t"\
		"movl	%[__r1] ,%%eax		\n\t"\
		"movl	%[__r17],%%ebx		\n\t"\
		"movl	%[__r9] ,%%ecx		\n\t"\
		"movl	%[__r25],%%edx		\n\t"\
		"							\n\t"\
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"movaps	    (%%eax),%%xmm2	\n\t"\
		"movaps	0x10(%%eax),%%xmm3	\n\t"\
		"							\n\t"\
		"addpd	    (%%ebx),%%xmm0	\n\t"\
		"addpd	0x10(%%ebx),%%xmm1	\n\t"\
		"subpd	    (%%ebx),%%xmm2	\n\t"\
		"subpd	0x10(%%ebx),%%xmm3	\n\t"\
		"							\n\t"\
		"movaps	    (%%ecx),%%xmm4	\n\t"\
		"movaps	0x10(%%ecx),%%xmm5	\n\t"\
		"movaps	    (%%ecx),%%xmm6	\n\t"\
		"movaps	0x10(%%ecx),%%xmm7	\n\t"\
		"							\n\t"\
		"addpd	    (%%edx),%%xmm4	\n\t"\
		"addpd	0x10(%%edx),%%xmm5	\n\t"\
		"subpd	    (%%edx),%%xmm6	\n\t"\
		"subpd	0x10(%%edx),%%xmm7	\n\t"\
		"							\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"movaps	%%xmm0,     (%%ebx)	\n\t"\
		"movaps	%%xmm1,0x010(%%ebx)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"movaps	%%xmm4,     (%%eax)	\n\t"\
		"movaps	%%xmm5,0x010(%%eax)	\n\t"\
		"							\n\t"\
		"subpd	%%xmm7,%%xmm2		\n\t"\
		"subpd	%%xmm6,%%xmm3		\n\t"\
		"movaps	%%xmm2,     (%%ecx)	\n\t"\
		"movaps	%%xmm3,0x010(%%edx)	\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	%%xmm7,     (%%edx)	\n\t"\
		"movaps	%%xmm6,0x010(%%ecx)	\n\t"\
		"/* SSE2_RADIX4_DIF_IN_PLACE(r5 , r21, r13, r29): */\n\t"\
		"movl	%[__r5] ,%%eax		\n\t"\
		"movl	%[__r21],%%ebx		\n\t"\
		"movl	%[__r13],%%ecx		\n\t"\
		"movl	%[__r29],%%edx		\n\t"\
		"							\n\t"\
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"movaps	    (%%eax),%%xmm2	\n\t"\
		"movaps	0x10(%%eax),%%xmm3	\n\t"\
		"							\n\t"\
		"addpd	    (%%ebx),%%xmm0	\n\t"\
		"addpd	0x10(%%ebx),%%xmm1	\n\t"\
		"subpd	    (%%ebx),%%xmm2	\n\t"\
		"subpd	0x10(%%ebx),%%xmm3	\n\t"\
		"							\n\t"\
		"movaps	    (%%ecx),%%xmm4	\n\t"\
		"movaps	0x10(%%ecx),%%xmm5	\n\t"\
		"movaps	    (%%ecx),%%xmm6	\n\t"\
		"movaps	0x10(%%ecx),%%xmm7	\n\t"\
		"							\n\t"\
		"addpd	    (%%edx),%%xmm4	\n\t"\
		"addpd	0x10(%%edx),%%xmm5	\n\t"\
		"subpd	    (%%edx),%%xmm6	\n\t"\
		"subpd	0x10(%%edx),%%xmm7	\n\t"\
		"							\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"movaps	%%xmm0,     (%%ebx)	\n\t"\
		"movaps	%%xmm1,0x010(%%ebx)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"movaps	%%xmm4,     (%%eax)	\n\t"\
		"movaps	%%xmm5,0x010(%%eax)	\n\t"\
		"							\n\t"\
		"subpd	%%xmm7,%%xmm2		\n\t"\
		"subpd	%%xmm6,%%xmm3		\n\t"\
		"movaps	%%xmm2,     (%%ecx)	\n\t"\
		"movaps	%%xmm3,0x010(%%edx)	\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	%%xmm7,     (%%edx)	\n\t"\
		"movaps	%%xmm6,0x010(%%ecx)	\n\t"\
		"/* SSE2_RADIX4_DIF_IN_PLACE(r3 , r19, r11, r27): */\n\t"\
		"movl	%[__r3] ,%%eax		\n\t"\
		"movl	%[__r19],%%ebx		\n\t"\
		"movl	%[__r11],%%ecx		\n\t"\
		"movl	%[__r27],%%edx		\n\t"\
		"							\n\t"\
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"movaps	    (%%eax),%%xmm2	\n\t"\
		"movaps	0x10(%%eax),%%xmm3	\n\t"\
		"							\n\t"\
		"addpd	    (%%ebx),%%xmm0	\n\t"\
		"addpd	0x10(%%ebx),%%xmm1	\n\t"\
		"subpd	    (%%ebx),%%xmm2	\n\t"\
		"subpd	0x10(%%ebx),%%xmm3	\n\t"\
		"							\n\t"\
		"movaps	    (%%ecx),%%xmm4	\n\t"\
		"movaps	0x10(%%ecx),%%xmm5	\n\t"\
		"movaps	    (%%ecx),%%xmm6	\n\t"\
		"movaps	0x10(%%ecx),%%xmm7	\n\t"\
		"							\n\t"\
		"addpd	    (%%edx),%%xmm4	\n\t"\
		"addpd	0x10(%%edx),%%xmm5	\n\t"\
		"subpd	    (%%edx),%%xmm6	\n\t"\
		"subpd	0x10(%%edx),%%xmm7	\n\t"\
		"							\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"movaps	%%xmm0,     (%%ebx)	\n\t"\
		"movaps	%%xmm1,0x010(%%ebx)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"movaps	%%xmm4,     (%%eax)	\n\t"\
		"movaps	%%xmm5,0x010(%%eax)	\n\t"\
		"							\n\t"\
		"subpd	%%xmm7,%%xmm2		\n\t"\
		"subpd	%%xmm6,%%xmm3		\n\t"\
		"movaps	%%xmm2,     (%%ecx)	\n\t"\
		"movaps	%%xmm3,0x010(%%edx)	\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	%%xmm7,     (%%edx)	\n\t"\
		"movaps	%%xmm6,0x010(%%ecx)	\n\t"\
		"/* SSE2_RADIX4_DIF_IN_PLACE(r7 , r23, r15, r31): */\n\t"\
		"movl	%[__r7] ,%%eax		\n\t"\
		"movl	%[__r23],%%ebx		\n\t"\
		"movl	%[__r15],%%ecx		\n\t"\
		"movl	%[__r31],%%edx		\n\t"\
		"							\n\t"\
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"movaps	    (%%eax),%%xmm2	\n\t"\
		"movaps	0x10(%%eax),%%xmm3	\n\t"\
		"							\n\t"\
		"addpd	    (%%ebx),%%xmm0	\n\t"\
		"addpd	0x10(%%ebx),%%xmm1	\n\t"\
		"subpd	    (%%ebx),%%xmm2	\n\t"\
		"subpd	0x10(%%ebx),%%xmm3	\n\t"\
		"							\n\t"\
		"movaps	    (%%ecx),%%xmm4	\n\t"\
		"movaps	0x10(%%ecx),%%xmm5	\n\t"\
		"movaps	    (%%ecx),%%xmm6	\n\t"\
		"movaps	0x10(%%ecx),%%xmm7	\n\t"\
		"							\n\t"\
		"addpd	    (%%edx),%%xmm4	\n\t"\
		"addpd	0x10(%%edx),%%xmm5	\n\t"\
		"subpd	    (%%edx),%%xmm6	\n\t"\
		"subpd	0x10(%%edx),%%xmm7	\n\t"\
		"							\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"movaps	%%xmm0,     (%%ebx)	\n\t"\
		"movaps	%%xmm1,0x010(%%ebx)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"movaps	%%xmm4,     (%%eax)	\n\t"\
		"movaps	%%xmm5,0x010(%%eax)	\n\t"\
		"							\n\t"\
		"subpd	%%xmm7,%%xmm2		\n\t"\
		"subpd	%%xmm6,%%xmm3		\n\t"\
		"movaps	%%xmm2,     (%%ecx)	\n\t"\
		"movaps	%%xmm3,0x010(%%edx)	\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	%%xmm7,     (%%edx)	\n\t"\
		"movaps	%%xmm6,0x010(%%ecx)	\n\t"\
		"							\n\t"\
		"movl	%[__add0],%%eax		\n\t"\
		"movl	%[__p1],%%ebx		\n\t"\
		"movl	%[__p2],%%ecx		\n\t"\
		"movl	%[__p3],%%edx		\n\t"\
		"movl	%[__p4],%%edi		\n\t"\
		"shll	$3,%%ebx			\n\t"\
		"shll	$3,%%ecx			\n\t"\
		"shll	$3,%%edx			\n\t"\
		"shll	$3,%%edi			\n\t"\
		"addl	%%eax,%%ebx			\n\t"\
		"addl	%%eax,%%ecx			\n\t"\
		"addl	%%eax,%%edx			\n\t"\
		"							\n\t"\
		"movl	%[__r1],%%esi		\n\t"\
		"							\n\t"\
		"movaps	     (%%esi),%%xmm0	\n\t"\
		"movaps	0x010(%%esi),%%xmm1	\n\t"\
		"movaps	0x040(%%esi),%%xmm2	\n\t"\
		"movaps	0x050(%%esi),%%xmm3	\n\t"\
		"							\n\t"\
		"subpd	0x040(%%esi),%%xmm0	\n\t"\
		"subpd	0x050(%%esi),%%xmm1	\n\t"\
		"addpd	     (%%esi),%%xmm2	\n\t"\
		"addpd	0x010(%%esi),%%xmm3	\n\t"\
		"							\n\t"\
		"movaps	0x020(%%esi),%%xmm4	\n\t"\
		"movaps	0x030(%%esi),%%xmm5	\n\t"\
		"movaps	0x060(%%esi),%%xmm6	\n\t"\
		"movaps	0x070(%%esi),%%xmm7	\n\t"\
		"							\n\t"\
		"subpd	0x060(%%esi),%%xmm4	\n\t"\
		"subpd	0x070(%%esi),%%xmm5	\n\t"\
		"addpd	0x020(%%esi),%%xmm6	\n\t"\
		"addpd	0x030(%%esi),%%xmm7	\n\t"\
		"							\n\t"\
		"subpd	%%xmm6,%%xmm2		\n\t"\
		"subpd	%%xmm7,%%xmm3		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"movaps	%%xmm2,    (%%ebx)	\n\t"\
		"movaps	%%xmm3,0x10(%%ebx)	\n\t"\
		"addpd	%%xmm2,		%%xmm6	\n\t"\
		"addpd	%%xmm3,		%%xmm7	\n\t"\
		"movaps	%%xmm6,    (%%eax)	\n\t"\
		"movaps	%%xmm7,0x10(%%eax)	\n\t"\
		"							\n\t"\
		"subpd	%%xmm5,		%%xmm0	\n\t"\
		"subpd	%%xmm4,		%%xmm1	\n\t"\
		"addpd	%%xmm5,		%%xmm5	\n\t"\
		"addpd	%%xmm4,		%%xmm4	\n\t"\
		"movaps	%%xmm0,    (%%ecx)	\n\t"\
		"movaps	%%xmm1,0x10(%%edx)	\n\t"\
		"addpd	%%xmm0,		%%xmm5	\n\t"\
		"addpd	%%xmm1,		%%xmm4	\n\t"\
		"movaps	%%xmm5,    (%%edx)	\n\t"\
		"movaps	%%xmm4,0x10(%%ecx)	\n\t"\
		"							\n\t"\
		"addl	%%edi,%%eax			\n\t"\
		"addl	%%edi,%%ebx			\n\t"\
		"addl	%%edi,%%ecx			\n\t"\
		"addl	%%edi,%%edx			\n\t"\
		"							\n\t"\
		"movl	%[__r17],%%esi		\n\t"\
		"							\n\t"\
		"movaps	     (%%esi),%%xmm0	\n\t"\
		"movaps	0x010(%%esi),%%xmm1	\n\t"\
		"movaps	0x040(%%esi),%%xmm2	\n\t"\
		"movaps	0x050(%%esi),%%xmm3	\n\t"\
		"							\n\t"\
		"subpd	0x050(%%esi),%%xmm0	\n\t"\
		"subpd	0x040(%%esi),%%xmm1	\n\t"\
		"addpd	0x010(%%esi),%%xmm2	\n\t"\
		"addpd	     (%%esi),%%xmm3	\n\t"\
		"							\n\t"\
		"movaps	0x020(%%esi),%%xmm4	\n\t"\
		"movaps	0x030(%%esi),%%xmm5	\n\t"\
		"movaps	0x060(%%esi),%%xmm6	\n\t"\
		"movaps	0x070(%%esi),%%xmm7	\n\t"\
		"							\n\t"\
		"subpd	0x030(%%esi),%%xmm4	\n\t"\
		"addpd	0x020(%%esi),%%xmm5	\n\t"\
		"addpd	0x070(%%esi),%%xmm6	\n\t"\
		"subpd	0x060(%%esi),%%xmm7	\n\t"\
		"							\n\t"\
		"movl	%[__isrt2],%%esi	\n\t"\
		"mulpd	(%%esi),%%xmm4		\n\t"\
		"mulpd	(%%esi),%%xmm5		\n\t"\
		"mulpd	(%%esi),%%xmm6		\n\t"\
		"mulpd	(%%esi),%%xmm7		\n\t"\
		"							\n\t"\
		"subpd	%%xmm6,%%xmm4		\n\t"\
		"subpd	%%xmm7,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm6		\n\t"\
		"addpd	%%xmm5,%%xmm7		\n\t"\
		"							\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm2		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"							\n\t"\
		"movaps	%%xmm0,    (%%ebx)	\n\t"\
		"movaps	%%xmm2,0x10(%%ebx)	\n\t"\
		"addpd	%%xmm0,		%%xmm4	\n\t"\
		"addpd	%%xmm2,		%%xmm5	\n\t"\
		"movaps	%%xmm4,    (%%eax)	\n\t"\
		"movaps	%%xmm5,0x10(%%eax)	\n\t"\
		"							\n\t"\
		"subpd	%%xmm7,		%%xmm3	\n\t"\
		"subpd	%%xmm6,		%%xmm1	\n\t"\
		"addpd	%%xmm7,		%%xmm7	\n\t"\
		"addpd	%%xmm6,		%%xmm6	\n\t"\
		"movaps	%%xmm3,    (%%ecx)	\n\t"\
		"movaps	%%xmm1,0x10(%%edx)	\n\t"\
		"addpd	%%xmm3,		%%xmm7	\n\t"\
		"addpd	%%xmm1,		%%xmm6	\n\t"\
		"movaps	%%xmm7,    (%%edx)	\n\t"\
		"movaps	%%xmm6,0x10(%%ecx)	\n\t"\
		"							\n\t"\
		"addl	%%edi,%%eax			\n\t"\
		"addl	%%edi,%%ebx			\n\t"\
		"addl	%%edi,%%ecx			\n\t"\
		"addl	%%edi,%%edx			\n\t"\
		"							\n\t"\
		"movl	%[__r9] ,%%esi		\n\t"\
		"movl	%[__cc0],%%edi		\n\t"\
		"movaps	0x020(%%esi),%%xmm4	\n\t"\
		"movaps	0x060(%%esi),%%xmm6	\n\t"\
		"movaps	0x030(%%esi),%%xmm5	\n\t"\
		"movaps	0x070(%%esi),%%xmm7	\n\t"\
		"movaps	%%xmm4,%%xmm0		\n\t"\
		"movaps	%%xmm6,%%xmm2		\n\t"\
		"movaps	%%xmm5,%%xmm1		\n\t"\
		"movaps	%%xmm7,%%xmm3		\n\t"\
		"							\n\t"\
		"							\n\t"\
		"mulpd	     (%%edi),%%xmm4	\n\t"\
		"mulpd	0x010(%%edi),%%xmm6	\n\t"\
		"mulpd	0x010(%%edi),%%xmm1	\n\t"\
		"mulpd	     (%%edi),%%xmm3	\n\t"\
		"mulpd	     (%%edi),%%xmm5	\n\t"\
		"mulpd	0x010(%%edi),%%xmm7	\n\t"\
		"mulpd	0x010(%%edi),%%xmm0	\n\t"\
		"mulpd	     (%%edi),%%xmm2	\n\t"\
		"subpd	%%xmm1,%%xmm4		\n\t"\
		"subpd	%%xmm3,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm5		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"							\n\t"\
		"subpd	%%xmm6,%%xmm4		\n\t"\
		"subpd	%%xmm7,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm6		\n\t"\
		"addpd	%%xmm5,%%xmm7		\n\t"\
		"							\n\t"\
		"movl	%[__isrt2],%%edi	\n\t"\
		"movaps	0x040(%%esi),%%xmm2	\n\t"\
		"movaps	0x050(%%esi),%%xmm3	\n\t"\
		"subpd	0x050(%%esi),%%xmm2	\n\t"\
		"addpd	0x040(%%esi),%%xmm3	\n\t"\
		"mulpd	(%%edi),%%xmm2		\n\t"\
		"mulpd	(%%edi),%%xmm3		\n\t"\
		"							\n\t"\
		"movaps	     (%%esi),%%xmm0	\n\t"\
		"movaps	0x010(%%esi),%%xmm1	\n\t"\
		"							\n\t"\
		"movl	%[__p4],%%edi		\n\t"\
		"shll	$3,%%edi			\n\t"\
		"							\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"addpd	     (%%esi),%%xmm2	\n\t"\
		"addpd	0x010(%%esi),%%xmm3	\n\t"\
		"							\n\t"\
		"subpd	%%xmm6,%%xmm2		\n\t"\
		"subpd	%%xmm7,%%xmm3		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"movaps	%%xmm2,    (%%ebx)	\n\t"\
		"movaps	%%xmm3,0x10(%%ebx)	\n\t"\
		"addpd	%%xmm2,		%%xmm6	\n\t"\
		"addpd	%%xmm3,		%%xmm7	\n\t"\
		"movaps	%%xmm6,    (%%eax)	\n\t"\
		"movaps	%%xmm7,0x10(%%eax)	\n\t"\
		"							\n\t"\
		"subpd	%%xmm5,		%%xmm0	\n\t"\
		"subpd	%%xmm4,		%%xmm1	\n\t"\
		"addpd	%%xmm5,		%%xmm5	\n\t"\
		"addpd	%%xmm4,		%%xmm4	\n\t"\
		"movaps	%%xmm0,    (%%ecx)	\n\t"\
		"movaps	%%xmm1,0x10(%%edx)	\n\t"\
		"addpd	%%xmm0,		%%xmm5	\n\t"\
		"addpd	%%xmm1,		%%xmm4	\n\t"\
		"movaps	%%xmm5,    (%%edx)	\n\t"\
		"movaps	%%xmm4,0x10(%%ecx)	\n\t"\
		"							\n\t"\
		"addl	%%edi,%%eax			\n\t"\
		"addl	%%edi,%%ebx			\n\t"\
		"addl	%%edi,%%ecx			\n\t"\
		"addl	%%edi,%%edx			\n\t"\
		"							\n\t"\
		"movl	%[__r25],%%esi		\n\t"\
		"movl	%[__cc0],%%edi		\n\t"\
		"movaps	0x020(%%esi),%%xmm4	\n\t"\
		"movaps	0x060(%%esi),%%xmm6	\n\t"\
		"movaps	0x030(%%esi),%%xmm5	\n\t"\
		"movaps	0x070(%%esi),%%xmm7	\n\t"\
		"movaps	%%xmm4,%%xmm0		\n\t"\
		"movaps	%%xmm6,%%xmm2		\n\t"\
		"movaps	%%xmm5,%%xmm1		\n\t"\
		"movaps	%%xmm7,%%xmm3		\n\t"\
		"							\n\t"\
		"							\n\t"\
		"mulpd	0x010(%%edi),%%xmm4	\n\t"\
		"mulpd	     (%%edi),%%xmm6	\n\t"\
		"mulpd	     (%%edi),%%xmm1	\n\t"\
		"mulpd	0x010(%%edi),%%xmm3	\n\t"\
		"mulpd	0x010(%%edi),%%xmm5	\n\t"\
		"mulpd	     (%%edi),%%xmm7	\n\t"\
		"mulpd	     (%%edi),%%xmm0	\n\t"\
		"mulpd	0x010(%%edi),%%xmm2	\n\t"\
		"subpd	%%xmm1,%%xmm4		\n\t"\
		"subpd	%%xmm3,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm5		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"							\n\t"\
		"subpd	%%xmm6,%%xmm4		\n\t"\
		"subpd	%%xmm7,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm6		\n\t"\
		"addpd	%%xmm5,%%xmm7		\n\t"\
		"							\n\t"\
		"movl	%[__isrt2],%%edi	\n\t"\
		"movaps	0x040(%%esi),%%xmm2	\n\t"\
		"movaps	0x050(%%esi),%%xmm3	\n\t"\
		"addpd	0x050(%%esi),%%xmm2	\n\t"\
		"subpd	0x040(%%esi),%%xmm3	\n\t"\
		"mulpd	(%%edi),%%xmm2		\n\t"\
		"mulpd	(%%edi),%%xmm3		\n\t"\
		"							\n\t"\
		"movaps	     (%%esi),%%xmm0	\n\t"\
		"movaps	0x010(%%esi),%%xmm1	\n\t"\
		"							\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"addpd	     (%%esi),%%xmm2	\n\t"\
		"addpd	0x010(%%esi),%%xmm3	\n\t"\
		"							\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"movaps	%%xmm0,    (%%ebx)	\n\t"\
		"movaps	%%xmm1,0x10(%%ebx)	\n\t"\
		"addpd	%%xmm0,		%%xmm4	\n\t"\
		"addpd	%%xmm1,		%%xmm5	\n\t"\
		"movaps	%%xmm4,    (%%eax)	\n\t"\
		"movaps	%%xmm5,0x10(%%eax)	\n\t"\
		"							\n\t"\
		"subpd	%%xmm7,		%%xmm2	\n\t"\
		"subpd	%%xmm6,		%%xmm3	\n\t"\
		"addpd	%%xmm7,		%%xmm7	\n\t"\
		"addpd	%%xmm6,		%%xmm6	\n\t"\
		"movaps	%%xmm2,    (%%ecx)	\n\t"\
		"movaps	%%xmm3,0x10(%%edx)	\n\t"\
		"addpd	%%xmm2,		%%xmm7	\n\t"\
		"addpd	%%xmm3,		%%xmm6	\n\t"\
		"movaps	%%xmm7,    (%%edx)	\n\t"\
		"movaps	%%xmm6,0x10(%%ecx)	\n\t"\
	"popl %%ebx	\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__p1] "m" (Xp1)\
		 ,[__p2] "m" (Xp2)\
		 ,[__p3] "m" (Xp3)\
		 ,[__p4] "m" (Xp4)\
		 ,[__r1] "m" (Xr1)\
		 ,[__r3] "m" (Xr3)\
		 ,[__r5] "m" (Xr5)\
		 ,[__r7] "m" (Xr7)\
		 ,[__r9] "m" (Xr9)\
		 ,[__r11] "m" (Xr11)\
		 ,[__r13] "m" (Xr13)\
		 ,[__r15] "m" (Xr15)\
		 ,[__r17] "m" (Xr17)\
		 ,[__r19] "m" (Xr19)\
		 ,[__r21] "m" (Xr21)\
		 ,[__r23] "m" (Xr23)\
		 ,[__r25] "m" (Xr25)\
		 ,[__r27] "m" (Xr27)\
		 ,[__r29] "m" (Xr29)\
		 ,[__r31] "m" (Xr31)\
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__cc0] "m" (Xcc0)\
		: "cc","memory","eax",/*"ebx",*/"ecx","edx","edi","esi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"	/* Clobbered registers */\
	);\
	}

#endif	/* radix16_ditN_cy_dif1_gcc_h_included */

