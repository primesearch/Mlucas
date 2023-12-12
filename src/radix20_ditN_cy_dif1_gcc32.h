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
#ifndef radix20_ditN_cy_dif1_gcc_h_included
#define radix20_ditN_cy_dif1_gcc_h_included

	#define	SSE2_RADIX20_DIT_NOTWIDDLE(Xadd0,Xp01,Xp04,Xr00,Xr10,Xr20,Xr30,Xcc1,Xa0,Xa1,Xa2,Xa3,Xa4,Xb0,Xb1,Xb2,Xb3,Xb4,Xc0,Xc1,Xc2,Xc3,Xc4,Xd0,Xd1,Xd2,Xd3,Xd4)\
	{\
	__asm__ volatile (\
	"pushl %%ebx	\n\t"/* Explicit save/restore of PIC register */\
		"/* Outputs in SSE2 modes are temps 2*5*16 = 10*16 = 0x0a0 bytes apart: */\n\t"\
		"movl	%[__add0],%%eax	/* Must use eax as base address throughout, since that is preserved in SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_C */\n\t"\
		"movl	%[__add0],%%edx		\n\t"\
		"movl	%[__p01],%%esi	/* esi will store power-of-2 multiples of p01 throughout */\n\t"\
		"shll	$3,%%esi		/* Pointer offset for floating doubles */\n\t"\
		"addl	%%esi,%%edx			\n\t"\
		"movl	%%edx,%%ebx		/* add1 = add0+p01 */\n\t"\
		"addl	%%esi,%%edx			\n\t"\
		"movl	%%edx,%%ecx		/* add3 = add0+p02 */\n\t"\
		"addl	%%esi,%%edx		/* add2 = add0+p03 */\n\t"\
		"/*	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_C(%%eax,ebx,edx,ecx, 0x0a0, 0x140, r00)	*/\n\t"\
		"movaps	0x10(%%eax),%%xmm3	\n\t"\
		"movaps	    (%%eax),%%xmm2	\n\t"\
		"movaps	    (%%ebx),%%xmm0	\n\t"\
		"movaps	0x10(%%ebx),%%xmm1	\n\t"\
		"movaps	    (%%ecx),%%xmm6	\n\t"\
		"movaps	0x10(%%ecx),%%xmm7	\n\t"\
		"movaps	0x10(%%edx),%%xmm5	\n\t"\
		"movaps	    (%%edx),%%xmm4	\n\t"\
		"movl	%[__r00],%%esi		\n\t"\
		"movl	$0x0a0,%%edi		\n\t"\
		"addl	%%esi,%%edi			\n\t"\
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
		"/* Finish radix-4 butterfly and store results into temporary-array slots: */\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm2		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"subpd	%%xmm6,%%xmm3		\n\t"\
		"movaps	%%xmm0,0x140(%%esi)	\n\t"\
		"movaps	%%xmm2,     (%%edi)	\n\t/* swap output location of xmm2 <-> xmm7 */"\
		"movaps	%%xmm1,0x150(%%esi)	\n\t"\
		"movaps	%%xmm3,0x150(%%edi)	\n\t/* swap output location of xmm3 <-> xmm6 */"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"addpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	%%xmm4,     (%%esi)	\n\t"\
		"movaps	%%xmm7,0x140(%%edi)	\n\t"\
		"movaps	%%xmm5,0x010(%%esi)	\n\t"\
		"movaps	%%xmm6,0x010(%%edi)	\n\t"\
		"							\n\t"\
		"/*	add0,1,2,3 = &a[j1+p04]+p3,2,1,0 */\n\t"\
		"movl	%[__p04],%%esi		\n\t"\
		"shll	$3,%%esi			\n\t"\
		"addl	%%esi,%%eax		/* &a[j1+p04] */\n\t"\
		"addl	%%esi,%%ebx			\n\t"\
		"addl	%%esi,%%ecx			\n\t"\
		"addl	%%esi,%%edx			\n\t"\
		"/*	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_C(%%edx,ecx,ebx,eax, 0x0a0, 0x140, r02)	*/\n\t"\
		"movl	%[__r00],%%esi		\n\t"\
		"movaps	0x10(%%edx),%%xmm3	\n\t"\
		"movaps	    (%%edx),%%xmm2	\n\t"\
		"movaps	    (%%ecx),%%xmm0	\n\t"\
		"movaps	0x10(%%ecx),%%xmm1	\n\t"\
		"movaps	    (%%ebx),%%xmm6	\n\t"\
		"movaps	0x10(%%ebx),%%xmm7	\n\t"\
		"movaps	0x10(%%eax),%%xmm5	\n\t"\
		"movaps	    (%%eax),%%xmm4	\n\t"\
		"addl	$0x020,%%esi	/* r02 */\n\t"\
		"movl	$0x0a0,%%edi		\n\t"\
		"addl	%%esi,%%edi			\n\t"\
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
		"/* Finish radix-4 butterfly and store results into temporary-array slots: */\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm2		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"subpd	%%xmm6,%%xmm3		\n\t"\
		"movaps	%%xmm0,0x140(%%esi)	\n\t"\
		"movaps	%%xmm2,0x140(%%edi)	\n\t"\
		"movaps	%%xmm1,0x150(%%esi)	\n\t"\
		"movaps	%%xmm3,0x010(%%edi)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"addpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	%%xmm4,     (%%esi)	\n\t"\
		"movaps	%%xmm7,     (%%edi)	\n\t"\
		"movaps	%%xmm5,0x010(%%esi)	\n\t"\
		"movaps	%%xmm6,0x150(%%edi)	\n\t"\
		"							\n\t"\
		"/*	add0,1,2,3 = &a[j1+p08]+p1,0,2,3 */\n\t"\
		"movl	%[__p04],%%esi		\n\t"\
		"shll	$3,%%esi			\n\t"\
		"addl	%%esi,%%eax		/* &a[j1+p08] */\n\t"\
		"addl	%%esi,%%ebx			\n\t"\
		"addl	%%esi,%%ecx			\n\t"\
		"addl	%%esi,%%edx			\n\t"\
		"/*	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_C(%%ebx,eax,ecx,edx, 0x0a0, 0x140, r04)	*/\n\t"\
		"movl	%[__r00],%%esi		\n\t"\
		"movaps	0x10(%%ebx),%%xmm3	\n\t"\
		"movaps	    (%%ebx),%%xmm2	\n\t"\
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"movaps	    (%%ecx),%%xmm6	\n\t"\
		"movaps	0x10(%%ecx),%%xmm7	\n\t"\
		"movaps	0x10(%%edx),%%xmm5	\n\t"\
		"movaps	    (%%edx),%%xmm4	\n\t"\
		"addl	$0x040,%%esi	/* r04 */\n\t"\
		"movl	$0x0a0,%%edi		\n\t"\
		"addl	%%esi,%%edi			\n\t"\
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
		"/* Finish radix-4 butterfly and store results into temporary-array slots: */\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm2		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"subpd	%%xmm6,%%xmm3		\n\t"\
		"movaps	%%xmm0,0x140(%%esi)	\n\t"\
		"movaps	%%xmm2,0x140(%%edi)	\n\t"\
		"movaps	%%xmm1,0x150(%%esi)	\n\t"\
		"movaps	%%xmm3,0x010(%%edi)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"addpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	%%xmm4,     (%%esi)	\n\t"\
		"movaps	%%xmm7,     (%%edi)	\n\t"\
		"movaps	%%xmm5,0x010(%%esi)	\n\t"\
		"movaps	%%xmm6,0x150(%%edi)	\n\t"\
		"							\n\t"\
		"/*	add0,1,2,3 = &a[j1+p12]+p2,3,0,1 */\n\t"\
		"movl	%[__p04],%%esi		\n\t"\
		"shll	$3,%%esi			\n\t"\
		"addl	%%esi,%%eax		/* &a[j1+p12] */\n\t"\
		"addl	%%esi,%%ebx			\n\t"\
		"addl	%%esi,%%ecx			\n\t"\
		"addl	%%esi,%%edx			\n\t"\
		"/*	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_C(%%ecx,edx,eax,ebx, 0x0a0, 0x140, r06)	*/\n\t"\
		"movl	%[__r00],%%esi		\n\t"\
		"movaps	0x10(%%ecx),%%xmm3	\n\t"\
		"movaps	    (%%ecx),%%xmm2	\n\t"\
		"movaps	    (%%edx),%%xmm0	\n\t"\
		"movaps	0x10(%%edx),%%xmm1	\n\t"\
		"movaps	    (%%eax),%%xmm6	\n\t"\
		"movaps	0x10(%%eax),%%xmm7	\n\t"\
		"movaps	0x10(%%ebx),%%xmm5	\n\t"\
		"movaps	    (%%ebx),%%xmm4	\n\t"\
		"addl	$0x060,%%esi	/* r06 */\n\t"\
		"movl	$0x0a0,%%edi		\n\t"\
		"addl	%%esi,%%edi			\n\t"\
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
		"/* Finish radix-4 butterfly and store results into temporary-array slots: */\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm2		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"subpd	%%xmm6,%%xmm3		\n\t"\
		"movaps	%%xmm0,0x140(%%esi)	\n\t"\
		"movaps	%%xmm2,0x140(%%edi)	\n\t"\
		"movaps	%%xmm1,0x150(%%esi)	\n\t"\
		"movaps	%%xmm3,0x010(%%edi)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"addpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	%%xmm4,     (%%esi)	\n\t"\
		"movaps	%%xmm7,     (%%edi)	\n\t"\
		"movaps	%%xmm5,0x010(%%esi)	\n\t"\
		"movaps	%%xmm6,0x150(%%edi)	\n\t"\
		"							\n\t"\
		"/*	add0,1,2,3 = &a[j1+p16]+p0,1,3,2 */\n\t"\
		"movl	%[__p04],%%esi		\n\t"\
		"shll	$3,%%esi			\n\t"\
		"addl	%%esi,%%eax		/* &a[j1+p16] */\n\t"\
		"addl	%%esi,%%ebx			\n\t"\
		"addl	%%esi,%%ecx			\n\t"\
		"addl	%%esi,%%edx			\n\t"\
		"/*	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_C(%%eax,ebx,edx,ecx, 0x0a0, 0x140, r08)	*/\n\t"\
		"movl	%[__r00],%%esi		\n\t"\
		"movaps	0x10(%%eax),%%xmm3	\n\t"\
		"movaps	    (%%eax),%%xmm2	\n\t"\
		"movaps	    (%%ebx),%%xmm0	\n\t"\
		"movaps	0x10(%%ebx),%%xmm1	\n\t"\
		"movaps	    (%%edx),%%xmm6	\n\t"\
		"movaps	0x10(%%edx),%%xmm7	\n\t"\
		"movaps	0x10(%%ecx),%%xmm5	\n\t"\
		"movaps	    (%%ecx),%%xmm4	\n\t"\
		"addl	$0x080,%%esi	/* r08 */\n\t"\
		"movl	$0x0a0,%%edi		\n\t"\
		"addl	%%esi,%%edi			\n\t"\
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
		"/* Finish radix-4 butterfly and store results into temporary-array slots: */\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm2		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"subpd	%%xmm6,%%xmm3		\n\t"\
		"movaps	%%xmm0,0x140(%%esi)	\n\t"\
		"movaps	%%xmm2,0x140(%%edi)	\n\t"\
		"movaps	%%xmm1,0x150(%%esi)	\n\t"\
		"movaps	%%xmm3,0x010(%%edi)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"addpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	%%xmm4,     (%%esi)	\n\t"\
		"movaps	%%xmm7,     (%%edi)	\n\t"\
		"movaps	%%xmm5,0x010(%%esi)	\n\t"\
		"movaps	%%xmm6,0x150(%%edi)	\n\t"\
		"							\n\t"\
		"/* Radix-5 DFT uses adjacenttemps, i.e. stride = 2*16 bytes: */\n\t"\
		"/*	SSE2_RADIX_05_DFT_0TWIDDLE(r00,r02,r04,r06,r08,cc1,Xa0,Xa1,Xa2,Xa3,Xa4)	*/\n\t"\
		"movl	%[__r00],%%esi		\n\t"\
		"movl	%%esi,%%eax			\n\t"\
		"movl	%%esi,%%ebx			\n\t"\
		"movl	%%esi,%%ecx			\n\t"\
		"movl	%%esi,%%edx			\n\t"\
		"addl	$0x020,%%eax		\n\t"\
		"addl	$0x040,%%ebx		\n\t"\
		"addl	$0x060,%%ecx		\n\t"\
		"addl	$0x080,%%edx		\n\t"\
		"movl	%[__a0],%%edi		\n\t"\
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"movaps	    (%%ebx),%%xmm2	\n\t"\
		"movaps	0x10(%%ebx),%%xmm3	\n\t"\
		"movaps	    (%%ecx),%%xmm4	\n\t"\
		"movaps	0x10(%%ecx),%%xmm5	\n\t"\
		"movaps	    (%%edx),%%xmm6	\n\t"\
		"movaps	0x10(%%edx),%%xmm7	\n\t"\
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
		"movl	%[__cc1],%%eax		\n\t"\
		"subpd	%%xmm4,%%xmm6		\n\t"\
		"subpd	%%xmm5,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm5		\n\t"\
		"addpd	    (%%esi),%%xmm4	\n\t"\
		"addpd	0x10(%%esi),%%xmm5	\n\t"\
		"movaps	%%xmm4,    (%%edi)	\n\t"\
		"movaps	%%xmm5,0x10(%%edi)	\n\t"\
		"mulpd	0x10(%%eax),%%xmm6	\n\t"\
		"mulpd	0x10(%%eax),%%xmm7	\n\t"\
		"subpd	    (%%esi),%%xmm4	\n\t"\
		"subpd	0x10(%%esi),%%xmm5	\n\t"\
		"mulpd	    (%%eax),%%xmm4	\n\t"\
		"mulpd	    (%%eax),%%xmm5	\n\t"\
		"addpd	    (%%edi),%%xmm4	\n\t"\
		"addpd	0x10(%%edi),%%xmm5	\n\t"\
		"subpd	%%xmm6,%%xmm4		\n\t"\
		"subpd	%%xmm7,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm6		\n\t"\
		"addpd	%%xmm5,%%xmm7		\n\t"\
		"movaps	%%xmm4,    (%%esi)	\n\t"\
		"movaps	%%xmm5,0x10(%%esi)	\n\t"\
		"movaps	%%xmm0,%%xmm4		\n\t"\
		"movaps	%%xmm1,%%xmm5		\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"mulpd	0x20(%%eax),%%xmm0	\n\t"\
		"mulpd	0x20(%%eax),%%xmm1	\n\t"\
		"mulpd	0x30(%%eax),%%xmm2	\n\t"\
		"mulpd	0x30(%%eax),%%xmm3	\n\t"\
		"mulpd	0x40(%%eax),%%xmm4	\n\t"\
		"mulpd	0x40(%%eax),%%xmm5	\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"movaps	    (%%esi),%%xmm4	\n\t"\
		"movaps	0x10(%%esi),%%xmm5	\n\t"\
		"movl	%[__a1],%%eax		\n\t"\
		"movl	%[__a4],%%edx		\n\t"\
		"subpd	%%xmm3,%%xmm6		\n\t"\
		"subpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"movaps	%%xmm6,    (%%eax)	\n\t"\
		"movaps	%%xmm7,0x10(%%edx)	\n\t"\
		"addpd	%%xmm6,%%xmm3		\n\t"\
		"addpd	%%xmm7,%%xmm2		\n\t"\
		"movaps	%%xmm3,    (%%edx)	\n\t"\
		"movaps	%%xmm2,0x10(%%eax)	\n\t"\
		"movl	%[__a2],%%ebx		\n\t"\
		"movl	%[__a3],%%ecx		\n\t"\
		"subpd	%%xmm1,%%xmm4		\n\t"\
		"subpd	%%xmm0,%%xmm5		\n\t"\
		"addpd	%%xmm1,%%xmm1		\n\t"\
		"addpd	%%xmm0,%%xmm0		\n\t"\
		"movaps	%%xmm4,    (%%ebx)	\n\t"\
		"movaps	%%xmm5,0x10(%%ecx)	\n\t"\
		"addpd	%%xmm4,%%xmm1		\n\t"\
		"addpd	%%xmm5,%%xmm0		\n\t"\
		"movaps	%%xmm1,    (%%ecx)	\n\t"\
		"movaps	%%xmm0,0x10(%%ebx)	\n\t"\
		"							\n\t"\
		"/*	SSE2_RADIX_05_DFT_0TWIDDLE(r10,r12,r14,r16,r18,cc1,Xb0,Xb1,Xb2,Xb3,Xb4)	*/\n\t"\
		"movl	%[__r10],%%esi		\n\t"\
		"movl	%%esi,%%eax			\n\t"\
		"movl	%%esi,%%ebx			\n\t"\
		"movl	%%esi,%%ecx			\n\t"\
		"movl	%%esi,%%edx			\n\t"\
		"addl	$0x020,%%eax		\n\t"\
		"addl	$0x040,%%ebx		\n\t"\
		"addl	$0x060,%%ecx		\n\t"\
		"addl	$0x080,%%edx		\n\t"\
		"movl	%[__b0],%%edi		\n\t"\
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"movaps	    (%%ebx),%%xmm2	\n\t"\
		"movaps	0x10(%%ebx),%%xmm3	\n\t"\
		"movaps	    (%%ecx),%%xmm4	\n\t"\
		"movaps	0x10(%%ecx),%%xmm5	\n\t"\
		"movaps	    (%%edx),%%xmm6	\n\t"\
		"movaps	0x10(%%edx),%%xmm7	\n\t"\
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
		"movl	%[__cc1],%%eax		\n\t"\
		"subpd	%%xmm4,%%xmm6		\n\t"\
		"subpd	%%xmm5,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm5		\n\t"\
		"addpd	    (%%esi),%%xmm4	\n\t"\
		"addpd	0x10(%%esi),%%xmm5	\n\t"\
		"movaps	%%xmm4,    (%%edi)	\n\t"\
		"movaps	%%xmm5,0x10(%%edi)	\n\t"\
		"mulpd	0x10(%%eax),%%xmm6	\n\t"\
		"mulpd	0x10(%%eax),%%xmm7	\n\t"\
		"subpd	    (%%esi),%%xmm4	\n\t"\
		"subpd	0x10(%%esi),%%xmm5	\n\t"\
		"mulpd	    (%%eax),%%xmm4	\n\t"\
		"mulpd	    (%%eax),%%xmm5	\n\t"\
		"addpd	    (%%edi),%%xmm4	\n\t"\
		"addpd	0x10(%%edi),%%xmm5	\n\t"\
		"subpd	%%xmm6,%%xmm4		\n\t"\
		"subpd	%%xmm7,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm6		\n\t"\
		"addpd	%%xmm5,%%xmm7		\n\t"\
		"movaps	%%xmm4,    (%%esi)	\n\t"\
		"movaps	%%xmm5,0x10(%%esi)	\n\t"\
		"movaps	%%xmm0,%%xmm4		\n\t"\
		"movaps	%%xmm1,%%xmm5		\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"mulpd	0x20(%%eax),%%xmm0	\n\t"\
		"mulpd	0x20(%%eax),%%xmm1	\n\t"\
		"mulpd	0x30(%%eax),%%xmm2	\n\t"\
		"mulpd	0x30(%%eax),%%xmm3	\n\t"\
		"mulpd	0x40(%%eax),%%xmm4	\n\t"\
		"mulpd	0x40(%%eax),%%xmm5	\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"movaps	    (%%esi),%%xmm4	\n\t"\
		"movaps	0x10(%%esi),%%xmm5	\n\t"\
		"movl	%[__b1],%%eax		\n\t"\
		"movl	%[__b4],%%edx		\n\t"\
		"subpd	%%xmm3,%%xmm6		\n\t"\
		"subpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"movaps	%%xmm6,    (%%eax)	\n\t"\
		"movaps	%%xmm7,0x10(%%edx)	\n\t"\
		"addpd	%%xmm6,%%xmm3		\n\t"\
		"addpd	%%xmm7,%%xmm2		\n\t"\
		"movaps	%%xmm3,    (%%edx)	\n\t"\
		"movaps	%%xmm2,0x10(%%eax)	\n\t"\
		"movl	%[__b2],%%ebx		\n\t"\
		"movl	%[__b3],%%ecx		\n\t"\
		"subpd	%%xmm1,%%xmm4		\n\t"\
		"subpd	%%xmm0,%%xmm5		\n\t"\
		"addpd	%%xmm1,%%xmm1		\n\t"\
		"addpd	%%xmm0,%%xmm0		\n\t"\
		"movaps	%%xmm4,    (%%ebx)	\n\t"\
		"movaps	%%xmm5,0x10(%%ecx)	\n\t"\
		"addpd	%%xmm4,%%xmm1		\n\t"\
		"addpd	%%xmm5,%%xmm0		\n\t"\
		"movaps	%%xmm1,    (%%ecx)	\n\t"\
		"movaps	%%xmm0,0x10(%%ebx)	\n\t"\
		"							\n\t"\
		"/*	SSE2_RADIX_05_DFT_0TWIDDLE(r20,r22,r24,r26,r28,cc1,Xc0,Xc1,Xc2,Xc3,Xc4)	*/\n\t"\
		"movl	%[__r20],%%esi		\n\t"\
		"movl	%%esi,%%eax			\n\t"\
		"movl	%%esi,%%ebx			\n\t"\
		"movl	%%esi,%%ecx			\n\t"\
		"movl	%%esi,%%edx			\n\t"\
		"addl	$0x020,%%eax		\n\t"\
		"addl	$0x040,%%ebx		\n\t"\
		"addl	$0x060,%%ecx		\n\t"\
		"addl	$0x080,%%edx		\n\t"\
		"movl	%[__c0],%%edi		\n\t"\
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"movaps	    (%%ebx),%%xmm2	\n\t"\
		"movaps	0x10(%%ebx),%%xmm3	\n\t"\
		"movaps	    (%%ecx),%%xmm4	\n\t"\
		"movaps	0x10(%%ecx),%%xmm5	\n\t"\
		"movaps	    (%%edx),%%xmm6	\n\t"\
		"movaps	0x10(%%edx),%%xmm7	\n\t"\
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
		"movl	%[__cc1],%%eax		\n\t"\
		"subpd	%%xmm4,%%xmm6		\n\t"\
		"subpd	%%xmm5,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm5		\n\t"\
		"addpd	    (%%esi),%%xmm4	\n\t"\
		"addpd	0x10(%%esi),%%xmm5	\n\t"\
		"movaps	%%xmm4,    (%%edi)	\n\t"\
		"movaps	%%xmm5,0x10(%%edi)	\n\t"\
		"mulpd	0x10(%%eax),%%xmm6	\n\t"\
		"mulpd	0x10(%%eax),%%xmm7	\n\t"\
		"subpd	    (%%esi),%%xmm4	\n\t"\
		"subpd	0x10(%%esi),%%xmm5	\n\t"\
		"mulpd	    (%%eax),%%xmm4	\n\t"\
		"mulpd	    (%%eax),%%xmm5	\n\t"\
		"addpd	    (%%edi),%%xmm4	\n\t"\
		"addpd	0x10(%%edi),%%xmm5	\n\t"\
		"subpd	%%xmm6,%%xmm4		\n\t"\
		"subpd	%%xmm7,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm6		\n\t"\
		"addpd	%%xmm5,%%xmm7		\n\t"\
		"movaps	%%xmm4,    (%%esi)	\n\t"\
		"movaps	%%xmm5,0x10(%%esi)	\n\t"\
		"movaps	%%xmm0,%%xmm4		\n\t"\
		"movaps	%%xmm1,%%xmm5		\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"mulpd	0x20(%%eax),%%xmm0	\n\t"\
		"mulpd	0x20(%%eax),%%xmm1	\n\t"\
		"mulpd	0x30(%%eax),%%xmm2	\n\t"\
		"mulpd	0x30(%%eax),%%xmm3	\n\t"\
		"mulpd	0x40(%%eax),%%xmm4	\n\t"\
		"mulpd	0x40(%%eax),%%xmm5	\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"movaps	    (%%esi),%%xmm4	\n\t"\
		"movaps	0x10(%%esi),%%xmm5	\n\t"\
		"movl	%[__c1],%%eax		\n\t"\
		"movl	%[__c4],%%edx		\n\t"\
		"subpd	%%xmm3,%%xmm6		\n\t"\
		"subpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"movaps	%%xmm6,    (%%eax)	\n\t"\
		"movaps	%%xmm7,0x10(%%edx)	\n\t"\
		"addpd	%%xmm6,%%xmm3		\n\t"\
		"addpd	%%xmm7,%%xmm2		\n\t"\
		"movaps	%%xmm3,    (%%edx)	\n\t"\
		"movaps	%%xmm2,0x10(%%eax)	\n\t"\
		"movl	%[__c2],%%ebx		\n\t"\
		"movl	%[__c3],%%ecx		\n\t"\
		"subpd	%%xmm1,%%xmm4		\n\t"\
		"subpd	%%xmm0,%%xmm5		\n\t"\
		"addpd	%%xmm1,%%xmm1		\n\t"\
		"addpd	%%xmm0,%%xmm0		\n\t"\
		"movaps	%%xmm4,    (%%ebx)	\n\t"\
		"movaps	%%xmm5,0x10(%%ecx)	\n\t"\
		"addpd	%%xmm4,%%xmm1		\n\t"\
		"addpd	%%xmm5,%%xmm0		\n\t"\
		"movaps	%%xmm1,    (%%ecx)	\n\t"\
		"movaps	%%xmm0,0x10(%%ebx)	\n\t"\
		"							\n\t"\
		"/*	SSE2_RADIX_05_DFT_0TWIDDLE(r30,r32,r34,r36,r38,cc1,Xd0,Xd1,Xd2,Xd3,Xd4)	*/\n\t"\
		"movl	%[__r30],%%esi		\n\t"\
		"movl	%%esi,%%eax			\n\t"\
		"movl	%%esi,%%ebx			\n\t"\
		"movl	%%esi,%%ecx			\n\t"\
		"movl	%%esi,%%edx			\n\t"\
		"addl	$0x020,%%eax		\n\t"\
		"addl	$0x040,%%ebx		\n\t"\
		"addl	$0x060,%%ecx		\n\t"\
		"addl	$0x080,%%edx		\n\t"\
		"movl	%[__d0],%%edi		\n\t"\
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"movaps	    (%%ebx),%%xmm2	\n\t"\
		"movaps	0x10(%%ebx),%%xmm3	\n\t"\
		"movaps	    (%%ecx),%%xmm4	\n\t"\
		"movaps	0x10(%%ecx),%%xmm5	\n\t"\
		"movaps	    (%%edx),%%xmm6	\n\t"\
		"movaps	0x10(%%edx),%%xmm7	\n\t"\
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
		"movl	%[__cc1],%%eax		\n\t"\
		"subpd	%%xmm4,%%xmm6		\n\t"\
		"subpd	%%xmm5,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm5		\n\t"\
		"addpd	    (%%esi),%%xmm4	\n\t"\
		"addpd	0x10(%%esi),%%xmm5	\n\t"\
		"movaps	%%xmm4,    (%%edi)	\n\t"\
		"movaps	%%xmm5,0x10(%%edi)	\n\t"\
		"mulpd	0x10(%%eax),%%xmm6	\n\t"\
		"mulpd	0x10(%%eax),%%xmm7	\n\t"\
		"subpd	    (%%esi),%%xmm4	\n\t"\
		"subpd	0x10(%%esi),%%xmm5	\n\t"\
		"mulpd	    (%%eax),%%xmm4	\n\t"\
		"mulpd	    (%%eax),%%xmm5	\n\t"\
		"addpd	    (%%edi),%%xmm4	\n\t"\
		"addpd	0x10(%%edi),%%xmm5	\n\t"\
		"subpd	%%xmm6,%%xmm4		\n\t"\
		"subpd	%%xmm7,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm6		\n\t"\
		"addpd	%%xmm5,%%xmm7		\n\t"\
		"movaps	%%xmm4,    (%%esi)	\n\t"\
		"movaps	%%xmm5,0x10(%%esi)	\n\t"\
		"movaps	%%xmm0,%%xmm4		\n\t"\
		"movaps	%%xmm1,%%xmm5		\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"mulpd	0x20(%%eax),%%xmm0	\n\t"\
		"mulpd	0x20(%%eax),%%xmm1	\n\t"\
		"mulpd	0x30(%%eax),%%xmm2	\n\t"\
		"mulpd	0x30(%%eax),%%xmm3	\n\t"\
		"mulpd	0x40(%%eax),%%xmm4	\n\t"\
		"mulpd	0x40(%%eax),%%xmm5	\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"movaps	    (%%esi),%%xmm4	\n\t"\
		"movaps	0x10(%%esi),%%xmm5	\n\t"\
		"movl	%[__d1],%%eax		\n\t"\
		"movl	%[__d4],%%edx		\n\t"\
		"subpd	%%xmm3,%%xmm6		\n\t"\
		"subpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"movaps	%%xmm6,    (%%eax)	\n\t"\
		"movaps	%%xmm7,0x10(%%edx)	\n\t"\
		"addpd	%%xmm6,%%xmm3		\n\t"\
		"addpd	%%xmm7,%%xmm2		\n\t"\
		"movaps	%%xmm3,    (%%edx)	\n\t"\
		"movaps	%%xmm2,0x10(%%eax)	\n\t"\
		"movl	%[__d2],%%ebx		\n\t"\
		"movl	%[__d3],%%ecx		\n\t"\
		"subpd	%%xmm1,%%xmm4		\n\t"\
		"subpd	%%xmm0,%%xmm5		\n\t"\
		"addpd	%%xmm1,%%xmm1		\n\t"\
		"addpd	%%xmm0,%%xmm0		\n\t"\
		"movaps	%%xmm4,    (%%ebx)	\n\t"\
		"movaps	%%xmm5,0x10(%%ecx)	\n\t"\
		"addpd	%%xmm4,%%xmm1		\n\t"\
		"addpd	%%xmm5,%%xmm0		\n\t"\
		"movaps	%%xmm1,    (%%ecx)	\n\t"\
		"movaps	%%xmm0,0x10(%%ebx)	\n\t"\
	"popl %%ebx	\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__p01] "m" (Xp01)\
		 ,[__p04] "m" (Xp04)\
		 ,[__r00] "m" (Xr00)\
		 ,[__r10] "m" (Xr10)\
		 ,[__r20] "m" (Xr20)\
		 ,[__r30] "m" (Xr30)\
		 ,[__cc1] "m" (Xcc1)\
		 ,[__a0] "m" (Xa0)\
		 ,[__a1] "m" (Xa1)\
		 ,[__a2] "m" (Xa2)\
		 ,[__a3] "m" (Xa3)\
		 ,[__a4] "m" (Xa4)\
		 ,[__b0] "m" (Xb0)\
		 ,[__b1] "m" (Xb1)\
		 ,[__b2] "m" (Xb2)\
		 ,[__b3] "m" (Xb3)\
		 ,[__b4] "m" (Xb4)\
		 ,[__c0] "m" (Xc0)\
		 ,[__c1] "m" (Xc1)\
		 ,[__c2] "m" (Xc2)\
		 ,[__c3] "m" (Xc3)\
		 ,[__c4] "m" (Xc4)\
		 ,[__d0] "m" (Xd0)\
		 ,[__d1] "m" (Xd1)\
		 ,[__d2] "m" (Xd2)\
		 ,[__d3] "m" (Xd3)\
		 ,[__d4] "m" (Xd4)\
		: "cc","memory","eax",/*"ebx",*/"ecx","edx","edi","esi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"	/* Clobbered registers */\
	);\
	}


	#define	SSE2_RADIX20_DIF_NOTWIDDLE(Xadd0,Xp01,Xp04,Xp08,Xp16,Xa0,Xa1,Xa2,Xa3,Xa4,Xb0,Xb1,Xb2,Xb3,Xb4,Xc0,Xc1,Xc2,Xc3,Xc4,Xd0,Xd1,Xd2,Xd3,Xd4,Xcc1,Xr00,Xr10,Xr20,Xr30)\
	{\
	__asm__ volatile (\
	"pushl %%ebx	\n\t"/* Explicit save/restore of PIC register */\
		"/*	SSE2_RADIX_05_DFT_0TWIDDLE(Xa0,Xa1,Xa2,Xa3,Xa4,cc1,r00,r02,r04,r08,r06) */\n\t"\
		"movl	%[__a0],%%esi		\n\t"\
		"movl	%[__a1],%%eax		\n\t"\
		"movl	%[__a2],%%ebx		\n\t"\
		"movl	%[__a3],%%ecx		\n\t"\
		"movl	%[__a4],%%edx		\n\t"\
		"movl	%[__r00],%%edi		\n\t"\
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"movaps	    (%%ebx),%%xmm2	\n\t"\
		"movaps	0x10(%%ebx),%%xmm3	\n\t"\
		"movaps	    (%%ecx),%%xmm4	\n\t"\
		"movaps	0x10(%%ecx),%%xmm5	\n\t"\
		"movaps	    (%%edx),%%xmm6	\n\t"\
		"movaps	0x10(%%edx),%%xmm7	\n\t"\
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
		"movl	%[__cc1],%%eax		\n\t"\
		"subpd	%%xmm4,%%xmm6		\n\t"\
		"subpd	%%xmm5,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm5		\n\t"\
		"addpd	    (%%esi),%%xmm4	\n\t"\
		"addpd	0x10(%%esi),%%xmm5	\n\t"\
		"movaps	%%xmm4,    (%%edi)	\n\t"\
		"movaps	%%xmm5,0x10(%%edi)	\n\t"\
		"mulpd	0x10(%%eax),%%xmm6	\n\t"\
		"mulpd	0x10(%%eax),%%xmm7	\n\t"\
		"subpd	    (%%esi),%%xmm4	\n\t"\
		"subpd	0x10(%%esi),%%xmm5	\n\t"\
		"mulpd	    (%%eax),%%xmm4	\n\t"\
		"mulpd	    (%%eax),%%xmm5	\n\t"\
		"addpd	    (%%edi),%%xmm4	\n\t"\
		"addpd	0x10(%%edi),%%xmm5	\n\t"\
		"subpd	%%xmm6,%%xmm4		\n\t"\
		"subpd	%%xmm7,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm6		\n\t"\
		"addpd	%%xmm5,%%xmm7		\n\t"\
		"movaps	%%xmm4,    (%%esi)	\n\t"\
		"movaps	%%xmm5,0x10(%%esi)	\n\t"\
		"movaps	%%xmm0,%%xmm4		\n\t"\
		"movaps	%%xmm1,%%xmm5		\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"mulpd	0x20(%%eax),%%xmm0	\n\t"\
		"mulpd	0x20(%%eax),%%xmm1	\n\t"\
		"mulpd	0x30(%%eax),%%xmm2	\n\t"\
		"mulpd	0x30(%%eax),%%xmm3	\n\t"\
		"mulpd	0x40(%%eax),%%xmm4	\n\t"\
		"mulpd	0x40(%%eax),%%xmm5	\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"movaps	    (%%esi),%%xmm4	\n\t"\
		"movaps	0x10(%%esi),%%xmm5	\n\t"\
		"movl	%%edi,%%eax			\n\t"\
		"movl	%%edi,%%edx			\n\t"\
		"addl	$0x020,%%eax		\n\t"\
		"addl	$0x060,%%edx		\n\t/* o3 and o4 outputs swapped */"\
		"subpd	%%xmm3,%%xmm6		\n\t"\
		"subpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"movaps	%%xmm6,    (%%eax)	\n\t"\
		"movaps	%%xmm7,0x10(%%edx)	\n\t"\
		"addpd	%%xmm6,%%xmm3		\n\t"\
		"addpd	%%xmm7,%%xmm2		\n\t"\
		"movaps	%%xmm3,    (%%edx)	\n\t"\
		"movaps	%%xmm2,0x10(%%eax)	\n\t"\
		"movl	%%edi,%%ebx			\n\t"\
		"movl	%%edi,%%ecx			\n\t"\
		"addl	$0x040,%%ebx		\n\t"\
		"addl	$0x080,%%ecx		\n\t/* o3 and o4 outputs swapped */"\
		"subpd	%%xmm1,%%xmm4		\n\t"\
		"subpd	%%xmm0,%%xmm5		\n\t"\
		"addpd	%%xmm1,%%xmm1		\n\t"\
		"addpd	%%xmm0,%%xmm0		\n\t"\
		"movaps	%%xmm4,    (%%ebx)	\n\t"\
		"movaps	%%xmm5,0x10(%%ecx)	\n\t"\
		"addpd	%%xmm4,%%xmm1		\n\t"\
		"addpd	%%xmm5,%%xmm0		\n\t"\
		"movaps	%%xmm1,    (%%ecx)	\n\t"\
		"movaps	%%xmm0,0x10(%%ebx)	\n\t"\
		"							\n\t"\
		"/*	SSE2_RADIX_05_DFT_0TWIDDLE(Xb0,Xb1,Xb2,Xb3,Xb4,cc1,r10,r12,r14,r18,r16)	*/\n\t"\
		"movl	%[__b0],%%esi		\n\t"\
		"movl	%[__b1],%%eax		\n\t"\
		"movl	%[__b2],%%ebx		\n\t"\
		"movl	%[__b3],%%ecx		\n\t"\
		"movl	%[__b4],%%edx		\n\t"\
		"movl	%[__r10],%%edi		\n\t"\
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"movaps	    (%%ebx),%%xmm2	\n\t"\
		"movaps	0x10(%%ebx),%%xmm3	\n\t"\
		"movaps	    (%%ecx),%%xmm4	\n\t"\
		"movaps	0x10(%%ecx),%%xmm5	\n\t"\
		"movaps	    (%%edx),%%xmm6	\n\t"\
		"movaps	0x10(%%edx),%%xmm7	\n\t"\
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
		"movl	%[__cc1],%%eax		\n\t"\
		"subpd	%%xmm4,%%xmm6		\n\t"\
		"subpd	%%xmm5,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm5		\n\t"\
		"addpd	    (%%esi),%%xmm4	\n\t"\
		"addpd	0x10(%%esi),%%xmm5	\n\t"\
		"movaps	%%xmm4,    (%%edi)	\n\t"\
		"movaps	%%xmm5,0x10(%%edi)	\n\t"\
		"mulpd	0x10(%%eax),%%xmm6	\n\t"\
		"mulpd	0x10(%%eax),%%xmm7	\n\t"\
		"subpd	    (%%esi),%%xmm4	\n\t"\
		"subpd	0x10(%%esi),%%xmm5	\n\t"\
		"mulpd	    (%%eax),%%xmm4	\n\t"\
		"mulpd	    (%%eax),%%xmm5	\n\t"\
		"addpd	    (%%edi),%%xmm4	\n\t"\
		"addpd	0x10(%%edi),%%xmm5	\n\t"\
		"subpd	%%xmm6,%%xmm4		\n\t"\
		"subpd	%%xmm7,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm6		\n\t"\
		"addpd	%%xmm5,%%xmm7		\n\t"\
		"movaps	%%xmm4,    (%%esi)	\n\t"\
		"movaps	%%xmm5,0x10(%%esi)	\n\t"\
		"movaps	%%xmm0,%%xmm4		\n\t"\
		"movaps	%%xmm1,%%xmm5		\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"mulpd	0x20(%%eax),%%xmm0	\n\t"\
		"mulpd	0x20(%%eax),%%xmm1	\n\t"\
		"mulpd	0x30(%%eax),%%xmm2	\n\t"\
		"mulpd	0x30(%%eax),%%xmm3	\n\t"\
		"mulpd	0x40(%%eax),%%xmm4	\n\t"\
		"mulpd	0x40(%%eax),%%xmm5	\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"movaps	    (%%esi),%%xmm4	\n\t"\
		"movaps	0x10(%%esi),%%xmm5	\n\t"\
		"movl	%%edi,%%eax			\n\t"\
		"movl	%%edi,%%edx			\n\t"\
		"addl	$0x020,%%eax		\n\t"\
		"addl	$0x060,%%edx		\n\t/* o3 and o4 outputs swapped */"\
		"subpd	%%xmm3,%%xmm6		\n\t"\
		"subpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"movaps	%%xmm6,    (%%eax)	\n\t"\
		"movaps	%%xmm7,0x10(%%edx)	\n\t"\
		"addpd	%%xmm6,%%xmm3		\n\t"\
		"addpd	%%xmm7,%%xmm2		\n\t"\
		"movaps	%%xmm3,    (%%edx)	\n\t"\
		"movaps	%%xmm2,0x10(%%eax)	\n\t"\
		"movl	%%edi,%%ebx			\n\t"\
		"movl	%%edi,%%ecx			\n\t"\
		"addl	$0x040,%%ebx		\n\t"\
		"addl	$0x080,%%ecx		\n\t/* o3 and o4 outputs swapped */"\
		"subpd	%%xmm1,%%xmm4		\n\t"\
		"subpd	%%xmm0,%%xmm5		\n\t"\
		"addpd	%%xmm1,%%xmm1		\n\t"\
		"addpd	%%xmm0,%%xmm0		\n\t"\
		"movaps	%%xmm4,    (%%ebx)	\n\t"\
		"movaps	%%xmm5,0x10(%%ecx)	\n\t"\
		"addpd	%%xmm4,%%xmm1		\n\t"\
		"addpd	%%xmm5,%%xmm0		\n\t"\
		"movaps	%%xmm1,    (%%ecx)	\n\t"\
		"movaps	%%xmm0,0x10(%%ebx)	\n\t"\
		"							\n\t"\
		"/*	SSE2_RADIX_05_DFT_0TWIDDLE(Xc0,Xc1,Xc2,Xc3,Xc4,cc1,r20,r22,r24,r28,r26)	*/\n\t"\
		"movl	%[__c0],%%esi		\n\t"\
		"movl	%[__c1],%%eax		\n\t"\
		"movl	%[__c2],%%ebx		\n\t"\
		"movl	%[__c3],%%ecx		\n\t"\
		"movl	%[__c4],%%edx		\n\t"\
		"movl	%[__r20],%%edi		\n\t"\
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"movaps	    (%%ebx),%%xmm2	\n\t"\
		"movaps	0x10(%%ebx),%%xmm3	\n\t"\
		"movaps	    (%%ecx),%%xmm4	\n\t"\
		"movaps	0x10(%%ecx),%%xmm5	\n\t"\
		"movaps	    (%%edx),%%xmm6	\n\t"\
		"movaps	0x10(%%edx),%%xmm7	\n\t"\
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
		"movl	%[__cc1],%%eax		\n\t"\
		"subpd	%%xmm4,%%xmm6		\n\t"\
		"subpd	%%xmm5,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm5		\n\t"\
		"addpd	    (%%esi),%%xmm4	\n\t"\
		"addpd	0x10(%%esi),%%xmm5	\n\t"\
		"movaps	%%xmm4,    (%%edi)	\n\t"\
		"movaps	%%xmm5,0x10(%%edi)	\n\t"\
		"mulpd	0x10(%%eax),%%xmm6	\n\t"\
		"mulpd	0x10(%%eax),%%xmm7	\n\t"\
		"subpd	    (%%esi),%%xmm4	\n\t"\
		"subpd	0x10(%%esi),%%xmm5	\n\t"\
		"mulpd	    (%%eax),%%xmm4	\n\t"\
		"mulpd	    (%%eax),%%xmm5	\n\t"\
		"addpd	    (%%edi),%%xmm4	\n\t"\
		"addpd	0x10(%%edi),%%xmm5	\n\t"\
		"subpd	%%xmm6,%%xmm4		\n\t"\
		"subpd	%%xmm7,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm6		\n\t"\
		"addpd	%%xmm5,%%xmm7		\n\t"\
		"movaps	%%xmm4,    (%%esi)	\n\t"\
		"movaps	%%xmm5,0x10(%%esi)	\n\t"\
		"movaps	%%xmm0,%%xmm4		\n\t"\
		"movaps	%%xmm1,%%xmm5		\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"mulpd	0x20(%%eax),%%xmm0	\n\t"\
		"mulpd	0x20(%%eax),%%xmm1	\n\t"\
		"mulpd	0x30(%%eax),%%xmm2	\n\t"\
		"mulpd	0x30(%%eax),%%xmm3	\n\t"\
		"mulpd	0x40(%%eax),%%xmm4	\n\t"\
		"mulpd	0x40(%%eax),%%xmm5	\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"movaps	    (%%esi),%%xmm4	\n\t"\
		"movaps	0x10(%%esi),%%xmm5	\n\t"\
		"movl	%%edi,%%eax			\n\t"\
		"movl	%%edi,%%edx			\n\t"\
		"addl	$0x020,%%eax		\n\t"\
		"addl	$0x060,%%edx		\n\t/* o3 and o4 outputs swapped */"\
		"subpd	%%xmm3,%%xmm6		\n\t"\
		"subpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"movaps	%%xmm6,    (%%eax)	\n\t"\
		"movaps	%%xmm7,0x10(%%edx)	\n\t"\
		"addpd	%%xmm6,%%xmm3		\n\t"\
		"addpd	%%xmm7,%%xmm2		\n\t"\
		"movaps	%%xmm3,    (%%edx)	\n\t"\
		"movaps	%%xmm2,0x10(%%eax)	\n\t"\
		"movl	%%edi,%%ebx			\n\t"\
		"movl	%%edi,%%ecx			\n\t"\
		"addl	$0x040,%%ebx		\n\t"\
		"addl	$0x080,%%ecx		\n\t/* o3 and o4 outputs swapped */"\
		"subpd	%%xmm1,%%xmm4		\n\t"\
		"subpd	%%xmm0,%%xmm5		\n\t"\
		"addpd	%%xmm1,%%xmm1		\n\t"\
		"addpd	%%xmm0,%%xmm0		\n\t"\
		"movaps	%%xmm4,    (%%ebx)	\n\t"\
		"movaps	%%xmm5,0x10(%%ecx)	\n\t"\
		"addpd	%%xmm4,%%xmm1		\n\t"\
		"addpd	%%xmm5,%%xmm0		\n\t"\
		"movaps	%%xmm1,    (%%ecx)	\n\t"\
		"movaps	%%xmm0,0x10(%%ebx)	\n\t"\
		"							\n\t"\
		"/*	SSE2_RADIX_05_DFT_0TWIDDLE(Xd0,Xd1,Xd2,Xd3,Xd4,cc1,r30,r32,r34,r38,r36)	*/\n\t"\
		"movl	%[__d0],%%esi		\n\t"\
		"movl	%[__d1],%%eax		\n\t"\
		"movl	%[__d2],%%ebx		\n\t"\
		"movl	%[__d3],%%ecx		\n\t"\
		"movl	%[__d4],%%edx		\n\t"\
		"movl	%[__r30],%%edi		\n\t"\
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"movaps	    (%%ebx),%%xmm2	\n\t"\
		"movaps	0x10(%%ebx),%%xmm3	\n\t"\
		"movaps	    (%%ecx),%%xmm4	\n\t"\
		"movaps	0x10(%%ecx),%%xmm5	\n\t"\
		"movaps	    (%%edx),%%xmm6	\n\t"\
		"movaps	0x10(%%edx),%%xmm7	\n\t"\
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
		"movl	%[__cc1],%%eax		\n\t"\
		"subpd	%%xmm4,%%xmm6		\n\t"\
		"subpd	%%xmm5,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm5		\n\t"\
		"addpd	    (%%esi),%%xmm4	\n\t"\
		"addpd	0x10(%%esi),%%xmm5	\n\t"\
		"movaps	%%xmm4,    (%%edi)	\n\t"\
		"movaps	%%xmm5,0x10(%%edi)	\n\t"\
		"mulpd	0x10(%%eax),%%xmm6	\n\t"\
		"mulpd	0x10(%%eax),%%xmm7	\n\t"\
		"subpd	    (%%esi),%%xmm4	\n\t"\
		"subpd	0x10(%%esi),%%xmm5	\n\t"\
		"mulpd	    (%%eax),%%xmm4	\n\t"\
		"mulpd	    (%%eax),%%xmm5	\n\t"\
		"addpd	    (%%edi),%%xmm4	\n\t"\
		"addpd	0x10(%%edi),%%xmm5	\n\t"\
		"subpd	%%xmm6,%%xmm4		\n\t"\
		"subpd	%%xmm7,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm6		\n\t"\
		"addpd	%%xmm5,%%xmm7		\n\t"\
		"movaps	%%xmm4,    (%%esi)	\n\t"\
		"movaps	%%xmm5,0x10(%%esi)	\n\t"\
		"movaps	%%xmm0,%%xmm4		\n\t"\
		"movaps	%%xmm1,%%xmm5		\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"mulpd	0x20(%%eax),%%xmm0	\n\t"\
		"mulpd	0x20(%%eax),%%xmm1	\n\t"\
		"mulpd	0x30(%%eax),%%xmm2	\n\t"\
		"mulpd	0x30(%%eax),%%xmm3	\n\t"\
		"mulpd	0x40(%%eax),%%xmm4	\n\t"\
		"mulpd	0x40(%%eax),%%xmm5	\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"movaps	    (%%esi),%%xmm4	\n\t"\
		"movaps	0x10(%%esi),%%xmm5	\n\t"\
		"movl	%%edi,%%eax			\n\t"\
		"movl	%%edi,%%edx			\n\t"\
		"addl	$0x020,%%eax		\n\t"\
		"addl	$0x060,%%edx		\n\t/* o3 and o4 outputs swapped */"\
		"subpd	%%xmm3,%%xmm6		\n\t"\
		"subpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"movaps	%%xmm6,    (%%eax)	\n\t"\
		"movaps	%%xmm7,0x10(%%edx)	\n\t"\
		"addpd	%%xmm6,%%xmm3		\n\t"\
		"addpd	%%xmm7,%%xmm2		\n\t"\
		"movaps	%%xmm3,    (%%edx)	\n\t"\
		"movaps	%%xmm2,0x10(%%eax)	\n\t"\
		"movl	%%edi,%%ebx			\n\t"\
		"movl	%%edi,%%ecx			\n\t"\
		"addl	$0x040,%%ebx		\n\t"\
		"addl	$0x080,%%ecx		\n\t/* o3 and o4 outputs swapped */"\
		"subpd	%%xmm1,%%xmm4		\n\t"\
		"subpd	%%xmm0,%%xmm5		\n\t"\
		"addpd	%%xmm1,%%xmm1		\n\t"\
		"addpd	%%xmm0,%%xmm0		\n\t"\
		"movaps	%%xmm4,    (%%ebx)	\n\t"\
		"movaps	%%xmm5,0x10(%%ecx)	\n\t"\
		"addpd	%%xmm4,%%xmm1		\n\t"\
		"addpd	%%xmm5,%%xmm0		\n\t"\
		"movaps	%%xmm1,    (%%ecx)	\n\t"\
		"movaps	%%xmm0,0x10(%%ebx)	\n\t"\
		"							\n\t"\
		"/*	add0 = &a[j1    ]; 	add1= add0+p01;	add2 = add0+p03;	add3 = add0+p02;	*/\n\t"\
		"movl	%[__add0],%%eax	/* Must use eax as base address throughout, since that is preserved in SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_B */\n\t"\
		"movl	%[__add0],%%edx		\n\t"\
		"movl	%[__p01],%%esi	/* esi will store power-of-2 multiples of p01 throughout */\n\t"\
		"shll	$3,%%esi		/* Pointer offset for floating doubles */\n\t"\
		"addl	%%esi,%%edx			\n\t"\
		"movl	%%edx,%%ebx		/* add0+p01 */\n\t"\
		"addl	%%esi,%%edx			\n\t"\
		"movl	%%edx,%%ecx		/* add0+p02 */\n\t"\
		"addl	%%esi,%%edx		/* add0+p03 */\n\t"\
		"/*	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_C(r00, 0x0a0, 0x140, eax,ebx,edx,ecx) */\n\t"\
		"movl	%[__r00],%%edi		\n\t"\
		"movl	$0x0a0,%%esi		\n\t"\
		"addl	%%edi,%%esi			\n\t"\
		"movaps	    (%%edi),%%xmm4	\n\t"\
		"movaps	    (%%esi),%%xmm6	\n\t"\
		"movaps	0x10(%%edi),%%xmm5	\n\t"\
		"movaps	0x10(%%esi),%%xmm7	\n\t"\
		"addl	$0x140,%%edi		\n\t"\
		"addl	$0x140,%%esi		\n\t"\
		"movaps	    (%%edi),%%xmm0	\n\t"\
		"movaps	    (%%esi),%%xmm2	\n\t"\
		"movaps	0x10(%%edi),%%xmm1	\n\t"\
		"movaps	0x10(%%esi),%%xmm3	\n\t"\
		"subpd	%%xmm0,%%xmm4		\n\t"\
		"subpd	%%xmm2,%%xmm6		\n\t"\
		"subpd	%%xmm1,%%xmm5		\n\t"\
		"subpd	%%xmm3,%%xmm7		\n\t"\
		"addpd	%%xmm0,%%xmm0		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"addpd	%%xmm1,%%xmm1		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm4,%%xmm0		\n\t"\
		"addpd	%%xmm6,%%xmm2		\n\t"\
		"addpd	%%xmm5,%%xmm1		\n\t"\
		"addpd	%%xmm7,%%xmm3		\n\t"\
		"/* Finish radix-4 butterfly and store results into main-array slots: */\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm4		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"subpd	%%xmm6,%%xmm5		\n\t"\
		"movaps	%%xmm0,     (%%ebx)	\n\t"\
		"movaps	%%xmm4,     (%%edx)	\n\t"\
		"movaps	%%xmm1,0x010(%%ebx)	\n\t"\
		"movaps	%%xmm5,0x010(%%ecx)	\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm4,%%xmm7		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"addpd	%%xmm5,%%xmm6		\n\t"\
		"movaps	%%xmm2,     (%%eax)	\n\t"\
		"movaps	%%xmm7,     (%%ecx)	\n\t"\
		"movaps	%%xmm3,0x010(%%eax)	\n\t"\
		"movaps	%%xmm6,0x010(%%edx)	\n\t"\
		"							\n\t"\
		"/*	add0,1,2,3 = &a[j1+p16]+p0,1,3,2 */\n\t"\
		"movl	%[__p16],%%esi		\n\t"\
		"shll	$3,%%esi			\n\t"\
		"addl	%%esi,%%eax		/* &a[j1+p16] */\n\t"\
		"addl	%%esi,%%ebx			\n\t"\
		"addl	%%esi,%%ecx			\n\t"\
		"addl	%%esi,%%edx			\n\t"\
		"/*	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_C(r02, 0x0a0, 0x140, eax,ebx,edx,ecx) */\n\t"\
		"subl	$0x120,%%edi	/* r02 */\n\t"\
		"movl	$0x0a0,%%esi		\n\t"\
		"addl	%%edi,%%esi			\n\t"\
		"movaps	    (%%edi),%%xmm4	\n\t"\
		"movaps	    (%%esi),%%xmm6	\n\t"\
		"movaps	0x10(%%edi),%%xmm5	\n\t"\
		"movaps	0x10(%%esi),%%xmm7	\n\t"\
		"addl	$0x140,%%edi		\n\t"\
		"addl	$0x140,%%esi		\n\t"\
		"movaps	    (%%edi),%%xmm0	\n\t"\
		"movaps	    (%%esi),%%xmm2	\n\t"\
		"movaps	0x10(%%edi),%%xmm1	\n\t"\
		"movaps	0x10(%%esi),%%xmm3	\n\t"\
		"subpd	%%xmm0,%%xmm4		\n\t"\
		"subpd	%%xmm2,%%xmm6		\n\t"\
		"subpd	%%xmm1,%%xmm5		\n\t"\
		"subpd	%%xmm3,%%xmm7		\n\t"\
		"addpd	%%xmm0,%%xmm0		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"addpd	%%xmm1,%%xmm1		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm4,%%xmm0		\n\t"\
		"addpd	%%xmm6,%%xmm2		\n\t"\
		"addpd	%%xmm5,%%xmm1		\n\t"\
		"addpd	%%xmm7,%%xmm3		\n\t"\
		"/* Finish radix-4 butterfly and store results into main-array slots: */\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm4		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"subpd	%%xmm6,%%xmm5		\n\t"\
		"movaps	%%xmm0,     (%%ebx)	\n\t"\
		"movaps	%%xmm4,     (%%edx)	\n\t"\
		"movaps	%%xmm1,0x010(%%ebx)	\n\t"\
		"movaps	%%xmm5,0x010(%%ecx)	\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm4,%%xmm7		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"addpd	%%xmm5,%%xmm6		\n\t"\
		"movaps	%%xmm2,     (%%eax)	\n\t"\
		"movaps	%%xmm7,     (%%ecx)	\n\t"\
		"movaps	%%xmm3,0x010(%%eax)	\n\t"\
		"movaps	%%xmm6,0x010(%%edx)	\n\t"\
		"							\n\t"\
		"/*	add0,1,2,3 = &a[j1+p12]+p2,3,0,1 */\n\t"\
		"movl	%[__p04],%%esi		\n\t"\
		"shll	$3,%%esi			\n\t"\
		"subl	%%esi,%%eax		/* &a[j1+p12] */\n\t"\
		"subl	%%esi,%%ebx			\n\t"\
		"subl	%%esi,%%ecx			\n\t"\
		"subl	%%esi,%%edx			\n\t"\
		"/*	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_C(r04, 0x0a0, 0x140, ecx,edx,eax,ebx) */\n\t"\
		"subl	$0x120,%%edi	/* r04 */\n\t"\
		"movl	$0x0a0,%%esi		\n\t"\
		"addl	%%edi,%%esi			\n\t"\
		"movaps	    (%%edi),%%xmm4	\n\t"\
		"movaps	    (%%esi),%%xmm6	\n\t"\
		"movaps	0x10(%%edi),%%xmm5	\n\t"\
		"movaps	0x10(%%esi),%%xmm7	\n\t"\
		"addl	$0x140,%%edi		\n\t"\
		"addl	$0x140,%%esi		\n\t"\
		"movaps	    (%%edi),%%xmm0	\n\t"\
		"movaps	    (%%esi),%%xmm2	\n\t"\
		"movaps	0x10(%%edi),%%xmm1	\n\t"\
		"movaps	0x10(%%esi),%%xmm3	\n\t"\
		"subpd	%%xmm0,%%xmm4		\n\t"\
		"subpd	%%xmm2,%%xmm6		\n\t"\
		"subpd	%%xmm1,%%xmm5		\n\t"\
		"subpd	%%xmm3,%%xmm7		\n\t"\
		"addpd	%%xmm0,%%xmm0		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"addpd	%%xmm1,%%xmm1		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm4,%%xmm0		\n\t"\
		"addpd	%%xmm6,%%xmm2		\n\t"\
		"addpd	%%xmm5,%%xmm1		\n\t"\
		"addpd	%%xmm7,%%xmm3		\n\t"\
		"/* Finish radix-4 butterfly and store results into main-array slots: */\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm4		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"subpd	%%xmm6,%%xmm5		\n\t"\
		"movaps	%%xmm0,     (%%edx)	\n\t"\
		"movaps	%%xmm4,     (%%eax)	\n\t"\
		"movaps	%%xmm1,0x010(%%edx)	\n\t"\
		"movaps	%%xmm5,0x010(%%ebx)	\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm4,%%xmm7		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"addpd	%%xmm5,%%xmm6		\n\t"\
		"movaps	%%xmm2,     (%%ecx)	\n\t"\
		"movaps	%%xmm7,     (%%ebx)	\n\t"\
		"movaps	%%xmm3,0x010(%%ecx)	\n\t"\
		"movaps	%%xmm6,0x010(%%eax)	\n\t"\
		"							\n\t"\
		"/*	add0,1,2,3 = &a[j1+p04]+p3,2,1,0 */\n\t"\
		"movl	%[__p08],%%esi		\n\t"\
		"shll	$3,%%esi			\n\t"\
		"subl	%%esi,%%eax		/* &a[j1+p04] */\n\t"\
		"subl	%%esi,%%ebx			\n\t"\
		"subl	%%esi,%%ecx			\n\t"\
		"subl	%%esi,%%edx			\n\t"\
		"/*	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_C(r06, 0x0a0, 0x140, edx,ecx,ebx,eax) */\n\t"\
		"subl	$0x120,%%edi	/* r06 */\n\t"\
		"movl	$0x0a0,%%esi		\n\t"\
		"addl	%%edi,%%esi			\n\t"\
		"movaps	    (%%edi),%%xmm4	\n\t"\
		"movaps	    (%%esi),%%xmm6	\n\t"\
		"movaps	0x10(%%edi),%%xmm5	\n\t"\
		"movaps	0x10(%%esi),%%xmm7	\n\t"\
		"addl	$0x140,%%edi		\n\t"\
		"addl	$0x140,%%esi		\n\t"\
		"movaps	    (%%edi),%%xmm0	\n\t"\
		"movaps	    (%%esi),%%xmm2	\n\t"\
		"movaps	0x10(%%edi),%%xmm1	\n\t"\
		"movaps	0x10(%%esi),%%xmm3	\n\t"\
		"subpd	%%xmm0,%%xmm4		\n\t"\
		"subpd	%%xmm2,%%xmm6		\n\t"\
		"subpd	%%xmm1,%%xmm5		\n\t"\
		"subpd	%%xmm3,%%xmm7		\n\t"\
		"addpd	%%xmm0,%%xmm0		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"addpd	%%xmm1,%%xmm1		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm4,%%xmm0		\n\t"\
		"addpd	%%xmm6,%%xmm2		\n\t"\
		"addpd	%%xmm5,%%xmm1		\n\t"\
		"addpd	%%xmm7,%%xmm3		\n\t"\
		"/* Finish radix-4 butterfly and store results into main-array slots: */\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm4		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"subpd	%%xmm6,%%xmm5		\n\t"\
		"movaps	%%xmm0,     (%%ecx)	\n\t"\
		"movaps	%%xmm4,     (%%ebx)	\n\t"\
		"movaps	%%xmm1,0x010(%%ecx)	\n\t"\
		"movaps	%%xmm5,0x010(%%eax)	\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm4,%%xmm7		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"addpd	%%xmm5,%%xmm6		\n\t"\
		"movaps	%%xmm2,     (%%edx)	\n\t"\
		"movaps	%%xmm7,     (%%eax)	\n\t"\
		"movaps	%%xmm3,0x010(%%edx)	\n\t"\
		"movaps	%%xmm6,0x010(%%ebx)	\n\t"\
		"							\n\t"\
		"/*	add0,1,2,3 = &a[j1+p08]+p1,0,2,3 */\n\t"\
		"movl	%[__p04],%%esi		\n\t"\
		"shll	$3,%%esi			\n\t"\
		"addl	%%esi,%%eax		/* &a[j1+p08] */\n\t"\
		"addl	%%esi,%%ebx			\n\t"\
		"addl	%%esi,%%ecx			\n\t"\
		"addl	%%esi,%%edx			\n\t"\
		"/*	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_C(r08, 0x0a0, 0x140, ebx,eax,ecx,edx) */\n\t"\
		"subl	$0x120,%%edi	/* r08 */\n\t"\
		"movl	$0x0a0,%%esi		\n\t"\
		"addl	%%edi,%%esi			\n\t"\
		"movaps	    (%%edi),%%xmm4	\n\t"\
		"movaps	    (%%esi),%%xmm6	\n\t"\
		"movaps	0x10(%%edi),%%xmm5	\n\t"\
		"movaps	0x10(%%esi),%%xmm7	\n\t"\
		"addl	$0x140,%%edi		\n\t"\
		"addl	$0x140,%%esi		\n\t"\
		"movaps	    (%%edi),%%xmm0	\n\t"\
		"movaps	    (%%esi),%%xmm2	\n\t"\
		"movaps	0x10(%%edi),%%xmm1	\n\t"\
		"movaps	0x10(%%esi),%%xmm3	\n\t"\
		"subpd	%%xmm0,%%xmm4		\n\t"\
		"subpd	%%xmm2,%%xmm6		\n\t"\
		"subpd	%%xmm1,%%xmm5		\n\t"\
		"subpd	%%xmm3,%%xmm7		\n\t"\
		"addpd	%%xmm0,%%xmm0		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"addpd	%%xmm1,%%xmm1		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm4,%%xmm0		\n\t"\
		"addpd	%%xmm6,%%xmm2		\n\t"\
		"addpd	%%xmm5,%%xmm1		\n\t"\
		"addpd	%%xmm7,%%xmm3		\n\t"\
		"/* Finish radix-4 butterfly and store results into main-array slots: */\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm4		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"subpd	%%xmm6,%%xmm5		\n\t"\
		"movaps	%%xmm0,     (%%eax)	\n\t"\
		"movaps	%%xmm4,     (%%ecx)	\n\t"\
		"movaps	%%xmm1,0x010(%%eax)	\n\t"\
		"movaps	%%xmm5,0x010(%%edx)	\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm4,%%xmm7		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"addpd	%%xmm5,%%xmm6		\n\t"\
		"movaps	%%xmm2,     (%%ebx)	\n\t"\
		"movaps	%%xmm7,     (%%edx)	\n\t"\
		"movaps	%%xmm3,0x010(%%ebx)	\n\t"\
		"movaps	%%xmm6,0x010(%%ecx)	\n\t"\
	"popl %%ebx	\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__p01] "m" (Xp01)\
		 ,[__p04] "m" (Xp04)\
		 ,[__p08] "m" (Xp08)\
		 ,[__p16] "m" (Xp16)\
		 ,[__a0] "m" (Xa0)\
		 ,[__a1] "m" (Xa1)\
		 ,[__a2] "m" (Xa2)\
		 ,[__a3] "m" (Xa3)\
		 ,[__a4] "m" (Xa4)\
		 ,[__b0] "m" (Xb0)\
		 ,[__b1] "m" (Xb1)\
		 ,[__b2] "m" (Xb2)\
		 ,[__b3] "m" (Xb3)\
		 ,[__b4] "m" (Xb4)\
		 ,[__c0] "m" (Xc0)\
		 ,[__c1] "m" (Xc1)\
		 ,[__c2] "m" (Xc2)\
		 ,[__c3] "m" (Xc3)\
		 ,[__c4] "m" (Xc4)\
		 ,[__d0] "m" (Xd0)\
		 ,[__d1] "m" (Xd1)\
		 ,[__d2] "m" (Xd2)\
		 ,[__d3] "m" (Xd3)\
		 ,[__d4] "m" (Xd4)\
		 ,[__cc1] "m" (Xcc1)\
		 ,[__r00] "m" (Xr00)\
		 ,[__r10] "m" (Xr10)\
		 ,[__r20] "m" (Xr20)\
		 ,[__r30] "m" (Xr30)\
		: "cc","memory","eax",/*"ebx",*/"ecx","edx","edi","esi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"	/* Clobbered registers */\
	);\
	}

#endif	/* radix20_ditN_cy_dif1_gcc_h_included */

