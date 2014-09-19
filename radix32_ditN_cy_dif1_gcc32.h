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
#ifndef radix32_ditN_cy_dif1_gcc_h_included
#define radix32_ditN_cy_dif1_gcc_h_included

	#define SSE2_RADIX32_DIT_NOTWIDDLE(Xadd,Xarr_offsets, Xr00, Xisrt2)\
	{\
	__asm__ volatile (\
	"pushl %%ebx	\n\t"/* Explicit save/restore of PIC register */\
	/*...Block 1: */\
		"movl	%[__r00],%%esi	\n\t"\
		"movl	%[__arr_offsets],%%edi	\n\t"\
		/* Can't simply addl these to a 64-bit base array address since the offsets are 4-bit array data: */\
		"movl	0x0(%%edi),%%eax	\n\t"/* p00 */\
		"movl	0x4(%%edi),%%ebx	\n\t"/* p01 */\
		"movl	0x8(%%edi),%%ecx	\n\t"/* p02 */\
		"movl	0xc(%%edi),%%edx	\n\t"/* p03 */\
		"movl	%[__add],%%edi		\n\t"\
		"addl	%%edi,%%eax	\n\t"/* add0 + p00 */\
		"addl	%%edi,%%ebx	\n\t"/* add0 + p01 */\
		"addl	%%edi,%%ecx	\n\t"/* add0 + p02 */\
		"addl	%%edi,%%edx	\n\t"/* add0 + p03 */\
		"\n\t"/* SSE2_RADIX4_DIT_0TWIDDLE(add0, add1, add2, add3, r00): */\
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
		"\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm2		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"subpd	%%xmm6,%%xmm3		\n\t"\
		"movaps	%%xmm0,0x40(%%esi)	\n\t"\
		"movaps	%%xmm2,0x60(%%esi)	\n\t"\
		"movaps	%%xmm1,0x50(%%esi)	\n\t"\
		"movaps	%%xmm3,0x30(%%esi)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"addpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	%%xmm4,    (%%esi)	\n\t"\
		"movaps	%%xmm7,0x20(%%esi)	\n\t"\
		"movaps	%%xmm5,0x10(%%esi)	\n\t"\
		"movaps	%%xmm6,0x70(%%esi)	\n\t"\
		"\n\t"/* SSE2_RADIX4_DIT_0TWIDDLE_2NDOFTWO(add4, add5, add6, add7, r08): */\
		"addl	$0x80,%%esi			\n\t"\
		"movl	%[__arr_offsets],%%edi	\n\t"\
		"movl	0x10(%%edi),%%eax	\n\t"/* p04 */\
		"movl	0x14(%%edi),%%ebx	\n\t"/* p05 */\
		"movl	0x18(%%edi),%%ecx	\n\t"/* p06 */\
		"movl	0x1c(%%edi),%%edx	\n\t"/* p07 */\
		"movl	%[__add],%%edi		\n\t"\
		"addl	%%edi,%%eax	\n\t"\
		"addl	%%edi,%%ebx	\n\t"\
		"addl	%%edi,%%ecx	\n\t"\
		"addl	%%edi,%%edx	\n\t"\
		"\n\t"\
		"movaps	    (%%eax),%%xmm2	\n\t"\
		"movaps	    (%%ecx),%%xmm6	\n\t"\
		"movaps	0x10(%%eax),%%xmm3	\n\t"\
		"movaps	0x10(%%ecx),%%xmm7	\n\t"\
		"movaps	    (%%ebx),%%xmm0	\n\t"\
		"movaps	    (%%edx),%%xmm4	\n\t"\
		"movaps	0x10(%%ebx),%%xmm1	\n\t"\
		"movaps	0x10(%%edx),%%xmm5	\n\t"\
		"\n\t"\
		"movl	%[__isrt2],%%edi	\n\t"\
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
		"movaps	%%xmm0,0x40(%%esi)	\n\t"\
		"movaps	%%xmm1,0x50(%%esi)	\n\t"\
		"subpd	%%xmm7,%%xmm2		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"subpd	%%xmm6,%%xmm3		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"movaps	%%xmm4,    (%%esi)	\n\t"\
		"addpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	%%xmm5,0x10(%%esi)	\n\t"\
		"\n\t"\
		"movaps	(%%edi),%%xmm5	\n\t"/* ISRT2 */\
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
		"movaps	%%xmm3,0x30(%%esi)	\n\t"\
		"movaps	%%xmm6,0x70(%%esi)	\n\t"\
		"movaps	%%xmm0,0x20(%%esi)	\n\t"\
		"movaps	%%xmm1,0x60(%%esi)	\n\t"\
		"\n\t"/* SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS(r00,r02,r04,r06,r08,r0A,r0C,r0E): */\
		"movaps		-0x80(%%esi),%%xmm0	\n\t"/* r00 */\
		"movaps		-0x40(%%esi),%%xmm4	\n\t"/* r04 */\
		"movaps		-0x70(%%esi),%%xmm1	\n\t"/* r01 */\
		"movaps		-0x30(%%esi),%%xmm5	\n\t"/* r05 */\
		"movaps		     (%%esi),%%xmm2	\n\t"/* r08 */\
		"movaps		 0x50(%%esi),%%xmm7	\n\t"/* r0D */\
		"movaps		 0x10(%%esi),%%xmm3	\n\t"/* r09 */\
		"movaps		 0x40(%%esi),%%xmm6	\n\t"/* r0C */\
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
		"movaps		%%xmm0,     (%%esi)\n\t"\
		"movaps		%%xmm4, 0x40(%%esi)\n\t"\
		"movaps		%%xmm1, 0x10(%%esi)\n\t"\
		"movaps		%%xmm5,-0x30(%%esi)\n\t"\
		"movaps		%%xmm2,-0x80(%%esi)\n\t"\
		"movaps		%%xmm7,-0x40(%%esi)\n\t"\
		"movaps		%%xmm3,-0x70(%%esi)\n\t"\
		"movaps		%%xmm6, 0x50(%%esi)\n\t"\
		"\n\t"\
		"movaps		-0x60(%%esi),%%xmm0	\n\t"/* r02 */\
		"movaps		-0x20(%%esi),%%xmm4	\n\t"/* r06 */\
		"movaps		-0x50(%%esi),%%xmm1	\n\t"/* r03 */\
		"movaps		-0x10(%%esi),%%xmm5	\n\t"/* r07 */\
		"movaps		 0x20(%%esi),%%xmm2	\n\t"/* r0A */\
		"movaps		 0x70(%%esi),%%xmm7	\n\t"/* r0F */\
		"movaps		 0x30(%%esi),%%xmm3	\n\t"/* r0B */\
		"movaps		 0x60(%%esi),%%xmm6	\n\t"/* r0E */\
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
		"movaps		%%xmm0, 0x20(%%esi)\n\t"\
		"movaps		%%xmm4, 0x60(%%esi)\n\t"\
		"movaps		%%xmm1, 0x30(%%esi)\n\t"\
		"movaps		%%xmm5,-0x10(%%esi)\n\t"\
		"movaps		%%xmm2,-0x60(%%esi)\n\t"\
		"movaps		%%xmm7,-0x20(%%esi)\n\t"\
		"movaps		%%xmm3,-0x50(%%esi)\n\t"\
		"movaps		%%xmm6, 0x70(%%esi)\n\t"\
		"\n\t"\
	/*...Block 2: */\
		"addl	$0x80,%%esi			\n\t"\
		"movl	%[__arr_offsets],%%edi	\n\t"\
		"movl	0x20(%%edi),%%eax	\n\t"/* p08 */\
		"movl	0x24(%%edi),%%ebx	\n\t"/* p09 */\
		"movl	0x28(%%edi),%%ecx	\n\t"/* p0a */\
		"movl	0x2c(%%edi),%%edx	\n\t"/* p0b */\
		"movl	%[__add],%%edi		\n\t"\
		"addl	%%edi,%%eax	\n\t"\
		"addl	%%edi,%%ebx	\n\t"\
		"addl	%%edi,%%ecx	\n\t"\
		"addl	%%edi,%%edx	\n\t"\
		"\n\t"\
		"\n\t"/* SSE2_RADIX4_DIT_0TWIDDLE(add0, add1, add2, add3, r10): */\
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
		"\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm2		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"subpd	%%xmm6,%%xmm3		\n\t"\
		"movaps	%%xmm0,0x40(%%esi)	\n\t"\
		"movaps	%%xmm2,0x60(%%esi)	\n\t"\
		"movaps	%%xmm1,0x50(%%esi)	\n\t"\
		"movaps	%%xmm3,0x30(%%esi)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"addpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	%%xmm4,    (%%esi)	\n\t"\
		"movaps	%%xmm7,0x20(%%esi)	\n\t"\
		"movaps	%%xmm5,0x10(%%esi)	\n\t"\
		"movaps	%%xmm6,0x70(%%esi)	\n\t"\
		"\n\t"/* SSE2_RADIX4_DIT_0TWIDDLE_2NDOFTWO(add4, add5, add6, add7, r18): */\
		"addl	$0x80,%%esi			\n\t"\
		"movl	%[__arr_offsets],%%edi	\n\t"\
		"movl	0x30(%%edi),%%eax	\n\t"/* p0c */\
		"movl	0x34(%%edi),%%ebx	\n\t"/* p0d */\
		"movl	0x38(%%edi),%%ecx	\n\t"/* p0e */\
		"movl	0x3c(%%edi),%%edx	\n\t"/* p0f */\
		"movl	%[__add],%%edi		\n\t"\
		"addl	%%edi,%%eax	\n\t"\
		"addl	%%edi,%%ebx	\n\t"\
		"addl	%%edi,%%ecx	\n\t"\
		"addl	%%edi,%%edx	\n\t"\
		"\n\t"\
		"movaps	    (%%eax),%%xmm2	\n\t"\
		"movaps	    (%%ecx),%%xmm6	\n\t"\
		"movaps	0x10(%%eax),%%xmm3	\n\t"\
		"movaps	0x10(%%ecx),%%xmm7	\n\t"\
		"movaps	    (%%ebx),%%xmm0	\n\t"\
		"movaps	    (%%edx),%%xmm4	\n\t"\
		"movaps	0x10(%%ebx),%%xmm1	\n\t"\
		"movaps	0x10(%%edx),%%xmm5	\n\t"\
		"\n\t"\
		"movl	%[__isrt2],%%edi	\n\t"\
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
		"movaps	%%xmm0,0x40(%%esi)	\n\t"\
		"movaps	%%xmm1,0x50(%%esi)	\n\t"\
		"subpd	%%xmm7,%%xmm2		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"subpd	%%xmm6,%%xmm3		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"movaps	%%xmm4,    (%%esi)	\n\t"\
		"addpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	%%xmm5,0x10(%%esi)	\n\t"\
		"\n\t"\
		"movaps	(%%edi),%%xmm5	\n\t"/* ISRT2 */\
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
		"movaps	%%xmm3,0x30(%%esi)	\n\t"\
		"movaps	%%xmm6,0x70(%%esi)	\n\t"\
		"movaps	%%xmm0,0x20(%%esi)	\n\t"\
		"movaps	%%xmm1,0x60(%%esi)	\n\t"\
		"\n\t"/* SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS(r10,r12,r14,r16,r18,r1A,r1C,r1E): */\
		"movaps		-0x80(%%esi),%%xmm0	\n\t"/* r10 */\
		"movaps		-0x40(%%esi),%%xmm4	\n\t"/* r14 */\
		"movaps		-0x70(%%esi),%%xmm1	\n\t"/* r11 */\
		"movaps		-0x30(%%esi),%%xmm5	\n\t"/* r15 */\
		"movaps		     (%%esi),%%xmm2	\n\t"/* r18 */\
		"movaps		 0x50(%%esi),%%xmm7	\n\t"/* r1D */\
		"movaps		 0x10(%%esi),%%xmm3	\n\t"/* r19 */\
		"movaps		 0x40(%%esi),%%xmm6	\n\t"/* r1C */\
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
		"movaps		%%xmm0,     (%%esi)\n\t"\
		"movaps		%%xmm4, 0x40(%%esi)\n\t"\
		"movaps		%%xmm1, 0x10(%%esi)\n\t"\
		"movaps		%%xmm5,-0x30(%%esi)\n\t"\
		"movaps		%%xmm2,-0x80(%%esi)\n\t"\
		"movaps		%%xmm7,-0x40(%%esi)\n\t"\
		"movaps		%%xmm3,-0x70(%%esi)\n\t"\
		"movaps		%%xmm6, 0x50(%%esi)\n\t"\
		"\n\t"\
		"movaps		-0x60(%%esi),%%xmm0	\n\t"/* r12 */\
		"movaps		-0x20(%%esi),%%xmm4	\n\t"/* r16 */\
		"movaps		-0x50(%%esi),%%xmm1	\n\t"/* r13 */\
		"movaps		-0x10(%%esi),%%xmm5	\n\t"/* r17 */\
		"movaps		 0x20(%%esi),%%xmm2	\n\t"/* r1A */\
		"movaps		 0x70(%%esi),%%xmm7	\n\t"/* r1F */\
		"movaps		 0x30(%%esi),%%xmm3	\n\t"/* r1B */\
		"movaps		 0x60(%%esi),%%xmm6	\n\t"/* r1E */\
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
		"movaps		%%xmm0, 0x20(%%esi)\n\t"\
		"movaps		%%xmm4, 0x60(%%esi)\n\t"\
		"movaps		%%xmm1, 0x30(%%esi)\n\t"\
		"movaps		%%xmm5,-0x10(%%esi)\n\t"\
		"movaps		%%xmm2,-0x60(%%esi)\n\t"\
		"movaps		%%xmm7,-0x20(%%esi)\n\t"\
		"movaps		%%xmm3,-0x50(%%esi)\n\t"\
		"movaps		%%xmm6, 0x70(%%esi)\n\t"\
		"\n\t"\
	/*...Block 3: */\
		"addl	$0x80,%%esi			\n\t"\
		"movl	%[__arr_offsets],%%edi	\n\t"\
		"movl	0x40(%%edi),%%eax	\n\t"/* p10 */\
		"movl	0x44(%%edi),%%ebx	\n\t"/* p11 */\
		"movl	0x48(%%edi),%%ecx	\n\t"/* p12 */\
		"movl	0x4c(%%edi),%%edx	\n\t"/* p13 */\
		"movl	%[__add],%%edi		\n\t"\
		"addl	%%edi,%%eax	\n\t"\
		"addl	%%edi,%%ebx	\n\t"\
		"addl	%%edi,%%ecx	\n\t"\
		"addl	%%edi,%%edx	\n\t"\
		"\n\t"\
		"\n\t"/* SSE2_RADIX4_DIT_0TWIDDLE(add0, add1, add2, add3, r20): */\
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
		"\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm2		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"subpd	%%xmm6,%%xmm3		\n\t"\
		"movaps	%%xmm0,0x40(%%esi)	\n\t"\
		"movaps	%%xmm2,0x60(%%esi)	\n\t"\
		"movaps	%%xmm1,0x50(%%esi)	\n\t"\
		"movaps	%%xmm3,0x30(%%esi)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"addpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	%%xmm4,    (%%esi)	\n\t"\
		"movaps	%%xmm7,0x20(%%esi)	\n\t"\
		"movaps	%%xmm5,0x10(%%esi)	\n\t"\
		"movaps	%%xmm6,0x70(%%esi)	\n\t"\
		"\n\t"/* SSE2_RADIX4_DIT_0TWIDDLE_2NDOFTWO(add4, add5, add6, add7, r28): */\
		"addl	$0x80,%%esi			\n\t"\
		"movl	%[__arr_offsets],%%edi	\n\t"\
		"movl	0x50(%%edi),%%eax	\n\t"/* p14 */\
		"movl	0x54(%%edi),%%ebx	\n\t"/* p15 */\
		"movl	0x58(%%edi),%%ecx	\n\t"/* p16 */\
		"movl	0x5c(%%edi),%%edx	\n\t"/* p17 */\
		"movl	%[__add],%%edi		\n\t"\
		"addl	%%edi,%%eax	\n\t"\
		"addl	%%edi,%%ebx	\n\t"\
		"addl	%%edi,%%ecx	\n\t"\
		"addl	%%edi,%%edx	\n\t"\
		"\n\t"\
		"movaps	    (%%eax),%%xmm2	\n\t"\
		"movaps	    (%%ecx),%%xmm6	\n\t"\
		"movaps	0x10(%%eax),%%xmm3	\n\t"\
		"movaps	0x10(%%ecx),%%xmm7	\n\t"\
		"movaps	    (%%ebx),%%xmm0	\n\t"\
		"movaps	    (%%edx),%%xmm4	\n\t"\
		"movaps	0x10(%%ebx),%%xmm1	\n\t"\
		"movaps	0x10(%%edx),%%xmm5	\n\t"\
		"\n\t"\
		"movl	%[__isrt2],%%edi	\n\t"\
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
		"movaps	%%xmm0,0x40(%%esi)	\n\t"\
		"movaps	%%xmm1,0x50(%%esi)	\n\t"\
		"subpd	%%xmm7,%%xmm2		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"subpd	%%xmm6,%%xmm3		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"movaps	%%xmm4,    (%%esi)	\n\t"\
		"addpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	%%xmm5,0x10(%%esi)	\n\t"\
		"\n\t"\
		"movaps	(%%edi),%%xmm5	\n\t"/* ISRT2 */\
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
		"movaps	%%xmm3,0x30(%%esi)	\n\t"\
		"movaps	%%xmm6,0x70(%%esi)	\n\t"\
		"movaps	%%xmm0,0x20(%%esi)	\n\t"\
		"movaps	%%xmm1,0x60(%%esi)	\n\t"\
		"\n\t"/* SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS(r20,r22,r24,r26,r28,r2A,r2C,r2E): */\
		"movaps		-0x80(%%esi),%%xmm0	\n\t"/* r20 */\
		"movaps		-0x40(%%esi),%%xmm4	\n\t"/* r24 */\
		"movaps		-0x70(%%esi),%%xmm1	\n\t"/* r21 */\
		"movaps		-0x30(%%esi),%%xmm5	\n\t"/* r25 */\
		"movaps		     (%%esi),%%xmm2	\n\t"/* r28 */\
		"movaps		 0x50(%%esi),%%xmm7	\n\t"/* r2D */\
		"movaps		 0x10(%%esi),%%xmm3	\n\t"/* r29 */\
		"movaps		 0x40(%%esi),%%xmm6	\n\t"/* r2C */\
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
		"movaps		%%xmm0,     (%%esi)\n\t"\
		"movaps		%%xmm4, 0x40(%%esi)\n\t"\
		"movaps		%%xmm1, 0x10(%%esi)\n\t"\
		"movaps		%%xmm5,-0x30(%%esi)\n\t"\
		"movaps		%%xmm2,-0x80(%%esi)\n\t"\
		"movaps		%%xmm7,-0x40(%%esi)\n\t"\
		"movaps		%%xmm3,-0x70(%%esi)\n\t"\
		"movaps		%%xmm6, 0x50(%%esi)\n\t"\
		"\n\t"\
		"movaps		-0x60(%%esi),%%xmm0	\n\t"/* r22 */\
		"movaps		-0x20(%%esi),%%xmm4	\n\t"/* r26 */\
		"movaps		-0x50(%%esi),%%xmm1	\n\t"/* r23 */\
		"movaps		-0x10(%%esi),%%xmm5	\n\t"/* r27 */\
		"movaps		 0x20(%%esi),%%xmm2	\n\t"/* r2A */\
		"movaps		 0x70(%%esi),%%xmm7	\n\t"/* r2F */\
		"movaps		 0x30(%%esi),%%xmm3	\n\t"/* r2B */\
		"movaps		 0x60(%%esi),%%xmm6	\n\t"/* r2E */\
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
		"movaps		%%xmm0, 0x20(%%esi)\n\t"\
		"movaps		%%xmm4, 0x60(%%esi)\n\t"\
		"movaps		%%xmm1, 0x30(%%esi)\n\t"\
		"movaps		%%xmm5,-0x10(%%esi)\n\t"\
		"movaps		%%xmm2,-0x60(%%esi)\n\t"\
		"movaps		%%xmm7,-0x20(%%esi)\n\t"\
		"movaps		%%xmm3,-0x50(%%esi)\n\t"\
		"movaps		%%xmm6, 0x70(%%esi)\n\t"\
		"\n\t"\
	/*...Block 4: */\
		"addl	$0x80,%%esi			\n\t"\
		"movl	%[__arr_offsets],%%edi	\n\t"\
		"movl	0x60(%%edi),%%eax	\n\t"/* p18 */\
		"movl	0x64(%%edi),%%ebx	\n\t"/* p19 */\
		"movl	0x68(%%edi),%%ecx	\n\t"/* p1a */\
		"movl	0x6c(%%edi),%%edx	\n\t"/* p1b */\
		"movl	%[__add],%%edi		\n\t"\
		"addl	%%edi,%%eax	\n\t"\
		"addl	%%edi,%%ebx	\n\t"\
		"addl	%%edi,%%ecx	\n\t"\
		"addl	%%edi,%%edx	\n\t"\
		"\n\t"\
		"\n\t"/* SSE2_RADIX4_DIT_0TWIDDLE(add0, add1, add2, add3, r30): */\
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
		"\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm2		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"subpd	%%xmm6,%%xmm3		\n\t"\
		"movaps	%%xmm0,0x40(%%esi)	\n\t"\
		"movaps	%%xmm2,0x60(%%esi)	\n\t"\
		"movaps	%%xmm1,0x50(%%esi)	\n\t"\
		"movaps	%%xmm3,0x30(%%esi)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"addpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	%%xmm4,    (%%esi)	\n\t"\
		"movaps	%%xmm7,0x20(%%esi)	\n\t"\
		"movaps	%%xmm5,0x10(%%esi)	\n\t"\
		"movaps	%%xmm6,0x70(%%esi)	\n\t"\
		"\n\t"/* SSE2_RADIX4_DIT_0TWIDDLE_2NDOFTWO(add4, add5, add6, add7, r38): */\
		"addl	$0x80,%%esi			\n\t"\
		"movl	%[__arr_offsets],%%edi	\n\t"\
		"movl	0x70(%%edi),%%eax	\n\t"/* p1c */\
		"movl	0x74(%%edi),%%ebx	\n\t"/* p1d */\
		"movl	0x78(%%edi),%%ecx	\n\t"/* p1e */\
		"movl	0x7c(%%edi),%%edx	\n\t"/* p1f */\
		"movl	%[__add],%%edi		\n\t"\
		"addl	%%edi,%%eax	\n\t"\
		"addl	%%edi,%%ebx	\n\t"\
		"addl	%%edi,%%ecx	\n\t"\
		"addl	%%edi,%%edx	\n\t"\
		"\n\t"\
		"movaps	    (%%eax),%%xmm2	\n\t"\
		"movaps	    (%%ecx),%%xmm6	\n\t"\
		"movaps	0x10(%%eax),%%xmm3	\n\t"\
		"movaps	0x10(%%ecx),%%xmm7	\n\t"\
		"movaps	    (%%ebx),%%xmm0	\n\t"\
		"movaps	    (%%edx),%%xmm4	\n\t"\
		"movaps	0x10(%%ebx),%%xmm1	\n\t"\
		"movaps	0x10(%%edx),%%xmm5	\n\t"\
		"\n\t"\
		"movl	%[__isrt2],%%edi	\n\t"\
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
		"movaps	%%xmm0,0x40(%%esi)	\n\t"\
		"movaps	%%xmm1,0x50(%%esi)	\n\t"\
		"subpd	%%xmm7,%%xmm2		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"subpd	%%xmm6,%%xmm3		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"movaps	%%xmm4,    (%%esi)	\n\t"\
		"addpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	%%xmm5,0x10(%%esi)	\n\t"\
		"\n\t"\
		"movaps	(%%edi),%%xmm5	\n\t"/* ISRT2 */\
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
		"movaps	%%xmm3,0x30(%%esi)	\n\t"\
		"movaps	%%xmm6,0x70(%%esi)	\n\t"\
		"movaps	%%xmm0,0x20(%%esi)	\n\t"\
		"movaps	%%xmm1,0x60(%%esi)	\n\t"\
		"\n\t"/* SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS(r30,r32,r34,r36,r38,r3A,r3C,r3E): */\
		"movaps		-0x80(%%esi),%%xmm0	\n\t"/* r30 */\
		"movaps		-0x40(%%esi),%%xmm4	\n\t"/* r34 */\
		"movaps		-0x70(%%esi),%%xmm1	\n\t"/* r31 */\
		"movaps		-0x30(%%esi),%%xmm5	\n\t"/* r35 */\
		"movaps		     (%%esi),%%xmm2	\n\t"/* r38 */\
		"movaps		 0x50(%%esi),%%xmm7	\n\t"/* r3D */\
		"movaps		 0x10(%%esi),%%xmm3	\n\t"/* r39 */\
		"movaps		 0x40(%%esi),%%xmm6	\n\t"/* r3C */\
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
		"movaps		%%xmm0,     (%%esi)\n\t"\
		"movaps		%%xmm4, 0x40(%%esi)\n\t"\
		"movaps		%%xmm1, 0x10(%%esi)\n\t"\
		"movaps		%%xmm5,-0x30(%%esi)\n\t"\
		"movaps		%%xmm2,-0x80(%%esi)\n\t"\
		"movaps		%%xmm7,-0x40(%%esi)\n\t"\
		"movaps		%%xmm3,-0x70(%%esi)\n\t"\
		"movaps		%%xmm6, 0x50(%%esi)\n\t"\
		"\n\t"\
		"movaps		-0x60(%%esi),%%xmm0	\n\t"/* r32 */\
		"movaps		-0x20(%%esi),%%xmm4	\n\t"/* r36 */\
		"movaps		-0x50(%%esi),%%xmm1	\n\t"/* r33 */\
		"movaps		-0x10(%%esi),%%xmm5	\n\t"/* r37 */\
		"movaps		 0x20(%%esi),%%xmm2	\n\t"/* r3A */\
		"movaps		 0x70(%%esi),%%xmm7	\n\t"/* r3F */\
		"movaps		 0x30(%%esi),%%xmm3	\n\t"/* r3B */\
		"movaps		 0x60(%%esi),%%xmm6	\n\t"/* r3E */\
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
		"movaps		%%xmm0, 0x20(%%esi)\n\t"\
		"movaps		%%xmm4, 0x60(%%esi)\n\t"\
		"movaps		%%xmm1, 0x30(%%esi)\n\t"\
		"movaps		%%xmm5,-0x10(%%esi)\n\t"\
		"movaps		%%xmm2,-0x60(%%esi)\n\t"\
		"movaps		%%xmm7,-0x20(%%esi)\n\t"\
		"movaps		%%xmm3,-0x50(%%esi)\n\t"\
		"movaps		%%xmm6, 0x70(%%esi)\n\t"\
		"\n\t"\
	"\n\t"/*...and now do eight radix-4 transforms, including the internal twiddle factors:	*/\
		"\n\t"/*...Block 1: r00,r10,r20,r30	*/\
		"\n\t"\
		"movl	%[__isrt2],%%edi	\n\t"\
		"leal	0x10(%%edi),%%esi	\n\t/* cc0 */"\
		"movl	%[__r00],%%eax	\n\t"\
		"movl	%%eax,%%ebx		\n\t"\
		"movl	%%eax,%%ecx		\n\t"\
		"movl	%%eax,%%edx		\n\t"\
		"addl	$0x100,%%ebx	\n\t"\
		"addl	$0x200,%%ecx	\n\t"\
		"addl	$0x300,%%edx	\n\t"\
		"\n\t"\
		"movaps	    (%%eax),%%xmm0		\n\t"\
		"movaps	0x10(%%eax),%%xmm1		\n\t"\
		"movaps	    (%%ebx),%%xmm2		\n\t"\
		"movaps	0x10(%%ebx),%%xmm3		\n\t"\
		"subpd	    (%%ebx),%%xmm0		\n\t"\
		"subpd	0x10(%%ebx),%%xmm1		\n\t"\
		"addpd	    (%%eax),%%xmm2		\n\t"\
		"addpd	0x10(%%eax),%%xmm3		\n\t"\
		"movaps	    (%%ecx),%%xmm4		\n\t"\
		"movaps	0x10(%%ecx),%%xmm5		\n\t"\
		"movaps	    (%%edx),%%xmm6		\n\t"\
		"movaps	0x10(%%edx),%%xmm7		\n\t"\
		"subpd	    (%%edx),%%xmm4		\n\t"\
		"subpd	0x10(%%edx),%%xmm5		\n\t"\
		"addpd	    (%%ecx),%%xmm6		\n\t"\
		"addpd	0x10(%%ecx),%%xmm7		\n\t"\
		"\n\t"\
		"\n\t"/* SSE2_RADIX4_DIT_LASTPAIR_IN_PLACE_STRIDE(xmm2,xmm3,xmm0,xmm1,xmm6,xmm7,xmm4,xmm5): swap xmm[01]<->[23], xmm[45]<->[67] */\
		"addpd	%%xmm6,%%xmm2		\n\t"\
		"addpd	%%xmm7,%%xmm3		\n\t"\
		"movaps	%%xmm2,    (%%eax)	\n\t"\
		"movaps	%%xmm3,0x10(%%eax)	\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"subpd	%%xmm6,%%xmm2		\n\t"\
		"subpd	%%xmm7,%%xmm3		\n\t"\
		"movaps	%%xmm2,    (%%ecx)	\n\t"\
		"movaps	%%xmm3,0x10(%%ecx)	\n\t"\
		"\n\t"\
		"addpd	%%xmm5,%%xmm0		\n\t"\
		"subpd	%%xmm4,%%xmm1		\n\t"\
		"movaps	%%xmm0,    (%%ebx)	\n\t"\
		"movaps	%%xmm1,0x10(%%ebx)	\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"subpd	%%xmm5,%%xmm0		\n\t"\
		"addpd	%%xmm4,%%xmm1		\n\t"\
		"movaps	%%xmm0,    (%%edx)	\n\t"\
		"movaps	%%xmm1,0x10(%%edx)	\n\t"\
		"\n\t"\
		"\n\t"/*...Block 5: r08,r18,r28,r38	*/\
		"\n\t"\
		"addl	$0x80,%%eax		\n\t"/* r08 */\
		"addl	$0x80,%%ebx		\n\t"\
		"addl	$0x80,%%ecx		\n\t"\
		"addl	$0x80,%%edx		\n\t"\
		"movaps	(%%edi),%%xmm2	\n\t"/* isrt2 */\
		"\n\t"\
		"movaps	    (%%ecx),%%xmm4	\n\t"\
		"movaps	0x10(%%ecx),%%xmm5	\n\t"\
		"movaps	    (%%edx),%%xmm0	\n\t"\
		"movaps	0x10(%%edx),%%xmm1	\n\t"\
		"addpd	0x10(%%ecx),%%xmm4	\n\t"\
		"subpd	    (%%ecx),%%xmm5	\n\t"\
		"subpd	0x10(%%edx),%%xmm0	\n\t"\
		"addpd	    (%%edx),%%xmm1	\n\t"\
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
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"movaps	    (%%ebx),%%xmm2	\n\t"\
		"movaps	0x10(%%ebx),%%xmm3	\n\t"\
		"subpd	0x10(%%ebx),%%xmm0	\n\t"\
		"subpd	    (%%ebx),%%xmm1	\n\t"\
		"addpd	    (%%eax),%%xmm3	\n\t"\
		"addpd	0x10(%%eax),%%xmm2	\n\t"\
		"\n\t"\
		"\n\t"/* SSE2_RADIX4_DIT_LASTPAIR_IN_PLACE_STRIDE(xmm3,xmm1,xmm0,xmm2,xmm4,xmm5,xmm6,xmm7): swap xmm0123<->3102 */\
		"addpd	%%xmm4,%%xmm3		\n\t"\
		"addpd	%%xmm5,%%xmm1		\n\t"\
		"movaps	%%xmm3,    (%%eax)	\n\t"\
		"movaps	%%xmm1,0x10(%%eax)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"subpd	%%xmm4,%%xmm3		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"movaps	%%xmm3,    (%%ecx)	\n\t"\
		"movaps	%%xmm1,0x10(%%ecx)	\n\t"\
		"\n\t"\
		"addpd	%%xmm7,%%xmm0		\n\t"\
		"subpd	%%xmm6,%%xmm2		\n\t"\
		"movaps	%%xmm0,    (%%ebx)	\n\t"\
		"movaps	%%xmm2,0x10(%%ebx)	\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"subpd	%%xmm7,%%xmm0		\n\t"\
		"addpd	%%xmm6,%%xmm2		\n\t"\
		"movaps	%%xmm0,    (%%edx)	\n\t"\
		"movaps	%%xmm2,0x10(%%edx)	\n\t"\
		"\n\t"\
		"\n\t"/*...Block 3: r04,r14,r24,r34	*/\
		"\n\t"\
		"subl	$0x40,%%eax 	\n\t"/* r04 */\
		"subl	$0x40,%%ebx		\n\t"\
		"subl	$0x40,%%ecx		\n\t"\
		"subl	$0x40,%%edx		\n\t"\
		"\n\t"\
		"movaps	    (%%ecx),%%xmm4	\n\t"\
		"movaps	    (%%edx),%%xmm0	\n\t"\
		"movaps	0x10(%%ecx),%%xmm5	\n\t"\
		"movaps	0x10(%%edx),%%xmm1	\n\t"\
		"movaps	    (%%ecx),%%xmm6	\n\t"\
		"movaps	    (%%edx),%%xmm2	\n\t"\
		"movaps	0x10(%%ecx),%%xmm7	\n\t"\
		"movaps	0x10(%%edx),%%xmm3	\n\t"\
		"\n\t"\
		"mulpd	    (%%esi),%%xmm4	\n\t"\
		"mulpd	0x10(%%esi),%%xmm0	\n\t"\
		"mulpd	    (%%esi),%%xmm5	\n\t"\
		"mulpd	0x10(%%esi),%%xmm1	\n\t"\
		"mulpd	0x10(%%esi),%%xmm6	\n\t"\
		"mulpd	    (%%esi),%%xmm2	\n\t"\
		"mulpd	0x10(%%esi),%%xmm7	\n\t"\
		"mulpd	    (%%esi),%%xmm3	\n\t"\
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
		"movaps	    (%%ebx),%%xmm2	\n\t"\
		"movaps	0x10(%%ebx),%%xmm3	\n\t"\
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"addpd	0x10(%%ebx),%%xmm2	\n\t"\
		"subpd	    (%%ebx),%%xmm3	\n\t"\
		"mulpd	    (%%edi),%%xmm2		\n\t"\
		"mulpd	    (%%edi),%%xmm3		\n\t"\
		"\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"\n\t"\
		"\n\t"/* SSE2_RADIX4_DIT_LASTPAIR_IN_PLACE_STRIDE(xmm2,xmm3,xmm0,xmm1,xmm4,xmm5,xmm6,xmm7): swap xmm[01]<->[23] */\
		"addpd	%%xmm4,%%xmm2		\n\t"\
		"addpd	%%xmm5,%%xmm3		\n\t"\
		"movaps	%%xmm2,    (%%eax)	\n\t"\
		"movaps	%%xmm3,0x10(%%eax)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"subpd	%%xmm4,%%xmm2		\n\t"\
		"subpd	%%xmm5,%%xmm3		\n\t"\
		"movaps	%%xmm2,    (%%ecx)	\n\t"\
		"movaps	%%xmm3,0x10(%%ecx)	\n\t"\
		"\n\t"\
		"addpd	%%xmm7,%%xmm0		\n\t"\
		"subpd	%%xmm6,%%xmm1		\n\t"\
		"movaps	%%xmm0,    (%%ebx)	\n\t"\
		"movaps	%%xmm1,0x10(%%ebx)	\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"subpd	%%xmm7,%%xmm0		\n\t"\
		"addpd	%%xmm6,%%xmm1		\n\t"\
		"movaps	%%xmm0,    (%%edx)	\n\t"\
		"movaps	%%xmm1,0x10(%%edx)	\n\t"\
		"\n\t"\
		"\n\t"/*...Block 7: r0C,r1C,r2C,r3C	*/\
		"\n\t"\
		"addl	$0x80,%%eax 	\n\t"/* r0C */\
		"addl	$0x80,%%ebx		\n\t"\
		"addl	$0x80,%%ecx		\n\t"\
		"addl	$0x80,%%edx		\n\t"\
		"\n\t"\
		"movaps	    (%%ecx),%%xmm4	\n\t"\
		"movaps	    (%%edx),%%xmm0	\n\t"\
		"movaps	0x10(%%ecx),%%xmm5	\n\t"\
		"movaps	0x10(%%edx),%%xmm1	\n\t"\
		"movaps	    (%%ecx),%%xmm6	\n\t"\
		"movaps	    (%%edx),%%xmm2	\n\t"\
		"movaps	0x10(%%ecx),%%xmm7	\n\t"\
		"movaps	0x10(%%edx),%%xmm3	\n\t"\
		"\n\t"\
		"mulpd	0x10(%%esi),%%xmm4	\n\t"\
		"mulpd	    (%%esi),%%xmm0	\n\t"\
		"mulpd	0x10(%%esi),%%xmm5	\n\t"\
		"mulpd	    (%%esi),%%xmm1	\n\t"\
		"mulpd	    (%%esi),%%xmm6	\n\t"\
		"mulpd	0x10(%%esi),%%xmm2	\n\t"\
		"mulpd	    (%%esi),%%xmm7	\n\t"\
		"mulpd	0x10(%%esi),%%xmm3	\n\t"\
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
		"movaps	    (%%ebx),%%xmm2	\n\t"\
		"movaps	0x10(%%ebx),%%xmm3	\n\t"\
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"subpd	0x10(%%ebx),%%xmm2	\n\t"\
		"addpd	    (%%ebx),%%xmm3	\n\t"\
		"mulpd	    (%%edi),%%xmm2		\n\t"\
		"mulpd	    (%%edi),%%xmm3		\n\t"\
		"\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"\n\t"\
		"\n\t"/* SSE2_RADIX4_DIT_LASTPAIR_IN_PLACE_STRIDE(xmm0,xmm1,xmm2,xmm3,xmm6,xmm7,xmm4,xmm5): swap xmm[45]<->[67] */\
		"addpd	%%xmm6,%%xmm0		\n\t"\
		"addpd	%%xmm7,%%xmm1		\n\t"\
		"movaps	%%xmm0,    (%%eax)	\n\t"\
		"movaps	%%xmm1,0x10(%%eax)	\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"subpd	%%xmm6,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm1		\n\t"\
		"movaps	%%xmm0,    (%%ecx)	\n\t"\
		"movaps	%%xmm1,0x10(%%ecx)	\n\t"\
		"\n\t"\
		"addpd	%%xmm5,%%xmm2		\n\t"\
		"subpd	%%xmm4,%%xmm3		\n\t"\
		"movaps	%%xmm2,    (%%ebx)	\n\t"\
		"movaps	%%xmm3,0x10(%%ebx)	\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"subpd	%%xmm5,%%xmm2		\n\t"\
		"addpd	%%xmm4,%%xmm3		\n\t"\
		"movaps	%%xmm2,    (%%edx)	\n\t"\
		"movaps	%%xmm3,0x10(%%edx)	\n\t"\
		"\n\t"\
		"\n\t"/*...Block 2: r02,r12,r22,r32	*/\
		"\n\t"\
		"subl	$0xa0,%%eax 	\n\t"/* r02 */\
		"subl	$0xa0,%%ebx			\n\t"\
		"subl	$0xa0,%%ecx			\n\t"\
		"subl	$0xa0,%%edx			\n\t"\
		"addl	$0x30,%%edi \n\t"/* cc1 */\
		"addl	$0x40,%%esi \n\t"/* cc3 */\
		"\n\t"\
		"movaps	    (%%ecx),%%xmm4	\n\t"\
		"movaps	    (%%edx),%%xmm0	\n\t"\
		"movaps	0x10(%%ecx),%%xmm5	\n\t"\
		"movaps	0x10(%%edx),%%xmm1	\n\t"\
		"movaps	    (%%ecx),%%xmm6	\n\t"\
		"movaps	    (%%edx),%%xmm2	\n\t"\
		"movaps	0x10(%%ecx),%%xmm7	\n\t"\
		"movaps	0x10(%%edx),%%xmm3	\n\t"\
		"\n\t"\
		"mulpd	    (%%edi),%%xmm4	\n\t"\
		"mulpd	    (%%esi),%%xmm0	\n\t"\
		"mulpd	    (%%edi),%%xmm5	\n\t"\
		"mulpd	    (%%esi),%%xmm1	\n\t"\
		"mulpd	0x10(%%edi),%%xmm6	\n\t"\
		"mulpd	0x10(%%esi),%%xmm2	\n\t"\
		"mulpd	0x10(%%edi),%%xmm7	\n\t"\
		"mulpd	0x10(%%esi),%%xmm3	\n\t"\
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
		"subl	$0x40,%%esi \n\t"/* cc0 */\
		"movaps	    (%%ebx),%%xmm2	\n\t"\
		"movaps	0x10(%%ebx),%%xmm3	\n\t"\
		"movaps	    (%%ebx),%%xmm0	\n\t"\
		"movaps	0x10(%%ebx),%%xmm1	\n\t"\
		"mulpd	    (%%esi),%%xmm2	\n\t"\
		"mulpd	0x10(%%esi),%%xmm1	\n\t"\
		"mulpd	    (%%esi),%%xmm3	\n\t"\
		"mulpd	0x10(%%esi),%%xmm0	\n\t"\
		"addpd	%%xmm1,%%xmm2		\n\t"\
		"subpd	%%xmm0,%%xmm3		\n\t"\
		"\n\t"\
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"\n\t"\
		"\n\t"/* SSE2_RADIX4_DIT_LASTPAIR_IN_PLACE_STRIDE(xmm2,xmm3,xmm0,xmm1,xmm6,xmm7,xmm4,xmm5): swap xmm[01]<->[23], xmm[45]<->[67] */\
		"addpd	%%xmm6,%%xmm2		\n\t"\
		"addpd	%%xmm7,%%xmm3		\n\t"\
		"movaps	%%xmm2,    (%%eax)	\n\t"\
		"movaps	%%xmm3,0x10(%%eax)	\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"subpd	%%xmm6,%%xmm2		\n\t"\
		"subpd	%%xmm7,%%xmm3		\n\t"\
		"movaps	%%xmm2,    (%%ecx)	\n\t"\
		"movaps	%%xmm3,0x10(%%ecx)	\n\t"\
		"\n\t"\
		"addpd	%%xmm5,%%xmm0		\n\t"\
		"subpd	%%xmm4,%%xmm1		\n\t"\
		"movaps	%%xmm0,    (%%ebx)	\n\t"\
		"movaps	%%xmm1,0x10(%%ebx)	\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"subpd	%%xmm5,%%xmm0		\n\t"\
		"addpd	%%xmm4,%%xmm1		\n\t"\
		"movaps	%%xmm0,    (%%edx)	\n\t"\
		"movaps	%%xmm1,0x10(%%edx)	\n\t"\
		"\n\t"\
		"\n\t"/*...Block 6: r0A,r1A,r2A,r3A	*/\
		"\n\t"\
		"addl	$0x80,%%eax 	\n\t"/* r0A */\
		"addl	$0x80,%%ebx			\n\t"\
		"addl	$0x80,%%ecx			\n\t"\
		"addl	$0x80,%%edx			\n\t"\
		"addl	$0x40,%%esi \n\t"/* cc3 */\
		"\n\t"\
		"movaps	    (%%ecx),%%xmm4	\n\t"\
		"movaps	    (%%edx),%%xmm0	\n\t"\
		"movaps	0x10(%%ecx),%%xmm5	\n\t"\
		"movaps	0x10(%%edx),%%xmm1	\n\t"\
		"movaps	    (%%ecx),%%xmm6	\n\t"\
		"movaps	    (%%edx),%%xmm2	\n\t"\
		"movaps	0x10(%%ecx),%%xmm7	\n\t"\
		"movaps	0x10(%%edx),%%xmm3	\n\t"\
		"\n\t"\
		"mulpd	0x10(%%esi),%%xmm4	\n\t"\
		"mulpd	    (%%edi),%%xmm0	\n\t"\
		"mulpd	0x10(%%esi),%%xmm5	\n\t"\
		"mulpd	    (%%edi),%%xmm1	\n\t"\
		"mulpd	    (%%esi),%%xmm6	\n\t"\
		"mulpd	0x10(%%edi),%%xmm2	\n\t"\
		"mulpd	    (%%esi),%%xmm7	\n\t"\
		"mulpd	0x10(%%edi),%%xmm3	\n\t"\
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
		"subl	$0x40,%%esi \n\t"/* cc0 */\
		"movaps	    (%%ebx),%%xmm2	\n\t"\
		"movaps	0x10(%%ebx),%%xmm3	\n\t"\
		"movaps	    (%%ebx),%%xmm0	\n\t"\
		"movaps	0x10(%%ebx),%%xmm1	\n\t"\
		"mulpd	0x10(%%esi),%%xmm2	\n\t"\
		"mulpd	    (%%esi),%%xmm1	\n\t"\
		"mulpd	0x10(%%esi),%%xmm3	\n\t"\
		"mulpd	    (%%esi),%%xmm0	\n\t"\
		"subpd	%%xmm1,%%xmm2		\n\t"\
		"addpd	%%xmm0,%%xmm3		\n\t"\
		"\n\t"\
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"\n\t"\
		"\n\t"/* SSE2_RADIX4_DIT_LASTPAIR_IN_PLACE_STRIDE(xmm0,xmm1,xmm2,xmm3,xmm6,xmm7,xmm4,xmm5): swap xmm[45]<->[67] */\
		"addpd	%%xmm6,%%xmm0		\n\t"\
		"addpd	%%xmm7,%%xmm1		\n\t"\
		"movaps	%%xmm0,    (%%eax)	\n\t"\
		"movaps	%%xmm1,0x10(%%eax)	\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"subpd	%%xmm6,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm1		\n\t"\
		"movaps	%%xmm0,    (%%ecx)	\n\t"\
		"movaps	%%xmm1,0x10(%%ecx)	\n\t"\
		"\n\t"\
		"addpd	%%xmm5,%%xmm2		\n\t"\
		"subpd	%%xmm4,%%xmm3		\n\t"\
		"movaps	%%xmm2,    (%%ebx)	\n\t"\
		"movaps	%%xmm3,0x10(%%ebx)	\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"subpd	%%xmm5,%%xmm2		\n\t"\
		"addpd	%%xmm4,%%xmm3		\n\t"\
		"movaps	%%xmm2,    (%%edx)	\n\t"\
		"movaps	%%xmm3,0x10(%%edx)	\n\t"\
		"\n\t"\
		"\n\t"/*...Block 4: r06,r16,r26,r36	*/\
		"\n\t"\
		"subl	$0x40,%%eax 	\n\t"/* r06 */\
		"subl	$0x40,%%ebx			\n\t"\
		"subl	$0x40,%%ecx			\n\t"\
		"subl	$0x40,%%edx			\n\t"\
		"addl	$0x40,%%esi \n\t"/* cc3 */\
		"\n\t"\
		"movaps	    (%%ecx),%%xmm4	\n\t"\
		"movaps	    (%%edx),%%xmm0	\n\t"\
		"movaps	0x10(%%ecx),%%xmm5	\n\t"\
		"movaps	0x10(%%edx),%%xmm1	\n\t"\
		"movaps	    (%%ecx),%%xmm6	\n\t"\
		"movaps	    (%%edx),%%xmm2	\n\t"\
		"movaps	0x10(%%ecx),%%xmm7	\n\t"\
		"movaps	0x10(%%edx),%%xmm3	\n\t"\
		"\n\t"\
		"mulpd	    (%%esi),%%xmm4	\n\t"\
		"mulpd	0x10(%%edi),%%xmm0	\n\t"\
		"mulpd	    (%%esi),%%xmm5	\n\t"\
		"mulpd	0x10(%%edi),%%xmm1	\n\t"\
		"mulpd	0x10(%%esi),%%xmm6	\n\t"\
		"mulpd	    (%%edi),%%xmm2	\n\t"\
		"mulpd	0x10(%%esi),%%xmm7	\n\t"\
		"mulpd	    (%%edi),%%xmm3	\n\t"\
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
		"subl	$0x40,%%esi \n\t"/* cc0 */\
		"movaps	    (%%ebx),%%xmm2	\n\t"\
		"movaps	0x10(%%ebx),%%xmm3	\n\t"\
		"movaps	    (%%ebx),%%xmm0	\n\t"\
		"movaps	0x10(%%ebx),%%xmm1	\n\t"\
		"mulpd	0x10(%%esi),%%xmm2	\n\t"\
		"mulpd	    (%%esi),%%xmm1	\n\t"\
		"mulpd	0x10(%%esi),%%xmm3	\n\t"\
		"mulpd	    (%%esi),%%xmm0	\n\t"\
		"addpd	%%xmm1,%%xmm2		\n\t"\
		"subpd	%%xmm0,%%xmm3		\n\t"\
		"\n\t"\
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"\n\t"\
		"\n\t"/* SSE2_RADIX4_DIT_LASTPAIR_IN_PLACE_STRIDE(xmm2,xmm3,xmm0,xmm1,xmm6,xmm7,xmm4,xmm5): swap xmm[01]<->[23], xmm[45]<->[67] */\
		"addpd	%%xmm6,%%xmm2		\n\t"\
		"addpd	%%xmm7,%%xmm3		\n\t"\
		"movaps	%%xmm2,    (%%eax)	\n\t"\
		"movaps	%%xmm3,0x10(%%eax)	\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"subpd	%%xmm6,%%xmm2		\n\t"\
		"subpd	%%xmm7,%%xmm3		\n\t"\
		"movaps	%%xmm2,    (%%ecx)	\n\t"\
		"movaps	%%xmm3,0x10(%%ecx)	\n\t"\
		"\n\t"\
		"addpd	%%xmm5,%%xmm0		\n\t"\
		"subpd	%%xmm4,%%xmm1		\n\t"\
		"movaps	%%xmm0,    (%%ebx)	\n\t"\
		"movaps	%%xmm1,0x10(%%ebx)	\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"subpd	%%xmm5,%%xmm0		\n\t"\
		"addpd	%%xmm4,%%xmm1		\n\t"\
		"movaps	%%xmm0,    (%%edx)	\n\t"\
		"movaps	%%xmm1,0x10(%%edx)	\n\t"\
		"\n\t"\
		"\n\t"/*...Block 8: r0E,r1E,r2E,r3E	*/\
		"\n\t"\
		"addl	$0x80,%%eax 	\n\t"/* r0E */\
		"addl	$0x80,%%ebx		\n\t"\
		"addl	$0x80,%%ecx		\n\t"\
		"addl	$0x80,%%edx		\n\t"\
		"addl	$0x40,%%esi \n\t"/* cc3 */\
		"\n\t"\
		"movaps	    (%%ecx),%%xmm4	\n\t"\
		"movaps	    (%%edx),%%xmm0	\n\t"\
		"movaps	0x10(%%ecx),%%xmm5	\n\t"\
		"movaps	0x10(%%edx),%%xmm1	\n\t"\
		"movaps	    (%%ecx),%%xmm6	\n\t"\
		"movaps	    (%%edx),%%xmm2	\n\t"\
		"movaps	0x10(%%ecx),%%xmm7	\n\t"\
		"movaps	0x10(%%edx),%%xmm3	\n\t"\
		"\n\t"\
		"mulpd	0x10(%%edi),%%xmm4	\n\t"\
		"mulpd	0x10(%%esi),%%xmm0	\n\t"\
		"mulpd	0x10(%%edi),%%xmm5	\n\t"\
		"mulpd	0x10(%%esi),%%xmm1	\n\t"\
		"mulpd	    (%%edi),%%xmm6	\n\t"\
		"mulpd	    (%%esi),%%xmm2	\n\t"\
		"mulpd	    (%%edi),%%xmm7	\n\t"\
		"mulpd	    (%%esi),%%xmm3	\n\t"\
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
		"subl	$0x40,%%esi \n\t"/* cc0 */\
		"movaps	    (%%ebx),%%xmm2	\n\t"\
		"movaps	0x10(%%ebx),%%xmm3	\n\t"\
		"movaps	    (%%ebx),%%xmm0	\n\t"\
		"movaps	0x10(%%ebx),%%xmm1	\n\t"\
		"mulpd	    (%%esi),%%xmm2	\n\t"\
		"mulpd	0x10(%%esi),%%xmm1	\n\t"\
		"mulpd	    (%%esi),%%xmm3	\n\t"\
		"mulpd	0x10(%%esi),%%xmm0	\n\t"\
		"subpd	%%xmm1,%%xmm2		\n\t"\
		"addpd	%%xmm0,%%xmm3		\n\t"\
		"\n\t"\
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"\n\t"\
		"\n\t"/* SSE2_RADIX4_DIT_LASTPAIR_IN_PLACE_STRIDE(xmm0,xmm1,xmm2,xmm3,xmm6,xmm7,xmm4,xmm5): swap xmm[45]<->[67] */\
		"addpd	%%xmm6,%%xmm0		\n\t"\
		"addpd	%%xmm7,%%xmm1		\n\t"\
		"movaps	%%xmm0,    (%%eax)	\n\t"\
		"movaps	%%xmm1,0x10(%%eax)	\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"subpd	%%xmm6,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm1		\n\t"\
		"movaps	%%xmm0,    (%%ecx)	\n\t"\
		"movaps	%%xmm1,0x10(%%ecx)	\n\t"\
		"\n\t"\
		"addpd	%%xmm5,%%xmm2		\n\t"\
		"subpd	%%xmm4,%%xmm3		\n\t"\
		"movaps	%%xmm2,    (%%ebx)	\n\t"\
		"movaps	%%xmm3,0x10(%%ebx)	\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"subpd	%%xmm5,%%xmm2		\n\t"\
		"addpd	%%xmm4,%%xmm3		\n\t"\
		"movaps	%%xmm2,    (%%edx)	\n\t"\
		"movaps	%%xmm3,0x10(%%edx)	\n\t"\
	"popl %%ebx	\n\t"\
		:					/* outputs: none */\
		: [__add] "m" (Xadd)	/* All inputs from memory addresses here */\
		 ,[__arr_offsets] "m" (Xarr_offsets)\
		 ,[__r00] "m" (Xr00)\
		 ,[__isrt2] "m" (Xisrt2)\
		: "cc","memory","eax",/*"ebx",*/"ecx","edx","edi","esi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"	/* Clobbered registers */\
	);\
	}

	#define SSE2_RADIX32_DIF_NOTWIDDLE(Xadd,Xarr_offsets, Xr00, Xisrt2)\
	{\
	__asm__ volatile (\
	"pushl %%ebx	\n\t"/* Explicit save/restore of PIC register */\
		"\n\t"/* SSE2_RADIX4_DIF_IN_PLACE(r00,r20,r10,r30): */\
		"movl	%[__isrt2],%%edi	\n\t"\
		"movl	%[__r00],%%eax	\n\t"\
		"movl	$0x200,%%ebx	\n\t"\
		"movl	$0x100,%%ecx	\n\t"\
		"movl	$0x300,%%edx	\n\t"\
		"addl	%%eax,%%ebx		\n\t"\
		"addl	%%eax,%%ecx		\n\t"\
		"addl	%%eax,%%edx		\n\t"\
		"\n\t"\
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"movaps	    (%%eax),%%xmm2	\n\t"\
		"movaps	0x10(%%eax),%%xmm3	\n\t"\
		"addpd	    (%%ebx),%%xmm0	\n\t"\
		"addpd	0x10(%%ebx),%%xmm1	\n\t"\
		"subpd	    (%%ebx),%%xmm2	\n\t"\
		"subpd	0x10(%%ebx),%%xmm3	\n\t"\
		"\n\t"\
		"movaps	    (%%ecx),%%xmm4	\n\t"\
		"movaps	0x10(%%ecx),%%xmm5	\n\t"\
		"movaps	    (%%ecx),%%xmm6	\n\t"\
		"movaps	0x10(%%ecx),%%xmm7	\n\t"\
		"addpd	    (%%edx),%%xmm4	\n\t"\
		"addpd	0x10(%%edx),%%xmm5	\n\t"\
		"subpd	    (%%edx),%%xmm6	\n\t"\
		"subpd	0x10(%%edx),%%xmm7	\n\t"\
		"\n\t"/* Finish radix-4 butterfly and store results into temporary-array slots: */\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"movaps	%%xmm0,    (%%ebx)	\n\t"\
		"movaps	%%xmm1,0x10(%%ebx)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"movaps	%%xmm4,    (%%eax)	\n\t"\
		"movaps	%%xmm5,0x10(%%eax)	\n\t"\
		"\n\t"\
		"subpd	%%xmm7,%%xmm2		\n\t"\
		"subpd	%%xmm6,%%xmm3		\n\t"\
		"movaps	%%xmm2,    (%%ecx)	\n\t"\
		"movaps	%%xmm3,0x10(%%edx)	\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	%%xmm7,    (%%edx)	\n\t"\
		"movaps	%%xmm6,0x10(%%ecx)	\n\t"\
		"\n\t"\
		"\n\t"/* SSE2_RADIX4_DIF_IN_PLACE_2NDOFTWO(r08,r28,r18,r38): */\
		"addl	$0x80,%%eax		\n\t"\
		"addl	$0x80,%%ebx		\n\t"\
		"addl	$0x80,%%ecx		\n\t"\
		"addl	$0x80,%%edx		\n\t"\
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"movaps	    (%%eax),%%xmm2	\n\t"\
		"movaps	0x10(%%eax),%%xmm3	\n\t"\
		"addpd	    (%%ebx),%%xmm0	\n\t"\
		"addpd	0x10(%%ebx),%%xmm1	\n\t"\
		"subpd	    (%%ebx),%%xmm2	\n\t"\
		"subpd	0x10(%%ebx),%%xmm3	\n\t"\
		"movaps	    (%%ecx),%%xmm4	\n\t"\
		"movaps	0x10(%%ecx),%%xmm5	\n\t"\
		"movaps	    (%%ecx),%%xmm6	\n\t"\
		"movaps	0x10(%%ecx),%%xmm7	\n\t"\
		"addpd	    (%%edx),%%xmm6	\n\t"\
		"addpd	0x10(%%edx),%%xmm7	\n\t"\
		"subpd	    (%%edx),%%xmm4	\n\t"\
		"subpd	0x10(%%edx),%%xmm5	\n\t"\
		"\n\t"/* Finish radix-4 butterfly and store results into temporary-array slots: */\
		"subpd	%%xmm6,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm1		\n\t"\
		"subpd	%%xmm5,%%xmm2		\n\t"\
		"movaps	%%xmm0,    (%%ebx)	\n\t"\
		"movaps	%%xmm1,0x10(%%ebx)	\n\t"\
		"subpd	%%xmm4,%%xmm3		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm0,%%xmm6		\n\t"\
		"addpd	%%xmm2,%%xmm5		\n\t"\
		"addpd	%%xmm1,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm4		\n\t"\
		"movaps	%%xmm6,    (%%eax)	\n\t"\
		"movaps	%%xmm7,0x10(%%eax)	\n\t"\
		"\n\t"\
		"movaps	(%%edi),%%xmm0	\n\t"/* isrt2 */\
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
		"movaps	%%xmm2,    (%%ecx)	\n\t"\
		"movaps	%%xmm5,    (%%edx)	\n\t"\
		"movaps	%%xmm6,0x10(%%ecx)	\n\t"\
		"movaps	%%xmm7,0x10(%%edx)	\n\t"\
		"\n\t"/* SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS(r00,r10,r20,r30,r08,r18,r28,r38): */\
		"subl	$0x80,%%eax	\n\t"/* r00 */\
		"subl	$0x200,%%ebx	\n\t"/* r08 */\
		"addl	$0x80,%%ecx	\n\t"/* r20 */\
		"subl	$0x100,%%edx	\n\t"/* r28 */\
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	    (%%ecx),%%xmm4	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"movaps	0x10(%%ecx),%%xmm5	\n\t"\
		"movaps	    (%%ebx),%%xmm2	\n\t"\
		"movaps	0x10(%%edx),%%xmm7	\n\t"\
		"movaps	0x10(%%ebx),%%xmm3	\n\t"\
		"movaps	    (%%edx),%%xmm6	\n\t"\
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
		"movaps	%%xmm0,    (%%ebx)	\n\t"\
		"movaps	%%xmm4,    (%%ecx)	\n\t"\
		"movaps	%%xmm1,0x10(%%ebx)	\n\t"\
		"movaps	%%xmm5,0x10(%%edx)	\n\t"\
		"movaps	%%xmm2,    (%%eax)	\n\t"\
		"movaps	%%xmm7,    (%%edx)	\n\t"\
		"movaps	%%xmm3,0x10(%%eax)	\n\t"\
		"movaps	%%xmm6,0x10(%%ecx)	\n\t"\
		"\n\t"\
		"addl	$0x100,%%eax	\n\t"/* r10 */\
		"addl	$0x100,%%ebx	\n\t"/* r18 */\
		"addl	$0x100,%%ecx	\n\t"/* r30 */\
		"addl	$0x100,%%edx	\n\t"/* r38 */\
		"\n\t"\
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	    (%%ecx),%%xmm4	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"movaps	0x10(%%ecx),%%xmm5	\n\t"\
		"movaps	    (%%ebx),%%xmm2	\n\t"\
		"movaps	0x10(%%edx),%%xmm7	\n\t"\
		"movaps	0x10(%%ebx),%%xmm3	\n\t"\
		"movaps	    (%%edx),%%xmm6	\n\t"\
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
		"movaps	%%xmm0,    (%%ebx)	\n\t"\
		"movaps	%%xmm4,    (%%ecx)	\n\t"\
		"movaps	%%xmm1,0x10(%%ebx)	\n\t"\
		"movaps	%%xmm5,0x10(%%edx)	\n\t"\
		"movaps	%%xmm2,    (%%eax)	\n\t"\
		"movaps	%%xmm7,    (%%edx)	\n\t"\
		"movaps	%%xmm3,0x10(%%eax)	\n\t"\
		"movaps	%%xmm6,0x10(%%ecx)	\n\t"\
		"\n\t"\
		"\n\t"/* SSE2_RADIX4_DIF_IN_PLACE(r04,r24,r14,r34): */\
		"\n\t"\
		"subl	$0xC0,%%eax	\n\t"/* r04 */\
		"addl	$0xC0,%%ebx	\n\t"/* r24 */\
		"subl	$0x1C0,%%ecx	\n\t"/* r14 */\
		"subl	$0x40,%%edx	\n\t"/* r34 */\
		"\n\t"\
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"movaps	    (%%eax),%%xmm2	\n\t"\
		"movaps	0x10(%%eax),%%xmm3	\n\t"\
		"addpd	    (%%ebx),%%xmm0	\n\t"\
		"addpd	0x10(%%ebx),%%xmm1	\n\t"\
		"subpd	    (%%ebx),%%xmm2	\n\t"\
		"subpd	0x10(%%ebx),%%xmm3	\n\t"\
		"\n\t"\
		"movaps	    (%%ecx),%%xmm4	\n\t"\
		"movaps	0x10(%%ecx),%%xmm5	\n\t"\
		"movaps	    (%%ecx),%%xmm6	\n\t"\
		"movaps	0x10(%%ecx),%%xmm7	\n\t"\
		"addpd	    (%%edx),%%xmm4	\n\t"\
		"addpd	0x10(%%edx),%%xmm5	\n\t"\
		"subpd	    (%%edx),%%xmm6	\n\t"\
		"subpd	0x10(%%edx),%%xmm7	\n\t"\
		"\n\t"/* Finish radix-4 butterfly and store results into temporary-array slots: */\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"movaps	%%xmm0,    (%%ebx)	\n\t"\
		"movaps	%%xmm1,0x10(%%ebx)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"movaps	%%xmm4,    (%%eax)	\n\t"\
		"movaps	%%xmm5,0x10(%%eax)	\n\t"\
		"\n\t"\
		"subpd	%%xmm7,%%xmm2		\n\t"\
		"subpd	%%xmm6,%%xmm3		\n\t"\
		"movaps	%%xmm2,    (%%ecx)	\n\t"\
		"movaps	%%xmm3,0x10(%%edx)	\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	%%xmm7,    (%%edx)	\n\t"\
		"movaps	%%xmm6,0x10(%%ecx)	\n\t"\
		"\n\t"\
		"\n\t"/* SSE2_RADIX4_DIF_IN_PLACE_2NDOFTWO(r0C,r2C,r1C,r3C): */\
		"addl	$0x80,%%eax		\n\t"\
		"addl	$0x80,%%ebx		\n\t"\
		"addl	$0x80,%%ecx		\n\t"\
		"addl	$0x80,%%edx		\n\t"\
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"movaps	    (%%eax),%%xmm2	\n\t"\
		"movaps	0x10(%%eax),%%xmm3	\n\t"\
		"addpd	    (%%ebx),%%xmm0	\n\t"\
		"addpd	0x10(%%ebx),%%xmm1	\n\t"\
		"subpd	    (%%ebx),%%xmm2	\n\t"\
		"subpd	0x10(%%ebx),%%xmm3	\n\t"\
		"movaps	    (%%ecx),%%xmm4	\n\t"\
		"movaps	0x10(%%ecx),%%xmm5	\n\t"\
		"movaps	    (%%ecx),%%xmm6	\n\t"\
		"movaps	0x10(%%ecx),%%xmm7	\n\t"\
		"addpd	    (%%edx),%%xmm6	\n\t"\
		"addpd	0x10(%%edx),%%xmm7	\n\t"\
		"subpd	    (%%edx),%%xmm4	\n\t"\
		"subpd	0x10(%%edx),%%xmm5	\n\t"\
		"\n\t"/* Finish radix-4 butterfly and store results into temporary-array slots: */\
		"subpd	%%xmm6,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm1		\n\t"\
		"subpd	%%xmm5,%%xmm2		\n\t"\
		"movaps	%%xmm0,    (%%ebx)	\n\t"\
		"movaps	%%xmm1,0x10(%%ebx)	\n\t"\
		"subpd	%%xmm4,%%xmm3		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm0,%%xmm6		\n\t"\
		"addpd	%%xmm2,%%xmm5		\n\t"\
		"addpd	%%xmm1,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm4		\n\t"\
		"movaps	%%xmm6,    (%%eax)	\n\t"\
		"movaps	%%xmm7,0x10(%%eax)	\n\t"\
		"\n\t"\
		"movaps	(%%edi),%%xmm0	\n\t"/* isrt2 */\
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
		"movaps	%%xmm2,    (%%ecx)	\n\t"\
		"movaps	%%xmm5,    (%%edx)	\n\t"\
		"movaps	%%xmm6,0x10(%%ecx)	\n\t"\
		"movaps	%%xmm7,0x10(%%edx)	\n\t"\
		"\n\t"/* SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS(r04,r14,r24,r34,r0C,r1C,r2C,r3C): */\
		"subl	$0x80,%%eax	\n\t"/* r04 */\
		"subl	$0x200,%%ebx	\n\t"/* r0C */\
		"addl	$0x80,%%ecx	\n\t"/* r24 */\
		"subl	$0x100,%%edx	\n\t"/* r2C */\
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	    (%%ecx),%%xmm4	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"movaps	0x10(%%ecx),%%xmm5	\n\t"\
		"movaps	    (%%ebx),%%xmm2	\n\t"\
		"movaps	0x10(%%edx),%%xmm7	\n\t"\
		"movaps	0x10(%%ebx),%%xmm3	\n\t"\
		"movaps	    (%%edx),%%xmm6	\n\t"\
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
		"movaps	%%xmm0,    (%%ebx)	\n\t"\
		"movaps	%%xmm4,    (%%ecx)	\n\t"\
		"movaps	%%xmm1,0x10(%%ebx)	\n\t"\
		"movaps	%%xmm5,0x10(%%edx)	\n\t"\
		"movaps	%%xmm2,    (%%eax)	\n\t"\
		"movaps	%%xmm7,    (%%edx)	\n\t"\
		"movaps	%%xmm3,0x10(%%eax)	\n\t"\
		"movaps	%%xmm6,0x10(%%ecx)	\n\t"\
		"\n\t"\
		"addl	$0x100,%%eax	\n\t"/* r14 */\
		"addl	$0x100,%%ebx	\n\t"/* r1C */\
		"addl	$0x100,%%ecx	\n\t"/* r34 */\
		"addl	$0x100,%%edx	\n\t"/* r3C */\
		"\n\t"\
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	    (%%ecx),%%xmm4	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"movaps	0x10(%%ecx),%%xmm5	\n\t"\
		"movaps	    (%%ebx),%%xmm2	\n\t"\
		"movaps	0x10(%%edx),%%xmm7	\n\t"\
		"movaps	0x10(%%ebx),%%xmm3	\n\t"\
		"movaps	    (%%edx),%%xmm6	\n\t"\
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
		"movaps	%%xmm0,    (%%ebx)	\n\t"\
		"movaps	%%xmm4,    (%%ecx)	\n\t"\
		"movaps	%%xmm1,0x10(%%ebx)	\n\t"\
		"movaps	%%xmm5,0x10(%%edx)	\n\t"\
		"movaps	%%xmm2,    (%%eax)	\n\t"\
		"movaps	%%xmm7,    (%%edx)	\n\t"\
		"movaps	%%xmm3,0x10(%%eax)	\n\t"\
		"movaps	%%xmm6,0x10(%%ecx)	\n\t"\
		"\n\t"\
		"\n\t"/* SSE2_RADIX4_DIF_IN_PLACE(r02,r22,r12,r32): */\
		"\n\t"\
		"subl	$0x120,%%eax	\n\t"/* r02 */\
		"addl	$0x60,%%ebx	\n\t"/* r22 */\
		"subl	$0x220,%%ecx	\n\t"/* r12 */\
		"subl	$0xa0,%%edx	\n\t"/* r32 */\
		"\n\t"\
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"movaps	    (%%eax),%%xmm2	\n\t"\
		"movaps	0x10(%%eax),%%xmm3	\n\t"\
		"addpd	    (%%ebx),%%xmm0	\n\t"\
		"addpd	0x10(%%ebx),%%xmm1	\n\t"\
		"subpd	    (%%ebx),%%xmm2	\n\t"\
		"subpd	0x10(%%ebx),%%xmm3	\n\t"\
		"\n\t"\
		"movaps	    (%%ecx),%%xmm4	\n\t"\
		"movaps	0x10(%%ecx),%%xmm5	\n\t"\
		"movaps	    (%%ecx),%%xmm6	\n\t"\
		"movaps	0x10(%%ecx),%%xmm7	\n\t"\
		"addpd	    (%%edx),%%xmm4	\n\t"\
		"addpd	0x10(%%edx),%%xmm5	\n\t"\
		"subpd	    (%%edx),%%xmm6	\n\t"\
		"subpd	0x10(%%edx),%%xmm7	\n\t"\
		"\n\t"/* Finish radix-4 butterfly and store results into temporary-array slots: */\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"movaps	%%xmm0,    (%%ebx)	\n\t"\
		"movaps	%%xmm1,0x10(%%ebx)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"movaps	%%xmm4,    (%%eax)	\n\t"\
		"movaps	%%xmm5,0x10(%%eax)	\n\t"\
		"\n\t"\
		"subpd	%%xmm7,%%xmm2		\n\t"\
		"subpd	%%xmm6,%%xmm3		\n\t"\
		"movaps	%%xmm2,    (%%ecx)	\n\t"\
		"movaps	%%xmm3,0x10(%%edx)	\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	%%xmm7,    (%%edx)	\n\t"\
		"movaps	%%xmm6,0x10(%%ecx)	\n\t"\
		"\n\t"\
		"\n\t"/* SSE2_RADIX4_DIF_IN_PLACE_2NDOFTWO(r0A,r2A,r1A,r3A): */\
		"addl	$0x80,%%eax		\n\t"\
		"addl	$0x80,%%ebx		\n\t"\
		"addl	$0x80,%%ecx		\n\t"\
		"addl	$0x80,%%edx		\n\t"\
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"movaps	    (%%eax),%%xmm2	\n\t"\
		"movaps	0x10(%%eax),%%xmm3	\n\t"\
		"addpd	    (%%ebx),%%xmm0	\n\t"\
		"addpd	0x10(%%ebx),%%xmm1	\n\t"\
		"subpd	    (%%ebx),%%xmm2	\n\t"\
		"subpd	0x10(%%ebx),%%xmm3	\n\t"\
		"movaps	    (%%ecx),%%xmm4	\n\t"\
		"movaps	0x10(%%ecx),%%xmm5	\n\t"\
		"movaps	    (%%ecx),%%xmm6	\n\t"\
		"movaps	0x10(%%ecx),%%xmm7	\n\t"\
		"addpd	    (%%edx),%%xmm6	\n\t"\
		"addpd	0x10(%%edx),%%xmm7	\n\t"\
		"subpd	    (%%edx),%%xmm4	\n\t"\
		"subpd	0x10(%%edx),%%xmm5	\n\t"\
		"\n\t"/* Finish radix-4 butterfly and store results into temporary-array slots: */\
		"subpd	%%xmm6,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm1		\n\t"\
		"subpd	%%xmm5,%%xmm2		\n\t"\
		"movaps	%%xmm0,    (%%ebx)	\n\t"\
		"movaps	%%xmm1,0x10(%%ebx)	\n\t"\
		"subpd	%%xmm4,%%xmm3		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm0,%%xmm6		\n\t"\
		"addpd	%%xmm2,%%xmm5		\n\t"\
		"addpd	%%xmm1,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm4		\n\t"\
		"movaps	%%xmm6,    (%%eax)	\n\t"\
		"movaps	%%xmm7,0x10(%%eax)	\n\t"\
		"\n\t"\
		"movaps	(%%edi),%%xmm0	\n\t"/* isrt2 */\
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
		"movaps	%%xmm2,    (%%ecx)	\n\t"\
		"movaps	%%xmm5,    (%%edx)	\n\t"\
		"movaps	%%xmm6,0x10(%%ecx)	\n\t"\
		"movaps	%%xmm7,0x10(%%edx)	\n\t"\
		"\n\t"/* SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS(r02,r12,r22,r32,r0A,r1A,r2A,r3A): */\
		"subl	$0x80,%%eax	\n\t"/* r02 */\
		"subl	$0x200,%%ebx	\n\t"/* r0A */\
		"addl	$0x80,%%ecx	\n\t"/* r22 */\
		"subl	$0x100,%%edx	\n\t"/* r2A */\
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	    (%%ecx),%%xmm4	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"movaps	0x10(%%ecx),%%xmm5	\n\t"\
		"movaps	    (%%ebx),%%xmm2	\n\t"\
		"movaps	0x10(%%edx),%%xmm7	\n\t"\
		"movaps	0x10(%%ebx),%%xmm3	\n\t"\
		"movaps	    (%%edx),%%xmm6	\n\t"\
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
		"movaps	%%xmm0,    (%%ebx)	\n\t"\
		"movaps	%%xmm4,    (%%ecx)	\n\t"\
		"movaps	%%xmm1,0x10(%%ebx)	\n\t"\
		"movaps	%%xmm5,0x10(%%edx)	\n\t"\
		"movaps	%%xmm2,    (%%eax)	\n\t"\
		"movaps	%%xmm7,    (%%edx)	\n\t"\
		"movaps	%%xmm3,0x10(%%eax)	\n\t"\
		"movaps	%%xmm6,0x10(%%ecx)	\n\t"\
		"\n\t"\
		"addl	$0x100,%%eax	\n\t"/* r12 */\
		"addl	$0x100,%%ebx	\n\t"/* r1A */\
		"addl	$0x100,%%ecx	\n\t"/* r32 */\
		"addl	$0x100,%%edx	\n\t"/* r3A */\
		"\n\t"\
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	    (%%ecx),%%xmm4	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"movaps	0x10(%%ecx),%%xmm5	\n\t"\
		"movaps	    (%%ebx),%%xmm2	\n\t"\
		"movaps	0x10(%%edx),%%xmm7	\n\t"\
		"movaps	0x10(%%ebx),%%xmm3	\n\t"\
		"movaps	    (%%edx),%%xmm6	\n\t"\
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
		"movaps	%%xmm0,    (%%ebx)	\n\t"\
		"movaps	%%xmm4,    (%%ecx)	\n\t"\
		"movaps	%%xmm1,0x10(%%ebx)	\n\t"\
		"movaps	%%xmm5,0x10(%%edx)	\n\t"\
		"movaps	%%xmm2,    (%%eax)	\n\t"\
		"movaps	%%xmm7,    (%%edx)	\n\t"\
		"movaps	%%xmm3,0x10(%%eax)	\n\t"\
		"movaps	%%xmm6,0x10(%%ecx)	\n\t"\
		"\n\t"\
		"\n\t"/* SSE2_RADIX4_DIF_IN_PLACE(r06,r26,r16,r36): */\
		"\n\t"\
		"subl	$0xC0,%%eax	\n\t"/* r06 */\
		"addl	$0xC0,%%ebx	\n\t"/* r26 */\
		"subl	$0x1C0,%%ecx	\n\t"/* r16 */\
		"subl	$0x40,%%edx	\n\t"/* r36 */\
		"\n\t"\
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"movaps	    (%%eax),%%xmm2	\n\t"\
		"movaps	0x10(%%eax),%%xmm3	\n\t"\
		"addpd	    (%%ebx),%%xmm0	\n\t"\
		"addpd	0x10(%%ebx),%%xmm1	\n\t"\
		"subpd	    (%%ebx),%%xmm2	\n\t"\
		"subpd	0x10(%%ebx),%%xmm3	\n\t"\
		"\n\t"\
		"movaps	    (%%ecx),%%xmm4	\n\t"\
		"movaps	0x10(%%ecx),%%xmm5	\n\t"\
		"movaps	    (%%ecx),%%xmm6	\n\t"\
		"movaps	0x10(%%ecx),%%xmm7	\n\t"\
		"addpd	    (%%edx),%%xmm4	\n\t"\
		"addpd	0x10(%%edx),%%xmm5	\n\t"\
		"subpd	    (%%edx),%%xmm6	\n\t"\
		"subpd	0x10(%%edx),%%xmm7	\n\t"\
		"\n\t"/* Finish radix-4 butterfly and store results into temporary-array slots: */\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"movaps	%%xmm0,    (%%ebx)	\n\t"\
		"movaps	%%xmm1,0x10(%%ebx)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"movaps	%%xmm4,    (%%eax)	\n\t"\
		"movaps	%%xmm5,0x10(%%eax)	\n\t"\
		"\n\t"\
		"subpd	%%xmm7,%%xmm2		\n\t"\
		"subpd	%%xmm6,%%xmm3		\n\t"\
		"movaps	%%xmm2,    (%%ecx)	\n\t"\
		"movaps	%%xmm3,0x10(%%edx)	\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	%%xmm7,    (%%edx)	\n\t"\
		"movaps	%%xmm6,0x10(%%ecx)	\n\t"\
		"\n\t"\
		"\n\t"/* SSE2_RADIX4_DIF_IN_PLACE_2NDOFTWO(r0E,r2E,r1E,r3E): */\
		"addl	$0x80,%%eax		\n\t"\
		"addl	$0x80,%%ebx		\n\t"\
		"addl	$0x80,%%ecx		\n\t"\
		"addl	$0x80,%%edx		\n\t"\
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"movaps	    (%%eax),%%xmm2	\n\t"\
		"movaps	0x10(%%eax),%%xmm3	\n\t"\
		"addpd	    (%%ebx),%%xmm0	\n\t"\
		"addpd	0x10(%%ebx),%%xmm1	\n\t"\
		"subpd	    (%%ebx),%%xmm2	\n\t"\
		"subpd	0x10(%%ebx),%%xmm3	\n\t"\
		"movaps	    (%%ecx),%%xmm4	\n\t"\
		"movaps	0x10(%%ecx),%%xmm5	\n\t"\
		"movaps	    (%%ecx),%%xmm6	\n\t"\
		"movaps	0x10(%%ecx),%%xmm7	\n\t"\
		"addpd	    (%%edx),%%xmm6	\n\t"\
		"addpd	0x10(%%edx),%%xmm7	\n\t"\
		"subpd	    (%%edx),%%xmm4	\n\t"\
		"subpd	0x10(%%edx),%%xmm5	\n\t"\
		"\n\t"/* Finish radix-4 butterfly and store results into temporary-array slots: */\
		"subpd	%%xmm6,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm1		\n\t"\
		"subpd	%%xmm5,%%xmm2		\n\t"\
		"movaps	%%xmm0,    (%%ebx)	\n\t"\
		"movaps	%%xmm1,0x10(%%ebx)	\n\t"\
		"subpd	%%xmm4,%%xmm3		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm0,%%xmm6		\n\t"\
		"addpd	%%xmm2,%%xmm5		\n\t"\
		"addpd	%%xmm1,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm4		\n\t"\
		"movaps	%%xmm6,    (%%eax)	\n\t"\
		"movaps	%%xmm7,0x10(%%eax)	\n\t"\
		"\n\t"\
		"movaps	(%%edi),%%xmm0	\n\t"/* isrt2 */\
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
		"movaps	%%xmm2,    (%%ecx)	\n\t"\
		"movaps	%%xmm5,    (%%edx)	\n\t"\
		"movaps	%%xmm6,0x10(%%ecx)	\n\t"\
		"movaps	%%xmm7,0x10(%%edx)	\n\t"\
		"\n\t"/* SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS(r06,r16,r26,r36,r0E,r1E,r2E,r3E): */\
		"subl	$0x80,%%eax	\n\t"/* r02 */\
		"subl	$0x200,%%ebx	\n\t"/* r0A */\
		"addl	$0x80,%%ecx	\n\t"/* r22 */\
		"subl	$0x100,%%edx	\n\t"/* r2A */\
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	    (%%ecx),%%xmm4	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"movaps	0x10(%%ecx),%%xmm5	\n\t"\
		"movaps	    (%%ebx),%%xmm2	\n\t"\
		"movaps	0x10(%%edx),%%xmm7	\n\t"\
		"movaps	0x10(%%ebx),%%xmm3	\n\t"\
		"movaps	    (%%edx),%%xmm6	\n\t"\
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
		"movaps	%%xmm0,    (%%ebx)	\n\t"\
		"movaps	%%xmm4,    (%%ecx)	\n\t"\
		"movaps	%%xmm1,0x10(%%ebx)	\n\t"\
		"movaps	%%xmm5,0x10(%%edx)	\n\t"\
		"movaps	%%xmm2,    (%%eax)	\n\t"\
		"movaps	%%xmm7,    (%%edx)	\n\t"\
		"movaps	%%xmm3,0x10(%%eax)	\n\t"\
		"movaps	%%xmm6,0x10(%%ecx)	\n\t"\
		"\n\t"\
		"addl	$0x100,%%eax	\n\t"/* r12 */\
		"addl	$0x100,%%ebx	\n\t"/* r1A */\
		"addl	$0x100,%%ecx	\n\t"/* r32 */\
		"addl	$0x100,%%edx	\n\t"/* r3A */\
		"\n\t"\
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	    (%%ecx),%%xmm4	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"movaps	0x10(%%ecx),%%xmm5	\n\t"\
		"movaps	    (%%ebx),%%xmm2	\n\t"\
		"movaps	0x10(%%edx),%%xmm7	\n\t"\
		"movaps	0x10(%%ebx),%%xmm3	\n\t"\
		"movaps	    (%%edx),%%xmm6	\n\t"\
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
		"movaps	%%xmm0,    (%%ebx)	\n\t"\
		"movaps	%%xmm4,    (%%ecx)	\n\t"\
		"movaps	%%xmm1,0x10(%%ebx)	\n\t"\
		"movaps	%%xmm5,0x10(%%edx)	\n\t"\
		"movaps	%%xmm2,    (%%eax)	\n\t"\
		"movaps	%%xmm7,    (%%edx)	\n\t"\
		"movaps	%%xmm3,0x10(%%eax)	\n\t"\
		"movaps	%%xmm6,0x10(%%ecx)	\n\t"\
		"\n\t"\
	/**********************************************************************************/\
	/*...and now do eight radix-4 transforms, including the internal twiddle factors: */\
	/**********************************************************************************/\
		"\n\t"\
	/*...Block 1: t00,t10,t20,t30 in r00,04,02,06 - note swapped middle 2 indices! */\
		"movl	%[__r00],%%esi	\n\t"\
		"movl	%[__arr_offsets],%%edi	\n\t"\
		"movl	0x0(%%edi),%%eax	\n\t"/* p00 */\
		"movl	0x4(%%edi),%%ebx	\n\t"/* p01 */\
		"movl	0x8(%%edi),%%ecx	\n\t"/* p02 */\
		"movl	0xc(%%edi),%%edx	\n\t"/* p03 */\
		"movl	%[__add],%%edi		\n\t"\
		"addl	%%edi,%%eax	\n\t"\
		"addl	%%edi,%%ebx	\n\t"\
		"addl	%%edi,%%ecx	\n\t"\
		"addl	%%edi,%%edx	\n\t"\
		"movaps	    (%%esi),%%xmm0	\n\t"/* t00 */\
		"movaps	0x20(%%esi),%%xmm4	\n\t"/* t20 */\
		"movaps	0x10(%%esi),%%xmm1	\n\t"/* t01 */\
		"movaps	0x30(%%esi),%%xmm5	\n\t"/* t21 */\
		"movaps	0x40(%%esi),%%xmm2	\n\t"/* t10 */\
		"movaps	0x60(%%esi),%%xmm6	\n\t"/* t30 */\
		"movaps	0x50(%%esi),%%xmm3	\n\t"/* t11 */\
		"movaps	0x70(%%esi),%%xmm7	\n\t"/* t31 */\
		"subpd	0x40(%%esi),%%xmm0	\n\t"/* t10=t00-rt */\
		"subpd	0x60(%%esi),%%xmm4	\n\t"/* t30=t20-rt */\
		"subpd	0x50(%%esi),%%xmm1	\n\t"/* t11=t01-it */\
		"subpd	0x70(%%esi),%%xmm5	\n\t"/* t31=t21-it */\
		"addpd	    (%%esi),%%xmm2	\n\t"/* t00=t00+rt */\
		"addpd	0x20(%%esi),%%xmm6	\n\t"/* t20=t20+rt */\
		"addpd	0x10(%%esi),%%xmm3	\n\t"/* t01=t01+it */\
		"addpd	0x30(%%esi),%%xmm7	\n\t"/* t21=t21+it */\
		"subpd	%%xmm6,%%xmm2		\n\t"/* t00 <- t00-t20 */\
		"subpd	%%xmm5,%%xmm0		\n\t"/* t10 <- t10-t31 */\
		"subpd	%%xmm7,%%xmm3		\n\t"/* t01 <- t01-t21 */\
		"subpd	%%xmm4,%%xmm1		\n\t"/* t11 <- t11-t30 */\
		"addpd	%%xmm6,%%xmm6		\n\t"/*          2*t20 */\
		"addpd	%%xmm5,%%xmm5		\n\t"/*          2*t31 */\
		"addpd	%%xmm7,%%xmm7		\n\t"/*          2*t21 */\
		"addpd	%%xmm4,%%xmm4		\n\t"/*          2*t30 */\
		"movaps	%%xmm2,    (%%ebx)	\n\t"/* a(jt+p1 ) */\
		"movaps	%%xmm0,    (%%ecx)	\n\t"/* a(jt+p2 ) */\
		"movaps	%%xmm3,0x10(%%ebx)	\n\t"/* a(jp+p1 ) */\
		"movaps	%%xmm1,0x10(%%edx)	\n\t"/* a(jp+p3 ) */\
		"addpd	%%xmm2,%%xmm6		\n\t"/* t20 <- t00+t20 */\
		"addpd	%%xmm0,%%xmm5		\n\t"/* t31 <- t10+t31 */\
		"addpd	%%xmm3,%%xmm7		\n\t"/* t21 <- t01+t21 */\
		"addpd	%%xmm1,%%xmm4		\n\t"/* t30 <- t11+t30 */\
		"movaps	%%xmm6,    (%%eax)	\n\t"/* a(jt+p0 ) */\
		"movaps	%%xmm5,    (%%edx)	\n\t"/* a(jt+p3 ) */\
		"movaps	%%xmm7,0x10(%%eax)	\n\t"/* a(jp+p0 ) */\
		"movaps	%%xmm4,0x10(%%ecx)	\n\t"/* a(jp+p2 ) */\
		"\n\t"/*...Block 5: t08,t18,t28,t38	*/\
		"addl	$0x80,%%esi			\n\t"/* r08 */\
		"movl	%[__arr_offsets],%%edi	\n\t"\
		"movl	0x10(%%edi),%%eax	\n\t"/* p04 */\
		"movl	0x14(%%edi),%%ebx	\n\t"/* p05 */\
		"movl	0x18(%%edi),%%ecx	\n\t"/* p06 */\
		"movl	0x1c(%%edi),%%edx	\n\t"/* p07 */\
		"movl	%[__add],%%edi		\n\t"\
		"addl	%%edi,%%eax	\n\t"\
		"addl	%%edi,%%ebx	\n\t"\
		"addl	%%edi,%%ecx	\n\t"\
		"addl	%%edi,%%edx	\n\t"\
		"movl	%[__isrt2],%%edi	\n\t"\
		"\n\t"\
		"movaps	0x20(%%esi),%%xmm4	\n\t"/* t28 */\
		"movaps	0x30(%%esi),%%xmm5	\n\t"/* t29 */\
		"movaps	0x60(%%esi),%%xmm6	\n\t"/* t38 */\
		"movaps	0x70(%%esi),%%xmm7	\n\t"/* t39 */\
		"movaps	    (%%esi),%%xmm0	\n\t"/* t08 */\
		"subpd	0x30(%%esi),%%xmm4	\n\t"/* t28-t29 */\
		"movaps	0x10(%%esi),%%xmm1	\n\t"/* t09 */\
		"addpd	0x20(%%esi),%%xmm5	\n\t"/* t29+t28 */\
		"movaps	0x40(%%esi),%%xmm2	\n\t"/* t18 */\
		"mulpd	    (%%edi),%%xmm4		\n\t"/* t28 = (t28-t29)*ISRT2 */\
		"movaps	0x50(%%esi),%%xmm3	\n\t"/* t19 */\
		"mulpd	    (%%edi),%%xmm5		\n\t"/* t29 = (t29+t28)*ISRT2 */\
		"subpd	0x50(%%esi),%%xmm0	\n\t"/* t08=t08-t19*/\
		"addpd	0x70(%%esi),%%xmm6	\n\t"/* t38+t39 */\
		"subpd	0x40(%%esi),%%xmm1	\n\t"/* t19=t09-t18*/\
		"subpd	0x60(%%esi),%%xmm7	\n\t"/* t39-t38 */\
		"addpd	0x10(%%esi),%%xmm2	\n\t"/* t09=t18+t09*/\
		"mulpd	    (%%edi),%%xmm6		\n\t"/*  rt = (t38+t39)*ISRT2 */\
		"addpd	    (%%esi),%%xmm3	\n\t"/* t18=t19+t08*/\
		"mulpd	    (%%edi),%%xmm7		\n\t"/*  it = (t39-t38)*ISRT2 */\
		"\n\t"\
		"subpd	%%xmm6,%%xmm4		\n\t"/* t28=t28-rt */\
		"subpd	%%xmm7,%%xmm5		\n\t"/* t29=t29-it */\
		"addpd	%%xmm6,%%xmm6		\n\t"/*      2* rt */\
		"addpd	%%xmm7,%%xmm7		\n\t"/*      2* it */\
		"addpd	%%xmm4,%%xmm6		\n\t"/* t38=t28+rt */\
		"addpd	%%xmm5,%%xmm7		\n\t"/* t39=t29+it */\
		"subpd	%%xmm4,%%xmm0		\n\t"/* t08-t28 */\
		"subpd	%%xmm5,%%xmm2		\n\t"/* t09-t29 */\
		"addpd	%%xmm4,%%xmm4		\n\t"/*   2*t28 */\
		"addpd	%%xmm5,%%xmm5		\n\t"/*   2*t29 */\
		"\n\t"\
		"subpd	%%xmm7,%%xmm3		\n\t"/* t18-t39 */\
		"subpd	%%xmm6,%%xmm1		\n\t"/* t19-t38 */\
		"addpd	%%xmm7,%%xmm7		\n\t"/*   2*t39 */\
		"addpd	%%xmm6,%%xmm6		\n\t"/*   2*t38 */\
		"movaps	%%xmm0,    (%%ebx)	\n\t"/* a(jt+p1 ) */\
		"movaps	%%xmm2,0x10(%%ebx)	\n\t"/* a(jp+p1 ) */\
		"movaps	%%xmm3,    (%%ecx)	\n\t"/* a(jt+p2 ) */\
		"movaps	%%xmm1,0x10(%%edx)	\n\t"/* a(jp+p3 ) */\
		"addpd	%%xmm0,%%xmm4		\n\t"/* t08+t28 */\
		"addpd	%%xmm2,%%xmm5		\n\t"/* t09+t29 */\
		"addpd	%%xmm3,%%xmm7		\n\t"/* t18+t39 */\
		"addpd	%%xmm1,%%xmm6		\n\t"/* t19+t38 */\
		"movaps	%%xmm4,    (%%eax)	\n\t"/* a(jt+p0 ) */\
		"movaps	%%xmm5,0x10(%%eax)	\n\t"/* a(jp+p0 ) */\
		"movaps	%%xmm7,    (%%edx)	\n\t"/* a(jt+p3 ) */\
		"movaps	%%xmm6,0x10(%%ecx)	\n\t"/* a(jp+p2 ) */\
		"\n\t"/*...Block 3: t04,t14,t24,t34	*/\
		"addl	$0x180,%%esi		\n\t"/* r20 */\
		"movl	%[__arr_offsets],%%edi	\n\t"\
		"movl	0x20(%%edi),%%eax	\n\t"/* p08 */\
		"movl	0x24(%%edi),%%ebx	\n\t"/* p09 */\
		"movl	0x28(%%edi),%%ecx	\n\t"/* p0a */\
		"movl	0x2c(%%edi),%%edx	\n\t"/* p0b */\
		"movl	%[__add],%%edi		\n\t"\
		"addl	%%edi,%%eax	\n\t"\
		"addl	%%edi,%%ebx	\n\t"\
		"addl	%%edi,%%ecx	\n\t"\
		"addl	%%edi,%%edx	\n\t"\
		"movl	%[__isrt2],%%edi	\n\t"\
		"addl	$0x10,%%edi			\n\t"/* cc0 */\
		"\n\t"\
		"movaps	0x20(%%esi),%%xmm4	\n\t"/* t24 */\
		"movaps	0x60(%%esi),%%xmm6	\n\t"/* t34 */\
		"movaps	0x30(%%esi),%%xmm5	\n\t"/* t25 */\
		"movaps	0x70(%%esi),%%xmm7	\n\t"/* t35 */\
		"movaps	0x20(%%esi),%%xmm0	\n\t"/* copy t24 */\
		"movaps	0x60(%%esi),%%xmm2	\n\t"/* copy t34 */\
		"movaps	0x30(%%esi),%%xmm1	\n\t"/* copy t25 */\
		"movaps	0x70(%%esi),%%xmm3	\n\t"/* copy t35 */\
		"mulpd	    (%%edi),%%xmm4	\n\t"/* t24*c */\
		"mulpd	0x10(%%edi),%%xmm6	\n\t"/* t34*s */\
		"mulpd	0x10(%%edi),%%xmm1	\n\t"/* t25*s */\
		"mulpd	    (%%edi),%%xmm3	\n\t"/* t35*c */\
		"mulpd	    (%%edi),%%xmm5	\n\t"/* t25*c */\
		"mulpd	0x10(%%edi),%%xmm7	\n\t"/* t35*s */\
		"mulpd	0x10(%%edi),%%xmm0	\n\t"/* t24*s */\
		"mulpd	    (%%edi),%%xmm2	\n\t"/* t34*c */\
		"subpd	%%xmm1,%%xmm4		\n\t"/* ~t24 */\
		"subpd	%%xmm3,%%xmm6		\n\t"/* rt */\
		"addpd	%%xmm0,%%xmm5		\n\t"/* ~t25 */\
		"addpd	%%xmm2,%%xmm7		\n\t"/* it */\
		"subl	$0x10,%%edi			\n\t"/* isrt2 */\
		"\n\t"\
		"movaps	0x40(%%esi),%%xmm2	\n\t"/* t14 */\
		"subpd	%%xmm6,%%xmm4		\n\t"/*~t34=t24-rt */\
		"movaps	0x50(%%esi),%%xmm3	\n\t"/* t15 */\
		"subpd	%%xmm7,%%xmm5		\n\t"/*~t35=t25-it */\
		"subpd	0x50(%%esi),%%xmm2	\n\t"/* t14-t15 */\
		"addpd	%%xmm6,%%xmm6		\n\t"/*      2* rt */\
		"addpd	0x40(%%esi),%%xmm3	\n\t"/* t15+t14 */\
		"addpd	%%xmm7,%%xmm7		\n\t"/*      2* it */\
		"mulpd	    (%%edi),%%xmm2		\n\t"/* rt = (t14-t15)*ISRT2 */\
		"addpd	%%xmm4,%%xmm6		\n\t"/*~t24=t24+rt */\
		"mulpd	    (%%edi),%%xmm3		\n\t"/* it = (t15+t14)*ISRT2 */\
		"addpd	%%xmm5,%%xmm7		\n\t"/*~t25=t25+it */\
		"movaps	    (%%esi),%%xmm0	\n\t"/* t04 */\
		"movaps	0x10(%%esi),%%xmm1	\n\t"/* t05 */\
		"subpd	%%xmm2,%%xmm0		\n\t"/*~t14=t04-rt */\
		"subpd	%%xmm3,%%xmm1		\n\t"/*~t15=t05-it */\
		"addpd	    (%%esi),%%xmm2	\n\t"/*~t04=rt +t04*/\
		"addpd	0x10(%%esi),%%xmm3	\n\t"/*~t05=it +t05*/\
		"subpd	%%xmm6,%%xmm2		\n\t"/* t04-t24 */\
		"subpd	%%xmm5,%%xmm0		\n\t"/* t14-t35 */\
		"subpd	%%xmm7,%%xmm3		\n\t"/* t05-t25 */\
		"subpd	%%xmm4,%%xmm1		\n\t"/* t15-t34 */\
		"addpd	%%xmm6,%%xmm6		\n\t"/*   2*t24 */\
		"addpd	%%xmm5,%%xmm5		\n\t"/*          2*t35 */\
		"addpd	%%xmm7,%%xmm7		\n\t"/*   2*t25 */\
		"addpd	%%xmm4,%%xmm4		\n\t"/*          2*t34 */\
		"movaps	%%xmm2,    (%%ebx)	\n\t"/* a(jt+p1 ) */\
		"movaps	%%xmm0,    (%%ecx)	\n\t"/* a(jt+p2 ) */\
		"movaps	%%xmm3,0x10(%%ebx)	\n\t"/* a(jp+p1 ) */\
		"movaps	%%xmm1,0x10(%%edx)	\n\t"/* a(jp+p3 ) */\
		"addpd	%%xmm2,%%xmm6		\n\t"/* t04+t24 */\
		"addpd	%%xmm0,%%xmm5		\n\t"/* t14+t35 */\
		"addpd	%%xmm3,%%xmm7		\n\t"/* t05+t25 */\
		"addpd	%%xmm1,%%xmm4		\n\t"/* t15+t34 */\
		"movaps	%%xmm6,    (%%eax)	\n\t"/* a(jt+p0 ) */\
		"movaps	%%xmm5,    (%%edx)	\n\t"/* a(jt+p3 ) */\
		"movaps	%%xmm7,0x10(%%eax)	\n\t"/* a(jp+p0 ) */\
		"movaps	%%xmm4,0x10(%%ecx)	\n\t"/* a(jp+p2 ) */\
		"\n\t"/*...Block 7: t0C,t1C,t2C,t3C	*/\
		"addl	$0x80,%%esi			\n\t"/* r28 */\
		"movl	%[__arr_offsets],%%edi	\n\t"\
		"movl	0x30(%%edi),%%eax	\n\t"/* p0c */\
		"movl	0x34(%%edi),%%ebx	\n\t"/* p0d */\
		"movl	0x38(%%edi),%%ecx	\n\t"/* p0e */\
		"movl	0x3c(%%edi),%%edx	\n\t"/* p0f */\
		"movl	%[__add],%%edi		\n\t"\
		"addl	%%edi,%%eax	\n\t"\
		"addl	%%edi,%%ebx	\n\t"\
		"addl	%%edi,%%ecx	\n\t"\
		"addl	%%edi,%%edx	\n\t"\
		"movl	%[__isrt2],%%edi	\n\t"\
		"addl	$0x10,%%edi			\n\t"/* cc0 */\
		"\n\t"\
		"movaps	0x20(%%esi),%%xmm4	\n\t"/* t2C */\
		"movaps	0x60(%%esi),%%xmm6	\n\t"/* t3C */\
		"movaps	0x30(%%esi),%%xmm5	\n\t"/* t2D */\
		"movaps	0x70(%%esi),%%xmm7	\n\t"/* t3D */\
		"movaps	0x20(%%esi),%%xmm0	\n\t"/* copy t2C */\
		"movaps	0x60(%%esi),%%xmm2	\n\t"/* copy t3C */\
		"movaps	0x30(%%esi),%%xmm1	\n\t"/* copy t2D */\
		"movaps	0x70(%%esi),%%xmm3	\n\t"/* copy t3D */\
		"mulpd	0x10(%%edi),%%xmm4	\n\t"/* t2C*s */\
		"mulpd	    (%%edi),%%xmm6	\n\t"/* t3C*c */\
		"mulpd	    (%%edi),%%xmm1	\n\t"/* t2D*c */\
		"mulpd	0x10(%%edi),%%xmm3	\n\t"/* t3D*s */\
		"mulpd	0x10(%%edi),%%xmm5	\n\t"/* t2D*s */\
		"mulpd	    (%%edi),%%xmm7	\n\t"/* t3D*c */\
		"mulpd	    (%%edi),%%xmm0	\n\t"/* t2C*c */\
		"mulpd	0x10(%%edi),%%xmm2	\n\t"/* t3C*s */\
		"subpd	%%xmm1,%%xmm4		\n\t"/* ~t24 */\
		"subpd	%%xmm3,%%xmm6		\n\t"/* rt */\
		"addpd	%%xmm0,%%xmm5		\n\t"/* ~t25 */\
		"addpd	%%xmm2,%%xmm7		\n\t"/* it */\
		"subl	$0x10,%%edi			\n\t"/* isrt2 */\
		"\n\t"\
		"movaps	0x40(%%esi),%%xmm2	\n\t"/* t14 */\
		"subpd	%%xmm6,%%xmm4		\n\t"/*~t2C=t2C-rt */\
		"movaps	0x50(%%esi),%%xmm3	\n\t"/* t1D */\
		"subpd	%%xmm7,%%xmm5		\n\t"/*~t2D=t2D-it */\
		"addpd	0x50(%%esi),%%xmm2	\n\t"/* t1C+t1D */\
		"addpd	%%xmm6,%%xmm6		\n\t"/*      2* rt */\
		"subpd	0x40(%%esi),%%xmm3	\n\t"/* t1D-t1C */\
		"addpd	%%xmm7,%%xmm7		\n\t"/*      2* it */\
		"mulpd	    (%%edi),%%xmm2		\n\t"/* rt = (t1C+t1D)*ISRT2 */\
		"addpd	%%xmm4,%%xmm6		\n\t"/*~t3C=t2C+rt */\
		"mulpd	    (%%edi),%%xmm3		\n\t"/* it = (t1D-t1C)*ISRT2 */\
		"addpd	%%xmm5,%%xmm7		\n\t"/*~t3D=t2D+it */\
		"movaps	    (%%esi),%%xmm0	\n\t"/* t0C */\
		"movaps	0x10(%%esi),%%xmm1	\n\t"/* t0D */\
		"subpd	%%xmm2,%%xmm0		\n\t"/*~t0C=t0C-rt */\
		"subpd	%%xmm3,%%xmm1		\n\t"/*~t0D=t0D-it */\
		"addpd	    (%%esi),%%xmm2	\n\t"/*~t1C=rt +t0C*/\
		"addpd	0x10(%%esi),%%xmm3	\n\t"/*~t1D=it +t0D*/\
		"subpd	%%xmm4,%%xmm0		\n\t"/* t0C-t2C */\
		"subpd	%%xmm7,%%xmm2		\n\t"/* t1C-t3D */\
		"subpd	%%xmm5,%%xmm1		\n\t"/* t0D-t2D */\
		"subpd	%%xmm6,%%xmm3		\n\t"/* t1D-t3C */\
		"addpd	%%xmm4,%%xmm4		\n\t"/*   2*t2C */\
		"addpd	%%xmm7,%%xmm7		\n\t"/*   2*t3D */\
		"addpd	%%xmm5,%%xmm5		\n\t"/*   2*t2D */\
		"addpd	%%xmm6,%%xmm6		\n\t"/*   2*t3C */\
		"movaps	%%xmm0,    (%%ebx)	\n\t"/* a(jt+p1 ) */\
		"movaps	%%xmm2,    (%%ecx)	\n\t"/* a(jt+p2 ) */\
		"movaps	%%xmm1,0x10(%%ebx)	\n\t"/* a(jp+p1 ) */\
		"movaps	%%xmm3,0x10(%%edx)	\n\t"/* a(jp+p3 ) */\
		"addpd	%%xmm0,%%xmm4		\n\t"/* t0C+t2C */\
		"addpd	%%xmm2,%%xmm7		\n\t"/* t1C+t3D */\
		"addpd	%%xmm1,%%xmm5		\n\t"/* t0D+t2D */\
		"addpd	%%xmm3,%%xmm6		\n\t"/* t1D+t3C */\
		"movaps	%%xmm4,    (%%eax)	\n\t"/* a(jt+p0 ) */\
		"movaps	%%xmm7,    (%%edx)	\n\t"/* a(jt+p3 ) */\
		"movaps	%%xmm5,0x10(%%eax)	\n\t"/* a(jp+p0 ) */\
		"movaps	%%xmm6,0x10(%%ecx)	\n\t"/* a(jp+p2 ) */\
		"\n\t"/*...Block 2: t02,t12,t22,t32	*/\
		"subl	$0x180,%%esi		\n\t"/* r10 */\
		"movl	%[__arr_offsets],%%edi	\n\t"\
		"movl	0x40(%%edi),%%eax	\n\t"/* p10 */\
		"movl	0x44(%%edi),%%ebx	\n\t"/* p11 */\
		"movl	0x48(%%edi),%%ecx	\n\t"/* p12 */\
		"movl	0x4c(%%edi),%%edx	\n\t"/* p13 */\
		"movl	%[__add],%%edi		\n\t"\
		"addl	%%edi,%%eax	\n\t"\
		"addl	%%edi,%%ebx	\n\t"\
		"addl	%%edi,%%ecx	\n\t"\
		"addl	%%edi,%%edx	\n\t"\
		"movl	%[__isrt2],%%edi	\n\t"\
		"addl	$0x10,%%edi			\n\t"/* cc0 */\
		"\n\t"\
		"movaps	0x20(%%esi),%%xmm4	\n\t"/* t22 */\
		"movaps	0x60(%%esi),%%xmm6	\n\t"/* t32 */\
		"movaps	0x30(%%esi),%%xmm5	\n\t"/* t23 */\
		"movaps	0x70(%%esi),%%xmm7	\n\t"/* t33 */\
		"movaps	0x20(%%esi),%%xmm0	\n\t"/* copy t22 */\
		"movaps	0x60(%%esi),%%xmm2	\n\t"/* copy t32 */\
		"movaps	0x30(%%esi),%%xmm1	\n\t"/* copy t23 */\
		"movaps	0x70(%%esi),%%xmm3	\n\t"/* copy t33 */\
		"mulpd	0x20(%%edi),%%xmm4	\n\t"/* t22*c32_1 */\
		"mulpd	0x40(%%edi),%%xmm6	\n\t"/* t32*c32_3 */\
		"mulpd	0x30(%%edi),%%xmm1	\n\t"/* t23*s32_1 */\
		"mulpd	0x50(%%edi),%%xmm3	\n\t"/* t33*s32_3 */\
		"mulpd	0x20(%%edi),%%xmm5	\n\t"/* t23*c32_1 */\
		"mulpd	0x40(%%edi),%%xmm7	\n\t"/* t33*c32_3 */\
		"mulpd	0x30(%%edi),%%xmm0	\n\t"/* t22*s32_1 */\
		"mulpd	0x50(%%edi),%%xmm2	\n\t"/* t32*s32_3 */\
		"subpd	%%xmm1,%%xmm4		\n\t"/* ~t22 */\
		"subpd	%%xmm3,%%xmm6		\n\t"/* rt */\
		"addpd	%%xmm0,%%xmm5		\n\t"/* ~t23 */\
		"addpd	%%xmm2,%%xmm7		\n\t"/* it */\
		"\n\t"\
		"movaps	0x40(%%esi),%%xmm2	\n\t"/* t12 */\
		"movaps	0x50(%%esi),%%xmm0	\n\t"/* t13 */\
		"movaps	0x40(%%esi),%%xmm1	\n\t"/* copy t12 */\
		"movaps	0x50(%%esi),%%xmm3	\n\t"/* copy t13 */\
		"subpd	%%xmm6,%%xmm4		\n\t"/*~t32=t22-rt */\
		"mulpd	    (%%edi),%%xmm2	\n\t"/* t12*c */\
		"subpd	%%xmm7,%%xmm5		\n\t"/*~t33=t23-it */\
		"mulpd	0x10(%%edi),%%xmm0	\n\t"/* t13*s */\
		"addpd	%%xmm6,%%xmm6		\n\t"/*      2* rt */\
		"mulpd	    (%%edi),%%xmm3	\n\t"/* t13*c */\
		"addpd	%%xmm7,%%xmm7		\n\t"/*      2* it */\
		"mulpd	0x10(%%edi),%%xmm1	\n\t"/* t12*s */\
		"addpd	%%xmm4,%%xmm6		\n\t"/*~t22=t22+rt */\
		"subpd	%%xmm0,%%xmm2		\n\t"/* rt */\
		"addpd	%%xmm5,%%xmm7		\n\t"/*~t23=t23+it */\
		"addpd	%%xmm1,%%xmm3		\n\t"/* it */\
		"\n\t"\
		"movaps	    (%%esi),%%xmm0	\n\t"/* t02 */\
		"movaps	0x10(%%esi),%%xmm1	\n\t"/* t03 */\
		"subpd	%%xmm2,%%xmm0		\n\t"/*~t12=t02-rt */\
		"subpd	%%xmm3,%%xmm1		\n\t"/*~t13=t03-it */\
		"addpd	    (%%esi),%%xmm2	\n\t"/*~t02=rt+t02 */\
		"addpd	0x10(%%esi),%%xmm3	\n\t"/*~t03=it+t03 */\
		"subpd	%%xmm6,%%xmm2		\n\t"/* t02-t22 */\
		"subpd	%%xmm5,%%xmm0		\n\t"/* t12-t33 */\
		"subpd	%%xmm7,%%xmm3		\n\t"/* t03-t23 */\
		"subpd	%%xmm4,%%xmm1		\n\t"/* t13-t32 */\
		"addpd	%%xmm6,%%xmm6		\n\t"/*   2*t22 */\
		"addpd	%%xmm5,%%xmm5		\n\t"/*   2*t33 */\
		"addpd	%%xmm7,%%xmm7		\n\t"/*   2*t23 */\
		"addpd	%%xmm4,%%xmm4		\n\t"/*   2*t32 */\
		"movaps	%%xmm2,    (%%ebx)	\n\t"/* a(jt+p1 ) */\
		"movaps	%%xmm0,    (%%ecx)	\n\t"/* a(jt+p2 ) */\
		"movaps	%%xmm3,0x10(%%ebx)	\n\t"/* a(jp+p1 ) */\
		"movaps	%%xmm1,0x10(%%edx)	\n\t"/* a(jp+p3 ) */\
		"addpd	%%xmm2,%%xmm6		\n\t"/* t02+t22 */\
		"addpd	%%xmm0,%%xmm5		\n\t"/* t12+t33 */\
		"addpd	%%xmm3,%%xmm7		\n\t"/* t03+t23 */\
		"addpd	%%xmm1,%%xmm4		\n\t"/* t13+t32 */\
		"movaps	%%xmm6,    (%%eax)	\n\t"/* a(jt+p0 ) */\
		"movaps	%%xmm5,    (%%edx)	\n\t"/* a(jt+p3 ) */\
		"movaps	%%xmm7,0x10(%%eax)	\n\t"/* a(jp+p0 ) */\
		"movaps	%%xmm4,0x10(%%ecx)	\n\t"/* a(jp+p2 ) */\
		"\n\t"/*...Block 6: t0A,t1A,t2A,t3A	*/\
		"addl	$0x80,%%esi			\n\t"/* r18 */\
		"movl	%[__arr_offsets],%%edi	\n\t"\
		"movl	0x50(%%edi),%%eax	\n\t"/* p14 */\
		"movl	0x54(%%edi),%%ebx	\n\t"/* p15 */\
		"movl	0x58(%%edi),%%ecx	\n\t"/* p16 */\
		"movl	0x5c(%%edi),%%edx	\n\t"/* p17 */\
		"movl	%[__add],%%edi		\n\t"\
		"addl	%%edi,%%eax	\n\t"\
		"addl	%%edi,%%ebx	\n\t"\
		"addl	%%edi,%%ecx	\n\t"\
		"addl	%%edi,%%edx	\n\t"\
		"movl	%[__isrt2],%%edi	\n\t"\
		"addl	$0x10,%%edi			\n\t"/* cc0 */\
		"\n\t"\
		"movaps	0x20(%%esi),%%xmm4	\n\t"/* t2A */\
		"movaps	0x60(%%esi),%%xmm6	\n\t"/* t3A */\
		"movaps	0x30(%%esi),%%xmm5	\n\t"/* t2B */\
		"movaps	0x70(%%esi),%%xmm7	\n\t"/* t3B */\
		"movaps	0x20(%%esi),%%xmm0	\n\t"/* copy t2A */\
		"movaps	0x60(%%esi),%%xmm2	\n\t"/* copy t3A */\
		"movaps	0x30(%%esi),%%xmm1	\n\t"/* copy t2B */\
		"movaps	0x70(%%esi),%%xmm3	\n\t"/* copy t3B */\
		"mulpd	0x50(%%edi),%%xmm4	\n\t"/* t2A*s32_3 */\
		"mulpd	0x20(%%edi),%%xmm6	\n\t"/* t3A*c32_1 */\
		"mulpd	0x40(%%edi),%%xmm1	\n\t"/* t2B*c32_3 */\
		"mulpd	0x30(%%edi),%%xmm3	\n\t"/* t3B*s32_1 */\
		"mulpd	0x50(%%edi),%%xmm5	\n\t"/* t2B*s32_3 */\
		"mulpd	0x20(%%edi),%%xmm7	\n\t"/* t3B*c32_1 */\
		"mulpd	0x40(%%edi),%%xmm0	\n\t"/* t2A*c32_3 */\
		"mulpd	0x30(%%edi),%%xmm2	\n\t"/* t3A*s32_1 */\
		"subpd	%%xmm1,%%xmm4		\n\t"/* ~t2A */\
		"addpd	%%xmm3,%%xmm6		\n\t"/* rt */\
		"addpd	%%xmm0,%%xmm5		\n\t"/* ~t2B */\
		"subpd	%%xmm2,%%xmm7		\n\t"/* it */\
		"\n\t"\
		"movaps	0x40(%%esi),%%xmm2	\n\t"/* t1A */\
		"movaps	0x50(%%esi),%%xmm0	\n\t"/* t1B */\
		"movaps	0x40(%%esi),%%xmm1	\n\t"/* copy t1A */\
		"movaps	0x50(%%esi),%%xmm3	\n\t"/* copy t1B */\
		"subpd	%%xmm6,%%xmm4		\n\t"/*~t2A=t2A-rt */\
		"mulpd	0x10(%%edi),%%xmm2	\n\t"/* t1A*s */\
		"subpd	%%xmm7,%%xmm5		\n\t"/*~t2B=t2B-it */\
		"mulpd	    (%%edi),%%xmm0	\n\t"/* t1B*c */\
		"addpd	%%xmm6,%%xmm6		\n\t"/*      2* rt */\
		"mulpd	0x10(%%edi),%%xmm3	\n\t"/* t1B*s */\
		"addpd	%%xmm7,%%xmm7		\n\t"/*      2* it */\
		"mulpd	    (%%edi),%%xmm1	\n\t"/* t1A*c */\
		"addpd	%%xmm4,%%xmm6		\n\t"/*~t3A=t2A+rt */\
		"addpd	%%xmm0,%%xmm2		\n\t"/* rt */\
		"addpd	%%xmm5,%%xmm7		\n\t"/*~t3B=t2B+it */\
		"subpd	%%xmm1,%%xmm3		\n\t"/* it */\
		"\n\t"\
		"movaps	    (%%esi),%%xmm0	\n\t"/* t0A */\
		"movaps	0x10(%%esi),%%xmm1	\n\t"/* t0B */\
		"subpd	%%xmm2,%%xmm0		\n\t"/*~t0A=t0A-rt */\
		"subpd	%%xmm3,%%xmm1		\n\t"/*~t0B=t0B-it */\
		"addpd	    (%%esi),%%xmm2	\n\t"/*~t1A=rt+t0A */\
		"addpd	0x10(%%esi),%%xmm3	\n\t"/*~t1B=it+t0B */\
		"subpd	%%xmm4,%%xmm0		\n\t"/* t0A-t2A */\
		"subpd	%%xmm7,%%xmm2		\n\t"/* t1A-t3B */\
		"subpd	%%xmm5,%%xmm1		\n\t"/* t0B-t2B */\
		"subpd	%%xmm6,%%xmm3		\n\t"/* t1B-t3A */\
		"addpd	%%xmm4,%%xmm4		\n\t"/*   2*t2A */\
		"addpd	%%xmm7,%%xmm7		\n\t"/*   2*t3B */\
		"addpd	%%xmm5,%%xmm5		\n\t"/*   2*t2B */\
		"addpd	%%xmm6,%%xmm6		\n\t"/*   2*t3A */\
		"movaps	%%xmm0,    (%%ebx)	\n\t"/* a(jt+p1 ) */\
		"movaps	%%xmm2,    (%%ecx)	\n\t"/* a(jt+p2 ) */\
		"movaps	%%xmm1,0x10(%%ebx)	\n\t"/* a(jp+p1 ) */\
		"movaps	%%xmm3,0x10(%%edx)	\n\t"/* a(jp+p3 ) */\
		"addpd	%%xmm0,%%xmm4		\n\t"/* t0A+t2A */\
		"addpd	%%xmm2,%%xmm7		\n\t"/* t1A+t3B */\
		"addpd	%%xmm1,%%xmm5		\n\t"/* t0B+t2B */\
		"addpd	%%xmm3,%%xmm6		\n\t"/* t1B+t3A */\
		"movaps	%%xmm4,    (%%eax)	\n\t"/* a(jt+p0 ) */\
		"movaps	%%xmm7,    (%%edx)	\n\t"/* a(jt+p3 ) */\
		"movaps	%%xmm5,0x10(%%eax)	\n\t"/* a(jp+p0 ) */\
		"movaps	%%xmm6,0x10(%%ecx)	\n\t"/* a(jp+p2 ) */\
		"\n\t"/*...Block 4: t06,t16,t26,t36	*/\
		"addl	$0x180,%%esi		\n\t"/* r30 */\
		"movl	%[__arr_offsets],%%edi	\n\t"\
		"movl	0x60(%%edi),%%eax	\n\t"/* p18 */\
		"movl	0x64(%%edi),%%ebx	\n\t"/* p19 */\
		"movl	0x68(%%edi),%%ecx	\n\t"/* p1a */\
		"movl	0x6c(%%edi),%%edx	\n\t"/* p1b */\
		"movl	%[__add],%%edi		\n\t"\
		"addl	%%edi,%%eax	\n\t"\
		"addl	%%edi,%%ebx	\n\t"\
		"addl	%%edi,%%ecx	\n\t"\
		"addl	%%edi,%%edx	\n\t"\
		"movl	%[__isrt2],%%edi	\n\t"\
		"addl	$0x10,%%edi			\n\t"/* cc0 */\
		"\n\t"\
		"movaps	0x20(%%esi),%%xmm4	\n\t"/* t26 */\
		"movaps	0x60(%%esi),%%xmm6	\n\t"/* t36 */\
		"movaps	0x30(%%esi),%%xmm5	\n\t"/* t27 */\
		"movaps	0x70(%%esi),%%xmm7	\n\t"/* t37 */\
		"movaps	0x20(%%esi),%%xmm0	\n\t"/* copy t26 */\
		"movaps	0x60(%%esi),%%xmm2	\n\t"/* copy t36 */\
		"movaps	0x30(%%esi),%%xmm1	\n\t"/* copy t27 */\
		"movaps	0x70(%%esi),%%xmm3	\n\t"/* copy t37 */\
		"mulpd	0x40(%%edi),%%xmm4	\n\t"/* t26*s32_3 */\
		"mulpd	0x30(%%edi),%%xmm6	\n\t"/* t36*s32_1 */\
		"mulpd	0x50(%%edi),%%xmm1	\n\t"/* t27*s32_3 */\
		"mulpd	0x20(%%edi),%%xmm3	\n\t"/* t37*c32_1 */\
		"mulpd	0x40(%%edi),%%xmm5	\n\t"/* t27*c32_3 */\
		"mulpd	0x30(%%edi),%%xmm7	\n\t"/* t37*s32_1 */\
		"mulpd	0x50(%%edi),%%xmm0	\n\t"/* t26*s32_3 */\
		"mulpd	0x20(%%edi),%%xmm2	\n\t"/* t36*c32_1 */\
		"subpd	%%xmm1,%%xmm4		\n\t"/* ~t26 */\
		"addpd	%%xmm3,%%xmm6		\n\t"/* rt */\
		"addpd	%%xmm0,%%xmm5		\n\t"/* ~t27 */\
		"subpd	%%xmm2,%%xmm7		\n\t"/* it */\
		"\n\t"\
		"movaps	0x40(%%esi),%%xmm2	\n\t"/* t16 */\
		"movaps	0x50(%%esi),%%xmm0	\n\t"/* t17 */\
		"movaps	0x40(%%esi),%%xmm1	\n\t"/* copy t16 */\
		"movaps	0x50(%%esi),%%xmm3	\n\t"/* copy t17 */\
		"subpd	%%xmm6,%%xmm4		\n\t"/*~t26=t26-rt */\
		"mulpd	0x10(%%edi),%%xmm2	\n\t"/* t16*s */\
		"subpd	%%xmm7,%%xmm5		\n\t"/*~t27=t27-it */\
		"mulpd	    (%%edi),%%xmm0	\n\t"/* t17*c */\
		"addpd	%%xmm6,%%xmm6		\n\t"/*      2* rt */\
		"mulpd	0x10(%%edi),%%xmm3	\n\t"/* t17*s */\
		"addpd	%%xmm7,%%xmm7		\n\t"/*      2* it */\
		"mulpd	    (%%edi),%%xmm1	\n\t"/* t16*c */\
		"addpd	%%xmm4,%%xmm6		\n\t"/*~t36=t26+rt */\
		"subpd	%%xmm0,%%xmm2		\n\t"/* rt */\
		"addpd	%%xmm5,%%xmm7		\n\t"/*~t37=t27+it */\
		"addpd	%%xmm1,%%xmm3		\n\t"/* it */\
		"\n\t"\
		"movaps	    (%%esi),%%xmm0	\n\t"/* t06 */\
		"movaps	0x10(%%esi),%%xmm1	\n\t"/* t07 */\
		"subpd	%%xmm2,%%xmm0		\n\t"/*~t16=t06-rt */\
		"subpd	%%xmm3,%%xmm1		\n\t"/*~t17=t07-it */\
		"addpd	    (%%esi),%%xmm2	\n\t"/*~t06=rt+t06 */\
		"addpd	0x10(%%esi),%%xmm3	\n\t"/*~t07=it+t07 */\
		"subpd	%%xmm4,%%xmm2		\n\t"/* t06-t26 */\
		"subpd	%%xmm7,%%xmm0		\n\t"/* t16-t37 */\
		"subpd	%%xmm5,%%xmm3		\n\t"/* t07-t27 */\
		"subpd	%%xmm6,%%xmm1		\n\t"/* t17-t36 */\
		"addpd	%%xmm4,%%xmm4		\n\t"/*   2*t26 */\
		"addpd	%%xmm7,%%xmm7		\n\t"/*   2*t37 */\
		"addpd	%%xmm5,%%xmm5		\n\t"/*   2*t27 */\
		"addpd	%%xmm6,%%xmm6		\n\t"/*   2*t36 */\
		"movaps	%%xmm2,    (%%ebx)	\n\t"/* a(jt+p1 ) */\
		"movaps	%%xmm0,    (%%ecx)	\n\t"/* a(jt+p2 ) */\
		"movaps	%%xmm3,0x10(%%ebx)	\n\t"/* a(jp+p1 ) */\
		"movaps	%%xmm1,0x10(%%edx)	\n\t"/* a(jp+p3 ) */\
		"addpd	%%xmm2,%%xmm4		\n\t"/* t06+t26 */\
		"addpd	%%xmm0,%%xmm7		\n\t"/* t16+t37 */\
		"addpd	%%xmm3,%%xmm5		\n\t"/* t07+t27 */\
		"addpd	%%xmm1,%%xmm6		\n\t"/* t17+t36 */\
		"movaps	%%xmm4,    (%%eax)	\n\t"/* a(jt+p0 ) */\
		"movaps	%%xmm7,    (%%edx)	\n\t"/* a(jt+p3 ) */\
		"movaps	%%xmm5,0x10(%%eax)	\n\t"/* a(jp+p0 ) */\
		"movaps	%%xmm6,0x10(%%ecx)	\n\t"/* a(jp+p2 ) */\
		"\n\t"/*...Block 8: t0E,t1E,t2E,t3E	*/\
		"addl	$0x80,%%esi			\n\t"/* r38 */\
		"movl	%[__arr_offsets],%%edi	\n\t"\
		"movl	0x70(%%edi),%%eax	\n\t"/* p1c */\
		"movl	0x74(%%edi),%%ebx	\n\t"/* p1d */\
		"movl	0x78(%%edi),%%ecx	\n\t"/* p1e */\
		"movl	0x7c(%%edi),%%edx	\n\t"/* p1f */\
		"movl	%[__add],%%edi		\n\t"\
		"addl	%%edi,%%eax	\n\t"\
		"addl	%%edi,%%ebx	\n\t"\
		"addl	%%edi,%%ecx	\n\t"\
		"addl	%%edi,%%edx	\n\t"\
		"movl	%[__isrt2],%%edi	\n\t"\
		"addl	$0x10,%%edi			\n\t"/* cc0 */\
		"\n\t"\
		"movaps	0x20(%%esi),%%xmm4	\n\t"/* t2E */\
		"movaps	0x60(%%esi),%%xmm6	\n\t"/* t3E */\
		"movaps	0x30(%%esi),%%xmm5	\n\t"/* t2F */\
		"movaps	0x70(%%esi),%%xmm7	\n\t"/* t3F */\
		"movaps	0x20(%%esi),%%xmm0	\n\t"/* copy t2E */\
		"movaps	0x60(%%esi),%%xmm2	\n\t"/* copy t3E */\
		"movaps	0x30(%%esi),%%xmm1	\n\t"/* copy t2F */\
		"movaps	0x70(%%esi),%%xmm3	\n\t"/* copy t3F */\
		"mulpd	0x30(%%edi),%%xmm4	\n\t"/* t2E*s32_1 */\
		"mulpd	0x50(%%edi),%%xmm6	\n\t"/* t3E*c32_3 */\
		"mulpd	0x20(%%edi),%%xmm1	\n\t"/* t2F*c32_1 */\
		"mulpd	0x40(%%edi),%%xmm3	\n\t"/* t3F*s32_3 */\
		"mulpd	0x30(%%edi),%%xmm5	\n\t"/* t2F*s32_1 */\
		"mulpd	0x50(%%edi),%%xmm7	\n\t"/* t3F*c32_3 */\
		"mulpd	0x20(%%edi),%%xmm0	\n\t"/* t2E*c32_1 */\
		"mulpd	0x40(%%edi),%%xmm2	\n\t"/* t3E*s32_3 */\
		"subpd	%%xmm1,%%xmm4		\n\t"/* ~t2E */\
		"subpd	%%xmm3,%%xmm6		\n\t"/* rt */\
		"addpd	%%xmm0,%%xmm5		\n\t"/* ~t2F */\
		"addpd	%%xmm2,%%xmm7		\n\t"/* it */\
		"\n\t"\
		"movaps	0x40(%%esi),%%xmm2	\n\t"/* t1E */\
		"movaps	0x50(%%esi),%%xmm0	\n\t"/* t1F */\
		"movaps	0x40(%%esi),%%xmm1	\n\t"/* copy t1E */\
		"movaps	0x50(%%esi),%%xmm3	\n\t"/* copy t1F */\
		"subpd	%%xmm6,%%xmm4		\n\t"/*~t2E=t2E-rt */\
		"mulpd	    (%%edi),%%xmm2	\n\t"/* t1E*c */\
		"subpd	%%xmm7,%%xmm5		\n\t"/*~t2F=t2F-it */\
		"mulpd	0x10(%%edi),%%xmm0	\n\t"/* t1F*s */\
		"addpd	%%xmm6,%%xmm6		\n\t"/*      2* rt */\
		"mulpd	    (%%edi),%%xmm3	\n\t"/* t1F*c */\
		"addpd	%%xmm7,%%xmm7		\n\t"/*      2* it */\
		"mulpd	0x10(%%edi),%%xmm1	\n\t"/* t1E*s */\
		"addpd	%%xmm4,%%xmm6		\n\t"/*~t3E=t2E+rt */\
		"addpd	%%xmm0,%%xmm2		\n\t"/* rt */\
		"addpd	%%xmm5,%%xmm7		\n\t"/*~t3F=t2F+it */\
		"subpd	%%xmm1,%%xmm3		\n\t"/* it */\
		"\n\t"\
		"movaps	    (%%esi),%%xmm0	\n\t"/* t0E */\
		"movaps	0x10(%%esi),%%xmm1	\n\t"/* t0F */\
		"subpd	%%xmm2,%%xmm0		\n\t"/*~t0E=t0E-rt */\
		"subpd	%%xmm3,%%xmm1		\n\t"/*~t0F=t0F-it */\
		"addpd	    (%%esi),%%xmm2	\n\t"/*~t1E=rt+t0E */\
		"addpd	0x10(%%esi),%%xmm3	\n\t"/*~t1F=it+t0F */\
		"subpd	%%xmm4,%%xmm0		\n\t"/* t0E-t2E */\
		"subpd	%%xmm7,%%xmm2		\n\t"/* t1E-t3F */\
		"subpd	%%xmm5,%%xmm1		\n\t"/* t0F-t2F */\
		"subpd	%%xmm6,%%xmm3		\n\t"/* t1F-t3E */\
		"addpd	%%xmm4,%%xmm4		\n\t"/*   2*t2E */\
		"addpd	%%xmm7,%%xmm7		\n\t"/*   2*t3F */\
		"addpd	%%xmm5,%%xmm5		\n\t"/*   2*t2F */\
		"addpd	%%xmm6,%%xmm6		\n\t"/*   2*t3E */\
		"movaps	%%xmm0,    (%%ebx)	\n\t"/* a(jt+p1 ) */\
		"movaps	%%xmm2,    (%%ecx)	\n\t"/* a(jt+p2 ) */\
		"movaps	%%xmm1,0x10(%%ebx)	\n\t"/* a(jp+p1 ) */\
		"movaps	%%xmm3,0x10(%%edx)	\n\t"/* a(jp+p3 ) */\
		"addpd	%%xmm0,%%xmm4		\n\t"/* t0E+t2E */\
		"addpd	%%xmm2,%%xmm7		\n\t"/* t1E+t3F */\
		"addpd	%%xmm1,%%xmm5		\n\t"/* t0F+t2F */\
		"addpd	%%xmm3,%%xmm6		\n\t"/* t1F+t3E */\
		"movaps	%%xmm4,    (%%eax)	\n\t"/* a(jt+p0 ) */\
		"movaps	%%xmm7,    (%%edx)	\n\t"/* a(jt+p3 ) */\
		"movaps	%%xmm5,0x10(%%eax)	\n\t"/* a(jp+p0 ) */\
		"movaps	%%xmm6,0x10(%%ecx)	\n\t"/* a(jp+p2 ) */\
	"popl %%ebx	\n\t"\
		:					/* outputs: none */\
		: [__add] "m" (Xadd)	/* All inputs from memory addresses here */\
		 ,[__arr_offsets] "m" (Xarr_offsets)\
		 ,[__r00] "m" (Xr00)\
		 ,[__isrt2] "m" (Xisrt2)\
		: "cc","memory","eax",/*"ebx",*/"ecx","edx","edi","esi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"	/* Clobbered registers */\
	);\
	}

#endif	/* radix32_ditN_cy_dif1_gcc_h_included */

