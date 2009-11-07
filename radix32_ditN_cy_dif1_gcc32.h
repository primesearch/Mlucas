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
		"movl	%[__r00],%%esi	\n\t"\
		"movl	%[__add],%%eax	\n\t"\
		"movl	%[__p01],%%ebx	\n\t"\
		"movl	%[__p02],%%ecx	\n\t"\
		"movl	%[__p03],%%edx	\n\t"\
		"shll	$3,%%ebx		\n\t"\
		"shll	$3,%%ecx		\n\t"\
		"shll	$3,%%edx		\n\t"\
		"addl	%%eax,%%ebx		\n\t"\
		"addl	%%eax,%%ecx		\n\t"\
		"addl	%%eax,%%edx		\n\t"\
		"/* SSE2_RADIX4_DIT_0TWIDDLE(add0, add1, add2, add3, r00): */\n\t"\
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
		"/* SSE2_RADIX4_DIT_0TWIDDLE_2NDOFTWO(add4, add5, add6, add7, r08): */\n\t"\
		"addl	$0x80,%%esi			\n\t"\
		"movl	%[__p04],%%edi		\n\t"\
		"shll	$3,%%edi			\n\t"\
		"addl	%%edi,%%eax			\n\t"\
		"addl	%%edi,%%ebx			\n\t"\
		"addl	%%edi,%%ecx			\n\t"\
		"addl	%%edi,%%edx			\n\t"\
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
		"movaps	%%xmm0,0x040(%%esi)	\n\t"\
		"movaps	%%xmm1,0x050(%%esi)	\n\t"\
		"subpd	%%xmm7,%%xmm2		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"subpd	%%xmm6,%%xmm3		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"movaps	%%xmm4,     (%%esi)	\n\t"\
		"addpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	%%xmm5,0x010(%%esi)	\n\t"\
		"\n\t"\
		"movaps	(%%edi),%%xmm5	/* ISRT2 */\n\t"\
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
		"movaps	%%xmm3,0x030(%%esi)	\n\t"\
		"movaps	%%xmm6,0x070(%%esi)	\n\t"\
		"movaps	%%xmm0,0x020(%%esi)	\n\t"\
		"movaps	%%xmm1,0x060(%%esi)	\n\t"\
		"/* SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS(r00,r02,r04,r06,r08,r0A,r0C,r0E): */\n\t"\
		"movaps		-0x080(%%esi),%%xmm0	/* r00 */\n\t"\
		"movaps		-0x040(%%esi),%%xmm4	/* r04 */\n\t"\
		"movaps		-0x070(%%esi),%%xmm1	/* r01 */\n\t"\
		"movaps		-0x030(%%esi),%%xmm5	/* r05 */\n\t"\
		"movaps		      (%%esi),%%xmm2	/* r08 */\n\t"\
		"movaps		 0x050(%%esi),%%xmm7	/* r0D */\n\t"\
		"movaps		 0x010(%%esi),%%xmm3	/* r09 */\n\t"\
		"movaps		 0x040(%%esi),%%xmm6	/* r0C */\n\t"\
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
		"movaps		%%xmm0,      (%%esi)\n\t"\
		"movaps		%%xmm4, 0x040(%%esi)\n\t"\
		"movaps		%%xmm1, 0x010(%%esi)\n\t"\
		"movaps		%%xmm5,-0x030(%%esi)\n\t"\
		"movaps		%%xmm2,-0x080(%%esi)\n\t"\
		"movaps		%%xmm7,-0x040(%%esi)\n\t"\
		"movaps		%%xmm3,-0x070(%%esi)\n\t"\
		"movaps		%%xmm6, 0x050(%%esi)\n\t"\
		"\n\t"\
		"movaps		-0x060(%%esi),%%xmm0	/* r02 */\n\t"\
		"movaps		-0x020(%%esi),%%xmm4	/* r06 */\n\t"\
		"movaps		-0x050(%%esi),%%xmm1	/* r03 */\n\t"\
		"movaps		-0x010(%%esi),%%xmm5	/* r07 */\n\t"\
		"movaps		 0x020(%%esi),%%xmm2	/* r0A */\n\t"\
		"movaps		 0x070(%%esi),%%xmm7	/* r0F */\n\t"\
		"movaps		 0x030(%%esi),%%xmm3	/* r0B */\n\t"\
		"movaps		 0x060(%%esi),%%xmm6	/* r0E */\n\t"\
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
		"movaps		%%xmm0, 0x020(%%esi)\n\t"\
		"movaps		%%xmm4, 0x060(%%esi)\n\t"\
		"movaps		%%xmm1, 0x030(%%esi)\n\t"\
		"movaps		%%xmm5,-0x010(%%esi)\n\t"\
		"movaps		%%xmm2,-0x060(%%esi)\n\t"\
		"movaps		%%xmm7,-0x020(%%esi)\n\t"\
		"movaps		%%xmm3,-0x050(%%esi)\n\t"\
		"movaps		%%xmm6, 0x070(%%esi)\n\t"\
		"\n\t"\
		"/*...Block 2: */\n\t"\
		"addl	$0x80,%%esi			\n\t"\
		"movl	%[__p04],%%edi		\n\t"\
		"shll	$3,%%edi			\n\t"\
		"addl	%%edi,%%eax			\n\t"\
		"addl	%%edi,%%ebx			\n\t"\
		"addl	%%edi,%%ecx			\n\t"\
		"addl	%%edi,%%edx			\n\t"\
		"\n\t"\
		"/* SSE2_RADIX4_DIT_0TWIDDLE(add0, add1, add2, add3, r10): */\n\t"\
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
		"/* SSE2_RADIX4_DIT_0TWIDDLE_2NDOFTWO(add4, add5, add6, add7, r18): */\n\t"\
		"addl	$0x80,%%esi			\n\t"\
		"movl	%[__p04],%%edi		\n\t"\
		"shll	$3,%%edi			\n\t"\
		"addl	%%edi,%%eax			\n\t"\
		"addl	%%edi,%%ebx			\n\t"\
		"addl	%%edi,%%ecx			\n\t"\
		"addl	%%edi,%%edx			\n\t"\
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
		"movaps	%%xmm0,0x040(%%esi)	\n\t"\
		"movaps	%%xmm1,0x050(%%esi)	\n\t"\
		"subpd	%%xmm7,%%xmm2		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"subpd	%%xmm6,%%xmm3		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"movaps	%%xmm4,     (%%esi)	\n\t"\
		"addpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	%%xmm5,0x010(%%esi)	\n\t"\
		"\n\t"\
		"movaps	(%%edi),%%xmm5	/* ISRT2 */\n\t"\
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
		"movaps	%%xmm3,0x030(%%esi)	\n\t"\
		"movaps	%%xmm6,0x070(%%esi)	\n\t"\
		"movaps	%%xmm0,0x020(%%esi)	\n\t"\
		"movaps	%%xmm1,0x060(%%esi)	\n\t"\
		"/* SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS(r10,r12,r14,r16,r18,r1A,r1C,r1E): */\n\t"\
		"movaps		-0x080(%%esi),%%xmm0	/* r10 */\n\t"\
		"movaps		-0x040(%%esi),%%xmm4	/* r14 */\n\t"\
		"movaps		-0x070(%%esi),%%xmm1	/* r11 */\n\t"\
		"movaps		-0x030(%%esi),%%xmm5	/* r15 */\n\t"\
		"movaps		      (%%esi),%%xmm2	/* r18 */\n\t"\
		"movaps		 0x050(%%esi),%%xmm7	/* r1D */\n\t"\
		"movaps		 0x010(%%esi),%%xmm3	/* r19 */\n\t"\
		"movaps		 0x040(%%esi),%%xmm6	/* r1C */\n\t"\
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
		"movaps		%%xmm0,      (%%esi)\n\t"\
		"movaps		%%xmm4, 0x040(%%esi)\n\t"\
		"movaps		%%xmm1, 0x010(%%esi)\n\t"\
		"movaps		%%xmm5,-0x030(%%esi)\n\t"\
		"movaps		%%xmm2,-0x080(%%esi)\n\t"\
		"movaps		%%xmm7,-0x040(%%esi)\n\t"\
		"movaps		%%xmm3,-0x070(%%esi)\n\t"\
		"movaps		%%xmm6, 0x050(%%esi)\n\t"\
		"\n\t"\
		"movaps		-0x060(%%esi),%%xmm0	/* r12 */\n\t"\
		"movaps		-0x020(%%esi),%%xmm4	/* r16 */\n\t"\
		"movaps		-0x050(%%esi),%%xmm1	/* r13 */\n\t"\
		"movaps		-0x010(%%esi),%%xmm5	/* r17 */\n\t"\
		"movaps		 0x020(%%esi),%%xmm2	/* r1A */\n\t"\
		"movaps		 0x070(%%esi),%%xmm7	/* r1F */\n\t"\
		"movaps		 0x030(%%esi),%%xmm3	/* r1B */\n\t"\
		"movaps		 0x060(%%esi),%%xmm6	/* r1E */\n\t"\
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
		"movaps		%%xmm0, 0x020(%%esi)\n\t"\
		"movaps		%%xmm4, 0x060(%%esi)\n\t"\
		"movaps		%%xmm1, 0x030(%%esi)\n\t"\
		"movaps		%%xmm5,-0x010(%%esi)\n\t"\
		"movaps		%%xmm2,-0x060(%%esi)\n\t"\
		"movaps		%%xmm7,-0x020(%%esi)\n\t"\
		"movaps		%%xmm3,-0x050(%%esi)\n\t"\
		"movaps		%%xmm6, 0x070(%%esi)\n\t"\
		"\n\t"\
		"/*...Block 3: */\n\t"\
		"addl	$0x80,%%esi			\n\t"\
		"movl	%[__p04],%%edi		\n\t"\
		"shll	$3,%%edi			\n\t"\
		"addl	%%edi,%%eax			\n\t"\
		"addl	%%edi,%%ebx			\n\t"\
		"addl	%%edi,%%ecx			\n\t"\
		"addl	%%edi,%%edx			\n\t"\
		"\n\t"\
		"/* SSE2_RADIX4_DIT_0TWIDDLE(add0, add1, add2, add3, r20): */\n\t"\
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
		"/* SSE2_RADIX4_DIT_0TWIDDLE_2NDOFTWO(add4, add5, add6, add7, r28): */\n\t"\
		"addl	$0x80,%%esi			\n\t"\
		"movl	%[__p04],%%edi		\n\t"\
		"shll	$3,%%edi			\n\t"\
		"addl	%%edi,%%eax			\n\t"\
		"addl	%%edi,%%ebx			\n\t"\
		"addl	%%edi,%%ecx			\n\t"\
		"addl	%%edi,%%edx			\n\t"\
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
		"movaps	%%xmm0,0x040(%%esi)	\n\t"\
		"movaps	%%xmm1,0x050(%%esi)	\n\t"\
		"subpd	%%xmm7,%%xmm2		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"subpd	%%xmm6,%%xmm3		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"movaps	%%xmm4,     (%%esi)	\n\t"\
		"addpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	%%xmm5,0x010(%%esi)	\n\t"\
		"\n\t"\
		"movaps	(%%edi),%%xmm5	/* ISRT2 */\n\t"\
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
		"movaps	%%xmm3,0x030(%%esi)	\n\t"\
		"movaps	%%xmm6,0x070(%%esi)	\n\t"\
		"movaps	%%xmm0,0x020(%%esi)	\n\t"\
		"movaps	%%xmm1,0x060(%%esi)	\n\t"\
		"/* SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS(r20,r22,r24,r26,r28,r2A,r2C,r2E): */\n\t"\
		"movaps		-0x080(%%esi),%%xmm0	/* r20 */\n\t"\
		"movaps		-0x040(%%esi),%%xmm4	/* r24 */\n\t"\
		"movaps		-0x070(%%esi),%%xmm1	/* r21 */\n\t"\
		"movaps		-0x030(%%esi),%%xmm5	/* r25 */\n\t"\
		"movaps		      (%%esi),%%xmm2	/* r28 */\n\t"\
		"movaps		 0x050(%%esi),%%xmm7	/* r2D */\n\t"\
		"movaps		 0x010(%%esi),%%xmm3	/* r29 */\n\t"\
		"movaps		 0x040(%%esi),%%xmm6	/* r2C */\n\t"\
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
		"movaps		%%xmm0,      (%%esi)\n\t"\
		"movaps		%%xmm4, 0x040(%%esi)\n\t"\
		"movaps		%%xmm1, 0x010(%%esi)\n\t"\
		"movaps		%%xmm5,-0x030(%%esi)\n\t"\
		"movaps		%%xmm2,-0x080(%%esi)\n\t"\
		"movaps		%%xmm7,-0x040(%%esi)\n\t"\
		"movaps		%%xmm3,-0x070(%%esi)\n\t"\
		"movaps		%%xmm6, 0x050(%%esi)\n\t"\
		"\n\t"\
		"movaps		-0x060(%%esi),%%xmm0	/* r22 */\n\t"\
		"movaps		-0x020(%%esi),%%xmm4	/* r26 */\n\t"\
		"movaps		-0x050(%%esi),%%xmm1	/* r23 */\n\t"\
		"movaps		-0x010(%%esi),%%xmm5	/* r27 */\n\t"\
		"movaps		 0x020(%%esi),%%xmm2	/* r2A */\n\t"\
		"movaps		 0x070(%%esi),%%xmm7	/* r2F */\n\t"\
		"movaps		 0x030(%%esi),%%xmm3	/* r2B */\n\t"\
		"movaps		 0x060(%%esi),%%xmm6	/* r2E */\n\t"\
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
		"movaps		%%xmm0, 0x020(%%esi)\n\t"\
		"movaps		%%xmm4, 0x060(%%esi)\n\t"\
		"movaps		%%xmm1, 0x030(%%esi)\n\t"\
		"movaps		%%xmm5,-0x010(%%esi)\n\t"\
		"movaps		%%xmm2,-0x060(%%esi)\n\t"\
		"movaps		%%xmm7,-0x020(%%esi)\n\t"\
		"movaps		%%xmm3,-0x050(%%esi)\n\t"\
		"movaps		%%xmm6, 0x070(%%esi)\n\t"\
		"\n\t"\
		"/*...Block 4: */\n\t"\
		"addl	$0x80,%%esi			\n\t"\
		"movl	%[__p04],%%edi		\n\t"\
		"shll	$3,%%edi			\n\t"\
		"addl	%%edi,%%eax			\n\t"\
		"addl	%%edi,%%ebx			\n\t"\
		"addl	%%edi,%%ecx			\n\t"\
		"addl	%%edi,%%edx			\n\t"\
		"\n\t"\
		"/* SSE2_RADIX4_DIT_0TWIDDLE(add0, add1, add2, add3, r30): */\n\t"\
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
		"/* SSE2_RADIX4_DIT_0TWIDDLE_2NDOFTWO(add4, add5, add6, add7, r38): */\n\t"\
		"addl	$0x80,%%esi			\n\t"\
		"movl	%[__p04],%%edi		\n\t"\
		"shll	$3,%%edi			\n\t"\
		"addl	%%edi,%%eax			\n\t"\
		"addl	%%edi,%%ebx			\n\t"\
		"addl	%%edi,%%ecx			\n\t"\
		"addl	%%edi,%%edx			\n\t"\
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
		"movaps	%%xmm0,0x040(%%esi)	\n\t"\
		"movaps	%%xmm1,0x050(%%esi)	\n\t"\
		"subpd	%%xmm7,%%xmm2		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"subpd	%%xmm6,%%xmm3		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"movaps	%%xmm4,     (%%esi)	\n\t"\
		"addpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	%%xmm5,0x010(%%esi)	\n\t"\
		"\n\t"\
		"movaps	(%%edi),%%xmm5	/* ISRT2 */\n\t"\
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
		"movaps	%%xmm3,0x030(%%esi)	\n\t"\
		"movaps	%%xmm6,0x070(%%esi)	\n\t"\
		"movaps	%%xmm0,0x020(%%esi)	\n\t"\
		"movaps	%%xmm1,0x060(%%esi)	\n\t"\
		"/* SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS(r30,r32,r34,r36,r38,r3A,r3C,r3E): */\n\t"\
		"movaps		-0x080(%%esi),%%xmm0	/* r30 */\n\t"\
		"movaps		-0x040(%%esi),%%xmm4	/* r34 */\n\t"\
		"movaps		-0x070(%%esi),%%xmm1	/* r31 */\n\t"\
		"movaps		-0x030(%%esi),%%xmm5	/* r35 */\n\t"\
		"movaps		      (%%esi),%%xmm2	/* r38 */\n\t"\
		"movaps		 0x050(%%esi),%%xmm7	/* r3D */\n\t"\
		"movaps		 0x010(%%esi),%%xmm3	/* r39 */\n\t"\
		"movaps		 0x040(%%esi),%%xmm6	/* r3C */\n\t"\
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
		"movaps		%%xmm0,      (%%esi)\n\t"\
		"movaps		%%xmm4, 0x040(%%esi)\n\t"\
		"movaps		%%xmm1, 0x010(%%esi)\n\t"\
		"movaps		%%xmm5,-0x030(%%esi)\n\t"\
		"movaps		%%xmm2,-0x080(%%esi)\n\t"\
		"movaps		%%xmm7,-0x040(%%esi)\n\t"\
		"movaps		%%xmm3,-0x070(%%esi)\n\t"\
		"movaps		%%xmm6, 0x050(%%esi)\n\t"\
		"\n\t"\
		"movaps		-0x060(%%esi),%%xmm0	/* r32 */\n\t"\
		"movaps		-0x020(%%esi),%%xmm4	/* r36 */\n\t"\
		"movaps		-0x050(%%esi),%%xmm1	/* r33 */\n\t"\
		"movaps		-0x010(%%esi),%%xmm5	/* r37 */\n\t"\
		"movaps		 0x020(%%esi),%%xmm2	/* r3A */\n\t"\
		"movaps		 0x070(%%esi),%%xmm7	/* r3F */\n\t"\
		"movaps		 0x030(%%esi),%%xmm3	/* r3B */\n\t"\
		"movaps		 0x060(%%esi),%%xmm6	/* r3E */\n\t"\
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
		"movaps		%%xmm0, 0x020(%%esi)\n\t"\
		"movaps		%%xmm4, 0x060(%%esi)\n\t"\
		"movaps		%%xmm1, 0x030(%%esi)\n\t"\
		"movaps		%%xmm5,-0x010(%%esi)\n\t"\
		"movaps		%%xmm2,-0x060(%%esi)\n\t"\
		"movaps		%%xmm7,-0x020(%%esi)\n\t"\
		"movaps		%%xmm3,-0x050(%%esi)\n\t"\
		"movaps		%%xmm6, 0x070(%%esi)\n\t"\
		"\n\t"\
	"/*...and now do eight radix-4 transforms, including the internal twiddle factors:	*/\n\t"\
		"/*...Block 1: r00,r10,r20,r30	*/\n\t"\
		"\n\t"\
		"movl	%[__isrt2],%%edi	\n\t"\
		"movl	%[__cc0],%%esi	\n\t"\
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
		"\n\t"\
		"subpd	    (%%ebx),%%xmm0		\n\t"\
		"subpd	0x10(%%ebx),%%xmm1		\n\t"\
		"addpd	    (%%eax),%%xmm2		\n\t"\
		"addpd	0x10(%%eax),%%xmm3		\n\t"\
		"\n\t"\
		"movaps	    (%%ecx),%%xmm4		\n\t"\
		"movaps	0x10(%%ecx),%%xmm5		\n\t"\
		"movaps	    (%%edx),%%xmm6		\n\t"\
		"movaps	0x10(%%edx),%%xmm7		\n\t"\
		"\n\t"\
		"subpd	    (%%edx),%%xmm4		\n\t"\
		"subpd	0x10(%%edx),%%xmm5		\n\t"\
		"addpd	    (%%ecx),%%xmm6		\n\t"\
		"addpd	0x10(%%ecx),%%xmm7		\n\t"\
		"\n\t"\
		"/* SSE2_RADIX4_DIT_LASTPAIR_IN_PLACE_STRIDE(xmm2,xmm3,xmm0,xmm1,xmm6,xmm7,xmm4,xmm5): swap xmm[01]<->[23], xmm[45]<->[67] */\n\t"\
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
		"/*...Block 5: r08,r18,r28,r38	*/\n\t"\
		"\n\t"\
		"addl	$0x80,%%eax		/* r08 */\n\t"\
		"addl	$0x80,%%ebx		\n\t"\
		"addl	$0x80,%%ecx		\n\t"\
		"addl	$0x80,%%edx		\n\t"\
		"movaps	(%%edi),%%xmm2	/* isrt2 */\n\t"\
		"\n\t"\
		"movaps	    (%%ecx),%%xmm4	\n\t"\
		"movaps	0x10(%%ecx),%%xmm5	\n\t"\
		"movaps	    (%%edx),%%xmm0	\n\t"\
		"movaps	0x10(%%edx),%%xmm1	\n\t"\
		"\n\t"\
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
		"\n\t"\
		"subpd	0x10(%%ebx),%%xmm0	\n\t"\
		"subpd	    (%%ebx),%%xmm1	\n\t"\
		"addpd	    (%%eax),%%xmm3	\n\t"\
		"addpd	0x10(%%eax),%%xmm2	\n\t"\
		"\n\t"\
		"/* SSE2_RADIX4_DIT_LASTPAIR_IN_PLACE_STRIDE(xmm3,xmm1,xmm0,xmm2,xmm4,xmm5,xmm6,xmm7): swap xmm0123<->3102 */\n\t"\
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
		"/*...Block 3: r04,r14,r24,r34	*/\n\t"\
		"\n\t"\
		"subl	$0x40,%%eax 	/* r04 */\n\t"\
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
		"mulpd	(%%edi),%%xmm2		\n\t"\
		"mulpd	(%%edi),%%xmm3		\n\t"\
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
		"/*...Block 7: r0C,r1C,r2C,r3C	*/\n\t"\
		"\n\t"\
		"addl	$0x80,%%eax 	/* r0C */\n\t"\
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
		"mulpd	(%%edi),%%xmm2		\n\t"\
		"mulpd	(%%edi),%%xmm3		\n\t"\
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
		"/*...Block 2: r02,r12,r22,r32	*/\n\t"\
		"\n\t"\
		"subl	$0xa0,%%eax 	/* r02 */\n\t"\
		"subl	$0xa0,%%ebx			\n\t"\
		"subl	$0xa0,%%ecx			\n\t"\
		"subl	$0xa0,%%edx			\n\t"\
		"addl	$0x30,%%edi /* cc1 */\n\t"\
		"addl	$0x40,%%esi /* cc3 */\n\t"\
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
		"subl	$0x40,%%esi /* cc0 */\n\t"\
		"movaps	    (%%ebx),%%xmm2	\n\t"\
		"movaps	0x10(%%ebx),%%xmm3	\n\t"\
		"movaps	    (%%ebx),%%xmm0	\n\t"\
		"movaps	0x10(%%ebx),%%xmm1	\n\t"\
		"\n\t"\
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
		"/* SSE2_RADIX4_DIT_LASTPAIR_IN_PLACE_STRIDE(xmm2,xmm3,xmm0,xmm1,xmm6,xmm7,xmm4,xmm5): swap xmm[01]<->[23], xmm[45]<->[67] */\n\t"\
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
		"/*...Block 6: r0A,r1A,r2A,r3A	*/\n\t"\
		"\n\t"\
		"addl	$0x80,%%eax 	/* r0A */\n\t"\
		"addl	$0x80,%%ebx			\n\t"\
		"addl	$0x80,%%ecx			\n\t"\
		"addl	$0x80,%%edx			\n\t"\
		"addl	$0x40,%%esi /* cc3 */\n\t"\
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
		"subl	$0x40,%%esi /* cc0 */\n\t"\
		"movaps	    (%%ebx),%%xmm2	\n\t"\
		"movaps	0x10(%%ebx),%%xmm3	\n\t"\
		"movaps	    (%%ebx),%%xmm0	\n\t"\
		"movaps	0x10(%%ebx),%%xmm1	\n\t"\
		"\n\t"\
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
		"/* SSE2_RADIX4_DIT_LASTPAIR_IN_PLACE_STRIDE(xmm0,xmm1,xmm2,xmm3,xmm6,xmm7,xmm4,xmm5): swap xmm[45]<->[67] */\n\t"\
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
		"/*...Block 4: r06,r16,r26,r36	*/\n\t"\
		"\n\t"\
		"subl	$0x40,%%eax 	/* r06 */\n\t"\
		"subl	$0x40,%%ebx			\n\t"\
		"subl	$0x40,%%ecx			\n\t"\
		"subl	$0x40,%%edx			\n\t"\
		"addl	$0x40,%%esi /* cc3 */\n\t"\
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
		"subl	$0x40,%%esi /* cc0 */\n\t"\
		"movaps	    (%%ebx),%%xmm2	\n\t"\
		"movaps	0x10(%%ebx),%%xmm3	\n\t"\
		"movaps	    (%%ebx),%%xmm0	\n\t"\
		"movaps	0x10(%%ebx),%%xmm1	\n\t"\
		"\n\t"\
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
		"/* SSE2_RADIX4_DIT_LASTPAIR_IN_PLACE_STRIDE(xmm2,xmm3,xmm0,xmm1,xmm6,xmm7,xmm4,xmm5): swap xmm[01]<->[23], xmm[45]<->[67] */\n\t"\
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
		"/*...Block 8: r0E,r1E,r2E,r3E	*/\n\t"\
		"\n\t"\
		"addl	$0x80,%%eax 	/* r0E */\n\t"\
		"addl	$0x80,%%ebx		\n\t"\
		"addl	$0x80,%%ecx		\n\t"\
		"addl	$0x80,%%edx		\n\t"\
		"addl	$0x40,%%esi /* cc3 */\n\t"\
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
		"subl	$0x40,%%esi /* cc0 */\n\t"\
		"movaps	    (%%ebx),%%xmm2	\n\t"\
		"movaps	0x10(%%ebx),%%xmm3	\n\t"\
		"movaps	    (%%ebx),%%xmm0	\n\t"\
		"movaps	0x10(%%ebx),%%xmm1	\n\t"\
		"\n\t"\
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
		"/* SSE2_RADIX4_DIT_LASTPAIR_IN_PLACE_STRIDE(xmm0,xmm1,xmm2,xmm3,xmm6,xmm7,xmm4,xmm5): swap xmm[45]<->[67] */\n\t"\
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
		: "eax","ebx","ecx","edx","edi","esi"		/* Clobbered registers */\
	);\
	}


	#define SSE2_RADIX32_DIF_NOTWIDDLE(Xadd,Xp01,Xp02,Xp03,Xp04,Xp08,Xp10,Xp18,Xr00,Xisrt2,Xcc0)\
	{\
	__asm__ volatile (\
		"/* SSE2_RADIX4_DIF_IN_PLACE(r00,r20,r10,r30): */\n\t"\
		"\n\t"\
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
		"\n\t"\
		"addpd	    (%%ebx),%%xmm0	\n\t"\
		"addpd	0x10(%%ebx),%%xmm1	\n\t"\
		"subpd	    (%%ebx),%%xmm2	\n\t"\
		"subpd	0x10(%%ebx),%%xmm3	\n\t"\
		"\n\t"\
		"movaps	    (%%ecx),%%xmm4	\n\t"\
		"movaps	0x10(%%ecx),%%xmm5	\n\t"\
		"movaps	    (%%ecx),%%xmm6	\n\t"\
		"movaps	0x10(%%ecx),%%xmm7	\n\t"\
		"\n\t"\
		"addpd	    (%%edx),%%xmm4	\n\t"\
		"addpd	0x10(%%edx),%%xmm5	\n\t"\
		"subpd	    (%%edx),%%xmm6	\n\t"\
		"subpd	0x10(%%edx),%%xmm7	\n\t"\
		"/* Finish radix-4 butterfly and store results into temporary-array slots: */\n\t"\
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
		"\n\t"\
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
		"\n\t"\
		"/* SSE2_RADIX4_DIF_IN_PLACE_2NDOFTWO(r08,r28,r18,r38): */\n\t"\
		"addl	$0x080,%%eax		\n\t"\
		"addl	$0x080,%%ebx		\n\t"\
		"addl	$0x080,%%ecx		\n\t"\
		"addl	$0x080,%%edx		\n\t"\
		"movaps	     (%%eax),%%xmm0	\n\t"\
		"movaps	0x010(%%eax),%%xmm1	\n\t"\
		"movaps	     (%%eax),%%xmm2	\n\t"\
		"movaps	0x010(%%eax),%%xmm3	\n\t"\
		"addpd	     (%%ebx),%%xmm0	\n\t"\
		"addpd	0x010(%%ebx),%%xmm1	\n\t"\
		"subpd	     (%%ebx),%%xmm2	\n\t"\
		"subpd	0x010(%%ebx),%%xmm3	\n\t"\
		"movaps	     (%%ecx),%%xmm4	\n\t"\
		"movaps	0x010(%%ecx),%%xmm5	\n\t"\
		"movaps	     (%%ecx),%%xmm6	\n\t"\
		"movaps	0x010(%%ecx),%%xmm7	\n\t"\
		"addpd	     (%%edx),%%xmm6	\n\t"\
		"addpd	0x010(%%edx),%%xmm7	\n\t"\
		"subpd	     (%%edx),%%xmm4	\n\t"\
		"subpd	0x010(%%edx),%%xmm5	\n\t"\
		"/* Finish radix-4 butterfly and store results into temporary-array slots: */\n\t"\
		"subpd	%%xmm6,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm1		\n\t"\
		"subpd	%%xmm5,%%xmm2		\n\t"\
		"movaps	%%xmm0,     (%%ebx)	\n\t"\
		"movaps	%%xmm1,0x010(%%ebx)	\n\t"\
		"subpd	%%xmm4,%%xmm3		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm0,%%xmm6		\n\t"\
		"addpd	%%xmm2,%%xmm5		\n\t"\
		"addpd	%%xmm1,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm4		\n\t"\
		"movaps	%%xmm6,     (%%eax)	\n\t"\
		"movaps	%%xmm7,0x010(%%eax)	\n\t"\
		"\n\t"\
		"movaps	(%%edi),%%xmm0	/* isrt2 */\n\t"\
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
		"movaps	%%xmm2,     (%%ecx)	\n\t"\
		"movaps	%%xmm5,     (%%edx)	\n\t"\
		"movaps	%%xmm6,0x010(%%ecx)	\n\t"\
		"movaps	%%xmm7,0x010(%%edx)	\n\t"\
		"/* SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS(r00,r10,r20,r30,r08,r18,r28,r38): */\n\t"\
		"subl	$0x080,%%eax	/* r00 */\n\t"\
		"subl	$0x200,%%ebx	/* r08 */\n\t"\
		"addl	$0x080,%%ecx	/* r20 */\n\t"\
		"subl	$0x100,%%edx	/* r28 */\n\t"\
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
		"addl	$0x100,%%eax	/* r10 */	\n\t"\
		"addl	$0x100,%%ebx	/* r18 */	\n\t"\
		"addl	$0x100,%%ecx	/* r30 */	\n\t"\
		"addl	$0x100,%%edx	/* r38 */	\n\t"\
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
		"/* SSE2_RADIX4_DIF_IN_PLACE(r04,r24,r14,r34): */\n\t"\
		"\n\t"\
		"subl	$0x0C0,%%eax	/* r04 */	\n\t"\
		"addl	$0x0C0,%%ebx	/* r24 */	\n\t"\
		"subl	$0x1C0,%%ecx	/* r14 */	\n\t"\
		"subl	$0x040,%%edx	/* r34 */	\n\t"\
		"\n\t"\
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"movaps	    (%%eax),%%xmm2	\n\t"\
		"movaps	0x10(%%eax),%%xmm3	\n\t"\
		"\n\t"\
		"addpd	    (%%ebx),%%xmm0	\n\t"\
		"addpd	0x10(%%ebx),%%xmm1	\n\t"\
		"subpd	    (%%ebx),%%xmm2	\n\t"\
		"subpd	0x10(%%ebx),%%xmm3	\n\t"\
		"\n\t"\
		"movaps	    (%%ecx),%%xmm4	\n\t"\
		"movaps	0x10(%%ecx),%%xmm5	\n\t"\
		"movaps	    (%%ecx),%%xmm6	\n\t"\
		"movaps	0x10(%%ecx),%%xmm7	\n\t"\
		"\n\t"\
		"addpd	    (%%edx),%%xmm4	\n\t"\
		"addpd	0x10(%%edx),%%xmm5	\n\t"\
		"subpd	    (%%edx),%%xmm6	\n\t"\
		"subpd	0x10(%%edx),%%xmm7	\n\t"\
		"/* Finish radix-4 butterfly and store results into temporary-array slots: */\n\t"\
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
		"\n\t"\
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
		"\n\t"\
		"/* SSE2_RADIX4_DIF_IN_PLACE_2NDOFTWO(r0C,r2C,r1C,r3C): */\n\t"\
		"addl	$0x080,%%eax		\n\t"\
		"addl	$0x080,%%ebx		\n\t"\
		"addl	$0x080,%%ecx		\n\t"\
		"addl	$0x080,%%edx		\n\t"\
		"movaps	     (%%eax),%%xmm0	\n\t"\
		"movaps	0x010(%%eax),%%xmm1	\n\t"\
		"movaps	     (%%eax),%%xmm2	\n\t"\
		"movaps	0x010(%%eax),%%xmm3	\n\t"\
		"addpd	     (%%ebx),%%xmm0	\n\t"\
		"addpd	0x010(%%ebx),%%xmm1	\n\t"\
		"subpd	     (%%ebx),%%xmm2	\n\t"\
		"subpd	0x010(%%ebx),%%xmm3	\n\t"\
		"movaps	     (%%ecx),%%xmm4	\n\t"\
		"movaps	0x010(%%ecx),%%xmm5	\n\t"\
		"movaps	     (%%ecx),%%xmm6	\n\t"\
		"movaps	0x010(%%ecx),%%xmm7	\n\t"\
		"addpd	     (%%edx),%%xmm6	\n\t"\
		"addpd	0x010(%%edx),%%xmm7	\n\t"\
		"subpd	     (%%edx),%%xmm4	\n\t"\
		"subpd	0x010(%%edx),%%xmm5	\n\t"\
		"/* Finish radix-4 butterfly and store results into temporary-array slots: */\n\t"\
		"subpd	%%xmm6,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm1		\n\t"\
		"subpd	%%xmm5,%%xmm2		\n\t"\
		"movaps	%%xmm0,     (%%ebx)	\n\t"\
		"movaps	%%xmm1,0x010(%%ebx)	\n\t"\
		"subpd	%%xmm4,%%xmm3		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm0,%%xmm6		\n\t"\
		"addpd	%%xmm2,%%xmm5		\n\t"\
		"addpd	%%xmm1,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm4		\n\t"\
		"movaps	%%xmm6,     (%%eax)	\n\t"\
		"movaps	%%xmm7,0x010(%%eax)	\n\t"\
		"\n\t"\
		"movaps	(%%edi),%%xmm0	/* isrt2 */\n\t"\
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
		"movaps	%%xmm2,     (%%ecx)	\n\t"\
		"movaps	%%xmm5,     (%%edx)	\n\t"\
		"movaps	%%xmm6,0x010(%%ecx)	\n\t"\
		"movaps	%%xmm7,0x010(%%edx)	\n\t"\
		"/* SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS(r04,r14,r24,r34,r0C,r1C,r2C,r3C): */\n\t"\
		"subl	$0x080,%%eax	/* r04 */\n\t"\
		"subl	$0x200,%%ebx	/* r0C */\n\t"\
		"addl	$0x080,%%ecx	/* r24 */\n\t"\
		"subl	$0x100,%%edx	/* r2C */\n\t"\
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
		"addl	$0x100,%%eax	/* r14 */	\n\t"\
		"addl	$0x100,%%ebx	/* r1C */	\n\t"\
		"addl	$0x100,%%ecx	/* r34 */	\n\t"\
		"addl	$0x100,%%edx	/* r3C */	\n\t"\
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
		"/* SSE2_RADIX4_DIF_IN_PLACE(r02,r22,r12,r32): */\n\t"\
		"\n\t"\
		"subl	$0x120,%%eax	/* r02 */	\n\t"\
		"addl	$0x060,%%ebx	/* r22 */	\n\t"\
		"subl	$0x220,%%ecx	/* r12 */	\n\t"\
		"subl	$0x0a0,%%edx	/* r32 */	\n\t"\
		"\n\t"\
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"movaps	    (%%eax),%%xmm2	\n\t"\
		"movaps	0x10(%%eax),%%xmm3	\n\t"\
		"\n\t"\
		"addpd	    (%%ebx),%%xmm0	\n\t"\
		"addpd	0x10(%%ebx),%%xmm1	\n\t"\
		"subpd	    (%%ebx),%%xmm2	\n\t"\
		"subpd	0x10(%%ebx),%%xmm3	\n\t"\
		"\n\t"\
		"movaps	    (%%ecx),%%xmm4	\n\t"\
		"movaps	0x10(%%ecx),%%xmm5	\n\t"\
		"movaps	    (%%ecx),%%xmm6	\n\t"\
		"movaps	0x10(%%ecx),%%xmm7	\n\t"\
		"\n\t"\
		"addpd	    (%%edx),%%xmm4	\n\t"\
		"addpd	0x10(%%edx),%%xmm5	\n\t"\
		"subpd	    (%%edx),%%xmm6	\n\t"\
		"subpd	0x10(%%edx),%%xmm7	\n\t"\
		"/* Finish radix-4 butterfly and store results into temporary-array slots: */\n\t"\
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
		"\n\t"\
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
		"\n\t"\
		"/* SSE2_RADIX4_DIF_IN_PLACE_2NDOFTWO(r0A,r2A,r1A,r3A): */\n\t"\
		"addl	$0x080,%%eax		\n\t"\
		"addl	$0x080,%%ebx		\n\t"\
		"addl	$0x080,%%ecx		\n\t"\
		"addl	$0x080,%%edx		\n\t"\
		"movaps	     (%%eax),%%xmm0	\n\t"\
		"movaps	0x010(%%eax),%%xmm1	\n\t"\
		"movaps	     (%%eax),%%xmm2	\n\t"\
		"movaps	0x010(%%eax),%%xmm3	\n\t"\
		"addpd	     (%%ebx),%%xmm0	\n\t"\
		"addpd	0x010(%%ebx),%%xmm1	\n\t"\
		"subpd	     (%%ebx),%%xmm2	\n\t"\
		"subpd	0x010(%%ebx),%%xmm3	\n\t"\
		"movaps	     (%%ecx),%%xmm4	\n\t"\
		"movaps	0x010(%%ecx),%%xmm5	\n\t"\
		"movaps	     (%%ecx),%%xmm6	\n\t"\
		"movaps	0x010(%%ecx),%%xmm7	\n\t"\
		"addpd	     (%%edx),%%xmm6	\n\t"\
		"addpd	0x010(%%edx),%%xmm7	\n\t"\
		"subpd	     (%%edx),%%xmm4	\n\t"\
		"subpd	0x010(%%edx),%%xmm5	\n\t"\
		"/* Finish radix-4 butterfly and store results into temporary-array slots: */\n\t"\
		"subpd	%%xmm6,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm1		\n\t"\
		"subpd	%%xmm5,%%xmm2		\n\t"\
		"movaps	%%xmm0,     (%%ebx)	\n\t"\
		"movaps	%%xmm1,0x010(%%ebx)	\n\t"\
		"subpd	%%xmm4,%%xmm3		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm0,%%xmm6		\n\t"\
		"addpd	%%xmm2,%%xmm5		\n\t"\
		"addpd	%%xmm1,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm4		\n\t"\
		"movaps	%%xmm6,     (%%eax)	\n\t"\
		"movaps	%%xmm7,0x010(%%eax)	\n\t"\
		"\n\t"\
		"movaps	(%%edi),%%xmm0	/* isrt2 */\n\t"\
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
		"movaps	%%xmm2,     (%%ecx)	\n\t"\
		"movaps	%%xmm5,     (%%edx)	\n\t"\
		"movaps	%%xmm6,0x010(%%ecx)	\n\t"\
		"movaps	%%xmm7,0x010(%%edx)	\n\t"\
		"/* SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS(r02,r12,r22,r32,r0A,r1A,r2A,r3A): */\n\t"\
		"subl	$0x080,%%eax	/* r02 */\n\t"\
		"subl	$0x200,%%ebx	/* r0A */\n\t"\
		"addl	$0x080,%%ecx	/* r22 */\n\t"\
		"subl	$0x100,%%edx	/* r2A */\n\t"\
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
		"addl	$0x100,%%eax	/* r12 */	\n\t"\
		"addl	$0x100,%%ebx	/* r1A */	\n\t"\
		"addl	$0x100,%%ecx	/* r32 */	\n\t"\
		"addl	$0x100,%%edx	/* r3A */	\n\t"\
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
		"/* SSE2_RADIX4_DIF_IN_PLACE(r06,r26,r16,r36): */\n\t"\
		"\n\t"\
		"subl	$0x0C0,%%eax	/* r06 */	\n\t"\
		"addl	$0x0C0,%%ebx	/* r26 */	\n\t"\
		"subl	$0x1C0,%%ecx	/* r16 */	\n\t"\
		"subl	$0x040,%%edx	/* r36 */	\n\t"\
		"\n\t"\
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"movaps	    (%%eax),%%xmm2	\n\t"\
		"movaps	0x10(%%eax),%%xmm3	\n\t"\
		"\n\t"\
		"addpd	    (%%ebx),%%xmm0	\n\t"\
		"addpd	0x10(%%ebx),%%xmm1	\n\t"\
		"subpd	    (%%ebx),%%xmm2	\n\t"\
		"subpd	0x10(%%ebx),%%xmm3	\n\t"\
		"\n\t"\
		"movaps	    (%%ecx),%%xmm4	\n\t"\
		"movaps	0x10(%%ecx),%%xmm5	\n\t"\
		"movaps	    (%%ecx),%%xmm6	\n\t"\
		"movaps	0x10(%%ecx),%%xmm7	\n\t"\
		"\n\t"\
		"addpd	    (%%edx),%%xmm4	\n\t"\
		"addpd	0x10(%%edx),%%xmm5	\n\t"\
		"subpd	    (%%edx),%%xmm6	\n\t"\
		"subpd	0x10(%%edx),%%xmm7	\n\t"\
		"/* Finish radix-4 butterfly and store results into temporary-array slots: */\n\t"\
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
		"\n\t"\
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
		"\n\t"\
		"/* SSE2_RADIX4_DIF_IN_PLACE_2NDOFTWO(r0E,r2E,r1E,r3E): */\n\t"\
		"addl	$0x080,%%eax		\n\t"\
		"addl	$0x080,%%ebx		\n\t"\
		"addl	$0x080,%%ecx		\n\t"\
		"addl	$0x080,%%edx		\n\t"\
		"movaps	     (%%eax),%%xmm0	\n\t"\
		"movaps	0x010(%%eax),%%xmm1	\n\t"\
		"movaps	     (%%eax),%%xmm2	\n\t"\
		"movaps	0x010(%%eax),%%xmm3	\n\t"\
		"addpd	     (%%ebx),%%xmm0	\n\t"\
		"addpd	0x010(%%ebx),%%xmm1	\n\t"\
		"subpd	     (%%ebx),%%xmm2	\n\t"\
		"subpd	0x010(%%ebx),%%xmm3	\n\t"\
		"movaps	     (%%ecx),%%xmm4	\n\t"\
		"movaps	0x010(%%ecx),%%xmm5	\n\t"\
		"movaps	     (%%ecx),%%xmm6	\n\t"\
		"movaps	0x010(%%ecx),%%xmm7	\n\t"\
		"addpd	     (%%edx),%%xmm6	\n\t"\
		"addpd	0x010(%%edx),%%xmm7	\n\t"\
		"subpd	     (%%edx),%%xmm4	\n\t"\
		"subpd	0x010(%%edx),%%xmm5	\n\t"\
		"/* Finish radix-4 butterfly and store results into temporary-array slots: */\n\t"\
		"subpd	%%xmm6,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm1		\n\t"\
		"subpd	%%xmm5,%%xmm2		\n\t"\
		"movaps	%%xmm0,     (%%ebx)	\n\t"\
		"movaps	%%xmm1,0x010(%%ebx)	\n\t"\
		"subpd	%%xmm4,%%xmm3		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm0,%%xmm6		\n\t"\
		"addpd	%%xmm2,%%xmm5		\n\t"\
		"addpd	%%xmm1,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm4		\n\t"\
		"movaps	%%xmm6,     (%%eax)	\n\t"\
		"movaps	%%xmm7,0x010(%%eax)	\n\t"\
		"\n\t"\
		"movaps	(%%edi),%%xmm0	/* isrt2 */\n\t"\
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
		"movaps	%%xmm2,     (%%ecx)	\n\t"\
		"movaps	%%xmm5,     (%%edx)	\n\t"\
		"movaps	%%xmm6,0x010(%%ecx)	\n\t"\
		"movaps	%%xmm7,0x010(%%edx)	\n\t"\
		"/* SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS(r06,r16,r26,r36,r0E,r1E,r2E,r3E): */\n\t"\
		"subl	$0x080,%%eax	/* r02 */\n\t"\
		"subl	$0x200,%%ebx	/* r0A */\n\t"\
		"addl	$0x080,%%ecx	/* r22 */\n\t"\
		"subl	$0x100,%%edx	/* r2A */\n\t"\
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
		"addl	$0x100,%%eax	/* r12 */	\n\t"\
		"addl	$0x100,%%ebx	/* r1A */	\n\t"\
		"addl	$0x100,%%ecx	/* r32 */	\n\t"\
		"addl	$0x100,%%edx	/* r3A */	\n\t"\
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
		"/*...and now do eight radix-4 transforms, including the internal twiddle factors: */\n\t"\
		"/*...Block 1: t00,t10,t20,t30 in r00,04,02,06 - note swapped middle 2 indices! */\n\t"\
		"\n\t"\
		"movl	%[__r00],%%esi	\n\t"\
		"movl	%[__add],%%eax	\n\t"\
		"movl	%[__p01],%%ebx	\n\t"\
		"movl	%[__p02],%%ecx	\n\t"\
		"movl	%[__p03],%%edx	\n\t"\
		"shll	$3,%%ebx		\n\t"\
		"shll	$3,%%ecx		\n\t"\
		"shll	$3,%%edx		\n\t"\
		"addl	%%eax,%%ebx		\n\t"\
		"addl	%%eax,%%ecx		\n\t"\
		"addl	%%eax,%%edx		\n\t"\
		"movaps	    (%%esi),%%xmm0	/* t00 */\n\t"\
		"movaps	0x20(%%esi),%%xmm4	/* t20 */\n\t"\
		"movaps	0x10(%%esi),%%xmm1	/* t01 */\n\t"\
		"movaps	0x30(%%esi),%%xmm5	/* t21 */\n\t"\
		"movaps	0x40(%%esi),%%xmm2	/* t10 */\n\t"\
		"movaps	0x60(%%esi),%%xmm6	/* t30 */\n\t"\
		"movaps	0x50(%%esi),%%xmm3	/* t11 */\n\t"\
		"movaps	0x70(%%esi),%%xmm7	/* t31 */\n\t"\
		"subpd	0x40(%%esi),%%xmm0	/* t10=t00-rt */\n\t"\
		"subpd	0x60(%%esi),%%xmm4	/* t30=t20-rt */\n\t"\
		"subpd	0x50(%%esi),%%xmm1	/* t11=t01-it */\n\t"\
		"subpd	0x70(%%esi),%%xmm5	/* t31=t21-it */\n\t"\
		"addpd	    (%%esi),%%xmm2	/* t00=t00+rt */\n\t"\
		"addpd	0x20(%%esi),%%xmm6	/* t20=t20+rt */\n\t"\
		"addpd	0x10(%%esi),%%xmm3	/* t01=t01+it */\n\t"\
		"addpd	0x30(%%esi),%%xmm7	/* t21=t21+it */\n\t"\
		"subpd	%%xmm6,%%xmm2		/* t00 <- t00-t20 */\n\t"\
		"subpd	%%xmm5,%%xmm0		/* t10 <- t10-t31 */\n\t"\
		"subpd	%%xmm7,%%xmm3		/* t01 <- t01-t21 */\n\t"\
		"subpd	%%xmm4,%%xmm1		/* t11 <- t11-t30 */\n\t"\
		"addpd	%%xmm6,%%xmm6		/*          2*t20 */\n\t"\
		"addpd	%%xmm5,%%xmm5		/*          2*t31 */\n\t"\
		"addpd	%%xmm7,%%xmm7		/*          2*t21 */\n\t"\
		"addpd	%%xmm4,%%xmm4		/*          2*t30 */\n\t"\
		"movaps	%%xmm2,    (%%ebx)	/* a(jt+p1 ) */\n\t"\
		"movaps	%%xmm0,    (%%ecx)	/* a(jt+p2 ) */\n\t"\
		"movaps	%%xmm3,0x10(%%ebx)	/* a(jp+p1 ) */\n\t"\
		"movaps	%%xmm1,0x10(%%edx)	/* a(jp+p3 ) */\n\t"\
		"addpd	%%xmm2,%%xmm6		/* t20 <- t00+t20 */\n\t"\
		"addpd	%%xmm0,%%xmm5		/* t31 <- t10+t31 */\n\t"\
		"addpd	%%xmm3,%%xmm7		/* t21 <- t01+t21 */\n\t"\
		"addpd	%%xmm1,%%xmm4		/* t30 <- t11+t30 */\n\t"\
		"movaps	%%xmm6,    (%%eax)	/* a(jt+p0 ) */\n\t"\
		"movaps	%%xmm5,    (%%edx)	/* a(jt+p3 ) */\n\t"\
		"movaps	%%xmm7,0x10(%%eax)	/* a(jp+p0 ) */\n\t"\
		"movaps	%%xmm4,0x10(%%ecx)	/* a(jp+p2 ) */\n\t"\
		"/*...Block 5: t08,t18,t28,t38	*/\n\t"\
		"addl	$0x80,%%esi			/* r08 */\n\t"\
		"movl	%[__p04],%%edi		\n\t"\
		"shll	$3,%%edi			\n\t"\
		"addl	%%edi,%%eax			/* ap04 */\n\t"\
		"addl	%%edi,%%ebx			\n\t"\
		"addl	%%edi,%%ecx			\n\t"\
		"addl	%%edi,%%edx			\n\t"\
		"movl	%[__isrt2],%%edi	\n\t"\
		"\n\t"\
		"movaps	0x20(%%esi),%%xmm4	/* t28 */\n\t"\
		"movaps	0x30(%%esi),%%xmm5	/* t29 */\n\t"\
		"movaps	0x60(%%esi),%%xmm6	/* t38 */\n\t"\
		"movaps	0x70(%%esi),%%xmm7	/* t39 */\n\t"\
		"movaps	    (%%esi),%%xmm0	/* t08 */\n\t"\
		"subpd	0x30(%%esi),%%xmm4	/* t28-t29 */\n\t"\
		"movaps	0x10(%%esi),%%xmm1	/* t09 */\n\t"\
		"addpd	0x20(%%esi),%%xmm5	/* t29+t28 */\n\t"\
		"movaps	0x40(%%esi),%%xmm2	/* t18 */\n\t"\
		"mulpd	(%%edi),%%xmm4		/* t28 = (t28-t29)*ISRT2 */\n\t"\
		"movaps	0x50(%%esi),%%xmm3	/* t19 */\n\t"\
		"mulpd	(%%edi),%%xmm5		/* t29 = (t29+t28)*ISRT2 */\n\t"\
		"subpd	0x50(%%esi),%%xmm0	/* t08=t08-t19*/\n\t"\
		"addpd	0x70(%%esi),%%xmm6	/* t38+t39 */\n\t"\
		"subpd	0x40(%%esi),%%xmm1	/* t19=t09-t18*/\n\t"\
		"subpd	0x60(%%esi),%%xmm7	/* t39-t38 */\n\t"\
		"addpd	0x10(%%esi),%%xmm2	/* t09=t18+t09*/\n\t"\
		"mulpd	(%%edi),%%xmm6		/*  rt = (t38+t39)*ISRT2 */\n\t"\
		"addpd	    (%%esi),%%xmm3	/* t18=t19+t08*/\n\t"\
		"mulpd	(%%edi),%%xmm7		/*  it = (t39-t38)*ISRT2 */\n\t"\
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
		"movaps	%%xmm0,    (%%ebx)	/* a(jt+p1 ) */\n\t"\
		"movaps	%%xmm2,0x10(%%ebx)	/* a(jp+p1 ) */\n\t"\
		"movaps	%%xmm3,    (%%ecx)	/* a(jt+p2 ) */\n\t"\
		"movaps	%%xmm1,0x10(%%edx)	/* a(jp+p3 ) */\n\t"\
		"addpd	%%xmm0,%%xmm4		/* t08+t28 */\n\t"\
		"addpd	%%xmm2,%%xmm5		/* t09+t29 */\n\t"\
		"addpd	%%xmm3,%%xmm7		/* t18+t39 */\n\t"\
		"addpd	%%xmm1,%%xmm6		/* t19+t38 */\n\t"\
		"movaps	%%xmm4,    (%%eax)	/* a(jt+p0 ) */\n\t"\
		"movaps	%%xmm5,0x10(%%eax)	/* a(jp+p0 ) */\n\t"\
		"movaps	%%xmm7,    (%%edx)	/* a(jt+p3 ) */\n\t"\
		"movaps	%%xmm6,0x10(%%ecx)	/* a(jp+p2 ) */\n\t"\
		"/*...Block 3: t04,t14,t24,t34	*/\n\t"\
		"addl	$0x180,%%esi		/* r20 */\n\t"\
		"movl	%[__p04],%%edi		\n\t"\
		"shll	$3,%%edi			\n\t"\
		"addl	%%edi,%%eax			/* ap08 */\n\t"\
		"addl	%%edi,%%ebx			\n\t"\
		"addl	%%edi,%%ecx			\n\t"\
		"addl	%%edi,%%edx			\n\t"\
		"movl	%[__cc0],%%edi	\n\t"\
		"\n\t"\
		"movaps	0x20(%%esi),%%xmm4	/* t24 */\n\t"\
		"movaps	0x60(%%esi),%%xmm6	/* t34 */\n\t"\
		"movaps	0x30(%%esi),%%xmm5	/* t25 */\n\t"\
		"movaps	0x70(%%esi),%%xmm7	/* t35 */\n\t"\
		"movaps	0x20(%%esi),%%xmm0	/* copy t24 */\n\t"\
		"movaps	0x60(%%esi),%%xmm2	/* copy t34 */\n\t"\
		"movaps	0x30(%%esi),%%xmm1	/* copy t25 */\n\t"\
		"movaps	0x70(%%esi),%%xmm3	/* copy t35 */\n\t"\
		"mulpd	    (%%edi),%%xmm4	/* t24*c */\n\t"\
		"mulpd	0x10(%%edi),%%xmm6	/* t34*s */\n\t"\
		"mulpd	0x10(%%edi),%%xmm1	/* t25*s */\n\t"\
		"mulpd	    (%%edi),%%xmm3	/* t35*c */\n\t"\
		"mulpd	    (%%edi),%%xmm5	/* t25*c */\n\t"\
		"mulpd	0x10(%%edi),%%xmm7	/* t35*s */\n\t"\
		"mulpd	0x10(%%edi),%%xmm0	/* t24*s */\n\t"\
		"mulpd	    (%%edi),%%xmm2	/* t34*c */\n\t"\
		"subpd	%%xmm1,%%xmm4		/* ~t24 */\n\t"\
		"subpd	%%xmm3,%%xmm6		/* rt */\n\t"\
		"addpd	%%xmm0,%%xmm5		/* ~t25 */\n\t"\
		"addpd	%%xmm2,%%xmm7		/* it */\n\t"\
		"subl	$0x10,%%edi			/* isrt2 */\n\t"\
		"\n\t"\
		"movaps	0x40(%%esi),%%xmm2	/* t14 */\n\t"\
		"subpd	%%xmm6,%%xmm4		/*~t34=t24-rt */\n\t"\
		"movaps	0x50(%%esi),%%xmm3	/* t15 */\n\t"\
		"subpd	%%xmm7,%%xmm5		/*~t35=t25-it */\n\t"\
		"subpd	0x50(%%esi),%%xmm2	/* t14-t15 */\n\t"\
		"addpd	%%xmm6,%%xmm6		/*      2* rt */\n\t"\
		"addpd	0x40(%%esi),%%xmm3	/* t15+t14 */\n\t"\
		"addpd	%%xmm7,%%xmm7		/*      2* it */\n\t"\
		"mulpd	(%%edi),%%xmm2		/* rt = (t14-t15)*ISRT2 */\n\t"\
		"addpd	%%xmm4,%%xmm6		/*~t24=t24+rt */\n\t"\
		"mulpd	(%%edi),%%xmm3		/* it = (t15+t14)*ISRT2 */\n\t"\
		"addpd	%%xmm5,%%xmm7		/*~t25=t25+it */\n\t"\
		"movaps	    (%%esi),%%xmm0	/* t04 */\n\t"\
		"movaps	0x10(%%esi),%%xmm1	/* t05 */\n\t"\
		"subpd	%%xmm2,%%xmm0		/*~t14=t04-rt */\n\t"\
		"subpd	%%xmm3,%%xmm1		/*~t15=t05-it */\n\t"\
		"addpd	    (%%esi),%%xmm2	/*~t04=rt +t04*/\n\t"\
		"addpd	0x10(%%esi),%%xmm3	/*~t05=it +t05*/\n\t"\
		"subpd	%%xmm6,%%xmm2		/* t04-t24 */\n\t"\
		"subpd	%%xmm5,%%xmm0		/* t14-t35 */\n\t"\
		"subpd	%%xmm7,%%xmm3		/* t05-t25 */\n\t"\
		"subpd	%%xmm4,%%xmm1		/* t15-t34 */\n\t"\
		"addpd	%%xmm6,%%xmm6		/*   2*t24 */\n\t"\
		"addpd	%%xmm5,%%xmm5		/*          2*t35 */\n\t"\
		"addpd	%%xmm7,%%xmm7		/*   2*t25 */\n\t"\
		"addpd	%%xmm4,%%xmm4		/*          2*t34 */\n\t"\
		"movaps	%%xmm2,    (%%ebx)	/* a(jt+p1 ) */\n\t"\
		"movaps	%%xmm0,    (%%ecx)	/* a(jt+p2 ) */\n\t"\
		"movaps	%%xmm3,0x10(%%ebx)	/* a(jp+p1 ) */\n\t"\
		"movaps	%%xmm1,0x10(%%edx)	/* a(jp+p3 ) */\n\t"\
		"addpd	%%xmm2,%%xmm6		/* t04+t24 */\n\t"\
		"addpd	%%xmm0,%%xmm5		/* t14+t35 */\n\t"\
		"addpd	%%xmm3,%%xmm7		/* t05+t25 */\n\t"\
		"addpd	%%xmm1,%%xmm4		/* t15+t34 */\n\t"\
		"movaps	%%xmm6,    (%%eax)	/* a(jt+p0 ) */\n\t"\
		"movaps	%%xmm5,    (%%edx)	/* a(jt+p3 ) */\n\t"\
		"movaps	%%xmm7,0x10(%%eax)	/* a(jp+p0 ) */\n\t"\
		"movaps	%%xmm4,0x10(%%ecx)	/* a(jp+p2 ) */\n\t"\
		"/*...Block 7: t0C,t1C,t2C,t3C	*/\n\t"\
		"addl	$0x80,%%esi			/* r28 */\n\t"\
		"movl	%[__p04],%%edi		\n\t"\
		"shll	$3,%%edi			\n\t"\
		"addl	%%edi,%%eax			/* ap0C */\n\t"\
		"addl	%%edi,%%ebx			\n\t"\
		"addl	%%edi,%%ecx			\n\t"\
		"addl	%%edi,%%edx			\n\t"\
		"movl	%[__cc0],%%edi	\n\t"\
		"\n\t"\
		"movaps	0x20(%%esi),%%xmm4	/* t2C */\n\t"\
		"movaps	0x60(%%esi),%%xmm6	/* t3C */\n\t"\
		"movaps	0x30(%%esi),%%xmm5	/* t2D */\n\t"\
		"movaps	0x70(%%esi),%%xmm7	/* t3D */\n\t"\
		"movaps	0x20(%%esi),%%xmm0	/* copy t2C */\n\t"\
		"movaps	0x60(%%esi),%%xmm2	/* copy t3C */\n\t"\
		"movaps	0x30(%%esi),%%xmm1	/* copy t2D */\n\t"\
		"movaps	0x70(%%esi),%%xmm3	/* copy t3D */\n\t"\
		"mulpd	0x10(%%edi),%%xmm4	/* t2C*s */\n\t"\
		"mulpd	    (%%edi),%%xmm6	/* t3C*c */\n\t"\
		"mulpd	    (%%edi),%%xmm1	/* t2D*c */\n\t"\
		"mulpd	0x10(%%edi),%%xmm3	/* t3D*s */\n\t"\
		"mulpd	0x10(%%edi),%%xmm5	/* t2D*s */\n\t"\
		"mulpd	    (%%edi),%%xmm7	/* t3D*c */\n\t"\
		"mulpd	    (%%edi),%%xmm0	/* t2C*c */\n\t"\
		"mulpd	0x10(%%edi),%%xmm2	/* t3C*s */\n\t"\
		"subpd	%%xmm1,%%xmm4		/* ~t24 */\n\t"\
		"subpd	%%xmm3,%%xmm6		/* rt */\n\t"\
		"addpd	%%xmm0,%%xmm5		/* ~t25 */\n\t"\
		"addpd	%%xmm2,%%xmm7		/* it */\n\t"\
		"subl	$0x10,%%edi			/* isrt2 */\n\t"\
		"\n\t"\
		"movaps	0x40(%%esi),%%xmm2	/* t14 */\n\t"\
		"subpd	%%xmm6,%%xmm4		/*~t2C=t2C-rt */\n\t"\
		"movaps	0x50(%%esi),%%xmm3	/* t1D */\n\t"\
		"subpd	%%xmm7,%%xmm5		/*~t2D=t2D-it */\n\t"\
		"addpd	0x50(%%esi),%%xmm2	/* t1C+t1D */\n\t"\
		"addpd	%%xmm6,%%xmm6		/*      2* rt */\n\t"\
		"subpd	0x40(%%esi),%%xmm3	/* t1D-t1C */\n\t"\
		"addpd	%%xmm7,%%xmm7		/*      2* it */\n\t"\
		"mulpd	(%%edi),%%xmm2		/* rt = (t1C+t1D)*ISRT2 */\n\t"\
		"addpd	%%xmm4,%%xmm6		/*~t3C=t2C+rt */\n\t"\
		"mulpd	(%%edi),%%xmm3		/* it = (t1D-t1C)*ISRT2 */\n\t"\
		"addpd	%%xmm5,%%xmm7		/*~t3D=t2D+it */\n\t"\
		"movaps	    (%%esi),%%xmm0	/* t0C */\n\t"\
		"movaps	0x10(%%esi),%%xmm1	/* t0D */\n\t"\
		"subpd	%%xmm2,%%xmm0		/*~t0C=t0C-rt */\n\t"\
		"subpd	%%xmm3,%%xmm1		/*~t0D=t0D-it */\n\t"\
		"addpd	    (%%esi),%%xmm2	/*~t1C=rt +t0C*/\n\t"\
		"addpd	0x10(%%esi),%%xmm3	/*~t1D=it +t0D*/\n\t"\
		"subpd	%%xmm4,%%xmm0		/* t0C-t2C */\n\t"\
		"subpd	%%xmm7,%%xmm2		/* t1C-t3D */\n\t"\
		"subpd	%%xmm5,%%xmm1		/* t0D-t2D */\n\t"\
		"subpd	%%xmm6,%%xmm3		/* t1D-t3C */\n\t"\
		"addpd	%%xmm4,%%xmm4		/*   2*t2C */\n\t"\
		"addpd	%%xmm7,%%xmm7		/*   2*t3D */\n\t"\
		"addpd	%%xmm5,%%xmm5		/*   2*t2D */\n\t"\
		"addpd	%%xmm6,%%xmm6		/*   2*t3C */\n\t"\
		"movaps	%%xmm0,    (%%ebx)	/* a(jt+p1 ) */\n\t"\
		"movaps	%%xmm2,    (%%ecx)	/* a(jt+p2 ) */\n\t"\
		"movaps	%%xmm1,0x10(%%ebx)	/* a(jp+p1 ) */\n\t"\
		"movaps	%%xmm3,0x10(%%edx)	/* a(jp+p3 ) */\n\t"\
		"addpd	%%xmm0,%%xmm4		/* t0C+t2C */\n\t"\
		"addpd	%%xmm2,%%xmm7		/* t1C+t3D */\n\t"\
		"addpd	%%xmm1,%%xmm5		/* t0D+t2D */\n\t"\
		"addpd	%%xmm3,%%xmm6		/* t1D+t3C */\n\t"\
		"movaps	%%xmm4,    (%%eax)	/* a(jt+p0 ) */\n\t"\
		"movaps	%%xmm7,    (%%edx)	/* a(jt+p3 ) */\n\t"\
		"movaps	%%xmm5,0x10(%%eax)	/* a(jp+p0 ) */\n\t"\
		"movaps	%%xmm6,0x10(%%ecx)	/* a(jp+p2 ) */\n\t"\
		"/*...Block 2: t02,t12,t22,t32	*/\n\t"\
		"subl	$0x180,%%esi		/* r10 */\n\t"\
		"movl	%[__p04],%%edi		\n\t"\
		"shll	$3,%%edi			\n\t"\
		"addl	%%edi,%%eax			/* ap10 */\n\t"\
		"addl	%%edi,%%ebx			\n\t"\
		"addl	%%edi,%%ecx			\n\t"\
		"addl	%%edi,%%edx			\n\t"\
		"movl	%[__cc0],%%edi	\n\t"\
		"\n\t"\
		"movaps	0x20(%%esi),%%xmm4	/* t22 */\n\t"\
		"movaps	0x60(%%esi),%%xmm6	/* t32 */\n\t"\
		"movaps	0x30(%%esi),%%xmm5	/* t23 */\n\t"\
		"movaps	0x70(%%esi),%%xmm7	/* t33 */\n\t"\
		"movaps	0x20(%%esi),%%xmm0	/* copy t22 */\n\t"\
		"movaps	0x60(%%esi),%%xmm2	/* copy t32 */\n\t"\
		"movaps	0x30(%%esi),%%xmm1	/* copy t23 */\n\t"\
		"movaps	0x70(%%esi),%%xmm3	/* copy t33 */\n\t"\
		"mulpd	0x20(%%edi),%%xmm4	/* t22*c32_1 */\n\t"\
		"mulpd	0x40(%%edi),%%xmm6	/* t32*c32_3 */\n\t"\
		"mulpd	0x30(%%edi),%%xmm1	/* t23*s32_1 */\n\t"\
		"mulpd	0x50(%%edi),%%xmm3	/* t33*s32_3 */\n\t"\
		"mulpd	0x20(%%edi),%%xmm5	/* t23*c32_1 */\n\t"\
		"mulpd	0x40(%%edi),%%xmm7	/* t33*c32_3 */\n\t"\
		"mulpd	0x30(%%edi),%%xmm0	/* t22*s32_1 */\n\t"\
		"mulpd	0x50(%%edi),%%xmm2	/* t32*s32_3 */\n\t"\
		"subpd	%%xmm1,%%xmm4		/* ~t22 */\n\t"\
		"subpd	%%xmm3,%%xmm6		/* rt */\n\t"\
		"addpd	%%xmm0,%%xmm5		/* ~t23 */\n\t"\
		"addpd	%%xmm2,%%xmm7		/* it */\n\t"\
		"\n\t"\
		"movaps	0x40(%%esi),%%xmm2	/* t12 */\n\t"\
		"movaps	0x50(%%esi),%%xmm0	/* t13 */\n\t"\
		"movaps	0x40(%%esi),%%xmm1	/* copy t12 */\n\t"\
		"movaps	0x50(%%esi),%%xmm3	/* copy t13 */\n\t"\
		"subpd	%%xmm6,%%xmm4		/*~t32=t22-rt */\n\t"\
		"mulpd	    (%%edi),%%xmm2	/* t12*c */\n\t"\
		"subpd	%%xmm7,%%xmm5		/*~t33=t23-it */\n\t"\
		"mulpd	0x10(%%edi),%%xmm0	/* t13*s */\n\t"\
		"addpd	%%xmm6,%%xmm6		/*      2* rt */\n\t"\
		"mulpd	    (%%edi),%%xmm3	/* t13*c */\n\t"\
		"addpd	%%xmm7,%%xmm7		/*      2* it */\n\t"\
		"mulpd	0x10(%%edi),%%xmm1	/* t12*s */\n\t"\
		"addpd	%%xmm4,%%xmm6		/*~t22=t22+rt */\n\t"\
		"subpd	%%xmm0,%%xmm2		/* rt */\n\t"\
		"addpd	%%xmm5,%%xmm7		/*~t23=t23+it */\n\t"\
		"addpd	%%xmm1,%%xmm3		/* it */\n\t"\
		"\n\t"\
		"movaps	    (%%esi),%%xmm0	/* t02 */\n\t"\
		"movaps	0x10(%%esi),%%xmm1	/* t03 */\n\t"\
		"subpd	%%xmm2,%%xmm0		/*~t12=t02-rt */\n\t"\
		"subpd	%%xmm3,%%xmm1		/*~t13=t03-it */\n\t"\
		"addpd	    (%%esi),%%xmm2	/*~t02=rt+t02 */\n\t"\
		"addpd	0x10(%%esi),%%xmm3	/*~t03=it+t03 */\n\t"\
		"subpd	%%xmm6,%%xmm2		/* t02-t22 */\n\t"\
		"subpd	%%xmm5,%%xmm0		/* t12-t33 */\n\t"\
		"subpd	%%xmm7,%%xmm3		/* t03-t23 */\n\t"\
		"subpd	%%xmm4,%%xmm1		/* t13-t32 */\n\t"\
		"addpd	%%xmm6,%%xmm6		/*   2*t22 */\n\t"\
		"addpd	%%xmm5,%%xmm5		/*   2*t33 */\n\t"\
		"addpd	%%xmm7,%%xmm7		/*   2*t23 */\n\t"\
		"addpd	%%xmm4,%%xmm4		/*   2*t32 */\n\t"\
		"movaps	%%xmm2,    (%%ebx)	/* a(jt+p1 ) */\n\t"\
		"movaps	%%xmm0,    (%%ecx)	/* a(jt+p2 ) */\n\t"\
		"movaps	%%xmm3,0x10(%%ebx)	/* a(jp+p1 ) */\n\t"\
		"movaps	%%xmm1,0x10(%%edx)	/* a(jp+p3 ) */\n\t"\
		"addpd	%%xmm2,%%xmm6		/* t02+t22 */\n\t"\
		"addpd	%%xmm0,%%xmm5		/* t12+t33 */\n\t"\
		"addpd	%%xmm3,%%xmm7		/* t03+t23 */\n\t"\
		"addpd	%%xmm1,%%xmm4		/* t13+t32 */\n\t"\
		"movaps	%%xmm6,    (%%eax)	/* a(jt+p0 ) */\n\t"\
		"movaps	%%xmm5,    (%%edx)	/* a(jt+p3 ) */\n\t"\
		"movaps	%%xmm7,0x10(%%eax)	/* a(jp+p0 ) */\n\t"\
		"movaps	%%xmm4,0x10(%%ecx)	/* a(jp+p2 ) */\n\t"\
		"/*...Block 6: t0A,t1A,t2A,t3A	*/\n\t"\
		"addl	$0x80,%%esi			/* r18 */\n\t"\
		"movl	%[__p04],%%edi		\n\t"\
		"shll	$3,%%edi			\n\t"\
		"addl	%%edi,%%eax			/* ap14 */\n\t"\
		"addl	%%edi,%%ebx			\n\t"\
		"addl	%%edi,%%ecx			\n\t"\
		"addl	%%edi,%%edx			\n\t"\
		"movl	%[__cc0],%%edi	\n\t"\
		"\n\t"\
		"movaps	0x20(%%esi),%%xmm4	/* t2A */\n\t"\
		"movaps	0x60(%%esi),%%xmm6	/* t3A */\n\t"\
		"movaps	0x30(%%esi),%%xmm5	/* t2B */\n\t"\
		"movaps	0x70(%%esi),%%xmm7	/* t3B */\n\t"\
		"movaps	0x20(%%esi),%%xmm0	/* copy t2A */\n\t"\
		"movaps	0x60(%%esi),%%xmm2	/* copy t3A */\n\t"\
		"movaps	0x30(%%esi),%%xmm1	/* copy t2B */\n\t"\
		"movaps	0x70(%%esi),%%xmm3	/* copy t3B */\n\t"\
		"mulpd	0x50(%%edi),%%xmm4	/* t2A*s32_3 */\n\t"\
		"mulpd	0x20(%%edi),%%xmm6	/* t3A*c32_1 */\n\t"\
		"mulpd	0x40(%%edi),%%xmm1	/* t2B*c32_3 */\n\t"\
		"mulpd	0x30(%%edi),%%xmm3	/* t3B*s32_1 */\n\t"\
		"mulpd	0x50(%%edi),%%xmm5	/* t2B*s32_3 */\n\t"\
		"mulpd	0x20(%%edi),%%xmm7	/* t3B*c32_1 */\n\t"\
		"mulpd	0x40(%%edi),%%xmm0	/* t2A*c32_3 */\n\t"\
		"mulpd	0x30(%%edi),%%xmm2	/* t3A*s32_1 */\n\t"\
		"subpd	%%xmm1,%%xmm4		/* ~t2A */\n\t"\
		"addpd	%%xmm3,%%xmm6		/* rt */\n\t"\
		"addpd	%%xmm0,%%xmm5		/* ~t2B */\n\t"\
		"subpd	%%xmm2,%%xmm7		/* it */\n\t"\
		"\n\t"\
		"movaps	0x40(%%esi),%%xmm2	/* t1A */\n\t"\
		"movaps	0x50(%%esi),%%xmm0	/* t1B */\n\t"\
		"movaps	0x40(%%esi),%%xmm1	/* copy t1A */\n\t"\
		"movaps	0x50(%%esi),%%xmm3	/* copy t1B */\n\t"\
		"subpd	%%xmm6,%%xmm4		/*~t2A=t2A-rt */\n\t"\
		"mulpd	0x10(%%edi),%%xmm2	/* t1A*s */\n\t"\
		"subpd	%%xmm7,%%xmm5		/*~t2B=t2B-it */\n\t"\
		"mulpd	    (%%edi),%%xmm0	/* t1B*c */\n\t"\
		"addpd	%%xmm6,%%xmm6		/*      2* rt */\n\t"\
		"mulpd	0x10(%%edi),%%xmm3	/* t1B*s */\n\t"\
		"addpd	%%xmm7,%%xmm7		/*      2* it */\n\t"\
		"mulpd	    (%%edi),%%xmm1	/* t1A*c */\n\t"\
		"addpd	%%xmm4,%%xmm6		/*~t3A=t2A+rt */\n\t"\
		"addpd	%%xmm0,%%xmm2		/* rt */\n\t"\
		"addpd	%%xmm5,%%xmm7		/*~t3B=t2B+it */\n\t"\
		"subpd	%%xmm1,%%xmm3		/* it */\n\t"\
		"\n\t"\
		"movaps	    (%%esi),%%xmm0	/* t0A */\n\t"\
		"movaps	0x10(%%esi),%%xmm1	/* t0B */\n\t"\
		"subpd	%%xmm2,%%xmm0		/*~t0A=t0A-rt */\n\t"\
		"subpd	%%xmm3,%%xmm1		/*~t0B=t0B-it */\n\t"\
		"addpd	    (%%esi),%%xmm2	/*~t1A=rt+t0A */\n\t"\
		"addpd	0x10(%%esi),%%xmm3	/*~t1B=it+t0B */\n\t"\
		"subpd	%%xmm4,%%xmm0		/* t0A-t2A */\n\t"\
		"subpd	%%xmm7,%%xmm2		/* t1A-t3B */\n\t"\
		"subpd	%%xmm5,%%xmm1		/* t0B-t2B */\n\t"\
		"subpd	%%xmm6,%%xmm3		/* t1B-t3A */\n\t"\
		"addpd	%%xmm4,%%xmm4		/*   2*t2A */\n\t"\
		"addpd	%%xmm7,%%xmm7		/*   2*t3B */\n\t"\
		"addpd	%%xmm5,%%xmm5		/*   2*t2B */\n\t"\
		"addpd	%%xmm6,%%xmm6		/*   2*t3A */\n\t"\
		"movaps	%%xmm0,    (%%ebx)	/* a(jt+p1 ) */\n\t"\
		"movaps	%%xmm2,    (%%ecx)	/* a(jt+p2 ) */\n\t"\
		"movaps	%%xmm1,0x10(%%ebx)	/* a(jp+p1 ) */\n\t"\
		"movaps	%%xmm3,0x10(%%edx)	/* a(jp+p3 ) */\n\t"\
		"addpd	%%xmm0,%%xmm4		/* t0A+t2A */\n\t"\
		"addpd	%%xmm2,%%xmm7		/* t1A+t3B */\n\t"\
		"addpd	%%xmm1,%%xmm5		/* t0B+t2B */\n\t"\
		"addpd	%%xmm3,%%xmm6		/* t1B+t3A */\n\t"\
		"movaps	%%xmm4,    (%%eax)	/* a(jt+p0 ) */\n\t"\
		"movaps	%%xmm7,    (%%edx)	/* a(jt+p3 ) */\n\t"\
		"movaps	%%xmm5,0x10(%%eax)	/* a(jp+p0 ) */\n\t"\
		"movaps	%%xmm6,0x10(%%ecx)	/* a(jp+p2 ) */\n\t"\
		"/*...Block 4: t06,t16,t26,t36	*/\n\t"\
		"addl	$0x180,%%esi		/* r30 */\n\t"\
		"movl	%[__p04],%%edi		\n\t"\
		"shll	$3,%%edi			\n\t"\
		"addl	%%edi,%%eax			/* ap18 */\n\t"\
		"addl	%%edi,%%ebx			\n\t"\
		"addl	%%edi,%%ecx			\n\t"\
		"addl	%%edi,%%edx			\n\t"\
		"movl	%[__cc0],%%edi	\n\t"\
		"\n\t"\
		"movaps	0x20(%%esi),%%xmm4	/* t26 */\n\t"\
		"movaps	0x60(%%esi),%%xmm6	/* t36 */\n\t"\
		"movaps	0x30(%%esi),%%xmm5	/* t27 */\n\t"\
		"movaps	0x70(%%esi),%%xmm7	/* t37 */\n\t"\
		"movaps	0x20(%%esi),%%xmm0	/* copy t26 */\n\t"\
		"movaps	0x60(%%esi),%%xmm2	/* copy t36 */\n\t"\
		"movaps	0x30(%%esi),%%xmm1	/* copy t27 */\n\t"\
		"movaps	0x70(%%esi),%%xmm3	/* copy t37 */\n\t"\
		"mulpd	0x40(%%edi),%%xmm4	/* t26*s32_3 */\n\t"\
		"mulpd	0x30(%%edi),%%xmm6	/* t36*s32_1 */\n\t"\
		"mulpd	0x50(%%edi),%%xmm1	/* t27*s32_3 */\n\t"\
		"mulpd	0x20(%%edi),%%xmm3	/* t37*c32_1 */\n\t"\
		"mulpd	0x40(%%edi),%%xmm5	/* t27*c32_3 */\n\t"\
		"mulpd	0x30(%%edi),%%xmm7	/* t37*s32_1 */\n\t"\
		"mulpd	0x50(%%edi),%%xmm0	/* t26*s32_3 */\n\t"\
		"mulpd	0x20(%%edi),%%xmm2	/* t36*c32_1 */\n\t"\
		"subpd	%%xmm1,%%xmm4		/* ~t26 */\n\t"\
		"addpd	%%xmm3,%%xmm6		/* rt */\n\t"\
		"addpd	%%xmm0,%%xmm5		/* ~t27 */\n\t"\
		"subpd	%%xmm2,%%xmm7		/* it */\n\t"\
		"\n\t"\
		"movaps	0x40(%%esi),%%xmm2	/* t16 */\n\t"\
		"movaps	0x50(%%esi),%%xmm0	/* t17 */\n\t"\
		"movaps	0x40(%%esi),%%xmm1	/* copy t16 */\n\t"\
		"movaps	0x50(%%esi),%%xmm3	/* copy t17 */\n\t"\
		"subpd	%%xmm6,%%xmm4		/*~t26=t26-rt */\n\t"\
		"mulpd	0x10(%%edi),%%xmm2	/* t16*s */\n\t"\
		"subpd	%%xmm7,%%xmm5		/*~t27=t27-it */\n\t"\
		"mulpd	    (%%edi),%%xmm0	/* t17*c */\n\t"\
		"addpd	%%xmm6,%%xmm6		/*      2* rt */\n\t"\
		"mulpd	0x10(%%edi),%%xmm3	/* t17*s */\n\t"\
		"addpd	%%xmm7,%%xmm7		/*      2* it */\n\t"\
		"mulpd	    (%%edi),%%xmm1	/* t16*c */\n\t"\
		"addpd	%%xmm4,%%xmm6		/*~t36=t26+rt */\n\t"\
		"subpd	%%xmm0,%%xmm2		/* rt */\n\t"\
		"addpd	%%xmm5,%%xmm7		/*~t37=t27+it */\n\t"\
		"addpd	%%xmm1,%%xmm3		/* it */\n\t"\
		"\n\t"\
		"movaps	    (%%esi),%%xmm0	/* t06 */\n\t"\
		"movaps	0x10(%%esi),%%xmm1	/* t07 */\n\t"\
		"subpd	%%xmm2,%%xmm0		/*~t16=t06-rt */\n\t"\
		"subpd	%%xmm3,%%xmm1		/*~t17=t07-it */\n\t"\
		"addpd	    (%%esi),%%xmm2	/*~t06=rt+t06 */\n\t"\
		"addpd	0x10(%%esi),%%xmm3	/*~t07=it+t07 */\n\t"\
		"subpd	%%xmm4,%%xmm2		/* t06-t26 */\n\t"\
		"subpd	%%xmm7,%%xmm0		/* t16-t37 */\n\t"\
		"subpd	%%xmm5,%%xmm3		/* t07-t27 */\n\t"\
		"subpd	%%xmm6,%%xmm1		/* t17-t36 */\n\t"\
		"addpd	%%xmm4,%%xmm4		/*   2*t26 */\n\t"\
		"addpd	%%xmm7,%%xmm7		/*   2*t37 */\n\t"\
		"addpd	%%xmm5,%%xmm5		/*   2*t27 */\n\t"\
		"addpd	%%xmm6,%%xmm6		/*   2*t36 */\n\t"\
		"movaps	%%xmm2,    (%%ebx)	/* a(jt+p1 ) */\n\t"\
		"movaps	%%xmm0,    (%%ecx)	/* a(jt+p2 ) */\n\t"\
		"movaps	%%xmm3,0x10(%%ebx)	/* a(jp+p1 ) */\n\t"\
		"movaps	%%xmm1,0x10(%%edx)	/* a(jp+p3 ) */\n\t"\
		"addpd	%%xmm2,%%xmm4		/* t06+t26 */\n\t"\
		"addpd	%%xmm0,%%xmm7		/* t16+t37 */\n\t"\
		"addpd	%%xmm3,%%xmm5		/* t07+t27 */\n\t"\
		"addpd	%%xmm1,%%xmm6		/* t17+t36 */\n\t"\
		"movaps	%%xmm4,    (%%eax)	/* a(jt+p0 ) */\n\t"\
		"movaps	%%xmm7,    (%%edx)	/* a(jt+p3 ) */\n\t"\
		"movaps	%%xmm5,0x10(%%eax)	/* a(jp+p0 ) */\n\t"\
		"movaps	%%xmm6,0x10(%%ecx)	/* a(jp+p2 ) */\n\t"\
		"/*...Block 8: t0E,t1E,t2E,t3E	*/\n\t"\
		"addl	$0x80,%%esi			/* r38 */\n\t"\
		"movl	%[__p04],%%edi		\n\t"\
		"shll	$3,%%edi			\n\t"\
		"addl	%%edi,%%eax			/* ap1C */\n\t"\
		"addl	%%edi,%%ebx			\n\t"\
		"addl	%%edi,%%ecx			\n\t"\
		"addl	%%edi,%%edx			\n\t"\
		"movl	%[__cc0],%%edi	\n\t"\
		"\n\t"\
		"movaps	0x20(%%esi),%%xmm4	/* t2E */\n\t"\
		"movaps	0x60(%%esi),%%xmm6	/* t3E */\n\t"\
		"movaps	0x30(%%esi),%%xmm5	/* t2F */\n\t"\
		"movaps	0x70(%%esi),%%xmm7	/* t3F */\n\t"\
		"movaps	0x20(%%esi),%%xmm0	/* copy t2E */\n\t"\
		"movaps	0x60(%%esi),%%xmm2	/* copy t3E */\n\t"\
		"movaps	0x30(%%esi),%%xmm1	/* copy t2F */\n\t"\
		"movaps	0x70(%%esi),%%xmm3	/* copy t3F */\n\t"\
		"mulpd	0x30(%%edi),%%xmm4	/* t2E*s32_1 */\n\t"\
		"mulpd	0x50(%%edi),%%xmm6	/* t3E*c32_3 */\n\t"\
		"mulpd	0x20(%%edi),%%xmm1	/* t2F*c32_1 */\n\t"\
		"mulpd	0x40(%%edi),%%xmm3	/* t3F*s32_3 */\n\t"\
		"mulpd	0x30(%%edi),%%xmm5	/* t2F*s32_1 */\n\t"\
		"mulpd	0x50(%%edi),%%xmm7	/* t3F*c32_3 */\n\t"\
		"mulpd	0x20(%%edi),%%xmm0	/* t2E*c32_1 */\n\t"\
		"mulpd	0x40(%%edi),%%xmm2	/* t3E*s32_3 */\n\t"\
		"subpd	%%xmm1,%%xmm4		/* ~t2E */\n\t"\
		"subpd	%%xmm3,%%xmm6		/* rt */\n\t"\
		"addpd	%%xmm0,%%xmm5		/* ~t2F */\n\t"\
		"addpd	%%xmm2,%%xmm7		/* it */\n\t"\
		"\n\t"\
		"movaps	0x40(%%esi),%%xmm2	/* t1E */\n\t"\
		"movaps	0x50(%%esi),%%xmm0	/* t1F */\n\t"\
		"movaps	0x40(%%esi),%%xmm1	/* copy t1E */\n\t"\
		"movaps	0x50(%%esi),%%xmm3	/* copy t1F */\n\t"\
		"subpd	%%xmm6,%%xmm4		/*~t2E=t2E-rt */\n\t"\
		"mulpd	    (%%edi),%%xmm2	/* t1E*c */\n\t"\
		"subpd	%%xmm7,%%xmm5		/*~t2F=t2F-it */\n\t"\
		"mulpd	0x10(%%edi),%%xmm0	/* t1F*s */\n\t"\
		"addpd	%%xmm6,%%xmm6		/*      2* rt */\n\t"\
		"mulpd	    (%%edi),%%xmm3	/* t1F*c */\n\t"\
		"addpd	%%xmm7,%%xmm7		/*      2* it */\n\t"\
		"mulpd	0x10(%%edi),%%xmm1	/* t1E*s */\n\t"\
		"addpd	%%xmm4,%%xmm6		/*~t3E=t2E+rt */\n\t"\
		"addpd	%%xmm0,%%xmm2		/* rt */\n\t"\
		"addpd	%%xmm5,%%xmm7		/*~t3F=t2F+it */\n\t"\
		"subpd	%%xmm1,%%xmm3		/* it */\n\t"\
		"\n\t"\
		"movaps	    (%%esi),%%xmm0	/* t0E */\n\t"\
		"movaps	0x10(%%esi),%%xmm1	/* t0F */\n\t"\
		"subpd	%%xmm2,%%xmm0		/*~t0E=t0E-rt */\n\t"\
		"subpd	%%xmm3,%%xmm1		/*~t0F=t0F-it */\n\t"\
		"addpd	    (%%esi),%%xmm2	/*~t1E=rt+t0E */\n\t"\
		"addpd	0x10(%%esi),%%xmm3	/*~t1F=it+t0F */\n\t"\
		"subpd	%%xmm4,%%xmm0		/* t0E-t2E */\n\t"\
		"subpd	%%xmm7,%%xmm2		/* t1E-t3F */\n\t"\
		"subpd	%%xmm5,%%xmm1		/* t0F-t2F */\n\t"\
		"subpd	%%xmm6,%%xmm3		/* t1F-t3E */\n\t"\
		"addpd	%%xmm4,%%xmm4		/*   2*t2E */\n\t"\
		"addpd	%%xmm7,%%xmm7		/*   2*t3F */\n\t"\
		"addpd	%%xmm5,%%xmm5		/*   2*t2F */\n\t"\
		"addpd	%%xmm6,%%xmm6		/*   2*t3E */\n\t"\
		"movaps	%%xmm0,    (%%ebx)	/* a(jt+p1 ) */\n\t"\
		"movaps	%%xmm2,    (%%ecx)	/* a(jt+p2 ) */\n\t"\
		"movaps	%%xmm1,0x10(%%ebx)	/* a(jp+p1 ) */\n\t"\
		"movaps	%%xmm3,0x10(%%edx)	/* a(jp+p3 ) */\n\t"\
		"addpd	%%xmm0,%%xmm4		/* t0E+t2E */\n\t"\
		"addpd	%%xmm2,%%xmm7		/* t1E+t3F */\n\t"\
		"addpd	%%xmm1,%%xmm5		/* t0F+t2F */\n\t"\
		"addpd	%%xmm3,%%xmm6		/* t1F+t3E */\n\t"\
		"movaps	%%xmm4,    (%%eax)	/* a(jt+p0 ) */\n\t"\
		"movaps	%%xmm7,    (%%edx)	/* a(jt+p3 ) */\n\t"\
		"movaps	%%xmm5,0x10(%%eax)	/* a(jp+p0 ) */\n\t"\
		"movaps	%%xmm6,0x10(%%ecx)	/* a(jp+p2 ) */\n\t"\
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
		: "eax","ebx","ecx","edx","edi","esi"		/* Clobbered registers */\
	);\
	}

#endif	/* radix32_ditN_cy_dif1_gcc_h_included */

