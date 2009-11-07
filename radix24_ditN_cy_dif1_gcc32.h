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
			"movl	%[__p04],%%edi	\n\t"\
			"movl	%[__add],%%eax	/* Use eax as base address throughout */\n\t"\
			"shll	$3,%%edi		\n\t"\
			"movl	%[__p01],%%ebx	\n\t"\
			"movl	%[__p02],%%ecx	\n\t"\
			"movl	%[__p03],%%edx	\n\t"\
			"addl	%%edi,%%eax		/* eax <- add0+p04 */\n\t"\
			"shll	$3,%%ebx		/* Pointer offset for floating doubles */\n\t"\
			"shll	$3,%%ecx		\n\t"\
			"shll	$3,%%edx		\n\t"\
			"addl	%%eax,%%ebx		/* eax <- add0+p05 */\n\t"\
			"addl	%%eax,%%ecx		/* ebx <- add0+p06 */\n\t"\
			"addl	%%eax,%%edx		/* ecx <- add0+p07 */\n\t"\
			"movl	%[__out],%%esi	/* s1p00r */\n\t"\
			"/* MSVC macro assumes add8+p[7,6,5,4] in eax,ebx,ecx,edx, but here get add0+p[4,5,6,7], so replace eax <-> edx and ebx <-> ecx: */\n\t"\
			"/* Do the p0,p4 combo: */\n\t"\
			"movaps	    (%%edx),%%xmm0			\n\t"\
			"movaps	0x10(%%edx),%%xmm1			\n\t"\
			"movaps	%%xmm0,%%xmm2				\n\t"\
			"movaps	%%xmm1,%%xmm3				\n\t"\
			"addpd	    (%%ecx),%%xmm2			\n\t"\
			"addpd	0x10(%%ecx),%%xmm3			\n\t"\
			"subpd	    (%%ecx),%%xmm0			\n\t"\
			"subpd	0x10(%%ecx),%%xmm1			\n\t"\
			"\n\t"\
			"movaps	    (%%ebx),%%xmm4			\n\t"\
			"movaps	0x10(%%ebx),%%xmm5			\n\t"\
			"movaps	%%xmm4,%%xmm6				\n\t"\
			"movaps	%%xmm5,%%xmm7				\n\t"\
			"addpd	    (%%eax),%%xmm6			\n\t"\
			"addpd	0x10(%%eax),%%xmm7			\n\t"\
			"subpd	    (%%eax),%%xmm4			\n\t"\
			"subpd	0x10(%%eax),%%xmm5			\n\t"\
			"/* Copy t6r,i into main-array slot add6 */\n\t"\
			"movaps	%%xmm6,0xc0(%%esi)			\n\t"\
			"movaps	%%xmm7,0xd0(%%esi)			\n\t"\
			"/* Copy t7r,i into main-array slot add7 */\n\t"\
			"movaps	%%xmm4,0xe0(%%esi)			\n\t"\
			"movaps	%%xmm5,0xf0(%%esi)			\n\t"\
			"\n\t"\
			"/* Move outputs t5r,i into a(j1,j2+p5), first doing the addsub and mul by ISRT2: */\n\t"\
			"/* Move outputs t7r,i into a(j1,j2+p7), first doing the addsub and mul by ISRT2: */\n\t"\
			"\n\t"\
			"addpd	%%xmm2,%%xmm6				\n\t"\
			"addpd	%%xmm3,%%xmm7				\n\t"\
			"subpd	0xc0(%%esi),%%xmm2			\n\t"\
			"subpd	0xd0(%%esi),%%xmm3			\n\t"\
			"/* Move t6r,i into a(j1,j2+p6) */	\n\t"\
			"movaps	%%xmm2,0xc0(%%esi)			\n\t"\
			"movaps	%%xmm3,0xd0(%%esi)			\n\t"\
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
			"movaps	0x300(%%esi),%%xmm1	/* isrt2 */	\n\t"\
			"subpd	%%xmm3,%%xmm2				\n\t"\
			"mulpd	%%xmm1,%%xmm5				\n\t"\
			"mulpd	%%xmm1,%%xmm2				\n\t"\
			"movaps	%%xmm0,%%xmm3				\n\t"\
			"\n\t"\
			"movaps	%%xmm5,0xa0(%%esi)			\n\t"\
			"movaps	%%xmm4,%%xmm5				\n\t"\
			"addpd	%%xmm4,%%xmm0				\n\t"\
			"movaps	%%xmm2,0xb0(%%esi)			\n\t"\
			"subpd	%%xmm5,%%xmm3				\n\t"\
			"mulpd	%%xmm1,%%xmm0				\n\t"\
			"mulpd	%%xmm1,%%xmm3				\n\t"\
			"movaps	%%xmm0,0xe0(%%esi)			\n\t"\
			"movaps	%%xmm3,0xf0(%%esi)			\n\t"\
			"\n\t"\
			"/* MSVC macro assumes add0+p[0,1,3,2] in eax,ebx,ecx,edx, but here get add0+p[0,1,2,3], so replace ecx <-> edx: */\n\t"\
			"subl	%%edi,%%eax		/* add2 = add0     */\n\t"\
			"subl	%%edi,%%ebx		/* add3 = add0+p01 */\n\t"\
			"subl	%%edi,%%ecx		/* add1 = add0+p02 */\n\t"\
			"subl	%%edi,%%edx		/* add0 = add0+p03 */\n\t"\
			"\n\t"\
			"movaps	    (%%eax),%%xmm0			\n\t"\
			"movaps	0x10(%%eax),%%xmm1			\n\t"\
			"movaps	%%xmm0,%%xmm2				\n\t"\
			"movaps	%%xmm1,%%xmm3				\n\t"\
			"addpd	    (%%ebx),%%xmm2			\n\t"\
			"addpd	0x10(%%ebx),%%xmm3			\n\t"\
			"subpd	    (%%ebx),%%xmm0			\n\t"\
			"subpd	0x10(%%ebx),%%xmm1			\n\t"\
			"\n\t"\
			"/* Move   t4r,i into temp(j1+p0) in anticipation of final outputs (t0+t4)r,i which will go there: */\n\t"\
			"movaps	%%xmm6,    (%%esi)			\n\t"\
			"movaps	%%xmm7,0x10(%%esi)			\n\t"\
			"addpd	%%xmm6,%%xmm6				\n\t"\
			"addpd	%%xmm7,%%xmm7				\n\t"\
			"/* Move 2*t4r,i into temp(j1+p4) in anticipation of final outputs (t0-t4)r,i which will go there: */\n\t"\
			"movaps	%%xmm6,0x80(%%esi)			\n\t"\
			"movaps	%%xmm7,0x90(%%esi)			\n\t"\
			"\n\t"\
			"movaps	    (%%edx),%%xmm4			\n\t"\
			"movaps	0x10(%%edx),%%xmm5			\n\t"\
			"movaps	%%xmm4,%%xmm6				\n\t"\
			"movaps	%%xmm5,%%xmm7				\n\t"\
			"addpd	    (%%ecx),%%xmm6			\n\t"\
			"addpd	0x10(%%ecx),%%xmm7			\n\t"\
			"subpd	    (%%ecx),%%xmm4			\n\t"\
			"subpd	0x10(%%ecx),%%xmm5			\n\t"\
			"/* Copy t5,6 into temp(j1+p2) */	\n\t"\
			"movaps	%%xmm6,0x40(%%esi)			\n\t"\
			"movaps	%%xmm7,0x50(%%esi)			\n\t"\
			"addpd	%%xmm2,%%xmm6				\n\t"\
			"addpd	%%xmm3,%%xmm7				\n\t"\
			"subpd	0x40(%%esi),%%xmm2			\n\t"\
			"subpd	0x50(%%esi),%%xmm3			\n\t"\
			"\n\t"\
			"/* Compute and dump first 2 outputs now, in order to free up 2 registers: */\n\t"\
			"addpd	    (%%esi),%%xmm6			\n\t"\
			"addpd	0x10(%%esi),%%xmm7			\n\t"\
			"movaps	%%xmm6,    (%%esi)			\n\t"\
			"movaps	%%xmm7,0x10(%%esi)			\n\t"\
			"\n\t"\
			"subpd	0x80(%%esi),%%xmm6			\n\t"\
			"subpd	0x90(%%esi),%%xmm7			\n\t"\
			"movaps	%%xmm6,0x80(%%esi)			\n\t"\
			"movaps	%%xmm7,0x90(%%esi)			\n\t"\
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
			"addpd	0xd0(%%esi),%%xmm2			\n\t"\
			"subpd	0xc0(%%esi),%%xmm3			\n\t"\
			"subpd	0xd0(%%esi),%%xmm6			\n\t"\
			"addpd	0xc0(%%esi),%%xmm7			\n\t"\
			"movaps	%%xmm2,0xc0(%%esi)			\n\t"\
			"movaps	%%xmm3,0xd0(%%esi)			\n\t"\
			"movaps	%%xmm6,0x40(%%esi)			\n\t"\
			"movaps	%%xmm7,0x50(%%esi)			\n\t"\
			"\n\t"\
			"movaps	%%xmm5,%%xmm2				\n\t"\
			"movaps	%%xmm0,%%xmm6				\n\t"\
			"movaps	%%xmm1,%%xmm3				\n\t"\
			"movaps	%%xmm4,%%xmm7				\n\t"\
			"addpd	0xa0(%%esi),%%xmm5			\n\t"\
			"subpd	0xf0(%%esi),%%xmm0			\n\t"\
			"subpd	0xb0(%%esi),%%xmm1			\n\t"\
			"subpd	0xe0(%%esi),%%xmm4			\n\t"\
			"subpd	0xa0(%%esi),%%xmm2			\n\t"\
			"addpd	0xf0(%%esi),%%xmm6			\n\t"\
			"addpd	0xb0(%%esi),%%xmm3			\n\t"\
			"addpd	0xe0(%%esi),%%xmm7			\n\t"\
			"movaps	%%xmm5,0xe0(%%esi)			\n\t"\
			"movaps	%%xmm0,0xa0(%%esi)			\n\t"\
			"movaps	%%xmm1,0xf0(%%esi)			\n\t"\
			"movaps	%%xmm4,0xb0(%%esi)			\n\t"\
			"movaps	%%xmm2,0x60(%%esi)			\n\t"\
			"movaps	%%xmm6,0x20(%%esi)			\n\t"\
			"movaps	%%xmm3,0x70(%%esi)			\n\t"\
			"movaps	%%xmm7,0x30(%%esi)			\n\t"\
			"\n\t"\
		"/* SSE2_RADIX8_DIT_0TWIDDLE([add0+p08]+p[5,4,6,7,1,0,2,3], s1p08r) */\n\t"\
			"movl	%[__p08],%%edi	\n\t"\
			"shll	$3,%%edi		/* Pointer offset for floating doubles */\n\t"\
			"addl	%%edi,%%eax		/* edx <- add8     */\n\t"\
			"addl	%%edi,%%ebx		/* ecx <- add8+p01 */\n\t"\
			"addl	%%edi,%%ecx		/* ebx <- add8+p02 */\n\t"\
			"addl	%%edi,%%edx		/* eax <- add8+p03 */\n\t"\
			"addl	$0x100,%%esi	/* s1p08r */\n\t"\
			"/* MSVC macro assumes add8+p[0,1,2,3] in eax,ebx,ecx,edx, but here get add+p[1,0,2,3], so replace eax <-> ebx: */\n\t"\
			"/* Do the p0,p4 combo: */\n\t"\
			"movaps	    (%%ebx),%%xmm0			\n\t"\
			"movaps	0x10(%%ebx),%%xmm1			\n\t"\
			"movaps	%%xmm0,%%xmm2				\n\t"\
			"movaps	%%xmm1,%%xmm3				\n\t"\
			"addpd	    (%%eax),%%xmm2			\n\t"\
			"addpd	0x10(%%eax),%%xmm3			\n\t"\
			"subpd	    (%%eax),%%xmm0			\n\t"\
			"subpd	0x10(%%eax),%%xmm1			\n\t"\
			"\n\t"\
			"movaps	    (%%ecx),%%xmm4			\n\t"\
			"movaps	0x10(%%ecx),%%xmm5			\n\t"\
			"movaps	%%xmm4,%%xmm6				\n\t"\
			"movaps	%%xmm5,%%xmm7				\n\t"\
			"addpd	    (%%edx),%%xmm6			\n\t"\
			"addpd	0x10(%%edx),%%xmm7			\n\t"\
			"subpd	    (%%edx),%%xmm4			\n\t"\
			"subpd	0x10(%%edx),%%xmm5			\n\t"\
			"/* Copy t6r,i into main-array slot add6 */\n\t"\
			"movaps	%%xmm6,0xc0(%%esi)			\n\t"\
			"movaps	%%xmm7,0xd0(%%esi)			\n\t"\
			"/* Copy t7r,i into main-array slot add7 */\n\t"\
			"movaps	%%xmm4,0xe0(%%esi)			\n\t"\
			"movaps	%%xmm5,0xf0(%%esi)			\n\t"\
			"\n\t"\
			"/* Move outputs t5r,i into a(j1,j2+p5), first doing the addsub and mul by ISRT2: */\n\t"\
			"/* Move outputs t7r,i into a(j1,j2+p7), first doing the addsub and mul by ISRT2: */\n\t"\
			"\n\t"\
			"addpd	%%xmm2,%%xmm6				\n\t"\
			"addpd	%%xmm3,%%xmm7				\n\t"\
			"subpd	0xc0(%%esi),%%xmm2			\n\t"\
			"subpd	0xd0(%%esi),%%xmm3			\n\t"\
			"/* Move t6r,i into a(j1,j2+p6) */	\n\t"\
			"movaps	%%xmm2,0xc0(%%esi)			\n\t"\
			"movaps	%%xmm3,0xd0(%%esi)			\n\t"\
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
			"movaps	0x200(%%esi),%%xmm1	/* isrt2 */	\n\t"\
			"subpd	%%xmm3,%%xmm2				\n\t"\
			"mulpd	%%xmm1,%%xmm5				\n\t"\
			"mulpd	%%xmm1,%%xmm2				\n\t"\
			"movaps	%%xmm0,%%xmm3				\n\t"\
			"\n\t"\
			"movaps	%%xmm5,0xa0(%%esi)			\n\t"\
			"movaps	%%xmm4,%%xmm5				\n\t"\
			"addpd	%%xmm4,%%xmm0				\n\t"\
			"movaps	%%xmm2,0xb0(%%esi)			\n\t"\
			"subpd	%%xmm5,%%xmm3				\n\t"\
			"mulpd	%%xmm1,%%xmm0				\n\t"\
			"mulpd	%%xmm1,%%xmm3				\n\t"\
			"movaps	%%xmm0,0xe0(%%esi)			\n\t"\
			"movaps	%%xmm3,0xf0(%%esi)			\n\t"\
			"\n\t"\
			"/* MSVC macro assumes add8+p[5,4,6,7] in eax,ebx,ecx,edx, but here get add8+p[4,5,6,7], so replace eax <-> ebx: */\n\t"\
			"movl	%[__p04],%%edi	\n\t"\
			"shll	$3,%%edi		/* Pointer offset for floating doubles */\n\t"\
			"addl	%%edi,%%eax		/* add0 = add8+p04 */\n\t"\
			"addl	%%edi,%%ebx		/* add1 = add8+p05 */\n\t"\
			"addl	%%edi,%%ecx		/* add3 = add8+p06 */\n\t"\
			"addl	%%edi,%%edx		/* add2 = add8+p07 */\n\t"\
			"\n\t"\
			"movaps	    (%%ebx),%%xmm0			\n\t"\
			"movaps	0x10(%%ebx),%%xmm1			\n\t"\
			"movaps	%%xmm0,%%xmm2				\n\t"\
			"movaps	%%xmm1,%%xmm3				\n\t"\
			"addpd	    (%%eax),%%xmm2			\n\t"\
			"addpd	0x10(%%eax),%%xmm3			\n\t"\
			"subpd	    (%%eax),%%xmm0			\n\t"\
			"subpd	0x10(%%eax),%%xmm1			\n\t"\
			"\n\t"\
			"/* Move   t4r,i into temp(j1+p0) in anticipation of final outputs (t0+t4)r,i which will go there: */\n\t"\
			"movaps	%%xmm6,    (%%esi)			\n\t"\
			"movaps	%%xmm7,0x10(%%esi)			\n\t"\
			"addpd	%%xmm6,%%xmm6				\n\t"\
			"addpd	%%xmm7,%%xmm7				\n\t"\
			"/* Move 2*t4r,i into temp(j1+p4) in anticipation of final outputs (t0-t4)r,i which will go there: */\n\t"\
			"movaps	%%xmm6,0x80(%%esi)			\n\t"\
			"movaps	%%xmm7,0x90(%%esi)			\n\t"\
			"\n\t"\
			"movaps	    (%%ecx),%%xmm4			\n\t"\
			"movaps	0x10(%%ecx),%%xmm5			\n\t"\
			"movaps	%%xmm4,%%xmm6				\n\t"\
			"movaps	%%xmm5,%%xmm7				\n\t"\
			"addpd	    (%%edx),%%xmm6			\n\t"\
			"addpd	0x10(%%edx),%%xmm7			\n\t"\
			"subpd	    (%%edx),%%xmm4			\n\t"\
			"subpd	0x10(%%edx),%%xmm5			\n\t"\
			"/* Copy t5,6 into temp(j1+p2) */	\n\t"\
			"movaps	%%xmm6,0x40(%%esi)			\n\t"\
			"movaps	%%xmm7,0x50(%%esi)			\n\t"\
			"addpd	%%xmm2,%%xmm6				\n\t"\
			"addpd	%%xmm3,%%xmm7				\n\t"\
			"subpd	0x40(%%esi),%%xmm2			\n\t"\
			"subpd	0x50(%%esi),%%xmm3			\n\t"\
			"\n\t"\
			"/* Compute and dump first 2 outputs now, in order to free up 2 registers: */\n\t"\
			"addpd	    (%%esi),%%xmm6			\n\t"\
			"addpd	0x10(%%esi),%%xmm7			\n\t"\
			"movaps	%%xmm6,    (%%esi)			\n\t"\
			"movaps	%%xmm7,0x10(%%esi)			\n\t"\
			"\n\t"\
			"subpd	0x80(%%esi),%%xmm6			\n\t"\
			"subpd	0x90(%%esi),%%xmm7			\n\t"\
			"movaps	%%xmm6,0x80(%%esi)			\n\t"\
			"movaps	%%xmm7,0x90(%%esi)			\n\t"\
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
			"addpd	0xd0(%%esi),%%xmm2			\n\t"\
			"subpd	0xc0(%%esi),%%xmm3			\n\t"\
			"subpd	0xd0(%%esi),%%xmm6			\n\t"\
			"addpd	0xc0(%%esi),%%xmm7			\n\t"\
			"movaps	%%xmm2,0xc0(%%esi)			\n\t"\
			"movaps	%%xmm3,0xd0(%%esi)			\n\t"\
			"movaps	%%xmm6,0x40(%%esi)			\n\t"\
			"movaps	%%xmm7,0x50(%%esi)			\n\t"\
			"\n\t"\
			"movaps	%%xmm5,%%xmm2				\n\t"\
			"movaps	%%xmm0,%%xmm6				\n\t"\
			"movaps	%%xmm1,%%xmm3				\n\t"\
			"movaps	%%xmm4,%%xmm7				\n\t"\
			"addpd	0xa0(%%esi),%%xmm5			\n\t"\
			"subpd	0xf0(%%esi),%%xmm0			\n\t"\
			"subpd	0xb0(%%esi),%%xmm1			\n\t"\
			"subpd	0xe0(%%esi),%%xmm4			\n\t"\
			"subpd	0xa0(%%esi),%%xmm2			\n\t"\
			"addpd	0xf0(%%esi),%%xmm6			\n\t"\
			"addpd	0xb0(%%esi),%%xmm3			\n\t"\
			"addpd	0xe0(%%esi),%%xmm7			\n\t"\
			"movaps	%%xmm5,0xe0(%%esi)			\n\t"\
			"movaps	%%xmm0,0xa0(%%esi)			\n\t"\
			"movaps	%%xmm1,0xf0(%%esi)			\n\t"\
			"movaps	%%xmm4,0xb0(%%esi)			\n\t"\
			"movaps	%%xmm2,0x60(%%esi)			\n\t"\
			"movaps	%%xmm6,0x20(%%esi)			\n\t"\
			"movaps	%%xmm3,0x70(%%esi)			\n\t"\
			"movaps	%%xmm7,0x30(%%esi)			\n\t"\
			"\n\t"\
		"/* SSE2_RADIX8_DIT_0TWIDDLE([add0+p16]+p[2,3,0,1,4,5,7,6], s1p16r) */\n\t"\
			"movl	%[__p08],%%edi	\n\t"\
			"shll	$3,%%edi		/* Pointer offset for floating doubles */\n\t"\
			"subl	%%edi,%%eax		/* edx <- add0 +p04 */\n\t"\
			"subl	%%edi,%%ebx		/* ecx <- add0 +p05 */\n\t"\
			"subl	%%edi,%%ecx		/* ebx <- add0 +p06 */\n\t"\
			"subl	%%edi,%%edx		/* eax <- add0 +p07 */\n\t"\
			"movl	%[__p16],%%edi	\n\t"\
			"shll	$3,%%edi		/* Pointer offset for floating doubles */\n\t"\
			"addl	%%edi,%%eax		/* edx <- add16+p04 */\n\t"\
			"addl	%%edi,%%ebx		/* ecx <- add16+p05 */\n\t"\
			"addl	%%edi,%%ecx		/* ebx <- add16+p06 */\n\t"\
			"addl	%%edi,%%edx		/* eax <- add16+p07 */\n\t"\
			"addl	$0x100,%%esi	/* s1p16r */\n\t"\
			"/* MSVC macro assumes add16+p[4,5,7,6] in eax,ebx,ecx,edx, but here get add16+p[4,5,6,7], so replace ecx <-> edx: */\n\t"\
			"/* Do the p0,p4 combo: */\n\t"\
			"movaps	    (%%eax),%%xmm0			\n\t"\
			"movaps	0x10(%%eax),%%xmm1			\n\t"\
			"movaps	%%xmm0,%%xmm2				\n\t"\
			"movaps	%%xmm1,%%xmm3				\n\t"\
			"addpd	    (%%ebx),%%xmm2			\n\t"\
			"addpd	0x10(%%ebx),%%xmm3			\n\t"\
			"subpd	    (%%ebx),%%xmm0			\n\t"\
			"subpd	0x10(%%ebx),%%xmm1			\n\t"\
			"\n\t"\
			"movaps	    (%%edx),%%xmm4			\n\t"\
			"movaps	0x10(%%edx),%%xmm5			\n\t"\
			"movaps	%%xmm4,%%xmm6				\n\t"\
			"movaps	%%xmm5,%%xmm7				\n\t"\
			"addpd	    (%%ecx),%%xmm6			\n\t"\
			"addpd	0x10(%%ecx),%%xmm7			\n\t"\
			"subpd	    (%%ecx),%%xmm4			\n\t"\
			"subpd	0x10(%%ecx),%%xmm5			\n\t"\
			"/* Copy t6r,i into main-array slot add6 */\n\t"\
			"movaps	%%xmm6,0xc0(%%esi)			\n\t"\
			"movaps	%%xmm7,0xd0(%%esi)			\n\t"\
			"/* Copy t7r,i into main-array slot add7 */\n\t"\
			"movaps	%%xmm4,0xe0(%%esi)			\n\t"\
			"movaps	%%xmm5,0xf0(%%esi)			\n\t"\
			"\n\t"\
			"/* Move outputs t5r,i into a(j1,j2+p5), first doing the addsub and mul by ISRT2: */\n\t"\
			"/* Move outputs t7r,i into a(j1,j2+p7), first doing the addsub and mul by ISRT2: */\n\t"\
			"\n\t"\
			"addpd	%%xmm2,%%xmm6				\n\t"\
			"addpd	%%xmm3,%%xmm7				\n\t"\
			"subpd	0xc0(%%esi),%%xmm2			\n\t"\
			"subpd	0xd0(%%esi),%%xmm3			\n\t"\
			"/* Move t6r,i into a(j1,j2+p6) */	\n\t"\
			"movaps	%%xmm2,0xc0(%%esi)			\n\t"\
			"movaps	%%xmm3,0xd0(%%esi)			\n\t"\
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
			"movaps	0x100(%%esi),%%xmm1	/* isrt2 */	\n\t"\
			"subpd	%%xmm3,%%xmm2				\n\t"\
			"mulpd	%%xmm1,%%xmm5				\n\t"\
			"mulpd	%%xmm1,%%xmm2				\n\t"\
			"movaps	%%xmm0,%%xmm3				\n\t"\
			"\n\t"\
			"movaps	%%xmm5,0xa0(%%esi)			\n\t"\
			"movaps	%%xmm4,%%xmm5				\n\t"\
			"addpd	%%xmm4,%%xmm0				\n\t"\
			"movaps	%%xmm2,0xb0(%%esi)			\n\t"\
			"subpd	%%xmm5,%%xmm3				\n\t"\
			"mulpd	%%xmm1,%%xmm0				\n\t"\
			"mulpd	%%xmm1,%%xmm3				\n\t"\
			"movaps	%%xmm0,0xe0(%%esi)			\n\t"\
			"movaps	%%xmm3,0xf0(%%esi)			\n\t"\
			"\n\t"\
			"/* MSVC macro assumes add16+p[2,3,0,1] in eax,ebx,ecx,edx, but here get add16+p[0,1,2,3], so replace eax <-> ecx, ebx <-> edx: */\n\t"\
			"movl	%[__p04],%%edi	\n\t"\
			"shll	$3,%%edi		/* Pointer offset for floating doubles */\n\t"\
			"subl	%%edi,%%eax		/* add0 = add8+p04 */\n\t"\
			"subl	%%edi,%%ebx		/* add1 = add8+p05 */\n\t"\
			"subl	%%edi,%%ecx		/* add3 = add8+p06 */\n\t"\
			"subl	%%edi,%%edx		/* add2 = add8+p07 */\n\t"\
			"\n\t"\
			"movaps	    (%%ecx),%%xmm0			\n\t"\
			"movaps	0x10(%%ecx),%%xmm1			\n\t"\
			"movaps	%%xmm0,%%xmm2				\n\t"\
			"movaps	%%xmm1,%%xmm3				\n\t"\
			"addpd	    (%%edx),%%xmm2			\n\t"\
			"addpd	0x10(%%edx),%%xmm3			\n\t"\
			"subpd	    (%%edx),%%xmm0			\n\t"\
			"subpd	0x10(%%edx),%%xmm1			\n\t"\
			"\n\t"\
			"/* Move   t4r,i into temp(j1+p0) in anticipation of final outputs (t0+t4)r,i which will go there: */\n\t"\
			"movaps	%%xmm6,    (%%esi)			\n\t"\
			"movaps	%%xmm7,0x10(%%esi)			\n\t"\
			"addpd	%%xmm6,%%xmm6				\n\t"\
			"addpd	%%xmm7,%%xmm7				\n\t"\
			"/* Move 2*t4r,i into temp(j1+p4) in anticipation of final outputs (t0-t4)r,i which will go there: */\n\t"\
			"movaps	%%xmm6,0x80(%%esi)			\n\t"\
			"movaps	%%xmm7,0x90(%%esi)			\n\t"\
			"\n\t"\
			"movaps	    (%%eax),%%xmm4			\n\t"\
			"movaps	0x10(%%eax),%%xmm5			\n\t"\
			"movaps	%%xmm4,%%xmm6				\n\t"\
			"movaps	%%xmm5,%%xmm7				\n\t"\
			"addpd	    (%%ebx),%%xmm6			\n\t"\
			"addpd	0x10(%%ebx),%%xmm7			\n\t"\
			"subpd	    (%%ebx),%%xmm4			\n\t"\
			"subpd	0x10(%%ebx),%%xmm5			\n\t"\
			"/* Copy t5,6 into temp(j1+p2) */	\n\t"\
			"movaps	%%xmm6,0x40(%%esi)			\n\t"\
			"movaps	%%xmm7,0x50(%%esi)			\n\t"\
			"addpd	%%xmm2,%%xmm6				\n\t"\
			"addpd	%%xmm3,%%xmm7				\n\t"\
			"subpd	0x40(%%esi),%%xmm2			\n\t"\
			"subpd	0x50(%%esi),%%xmm3			\n\t"\
			"\n\t"\
			"/* Compute and dump first 2 outputs now, in order to free up 2 registers: */\n\t"\
			"addpd	    (%%esi),%%xmm6			\n\t"\
			"addpd	0x10(%%esi),%%xmm7			\n\t"\
			"movaps	%%xmm6,    (%%esi)			\n\t"\
			"movaps	%%xmm7,0x10(%%esi)			\n\t"\
			"\n\t"\
			"subpd	0x80(%%esi),%%xmm6			\n\t"\
			"subpd	0x90(%%esi),%%xmm7			\n\t"\
			"movaps	%%xmm6,0x80(%%esi)			\n\t"\
			"movaps	%%xmm7,0x90(%%esi)			\n\t"\
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
			"addpd	0xd0(%%esi),%%xmm2			\n\t"\
			"subpd	0xc0(%%esi),%%xmm3			\n\t"\
			"subpd	0xd0(%%esi),%%xmm6			\n\t"\
			"addpd	0xc0(%%esi),%%xmm7			\n\t"\
			"movaps	%%xmm2,0xc0(%%esi)			\n\t"\
			"movaps	%%xmm3,0xd0(%%esi)			\n\t"\
			"movaps	%%xmm6,0x40(%%esi)			\n\t"\
			"movaps	%%xmm7,0x50(%%esi)			\n\t"\
			"\n\t"\
			"movaps	%%xmm5,%%xmm2				\n\t"\
			"movaps	%%xmm0,%%xmm6				\n\t"\
			"movaps	%%xmm1,%%xmm3				\n\t"\
			"movaps	%%xmm4,%%xmm7				\n\t"\
			"addpd	0xa0(%%esi),%%xmm5			\n\t"\
			"subpd	0xf0(%%esi),%%xmm0			\n\t"\
			"subpd	0xb0(%%esi),%%xmm1			\n\t"\
			"subpd	0xe0(%%esi),%%xmm4			\n\t"\
			"subpd	0xa0(%%esi),%%xmm2			\n\t"\
			"addpd	0xf0(%%esi),%%xmm6			\n\t"\
			"addpd	0xb0(%%esi),%%xmm3			\n\t"\
			"addpd	0xe0(%%esi),%%xmm7			\n\t"\
			"movaps	%%xmm5,0xe0(%%esi)			\n\t"\
			"movaps	%%xmm0,0xa0(%%esi)			\n\t"\
			"movaps	%%xmm1,0xf0(%%esi)			\n\t"\
			"movaps	%%xmm4,0xb0(%%esi)			\n\t"\
			"movaps	%%xmm2,0x60(%%esi)			\n\t"\
			"movaps	%%xmm6,0x20(%%esi)			\n\t"\
			"movaps	%%xmm3,0x70(%%esi)			\n\t"\
			"movaps	%%xmm7,0x30(%%esi)			\n\t"\
			"\n\t"\
			"movl	%[__cc3],%%edx				\n\t"\
			"movl	 %%esi,%%eax	/* cpy s1p16r */\n\t"\
			"movl	 %%esi,%%ebx	/* cpy s1p16r */\n\t"\
			"movl	 %%esi,%%ecx	/* cpy s1p16r */\n\t"\
			"\n\t"\
		"/* SSE2_RADIX_03_DFT(s1p00r,s1p08r,s1p16r,cc3,s1p00r,s1p16r,s1p08r) */\n\t"\
			"subl	$0x200,%%eax	/* s1p00r */\n\t"\
			"subl	$0x100,%%ebx	/* s1p08r */\n\t"\
			"\n\t"\
			"movaps	    (%%ebx),%%xmm2			\n\t"\
			"movaps	0x10(%%ebx),%%xmm3			\n\t"\
			"movaps	    (%%eax),%%xmm0			\n\t"\
			"movaps	0x10(%%eax),%%xmm1			\n\t"\
			"movaps	    (%%ecx),%%xmm6			\n\t"\
			"movaps	0x10(%%ecx),%%xmm7			\n\t"\
			"movaps	%%xmm2,%%xmm4				\n\t"\
			"movaps	%%xmm3,%%xmm5				\n\t"\
			"\n\t"\
			"addpd	%%xmm6,%%xmm2				\n\t"\
			"addpd	%%xmm7,%%xmm3				\n\t"\
			"subpd	%%xmm6,%%xmm4				\n\t"\
			"subpd	%%xmm7,%%xmm5				\n\t"\
			"addpd	%%xmm2,%%xmm0				\n\t"\
			"addpd	%%xmm3,%%xmm1				\n\t"\
			"movaps	    (%%edx),%%xmm6			\n\t"\
			"movaps	0x10(%%edx),%%xmm7			\n\t"\
			"movaps	%%xmm0,     (%%eax)			\n\t"\
			"movaps	%%xmm1,0x010(%%eax)			\n\t"\
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
			"movaps	%%xmm2,     (%%ecx)	/*b<>c*/\n\t"\
			"movaps	%%xmm3,0x010(%%ecx)			\n\t"\
			"movaps	%%xmm0,     (%%ebx)			\n\t"\
			"movaps	%%xmm1,0x010(%%ebx)			\n\t"\
			"\n\t"\
		"/* SSE2_RADIX_03_DFT(s1p01r,s1p09r,s1p17r,cc3,s1p09r,s1p01r,s1p17r) */\n\t"\
			"addl	$0x020,%%eax	/* s1p01r */\n\t"\
			"addl	$0x020,%%ebx	/* s1p09r */\n\t"\
			"addl	$0x020,%%ecx	/* s1p17r */\n\t"\
			"\n\t"\
			"movaps	    (%%ebx),%%xmm2			\n\t"\
			"movaps	0x10(%%ebx),%%xmm3			\n\t"\
			"movaps	    (%%eax),%%xmm0			\n\t"\
			"movaps	0x10(%%eax),%%xmm1			\n\t"\
			"movaps	    (%%ecx),%%xmm6			\n\t"\
			"movaps	0x10(%%ecx),%%xmm7			\n\t"\
			"movaps	%%xmm2,%%xmm4				\n\t"\
			"movaps	%%xmm3,%%xmm5				\n\t"\
			"\n\t"\
			"addpd	%%xmm6,%%xmm2				\n\t"\
			"addpd	%%xmm7,%%xmm3				\n\t"\
			"subpd	%%xmm6,%%xmm4				\n\t"\
			"subpd	%%xmm7,%%xmm5				\n\t"\
			"addpd	%%xmm2,%%xmm0				\n\t"\
			"addpd	%%xmm3,%%xmm1				\n\t"\
			"movaps	    (%%edx),%%xmm6			\n\t"\
			"movaps	0x10(%%edx),%%xmm7			\n\t"\
			"movaps	%%xmm0,     (%%ebx)	/*a<>b*/\n\t"\
			"movaps	%%xmm1,0x010(%%ebx)			\n\t"\
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
			"movaps	%%xmm2,     (%%eax)			\n\t"\
			"movaps	%%xmm3,0x010(%%eax)			\n\t"\
			"movaps	%%xmm0,     (%%ecx)			\n\t"\
			"movaps	%%xmm1,0x010(%%ecx)			\n\t"\
			"\n\t"\
		"/* SSE2_RADIX_03_DFT(s1p02r,s1p10r,s1p18r,cc3,s1p18r,s1p10r,s1p02r) */\n\t"\
			"addl	$0x020,%%eax	\n\t"\
			"addl	$0x020,%%ebx	\n\t"\
			"addl	$0x020,%%ecx	\n\t"\
			"\n\t"\
			"movaps	    (%%ebx),%%xmm2			\n\t"\
			"movaps	0x10(%%ebx),%%xmm3			\n\t"\
			"movaps	    (%%eax),%%xmm0			\n\t"\
			"movaps	0x10(%%eax),%%xmm1			\n\t"\
			"movaps	    (%%ecx),%%xmm6			\n\t"\
			"movaps	0x10(%%ecx),%%xmm7			\n\t"\
			"movaps	%%xmm2,%%xmm4				\n\t"\
			"movaps	%%xmm3,%%xmm5				\n\t"\
			"\n\t"\
			"addpd	%%xmm6,%%xmm2				\n\t"\
			"addpd	%%xmm7,%%xmm3				\n\t"\
			"subpd	%%xmm6,%%xmm4				\n\t"\
			"subpd	%%xmm7,%%xmm5				\n\t"\
			"addpd	%%xmm2,%%xmm0				\n\t"\
			"addpd	%%xmm3,%%xmm1				\n\t"\
			"movaps	    (%%edx),%%xmm6			\n\t"\
			"movaps	0x10(%%edx),%%xmm7			\n\t"\
			"movaps	%%xmm0,     (%%ecx)	/*a<>c*/\n\t"\
			"movaps	%%xmm1,0x010(%%ecx)			\n\t"\
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
			"movaps	%%xmm2,     (%%ebx)			\n\t"\
			"movaps	%%xmm3,0x010(%%ebx)			\n\t"\
			"movaps	%%xmm0,     (%%eax)			\n\t"\
			"movaps	%%xmm1,0x010(%%eax)			\n\t"\
			"\n\t"\
		"/* SSE2_RADIX_03_DFT(s1p03r,s1p11r,s1p19r,cc3,s1p03r,s1p19r,s1p11r) */\n\t"\
			"addl	$0x020,%%eax	\n\t"\
			"addl	$0x020,%%ebx	\n\t"\
			"addl	$0x020,%%ecx	\n\t"\
			"\n\t"\
			"movaps	    (%%ebx),%%xmm2			\n\t"\
			"movaps	0x10(%%ebx),%%xmm3			\n\t"\
			"movaps	    (%%eax),%%xmm0			\n\t"\
			"movaps	0x10(%%eax),%%xmm1			\n\t"\
			"movaps	    (%%ecx),%%xmm6			\n\t"\
			"movaps	0x10(%%ecx),%%xmm7			\n\t"\
			"movaps	%%xmm2,%%xmm4				\n\t"\
			"movaps	%%xmm3,%%xmm5				\n\t"\
			"\n\t"\
			"addpd	%%xmm6,%%xmm2				\n\t"\
			"addpd	%%xmm7,%%xmm3				\n\t"\
			"subpd	%%xmm6,%%xmm4				\n\t"\
			"subpd	%%xmm7,%%xmm5				\n\t"\
			"addpd	%%xmm2,%%xmm0				\n\t"\
			"addpd	%%xmm3,%%xmm1				\n\t"\
			"movaps	    (%%edx),%%xmm6			\n\t"\
			"movaps	0x10(%%edx),%%xmm7			\n\t"\
			"movaps	%%xmm0,     (%%eax)	/*b<>c*/\n\t"\
			"movaps	%%xmm1,0x010(%%eax)			\n\t"\
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
			"movaps	%%xmm2,     (%%ecx)			\n\t"\
			"movaps	%%xmm3,0x010(%%ecx)			\n\t"\
			"movaps	%%xmm0,     (%%ebx)			\n\t"\
			"movaps	%%xmm1,0x010(%%ebx)			\n\t"\
			"\n\t"\
		"/* SSE2_RADIX_03_DFT(s1p04r,s1p12r,s1p20r,cc3,s1p12r,s1p04r,s1p20r) */\n\t"\
			"addl	$0x020,%%eax	\n\t"\
			"addl	$0x020,%%ebx	\n\t"\
			"addl	$0x020,%%ecx	\n\t"\
			"\n\t"\
			"movaps	    (%%ebx),%%xmm2			\n\t"\
			"movaps	0x10(%%ebx),%%xmm3			\n\t"\
			"movaps	    (%%eax),%%xmm0			\n\t"\
			"movaps	0x10(%%eax),%%xmm1			\n\t"\
			"movaps	    (%%ecx),%%xmm6			\n\t"\
			"movaps	0x10(%%ecx),%%xmm7			\n\t"\
			"movaps	%%xmm2,%%xmm4				\n\t"\
			"movaps	%%xmm3,%%xmm5				\n\t"\
			"\n\t"\
			"addpd	%%xmm6,%%xmm2				\n\t"\
			"addpd	%%xmm7,%%xmm3				\n\t"\
			"subpd	%%xmm6,%%xmm4				\n\t"\
			"subpd	%%xmm7,%%xmm5				\n\t"\
			"addpd	%%xmm2,%%xmm0				\n\t"\
			"addpd	%%xmm3,%%xmm1				\n\t"\
			"movaps	    (%%edx),%%xmm6			\n\t"\
			"movaps	0x10(%%edx),%%xmm7			\n\t"\
			"movaps	%%xmm0,     (%%ebx)	/*a<>b*/\n\t"\
			"movaps	%%xmm1,0x010(%%ebx)			\n\t"\
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
			"movaps	%%xmm2,     (%%eax)			\n\t"\
			"movaps	%%xmm3,0x010(%%eax)			\n\t"\
			"movaps	%%xmm0,     (%%ecx)			\n\t"\
			"movaps	%%xmm1,0x010(%%ecx)			\n\t"\
			"\n\t"\
		"/* SSE2_RADIX_03_DFT(s1p05r,s1p13r,s1p21r,cc3,s1p21r,s1p13r,s1p05r) */\n\t"\
			"addl	$0x020,%%eax	\n\t"\
			"addl	$0x020,%%ebx	\n\t"\
			"addl	$0x020,%%ecx	\n\t"\
			"\n\t"\
			"movaps	    (%%ebx),%%xmm2			\n\t"\
			"movaps	0x10(%%ebx),%%xmm3			\n\t"\
			"movaps	    (%%eax),%%xmm0			\n\t"\
			"movaps	0x10(%%eax),%%xmm1			\n\t"\
			"movaps	    (%%ecx),%%xmm6			\n\t"\
			"movaps	0x10(%%ecx),%%xmm7			\n\t"\
			"movaps	%%xmm2,%%xmm4				\n\t"\
			"movaps	%%xmm3,%%xmm5				\n\t"\
			"\n\t"\
			"addpd	%%xmm6,%%xmm2				\n\t"\
			"addpd	%%xmm7,%%xmm3				\n\t"\
			"subpd	%%xmm6,%%xmm4				\n\t"\
			"subpd	%%xmm7,%%xmm5				\n\t"\
			"addpd	%%xmm2,%%xmm0				\n\t"\
			"addpd	%%xmm3,%%xmm1				\n\t"\
			"movaps	    (%%edx),%%xmm6			\n\t"\
			"movaps	0x10(%%edx),%%xmm7			\n\t"\
			"movaps	%%xmm0,     (%%ecx)	/*a<>c*/\n\t"\
			"movaps	%%xmm1,0x010(%%ecx)			\n\t"\
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
			"movaps	%%xmm2,     (%%ebx)			\n\t"\
			"movaps	%%xmm3,0x010(%%ebx)			\n\t"\
			"movaps	%%xmm0,     (%%eax)			\n\t"\
			"movaps	%%xmm1,0x010(%%eax)			\n\t"\
			"\n\t"\
		"/* SSE2_RADIX_03_DFT(s1p06r,s1p14r,s1p22r,cc3,s1p06r,s1p22r,s1p14r) */\n\t"\
			"addl	$0x020,%%eax	\n\t"\
			"addl	$0x020,%%ebx	\n\t"\
			"addl	$0x020,%%ecx	\n\t"\
			"\n\t"\
			"movaps	    (%%ebx),%%xmm2			\n\t"\
			"movaps	0x10(%%ebx),%%xmm3			\n\t"\
			"movaps	    (%%eax),%%xmm0			\n\t"\
			"movaps	0x10(%%eax),%%xmm1			\n\t"\
			"movaps	    (%%ecx),%%xmm6			\n\t"\
			"movaps	0x10(%%ecx),%%xmm7			\n\t"\
			"movaps	%%xmm2,%%xmm4				\n\t"\
			"movaps	%%xmm3,%%xmm5				\n\t"\
			"\n\t"\
			"addpd	%%xmm6,%%xmm2				\n\t"\
			"addpd	%%xmm7,%%xmm3				\n\t"\
			"subpd	%%xmm6,%%xmm4				\n\t"\
			"subpd	%%xmm7,%%xmm5				\n\t"\
			"addpd	%%xmm2,%%xmm0				\n\t"\
			"addpd	%%xmm3,%%xmm1				\n\t"\
			"movaps	    (%%edx),%%xmm6			\n\t"\
			"movaps	0x10(%%edx),%%xmm7			\n\t"\
			"movaps	%%xmm0,     (%%eax)	/*b<>c*/\n\t"\
			"movaps	%%xmm1,0x010(%%eax)			\n\t"\
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
			"movaps	%%xmm2,     (%%ecx)			\n\t"\
			"movaps	%%xmm3,0x010(%%ecx)			\n\t"\
			"movaps	%%xmm0,     (%%ebx)			\n\t"\
			"movaps	%%xmm1,0x010(%%ebx)			\n\t"\
			"\n\t"\
		"/* SSE2_RADIX_03_DFT(s1p07r,s1p15r,s1p23r,cc3,s1p15r,s1p07r,s1p23r) */\n\t"\
			"addl	$0x020,%%eax	\n\t"\
			"addl	$0x020,%%ebx	\n\t"\
			"addl	$0x020,%%ecx	\n\t"\
			"\n\t"\
			"movaps	    (%%ebx),%%xmm2			\n\t"\
			"movaps	0x10(%%ebx),%%xmm3			\n\t"\
			"movaps	    (%%eax),%%xmm0			\n\t"\
			"movaps	0x10(%%eax),%%xmm1			\n\t"\
			"movaps	    (%%ecx),%%xmm6			\n\t"\
			"movaps	0x10(%%ecx),%%xmm7			\n\t"\
			"movaps	%%xmm2,%%xmm4				\n\t"\
			"movaps	%%xmm3,%%xmm5				\n\t"\
			"\n\t"\
			"addpd	%%xmm6,%%xmm2				\n\t"\
			"addpd	%%xmm7,%%xmm3				\n\t"\
			"subpd	%%xmm6,%%xmm4				\n\t"\
			"subpd	%%xmm7,%%xmm5				\n\t"\
			"addpd	%%xmm2,%%xmm0				\n\t"\
			"addpd	%%xmm3,%%xmm1				\n\t"\
			"movaps	    (%%edx),%%xmm6			\n\t"\
			"movaps	0x10(%%edx),%%xmm7			\n\t"\
			"movaps	%%xmm0,     (%%ebx)	/*a<>b*/\n\t"\
			"movaps	%%xmm1,0x010(%%ebx)			\n\t"\
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
			"movaps	%%xmm2,     (%%eax)			\n\t"\
			"movaps	%%xmm3,0x010(%%eax)			\n\t"\
			"movaps	%%xmm0,     (%%ecx)			\n\t"\
			"movaps	%%xmm1,0x010(%%ecx)			\n\t"\
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
			: "eax","ebx","ecx","edx","edi","esi"		/* Clobbered registers */\
		);\
		}

	#define	SSE2_RADIX24_DIF_NOTWIDDLE(Xadd,Xp01,Xp02,Xp03,Xp04,Xp05,Xp06,Xp07,Xp08,Xp16,Xout,Xisrt2,Xcc3)\
	{\
	__asm__ volatile (\
			"movl	%[__out],%%eax	/* s1p00r */\n\t"\
			"movl	 %%eax,%%ebx	/* s1p00r */\n\t"\
			"movl	 %%eax,%%ecx	/* s1p08r */\n\t"\
			"movl	%[__cc3],%%edx				\n\t"\
			"addl	$0x100,%%ebx	/* s1p08r */\n\t"\
			"addl	$0x200,%%ecx	/* s1p16r */\n\t"\
			"\n\t"\
		"/* On the DIF side it's always an *input* pair to the radix-3 DFT that gets swapped: */\n\t"\
		"/* SSE2_RADIX_03_DFT(s1p00r,s1p16r,s1p08r,cc3,s1p00r,s1p08r,s1p16r) */\n\t"\
			"movaps	    (%%ecx),%%xmm2	/*b<>c*/\n\t"\
			"movaps	0x10(%%ecx),%%xmm3			\n\t"\
			"movaps	    (%%eax),%%xmm0			\n\t"\
			"movaps	0x10(%%eax),%%xmm1			\n\t"\
			"movaps	    (%%ebx),%%xmm6			\n\t"\
			"movaps	0x10(%%ebx),%%xmm7			\n\t"\
			"movaps	%%xmm2,%%xmm4				\n\t"\
			"movaps	%%xmm3,%%xmm5				\n\t"\
			"\n\t"\
			"addpd	%%xmm6,%%xmm2				\n\t"\
			"addpd	%%xmm7,%%xmm3				\n\t"\
			"subpd	%%xmm6,%%xmm4				\n\t"\
			"subpd	%%xmm7,%%xmm5				\n\t"\
			"addpd	%%xmm2,%%xmm0				\n\t"\
			"addpd	%%xmm3,%%xmm1				\n\t"\
			"movaps	    (%%edx),%%xmm6			\n\t"\
			"movaps	0x10(%%edx),%%xmm7			\n\t"\
			"movaps	%%xmm0,     (%%eax)			\n\t"\
			"movaps	%%xmm1,0x010(%%eax)			\n\t"\
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
			"movaps	%%xmm2,     (%%ebx)			\n\t"\
			"movaps	%%xmm3,0x010(%%ebx)			\n\t"\
			"movaps	%%xmm0,     (%%ecx)			\n\t"\
			"movaps	%%xmm1,0x010(%%ecx)			\n\t"\
			"\n\t"\
		"/* SSE2_RADIX_03_DFT(s1p09r,s1p01r,s1p17r,cc3,s1p01r,s1p09r,s1p17r) */\n\t"\
			"addl	$0x020,%%eax	/* s1p01r */\n\t"\
			"addl	$0x020,%%ebx	/* s1p09r */\n\t"\
			"addl	$0x020,%%ecx	/* s1p17r */\n\t"\
			"\n\t"\
			"movaps	    (%%eax),%%xmm2	/*a<>b*/\n\t"\
			"movaps	0x10(%%eax),%%xmm3			\n\t"\
			"movaps	    (%%ebx),%%xmm0			\n\t"\
			"movaps	0x10(%%ebx),%%xmm1			\n\t"\
			"movaps	    (%%ecx),%%xmm6			\n\t"\
			"movaps	0x10(%%ecx),%%xmm7			\n\t"\
			"movaps	%%xmm2,%%xmm4				\n\t"\
			"movaps	%%xmm3,%%xmm5				\n\t"\
			"\n\t"\
			"addpd	%%xmm6,%%xmm2				\n\t"\
			"addpd	%%xmm7,%%xmm3				\n\t"\
			"subpd	%%xmm6,%%xmm4				\n\t"\
			"subpd	%%xmm7,%%xmm5				\n\t"\
			"addpd	%%xmm2,%%xmm0				\n\t"\
			"addpd	%%xmm3,%%xmm1				\n\t"\
			"movaps	    (%%edx),%%xmm6			\n\t"\
			"movaps	0x10(%%edx),%%xmm7			\n\t"\
			"movaps	%%xmm0,     (%%eax)			\n\t"\
			"movaps	%%xmm1,0x010(%%eax)			\n\t"\
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
			"movaps	%%xmm2,     (%%ebx)			\n\t"\
			"movaps	%%xmm3,0x010(%%ebx)			\n\t"\
			"movaps	%%xmm0,     (%%ecx)			\n\t"\
			"movaps	%%xmm1,0x010(%%ecx)			\n\t"\
			"\n\t"\
		"/* SSE2_RADIX_03_DFT(s1p18r,s1p10r,s1p02r,cc3,s1p02r,s1p10r,s1p18r) */\n\t"\
			"addl	$0x020,%%eax	/* s1p02r */\n\t"\
			"addl	$0x020,%%ebx	/* s1p10r */\n\t"\
			"addl	$0x020,%%ecx	/* s1p18r */\n\t"\
			"\n\t"\
			"movaps	    (%%ebx),%%xmm2	/*a<>c*/\n\t"\
			"movaps	0x10(%%ebx),%%xmm3			\n\t"\
			"movaps	    (%%ecx),%%xmm0			\n\t"\
			"movaps	0x10(%%ecx),%%xmm1			\n\t"\
			"movaps	    (%%eax),%%xmm6			\n\t"\
			"movaps	0x10(%%eax),%%xmm7			\n\t"\
			"movaps	%%xmm2,%%xmm4				\n\t"\
			"movaps	%%xmm3,%%xmm5				\n\t"\
			"\n\t"\
			"addpd	%%xmm6,%%xmm2				\n\t"\
			"addpd	%%xmm7,%%xmm3				\n\t"\
			"subpd	%%xmm6,%%xmm4				\n\t"\
			"subpd	%%xmm7,%%xmm5				\n\t"\
			"addpd	%%xmm2,%%xmm0				\n\t"\
			"addpd	%%xmm3,%%xmm1				\n\t"\
			"movaps	    (%%edx),%%xmm6			\n\t"\
			"movaps	0x10(%%edx),%%xmm7			\n\t"\
			"movaps	%%xmm0,     (%%eax)			\n\t"\
			"movaps	%%xmm1,0x010(%%eax)			\n\t"\
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
			"movaps	%%xmm2,     (%%ebx)			\n\t"\
			"movaps	%%xmm3,0x010(%%ebx)			\n\t"\
			"movaps	%%xmm0,     (%%ecx)			\n\t"\
			"movaps	%%xmm1,0x010(%%ecx)			\n\t"\
			"\n\t"\
		"/* SSE2_RADIX_03_DFT(s1p03r,s1p19r,s1p11r,cc3,s1p03r,s1p11r,s1p19r) */\n\t"\
			"addl	$0x020,%%eax	/* s1p03r */\n\t"\
			"addl	$0x020,%%ebx	/* s1p11r */\n\t"\
			"addl	$0x020,%%ecx	/* s1p19r */\n\t"\
			"\n\t"\
			"movaps	    (%%ecx),%%xmm2	/*b<>c*/\n\t"\
			"movaps	0x10(%%ecx),%%xmm3			\n\t"\
			"movaps	    (%%eax),%%xmm0			\n\t"\
			"movaps	0x10(%%eax),%%xmm1			\n\t"\
			"movaps	    (%%ebx),%%xmm6			\n\t"\
			"movaps	0x10(%%ebx),%%xmm7			\n\t"\
			"movaps	%%xmm2,%%xmm4				\n\t"\
			"movaps	%%xmm3,%%xmm5				\n\t"\
			"\n\t"\
			"addpd	%%xmm6,%%xmm2				\n\t"\
			"addpd	%%xmm7,%%xmm3				\n\t"\
			"subpd	%%xmm6,%%xmm4				\n\t"\
			"subpd	%%xmm7,%%xmm5				\n\t"\
			"addpd	%%xmm2,%%xmm0				\n\t"\
			"addpd	%%xmm3,%%xmm1				\n\t"\
			"movaps	    (%%edx),%%xmm6			\n\t"\
			"movaps	0x10(%%edx),%%xmm7			\n\t"\
			"movaps	%%xmm0,     (%%eax)			\n\t"\
			"movaps	%%xmm1,0x010(%%eax)			\n\t"\
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
			"movaps	%%xmm2,     (%%ebx)			\n\t"\
			"movaps	%%xmm3,0x010(%%ebx)			\n\t"\
			"movaps	%%xmm0,     (%%ecx)			\n\t"\
			"movaps	%%xmm1,0x010(%%ecx)			\n\t"\
			"\n\t"\
		"/* SSE2_RADIX_03_DFT(s1p12r,s1p04r,s1p20r,cc3,s1p04r,s1p12r,s1p20r) */\n\t"\
			"addl	$0x020,%%eax	/* s1p04r */\n\t"\
			"addl	$0x020,%%ebx	/* s1p12r */\n\t"\
			"addl	$0x020,%%ecx	/* s1p20r */\n\t"\
			"\n\t"\
			"movaps	    (%%eax),%%xmm2	/*a<>b*/\n\t"\
			"movaps	0x10(%%eax),%%xmm3			\n\t"\
			"movaps	    (%%ebx),%%xmm0			\n\t"\
			"movaps	0x10(%%ebx),%%xmm1			\n\t"\
			"movaps	    (%%ecx),%%xmm6			\n\t"\
			"movaps	0x10(%%ecx),%%xmm7			\n\t"\
			"movaps	%%xmm2,%%xmm4				\n\t"\
			"movaps	%%xmm3,%%xmm5				\n\t"\
			"\n\t"\
			"addpd	%%xmm6,%%xmm2				\n\t"\
			"addpd	%%xmm7,%%xmm3				\n\t"\
			"subpd	%%xmm6,%%xmm4				\n\t"\
			"subpd	%%xmm7,%%xmm5				\n\t"\
			"addpd	%%xmm2,%%xmm0				\n\t"\
			"addpd	%%xmm3,%%xmm1				\n\t"\
			"movaps	    (%%edx),%%xmm6			\n\t"\
			"movaps	0x10(%%edx),%%xmm7			\n\t"\
			"movaps	%%xmm0,     (%%eax)			\n\t"\
			"movaps	%%xmm1,0x010(%%eax)			\n\t"\
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
			"movaps	%%xmm2,     (%%ebx)			\n\t"\
			"movaps	%%xmm3,0x010(%%ebx)			\n\t"\
			"movaps	%%xmm0,     (%%ecx)			\n\t"\
			"movaps	%%xmm1,0x010(%%ecx)			\n\t"\
			"\n\t"\
		"/* SSE2_RADIX_03_DFT(s1p21r,s1p13r,s1p05r,cc3,s1p05r,s1p13r,s1p21r) */\n\t"\
			"addl	$0x020,%%eax	/* s1p05r */\n\t"\
			"addl	$0x020,%%ebx	/* s1p13r */\n\t"\
			"addl	$0x020,%%ecx	/* s1p21r */\n\t"\
			"\n\t"\
			"movaps	    (%%ebx),%%xmm2	/*a<>c*/\n\t"\
			"movaps	0x10(%%ebx),%%xmm3			\n\t"\
			"movaps	    (%%ecx),%%xmm0			\n\t"\
			"movaps	0x10(%%ecx),%%xmm1			\n\t"\
			"movaps	    (%%eax),%%xmm6			\n\t"\
			"movaps	0x10(%%eax),%%xmm7			\n\t"\
			"movaps	%%xmm2,%%xmm4				\n\t"\
			"movaps	%%xmm3,%%xmm5				\n\t"\
			"\n\t"\
			"addpd	%%xmm6,%%xmm2				\n\t"\
			"addpd	%%xmm7,%%xmm3				\n\t"\
			"subpd	%%xmm6,%%xmm4				\n\t"\
			"subpd	%%xmm7,%%xmm5				\n\t"\
			"addpd	%%xmm2,%%xmm0				\n\t"\
			"addpd	%%xmm3,%%xmm1				\n\t"\
			"movaps	    (%%edx),%%xmm6			\n\t"\
			"movaps	0x10(%%edx),%%xmm7			\n\t"\
			"movaps	%%xmm0,     (%%eax)			\n\t"\
			"movaps	%%xmm1,0x010(%%eax)			\n\t"\
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
			"movaps	%%xmm2,     (%%ebx)			\n\t"\
			"movaps	%%xmm3,0x010(%%ebx)			\n\t"\
			"movaps	%%xmm0,     (%%ecx)			\n\t"\
			"movaps	%%xmm1,0x010(%%ecx)			\n\t"\
			"\n\t"\
		"/* SSE2_RADIX_03_DFT(s1p06r,s1p22r,s1p14r,cc3,s1p06r,s1p14r,s1p22r) */\n\t"\
			"addl	$0x020,%%eax	/* s1p06r */\n\t"\
			"addl	$0x020,%%ebx	/* s1p14r */\n\t"\
			"addl	$0x020,%%ecx	/* s1p22r */\n\t"\
			"\n\t"\
			"movaps	    (%%ecx),%%xmm2	/*b<>c*/\n\t"\
			"movaps	0x10(%%ecx),%%xmm3			\n\t"\
			"movaps	    (%%eax),%%xmm0			\n\t"\
			"movaps	0x10(%%eax),%%xmm1			\n\t"\
			"movaps	    (%%ebx),%%xmm6			\n\t"\
			"movaps	0x10(%%ebx),%%xmm7			\n\t"\
			"movaps	%%xmm2,%%xmm4				\n\t"\
			"movaps	%%xmm3,%%xmm5				\n\t"\
			"\n\t"\
			"addpd	%%xmm6,%%xmm2				\n\t"\
			"addpd	%%xmm7,%%xmm3				\n\t"\
			"subpd	%%xmm6,%%xmm4				\n\t"\
			"subpd	%%xmm7,%%xmm5				\n\t"\
			"addpd	%%xmm2,%%xmm0				\n\t"\
			"addpd	%%xmm3,%%xmm1				\n\t"\
			"movaps	    (%%edx),%%xmm6			\n\t"\
			"movaps	0x10(%%edx),%%xmm7			\n\t"\
			"movaps	%%xmm0,     (%%eax)			\n\t"\
			"movaps	%%xmm1,0x010(%%eax)			\n\t"\
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
			"movaps	%%xmm2,     (%%ebx)			\n\t"\
			"movaps	%%xmm3,0x010(%%ebx)			\n\t"\
			"movaps	%%xmm0,     (%%ecx)			\n\t"\
			"movaps	%%xmm1,0x010(%%ecx)			\n\t"\
			"\n\t"\
		"/* SSE2_RADIX_03_DFT(s1p15r,s1p07r,s1p23r,cc3,s1p07r,s1p15r,s1p23r) */\n\t"\
			"addl	$0x020,%%eax	/* s1p07r */\n\t"\
			"addl	$0x020,%%ebx	/* s1p15r */\n\t"\
			"addl	$0x020,%%ecx	/* s1p23r */\n\t"\
			"\n\t"\
			"movaps	    (%%eax),%%xmm2	/*a<>b*/\n\t"\
			"movaps	0x10(%%eax),%%xmm3			\n\t"\
			"movaps	    (%%ebx),%%xmm0			\n\t"\
			"movaps	0x10(%%ebx),%%xmm1			\n\t"\
			"movaps	    (%%ecx),%%xmm6			\n\t"\
			"movaps	0x10(%%ecx),%%xmm7			\n\t"\
			"movaps	%%xmm2,%%xmm4				\n\t"\
			"movaps	%%xmm3,%%xmm5				\n\t"\
			"\n\t"\
			"addpd	%%xmm6,%%xmm2				\n\t"\
			"addpd	%%xmm7,%%xmm3				\n\t"\
			"subpd	%%xmm6,%%xmm4				\n\t"\
			"subpd	%%xmm7,%%xmm5				\n\t"\
			"addpd	%%xmm2,%%xmm0				\n\t"\
			"addpd	%%xmm3,%%xmm1				\n\t"\
			"movaps	    (%%edx),%%xmm6			\n\t"\
			"movaps	0x10(%%edx),%%xmm7			\n\t"\
			"movaps	%%xmm0,     (%%eax)			\n\t"\
			"movaps	%%xmm1,0x010(%%eax)			\n\t"\
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
			"movaps	%%xmm2,     (%%ebx)			\n\t"\
			"movaps	%%xmm3,0x010(%%ebx)			\n\t"\
			"movaps	%%xmm0,     (%%ecx)			\n\t"\
			"movaps	%%xmm1,0x010(%%ecx)			\n\t"\
			"\n\t"\
		"movl	%[__isrt2],%%esi			\n\t"\
			"\n\t"\
		"/* For the radix-8 DIF DFTs, the input offsets always have the same pattern; outputs are permuted */\n\t"\
		"/* SSE2_RADIX8_DIF_0TWIDDLE( i[0-7] = s1p00r + 0x[0a4e82c6]0, o[0-7] = add0 + p[01235476]) */\n\t"\
			"movl	%[__out],%%eax	/* i0 = s1p00r */\n\t"\
			"movl	$0x40  ,%%ebx	/* i2 */	\n\t"\
			"movl	$0x80  ,%%ecx	/* i4 */	\n\t"\
			"movl	$0xc0  ,%%edx	/* i6 */	\n\t"\
			"addl	%%eax,%%ebx					\n\t"\
			"addl	%%eax,%%ecx					\n\t"\
			"addl	%%eax,%%edx					\n\t"\
			"\n\t"\
			"/* Do the p0,p4 combo: */			\n\t"\
			"\n\t"\
			"movaps	    (%%eax),%%xmm0			\n\t"\
			"movaps	    (%%ecx),%%xmm4			\n\t"\
			"movaps	0x10(%%eax),%%xmm1			\n\t"\
			"movaps	0x10(%%ecx),%%xmm5			\n\t"\
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
			"movaps	    (%%ebx),%%xmm4			\n\t"\
			"movaps	0x10(%%ebx),%%xmm5			\n\t"\
			"movaps	%%xmm4,%%xmm6				\n\t"\
			"movaps	%%xmm5,%%xmm7				\n\t"\
			"\n\t"\
			"addpd	     (%%edx),%%xmm4			\n\t"\
			"addpd	0x010(%%edx),%%xmm5			\n\t"\
			"subpd	     (%%edx),%%xmm6			\n\t"\
			"subpd	0x010(%%edx),%%xmm7			\n\t"\
			"\n\t"\
			"subpd	%%xmm4,%%xmm0				\n\t"\
			"subpd	%%xmm7,%%xmm2				\n\t"\
			"subpd	%%xmm5,%%xmm1				\n\t"\
			"subpd	%%xmm6,%%xmm3				\n\t"\
			"movaps	%%xmm0,     (%%ecx)			\n\t"\
			"movaps	%%xmm2,     (%%ebx)			\n\t"\
			"movaps	%%xmm1,0x010(%%ecx)			\n\t"\
			"movaps	%%xmm3,0x010(%%edx)			\n\t"\
			"addpd	%%xmm4,%%xmm4				\n\t"\
			"addpd	%%xmm7,%%xmm7				\n\t"\
			"addpd	%%xmm5,%%xmm5				\n\t"\
			"addpd	%%xmm6,%%xmm6				\n\t"\
			"addpd	%%xmm0,%%xmm4				\n\t"\
			"addpd	%%xmm2,%%xmm7				\n\t"\
			"addpd	%%xmm1,%%xmm5				\n\t"\
			"addpd	%%xmm3,%%xmm6				\n\t"\
			"movaps	%%xmm4,     (%%eax)			\n\t"\
			"movaps	%%xmm7,     (%%edx)			\n\t"\
			"movaps	%%xmm5,0x010(%%eax)			\n\t"\
			"movaps	%%xmm6,0x010(%%ebx)			\n\t"\
			"\n\t"\
			"/* Inline of SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO_B(r8,c1): Do 2nd of 2 radix-4 subtransforms, but now keep outputs in registers: */\n\t"\
			"\n\t"\
			"movl	$0xe0,%%ebx		/* i3 */	\n\t"\
			"movl	$0x20,%%ecx		/* i5 */	\n\t"\
			"movl	$0x60,%%edx		/* i7 */	\n\t"\
			"addl	%%eax,%%ebx					\n\t"\
			"addl	%%eax,%%ecx					\n\t"\
			"addl	%%eax,%%edx					\n\t"\
			"addl	$0xa0,%%eax		/* i1 */	\n\t"\
			"\n\t"\
			"/* Do the p1,p5 combo: */			\n\t"\
			"\n\t"\
			"movaps	    (%%eax),%%xmm0			\n\t"\
			"movaps	    (%%ecx),%%xmm4			\n\t"\
			"movaps	0x10(%%eax),%%xmm1			\n\t"\
			"movaps	0x10(%%ecx),%%xmm5			\n\t"\
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
			"movaps	    (%%ebx),%%xmm4			\n\t"\
			"movaps	0x10(%%ebx),%%xmm5			\n\t"\
			"movaps	%%xmm4,%%xmm6				\n\t"\
			"movaps	%%xmm5,%%xmm7				\n\t"\
			"\n\t"\
			"subpd	    (%%edx),%%xmm4			\n\t"\
			"subpd	0x10(%%edx),%%xmm5			\n\t"\
			"addpd	    (%%edx),%%xmm6			\n\t"\
			"addpd	0x10(%%edx),%%xmm7			\n\t"\
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
			"movaps	%%xmm6,    (%%eax)			\n\t"\
			"movaps	%%xmm7,0x10(%%eax)			\n\t"\
			"\n\t"\
			"movaps	%%xmm2,%%xmm6				\n\t"\
			"movaps	%%xmm5,%%xmm7				\n\t"\
			"subpd	%%xmm4,%%xmm2				\n\t"\
			"subpd	%%xmm3,%%xmm5				\n\t"\
			"addpd	%%xmm6,%%xmm4				\n\t"\
			"addpd	%%xmm7,%%xmm3				\n\t"\
			"\n\t"\
			"movaps	(%%esi),%%xmm6				\n\t"\
			"mulpd	%%xmm6,%%xmm2				\n\t"\
			"mulpd	%%xmm6,%%xmm5				\n\t"\
			"mulpd	%%xmm6,%%xmm4				\n\t"\
			"mulpd	%%xmm6,%%xmm3				\n\t"\
			"\n\t"\
			"/* Inline the contents of SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_A(r0): Combine the 2 radix-4 subtransforms and write outputs to temp array: */\n\t"\
			"/***** t0,1,2,3,4,5,6,7 in __i0,2,4,6 (eax,ebx,ecx,edx) *****/\n\t"\
			"/***** t8,9,a,b,c,d,e,f in __i1,xmm2,4,0,1,5,3; xmm6,7 FREE */\n\t"\
			"subl	$0xa0,%%eax		/* i0 */	\n\t"\
			"movl	%%eax , %%edi				\n\t"\
			"\n\t"\
			"movl	%[__p05],%%eax	\n\t"\
			"movl	%[__p04],%%ebx	\n\t"\
			"movl	%[__p07],%%ecx	\n\t"\
			"movl	%[__p06],%%edx	\n\t"\
			"shll	$3,%%eax		/* Pointer offset for floating doubles */\n\t"\
			"shll	$3,%%ebx		\n\t"\
			"shll	$3,%%ecx		\n\t"\
			"shll	$3,%%edx		\n\t"\
			"addl	%[__add],%%eax	/* o4 */\n\t"\
			"addl	%[__add],%%ebx	/* o5 */\n\t"\
			"addl	%[__add],%%ecx	/* o6 */\n\t"\
			"addl	%[__add],%%edx	/* o7 */\n\t"\
			"\n\t"\
			"movaps	0x40(%%edi),%%xmm6		\n\t"\
			"movaps	0x50(%%edi),%%xmm7		\n\t"\
			"subpd   %%xmm2,%%xmm6				\n\t"\
			"subpd   %%xmm4,%%xmm7				\n\t"\
			"addpd   %%xmm2,%%xmm2				\n\t"\
			"addpd   %%xmm4,%%xmm4				\n\t"\
			"addpd   %%xmm6,%%xmm2				\n\t"\
			"addpd   %%xmm7,%%xmm4				\n\t"\
			"\n\t"\
			"movaps	%%xmm6,    (%%ebx)			\n\t"\
			"movaps	%%xmm7,0x10(%%ebx)			\n\t"\
			"movaps	%%xmm2,    (%%eax)			\n\t"\
			"movaps	%%xmm4,0x10(%%eax)			\n\t"\
			"\n\t"\
			"/* Can't do 2nd set side-by-side with 1st, since only 2 registers (xmm6,7) available for temp storage: */\n\t"\
			"movaps	0xc0(%%edi),%%xmm6		\n\t"\
			"movaps	0xd0(%%edi),%%xmm7		\n\t"\
			"subpd   %%xmm3,%%xmm6				\n\t"\
			"subpd   %%xmm5,%%xmm7				\n\t"\
			"addpd   %%xmm3,%%xmm3				\n\t"\
			"addpd   %%xmm5,%%xmm5				\n\t"\
			"addpd   %%xmm6,%%xmm3				\n\t"\
			"addpd   %%xmm7,%%xmm5				\n\t"\
			"\n\t"\
			"movaps	%%xmm6,    (%%ecx)			\n\t"\
			"movaps	%%xmm7,0x10(%%edx)			\n\t"\
			"movaps	%%xmm3,    (%%edx)			\n\t"\
			"movaps	%%xmm5,0x10(%%ecx)			\n\t"\
			"\n\t"\
			"subl	%[__add],%%ebx	/* p4 */\n\t"\
			"subl	%%ebx,%%eax	/* o1 = add+p1 */\n\t"\
			"subl	%%ebx,%%ecx	/* o2 = add+p3 */\n\t"\
			"subl	%%ebx,%%edx	/* o3 = add+p2 */\n\t"\
			"movl	%[__add],%%ebx	/* o0 */\n\t"\
			"\n\t"\
			"movaps	    (%%edi),%%xmm6			\n\t"\
			"movaps	0x80(%%edi),%%xmm4		\n\t"\
			"movaps	0x10(%%edi),%%xmm7			\n\t"\
			"movaps	0x90(%%edi),%%xmm5		\n\t"\
			"movaps	0xa0(%%edi),%%xmm2		\n\t"\
			"movaps	0xb0(%%edi),%%xmm3		\n\t"\
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
			"movaps	%%xmm6,    (%%eax)	/* a<>b */\n\t"\
			"movaps	%%xmm4,    (%%edx)	/* c<>d */\n\t"\
			"movaps	%%xmm7,0x10(%%eax)			\n\t"\
			"movaps	%%xmm5,0x10(%%ecx)			\n\t"\
			"movaps	%%xmm2,    (%%ebx)			\n\t"\
			"movaps	%%xmm1,    (%%ecx)			\n\t"\
			"movaps	%%xmm3,0x10(%%ebx)			\n\t"\
			"movaps	%%xmm0,0x10(%%edx)			\n\t"\
			"\n\t"\
		"/* SSE2_RADIX8_DIF_0TWIDDLE( i[0-7] = s1p16r + 0x[0a4e82c6]0, o[0-7] = add8 + p[54762310]) */\n\t"\
			"movl	%%edi,%%eax	/* s1p00r */\n\t"\
			"addl	$0x200 ,%%eax	/* i0 = s1p16r */	\n\t"\
			"movl	$0x40  ,%%ebx	/* i2 */\n\t"\
			"movl	$0x80  ,%%ecx	/* i4 */\n\t"\
			"movl	$0xc0  ,%%edx	/* i6 */\n\t"\
			"addl	%%eax,%%ebx				\n\t"\
			"addl	%%eax,%%ecx				\n\t"\
			"addl	%%eax,%%edx				\n\t"\
			"\n\t"\
			"/* Do the p0,p4 combo: */		\n\t"\
			"\n\t"\
			"movaps	    (%%eax),%%xmm0		\n\t"\
			"movaps	    (%%ecx),%%xmm4		\n\t"\
			"movaps	0x10(%%eax),%%xmm1		\n\t"\
			"movaps	0x10(%%ecx),%%xmm5		\n\t"\
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
			"movaps	    (%%ebx),%%xmm4		\n\t"\
			"movaps	0x10(%%ebx),%%xmm5		\n\t"\
			"movaps	%%xmm4,%%xmm6			\n\t"\
			"movaps	%%xmm5,%%xmm7			\n\t"\
			"\n\t"\
			"addpd	     (%%edx),%%xmm4		\n\t"\
			"addpd	0x010(%%edx),%%xmm5		\n\t"\
			"subpd	     (%%edx),%%xmm6		\n\t"\
			"subpd	0x010(%%edx),%%xmm7		\n\t"\
			"\n\t"\
			"subpd	%%xmm4,%%xmm0			\n\t"\
			"subpd	%%xmm7,%%xmm2			\n\t"\
			"subpd	%%xmm5,%%xmm1			\n\t"\
			"subpd	%%xmm6,%%xmm3			\n\t"\
			"movaps	%%xmm0,     (%%ecx)		\n\t"\
			"movaps	%%xmm2,     (%%ebx)		\n\t"\
			"movaps	%%xmm1,0x010(%%ecx)		\n\t"\
			"movaps	%%xmm3,0x010(%%edx)		\n\t"\
			"addpd	%%xmm4,%%xmm4			\n\t"\
			"addpd	%%xmm7,%%xmm7			\n\t"\
			"addpd	%%xmm5,%%xmm5			\n\t"\
			"addpd	%%xmm6,%%xmm6			\n\t"\
			"addpd	%%xmm0,%%xmm4			\n\t"\
			"addpd	%%xmm2,%%xmm7			\n\t"\
			"addpd	%%xmm1,%%xmm5			\n\t"\
			"addpd	%%xmm3,%%xmm6			\n\t"\
			"movaps	%%xmm4,     (%%eax)		\n\t"\
			"movaps	%%xmm7,     (%%edx)		\n\t"\
			"movaps	%%xmm5,0x010(%%eax)		\n\t"\
			"movaps	%%xmm6,0x010(%%ebx)		\n\t"\
			"\n\t"\
			"/* Inline of SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO_B(r8,c1): Do 2nd of 2 radix-4 subtransforms, but now keep outputs in registers: */\n\t"\
			"\n\t"\
			"movl	$0xe0,%%ebx		/* i3 */\n\t"\
			"movl	$0x20,%%ecx		/* i5 */\n\t"\
			"movl	$0x60,%%edx		/* i7 */\n\t"\
			"addl	%%eax,%%ebx				\n\t"\
			"addl	%%eax,%%ecx				\n\t"\
			"addl	%%eax,%%edx				\n\t"\
			"addl	$0xa0,%%eax		/* i1 */\n\t"\
			"\n\t"\
			"/* Do the p1,p5 combo: */		\n\t"\
			"\n\t"\
			"movaps	    (%%eax),%%xmm0		\n\t"\
			"movaps	    (%%ecx),%%xmm4		\n\t"\
			"movaps	0x10(%%eax),%%xmm1		\n\t"\
			"movaps	0x10(%%ecx),%%xmm5		\n\t"\
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
			"movaps	    (%%ebx),%%xmm4		\n\t"\
			"movaps	0x10(%%ebx),%%xmm5		\n\t"\
			"movaps	%%xmm4,%%xmm6			\n\t"\
			"movaps	%%xmm5,%%xmm7			\n\t"\
			"\n\t"\
			"subpd	    (%%edx),%%xmm4		\n\t"\
			"subpd	0x10(%%edx),%%xmm5		\n\t"\
			"addpd	    (%%edx),%%xmm6		\n\t"\
			"addpd	0x10(%%edx),%%xmm7		\n\t"\
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
			"movaps	%%xmm6,    (%%eax)		\n\t"\
			"movaps	%%xmm7,0x10(%%eax)		\n\t"\
			"\n\t"\
			"movaps	%%xmm2,%%xmm6			\n\t"\
			"movaps	%%xmm5,%%xmm7			\n\t"\
			"subpd	%%xmm4,%%xmm2			\n\t"\
			"subpd	%%xmm3,%%xmm5			\n\t"\
			"addpd	%%xmm6,%%xmm4			\n\t"\
			"addpd	%%xmm7,%%xmm3			\n\t"\
			"\n\t"\
			"movaps	(%%esi),%%xmm6			\n\t"\
			"mulpd	%%xmm6,%%xmm2			\n\t"\
			"mulpd	%%xmm6,%%xmm5			\n\t"\
			"mulpd	%%xmm6,%%xmm4			\n\t"\
			"mulpd	%%xmm6,%%xmm3			\n\t"\
			"\n\t"\
			"/* Inline the contents of SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_A(r0): Combine the 2 radix-4 subtransforms and write outputs to temp array: */\n\t"\
			"/***** t0,1,2,3,4,5,6,7 in __i0,2,4,6 (eax,ebx,ecx,edx) *****/\n\t"\
			"/***** t8,9,a,b,c,d,e,f in __i1,xmm2,4,0,1,5,3; xmm6,7 FREE */\n\t"\
			"subl	$0xa0,%%eax		/* i0 */	\n\t"\
			"movl	%%eax , %%edi			\n\t"\
			"\n\t"\
			"movl	%[__p08],%%edx			\n\t"\
			"shll	$3,%%edx		/* Pointer offset for floating doubles */\n\t"\
			"movl	%[__p02],%%eax			\n\t"\
			"movl	%[__p03],%%ebx			\n\t"\
			"movl	%[__p01],%%ecx			\n\t"\
			"addl	%[__add],%%edx	/* o7 */\n\t"\
			"shll	$3,%%eax		/* Pointer offset for floating doubles */\n\t"\
			"shll	$3,%%ebx				\n\t"\
			"shll	$3,%%ecx				\n\t"\
			"addl	%%edx,%%eax	/* o4 */	\n\t"\
			"addl	%%edx,%%ebx	/* o5 */	\n\t"\
			"addl	%%edx,%%ecx	/* o6 */	\n\t"\
			"\n\t"\
			"movaps	0x40(%%edi),%%xmm6		\n\t"\
			"movaps	0x50(%%edi),%%xmm7		\n\t"\
			"subpd   %%xmm2,%%xmm6			\n\t"\
			"subpd   %%xmm4,%%xmm7			\n\t"\
			"addpd   %%xmm2,%%xmm2			\n\t"\
			"addpd   %%xmm4,%%xmm4			\n\t"\
			"addpd   %%xmm6,%%xmm2			\n\t"\
			"addpd   %%xmm7,%%xmm4			\n\t"\
			"\n\t"\
			"movaps	%%xmm6,    (%%ebx)		\n\t"\
			"movaps	%%xmm7,0x10(%%ebx)		\n\t"\
			"movaps	%%xmm2,    (%%eax)		\n\t"\
			"movaps	%%xmm4,0x10(%%eax)		\n\t"\
			"\n\t"\
			"/* Can't do 2nd set side-by-side with 1st, since only 2 registers (xmm6,7) available for temp storage: */\n\t"\
			"movaps	0xc0(%%edi),%%xmm6		\n\t"\
			"movaps	0xd0(%%edi),%%xmm7		\n\t"\
			"subpd   %%xmm3,%%xmm6			\n\t"\
			"subpd   %%xmm5,%%xmm7			\n\t"\
			"addpd   %%xmm3,%%xmm3			\n\t"\
			"addpd   %%xmm5,%%xmm5			\n\t"\
			"addpd   %%xmm6,%%xmm3			\n\t"\
			"addpd   %%xmm7,%%xmm5			\n\t"\
			"\n\t"\
			"movaps	%%xmm6,    (%%ecx)		\n\t"\
			"movaps	%%xmm7,0x10(%%edx)		\n\t"\
			"movaps	%%xmm3,    (%%edx)		\n\t"\
			"movaps	%%xmm5,0x10(%%ecx)		\n\t"\
			"\n\t"\
			"movl	%[__p04],%%esi			\n\t"\
			"shll	$3,%%esi		/* Pointer offset for floating doubles */\n\t"\
			"addl	%%esi,%%eax	/* o3 = add8+p6 */\n\t"\
			"addl	%%esi,%%ebx	/* o2 = add8+p7 */\n\t"\
			"addl	%%esi,%%ecx	/* o0 = add8+p5 */\n\t"\
			"addl	%%esi,%%edx	/* o1 = add8+p4 */\n\t"\
			"\n\t"\
		"movl	%[__isrt2],%%esi	/* Restore isrt2 ptr */\n\t"\
			"\n\t"\
			"movaps	    (%%edi),%%xmm6		\n\t"\
			"movaps	0x80(%%edi),%%xmm4		\n\t"\
			"movaps	0x10(%%edi),%%xmm7		\n\t"\
			"movaps	0x90(%%edi),%%xmm5		\n\t"\
			"movaps	0xa0(%%edi),%%xmm2		\n\t"\
			"movaps	0xb0(%%edi),%%xmm3		\n\t"\
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
			"movaps	%%xmm6,    (%%edx)	/* abcd -> cdba */\n\t"\
			"movaps	%%xmm4,    (%%ebx)		\n\t"\
			"movaps	%%xmm7,0x10(%%edx)		\n\t"\
			"movaps	%%xmm5,0x10(%%eax)		\n\t"\
			"movaps	%%xmm2,    (%%ecx)		\n\t"\
			"movaps	%%xmm1,    (%%eax)		\n\t"\
			"movaps	%%xmm3,0x10(%%ecx)		\n\t"\
			"movaps	%%xmm0,0x10(%%ebx)		\n\t"\
			"\n\t"\
		"/* SSE2_RADIX8_DIF_0TWIDDLE( i[0-7] = s1p08r + 0x[0a4e82c6]0, o[0-7] = add16+ p[23107645]) */\n\t"\
			"movl	%%edi,%%eax	/* s1p00r */\n\t"\
			"subl	$0x100 ,%%eax	/* i0 = s1p08r */	\n\t"\
			"movl	$0x40  ,%%ebx	/* i2 */\n\t"\
			"movl	$0x80  ,%%ecx	/* i4 */\n\t"\
			"movl	$0xc0  ,%%edx	/* i6 */\n\t"\
			"addl	%%eax,%%ebx				\n\t"\
			"addl	%%eax,%%ecx				\n\t"\
			"addl	%%eax,%%edx				\n\t"\
			"\n\t"\
			"/* Do the p0,p4 combo: */		\n\t"\
			"\n\t"\
			"movaps	    (%%eax),%%xmm0		\n\t"\
			"movaps	    (%%ecx),%%xmm4		\n\t"\
			"movaps	0x10(%%eax),%%xmm1		\n\t"\
			"movaps	0x10(%%ecx),%%xmm5		\n\t"\
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
			"movaps	    (%%ebx),%%xmm4		\n\t"\
			"movaps	0x10(%%ebx),%%xmm5		\n\t"\
			"movaps	%%xmm4,%%xmm6			\n\t"\
			"movaps	%%xmm5,%%xmm7			\n\t"\
			"\n\t"\
			"addpd	     (%%edx),%%xmm4		\n\t"\
			"addpd	0x010(%%edx),%%xmm5		\n\t"\
			"subpd	     (%%edx),%%xmm6		\n\t"\
			"subpd	0x010(%%edx),%%xmm7		\n\t"\
			"\n\t"\
			"subpd	%%xmm4,%%xmm0			\n\t"\
			"subpd	%%xmm7,%%xmm2			\n\t"\
			"subpd	%%xmm5,%%xmm1			\n\t"\
			"subpd	%%xmm6,%%xmm3			\n\t"\
			"movaps	%%xmm0,     (%%ecx)		\n\t"\
			"movaps	%%xmm2,     (%%ebx)		\n\t"\
			"movaps	%%xmm1,0x010(%%ecx)		\n\t"\
			"movaps	%%xmm3,0x010(%%edx)		\n\t"\
			"addpd	%%xmm4,%%xmm4			\n\t"\
			"addpd	%%xmm7,%%xmm7			\n\t"\
			"addpd	%%xmm5,%%xmm5			\n\t"\
			"addpd	%%xmm6,%%xmm6			\n\t"\
			"addpd	%%xmm0,%%xmm4			\n\t"\
			"addpd	%%xmm2,%%xmm7			\n\t"\
			"addpd	%%xmm1,%%xmm5			\n\t"\
			"addpd	%%xmm3,%%xmm6			\n\t"\
			"movaps	%%xmm4,     (%%eax)		\n\t"\
			"movaps	%%xmm7,     (%%edx)		\n\t"\
			"movaps	%%xmm5,0x010(%%eax)		\n\t"\
			"movaps	%%xmm6,0x010(%%ebx)		\n\t"\
			"\n\t"\
			"/* Inline of SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO_B(r8,c1): Do 2nd of 2 radix-4 subtransforms, but now keep outputs in registers: */\n\t"\
			"\n\t"\
			"movl	$0xe0,%%ebx		/* i3 */\n\t"\
			"movl	$0x20,%%ecx		/* i5 */\n\t"\
			"movl	$0x60,%%edx		/* i7 */\n\t"\
			"addl	%%eax,%%ebx				\n\t"\
			"addl	%%eax,%%ecx				\n\t"\
			"addl	%%eax,%%edx				\n\t"\
			"addl	$0xa0,%%eax		/* i1 */\n\t"\
			"\n\t"\
			"/* Do the p1,p5 combo: */		\n\t"\
			"\n\t"\
			"movaps	    (%%eax),%%xmm0		\n\t"\
			"movaps	    (%%ecx),%%xmm4		\n\t"\
			"movaps	0x10(%%eax),%%xmm1		\n\t"\
			"movaps	0x10(%%ecx),%%xmm5		\n\t"\
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
			"movaps	    (%%ebx),%%xmm4		\n\t"\
			"movaps	0x10(%%ebx),%%xmm5		\n\t"\
			"movaps	%%xmm4,%%xmm6			\n\t"\
			"movaps	%%xmm5,%%xmm7			\n\t"\
			"\n\t"\
			"subpd	    (%%edx),%%xmm4		\n\t"\
			"subpd	0x10(%%edx),%%xmm5		\n\t"\
			"addpd	    (%%edx),%%xmm6		\n\t"\
			"addpd	0x10(%%edx),%%xmm7		\n\t"\
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
			"movaps	%%xmm6,    (%%eax)		\n\t"\
			"movaps	%%xmm7,0x10(%%eax)		\n\t"\
			"\n\t"\
			"movaps	%%xmm2,%%xmm6			\n\t"\
			"movaps	%%xmm5,%%xmm7			\n\t"\
			"subpd	%%xmm4,%%xmm2			\n\t"\
			"subpd	%%xmm3,%%xmm5			\n\t"\
			"addpd	%%xmm6,%%xmm4			\n\t"\
			"addpd	%%xmm7,%%xmm3			\n\t"\
			"\n\t"\
			"movaps	(%%esi),%%xmm6			\n\t"\
			"mulpd	%%xmm6,%%xmm2			\n\t"\
			"mulpd	%%xmm6,%%xmm5			\n\t"\
			"mulpd	%%xmm6,%%xmm4			\n\t"\
			"mulpd	%%xmm6,%%xmm3			\n\t"\
			"\n\t"\
			"/* Inline the contents of SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_A(r0): Combine the 2 radix-4 subtransforms and write outputs to temp array: */\n\t"\
			"/***** t0,1,2,3,4,5,6,7 in __i0,2,4,6 (eax,ebx,ecx,edx) *****/\n\t"\
			"/***** t8,9,a,b,c,d,e,f in __i1,xmm2,4,0,1,5,3; xmm6,7 FREE */\n\t"\
			"subl	$0xa0,%%eax		/* i0 */	\n\t"\
			"movl	%%eax , %%edi			\n\t"\
			"\n\t"\
			"movl	%[__p04],%%ecx			\n\t"\
			"addl	%[__p16],%%ecx			\n\t"\
			"movl	%[__p03],%%eax			\n\t"\
			"movl	%[__p02],%%ebx			\n\t"\
			"movl	%[__p01],%%edx			\n\t"\
			"shll	$3,%%ecx				\n\t"\
			"shll	$3,%%eax		/* Pointer offset for floating doubles */\n\t"\
			"shll	$3,%%ebx				\n\t"\
			"shll	$3,%%edx				\n\t"\
			"addl	%[__add],%%ecx	/* o6 = add16+p4 */\n\t"\
			"addl	%%ecx,%%eax		/* o4 = add16+p7 */\n\t"\
			"addl	%%ecx,%%ebx		/* o5 = add16+p6 */\n\t"\
			"addl	%%ecx,%%edx		/* o7 = add16+p5 */\n\t"\
			"\n\t"\
			"movaps	0x40(%%edi),%%xmm6		\n\t"\
			"movaps	0x50(%%edi),%%xmm7		\n\t"\
			"subpd   %%xmm2,%%xmm6			\n\t"\
			"subpd   %%xmm4,%%xmm7			\n\t"\
			"addpd   %%xmm2,%%xmm2			\n\t"\
			"addpd   %%xmm4,%%xmm4			\n\t"\
			"addpd   %%xmm6,%%xmm2			\n\t"\
			"addpd   %%xmm7,%%xmm4			\n\t"\
			"\n\t"\
			"movaps	%%xmm6,    (%%ebx)		\n\t"\
			"movaps	%%xmm7,0x10(%%ebx)		\n\t"\
			"movaps	%%xmm2,    (%%eax)		\n\t"\
			"movaps	%%xmm4,0x10(%%eax)		\n\t"\
			"\n\t"\
			"/* Can't do 2nd set side-by-side with 1st, since only 2 registers (xmm6,7) available for temp storage: */\n\t"\
			"movaps	0xc0(%%edi),%%xmm6		\n\t"\
			"movaps	0xd0(%%edi),%%xmm7		\n\t"\
			"subpd   %%xmm3,%%xmm6			\n\t"\
			"subpd   %%xmm5,%%xmm7			\n\t"\
			"addpd   %%xmm3,%%xmm3			\n\t"\
			"addpd   %%xmm5,%%xmm5			\n\t"\
			"addpd   %%xmm6,%%xmm3			\n\t"\
			"addpd   %%xmm7,%%xmm5			\n\t"\
			"\n\t"\
			"movaps	%%xmm6,    (%%ecx)		\n\t"\
			"movaps	%%xmm7,0x10(%%edx)		\n\t"\
			"movaps	%%xmm3,    (%%edx)		\n\t"\
			"movaps	%%xmm5,0x10(%%ecx)		\n\t"\
			"\n\t"\
			"movl	%[__p04],%%esi		/* overwrite __isrt2 ptr */\n\t"\
			"shll	$3,%%esi		/* Pointer offset for floating doubles */\n\t"\
			"subl	%%esi,%%eax	/* o1 = add16+p3 */\n\t"\
			"subl	%%esi,%%ebx	/* o0 = add16+p2 */\n\t"\
			"subl	%%esi,%%ecx	/* o3 = add16+p0 */\n\t"\
			"subl	%%esi,%%edx	/* o2 = add16+p1 */\n\t"\
			"\n\t"\
			"movaps	    (%%edi),%%xmm6		\n\t"\
			"movaps	0x80(%%edi),%%xmm4		\n\t"\
			"movaps	0x10(%%edi),%%xmm7		\n\t"\
			"movaps	0x90(%%edi),%%xmm5		\n\t"\
			"movaps	0xa0(%%edi),%%xmm2		\n\t"\
			"movaps	0xb0(%%edi),%%xmm3		\n\t"\
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
			"movaps	%%xmm6,    (%%eax)	/* abcd -> badc */\n\t"\
			"movaps	%%xmm4,    (%%edx)		\n\t"\
			"movaps	%%xmm7,0x10(%%eax)		\n\t"\
			"movaps	%%xmm5,0x10(%%ecx)		\n\t"\
			"movaps	%%xmm2,    (%%ebx)		\n\t"\
			"movaps	%%xmm1,    (%%ecx)		\n\t"\
			"movaps	%%xmm3,0x10(%%ebx)		\n\t"\
			"movaps	%%xmm0,0x10(%%edx)		\n\t"\
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
			: "eax","ebx","ecx","edx","edi","esi"		/* Clobbered registers */\
		);\
		}

#endif	/* radix24_ditN_cy_dif1_gcc_h_included */

