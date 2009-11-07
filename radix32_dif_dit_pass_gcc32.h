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
#ifndef radix32_dif_dit_pass_gcc_h_included
#define radix32_dif_dit_pass_gcc_h_included

	/*
	For GCC-macro version of this, use that isrt2 + 0x010,0x030,0x050 = cc0,cc1,cc3,
	and isrt2 + 0x070,0x0f0,0x170,0x1f0,0x270,0x2f0,0x370,0x3f0 = c00,04,02,06,01,05,03,07
	in order to reduce number of args to <= the GCC-allowed maximum of 30:
	*/
	#define SSE2_RADIX32_DIF_TWIDDLE(Xadd0,Xp01,Xp02,Xp03,Xp04,Xp08,Xp0C,Xp10,Xp18,Xr00)\
	{\
	__asm__ volatile (\
		"/*...Block 1: */	\n\t"\
		"movl	%[__add0],%%eax	\n\t"\
		"movl	%[__p08],%%ebx	\n\t"\
		"movl	%[__p10],%%ecx	\n\t"\
		"movl	%[__p18],%%edx	\n\t"\
		"shll	$3,%%ebx	\n\t"\
		"shll	$3,%%ecx	\n\t"\
		"shll	$3,%%edx	\n\t"\
		"addl	%%eax,%%ebx	\n\t"\
		"addl	%%eax,%%ecx	\n\t"\
		"addl	%%eax,%%edx	\n\t"\
		"/* SSE2_RADIX4_DIF_4TWIDDLE_A(r00,c10) */	\n\t"\
		"movl	%[__r00],%%esi	\n\t"\
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	    (%%ecx),%%xmm4	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"movaps	0x10(%%ecx),%%xmm5	\n\t"\
		"movaps	0x490(%%esi),%%xmm2	/* c10 */	\n\t"\
		"movaps	0x4a0(%%esi),%%xmm3	\n\t"\
		"movaps	%%xmm4,%%xmm6	\n\t"\
		"movaps	%%xmm5,%%xmm7	\n\t"\
		"mulpd	%%xmm2,%%xmm4	\n\t"\
		"mulpd	%%xmm2,%%xmm5	\n\t"\
		"mulpd	%%xmm3,%%xmm6	\n\t"\
		"movaps	%%xmm0,%%xmm2	\n\t"\
		"mulpd	%%xmm3,%%xmm7	\n\t"\
		"addpd	%%xmm6,%%xmm5	\n\t"\
		"movaps	%%xmm1,%%xmm3	\n\t"\
		"subpd	%%xmm7,%%xmm4	\n\t"\
		"addpd	%%xmm4,%%xmm0	\n\t"\
		"addpd	%%xmm5,%%xmm1	\n\t"\
		"subpd	%%xmm4,%%xmm2	\n\t"\
		"subpd	%%xmm5,%%xmm3	\n\t"\
		"movaps	    (%%edx),%%xmm4	\n\t"\
		"movaps	0x10(%%edx),%%xmm5	\n\t"\
		"movaps	    (%%edx),%%xmm6	\n\t"\
		"movaps	0x10(%%edx),%%xmm7	\n\t"\
		"mulpd	0x4d0(%%esi),%%xmm4	/* c18 */	\n\t"\
		"mulpd	0x4d0(%%esi),%%xmm5	\n\t"\
		"mulpd	0x4e0(%%esi),%%xmm6	\n\t"\
		"mulpd	0x4e0(%%esi),%%xmm7	\n\t"\
		"addpd	%%xmm6,%%xmm5	\n\t"\
		"subpd	%%xmm7,%%xmm4	\n\t"\
		"movaps	%%xmm5,0x010(%%esi)	\n\t"\
		"movaps	%%xmm4,     (%%esi)	/* r00 */	\n\t"\
		"movaps	    (%%ebx),%%xmm4	\n\t"\
		"movaps	0x10(%%ebx),%%xmm5	\n\t"\
		"movaps	    (%%ebx),%%xmm6	\n\t"\
		"movaps	0x10(%%ebx),%%xmm7	\n\t"\
		"mulpd	0x4b0(%%esi),%%xmm4	/* c08 */	\n\t"\
		"mulpd	0x4b0(%%esi),%%xmm5	\n\t"\
		"mulpd	0x4c0(%%esi),%%xmm6	\n\t"\
		"mulpd	0x4c0(%%esi),%%xmm7	\n\t"\
		"addpd	%%xmm6,%%xmm5	\n\t"\
		"subpd	%%xmm7,%%xmm4	\n\t"\
		"movaps	%%xmm5,%%xmm7	\n\t"\
		"movaps	%%xmm4,%%xmm6	\n\t"\
		"subpd	    (%%esi),%%xmm4	/* r00 */	\n\t"\
		"subpd	0x10(%%esi),%%xmm5	\n\t"\
		"addpd	    (%%esi),%%xmm6	\n\t"\
		"addpd	0x10(%%esi),%%xmm7	\n\t"\
		"subpd	%%xmm6,%%xmm0	\n\t"\
		"subpd	%%xmm5,%%xmm2	\n\t"\
		"subpd	%%xmm7,%%xmm1	\n\t"\
		"subpd	%%xmm4,%%xmm3	\n\t"\
		"movaps	%%xmm0,0x040(%%esi)	\n\t"\
		"movaps	%%xmm2,0x020(%%esi)	\n\t"\
		"movaps	%%xmm1,0x050(%%esi)	\n\t"\
		"movaps	%%xmm3,0x070(%%esi)	\n\t"\
		"addpd	%%xmm6,%%xmm6	\n\t"\
		"addpd	%%xmm5,%%xmm5	\n\t"\
		"addpd	%%xmm7,%%xmm7	\n\t"\
		"addpd	%%xmm4,%%xmm4	\n\t"\
		"addpd	%%xmm0,%%xmm6	\n\t"\
		"addpd	%%xmm2,%%xmm5	\n\t"\
		"addpd	%%xmm1,%%xmm7	\n\t"\
		"addpd	%%xmm3,%%xmm4	\n\t"\
		"movaps	%%xmm6,     (%%esi)	\n\t"\
		"movaps	%%xmm5,0x060(%%esi)	\n\t"\
		"movaps	%%xmm7,0x010(%%esi)	\n\t"\
		"movaps	%%xmm4,0x030(%%esi)	\n\t"\
		"/***************************************/	\n\t"\
		"movl	%[__p04],%%edi	\n\t"\
		"shll	$3,%%edi	\n\t"\
		"addl	%%edi,%%eax	/* &a[j1+p4] */	\n\t"\
		"addl	%%edi,%%ebx	\n\t"\
		"addl	%%edi,%%ecx	\n\t"\
		"addl	%%edi,%%edx	\n\t"\
		"/* SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO_B(r08,c04) */	\n\t"\
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	    (%%ecx),%%xmm4	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"movaps	0x10(%%ecx),%%xmm5	\n\t"\
		"movaps	0x4f0(%%esi),%%xmm6	/* c04 */	\n\t"\
		"movaps	0x500(%%esi),%%xmm7	\n\t"\
		"movaps	%%xmm0,%%xmm2	\n\t"\
		"movaps	%%xmm1,%%xmm3	\n\t"\
		"mulpd	%%xmm6,%%xmm0	\n\t"\
		"mulpd	%%xmm6,%%xmm1	\n\t"\
		"mulpd	%%xmm7,%%xmm2	\n\t"\
		"mulpd	%%xmm7,%%xmm3	\n\t"\
		"movaps	%%xmm4,%%xmm6	\n\t"\
		"addpd	%%xmm2,%%xmm1	\n\t"\
		"movaps	%%xmm5,%%xmm7	\n\t"\
		"mulpd	0x510(%%esi),%%xmm4	/* c14 */	\n\t"\
		"subpd	%%xmm3,%%xmm0	\n\t"\
		"mulpd	0x510(%%esi),%%xmm5	\n\t"\
		"mulpd	0x520(%%esi),%%xmm6	\n\t"\
		"movaps	%%xmm0,%%xmm2	\n\t"\
		"mulpd	0x520(%%esi),%%xmm7	\n\t"\
		"addpd	%%xmm6,%%xmm5	\n\t"\
		"movaps	%%xmm1,%%xmm3	\n\t"\
		"subpd	%%xmm7,%%xmm4	\n\t"\
		"addpd	%%xmm4,%%xmm0	\n\t"\
		"addpd	%%xmm5,%%xmm1	\n\t"\
		"subpd	%%xmm4,%%xmm2	\n\t"\
		"subpd	%%xmm5,%%xmm3	\n\t"\
		"movaps	    (%%edx),%%xmm4	\n\t"\
		"movaps	0x10(%%edx),%%xmm5	\n\t"\
		"movaps	    (%%edx),%%xmm6	\n\t"\
		"movaps	0x10(%%edx),%%xmm7	\n\t"\
		"mulpd	0x550(%%esi),%%xmm4	/* c1C */	\n\t"\
		"mulpd	0x550(%%esi),%%xmm5	\n\t"\
		"mulpd	0x560(%%esi),%%xmm6	\n\t"\
		"mulpd	0x560(%%esi),%%xmm7	\n\t"\
		"addl	$0x80,%%esi	/* r08 */	\n\t"\
		"addpd	%%xmm6,%%xmm5	\n\t"\
		"subpd	%%xmm7,%%xmm4	\n\t"\
		"movaps	%%xmm5,0x010(%%esi)	\n\t"\
		"movaps	%%xmm4,     (%%esi)	\n\t"\
		"movaps	    (%%ebx),%%xmm4	\n\t"\
		"movaps	0x10(%%ebx),%%xmm5	\n\t"\
		"movaps	    (%%ebx),%%xmm6	\n\t"\
		"movaps	0x10(%%ebx),%%xmm7	\n\t"\
		"mulpd	0x4b0(%%esi),%%xmm4	/* c0C */	\n\t"\
		"mulpd	0x4b0(%%esi),%%xmm5	\n\t"\
		"mulpd	0x4c0(%%esi),%%xmm6	\n\t"\
		"mulpd	0x4c0(%%esi),%%xmm7	\n\t"\
		"addpd	%%xmm6,%%xmm5	\n\t"\
		"subpd	%%xmm7,%%xmm4	\n\t"\
		"movaps	%%xmm5,%%xmm7	\n\t"\
		"movaps	%%xmm4,%%xmm6	\n\t"\
		"subpd	     (%%esi),%%xmm4	\n\t"\
		"subpd	0x010(%%esi),%%xmm5	\n\t"\
		"addpd	     (%%esi),%%xmm6	\n\t"\
		"addpd	0x010(%%esi),%%xmm7	\n\t"\
		"subpd	%%xmm6,%%xmm0	\n\t"\
		"subpd	%%xmm7,%%xmm1	\n\t"\
		"movaps	%%xmm0,0x040(%%esi)	\n\t"\
		"subpd	%%xmm5,%%xmm2	\n\t"\
		"movaps	%%xmm1,0x050(%%esi)	\n\t"\
		"subpd	%%xmm4,%%xmm3	\n\t"\
		"addpd	%%xmm6,%%xmm6	\n\t"\
		"addpd	%%xmm5,%%xmm5	\n\t"\
		"addpd	%%xmm7,%%xmm7	\n\t"\
		"addpd	%%xmm4,%%xmm4	\n\t"\
		"addpd	%%xmm0,%%xmm6	\n\t"\
		"addpd	%%xmm2,%%xmm5	\n\t"\
		"addpd	%%xmm1,%%xmm7	\n\t"\
		"addpd	%%xmm3,%%xmm4	\n\t"\
		"movaps	%%xmm6,     (%%esi)	\n\t"\
		"movaps	%%xmm7,0x010(%%esi)	\n\t"\
		"movaps	0x380(%%esi),%%xmm0	/* isrt2 */	\n\t"\
		"movaps	%%xmm2,%%xmm6	\n\t"\
		"movaps	%%xmm5,%%xmm7	\n\t"\
		"subpd	%%xmm4,%%xmm2	\n\t"\
		"subpd	%%xmm3,%%xmm5	\n\t"\
		"addpd	%%xmm4,%%xmm6	\n\t"\
		"addpd	%%xmm3,%%xmm7	\n\t"\
		"mulpd	%%xmm0,%%xmm2	\n\t"\
		"mulpd	%%xmm0,%%xmm5	\n\t"\
		"mulpd	%%xmm0,%%xmm6	\n\t"\
		"mulpd	%%xmm0,%%xmm7	\n\t"\
		"movaps	%%xmm2,0x020(%%esi)	\n\t"\
		"movaps	%%xmm5,0x060(%%esi)	\n\t"\
		"movaps	%%xmm6,0x030(%%esi)	\n\t"\
		"movaps	%%xmm7,0x070(%%esi)	\n\t"\
		"/***************************************/	\n\t"\
		"	\n\t"\
		"/* SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r00) */	\n\t"\
		"subl	$0x80,%%esi	/* r00 */	\n\t"\
		"movaps	    (%%esi),%%xmm0	\n\t"\
		"movaps	0x40(%%esi),%%xmm4	\n\t"\
		"movaps	0x10(%%esi),%%xmm1	\n\t"\
		"movaps	0x50(%%esi),%%xmm5	\n\t"\
		"movaps	0x80(%%esi),%%xmm2	\n\t"\
		"movaps	0xd0(%%esi),%%xmm7	\n\t"\
		"movaps	0x90(%%esi),%%xmm3	\n\t"\
		"movaps	0xc0(%%esi),%%xmm6	\n\t"\
		"subpd   %%xmm2,%%xmm0	\n\t"\
		"subpd   %%xmm7,%%xmm4	\n\t"\
		"subpd   %%xmm3,%%xmm1	\n\t"\
		"subpd   %%xmm6,%%xmm5	\n\t"\
		"addpd   %%xmm2,%%xmm2	\n\t"\
		"addpd   %%xmm7,%%xmm7	\n\t"\
		"addpd   %%xmm3,%%xmm3	\n\t"\
		"addpd   %%xmm6,%%xmm6	\n\t"\
		"addpd   %%xmm0,%%xmm2	\n\t"\
		"addpd   %%xmm4,%%xmm7	\n\t"\
		"addpd   %%xmm1,%%xmm3	\n\t"\
		"addpd   %%xmm5,%%xmm6	\n\t"\
		"movaps	%%xmm0,0x80(%%esi)	\n\t"\
		"movaps	%%xmm4,0x40(%%esi)	\n\t"\
		"movaps	%%xmm1,0x90(%%esi)	\n\t"\
		"movaps	%%xmm5,0xd0(%%esi)	\n\t"\
		"movaps	%%xmm2,    (%%esi)	\n\t"\
		"movaps	%%xmm7,0xc0(%%esi)	\n\t"\
		"movaps	%%xmm3,0x10(%%esi)	\n\t"\
		"movaps	%%xmm6,0x50(%%esi)	\n\t"\
		"movaps	0x20(%%esi),%%xmm0	\n\t"\
		"movaps	0x60(%%esi),%%xmm4	\n\t"\
		"movaps	0x30(%%esi),%%xmm1	\n\t"\
		"movaps	0x70(%%esi),%%xmm5	\n\t"\
		"movaps	0xa0(%%esi),%%xmm2	\n\t"\
		"movaps	0xf0(%%esi),%%xmm7	\n\t"\
		"movaps	0xb0(%%esi),%%xmm3	\n\t"\
		"movaps	0xe0(%%esi),%%xmm6	\n\t"\
		"subpd   %%xmm2,%%xmm0	\n\t"\
		"subpd   %%xmm7,%%xmm4	\n\t"\
		"subpd   %%xmm3,%%xmm1	\n\t"\
		"subpd   %%xmm6,%%xmm5	\n\t"\
		"addpd   %%xmm2,%%xmm2	\n\t"\
		"addpd   %%xmm7,%%xmm7	\n\t"\
		"addpd   %%xmm3,%%xmm3	\n\t"\
		"addpd   %%xmm6,%%xmm6	\n\t"\
		"addpd   %%xmm0,%%xmm2	\n\t"\
		"addpd   %%xmm4,%%xmm7	\n\t"\
		"addpd   %%xmm1,%%xmm3	\n\t"\
		"addpd   %%xmm5,%%xmm6	\n\t"\
		"movaps	%%xmm0,0xa0(%%esi)	\n\t"\
		"movaps	%%xmm4,0x60(%%esi)	\n\t"\
		"movaps	%%xmm1,0xb0(%%esi)	\n\t"\
		"movaps	%%xmm5,0xf0(%%esi)	\n\t"\
		"movaps	%%xmm2,0x20(%%esi)	\n\t"\
		"movaps	%%xmm7,0xe0(%%esi)	\n\t"\
		"movaps	%%xmm3,0x30(%%esi)	\n\t"\
		"movaps	%%xmm6,0x70(%%esi)	\n\t"\
		"/***************************************/	\n\t"\
		"	\n\t"\
		"/*...Block 2: */	\n\t"\
		"subl	%%edi,%%eax	/* &a[j1] */	\n\t"\
		"subl	%%edi,%%ebx	\n\t"\
		"subl	%%edi,%%ecx	\n\t"\
		"subl	%%edi,%%edx	\n\t"\
		"movl	%[__p02],%%edi	/* Do things this way [rather than repeatedly adding p1] since array-padding scheme means p2 == p+p1 not guaranteed. */	\n\t"\
		"shll	$3,%%edi	\n\t"\
		"addl	%%edi,%%eax	/* &a[j1+p2] */	\n\t"\
		"addl	%%edi,%%ebx	\n\t"\
		"addl	%%edi,%%ecx	\n\t"\
		"addl	%%edi,%%edx	\n\t"\
		"/* SSE2_RADIX4_DIF_4TWIDDLE_B(r10,c02) */	\n\t"\
		"addl	$0x100,%%esi	/* r10 */	\n\t"\
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	    (%%ecx),%%xmm4	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"movaps	0x10(%%ecx),%%xmm5	\n\t"\
		"movaps	0x470(%%esi),%%xmm6	/* c02 */	\n\t"\
		"movaps	0x480(%%esi),%%xmm7	\n\t"\
		"movaps	%%xmm0,%%xmm2	\n\t"\
		"movaps	%%xmm1,%%xmm3	\n\t"\
		"mulpd   %%xmm6,%%xmm0	\n\t"\
		"mulpd   %%xmm6,%%xmm1	\n\t"\
		"mulpd   %%xmm7,%%xmm2	\n\t"\
		"mulpd   %%xmm7,%%xmm3	\n\t"\
		"movaps	%%xmm4,%%xmm6	\n\t"\
		"addpd   %%xmm2,%%xmm1	\n\t"\
		"movaps	%%xmm5,%%xmm7	\n\t"\
		"mulpd   0x490(%%esi),%%xmm4	/* c12 */	\n\t"\
		"subpd   %%xmm3,%%xmm0	\n\t"\
		"mulpd   0x490(%%esi),%%xmm5	\n\t"\
		"mulpd   0x4a0(%%esi),%%xmm6	\n\t"\
		"movaps	%%xmm0,%%xmm2	\n\t"\
		"mulpd	0x4a0(%%esi),%%xmm7	\n\t"\
		"addpd	%%xmm6,%%xmm5	\n\t"\
		"movaps	%%xmm1,%%xmm3	\n\t"\
		"subpd	%%xmm7,%%xmm4	\n\t"\
		"addpd	%%xmm4,%%xmm0	\n\t"\
		"addpd	%%xmm5,%%xmm1	\n\t"\
		"subpd	%%xmm4,%%xmm2	\n\t"\
		"subpd	%%xmm5,%%xmm3	\n\t"\
		"movaps	    (%%edx),%%xmm4	\n\t"\
		"movaps	0x10(%%edx),%%xmm5	\n\t"\
		"movaps	    (%%edx),%%xmm6	\n\t"\
		"movaps	0x10(%%edx),%%xmm7	\n\t"\
		"mulpd	0x4d0(%%esi),%%xmm4	/* c1A */	\n\t"\
		"mulpd	0x4d0(%%esi),%%xmm5	\n\t"\
		"mulpd	0x4e0(%%esi),%%xmm6	\n\t"\
		"mulpd	0x4e0(%%esi),%%xmm7	\n\t"\
		"addpd	%%xmm6,%%xmm5	\n\t"\
		"subpd	%%xmm7,%%xmm4	\n\t"\
		"movaps	%%xmm5,0x010(%%esi)	\n\t"\
		"movaps	%%xmm4,     (%%esi)	\n\t"\
		"movaps	    (%%ebx),%%xmm4	\n\t"\
		"movaps	0x10(%%ebx),%%xmm5	\n\t"\
		"movaps	    (%%ebx),%%xmm6	\n\t"\
		"movaps	0x10(%%ebx),%%xmm7	\n\t"\
		"mulpd	0x4b0(%%esi),%%xmm4	/* c0A */	\n\t"\
		"mulpd	0x4b0(%%esi),%%xmm5	\n\t"\
		"mulpd	0x4c0(%%esi),%%xmm6	\n\t"\
		"mulpd	0x4c0(%%esi),%%xmm7	\n\t"\
		"addpd	%%xmm6,%%xmm5	\n\t"\
		"subpd	%%xmm7,%%xmm4	\n\t"\
		"movaps	%%xmm5,%%xmm7	\n\t"\
		"movaps	%%xmm4,%%xmm6	\n\t"\
		"subpd	    (%%esi),%%xmm4	\n\t"\
		"subpd	0x10(%%esi),%%xmm5	\n\t"\
		"addpd	    (%%esi),%%xmm6	\n\t"\
		"addpd	0x10(%%esi),%%xmm7	\n\t"\
		"subpd	%%xmm6,%%xmm0	\n\t"\
		"subpd	%%xmm5,%%xmm2	\n\t"\
		"subpd	%%xmm7,%%xmm1	\n\t"\
		"subpd	%%xmm4,%%xmm3	\n\t"\
		"movaps	%%xmm0,0x040(%%esi)	\n\t"\
		"movaps	%%xmm2,0x020(%%esi)	\n\t"\
		"movaps	%%xmm1,0x050(%%esi)	\n\t"\
		"movaps	%%xmm3,0x070(%%esi)	\n\t"\
		"addpd	%%xmm6,%%xmm6	\n\t"\
		"addpd	%%xmm5,%%xmm5	\n\t"\
		"addpd	%%xmm7,%%xmm7	\n\t"\
		"addpd	%%xmm4,%%xmm4	\n\t"\
		"addpd	%%xmm0,%%xmm6	\n\t"\
		"addpd	%%xmm2,%%xmm5	\n\t"\
		"addpd	%%xmm1,%%xmm7	\n\t"\
		"addpd	%%xmm3,%%xmm4	\n\t"\
		"movaps	%%xmm6,     (%%esi)	\n\t"\
		"movaps	%%xmm5,0x060(%%esi)	\n\t"\
		"movaps	%%xmm7,0x010(%%esi)	\n\t"\
		"movaps	%%xmm4,0x030(%%esi)	\n\t"\
		"/***************************************/	\n\t"\
		"movl	%[__p04],%%edi	\n\t"\
		"shll	$3,%%edi	\n\t"\
		"addl	%%edi,%%eax	/* &a[j1+p6] */	\n\t"\
		"addl	%%edi,%%ebx	\n\t"\
		"addl	%%edi,%%ecx	\n\t"\
		"addl	%%edi,%%edx	\n\t"\
		"/* SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO_B(r18,c06) */	\n\t"\
		"/* esi contains r10 */	\n\t"\
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	    (%%ecx),%%xmm4	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"movaps	0x10(%%ecx),%%xmm5	\n\t"\
		"movaps	0x4f0(%%esi),%%xmm6	/* c06 */	\n\t"\
		"movaps	0x500(%%esi),%%xmm7	\n\t"\
		"movaps	%%xmm0,%%xmm2	\n\t"\
		"movaps	%%xmm1,%%xmm3	\n\t"\
		"mulpd	%%xmm6,%%xmm0	\n\t"\
		"mulpd	%%xmm6,%%xmm1	\n\t"\
		"mulpd	%%xmm7,%%xmm2	\n\t"\
		"mulpd	%%xmm7,%%xmm3	\n\t"\
		"movaps	%%xmm4,%%xmm6	\n\t"\
		"addpd	%%xmm2,%%xmm1	\n\t"\
		"movaps	%%xmm5,%%xmm7	\n\t"\
		"mulpd	0x510(%%esi),%%xmm4	/* c16 */	\n\t"\
		"subpd	%%xmm3,%%xmm0	\n\t"\
		"mulpd	0x510(%%esi),%%xmm5	\n\t"\
		"mulpd	0x520(%%esi),%%xmm6	\n\t"\
		"movaps	%%xmm0,%%xmm2	\n\t"\
		"mulpd	0x520(%%esi),%%xmm7	\n\t"\
		"addpd	%%xmm6,%%xmm5	\n\t"\
		"movaps	%%xmm1,%%xmm3	\n\t"\
		"subpd	%%xmm7,%%xmm4	\n\t"\
		"addpd	%%xmm4,%%xmm0	\n\t"\
		"addpd	%%xmm5,%%xmm1	\n\t"\
		"subpd	%%xmm4,%%xmm2	\n\t"\
		"subpd	%%xmm5,%%xmm3	\n\t"\
		"movaps	    (%%edx),%%xmm4	\n\t"\
		"movaps	0x10(%%edx),%%xmm5	\n\t"\
		"movaps	    (%%edx),%%xmm6	\n\t"\
		"movaps	0x10(%%edx),%%xmm7	\n\t"\
		"mulpd	0x550(%%esi),%%xmm4	/* c1E */	\n\t"\
		"mulpd	0x550(%%esi),%%xmm5	\n\t"\
		"mulpd	0x560(%%esi),%%xmm6	\n\t"\
		"mulpd	0x560(%%esi),%%xmm7	\n\t"\
		"addl	$0x80,%%esi	/* r18 */	\n\t"\
		"addpd	%%xmm6,%%xmm5	\n\t"\
		"subpd	%%xmm7,%%xmm4	\n\t"\
		"movaps	%%xmm5,0x010(%%esi)	\n\t"\
		"movaps	%%xmm4,     (%%esi)	\n\t"\
		"movaps	    (%%ebx),%%xmm4	\n\t"\
		"movaps	0x10(%%ebx),%%xmm5	\n\t"\
		"movaps	    (%%ebx),%%xmm6	\n\t"\
		"movaps	0x10(%%ebx),%%xmm7	\n\t"\
		"mulpd	0x4b0(%%esi),%%xmm4	/* c0E */	\n\t"\
		"mulpd	0x4b0(%%esi),%%xmm5	\n\t"\
		"mulpd	0x4c0(%%esi),%%xmm6	\n\t"\
		"mulpd	0x4c0(%%esi),%%xmm7	\n\t"\
		"addpd	%%xmm6,%%xmm5	\n\t"\
		"subpd	%%xmm7,%%xmm4	\n\t"\
		"movaps	%%xmm5,%%xmm7	\n\t"\
		"movaps	%%xmm4,%%xmm6	\n\t"\
		"subpd	     (%%esi),%%xmm4	\n\t"\
		"subpd	0x010(%%esi),%%xmm5	\n\t"\
		"addpd	     (%%esi),%%xmm6	\n\t"\
		"addpd	0x010(%%esi),%%xmm7	\n\t"\
		"subpd	%%xmm6,%%xmm0	\n\t"\
		"subpd	%%xmm7,%%xmm1	\n\t"\
		"movaps	%%xmm0,0x040(%%esi)	\n\t"\
		"subpd	%%xmm5,%%xmm2	\n\t"\
		"movaps	%%xmm1,0x050(%%esi)	\n\t"\
		"subpd	%%xmm4,%%xmm3	\n\t"\
		"addpd	%%xmm6,%%xmm6	\n\t"\
		"addpd	%%xmm5,%%xmm5	\n\t"\
		"addpd	%%xmm7,%%xmm7	\n\t"\
		"addpd	%%xmm4,%%xmm4	\n\t"\
		"addpd	%%xmm0,%%xmm6	\n\t"\
		"addpd	%%xmm2,%%xmm5	\n\t"\
		"addpd	%%xmm1,%%xmm7	\n\t"\
		"addpd	%%xmm3,%%xmm4	\n\t"\
		"movaps	%%xmm6,     (%%esi)	\n\t"\
		"movaps	%%xmm7,0x010(%%esi)	\n\t"\
		"movaps	0x280(%%esi),%%xmm0	/* isrt2 */	\n\t"\
		"movaps	%%xmm2,%%xmm6	\n\t"\
		"movaps	%%xmm5,%%xmm7	\n\t"\
		"subpd	%%xmm4,%%xmm2	\n\t"\
		"subpd	%%xmm3,%%xmm5	\n\t"\
		"addpd	%%xmm4,%%xmm6	\n\t"\
		"addpd	%%xmm3,%%xmm7	\n\t"\
		"mulpd	%%xmm0,%%xmm2	\n\t"\
		"mulpd	%%xmm0,%%xmm5	\n\t"\
		"mulpd	%%xmm0,%%xmm6	\n\t"\
		"mulpd	%%xmm0,%%xmm7	\n\t"\
		"movaps	%%xmm2,0x020(%%esi)	\n\t"\
		"movaps	%%xmm5,0x060(%%esi)	\n\t"\
		"movaps	%%xmm6,0x030(%%esi)	\n\t"\
		"movaps	%%xmm7,0x070(%%esi)	\n\t"\
		"/***************************************/	\n\t"\
		"	\n\t"\
		"/* SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r10) */	\n\t"\
		"subl		$0x80,%%esi	/* r10 */	\n\t"\
		"movaps		    (%%esi),%%xmm0	\n\t"\
		"movaps		0x40(%%esi),%%xmm4	\n\t"\
		"movaps		0x10(%%esi),%%xmm1	\n\t"\
		"movaps		0x50(%%esi),%%xmm5	\n\t"\
		"movaps		0x80(%%esi),%%xmm2	\n\t"\
		"movaps		0xd0(%%esi),%%xmm7	\n\t"\
		"movaps		0x90(%%esi),%%xmm3	\n\t"\
		"movaps		0xc0(%%esi),%%xmm6	\n\t"\
		"subpd   	%%xmm2,%%xmm0	\n\t"\
		"subpd   	%%xmm7,%%xmm4	\n\t"\
		"subpd   	%%xmm3,%%xmm1	\n\t"\
		"subpd   	%%xmm6,%%xmm5	\n\t"\
		"addpd   	%%xmm2,%%xmm2	\n\t"\
		"addpd   	%%xmm7,%%xmm7	\n\t"\
		"addpd   	%%xmm3,%%xmm3	\n\t"\
		"addpd   	%%xmm6,%%xmm6	\n\t"\
		"addpd   	%%xmm0,%%xmm2	\n\t"\
		"addpd   	%%xmm4,%%xmm7	\n\t"\
		"addpd   	%%xmm1,%%xmm3	\n\t"\
		"addpd   	%%xmm5,%%xmm6	\n\t"\
		"movaps		%%xmm0,0x80(%%esi)	\n\t"\
		"movaps		%%xmm4,0x40(%%esi)	\n\t"\
		"movaps		%%xmm1,0x90(%%esi)	\n\t"\
		"movaps		%%xmm5,0xd0(%%esi)	\n\t"\
		"movaps		%%xmm2,    (%%esi)	\n\t"\
		"movaps		%%xmm7,0xc0(%%esi)	\n\t"\
		"movaps		%%xmm3,0x10(%%esi)	\n\t"\
		"movaps		%%xmm6,0x50(%%esi)	\n\t"\
		"movaps		0x20(%%esi),%%xmm0	\n\t"\
		"movaps		0x60(%%esi),%%xmm4	\n\t"\
		"movaps		0x30(%%esi),%%xmm1	\n\t"\
		"movaps		0x70(%%esi),%%xmm5	\n\t"\
		"movaps		0xa0(%%esi),%%xmm2	\n\t"\
		"movaps		0xf0(%%esi),%%xmm7	\n\t"\
		"movaps		0xb0(%%esi),%%xmm3	\n\t"\
		"movaps		0xe0(%%esi),%%xmm6	\n\t"\
		"subpd   	%%xmm2,%%xmm0	\n\t"\
		"subpd   	%%xmm7,%%xmm4	\n\t"\
		"subpd   	%%xmm3,%%xmm1	\n\t"\
		"subpd   	%%xmm6,%%xmm5	\n\t"\
		"addpd   	%%xmm2,%%xmm2	\n\t"\
		"addpd   	%%xmm7,%%xmm7	\n\t"\
		"addpd   	%%xmm3,%%xmm3	\n\t"\
		"addpd   	%%xmm6,%%xmm6	\n\t"\
		"addpd   	%%xmm0,%%xmm2	\n\t"\
		"addpd   	%%xmm4,%%xmm7	\n\t"\
		"addpd   	%%xmm1,%%xmm3	\n\t"\
		"addpd   	%%xmm5,%%xmm6	\n\t"\
		"movaps		%%xmm0,0xa0(%%esi)	\n\t"\
		"movaps		%%xmm4,0x60(%%esi)	\n\t"\
		"movaps		%%xmm1,0xb0(%%esi)	\n\t"\
		"movaps		%%xmm5,0xf0(%%esi)	\n\t"\
		"movaps		%%xmm2,0x20(%%esi)	\n\t"\
		"movaps		%%xmm7,0xe0(%%esi)	\n\t"\
		"movaps		%%xmm3,0x30(%%esi)	\n\t"\
		"movaps		%%xmm6,0x70(%%esi)	\n\t"\
		"/***************************************/	\n\t"\
		"	\n\t"\
		"/*...Block 3: */	\n\t"\
		"movl	%[__add0],%%eax	\n\t"\
		"movl	%[__p01],%%edi	/* Do things this way [rather than repeatedly adding p1] since array-padding scheme means p2 == p+p1 not guaranteed. */	\n\t"\
		"movl	%[__p08],%%ebx	/* Can't get these via simple load-one-and-shift-as-needed due to array padding scheme */	\n\t"\
		"movl	%[__p10],%%ecx	\n\t"\
		"movl	%[__p18],%%edx	\n\t"\
		"shll	$3,%%edi	\n\t"\
		"shll	$3,%%ebx	\n\t"\
		"shll	$3,%%ecx	\n\t"\
		"shll	$3,%%edx	\n\t"\
		"addl	%%edi,%%eax	/* &a[j1+p1] */	\n\t"\
		"addl	%%eax,%%ebx	\n\t"\
		"addl	%%eax,%%ecx	\n\t"\
		"addl	%%eax,%%edx	\n\t"\
		"/* SSE2_RADIX4_DIF_4TWIDDLE_B(r20,c01) */	\n\t"\
		"addl	$0x100,%%esi	/* r20 */	\n\t"\
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	    (%%ecx),%%xmm4	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"movaps	0x10(%%ecx),%%xmm5	\n\t"\
		"movaps	0x470(%%esi),%%xmm6	/* c01 */	\n\t"\
		"movaps	0x480(%%esi),%%xmm7	\n\t"\
		"movaps	%%xmm0,%%xmm2	\n\t"\
		"movaps	%%xmm1,%%xmm3	\n\t"\
		"mulpd   %%xmm6,%%xmm0	\n\t"\
		"mulpd   %%xmm6,%%xmm1	\n\t"\
		"mulpd   %%xmm7,%%xmm2	\n\t"\
		"mulpd   %%xmm7,%%xmm3	\n\t"\
		"movaps	%%xmm4,%%xmm6	\n\t"\
		"addpd   %%xmm2,%%xmm1	\n\t"\
		"movaps	%%xmm5,%%xmm7	\n\t"\
		"mulpd   0x490(%%esi),%%xmm4	/* c11 */	\n\t"\
		"subpd   %%xmm3,%%xmm0	\n\t"\
		"mulpd   0x490(%%esi),%%xmm5	\n\t"\
		"mulpd   0x4a0(%%esi),%%xmm6	\n\t"\
		"movaps	%%xmm0,%%xmm2	\n\t"\
		"mulpd	0x4a0(%%esi),%%xmm7	\n\t"\
		"addpd	%%xmm6,%%xmm5	\n\t"\
		"movaps	%%xmm1,%%xmm3	\n\t"\
		"subpd	%%xmm7,%%xmm4	\n\t"\
		"addpd	%%xmm4,%%xmm0	\n\t"\
		"addpd	%%xmm5,%%xmm1	\n\t"\
		"subpd	%%xmm4,%%xmm2	\n\t"\
		"subpd	%%xmm5,%%xmm3	\n\t"\
		"movaps	    (%%edx),%%xmm4	\n\t"\
		"movaps	0x10(%%edx),%%xmm5	\n\t"\
		"movaps	    (%%edx),%%xmm6	\n\t"\
		"movaps	0x10(%%edx),%%xmm7	\n\t"\
		"mulpd	0x4d0(%%esi),%%xmm4	/* c19 */	\n\t"\
		"mulpd	0x4d0(%%esi),%%xmm5	\n\t"\
		"mulpd	0x4e0(%%esi),%%xmm6	\n\t"\
		"mulpd	0x4e0(%%esi),%%xmm7	\n\t"\
		"addpd	%%xmm6,%%xmm5	\n\t"\
		"subpd	%%xmm7,%%xmm4	\n\t"\
		"movaps	%%xmm5,0x010(%%esi)	\n\t"\
		"movaps	%%xmm4,     (%%esi)	\n\t"\
		"movaps	    (%%ebx),%%xmm4	\n\t"\
		"movaps	0x10(%%ebx),%%xmm5	\n\t"\
		"movaps	    (%%ebx),%%xmm6	\n\t"\
		"movaps	0x10(%%ebx),%%xmm7	\n\t"\
		"mulpd	0x4b0(%%esi),%%xmm4	/* c01 */	\n\t"\
		"mulpd	0x4b0(%%esi),%%xmm5	\n\t"\
		"mulpd	0x4c0(%%esi),%%xmm6	\n\t"\
		"mulpd	0x4c0(%%esi),%%xmm7	\n\t"\
		"addpd	%%xmm6,%%xmm5	\n\t"\
		"subpd	%%xmm7,%%xmm4	\n\t"\
		"movaps	%%xmm5,%%xmm7	\n\t"\
		"movaps	%%xmm4,%%xmm6	\n\t"\
		"subpd	    (%%esi),%%xmm4	\n\t"\
		"subpd	0x10(%%esi),%%xmm5	\n\t"\
		"addpd	    (%%esi),%%xmm6	\n\t"\
		"addpd	0x10(%%esi),%%xmm7	\n\t"\
		"subpd	%%xmm6,%%xmm0	\n\t"\
		"subpd	%%xmm5,%%xmm2	\n\t"\
		"subpd	%%xmm7,%%xmm1	\n\t"\
		"subpd	%%xmm4,%%xmm3	\n\t"\
		"movaps	%%xmm0,0x040(%%esi)	\n\t"\
		"movaps	%%xmm2,0x020(%%esi)	\n\t"\
		"movaps	%%xmm1,0x050(%%esi)	\n\t"\
		"movaps	%%xmm3,0x070(%%esi)	\n\t"\
		"addpd	%%xmm6,%%xmm6	\n\t"\
		"addpd	%%xmm5,%%xmm5	\n\t"\
		"addpd	%%xmm7,%%xmm7	\n\t"\
		"addpd	%%xmm4,%%xmm4	\n\t"\
		"addpd	%%xmm0,%%xmm6	\n\t"\
		"addpd	%%xmm2,%%xmm5	\n\t"\
		"addpd	%%xmm1,%%xmm7	\n\t"\
		"addpd	%%xmm3,%%xmm4	\n\t"\
		"movaps	%%xmm6,     (%%esi)	\n\t"\
		"movaps	%%xmm5,0x060(%%esi)	\n\t"\
		"movaps	%%xmm7,0x010(%%esi)	\n\t"\
		"movaps	%%xmm4,0x030(%%esi)	\n\t"\
		"/***************************************/	\n\t"\
		"movl	%[__p04],%%edi	\n\t"\
		"shll	$3,%%edi	\n\t"\
		"addl	%%edi,%%eax	/* &a[j1+p5] */	\n\t"\
		"addl	%%edi,%%ebx	\n\t"\
		"addl	%%edi,%%ecx	\n\t"\
		"addl	%%edi,%%edx	\n\t"\
		"/* SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO_B(r28,c05) */	\n\t"\
		"/* esi contains r20 */	\n\t"\
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	    (%%ecx),%%xmm4	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"movaps	0x10(%%ecx),%%xmm5	\n\t"\
		"movaps	0x4f0(%%esi),%%xmm6	/* c05 */	\n\t"\
		"movaps	0x500(%%esi),%%xmm7	\n\t"\
		"movaps	%%xmm0,%%xmm2	\n\t"\
		"movaps	%%xmm1,%%xmm3	\n\t"\
		"mulpd	%%xmm6,%%xmm0	\n\t"\
		"mulpd	%%xmm6,%%xmm1	\n\t"\
		"mulpd	%%xmm7,%%xmm2	\n\t"\
		"mulpd	%%xmm7,%%xmm3	\n\t"\
		"movaps	%%xmm4,%%xmm6	\n\t"\
		"addpd	%%xmm2,%%xmm1	\n\t"\
		"movaps	%%xmm5,%%xmm7	\n\t"\
		"mulpd	0x510(%%esi),%%xmm4	/* c15 */	\n\t"\
		"subpd	%%xmm3,%%xmm0	\n\t"\
		"mulpd	0x510(%%esi),%%xmm5	\n\t"\
		"mulpd	0x520(%%esi),%%xmm6	\n\t"\
		"movaps	%%xmm0,%%xmm2	\n\t"\
		"mulpd	0x520(%%esi),%%xmm7	\n\t"\
		"addpd	%%xmm6,%%xmm5	\n\t"\
		"movaps	%%xmm1,%%xmm3	\n\t"\
		"subpd	%%xmm7,%%xmm4	\n\t"\
		"addpd	%%xmm4,%%xmm0	\n\t"\
		"addpd	%%xmm5,%%xmm1	\n\t"\
		"subpd	%%xmm4,%%xmm2	\n\t"\
		"subpd	%%xmm5,%%xmm3	\n\t"\
		"movaps	    (%%edx),%%xmm4	\n\t"\
		"movaps	0x10(%%edx),%%xmm5	\n\t"\
		"movaps	    (%%edx),%%xmm6	\n\t"\
		"movaps	0x10(%%edx),%%xmm7	\n\t"\
		"mulpd	0x550(%%esi),%%xmm4	/* c1D */	\n\t"\
		"mulpd	0x550(%%esi),%%xmm5	\n\t"\
		"mulpd	0x560(%%esi),%%xmm6	\n\t"\
		"mulpd	0x560(%%esi),%%xmm7	\n\t"\
		"addl	$0x80,%%esi	/* r28 */	\n\t"\
		"addpd	%%xmm6,%%xmm5	\n\t"\
		"subpd	%%xmm7,%%xmm4	\n\t"\
		"movaps	%%xmm5,0x010(%%esi)	\n\t"\
		"movaps	%%xmm4,     (%%esi)	\n\t"\
		"movaps	    (%%ebx),%%xmm4	\n\t"\
		"movaps	0x10(%%ebx),%%xmm5	\n\t"\
		"movaps	    (%%ebx),%%xmm6	\n\t"\
		"movaps	0x10(%%ebx),%%xmm7	\n\t"\
		"mulpd	0x4b0(%%esi),%%xmm4	/* c0D */	\n\t"\
		"mulpd	0x4b0(%%esi),%%xmm5	\n\t"\
		"mulpd	0x4c0(%%esi),%%xmm6	\n\t"\
		"mulpd	0x4c0(%%esi),%%xmm7	\n\t"\
		"addpd	%%xmm6,%%xmm5	\n\t"\
		"subpd	%%xmm7,%%xmm4	\n\t"\
		"movaps	%%xmm5,%%xmm7	\n\t"\
		"movaps	%%xmm4,%%xmm6	\n\t"\
		"subpd	     (%%esi),%%xmm4	\n\t"\
		"subpd	0x010(%%esi),%%xmm5	\n\t"\
		"addpd	     (%%esi),%%xmm6	\n\t"\
		"addpd	0x010(%%esi),%%xmm7	\n\t"\
		"subpd	%%xmm6,%%xmm0	\n\t"\
		"subpd	%%xmm7,%%xmm1	\n\t"\
		"movaps	%%xmm0,0x040(%%esi)	\n\t"\
		"subpd	%%xmm5,%%xmm2	\n\t"\
		"movaps	%%xmm1,0x050(%%esi)	\n\t"\
		"subpd	%%xmm4,%%xmm3	\n\t"\
		"addpd	%%xmm6,%%xmm6	\n\t"\
		"addpd	%%xmm5,%%xmm5	\n\t"\
		"addpd	%%xmm7,%%xmm7	\n\t"\
		"addpd	%%xmm4,%%xmm4	\n\t"\
		"addpd	%%xmm0,%%xmm6	\n\t"\
		"addpd	%%xmm2,%%xmm5	\n\t"\
		"addpd	%%xmm1,%%xmm7	\n\t"\
		"addpd	%%xmm3,%%xmm4	\n\t"\
		"movaps	%%xmm6,     (%%esi)	\n\t"\
		"movaps	%%xmm7,0x010(%%esi)	\n\t"\
		"movaps	0x180(%%esi),%%xmm0	/* isrt2 */	\n\t"\
		"movaps	%%xmm2,%%xmm6	\n\t"\
		"movaps	%%xmm5,%%xmm7	\n\t"\
		"subpd	%%xmm4,%%xmm2	\n\t"\
		"subpd	%%xmm3,%%xmm5	\n\t"\
		"addpd	%%xmm4,%%xmm6	\n\t"\
		"addpd	%%xmm3,%%xmm7	\n\t"\
		"mulpd	%%xmm0,%%xmm2	\n\t"\
		"mulpd	%%xmm0,%%xmm5	\n\t"\
		"mulpd	%%xmm0,%%xmm6	\n\t"\
		"mulpd	%%xmm0,%%xmm7	\n\t"\
		"movaps	%%xmm2,0x020(%%esi)	\n\t"\
		"movaps	%%xmm5,0x060(%%esi)	\n\t"\
		"movaps	%%xmm6,0x030(%%esi)	\n\t"\
		"movaps	%%xmm7,0x070(%%esi)	\n\t"\
		"/***************************************/	\n\t"\
		"	\n\t"\
		"/* SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r20) */	\n\t"\
		"subl		$0x80,%%esi	/* r20 */	\n\t"\
		"movaps		    (%%esi),%%xmm0	\n\t"\
		"movaps		0x40(%%esi),%%xmm4	\n\t"\
		"movaps		0x10(%%esi),%%xmm1	\n\t"\
		"movaps		0x50(%%esi),%%xmm5	\n\t"\
		"movaps		0x80(%%esi),%%xmm2	\n\t"\
		"movaps		0xd0(%%esi),%%xmm7	\n\t"\
		"movaps		0x90(%%esi),%%xmm3	\n\t"\
		"movaps		0xc0(%%esi),%%xmm6	\n\t"\
		"subpd   	%%xmm2,%%xmm0	\n\t"\
		"subpd   	%%xmm7,%%xmm4	\n\t"\
		"subpd   	%%xmm3,%%xmm1	\n\t"\
		"subpd   	%%xmm6,%%xmm5	\n\t"\
		"addpd   	%%xmm2,%%xmm2	\n\t"\
		"addpd   	%%xmm7,%%xmm7	\n\t"\
		"addpd   	%%xmm3,%%xmm3	\n\t"\
		"addpd   	%%xmm6,%%xmm6	\n\t"\
		"addpd   	%%xmm0,%%xmm2	\n\t"\
		"addpd   	%%xmm4,%%xmm7	\n\t"\
		"addpd   	%%xmm1,%%xmm3	\n\t"\
		"addpd   	%%xmm5,%%xmm6	\n\t"\
		"movaps		%%xmm0,0x80(%%esi)	\n\t"\
		"movaps		%%xmm4,0x40(%%esi)	\n\t"\
		"movaps		%%xmm1,0x90(%%esi)	\n\t"\
		"movaps		%%xmm5,0xd0(%%esi)	\n\t"\
		"movaps		%%xmm2,    (%%esi)	\n\t"\
		"movaps		%%xmm7,0xc0(%%esi)	\n\t"\
		"movaps		%%xmm3,0x10(%%esi)	\n\t"\
		"movaps		%%xmm6,0x50(%%esi)	\n\t"\
		"movaps		0x20(%%esi),%%xmm0	\n\t"\
		"movaps		0x60(%%esi),%%xmm4	\n\t"\
		"movaps		0x30(%%esi),%%xmm1	\n\t"\
		"movaps		0x70(%%esi),%%xmm5	\n\t"\
		"movaps		0xa0(%%esi),%%xmm2	\n\t"\
		"movaps		0xf0(%%esi),%%xmm7	\n\t"\
		"movaps		0xb0(%%esi),%%xmm3	\n\t"\
		"movaps		0xe0(%%esi),%%xmm6	\n\t"\
		"subpd   	%%xmm2,%%xmm0	\n\t"\
		"subpd   	%%xmm7,%%xmm4	\n\t"\
		"subpd   	%%xmm3,%%xmm1	\n\t"\
		"subpd   	%%xmm6,%%xmm5	\n\t"\
		"addpd   	%%xmm2,%%xmm2	\n\t"\
		"addpd   	%%xmm7,%%xmm7	\n\t"\
		"addpd   	%%xmm3,%%xmm3	\n\t"\
		"addpd   	%%xmm6,%%xmm6	\n\t"\
		"addpd   	%%xmm0,%%xmm2	\n\t"\
		"addpd   	%%xmm4,%%xmm7	\n\t"\
		"addpd   	%%xmm1,%%xmm3	\n\t"\
		"addpd   	%%xmm5,%%xmm6	\n\t"\
		"movaps		%%xmm0,0xa0(%%esi)	\n\t"\
		"movaps		%%xmm4,0x60(%%esi)	\n\t"\
		"movaps		%%xmm1,0xb0(%%esi)	\n\t"\
		"movaps		%%xmm5,0xf0(%%esi)	\n\t"\
		"movaps		%%xmm2,0x20(%%esi)	\n\t"\
		"movaps		%%xmm7,0xe0(%%esi)	\n\t"\
		"movaps		%%xmm3,0x30(%%esi)	\n\t"\
		"movaps		%%xmm6,0x70(%%esi)	\n\t"\
		"/***************************************/	\n\t"\
		"	\n\t"\
		"/*...Block 4: */	\n\t"\
		"movl	%[__add0],%%eax	\n\t"\
		"movl	%[__p03],%%edi	/* Do things this way [rather than repeatedly adding p1] since array-padding scheme means p2 == p+p1 not guaranteed. */	\n\t"\
		"movl	%[__p08],%%ebx	/* Can't get these via simple load-one-and-shift-as-needed due to array padding scheme */	\n\t"\
		"movl	%[__p10],%%ecx	\n\t"\
		"movl	%[__p18],%%edx	\n\t"\
		"shll	$3,%%edi	\n\t"\
		"shll	$3,%%ebx	\n\t"\
		"shll	$3,%%ecx	\n\t"\
		"shll	$3,%%edx	\n\t"\
		"addl	%%edi,%%eax	/* &a[j1+p3] */	\n\t"\
		"addl	%%eax,%%ebx	\n\t"\
		"addl	%%eax,%%ecx	\n\t"\
		"addl	%%eax,%%edx	\n\t"\
		"/* SSE2_RADIX4_DIF_4TWIDDLE_B(r30,c03) */	\n\t"\
		"addl	$0x100,%%esi	/* r30 */	\n\t"\
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	    (%%ecx),%%xmm4	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"movaps	0x10(%%ecx),%%xmm5	\n\t"\
		"movaps	0x470(%%esi),%%xmm6	/* c03 */	\n\t"\
		"movaps	0x480(%%esi),%%xmm7	\n\t"\
		"movaps	%%xmm0,%%xmm2	\n\t"\
		"movaps	%%xmm1,%%xmm3	\n\t"\
		"mulpd   %%xmm6,%%xmm0	\n\t"\
		"mulpd   %%xmm6,%%xmm1	\n\t"\
		"mulpd   %%xmm7,%%xmm2	\n\t"\
		"mulpd   %%xmm7,%%xmm3	\n\t"\
		"movaps	%%xmm4,%%xmm6	\n\t"\
		"addpd   %%xmm2,%%xmm1	\n\t"\
		"movaps	%%xmm5,%%xmm7	\n\t"\
		"mulpd   0x490(%%esi),%%xmm4	/* c13 */	\n\t"\
		"subpd   %%xmm3,%%xmm0	\n\t"\
		"mulpd   0x490(%%esi),%%xmm5	\n\t"\
		"mulpd   0x4a0(%%esi),%%xmm6	\n\t"\
		"movaps	%%xmm0,%%xmm2	\n\t"\
		"mulpd	0x4a0(%%esi),%%xmm7	\n\t"\
		"addpd	%%xmm6,%%xmm5	\n\t"\
		"movaps	%%xmm1,%%xmm3	\n\t"\
		"subpd	%%xmm7,%%xmm4	\n\t"\
		"addpd	%%xmm4,%%xmm0	\n\t"\
		"addpd	%%xmm5,%%xmm1	\n\t"\
		"subpd	%%xmm4,%%xmm2	\n\t"\
		"subpd	%%xmm5,%%xmm3	\n\t"\
		"movaps	    (%%edx),%%xmm4	\n\t"\
		"movaps	0x10(%%edx),%%xmm5	\n\t"\
		"movaps	    (%%edx),%%xmm6	\n\t"\
		"movaps	0x10(%%edx),%%xmm7	\n\t"\
		"mulpd	0x4d0(%%esi),%%xmm4	/* c1B */	\n\t"\
		"mulpd	0x4d0(%%esi),%%xmm5	\n\t"\
		"mulpd	0x4e0(%%esi),%%xmm6	\n\t"\
		"mulpd	0x4e0(%%esi),%%xmm7	\n\t"\
		"addpd	%%xmm6,%%xmm5	\n\t"\
		"subpd	%%xmm7,%%xmm4	\n\t"\
		"movaps	%%xmm5,0x010(%%esi)	\n\t"\
		"movaps	%%xmm4,     (%%esi)	\n\t"\
		"movaps	    (%%ebx),%%xmm4	\n\t"\
		"movaps	0x10(%%ebx),%%xmm5	\n\t"\
		"movaps	    (%%ebx),%%xmm6	\n\t"\
		"movaps	0x10(%%ebx),%%xmm7	\n\t"\
		"mulpd	0x4b0(%%esi),%%xmm4	/* c0B */	\n\t"\
		"mulpd	0x4b0(%%esi),%%xmm5	\n\t"\
		"mulpd	0x4c0(%%esi),%%xmm6	\n\t"\
		"mulpd	0x4c0(%%esi),%%xmm7	\n\t"\
		"addpd	%%xmm6,%%xmm5	\n\t"\
		"subpd	%%xmm7,%%xmm4	\n\t"\
		"movaps	%%xmm5,%%xmm7	\n\t"\
		"movaps	%%xmm4,%%xmm6	\n\t"\
		"subpd	    (%%esi),%%xmm4	\n\t"\
		"subpd	0x10(%%esi),%%xmm5	\n\t"\
		"addpd	    (%%esi),%%xmm6	\n\t"\
		"addpd	0x10(%%esi),%%xmm7	\n\t"\
		"subpd	%%xmm6,%%xmm0	\n\t"\
		"subpd	%%xmm5,%%xmm2	\n\t"\
		"subpd	%%xmm7,%%xmm1	\n\t"\
		"subpd	%%xmm4,%%xmm3	\n\t"\
		"movaps	%%xmm0,0x040(%%esi)	\n\t"\
		"movaps	%%xmm2,0x020(%%esi)	\n\t"\
		"movaps	%%xmm1,0x050(%%esi)	\n\t"\
		"movaps	%%xmm3,0x070(%%esi)	\n\t"\
		"addpd	%%xmm6,%%xmm6	\n\t"\
		"addpd	%%xmm5,%%xmm5	\n\t"\
		"addpd	%%xmm7,%%xmm7	\n\t"\
		"addpd	%%xmm4,%%xmm4	\n\t"\
		"addpd	%%xmm0,%%xmm6	\n\t"\
		"addpd	%%xmm2,%%xmm5	\n\t"\
		"addpd	%%xmm1,%%xmm7	\n\t"\
		"addpd	%%xmm3,%%xmm4	\n\t"\
		"movaps	%%xmm6,     (%%esi)	\n\t"\
		"movaps	%%xmm5,0x060(%%esi)	\n\t"\
		"movaps	%%xmm7,0x010(%%esi)	\n\t"\
		"movaps	%%xmm4,0x030(%%esi)	\n\t"\
		"/***************************************/	\n\t"\
		"movl	%[__p04],%%edi	\n\t"\
		"shll	$3,%%edi	\n\t"\
		"addl	%%edi,%%eax	/* &a[j1+p7] */	\n\t"\
		"addl	%%edi,%%ebx	\n\t"\
		"addl	%%edi,%%ecx	\n\t"\
		"addl	%%edi,%%edx	\n\t"\
		"/* SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO_B(r38,c07) */	\n\t"\
		"/* esi contains r30 */	\n\t"\
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	    (%%ecx),%%xmm4	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"movaps	0x10(%%ecx),%%xmm5	\n\t"\
		"movaps	0x4f0(%%esi),%%xmm6	/* c07 */	\n\t"\
		"movaps	0x500(%%esi),%%xmm7	\n\t"\
		"movaps	%%xmm0,%%xmm2	\n\t"\
		"movaps	%%xmm1,%%xmm3	\n\t"\
		"mulpd	%%xmm6,%%xmm0	\n\t"\
		"mulpd	%%xmm6,%%xmm1	\n\t"\
		"mulpd	%%xmm7,%%xmm2	\n\t"\
		"mulpd	%%xmm7,%%xmm3	\n\t"\
		"movaps	%%xmm4,%%xmm6	\n\t"\
		"addpd	%%xmm2,%%xmm1	\n\t"\
		"movaps	%%xmm5,%%xmm7	\n\t"\
		"mulpd	0x510(%%esi),%%xmm4	/* c17 */	\n\t"\
		"subpd	%%xmm3,%%xmm0	\n\t"\
		"mulpd	0x510(%%esi),%%xmm5	\n\t"\
		"mulpd	0x520(%%esi),%%xmm6	\n\t"\
		"movaps	%%xmm0,%%xmm2	\n\t"\
		"mulpd	0x520(%%esi),%%xmm7	\n\t"\
		"addpd	%%xmm6,%%xmm5	\n\t"\
		"movaps	%%xmm1,%%xmm3	\n\t"\
		"subpd	%%xmm7,%%xmm4	\n\t"\
		"addpd	%%xmm4,%%xmm0	\n\t"\
		"addpd	%%xmm5,%%xmm1	\n\t"\
		"subpd	%%xmm4,%%xmm2	\n\t"\
		"subpd	%%xmm5,%%xmm3	\n\t"\
		"movaps	    (%%edx),%%xmm4	\n\t"\
		"movaps	0x10(%%edx),%%xmm5	\n\t"\
		"movaps	    (%%edx),%%xmm6	\n\t"\
		"movaps	0x10(%%edx),%%xmm7	\n\t"\
		"mulpd	0x550(%%esi),%%xmm4	/* c1F */	\n\t"\
		"mulpd	0x550(%%esi),%%xmm5	\n\t"\
		"mulpd	0x560(%%esi),%%xmm6	\n\t"\
		"mulpd	0x560(%%esi),%%xmm7	\n\t"\
		"addl	$0x80,%%esi	/* r38 */	\n\t"\
		"addpd	%%xmm6,%%xmm5	\n\t"\
		"subpd	%%xmm7,%%xmm4	\n\t"\
		"movaps	%%xmm5,0x010(%%esi)	\n\t"\
		"movaps	%%xmm4,     (%%esi)	\n\t"\
		"movaps	    (%%ebx),%%xmm4	\n\t"\
		"movaps	0x10(%%ebx),%%xmm5	\n\t"\
		"movaps	    (%%ebx),%%xmm6	\n\t"\
		"movaps	0x10(%%ebx),%%xmm7	\n\t"\
		"mulpd	0x4b0(%%esi),%%xmm4	/* c0F */	\n\t"\
		"mulpd	0x4b0(%%esi),%%xmm5	\n\t"\
		"mulpd	0x4c0(%%esi),%%xmm6	\n\t"\
		"mulpd	0x4c0(%%esi),%%xmm7	\n\t"\
		"addpd	%%xmm6,%%xmm5	\n\t"\
		"subpd	%%xmm7,%%xmm4	\n\t"\
		"movaps	%%xmm5,%%xmm7	\n\t"\
		"movaps	%%xmm4,%%xmm6	\n\t"\
		"subpd	     (%%esi),%%xmm4	\n\t"\
		"subpd	0x010(%%esi),%%xmm5	\n\t"\
		"addpd	     (%%esi),%%xmm6	\n\t"\
		"addpd	0x010(%%esi),%%xmm7	\n\t"\
		"subpd	%%xmm6,%%xmm0	\n\t"\
		"subpd	%%xmm7,%%xmm1	\n\t"\
		"movaps	%%xmm0,0x040(%%esi)	\n\t"\
		"subpd	%%xmm5,%%xmm2	\n\t"\
		"movaps	%%xmm1,0x050(%%esi)	\n\t"\
		"subpd	%%xmm4,%%xmm3	\n\t"\
		"addpd	%%xmm6,%%xmm6	\n\t"\
		"addpd	%%xmm5,%%xmm5	\n\t"\
		"addpd	%%xmm7,%%xmm7	\n\t"\
		"addpd	%%xmm4,%%xmm4	\n\t"\
		"addpd	%%xmm0,%%xmm6	\n\t"\
		"addpd	%%xmm2,%%xmm5	\n\t"\
		"addpd	%%xmm1,%%xmm7	\n\t"\
		"addpd	%%xmm3,%%xmm4	\n\t"\
		"movaps	%%xmm6,     (%%esi)	\n\t"\
		"movaps	%%xmm7,0x010(%%esi)	\n\t"\
		"movaps	0x080(%%esi),%%xmm0	/* isrt2 */	\n\t"\
		"movaps	%%xmm2,%%xmm6	\n\t"\
		"movaps	%%xmm5,%%xmm7	\n\t"\
		"subpd	%%xmm4,%%xmm2	\n\t"\
		"subpd	%%xmm3,%%xmm5	\n\t"\
		"addpd	%%xmm4,%%xmm6	\n\t"\
		"addpd	%%xmm3,%%xmm7	\n\t"\
		"mulpd	%%xmm0,%%xmm2	\n\t"\
		"mulpd	%%xmm0,%%xmm5	\n\t"\
		"mulpd	%%xmm0,%%xmm6	\n\t"\
		"mulpd	%%xmm0,%%xmm7	\n\t"\
		"movaps	%%xmm2,0x020(%%esi)	\n\t"\
		"movaps	%%xmm5,0x060(%%esi)	\n\t"\
		"movaps	%%xmm6,0x030(%%esi)	\n\t"\
		"movaps	%%xmm7,0x070(%%esi)	\n\t"\
		"/***************************************/	\n\t"\
		"	\n\t"\
		"/* SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r30) */	\n\t"\
		"subl		$0x80,%%esi	/* r30 */	\n\t"\
		"movaps		    (%%esi),%%xmm0	\n\t"\
		"movaps		0x40(%%esi),%%xmm4	\n\t"\
		"movaps		0x10(%%esi),%%xmm1	\n\t"\
		"movaps		0x50(%%esi),%%xmm5	\n\t"\
		"movaps		0x80(%%esi),%%xmm2	\n\t"\
		"movaps		0xd0(%%esi),%%xmm7	\n\t"\
		"movaps		0x90(%%esi),%%xmm3	\n\t"\
		"movaps		0xc0(%%esi),%%xmm6	\n\t"\
		"subpd   	%%xmm2,%%xmm0	\n\t"\
		"subpd   	%%xmm7,%%xmm4	\n\t"\
		"subpd   	%%xmm3,%%xmm1	\n\t"\
		"subpd   	%%xmm6,%%xmm5	\n\t"\
		"addpd   	%%xmm2,%%xmm2	\n\t"\
		"addpd   	%%xmm7,%%xmm7	\n\t"\
		"addpd   	%%xmm3,%%xmm3	\n\t"\
		"addpd   	%%xmm6,%%xmm6	\n\t"\
		"addpd   	%%xmm0,%%xmm2	\n\t"\
		"addpd   	%%xmm4,%%xmm7	\n\t"\
		"addpd   	%%xmm1,%%xmm3	\n\t"\
		"addpd   	%%xmm5,%%xmm6	\n\t"\
		"movaps		%%xmm0,0x80(%%esi)	\n\t"\
		"movaps		%%xmm4,0x40(%%esi)	\n\t"\
		"movaps		%%xmm1,0x90(%%esi)	\n\t"\
		"movaps		%%xmm5,0xd0(%%esi)	\n\t"\
		"movaps		%%xmm2,    (%%esi)	\n\t"\
		"movaps		%%xmm7,0xc0(%%esi)	\n\t"\
		"movaps		%%xmm3,0x10(%%esi)	\n\t"\
		"movaps		%%xmm6,0x50(%%esi)	\n\t"\
		"movaps		0x20(%%esi),%%xmm0	\n\t"\
		"movaps		0x60(%%esi),%%xmm4	\n\t"\
		"movaps		0x30(%%esi),%%xmm1	\n\t"\
		"movaps		0x70(%%esi),%%xmm5	\n\t"\
		"movaps		0xa0(%%esi),%%xmm2	\n\t"\
		"movaps		0xf0(%%esi),%%xmm7	\n\t"\
		"movaps		0xb0(%%esi),%%xmm3	\n\t"\
		"movaps		0xe0(%%esi),%%xmm6	\n\t"\
		"subpd   	%%xmm2,%%xmm0	\n\t"\
		"subpd   	%%xmm7,%%xmm4	\n\t"\
		"subpd   	%%xmm3,%%xmm1	\n\t"\
		"subpd   	%%xmm6,%%xmm5	\n\t"\
		"addpd   	%%xmm2,%%xmm2	\n\t"\
		"addpd   	%%xmm7,%%xmm7	\n\t"\
		"addpd   	%%xmm3,%%xmm3	\n\t"\
		"addpd   	%%xmm6,%%xmm6	\n\t"\
		"addpd   	%%xmm0,%%xmm2	\n\t"\
		"addpd   	%%xmm4,%%xmm7	\n\t"\
		"addpd   	%%xmm1,%%xmm3	\n\t"\
		"addpd   	%%xmm5,%%xmm6	\n\t"\
		"movaps		%%xmm0,0xa0(%%esi)	\n\t"\
		"movaps		%%xmm4,0x60(%%esi)	\n\t"\
		"movaps		%%xmm1,0xb0(%%esi)	\n\t"\
		"movaps		%%xmm5,0xf0(%%esi)	\n\t"\
		"movaps		%%xmm2,0x20(%%esi)	\n\t"\
		"movaps		%%xmm7,0xe0(%%esi)	\n\t"\
		"movaps		%%xmm3,0x30(%%esi)	\n\t"\
		"movaps		%%xmm6,0x70(%%esi)	\n\t"\
		"	\n\t"\
		"/**********************************************************************************/	\n\t"\
		"/*...and now do eight radix-4 transforms, including the internal twiddle factors: */	\n\t"\
		"/**********************************************************************************/	\n\t"\
		"/*...Block 1: t00,t10,t20,t30	*/	\n\t"\
		"movl	%[__add0],%%eax	/* &a[j1] */	\n\t"\
		"movl	%[__p01],%%ebx	\n\t"\
		"movl	%[__p02],%%ecx	\n\t"\
		"movl	%[__p03],%%edx	\n\t"\
		"shll	$3,%%ebx	\n\t"\
		"shll	$3,%%ecx	\n\t"\
		"shll	$3,%%edx	\n\t"\
		"addl	%%eax,%%ebx	\n\t"\
		"addl	%%eax,%%ecx	\n\t"\
		"addl	%%eax,%%edx	\n\t"\
		"movl	%[__r00],%%esi	\n\t"\
		"movaps	     (%%esi),%%xmm0	\n\t"\
		"movaps	0x200(%%esi),%%xmm4	\n\t"\
		"movaps	0x010(%%esi),%%xmm1	\n\t"\
		"movaps	0x210(%%esi),%%xmm5	\n\t"\
		"movaps	0x100(%%esi),%%xmm2	\n\t"\
		"movaps	0x300(%%esi),%%xmm6	\n\t"\
		"movaps	0x110(%%esi),%%xmm3	\n\t"\
		"movaps	0x310(%%esi),%%xmm7	\n\t"\
		"subpd	0x100(%%esi),%%xmm0	\n\t"\
		"subpd	0x300(%%esi),%%xmm4	\n\t"\
		"subpd	0x110(%%esi),%%xmm1	\n\t"\
		"subpd	0x310(%%esi),%%xmm5	\n\t"\
		"addpd	     (%%esi),%%xmm2	\n\t"\
		"addpd	0x200(%%esi),%%xmm6	\n\t"\
		"addpd	0x010(%%esi),%%xmm3	\n\t"\
		"addpd	0x210(%%esi),%%xmm7	\n\t"\
		"subpd	%%xmm6,%%xmm2	\n\t"\
		"subpd	%%xmm5,%%xmm0	\n\t"\
		"subpd	%%xmm7,%%xmm3	\n\t"\
		"subpd	%%xmm4,%%xmm1	\n\t"\
		"addpd	%%xmm6,%%xmm6	\n\t"\
		"addpd	%%xmm5,%%xmm5	\n\t"\
		"addpd	%%xmm7,%%xmm7	\n\t"\
		"addpd	%%xmm4,%%xmm4	\n\t"\
		"movaps	%%xmm2,    (%%ebx)	\n\t"\
		"movaps	%%xmm0,    (%%ecx)	\n\t"\
		"movaps	%%xmm3,0x10(%%ebx)	\n\t"\
		"movaps	%%xmm1,0x10(%%edx)	\n\t"\
		"addpd	%%xmm2,%%xmm6	\n\t"\
		"addpd	%%xmm0,%%xmm5	\n\t"\
		"addpd	%%xmm3,%%xmm7	\n\t"\
		"addpd	%%xmm1,%%xmm4	\n\t"\
		"movaps	%%xmm6,    (%%eax)	\n\t"\
		"movaps	%%xmm5,    (%%edx)	\n\t"\
		"movaps	%%xmm7,0x10(%%eax)	\n\t"\
		"movaps	%%xmm4,0x10(%%ecx)	\n\t"\
		"	\n\t"\
		"/*...Block 5: t08,t18,t28,t38	*/	\n\t"\
		"addl	$0x80,%%esi	/* r08 */	\n\t"\
		"subl	%%eax,%%ebx	/* p01 << 3 */	\n\t"\
		"subl	%%eax,%%ecx	/* p02 << 3 */	\n\t"\
		"subl	%%eax,%%edx	/* p03 << 3 */	\n\t"\
		"movl	%[__p04],%%edi	\n\t"\
		"shll	$3,%%edi	\n\t"\
		"addl	%%edi,%%eax	/* &a[j1+p04] */	\n\t"\
		"addl	%%eax,%%ebx	\n\t"\
		"addl	%%eax,%%ecx	\n\t"\
		"addl	%%eax,%%edx	\n\t"\
		"movaps	0x380(%%esi),%%xmm3	/* isrt2 */	\n\t"\
		"movaps	0x200(%%esi),%%xmm4	\n\t"\
		"movaps	0x210(%%esi),%%xmm5	\n\t"\
		"movaps	0x300(%%esi),%%xmm6	\n\t"\
		"movaps	0x310(%%esi),%%xmm7	\n\t"\
		"mulpd	%%xmm3,%%xmm4	\n\t"\
		"movaps	     (%%esi),%%xmm0	\n\t"\
		"mulpd	%%xmm3,%%xmm5	\n\t"\
		"movaps	0x010(%%esi),%%xmm1	\n\t"\
		"mulpd	%%xmm3,%%xmm6	\n\t"\
		"movaps	0x100(%%esi),%%xmm2	\n\t"\
		"mulpd	%%xmm3,%%xmm7	\n\t"\
		"movaps	0x110(%%esi),%%xmm3	\n\t"\
		"subpd	%%xmm3,%%xmm0	\n\t"\
		"subpd	%%xmm5,%%xmm4	\n\t"\
		"subpd	%%xmm2,%%xmm1	\n\t"\
		"subpd	%%xmm6,%%xmm7	\n\t"\
		"addpd	%%xmm3,%%xmm3	\n\t"\
		"addpd	%%xmm5,%%xmm5	\n\t"\
		"addpd	%%xmm2,%%xmm2	\n\t"\
		"addpd	%%xmm6,%%xmm6	\n\t"\
		"addpd	%%xmm0,%%xmm3	\n\t"\
		"addpd	%%xmm4,%%xmm5	\n\t"\
		"addpd	%%xmm1,%%xmm2	\n\t"\
		"addpd	%%xmm7,%%xmm6	\n\t"\
		"subpd	%%xmm6,%%xmm4	\n\t"\
		"subpd	%%xmm7,%%xmm5	\n\t"\
		"addpd	%%xmm6,%%xmm6	\n\t"\
		"addpd	%%xmm7,%%xmm7	\n\t"\
		"addpd	%%xmm4,%%xmm6	\n\t"\
		"addpd	%%xmm5,%%xmm7	\n\t"\
		"subpd	%%xmm4,%%xmm0	\n\t"\
		"subpd	%%xmm7,%%xmm3	\n\t"\
		"subpd	%%xmm5,%%xmm2	\n\t"\
		"subpd	%%xmm6,%%xmm1	\n\t"\
		"addpd	%%xmm4,%%xmm4	\n\t"\
		"addpd	%%xmm7,%%xmm7	\n\t"\
		"addpd	%%xmm5,%%xmm5	\n\t"\
		"addpd	%%xmm6,%%xmm6	\n\t"\
		"movaps	%%xmm0,    (%%ebx)	\n\t"\
		"movaps	%%xmm3,    (%%ecx)	\n\t"\
		"movaps	%%xmm2,0x10(%%ebx)	\n\t"\
		"movaps	%%xmm1,0x10(%%edx)	\n\t"\
		"addpd	%%xmm0,%%xmm4	\n\t"\
		"addpd	%%xmm3,%%xmm7	\n\t"\
		"addpd	%%xmm2,%%xmm5	\n\t"\
		"addpd	%%xmm1,%%xmm6	\n\t"\
		"movaps	%%xmm4,    (%%eax)	\n\t"\
		"movaps	%%xmm7,    (%%edx)	\n\t"\
		"movaps	%%xmm5,0x10(%%eax)	\n\t"\
		"movaps	%%xmm6,0x10(%%ecx)	\n\t"\
		"	\n\t"\
		"/*...Block 3: t04,t14,t24,t34	*/	\n\t"\
		"subl	$0x40,%%esi	/* r04 */	\n\t"\
		"subl	%%eax,%%ebx	/* p01 << 3 */	\n\t"\
		"subl	%%eax,%%ecx	/* p02 << 3 */	\n\t"\
		"subl	%%eax,%%edx	/* p03 << 3 */	\n\t"\
		"subl	%%edi,%%eax	/* &a[j1] */	\n\t"\
		"movl	%[__p08],%%edi	\n\t"\
		"shll	$3,%%edi	\n\t"\
		"addl	%%edi,%%eax	/* &a[j1+p08] */	\n\t"\
		"addl	%%eax,%%ebx	\n\t"\
		"addl	%%eax,%%ecx	\n\t"\
		"addl	%%eax,%%edx	\n\t"\
		"movaps	0x200(%%esi),%%xmm4	\n\t"\
		"movaps	0x210(%%esi),%%xmm5	\n\t"\
		"movaps	0x3d0(%%esi),%%xmm3	/* cc0 */	\n\t"\
		"movaps	0x3e0(%%esi),%%xmm2	\n\t"\
		"movaps	%%xmm4,%%xmm6	\n\t"\
		"movaps	%%xmm5,%%xmm7	\n\t"\
		"mulpd	%%xmm3,%%xmm4	\n\t"\
		"mulpd	%%xmm3,%%xmm5	\n\t"\
		"mulpd	%%xmm2,%%xmm6	\n\t"\
		"movaps	0x300(%%esi),%%xmm0	\n\t"\
		"mulpd	%%xmm2,%%xmm7	\n\t"\
		"movaps	0x310(%%esi),%%xmm1	\n\t"\
		"addpd	%%xmm6,%%xmm5	\n\t"\
		"movaps	%%xmm0,%%xmm6	\n\t"\
		"subpd	%%xmm7,%%xmm4	\n\t"\
		"movaps	%%xmm1,%%xmm7	\n\t"\
		"mulpd	%%xmm2,%%xmm6	\n\t"\
		"mulpd	%%xmm2,%%xmm7	\n\t"\
		"mulpd	%%xmm3,%%xmm0	\n\t"\
		"mulpd	%%xmm3,%%xmm1	\n\t"\
		"addpd	%%xmm0,%%xmm7	\n\t"\
		"subpd	%%xmm1,%%xmm6	\n\t"\
		"movaps	%%xmm4,%%xmm2	\n\t"\
		"movaps	%%xmm5,%%xmm3	\n\t"\
		"subpd	%%xmm6,%%xmm4	\n\t"\
		"subpd	%%xmm7,%%xmm5	\n\t"\
		"addpd	%%xmm2,%%xmm6	\n\t"\
		"addpd	%%xmm3,%%xmm7	\n\t"\
		"movaps	0x100(%%esi),%%xmm2	\n\t"\
		"movaps	0x110(%%esi),%%xmm3	\n\t"\
		"movaps	0x3c0(%%esi),%%xmm1	/* isrt2 */	\n\t"\
		"movaps	%%xmm2,%%xmm0	\n\t"\
		"subpd	%%xmm3,%%xmm2	\n\t"\
		"addpd	%%xmm0,%%xmm3	\n\t"\
		"mulpd	%%xmm1,%%xmm2	\n\t"\
		"mulpd	%%xmm1,%%xmm3	\n\t"\
		"movaps	     (%%esi),%%xmm0	\n\t"\
		"movaps	0x010(%%esi),%%xmm1	\n\t"\
		"subpd	%%xmm2,%%xmm0	\n\t"\
		"subpd	%%xmm3,%%xmm1	\n\t"\
		"addpd	%%xmm2,%%xmm2	\n\t"\
		"addpd	%%xmm3,%%xmm3	\n\t"\
		"addpd	%%xmm0,%%xmm2	\n\t"\
		"addpd	%%xmm1,%%xmm3	\n\t"\
		"subpd	%%xmm6,%%xmm2	\n\t"\
		"subpd	%%xmm5,%%xmm0	\n\t"\
		"subpd	%%xmm7,%%xmm3	\n\t"\
		"subpd	%%xmm4,%%xmm1	\n\t"\
		"addpd	%%xmm6,%%xmm6	\n\t"\
		"addpd	%%xmm5,%%xmm5	\n\t"\
		"addpd	%%xmm7,%%xmm7	\n\t"\
		"addpd	%%xmm4,%%xmm4	\n\t"\
		"movaps	%%xmm2,    (%%ebx)	\n\t"\
		"movaps	%%xmm0,    (%%ecx)	\n\t"\
		"movaps	%%xmm3,0x10(%%ebx)	\n\t"\
		"movaps	%%xmm1,0x10(%%edx)	\n\t"\
		"addpd	%%xmm2,%%xmm6	\n\t"\
		"addpd	%%xmm0,%%xmm5	\n\t"\
		"addpd	%%xmm3,%%xmm7	\n\t"\
		"addpd	%%xmm1,%%xmm4	\n\t"\
		"movaps	%%xmm6,    (%%eax)	\n\t"\
		"movaps	%%xmm5,    (%%edx)	\n\t"\
		"movaps	%%xmm7,0x10(%%eax)	\n\t"\
		"movaps	%%xmm4,0x10(%%ecx)	\n\t"\
		"	\n\t"\
		"/*...Block 7: t0C,t1C,t2C,t3C	*/	\n\t"\
		"addl	$0x80,%%esi	/* r0C */	\n\t"\
		"subl	%%eax,%%ebx	/* p01 << 3 */	\n\t"\
		"subl	%%eax,%%ecx	/* p02 << 3 */	\n\t"\
		"subl	%%eax,%%edx	/* p03 << 3 */	\n\t"\
		"subl	%%edi,%%eax	/* &a[j1] */	\n\t"\
		"movl	%[__p0C],%%edi	\n\t"\
		"shll	$3,%%edi	\n\t"\
		"addl	%%edi,%%eax	/* &a[j1+p0C] */	\n\t"\
		"addl	%%eax,%%ebx	\n\t"\
		"addl	%%eax,%%ecx	\n\t"\
		"addl	%%eax,%%edx	\n\t"\
		"movaps	0x200(%%esi),%%xmm4	\n\t"\
		"movaps	0x210(%%esi),%%xmm5	\n\t"\
		"movaps	0x350(%%esi),%%xmm2	/* cc0 */	\n\t"\
		"movaps	0x360(%%esi),%%xmm3	\n\t"\
		"movaps	%%xmm4,%%xmm6	\n\t"\
		"movaps	%%xmm5,%%xmm7	\n\t"\
		"mulpd	%%xmm3,%%xmm4	\n\t"\
		"mulpd	%%xmm3,%%xmm5	\n\t"\
		"mulpd	%%xmm2,%%xmm6	\n\t"\
		"movaps	0x300(%%esi),%%xmm0	\n\t"\
		"mulpd	%%xmm2,%%xmm7	\n\t"\
		"movaps	0x310(%%esi),%%xmm1	\n\t"\
		"addpd	%%xmm6,%%xmm5	\n\t"\
		"movaps	%%xmm0,%%xmm6	\n\t"\
		"subpd	%%xmm7,%%xmm4	\n\t"\
		"movaps	%%xmm1,%%xmm7	\n\t"\
		"mulpd	%%xmm2,%%xmm6	\n\t"\
		"mulpd	%%xmm2,%%xmm7	\n\t"\
		"mulpd	%%xmm3,%%xmm0	\n\t"\
		"mulpd	%%xmm3,%%xmm1	\n\t"\
		"addpd	%%xmm0,%%xmm7	\n\t"\
		"subpd	%%xmm1,%%xmm6	\n\t"\
		"movaps	%%xmm4,%%xmm2	\n\t"\
		"movaps	%%xmm5,%%xmm3	\n\t"\
		"subpd	%%xmm6,%%xmm4	\n\t"\
		"subpd	%%xmm7,%%xmm5	\n\t"\
		"addpd	%%xmm2,%%xmm6	\n\t"\
		"addpd	%%xmm3,%%xmm7	\n\t"\
		"movaps	0x100(%%esi),%%xmm2	\n\t"\
		"movaps	0x110(%%esi),%%xmm3	\n\t"\
		"movaps	0x340(%%esi),%%xmm1	/* isrt2 */	\n\t"\
		"movaps	%%xmm2,%%xmm0	\n\t"\
		"addpd	%%xmm3,%%xmm2	\n\t"\
		"subpd	%%xmm0,%%xmm3	\n\t"\
		"mulpd	%%xmm1,%%xmm2	\n\t"\
		"mulpd	%%xmm1,%%xmm3	\n\t"\
		"movaps	     (%%esi),%%xmm0	\n\t"\
		"movaps	0x010(%%esi),%%xmm1	\n\t"\
		"subpd	%%xmm2,%%xmm0	\n\t"\
		"subpd	%%xmm3,%%xmm1	\n\t"\
		"addpd	%%xmm2,%%xmm2	\n\t"\
		"addpd	%%xmm3,%%xmm3	\n\t"\
		"addpd	%%xmm0,%%xmm2	\n\t"\
		"addpd	%%xmm1,%%xmm3	\n\t"\
		"subpd	%%xmm4,%%xmm0	\n\t"\
		"subpd	%%xmm7,%%xmm2	\n\t"\
		"subpd	%%xmm5,%%xmm1	\n\t"\
		"subpd	%%xmm6,%%xmm3	\n\t"\
		"addpd	%%xmm4,%%xmm4	\n\t"\
		"addpd	%%xmm7,%%xmm7	\n\t"\
		"addpd	%%xmm5,%%xmm5	\n\t"\
		"addpd	%%xmm6,%%xmm6	\n\t"\
		"movaps	%%xmm0,    (%%ebx)	\n\t"\
		"movaps	%%xmm2,    (%%ecx)	\n\t"\
		"movaps	%%xmm1,0x10(%%ebx)	\n\t"\
		"movaps	%%xmm3,0x10(%%edx)	\n\t"\
		"addpd	%%xmm0,%%xmm4	\n\t"\
		"addpd	%%xmm2,%%xmm7	\n\t"\
		"addpd	%%xmm1,%%xmm5	\n\t"\
		"addpd	%%xmm3,%%xmm6	\n\t"\
		"movaps	%%xmm4,    (%%eax)	\n\t"\
		"movaps	%%xmm7,    (%%edx)	\n\t"\
		"movaps	%%xmm5,0x10(%%eax)	\n\t"\
		"movaps	%%xmm6,0x10(%%ecx)	\n\t"\
		"	\n\t"\
		"/*...Block 2: t02,t12,t22,t32	*/	\n\t"\
		"subl	$0xa0,%%esi	/* r02 */	\n\t"\
		"movl	%[__p10],%%edi	\n\t"\
		"subl	%%eax,%%ebx	\n\t"\
		"subl	%%eax,%%ecx	\n\t"\
		"subl	%%eax,%%edx	\n\t"\
		"movl	%[__add0],%%eax	/* &a[j1] */	\n\t"\
		"shll	$3,%%edi	\n\t"\
		"addl	%%edi,%%eax	/* &a[j1+p10) */	\n\t"\
		"addl	%%eax,%%ebx	\n\t"\
		"addl	%%eax,%%ecx	\n\t"\
		"addl	%%eax,%%edx	\n\t"\
		"movaps	0x200(%%esi),%%xmm4	\n\t"\
		"movaps	0x210(%%esi),%%xmm5	\n\t"\
		"movaps	0x410(%%esi),%%xmm2	/* cc1 */	\n\t"\
		"movaps	0x420(%%esi),%%xmm3	\n\t"\
		"movaps	%%xmm4,%%xmm6	\n\t"\
		"movaps	%%xmm5,%%xmm7	\n\t"\
		"mulpd	%%xmm2,%%xmm4	\n\t"\
		"mulpd	%%xmm2,%%xmm5	\n\t"\
		"mulpd	%%xmm3,%%xmm6	\n\t"\
		"movaps	0x300(%%esi),%%xmm0	\n\t"\
		"mulpd	%%xmm3,%%xmm7	\n\t"\
		"movaps	0x310(%%esi),%%xmm1	\n\t"\
		"movaps	0x430(%%esi),%%xmm2	/* cc3 */	\n\t"\
		"movaps	0x440(%%esi),%%xmm3	\n\t"\
		"addpd	%%xmm6,%%xmm5	\n\t"\
		"movaps	%%xmm0,%%xmm6	\n\t"\
		"subpd	%%xmm7,%%xmm4	\n\t"\
		"movaps	%%xmm1,%%xmm7	\n\t"\
		"mulpd	%%xmm2,%%xmm6	\n\t"\
		"mulpd	%%xmm2,%%xmm7	\n\t"\
		"mulpd	%%xmm3,%%xmm0	\n\t"\
		"mulpd	%%xmm3,%%xmm1	\n\t"\
		"addpd	%%xmm0,%%xmm7	\n\t"\
		"subpd	%%xmm1,%%xmm6	\n\t"\
		"movaps	%%xmm4,%%xmm2	\n\t"\
		"movaps	%%xmm5,%%xmm3	\n\t"\
		"subpd	%%xmm6,%%xmm4	\n\t"\
		"subpd	%%xmm7,%%xmm5	\n\t"\
		"addpd	%%xmm2,%%xmm6	\n\t"\
		"addpd	%%xmm3,%%xmm7	\n\t"\
		"movaps	0x100(%%esi),%%xmm1	\n\t"\
		"movaps	0x110(%%esi),%%xmm3	\n\t"\
		"movaps	0x400(%%esi),%%xmm0	/* ss1 */	\n\t"\
		"movaps	%%xmm1,%%xmm2	\n\t"\
		"mulpd	%%xmm0,%%xmm1	\n\t"\
		"mulpd	%%xmm3,%%xmm0	\n\t"\
		"mulpd	0x3f0(%%esi),%%xmm2		/* cc1 */	\n\t"\
		"mulpd	0x3f0(%%esi),%%xmm3	\n\t"\
		"subpd	%%xmm0,%%xmm2	\n\t"\
		"addpd	%%xmm1,%%xmm3	\n\t"\
		"movaps	     (%%esi),%%xmm0	\n\t"\
		"movaps	0x010(%%esi),%%xmm1	\n\t"\
		"subpd	%%xmm2,%%xmm0	\n\t"\
		"subpd	%%xmm3,%%xmm1	\n\t"\
		"addpd	%%xmm2,%%xmm2	\n\t"\
		"addpd	%%xmm3,%%xmm3	\n\t"\
		"addpd	%%xmm0,%%xmm2	\n\t"\
		"addpd	%%xmm1,%%xmm3	\n\t"\
		"subpd	%%xmm6,%%xmm2	\n\t"\
		"subpd	%%xmm5,%%xmm0	\n\t"\
		"subpd	%%xmm7,%%xmm3	\n\t"\
		"subpd	%%xmm4,%%xmm1	\n\t"\
		"addpd	%%xmm6,%%xmm6	\n\t"\
		"addpd	%%xmm5,%%xmm5	\n\t"\
		"addpd	%%xmm7,%%xmm7	\n\t"\
		"addpd	%%xmm4,%%xmm4	\n\t"\
		"movaps	%%xmm2,    (%%ebx)	\n\t"\
		"movaps	%%xmm0,    (%%ecx)	\n\t"\
		"movaps	%%xmm3,0x10(%%ebx)	\n\t"\
		"movaps	%%xmm1,0x10(%%edx)	\n\t"\
		"addpd	%%xmm2,%%xmm6	\n\t"\
		"addpd	%%xmm0,%%xmm5	\n\t"\
		"addpd	%%xmm3,%%xmm7	\n\t"\
		"addpd	%%xmm1,%%xmm4	\n\t"\
		"movaps	%%xmm6,    (%%eax)	\n\t"\
		"movaps	%%xmm5,    (%%edx)	\n\t"\
		"movaps	%%xmm7,0x10(%%eax)	\n\t"\
		"movaps	%%xmm4,0x10(%%ecx)	\n\t"\
		"	\n\t"\
		"/*...Block 6: t0A,t1A,t2A,t3A	*/	\n\t"\
		"addl	$0x80,%%esi	/* r0A */	\n\t"\
		"movl	%[__p04],%%edi	\n\t"\
		"shll	$3,%%edi	\n\t"\
		"addl	%%edi,%%eax	/* &a[j1+p14] */	\n\t"\
		"addl	%%edi,%%ebx	\n\t"\
		"addl	%%edi,%%ecx	\n\t"\
		"addl	%%edi,%%edx	\n\t"\
		"movaps	0x200(%%esi),%%xmm4	\n\t"\
		"movaps	0x210(%%esi),%%xmm5	\n\t"\
		"movaps	0x3b0(%%esi),%%xmm3	/* cc3 */	\n\t"\
		"movaps	0x3c0(%%esi),%%xmm2	\n\t"\
		"movaps	%%xmm4,%%xmm6	\n\t"\
		"movaps	%%xmm5,%%xmm7	\n\t"\
		"mulpd	%%xmm2,%%xmm4	\n\t"\
		"mulpd	%%xmm2,%%xmm5	\n\t"\
		"mulpd	%%xmm3,%%xmm6	\n\t"\
		"movaps	0x300(%%esi),%%xmm0	\n\t"\
		"mulpd	%%xmm3,%%xmm7	\n\t"\
		"movaps	0x310(%%esi),%%xmm1	\n\t"\
		"movaps	0x390(%%esi),%%xmm2	/* cc1 */	\n\t"\
		"movaps	0x3a0(%%esi),%%xmm3	\n\t"\
		"addpd	%%xmm6,%%xmm5	\n\t"\
		"movaps	%%xmm0,%%xmm6	\n\t"\
		"subpd	%%xmm7,%%xmm4	\n\t"\
		"movaps	%%xmm1,%%xmm7	\n\t"\
		"mulpd	%%xmm2,%%xmm6	\n\t"\
		"mulpd	%%xmm2,%%xmm7	\n\t"\
		"mulpd	%%xmm3,%%xmm0	\n\t"\
		"mulpd	%%xmm3,%%xmm1	\n\t"\
		"subpd	%%xmm0,%%xmm7	\n\t"\
		"addpd	%%xmm1,%%xmm6	\n\t"\
		"movaps	%%xmm4,%%xmm2	\n\t"\
		"movaps	%%xmm5,%%xmm3	\n\t"\
		"subpd	%%xmm6,%%xmm4	\n\t"\
		"subpd	%%xmm7,%%xmm5	\n\t"\
		"addpd	%%xmm2,%%xmm6	\n\t"\
		"addpd	%%xmm3,%%xmm7	\n\t"\
		"movaps	0x100(%%esi),%%xmm1	\n\t"\
		"movaps	0x110(%%esi),%%xmm3	\n\t"\
		"movaps	0x370(%%esi),%%xmm0		/* cc0 */	\n\t"\
		"movaps	%%xmm1,%%xmm2	\n\t"\
		"mulpd	%%xmm0,%%xmm1	\n\t"\
		"mulpd	%%xmm3,%%xmm0	\n\t"\
		"mulpd	0x380(%%esi),%%xmm2	/* ss0 */	\n\t"\
		"mulpd	0x380(%%esi),%%xmm3	\n\t"\
		"addpd	%%xmm0,%%xmm2	\n\t"\
		"subpd	%%xmm1,%%xmm3	\n\t"\
		"movaps	     (%%esi),%%xmm0	\n\t"\
		"movaps	0x010(%%esi),%%xmm1	\n\t"\
		"subpd	%%xmm2,%%xmm0	\n\t"\
		"subpd	%%xmm3,%%xmm1	\n\t"\
		"addpd	%%xmm2,%%xmm2	\n\t"\
		"addpd	%%xmm3,%%xmm3	\n\t"\
		"addpd	%%xmm0,%%xmm2	\n\t"\
		"addpd	%%xmm1,%%xmm3	\n\t"\
		"subpd	%%xmm4,%%xmm0	\n\t"\
		"subpd	%%xmm7,%%xmm2	\n\t"\
		"subpd	%%xmm5,%%xmm1	\n\t"\
		"subpd	%%xmm6,%%xmm3	\n\t"\
		"addpd	%%xmm4,%%xmm4	\n\t"\
		"addpd	%%xmm7,%%xmm7	\n\t"\
		"addpd	%%xmm5,%%xmm5	\n\t"\
		"addpd	%%xmm6,%%xmm6	\n\t"\
		"movaps	%%xmm0,    (%%ebx)	\n\t"\
		"movaps	%%xmm2,    (%%ecx)	\n\t"\
		"movaps	%%xmm1,0x10(%%ebx)	\n\t"\
		"movaps	%%xmm3,0x10(%%edx)	\n\t"\
		"addpd	%%xmm0,	%%xmm4	\n\t"\
		"addpd	%%xmm2,	%%xmm7	\n\t"\
		"addpd	%%xmm1,	%%xmm5	\n\t"\
		"addpd	%%xmm3,	%%xmm6	\n\t"\
		"movaps	%%xmm4,    (%%eax)	\n\t"\
		"movaps	%%xmm7,    (%%edx)	\n\t"\
		"movaps	%%xmm5,0x10(%%eax)	\n\t"\
		"movaps	%%xmm6,0x10(%%ecx)	\n\t"\
		"	\n\t"\
		"/*...Block 4: t06,t16,t26,t36	*/	\n\t"\
		"subl	$0x40,%%esi	/* r06 */	\n\t"\
		"movl	%[__p18],%%edi	\n\t"\
		"subl	%%eax,%%ebx	\n\t"\
		"subl	%%eax,%%ecx	\n\t"\
		"subl	%%eax,%%edx	\n\t"\
		"movl	%[__add0],%%eax	/* &a[j1] */	\n\t"\
		"shll	$3,%%edi	\n\t"\
		"addl	%%edi,%%eax	/* &a[j1+p18] */	\n\t"\
		"addl	%%eax,%%ebx	\n\t"\
		"addl	%%eax,%%ecx	\n\t"\
		"addl	%%eax,%%edx	\n\t"\
		"movaps	0x200(%%esi),%%xmm4	\n\t"\
		"movaps	0x210(%%esi),%%xmm5	\n\t"\
		"movaps	0x3f0(%%esi),%%xmm2	/* cc3 */	\n\t"\
		"movaps	0x400(%%esi),%%xmm3	\n\t"\
		"movaps	%%xmm4,%%xmm6	\n\t"\
		"movaps	%%xmm5,%%xmm7	\n\t"\
		"mulpd	%%xmm2,%%xmm4	\n\t"\
		"mulpd	%%xmm2,%%xmm5	\n\t"\
		"mulpd	%%xmm3,%%xmm6	\n\t"\
		"movaps	0x300(%%esi),%%xmm0	\n\t"\
		"mulpd	%%xmm3,%%xmm7	\n\t"\
		"movaps	0x310(%%esi),%%xmm1	\n\t"\
		"movaps	0x3d0(%%esi),%%xmm3	/* cc1 */	\n\t"\
		"movaps	0x3e0(%%esi),%%xmm2	\n\t"\
		"addpd	%%xmm6,%%xmm5	\n\t"\
		"movaps	%%xmm0,%%xmm6	\n\t"\
		"subpd	%%xmm7,%%xmm4	\n\t"\
		"movaps	%%xmm1,%%xmm7	\n\t"\
		"mulpd	%%xmm2,%%xmm6	\n\t"\
		"mulpd	%%xmm2,%%xmm7	\n\t"\
		"mulpd	%%xmm3,%%xmm0	\n\t"\
		"mulpd	%%xmm3,%%xmm1	\n\t"\
		"subpd	%%xmm0,%%xmm7	\n\t"\
		"addpd	%%xmm1,%%xmm6	\n\t"\
		"movaps	%%xmm4,%%xmm2	\n\t"\
		"movaps	%%xmm5,%%xmm3	\n\t"\
		"subpd	%%xmm6,%%xmm4	\n\t"\
		"subpd	%%xmm7,%%xmm5	\n\t"\
		"addpd	%%xmm2,%%xmm6	\n\t"\
		"addpd	%%xmm3,%%xmm7	\n\t"\
		"movaps	0x100(%%esi),%%xmm1	\n\t"\
		"movaps	0x110(%%esi),%%xmm3	\n\t"\
		"movaps	0x3b0(%%esi),%%xmm0		/* cc0 */	\n\t"\
		"movaps	%%xmm1,%%xmm2	\n\t"\
		"mulpd	%%xmm0,%%xmm1	\n\t"\
		"mulpd	%%xmm3,%%xmm0	\n\t"\
		"mulpd	0x3c0(%%esi),%%xmm2		/* ss0 */	\n\t"\
		"mulpd	0x3c0(%%esi),%%xmm3	\n\t"\
		"subpd	%%xmm0,%%xmm2	\n\t"\
		"addpd	%%xmm1,%%xmm3	\n\t"\
		"movaps	     (%%esi),%%xmm0	\n\t"\
		"movaps	0x010(%%esi),%%xmm1	\n\t"\
		"subpd	%%xmm2,%%xmm0	\n\t"\
		"subpd	%%xmm3,%%xmm1	\n\t"\
		"addpd	%%xmm2,%%xmm2	\n\t"\
		"addpd	%%xmm3,%%xmm3	\n\t"\
		"addpd	%%xmm0,%%xmm2	\n\t"\
		"addpd	%%xmm1,%%xmm3	\n\t"\
		"subpd	%%xmm4,%%xmm2	\n\t"\
		"subpd	%%xmm7,%%xmm0	\n\t"\
		"subpd	%%xmm5,%%xmm3	\n\t"\
		"subpd	%%xmm6,%%xmm1	\n\t"\
		"addpd	%%xmm4,%%xmm4	\n\t"\
		"addpd	%%xmm7,%%xmm7	\n\t"\
		"addpd	%%xmm5,%%xmm5	\n\t"\
		"addpd	%%xmm6,%%xmm6	\n\t"\
		"movaps	%%xmm2,    (%%ebx)	\n\t"\
		"movaps	%%xmm0,    (%%ecx)	\n\t"\
		"movaps	%%xmm3,0x10(%%ebx)	\n\t"\
		"movaps	%%xmm1,0x10(%%edx)	\n\t"\
		"addpd	%%xmm2,	%%xmm4	\n\t"\
		"addpd	%%xmm0,	%%xmm7	\n\t"\
		"addpd	%%xmm3,	%%xmm5	\n\t"\
		"addpd	%%xmm1,	%%xmm6	\n\t"\
		"movaps	%%xmm4,    (%%eax)	\n\t"\
		"movaps	%%xmm7,    (%%edx)	\n\t"\
		"movaps	%%xmm5,0x10(%%eax)	\n\t"\
		"movaps	%%xmm6,0x10(%%ecx)	\n\t"\
		"	\n\t"\
		"/*...Block 8: t0E,t1E,t2E,t3E	*/	\n\t"\
		"addl	$0x80,%%esi	/* r0E */	\n\t"\
		"subl	%%eax,%%ebx	/* p01 << 3 */	\n\t"\
		"subl	%%eax,%%ecx	/* p02 << 3 */	\n\t"\
		"subl	%%eax,%%edx	/* p03 << 3 */	\n\t"\
		"movl	%[__p04],%%edi	\n\t"\
		"shll	$3,%%edi	\n\t"\
		"addl	%%edi,%%eax	/* &a[j1+p1C] */	\n\t"\
		"addl	%%eax,%%ebx	\n\t"\
		"addl	%%eax,%%ecx	\n\t"\
		"addl	%%eax,%%edx	\n\t"\
		"movaps	0x200(%%esi),%%xmm4	\n\t"\
		"movaps	0x210(%%esi),%%xmm5	\n\t"\
		"movaps	0x350(%%esi),%%xmm3	/* cc1 */	\n\t"\
		"movaps	0x360(%%esi),%%xmm2	\n\t"\
		"movaps	%%xmm4,%%xmm6	\n\t"\
		"movaps	%%xmm5,%%xmm7	\n\t"\
		"mulpd	%%xmm2,%%xmm4	\n\t"\
		"mulpd	%%xmm2,%%xmm5	\n\t"\
		"mulpd	%%xmm3,%%xmm6	\n\t"\
		"movaps	0x300(%%esi),%%xmm0	\n\t"\
		"mulpd	%%xmm3,%%xmm7	\n\t"\
		"movaps	0x310(%%esi),%%xmm1	\n\t"\
		"movaps	0x370(%%esi),%%xmm3	/* cc3 */	\n\t"\
		"movaps	0x380(%%esi),%%xmm2	\n\t"\
		"addpd	%%xmm6,%%xmm5	\n\t"\
		"movaps	%%xmm0,%%xmm6	\n\t"\
		"subpd	%%xmm7,%%xmm4	\n\t"\
		"movaps	%%xmm1,%%xmm7	\n\t"\
		"mulpd	%%xmm2,%%xmm6	\n\t"\
		"mulpd	%%xmm2,%%xmm7	\n\t"\
		"mulpd	%%xmm3,%%xmm0	\n\t"\
		"mulpd	%%xmm3,%%xmm1	\n\t"\
		"addpd	%%xmm0,%%xmm7	\n\t"\
		"subpd	%%xmm1,%%xmm6	\n\t"\
		"movaps	%%xmm4,%%xmm2	\n\t"\
		"movaps	%%xmm5,%%xmm3	\n\t"\
		"subpd	%%xmm6,%%xmm4	\n\t"\
		"subpd	%%xmm7,%%xmm5	\n\t"\
		"addpd	%%xmm2,%%xmm6	\n\t"\
		"addpd	%%xmm3,%%xmm7	\n\t"\
		"movaps	0x100(%%esi),%%xmm1	\n\t"\
		"movaps	0x110(%%esi),%%xmm3	\n\t"\
		"movaps	0x340(%%esi),%%xmm0	/* ss0 */	\n\t"\
		"movaps	%%xmm1,%%xmm2	\n\t"\
		"mulpd	%%xmm0,%%xmm1	\n\t"\
		"mulpd	%%xmm3,%%xmm0	\n\t"\
		"mulpd	0x330(%%esi),%%xmm2	/* cc0 */	\n\t"\
		"mulpd	0x330(%%esi),%%xmm3	\n\t"\
		"addpd	%%xmm0,%%xmm2	\n\t"\
		"subpd	%%xmm1,%%xmm3	\n\t"\
		"movaps	     (%%esi),%%xmm0	\n\t"\
		"movaps	0x010(%%esi),%%xmm1	\n\t"\
		"subpd	%%xmm2,%%xmm0	\n\t"\
		"subpd	%%xmm3,%%xmm1	\n\t"\
		"addpd	%%xmm2,%%xmm2	\n\t"\
		"addpd	%%xmm3,%%xmm3	\n\t"\
		"addpd	%%xmm0,%%xmm2	\n\t"\
		"addpd	%%xmm1,%%xmm3	\n\t"\
		"subpd	%%xmm4,%%xmm0	\n\t"\
		"subpd	%%xmm7,%%xmm2	\n\t"\
		"subpd	%%xmm5,%%xmm1	\n\t"\
		"subpd	%%xmm6,%%xmm3	\n\t"\
		"addpd	%%xmm4,%%xmm4	\n\t"\
		"addpd	%%xmm7,%%xmm7	\n\t"\
		"addpd	%%xmm5,%%xmm5	\n\t"\
		"addpd	%%xmm6,%%xmm6	\n\t"\
		"movaps	%%xmm0,    (%%ebx)	\n\t"\
		"movaps	%%xmm2,    (%%ecx)	\n\t"\
		"movaps	%%xmm1,0x10(%%ebx)	\n\t"\
		"movaps	%%xmm3,0x10(%%edx)	\n\t"\
		"addpd	%%xmm0,	%%xmm4	\n\t"\
		"addpd	%%xmm2,	%%xmm7	\n\t"\
		"addpd	%%xmm1,	%%xmm5	\n\t"\
		"addpd	%%xmm3,	%%xmm6	\n\t"\
		"movaps	%%xmm4,    (%%eax)	\n\t"\
		"movaps	%%xmm7,    (%%edx)	\n\t"\
		"movaps	%%xmm5,0x10(%%eax)	\n\t"\
		"movaps	%%xmm6,0x10(%%ecx)	\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__p01] "m" (Xp01)\
		 ,[__p02] "m" (Xp02)\
		 ,[__p03] "m" (Xp03)\
		 ,[__p04] "m" (Xp04)\
		 ,[__p08] "m" (Xp08)\
		 ,[__p0C] "m" (Xp0C)\
		 ,[__p10] "m" (Xp10)\
		 ,[__p18] "m" (Xp18)\
		 ,[__r00] "m" (Xr00)\
		: "eax","ebx","ecx","edx","edi","esi"		/* Clobbered registers */\
	);\
	}


	/*
	For GCC-macro version of this, use that isrt2 + 0x010,0x030,0x050 = cc0,cc1,cc3,
	and isrt2 + 0x070,0x0f0,0x170,0x1f0,0x270,0x2f0,0x370,0x3f0 = c00,04,02,06,01,05,03,07
	in order to reduce number of args to <= the GCC-allowed maximum of 30:
	*/
	#define SSE2_RADIX32_DIT_TWIDDLE(Xadd0,Xp01,Xp02,Xp03,Xp04,Xp05,Xp06,Xp07,Xp08,Xp10,Xp18,Xr00,Xr02,Xr04,Xr06,Xr08,Xr0A,Xr0C,Xr0E,Xr10,Xr12,Xr14,Xr16,Xr18,Xr1A,Xr1C,Xr1E,Xr20,Xr30,Xisrt2)\
	{\
	__asm__ volatile (\
		"movl	%[__add0],%%eax\n\t"\
		"movl	%[__p05],%%ebx\n\t"\
		"movl	%[__p06],%%ecx\n\t"\
		"movl	%[__p07],%%edx\n\t"\
		"movl	%[__p04],%%edi	/* edi will store copy of p4 throughout */\n\t"\
		"shll	$3,%%ebx\n\t"\
		"shll	$3,%%ecx\n\t"\
		"shll	$3,%%edx\n\t"\
		"shll	$3,%%edi\n\t"\
		"addl	%%eax,%%ebx\n\t"\
		"addl	%%eax,%%ecx\n\t"\
		"addl	%%eax,%%edx\n\t"\
		"addl	%%edi,%%eax\n\t"\
		"/* SSE2_RADIX8_DIT_0TWIDDLE(r00): */\n\t"\
		"movl	%[__r00]	,%%esi 		\n\t"\
		"movaps	    (%%eax)	,%%xmm0		\n\t"\
		"movaps	0x10(%%eax)	,%%xmm1		\n\t"\
		"movaps	%%xmm0		,%%xmm2		\n\t"\
		"movaps	%%xmm1		,%%xmm3		\n\t"\
		"addpd	    (%%ebx)	,%%xmm2		\n\t"\
		"addpd	0x10(%%ebx)	,%%xmm3		\n\t"\
		"subpd	    (%%ebx)	,%%xmm0		\n\t"\
		"subpd	0x10(%%ebx)	,%%xmm1		\n\t"\
		"movaps	    (%%ecx)	,%%xmm4		\n\t"\
		"movaps	0x10(%%ecx)	,%%xmm5		\n\t"\
		"movaps	%%xmm4		,%%xmm6		\n\t"\
		"movaps	%%xmm5		,%%xmm7		\n\t"\
		"addpd	    (%%edx)	,%%xmm6		\n\t"\
		"addpd	0x10(%%edx)	,%%xmm7		\n\t"\
		"subpd	    (%%edx)	,%%xmm4		\n\t"\
		"subpd	0x10(%%edx)	,%%xmm5		\n\t"\
		"movaps	%%xmm6		,0xc0(%%esi)\n\t"\
		"movaps	%%xmm7		,0xd0(%%esi)\n\t"\
		"addpd	%%xmm2		,%%xmm6		\n\t"\
		"addpd	%%xmm3		,%%xmm7		\n\t"\
		"subpd	0xc0(%%esi)	,%%xmm2		\n\t"\
		"subpd	0xd0(%%esi)	,%%xmm3		\n\t"\
		"movaps	%%xmm2		,0xc0(%%esi)\n\t"\
		"movaps	%%xmm3		,0xd0(%%esi)\n\t"\
		"movaps	%%xmm4		,%%xmm2		\n\t"\
		"movaps	%%xmm5		,%%xmm3		\n\t"\
		"addpd	%%xmm0		,%%xmm5		\n\t"\
		"subpd	%%xmm3		,%%xmm0		\n\t"\
		"addpd	%%xmm1		,%%xmm4		\n\t"\
		"subpd	%%xmm2		,%%xmm1		\n\t"\
		"movl	%[__isrt2]	,%%esi		\n\t"\
		"movaps	%%xmm5		,%%xmm2		\n\t"\
		"movaps	%%xmm1		,%%xmm3		\n\t"\
		"addpd	%%xmm1		,%%xmm5		\n\t"\
		"movaps	(%%esi)		,%%xmm1		\n\t"\
		"subpd	%%xmm3		,%%xmm2		\n\t"\
		"movl	%[__r00]	,%%esi		\n\t"\
		"mulpd	%%xmm1		,%%xmm5		\n\t"\
		"mulpd	%%xmm1		,%%xmm2		\n\t"\
		"movaps	%%xmm0		,%%xmm3		\n\t"\
		"movaps	%%xmm5		,0xa0(%%esi)\n\t"\
		"movaps	%%xmm4		,%%xmm5		\n\t"\
		"addpd	%%xmm4		,%%xmm0		\n\t"\
		"movaps	%%xmm2		,0xb0(%%esi)\n\t"\
		"subpd	%%xmm5		,%%xmm3		\n\t"\
		"mulpd	%%xmm1		,%%xmm0		\n\t"\
		"mulpd	%%xmm1		,%%xmm3		\n\t"\
		"movaps	%%xmm0		,0xe0(%%esi)\n\t"\
		"movaps	%%xmm3		,0xf0(%%esi)\n\t"\
		"subl	%%edi		,%%eax		\n\t"\
		"subl	%%edi		,%%ebx		\n\t"\
		"subl	%%edi		,%%ecx		\n\t"\
		"subl	%%edi		,%%edx		\n\t"\
		"movaps	    (%%eax)	,%%xmm0		\n\t"\
		"movaps	0x10(%%eax)	,%%xmm1		\n\t"\
		"movaps	%%xmm0		,%%xmm2		\n\t"\
		"movaps	%%xmm1		,%%xmm3		\n\t"\
		"addpd	    (%%ebx)	,%%xmm2		\n\t"\
		"addpd	0x10(%%ebx)	,%%xmm3		\n\t"\
		"subpd	    (%%ebx)	,%%xmm0		\n\t"\
		"subpd	0x10(%%ebx)	,%%xmm1		\n\t"\
		"movaps	%%xmm6		,    (%%esi)\n\t"\
		"movaps	%%xmm7		,0x10(%%esi)\n\t"\
		"addpd	%%xmm6		,%%xmm6		\n\t"\
		"addpd	%%xmm7		,%%xmm7		\n\t"\
		"movaps	%%xmm6		,0x80(%%esi)\n\t"\
		"movaps	%%xmm7		,0x90(%%esi)\n\t"\
		"movaps	    (%%ecx)	,%%xmm4		\n\t"\
		"movaps	0x10(%%ecx)	,%%xmm5		\n\t"\
		"movaps	%%xmm4		,%%xmm6		\n\t"\
		"movaps	%%xmm5		,%%xmm7		\n\t"\
		"addpd	    (%%edx)	,%%xmm6		\n\t"\
		"addpd	0x10(%%edx)	,%%xmm7		\n\t"\
		"subpd	    (%%edx)	,%%xmm4		\n\t"\
		"subpd	0x10(%%edx)	,%%xmm5		\n\t"\
		"movaps	%%xmm6		,0x40(%%esi)\n\t"\
		"movaps	%%xmm7		,0x50(%%esi)\n\t"\
		"addpd	%%xmm2		,%%xmm6		\n\t"\
		"addpd	%%xmm3		,%%xmm7		\n\t"\
		"subpd	0x40(%%esi)	,%%xmm2		\n\t"\
		"subpd	0x50(%%esi)	,%%xmm3		\n\t"\
		"addpd	    (%%esi)	,%%xmm6		\n\t"\
		"addpd	0x10(%%esi)	,%%xmm7		\n\t"\
		"movaps	%%xmm6		,    (%%esi)\n\t"\
		"movaps	%%xmm7		,0x10(%%esi)\n\t"\
		"subpd	0x80(%%esi)	,%%xmm6		\n\t"\
		"subpd	0x90(%%esi)	,%%xmm7		\n\t"\
		"movaps	%%xmm6		,0x80(%%esi)\n\t"\
		"movaps	%%xmm7		,0x90(%%esi)\n\t"\
		"movaps	%%xmm4		,%%xmm6		\n\t"\
		"movaps	%%xmm5		,%%xmm7		\n\t"\
		"addpd	%%xmm0		,%%xmm5		\n\t"\
		"subpd	%%xmm7		,%%xmm0		\n\t"\
		"addpd	%%xmm1		,%%xmm4		\n\t"\
		"subpd	%%xmm6		,%%xmm1		\n\t"\
		"movaps	%%xmm5		,%%xmm6		\n\t"\
		"movaps	%%xmm1		,%%xmm7		\n\t"\
		"addpd	0xa0(%%esi)	,%%xmm5		\n\t"\
		"subpd	0xb0(%%esi)	,%%xmm1		\n\t"\
		"subpd	0xa0(%%esi)	,%%xmm6		\n\t"\
		"addpd	0xb0(%%esi)	,%%xmm7		\n\t"\
		"movaps	%%xmm5		,0x20(%%esi)\n\t"\
		"movaps	%%xmm1		,0x30(%%esi)\n\t"\
		"movaps	%%xmm6		,0xa0(%%esi)\n\t"\
		"movaps	%%xmm7		,0xb0(%%esi)\n\t"\
		"movaps	%%xmm2		,%%xmm6		\n\t"\
		"movaps	%%xmm3		,%%xmm7		\n\t"\
		"addpd	0xd0(%%esi)	,%%xmm2		\n\t"\
		"subpd	0xc0(%%esi)	,%%xmm3		\n\t"\
		"subpd	0xd0(%%esi)	,%%xmm6		\n\t"\
		"addpd	0xc0(%%esi)	,%%xmm7		\n\t"\
		"movaps	%%xmm2		,0x40(%%esi)\n\t"\
		"movaps	%%xmm3		,0x50(%%esi)\n\t"\
		"movaps	%%xmm6		,0xc0(%%esi)\n\t"\
		"movaps	%%xmm7		,0xd0(%%esi)\n\t"\
		"movaps	%%xmm0		,%%xmm6		\n\t"\
		"movaps	%%xmm4		,%%xmm7		\n\t"\
		"subpd	0xf0(%%esi)	,%%xmm0		\n\t"\
		"subpd	0xe0(%%esi)	,%%xmm4		\n\t"\
		"addpd	0xf0(%%esi)	,%%xmm6		\n\t"\
		"addpd	0xe0(%%esi)	,%%xmm7		\n\t"\
		"movaps	%%xmm0		,0x60(%%esi)\n\t"\
		"movaps	%%xmm4		,0x70(%%esi)\n\t"\
		"movaps	%%xmm6		,0xe0(%%esi)\n\t"\
		"movaps	%%xmm7		,0xf0(%%esi)\n\t"\
		"\n\t"\
		"movl	%[__add0],%%eax\n\t"\
		"movl	%[__p08],%%esi\n\t"\
		"shll	$3,%%esi\n\t"\
		"addl	%%esi,%%eax	/* add1 = add0+p8 */\n\t"\
		"movl	%[__p05],%%ebx\n\t"\
		"movl	%[__p06],%%ecx\n\t"\
		"movl	%[__p07],%%edx\n\t"\
		"movl	%[__p04],%%edi\n\t"\
		"shll	$3,%%ebx\n\t"\
		"shll	$3,%%ecx\n\t"\
		"shll	$3,%%edx\n\t"\
		"shll	$3,%%edi\n\t"\
		"addl	%%eax,%%ebx\n\t"\
		"addl	%%eax,%%ecx\n\t"\
		"addl	%%eax,%%edx\n\t"\
		"addl	%%edi,%%eax\n\t"\
		"/* SSE2_RADIX8_DIT_0TWIDDLE(r10): */\n\t"\
		"movl	%[__r10]	,%%esi 		\n\t"\
		"movaps	    (%%eax)	,%%xmm0		\n\t"\
		"movaps	0x10(%%eax)	,%%xmm1		\n\t"\
		"movaps	%%xmm0		,%%xmm2		\n\t"\
		"movaps	%%xmm1		,%%xmm3		\n\t"\
		"addpd	    (%%ebx)	,%%xmm2		\n\t"\
		"addpd	0x10(%%ebx)	,%%xmm3		\n\t"\
		"subpd	    (%%ebx)	,%%xmm0		\n\t"\
		"subpd	0x10(%%ebx)	,%%xmm1		\n\t"\
		"movaps	    (%%ecx)	,%%xmm4		\n\t"\
		"movaps	0x10(%%ecx)	,%%xmm5		\n\t"\
		"movaps	%%xmm4		,%%xmm6		\n\t"\
		"movaps	%%xmm5		,%%xmm7		\n\t"\
		"addpd	    (%%edx)	,%%xmm6		\n\t"\
		"addpd	0x10(%%edx)	,%%xmm7		\n\t"\
		"subpd	    (%%edx)	,%%xmm4		\n\t"\
		"subpd	0x10(%%edx)	,%%xmm5		\n\t"\
		"movaps	%%xmm6		,0xc0(%%esi)\n\t"\
		"movaps	%%xmm7		,0xd0(%%esi)\n\t"\
		"addpd	%%xmm2		,%%xmm6		\n\t"\
		"addpd	%%xmm3		,%%xmm7		\n\t"\
		"subpd	0xc0(%%esi)	,%%xmm2		\n\t"\
		"subpd	0xd0(%%esi)	,%%xmm3		\n\t"\
		"movaps	%%xmm2		,0xc0(%%esi)\n\t"\
		"movaps	%%xmm3		,0xd0(%%esi)\n\t"\
		"movaps	%%xmm4		,%%xmm2		\n\t"\
		"movaps	%%xmm5		,%%xmm3		\n\t"\
		"addpd	%%xmm0		,%%xmm5		\n\t"\
		"subpd	%%xmm3		,%%xmm0		\n\t"\
		"addpd	%%xmm1		,%%xmm4		\n\t"\
		"subpd	%%xmm2		,%%xmm1		\n\t"\
		"movl	%[__isrt2]	,%%esi		\n\t"\
		"movaps	%%xmm5		,%%xmm2		\n\t"\
		"movaps	%%xmm1		,%%xmm3		\n\t"\
		"addpd	%%xmm1		,%%xmm5		\n\t"\
		"movaps	(%%esi)		,%%xmm1		\n\t"\
		"subpd	%%xmm3		,%%xmm2		\n\t"\
		"movl	%[__r10]	,%%esi		\n\t"\
		"mulpd	%%xmm1		,%%xmm5		\n\t"\
		"mulpd	%%xmm1		,%%xmm2		\n\t"\
		"movaps	%%xmm0		,%%xmm3		\n\t"\
		"movaps	%%xmm5		,0xa0(%%esi)\n\t"\
		"movaps	%%xmm4		,%%xmm5		\n\t"\
		"addpd	%%xmm4		,%%xmm0		\n\t"\
		"movaps	%%xmm2		,0xb0(%%esi)\n\t"\
		"subpd	%%xmm5		,%%xmm3		\n\t"\
		"mulpd	%%xmm1		,%%xmm0		\n\t"\
		"mulpd	%%xmm1		,%%xmm3		\n\t"\
		"movaps	%%xmm0		,0xe0(%%esi)\n\t"\
		"movaps	%%xmm3		,0xf0(%%esi)\n\t"\
		"subl	%%edi		,%%eax		\n\t"\
		"subl	%%edi		,%%ebx		\n\t"\
		"subl	%%edi		,%%ecx		\n\t"\
		"subl	%%edi		,%%edx		\n\t"\
		"movaps	    (%%eax)	,%%xmm0		\n\t"\
		"movaps	0x10(%%eax)	,%%xmm1		\n\t"\
		"movaps	%%xmm0		,%%xmm2		\n\t"\
		"movaps	%%xmm1		,%%xmm3		\n\t"\
		"addpd	    (%%ebx)	,%%xmm2		\n\t"\
		"addpd	0x10(%%ebx)	,%%xmm3		\n\t"\
		"subpd	    (%%ebx)	,%%xmm0		\n\t"\
		"subpd	0x10(%%ebx)	,%%xmm1		\n\t"\
		"movaps	%%xmm6		,    (%%esi)\n\t"\
		"movaps	%%xmm7		,0x10(%%esi)\n\t"\
		"addpd	%%xmm6		,%%xmm6		\n\t"\
		"addpd	%%xmm7		,%%xmm7		\n\t"\
		"movaps	%%xmm6		,0x80(%%esi)\n\t"\
		"movaps	%%xmm7		,0x90(%%esi)\n\t"\
		"movaps	    (%%ecx)	,%%xmm4		\n\t"\
		"movaps	0x10(%%ecx)	,%%xmm5		\n\t"\
		"movaps	%%xmm4		,%%xmm6		\n\t"\
		"movaps	%%xmm5		,%%xmm7		\n\t"\
		"addpd	    (%%edx)	,%%xmm6		\n\t"\
		"addpd	0x10(%%edx)	,%%xmm7		\n\t"\
		"subpd	    (%%edx)	,%%xmm4		\n\t"\
		"subpd	0x10(%%edx)	,%%xmm5		\n\t"\
		"movaps	%%xmm6		,0x40(%%esi)\n\t"\
		"movaps	%%xmm7		,0x50(%%esi)\n\t"\
		"addpd	%%xmm2		,%%xmm6		\n\t"\
		"addpd	%%xmm3		,%%xmm7		\n\t"\
		"subpd	0x40(%%esi)	,%%xmm2		\n\t"\
		"subpd	0x50(%%esi)	,%%xmm3		\n\t"\
		"addpd	    (%%esi)	,%%xmm6		\n\t"\
		"addpd	0x10(%%esi)	,%%xmm7		\n\t"\
		"movaps	%%xmm6		,    (%%esi)\n\t"\
		"movaps	%%xmm7		,0x10(%%esi)\n\t"\
		"subpd	0x80(%%esi)	,%%xmm6		\n\t"\
		"subpd	0x90(%%esi)	,%%xmm7		\n\t"\
		"movaps	%%xmm6		,0x80(%%esi)\n\t"\
		"movaps	%%xmm7		,0x90(%%esi)\n\t"\
		"movaps	%%xmm4		,%%xmm6		\n\t"\
		"movaps	%%xmm5		,%%xmm7		\n\t"\
		"addpd	%%xmm0		,%%xmm5		\n\t"\
		"subpd	%%xmm7		,%%xmm0		\n\t"\
		"addpd	%%xmm1		,%%xmm4		\n\t"\
		"subpd	%%xmm6		,%%xmm1		\n\t"\
		"movaps	%%xmm5		,%%xmm6		\n\t"\
		"movaps	%%xmm1		,%%xmm7		\n\t"\
		"addpd	0xa0(%%esi)	,%%xmm5		\n\t"\
		"subpd	0xb0(%%esi)	,%%xmm1		\n\t"\
		"subpd	0xa0(%%esi)	,%%xmm6		\n\t"\
		"addpd	0xb0(%%esi)	,%%xmm7		\n\t"\
		"movaps	%%xmm5		,0x20(%%esi)\n\t"\
		"movaps	%%xmm1		,0x30(%%esi)\n\t"\
		"movaps	%%xmm6		,0xa0(%%esi)\n\t"\
		"movaps	%%xmm7		,0xb0(%%esi)\n\t"\
		"movaps	%%xmm2		,%%xmm6		\n\t"\
		"movaps	%%xmm3		,%%xmm7		\n\t"\
		"addpd	0xd0(%%esi)	,%%xmm2		\n\t"\
		"subpd	0xc0(%%esi)	,%%xmm3		\n\t"\
		"subpd	0xd0(%%esi)	,%%xmm6		\n\t"\
		"addpd	0xc0(%%esi)	,%%xmm7		\n\t"\
		"movaps	%%xmm2		,0x40(%%esi)\n\t"\
		"movaps	%%xmm3		,0x50(%%esi)\n\t"\
		"movaps	%%xmm6		,0xc0(%%esi)\n\t"\
		"movaps	%%xmm7		,0xd0(%%esi)\n\t"\
		"movaps	%%xmm0		,%%xmm6		\n\t"\
		"movaps	%%xmm4		,%%xmm7		\n\t"\
		"subpd	0xf0(%%esi)	,%%xmm0		\n\t"\
		"subpd	0xe0(%%esi)	,%%xmm4		\n\t"\
		"addpd	0xf0(%%esi)	,%%xmm6		\n\t"\
		"addpd	0xe0(%%esi)	,%%xmm7		\n\t"\
		"movaps	%%xmm0		,0x60(%%esi)\n\t"\
		"movaps	%%xmm4		,0x70(%%esi)\n\t"\
		"movaps	%%xmm6		,0xe0(%%esi)\n\t"\
		"movaps	%%xmm7		,0xf0(%%esi)\n\t"\
		"\n\t"\
		"movl	%[__add0],%%eax\n\t"\
		"movl	%[__p10],%%esi\n\t"\
		"shll	$3,%%esi\n\t"\
		"addl	%%esi,%%eax	/* add2 = add0+p10*/\n\t"\
		"movl	%[__p05],%%ebx\n\t"\
		"movl	%[__p06],%%ecx\n\t"\
		"movl	%[__p07],%%edx\n\t"\
		"movl	%[__p04],%%edi\n\t"\
		"shll	$3,%%ebx\n\t"\
		"shll	$3,%%ecx\n\t"\
		"shll	$3,%%edx\n\t"\
		"shll	$3,%%edi\n\t"\
		"addl	%%eax,%%ebx\n\t"\
		"addl	%%eax,%%ecx\n\t"\
		"addl	%%eax,%%edx\n\t"\
		"addl	%%edi,%%eax\n\t"\
		"/* SSE2_RADIX8_DIT_0TWIDDLE(r20): */\n\t"\
		"movl	%[__r20]	,%%esi 		\n\t"\
		"movaps	    (%%eax)	,%%xmm0		\n\t"\
		"movaps	0x10(%%eax)	,%%xmm1		\n\t"\
		"movaps	%%xmm0		,%%xmm2		\n\t"\
		"movaps	%%xmm1		,%%xmm3		\n\t"\
		"addpd	    (%%ebx)	,%%xmm2		\n\t"\
		"addpd	0x10(%%ebx)	,%%xmm3		\n\t"\
		"subpd	    (%%ebx)	,%%xmm0		\n\t"\
		"subpd	0x10(%%ebx)	,%%xmm1		\n\t"\
		"movaps	    (%%ecx)	,%%xmm4		\n\t"\
		"movaps	0x10(%%ecx)	,%%xmm5		\n\t"\
		"movaps	%%xmm4		,%%xmm6		\n\t"\
		"movaps	%%xmm5		,%%xmm7		\n\t"\
		"addpd	    (%%edx)	,%%xmm6		\n\t"\
		"addpd	0x10(%%edx)	,%%xmm7		\n\t"\
		"subpd	    (%%edx)	,%%xmm4		\n\t"\
		"subpd	0x10(%%edx)	,%%xmm5		\n\t"\
		"movaps	%%xmm6		,0xc0(%%esi)\n\t"\
		"movaps	%%xmm7		,0xd0(%%esi)\n\t"\
		"addpd	%%xmm2		,%%xmm6		\n\t"\
		"addpd	%%xmm3		,%%xmm7		\n\t"\
		"subpd	0xc0(%%esi)	,%%xmm2		\n\t"\
		"subpd	0xd0(%%esi)	,%%xmm3		\n\t"\
		"movaps	%%xmm2		,0xc0(%%esi)\n\t"\
		"movaps	%%xmm3		,0xd0(%%esi)\n\t"\
		"movaps	%%xmm4		,%%xmm2		\n\t"\
		"movaps	%%xmm5		,%%xmm3		\n\t"\
		"addpd	%%xmm0		,%%xmm5		\n\t"\
		"subpd	%%xmm3		,%%xmm0		\n\t"\
		"addpd	%%xmm1		,%%xmm4		\n\t"\
		"subpd	%%xmm2		,%%xmm1		\n\t"\
		"movl	%[__isrt2]	,%%esi		\n\t"\
		"movaps	%%xmm5		,%%xmm2		\n\t"\
		"movaps	%%xmm1		,%%xmm3		\n\t"\
		"addpd	%%xmm1		,%%xmm5		\n\t"\
		"movaps	(%%esi)		,%%xmm1		\n\t"\
		"subpd	%%xmm3		,%%xmm2		\n\t"\
		"movl	%[__r20]	,%%esi		\n\t"\
		"mulpd	%%xmm1		,%%xmm5		\n\t"\
		"mulpd	%%xmm1		,%%xmm2		\n\t"\
		"movaps	%%xmm0		,%%xmm3		\n\t"\
		"movaps	%%xmm5		,0xa0(%%esi)\n\t"\
		"movaps	%%xmm4		,%%xmm5		\n\t"\
		"addpd	%%xmm4		,%%xmm0		\n\t"\
		"movaps	%%xmm2		,0xb0(%%esi)\n\t"\
		"subpd	%%xmm5		,%%xmm3		\n\t"\
		"mulpd	%%xmm1		,%%xmm0		\n\t"\
		"mulpd	%%xmm1		,%%xmm3		\n\t"\
		"movaps	%%xmm0		,0xe0(%%esi)\n\t"\
		"movaps	%%xmm3		,0xf0(%%esi)\n\t"\
		"subl	%%edi		,%%eax		\n\t"\
		"subl	%%edi		,%%ebx		\n\t"\
		"subl	%%edi		,%%ecx		\n\t"\
		"subl	%%edi		,%%edx		\n\t"\
		"movaps	    (%%eax)	,%%xmm0		\n\t"\
		"movaps	0x10(%%eax)	,%%xmm1		\n\t"\
		"movaps	%%xmm0		,%%xmm2		\n\t"\
		"movaps	%%xmm1		,%%xmm3		\n\t"\
		"addpd	    (%%ebx)	,%%xmm2		\n\t"\
		"addpd	0x10(%%ebx)	,%%xmm3		\n\t"\
		"subpd	    (%%ebx)	,%%xmm0		\n\t"\
		"subpd	0x10(%%ebx)	,%%xmm1		\n\t"\
		"movaps	%%xmm6		,    (%%esi)\n\t"\
		"movaps	%%xmm7		,0x10(%%esi)\n\t"\
		"addpd	%%xmm6		,%%xmm6		\n\t"\
		"addpd	%%xmm7		,%%xmm7		\n\t"\
		"movaps	%%xmm6		,0x80(%%esi)\n\t"\
		"movaps	%%xmm7		,0x90(%%esi)\n\t"\
		"movaps	    (%%ecx)	,%%xmm4		\n\t"\
		"movaps	0x10(%%ecx)	,%%xmm5		\n\t"\
		"movaps	%%xmm4		,%%xmm6		\n\t"\
		"movaps	%%xmm5		,%%xmm7		\n\t"\
		"addpd	    (%%edx)	,%%xmm6		\n\t"\
		"addpd	0x10(%%edx)	,%%xmm7		\n\t"\
		"subpd	    (%%edx)	,%%xmm4		\n\t"\
		"subpd	0x10(%%edx)	,%%xmm5		\n\t"\
		"movaps	%%xmm6		,0x40(%%esi)\n\t"\
		"movaps	%%xmm7		,0x50(%%esi)\n\t"\
		"addpd	%%xmm2		,%%xmm6		\n\t"\
		"addpd	%%xmm3		,%%xmm7		\n\t"\
		"subpd	0x40(%%esi)	,%%xmm2		\n\t"\
		"subpd	0x50(%%esi)	,%%xmm3		\n\t"\
		"addpd	    (%%esi)	,%%xmm6		\n\t"\
		"addpd	0x10(%%esi)	,%%xmm7		\n\t"\
		"movaps	%%xmm6		,    (%%esi)\n\t"\
		"movaps	%%xmm7		,0x10(%%esi)\n\t"\
		"subpd	0x80(%%esi)	,%%xmm6		\n\t"\
		"subpd	0x90(%%esi)	,%%xmm7		\n\t"\
		"movaps	%%xmm6		,0x80(%%esi)\n\t"\
		"movaps	%%xmm7		,0x90(%%esi)\n\t"\
		"movaps	%%xmm4		,%%xmm6		\n\t"\
		"movaps	%%xmm5		,%%xmm7		\n\t"\
		"addpd	%%xmm0		,%%xmm5		\n\t"\
		"subpd	%%xmm7		,%%xmm0		\n\t"\
		"addpd	%%xmm1		,%%xmm4		\n\t"\
		"subpd	%%xmm6		,%%xmm1		\n\t"\
		"movaps	%%xmm5		,%%xmm6		\n\t"\
		"movaps	%%xmm1		,%%xmm7		\n\t"\
		"addpd	0xa0(%%esi)	,%%xmm5		\n\t"\
		"subpd	0xb0(%%esi)	,%%xmm1		\n\t"\
		"subpd	0xa0(%%esi)	,%%xmm6		\n\t"\
		"addpd	0xb0(%%esi)	,%%xmm7		\n\t"\
		"movaps	%%xmm5		,0x20(%%esi)\n\t"\
		"movaps	%%xmm1		,0x30(%%esi)\n\t"\
		"movaps	%%xmm6		,0xa0(%%esi)\n\t"\
		"movaps	%%xmm7		,0xb0(%%esi)\n\t"\
		"movaps	%%xmm2		,%%xmm6		\n\t"\
		"movaps	%%xmm3		,%%xmm7		\n\t"\
		"addpd	0xd0(%%esi)	,%%xmm2		\n\t"\
		"subpd	0xc0(%%esi)	,%%xmm3		\n\t"\
		"subpd	0xd0(%%esi)	,%%xmm6		\n\t"\
		"addpd	0xc0(%%esi)	,%%xmm7		\n\t"\
		"movaps	%%xmm2		,0x40(%%esi)\n\t"\
		"movaps	%%xmm3		,0x50(%%esi)\n\t"\
		"movaps	%%xmm6		,0xc0(%%esi)\n\t"\
		"movaps	%%xmm7		,0xd0(%%esi)\n\t"\
		"movaps	%%xmm0		,%%xmm6		\n\t"\
		"movaps	%%xmm4		,%%xmm7		\n\t"\
		"subpd	0xf0(%%esi)	,%%xmm0		\n\t"\
		"subpd	0xe0(%%esi)	,%%xmm4		\n\t"\
		"addpd	0xf0(%%esi)	,%%xmm6		\n\t"\
		"addpd	0xe0(%%esi)	,%%xmm7		\n\t"\
		"movaps	%%xmm0		,0x60(%%esi)\n\t"\
		"movaps	%%xmm4		,0x70(%%esi)\n\t"\
		"movaps	%%xmm6		,0xe0(%%esi)\n\t"\
		"movaps	%%xmm7		,0xf0(%%esi)\n\t"\
		"\n\t"\
		"movl	%[__add0],%%eax\n\t"\
		"movl	%[__p18],%%esi\n\t"\
		"shll	$3,%%esi\n\t"\
		"addl	%%esi,%%eax	/* add3 = add0+p18*/\n\t"\
		"movl	%[__p05],%%ebx\n\t"\
		"movl	%[__p06],%%ecx\n\t"\
		"movl	%[__p07],%%edx\n\t"\
		"movl	%[__p04],%%edi\n\t"\
		"shll	$3,%%ebx\n\t"\
		"shll	$3,%%ecx\n\t"\
		"shll	$3,%%edx\n\t"\
		"shll	$3,%%edi\n\t"\
		"addl	%%eax,%%ebx\n\t"\
		"addl	%%eax,%%ecx\n\t"\
		"addl	%%eax,%%edx\n\t"\
		"addl	%%edi,%%eax\n\t"\
		"/* SSE2_RADIX8_DIT_0TWIDDLE(r30): */\n\t"\
		"movl	%[__r30]	,%%esi 		\n\t"\
		"movaps	    (%%eax)	,%%xmm0		\n\t"\
		"movaps	0x10(%%eax)	,%%xmm1		\n\t"\
		"movaps	%%xmm0		,%%xmm2		\n\t"\
		"movaps	%%xmm1		,%%xmm3		\n\t"\
		"addpd	    (%%ebx)	,%%xmm2		\n\t"\
		"addpd	0x10(%%ebx)	,%%xmm3		\n\t"\
		"subpd	    (%%ebx)	,%%xmm0		\n\t"\
		"subpd	0x10(%%ebx)	,%%xmm1		\n\t"\
		"movaps	    (%%ecx)	,%%xmm4		\n\t"\
		"movaps	0x10(%%ecx)	,%%xmm5		\n\t"\
		"movaps	%%xmm4		,%%xmm6		\n\t"\
		"movaps	%%xmm5		,%%xmm7		\n\t"\
		"addpd	    (%%edx)	,%%xmm6		\n\t"\
		"addpd	0x10(%%edx)	,%%xmm7		\n\t"\
		"subpd	    (%%edx)	,%%xmm4		\n\t"\
		"subpd	0x10(%%edx)	,%%xmm5		\n\t"\
		"movaps	%%xmm6		,0xc0(%%esi)\n\t"\
		"movaps	%%xmm7		,0xd0(%%esi)\n\t"\
		"addpd	%%xmm2		,%%xmm6		\n\t"\
		"addpd	%%xmm3		,%%xmm7		\n\t"\
		"subpd	0xc0(%%esi)	,%%xmm2		\n\t"\
		"subpd	0xd0(%%esi)	,%%xmm3		\n\t"\
		"movaps	%%xmm2		,0xc0(%%esi)\n\t"\
		"movaps	%%xmm3		,0xd0(%%esi)\n\t"\
		"movaps	%%xmm4		,%%xmm2		\n\t"\
		"movaps	%%xmm5		,%%xmm3		\n\t"\
		"addpd	%%xmm0		,%%xmm5		\n\t"\
		"subpd	%%xmm3		,%%xmm0		\n\t"\
		"addpd	%%xmm1		,%%xmm4		\n\t"\
		"subpd	%%xmm2		,%%xmm1		\n\t"\
		"movl	%[__isrt2]	,%%esi		\n\t"\
		"movaps	%%xmm5		,%%xmm2		\n\t"\
		"movaps	%%xmm1		,%%xmm3		\n\t"\
		"addpd	%%xmm1		,%%xmm5		\n\t"\
		"movaps	(%%esi)		,%%xmm1		\n\t"\
		"subpd	%%xmm3		,%%xmm2		\n\t"\
		"movl	%[__r30]	,%%esi		\n\t"\
		"mulpd	%%xmm1		,%%xmm5		\n\t"\
		"mulpd	%%xmm1		,%%xmm2		\n\t"\
		"movaps	%%xmm0		,%%xmm3		\n\t"\
		"movaps	%%xmm5		,0xa0(%%esi)\n\t"\
		"movaps	%%xmm4		,%%xmm5		\n\t"\
		"addpd	%%xmm4		,%%xmm0		\n\t"\
		"movaps	%%xmm2		,0xb0(%%esi)\n\t"\
		"subpd	%%xmm5		,%%xmm3		\n\t"\
		"mulpd	%%xmm1		,%%xmm0		\n\t"\
		"mulpd	%%xmm1		,%%xmm3		\n\t"\
		"movaps	%%xmm0		,0xe0(%%esi)\n\t"\
		"movaps	%%xmm3		,0xf0(%%esi)\n\t"\
		"subl	%%edi		,%%eax		\n\t"\
		"subl	%%edi		,%%ebx		\n\t"\
		"subl	%%edi		,%%ecx		\n\t"\
		"subl	%%edi		,%%edx		\n\t"\
		"movaps	    (%%eax)	,%%xmm0		\n\t"\
		"movaps	0x10(%%eax)	,%%xmm1		\n\t"\
		"movaps	%%xmm0		,%%xmm2		\n\t"\
		"movaps	%%xmm1		,%%xmm3		\n\t"\
		"addpd	    (%%ebx)	,%%xmm2		\n\t"\
		"addpd	0x10(%%ebx)	,%%xmm3		\n\t"\
		"subpd	    (%%ebx)	,%%xmm0		\n\t"\
		"subpd	0x10(%%ebx)	,%%xmm1		\n\t"\
		"movaps	%%xmm6		,    (%%esi)\n\t"\
		"movaps	%%xmm7		,0x10(%%esi)\n\t"\
		"addpd	%%xmm6		,%%xmm6		\n\t"\
		"addpd	%%xmm7		,%%xmm7		\n\t"\
		"movaps	%%xmm6		,0x80(%%esi)\n\t"\
		"movaps	%%xmm7		,0x90(%%esi)\n\t"\
		"movaps	    (%%ecx)	,%%xmm4		\n\t"\
		"movaps	0x10(%%ecx)	,%%xmm5		\n\t"\
		"movaps	%%xmm4		,%%xmm6		\n\t"\
		"movaps	%%xmm5		,%%xmm7		\n\t"\
		"addpd	    (%%edx)	,%%xmm6		\n\t"\
		"addpd	0x10(%%edx)	,%%xmm7		\n\t"\
		"subpd	    (%%edx)	,%%xmm4		\n\t"\
		"subpd	0x10(%%edx)	,%%xmm5		\n\t"\
		"movaps	%%xmm6		,0x40(%%esi)\n\t"\
		"movaps	%%xmm7		,0x50(%%esi)\n\t"\
		"addpd	%%xmm2		,%%xmm6		\n\t"\
		"addpd	%%xmm3		,%%xmm7		\n\t"\
		"subpd	0x40(%%esi)	,%%xmm2		\n\t"\
		"subpd	0x50(%%esi)	,%%xmm3		\n\t"\
		"addpd	    (%%esi)	,%%xmm6		\n\t"\
		"addpd	0x10(%%esi)	,%%xmm7		\n\t"\
		"movaps	%%xmm6		,    (%%esi)\n\t"\
		"movaps	%%xmm7		,0x10(%%esi)\n\t"\
		"subpd	0x80(%%esi)	,%%xmm6		\n\t"\
		"subpd	0x90(%%esi)	,%%xmm7		\n\t"\
		"movaps	%%xmm6		,0x80(%%esi)\n\t"\
		"movaps	%%xmm7		,0x90(%%esi)\n\t"\
		"movaps	%%xmm4		,%%xmm6		\n\t"\
		"movaps	%%xmm5		,%%xmm7		\n\t"\
		"addpd	%%xmm0		,%%xmm5		\n\t"\
		"subpd	%%xmm7		,%%xmm0		\n\t"\
		"addpd	%%xmm1		,%%xmm4		\n\t"\
		"subpd	%%xmm6		,%%xmm1		\n\t"\
		"movaps	%%xmm5		,%%xmm6		\n\t"\
		"movaps	%%xmm1		,%%xmm7		\n\t"\
		"addpd	0xa0(%%esi)	,%%xmm5		\n\t"\
		"subpd	0xb0(%%esi)	,%%xmm1		\n\t"\
		"subpd	0xa0(%%esi)	,%%xmm6		\n\t"\
		"addpd	0xb0(%%esi)	,%%xmm7		\n\t"\
		"movaps	%%xmm5		,0x20(%%esi)\n\t"\
		"movaps	%%xmm1		,0x30(%%esi)\n\t"\
		"movaps	%%xmm6		,0xa0(%%esi)\n\t"\
		"movaps	%%xmm7		,0xb0(%%esi)\n\t"\
		"movaps	%%xmm2		,%%xmm6		\n\t"\
		"movaps	%%xmm3		,%%xmm7		\n\t"\
		"addpd	0xd0(%%esi)	,%%xmm2		\n\t"\
		"subpd	0xc0(%%esi)	,%%xmm3		\n\t"\
		"subpd	0xd0(%%esi)	,%%xmm6		\n\t"\
		"addpd	0xc0(%%esi)	,%%xmm7		\n\t"\
		"movaps	%%xmm2		,0x40(%%esi)\n\t"\
		"movaps	%%xmm3		,0x50(%%esi)\n\t"\
		"movaps	%%xmm6		,0xc0(%%esi)\n\t"\
		"movaps	%%xmm7		,0xd0(%%esi)\n\t"\
		"movaps	%%xmm0		,%%xmm6		\n\t"\
		"movaps	%%xmm4		,%%xmm7		\n\t"\
		"subpd	0xf0(%%esi)	,%%xmm0		\n\t"\
		"subpd	0xe0(%%esi)	,%%xmm4		\n\t"\
		"addpd	0xf0(%%esi)	,%%xmm6		\n\t"\
		"addpd	0xe0(%%esi)	,%%xmm7		\n\t"\
		"movaps	%%xmm0		,0x60(%%esi)\n\t"\
		"movaps	%%xmm4		,0x70(%%esi)\n\t"\
		"movaps	%%xmm6		,0xe0(%%esi)\n\t"\
		"movaps	%%xmm7		,0xf0(%%esi)\n\t"\
		"\n\t"\
		"/*...and now do eight radix-4 transforms, including the internal twiddle factors: */\n\t"\
		"\n\t"\
		"/*...Block 1: t00,t10,t20,t30	*/\n\t"\
		"movl	%[__add0],%%eax\n\t"\
		"movl	%[__p08],%%ebx\n\t"\
		"movl	%[__r00],%%ecx\n\t"\
		"movl	%[__r10],%%edx\n\t"\
		"movl	%[__p10],%%edi	/* edi will store copy of p10 throughout */\n\t"\
		"shll	$3,%%ebx\n\t"\
		"shll	$3,%%edi\n\t"\
		"addl	%%eax	,%%ebx	/* add1 = add0+p8 */\n\t"\
		"\n\t"\
		"movaps	     (%%edx),%%xmm2	/* t10 */\n\t"\
		"movaps	0x200(%%edx),%%xmm4	/* t30 */\n\t"\
		"movaps	0x010(%%edx),%%xmm3	/* t11 */\n\t"\
		"movaps	0x210(%%edx),%%xmm5	/* t31 */\n\t"\
		"movaps	     (%%ecx),%%xmm0	/* t00 */\n\t"\
		"movaps	0x200(%%ecx),%%xmm6	/* t20 */\n\t"\
		"movaps	0x010(%%ecx),%%xmm1	/* t01 */\n\t"\
		"movaps	0x210(%%ecx),%%xmm7	/* t21 */\n\t"\
		"\n\t"\
		"subpd	%%xmm2,%%xmm0	/*~t10=t00-t10*/\n\t"\
		"subpd	%%xmm4,%%xmm6	/*~t30=t20-t30*/\n\t"\
		"subpd	%%xmm3,%%xmm1	/*~t11=t01-t11*/\n\t"\
		"subpd	%%xmm5,%%xmm7	/*~t31=t21-t31*/\n\t"\
		"addpd	%%xmm2,%%xmm2	/*       2*t10*/\n\t"\
		"addpd	%%xmm4,%%xmm4	/*       2*t30*/\n\t"\
		"addpd	%%xmm3,%%xmm3	/*       2*t11*/\n\t"\
		"addpd	%%xmm5,%%xmm5	/*       2*t31*/\n\t"\
		"addpd	%%xmm0,%%xmm2	/*~t00=t00+t10*/\n\t"\
		"addpd	%%xmm6,%%xmm4	/*~t20=t20+t30*/\n\t"\
		"addpd	%%xmm1,%%xmm3	/*~t01=t01+t11*/\n\t"\
		"addpd	%%xmm7,%%xmm5	/*~t21=t21+t31*/\n\t"\
		"\n\t"\
		"/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c00): */\n\t"\
		"movl	%[__isrt2],%%esi\n\t"\
		"subpd	%%xmm4,%%xmm2\n\t"\
		"subpd	%%xmm7,%%xmm0\n\t"\
		"subpd	%%xmm5,%%xmm3\n\t"\
		"subpd	%%xmm6,%%xmm1\n\t"\
		"addpd	%%xmm4,%%xmm4\n\t"\
		"addpd	%%xmm7,%%xmm7\n\t"\
		"addpd	%%xmm5,%%xmm5\n\t"\
		"addpd	%%xmm6,%%xmm6\n\t"\
		"addpd	%%xmm2,%%xmm4\n\t"\
		"addpd	%%xmm0,%%xmm7\n\t"\
		"addpd	%%xmm3,%%xmm5\n\t"\
		"addpd	%%xmm1,%%xmm6\n\t"\
		"addl	$0x070,%%esi	/* c00 */\n\t"\
		"movaps	%%xmm2,     (%%ecx)\n\t"\
		"movaps	%%xmm0,0x010(%%edx)\n\t"\
		"movaps	%%xmm3,0x010(%%ecx)\n\t"\
		"movaps	%%xmm6,     (%%edx)\n\t"\
		"movaps	%%xmm4,%%xmm2\n\t"\
		"movaps	%%xmm7,%%xmm0\n\t"\
		"movaps	%%xmm5,%%xmm3\n\t"\
		"movaps	%%xmm1,%%xmm6\n\t"\
		"mulpd	     (%%esi),%%xmm4\n\t"\
		"mulpd	0x020(%%esi),%%xmm7\n\t"\
		"mulpd	     (%%esi),%%xmm5\n\t"\
		"mulpd	0x020(%%esi),%%xmm1\n\t"\
		"mulpd	0x010(%%esi),%%xmm2\n\t"\
		"mulpd	0x030(%%esi),%%xmm0\n\t"\
		"mulpd	0x010(%%esi),%%xmm3\n\t"\
		"mulpd	0x030(%%esi),%%xmm6\n\t"\
		"subpd	%%xmm2,%%xmm5\n\t"\
		"subpd	%%xmm0,%%xmm1\n\t"\
		"addpd	%%xmm3,%%xmm4\n\t"\
		"addpd	%%xmm6,%%xmm7\n\t"\
		"movaps	%%xmm5,0x10(%%eax)\n\t"\
		"movaps	%%xmm1,0x10(%%ebx)\n\t"\
		"movaps	%%xmm4,    (%%eax)\n\t"\
		"movaps	%%xmm7,    (%%ebx)\n\t"\
		"addl	%%edi,%%eax\n\t"\
		"addl	%%edi,%%ebx\n\t"\
		"addl	$0x40,%%esi\n\t"\
		"movaps	     (%%ecx),%%xmm4\n\t"\
		"movaps	0x010(%%edx),%%xmm0\n\t"\
		"movaps	0x010(%%ecx),%%xmm5\n\t"\
		"movaps	     (%%edx),%%xmm6\n\t"\
		"movaps	%%xmm4,%%xmm2\n\t"\
		"movaps	%%xmm0,%%xmm1\n\t"\
		"movaps	%%xmm5,%%xmm3\n\t"\
		"movaps	%%xmm6,%%xmm7\n\t"\
		"mulpd	     (%%esi),%%xmm4\n\t"\
		"mulpd	0x020(%%esi),%%xmm0\n\t"\
		"mulpd	     (%%esi),%%xmm5\n\t"\
		"mulpd	0x020(%%esi),%%xmm6\n\t"\
		"mulpd	0x010(%%esi),%%xmm2\n\t"\
		"mulpd	0x030(%%esi),%%xmm1\n\t"\
		"mulpd	0x010(%%esi),%%xmm3\n\t"\
		"mulpd	0x030(%%esi),%%xmm7\n\t"\
		"subpd	%%xmm2,%%xmm5\n\t"\
		"subpd	%%xmm1,%%xmm6\n\t"\
		"addpd	%%xmm3,%%xmm4\n\t"\
		"addpd	%%xmm7,%%xmm0\n\t"\
		"movaps	%%xmm5,0x10(%%eax)\n\t"\
		"movaps	%%xmm6,0x10(%%ebx)\n\t"\
		"movaps	%%xmm4,    (%%eax)\n\t"\
		"movaps	%%xmm0,    (%%ebx)\n\t"\
		"\n\t"\
		"/*...Block 2: t02,t12,t22,t32	*/\n\t"\
		"movl	%[__add0],%%eax\n\t"\
		"movl	%[__p01],%%esi\n\t"\
		"movl	%[__p08],%%ebx\n\t"\
		"shll	$3,%%esi\n\t"\
		"shll	$3,%%ebx\n\t"\
		"addl	%%esi	,%%eax	/* add0 = &a[j1+p1] */\n\t"\
		"addl	%%eax	,%%ebx	/* add1 = add0+p8 */\n\t"\
		"movl	%[__r02],%%ecx\n\t"\
		"movl	%[__r12],%%edx\n\t"\
		"movl	%[__isrt2],%%esi\n\t"\
		"addl	$0x200,%%ecx\n\t"\
		"addl	$0x200,%%edx\n\t"\
		"\n\t"\
		"addl	$0x030,%%esi	/* cc1 */\n\t"\
		"movaps	     (%%ecx),%%xmm4	/* t22 */\n\t"\
		"movaps	0x010(%%ecx),%%xmm5	/* t23 */\n\t"\
		"movaps	     (%%esi),%%xmm2	/* c32_1 */\n\t"\
		"movaps	0x010(%%esi),%%xmm3	/* s32_1 */\n\t"\
		"movaps	%%xmm4,%%xmm6	/* xmm2 <- cpy t22 */\n\t"\
		"movaps	%%xmm5,%%xmm7	/* xmm3 <- cpy t23 */\n\t"\
		"\n\t"\
		"mulpd	%%xmm2,%%xmm4	/* t22*c32_1 */\n\t"\
		"mulpd	%%xmm2,%%xmm5	/* t23*c32_1 */\n\t"\
		"mulpd	%%xmm3,%%xmm6	/* t22*s32_1 */\n\t"\
		"movaps	     (%%edx),%%xmm0	/* t32 */\n\t"\
		"mulpd	%%xmm3,%%xmm7	/* t23*s32_1 */\n\t"\
		"movaps	0x010(%%edx),%%xmm1	/* t33 */\n\t"\
		"addl	$0x020,%%esi	/* cc3 */\n\t"\
		"movaps	     (%%esi),%%xmm2	/* c32_3 */\n\t"\
		"movaps	0x010(%%esi),%%xmm3	/* s32_3 */\n\t"\
		"subpd	%%xmm6,%%xmm5	/* xmm5 <-~t23 */\n\t"\
		"movaps	%%xmm0,%%xmm6	/* xmm6 <- cpy t32 */\n\t"\
		"addpd	%%xmm7,%%xmm4	/* xmm4 <-~t22 */\n\t"\
		"movaps	%%xmm1,%%xmm7	/* xmm7 <- cpy t33 */\n\t"\
		"\n\t"\
		"mulpd	%%xmm2,%%xmm0	/* t32*c32_3 */\n\t"\
		"mulpd	%%xmm2,%%xmm1	/* t33*c32_3 */\n\t"\
		"mulpd	%%xmm3,%%xmm6	/* t32*s32_3 */\n\t"\
		"mulpd	%%xmm3,%%xmm7	/* t33*s32_3 */\n\t"\
		"subpd	%%xmm6,%%xmm1	/* xmm1 <- it */\n\t"\
		"addpd	%%xmm7,%%xmm0	/* xmm0 <- rt */\n\t"\
		"\n\t"\
		"movaps	%%xmm5,%%xmm7	/* xmm7 <- cpy~t23*/\n\t"\
		"movaps	%%xmm4,%%xmm6	/* xmm6 <- cpy~t22*/\n\t"\
		"\n\t"\
		"addpd	%%xmm0,%%xmm4	/* ~t22 <- t22+rt */\n\t"\
		"addpd	%%xmm1,%%xmm5	/* ~t23 <- t23+it */\n\t"\
		"subpd	%%xmm0,%%xmm6	/* ~t32 <- t22-rt */\n\t"\
		"subpd	%%xmm1,%%xmm7	/* ~t33 <- t23-it */\n\t"\
		"\n\t"\
		"subl	$0x200,%%ecx\n\t"\
		"subl	$0x200,%%edx\n\t"\
		"subl	$0x040,%%esi	/* cc0 */\n\t"\
		"movaps	     (%%edx),%%xmm1	/* t12 */\n\t"\
		"movaps	0x010(%%edx),%%xmm3	/* t13 */\n\t"\
		"movaps	0x010(%%esi),%%xmm2	/* s */\n\t"\
		"movaps	%%xmm1,%%xmm0	/* cpy t12 */\n\t"\
		"mulpd	%%xmm2,%%xmm1	/* t12*s */\n\t"\
		"mulpd	%%xmm3,%%xmm2	/* t13*s */\n\t"\
		"mulpd	(%%esi),%%xmm0	/* t12*c */\n\t"\
		"mulpd	(%%esi),%%xmm3	/* t13*c */\n\t"\
		"\n\t"\
		"addpd	%%xmm0,%%xmm2	/* rt =t12*c + t13*s */\n\t"\
		"subpd	%%xmm1,%%xmm3	/* it =t13*c - t12*s */\n\t"\
		"\n\t"\
		"movaps	     (%%ecx),%%xmm0	/* t02 */\n\t"\
		"movaps	0x010(%%ecx),%%xmm1	/* t03 */\n\t"\
		"subpd	%%xmm2,%%xmm0	/*~t12 <- t02- rt */\n\t"\
		"subpd	%%xmm3,%%xmm1	/*~t13 <- t03- it */\n\t"\
		"addpd	%%xmm2,%%xmm2	/*          2* rt */\n\t"\
		"addpd	%%xmm3,%%xmm3	/*          2* it */\n\t"\
		"addpd	%%xmm0,%%xmm2	/*~t02 <- t02+ rt */\n\t"\
		"addpd	%%xmm1,%%xmm3	/*~t03 <- t03+ it */\n\t"\
		"\n\t"\
		"/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c01): */\n\t"\
		"addl	$0x260,%%esi	/* c01 */\n\t"\
		"subpd	%%xmm4,%%xmm2\n\t"\
		"subpd	%%xmm7,%%xmm0\n\t"\
		"subpd	%%xmm5,%%xmm3\n\t"\
		"subpd	%%xmm6,%%xmm1\n\t"\
		"addpd	%%xmm4,%%xmm4\n\t"\
		"addpd	%%xmm7,%%xmm7\n\t"\
		"addpd	%%xmm5,%%xmm5\n\t"\
		"addpd	%%xmm6,%%xmm6\n\t"\
		"addpd	%%xmm2,%%xmm4\n\t"\
		"addpd	%%xmm0,%%xmm7\n\t"\
		"addpd	%%xmm3,%%xmm5\n\t"\
		"addpd	%%xmm1,%%xmm6\n\t"\
		"movaps	%%xmm2,     (%%ecx)\n\t"\
		"movaps	%%xmm0,0x010(%%edx)\n\t"\
		"movaps	%%xmm3,0x010(%%ecx)\n\t"\
		"movaps	%%xmm6,     (%%edx)\n\t"\
		"movaps	%%xmm4,%%xmm2\n\t"\
		"movaps	%%xmm7,%%xmm0\n\t"\
		"movaps	%%xmm5,%%xmm3\n\t"\
		"movaps	%%xmm1,%%xmm6\n\t"\
		"mulpd	     (%%esi),%%xmm4\n\t"\
		"mulpd	0x020(%%esi),%%xmm7\n\t"\
		"mulpd	     (%%esi),%%xmm5\n\t"\
		"mulpd	0x020(%%esi),%%xmm1\n\t"\
		"mulpd	0x010(%%esi),%%xmm2\n\t"\
		"mulpd	0x030(%%esi),%%xmm0\n\t"\
		"mulpd	0x010(%%esi),%%xmm3\n\t"\
		"mulpd	0x030(%%esi),%%xmm6\n\t"\
		"subpd	%%xmm2,%%xmm5\n\t"\
		"subpd	%%xmm0,%%xmm1\n\t"\
		"addpd	%%xmm3,%%xmm4\n\t"\
		"addpd	%%xmm6,%%xmm7\n\t"\
		"movaps	%%xmm5,0x10(%%eax)\n\t"\
		"movaps	%%xmm1,0x10(%%ebx)\n\t"\
		"movaps	%%xmm4,    (%%eax)\n\t"\
		"movaps	%%xmm7,    (%%ebx)\n\t"\
		"addl	%%edi,%%eax\n\t"\
		"addl	%%edi,%%ebx\n\t"\
		"addl	$0x40,%%esi\n\t"\
		"movaps	     (%%ecx),%%xmm4\n\t"\
		"movaps	0x010(%%edx),%%xmm0\n\t"\
		"movaps	0x010(%%ecx),%%xmm5\n\t"\
		"movaps	     (%%edx),%%xmm6\n\t"\
		"movaps	%%xmm4,%%xmm2\n\t"\
		"movaps	%%xmm0,%%xmm1\n\t"\
		"movaps	%%xmm5,%%xmm3\n\t"\
		"movaps	%%xmm6,%%xmm7\n\t"\
		"mulpd	     (%%esi),%%xmm4\n\t"\
		"mulpd	0x020(%%esi),%%xmm0\n\t"\
		"mulpd	     (%%esi),%%xmm5\n\t"\
		"mulpd	0x020(%%esi),%%xmm6\n\t"\
		"mulpd	0x010(%%esi),%%xmm2\n\t"\
		"mulpd	0x030(%%esi),%%xmm1\n\t"\
		"mulpd	0x010(%%esi),%%xmm3\n\t"\
		"mulpd	0x030(%%esi),%%xmm7\n\t"\
		"subpd	%%xmm2,%%xmm5\n\t"\
		"subpd	%%xmm1,%%xmm6\n\t"\
		"addpd	%%xmm3,%%xmm4\n\t"\
		"addpd	%%xmm7,%%xmm0\n\t"\
		"movaps	%%xmm5,0x10(%%eax)\n\t"\
		"movaps	%%xmm6,0x10(%%ebx)\n\t"\
		"movaps	%%xmm4,    (%%eax)\n\t"\
		"movaps	%%xmm0,    (%%ebx)\n\t"\
		"\n\t"\
		"/*...Block 3: t04,t14,t24,t34	*/\n\t"\
		"movl	%[__add0],%%eax\n\t"\
		"movl	%[__p02],%%esi\n\t"\
		"movl	%[__p08],%%ebx\n\t"\
		"shll	$3,%%esi\n\t"\
		"shll	$3,%%ebx\n\t"\
		"addl	%%esi	,%%eax	/* add0 = &a[j1+p2] */\n\t"\
		"addl	%%eax	,%%ebx	/* add1 = add0+p8 */\n\t"\
		"movl	%[__r04],%%ecx\n\t"\
		"movl	%[__r14],%%edx\n\t"\
		"movl	%[__isrt2],%%esi\n\t"\
		"addl	$0x200,%%ecx	/* t24 */\n\t"\
		"addl	$0x200,%%edx	/* t25 */\n\t"\
		"addl	$0x010,%%esi	/* cc0 */\n\t"\
		"\n\t"\
		"movaps	     (%%ecx),%%xmm4	/* t24 */\n\t"\
		"movaps	0x010(%%ecx),%%xmm5	/* t25 */\n\t"\
		"movaps	     (%%esi),%%xmm2	/* c */\n\t"\
		"movaps	0x010(%%esi),%%xmm3	/* s */\n\t"\
		"movaps	%%xmm4,%%xmm6	/* xmm2 <- cpy t24 */\n\t"\
		"movaps	%%xmm5,%%xmm7	/* xmm3 <- cpy t25 */\n\t"\
		"\n\t"\
		"mulpd	%%xmm2,%%xmm4	/* t24*c */\n\t"\
		"mulpd	%%xmm2,%%xmm5	/* t25*c */\n\t"\
		"mulpd	%%xmm3,%%xmm6	/* t24*s */\n\t"\
		"movaps	     (%%edx),%%xmm0	/* t34 */\n\t"\
		"mulpd	%%xmm3,%%xmm7	/* t25*s */\n\t"\
		"movaps	0x010(%%edx),%%xmm1	/* t35 */\n\t"\
		"subpd	%%xmm6,%%xmm5	/* xmm1 <-~t25 */\n\t"\
		"movaps	%%xmm0,%%xmm6	/* xmm6 <- cpy t34 */\n\t"\
		"addpd	%%xmm7,%%xmm4	/* xmm0 <-~t24 */\n\t"\
		"movaps	%%xmm1,%%xmm7	/* xmm7 <- cpy t35 */\n\t"\
		"\n\t"\
		"mulpd	%%xmm3,%%xmm0	/* t34*s */\n\t"\
		"mulpd	%%xmm3,%%xmm1	/* t35*s */\n\t"\
		"mulpd	%%xmm2,%%xmm6	/* t34*c */\n\t"\
		"mulpd	%%xmm2,%%xmm7	/* t35*c */\n\t"\
		"subpd	%%xmm6,%%xmm1	/* xmm5 <- it */\n\t"\
		"addpd	%%xmm7,%%xmm0	/* xmm4 <- rt */\n\t"\
		"\n\t"\
		"movaps	%%xmm5,%%xmm7	/* xmm3 <- cpy~t25*/\n\t"\
		"movaps	%%xmm4,%%xmm6	/* xmm2 <- cpy~t24*/\n\t"\
		"\n\t"\
		"addpd	%%xmm0,%%xmm4	/* ~t24 <- t24+rt */\n\t"\
		"addpd	%%xmm1,%%xmm5	/* ~t25 <- t25+it */\n\t"\
		"subpd	%%xmm0,%%xmm6	/* ~t34 <- t24-rt */\n\t"\
		"subpd	%%xmm1,%%xmm7	/* ~t35 <- t25-it */\n\t"\
		"\n\t"\
		"subl	$0x200,%%ecx\n\t"\
		"subl	$0x200,%%edx\n\t"\
		"subl	$0x10,%%esi	/* isrt2 */\n\t"\
		"movaps	     (%%edx),%%xmm2	/* t14 */\n\t"\
		"movaps	0x010(%%edx),%%xmm3	/* t15 */\n\t"\
		"movaps	(%%esi),%%xmm1	/* isrt2 */\n\t"\
		"movaps	%%xmm3,%%xmm0	/* cpy t15 */\n\t"\
		"subpd	%%xmm2,%%xmm3	/*~t15=t15-t14 */\n\t"\
		"addpd	%%xmm0,%%xmm2	/*~t14=t14+t15 */\n\t"\
		"mulpd	%%xmm1,%%xmm2	/* rt */\n\t"\
		"mulpd	%%xmm1,%%xmm3	/* it */\n\t"\
		"\n\t"\
		"movaps	     (%%ecx),%%xmm0	/* t04 */\n\t"\
		"movaps	0x010(%%ecx),%%xmm1	/* t05 */\n\t"\
		"subpd	%%xmm2,%%xmm0	/*~t14 <- t04- rt */\n\t"\
		"subpd	%%xmm3,%%xmm1	/*~t15 <- t05- it */\n\t"\
		"addpd	%%xmm2,%%xmm2	/*          2* rt */\n\t"\
		"addpd	%%xmm3,%%xmm3	/*          2* it */\n\t"\
		"addpd	%%xmm0,%%xmm2	/*~t04 <- t04+ rt */\n\t"\
		"addpd	%%xmm1,%%xmm3	/*~t05 <- t05+ it */\n\t"\
		"\n\t"\
		"/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c02): */\n\t"\
		"addl	$0x170,%%esi	/* c02 */\n\t"\
		"subpd	%%xmm4,%%xmm2\n\t"\
		"subpd	%%xmm7,%%xmm0\n\t"\
		"subpd	%%xmm5,%%xmm3\n\t"\
		"subpd	%%xmm6,%%xmm1\n\t"\
		"addpd	%%xmm4,%%xmm4\n\t"\
		"addpd	%%xmm7,%%xmm7\n\t"\
		"addpd	%%xmm5,%%xmm5\n\t"\
		"addpd	%%xmm6,%%xmm6\n\t"\
		"addpd	%%xmm2,%%xmm4\n\t"\
		"addpd	%%xmm0,%%xmm7\n\t"\
		"addpd	%%xmm3,%%xmm5\n\t"\
		"addpd	%%xmm1,%%xmm6\n\t"\
		"movaps	%%xmm2,     (%%ecx)\n\t"\
		"movaps	%%xmm0,0x010(%%edx)\n\t"\
		"movaps	%%xmm3,0x010(%%ecx)\n\t"\
		"movaps	%%xmm6,     (%%edx)\n\t"\
		"movaps	%%xmm4,%%xmm2\n\t"\
		"movaps	%%xmm7,%%xmm0\n\t"\
		"movaps	%%xmm5,%%xmm3\n\t"\
		"movaps	%%xmm1,%%xmm6\n\t"\
		"mulpd	     (%%esi),%%xmm4\n\t"\
		"mulpd	0x020(%%esi),%%xmm7\n\t"\
		"mulpd	     (%%esi),%%xmm5\n\t"\
		"mulpd	0x020(%%esi),%%xmm1\n\t"\
		"mulpd	0x010(%%esi),%%xmm2\n\t"\
		"mulpd	0x030(%%esi),%%xmm0\n\t"\
		"mulpd	0x010(%%esi),%%xmm3\n\t"\
		"mulpd	0x030(%%esi),%%xmm6\n\t"\
		"subpd	%%xmm2,%%xmm5\n\t"\
		"subpd	%%xmm0,%%xmm1\n\t"\
		"addpd	%%xmm3,%%xmm4\n\t"\
		"addpd	%%xmm6,%%xmm7\n\t"\
		"movaps	%%xmm5,0x10(%%eax)\n\t"\
		"movaps	%%xmm1,0x10(%%ebx)\n\t"\
		"movaps	%%xmm4,    (%%eax)\n\t"\
		"movaps	%%xmm7,    (%%ebx)\n\t"\
		"addl	%%edi,%%eax\n\t"\
		"addl	%%edi,%%ebx\n\t"\
		"addl	$0x40,%%esi\n\t"\
		"movaps	     (%%ecx),%%xmm4\n\t"\
		"movaps	0x010(%%edx),%%xmm0\n\t"\
		"movaps	0x010(%%ecx),%%xmm5\n\t"\
		"movaps	     (%%edx),%%xmm6\n\t"\
		"movaps	%%xmm4,%%xmm2\n\t"\
		"movaps	%%xmm0,%%xmm1\n\t"\
		"movaps	%%xmm5,%%xmm3\n\t"\
		"movaps	%%xmm6,%%xmm7\n\t"\
		"mulpd	     (%%esi),%%xmm4\n\t"\
		"mulpd	0x020(%%esi),%%xmm0\n\t"\
		"mulpd	     (%%esi),%%xmm5\n\t"\
		"mulpd	0x020(%%esi),%%xmm6\n\t"\
		"mulpd	0x010(%%esi),%%xmm2\n\t"\
		"mulpd	0x030(%%esi),%%xmm1\n\t"\
		"mulpd	0x010(%%esi),%%xmm3\n\t"\
		"mulpd	0x030(%%esi),%%xmm7\n\t"\
		"subpd	%%xmm2,%%xmm5\n\t"\
		"subpd	%%xmm1,%%xmm6\n\t"\
		"addpd	%%xmm3,%%xmm4\n\t"\
		"addpd	%%xmm7,%%xmm0\n\t"\
		"movaps	%%xmm5,0x10(%%eax)\n\t"\
		"movaps	%%xmm6,0x10(%%ebx)\n\t"\
		"movaps	%%xmm4,    (%%eax)\n\t"\
		"movaps	%%xmm0,    (%%ebx)\n\t"\
		"\n\t"\
		"/*...Block 4: t06,t16,t26,t36	*/\n\t"\
		"movl	%[__add0],%%eax\n\t"\
		"movl	%[__p03],%%esi\n\t"\
		"movl	%[__p08],%%ebx\n\t"\
		"shll	$3,%%esi\n\t"\
		"shll	$3,%%ebx\n\t"\
		"addl	%%esi	,%%eax	/* add0 = &a[j1+p3] */\n\t"\
		"addl	%%eax	,%%ebx	/* add1 = add0+p8 */\n\t"\
		"movl	%[__r06],%%ecx\n\t"\
		"movl	%[__r16],%%edx\n\t"\
		"movl	%[__isrt2],%%esi\n\t"\
		"addl	$0x200,%%ecx\n\t"\
		"addl	$0x200,%%edx\n\t"\
		"addl	$0x050,%%esi	/* cc3 */\n\t"\
		"\n\t"\
		"movaps	     (%%ecx),%%xmm4	/* t26 */\n\t"\
		"movaps	0x010(%%ecx),%%xmm5	/* t27 */\n\t"\
		"movaps	     (%%esi),%%xmm2	/* c32_3 */\n\t"\
		"movaps	0x010(%%esi),%%xmm3	/* s32_3 */\n\t"\
		"movaps	%%xmm4,%%xmm6	/* xmm2 <- cpy t26 */\n\t"\
		"movaps	%%xmm5,%%xmm7	/* xmm3 <- cpy t27 */\n\t"\
		"\n\t"\
		"mulpd	%%xmm2,%%xmm4	/* t26*c32_3 */\n\t"\
		"mulpd	%%xmm2,%%xmm5	/* t27*c32_3 */\n\t"\
		"mulpd	%%xmm3,%%xmm6	/* t26*s32_3 */\n\t"\
		"movaps	     (%%edx),%%xmm0	/* t36 */\n\t"\
		"mulpd	%%xmm3,%%xmm7	/* t27*s32_3 */\n\t"\
		"movaps	0x010(%%edx),%%xmm1	/* t37 */\n\t"\
		"subl	$0x20,%%esi	/* cc1 */\n\t"\
		"movaps	     (%%esi),%%xmm3	/* c32_1 */\n\t"\
		"movaps	0x010(%%esi),%%xmm2	/* s32_1 */\n\t"\
		"subpd	%%xmm6,%%xmm5	/* xmm5 <-~t27 */\n\t"\
		"movaps	%%xmm0,%%xmm6	/* xmm6 <- cpy t36 */\n\t"\
		"addpd	%%xmm7,%%xmm4	/* xmm4 <-~t26 */\n\t"\
		"movaps	%%xmm1,%%xmm7	/* xmm7 <- cpy t37 */\n\t"\
		"\n\t"\
		"mulpd	%%xmm2,%%xmm0	/* t36*s32_1 */\n\t"\
		"mulpd	%%xmm2,%%xmm1	/* t37*s32_1 */\n\t"\
		"mulpd	%%xmm3,%%xmm6	/* t36*c32_1 */\n\t"\
		"mulpd	%%xmm3,%%xmm7	/* t37*c32_1 */\n\t"\
		"addpd	%%xmm6,%%xmm1	/* xmm1 <- it */\n\t"\
		"subpd	%%xmm7,%%xmm0	/* xmm0 <- rt */\n\t"\
		"\n\t"\
		"movaps	%%xmm5,%%xmm7	/* xmm7 <- cpy~t27*/\n\t"\
		"movaps	%%xmm4,%%xmm6	/* xmm6 <- cpy~t26*/\n\t"\
		"\n\t"\
		"addpd	%%xmm0,%%xmm6	/* ~t36 <- t26+rt */\n\t"\
		"addpd	%%xmm1,%%xmm7	/* ~t37 <- t27+it */\n\t"\
		"subpd	%%xmm0,%%xmm4	/* ~t26 <- t26-rt */\n\t"\
		"subpd	%%xmm1,%%xmm5	/* ~t27 <- t27-it */\n\t"\
		"\n\t"\
		"subl	$0x200,%%ecx\n\t"\
		"subl	$0x200,%%edx\n\t"\
		"subl	$0x20,%%esi	/* cc0 */\n\t"\
		"movaps	     (%%edx),%%xmm2	/* t16 */\n\t"\
		"movaps	0x010(%%edx),%%xmm0	/* t17 */\n\t"\
		"movaps	0x010(%%esi),%%xmm3	/* s */\n\t"\
		"movaps	%%xmm2,%%xmm1	/* cpy t16 */\n\t"\
		"mulpd	%%xmm3,%%xmm2	/* t16*s */\n\t"\
		"mulpd	%%xmm0,%%xmm3	/* s*t17 */\n\t"\
		"mulpd	(%%esi),%%xmm1	/* t16*c */\n\t"\
		"mulpd	(%%esi),%%xmm0	/* t17*c */\n\t"\
		"\n\t"\
		"addpd	%%xmm0,%%xmm2	/* rt =t16*s - t17*c */\n\t"\
		"subpd	%%xmm1,%%xmm3	/* it =t17*s + t16*c */\n\t"\
		"\n\t"\
		"movaps	     (%%ecx),%%xmm0	/* t06 */\n\t"\
		"movaps	0x010(%%ecx),%%xmm1	/* t07 */\n\t"\
		"subpd	%%xmm2,%%xmm0	/*~t16 <- t06- rt */\n\t"\
		"subpd	%%xmm3,%%xmm1	/*~t17 <- t07- it */\n\t"\
		"addpd	%%xmm2,%%xmm2	/*          2* rt */\n\t"\
		"addpd	%%xmm3,%%xmm3	/*          2* it */\n\t"\
		"addpd	%%xmm0,%%xmm2	/*~t06 <- t06+ rt */\n\t"\
		"addpd	%%xmm1,%%xmm3	/*~t07 <- t07+ it */\n\t"\
		"\n\t"\
		"/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c03): */\n\t"\
		"addl	$0x360,%%esi	/* c03 */\n\t"\
		"subpd	%%xmm4,%%xmm2\n\t"\
		"subpd	%%xmm7,%%xmm0\n\t"\
		"subpd	%%xmm5,%%xmm3\n\t"\
		"subpd	%%xmm6,%%xmm1\n\t"\
		"addpd	%%xmm4,%%xmm4\n\t"\
		"addpd	%%xmm7,%%xmm7\n\t"\
		"addpd	%%xmm5,%%xmm5\n\t"\
		"addpd	%%xmm6,%%xmm6\n\t"\
		"addpd	%%xmm2,%%xmm4\n\t"\
		"addpd	%%xmm0,%%xmm7\n\t"\
		"addpd	%%xmm3,%%xmm5\n\t"\
		"addpd	%%xmm1,%%xmm6\n\t"\
		"movaps	%%xmm2,     (%%ecx)\n\t"\
		"movaps	%%xmm0,0x010(%%edx)\n\t"\
		"movaps	%%xmm3,0x010(%%ecx)\n\t"\
		"movaps	%%xmm6,     (%%edx)\n\t"\
		"movaps	%%xmm4,%%xmm2\n\t"\
		"movaps	%%xmm7,%%xmm0\n\t"\
		"movaps	%%xmm5,%%xmm3\n\t"\
		"movaps	%%xmm1,%%xmm6\n\t"\
		"mulpd	     (%%esi),%%xmm4\n\t"\
		"mulpd	0x020(%%esi),%%xmm7\n\t"\
		"mulpd	     (%%esi),%%xmm5\n\t"\
		"mulpd	0x020(%%esi),%%xmm1\n\t"\
		"mulpd	0x010(%%esi),%%xmm2\n\t"\
		"mulpd	0x030(%%esi),%%xmm0\n\t"\
		"mulpd	0x010(%%esi),%%xmm3\n\t"\
		"mulpd	0x030(%%esi),%%xmm6\n\t"\
		"subpd	%%xmm2,%%xmm5\n\t"\
		"subpd	%%xmm0,%%xmm1\n\t"\
		"addpd	%%xmm3,%%xmm4\n\t"\
		"addpd	%%xmm6,%%xmm7\n\t"\
		"movaps	%%xmm5,0x10(%%eax)\n\t"\
		"movaps	%%xmm1,0x10(%%ebx)\n\t"\
		"movaps	%%xmm4,    (%%eax)\n\t"\
		"movaps	%%xmm7,    (%%ebx)\n\t"\
		"addl	%%edi,%%eax\n\t"\
		"addl	%%edi,%%ebx\n\t"\
		"addl	$0x40,%%esi	/* c0B */\n\t"\
		"movaps	     (%%ecx),%%xmm4\n\t"\
		"movaps	0x010(%%edx),%%xmm0\n\t"\
		"movaps	0x010(%%ecx),%%xmm5\n\t"\
		"movaps	     (%%edx),%%xmm6\n\t"\
		"movaps	%%xmm4,%%xmm2\n\t"\
		"movaps	%%xmm0,%%xmm1\n\t"\
		"movaps	%%xmm5,%%xmm3\n\t"\
		"movaps	%%xmm6,%%xmm7\n\t"\
		"mulpd	     (%%esi),%%xmm4\n\t"\
		"mulpd	0x020(%%esi),%%xmm0\n\t"\
		"mulpd	     (%%esi),%%xmm5\n\t"\
		"mulpd	0x020(%%esi),%%xmm6\n\t"\
		"mulpd	0x010(%%esi),%%xmm2\n\t"\
		"mulpd	0x030(%%esi),%%xmm1\n\t"\
		"mulpd	0x010(%%esi),%%xmm3\n\t"\
		"mulpd	0x030(%%esi),%%xmm7\n\t"\
		"subpd	%%xmm2,%%xmm5\n\t"\
		"subpd	%%xmm1,%%xmm6\n\t"\
		"addpd	%%xmm3,%%xmm4\n\t"\
		"addpd	%%xmm7,%%xmm0\n\t"\
		"movaps	%%xmm5,0x10(%%eax)\n\t"\
		"movaps	%%xmm6,0x10(%%ebx)\n\t"\
		"movaps	%%xmm4,    (%%eax)\n\t"\
		"movaps	%%xmm0,    (%%ebx)\n\t"\
		"\n\t"\
		"/*...Block 5: t08,t18,t28,t38	*/\n\t"\
		"movl	%[__add0],%%eax\n\t"\
		"movl	%[__p04],%%esi\n\t"\
		"movl	%[__p08],%%ebx\n\t"\
		"shll	$3,%%esi\n\t"\
		"shll	$3,%%ebx\n\t"\
		"addl	%%esi	,%%eax	/* add0 = &a[j1+p4] */\n\t"\
		"addl	%%eax	,%%ebx	/* add1 = add0+p8 */\n\t"\
		"movl	%[__r08],%%ecx\n\t"\
		"movl	%[__r18],%%edx\n\t"\
		"movl	%[__isrt2],%%esi\n\t"\
		"addl	$0x200,%%ecx\n\t"\
		"addl	$0x200,%%edx\n\t"\
		"movaps	(%%esi),%%xmm2	/* isrt2 */\n\t"\
		"\n\t"\
		"movaps	     (%%ecx),%%xmm4	/* t28 */\n\t"\
		"movaps	0x010(%%ecx),%%xmm5	/* t29 */\n\t"\
		"movaps	     (%%edx),%%xmm6	/* t38 */\n\t"\
		"movaps	0x010(%%edx),%%xmm7	/* t39 */\n\t"\
		"subl	$0x200,%%ecx\n\t"\
		"subl	$0x200,%%edx\n\t"\
		"mulpd	%%xmm2,%%xmm4\n\t"\
		"mulpd	%%xmm2,%%xmm5\n\t"\
		"mulpd	%%xmm2,%%xmm6\n\t"\
		"mulpd	%%xmm2,%%xmm7\n\t"\
		"\n\t"\
		"subpd	%%xmm4,%%xmm5	/*~t29=t29-t28*/\n\t"\
		"movaps	     (%%ecx),%%xmm0	/* t08 */\n\t"\
		"subpd	%%xmm7,%%xmm6	/* rt =t38-t39*/\n\t"\
		"movaps	0x010(%%edx),%%xmm2	/* t19 */\n\t"\
		"addpd	%%xmm4,%%xmm4	/*       2*t28*/\n\t"\
		"movaps	0x010(%%ecx),%%xmm3	/* t09 */\n\t"\
		"addpd	%%xmm7,%%xmm7	/*       2*t39*/\n\t"\
		"movaps	     (%%edx),%%xmm1	/* t18 */\n\t"\
		"addpd	%%xmm5,%%xmm4	/*~t28=t28+t29*/\n\t"\
		"addpd	%%xmm6,%%xmm7	/* it =t39+t38*/\n\t"\
		"\n\t"\
		"subpd	%%xmm6,%%xmm4	/*~t28=t28-rt */\n\t"\
		"subpd	%%xmm2,%%xmm0	/*~t18=t08-t19*/\n\t"\
		"subpd	%%xmm7,%%xmm5	/*~t29=t29-it */\n\t"\
		"subpd	%%xmm1,%%xmm3	/*~t09=t09-t18*/\n\t"\
		"addpd	%%xmm6,%%xmm6	/*       2*rt */\n\t"\
		"addpd	%%xmm2,%%xmm2	/*       2*t08*/\n\t"\
		"addpd	%%xmm7,%%xmm7	/*       2*it */\n\t"\
		"addpd	%%xmm1,%%xmm1	/*       2*t09*/\n\t"\
		"addpd	%%xmm4,%%xmm6	/*~t38=t28+rt */\n\t"\
		"addpd	%%xmm0,%%xmm2	/*~t08=t19+t08*/\n\t"\
		"addpd	%%xmm5,%%xmm7	/*~t39=t29+it */\n\t"\
		"addpd	%%xmm3,%%xmm1	/*~t19=t18+t09*/\n\t"\
		"\n\t"\
		"/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c04): */\n\t"\
		"addl	$0x0f0,%%esi	/* c04 */\n\t"\
		"subpd	%%xmm4,%%xmm2\n\t"\
		"subpd	%%xmm7,%%xmm0\n\t"\
		"subpd	%%xmm5,%%xmm3\n\t"\
		"subpd	%%xmm6,%%xmm1\n\t"\
		"addpd	%%xmm4,%%xmm4\n\t"\
		"addpd	%%xmm7,%%xmm7\n\t"\
		"addpd	%%xmm5,%%xmm5\n\t"\
		"addpd	%%xmm6,%%xmm6\n\t"\
		"addpd	%%xmm2,%%xmm4\n\t"\
		"addpd	%%xmm0,%%xmm7\n\t"\
		"addpd	%%xmm3,%%xmm5\n\t"\
		"addpd	%%xmm1,%%xmm6\n\t"\
		"movaps	%%xmm2,     (%%ecx)\n\t"\
		"movaps	%%xmm0,0x010(%%edx)\n\t"\
		"movaps	%%xmm3,0x010(%%ecx)\n\t"\
		"movaps	%%xmm6,     (%%edx)\n\t"\
		"movaps	%%xmm4,%%xmm2\n\t"\
		"movaps	%%xmm7,%%xmm0\n\t"\
		"movaps	%%xmm5,%%xmm3\n\t"\
		"movaps	%%xmm1,%%xmm6\n\t"\
		"mulpd	     (%%esi),%%xmm4\n\t"\
		"mulpd	0x020(%%esi),%%xmm7\n\t"\
		"mulpd	     (%%esi),%%xmm5\n\t"\
		"mulpd	0x020(%%esi),%%xmm1\n\t"\
		"mulpd	0x010(%%esi),%%xmm2\n\t"\
		"mulpd	0x030(%%esi),%%xmm0\n\t"\
		"mulpd	0x010(%%esi),%%xmm3\n\t"\
		"mulpd	0x030(%%esi),%%xmm6\n\t"\
		"subpd	%%xmm2,%%xmm5\n\t"\
		"subpd	%%xmm0,%%xmm1\n\t"\
		"addpd	%%xmm3,%%xmm4\n\t"\
		"addpd	%%xmm6,%%xmm7\n\t"\
		"movaps	%%xmm5,0x10(%%eax)\n\t"\
		"movaps	%%xmm1,0x10(%%ebx)\n\t"\
		"movaps	%%xmm4,    (%%eax)\n\t"\
		"movaps	%%xmm7,    (%%ebx)\n\t"\
		"addl	%%edi,%%eax\n\t"\
		"addl	%%edi,%%ebx\n\t"\
		"addl	$0x40,%%esi	/* c0C */\n\t"\
		"movaps	     (%%ecx),%%xmm4\n\t"\
		"movaps	0x010(%%edx),%%xmm0\n\t"\
		"movaps	0x010(%%ecx),%%xmm5\n\t"\
		"movaps	     (%%edx),%%xmm6\n\t"\
		"movaps	%%xmm4,%%xmm2\n\t"\
		"movaps	%%xmm0,%%xmm1\n\t"\
		"movaps	%%xmm5,%%xmm3\n\t"\
		"movaps	%%xmm6,%%xmm7\n\t"\
		"mulpd	     (%%esi),%%xmm4\n\t"\
		"mulpd	0x020(%%esi),%%xmm0\n\t"\
		"mulpd	     (%%esi),%%xmm5\n\t"\
		"mulpd	0x020(%%esi),%%xmm6\n\t"\
		"mulpd	0x010(%%esi),%%xmm2\n\t"\
		"mulpd	0x030(%%esi),%%xmm1\n\t"\
		"mulpd	0x010(%%esi),%%xmm3\n\t"\
		"mulpd	0x030(%%esi),%%xmm7\n\t"\
		"subpd	%%xmm2,%%xmm5\n\t"\
		"subpd	%%xmm1,%%xmm6\n\t"\
		"addpd	%%xmm3,%%xmm4\n\t"\
		"addpd	%%xmm7,%%xmm0\n\t"\
		"movaps	%%xmm5,0x10(%%eax)\n\t"\
		"movaps	%%xmm6,0x10(%%ebx)\n\t"\
		"movaps	%%xmm4,    (%%eax)\n\t"\
		"movaps	%%xmm0,    (%%ebx)\n\t"\
		"\n\t"\
		"/*...Block 6: t0A,t1A,t2A,t3A	*/\n\t"\
		"movl	%[__add0],%%eax\n\t"\
		"movl	%[__p05],%%esi\n\t"\
		"movl	%[__p08],%%ebx\n\t"\
		"shll	$3,%%esi\n\t"\
		"shll	$3,%%ebx\n\t"\
		"addl	%%esi	,%%eax	/* add0 = &a[j1+p5] */\n\t"\
		"addl	%%eax	,%%ebx	/* add1 = add0+p8 */\n\t"\
		"movl	%[__r0A],%%ecx\n\t"\
		"movl	%[__r1A],%%edx\n\t"\
		"movl	%[__isrt2],%%esi\n\t"\
		"addl	$0x200,%%ecx\n\t"\
		"addl	$0x200,%%edx\n\t"\
		"addl	$0x050,%%esi	/* cc3 */\n\t"\
		"\n\t"\
		"movaps	     (%%ecx),%%xmm4	/* t2A */\n\t"\
		"movaps	0x010(%%ecx),%%xmm5	/* t2B */\n\t"\
		"movaps	     (%%esi),%%xmm3	/* c32_3 */\n\t"\
		"movaps	0x010(%%esi),%%xmm2	/* s32_3 */\n\t"\
		"movaps	%%xmm4,%%xmm6	/* xmm2 <- cpy t2A */\n\t"\
		"movaps	%%xmm5,%%xmm7	/* xmm3 <- cpy t2B */\n\t"\
		"\n\t"\
		"mulpd	%%xmm2,%%xmm4	/* t2A*s32_3 */\n\t"\
		"mulpd	%%xmm2,%%xmm5	/* t2B*s32_3 */\n\t"\
		"mulpd	%%xmm3,%%xmm6	/* t2A*c32_3 */\n\t"\
		"movaps	     (%%edx),%%xmm0	/* t3A */\n\t"\
		"mulpd	%%xmm3,%%xmm7	/* t2B*c32_3 */\n\t"\
		"movaps	0x010(%%edx),%%xmm1	/* t3B */\n\t"\
		"subl	$0x20,%%esi	/* cc1 */\n\t"\
		"movaps	     (%%esi),%%xmm2	/* c32_1 */\n\t"\
		"movaps	0x010(%%esi),%%xmm3	/* s32_1 */\n\t"\
		"subpd	%%xmm6,%%xmm5	/* xmm5 <-~t2B */\n\t"\
		"movaps	%%xmm0,%%xmm6	/* xmm6 <- cpy t3A */\n\t"\
		"addpd	%%xmm7,%%xmm4	/* xmm4 <-~t2A */\n\t"\
		"movaps	%%xmm1,%%xmm7	/* xmm7 <- cpy t3B */\n\t"\
		"\n\t"\
		"mulpd	%%xmm2,%%xmm0	/* t3A*c32_1 */\n\t"\
		"mulpd	%%xmm2,%%xmm1	/* t3B*c32_1 */\n\t"\
		"mulpd	%%xmm3,%%xmm6	/* t3A*s32_1 */\n\t"\
		"mulpd	%%xmm3,%%xmm7	/* t3B*s32_1 */\n\t"\
		"addpd	%%xmm6,%%xmm1	/* xmm1 <- it */\n\t"\
		"subpd	%%xmm7,%%xmm0	/* xmm0 <- rt */\n\t"\
		"\n\t"\
		"movaps	%%xmm5,%%xmm7	/* xmm7 <- cpy~t2B*/\n\t"\
		"movaps	%%xmm4,%%xmm6	/* xmm6 <- cpy~t2A*/\n\t"\
		"\n\t"\
		"addpd	%%xmm0,%%xmm6	/* ~t3A <- t2A+rt */\n\t"\
		"addpd	%%xmm1,%%xmm7	/* ~t3B <- t2B+it */\n\t"\
		"subpd	%%xmm0,%%xmm4	/* ~t2A <- t2A-rt */\n\t"\
		"subpd	%%xmm1,%%xmm5	/* ~t2B <- t2B-it */\n\t"\
		"\n\t"\
		"subl	$0x200,%%ecx\n\t"\
		"subl	$0x200,%%edx\n\t"\
		"subl	$0x20,%%esi	/* cc0 */\n\t"\
		"movaps	     (%%edx),%%xmm0	/* t1A */\n\t"\
		"movaps	0x010(%%edx),%%xmm2	/* t1B */\n\t"\
		"movaps	0x010(%%esi),%%xmm1	/* s */\n\t"\
		"movaps	%%xmm0,%%xmm3	/* cpy t1A */\n\t"\
		"mulpd	%%xmm1,%%xmm0	/* t1A*s */\n\t"\
		"mulpd	%%xmm2,%%xmm1	/* s*t1B */\n\t"\
		"mulpd	(%%esi),%%xmm3	/* t1A*c */\n\t"\
		"mulpd	(%%esi),%%xmm2	/* t1B*c */\n\t"\
		"\n\t"\
		"subpd	%%xmm2,%%xmm0	/* rt =t1A*s - t1B*c */\n\t"\
		"addpd	%%xmm3,%%xmm1	/* it =t1B*s + t1A*c */\n\t"\
		"\n\t"\
		"movaps	     (%%ecx),%%xmm2	/* t0A */\n\t"\
		"movaps	0x010(%%ecx),%%xmm3	/* t0B */\n\t"\
		"subpd	%%xmm0,%%xmm2	/*~t0A <- t0A- rt */\n\t"\
		"subpd	%%xmm1,%%xmm3	/*~t0B <- t0B- it */\n\t"\
		"addpd	%%xmm0,%%xmm0	/*          2* rt */\n\t"\
		"addpd	%%xmm1,%%xmm1	/*          2* it */\n\t"\
		"addpd	%%xmm2,%%xmm0	/*~t1A <- t0A+ rt */\n\t"\
		"addpd	%%xmm3,%%xmm1	/*~t1B <- t0B+ it */\n\t"\
		"\n\t"\
		"/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c05): */\n\t"\
		"addl	$0x2e0,%%esi	/* c05 */\n\t"\
		"subpd	%%xmm4,%%xmm2\n\t"\
		"subpd	%%xmm7,%%xmm0\n\t"\
		"subpd	%%xmm5,%%xmm3\n\t"\
		"subpd	%%xmm6,%%xmm1\n\t"\
		"addpd	%%xmm4,%%xmm4\n\t"\
		"addpd	%%xmm7,%%xmm7\n\t"\
		"addpd	%%xmm5,%%xmm5\n\t"\
		"addpd	%%xmm6,%%xmm6\n\t"\
		"addpd	%%xmm2,%%xmm4\n\t"\
		"addpd	%%xmm0,%%xmm7\n\t"\
		"addpd	%%xmm3,%%xmm5\n\t"\
		"addpd	%%xmm1,%%xmm6\n\t"\
		"movaps	%%xmm2,     (%%ecx)\n\t"\
		"movaps	%%xmm0,0x010(%%edx)\n\t"\
		"movaps	%%xmm3,0x010(%%ecx)\n\t"\
		"movaps	%%xmm6,     (%%edx)\n\t"\
		"movaps	%%xmm4,%%xmm2\n\t"\
		"movaps	%%xmm7,%%xmm0\n\t"\
		"movaps	%%xmm5,%%xmm3\n\t"\
		"movaps	%%xmm1,%%xmm6\n\t"\
		"mulpd	     (%%esi),%%xmm4\n\t"\
		"mulpd	0x020(%%esi),%%xmm7\n\t"\
		"mulpd	     (%%esi),%%xmm5\n\t"\
		"mulpd	0x020(%%esi),%%xmm1\n\t"\
		"mulpd	0x010(%%esi),%%xmm2\n\t"\
		"mulpd	0x030(%%esi),%%xmm0\n\t"\
		"mulpd	0x010(%%esi),%%xmm3\n\t"\
		"mulpd	0x030(%%esi),%%xmm6\n\t"\
		"subpd	%%xmm2,%%xmm5\n\t"\
		"subpd	%%xmm0,%%xmm1\n\t"\
		"addpd	%%xmm3,%%xmm4\n\t"\
		"addpd	%%xmm6,%%xmm7\n\t"\
		"movaps	%%xmm5,0x10(%%eax)\n\t"\
		"movaps	%%xmm1,0x10(%%ebx)\n\t"\
		"movaps	%%xmm4,    (%%eax)\n\t"\
		"movaps	%%xmm7,    (%%ebx)\n\t"\
		"addl	%%edi,%%eax\n\t"\
		"addl	%%edi,%%ebx\n\t"\
		"addl	$0x40,%%esi	/* c0D */\n\t"\
		"movaps	     (%%ecx),%%xmm4\n\t"\
		"movaps	0x010(%%edx),%%xmm0\n\t"\
		"movaps	0x010(%%ecx),%%xmm5\n\t"\
		"movaps	     (%%edx),%%xmm6\n\t"\
		"movaps	%%xmm4,%%xmm2\n\t"\
		"movaps	%%xmm0,%%xmm1\n\t"\
		"movaps	%%xmm5,%%xmm3\n\t"\
		"movaps	%%xmm6,%%xmm7\n\t"\
		"mulpd	     (%%esi),%%xmm4\n\t"\
		"mulpd	0x020(%%esi),%%xmm0\n\t"\
		"mulpd	     (%%esi),%%xmm5\n\t"\
		"mulpd	0x020(%%esi),%%xmm6\n\t"\
		"mulpd	0x010(%%esi),%%xmm2\n\t"\
		"mulpd	0x030(%%esi),%%xmm1\n\t"\
		"mulpd	0x010(%%esi),%%xmm3\n\t"\
		"mulpd	0x030(%%esi),%%xmm7\n\t"\
		"subpd	%%xmm2,%%xmm5\n\t"\
		"subpd	%%xmm1,%%xmm6\n\t"\
		"addpd	%%xmm3,%%xmm4\n\t"\
		"addpd	%%xmm7,%%xmm0\n\t"\
		"movaps	%%xmm5,0x10(%%eax)\n\t"\
		"movaps	%%xmm6,0x10(%%ebx)\n\t"\
		"movaps	%%xmm4,    (%%eax)\n\t"\
		"movaps	%%xmm0,    (%%ebx)\n\t"\
		"\n\t"\
		"/*...Block 7: t0C,t1C,t2C,t3C	*/\n\t"\
		"movl	%[__add0],%%eax\n\t"\
		"movl	%[__p06],%%esi\n\t"\
		"movl	%[__p08],%%ebx\n\t"\
		"shll	$3,%%esi\n\t"\
		"shll	$3,%%ebx\n\t"\
		"addl	%%esi	,%%eax	/* add0 = &a[j1+p6] */\n\t"\
		"addl	%%eax	,%%ebx	/* add1 = add0+p8 */\n\t"\
		"movl	%[__r0C],%%ecx\n\t"\
		"movl	%[__r1C],%%edx\n\t"\
		"movl	%[__isrt2],%%esi\n\t"\
		"addl	$0x200,%%ecx\n\t"\
		"addl	$0x200,%%edx\n\t"\
		"addl	$0x010,%%esi	/* cc0 */\n\t"\
		"\n\t"\
		"movaps	     (%%ecx),%%xmm4	/* t2C */\n\t"\
		"movaps	0x010(%%ecx),%%xmm5	/* t2D */\n\t"\
		"movaps	     (%%esi),%%xmm3	/* c */\n\t"\
		"movaps	0x010(%%esi),%%xmm2	/* s */\n\t"\
		"movaps	%%xmm4,%%xmm6	/* xmm2 <- cpy t2C */\n\t"\
		"movaps	%%xmm5,%%xmm7	/* xmm3 <- cpy t2D */\n\t"\
		"\n\t"\
		"mulpd	%%xmm2,%%xmm4	/* t2C*s */\n\t"\
		"mulpd	%%xmm2,%%xmm5	/* t2D*s */\n\t"\
		"mulpd	%%xmm3,%%xmm6	/* t2C*c */\n\t"\
		"movaps	     (%%edx),%%xmm0	/* t3C */\n\t"\
		"mulpd	%%xmm3,%%xmm7	/* t2D*c */\n\t"\
		"movaps	0x010(%%edx),%%xmm1	/* t3D */\n\t"\
		"subpd	%%xmm6,%%xmm5	/* xmm5 <-~t2D */\n\t"\
		"movaps	%%xmm0,%%xmm6	/* xmm6 <- cpy t3C */\n\t"\
		"addpd	%%xmm7,%%xmm4	/* xmm4 <-~t2C */\n\t"\
		"movaps	%%xmm1,%%xmm7	/* xmm7 <- cpy t3D */\n\t"\
		"\n\t"\
		"mulpd	%%xmm3,%%xmm0	/* t3C*c */\n\t"\
		"mulpd	%%xmm3,%%xmm1	/* t3D*c */\n\t"\
		"mulpd	%%xmm2,%%xmm6	/* t3C*s */\n\t"\
		"mulpd	%%xmm2,%%xmm7	/* t3D*s */\n\t"\
		"subpd	%%xmm6,%%xmm1	/* xmm1 <- it */\n\t"\
		"addpd	%%xmm7,%%xmm0	/* xmm0 <- rt */\n\t"\
		"\n\t"\
		"movaps	%%xmm5,%%xmm7	/* xmm7 <- cpy~t2D*/\n\t"\
		"movaps	%%xmm4,%%xmm6	/* xmm6 <- cpy~t2C*/\n\t"\
		"\n\t"\
		"addpd	%%xmm0,%%xmm6	/* ~t3C <- t2C+rt */\n\t"\
		"addpd	%%xmm1,%%xmm7	/* ~t3D <- t2D+it */\n\t"\
		"subpd	%%xmm0,%%xmm4	/* ~t2C <- t2C-rt */\n\t"\
		"subpd	%%xmm1,%%xmm5	/* ~t2D <- t2D-it */\n\t"\
		"\n\t"\
		"subl	$0x200,%%ecx\n\t"\
		"subl	$0x200,%%edx\n\t"\
		"subl	$0x10,%%esi	/* isrt2 */\n\t"\
		"movaps	     (%%edx),%%xmm0	/* t1C */\n\t"\
		"movaps	0x010(%%edx),%%xmm1	/* t1D */\n\t"\
		"movaps	(%%esi),%%xmm3	/* isrt2 */\n\t"\
		"movaps	%%xmm0,%%xmm2	/* cpy t1C */\n\t"\
		"subpd	%%xmm1,%%xmm0	/*~t1C=t1C-t1D */\n\t"\
		"addpd	%%xmm2,%%xmm1	/*~t1D=t1D+t1C */\n\t"\
		"mulpd	%%xmm3,%%xmm0	/* it */\n\t"\
		"mulpd	%%xmm3,%%xmm1	/* rt */\n\t"\
		"\n\t"\
		"movaps	     (%%ecx),%%xmm2	/* t0C */\n\t"\
		"movaps	0x010(%%ecx),%%xmm3	/* t0D */\n\t"\
		"subpd	%%xmm0,%%xmm2	/*~t0C <- t0C- rt */\n\t"\
		"subpd	%%xmm1,%%xmm3	/*~t0D <- t0D- it */\n\t"\
		"addpd	%%xmm0,%%xmm0	/*          2* rt */\n\t"\
		"addpd	%%xmm1,%%xmm1	/*          2* it */\n\t"\
		"addpd	%%xmm2,%%xmm0	/*~t1C <- t0C+ rt */\n\t"\
		"addpd	%%xmm3,%%xmm1	/*~t1D <- t0D+ it */\n\t"\
		"\n\t"\
		"/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c06): */\n\t"\
		"addl	$0x1f0,%%esi	/* c06 */\n\t"\
		"subpd	%%xmm4,%%xmm2\n\t"\
		"subpd	%%xmm7,%%xmm0\n\t"\
		"subpd	%%xmm5,%%xmm3\n\t"\
		"subpd	%%xmm6,%%xmm1\n\t"\
		"addpd	%%xmm4,%%xmm4\n\t"\
		"addpd	%%xmm7,%%xmm7\n\t"\
		"addpd	%%xmm5,%%xmm5\n\t"\
		"addpd	%%xmm6,%%xmm6\n\t"\
		"addpd	%%xmm2,%%xmm4\n\t"\
		"addpd	%%xmm0,%%xmm7\n\t"\
		"addpd	%%xmm3,%%xmm5\n\t"\
		"addpd	%%xmm1,%%xmm6\n\t"\
		"movaps	%%xmm2,     (%%ecx)\n\t"\
		"movaps	%%xmm0,0x010(%%edx)\n\t"\
		"movaps	%%xmm3,0x010(%%ecx)\n\t"\
		"movaps	%%xmm6,     (%%edx)\n\t"\
		"movaps	%%xmm4,%%xmm2\n\t"\
		"movaps	%%xmm7,%%xmm0\n\t"\
		"movaps	%%xmm5,%%xmm3\n\t"\
		"movaps	%%xmm1,%%xmm6\n\t"\
		"mulpd	     (%%esi),%%xmm4\n\t"\
		"mulpd	0x020(%%esi),%%xmm7\n\t"\
		"mulpd	     (%%esi),%%xmm5\n\t"\
		"mulpd	0x020(%%esi),%%xmm1\n\t"\
		"mulpd	0x010(%%esi),%%xmm2\n\t"\
		"mulpd	0x030(%%esi),%%xmm0\n\t"\
		"mulpd	0x010(%%esi),%%xmm3\n\t"\
		"mulpd	0x030(%%esi),%%xmm6\n\t"\
		"subpd	%%xmm2,%%xmm5\n\t"\
		"subpd	%%xmm0,%%xmm1\n\t"\
		"addpd	%%xmm3,%%xmm4\n\t"\
		"addpd	%%xmm6,%%xmm7\n\t"\
		"movaps	%%xmm5,0x10(%%eax)\n\t"\
		"movaps	%%xmm1,0x10(%%ebx)\n\t"\
		"movaps	%%xmm4,    (%%eax)\n\t"\
		"movaps	%%xmm7,    (%%ebx)\n\t"\
		"addl	%%edi,%%eax\n\t"\
		"addl	%%edi,%%ebx\n\t"\
		"addl	$0x40,%%esi	/* c0E */\n\t"\
		"movaps	     (%%ecx),%%xmm4\n\t"\
		"movaps	0x010(%%edx),%%xmm0\n\t"\
		"movaps	0x010(%%ecx),%%xmm5\n\t"\
		"movaps	     (%%edx),%%xmm6\n\t"\
		"movaps	%%xmm4,%%xmm2\n\t"\
		"movaps	%%xmm0,%%xmm1\n\t"\
		"movaps	%%xmm5,%%xmm3\n\t"\
		"movaps	%%xmm6,%%xmm7\n\t"\
		"mulpd	     (%%esi),%%xmm4\n\t"\
		"mulpd	0x020(%%esi),%%xmm0\n\t"\
		"mulpd	     (%%esi),%%xmm5\n\t"\
		"mulpd	0x020(%%esi),%%xmm6\n\t"\
		"mulpd	0x010(%%esi),%%xmm2\n\t"\
		"mulpd	0x030(%%esi),%%xmm1\n\t"\
		"mulpd	0x010(%%esi),%%xmm3\n\t"\
		"mulpd	0x030(%%esi),%%xmm7\n\t"\
		"subpd	%%xmm2,%%xmm5\n\t"\
		"subpd	%%xmm1,%%xmm6\n\t"\
		"addpd	%%xmm3,%%xmm4\n\t"\
		"addpd	%%xmm7,%%xmm0\n\t"\
		"movaps	%%xmm5,0x10(%%eax)\n\t"\
		"movaps	%%xmm6,0x10(%%ebx)\n\t"\
		"movaps	%%xmm4,    (%%eax)\n\t"\
		"movaps	%%xmm0,    (%%ebx)\n\t"\
		"\n\t"\
		"/*...Block 8: t0E,t1E,t2E,t3E	*/\n\t"\
		"movl	%[__add0],%%eax\n\t"\
		"movl	%[__p07],%%esi\n\t"\
		"movl	%[__p08],%%ebx\n\t"\
		"shll	$3,%%esi\n\t"\
		"shll	$3	,%%ebx	/* Pointer offset for floating doubles */\n\t"\
		"addl	%%esi	,%%eax	/* add0 = &a[j1+p7] */\n\t"\
		"addl	%%eax	,%%ebx	/* add1 = add0+p8 */\n\t"\
		"movl	%[__r0E],%%ecx\n\t"\
		"movl	%[__r1E],%%edx\n\t"\
		"movl	%[__isrt2],%%esi\n\t"\
		"addl	$0x200,%%ecx\n\t"\
		"addl	$0x200,%%edx\n\t"\
		"addl	$0x030,%%esi	/* cc1 */\n\t"\
		"\n\t"\
		"movaps	     (%%ecx),%%xmm4	/* t2E */\n\t"\
		"movaps	0x010(%%ecx),%%xmm5	/* t2F */\n\t"\
		"movaps	     (%%esi),%%xmm3	/* c32_1 */\n\t"\
		"movaps	0x010(%%esi),%%xmm2	/* s32_1 */\n\t"\
		"movaps	%%xmm4,%%xmm6	/* xmm2 <- cpy t2E */\n\t"\
		"movaps	%%xmm5,%%xmm7	/* xmm3 <- cpy t2F */\n\t"\
		"\n\t"\
		"mulpd	%%xmm2,%%xmm4	/* t2E*s32_1 */\n\t"\
		"mulpd	%%xmm2,%%xmm5	/* t2F*s32_1 */\n\t"\
		"mulpd	%%xmm3,%%xmm6	/* t2E*c32_1 */\n\t"\
		"movaps	     (%%edx),%%xmm0	/* t3E */\n\t"\
		"mulpd	%%xmm3,%%xmm7	/* t2F*c32_1 */\n\t"\
		"movaps	0x010(%%edx),%%xmm1	/* t3F */\n\t"\
		"addl	$0x20,%%esi	/* cc3 */\n\t"\
		"movaps	     (%%esi),%%xmm3	/* c32_3 */\n\t"\
		"movaps	0x010(%%esi),%%xmm2	/* s32_3 */\n\t"\
		"subpd	%%xmm6,%%xmm5	/* xmm5 <-~t2F */\n\t"\
		"movaps	%%xmm0,%%xmm6	/* xmm6 <- cpy t3E */\n\t"\
		"addpd	%%xmm7,%%xmm4	/* xmm4 <-~t2E */\n\t"\
		"movaps	%%xmm1,%%xmm7	/* xmm7 <- cpy t3F */\n\t"\
		"\n\t"\
		"mulpd	%%xmm2,%%xmm0	/* t3E*s32_3 */\n\t"\
		"mulpd	%%xmm2,%%xmm1	/* t3F*s32_3 */\n\t"\
		"mulpd	%%xmm3,%%xmm6	/* t3E*c32_3 */\n\t"\
		"mulpd	%%xmm3,%%xmm7	/* t3F*c32_3 */\n\t"\
		"subpd	%%xmm6,%%xmm1	/* xmm1 <- it */\n\t"\
		"addpd	%%xmm7,%%xmm0	/* xmm0 <- rt */\n\t"\
		"\n\t"\
		"movaps	%%xmm5,%%xmm7	/* xmm7 <- cpy~t2F*/\n\t"\
		"movaps	%%xmm4,%%xmm6	/* xmm6 <- cpy~t2E*/\n\t"\
		"\n\t"\
		"subpd	%%xmm0,%%xmm4	/* ~t2E <- t2E-rt */\n\t"\
		"subpd	%%xmm1,%%xmm5	/* ~t2F <- t2F-it */\n\t"\
		"addpd	%%xmm0,%%xmm6	/* ~t3E <- t2E+rt */\n\t"\
		"addpd	%%xmm1,%%xmm7	/* ~t3F <- t2F+it */\n\t"\
		"\n\t"\
		"subl	$0x200,%%ecx\n\t"\
		"subl	$0x200,%%edx\n\t"\
		"subl	$0x40,%%esi	 /* cc0 */\n\t"\
		"movaps	     (%%edx),%%xmm3	/* t1E */\n\t"\
		"movaps	0x010(%%edx),%%xmm1	/* t1F */\n\t"\
		"movaps	0x010(%%esi),%%xmm2	/* s */\n\t"\
		"movaps	%%xmm3,%%xmm0	/* cpy t1E */\n\t"\
		"mulpd	%%xmm2,%%xmm3	/* t1E*s */\n\t"\
		"mulpd	%%xmm1,%%xmm2	/* t1F*s */\n\t"\
		"mulpd	(%%esi),%%xmm0	/* t1E*c */\n\t"\
		"mulpd	(%%esi),%%xmm1	/* t1F*c */\n\t"\
		"\n\t"\
		"subpd	%%xmm2,%%xmm0	/* rt =t1E*c - t1F*s */\n\t"\
		"addpd	%%xmm3,%%xmm1	/* it =t1F*c + t1E*s */\n\t"\
		"\n\t"\
		"movaps	     (%%ecx),%%xmm2	/* t0E */\n\t"\
		"movaps	0x010(%%ecx),%%xmm3	/* t0F */\n\t"\
		"subpd	%%xmm0,%%xmm2	/*~t0E <- t0E- rt */\n\t"\
		"subpd	%%xmm1,%%xmm3	/*~t0F <- t0F- it */\n\t"\
		"addpd	%%xmm0,%%xmm0	/*          2* rt */\n\t"\
		"addpd	%%xmm1,%%xmm1	/*          2* it */\n\t"\
		"addpd	%%xmm2,%%xmm0	/*~t1E <- t0E+ rt */\n\t"\
		"addpd	%%xmm3,%%xmm1	/*~t1F <- t0F+ it */\n\t"\
		"/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c07): */\n\t"\
		"addl	$0x3e0,%%esi	/* c07 */\n\t"\
		"subpd	%%xmm4,%%xmm2\n\t"\
		"subpd	%%xmm7,%%xmm0\n\t"\
		"subpd	%%xmm5,%%xmm3\n\t"\
		"subpd	%%xmm6,%%xmm1\n\t"\
		"addpd	%%xmm4,%%xmm4\n\t"\
		"addpd	%%xmm7,%%xmm7\n\t"\
		"addpd	%%xmm5,%%xmm5\n\t"\
		"addpd	%%xmm6,%%xmm6\n\t"\
		"addpd	%%xmm2,%%xmm4\n\t"\
		"addpd	%%xmm0,%%xmm7\n\t"\
		"addpd	%%xmm3,%%xmm5\n\t"\
		"addpd	%%xmm1,%%xmm6\n\t"\
		"movaps	%%xmm2,     (%%ecx)\n\t"\
		"movaps	%%xmm0,0x010(%%edx)\n\t"\
		"movaps	%%xmm3,0x010(%%ecx)\n\t"\
		"movaps	%%xmm6,     (%%edx)\n\t"\
		"movaps	%%xmm4,%%xmm2\n\t"\
		"movaps	%%xmm7,%%xmm0\n\t"\
		"movaps	%%xmm5,%%xmm3\n\t"\
		"movaps	%%xmm1,%%xmm6\n\t"\
		"mulpd	     (%%esi),%%xmm4\n\t"\
		"mulpd	0x020(%%esi),%%xmm7\n\t"\
		"mulpd	     (%%esi),%%xmm5\n\t"\
		"mulpd	0x020(%%esi),%%xmm1\n\t"\
		"mulpd	0x010(%%esi),%%xmm2\n\t"\
		"mulpd	0x030(%%esi),%%xmm0\n\t"\
		"mulpd	0x010(%%esi),%%xmm3\n\t"\
		"mulpd	0x030(%%esi),%%xmm6\n\t"\
		"subpd	%%xmm2,%%xmm5\n\t"\
		"subpd	%%xmm0,%%xmm1\n\t"\
		"addpd	%%xmm3,%%xmm4\n\t"\
		"addpd	%%xmm6,%%xmm7\n\t"\
		"movaps	%%xmm5,0x10(%%eax)\n\t"\
		"movaps	%%xmm1,0x10(%%ebx)\n\t"\
		"movaps	%%xmm4,    (%%eax)\n\t"\
		"movaps	%%xmm7,    (%%ebx)\n\t"\
		"addl	%%edi,%%eax\n\t"\
		"addl	%%edi,%%ebx\n\t"\
		"addl	$0x40,%%esi	/* c0F */\n\t"\
		"movaps	     (%%ecx),%%xmm4\n\t"\
		"movaps	0x010(%%edx),%%xmm0\n\t"\
		"movaps	0x010(%%ecx),%%xmm5\n\t"\
		"movaps	     (%%edx),%%xmm6\n\t"\
		"movaps	%%xmm4,%%xmm2\n\t"\
		"movaps	%%xmm0,%%xmm1\n\t"\
		"movaps	%%xmm5,%%xmm3\n\t"\
		"movaps	%%xmm6,%%xmm7\n\t"\
		"mulpd	     (%%esi),%%xmm4\n\t"\
		"mulpd	0x020(%%esi),%%xmm0\n\t"\
		"mulpd	     (%%esi),%%xmm5\n\t"\
		"mulpd	0x020(%%esi),%%xmm6\n\t"\
		"mulpd	0x010(%%esi),%%xmm2\n\t"\
		"mulpd	0x030(%%esi),%%xmm1\n\t"\
		"mulpd	0x010(%%esi),%%xmm3\n\t"\
		"mulpd	0x030(%%esi),%%xmm7\n\t"\
		"subpd	%%xmm2,%%xmm5\n\t"\
		"subpd	%%xmm1,%%xmm6\n\t"\
		"addpd	%%xmm3,%%xmm4\n\t"\
		"addpd	%%xmm7,%%xmm0\n\t"\
		"movaps	%%xmm5,0x10(%%eax)\n\t"\
		"movaps	%%xmm6,0x10(%%ebx)\n\t"\
		"movaps	%%xmm4,    (%%eax)\n\t"\
		"movaps	%%xmm0,    (%%ebx)\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__p01] "m" (Xp01)\
		 ,[__p02] "m" (Xp02)\
		 ,[__p03] "m" (Xp03)\
		 ,[__p04] "m" (Xp04)\
		 ,[__p05] "m" (Xp05)\
		 ,[__p06] "m" (Xp06)\
		 ,[__p07] "m" (Xp07)\
		 ,[__p08] "m" (Xp08)\
		 ,[__p10] "m" (Xp10)\
		 ,[__p18] "m" (Xp18)\
		 ,[__r00] "m" (Xr00)\
		 ,[__r02] "m" (Xr02)\
		 ,[__r04] "m" (Xr04)\
		 ,[__r06] "m" (Xr06)\
		 ,[__r08] "m" (Xr08)\
		 ,[__r0A] "m" (Xr0A)\
		 ,[__r0C] "m" (Xr0C)\
		 ,[__r0E] "m" (Xr0E)\
		 ,[__r10] "m" (Xr10)\
		 ,[__r12] "m" (Xr12)\
		 ,[__r14] "m" (Xr14)\
		 ,[__r16] "m" (Xr16)\
		 ,[__r18] "m" (Xr18)\
		 ,[__r1A] "m" (Xr1A)\
		 ,[__r1C] "m" (Xr1C)\
		 ,[__r1E] "m" (Xr1E)\
		 ,[__r20] "m" (Xr20)\
		 ,[__r30] "m" (Xr30)\
		 ,[__isrt2] "m" (Xisrt2)\
		: "eax","ebx","ecx","edx","edi","esi"		/* Clobbered registers */\
	);\
	}

#endif	/* radix32_dif_dit_pass_gcc_h_included */

