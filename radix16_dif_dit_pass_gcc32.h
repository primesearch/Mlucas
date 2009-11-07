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
#ifndef radix16_dif_dit_pass_gcc_h_included
#define radix16_dif_dit_pass_gcc_h_included

	#define SSE2_RADIX16_DIT_TWIDDLE(Xadd0,Xp1,Xp2,Xp3,Xp4,Xp8,Xr1,Xr3,Xr5,Xr7,Xr9,Xr11,Xr13,Xr15,Xr17,Xr25,Xisrt2,Xcc0,Xc0,Xc1,Xc2,Xc3)\
	{\
	__asm__ volatile (\
		"movl	%[__add0]	,%%eax\n\t"\
		"movl	%[__p1]		,%%ebx\n\t"\
		"movl	%[__p2]		,%%ecx\n\t"\
		"movl	%[__p3]		,%%edx\n\t"\
		"movl	%[__p4]		,%%edi\n\t"\
		"shll	$3			,%%ebx\n\t"\
		"shll	$3			,%%ecx\n\t"\
		"shll	$3			,%%edx\n\t"\
		"shll	$3			,%%edi\n\t"\
		"addl	%%eax		,%%ebx\n\t"\
		"addl	%%eax		,%%ecx\n\t"\
		"addl	%%eax		,%%edx\n\t"\
		"/* SSE2_RADIX4_DIT_0TWIDDLE_B(r1) */\n\t"\
		"movaps	    (%%eax),%%xmm2			\n\t	movaps	    (%%ecx),%%xmm6\n\t"\
		"movaps	0x10(%%eax),%%xmm3			\n\t	movaps	0x10(%%ecx),%%xmm7\n\t"\
		"movaps	    (%%ebx),%%xmm0			\n\t	movaps	    (%%edx),%%xmm4\n\t"\
		"movaps	0x10(%%ebx),%%xmm1			\n\t	movaps	0x10(%%edx),%%xmm5\n\t"\
		"\n\t"\
		"subpd	%%xmm0,%%xmm2				\n\t	subpd	%%xmm4,%%xmm6\n\t"\
		"subpd	%%xmm1,%%xmm3				\n\t	subpd	%%xmm5,%%xmm7\n\t"\
		"addpd	%%xmm0,%%xmm0				\n\t	addpd	%%xmm4,%%xmm4\n\t"\
		"addpd	%%xmm1,%%xmm1				\n\t	addpd	%%xmm5,%%xmm5\n\t"\
		"addpd	%%xmm2,%%xmm0				\n\t	addpd	%%xmm6,%%xmm4\n\t"\
		"addpd	%%xmm3,%%xmm1				\n\t	addpd	%%xmm7,%%xmm5\n\t"\
		"\n\t"\
		"movl	%[__r1]	,%%esi\n\t"\
		"subpd	%%xmm4	,%%xmm0				\n\t	subpd	%%xmm7,%%xmm2\n\t"\
		"subpd	%%xmm5	,%%xmm1				\n\t	subpd	%%xmm6,%%xmm3\n\t"\
		"movaps	%%xmm0	,0x040(%%esi)		\n\t	movaps	%%xmm2,0x060(%%esi)\n\t"\
		"movaps	%%xmm1	,0x050(%%esi)		\n\t	movaps	%%xmm3,0x030(%%esi)\n\t"\
		"addpd	%%xmm4	,%%xmm4				\n\t	addpd	%%xmm7,%%xmm7\n\t"\
		"addpd	%%xmm5	,%%xmm5				\n\t	addpd	%%xmm6,%%xmm6\n\t"\
		"addpd	%%xmm0	,%%xmm4				\n\t	addpd	%%xmm2,%%xmm7\n\t"\
		"addpd	%%xmm1	,%%xmm5				\n\t	addpd	%%xmm3,%%xmm6\n\t"\
		"movaps	%%xmm4	,     (%%esi)		\n\t	movaps	%%xmm7,0x020(%%esi)\n\t"\
		"movaps	%%xmm5	,0x010(%%esi)		\n\t	movaps	%%xmm6,0x070(%%esi)\n\t"\
		"\n\t"\
		"addl	%%edi	,%%eax\n\t"\
		"addl	%%edi	,%%ebx\n\t"\
		"addl	%%edi	,%%ecx\n\t"\
		"addl	%%edi	,%%edx\n\t"\
		"/* SSE2_RADIX4_DIT_0TWIDDLE_B(r9) */\n\t"\
		"movaps	    (%%eax),%%xmm2			\n\t	movaps	    (%%ecx),%%xmm6\n\t"\
		"movaps	0x10(%%eax),%%xmm3			\n\t	movaps	0x10(%%ecx),%%xmm7\n\t"\
		"movaps	    (%%ebx),%%xmm0			\n\t	movaps	    (%%edx),%%xmm4\n\t"\
		"movaps	0x10(%%ebx),%%xmm1			\n\t	movaps	0x10(%%edx),%%xmm5\n\t"\
		"\n\t"\
		"subpd	%%xmm0,%%xmm2				\n\t	subpd	%%xmm4,%%xmm6\n\t"\
		"subpd	%%xmm1,%%xmm3				\n\t	subpd	%%xmm5,%%xmm7\n\t"\
		"addpd	%%xmm0,%%xmm0				\n\t	addpd	%%xmm4,%%xmm4\n\t"\
		"addpd	%%xmm1,%%xmm1				\n\t	addpd	%%xmm5,%%xmm5\n\t"\
		"addpd	%%xmm2,%%xmm0				\n\t	addpd	%%xmm6,%%xmm4\n\t"\
		"addpd	%%xmm3,%%xmm1				\n\t	addpd	%%xmm7,%%xmm5\n\t"\
		"\n\t"\
		"movl	%[__r9]	,%%esi\n\t"\
		"subpd	%%xmm4	,%%xmm0				\n\t	subpd	%%xmm7,%%xmm2\n\t"\
		"subpd	%%xmm5	,%%xmm1				\n\t	subpd	%%xmm6,%%xmm3\n\t"\
		"movaps	%%xmm0	,0x040(%%esi)		\n\t	movaps	%%xmm2,0x060(%%esi)\n\t"\
		"movaps	%%xmm1	,0x050(%%esi)		\n\t	movaps	%%xmm3,0x030(%%esi)\n\t"\
		"addpd	%%xmm4	,%%xmm4				\n\t	addpd	%%xmm7,%%xmm7\n\t"\
		"addpd	%%xmm5	,%%xmm5				\n\t	addpd	%%xmm6,%%xmm6\n\t"\
		"addpd	%%xmm0	,%%xmm4				\n\t	addpd	%%xmm2,%%xmm7\n\t"\
		"addpd	%%xmm1	,%%xmm5				\n\t	addpd	%%xmm3,%%xmm6\n\t"\
		"movaps	%%xmm4	,     (%%esi)		\n\t	movaps	%%xmm7,0x020(%%esi)\n\t"\
		"movaps	%%xmm5	,0x010(%%esi)		\n\t	movaps	%%xmm6,0x070(%%esi)\n\t"\
		"\n\t"\
		"addl	%%edi	,%%eax\n\t"\
		"addl	%%edi	,%%ebx\n\t"\
		"addl	%%edi	,%%ecx\n\t"\
		"addl	%%edi	,%%edx\n\t"\
		"/* SSE2_RADIX4_DIT_0TWIDDLE_B(r17) */\n\t"\
		"movaps	    (%%eax),%%xmm2			\n\t	movaps	    (%%ecx),%%xmm6\n\t"\
		"movaps	0x10(%%eax),%%xmm3			\n\t	movaps	0x10(%%ecx),%%xmm7\n\t"\
		"movaps	    (%%ebx),%%xmm0			\n\t	movaps	    (%%edx),%%xmm4\n\t"\
		"movaps	0x10(%%ebx),%%xmm1			\n\t	movaps	0x10(%%edx),%%xmm5\n\t"\
		"\n\t"\
		"subpd	%%xmm0,%%xmm2				\n\t	subpd	%%xmm4,%%xmm6\n\t"\
		"subpd	%%xmm1,%%xmm3				\n\t	subpd	%%xmm5,%%xmm7\n\t"\
		"addpd	%%xmm0,%%xmm0				\n\t	addpd	%%xmm4,%%xmm4\n\t"\
		"addpd	%%xmm1,%%xmm1				\n\t	addpd	%%xmm5,%%xmm5\n\t"\
		"addpd	%%xmm2,%%xmm0				\n\t	addpd	%%xmm6,%%xmm4\n\t"\
		"addpd	%%xmm3,%%xmm1				\n\t	addpd	%%xmm7,%%xmm5\n\t"\
		"\n\t"\
		"movl	%[__r17],%%esi\n\t"\
		"subpd	%%xmm4	,%%xmm0				\n\t	subpd	%%xmm7,%%xmm2\n\t"\
		"subpd	%%xmm5	,%%xmm1				\n\t	subpd	%%xmm6,%%xmm3\n\t"\
		"movaps	%%xmm0	,0x040(%%esi)		\n\t	movaps	%%xmm2,0x060(%%esi)\n\t"\
		"movaps	%%xmm1	,0x050(%%esi)		\n\t	movaps	%%xmm3,0x030(%%esi)\n\t"\
		"addpd	%%xmm4	,%%xmm4				\n\t	addpd	%%xmm7,%%xmm7\n\t"\
		"addpd	%%xmm5	,%%xmm5				\n\t	addpd	%%xmm6,%%xmm6\n\t"\
		"addpd	%%xmm0	,%%xmm4				\n\t	addpd	%%xmm2,%%xmm7\n\t"\
		"addpd	%%xmm1	,%%xmm5				\n\t	addpd	%%xmm3,%%xmm6\n\t"\
		"movaps	%%xmm4	,     (%%esi)		\n\t	movaps	%%xmm7,0x020(%%esi)\n\t"\
		"movaps	%%xmm5	,0x010(%%esi)		\n\t	movaps	%%xmm6,0x070(%%esi)\n\t"\
		"\n\t"\
		"addl	%%edi	,%%eax\n\t"\
		"addl	%%edi	,%%ebx\n\t"\
		"addl	%%edi	,%%ecx\n\t"\
		"addl	%%edi	,%%edx\n\t"\
		"/* SSE2_RADIX4_DIT_0TWIDDLE_B(r25) */\n\t"\
		"movaps	    (%%eax),%%xmm2			\n\t	movaps	    (%%ecx),%%xmm6\n\t"\
		"movaps	0x10(%%eax),%%xmm3			\n\t	movaps	0x10(%%ecx),%%xmm7\n\t"\
		"movaps	    (%%ebx),%%xmm0			\n\t	movaps	    (%%edx),%%xmm4\n\t"\
		"movaps	0x10(%%ebx),%%xmm1			\n\t	movaps	0x10(%%edx),%%xmm5\n\t"\
		"\n\t"\
		"subpd	%%xmm0,%%xmm2				\n\t	subpd	%%xmm4,%%xmm6\n\t"\
		"subpd	%%xmm1,%%xmm3				\n\t	subpd	%%xmm5,%%xmm7\n\t"\
		"addpd	%%xmm0,%%xmm0				\n\t	addpd	%%xmm4,%%xmm4\n\t"\
		"addpd	%%xmm1,%%xmm1				\n\t	addpd	%%xmm5,%%xmm5\n\t"\
		"addpd	%%xmm2,%%xmm0				\n\t	addpd	%%xmm6,%%xmm4\n\t"\
		"addpd	%%xmm3,%%xmm1				\n\t	addpd	%%xmm7,%%xmm5\n\t"\
		"\n\t"\
		"movl	%[__r25],%%esi\n\t"\
		"subpd	%%xmm4	,%%xmm0				\n\t	subpd	%%xmm7,%%xmm2\n\t"\
		"subpd	%%xmm5	,%%xmm1				\n\t	subpd	%%xmm6,%%xmm3\n\t"\
		"movaps	%%xmm0	,0x040(%%esi)		\n\t	movaps	%%xmm2,0x060(%%esi)\n\t"\
		"movaps	%%xmm1	,0x050(%%esi)		\n\t	movaps	%%xmm3,0x030(%%esi)\n\t"\
		"addpd	%%xmm4	,%%xmm4				\n\t	addpd	%%xmm7,%%xmm7\n\t"\
		"addpd	%%xmm5	,%%xmm5				\n\t	addpd	%%xmm6,%%xmm6\n\t"\
		"addpd	%%xmm0	,%%xmm4				\n\t	addpd	%%xmm2,%%xmm7\n\t"\
		"addpd	%%xmm1	,%%xmm5				\n\t	addpd	%%xmm3,%%xmm6\n\t"\
		"movaps	%%xmm4	,     (%%esi)		\n\t	movaps	%%xmm7,0x020(%%esi)\n\t"\
		"movaps	%%xmm5	,0x010(%%esi)		\n\t	movaps	%%xmm6,0x070(%%esi)\n\t"\
		"\n\t"\
		"/****************************************************************************************************\n\t"\
		"!...and now do four more radix-4 transforms, including the internal and external twiddle factors:   !\n\t"\
		"****************************************************************************************************/\n\t"\
		"\n\t"\
	"/*...Block 1: 24 MOVapd, 26 ADD/SUBpd, 12 MULpd */\n\t"\
		"\n\t"\
		"movl	%[__add0]	,%%eax\n\t"\
		"movl	%[__p4]	,%%ebx\n\t"\
		"movl	%[__r1]	,%%ecx\n\t"\
		"movl	%[__r9]	,%%edx\n\t"\
		"movl	%[__p8]	,%%edi\n\t"\
		"shll	$3		,%%ebx\n\t"\
		"shll	$3		,%%edi\n\t"\
		"addl	%%eax	,%%ebx\n\t"\
		"\n\t"\
		"movaps	     (%%edx)	,%%xmm2		\n\t	movaps	0x100(%%edx),%%xmm4\n\t"\
		"movaps	0x010(%%edx)	,%%xmm3		\n\t	movaps	0x110(%%edx),%%xmm5\n\t"\
		"movaps	     (%%ecx)	,%%xmm0		\n\t	movaps	0x100(%%ecx),%%xmm6\n\t"\
		"movaps	0x010(%%ecx)	,%%xmm1		\n\t	movaps	0x110(%%ecx),%%xmm7\n\t"\
		"\n\t"\
		"subpd	%%xmm2,%%xmm0				\n\t	subpd	%%xmm4,%%xmm6\n\t"\
		"subpd	%%xmm3,%%xmm1				\n\t	subpd	%%xmm5,%%xmm7\n\t"\
		"addpd	%%xmm2,%%xmm2				\n\t	addpd	%%xmm4,%%xmm4\n\t"\
		"addpd	%%xmm3,%%xmm3				\n\t	addpd	%%xmm5,%%xmm5\n\t"\
		"addpd	%%xmm0,%%xmm2				\n\t	addpd	%%xmm6,%%xmm4\n\t"\
		"addpd	%%xmm1,%%xmm3				\n\t	addpd	%%xmm7,%%xmm5\n\t"\
		"\n\t"\
		"/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c0) */\n\t"\
		"movl	%[__c0], %%esi\n\t"\
		"\n\t"\
		"subpd	%%xmm4,%%xmm2				\n\t	subpd	%%xmm7,%%xmm0\n\t"\
		"subpd	%%xmm5,%%xmm3				\n\t	subpd	%%xmm6,%%xmm1\n\t"\
		"addpd	%%xmm4,%%xmm4				\n\t	addpd	%%xmm7,%%xmm7\n\t"\
		"addpd	%%xmm5,%%xmm5				\n\t	addpd	%%xmm6,%%xmm6\n\t"\
		"addpd	%%xmm2,%%xmm4				\n\t	addpd	%%xmm0,%%xmm7\n\t"\
		"addpd	%%xmm3,%%xmm5				\n\t	addpd	%%xmm1,%%xmm6\n\t"\
		"movaps	%%xmm2,     (%%ecx)			\n\t	movaps	%%xmm0,0x010(%%edx)\n\t"\
		"movaps	%%xmm3,0x010(%%ecx)			\n\t	movaps	%%xmm6,     (%%edx)\n\t"\
		"movaps	%%xmm4,%%xmm2				\n\t	movaps	%%xmm7,%%xmm0\n\t"\
		"movaps	%%xmm5,%%xmm3				\n\t	movaps	%%xmm1,%%xmm6\n\t"\
		"mulpd	    (%%esi),%%xmm4			\n\t	mulpd	0x20(%%esi),%%xmm7\n\t"\
		"mulpd	    (%%esi),%%xmm5			\n\t	mulpd	0x20(%%esi),%%xmm1\n\t"\
		"mulpd	0x10(%%esi),%%xmm2			\n\t	mulpd	0x30(%%esi),%%xmm0\n\t"\
		"mulpd	0x10(%%esi),%%xmm3			\n\t	mulpd	0x30(%%esi),%%xmm6\n\t"\
		"subpd	%%xmm2,%%xmm5				\n\t	subpd	%%xmm0,%%xmm1\n\t"\
		"addpd	%%xmm3,%%xmm4				\n\t	addpd	%%xmm6,%%xmm7\n\t"\
		"\n\t"\
		"movaps	%%xmm5,0x10(%%eax)			\n\t	movaps	%%xmm1,0x10(%%ebx)\n\t"\
		"movaps	%%xmm4,    (%%eax)			\n\t	movaps	%%xmm7,    (%%ebx)\n\t"\
		"\n\t"\
		"addl	%%edi	,%%eax				\n\t	addl	 %%edi,%%ebx\n\t"\
		"addl	$0x40,%%esi	/* c2 */\n\t"\
		"\n\t"\
		"movaps	     (%%ecx),%%xmm4			\n\t	movaps	0x010(%%edx),%%xmm0\n\t"\
		"movaps	0x010(%%ecx),%%xmm5			\n\t	movaps	     (%%edx),%%xmm6\n\t"\
		"movaps	%%xmm4,%%xmm2				\n\t	movaps	%%xmm0,%%xmm1\n\t"\
		"movaps	%%xmm5,%%xmm3				\n\t	movaps	%%xmm6,%%xmm7\n\t"\
		"mulpd	    (%%esi),%%xmm4			\n\t	mulpd	0x20(%%esi),%%xmm0\n\t"\
		"mulpd	    (%%esi),%%xmm5			\n\t	mulpd	0x20(%%esi),%%xmm6\n\t"\
		"mulpd	0x10(%%esi),%%xmm2			\n\t	mulpd	0x30(%%esi),%%xmm1\n\t"\
		"mulpd	0x10(%%esi),%%xmm3			\n\t	mulpd	0x30(%%esi),%%xmm7\n\t"\
		"subpd	%%xmm2,%%xmm5				\n\t	subpd	%%xmm1,%%xmm6\n\t"\
		"addpd	%%xmm3,%%xmm4				\n\t	addpd	%%xmm7,%%xmm0\n\t"\
		"\n\t"\
		"movaps	%%xmm5,0x10(%%eax)			\n\t	movaps	%%xmm6,0x10(%%ebx)\n\t"\
		"movaps	%%xmm4,    (%%eax)			\n\t	movaps	%%xmm0,    (%%ebx)\n\t"\
		"\n\t"\
	"/*...Block 2: 34 MOVapd, 36 ADD/SUBpd, 26 MULpd */\n\t"\
		"\n\t"\
		"movl	%[__add0]	,%%eax\n\t"\
		"movl	%[__p1]	,%%esi\n\t"\
		"movl	%[__p4]	,%%ebx\n\t"\
		"shll	$3		,%%esi\n\t"\
		"shll	$3		,%%ebx\n\t"\
		"addl	%%esi		,%%eax\n\t"\
		"addl	%%eax		,%%ebx\n\t"\
		"movl	%[__r3]	,%%ecx\n\t"\
		"movl	%[__r11]	,%%edx\n\t"\
		"addl	$0x100	,%%ecx\n\t"\
		"addl	$0x100	,%%edx\n\t"\
		"movl	%[__cc0]	,%%esi\n\t"\
		"\n\t"\
		"movaps	     (%%ecx)	,%%xmm4\n\t"\
		"movaps	0x010(%%ecx)	,%%xmm5\n\t"\
		"movaps	     (%%esi)	,%%xmm2\n\t"\
		"movaps	0x010(%%esi)	,%%xmm3\n\t"\
		"movaps	%%xmm4,%%xmm6\n\t"\
		"movaps	%%xmm5,%%xmm7\n\t"\
		"\n\t"\
		"mulpd	%%xmm2,%%xmm4\n\t"\
		"mulpd	%%xmm2,%%xmm5\n\t"\
		"mulpd	%%xmm3,%%xmm6				\n\t	movaps	     (%%edx),%%xmm0\n\t"\
		"mulpd	%%xmm3,%%xmm7				\n\t	movaps	0x010(%%edx),%%xmm1\n\t"\
		"subpd	%%xmm6,%%xmm5				\n\t	movaps	%%xmm0,%%xmm6\n\t"\
		"addpd	%%xmm7,%%xmm4				\n\t	movaps	%%xmm1,%%xmm7\n\t"\
		"\n\t"\
		"									\n\t	mulpd	%%xmm3,%%xmm0\n\t"\
		"									\n\t	mulpd	%%xmm3,%%xmm1\n\t"\
		"									\n\t	mulpd	%%xmm2,%%xmm6\n\t"\
		"									\n\t	mulpd	%%xmm2,%%xmm7\n\t"\
		"									\n\t	subpd	%%xmm6,%%xmm1\n\t"\
		"									\n\t	addpd	%%xmm7,%%xmm0\n\t"\
		"\n\t"\
		"movaps	%%xmm5,%%xmm7\n\t"\
		"movaps	%%xmm4,%%xmm6\n\t"\
		"\n\t"\
		"addpd	%%xmm0,%%xmm4\n\t"\
		"addpd	%%xmm1,%%xmm5\n\t"\
		"subpd	%%xmm0,%%xmm6\n\t"\
		"subpd	%%xmm1,%%xmm7\n\t"\
		"\n\t"\
		"subl	$0x100	,%%ecx\n\t"\
		"subl	$0x100	,%%edx\n\t"\
		"movl	%[__isrt2],%%esi\n\t"\
		"movaps	     (%%edx),%%xmm2\n\t"\
		"movaps	0x010(%%edx),%%xmm3\n\t"\
		"movaps	(%%esi),%%xmm1\n\t"\
		"movaps	%%xmm3,%%xmm0\n\t"\
		"subpd	%%xmm2,%%xmm3\n\t"\
		"addpd	%%xmm0,%%xmm2\n\t"\
		"mulpd	%%xmm1,%%xmm2\n\t"\
		"mulpd	%%xmm1,%%xmm3\n\t"\
		"\n\t"\
		"movaps	     (%%ecx),%%xmm0\n\t"\
		"movaps	0x010(%%ecx),%%xmm1\n\t"\
		"subpd	%%xmm2,%%xmm0\n\t"\
		"subpd	%%xmm3,%%xmm1\n\t"\
		"addpd	%%xmm2,%%xmm2\n\t"\
		"addpd	%%xmm3,%%xmm3\n\t"\
		"addpd	%%xmm0,%%xmm2\n\t"\
		"addpd	%%xmm1,%%xmm3\n\t"\
		"\n\t"\
		"/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c1) */\n\t"\
		"movl	%[__c1], %%esi\n\t"\
		"\n\t"\
		"subpd	%%xmm4,%%xmm2				\n\t	subpd	%%xmm7,%%xmm0\n\t"\
		"subpd	%%xmm5,%%xmm3				\n\t	subpd	%%xmm6,%%xmm1\n\t"\
		"addpd	%%xmm4,%%xmm4				\n\t	addpd	%%xmm7,%%xmm7\n\t"\
		"addpd	%%xmm5,%%xmm5				\n\t	addpd	%%xmm6,%%xmm6\n\t"\
		"addpd	%%xmm2,%%xmm4				\n\t	addpd	%%xmm0,%%xmm7\n\t"\
		"addpd	%%xmm3,%%xmm5				\n\t	addpd	%%xmm1,%%xmm6\n\t"\
		"movaps	%%xmm2,     (%%ecx)			\n\t	movaps	%%xmm0,0x010(%%edx)\n\t"\
		"movaps	%%xmm3,0x010(%%ecx)			\n\t	movaps	%%xmm6,     (%%edx)\n\t"\
		"movaps	%%xmm4,%%xmm2				\n\t	movaps	%%xmm7,%%xmm0\n\t"\
		"movaps	%%xmm5,%%xmm3				\n\t	movaps	%%xmm1,%%xmm6\n\t"\
		"mulpd	    (%%esi),%%xmm4			\n\t	mulpd	0x20(%%esi),%%xmm7\n\t"\
		"mulpd	    (%%esi),%%xmm5			\n\t	mulpd	0x20(%%esi),%%xmm1\n\t"\
		"mulpd	0x10(%%esi),%%xmm2			\n\t	mulpd	0x30(%%esi),%%xmm0\n\t"\
		"mulpd	0x10(%%esi),%%xmm3			\n\t	mulpd	0x30(%%esi),%%xmm6\n\t"\
		"subpd	%%xmm2,%%xmm5				\n\t	subpd	%%xmm0,%%xmm1\n\t"\
		"addpd	%%xmm3,%%xmm4				\n\t	addpd	%%xmm6,%%xmm7\n\t"\
		"\n\t"\
		"movaps	%%xmm5,0x10(%%eax)			\n\t	movaps	%%xmm1,0x10(%%ebx)\n\t"\
		"movaps	%%xmm4,    (%%eax)			\n\t	movaps	%%xmm7,    (%%ebx)\n\t"\
		"\n\t"\
		"addl	%%edi	,%%eax				\n\t	addl	 %%edi,%%ebx\n\t"\
		"addl	$0x40,%%esi	/* c2 */\n\t"\
		"\n\t"\
		"movaps	     (%%ecx),%%xmm4			\n\t	movaps	0x010(%%edx),%%xmm0\n\t"\
		"movaps	0x010(%%ecx),%%xmm5			\n\t	movaps	     (%%edx),%%xmm6\n\t"\
		"movaps	%%xmm4,%%xmm2				\n\t	movaps	%%xmm0,%%xmm1\n\t"\
		"movaps	%%xmm5,%%xmm3				\n\t	movaps	%%xmm6,%%xmm7\n\t"\
		"mulpd	    (%%esi),%%xmm4			\n\t	mulpd	0x20(%%esi),%%xmm0\n\t"\
		"mulpd	    (%%esi),%%xmm5			\n\t	mulpd	0x20(%%esi),%%xmm6\n\t"\
		"mulpd	0x10(%%esi),%%xmm2			\n\t	mulpd	0x30(%%esi),%%xmm1\n\t"\
		"mulpd	0x10(%%esi),%%xmm3			\n\t	mulpd	0x30(%%esi),%%xmm7\n\t"\
		"subpd	%%xmm2,%%xmm5				\n\t	subpd	%%xmm1,%%xmm6\n\t"\
		"addpd	%%xmm3,%%xmm4				\n\t	addpd	%%xmm7,%%xmm0\n\t"\
		"\n\t"\
		"movaps	%%xmm5,0x10(%%eax)			\n\t	movaps	%%xmm6,0x10(%%ebx)\n\t"\
		"movaps	%%xmm4,    (%%eax)			\n\t	movaps	%%xmm0,    (%%ebx)\n\t"\
		"\n\t"\
	"/*...Block 3: 31 MOVapd, 32 ADD/SUBpd, 20 MULpd */\n\t"\
		"\n\t"\
		"movl	%[__add0]	,%%eax\n\t"\
		"movl	%[__p2]	,%%esi\n\t"\
		"movl	%[__p4]	,%%ebx\n\t"\
		"shll	$3	,%%esi\n\t"\
		"shll	$3	,%%ebx\n\t"\
		"addl	%%esi	,%%eax\n\t"\
		"addl	%%eax	,%%ebx\n\t"\
		"movl	%[__r5]	,%%ecx\n\t"\
		"movl	%[__r13]	,%%edx\n\t"\
		"addl	$0x100	,%%ecx\n\t"\
		"addl	$0x100	,%%edx\n\t"\
		"movl	%[__isrt2],%%esi\n\t"\
		"movaps	(%%esi)	,%%xmm2	/* isrt2 */\n\t"\
		"\n\t"\
		"movaps	     (%%ecx),%%xmm4\n\t"\
		"movaps	0x010(%%ecx),%%xmm5\n\t"\
		"movaps	     (%%edx),%%xmm6\n\t"\
		"movaps	0x010(%%edx),%%xmm7\n\t"\
		"subl	$0x100,%%ecx\n\t"\
		"subl	$0x100,%%edx\n\t"\
		"mulpd	%%xmm2,%%xmm4\n\t"\
		"mulpd	%%xmm2,%%xmm5\n\t"\
		"mulpd	%%xmm2,%%xmm6\n\t"\
		"mulpd	%%xmm2,%%xmm7\n\t"\
		"\n\t"\
		"subpd	%%xmm4,%%xmm5				\n\t	movaps	     (%%ecx),%%xmm0\n\t"\
		"subpd	%%xmm7,%%xmm6				\n\t	movaps	0x010(%%edx),%%xmm2\n\t"\
		"addpd	%%xmm4,%%xmm4				\n\t	movaps	0x010(%%ecx),%%xmm3\n\t"\
		"addpd	%%xmm7,%%xmm7				\n\t	movaps	     (%%edx),%%xmm1\n\t"\
		"addpd	%%xmm5,%%xmm4\n\t"\
		"addpd	%%xmm6,%%xmm7\n\t"\
		"\n\t"\
		"subpd	%%xmm6,%%xmm4				\n\t	subpd	%%xmm2,%%xmm0\n\t"\
		"subpd	%%xmm7,%%xmm5				\n\t	subpd	%%xmm1,%%xmm3\n\t"\
		"addpd	%%xmm6,%%xmm6				\n\t	addpd	%%xmm2,%%xmm2\n\t"\
		"addpd	%%xmm7,%%xmm7				\n\t	addpd	%%xmm1,%%xmm1\n\t"\
		"addpd	%%xmm4,%%xmm6				\n\t	addpd	%%xmm0,%%xmm2\n\t"\
		"addpd	%%xmm5,%%xmm7				\n\t	addpd	%%xmm3,%%xmm1\n\t"\
		"\n\t"\
		"/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c2) */\n\t"\
		"movl	%[__c2], %%esi\n\t"\
		"\n\t"\
		"subpd	%%xmm4,%%xmm2				\n\t	subpd	%%xmm7,%%xmm0\n\t"\
		"subpd	%%xmm5,%%xmm3				\n\t	subpd	%%xmm6,%%xmm1\n\t"\
		"addpd	%%xmm4,%%xmm4				\n\t	addpd	%%xmm7,%%xmm7\n\t"\
		"addpd	%%xmm5,%%xmm5				\n\t	addpd	%%xmm6,%%xmm6\n\t"\
		"addpd	%%xmm2,%%xmm4				\n\t	addpd	%%xmm0,%%xmm7\n\t"\
		"addpd	%%xmm3,%%xmm5				\n\t	addpd	%%xmm1,%%xmm6\n\t"\
		"movaps	%%xmm2,     (%%ecx)			\n\t	movaps	%%xmm0,0x010(%%edx)\n\t"\
		"movaps	%%xmm3,0x010(%%ecx)			\n\t	movaps	%%xmm6,     (%%edx)\n\t"\
		"movaps	%%xmm4,%%xmm2				\n\t	movaps	%%xmm7,%%xmm0\n\t"\
		"movaps	%%xmm5,%%xmm3				\n\t	movaps	%%xmm1,%%xmm6\n\t"\
		"mulpd	    (%%esi),%%xmm4			\n\t	mulpd	0x20(%%esi),%%xmm7\n\t"\
		"mulpd	    (%%esi),%%xmm5			\n\t	mulpd	0x20(%%esi),%%xmm1\n\t"\
		"mulpd	0x10(%%esi),%%xmm2			\n\t	mulpd	0x30(%%esi),%%xmm0\n\t"\
		"mulpd	0x10(%%esi),%%xmm3			\n\t	mulpd	0x30(%%esi),%%xmm6\n\t"\
		"subpd	%%xmm2,%%xmm5				\n\t	subpd	%%xmm0,%%xmm1\n\t"\
		"addpd	%%xmm3,%%xmm4				\n\t	addpd	%%xmm6,%%xmm7\n\t"\
		"\n\t"\
		"movaps	%%xmm5,0x10(%%eax)			\n\t	movaps	%%xmm1,0x10(%%ebx)\n\t"\
		"movaps	%%xmm4,    (%%eax)			\n\t	movaps	%%xmm7,    (%%ebx)\n\t"\
		"\n\t"\
		"addl	%%edi	,%%eax				\n\t	addl	 %%edi,%%ebx\n\t"\
		"addl	$0x40,%%esi	/* c2 */\n\t"\
		"\n\t"\
		"movaps	     (%%ecx),%%xmm4			\n\t	movaps	0x010(%%edx),%%xmm0\n\t"\
		"movaps	0x010(%%ecx),%%xmm5			\n\t	movaps	     (%%edx),%%xmm6\n\t"\
		"movaps	%%xmm4,%%xmm2				\n\t	movaps	%%xmm0,%%xmm1\n\t"\
		"movaps	%%xmm5,%%xmm3				\n\t	movaps	%%xmm6,%%xmm7\n\t"\
		"mulpd	    (%%esi),%%xmm4			\n\t	mulpd	0x20(%%esi),%%xmm0\n\t"\
		"mulpd	    (%%esi),%%xmm5			\n\t	mulpd	0x20(%%esi),%%xmm6\n\t"\
		"mulpd	0x10(%%esi),%%xmm2			\n\t	mulpd	0x30(%%esi),%%xmm1\n\t"\
		"mulpd	0x10(%%esi),%%xmm3			\n\t	mulpd	0x30(%%esi),%%xmm7\n\t"\
		"subpd	%%xmm2,%%xmm5				\n\t	subpd	%%xmm1,%%xmm6\n\t"\
		"addpd	%%xmm3,%%xmm4				\n\t	addpd	%%xmm7,%%xmm0\n\t"\
		"\n\t"\
		"movaps	%%xmm5,0x10(%%eax)			\n\t	movaps	%%xmm6,0x10(%%ebx)\n\t"\
		"movaps	%%xmm4,    (%%eax)			\n\t	movaps	%%xmm0,    (%%ebx)\n\t"\
		"\n\t"\
	"/*...Block 4: 34 MOVapd, 36 ADD/SUBpd, 26 MULpd */\n\t"\
		"\n\t"\
		"movl	%[__add0]	,%%eax\n\t"\
		"movl	%[__p3]	,%%esi\n\t"\
		"movl	%[__p4]	,%%ebx\n\t"\
		"shll	$3	,%%esi\n\t"\
		"shll	$3	,%%ebx\n\t"\
		"addl	%%esi	,%%eax\n\t"\
		"addl	%%eax	,%%ebx\n\t"\
		"movl	%[__r7]	,%%ecx\n\t"\
		"movl	%[__r15]	,%%edx\n\t"\
		"addl	$0x100	,%%ecx\n\t"\
		"addl	$0x100	,%%edx\n\t"\
		"movl	%[__cc0]	,%%esi\n\t"\
		"\n\t"\
		"movaps	     (%%ecx),%%xmm4\n\t"\
		"movaps	0x010(%%ecx),%%xmm5\n\t"\
		"movaps	     (%%esi),%%xmm3\n\t"\
		"movaps	0x010(%%esi),%%xmm2\n\t"\
		"movaps	%%xmm4,%%xmm6\n\t"\
		"movaps	%%xmm5,%%xmm7\n\t"\
		"\n\t"\
		"mulpd	%%xmm2,%%xmm4\n\t"\
		"mulpd	%%xmm2,%%xmm5\n\t"\
		"mulpd	%%xmm3,%%xmm6				\n\t	movaps	     (%%edx),%%xmm0\n\t"\
		"mulpd	%%xmm3,%%xmm7				\n\t	movaps	0x010(%%edx),%%xmm1\n\t"\
		"subpd	%%xmm6,%%xmm5				\n\t	movaps	%%xmm0,%%xmm6\n\t"\
		"addpd	%%xmm7,%%xmm4				\n\t	movaps	%%xmm1,%%xmm7\n\t"\
		"\n\t"\
		"									\n\t	mulpd	%%xmm3,%%xmm0\n\t"\
		"									\n\t	mulpd	%%xmm3,%%xmm1\n\t"\
		"									\n\t	mulpd	%%xmm2,%%xmm6\n\t"\
		"									\n\t	mulpd	%%xmm2,%%xmm7\n\t"\
		"									\n\t	subpd	%%xmm6,%%xmm1\n\t"\
		"									\n\t	addpd	%%xmm7,%%xmm0\n\t"\
		"\n\t"\
		"movaps	%%xmm5,%%xmm7\n\t"\
		"movaps	%%xmm4,%%xmm6\n\t"\
		"\n\t"\
		"addpd	%%xmm0,%%xmm6\n\t"\
		"addpd	%%xmm1,%%xmm7\n\t"\
		"subpd	%%xmm0,%%xmm4\n\t"\
		"subpd	%%xmm1,%%xmm5\n\t"\
		"\n\t"\
		"subl	$0x100	,%%ecx\n\t"\
		"subl	$0x100	,%%edx\n\t"\
		"movl	%[__isrt2],%%esi\n\t"\
		"movaps	     (%%edx)	,%%xmm0\n\t"\
		"movaps	0x010(%%edx)	,%%xmm1\n\t"\
		"movaps	(%%esi)	,%%xmm3\n\t"\
		"movaps	%%xmm0	,%%xmm2\n\t"\
		"subpd	%%xmm1	,%%xmm0\n\t"\
		"addpd	%%xmm2	,%%xmm1\n\t"\
		"mulpd	%%xmm3	,%%xmm0\n\t"\
		"mulpd	%%xmm3	,%%xmm1\n\t"\
		"\n\t"\
		"movaps	     (%%ecx)	,%%xmm2\n\t"\
		"movaps	0x010(%%ecx)	,%%xmm3\n\t"\
		"subpd	%%xmm0	,%%xmm2\n\t"\
		"subpd	%%xmm1	,%%xmm3\n\t"\
		"addpd	%%xmm0	,%%xmm0\n\t"\
		"addpd	%%xmm1	,%%xmm1\n\t"\
		"addpd	%%xmm2	,%%xmm0\n\t"\
		"addpd	%%xmm3	,%%xmm1\n\t"\
		"\n\t"\
		"/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c3) */\n\t"\
		"movl	%[__c3], %%esi\n\t"\
		"\n\t"\
		"subpd	%%xmm4,%%xmm2				\n\t	subpd	%%xmm7,%%xmm0\n\t"\
		"subpd	%%xmm5,%%xmm3				\n\t	subpd	%%xmm6,%%xmm1\n\t"\
		"addpd	%%xmm4,%%xmm4				\n\t	addpd	%%xmm7,%%xmm7\n\t"\
		"addpd	%%xmm5,%%xmm5				\n\t	addpd	%%xmm6,%%xmm6\n\t"\
		"addpd	%%xmm2,%%xmm4				\n\t	addpd	%%xmm0,%%xmm7\n\t"\
		"addpd	%%xmm3,%%xmm5				\n\t	addpd	%%xmm1,%%xmm6\n\t"\
		"movaps	%%xmm2,     (%%ecx)			\n\t	movaps	%%xmm0,0x010(%%edx)\n\t"\
		"movaps	%%xmm3,0x010(%%ecx)			\n\t	movaps	%%xmm6,     (%%edx)\n\t"\
		"movaps	%%xmm4,%%xmm2				\n\t	movaps	%%xmm7,%%xmm0\n\t"\
		"movaps	%%xmm5,%%xmm3				\n\t	movaps	%%xmm1,%%xmm6\n\t"\
		"mulpd	    (%%esi),%%xmm4			\n\t	mulpd	0x20(%%esi),%%xmm7\n\t"\
		"mulpd	    (%%esi),%%xmm5			\n\t	mulpd	0x20(%%esi),%%xmm1\n\t"\
		"mulpd	0x10(%%esi),%%xmm2			\n\t	mulpd	0x30(%%esi),%%xmm0\n\t"\
		"mulpd	0x10(%%esi),%%xmm3			\n\t	mulpd	0x30(%%esi),%%xmm6\n\t"\
		"subpd	%%xmm2,%%xmm5				\n\t	subpd	%%xmm0,%%xmm1\n\t"\
		"addpd	%%xmm3,%%xmm4				\n\t	addpd	%%xmm6,%%xmm7\n\t"\
		"\n\t"\
		"movaps	%%xmm5,0x10(%%eax)			\n\t	movaps	%%xmm1,0x10(%%ebx)\n\t"\
		"movaps	%%xmm4,    (%%eax)			\n\t	movaps	%%xmm7,    (%%ebx)\n\t"\
		"\n\t"\
		"addl	%%edi	,%%eax				\n\t	addl	 %%edi,%%ebx\n\t"\
		"addl	$0x40,%%esi	/* c2 */\n\t"\
		"\n\t"\
		"movaps	     (%%ecx),%%xmm4			\n\t	movaps	0x010(%%edx),%%xmm0\n\t"\
		"movaps	0x010(%%ecx),%%xmm5			\n\t	movaps	     (%%edx),%%xmm6\n\t"\
		"movaps	%%xmm4,%%xmm2				\n\t	movaps	%%xmm0,%%xmm1\n\t"\
		"movaps	%%xmm5,%%xmm3				\n\t	movaps	%%xmm6,%%xmm7\n\t"\
		"mulpd	    (%%esi),%%xmm4			\n\t	mulpd	0x20(%%esi),%%xmm0\n\t"\
		"mulpd	    (%%esi),%%xmm5			\n\t	mulpd	0x20(%%esi),%%xmm6\n\t"\
		"mulpd	0x10(%%esi),%%xmm2			\n\t	mulpd	0x30(%%esi),%%xmm1\n\t"\
		"mulpd	0x10(%%esi),%%xmm3			\n\t	mulpd	0x30(%%esi),%%xmm7\n\t"\
		"subpd	%%xmm2,%%xmm5				\n\t	subpd	%%xmm1,%%xmm6\n\t"\
		"addpd	%%xmm3,%%xmm4				\n\t	addpd	%%xmm7,%%xmm0\n\t"\
		"\n\t"\
		"movaps	%%xmm5,0x10(%%eax)			\n\t	movaps	%%xmm6,0x10(%%ebx)\n\t"\
		"movaps	%%xmm4,    (%%eax)			\n\t	movaps	%%xmm0,    (%%ebx)\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__p1] "m" (Xp1)\
		 ,[__p2] "m" (Xp2)\
		 ,[__p3] "m" (Xp3)\
		 ,[__p4] "m" (Xp4)\
		 ,[__p8] "m" (Xp8)\
		 ,[__r1] "m" (Xr1)\
		 ,[__r3] "m" (Xr3)\
		 ,[__r5] "m" (Xr5)\
		 ,[__r7] "m" (Xr7)\
		 ,[__r9] "m" (Xr9)\
		 ,[__r11] "m" (Xr11)\
		 ,[__r13] "m" (Xr13)\
		 ,[__r15] "m" (Xr15)\
		 ,[__r17] "m" (Xr17)\
		 ,[__r25] "m" (Xr25)\
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__cc0] "m" (Xcc0)\
		 ,[__c0] "m" (Xc0)\
		 ,[__c1] "m" (Xc1)\
		 ,[__c2] "m" (Xc2)\
		 ,[__c3] "m" (Xc3)\
		: "eax","ebx","ecx","edx","edi","esi"	/* Clobbered registers */\
	);\
	}


	#define SSE2_RADIX16_DIF_TWIDDLE(Xadd0,Xp1,Xp2,Xp3,Xp4,Xp8,Xp12,Xr1,Xr3,Xr5,Xr7,Xr9,Xr17,Xr25,Xisrt2,Xcc0,Xc0,Xc1,Xc2,Xc3)\
	{\
	__asm__ volatile (\
		"movl	%[__add0]	,%%eax			\n\t"\
		"movl	%[__p4]		,%%ebx			\n\t"\
		"movl	%[__p8]		,%%ecx			\n\t"\
		"movl	%[__p12]	,%%edx			\n\t"\
		"shll	$3			,%%ebx			\n\t"\
		"shll	$3			,%%ecx			\n\t"\
		"shll	$3			,%%edx			\n\t"\
		"addl	%%eax		,%%ebx			\n\t"\
		"addl	%%eax		,%%ecx			\n\t"\
		"addl	%%eax		,%%edx			\n\t"\
		"/* SSE2_RADIX4_DIF_4TWIDDLE_B(r1 ,c0) */\n\t"\
		"/* Do the p0,p8 combo: */			\n\t"\
		"movl	%[__c0]		,%%esi 			\n\t"\
		"movaps	    (%%eax)	,%%xmm0			\n\t		movaps	    (%%ecx)	,%%xmm4			\n\t"\
		"movaps	0x10(%%eax)	,%%xmm1			\n\t		movaps	0x10(%%ecx)	,%%xmm5			\n\t"\
		"movaps	    (%%esi)	,%%xmm6			\n\t"\
		"movaps	0x10(%%esi)	,%%xmm7			\n\t"\
		"movaps	%%xmm0		,%%xmm2			\n\t"\
		"movaps	%%xmm1		,%%xmm3			\n\t"\
		"mulpd   %%xmm6		,%%xmm0			\n\t"\
		"mulpd   %%xmm6		,%%xmm1			\n\t"\
		"mulpd   %%xmm7		,%%xmm2			\n\t"\
		"mulpd   %%xmm7		,%%xmm3			\n\t"\
		"												movaps	%%xmm4		,%%xmm6			\n\t"\
		"addpd   %%xmm2		,%%xmm1			\n\t		movaps	%%xmm5		,%%xmm7			\n\t"\
		"												mulpd   0x20(%%esi) ,%%xmm4			\n\t"\
		"subpd   %%xmm3		,%%xmm0			\n\t		mulpd   0x20(%%esi)	,%%xmm5			\n\t"\
		"												mulpd   0x30(%%esi)	,%%xmm6			\n\t"\
		"movaps	%%xmm0		,%%xmm2			\n\t		mulpd	0x30(%%esi)	,%%xmm7			\n\t"\
		"												addpd   %%xmm6	    ,%%xmm5			\n\t"\
		"movaps	%%xmm1		,%%xmm3			\n\t		subpd	%%xmm7		,%%xmm4			\n\t"\
		"addpd	%%xmm4		,%%xmm0			\n\t"\
		"addpd	%%xmm5		,%%xmm1			\n\t"\
		"subpd	%%xmm4		,%%xmm2			\n\t"\
		"subpd	%%xmm5		,%%xmm3			\n\t"\
		"/* Do the p4,12 combo: */			\n\t"\
		"addl	$0x60		,%%esi 			\n\t"\
		"movaps	    (%%edx)	,%%xmm4			\n\t"\
		"movaps	0x10(%%edx)	,%%xmm5			\n\t"\
		"movaps	    (%%edx)	,%%xmm6			\n\t"\
		"movaps	0x10(%%edx)	,%%xmm7			\n\t"\
		"mulpd	    (%%esi)	,%%xmm4			\n\t"\
		"mulpd	    (%%esi)	,%%xmm5			\n\t"\
		"mulpd	0x10(%%esi)	,%%xmm6			\n\t"\
		"mulpd	0x10(%%esi)	,%%xmm7			\n\t"\
		"movl	%[__r1]		,%%esi 			\n\t"\
		"addpd	%%xmm6		,%%xmm5			\n\t"\
		"subpd	%%xmm7		,%%xmm4			\n\t"\
		"movaps	%%xmm5		,0x010(%%esi)	\n\t"\
		"movaps	%%xmm4		,     (%%esi)	\n\t"\
		"movl	%[__c0]		,%%esi 			\n\t"\
		"addl	$0x40		,%%esi 			\n\t"\
		"movaps	    (%%ebx)	,%%xmm4			\n\t"\
		"movaps	0x10(%%ebx)	,%%xmm5			\n\t"\
		"movaps	    (%%ebx)	,%%xmm6			\n\t"\
		"movaps	0x10(%%ebx)	,%%xmm7			\n\t"\
		"mulpd	    (%%esi)	,%%xmm4			\n\t"\
		"mulpd	    (%%esi)	,%%xmm5			\n\t"\
		"mulpd	0x10(%%esi)	,%%xmm6			\n\t"\
		"mulpd	0x10(%%esi)	,%%xmm7			\n\t"\
		"movl	%[__r1]		,%%esi 			\n\t"\
		"addpd	%%xmm6		,%%xmm5			\n\t"\
		"subpd	%%xmm7		,%%xmm4			\n\t"\
		"movaps	%%xmm5		,%%xmm7			\n\t"\
		"movaps	%%xmm4		,%%xmm6			\n\t"\
		"subpd	     (%%esi),%%xmm4			\n\t"\
		"subpd	0x010(%%esi),%%xmm5			\n\t"\
		"addpd	     (%%esi),%%xmm6			\n\t"\
		"addpd	0x010(%%esi),%%xmm7			\n\t"\
		"/* Finish radix-4 butterfly and store results into temporary-array slots: */\n\t"\
		"subpd	%%xmm6,%%xmm0				\n\t			subpd	%%xmm5	,%%xmm2				\n\t"\
		"subpd	%%xmm7,%%xmm1				\n\t			subpd	%%xmm4	,%%xmm3				\n\t"\
		"movaps	%%xmm0,0x040(%%esi)			\n\t			movaps	%%xmm2	,0x020(%%esi)		\n\t"\
		"movaps	%%xmm1,0x050(%%esi)			\n\t			movaps	%%xmm3	,0x070(%%esi)		\n\t"\
		"addpd	%%xmm6,%%xmm6				\n\t			addpd	%%xmm5	,%%xmm5				\n\t"\
		"addpd	%%xmm7,%%xmm7				\n\t			addpd	%%xmm4	,%%xmm4				\n\t"\
		"addpd	%%xmm0,%%xmm6				\n\t			addpd	%%xmm2	,%%xmm5				\n\t"\
		"addpd	%%xmm1,%%xmm7				\n\t			addpd	%%xmm3	,%%xmm4				\n\t"\
		"movaps	%%xmm6,     (%%esi)			\n\t			movaps	%%xmm5	,0x060(%%esi)		\n\t"\
		"movaps	%%xmm7,0x010(%%esi)			\n\t			movaps	%%xmm4	,0x030(%%esi)		\n\t"\
		"\n\t"\
		"movl	%[__p2]		,%%edi			\n\t"\
		"shll	$3			,%%edi			\n\t"\
		"addl	%%edi		,%%eax			\n\t"\
		"addl	%%edi		,%%ebx			\n\t"\
		"addl	%%edi		,%%ecx			\n\t"\
		"addl	%%edi		,%%edx			\n\t"\
		"/* SSE2_RADIX4_DIF_4TWIDDLE_B(r9 ,c2) */\n\t"\
		"/* Do the p0,p8 combo: */			\n\t"\
		"movl	%[__c2]		,%%esi 			\n\t"\
		"movaps	    (%%eax)	,%%xmm0			\n\t			movaps	    (%%ecx)	,%%xmm4			\n\t"\
		"movaps	0x10(%%eax)	,%%xmm1			\n\t			movaps	0x10(%%ecx)	,%%xmm5			\n\t"\
		"movaps	    (%%esi)	,%%xmm6			\n\t"\
		"movaps	0x10(%%esi)	,%%xmm7			\n\t"\
		"movaps	%%xmm0		,%%xmm2			\n\t"\
		"movaps	%%xmm1		,%%xmm3			\n\t"\
		"mulpd   %%xmm6		,%%xmm0			\n\t"\
		"mulpd   %%xmm6		,%%xmm1			\n\t"\
		"mulpd   %%xmm7		,%%xmm2			\n\t"\
		"mulpd   %%xmm7		,%%xmm3			\n\t"\
		"													movaps	%%xmm4		,%%xmm6			\n\t"\
		"addpd   %%xmm2		,%%xmm1			\n\t			movaps	%%xmm5		,%%xmm7			\n\t"\
		"													mulpd   0x20(%%esi) ,%%xmm4			\n\t"\
		"subpd   %%xmm3		,%%xmm0			\n\t			mulpd   0x20(%%esi)	,%%xmm5			\n\t"\
		"													mulpd   0x30(%%esi)	,%%xmm6			\n\t"\
		"movaps	%%xmm0		,%%xmm2			\n\t			mulpd	0x30(%%esi)	,%%xmm7			\n\t"\
		"													addpd   %%xmm6	    ,%%xmm5			\n\t"\
		"movaps	%%xmm1		,%%xmm3			\n\t			subpd	%%xmm7		,%%xmm4			\n\t"\
		"addpd	%%xmm4		,%%xmm0			\n\t"\
		"addpd	%%xmm5		,%%xmm1			\n\t"\
		"subpd	%%xmm4		,%%xmm2			\n\t"\
		"subpd	%%xmm5		,%%xmm3			\n\t"\
		"/* Do the p4,12 combo: */			\n\t"\
		"addl	$0x60		,%%esi 			\n\t"\
		"movaps	    (%%edx)	,%%xmm4			\n\t"\
		"movaps	0x10(%%edx)	,%%xmm5			\n\t"\
		"movaps	    (%%edx)	,%%xmm6			\n\t"\
		"movaps	0x10(%%edx)	,%%xmm7			\n\t"\
		"mulpd	    (%%esi)	,%%xmm4			\n\t"\
		"mulpd	    (%%esi)	,%%xmm5			\n\t"\
		"mulpd	0x10(%%esi)	,%%xmm6			\n\t"\
		"mulpd	0x10(%%esi)	,%%xmm7			\n\t"\
		"movl	%[__r9]		,%%esi 			\n\t"\
		"addpd	%%xmm6		,%%xmm5			\n\t"\
		"subpd	%%xmm7		,%%xmm4			\n\t"\
		"movaps	%%xmm5		,0x010(%%esi)	\n\t"\
		"movaps	%%xmm4		,     (%%esi)	\n\t"\
		"movl	%[__c2]		,%%esi 			\n\t"\
		"addl	$0x40		,%%esi 			\n\t"\
		"movaps	    (%%ebx)	,%%xmm4			\n\t"\
		"movaps	0x10(%%ebx)	,%%xmm5			\n\t"\
		"movaps	    (%%ebx)	,%%xmm6			\n\t"\
		"movaps	0x10(%%ebx)	,%%xmm7			\n\t"\
		"mulpd	    (%%esi)	,%%xmm4			\n\t"\
		"mulpd	    (%%esi)	,%%xmm5			\n\t"\
		"mulpd	0x10(%%esi)	,%%xmm6			\n\t"\
		"mulpd	0x10(%%esi)	,%%xmm7			\n\t"\
		"movl	%[__r9]		,%%esi 			\n\t"\
		"addpd	%%xmm6		,%%xmm5			\n\t"\
		"subpd	%%xmm7		,%%xmm4			\n\t"\
		"movaps	%%xmm5		,%%xmm7			\n\t"\
		"movaps	%%xmm4		,%%xmm6			\n\t"\
		"subpd	     (%%esi),%%xmm4			\n\t"\
		"subpd	0x010(%%esi),%%xmm5			\n\t"\
		"addpd	     (%%esi),%%xmm6			\n\t"\
		"addpd	0x010(%%esi),%%xmm7			\n\t"\
		"/* Finish radix-4 butterfly and store results into temporary-array slots: */	\n\t"\
		"subpd	%%xmm6,%%xmm0				\n\t			subpd	%%xmm5	,%%xmm2				\n\t"\
		"subpd	%%xmm7,%%xmm1				\n\t			subpd	%%xmm4	,%%xmm3				\n\t"\
		"movaps	%%xmm0,0x040(%%esi)			\n\t			movaps	%%xmm2	,0x020(%%esi)		\n\t"\
		"movaps	%%xmm1,0x050(%%esi)			\n\t			movaps	%%xmm3	,0x070(%%esi)		\n\t"\
		"addpd	%%xmm6,%%xmm6				\n\t			addpd	%%xmm5	,%%xmm5				\n\t"\
		"addpd	%%xmm7,%%xmm7				\n\t			addpd	%%xmm4	,%%xmm4				\n\t"\
		"addpd	%%xmm0,%%xmm6				\n\t			addpd	%%xmm2	,%%xmm5				\n\t"\
		"addpd	%%xmm1,%%xmm7				\n\t			addpd	%%xmm3	,%%xmm4				\n\t"\
		"movaps	%%xmm6,     (%%esi)			\n\t			movaps	%%xmm5	,0x060(%%esi)		\n\t"\
		"movaps	%%xmm7,0x010(%%esi)			\n\t			movaps	%%xmm4	,0x030(%%esi)		\n\t"\
		"\n\t"\
		"subl	%%edi		,%%eax			\n\t"\
		"subl	%%edi		,%%ebx			\n\t"\
		"subl	%%edi		,%%ecx			\n\t"\
		"subl	%%edi		,%%edx			\n\t"\
		"movl	%[__p1]		,%%edi			\n\t"\
		"shll	$3			,%%edi			\n\t"\
		"addl	%%edi		,%%eax			\n\t"\
		"addl	%%edi		,%%ebx			\n\t"\
		"addl	%%edi		,%%ecx			\n\t"\
		"addl	%%edi		,%%edx			\n\t"\
		"/* SSE2_RADIX4_DIF_4TWIDDLE_B(r17,c1) */	\n\t"\
		"/* Do the p0,p8 combo: */			\n\t"\
		"movl	%[__c1]		,%%esi 			\n\t"\
		"movaps	    (%%eax)	,%%xmm0			\n\t			movaps	    (%%ecx)	,%%xmm4			\n\t"\
		"movaps	0x10(%%eax)	,%%xmm1			\n\t			movaps	0x10(%%ecx)	,%%xmm5			\n\t"\
		"movaps	    (%%esi)	,%%xmm6			\n\t"\
		"movaps	0x10(%%esi)	,%%xmm7			\n\t"\
		"movaps	%%xmm0		,%%xmm2			\n\t"\
		"movaps	%%xmm1		,%%xmm3			\n\t"\
		"mulpd   %%xmm6		,%%xmm0			\n\t"\
		"mulpd   %%xmm6		,%%xmm1			\n\t"\
		"mulpd   %%xmm7		,%%xmm2			\n\t"\
		"mulpd   %%xmm7		,%%xmm3			\n\t"\
		"													movaps	%%xmm4		,%%xmm6			\n\t"\
		"addpd   %%xmm2		,%%xmm1			\n\t			movaps	%%xmm5		,%%xmm7			\n\t"\
		"													mulpd   0x20(%%esi) ,%%xmm4			\n\t"\
		"subpd   %%xmm3		,%%xmm0			\n\t			mulpd   0x20(%%esi)	,%%xmm5			\n\t"\
		"													mulpd   0x30(%%esi)	,%%xmm6			\n\t"\
		"movaps	%%xmm0		,%%xmm2			\n\t			mulpd	0x30(%%esi)	,%%xmm7			\n\t"\
		"													addpd   %%xmm6	    ,%%xmm5			\n\t"\
		"movaps	%%xmm1		,%%xmm3			\n\t			subpd	%%xmm7		,%%xmm4			\n\t"\
		"addpd	%%xmm4		,%%xmm0			\n\t"\
		"addpd	%%xmm5		,%%xmm1			\n\t"\
		"subpd	%%xmm4		,%%xmm2			\n\t"\
		"subpd	%%xmm5		,%%xmm3			\n\t"\
		"/* Do the p4,12 combo: */			\n\t"\
		"addl	$0x60		,%%esi 			\n\t"\
		"movaps	    (%%edx)	,%%xmm4			\n\t"\
		"movaps	0x10(%%edx)	,%%xmm5			\n\t"\
		"movaps	    (%%edx)	,%%xmm6			\n\t"\
		"movaps	0x10(%%edx)	,%%xmm7			\n\t"\
		"mulpd	    (%%esi)	,%%xmm4			\n\t"\
		"mulpd	    (%%esi)	,%%xmm5			\n\t"\
		"mulpd	0x10(%%esi)	,%%xmm6			\n\t"\
		"mulpd	0x10(%%esi)	,%%xmm7			\n\t"\
		"movl	%[__r17]	,%%esi 			\n\t"\
		"addpd	%%xmm6		,%%xmm5			\n\t"\
		"subpd	%%xmm7		,%%xmm4			\n\t"\
		"movaps	%%xmm5		,0x010(%%esi)	\n\t"\
		"movaps	%%xmm4		,     (%%esi)	\n\t"\
		"movl	%[__c1]		,%%esi 			\n\t"\
		"addl	$0x40		,%%esi 			\n\t"\
		"movaps	    (%%ebx)	,%%xmm4			\n\t"\
		"movaps	0x10(%%ebx)	,%%xmm5			\n\t"\
		"movaps	    (%%ebx)	,%%xmm6			\n\t"\
		"movaps	0x10(%%ebx)	,%%xmm7			\n\t"\
		"mulpd	    (%%esi)	,%%xmm4			\n\t"\
		"mulpd	    (%%esi)	,%%xmm5			\n\t"\
		"mulpd	0x10(%%esi)	,%%xmm6			\n\t"\
		"mulpd	0x10(%%esi)	,%%xmm7			\n\t"\
		"movl	%[__r17]	,%%esi 			\n\t"\
		"addpd	%%xmm6		,%%xmm5			\n\t"\
		"subpd	%%xmm7		,%%xmm4			\n\t"\
		"movaps	%%xmm5		,%%xmm7			\n\t"\
		"movaps	%%xmm4		,%%xmm6			\n\t"\
		"subpd	     (%%esi),%%xmm4			\n\t"\
		"subpd	0x010(%%esi),%%xmm5			\n\t"\
		"addpd	     (%%esi),%%xmm6			\n\t"\
		"addpd	0x010(%%esi),%%xmm7			\n\t"\
		"/* Finish radix-4 butterfly and store results into temporary-array slots: */	\n\t"\
		"subpd	%%xmm6,%%xmm0				\n\t			subpd	%%xmm5	,%%xmm2				\n\t"\
		"subpd	%%xmm7,%%xmm1				\n\t			subpd	%%xmm4	,%%xmm3				\n\t"\
		"movaps	%%xmm0,0x040(%%esi)			\n\t			movaps	%%xmm2	,0x020(%%esi)		\n\t"\
		"movaps	%%xmm1,0x050(%%esi)			\n\t			movaps	%%xmm3	,0x070(%%esi)		\n\t"\
		"addpd	%%xmm6,%%xmm6				\n\t			addpd	%%xmm5	,%%xmm5				\n\t"\
		"addpd	%%xmm7,%%xmm7				\n\t			addpd	%%xmm4	,%%xmm4				\n\t"\
		"addpd	%%xmm0,%%xmm6				\n\t			addpd	%%xmm2	,%%xmm5				\n\t"\
		"addpd	%%xmm1,%%xmm7				\n\t			addpd	%%xmm3	,%%xmm4				\n\t"\
		"movaps	%%xmm6,     (%%esi)			\n\t			movaps	%%xmm5	,0x060(%%esi)		\n\t"\
		"movaps	%%xmm7,0x010(%%esi)			\n\t			movaps	%%xmm4	,0x030(%%esi)		\n\t"\
		"\n\t"\
		"subl	%%edi		,%%eax			\n\t"\
		"subl	%%edi		,%%ebx			\n\t"\
		"subl	%%edi		,%%ecx			\n\t"\
		"subl	%%edi		,%%edx			\n\t"\
		"movl	%[__p3]		,%%edi			\n\t"\
		"shll	$3			,%%edi			\n\t"\
		"addl	%%edi		,%%eax			\n\t"\
		"addl	%%edi		,%%ebx			\n\t"\
		"addl	%%edi		,%%ecx			\n\t"\
		"addl	%%edi		,%%edx			\n\t"\
		"/* SSE2_RADIX4_DIF_4TWIDDLE_B(r25,c3) */\n\t"\
		"/* Do the p0,p8 combo: */			\n\t"\
		"movl	%[__c3]		,%%esi 			\n\t"\
		"movaps	    (%%eax)	,%%xmm0			\n\t			movaps	    (%%ecx)	,%%xmm4			\n\t"\
		"movaps	0x10(%%eax)	,%%xmm1			\n\t			movaps	0x10(%%ecx)	,%%xmm5			\n\t"\
		"movaps	    (%%esi)	,%%xmm6			\n\t"\
		"movaps	0x10(%%esi)	,%%xmm7			\n\t"\
		"movaps	%%xmm0		,%%xmm2			\n\t"\
		"movaps	%%xmm1		,%%xmm3			\n\t"\
		"mulpd   %%xmm6		,%%xmm0			\n\t"\
		"mulpd   %%xmm6		,%%xmm1			\n\t"\
		"mulpd   %%xmm7		,%%xmm2			\n\t"\
		"mulpd   %%xmm7		,%%xmm3			\n\t"\
		"													movaps	%%xmm4		,%%xmm6			\n\t"\
		"addpd   %%xmm2		,%%xmm1			\n\t			movaps	%%xmm5		,%%xmm7			\n\t"\
		"													mulpd   0x20(%%esi) ,%%xmm4			\n\t"\
		"subpd   %%xmm3		,%%xmm0			\n\t			mulpd   0x20(%%esi)	,%%xmm5			\n\t"\
		"													mulpd   0x30(%%esi)	,%%xmm6			\n\t"\
		"movaps	%%xmm0		,%%xmm2			\n\t			mulpd	0x30(%%esi)	,%%xmm7			\n\t"\
		"													addpd   %%xmm6	    ,%%xmm5			\n\t"\
		"movaps	%%xmm1		,%%xmm3			\n\t			subpd	%%xmm7		,%%xmm4			\n\t"\
		"addpd	%%xmm4		,%%xmm0			\n\t"\
		"addpd	%%xmm5		,%%xmm1			\n\t"\
		"subpd	%%xmm4		,%%xmm2			\n\t"\
		"subpd	%%xmm5		,%%xmm3			\n\t"\
		"/* Do the p4,12 combo: */			\n\t"\
		"addl	$0x60		,%%esi 			\n\t"\
		"movaps	    (%%edx)	,%%xmm4			\n\t"\
		"movaps	0x10(%%edx)	,%%xmm5			\n\t"\
		"movaps	    (%%edx)	,%%xmm6			\n\t"\
		"movaps	0x10(%%edx)	,%%xmm7			\n\t"\
		"mulpd	    (%%esi)	,%%xmm4			\n\t"\
		"mulpd	    (%%esi)	,%%xmm5			\n\t"\
		"mulpd	0x10(%%esi)	,%%xmm6			\n\t"\
		"mulpd	0x10(%%esi)	,%%xmm7			\n\t"\
		"movl	%[__r25]	,%%esi 			\n\t"\
		"addpd	%%xmm6		,%%xmm5			\n\t"\
		"subpd	%%xmm7		,%%xmm4			\n\t"\
		"movaps	%%xmm5		,0x010(%%esi)	\n\t"\
		"movaps	%%xmm4		,     (%%esi)	\n\t"\
		"movl	%[__c3]		,%%esi 			\n\t"\
		"addl	$0x40		,%%esi 			\n\t"\
		"movaps	    (%%ebx)	,%%xmm4			\n\t"\
		"movaps	0x10(%%ebx)	,%%xmm5			\n\t"\
		"movaps	    (%%ebx)	,%%xmm6			\n\t"\
		"movaps	0x10(%%ebx)	,%%xmm7			\n\t"\
		"mulpd	    (%%esi)	,%%xmm4			\n\t"\
		"mulpd	    (%%esi)	,%%xmm5			\n\t"\
		"mulpd	0x10(%%esi)	,%%xmm6			\n\t"\
		"mulpd	0x10(%%esi)	,%%xmm7			\n\t"\
		"movl	%[__r25]	,%%esi 			\n\t"\
		"addpd	%%xmm6		,%%xmm5			\n\t"\
		"subpd	%%xmm7		,%%xmm4			\n\t"\
		"movaps	%%xmm5		,%%xmm7			\n\t"\
		"movaps	%%xmm4		,%%xmm6			\n\t"\
		"subpd	     (%%esi),%%xmm4			\n\t"\
		"subpd	0x010(%%esi),%%xmm5			\n\t"\
		"addpd	     (%%esi),%%xmm6			\n\t"\
		"addpd	0x010(%%esi),%%xmm7			\n\t"\
		"/* Finish radix-4 butterfly and store results into temporary-array slots: */			\n\t"\
		"subpd	%%xmm6,%%xmm0				\n\t			subpd	%%xmm5	,%%xmm2				\n\t"\
		"subpd	%%xmm7,%%xmm1				\n\t			subpd	%%xmm4	,%%xmm3				\n\t"\
		"movaps	%%xmm0,0x040(%%esi)			\n\t			movaps	%%xmm2	,0x020(%%esi)		\n\t"\
		"movaps	%%xmm1,0x050(%%esi)			\n\t			movaps	%%xmm3	,0x070(%%esi)		\n\t"\
		"addpd	%%xmm6,%%xmm6				\n\t			addpd	%%xmm5	,%%xmm5				\n\t"\
		"addpd	%%xmm7,%%xmm7				\n\t			addpd	%%xmm4	,%%xmm4				\n\t"\
		"addpd	%%xmm0,%%xmm6				\n\t			addpd	%%xmm2	,%%xmm5				\n\t"\
		"addpd	%%xmm1,%%xmm7				\n\t			addpd	%%xmm3	,%%xmm4				\n\t"\
		"movaps	%%xmm6,     (%%esi)			\n\t			movaps	%%xmm5	,0x060(%%esi)		\n\t"\
		"movaps	%%xmm7,0x010(%%esi)			\n\t			movaps	%%xmm4	,0x030(%%esi)		\n\t"\
		"/*************************************************************************************/\n\t"\
		"/*  And now do four more radix-4 transforms, including the internal twiddle factors: */\n\t"\
		"/*************************************************************************************/\n\t"\
		"/* Block 1: */\n\t"\
		"movl	%[__add0]	,%%eax			\n\t"\
		"movl	%[__p1]		,%%ebx			\n\t"\
		"movl	%[__p2]		,%%ecx			\n\t"\
		"movl	%[__p3]		,%%edx			\n\t"\
		"shll	$3			,%%ebx			\n\t"\
		"shll	$3			,%%ecx			\n\t"\
		"shll	$3			,%%edx			\n\t"\
		"addl	%%eax		,%%ebx			\n\t"\
		"addl	%%eax		,%%ecx			\n\t"\
		"addl	%%eax		,%%edx			\n\t"\
		"\n\t"\
		"movl	%[__r1]		,%%edi			\n\t"\
		"movaps	     (%%edi),%%xmm0			\n\t		movaps	0x100(%%edi),%%xmm4			\n\t"\
		"movaps	0x010(%%edi),%%xmm1			\n\t		movaps	0x110(%%edi),%%xmm5			\n\t"\
		"movaps	0x080(%%edi),%%xmm2			\n\t		movaps	0x180(%%edi),%%xmm6			\n\t"\
		"movaps	0x090(%%edi),%%xmm3			\n\t		movaps	0x190(%%edi),%%xmm7			\n\t"\
		"\n\t"\
		"subpd	%%xmm2		,%%xmm0			\n\t		subpd	%%xmm6		,%%xmm4			\n\t"\
		"subpd	%%xmm3		,%%xmm1			\n\t		subpd	%%xmm7		,%%xmm5			\n\t"\
		"addpd	%%xmm2		,%%xmm2			\n\t		addpd	%%xmm6		,%%xmm6			\n\t"\
		"addpd	%%xmm3		,%%xmm3			\n\t		addpd	%%xmm7		,%%xmm7			\n\t"\
		"addpd	%%xmm0		,%%xmm2			\n\t		addpd	%%xmm4		,%%xmm6			\n\t"\
		"addpd	%%xmm1		,%%xmm3			\n\t		addpd	%%xmm5		,%%xmm7			\n\t"\
		"\n\t"\
		"subpd	%%xmm6		,%%xmm2			\n\t"\
		"subpd	%%xmm7		,%%xmm3			\n\t"\
		"addpd	%%xmm6		,%%xmm6			\n\t"\
		"addpd	%%xmm7		,%%xmm7			\n\t"\
		"movaps	%%xmm2		,    (%%ebx)	\n\t"\
		"movaps	%%xmm3		,0x10(%%ebx)	\n\t"\
		"addpd	%%xmm2		,	%%xmm6		\n\t"\
		"addpd	%%xmm3		,	%%xmm7		\n\t"\
		"movaps	%%xmm6		,    (%%eax)	\n\t"\
		"movaps	%%xmm7		,0x10(%%eax)	\n\t"\
		"\n\t"\
		"subpd	%%xmm5		,	%%xmm0		\n\t"\
		"subpd	%%xmm4		,	%%xmm1		\n\t"\
		"addpd	%%xmm5		,	%%xmm5		\n\t"\
		"addpd	%%xmm4		,	%%xmm4		\n\t"\
		"movaps	%%xmm0		,    (%%ecx)	\n\t"\
		"movaps	%%xmm1		,0x10(%%edx)	\n\t"\
		"addpd	%%xmm0		,	%%xmm5		\n\t"\
		"addpd	%%xmm1		,	%%xmm4		\n\t"\
		"movaps	%%xmm5		,    (%%edx)	\n\t"\
		"movaps	%%xmm4		,0x10(%%ecx)	\n\t"\
		"\n\t"\
		"/* Block 3: */\n\t"\
		"movl	%[__p4]		,%%edi			\n\t"\
		"shll	$3			,%%edi			\n\t"\
		"addl	%%edi		,%%eax			\n\t"\
		"addl	%%edi		,%%ebx			\n\t"\
		"addl	%%edi		,%%ecx			\n\t"\
		"addl	%%edi		,%%edx			\n\t"\
		"\n\t"\
		"movl	%[__r5]		,%%edi			\n\t"\
		"movl	%[__isrt2]	,%%esi			\n\t"\
		"movaps	(%%esi)		,%%xmm3			/* isrt2 */	\n\t"\
		"									\n\t		movaps	0x100(%%edi),%%xmm4			\n\t"\
		"									\n\t		movaps	0x110(%%edi),%%xmm5			\n\t"\
		"									\n\t		movaps	0x180(%%edi),%%xmm6			\n\t"\
		"									\n\t		movaps	0x190(%%edi),%%xmm7			\n\t"\
		"									\n\t		mulpd	%%xmm3		,%%xmm4			\n\t"\
		"movaps	     (%%edi),%%xmm0			\n\t		mulpd	%%xmm3		,%%xmm5			\n\t"\
		"movaps	0x010(%%edi),%%xmm1			\n\t		mulpd	%%xmm3		,%%xmm6			\n\t"\
		"movaps	0x080(%%edi),%%xmm2			\n\t		mulpd	%%xmm3		,%%xmm7			\n\t"\
		"movaps	0x090(%%edi),%%xmm3			\n\t"\
		"\n\t"\
		"subpd	%%xmm3		,%%xmm0			\n\t		subpd	%%xmm5		,%%xmm4			\n\t"\
		"subpd	%%xmm2		,%%xmm1			\n\t		subpd	%%xmm6		,%%xmm7			\n\t"\
		"addpd	%%xmm3		,%%xmm3			\n\t		addpd	%%xmm5		,%%xmm5			\n\t"\
		"addpd	%%xmm2		,%%xmm2			\n\t		addpd	%%xmm6		,%%xmm6			\n\t"\
		"addpd	%%xmm0		,%%xmm3			\n\t		addpd	%%xmm4		,%%xmm5			\n\t"\
		"addpd	%%xmm1		,%%xmm2			\n\t		addpd	%%xmm7		,%%xmm6			\n\t"\
		"\n\t"\
		"subpd	%%xmm6		,%%xmm4			\n\t"\
		"subpd	%%xmm7		,%%xmm5			\n\t"\
		"addpd	%%xmm6		,%%xmm6			\n\t"\
		"addpd	%%xmm7		,%%xmm7			\n\t"\
		"addpd	%%xmm4		,%%xmm6			\n\t"\
		"addpd	%%xmm5		,%%xmm7			\n\t"\
		"\n\t"\
		"subpd	%%xmm4		,%%xmm0			\n\t		subpd	%%xmm7		,%%xmm3			\n\t"\
		"subpd	%%xmm5		,%%xmm2			\n\t		subpd	%%xmm6		,%%xmm1			\n\t"\
		"addpd	%%xmm4		,%%xmm4			\n\t		addpd	%%xmm7		,%%xmm7			\n\t"\
		"addpd	%%xmm5		,%%xmm5			\n\t		addpd	%%xmm6		,%%xmm6			\n\t"\
		"movaps	%%xmm0		,    (%%ebx)	\n\t		movaps	%%xmm3		,    (%%ecx)	\n\t"\
		"movaps	%%xmm2		,0x10(%%ebx)	\n\t		movaps	%%xmm1		,0x10(%%edx)	\n\t"\
		"addpd	%%xmm0		,	%%xmm4		\n\t		addpd	%%xmm3		,	%%xmm7		\n\t"\
		"addpd	%%xmm2		,	%%xmm5		\n\t		addpd	%%xmm1		,	%%xmm6		\n\t"\
		"movaps	%%xmm4		,    (%%eax)	\n\t		movaps	%%xmm7		,    (%%edx)	\n\t"\
		"movaps	%%xmm5		,0x10(%%eax)	\n\t		movaps	%%xmm6		,0x10(%%ecx)	\n\t"\
		"\n\t"\
		"/* Block 2: */\n\t"\
		"movl	%[__p4]		,%%edi			\n\t"\
		"shll	$3			,%%edi			\n\t"\
		"addl	%%edi		,%%eax			\n\t"\
		"addl	%%edi		,%%ebx			\n\t"\
		"addl	%%edi		,%%ecx			\n\t"\
		"addl	%%edi		,%%edx			\n\t"\
		"\n\t"\
		"movl	%[__r3]		,%%edi			\n\t"\
		"movl	%[__cc0]	,%%esi			\n\t"\
		"movaps	0x100(%%edi),%%xmm4			\n\t"\
		"movaps	0x110(%%edi),%%xmm5			\n\t"\
		"movaps	     (%%esi),%%xmm3			\n\t"\
		"movaps	0x010(%%esi),%%xmm2			\n\t"\
		"movaps	%%xmm4		,%%xmm6			\n\t"\
		"movaps	%%xmm5		,%%xmm7			\n\t"\
		"\n\t"\
		"mulpd	%%xmm3		,%%xmm4			\n\t"\
		"mulpd	%%xmm3		,%%xmm5			\n\t"\
		"mulpd	%%xmm2		,%%xmm6			\n\t		movaps	0x180(%%edi),%%xmm0			\n\t"\
		"mulpd	%%xmm2		,%%xmm7			\n\t		movaps	0x190(%%edi),%%xmm1			\n\t"\
		"addpd	%%xmm6		,%%xmm5			\n\t		movaps	%%xmm0		,%%xmm6			\n\t"\
		"subpd	%%xmm7		,%%xmm4			\n\t		movaps	%%xmm1		,%%xmm7			\n\t"\
		"\n\t"\
		"									\n\t		mulpd	%%xmm2		,%%xmm6			\n\t"\
		"									\n\t		mulpd	%%xmm2		,%%xmm7			\n\t"\
		"									\n\t		mulpd	%%xmm3		,%%xmm0			\n\t"\
		"									\n\t		mulpd	%%xmm3		,%%xmm1			\n\t"\
		"									\n\t		addpd	%%xmm0		,%%xmm7			\n\t"\
		"									\n\t		subpd	%%xmm1		,%%xmm6			\n\t"\
		"\n\t"\
		"movaps	%%xmm4		,%%xmm2			\n\t"\
		"movaps	%%xmm5		,%%xmm3			\n\t"\
		"subpd	%%xmm6		,%%xmm4			\n\t"\
		"subpd	%%xmm7		,%%xmm5			\n\t"\
		"addpd	%%xmm2		,%%xmm6			\n\t"\
		"addpd	%%xmm3		,%%xmm7			\n\t"\
		"\n\t"\
		"movl	%[__isrt2]	,%%esi			\n\t"\
		"movaps	0x080(%%edi),%%xmm2			\n\t"\
		"movaps	0x090(%%edi),%%xmm3			\n\t"\
		"movaps	(%%esi)		,%%xmm1			/* isrt2 */	\n\t"\
		"movaps	%%xmm2		,%%xmm0			\n\t"\
		"subpd	%%xmm3		,%%xmm2			\n\t"\
		"addpd	%%xmm0		,%%xmm3			\n\t"\
		"mulpd	%%xmm1		,%%xmm2			\n\t"\
		"mulpd	%%xmm1		,%%xmm3			\n\t"\
		"\n\t"\
		"movaps	     (%%edi),%%xmm0			\n\t"\
		"movaps	0x010(%%edi),%%xmm1			\n\t"\
		"\n\t"\
		"subpd	%%xmm2		,%%xmm0			\n\t"\
		"subpd	%%xmm3		,%%xmm1			\n\t"\
		"addpd	%%xmm2		,%%xmm2			\n\t"\
		"addpd	%%xmm3		,%%xmm3			\n\t"\
		"addpd	%%xmm0		,%%xmm2			\n\t"\
		"addpd	%%xmm1		,%%xmm3			\n\t"\
		"\n\t"\
		"subpd	%%xmm6		,%%xmm2			\n\t		subpd	%%xmm5		,%%xmm0			\n\t"\
		"subpd	%%xmm7		,%%xmm3			\n\t		subpd	%%xmm4		,%%xmm1			\n\t"\
		"addpd	%%xmm6		,%%xmm6			\n\t		addpd	%%xmm5		,%%xmm5			\n\t"\
		"addpd	%%xmm7		,%%xmm7			\n\t		addpd	%%xmm4		,%%xmm4			\n\t"\
		"movaps	%%xmm2		,    (%%ebx)	\n\t		movaps	%%xmm0		,    (%%ecx)	\n\t"\
		"movaps	%%xmm3		,0x10(%%ebx)	\n\t		movaps	%%xmm1		,0x10(%%edx)	\n\t"\
		"addpd	%%xmm2		,	%%xmm6		\n\t		addpd	%%xmm0		,	%%xmm5		\n\t"\
		"addpd	%%xmm3		,	%%xmm7		\n\t		addpd	%%xmm1		,	%%xmm4		\n\t"\
		"movaps	%%xmm6		,    (%%eax)	\n\t		movaps	%%xmm5		,    (%%edx)	\n\t"\
		"movaps	%%xmm7		,0x10(%%eax)	\n\t		movaps	%%xmm4		,0x10(%%ecx)	\n\t"\
		"/* Block 4: */\n\t"\
		"movl	%[__p4]		,%%edi			\n\t"\
		"shll	$3			,%%edi			\n\t"\
		"addl	%%edi		,%%eax			\n\t"\
		"addl	%%edi		,%%ebx			\n\t"\
		"addl	%%edi		,%%ecx			\n\t"\
		"addl	%%edi		,%%edx			\n\t"\
		"\n\t"\
		"movl	%[__r7]		,%%edi			\n\t"\
		"movl	%[__cc0]	,%%esi			\n\t"\
		"movaps	0x100(%%edi),%%xmm4			\n\t"\
		"movaps	0x110(%%edi),%%xmm5			\n\t"\
		"movaps	     (%%esi),%%xmm2			\n\t"\
		"movaps	0x010(%%esi),%%xmm3			\n\t"\
		"movaps	%%xmm4		,%%xmm6			\n\t"\
		"movaps	%%xmm5		,%%xmm7			\n\t"\
		"\n\t"\
		"mulpd	%%xmm3		,%%xmm4			\n\t"\
		"mulpd	%%xmm3		,%%xmm5			\n\t"\
		"mulpd	%%xmm2		,%%xmm6			\n\t		movaps	0x180(%%edi),%%xmm0			\n\t"\
		"mulpd	%%xmm2		,%%xmm7			\n\t		movaps	0x190(%%edi),%%xmm1			\n\t"\
		"addpd	%%xmm6		,%%xmm5			\n\t		movaps	%%xmm0		,%%xmm6			\n\t"\
		"subpd	%%xmm7		,%%xmm4			\n\t		movaps	%%xmm1		,%%xmm7			\n\t"\
		"\n\t"\
		"									\n\t		mulpd	%%xmm2		,%%xmm6			\n\t"\
		"									\n\t		mulpd	%%xmm2		,%%xmm7			\n\t"\
		"									\n\t		mulpd	%%xmm3		,%%xmm0			\n\t"\
		"									\n\t		mulpd	%%xmm3		,%%xmm1			\n\t"\
		"									\n\t		addpd	%%xmm0		,%%xmm7			\n\t"\
		"									\n\t		subpd	%%xmm1		,%%xmm6			\n\t"\
		"\n\t"\
		"movaps	%%xmm4		,%%xmm2			\n\t"\
		"movaps	%%xmm5		,%%xmm3			\n\t"\
		"subpd	%%xmm6		,%%xmm4			\n\t"\
		"subpd	%%xmm7		,%%xmm5			\n\t"\
		"addpd	%%xmm2		,%%xmm6			\n\t"\
		"addpd	%%xmm3		,%%xmm7			\n\t"\
		"\n\t"\
		"movl	%[__isrt2]	,%%esi			\n\t"\
		"movaps	0x080(%%edi),%%xmm2			\n\t"\
		"movaps	0x090(%%edi),%%xmm3			\n\t"\
		"movaps	(%%esi)		,%%xmm1			/* isrt2 */	\n\t"\
		"movaps	%%xmm2		,%%xmm0			\n\t"\
		"addpd	%%xmm3		,%%xmm2			\n\t"\
		"subpd	%%xmm0		,%%xmm3			\n\t"\
		"mulpd	%%xmm1		,%%xmm2			\n\t"\
		"mulpd	%%xmm1		,%%xmm3			\n\t"\
		"\n\t"\
		"movaps	     (%%edi),%%xmm0			\n\t"\
		"movaps	0x010(%%edi),%%xmm1			\n\t"\
		"\n\t"\
		"subpd	%%xmm2		,%%xmm0			\n\t"\
		"subpd	%%xmm3		,%%xmm1			\n\t"\
		"addpd	%%xmm2		,%%xmm2			\n\t"\
		"addpd	%%xmm3		,%%xmm3			\n\t"\
		"addpd	%%xmm0		,%%xmm2			\n\t"\
		"addpd	%%xmm1		,%%xmm3			\n\t"\
		"\n\t"\
		"subpd	%%xmm4		,%%xmm0			\n\t		subpd	%%xmm7		,%%xmm2			\n\t"\
		"subpd	%%xmm5		,%%xmm1			\n\t		subpd	%%xmm6		,%%xmm3			\n\t"\
		"addpd	%%xmm4		,%%xmm4			\n\t		addpd	%%xmm7		,%%xmm7			\n\t"\
		"addpd	%%xmm5		,%%xmm5			\n\t		addpd	%%xmm6		,%%xmm6			\n\t"\
		"movaps	%%xmm0		,    (%%ebx)	\n\t		movaps	%%xmm2		,    (%%ecx)	\n\t"\
		"movaps	%%xmm1		,0x10(%%ebx)	\n\t		movaps	%%xmm3		,0x10(%%edx)	\n\t"\
		"addpd	%%xmm0		,	%%xmm4		\n\t		addpd	%%xmm2		,%%xmm7			\n\t"\
		"addpd	%%xmm1		,	%%xmm5		\n\t		addpd	%%xmm3		,%%xmm6			\n\t"\
		"movaps	%%xmm4		,    (%%eax)	\n\t		movaps	%%xmm7		,    (%%edx)	\n\t"\
		"movaps	%%xmm5		,0x10(%%eax)	\n\t		movaps	%%xmm6		,0x10(%%ecx)	\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__p1] "m" (Xp1)\
		 ,[__p2] "m" (Xp2)\
		 ,[__p3] "m" (Xp3)\
		 ,[__p4] "m" (Xp4)\
		 ,[__p8] "m" (Xp8)\
		 ,[__p12] "m" (Xp12)\
		 ,[__r1] "m" (Xr1)\
		 ,[__r3] "m" (Xr3)\
		 ,[__r5] "m" (Xr5)\
		 ,[__r7] "m" (Xr7)\
		 ,[__r9] "m" (Xr9)\
		 ,[__r17] "m" (Xr17)\
		 ,[__r25] "m" (Xr25)\
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__cc0] "m" (Xcc0)\
		 ,[__c0] "m" (Xc0)\
		 ,[__c1] "m" (Xc1)\
		 ,[__c2] "m" (Xc2)\
		 ,[__c3] "m" (Xc3)\
		: "eax","ebx","ecx","edx","edi","esi"		/* Clobbered registers */\
	);\
	}

#endif	/* radix16_dif_dit_pass_gcc_h_included */

