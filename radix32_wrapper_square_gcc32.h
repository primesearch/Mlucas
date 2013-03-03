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
#ifndef radix32_wrapper_square_gcc_h_included
#define radix32_wrapper_square_gcc_h_included

	#define SSE2_RADIX32_WRAPPER_DIF(Xadd0,Xadd1,Xr00,Xr10,Xr20,Xr30,Xisrt2,Xcc0,Xc00,Xc01,Xc02,Xc03,Xc05,Xc07)\
	{\
	__asm__ volatile (\
		"/************************************************************************/\n\t"\
		"/* Forward DIF radix-32 pass on the interleaved block1 and block2 data: */\n\t"\
		"/************************************************************************/\n\t"\
		"/*...Block 0:	*/\n\t"\
		"movl		%[__add0]		,%%eax\n\t"\
		"movl		%[__add1]		,%%ebx\n\t"\
		"/*****	SSE2_RADIX4_DIF_4WRAPPER(c00,c08,c10,c18,r00) *****/\n\t"\
		"movl		%[__r00]			,%%ecx\n\t"\
		"movl		%[__c00]			,%%edx\n\t"\
		"movaps		      (%%eax)	,%%xmm6\n\t"\
		"movaps		      (%%eax)	,%%xmm0\n\t"\
		"unpckhpd	      (%%ebx)	,%%xmm6\n\t"\
		"unpcklpd	      (%%ebx)	,%%xmm0\n\t"\
		"movaps		%%xmm6		, 0x200(%%ecx)\n\t"\
		"movaps		 0x010(%%eax)	,%%xmm7\n\t"\
		"movaps		 0x010(%%eax)	,%%xmm1\n\t"\
		"unpckhpd	 0x010(%%ebx)	,%%xmm7\n\t"\
		"unpcklpd	 0x010(%%ebx)	,%%xmm1\n\t"\
		"movaps		%%xmm7		, 0x210(%%ecx)\n\t"\
		"movaps		%%xmm0		,%%xmm2\n\t"\
		"movaps		%%xmm1		,%%xmm3\n\t"\
		"mulpd		      (%%edx)	,%%xmm0\n\t"\
		"mulpd		      (%%edx)	,%%xmm1\n\t"\
		"mulpd		 0x010(%%edx)	,%%xmm2\n\t"\
		"mulpd		 0x010(%%edx)	,%%xmm3\n\t"\
		"addpd		%%xmm2		,%%xmm1\n\t"\
		"subpd		%%xmm3		,%%xmm0\n\t"\
		"movaps		%%xmm0		,%%xmm2\n\t"\
		"movaps		%%xmm1		,%%xmm3\n\t"\
		"addl		$0x020		,%%edx\n\t"\
		"movaps		 0x100(%%eax)	,%%xmm6\n\t"\
		"movaps		 0x100(%%eax)	,%%xmm4\n\t"\
		"unpckhpd	 0x100(%%ebx)	,%%xmm6\n\t"\
		"unpcklpd	 0x100(%%ebx)	,%%xmm4\n\t"\
		"movaps		%%xmm6		, 0x220(%%ecx)\n\t"\
		"movaps		 0x110(%%eax)	,%%xmm7\n\t"\
		"movaps		 0x110(%%eax)	,%%xmm5\n\t"\
		"unpckhpd	 0x110(%%ebx)	,%%xmm7\n\t"\
		"unpcklpd	 0x110(%%ebx)	,%%xmm5\n\t"\
		"movaps		%%xmm7		, 0x230(%%ecx)\n\t"\
		"movaps		%%xmm4		,%%xmm6\n\t"\
		"movaps		%%xmm5		,%%xmm7\n\t"\
		"mulpd		      (%%edx)	,%%xmm4\n\t"\
		"mulpd		      (%%edx)	,%%xmm5\n\t"\
		"mulpd		 0x010(%%edx)	,%%xmm6\n\t"\
		"mulpd		 0x010(%%edx)	,%%xmm7\n\t"\
		"addpd		%%xmm6		,%%xmm5\n\t"\
		"subpd		%%xmm7		,%%xmm4\n\t"\
		"addpd		%%xmm4		,%%xmm0\n\t"\
		"addpd		%%xmm5		,%%xmm1\n\t"\
		"subpd		%%xmm4		,%%xmm2\n\t"\
		"subpd		%%xmm5		,%%xmm3\n\t"\
		"addl		$0x040		,%%edx\n\t"\
		"movaps		 0x180(%%eax)	,%%xmm6\n\t"\
		"movaps		 0x180(%%eax)	,%%xmm4\n\t"\
		"unpckhpd	 0x180(%%ebx)	,%%xmm6\n\t"\
		"unpcklpd	 0x180(%%ebx)	,%%xmm4\n\t"\
		"movaps		%%xmm6		, 0x260(%%ecx)\n\t"\
		"movaps		 0x190(%%eax)	,%%xmm7\n\t"\
		"movaps		 0x190(%%eax)	,%%xmm5\n\t"\
		"unpckhpd	 0x190(%%ebx)	,%%xmm7\n\t"\
		"unpcklpd	 0x190(%%ebx)	,%%xmm5\n\t"\
		"movaps		%%xmm7		, 0x270(%%ecx)\n\t"\
		"movaps		%%xmm4		,%%xmm6\n\t"\
		"movaps		%%xmm5		,%%xmm7\n\t"\
		"mulpd		      (%%edx)	,%%xmm4\n\t"\
		"mulpd		      (%%edx)	,%%xmm5\n\t"\
		"mulpd		 0x010(%%edx)	,%%xmm6\n\t"\
		"mulpd		 0x010(%%edx)	,%%xmm7\n\t"\
		"addpd		%%xmm6		,%%xmm5\n\t"\
		"subpd		%%xmm7		,%%xmm4\n\t"\
		"movaps		%%xmm5		, 0x010(%%ecx)\n\t"\
		"movaps		%%xmm4		,      (%%ecx)\n\t"\
		"subl		$0x020		,%%edx\n\t"\
		"movaps		 0x080(%%eax)	,%%xmm6\n\t"\
		"movaps		 0x080(%%eax)	,%%xmm4\n\t"\
		"unpckhpd	 0x080(%%ebx)	,%%xmm6\n\t"\
		"unpcklpd	 0x080(%%ebx)	,%%xmm4\n\t"\
		"movaps		%%xmm6		, 0x240(%%ecx)\n\t"\
		"movaps		 0x090(%%eax)	,%%xmm7\n\t"\
		"movaps		 0x090(%%eax)	,%%xmm5\n\t"\
		"unpckhpd	 0x090(%%ebx)	,%%xmm7\n\t"\
		"unpcklpd	 0x090(%%ebx)	,%%xmm5\n\t"\
		"movaps		%%xmm7		, 0x250(%%ecx)\n\t"\
		"movaps		%%xmm4		,%%xmm6\n\t"\
		"movaps		%%xmm5		,%%xmm7\n\t"\
		"mulpd		      (%%edx)	,%%xmm4\n\t"\
		"mulpd		      (%%edx)	,%%xmm5\n\t"\
		"mulpd		 0x010(%%edx)	,%%xmm6\n\t"\
		"mulpd		 0x010(%%edx)	,%%xmm7\n\t"\
		"addpd		%%xmm6		,%%xmm5\n\t"\
		"subpd		%%xmm7		,%%xmm4\n\t"\
		"movaps		%%xmm5		,%%xmm7\n\t"\
		"movaps		%%xmm4		,%%xmm6\n\t"\
		"subpd		      (%%ecx)	,%%xmm4\n\t"\
		"subpd		 0x010(%%ecx)	,%%xmm5\n\t"\
		"addpd		      (%%ecx)	,%%xmm6\n\t"\
		"addpd		 0x010(%%ecx)	,%%xmm7\n\t"\
		"subpd		%%xmm6		,%%xmm0\n\t"\
		"subpd		%%xmm7		,%%xmm1\n\t"\
		"movaps		%%xmm0		, 0x040(%%ecx)\n\t"\
		"movaps		%%xmm1		, 0x050(%%ecx)\n\t"\
		"addpd		%%xmm6		,%%xmm6\n\t"\
		"addpd		%%xmm7		,%%xmm7\n\t"\
		"addpd		%%xmm0		,%%xmm6\n\t"\
		"addpd		%%xmm1		,%%xmm7\n\t"\
		"movaps		%%xmm6		,      (%%ecx)\n\t"\
		"movaps		%%xmm7		, 0x010(%%ecx)\n\t"\
		"subpd		%%xmm5		,%%xmm2\n\t"\
		"subpd		%%xmm4		,%%xmm3\n\t"\
		"movaps		%%xmm2		, 0x020(%%ecx)\n\t"\
		"movaps		%%xmm3		, 0x070(%%ecx)\n\t"\
		"addpd		%%xmm5		,%%xmm5\n\t"\
		"addpd		%%xmm4		,%%xmm4\n\t"\
		"addpd		%%xmm2		,%%xmm5\n\t"\
		"addpd		%%xmm3		,%%xmm4\n\t"\
		"movaps		%%xmm5		, 0x060(%%ecx)\n\t"\
		"movaps		%%xmm4		, 0x030(%%ecx)\n\t"\
		"/*****	SSE2_RADIX4_DIF_4WRAPPER_2NDOFTWO(c04,c0C,c14,c1C,r08)	*****/\n\t"\
		"addl		$0x080		,%%ecx\n\t"\
		"addl		$0x040		,%%edx\n\t"\
		"movaps		 0x40 (%%eax)	,%%xmm6\n\t"\
		"movaps		 0x40 (%%eax)	,%%xmm0\n\t"\
		"unpckhpd	 0x40 (%%ebx)	,%%xmm6\n\t"\
		"unpcklpd	 0x40 (%%ebx)	,%%xmm0\n\t"\
		"movaps		%%xmm6		, 0x200(%%ecx)\n\t"\
		"movaps		 0x50 (%%eax)	,%%xmm7\n\t"\
		"movaps		 0x50 (%%eax)	,%%xmm1\n\t"\
		"unpckhpd	 0x50 (%%ebx)	,%%xmm7\n\t"\
		"unpcklpd	 0x50 (%%ebx)	,%%xmm1\n\t"\
		"movaps		%%xmm7		, 0x210(%%ecx)\n\t"\
		"movaps		%%xmm0		,%%xmm2\n\t"\
		"movaps		%%xmm1		,%%xmm3\n\t"\
		"mulpd		      (%%edx)	,%%xmm0\n\t"\
		"mulpd		      (%%edx)	,%%xmm1\n\t"\
		"mulpd		 0x010(%%edx)	,%%xmm2\n\t"\
		"mulpd		 0x010(%%edx)	,%%xmm3\n\t"\
		"addpd		%%xmm2		,%%xmm1\n\t"\
		"subpd		%%xmm3		,%%xmm0\n\t"\
		"movaps		%%xmm0		,%%xmm2\n\t"\
		"movaps		%%xmm1		,%%xmm3\n\t"\
		"addl		$0x020		,%%edx\n\t"\
		"movaps		 0x140(%%eax)	,%%xmm6\n\t"\
		"movaps		 0x140(%%eax)	,%%xmm4\n\t"\
		"unpckhpd	 0x140(%%ebx)	,%%xmm6\n\t"\
		"unpcklpd	 0x140(%%ebx)	,%%xmm4\n\t"\
		"movaps		%%xmm6		, 0x220(%%ecx)\n\t"\
		"movaps		 0x150(%%eax)	,%%xmm7\n\t"\
		"movaps		 0x150(%%eax)	,%%xmm5\n\t"\
		"unpckhpd	 0x150(%%ebx)	,%%xmm7\n\t"\
		"unpcklpd	 0x150(%%ebx)	,%%xmm5\n\t"\
		"movaps		%%xmm7		, 0x230(%%ecx)\n\t"\
		"movaps		%%xmm4		,%%xmm6\n\t"\
		"movaps		%%xmm5		,%%xmm7\n\t"\
		"mulpd		      (%%edx)	,%%xmm4\n\t"\
		"mulpd		      (%%edx)	,%%xmm5\n\t"\
		"mulpd		 0x010(%%edx)	,%%xmm6\n\t"\
		"mulpd		 0x010(%%edx)	,%%xmm7\n\t"\
		"addpd		%%xmm6		,%%xmm5\n\t"\
		"subpd		%%xmm7		,%%xmm4\n\t"\
		"addpd		%%xmm4		,%%xmm0\n\t"\
		"addpd		%%xmm5		,%%xmm1\n\t"\
		"subpd		%%xmm4		,%%xmm2\n\t"\
		"subpd		%%xmm5		,%%xmm3\n\t"\
		"addl		$0x040		,%%edx\n\t"\
		"movaps		 0x1c0(%%eax)	,%%xmm6\n\t"\
		"movaps		 0x1c0(%%eax)	,%%xmm4\n\t"\
		"unpckhpd	 0x1c0(%%ebx)	,%%xmm6\n\t"\
		"unpcklpd	 0x1c0(%%ebx)	,%%xmm4\n\t"\
		"movaps		%%xmm6		, 0x260(%%ecx)\n\t"\
		"movaps		 0x1d0(%%eax)	,%%xmm7\n\t"\
		"movaps		 0x1d0(%%eax)	,%%xmm5\n\t"\
		"unpckhpd	 0x1d0(%%ebx)	,%%xmm7\n\t"\
		"unpcklpd	 0x1d0(%%ebx)	,%%xmm5\n\t"\
		"movaps		%%xmm7		, 0x270(%%ecx)\n\t"\
		"movaps		%%xmm4		,%%xmm6\n\t"\
		"movaps		%%xmm5		,%%xmm7\n\t"\
		"mulpd		      (%%edx)	,%%xmm4\n\t"\
		"mulpd		      (%%edx)	,%%xmm5\n\t"\
		"mulpd		 0x010(%%edx)	,%%xmm6\n\t"\
		"mulpd		 0x010(%%edx)	,%%xmm7\n\t"\
		"addpd		%%xmm6		,%%xmm5\n\t"\
		"subpd		%%xmm7		,%%xmm4\n\t"\
		"movaps		%%xmm5		, 0x010(%%ecx)\n\t"\
		"movaps		%%xmm4		,      (%%ecx)\n\t"\
		"subl		$0x020		,%%edx\n\t"\
		"movaps		 0x0c0(%%eax)	,%%xmm6\n\t"\
		"movaps		 0x0c0(%%eax)	,%%xmm4\n\t"\
		"unpckhpd	 0x0c0(%%ebx)	,%%xmm6\n\t"\
		"unpcklpd	 0x0c0(%%ebx)	,%%xmm4\n\t"\
		"movaps		%%xmm6		, 0x240(%%ecx)\n\t"\
		"movaps		 0x0d0(%%eax)	,%%xmm7\n\t"\
		"movaps		 0x0d0(%%eax)	,%%xmm5\n\t"\
		"unpckhpd	 0x0d0(%%ebx)	,%%xmm7\n\t"\
		"unpcklpd	 0x0d0(%%ebx)	,%%xmm5\n\t"\
		"movaps		%%xmm7		, 0x250(%%ecx)\n\t"\
		"movaps		%%xmm4		,%%xmm6\n\t"\
		"movaps		%%xmm5		,%%xmm7\n\t"\
		"mulpd		      (%%edx)	,%%xmm4\n\t"\
		"mulpd		      (%%edx)	,%%xmm5\n\t"\
		"mulpd		 0x010(%%edx)	,%%xmm6\n\t"\
		"mulpd		 0x010(%%edx)	,%%xmm7\n\t"\
		"addpd		%%xmm6		,%%xmm5\n\t"\
		"subpd		%%xmm7		,%%xmm4\n\t"\
		"movaps		%%xmm5		,%%xmm7\n\t"\
		"movaps		%%xmm4		,%%xmm6\n\t"\
		"subpd		      (%%ecx)	,%%xmm4\n\t"\
		"subpd		 0x010(%%ecx)	,%%xmm5\n\t"\
		"addpd		      (%%ecx)	,%%xmm6\n\t"\
		"addpd		 0x010(%%ecx)	,%%xmm7\n\t"\
		"subpd		%%xmm6		,%%xmm0\n\t"\
		"subpd		%%xmm7		,%%xmm1\n\t"\
		"movaps		%%xmm0		, 0x040(%%ecx)\n\t"\
		"movaps		%%xmm1		, 0x050(%%ecx)\n\t"\
		"addpd		%%xmm6		,%%xmm6\n\t"\
		"addpd		%%xmm7		,%%xmm7\n\t"\
		"addpd		%%xmm0		,%%xmm6\n\t"\
		"addpd		%%xmm1		,%%xmm7\n\t"\
		"movaps		%%xmm6		,      (%%ecx)\n\t"\
		"movaps		%%xmm7		, 0x010(%%ecx)\n\t"\
		"movl		%[__isrt2]		,%%edx\n\t"\
		"subpd		%%xmm5		,%%xmm2\n\t"\
		"subpd		%%xmm4		,%%xmm3\n\t"\
		"addpd		%%xmm5		,%%xmm5\n\t"\
		"addpd		%%xmm4		,%%xmm4\n\t"\
		"addpd		%%xmm2		,%%xmm5\n\t"\
		"addpd		%%xmm3		,%%xmm4\n\t"\
		"movaps		%%xmm2		,%%xmm6\n\t"\
		"movaps		%%xmm5		,%%xmm7\n\t"\
		"subpd		%%xmm4		,%%xmm2\n\t"\
		"subpd		%%xmm3		,%%xmm5\n\t"\
		"addpd		%%xmm4		,%%xmm6\n\t"\
		"addpd		%%xmm3		,%%xmm7\n\t"\
		"mulpd		      (%%edx)	,%%xmm2\n\t"\
		"mulpd		      (%%edx)	,%%xmm5\n\t"\
		"mulpd		      (%%edx)	,%%xmm6\n\t"\
		"mulpd		      (%%edx)	,%%xmm7\n\t"\
		"movaps		%%xmm2		, 0x020(%%ecx)\n\t"\
		"movaps		%%xmm5		, 0x060(%%ecx)\n\t"\
		"movaps		%%xmm6		, 0x030(%%ecx)\n\t"\
		"movaps		%%xmm7		, 0x070(%%ecx)\n\t"\
		"/*****	SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r00)	*****/\n\t"\
		"movl		%[__r00]		,%%esi\n\t"\
		"movaps		      (%%esi)	,%%xmm0\n\t"\
		"movaps		 0x40 (%%esi)	,%%xmm4\n\t"\
		"movaps		 0x010(%%esi)	,%%xmm1\n\t"\
		"movaps		 0x50 (%%esi)	,%%xmm5\n\t"\
		"movaps		 0x80 (%%esi)	,%%xmm2\n\t"\
		"movaps		 0xd0 (%%esi)	,%%xmm7\n\t"\
		"movaps		 0x90 (%%esi)	,%%xmm3\n\t"\
		"movaps		 0xc0 (%%esi)	,%%xmm6\n\t"\
		"subpd		%%xmm2		,%%xmm0\n\t"\
		"subpd		%%xmm7		,%%xmm4\n\t"\
		"subpd		%%xmm3		,%%xmm1\n\t"\
		"subpd		%%xmm6		,%%xmm5\n\t"\
		"addpd		%%xmm2		,%%xmm2\n\t"\
		"addpd		%%xmm7		,%%xmm7\n\t"\
		"addpd		%%xmm3		,%%xmm3\n\t"\
		"addpd		%%xmm6		,%%xmm6\n\t"\
		"addpd		%%xmm0		,%%xmm2\n\t"\
		"addpd		%%xmm4		,%%xmm7\n\t"\
		"addpd		%%xmm1		,%%xmm3\n\t"\
		"addpd		%%xmm5		,%%xmm6\n\t"\
		"movaps		%%xmm0		, 0x080(%%esi)\n\t"\
		"movaps		%%xmm4		, 0x040(%%esi)\n\t"\
		"movaps		%%xmm1		, 0x090(%%esi)\n\t"\
		"movaps		%%xmm5		, 0x0d0(%%esi)\n\t"\
		"movaps		%%xmm2		,      (%%esi)\n\t"\
		"movaps		%%xmm7		, 0x0c0(%%esi)\n\t"\
		"movaps		%%xmm3		, 0x010(%%esi)\n\t"\
		"movaps		%%xmm6		, 0x050(%%esi)\n\t"\
		"movaps		 0x20 (%%esi)	,%%xmm0\n\t"\
		"movaps		 0x60 (%%esi)	,%%xmm4\n\t"\
		"movaps		 0x30 (%%esi)	,%%xmm1\n\t"\
		"movaps		 0x70 (%%esi)	,%%xmm5\n\t"\
		"movaps		 0xa0 (%%esi)	,%%xmm2\n\t"\
		"movaps		 0xf0 (%%esi)	,%%xmm7\n\t"\
		"movaps		 0xb0 (%%esi)	,%%xmm3\n\t"\
		"movaps		 0xe0 (%%esi)	,%%xmm6\n\t"\
		"subpd		%%xmm2		,%%xmm0\n\t"\
		"subpd		%%xmm7		,%%xmm4\n\t"\
		"subpd		%%xmm3		,%%xmm1\n\t"\
		"subpd		%%xmm6		,%%xmm5\n\t"\
		"addpd		%%xmm2		,%%xmm2\n\t"\
		"addpd		%%xmm7		,%%xmm7\n\t"\
		"addpd		%%xmm3		,%%xmm3\n\t"\
		"addpd		%%xmm6		,%%xmm6\n\t"\
		"addpd		%%xmm0		,%%xmm2\n\t"\
		"addpd		%%xmm4		,%%xmm7\n\t"\
		"addpd		%%xmm1		,%%xmm3\n\t"\
		"addpd		%%xmm5		,%%xmm6\n\t"\
		"movaps		%%xmm0		, 0x0a0(%%esi)\n\t"\
		"movaps		%%xmm4		, 0x060(%%esi)\n\t"\
		"movaps		%%xmm1		, 0x0b0(%%esi)\n\t"\
		"movaps		%%xmm5		, 0x0f0(%%esi)\n\t"\
		"movaps		%%xmm2		, 0x020(%%esi)\n\t"\
		"movaps		%%xmm7		, 0x0e0(%%esi)\n\t"\
		"movaps		%%xmm3		, 0x030(%%esi)\n\t"\
		"movaps		%%xmm6		, 0x070(%%esi)\n\t"\
		"/*...Block 2:	*/\n\t"\
		"addl		$0x20		,%%eax\n\t"\
		"addl		$0x20		,%%ebx\n\t"\
		"/*****	SSE2_RADIX4_DIF_4WRAPPER(c02,c0A,c12,c1A,r10)	*****/\n\t"\
		"movl		%[__r10]			,%%ecx\n\t"\
		"movl		%[__c02]			,%%edx\n\t"\
		"movaps		      (%%eax)	,%%xmm6\n\t"\
		"movaps		      (%%eax)	,%%xmm0\n\t"\
		"unpckhpd	      (%%ebx)	,%%xmm6\n\t"\
		"unpcklpd	      (%%ebx)	,%%xmm0\n\t"\
		"movaps		%%xmm6		, 0x200(%%ecx)\n\t"\
		"movaps		 0x010(%%eax)	,%%xmm7\n\t"\
		"movaps		 0x010(%%eax)	,%%xmm1\n\t"\
		"unpckhpd	 0x010(%%ebx)	,%%xmm7\n\t"\
		"unpcklpd	 0x010(%%ebx)	,%%xmm1\n\t"\
		"movaps		%%xmm7		, 0x210(%%ecx)\n\t"\
		"movaps		%%xmm0		,%%xmm2\n\t"\
		"movaps		%%xmm1		,%%xmm3\n\t"\
		"mulpd		      (%%edx)	,%%xmm0\n\t"\
		"mulpd		      (%%edx)	,%%xmm1\n\t"\
		"mulpd		 0x010(%%edx)	,%%xmm2\n\t"\
		"mulpd		 0x010(%%edx)	,%%xmm3\n\t"\
		"addpd		%%xmm2		,%%xmm1\n\t"\
		"subpd		%%xmm3		,%%xmm0\n\t"\
		"movaps		%%xmm0		,%%xmm2\n\t"\
		"movaps		%%xmm1		,%%xmm3\n\t"\
		"addl		$0x020		,%%edx\n\t"\
		"movaps		 0x100(%%eax)	,%%xmm6\n\t"\
		"movaps		 0x100(%%eax)	,%%xmm4\n\t"\
		"unpckhpd	 0x100(%%ebx)	,%%xmm6\n\t"\
		"unpcklpd	 0x100(%%ebx)	,%%xmm4\n\t"\
		"movaps		%%xmm6		, 0x220(%%ecx)\n\t"\
		"movaps		 0x110(%%eax)	,%%xmm7\n\t"\
		"movaps		 0x110(%%eax)	,%%xmm5\n\t"\
		"unpckhpd	 0x110(%%ebx)	,%%xmm7\n\t"\
		"unpcklpd	 0x110(%%ebx)	,%%xmm5\n\t"\
		"movaps		%%xmm7		, 0x230(%%ecx)\n\t"\
		"movaps		%%xmm4		,%%xmm6\n\t"\
		"movaps		%%xmm5		,%%xmm7\n\t"\
		"mulpd		      (%%edx)	,%%xmm4\n\t"\
		"mulpd		      (%%edx)	,%%xmm5\n\t"\
		"mulpd		 0x010(%%edx)	,%%xmm6\n\t"\
		"mulpd		 0x010(%%edx)	,%%xmm7\n\t"\
		"addpd		%%xmm6		,%%xmm5\n\t"\
		"subpd		%%xmm7		,%%xmm4\n\t"\
		"addpd		%%xmm4		,%%xmm0\n\t"\
		"addpd		%%xmm5		,%%xmm1\n\t"\
		"subpd		%%xmm4		,%%xmm2\n\t"\
		"subpd		%%xmm5		,%%xmm3\n\t"\
		"addl		$0x040		,%%edx\n\t"\
		"movaps		 0x180(%%eax)	,%%xmm6\n\t"\
		"movaps		 0x180(%%eax)	,%%xmm4\n\t"\
		"unpckhpd	 0x180(%%ebx)	,%%xmm6\n\t"\
		"unpcklpd	 0x180(%%ebx)	,%%xmm4\n\t"\
		"movaps		%%xmm6		, 0x260(%%ecx)\n\t"\
		"movaps		 0x190(%%eax)	,%%xmm7\n\t"\
		"movaps		 0x190(%%eax)	,%%xmm5\n\t"\
		"unpckhpd	 0x190(%%ebx)	,%%xmm7\n\t"\
		"unpcklpd	 0x190(%%ebx)	,%%xmm5\n\t"\
		"movaps		%%xmm7		, 0x270(%%ecx)\n\t"\
		"movaps		%%xmm4		,%%xmm6\n\t"\
		"movaps		%%xmm5		,%%xmm7\n\t"\
		"mulpd		      (%%edx)	,%%xmm4\n\t"\
		"mulpd		      (%%edx)	,%%xmm5\n\t"\
		"mulpd		 0x010(%%edx)	,%%xmm6\n\t"\
		"mulpd		 0x010(%%edx)	,%%xmm7\n\t"\
		"addpd		%%xmm6		,%%xmm5\n\t"\
		"subpd		%%xmm7		,%%xmm4\n\t"\
		"movaps		%%xmm5		, 0x010(%%ecx)\n\t"\
		"movaps		%%xmm4		,      (%%ecx)\n\t"\
		"subl		$0x020		,%%edx\n\t"\
		"movaps		 0x080(%%eax)	,%%xmm6\n\t"\
		"movaps		 0x080(%%eax)	,%%xmm4\n\t"\
		"unpckhpd	 0x080(%%ebx)	,%%xmm6\n\t"\
		"unpcklpd	 0x080(%%ebx)	,%%xmm4\n\t"\
		"movaps		%%xmm6		, 0x240(%%ecx)\n\t"\
		"movaps		 0x090(%%eax)	,%%xmm7\n\t"\
		"movaps		 0x090(%%eax)	,%%xmm5\n\t"\
		"unpckhpd	 0x090(%%ebx)	,%%xmm7\n\t"\
		"unpcklpd	 0x090(%%ebx)	,%%xmm5\n\t"\
		"movaps		%%xmm7		, 0x250(%%ecx)\n\t"\
		"movaps		%%xmm4		,%%xmm6\n\t"\
		"movaps		%%xmm5		,%%xmm7\n\t"\
		"mulpd		      (%%edx)	,%%xmm4\n\t"\
		"mulpd		      (%%edx)	,%%xmm5\n\t"\
		"mulpd		 0x010(%%edx)	,%%xmm6\n\t"\
		"mulpd		 0x010(%%edx)	,%%xmm7\n\t"\
		"addpd		%%xmm6		,%%xmm5\n\t"\
		"subpd		%%xmm7		,%%xmm4\n\t"\
		"movaps		%%xmm5		,%%xmm7\n\t"\
		"movaps		%%xmm4		,%%xmm6\n\t"\
		"subpd		      (%%ecx)	,%%xmm4\n\t"\
		"subpd		 0x010(%%ecx)	,%%xmm5\n\t"\
		"addpd		      (%%ecx)	,%%xmm6\n\t"\
		"addpd		 0x010(%%ecx)	,%%xmm7\n\t"\
		"subpd		%%xmm6		,%%xmm0\n\t"\
		"subpd		%%xmm7		,%%xmm1\n\t"\
		"movaps		%%xmm0		, 0x040(%%ecx)\n\t"\
		"movaps		%%xmm1		, 0x050(%%ecx)\n\t"\
		"addpd		%%xmm6		,%%xmm6\n\t"\
		"addpd		%%xmm7		,%%xmm7\n\t"\
		"addpd		%%xmm0		,%%xmm6\n\t"\
		"addpd		%%xmm1		,%%xmm7\n\t"\
		"movaps		%%xmm6		,      (%%ecx)\n\t"\
		"movaps		%%xmm7		, 0x010(%%ecx)\n\t"\
		"subpd		%%xmm5		,%%xmm2\n\t"\
		"subpd		%%xmm4		,%%xmm3\n\t"\
		"movaps		%%xmm2		, 0x020(%%ecx)\n\t"\
		"movaps		%%xmm3		, 0x070(%%ecx)\n\t"\
		"addpd		%%xmm5		,%%xmm5\n\t"\
		"addpd		%%xmm4		,%%xmm4\n\t"\
		"addpd		%%xmm2		,%%xmm5\n\t"\
		"addpd		%%xmm3		,%%xmm4\n\t"\
		"movaps		%%xmm5		, 0x060(%%ecx)\n\t"\
		"movaps		%%xmm4		, 0x030(%%ecx)\n\t"\
		"/*****	SSE2_RADIX4_DIF_4WRAPPER_2NDOFTWO(c06,c0E,c16,c1E,r18)	*****/\n\t"\
		"addl		$0x080		,%%ecx\n\t"\
		"addl		$0x040		,%%edx\n\t"\
		"movaps		 0x40 (%%eax)	,%%xmm6\n\t"\
		"movaps		 0x40 (%%eax)	,%%xmm0\n\t"\
		"unpckhpd	 0x40 (%%ebx)	,%%xmm6\n\t"\
		"unpcklpd	 0x40 (%%ebx)	,%%xmm0\n\t"\
		"movaps		%%xmm6		, 0x200(%%ecx)\n\t"\
		"movaps		 0x50 (%%eax)	,%%xmm7\n\t"\
		"movaps		 0x50 (%%eax)	,%%xmm1\n\t"\
		"unpckhpd	 0x50 (%%ebx)	,%%xmm7\n\t"\
		"unpcklpd	 0x50 (%%ebx)	,%%xmm1\n\t"\
		"movaps		%%xmm7		, 0x210(%%ecx)\n\t"\
		"movaps		%%xmm0		,%%xmm2\n\t"\
		"movaps		%%xmm1		,%%xmm3\n\t"\
		"mulpd		      (%%edx)	,%%xmm0\n\t"\
		"mulpd		      (%%edx)	,%%xmm1\n\t"\
		"mulpd		 0x010(%%edx)	,%%xmm2\n\t"\
		"mulpd		 0x010(%%edx)	,%%xmm3\n\t"\
		"addpd		%%xmm2		,%%xmm1\n\t"\
		"subpd		%%xmm3		,%%xmm0\n\t"\
		"movaps		%%xmm0		,%%xmm2\n\t"\
		"movaps		%%xmm1		,%%xmm3\n\t"\
		"addl		$0x020		,%%edx\n\t"\
		"movaps		 0x140(%%eax)	,%%xmm6\n\t"\
		"movaps		 0x140(%%eax)	,%%xmm4\n\t"\
		"unpckhpd	 0x140(%%ebx)	,%%xmm6\n\t"\
		"unpcklpd	 0x140(%%ebx)	,%%xmm4\n\t"\
		"movaps		%%xmm6		, 0x220(%%ecx)\n\t"\
		"movaps		 0x150(%%eax)	,%%xmm7\n\t"\
		"movaps		 0x150(%%eax)	,%%xmm5\n\t"\
		"unpckhpd	 0x150(%%ebx)	,%%xmm7\n\t"\
		"unpcklpd	 0x150(%%ebx)	,%%xmm5\n\t"\
		"movaps		%%xmm7		, 0x230(%%ecx)\n\t"\
		"movaps		%%xmm4		,%%xmm6\n\t"\
		"movaps		%%xmm5		,%%xmm7\n\t"\
		"mulpd		      (%%edx)	,%%xmm4\n\t"\
		"mulpd		      (%%edx)	,%%xmm5\n\t"\
		"mulpd		 0x010(%%edx)	,%%xmm6\n\t"\
		"mulpd		 0x010(%%edx)	,%%xmm7\n\t"\
		"addpd		%%xmm6		,%%xmm5\n\t"\
		"subpd		%%xmm7		,%%xmm4\n\t"\
		"addpd		%%xmm4		,%%xmm0\n\t"\
		"addpd		%%xmm5		,%%xmm1\n\t"\
		"subpd		%%xmm4		,%%xmm2\n\t"\
		"subpd		%%xmm5		,%%xmm3\n\t"\
		"addl		$0x040		,%%edx\n\t"\
		"movaps		 0x1c0(%%eax)	,%%xmm6\n\t"\
		"movaps		 0x1c0(%%eax)	,%%xmm4\n\t"\
		"unpckhpd	 0x1c0(%%ebx)	,%%xmm6\n\t"\
		"unpcklpd	 0x1c0(%%ebx)	,%%xmm4\n\t"\
		"movaps		%%xmm6		, 0x260(%%ecx)\n\t"\
		"movaps		 0x1d0(%%eax)	,%%xmm7\n\t"\
		"movaps		 0x1d0(%%eax)	,%%xmm5\n\t"\
		"unpckhpd	 0x1d0(%%ebx)	,%%xmm7\n\t"\
		"unpcklpd	 0x1d0(%%ebx)	,%%xmm5\n\t"\
		"movaps		%%xmm7		, 0x270(%%ecx)\n\t"\
		"movaps		%%xmm4		,%%xmm6\n\t"\
		"movaps		%%xmm5		,%%xmm7\n\t"\
		"mulpd		      (%%edx)	,%%xmm4\n\t"\
		"mulpd		      (%%edx)	,%%xmm5\n\t"\
		"mulpd		 0x010(%%edx)	,%%xmm6\n\t"\
		"mulpd		 0x010(%%edx)	,%%xmm7\n\t"\
		"addpd		%%xmm6		,%%xmm5\n\t"\
		"subpd		%%xmm7		,%%xmm4\n\t"\
		"movaps		%%xmm5		, 0x010(%%ecx)\n\t"\
		"movaps		%%xmm4		,      (%%ecx)\n\t"\
		"subl		$0x020		,%%edx\n\t"\
		"movaps		 0x0c0(%%eax)	,%%xmm6\n\t"\
		"movaps		 0x0c0(%%eax)	,%%xmm4\n\t"\
		"unpckhpd	 0x0c0(%%ebx)	,%%xmm6\n\t"\
		"unpcklpd	 0x0c0(%%ebx)	,%%xmm4\n\t"\
		"movaps		%%xmm6		, 0x240(%%ecx)\n\t"\
		"movaps		 0x0d0(%%eax)	,%%xmm7\n\t"\
		"movaps		 0x0d0(%%eax)	,%%xmm5\n\t"\
		"unpckhpd	 0x0d0(%%ebx)	,%%xmm7\n\t"\
		"unpcklpd	 0x0d0(%%ebx)	,%%xmm5\n\t"\
		"movaps		%%xmm7		, 0x250(%%ecx)\n\t"\
		"movaps		%%xmm4		,%%xmm6\n\t"\
		"movaps		%%xmm5		,%%xmm7\n\t"\
		"mulpd		      (%%edx)	,%%xmm4\n\t"\
		"mulpd		      (%%edx)	,%%xmm5\n\t"\
		"mulpd		 0x010(%%edx)	,%%xmm6\n\t"\
		"mulpd		 0x010(%%edx)	,%%xmm7\n\t"\
		"addpd		%%xmm6		,%%xmm5\n\t"\
		"subpd		%%xmm7		,%%xmm4\n\t"\
		"movaps		%%xmm5		,%%xmm7\n\t"\
		"movaps		%%xmm4		,%%xmm6\n\t"\
		"subpd		      (%%ecx)	,%%xmm4\n\t"\
		"subpd		 0x010(%%ecx)	,%%xmm5\n\t"\
		"addpd		      (%%ecx)	,%%xmm6\n\t"\
		"addpd		 0x010(%%ecx)	,%%xmm7\n\t"\
		"subpd		%%xmm6		,%%xmm0\n\t"\
		"subpd		%%xmm7		,%%xmm1\n\t"\
		"movaps		%%xmm0		, 0x040(%%ecx)\n\t"\
		"movaps		%%xmm1		, 0x050(%%ecx)\n\t"\
		"addpd		%%xmm6		,%%xmm6\n\t"\
		"addpd		%%xmm7		,%%xmm7\n\t"\
		"addpd		%%xmm0		,%%xmm6\n\t"\
		"addpd		%%xmm1		,%%xmm7\n\t"\
		"movaps		%%xmm6		,      (%%ecx)\n\t"\
		"movaps		%%xmm7		, 0x010(%%ecx)\n\t"\
		"movl		%[__isrt2]	,%%edx\n\t"\
		"subpd		%%xmm5		,%%xmm2\n\t"\
		"subpd		%%xmm4		,%%xmm3\n\t"\
		"addpd		%%xmm5		,%%xmm5\n\t"\
		"addpd		%%xmm4		,%%xmm4\n\t"\
		"addpd		%%xmm2		,%%xmm5\n\t"\
		"addpd		%%xmm3		,%%xmm4\n\t"\
		"movaps		%%xmm2		,%%xmm6\n\t"\
		"movaps		%%xmm5		,%%xmm7\n\t"\
		"subpd		%%xmm4		,%%xmm2\n\t"\
		"subpd		%%xmm3		,%%xmm5\n\t"\
		"addpd		%%xmm4		,%%xmm6\n\t"\
		"addpd		%%xmm3		,%%xmm7\n\t"\
		"mulpd		      (%%edx)	,%%xmm2\n\t"\
		"mulpd		      (%%edx)	,%%xmm5\n\t"\
		"mulpd		      (%%edx)	,%%xmm6\n\t"\
		"mulpd		      (%%edx)	,%%xmm7\n\t"\
		"movaps		%%xmm2		, 0x020(%%ecx)\n\t"\
		"movaps		%%xmm5		, 0x060(%%ecx)\n\t"\
		"movaps		%%xmm6		, 0x030(%%ecx)\n\t"\
		"movaps		%%xmm7		, 0x070(%%ecx)\n\t"\
		"/*****	SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r10)	*****/\n\t"\
		"movl		%[__r10]		,%%esi\n\t"\
		"movaps		      (%%esi)	,%%xmm0\n\t"\
		"movaps		 0x40 (%%esi)	,%%xmm4\n\t"\
		"movaps		 0x010(%%esi)	,%%xmm1\n\t"\
		"movaps		 0x50 (%%esi)	,%%xmm5\n\t"\
		"movaps		 0x80 (%%esi)	,%%xmm2\n\t"\
		"movaps		 0xd0 (%%esi)	,%%xmm7\n\t"\
		"movaps		 0x90 (%%esi)	,%%xmm3\n\t"\
		"movaps		 0xc0 (%%esi)	,%%xmm6\n\t"\
		"subpd		%%xmm2		,%%xmm0\n\t"\
		"subpd		%%xmm7		,%%xmm4\n\t"\
		"subpd		%%xmm3		,%%xmm1\n\t"\
		"subpd		%%xmm6		,%%xmm5\n\t"\
		"addpd		%%xmm2		,%%xmm2\n\t"\
		"addpd		%%xmm7		,%%xmm7\n\t"\
		"addpd		%%xmm3		,%%xmm3\n\t"\
		"addpd		%%xmm6		,%%xmm6\n\t"\
		"addpd		%%xmm0		,%%xmm2\n\t"\
		"addpd		%%xmm4		,%%xmm7\n\t"\
		"addpd		%%xmm1		,%%xmm3\n\t"\
		"addpd		%%xmm5		,%%xmm6\n\t"\
		"movaps		%%xmm0		, 0x080(%%esi)\n\t"\
		"movaps		%%xmm4		, 0x040(%%esi)\n\t"\
		"movaps		%%xmm1		, 0x090(%%esi)\n\t"\
		"movaps		%%xmm5		, 0x0d0(%%esi)\n\t"\
		"movaps		%%xmm2		,      (%%esi)\n\t"\
		"movaps		%%xmm7		, 0x0c0(%%esi)\n\t"\
		"movaps		%%xmm3		, 0x010(%%esi)\n\t"\
		"movaps		%%xmm6		, 0x050(%%esi)\n\t"\
		"movaps		 0x20 (%%esi)	,%%xmm0\n\t"\
		"movaps		 0x60 (%%esi)	,%%xmm4\n\t"\
		"movaps		 0x30 (%%esi)	,%%xmm1\n\t"\
		"movaps		 0x70 (%%esi)	,%%xmm5\n\t"\
		"movaps		 0xa0 (%%esi)	,%%xmm2\n\t"\
		"movaps		 0xf0 (%%esi)	,%%xmm7\n\t"\
		"movaps		 0xb0 (%%esi)	,%%xmm3\n\t"\
		"movaps		 0xe0 (%%esi)	,%%xmm6\n\t"\
		"subpd		%%xmm2		,%%xmm0\n\t"\
		"subpd		%%xmm7		,%%xmm4\n\t"\
		"subpd		%%xmm3		,%%xmm1\n\t"\
		"subpd		%%xmm6		,%%xmm5\n\t"\
		"addpd		%%xmm2		,%%xmm2\n\t"\
		"addpd		%%xmm7		,%%xmm7\n\t"\
		"addpd		%%xmm3		,%%xmm3\n\t"\
		"addpd		%%xmm6		,%%xmm6\n\t"\
		"addpd		%%xmm0		,%%xmm2\n\t"\
		"addpd		%%xmm4		,%%xmm7\n\t"\
		"addpd		%%xmm1		,%%xmm3\n\t"\
		"addpd		%%xmm5		,%%xmm6\n\t"\
		"movaps		%%xmm0		, 0x0a0(%%esi)\n\t"\
		"movaps		%%xmm4		, 0x060(%%esi)\n\t"\
		"movaps		%%xmm1		, 0x0b0(%%esi)\n\t"\
		"movaps		%%xmm5		, 0x0f0(%%esi)\n\t"\
		"movaps		%%xmm2		, 0x020(%%esi)\n\t"\
		"movaps		%%xmm7		, 0x0e0(%%esi)\n\t"\
		"movaps		%%xmm3		, 0x030(%%esi)\n\t"\
		"movaps		%%xmm6		, 0x070(%%esi)\n\t"\
		"/******************************************************************************/\n\t"\
		"/* Each of the above 2 radix-4 blocks is unique, but the ensuing 2 blocks     */\n\t"\
		"/* (operating on the odd-indexed elements from the unpck*pd commands which we */\n\t"\
		"/******************************************************************************/\n\t"\
		"/*...Block 3:	*/\n\t"\
		"/*	SSE2_RADIX4_DIF_4TWIDDLE         (r20,r24,r22,r26,r20,c01) */\n\t"\
		"addl		$0x080		,%%ecx\n\t"\
		"movl		%[__c01]		,%%ebx\n\t"\
		"movl		%%ecx		,%%eax\n\t"\
		"addl		$0x020		,%%ecx\n\t"\
		"movaps		      (%%eax)	,%%xmm0\n\t"\
		"movaps		      (%%ecx)	,%%xmm4\n\t"\
		"movaps		 0x010(%%eax)	,%%xmm1\n\t"\
		"movaps		 0x010(%%ecx)	,%%xmm5\n\t"\
		"movaps		      (%%ebx)	,%%xmm6\n\t"\
		"movaps		 0x010(%%ebx)	,%%xmm7\n\t"\
		"movaps		%%xmm0		,%%xmm2\n\t"\
		"movaps		%%xmm1		,%%xmm3\n\t"\
		"mulpd		%%xmm6		,%%xmm0\n\t"\
		"mulpd		%%xmm6		,%%xmm1\n\t"\
		"mulpd		%%xmm7		,%%xmm2\n\t"\
		"mulpd		%%xmm7		,%%xmm3\n\t"\
		"movaps		%%xmm4		,%%xmm6\n\t"\
		"addpd		%%xmm2		,%%xmm1\n\t"\
		"movaps		%%xmm5		,%%xmm7\n\t"\
		"mulpd		 0x20 (%%ebx)	,%%xmm4\n\t"\
		"subpd		%%xmm3		,%%xmm0\n\t"\
		"mulpd		 0x20 (%%ebx)	,%%xmm5\n\t"\
		"mulpd		 0x30 (%%ebx)	,%%xmm6\n\t"\
		"movaps		%%xmm0		,%%xmm2\n\t"\
		"mulpd		 0x30 (%%ebx)	,%%xmm7\n\t"\
		"addpd		%%xmm6		,%%xmm5\n\t"\
		"movaps		%%xmm1		,%%xmm3\n\t"\
		"subpd		%%xmm7		,%%xmm4\n\t"\
		"addl		$0x040		,%%ecx\n\t"\
		"addl		$0x60		,%%ebx\n\t"\
		"movaps		      (%%ecx)	,%%xmm6\n\t"\
		"movaps		 0x010(%%ecx)	,%%xmm7\n\t"\
		"addpd		%%xmm4		,%%xmm0\n\t"\
		"addpd		%%xmm5		,%%xmm1\n\t"\
		"subpd		%%xmm4		,%%xmm2\n\t"\
		"subpd		%%xmm5		,%%xmm3\n\t"\
		"movaps		%%xmm6		,%%xmm4\n\t"\
		"movaps		%%xmm7		,%%xmm5\n\t"\
		"mulpd		      (%%ebx)	,%%xmm4\n\t"\
		"mulpd		      (%%ebx)	,%%xmm5\n\t"\
		"mulpd		 0x010(%%ebx)	,%%xmm6\n\t"\
		"mulpd		 0x010(%%ebx)	,%%xmm7\n\t"\
		"movl		%%eax		,%%edx\n\t"\
		"addpd		%%xmm6		,%%xmm5\n\t"\
		"subpd		%%xmm7		,%%xmm4\n\t"\
		"movaps		%%xmm5		, 0x010(%%edx)\n\t"\
		"movaps		%%xmm4		,      (%%edx)\n\t"\
		"addl		$0x040		,%%eax\n\t"\
		"subl		$0x20		,%%ebx\n\t"\
		"movaps		      (%%eax)	,%%xmm4\n\t"\
		"movaps		 0x010(%%eax)	,%%xmm5\n\t"\
		"movaps		%%xmm4		,%%xmm6\n\t"\
		"movaps		%%xmm5		,%%xmm7\n\t"\
		"mulpd		      (%%ebx)	,%%xmm4\n\t"\
		"mulpd		      (%%ebx)	,%%xmm5\n\t"\
		"mulpd		 0x010(%%ebx)	,%%xmm6\n\t"\
		"mulpd		 0x010(%%ebx)	,%%xmm7\n\t"\
		"addpd		%%xmm6		,%%xmm5\n\t"\
		"subpd		%%xmm7		,%%xmm4\n\t"\
		"movaps		%%xmm5		,%%xmm7\n\t"\
		"movaps		%%xmm4		,%%xmm6\n\t"\
		"subpd		      (%%edx)	,%%xmm4\n\t"\
		"subpd		 0x010(%%edx)	,%%xmm5\n\t"\
		"addpd		      (%%edx)	,%%xmm6\n\t"\
		"addpd		 0x010(%%edx)	,%%xmm7\n\t"\
		"subpd		%%xmm6		,%%xmm0\n\t"\
		"subpd		%%xmm5		,%%xmm2\n\t"\
		"subpd		%%xmm7		,%%xmm1\n\t"\
		"subpd		%%xmm4		,%%xmm3\n\t"\
		"movaps		%%xmm0		, 0x040(%%edx)\n\t"\
		"movaps		%%xmm2		, 0x020(%%edx)\n\t"\
		"movaps		%%xmm1		, 0x050(%%edx)\n\t"\
		"movaps		%%xmm3		, 0x070(%%edx)\n\t"\
		"addpd		%%xmm6		,%%xmm6\n\t"\
		"addpd		%%xmm5		,%%xmm5\n\t"\
		"addpd		%%xmm7		,%%xmm7\n\t"\
		"addpd		%%xmm4		,%%xmm4\n\t"\
		"addpd		%%xmm0		,%%xmm6\n\t"\
		"addpd		%%xmm2		,%%xmm5\n\t"\
		"addpd		%%xmm1		,%%xmm7\n\t"\
		"addpd		%%xmm3		,%%xmm4\n\t"\
		"movaps		%%xmm6		,      (%%edx)\n\t"\
		"movaps		%%xmm5		, 0x060(%%edx)\n\t"\
		"movaps		%%xmm7		, 0x010(%%edx)\n\t"\
		"movaps		%%xmm4		, 0x030(%%edx)\n\t"\
		"/*****	SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO(r28,r2C,r2A,r2E,r28,c05)	*****/\n\t"\
		"movl		%[__c05]		,%%ebx\n\t"\
		"addl		$0x040		,%%eax\n\t"\
		"addl		$0x040		,%%ecx\n\t"\
		"movl		%%eax		,%%edx\n\t"\
		"movaps		      (%%eax)	,%%xmm0\n\t"\
		"movaps		      (%%ecx)	,%%xmm4\n\t"\
		"movaps		 0x010(%%eax)	,%%xmm1\n\t"\
		"movaps		 0x010(%%ecx)	,%%xmm5\n\t"\
		"movaps		      (%%ebx)	,%%xmm6\n\t"\
		"movaps		 0x010(%%ebx)	,%%xmm7\n\t"\
		"movaps		%%xmm0		,%%xmm2\n\t"\
		"movaps		%%xmm1		,%%xmm3\n\t"\
		"mulpd		%%xmm6		,%%xmm0\n\t"\
		"mulpd		%%xmm6		,%%xmm1\n\t"\
		"mulpd		%%xmm7		,%%xmm2\n\t"\
		"mulpd		%%xmm7		,%%xmm3\n\t"\
		"movaps		%%xmm4		,%%xmm6\n\t"\
		"addpd		%%xmm2		,%%xmm1\n\t"\
		"movaps		%%xmm5		,%%xmm7\n\t"\
		"mulpd		 0x20 (%%ebx)	,%%xmm4\n\t"\
		"subpd		%%xmm3		,%%xmm0\n\t"\
		"mulpd		 0x20 (%%ebx)	,%%xmm5\n\t"\
		"mulpd		 0x30 (%%ebx)	,%%xmm6\n\t"\
		"movaps		%%xmm0		,%%xmm2\n\t"\
		"mulpd		 0x30 (%%ebx)	,%%xmm7\n\t"\
		"addpd		%%xmm6		,%%xmm5\n\t"\
		"movaps		%%xmm1		,%%xmm3\n\t"\
		"subpd		%%xmm7		,%%xmm4\n\t"\
		"addpd		%%xmm4		,%%xmm0\n\t"\
		"addpd		%%xmm5		,%%xmm1\n\t"\
		"subpd		%%xmm4		,%%xmm2\n\t"\
		"subpd		%%xmm5		,%%xmm3\n\t"\
		"addl		$0x60		,%%ebx\n\t"\
		"addl		$0x040		,%%eax\n\t"\
		"addl		$0x040		,%%ecx\n\t"\
		"movaps		      (%%ecx)	,%%xmm4\n\t"\
		"movaps		 0x010(%%ecx)	,%%xmm5\n\t"\
		"movaps		      (%%ecx)	,%%xmm6\n\t"\
		"movaps		 0x010(%%ecx)	,%%xmm7\n\t"\
		"mulpd		      (%%ebx)	,%%xmm4\n\t"\
		"mulpd		      (%%ebx)	,%%xmm5\n\t"\
		"mulpd		 0x010(%%ebx)	,%%xmm6\n\t"\
		"mulpd		 0x010(%%ebx)	,%%xmm7\n\t"\
		"addpd		%%xmm6		,%%xmm5\n\t"\
		"subpd		%%xmm7		,%%xmm4\n\t"\
		"movaps		%%xmm5		, 0x010(%%edx)\n\t"\
		"movaps		%%xmm4		,      (%%edx)\n\t"\
		"subl		$0x20		,%%ebx\n\t"\
		"movaps		      (%%eax)	,%%xmm4\n\t"\
		"movaps		 0x010(%%eax)	,%%xmm5\n\t"\
		"movaps		      (%%eax)	,%%xmm6\n\t"\
		"movaps		 0x010(%%eax)	,%%xmm7\n\t"\
		"mulpd		      (%%ebx)	,%%xmm4\n\t"\
		"mulpd		      (%%ebx)	,%%xmm5\n\t"\
		"mulpd		 0x010(%%ebx)	,%%xmm6\n\t"\
		"mulpd		 0x010(%%ebx)	,%%xmm7\n\t"\
		"addpd		%%xmm6		,%%xmm5\n\t"\
		"subpd		%%xmm7		,%%xmm4\n\t"\
		"movaps		%%xmm5		,%%xmm7\n\t"\
		"movaps		%%xmm4		,%%xmm6\n\t"\
		"subpd		      (%%edx)	,%%xmm4\n\t"\
		"subpd		 0x010(%%edx)	,%%xmm5\n\t"\
		"addpd		      (%%edx)	,%%xmm6\n\t"\
		"addpd		 0x010(%%edx)	,%%xmm7\n\t"\
		"subpd		%%xmm6		,%%xmm0\n\t"\
		"movl		%[__isrt2]		,%%ecx\n\t"\
		"subpd		%%xmm7		,%%xmm1\n\t"\
		"movaps		%%xmm0		, 0x040(%%edx)\n\t"\
		"subpd		%%xmm5		,%%xmm2\n\t"\
		"movaps		%%xmm1		, 0x050(%%edx)\n\t"\
		"subpd		%%xmm4		,%%xmm3\n\t"\
		"addpd		%%xmm6		,%%xmm6\n\t"\
		"addpd		%%xmm5		,%%xmm5\n\t"\
		"addpd		%%xmm7		,%%xmm7\n\t"\
		"addpd		%%xmm4		,%%xmm4\n\t"\
		"addpd		%%xmm0		,%%xmm6\n\t"\
		"addpd		%%xmm2		,%%xmm5\n\t"\
		"addpd		%%xmm1		,%%xmm7\n\t"\
		"addpd		%%xmm3		,%%xmm4\n\t"\
		"movaps		%%xmm6		,      (%%edx)\n\t"\
		"movaps		%%xmm7		, 0x010(%%edx)\n\t"\
		"movaps		      (%%ecx)	,%%xmm0\n\t"\
		"movaps		%%xmm2		,%%xmm6\n\t"\
		"movaps		%%xmm5		,%%xmm7\n\t"\
		"subpd		%%xmm4		,%%xmm2\n\t"\
		"subpd		%%xmm3		,%%xmm5\n\t"\
		"addpd		%%xmm4		,%%xmm6\n\t"\
		"addpd		%%xmm3		,%%xmm7\n\t"\
		"mulpd		%%xmm0		,%%xmm2\n\t"\
		"mulpd		%%xmm0		,%%xmm5\n\t"\
		"mulpd		%%xmm0		,%%xmm6\n\t"\
		"mulpd		%%xmm0		,%%xmm7\n\t"\
		"movaps		%%xmm2		, 0x020(%%edx)\n\t"\
		"movaps		%%xmm5		, 0x060(%%edx)\n\t"\
		"movaps		%%xmm6		, 0x030(%%edx)\n\t"\
		"movaps		%%xmm7		, 0x070(%%edx)\n\t"\
		"/*****	SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r20)	*****/\n\t"\
		"movl		%[__r20]		,%%esi\n\t"\
		"movaps		      (%%esi)	,%%xmm0\n\t"\
		"movaps		 0x40 (%%esi)	,%%xmm4\n\t"\
		"movaps		 0x010(%%esi)	,%%xmm1\n\t"\
		"movaps		 0x50 (%%esi)	,%%xmm5\n\t"\
		"movaps		 0x80 (%%esi)	,%%xmm2\n\t"\
		"movaps		 0xd0 (%%esi)	,%%xmm7\n\t"\
		"movaps		 0x90 (%%esi)	,%%xmm3\n\t"\
		"movaps		 0xc0 (%%esi)	,%%xmm6\n\t"\
		"subpd		%%xmm2		,%%xmm0\n\t"\
		"subpd		%%xmm7		,%%xmm4\n\t"\
		"subpd		%%xmm3		,%%xmm1\n\t"\
		"subpd		%%xmm6		,%%xmm5\n\t"\
		"addpd		%%xmm2		,%%xmm2\n\t"\
		"addpd		%%xmm7		,%%xmm7\n\t"\
		"addpd		%%xmm3		,%%xmm3\n\t"\
		"addpd		%%xmm6		,%%xmm6\n\t"\
		"addpd		%%xmm0		,%%xmm2\n\t"\
		"addpd		%%xmm4		,%%xmm7\n\t"\
		"addpd		%%xmm1		,%%xmm3\n\t"\
		"addpd		%%xmm5		,%%xmm6\n\t"\
		"movaps		%%xmm0		, 0x080(%%esi)\n\t"\
		"movaps		%%xmm4		, 0x040(%%esi)\n\t"\
		"movaps		%%xmm1		, 0x090(%%esi)\n\t"\
		"movaps		%%xmm5		, 0x0d0(%%esi)\n\t"\
		"movaps		%%xmm2		,      (%%esi)\n\t"\
		"movaps		%%xmm7		, 0x0c0(%%esi)\n\t"\
		"movaps		%%xmm3		, 0x010(%%esi)\n\t"\
		"movaps		%%xmm6		, 0x050(%%esi)\n\t"\
		"movaps		 0x20 (%%esi)	,%%xmm0\n\t"\
		"movaps		 0x60 (%%esi)	,%%xmm4\n\t"\
		"movaps		 0x30 (%%esi)	,%%xmm1\n\t"\
		"movaps		 0x70 (%%esi)	,%%xmm5\n\t"\
		"movaps		 0xa0 (%%esi)	,%%xmm2\n\t"\
		"movaps		 0xf0 (%%esi)	,%%xmm7\n\t"\
		"movaps		 0xb0 (%%esi)	,%%xmm3\n\t"\
		"movaps		 0xe0 (%%esi)	,%%xmm6\n\t"\
		"subpd		%%xmm2		,%%xmm0\n\t"\
		"subpd		%%xmm7		,%%xmm4\n\t"\
		"subpd		%%xmm3		,%%xmm1\n\t"\
		"subpd		%%xmm6		,%%xmm5\n\t"\
		"addpd		%%xmm2		,%%xmm2\n\t"\
		"addpd		%%xmm7		,%%xmm7\n\t"\
		"addpd		%%xmm3		,%%xmm3\n\t"\
		"addpd		%%xmm6		,%%xmm6\n\t"\
		"addpd		%%xmm0		,%%xmm2\n\t"\
		"addpd		%%xmm4		,%%xmm7\n\t"\
		"addpd		%%xmm1		,%%xmm3\n\t"\
		"addpd		%%xmm5		,%%xmm6\n\t"\
		"movaps		%%xmm0		, 0x0a0(%%esi)\n\t"\
		"movaps		%%xmm4		, 0x060(%%esi)\n\t"\
		"movaps		%%xmm1		, 0x0b0(%%esi)\n\t"\
		"movaps		%%xmm5		, 0x0f0(%%esi)\n\t"\
		"movaps		%%xmm2		, 0x020(%%esi)\n\t"\
		"movaps		%%xmm7		, 0x0e0(%%esi)\n\t"\
		"movaps		%%xmm3		, 0x030(%%esi)\n\t"\
		"movaps		%%xmm6		, 0x070(%%esi)\n\t"\
		"/*...Block 4:	*/\n\t"\
		"/*****	SSE2_RADIX4_DIF_4TWIDDLE         (r30,r34,r32,r36,r30,c03)	*****/\n\t"\
		"movl		%[__c03]		,%%ebx\n\t"\
		"movl		%[__r30]		,%%eax\n\t"\
		"movl		%%eax		,%%ecx\n\t"\
		"addl		$0x020		,%%ecx\n\t"\
		"movaps		      (%%eax)	,%%xmm0\n\t"\
		"movaps		      (%%ecx)	,%%xmm4\n\t"\
		"movaps		 0x010(%%eax)	,%%xmm1\n\t"\
		"movaps		 0x010(%%ecx)	,%%xmm5\n\t"\
		"movaps		      (%%ebx)	,%%xmm6\n\t"\
		"movaps		 0x010(%%ebx)	,%%xmm7\n\t"\
		"movaps		%%xmm0		,%%xmm2\n\t"\
		"movaps		%%xmm1		,%%xmm3\n\t"\
		"mulpd		%%xmm6		,%%xmm0\n\t"\
		"mulpd		%%xmm6		,%%xmm1\n\t"\
		"mulpd		%%xmm7		,%%xmm2\n\t"\
		"mulpd		%%xmm7		,%%xmm3\n\t"\
		"movaps		%%xmm4		,%%xmm6\n\t"\
		"addpd		%%xmm2		,%%xmm1\n\t"\
		"movaps		%%xmm5		,%%xmm7\n\t"\
		"mulpd		 0x20 (%%ebx)	,%%xmm4\n\t"\
		"subpd		%%xmm3		,%%xmm0\n\t"\
		"mulpd		 0x20 (%%ebx)	,%%xmm5\n\t"\
		"mulpd		 0x30 (%%ebx)	,%%xmm6\n\t"\
		"movaps		%%xmm0		,%%xmm2\n\t"\
		"mulpd		 0x30 (%%ebx)	,%%xmm7\n\t"\
		"addpd		%%xmm6		,%%xmm5\n\t"\
		"movaps		%%xmm1		,%%xmm3\n\t"\
		"subpd		%%xmm7		,%%xmm4\n\t"\
		"addl		$0x040		,%%ecx\n\t"\
		"addl		$0x60		,%%ebx\n\t"\
		"movaps		      (%%ecx)	,%%xmm6\n\t"\
		"movaps		 0x010(%%ecx)	,%%xmm7\n\t"\
		"addpd		%%xmm4		,%%xmm0\n\t"\
		"addpd		%%xmm5		,%%xmm1\n\t"\
		"subpd		%%xmm4		,%%xmm2\n\t"\
		"subpd		%%xmm5		,%%xmm3\n\t"\
		"movaps		%%xmm6		,%%xmm4\n\t"\
		"movaps		%%xmm7		,%%xmm5\n\t"\
		"mulpd		      (%%ebx)	,%%xmm4\n\t"\
		"mulpd		      (%%ebx)	,%%xmm5\n\t"\
		"mulpd		 0x010(%%ebx)	,%%xmm6\n\t"\
		"mulpd		 0x010(%%ebx)	,%%xmm7\n\t"\
		"movl		%%eax		,%%edx\n\t"\
		"addpd		%%xmm6		,%%xmm5\n\t"\
		"subpd		%%xmm7		,%%xmm4\n\t"\
		"movaps		%%xmm5		, 0x010(%%edx)\n\t"\
		"movaps		%%xmm4		,      (%%edx)\n\t"\
		"addl		$0x040		,%%eax\n\t"\
		"subl		$0x20		,%%ebx\n\t"\
		"movaps		      (%%eax)	,%%xmm4\n\t"\
		"movaps		 0x010(%%eax)	,%%xmm5\n\t"\
		"movaps		%%xmm4		,%%xmm6\n\t"\
		"movaps		%%xmm5		,%%xmm7\n\t"\
		"mulpd		      (%%ebx)	,%%xmm4\n\t"\
		"mulpd		      (%%ebx)	,%%xmm5\n\t"\
		"mulpd		 0x010(%%ebx)	,%%xmm6\n\t"\
		"mulpd		 0x010(%%ebx)	,%%xmm7\n\t"\
		"addpd		%%xmm6		,%%xmm5\n\t"\
		"subpd		%%xmm7		,%%xmm4\n\t"\
		"movaps		%%xmm5		,%%xmm7\n\t"\
		"movaps		%%xmm4		,%%xmm6\n\t"\
		"subpd		      (%%edx)	,%%xmm4\n\t"\
		"subpd		 0x010(%%edx)	,%%xmm5\n\t"\
		"addpd		      (%%edx)	,%%xmm6\n\t"\
		"addpd		 0x010(%%edx)	,%%xmm7\n\t"\
		"subpd		%%xmm6		,%%xmm0\n\t"\
		"subpd		%%xmm5		,%%xmm2\n\t"\
		"subpd		%%xmm7		,%%xmm1\n\t"\
		"subpd		%%xmm4		,%%xmm3\n\t"\
		"movaps		%%xmm0		, 0x040(%%edx)\n\t"\
		"movaps		%%xmm2		, 0x020(%%edx)\n\t"\
		"movaps		%%xmm1		, 0x050(%%edx)\n\t"\
		"movaps		%%xmm3		, 0x070(%%edx)\n\t"\
		"addpd		%%xmm6		,%%xmm6\n\t"\
		"addpd		%%xmm5		,%%xmm5\n\t"\
		"addpd		%%xmm7		,%%xmm7\n\t"\
		"addpd		%%xmm4		,%%xmm4\n\t"\
		"addpd		%%xmm0		,%%xmm6\n\t"\
		"addpd		%%xmm2		,%%xmm5\n\t"\
		"addpd		%%xmm1		,%%xmm7\n\t"\
		"addpd		%%xmm3		,%%xmm4\n\t"\
		"movaps		%%xmm6		,      (%%edx)\n\t"\
		"movaps		%%xmm5		, 0x060(%%edx)\n\t"\
		"movaps		%%xmm7		, 0x010(%%edx)\n\t"\
		"movaps		%%xmm4		, 0x030(%%edx)\n\t"\
		"/*****	SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO(r38,r3C,r3A,r3E,r38,c07)	*****/\n\t"\
		"movl		%[__c07]		,%%ebx\n\t"\
		"addl		$0x040		,%%eax\n\t"\
		"addl		$0x040		,%%ecx\n\t"\
		"movl		%%eax		,%%edx\n\t"\
		"movaps		      (%%eax)	,%%xmm0\n\t"\
		"movaps		      (%%ecx)	,%%xmm4\n\t"\
		"movaps		 0x010(%%eax)	,%%xmm1\n\t"\
		"movaps		 0x010(%%ecx)	,%%xmm5\n\t"\
		"movaps		      (%%ebx)	,%%xmm6\n\t"\
		"movaps		 0x010(%%ebx)	,%%xmm7\n\t"\
		"movaps		%%xmm0		,%%xmm2\n\t"\
		"movaps		%%xmm1		,%%xmm3\n\t"\
		"mulpd		%%xmm6		,%%xmm0\n\t"\
		"mulpd		%%xmm6		,%%xmm1\n\t"\
		"mulpd		%%xmm7		,%%xmm2\n\t"\
		"mulpd		%%xmm7		,%%xmm3\n\t"\
		"movaps		%%xmm4		,%%xmm6\n\t"\
		"addpd		%%xmm2		,%%xmm1\n\t"\
		"movaps		%%xmm5		,%%xmm7\n\t"\
		"mulpd		 0x20 (%%ebx)	,%%xmm4\n\t"\
		"subpd		%%xmm3		,%%xmm0\n\t"\
		"mulpd		 0x20 (%%ebx)	,%%xmm5\n\t"\
		"mulpd		 0x30 (%%ebx)	,%%xmm6\n\t"\
		"movaps		%%xmm0		,%%xmm2\n\t"\
		"mulpd		 0x30 (%%ebx)	,%%xmm7\n\t"\
		"addpd		%%xmm6		,%%xmm5\n\t"\
		"movaps		%%xmm1		,%%xmm3\n\t"\
		"subpd		%%xmm7		,%%xmm4\n\t"\
		"addpd		%%xmm4		,%%xmm0\n\t"\
		"addpd		%%xmm5		,%%xmm1\n\t"\
		"subpd		%%xmm4		,%%xmm2\n\t"\
		"subpd		%%xmm5		,%%xmm3\n\t"\
		"addl		$0x60		,%%ebx\n\t"\
		"addl		$0x040		,%%eax\n\t"\
		"addl		$0x040		,%%ecx\n\t"\
		"movaps		      (%%ecx)	,%%xmm4\n\t"\
		"movaps		 0x010(%%ecx)	,%%xmm5\n\t"\
		"movaps		      (%%ecx)	,%%xmm6\n\t"\
		"movaps		 0x010(%%ecx)	,%%xmm7\n\t"\
		"mulpd		      (%%ebx)	,%%xmm4\n\t"\
		"mulpd		      (%%ebx)	,%%xmm5\n\t"\
		"mulpd		 0x010(%%ebx)	,%%xmm6\n\t"\
		"mulpd		 0x010(%%ebx)	,%%xmm7\n\t"\
		"addpd		%%xmm6		,%%xmm5\n\t"\
		"subpd		%%xmm7		,%%xmm4\n\t"\
		"movaps		%%xmm5		, 0x010(%%edx)\n\t"\
		"movaps		%%xmm4		,      (%%edx)\n\t"\
		"subl		$0x20		,%%ebx\n\t"\
		"movaps		      (%%eax)	,%%xmm4\n\t"\
		"movaps		 0x010(%%eax)	,%%xmm5\n\t"\
		"movaps		      (%%eax)	,%%xmm6\n\t"\
		"movaps		 0x010(%%eax)	,%%xmm7\n\t"\
		"mulpd		      (%%ebx)	,%%xmm4\n\t"\
		"mulpd		      (%%ebx)	,%%xmm5\n\t"\
		"mulpd		 0x010(%%ebx)	,%%xmm6\n\t"\
		"mulpd		 0x010(%%ebx)	,%%xmm7\n\t"\
		"addpd		%%xmm6		,%%xmm5\n\t"\
		"subpd		%%xmm7		,%%xmm4\n\t"\
		"movaps		%%xmm5		,%%xmm7\n\t"\
		"movaps		%%xmm4		,%%xmm6\n\t"\
		"subpd		      (%%edx)	,%%xmm4\n\t"\
		"subpd		 0x010(%%edx)	,%%xmm5\n\t"\
		"addpd		      (%%edx)	,%%xmm6\n\t"\
		"addpd		 0x010(%%edx)	,%%xmm7\n\t"\
		"subpd		%%xmm6		,%%xmm0\n\t"\
		"movl		%[__isrt2]		,%%ecx\n\t"\
		"subpd		%%xmm7		,%%xmm1\n\t"\
		"movaps		%%xmm0		, 0x040(%%edx)\n\t"\
		"subpd		%%xmm5		,%%xmm2\n\t"\
		"movaps		%%xmm1		, 0x050(%%edx)\n\t"\
		"subpd		%%xmm4		,%%xmm3\n\t"\
		"addpd		%%xmm6		,%%xmm6\n\t"\
		"addpd		%%xmm5		,%%xmm5\n\t"\
		"addpd		%%xmm7		,%%xmm7\n\t"\
		"addpd		%%xmm4		,%%xmm4\n\t"\
		"addpd		%%xmm0		,%%xmm6\n\t"\
		"addpd		%%xmm2		,%%xmm5\n\t"\
		"addpd		%%xmm1		,%%xmm7\n\t"\
		"addpd		%%xmm3		,%%xmm4\n\t"\
		"movaps		%%xmm6		,      (%%edx)\n\t"\
		"movaps		%%xmm7		, 0x010(%%edx)\n\t"\
		"movaps		      (%%ecx)	,%%xmm0\n\t"\
		"movaps		%%xmm2		,%%xmm6\n\t"\
		"movaps		%%xmm5		,%%xmm7\n\t"\
		"subpd		%%xmm4		,%%xmm2\n\t"\
		"subpd		%%xmm3		,%%xmm5\n\t"\
		"addpd		%%xmm4		,%%xmm6\n\t"\
		"addpd		%%xmm3		,%%xmm7\n\t"\
		"mulpd		%%xmm0		,%%xmm2\n\t"\
		"mulpd		%%xmm0		,%%xmm5\n\t"\
		"mulpd		%%xmm0		,%%xmm6\n\t"\
		"mulpd		%%xmm0		,%%xmm7\n\t"\
		"movaps		%%xmm2		, 0x020(%%edx)\n\t"\
		"movaps		%%xmm5		, 0x060(%%edx)\n\t"\
		"movaps		%%xmm6		, 0x030(%%edx)\n\t"\
		"movaps		%%xmm7		, 0x070(%%edx)\n\t"\
		"/*****	SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r30)	*****/\n\t"\
		"movl		%[__r30]		,%%esi\n\t"\
		"movaps		      (%%esi)	,%%xmm0\n\t"\
		"movaps		 0x40 (%%esi)	,%%xmm4\n\t"\
		"movaps		 0x010(%%esi)	,%%xmm1\n\t"\
		"movaps		 0x50 (%%esi)	,%%xmm5\n\t"\
		"movaps		 0x80 (%%esi)	,%%xmm2\n\t"\
		"movaps		 0xd0 (%%esi)	,%%xmm7\n\t"\
		"movaps		 0x90 (%%esi)	,%%xmm3\n\t"\
		"movaps		 0xc0 (%%esi)	,%%xmm6\n\t"\
		"subpd		%%xmm2		,%%xmm0\n\t"\
		"subpd		%%xmm7		,%%xmm4\n\t"\
		"subpd		%%xmm3		,%%xmm1\n\t"\
		"subpd		%%xmm6		,%%xmm5\n\t"\
		"addpd		%%xmm2		,%%xmm2\n\t"\
		"addpd		%%xmm7		,%%xmm7\n\t"\
		"addpd		%%xmm3		,%%xmm3\n\t"\
		"addpd		%%xmm6		,%%xmm6\n\t"\
		"addpd		%%xmm0		,%%xmm2\n\t"\
		"addpd		%%xmm4		,%%xmm7\n\t"\
		"addpd		%%xmm1		,%%xmm3\n\t"\
		"addpd		%%xmm5		,%%xmm6\n\t"\
		"movaps		%%xmm0		, 0x080(%%esi)\n\t"\
		"movaps		%%xmm4		, 0x040(%%esi)\n\t"\
		"movaps		%%xmm1		, 0x090(%%esi)\n\t"\
		"movaps		%%xmm5		, 0x0d0(%%esi)\n\t"\
		"movaps		%%xmm2		,      (%%esi)\n\t"\
		"movaps		%%xmm7		, 0x0c0(%%esi)\n\t"\
		"movaps		%%xmm3		, 0x010(%%esi)\n\t"\
		"movaps		%%xmm6		, 0x050(%%esi)\n\t"\
		"movaps		 0x20 (%%esi)	,%%xmm0\n\t"\
		"movaps		 0x60 (%%esi)	,%%xmm4\n\t"\
		"movaps		 0x30 (%%esi)	,%%xmm1\n\t"\
		"movaps		 0x70 (%%esi)	,%%xmm5\n\t"\
		"movaps		 0xa0 (%%esi)	,%%xmm2\n\t"\
		"movaps		 0xf0 (%%esi)	,%%xmm7\n\t"\
		"movaps		 0xb0 (%%esi)	,%%xmm3\n\t"\
		"movaps		 0xe0 (%%esi)	,%%xmm6\n\t"\
		"subpd		%%xmm2		,%%xmm0\n\t"\
		"subpd		%%xmm7		,%%xmm4\n\t"\
		"subpd		%%xmm3		,%%xmm1\n\t"\
		"subpd		%%xmm6		,%%xmm5\n\t"\
		"addpd		%%xmm2		,%%xmm2\n\t"\
		"addpd		%%xmm7		,%%xmm7\n\t"\
		"addpd		%%xmm3		,%%xmm3\n\t"\
		"addpd		%%xmm6		,%%xmm6\n\t"\
		"addpd		%%xmm0		,%%xmm2\n\t"\
		"addpd		%%xmm4		,%%xmm7\n\t"\
		"addpd		%%xmm1		,%%xmm3\n\t"\
		"addpd		%%xmm5		,%%xmm6\n\t"\
		"movaps		%%xmm0		, 0x0a0(%%esi)\n\t"\
		"movaps		%%xmm4		, 0x060(%%esi)\n\t"\
		"movaps		%%xmm1		, 0x0b0(%%esi)\n\t"\
		"movaps		%%xmm5		, 0x0f0(%%esi)\n\t"\
		"movaps		%%xmm2		, 0x020(%%esi)\n\t"\
		"movaps		%%xmm7		, 0x0e0(%%esi)\n\t"\
		"movaps		%%xmm3		, 0x030(%%esi)\n\t"\
		"movaps		%%xmm6		, 0x070(%%esi)\n\t"\
		"/**********************************************************************************/\n\t"\
		"/*...and now do eight radix-4 transforms, including the internal twiddle factots: */\n\t"\
		"/**********************************************************************************/\n\t"\
		"/*...Block 1: t00,t10,t20,t30	*/\n\t"\
		"movl		%[__r00]		,%%eax\n\t"\
		"movl		%[__r10]		,%%ebx\n\t"\
		"movl		%[__r20]		,%%ecx\n\t"\
		"movl		%[__r30]		,%%edx\n\t"\
		"movaps		      (%%eax)	,%%xmm0\n\t"\
		"movaps		 0x010(%%eax)	,%%xmm1\n\t"\
		"movaps		      (%%ebx)	,%%xmm2\n\t"\
		"movaps		 0x010(%%ebx)	,%%xmm3\n\t"\
		"subpd		      (%%ebx)	,%%xmm0\n\t"\
		"subpd		 0x010(%%ebx)	,%%xmm1\n\t"\
		"addpd		      (%%eax)	,%%xmm2\n\t"\
		"addpd		 0x010(%%eax)	,%%xmm3\n\t"\
		"movaps		      (%%ecx)	,%%xmm4\n\t"\
		"movaps		 0x010(%%ecx)	,%%xmm5\n\t"\
		"movaps		      (%%edx)	,%%xmm6\n\t"\
		"movaps		 0x010(%%edx)	,%%xmm7\n\t"\
		"subpd		      (%%edx)	,%%xmm4\n\t"\
		"subpd		 0x010(%%edx)	,%%xmm5\n\t"\
		"addpd		      (%%ecx)	,%%xmm6\n\t"\
		"addpd		 0x010(%%ecx)	,%%xmm7\n\t"\
		"subpd		%%xmm6		,%%xmm2\n\t"\
		"subpd		%%xmm7		,%%xmm3\n\t"\
		"addpd		%%xmm6		,%%xmm6\n\t"\
		"addpd		%%xmm7		,%%xmm7\n\t"\
		"movaps		%%xmm2		,      (%%ecx)\n\t"\
		"movaps		%%xmm3		, 0x010(%%ecx)\n\t"\
		"addpd		%%xmm2		,%%xmm6\n\t"\
		"addpd		%%xmm3		,%%xmm7\n\t"\
		"movaps		%%xmm6		,      (%%eax)\n\t"\
		"movaps		%%xmm7		, 0x010(%%eax)\n\t"\
		"subpd		%%xmm5		,%%xmm0\n\t"\
		"subpd		%%xmm4		,%%xmm1\n\t"\
		"addpd		%%xmm5		,%%xmm5\n\t"\
		"addpd		%%xmm4		,%%xmm4\n\t"\
		"movaps		%%xmm0		,      (%%ebx)\n\t"\
		"movaps		%%xmm1		, 0x010(%%edx)\n\t"\
		"addpd		%%xmm0		,%%xmm5\n\t"\
		"addpd		%%xmm1		,%%xmm4\n\t"\
		"movaps		%%xmm5		,      (%%edx)\n\t"\
		"movaps		%%xmm4		, 0x010(%%ebx)\n\t"\
		"/*...Block 5: t08,t18,t28,t38	*/\n\t"\
		"addl		$0x080		,%%eax\n\t"\
		"addl		$0x080		,%%ebx\n\t"\
		"addl		$0x080		,%%ecx\n\t"\
		"addl		$0x080		,%%edx\n\t"\
		"movl		%[__isrt2]		,%%esi\n\t"\
		"movaps		      (%%eax)	,%%xmm0\n\t"\
		"movaps		 0x010(%%eax)	,%%xmm1\n\t"\
		"movaps		      (%%ebx)	,%%xmm2\n\t"\
		"movaps		 0x010(%%ebx)	,%%xmm3\n\t"\
		"subpd		 0x010(%%ebx)	,%%xmm0\n\t"\
		"subpd		      (%%ebx)	,%%xmm1\n\t"\
		"addpd		 0x010(%%eax)	,%%xmm2\n\t"\
		"addpd		      (%%eax)	,%%xmm3\n\t"\
		"movaps		      (%%ecx)	,%%xmm4\n\t"\
		"movaps		 0x010(%%ecx)	,%%xmm5\n\t"\
		"movaps		      (%%edx)	,%%xmm6\n\t"\
		"movaps		 0x010(%%edx)	,%%xmm7\n\t"\
		"subpd		 0x010(%%ecx)	,%%xmm4\n\t"\
		"addpd		      (%%ecx)	,%%xmm5\n\t"\
		"mulpd		      (%%esi)	,%%xmm4\n\t"\
		"mulpd		      (%%esi)	,%%xmm5\n\t"\
		"addpd		 0x010(%%edx)	,%%xmm6\n\t"\
		"subpd		      (%%edx)	,%%xmm7\n\t"\
		"mulpd		      (%%esi)	,%%xmm6\n\t"\
		"mulpd		      (%%esi)	,%%xmm7\n\t"\
		"subpd		%%xmm6		,%%xmm4\n\t"\
		"subpd		%%xmm7		,%%xmm5\n\t"\
		"addpd		%%xmm6		,%%xmm6\n\t"\
		"addpd		%%xmm7		,%%xmm7\n\t"\
		"addpd		%%xmm4		,%%xmm6\n\t"\
		"addpd		%%xmm5		,%%xmm7\n\t"\
		"subpd		%%xmm4		,%%xmm0\n\t"\
		"subpd		%%xmm5		,%%xmm2\n\t"\
		"addpd		%%xmm4		,%%xmm4\n\t"\
		"addpd		%%xmm5		,%%xmm5\n\t"\
		"movaps		%%xmm0		,      (%%ecx)\n\t"\
		"movaps		%%xmm2		, 0x010(%%ecx)\n\t"\
		"addpd		%%xmm0		,%%xmm4\n\t"\
		"addpd		%%xmm2		,%%xmm5\n\t"\
		"movaps		%%xmm4		,      (%%eax)\n\t"\
		"movaps		%%xmm5		, 0x010(%%eax)\n\t"\
		"subpd		%%xmm7		,%%xmm3\n\t"\
		"subpd		%%xmm6		,%%xmm1\n\t"\
		"addpd		%%xmm7		,%%xmm7\n\t"\
		"addpd		%%xmm6		,%%xmm6\n\t"\
		"movaps		%%xmm3		,      (%%ebx)\n\t"\
		"movaps		%%xmm1		, 0x010(%%edx)\n\t"\
		"addpd		%%xmm3		,%%xmm7\n\t"\
		"addpd		%%xmm1		,%%xmm6\n\t"\
		"movaps		%%xmm7		,      (%%edx)\n\t"\
		"movaps		%%xmm6		, 0x010(%%ebx)\n\t"\
		"/*...Block 3: t04,t14,t24,t34	*/\n\t"\
		"subl		$0x040		,%%eax\n\t"\
		"subl		$0x040		,%%ebx\n\t"\
		"subl		$0x040		,%%ecx\n\t"\
		"subl		$0x040		,%%edx\n\t"\
		"movl		%[__cc0]		,%%edi\n\t"\
		"movaps		      (%%ecx)	,%%xmm4\n\t"\
		"movaps		      (%%edx)	,%%xmm6\n\t"\
		"movaps		 0x010(%%ecx)	,%%xmm5\n\t"\
		"movaps		 0x010(%%edx)	,%%xmm7\n\t"\
		"movaps		      (%%ecx)	,%%xmm0\n\t"\
		"movaps		      (%%edx)	,%%xmm2\n\t"\
		"movaps		 0x010(%%ecx)	,%%xmm1\n\t"\
		"movaps		 0x010(%%edx)	,%%xmm3\n\t"\
		"mulpd		      (%%edi)	,%%xmm4\n\t"\
		"mulpd		 0x010(%%edi)	,%%xmm6\n\t"\
		"mulpd		 0x010(%%edi)	,%%xmm1\n\t"\
		"mulpd		      (%%edi)	,%%xmm3\n\t"\
		"mulpd		      (%%edi)	,%%xmm5\n\t"\
		"mulpd		 0x010(%%edi)	,%%xmm7\n\t"\
		"mulpd		 0x010(%%edi)	,%%xmm0\n\t"\
		"mulpd		      (%%edi)	,%%xmm2\n\t"\
		"subpd		%%xmm1		,%%xmm4\n\t"\
		"subpd		%%xmm3		,%%xmm6\n\t"\
		"addpd		%%xmm0		,%%xmm5\n\t"\
		"addpd		%%xmm2		,%%xmm7\n\t"\
		"subpd		%%xmm6		,%%xmm4\n\t"\
		"subpd		%%xmm7		,%%xmm5\n\t"\
		"addpd		%%xmm6		,%%xmm6\n\t"\
		"addpd		%%xmm7		,%%xmm7\n\t"\
		"addpd		%%xmm4		,%%xmm6\n\t"\
		"addpd		%%xmm5		,%%xmm7\n\t"\
		"movaps		      (%%ebx)	,%%xmm2\n\t"\
		"movaps		 0x010(%%ebx)	,%%xmm3\n\t"\
		"subpd		 0x010(%%ebx)	,%%xmm2\n\t"\
		"addpd		      (%%ebx)	,%%xmm3\n\t"\
		"mulpd		      (%%esi)	,%%xmm2\n\t"\
		"mulpd		      (%%esi)	,%%xmm3\n\t"\
		"movaps		      (%%eax)	,%%xmm0\n\t"\
		"movaps		 0x010(%%eax)	,%%xmm1\n\t"\
		"subpd		%%xmm2		,%%xmm0\n\t"\
		"subpd		%%xmm3		,%%xmm1\n\t"\
		"addpd		      (%%eax)	,%%xmm2\n\t"\
		"addpd		 0x010(%%eax)	,%%xmm3\n\t"\
		"subpd		%%xmm6		,%%xmm2\n\t"\
		"subpd		%%xmm7		,%%xmm3\n\t"\
		"addpd		%%xmm6		,%%xmm6\n\t"\
		"addpd		%%xmm7		,%%xmm7\n\t"\
		"movaps		%%xmm2		,      (%%ecx)\n\t"\
		"movaps		%%xmm3		, 0x010(%%ecx)\n\t"\
		"addpd		%%xmm2		,%%xmm6\n\t"\
		"addpd		%%xmm3		,%%xmm7\n\t"\
		"movaps		%%xmm6		,      (%%eax)\n\t"\
		"movaps		%%xmm7		, 0x010(%%eax)\n\t"\
		"subpd		%%xmm5		,%%xmm0\n\t"\
		"subpd		%%xmm4		,%%xmm1\n\t"\
		"addpd		%%xmm5		,%%xmm5\n\t"\
		"addpd		%%xmm4		,%%xmm4\n\t"\
		"movaps		%%xmm0		,      (%%ebx)\n\t"\
		"movaps		%%xmm1		, 0x010(%%edx)\n\t"\
		"addpd		%%xmm0		,%%xmm5\n\t"\
		"addpd		%%xmm1		,%%xmm4\n\t"\
		"movaps		%%xmm5		,      (%%edx)\n\t"\
		"movaps		%%xmm4		, 0x010(%%ebx)\n\t"\
		"/*...Block 7: t0C,t1C,t2C,t3C	*/\n\t"\
		"addl		$0x080		,%%eax\n\t"\
		"addl		$0x080		,%%ebx\n\t"\
		"addl		$0x080		,%%ecx\n\t"\
		"addl		$0x080		,%%edx\n\t"\
		"movaps		      (%%ecx)	,%%xmm4\n\t"\
		"movaps		      (%%edx)	,%%xmm6\n\t"\
		"movaps		 0x010(%%ecx)	,%%xmm5\n\t"\
		"movaps		 0x010(%%edx)	,%%xmm7\n\t"\
		"movaps		      (%%ecx)	,%%xmm0\n\t"\
		"movaps		      (%%edx)	,%%xmm2\n\t"\
		"movaps		 0x010(%%ecx)	,%%xmm1\n\t"\
		"movaps		 0x010(%%edx)	,%%xmm3\n\t"\
		"mulpd		 0x010(%%edi)	,%%xmm4\n\t"\
		"mulpd		      (%%edi)	,%%xmm6\n\t"\
		"mulpd		      (%%edi)	,%%xmm1\n\t"\
		"mulpd		 0x010(%%edi)	,%%xmm3\n\t"\
		"mulpd		 0x010(%%edi)	,%%xmm5\n\t"\
		"mulpd		      (%%edi)	,%%xmm7\n\t"\
		"mulpd		      (%%edi)	,%%xmm0\n\t"\
		"mulpd		 0x010(%%edi)	,%%xmm2\n\t"\
		"subpd		%%xmm1		,%%xmm4\n\t"\
		"subpd		%%xmm3		,%%xmm6\n\t"\
		"addpd		%%xmm0		,%%xmm5\n\t"\
		"addpd		%%xmm2		,%%xmm7\n\t"\
		"subpd		%%xmm6		,%%xmm4\n\t"\
		"subpd		%%xmm7		,%%xmm5\n\t"\
		"addpd		%%xmm6		,%%xmm6\n\t"\
		"addpd		%%xmm7		,%%xmm7\n\t"\
		"addpd		%%xmm4		,%%xmm6\n\t"\
		"addpd		%%xmm5		,%%xmm7\n\t"\
		"movaps		      (%%ebx)	,%%xmm2\n\t"\
		"movaps		 0x010(%%ebx)	,%%xmm3\n\t"\
		"addpd		 0x010(%%ebx)	,%%xmm2\n\t"\
		"subpd		      (%%ebx)	,%%xmm3\n\t"\
		"mulpd		      (%%esi)	,%%xmm2\n\t"\
		"mulpd		      (%%esi)	,%%xmm3\n\t"\
		"movaps		      (%%eax)	,%%xmm0\n\t"\
		"movaps		 0x010(%%eax)	,%%xmm1\n\t"\
		"subpd		%%xmm2		,%%xmm0\n\t"\
		"subpd		%%xmm3		,%%xmm1\n\t"\
		"addpd		      (%%eax)	,%%xmm2\n\t"\
		"addpd		 0x010(%%eax)	,%%xmm3\n\t"\
		"subpd		%%xmm4		,%%xmm0\n\t"\
		"subpd		%%xmm5		,%%xmm1\n\t"\
		"addpd		%%xmm4		,%%xmm4\n\t"\
		"addpd		%%xmm5		,%%xmm5\n\t"\
		"movaps		%%xmm0		,      (%%ecx)\n\t"\
		"movaps		%%xmm1		, 0x010(%%ecx)\n\t"\
		"addpd		%%xmm0		,%%xmm4\n\t"\
		"addpd		%%xmm1		,%%xmm5\n\t"\
		"movaps		%%xmm4		,      (%%eax)\n\t"\
		"movaps		%%xmm5		, 0x010(%%eax)\n\t"\
		"subpd		%%xmm7		,%%xmm2\n\t"\
		"subpd		%%xmm6		,%%xmm3\n\t"\
		"addpd		%%xmm7		,%%xmm7\n\t"\
		"addpd		%%xmm6		,%%xmm6\n\t"\
		"movaps		%%xmm2		,      (%%ebx)\n\t"\
		"movaps		%%xmm3		, 0x010(%%edx)\n\t"\
		"addpd		%%xmm2		,%%xmm7\n\t"\
		"addpd		%%xmm3		,%%xmm6\n\t"\
		"movaps		%%xmm7		,      (%%edx)\n\t"\
		"movaps		%%xmm6		, 0x010(%%ebx)\n\t"\
		"/*...Block 2: t02,t12,t22,t32	*/\n\t"\
		"subl		$0x0a0		,%%eax\n\t"\
		"subl		$0x0a0		,%%ebx\n\t"\
		"subl		$0x0a0		,%%ecx\n\t"\
		"subl		$0x0a0		,%%edx\n\t"\
		"addl		$0x30		,%%esi\n\t"\
		"movaps		      (%%ecx)	,%%xmm4\n\t"\
		"movaps		      (%%edx)	,%%xmm6\n\t"\
		"movaps		 0x010(%%ecx)	,%%xmm5\n\t"\
		"movaps		 0x010(%%edx)	,%%xmm7\n\t"\
		"movaps		      (%%ecx)	,%%xmm0\n\t"\
		"movaps		      (%%edx)	,%%xmm2\n\t"\
		"movaps		 0x010(%%ecx)	,%%xmm1\n\t"\
		"movaps		 0x010(%%edx)	,%%xmm3\n\t"\
		"mulpd		      (%%esi)	,%%xmm4\n\t"\
		"mulpd		 0x20 (%%esi)	,%%xmm6\n\t"\
		"mulpd		 0x010(%%esi)	,%%xmm1\n\t"\
		"mulpd		 0x30 (%%esi)	,%%xmm3\n\t"\
		"mulpd		      (%%esi)	,%%xmm5\n\t"\
		"mulpd		 0x20 (%%esi)	,%%xmm7\n\t"\
		"mulpd		 0x010(%%esi)	,%%xmm0\n\t"\
		"mulpd		 0x30 (%%esi)	,%%xmm2\n\t"\
		"subpd		%%xmm1		,%%xmm4\n\t"\
		"subpd		%%xmm3		,%%xmm6\n\t"\
		"addpd		%%xmm0		,%%xmm5\n\t"\
		"addpd		%%xmm2		,%%xmm7\n\t"\
		"subpd		%%xmm6		,%%xmm4\n\t"\
		"subpd		%%xmm7		,%%xmm5\n\t"\
		"addpd		%%xmm6		,%%xmm6\n\t"\
		"addpd		%%xmm7		,%%xmm7\n\t"\
		"addpd		%%xmm4		,%%xmm6\n\t"\
		"addpd		%%xmm5		,%%xmm7\n\t"\
		"movaps		      (%%ebx)	,%%xmm2\n\t"\
		"movaps		 0x010(%%ebx)	,%%xmm0\n\t"\
		"movaps		      (%%ebx)	,%%xmm1\n\t"\
		"movaps		 0x010(%%ebx)	,%%xmm3\n\t"\
		"mulpd		      (%%edi)	,%%xmm2\n\t"\
		"mulpd		 0x010(%%edi)	,%%xmm0\n\t"\
		"mulpd		      (%%edi)	,%%xmm3\n\t"\
		"mulpd		 0x010(%%edi)	,%%xmm1\n\t"\
		"subpd		%%xmm0		,%%xmm2\n\t"\
		"addpd		%%xmm1		,%%xmm3\n\t"\
		"movaps		      (%%eax)	,%%xmm0\n\t"\
		"movaps		 0x010(%%eax)	,%%xmm1\n\t"\
		"subpd		%%xmm2		,%%xmm0\n\t"\
		"subpd		%%xmm3		,%%xmm1\n\t"\
		"addpd		      (%%eax)	,%%xmm2\n\t"\
		"addpd		 0x010(%%eax)	,%%xmm3\n\t"\
		"subpd		%%xmm6		,%%xmm2\n\t"\
		"subpd		%%xmm7		,%%xmm3\n\t"\
		"addpd		%%xmm6		,%%xmm6\n\t"\
		"addpd		%%xmm7		,%%xmm7\n\t"\
		"movaps		%%xmm2		,      (%%ecx)\n\t"\
		"movaps		%%xmm3		, 0x010(%%ecx)\n\t"\
		"addpd		%%xmm2		,%%xmm6\n\t"\
		"addpd		%%xmm3		,%%xmm7\n\t"\
		"movaps		%%xmm6		,      (%%eax)\n\t"\
		"movaps		%%xmm7		, 0x010(%%eax)\n\t"\
		"subpd		%%xmm5		,%%xmm0\n\t"\
		"subpd		%%xmm4		,%%xmm1\n\t"\
		"addpd		%%xmm5		,%%xmm5\n\t"\
		"addpd		%%xmm4		,%%xmm4\n\t"\
		"movaps		%%xmm0		,      (%%ebx)\n\t"\
		"movaps		%%xmm1		, 0x010(%%edx)\n\t"\
		"addpd		%%xmm0		,%%xmm5\n\t"\
		"addpd		%%xmm1		,%%xmm4\n\t"\
		"movaps		%%xmm5		,      (%%edx)\n\t"\
		"movaps		%%xmm4		, 0x010(%%ebx)\n\t"\
		"/*...Block 6: t0A,t1A,t2A,t3A	*/\n\t"\
		"addl		$0x080		,%%eax\n\t"\
		"addl		$0x080		,%%ebx\n\t"\
		"addl		$0x080		,%%ecx\n\t"\
		"addl		$0x080		,%%edx\n\t"\
		"movaps		      (%%ecx)	,%%xmm4\n\t"\
		"movaps		      (%%edx)	,%%xmm6\n\t"\
		"movaps		 0x010(%%ecx)	,%%xmm5\n\t"\
		"movaps		 0x010(%%edx)	,%%xmm7\n\t"\
		"movaps		      (%%ecx)	,%%xmm0\n\t"\
		"movaps		      (%%edx)	,%%xmm2\n\t"\
		"movaps		 0x010(%%ecx)	,%%xmm1\n\t"\
		"movaps		 0x010(%%edx)	,%%xmm3\n\t"\
		"mulpd		 0x30 (%%esi)	,%%xmm4\n\t"\
		"mulpd		      (%%esi)	,%%xmm6\n\t"\
		"mulpd		 0x20 (%%esi)	,%%xmm1\n\t"\
		"mulpd		 0x010(%%esi)	,%%xmm3\n\t"\
		"mulpd		 0x30 (%%esi)	,%%xmm5\n\t"\
		"mulpd		      (%%esi)	,%%xmm7\n\t"\
		"mulpd		 0x20 (%%esi)	,%%xmm0\n\t"\
		"mulpd		 0x010(%%esi)	,%%xmm2\n\t"\
		"subpd		%%xmm1		,%%xmm4\n\t"\
		"addpd		%%xmm3		,%%xmm6\n\t"\
		"addpd		%%xmm0		,%%xmm5\n\t"\
		"subpd		%%xmm2		,%%xmm7\n\t"\
		"subpd		%%xmm6		,%%xmm4\n\t"\
		"subpd		%%xmm7		,%%xmm5\n\t"\
		"addpd		%%xmm6		,%%xmm6\n\t"\
		"addpd		%%xmm7		,%%xmm7\n\t"\
		"addpd		%%xmm4		,%%xmm6\n\t"\
		"addpd		%%xmm5		,%%xmm7\n\t"\
		"movaps		      (%%ebx)	,%%xmm2\n\t"\
		"movaps		 0x010(%%ebx)	,%%xmm0\n\t"\
		"movaps		      (%%ebx)	,%%xmm1\n\t"\
		"movaps		 0x010(%%ebx)	,%%xmm3\n\t"\
		"mulpd		 0x010(%%edi)	,%%xmm2\n\t"\
		"mulpd		      (%%edi)	,%%xmm0\n\t"\
		"mulpd		 0x010(%%edi)	,%%xmm3\n\t"\
		"mulpd		      (%%edi)	,%%xmm1\n\t"\
		"addpd		%%xmm0		,%%xmm2\n\t"\
		"subpd		%%xmm1		,%%xmm3\n\t"\
		"movaps		      (%%eax)	,%%xmm0\n\t"\
		"movaps		 0x010(%%eax)	,%%xmm1\n\t"\
		"subpd		%%xmm2		,%%xmm0\n\t"\
		"subpd		%%xmm3		,%%xmm1\n\t"\
		"addpd		      (%%eax)	,%%xmm2\n\t"\
		"addpd		 0x010(%%eax)	,%%xmm3\n\t"\
		"subpd		%%xmm4		,%%xmm0\n\t"\
		"subpd		%%xmm5		,%%xmm1\n\t"\
		"addpd		%%xmm4		,%%xmm4\n\t"\
		"addpd		%%xmm5		,%%xmm5\n\t"\
		"movaps		%%xmm0		,      (%%ecx)\n\t"\
		"movaps		%%xmm1		, 0x010(%%ecx)\n\t"\
		"addpd		%%xmm0		,%%xmm4\n\t"\
		"addpd		%%xmm1		,%%xmm5\n\t"\
		"movaps		%%xmm4		,      (%%eax)\n\t"\
		"movaps		%%xmm5		, 0x010(%%eax)\n\t"\
		"subpd		%%xmm7		,%%xmm2\n\t"\
		"subpd		%%xmm6		,%%xmm3\n\t"\
		"addpd		%%xmm7		,%%xmm7\n\t"\
		"addpd		%%xmm6		,%%xmm6\n\t"\
		"movaps		%%xmm2		,      (%%ebx)\n\t"\
		"movaps		%%xmm3		, 0x010(%%edx)\n\t"\
		"addpd		%%xmm2		,%%xmm7\n\t"\
		"addpd		%%xmm3		,%%xmm6\n\t"\
		"movaps		%%xmm7		,      (%%edx)\n\t"\
		"movaps		%%xmm6		, 0x010(%%ebx)\n\t"\
		"/*...Block 4: t06,t16,t26,t36	*/\n\t"\
		"subl		$0x040		,%%eax\n\t"\
		"subl		$0x040		,%%ebx\n\t"\
		"subl		$0x040		,%%ecx\n\t"\
		"subl		$0x040		,%%edx\n\t"\
		"movaps		      (%%ecx)	,%%xmm4\n\t"\
		"movaps		      (%%edx)	,%%xmm6\n\t"\
		"movaps		 0x010(%%ecx)	,%%xmm5\n\t"\
		"movaps		 0x010(%%edx)	,%%xmm7\n\t"\
		"movaps		      (%%ecx)	,%%xmm0\n\t"\
		"movaps		      (%%edx)	,%%xmm2\n\t"\
		"movaps		 0x010(%%ecx)	,%%xmm1\n\t"\
		"movaps		 0x010(%%edx)	,%%xmm3\n\t"\
		"mulpd		 0x20 (%%esi)	,%%xmm4\n\t"\
		"mulpd		 0x010(%%esi)	,%%xmm6\n\t"\
		"mulpd		 0x30 (%%esi)	,%%xmm1\n\t"\
		"mulpd		      (%%esi)	,%%xmm3\n\t"\
		"mulpd		 0x20 (%%esi)	,%%xmm5\n\t"\
		"mulpd		 0x010(%%esi)	,%%xmm7\n\t"\
		"mulpd		 0x30 (%%esi)	,%%xmm0\n\t"\
		"mulpd		      (%%esi)	,%%xmm2\n\t"\
		"subpd		%%xmm1		,%%xmm4\n\t"\
		"addpd		%%xmm3		,%%xmm6\n\t"\
		"addpd		%%xmm0		,%%xmm5\n\t"\
		"subpd		%%xmm2		,%%xmm7\n\t"\
		"subpd		%%xmm6		,%%xmm4\n\t"\
		"subpd		%%xmm7		,%%xmm5\n\t"\
		"addpd		%%xmm6		,%%xmm6\n\t"\
		"addpd		%%xmm7		,%%xmm7\n\t"\
		"addpd		%%xmm4		,%%xmm6\n\t"\
		"addpd		%%xmm5		,%%xmm7\n\t"\
		"movaps		      (%%ebx)	,%%xmm2\n\t"\
		"movaps		 0x010(%%ebx)	,%%xmm0\n\t"\
		"movaps		      (%%ebx)	,%%xmm1\n\t"\
		"movaps		 0x010(%%ebx)	,%%xmm3\n\t"\
		"mulpd		 0x010(%%edi)	,%%xmm2\n\t"\
		"mulpd		      (%%edi)	,%%xmm0\n\t"\
		"mulpd		 0x010(%%edi)	,%%xmm3\n\t"\
		"mulpd		      (%%edi)	,%%xmm1\n\t"\
		"subpd		%%xmm0		,%%xmm2\n\t"\
		"addpd		%%xmm1		,%%xmm3\n\t"\
		"movaps		      (%%eax)	,%%xmm0\n\t"\
		"movaps		 0x010(%%eax)	,%%xmm1\n\t"\
		"subpd		%%xmm2		,%%xmm0\n\t"\
		"subpd		%%xmm3		,%%xmm1\n\t"\
		"addpd		      (%%eax)	,%%xmm2\n\t"\
		"addpd		 0x010(%%eax)	,%%xmm3\n\t"\
		"subpd		%%xmm4		,%%xmm2\n\t"\
		"subpd		%%xmm5		,%%xmm3\n\t"\
		"addpd		%%xmm4		,%%xmm4\n\t"\
		"addpd		%%xmm5		,%%xmm5\n\t"\
		"movaps		%%xmm2		,      (%%ecx)\n\t"\
		"movaps		%%xmm3		, 0x010(%%ecx)\n\t"\
		"addpd		%%xmm2		,%%xmm4\n\t"\
		"addpd		%%xmm3		,%%xmm5\n\t"\
		"movaps		%%xmm4		,      (%%eax)\n\t"\
		"movaps		%%xmm5		, 0x010(%%eax)\n\t"\
		"subpd		%%xmm7		,%%xmm0\n\t"\
		"subpd		%%xmm6		,%%xmm1\n\t"\
		"addpd		%%xmm7		,%%xmm7\n\t"\
		"addpd		%%xmm6		,%%xmm6\n\t"\
		"movaps		%%xmm0		,      (%%ebx)\n\t"\
		"movaps		%%xmm1		, 0x010(%%edx)\n\t"\
		"addpd		%%xmm0		,%%xmm7\n\t"\
		"addpd		%%xmm1		,%%xmm6\n\t"\
		"movaps		%%xmm7		,      (%%edx)\n\t"\
		"movaps		%%xmm6		, 0x010(%%ebx)\n\t"\
		"/*...Block 8: t0E,t1E,t2E,t3E	*/\n\t"\
		"addl		$0x080		,%%eax\n\t"\
		"addl		$0x080		,%%ebx\n\t"\
		"addl		$0x080		,%%ecx\n\t"\
		"addl		$0x080		,%%edx\n\t"\
		"movaps		      (%%ecx)	,%%xmm4\n\t"\
		"movaps		      (%%edx)	,%%xmm6\n\t"\
		"movaps		 0x010(%%ecx)	,%%xmm5\n\t"\
		"movaps		 0x010(%%edx)	,%%xmm7\n\t"\
		"movaps		      (%%ecx)	,%%xmm0\n\t"\
		"movaps		      (%%edx)	,%%xmm2\n\t"\
		"movaps		 0x010(%%ecx)	,%%xmm1\n\t"\
		"movaps		 0x010(%%edx)	,%%xmm3\n\t"\
		"mulpd		 0x010(%%esi)	,%%xmm4\n\t"\
		"mulpd		 0x30 (%%esi)	,%%xmm6\n\t"\
		"mulpd		      (%%esi)	,%%xmm1\n\t"\
		"mulpd		 0x20 (%%esi)	,%%xmm3\n\t"\
		"mulpd		 0x010(%%esi)	,%%xmm5\n\t"\
		"mulpd		 0x30 (%%esi)	,%%xmm7\n\t"\
		"mulpd		      (%%esi)	,%%xmm0\n\t"\
		"mulpd		 0x20 (%%esi)	,%%xmm2\n\t"\
		"subpd		%%xmm1		,%%xmm4\n\t"\
		"subpd		%%xmm3		,%%xmm6\n\t"\
		"addpd		%%xmm0		,%%xmm5\n\t"\
		"addpd		%%xmm2		,%%xmm7\n\t"\
		"subpd		%%xmm6		,%%xmm4\n\t"\
		"subpd		%%xmm7		,%%xmm5\n\t"\
		"addpd		%%xmm6		,%%xmm6\n\t"\
		"addpd		%%xmm7		,%%xmm7\n\t"\
		"addpd		%%xmm4		,%%xmm6\n\t"\
		"addpd		%%xmm5		,%%xmm7\n\t"\
		"movaps		      (%%ebx)	,%%xmm2\n\t"\
		"movaps		 0x010(%%ebx)	,%%xmm0\n\t"\
		"movaps		      (%%ebx)	,%%xmm1\n\t"\
		"movaps		 0x010(%%ebx)	,%%xmm3\n\t"\
		"mulpd		      (%%edi)	,%%xmm2\n\t"\
		"mulpd		 0x010(%%edi)	,%%xmm0\n\t"\
		"mulpd		      (%%edi)	,%%xmm3\n\t"\
		"mulpd		 0x010(%%edi)	,%%xmm1\n\t"\
		"addpd		%%xmm0		,%%xmm2\n\t"\
		"subpd		%%xmm1		,%%xmm3\n\t"\
		"movaps		      (%%eax)	,%%xmm0\n\t"\
		"movaps		 0x010(%%eax)	,%%xmm1\n\t"\
		"subpd		%%xmm2		,%%xmm0\n\t"\
		"subpd		%%xmm3		,%%xmm1\n\t"\
		"addpd		      (%%eax)	,%%xmm2\n\t"\
		"addpd		 0x010(%%eax)	,%%xmm3\n\t"\
		"subpd		%%xmm4		,%%xmm0\n\t"\
		"subpd		%%xmm5		,%%xmm1\n\t"\
		"addpd		%%xmm4		,%%xmm4\n\t"\
		"addpd		%%xmm5		,%%xmm5\n\t"\
		"movaps		%%xmm0		,      (%%ecx)\n\t"\
		"movaps		%%xmm1		, 0x010(%%ecx)\n\t"\
		"addpd		%%xmm0		,%%xmm4\n\t"\
		"addpd		%%xmm1		,%%xmm5\n\t"\
		"movaps		%%xmm4		,      (%%eax)\n\t"\
		"movaps		%%xmm5		, 0x010(%%eax)\n\t"\
		"subpd		%%xmm7		,%%xmm2\n\t"\
		"subpd		%%xmm6		,%%xmm3\n\t"\
		"addpd		%%xmm7		,%%xmm7\n\t"\
		"addpd		%%xmm6		,%%xmm6\n\t"\
		"movaps		%%xmm2		,      (%%ebx)\n\t"\
		"movaps		%%xmm3		, 0x010(%%edx)\n\t"\
		"addpd		%%xmm2		,%%xmm7\n\t"\
		"addpd		%%xmm3		,%%xmm6\n\t"\
		"movaps		%%xmm7		,      (%%edx)\n\t"\
		"movaps		%%xmm6		, 0x010(%%ebx)\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__add1] "m" (Xadd1)\
		 ,[__r00] "m" (Xr00)\
		 ,[__r10] "m" (Xr10)\
		 ,[__r20] "m" (Xr20)\
		 ,[__r30] "m" (Xr30)\
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__cc0] "m" (Xcc0)\
		 ,[__c00] "m" (Xc00)\
		 ,[__c01] "m" (Xc01)\
		 ,[__c02] "m" (Xc02)\
		 ,[__c03] "m" (Xc03)\
		 ,[__c05] "m" (Xc05)\
		 ,[__c07] "m" (Xc07)\
		: "eax","ebx","ecx","edx","edi","esi"		/* Clobbered registers */\
	);\
	}

	#define SSE2_RADIX32_WRAPPER_DIT(Xadd0,Xadd1,Xisrt2,Xr00,Xr08,Xr10,Xr20,Xr28,Xr30,Xc00,Xc01,Xc02,Xc03,Xc04,Xc05,Xc06,Xc07,Xc08,Xc0A,Xc0C,Xc0E,Xc10,Xc12,Xc14,Xc16,Xc18,Xc1A,Xc1C,Xc1E)\
	{\
	__asm__ volatile (\
		"/************************************************************************/\n\t"\
		"/*...And do an inverse DIT radix-32 pass on the squared-data blocks.	*/\n\t"\
		"/************************************************************************/\n\t"\
		"/*...Block 1: */\n\t"\
		"movl		%[__isrt2]		,%%esi\n\t"\
		"movl		%[__r00]		,%%eax\n\t"\
		"movl		%%eax		,%%ebx\n\t"\
		"movl		%%eax		,%%ecx\n\t"\
		"movl		%%eax		,%%edx\n\t"\
		"addl		$0x200		,%%ebx\n\t"\
		"addl		$0x100		,%%ecx\n\t"\
		"addl		$0x300		,%%edx\n\t"\
		"/*****	SSE2_RADIX4_DIT_IN_PLACE()	*****/\n\t"\
		"movaps		      (%%eax)	,%%xmm0\n\t"\
		"movaps		 0x010(%%eax)	,%%xmm1\n\t"\
		"movaps		      (%%eax)	,%%xmm2\n\t"\
		"movaps		 0x010(%%eax)	,%%xmm3\n\t"\
		"addpd		      (%%ebx)	,%%xmm0\n\t"\
		"addpd		 0x010(%%ebx)	,%%xmm1\n\t"\
		"subpd		      (%%ebx)	,%%xmm2\n\t"\
		"subpd		 0x010(%%ebx)	,%%xmm3\n\t"\
		"movaps		      (%%ecx)	,%%xmm4\n\t"\
		"movaps		 0x010(%%ecx)	,%%xmm5\n\t"\
		"movaps		      (%%ecx)	,%%xmm6\n\t"\
		"movaps		 0x010(%%ecx)	,%%xmm7\n\t"\
		"addpd		      (%%edx)	,%%xmm4\n\t"\
		"addpd		 0x010(%%edx)	,%%xmm5\n\t"\
		"subpd		      (%%edx)	,%%xmm6\n\t"\
		"subpd		 0x010(%%edx)	,%%xmm7\n\t"\
		"subpd		%%xmm4		,%%xmm0\n\t"\
		"subpd		%%xmm5		,%%xmm1\n\t"\
		"movaps		%%xmm0		,      (%%ebx)\n\t"\
		"movaps		%%xmm1		, 0x010(%%ebx)\n\t"\
		"addpd		%%xmm4		,%%xmm4\n\t"\
		"addpd		%%xmm5		,%%xmm5\n\t"\
		"addpd		%%xmm0		,%%xmm4\n\t"\
		"addpd		%%xmm1		,%%xmm5\n\t"\
		"movaps		%%xmm4		,      (%%eax)\n\t"\
		"movaps		%%xmm5		, 0x010(%%eax)\n\t"\
		"subpd		%%xmm7		,%%xmm2\n\t"\
		"subpd		%%xmm6		,%%xmm3\n\t"\
		"movaps		%%xmm2		,      (%%edx)\n\t"\
		"movaps		%%xmm3		, 0x010(%%ecx)\n\t"\
		"addpd		%%xmm7		,%%xmm7\n\t"\
		"addpd		%%xmm6		,%%xmm6\n\t"\
		"addpd		%%xmm2		,%%xmm7\n\t"\
		"addpd		%%xmm3		,%%xmm6\n\t"\
		"movaps		%%xmm7		,      (%%ecx)\n\t"\
		"movaps		%%xmm6		, 0x010(%%edx)\n\t"\
		"addl		$0x080		,%%eax\n\t"\
		"addl		$0x080		,%%ebx\n\t"\
		"addl		$0x080		,%%ecx\n\t"\
		"addl		$0x080		,%%edx\n\t"\
		"movaps		      (%%eax)	,%%xmm0\n\t"\
		"movaps		 0x010(%%eax)	,%%xmm1\n\t"\
		"movaps		      (%%eax)	,%%xmm2\n\t"\
		"movaps		 0x010(%%eax)	,%%xmm3\n\t"\
		"addpd		      (%%ebx)	,%%xmm0\n\t"\
		"addpd		 0x010(%%ebx)	,%%xmm1\n\t"\
		"subpd		      (%%ebx)	,%%xmm2\n\t"\
		"subpd		 0x010(%%ebx)	,%%xmm3\n\t"\
		"movaps		      (%%ecx)	,%%xmm4\n\t"\
		"movaps		 0x010(%%ecx)	,%%xmm5\n\t"\
		"movaps		      (%%ecx)	,%%xmm6\n\t"\
		"movaps		 0x010(%%ecx)	,%%xmm7\n\t"\
		"addpd		      (%%edx)	,%%xmm4\n\t"\
		"addpd		 0x010(%%edx)	,%%xmm5\n\t"\
		"subpd		      (%%edx)	,%%xmm6\n\t"\
		"subpd		 0x010(%%edx)	,%%xmm7\n\t"\
		"subpd		%%xmm4		,%%xmm0\n\t"\
		"subpd		%%xmm5		,%%xmm1\n\t"\
		"movaps		%%xmm0		,      (%%ebx)\n\t"\
		"movaps		%%xmm1		, 0x010(%%ebx)\n\t"\
		"addpd		%%xmm4		,%%xmm4\n\t"\
		"addpd		%%xmm5		,%%xmm5\n\t"\
		"addpd		%%xmm0		,%%xmm4\n\t"\
		"addpd		%%xmm1		,%%xmm5\n\t"\
		"movaps		%%xmm4		,      (%%eax)\n\t"\
		"movaps		%%xmm5		, 0x010(%%eax)\n\t"\
		"subpd		%%xmm7		,%%xmm2\n\t"\
		"subpd		%%xmm6		,%%xmm3\n\t"\
		"addpd		%%xmm7		,%%xmm7\n\t"\
		"addpd		%%xmm6		,%%xmm6\n\t"\
		"addpd		%%xmm2		,%%xmm7\n\t"\
		"addpd		%%xmm3		,%%xmm6\n\t"\
		"movaps		%%xmm3		,%%xmm0\n\t"\
		"movaps		%%xmm6		,%%xmm1\n\t"\
		"subpd		%%xmm7		,%%xmm3\n\t"\
		"subpd		%%xmm2		,%%xmm6\n\t"\
		"addpd		%%xmm7		,%%xmm0\n\t"\
		"addpd		%%xmm2		,%%xmm1\n\t"\
		"mulpd		      (%%esi)	,%%xmm3\n\t"\
		"mulpd		      (%%esi)	,%%xmm6\n\t"\
		"mulpd		      (%%esi)	,%%xmm0\n\t"\
		"mulpd		      (%%esi)	,%%xmm1\n\t"\
		"movaps		%%xmm3		, 0x010(%%ecx)\n\t"\
		"movaps		%%xmm6		, 0x010(%%edx)\n\t"\
		"movaps		%%xmm0		,      (%%ecx)\n\t"\
		"movaps		%%xmm1		,      (%%edx)\n\t"\
		"/*****	SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS(r00,r10,r20,r30,r08,r18,r28,r38)	*****/\n\t"\
		"movaps		-0x080(%%eax)	,%%xmm0\n\t"\
		"movaps		-0x080(%%ebx)	,%%xmm4\n\t"\
		"movaps		-0x070(%%eax)	,%%xmm1\n\t"\
		"movaps		-0x070(%%ebx)	,%%xmm5\n\t"\
		"movaps		      (%%eax)	,%%xmm2\n\t"\
		"movaps		 0x010(%%ebx)	,%%xmm7\n\t"\
		"movaps		 0x010(%%eax)	,%%xmm3\n\t"\
		"movaps		      (%%ebx)	,%%xmm6\n\t"\
		"subpd		%%xmm2		,%%xmm0\n\t"\
		"subpd		%%xmm7		,%%xmm4\n\t"\
		"subpd		%%xmm3		,%%xmm1\n\t"\
		"subpd		%%xmm6		,%%xmm5\n\t"\
		"addpd		%%xmm2		,%%xmm2\n\t"\
		"addpd		%%xmm7		,%%xmm7\n\t"\
		"addpd		%%xmm3		,%%xmm3\n\t"\
		"addpd		%%xmm6		,%%xmm6\n\t"\
		"addpd		%%xmm0		,%%xmm2\n\t"\
		"addpd		%%xmm4		,%%xmm7\n\t"\
		"addpd		%%xmm1		,%%xmm3\n\t"\
		"addpd		%%xmm5		,%%xmm6\n\t"\
		"movaps		%%xmm0		,      (%%eax)\n\t"\
		"movaps		%%xmm4		,      (%%ebx)\n\t"\
		"movaps		%%xmm1		, 0x010(%%eax)\n\t"\
		"movaps		%%xmm5		,-0x070(%%ebx)\n\t"\
		"movaps		%%xmm2		,-0x080(%%eax)\n\t"\
		"movaps		%%xmm7		,-0x080(%%ebx)\n\t"\
		"movaps		%%xmm3		,-0x070(%%eax)\n\t"\
		"movaps		%%xmm6		, 0x010(%%ebx)\n\t"\
		"movaps		-0x080(%%ecx)	,%%xmm0\n\t"\
		"movaps		-0x080(%%edx)	,%%xmm4\n\t"\
		"movaps		-0x070(%%ecx)	,%%xmm1\n\t"\
		"movaps		-0x070(%%edx)	,%%xmm5\n\t"\
		"movaps		      (%%ecx)	,%%xmm2\n\t"\
		"movaps		 0x010(%%edx)	,%%xmm7\n\t"\
		"movaps		 0x010(%%ecx)	,%%xmm3\n\t"\
		"movaps		      (%%edx)	,%%xmm6\n\t"\
		"subpd		%%xmm2		,%%xmm0\n\t"\
		"subpd		%%xmm7		,%%xmm4\n\t"\
		"subpd		%%xmm3		,%%xmm1\n\t"\
		"subpd		%%xmm6		,%%xmm5\n\t"\
		"addpd		%%xmm2		,%%xmm2\n\t"\
		"addpd		%%xmm7		,%%xmm7\n\t"\
		"addpd		%%xmm3		,%%xmm3\n\t"\
		"addpd		%%xmm6		,%%xmm6\n\t"\
		"addpd		%%xmm0		,%%xmm2\n\t"\
		"addpd		%%xmm4		,%%xmm7\n\t"\
		"addpd		%%xmm1		,%%xmm3\n\t"\
		"addpd		%%xmm5		,%%xmm6\n\t"\
		"movaps		%%xmm0		,      (%%ecx)\n\t"\
		"movaps		%%xmm4		,      (%%edx)\n\t"\
		"movaps		%%xmm1		, 0x010(%%ecx)\n\t"\
		"movaps		%%xmm5		,-0x070(%%edx)\n\t"\
		"movaps		%%xmm2		,-0x080(%%ecx)\n\t"\
		"movaps		%%xmm7		,-0x080(%%edx)\n\t"\
		"movaps		%%xmm3		,-0x070(%%ecx)\n\t"\
		"movaps		%%xmm6		, 0x010(%%edx)\n\t"\
		"/*...Block 2:	*/\n\t"\
		"subl		$0x040		,%%eax\n\t"\
		"subl		$0x040		,%%ebx\n\t"\
		"subl		$0x040		,%%ecx\n\t"\
		"subl		$0x040		,%%edx\n\t"\
		"movaps		      (%%eax)	,%%xmm0\n\t"\
		"movaps		 0x010(%%eax)	,%%xmm1\n\t"\
		"movaps		      (%%eax)	,%%xmm2\n\t"\
		"movaps		 0x010(%%eax)	,%%xmm3\n\t"\
		"addpd		      (%%ebx)	,%%xmm0\n\t"\
		"addpd		 0x010(%%ebx)	,%%xmm1\n\t"\
		"subpd		      (%%ebx)	,%%xmm2\n\t"\
		"subpd		 0x010(%%ebx)	,%%xmm3\n\t"\
		"movaps		      (%%ecx)	,%%xmm4\n\t"\
		"movaps		 0x010(%%ecx)	,%%xmm5\n\t"\
		"movaps		      (%%ecx)	,%%xmm6\n\t"\
		"movaps		 0x010(%%ecx)	,%%xmm7\n\t"\
		"addpd		      (%%edx)	,%%xmm4\n\t"\
		"addpd		 0x010(%%edx)	,%%xmm5\n\t"\
		"subpd		      (%%edx)	,%%xmm6\n\t"\
		"subpd		 0x010(%%edx)	,%%xmm7\n\t"\
		"subpd		%%xmm4		,%%xmm0\n\t"\
		"subpd		%%xmm5		,%%xmm1\n\t"\
		"movaps		%%xmm0		,      (%%ebx)\n\t"\
		"movaps		%%xmm1		, 0x010(%%ebx)\n\t"\
		"addpd		%%xmm4		,%%xmm4\n\t"\
		"addpd		%%xmm5		,%%xmm5\n\t"\
		"addpd		%%xmm0		,%%xmm4\n\t"\
		"addpd		%%xmm1		,%%xmm5\n\t"\
		"movaps		%%xmm4		,      (%%eax)\n\t"\
		"movaps		%%xmm5		, 0x010(%%eax)\n\t"\
		"subpd		%%xmm7		,%%xmm2\n\t"\
		"subpd		%%xmm6		,%%xmm3\n\t"\
		"movaps		%%xmm2		,      (%%edx)\n\t"\
		"movaps		%%xmm3		, 0x010(%%ecx)\n\t"\
		"addpd		%%xmm7		,%%xmm7\n\t"\
		"addpd		%%xmm6		,%%xmm6\n\t"\
		"addpd		%%xmm2		,%%xmm7\n\t"\
		"addpd		%%xmm3		,%%xmm6\n\t"\
		"movaps		%%xmm7		,      (%%ecx)\n\t"\
		"movaps		%%xmm6		, 0x010(%%edx)\n\t"\
		"addl		$0x080		,%%eax\n\t"\
		"addl		$0x080		,%%ebx\n\t"\
		"addl		$0x080		,%%ecx\n\t"\
		"addl		$0x080		,%%edx\n\t"\
		"movaps		      (%%eax)	,%%xmm0\n\t"\
		"movaps		 0x010(%%eax)	,%%xmm1\n\t"\
		"movaps		      (%%eax)	,%%xmm2\n\t"\
		"movaps		 0x010(%%eax)	,%%xmm3\n\t"\
		"addpd		      (%%ebx)	,%%xmm0\n\t"\
		"addpd		 0x010(%%ebx)	,%%xmm1\n\t"\
		"subpd		      (%%ebx)	,%%xmm2\n\t"\
		"subpd		 0x010(%%ebx)	,%%xmm3\n\t"\
		"movaps		      (%%ecx)	,%%xmm4\n\t"\
		"movaps		 0x010(%%ecx)	,%%xmm5\n\t"\
		"movaps		      (%%ecx)	,%%xmm6\n\t"\
		"movaps		 0x010(%%ecx)	,%%xmm7\n\t"\
		"addpd		      (%%edx)	,%%xmm4\n\t"\
		"addpd		 0x010(%%edx)	,%%xmm5\n\t"\
		"subpd		      (%%edx)	,%%xmm6\n\t"\
		"subpd		 0x010(%%edx)	,%%xmm7\n\t"\
		"subpd		%%xmm4		,%%xmm0\n\t"\
		"subpd		%%xmm5		,%%xmm1\n\t"\
		"movaps		%%xmm0		,      (%%ebx)\n\t"\
		"movaps		%%xmm1		, 0x010(%%ebx)\n\t"\
		"addpd		%%xmm4		,%%xmm4\n\t"\
		"addpd		%%xmm5		,%%xmm5\n\t"\
		"addpd		%%xmm0		,%%xmm4\n\t"\
		"addpd		%%xmm1		,%%xmm5\n\t"\
		"movaps		%%xmm4		,      (%%eax)\n\t"\
		"movaps		%%xmm5		, 0x010(%%eax)\n\t"\
		"subpd		%%xmm7		,%%xmm2\n\t"\
		"subpd		%%xmm6		,%%xmm3\n\t"\
		"addpd		%%xmm7		,%%xmm7\n\t"\
		"addpd		%%xmm6		,%%xmm6\n\t"\
		"addpd		%%xmm2		,%%xmm7\n\t"\
		"addpd		%%xmm3		,%%xmm6\n\t"\
		"movaps		%%xmm3		,%%xmm0\n\t"\
		"movaps		%%xmm6		,%%xmm1\n\t"\
		"subpd		%%xmm7		,%%xmm3\n\t"\
		"subpd		%%xmm2		,%%xmm6\n\t"\
		"addpd		%%xmm7		,%%xmm0\n\t"\
		"addpd		%%xmm2		,%%xmm1\n\t"\
		"mulpd		      (%%esi)	,%%xmm3\n\t"\
		"mulpd		      (%%esi)	,%%xmm6\n\t"\
		"mulpd		      (%%esi)	,%%xmm0\n\t"\
		"mulpd		      (%%esi)	,%%xmm1\n\t"\
		"movaps		%%xmm3		, 0x010(%%ecx)\n\t"\
		"movaps		%%xmm6		, 0x010(%%edx)\n\t"\
		"movaps		%%xmm0		,      (%%ecx)\n\t"\
		"movaps		%%xmm1		,      (%%edx)\n\t"\
		"/*****	SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS(r04,r14,r24,r34,r0C,r1C,r2C,r3C)	*****/\n\t"\
		"movaps		-0x080(%%eax)	,%%xmm0\n\t"\
		"movaps		-0x080(%%ebx)	,%%xmm4\n\t"\
		"movaps		-0x070(%%eax)	,%%xmm1\n\t"\
		"movaps		-0x070(%%ebx)	,%%xmm5\n\t"\
		"movaps		      (%%eax)	,%%xmm2\n\t"\
		"movaps		 0x010(%%ebx)	,%%xmm7\n\t"\
		"movaps		 0x010(%%eax)	,%%xmm3\n\t"\
		"movaps		      (%%ebx)	,%%xmm6\n\t"\
		"subpd		%%xmm2		,%%xmm0\n\t"\
		"subpd		%%xmm7		,%%xmm4\n\t"\
		"subpd		%%xmm3		,%%xmm1\n\t"\
		"subpd		%%xmm6		,%%xmm5\n\t"\
		"addpd		%%xmm2		,%%xmm2\n\t"\
		"addpd		%%xmm7		,%%xmm7\n\t"\
		"addpd		%%xmm3		,%%xmm3\n\t"\
		"addpd		%%xmm6		,%%xmm6\n\t"\
		"addpd		%%xmm0		,%%xmm2\n\t"\
		"addpd		%%xmm4		,%%xmm7\n\t"\
		"addpd		%%xmm1		,%%xmm3\n\t"\
		"addpd		%%xmm5		,%%xmm6\n\t"\
		"movaps		%%xmm0		,      (%%eax)\n\t"\
		"movaps		%%xmm4		,      (%%ebx)\n\t"\
		"movaps		%%xmm1		, 0x010(%%eax)\n\t"\
		"movaps		%%xmm5		,-0x070(%%ebx)\n\t"\
		"movaps		%%xmm2		,-0x080(%%eax)\n\t"\
		"movaps		%%xmm7		,-0x080(%%ebx)\n\t"\
		"movaps		%%xmm3		,-0x070(%%eax)\n\t"\
		"movaps		%%xmm6		, 0x010(%%ebx)\n\t"\
		"movaps		-0x080(%%ecx)	,%%xmm0\n\t"\
		"movaps		-0x080(%%edx)	,%%xmm4\n\t"\
		"movaps		-0x070(%%ecx)	,%%xmm1\n\t"\
		"movaps		-0x070(%%edx)	,%%xmm5\n\t"\
		"movaps		      (%%ecx)	,%%xmm2\n\t"\
		"movaps		 0x010(%%edx)	,%%xmm7\n\t"\
		"movaps		 0x010(%%ecx)	,%%xmm3\n\t"\
		"movaps		      (%%edx)	,%%xmm6\n\t"\
		"subpd		%%xmm2		,%%xmm0\n\t"\
		"subpd		%%xmm7		,%%xmm4\n\t"\
		"subpd		%%xmm3		,%%xmm1\n\t"\
		"subpd		%%xmm6		,%%xmm5\n\t"\
		"addpd		%%xmm2		,%%xmm2\n\t"\
		"addpd		%%xmm7		,%%xmm7\n\t"\
		"addpd		%%xmm3		,%%xmm3\n\t"\
		"addpd		%%xmm6		,%%xmm6\n\t"\
		"addpd		%%xmm0		,%%xmm2\n\t"\
		"addpd		%%xmm4		,%%xmm7\n\t"\
		"addpd		%%xmm1		,%%xmm3\n\t"\
		"addpd		%%xmm5		,%%xmm6\n\t"\
		"movaps		%%xmm0		,      (%%ecx)\n\t"\
		"movaps		%%xmm4		,      (%%edx)\n\t"\
		"movaps		%%xmm1		, 0x010(%%ecx)\n\t"\
		"movaps		%%xmm5		,-0x070(%%edx)\n\t"\
		"movaps		%%xmm2		,-0x080(%%ecx)\n\t"\
		"movaps		%%xmm7		,-0x080(%%edx)\n\t"\
		"movaps		%%xmm3		,-0x070(%%ecx)\n\t"\
		"movaps		%%xmm6		, 0x010(%%edx)\n\t"\
		"/*...Block 3:	*/\n\t"\
		"subl		$0x0a0		,%%eax\n\t"\
		"subl		$0x0a0		,%%ebx\n\t"\
		"subl		$0x0a0		,%%ecx\n\t"\
		"subl		$0x0a0		,%%edx\n\t"\
		"movaps		      (%%eax)	,%%xmm0\n\t"\
		"movaps		 0x010(%%eax)	,%%xmm1\n\t"\
		"movaps		      (%%eax)	,%%xmm2\n\t"\
		"movaps		 0x010(%%eax)	,%%xmm3\n\t"\
		"addpd		      (%%ebx)	,%%xmm0\n\t"\
		"addpd		 0x010(%%ebx)	,%%xmm1\n\t"\
		"subpd		      (%%ebx)	,%%xmm2\n\t"\
		"subpd		 0x010(%%ebx)	,%%xmm3\n\t"\
		"movaps		      (%%ecx)	,%%xmm4\n\t"\
		"movaps		 0x010(%%ecx)	,%%xmm5\n\t"\
		"movaps		      (%%ecx)	,%%xmm6\n\t"\
		"movaps		 0x010(%%ecx)	,%%xmm7\n\t"\
		"addpd		      (%%edx)	,%%xmm4\n\t"\
		"addpd		 0x010(%%edx)	,%%xmm5\n\t"\
		"subpd		      (%%edx)	,%%xmm6\n\t"\
		"subpd		 0x010(%%edx)	,%%xmm7\n\t"\
		"subpd		%%xmm4		,%%xmm0\n\t"\
		"subpd		%%xmm5		,%%xmm1\n\t"\
		"movaps		%%xmm0		,      (%%ebx)\n\t"\
		"movaps		%%xmm1		, 0x010(%%ebx)\n\t"\
		"addpd		%%xmm4		,%%xmm4\n\t"\
		"addpd		%%xmm5		,%%xmm5\n\t"\
		"addpd		%%xmm0		,%%xmm4\n\t"\
		"addpd		%%xmm1		,%%xmm5\n\t"\
		"movaps		%%xmm4		,      (%%eax)\n\t"\
		"movaps		%%xmm5		, 0x010(%%eax)\n\t"\
		"subpd		%%xmm7		,%%xmm2\n\t"\
		"subpd		%%xmm6		,%%xmm3\n\t"\
		"movaps		%%xmm2		,      (%%edx)\n\t"\
		"movaps		%%xmm3		, 0x010(%%ecx)\n\t"\
		"addpd		%%xmm7		,%%xmm7\n\t"\
		"addpd		%%xmm6		,%%xmm6\n\t"\
		"addpd		%%xmm2		,%%xmm7\n\t"\
		"addpd		%%xmm3		,%%xmm6\n\t"\
		"movaps		%%xmm7		,      (%%ecx)\n\t"\
		"movaps		%%xmm6		, 0x010(%%edx)\n\t"\
		"addl		$0x080		,%%eax\n\t"\
		"addl		$0x080		,%%ebx\n\t"\
		"addl		$0x080		,%%ecx\n\t"\
		"addl		$0x080		,%%edx\n\t"\
		"movaps		      (%%eax)	,%%xmm0\n\t"\
		"movaps		 0x010(%%eax)	,%%xmm1\n\t"\
		"movaps		      (%%eax)	,%%xmm2\n\t"\
		"movaps		 0x010(%%eax)	,%%xmm3\n\t"\
		"addpd		      (%%ebx)	,%%xmm0\n\t"\
		"addpd		 0x010(%%ebx)	,%%xmm1\n\t"\
		"subpd		      (%%ebx)	,%%xmm2\n\t"\
		"subpd		 0x010(%%ebx)	,%%xmm3\n\t"\
		"movaps		      (%%ecx)	,%%xmm4\n\t"\
		"movaps		 0x010(%%ecx)	,%%xmm5\n\t"\
		"movaps		      (%%ecx)	,%%xmm6\n\t"\
		"movaps		 0x010(%%ecx)	,%%xmm7\n\t"\
		"addpd		      (%%edx)	,%%xmm4\n\t"\
		"addpd		 0x010(%%edx)	,%%xmm5\n\t"\
		"subpd		      (%%edx)	,%%xmm6\n\t"\
		"subpd		 0x010(%%edx)	,%%xmm7\n\t"\
		"subpd		%%xmm4		,%%xmm0\n\t"\
		"subpd		%%xmm5		,%%xmm1\n\t"\
		"movaps		%%xmm0		,      (%%ebx)\n\t"\
		"movaps		%%xmm1		, 0x010(%%ebx)\n\t"\
		"addpd		%%xmm4		,%%xmm4\n\t"\
		"addpd		%%xmm5		,%%xmm5\n\t"\
		"addpd		%%xmm0		,%%xmm4\n\t"\
		"addpd		%%xmm1		,%%xmm5\n\t"\
		"movaps		%%xmm4		,      (%%eax)\n\t"\
		"movaps		%%xmm5		, 0x010(%%eax)\n\t"\
		"subpd		%%xmm7		,%%xmm2\n\t"\
		"subpd		%%xmm6		,%%xmm3\n\t"\
		"addpd		%%xmm7		,%%xmm7\n\t"\
		"addpd		%%xmm6		,%%xmm6\n\t"\
		"addpd		%%xmm2		,%%xmm7\n\t"\
		"addpd		%%xmm3		,%%xmm6\n\t"\
		"movaps		%%xmm3		,%%xmm0\n\t"\
		"movaps		%%xmm6		,%%xmm1\n\t"\
		"subpd		%%xmm7		,%%xmm3\n\t"\
		"subpd		%%xmm2		,%%xmm6\n\t"\
		"addpd		%%xmm7		,%%xmm0\n\t"\
		"addpd		%%xmm2		,%%xmm1\n\t"\
		"mulpd		      (%%esi)	,%%xmm3\n\t"\
		"mulpd		      (%%esi)	,%%xmm6\n\t"\
		"mulpd		      (%%esi)	,%%xmm0\n\t"\
		"mulpd		      (%%esi)	,%%xmm1\n\t"\
		"movaps		%%xmm3		, 0x010(%%ecx)\n\t"\
		"movaps		%%xmm6		, 0x010(%%edx)\n\t"\
		"movaps		%%xmm0		,      (%%ecx)\n\t"\
		"movaps		%%xmm1		,      (%%edx)\n\t"\
		"/*****	SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS(r02,r12,r22,r32,r0A,r1A,r2A,r3A)	*****/\n\t"\
		"movaps		-0x080(%%eax)	,%%xmm0\n\t"\
		"movaps		-0x080(%%ebx)	,%%xmm4\n\t"\
		"movaps		-0x070(%%eax)	,%%xmm1\n\t"\
		"movaps		-0x070(%%ebx)	,%%xmm5\n\t"\
		"movaps		      (%%eax)	,%%xmm2\n\t"\
		"movaps		 0x010(%%ebx)	,%%xmm7\n\t"\
		"movaps		 0x010(%%eax)	,%%xmm3\n\t"\
		"movaps		      (%%ebx)	,%%xmm6\n\t"\
		"subpd		%%xmm2		,%%xmm0\n\t"\
		"subpd		%%xmm7		,%%xmm4\n\t"\
		"subpd		%%xmm3		,%%xmm1\n\t"\
		"subpd		%%xmm6		,%%xmm5\n\t"\
		"addpd		%%xmm2		,%%xmm2\n\t"\
		"addpd		%%xmm7		,%%xmm7\n\t"\
		"addpd		%%xmm3		,%%xmm3\n\t"\
		"addpd		%%xmm6		,%%xmm6\n\t"\
		"addpd		%%xmm0		,%%xmm2\n\t"\
		"addpd		%%xmm4		,%%xmm7\n\t"\
		"addpd		%%xmm1		,%%xmm3\n\t"\
		"addpd		%%xmm5		,%%xmm6\n\t"\
		"movaps		%%xmm0		,      (%%eax)\n\t"\
		"movaps		%%xmm4		,      (%%ebx)\n\t"\
		"movaps		%%xmm1		, 0x010(%%eax)\n\t"\
		"movaps		%%xmm5		,-0x070(%%ebx)\n\t"\
		"movaps		%%xmm2		,-0x080(%%eax)\n\t"\
		"movaps		%%xmm7		,-0x080(%%ebx)\n\t"\
		"movaps		%%xmm3		,-0x070(%%eax)\n\t"\
		"movaps		%%xmm6		, 0x010(%%ebx)\n\t"\
		"movaps		-0x080(%%ecx)	,%%xmm0\n\t"\
		"movaps		-0x080(%%edx)	,%%xmm4\n\t"\
		"movaps		-0x070(%%ecx)	,%%xmm1\n\t"\
		"movaps		-0x070(%%edx)	,%%xmm5\n\t"\
		"movaps		      (%%ecx)	,%%xmm2\n\t"\
		"movaps		 0x010(%%edx)	,%%xmm7\n\t"\
		"movaps		 0x010(%%ecx)	,%%xmm3\n\t"\
		"movaps		      (%%edx)	,%%xmm6\n\t"\
		"subpd		%%xmm2		,%%xmm0\n\t"\
		"subpd		%%xmm7		,%%xmm4\n\t"\
		"subpd		%%xmm3		,%%xmm1\n\t"\
		"subpd		%%xmm6		,%%xmm5\n\t"\
		"addpd		%%xmm2		,%%xmm2\n\t"\
		"addpd		%%xmm7		,%%xmm7\n\t"\
		"addpd		%%xmm3		,%%xmm3\n\t"\
		"addpd		%%xmm6		,%%xmm6\n\t"\
		"addpd		%%xmm0		,%%xmm2\n\t"\
		"addpd		%%xmm4		,%%xmm7\n\t"\
		"addpd		%%xmm1		,%%xmm3\n\t"\
		"addpd		%%xmm5		,%%xmm6\n\t"\
		"movaps		%%xmm0		,      (%%ecx)\n\t"\
		"movaps		%%xmm4		,      (%%edx)\n\t"\
		"movaps		%%xmm1		, 0x010(%%ecx)\n\t"\
		"movaps		%%xmm5		,-0x070(%%edx)\n\t"\
		"movaps		%%xmm2		,-0x080(%%ecx)\n\t"\
		"movaps		%%xmm7		,-0x080(%%edx)\n\t"\
		"movaps		%%xmm3		,-0x070(%%ecx)\n\t"\
		"movaps		%%xmm6		, 0x010(%%edx)\n\t"\
		"/*...Block 4:	*/\n\t"\
		"subl		$0x040		,%%eax\n\t"\
		"subl		$0x040		,%%ebx\n\t"\
		"subl		$0x040		,%%ecx\n\t"\
		"subl		$0x040		,%%edx\n\t"\
		"movaps		      (%%eax)	,%%xmm0\n\t"\
		"movaps		 0x010(%%eax)	,%%xmm1\n\t"\
		"movaps		      (%%eax)	,%%xmm2\n\t"\
		"movaps		 0x010(%%eax)	,%%xmm3\n\t"\
		"addpd		      (%%ebx)	,%%xmm0\n\t"\
		"addpd		 0x010(%%ebx)	,%%xmm1\n\t"\
		"subpd		      (%%ebx)	,%%xmm2\n\t"\
		"subpd		 0x010(%%ebx)	,%%xmm3\n\t"\
		"movaps		      (%%ecx)	,%%xmm4\n\t"\
		"movaps		 0x010(%%ecx)	,%%xmm5\n\t"\
		"movaps		      (%%ecx)	,%%xmm6\n\t"\
		"movaps		 0x010(%%ecx)	,%%xmm7\n\t"\
		"addpd		      (%%edx)	,%%xmm4\n\t"\
		"addpd		 0x010(%%edx)	,%%xmm5\n\t"\
		"subpd		      (%%edx)	,%%xmm6\n\t"\
		"subpd		 0x010(%%edx)	,%%xmm7\n\t"\
		"subpd		%%xmm4		,%%xmm0\n\t"\
		"subpd		%%xmm5		,%%xmm1\n\t"\
		"movaps		%%xmm0		,      (%%ebx)\n\t"\
		"movaps		%%xmm1		, 0x010(%%ebx)\n\t"\
		"addpd		%%xmm4		,%%xmm4\n\t"\
		"addpd		%%xmm5		,%%xmm5\n\t"\
		"addpd		%%xmm0		,%%xmm4\n\t"\
		"addpd		%%xmm1		,%%xmm5\n\t"\
		"movaps		%%xmm4		,      (%%eax)\n\t"\
		"movaps		%%xmm5		, 0x010(%%eax)\n\t"\
		"subpd		%%xmm7		,%%xmm2\n\t"\
		"subpd		%%xmm6		,%%xmm3\n\t"\
		"movaps		%%xmm2		,      (%%edx)\n\t"\
		"movaps		%%xmm3		, 0x010(%%ecx)\n\t"\
		"addpd		%%xmm7		,%%xmm7\n\t"\
		"addpd		%%xmm6		,%%xmm6\n\t"\
		"addpd		%%xmm2		,%%xmm7\n\t"\
		"addpd		%%xmm3		,%%xmm6\n\t"\
		"movaps		%%xmm7		,      (%%ecx)\n\t"\
		"movaps		%%xmm6		, 0x010(%%edx)\n\t"\
		"addl		$0x080		,%%eax\n\t"\
		"addl		$0x080		,%%ebx\n\t"\
		"addl		$0x080		,%%ecx\n\t"\
		"addl		$0x080		,%%edx\n\t"\
		"movaps		      (%%eax)	,%%xmm0\n\t"\
		"movaps		 0x010(%%eax)	,%%xmm1\n\t"\
		"movaps		      (%%eax)	,%%xmm2\n\t"\
		"movaps		 0x010(%%eax)	,%%xmm3\n\t"\
		"addpd		      (%%ebx)	,%%xmm0\n\t"\
		"addpd		 0x010(%%ebx)	,%%xmm1\n\t"\
		"subpd		      (%%ebx)	,%%xmm2\n\t"\
		"subpd		 0x010(%%ebx)	,%%xmm3\n\t"\
		"movaps		      (%%ecx)	,%%xmm4\n\t"\
		"movaps		 0x010(%%ecx)	,%%xmm5\n\t"\
		"movaps		      (%%ecx)	,%%xmm6\n\t"\
		"movaps		 0x010(%%ecx)	,%%xmm7\n\t"\
		"addpd		      (%%edx)	,%%xmm4\n\t"\
		"addpd		 0x010(%%edx)	,%%xmm5\n\t"\
		"subpd		      (%%edx)	,%%xmm6\n\t"\
		"subpd		 0x010(%%edx)	,%%xmm7\n\t"\
		"subpd		%%xmm4		,%%xmm0\n\t"\
		"subpd		%%xmm5		,%%xmm1\n\t"\
		"movaps		%%xmm0		,      (%%ebx)\n\t"\
		"movaps		%%xmm1		, 0x010(%%ebx)\n\t"\
		"addpd		%%xmm4		,%%xmm4\n\t"\
		"addpd		%%xmm5		,%%xmm5\n\t"\
		"addpd		%%xmm0		,%%xmm4\n\t"\
		"addpd		%%xmm1		,%%xmm5\n\t"\
		"movaps		%%xmm4		,      (%%eax)\n\t"\
		"movaps		%%xmm5		, 0x010(%%eax)\n\t"\
		"subpd		%%xmm7		,%%xmm2\n\t"\
		"subpd		%%xmm6		,%%xmm3\n\t"\
		"addpd		%%xmm7		,%%xmm7\n\t"\
		"addpd		%%xmm6		,%%xmm6\n\t"\
		"addpd		%%xmm2		,%%xmm7\n\t"\
		"addpd		%%xmm3		,%%xmm6\n\t"\
		"movaps		%%xmm3		,%%xmm0\n\t"\
		"movaps		%%xmm6		,%%xmm1\n\t"\
		"subpd		%%xmm7		,%%xmm3\n\t"\
		"subpd		%%xmm2		,%%xmm6\n\t"\
		"addpd		%%xmm7		,%%xmm0\n\t"\
		"addpd		%%xmm2		,%%xmm1\n\t"\
		"mulpd		      (%%esi)	,%%xmm3\n\t"\
		"mulpd		      (%%esi)	,%%xmm6\n\t"\
		"mulpd		      (%%esi)	,%%xmm0\n\t"\
		"mulpd		      (%%esi)	,%%xmm1\n\t"\
		"movaps		%%xmm3		, 0x010(%%ecx)\n\t"\
		"movaps		%%xmm6		, 0x010(%%edx)\n\t"\
		"movaps		%%xmm0		,      (%%ecx)\n\t"\
		"movaps		%%xmm1		,      (%%edx)\n\t"\
		"/*****	SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS(r06,r16,r26,r36,r0E,r1E,r2E,r3E)	*****/\n\t"\
		"movaps		-0x080(%%eax)	,%%xmm0\n\t"\
		"movaps		-0x080(%%ebx)	,%%xmm4\n\t"\
		"movaps		-0x070(%%eax)	,%%xmm1\n\t"\
		"movaps		-0x070(%%ebx)	,%%xmm5\n\t"\
		"movaps		      (%%eax)	,%%xmm2\n\t"\
		"movaps		 0x010(%%ebx)	,%%xmm7\n\t"\
		"movaps		 0x010(%%eax)	,%%xmm3\n\t"\
		"movaps		      (%%ebx)	,%%xmm6\n\t"\
		"subpd		%%xmm2		,%%xmm0\n\t"\
		"subpd		%%xmm7		,%%xmm4\n\t"\
		"subpd		%%xmm3		,%%xmm1\n\t"\
		"subpd		%%xmm6		,%%xmm5\n\t"\
		"addpd		%%xmm2		,%%xmm2\n\t"\
		"addpd		%%xmm7		,%%xmm7\n\t"\
		"addpd		%%xmm3		,%%xmm3\n\t"\
		"addpd		%%xmm6		,%%xmm6\n\t"\
		"addpd		%%xmm0		,%%xmm2\n\t"\
		"addpd		%%xmm4		,%%xmm7\n\t"\
		"addpd		%%xmm1		,%%xmm3\n\t"\
		"addpd		%%xmm5		,%%xmm6\n\t"\
		"movaps		%%xmm0		,      (%%eax)\n\t"\
		"movaps		%%xmm4		,      (%%ebx)\n\t"\
		"movaps		%%xmm1		, 0x010(%%eax)\n\t"\
		"movaps		%%xmm5		,-0x070(%%ebx)\n\t"\
		"movaps		%%xmm2		,-0x080(%%eax)\n\t"\
		"movaps		%%xmm7		,-0x080(%%ebx)\n\t"\
		"movaps		%%xmm3		,-0x070(%%eax)\n\t"\
		"movaps		%%xmm6		, 0x010(%%ebx)\n\t"\
		"movaps		-0x080(%%ecx)	,%%xmm0\n\t"\
		"movaps		-0x080(%%edx)	,%%xmm4\n\t"\
		"movaps		-0x070(%%ecx)	,%%xmm1\n\t"\
		"movaps		-0x070(%%edx)	,%%xmm5\n\t"\
		"movaps		      (%%ecx)	,%%xmm2\n\t"\
		"movaps		 0x010(%%edx)	,%%xmm7\n\t"\
		"movaps		 0x010(%%ecx)	,%%xmm3\n\t"\
		"movaps		      (%%edx)	,%%xmm6\n\t"\
		"subpd		%%xmm2		,%%xmm0\n\t"\
		"subpd		%%xmm7		,%%xmm4\n\t"\
		"subpd		%%xmm3		,%%xmm1\n\t"\
		"subpd		%%xmm6		,%%xmm5\n\t"\
		"addpd		%%xmm2		,%%xmm2\n\t"\
		"addpd		%%xmm7		,%%xmm7\n\t"\
		"addpd		%%xmm3		,%%xmm3\n\t"\
		"addpd		%%xmm6		,%%xmm6\n\t"\
		"addpd		%%xmm0		,%%xmm2\n\t"\
		"addpd		%%xmm4		,%%xmm7\n\t"\
		"addpd		%%xmm1		,%%xmm3\n\t"\
		"addpd		%%xmm5		,%%xmm6\n\t"\
		"movaps		%%xmm0		,      (%%ecx)\n\t"\
		"movaps		%%xmm4		,      (%%edx)\n\t"\
		"movaps		%%xmm1		, 0x010(%%ecx)\n\t"\
		"movaps		%%xmm5		,-0x070(%%edx)\n\t"\
		"movaps		%%xmm2		,-0x080(%%ecx)\n\t"\
		"movaps		%%xmm7		,-0x080(%%edx)\n\t"\
		"movaps		%%xmm3		,-0x070(%%ecx)\n\t"\
		"movaps		%%xmm6		, 0x010(%%edx)\n\t"\
		"/***************************************************************************************/\n\t"\
		"/* Now do eight more radix-4 transforms, including the internal and external twiddles: */\n\t"\
		"/***************************************************************************************/\n\t"\
		"/************************************************/\n\t"\
		"/*...Block 2: t02,t12,t22,t32 -> r10,14,12,16:	*/\n\t"\
		"/************************************************/\n\t"\
		"movl		%[__isrt2]		,%%esi\n\t"\
		"movl		%[__r10]		,%%eax\n\t"\
		"movl		%[__add1]		,%%ebx\n\t"\
		"movl		%%esi			,%%ecx\n\t"\
		"movl		%%esi			,%%edx\n\t"\
		"movl		%[__c01]		,%%edi\n\t"\
		"addl		$0x010			,%%ecx	/* cc0 */\n\t"\
		"addl		$0x030			,%%edx	/* cc1 */\n\t"\
		"movaps		 0x020(%%eax)	,%%xmm4\n\t"\
		"movaps		 0x060(%%eax)	,%%xmm0\n\t"\
		"movaps		 0x030(%%eax)	,%%xmm5\n\t"\
		"movaps		 0x070(%%eax)	,%%xmm1\n\t"\
		"movaps		 0x020(%%eax)	,%%xmm6\n\t"\
		"movaps		 0x060(%%eax)	,%%xmm2\n\t"\
		"movaps		 0x030(%%eax)	,%%xmm7\n\t"\
		"movaps		 0x070(%%eax)	,%%xmm3\n\t"\
		"mulpd		      (%%edx)	,%%xmm4\n\t"\
		"mulpd		 0x20 (%%edx)	,%%xmm0\n\t"\
		"mulpd		      (%%edx)	,%%xmm5\n\t"\
		"mulpd		 0x20 (%%edx)	,%%xmm1\n\t"\
		"mulpd		 0x010(%%edx)	,%%xmm6\n\t"\
		"mulpd		 0x30 (%%edx)	,%%xmm2\n\t"\
		"mulpd		 0x010(%%edx)	,%%xmm7\n\t"\
		"mulpd		 0x30 (%%edx)	,%%xmm3\n\t"\
		"subpd		%%xmm6		,%%xmm5\n\t"\
		"subpd		%%xmm2		,%%xmm1\n\t"\
		"addpd		%%xmm7		,%%xmm4\n\t"\
		"addpd		%%xmm3		,%%xmm0\n\t"\
		"movaps		%%xmm5		,%%xmm7\n\t"\
		"movaps		%%xmm4		,%%xmm6\n\t"\
		"subpd		%%xmm0		,%%xmm4\n\t"\
		"subpd		%%xmm1		,%%xmm5\n\t"\
		"addpd		%%xmm0		,%%xmm6\n\t"\
		"addpd		%%xmm1		,%%xmm7\n\t"\
		"movaps		 0x040(%%eax)	,%%xmm2\n\t"\
		"movaps		 0x050(%%eax)	,%%xmm3\n\t"\
		"movaps		 0x040(%%eax)	,%%xmm0\n\t"\
		"movaps		 0x050(%%eax)	,%%xmm1\n\t"\
		"mulpd		      (%%ecx)	,%%xmm2\n\t"\
		"mulpd		 0x010(%%ecx)	,%%xmm1\n\t"\
		"mulpd		      (%%ecx)	,%%xmm3\n\t"\
		"mulpd		 0x010(%%ecx)	,%%xmm0\n\t"\
		"addpd		%%xmm1		,%%xmm2\n\t"\
		"subpd		%%xmm0		,%%xmm3\n\t"\
		"movaps		      (%%eax)	,%%xmm0\n\t"\
		"movaps		 0x010(%%eax)	,%%xmm1\n\t"\
		"subpd		%%xmm2		,%%xmm0\n\t"\
		"subpd		%%xmm3		,%%xmm1\n\t"\
		"addpd		%%xmm2		,%%xmm2\n\t"\
		"addpd		%%xmm3		,%%xmm3\n\t"\
		"addpd		%%xmm0		,%%xmm2\n\t"\
		"addpd		%%xmm1		,%%xmm3\n\t"\
		"subpd		%%xmm6		,%%xmm2\n\t"\
		"subpd		%%xmm7		,%%xmm3\n\t"\
		"addpd		%%xmm6		,%%xmm6\n\t"\
		"addpd		%%xmm7		,%%xmm7\n\t"\
		"addpd		%%xmm2		,%%xmm6\n\t"\
		"addpd		%%xmm3		,%%xmm7\n\t"\
		"movaps		%%xmm2		, 0x020(%%eax)\n\t"\
		"movaps		%%xmm3		, 0x030(%%eax)\n\t"\
		"movaps		%%xmm6		,%%xmm2\n\t"\
		"movaps		%%xmm7		,%%xmm3\n\t"\
		"mulpd		      (%%edi)	,%%xmm6\n\t"\
		"mulpd		      (%%edi)	,%%xmm7\n\t"\
		"mulpd		 0x010(%%edi)	,%%xmm2\n\t"\
		"mulpd		 0x010(%%edi)	,%%xmm3\n\t"\
		"subpd		%%xmm2		,%%xmm7\n\t"\
		"addpd		%%xmm3		,%%xmm6\n\t"\
		"movaps		%%xmm7		, 0x010(%%ebx)\n\t"\
		"movaps		%%xmm6		,      (%%ebx)\n\t"\
		"movaps		 0x020(%%eax)	,%%xmm6\n\t"\
		"movaps		 0x030(%%eax)	,%%xmm7\n\t"\
		"movaps		%%xmm6		,%%xmm2\n\t"\
		"movaps		%%xmm7		,%%xmm3\n\t"\
		"mulpd		 0x20 (%%edi)	,%%xmm6\n\t"\
		"mulpd		 0x20 (%%edi)	,%%xmm7\n\t"\
		"mulpd		 0x30 (%%edi)	,%%xmm2\n\t"\
		"mulpd		 0x30 (%%edi)	,%%xmm3\n\t"\
		"subpd		%%xmm2		,%%xmm7\n\t"\
		"addpd		%%xmm3		,%%xmm6\n\t"\
		"movaps		%%xmm7		, 0x110(%%ebx)\n\t"\
		"movaps		%%xmm6		, 0x100(%%ebx)\n\t"\
		"addl		$0x040		,%%edi\n\t"\
		"subpd		%%xmm5		,%%xmm0\n\t"\
		"subpd		%%xmm4		,%%xmm1\n\t"\
		"addpd		%%xmm5		,%%xmm5\n\t"\
		"addpd		%%xmm4		,%%xmm4\n\t"\
		"addpd		%%xmm0		,%%xmm5\n\t"\
		"addpd		%%xmm1		,%%xmm4\n\t"\
		"movaps		%%xmm5		,%%xmm2\n\t"\
		"movaps		%%xmm1		,%%xmm3\n\t"\
		"mulpd		      (%%edi)	,%%xmm5\n\t"\
		"mulpd		      (%%edi)	,%%xmm1\n\t"\
		"mulpd		 0x010(%%edi)	,%%xmm2\n\t"\
		"mulpd		 0x010(%%edi)	,%%xmm3\n\t"\
		"subpd		%%xmm2		,%%xmm1\n\t"\
		"addpd		%%xmm3		,%%xmm5\n\t"\
		"movaps		%%xmm1		, 0x090(%%ebx)\n\t"\
		"movaps		%%xmm5		, 0x080(%%ebx)\n\t"\
		"movaps		%%xmm0		,%%xmm2\n\t"\
		"movaps		%%xmm4		,%%xmm3\n\t"\
		"mulpd		 0x20 (%%edi)	,%%xmm0\n\t"\
		"mulpd		 0x20 (%%edi)	,%%xmm4\n\t"\
		"mulpd		 0x30 (%%edi)	,%%xmm2\n\t"\
		"mulpd		 0x30 (%%edi)	,%%xmm3\n\t"\
		"subpd		%%xmm2		,%%xmm4\n\t"\
		"addpd		%%xmm3		,%%xmm0\n\t"\
		"movaps		%%xmm4		, 0x190(%%ebx)\n\t"\
		"movaps		%%xmm0		, 0x180(%%ebx)\n\t"\
		"/************************************************/\n\t"\
		"/*...Block 6: t0A,t1A,t2A,t3A -> r18,1C,1A,1E:	*/\n\t"\
		"/************************************************/\n\t"\
		"addl		$0x080		,%%eax\n\t"\
		"movaps		 0x020(%%eax)	,%%xmm4\n\t"\
		"movaps		 0x060(%%eax)	,%%xmm0\n\t"\
		"movaps		 0x030(%%eax)	,%%xmm5\n\t"\
		"movaps		 0x070(%%eax)	,%%xmm1\n\t"\
		"movaps		 0x020(%%eax)	,%%xmm6\n\t"\
		"movaps		 0x060(%%eax)	,%%xmm2\n\t"\
		"movaps		 0x030(%%eax)	,%%xmm7\n\t"\
		"movaps		 0x070(%%eax)	,%%xmm3\n\t"\
		"mulpd		 0x30 (%%edx)	,%%xmm4\n\t"\
		"mulpd		      (%%edx)	,%%xmm0\n\t"\
		"mulpd		 0x30 (%%edx)	,%%xmm5\n\t"\
		"mulpd		      (%%edx)	,%%xmm1\n\t"\
		"mulpd		 0x20 (%%edx)	,%%xmm6\n\t"\
		"mulpd		 0x010(%%edx)	,%%xmm2\n\t"\
		"mulpd		 0x20 (%%edx)	,%%xmm7\n\t"\
		"mulpd		 0x010(%%edx)	,%%xmm3\n\t"\
		"subpd		%%xmm6		,%%xmm5\n\t"\
		"addpd		%%xmm2		,%%xmm1\n\t"\
		"addpd		%%xmm7		,%%xmm4\n\t"\
		"subpd		%%xmm3		,%%xmm0\n\t"\
		"movaps		%%xmm5		,%%xmm7\n\t"\
		"movaps		%%xmm4		,%%xmm6\n\t"\
		"addpd		%%xmm0		,%%xmm4\n\t"\
		"addpd		%%xmm1		,%%xmm5\n\t"\
		"subpd		%%xmm0		,%%xmm6\n\t"\
		"subpd		%%xmm1		,%%xmm7\n\t"\
		"movaps		 0x040(%%eax)	,%%xmm2\n\t"\
		"movaps		 0x050(%%eax)	,%%xmm3\n\t"\
		"movaps		 0x040(%%eax)	,%%xmm0\n\t"\
		"movaps		 0x050(%%eax)	,%%xmm1\n\t"\
		"mulpd		 0x010(%%ecx)	,%%xmm2\n\t"\
		"mulpd		      (%%ecx)	,%%xmm1\n\t"\
		"mulpd		 0x010(%%ecx)	,%%xmm3\n\t"\
		"mulpd		      (%%ecx)	,%%xmm0\n\t"\
		"subpd		%%xmm1		,%%xmm2\n\t"\
		"addpd		%%xmm0		,%%xmm3\n\t"\
		"movaps		      (%%eax)	,%%xmm0\n\t"\
		"movaps		 0x010(%%eax)	,%%xmm1\n\t"\
		"subpd		%%xmm2		,%%xmm0\n\t"\
		"subpd		%%xmm3		,%%xmm1\n\t"\
		"addpd		%%xmm2		,%%xmm2\n\t"\
		"addpd		%%xmm3		,%%xmm3\n\t"\
		"addpd		%%xmm0		,%%xmm2\n\t"\
		"addpd		%%xmm1		,%%xmm3\n\t"\
		"addl		$0x040		,%%edi\n\t"\
		"subpd		%%xmm6		,%%xmm0\n\t"\
		"subpd		%%xmm7		,%%xmm1\n\t"\
		"addpd		%%xmm6		,%%xmm6\n\t"\
		"addpd		%%xmm7		,%%xmm7\n\t"\
		"addpd		%%xmm0		,%%xmm6\n\t"\
		"addpd		%%xmm1		,%%xmm7\n\t"\
		"movaps		%%xmm0		, 0x020(%%eax)\n\t"\
		"movaps		%%xmm1		, 0x030(%%eax)\n\t"\
		"movaps		%%xmm6		,%%xmm0\n\t"\
		"movaps		%%xmm7		,%%xmm1\n\t"\
		"mulpd		      (%%edi)	,%%xmm6\n\t"\
		"mulpd		      (%%edi)	,%%xmm7\n\t"\
		"mulpd		 0x010(%%edi)	,%%xmm0\n\t"\
		"mulpd		 0x010(%%edi)	,%%xmm1\n\t"\
		"subpd		%%xmm0		,%%xmm7\n\t"\
		"addpd		%%xmm1		,%%xmm6\n\t"\
		"movaps		%%xmm7		, 0x050(%%ebx)\n\t"\
		"movaps		%%xmm6		, 0x040(%%ebx)\n\t"\
		"movaps		 0x020(%%eax)	,%%xmm6\n\t"\
		"movaps		 0x030(%%eax)	,%%xmm7\n\t"\
		"movaps		%%xmm6		,%%xmm0\n\t"\
		"movaps		%%xmm7		,%%xmm1\n\t"\
		"mulpd		 0x20 (%%edi)	,%%xmm6\n\t"\
		"mulpd		 0x20 (%%edi)	,%%xmm7\n\t"\
		"mulpd		 0x30 (%%edi)	,%%xmm0\n\t"\
		"mulpd		 0x30 (%%edi)	,%%xmm1\n\t"\
		"subpd		%%xmm0		,%%xmm7\n\t"\
		"addpd		%%xmm1		,%%xmm6\n\t"\
		"movaps		%%xmm7		, 0x150(%%ebx)\n\t"\
		"movaps		%%xmm6		, 0x140(%%ebx)\n\t"\
		"addl		$0x040		,%%edi\n\t"\
		"subpd		%%xmm5		,%%xmm2\n\t"\
		"subpd		%%xmm4		,%%xmm3\n\t"\
		"addpd		%%xmm5		,%%xmm5\n\t"\
		"addpd		%%xmm4		,%%xmm4\n\t"\
		"addpd		%%xmm2		,%%xmm5\n\t"\
		"addpd		%%xmm3		,%%xmm4\n\t"\
		"movaps		%%xmm5		,%%xmm0\n\t"\
		"movaps		%%xmm3		,%%xmm1\n\t"\
		"mulpd		      (%%edi)	,%%xmm5\n\t"\
		"mulpd		      (%%edi)	,%%xmm3\n\t"\
		"mulpd		 0x010(%%edi)	,%%xmm0\n\t"\
		"mulpd		 0x010(%%edi)	,%%xmm1\n\t"\
		"subpd		%%xmm0		,%%xmm3\n\t"\
		"addpd		%%xmm1		,%%xmm5\n\t"\
		"movaps		%%xmm3		, 0x0d0(%%ebx)\n\t"\
		"movaps		%%xmm5		, 0x0c0(%%ebx)\n\t"\
		"movaps		%%xmm2		,%%xmm0\n\t"\
		"movaps		%%xmm4		,%%xmm1\n\t"\
		"mulpd		 0x20 (%%edi)	,%%xmm2\n\t"\
		"mulpd		 0x20 (%%edi)	,%%xmm4\n\t"\
		"mulpd		 0x30 (%%edi)	,%%xmm0\n\t"\
		"mulpd		 0x30 (%%edi)	,%%xmm1\n\t"\
		"subpd		%%xmm0		,%%xmm4\n\t"\
		"addpd		%%xmm1		,%%xmm2\n\t"\
		"movaps		%%xmm4		, 0x1d0(%%ebx)\n\t"\
		"movaps		%%xmm2		, 0x1c0(%%ebx)\n\t"\
		"/************************************************/\n\t"\
		"/*...Block 4: t06,t16,t26,t36 -> r30,34,32,36:	*/\n\t"\
		"/************************************************/\n\t"\
		"addl		$0x180		,%%eax\n\t"\
		"movaps		 0x020(%%eax)	,%%xmm4\n\t"\
		"movaps		 0x060(%%eax)	,%%xmm0\n\t"\
		"movaps		 0x030(%%eax)	,%%xmm5\n\t"\
		"movaps		 0x070(%%eax)	,%%xmm1\n\t"\
		"movaps		 0x020(%%eax)	,%%xmm6\n\t"\
		"movaps		 0x060(%%eax)	,%%xmm2\n\t"\
		"movaps		 0x030(%%eax)	,%%xmm7\n\t"\
		"movaps		 0x070(%%eax)	,%%xmm3\n\t"\
		"mulpd		 0x20 (%%edx)	,%%xmm4\n\t"\
		"mulpd		 0x010(%%edx)	,%%xmm0\n\t"\
		"mulpd		 0x20 (%%edx)	,%%xmm5\n\t"\
		"mulpd		 0x010(%%edx)	,%%xmm1\n\t"\
		"mulpd		 0x30 (%%edx)	,%%xmm6\n\t"\
		"mulpd		      (%%edx)	,%%xmm2\n\t"\
		"mulpd		 0x30 (%%edx)	,%%xmm7\n\t"\
		"mulpd		      (%%edx)	,%%xmm3\n\t"\
		"subpd		%%xmm6		,%%xmm5\n\t"\
		"addpd		%%xmm2		,%%xmm1\n\t"\
		"addpd		%%xmm7		,%%xmm4\n\t"\
		"subpd		%%xmm3		,%%xmm0\n\t"\
		"movaps		%%xmm5		,%%xmm7\n\t"\
		"movaps		%%xmm4		,%%xmm6\n\t"\
		"addpd		%%xmm0		,%%xmm4\n\t"\
		"addpd		%%xmm1		,%%xmm5\n\t"\
		"subpd		%%xmm0		,%%xmm6\n\t"\
		"subpd		%%xmm1		,%%xmm7\n\t"\
		"movaps		 0x040(%%eax)	,%%xmm2\n\t"\
		"movaps		 0x050(%%eax)	,%%xmm3\n\t"\
		"movaps		 0x040(%%eax)	,%%xmm0\n\t"\
		"movaps		 0x050(%%eax)	,%%xmm1\n\t"\
		"mulpd		 0x010(%%ecx)	,%%xmm2\n\t"\
		"mulpd		      (%%ecx)	,%%xmm1\n\t"\
		"mulpd		 0x010(%%ecx)	,%%xmm3\n\t"\
		"mulpd		      (%%ecx)	,%%xmm0\n\t"\
		"addpd		%%xmm1		,%%xmm2\n\t"\
		"subpd		%%xmm0		,%%xmm3\n\t"\
		"movaps		      (%%eax)	,%%xmm0\n\t"\
		"movaps		 0x010(%%eax)	,%%xmm1\n\t"\
		"subpd		%%xmm2		,%%xmm0\n\t"\
		"subpd		%%xmm3		,%%xmm1\n\t"\
		"addpd		%%xmm2		,%%xmm2\n\t"\
		"addpd		%%xmm3		,%%xmm3\n\t"\
		"addpd		%%xmm0		,%%xmm2\n\t"\
		"addpd		%%xmm1		,%%xmm3\n\t"\
		"addl		$0x040		,%%edi\n\t"\
		"subpd		%%xmm6		,%%xmm2\n\t"\
		"subpd		%%xmm7		,%%xmm3\n\t"\
		"addpd		%%xmm6		,%%xmm6\n\t"\
		"addpd		%%xmm7		,%%xmm7\n\t"\
		"addpd		%%xmm2		,%%xmm6\n\t"\
		"addpd		%%xmm3		,%%xmm7\n\t"\
		"movaps		%%xmm2		, 0x020(%%eax)\n\t"\
		"movaps		%%xmm3		, 0x030(%%eax)\n\t"\
		"movaps		%%xmm6		,%%xmm2\n\t"\
		"movaps		%%xmm7		,%%xmm3\n\t"\
		"mulpd		      (%%edi)	,%%xmm6\n\t"\
		"mulpd		      (%%edi)	,%%xmm7\n\t"\
		"mulpd		 0x010(%%edi)	,%%xmm2\n\t"\
		"mulpd		 0x010(%%edi)	,%%xmm3\n\t"\
		"subpd		%%xmm2		,%%xmm7\n\t"\
		"addpd		%%xmm3		,%%xmm6\n\t"\
		"movaps		%%xmm7		, 0x030(%%ebx)\n\t"\
		"movaps		%%xmm6		, 0x020(%%ebx)\n\t"\
		"movaps		 0x020(%%eax)	,%%xmm6\n\t"\
		"movaps		 0x030(%%eax)	,%%xmm7\n\t"\
		"movaps		%%xmm6		,%%xmm2\n\t"\
		"movaps		%%xmm7		,%%xmm3\n\t"\
		"mulpd		 0x20 (%%edi)	,%%xmm6\n\t"\
		"mulpd		 0x20 (%%edi)	,%%xmm7\n\t"\
		"mulpd		 0x30 (%%edi)	,%%xmm2\n\t"\
		"mulpd		 0x30 (%%edi)	,%%xmm3\n\t"\
		"subpd		%%xmm2		,%%xmm7\n\t"\
		"addpd		%%xmm3		,%%xmm6\n\t"\
		"movaps		%%xmm7		, 0x130(%%ebx)\n\t"\
		"movaps		%%xmm6		, 0x120(%%ebx)\n\t"\
		"addl		$0x040		,%%edi\n\t"\
		"subpd		%%xmm5		,%%xmm0\n\t"\
		"subpd		%%xmm4		,%%xmm1\n\t"\
		"addpd		%%xmm5		,%%xmm5\n\t"\
		"addpd		%%xmm4		,%%xmm4\n\t"\
		"addpd		%%xmm0		,%%xmm5\n\t"\
		"addpd		%%xmm1		,%%xmm4\n\t"\
		"movaps		%%xmm5		,%%xmm2\n\t"\
		"movaps		%%xmm1		,%%xmm3\n\t"\
		"mulpd		      (%%edi)	,%%xmm5\n\t"\
		"mulpd		      (%%edi)	,%%xmm1\n\t"\
		"mulpd		 0x010(%%edi)	,%%xmm2\n\t"\
		"mulpd		 0x010(%%edi)	,%%xmm3\n\t"\
		"subpd		%%xmm2		,%%xmm1\n\t"\
		"addpd		%%xmm3		,%%xmm5\n\t"\
		"movaps		%%xmm1		, 0x0b0(%%ebx)\n\t"\
		"movaps		%%xmm5		, 0x0a0(%%ebx)\n\t"\
		"movaps		%%xmm0		,%%xmm2\n\t"\
		"movaps		%%xmm4		,%%xmm3\n\t"\
		"mulpd		 0x20 (%%edi)	,%%xmm0\n\t"\
		"mulpd		 0x20 (%%edi)	,%%xmm4\n\t"\
		"mulpd		 0x30 (%%edi)	,%%xmm2\n\t"\
		"mulpd		 0x30 (%%edi)	,%%xmm3\n\t"\
		"subpd		%%xmm2		,%%xmm4\n\t"\
		"addpd		%%xmm3		,%%xmm0\n\t"\
		"movaps		%%xmm4		, 0x1b0(%%ebx)\n\t"\
		"movaps		%%xmm0		, 0x1a0(%%ebx)\n\t"\
		"/************************************************/\n\t"\
		"/*...Block 8: t0E,t1E,t2E,t3E -> r38,3C,3A,3E:	*/\n\t"\
		"/************************************************/\n\t"\
		"addl		$0x080		,%%eax\n\t"\
		"movaps		 0x020(%%eax)	,%%xmm4\n\t"\
		"movaps		 0x060(%%eax)	,%%xmm0\n\t"\
		"movaps		 0x030(%%eax)	,%%xmm5\n\t"\
		"movaps		 0x070(%%eax)	,%%xmm1\n\t"\
		"movaps		 0x020(%%eax)	,%%xmm6\n\t"\
		"movaps		 0x060(%%eax)	,%%xmm2\n\t"\
		"movaps		 0x030(%%eax)	,%%xmm7\n\t"\
		"movaps		 0x070(%%eax)	,%%xmm3\n\t"\
		"mulpd		 0x010(%%edx)	,%%xmm4\n\t"\
		"mulpd		 0x30 (%%edx)	,%%xmm0\n\t"\
		"mulpd		 0x010(%%edx)	,%%xmm5\n\t"\
		"mulpd		 0x30 (%%edx)	,%%xmm1\n\t"\
		"mulpd		      (%%edx)	,%%xmm6\n\t"\
		"mulpd		 0x20 (%%edx)	,%%xmm2\n\t"\
		"mulpd		      (%%edx)	,%%xmm7\n\t"\
		"mulpd		 0x20 (%%edx)	,%%xmm3\n\t"\
		"subpd		%%xmm6		,%%xmm5\n\t"\
		"subpd		%%xmm2		,%%xmm1\n\t"\
		"addpd		%%xmm7		,%%xmm4\n\t"\
		"addpd		%%xmm3		,%%xmm0\n\t"\
		"movaps		%%xmm5		,%%xmm7\n\t"\
		"movaps		%%xmm4		,%%xmm6\n\t"\
		"addpd		%%xmm0		,%%xmm4\n\t"\
		"addpd		%%xmm1		,%%xmm5\n\t"\
		"subpd		%%xmm0		,%%xmm6\n\t"\
		"subpd		%%xmm1		,%%xmm7\n\t"\
		"movaps		 0x040(%%eax)	,%%xmm2\n\t"\
		"movaps		 0x050(%%eax)	,%%xmm3\n\t"\
		"movaps		 0x040(%%eax)	,%%xmm0\n\t"\
		"movaps		 0x050(%%eax)	,%%xmm1\n\t"\
		"mulpd		      (%%ecx)	,%%xmm2\n\t"\
		"mulpd		 0x010(%%ecx)	,%%xmm1\n\t"\
		"mulpd		      (%%ecx)	,%%xmm3\n\t"\
		"mulpd		 0x010(%%ecx)	,%%xmm0\n\t"\
		"subpd		%%xmm1		,%%xmm2\n\t"\
		"addpd		%%xmm0		,%%xmm3\n\t"\
		"movaps		      (%%eax)	,%%xmm0\n\t"\
		"movaps		 0x010(%%eax)	,%%xmm1\n\t"\
		"subpd		%%xmm2		,%%xmm0\n\t"\
		"subpd		%%xmm3		,%%xmm1\n\t"\
		"addpd		%%xmm2		,%%xmm2\n\t"\
		"addpd		%%xmm3		,%%xmm3\n\t"\
		"addpd		%%xmm0		,%%xmm2\n\t"\
		"addpd		%%xmm1		,%%xmm3\n\t"\
		"addl		$0x040		,%%edi\n\t"\
		"subpd		%%xmm6		,%%xmm0\n\t"\
		"subpd		%%xmm7		,%%xmm1\n\t"\
		"addpd		%%xmm6		,%%xmm6\n\t"\
		"addpd		%%xmm7		,%%xmm7\n\t"\
		"addpd		%%xmm0		,%%xmm6\n\t"\
		"addpd		%%xmm1		,%%xmm7\n\t"\
		"movaps		%%xmm0		, 0x020(%%eax)\n\t"\
		"movaps		%%xmm1		, 0x030(%%eax)\n\t"\
		"movaps		%%xmm6		,%%xmm0\n\t"\
		"movaps		%%xmm7		,%%xmm1\n\t"\
		"mulpd		      (%%edi)	,%%xmm6\n\t"\
		"mulpd		      (%%edi)	,%%xmm7\n\t"\
		"mulpd		 0x010(%%edi)	,%%xmm0\n\t"\
		"mulpd		 0x010(%%edi)	,%%xmm1\n\t"\
		"subpd		%%xmm0		,%%xmm7\n\t"\
		"addpd		%%xmm1		,%%xmm6\n\t"\
		"movaps		%%xmm7		, 0x070(%%ebx)\n\t"\
		"movaps		%%xmm6		, 0x060(%%ebx)\n\t"\
		"movaps		 0x020(%%eax)	,%%xmm6\n\t"\
		"movaps		 0x030(%%eax)	,%%xmm7\n\t"\
		"movaps		%%xmm6		,%%xmm0\n\t"\
		"movaps		%%xmm7		,%%xmm1\n\t"\
		"mulpd		 0x20 (%%edi)	,%%xmm6\n\t"\
		"mulpd		 0x20 (%%edi)	,%%xmm7\n\t"\
		"mulpd		 0x30 (%%edi)	,%%xmm0\n\t"\
		"mulpd		 0x30 (%%edi)	,%%xmm1\n\t"\
		"subpd		%%xmm0		,%%xmm7\n\t"\
		"addpd		%%xmm1		,%%xmm6\n\t"\
		"movaps		%%xmm7		, 0x170(%%ebx)\n\t"\
		"movaps		%%xmm6		, 0x160(%%ebx)\n\t"\
		"addl		$0x040		,%%edi\n\t"\
		"subpd		%%xmm5		,%%xmm2\n\t"\
		"subpd		%%xmm4		,%%xmm3\n\t"\
		"addpd		%%xmm5		,%%xmm5\n\t"\
		"addpd		%%xmm4		,%%xmm4\n\t"\
		"addpd		%%xmm2		,%%xmm5\n\t"\
		"addpd		%%xmm3		,%%xmm4\n\t"\
		"movaps		%%xmm5		,%%xmm0\n\t"\
		"movaps		%%xmm3		,%%xmm1\n\t"\
		"mulpd		      (%%edi)	,%%xmm5\n\t"\
		"mulpd		      (%%edi)	,%%xmm3\n\t"\
		"mulpd		 0x010(%%edi)	,%%xmm0\n\t"\
		"mulpd		 0x010(%%edi)	,%%xmm1\n\t"\
		"subpd		%%xmm0		,%%xmm3\n\t"\
		"addpd		%%xmm1		,%%xmm5\n\t"\
		"movaps		%%xmm3		, 0x0f0(%%ebx)\n\t"\
		"movaps		%%xmm5		, 0x0e0(%%ebx)\n\t"\
		"movaps		%%xmm2		,%%xmm0\n\t"\
		"movaps		%%xmm4		,%%xmm1\n\t"\
		"mulpd		 0x20 (%%edi)	,%%xmm2\n\t"\
		"mulpd		 0x20 (%%edi)	,%%xmm4\n\t"\
		"mulpd		 0x30 (%%edi)	,%%xmm0\n\t"\
		"mulpd		 0x30 (%%edi)	,%%xmm1\n\t"\
		"subpd		%%xmm0		,%%xmm4\n\t"\
		"addpd		%%xmm1		,%%xmm2\n\t"\
		"movaps		%%xmm4		, 0x1f0(%%ebx)\n\t"\
		"movaps		%%xmm2		, 0x1e0(%%ebx)\n\t"\
		"/************************************************/\n\t"\
		"/*...Block 1: t00,t10,t20,t30 -> r00,04,02,06:	*/\n\t"\
		"/************************************************/\n\t"\
		"movl		%[__r00]		,%%edx\n\t"\
		"movaps		      (%%edx)	,%%xmm0\n\t"\
		"movaps		 0x010(%%edx)	,%%xmm1\n\t"\
		"movaps		 0x040(%%edx)	,%%xmm2\n\t"\
		"movaps		 0x050(%%edx)	,%%xmm3\n\t"\
		"subpd		 0x040(%%edx)	,%%xmm0\n\t"\
		"subpd		 0x050(%%edx)	,%%xmm1\n\t"\
		"addpd		      (%%edx)	,%%xmm2\n\t"\
		"addpd		 0x010(%%edx)	,%%xmm3\n\t"\
		"movaps		 0x020(%%edx)	,%%xmm4\n\t"\
		"movaps		 0x030(%%edx)	,%%xmm5\n\t"\
		"movaps		 0x060(%%edx)	,%%xmm6\n\t"\
		"movaps		 0x070(%%edx)	,%%xmm7\n\t"\
		"subpd		 0x060(%%edx)	,%%xmm4\n\t"\
		"subpd		 0x070(%%edx)	,%%xmm5\n\t"\
		"addpd		 0x020(%%edx)	,%%xmm6\n\t"\
		"addpd		 0x030(%%edx)	,%%xmm7\n\t"\
		"movl		%[__add0]		,%%eax\n\t"\
		"movl		%[__c10]		,%%ecx\n\t"\
		"addpd		%%xmm6		,%%xmm2\n\t"\
		"addpd		%%xmm7		,%%xmm3\n\t"\
		"movaps		%%xmm2		,      (%%edx)\n\t"\
		"movaps		%%xmm3		, 0x010(%%edx)\n\t"\
		"addpd		%%xmm6		,%%xmm6\n\t"\
		"addpd		%%xmm7		,%%xmm7\n\t"\
		"subpd		%%xmm6		,%%xmm2\n\t"\
		"subpd		%%xmm7		,%%xmm3\n\t"\
		"movaps		%%xmm2		,%%xmm6\n\t"\
		"movaps		%%xmm3		,%%xmm7\n\t"\
		"mulpd		      (%%ecx)	,%%xmm2\n\t"\
		"mulpd		      (%%ecx)	,%%xmm3\n\t"\
		"mulpd		 0x010(%%ecx)	,%%xmm6\n\t"\
		"mulpd		 0x010(%%ecx)	,%%xmm7\n\t"\
		"subpd		%%xmm6		,%%xmm3\n\t"\
		"addpd		%%xmm7		,%%xmm2\n\t"\
		"movl		%[__add1]		,%%ebx\n\t"\
		"movaps		%%xmm3		,%%xmm7\n\t"\
		"movaps		%%xmm2		,%%xmm6\n\t"\
		"unpckhpd	 0x110(%%ebx)	,%%xmm7\n\t"\
		"unpcklpd	 0x110(%%ebx)	,%%xmm3\n\t"\
		"movaps		%%xmm7		, 0x110(%%ebx)\n\t"\
		"unpckhpd	 0x100(%%ebx)	,%%xmm6\n\t"\
		"unpcklpd	 0x100(%%ebx)	,%%xmm2\n\t"\
		"movaps		%%xmm6		, 0x100(%%ebx)\n\t"\
		"movaps		%%xmm3		, 0x110(%%eax)\n\t"\
		"movaps		%%xmm2		, 0x100(%%eax)\n\t"\
		"movaps		 0x010(%%edx)	,%%xmm3\n\t"\
		"movaps		      (%%edx)	,%%xmm2\n\t"\
		"movaps		%%xmm3		,%%xmm7\n\t"\
		"movaps		%%xmm2		,%%xmm6\n\t"\
		"unpckhpd	 0x010(%%ebx)	,%%xmm7\n\t"\
		"unpcklpd	 0x010(%%ebx)	,%%xmm3\n\t"\
		"movaps		%%xmm7		, 0x010(%%ebx)\n\t"\
		"unpckhpd	      (%%ebx)	,%%xmm6\n\t"\
		"unpcklpd	      (%%ebx)	,%%xmm2\n\t"\
		"movaps		%%xmm6		,      (%%ebx)\n\t"\
		"movaps		%%xmm3		, 0x010(%%eax)\n\t"\
		"movaps		%%xmm2		,      (%%eax)\n\t"\
		"movl		%[__c08]		,%%ecx\n\t"\
		"addpd		%%xmm5		,%%xmm0\n\t"\
		"subpd		%%xmm4		,%%xmm1\n\t"\
		"movaps		%%xmm0		,%%xmm2\n\t"\
		"movaps		%%xmm1		,%%xmm3\n\t"\
		"addpd		%%xmm5		,%%xmm5\n\t"\
		"addpd		%%xmm4		,%%xmm4\n\t"\
		"movaps		%%xmm0		,%%xmm6\n\t"\
		"movaps		%%xmm1		,%%xmm7\n\t"\
		"mulpd		      (%%ecx)	,%%xmm2\n\t"\
		"mulpd		      (%%ecx)	,%%xmm3\n\t"\
		"mulpd		 0x010(%%ecx)	,%%xmm6\n\t"\
		"mulpd		 0x010(%%ecx)	,%%xmm7\n\t"\
		"subpd		%%xmm6		,%%xmm3\n\t"\
		"addpd		%%xmm7		,%%xmm2\n\t"\
		"movaps		%%xmm3		,%%xmm7\n\t"\
		"movaps		%%xmm2		,%%xmm6\n\t"\
		"unpckhpd	 0x090(%%ebx)	,%%xmm7\n\t"\
		"unpcklpd	 0x090(%%ebx)	,%%xmm3\n\t"\
		"movaps		%%xmm7		, 0x090(%%ebx)\n\t"\
		"unpckhpd	 0x080(%%ebx)	,%%xmm6\n\t"\
		"unpcklpd	 0x080(%%ebx)	,%%xmm2\n\t"\
		"movaps		%%xmm6		, 0x080(%%ebx)\n\t"\
		"movaps		%%xmm3		, 0x090(%%eax)\n\t"\
		"movaps		%%xmm2		, 0x080(%%eax)\n\t"\
		"movl		%[__c18]		,%%ecx\n\t"\
		"subpd		%%xmm5		,%%xmm0\n\t"\
		"addpd		%%xmm4		,%%xmm1\n\t"\
		"movaps		%%xmm0		,%%xmm6\n\t"\
		"movaps		%%xmm1		,%%xmm7\n\t"\
		"mulpd		      (%%ecx)	,%%xmm0\n\t"\
		"mulpd		      (%%ecx)	,%%xmm1\n\t"\
		"mulpd		 0x010(%%ecx)	,%%xmm6\n\t"\
		"mulpd		 0x010(%%ecx)	,%%xmm7\n\t"\
		"subpd		%%xmm6		,%%xmm1\n\t"\
		"addpd		%%xmm7		,%%xmm0\n\t"\
		"movaps		%%xmm1		,%%xmm7\n\t"\
		"movaps		%%xmm0		,%%xmm6\n\t"\
		"unpckhpd	 0x190(%%ebx)	,%%xmm7\n\t"\
		"unpcklpd	 0x190(%%ebx)	,%%xmm1\n\t"\
		"movaps		%%xmm7		, 0x190(%%ebx)\n\t"\
		"unpckhpd	 0x180(%%ebx)	,%%xmm6\n\t"\
		"unpcklpd	 0x180(%%ebx)	,%%xmm0\n\t"\
		"movaps		%%xmm6		, 0x180(%%ebx)\n\t"\
		"movaps		%%xmm1		, 0x190(%%eax)\n\t"\
		"movaps		%%xmm0		, 0x180(%%eax)\n\t"\
		"/************************************************/\n\t"\
		"/*...Block 5: t08,t18,t28,t38 -> r08,0C,0A,0E	*/\n\t"\
		"/************************************************/\n\t"\
		"movl		%[__r08]		,%%edx\n\t"\
		"movaps		 (%%esi)	,%%xmm2\n\t	/* isrt2 */"\
		"movaps		 0x020(%%edx)	,%%xmm4\n\t"\
		"movaps		 0x030(%%edx)	,%%xmm5\n\t"\
		"movaps		 0x060(%%edx)	,%%xmm0\n\t"\
		"movaps		 0x070(%%edx)	,%%xmm1\n\t"\
		"addpd		 0x030(%%edx)	,%%xmm4\n\t"\
		"subpd		 0x020(%%edx)	,%%xmm5\n\t"\
		"subpd		 0x070(%%edx)	,%%xmm0\n\t"\
		"addpd		 0x060(%%edx)	,%%xmm1\n\t"\
		"mulpd		%%xmm2		,%%xmm4\n\t"\
		"mulpd		%%xmm2		,%%xmm5\n\t"\
		"mulpd		%%xmm2		,%%xmm0\n\t"\
		"mulpd		%%xmm2		,%%xmm1\n\t"\
		"movaps		%%xmm4		,%%xmm6\n\t"\
		"movaps		%%xmm5		,%%xmm7\n\t"\
		"subpd		%%xmm0		,%%xmm4\n\t"\
		"subpd		%%xmm1		,%%xmm5\n\t"\
		"addpd		%%xmm0		,%%xmm6\n\t"\
		"addpd		%%xmm1		,%%xmm7\n\t"\
		"movaps		      (%%edx)	,%%xmm0\n\t"\
		"movaps		 0x010(%%edx)	,%%xmm1\n\t"\
		"movaps		 0x040(%%edx)	,%%xmm2\n\t"\
		"movaps		 0x050(%%edx)	,%%xmm3\n\t"\
		"subpd		 0x050(%%edx)	,%%xmm0\n\t"\
		"subpd		 0x040(%%edx)	,%%xmm1\n\t"\
		"addpd		      (%%edx)	,%%xmm3\n\t"\
		"addpd		 0x010(%%edx)	,%%xmm2\n\t"\
		"movl		%[__add0]		,%%eax\n\t"\
		"movl		%[__c04]		,%%ecx\n\t"\
		"subpd		%%xmm4		,%%xmm3\n\t"\
		"subpd		%%xmm5		,%%xmm1\n\t"\
		"addpd		%%xmm4		,%%xmm4\n\t"\
		"addpd		%%xmm5		,%%xmm5\n\t"\
		"addpd		%%xmm3		,%%xmm4\n\t"\
		"addpd		%%xmm1		,%%xmm5\n\t"\
		"movaps		%%xmm3		,      (%%edx)\n\t"\
		"movaps		%%xmm1		, 0x010(%%edx)\n\t"\
		"movaps		%%xmm4		,%%xmm3\n\t"\
		"movaps		%%xmm5		,%%xmm1\n\t"\
		"mulpd		      (%%ecx)	,%%xmm4\n\t"\
		"mulpd		      (%%ecx)	,%%xmm5\n\t"\
		"mulpd		 0x010(%%ecx)	,%%xmm3\n\t"\
		"mulpd		 0x010(%%ecx)	,%%xmm1\n\t"\
		"subpd		%%xmm3		,%%xmm5\n\t"\
		"addpd		%%xmm1		,%%xmm4\n\t"\
		"movl		%[__add1]		,%%ebx\n\t"\
		"movaps		%%xmm5		,%%xmm3\n\t"\
		"movaps		%%xmm4		,%%xmm1\n\t"\
		"unpckhpd	 0x050(%%ebx)	,%%xmm3\n\t"\
		"unpcklpd	 0x050(%%ebx)	,%%xmm5\n\t"\
		"movaps		%%xmm3		, 0x050(%%ebx)\n\t"\
		"unpckhpd	 0x040(%%ebx)	,%%xmm1\n\t"\
		"unpcklpd	 0x040(%%ebx)	,%%xmm4\n\t"\
		"movaps		%%xmm1		, 0x040(%%ebx)\n\t"\
		"movaps		%%xmm5		, 0x050(%%eax)\n\t"\
		"movaps		%%xmm4		, 0x040(%%eax)\n\t"\
		"movl		%[__c14]		,%%ecx\n\t"\
		"movaps		      (%%edx)	,%%xmm4\n\t"\
		"movaps		 0x010(%%edx)	,%%xmm5\n\t"\
		"movaps		%%xmm4		,%%xmm3\n\t"\
		"movaps		%%xmm5		,%%xmm1\n\t"\
		"mulpd		      (%%ecx)	,%%xmm4\n\t"\
		"mulpd		      (%%ecx)	,%%xmm5\n\t"\
		"mulpd		 0x010(%%ecx)	,%%xmm3\n\t"\
		"mulpd		 0x010(%%ecx)	,%%xmm1\n\t"\
		"subpd		%%xmm3		,%%xmm5\n\t"\
		"addpd		%%xmm1		,%%xmm4\n\t"\
		"movaps		%%xmm5		,%%xmm3\n\t"\
		"movaps		%%xmm4		,%%xmm1\n\t"\
		"unpckhpd	 0x150(%%ebx)	,%%xmm3\n\t"\
		"unpcklpd	 0x150(%%ebx)	,%%xmm5\n\t"\
		"movaps		%%xmm3		, 0x150(%%ebx)\n\t"\
		"unpckhpd	 0x140(%%ebx)	,%%xmm1\n\t"\
		"unpcklpd	 0x140(%%ebx)	,%%xmm4\n\t"\
		"movaps		%%xmm1		, 0x140(%%ebx)\n\t"\
		"movaps		%%xmm5		, 0x150(%%eax)\n\t"\
		"movaps		%%xmm4		, 0x140(%%eax)\n\t"\
		"movl		%[__c0C]		,%%ecx\n\t"\
		"subpd		%%xmm7		,%%xmm0\n\t"\
		"subpd		%%xmm6		,%%xmm2\n\t"\
		"addpd		%%xmm7		,%%xmm7\n\t"\
		"addpd		%%xmm6		,%%xmm6\n\t"\
		"addpd		%%xmm0		,%%xmm7\n\t"\
		"addpd		%%xmm2		,%%xmm6\n\t"\
		"movaps		%%xmm7		,%%xmm4\n\t"\
		"movaps		%%xmm2		,%%xmm5\n\t"\
		"mulpd		      (%%ecx)	,%%xmm7\n\t"\
		"mulpd		      (%%ecx)	,%%xmm2\n\t"\
		"mulpd		 0x010(%%ecx)	,%%xmm4\n\t"\
		"mulpd		 0x010(%%ecx)	,%%xmm5\n\t"\
		"subpd		%%xmm4		,%%xmm2\n\t"\
		"addpd		%%xmm5		,%%xmm7\n\t"\
		"movaps		%%xmm2		,%%xmm5\n\t"\
		"movaps		%%xmm7		,%%xmm4\n\t"\
		"unpckhpd	 0x0d0(%%ebx)	,%%xmm5\n\t"\
		"unpcklpd	 0x0d0(%%ebx)	,%%xmm2\n\t"\
		"movaps		%%xmm5		, 0x0d0(%%ebx)\n\t"\
		"unpckhpd	 0x0c0(%%ebx)	,%%xmm4\n\t"\
		"unpcklpd	 0x0c0(%%ebx)	,%%xmm7\n\t"\
		"movaps		%%xmm4		, 0x0c0(%%ebx)\n\t"\
		"movaps		%%xmm2		, 0x0d0(%%eax)\n\t"\
		"movaps		%%xmm7		, 0x0c0(%%eax)\n\t"\
		"movl		%[__c1C]		,%%ecx\n\t"\
		"movaps		%%xmm0		,%%xmm4\n\t"\
		"movaps		%%xmm6		,%%xmm5\n\t"\
		"mulpd		      (%%ecx)	,%%xmm0\n\t"\
		"mulpd		      (%%ecx)	,%%xmm6\n\t"\
		"mulpd		 0x010(%%ecx)	,%%xmm4\n\t"\
		"mulpd		 0x010(%%ecx)	,%%xmm5\n\t"\
		"subpd		%%xmm4		,%%xmm6\n\t"\
		"addpd		%%xmm5		,%%xmm0\n\t"\
		"movaps		%%xmm6		,%%xmm5\n\t"\
		"movaps		%%xmm0		,%%xmm4\n\t"\
		"unpckhpd	 0x1d0(%%ebx)	,%%xmm5\n\t"\
		"unpcklpd	 0x1d0(%%ebx)	,%%xmm6\n\t"\
		"movaps		%%xmm5		, 0x1d0(%%ebx)\n\t"\
		"unpckhpd	 0x1c0(%%ebx)	,%%xmm4\n\t"\
		"unpcklpd	 0x1c0(%%ebx)	,%%xmm0\n\t"\
		"movaps		%%xmm4		, 0x1c0(%%ebx)\n\t"\
		"movaps		%%xmm6		, 0x1d0(%%eax)\n\t"\
		"movaps		%%xmm0		, 0x1c0(%%eax)\n\t"\
		"/************************************************/\n\t"\
		"/*...Block 3: t04,t14,t24,t34 -> r20,24,22,26	*/\n\t"\
		"/************************************************/\n\t"\
		"movl		%[__r20]		,%%edx\n\t"\
		"movl		%%esi			,%%ecx\n\t"\
		"addl		$0x010			,%%ecx	/* cc0 */\n\t"\
		"movaps		 0x020(%%edx)	,%%xmm4\n\t"\
		"movaps		 0x060(%%edx)	,%%xmm0\n\t"\
		"movaps		 0x030(%%edx)	,%%xmm5\n\t"\
		"movaps		 0x070(%%edx)	,%%xmm1\n\t"\
		"movaps		 0x020(%%edx)	,%%xmm6\n\t"\
		"movaps		 0x060(%%edx)	,%%xmm2\n\t"\
		"movaps		 0x030(%%edx)	,%%xmm7\n\t"\
		"movaps		 0x070(%%edx)	,%%xmm3\n\t"\
		"mulpd		      (%%ecx)	,%%xmm4\n\t"\
		"mulpd		 0x010(%%ecx)	,%%xmm0\n\t"\
		"mulpd		      (%%ecx)	,%%xmm5\n\t"\
		"mulpd		 0x010(%%ecx)	,%%xmm1\n\t"\
		"mulpd		 0x010(%%ecx)	,%%xmm6\n\t"\
		"mulpd		      (%%ecx)	,%%xmm2\n\t"\
		"mulpd		 0x010(%%ecx)	,%%xmm7\n\t"\
		"mulpd		      (%%ecx)	,%%xmm3\n\t"\
		"subpd		%%xmm6		,%%xmm5\n\t"\
		"subpd		%%xmm2		,%%xmm1\n\t"\
		"addpd		%%xmm7		,%%xmm4\n\t"\
		"addpd		%%xmm3		,%%xmm0\n\t"\
		"movaps		%%xmm5		,%%xmm7\n\t"\
		"movaps		%%xmm4		,%%xmm6\n\t"\
		"addpd		%%xmm0		,%%xmm4\n\t"\
		"addpd		%%xmm1		,%%xmm5\n\t"\
		"subpd		%%xmm0		,%%xmm6\n\t"\
		"subpd		%%xmm1		,%%xmm7\n\t"\
		"movaps		 0x040(%%edx)	,%%xmm2\n\t"\
		"movaps		 0x050(%%edx)	,%%xmm3\n\t"\
		"movaps		      (%%edx)	,%%xmm0\n\t"\
		"movaps		 0x010(%%edx)	,%%xmm1\n\t"\
		"addpd		 0x050(%%edx)	,%%xmm2\n\t"\
		"subpd		 0x040(%%edx)	,%%xmm3\n\t"\
		"mulpd		      (%%esi)	,%%xmm2\n\t"\
		"mulpd		      (%%esi)	,%%xmm3\n\t"\
		"subpd		%%xmm2		,%%xmm0\n\t"\
		"subpd		%%xmm3		,%%xmm1\n\t"\
		"addpd		%%xmm2		,%%xmm2\n\t"\
		"addpd		%%xmm3		,%%xmm3\n\t"\
		"addpd		%%xmm0		,%%xmm2\n\t"\
		"addpd		%%xmm1		,%%xmm3\n\t"\
		"movl		%[__add0]		,%%eax\n\t"\
		"movl		%[__c02]		,%%ecx\n\t"\
		"subpd		%%xmm4		,%%xmm2\n\t"\
		"subpd		%%xmm5		,%%xmm3\n\t"\
		"addpd		%%xmm4		,%%xmm4\n\t"\
		"addpd		%%xmm5		,%%xmm5\n\t"\
		"addpd		%%xmm2		,%%xmm4\n\t"\
		"addpd		%%xmm3		,%%xmm5\n\t"\
		"movaps		%%xmm2		,      (%%edx)\n\t"\
		"movaps		%%xmm3		, 0x010(%%edx)\n\t"\
		"movaps		%%xmm4		,%%xmm2\n\t"\
		"movaps		%%xmm5		,%%xmm3\n\t"\
		"mulpd		      (%%ecx)	,%%xmm4\n\t"\
		"mulpd		      (%%ecx)	,%%xmm5\n\t"\
		"mulpd		 0x010(%%ecx)	,%%xmm2\n\t"\
		"mulpd		 0x010(%%ecx)	,%%xmm3\n\t"\
		"subpd		%%xmm2		,%%xmm5\n\t"\
		"addpd		%%xmm3		,%%xmm4\n\t"\
		"movl		%[__add1]		,%%ebx\n\t"\
		"movaps		%%xmm5		,%%xmm3\n\t"\
		"movaps		%%xmm4		,%%xmm2\n\t"\
		"unpckhpd	 0x030(%%ebx)	,%%xmm3\n\t"\
		"unpcklpd	 0x030(%%ebx)	,%%xmm5\n\t"\
		"movaps		%%xmm3		, 0x030(%%ebx)\n\t"\
		"unpckhpd	 0x020(%%ebx)	,%%xmm2\n\t"\
		"unpcklpd	 0x020(%%ebx)	,%%xmm4\n\t"\
		"movaps		%%xmm2		, 0x020(%%ebx)\n\t"\
		"movaps		%%xmm5		, 0x030(%%eax)\n\t"\
		"movaps		%%xmm4		, 0x020(%%eax)\n\t"\
		"movl		%[__c12]		,%%ecx\n\t"\
		"movaps		      (%%edx)	,%%xmm4\n\t"\
		"movaps		 0x010(%%edx)	,%%xmm5\n\t"\
		"movaps		%%xmm4		,%%xmm2\n\t"\
		"movaps		%%xmm5		,%%xmm3\n\t"\
		"mulpd		      (%%ecx)	,%%xmm4\n\t"\
		"mulpd		      (%%ecx)	,%%xmm5\n\t"\
		"mulpd		 0x010(%%ecx)	,%%xmm2\n\t"\
		"mulpd		 0x010(%%ecx)	,%%xmm3\n\t"\
		"subpd		%%xmm2		,%%xmm5\n\t"\
		"addpd		%%xmm3		,%%xmm4\n\t"\
		"movaps		%%xmm5		,%%xmm3\n\t"\
		"movaps		%%xmm4		,%%xmm2\n\t"\
		"unpckhpd	 0x130(%%ebx)	,%%xmm3\n\t"\
		"unpcklpd	 0x130(%%ebx)	,%%xmm5\n\t"\
		"movaps		%%xmm3		, 0x130(%%ebx)\n\t"\
		"unpckhpd	 0x120(%%ebx)	,%%xmm2\n\t"\
		"unpcklpd	 0x120(%%ebx)	,%%xmm4\n\t"\
		"movaps		%%xmm2		, 0x120(%%ebx)\n\t"\
		"movaps		%%xmm5		, 0x130(%%eax)\n\t"\
		"movaps		%%xmm4		, 0x120(%%eax)\n\t"\
		"movl		%[__c0A]		,%%ecx\n\t"\
		"subpd		%%xmm7		,%%xmm0\n\t"\
		"subpd		%%xmm6		,%%xmm1\n\t"\
		"addpd		%%xmm7		,%%xmm7\n\t"\
		"addpd		%%xmm6		,%%xmm6\n\t"\
		"addpd		%%xmm0		,%%xmm7\n\t"\
		"addpd		%%xmm1		,%%xmm6\n\t"\
		"movaps		%%xmm7		,%%xmm4\n\t"\
		"movaps		%%xmm1		,%%xmm5\n\t"\
		"mulpd		      (%%ecx)	,%%xmm7\n\t"\
		"mulpd		      (%%ecx)	,%%xmm1\n\t"\
		"mulpd		 0x010(%%ecx)	,%%xmm4\n\t"\
		"mulpd		 0x010(%%ecx)	,%%xmm5\n\t"\
		"subpd		%%xmm4		,%%xmm1\n\t"\
		"addpd		%%xmm5		,%%xmm7\n\t"\
		"movaps		%%xmm1		,%%xmm5\n\t"\
		"movaps		%%xmm7		,%%xmm4\n\t"\
		"unpckhpd	 0x0b0(%%ebx)	,%%xmm5\n\t"\
		"unpcklpd	 0x0b0(%%ebx)	,%%xmm1\n\t"\
		"movaps		%%xmm5		, 0x0b0(%%ebx)\n\t"\
		"unpckhpd	 0x0a0(%%ebx)	,%%xmm4\n\t"\
		"unpcklpd	 0x0a0(%%ebx)	,%%xmm7\n\t"\
		"movaps		%%xmm4		, 0x0a0(%%ebx)\n\t"\
		"movaps		%%xmm1		, 0x0b0(%%eax)\n\t"\
		"movaps		%%xmm7		, 0x0a0(%%eax)\n\t"\
		"movl		%[__c1A]		,%%ecx\n\t"\
		"movaps		%%xmm0		,%%xmm4\n\t"\
		"movaps		%%xmm6		,%%xmm5\n\t"\
		"mulpd		      (%%ecx)	,%%xmm0\n\t"\
		"mulpd		      (%%ecx)	,%%xmm6\n\t"\
		"mulpd		 0x010(%%ecx)	,%%xmm4\n\t"\
		"mulpd		 0x010(%%ecx)	,%%xmm5\n\t"\
		"subpd		%%xmm4		,%%xmm6\n\t"\
		"addpd		%%xmm5		,%%xmm0\n\t"\
		"movaps		%%xmm6		,%%xmm5\n\t"\
		"movaps		%%xmm0		,%%xmm4\n\t"\
		"unpckhpd	 0x1b0(%%ebx)	,%%xmm5\n\t"\
		"unpcklpd	 0x1b0(%%ebx)	,%%xmm6\n\t"\
		"movaps		%%xmm5		, 0x1b0(%%ebx)\n\t"\
		"unpckhpd	 0x1a0(%%ebx)	,%%xmm4\n\t"\
		"unpcklpd	 0x1a0(%%ebx)	,%%xmm0\n\t"\
		"movaps		%%xmm4		, 0x1a0(%%ebx)\n\t"\
		"movaps		%%xmm6		, 0x1b0(%%eax)\n\t"\
		"movaps		%%xmm0		, 0x1a0(%%eax)\n\t"\
		"/************************************************/\n\t"\
		"/*...Block 7: t0C,t1C,t2C,t3C -> r28,2C,2A,2E	*/\n\t"\
		"/************************************************/\n\t"\
		"movl		%[__r28]		,%%edx\n\t"\
		"movl		%%esi			,%%ecx\n\t"\
		"addl		$0x010			,%%ecx	/* cc0 */\n\t"\
		"movaps		 0x020(%%edx)	,%%xmm4\n\t"\
		"movaps		 0x060(%%edx)	,%%xmm0\n\t"\
		"movaps		 0x030(%%edx)	,%%xmm5\n\t"\
		"movaps		 0x070(%%edx)	,%%xmm1\n\t"\
		"movaps		 0x020(%%edx)	,%%xmm6\n\t"\
		"movaps		 0x060(%%edx)	,%%xmm2\n\t"\
		"movaps		 0x030(%%edx)	,%%xmm7\n\t"\
		"movaps		 0x070(%%edx)	,%%xmm3\n\t"\
		"mulpd		 0x010(%%ecx)	,%%xmm4\n\t"\
		"mulpd		      (%%ecx)	,%%xmm0\n\t"\
		"mulpd		 0x010(%%ecx)	,%%xmm5\n\t"\
		"mulpd		      (%%ecx)	,%%xmm1\n\t"\
		"mulpd		      (%%ecx)	,%%xmm6\n\t"\
		"mulpd		 0x010(%%ecx)	,%%xmm2\n\t"\
		"mulpd		      (%%ecx)	,%%xmm7\n\t"\
		"mulpd		 0x010(%%ecx)	,%%xmm3\n\t"\
		"subpd		%%xmm6		,%%xmm5\n\t"\
		"subpd		%%xmm2		,%%xmm1\n\t"\
		"addpd		%%xmm7		,%%xmm4\n\t"\
		"addpd		%%xmm3		,%%xmm0\n\t"\
		"movaps		%%xmm5		,%%xmm7\n\t"\
		"movaps		%%xmm4		,%%xmm6\n\t"\
		"addpd		%%xmm0		,%%xmm4\n\t"\
		"addpd		%%xmm1		,%%xmm5\n\t"\
		"subpd		%%xmm0		,%%xmm6\n\t"\
		"subpd		%%xmm1		,%%xmm7\n\t"\
		"movaps		 0x040(%%edx)	,%%xmm2\n\t"\
		"movaps		 0x050(%%edx)	,%%xmm3\n\t"\
		"movaps		      (%%edx)	,%%xmm0\n\t"\
		"movaps		 0x010(%%edx)	,%%xmm1\n\t"\
		"subpd		 0x050(%%edx)	,%%xmm2\n\t"\
		"addpd		 0x040(%%edx)	,%%xmm3\n\t"\
		"mulpd		      (%%esi)	,%%xmm2\n\t"\
		"mulpd		      (%%esi)	,%%xmm3\n\t"\
		"subpd		%%xmm2		,%%xmm0\n\t"\
		"subpd		%%xmm3		,%%xmm1\n\t"\
		"addpd		%%xmm2		,%%xmm2\n\t"\
		"addpd		%%xmm3		,%%xmm3\n\t"\
		"addpd		%%xmm0		,%%xmm2\n\t"\
		"addpd		%%xmm1		,%%xmm3\n\t"\
		"movl		%[__add0]		,%%eax\n\t"\
		"movl		%[__c06]		,%%ecx\n\t"\
		"subpd		%%xmm6		,%%xmm0\n\t"\
		"subpd		%%xmm7		,%%xmm1\n\t"\
		"addpd		%%xmm6		,%%xmm6\n\t"\
		"addpd		%%xmm7		,%%xmm7\n\t"\
		"addpd		%%xmm0		,%%xmm6\n\t"\
		"addpd		%%xmm1		,%%xmm7\n\t"\
		"movaps		%%xmm0		,      (%%edx)\n\t"\
		"movaps		%%xmm1		, 0x010(%%edx)\n\t"\
		"movaps		%%xmm6		,%%xmm0\n\t"\
		"movaps		%%xmm7		,%%xmm1\n\t"\
		"mulpd		      (%%ecx)	,%%xmm6\n\t"\
		"mulpd		      (%%ecx)	,%%xmm7\n\t"\
		"mulpd		 0x010(%%ecx)	,%%xmm0\n\t"\
		"mulpd		 0x010(%%ecx)	,%%xmm1\n\t"\
		"subpd		%%xmm0		,%%xmm7\n\t"\
		"addpd		%%xmm1		,%%xmm6\n\t"\
		"movl		%[__add1]		,%%ebx\n\t"\
		"movaps		%%xmm7		,%%xmm1\n\t"\
		"movaps		%%xmm6		,%%xmm0\n\t"\
		"unpckhpd	 0x070(%%ebx)	,%%xmm1\n\t"\
		"unpcklpd	 0x070(%%ebx)	,%%xmm7\n\t"\
		"movaps		%%xmm1		, 0x070(%%ebx)\n\t"\
		"unpckhpd	 0x060(%%ebx)	,%%xmm0\n\t"\
		"unpcklpd	 0x060(%%ebx)	,%%xmm6\n\t"\
		"movaps		%%xmm0		, 0x060(%%ebx)\n\t"\
		"movaps		%%xmm7		, 0x070(%%eax)\n\t"\
		"movaps		%%xmm6		, 0x060(%%eax)\n\t"\
		"movl		%[__c16]		,%%ecx\n\t"\
		"movaps		      (%%edx)	,%%xmm6\n\t"\
		"movaps		 0x010(%%edx)	,%%xmm7\n\t"\
		"movaps		%%xmm6		,%%xmm0\n\t"\
		"movaps		%%xmm7		,%%xmm1\n\t"\
		"mulpd		      (%%ecx)	,%%xmm6\n\t"\
		"mulpd		      (%%ecx)	,%%xmm7\n\t"\
		"mulpd		 0x010(%%ecx)	,%%xmm0\n\t"\
		"mulpd		 0x010(%%ecx)	,%%xmm1\n\t"\
		"subpd		%%xmm0		,%%xmm7\n\t"\
		"addpd		%%xmm1		,%%xmm6\n\t"\
		"movaps		%%xmm7		,%%xmm1\n\t"\
		"movaps		%%xmm6		,%%xmm0\n\t"\
		"unpckhpd	 0x170(%%ebx)	,%%xmm1\n\t"\
		"unpcklpd	 0x170(%%ebx)	,%%xmm7\n\t"\
		"movaps		%%xmm1		, 0x170(%%ebx)\n\t"\
		"unpckhpd	 0x160(%%ebx)	,%%xmm0\n\t"\
		"unpcklpd	 0x160(%%ebx)	,%%xmm6\n\t"\
		"movaps		%%xmm0		, 0x160(%%ebx)\n\t"\
		"movaps		%%xmm7		, 0x170(%%eax)\n\t"\
		"movaps		%%xmm6		, 0x160(%%eax)\n\t"\
		"movl		%[__c0E]		,%%ecx\n\t"\
		"subpd		%%xmm5		,%%xmm2\n\t"\
		"subpd		%%xmm4		,%%xmm3\n\t"\
		"addpd		%%xmm5		,%%xmm5\n\t"\
		"addpd		%%xmm4		,%%xmm4\n\t"\
		"addpd		%%xmm2		,%%xmm5\n\t"\
		"addpd		%%xmm3		,%%xmm4\n\t"\
		"movaps		%%xmm5		,%%xmm0\n\t"\
		"movaps		%%xmm3		,%%xmm1\n\t"\
		"mulpd		      (%%ecx)	,%%xmm5\n\t"\
		"mulpd		      (%%ecx)	,%%xmm3\n\t"\
		"mulpd		 0x010(%%ecx)	,%%xmm0\n\t"\
		"mulpd		 0x010(%%ecx)	,%%xmm1\n\t"\
		"subpd		%%xmm0		,%%xmm3\n\t"\
		"addpd		%%xmm1		,%%xmm5\n\t"\
		"movaps		%%xmm3		,%%xmm1\n\t"\
		"movaps		%%xmm5		,%%xmm0\n\t"\
		"unpckhpd	 0x0f0(%%ebx)	,%%xmm1\n\t"\
		"unpcklpd	 0x0f0(%%ebx)	,%%xmm3\n\t"\
		"movaps		%%xmm1		, 0x0f0(%%ebx)\n\t"\
		"unpckhpd	 0x0e0(%%ebx)	,%%xmm0\n\t"\
		"unpcklpd	 0x0e0(%%ebx)	,%%xmm5\n\t"\
		"movaps		%%xmm0		, 0x0e0(%%ebx)\n\t"\
		"movaps		%%xmm3		, 0x0f0(%%eax)\n\t"\
		"movaps		%%xmm5		, 0x0e0(%%eax)\n\t"\
		"movl		%[__c1E]		,%%ecx\n\t"\
		"movaps		%%xmm2		,%%xmm0\n\t"\
		"movaps		%%xmm4		,%%xmm1\n\t"\
		"mulpd		      (%%ecx)	,%%xmm2\n\t"\
		"mulpd		      (%%ecx)	,%%xmm4\n\t"\
		"mulpd		 0x010(%%ecx)	,%%xmm0\n\t"\
		"mulpd		 0x010(%%ecx)	,%%xmm1\n\t"\
		"subpd		%%xmm0		,%%xmm4\n\t"\
		"addpd		%%xmm1		,%%xmm2\n\t"\
		"movaps		%%xmm4		,%%xmm1\n\t"\
		"movaps		%%xmm2		,%%xmm0\n\t"\
		"unpckhpd	 0x1f0(%%ebx)	,%%xmm1\n\t"\
		"unpcklpd	 0x1f0(%%ebx)	,%%xmm4\n\t"\
		"movaps		%%xmm1		, 0x1f0(%%ebx)\n\t"\
		"unpckhpd	 0x1e0(%%ebx)	,%%xmm0\n\t"\
		"unpcklpd	 0x1e0(%%ebx)	,%%xmm2\n\t"\
		"movaps		%%xmm0		, 0x1e0(%%ebx)\n\t"\
		"movaps		%%xmm4		, 0x1f0(%%eax)\n\t"\
		"movaps		%%xmm2		, 0x1e0(%%eax)\n\t"\
		:					/* outputs: none */\
		: [__add0 ] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__add1] "m" (Xadd1)\
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__r00] "m" (Xr00)\
		 ,[__r08] "m" (Xr08)\
		 ,[__r10] "m" (Xr10)\
		 ,[__r20] "m" (Xr20)\
		 ,[__r28] "m" (Xr28)\
		 ,[__r30] "m" (Xr30)\
		 ,[__c00] "m" (Xc00)\
		 ,[__c01] "m" (Xc01)\
		 ,[__c02] "m" (Xc02)\
		 ,[__c03] "m" (Xc03)\
		 ,[__c04] "m" (Xc04)\
		 ,[__c05] "m" (Xc05)\
		 ,[__c06] "m" (Xc06)\
		 ,[__c07] "m" (Xc07)\
		 ,[__c08] "m" (Xc08)\
		 ,[__c0A] "m" (Xc0A)\
		 ,[__c0C] "m" (Xc0C)\
		 ,[__c0E] "m" (Xc0E)\
		 ,[__c10] "m" (Xc10)\
		 ,[__c12] "m" (Xc12)\
		 ,[__c14] "m" (Xc14)\
		 ,[__c16] "m" (Xc16)\
		 ,[__c18] "m" (Xc18)\
		 ,[__c1A] "m" (Xc1A)\
		 ,[__c1C] "m" (Xc1C)\
		 ,[__c1E] "m" (Xc1E)\
		: "eax","ebx","ecx","edx","edi","esi"		/* Clobbered registers */\
	);\
	}

#endif	/* radix32_wrapper_square_gcc_h_included */

