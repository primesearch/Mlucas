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
#ifndef radix8_dif_dit_pass_gcc_h_included
#define radix8_dif_dit_pass_gcc_h_included

	#define SSE2_RADIX8_DIF_TWIDDLE(Xadd0,Xadd1,Xadd2,Xadd3,Xadd4,Xadd5,Xadd6,Xadd7,Xisrt2,Xc1,Xc2,Xc3,Xc4,Xc5,Xc6,Xc7)\
	{\
	__asm__ volatile (\
		"movl		%[add0]	,%%eax			\n\t"\
		"movl		%[add4]	,%%esi			\n\t"\
		"movl		%[c4]	,%%ecx			\n\t"\
		"movaps		    (%%eax)	,%%xmm0		\n\t"\
		"movaps		0x10(%%eax)	,%%xmm1		\n\t"\
		"movaps		    (%%eax)	,%%xmm6		\n\t"\
		"movaps		0x10(%%eax)	,%%xmm7		\n\t"\
		"movaps		    (%%esi)	,%%xmm2		\n\t"\
		"movaps		0x10(%%esi)	,%%xmm3		\n\t"\
		"movaps		    (%%esi)	,%%xmm4		\n\t"\
		"movaps		0x10(%%esi)	,%%xmm5		\n\t"\
		"mulpd		    (%%ecx)	,%%xmm2		\n\t"\
		"mulpd		    (%%ecx)	,%%xmm3		\n\t"\
		"mulpd		0x10(%%ecx)	,%%xmm4		\n\t"\
		"mulpd		0x10(%%ecx)	,%%xmm5		\n\t"\
		"subpd		%%xmm5		,%%xmm2		\n\t"\
		"addpd		%%xmm4		,%%xmm3		\n\t"\
		"addpd		%%xmm2		,%%xmm0		\n\t"\
		"addpd		%%xmm3		,%%xmm1		\n\t"\
		"subpd		%%xmm2		,%%xmm6		\n\t"\
		"subpd		%%xmm3		,%%xmm7		\n\t"\
		"movaps		%%xmm0		,    (%%eax)	\n\t"\
		"movaps		%%xmm1		,0x10(%%eax)	\n\t"\
		"movaps		%%xmm6		,    (%%esi)	\n\t"\
		"movaps		%%xmm7		,0x10(%%esi)	\n\t"\
		"movl		%[add2]	,%%eax			\n\t"\
		"movl		%[add6]	,%%esi			\n\t"\
		"movl		%[c2]	,%%ecx			\n\t"\
		"movl		%[c6]	,%%edx			\n\t"\
		"movaps		    (%%eax)	,%%xmm0		\n\t"\
		"movaps		0x10(%%eax)	,%%xmm2		\n\t"\
		"movaps		    (%%eax)	,%%xmm1		\n\t"\
		"movaps		0x10(%%eax)	,%%xmm3		\n\t"\
		"mulpd		    (%%ecx)	,%%xmm0		\n\t"\
		"mulpd		0x10(%%ecx)	,%%xmm2		\n\t"\
		"mulpd		0x10(%%ecx)	,%%xmm1		\n\t"\
		"mulpd		    (%%ecx)	,%%xmm3		\n\t"\
		"subpd		%%xmm2		,%%xmm0		\n\t"\
		"addpd		%%xmm3		,%%xmm1		\n\t"\
		"movaps		    (%%esi)	,%%xmm2		\n\t"\
		"movaps		0x10(%%esi)	,%%xmm3		\n\t"\
		"movaps		    (%%esi)	,%%xmm4		\n\t"\
		"movaps		0x10(%%esi)	,%%xmm5		\n\t"\
		"mulpd		    (%%edx)	,%%xmm2		\n\t"\
		"mulpd		0x10(%%edx)	,%%xmm3		\n\t"\
		"mulpd		0x10(%%edx)	,%%xmm4		\n\t"\
		"mulpd		    (%%edx)	,%%xmm5		\n\t"\
		"subpd		%%xmm3		,%%xmm2		\n\t"\
		"addpd		%%xmm5		,%%xmm4		\n\t"\
		"movaps		%%xmm2		,%%xmm3		\n\t"\
		"movaps		%%xmm4		,%%xmm5		\n\t"\
		"addpd		%%xmm0		,%%xmm2		\n\t"\
		"subpd		%%xmm3		,%%xmm0		\n\t"\
		"addpd		%%xmm1		,%%xmm4		\n\t"\
		"subpd		%%xmm5		,%%xmm1		\n\t"\
		"movaps		%%xmm2		,    (%%eax)	\n\t"\
		"movaps		%%xmm4		,0x10(%%eax)	\n\t"\
		"movaps		%%xmm0		,    (%%esi)	\n\t"\
		"movaps		%%xmm1		,0x10(%%esi)	\n\t"\
		"movl		%[add1]	,%%eax			\n\t"\
		"movl		%[add5]	,%%esi			\n\t"\
		"movl		%[c1]	,%%ecx			\n\t"\
		"movl		%[c5]	,%%edx			\n\t"\
		"movaps		    (%%eax)	,%%xmm0		\n\t"\
		"movaps		0x10(%%eax)	,%%xmm2		\n\t"\
		"movaps		    (%%eax)	,%%xmm1		\n\t"\
		"movaps		0x10(%%eax)	,%%xmm3		\n\t"\
		"mulpd		    (%%ecx)	,%%xmm0		\n\t"\
		"mulpd		0x10(%%ecx)	,%%xmm2		\n\t"\
		"mulpd		0x10(%%ecx)	,%%xmm1		\n\t"\
		"mulpd		    (%%ecx)	,%%xmm3		\n\t"\
		"subpd		%%xmm2		,%%xmm0		\n\t"\
		"addpd		%%xmm3		,%%xmm1		\n\t"\
		"movaps		    (%%esi)	,%%xmm2		\n\t"\
		"movaps		0x10(%%esi)	,%%xmm3		\n\t"\
		"movaps		    (%%esi)	,%%xmm4		\n\t"\
		"movaps		0x10(%%esi)	,%%xmm5		\n\t"\
		"mulpd		    (%%edx)	,%%xmm2		\n\t"\
		"mulpd		0x10(%%edx)	,%%xmm3		\n\t"\
		"mulpd		0x10(%%edx)	,%%xmm4		\n\t"\
		"mulpd		    (%%edx)	,%%xmm5		\n\t"\
		"subpd		%%xmm3		,%%xmm2		\n\t"\
		"addpd		%%xmm5		,%%xmm4		\n\t"\
		"movaps		%%xmm2		,%%xmm3		\n\t"\
		"movaps		%%xmm4		,%%xmm5		\n\t"\
		"addpd		%%xmm0		,%%xmm2		\n\t"\
		"subpd		%%xmm3		,%%xmm0		\n\t"\
		"addpd		%%xmm1		,%%xmm4		\n\t"\
		"subpd		%%xmm5		,%%xmm1		\n\t"\
		"movaps		%%xmm2		,    (%%eax)	\n\t"\
		"movaps		%%xmm4		,0x10(%%eax)	\n\t"\
		"movaps		%%xmm0		,    (%%esi)	\n\t"\
		"movaps		%%xmm1		,0x10(%%esi)	\n\t"\
		"movl		%[add3]	,%%eax			\n\t"\
		"movl		%[add7]	,%%esi			\n\t"\
		"movl		%[c3]	,%%ecx			\n\t"\
		"movl		%[c7]	,%%edx			\n\t"\
		"movaps		    (%%eax)	,%%xmm0		\n\t"\
		"movaps		0x10(%%eax)	,%%xmm2		\n\t"\
		"movaps		    (%%eax)	,%%xmm1		\n\t"\
		"movaps		0x10(%%eax)	,%%xmm3		\n\t"\
		"mulpd		    (%%ecx)	,%%xmm0		\n\t"\
		"mulpd		0x10(%%ecx)	,%%xmm2		\n\t"\
		"mulpd		0x10(%%ecx)	,%%xmm1		\n\t"\
		"mulpd		    (%%ecx)	,%%xmm3		\n\t"\
		"subpd		%%xmm2		,%%xmm0		\n\t"\
		"addpd		%%xmm3		,%%xmm1		\n\t"\
		"movaps		    (%%esi)	,%%xmm2		\n\t"\
		"movaps		0x10(%%esi)	,%%xmm3		\n\t"\
		"movaps		    (%%esi)	,%%xmm4		\n\t"\
		"movaps		0x10(%%esi)	,%%xmm5		\n\t"\
		"mulpd		    (%%edx)	,%%xmm2		\n\t"\
		"mulpd		0x10(%%edx)	,%%xmm3		\n\t"\
		"mulpd		0x10(%%edx)	,%%xmm4		\n\t"\
		"mulpd		    (%%edx)	,%%xmm5		\n\t"\
		"subpd		%%xmm3		,%%xmm2		\n\t"\
		"addpd		%%xmm5		,%%xmm4		\n\t"\
		"movaps		%%xmm2		,%%xmm3		\n\t"\
		"movaps		%%xmm4		,%%xmm5		\n\t"\
		"addpd		%%xmm0		,%%xmm2		\n\t"\
		"subpd		%%xmm3		,%%xmm0		\n\t"\
		"addpd		%%xmm1		,%%xmm4		\n\t"\
		"subpd		%%xmm5		,%%xmm1		\n\t"\
		"movaps		%%xmm2		,    (%%eax)	\n\t"\
		"movaps		%%xmm4		,0x10(%%eax)	\n\t"\
		"movaps		%%xmm0		,    (%%esi)	\n\t"\
		"movaps		%%xmm1		,0x10(%%esi)	\n\t"\
		"movl		%[add0]	,%%eax			\n\t"\
		"movl		%[add2]	,%%esi			\n\t"\
		"movaps		    (%%eax)	,%%xmm0		\n\t"\
		"movaps		0x10(%%eax)	,%%xmm1		\n\t"\
		"movaps		%%xmm0		,%%xmm4		\n\t"\
		"movaps		%%xmm1		,%%xmm5		\n\t"\
		"addpd		    (%%esi)	,%%xmm0		\n\t"\
		"subpd		    (%%esi)	,%%xmm4		\n\t"\
		"addpd		0x10(%%esi)	,%%xmm1		\n\t"\
		"subpd		0x10(%%esi)	,%%xmm5		\n\t"\
		"movaps		%%xmm0		,    (%%eax)	\n\t"\
		"movaps		%%xmm1		,0x10(%%eax)	\n\t"\
		"movaps		%%xmm4		,    (%%esi)	\n\t"\
		"movaps		%%xmm5		,0x10(%%esi)	\n\t"\
		"movl		%[add1]	,%%ecx			\n\t"\
		"movl		%[add3]	,%%edx			\n\t"\
		"movaps		    (%%ecx)	,%%xmm2		\n\t"\
		"movaps		0x10(%%ecx)	,%%xmm3		\n\t"\
		"movaps		%%xmm2		,%%xmm6		\n\t"\
		"movaps		%%xmm3		,%%xmm7		\n\t"\
		"addpd		    (%%edx)	,%%xmm2		\n\t"\
		"subpd		    (%%edx)	,%%xmm6		\n\t"\
		"addpd		0x10(%%edx)	,%%xmm3		\n\t"\
		"subpd		0x10(%%edx)	,%%xmm7		\n\t"\
		"subpd		%%xmm2		,%%xmm0		\n\t"\
		"subpd		%%xmm3		,%%xmm1		\n\t"\
		"subpd		%%xmm7		,%%xmm4		\n\t"\
		"subpd		%%xmm6		,%%xmm5		\n\t"\
		"addpd		    (%%eax)	,%%xmm2		\n\t"\
		"addpd		0x10(%%eax)	,%%xmm3		\n\t"\
		"addpd		    (%%esi)	,%%xmm7		\n\t"\
		"addpd		0x10(%%esi)	,%%xmm6		\n\t"\
		"movaps		%%xmm2		,    (%%eax)	\n\t"\
		"movaps		%%xmm3		,0x10(%%eax)	\n\t"\
		"movaps		%%xmm4		,    (%%esi)	\n\t"\
		"movaps		%%xmm6		,0x10(%%esi)	\n\t"\
		"movaps		%%xmm0		,    (%%ecx)	\n\t"\
		"movaps		%%xmm1		,0x10(%%ecx)	\n\t"\
		"movaps		%%xmm7		,    (%%edx)	\n\t"\
		"movaps		%%xmm5		,0x10(%%edx)	\n\t"\
		"movl		%[add4]	,%%eax			\n\t"\
		"movl		%[add6]	,%%esi			\n\t"\
		"movaps		    (%%eax)	,%%xmm0		\n\t"\
		"movaps		0x10(%%eax)	,%%xmm1		\n\t"\
		"movaps		%%xmm0		,%%xmm4		\n\t"\
		"movaps		%%xmm1		,%%xmm5		\n\t"\
		"subpd		0x10(%%esi)	,%%xmm0		\n\t"\
		"addpd		0x10(%%esi)	,%%xmm4		\n\t"\
		"addpd		    (%%esi)	,%%xmm1		\n\t"\
		"subpd		    (%%esi)	,%%xmm5		\n\t"\
		"movaps		%%xmm0		,    (%%eax)	\n\t"\
		"movaps		%%xmm1		,0x10(%%eax)	\n\t"\
		"movaps		%%xmm4		,    (%%esi)	\n\t"\
		"movaps		%%xmm5		,0x10(%%esi)	\n\t"\
		"movl		%[add5]	,%%ecx			\n\t"\
		"movl		%[add7]	,%%edx			\n\t"\
		"movaps		    (%%ecx)	,%%xmm2		\n\t"\
		"movaps		0x10(%%ecx)	,%%xmm3		\n\t"\
		"movaps		%%xmm2		,%%xmm6		\n\t"\
		"movaps		%%xmm3		,%%xmm7		\n\t"\
		"subpd		0x10(%%edx)	,%%xmm2		\n\t"\
		"addpd		0x10(%%edx)	,%%xmm6		\n\t"\
		"addpd		    (%%edx)	,%%xmm3		\n\t"\
		"subpd		    (%%edx)	,%%xmm7		\n\t"\
		"movaps		%%xmm4		,    (%%edx)	\n\t"\
		"movaps		%%xmm5		,0x10(%%edx)	\n\t"\
		"movl		%[isrt2]	,%%eax		\n\t"\
		"movaps		%%xmm2		,%%xmm5		\n\t"\
		"subpd		%%xmm3		,%%xmm2		\n\t"\
		"addpd		%%xmm3		,%%xmm5		\n\t"\
		"mulpd		    (%%eax)	,%%xmm2		\n\t"\
		"mulpd		    (%%eax)	,%%xmm5		\n\t"\
		"movaps		0x10(%%edx)	,%%xmm3		\n\t"\
		"movaps		%%xmm7		,%%xmm4		\n\t"\
		"addpd		%%xmm6		,%%xmm4		\n\t"\
		"subpd		%%xmm6		,%%xmm7		\n\t"\
		"mulpd		(%%eax)		,%%xmm4		\n\t"\
		"mulpd		(%%eax)		,%%xmm7		\n\t"\
		"movaps		(%%edx)		,%%xmm6		\n\t"\
		"movl		%[add4]	,%%eax			\n\t"\
		"subpd		%%xmm2		,%%xmm0		\n\t"\
		"subpd		%%xmm5		,%%xmm1		\n\t"\
		"subpd		%%xmm4		,%%xmm6		\n\t"\
		"subpd		%%xmm7		,%%xmm3		\n\t"\
		"addpd		    (%%eax)	,%%xmm2		\n\t"\
		"addpd		0x10(%%eax)	,%%xmm5		\n\t"\
		"addpd		    (%%esi)	,%%xmm4		\n\t"\
		"addpd		0x10(%%esi)	,%%xmm7		\n\t"\
		"movaps		%%xmm2		,    (%%eax)	\n\t"\
		"movaps		%%xmm5		,0x10(%%eax)	\n\t"\
		"movaps		%%xmm6		,    (%%esi)	\n\t"\
		"movaps		%%xmm3		,0x10(%%esi)	\n\t"\
		"movaps		%%xmm0		,    (%%ecx)	\n\t"\
		"movaps		%%xmm1		,0x10(%%ecx)	\n\t"\
		"movaps		%%xmm4		,    (%%edx)	\n\t"\
		"movaps		%%xmm7		,0x10(%%edx)	\n\t"\
		:					/* outputs: none */\
		: [add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[add1] "m" (Xadd1)\
		 ,[add2] "m" (Xadd2)\
		 ,[add3] "m" (Xadd3)\
		 ,[add4] "m" (Xadd4)\
		 ,[add5] "m" (Xadd5)\
		 ,[add6] "m" (Xadd6)\
		 ,[add7] "m" (Xadd7)\
		 ,[isrt2] "m" (Xisrt2)\
		 ,[c1] "m" (Xc1)\
		 ,[c2] "m" (Xc2)\
		 ,[c3] "m" (Xc3)\
		 ,[c4] "m" (Xc4)\
		 ,[c5] "m" (Xc5)\
		 ,[c6] "m" (Xc6)\
		 ,[c7] "m" (Xc7)\
		: "cc","memory","eax","esi","ecx","edx","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"	/* Clobbered registers */\
	);\
	}

	#define SSE2_RADIX8_DIT_TWIDDLE(Xadd0,Xadd1,Xadd2,Xadd3,Xadd4,Xadd5,Xadd6,Xadd7,Xisrt2,Xc1,Xc2,Xc3,Xc4,Xc5,Xc6,Xc7)\
	{\
	__asm__ volatile (\
		"movl		%[add4]		,%%eax		\n\t"\
		"movl		%[add5]		,%%esi		\n\t"\
		"movaps		(%%eax)		,%%xmm0		\n\t"\
		"movaps		0x10(%%eax)	,%%xmm1		\n\t"\
		"movaps		%%xmm0		,%%xmm2		\n\t"\
		"movaps		%%xmm1		,%%xmm3		\n\t"\
		"addpd		(%%esi)		,%%xmm2		\n\t"\
		"addpd		0x10(%%esi)	,%%xmm3		\n\t"\
		"subpd		(%%esi)		,%%xmm0		\n\t"\
		"subpd		0x10(%%esi)	,%%xmm1		\n\t"\
		"movl		%[add6]		,%%ecx		\n\t"\
		"movl		%[add7]		,%%edx		\n\t"\
		"movaps		(%%ecx)		,%%xmm4		\n\t"\
		"movaps		0x10(%%ecx)	,%%xmm5		\n\t"\
		"movaps		%%xmm4		,%%xmm6		\n\t"\
		"movaps		%%xmm5		,%%xmm7		\n\t"\
		"addpd		(%%edx)		,%%xmm6		\n\t"\
		"addpd		0x10(%%edx)	,%%xmm7		\n\t"\
		"subpd		(%%edx)		,%%xmm4		\n\t"\
		"subpd		0x10(%%edx)	,%%xmm5		\n\t"\
		"movaps		%%xmm6		,    (%%ecx)\n\t"\
		"movaps		%%xmm7		,0x10(%%ecx)\n\t"\
		"movaps		%%xmm4		,    (%%edx)\n\t"\
		"movaps		%%xmm5		,0x10(%%edx)\n\t"\
		"addpd		%%xmm2		,%%xmm6		\n\t"\
		"addpd		%%xmm3		,%%xmm7		\n\t"\
		"subpd		    (%%ecx)	,%%xmm2		\n\t"\
		"subpd		0x10(%%ecx)	,%%xmm3		\n\t"\
		"movaps		%%xmm2		,    (%%ecx)\n\t"\
		"movaps		%%xmm3		,0x10(%%ecx)\n\t"\
		"movaps		%%xmm4		,%%xmm2		\n\t"\
		"movaps		%%xmm5		,%%xmm3		\n\t"\
		"addpd		%%xmm0		,%%xmm5		\n\t"\
		"subpd		%%xmm3		,%%xmm0		\n\t"\
		"addpd		%%xmm1		,%%xmm4		\n\t"\
		"subpd		%%xmm2		,%%xmm1		\n\t"\
		"movaps		%%xmm5		,%%xmm2		\n\t"\
		"movaps		%%xmm1		,%%xmm3		\n\t"\
		"addpd		%%xmm1		,%%xmm5		\n\t"\
		"movl		%[isrt2]	,%%ecx		\n\t"\
		"movaps		(%%ecx)		,%%xmm1		\n\t"\
		"subpd		%%xmm3		,%%xmm2		\n\t"\
		"mulpd		%%xmm1		,%%xmm5		\n\t"\
		"mulpd		%%xmm1		,%%xmm2		\n\t"\
		"movaps		%%xmm0		,%%xmm3		\n\t"\
		"movaps		%%xmm5		,    (%%esi)\n\t"\
		"movaps		%%xmm4		,%%xmm5		\n\t"\
		"addpd		%%xmm4		,%%xmm0		\n\t"\
		"movaps		%%xmm2		,0x10(%%esi)\n\t"\
		"subpd		%%xmm5		,%%xmm3		\n\t"\
		"mulpd		%%xmm1		,%%xmm0		\n\t"\
		"mulpd		%%xmm1		,%%xmm3		\n\t"\
		"movaps		%%xmm0		,    (%%edx)\n\t"\
		"movaps		%%xmm3		,0x10(%%edx)\n\t"\
		"movl		%[add0]		,%%eax		\n\t"\
		"movl		%[add1]		,%%esi		\n\t"\
		"movaps		(%%eax)		,%%xmm0		\n\t"\
		"movaps		0x10(%%eax)	,%%xmm1		\n\t"\
		"movaps		%%xmm0		,%%xmm2		\n\t"\
		"movaps		%%xmm1		,%%xmm3		\n\t"\
		"addpd		(%%esi)		,%%xmm2		\n\t"\
		"addpd		0x10(%%esi)	,%%xmm3		\n\t"\
		"subpd		(%%esi)		,%%xmm0		\n\t"\
		"subpd		0x10(%%esi)	,%%xmm1		\n\t"\
		"movaps		%%xmm6		,    (%%eax)\n\t"\
		"movaps		%%xmm7		,0x10(%%eax)\n\t"\
		"addpd		%%xmm6		,%%xmm6		\n\t"\
		"addpd		%%xmm7		,%%xmm7		\n\t"\
		"movl		%[add4]		,%%ecx		\n\t"\
		"movaps		%%xmm6		,    (%%ecx)\n\t"\
		"movaps		%%xmm7		,0x10(%%ecx)\n\t"\
		"movl		%[add2]		,%%ecx		\n\t"\
		"movl		%[add3]		,%%edx		\n\t"\
		"movaps		(%%ecx)		,%%xmm4		\n\t"\
		"movaps		0x10(%%ecx)	,%%xmm5		\n\t"\
		"movaps		%%xmm4		,%%xmm6		\n\t"\
		"movaps		%%xmm5		,%%xmm7		\n\t"\
		"addpd		(%%edx)		,%%xmm6		\n\t"\
		"addpd		0x10(%%edx)	,%%xmm7		\n\t"\
		"subpd		(%%edx)		,%%xmm4		\n\t"\
		"subpd		0x10(%%edx)	,%%xmm5		\n\t"\
		"movaps		%%xmm6		,    (%%ecx)\n\t"\
		"movaps		%%xmm7		,0x10(%%ecx)\n\t"\
		"addpd		%%xmm2		,%%xmm6		\n\t"\
		"addpd		%%xmm3		,%%xmm7		\n\t"\
		"subpd		    (%%ecx)	,%%xmm2		\n\t"\
		"subpd		0x10(%%ecx)	,%%xmm3		\n\t"\
		"addpd		    (%%eax)	,%%xmm6		\n\t"\
		"addpd		0x10(%%eax)	,%%xmm7		\n\t"\
		"movaps		%%xmm6		,    (%%eax)\n\t"\
		"movaps		%%xmm7		,0x10(%%eax)\n\t"\
		"movl		%[add4]		,%%esi		\n\t"\
		"subpd		    (%%esi)	,%%xmm6		\n\t"\
		"subpd		0x10(%%esi)	,%%xmm7		\n\t"\
		"movaps		%%xmm6		,    (%%esi)\n\t"\
		"movaps		%%xmm7		,0x10(%%esi)\n\t"\
		"movaps		%%xmm6		,    (%%edx)\n\t"\
		"movaps		%%xmm7		,0x10(%%edx)\n\t"\
		"movaps		%%xmm4		,%%xmm6		\n\t"\
		"movaps		%%xmm5		,%%xmm7		\n\t"\
		"addpd		%%xmm0		,%%xmm5		\n\t"\
		"subpd		%%xmm7		,%%xmm0		\n\t"\
		"addpd		%%xmm1		,%%xmm4		\n\t"\
		"subpd		%%xmm6		,%%xmm1		\n\t"\
		"movl		%[c4]		,%%ecx		\n\t"\
		"movaps		    (%%esi)	,%%xmm6		\n\t"\
		"movaps		0x10(%%esi)	,%%xmm7		\n\t"\
		"mulpd		0x10(%%ecx)	,%%xmm6		\n\t"\
		"mulpd		0x10(%%ecx)	,%%xmm7		\n\t"\
		"movaps		%%xmm6		,    (%%esi)\n\t"\
		"movaps		%%xmm7		,0x10(%%esi)\n\t"\
		"movaps		    (%%edx)	,%%xmm6		\n\t"\
		"movaps		0x10(%%edx)	,%%xmm7		\n\t"\
		"mulpd		    (%%ecx)	,%%xmm6		\n\t"\
		"mulpd		    (%%ecx)	,%%xmm7		\n\t"\
		"addpd		0x10(%%esi)	,%%xmm6		\n\t"\
		"subpd		    (%%esi)	,%%xmm7		\n\t"\
		"movaps		%%xmm6		,    (%%esi)\n\t"\
		"movaps		%%xmm7		,0x10(%%esi)\n\t"\
		"movl		%[add1]		,%%eax		\n\t"\
		"movl		%[add5]		,%%esi		\n\t"\
		"movl		%[c1]		,%%ecx		\n\t"\
		"movl		%[c5]		,%%edx		\n\t"\
		"movaps		%%xmm5		,%%xmm6		\n\t"\
		"movaps		%%xmm1		,%%xmm7		\n\t"\
		"addpd		    (%%esi)	,%%xmm5		\n\t"\
		"subpd		0x10(%%esi)	,%%xmm1		\n\t"\
		"subpd		    (%%esi)	,%%xmm6		\n\t"\
		"addpd		0x10(%%esi)	,%%xmm7		\n\t"\
		"movaps		%%xmm6		,    (%%esi)\n\t"\
		"movaps		%%xmm7		,0x10(%%esi)\n\t"\
		"movaps		%%xmm5		,%%xmm6		\n\t"\
		"movaps		%%xmm1		,%%xmm7		\n\t"\
		"mulpd		    (%%ecx)	,%%xmm5		\n\t"\
		"mulpd		    (%%ecx)	,%%xmm1		\n\t"\
		"mulpd		0x10(%%ecx)	,%%xmm6		\n\t"\
		"mulpd		0x10(%%ecx)	,%%xmm7		\n\t"\
		"subpd		%%xmm6		,%%xmm1		\n\t"\
		"addpd		%%xmm7		,%%xmm5		\n\t"\
		"movaps		%%xmm1		,0x10(%%eax)\n\t"\
		"movaps		%%xmm5		,    (%%eax)\n\t"\
		"movaps		    (%%esi)	,%%xmm5		\n\t"\
		"movaps		0x10(%%esi)	,%%xmm1		\n\t"\
		"movaps		%%xmm5		,%%xmm6		\n\t"\
		"movaps		%%xmm1		,%%xmm7		\n\t"\
		"mulpd		    (%%edx)	,%%xmm5		\n\t"\
		"mulpd		    (%%edx)	,%%xmm1		\n\t"\
		"mulpd		0x10(%%edx)	,%%xmm6		\n\t"\
		"mulpd		0x10(%%edx)	,%%xmm7		\n\t"\
		"subpd		%%xmm6		,%%xmm1		\n\t"\
		"addpd		%%xmm7		,%%xmm5		\n\t"\
		"movaps		%%xmm1		,0x10(%%esi)\n\t"\
		"movaps		%%xmm5		,    (%%esi)\n\t"\
		"movl		%[add2]		,%%eax		\n\t"\
		"movl		%[add6]		,%%esi		\n\t"\
		"movl		%[c2]		,%%ecx		\n\t"\
		"movl		%[c6]		,%%edx		\n\t"\
		"movaps		%%xmm2		,%%xmm6		\n\t"\
		"movaps		%%xmm3		,%%xmm7		\n\t"\
		"addpd		0x10(%%esi)	,%%xmm2		\n\t"\
		"subpd		    (%%esi)	,%%xmm3		\n\t"\
		"subpd		0x10(%%esi)	,%%xmm6		\n\t"\
		"addpd		    (%%esi)	,%%xmm7		\n\t"\
		"movaps		%%xmm2		,%%xmm1		\n\t"\
		"movaps		%%xmm3		,%%xmm5		\n\t"\
		"mulpd		    (%%ecx)	,%%xmm2		\n\t"\
		"mulpd		    (%%ecx)	,%%xmm3		\n\t"\
		"mulpd		0x10(%%ecx)	,%%xmm1		\n\t"\
		"mulpd		0x10(%%ecx)	,%%xmm5		\n\t"\
		"subpd		%%xmm1		,%%xmm3		\n\t"\
		"addpd		%%xmm5		,%%xmm2		\n\t"\
		"movaps		%%xmm3		,0x10(%%eax)\n\t"\
		"movaps		%%xmm2		,    (%%eax)\n\t"\
		"movaps		%%xmm6		,%%xmm1		\n\t"\
		"movaps		%%xmm7		,%%xmm5		\n\t"\
		"mulpd		    (%%edx)	,%%xmm6		\n\t"\
		"mulpd		    (%%edx)	,%%xmm7		\n\t"\
		"mulpd		0x10(%%edx)	,%%xmm1		\n\t"\
		"mulpd		0x10(%%edx)	,%%xmm5		\n\t"\
		"subpd		%%xmm1		,%%xmm7		\n\t"\
		"addpd		%%xmm5		,%%xmm6		\n\t"\
		"movaps		%%xmm7		,0x10(%%esi)\n\t"\
		"movaps		%%xmm6		,    (%%esi)\n\t"\
		"movl		%[add3]		,%%eax		\n\t"\
		"movl		%[add7]		,%%esi		\n\t"\
		"movl		%[c3]		,%%ecx		\n\t"\
		"movl		%[c7]		,%%edx		\n\t"\
		"movaps		%%xmm0		,%%xmm6		\n\t"\
		"movaps		%%xmm4		,%%xmm7		\n\t"\
		"subpd		0x10(%%esi)	,%%xmm0		\n\t"\
		"subpd		    (%%esi)	,%%xmm4		\n\t"\
		"addpd		0x10(%%esi)	,%%xmm6		\n\t"\
		"addpd		    (%%esi)	,%%xmm7		\n\t"\
		"movaps		%%xmm0		,%%xmm1		\n\t"\
		"movaps		%%xmm4		,%%xmm5		\n\t"\
		"mulpd		    (%%ecx)	,%%xmm0		\n\t"\
		"mulpd		    (%%ecx)	,%%xmm4		\n\t"\
		"mulpd		0x10(%%ecx)	,%%xmm1		\n\t"\
		"mulpd		0x10(%%ecx)	,%%xmm5		\n\t"\
		"subpd		%%xmm1		,%%xmm4		\n\t"\
		"addpd		%%xmm5		,%%xmm0		\n\t"\
		"movaps		%%xmm4		,0x10(%%eax)\n\t"\
		"movaps		%%xmm0		,    (%%eax)\n\t"\
		"movaps		%%xmm6		,%%xmm1		\n\t"\
		"movaps		%%xmm7		,%%xmm5		\n\t"\
		"mulpd		    (%%edx)	,%%xmm6		\n\t"\
		"mulpd		    (%%edx)	,%%xmm7		\n\t"\
		"mulpd		0x10(%%edx)	,%%xmm1		\n\t"\
		"mulpd		0x10(%%edx)	,%%xmm5		\n\t"\
		"subpd		%%xmm1		,%%xmm7		\n\t"\
		"addpd		%%xmm5		,%%xmm6		\n\t"\
		"movaps		%%xmm7		,0x10(%%esi)\n\t"\
		"movaps		%%xmm6		,    (%%esi)\n\t"\
		:					/* outputs: none */\
		: [add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[add1] "m" (Xadd1)\
		 ,[add2] "m" (Xadd2)\
		 ,[add3] "m" (Xadd3)\
		 ,[add4] "m" (Xadd4)\
		 ,[add5] "m" (Xadd5)\
		 ,[add6] "m" (Xadd6)\
		 ,[add7] "m" (Xadd7)\
		 ,[isrt2] "m" (Xisrt2)\
		 ,[c1] "m" (Xc1)\
		 ,[c2] "m" (Xc2)\
		 ,[c3] "m" (Xc3)\
		 ,[c4] "m" (Xc4)\
		 ,[c5] "m" (Xc5)\
		 ,[c6] "m" (Xc6)\
		 ,[c7] "m" (Xc7)\
		: "cc","memory","eax","esi","ecx","edx","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"	/* Clobbered registers */\
	);\
	}

#endif	/* radix8_dif_dit_pass_gcc_h_included */

