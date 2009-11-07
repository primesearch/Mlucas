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
#ifndef radix8_dif_dit_pass_gcc_h_included
#define radix8_dif_dit_pass_gcc_h_included

	#define SSE2_RADIX8_DIF_TWIDDLE(Xadd0,Xadd1,Xadd2,Xadd3,Xadd4,Xadd5,Xadd6,Xadd7,Xisrt2,Xc1,Xc2,Xc3,Xc4,Xc5,Xc6,Xc7)\
	{\
	__asm__ volatile (\
		"movq		%[add0]	,%%rax			\n\t"\
		"movq		%[add4]	,%%r15			\n\t"\
		"movq		%[c4]	,%%rcx			\n\t"\
		"movaps		    (%%rax)	,%%xmm0		\n\t"\
		"movaps		0x10(%%rax)	,%%xmm1		\n\t"\
		"movaps		    (%%rax)	,%%xmm6		\n\t"\
		"movaps		0x10(%%rax)	,%%xmm7		\n\t"\
		"movaps		    (%%r15)	,%%xmm2		\n\t"\
		"movaps		0x10(%%r15)	,%%xmm3		\n\t"\
		"movaps		    (%%r15)	,%%xmm4		\n\t"\
		"movaps		0x10(%%r15)	,%%xmm5		\n\t"\
		"mulpd		    (%%rcx)	,%%xmm2		\n\t"\
		"mulpd		    (%%rcx)	,%%xmm3		\n\t"\
		"mulpd		0x10(%%rcx)	,%%xmm4		\n\t"\
		"mulpd		0x10(%%rcx)	,%%xmm5		\n\t"\
		"subpd		%%xmm5		,%%xmm2		\n\t"\
		"addpd		%%xmm4		,%%xmm3		\n\t"\
		"addpd		%%xmm2		,%%xmm0		\n\t"\
		"addpd		%%xmm3		,%%xmm1		\n\t"\
		"subpd		%%xmm2		,%%xmm6		\n\t"\
		"subpd		%%xmm3		,%%xmm7		\n\t"\
		"movaps		%%xmm0		,    (%%rax)	\n\t"\
		"movaps		%%xmm1		,0x10(%%rax)	\n\t"\
		"movaps		%%xmm6		,    (%%r15)	\n\t"\
		"movaps		%%xmm7		,0x10(%%r15)	\n\t"\
		"movq		%[add2]	,%%rax			\n\t"\
		"movq		%[add6]	,%%r15			\n\t"\
		"movq		%[c2]	,%%rcx			\n\t"\
		"movq		%[c6]	,%%rdx			\n\t"\
		"movaps		    (%%rax)	,%%xmm0		\n\t"\
		"movaps		0x10(%%rax)	,%%xmm2		\n\t"\
		"movaps		    (%%rax)	,%%xmm1		\n\t"\
		"movaps		0x10(%%rax)	,%%xmm3		\n\t"\
		"mulpd		    (%%rcx)	,%%xmm0		\n\t"\
		"mulpd		0x10(%%rcx)	,%%xmm2		\n\t"\
		"mulpd		0x10(%%rcx)	,%%xmm1		\n\t"\
		"mulpd		    (%%rcx)	,%%xmm3		\n\t"\
		"subpd		%%xmm2		,%%xmm0		\n\t"\
		"addpd		%%xmm3		,%%xmm1		\n\t"\
		"movaps		    (%%r15)	,%%xmm2		\n\t"\
		"movaps		0x10(%%r15)	,%%xmm3		\n\t"\
		"movaps		    (%%r15)	,%%xmm4		\n\t"\
		"movaps		0x10(%%r15)	,%%xmm5		\n\t"\
		"mulpd		    (%%rdx)	,%%xmm2		\n\t"\
		"mulpd		0x10(%%rdx)	,%%xmm3		\n\t"\
		"mulpd		0x10(%%rdx)	,%%xmm4		\n\t"\
		"mulpd		    (%%rdx)	,%%xmm5		\n\t"\
		"subpd		%%xmm3		,%%xmm2		\n\t"\
		"addpd		%%xmm5		,%%xmm4		\n\t"\
		"movaps		%%xmm2		,%%xmm3		\n\t"\
		"movaps		%%xmm4		,%%xmm5		\n\t"\
		"addpd		%%xmm0		,%%xmm2		\n\t"\
		"subpd		%%xmm3		,%%xmm0		\n\t"\
		"addpd		%%xmm1		,%%xmm4		\n\t"\
		"subpd		%%xmm5		,%%xmm1		\n\t"\
		"movaps		%%xmm2		,    (%%rax)	\n\t"\
		"movaps		%%xmm4		,0x10(%%rax)	\n\t"\
		"movaps		%%xmm0		,    (%%r15)	\n\t"\
		"movaps		%%xmm1		,0x10(%%r15)	\n\t"\
		"movq		%[add1]	,%%rax			\n\t"\
		"movq		%[add5]	,%%r15			\n\t"\
		"movq		%[c1]	,%%rcx			\n\t"\
		"movq		%[c5]	,%%rdx			\n\t"\
		"movaps		    (%%rax)	,%%xmm0		\n\t"\
		"movaps		0x10(%%rax)	,%%xmm2		\n\t"\
		"movaps		    (%%rax)	,%%xmm1		\n\t"\
		"movaps		0x10(%%rax)	,%%xmm3		\n\t"\
		"mulpd		    (%%rcx)	,%%xmm0		\n\t"\
		"mulpd		0x10(%%rcx)	,%%xmm2		\n\t"\
		"mulpd		0x10(%%rcx)	,%%xmm1		\n\t"\
		"mulpd		    (%%rcx)	,%%xmm3		\n\t"\
		"subpd		%%xmm2		,%%xmm0		\n\t"\
		"addpd		%%xmm3		,%%xmm1		\n\t"\
		"movaps		    (%%r15)	,%%xmm2		\n\t"\
		"movaps		0x10(%%r15)	,%%xmm3		\n\t"\
		"movaps		    (%%r15)	,%%xmm4		\n\t"\
		"movaps		0x10(%%r15)	,%%xmm5		\n\t"\
		"mulpd		    (%%rdx)	,%%xmm2		\n\t"\
		"mulpd		0x10(%%rdx)	,%%xmm3		\n\t"\
		"mulpd		0x10(%%rdx)	,%%xmm4		\n\t"\
		"mulpd		    (%%rdx)	,%%xmm5		\n\t"\
		"subpd		%%xmm3		,%%xmm2		\n\t"\
		"addpd		%%xmm5		,%%xmm4		\n\t"\
		"movaps		%%xmm2		,%%xmm3		\n\t"\
		"movaps		%%xmm4		,%%xmm5		\n\t"\
		"addpd		%%xmm0		,%%xmm2		\n\t"\
		"subpd		%%xmm3		,%%xmm0		\n\t"\
		"addpd		%%xmm1		,%%xmm4		\n\t"\
		"subpd		%%xmm5		,%%xmm1		\n\t"\
		"movaps		%%xmm2		,    (%%rax)	\n\t"\
		"movaps		%%xmm4		,0x10(%%rax)	\n\t"\
		"movaps		%%xmm0		,    (%%r15)	\n\t"\
		"movaps		%%xmm1		,0x10(%%r15)	\n\t"\
		"movq		%[add3]	,%%rax			\n\t"\
		"movq		%[add7]	,%%r15			\n\t"\
		"movq		%[c3]	,%%rcx			\n\t"\
		"movq		%[c7]	,%%rdx			\n\t"\
		"movaps		    (%%rax)	,%%xmm0		\n\t"\
		"movaps		0x10(%%rax)	,%%xmm2		\n\t"\
		"movaps		    (%%rax)	,%%xmm1		\n\t"\
		"movaps		0x10(%%rax)	,%%xmm3		\n\t"\
		"mulpd		    (%%rcx)	,%%xmm0		\n\t"\
		"mulpd		0x10(%%rcx)	,%%xmm2		\n\t"\
		"mulpd		0x10(%%rcx)	,%%xmm1		\n\t"\
		"mulpd		    (%%rcx)	,%%xmm3		\n\t"\
		"subpd		%%xmm2		,%%xmm0		\n\t"\
		"addpd		%%xmm3		,%%xmm1		\n\t"\
		"movaps		    (%%r15)	,%%xmm2		\n\t"\
		"movaps		0x10(%%r15)	,%%xmm3		\n\t"\
		"movaps		    (%%r15)	,%%xmm4		\n\t"\
		"movaps		0x10(%%r15)	,%%xmm5		\n\t"\
		"mulpd		    (%%rdx)	,%%xmm2		\n\t"\
		"mulpd		0x10(%%rdx)	,%%xmm3		\n\t"\
		"mulpd		0x10(%%rdx)	,%%xmm4		\n\t"\
		"mulpd		    (%%rdx)	,%%xmm5		\n\t"\
		"subpd		%%xmm3		,%%xmm2		\n\t"\
		"addpd		%%xmm5		,%%xmm4		\n\t"\
		"movaps		%%xmm2		,%%xmm3		\n\t"\
		"movaps		%%xmm4		,%%xmm5		\n\t"\
		"addpd		%%xmm0		,%%xmm2		\n\t"\
		"subpd		%%xmm3		,%%xmm0		\n\t"\
		"addpd		%%xmm1		,%%xmm4		\n\t"\
		"subpd		%%xmm5		,%%xmm1		\n\t"\
		"movaps		%%xmm2		,    (%%rax)	\n\t"\
		"movaps		%%xmm4		,0x10(%%rax)	\n\t"\
		"movaps		%%xmm0		,    (%%r15)	\n\t"\
		"movaps		%%xmm1		,0x10(%%r15)	\n\t"\
		"movq		%[add0]	,%%rax			\n\t"\
		"movq		%[add2]	,%%r15			\n\t"\
		"movaps		    (%%rax)	,%%xmm0		\n\t"\
		"movaps		0x10(%%rax)	,%%xmm1		\n\t"\
		"movaps		%%xmm0		,%%xmm4		\n\t"\
		"movaps		%%xmm1		,%%xmm5		\n\t"\
		"addpd		    (%%r15)	,%%xmm0		\n\t"\
		"subpd		    (%%r15)	,%%xmm4		\n\t"\
		"addpd		0x10(%%r15)	,%%xmm1		\n\t"\
		"subpd		0x10(%%r15)	,%%xmm5		\n\t"\
		"movaps		%%xmm0		,    (%%rax)	\n\t"\
		"movaps		%%xmm1		,0x10(%%rax)	\n\t"\
		"movaps		%%xmm4		,    (%%r15)	\n\t"\
		"movaps		%%xmm5		,0x10(%%r15)	\n\t"\
		"movq		%[add1]	,%%rcx			\n\t"\
		"movq		%[add3]	,%%rdx			\n\t"\
		"movaps		    (%%rcx)	,%%xmm2		\n\t"\
		"movaps		0x10(%%rcx)	,%%xmm3		\n\t"\
		"movaps		%%xmm2		,%%xmm6		\n\t"\
		"movaps		%%xmm3		,%%xmm7		\n\t"\
		"addpd		    (%%rdx)	,%%xmm2		\n\t"\
		"subpd		    (%%rdx)	,%%xmm6		\n\t"\
		"addpd		0x10(%%rdx)	,%%xmm3		\n\t"\
		"subpd		0x10(%%rdx)	,%%xmm7		\n\t"\
		"subpd		%%xmm2		,%%xmm0		\n\t"\
		"subpd		%%xmm3		,%%xmm1		\n\t"\
		"subpd		%%xmm7		,%%xmm4		\n\t"\
		"subpd		%%xmm6		,%%xmm5		\n\t"\
		"addpd		    (%%rax)	,%%xmm2		\n\t"\
		"addpd		0x10(%%rax)	,%%xmm3		\n\t"\
		"addpd		    (%%r15)	,%%xmm7		\n\t"\
		"addpd		0x10(%%r15)	,%%xmm6		\n\t"\
		"movaps		%%xmm2		,    (%%rax)	\n\t"\
		"movaps		%%xmm3		,0x10(%%rax)	\n\t"\
		"movaps		%%xmm4		,    (%%r15)	\n\t"\
		"movaps		%%xmm6		,0x10(%%r15)	\n\t"\
		"movaps		%%xmm0		,    (%%rcx)	\n\t"\
		"movaps		%%xmm1		,0x10(%%rcx)	\n\t"\
		"movaps		%%xmm7		,    (%%rdx)	\n\t"\
		"movaps		%%xmm5		,0x10(%%rdx)	\n\t"\
		"movq		%[add4]	,%%rax			\n\t"\
		"movq		%[add6]	,%%r15			\n\t"\
		"movaps		    (%%rax)	,%%xmm0		\n\t"\
		"movaps		0x10(%%rax)	,%%xmm1		\n\t"\
		"movaps		%%xmm0		,%%xmm4		\n\t"\
		"movaps		%%xmm1		,%%xmm5		\n\t"\
		"subpd		0x10(%%r15)	,%%xmm0		\n\t"\
		"addpd		0x10(%%r15)	,%%xmm4		\n\t"\
		"addpd		    (%%r15)	,%%xmm1		\n\t"\
		"subpd		    (%%r15)	,%%xmm5		\n\t"\
		"movaps		%%xmm0		,    (%%rax)	\n\t"\
		"movaps		%%xmm1		,0x10(%%rax)	\n\t"\
		"movaps		%%xmm4		,    (%%r15)	\n\t"\
		"movaps		%%xmm5		,0x10(%%r15)	\n\t"\
		"movq		%[add5]	,%%rcx			\n\t"\
		"movq		%[add7]	,%%rdx			\n\t"\
		"movaps		    (%%rcx)	,%%xmm2		\n\t"\
		"movaps		0x10(%%rcx)	,%%xmm3		\n\t"\
		"movaps		%%xmm2		,%%xmm6		\n\t"\
		"movaps		%%xmm3		,%%xmm7		\n\t"\
		"subpd		0x10(%%rdx)	,%%xmm2		\n\t"\
		"addpd		0x10(%%rdx)	,%%xmm6		\n\t"\
		"addpd		    (%%rdx)	,%%xmm3		\n\t"\
		"subpd		    (%%rdx)	,%%xmm7		\n\t"\
		"movaps		%%xmm4		,    (%%rdx)	\n\t"\
		"movaps		%%xmm5		,0x10(%%rdx)	\n\t"\
		"movq		%[isrt2]	,%%rax		\n\t"\
		"movaps		%%xmm2		,%%xmm5		\n\t"\
		"subpd		%%xmm3		,%%xmm2		\n\t"\
		"addpd		%%xmm3		,%%xmm5		\n\t"\
		"mulpd		    (%%rax)	,%%xmm2		\n\t"\
		"mulpd		    (%%rax)	,%%xmm5		\n\t"\
		"movaps		0x10(%%rdx)	,%%xmm3		\n\t"\
		"movaps		%%xmm7		,%%xmm4		\n\t"\
		"addpd		%%xmm6		,%%xmm4		\n\t"\
		"subpd		%%xmm6		,%%xmm7		\n\t"\
		"mulpd		(%%rax)		,%%xmm4		\n\t"\
		"mulpd		(%%rax)		,%%xmm7		\n\t"\
		"movaps		(%%rdx)		,%%xmm6		\n\t"\
		"movq		%[add4]	,%%rax			\n\t"\
		"subpd		%%xmm2		,%%xmm0		\n\t"\
		"subpd		%%xmm5		,%%xmm1		\n\t"\
		"subpd		%%xmm4		,%%xmm6		\n\t"\
		"subpd		%%xmm7		,%%xmm3		\n\t"\
		"addpd		    (%%rax)	,%%xmm2		\n\t"\
		"addpd		0x10(%%rax)	,%%xmm5		\n\t"\
		"addpd		    (%%r15)	,%%xmm4		\n\t"\
		"addpd		0x10(%%r15)	,%%xmm7		\n\t"\
		"movaps		%%xmm2		,    (%%rax)	\n\t"\
		"movaps		%%xmm5		,0x10(%%rax)	\n\t"\
		"movaps		%%xmm6		,    (%%r15)	\n\t"\
		"movaps		%%xmm3		,0x10(%%r15)	\n\t"\
		"movaps		%%xmm0		,    (%%rcx)	\n\t"\
		"movaps		%%xmm1		,0x10(%%rcx)	\n\t"\
		"movaps		%%xmm4		,    (%%rdx)	\n\t"\
		"movaps		%%xmm7		,0x10(%%rdx)	\n\t"\
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
		: "rax","r15","rcx","rdx","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
	);\
	}

	#define SSE2_RADIX8_DIT_TWIDDLE(Xadd0,Xadd1,Xadd2,Xadd3,Xadd4,Xadd5,Xadd6,Xadd7,Xisrt2,Xc1,Xc2,Xc3,Xc4,Xc5,Xc6,Xc7)\
	{\
	__asm__ volatile (\
		"movq		%[add4]		,%%rax		\n\t"\
		"movq		%[add5]		,%%r15		\n\t"\
		"movaps		(%%rax)		,%%xmm0		\n\t"\
		"movaps		0x10(%%rax)	,%%xmm1		\n\t"\
		"movaps		%%xmm0		,%%xmm2		\n\t"\
		"movaps		%%xmm1		,%%xmm3		\n\t"\
		"addpd		(%%r15)		,%%xmm2		\n\t"\
		"addpd		0x10(%%r15)	,%%xmm3		\n\t"\
		"subpd		(%%r15)		,%%xmm0		\n\t"\
		"subpd		0x10(%%r15)	,%%xmm1		\n\t"\
		"movq		%[add6]		,%%rcx		\n\t"\
		"movq		%[add7]		,%%rdx		\n\t"\
		"movaps		(%%rcx)		,%%xmm4		\n\t"\
		"movaps		0x10(%%rcx)	,%%xmm5		\n\t"\
		"movaps		%%xmm4		,%%xmm6		\n\t"\
		"movaps		%%xmm5		,%%xmm7		\n\t"\
		"addpd		(%%rdx)		,%%xmm6		\n\t"\
		"addpd		0x10(%%rdx)	,%%xmm7		\n\t"\
		"subpd		(%%rdx)		,%%xmm4		\n\t"\
		"subpd		0x10(%%rdx)	,%%xmm5		\n\t"\
		"movaps		%%xmm6		,    (%%rcx)\n\t"\
		"movaps		%%xmm7		,0x10(%%rcx)\n\t"\
		"movaps		%%xmm4		,    (%%rdx)\n\t"\
		"movaps		%%xmm5		,0x10(%%rdx)\n\t"\
		"addpd		%%xmm2		,%%xmm6		\n\t"\
		"addpd		%%xmm3		,%%xmm7		\n\t"\
		"subpd		    (%%rcx)	,%%xmm2		\n\t"\
		"subpd		0x10(%%rcx)	,%%xmm3		\n\t"\
		"movaps		%%xmm2		,    (%%rcx)\n\t"\
		"movaps		%%xmm3		,0x10(%%rcx)\n\t"\
		"movaps		%%xmm4		,%%xmm2		\n\t"\
		"movaps		%%xmm5		,%%xmm3		\n\t"\
		"addpd		%%xmm0		,%%xmm5		\n\t"\
		"subpd		%%xmm3		,%%xmm0		\n\t"\
		"addpd		%%xmm1		,%%xmm4		\n\t"\
		"subpd		%%xmm2		,%%xmm1		\n\t"\
		"movaps		%%xmm5		,%%xmm2		\n\t"\
		"movaps		%%xmm1		,%%xmm3		\n\t"\
		"addpd		%%xmm1		,%%xmm5		\n\t"\
		"movq		%[isrt2]	,%%rcx		\n\t"\
		"movaps		(%%rcx)		,%%xmm1		\n\t"\
		"subpd		%%xmm3		,%%xmm2		\n\t"\
		"mulpd		%%xmm1		,%%xmm5		\n\t"\
		"mulpd		%%xmm1		,%%xmm2		\n\t"\
		"movaps		%%xmm0		,%%xmm3		\n\t"\
		"movaps		%%xmm5		,    (%%r15)\n\t"\
		"movaps		%%xmm4		,%%xmm5		\n\t"\
		"addpd		%%xmm4		,%%xmm0		\n\t"\
		"movaps		%%xmm2		,0x10(%%r15)\n\t"\
		"subpd		%%xmm5		,%%xmm3		\n\t"\
		"mulpd		%%xmm1		,%%xmm0		\n\t"\
		"mulpd		%%xmm1		,%%xmm3		\n\t"\
		"movaps		%%xmm0		,    (%%rdx)\n\t"\
		"movaps		%%xmm3		,0x10(%%rdx)\n\t"\
		"movq		%[add0]		,%%rax		\n\t"\
		"movq		%[add1]		,%%r15		\n\t"\
		"movaps		(%%rax)		,%%xmm0		\n\t"\
		"movaps		0x10(%%rax)	,%%xmm1		\n\t"\
		"movaps		%%xmm0		,%%xmm2		\n\t"\
		"movaps		%%xmm1		,%%xmm3		\n\t"\
		"addpd		(%%r15)		,%%xmm2		\n\t"\
		"addpd		0x10(%%r15)	,%%xmm3		\n\t"\
		"subpd		(%%r15)		,%%xmm0		\n\t"\
		"subpd		0x10(%%r15)	,%%xmm1		\n\t"\
		"movaps		%%xmm6		,    (%%rax)\n\t"\
		"movaps		%%xmm7		,0x10(%%rax)\n\t"\
		"addpd		%%xmm6		,%%xmm6		\n\t"\
		"addpd		%%xmm7		,%%xmm7		\n\t"\
		"movq		%[add4]		,%%rcx		\n\t"\
		"movaps		%%xmm6		,    (%%rcx)\n\t"\
		"movaps		%%xmm7		,0x10(%%rcx)\n\t"\
		"movq		%[add2]		,%%rcx		\n\t"\
		"movq		%[add3]		,%%rdx		\n\t"\
		"movaps		(%%rcx)		,%%xmm4		\n\t"\
		"movaps		0x10(%%rcx)	,%%xmm5		\n\t"\
		"movaps		%%xmm4		,%%xmm6		\n\t"\
		"movaps		%%xmm5		,%%xmm7		\n\t"\
		"addpd		(%%rdx)		,%%xmm6		\n\t"\
		"addpd		0x10(%%rdx)	,%%xmm7		\n\t"\
		"subpd		(%%rdx)		,%%xmm4		\n\t"\
		"subpd		0x10(%%rdx)	,%%xmm5		\n\t"\
		"movaps		%%xmm6		,    (%%rcx)\n\t"\
		"movaps		%%xmm7		,0x10(%%rcx)\n\t"\
		"addpd		%%xmm2		,%%xmm6		\n\t"\
		"addpd		%%xmm3		,%%xmm7		\n\t"\
		"subpd		    (%%rcx)	,%%xmm2		\n\t"\
		"subpd		0x10(%%rcx)	,%%xmm3		\n\t"\
		"addpd		    (%%rax)	,%%xmm6		\n\t"\
		"addpd		0x10(%%rax)	,%%xmm7		\n\t"\
		"movaps		%%xmm6		,    (%%rax)\n\t"\
		"movaps		%%xmm7		,0x10(%%rax)\n\t"\
		"movq		%[add4]		,%%r15		\n\t"\
		"subpd		    (%%r15)	,%%xmm6		\n\t"\
		"subpd		0x10(%%r15)	,%%xmm7		\n\t"\
		"movaps		%%xmm6		,    (%%r15)\n\t"\
		"movaps		%%xmm7		,0x10(%%r15)\n\t"\
		"movaps		%%xmm6		,    (%%rdx)\n\t"\
		"movaps		%%xmm7		,0x10(%%rdx)\n\t"\
		"movaps		%%xmm4		,%%xmm6		\n\t"\
		"movaps		%%xmm5		,%%xmm7		\n\t"\
		"addpd		%%xmm0		,%%xmm5		\n\t"\
		"subpd		%%xmm7		,%%xmm0		\n\t"\
		"addpd		%%xmm1		,%%xmm4		\n\t"\
		"subpd		%%xmm6		,%%xmm1		\n\t"\
		"movq		%[c4]		,%%rcx		\n\t"\
		"movaps		    (%%r15)	,%%xmm6		\n\t"\
		"movaps		0x10(%%r15)	,%%xmm7		\n\t"\
		"mulpd		0x10(%%rcx)	,%%xmm6		\n\t"\
		"mulpd		0x10(%%rcx)	,%%xmm7		\n\t"\
		"movaps		%%xmm6		,    (%%r15)\n\t"\
		"movaps		%%xmm7		,0x10(%%r15)\n\t"\
		"movaps		    (%%rdx)	,%%xmm6		\n\t"\
		"movaps		0x10(%%rdx)	,%%xmm7		\n\t"\
		"mulpd		    (%%rcx)	,%%xmm6		\n\t"\
		"mulpd		    (%%rcx)	,%%xmm7		\n\t"\
		"addpd		0x10(%%r15)	,%%xmm6		\n\t"\
		"subpd		    (%%r15)	,%%xmm7		\n\t"\
		"movaps		%%xmm6		,    (%%r15)\n\t"\
		"movaps		%%xmm7		,0x10(%%r15)\n\t"\
		"movq		%[add1]		,%%rax		\n\t"\
		"movq		%[add5]		,%%r15		\n\t"\
		"movq		%[c1]		,%%rcx		\n\t"\
		"movq		%[c5]		,%%rdx		\n\t"\
		"movaps		%%xmm5		,%%xmm6		\n\t"\
		"movaps		%%xmm1		,%%xmm7		\n\t"\
		"addpd		    (%%r15)	,%%xmm5		\n\t"\
		"subpd		0x10(%%r15)	,%%xmm1		\n\t"\
		"subpd		    (%%r15)	,%%xmm6		\n\t"\
		"addpd		0x10(%%r15)	,%%xmm7		\n\t"\
		"movaps		%%xmm6		,    (%%r15)\n\t"\
		"movaps		%%xmm7		,0x10(%%r15)\n\t"\
		"movaps		%%xmm5		,%%xmm6		\n\t"\
		"movaps		%%xmm1		,%%xmm7		\n\t"\
		"mulpd		    (%%rcx)	,%%xmm5		\n\t"\
		"mulpd		    (%%rcx)	,%%xmm1		\n\t"\
		"mulpd		0x10(%%rcx)	,%%xmm6		\n\t"\
		"mulpd		0x10(%%rcx)	,%%xmm7		\n\t"\
		"subpd		%%xmm6		,%%xmm1		\n\t"\
		"addpd		%%xmm7		,%%xmm5		\n\t"\
		"movaps		%%xmm1		,0x10(%%rax)\n\t"\
		"movaps		%%xmm5		,    (%%rax)\n\t"\
		"movaps		    (%%r15)	,%%xmm5		\n\t"\
		"movaps		0x10(%%r15)	,%%xmm1		\n\t"\
		"movaps		%%xmm5		,%%xmm6		\n\t"\
		"movaps		%%xmm1		,%%xmm7		\n\t"\
		"mulpd		    (%%rdx)	,%%xmm5		\n\t"\
		"mulpd		    (%%rdx)	,%%xmm1		\n\t"\
		"mulpd		0x10(%%rdx)	,%%xmm6		\n\t"\
		"mulpd		0x10(%%rdx)	,%%xmm7		\n\t"\
		"subpd		%%xmm6		,%%xmm1		\n\t"\
		"addpd		%%xmm7		,%%xmm5		\n\t"\
		"movaps		%%xmm1		,0x10(%%r15)\n\t"\
		"movaps		%%xmm5		,    (%%r15)\n\t"\
		"movq		%[add2]		,%%rax		\n\t"\
		"movq		%[add6]		,%%r15		\n\t"\
		"movq		%[c2]		,%%rcx		\n\t"\
		"movq		%[c6]		,%%rdx		\n\t"\
		"movaps		%%xmm2		,%%xmm6		\n\t"\
		"movaps		%%xmm3		,%%xmm7		\n\t"\
		"addpd		0x10(%%r15)	,%%xmm2		\n\t"\
		"subpd		    (%%r15)	,%%xmm3		\n\t"\
		"subpd		0x10(%%r15)	,%%xmm6		\n\t"\
		"addpd		    (%%r15)	,%%xmm7		\n\t"\
		"movaps		%%xmm2		,%%xmm1		\n\t"\
		"movaps		%%xmm3		,%%xmm5		\n\t"\
		"mulpd		    (%%rcx)	,%%xmm2		\n\t"\
		"mulpd		    (%%rcx)	,%%xmm3		\n\t"\
		"mulpd		0x10(%%rcx)	,%%xmm1		\n\t"\
		"mulpd		0x10(%%rcx)	,%%xmm5		\n\t"\
		"subpd		%%xmm1		,%%xmm3		\n\t"\
		"addpd		%%xmm5		,%%xmm2		\n\t"\
		"movaps		%%xmm3		,0x10(%%rax)\n\t"\
		"movaps		%%xmm2		,    (%%rax)\n\t"\
		"movaps		%%xmm6		,%%xmm1		\n\t"\
		"movaps		%%xmm7		,%%xmm5		\n\t"\
		"mulpd		    (%%rdx)	,%%xmm6		\n\t"\
		"mulpd		    (%%rdx)	,%%xmm7		\n\t"\
		"mulpd		0x10(%%rdx)	,%%xmm1		\n\t"\
		"mulpd		0x10(%%rdx)	,%%xmm5		\n\t"\
		"subpd		%%xmm1		,%%xmm7		\n\t"\
		"addpd		%%xmm5		,%%xmm6		\n\t"\
		"movaps		%%xmm7		,0x10(%%r15)\n\t"\
		"movaps		%%xmm6		,    (%%r15)\n\t"\
		"movq		%[add3]		,%%rax		\n\t"\
		"movq		%[add7]		,%%r15		\n\t"\
		"movq		%[c3]		,%%rcx		\n\t"\
		"movq		%[c7]		,%%rdx		\n\t"\
		"movaps		%%xmm0		,%%xmm6		\n\t"\
		"movaps		%%xmm4		,%%xmm7		\n\t"\
		"subpd		0x10(%%r15)	,%%xmm0		\n\t"\
		"subpd		    (%%r15)	,%%xmm4		\n\t"\
		"addpd		0x10(%%r15)	,%%xmm6		\n\t"\
		"addpd		    (%%r15)	,%%xmm7		\n\t"\
		"movaps		%%xmm0		,%%xmm1		\n\t"\
		"movaps		%%xmm4		,%%xmm5		\n\t"\
		"mulpd		    (%%rcx)	,%%xmm0		\n\t"\
		"mulpd		    (%%rcx)	,%%xmm4		\n\t"\
		"mulpd		0x10(%%rcx)	,%%xmm1		\n\t"\
		"mulpd		0x10(%%rcx)	,%%xmm5		\n\t"\
		"subpd		%%xmm1		,%%xmm4		\n\t"\
		"addpd		%%xmm5		,%%xmm0		\n\t"\
		"movaps		%%xmm4		,0x10(%%rax)\n\t"\
		"movaps		%%xmm0		,    (%%rax)\n\t"\
		"movaps		%%xmm6		,%%xmm1		\n\t"\
		"movaps		%%xmm7		,%%xmm5		\n\t"\
		"mulpd		    (%%rdx)	,%%xmm6		\n\t"\
		"mulpd		    (%%rdx)	,%%xmm7		\n\t"\
		"mulpd		0x10(%%rdx)	,%%xmm1		\n\t"\
		"mulpd		0x10(%%rdx)	,%%xmm5		\n\t"\
		"subpd		%%xmm1		,%%xmm7		\n\t"\
		"addpd		%%xmm5		,%%xmm6		\n\t"\
		"movaps		%%xmm7		,0x10(%%r15)\n\t"\
		"movaps		%%xmm6		,    (%%r15)\n\t"\
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
		: "rax","r15","rcx","rdx","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
	);\
	}

#endif	/* radix8_dif_dit_pass_gcc_h_included */

