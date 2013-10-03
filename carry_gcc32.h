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
#ifndef carry_gcc_h_included
#define carry_gcc_h_included

	/************** See the Visual-studio-style 32-bit analogs of these in carry.h for commented versions: **********/

	/*************************************************************/
	/**************** FERMAT  -MOD CARRY MACROS ******************/
	/*************************************************************/

	#define SSE2_fermat_carry_norm_pow2_errcheck(Xdata,Xcy,Xnrt_bits,Xnrtm1,Xidx_offset,Xidx_incr,Xhalf_arr,Xsign_mask,Xadd1,Xadd2)\
	{\
	__asm__ volatile (\
	"pushl %%ebx	\n\t"/* Explicit save/restore of PIC register */\
		"movl	%[__idx_offset],%%esi	\n\t"\
		"movl		%[__nrt_bits],%%ecx	\n\t"\
		"movl		%%esi,%%eax			\n\t"\
		"shrl		$1,%%eax			\n\t"\
		"movl		%%eax,%%ebx			\n\t"\
		"andl		%[__nrtm1],%%eax	\n\t"\
		"shrl		%%cl,%%ebx			\n\t"\
		"shll		$4,%%eax			\n\t"\
		"shll		$4,%%ebx			\n\t"\
		"addl		%[__add1],%%eax		\n\t"\
		"addl		%[__add2],%%ebx		\n\t"\
		"movaps		(%%eax),%%xmm0		\n\t"\
		"movaps		(%%ebx),%%xmm1		\n\t"\
		"movl		%%esi,%%eax			\n\t"\
		"movaps		%%xmm1,%%xmm2		\n\t"\
		"shufpd	$1,	%%xmm2,%%xmm2		\n\t"\
		"mulpd		%%xmm0,%%xmm1		\n\t"\
		"mulpd		%%xmm0,%%xmm2		\n\t"\
		"addl		$2,%%eax			\n\t"\
		"shrl		$1,%%eax			\n\t"\
		"movl		%%eax,%%ebx			\n\t"\
		"andl		%[__nrtm1],%%eax	\n\t"\
		"shrl		%%cl,%%ebx			\n\t"\
		"shll		$4,%%eax			\n\t"\
		"shll		$4,%%ebx			\n\t"\
		"addl		%[__add1],%%eax		\n\t"\
		"addl		%[__add2],%%ebx		\n\t"\
		"movaps		(%%eax),%%xmm0		\n\t"\
		"movaps		(%%ebx),%%xmm3		\n\t"\
		"movl		%%esi,%%eax			\n\t"\
		"movaps		%%xmm3,%%xmm4		\n\t"\
		"shufpd	$1,	%%xmm4,%%xmm4		\n\t"\
		"mulpd		%%xmm0,%%xmm3		\n\t"\
		"mulpd		%%xmm0,%%xmm4		\n\t"\
		"movaps		%%xmm1,%%xmm0		\n\t"\
		"unpcklpd	%%xmm3,%%xmm0		\n\t"\
		"unpckhpd	%%xmm3,%%xmm1		\n\t"\
		"subpd		%%xmm1,%%xmm0		\n\t"\
		"movaps		%%xmm2,%%xmm1		\n\t"\
		"unpcklpd	%%xmm4,%%xmm1		\n\t"\
		"unpckhpd	%%xmm4,%%xmm2		\n\t"\
		"addpd		%%xmm2,%%xmm1		\n\t"\
		"movl		%[__half_arr],%%ecx	\n\t"\
		"movl		%[__data],%%edx		\n\t"\
		"movaps		     (%%edx),%%xmm4	\n\t"\
		"movaps		 0x10(%%edx),%%xmm2	\n\t"\
		"movaps		0x020(%%ecx),%%xmm5	\n\t"\
		"mulpd		%%xmm5,%%xmm4		\n\t"\
		"mulpd		%%xmm5,%%xmm2		\n\t"\
		"movaps		%%xmm4,%%xmm5		\n\t"\
		"movaps		%%xmm2,%%xmm3		\n\t"\
		"mulpd		%%xmm0,%%xmm4		\n\t"\
		"mulpd		%%xmm1,%%xmm3		\n\t"\
		"mulpd		%%xmm0,%%xmm2		\n\t"\
		"mulpd		%%xmm1,%%xmm5		\n\t"\
		"addpd		%%xmm3,%%xmm4		\n\t"\
		"subpd		%%xmm5,%%xmm2		\n\t"\
		"movl		%[__cy],%%ebx		\n\t"\
		"movaps		%%xmm4,%%xmm5		\n\t"\
		"shufpd	$0,	%%xmm2,%%xmm4		\n\t"\
		"shufpd	$3,	%%xmm2,%%xmm5		\n\t"\
		"addpd		     (%%ebx),%%xmm4	\n\t"\
		"movaps		-0x20(%%ecx),%%xmm6	\n\t"\
		"movaps		-0x10(%%ecx),%%xmm7	\n\t"\
		"movaps		%%xmm4,%%xmm2		\n\t"\
		"addpd		%%xmm7,%%xmm4		\n\t"\
		"subpd		%%xmm7,%%xmm4		\n\t"\
		"movl		%[__sign_mask],%%eax\n\t"\
		"subpd		%%xmm4,%%xmm2		\n\t"\
		"andpd		     (%%eax),%%xmm2	\n\t"\
		"maxpd		%%xmm6,%%xmm2		\n\t"\
		"movaps		%%xmm2,%%xmm6		\n\t"\
		"movaps		%%xmm4,%%xmm2		\n\t"\
		"mulpd		0x10(%%ecx),%%xmm2	\n\t"\
		"addpd		%%xmm7,%%xmm2		\n\t"\
		"subpd		%%xmm7,%%xmm2		\n\t"\
		"movaps		%%xmm2,%%xmm3		\n\t"\
		"mulpd		    (%%ecx),%%xmm3	\n\t"\
		"subpd		%%xmm3,%%xmm4		\n\t"\
		"addpd		%%xmm2,%%xmm5		\n\t"\
		"movaps		%%xmm5,%%xmm2		\n\t"\
		"addpd		%%xmm7,%%xmm5		\n\t"\
		"subpd		%%xmm7,%%xmm5		\n\t"\
		"subpd		%%xmm5,%%xmm2		\n\t"\
		"andpd		     (%%eax),%%xmm2	\n\t"\
		"maxpd		%%xmm6,%%xmm2		\n\t"\
		"movaps		%%xmm2,%%xmm6		\n\t"\
		"movaps		%%xmm5,%%xmm2		\n\t"\
		"mulpd		 0x10(%%ecx),%%xmm2	\n\t"\
		"addpd		%%xmm7,%%xmm2		\n\t"\
		"subpd		%%xmm7,%%xmm2		\n\t"\
		"movaps		%%xmm2,%%xmm3		\n\t"\
		"mulpd		     (%%ecx),%%xmm3	\n\t"\
		"subpd		%%xmm3,%%xmm5		\n\t"\
		"movaps		%%xmm2,(%%ebx)		\n\t"\
		"movaps		%%xmm4,%%xmm2		\n\t"\
		"shufpd	$0,	%%xmm5,%%xmm4		\n\t"\
		"shufpd	$3,	%%xmm5,%%xmm2		\n\t"\
		"movaps		%%xmm4,%%xmm5		\n\t"\
		"movaps		%%xmm2,%%xmm3		\n\t"\
		"movaps		%%xmm6,-0x20(%%ecx)	\n\t"\
		"mulpd		%%xmm0,%%xmm4		\n\t"\
		"mulpd		%%xmm1,%%xmm3		\n\t"\
		"mulpd		%%xmm0,%%xmm2		\n\t"\
		"mulpd		%%xmm1,%%xmm5		\n\t"\
		"subpd		%%xmm3,%%xmm4		\n\t"\
		"addpd		%%xmm2,%%xmm5		\n\t"\
		"movaps		%%xmm4,    (%%edx)	\n\t"\
		"movaps		%%xmm5,0x10(%%edx)	\n\t"\
		"addl	%[__idx_incr],%%esi		\n\t"\
		"movl	%%esi, %[__idx_offset]	/* Store incremented idx_offset */	\n\t"\
	"popl %%ebx	\n\t"\
	:						/* outputs: none */\
	:	[__data]		"m" (Xdata)	/* All inputs from memory addresses here */\
	,	[__cy]			"m" (Xcy)\
	,	[__nrt_bits]	"m" (Xnrt_bits)\
	,	[__nrtm1]		"m" (Xnrtm1)\
	,	[__idx_offset]	"m" (Xidx_offset)\
	,	[__idx_incr]	"m" (Xidx_incr)\
	,	[__half_arr]	"m" (Xhalf_arr)\
	,	[__sign_mask]	"m" (Xsign_mask)\
	,	[__add1]		"m" (Xadd1)\
	,	[__add2]		"m" (Xadd2)\
		: "cc","memory","eax","ecx","edx","esi","edi"	/* Clobbered registers */\
	);\
	}

	#define SSE2_fermat_carry_norm_errcheck(Xdata,Xcy,Xnrt_bits,Xnrtm1,Xidx_offset,Xidx_incr,Xodd_radix,Xhalf_arr,Xsign_mask,Xadd1,Xadd2,Xoffset0,Xoffset1)\
	{\
	__asm__ volatile (\
	"pushl %%ebx	\n\t"/* Explicit save/restore of PIC register */\
		"movl	%[__idx_offset],%%esi	\n\t"\
		"movl %[__odd_radix],%%edi	\n\t"\
		"movl	%[__nrt_bits],%%ecx		\n\t"\
		"movl		%%esi,%%eax			\n\t"\
		"shrl		$1,%%eax			\n\t"\
		"movl		%%eax,%%ebx			\n\t"\
		"andl		%[__nrtm1],%%eax	\n\t"\
		"shrl		%%cl,%%ebx			\n\t"\
		"shll		$4,%%eax			\n\t"\
		"shll		$4,%%ebx			\n\t"\
		"shll		$4,%%edi			\n\t"\
		"addl		%[__add1],%%eax		\n\t"\
		"addl		%[__add2],%%ebx		\n\t"\
		"movaps		(%%eax),%%xmm0		\n\t"\
		"movaps		(%%ebx),%%xmm1		\n\t"\
		"movl		%%esi,%%eax			\n\t"\
		"movaps		%%xmm1,%%xmm2		\n\t"\
		"shufpd	$1,	%%xmm2,%%xmm2		\n\t"\
		"mulpd		%%xmm0,%%xmm1		\n\t"\
		"mulpd		%%xmm0,%%xmm2		\n\t"\
		"addl		$2,%%eax			\n\t"\
		"shrl		$1,%%eax			\n\t"\
		"movl		%%eax,%%ebx			\n\t"\
		"andl		%[__nrtm1],%%eax	\n\t"\
		"shrl		%%cl,%%ebx			\n\t"\
		"shll		$4,%%eax			\n\t"\
		"shll		$4,%%ebx			\n\t"\
		"addl		%[__add1],%%eax		\n\t"\
		"addl		%[__add2],%%ebx		\n\t"\
		"movaps		(%%eax),%%xmm0		\n\t"\
		"movaps		(%%ebx),%%xmm3		\n\t"\
		"movl		%%esi,%%eax			\n\t"\
		"movaps		%%xmm3,%%xmm4		\n\t"\
		"shufpd	$1,	%%xmm4,%%xmm4		\n\t"\
		"mulpd		%%xmm0,%%xmm3		\n\t"\
		"mulpd		%%xmm0,%%xmm4		\n\t"\
		"movaps		%%xmm1,%%xmm0		\n\t"\
		"unpcklpd	%%xmm3,%%xmm0		\n\t"\
		"unpckhpd	%%xmm3,%%xmm1		\n\t"\
		"subpd		%%xmm1,%%xmm0		\n\t"\
		"movaps		%%xmm2,%%xmm1		\n\t"\
		"unpcklpd	%%xmm4,%%xmm1		\n\t"\
		"unpckhpd	%%xmm4,%%xmm2		\n\t"\
		"addpd		%%xmm2,%%xmm1		\n\t"\
		"movl		%[__half_arr],%%ecx	\n\t"\
		"movl		%[__data],%%edx		\n\t"\
		"movaps		     (%%edx),%%xmm4	\n\t"\
		"movaps		 0x10(%%edx),%%xmm2	\n\t"\
		"addl		%[__offset0],%%ecx	\n\t"\
		"/* movaps	(%%ecx,%%edi),%%xmm5	*/\n\t"\
		"addl	%%ecx,%%edi	\n\t"\
		"movaps	(%%ecx),%%xmm5	\n\t"\
		"subl	%%ecx,%%edi	\n\t"\
		"subl		%[__offset0],%%ecx	\n\t"\
		"mulpd		%%xmm5,%%xmm4		\n\t"\
		"mulpd		%%xmm5,%%xmm2		\n\t"\
		"movaps		%%xmm4,%%xmm5		\n\t"\
		"movaps		%%xmm2,%%xmm3		\n\t"\
		"mulpd		%%xmm0,%%xmm4		\n\t"\
		"mulpd		%%xmm1,%%xmm3		\n\t"\
		"mulpd		%%xmm0,%%xmm2		\n\t"\
		"mulpd		%%xmm1,%%xmm5		\n\t"\
		"addpd		%%xmm3,%%xmm4		\n\t"\
		"subpd		%%xmm5,%%xmm2		\n\t"\
		"movl		%[__cy],%%ebx		\n\t"\
		"movaps		%%xmm4,%%xmm5		\n\t"\
		"shufpd	$0,	%%xmm2,%%xmm4		\n\t"\
		"shufpd	$3,	%%xmm2,%%xmm5		\n\t"\
		"addpd		     (%%ebx),%%xmm4	\n\t"\
		"movaps		-0x20(%%ecx),%%xmm6	\n\t"\
		"movaps		-0x10(%%ecx),%%xmm7	\n\t"\
		"addl	   %[__offset0],%%ecx	\n\t"\
		"movaps		%%xmm4,%%xmm2		\n\t"\
		"shll	   $1,%%edi				\n\t"\
		"addpd		%%xmm7,%%xmm4		\n\t"\
		"subpd		%%xmm7,%%xmm4		\n\t"\
		"movl		%[__sign_mask],%%eax\n\t"\
		"subpd		%%xmm4,%%xmm2		\n\t"\
		"andpd		(%%eax),%%xmm2		\n\t"\
		"maxpd		%%xmm6,%%xmm2		\n\t"\
		"movaps		%%xmm2,%%xmm6		\n\t"\
		"addl		%%edi,%%ecx			\n\t"\
		"shrl		$1,%%edi			\n\t"\
		"movaps		%%xmm4,%%xmm2		\n\t"\
		"mulpd  (%%ecx,%%edi),%%xmm2	\n\t"\
		"addpd		%%xmm7,%%xmm2		\n\t"\
		"subpd		%%xmm7,%%xmm2		\n\t"\
		"movaps		%%xmm2,%%xmm3		\n\t"\
		"mulpd		(%%ecx),%%xmm3		\n\t"\
		"subl		%[__offset0],%%ecx	\n\t"\
		"subpd		%%xmm3,%%xmm4		\n\t"\
		"addpd		%%xmm2,%%xmm5		\n\t"\
		"movaps		%%xmm5,%%xmm2		\n\t"\
		"addpd		%%xmm7,%%xmm5		\n\t"\
		"subpd		%%xmm7,%%xmm5		\n\t"\
		"subpd		%%xmm5,%%xmm2		\n\t"\
		"andpd		(%%eax),%%xmm2		\n\t"\
		"maxpd		%%xmm6,%%xmm2		\n\t"\
		"movaps		%%xmm2,%%xmm6		\n\t"\
		"movaps		%%xmm5,%%xmm2		\n\t"\
		"addl		%[__offset1],%%ecx	\n\t"\
		"mulpd  (%%ecx,%%edi),%%xmm2	\n\t"\
		"addpd		%%xmm7,%%xmm2		\n\t"\
		"subpd		%%xmm7,%%xmm2		\n\t"\
		"shll		$1,%%edi			\n\t"\
		"movaps		%%xmm2,%%xmm3		\n\t"\
		"mulpd		(%%ecx),%%xmm3		\n\t"\
		"subl		%[__offset1],%%ecx	\n\t"\
		"subpd		%%xmm3,%%xmm5		\n\t"\
		"movaps		%%xmm2,(%%ebx)		\n\t"\
		"movaps		%%xmm4,%%xmm2		\n\t"\
		"shufpd	$0,	%%xmm5,%%xmm4		\n\t"\
		"shufpd	$3,	%%xmm5,%%xmm2		\n\t"\
		"movaps		%%xmm4,%%xmm5		\n\t"\
		"movaps		%%xmm2,%%xmm3		\n\t"\
		"subl		%%edi,%%ecx			\n\t"\
		"mulpd		%%xmm0,%%xmm4		\n\t"\
		"mulpd		%%xmm1,%%xmm3		\n\t"\
		"movaps	 	%%xmm6,-0x20(%%ecx)	\n\t"\
		"addl		%[__offset0],%%ecx	\n\t"\
		"mulpd		%%xmm0,%%xmm2		\n\t"\
		"mulpd		%%xmm1,%%xmm5		\n\t"\
		"movaps		(%%ecx),%%xmm0		\n\t"\
		"subpd		%%xmm3,%%xmm4		\n\t"\
		"addpd		%%xmm2,%%xmm5		\n\t"\
		"mulpd		%%xmm0,%%xmm4		\n\t"\
		"mulpd		%%xmm0,%%xmm5		\n\t"\
		"movaps		%%xmm4,    (%%edx)	\n\t"\
		"movaps		%%xmm5,0x10(%%edx)	\n\t"\
		"addl	%[__idx_incr],%%esi		\n\t"\
		"movl	%%esi, %[__idx_offset]	/* Store incremented idx_offset */	\n\t"\
	"popl %%ebx	\n\t"\
	:						/* outputs: none */\
	:	[__data]		"m" (Xdata)	/* All inputs from memory addresses here */\
	,	[__cy]			"m" (Xcy)\
	,	[__nrt_bits]	"m" (Xnrt_bits)\
	,	[__nrtm1]		"m" (Xnrtm1)\
	,	[__idx_offset]	"m" (Xidx_offset)\
	,	[__idx_incr]	"m" (Xidx_incr)\
	,	[__odd_radix]   "m" (Xodd_radix)\
	,	[__half_arr]	"m" (Xhalf_arr)\
	,	[__sign_mask]	"m" (Xsign_mask)\
	,	[__add1]		"m" (Xadd1)\
	,	[__add2]		"m" (Xadd2)\
	,	[__offset0]		"m" (Xoffset0)\
	,	[__offset1]		"m" (Xoffset1)\
		: "cc","memory","eax","ecx","edx","esi","edi"   /* Clobbered registers */\
	);\
	}

	/*************************************************************/
	/**************** MERSENNE-MOD CARRY MACROS ******************/
	/*************************************************************/

	/********* Packed 32-bit-int version of SSE2_cmplx_carry_norm_pow2_errcheck0_2x:***********/
	#define SSE2_cmplx_carry_norm_pow2_errcheck0_2B(Xdata,XwtA,XwtB,XwtC,XcyA,XcyB,Xbjmod_0,Xhalf_arr,Xi,Xn_minus_silp1,Xn_minus_sil,Xsign_mask,Xsinwt,Xsinwtm1,Xsse_bw,Xsse_nm1,Xsse_sw)\
	{\
	__asm__ volatile (\
	"pushl %%ebx	\n\t"/* Explicit save/restore of PIC register */\
	"/***************Unpack the data:*************************/\n\t"\
		"movl	%[__data]	,%%eax	\n\t"\
		"movaps		    (%%eax)	,%%xmm1	\n\t	movaps		0x40(%%eax)	,%%xmm5	\n\t"\
		"movaps		    (%%eax)	,%%xmm2	\n\t	movaps		0x40(%%eax)	,%%xmm6	\n\t"\
		"unpcklpd	0x20(%%eax)	,%%xmm1	\n\t	unpcklpd	0x60(%%eax)	,%%xmm5	\n\t"\
		"unpckhpd	0x20(%%eax)	,%%xmm2	\n\t	unpckhpd	0x60(%%eax)	,%%xmm6	\n\t"\
		"movaps		%%xmm2, 0x20(%%eax)	\n\t	movaps		%%xmm6, 0x60(%%eax)	\n\t"\
		"\n\t"\
		"movaps		0x10(%%eax)	,%%xmm2	\n\t	movaps		0x50(%%eax)	,%%xmm6	\n\t"\
		"movaps		0x10(%%eax)	,%%xmm3	\n\t	movaps		0x50(%%eax)	,%%xmm7	\n\t"\
		"unpcklpd	0x30(%%eax)	,%%xmm2	\n\t	unpcklpd	0x70(%%eax)	,%%xmm6	\n\t"\
		"unpckhpd	0x30(%%eax)	,%%xmm3	\n\t	unpckhpd	0x70(%%eax)	,%%xmm7	\n\t"\
		"movaps		%%xmm2, 0x10(%%eax)	\n\t	movaps		%%xmm6, 0x50(%%eax)	\n\t"\
		"movaps		%%xmm3, 0x30(%%eax)	\n\t	movaps		%%xmm7, 0x70(%%eax)	\n\t"\
	"/**********************************************/\n\t"\
	"/*          Real      parts                   */\n\t"\
	"/**********************************************/\n\t"\
		"movl	%[__bjmod_0],	%%eax	\n\t"\
		"movaps		(%%eax)	,	%%xmm0	\n\t"\
		"movl	%[__sse_sw]	,	%%ebx	\n\t"\
		"movaps		(%%ebx)	,	%%xmm7	\n\t"\
		"psubd		%%xmm0	,	%%xmm7	\n\t"\
		"movmskps	%%xmm7	,	%%esi	\n\t"\
		"movl	%[__i]	,%%ecx			\n\t"\
		"andl	$0xfffffffe,%%esi		\n\t"\
		"addl	%%ecx	,%%esi			\n\t"\
		"shll	$24		,%%esi			\n\t"\
		"movaps		%%xmm0	,%%xmm7		\n\t"\
		"movl	%[__n_minus_sil],%%ecx	\n\t"\
		"movd	%%ecx	,%%xmm2			\n\t"\
		"pshufd	$0,	%%xmm2	,%%xmm2		\n\t"\
		"psubd		%%xmm0	,%%xmm2		\n\t"\
		"movmskps	%%xmm2	,%%ecx		\n\t"\
		"shll	$16		,%%ecx			\n\t"\
		"addl	%%ecx	,%%esi			\n\t"\
		"movl	%[__sinwt]	,%%edx		\n\t"\
		"movd	%%edx	,%%xmm3			\n\t"\
		"pshufd	$0,	%%xmm3	,%%xmm3		\n\t"\
		"psubd		%%xmm3	,%%xmm7		\n\t"\
		"movmskps	%%xmm7	,%%edx		\n\t"\
		"shll	$8		,%%edx			\n\t"\
		"addl	%%edx	,%%esi			\n\t"\
		"\n\t"\
		"movl	%[__half_arr]	,%%eax		\n\t"\
		"movl	%[__wtA]	,%%ebx		\n\t"\
		"movl	%[__wtB]	,%%ecx		\n\t"\
		"movaps		-0x10(%%eax),%%xmm4	\n\t"\
		"movaps		     (%%ebx),%%xmm2	\n\t	movaps		 0x10(%%ebx),%%xmm6	\n\t"\
		"movhpd		     (%%ecx),%%xmm3	\n\t	movhpd		-0x10(%%ecx),%%xmm7	\n\t"\
		"movlpd		 0x08(%%ecx),%%xmm3	\n\t	movlpd		-0x08(%%ecx),%%xmm7	\n\t"\
		"\n\t"\
		"subl	$0x20	,%%ecx			\n\t"\
		"movl	%%ecx	,%[__wtB]		\n\t"\
		"\n\t"\
		"shld	$20,	%%esi	,%%edi	\n\t	shld	$18,	%%esi	,%%ebx	\n\t"\
		"andl	$0x00000030	,%%edi		\n\t	andl	$0x00000030	,%%ebx		\n\t"\
		"shld	$28,	%%esi,	%%edx	\n\t	shld	$26,	%%esi,	%%ecx	\n\t"\
		"andl	$0x00000030	,%%edx		\n\t	andl	$0x00000030	,%%ecx		\n\t"\
		"addl	%%eax	,%%edi			\n\t	addl	%%eax	,%%ebx			\n\t"\
		"addl	%%eax	,%%edx			\n\t	addl	%%eax	,%%ecx			\n\t"\
		"\n\t"\
		"mulpd	 	0x100(%%eax),%%xmm2	\n\t	mulpd	 	0x100(%%eax),%%xmm6	\n\t"\
		"mulpd	 	0x110(%%eax),%%xmm3	\n\t	mulpd	 	0x110(%%eax),%%xmm7	\n\t"\
		"mulpd	 	     (%%edi),%%xmm2	\n\t	mulpd	 	     (%%ebx),%%xmm6	\n\t"\
		"mulpd	 	0x040(%%edx),%%xmm3	\n\t	mulpd	 	0x040(%%ecx),%%xmm7	\n\t"\
		"\n\t"\
		"movl	%[__cyA]	,%%ecx		\n\t	movl	%[__cyB]	,%%edx		\n\t"\
		"mulpd		%%xmm3	,%%xmm1		\n\t	mulpd		%%xmm7	,%%xmm5		\n\t"\
		"addpd		(%%ecx)	,%%xmm1		\n\t	addpd		(%%edx)	,%%xmm5		\n\t"\
		"movaps		%%xmm1	,%%xmm3		\n\t	movaps		%%xmm5	,%%xmm7		\n\t"\
		"addpd		%%xmm4	,%%xmm3		\n\t	addpd		%%xmm4	,%%xmm7		\n\t"\
		"subpd		%%xmm4	,%%xmm3		\n\t	subpd		%%xmm4	,%%xmm7		\n\t"\
		"\n\t"\
		"movl	%[__sign_mask],%%ebx	\n\t"\
		"subpd			%%xmm3	,%%xmm1	\n\t	subpd			%%xmm7	,%%xmm5	\n\t"\
		"andpd			(%%ebx)	,%%xmm1	\n\t	andpd			(%%ebx)	,%%xmm5	\n\t"\
		"maxpd			%%xmm5	,%%xmm1	\n\t"\
		"maxpd		-0x20(%%eax),%%xmm1	\n\t"\
		"movaps		%%xmm1,-0x20(%%eax)	\n\t"\
		"\n\t"\
		"shld	$12		,%%esi	,%%edi	\n\t	shld	$10		,%%esi	,%%ebx	\n\t"\
		"andl	$0x00000030	,%%edi		\n\t	andl	$0x00000030	,%%ebx		\n\t"\
		"addl	%%eax	,%%edi			\n\t	addl	%%eax	,%%ebx			\n\t"\
		"movaps		%%xmm3	,%%xmm1		\n\t	movaps		%%xmm7	,%%xmm5		\n\t"\
		"mulpd	 	0xc0(%%edi)	,%%xmm3	\n\t	mulpd	 	0xc0(%%ebx)	,%%xmm7	\n\t"\
		"addpd		%%xmm4	,%%xmm3		\n\t	addpd		%%xmm4	,%%xmm7		\n\t"\
		"subpd		%%xmm4	,%%xmm3		\n\t	subpd		%%xmm4	,%%xmm7		\n\t"\
		"movaps		%%xmm3	,(%%ecx)	\n\t	movaps		%%xmm7	,(%%edx)	\n\t"\
		"\n\t"\
		"\n\t"\
		"movl	%[__data]	,%%eax		\n\t"\
		"mulpd	 0x80(%%edi)	,%%xmm3	\n\t	mulpd	 0x80(%%ebx)	,%%xmm7	\n\t"\
		"subpd		%%xmm3	,	%%xmm1	\n\t	subpd		%%xmm7	,	%%xmm5	\n\t"\
		"mulpd		%%xmm2	,	%%xmm1	\n\t	mulpd		%%xmm6	,	%%xmm5	\n\t"\
		"movaps		%%xmm1	,    (%%eax)\n\t	movaps		%%xmm5	,0x40(%%eax)\n\t"\
		"\n\t"\
		"\n\t"\
		"movl	%[__sse_bw]	,%%eax		\n\t"\
		"movl	%[__sse_nm1]	,%%ebx	\n\t"\
		"paddd		(%%eax)	,%%xmm0		\n\t"\
		"pand		(%%ebx)	,%%xmm0		\n\t"\
	"/**********************************************/\n\t"\
	"/*          Imaginary parts               */\n\t"\
	"/**********************************************/\n\t"\
		"movl	%[__sse_sw]	,%%edx		\n\t"\
		"movaps	(%%edx)	,%%xmm1			\n\t"\
		"psubd	%%xmm0	,%%xmm1			\n\t"\
		"movmskps	%%xmm1	,%%esi		\n\t"\
		"shll	$24	,%%esi				\n\t"\
		"movaps	%%xmm0	,%%xmm1			\n\t"\
		"movl	%[__n_minus_silp1],%%ecx\n\t"\
		"movd	%%ecx	,%%xmm2			\n\t"\
		"pshufd	$0	,%%xmm2	,%%xmm2		\n\t"\
		"psubd	%%xmm0	,%%xmm2			\n\t"\
		"movmskps	%%xmm2	,%%ecx		\n\t"\
		"shll	$16	,%%ecx				\n\t"\
		"addl	%%ecx	,%%esi			\n\t"\
		"movl	%[__sinwtm1]	,%%edx	\n\t"\
		"movd	%%edx	,%%xmm3			\n\t"\
		"pshufd	$0	,%%xmm3	,%%xmm3		\n\t"\
		"psubd	%%xmm3	,%%xmm1			\n\t"\
		"movmskps	%%xmm1	,%%edx		\n\t"\
		"shll	$8	,%%edx				\n\t"\
		"addl	%%edx	,%%esi			\n\t"\
		"movl	%[__data]	,%%eax		\n\t"\
		"movaps	 0x10(%%eax)	,%%xmm1	\n\t"\
		"movaps	 0x50(%%eax)	,%%xmm5	\n\t"\
		"\n\t"\
		"movl	%[__half_arr]	,%%eax	\n\t"\
		"movl	%[__wtA]	,%%ebx		\n\t"\
		"movl	%[__wtC]	,%%ecx		\n\t"\
		"\n\t"\
		"movaps	     (%%ebx)	,%%xmm2	\n\t	movaps	 0x10(%%ebx)	,%%xmm6	\n\t"\
		"movhpd	     (%%ecx)	,%%xmm3	\n\t	movhpd	-0x10(%%ecx)	,%%xmm7	\n\t"\
		"movlpd	 0x08(%%ecx)	,%%xmm3	\n\t	movlpd	-0x08(%%ecx)	,%%xmm7	\n\t"\
		"\n\t"\
		"addl	$0x20	,%%ebx			\n\t"\
		"subl	$0x20	,%%ecx			\n\t"\
		"movl	%%ebx	,%[__wtA]		\n\t"\
		"movl	%%ecx	,%[__wtC]		\n\t"\
		"shld	$20		,%%esi	,%%edi	\n\t	shld	$18		,%%esi	,%%ebx	\n\t"\
		"andl	$0x00000030	,%%edi		\n\t	andl	$0x00000030	,%%ebx		\n\t"\
		"shld	$28,	%%esi,	%%edx	\n\t	shld	$26,	%%esi,	%%ecx	\n\t"\
		"andl	$0x00000030	,%%edx		\n\t	andl	$0x00000030	,%%ecx		\n\t"\
		"addl	%%eax	,%%edi			\n\t	addl	%%eax	,%%ebx			\n\t"\
		"addl	%%eax	,%%edx			\n\t	addl	%%eax	,%%ecx			\n\t"\
		"\n\t"\
		"mulpd	 0x120(%%eax)	,%%xmm2	\n\t	mulpd	 0x120(%%eax)	,%%xmm6	\n\t"\
		"mulpd	 0x130(%%eax)	,%%xmm3	\n\t	mulpd	 0x130(%%eax)	,%%xmm7	\n\t"\
		"mulpd	      (%%edi)	,%%xmm2	\n\t	mulpd	      (%%ebx)	,%%xmm6	\n\t"\
		"mulpd	 0x040(%%edx)	,%%xmm3	\n\t	mulpd	 0x040(%%ecx)	,%%xmm7	\n\t"\
		"\n\t"\
		"movl	%[__cyA]	,%%ecx		\n\t	movl	%[__cyB]	,%%edx		\n\t"\
		"mulpd	%%xmm3	,%%xmm1			\n\t	mulpd	%%xmm7	,%%xmm5			\n\t"\
		"addpd	(%%ecx)	,%%xmm1			\n\t	addpd	(%%edx)	,%%xmm5			\n\t"\
		"movaps	%%xmm1	,%%xmm3			\n\t	movaps	%%xmm5	,%%xmm7			\n\t"\
		"addpd	%%xmm4	,%%xmm3			\n\t	addpd	%%xmm4	,%%xmm7			\n\t"\
		"subpd	%%xmm4	,%%xmm3			\n\t	subpd	%%xmm4	,%%xmm7			\n\t"\
		"\n\t"\
		"movl	%[__sign_mask],%%ebx	\n\t"\
		"subpd	%%xmm3	,%%xmm1			\n\t	subpd	%%xmm7	,%%xmm5			\n\t"\
		"andpd	(%%ebx)	,%%xmm1			\n\t	andpd	(%%ebx)	,%%xmm5			\n\t"\
		"maxpd	%%xmm5	,%%xmm1			\n\t"\
		"maxpd	-0x20(%%eax)	,%%xmm1	\n\t"\
		"movaps	%%xmm1	,-0x20(%%eax)	\n\t"\
		"\n\t"\
		"shld	$12		,%%esi	,%%edi	\n\t	shld	$10		,%%esi	,%%ebx	\n\t"\
		"andl	$0x00000030	,%%edi		\n\t	andl	$0x00000030	,%%ebx		\n\t"\
		"addl	%%eax	,%%edi			\n\t	addl	%%eax	,%%ebx			\n\t"\
		"movaps	%%xmm3	,%%xmm1			\n\t	movaps	%%xmm7	,%%xmm5			\n\t"\
		"mulpd	 0xc0(%%edi)	,%%xmm3	\n\t	mulpd	 0xc0(%%ebx)	,%%xmm7	\n\t"\
		"addpd	%%xmm4	,%%xmm3			\n\t	addpd	%%xmm4	,%%xmm7			\n\t"\
		"subpd	%%xmm4	,%%xmm3			\n\t	subpd	%%xmm4	,%%xmm7			\n\t"\
		"movaps	%%xmm3	,(%%ecx)		\n\t	movaps	%%xmm7	,(%%edx)		\n\t"\
		"\n\t"\
		"\n\t"\
		"movl	%[__data]	,%%eax		\n\t"\
		"mulpd	 0x80(%%edi)	,%%xmm3	\n\t	mulpd	 0x80(%%ebx)	,%%xmm7	\n\t"\
		"subpd	%%xmm3	,%%xmm1			\n\t	subpd	%%xmm7	,%%xmm5			\n\t"\
		"mulpd	%%xmm2	,%%xmm1			\n\t	mulpd	%%xmm6	,%%xmm5			\n\t"\
		"movaps	%%xmm1	, 0x10(%%eax)	\n\t	movaps	%%xmm5	, 0x50(%%eax)	\n\t"\
		"\n\t"\
		"\n\t"\
		"movl	%[__sse_bw]	,%%eax		\n\t"\
		"movl	%[__sse_nm1]	,%%ebx	\n\t"\
		"paddd	(%%eax)	,%%xmm0			\n\t"\
		"pand	(%%ebx)	,%%xmm0			\n\t"\
		"movl	%[__bjmod_0],%%ecx		\n\t"\
		"movaps	%%xmm0,(%%ecx)			\n\t"\
	"popl %%ebx	\n\t"\
		:					/* outputs: none */\
		: [__data]		"m" (Xdata)	/* All inputs from memory addresses here */\
		, [__wtA]		"m" (XwtA)		\
		, [__wtB]		"m" (XwtB)		\
		, [__wtC]		"m" (XwtC)		\
		, [__cyA]		"m" (XcyA)		\
		, [__cyB]		"m" (XcyB)		\
		, [__bjmod_0]	"m" (Xbjmod_0)		\
		, [__half_arr]	"m" (Xhalf_arr)		\
		, [__i]			"m" (Xi)			\
		, [__n_minus_silp1] "m" (Xn_minus_silp1)\
		, [__n_minus_sil]	"m" (Xn_minus_sil)	\
		, [__sign_mask]	"m" (Xsign_mask)	\
		, [__sinwt]		"m" (Xsinwt)		\
		, [__sinwtm1]		"m" (Xsinwtm1)		\
		, [__sse_bw]		"m" (Xsse_bw)		\
		, [__sse_nm1]		"m" (Xsse_nm1)		\
		, [__sse_sw]		"m" (Xsse_sw)		\
		: "cc","memory","eax","ecx","edx","edi","esi"		/* Clobbered registers */\
	);\
	}

	/********* Packed 32-bit-int version of SSE2_cmplx_carry_norm_pow2_errcheck1_2x:***********/
	#define SSE2_cmplx_carry_norm_pow2_errcheck1_2B(Xdata,XwtA,XwtB,XwtC,XcyA,XcyB,Xbjmod_0,Xhalf_arr,Xn_minus_silp1,Xn_minus_sil,Xsign_mask,Xsinwt,Xsinwtm1,Xsse_bw,Xsse_nm1,Xsse_sw)\
	{\
	__asm__ volatile (\
	"pushl %%ebx	\n\t"/* Explicit save/restore of PIC register */\
	"/***************Unpack the data:*************************/\n\t"\
		"movl	%[__data]	,%%eax	\n\t"\
		"movaps		    (%%eax)	,%%xmm1	\n\t	movaps		0x40(%%eax)	,%%xmm5	\n\t"\
		"movaps		    (%%eax)	,%%xmm2	\n\t	movaps		0x40(%%eax)	,%%xmm6	\n\t"\
		"unpcklpd	0x20(%%eax)	,%%xmm1	\n\t	unpcklpd	0x60(%%eax)	,%%xmm5	\n\t"\
		"unpckhpd	0x20(%%eax)	,%%xmm2	\n\t	unpckhpd	0x60(%%eax)	,%%xmm6	\n\t"\
		"movaps		%%xmm2, 0x20(%%eax)	\n\t	movaps		%%xmm6, 0x60(%%eax)	\n\t"\
		"\n\t"\
		"movaps		0x10(%%eax)	,%%xmm2	\n\t	movaps		0x50(%%eax)	,%%xmm6	\n\t"\
		"movaps		0x10(%%eax)	,%%xmm3	\n\t	movaps		0x50(%%eax)	,%%xmm7	\n\t"\
		"unpcklpd	0x30(%%eax)	,%%xmm2	\n\t	unpcklpd	0x70(%%eax)	,%%xmm6	\n\t"\
		"unpckhpd	0x30(%%eax)	,%%xmm3	\n\t	unpckhpd	0x70(%%eax)	,%%xmm7	\n\t"\
		"movaps		%%xmm2, 0x10(%%eax)	\n\t	movaps		%%xmm6, 0x50(%%eax)	\n\t"\
		"movaps		%%xmm3, 0x30(%%eax)	\n\t	movaps		%%xmm7, 0x70(%%eax)	\n\t"\
	"/**********************************************/\n\t"\
	"/*          Real      parts                   */\n\t"\
	"/**********************************************/\n\t"\
		"movl	%[__bjmod_0],	%%eax	\n\t"\
		"movaps		(%%eax)	,	%%xmm0	\n\t"\
		"movl	%[__sse_sw]	,	%%ebx	\n\t"\
		"movaps		(%%ebx)	,	%%xmm7	\n\t"\
		"psubd		%%xmm0	,	%%xmm7	\n\t"\
		"movmskps	%%xmm7	,	%%esi	\n\t"\
		"\n\t"\
		"shll	$24		,%%esi			\n\t"\
		"movaps		%%xmm0	,%%xmm7		\n\t"\
		"movl	%[__n_minus_sil],%%ecx	\n\t"\
		"movd	%%ecx	,%%xmm2			\n\t"\
		"pshufd	$0,	%%xmm2	,%%xmm2		\n\t"\
		"psubd		%%xmm0	,%%xmm2		\n\t"\
		"movmskps	%%xmm2	,%%ecx		\n\t"\
		"shll	$16		,%%ecx			\n\t"\
		"addl	%%ecx	,%%esi			\n\t"\
		"movl	%[__sinwt]	,%%edx		\n\t"\
		"movd	%%edx	,%%xmm3			\n\t"\
		"pshufd	$0,	%%xmm3	,%%xmm3		\n\t"\
		"psubd		%%xmm3	,%%xmm7		\n\t"\
		"movmskps	%%xmm7	,%%edx		\n\t"\
		"shll	$8		,%%edx			\n\t"\
		"addl	%%edx	,%%esi			\n\t"\
		"\n\t"\
		"movl	%[__half_arr]	,%%eax	\n\t"\
		"movl	%[__wtA]	,%%ebx		\n\t"\
		"movl	%[__wtB]	,%%ecx		\n\t"\
		"movaps		-0x10(%%eax),%%xmm4	\n\t"\
		"movaps		     (%%ebx),%%xmm2	\n\t	movaps		 0x10(%%ebx),%%xmm6	\n\t"\
		"movhpd		     (%%ecx),%%xmm3	\n\t	movhpd		-0x10(%%ecx),%%xmm7	\n\t"\
		"movlpd		 0x08(%%ecx),%%xmm3	\n\t	movlpd		-0x08(%%ecx),%%xmm7	\n\t"\
		"\n\t"\
		"subl	$0x20	,%%ecx			\n\t"\
		"movl	%%ecx	,%[__wtB]		\n\t"\
		"\n\t"\
		"shld	$20,	%%esi	,%%edi	\n\t	shld	$18,	%%esi	,%%ebx	\n\t"\
		"andl	$0x00000030	,%%edi		\n\t	andl	$0x00000030	,%%ebx		\n\t"\
		"shld	$28,	%%esi,	%%edx	\n\t	shld	$26,	%%esi,	%%ecx	\n\t"\
		"andl	$0x00000030	,%%edx		\n\t	andl	$0x00000030	,%%ecx		\n\t"\
		"addl	%%eax	,%%edi			\n\t	addl	%%eax	,%%ebx			\n\t"\
		"addl	%%eax	,%%edx			\n\t	addl	%%eax	,%%ecx			\n\t"\
		"\n\t"\
		"mulpd	 	0x100(%%eax),%%xmm2	\n\t	mulpd	 	0x100(%%eax),%%xmm6	\n\t"\
		"mulpd	 	0x110(%%eax),%%xmm3	\n\t	mulpd	 	0x110(%%eax),%%xmm7	\n\t"\
		"mulpd	 	     (%%edi),%%xmm2	\n\t	mulpd	 	     (%%ebx),%%xmm6	\n\t"\
		"mulpd	 	0x040(%%edx),%%xmm3	\n\t	mulpd	 	0x040(%%ecx),%%xmm7	\n\t"\
		"\n\t"\
		"movl	%[__cyA]	,%%ecx		\n\t	movl	%[__cyB]	,%%edx		\n\t"\
		"mulpd		%%xmm3	,%%xmm1		\n\t	mulpd		%%xmm7	,%%xmm5		\n\t"\
		"addpd		(%%ecx)	,%%xmm1		\n\t	addpd		(%%edx)	,%%xmm5		\n\t"\
		"movaps		%%xmm1	,%%xmm3		\n\t	movaps		%%xmm5	,%%xmm7		\n\t"\
		"addpd		%%xmm4	,%%xmm3		\n\t	addpd		%%xmm4	,%%xmm7		\n\t"\
		"subpd		%%xmm4	,%%xmm3		\n\t	subpd		%%xmm4	,%%xmm7		\n\t"\
		"\n\t"\
		"movl	%[__sign_mask],%%ebx	\n\t"\
		"subpd			%%xmm3	,%%xmm1	\n\t	subpd			%%xmm7	,%%xmm5	\n\t"\
		"andpd			(%%ebx)	,%%xmm1	\n\t	andpd			(%%ebx)	,%%xmm5	\n\t"\
		"maxpd			%%xmm5	,%%xmm1	\n\t"\
		"maxpd		-0x20(%%eax),%%xmm1	\n\t"\
		"movaps		%%xmm1,-0x20(%%eax)	\n\t"\
		"\n\t"\
		"shld	$12		,%%esi	,%%edi	\n\t	shld	$10		,%%esi	,%%ebx	\n\t"\
		"andl	$0x00000030	,%%edi		\n\t	andl	$0x00000030	,%%ebx		\n\t"\
		"addl	%%eax	,%%edi			\n\t	addl	%%eax	,%%ebx			\n\t"\
		"movaps		%%xmm3	,%%xmm1		\n\t	movaps		%%xmm7	,%%xmm5		\n\t"\
		"mulpd	 	0xc0(%%edi)	,%%xmm3	\n\t	mulpd	 	0xc0(%%ebx)	,%%xmm7	\n\t"\
		"addpd		%%xmm4	,%%xmm3		\n\t	addpd		%%xmm4	,%%xmm7		\n\t"\
		"subpd		%%xmm4	,%%xmm3		\n\t	subpd		%%xmm4	,%%xmm7		\n\t"\
		"movaps		%%xmm3	,(%%ecx)	\n\t	movaps		%%xmm7	,(%%edx)	\n\t"\
		"\n\t"\
		"\n\t"\
		"movl	%[__data]	,%%eax		\n\t"\
		"mulpd	 0x80(%%edi)	,%%xmm3	\n\t	mulpd	 0x80(%%ebx)	,%%xmm7	\n\t"\
		"subpd		%%xmm3	,	%%xmm1	\n\t	subpd		%%xmm7	,	%%xmm5	\n\t"\
		"mulpd		%%xmm2	,	%%xmm1	\n\t	mulpd		%%xmm6	,	%%xmm5	\n\t"\
		"movaps		%%xmm1	,    (%%eax)\n\t	movaps		%%xmm5	,0x40(%%eax)\n\t"\
		"\n\t"\
		"\n\t"\
		"movl	%[__sse_bw]	,%%eax		\n\t"\
		"movl	%[__sse_nm1]	,%%ebx	\n\t"\
		"paddd		(%%eax)	,%%xmm0		\n\t"\
		"pand		(%%ebx)	,%%xmm0		\n\t"\
	"/**********************************************/\n\t"\
	"/*          Imaginary parts               */\n\t"\
	"/**********************************************/\n\t"\
		"movl	%[__sse_sw]	,%%edx		\n\t"\
		"movaps	(%%edx)	,%%xmm1			\n\t"\
		"psubd	%%xmm0	,%%xmm1			\n\t"\
		"movmskps	%%xmm1	,%%esi		\n\t"\
		"shll	$24	,%%esi				\n\t"\
		"movaps	%%xmm0	,%%xmm1			\n\t"\
		"movl	%[__n_minus_silp1],%%ecx\n\t"\
		"movd	%%ecx	,%%xmm2			\n\t"\
		"pshufd	$0	,%%xmm2	,%%xmm2		\n\t"\
		"psubd	%%xmm0	,%%xmm2			\n\t"\
		"movmskps	%%xmm2	,%%ecx		\n\t"\
		"shll	$16	,%%ecx				\n\t"\
		"addl	%%ecx	,%%esi			\n\t"\
		"movl	%[__sinwtm1]	,%%edx	\n\t"\
		"movd	%%edx	,%%xmm3			\n\t"\
		"pshufd	$0	,%%xmm3	,%%xmm3		\n\t"\
		"psubd	%%xmm3	,%%xmm1			\n\t"\
		"movmskps	%%xmm1	,%%edx		\n\t"\
		"shll	$8	,%%edx				\n\t"\
		"addl	%%edx	,%%esi			\n\t"\
		"movl	%[__data]	,%%eax		\n\t"\
		"movaps	 0x10(%%eax)	,%%xmm1	\n\t"\
		"movaps	 0x50(%%eax)	,%%xmm5	\n\t"\
		"\n\t"\
		"movl	%[__half_arr]	,%%eax	\n\t"\
		"movl	%[__wtA]	,%%ebx		\n\t"\
		"movl	%[__wtC]	,%%ecx		\n\t"\
		"\n\t"\
		"movaps	     (%%ebx)	,%%xmm2	\n\t	movaps	 0x10(%%ebx)	,%%xmm6	\n\t"\
		"movhpd	     (%%ecx)	,%%xmm3	\n\t	movhpd	-0x10(%%ecx)	,%%xmm7	\n\t"\
		"movlpd	 0x08(%%ecx)	,%%xmm3	\n\t	movlpd	-0x08(%%ecx)	,%%xmm7	\n\t"\
		"\n\t"\
		"addl	$0x20	,%%ebx			\n\t"\
		"subl	$0x20	,%%ecx			\n\t"\
		"movl	%%ebx	,%[__wtA]		\n\t"\
		"movl	%%ecx	,%[__wtC]		\n\t"\
		"shld	$20		,%%esi	,%%edi	\n\t	shld	$18		,%%esi	,%%ebx	\n\t"\
		"andl	$0x00000030	,%%edi		\n\t	andl	$0x00000030	,%%ebx		\n\t"\
		"shld	$28,	%%esi,	%%edx	\n\t	shld	$26,	%%esi,	%%ecx	\n\t"\
		"andl	$0x00000030	,%%edx		\n\t	andl	$0x00000030	,%%ecx		\n\t"\
		"addl	%%eax	,%%edi			\n\t	addl	%%eax	,%%ebx			\n\t"\
		"addl	%%eax	,%%edx			\n\t	addl	%%eax	,%%ecx			\n\t"\
		"\n\t"\
		"mulpd	 0x120(%%eax)	,%%xmm2	\n\t	mulpd	 0x120(%%eax)	,%%xmm6	\n\t"\
		"mulpd	 0x130(%%eax)	,%%xmm3	\n\t	mulpd	 0x130(%%eax)	,%%xmm7	\n\t"\
		"mulpd	      (%%edi)	,%%xmm2	\n\t	mulpd	      (%%ebx)	,%%xmm6	\n\t"\
		"mulpd	 0x040(%%edx)	,%%xmm3	\n\t	mulpd	 0x040(%%ecx)	,%%xmm7	\n\t"\
		"\n\t"\
		"movl	%[__cyA]	,%%ecx		\n\t	movl	%[__cyB]	,%%edx		\n\t"\
		"mulpd	%%xmm3	,%%xmm1			\n\t	mulpd	%%xmm7	,%%xmm5			\n\t"\
		"addpd	(%%ecx)	,%%xmm1			\n\t	addpd	(%%edx)	,%%xmm5			\n\t"\
		"movaps	%%xmm1	,%%xmm3			\n\t	movaps	%%xmm5	,%%xmm7			\n\t"\
		"addpd	%%xmm4	,%%xmm3			\n\t	addpd	%%xmm4	,%%xmm7			\n\t"\
		"subpd	%%xmm4	,%%xmm3			\n\t	subpd	%%xmm4	,%%xmm7			\n\t"\
		"\n\t"\
		"movl	%[__sign_mask],%%ebx	\n\t"\
		"subpd	%%xmm3	,%%xmm1			\n\t	subpd	%%xmm7	,%%xmm5			\n\t"\
		"andpd	(%%ebx)	,%%xmm1			\n\t	andpd	(%%ebx)	,%%xmm5			\n\t"\
		"maxpd	%%xmm5	,%%xmm1			\n\t"\
		"maxpd	-0x20(%%eax)	,%%xmm1	\n\t"\
		"movaps	%%xmm1	,-0x20(%%eax)	\n\t"\
		"\n\t"\
		"shld	$12		,%%esi	,%%edi	\n\t	shld	$10		,%%esi	,%%ebx	\n\t"\
		"andl	$0x00000030	,%%edi		\n\t	andl	$0x00000030	,%%ebx		\n\t"\
		"addl	%%eax	,%%edi			\n\t	addl	%%eax	,%%ebx			\n\t"\
		"movaps	%%xmm3	,%%xmm1			\n\t	movaps	%%xmm7	,%%xmm5			\n\t"\
		"mulpd	 0xc0(%%edi)	,%%xmm3	\n\t	mulpd	 0xc0(%%ebx)	,%%xmm7	\n\t"\
		"addpd	%%xmm4	,%%xmm3			\n\t	addpd	%%xmm4	,%%xmm7			\n\t"\
		"subpd	%%xmm4	,%%xmm3			\n\t	subpd	%%xmm4	,%%xmm7			\n\t"\
		"movaps	%%xmm3	,(%%ecx)		\n\t	movaps	%%xmm7	,(%%edx)		\n\t"\
		"\n\t"\
		"\n\t"\
		"movl	%[__data]	,%%eax		\n\t"\
		"mulpd	 0x80(%%edi)	,%%xmm3	\n\t	mulpd	 0x80(%%ebx)	,%%xmm7	\n\t"\
		"subpd	%%xmm3	,%%xmm1			\n\t	subpd	%%xmm7	,%%xmm5			\n\t"\
		"mulpd	%%xmm2	,%%xmm1			\n\t	mulpd	%%xmm6	,%%xmm5			\n\t"\
		"movaps	%%xmm1	, 0x10(%%eax)	\n\t	movaps	%%xmm5	, 0x50(%%eax)	\n\t"\
		"\n\t"\
		"\n\t"\
		"movl	%[__sse_bw]	,%%eax		\n\t"\
		"movl	%[__sse_nm1]	,%%ebx	\n\t"\
		"paddd	(%%eax)	,%%xmm0			\n\t"\
		"pand	(%%ebx)	,%%xmm0			\n\t"\
		"movl	%[__bjmod_0],%%ecx		\n\t"\
		"movaps	%%xmm0,(%%ecx)			\n\t"\
	"popl %%ebx	\n\t"\
		:					/* outputs: none */\
		: [__data]		"m" (Xdata)	/* All inputs from memory addresses here */\
		, [__wtA]		"m" (XwtA)		\
		, [__wtB]		"m" (XwtB)		\
		, [__wtC]		"m" (XwtC)		\
		, [__cyA]		"m" (XcyA)		\
		, [__cyB]		"m" (XcyB)		\
		, [__bjmod_0]	"m" (Xbjmod_0)		\
		, [__half_arr]	"m" (Xhalf_arr)		\
		, [__n_minus_silp1] "m" (Xn_minus_silp1)\
		, [__n_minus_sil]	"m" (Xn_minus_sil)	\
		, [__sign_mask]	"m" (Xsign_mask)	\
		, [__sinwt]		"m" (Xsinwt)		\
		, [__sinwtm1]		"m" (Xsinwtm1)		\
		, [__sse_bw]		"m" (Xsse_bw)		\
		, [__sse_nm1]		"m" (Xsse_nm1)		\
		, [__sse_sw]		"m" (Xsse_sw)		\
		: "cc","memory","eax","ecx","edx","edi","esi"		/* Clobbered registers */\
	);\
	}

	/********* Packed 32-bit-int version of SSE2_cmplx_carry_norm_pow2_errcheck2_2x:***********/
	#define SSE2_cmplx_carry_norm_pow2_errcheck2_2B(Xdata,XwtA,XwtB,XcyA,XcyB,Xbjmod_0,Xhalf_arr,Xn_minus_silp1,Xn_minus_sil,Xsign_mask,Xsinwt,Xsinwtm1,Xsse_bw,Xsse_nm1,Xsse_sw)\
	{\
	__asm__ volatile (\
	"pushl %%ebx	\n\t"/* Explicit save/restore of PIC register */\
	"/**********************************************/\n\t"\
	"/*          Real      parts                   */\n\t"\
	"/**********************************************/\n\t"\
		"movl	%[__bjmod_0],	%%eax	\n\t"\
		"movaps		(%%eax)	,	%%xmm0	\n\t"\
		"movl	%[__sse_sw]	,	%%ebx	\n\t"\
		"movaps		(%%ebx)	,	%%xmm1	\n\t"\
		"psubd		%%xmm0	,	%%xmm1	\n\t"\
		"movmskps	%%xmm1	,	%%esi	\n\t"\
		"\n\t"\
		"shll	$24		,%%esi			\n\t"\
		"movaps		%%xmm0	,%%xmm1		\n\t"\
		"movl	%[__n_minus_sil],%%ecx	\n\t"\
		"movd	%%ecx	,%%xmm2			\n\t"\
		"pshufd	$0,	%%xmm2	,%%xmm2		\n\t"\
		"psubd		%%xmm0	,%%xmm2		\n\t"\
		"movmskps	%%xmm2	,%%ecx		\n\t"\
		"shll	$16		,%%ecx			\n\t"\
		"addl	%%ecx	,%%esi			\n\t"\
		"movl	%[__sinwt]	,%%edx		\n\t"\
		"movd	%%edx	,%%xmm3			\n\t"\
		"pshufd	$0,	%%xmm3	,%%xmm3		\n\t"\
		"psubd		%%xmm3	,%%xmm1		\n\t"\
		"movmskps	%%xmm1	,%%edx		\n\t"\
		"shll	$8		,%%edx			\n\t"\
		"addl	%%edx	,%%esi			\n\t"\
		"\n\t"\
		"mov	%[__data],%%eax			\n\t"\
		"movaps	 0x20(%%eax),%%xmm1		\n\t	movaps		 0x60(%%eax),%%xmm5\n\t"\
		"\n\t"\
		"movl	%[__half_arr]	,%%eax	\n\t"\
		"movl	%[__wtA]	,%%ebx		\n\t"\
		"movl	%[__wtB]	,%%ecx		\n\t"\
		"movaps		-0x10(%%eax),%%xmm4	\n\t"\
		"movaps		     (%%ebx),%%xmm2	\n\t	movaps		 0x10(%%ebx),%%xmm6	\n\t"\
		"movhpd		     (%%ecx),%%xmm3	\n\t	movhpd		-0x10(%%ecx),%%xmm7	\n\t"\
		"movlpd		 0x08(%%ecx),%%xmm3	\n\t	movlpd		-0x08(%%ecx),%%xmm7	\n\t"\
		"\n\t"\
		"shld	$20,	%%esi	,%%edi	\n\t	shld	$18,	%%esi	,%%ebx	\n\t"\
		"andl	$0x00000030	,%%edi		\n\t	andl	$0x00000030	,%%ebx		\n\t"\
		"shld	$28,	%%esi,	%%edx	\n\t	shld	$26,	%%esi,	%%ecx	\n\t"\
		"andl	$0x00000030	,%%edx		\n\t	andl	$0x00000030	,%%ecx		\n\t"\
		"addl	%%eax	,%%edi			\n\t	addl	%%eax	,%%ebx			\n\t"\
		"addl	%%eax	,%%edx			\n\t	addl	%%eax	,%%ecx			\n\t"\
		"\n\t"\
		"mulpd	 	0x100(%%eax),%%xmm2	\n\t	mulpd	 	0x100(%%eax),%%xmm6	\n\t"\
		"mulpd	 	0x110(%%eax),%%xmm3	\n\t	mulpd	 	0x110(%%eax),%%xmm7	\n\t"\
		"mulpd	 	     (%%edi),%%xmm2	\n\t	mulpd	 	     (%%ebx),%%xmm6	\n\t"\
		"mulpd	 	0x040(%%edx),%%xmm3	\n\t	mulpd	 	0x040(%%ecx),%%xmm7	\n\t"\
		"\n\t"\
		"movl	%[__cyA]	,%%ecx		\n\t	movl	%[__cyB]	,%%edx		\n\t"\
		"mulpd		%%xmm3	,%%xmm1		\n\t	mulpd		%%xmm7	,%%xmm5		\n\t"\
		"addpd		(%%ecx)	,%%xmm1		\n\t	addpd		(%%edx)	,%%xmm5		\n\t"\
		"movaps		%%xmm1	,%%xmm3		\n\t	movaps		%%xmm5	,%%xmm7		\n\t"\
		"addpd		%%xmm4	,%%xmm3		\n\t	addpd		%%xmm4	,%%xmm7		\n\t"\
		"subpd		%%xmm4	,%%xmm3		\n\t	subpd		%%xmm4	,%%xmm7		\n\t"\
		"\n\t"\
		"movl	%[__sign_mask],%%ebx	\n\t"\
		"subpd			%%xmm3	,%%xmm1	\n\t	subpd			%%xmm7	,%%xmm5	\n\t"\
		"andpd			(%%ebx)	,%%xmm1	\n\t	andpd			(%%ebx)	,%%xmm5	\n\t"\
		"maxpd			%%xmm5	,%%xmm1	\n\t"\
		"maxpd		-0x20(%%eax),%%xmm1	\n\t"\
		"movaps		%%xmm1,-0x20(%%eax)	\n\t"\
		"\n\t"\
		"shld	$12		,%%esi	,%%edi	\n\t	shld	$10		,%%esi	,%%ebx	\n\t"\
		"andl	$0x00000030	,%%edi		\n\t	andl	$0x00000030	,%%ebx		\n\t"\
		"addl	%%eax	,%%edi			\n\t	addl	%%eax	,%%ebx			\n\t"\
		"movaps		%%xmm3	,%%xmm1		\n\t	movaps		%%xmm7	,%%xmm5		\n\t"\
		"mulpd	 	0xc0(%%edi)	,%%xmm3	\n\t	mulpd	 	0xc0(%%ebx)	,%%xmm7	\n\t"\
		"addpd		%%xmm4	,%%xmm3		\n\t	addpd		%%xmm4	,%%xmm7		\n\t"\
		"subpd		%%xmm4	,%%xmm3		\n\t	subpd		%%xmm4	,%%xmm7		\n\t"\
		"movaps		%%xmm3	,(%%ecx)	\n\t	movaps		%%xmm7	,(%%edx)	\n\t"\
		"\n\t"\
		"\n\t"\
		"movl	%[__data]	,%%eax		\n\t"\
		"mulpd	 0x80(%%edi)	,%%xmm3	\n\t	mulpd	 0x80(%%ebx)	,%%xmm7	\n\t"\
		"subpd		%%xmm3	,	%%xmm1	\n\t	subpd		%%xmm7	,	%%xmm5	\n\t"\
		"mulpd		%%xmm2	,	%%xmm1	\n\t	mulpd		%%xmm6	,	%%xmm5	\n\t"\
		"movaps		%%xmm1	,0x20(%%eax)\n\t	movaps		%%xmm5	,0x60(%%eax)\n\t"\
		"\n\t"\
		"\n\t"\
		"movl	%[__sse_bw]	,%%eax		\n\t"\
		"movl	%[__sse_nm1]	,%%ebx	\n\t"\
		"paddd		(%%eax)	,%%xmm0		\n\t"\
		"pand		(%%ebx)	,%%xmm0		\n\t"\
	"/**********************************************/\n\t"\
	"/*          Imaginary parts               */\n\t"\
	"/**********************************************/\n\t"\
		"movl	%[__sse_sw]	,%%edx		\n\t"\
		"movaps	(%%edx)	,%%xmm1			\n\t"\
		"psubd	%%xmm0	,%%xmm1			\n\t"\
		"movmskps	%%xmm1	,%%esi		\n\t"\
		"shll	$24	,%%esi				\n\t"\
		"movaps	%%xmm0	,%%xmm1			\n\t"\
		"movl	%[__n_minus_silp1],%%ecx\n\t"\
		"movd	%%ecx	,%%xmm2			\n\t"\
		"pshufd	$0	,%%xmm2	,%%xmm2		\n\t"\
		"psubd	%%xmm0	,%%xmm2			\n\t"\
		"movmskps	%%xmm2	,%%ecx		\n\t"\
		"shll	$16	,%%ecx				\n\t"\
		"addl	%%ecx	,%%esi			\n\t"\
		"movl	%[__sinwtm1]	,%%edx	\n\t"\
		"movd	%%edx	,%%xmm3			\n\t"\
		"pshufd	$0	,%%xmm3	,%%xmm3		\n\t"\
		"psubd	%%xmm3	,%%xmm1			\n\t"\
		"movmskps	%%xmm1	,%%edx		\n\t"\
		"shll	$8	,%%edx				\n\t"\
		"addl	%%edx	,%%esi			\n\t"\
		"movl	%[__data]	,%%eax		\n\t"\
		"movaps	 0x30(%%eax)	,%%xmm1	\n\t"\
		"movaps	 0x70(%%eax)	,%%xmm5	\n\t"\
		"\n\t"\
		"movl	%[__half_arr]	,%%eax	\n\t"\
		"movl	%[__wtA]	,%%ebx		\n\t"\
		"movl	%[__wtB]	,%%ecx		\n\t"\
		"\n\t"\
		"movaps	     (%%ebx)	,%%xmm2	\n\t	movaps	 0x10(%%ebx)	,%%xmm6	\n\t"\
		"movhpd	     (%%ecx)	,%%xmm3	\n\t	movhpd	-0x10(%%ecx)	,%%xmm7	\n\t"\
		"movlpd	 0x08(%%ecx)	,%%xmm3	\n\t	movlpd	-0x08(%%ecx)	,%%xmm7	\n\t"\
		"\n\t"\
		"addl	$0x20	,%%ebx			\n\t"\
		"subl	$0x20	,%%ecx			\n\t"\
		"movl	%%ebx	,%[__wtA]		\n\t"\
		"movl	%%ecx	,%[__wtB]		\n\t"\
		"shld	$20		,%%esi	,%%edi	\n\t	shld	$18		,%%esi	,%%ebx	\n\t"\
		"andl	$0x00000030	,%%edi		\n\t	andl	$0x00000030	,%%ebx		\n\t"\
		"shld	$28,	%%esi,	%%edx	\n\t	shld	$26,	%%esi,	%%ecx	\n\t"\
		"andl	$0x00000030	,%%edx		\n\t	andl	$0x00000030	,%%ecx		\n\t"\
		"addl	%%eax	,%%edi			\n\t	addl	%%eax	,%%ebx			\n\t"\
		"addl	%%eax	,%%edx			\n\t	addl	%%eax	,%%ecx			\n\t"\
		"\n\t"\
		"mulpd	 0x120(%%eax)	,%%xmm2	\n\t	mulpd	 0x120(%%eax)	,%%xmm6	\n\t"\
		"mulpd	 0x130(%%eax)	,%%xmm3	\n\t	mulpd	 0x130(%%eax)	,%%xmm7	\n\t"\
		"mulpd	      (%%edi)	,%%xmm2	\n\t	mulpd	      (%%ebx)	,%%xmm6	\n\t"\
		"mulpd	 0x040(%%edx)	,%%xmm3	\n\t	mulpd	 0x040(%%ecx)	,%%xmm7	\n\t"\
		"\n\t"\
		"movl	%[__cyA]	,%%ecx		\n\t	movl	%[__cyB]	,%%edx		\n\t"\
		"mulpd	%%xmm3	,%%xmm1			\n\t	mulpd	%%xmm7	,%%xmm5			\n\t"\
		"addpd	(%%ecx)	,%%xmm1			\n\t	addpd	(%%edx)	,%%xmm5			\n\t"\
		"movaps	%%xmm1	,%%xmm3			\n\t	movaps	%%xmm5	,%%xmm7			\n\t"\
		"addpd	%%xmm4	,%%xmm3			\n\t	addpd	%%xmm4	,%%xmm7			\n\t"\
		"subpd	%%xmm4	,%%xmm3			\n\t	subpd	%%xmm4	,%%xmm7			\n\t"\
		"\n\t"\
		"movl	%[__sign_mask],%%ebx	\n\t"\
		"subpd	%%xmm3	,%%xmm1			\n\t	subpd	%%xmm7	,%%xmm5			\n\t"\
		"andpd	(%%ebx)	,%%xmm1			\n\t	andpd	(%%ebx)	,%%xmm5			\n\t"\
		"maxpd	%%xmm5	,%%xmm1			\n\t"\
		"maxpd	-0x20(%%eax)	,%%xmm1	\n\t"\
		"movaps	%%xmm1	,-0x20(%%eax)	\n\t"\
		"\n\t"\
		"shld	$12		,%%esi	,%%edi	\n\t	shld	$10		,%%esi	,%%ebx	\n\t"\
		"andl	$0x00000030	,%%edi		\n\t	andl	$0x00000030	,%%ebx		\n\t"\
		"addl	%%eax	,%%edi			\n\t	addl	%%eax	,%%ebx			\n\t"\
		"movaps	%%xmm3	,%%xmm1			\n\t	movaps	%%xmm7	,%%xmm5			\n\t"\
		"mulpd	 0xc0(%%edi)	,%%xmm3	\n\t	mulpd	 0xc0(%%ebx)	,%%xmm7	\n\t"\
		"addpd	%%xmm4	,%%xmm3			\n\t	addpd	%%xmm4	,%%xmm7			\n\t"\
		"subpd	%%xmm4	,%%xmm3			\n\t	subpd	%%xmm4	,%%xmm7			\n\t"\
		"movaps	%%xmm3	,(%%ecx)		\n\t	movaps	%%xmm7	,(%%edx)		\n\t"\
		"\n\t"\
		"\n\t"\
		"movl	%[__data]	,%%eax		\n\t"\
		"mulpd	 0x80(%%edi)	,%%xmm3	\n\t	mulpd	 0x80(%%ebx)	,%%xmm7	\n\t"\
		"subpd	%%xmm3	,%%xmm1			\n\t	subpd	%%xmm7	,%%xmm5			\n\t"\
		"mulpd	%%xmm2	,%%xmm1			\n\t	mulpd	%%xmm6	,%%xmm5			\n\t"\
		"movaps	%%xmm1	, 0x30(%%eax)	\n\t	movaps	%%xmm5	, 0x70(%%eax)	\n\t"\
		"\n\t"\
		"\n\t"\
		"movl	%[__sse_bw]	,%%eax		\n\t"\
		"movl	%[__sse_nm1]	,%%ebx	\n\t"\
		"paddd	(%%eax)	,%%xmm0			\n\t"\
		"pand	(%%ebx)	,%%xmm0			\n\t"\
		"movl	%[__bjmod_0],%%ecx		\n\t"\
		"movaps	%%xmm0,(%%ecx)			\n\t"\
	"/**********************************************/\n\t"\
	"/*          Repack the data:              */\n\t"\
	"/**********************************************/\n\t"\
		"mov	%[__data],%%eax			\n\t"\
		"movaps		0x10(%%eax)	,%%xmm1	\n\t	movaps		0x50(%%eax)	,%%xmm5	\n\t"\
		"movaps		    (%%eax)	,%%xmm0	\n\t	movaps		0x40(%%eax)	,%%xmm4	\n\t"\
		"movaps		%%xmm1		,%%xmm3	\n\t	movaps		%%xmm5		,%%xmm7	\n\t"\
		"movaps		%%xmm0		,%%xmm2	\n\t	movaps		%%xmm4		,%%xmm6	\n\t"\
		"unpckhpd	0x30(%%eax)	,%%xmm3	\n\t	unpckhpd	0x70(%%eax)	,%%xmm7	\n\t"\
		"unpcklpd	0x30(%%eax)	,%%xmm1	\n\t	unpcklpd	0x70(%%eax)	,%%xmm5	\n\t"\
		"movaps		%%xmm3,0x30(%%eax)	\n\t	movaps		%%xmm7,0x70(%%eax)	\n\t"\
		"unpckhpd	0x20(%%eax)	,%%xmm2	\n\t	unpckhpd	0x60(%%eax)	,%%xmm6	\n\t"\
		"unpcklpd	0x20(%%eax)	,%%xmm0	\n\t	unpcklpd	0x60(%%eax)	,%%xmm4	\n\t"\
		"movaps		%%xmm2,0x20(%%eax)	\n\t	movaps		%%xmm6,0x60(%%eax)	\n\t"\
		"movaps		%%xmm1,0x10(%%eax)	\n\t	movaps		%%xmm5,0x50(%%eax)	\n\t"\
		"movaps		%%xmm0,    (%%eax)	\n\t	movaps		%%xmm4,0x40(%%eax)	\n\t"\
	"popl %%ebx	\n\t"\
		:					/* outputs: none */\
		: [__data]		"m" (Xdata)	/* All inputs from memory addresses here */\
		, [__wtA]		"m" (XwtA)		\
		, [__wtB]		"m" (XwtB)		\
		, [__cyA]		"m" (XcyA)		\
		, [__cyB]		"m" (XcyB)		\
		, [__bjmod_0]	"m" (Xbjmod_0)		\
		, [__half_arr]	"m" (Xhalf_arr)		\
		, [__i]			"m" (Xi)			\
		, [__n_minus_silp1] "m" (Xn_minus_silp1)\
		, [__n_minus_sil]	"m" (Xn_minus_sil)	\
		, [__sign_mask]	"m" (Xsign_mask)	\
		, [__sinwt]		"m" (Xsinwt)		\
		, [__sinwtm1]		"m" (Xsinwtm1)		\
		, [__sse_bw]		"m" (Xsse_bw)		\
		, [__sse_nm1]		"m" (Xsse_nm1)		\
		, [__sse_sw]		"m" (Xsse_sw)		\
		: "cc","memory","eax","ecx","edx","edi","esi"		/* Clobbered registers */\
	);\
	}

	/******************************************************************************************************************************************************************/
	/********** No-ROE-Check versions of the latter 2 of the above 3 macros - the first is called too infrequently to bother with a special non-ROE version: **********/
	/******************************************************************************************************************************************************************/

	#define SSE2_cmplx_carry_norm_pow2_nocheck1_2B(Xdata,XwtA,XwtB,XwtC,XcyA,XcyB,Xbjmod_0,Xhalf_arr,Xn_minus_silp1,Xn_minus_sil,Xsinwt,Xsinwtm1,Xsse_bw,Xsse_nm1,Xsse_sw)\
	{\
	__asm__ volatile (\
	"pushl %%ebx	\n\t"/* Explicit save/restore of PIC register */\
	"/***************Unpack the data:*************************/\n\t"\
		"movl	%[__data]	,%%eax	\n\t"\
		"movaps		    (%%eax)	,%%xmm1	\n\t	movaps		0x40(%%eax)	,%%xmm5	\n\t"\
		"movaps		    (%%eax)	,%%xmm2	\n\t	movaps		0x40(%%eax)	,%%xmm6	\n\t"\
		"unpcklpd	0x20(%%eax)	,%%xmm1	\n\t	unpcklpd	0x60(%%eax)	,%%xmm5	\n\t"\
		"unpckhpd	0x20(%%eax)	,%%xmm2	\n\t	unpckhpd	0x60(%%eax)	,%%xmm6	\n\t"\
		"movaps		%%xmm2, 0x20(%%eax)	\n\t	movaps		%%xmm6, 0x60(%%eax)	\n\t"\
		"\n\t"\
		"movaps		0x10(%%eax)	,%%xmm2	\n\t	movaps		0x50(%%eax)	,%%xmm6	\n\t"\
		"movaps		0x10(%%eax)	,%%xmm3	\n\t	movaps		0x50(%%eax)	,%%xmm7	\n\t"\
		"unpcklpd	0x30(%%eax)	,%%xmm2	\n\t	unpcklpd	0x70(%%eax)	,%%xmm6	\n\t"\
		"unpckhpd	0x30(%%eax)	,%%xmm3	\n\t	unpckhpd	0x70(%%eax)	,%%xmm7	\n\t"\
		"movaps		%%xmm2, 0x10(%%eax)	\n\t	movaps		%%xmm6, 0x50(%%eax)	\n\t"\
		"movaps		%%xmm3, 0x30(%%eax)	\n\t	movaps		%%xmm7, 0x70(%%eax)	\n\t"\
	"/**********************************************/\n\t"\
	"/*          Real      parts                   */\n\t"\
	"/**********************************************/\n\t"\
		"movl	%[__bjmod_0],	%%eax	\n\t"\
		"movaps		(%%eax)	,	%%xmm0	\n\t"\
		"movl	%[__sse_sw]	,	%%ebx	\n\t"\
		"movaps		(%%ebx)	,	%%xmm7	\n\t"\
		"psubd		%%xmm0	,	%%xmm7	\n\t"\
		"movmskps	%%xmm7	,	%%esi	\n\t"\
		"\n\t"\
		"shll	$24		,%%esi			\n\t"\
		"movaps		%%xmm0	,%%xmm7		\n\t"\
		"movl	%[__n_minus_sil],%%ecx	\n\t"\
		"movd	%%ecx	,%%xmm2			\n\t"\
		"pshufd	$0,	%%xmm2	,%%xmm2		\n\t"\
		"psubd		%%xmm0	,%%xmm2		\n\t"\
		"movmskps	%%xmm2	,%%ecx		\n\t"\
		"shll	$16		,%%ecx			\n\t"\
		"addl	%%ecx	,%%esi			\n\t"\
		"movl	%[__sinwt]	,%%edx		\n\t"\
		"movd	%%edx	,%%xmm3			\n\t"\
		"pshufd	$0,	%%xmm3	,%%xmm3		\n\t"\
		"psubd		%%xmm3	,%%xmm7		\n\t"\
		"movmskps	%%xmm7	,%%edx		\n\t"\
		"shll	$8		,%%edx			\n\t"\
		"addl	%%edx	,%%esi			\n\t"\
		"\n\t"\
		"movl	%[__half_arr]	,%%eax	\n\t"\
		"movl	%[__wtA]	,%%ebx		\n\t"\
		"movl	%[__wtB]	,%%ecx		\n\t"\
		"movaps		-0x10(%%eax),%%xmm4	\n\t"\
		"movaps		     (%%ebx),%%xmm2	\n\t	movaps		 0x10(%%ebx),%%xmm6	\n\t"\
		"movhpd		     (%%ecx),%%xmm3	\n\t	movhpd		-0x10(%%ecx),%%xmm7	\n\t"\
		"movlpd		 0x08(%%ecx),%%xmm3	\n\t	movlpd		-0x08(%%ecx),%%xmm7	\n\t"\
		"\n\t"\
		"subl	$0x20	,%%ecx			\n\t"\
		"movl	%%ecx	,%[__wtB]		\n\t"\
		"\n\t"\
		"shld	$20,	%%esi	,%%edi	\n\t	shld	$18,	%%esi	,%%ebx	\n\t"\
		"andl	$0x00000030	,%%edi		\n\t	andl	$0x00000030	,%%ebx		\n\t"\
		"shld	$28,	%%esi,	%%edx	\n\t	shld	$26,	%%esi,	%%ecx	\n\t"\
		"andl	$0x00000030	,%%edx		\n\t	andl	$0x00000030	,%%ecx		\n\t"\
		"addl	%%eax	,%%edi			\n\t	addl	%%eax	,%%ebx			\n\t"\
		"addl	%%eax	,%%edx			\n\t	addl	%%eax	,%%ecx			\n\t"\
		"\n\t"\
		"mulpd	 	0x100(%%eax),%%xmm2	\n\t	mulpd	 	0x100(%%eax),%%xmm6	\n\t"\
		"mulpd	 	0x110(%%eax),%%xmm3	\n\t	mulpd	 	0x110(%%eax),%%xmm7	\n\t"\
		"mulpd	 	     (%%edi),%%xmm2	\n\t	mulpd	 	     (%%ebx),%%xmm6	\n\t"\
		"mulpd	 	0x040(%%edx),%%xmm3	\n\t	mulpd	 	0x040(%%ecx),%%xmm7	\n\t"\
		"\n\t"\
		"movl	%[__cyA]	,%%ecx		\n\t	movl	%[__cyB]	,%%edx		\n\t"\
		"mulpd		%%xmm3	,%%xmm1		\n\t	mulpd		%%xmm7	,%%xmm5		\n\t"\
		"addpd		(%%ecx)	,%%xmm1		\n\t	addpd		(%%edx)	,%%xmm5		\n\t"\
		"movaps		%%xmm1	,%%xmm3		\n\t	movaps		%%xmm5	,%%xmm7		\n\t"\
		"addpd		%%xmm4	,%%xmm3		\n\t	addpd		%%xmm4	,%%xmm7		\n\t"\
		"subpd		%%xmm4	,%%xmm3		\n\t	subpd		%%xmm4	,%%xmm7		\n\t"\
		"\n\t"\
		"/* NO ROE HERE */		\n\t"\
		"\n\t"\
		"shld	$12		,%%esi	,%%edi	\n\t	shld	$10		,%%esi	,%%ebx	\n\t"\
		"andl	$0x00000030	,%%edi		\n\t	andl	$0x00000030	,%%ebx		\n\t"\
		"addl	%%eax	,%%edi			\n\t	addl	%%eax	,%%ebx			\n\t"\
		"movaps		%%xmm3	,%%xmm1		\n\t	movaps		%%xmm7	,%%xmm5		\n\t"\
		"mulpd	 	0xc0(%%edi)	,%%xmm3	\n\t	mulpd	 	0xc0(%%ebx)	,%%xmm7	\n\t"\
		"addpd		%%xmm4	,%%xmm3		\n\t	addpd		%%xmm4	,%%xmm7		\n\t"\
		"subpd		%%xmm4	,%%xmm3		\n\t	subpd		%%xmm4	,%%xmm7		\n\t"\
		"movaps		%%xmm3	,(%%ecx)	\n\t	movaps		%%xmm7	,(%%edx)	\n\t"\
		"\n\t"\
		"\n\t"\
		"movl	%[__data]	,%%eax		\n\t"\
		"mulpd	 0x80(%%edi)	,%%xmm3	\n\t	mulpd	 0x80(%%ebx)	,%%xmm7	\n\t"\
		"subpd		%%xmm3	,	%%xmm1	\n\t	subpd		%%xmm7	,	%%xmm5	\n\t"\
		"mulpd		%%xmm2	,	%%xmm1	\n\t	mulpd		%%xmm6	,	%%xmm5	\n\t"\
		"movaps		%%xmm1	,    (%%eax)\n\t	movaps		%%xmm5	,0x40(%%eax)\n\t"\
		"\n\t"\
		"\n\t"\
		"movl	%[__sse_bw]	,%%eax		\n\t"\
		"movl	%[__sse_nm1]	,%%ebx	\n\t"\
		"paddd		(%%eax)	,%%xmm0		\n\t"\
		"pand		(%%ebx)	,%%xmm0		\n\t"\
	"/**********************************************/\n\t"\
	"/*          Imaginary parts               */\n\t"\
	"/**********************************************/\n\t"\
		"movl	%[__sse_sw]	,%%edx		\n\t"\
		"movaps	(%%edx)	,%%xmm1			\n\t"\
		"psubd	%%xmm0	,%%xmm1			\n\t"\
		"movmskps	%%xmm1	,%%esi		\n\t"\
		"shll	$24	,%%esi				\n\t"\
		"movaps	%%xmm0	,%%xmm1			\n\t"\
		"movl	%[__n_minus_silp1],%%ecx\n\t"\
		"movd	%%ecx	,%%xmm2			\n\t"\
		"pshufd	$0	,%%xmm2	,%%xmm2		\n\t"\
		"psubd	%%xmm0	,%%xmm2			\n\t"\
		"movmskps	%%xmm2	,%%ecx		\n\t"\
		"shll	$16	,%%ecx				\n\t"\
		"addl	%%ecx	,%%esi			\n\t"\
		"movl	%[__sinwtm1]	,%%edx	\n\t"\
		"movd	%%edx	,%%xmm3			\n\t"\
		"pshufd	$0	,%%xmm3	,%%xmm3		\n\t"\
		"psubd	%%xmm3	,%%xmm1			\n\t"\
		"movmskps	%%xmm1	,%%edx		\n\t"\
		"shll	$8	,%%edx				\n\t"\
		"addl	%%edx	,%%esi			\n\t"\
		"movl	%[__data]	,%%eax		\n\t"\
		"movaps	 0x10(%%eax)	,%%xmm1	\n\t"\
		"movaps	 0x50(%%eax)	,%%xmm5	\n\t"\
		"\n\t"\
		"movl	%[__half_arr]	,%%eax	\n\t"\
		"movl	%[__wtA]	,%%ebx		\n\t"\
		"movl	%[__wtC]	,%%ecx		\n\t"\
		"\n\t"\
		"movaps	     (%%ebx)	,%%xmm2	\n\t	movaps	 0x10(%%ebx)	,%%xmm6	\n\t"\
		"movhpd	     (%%ecx)	,%%xmm3	\n\t	movhpd	-0x10(%%ecx)	,%%xmm7	\n\t"\
		"movlpd	 0x08(%%ecx)	,%%xmm3	\n\t	movlpd	-0x08(%%ecx)	,%%xmm7	\n\t"\
		"\n\t"\
		"addl	$0x20	,%%ebx			\n\t"\
		"subl	$0x20	,%%ecx			\n\t"\
		"movl	%%ebx	,%[__wtA]		\n\t"\
		"movl	%%ecx	,%[__wtC]		\n\t"\
		"shld	$20		,%%esi	,%%edi	\n\t	shld	$18		,%%esi	,%%ebx	\n\t"\
		"andl	$0x00000030	,%%edi		\n\t	andl	$0x00000030	,%%ebx		\n\t"\
		"shld	$28,	%%esi,	%%edx	\n\t	shld	$26,	%%esi,	%%ecx	\n\t"\
		"andl	$0x00000030	,%%edx		\n\t	andl	$0x00000030	,%%ecx		\n\t"\
		"addl	%%eax	,%%edi			\n\t	addl	%%eax	,%%ebx			\n\t"\
		"addl	%%eax	,%%edx			\n\t	addl	%%eax	,%%ecx			\n\t"\
		"\n\t"\
		"mulpd	 0x120(%%eax)	,%%xmm2	\n\t	mulpd	 0x120(%%eax)	,%%xmm6	\n\t"\
		"mulpd	 0x130(%%eax)	,%%xmm3	\n\t	mulpd	 0x130(%%eax)	,%%xmm7	\n\t"\
		"mulpd	      (%%edi)	,%%xmm2	\n\t	mulpd	      (%%ebx)	,%%xmm6	\n\t"\
		"mulpd	 0x040(%%edx)	,%%xmm3	\n\t	mulpd	 0x040(%%ecx)	,%%xmm7	\n\t"\
		"\n\t"\
		"movl	%[__cyA]	,%%ecx		\n\t	movl	%[__cyB]	,%%edx		\n\t"\
		"mulpd	%%xmm3	,%%xmm1			\n\t	mulpd	%%xmm7	,%%xmm5			\n\t"\
		"addpd	(%%ecx)	,%%xmm1			\n\t	addpd	(%%edx)	,%%xmm5			\n\t"\
		"movaps	%%xmm1	,%%xmm3			\n\t	movaps	%%xmm5	,%%xmm7			\n\t"\
		"addpd	%%xmm4	,%%xmm3			\n\t	addpd	%%xmm4	,%%xmm7			\n\t"\
		"subpd	%%xmm4	,%%xmm3			\n\t	subpd	%%xmm4	,%%xmm7			\n\t"\
		"\n\t"\
		"/* NO ROE HERE */		\n\t"\
		"\n\t"\
		"shld	$12		,%%esi	,%%edi	\n\t	shld	$10		,%%esi	,%%ebx	\n\t"\
		"andl	$0x00000030	,%%edi		\n\t	andl	$0x00000030	,%%ebx		\n\t"\
		"addl	%%eax	,%%edi			\n\t	addl	%%eax	,%%ebx			\n\t"\
		"movaps	%%xmm3	,%%xmm1			\n\t	movaps	%%xmm7	,%%xmm5			\n\t"\
		"mulpd	 0xc0(%%edi)	,%%xmm3	\n\t	mulpd	 0xc0(%%ebx)	,%%xmm7	\n\t"\
		"addpd	%%xmm4	,%%xmm3			\n\t	addpd	%%xmm4	,%%xmm7			\n\t"\
		"subpd	%%xmm4	,%%xmm3			\n\t	subpd	%%xmm4	,%%xmm7			\n\t"\
		"movaps	%%xmm3	,(%%ecx)		\n\t	movaps	%%xmm7	,(%%edx)		\n\t"\
		"\n\t"\
		"\n\t"\
		"movl	%[__data]	,%%eax		\n\t"\
		"mulpd	 0x80(%%edi)	,%%xmm3	\n\t	mulpd	 0x80(%%ebx)	,%%xmm7	\n\t"\
		"subpd	%%xmm3	,%%xmm1			\n\t	subpd	%%xmm7	,%%xmm5			\n\t"\
		"mulpd	%%xmm2	,%%xmm1			\n\t	mulpd	%%xmm6	,%%xmm5			\n\t"\
		"movaps	%%xmm1	, 0x10(%%eax)	\n\t	movaps	%%xmm5	, 0x50(%%eax)	\n\t"\
		"\n\t"\
		"\n\t"\
		"movl	%[__sse_bw]	,%%eax		\n\t"\
		"movl	%[__sse_nm1]	,%%ebx	\n\t"\
		"paddd	(%%eax)	,%%xmm0			\n\t"\
		"pand	(%%ebx)	,%%xmm0			\n\t"\
		"movl	%[__bjmod_0],%%ecx		\n\t"\
		"movaps	%%xmm0,(%%ecx)			\n\t"\
	"popl %%ebx	\n\t"\
			:					/* outputs: none */\
		: [__data]		"m" (Xdata)	/* All inputs from memory addresses here */\
		, [__wtA]		"m" (XwtA)		\
		, [__wtB]		"m" (XwtB)		\
		, [__wtC]		"m" (XwtC)		\
		, [__cyA]		"m" (XcyA)		\
		, [__cyB]		"m" (XcyB)		\
		, [__bjmod_0]	"m" (Xbjmod_0)		\
		, [__half_arr]	"m" (Xhalf_arr)		\
		, [__n_minus_silp1] "m" (Xn_minus_silp1)\
		, [__n_minus_sil]	"m" (Xn_minus_sil)	\
		, [__sinwt]		"m" (Xsinwt)		\
		, [__sinwtm1]		"m" (Xsinwtm1)		\
		, [__sse_bw]		"m" (Xsse_bw)		\
		, [__sse_nm1]		"m" (Xsse_nm1)		\
		, [__sse_sw]		"m" (Xsse_sw)		\
		: "cc","memory","eax","ecx","edx","edi","esi"		/* Clobbered registers */\
	);\
	}

	#define SSE2_cmplx_carry_norm_pow2_nocheck2_2B(Xdata,XwtA,XwtB,XcyA,XcyB,Xbjmod_0,Xhalf_arr,Xn_minus_silp1,Xn_minus_sil,Xsinwt,Xsinwtm1,Xsse_bw,Xsse_nm1,Xsse_sw)\
	{\
	__asm__ volatile (\
	"pushl %%ebx	\n\t"/* Explicit save/restore of PIC register */\
	"/**********************************************/\n\t"\
	"/*          Real      parts                   */\n\t"\
	"/**********************************************/\n\t"\
		"movl	%[__bjmod_0],	%%eax	\n\t"\
		"movaps		(%%eax)	,	%%xmm0	\n\t"\
		"movl	%[__sse_sw]	,	%%ebx	\n\t"\
		"movaps		(%%ebx)	,	%%xmm1	\n\t"\
		"psubd		%%xmm0	,	%%xmm1	\n\t"\
		"movmskps	%%xmm1	,	%%esi	\n\t"\
		"\n\t"\
		"shll	$24		,%%esi			\n\t"\
		"movaps		%%xmm0	,%%xmm1		\n\t"\
		"movl	%[__n_minus_sil],%%ecx	\n\t"\
		"movd	%%ecx	,%%xmm2			\n\t"\
		"pshufd	$0,	%%xmm2	,%%xmm2		\n\t"\
		"psubd		%%xmm0	,%%xmm2		\n\t"\
		"movmskps	%%xmm2	,%%ecx		\n\t"\
		"shll	$16		,%%ecx			\n\t"\
		"addl	%%ecx	,%%esi			\n\t"\
		"movl	%[__sinwt]	,%%edx		\n\t"\
		"movd	%%edx	,%%xmm3			\n\t"\
		"pshufd	$0,	%%xmm3	,%%xmm3		\n\t"\
		"psubd		%%xmm3	,%%xmm1		\n\t"\
		"movmskps	%%xmm1	,%%edx		\n\t"\
		"shll	$8		,%%edx			\n\t"\
		"addl	%%edx	,%%esi			\n\t"\
		"\n\t"\
		"mov	%[__data],%%eax			\n\t"\
		"movaps	 0x20(%%eax),%%xmm1		\n\t	movaps		 0x60(%%eax),%%xmm5\n\t"\
		"\n\t"\
		"movl	%[__half_arr]	,%%eax	\n\t"\
		"movl	%[__wtA]	,%%ebx		\n\t"\
		"movl	%[__wtB]	,%%ecx		\n\t"\
		"movaps		-0x10(%%eax),%%xmm4	\n\t"\
		"movaps		     (%%ebx),%%xmm2	\n\t	movaps		 0x10(%%ebx),%%xmm6	\n\t"\
		"movhpd		     (%%ecx),%%xmm3	\n\t	movhpd		-0x10(%%ecx),%%xmm7	\n\t"\
		"movlpd		 0x08(%%ecx),%%xmm3	\n\t	movlpd		-0x08(%%ecx),%%xmm7	\n\t"\
		"\n\t"\
		"shld	$20,	%%esi	,%%edi	\n\t	shld	$18,	%%esi	,%%ebx	\n\t"\
		"andl	$0x00000030	,%%edi		\n\t	andl	$0x00000030	,%%ebx		\n\t"\
		"shld	$28,	%%esi,	%%edx	\n\t	shld	$26,	%%esi,	%%ecx	\n\t"\
		"andl	$0x00000030	,%%edx		\n\t	andl	$0x00000030	,%%ecx		\n\t"\
		"addl	%%eax	,%%edi			\n\t	addl	%%eax	,%%ebx			\n\t"\
		"addl	%%eax	,%%edx			\n\t	addl	%%eax	,%%ecx			\n\t"\
		"\n\t"\
		"mulpd	 	0x100(%%eax),%%xmm2	\n\t	mulpd	 	0x100(%%eax),%%xmm6	\n\t"\
		"mulpd	 	0x110(%%eax),%%xmm3	\n\t	mulpd	 	0x110(%%eax),%%xmm7	\n\t"\
		"mulpd	 	     (%%edi),%%xmm2	\n\t	mulpd	 	     (%%ebx),%%xmm6	\n\t"\
		"mulpd	 	0x040(%%edx),%%xmm3	\n\t	mulpd	 	0x040(%%ecx),%%xmm7	\n\t"\
		"\n\t"\
		"movl	%[__cyA]	,%%ecx		\n\t	movl	%[__cyB]	,%%edx		\n\t"\
		"mulpd		%%xmm3	,%%xmm1		\n\t	mulpd		%%xmm7	,%%xmm5		\n\t"\
		"addpd		(%%ecx)	,%%xmm1		\n\t	addpd		(%%edx)	,%%xmm5		\n\t"\
		"movaps		%%xmm1	,%%xmm3		\n\t	movaps		%%xmm5	,%%xmm7		\n\t"\
		"addpd		%%xmm4	,%%xmm3		\n\t	addpd		%%xmm4	,%%xmm7		\n\t"\
		"subpd		%%xmm4	,%%xmm3		\n\t	subpd		%%xmm4	,%%xmm7		\n\t"\
		"\n\t"\
		"/* NO ROE HERE */		\n\t"\
		"\n\t"\
		"shld	$12		,%%esi	,%%edi	\n\t	shld	$10		,%%esi	,%%ebx	\n\t"\
		"andl	$0x00000030	,%%edi		\n\t	andl	$0x00000030	,%%ebx		\n\t"\
		"addl	%%eax	,%%edi			\n\t	addl	%%eax	,%%ebx			\n\t"\
		"movaps		%%xmm3	,%%xmm1		\n\t	movaps		%%xmm7	,%%xmm5		\n\t"\
		"mulpd	 	0xc0(%%edi)	,%%xmm3	\n\t	mulpd	 	0xc0(%%ebx)	,%%xmm7	\n\t"\
		"addpd		%%xmm4	,%%xmm3		\n\t	addpd		%%xmm4	,%%xmm7		\n\t"\
		"subpd		%%xmm4	,%%xmm3		\n\t	subpd		%%xmm4	,%%xmm7		\n\t"\
		"movaps		%%xmm3	,(%%ecx)	\n\t	movaps		%%xmm7	,(%%edx)	\n\t"\
		"\n\t"\
		"\n\t"\
		"movl	%[__data]	,%%eax		\n\t"\
		"mulpd	 0x80(%%edi)	,%%xmm3	\n\t	mulpd	 0x80(%%ebx)	,%%xmm7	\n\t"\
		"subpd		%%xmm3	,	%%xmm1	\n\t	subpd		%%xmm7	,	%%xmm5	\n\t"\
		"mulpd		%%xmm2	,	%%xmm1	\n\t	mulpd		%%xmm6	,	%%xmm5	\n\t"\
		"movaps		%%xmm1	,0x20(%%eax)\n\t	movaps		%%xmm5	,0x60(%%eax)\n\t"\
		"\n\t"\
		"\n\t"\
		"movl	%[__sse_bw]	,%%eax		\n\t"\
		"movl	%[__sse_nm1]	,%%ebx	\n\t"\
		"paddd		(%%eax)	,%%xmm0		\n\t"\
		"pand		(%%ebx)	,%%xmm0		\n\t"\
	"/**********************************************/\n\t"\
	"/*          Imaginary parts               */\n\t"\
	"/**********************************************/\n\t"\
		"movl	%[__sse_sw]	,%%edx		\n\t"\
		"movaps	(%%edx)	,%%xmm1			\n\t"\
		"psubd	%%xmm0	,%%xmm1			\n\t"\
		"movmskps	%%xmm1	,%%esi		\n\t"\
		"shll	$24	,%%esi				\n\t"\
		"movaps	%%xmm0	,%%xmm1			\n\t"\
		"movl	%[__n_minus_silp1],%%ecx\n\t"\
		"movd	%%ecx	,%%xmm2			\n\t"\
		"pshufd	$0	,%%xmm2	,%%xmm2		\n\t"\
		"psubd	%%xmm0	,%%xmm2			\n\t"\
		"movmskps	%%xmm2	,%%ecx		\n\t"\
		"shll	$16	,%%ecx				\n\t"\
		"addl	%%ecx	,%%esi			\n\t"\
		"movl	%[__sinwtm1]	,%%edx	\n\t"\
		"movd	%%edx	,%%xmm3			\n\t"\
		"pshufd	$0	,%%xmm3	,%%xmm3		\n\t"\
		"psubd	%%xmm3	,%%xmm1			\n\t"\
		"movmskps	%%xmm1	,%%edx		\n\t"\
		"shll	$8	,%%edx				\n\t"\
		"addl	%%edx	,%%esi			\n\t"\
		"movl	%[__data]	,%%eax		\n\t"\
		"movaps	 0x30(%%eax)	,%%xmm1	\n\t"\
		"movaps	 0x70(%%eax)	,%%xmm5	\n\t"\
		"\n\t"\
		"movl	%[__half_arr]	,%%eax	\n\t"\
		"movl	%[__wtA]	,%%ebx		\n\t"\
		"movl	%[__wtB]	,%%ecx		\n\t"\
		"\n\t"\
		"movaps	     (%%ebx)	,%%xmm2	\n\t	movaps	 0x10(%%ebx)	,%%xmm6	\n\t"\
		"movhpd	     (%%ecx)	,%%xmm3	\n\t	movhpd	-0x10(%%ecx)	,%%xmm7	\n\t"\
		"movlpd	 0x08(%%ecx)	,%%xmm3	\n\t	movlpd	-0x08(%%ecx)	,%%xmm7	\n\t"\
		"\n\t"\
		"addl	$0x20	,%%ebx			\n\t"\
		"subl	$0x20	,%%ecx			\n\t"\
		"movl	%%ebx	,%[__wtA]		\n\t"\
		"movl	%%ecx	,%[__wtB]		\n\t"\
		"shld	$20		,%%esi	,%%edi	\n\t	shld	$18		,%%esi	,%%ebx	\n\t"\
		"andl	$0x00000030	,%%edi		\n\t	andl	$0x00000030	,%%ebx		\n\t"\
		"shld	$28,	%%esi,	%%edx	\n\t	shld	$26,	%%esi,	%%ecx	\n\t"\
		"andl	$0x00000030	,%%edx		\n\t	andl	$0x00000030	,%%ecx		\n\t"\
		"addl	%%eax	,%%edi			\n\t	addl	%%eax	,%%ebx			\n\t"\
		"addl	%%eax	,%%edx			\n\t	addl	%%eax	,%%ecx			\n\t"\
		"\n\t"\
		"mulpd	 0x120(%%eax)	,%%xmm2	\n\t	mulpd	 0x120(%%eax)	,%%xmm6	\n\t"\
		"mulpd	 0x130(%%eax)	,%%xmm3	\n\t	mulpd	 0x130(%%eax)	,%%xmm7	\n\t"\
		"mulpd	      (%%edi)	,%%xmm2	\n\t	mulpd	      (%%ebx)	,%%xmm6	\n\t"\
		"mulpd	 0x040(%%edx)	,%%xmm3	\n\t	mulpd	 0x040(%%ecx)	,%%xmm7	\n\t"\
		"\n\t"\
		"movl	%[__cyA]	,%%ecx		\n\t	movl	%[__cyB]	,%%edx		\n\t"\
		"mulpd	%%xmm3	,%%xmm1			\n\t	mulpd	%%xmm7	,%%xmm5			\n\t"\
		"addpd	(%%ecx)	,%%xmm1			\n\t	addpd	(%%edx)	,%%xmm5			\n\t"\
		"movaps	%%xmm1	,%%xmm3			\n\t	movaps	%%xmm5	,%%xmm7			\n\t"\
		"addpd	%%xmm4	,%%xmm3			\n\t	addpd	%%xmm4	,%%xmm7			\n\t"\
		"subpd	%%xmm4	,%%xmm3			\n\t	subpd	%%xmm4	,%%xmm7			\n\t"\
		"\n\t"\
		"/* NO ROE HERE */		\n\t"\
		"\n\t"\
		"shld	$12		,%%esi	,%%edi	\n\t	shld	$10		,%%esi	,%%ebx	\n\t"\
		"andl	$0x00000030	,%%edi		\n\t	andl	$0x00000030	,%%ebx		\n\t"\
		"addl	%%eax	,%%edi			\n\t	addl	%%eax	,%%ebx			\n\t"\
		"movaps	%%xmm3	,%%xmm1			\n\t	movaps	%%xmm7	,%%xmm5			\n\t"\
		"mulpd	 0xc0(%%edi)	,%%xmm3	\n\t	mulpd	 0xc0(%%ebx)	,%%xmm7	\n\t"\
		"addpd	%%xmm4	,%%xmm3			\n\t	addpd	%%xmm4	,%%xmm7			\n\t"\
		"subpd	%%xmm4	,%%xmm3			\n\t	subpd	%%xmm4	,%%xmm7			\n\t"\
		"movaps	%%xmm3	,(%%ecx)		\n\t	movaps	%%xmm7	,(%%edx)		\n\t"\
		"\n\t"\
		"\n\t"\
		"movl	%[__data]	,%%eax		\n\t"\
		"mulpd	 0x80(%%edi)	,%%xmm3	\n\t	mulpd	 0x80(%%ebx)	,%%xmm7	\n\t"\
		"subpd	%%xmm3	,%%xmm1			\n\t	subpd	%%xmm7	,%%xmm5			\n\t"\
		"mulpd	%%xmm2	,%%xmm1			\n\t	mulpd	%%xmm6	,%%xmm5			\n\t"\
		"movaps	%%xmm1	, 0x30(%%eax)	\n\t	movaps	%%xmm5	, 0x70(%%eax)	\n\t"\
		"\n\t"\
		"\n\t"\
		"movl	%[__sse_bw]	,%%eax		\n\t"\
		"movl	%[__sse_nm1]	,%%ebx	\n\t"\
		"paddd	(%%eax)	,%%xmm0			\n\t"\
		"pand	(%%ebx)	,%%xmm0			\n\t"\
		"movl	%[__bjmod_0],%%ecx		\n\t"\
		"movaps	%%xmm0,(%%ecx)			\n\t"\
	"/**********************************************/\n\t"\
	"/*          Repack the data:              */\n\t"\
	"/**********************************************/\n\t"\
		"mov	%[__data],%%eax			\n\t"\
		"movaps		0x10(%%eax)	,%%xmm1	\n\t	movaps		0x50(%%eax)	,%%xmm5	\n\t"\
		"movaps		    (%%eax)	,%%xmm0	\n\t	movaps		0x40(%%eax)	,%%xmm4	\n\t"\
		"movaps		%%xmm1		,%%xmm3	\n\t	movaps		%%xmm5		,%%xmm7	\n\t"\
		"movaps		%%xmm0		,%%xmm2	\n\t	movaps		%%xmm4		,%%xmm6	\n\t"\
		"unpckhpd	0x30(%%eax)	,%%xmm3	\n\t	unpckhpd	0x70(%%eax)	,%%xmm7	\n\t"\
		"unpcklpd	0x30(%%eax)	,%%xmm1	\n\t	unpcklpd	0x70(%%eax)	,%%xmm5	\n\t"\
		"movaps		%%xmm3,0x30(%%eax)	\n\t	movaps		%%xmm7,0x70(%%eax)	\n\t"\
		"unpckhpd	0x20(%%eax)	,%%xmm2	\n\t	unpckhpd	0x60(%%eax)	,%%xmm6	\n\t"\
		"unpcklpd	0x20(%%eax)	,%%xmm0	\n\t	unpcklpd	0x60(%%eax)	,%%xmm4	\n\t"\
		"movaps		%%xmm2,0x20(%%eax)	\n\t	movaps		%%xmm6,0x60(%%eax)	\n\t"\
		"movaps		%%xmm1,0x10(%%eax)	\n\t	movaps		%%xmm5,0x50(%%eax)	\n\t"\
		"movaps		%%xmm0,    (%%eax)	\n\t	movaps		%%xmm4,0x40(%%eax)	\n\t"\
	"popl %%ebx	\n\t"\
		:					/* outputs: none */\
		: [__data]		"m" (Xdata)	/* All inputs from memory addresses here */\
		, [__wtA]		"m" (XwtA)		\
		, [__wtB]		"m" (XwtB)		\
		, [__cyA]		"m" (XcyA)		\
		, [__cyB]		"m" (XcyB)		\
		, [__bjmod_0]	"m" (Xbjmod_0)		\
		, [__half_arr]	"m" (Xhalf_arr)		\
		, [__n_minus_silp1] "m" (Xn_minus_silp1)\
		, [__n_minus_sil]	"m" (Xn_minus_sil)	\
		, [__sinwt]		"m" (Xsinwt)		\
		, [__sinwtm1]		"m" (Xsinwtm1)		\
		, [__sse_bw]		"m" (Xsse_bw)		\
		, [__sse_nm1]		"m" (Xsse_nm1)		\
		, [__sse_sw]		"m" (Xsse_sw)		\
		: "cc","memory","eax","ecx","edx","edi","esi"		/* Clobbered registers */\
	);\
	}


	/***************************************************************************************************************************************************/
	/********* Non-power-of-2-FFT versions of SSE2_cmplx_carry_norm_pow2_errcheck0_2B,1_2B,2_2B (only give sans-error-check version of latter 2: *******/
	/***************************************************************************************************************************************************/

	#define SSE2_cmplx_carry_norm_errcheck0_2B(Xdata,XwtA,XwtB,XwtC,XcyA,XcyB,Xbjmod_0,Xhalf_arr,Xi,Xn_minus_silp1,Xn_minus_sil,Xsign_mask,Xsinwt,Xsinwtm1,Xsse_bw,Xsse_n,Xsse_sw)\
	{\
	__asm__ volatile (\
	"pushl %%ebx	\n\t"/* Explicit save/restore of PIC register */\
	"/***************Unpack the data:*************************/\n\t"\
		"movl	%[__data]	,%%eax	\n\t"\
		"movaps		    (%%eax)	,%%xmm1	\n\t	movaps		0x40(%%eax)	,%%xmm5	\n\t"\
		"movaps		    (%%eax)	,%%xmm2	\n\t	movaps		0x40(%%eax)	,%%xmm6	\n\t"\
		"unpcklpd	0x20(%%eax)	,%%xmm1	\n\t	unpcklpd	0x60(%%eax)	,%%xmm5	\n\t"\
		"unpckhpd	0x20(%%eax)	,%%xmm2	\n\t	unpckhpd	0x60(%%eax)	,%%xmm6	\n\t"\
		"movaps		%%xmm2, 0x20(%%eax)	\n\t	movaps		%%xmm6, 0x60(%%eax)	\n\t"\
		"\n\t"\
		"movaps		0x10(%%eax)	,%%xmm2	\n\t	movaps		0x50(%%eax)	,%%xmm6	\n\t"\
		"movaps		0x10(%%eax)	,%%xmm3	\n\t	movaps		0x50(%%eax)	,%%xmm7	\n\t"\
		"unpcklpd	0x30(%%eax)	,%%xmm2	\n\t	unpcklpd	0x70(%%eax)	,%%xmm6	\n\t"\
		"unpckhpd	0x30(%%eax)	,%%xmm3	\n\t	unpckhpd	0x70(%%eax)	,%%xmm7	\n\t"\
		"movaps		%%xmm2, 0x10(%%eax)	\n\t	movaps		%%xmm6, 0x50(%%eax)	\n\t"\
		"movaps		%%xmm3, 0x30(%%eax)	\n\t	movaps		%%xmm7, 0x70(%%eax)	\n\t"\
	"/**********************************************/\n\t"\
	"/*          Real      parts                   */\n\t"\
	"/**********************************************/\n\t"\
		"movl	%[__bjmod_0],	%%eax	\n\t"\
		"movaps		(%%eax)	,	%%xmm0	\n\t"\
		"movl	%[__sse_sw]	,	%%ebx	\n\t"\
		"movaps		(%%ebx)	,	%%xmm7	\n\t"\
		"psubd		%%xmm0	,	%%xmm7	\n\t"\
		"movmskps	%%xmm7	,	%%esi	\n\t"\
		"movl	%[__i]	,%%ecx			\n\t"\
		"andl	$0xfffffffe,%%esi		\n\t"\
		"addl	%%ecx	,%%esi			\n\t"\
		"shll	$24		,%%esi			\n\t"\
		"movaps		%%xmm0	,%%xmm7		\n\t"\
		"movl	%[__n_minus_sil],%%ecx	\n\t"\
		"movd	%%ecx	,%%xmm2			\n\t"\
		"pshufd	$0,	%%xmm2	,%%xmm2		\n\t"\
		"psubd		%%xmm0	,%%xmm2		\n\t"\
		"movmskps	%%xmm2	,%%ecx		\n\t"\
		"shll	$16		,%%ecx			\n\t"\
		"addl	%%ecx	,%%esi			\n\t"\
		"movl	%[__sinwt]	,%%edx		\n\t"\
		"movd	%%edx	,%%xmm3			\n\t"\
		"pshufd	$0,	%%xmm3	,%%xmm3		\n\t"\
		"psubd		%%xmm3	,%%xmm7		\n\t"\
		"movmskps	%%xmm7	,%%edx		\n\t"\
		"shll	$8		,%%edx			\n\t"\
		"addl	%%edx	,%%esi			\n\t"\
		"\n\t"\
		"movl	%[__half_arr]	,%%eax	\n\t"\
		"movl	%[__wtA]	,%%ebx		\n\t"\
		"movl	%[__wtB]	,%%ecx		\n\t"\
		"movaps		-0x10(%%eax),%%xmm4	\n\t"\
		"movaps		     (%%ebx),%%xmm2	\n\t	movaps		 0x10(%%ebx),%%xmm6	\n\t"\
		"movhpd		     (%%ecx),%%xmm3	\n\t	movhpd		-0x10(%%ecx),%%xmm7	\n\t"\
		"movlpd		 0x08(%%ecx),%%xmm3	\n\t	movlpd		-0x08(%%ecx),%%xmm7	\n\t"\
		"\n\t"\
		"subl	$0x20	,%%ecx			\n\t"\
		"movl	%%ecx	,%[__wtB]		\n\t"\
		"\n\t"\
		"shld	$20,	%%esi	,%%edi	\n\t	shld	$18,	%%esi	,%%ebx	\n\t"\
		"andl	$0x00000030	,%%edi		\n\t	andl	$0x00000030	,%%ebx		\n\t"\
		"shld	$28,	%%esi,	%%edx	\n\t	shld	$26,	%%esi,	%%ecx	\n\t"\
		"andl	$0x00000030	,%%edx		\n\t	andl	$0x00000030	,%%ecx		\n\t"\
		"addl	%%eax	,%%edi			\n\t	addl	%%eax	,%%ebx			\n\t"\
		"addl	%%eax	,%%edx			\n\t	addl	%%eax	,%%ecx			\n\t"\
		"\n\t"\
		"mulpd	 	0x100(%%eax),%%xmm2	\n\t	mulpd	 	0x100(%%eax),%%xmm6	\n\t"\
		"mulpd	 	0x110(%%eax),%%xmm3	\n\t	mulpd	 	0x110(%%eax),%%xmm7	\n\t"\
		"mulpd	 	     (%%edi),%%xmm2	\n\t	mulpd	 	     (%%ebx),%%xmm6	\n\t"\
		"mulpd	 	0x040(%%edx),%%xmm3	\n\t	mulpd	 	0x040(%%ecx),%%xmm7	\n\t"\
		"\n\t"\
		"movl	%[__cyA]	,%%ecx		\n\t	movl	%[__cyB]	,%%edx		\n\t"\
		"mulpd		%%xmm3	,%%xmm1		\n\t	mulpd		%%xmm7	,%%xmm5		\n\t"\
		"addpd		(%%ecx)	,%%xmm1		\n\t	addpd		(%%edx)	,%%xmm5		\n\t"\
		"movaps		%%xmm1	,%%xmm3		\n\t	movaps		%%xmm5	,%%xmm7		\n\t"\
		"addpd		%%xmm4	,%%xmm3		\n\t	addpd		%%xmm4	,%%xmm7		\n\t"\
		"subpd		%%xmm4	,%%xmm3		\n\t	subpd		%%xmm4	,%%xmm7		\n\t"\
		"\n\t"\
		"movl	%[__sign_mask],%%ebx	\n\t"\
		"subpd			%%xmm3	,%%xmm1	\n\t	subpd			%%xmm7	,%%xmm5	\n\t"\
		"andpd			(%%ebx)	,%%xmm1	\n\t	andpd			(%%ebx)	,%%xmm5	\n\t"\
		"maxpd			%%xmm5	,%%xmm1	\n\t"\
		"maxpd		-0x20(%%eax),%%xmm1	\n\t"\
		"movaps		%%xmm1,-0x20(%%eax)	\n\t"\
		"\n\t"\
		"shld	$12		,%%esi	,%%edi	\n\t	shld	$10		,%%esi	,%%ebx	\n\t"\
		"andl	$0x00000030	,%%edi		\n\t	andl	$0x00000030	,%%ebx		\n\t"\
		"addl	%%eax	,%%edi			\n\t	addl	%%eax	,%%ebx			\n\t"\
		"movaps		%%xmm3	,%%xmm1		\n\t	movaps		%%xmm7	,%%xmm5		\n\t"\
		"mulpd	 	0xc0(%%edi)	,%%xmm3	\n\t	mulpd	 	0xc0(%%ebx)	,%%xmm7	\n\t"\
		"addpd		%%xmm4	,%%xmm3		\n\t	addpd		%%xmm4	,%%xmm7		\n\t"\
		"subpd		%%xmm4	,%%xmm3		\n\t	subpd		%%xmm4	,%%xmm7		\n\t"\
		"movaps		%%xmm3	,(%%ecx)	\n\t	movaps		%%xmm7	,(%%edx)	\n\t"\
		"\n\t"\
		"\n\t"\
		"movl	%[__data]	,%%eax		\n\t"\
		"mulpd	 0x80(%%edi)	,%%xmm3	\n\t	mulpd	 0x80(%%ebx)	,%%xmm7	\n\t"\
		"subpd		%%xmm3	,	%%xmm1	\n\t	subpd		%%xmm7	,	%%xmm5	\n\t"\
		"mulpd		%%xmm2	,	%%xmm1	\n\t	mulpd		%%xmm6	,	%%xmm5	\n\t"\
		"movaps		%%xmm1	,    (%%eax)\n\t	movaps		%%xmm5	,0x40(%%eax)\n\t"\
		"\n\t"\
		"movl	%[__sse_n]	,%%ebx		\n\t"\
		"movaps		(%%ebx)	,%%xmm2		\n\t"\
		"movl	%[__sse_bw]	,%%eax		\n\t"\
		"paddd		(%%eax)	,%%xmm0		\n\t"\
		"movaps		%%xmm0	,%%xmm1		\n\t"\
		"pcmpgtd	%%xmm2	,%%xmm1		\n\t"\
		"pand		(%%ebx)	,%%xmm1		\n\t"\
		"psubd		%%xmm1	,%%xmm0		\n\t"\
	"/**********************************************/\n\t"\
	"/*          Imaginary parts               */\n\t"\
	"/**********************************************/\n\t"\
		"movl	%[__sse_sw]	,%%edx		\n\t"\
		"movaps	(%%edx)	,%%xmm1			\n\t"\
		"psubd	%%xmm0	,%%xmm1			\n\t"\
		"movmskps	%%xmm1	,%%esi		\n\t"\
		"shll	$24	,%%esi				\n\t"\
		"movaps	%%xmm0	,%%xmm1			\n\t"\
		"movl	%[__n_minus_silp1],%%ecx\n\t"\
		"movd	%%ecx	,%%xmm2			\n\t"\
		"pshufd	$0	,%%xmm2	,%%xmm2		\n\t"\
		"psubd	%%xmm0	,%%xmm2			\n\t"\
		"movmskps	%%xmm2	,%%ecx		\n\t"\
		"shll	$16	,%%ecx				\n\t"\
		"addl	%%ecx	,%%esi			\n\t"\
		"movl	%[__sinwtm1]	,%%edx	\n\t"\
		"movd	%%edx	,%%xmm3			\n\t"\
		"pshufd	$0	,%%xmm3	,%%xmm3		\n\t"\
		"psubd	%%xmm3	,%%xmm1			\n\t"\
		"movmskps	%%xmm1	,%%edx		\n\t"\
		"shll	$8	,%%edx				\n\t"\
		"addl	%%edx	,%%esi			\n\t"\
		"movl	%[__data]	,%%eax		\n\t"\
		"movaps	 0x10(%%eax)	,%%xmm1	\n\t"\
		"movaps	 0x50(%%eax)	,%%xmm5	\n\t"\
		"\n\t"\
		"movl	%[__half_arr]	,%%eax	\n\t"\
		"movl	%[__wtA]	,%%ebx		\n\t"\
		"movl	%[__wtC]	,%%ecx		\n\t"\
		"\n\t"\
		"movaps	     (%%ebx)	,%%xmm2	\n\t	movaps	 0x10(%%ebx)	,%%xmm6	\n\t"\
		"movhpd	     (%%ecx)	,%%xmm3	\n\t	movhpd	-0x10(%%ecx)	,%%xmm7	\n\t"\
		"movlpd	 0x08(%%ecx)	,%%xmm3	\n\t	movlpd	-0x08(%%ecx)	,%%xmm7	\n\t"\
		"\n\t"\
		"addl	$0x20	,%%ebx			\n\t"\
		"subl	$0x20	,%%ecx			\n\t"\
		"movl	%%ebx	,%[__wtA]		\n\t"\
		"movl	%%ecx	,%[__wtC]		\n\t"\
		"shld	$20		,%%esi	,%%edi	\n\t	shld	$18		,%%esi	,%%ebx	\n\t"\
		"andl	$0x00000030	,%%edi		\n\t	andl	$0x00000030	,%%ebx		\n\t"\
		"shld	$28,	%%esi,	%%edx	\n\t	shld	$26,	%%esi,	%%ecx	\n\t"\
		"andl	$0x00000030	,%%edx		\n\t	andl	$0x00000030	,%%ecx		\n\t"\
		"addl	%%eax	,%%edi			\n\t	addl	%%eax	,%%ebx			\n\t"\
		"addl	%%eax	,%%edx			\n\t	addl	%%eax	,%%ecx			\n\t"\
		"\n\t"\
		"mulpd	 0x120(%%eax)	,%%xmm2	\n\t	mulpd	 0x120(%%eax)	,%%xmm6	\n\t"\
		"mulpd	 0x130(%%eax)	,%%xmm3	\n\t	mulpd	 0x130(%%eax)	,%%xmm7	\n\t"\
		"mulpd	      (%%edi)	,%%xmm2	\n\t	mulpd	      (%%ebx)	,%%xmm6	\n\t"\
		"mulpd	 0x040(%%edx)	,%%xmm3	\n\t	mulpd	 0x040(%%ecx)	,%%xmm7	\n\t"\
		"\n\t"\
		"movl	%[__cyA]	,%%ecx		\n\t	movl	%[__cyB]	,%%edx		\n\t"\
		"mulpd	%%xmm3	,%%xmm1			\n\t	mulpd	%%xmm7	,%%xmm5			\n\t"\
		"addpd	(%%ecx)	,%%xmm1			\n\t	addpd	(%%edx)	,%%xmm5			\n\t"\
		"movaps	%%xmm1	,%%xmm3			\n\t	movaps	%%xmm5	,%%xmm7			\n\t"\
		"addpd	%%xmm4	,%%xmm3			\n\t	addpd	%%xmm4	,%%xmm7			\n\t"\
		"subpd	%%xmm4	,%%xmm3			\n\t	subpd	%%xmm4	,%%xmm7			\n\t"\
		"\n\t"\
		"movl	%[__sign_mask],%%ebx	\n\t"\
		"subpd	%%xmm3	,%%xmm1			\n\t	subpd	%%xmm7	,%%xmm5			\n\t"\
		"andpd	(%%ebx)	,%%xmm1			\n\t	andpd	(%%ebx)	,%%xmm5			\n\t"\
		"maxpd	%%xmm5	,%%xmm1			\n\t"\
		"maxpd	-0x20(%%eax)	,%%xmm1	\n\t"\
		"movaps	%%xmm1	,-0x20(%%eax)	\n\t"\
		"\n\t"\
		"shld	$12		,%%esi	,%%edi	\n\t	shld	$10		,%%esi	,%%ebx	\n\t"\
		"andl	$0x00000030	,%%edi		\n\t	andl	$0x00000030	,%%ebx		\n\t"\
		"addl	%%eax	,%%edi			\n\t	addl	%%eax	,%%ebx			\n\t"\
		"movaps	%%xmm3	,%%xmm1			\n\t	movaps	%%xmm7	,%%xmm5			\n\t"\
		"mulpd	 0xc0(%%edi)	,%%xmm3	\n\t	mulpd	 0xc0(%%ebx)	,%%xmm7	\n\t"\
		"addpd	%%xmm4	,%%xmm3			\n\t	addpd	%%xmm4	,%%xmm7			\n\t"\
		"subpd	%%xmm4	,%%xmm3			\n\t	subpd	%%xmm4	,%%xmm7			\n\t"\
		"movaps	%%xmm3	,(%%ecx)		\n\t	movaps	%%xmm7	,(%%edx)		\n\t"\
		"\n\t"\
		"\n\t"\
		"movl	%[__data]	,%%eax		\n\t"\
		"mulpd	 0x80(%%edi)	,%%xmm3	\n\t	mulpd	 0x80(%%ebx)	,%%xmm7	\n\t"\
		"subpd	%%xmm3	,%%xmm1			\n\t	subpd	%%xmm7	,%%xmm5			\n\t"\
		"mulpd	%%xmm2	,%%xmm1			\n\t	mulpd	%%xmm6	,%%xmm5			\n\t"\
		"movaps	%%xmm1	, 0x10(%%eax)	\n\t	movaps	%%xmm5	, 0x50(%%eax)	\n\t"\
		"\n\t"\
		"\n\t"\
		"movl	%[__sse_n]	,%%ebx		\n\t"\
		"movaps		(%%ebx)	,%%xmm2		\n\t"\
		"movl	%[__sse_bw]	,%%eax		\n\t"\
		"paddd		(%%eax)	,%%xmm0		\n\t"\
		"movaps		%%xmm0	,%%xmm1		\n\t"\
		"pcmpgtd	%%xmm2	,%%xmm1		\n\t"\
		"pand		(%%ebx)	,%%xmm1		\n\t"\
		"psubd		%%xmm1	,%%xmm0		\n\t"\
		"movl	%[__bjmod_0],%%ecx		\n\t"\
		"movaps	%%xmm0,(%%ecx)			\n\t"\
	"popl %%ebx	\n\t"\
		:					/* outputs: none */\
		: [__data]		"m" (Xdata)	/* All inputs from memory addresses here */\
		, [__wtA]		"m" (XwtA)		\
		, [__wtB]		"m" (XwtB)		\
		, [__wtC]		"m" (XwtC)		\
		, [__cyA]		"m" (XcyA)		\
		, [__cyB]		"m" (XcyB)		\
		, [__bjmod_0]	"m" (Xbjmod_0)		\
		, [__half_arr]	"m" (Xhalf_arr)		\
		, [__i]			"m" (Xi)			\
		, [__n_minus_silp1] "m" (Xn_minus_silp1)\
		, [__n_minus_sil]	"m" (Xn_minus_sil)	\
		, [__sign_mask]	"m" (Xsign_mask)	\
		, [__sinwt]		"m" (Xsinwt)		\
		, [__sinwtm1]		"m" (Xsinwtm1)		\
		, [__sse_bw]		"m" (Xsse_bw)		\
		, [__sse_n]		"m" (Xsse_n)		\
		, [__sse_sw]		"m" (Xsse_sw)		\
		: "cc","memory","eax","ecx","edx","edi","esi"		/* Clobbered registers */\
	);\
	}

	#define SSE2_cmplx_carry_norm_errcheck1_2B(Xdata,XwtA,XwtB,XwtC,XcyA,XcyB,Xbjmod_0,Xhalf_arr,Xn_minus_silp1,Xn_minus_sil,Xsign_mask,Xsinwt,Xsinwtm1,Xsse_bw,Xsse_n,Xsse_sw)\
	{\
	__asm__ volatile (\
	"pushl %%ebx	\n\t"/* Explicit save/restore of PIC register */\
	"/***************Unpack the data:*************************/\n\t"\
		"movl	%[__data]	,%%eax	\n\t"\
		"movaps		    (%%eax)	,%%xmm1	\n\t	movaps		0x40(%%eax)	,%%xmm5	\n\t"\
		"movaps		    (%%eax)	,%%xmm2	\n\t	movaps		0x40(%%eax)	,%%xmm6	\n\t"\
		"unpcklpd	0x20(%%eax)	,%%xmm1	\n\t	unpcklpd	0x60(%%eax)	,%%xmm5	\n\t"\
		"unpckhpd	0x20(%%eax)	,%%xmm2	\n\t	unpckhpd	0x60(%%eax)	,%%xmm6	\n\t"\
		"movaps		%%xmm2, 0x20(%%eax)	\n\t	movaps		%%xmm6, 0x60(%%eax)	\n\t"\
		"\n\t"\
		"movaps		0x10(%%eax)	,%%xmm2	\n\t	movaps		0x50(%%eax)	,%%xmm6	\n\t"\
		"movaps		0x10(%%eax)	,%%xmm3	\n\t	movaps		0x50(%%eax)	,%%xmm7	\n\t"\
		"unpcklpd	0x30(%%eax)	,%%xmm2	\n\t	unpcklpd	0x70(%%eax)	,%%xmm6	\n\t"\
		"unpckhpd	0x30(%%eax)	,%%xmm3	\n\t	unpckhpd	0x70(%%eax)	,%%xmm7	\n\t"\
		"movaps		%%xmm2, 0x10(%%eax)	\n\t	movaps		%%xmm6, 0x50(%%eax)	\n\t"\
		"movaps		%%xmm3, 0x30(%%eax)	\n\t	movaps		%%xmm7, 0x70(%%eax)	\n\t"\
	"/**********************************************/\n\t"\
	"/*          Real      parts                   */\n\t"\
	"/**********************************************/\n\t"\
		"movl	%[__bjmod_0],	%%eax	\n\t"\
		"movaps		(%%eax)	,	%%xmm0	\n\t"\
		"movl	%[__sse_sw]	,	%%ebx	\n\t"\
		"movaps		(%%ebx)	,	%%xmm7	\n\t"\
		"psubd		%%xmm0	,	%%xmm7	\n\t"\
		"movmskps	%%xmm7	,	%%esi	\n\t"\
		"\n\t"\
		"shll	$24		,%%esi			\n\t"\
		"movaps		%%xmm0	,%%xmm7		\n\t"\
		"movl	%[__n_minus_sil],%%ecx	\n\t"\
		"movd	%%ecx	,%%xmm2			\n\t"\
		"pshufd	$0,	%%xmm2	,%%xmm2		\n\t"\
		"psubd		%%xmm0	,%%xmm2		\n\t"\
		"movmskps	%%xmm2	,%%ecx		\n\t"\
		"shll	$16		,%%ecx			\n\t"\
		"addl	%%ecx	,%%esi			\n\t"\
		"movl	%[__sinwt]	,%%edx		\n\t"\
		"movd	%%edx	,%%xmm3			\n\t"\
		"pshufd	$0,	%%xmm3	,%%xmm3		\n\t"\
		"psubd		%%xmm3	,%%xmm7		\n\t"\
		"movmskps	%%xmm7	,%%edx		\n\t"\
		"shll	$8		,%%edx			\n\t"\
		"addl	%%edx	,%%esi			\n\t"\
		"\n\t"\
		"movl	%[__half_arr]	,%%eax	\n\t"\
		"movl	%[__wtA]	,%%ebx		\n\t"\
		"movl	%[__wtB]	,%%ecx		\n\t"\
		"movaps		-0x10(%%eax),%%xmm4	\n\t"\
		"movaps		     (%%ebx),%%xmm2	\n\t	movaps		 0x10(%%ebx),%%xmm6	\n\t"\
		"movhpd		     (%%ecx),%%xmm3	\n\t	movhpd		-0x10(%%ecx),%%xmm7	\n\t"\
		"movlpd		 0x08(%%ecx),%%xmm3	\n\t	movlpd		-0x08(%%ecx),%%xmm7	\n\t"\
		"\n\t"\
		"subl	$0x20	,%%ecx			\n\t"\
		"movl	%%ecx	,%[__wtB]		\n\t"\
		"\n\t"\
		"shld	$20,	%%esi	,%%edi	\n\t	shld	$18,	%%esi	,%%ebx	\n\t"\
		"andl	$0x00000030	,%%edi		\n\t	andl	$0x00000030	,%%ebx		\n\t"\
		"shld	$28,	%%esi,	%%edx	\n\t	shld	$26,	%%esi,	%%ecx	\n\t"\
		"andl	$0x00000030	,%%edx		\n\t	andl	$0x00000030	,%%ecx		\n\t"\
		"addl	%%eax	,%%edi			\n\t	addl	%%eax	,%%ebx			\n\t"\
		"addl	%%eax	,%%edx			\n\t	addl	%%eax	,%%ecx			\n\t"\
		"\n\t"\
		"mulpd	 	0x100(%%eax),%%xmm2	\n\t	mulpd	 	0x100(%%eax),%%xmm6	\n\t"\
		"mulpd	 	0x110(%%eax),%%xmm3	\n\t	mulpd	 	0x110(%%eax),%%xmm7	\n\t"\
		"mulpd	 	     (%%edi),%%xmm2	\n\t	mulpd	 	     (%%ebx),%%xmm6	\n\t"\
		"mulpd	 	0x040(%%edx),%%xmm3	\n\t	mulpd	 	0x040(%%ecx),%%xmm7	\n\t"\
		"\n\t"\
		"movl	%[__cyA]	,%%ecx		\n\t	movl	%[__cyB]	,%%edx		\n\t"\
		"mulpd		%%xmm3	,%%xmm1		\n\t	mulpd		%%xmm7	,%%xmm5		\n\t"\
		"addpd		(%%ecx)	,%%xmm1		\n\t	addpd		(%%edx)	,%%xmm5		\n\t"\
		"movaps		%%xmm1	,%%xmm3		\n\t	movaps		%%xmm5	,%%xmm7		\n\t"\
		"addpd		%%xmm4	,%%xmm3		\n\t	addpd		%%xmm4	,%%xmm7		\n\t"\
		"subpd		%%xmm4	,%%xmm3		\n\t	subpd		%%xmm4	,%%xmm7		\n\t"\
		"\n\t"\
		"movl	%[__sign_mask],%%ebx	\n\t"\
		"subpd			%%xmm3	,%%xmm1	\n\t	subpd			%%xmm7	,%%xmm5	\n\t"\
		"andpd			(%%ebx)	,%%xmm1	\n\t	andpd			(%%ebx)	,%%xmm5	\n\t"\
		"maxpd			%%xmm5	,%%xmm1	\n\t"\
		"maxpd		-0x20(%%eax),%%xmm1	\n\t"\
		"movaps		%%xmm1,-0x20(%%eax)	\n\t"\
		"\n\t"\
		"shld	$12		,%%esi	,%%edi	\n\t	shld	$10		,%%esi	,%%ebx	\n\t"\
		"andl	$0x00000030	,%%edi		\n\t	andl	$0x00000030	,%%ebx		\n\t"\
		"addl	%%eax	,%%edi			\n\t	addl	%%eax	,%%ebx			\n\t"\
		"movaps		%%xmm3	,%%xmm1		\n\t	movaps		%%xmm7	,%%xmm5		\n\t"\
		"mulpd	 	0xc0(%%edi)	,%%xmm3	\n\t	mulpd	 	0xc0(%%ebx)	,%%xmm7	\n\t"\
		"addpd		%%xmm4	,%%xmm3		\n\t	addpd		%%xmm4	,%%xmm7		\n\t"\
		"subpd		%%xmm4	,%%xmm3		\n\t	subpd		%%xmm4	,%%xmm7		\n\t"\
		"movaps		%%xmm3	,(%%ecx)	\n\t	movaps		%%xmm7	,(%%edx)	\n\t"\
		"\n\t"\
		"\n\t"\
		"movl	%[__data]	,%%eax		\n\t"\
		"mulpd	 0x80(%%edi)	,%%xmm3	\n\t	mulpd	 0x80(%%ebx)	,%%xmm7	\n\t"\
		"subpd		%%xmm3	,	%%xmm1	\n\t	subpd		%%xmm7	,	%%xmm5	\n\t"\
		"mulpd		%%xmm2	,	%%xmm1	\n\t	mulpd		%%xmm6	,	%%xmm5	\n\t"\
		"movaps		%%xmm1	,    (%%eax)\n\t	movaps		%%xmm5	,0x40(%%eax)\n\t"\
		"\n\t"\
		"movl	%[__sse_n]	,%%ebx		\n\t"\
		"movaps		(%%ebx)	,%%xmm2		\n\t"\
		"movl	%[__sse_bw]	,%%eax		\n\t"\
		"paddd		(%%eax)	,%%xmm0		\n\t"\
		"movaps		%%xmm0	,%%xmm1		\n\t"\
		"pcmpgtd	%%xmm2	,%%xmm1		\n\t"\
		"pand		(%%ebx)	,%%xmm1		\n\t"\
		"psubd		%%xmm1	,%%xmm0		\n\t"\
	"/**********************************************/\n\t"\
	"/*          Imaginary parts               */\n\t"\
	"/**********************************************/\n\t"\
		"movl	%[__sse_sw]	,%%edx		\n\t"\
		"movaps	(%%edx)	,%%xmm1			\n\t"\
		"psubd	%%xmm0	,%%xmm1			\n\t"\
		"movmskps	%%xmm1	,%%esi		\n\t"\
		"shll	$24	,%%esi				\n\t"\
		"movaps	%%xmm0	,%%xmm1			\n\t"\
		"movl	%[__n_minus_silp1],%%ecx\n\t"\
		"movd	%%ecx	,%%xmm2			\n\t"\
		"pshufd	$0	,%%xmm2	,%%xmm2		\n\t"\
		"psubd	%%xmm0	,%%xmm2			\n\t"\
		"movmskps	%%xmm2	,%%ecx		\n\t"\
		"shll	$16	,%%ecx				\n\t"\
		"addl	%%ecx	,%%esi			\n\t"\
		"movl	%[__sinwtm1]	,%%edx	\n\t"\
		"movd	%%edx	,%%xmm3			\n\t"\
		"pshufd	$0	,%%xmm3	,%%xmm3		\n\t"\
		"psubd	%%xmm3	,%%xmm1			\n\t"\
		"movmskps	%%xmm1	,%%edx		\n\t"\
		"shll	$8	,%%edx				\n\t"\
		"addl	%%edx	,%%esi			\n\t"\
		"movl	%[__data]	,%%eax		\n\t"\
		"movaps	 0x10(%%eax)	,%%xmm1	\n\t"\
		"movaps	 0x50(%%eax)	,%%xmm5	\n\t"\
		"\n\t"\
		"movl	%[__half_arr]	,%%eax	\n\t"\
		"movl	%[__wtA]	,%%ebx		\n\t"\
		"movl	%[__wtC]	,%%ecx		\n\t"\
		"\n\t"\
		"movaps	     (%%ebx)	,%%xmm2	\n\t	movaps	 0x10(%%ebx)	,%%xmm6	\n\t"\
		"movhpd	     (%%ecx)	,%%xmm3	\n\t	movhpd	-0x10(%%ecx)	,%%xmm7	\n\t"\
		"movlpd	 0x08(%%ecx)	,%%xmm3	\n\t	movlpd	-0x08(%%ecx)	,%%xmm7	\n\t"\
		"\n\t"\
		"addl	$0x20	,%%ebx			\n\t"\
		"subl	$0x20	,%%ecx			\n\t"\
		"movl	%%ebx	,%[__wtA]		\n\t"\
		"movl	%%ecx	,%[__wtC]		\n\t"\
		"shld	$20		,%%esi	,%%edi	\n\t	shld	$18		,%%esi	,%%ebx	\n\t"\
		"andl	$0x00000030	,%%edi		\n\t	andl	$0x00000030	,%%ebx		\n\t"\
		"shld	$28,	%%esi,	%%edx	\n\t	shld	$26,	%%esi,	%%ecx	\n\t"\
		"andl	$0x00000030	,%%edx		\n\t	andl	$0x00000030	,%%ecx		\n\t"\
		"addl	%%eax	,%%edi			\n\t	addl	%%eax	,%%ebx			\n\t"\
		"addl	%%eax	,%%edx			\n\t	addl	%%eax	,%%ecx			\n\t"\
		"\n\t"\
		"mulpd	 0x120(%%eax)	,%%xmm2	\n\t	mulpd	 0x120(%%eax)	,%%xmm6	\n\t"\
		"mulpd	 0x130(%%eax)	,%%xmm3	\n\t	mulpd	 0x130(%%eax)	,%%xmm7	\n\t"\
		"mulpd	      (%%edi)	,%%xmm2	\n\t	mulpd	      (%%ebx)	,%%xmm6	\n\t"\
		"mulpd	 0x040(%%edx)	,%%xmm3	\n\t	mulpd	 0x040(%%ecx)	,%%xmm7	\n\t"\
		"\n\t"\
		"movl	%[__cyA]	,%%ecx		\n\t	movl	%[__cyB]	,%%edx		\n\t"\
		"mulpd	%%xmm3	,%%xmm1			\n\t	mulpd	%%xmm7	,%%xmm5			\n\t"\
		"addpd	(%%ecx)	,%%xmm1			\n\t	addpd	(%%edx)	,%%xmm5			\n\t"\
		"movaps	%%xmm1	,%%xmm3			\n\t	movaps	%%xmm5	,%%xmm7			\n\t"\
		"addpd	%%xmm4	,%%xmm3			\n\t	addpd	%%xmm4	,%%xmm7			\n\t"\
		"subpd	%%xmm4	,%%xmm3			\n\t	subpd	%%xmm4	,%%xmm7			\n\t"\
		"\n\t"\
		"movl	%[__sign_mask],%%ebx	\n\t"\
		"subpd			%%xmm3	,%%xmm1	\n\t	subpd			%%xmm7	,%%xmm5	\n\t"\
		"andpd			(%%ebx)	,%%xmm1	\n\t	andpd			(%%ebx)	,%%xmm5	\n\t"\
		"maxpd			%%xmm5	,%%xmm1	\n\t"\
		"maxpd		-0x20(%%eax),%%xmm1	\n\t"\
		"movaps		%%xmm1,-0x20(%%eax)	\n\t"\
		"\n\t"\
		"shld	$12		,%%esi	,%%edi	\n\t	shld	$10		,%%esi	,%%ebx	\n\t"\
		"andl	$0x00000030	,%%edi		\n\t	andl	$0x00000030	,%%ebx		\n\t"\
		"addl	%%eax	,%%edi			\n\t	addl	%%eax	,%%ebx			\n\t"\
		"movaps	%%xmm3	,%%xmm1			\n\t	movaps	%%xmm7	,%%xmm5			\n\t"\
		"mulpd	 0xc0(%%edi)	,%%xmm3	\n\t	mulpd	 0xc0(%%ebx)	,%%xmm7	\n\t"\
		"addpd	%%xmm4	,%%xmm3			\n\t	addpd	%%xmm4	,%%xmm7			\n\t"\
		"subpd	%%xmm4	,%%xmm3			\n\t	subpd	%%xmm4	,%%xmm7			\n\t"\
		"movaps	%%xmm3	,(%%ecx)		\n\t	movaps	%%xmm7	,(%%edx)		\n\t"\
		"\n\t"\
		"movl	%[__data]	,%%eax		\n\t"\
		"mulpd	 0x80(%%edi)	,%%xmm3	\n\t	mulpd	 0x80(%%ebx)	,%%xmm7	\n\t"\
		"subpd	%%xmm3	,%%xmm1			\n\t	subpd	%%xmm7	,%%xmm5			\n\t"\
		"mulpd	%%xmm2	,%%xmm1			\n\t	mulpd	%%xmm6	,%%xmm5			\n\t"\
		"movaps	%%xmm1	, 0x10(%%eax)	\n\t	movaps	%%xmm5	, 0x50(%%eax)	\n\t"\
		"\n\t"\
		"movl	%[__sse_n]	,%%ebx		\n\t"\
		"movaps		(%%ebx)	,%%xmm2		\n\t"\
		"movl	%[__sse_bw]	,%%eax		\n\t"\
		"paddd		(%%eax)	,%%xmm0		\n\t"\
		"movaps		%%xmm0	,%%xmm1		\n\t"\
		"pcmpgtd	%%xmm2	,%%xmm1		\n\t"\
		"pand		(%%ebx)	,%%xmm1		\n\t"\
		"psubd		%%xmm1	,%%xmm0		\n\t"\
		"movl	%[__bjmod_0],%%ecx		\n\t"\
		"movaps	%%xmm0,(%%ecx)			\n\t"\
	"popl %%ebx	\n\t"\
		:					/* outputs: none */\
		: [__data]		"m" (Xdata)	/* All inputs from memory addresses here */\
		, [__wtA]		"m" (XwtA)		\
		, [__wtB]		"m" (XwtB)		\
		, [__wtC]		"m" (XwtC)		\
		, [__cyA]		"m" (XcyA)		\
		, [__cyB]		"m" (XcyB)		\
		, [__bjmod_0]	"m" (Xbjmod_0)		\
		, [__half_arr]	"m" (Xhalf_arr)		\
		, [__n_minus_silp1] "m" (Xn_minus_silp1)\
		, [__n_minus_sil]	"m" (Xn_minus_sil)	\
		, [__sign_mask]	"m" (Xsign_mask)	\
		, [__sinwt]		"m" (Xsinwt)		\
		, [__sinwtm1]		"m" (Xsinwtm1)		\
		, [__sse_bw]		"m" (Xsse_bw)		\
		, [__sse_n]		"m" (Xsse_n)		\
		, [__sse_sw]		"m" (Xsse_sw)		\
		: "cc","memory","eax","ecx","edx","edi","esi"		/* Clobbered registers */\
	);\
	}

	#define SSE2_cmplx_carry_norm_errcheck2_2B(Xdata,XwtA,XwtB,XcyA,XcyB,Xbjmod_0,Xhalf_arr,Xn_minus_silp1,Xn_minus_sil,Xsign_mask,Xsinwt,Xsinwtm1,Xsse_bw,Xsse_n,Xsse_sw)\
	{\
	__asm__ volatile (\
	"pushl %%ebx	\n\t"/* Explicit save/restore of PIC register */\
	"/**********************************************/\n\t"\
	"/*          Real      parts                   */\n\t"\
	"/**********************************************/\n\t"\
		"movl	%[__bjmod_0],	%%eax	\n\t"\
		"movaps		(%%eax)	,	%%xmm0	\n\t"\
		"movl	%[__sse_sw]	,	%%ebx	\n\t"\
		"movaps		(%%ebx)	,	%%xmm1	\n\t"\
		"psubd		%%xmm0	,	%%xmm1	\n\t"\
		"movmskps	%%xmm1	,	%%esi	\n\t"\
		"\n\t"\
		"shll	$24		,%%esi			\n\t"\
		"movaps		%%xmm0	,%%xmm1		\n\t"\
		"movl	%[__n_minus_sil],%%ecx	\n\t"\
		"movd	%%ecx	,%%xmm2			\n\t"\
		"pshufd	$0,	%%xmm2	,%%xmm2		\n\t"\
		"psubd		%%xmm0	,%%xmm2		\n\t"\
		"movmskps	%%xmm2	,%%ecx		\n\t"\
		"shll	$16		,%%ecx			\n\t"\
		"addl	%%ecx	,%%esi			\n\t"\
		"movl	%[__sinwt]	,%%edx		\n\t"\
		"movd	%%edx	,%%xmm3			\n\t"\
		"pshufd	$0,	%%xmm3	,%%xmm3		\n\t"\
		"psubd		%%xmm3	,%%xmm1		\n\t"\
		"movmskps	%%xmm1	,%%edx		\n\t"\
		"shll	$8		,%%edx			\n\t"\
		"addl	%%edx	,%%esi			\n\t"\
		"\n\t"\
		"mov	%[__data],%%eax			\n\t"\
		"movaps	 0x20(%%eax),%%xmm1		\n\t	movaps		 0x60(%%eax),%%xmm5\n\t"\
		"\n\t"\
		"movl	%[__half_arr]	,%%eax	\n\t"\
		"movl	%[__wtA]	,%%ebx		\n\t"\
		"movl	%[__wtB]	,%%ecx		\n\t"\
		"movaps		-0x10(%%eax),%%xmm4	\n\t"\
		"movaps		     (%%ebx),%%xmm2	\n\t	movaps		 0x10(%%ebx),%%xmm6	\n\t"\
		"movhpd		     (%%ecx),%%xmm3	\n\t	movhpd		-0x10(%%ecx),%%xmm7	\n\t"\
		"movlpd		 0x08(%%ecx),%%xmm3	\n\t	movlpd		-0x08(%%ecx),%%xmm7	\n\t"\
		"\n\t"\
		"shld	$20,	%%esi	,%%edi	\n\t	shld	$18,	%%esi	,%%ebx	\n\t"\
		"andl	$0x00000030	,%%edi		\n\t	andl	$0x00000030	,%%ebx		\n\t"\
		"shld	$28,	%%esi,	%%edx	\n\t	shld	$26,	%%esi,	%%ecx	\n\t"\
		"andl	$0x00000030	,%%edx		\n\t	andl	$0x00000030	,%%ecx		\n\t"\
		"addl	%%eax	,%%edi			\n\t	addl	%%eax	,%%ebx			\n\t"\
		"addl	%%eax	,%%edx			\n\t	addl	%%eax	,%%ecx			\n\t"\
		"\n\t"\
		"mulpd	 	0x100(%%eax),%%xmm2	\n\t	mulpd	 	0x100(%%eax),%%xmm6	\n\t"\
		"mulpd	 	0x110(%%eax),%%xmm3	\n\t	mulpd	 	0x110(%%eax),%%xmm7	\n\t"\
		"mulpd	 	     (%%edi),%%xmm2	\n\t	mulpd	 	     (%%ebx),%%xmm6	\n\t"\
		"mulpd	 	0x040(%%edx),%%xmm3	\n\t	mulpd	 	0x040(%%ecx),%%xmm7	\n\t"\
		"\n\t"\
		"movl	%[__cyA]	,%%ecx		\n\t	movl	%[__cyB]	,%%edx		\n\t"\
		"mulpd		%%xmm3	,%%xmm1		\n\t	mulpd		%%xmm7	,%%xmm5		\n\t"\
		"addpd		(%%ecx)	,%%xmm1		\n\t	addpd		(%%edx)	,%%xmm5		\n\t"\
		"movaps		%%xmm1	,%%xmm3		\n\t	movaps		%%xmm5	,%%xmm7		\n\t"\
		"addpd		%%xmm4	,%%xmm3		\n\t	addpd		%%xmm4	,%%xmm7		\n\t"\
		"subpd		%%xmm4	,%%xmm3		\n\t	subpd		%%xmm4	,%%xmm7		\n\t"\
		"\n\t"\
		"movl	%[__sign_mask],%%ebx	\n\t"\
		"subpd			%%xmm3	,%%xmm1	\n\t	subpd			%%xmm7	,%%xmm5	\n\t"\
		"andpd			(%%ebx)	,%%xmm1	\n\t	andpd			(%%ebx)	,%%xmm5	\n\t"\
		"maxpd			%%xmm5	,%%xmm1	\n\t"\
		"maxpd		-0x20(%%eax),%%xmm1	\n\t"\
		"movaps		%%xmm1,-0x20(%%eax)	\n\t"\
		"\n\t"\
		"shld	$12		,%%esi	,%%edi	\n\t	shld	$10		,%%esi	,%%ebx	\n\t"\
		"andl	$0x00000030	,%%edi		\n\t	andl	$0x00000030	,%%ebx		\n\t"\
		"addl	%%eax	,%%edi			\n\t	addl	%%eax	,%%ebx			\n\t"\
		"movaps		%%xmm3	,%%xmm1		\n\t	movaps		%%xmm7	,%%xmm5		\n\t"\
		"mulpd	 	0xc0(%%edi)	,%%xmm3	\n\t	mulpd	 	0xc0(%%ebx)	,%%xmm7	\n\t"\
		"addpd		%%xmm4	,%%xmm3		\n\t	addpd		%%xmm4	,%%xmm7		\n\t"\
		"subpd		%%xmm4	,%%xmm3		\n\t	subpd		%%xmm4	,%%xmm7		\n\t"\
		"movaps		%%xmm3	,(%%ecx)	\n\t	movaps		%%xmm7	,(%%edx)	\n\t"\
		"\n\t"\
		"\n\t"\
		"movl	%[__data]	,%%eax		\n\t"\
		"mulpd	 0x80(%%edi)	,%%xmm3	\n\t	mulpd	 0x80(%%ebx)	,%%xmm7	\n\t"\
		"subpd		%%xmm3	,	%%xmm1	\n\t	subpd		%%xmm7	,	%%xmm5	\n\t"\
		"mulpd		%%xmm2	,	%%xmm1	\n\t	mulpd		%%xmm6	,	%%xmm5	\n\t"\
		"movaps		%%xmm1	,0x20(%%eax)\n\t	movaps		%%xmm5	,0x60(%%eax)\n\t"\
		"\n\t"\
		"movl	%[__sse_n]	,%%ebx		\n\t"\
		"movaps		(%%ebx)	,%%xmm2		\n\t"\
		"movl	%[__sse_bw]	,%%eax		\n\t"\
		"paddd		(%%eax)	,%%xmm0		\n\t"\
		"movaps		%%xmm0	,%%xmm1		\n\t"\
		"pcmpgtd	%%xmm2	,%%xmm1		\n\t"\
		"pand		(%%ebx)	,%%xmm1		\n\t"\
		"psubd		%%xmm1	,%%xmm0		\n\t"\
	"/**********************************************/\n\t"\
	"/*          Imaginary parts               */\n\t"\
	"/**********************************************/\n\t"\
		"movl	%[__sse_sw]	,%%edx		\n\t"\
		"movaps	(%%edx)	,%%xmm1			\n\t"\
		"psubd	%%xmm0	,%%xmm1			\n\t"\
		"movmskps	%%xmm1	,%%esi		\n\t"\
		"shll	$24	,%%esi				\n\t"\
		"movaps	%%xmm0	,%%xmm1			\n\t"\
		"movl	%[__n_minus_silp1],%%ecx\n\t"\
		"movd	%%ecx	,%%xmm2			\n\t"\
		"pshufd	$0	,%%xmm2	,%%xmm2		\n\t"\
		"psubd	%%xmm0	,%%xmm2			\n\t"\
		"movmskps	%%xmm2	,%%ecx		\n\t"\
		"shll	$16	,%%ecx				\n\t"\
		"addl	%%ecx	,%%esi			\n\t"\
		"movl	%[__sinwtm1]	,%%edx	\n\t"\
		"movd	%%edx	,%%xmm3			\n\t"\
		"pshufd	$0	,%%xmm3	,%%xmm3		\n\t"\
		"psubd	%%xmm3	,%%xmm1			\n\t"\
		"movmskps	%%xmm1	,%%edx		\n\t"\
		"shll	$8	,%%edx				\n\t"\
		"addl	%%edx	,%%esi			\n\t"\
		"movl	%[__data]	,%%eax		\n\t"\
		"movaps	 0x30(%%eax)	,%%xmm1	\n\t"\
		"movaps	 0x70(%%eax)	,%%xmm5	\n\t"\
		"\n\t"\
		"movl	%[__half_arr]	,%%eax	\n\t"\
		"movl	%[__wtA]	,%%ebx		\n\t"\
		"movl	%[__wtB]	,%%ecx		\n\t"\
		"\n\t"\
		"movaps	     (%%ebx)	,%%xmm2	\n\t	movaps	 0x10(%%ebx)	,%%xmm6	\n\t"\
		"movhpd	     (%%ecx)	,%%xmm3	\n\t	movhpd	-0x10(%%ecx)	,%%xmm7	\n\t"\
		"movlpd	 0x08(%%ecx)	,%%xmm3	\n\t	movlpd	-0x08(%%ecx)	,%%xmm7	\n\t"\
		"\n\t"\
		"addl	$0x20	,%%ebx			\n\t"\
		"subl	$0x20	,%%ecx			\n\t"\
		"movl	%%ebx	,%[__wtA]		\n\t"\
		"movl	%%ecx	,%[__wtB]		\n\t"\
		"shld	$20		,%%esi	,%%edi	\n\t	shld	$18		,%%esi	,%%ebx	\n\t"\
		"andl	$0x00000030	,%%edi		\n\t	andl	$0x00000030	,%%ebx		\n\t"\
		"shld	$28,	%%esi,	%%edx	\n\t	shld	$26,	%%esi,	%%ecx	\n\t"\
		"andl	$0x00000030	,%%edx		\n\t	andl	$0x00000030	,%%ecx		\n\t"\
		"addl	%%eax	,%%edi			\n\t	addl	%%eax	,%%ebx			\n\t"\
		"addl	%%eax	,%%edx			\n\t	addl	%%eax	,%%ecx			\n\t"\
		"\n\t"\
		"mulpd	 0x120(%%eax)	,%%xmm2	\n\t	mulpd	 0x120(%%eax)	,%%xmm6	\n\t"\
		"mulpd	 0x130(%%eax)	,%%xmm3	\n\t	mulpd	 0x130(%%eax)	,%%xmm7	\n\t"\
		"mulpd	      (%%edi)	,%%xmm2	\n\t	mulpd	      (%%ebx)	,%%xmm6	\n\t"\
		"mulpd	 0x040(%%edx)	,%%xmm3	\n\t	mulpd	 0x040(%%ecx)	,%%xmm7	\n\t"\
		"\n\t"\
		"movl	%[__cyA]	,%%ecx		\n\t	movl	%[__cyB]	,%%edx		\n\t"\
		"mulpd	%%xmm3	,%%xmm1			\n\t	mulpd	%%xmm7	,%%xmm5			\n\t"\
		"addpd	(%%ecx)	,%%xmm1			\n\t	addpd	(%%edx)	,%%xmm5			\n\t"\
		"movaps	%%xmm1	,%%xmm3			\n\t	movaps	%%xmm5	,%%xmm7			\n\t"\
		"addpd	%%xmm4	,%%xmm3			\n\t	addpd	%%xmm4	,%%xmm7			\n\t"\
		"subpd	%%xmm4	,%%xmm3			\n\t	subpd	%%xmm4	,%%xmm7			\n\t"\
		"\n\t"\
		"movl	%[__sign_mask],%%ebx	\n\t"\
		"subpd			%%xmm3	,%%xmm1	\n\t	subpd			%%xmm7	,%%xmm5	\n\t"\
		"andpd			(%%ebx)	,%%xmm1	\n\t	andpd			(%%ebx)	,%%xmm5	\n\t"\
		"maxpd			%%xmm5	,%%xmm1	\n\t"\
		"maxpd		-0x20(%%eax),%%xmm1	\n\t"\
		"movaps		%%xmm1,-0x20(%%eax)	\n\t"\
		"\n\t"\
		"shld	$12		,%%esi	,%%edi	\n\t	shld	$10		,%%esi	,%%ebx	\n\t"\
		"andl	$0x00000030	,%%edi		\n\t	andl	$0x00000030	,%%ebx		\n\t"\
		"addl	%%eax	,%%edi			\n\t	addl	%%eax	,%%ebx			\n\t"\
		"movaps	%%xmm3	,%%xmm1			\n\t	movaps	%%xmm7	,%%xmm5			\n\t"\
		"mulpd	 0xc0(%%edi)	,%%xmm3	\n\t	mulpd	 0xc0(%%ebx)	,%%xmm7	\n\t"\
		"addpd	%%xmm4	,%%xmm3			\n\t	addpd	%%xmm4	,%%xmm7			\n\t"\
		"subpd	%%xmm4	,%%xmm3			\n\t	subpd	%%xmm4	,%%xmm7			\n\t"\
		"movaps	%%xmm3	,(%%ecx)		\n\t	movaps	%%xmm7	,(%%edx)		\n\t"\
		"\n\t"\
		"\n\t"\
		"movl	%[__data]	,%%eax		\n\t"\
		"mulpd	 0x80(%%edi)	,%%xmm3	\n\t	mulpd	 0x80(%%ebx)	,%%xmm7	\n\t"\
		"subpd	%%xmm3	,%%xmm1			\n\t	subpd	%%xmm7	,%%xmm5			\n\t"\
		"mulpd	%%xmm2	,%%xmm1			\n\t	mulpd	%%xmm6	,%%xmm5			\n\t"\
		"movaps	%%xmm1	, 0x30(%%eax)	\n\t	movaps	%%xmm5	, 0x70(%%eax)	\n\t"\
		"\n\t"\
		"movl	%[__sse_n]	,%%ebx		\n\t"\
		"movaps		(%%ebx)	,%%xmm2		\n\t"\
		"movl	%[__sse_bw]	,%%eax		\n\t"\
		"paddd		(%%eax)	,%%xmm0		\n\t"\
		"movaps		%%xmm0	,%%xmm1		\n\t"\
		"pcmpgtd	%%xmm2	,%%xmm1		\n\t"\
		"pand		(%%ebx)	,%%xmm1		\n\t"\
		"psubd		%%xmm1	,%%xmm0		\n\t"\
		"movl	%[__bjmod_0],%%ecx		\n\t"\
		"movaps	%%xmm0,(%%ecx)			\n\t"\
	"/**********************************************/\n\t"\
	"/*          Repack the data:              */\n\t"\
	"/**********************************************/\n\t"\
		"mov	%[__data],%%eax			\n\t"\
		"movaps		0x10(%%eax)	,%%xmm1	\n\t	movaps		0x50(%%eax)	,%%xmm5	\n\t"\
		"movaps		    (%%eax)	,%%xmm0	\n\t	movaps		0x40(%%eax)	,%%xmm4	\n\t"\
		"movaps		%%xmm1		,%%xmm3	\n\t	movaps		%%xmm5		,%%xmm7	\n\t"\
		"movaps		%%xmm0		,%%xmm2	\n\t	movaps		%%xmm4		,%%xmm6	\n\t"\
		"unpckhpd	0x30(%%eax)	,%%xmm3	\n\t	unpckhpd	0x70(%%eax)	,%%xmm7	\n\t"\
		"unpcklpd	0x30(%%eax)	,%%xmm1	\n\t	unpcklpd	0x70(%%eax)	,%%xmm5	\n\t"\
		"movaps		%%xmm3,0x30(%%eax)	\n\t	movaps		%%xmm7,0x70(%%eax)	\n\t"\
		"unpckhpd	0x20(%%eax)	,%%xmm2	\n\t	unpckhpd	0x60(%%eax)	,%%xmm6	\n\t"\
		"unpcklpd	0x20(%%eax)	,%%xmm0	\n\t	unpcklpd	0x60(%%eax)	,%%xmm4	\n\t"\
		"movaps		%%xmm2,0x20(%%eax)	\n\t	movaps		%%xmm6,0x60(%%eax)	\n\t"\
		"movaps		%%xmm1,0x10(%%eax)	\n\t	movaps		%%xmm5,0x50(%%eax)	\n\t"\
		"movaps		%%xmm0,    (%%eax)	\n\t	movaps		%%xmm4,0x40(%%eax)	\n\t"\
	"popl %%ebx	\n\t"\
		:					/* outputs: none */\
		: [__data]		"m" (Xdata)	/* All inputs from memory addresses here */\
		, [__wtA]		"m" (XwtA)		\
		, [__wtB]		"m" (XwtB)		\
		, [__cyA]		"m" (XcyA)		\
		, [__cyB]		"m" (XcyB)		\
		, [__bjmod_0]	"m" (Xbjmod_0)		\
		, [__half_arr]	"m" (Xhalf_arr)		\
		, [__n_minus_silp1] "m" (Xn_minus_silp1)\
		, [__n_minus_sil]	"m" (Xn_minus_sil)	\
		, [__sign_mask]	"m" (Xsign_mask)	\
		, [__sinwt]		"m" (Xsinwt)		\
		, [__sinwtm1]		"m" (Xsinwtm1)		\
		, [__sse_bw]		"m" (Xsse_bw)		\
		, [__sse_n]		"m" (Xsse_n)		\
		, [__sse_sw]		"m" (Xsse_sw)		\
		: "cc","memory","eax","ecx","edx","edi","esi"		/* Clobbered registers */\
	);\
	}

	/******************************************************************************************************************************************************************/
	/********** No-ROE-Check versions of the latter 2 of the above 3 macros - the first is called too infrequently to bother with a special non-ROE version: **********/
	/******************************************************************************************************************************************************************/

	#define SSE2_cmplx_carry_norm_nocheck1_2B(Xdata,XwtA,XwtB,XwtC,XcyA,XcyB,Xbjmod_0,Xhalf_arr,Xn_minus_silp1,Xn_minus_sil,Xsinwt,Xsinwtm1,Xsse_bw,Xsse_n,Xsse_sw)\
	{\
	__asm__ volatile (\
	"pushl %%ebx	\n\t"/* Explicit save/restore of PIC register */\
	"/***************Unpack the data:*************************/\n\t"\
		"movl	%[__data]	,%%eax	\n\t"\
		"movaps		    (%%eax)	,%%xmm1	\n\t	movaps		0x40(%%eax)	,%%xmm5	\n\t"\
		"movaps		    (%%eax)	,%%xmm2	\n\t	movaps		0x40(%%eax)	,%%xmm6	\n\t"\
		"unpcklpd	0x20(%%eax)	,%%xmm1	\n\t	unpcklpd	0x60(%%eax)	,%%xmm5	\n\t"\
		"unpckhpd	0x20(%%eax)	,%%xmm2	\n\t	unpckhpd	0x60(%%eax)	,%%xmm6	\n\t"\
		"movaps		%%xmm2, 0x20(%%eax)	\n\t	movaps		%%xmm6, 0x60(%%eax)	\n\t"\
		"\n\t"\
		"movaps		0x10(%%eax)	,%%xmm2	\n\t	movaps		0x50(%%eax)	,%%xmm6	\n\t"\
		"movaps		0x10(%%eax)	,%%xmm3	\n\t	movaps		0x50(%%eax)	,%%xmm7	\n\t"\
		"unpcklpd	0x30(%%eax)	,%%xmm2	\n\t	unpcklpd	0x70(%%eax)	,%%xmm6	\n\t"\
		"unpckhpd	0x30(%%eax)	,%%xmm3	\n\t	unpckhpd	0x70(%%eax)	,%%xmm7	\n\t"\
		"movaps		%%xmm2, 0x10(%%eax)	\n\t	movaps		%%xmm6, 0x50(%%eax)	\n\t"\
		"movaps		%%xmm3, 0x30(%%eax)	\n\t	movaps		%%xmm7, 0x70(%%eax)	\n\t"\
	"/**********************************************/\n\t"\
	"/*          Real      parts                   */\n\t"\
	"/**********************************************/\n\t"\
		"movl	%[__bjmod_0],	%%eax	\n\t"\
		"movaps		(%%eax)	,	%%xmm0	\n\t"\
		"movl	%[__sse_sw]	,	%%ebx	\n\t"\
		"movaps		(%%ebx)	,	%%xmm7	\n\t"\
		"psubd		%%xmm0	,	%%xmm7	\n\t"\
		"movmskps	%%xmm7	,	%%esi	\n\t"\
		"\n\t"\
		"shll	$24		,%%esi			\n\t"\
		"movaps		%%xmm0	,%%xmm7		\n\t"\
		"movl	%[__n_minus_sil],%%ecx	\n\t"\
		"movd	%%ecx	,%%xmm2			\n\t"\
		"pshufd	$0,	%%xmm2	,%%xmm2		\n\t"\
		"psubd		%%xmm0	,%%xmm2		\n\t"\
		"movmskps	%%xmm2	,%%ecx		\n\t"\
		"shll	$16		,%%ecx			\n\t"\
		"addl	%%ecx	,%%esi			\n\t"\
		"movl	%[__sinwt]	,%%edx		\n\t"\
		"movd	%%edx	,%%xmm3			\n\t"\
		"pshufd	$0,	%%xmm3	,%%xmm3		\n\t"\
		"psubd		%%xmm3	,%%xmm7		\n\t"\
		"movmskps	%%xmm7	,%%edx		\n\t"\
		"shll	$8		,%%edx			\n\t"\
		"addl	%%edx	,%%esi			\n\t"\
		"\n\t"\
		"movl	%[__half_arr]	,%%eax	\n\t"\
		"movl	%[__wtA]	,%%ebx		\n\t"\
		"movl	%[__wtB]	,%%ecx		\n\t"\
		"movaps		-0x10(%%eax),%%xmm4	\n\t"\
		"movaps		     (%%ebx),%%xmm2	\n\t	movaps		 0x10(%%ebx),%%xmm6	\n\t"\
		"movhpd		     (%%ecx),%%xmm3	\n\t	movhpd		-0x10(%%ecx),%%xmm7	\n\t"\
		"movlpd		 0x08(%%ecx),%%xmm3	\n\t	movlpd		-0x08(%%ecx),%%xmm7	\n\t"\
		"\n\t"\
		"subl	$0x20	,%%ecx			\n\t"\
		"movl	%%ecx	,%[__wtB]		\n\t"\
		"\n\t"\
		"shld	$20,	%%esi	,%%edi	\n\t	shld	$18,	%%esi	,%%ebx	\n\t"\
		"andl	$0x00000030	,%%edi		\n\t	andl	$0x00000030	,%%ebx		\n\t"\
		"shld	$28,	%%esi,	%%edx	\n\t	shld	$26,	%%esi,	%%ecx	\n\t"\
		"andl	$0x00000030	,%%edx		\n\t	andl	$0x00000030	,%%ecx		\n\t"\
		"addl	%%eax	,%%edi			\n\t	addl	%%eax	,%%ebx			\n\t"\
		"addl	%%eax	,%%edx			\n\t	addl	%%eax	,%%ecx			\n\t"\
		"\n\t"\
		"mulpd	 	0x100(%%eax),%%xmm2	\n\t	mulpd	 	0x100(%%eax),%%xmm6	\n\t"\
		"mulpd	 	0x110(%%eax),%%xmm3	\n\t	mulpd	 	0x110(%%eax),%%xmm7	\n\t"\
		"mulpd	 	     (%%edi),%%xmm2	\n\t	mulpd	 	     (%%ebx),%%xmm6	\n\t"\
		"mulpd	 	0x040(%%edx),%%xmm3	\n\t	mulpd	 	0x040(%%ecx),%%xmm7	\n\t"\
		"\n\t"\
		"movl	%[__cyA]	,%%ecx		\n\t	movl	%[__cyB]	,%%edx		\n\t"\
		"mulpd		%%xmm3	,%%xmm1		\n\t	mulpd		%%xmm7	,%%xmm5		\n\t"\
		"addpd		(%%ecx)	,%%xmm1		\n\t	addpd		(%%edx)	,%%xmm5		\n\t"\
		"movaps		%%xmm1	,%%xmm3		\n\t	movaps		%%xmm5	,%%xmm7		\n\t"\
		"addpd		%%xmm4	,%%xmm3		\n\t	addpd		%%xmm4	,%%xmm7		\n\t"\
		"subpd		%%xmm4	,%%xmm3		\n\t	subpd		%%xmm4	,%%xmm7		\n\t"\
		"\n\t"\
		"/* NO ROE HERE */				\n\t"\
		"\n\t"\
		"shld	$12		,%%esi	,%%edi	\n\t	shld	$10		,%%esi	,%%ebx	\n\t"\
		"andl	$0x00000030	,%%edi		\n\t	andl	$0x00000030	,%%ebx		\n\t"\
		"addl	%%eax	,%%edi			\n\t	addl	%%eax	,%%ebx			\n\t"\
		"movaps		%%xmm3	,%%xmm1		\n\t	movaps		%%xmm7	,%%xmm5		\n\t"\
		"mulpd	 	0xc0(%%edi)	,%%xmm3	\n\t	mulpd	 	0xc0(%%ebx)	,%%xmm7	\n\t"\
		"addpd		%%xmm4	,%%xmm3		\n\t	addpd		%%xmm4	,%%xmm7		\n\t"\
		"subpd		%%xmm4	,%%xmm3		\n\t	subpd		%%xmm4	,%%xmm7		\n\t"\
		"movaps		%%xmm3	,(%%ecx)	\n\t	movaps		%%xmm7	,(%%edx)	\n\t"\
		"\n\t"\
		"\n\t"\
		"movl	%[__data]	,%%eax		\n\t"\
		"mulpd	 0x80(%%edi)	,%%xmm3	\n\t	mulpd	 0x80(%%ebx)	,%%xmm7	\n\t"\
		"subpd		%%xmm3	,	%%xmm1	\n\t	subpd		%%xmm7	,	%%xmm5	\n\t"\
		"mulpd		%%xmm2	,	%%xmm1	\n\t	mulpd		%%xmm6	,	%%xmm5	\n\t"\
		"movaps		%%xmm1	,    (%%eax)\n\t	movaps		%%xmm5	,0x40(%%eax)\n\t"\
		"\n\t"\
		"movl	%[__sse_n]	,%%ebx		\n\t"\
		"movaps		(%%ebx)	,%%xmm2		\n\t"\
		"movl	%[__sse_bw]	,%%eax		\n\t"\
		"paddd		(%%eax)	,%%xmm0		\n\t"\
		"movaps		%%xmm0	,%%xmm1		\n\t"\
		"pcmpgtd	%%xmm2	,%%xmm1		\n\t"\
		"pand		(%%ebx)	,%%xmm1		\n\t"\
		"psubd		%%xmm1	,%%xmm0		\n\t"\
	"/**********************************************/\n\t"\
	"/*          Imaginary parts               */\n\t"\
	"/**********************************************/\n\t"\
		"movl	%[__sse_sw]	,%%edx		\n\t"\
		"movaps	(%%edx)	,%%xmm1			\n\t"\
		"psubd	%%xmm0	,%%xmm1			\n\t"\
		"movmskps	%%xmm1	,%%esi		\n\t"\
		"shll	$24	,%%esi				\n\t"\
		"movaps	%%xmm0	,%%xmm1			\n\t"\
		"movl	%[__n_minus_silp1],%%ecx\n\t"\
		"movd	%%ecx	,%%xmm2			\n\t"\
		"pshufd	$0	,%%xmm2	,%%xmm2		\n\t"\
		"psubd	%%xmm0	,%%xmm2			\n\t"\
		"movmskps	%%xmm2	,%%ecx		\n\t"\
		"shll	$16	,%%ecx				\n\t"\
		"addl	%%ecx	,%%esi			\n\t"\
		"movl	%[__sinwtm1]	,%%edx	\n\t"\
		"movd	%%edx	,%%xmm3			\n\t"\
		"pshufd	$0	,%%xmm3	,%%xmm3		\n\t"\
		"psubd	%%xmm3	,%%xmm1			\n\t"\
		"movmskps	%%xmm1	,%%edx		\n\t"\
		"shll	$8	,%%edx				\n\t"\
		"addl	%%edx	,%%esi			\n\t"\
		"movl	%[__data]	,%%eax		\n\t"\
		"movaps	 0x10(%%eax)	,%%xmm1	\n\t"\
		"movaps	 0x50(%%eax)	,%%xmm5	\n\t"\
		"\n\t"\
		"movl	%[__half_arr]	,%%eax	\n\t"\
		"movl	%[__wtA]	,%%ebx		\n\t"\
		"movl	%[__wtC]	,%%ecx		\n\t"\
		"\n\t"\
		"movaps	     (%%ebx)	,%%xmm2	\n\t	movaps	 0x10(%%ebx)	,%%xmm6	\n\t"\
		"movhpd	     (%%ecx)	,%%xmm3	\n\t	movhpd	-0x10(%%ecx)	,%%xmm7	\n\t"\
		"movlpd	 0x08(%%ecx)	,%%xmm3	\n\t	movlpd	-0x08(%%ecx)	,%%xmm7	\n\t"\
		"\n\t"\
		"addl	$0x20	,%%ebx			\n\t"\
		"subl	$0x20	,%%ecx			\n\t"\
		"movl	%%ebx	,%[__wtA]		\n\t"\
		"movl	%%ecx	,%[__wtC]		\n\t"\
		"shld	$20		,%%esi	,%%edi	\n\t	shld	$18		,%%esi	,%%ebx	\n\t"\
		"andl	$0x00000030	,%%edi		\n\t	andl	$0x00000030	,%%ebx		\n\t"\
		"shld	$28,	%%esi,	%%edx	\n\t	shld	$26,	%%esi,	%%ecx	\n\t"\
		"andl	$0x00000030	,%%edx		\n\t	andl	$0x00000030	,%%ecx		\n\t"\
		"addl	%%eax	,%%edi			\n\t	addl	%%eax	,%%ebx			\n\t"\
		"addl	%%eax	,%%edx			\n\t	addl	%%eax	,%%ecx			\n\t"\
		"\n\t"\
		"mulpd	 0x120(%%eax)	,%%xmm2	\n\t	mulpd	 0x120(%%eax)	,%%xmm6	\n\t"\
		"mulpd	 0x130(%%eax)	,%%xmm3	\n\t	mulpd	 0x130(%%eax)	,%%xmm7	\n\t"\
		"mulpd	      (%%edi)	,%%xmm2	\n\t	mulpd	      (%%ebx)	,%%xmm6	\n\t"\
		"mulpd	 0x040(%%edx)	,%%xmm3	\n\t	mulpd	 0x040(%%ecx)	,%%xmm7	\n\t"\
		"\n\t"\
		"movl	%[__cyA]	,%%ecx		\n\t	movl	%[__cyB]	,%%edx		\n\t"\
		"mulpd	%%xmm3	,%%xmm1			\n\t	mulpd	%%xmm7	,%%xmm5			\n\t"\
		"addpd	(%%ecx)	,%%xmm1			\n\t	addpd	(%%edx)	,%%xmm5			\n\t"\
		"movaps	%%xmm1	,%%xmm3			\n\t	movaps	%%xmm5	,%%xmm7			\n\t"\
		"addpd	%%xmm4	,%%xmm3			\n\t	addpd	%%xmm4	,%%xmm7			\n\t"\
		"subpd	%%xmm4	,%%xmm3			\n\t	subpd	%%xmm4	,%%xmm7			\n\t"\
		"\n\t"\
		"/* NO ROE HERE */				\n\t"\
		"\n\t"\
		"shld	$12		,%%esi	,%%edi	\n\t	shld	$10		,%%esi	,%%ebx	\n\t"\
		"andl	$0x00000030	,%%edi		\n\t	andl	$0x00000030	,%%ebx		\n\t"\
		"addl	%%eax	,%%edi			\n\t	addl	%%eax	,%%ebx			\n\t"\
		"movaps	%%xmm3	,%%xmm1			\n\t	movaps	%%xmm7	,%%xmm5			\n\t"\
		"mulpd	 0xc0(%%edi)	,%%xmm3	\n\t	mulpd	 0xc0(%%ebx)	,%%xmm7	\n\t"\
		"addpd	%%xmm4	,%%xmm3			\n\t	addpd	%%xmm4	,%%xmm7			\n\t"\
		"subpd	%%xmm4	,%%xmm3			\n\t	subpd	%%xmm4	,%%xmm7			\n\t"\
		"movaps	%%xmm3	,(%%ecx)		\n\t	movaps	%%xmm7	,(%%edx)		\n\t"\
		"\n\t"\
		"movl	%[__data]	,%%eax		\n\t"\
		"mulpd	 0x80(%%edi)	,%%xmm3	\n\t	mulpd	 0x80(%%ebx)	,%%xmm7	\n\t"\
		"subpd	%%xmm3	,%%xmm1			\n\t	subpd	%%xmm7	,%%xmm5			\n\t"\
		"mulpd	%%xmm2	,%%xmm1			\n\t	mulpd	%%xmm6	,%%xmm5			\n\t"\
		"movaps	%%xmm1	, 0x10(%%eax)	\n\t	movaps	%%xmm5	, 0x50(%%eax)	\n\t"\
		"\n\t"\
		"movl	%[__sse_n]	,%%ebx		\n\t"\
		"movaps		(%%ebx)	,%%xmm2		\n\t"\
		"movl	%[__sse_bw]	,%%eax		\n\t"\
		"paddd		(%%eax)	,%%xmm0		\n\t"\
		"movaps		%%xmm0	,%%xmm1		\n\t"\
		"pcmpgtd	%%xmm2	,%%xmm1		\n\t"\
		"pand		(%%ebx)	,%%xmm1		\n\t"\
		"psubd		%%xmm1	,%%xmm0		\n\t"\
		"movl	%[__bjmod_0],%%ecx		\n\t"\
		"movaps	%%xmm0,(%%ecx)			\n\t"\
	"popl %%ebx	\n\t"\
		:					/* outputs: none */\
		: [__data]		"m" (Xdata)	/* All inputs from memory addresses here */\
		, [__wtA]		"m" (XwtA)		\
		, [__wtB]		"m" (XwtB)		\
		, [__wtC]		"m" (XwtC)		\
		, [__cyA]		"m" (XcyA)		\
		, [__cyB]		"m" (XcyB)		\
		, [__bjmod_0]	"m" (Xbjmod_0)		\
		, [__half_arr]	"m" (Xhalf_arr)		\
		, [__n_minus_silp1] "m" (Xn_minus_silp1)\
		, [__n_minus_sil]	"m" (Xn_minus_sil)	\
		, [__sinwt]		"m" (Xsinwt)		\
		, [__sinwtm1]		"m" (Xsinwtm1)		\
		, [__sse_bw]		"m" (Xsse_bw)		\
		, [__sse_n]		"m" (Xsse_n)		\
		, [__sse_sw]		"m" (Xsse_sw)		\
		: "cc","memory","eax","ecx","edx","edi","esi"		/* Clobbered registers */\
	);\
	}

	#define SSE2_cmplx_carry_norm_nocheck2_2B(Xdata,XwtA,XwtB,XcyA,XcyB,Xbjmod_0,Xhalf_arr,Xn_minus_silp1,Xn_minus_sil,Xsinwt,Xsinwtm1,Xsse_bw,Xsse_n,Xsse_sw)\
	{\
	__asm__ volatile (\
	"pushl %%ebx	\n\t"/* Explicit save/restore of PIC register */\
	"/**********************************************/\n\t"\
	"/*          Real      parts                   */\n\t"\
	"/**********************************************/\n\t"\
		"movl	%[__bjmod_0],	%%eax	\n\t"\
		"movaps		(%%eax)	,	%%xmm0	\n\t"\
		"movl	%[__sse_sw]	,	%%ebx	\n\t"\
		"movaps		(%%ebx)	,	%%xmm1	\n\t"\
		"psubd		%%xmm0	,	%%xmm1	\n\t"\
		"movmskps	%%xmm1	,	%%esi	\n\t"\
		"\n\t"\
		"shll	$24		,%%esi			\n\t"\
		"movaps		%%xmm0	,%%xmm1		\n\t"\
		"movl	%[__n_minus_sil],%%ecx	\n\t"\
		"movd	%%ecx	,%%xmm2			\n\t"\
		"pshufd	$0,	%%xmm2	,%%xmm2		\n\t"\
		"psubd		%%xmm0	,%%xmm2		\n\t"\
		"movmskps	%%xmm2	,%%ecx		\n\t"\
		"shll	$16		,%%ecx			\n\t"\
		"addl	%%ecx	,%%esi			\n\t"\
		"movl	%[__sinwt]	,%%edx		\n\t"\
		"movd	%%edx	,%%xmm3			\n\t"\
		"pshufd	$0,	%%xmm3	,%%xmm3		\n\t"\
		"psubd		%%xmm3	,%%xmm1		\n\t"\
		"movmskps	%%xmm1	,%%edx		\n\t"\
		"shll	$8		,%%edx			\n\t"\
		"addl	%%edx	,%%esi			\n\t"\
		"\n\t"\
		"mov	%[__data],%%eax			\n\t"\
		"movaps	 0x20(%%eax),%%xmm1		\n\t	movaps		 0x60(%%eax),%%xmm5\n\t"\
		"\n\t"\
		"movl	%[__half_arr]	,%%eax	\n\t"\
		"movl	%[__wtA]	,%%ebx		\n\t"\
		"movl	%[__wtB]	,%%ecx		\n\t"\
		"movaps		-0x10(%%eax),%%xmm4	\n\t"\
		"movaps		     (%%ebx),%%xmm2	\n\t	movaps		 0x10(%%ebx),%%xmm6	\n\t"\
		"movhpd		     (%%ecx),%%xmm3	\n\t	movhpd		-0x10(%%ecx),%%xmm7	\n\t"\
		"movlpd		 0x08(%%ecx),%%xmm3	\n\t	movlpd		-0x08(%%ecx),%%xmm7	\n\t"\
		"\n\t"\
		"shld	$20,	%%esi	,%%edi	\n\t	shld	$18,	%%esi	,%%ebx	\n\t"\
		"andl	$0x00000030	,%%edi		\n\t	andl	$0x00000030	,%%ebx		\n\t"\
		"shld	$28,	%%esi,	%%edx	\n\t	shld	$26,	%%esi,	%%ecx	\n\t"\
		"andl	$0x00000030	,%%edx		\n\t	andl	$0x00000030	,%%ecx		\n\t"\
		"addl	%%eax	,%%edi			\n\t	addl	%%eax	,%%ebx			\n\t"\
		"addl	%%eax	,%%edx			\n\t	addl	%%eax	,%%ecx			\n\t"\
		"\n\t"\
		"mulpd	 	0x100(%%eax),%%xmm2	\n\t	mulpd	 	0x100(%%eax),%%xmm6	\n\t"\
		"mulpd	 	0x110(%%eax),%%xmm3	\n\t	mulpd	 	0x110(%%eax),%%xmm7	\n\t"\
		"mulpd	 	     (%%edi),%%xmm2	\n\t	mulpd	 	     (%%ebx),%%xmm6	\n\t"\
		"mulpd	 	0x040(%%edx),%%xmm3	\n\t	mulpd	 	0x040(%%ecx),%%xmm7	\n\t"\
		"\n\t"\
		"movl	%[__cyA]	,%%ecx		\n\t	movl	%[__cyB]	,%%edx		\n\t"\
		"mulpd		%%xmm3	,%%xmm1		\n\t	mulpd		%%xmm7	,%%xmm5		\n\t"\
		"addpd		(%%ecx)	,%%xmm1		\n\t	addpd		(%%edx)	,%%xmm5		\n\t"\
		"movaps		%%xmm1	,%%xmm3		\n\t	movaps		%%xmm5	,%%xmm7		\n\t"\
		"addpd		%%xmm4	,%%xmm3		\n\t	addpd		%%xmm4	,%%xmm7		\n\t"\
		"subpd		%%xmm4	,%%xmm3		\n\t	subpd		%%xmm4	,%%xmm7		\n\t"\
		"\n\t"\
		"/* NO ROE HERE */				\n\t"\
		"\n\t"\
		"shld	$12		,%%esi	,%%edi	\n\t	shld	$10		,%%esi	,%%ebx	\n\t"\
		"andl	$0x00000030	,%%edi		\n\t	andl	$0x00000030	,%%ebx		\n\t"\
		"addl	%%eax	,%%edi			\n\t	addl	%%eax	,%%ebx			\n\t"\
		"movaps		%%xmm3	,%%xmm1		\n\t	movaps		%%xmm7	,%%xmm5		\n\t"\
		"mulpd	 	0xc0(%%edi)	,%%xmm3	\n\t	mulpd	 	0xc0(%%ebx)	,%%xmm7	\n\t"\
		"addpd		%%xmm4	,%%xmm3		\n\t	addpd		%%xmm4	,%%xmm7		\n\t"\
		"subpd		%%xmm4	,%%xmm3		\n\t	subpd		%%xmm4	,%%xmm7		\n\t"\
		"movaps		%%xmm3	,(%%ecx)	\n\t	movaps		%%xmm7	,(%%edx)	\n\t"\
		"\n\t"\
		"\n\t"\
		"movl	%[__data]	,%%eax		\n\t"\
		"mulpd	 0x80(%%edi)	,%%xmm3	\n\t	mulpd	 0x80(%%ebx)	,%%xmm7	\n\t"\
		"subpd		%%xmm3	,	%%xmm1	\n\t	subpd		%%xmm7	,	%%xmm5	\n\t"\
		"mulpd		%%xmm2	,	%%xmm1	\n\t	mulpd		%%xmm6	,	%%xmm5	\n\t"\
		"movaps		%%xmm1	,0x20(%%eax)\n\t	movaps		%%xmm5	,0x60(%%eax)\n\t"\
		"\n\t"\
		"movl	%[__sse_n]	,%%ebx		\n\t"\
		"movaps		(%%ebx)	,%%xmm2		\n\t"\
		"movl	%[__sse_bw]	,%%eax		\n\t"\
		"paddd		(%%eax)	,%%xmm0		\n\t"\
		"movaps		%%xmm0	,%%xmm1		\n\t"\
		"pcmpgtd	%%xmm2	,%%xmm1		\n\t"\
		"pand		(%%ebx)	,%%xmm1		\n\t"\
		"psubd		%%xmm1	,%%xmm0		\n\t"\
	"/**********************************************/\n\t"\
	"/*          Imaginary parts               */\n\t"\
	"/**********************************************/\n\t"\
		"movl	%[__sse_sw]	,%%edx		\n\t"\
		"movaps	(%%edx)	,%%xmm1			\n\t"\
		"psubd	%%xmm0	,%%xmm1			\n\t"\
		"movmskps	%%xmm1	,%%esi		\n\t"\
		"shll	$24	,%%esi				\n\t"\
		"movaps	%%xmm0	,%%xmm1			\n\t"\
		"movl	%[__n_minus_silp1],%%ecx\n\t"\
		"movd	%%ecx	,%%xmm2			\n\t"\
		"pshufd	$0	,%%xmm2	,%%xmm2		\n\t"\
		"psubd	%%xmm0	,%%xmm2			\n\t"\
		"movmskps	%%xmm2	,%%ecx		\n\t"\
		"shll	$16	,%%ecx				\n\t"\
		"addl	%%ecx	,%%esi			\n\t"\
		"movl	%[__sinwtm1]	,%%edx	\n\t"\
		"movd	%%edx	,%%xmm3			\n\t"\
		"pshufd	$0	,%%xmm3	,%%xmm3		\n\t"\
		"psubd	%%xmm3	,%%xmm1			\n\t"\
		"movmskps	%%xmm1	,%%edx		\n\t"\
		"shll	$8	,%%edx				\n\t"\
		"addl	%%edx	,%%esi			\n\t"\
		"movl	%[__data]	,%%eax		\n\t"\
		"movaps	 0x30(%%eax)	,%%xmm1	\n\t"\
		"movaps	 0x70(%%eax)	,%%xmm5	\n\t"\
		"\n\t"\
		"movl	%[__half_arr]	,%%eax	\n\t"\
		"movl	%[__wtA]	,%%ebx		\n\t"\
		"movl	%[__wtB]	,%%ecx		\n\t"\
		"\n\t"\
		"movaps	     (%%ebx)	,%%xmm2	\n\t	movaps	 0x10(%%ebx)	,%%xmm6	\n\t"\
		"movhpd	     (%%ecx)	,%%xmm3	\n\t	movhpd	-0x10(%%ecx)	,%%xmm7	\n\t"\
		"movlpd	 0x08(%%ecx)	,%%xmm3	\n\t	movlpd	-0x08(%%ecx)	,%%xmm7	\n\t"\
		"\n\t"\
		"addl	$0x20	,%%ebx			\n\t"\
		"subl	$0x20	,%%ecx			\n\t"\
		"movl	%%ebx	,%[__wtA]		\n\t"\
		"movl	%%ecx	,%[__wtB]		\n\t"\
		"shld	$20		,%%esi	,%%edi	\n\t	shld	$18		,%%esi	,%%ebx	\n\t"\
		"andl	$0x00000030	,%%edi		\n\t	andl	$0x00000030	,%%ebx		\n\t"\
		"shld	$28,	%%esi,	%%edx	\n\t	shld	$26,	%%esi,	%%ecx	\n\t"\
		"andl	$0x00000030	,%%edx		\n\t	andl	$0x00000030	,%%ecx		\n\t"\
		"addl	%%eax	,%%edi			\n\t	addl	%%eax	,%%ebx			\n\t"\
		"addl	%%eax	,%%edx			\n\t	addl	%%eax	,%%ecx			\n\t"\
		"\n\t"\
		"mulpd	 0x120(%%eax)	,%%xmm2	\n\t	mulpd	 0x120(%%eax)	,%%xmm6	\n\t"\
		"mulpd	 0x130(%%eax)	,%%xmm3	\n\t	mulpd	 0x130(%%eax)	,%%xmm7	\n\t"\
		"mulpd	      (%%edi)	,%%xmm2	\n\t	mulpd	      (%%ebx)	,%%xmm6	\n\t"\
		"mulpd	 0x040(%%edx)	,%%xmm3	\n\t	mulpd	 0x040(%%ecx)	,%%xmm7	\n\t"\
		"\n\t"\
		"movl	%[__cyA]	,%%ecx		\n\t	movl	%[__cyB]	,%%edx		\n\t"\
		"mulpd	%%xmm3	,%%xmm1			\n\t	mulpd	%%xmm7	,%%xmm5			\n\t"\
		"addpd	(%%ecx)	,%%xmm1			\n\t	addpd	(%%edx)	,%%xmm5			\n\t"\
		"movaps	%%xmm1	,%%xmm3			\n\t	movaps	%%xmm5	,%%xmm7			\n\t"\
		"addpd	%%xmm4	,%%xmm3			\n\t	addpd	%%xmm4	,%%xmm7			\n\t"\
		"subpd	%%xmm4	,%%xmm3			\n\t	subpd	%%xmm4	,%%xmm7			\n\t"\
		"\n\t"\
		"/* NO ROE HERE */				\n\t"\
		"\n\t"\
		"shld	$12		,%%esi	,%%edi	\n\t	shld	$10		,%%esi	,%%ebx	\n\t"\
		"andl	$0x00000030	,%%edi		\n\t	andl	$0x00000030	,%%ebx		\n\t"\
		"addl	%%eax	,%%edi			\n\t	addl	%%eax	,%%ebx			\n\t"\
		"movaps	%%xmm3	,%%xmm1			\n\t	movaps	%%xmm7	,%%xmm5			\n\t"\
		"mulpd	 0xc0(%%edi)	,%%xmm3	\n\t	mulpd	 0xc0(%%ebx)	,%%xmm7	\n\t"\
		"addpd	%%xmm4	,%%xmm3			\n\t	addpd	%%xmm4	,%%xmm7			\n\t"\
		"subpd	%%xmm4	,%%xmm3			\n\t	subpd	%%xmm4	,%%xmm7			\n\t"\
		"movaps	%%xmm3	,(%%ecx)		\n\t	movaps	%%xmm7	,(%%edx)		\n\t"\
		"\n\t"\
		"\n\t"\
		"movl	%[__data]	,%%eax		\n\t"\
		"mulpd	 0x80(%%edi)	,%%xmm3	\n\t	mulpd	 0x80(%%ebx)	,%%xmm7	\n\t"\
		"subpd	%%xmm3	,%%xmm1			\n\t	subpd	%%xmm7	,%%xmm5			\n\t"\
		"mulpd	%%xmm2	,%%xmm1			\n\t	mulpd	%%xmm6	,%%xmm5			\n\t"\
		"movaps	%%xmm1	, 0x30(%%eax)	\n\t	movaps	%%xmm5	, 0x70(%%eax)	\n\t"\
		"\n\t"\
		"movl	%[__sse_n]	,%%ebx		\n\t"\
		"movaps		(%%ebx)	,%%xmm2		\n\t"\
		"movl	%[__sse_bw]	,%%eax		\n\t"\
		"paddd		(%%eax)	,%%xmm0		\n\t"\
		"movaps		%%xmm0	,%%xmm1		\n\t"\
		"pcmpgtd	%%xmm2	,%%xmm1		\n\t"\
		"pand		(%%ebx)	,%%xmm1		\n\t"\
		"psubd		%%xmm1	,%%xmm0		\n\t"\
		"movl	%[__bjmod_0],%%ecx		\n\t"\
		"movaps	%%xmm0,(%%ecx)			\n\t"\
	"/**********************************************/\n\t"\
	"/*          Repack the data:              */\n\t"\
	"/**********************************************/\n\t"\
		"mov	%[__data],%%eax			\n\t"\
		"movaps		0x10(%%eax)	,%%xmm1	\n\t	movaps		0x50(%%eax)	,%%xmm5	\n\t"\
		"movaps		    (%%eax)	,%%xmm0	\n\t	movaps		0x40(%%eax)	,%%xmm4	\n\t"\
		"movaps		%%xmm1		,%%xmm3	\n\t	movaps		%%xmm5		,%%xmm7	\n\t"\
		"movaps		%%xmm0		,%%xmm2	\n\t	movaps		%%xmm4		,%%xmm6	\n\t"\
		"unpckhpd	0x30(%%eax)	,%%xmm3	\n\t	unpckhpd	0x70(%%eax)	,%%xmm7	\n\t"\
		"unpcklpd	0x30(%%eax)	,%%xmm1	\n\t	unpcklpd	0x70(%%eax)	,%%xmm5	\n\t"\
		"movaps		%%xmm3,0x30(%%eax)	\n\t	movaps		%%xmm7,0x70(%%eax)	\n\t"\
		"unpckhpd	0x20(%%eax)	,%%xmm2	\n\t	unpckhpd	0x60(%%eax)	,%%xmm6	\n\t"\
		"unpcklpd	0x20(%%eax)	,%%xmm0	\n\t	unpcklpd	0x60(%%eax)	,%%xmm4	\n\t"\
		"movaps		%%xmm2,0x20(%%eax)	\n\t	movaps		%%xmm6,0x60(%%eax)	\n\t"\
		"movaps		%%xmm1,0x10(%%eax)	\n\t	movaps		%%xmm5,0x50(%%eax)	\n\t"\
		"movaps		%%xmm0,    (%%eax)	\n\t	movaps		%%xmm4,0x40(%%eax)	\n\t"\
	"popl %%ebx	\n\t"\
		:					/* outputs: none */\
		: [__data]		"m" (Xdata)	/* All inputs from memory addresses here */\
		, [__wtA]		"m" (XwtA)		\
		, [__wtB]		"m" (XwtB)		\
		, [__cyA]		"m" (XcyA)		\
		, [__cyB]		"m" (XcyB)		\
		, [__bjmod_0]	"m" (Xbjmod_0)		\
		, [__half_arr]	"m" (Xhalf_arr)		\
		, [__n_minus_silp1] "m" (Xn_minus_silp1)\
		, [__n_minus_sil]	"m" (Xn_minus_sil)	\
		, [__sinwt]		"m" (Xsinwt)		\
		, [__sinwtm1]		"m" (Xsinwtm1)		\
		, [__sse_bw]		"m" (Xsse_bw)		\
		, [__sse_n]		"m" (Xsse_n)		\
		, [__sse_sw]		"m" (Xsse_sw)		\
		: "cc","memory","eax","ecx","edx","edi","esi"		/* Clobbered registers */\
	);\
	}

#endif	/* carry_gcc_h_included */

