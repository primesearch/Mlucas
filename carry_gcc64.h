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

	/* Power-of-2-runlength Fermat-mod acyclic-transform carry macro. (No IBDWT needed for power-of-2 runlenghts).
	*/
	#define SSE2_fermat_carry_norm_pow2_errcheck(Xdata,Xcy,Xnrt_bits,Xnrtm1,Xidx_offset,Xidx_incr,Xhalf_arr,Xsign_mask,Xadd1,Xadd2)\
	{\
	__asm__ volatile (\
		"movslq	%[__idx_offset],%%rsi	\n\t"\
		"movslq		%[__nrt_bits],%%rcx	\n\t"\
		"movslq		%[__nrtm1],%%rdi	\n\t"\
		"movq		%%rsi,%%rax			\n\t"\
		"shrq		$1,%%rax			\n\t"\
		"movq		%%rax,%%rbx			\n\t"\
		"andq		%%rdi,%%rax	\n\t"\
		"shrq		%%cl,%%rbx			\n\t"\
		"shlq		$4,%%rax			\n\t"\
		"shlq		$4,%%rbx			\n\t"\
		"addq		%[__add1],%%rax		\n\t"\
		"addq		%[__add2],%%rbx		\n\t"\
		"movaps		(%%rax),%%xmm0		\n\t"\
		"movaps		(%%rbx),%%xmm1		\n\t"\
		"movq		%%rsi,%%rax			\n\t"\
		"movaps		%%xmm1,%%xmm2		\n\t"\
		"shufpd	$1,	%%xmm2,%%xmm2		\n\t"\
		"mulpd		%%xmm0,%%xmm1		\n\t"\
		"mulpd		%%xmm0,%%xmm2		\n\t"\
		"addq		$2,%%rax			\n\t"\
		"shrq		$1,%%rax			\n\t"\
		"movq		%%rax,%%rbx			\n\t"\
		"andq		%%rdi,%%rax	\n\t"\
		"shrq		%%cl,%%rbx			\n\t"\
		"shlq		$4,%%rax			\n\t"\
		"shlq		$4,%%rbx			\n\t"\
		"addq		%[__add1],%%rax		\n\t"\
		"addq		%[__add2],%%rbx		\n\t"\
		"movaps		(%%rax),%%xmm0		\n\t"\
		"movaps		(%%rbx),%%xmm3		\n\t"\
		"movq		%%rsi,%%rax			\n\t"\
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
		"movq		%[__half_arr],%%rcx	\n\t"\
		"movq		%[__data],%%rdx		\n\t"\
		"movaps		     (%%rdx),%%xmm4	\n\t"\
		"movaps		 0x10(%%rdx),%%xmm2	\n\t"\
		"movaps		0x020(%%rcx),%%xmm5	\n\t"\
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
		"movq		%[__cy],%%rbx		\n\t"\
		"movaps		%%xmm4,%%xmm5		\n\t"\
		"shufpd	$0,	%%xmm2,%%xmm4		\n\t"\
		"shufpd	$3,	%%xmm2,%%xmm5		\n\t"\
		"addpd		     (%%rbx),%%xmm4	\n\t"\
		"movaps		-0x20(%%rcx),%%xmm6	\n\t"\
		"movaps		-0x10(%%rcx),%%xmm7	\n\t"\
		"movaps		%%xmm4,%%xmm2		\n\t"\
		"addpd		%%xmm7,%%xmm4		\n\t"\
		"subpd		%%xmm7,%%xmm4		\n\t"\
		"movq		%[__sign_mask],%%rax\n\t"\
		"subpd		%%xmm4,%%xmm2		\n\t"\
		"andpd		     (%%rax),%%xmm2	\n\t"\
		"maxpd		%%xmm6,%%xmm2		\n\t"\
		"movaps		%%xmm2,%%xmm6		\n\t"\
		"movaps		%%xmm4,%%xmm2		\n\t"\
		"mulpd		0x10(%%rcx),%%xmm2	\n\t"\
		"addpd		%%xmm7,%%xmm2		\n\t"\
		"subpd		%%xmm7,%%xmm2		\n\t"\
		"movaps		%%xmm2,%%xmm3		\n\t"\
		"mulpd		    (%%rcx),%%xmm3	\n\t"\
		"subpd		%%xmm3,%%xmm4		\n\t"\
		"addpd		%%xmm2,%%xmm5		\n\t"\
		"movaps		%%xmm5,%%xmm2		\n\t"\
		"addpd		%%xmm7,%%xmm5		\n\t"\
		"subpd		%%xmm7,%%xmm5		\n\t"\
		"subpd		%%xmm5,%%xmm2		\n\t"\
		"andpd		     (%%rax),%%xmm2	\n\t"\
		"maxpd		%%xmm6,%%xmm2		\n\t"\
		"movaps		%%xmm2,%%xmm6		\n\t"\
		"movaps		%%xmm5,%%xmm2		\n\t"\
		"mulpd		 0x10(%%rcx),%%xmm2	\n\t"\
		"addpd		%%xmm7,%%xmm2		\n\t"\
		"subpd		%%xmm7,%%xmm2		\n\t"\
		"movaps		%%xmm2,%%xmm3		\n\t"\
		"mulpd		     (%%rcx),%%xmm3	\n\t"\
		"subpd		%%xmm3,%%xmm5		\n\t"\
		"movaps		%%xmm2,(%%rbx)		\n\t"\
		"movaps		%%xmm4,%%xmm2		\n\t"\
		"shufpd	$0,	%%xmm5,%%xmm4		\n\t"\
		"shufpd	$3,	%%xmm5,%%xmm2		\n\t"\
		"movaps		%%xmm4,%%xmm5		\n\t"\
		"movaps		%%xmm2,%%xmm3		\n\t"\
		"movaps		%%xmm6,-0x20(%%rcx)	\n\t"\
		"mulpd		%%xmm0,%%xmm4		\n\t"\
		"mulpd		%%xmm1,%%xmm3		\n\t"\
		"mulpd		%%xmm0,%%xmm2		\n\t"\
		"mulpd		%%xmm1,%%xmm5		\n\t"\
		"subpd		%%xmm3,%%xmm4		\n\t"\
		"addpd		%%xmm2,%%xmm5		\n\t"\
		"movaps		%%xmm4,    (%%rdx)	\n\t"\
		"movaps		%%xmm5,0x10(%%rdx)	\n\t"\
		"movslq	%[__idx_incr],%%rdi		\n\t"\
		"addq	%%rdi,%%rsi		\n\t"\
		"mov	%%esi, %[__idx_offset]	/* Store incremented idx_offset */	\n\t"\
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
		: "cc","memory","rax","rbx","rcx","rdx","rsi","rdi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"	/* Clobbered registers */\
	);\
	}

	/* Same power-of-2-transform deal as above, but use xmm8-15 to process 2 sets of carries side-by-side.
	Data/Carry #2 assumed offset by +0x20/0x10 from #1 (which are accessed via the [__data/__cy] pointers, resp.):
	*/
	#define SSE2_fermat_carry_norm_pow2_errcheck_X2(Xdata,Xcy,Xnrt_bits,Xnrtm1,Xidx_offset,Xidx_incr,Xhalf_arr,Xsign_mask,Xadd1,Xadd2)\
	{\
	__asm__ volatile (\
		"/* lcol -> rcol index analogs: [rsi,rax,rbx] -> [r10,r11,r12], [rcx,rdx,rdi] shared */\n\t"\
		"movslq	%[__idx_offset],%%rsi	\n\t		movslq	%[__idx_incr],%%r10		\n\t"\
		"movslq		%[__nrt_bits],%%rcx	\n\t		addq	%%rsi,%%r10				\n\t"\
		"movslq		%[__nrtm1],%%rdi	/* r10 contains idx_offset2, i.e. is the rcol-analog of rsi in lcol: */\n\t"\
		"movq		%%rsi,%%rax			\n\t		movq		%%r10,%%r11			\n\t"\
		"shrq		$1,%%rax			\n\t		shrq		$1,%%r11			\n\t"\
		"movq		%%rax,%%rbx			\n\t		movq		%%r11,%%r12			\n\t"\
		"andq		%%rdi,%%rax			\n\t		andq		%%rdi,%%r11			\n\t"\
		"shrq		%%cl,%%rbx			\n\t		shrq		%%cl,%%r12			\n\t"\
		"shlq		$4,%%rax			\n\t		shlq		$4,%%r11			\n\t"\
		"shlq		$4,%%rbx			\n\t		shlq		$4,%%r12			\n\t"\
		"addq		%[__add1],%%rax		\n\t		addq		%[__add1],%%r11		\n\t"\
		"addq		%[__add2],%%rbx		\n\t		addq		%[__add2],%%r12		\n\t"\
		"movaps		(%%rax),%%xmm0		\n\t		movaps		(%%r11),%%xmm8 		\n\t"\
		"movaps		(%%rbx),%%xmm1		\n\t		movaps		(%%r12),%%xmm9 		\n\t"\
		"movq		%%rsi,%%rax			\n\t		movq		%%r10,%%r11			\n\t"\
		"movaps		%%xmm1,%%xmm2		\n\t		movaps		%%xmm9 ,%%xmm10		\n\t"\
		"shufpd	$1,	%%xmm2,%%xmm2		\n\t		shufpd	$1,	%%xmm10,%%xmm10		\n\t"\
		"mulpd		%%xmm0,%%xmm1		\n\t		mulpd		%%xmm8 ,%%xmm9 		\n\t"\
		"mulpd		%%xmm0,%%xmm2		\n\t		mulpd		%%xmm8 ,%%xmm10		\n\t"\
		"addq		$2,%%rax			\n\t		addq		$2,%%r11			\n\t"\
		"shrq		$1,%%rax			\n\t		shrq		$1,%%r11			\n\t"\
		"movq		%%rax,%%rbx			\n\t		movq		%%r11,%%r12			\n\t"\
		"andq		%%rdi,%%rax			\n\t		andq		%%rdi,%%r11			\n\t"\
		"shrq		%%cl,%%rbx			\n\t		shrq		%%cl,%%r12			\n\t"\
		"shlq		$4,%%rax			\n\t		shlq		$4,%%r11			\n\t"\
		"shlq		$4,%%rbx			\n\t		shlq		$4,%%r12			\n\t"\
		"addq		%[__add1],%%rax		\n\t		addq		%[__add1],%%r11		\n\t"\
		"addq		%[__add2],%%rbx		\n\t		addq		%[__add2],%%r12		\n\t"\
		"movaps		(%%rax),%%xmm0		\n\t		movaps		(%%r11),%%xmm8 		\n\t"\
		"movaps		(%%rbx),%%xmm3		\n\t		movaps		(%%r12),%%xmm11		\n\t"\
		"movq		%%rsi,%%rax			\n\t		movq		%%r10,%%r11			\n\t"\
		"movaps		%%xmm3,%%xmm4		\n\t		movaps		%%xmm11,%%xmm12		\n\t"\
		"shufpd	$1,	%%xmm4,%%xmm4		\n\t		shufpd	$1,	%%xmm12,%%xmm12		\n\t"\
		"mulpd		%%xmm0,%%xmm3		\n\t		mulpd		%%xmm8 ,%%xmm11		\n\t"\
		"mulpd		%%xmm0,%%xmm4		\n\t		mulpd		%%xmm8 ,%%xmm12		\n\t"\
		"movaps		%%xmm1,%%xmm0		\n\t		movaps		%%xmm9 ,%%xmm8 		\n\t"\
		"unpcklpd	%%xmm3,%%xmm0		\n\t		unpcklpd	%%xmm11,%%xmm8 		\n\t"\
		"unpckhpd	%%xmm3,%%xmm1		\n\t		unpckhpd	%%xmm11,%%xmm9 		\n\t"\
		"subpd		%%xmm1,%%xmm0		\n\t		subpd		%%xmm9 ,%%xmm8 		\n\t"\
		"movaps		%%xmm2,%%xmm1		\n\t		movaps		%%xmm10,%%xmm9 		\n\t"\
		"unpcklpd	%%xmm4,%%xmm1		\n\t		unpcklpd	%%xmm12,%%xmm9 		\n\t"\
		"unpckhpd	%%xmm4,%%xmm2		\n\t		unpckhpd	%%xmm12,%%xmm10		\n\t"\
		"addpd		%%xmm2,%%xmm1		\n\t		addpd		%%xmm10,%%xmm9 		\n\t"\
		"movq		%[__half_arr],%%rcx	/* rcx shared, has same offset lcol/rcol: */\n\t"\
		"movq		%[__data],%%rdx		/* rdx shared, offset +0x20 in rcol: */		\n\t"\
		"movaps		    (%%rdx),%%xmm4	\n\t		movaps		0x20(%%rdx),%%xmm12	\n\t"\
		"movaps		0x10(%%rdx),%%xmm2	\n\t		movaps		0x30(%%rdx),%%xmm10	\n\t"\
		"movaps		0x20(%%rcx),%%xmm5	\n\t		movaps		0x20(%%rcx),%%xmm13	\n\t"\
		"mulpd		%%xmm5,%%xmm4		\n\t		mulpd		%%xmm13,%%xmm12		\n\t"\
		"mulpd		%%xmm5,%%xmm2		\n\t		mulpd		%%xmm13,%%xmm10		\n\t"\
		"movaps		%%xmm4,%%xmm5		\n\t		movaps		%%xmm12,%%xmm13		\n\t"\
		"movaps		%%xmm2,%%xmm3		\n\t		movaps		%%xmm10,%%xmm11		\n\t"\
		"mulpd		%%xmm0,%%xmm4		\n\t		mulpd		%%xmm8 ,%%xmm12		\n\t"\
		"mulpd		%%xmm1,%%xmm3		\n\t		mulpd		%%xmm9 ,%%xmm11		\n\t"\
		"mulpd		%%xmm0,%%xmm2		\n\t		mulpd		%%xmm8 ,%%xmm10		\n\t"\
		"mulpd		%%xmm1,%%xmm5		\n\t		mulpd		%%xmm9 ,%%xmm13		\n\t"\
		"addpd		%%xmm3,%%xmm4		\n\t		addpd		%%xmm11,%%xmm12		\n\t"\
		"subpd		%%xmm5,%%xmm2		\n\t		subpd		%%xmm13,%%xmm10		\n\t"\
		"movq		%[__cy],%%rbx	/* rbx -> rbx+0x10 (carry offset only half of data-offset) in rcol, shared from here */	\n\t"\
		"movaps		%%xmm4,%%xmm5		\n\t		movaps		%%xmm12,%%xmm13		\n\t"\
		"shufpd	$0,	%%xmm2,%%xmm4		\n\t		shufpd	$0,	%%xmm10,%%xmm12		\n\t"\
		"shufpd	$3,	%%xmm2,%%xmm5		\n\t		shufpd	$3,	%%xmm10,%%xmm13		\n\t"\
		"addpd		     (%%rbx),%%xmm4	\n\t		addpd		 0x10(%%rbx),%%xmm12\n\t"\
		"movaps		-0x20(%%rcx),%%xmm6	\n\t	movaps	-0x20(%%rcx),%%xmm14	/* Use 2 copies of maxerr, merge at end */\n\t"\
		"movaps		%%xmm4,%%xmm2		\n\t		movaps		%%xmm12,%%xmm10		\n\t"\
		"roundpd	$0,%%xmm4,%%xmm4	\n\t		roundpd	$0,%%xmm12,%%xmm12		\n\t"\
		"movq		%[__sign_mask],%%rax/* rax shared between lcol/rcol from here */\n\t"\
		"subpd		%%xmm4,%%xmm2		\n\t		subpd		%%xmm12,%%xmm10		\n\t"\
		"andpd		     (%%rax),%%xmm2	\n\t		andpd		     (%%rax),%%xmm10\n\t"\
		"maxpd		%%xmm6,%%xmm2		\n\t		maxpd		%%xmm14,%%xmm10		\n\t"\
		"movaps		%%xmm2,%%xmm6		\n\t		movaps		%%xmm10,%%xmm14		\n\t"\
		"movaps		%%xmm4,%%xmm2		\n\t		movaps		%%xmm12,%%xmm10		\n\t"\
		"mulpd		0x10(%%rcx),%%xmm2	\n\t		mulpd		0x10(%%rcx),%%xmm10	\n\t"\
		"roundpd	$0,%%xmm2,%%xmm2	\n\t		roundpd	$0,%%xmm10,%%xmm10		\n\t"\
		"movaps		%%xmm2,%%xmm3		\n\t		movaps		%%xmm10,%%xmm11		\n\t"\
		"mulpd		    (%%rcx),%%xmm3	\n\t		mulpd		    (%%rcx),%%xmm11	\n\t"\
		"subpd		%%xmm3,%%xmm4		\n\t		subpd		%%xmm11,%%xmm12		\n\t"\
		"addpd		%%xmm2,%%xmm5		\n\t		addpd		%%xmm10,%%xmm13		\n\t"\
		"movaps		%%xmm5,%%xmm2		\n\t		movaps		%%xmm13,%%xmm10		\n\t"\
		"roundpd	$0,%%xmm5,%%xmm5	\n\t		roundpd	$0,%%xmm13,%%xmm13		\n\t"\
		"subpd		%%xmm5,%%xmm2		\n\t		subpd		%%xmm13,%%xmm10		\n\t"\
		"andpd		     (%%rax),%%xmm2	\n\t		andpd		     (%%rax),%%xmm10\n\t"\
		"maxpd		%%xmm6,%%xmm2		\n\t		maxpd		%%xmm14,%%xmm10		\n\t"\
		"movaps		%%xmm2,%%xmm6		\n\t		movaps		%%xmm10,%%xmm14		\n\t"\
		"movaps		%%xmm5,%%xmm2		\n\t		movaps		%%xmm13,%%xmm10		\n\t"\
		"mulpd		 0x10(%%rcx),%%xmm2	\n\t		mulpd		 0x10(%%rcx),%%xmm10\n\t"\
		"roundpd	$0,%%xmm2,%%xmm2	\n\t		roundpd	$0,%%xmm10,%%xmm10		\n\t"\
		"movaps		%%xmm2,%%xmm3		\n\t		movaps		%%xmm10,%%xmm11		\n\t"\
		"mulpd		     (%%rcx),%%xmm3	\n\t		mulpd		     (%%rcx),%%xmm11\n\t"\
		"subpd		%%xmm3,%%xmm5		\n\t		subpd		%%xmm11,%%xmm13		\n\t"\
		"movaps		%%xmm2,(%%rbx)		\n\t		movaps		%%xmm10,0x10(%%rbx)	\n\t"\
		"movaps		%%xmm4,%%xmm2		\n\t		movaps		%%xmm12,%%xmm10		\n\t"\
		"shufpd	$0,	%%xmm5,%%xmm4		\n\t		shufpd	$0,	%%xmm13,%%xmm12		\n\t"\
		"shufpd	$3,	%%xmm5,%%xmm2		\n\t		shufpd	$3,	%%xmm13,%%xmm10		\n\t"\
		"movaps		%%xmm4,%%xmm5		\n\t		movaps		%%xmm12,%%xmm13		\n\t"\
		"movaps		%%xmm2,%%xmm3		\n\t		movaps		%%xmm10,%%xmm11		\n\t"\
		"mulpd		%%xmm0,%%xmm4		\n\t		mulpd		%%xmm8 ,%%xmm12		\n\t"\
		"mulpd		%%xmm1,%%xmm3		\n\t		mulpd		%%xmm9 ,%%xmm11		\n\t"\
		"mulpd		%%xmm0,%%xmm2		\n\t		mulpd		%%xmm8 ,%%xmm10		\n\t"\
		"mulpd		%%xmm1,%%xmm5		\n\t		mulpd		%%xmm9 ,%%xmm13		\n\t"\
		"subpd		%%xmm3,%%xmm4		\n\t		subpd		%%xmm11,%%xmm12		\n\t"\
		"addpd		%%xmm2,%%xmm5		\n\t		addpd		%%xmm10,%%xmm13		\n\t"\
		"movaps		%%xmm4,    (%%rdx)	\n\t		movaps		%%xmm12,0x20(%%rdx)	\n\t"\
		"movaps		%%xmm5,0x10(%%rdx)	\n\t		movaps		%%xmm13,0x30(%%rdx)	\n\t"\
		"/* Store larger of maxerr1,2: */	\n\t	movslq	%[__idx_incr],%%rdi		\n\t"\
		"maxpd		%%xmm14,%%xmm6		\n\t		addq	%%r10,%%rdi				\n\t"\
		"movaps		%%xmm6,-0x20(%%rcx)	\n\t		mov	%%edi, %[__idx_offset]	/* Store twice-incremented idx_offset */	\n\t"\
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
		: "cc","memory","rax","rbx","rcx","rdx","rsi","rdi","r10","r11","r12","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14"		/* Clobbered registers */\
	);\
	}

	/* Non-power-of-2-runlength Fermat-mod acyclic-transform/IBDWT carry macro.
	The array indices icycle0/1 declared int in caller but assumed to have been << 4 at time this macro called, thus can use as complex-address-offsets.
	*/
	#define SSE2_fermat_carry_norm_errcheck(Xdata,Xcy,Xnrt_bits,Xnrtm1,Xidx_offset,Xidx_incr,Xodd_radix,Xhalf_arr,Xsign_mask,Xadd1,Xadd2,Xicycle0,Xjcycle0)\
	{\
	__asm__ volatile (\
		"movslq	%[__idx_offset],%%rsi	\n\t"\
		"movslq %[__odd_radix],%%rdi	\n\t"\
		"movslq	%[__nrt_bits],%%rcx		\n\t"\
		"movslq %[__icycle0],%%r10		\n\t"\
		"movslq	%[__jcycle0],%%r11		\n\t"\
		"movq		%%rsi,%%rax			\n\t"\
		"shrq		$1,%%rax			\n\t"\
		"movq		%%rax,%%rbx			\n\t"\
		"andq		%[__nrtm1],%%rax	\n\t"\
		"shrq		%%cl,%%rbx			\n\t"\
		"shlq		$4,%%rax			\n\t"\
		"shlq		$4,%%rbx			\n\t"\
		"shlq		$4,%%rdi			\n\t"\
		"addq		%[__add1],%%rax		\n\t"\
		"addq		%[__add2],%%rbx		\n\t"\
		"movaps		(%%rax),%%xmm0		\n\t"\
		"movaps		(%%rbx),%%xmm1		\n\t"\
		"movq		%%rsi,%%rax			\n\t"\
		"movaps		%%xmm1,%%xmm2		\n\t"\
		"shufpd	$1,	%%xmm2,%%xmm2		\n\t"\
		"mulpd		%%xmm0,%%xmm1		\n\t"\
		"mulpd		%%xmm0,%%xmm2		\n\t"\
		"addq		$2,%%rax			\n\t"\
		"shrq		$1,%%rax			\n\t"\
		"movq		%%rax,%%rbx			\n\t"\
		"andq		%[__nrtm1],%%rax	\n\t"\
		"shrq		%%cl,%%rbx			\n\t"\
		"shlq		$4,%%rax			\n\t"\
		"shlq		$4,%%rbx			\n\t"\
		"addq		%[__add1],%%rax		\n\t"\
		"addq		%[__add2],%%rbx		\n\t"\
		"movaps		(%%rax),%%xmm0		\n\t"\
		"movaps		(%%rbx),%%xmm3		\n\t"\
		"movq		%%rsi,%%rax			\n\t"\
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
		"movq		%[__half_arr],%%rcx	\n\t"\
		"movq		%[__data],%%rdx		\n\t"\
		"movaps		     (%%rdx),%%xmm4	\n\t"\
		"movaps		 0x10(%%rdx),%%xmm2	\n\t"\
		"addq		%%r10,%%rcx	\n\t"\
		"movaps	(%%rcx,%%rdi),%%xmm5	\n\t"\
		"subq		%%r10,%%rcx	\n\t"\
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
		"movq		%[__cy],%%rbx		\n\t"\
		"movaps		%%xmm4,%%xmm5		\n\t"\
		"shufpd	$0,	%%xmm2,%%xmm4		\n\t"\
		"shufpd	$3,	%%xmm2,%%xmm5		\n\t"\
		"addpd		     (%%rbx),%%xmm4	\n\t"\
		"movaps		-0x20(%%rcx),%%xmm6	\n\t"\
		"movaps		-0x10(%%rcx),%%xmm7	\n\t"\
		"addq	   %%r10,%%rcx	\n\t"\
		"movaps		%%xmm4,%%xmm2		\n\t"\
		"shlq	   $1,%%rdi				\n\t"\
		"addpd		%%xmm7,%%xmm4		\n\t"\
		"subpd		%%xmm7,%%xmm4		\n\t"\
		"movq		%[__sign_mask],%%rax\n\t"\
		"subpd		%%xmm4,%%xmm2		\n\t"\
		"andpd		(%%rax),%%xmm2		\n\t"\
		"maxpd		%%xmm6,%%xmm2		\n\t"\
		"movaps		%%xmm2,%%xmm6		\n\t"\
		"addq		%%rdi,%%rcx			\n\t"\
		"shrq		$1,%%rdi			\n\t"\
		"movaps		%%xmm4,%%xmm2		\n\t"\
		"mulpd  (%%rcx,%%rdi),%%xmm2	\n\t"\
		"addpd		%%xmm7,%%xmm2		\n\t"\
		"subpd		%%xmm7,%%xmm2		\n\t"\
		"movaps		%%xmm2,%%xmm3		\n\t"\
		"mulpd		(%%rcx),%%xmm3		\n\t"\
		"subq		%%r10,%%rcx	\n\t"\
		"subpd		%%xmm3,%%xmm4		\n\t"\
		"addpd		%%xmm2,%%xmm5		\n\t"\
		"movaps		%%xmm5,%%xmm2		\n\t"\
		"addpd		%%xmm7,%%xmm5		\n\t"\
		"subpd		%%xmm7,%%xmm5		\n\t"\
		"subpd		%%xmm5,%%xmm2		\n\t"\
		"andpd		(%%rax),%%xmm2		\n\t"\
		"maxpd		%%xmm6,%%xmm2		\n\t"\
		"movaps		%%xmm2,%%xmm6		\n\t"\
		"movaps		%%xmm5,%%xmm2		\n\t"\
		"addq		%%r11,%%rcx	\n\t"\
		"mulpd  (%%rcx,%%rdi),%%xmm2	\n\t"\
		"addpd		%%xmm7,%%xmm2		\n\t"\
		"subpd		%%xmm7,%%xmm2		\n\t"\
		"shlq		$1,%%rdi			\n\t"\
		"movaps		%%xmm2,%%xmm3		\n\t"\
		"mulpd		(%%rcx),%%xmm3		\n\t"\
		"subq		%%r11,%%rcx	\n\t"\
		"subpd		%%xmm3,%%xmm5		\n\t"\
		"movaps		%%xmm2,(%%rbx)		\n\t"\
		"movaps		%%xmm4,%%xmm2		\n\t"\
		"shufpd	$0,	%%xmm5,%%xmm4		\n\t"\
		"shufpd	$3,	%%xmm5,%%xmm2		\n\t"\
		"movaps		%%xmm4,%%xmm5		\n\t"\
		"movaps		%%xmm2,%%xmm3		\n\t"\
		"subq		%%rdi,%%rcx			\n\t"\
		"mulpd		%%xmm0,%%xmm4		\n\t"\
		"mulpd		%%xmm1,%%xmm3		\n\t"\
		"movaps	 	%%xmm6,-0x20(%%rcx)	\n\t"\
		"addq		%%r10,%%rcx	\n\t"\
		"mulpd		%%xmm0,%%xmm2		\n\t"\
		"mulpd		%%xmm1,%%xmm5		\n\t"\
		"movaps		(%%rcx),%%xmm0		\n\t"\
		"subpd		%%xmm3,%%xmm4		\n\t"\
		"addpd		%%xmm2,%%xmm5		\n\t"\
		"mulpd		%%xmm0,%%xmm4		\n\t"\
		"mulpd		%%xmm0,%%xmm5		\n\t"\
		"movaps		%%xmm4,    (%%rdx)	\n\t"\
		"movaps		%%xmm5,0x10(%%rdx)	\n\t"\
		"addq	%[__idx_incr],%%rsi		\n\t"\
		"mov	%%esi, %[__idx_offset]	/* Store incremented idx_offset */	\n\t"\
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
	,	[__icycle0]		"m" (Xicycle0)\
	,	[__jcycle0]		"m" (Xjcycle0)\
		: "cc","memory","rax","rbx","rcx","rdx","rsi","rdi","r10","r11","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"   /* Clobbered registers */\
	);\
	}

	/* Same non-power-of-2-transform deal as above, but use xmm8-15 to process 2 sets of carries side-by-side.
	Data/Carry #2 assumed offset by +0x20/0x10 from #1 (which are accessed via the [__data/__cy] pointers, resp.)
	i/jcycle0 and i/jcycle1 are the address offsets needed for IBDWT array indexing for the 2 resp. carries.
	*/
	#define SSE2_fermat_carry_norm_errcheck_X2(Xdata,Xcy,Xnrt_bits,Xnrtm1,Xidx_offset,Xidx_incr,Xodd_radix,Xhalf_arr,Xsign_mask,Xadd1,Xadd2,Xicycle0,Xjcycle0,Xicycle1,Xjcycle1)\
	{\
	__asm__ volatile (\
		"/* lcol -> rcol index analogs: [rsi,rax,rbx] -> [r10,r11,r12], [rcx,rdx,rdi] shared */\n\t"\
		"movslq	%[__idx_offset],%%rsi	\n\t		movslq	%[__idx_incr],%%r10		\n\t"\
		"movslq		%[__nrt_bits],%%rcx	\n\t		addq	%%rsi,%%r10				\n\t"\
		"movslq		%[__nrtm1],%%rdi	/* r10 contains idx_offset2, i.e. is the rcol-analog of rsi in lcol: */\n\t"\
		"movq		%%rsi,%%rax			\n\t		movq		%%r10,%%r11			\n\t"\
		"shrq		$1,%%rax			\n\t		shrq		$1,%%r11			\n\t"\
		"movq		%%rax,%%rbx			\n\t		movq		%%r11,%%r12			\n\t"\
		"andq		%%rdi,%%rax			\n\t		andq		%%rdi,%%r11			\n\t"\
		"shrq		%%cl,%%rbx			\n\t		shrq		%%cl,%%r12			\n\t"\
		"shlq		$4,%%rax			\n\t		shlq		$4,%%r11			\n\t"\
		"shlq		$4,%%rbx			\n\t		shlq		$4,%%r12			\n\t"\
		"addq		%[__add1],%%rax		\n\t		addq		%[__add1],%%r11		\n\t"\
		"addq		%[__add2],%%rbx		\n\t		addq		%[__add2],%%r12		\n\t"\
		"movaps		(%%rax),%%xmm0		\n\t		movaps		(%%r11),%%xmm8 		\n\t"\
		"movaps		(%%rbx),%%xmm1		\n\t		movaps		(%%r12),%%xmm9 		\n\t"\
		"movq		%%rsi,%%rax			\n\t		movq		%%r10,%%r11			\n\t"\
		"movaps		%%xmm1,%%xmm2		\n\t		movaps		%%xmm9 ,%%xmm10		\n\t"\
		"shufpd	$1,	%%xmm2,%%xmm2		\n\t		shufpd	$1,	%%xmm10,%%xmm10		\n\t"\
		"mulpd		%%xmm0,%%xmm1		\n\t		mulpd		%%xmm8 ,%%xmm9 		\n\t"\
		"mulpd		%%xmm0,%%xmm2		\n\t		mulpd		%%xmm8 ,%%xmm10		\n\t"\
		"addq		$2,%%rax			\n\t		addq		$2,%%r11			\n\t"\
		"shrq		$1,%%rax			\n\t		shrq		$1,%%r11			\n\t"\
		"movq		%%rax,%%rbx			\n\t		movq		%%r11,%%r12			\n\t"\
		"andq		%%rdi,%%rax			\n\t		andq		%%rdi,%%r11			\n\t"\
		"shrq		%%cl,%%rbx			\n\t		shrq		%%cl,%%r12			\n\t"\
		"shlq		$4,%%rax			\n\t		shlq		$4,%%r11			\n\t"\
		"shlq		$4,%%rbx			\n\t		shlq		$4,%%r12			\n\t"\
		"addq		%[__add1],%%rax		\n\t		addq		%[__add1],%%r11		\n\t"\
		"addq		%[__add2],%%rbx		\n\t		addq		%[__add2],%%r12		\n\t"\
		"movaps		(%%rax),%%xmm0		\n\t		movaps		(%%r11),%%xmm8 		\n\t"\
		"movaps		(%%rbx),%%xmm3		\n\t		movaps		(%%r12),%%xmm11		\n\t"\
		"movq		%%rsi,%%rax			\n\t		movq		%%r10,%%r11			\n\t"\
		"movaps		%%xmm3,%%xmm4		\n\t		movaps		%%xmm11,%%xmm12		\n\t"\
		"shufpd	$1,	%%xmm4,%%xmm4		\n\t		shufpd	$1,	%%xmm12,%%xmm12		\n\t"\
		"mulpd		%%xmm0,%%xmm3		\n\t		mulpd		%%xmm8 ,%%xmm11		\n\t"\
		"mulpd		%%xmm0,%%xmm4		\n\t		mulpd		%%xmm8 ,%%xmm12		\n\t"\
		"movaps		%%xmm1,%%xmm0		\n\t		movaps		%%xmm9 ,%%xmm8 		\n\t"\
		"unpcklpd	%%xmm3,%%xmm0		\n\t		unpcklpd	%%xmm11,%%xmm8 		\n\t"\
		"unpckhpd	%%xmm3,%%xmm1		\n\t		unpckhpd	%%xmm11,%%xmm9 		\n\t"\
		"subpd		%%xmm1,%%xmm0		\n\t		subpd		%%xmm9 ,%%xmm8 		\n\t"\
		"movaps		%%xmm2,%%xmm1		\n\t		movaps		%%xmm10,%%xmm9 		\n\t"\
		"unpcklpd	%%xmm4,%%xmm1		\n\t		unpcklpd	%%xmm12,%%xmm9 		\n\t"\
		"unpckhpd	%%xmm4,%%xmm2		\n\t		unpckhpd	%%xmm12,%%xmm10		\n\t"\
		"addpd		%%xmm2,%%xmm1		\n\t		addpd		%%xmm10,%%xmm9 		\n\t"\
		"/* Store twice-incremented idx_offset to free up registers rdi,r10,r11: */	\n\t"\
		"											movslq	%[__idx_incr],%%rdi		\n\t"\
		"											addq	%%r10,%%rdi				\n\t"\
		"											mov	%%edi, %[__idx_offset]		\n\t"\
		"movslq	%[__icycle0],%%r8 		\n\t		movslq	%[__icycle1],%%r10		\n\t"\
		"movslq	%[__jcycle0],%%r9 		\n\t		movslq	%[__jcycle1],%%r11		\n\t"\
		"movslq	%[__odd_radix],%%rdi												\n\t"\
		"movq		%[__half_arr],%%rcx	/* Need separate rcol copy of rcx below */	\n\t"\
		"movq		%[__data],%%rdx		/* rdx shared, offset +0x20 in rcol: */		\n\t"\
		"shlq		$4,%%rdi			\n\t		movq		%%rcx,%%r12	/* rcol-copy for incr/decr: */\n\t"\
		"movaps		    (%%rdx),%%xmm4	\n\t		movaps		0x20(%%rdx),%%xmm12	\n\t"\
		"movaps		0x10(%%rdx),%%xmm2	\n\t		movaps		0x30(%%rdx),%%xmm10	\n\t"\
		"addq		%%r8 ,%%rcx			\n\t		addq		%%r10,%%r12			\n\t"\
		"movaps	(%%rcx,%%rdi),%%xmm5	\n\t		movaps	(%%r12,%%rdi),%%xmm13	\n\t"\
		"subq		%%r8 ,%%rcx			\n\t		subq		%%r10,%%r12	/* rcx == r12 again */\n\t"\
		"mulpd		%%xmm5,%%xmm4		\n\t		mulpd		%%xmm13,%%xmm12		\n\t"\
		"mulpd		%%xmm5,%%xmm2		\n\t		mulpd		%%xmm13,%%xmm10		\n\t"\
		"movaps		%%xmm4,%%xmm5		\n\t		movaps		%%xmm12,%%xmm13		\n\t"\
		"movaps		%%xmm2,%%xmm3		\n\t		movaps		%%xmm10,%%xmm11		\n\t"\
		"mulpd		%%xmm0,%%xmm4		\n\t		mulpd		%%xmm8 ,%%xmm12		\n\t"\
		"mulpd		%%xmm1,%%xmm3		\n\t		mulpd		%%xmm9 ,%%xmm11		\n\t"\
		"mulpd		%%xmm0,%%xmm2		\n\t		mulpd		%%xmm8 ,%%xmm10		\n\t"\
		"mulpd		%%xmm1,%%xmm5		\n\t		mulpd		%%xmm9 ,%%xmm13		\n\t"\
		"addpd		%%xmm3,%%xmm4		\n\t		addpd		%%xmm11,%%xmm12		\n\t"\
		"subpd		%%xmm5,%%xmm2		\n\t		subpd		%%xmm13,%%xmm10		\n\t"\
		"movq		%[__cy],%%rbx	/* rbx -> rbx+0x10 (carry offset only half of data-offset) in rcol, shared from here */	\n\t"\
		"movaps		%%xmm4,%%xmm5		\n\t		movaps		%%xmm12,%%xmm13		\n\t"\
		"shufpd	$0,	%%xmm2,%%xmm4		\n\t		shufpd	$0,	%%xmm10,%%xmm12		\n\t"\
		"shufpd	$3,	%%xmm2,%%xmm5		\n\t		shufpd	$3,	%%xmm10,%%xmm13		\n\t"\
		"addpd		     (%%rbx),%%xmm4	\n\t		addpd		 0x10(%%rbx),%%xmm12\n\t"\
		"movaps		-0x20(%%rcx),%%xmm6	/* maxerr, will make rcol-copy below, re-merge at end */\n\t"\
		"addq		%%r8 ,%%rcx			\n\t		addq		%%r10,%%r12			\n\t"\
		"movaps		%%xmm4,%%xmm2		\n\t		movaps		%%xmm12,%%xmm10		\n\t"\
		"shlq		$1,%%rdi			\n\t		movaps		%%xmm6,%%xmm14	/* rcol-copy of maxerr */\n\t"\
		"roundpd	$0,%%xmm4,%%xmm4	\n\t		roundpd	$0,%%xmm12,%%xmm12		\n\t"\
		"movq		%[__sign_mask],%%rax/* rax shared between lcol/rcol from here */\n\t"\
		"subpd		%%xmm4,%%xmm2		\n\t		subpd		%%xmm12,%%xmm10		\n\t"\
		"andpd		     (%%rax),%%xmm2	\n\t		andpd		     (%%rax),%%xmm10	\n\t"\
		"maxpd		%%xmm6,%%xmm2		\n\t		maxpd		%%xmm14,%%xmm10		\n\t"\
		"movaps		%%xmm2,%%xmm6		\n\t		movaps		%%xmm10,%%xmm14		\n\t"\
		"addq		%%rdi,%%rcx			\n\t		addq		%%rdi,%%r12			\n\t"\
		"shrq		$1,%%rdi														\n\t"\
		"movaps		%%xmm4,%%xmm2		\n\t		movaps		%%xmm12,%%xmm10		\n\t"\
		"mulpd	(%%rcx,%%rdi),%%xmm2	\n\t		mulpd	(%%r12,%%rdi),%%xmm10	\n\t"\
		"roundpd	$0,%%xmm2,%%xmm2	\n\t		roundpd	$0,%%xmm10,%%xmm10		\n\t"\
		"movaps		%%xmm2,%%xmm3		\n\t		movaps		%%xmm10,%%xmm11		\n\t"\
		"mulpd		    (%%rcx),%%xmm3	\n\t		mulpd		    (%%r12),%%xmm11	\n\t"\
		"subq		%%r8 ,%%rcx			\n\t		subq		%%r10,%%r12			\n\t"\
		"subpd		%%xmm3,%%xmm4		\n\t		subpd		%%xmm11,%%xmm12		\n\t"\
		"addpd		%%xmm2,%%xmm5		\n\t		addpd		%%xmm10,%%xmm13		\n\t"\
		"movaps		%%xmm5,%%xmm2		\n\t		movaps		%%xmm13,%%xmm10		\n\t"\
		"roundpd	$0,%%xmm5,%%xmm5	\n\t		roundpd	$0,%%xmm13,%%xmm13		\n\t"\
		"subpd		%%xmm5,%%xmm2		\n\t		subpd		%%xmm13,%%xmm10		\n\t"\
		"andpd		     (%%rax),%%xmm2	\n\t		andpd		     (%%rax),%%xmm10	\n\t"\
		"maxpd		%%xmm6,%%xmm2		\n\t		maxpd		%%xmm14,%%xmm10		\n\t"\
		"movaps		%%xmm2,%%xmm6		\n\t		movaps		%%xmm10,%%xmm14		\n\t"\
		"movaps		%%xmm5,%%xmm2		\n\t		movaps		%%xmm13,%%xmm10		\n\t"\
		"addq		%%r9 ,%%rcx			\n\t		addq		%%r11,%%r12			\n\t"\
		"maxpd		%%xmm14,%%xmm6		/* Save larger of maxerr1,2: */	\n\t"\
		"mulpd	(%%rcx,%%rdi),%%xmm2	\n\t		mulpd	(%%r12,%%rdi),%%xmm10	\n\t"\
		"roundpd	$0,%%xmm2,%%xmm2	\n\t		roundpd	$0,%%xmm10,%%xmm10		\n\t"\
		"shlq		$1,%%rdi														\n\t"\
		"movaps		%%xmm2,%%xmm3		\n\t		movaps		%%xmm10,%%xmm11		\n\t"\
		"mulpd		     (%%rcx),%%xmm3	\n\t		mulpd		     (%%r12),%%xmm11	\n\t"\
		"subq		%%r9 ,%%rcx			\n\t		subq		%%r11,%%r12			\n\t"\
		"subpd		%%xmm3,%%xmm5		\n\t		subpd		%%xmm11,%%xmm13		\n\t"\
		"movaps		%%xmm2,(%%rbx)		\n\t		movaps		%%xmm10,0x10(%%rbx)	\n\t"\
		"movaps		%%xmm4,%%xmm2		\n\t		movaps		%%xmm12,%%xmm10		\n\t"\
		"shufpd	$0,	%%xmm5,%%xmm4		\n\t		shufpd	$0,	%%xmm13,%%xmm12		\n\t"\
		"shufpd	$3,	%%xmm5,%%xmm2		\n\t		shufpd	$3,	%%xmm13,%%xmm10		\n\t"\
		"movaps		%%xmm4,%%xmm5		\n\t		movaps		%%xmm12,%%xmm13		\n\t"\
		"movaps		%%xmm2,%%xmm3		\n\t		movaps		%%xmm10,%%xmm11		\n\t"\
		"subq		%%rdi,%%rcx			\n\t		subq		%%rdi,%%r12	/* rcx == r12 again */\n\t"\
		"mulpd		%%xmm0,%%xmm4		\n\t		mulpd		%%xmm8 ,%%xmm12		\n\t"\
		"mulpd		%%xmm1,%%xmm3		\n\t		mulpd		%%xmm9 ,%%xmm11		\n\t"\
		"movaps		%%xmm6,-0x20(%%rcx)	/* Store maxerr: */	\n\t"\
		"addq		%%r8 ,%%rcx			\n\t		addq		%%r10,%%r12			\n\t"\
		"mulpd		%%xmm0,%%xmm2		\n\t		mulpd		%%xmm8 ,%%xmm10		\n\t"\
		"mulpd		%%xmm1,%%xmm5		\n\t		mulpd		%%xmm9 ,%%xmm13		\n\t"\
		"movaps		(%%rcx),%%xmm0		\n\t		movaps		(%%r12),%%xmm8 		\n\t"\
		"subpd		%%xmm3,%%xmm4		\n\t		subpd		%%xmm11,%%xmm12		\n\t"\
		"addpd		%%xmm2,%%xmm5		\n\t		addpd		%%xmm10,%%xmm13		\n\t"\
		"mulpd		%%xmm0,%%xmm4		\n\t		mulpd		%%xmm8 ,%%xmm12		\n\t"\
		"mulpd		%%xmm0,%%xmm5		\n\t		mulpd		%%xmm8 ,%%xmm13		\n\t"\
		"movaps		%%xmm4,    (%%rdx)	\n\t		movaps		%%xmm12,0x20(%%rdx)	\n\t"\
		"movaps		%%xmm5,0x10(%%rdx)	\n\t		movaps		%%xmm13,0x30(%%rdx)	\n\t"\
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
	,	[__icycle0]		"m" (Xicycle0)\
	,	[__jcycle0]		"m" (Xjcycle0)\
	,	[__icycle1]		"m" (Xicycle1)\
	,	[__jcycle1]		"m" (Xjcycle1)\
		: "cc","memory","rax","rbx","rcx","rdx","rsi","rdi","r8","r9","r10","r11","r12","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14"		/* Clobbered registers */\
	);\
	}

	/*************************************************************/
	/**************** MERSENNE-MOD CARRY MACROS ******************/
	/*************************************************************/

	/********* Packed 32-bit-int version of SSE2_cmplx_carry_norm_pow2_errcheck0_2x:***********/
	#define SSE2_cmplx_carry_norm_pow2_errcheck0_2B(Xdata,XwtA,XwtB,XwtC,XcyA,XcyB,Xbjmod_0,Xhalf_arr,Xi,Xn_minus_silp1,Xn_minus_sil,Xsign_mask,Xsinwt,Xsinwtm1,Xsse_bw,Xsse_nm1,Xsse_sw)\
	{\
	__asm__ volatile (\
	"/***************Unpack the data:*************************/\n\t"\
		"movq	%[__data]	,%%rax	\n\t"\
		"movaps		    (%%rax)	,%%xmm1	\n\t	movaps		0x40(%%rax)	,%%xmm5	\n\t"\
		"movaps		    (%%rax)	,%%xmm2	\n\t	movaps		0x40(%%rax)	,%%xmm6	\n\t"\
		"unpcklpd	0x20(%%rax)	,%%xmm1	\n\t	unpcklpd	0x60(%%rax)	,%%xmm5	\n\t"\
		"unpckhpd	0x20(%%rax)	,%%xmm2	\n\t	unpckhpd	0x60(%%rax)	,%%xmm6	\n\t"\
		"movaps		%%xmm2, 0x20(%%rax)	\n\t	movaps		%%xmm6, 0x60(%%rax)	\n\t"\
		"\n\t"\
		"movaps		0x10(%%rax)	,%%xmm2	\n\t	movaps		0x50(%%rax)	,%%xmm6	\n\t"\
		"movaps		0x10(%%rax)	,%%xmm3	\n\t	movaps		0x50(%%rax)	,%%xmm7	\n\t"\
		"unpcklpd	0x30(%%rax)	,%%xmm2	\n\t	unpcklpd	0x70(%%rax)	,%%xmm6	\n\t"\
		"unpckhpd	0x30(%%rax)	,%%xmm3	\n\t	unpckhpd	0x70(%%rax)	,%%xmm7	\n\t"\
		"movaps		%%xmm2, 0x10(%%rax)	\n\t	movaps		%%xmm6, 0x50(%%rax)	\n\t"\
		"movaps		%%xmm3, 0x30(%%rax)	\n\t	movaps		%%xmm7, 0x70(%%rax)	\n\t"\
	"/**********************************************/\n\t"\
	"/*          Real      parts                   */\n\t"\
	"/**********************************************/\n\t"\
		"movq	%[__bjmod_0],	%%rax	\n\t"\
		"movaps		(%%rax)	,	%%xmm0	\n\t"\
		"movq	%[__sse_sw]	,	%%rbx	\n\t"\
		"movaps		(%%rbx)	,	%%xmm7	\n\t"\
		"psubd		%%xmm0	,	%%xmm7	\n\t"\
		"movmskps	%%xmm7	,	%%rsi	\n\t"\
		"movslq	%[__i]	,%%rcx			\n\t"\
		"andq	$0xfffffffffffffffe,%%rsi		\n\t"\
		"addq	%%rcx	,%%rsi			\n\t"\
		"shlq	$24		,%%rsi			\n\t"\
		"movaps		%%xmm0	,%%xmm7		\n\t"\
		"movslq	%[__n_minus_sil],%%rcx	\n\t"\
		"movd	%%rcx	,%%xmm2			\n\t"\
		"pshufd	$0,	%%xmm2	,%%xmm2		\n\t"\
		"psubd		%%xmm0	,%%xmm2		\n\t"\
		"movmskps	%%xmm2	,%%rcx		\n\t"\
		"shlq	$16		,%%rcx			\n\t"\
		"addq	%%rcx	,%%rsi			\n\t"\
		"movq	%[__sinwt]	,%%rdx		\n\t"\
		"movd	%%rdx	,%%xmm3			\n\t"\
		"pshufd	$0,	%%xmm3	,%%xmm3		\n\t"\
		"psubd		%%xmm3	,%%xmm7		\n\t"\
		"movmskps	%%xmm7	,%%rdx		\n\t"\
		"shlq	$8		,%%rdx			\n\t"\
		"addq	%%rdx	,%%rsi			\n\t"\
		"\n\t"\
		"movq	%[__half_arr]	,%%rax		\n\t"\
		"movq	%[__wtA]	,%%rbx		\n\t"\
		"movq	%[__wtB]	,%%rcx		\n\t"\
		"movaps		     (%%rbx),%%xmm2	\n\t	movaps		 0x10(%%rbx),%%xmm6	\n\t"\
		"movhpd		     (%%rcx),%%xmm3	\n\t	movhpd		-0x10(%%rcx),%%xmm7	\n\t"\
		"movlpd		 0x08(%%rcx),%%xmm3	\n\t	movlpd		-0x08(%%rcx),%%xmm7	\n\t"\
		"\n\t"\
		"subq	$0x20	,%%rcx			\n\t"\
		"movq	%%rcx	,%[__wtB]		\n\t"\
		"\n\t"\
		"movq	%%rsi,	%%rdi	\n\t	movq	%%rsi,	%%rbx	\n\t"\
		"shrq	$12,	%%rdi	\n\t	shrq	$14,	%%rbx	\n\t"\
		"andq	$0x0000000000000030	,%%rdi		\n\t	andq	$0x0000000000000030	,%%rbx		\n\t"\
		"movq	%%rsi,	%%rdx	\n\t	movq	%%rsi,	%%rcx	\n\t"\
		"shrq	$4 ,	%%rdx	\n\t	shrq	$6 ,	%%rcx	\n\t"\
		"andq	$0x0000000000000030	,%%rdx		\n\t	andq	$0x0000000000000030	,%%rcx		\n\t"\
		"addq	%%rax	,%%rdi			\n\t	addq	%%rax	,%%rbx			\n\t"\
		"addq	%%rax	,%%rdx			\n\t	addq	%%rax	,%%rcx			\n\t"\
		"\n\t"\
		"mulpd	 	0x100(%%rax),%%xmm2	\n\t	mulpd	 	0x100(%%rax),%%xmm6	\n\t"\
		"mulpd	 	0x110(%%rax),%%xmm3	\n\t	mulpd	 	0x110(%%rax),%%xmm7	\n\t"\
		"mulpd	 	     (%%rdi),%%xmm2	\n\t	mulpd	 	     (%%rbx),%%xmm6	\n\t"\
		"mulpd	 	0x040(%%rdx),%%xmm3	\n\t	mulpd	 	0x040(%%rcx),%%xmm7	\n\t"\
		"\n\t"\
		"movq	%[__cyA]	,%%rcx		\n\t	movq	%[__cyB]	,%%rdx		\n\t"\
		"mulpd		%%xmm3	,%%xmm1		\n\t	mulpd		%%xmm7	,%%xmm5		\n\t"\
		"addpd		(%%rcx)	,%%xmm1		\n\t	addpd		(%%rdx)	,%%xmm5		\n\t"\
		"movaps		%%xmm1	,%%xmm3		\n\t	movaps		%%xmm5	,%%xmm7		\n\t"\
		"roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd		$0,%%xmm7,%%xmm7	\n\t"\
		"\n\t"\
		"movq	%[__sign_mask],%%rbx	\n\t"\
		"subpd			%%xmm3	,%%xmm1	\n\t	subpd			%%xmm7	,%%xmm5	\n\t"\
		"andpd			(%%rbx)	,%%xmm1	\n\t	andpd			(%%rbx)	,%%xmm5	\n\t"\
		"maxpd			%%xmm5	,%%xmm1	\n\t"\
		"maxpd		-0x20(%%rax),%%xmm1	\n\t"\
		"movaps		%%xmm1,-0x20(%%rax)	\n\t"\
		"\n\t"\
		"movq	%%rsi,	%%rdi	\n\t	movq	%%rsi,	%%rbx	\n\t"\
		"shrq	$20,	%%rdi	\n\t	shrq	$22,	%%rbx	\n\t"\
		"andq	$0x0000000000000030	,%%rdi		\n\t	andq	$0x0000000000000030	,%%rbx		\n\t"\
		"addq	%%rax	,%%rdi			\n\t	addq	%%rax	,%%rbx			\n\t"\
		"movaps		%%xmm3	,%%xmm1		\n\t	movaps		%%xmm7	,%%xmm5		\n\t"\
		"mulpd	 	0xc0(%%rdi)	,%%xmm3	\n\t	mulpd	 	0xc0(%%rbx)	,%%xmm7	\n\t"\
		"roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd		$0,%%xmm7,%%xmm7	\n\t"\
		"movaps		%%xmm3	,(%%rcx)	\n\t	movaps		%%xmm7	,(%%rdx)	\n\t"\
		"\n\t"\
		"\n\t"\
		"movq	%[__data]	,%%rax		\n\t"\
		"mulpd	 0x80(%%rdi)	,%%xmm3	\n\t	mulpd	 0x80(%%rbx)	,%%xmm7	\n\t"\
		"subpd		%%xmm3	,	%%xmm1	\n\t	subpd		%%xmm7	,	%%xmm5	\n\t"\
		"mulpd		%%xmm2	,	%%xmm1	\n\t	mulpd		%%xmm6	,	%%xmm5	\n\t"\
		"movaps		%%xmm1	,    (%%rax)\n\t	movaps		%%xmm5	,0x40(%%rax)\n\t"\
		"\n\t"\
		"\n\t"\
		"movq	%[__sse_bw]	,%%rax		\n\t"\
		"movq	%[__sse_nm1]	,%%rbx	\n\t"\
		"paddd		(%%rax)	,%%xmm0		\n\t"\
		"pand		(%%rbx)	,%%xmm0		\n\t"\
	"/**********************************************/\n\t"\
	"/*          Imaginary parts               */\n\t"\
	"/**********************************************/\n\t"\
		"movq	%[__sse_sw]	,%%rdx		\n\t"\
		"movaps	(%%rdx)	,%%xmm1			\n\t"\
		"psubd	%%xmm0	,%%xmm1			\n\t"\
		"movmskps	%%xmm1	,%%rsi		\n\t"\
		"shlq	$24	,%%rsi				\n\t"\
		"movaps	%%xmm0	,%%xmm1			\n\t"\
		"movslq	%[__n_minus_silp1],%%rcx\n\t"\
		"movd	%%rcx	,%%xmm2			\n\t"\
		"pshufd	$0	,%%xmm2	,%%xmm2		\n\t"\
		"psubd	%%xmm0	,%%xmm2			\n\t"\
		"movmskps	%%xmm2	,%%rcx		\n\t"\
		"shlq	$16	,%%rcx				\n\t"\
		"addq	%%rcx	,%%rsi			\n\t"\
		"movq	%[__sinwtm1]	,%%rdx	\n\t"\
		"movd	%%rdx	,%%xmm3			\n\t"\
		"pshufd	$0	,%%xmm3	,%%xmm3		\n\t"\
		"psubd	%%xmm3	,%%xmm1			\n\t"\
		"movmskps	%%xmm1	,%%rdx		\n\t"\
		"shlq	$8	,%%rdx				\n\t"\
		"addq	%%rdx	,%%rsi			\n\t"\
		"movq	%[__data]	,%%rax		\n\t"\
		"movaps	 0x10(%%rax)	,%%xmm1	\n\t"\
		"movaps	 0x50(%%rax)	,%%xmm5	\n\t"\
		"\n\t"\
		"movq	%[__half_arr]	,%%rax	\n\t"\
		"movq	%[__wtA]	,%%rbx		\n\t"\
		"movq	%[__wtC]	,%%rcx		\n\t"\
		"\n\t"\
		"movaps	     (%%rbx)	,%%xmm2	\n\t	movaps	 0x10(%%rbx)	,%%xmm6	\n\t"\
		"movhpd	     (%%rcx)	,%%xmm3	\n\t	movhpd	-0x10(%%rcx)	,%%xmm7	\n\t"\
		"movlpd	 0x08(%%rcx)	,%%xmm3	\n\t	movlpd	-0x08(%%rcx)	,%%xmm7	\n\t"\
		"\n\t"\
		"addq	$0x20	,%%rbx			\n\t"\
		"subq	$0x20	,%%rcx			\n\t"\
		"movq	%%rbx	,%[__wtA]		\n\t"\
		"movq	%%rcx	,%[__wtC]		\n\t"\
		"movq	%%rsi,	%%rdi	\n\t	movq	%%rsi,	%%rbx	\n\t"\
		"shrq	$12,	%%rdi	\n\t	shrq	$14,	%%rbx	\n\t"\
		"andq	$0x0000000000000030	,%%rdi		\n\t	andq	$0x0000000000000030	,%%rbx		\n\t"\
		"movq	%%rsi,	%%rdx	\n\t	movq	%%rsi,	%%rcx	\n\t"\
		"shrq	$4 ,	%%rdx	\n\t	shrq	$6 ,	%%rcx	\n\t"\
		"andq	$0x0000000000000030	,%%rdx		\n\t	andq	$0x0000000000000030	,%%rcx		\n\t"\
		"addq	%%rax	,%%rdi			\n\t	addq	%%rax	,%%rbx			\n\t"\
		"addq	%%rax	,%%rdx			\n\t	addq	%%rax	,%%rcx			\n\t"\
		"\n\t"\
		"mulpd	 0x120(%%rax)	,%%xmm2	\n\t	mulpd	 0x120(%%rax)	,%%xmm6	\n\t"\
		"mulpd	 0x130(%%rax)	,%%xmm3	\n\t	mulpd	 0x130(%%rax)	,%%xmm7	\n\t"\
		"mulpd	      (%%rdi)	,%%xmm2	\n\t	mulpd	      (%%rbx)	,%%xmm6	\n\t"\
		"mulpd	 0x040(%%rdx)	,%%xmm3	\n\t	mulpd	 0x040(%%rcx)	,%%xmm7	\n\t"\
		"\n\t"\
		"movq	%[__cyA]	,%%rcx		\n\t	movq	%[__cyB]	,%%rdx		\n\t"\
		"mulpd	%%xmm3	,%%xmm1			\n\t	mulpd	%%xmm7	,%%xmm5			\n\t"\
		"addpd	(%%rcx)	,%%xmm1			\n\t	addpd	(%%rdx)	,%%xmm5			\n\t"\
		"movaps	%%xmm1	,%%xmm3			\n\t	movaps	%%xmm5	,%%xmm7			\n\t"\
		"roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd		$0,%%xmm7,%%xmm7	\n\t"\
		"\n\t"\
		"movq	%[__sign_mask],%%rbx	\n\t"\
		"subpd	%%xmm3	,%%xmm1			\n\t	subpd	%%xmm7	,%%xmm5			\n\t"\
		"andpd	(%%rbx)	,%%xmm1			\n\t	andpd	(%%rbx)	,%%xmm5			\n\t"\
		"maxpd	%%xmm5	,%%xmm1			\n\t"\
		"maxpd	-0x20(%%rax)	,%%xmm1	\n\t"\
		"movaps	%%xmm1	,-0x20(%%rax)	\n\t"\
		"\n\t"\
		"movq	%%rsi,	%%rdi	\n\t	movq	%%rsi,	%%rbx	\n\t"\
		"shrq	$20,	%%rdi	\n\t	shrq	$22,	%%rbx	\n\t"\
		"andq	$0x0000000000000030	,%%rdi		\n\t	andq	$0x0000000000000030	,%%rbx		\n\t"\
		"addq	%%rax	,%%rdi			\n\t	addq	%%rax	,%%rbx			\n\t"\
		"movaps	%%xmm3	,%%xmm1			\n\t	movaps	%%xmm7	,%%xmm5			\n\t"\
		"mulpd	 0xc0(%%rdi)	,%%xmm3	\n\t	mulpd	 0xc0(%%rbx)	,%%xmm7	\n\t"\
		"roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd		$0,%%xmm7,%%xmm7	\n\t"\
		"movaps	%%xmm3	,(%%rcx)		\n\t	movaps	%%xmm7	,(%%rdx)		\n\t"\
		"\n\t"\
		"\n\t"\
		"movq	%[__data]	,%%rax		\n\t"\
		"mulpd	 0x80(%%rdi)	,%%xmm3	\n\t	mulpd	 0x80(%%rbx)	,%%xmm7	\n\t"\
		"subpd	%%xmm3	,%%xmm1			\n\t	subpd	%%xmm7	,%%xmm5			\n\t"\
		"mulpd	%%xmm2	,%%xmm1			\n\t	mulpd	%%xmm6	,%%xmm5			\n\t"\
		"movaps	%%xmm1	, 0x10(%%rax)	\n\t	movaps	%%xmm5	, 0x50(%%rax)	\n\t"\
		"\n\t"\
		"\n\t"\
		"movq	%[__sse_bw]	,%%rax		\n\t"\
		"movq	%[__sse_nm1]	,%%rbx	\n\t"\
		"paddd	(%%rax)	,%%xmm0			\n\t"\
		"pand	(%%rbx)	,%%xmm0			\n\t"\
		"movq	%[__bjmod_0],%%rcx		\n\t"\
		"movaps	%%xmm0,(%%rcx)			\n\t"\
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
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","xmm0","xmm1","xmm2","xmm3","xmm5","xmm6","xmm7"	/* Clobbered registers - Use of rondpd frees up xmm4 */\
	);\
	}

	/********* Packed 32-bit-int version of SSE2_cmplx_carry_norm_pow2_errcheck1_2x:***********/
	#define SSE2_cmplx_carry_norm_pow2_errcheck1_2B(Xdata,XwtA,XwtB,XwtC,XcyA,XcyB,Xbjmod_0,Xhalf_arr,Xn_minus_silp1,Xn_minus_sil,Xsign_mask,Xsinwt,Xsinwtm1,Xsse_bw,Xsse_nm1,Xsse_sw)\
	{\
	__asm__ volatile (\
	"/***************Unpack the data:*************************/\n\t"\
		"movq	%[__data]	,%%rax	\n\t"\
		"movaps		    (%%rax)	,%%xmm1	\n\t	movaps		0x40(%%rax)	,%%xmm5	\n\t"\
		"movaps		    (%%rax)	,%%xmm2	\n\t	movaps		0x40(%%rax)	,%%xmm6	\n\t"\
		"unpcklpd	0x20(%%rax)	,%%xmm1	\n\t	unpcklpd	0x60(%%rax)	,%%xmm5	\n\t"\
		"unpckhpd	0x20(%%rax)	,%%xmm2	\n\t	unpckhpd	0x60(%%rax)	,%%xmm6	\n\t"\
		"movaps		%%xmm2, 0x20(%%rax)	\n\t	movaps		%%xmm6, 0x60(%%rax)	\n\t"\
		"\n\t"\
		"movaps		0x10(%%rax)	,%%xmm2	\n\t	movaps		0x50(%%rax)	,%%xmm6	\n\t"\
		"movaps		0x10(%%rax)	,%%xmm3	\n\t	movaps		0x50(%%rax)	,%%xmm7	\n\t"\
		"unpcklpd	0x30(%%rax)	,%%xmm2	\n\t	unpcklpd	0x70(%%rax)	,%%xmm6	\n\t"\
		"unpckhpd	0x30(%%rax)	,%%xmm3	\n\t	unpckhpd	0x70(%%rax)	,%%xmm7	\n\t"\
		"movaps		%%xmm2, 0x10(%%rax)	\n\t	movaps		%%xmm6, 0x50(%%rax)	\n\t"\
		"movaps		%%xmm3, 0x30(%%rax)	\n\t	movaps		%%xmm7, 0x70(%%rax)	\n\t"\
	"/**********************************************/\n\t"\
	"/*          Real      parts                   */\n\t"\
	"/**********************************************/\n\t"\
		"movq	%[__bjmod_0],	%%rax	\n\t"\
		"movaps		(%%rax)	,	%%xmm0	\n\t"\
		"movq	%[__sse_sw]	,	%%rbx	\n\t"\
		"movaps		(%%rbx)	,	%%xmm7	\n\t"\
		"psubd		%%xmm0	,	%%xmm7	\n\t"\
		"movmskps	%%xmm7	,	%%rsi	\n\t"\
		"\n\t"\
		"shlq	$24		,%%rsi			\n\t"\
		"movaps		%%xmm0	,%%xmm7		\n\t"\
		"movslq	%[__n_minus_sil],%%rcx	\n\t"\
		"movd	%%rcx	,%%xmm2			\n\t"\
		"pshufd	$0,	%%xmm2	,%%xmm2		\n\t"\
		"psubd		%%xmm0	,%%xmm2		\n\t"\
		"movmskps	%%xmm2	,%%rcx		\n\t"\
		"shlq	$16		,%%rcx			\n\t"\
		"addq	%%rcx	,%%rsi			\n\t"\
		"movq	%[__sinwt]	,%%rdx		\n\t"\
		"movd	%%rdx	,%%xmm3			\n\t"\
		"pshufd	$0,	%%xmm3	,%%xmm3		\n\t"\
		"psubd		%%xmm3	,%%xmm7		\n\t"\
		"movmskps	%%xmm7	,%%rdx		\n\t"\
		"shlq	$8		,%%rdx			\n\t"\
		"addq	%%rdx	,%%rsi			\n\t"\
		"\n\t"\
		"movq	%[__half_arr]	,%%rax	\n\t"\
		"movq	%[__wtA]	,%%rbx		\n\t"\
		"movq	%[__wtB]	,%%rcx		\n\t"\
		"movaps		     (%%rbx),%%xmm2	\n\t	movaps		 0x10(%%rbx),%%xmm6	\n\t"\
		"movhpd		     (%%rcx),%%xmm3	\n\t	movhpd		-0x10(%%rcx),%%xmm7	\n\t"\
		"movlpd		 0x08(%%rcx),%%xmm3	\n\t	movlpd		-0x08(%%rcx),%%xmm7	\n\t"\
		"\n\t"\
		"subq	$0x20	,%%rcx			\n\t"\
		"movq	%%rcx	,%[__wtB]		\n\t"\
		"\n\t"\
		"movq	%%rsi,	%%rdi	\n\t	movq	%%rsi,	%%rbx	\n\t"\
		"shrq	$12,	%%rdi	\n\t	shrq	$14,	%%rbx	\n\t"\
		"andq	$0x0000000000000030	,%%rdi		\n\t	andq	$0x0000000000000030	,%%rbx		\n\t"\
		"movq	%%rsi,	%%rdx	\n\t	movq	%%rsi,	%%rcx	\n\t"\
		"shrq	$4 ,	%%rdx	\n\t	shrq	$6 ,	%%rcx	\n\t"\
		"andq	$0x0000000000000030	,%%rdx		\n\t	andq	$0x0000000000000030	,%%rcx		\n\t"\
		"addq	%%rax	,%%rdi			\n\t	addq	%%rax	,%%rbx			\n\t"\
		"addq	%%rax	,%%rdx			\n\t	addq	%%rax	,%%rcx			\n\t"\
		"\n\t"\
		"mulpd	 	0x100(%%rax),%%xmm2	\n\t	mulpd	 	0x100(%%rax),%%xmm6	\n\t"\
		"mulpd	 	0x110(%%rax),%%xmm3	\n\t	mulpd	 	0x110(%%rax),%%xmm7	\n\t"\
		"mulpd	 	     (%%rdi),%%xmm2	\n\t	mulpd	 	     (%%rbx),%%xmm6	\n\t"\
		"mulpd	 	0x040(%%rdx),%%xmm3	\n\t	mulpd	 	0x040(%%rcx),%%xmm7	\n\t"\
		"\n\t"\
		"movq	%[__cyA]	,%%rcx		\n\t	movq	%[__cyB]	,%%rdx		\n\t"\
		"mulpd		%%xmm3	,%%xmm1		\n\t	mulpd		%%xmm7	,%%xmm5		\n\t"\
		"addpd		(%%rcx)	,%%xmm1		\n\t	addpd		(%%rdx)	,%%xmm5		\n\t"\
		"movaps		%%xmm1	,%%xmm3		\n\t	movaps		%%xmm5	,%%xmm7		\n\t"\
		"roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd		$0,%%xmm7,%%xmm7	\n\t"\
		"\n\t"\
		"movq	%[__sign_mask],%%rbx	\n\t"\
		"subpd			%%xmm3	,%%xmm1	\n\t	subpd			%%xmm7	,%%xmm5	\n\t"\
		"andpd			(%%rbx)	,%%xmm1	\n\t	andpd			(%%rbx)	,%%xmm5	\n\t"\
		"maxpd			%%xmm5	,%%xmm1	\n\t"\
		"maxpd		-0x20(%%rax),%%xmm1	\n\t"\
		"movaps		%%xmm1,-0x20(%%rax)	\n\t"\
		"\n\t"\
		"movq	%%rsi,	%%rdi	\n\t	movq	%%rsi,	%%rbx	\n\t"\
		"shrq	$20,	%%rdi	\n\t	shrq	$22,	%%rbx	\n\t"\
		"andq	$0x0000000000000030	,%%rdi		\n\t	andq	$0x0000000000000030	,%%rbx		\n\t"\
		"addq	%%rax	,%%rdi			\n\t	addq	%%rax	,%%rbx			\n\t"\
		"movaps		%%xmm3	,%%xmm1		\n\t	movaps		%%xmm7	,%%xmm5		\n\t"\
		"mulpd	 	0xc0(%%rdi)	,%%xmm3	\n\t	mulpd	 	0xc0(%%rbx)	,%%xmm7	\n\t"\
		"roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd		$0,%%xmm7,%%xmm7	\n\t"\
		"movaps		%%xmm3	,(%%rcx)	\n\t	movaps		%%xmm7	,(%%rdx)	\n\t"\
		"\n\t"\
		"\n\t"\
		"movq	%[__data]	,%%rax		\n\t"\
		"mulpd	 0x80(%%rdi)	,%%xmm3	\n\t	mulpd	 0x80(%%rbx)	,%%xmm7	\n\t"\
		"subpd		%%xmm3	,	%%xmm1	\n\t	subpd		%%xmm7	,	%%xmm5	\n\t"\
		"mulpd		%%xmm2	,	%%xmm1	\n\t	mulpd		%%xmm6	,	%%xmm5	\n\t"\
		"movaps		%%xmm1	,    (%%rax)\n\t	movaps		%%xmm5	,0x40(%%rax)\n\t"\
		"\n\t"\
		"\n\t"\
		"movq	%[__sse_bw]	,%%rax		\n\t"\
		"movq	%[__sse_nm1]	,%%rbx	\n\t"\
		"paddd		(%%rax)	,%%xmm0		\n\t"\
		"pand		(%%rbx)	,%%xmm0		\n\t"\
	"/**********************************************/\n\t"\
	"/*          Imaginary parts               */\n\t"\
	"/**********************************************/\n\t"\
		"movq	%[__sse_sw]	,%%rdx		\n\t"\
		"movaps	(%%rdx)	,%%xmm1			\n\t"\
		"psubd	%%xmm0	,%%xmm1			\n\t"\
		"movmskps	%%xmm1	,%%rsi		\n\t"\
		"shlq	$24	,%%rsi				\n\t"\
		"movaps	%%xmm0	,%%xmm1			\n\t"\
		"movslq	%[__n_minus_silp1],%%rcx\n\t"\
		"movd	%%rcx	,%%xmm2			\n\t"\
		"pshufd	$0	,%%xmm2	,%%xmm2		\n\t"\
		"psubd	%%xmm0	,%%xmm2			\n\t"\
		"movmskps	%%xmm2	,%%rcx		\n\t"\
		"shlq	$16	,%%rcx				\n\t"\
		"addq	%%rcx	,%%rsi			\n\t"\
		"movq	%[__sinwtm1]	,%%rdx	\n\t"\
		"movd	%%rdx	,%%xmm3			\n\t"\
		"pshufd	$0	,%%xmm3	,%%xmm3		\n\t"\
		"psubd	%%xmm3	,%%xmm1			\n\t"\
		"movmskps	%%xmm1	,%%rdx		\n\t"\
		"shlq	$8	,%%rdx				\n\t"\
		"addq	%%rdx	,%%rsi			\n\t"\
		"movq	%[__data]	,%%rax		\n\t"\
		"movaps	 0x10(%%rax)	,%%xmm1	\n\t"\
		"movaps	 0x50(%%rax)	,%%xmm5	\n\t"\
		"\n\t"\
		"movq	%[__half_arr]	,%%rax	\n\t"\
		"movq	%[__wtA]	,%%rbx		\n\t"\
		"movq	%[__wtC]	,%%rcx		\n\t"\
		"\n\t"\
		"movaps	     (%%rbx)	,%%xmm2	\n\t	movaps	 0x10(%%rbx)	,%%xmm6	\n\t"\
		"movhpd	     (%%rcx)	,%%xmm3	\n\t	movhpd	-0x10(%%rcx)	,%%xmm7	\n\t"\
		"movlpd	 0x08(%%rcx)	,%%xmm3	\n\t	movlpd	-0x08(%%rcx)	,%%xmm7	\n\t"\
		"\n\t"\
		"addq	$0x20	,%%rbx			\n\t"\
		"subq	$0x20	,%%rcx			\n\t"\
		"movq	%%rbx	,%[__wtA]		\n\t"\
		"movq	%%rcx	,%[__wtC]		\n\t"\
		"movq	%%rsi,	%%rdi	\n\t	movq	%%rsi,	%%rbx	\n\t"\
		"shrq	$12,	%%rdi	\n\t	shrq	$14,	%%rbx	\n\t"\
		"andq	$0x0000000000000030	,%%rdi		\n\t	andq	$0x0000000000000030	,%%rbx		\n\t"\
		"movq	%%rsi,	%%rdx	\n\t	movq	%%rsi,	%%rcx	\n\t"\
		"shrq	$4 ,	%%rdx	\n\t	shrq	$6 ,	%%rcx	\n\t"\
		"andq	$0x0000000000000030	,%%rdx		\n\t	andq	$0x0000000000000030	,%%rcx		\n\t"\
		"addq	%%rax	,%%rdi			\n\t	addq	%%rax	,%%rbx			\n\t"\
		"addq	%%rax	,%%rdx			\n\t	addq	%%rax	,%%rcx			\n\t"\
		"\n\t"\
		"mulpd	 0x120(%%rax)	,%%xmm2	\n\t	mulpd	 0x120(%%rax)	,%%xmm6	\n\t"\
		"mulpd	 0x130(%%rax)	,%%xmm3	\n\t	mulpd	 0x130(%%rax)	,%%xmm7	\n\t"\
		"mulpd	      (%%rdi)	,%%xmm2	\n\t	mulpd	      (%%rbx)	,%%xmm6	\n\t"\
		"mulpd	 0x040(%%rdx)	,%%xmm3	\n\t	mulpd	 0x040(%%rcx)	,%%xmm7	\n\t"\
		"\n\t"\
		"movq	%[__cyA]	,%%rcx		\n\t	movq	%[__cyB]	,%%rdx		\n\t"\
		"mulpd	%%xmm3	,%%xmm1			\n\t	mulpd	%%xmm7	,%%xmm5			\n\t"\
		"addpd	(%%rcx)	,%%xmm1			\n\t	addpd	(%%rdx)	,%%xmm5			\n\t"\
		"movaps	%%xmm1	,%%xmm3			\n\t	movaps	%%xmm5	,%%xmm7			\n\t"\
		"roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd		$0,%%xmm7,%%xmm7	\n\t"\
		"\n\t"\
		"movq	%[__sign_mask],%%rbx	\n\t"\
		"subpd	%%xmm3	,%%xmm1			\n\t	subpd	%%xmm7	,%%xmm5			\n\t"\
		"andpd	(%%rbx)	,%%xmm1			\n\t	andpd	(%%rbx)	,%%xmm5			\n\t"\
		"maxpd	%%xmm5	,%%xmm1			\n\t"\
		"maxpd	-0x20(%%rax)	,%%xmm1	\n\t"\
		"movaps	%%xmm1	,-0x20(%%rax)	\n\t"\
		"\n\t"\
		"movq	%%rsi,	%%rdi	\n\t	movq	%%rsi,	%%rbx	\n\t"\
		"shrq	$20,	%%rdi	\n\t	shrq	$22,	%%rbx	\n\t"\
		"andq	$0x0000000000000030	,%%rdi		\n\t	andq	$0x0000000000000030	,%%rbx		\n\t"\
		"addq	%%rax	,%%rdi			\n\t	addq	%%rax	,%%rbx			\n\t"\
		"movaps	%%xmm3	,%%xmm1			\n\t	movaps	%%xmm7	,%%xmm5			\n\t"\
		"mulpd	 0xc0(%%rdi)	,%%xmm3	\n\t	mulpd	 0xc0(%%rbx)	,%%xmm7	\n\t"\
		"roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd		$0,%%xmm7,%%xmm7	\n\t"\
		"movaps	%%xmm3	,(%%rcx)		\n\t	movaps	%%xmm7	,(%%rdx)		\n\t"\
		"\n\t"\
		"\n\t"\
		"movq	%[__data]	,%%rax		\n\t"\
		"mulpd	 0x80(%%rdi)	,%%xmm3	\n\t	mulpd	 0x80(%%rbx)	,%%xmm7	\n\t"\
		"subpd	%%xmm3	,%%xmm1			\n\t	subpd	%%xmm7	,%%xmm5			\n\t"\
		"mulpd	%%xmm2	,%%xmm1			\n\t	mulpd	%%xmm6	,%%xmm5			\n\t"\
		"movaps	%%xmm1	, 0x10(%%rax)	\n\t	movaps	%%xmm5	, 0x50(%%rax)	\n\t"\
		"\n\t"\
		"\n\t"\
		"movq	%[__sse_bw]	,%%rax		\n\t"\
		"movq	%[__sse_nm1]	,%%rbx	\n\t"\
		"paddd	(%%rax)	,%%xmm0			\n\t"\
		"pand	(%%rbx)	,%%xmm0			\n\t"\
		"movq	%[__bjmod_0],%%rcx		\n\t"\
		"movaps	%%xmm0,(%%rcx)			\n\t"\
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
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","xmm0","xmm1","xmm2","xmm3","xmm5","xmm6","xmm7"		/* Clobbered registers */\
	);\
	}

	/********* Packed 32-bit-int version of SSE2_cmplx_carry_norm_pow2_errcheck2_2x:***********/
	#define SSE2_cmplx_carry_norm_pow2_errcheck2_2B(Xdata,XwtA,XwtB,XcyA,XcyB,Xbjmod_0,Xhalf_arr,Xn_minus_silp1,Xn_minus_sil,Xsign_mask,Xsinwt,Xsinwtm1,Xsse_bw,Xsse_nm1,Xsse_sw)\
	{\
	__asm__ volatile (\
	"/**********************************************/\n\t"\
	"/*          Real      parts                   */\n\t"\
	"/**********************************************/\n\t"\
		"movq	%[__bjmod_0],	%%rax	\n\t"\
		"movaps		(%%rax)	,	%%xmm0	\n\t"\
		"movq	%[__sse_sw]	,	%%rbx	\n\t"\
		"movaps		(%%rbx)	,	%%xmm1	\n\t"\
		"psubd		%%xmm0	,	%%xmm1	\n\t"\
		"movmskps	%%xmm1	,	%%rsi	\n\t"\
		"\n\t"\
		"shlq	$24		,%%rsi			\n\t"\
		"movaps		%%xmm0	,%%xmm1		\n\t"\
		"movslq	%[__n_minus_sil],%%rcx	\n\t"\
		"movd	%%rcx	,%%xmm2			\n\t"\
		"pshufd	$0,	%%xmm2	,%%xmm2		\n\t"\
		"psubd		%%xmm0	,%%xmm2		\n\t"\
		"movmskps	%%xmm2	,%%rcx		\n\t"\
		"shlq	$16		,%%rcx			\n\t"\
		"addq	%%rcx	,%%rsi			\n\t"\
		"movq	%[__sinwt]	,%%rdx		\n\t"\
		"movd	%%rdx	,%%xmm3			\n\t"\
		"pshufd	$0,	%%xmm3	,%%xmm3		\n\t"\
		"psubd		%%xmm3	,%%xmm1		\n\t"\
		"movmskps	%%xmm1	,%%rdx		\n\t"\
		"shlq	$8		,%%rdx			\n\t"\
		"addq	%%rdx	,%%rsi			\n\t"\
		"\n\t"\
		"mov	%[__data],%%rax			\n\t"\
		"movaps	 0x20(%%rax),%%xmm1		\n\t	movaps		 0x60(%%rax),%%xmm5\n\t"\
		"\n\t"\
		"movq	%[__half_arr]	,%%rax	\n\t"\
		"movq	%[__wtA]	,%%rbx		\n\t"\
		"movq	%[__wtB]	,%%rcx		\n\t"\
		"movaps		     (%%rbx),%%xmm2	\n\t	movaps		 0x10(%%rbx),%%xmm6	\n\t"\
		"movhpd		     (%%rcx),%%xmm3	\n\t	movhpd		-0x10(%%rcx),%%xmm7	\n\t"\
		"movlpd		 0x08(%%rcx),%%xmm3	\n\t	movlpd		-0x08(%%rcx),%%xmm7	\n\t"\
		"\n\t"\
		"movq	%%rsi,	%%rdi	\n\t	movq	%%rsi,	%%rbx	\n\t"\
		"shrq	$12,	%%rdi	\n\t	shrq	$14,	%%rbx	\n\t"\
		"andq	$0x0000000000000030	,%%rdi		\n\t	andq	$0x0000000000000030	,%%rbx		\n\t"\
		"movq	%%rsi,	%%rdx	\n\t	movq	%%rsi,	%%rcx	\n\t"\
		"shrq	$4 ,	%%rdx	\n\t	shrq	$6 ,	%%rcx	\n\t"\
		"andq	$0x0000000000000030	,%%rdx		\n\t	andq	$0x0000000000000030	,%%rcx		\n\t"\
		"addq	%%rax	,%%rdi			\n\t	addq	%%rax	,%%rbx			\n\t"\
		"addq	%%rax	,%%rdx			\n\t	addq	%%rax	,%%rcx			\n\t"\
		"\n\t"\
		"mulpd	 	0x100(%%rax),%%xmm2	\n\t	mulpd	 	0x100(%%rax),%%xmm6	\n\t"\
		"mulpd	 	0x110(%%rax),%%xmm3	\n\t	mulpd	 	0x110(%%rax),%%xmm7	\n\t"\
		"mulpd	 	     (%%rdi),%%xmm2	\n\t	mulpd	 	     (%%rbx),%%xmm6	\n\t"\
		"mulpd	 	0x040(%%rdx),%%xmm3	\n\t	mulpd	 	0x040(%%rcx),%%xmm7	\n\t"\
		"\n\t"\
		"movq	%[__cyA]	,%%rcx		\n\t	movq	%[__cyB]	,%%rdx		\n\t"\
		"mulpd		%%xmm3	,%%xmm1		\n\t	mulpd		%%xmm7	,%%xmm5		\n\t"\
		"addpd		(%%rcx)	,%%xmm1		\n\t	addpd		(%%rdx)	,%%xmm5		\n\t"\
		"movaps		%%xmm1	,%%xmm3		\n\t	movaps		%%xmm5	,%%xmm7		\n\t"\
		"roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd		$0,%%xmm7,%%xmm7	\n\t"\
		"\n\t"\
		"movq	%[__sign_mask],%%rbx	\n\t"\
		"subpd			%%xmm3	,%%xmm1	\n\t	subpd			%%xmm7	,%%xmm5	\n\t"\
		"andpd			(%%rbx)	,%%xmm1	\n\t	andpd			(%%rbx)	,%%xmm5	\n\t"\
		"maxpd			%%xmm5	,%%xmm1	\n\t"\
		"maxpd		-0x20(%%rax),%%xmm1	\n\t"\
		"movaps		%%xmm1,-0x20(%%rax)	\n\t"\
		"\n\t"\
		"movq	%%rsi,	%%rdi	\n\t	movq	%%rsi,	%%rbx	\n\t"\
		"shrq	$20,	%%rdi	\n\t	shrq	$22,	%%rbx	\n\t"\
		"andq	$0x0000000000000030	,%%rdi		\n\t	andq	$0x0000000000000030	,%%rbx		\n\t"\
		"addq	%%rax	,%%rdi			\n\t	addq	%%rax	,%%rbx			\n\t"\
		"movaps		%%xmm3	,%%xmm1		\n\t	movaps		%%xmm7	,%%xmm5		\n\t"\
		"mulpd	 	0xc0(%%rdi)	,%%xmm3	\n\t	mulpd	 	0xc0(%%rbx)	,%%xmm7	\n\t"\
		"roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd		$0,%%xmm7,%%xmm7	\n\t"\
		"movaps		%%xmm3	,(%%rcx)	\n\t	movaps		%%xmm7	,(%%rdx)	\n\t"\
		"\n\t"\
		"\n\t"\
		"movq	%[__data]	,%%rax		\n\t"\
		"mulpd	 0x80(%%rdi)	,%%xmm3	\n\t	mulpd	 0x80(%%rbx)	,%%xmm7	\n\t"\
		"subpd		%%xmm3	,	%%xmm1	\n\t	subpd		%%xmm7	,	%%xmm5	\n\t"\
		"mulpd		%%xmm2	,	%%xmm1	\n\t	mulpd		%%xmm6	,	%%xmm5	\n\t"\
		"movaps		%%xmm1	,0x20(%%rax)\n\t	movaps		%%xmm5	,0x60(%%rax)\n\t"\
		"\n\t"\
		"\n\t"\
		"movq	%[__sse_bw]	,%%rax		\n\t"\
		"movq	%[__sse_nm1]	,%%rbx	\n\t"\
		"paddd		(%%rax)	,%%xmm0		\n\t"\
		"pand		(%%rbx)	,%%xmm0		\n\t"\
	"/**********************************************/\n\t"\
	"/*          Imaginary parts               */\n\t"\
	"/**********************************************/\n\t"\
		"movq	%[__sse_sw]	,%%rdx		\n\t"\
		"movaps	(%%rdx)	,%%xmm1			\n\t"\
		"psubd	%%xmm0	,%%xmm1			\n\t"\
		"movmskps	%%xmm1	,%%rsi		\n\t"\
		"shlq	$24	,%%rsi				\n\t"\
		"movaps	%%xmm0	,%%xmm1			\n\t"\
		"movslq	%[__n_minus_silp1],%%rcx\n\t"\
		"movd	%%rcx	,%%xmm2			\n\t"\
		"pshufd	$0	,%%xmm2	,%%xmm2		\n\t"\
		"psubd	%%xmm0	,%%xmm2			\n\t"\
		"movmskps	%%xmm2	,%%rcx		\n\t"\
		"shlq	$16	,%%rcx				\n\t"\
		"addq	%%rcx	,%%rsi			\n\t"\
		"movq	%[__sinwtm1]	,%%rdx	\n\t"\
		"movd	%%rdx	,%%xmm3			\n\t"\
		"pshufd	$0	,%%xmm3	,%%xmm3		\n\t"\
		"psubd	%%xmm3	,%%xmm1			\n\t"\
		"movmskps	%%xmm1	,%%rdx		\n\t"\
		"shlq	$8	,%%rdx				\n\t"\
		"addq	%%rdx	,%%rsi			\n\t"\
		"movq	%[__data]	,%%rax		\n\t"\
		"movaps	 0x30(%%rax)	,%%xmm1	\n\t"\
		"movaps	 0x70(%%rax)	,%%xmm5	\n\t"\
		"\n\t"\
		"movq	%[__half_arr]	,%%rax	\n\t"\
		"movq	%[__wtA]	,%%rbx		\n\t"\
		"movq	%[__wtB]	,%%rcx		\n\t"\
		"\n\t"\
		"movaps	     (%%rbx)	,%%xmm2	\n\t	movaps	 0x10(%%rbx)	,%%xmm6	\n\t"\
		"movhpd	     (%%rcx)	,%%xmm3	\n\t	movhpd	-0x10(%%rcx)	,%%xmm7	\n\t"\
		"movlpd	 0x08(%%rcx)	,%%xmm3	\n\t	movlpd	-0x08(%%rcx)	,%%xmm7	\n\t"\
		"\n\t"\
		"addq	$0x20	,%%rbx			\n\t"\
		"subq	$0x20	,%%rcx			\n\t"\
		"movq	%%rbx	,%[__wtA]		\n\t"\
		"movq	%%rcx	,%[__wtB]		\n\t"\
		"movq	%%rsi,	%%rdi	\n\t	movq	%%rsi,	%%rbx	\n\t"\
		"shrq	$12,	%%rdi	\n\t	shrq	$14,	%%rbx	\n\t"\
		"andq	$0x0000000000000030	,%%rdi		\n\t	andq	$0x0000000000000030	,%%rbx		\n\t"\
		"movq	%%rsi,	%%rdx	\n\t	movq	%%rsi,	%%rcx	\n\t"\
		"shrq	$4 ,	%%rdx	\n\t	shrq	$6 ,	%%rcx	\n\t"\
		"andq	$0x0000000000000030	,%%rdx		\n\t	andq	$0x0000000000000030	,%%rcx		\n\t"\
		"addq	%%rax	,%%rdi			\n\t	addq	%%rax	,%%rbx			\n\t"\
		"addq	%%rax	,%%rdx			\n\t	addq	%%rax	,%%rcx			\n\t"\
		"\n\t"\
		"mulpd	 0x120(%%rax)	,%%xmm2	\n\t	mulpd	 0x120(%%rax)	,%%xmm6	\n\t"\
		"mulpd	 0x130(%%rax)	,%%xmm3	\n\t	mulpd	 0x130(%%rax)	,%%xmm7	\n\t"\
		"mulpd	      (%%rdi)	,%%xmm2	\n\t	mulpd	      (%%rbx)	,%%xmm6	\n\t"\
		"mulpd	 0x040(%%rdx)	,%%xmm3	\n\t	mulpd	 0x040(%%rcx)	,%%xmm7	\n\t"\
		"\n\t"\
		"movq	%[__cyA]	,%%rcx		\n\t	movq	%[__cyB]	,%%rdx		\n\t"\
		"mulpd	%%xmm3	,%%xmm1			\n\t	mulpd	%%xmm7	,%%xmm5			\n\t"\
		"addpd	(%%rcx)	,%%xmm1			\n\t	addpd	(%%rdx)	,%%xmm5			\n\t"\
		"movaps	%%xmm1	,%%xmm3			\n\t	movaps	%%xmm5	,%%xmm7			\n\t"\
		"roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd		$0,%%xmm7,%%xmm7	\n\t"\
		"\n\t"\
		"movq	%[__sign_mask],%%rbx	\n\t"\
		"subpd	%%xmm3	,%%xmm1			\n\t	subpd	%%xmm7	,%%xmm5			\n\t"\
		"andpd	(%%rbx)	,%%xmm1			\n\t	andpd	(%%rbx)	,%%xmm5			\n\t"\
		"maxpd	%%xmm5	,%%xmm1			\n\t"\
		"maxpd	-0x20(%%rax)	,%%xmm1	\n\t"\
		"movaps	%%xmm1	,-0x20(%%rax)	\n\t"\
		"\n\t"\
		"movq	%%rsi,	%%rdi	\n\t	movq	%%rsi,	%%rbx	\n\t"\
		"shrq	$20,	%%rdi	\n\t	shrq	$22,	%%rbx	\n\t"\
		"andq	$0x0000000000000030	,%%rdi		\n\t	andq	$0x0000000000000030	,%%rbx		\n\t"\
		"addq	%%rax	,%%rdi			\n\t	addq	%%rax	,%%rbx			\n\t"\
		"movaps	%%xmm3	,%%xmm1			\n\t	movaps	%%xmm7	,%%xmm5			\n\t"\
		"mulpd	 0xc0(%%rdi)	,%%xmm3	\n\t	mulpd	 0xc0(%%rbx)	,%%xmm7	\n\t"\
		"roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd		$0,%%xmm7,%%xmm7	\n\t"\
		"movaps	%%xmm3	,(%%rcx)		\n\t	movaps	%%xmm7	,(%%rdx)		\n\t"\
		"\n\t"\
		"\n\t"\
		"movq	%[__data]	,%%rax		\n\t"\
		"mulpd	 0x80(%%rdi)	,%%xmm3	\n\t	mulpd	 0x80(%%rbx)	,%%xmm7	\n\t"\
		"subpd	%%xmm3	,%%xmm1			\n\t	subpd	%%xmm7	,%%xmm5			\n\t"\
		"mulpd	%%xmm2	,%%xmm1			\n\t	mulpd	%%xmm6	,%%xmm5			\n\t"\
		"movaps	%%xmm1	, 0x30(%%rax)	\n\t	movaps	%%xmm5	, 0x70(%%rax)	\n\t"\
		"\n\t"\
		"\n\t"\
		"movq	%[__sse_bw]	,%%rax		\n\t"\
		"movq	%[__sse_nm1]	,%%rbx	\n\t"\
		"paddd	(%%rax)	,%%xmm0			\n\t"\
		"pand	(%%rbx)	,%%xmm0			\n\t"\
		"movq	%[__bjmod_0],%%rcx		\n\t"\
		"movaps	%%xmm0,(%%rcx)			\n\t"\
	"/**********************************************/\n\t"\
	"/*          Repack the data:              */\n\t"\
	"/**********************************************/\n\t"\
		"mov	%[__data],%%rax			\n\t"\
		"movaps		0x10(%%rax)	,%%xmm1	\n\t	movaps		0x50(%%rax)	,%%xmm5	\n\t"\
		"movaps		    (%%rax)	,%%xmm0	\n\t	movaps		0x40(%%rax)	,%%xmm4	\n\t"\
		"movaps		%%xmm1		,%%xmm3	\n\t	movaps		%%xmm5		,%%xmm7	\n\t"\
		"movaps		%%xmm0		,%%xmm2	\n\t	movaps		%%xmm4		,%%xmm6	\n\t"\
		"unpckhpd	0x30(%%rax)	,%%xmm3	\n\t	unpckhpd	0x70(%%rax)	,%%xmm7	\n\t"\
		"unpcklpd	0x30(%%rax)	,%%xmm1	\n\t	unpcklpd	0x70(%%rax)	,%%xmm5	\n\t"\
		"movaps		%%xmm3,0x30(%%rax)	\n\t	movaps		%%xmm7,0x70(%%rax)	\n\t"\
		"unpckhpd	0x20(%%rax)	,%%xmm2	\n\t	unpckhpd	0x60(%%rax)	,%%xmm6	\n\t"\
		"unpcklpd	0x20(%%rax)	,%%xmm0	\n\t	unpcklpd	0x60(%%rax)	,%%xmm4	\n\t"\
		"movaps		%%xmm2,0x20(%%rax)	\n\t	movaps		%%xmm6,0x60(%%rax)	\n\t"\
		"movaps		%%xmm1,0x10(%%rax)	\n\t	movaps		%%xmm5,0x50(%%rax)	\n\t"\
		"movaps		%%xmm0,    (%%rax)	\n\t	movaps		%%xmm4,0x40(%%rax)	\n\t"\
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
		, [__sse_nm1]		"m" (Xsse_nm1)		\
		, [__sse_sw]		"m" (Xsse_sw)		\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
	);\
	}

	/******************************************************************************************************************************************************************/
	/********** No-ROE-Check versions of the latter 2 of the above 3 macros - the first is called too infrequently to bother with a special non-ROE version: **********/
	/******************************************************************************************************************************************************************/

	#define SSE2_cmplx_carry_norm_pow2_nocheck1_2B(Xdata,XwtA,XwtB,XwtC,XcyA,XcyB,Xbjmod_0,Xhalf_arr,Xn_minus_silp1,Xn_minus_sil,Xsinwt,Xsinwtm1,Xsse_bw,Xsse_nm1,Xsse_sw)\
	{\
	__asm__ volatile (\
	"/***************Unpack the data:*************************/\n\t"\
		"movq	%[__data]	,%%rax	\n\t"\
		"movaps		    (%%rax)	,%%xmm1	\n\t	movaps		0x40(%%rax)	,%%xmm5	\n\t"\
		"movaps		    (%%rax)	,%%xmm2	\n\t	movaps		0x40(%%rax)	,%%xmm6	\n\t"\
		"unpcklpd	0x20(%%rax)	,%%xmm1	\n\t	unpcklpd	0x60(%%rax)	,%%xmm5	\n\t"\
		"unpckhpd	0x20(%%rax)	,%%xmm2	\n\t	unpckhpd	0x60(%%rax)	,%%xmm6	\n\t"\
		"movaps		%%xmm2, 0x20(%%rax)	\n\t	movaps		%%xmm6, 0x60(%%rax)	\n\t"\
		"\n\t"\
		"movaps		0x10(%%rax)	,%%xmm2	\n\t	movaps		0x50(%%rax)	,%%xmm6	\n\t"\
		"movaps		0x10(%%rax)	,%%xmm3	\n\t	movaps		0x50(%%rax)	,%%xmm7	\n\t"\
		"unpcklpd	0x30(%%rax)	,%%xmm2	\n\t	unpcklpd	0x70(%%rax)	,%%xmm6	\n\t"\
		"unpckhpd	0x30(%%rax)	,%%xmm3	\n\t	unpckhpd	0x70(%%rax)	,%%xmm7	\n\t"\
		"movaps		%%xmm2, 0x10(%%rax)	\n\t	movaps		%%xmm6, 0x50(%%rax)	\n\t"\
		"movaps		%%xmm3, 0x30(%%rax)	\n\t	movaps		%%xmm7, 0x70(%%rax)	\n\t"\
	"/**********************************************/\n\t"\
	"/*          Real      parts                   */\n\t"\
	"/**********************************************/\n\t"\
		"movq	%[__bjmod_0],	%%rax	\n\t"\
		"movaps		(%%rax)	,	%%xmm0	\n\t"\
		"movq	%[__sse_sw]	,	%%rbx	\n\t"\
		"movaps		(%%rbx)	,	%%xmm7	\n\t"\
		"psubd		%%xmm0	,	%%xmm7	\n\t"\
		"movmskps	%%xmm7	,	%%rsi	\n\t"\
		"\n\t"\
		"shlq	$24		,%%rsi			\n\t"\
		"movaps		%%xmm0	,%%xmm7		\n\t"\
		"movslq	%[__n_minus_sil],%%rcx	\n\t"\
		"movd	%%rcx	,%%xmm2			\n\t"\
		"pshufd	$0,	%%xmm2	,%%xmm2		\n\t"\
		"psubd		%%xmm0	,%%xmm2		\n\t"\
		"movmskps	%%xmm2	,%%rcx		\n\t"\
		"shlq	$16		,%%rcx			\n\t"\
		"addq	%%rcx	,%%rsi			\n\t"\
		"movq	%[__sinwt]	,%%rdx		\n\t"\
		"movd	%%rdx	,%%xmm3			\n\t"\
		"pshufd	$0,	%%xmm3	,%%xmm3		\n\t"\
		"psubd		%%xmm3	,%%xmm7		\n\t"\
		"movmskps	%%xmm7	,%%rdx		\n\t"\
		"shlq	$8		,%%rdx			\n\t"\
		"addq	%%rdx	,%%rsi			\n\t"\
		"\n\t"\
		"movq	%[__half_arr]	,%%rax	\n\t"\
		"movq	%[__wtA]	,%%rbx		\n\t"\
		"movq	%[__wtB]	,%%rcx		\n\t"\
		"movaps		     (%%rbx),%%xmm2	\n\t	movaps		 0x10(%%rbx),%%xmm6	\n\t"\
		"movhpd		     (%%rcx),%%xmm3	\n\t	movhpd		-0x10(%%rcx),%%xmm7	\n\t"\
		"movlpd		 0x08(%%rcx),%%xmm3	\n\t	movlpd		-0x08(%%rcx),%%xmm7	\n\t"\
		"\n\t"\
		"subq	$0x20	,%%rcx			\n\t"\
		"movq	%%rcx	,%[__wtB]		\n\t"\
		"\n\t"\
		"movq	%%rsi,	%%rdi	\n\t	movq	%%rsi,	%%rbx	\n\t"\
		"shrq	$12,	%%rdi	\n\t	shrq	$14,	%%rbx	\n\t"\
		"andq	$0x0000000000000030	,%%rdi		\n\t	andq	$0x0000000000000030	,%%rbx		\n\t"\
		"movq	%%rsi,	%%rdx	\n\t	movq	%%rsi,	%%rcx	\n\t"\
		"shrq	$4 ,	%%rdx	\n\t	shrq	$6 ,	%%rcx	\n\t"\
		"andq	$0x0000000000000030	,%%rdx		\n\t	andq	$0x0000000000000030	,%%rcx		\n\t"\
		"addq	%%rax	,%%rdi			\n\t	addq	%%rax	,%%rbx			\n\t"\
		"addq	%%rax	,%%rdx			\n\t	addq	%%rax	,%%rcx			\n\t"\
		"\n\t"\
		"mulpd	 	0x100(%%rax),%%xmm2	\n\t	mulpd	 	0x100(%%rax),%%xmm6	\n\t"\
		"mulpd	 	0x110(%%rax),%%xmm3	\n\t	mulpd	 	0x110(%%rax),%%xmm7	\n\t"\
		"mulpd	 	     (%%rdi),%%xmm2	\n\t	mulpd	 	     (%%rbx),%%xmm6	\n\t"\
		"mulpd	 	0x040(%%rdx),%%xmm3	\n\t	mulpd	 	0x040(%%rcx),%%xmm7	\n\t"\
		"\n\t"\
		"movq	%[__cyA]	,%%rcx		\n\t	movq	%[__cyB]	,%%rdx		\n\t"\
		"mulpd		%%xmm3	,%%xmm1		\n\t	mulpd		%%xmm7	,%%xmm5		\n\t"\
		"addpd		(%%rcx)	,%%xmm1		\n\t	addpd		(%%rdx)	,%%xmm5		\n\t"\
		"movaps		%%xmm1	,%%xmm3		\n\t	movaps		%%xmm5	,%%xmm7		\n\t"\
		"roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd		$0,%%xmm7,%%xmm7	\n\t"\
		"\n\t"\
		"/* NO ROE HERE */		\n\t"\
		"\n\t"\
		"movq	%%rsi,	%%rdi	\n\t	movq	%%rsi,	%%rbx	\n\t"\
		"shrq	$20,	%%rdi	\n\t	shrq	$22,	%%rbx	\n\t"\
		"andq	$0x0000000000000030	,%%rdi		\n\t	andq	$0x0000000000000030	,%%rbx		\n\t"\
		"addq	%%rax	,%%rdi			\n\t	addq	%%rax	,%%rbx			\n\t"\
		"movaps		%%xmm3	,%%xmm1		\n\t	movaps		%%xmm7	,%%xmm5		\n\t"\
		"mulpd	 	0xc0(%%rdi)	,%%xmm3	\n\t	mulpd	 	0xc0(%%rbx)	,%%xmm7	\n\t"\
		"roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd		$0,%%xmm7,%%xmm7	\n\t"\
		"movaps		%%xmm3	,(%%rcx)	\n\t	movaps		%%xmm7	,(%%rdx)	\n\t"\
		"\n\t"\
		"\n\t"\
		"movq	%[__data]	,%%rax		\n\t"\
		"mulpd	 0x80(%%rdi)	,%%xmm3	\n\t	mulpd	 0x80(%%rbx)	,%%xmm7	\n\t"\
		"subpd		%%xmm3	,	%%xmm1	\n\t	subpd		%%xmm7	,	%%xmm5	\n\t"\
		"mulpd		%%xmm2	,	%%xmm1	\n\t	mulpd		%%xmm6	,	%%xmm5	\n\t"\
		"movaps		%%xmm1	,    (%%rax)\n\t	movaps		%%xmm5	,0x40(%%rax)\n\t"\
		"\n\t"\
		"\n\t"\
		"movq	%[__sse_bw]	,%%rax		\n\t"\
		"movq	%[__sse_nm1]	,%%rbx	\n\t"\
		"paddd		(%%rax)	,%%xmm0		\n\t"\
		"pand		(%%rbx)	,%%xmm0		\n\t"\
	"/**********************************************/\n\t"\
	"/*          Imaginary parts               */\n\t"\
	"/**********************************************/\n\t"\
		"movq	%[__sse_sw]	,%%rdx		\n\t"\
		"movaps	(%%rdx)	,%%xmm1			\n\t"\
		"psubd	%%xmm0	,%%xmm1			\n\t"\
		"movmskps	%%xmm1	,%%rsi		\n\t"\
		"shlq	$24	,%%rsi				\n\t"\
		"movaps	%%xmm0	,%%xmm1			\n\t"\
		"movslq	%[__n_minus_silp1],%%rcx\n\t"\
		"movd	%%rcx	,%%xmm2			\n\t"\
		"pshufd	$0	,%%xmm2	,%%xmm2		\n\t"\
		"psubd	%%xmm0	,%%xmm2			\n\t"\
		"movmskps	%%xmm2	,%%rcx		\n\t"\
		"shlq	$16	,%%rcx				\n\t"\
		"addq	%%rcx	,%%rsi			\n\t"\
		"movq	%[__sinwtm1]	,%%rdx	\n\t"\
		"movd	%%rdx	,%%xmm3			\n\t"\
		"pshufd	$0	,%%xmm3	,%%xmm3		\n\t"\
		"psubd	%%xmm3	,%%xmm1			\n\t"\
		"movmskps	%%xmm1	,%%rdx		\n\t"\
		"shlq	$8	,%%rdx				\n\t"\
		"addq	%%rdx	,%%rsi			\n\t"\
		"movq	%[__data]	,%%rax		\n\t"\
		"movaps	 0x10(%%rax)	,%%xmm1	\n\t"\
		"movaps	 0x50(%%rax)	,%%xmm5	\n\t"\
		"\n\t"\
		"movq	%[__half_arr]	,%%rax	\n\t"\
		"movq	%[__wtA]	,%%rbx		\n\t"\
		"movq	%[__wtC]	,%%rcx		\n\t"\
		"\n\t"\
		"movaps	     (%%rbx)	,%%xmm2	\n\t	movaps	 0x10(%%rbx)	,%%xmm6	\n\t"\
		"movhpd	     (%%rcx)	,%%xmm3	\n\t	movhpd	-0x10(%%rcx)	,%%xmm7	\n\t"\
		"movlpd	 0x08(%%rcx)	,%%xmm3	\n\t	movlpd	-0x08(%%rcx)	,%%xmm7	\n\t"\
		"\n\t"\
		"addq	$0x20	,%%rbx			\n\t"\
		"subq	$0x20	,%%rcx			\n\t"\
		"movq	%%rbx	,%[__wtA]		\n\t"\
		"movq	%%rcx	,%[__wtC]		\n\t"\
		"movq	%%rsi,	%%rdi	\n\t	movq	%%rsi,	%%rbx	\n\t"\
		"shrq	$12,	%%rdi	\n\t	shrq	$14,	%%rbx	\n\t"\
		"andq	$0x0000000000000030	,%%rdi		\n\t	andq	$0x0000000000000030	,%%rbx		\n\t"\
		"movq	%%rsi,	%%rdx	\n\t	movq	%%rsi,	%%rcx	\n\t"\
		"shrq	$4 ,	%%rdx	\n\t	shrq	$6 ,	%%rcx	\n\t"\
		"andq	$0x0000000000000030	,%%rdx		\n\t	andq	$0x0000000000000030	,%%rcx		\n\t"\
		"addq	%%rax	,%%rdi			\n\t	addq	%%rax	,%%rbx			\n\t"\
		"addq	%%rax	,%%rdx			\n\t	addq	%%rax	,%%rcx			\n\t"\
		"\n\t"\
		"mulpd	 0x120(%%rax)	,%%xmm2	\n\t	mulpd	 0x120(%%rax)	,%%xmm6	\n\t"\
		"mulpd	 0x130(%%rax)	,%%xmm3	\n\t	mulpd	 0x130(%%rax)	,%%xmm7	\n\t"\
		"mulpd	      (%%rdi)	,%%xmm2	\n\t	mulpd	      (%%rbx)	,%%xmm6	\n\t"\
		"mulpd	 0x040(%%rdx)	,%%xmm3	\n\t	mulpd	 0x040(%%rcx)	,%%xmm7	\n\t"\
		"\n\t"\
		"movq	%[__cyA]	,%%rcx		\n\t	movq	%[__cyB]	,%%rdx		\n\t"\
		"mulpd	%%xmm3	,%%xmm1			\n\t	mulpd	%%xmm7	,%%xmm5			\n\t"\
		"addpd	(%%rcx)	,%%xmm1			\n\t	addpd	(%%rdx)	,%%xmm5			\n\t"\
		"movaps	%%xmm1	,%%xmm3			\n\t	movaps	%%xmm5	,%%xmm7			\n\t"\
		"roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd		$0,%%xmm7,%%xmm7	\n\t"\
		"\n\t"\
		"/* NO ROE HERE */		\n\t"\
		"\n\t"\
		"movq	%%rsi,	%%rdi	\n\t	movq	%%rsi,	%%rbx	\n\t"\
		"shrq	$20,	%%rdi	\n\t	shrq	$22,	%%rbx	\n\t"\
		"andq	$0x0000000000000030	,%%rdi		\n\t	andq	$0x0000000000000030	,%%rbx		\n\t"\
		"addq	%%rax	,%%rdi			\n\t	addq	%%rax	,%%rbx			\n\t"\
		"movaps	%%xmm3	,%%xmm1			\n\t	movaps	%%xmm7	,%%xmm5			\n\t"\
		"mulpd	 0xc0(%%rdi)	,%%xmm3	\n\t	mulpd	 0xc0(%%rbx)	,%%xmm7	\n\t"\
		"roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd		$0,%%xmm7,%%xmm7	\n\t"\
		"movaps	%%xmm3	,(%%rcx)		\n\t	movaps	%%xmm7	,(%%rdx)		\n\t"\
		"\n\t"\
		"\n\t"\
		"movq	%[__data]	,%%rax		\n\t"\
		"mulpd	 0x80(%%rdi)	,%%xmm3	\n\t	mulpd	 0x80(%%rbx)	,%%xmm7	\n\t"\
		"subpd	%%xmm3	,%%xmm1			\n\t	subpd	%%xmm7	,%%xmm5			\n\t"\
		"mulpd	%%xmm2	,%%xmm1			\n\t	mulpd	%%xmm6	,%%xmm5			\n\t"\
		"movaps	%%xmm1	, 0x10(%%rax)	\n\t	movaps	%%xmm5	, 0x50(%%rax)	\n\t"\
		"\n\t"\
		"\n\t"\
		"movq	%[__sse_bw]	,%%rax		\n\t"\
		"movq	%[__sse_nm1]	,%%rbx	\n\t"\
		"paddd	(%%rax)	,%%xmm0			\n\t"\
		"pand	(%%rbx)	,%%xmm0			\n\t"\
		"movq	%[__bjmod_0],%%rcx		\n\t"\
		"movaps	%%xmm0,(%%rcx)			\n\t"\
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
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","xmm0","xmm1","xmm2","xmm3","xmm5","xmm6","xmm7"		/* Clobbered registers */\
	);\
	}

	#define SSE2_cmplx_carry_norm_pow2_nocheck2_2B(Xdata,XwtA,XwtB,XcyA,XcyB,Xbjmod_0,Xhalf_arr,Xn_minus_silp1,Xn_minus_sil,Xsinwt,Xsinwtm1,Xsse_bw,Xsse_nm1,Xsse_sw)\
	{\
	__asm__ volatile (\
	"/**********************************************/\n\t"\
	"/*          Real      parts                   */\n\t"\
	"/**********************************************/\n\t"\
		"movq	%[__bjmod_0],	%%rax	\n\t"\
		"movaps		(%%rax)	,	%%xmm0	\n\t"\
		"movq	%[__sse_sw]	,	%%rbx	\n\t"\
		"movaps		(%%rbx)	,	%%xmm1	\n\t"\
		"psubd		%%xmm0	,	%%xmm1	\n\t"\
		"movmskps	%%xmm1	,	%%rsi	\n\t"\
		"\n\t"\
		"shlq	$24		,%%rsi			\n\t"\
		"movaps		%%xmm0	,%%xmm1		\n\t"\
		"movslq	%[__n_minus_sil],%%rcx	\n\t"\
		"movd	%%rcx	,%%xmm2			\n\t"\
		"pshufd	$0,	%%xmm2	,%%xmm2		\n\t"\
		"psubd		%%xmm0	,%%xmm2		\n\t"\
		"movmskps	%%xmm2	,%%rcx		\n\t"\
		"shlq	$16		,%%rcx			\n\t"\
		"addq	%%rcx	,%%rsi			\n\t"\
		"movq	%[__sinwt]	,%%rdx		\n\t"\
		"movd	%%rdx	,%%xmm3			\n\t"\
		"pshufd	$0,	%%xmm3	,%%xmm3		\n\t"\
		"psubd		%%xmm3	,%%xmm1		\n\t"\
		"movmskps	%%xmm1	,%%rdx		\n\t"\
		"shlq	$8		,%%rdx			\n\t"\
		"addq	%%rdx	,%%rsi			\n\t"\
		"\n\t"\
		"mov	%[__data],%%rax			\n\t"\
		"movaps	 0x20(%%rax),%%xmm1		\n\t	movaps		 0x60(%%rax),%%xmm5\n\t"\
		"\n\t"\
		"movq	%[__half_arr]	,%%rax	\n\t"\
		"movq	%[__wtA]	,%%rbx		\n\t"\
		"movq	%[__wtB]	,%%rcx		\n\t"\
		"movaps		     (%%rbx),%%xmm2	\n\t	movaps		 0x10(%%rbx),%%xmm6	\n\t"\
		"movhpd		     (%%rcx),%%xmm3	\n\t	movhpd		-0x10(%%rcx),%%xmm7	\n\t"\
		"movlpd		 0x08(%%rcx),%%xmm3	\n\t	movlpd		-0x08(%%rcx),%%xmm7	\n\t"\
		"\n\t"\
		"movq	%%rsi,	%%rdi	\n\t	movq	%%rsi,	%%rbx	\n\t"\
		"shrq	$12,	%%rdi	\n\t	shrq	$14,	%%rbx	\n\t"\
		"andq	$0x0000000000000030	,%%rdi		\n\t	andq	$0x0000000000000030	,%%rbx		\n\t"\
		"movq	%%rsi,	%%rdx	\n\t	movq	%%rsi,	%%rcx	\n\t"\
		"shrq	$4 ,	%%rdx	\n\t	shrq	$6 ,	%%rcx	\n\t"\
		"andq	$0x0000000000000030	,%%rdx		\n\t	andq	$0x0000000000000030	,%%rcx		\n\t"\
		"addq	%%rax	,%%rdi			\n\t	addq	%%rax	,%%rbx			\n\t"\
		"addq	%%rax	,%%rdx			\n\t	addq	%%rax	,%%rcx			\n\t"\
		"\n\t"\
		"mulpd	 	0x100(%%rax),%%xmm2	\n\t	mulpd	 	0x100(%%rax),%%xmm6	\n\t"\
		"mulpd	 	0x110(%%rax),%%xmm3	\n\t	mulpd	 	0x110(%%rax),%%xmm7	\n\t"\
		"mulpd	 	     (%%rdi),%%xmm2	\n\t	mulpd	 	     (%%rbx),%%xmm6	\n\t"\
		"mulpd	 	0x040(%%rdx),%%xmm3	\n\t	mulpd	 	0x040(%%rcx),%%xmm7	\n\t"\
		"\n\t"\
		"movq	%[__cyA]	,%%rcx		\n\t	movq	%[__cyB]	,%%rdx		\n\t"\
		"mulpd		%%xmm3	,%%xmm1		\n\t	mulpd		%%xmm7	,%%xmm5		\n\t"\
		"addpd		(%%rcx)	,%%xmm1		\n\t	addpd		(%%rdx)	,%%xmm5		\n\t"\
		"movaps		%%xmm1	,%%xmm3		\n\t	movaps		%%xmm5	,%%xmm7		\n\t"\
		"roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd		$0,%%xmm7,%%xmm7	\n\t"\
		"\n\t"\
		"/* NO ROE HERE */		\n\t"\
		"\n\t"\
		"movq	%%rsi,	%%rdi	\n\t	movq	%%rsi,	%%rbx	\n\t"\
		"shrq	$20,	%%rdi	\n\t	shrq	$22,	%%rbx	\n\t"\
		"andq	$0x0000000000000030	,%%rdi		\n\t	andq	$0x0000000000000030	,%%rbx		\n\t"\
		"addq	%%rax	,%%rdi			\n\t	addq	%%rax	,%%rbx			\n\t"\
		"movaps		%%xmm3	,%%xmm1		\n\t	movaps		%%xmm7	,%%xmm5		\n\t"\
		"mulpd	 	0xc0(%%rdi)	,%%xmm3	\n\t	mulpd	 	0xc0(%%rbx)	,%%xmm7	\n\t"\
		"roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd		$0,%%xmm7,%%xmm7	\n\t"\
		"movaps		%%xmm3	,(%%rcx)	\n\t	movaps		%%xmm7	,(%%rdx)	\n\t"\
		"\n\t"\
		"\n\t"\
		"movq	%[__data]	,%%rax		\n\t"\
		"mulpd	 0x80(%%rdi)	,%%xmm3	\n\t	mulpd	 0x80(%%rbx)	,%%xmm7	\n\t"\
		"subpd		%%xmm3	,	%%xmm1	\n\t	subpd		%%xmm7	,	%%xmm5	\n\t"\
		"mulpd		%%xmm2	,	%%xmm1	\n\t	mulpd		%%xmm6	,	%%xmm5	\n\t"\
		"movaps		%%xmm1	,0x20(%%rax)\n\t	movaps		%%xmm5	,0x60(%%rax)\n\t"\
		"\n\t"\
		"\n\t"\
		"movq	%[__sse_bw]	,%%rax		\n\t"\
		"movq	%[__sse_nm1]	,%%rbx	\n\t"\
		"paddd		(%%rax)	,%%xmm0		\n\t"\
		"pand		(%%rbx)	,%%xmm0		\n\t"\
	"/**********************************************/\n\t"\
	"/*          Imaginary parts               */\n\t"\
	"/**********************************************/\n\t"\
		"movq	%[__sse_sw]	,%%rdx		\n\t"\
		"movaps	(%%rdx)	,%%xmm1			\n\t"\
		"psubd	%%xmm0	,%%xmm1			\n\t"\
		"movmskps	%%xmm1	,%%rsi		\n\t"\
		"shlq	$24	,%%rsi				\n\t"\
		"movaps	%%xmm0	,%%xmm1			\n\t"\
		"movslq	%[__n_minus_silp1],%%rcx\n\t"\
		"movd	%%rcx	,%%xmm2			\n\t"\
		"pshufd	$0	,%%xmm2	,%%xmm2		\n\t"\
		"psubd	%%xmm0	,%%xmm2			\n\t"\
		"movmskps	%%xmm2	,%%rcx		\n\t"\
		"shlq	$16	,%%rcx				\n\t"\
		"addq	%%rcx	,%%rsi			\n\t"\
		"movq	%[__sinwtm1]	,%%rdx	\n\t"\
		"movd	%%rdx	,%%xmm3			\n\t"\
		"pshufd	$0	,%%xmm3	,%%xmm3		\n\t"\
		"psubd	%%xmm3	,%%xmm1			\n\t"\
		"movmskps	%%xmm1	,%%rdx		\n\t"\
		"shlq	$8	,%%rdx				\n\t"\
		"addq	%%rdx	,%%rsi			\n\t"\
		"movq	%[__data]	,%%rax		\n\t"\
		"movaps	 0x30(%%rax)	,%%xmm1	\n\t"\
		"movaps	 0x70(%%rax)	,%%xmm5	\n\t"\
		"\n\t"\
		"movq	%[__half_arr]	,%%rax	\n\t"\
		"movq	%[__wtA]	,%%rbx		\n\t"\
		"movq	%[__wtB]	,%%rcx		\n\t"\
		"\n\t"\
		"movaps	     (%%rbx)	,%%xmm2	\n\t	movaps	 0x10(%%rbx)	,%%xmm6	\n\t"\
		"movhpd	     (%%rcx)	,%%xmm3	\n\t	movhpd	-0x10(%%rcx)	,%%xmm7	\n\t"\
		"movlpd	 0x08(%%rcx)	,%%xmm3	\n\t	movlpd	-0x08(%%rcx)	,%%xmm7	\n\t"\
		"\n\t"\
		"addq	$0x20	,%%rbx			\n\t"\
		"subq	$0x20	,%%rcx			\n\t"\
		"movq	%%rbx	,%[__wtA]		\n\t"\
		"movq	%%rcx	,%[__wtB]		\n\t"\
		"movq	%%rsi,	%%rdi	\n\t	movq	%%rsi,	%%rbx	\n\t"\
		"shrq	$12,	%%rdi	\n\t	shrq	$14,	%%rbx	\n\t"\
		"andq	$0x0000000000000030	,%%rdi		\n\t	andq	$0x0000000000000030	,%%rbx		\n\t"\
		"movq	%%rsi,	%%rdx	\n\t	movq	%%rsi,	%%rcx	\n\t"\
		"shrq	$4 ,	%%rdx	\n\t	shrq	$6 ,	%%rcx	\n\t"\
		"andq	$0x0000000000000030	,%%rdx		\n\t	andq	$0x0000000000000030	,%%rcx		\n\t"\
		"addq	%%rax	,%%rdi			\n\t	addq	%%rax	,%%rbx			\n\t"\
		"addq	%%rax	,%%rdx			\n\t	addq	%%rax	,%%rcx			\n\t"\
		"\n\t"\
		"mulpd	 0x120(%%rax)	,%%xmm2	\n\t	mulpd	 0x120(%%rax)	,%%xmm6	\n\t"\
		"mulpd	 0x130(%%rax)	,%%xmm3	\n\t	mulpd	 0x130(%%rax)	,%%xmm7	\n\t"\
		"mulpd	      (%%rdi)	,%%xmm2	\n\t	mulpd	      (%%rbx)	,%%xmm6	\n\t"\
		"mulpd	 0x040(%%rdx)	,%%xmm3	\n\t	mulpd	 0x040(%%rcx)	,%%xmm7	\n\t"\
		"\n\t"\
		"movq	%[__cyA]	,%%rcx		\n\t	movq	%[__cyB]	,%%rdx		\n\t"\
		"mulpd	%%xmm3	,%%xmm1			\n\t	mulpd	%%xmm7	,%%xmm5			\n\t"\
		"addpd	(%%rcx)	,%%xmm1			\n\t	addpd	(%%rdx)	,%%xmm5			\n\t"\
		"movaps	%%xmm1	,%%xmm3			\n\t	movaps	%%xmm5	,%%xmm7			\n\t"\
		"roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd		$0,%%xmm7,%%xmm7	\n\t"\
		"\n\t"\
		"/* NO ROE HERE */		\n\t"\
		"\n\t"\
		"movq	%%rsi,	%%rdi	\n\t	movq	%%rsi,	%%rbx	\n\t"\
		"shrq	$20,	%%rdi	\n\t	shrq	$22,	%%rbx	\n\t"\
		"andq	$0x0000000000000030	,%%rdi		\n\t	andq	$0x0000000000000030	,%%rbx		\n\t"\
		"addq	%%rax	,%%rdi			\n\t	addq	%%rax	,%%rbx			\n\t"\
		"movaps	%%xmm3	,%%xmm1			\n\t	movaps	%%xmm7	,%%xmm5			\n\t"\
		"mulpd	 0xc0(%%rdi)	,%%xmm3	\n\t	mulpd	 0xc0(%%rbx)	,%%xmm7	\n\t"\
		"roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd		$0,%%xmm7,%%xmm7	\n\t"\
		"movaps	%%xmm3	,(%%rcx)		\n\t	movaps	%%xmm7	,(%%rdx)		\n\t"\
		"\n\t"\
		"\n\t"\
		"movq	%[__data]	,%%rax		\n\t"\
		"mulpd	 0x80(%%rdi)	,%%xmm3	\n\t	mulpd	 0x80(%%rbx)	,%%xmm7	\n\t"\
		"subpd	%%xmm3	,%%xmm1			\n\t	subpd	%%xmm7	,%%xmm5			\n\t"\
		"mulpd	%%xmm2	,%%xmm1			\n\t	mulpd	%%xmm6	,%%xmm5			\n\t"\
		"movaps	%%xmm1	, 0x30(%%rax)	\n\t	movaps	%%xmm5	, 0x70(%%rax)	\n\t"\
		"\n\t"\
		"\n\t"\
		"movq	%[__sse_bw]	,%%rax		\n\t"\
		"movq	%[__sse_nm1]	,%%rbx	\n\t"\
		"paddd	(%%rax)	,%%xmm0			\n\t"\
		"pand	(%%rbx)	,%%xmm0			\n\t"\
		"movq	%[__bjmod_0],%%rcx		\n\t"\
		"movaps	%%xmm0,(%%rcx)			\n\t"\
	"/**********************************************/\n\t"\
	"/*          Repack the data:              */\n\t"\
	"/**********************************************/\n\t"\
		"mov	%[__data],%%rax			\n\t"\
		"movaps		0x10(%%rax)	,%%xmm1	\n\t	movaps		0x50(%%rax)	,%%xmm5	\n\t"\
		"movaps		    (%%rax)	,%%xmm0	\n\t	movaps		0x40(%%rax)	,%%xmm4	\n\t"\
		"movaps		%%xmm1		,%%xmm3	\n\t	movaps		%%xmm5		,%%xmm7	\n\t"\
		"movaps		%%xmm0		,%%xmm2	\n\t	movaps		%%xmm4		,%%xmm6	\n\t"\
		"unpckhpd	0x30(%%rax)	,%%xmm3	\n\t	unpckhpd	0x70(%%rax)	,%%xmm7	\n\t"\
		"unpcklpd	0x30(%%rax)	,%%xmm1	\n\t	unpcklpd	0x70(%%rax)	,%%xmm5	\n\t"\
		"movaps		%%xmm3,0x30(%%rax)	\n\t	movaps		%%xmm7,0x70(%%rax)	\n\t"\
		"unpckhpd	0x20(%%rax)	,%%xmm2	\n\t	unpckhpd	0x60(%%rax)	,%%xmm6	\n\t"\
		"unpcklpd	0x20(%%rax)	,%%xmm0	\n\t	unpcklpd	0x60(%%rax)	,%%xmm4	\n\t"\
		"movaps		%%xmm2,0x20(%%rax)	\n\t	movaps		%%xmm6,0x60(%%rax)	\n\t"\
		"movaps		%%xmm1,0x10(%%rax)	\n\t	movaps		%%xmm5,0x50(%%rax)	\n\t"\
		"movaps		%%xmm0,    (%%rax)	\n\t	movaps		%%xmm4,0x40(%%rax)	\n\t"\
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
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
	);\
	}


	/***************************************************************************************************************************************************/
	/********* Non-power-of-2-FFT versions of SSE2_cmplx_carry_norm_pow2_errcheck0_2B,1_2B,2_2B (only give sans-error-check version of latter 2: *******/
	/***************************************************************************************************************************************************/

	#define SSE2_cmplx_carry_norm_errcheck0_2B(Xdata,XwtA,XwtB,XwtC,XcyA,XcyB,Xbjmod_0,Xhalf_arr,Xi,Xn_minus_silp1,Xn_minus_sil,Xsign_mask,Xsinwt,Xsinwtm1,Xsse_bw,Xsse_n,Xsse_sw)\
	{\
	__asm__ volatile (\
	"/***************Unpack the data:*************************/\n\t"\
		"movq	%[__data]	,%%rax	\n\t"\
		"movaps		    (%%rax)	,%%xmm1	\n\t	movaps		0x40(%%rax)	,%%xmm5	\n\t"\
		"movaps		    (%%rax)	,%%xmm2	\n\t	movaps		0x40(%%rax)	,%%xmm6	\n\t"\
		"unpcklpd	0x20(%%rax)	,%%xmm1	\n\t	unpcklpd	0x60(%%rax)	,%%xmm5	\n\t"\
		"unpckhpd	0x20(%%rax)	,%%xmm2	\n\t	unpckhpd	0x60(%%rax)	,%%xmm6	\n\t"\
		"movaps		%%xmm2, 0x20(%%rax)	\n\t	movaps		%%xmm6, 0x60(%%rax)	\n\t"\
		"\n\t"\
		"movaps		0x10(%%rax)	,%%xmm2	\n\t	movaps		0x50(%%rax)	,%%xmm6	\n\t"\
		"movaps		0x10(%%rax)	,%%xmm3	\n\t	movaps		0x50(%%rax)	,%%xmm7	\n\t"\
		"unpcklpd	0x30(%%rax)	,%%xmm2	\n\t	unpcklpd	0x70(%%rax)	,%%xmm6	\n\t"\
		"unpckhpd	0x30(%%rax)	,%%xmm3	\n\t	unpckhpd	0x70(%%rax)	,%%xmm7	\n\t"\
		"movaps		%%xmm2, 0x10(%%rax)	\n\t	movaps		%%xmm6, 0x50(%%rax)	\n\t"\
		"movaps		%%xmm3, 0x30(%%rax)	\n\t	movaps		%%xmm7, 0x70(%%rax)	\n\t"\
	"/**********************************************/\n\t"\
	"/*          Real      parts                   */\n\t"\
	"/**********************************************/\n\t"\
		"movq	%[__bjmod_0],	%%rax	\n\t"\
		"movaps		(%%rax)	,	%%xmm0	\n\t"\
		"movq	%[__sse_sw]	,	%%rbx	\n\t"\
		"movaps		(%%rbx)	,	%%xmm7	\n\t"\
		"psubd		%%xmm0	,	%%xmm7	\n\t"\
		"movmskps	%%xmm7	,	%%rsi	\n\t"\
		"movslq	%[__i]	,%%rcx			\n\t"\
		"andq	$0xfffffffffffffffe,%%rsi		\n\t"\
		"addq	%%rcx	,%%rsi			\n\t"\
		"shlq	$24		,%%rsi			\n\t"\
		"movaps		%%xmm0	,%%xmm7		\n\t"\
		"movslq	%[__n_minus_sil],%%rcx	\n\t"\
		"movd	%%rcx	,%%xmm2			\n\t"\
		"pshufd	$0,	%%xmm2	,%%xmm2		\n\t"\
		"psubd		%%xmm0	,%%xmm2		\n\t"\
		"movmskps	%%xmm2	,%%rcx		\n\t"\
		"shlq	$16		,%%rcx			\n\t"\
		"addq	%%rcx	,%%rsi			\n\t"\
		"movslq	%[__sinwt]	,%%rdx		\n\t"\
		"movd	%%rdx	,%%xmm3			\n\t"\
		"pshufd	$0,	%%xmm3	,%%xmm3		\n\t"\
		"psubd		%%xmm3	,%%xmm7		\n\t"\
		"movmskps	%%xmm7	,%%rdx		\n\t"\
		"shlq	$8		,%%rdx			\n\t"\
		"addq	%%rdx	,%%rsi			\n\t"\
		"\n\t"\
		"movq	%[__half_arr]	,%%rax	\n\t"\
		"movq	%[__wtA]	,%%rbx		\n\t"\
		"movq	%[__wtB]	,%%rcx		\n\t"\
		"movaps		     (%%rbx),%%xmm2	\n\t	movaps		 0x10(%%rbx),%%xmm6	\n\t"\
		"movhpd		     (%%rcx),%%xmm3	\n\t	movhpd		-0x10(%%rcx),%%xmm7	\n\t"\
		"movlpd		 0x08(%%rcx),%%xmm3	\n\t	movlpd		-0x08(%%rcx),%%xmm7	\n\t"\
		"\n\t"\
		"subq	$0x20	,%%rcx			\n\t"\
		"movq	%%rcx	,%[__wtB]		\n\t"\
		"\n\t"\
		"movq	%%rsi,	%%rdi	\n\t	movq	%%rsi,	%%rbx	\n\t"\
		"shrq	$12,	%%rdi	\n\t	shrq	$14,	%%rbx	\n\t"\
		"andq	$0x0000000000000030	,%%rdi		\n\t	andq	$0x0000000000000030	,%%rbx		\n\t"\
		"movq	%%rsi,	%%rdx	\n\t	movq	%%rsi,	%%rcx	\n\t"\
		"shrq	$4 ,	%%rdx	\n\t	shrq	$6 ,	%%rcx	\n\t"\
		"andq	$0x0000000000000030	,%%rdx		\n\t	andq	$0x0000000000000030	,%%rcx		\n\t"\
		"addq	%%rax	,%%rdi			\n\t	addq	%%rax	,%%rbx			\n\t"\
		"addq	%%rax	,%%rdx			\n\t	addq	%%rax	,%%rcx			\n\t"\
		"\n\t"\
		"mulpd	 	0x100(%%rax),%%xmm2	\n\t	mulpd	 	0x100(%%rax),%%xmm6	\n\t"\
		"mulpd	 	0x110(%%rax),%%xmm3	\n\t	mulpd	 	0x110(%%rax),%%xmm7	\n\t"\
		"mulpd	 	     (%%rdi),%%xmm2	\n\t	mulpd	 	     (%%rbx),%%xmm6	\n\t"\
		"mulpd	 	0x040(%%rdx),%%xmm3	\n\t	mulpd	 	0x040(%%rcx),%%xmm7	\n\t"\
		"\n\t"\
		"movq	%[__cyA]	,%%rcx		\n\t	movq	%[__cyB]	,%%rdx		\n\t"\
		"mulpd		%%xmm3	,%%xmm1		\n\t	mulpd		%%xmm7	,%%xmm5		\n\t"\
		"addpd		(%%rcx)	,%%xmm1		\n\t	addpd		(%%rdx)	,%%xmm5		\n\t"\
		"movaps		%%xmm1	,%%xmm3		\n\t	movaps		%%xmm5	,%%xmm7		\n\t"\
		"roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd		$0,%%xmm7,%%xmm7	\n\t"\
		"\n\t"\
		"movq	%[__sign_mask],%%rbx	\n\t"\
		"subpd			%%xmm3	,%%xmm1	\n\t	subpd			%%xmm7	,%%xmm5	\n\t"\
		"andpd			(%%rbx)	,%%xmm1	\n\t	andpd			(%%rbx)	,%%xmm5	\n\t"\
		"maxpd			%%xmm5	,%%xmm1	\n\t"\
		"maxpd		-0x20(%%rax),%%xmm1	\n\t"\
		"movaps		%%xmm1,-0x20(%%rax)	\n\t"\
		"\n\t"\
		"movq	%%rsi,	%%rdi	\n\t	movq	%%rsi,	%%rbx	\n\t"\
		"shrq	$20,	%%rdi	\n\t	shrq	$22,	%%rbx	\n\t"\
		"andq	$0x0000000000000030	,%%rdi		\n\t	andq	$0x0000000000000030	,%%rbx		\n\t"\
		"addq	%%rax	,%%rdi			\n\t	addq	%%rax	,%%rbx			\n\t"\
		"movaps		%%xmm3	,%%xmm1		\n\t	movaps		%%xmm7	,%%xmm5		\n\t"\
		"mulpd	 	0xc0(%%rdi)	,%%xmm3	\n\t	mulpd	 	0xc0(%%rbx)	,%%xmm7	\n\t"\
		"roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd		$0,%%xmm7,%%xmm7	\n\t"\
		"movaps		%%xmm3	,(%%rcx)	\n\t	movaps		%%xmm7	,(%%rdx)	\n\t"\
		"\n\t"\
		"\n\t"\
		"movq	%[__data]	,%%rax		\n\t"\
		"mulpd	 0x80(%%rdi)	,%%xmm3	\n\t	mulpd	 0x80(%%rbx)	,%%xmm7	\n\t"\
		"subpd		%%xmm3	,	%%xmm1	\n\t	subpd		%%xmm7	,	%%xmm5	\n\t"\
		"mulpd		%%xmm2	,	%%xmm1	\n\t	mulpd		%%xmm6	,	%%xmm5	\n\t"\
		"movaps		%%xmm1	,    (%%rax)\n\t	movaps		%%xmm5	,0x40(%%rax)\n\t"\
		"\n\t"\
		"movq	%[__sse_n]	,%%rbx		\n\t"\
		"movaps		(%%rbx)	,%%xmm2		\n\t"\
		"movq	%[__sse_bw]	,%%rax		\n\t"\
		"paddd		(%%rax)	,%%xmm0		\n\t"\
		"movaps		%%xmm0	,%%xmm1		\n\t"\
		"pcmpgtd	%%xmm2	,%%xmm1		\n\t"\
		"pand		(%%rbx)	,%%xmm1		\n\t"\
		"psubd		%%xmm1	,%%xmm0		\n\t"\
	"/**********************************************/\n\t"\
	"/*          Imaginary parts               */\n\t"\
	"/**********************************************/\n\t"\
		"movq	%[__sse_sw]	,%%rdx		\n\t"\
		"movaps	(%%rdx)	,%%xmm1			\n\t"\
		"psubd	%%xmm0	,%%xmm1			\n\t"\
		"movmskps	%%xmm1	,%%rsi		\n\t"\
		"shlq	$24	,%%rsi				\n\t"\
		"movaps	%%xmm0	,%%xmm1			\n\t"\
		"movslq	%[__n_minus_silp1],%%rcx\n\t"\
		"movd	%%rcx	,%%xmm2			\n\t"\
		"pshufd	$0	,%%xmm2	,%%xmm2		\n\t"\
		"psubd	%%xmm0	,%%xmm2			\n\t"\
		"movmskps	%%xmm2	,%%rcx		\n\t"\
		"shlq	$16	,%%rcx				\n\t"\
		"addq	%%rcx	,%%rsi			\n\t"\
		"movslq	%[__sinwtm1]	,%%rdx	\n\t"\
		"movd	%%rdx	,%%xmm3			\n\t"\
		"pshufd	$0	,%%xmm3	,%%xmm3		\n\t"\
		"psubd	%%xmm3	,%%xmm1			\n\t"\
		"movmskps	%%xmm1	,%%rdx		\n\t"\
		"shlq	$8	,%%rdx				\n\t"\
		"addq	%%rdx	,%%rsi			\n\t"\
		"movq	%[__data]	,%%rax		\n\t"\
		"movaps	 0x10(%%rax)	,%%xmm1	\n\t"\
		"movaps	 0x50(%%rax)	,%%xmm5	\n\t"\
		"\n\t"\
		"movq	%[__half_arr]	,%%rax	\n\t"\
		"movq	%[__wtA]	,%%rbx		\n\t"\
		"movq	%[__wtC]	,%%rcx		\n\t"\
		"\n\t"\
		"movaps	     (%%rbx)	,%%xmm2	\n\t	movaps	 0x10(%%rbx)	,%%xmm6	\n\t"\
		"movhpd	     (%%rcx)	,%%xmm3	\n\t	movhpd	-0x10(%%rcx)	,%%xmm7	\n\t"\
		"movlpd	 0x08(%%rcx)	,%%xmm3	\n\t	movlpd	-0x08(%%rcx)	,%%xmm7	\n\t"\
		"\n\t"\
		"addq	$0x20	,%%rbx			\n\t"\
		"subq	$0x20	,%%rcx			\n\t"\
		"movq	%%rbx	,%[__wtA]		\n\t"\
		"movq	%%rcx	,%[__wtC]		\n\t"\
		"movq	%%rsi,	%%rdi	\n\t	movq	%%rsi,	%%rbx	\n\t"\
		"shrq	$12,	%%rdi	\n\t	shrq	$14,	%%rbx	\n\t"\
		"andq	$0x0000000000000030	,%%rdi		\n\t	andq	$0x0000000000000030	,%%rbx		\n\t"\
		"movq	%%rsi,	%%rdx	\n\t	movq	%%rsi,	%%rcx	\n\t"\
		"shrq	$4 ,	%%rdx	\n\t	shrq	$6 ,	%%rcx	\n\t"\
		"andq	$0x0000000000000030	,%%rdx		\n\t	andq	$0x0000000000000030	,%%rcx		\n\t"\
		"addq	%%rax	,%%rdi			\n\t	addq	%%rax	,%%rbx			\n\t"\
		"addq	%%rax	,%%rdx			\n\t	addq	%%rax	,%%rcx			\n\t"\
		"\n\t"\
		"mulpd	 0x120(%%rax)	,%%xmm2	\n\t	mulpd	 0x120(%%rax)	,%%xmm6	\n\t"\
		"mulpd	 0x130(%%rax)	,%%xmm3	\n\t	mulpd	 0x130(%%rax)	,%%xmm7	\n\t"\
		"mulpd	      (%%rdi)	,%%xmm2	\n\t	mulpd	      (%%rbx)	,%%xmm6	\n\t"\
		"mulpd	 0x040(%%rdx)	,%%xmm3	\n\t	mulpd	 0x040(%%rcx)	,%%xmm7	\n\t"\
		"\n\t"\
		"movq	%[__cyA]	,%%rcx		\n\t	movq	%[__cyB]	,%%rdx		\n\t"\
		"mulpd	%%xmm3	,%%xmm1			\n\t	mulpd	%%xmm7	,%%xmm5			\n\t"\
		"addpd	(%%rcx)	,%%xmm1			\n\t	addpd	(%%rdx)	,%%xmm5			\n\t"\
		"movaps	%%xmm1	,%%xmm3			\n\t	movaps	%%xmm5	,%%xmm7			\n\t"\
		"roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd		$0,%%xmm7,%%xmm7	\n\t"\
		"\n\t"\
		"movq	%[__sign_mask],%%rbx	\n\t"\
		"subpd	%%xmm3	,%%xmm1			\n\t	subpd	%%xmm7	,%%xmm5			\n\t"\
		"andpd	(%%rbx)	,%%xmm1			\n\t	andpd	(%%rbx)	,%%xmm5			\n\t"\
		"maxpd	%%xmm5	,%%xmm1			\n\t"\
		"maxpd	-0x20(%%rax)	,%%xmm1	\n\t"\
		"movaps	%%xmm1	,-0x20(%%rax)	\n\t"\
		"\n\t"\
		"movq	%%rsi,	%%rdi	\n\t	movq	%%rsi,	%%rbx	\n\t"\
		"shrq	$20,	%%rdi	\n\t	shrq	$22,	%%rbx	\n\t"\
		"andq	$0x0000000000000030	,%%rdi		\n\t	andq	$0x0000000000000030	,%%rbx		\n\t"\
		"addq	%%rax	,%%rdi			\n\t	addq	%%rax	,%%rbx			\n\t"\
		"movaps	%%xmm3	,%%xmm1			\n\t	movaps	%%xmm7	,%%xmm5			\n\t"\
		"mulpd	 0xc0(%%rdi)	,%%xmm3	\n\t	mulpd	 0xc0(%%rbx)	,%%xmm7	\n\t"\
		"roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd		$0,%%xmm7,%%xmm7	\n\t"\
		"movaps	%%xmm3	,(%%rcx)		\n\t	movaps	%%xmm7	,(%%rdx)		\n\t"\
		"\n\t"\
		"\n\t"\
		"movq	%[__data]	,%%rax		\n\t"\
		"mulpd	 0x80(%%rdi)	,%%xmm3	\n\t	mulpd	 0x80(%%rbx)	,%%xmm7	\n\t"\
		"subpd	%%xmm3	,%%xmm1			\n\t	subpd	%%xmm7	,%%xmm5			\n\t"\
		"mulpd	%%xmm2	,%%xmm1			\n\t	mulpd	%%xmm6	,%%xmm5			\n\t"\
		"movaps	%%xmm1	, 0x10(%%rax)	\n\t	movaps	%%xmm5	, 0x50(%%rax)	\n\t"\
		"\n\t"\
		"\n\t"\
		"movq	%[__sse_n]	,%%rbx		\n\t"\
		"movaps		(%%rbx)	,%%xmm2		\n\t"\
		"movq	%[__sse_bw]	,%%rax		\n\t"\
		"paddd		(%%rax)	,%%xmm0		\n\t"\
		"movaps		%%xmm0	,%%xmm1		\n\t"\
		"pcmpgtd	%%xmm2	,%%xmm1		\n\t"\
		"pand		(%%rbx)	,%%xmm1		\n\t"\
		"psubd		%%xmm1	,%%xmm0		\n\t"\
		"movq	%[__bjmod_0],%%rcx		\n\t"\
		"movaps	%%xmm0,(%%rcx)			\n\t"\
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
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","xmm0","xmm1","xmm2","xmm3","xmm5","xmm6","xmm7"		/* Clobbered registers */\
	);\
	}

	#define SSE2_cmplx_carry_norm_errcheck1_2B(Xdata,XwtA,XwtB,XwtC,XcyA,XcyB,Xbjmod_0,Xhalf_arr,Xn_minus_silp1,Xn_minus_sil,Xsign_mask,Xsinwt,Xsinwtm1,Xsse_bw,Xsse_n,Xsse_sw)\
	{\
	__asm__ volatile (\
	"/***************Unpack the data:*************************/\n\t"\
		"movq	%[__data]	,%%rax	\n\t"\
		"movaps		    (%%rax)	,%%xmm1	\n\t	movaps		0x40(%%rax)	,%%xmm5	\n\t"\
		"movaps		    (%%rax)	,%%xmm2	\n\t	movaps		0x40(%%rax)	,%%xmm6	\n\t"\
		"unpcklpd	0x20(%%rax)	,%%xmm1	\n\t	unpcklpd	0x60(%%rax)	,%%xmm5	\n\t"\
		"unpckhpd	0x20(%%rax)	,%%xmm2	\n\t	unpckhpd	0x60(%%rax)	,%%xmm6	\n\t"\
		"movaps		%%xmm2, 0x20(%%rax)	\n\t	movaps		%%xmm6, 0x60(%%rax)	\n\t"\
		"\n\t"\
		"movaps		0x10(%%rax)	,%%xmm2	\n\t	movaps		0x50(%%rax)	,%%xmm6	\n\t"\
		"movaps		0x10(%%rax)	,%%xmm3	\n\t	movaps		0x50(%%rax)	,%%xmm7	\n\t"\
		"unpcklpd	0x30(%%rax)	,%%xmm2	\n\t	unpcklpd	0x70(%%rax)	,%%xmm6	\n\t"\
		"unpckhpd	0x30(%%rax)	,%%xmm3	\n\t	unpckhpd	0x70(%%rax)	,%%xmm7	\n\t"\
		"movaps		%%xmm2, 0x10(%%rax)	\n\t	movaps		%%xmm6, 0x50(%%rax)	\n\t"\
		"movaps		%%xmm3, 0x30(%%rax)	\n\t	movaps		%%xmm7, 0x70(%%rax)	\n\t"\
	"/**********************************************/\n\t"\
	"/*          Real      parts                   */\n\t"\
	"/**********************************************/\n\t"\
		"movq	%[__bjmod_0],	%%rax	\n\t"\
		"movaps		(%%rax)	,	%%xmm0	\n\t"\
		"movq	%[__sse_sw]	,	%%rbx	\n\t"\
		"movaps		(%%rbx)	,	%%xmm7	\n\t"\
		"psubd		%%xmm0	,	%%xmm7	\n\t"\
		"movmskps	%%xmm7	,	%%rsi	\n\t"\
		"\n\t"\
		"shlq	$24		,%%rsi			\n\t"\
		"movaps		%%xmm0	,%%xmm7		\n\t"\
		"movslq	%[__n_minus_sil],%%rcx	\n\t"\
		"movd	%%rcx	,%%xmm2			\n\t"\
		"pshufd	$0,	%%xmm2	,%%xmm2		\n\t"\
		"psubd		%%xmm0	,%%xmm2		\n\t"\
		"movmskps	%%xmm2	,%%rcx		\n\t"\
		"shlq	$16		,%%rcx			\n\t"\
		"addq	%%rcx	,%%rsi			\n\t"\
		"movslq	%[__sinwt]	,%%rdx		\n\t"\
		"movd	%%rdx	,%%xmm3			\n\t"\
		"pshufd	$0,	%%xmm3	,%%xmm3		\n\t"\
		"psubd		%%xmm3	,%%xmm7		\n\t"\
		"movmskps	%%xmm7	,%%rdx		\n\t"\
		"shlq	$8		,%%rdx			\n\t"\
		"addq	%%rdx	,%%rsi			\n\t"\
		"\n\t"\
		"movq	%[__half_arr]	,%%rax	\n\t"\
		"movq	%[__wtA]	,%%rbx		\n\t"\
		"movq	%[__wtB]	,%%rcx		\n\t"\
		"movaps		     (%%rbx),%%xmm2	\n\t	movaps		 0x10(%%rbx),%%xmm6	\n\t"\
		"movhpd		     (%%rcx),%%xmm3	\n\t	movhpd		-0x10(%%rcx),%%xmm7	\n\t"\
		"movlpd		 0x08(%%rcx),%%xmm3	\n\t	movlpd		-0x08(%%rcx),%%xmm7	\n\t"\
		"\n\t"\
		"subq	$0x20	,%%rcx			\n\t"\
		"movq	%%rcx	,%[__wtB]		\n\t"\
		"\n\t"\
		"movq	%%rsi,	%%rdi	\n\t	movq	%%rsi,	%%rbx	\n\t"\
		"shrq	$12,	%%rdi	\n\t	shrq	$14,	%%rbx	\n\t"\
		"andq	$0x0000000000000030	,%%rdi		\n\t	andq	$0x0000000000000030	,%%rbx		\n\t"\
		"movq	%%rsi,	%%rdx	\n\t	movq	%%rsi,	%%rcx	\n\t"\
		"shrq	$4 ,	%%rdx	\n\t	shrq	$6 ,	%%rcx	\n\t"\
		"andq	$0x0000000000000030	,%%rdx		\n\t	andq	$0x0000000000000030	,%%rcx		\n\t"\
		"addq	%%rax	,%%rdi			\n\t	addq	%%rax	,%%rbx			\n\t"\
		"addq	%%rax	,%%rdx			\n\t	addq	%%rax	,%%rcx			\n\t"\
		"\n\t"\
		"mulpd	 	0x100(%%rax),%%xmm2	\n\t	mulpd	 	0x100(%%rax),%%xmm6	\n\t"\
		"mulpd	 	0x110(%%rax),%%xmm3	\n\t	mulpd	 	0x110(%%rax),%%xmm7	\n\t"\
		"mulpd	 	     (%%rdi),%%xmm2	\n\t	mulpd	 	     (%%rbx),%%xmm6	\n\t"\
		"mulpd	 	0x040(%%rdx),%%xmm3	\n\t	mulpd	 	0x040(%%rcx),%%xmm7	\n\t"\
		"\n\t"\
		"movq	%[__cyA]	,%%rcx		\n\t	movq	%[__cyB]	,%%rdx		\n\t"\
		"mulpd		%%xmm3	,%%xmm1		\n\t	mulpd		%%xmm7	,%%xmm5		\n\t"\
		"addpd		(%%rcx)	,%%xmm1		\n\t	addpd		(%%rdx)	,%%xmm5		\n\t"\
		"movaps		%%xmm1	,%%xmm3		\n\t	movaps		%%xmm5	,%%xmm7		\n\t"\
		"roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd		$0,%%xmm7,%%xmm7	\n\t"\
		"\n\t"\
		"movq	%[__sign_mask],%%rbx	\n\t"\
		"subpd	%%xmm3	,%%xmm1			\n\t	subpd	%%xmm7	,%%xmm5			\n\t"\
		"andpd	(%%rbx)	,%%xmm1			\n\t	andpd	(%%rbx)	,%%xmm5			\n\t"\
		"maxpd	%%xmm5	,%%xmm1			\n\t"\
		"maxpd	-0x20(%%rax)	,%%xmm1	\n\t"\
		"movaps	%%xmm1	,-0x20(%%rax)	\n\t"\
		"\n\t"\
		"movq	%%rsi,	%%rdi	\n\t	movq	%%rsi,	%%rbx	\n\t"\
		"shrq	$20,	%%rdi	\n\t	shrq	$22,	%%rbx	\n\t"\
		"andq	$0x0000000000000030	,%%rdi		\n\t	andq	$0x0000000000000030	,%%rbx		\n\t"\
		"addq	%%rax	,%%rdi			\n\t	addq	%%rax	,%%rbx			\n\t"\
		"movaps		%%xmm3	,%%xmm1		\n\t	movaps		%%xmm7	,%%xmm5		\n\t"\
		"mulpd	 	0xc0(%%rdi)	,%%xmm3	\n\t	mulpd	 	0xc0(%%rbx)	,%%xmm7	\n\t"\
		"roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd		$0,%%xmm7,%%xmm7	\n\t"\
		"movaps		%%xmm3	,(%%rcx)	\n\t	movaps		%%xmm7	,(%%rdx)	\n\t"\
		"\n\t"\
		"\n\t"\
		"movq	%[__data]	,%%rax		\n\t"\
		"mulpd	 0x80(%%rdi)	,%%xmm3	\n\t	mulpd	 0x80(%%rbx)	,%%xmm7	\n\t"\
		"subpd		%%xmm3	,	%%xmm1	\n\t	subpd		%%xmm7	,	%%xmm5	\n\t"\
		"mulpd		%%xmm2	,	%%xmm1	\n\t	mulpd		%%xmm6	,	%%xmm5	\n\t"\
		"movaps		%%xmm1	,    (%%rax)\n\t	movaps		%%xmm5	,0x40(%%rax)\n\t"\
		"\n\t"\
		"movq	%[__sse_n]	,%%rbx		\n\t"\
		"movaps		(%%rbx)	,%%xmm2		\n\t"\
		"movq	%[__sse_bw]	,%%rax		\n\t"\
		"paddd		(%%rax)	,%%xmm0		\n\t"\
		"movaps		%%xmm0	,%%xmm1		\n\t"\
		"pcmpgtd	%%xmm2	,%%xmm1		\n\t"\
		"pand		(%%rbx)	,%%xmm1		\n\t"\
		"psubd		%%xmm1	,%%xmm0		\n\t"\
	"/**********************************************/\n\t"\
	"/*          Imaginary parts               */\n\t"\
	"/**********************************************/\n\t"\
		"movq	%[__sse_sw]	,%%rdx		\n\t"\
		"movaps	(%%rdx)	,%%xmm1			\n\t"\
		"psubd	%%xmm0	,%%xmm1			\n\t"\
		"movmskps	%%xmm1	,%%rsi		\n\t"\
		"shlq	$24	,%%rsi				\n\t"\
		"movaps	%%xmm0	,%%xmm1			\n\t"\
		"movslq	%[__n_minus_silp1],%%rcx\n\t"\
		"movd	%%rcx	,%%xmm2			\n\t"\
		"pshufd	$0	,%%xmm2	,%%xmm2		\n\t"\
		"psubd	%%xmm0	,%%xmm2			\n\t"\
		"movmskps	%%xmm2	,%%rcx		\n\t"\
		"shlq	$16	,%%rcx				\n\t"\
		"addq	%%rcx	,%%rsi			\n\t"\
		"movslq	%[__sinwtm1]	,%%rdx	\n\t"\
		"movd	%%rdx	,%%xmm3			\n\t"\
		"pshufd	$0	,%%xmm3	,%%xmm3		\n\t"\
		"psubd	%%xmm3	,%%xmm1			\n\t"\
		"movmskps	%%xmm1	,%%rdx		\n\t"\
		"shlq	$8	,%%rdx				\n\t"\
		"addq	%%rdx	,%%rsi			\n\t"\
		"movq	%[__data]	,%%rax		\n\t"\
		"movaps	 0x10(%%rax)	,%%xmm1	\n\t"\
		"movaps	 0x50(%%rax)	,%%xmm5	\n\t"\
		"\n\t"\
		"movq	%[__half_arr]	,%%rax	\n\t"\
		"movq	%[__wtA]	,%%rbx		\n\t"\
		"movq	%[__wtC]	,%%rcx		\n\t"\
		"\n\t"\
		"movaps	     (%%rbx)	,%%xmm2	\n\t	movaps	 0x10(%%rbx)	,%%xmm6	\n\t"\
		"movhpd	     (%%rcx)	,%%xmm3	\n\t	movhpd	-0x10(%%rcx)	,%%xmm7	\n\t"\
		"movlpd	 0x08(%%rcx)	,%%xmm3	\n\t	movlpd	-0x08(%%rcx)	,%%xmm7	\n\t"\
		"\n\t"\
		"addq	$0x20	,%%rbx			\n\t"\
		"subq	$0x20	,%%rcx			\n\t"\
		"movq	%%rbx	,%[__wtA]		\n\t"\
		"movq	%%rcx	,%[__wtC]		\n\t"\
		"movq	%%rsi,	%%rdi	\n\t	movq	%%rsi,	%%rbx	\n\t"\
		"shrq	$12,	%%rdi	\n\t	shrq	$14,	%%rbx	\n\t"\
		"andq	$0x0000000000000030	,%%rdi		\n\t	andq	$0x0000000000000030	,%%rbx		\n\t"\
		"movq	%%rsi,	%%rdx	\n\t	movq	%%rsi,	%%rcx	\n\t"\
		"shrq	$4 ,	%%rdx	\n\t	shrq	$6 ,	%%rcx	\n\t"\
		"andq	$0x0000000000000030	,%%rdx		\n\t	andq	$0x0000000000000030	,%%rcx		\n\t"\
		"addq	%%rax	,%%rdi			\n\t	addq	%%rax	,%%rbx			\n\t"\
		"addq	%%rax	,%%rdx			\n\t	addq	%%rax	,%%rcx			\n\t"\
		"\n\t"\
		"mulpd	 0x120(%%rax)	,%%xmm2	\n\t	mulpd	 0x120(%%rax)	,%%xmm6	\n\t"\
		"mulpd	 0x130(%%rax)	,%%xmm3	\n\t	mulpd	 0x130(%%rax)	,%%xmm7	\n\t"\
		"mulpd	      (%%rdi)	,%%xmm2	\n\t	mulpd	      (%%rbx)	,%%xmm6	\n\t"\
		"mulpd	 0x040(%%rdx)	,%%xmm3	\n\t	mulpd	 0x040(%%rcx)	,%%xmm7	\n\t"\
		"\n\t"\
		"movq	%[__cyA]	,%%rcx		\n\t	movq	%[__cyB]	,%%rdx		\n\t"\
		"mulpd	%%xmm3	,%%xmm1			\n\t	mulpd	%%xmm7	,%%xmm5			\n\t"\
		"addpd	(%%rcx)	,%%xmm1			\n\t	addpd	(%%rdx)	,%%xmm5			\n\t"\
		"movaps	%%xmm1	,%%xmm3			\n\t	movaps	%%xmm5	,%%xmm7			\n\t"\
		"roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd		$0,%%xmm7,%%xmm7	\n\t"\
		"\n\t"\
		"movq	%[__sign_mask],%%rbx	\n\t"\
		"subpd	%%xmm3	,%%xmm1			\n\t	subpd	%%xmm7	,%%xmm5			\n\t"\
		"andpd	(%%rbx)	,%%xmm1			\n\t	andpd	(%%rbx)	,%%xmm5			\n\t"\
		"maxpd	%%xmm5	,%%xmm1			\n\t"\
		"maxpd	-0x20(%%rax)	,%%xmm1	\n\t"\
		"movaps	%%xmm1	,-0x20(%%rax)	\n\t"\
		"\n\t"\
		"movq	%%rsi,	%%rdi	\n\t	movq	%%rsi,	%%rbx	\n\t"\
		"shrq	$20,	%%rdi	\n\t	shrq	$22,	%%rbx	\n\t"\
		"andq	$0x0000000000000030	,%%rdi		\n\t	andq	$0x0000000000000030	,%%rbx		\n\t"\
		"addq	%%rax	,%%rdi			\n\t	addq	%%rax	,%%rbx			\n\t"\
		"movaps	%%xmm3	,%%xmm1			\n\t	movaps	%%xmm7	,%%xmm5			\n\t"\
		"mulpd	 0xc0(%%rdi)	,%%xmm3	\n\t	mulpd	 0xc0(%%rbx)	,%%xmm7	\n\t"\
		"roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd		$0,%%xmm7,%%xmm7	\n\t"\
		"movaps	%%xmm3	,(%%rcx)		\n\t	movaps	%%xmm7	,(%%rdx)		\n\t"\
		"\n\t"\
		"movq	%[__data]	,%%rax		\n\t"\
		"mulpd	 0x80(%%rdi)	,%%xmm3	\n\t	mulpd	 0x80(%%rbx)	,%%xmm7	\n\t"\
		"subpd	%%xmm3	,%%xmm1			\n\t	subpd	%%xmm7	,%%xmm5			\n\t"\
		"mulpd	%%xmm2	,%%xmm1			\n\t	mulpd	%%xmm6	,%%xmm5			\n\t"\
		"movaps	%%xmm1	, 0x10(%%rax)	\n\t	movaps	%%xmm5	, 0x50(%%rax)	\n\t"\
		"\n\t"\
		"movq	%[__sse_n]	,%%rbx		\n\t"\
		"movaps		(%%rbx)	,%%xmm2		\n\t"\
		"movq	%[__sse_bw]	,%%rax		\n\t"\
		"paddd		(%%rax)	,%%xmm0		\n\t"\
		"movaps		%%xmm0	,%%xmm1		\n\t"\
		"pcmpgtd	%%xmm2	,%%xmm1		\n\t"\
		"pand		(%%rbx)	,%%xmm1		\n\t"\
		"psubd		%%xmm1	,%%xmm0		\n\t"\
		"movq	%[__bjmod_0],%%rcx		\n\t"\
		"movaps	%%xmm0,(%%rcx)			\n\t"\
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
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","xmm0","xmm1","xmm2","xmm3","xmm5","xmm6","xmm7"		/* Clobbered registers */\
	);\
	}

	#define SSE2_cmplx_carry_norm_errcheck2_2B(Xdata,XwtA,XwtB,XcyA,XcyB,Xbjmod_0,Xhalf_arr,Xn_minus_silp1,Xn_minus_sil,Xsign_mask,Xsinwt,Xsinwtm1,Xsse_bw,Xsse_n,Xsse_sw)\
	{\
	__asm__ volatile (\
	"/**********************************************/\n\t"\
	"/*          Real      parts                   */\n\t"\
	"/**********************************************/\n\t"\
		"movq	%[__bjmod_0],	%%rax	\n\t"\
		"movaps		(%%rax)	,	%%xmm0	\n\t"\
		"movq	%[__sse_sw]	,	%%rbx	\n\t"\
		"movaps		(%%rbx)	,	%%xmm1	\n\t"\
		"psubd		%%xmm0	,	%%xmm1	\n\t"\
		"movmskps	%%xmm1	,	%%rsi	\n\t"\
		"\n\t"\
		"shlq	$24		,%%rsi			\n\t"\
		"movaps		%%xmm0	,%%xmm1		\n\t"\
		"movslq	%[__n_minus_sil],%%rcx	\n\t"\
		"movd	%%rcx	,%%xmm2			\n\t"\
		"pshufd	$0,	%%xmm2	,%%xmm2		\n\t"\
		"psubd		%%xmm0	,%%xmm2		\n\t"\
		"movmskps	%%xmm2	,%%rcx		\n\t"\
		"shlq	$16		,%%rcx			\n\t"\
		"addq	%%rcx	,%%rsi			\n\t"\
		"movslq	%[__sinwt]	,%%rdx		\n\t"\
		"movd	%%rdx	,%%xmm3			\n\t"\
		"pshufd	$0,	%%xmm3	,%%xmm3		\n\t"\
		"psubd		%%xmm3	,%%xmm1		\n\t"\
		"movmskps	%%xmm1	,%%rdx		\n\t"\
		"shlq	$8		,%%rdx			\n\t"\
		"addq	%%rdx	,%%rsi			\n\t"\
		"\n\t"\
		"mov	%[__data],%%rax			\n\t"\
		"movaps	 0x20(%%rax),%%xmm1		\n\t	movaps		 0x60(%%rax),%%xmm5\n\t"\
		"\n\t"\
		"movq	%[__half_arr]	,%%rax	\n\t"\
		"movq	%[__wtA]	,%%rbx		\n\t"\
		"movq	%[__wtB]	,%%rcx		\n\t"\
		"movaps		     (%%rbx),%%xmm2	\n\t	movaps		 0x10(%%rbx),%%xmm6	\n\t"\
		"movhpd		     (%%rcx),%%xmm3	\n\t	movhpd		-0x10(%%rcx),%%xmm7	\n\t"\
		"movlpd		 0x08(%%rcx),%%xmm3	\n\t	movlpd		-0x08(%%rcx),%%xmm7	\n\t"\
		"\n\t"\
		"movq	%%rsi,	%%rdi	\n\t	movq	%%rsi,	%%rbx	\n\t"\
		"shrq	$12,	%%rdi	\n\t	shrq	$14,	%%rbx	\n\t"\
		"andq	$0x0000000000000030	,%%rdi		\n\t	andq	$0x0000000000000030	,%%rbx		\n\t"\
		"movq	%%rsi,	%%rdx	\n\t	movq	%%rsi,	%%rcx	\n\t"\
		"shrq	$4 ,	%%rdx	\n\t	shrq	$6 ,	%%rcx	\n\t"\
		"andq	$0x0000000000000030	,%%rdx		\n\t	andq	$0x0000000000000030	,%%rcx		\n\t"\
		"addq	%%rax	,%%rdi			\n\t	addq	%%rax	,%%rbx			\n\t"\
		"addq	%%rax	,%%rdx			\n\t	addq	%%rax	,%%rcx			\n\t"\
		"\n\t"\
		"mulpd	 	0x100(%%rax),%%xmm2	\n\t	mulpd	 	0x100(%%rax),%%xmm6	\n\t"\
		"mulpd	 	0x110(%%rax),%%xmm3	\n\t	mulpd	 	0x110(%%rax),%%xmm7	\n\t"\
		"mulpd	 	     (%%rdi),%%xmm2	\n\t	mulpd	 	     (%%rbx),%%xmm6	\n\t"\
		"mulpd	 	0x040(%%rdx),%%xmm3	\n\t	mulpd	 	0x040(%%rcx),%%xmm7	\n\t"\
		"\n\t"\
		"movq	%[__cyA]	,%%rcx		\n\t	movq	%[__cyB]	,%%rdx		\n\t"\
		"mulpd		%%xmm3	,%%xmm1		\n\t	mulpd		%%xmm7	,%%xmm5		\n\t"\
		"addpd		(%%rcx)	,%%xmm1		\n\t	addpd		(%%rdx)	,%%xmm5		\n\t"\
		"movaps		%%xmm1	,%%xmm3		\n\t	movaps		%%xmm5	,%%xmm7		\n\t"\
		"roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd		$0,%%xmm7,%%xmm7	\n\t"\
		"\n\t"\
		"/* NO ROE HERE */				\n\t"\
		"\n\t"\
		"movq	%%rsi,	%%rdi	\n\t	movq	%%rsi,	%%rbx	\n\t"\
		"shrq	$20,	%%rdi	\n\t	shrq	$22,	%%rbx	\n\t"\
		"andq	$0x0000000000000030	,%%rdi		\n\t	andq	$0x0000000000000030	,%%rbx		\n\t"\
		"addq	%%rax	,%%rdi			\n\t	addq	%%rax	,%%rbx			\n\t"\
		"movaps		%%xmm3	,%%xmm1		\n\t	movaps		%%xmm7	,%%xmm5		\n\t"\
		"mulpd	 	0xc0(%%rdi)	,%%xmm3	\n\t	mulpd	 	0xc0(%%rbx)	,%%xmm7	\n\t"\
		"roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd		$0,%%xmm7,%%xmm7	\n\t"\
		"movaps		%%xmm3	,(%%rcx)	\n\t	movaps		%%xmm7	,(%%rdx)	\n\t"\
		"\n\t"\
		"\n\t"\
		"movq	%[__data]	,%%rax		\n\t"\
		"mulpd	 0x80(%%rdi)	,%%xmm3	\n\t	mulpd	 0x80(%%rbx)	,%%xmm7	\n\t"\
		"subpd		%%xmm3	,	%%xmm1	\n\t	subpd		%%xmm7	,	%%xmm5	\n\t"\
		"mulpd		%%xmm2	,	%%xmm1	\n\t	mulpd		%%xmm6	,	%%xmm5	\n\t"\
		"movaps		%%xmm1	,0x20(%%rax)\n\t	movaps		%%xmm5	,0x60(%%rax)\n\t"\
		"\n\t"\
		"movq	%[__sse_n]	,%%rbx		\n\t"\
		"movaps		(%%rbx)	,%%xmm2		\n\t"\
		"movq	%[__sse_bw]	,%%rax		\n\t"\
		"paddd		(%%rax)	,%%xmm0		\n\t"\
		"movaps		%%xmm0	,%%xmm1		\n\t"\
		"pcmpgtd	%%xmm2	,%%xmm1		\n\t"\
		"pand		(%%rbx)	,%%xmm1		\n\t"\
		"psubd		%%xmm1	,%%xmm0		\n\t"\
	"/**********************************************/\n\t"\
	"/*          Imaginary parts               */\n\t"\
	"/**********************************************/\n\t"\
		"movq	%[__sse_sw]	,%%rdx		\n\t"\
		"movaps	(%%rdx)	,%%xmm1			\n\t"\
		"psubd	%%xmm0	,%%xmm1			\n\t"\
		"movmskps	%%xmm1	,%%rsi		\n\t"\
		"shlq	$24	,%%rsi				\n\t"\
		"movaps	%%xmm0	,%%xmm1			\n\t"\
		"movslq	%[__n_minus_silp1],%%rcx\n\t"\
		"movd	%%rcx	,%%xmm2			\n\t"\
		"pshufd	$0	,%%xmm2	,%%xmm2		\n\t"\
		"psubd	%%xmm0	,%%xmm2			\n\t"\
		"movmskps	%%xmm2	,%%rcx		\n\t"\
		"shlq	$16	,%%rcx				\n\t"\
		"addq	%%rcx	,%%rsi			\n\t"\
		"movslq	%[__sinwtm1]	,%%rdx	\n\t"\
		"movd	%%rdx	,%%xmm3			\n\t"\
		"pshufd	$0	,%%xmm3	,%%xmm3		\n\t"\
		"psubd	%%xmm3	,%%xmm1			\n\t"\
		"movmskps	%%xmm1	,%%rdx		\n\t"\
		"shlq	$8	,%%rdx				\n\t"\
		"addq	%%rdx	,%%rsi			\n\t"\
		"movq	%[__data]	,%%rax		\n\t"\
		"movaps	 0x30(%%rax)	,%%xmm1	\n\t"\
		"movaps	 0x70(%%rax)	,%%xmm5	\n\t"\
		"\n\t"\
		"movq	%[__half_arr]	,%%rax	\n\t"\
		"movq	%[__wtA]	,%%rbx		\n\t"\
		"movq	%[__wtB]	,%%rcx		\n\t"\
		"\n\t"\
		"movaps	     (%%rbx)	,%%xmm2	\n\t	movaps	 0x10(%%rbx)	,%%xmm6	\n\t"\
		"movhpd	     (%%rcx)	,%%xmm3	\n\t	movhpd	-0x10(%%rcx)	,%%xmm7	\n\t"\
		"movlpd	 0x08(%%rcx)	,%%xmm3	\n\t	movlpd	-0x08(%%rcx)	,%%xmm7	\n\t"\
		"\n\t"\
		"addq	$0x20	,%%rbx			\n\t"\
		"subq	$0x20	,%%rcx			\n\t"\
		"movq	%%rbx	,%[__wtA]		\n\t"\
		"movq	%%rcx	,%[__wtB]		\n\t"\
		"movq	%%rsi,	%%rdi	\n\t	movq	%%rsi,	%%rbx	\n\t"\
		"shrq	$12,	%%rdi	\n\t	shrq	$14,	%%rbx	\n\t"\
		"andq	$0x0000000000000030	,%%rdi		\n\t	andq	$0x0000000000000030	,%%rbx		\n\t"\
		"movq	%%rsi,	%%rdx	\n\t	movq	%%rsi,	%%rcx	\n\t"\
		"shrq	$4 ,	%%rdx	\n\t	shrq	$6 ,	%%rcx	\n\t"\
		"andq	$0x0000000000000030	,%%rdx		\n\t	andq	$0x0000000000000030	,%%rcx		\n\t"\
		"addq	%%rax	,%%rdi			\n\t	addq	%%rax	,%%rbx			\n\t"\
		"addq	%%rax	,%%rdx			\n\t	addq	%%rax	,%%rcx			\n\t"\
		"\n\t"\
		"mulpd	 0x120(%%rax)	,%%xmm2	\n\t	mulpd	 0x120(%%rax)	,%%xmm6	\n\t"\
		"mulpd	 0x130(%%rax)	,%%xmm3	\n\t	mulpd	 0x130(%%rax)	,%%xmm7	\n\t"\
		"mulpd	      (%%rdi)	,%%xmm2	\n\t	mulpd	      (%%rbx)	,%%xmm6	\n\t"\
		"mulpd	 0x040(%%rdx)	,%%xmm3	\n\t	mulpd	 0x040(%%rcx)	,%%xmm7	\n\t"\
		"\n\t"\
		"movq	%[__cyA]	,%%rcx		\n\t	movq	%[__cyB]	,%%rdx		\n\t"\
		"mulpd	%%xmm3	,%%xmm1			\n\t	mulpd	%%xmm7	,%%xmm5			\n\t"\
		"addpd	(%%rcx)	,%%xmm1			\n\t	addpd	(%%rdx)	,%%xmm5			\n\t"\
		"movaps	%%xmm1	,%%xmm3			\n\t	movaps	%%xmm5	,%%xmm7			\n\t"\
		"roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd		$0,%%xmm7,%%xmm7	\n\t"\
		"\n\t"\
		"/* NO ROE HERE */				\n\t"\
		"\n\t"\
		"movq	%%rsi,	%%rdi	\n\t	movq	%%rsi,	%%rbx	\n\t"\
		"shrq	$20,	%%rdi	\n\t	shrq	$22,	%%rbx	\n\t"\
		"andq	$0x0000000000000030	,%%rdi		\n\t	andq	$0x0000000000000030	,%%rbx		\n\t"\
		"addq	%%rax	,%%rdi			\n\t	addq	%%rax	,%%rbx			\n\t"\
		"movaps	%%xmm3	,%%xmm1			\n\t	movaps	%%xmm7	,%%xmm5			\n\t"\
		"mulpd	 0xc0(%%rdi)	,%%xmm3	\n\t	mulpd	 0xc0(%%rbx)	,%%xmm7	\n\t"\
		"roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd		$0,%%xmm7,%%xmm7	\n\t"\
		"movaps	%%xmm3	,(%%rcx)		\n\t	movaps	%%xmm7	,(%%rdx)		\n\t"\
		"\n\t"\
		"\n\t"\
		"movq	%[__data]	,%%rax		\n\t"\
		"mulpd	 0x80(%%rdi)	,%%xmm3	\n\t	mulpd	 0x80(%%rbx)	,%%xmm7	\n\t"\
		"subpd	%%xmm3	,%%xmm1			\n\t	subpd	%%xmm7	,%%xmm5			\n\t"\
		"mulpd	%%xmm2	,%%xmm1			\n\t	mulpd	%%xmm6	,%%xmm5			\n\t"\
		"movaps	%%xmm1	, 0x30(%%rax)	\n\t	movaps	%%xmm5	, 0x70(%%rax)	\n\t"\
		"\n\t"\
		"movq	%[__sse_n]	,%%rbx		\n\t"\
		"movaps		(%%rbx)	,%%xmm2		\n\t"\
		"movq	%[__sse_bw]	,%%rax		\n\t"\
		"paddd		(%%rax)	,%%xmm0		\n\t"\
		"movaps		%%xmm0	,%%xmm1		\n\t"\
		"pcmpgtd	%%xmm2	,%%xmm1		\n\t"\
		"pand		(%%rbx)	,%%xmm1		\n\t"\
		"psubd		%%xmm1	,%%xmm0		\n\t"\
		"movq	%[__bjmod_0],%%rcx		\n\t"\
		"movaps	%%xmm0,(%%rcx)			\n\t"\
	"/**********************************************/\n\t"\
	"/*          Repack the data:              */\n\t"\
	"/**********************************************/\n\t"\
		"mov	%[__data],%%rax			\n\t"\
		"movaps		0x10(%%rax)	,%%xmm1	\n\t	movaps		0x50(%%rax)	,%%xmm5	\n\t"\
		"movaps		    (%%rax)	,%%xmm0	\n\t	movaps		0x40(%%rax)	,%%xmm4	\n\t"\
		"movaps		%%xmm1		,%%xmm3	\n\t	movaps		%%xmm5		,%%xmm7	\n\t"\
		"movaps		%%xmm0		,%%xmm2	\n\t	movaps		%%xmm4		,%%xmm6	\n\t"\
		"unpckhpd	0x30(%%rax)	,%%xmm3	\n\t	unpckhpd	0x70(%%rax)	,%%xmm7	\n\t"\
		"unpcklpd	0x30(%%rax)	,%%xmm1	\n\t	unpcklpd	0x70(%%rax)	,%%xmm5	\n\t"\
		"movaps		%%xmm3,0x30(%%rax)	\n\t	movaps		%%xmm7,0x70(%%rax)	\n\t"\
		"unpckhpd	0x20(%%rax)	,%%xmm2	\n\t	unpckhpd	0x60(%%rax)	,%%xmm6	\n\t"\
		"unpcklpd	0x20(%%rax)	,%%xmm0	\n\t	unpcklpd	0x60(%%rax)	,%%xmm4	\n\t"\
		"movaps		%%xmm2,0x20(%%rax)	\n\t	movaps		%%xmm6,0x60(%%rax)	\n\t"\
		"movaps		%%xmm1,0x10(%%rax)	\n\t	movaps		%%xmm5,0x50(%%rax)	\n\t"\
		"movaps		%%xmm0,    (%%rax)	\n\t	movaps		%%xmm4,0x40(%%rax)	\n\t"\
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
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
	);\
	}

	/******************************************************************************************************************************************************************/
	/********** No-ROE-Check versions of the latter 2 of the above 3 macros - the first is called too infrequently to bother with a special non-ROE version: **********/
	/******************************************************************************************************************************************************************/

	#define SSE2_cmplx_carry_norm_nocheck1_2B(Xdata,XwtA,XwtB,XwtC,XcyA,XcyB,Xbjmod_0,Xhalf_arr,Xn_minus_silp1,Xn_minus_sil,Xsinwt,Xsinwtm1,Xsse_bw,Xsse_n,Xsse_sw)\
	{\
	__asm__ volatile (\
	"/***************Unpack the data:*************************/\n\t"\
		"movq	%[__data]	,%%rax	\n\t"\
		"movaps		    (%%rax)	,%%xmm1	\n\t	movaps		0x40(%%rax)	,%%xmm5	\n\t"\
		"movaps		    (%%rax)	,%%xmm2	\n\t	movaps		0x40(%%rax)	,%%xmm6	\n\t"\
		"unpcklpd	0x20(%%rax)	,%%xmm1	\n\t	unpcklpd	0x60(%%rax)	,%%xmm5	\n\t"\
		"unpckhpd	0x20(%%rax)	,%%xmm2	\n\t	unpckhpd	0x60(%%rax)	,%%xmm6	\n\t"\
		"movaps		%%xmm2, 0x20(%%rax)	\n\t	movaps		%%xmm6, 0x60(%%rax)	\n\t"\
		"\n\t"\
		"movaps		0x10(%%rax)	,%%xmm2	\n\t	movaps		0x50(%%rax)	,%%xmm6	\n\t"\
		"movaps		0x10(%%rax)	,%%xmm3	\n\t	movaps		0x50(%%rax)	,%%xmm7	\n\t"\
		"unpcklpd	0x30(%%rax)	,%%xmm2	\n\t	unpcklpd	0x70(%%rax)	,%%xmm6	\n\t"\
		"unpckhpd	0x30(%%rax)	,%%xmm3	\n\t	unpckhpd	0x70(%%rax)	,%%xmm7	\n\t"\
		"movaps		%%xmm2, 0x10(%%rax)	\n\t	movaps		%%xmm6, 0x50(%%rax)	\n\t"\
		"movaps		%%xmm3, 0x30(%%rax)	\n\t	movaps		%%xmm7, 0x70(%%rax)	\n\t"\
	"/**********************************************/\n\t"\
	"/*          Real      parts                   */\n\t"\
	"/**********************************************/\n\t"\
		"movq	%[__bjmod_0],	%%rax	\n\t"\
		"movaps		(%%rax)	,	%%xmm0	\n\t"\
		"movq	%[__sse_sw]	,	%%rbx	\n\t"\
		"movaps		(%%rbx)	,	%%xmm7	\n\t"\
		"psubd		%%xmm0	,	%%xmm7	\n\t"\
		"movmskps	%%xmm7	,	%%rsi	\n\t"\
		"\n\t"\
		"shlq	$24		,%%rsi			\n\t"\
		"movaps		%%xmm0	,%%xmm7		\n\t"\
		"movslq	%[__n_minus_sil],%%rcx	\n\t"\
		"movd	%%rcx	,%%xmm2			\n\t"\
		"pshufd	$0,	%%xmm2	,%%xmm2		\n\t"\
		"psubd		%%xmm0	,%%xmm2		\n\t"\
		"movmskps	%%xmm2	,%%rcx		\n\t"\
		"shlq	$16		,%%rcx			\n\t"\
		"addq	%%rcx	,%%rsi			\n\t"\
		"movslq	%[__sinwt]	,%%rdx		\n\t"\
		"movd	%%rdx	,%%xmm3			\n\t"\
		"pshufd	$0,	%%xmm3	,%%xmm3		\n\t"\
		"psubd		%%xmm3	,%%xmm7		\n\t"\
		"movmskps	%%xmm7	,%%rdx		\n\t"\
		"shlq	$8		,%%rdx			\n\t"\
		"addq	%%rdx	,%%rsi			\n\t"\
		"\n\t"\
		"movq	%[__half_arr]	,%%rax	\n\t"\
		"movq	%[__wtA]	,%%rbx		\n\t"\
		"movq	%[__wtB]	,%%rcx		\n\t"\
		"movaps		     (%%rbx),%%xmm2	\n\t	movaps		 0x10(%%rbx),%%xmm6	\n\t"\
		"movhpd		     (%%rcx),%%xmm3	\n\t	movhpd		-0x10(%%rcx),%%xmm7	\n\t"\
		"movlpd		 0x08(%%rcx),%%xmm3	\n\t	movlpd		-0x08(%%rcx),%%xmm7	\n\t"\
		"\n\t"\
		"subq	$0x20	,%%rcx			\n\t"\
		"movq	%%rcx	,%[__wtB]		\n\t"\
		"\n\t"\
		"movq	%%rsi,	%%rdi	\n\t	movq	%%rsi,	%%rbx	\n\t"\
		"shrq	$12,	%%rdi	\n\t	shrq	$14,	%%rbx	\n\t"\
		"andq	$0x0000000000000030	,%%rdi		\n\t	andq	$0x0000000000000030	,%%rbx		\n\t"\
		"movq	%%rsi,	%%rdx	\n\t	movq	%%rsi,	%%rcx	\n\t"\
		"shrq	$4 ,	%%rdx	\n\t	shrq	$6 ,	%%rcx	\n\t"\
		"andq	$0x0000000000000030	,%%rdx		\n\t	andq	$0x0000000000000030	,%%rcx		\n\t"\
		"addq	%%rax	,%%rdi			\n\t	addq	%%rax	,%%rbx			\n\t"\
		"addq	%%rax	,%%rdx			\n\t	addq	%%rax	,%%rcx			\n\t"\
		"\n\t"\
		"mulpd	 	0x100(%%rax),%%xmm2	\n\t	mulpd	 	0x100(%%rax),%%xmm6	\n\t"\
		"mulpd	 	0x110(%%rax),%%xmm3	\n\t	mulpd	 	0x110(%%rax),%%xmm7	\n\t"\
		"mulpd	 	     (%%rdi),%%xmm2	\n\t	mulpd	 	     (%%rbx),%%xmm6	\n\t"\
		"mulpd	 	0x040(%%rdx),%%xmm3	\n\t	mulpd	 	0x040(%%rcx),%%xmm7	\n\t"\
		"\n\t"\
		"movq	%[__cyA]	,%%rcx		\n\t	movq	%[__cyB]	,%%rdx		\n\t"\
		"mulpd		%%xmm3	,%%xmm1		\n\t	mulpd		%%xmm7	,%%xmm5		\n\t"\
		"addpd		(%%rcx)	,%%xmm1		\n\t	addpd		(%%rdx)	,%%xmm5		\n\t"\
		"movaps		%%xmm1	,%%xmm3		\n\t	movaps		%%xmm5	,%%xmm7		\n\t"\
		"roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd		$0,%%xmm7,%%xmm7	\n\t"\
		"\n\t"\
		"/* NO ROE HERE */				\n\t"\
		"\n\t"\
		"movq	%%rsi,	%%rdi	\n\t	movq	%%rsi,	%%rbx	\n\t"\
		"shrq	$20,	%%rdi	\n\t	shrq	$22,	%%rbx	\n\t"\
		"andq	$0x0000000000000030	,%%rdi		\n\t	andq	$0x0000000000000030	,%%rbx		\n\t"\
		"addq	%%rax	,%%rdi			\n\t	addq	%%rax	,%%rbx			\n\t"\
		"movaps		%%xmm3	,%%xmm1		\n\t	movaps		%%xmm7	,%%xmm5		\n\t"\
		"mulpd	 	0xc0(%%rdi)	,%%xmm3	\n\t	mulpd	 	0xc0(%%rbx)	,%%xmm7	\n\t"\
		"roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd		$0,%%xmm7,%%xmm7	\n\t"\
		"movaps		%%xmm3	,(%%rcx)	\n\t	movaps		%%xmm7	,(%%rdx)	\n\t"\
		"\n\t"\
		"\n\t"\
		"movq	%[__data]	,%%rax		\n\t"\
		"mulpd	 0x80(%%rdi)	,%%xmm3	\n\t	mulpd	 0x80(%%rbx)	,%%xmm7	\n\t"\
		"subpd		%%xmm3	,	%%xmm1	\n\t	subpd		%%xmm7	,	%%xmm5	\n\t"\
		"mulpd		%%xmm2	,	%%xmm1	\n\t	mulpd		%%xmm6	,	%%xmm5	\n\t"\
		"movaps		%%xmm1	,    (%%rax)\n\t	movaps		%%xmm5	,0x40(%%rax)\n\t"\
		"\n\t"\
		"movq	%[__sse_n]	,%%rbx		\n\t"\
		"movaps		(%%rbx)	,%%xmm2		\n\t"\
		"movq	%[__sse_bw]	,%%rax		\n\t"\
		"paddd		(%%rax)	,%%xmm0		\n\t"\
		"movaps		%%xmm0	,%%xmm1		\n\t"\
		"pcmpgtd	%%xmm2	,%%xmm1		\n\t"\
		"pand		(%%rbx)	,%%xmm1		\n\t"\
		"psubd		%%xmm1	,%%xmm0		\n\t"\
	"/**********************************************/\n\t"\
	"/*          Imaginary parts               */\n\t"\
	"/**********************************************/\n\t"\
		"movq	%[__sse_sw]	,%%rdx		\n\t"\
		"movaps	(%%rdx)	,%%xmm1			\n\t"\
		"psubd	%%xmm0	,%%xmm1			\n\t"\
		"movmskps	%%xmm1	,%%rsi		\n\t"\
		"shlq	$24	,%%rsi				\n\t"\
		"movaps	%%xmm0	,%%xmm1			\n\t"\
		"movslq	%[__n_minus_silp1],%%rcx\n\t"\
		"movd	%%rcx	,%%xmm2			\n\t"\
		"pshufd	$0	,%%xmm2	,%%xmm2		\n\t"\
		"psubd	%%xmm0	,%%xmm2			\n\t"\
		"movmskps	%%xmm2	,%%rcx		\n\t"\
		"shlq	$16	,%%rcx				\n\t"\
		"addq	%%rcx	,%%rsi			\n\t"\
		"movslq	%[__sinwtm1]	,%%rdx	\n\t"\
		"movd	%%rdx	,%%xmm3			\n\t"\
		"pshufd	$0	,%%xmm3	,%%xmm3		\n\t"\
		"psubd	%%xmm3	,%%xmm1			\n\t"\
		"movmskps	%%xmm1	,%%rdx		\n\t"\
		"shlq	$8	,%%rdx				\n\t"\
		"addq	%%rdx	,%%rsi			\n\t"\
		"movq	%[__data]	,%%rax		\n\t"\
		"movaps	 0x10(%%rax)	,%%xmm1	\n\t"\
		"movaps	 0x50(%%rax)	,%%xmm5	\n\t"\
		"\n\t"\
		"movq	%[__half_arr]	,%%rax	\n\t"\
		"movq	%[__wtA]	,%%rbx		\n\t"\
		"movq	%[__wtC]	,%%rcx		\n\t"\
		"\n\t"\
		"movaps	     (%%rbx)	,%%xmm2	\n\t	movaps	 0x10(%%rbx)	,%%xmm6	\n\t"\
		"movhpd	     (%%rcx)	,%%xmm3	\n\t	movhpd	-0x10(%%rcx)	,%%xmm7	\n\t"\
		"movlpd	 0x08(%%rcx)	,%%xmm3	\n\t	movlpd	-0x08(%%rcx)	,%%xmm7	\n\t"\
		"\n\t"\
		"addq	$0x20	,%%rbx			\n\t"\
		"subq	$0x20	,%%rcx			\n\t"\
		"movq	%%rbx	,%[__wtA]		\n\t"\
		"movq	%%rcx	,%[__wtC]		\n\t"\
		"movq	%%rsi,	%%rdi	\n\t	movq	%%rsi,	%%rbx	\n\t"\
		"shrq	$12,	%%rdi	\n\t	shrq	$14,	%%rbx	\n\t"\
		"andq	$0x0000000000000030	,%%rdi		\n\t	andq	$0x0000000000000030	,%%rbx		\n\t"\
		"movq	%%rsi,	%%rdx	\n\t	movq	%%rsi,	%%rcx	\n\t"\
		"shrq	$4 ,	%%rdx	\n\t	shrq	$6 ,	%%rcx	\n\t"\
		"andq	$0x0000000000000030	,%%rdx		\n\t	andq	$0x0000000000000030	,%%rcx		\n\t"\
		"addq	%%rax	,%%rdi			\n\t	addq	%%rax	,%%rbx			\n\t"\
		"addq	%%rax	,%%rdx			\n\t	addq	%%rax	,%%rcx			\n\t"\
		"\n\t"\
		"mulpd	 0x120(%%rax)	,%%xmm2	\n\t	mulpd	 0x120(%%rax)	,%%xmm6	\n\t"\
		"mulpd	 0x130(%%rax)	,%%xmm3	\n\t	mulpd	 0x130(%%rax)	,%%xmm7	\n\t"\
		"mulpd	      (%%rdi)	,%%xmm2	\n\t	mulpd	      (%%rbx)	,%%xmm6	\n\t"\
		"mulpd	 0x040(%%rdx)	,%%xmm3	\n\t	mulpd	 0x040(%%rcx)	,%%xmm7	\n\t"\
		"\n\t"\
		"movq	%[__cyA]	,%%rcx		\n\t	movq	%[__cyB]	,%%rdx		\n\t"\
		"mulpd	%%xmm3	,%%xmm1			\n\t	mulpd	%%xmm7	,%%xmm5			\n\t"\
		"addpd	(%%rcx)	,%%xmm1			\n\t	addpd	(%%rdx)	,%%xmm5			\n\t"\
		"movaps	%%xmm1	,%%xmm3			\n\t	movaps	%%xmm5	,%%xmm7			\n\t"\
		"roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd		$0,%%xmm7,%%xmm7	\n\t"\
		"\n\t"\
		"/* NO ROE HERE */				\n\t"\
		"\n\t"\
		"movq	%%rsi,	%%rdi	\n\t	movq	%%rsi,	%%rbx	\n\t"\
		"shrq	$20,	%%rdi	\n\t	shrq	$22,	%%rbx	\n\t"\
		"andq	$0x0000000000000030	,%%rdi		\n\t	andq	$0x0000000000000030	,%%rbx		\n\t"\
		"addq	%%rax	,%%rdi			\n\t	addq	%%rax	,%%rbx			\n\t"\
		"movaps	%%xmm3	,%%xmm1			\n\t	movaps	%%xmm7	,%%xmm5			\n\t"\
		"mulpd	 0xc0(%%rdi)	,%%xmm3	\n\t	mulpd	 0xc0(%%rbx)	,%%xmm7	\n\t"\
		"roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd		$0,%%xmm7,%%xmm7	\n\t"\
		"movaps	%%xmm3	,(%%rcx)		\n\t	movaps	%%xmm7	,(%%rdx)		\n\t"\
		"\n\t"\
		"movq	%[__data]	,%%rax		\n\t"\
		"mulpd	 0x80(%%rdi)	,%%xmm3	\n\t	mulpd	 0x80(%%rbx)	,%%xmm7	\n\t"\
		"subpd	%%xmm3	,%%xmm1			\n\t	subpd	%%xmm7	,%%xmm5			\n\t"\
		"mulpd	%%xmm2	,%%xmm1			\n\t	mulpd	%%xmm6	,%%xmm5			\n\t"\
		"movaps	%%xmm1	, 0x10(%%rax)	\n\t	movaps	%%xmm5	, 0x50(%%rax)	\n\t"\
		"\n\t"\
		"movq	%[__sse_n]	,%%rbx		\n\t"\
		"movaps		(%%rbx)	,%%xmm2		\n\t"\
		"movq	%[__sse_bw]	,%%rax		\n\t"\
		"paddd		(%%rax)	,%%xmm0		\n\t"\
		"movaps		%%xmm0	,%%xmm1		\n\t"\
		"pcmpgtd	%%xmm2	,%%xmm1		\n\t"\
		"pand		(%%rbx)	,%%xmm1		\n\t"\
		"psubd		%%xmm1	,%%xmm0		\n\t"\
		"movq	%[__bjmod_0],%%rcx		\n\t"\
		"movaps	%%xmm0,(%%rcx)			\n\t"\
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
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
	);\
	}

	#define SSE2_cmplx_carry_norm_nocheck2_2B(Xdata,XwtA,XwtB,XcyA,XcyB,Xbjmod_0,Xhalf_arr,Xn_minus_silp1,Xn_minus_sil,Xsinwt,Xsinwtm1,Xsse_bw,Xsse_n,Xsse_sw)\
	{\
	__asm__ volatile (\
	"/**********************************************/\n\t"\
	"/*          Real      parts                   */\n\t"\
	"/**********************************************/\n\t"\
		"movq	%[__bjmod_0],	%%rax	\n\t"\
		"movaps		(%%rax)	,	%%xmm0	\n\t"\
		"movq	%[__sse_sw]	,	%%rbx	\n\t"\
		"movaps		(%%rbx)	,	%%xmm1	\n\t"\
		"psubd		%%xmm0	,	%%xmm1	\n\t"\
		"movmskps	%%xmm1	,	%%rsi	\n\t"\
		"\n\t"\
		"shlq	$24		,%%rsi			\n\t"\
		"movaps		%%xmm0	,%%xmm1		\n\t"\
		"movslq	%[__n_minus_sil],%%rcx	\n\t"\
		"movd	%%rcx	,%%xmm2			\n\t"\
		"pshufd	$0,	%%xmm2	,%%xmm2		\n\t"\
		"psubd		%%xmm0	,%%xmm2		\n\t"\
		"movmskps	%%xmm2	,%%rcx		\n\t"\
		"shlq	$16		,%%rcx			\n\t"\
		"addq	%%rcx	,%%rsi			\n\t"\
		"movslq	%[__sinwt]	,%%rdx		\n\t"\
		"movd	%%rdx	,%%xmm3			\n\t"\
		"pshufd	$0,	%%xmm3	,%%xmm3		\n\t"\
		"psubd		%%xmm3	,%%xmm1		\n\t"\
		"movmskps	%%xmm1	,%%rdx		\n\t"\
		"shlq	$8		,%%rdx			\n\t"\
		"addq	%%rdx	,%%rsi			\n\t"\
		"\n\t"\
		"mov	%[__data],%%rax			\n\t"\
		"movaps	 0x20(%%rax),%%xmm1		\n\t	movaps		 0x60(%%rax),%%xmm5\n\t"\
		"\n\t"\
		"movq	%[__half_arr]	,%%rax	\n\t"\
		"movq	%[__wtA]	,%%rbx		\n\t"\
		"movq	%[__wtB]	,%%rcx		\n\t"\
		"movaps		     (%%rbx),%%xmm2	\n\t	movaps		 0x10(%%rbx),%%xmm6	\n\t"\
		"movhpd		     (%%rcx),%%xmm3	\n\t	movhpd		-0x10(%%rcx),%%xmm7	\n\t"\
		"movlpd		 0x08(%%rcx),%%xmm3	\n\t	movlpd		-0x08(%%rcx),%%xmm7	\n\t"\
		"\n\t"\
		"movq	%%rsi,	%%rdi	\n\t	movq	%%rsi,	%%rbx	\n\t"\
		"shrq	$12,	%%rdi	\n\t	shrq	$14,	%%rbx	\n\t"\
		"andq	$0x0000000000000030	,%%rdi		\n\t	andq	$0x0000000000000030	,%%rbx		\n\t"\
		"movq	%%rsi,	%%rdx	\n\t	movq	%%rsi,	%%rcx	\n\t"\
		"shrq	$4 ,	%%rdx	\n\t	shrq	$6 ,	%%rcx	\n\t"\
		"andq	$0x0000000000000030	,%%rdx		\n\t	andq	$0x0000000000000030	,%%rcx		\n\t"\
		"addq	%%rax	,%%rdi			\n\t	addq	%%rax	,%%rbx			\n\t"\
		"addq	%%rax	,%%rdx			\n\t	addq	%%rax	,%%rcx			\n\t"\
		"\n\t"\
		"mulpd	 	0x100(%%rax),%%xmm2	\n\t	mulpd	 	0x100(%%rax),%%xmm6	\n\t"\
		"mulpd	 	0x110(%%rax),%%xmm3	\n\t	mulpd	 	0x110(%%rax),%%xmm7	\n\t"\
		"mulpd	 	     (%%rdi),%%xmm2	\n\t	mulpd	 	     (%%rbx),%%xmm6	\n\t"\
		"mulpd	 	0x040(%%rdx),%%xmm3	\n\t	mulpd	 	0x040(%%rcx),%%xmm7	\n\t"\
		"\n\t"\
		"movq	%[__cyA]	,%%rcx		\n\t	movq	%[__cyB]	,%%rdx		\n\t"\
		"mulpd		%%xmm3	,%%xmm1		\n\t	mulpd		%%xmm7	,%%xmm5		\n\t"\
		"addpd		(%%rcx)	,%%xmm1		\n\t	addpd		(%%rdx)	,%%xmm5		\n\t"\
		"movaps		%%xmm1	,%%xmm3		\n\t	movaps		%%xmm5	,%%xmm7		\n\t"\
		"roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd		$0,%%xmm7,%%xmm7	\n\t"\
		"\n\t"\
		"/* NO ROE HERE */				\n\t"\
		"\n\t"\
		"movq	%%rsi,	%%rdi	\n\t	movq	%%rsi,	%%rbx	\n\t"\
		"shrq	$20,	%%rdi	\n\t	shrq	$22,	%%rbx	\n\t"\
		"andq	$0x0000000000000030	,%%rdi		\n\t	andq	$0x0000000000000030	,%%rbx		\n\t"\
		"addq	%%rax	,%%rdi			\n\t	addq	%%rax	,%%rbx			\n\t"\
		"movaps		%%xmm3	,%%xmm1		\n\t	movaps		%%xmm7	,%%xmm5		\n\t"\
		"mulpd	 	0xc0(%%rdi)	,%%xmm3	\n\t	mulpd	 	0xc0(%%rbx)	,%%xmm7	\n\t"\
		"roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd		$0,%%xmm7,%%xmm7	\n\t"\
		"movaps		%%xmm3	,(%%rcx)	\n\t	movaps		%%xmm7	,(%%rdx)	\n\t"\
		"\n\t"\
		"\n\t"\
		"movq	%[__data]	,%%rax		\n\t"\
		"mulpd	 0x80(%%rdi)	,%%xmm3	\n\t	mulpd	 0x80(%%rbx)	,%%xmm7	\n\t"\
		"subpd		%%xmm3	,	%%xmm1	\n\t	subpd		%%xmm7	,	%%xmm5	\n\t"\
		"mulpd		%%xmm2	,	%%xmm1	\n\t	mulpd		%%xmm6	,	%%xmm5	\n\t"\
		"movaps		%%xmm1	,0x20(%%rax)\n\t	movaps		%%xmm5	,0x60(%%rax)\n\t"\
		"\n\t"\
		"movq	%[__sse_n]	,%%rbx		\n\t"\
		"movaps		(%%rbx)	,%%xmm2		\n\t"\
		"movq	%[__sse_bw]	,%%rax		\n\t"\
		"paddd		(%%rax)	,%%xmm0		\n\t"\
		"movaps		%%xmm0	,%%xmm1		\n\t"\
		"pcmpgtd	%%xmm2	,%%xmm1		\n\t"\
		"pand		(%%rbx)	,%%xmm1		\n\t"\
		"psubd		%%xmm1	,%%xmm0		\n\t"\
	"/**********************************************/\n\t"\
	"/*          Imaginary parts               */\n\t"\
	"/**********************************************/\n\t"\
		"movq	%[__sse_sw]	,%%rdx		\n\t"\
		"movaps	(%%rdx)	,%%xmm1			\n\t"\
		"psubd	%%xmm0	,%%xmm1			\n\t"\
		"movmskps	%%xmm1	,%%rsi		\n\t"\
		"shlq	$24	,%%rsi				\n\t"\
		"movaps	%%xmm0	,%%xmm1			\n\t"\
		"movslq	%[__n_minus_silp1],%%rcx\n\t"\
		"movd	%%rcx	,%%xmm2			\n\t"\
		"pshufd	$0	,%%xmm2	,%%xmm2		\n\t"\
		"psubd	%%xmm0	,%%xmm2			\n\t"\
		"movmskps	%%xmm2	,%%rcx		\n\t"\
		"shlq	$16	,%%rcx				\n\t"\
		"addq	%%rcx	,%%rsi			\n\t"\
		"movslq	%[__sinwtm1]	,%%rdx	\n\t"\
		"movd	%%rdx	,%%xmm3			\n\t"\
		"pshufd	$0	,%%xmm3	,%%xmm3		\n\t"\
		"psubd	%%xmm3	,%%xmm1			\n\t"\
		"movmskps	%%xmm1	,%%rdx		\n\t"\
		"shlq	$8	,%%rdx				\n\t"\
		"addq	%%rdx	,%%rsi			\n\t"\
		"movq	%[__data]	,%%rax		\n\t"\
		"movaps	 0x30(%%rax)	,%%xmm1	\n\t"\
		"movaps	 0x70(%%rax)	,%%xmm5	\n\t"\
		"\n\t"\
		"movq	%[__half_arr]	,%%rax	\n\t"\
		"movq	%[__wtA]	,%%rbx		\n\t"\
		"movq	%[__wtB]	,%%rcx		\n\t"\
		"\n\t"\
		"movaps	     (%%rbx)	,%%xmm2	\n\t	movaps	 0x10(%%rbx)	,%%xmm6	\n\t"\
		"movhpd	     (%%rcx)	,%%xmm3	\n\t	movhpd	-0x10(%%rcx)	,%%xmm7	\n\t"\
		"movlpd	 0x08(%%rcx)	,%%xmm3	\n\t	movlpd	-0x08(%%rcx)	,%%xmm7	\n\t"\
		"\n\t"\
		"addq	$0x20	,%%rbx			\n\t"\
		"subq	$0x20	,%%rcx			\n\t"\
		"movq	%%rbx	,%[__wtA]		\n\t"\
		"movq	%%rcx	,%[__wtB]		\n\t"\
		"movq	%%rsi,	%%rdi	\n\t	movq	%%rsi,	%%rbx	\n\t"\
		"shrq	$12,	%%rdi	\n\t	shrq	$14,	%%rbx	\n\t"\
		"andq	$0x0000000000000030	,%%rdi		\n\t	andq	$0x0000000000000030	,%%rbx		\n\t"\
		"movq	%%rsi,	%%rdx	\n\t	movq	%%rsi,	%%rcx	\n\t"\
		"shrq	$4 ,	%%rdx	\n\t	shrq	$6 ,	%%rcx	\n\t"\
		"andq	$0x0000000000000030	,%%rdx		\n\t	andq	$0x0000000000000030	,%%rcx		\n\t"\
		"addq	%%rax	,%%rdi			\n\t	addq	%%rax	,%%rbx			\n\t"\
		"addq	%%rax	,%%rdx			\n\t	addq	%%rax	,%%rcx			\n\t"\
		"\n\t"\
		"mulpd	 0x120(%%rax)	,%%xmm2	\n\t	mulpd	 0x120(%%rax)	,%%xmm6	\n\t"\
		"mulpd	 0x130(%%rax)	,%%xmm3	\n\t	mulpd	 0x130(%%rax)	,%%xmm7	\n\t"\
		"mulpd	      (%%rdi)	,%%xmm2	\n\t	mulpd	      (%%rbx)	,%%xmm6	\n\t"\
		"mulpd	 0x040(%%rdx)	,%%xmm3	\n\t	mulpd	 0x040(%%rcx)	,%%xmm7	\n\t"\
		"\n\t"\
		"movq	%[__cyA]	,%%rcx		\n\t	movq	%[__cyB]	,%%rdx		\n\t"\
		"mulpd	%%xmm3	,%%xmm1			\n\t	mulpd	%%xmm7	,%%xmm5			\n\t"\
		"addpd	(%%rcx)	,%%xmm1			\n\t	addpd	(%%rdx)	,%%xmm5			\n\t"\
		"movaps	%%xmm1	,%%xmm3			\n\t	movaps	%%xmm5	,%%xmm7			\n\t"\
		"roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd		$0,%%xmm7,%%xmm7	\n\t"\
		"\n\t"\
		"/* NO ROE HERE */				\n\t"\
		"\n\t"\
		"movq	%%rsi,	%%rdi	\n\t	movq	%%rsi,	%%rbx	\n\t"\
		"shrq	$20,	%%rdi	\n\t	shrq	$22,	%%rbx	\n\t"\
		"andq	$0x0000000000000030	,%%rdi		\n\t	andq	$0x0000000000000030	,%%rbx		\n\t"\
		"addq	%%rax	,%%rdi			\n\t	addq	%%rax	,%%rbx			\n\t"\
		"movaps	%%xmm3	,%%xmm1			\n\t	movaps	%%xmm7	,%%xmm5			\n\t"\
		"mulpd	 0xc0(%%rdi)	,%%xmm3	\n\t	mulpd	 0xc0(%%rbx)	,%%xmm7	\n\t"\
		"roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd		$0,%%xmm7,%%xmm7	\n\t"\
		"movaps	%%xmm3	,(%%rcx)		\n\t	movaps	%%xmm7	,(%%rdx)		\n\t"\
		"\n\t"\
		"\n\t"\
		"movq	%[__data]	,%%rax		\n\t"\
		"mulpd	 0x80(%%rdi)	,%%xmm3	\n\t	mulpd	 0x80(%%rbx)	,%%xmm7	\n\t"\
		"subpd	%%xmm3	,%%xmm1			\n\t	subpd	%%xmm7	,%%xmm5			\n\t"\
		"mulpd	%%xmm2	,%%xmm1			\n\t	mulpd	%%xmm6	,%%xmm5			\n\t"\
		"movaps	%%xmm1	, 0x30(%%rax)	\n\t	movaps	%%xmm5	, 0x70(%%rax)	\n\t"\
		"\n\t"\
		"movq	%[__sse_n]	,%%rbx		\n\t"\
		"movaps		(%%rbx)	,%%xmm2		\n\t"\
		"movq	%[__sse_bw]	,%%rax		\n\t"\
		"paddd		(%%rax)	,%%xmm0		\n\t"\
		"movaps		%%xmm0	,%%xmm1		\n\t"\
		"pcmpgtd	%%xmm2	,%%xmm1		\n\t"\
		"pand		(%%rbx)	,%%xmm1		\n\t"\
		"psubd		%%xmm1	,%%xmm0		\n\t"\
		"movq	%[__bjmod_0],%%rcx		\n\t"\
		"movaps	%%xmm0,(%%rcx)			\n\t"\
	"/**********************************************/\n\t"\
	"/*          Repack the data:              */\n\t"\
	"/**********************************************/\n\t"\
		"mov	%[__data],%%rax			\n\t"\
		"movaps		0x10(%%rax)	,%%xmm1	\n\t	movaps		0x50(%%rax)	,%%xmm5	\n\t"\
		"movaps		    (%%rax)	,%%xmm0	\n\t	movaps		0x40(%%rax)	,%%xmm4	\n\t"\
		"movaps		%%xmm1		,%%xmm3	\n\t	movaps		%%xmm5		,%%xmm7	\n\t"\
		"movaps		%%xmm0		,%%xmm2	\n\t	movaps		%%xmm4		,%%xmm6	\n\t"\
		"unpckhpd	0x30(%%rax)	,%%xmm3	\n\t	unpckhpd	0x70(%%rax)	,%%xmm7	\n\t"\
		"unpcklpd	0x30(%%rax)	,%%xmm1	\n\t	unpcklpd	0x70(%%rax)	,%%xmm5	\n\t"\
		"movaps		%%xmm3,0x30(%%rax)	\n\t	movaps		%%xmm7,0x70(%%rax)	\n\t"\
		"unpckhpd	0x20(%%rax)	,%%xmm2	\n\t	unpckhpd	0x60(%%rax)	,%%xmm6	\n\t"\
		"unpcklpd	0x20(%%rax)	,%%xmm0	\n\t	unpcklpd	0x60(%%rax)	,%%xmm4	\n\t"\
		"movaps		%%xmm2,0x20(%%rax)	\n\t	movaps		%%xmm6,0x60(%%rax)	\n\t"\
		"movaps		%%xmm1,0x10(%%rax)	\n\t	movaps		%%xmm5,0x50(%%rax)	\n\t"\
		"movaps		%%xmm0,    (%%rax)	\n\t	movaps		%%xmm4,0x40(%%rax)	\n\t"\
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
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
	);\
	}

#endif	/* carry_gcc_h_included */

