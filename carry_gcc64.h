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
#ifndef carry_gcc_h_included
#define carry_gcc_h_included

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
		"movaps		-0x10(%%rax),%%xmm4	\n\t"\
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
		"addpd		%%xmm4	,%%xmm3		\n\t	addpd		%%xmm4	,%%xmm7		\n\t"\
		"subpd		%%xmm4	,%%xmm3		\n\t	subpd		%%xmm4	,%%xmm7		\n\t"\
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
		"addpd		%%xmm4	,%%xmm3		\n\t	addpd		%%xmm4	,%%xmm7		\n\t"\
		"subpd		%%xmm4	,%%xmm3		\n\t	subpd		%%xmm4	,%%xmm7		\n\t"\
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
		"addpd	%%xmm4	,%%xmm3			\n\t	addpd	%%xmm4	,%%xmm7			\n\t"\
		"subpd	%%xmm4	,%%xmm3			\n\t	subpd	%%xmm4	,%%xmm7			\n\t"\
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
		"addpd	%%xmm4	,%%xmm3			\n\t	addpd	%%xmm4	,%%xmm7			\n\t"\
		"subpd	%%xmm4	,%%xmm3			\n\t	subpd	%%xmm4	,%%xmm7			\n\t"\
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
		: "rax","rbx","rcx","rdx","rdi","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
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
		"movaps		-0x10(%%rax),%%xmm4	\n\t"\
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
		"addpd		%%xmm4	,%%xmm3		\n\t	addpd		%%xmm4	,%%xmm7		\n\t"\
		"subpd		%%xmm4	,%%xmm3		\n\t	subpd		%%xmm4	,%%xmm7		\n\t"\
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
		"addpd		%%xmm4	,%%xmm3		\n\t	addpd		%%xmm4	,%%xmm7		\n\t"\
		"subpd		%%xmm4	,%%xmm3		\n\t	subpd		%%xmm4	,%%xmm7		\n\t"\
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
		"addpd	%%xmm4	,%%xmm3			\n\t	addpd	%%xmm4	,%%xmm7			\n\t"\
		"subpd	%%xmm4	,%%xmm3			\n\t	subpd	%%xmm4	,%%xmm7			\n\t"\
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
		"addpd	%%xmm4	,%%xmm3			\n\t	addpd	%%xmm4	,%%xmm7			\n\t"\
		"subpd	%%xmm4	,%%xmm3			\n\t	subpd	%%xmm4	,%%xmm7			\n\t"\
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
		: "rax","rbx","rcx","rdx","rdi","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
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
		"movaps		-0x10(%%rax),%%xmm4	\n\t"\
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
		"addpd		%%xmm4	,%%xmm3		\n\t	addpd		%%xmm4	,%%xmm7		\n\t"\
		"subpd		%%xmm4	,%%xmm3		\n\t	subpd		%%xmm4	,%%xmm7		\n\t"\
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
		"addpd		%%xmm4	,%%xmm3		\n\t	addpd		%%xmm4	,%%xmm7		\n\t"\
		"subpd		%%xmm4	,%%xmm3		\n\t	subpd		%%xmm4	,%%xmm7		\n\t"\
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
		"addpd	%%xmm4	,%%xmm3			\n\t	addpd	%%xmm4	,%%xmm7			\n\t"\
		"subpd	%%xmm4	,%%xmm3			\n\t	subpd	%%xmm4	,%%xmm7			\n\t"\
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
		"addpd	%%xmm4	,%%xmm3			\n\t	addpd	%%xmm4	,%%xmm7			\n\t"\
		"subpd	%%xmm4	,%%xmm3			\n\t	subpd	%%xmm4	,%%xmm7			\n\t"\
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
		: "rax","rbx","rcx","rdx","rdi","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
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
		"movaps		-0x10(%%rax),%%xmm4	\n\t"\
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
		"addpd		%%xmm4	,%%xmm3		\n\t	addpd		%%xmm4	,%%xmm7		\n\t"\
		"subpd		%%xmm4	,%%xmm3		\n\t	subpd		%%xmm4	,%%xmm7		\n\t"\
		"\n\t"\
		"/* NO ROE HERE */		\n\t"\
		"\n\t"\
		"movq	%%rsi,	%%rdi	\n\t	movq	%%rsi,	%%rbx	\n\t"\
		"shrq	$20,	%%rdi	\n\t	shrq	$22,	%%rbx	\n\t"\
		"andq	$0x0000000000000030	,%%rdi		\n\t	andq	$0x0000000000000030	,%%rbx		\n\t"\
		"addq	%%rax	,%%rdi			\n\t	addq	%%rax	,%%rbx			\n\t"\
		"movaps		%%xmm3	,%%xmm1		\n\t	movaps		%%xmm7	,%%xmm5		\n\t"\
		"mulpd	 	0xc0(%%rdi)	,%%xmm3	\n\t	mulpd	 	0xc0(%%rbx)	,%%xmm7	\n\t"\
		"addpd		%%xmm4	,%%xmm3		\n\t	addpd		%%xmm4	,%%xmm7		\n\t"\
		"subpd		%%xmm4	,%%xmm3		\n\t	subpd		%%xmm4	,%%xmm7		\n\t"\
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
		"addpd	%%xmm4	,%%xmm3			\n\t	addpd	%%xmm4	,%%xmm7			\n\t"\
		"subpd	%%xmm4	,%%xmm3			\n\t	subpd	%%xmm4	,%%xmm7			\n\t"\
		"\n\t"\
		"/* NO ROE HERE */		\n\t"\
		"\n\t"\
		"movq	%%rsi,	%%rdi	\n\t	movq	%%rsi,	%%rbx	\n\t"\
		"shrq	$20,	%%rdi	\n\t	shrq	$22,	%%rbx	\n\t"\
		"andq	$0x0000000000000030	,%%rdi		\n\t	andq	$0x0000000000000030	,%%rbx		\n\t"\
		"addq	%%rax	,%%rdi			\n\t	addq	%%rax	,%%rbx			\n\t"\
		"movaps	%%xmm3	,%%xmm1			\n\t	movaps	%%xmm7	,%%xmm5			\n\t"\
		"mulpd	 0xc0(%%rdi)	,%%xmm3	\n\t	mulpd	 0xc0(%%rbx)	,%%xmm7	\n\t"\
		"addpd	%%xmm4	,%%xmm3			\n\t	addpd	%%xmm4	,%%xmm7			\n\t"\
		"subpd	%%xmm4	,%%xmm3			\n\t	subpd	%%xmm4	,%%xmm7			\n\t"\
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
		: "rax","rbx","rcx","rdx","rdi","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
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
		"movaps		-0x10(%%rax),%%xmm4	\n\t"\
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
		"addpd		%%xmm4	,%%xmm3		\n\t	addpd		%%xmm4	,%%xmm7		\n\t"\
		"subpd		%%xmm4	,%%xmm3		\n\t	subpd		%%xmm4	,%%xmm7		\n\t"\
		"\n\t"\
		"/* NO ROE HERE */		\n\t"\
		"\n\t"\
		"movq	%%rsi,	%%rdi	\n\t	movq	%%rsi,	%%rbx	\n\t"\
		"shrq	$20,	%%rdi	\n\t	shrq	$22,	%%rbx	\n\t"\
		"andq	$0x0000000000000030	,%%rdi		\n\t	andq	$0x0000000000000030	,%%rbx		\n\t"\
		"addq	%%rax	,%%rdi			\n\t	addq	%%rax	,%%rbx			\n\t"\
		"movaps		%%xmm3	,%%xmm1		\n\t	movaps		%%xmm7	,%%xmm5		\n\t"\
		"mulpd	 	0xc0(%%rdi)	,%%xmm3	\n\t	mulpd	 	0xc0(%%rbx)	,%%xmm7	\n\t"\
		"addpd		%%xmm4	,%%xmm3		\n\t	addpd		%%xmm4	,%%xmm7		\n\t"\
		"subpd		%%xmm4	,%%xmm3		\n\t	subpd		%%xmm4	,%%xmm7		\n\t"\
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
		"addpd	%%xmm4	,%%xmm3			\n\t	addpd	%%xmm4	,%%xmm7			\n\t"\
		"subpd	%%xmm4	,%%xmm3			\n\t	subpd	%%xmm4	,%%xmm7			\n\t"\
		"\n\t"\
		"/* NO ROE HERE */		\n\t"\
		"\n\t"\
		"movq	%%rsi,	%%rdi	\n\t	movq	%%rsi,	%%rbx	\n\t"\
		"shrq	$20,	%%rdi	\n\t	shrq	$22,	%%rbx	\n\t"\
		"andq	$0x0000000000000030	,%%rdi		\n\t	andq	$0x0000000000000030	,%%rbx		\n\t"\
		"addq	%%rax	,%%rdi			\n\t	addq	%%rax	,%%rbx			\n\t"\
		"movaps	%%xmm3	,%%xmm1			\n\t	movaps	%%xmm7	,%%xmm5			\n\t"\
		"mulpd	 0xc0(%%rdi)	,%%xmm3	\n\t	mulpd	 0xc0(%%rbx)	,%%xmm7	\n\t"\
		"addpd	%%xmm4	,%%xmm3			\n\t	addpd	%%xmm4	,%%xmm7			\n\t"\
		"subpd	%%xmm4	,%%xmm3			\n\t	subpd	%%xmm4	,%%xmm7			\n\t"\
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
		: "rax","rbx","rcx","rdx","rdi","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
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
		"movaps		-0x10(%%rax),%%xmm4	\n\t"\
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
		"addpd		%%xmm4	,%%xmm3		\n\t	addpd		%%xmm4	,%%xmm7		\n\t"\
		"subpd		%%xmm4	,%%xmm3		\n\t	subpd		%%xmm4	,%%xmm7		\n\t"\
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
		"addpd		%%xmm4	,%%xmm3		\n\t	addpd		%%xmm4	,%%xmm7		\n\t"\
		"subpd		%%xmm4	,%%xmm3		\n\t	subpd		%%xmm4	,%%xmm7		\n\t"\
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
		"addpd	%%xmm4	,%%xmm3			\n\t	addpd	%%xmm4	,%%xmm7			\n\t"\
		"subpd	%%xmm4	,%%xmm3			\n\t	subpd	%%xmm4	,%%xmm7			\n\t"\
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
		"addpd	%%xmm4	,%%xmm3			\n\t	addpd	%%xmm4	,%%xmm7			\n\t"\
		"subpd	%%xmm4	,%%xmm3			\n\t	subpd	%%xmm4	,%%xmm7			\n\t"\
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
		: "rax","rbx","rcx","rdx","rdi","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
	);\
	}

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
		"movaps		-0x10(%%rax),%%xmm4	\n\t"\
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
		"addpd		%%xmm4	,%%xmm3		\n\t	addpd		%%xmm4	,%%xmm7		\n\t"\
		"subpd		%%xmm4	,%%xmm3		\n\t	subpd		%%xmm4	,%%xmm7		\n\t"\
		"\n\t"\
		"/* NO ROE HERE */				\n\t"\
		"\n\t"\
		"movq	%%rsi,	%%rdi	\n\t	movq	%%rsi,	%%rbx	\n\t"\
		"shrq	$20,	%%rdi	\n\t	shrq	$22,	%%rbx	\n\t"\
		"andq	$0x0000000000000030	,%%rdi		\n\t	andq	$0x0000000000000030	,%%rbx		\n\t"\
		"addq	%%rax	,%%rdi			\n\t	addq	%%rax	,%%rbx			\n\t"\
		"movaps		%%xmm3	,%%xmm1		\n\t	movaps		%%xmm7	,%%xmm5		\n\t"\
		"mulpd	 	0xc0(%%rdi)	,%%xmm3	\n\t	mulpd	 	0xc0(%%rbx)	,%%xmm7	\n\t"\
		"addpd		%%xmm4	,%%xmm3		\n\t	addpd		%%xmm4	,%%xmm7		\n\t"\
		"subpd		%%xmm4	,%%xmm3		\n\t	subpd		%%xmm4	,%%xmm7		\n\t"\
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
		"addpd	%%xmm4	,%%xmm3			\n\t	addpd	%%xmm4	,%%xmm7			\n\t"\
		"subpd	%%xmm4	,%%xmm3			\n\t	subpd	%%xmm4	,%%xmm7			\n\t"\
		"\n\t"\
		"/* NO ROE HERE */				\n\t"\
		"\n\t"\
		"movq	%%rsi,	%%rdi	\n\t	movq	%%rsi,	%%rbx	\n\t"\
		"shrq	$20,	%%rdi	\n\t	shrq	$22,	%%rbx	\n\t"\
		"andq	$0x0000000000000030	,%%rdi		\n\t	andq	$0x0000000000000030	,%%rbx		\n\t"\
		"addq	%%rax	,%%rdi			\n\t	addq	%%rax	,%%rbx			\n\t"\
		"movaps	%%xmm3	,%%xmm1			\n\t	movaps	%%xmm7	,%%xmm5			\n\t"\
		"mulpd	 0xc0(%%rdi)	,%%xmm3	\n\t	mulpd	 0xc0(%%rbx)	,%%xmm7	\n\t"\
		"addpd	%%xmm4	,%%xmm3			\n\t	addpd	%%xmm4	,%%xmm7			\n\t"\
		"subpd	%%xmm4	,%%xmm3			\n\t	subpd	%%xmm4	,%%xmm7			\n\t"\
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
		: "rax","rbx","rcx","rdx","rdi","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
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
		"movaps		-0x10(%%rax),%%xmm4	\n\t"\
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
		"addpd		%%xmm4	,%%xmm3		\n\t	addpd		%%xmm4	,%%xmm7		\n\t"\
		"subpd		%%xmm4	,%%xmm3		\n\t	subpd		%%xmm4	,%%xmm7		\n\t"\
		"\n\t"\
		"/* NO ROE HERE */				\n\t"\
		"\n\t"\
		"movq	%%rsi,	%%rdi	\n\t	movq	%%rsi,	%%rbx	\n\t"\
		"shrq	$20,	%%rdi	\n\t	shrq	$22,	%%rbx	\n\t"\
		"andq	$0x0000000000000030	,%%rdi		\n\t	andq	$0x0000000000000030	,%%rbx		\n\t"\
		"addq	%%rax	,%%rdi			\n\t	addq	%%rax	,%%rbx			\n\t"\
		"movaps		%%xmm3	,%%xmm1		\n\t	movaps		%%xmm7	,%%xmm5		\n\t"\
		"mulpd	 	0xc0(%%rdi)	,%%xmm3	\n\t	mulpd	 	0xc0(%%rbx)	,%%xmm7	\n\t"\
		"addpd		%%xmm4	,%%xmm3		\n\t	addpd		%%xmm4	,%%xmm7		\n\t"\
		"subpd		%%xmm4	,%%xmm3		\n\t	subpd		%%xmm4	,%%xmm7		\n\t"\
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
		"addpd	%%xmm4	,%%xmm3			\n\t	addpd	%%xmm4	,%%xmm7			\n\t"\
		"subpd	%%xmm4	,%%xmm3			\n\t	subpd	%%xmm4	,%%xmm7			\n\t"\
		"\n\t"\
		"/* NO ROE HERE */				\n\t"\
		"\n\t"\
		"movq	%%rsi,	%%rdi	\n\t	movq	%%rsi,	%%rbx	\n\t"\
		"shrq	$20,	%%rdi	\n\t	shrq	$22,	%%rbx	\n\t"\
		"andq	$0x0000000000000030	,%%rdi		\n\t	andq	$0x0000000000000030	,%%rbx		\n\t"\
		"addq	%%rax	,%%rdi			\n\t	addq	%%rax	,%%rbx			\n\t"\
		"movaps	%%xmm3	,%%xmm1			\n\t	movaps	%%xmm7	,%%xmm5			\n\t"\
		"mulpd	 0xc0(%%rdi)	,%%xmm3	\n\t	mulpd	 0xc0(%%rbx)	,%%xmm7	\n\t"\
		"addpd	%%xmm4	,%%xmm3			\n\t	addpd	%%xmm4	,%%xmm7			\n\t"\
		"subpd	%%xmm4	,%%xmm3			\n\t	subpd	%%xmm4	,%%xmm7			\n\t"\
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
		: "rax","rbx","rcx","rdx","rdi","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
	);\
	}

#endif	/* carry_gcc_h_included */

