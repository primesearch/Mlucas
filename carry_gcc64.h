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
#if 0
COMPILER BUG NOTE: If encounter errors of the following kinds:
clang:
	fatal error: error in backend: Ran out of registers during register allocation!
	Please check your inline asm statement for invalid constraints:
gcc:
	error: cannot find a register in class ‘GENERAL_REGS’ while reloading ‘asm’
	error: ‘asm’ operand has impossible constraints

Check the compile optimization level - If 0, try upping to at east -O1.
#endif
/*******************************************************************************
   We now include this header file if it was not included before.
*******************************************************************************/
#ifndef carry_gcc_h_included
#define carry_gcc_h_included

	/************** See the Visual-studio-style 32-bit analogs of these in carry.h for commented versions: **********/

	/*************************************************************/
	/**************** FERMAT  -MOD CARRY MACROS ******************/
	/*************************************************************/

#ifdef USE_AVX

/*
Use MOVUPD (or 1-byte-shorter MOVUPS) in legacy 128-bit SSE form to load 2 doubles into lo128 without touching hi128;
Thus can do MOVUPD m128,xmm1 to fill lo128 of ymm1, then fill hi128 from xmm2/m128 via (with imm8 = 1):

	VINSERTF128 imm8,src2[xmm/m128],src1[ymm1],dest[ymm2]:
	imm8 = 0: dest.lo128 = src2, dest.hi128 = src1.hi128
	imm8 = 1: dest.lo128 = src1.lo128, dest.hi128 = src2

Once have 4 dcomplex roots loaded into 2 ymm as
ymm0.lo,hi = [c0,s0,c2,s2]
ymm1.lo,hi = [c1,s1,c3,s3]
Then interleave via

vunpcklpd ymm0,ymm1,ymmA
vunpckhpd ymm0,ymm1,ymmB

to get

ymmA = [c0,c1,c2,c3]
ymmB = [s0,s1,s2,s3]

Similarly for table2 [ = rn1 ] roots to get:

ymmC = [x0,x1,x2,x3]
ymmD = [y0,y1,y2,y3]

then do CMUL:

vmulpd	ymmA,ymmD,ymmE	// ymmE = c.y
vmulpd	ymmA,ymmC,ymmA	// ymmA = c.x

vmulpd	ymmB,ymmC,ymmC	// ymmC = s.x
vmulpd	ymmB,ymmD,ymmD	// ymmD = s.y

vsubpd	ymmA,ymmD,ymmA	// ymmA = c.x - s.y; ymmD free
vsubpd	ymmC,ymmE,ymmB	// ymmB = s.x + c.y; ymmC,E free
*/
  // FMA-based versions of selected macros in this file for Intel AVX2/FMA3
  #ifdef USE_AVX2

	/* Power-of-2-runlength Fermat-mod acyclic-transform/IBDWT carry macro.

	NOTE: The array indices i/j/k/lcycle declared int in caller but assumed to have been
	byte-shift-converted at time this macro called, thus can use as complex-address-offsets.
	*/
	#define SSE2_fermat_carry_norm_pow2_errcheck_X4(Xdata,Xbase_root,Xcmul_offset,Xcy_re,Xcy_im,Xhalf_arr,Xsign_mask, Xadd0,Xp1,Xp2,Xp3)\
	{\
	__asm__ volatile (\
		"movq		%[__add0],%%r14		\n\t"/* base address for 4 prefetches-from-main-data-array spread through this macro */\
		/* Base negacyclic roots at this address in [0,2,4,6]*0x20 (Re parts), [1,3,5,7]*0x20 (Imag parts) */\
		"movq	%[__base_root] ,%%rax			\n\t"	/* Won't need main-array again until output transpose, so re-use rax for base_root */\
		"movslq	%[__cmul_offset],%%rbx			\n\t"\
		"addq	%%rax,%%rbx	\n\t"	/* Index into complex const multipliers block, each applied to 4 sets of base roots */\
		/* Up-multiply quartet of negacyclic roots used in this macro invocation; store sets 2-4 back into mem, keep set 1 in ymm10,11 [that's why we do sets 1/2 after 3/4] */\
		"vmovaps	    (%%rbx),%%ymm10		\n\t	"	/* Multiply by exp(j*I*Pi/2)/RADIX, for j = 0-3 */\
		"vmovaps	0x20(%%rbx),%%ymm11		\n\t	"	/* c = Re(exp) in ymm0, s = Im(exp) in ymm1 */\
		"/* Sets 3/4: */"\
		"vmovaps	 0x80(%%rax),%%ymm0		\n\t	vmovaps	 0xc0(%%rax),%%ymm2			\n\t"	\
		"vmovaps	 0xa0(%%rax),%%ymm1		\n\t	vmovaps	 0xe0(%%rax),%%ymm3			\n\t"	\
		"vmovaps	%%ymm0,%%ymm4			\n\t	vmovaps		%%ymm2,%%ymm6			\n\t"	\
		"vmovaps	%%ymm1,%%ymm5			\n\t	vmovaps		%%ymm3,%%ymm7			\n\t"	\
		"vmulpd		%%ymm10,%%ymm0,%%ymm0	\n\t	vmulpd		%%ymm10,%%ymm2,%%ymm2	\n\t"	\
		"vmulpd		%%ymm11,%%ymm5,%%ymm5	\n\t	vmulpd		%%ymm11,%%ymm7,%%ymm7	\n\t"	\
		"vmulpd		%%ymm10,%%ymm1,%%ymm1	\n\t	vmulpd		%%ymm10,%%ymm3,%%ymm3	\n\t"	\
		"vmulpd		%%ymm11,%%ymm4,%%ymm4	\n\t	vmulpd		%%ymm11,%%ymm6,%%ymm6	\n\t"	\
		"vsubpd		%%ymm5 ,%%ymm0,%%ymm0	\n\t	vsubpd		%%ymm7 ,%%ymm2,%%ymm2	\n\t"	\
		"vaddpd		%%ymm4 ,%%ymm1,%%ymm1	\n\t	vaddpd		%%ymm6 ,%%ymm3,%%ymm3	\n\t"	\
		"vmovaps	%%ymm0 ,0x80(%%rax)		\n\t	vmovaps		%%ymm2 ,0xc0(%%rax)		\n\t"	/* Store result, overwriting input base root */\
		"vmovaps	%%ymm1 ,0xa0(%%rax)		\n\t	vmovaps		%%ymm3 ,0xe0(%%rax)		\n\t"	/* Im part */\
		"/* Sets 1/2: */"\
		"vmovaps	     (%%rax),%%ymm0		\n\t	vmovaps	 0x40(%%rax),%%ymm2			\n\t"	/* x = Re part of 1st base-root quartet */\
		"vmovaps	 0x20(%%rax),%%ymm1		\n\t	vmovaps	 0x60(%%rax),%%ymm3			\n\t"	/* y = Im part */\
		"vmovaps	%%ymm0,%%ymm4			\n\t	vmovaps		%%ymm2,%%ymm6			\n\t"	/* Copy x */\
		"vmovaps	%%ymm1,%%ymm5			\n\t	vmovaps		%%ymm3,%%ymm7			\n\t"	/* Copy y */\
		"vmulpd		%%ymm10,%%ymm0,%%ymm0	\n\t	vmulpd		%%ymm10,%%ymm2,%%ymm2	\n\t"	/* c.x */\
		"vmulpd		%%ymm11,%%ymm5,%%ymm5	\n\t	vmulpd		%%ymm11,%%ymm7,%%ymm7	\n\t"	/* s.y */\
		"vmulpd		%%ymm10,%%ymm1,%%ymm1	\n\t	vmulpd		%%ymm10,%%ymm3,%%ymm3	\n\t"	/* c.y */\
		"vmulpd		%%ymm11,%%ymm4,%%ymm4	\n\t	vmulpd		%%ymm11,%%ymm6,%%ymm6	\n\t"	/* s.x */\
		"vsubpd		%%ymm5 ,%%ymm0,%%ymm10	\n\t	vsubpd		%%ymm7 ,%%ymm2,%%ymm2	\n\t"	/* Out.re = c.x - s.y */\
		"vaddpd		%%ymm4 ,%%ymm1,%%ymm11	\n\t	vaddpd		%%ymm6 ,%%ymm3,%%ymm3	\n\t"	/* Out.im = c.y + s.x */\
		"											vmovaps		%%ymm2 ,0x40(%%rax)		\n\t"	/* Store result, overwriting input base root */\
		"											vmovaps		%%ymm3 ,0x60(%%rax)		\n\t"	/* Im part */\
		"movq		%[__data],%%rax		\n\t"\
		"movq		%[__half_arr],%%rdx		\n\t"\
		"vmovaps	0x40(%%rdx),%%ymm12	\n\t"	/* XMM12 = scale */\
	/* 4-way transpose of inputs (Re, Im parts separately): Inputs from r0/1,2/3,4/5.6/7. Outputs into ymm0-7: */\
	/* Real parts use ymm0,2,4,6, ymm8 as tmp-reg:					Imag parts use ymm1,3,5,7, ymm9 as tm-reg: */\
		"vmovaps	     (%%rax),%%ymm4						\n\t		vmovaps	0x020(%%rax),%%ymm5							\n\t"\
		"vmovaps	0x040(%%rax),%%ymm2						\n\t		vmovaps	0x060(%%rax),%%ymm3							\n\t"\
		"vshufpd	$15,%%ymm2,%%ymm4,%%ymm6				\n\t		vshufpd	$15,%%ymm3,%%ymm5,%%ymm7					\n\t"\
		"vshufpd	$0 ,%%ymm2,%%ymm4,%%ymm4				\n\t		vshufpd	$0 ,%%ymm3,%%ymm5,%%ymm5					\n\t"\
		"vmovaps	0x080(%%rax),%%ymm8						\n\t		vmovaps	0x0a0(%%rax),%%ymm9							\n\t"\
		"vmovaps	0x0c0(%%rax),%%ymm2						\n\t		vmovaps	0x0e0(%%rax),%%ymm3							\n\t"\
		"vshufpd	$15,%%ymm2,%%ymm8,%%ymm0				\n\t		vshufpd	$15,%%ymm3,%%ymm9,%%ymm1					\n\t"\
		"vshufpd	$0 ,%%ymm2,%%ymm8,%%ymm8				\n\t		vshufpd	$0 ,%%ymm3,%%ymm9,%%ymm9					\n\t"\
		"vperm2f128 $32,%%ymm0,%%ymm6,%%ymm2	/* Re B	*/	\n\t		vperm2f128 $32,%%ymm1,%%ymm7,%%ymm3		/* Im B	*/	\n\t"\
		"vperm2f128 $49,%%ymm0,%%ymm6,%%ymm6	/* Re D	*/	\n\t		vperm2f128 $49,%%ymm1,%%ymm7,%%ymm7		/* Im D	*/	\n\t"\
		"vperm2f128 $32,%%ymm8,%%ymm4,%%ymm0	/* Re A	*/	\n\t		vperm2f128 $32,%%ymm9,%%ymm5,%%ymm1 	/* Im A	*/	\n\t"\
		"vperm2f128 $49,%%ymm8,%%ymm4,%%ymm4	/* Re C	*/	\n\t		vperm2f128 $49,%%ymm9,%%ymm5,%%ymm5 	/* Im C	*/	\n\t"\
	/* Apply inverse-complex-runlength scaling factor to the data: */\
		"vmulpd		%%ymm12,%%ymm2,%%ymm2					\n\t		vmulpd		%%ymm12,%%ymm3,%%ymm3	\n\t"\
		"vmulpd		%%ymm12,%%ymm6,%%ymm6					\n\t		vmulpd		%%ymm12,%%ymm7,%%ymm7	\n\t"\
		"vmulpd		%%ymm12,%%ymm0,%%ymm0					\n\t		vmulpd		%%ymm12,%%ymm1,%%ymm1	\n\t"\
		"vmulpd		%%ymm12,%%ymm4,%%ymm4					\n\t		vmulpd		%%ymm12,%%ymm5,%%ymm5	\n\t"\
		/* Base negacyclic roots at this address in [0,2,4,6]*0x20 (Re parts), [1,3,5,7]*0x20 (Imag parts) */\
		"movq	%[__base_root] ,%%rax		\n\t"	/* Won't need main-array again until output transpose, so re-use rax for base_root */\
		/* half_arr[0,1,2,3] = [base*2, baseinv*2,wt_re*2,wt_im*2] */\
		"vmovaps	-0x40(%%rdx),%%ymm13	\n\t"	/* XMM13 = maxerr */\
		"/*vmovaps	-0x20(%%rdx),%%ymm15	*/\n\t"	/* rnd_const; prefer ROUNDPD in AVX mode, so ymm15 free */\
		"vmovaps	%%ymm13,%%ymm14			\n\t"	/* maxerr copy */\
		"movq		%[__cy_re],%%rbx		\n\t"\
		"movq		%[__cy_im],%%rcx		\n\t"\
	/* Do a-quartet: Data in ymm0,ymm1 */\
	"prefetcht0	(%%r14)	\n\t"\
		/* For a-quartet, needed negacyclic root already in ymm10/11: */\
		/* Data in ymm0,ymm1 */\
		"vmovaps	%%ymm0,%%ymm8			\n\t"	/* x copy */\
		"vmovaps	%%ymm1,%%ymm9			\n\t"	/* y copy */\
		/* Inverse negacyclic weight is (wt_re, -wt_im): */\
		"vmulpd		%%ymm10,%%ymm0,%%ymm0	\n\t"	/* wt_re*[x     ] */\
		"vmulpd		%%ymm11,%%ymm9,%%ymm9	\n\t"	/* wt_im*[y copy] */\
		"vmulpd		%%ymm10,%%ymm1,%%ymm1	\n\t"	/* wt_re*[y     ] */\
		"vmulpd		%%ymm11,%%ymm8,%%ymm8	\n\t"	/* wt_im*[x copy] */\
		"vaddpd		%%ymm9 ,%%ymm0,%%ymm0	\n\t"	/* [a0-3.re] = x*wt_re + y*wt_im */\
		"vsubpd		%%ymm8 ,%%ymm1,%%ymm1	\n\t"	/* [a0-3.im] = y*wt_re - x*wt_im */\
		/* normalize a-pair, compute carryout, compute ROE: */\
		"vaddpd		(%%rbx),%%ymm0,%%ymm0	\n\t"	/* [a0-3.re] + [cy0-3.re] */\
		"vaddpd		(%%rcx),%%ymm1,%%ymm1	\n\t"	/* [a0-3.im] + [cy0-3.im] */\
		"vmovaps	%%ymm0,%%ymm8			\n\t	vmovaps		%%ymm1,%%ymm9		\n\t"	/* copy x|y */\
		"vroundpd	$0,%%ymm0,%%ymm0		\n\t	vroundpd	$0,%%ymm1,%%ymm1	\n\t"	/* temp = DNINT(x|y) */\
		"movq		%[__sign_mask],%%rsi	\n\t"\
		"vmovaps	0x20(%%rdx),%%ymm12	\n\t"	/* [baseinv0-3] */\
		"vsubpd		%%ymm0 ,%%ymm8 ,%%ymm8 	\n\t	vsubpd		%%ymm1 ,%%ymm9 ,%%ymm9 	\n\t"	/* frac = [x - temp] */\
		"vandpd		(%%rsi),%%ymm8 ,%%ymm8 	\n\t	vandpd		(%%rsi),%%ymm9 ,%%ymm9 	\n\t"	/* frac = fabs(frac) */\
		"vmaxpd		%%ymm13,%%ymm8 ,%%ymm13	\n\t	vmaxpd		%%ymm14,%%ymm9 ,%%ymm14	\n\t"	/* if(frac > maxerr) maxerr=frac */\
		"vmovaps	%%ymm0 ,%%ymm8			\n\t	vmovaps		%%ymm1 ,%%ymm9			\n\t"	/* cpy temp */\
		"vmulpd		%%ymm12,%%ymm8 ,%%ymm8 	\n\t	vmulpd		%%ymm12,%%ymm9 ,%%ymm9 	\n\t"	/* temp*baseinv[0] ... inline the remaining +odd_radix offset in addressing */\
		"vmaxpd		%%ymm13,%%ymm14,%%ymm14	\n\t"	/* merge re|im maxerr vectors */\
		"vroundpd	$0,%%ymm8 ,%%ymm8 		\n\t	vroundpd	$0,%%ymm9,%%ymm9		\n\t"	/* [cy0-3.re] = DNINT(temp*baseinv[0]) */\
		"vmovaps	(%%rdx),%%ymm12	\n\t"	/* [base0-3] */\
		"vmovaps	%%ymm8,(%%rbx)			\n\t	vmovaps		%%ymm9,(%%rcx)			\n\t"	/* store [cy0-3.re|im] */\
		"vmulpd		%%ymm12,%%ymm8,%%ymm8	\n\t	vmulpd		%%ymm12,%%ymm9,%%ymm9	\n\t"	/* base[0]*[cy0-3.re|im] */\
		"vsubpd		%%ymm8 ,%%ymm0,%%ymm0	\n\t	vsubpd		%%ymm9 ,%%ymm1,%%ymm1	\n\t"	/* XMM0|1 = [a0-3.re|im] = temp - [cy0-3.re|im]*base[0] */\
		"vmovaps	%%ymm0 ,%%ymm8			\n\t	vmovaps		%%ymm1 ,%%ymm9			\n\t"	/* cpy x|y */\
		/* Forward acyclic-convo weight is (wt_re, +wt_im): */\
		"vmulpd		%%ymm10,%%ymm0,%%ymm0	\n\t"	/* wt_re*[x     ] */\
		"vmulpd		%%ymm11,%%ymm9,%%ymm9	\n\t"	/* wt_im*[y copy] */\
		"vmulpd		%%ymm10,%%ymm1,%%ymm1	\n\t"	/* wt_re*[y     ] */\
		"vmulpd		%%ymm11,%%ymm8,%%ymm8	\n\t"	/* wt_im*[x copy] */\
		"vsubpd		%%ymm9 ,%%ymm0,%%ymm0	\n\t"	/* [a0-3.re] = x*wt_re - y*wt_im */\
		"vaddpd		%%ymm8 ,%%ymm1,%%ymm1	\n\t"	/* [a0-3.im] = y*wt_re + x*wt_im */\
	/* Now do b-quartet: Data in ymm2,ymm3 */\
	"movslq		%[__p1],%%r15	\n\t"	\
	"prefetcht0	(%%r14,%%r15,8)	\n\t"\
		"vmovaps	%%ymm14,%%ymm13			\n\t"	/* maxerr copy */\
		/* Get next set of negacyclic roots: */\
		"vmovaps	0x40(%%rax),%%ymm10		\n\t"	/* c = Re part of 2nd base-root quartet */\
		"vmovaps	0x60(%%rax),%%ymm11		\n\t"	/* s = Im part */\
		/* Data in ymm2,ymm3 */\
		"vmovaps	%%ymm2,%%ymm8			\n\t"	/* x copy */\
		"vmovaps	%%ymm3,%%ymm9			\n\t"	/* y copy */\
		/* Inverse negacyclic weight is (wt_re, -wt_im): */\
		"vmulpd		%%ymm10,%%ymm2,%%ymm2	\n\t"	/* wt_re*[x     ] */\
		"vmulpd		%%ymm11,%%ymm9,%%ymm9	\n\t"	/* wt_im*[y copy] */\
		"vmulpd		%%ymm10,%%ymm3,%%ymm3	\n\t"	/* wt_re*[y     ] */\
		"vmulpd		%%ymm11,%%ymm8,%%ymm8	\n\t"	/* wt_im*[x copy] */\
		"vaddpd		%%ymm9 ,%%ymm2,%%ymm2	\n\t"	/* [a0-3.re] = x*wt_re + y*wt_im */\
		"vsubpd		%%ymm8 ,%%ymm3,%%ymm3	\n\t"	/* [a0-3.im] = y*wt_re - x*wt_im */\
		/* normalize, compute carryout, compute ROE: */\
		"vaddpd		(%%rbx),%%ymm2,%%ymm2	\n\t"	/* [a0-3.re] + [cy0-3.re] */\
		"vaddpd		(%%rcx),%%ymm3,%%ymm3	\n\t"	/* [a0-3.im] + [cy0-3.im] */\
		"vmovaps	%%ymm2,%%ymm8			\n\t	vmovaps		%%ymm3,%%ymm9		\n\t"	/* copy x|y */\
		"vroundpd	$0,%%ymm2,%%ymm2		\n\t	vroundpd	$0,%%ymm3,%%ymm3	\n\t"	/* temp = DNINT(x|y) */\
		"vmovaps	0x20(%%rdx),%%ymm12	\n\t"	/* [baseinv0-3] */\
		"vsubpd		%%ymm2 ,%%ymm8 ,%%ymm8 	\n\t	vsubpd		%%ymm3 ,%%ymm9 ,%%ymm9 	\n\t"	/* frac = [x - temp] */\
		"vandpd		(%%rsi),%%ymm8 ,%%ymm8 	\n\t	vandpd		(%%rsi),%%ymm9 ,%%ymm9 	\n\t"	/* frac = fabs(frac) */\
		"vmaxpd		%%ymm13,%%ymm8 ,%%ymm13	\n\t	vmaxpd		%%ymm14,%%ymm9 ,%%ymm14	\n\t"	/* if(frac > maxerr) maxerr=frac */\
		"vmovaps	%%ymm2 ,%%ymm8			\n\t	vmovaps		%%ymm3 ,%%ymm9			\n\t"	/* cpy temp */\
		"vmulpd		%%ymm12,%%ymm8 ,%%ymm8 	\n\t	vmulpd		%%ymm12,%%ymm9 ,%%ymm9 	\n\t"	/* temp*baseinv[0] ... inline the remaining +odd_radix offset in addressing */\
		"vmaxpd		%%ymm13,%%ymm14,%%ymm14	\n\t"	/* merge re|im maxerr vectors */\
		"vroundpd	$0,%%ymm8 ,%%ymm8 		\n\t	vroundpd	$0,%%ymm9,%%ymm9		\n\t"	/* [cy0-3.re] = DNINT(temp*baseinv[0]) */\
		"vmovaps	(%%rdx),%%ymm12	\n\t"	/* [base0-3] */\
		"vmovaps	%%ymm8,(%%rbx)			\n\t	vmovaps		%%ymm9,(%%rcx)			\n\t"	/* store [cy0-3.re|im] */\
		"vmulpd		%%ymm12,%%ymm8,%%ymm8	\n\t	vmulpd		%%ymm12,%%ymm9,%%ymm9	\n\t"	/* base[0]*[cy0-3.re|im] */\
		"vsubpd		%%ymm8 ,%%ymm2,%%ymm2	\n\t	vsubpd		%%ymm9 ,%%ymm3,%%ymm3	\n\t"	/* XMM0|1 = [a0-3.re|im] = temp - [cy0-3.re|im]*base[0] */\
		"vmovaps	%%ymm2 ,%%ymm8			\n\t	vmovaps		%%ymm3 ,%%ymm9			\n\t"	/* cpy x|y */\
		/* Forward acyclic-convo weight is (wt_re, +wt_im): */\
		"vmulpd		%%ymm10,%%ymm2,%%ymm2	\n\t"	/* wt_re*[x     ] */\
		"vmulpd		%%ymm11,%%ymm9,%%ymm9	\n\t"	/* wt_im*[y copy] */\
		"vmulpd		%%ymm10,%%ymm3,%%ymm3	\n\t"	/* wt_re*[y     ] */\
		"vmulpd		%%ymm11,%%ymm8,%%ymm8	\n\t"	/* wt_im*[x copy] */\
		"vsubpd		%%ymm9 ,%%ymm2,%%ymm2	\n\t"	/* [a0-3.re] = x*wt_re - y*wt_im */\
		"vaddpd		%%ymm8 ,%%ymm3,%%ymm3	\n\t"	/* [a0-3.im] = y*wt_re + x*wt_im */\
	/* Now do c-quartet: Data in ymm4,ymm5 */\
	"movslq		%[__p2],%%r15	\n\t"	\
	"prefetcht0	(%%r14,%%r15,8)	\n\t"\
		"vmovaps	%%ymm14,%%ymm13			\n\t"	/* maxerr copy */\
		/* Get next set of negacyclic roots: */\
		"vmovaps	0x80(%%rax),%%ymm10		\n\t"	/* c = Re part of 3rd base-root quartet */\
		"vmovaps	0xa0(%%rax),%%ymm11		\n\t"	/* s = Im part */\
		/* Data in ymm4,ymm5 */\
		"vmovaps	%%ymm4,%%ymm8			\n\t"	/* x copy */\
		"vmovaps	%%ymm5,%%ymm9			\n\t"	/* y copy */\
		/* Inverse negacyclic weight is (wt_re, -wt_im): */\
		"vmulpd		%%ymm10,%%ymm4,%%ymm4	\n\t"	/* wt_re*[x     ] */\
		"vmulpd		%%ymm11,%%ymm9,%%ymm9	\n\t"	/* wt_im*[y copy] */\
		"vmulpd		%%ymm10,%%ymm5,%%ymm5	\n\t"	/* wt_re*[y     ] */\
		"vmulpd		%%ymm11,%%ymm8,%%ymm8	\n\t"	/* wt_im*[x copy] */\
		"vaddpd		%%ymm9 ,%%ymm4,%%ymm4	\n\t"	/* [a0-3.re] = x*wt_re + y*wt_im */\
		"vsubpd		%%ymm8 ,%%ymm5,%%ymm5	\n\t"	/* [a0-3.im] = y*wt_re - x*wt_im */\
		/* normalize, compute carryout, compute ROE: */\
		"vaddpd		(%%rbx),%%ymm4,%%ymm4	\n\t"	/* [a0-3.re] + [cy0-3.re] */\
		"vaddpd		(%%rcx),%%ymm5,%%ymm5	\n\t"	/* [a0-3.im] + [cy0-3.im] */\
		"vmovaps	%%ymm4,%%ymm8			\n\t	vmovaps		%%ymm5,%%ymm9		\n\t"	/* copy x|y */\
		"vroundpd	$0,%%ymm4,%%ymm4		\n\t	vroundpd	$0,%%ymm5,%%ymm5	\n\t"	/* temp = DNINT(x|y) */\
		"vmovaps	0x20(%%rdx),%%ymm12	\n\t"	/* [baseinv0-3] */\
		"vsubpd		%%ymm4 ,%%ymm8 ,%%ymm8 	\n\t	vsubpd		%%ymm5 ,%%ymm9 ,%%ymm9 	\n\t"	/* frac = [x - temp] */\
		"vandpd		(%%rsi),%%ymm8 ,%%ymm8 	\n\t	vandpd		(%%rsi),%%ymm9 ,%%ymm9 	\n\t"	/* frac = fabs(frac) */\
		"vmaxpd		%%ymm13,%%ymm8 ,%%ymm13	\n\t	vmaxpd		%%ymm14,%%ymm9 ,%%ymm14	\n\t"	/* if(frac > maxerr) maxerr=frac */\
		"vmovaps	%%ymm4 ,%%ymm8			\n\t	vmovaps		%%ymm5 ,%%ymm9			\n\t"	/* cpy temp */\
		"vmulpd		%%ymm12,%%ymm8 ,%%ymm8 	\n\t	vmulpd		%%ymm12,%%ymm9 ,%%ymm9 	\n\t"	/* temp*baseinv[0] ... inline the remaining +odd_radix offset in addressing */\
		"vmaxpd		%%ymm13,%%ymm14,%%ymm14	\n\t"	/* merge re|im maxerr vectors */\
		"vroundpd	$0,%%ymm8 ,%%ymm8 		\n\t	vroundpd	$0,%%ymm9,%%ymm9		\n\t"	/* [cy0-3.re] = DNINT(temp*baseinv[0]) */\
		"vmovaps	(%%rdx),%%ymm12	\n\t"	/* [base0-3] */\
		"vmovaps	%%ymm8,(%%rbx)			\n\t	vmovaps		%%ymm9,(%%rcx)			\n\t"	/* store [cy0-3.re|im] */\
		"vmulpd		%%ymm12,%%ymm8,%%ymm8	\n\t	vmulpd		%%ymm12,%%ymm9,%%ymm9	\n\t"	/* base[0]*[cy0-3.re|im] */\
		"vsubpd		%%ymm8 ,%%ymm4,%%ymm4	\n\t	vsubpd		%%ymm9 ,%%ymm5,%%ymm5	\n\t"	/* XMM0|1 = [a0-3.re|im] = temp - [cy0-3.re|im]*base[0] */\
		"vmovaps	%%ymm4 ,%%ymm8			\n\t	vmovaps		%%ymm5 ,%%ymm9			\n\t"	/* cpy x|y */\
		/* Forward acyclic-convo weight is (wt_re, +wt_im): */\
		"vmulpd		%%ymm10,%%ymm4,%%ymm4	\n\t"	/* wt_re*[x     ] */\
		"vmulpd		%%ymm11,%%ymm9,%%ymm9	\n\t"	/* wt_im*[y copy] */\
		"vmulpd		%%ymm10,%%ymm5,%%ymm5	\n\t"	/* wt_re*[y     ] */\
		"vmulpd		%%ymm11,%%ymm8,%%ymm8	\n\t"	/* wt_im*[x copy] */\
		"vsubpd		%%ymm9 ,%%ymm4,%%ymm4	\n\t"	/* [a0-3.re] = x*wt_re - y*wt_im */\
		"vaddpd		%%ymm8 ,%%ymm5,%%ymm5	\n\t"	/* [a0-3.im] = y*wt_re + x*wt_im */\
	/* Now do d-quartet: Data in ymm6,ymm7 */\
	"movslq		%[__p3],%%r15	\n\t"	\
	"prefetcht0	(%%r14,%%r15,8)	\n\t"\
		"vmovaps	%%ymm14,%%ymm13			\n\t"	/* maxerr copy */\
		/* Get next set of negacyclic roots: */\
		"vmovaps	0xc0(%%rax),%%ymm10		\n\t"	/* c = Re part of 4th base-root quartet */\
		"vmovaps	0xe0(%%rax),%%ymm11		\n\t"	/* s = Im part */\
		/* Data in ymm6,ymm7 */\
		"vmovaps	%%ymm6,%%ymm8			\n\t"	/* x copy */\
		"vmovaps	%%ymm7,%%ymm9			\n\t"	/* y copy */\
		/* Inverse negacyclic weight is (wt_re, -wt_im): */\
		"vmulpd		%%ymm10,%%ymm6,%%ymm6	\n\t"	/* wt_re*[x     ] */\
		"vmulpd		%%ymm11,%%ymm9,%%ymm9	\n\t"	/* wt_im*[y copy] */\
		"vmulpd		%%ymm10,%%ymm7,%%ymm7	\n\t"	/* wt_re*[y     ] */\
		"vmulpd		%%ymm11,%%ymm8,%%ymm8	\n\t"	/* wt_im*[x copy] */\
		"vaddpd		%%ymm9 ,%%ymm6,%%ymm6	\n\t"	/* [a0-3.re] = x*wt_re + y*wt_im */\
		"vsubpd		%%ymm8 ,%%ymm7,%%ymm7	\n\t"	/* [a0-3.im] = y*wt_re - x*wt_im */\
		/* normalize, compute carryout, compute ROE: */\
		"vaddpd		(%%rbx),%%ymm6,%%ymm6	\n\t"	/* [a0-3.re] + [cy0-3.re] */\
		"vaddpd		(%%rcx),%%ymm7,%%ymm7	\n\t"	/* [a0-3.im] + [cy0-3.im] */\
		"vmovaps	%%ymm6,%%ymm8			\n\t	vmovaps		%%ymm7,%%ymm9		\n\t"	/* copy x|y */\
		"vroundpd	$0,%%ymm6,%%ymm6		\n\t	vroundpd	$0,%%ymm7,%%ymm7	\n\t"	/* temp = DNINT(x|y) */\
		"vmovaps	0x20(%%rdx),%%ymm12	\n\t"	/* [baseinv0-3] */\
		"vsubpd		%%ymm6 ,%%ymm8 ,%%ymm8 	\n\t	vsubpd		%%ymm7 ,%%ymm9 ,%%ymm9 	\n\t"	/* frac = [x - temp] */\
		"vandpd		(%%rsi),%%ymm8 ,%%ymm8 	\n\t	vandpd		(%%rsi),%%ymm9 ,%%ymm9 	\n\t"	/* frac = fabs(frac) */\
		"vmaxpd		%%ymm13,%%ymm8 ,%%ymm13	\n\t	vmaxpd		%%ymm14,%%ymm9 ,%%ymm14	\n\t"	/* if(frac > maxerr) maxerr=frac */\
		"vmovaps	%%ymm6 ,%%ymm8			\n\t	vmovaps		%%ymm7 ,%%ymm9			\n\t"	/* cpy temp */\
		"vmulpd		%%ymm12,%%ymm8 ,%%ymm8 	\n\t	vmulpd		%%ymm12,%%ymm9 ,%%ymm9 	\n\t"	/* temp*baseinv[0] ... inline the remaining +odd_radix offset in addressing */\
		"vmaxpd		%%ymm13,%%ymm14,%%ymm14	\n\t"	/* merge re|im maxerr vectors */\
		"vroundpd	$0,%%ymm8 ,%%ymm8 		\n\t	vroundpd	$0,%%ymm9,%%ymm9		\n\t"	/* [cy0-3.re] = DNINT(temp*baseinv[0]) */\
		"vmovaps	(%%rdx),%%ymm12	\n\t"	/* [base0-3] */\
		"vmovaps	%%ymm8,(%%rbx)			\n\t	vmovaps		%%ymm9,(%%rcx)			\n\t"	/* store [cy0-3.re|im] */\
		"vmulpd		%%ymm12,%%ymm8,%%ymm8	\n\t	vmulpd		%%ymm12,%%ymm9,%%ymm9	\n\t"	/* base[0]*[cy0-3.re|im] */\
		"vsubpd		%%ymm8 ,%%ymm6,%%ymm6	\n\t	vsubpd		%%ymm9 ,%%ymm7,%%ymm7	\n\t"	/* XMM0|1 = [a0-3.re|im] = temp - [cy0-3.re|im]*base[0] */\
		"vmovaps	%%ymm6 ,%%ymm8			\n\t	vmovaps		%%ymm7 ,%%ymm9			\n\t"	/* cpy x|y */\
		/* Forward acyclic-convo weight is (wt_re, +wt_im): */\
		"vmulpd		%%ymm10,%%ymm6,%%ymm6	\n\t"	/* wt_re*[x     ] */\
		"vmulpd		%%ymm11,%%ymm9,%%ymm9	\n\t"	/* wt_im*[y copy] */\
		"vmulpd		%%ymm10,%%ymm7,%%ymm7	\n\t"	/* wt_re*[y     ] */\
		"vmulpd		%%ymm11,%%ymm8,%%ymm8	\n\t"	/* wt_im*[x copy] */\
		"vsubpd		%%ymm9 ,%%ymm6,%%ymm6	\n\t"	/* [a0-3.re] = x*wt_re - y*wt_im */\
		"vaddpd		%%ymm8 ,%%ymm7,%%ymm7	\n\t"	/* [a0-3.im] = y*wt_re + x*wt_im */\
		/* Store maxerr: */\
		"movq		%[__half_arr],%%rdx		\n\t"\
		"vmovaps	%%ymm14,-0x40(%%rdx)	\n\t"\
	/* 4-way transpose of outputs (Re, Im parts separately): Inputs from ymm0-7. Outputs into r0/1,2/3,4/5.6/7: */	\
	/* Because default inputs for our 4 x 4 transpose macro (e.g. the one used at start of this carry macro) */\
	/* are into ymm4/2/8/2, munge inputs into that order, resolving name-conflicts via use of the now-available ymm8-15 for outputs: */\
		"movq		%[__data],%%rax			\n\t"\
		"vshufpd	$15,%%ymm2,%%ymm0,%%ymm10					\n\t		vshufpd	$15,%%ymm3,%%ymm1,%%ymm11						\n\t"\
		"vshufpd	$0 ,%%ymm2,%%ymm0,%%ymm0					\n\t		vshufpd	$0 ,%%ymm3,%%ymm1,%%ymm1						\n\t"\
		"vshufpd	$15,%%ymm6,%%ymm4,%%ymm12					\n\t		vshufpd	$15,%%ymm7,%%ymm5,%%ymm13						\n\t"\
		"vshufpd	$0 ,%%ymm6,%%ymm4,%%ymm4					\n\t		vshufpd	$0 ,%%ymm7,%%ymm5,%%ymm5						\n\t"\
		"vperm2f128 $32,%%ymm12,%%ymm10,%%ymm2 		/* Re B	*/	\n\t		vperm2f128 $32,%%ymm13,%%ymm11,%%ymm3		/* Im B	*/	\n\t"\
		"vperm2f128 $49,%%ymm12,%%ymm10,%%ymm10		/* Re D	*/	\n\t		vperm2f128 $49,%%ymm13,%%ymm11,%%ymm11		/* Im D	*/	\n\t"\
		"vperm2f128 $32,%%ymm4 ,%%ymm0 ,%%ymm12		/* Re A	*/	\n\t		vperm2f128 $32,%%ymm5 ,%%ymm1 ,%%ymm13 		/* Im A	*/	\n\t"\
		"vperm2f128 $49,%%ymm4 ,%%ymm0 ,%%ymm0 		/* Re C	*/	\n\t		vperm2f128 $49,%%ymm5 ,%%ymm1 ,%%ymm1		/* Im C	*/	\n\t"\
		/* And write 'em back to memory: */\
		"vmovaps	%%ymm12,     (%%rax)						\n\t		vmovaps	%%ymm13,0x020(%%rax)				\n\t"\
		"vmovaps	%%ymm2 ,0x040(%%rax)						\n\t		vmovaps	%%ymm3 ,0x060(%%rax)				\n\t"\
		"vmovaps	%%ymm0 ,0x080(%%rax)						\n\t		vmovaps	%%ymm1 ,0x0a0(%%rax)				\n\t"\
		"vmovaps	%%ymm10,0x0c0(%%rax)						\n\t		vmovaps	%%ymm11,0x0e0(%%rax)				\n\t"\
		:						/* outputs: none */\
		:	[__data]		"m" (Xdata)	/* All inputs from memory addresses here */\
		,	[__base_root]	"m" (Xbase_root)\
		,	[__cmul_offset] "m" (Xcmul_offset)\
		,	[__cy_re]		"m" (Xcy_re)\
		,	[__cy_im]		"m" (Xcy_im)\
		,	[__half_arr]	"m" (Xhalf_arr)\
		,	[__sign_mask]	"m" (Xsign_mask)\
		/* Prefetch: base address and 3 index offsets */\
		,	[__add0] "m" (Xadd0)\
		,	[__p1] "m" (Xp1)\
		,	[__p2] "m" (Xp2)\
		,	[__p3] "m" (Xp3)\
		: "cc","memory","rax","rbx","rcx","rdx","rsi","r14","r15","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14"	/* Clobbered registers */\
	);\
	}

	/* Non-power-of-2-runlength Fermat-mod acyclic-transform/IBDWT carry macro.

	Key differences vs pow2 version:
	- Use odd_radix as index offset into local storage for IBDWT weights and variable base/baseinv terms;
	- Apply inv/fwd IBDWT weights bookending the negacyclic weights;
	- Value of base/baseinv to be applied to output taken from odd_radix-length array, using same index as for selecting IBDWT weight.

	The array indices i/j/k/lcycle declared int in caller but assumed to have been byte-shift-converted at time this macro called,
	thus can use as complex-address-offsets.  Use bytewise literal offsets to save registers for several args here,as vvv-marked:
												                                           vvvvv The [1,2,3]-multiples of odd_radix assumed << l2_sz_vd on input */
	#define SSE2_fermat_carry_norm_errcheck_X4_hiacc(Xdata,Xbase_root,Xcmul_offset,Xcy_re,Xcy_im,Xodd_radix,Xodd_radm2,Xodd_radm3,Xhalf_arr,Xsign_mask,XicycleA,XicycleB,XicycleC,XicycleD, XjcycleA,XkcycleA,XlcycleA, Xadd0,Xp1,Xp2,Xp3)\
	{\
	__asm__ volatile (\
		"movq		%[__add0],%%rcx		\n\t"/* base address for 4 prefetches-from-main-data-array spread through this macro */\
		/* Base negacyclic roots at this address in [0,2,4,6]*0x20 (Re parts), [1,3,5,7]*0x20 (Imag parts) */\
		"movq	%[__base_root]  ,%%rax			\n\t"	/* Won't need main-array again until output transpose, so re-use rax for base_root */\
		"movslq	%[__cmul_offset],%%rbx			\n\t"\
		"addq	%%rax,%%rbx	\n\t"	/* Index into complex const multipliers block, each applied to 4 sets of base roots */\
		/* Up-multiply quartet of negacyclic roots used in this macro invocation; store sets 2-4 back into mem, keep set 1 in ymm10,11 [that's why we do sets 1/2 after 3/4] */\
		"vmovaps	    (%%rbx),%%ymm10		\n\t	"	/* Multiply by exp(j*I*Pi/2)/RADIX, for j = 0-3 */\
		"vmovaps	0x20(%%rbx),%%ymm11		\n\t	"	/* c = Re(exp) in ymm0, s = Im(exp) in ymm1 */\
		"/* Sets 3/4: */"\
		"vmovaps	 0x80(%%rax),%%ymm0		\n\t	vmovaps	 0xc0(%%rax),%%ymm2			\n\t"	\
		"vmovaps	 0xa0(%%rax),%%ymm1		\n\t	vmovaps	 0xe0(%%rax),%%ymm3			\n\t"	\
		"vmovaps	%%ymm0,%%ymm4			\n\t	vmovaps		%%ymm2,%%ymm6			\n\t"	\
		"vmovaps	%%ymm1,%%ymm5			\n\t	vmovaps		%%ymm3,%%ymm7			\n\t"	\
		"vmulpd		%%ymm10,%%ymm0,%%ymm0	\n\t	vmulpd		%%ymm10,%%ymm2,%%ymm2	\n\t"	\
		"vmulpd		%%ymm11,%%ymm5,%%ymm5	\n\t	vmulpd		%%ymm11,%%ymm7,%%ymm7	\n\t"	\
		"vmulpd		%%ymm10,%%ymm1,%%ymm1	\n\t	vmulpd		%%ymm10,%%ymm3,%%ymm3	\n\t"	\
		"vmulpd		%%ymm11,%%ymm4,%%ymm4	\n\t	vmulpd		%%ymm11,%%ymm6,%%ymm6	\n\t"	\
		"vsubpd		%%ymm5 ,%%ymm0,%%ymm0	\n\t	vsubpd		%%ymm7 ,%%ymm2,%%ymm2	\n\t"	\
		"vaddpd		%%ymm4 ,%%ymm1,%%ymm1	\n\t	vaddpd		%%ymm6 ,%%ymm3,%%ymm3	\n\t"	\
		"vmovaps	%%ymm0 ,0x80(%%rax)		\n\t	vmovaps		%%ymm2 ,0xc0(%%rax)		\n\t"	/* Store result, overwriting input base root */\
		"vmovaps	%%ymm1 ,0xa0(%%rax)		\n\t	vmovaps		%%ymm3 ,0xe0(%%rax)		\n\t"	/* Im part */\
		"/* Sets 1/2: */"\
		"vmovaps	     (%%rax),%%ymm0		\n\t	vmovaps	 0x40(%%rax),%%ymm2			\n\t"	/* x = Re part of 1st base-root quartet */\
		"vmovaps	 0x20(%%rax),%%ymm1		\n\t	vmovaps	 0x60(%%rax),%%ymm3			\n\t"	/* y = Im part */\
		"vmovaps	%%ymm0,%%ymm4			\n\t	vmovaps		%%ymm2,%%ymm6			\n\t"	/* Copy x */\
		"vmovaps	%%ymm1,%%ymm5			\n\t	vmovaps		%%ymm3,%%ymm7			\n\t"	/* Copy y */\
		"vmulpd		%%ymm10,%%ymm0,%%ymm0	\n\t	vmulpd		%%ymm10,%%ymm2,%%ymm2	\n\t"	/* c.x */\
		"vmulpd		%%ymm11,%%ymm5,%%ymm5	\n\t	vmulpd		%%ymm11,%%ymm7,%%ymm7	\n\t"	/* s.y */\
		"vmulpd		%%ymm10,%%ymm1,%%ymm1	\n\t	vmulpd		%%ymm10,%%ymm3,%%ymm3	\n\t"	/* c.y */\
		"vmulpd		%%ymm11,%%ymm4,%%ymm4	\n\t	vmulpd		%%ymm11,%%ymm6,%%ymm6	\n\t"	/* s.x */\
		"vsubpd		%%ymm5 ,%%ymm0,%%ymm10	\n\t	vsubpd		%%ymm7 ,%%ymm2,%%ymm2	\n\t"	/* Out.re = c.x - s.y */\
		"vaddpd		%%ymm4 ,%%ymm1,%%ymm11	\n\t	vaddpd		%%ymm6 ,%%ymm3,%%ymm3	\n\t"	/* Out.im = c.y + s.x */\
		"											vmovaps		%%ymm2 ,0x40(%%rax)		\n\t"	/* Store result, overwriting input base root */\
		"											vmovaps		%%ymm3 ,0x60(%%rax)		\n\t"	/* Im part */\
		"movq		%[__data],%%rax		\n\t"\
		"movq		%[__half_arr],%%rdx		\n\t"\
		"addq		$%c[__odd_radix],%%rdx				\n\t"	/* wt|wtinv|base|baseinv data offset by icycle array slots from resp. base addresses */\
		/* Multiply complex transform outputs [x,y] = [re,im] by inverse IBDWT weights, which include the 2/n scale factor: */\
		"movslq		%[__icycleA],%%rdi		\n\t"	\
		"movslq		%[__icycleB],%%r9 		\n\t"	\
		"movslq		%[__icycleC],%%r8 		\n\t"	\
		"movslq		%[__icycleD],%%r10		\n\t"	\
		"vmovaps	(%%rdx,%%rdi),%%ymm12	\n\t"	/* [wtinv0-3]A */\
		"vmovaps	(%%rdx,%%r9 ),%%ymm13	\n\t"	/* [wtinv0-3]B */\
		"vmovaps	(%%rdx,%%r8 ),%%ymm14	\n\t"	/* [wtinv0-3]C */\
		"vmovaps	(%%rdx,%%r10),%%ymm15	\n\t"	/* [wtinv0-3]D */\
	/* 4-way transpose of inputs (Re, Im parts separately): Inputs from r0/1,2/3,4/5.6/7. Outputs into ymm0-7: */\
	/* Real parts use ymm0,2,4,6, ymm8 as tmp-reg:					Imag parts use ymm1,3,5,7, ymm9 as tm-reg: */\
		"vmovaps	     (%%rax),%%ymm4						\n\t		vmovaps	0x020(%%rax),%%ymm5							\n\t"\
		"vmovaps	0x040(%%rax),%%ymm2						\n\t		vmovaps	0x060(%%rax),%%ymm3							\n\t"\
		"vmulpd		%%ymm12,%%ymm4,%%ymm4					\n\t		vmulpd		%%ymm12,%%ymm5,%%ymm5					\n\t"\
		"vmulpd		%%ymm13,%%ymm2,%%ymm2					\n\t		vmulpd		%%ymm13,%%ymm3,%%ymm3					\n\t"\
		"vshufpd	$15,%%ymm2,%%ymm4,%%ymm6				\n\t		vshufpd	$15,%%ymm3,%%ymm5,%%ymm7					\n\t"\
		"vshufpd	$0 ,%%ymm2,%%ymm4,%%ymm4				\n\t		vshufpd	$0 ,%%ymm3,%%ymm5,%%ymm5					\n\t"\
		"vmovaps	0x080(%%rax),%%ymm8						\n\t		vmovaps	0x0a0(%%rax),%%ymm9							\n\t"\
		"vmovaps	0x0c0(%%rax),%%ymm2						\n\t		vmovaps	0x0e0(%%rax),%%ymm3							\n\t"\
		"vmulpd		%%ymm14,%%ymm8,%%ymm8					\n\t		vmulpd		%%ymm14,%%ymm9,%%ymm9					\n\t"\
		"vmulpd		%%ymm15,%%ymm2,%%ymm2					\n\t		vmulpd		%%ymm15,%%ymm3,%%ymm3					\n\t"\
		"vshufpd	$15,%%ymm2,%%ymm8,%%ymm0				\n\t		vshufpd	$15,%%ymm3,%%ymm9,%%ymm1					\n\t"\
		"vshufpd	$0 ,%%ymm2,%%ymm8,%%ymm8				\n\t		vshufpd	$0 ,%%ymm3,%%ymm9,%%ymm9					\n\t"\
		"vperm2f128 $32,%%ymm0,%%ymm6,%%ymm2	/* Re B	*/	\n\t		vperm2f128 $32,%%ymm1,%%ymm7,%%ymm3		/* Im B	*/	\n\t"\
		"vperm2f128 $49,%%ymm0,%%ymm6,%%ymm6	/* Re D	*/	\n\t		vperm2f128 $49,%%ymm1,%%ymm7,%%ymm7		/* Im D	*/	\n\t"\
		"vperm2f128 $32,%%ymm8,%%ymm4,%%ymm0	/* Re A	*/	\n\t		vperm2f128 $32,%%ymm9,%%ymm5,%%ymm1 	/* Im A	*/	\n\t"\
		"vperm2f128 $49,%%ymm8,%%ymm4,%%ymm4	/* Re C	*/	\n\t		vperm2f128 $49,%%ymm9,%%ymm5,%%ymm5 	/* Im C	*/	\n\t"\
		"subq		$%c[__odd_radix],%%rdx				\n\t"\
		/* Base negacyclic roots at this address in [0,2,4,6]*0x20 (Re parts), [1,3,5,7]*0x20 (Imag parts) */\
		"movq		%[__sign_mask],%%rax	\n\t"\
		"vmovaps	(%%rax),%%ymm15	\n\t"	/* ymm15 free for rest of way; use to store sign_mask needed for floating ABS */\
		"movq	%[__base_root] ,%%rax		\n\t"	/* Won't need main-array again until output transpose, so re-use rax for base_root */\
		/* half_arr[0,1,2,3] = [base*2, baseinv*2,wt_re*2,wt_im*2] */\
		"vmovaps	-0x40(%%rdx),%%ymm13	\n\t"	/* XMM13 = maxerr */\
		"addq		%%rdi,%%rdx				\n\t"	/* icycle assumed already in left-shifted ptr-byte-offset form */\
		"vmovaps	%%ymm13,%%ymm14			\n\t"	/* maxerr copy */\
		"movq		%[__cy_re],%%rbx		\n\t"\
	/* Do a-quartet: Data in ymm0,ymm1 */\
	"prefetcht0	(%%rcx)	\n\t"\
		/* For a-quartet, needed negacyclic root already in ymm10/11: */\
		/* Data in ymm0,ymm1 */\
		"vmovaps	%%ymm0,%%ymm8			\n\t"	/* x copy */\
		"vmovaps	%%ymm1,%%ymm9			\n\t"	/* y copy */\
		/* Inverse negacyclic weight is (wt_re, -wt_im): */\
		"vmulpd		%%ymm10,%%ymm0,%%ymm0	\n\t"	/* wt_re*[x     ] */\
		"vmulpd		%%ymm11,%%ymm9,%%ymm9	\n\t"	/* wt_im*[y copy] */\
		"vmulpd		%%ymm10,%%ymm1,%%ymm1	\n\t"	/* wt_re*[y     ] */\
		"vmulpd		%%ymm11,%%ymm8,%%ymm8	\n\t"	/* wt_im*[x copy] */\
		"vaddpd		%%ymm9 ,%%ymm0,%%ymm0	\n\t"	/* [a0-3.re] = x*wt_re + y*wt_im */\
		"vsubpd		%%ymm8 ,%%ymm1,%%ymm1	\n\t"	/* [a0-3.im] = y*wt_re - x*wt_im */\
		/* normalize a-pair, compute carryout, compute ROE: */\
		"vaddpd		           (%%rbx),%%ymm0,%%ymm0	\n\t"	/* [a0-3.re] + [cy0-3.re] */\
		"vaddpd		%c[__cy_im](%%rbx),%%ymm1,%%ymm1	\n\t"	/* [a0-3.im] + [cy0-3.im] */\
		"vmovaps	%%ymm0,%%ymm8			\n\t	vmovaps		%%ymm1,%%ymm9		\n\t"	/* copy x|y */\
		"vroundpd	$0,%%ymm0,%%ymm0		\n\t	vroundpd	$0,%%ymm1,%%ymm1	\n\t"	/* temp = DNINT(x|y) */\
		"vmovaps	%c[__odd_radm3](%%rdx),%%ymm12	\n\t"	/* [baseinv0-3] */\
		"vsubpd		%%ymm0 ,%%ymm8 ,%%ymm8 	\n\t	vsubpd		%%ymm1 ,%%ymm9 ,%%ymm9 	\n\t"	/* frac = [x - temp] */\
		"vandpd		%%ymm15,%%ymm8 ,%%ymm8 	\n\t	vandpd		%%ymm15,%%ymm9 ,%%ymm9 	\n\t"	/* frac = fabs(frac) */\
		"vmaxpd		%%ymm13,%%ymm8 ,%%ymm13	\n\t	vmaxpd		%%ymm14,%%ymm9 ,%%ymm14	\n\t"	/* if(frac > maxerr) maxerr=frac */\
		"vmovaps	%%ymm0 ,%%ymm8			\n\t	vmovaps		%%ymm1 ,%%ymm9			\n\t"	/* cpy temp */\
		"vmulpd		%%ymm12,%%ymm8 ,%%ymm8 	\n\t	vmulpd		%%ymm12,%%ymm9 ,%%ymm9 	\n\t"	/* temp*baseinv[0] ... inline the remaining +odd_radix offset in addressing */\
		"vmaxpd		%%ymm13,%%ymm14,%%ymm14	\n\t"	/* merge re|im maxerr vectors */\
		"vroundpd	$0,%%ymm8 ,%%ymm8 		\n\t	vroundpd	$0,%%ymm9,%%ymm9		\n\t"	/* [cy0-3.re] = DNINT(temp*baseinv[0]) */\
		"vmovaps	%c[__odd_radm2](%%rdx),%%ymm12	\n\t"	/* [base0-3] */\
		"vmovaps	%%ymm8,(%%rbx)			\n\t	vmovaps		%%ymm9,%c[__cy_im](%%rbx)\n\t"	/* store [cy0-3.re|im] */\
		"vmulpd		%%ymm12,%%ymm8,%%ymm8	\n\t	vmulpd		%%ymm12,%%ymm9,%%ymm9	\n\t"	/* base[0]*[cy0-3.re|im] */\
		"vsubpd		%%ymm8 ,%%ymm0,%%ymm0	\n\t	vsubpd		%%ymm9 ,%%ymm1,%%ymm1	\n\t"	/* XMM0|1 = [a0-3.re|im] = temp - [cy0-3.re|im]*base[0] */\
		"vmovaps	%%ymm0 ,%%ymm8			\n\t	vmovaps		%%ymm1 ,%%ymm9			\n\t"	/* cpy x|y */\
		/* Forward acyclic-convo weight is (wt_re, +wt_im): */\
		"vmulpd		%%ymm10,%%ymm0,%%ymm0	\n\t"	/* wt_re*[x     ] */\
		"vmulpd		%%ymm11,%%ymm9,%%ymm9	\n\t"	/* wt_im*[y copy] */\
		"vmulpd		%%ymm10,%%ymm1,%%ymm1	\n\t"	/* wt_re*[y     ] */\
		"vmulpd		%%ymm11,%%ymm8,%%ymm8	\n\t"	/* wt_im*[x copy] */\
		"vsubpd		%%ymm9 ,%%ymm0,%%ymm0	\n\t"	/* [a0-3.re] = x*wt_re - y*wt_im */\
		"vaddpd		%%ymm8 ,%%ymm1,%%ymm1	\n\t"	/* [a0-3.im] = y*wt_re + x*wt_im */\
	/* Now do b-quartet: Data in ymm2,ymm3 */\
	"movslq		%[__p1],%%r15	\n\t"	\
	"prefetcht0	(%%rcx,%%r15,8)	\n\t"\
		"movq		%[__half_arr],%%rdx		\n\t"\
		"vmovaps	%%ymm14,%%ymm13			\n\t"	/* maxerr copy */\
		"movslq		%[__jcycleA],%%rdi		\n\t"\
		"addq		%%rdi,%%rdx				\n\t"	/* jcycle assumed already in left-shifted ptr-byte-offset form */\
		/* Get next set of negacyclic roots: */\
		"vmovaps	0x40(%%rax),%%ymm10		\n\t"	/* c = Re part of 2nd base-root quartet */\
		"vmovaps	0x60(%%rax),%%ymm11		\n\t"	/* s = Im part */\
		/* Data in ymm2,ymm3 */\
		"vmovaps	%%ymm2,%%ymm8			\n\t"	/* x copy */\
		"vmovaps	%%ymm3,%%ymm9			\n\t"	/* y copy */\
		/* Inverse negacyclic weight is (wt_re, -wt_im): */\
		"vmulpd		%%ymm10,%%ymm2,%%ymm2	\n\t"	/* wt_re*[x     ] */\
		"vmulpd		%%ymm11,%%ymm9,%%ymm9	\n\t"	/* wt_im*[y copy] */\
		"vmulpd		%%ymm10,%%ymm3,%%ymm3	\n\t"	/* wt_re*[y     ] */\
		"vmulpd		%%ymm11,%%ymm8,%%ymm8	\n\t"	/* wt_im*[x copy] */\
		"vaddpd		%%ymm9 ,%%ymm2,%%ymm2	\n\t"	/* [a0-3.re] = x*wt_re + y*wt_im */\
		"vsubpd		%%ymm8 ,%%ymm3,%%ymm3	\n\t"	/* [a0-3.im] = y*wt_re - x*wt_im */\
		/* normalize, compute carryout, compute ROE: */\
		"vaddpd		           (%%rbx),%%ymm2,%%ymm2	\n\t"	/* [a0-3.re] + [cy0-3.re] */\
		"vaddpd		%c[__cy_im](%%rbx),%%ymm3,%%ymm3	\n\t"	/* [a0-3.im] + [cy0-3.im] */\
		"vmovaps	%%ymm2,%%ymm8			\n\t	vmovaps		%%ymm3,%%ymm9		\n\t"	/* copy x|y */\
		"vroundpd	$0,%%ymm2,%%ymm2		\n\t	vroundpd	$0,%%ymm3,%%ymm3	\n\t"	/* temp = DNINT(x|y) */\
		"vmovaps	%c[__odd_radm3](%%rdx),%%ymm12	\n\t"	/* [baseinv0-3] */\
		"vsubpd		%%ymm2 ,%%ymm8 ,%%ymm8 	\n\t	vsubpd		%%ymm3 ,%%ymm9 ,%%ymm9 	\n\t"	/* frac = [x - temp] */\
		"vandpd		%%ymm15,%%ymm8 ,%%ymm8 	\n\t	vandpd		%%ymm15,%%ymm9 ,%%ymm9 	\n\t"	/* frac = fabs(frac) */\
		"vmaxpd		%%ymm13,%%ymm8 ,%%ymm13	\n\t	vmaxpd		%%ymm14,%%ymm9 ,%%ymm14	\n\t"	/* if(frac > maxerr) maxerr=frac */\
		"vmovaps	%%ymm2 ,%%ymm8			\n\t	vmovaps		%%ymm3 ,%%ymm9			\n\t"	/* cpy temp */\
		"vmulpd		%%ymm12,%%ymm8 ,%%ymm8 	\n\t	vmulpd		%%ymm12,%%ymm9 ,%%ymm9 	\n\t"	/* temp*baseinv[0] ... inline the remaining +odd_radix offset in addressing */\
		"vmaxpd		%%ymm13,%%ymm14,%%ymm14	\n\t"	/* merge re|im maxerr vectors */\
		"vroundpd	$0,%%ymm8 ,%%ymm8 		\n\t	vroundpd	$0,%%ymm9,%%ymm9		\n\t"	/* [cy0-3.re] = DNINT(temp*baseinv[0]) */\
		"vmovaps	%c[__odd_radm2](%%rdx),%%ymm12	\n\t"	/* [base0-3] */\
		"vmovaps	%%ymm8,(%%rbx)			\n\t	vmovaps		%%ymm9,%c[__cy_im](%%rbx)\n\t"	/* store [cy0-3.re|im] */\
		"vmulpd		%%ymm12,%%ymm8,%%ymm8	\n\t	vmulpd		%%ymm12,%%ymm9,%%ymm9	\n\t"	/* base[0]*[cy0-3.re|im] */\
		"vsubpd		%%ymm8 ,%%ymm2,%%ymm2	\n\t	vsubpd		%%ymm9 ,%%ymm3,%%ymm3	\n\t"	/* XMM0|1 = [a0-3.re|im] = temp - [cy0-3.re|im]*base[0] */\
		"vmovaps	%%ymm2 ,%%ymm8			\n\t	vmovaps		%%ymm3 ,%%ymm9			\n\t"	/* cpy x|y */\
		/* Forward acyclic-convo weight is (wt_re, +wt_im): */\
		"vmulpd		%%ymm10,%%ymm2,%%ymm2	\n\t"	/* wt_re*[x     ] */\
		"vmulpd		%%ymm11,%%ymm9,%%ymm9	\n\t"	/* wt_im*[y copy] */\
		"vmulpd		%%ymm10,%%ymm3,%%ymm3	\n\t"	/* wt_re*[y     ] */\
		"vmulpd		%%ymm11,%%ymm8,%%ymm8	\n\t"	/* wt_im*[x copy] */\
		"vsubpd		%%ymm9 ,%%ymm2,%%ymm2	\n\t"	/* [a0-3.re] = x*wt_re - y*wt_im */\
		"vaddpd		%%ymm8 ,%%ymm3,%%ymm3	\n\t"	/* [a0-3.im] = y*wt_re + x*wt_im */\
	/* Now do c-quartet: Data in ymm4,ymm5 */\
	"movslq		%[__p2],%%r15	\n\t"	\
	"prefetcht0	(%%rcx,%%r15,8)	\n\t"\
		"movq		%[__half_arr],%%rdx		\n\t"\
		"vmovaps	%%ymm14,%%ymm13			\n\t"	/* maxerr copy */\
		"movslq		%[__kcycleA],%%rdi		\n\t"\
		"addq		%%rdi,%%rdx				\n\t"	/* kcycle assumed already in left-shifted ptr-byte-offset form */\
		/* Get next set of negacyclic roots: */\
		"vmovaps	0x80(%%rax),%%ymm10		\n\t"	/* c = Re part of 3rd base-root quartet */\
		"vmovaps	0xa0(%%rax),%%ymm11		\n\t"	/* s = Im part */\
		/* Data in ymm4,ymm5 */\
		"vmovaps	%%ymm4,%%ymm8			\n\t"	/* x copy */\
		"vmovaps	%%ymm5,%%ymm9			\n\t"	/* y copy */\
		/* Inverse negacyclic weight is (wt_re, -wt_im): */\
		"vmulpd		%%ymm10,%%ymm4,%%ymm4	\n\t"	/* wt_re*[x     ] */\
		"vmulpd		%%ymm11,%%ymm9,%%ymm9	\n\t"	/* wt_im*[y copy] */\
		"vmulpd		%%ymm10,%%ymm5,%%ymm5	\n\t"	/* wt_re*[y     ] */\
		"vmulpd		%%ymm11,%%ymm8,%%ymm8	\n\t"	/* wt_im*[x copy] */\
		"vaddpd		%%ymm9 ,%%ymm4,%%ymm4	\n\t"	/* [a0-3.re] = x*wt_re + y*wt_im */\
		"vsubpd		%%ymm8 ,%%ymm5,%%ymm5	\n\t"	/* [a0-3.im] = y*wt_re - x*wt_im */\
		/* normalize, compute carryout, compute ROE: */\
		"vaddpd		           (%%rbx),%%ymm4,%%ymm4	\n\t"	/* [a0-3.re] + [cy0-3.re] */\
		"vaddpd		%c[__cy_im](%%rbx),%%ymm5,%%ymm5	\n\t"	/* [a0-3.im] + [cy0-3.im] */\
		"vmovaps	%%ymm4,%%ymm8			\n\t	vmovaps		%%ymm5,%%ymm9		\n\t"	/* copy x|y */\
		"vroundpd	$0,%%ymm4,%%ymm4		\n\t	vroundpd	$0,%%ymm5,%%ymm5	\n\t"	/* temp = DNINT(x|y) */\
		"vmovaps	%c[__odd_radm3](%%rdx),%%ymm12	\n\t"	/* [baseinv0-3] */\
		"vsubpd		%%ymm4 ,%%ymm8 ,%%ymm8 	\n\t	vsubpd		%%ymm5 ,%%ymm9 ,%%ymm9 	\n\t"	/* frac = [x - temp] */\
		"vandpd		%%ymm15,%%ymm8 ,%%ymm8 	\n\t	vandpd		%%ymm15,%%ymm9 ,%%ymm9 	\n\t"	/* frac = fabs(frac) */\
		"vmaxpd		%%ymm13,%%ymm8 ,%%ymm13	\n\t	vmaxpd		%%ymm14,%%ymm9 ,%%ymm14	\n\t"	/* if(frac > maxerr) maxerr=frac */\
		"vmovaps	%%ymm4 ,%%ymm8			\n\t	vmovaps		%%ymm5 ,%%ymm9			\n\t"	/* cpy temp */\
		"vmulpd		%%ymm12,%%ymm8 ,%%ymm8 	\n\t	vmulpd		%%ymm12,%%ymm9 ,%%ymm9 	\n\t"	/* temp*baseinv[0] ... inline the remaining +odd_radix offset in addressing */\
		"vmaxpd		%%ymm13,%%ymm14,%%ymm14	\n\t"	/* merge re|im maxerr vectors */\
		"vroundpd	$0,%%ymm8 ,%%ymm8 		\n\t	vroundpd	$0,%%ymm9,%%ymm9		\n\t"	/* [cy0-3.re] = DNINT(temp*baseinv[0]) */\
		"vmovaps	%c[__odd_radm2](%%rdx),%%ymm12	\n\t"	/* [base0-3] */\
		"vmovaps	%%ymm8,(%%rbx)			\n\t	vmovaps		%%ymm9,%c[__cy_im](%%rbx)\n\t"	/* store [cy0-3.re|im] */\
		"vmulpd		%%ymm12,%%ymm8,%%ymm8	\n\t	vmulpd		%%ymm12,%%ymm9,%%ymm9	\n\t"	/* base[0]*[cy0-3.re|im] */\
		"vsubpd		%%ymm8 ,%%ymm4,%%ymm4	\n\t	vsubpd		%%ymm9 ,%%ymm5,%%ymm5	\n\t"	/* XMM0|1 = [a0-3.re|im] = temp - [cy0-3.re|im]*base[0] */\
		"vmovaps	%%ymm4 ,%%ymm8			\n\t	vmovaps		%%ymm5 ,%%ymm9			\n\t"	/* cpy x|y */\
		/* Forward acyclic-convo weight is (wt_re, +wt_im): */\
		"vmulpd		%%ymm10,%%ymm4,%%ymm4	\n\t"	/* wt_re*[x     ] */\
		"vmulpd		%%ymm11,%%ymm9,%%ymm9	\n\t"	/* wt_im*[y copy] */\
		"vmulpd		%%ymm10,%%ymm5,%%ymm5	\n\t"	/* wt_re*[y     ] */\
		"vmulpd		%%ymm11,%%ymm8,%%ymm8	\n\t"	/* wt_im*[x copy] */\
		"vsubpd		%%ymm9 ,%%ymm4,%%ymm4	\n\t"	/* [a0-3.re] = x*wt_re - y*wt_im */\
		"vaddpd		%%ymm8 ,%%ymm5,%%ymm5	\n\t"	/* [a0-3.im] = y*wt_re + x*wt_im */\
	/* Now do d-quartet: Data in ymm6,ymm7 */\
	"movslq		%[__p3],%%r15	\n\t"	\
	"prefetcht0	(%%rcx,%%r15,8)	\n\t"\
		"movq		%[__half_arr],%%rdx		\n\t"\
		"vmovaps	%%ymm14,%%ymm13			\n\t"	/* maxerr copy */\
		"movslq		%[__lcycleA],%%rdi		\n\t"\
		"addq		%%rdi,%%rdx				\n\t"	/* lcycle assumed already in left-shifted ptr-byte-offset form */\
		/* Get next set of negacyclic roots: */\
		"vmovaps	0xc0(%%rax),%%ymm10		\n\t"	/* c = Re part of 4th base-root quartet */\
		"vmovaps	0xe0(%%rax),%%ymm11		\n\t"	/* s = Im part */\
		/* Data in ymm6,ymm7 */\
		"vmovaps	%%ymm6,%%ymm8			\n\t"	/* x copy */\
		"vmovaps	%%ymm7,%%ymm9			\n\t"	/* y copy */\
		/* Inverse negacyclic weight is (wt_re, -wt_im): */\
		"vmulpd		%%ymm10,%%ymm6,%%ymm6	\n\t"	/* wt_re*[x     ] */\
		"vmulpd		%%ymm11,%%ymm9,%%ymm9	\n\t"	/* wt_im*[y copy] */\
		"vmulpd		%%ymm10,%%ymm7,%%ymm7	\n\t"	/* wt_re*[y     ] */\
		"vmulpd		%%ymm11,%%ymm8,%%ymm8	\n\t"	/* wt_im*[x copy] */\
		"vaddpd		%%ymm9 ,%%ymm6,%%ymm6	\n\t"	/* [a0-3.re] = x*wt_re + y*wt_im */\
		"vsubpd		%%ymm8 ,%%ymm7,%%ymm7	\n\t"	/* [a0-3.im] = y*wt_re - x*wt_im */\
		/* normalize, compute carryout, compute ROE: */\
		"vaddpd		           (%%rbx),%%ymm6,%%ymm6	\n\t"	/* [a0-3.re] + [cy0-3.re] */\
		"vaddpd		%c[__cy_im](%%rbx),%%ymm7,%%ymm7	\n\t"	/* [a0-3.im] + [cy0-3.im] */\
		"vmovaps	%%ymm6,%%ymm8			\n\t	vmovaps		%%ymm7,%%ymm9		\n\t"	/* copy x|y */\
		"vroundpd	$0,%%ymm6,%%ymm6		\n\t	vroundpd	$0,%%ymm7,%%ymm7	\n\t"	/* temp = DNINT(x|y) */\
		"vmovaps	%c[__odd_radm3](%%rdx),%%ymm12	\n\t"	/* [baseinv0-3] */\
		"vsubpd		%%ymm6 ,%%ymm8 ,%%ymm8 	\n\t	vsubpd		%%ymm7 ,%%ymm9 ,%%ymm9 	\n\t"	/* frac = [x - temp] */\
		"vandpd		%%ymm15,%%ymm8 ,%%ymm8 	\n\t	vandpd		%%ymm15,%%ymm9 ,%%ymm9 	\n\t"	/* frac = fabs(frac) */\
		"vmaxpd		%%ymm13,%%ymm8 ,%%ymm13	\n\t	vmaxpd		%%ymm14,%%ymm9 ,%%ymm14	\n\t"	/* if(frac > maxerr) maxerr=frac */\
		"vmovaps	%%ymm6 ,%%ymm8			\n\t	vmovaps		%%ymm7 ,%%ymm9			\n\t"	/* cpy temp */\
		"vmulpd		%%ymm12,%%ymm8 ,%%ymm8 	\n\t	vmulpd		%%ymm12,%%ymm9 ,%%ymm9 	\n\t"	/* temp*baseinv[0] ... inline the remaining +odd_radix offset in addressing */\
		"vmaxpd		%%ymm13,%%ymm14,%%ymm14	\n\t"	/* merge re|im maxerr vectors */\
		"vroundpd	$0,%%ymm8 ,%%ymm8 		\n\t	vroundpd	$0,%%ymm9,%%ymm9		\n\t"	/* [cy0-3.re] = DNINT(temp*baseinv[0]) */\
		"vmovaps	%c[__odd_radm2](%%rdx),%%ymm12	\n\t"	/* [base0-3] */\
		"vmovaps	%%ymm8,(%%rbx)			\n\t	vmovaps		%%ymm9,%c[__cy_im](%%rbx)\n\t"	/* store [cy0-3.re|im] */\
		"vmulpd		%%ymm12,%%ymm8,%%ymm8	\n\t	vmulpd		%%ymm12,%%ymm9,%%ymm9	\n\t"	/* base[0]*[cy0-3.re|im] */\
		"vsubpd		%%ymm8 ,%%ymm6,%%ymm6	\n\t	vsubpd		%%ymm9 ,%%ymm7,%%ymm7	\n\t"	/* XMM0|1 = [a0-3.re|im] = temp - [cy0-3.re|im]*base[0] */\
		"vmovaps	%%ymm6 ,%%ymm8			\n\t	vmovaps		%%ymm7 ,%%ymm9			\n\t"	/* cpy x|y */\
		/* Forward acyclic-convo weight is (wt_re, +wt_im): */\
		"vmulpd		%%ymm10,%%ymm6,%%ymm6	\n\t"	/* wt_re*[x     ] */\
		"vmulpd		%%ymm11,%%ymm9,%%ymm9	\n\t"	/* wt_im*[y copy] */\
		"vmulpd		%%ymm10,%%ymm7,%%ymm7	\n\t"	/* wt_re*[y     ] */\
		"vmulpd		%%ymm11,%%ymm8,%%ymm8	\n\t"	/* wt_im*[x copy] */\
		"vsubpd		%%ymm9 ,%%ymm6,%%ymm6	\n\t"	/* [a0-3.re] = x*wt_re - y*wt_im */\
		"vaddpd		%%ymm8 ,%%ymm7,%%ymm7	\n\t"	/* [a0-3.im] = y*wt_re + x*wt_im */\
		/* Store maxerr: */\
		"movq		%[__half_arr],%%rdx		\n\t"\
		"vmovaps	%%ymm14,-0x40(%%rdx)	\n\t"\
	/* 4-way transpose of outputs (Re, Im parts separately): Inputs from ymm0-7. Outputs into r0/1,2/3,4/5.6/7: */	\
	/* Because default inputs for our 4 x 4 transpose macro (e.g. the one used at start of this carry macro) */\
	/* are into ymm4/2/8/2, munge inputs into that order, resolving name-conflicts via use of the now-available ymm8-15 for outputs: */\
		"movq		%[__data],%%rax			\n\t"\
		"vshufpd	$15,%%ymm2,%%ymm0,%%ymm10					\n\t		vshufpd	$15,%%ymm3,%%ymm1,%%ymm11						\n\t"\
		"vshufpd	$0 ,%%ymm2,%%ymm0,%%ymm0					\n\t		vshufpd	$0 ,%%ymm3,%%ymm1,%%ymm1						\n\t"\
		"vshufpd	$15,%%ymm6,%%ymm4,%%ymm12					\n\t		vshufpd	$15,%%ymm7,%%ymm5,%%ymm13						\n\t"\
		"vshufpd	$0 ,%%ymm6,%%ymm4,%%ymm4					\n\t		vshufpd	$0 ,%%ymm7,%%ymm5,%%ymm5						\n\t"\
		"vperm2f128 $32,%%ymm12,%%ymm10,%%ymm2 		/* Re B	*/	\n\t		vperm2f128 $32,%%ymm13,%%ymm11,%%ymm3		/* Im B	*/	\n\t"\
		"vperm2f128 $49,%%ymm12,%%ymm10,%%ymm10		/* Re D	*/	\n\t		vperm2f128 $49,%%ymm13,%%ymm11,%%ymm11		/* Im D	*/	\n\t"\
		"vperm2f128 $32,%%ymm4 ,%%ymm0 ,%%ymm12		/* Re A	*/	\n\t		vperm2f128 $32,%%ymm5 ,%%ymm1 ,%%ymm13 		/* Im A	*/	\n\t"\
		"vperm2f128 $49,%%ymm4 ,%%ymm0 ,%%ymm0 		/* Re C	*/	\n\t		vperm2f128 $49,%%ymm5 ,%%ymm1 ,%%ymm1		/* Im C	*/	\n\t"\
		/* Multiply normalized, re-permuted transform outputs by forward IBDWT weights: */\
		"movslq		%[__icycleA],%%rdi		\n\t"	\
		"vmovaps	(%%rdx,%%rdi),%%ymm4	\n\t"	/* [wt0-3]A */\
		"vmovaps	(%%rdx,%%r9 ),%%ymm5	\n\t"	/* [wt0-3]B */\
		"vmovaps	(%%rdx,%%r8 ),%%ymm6	\n\t"	/* [wt0-3]C */\
		"vmovaps	(%%rdx,%%r10),%%ymm7	\n\t"	/* [wt0-3]D */\
		"vmulpd		%%ymm4,%%ymm12,%%ymm12						\n\t		vmulpd		%%ymm4,%%ymm13,%%ymm13			\n\t"\
		"vmulpd		%%ymm5,%%ymm2 ,%%ymm2 						\n\t		vmulpd		%%ymm5,%%ymm3 ,%%ymm3 			\n\t"\
		"vmulpd		%%ymm6,%%ymm0 ,%%ymm0 						\n\t		vmulpd		%%ymm6,%%ymm1 ,%%ymm1 			\n\t"\
		"vmulpd		%%ymm7,%%ymm10,%%ymm10						\n\t		vmulpd		%%ymm7,%%ymm11,%%ymm11			\n\t"\
		/* And write 'em back to memory: */\
		"vmovaps	%%ymm12,     (%%rax)						\n\t		vmovaps	%%ymm13,0x020(%%rax)				\n\t"\
		"vmovaps	%%ymm2 ,0x040(%%rax)						\n\t		vmovaps	%%ymm3 ,0x060(%%rax)				\n\t"\
		"vmovaps	%%ymm0 ,0x080(%%rax)						\n\t		vmovaps	%%ymm1 ,0x0a0(%%rax)				\n\t"\
		"vmovaps	%%ymm10,0x0c0(%%rax)						\n\t		vmovaps	%%ymm11,0x0e0(%%rax)				\n\t"\
		:						/* outputs: none */\
		:	[__data]		"m" (Xdata)	/* All inputs from memory addresses here */\
		,	[__base_root]	"m" (Xbase_root)\
		,	[__cmul_offset] "m" (Xcmul_offset)\
		,	[__cy_re]		"m" (Xcy_re)\
		,	[__cy_im]	"e" (Xcy_im)	/* Use literal-byte-offset for this ome to save a reg */\
		/* [1,2,3]*odd_radix are the index offsets to the wtinv, base, and base_inv values, resp. - assumed << l2_sz_vd on input: */\
		,	[__odd_radix]   "e" (Xodd_radix)\
		,	[__odd_radm2]   "e" (Xodd_radm2)\
		,	[__odd_radm3]   "e" (Xodd_radm3)\
		,	[__half_arr]	"m" (Xhalf_arr)\
		,	[__sign_mask]	"m" (Xsign_mask)\
		/* Need quartet of ascending [modulo odd_radix] icycle indices for IBDWT weights: */\
		,	[__icycleA]		"m" (XicycleA)\
		,	[__icycleB]		"m" (XicycleB)\
		,	[__icycleC]		"m" (XicycleC)\
		,	[__icycleD]		"m" (XicycleD)\
		/* Need quartet of same-index [i,j,k,l]cycle indices for negacyclic weights and base/baseinv normalizations: */\
		,	[__jcycleA]		"m" (XjcycleA)\
		,	[__kcycleA]		"m" (XkcycleA)\
		,	[__lcycleA]		"m" (XlcycleA)\
		/* Prefetch: base address and 3 index offsets */\
		,	[__add0] "m" (Xadd0)\
		,	[__p1] "m" (Xp1)\
		,	[__p2] "m" (Xp2)\
		,	[__p3] "m" (Xp3)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","r8","r9","r10","r15","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"   /* Clobbered registers */\
	);\
	}

	#define SSE2_fermat_carry_init_loacc(Xbase_root)\
	{\
	__asm__ volatile (\
		"movq		%[__base_root] ,%%rax	\n\t	"	/* Base negacyclic roots at this address +8*0x20 (Re parts), +9*0x20 (Imag parts) */\
		"vmovaps	0x100(%%rax),%%ymm10	\n\t	"	/* Multiply by exp(j*I*Pi/2)/RADIX, for j = 0-3 */\
		"vmovaps	0x120(%%rax),%%ymm11	\n\t	"	/* c = Re(exp) in ymm0, s = Im(exp) in ymm1 */\
		"vmovaps	     (%%rax),%%ymm0		\n\t	vmovaps	 0x40(%%rax),%%ymm2			\n\t"	/* x = Re part of 1st base-root quartet */\
		"vmovaps	 0x20(%%rax),%%ymm1		\n\t	vmovaps	 0x60(%%rax),%%ymm3			\n\t"	/* y = Im part */\
		"vmovaps	%%ymm0,%%ymm4			\n\t	vmovaps		%%ymm2,%%ymm6			\n\t"	/* Copy x */\
		"vmovaps	%%ymm1,%%ymm5			\n\t	vmovaps		%%ymm3,%%ymm7			\n\t"	/* Copy y */\
		"vmulpd		%%ymm10,%%ymm0,%%ymm0	\n\t	vmulpd		%%ymm10,%%ymm2,%%ymm2	\n\t"	/* c.x */\
		"vmulpd		%%ymm10,%%ymm1,%%ymm1	\n\t	vmulpd		%%ymm10,%%ymm3,%%ymm3	\n\t"	/* c.y */\
	"vfnmadd231pd	%%ymm11,%%ymm5,%%ymm0 	\n\t vfnmadd231pd	%%ymm11,%%ymm7,%%ymm2	\n\t"/* Out.re = c.x - s.y */\
	" vfmadd231pd	%%ymm11,%%ymm4,%%ymm1 	\n\t  vfmadd231pd	%%ymm11,%%ymm6,%%ymm3	\n\t"/* Out.im = c.y + s.x */\
		"vmovaps	%%ymm0 ,    (%%rax)		\n\t	vmovaps		%%ymm2 ,0x40(%%rax)		\n\t"	/* Store result, overwriting input base root */\
		"vmovaps	%%ymm1 ,0x20(%%rax)		\n\t	vmovaps		%%ymm3 ,0x60(%%rax)		\n\t"	/* Im part */\
		"/* Process next 2 base-root quartets: */"\
		"vmovaps	 0x80(%%rax),%%ymm0		\n\t	vmovaps	 0xc0(%%rax),%%ymm2			\n\t"	\
		"vmovaps	 0xa0(%%rax),%%ymm1		\n\t	vmovaps	 0xe0(%%rax),%%ymm3			\n\t"	\
		"vmovaps	%%ymm0,%%ymm4			\n\t	vmovaps		%%ymm2,%%ymm6			\n\t"	\
		"vmovaps	%%ymm1,%%ymm5			\n\t	vmovaps		%%ymm3,%%ymm7			\n\t"	\
		"vmulpd		%%ymm10,%%ymm0,%%ymm0	\n\t	vmulpd		%%ymm10,%%ymm2,%%ymm2	\n\t"	\
		"vmulpd		%%ymm10,%%ymm1,%%ymm1	\n\t	vmulpd		%%ymm10,%%ymm3,%%ymm3	\n\t"	\
	"vfnmadd231pd	%%ymm11,%%ymm5,%%ymm0 	\n\t vfnmadd231pd	%%ymm11,%%ymm7,%%ymm2	\n\t"	\
	" vfmadd231pd	%%ymm11,%%ymm4,%%ymm1 	\n\t  vfmadd231pd	%%ymm11,%%ymm6,%%ymm3	\n\t"	\
		"vmovaps	%%ymm0 ,0x80(%%rax)		\n\t	vmovaps		%%ymm2 ,0xc0(%%rax)		\n\t"	/* Store result, overwriting input base root */\
		"vmovaps	%%ymm1 ,0xa0(%%rax)		\n\t	vmovaps		%%ymm3 ,0xe0(%%rax)		\n\t"	/* Im part */\
	:						/* outputs: none */\
	:	[__base_root]	"m" (Xbase_root)	/* All inputs from memory addresses here */\
		: "cc","memory","rax","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm10","xmm11"   /* Clobbered registers */\
	);\
	}

	/* Non-power-of-2-runlength Fermat-mod acyclic-transform/IBDWT carry macro.
	The array indices i/j/k/lcycle declared int in caller but assumed to have been byte-shift-converted at time this macro called,
	thus can use as complex-address-offsets.  Use bytewise literal offsets to save registers for several args here,as vvv-marked:
												                             vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
	#define SSE2_fermat_carry_norm_errcheck_X4_loacc(Xdata,Xbase_root,Xcy_re,Xcy_im,Xodd_radix,Xodd_radm2,Xodd_radm3,Xhalf_arr,Xsign_mask,XicycleA,XicycleB,XicycleC,XicycleD, XjcycleA,XkcycleA,XlcycleA, Xadd0,Xp1,Xp2,Xp3)\
	{\
	__asm__ volatile (\
		"movq		%[__add0],%%rcx		\n\t"/* base address for 4 prefetches-from-main-data-array spread through this macro */\
		"movq		%[__data],%%rax		\n\t"\
		"movq		%[__half_arr],%%rdx		\n\t"\
		"addq		$%c[__odd_radix],%%rdx	\n\t"	/* wt|wtinv|base|baseinv data offset by icycle array slots from resp. base addresses */\
		/* Multiply complex transform outputs [x,y] = [re,im] by inverse IBDWT weights, which include the 2/n scale factor: */\
		"movslq		%[__icycleA],%%r15		\n\t"	\
		"movslq		%[__icycleB],%%r9 		\n\t"	\
		"movslq		%[__icycleC],%%r8 		\n\t"	\
		"movslq		%[__icycleD],%%r10		\n\t"	\
		"vmovaps	(%%rdx,%%r15),%%ymm10	\n\t"	/* [wtinv0-3]A */\
		"vmovaps	(%%rdx,%%r9 ),%%ymm11	\n\t"	/* [wtinv0-3]B */\
		"vmovaps	(%%rdx,%%r8 ),%%ymm12	\n\t"	/* [wtinv0-3]C */\
		"vmovaps	(%%rdx,%%r10),%%ymm13	\n\t"	/* [wtinv0-3]D */\
	/* 4-way transpose of inputs (Re, Im parts separately): Inputs from r0/1,2/3,4/5.6/7. Outputs into ymm0-7: */\
	/* Real parts use ymm0,2,4,6, ymm8 as tmp-reg:					Imag parts use ymm1,3,5,7, ymm9 as tm-reg: */\
		"vmovaps	     (%%rax),%%ymm4						\n\t		vmovaps	0x020(%%rax),%%ymm5							\n\t"\
		"vmovaps	0x040(%%rax),%%ymm2						\n\t		vmovaps	0x060(%%rax),%%ymm3							\n\t"\
		"vmulpd		%%ymm10,%%ymm4,%%ymm4					\n\t		vmulpd		%%ymm10,%%ymm5,%%ymm5					\n\t"\
		"vmulpd		%%ymm11,%%ymm2,%%ymm2					\n\t		vmulpd		%%ymm11,%%ymm3,%%ymm3					\n\t"\
		"vshufpd	$15,%%ymm2,%%ymm4,%%ymm6				\n\t		vshufpd	$15,%%ymm3,%%ymm5,%%ymm7					\n\t"\
		"vshufpd	$0 ,%%ymm2,%%ymm4,%%ymm4				\n\t		vshufpd	$0 ,%%ymm3,%%ymm5,%%ymm5					\n\t"\
		"vmovaps	0x080(%%rax),%%ymm8						\n\t		vmovaps	0x0a0(%%rax),%%ymm9							\n\t"\
		"vmovaps	0x0c0(%%rax),%%ymm2						\n\t		vmovaps	0x0e0(%%rax),%%ymm3							\n\t"\
		"vmulpd		%%ymm12,%%ymm8,%%ymm8					\n\t		vmulpd		%%ymm12,%%ymm9,%%ymm9					\n\t"\
		"vmulpd		%%ymm13,%%ymm2,%%ymm2					\n\t		vmulpd		%%ymm13,%%ymm3,%%ymm3					\n\t"\
		"vshufpd	$15,%%ymm2,%%ymm8,%%ymm0				\n\t		vshufpd	$15,%%ymm3,%%ymm9,%%ymm1					\n\t"\
		"vshufpd	$0 ,%%ymm2,%%ymm8,%%ymm8				\n\t		vshufpd	$0 ,%%ymm3,%%ymm9,%%ymm9					\n\t"\
		"vperm2f128 $32,%%ymm0,%%ymm6,%%ymm2	/* Re B	*/	\n\t		vperm2f128 $32,%%ymm1,%%ymm7,%%ymm3		/* Im B	*/	\n\t"\
		"vperm2f128 $49,%%ymm0,%%ymm6,%%ymm6	/* Re D	*/	\n\t		vperm2f128 $49,%%ymm1,%%ymm7,%%ymm7		/* Im D	*/	\n\t"\
		"vperm2f128 $32,%%ymm8,%%ymm4,%%ymm0	/* Re A	*/	\n\t		vperm2f128 $32,%%ymm9,%%ymm5,%%ymm1 	/* Im A	*/	\n\t"\
		"vperm2f128 $49,%%ymm8,%%ymm4,%%ymm4	/* Re C	*/	\n\t		vperm2f128 $49,%%ymm9,%%ymm5,%%ymm5 	/* Im C	*/	\n\t"\
		"subq		$%c[__odd_radix],%%rdx				\n\t"\
		/* Base negacyclic roots at this address in [0,2,4,6]*0x20 (Re parts), [1,3,5,7]*0x20 (Imag parts) */\
		"movq		%[__sign_mask],%%rax	\n\t"\
		"vmovaps	(%%rax),%%ymm15	\n\t"	/* ymm15 free for rest of way; use to store sign_mask needed for floating ABS */\
		"movq	%[__base_root] ,%%rax		\n\t"	/* Won't need main-array again until output transpose, so re-use rax for base_root */\
		/* half_arr[0,1,2,3] = [base*2, baseinv*2,wt_re*2,wt_im*2] */\
		"vmovaps	-0x40(%%rdx),%%ymm13	\n\t"	/* XMM13 = maxerr */\
		"addq		%%r15,%%rdx				\n\t"	/* icycle assumed already in left-shifted ptr-byte-offset form */\
		/* Get next set of negacyclic roots: */\
		"vmovaps	    (%%rax),%%ymm10		\n\t"	/* c = Re part of 1st base-root quartet */\
		"vmovaps	0x20(%%rax),%%ymm11		\n\t"	/* s = Im part */\
	/* Do a-quartet: Data in ymm0,ymm1 */\
	"prefetcht0	(%%rcx)	\n\t"\
		"vmovaps	%%ymm13,%%ymm14			\n\t"	/* maxerr copy */\
		"movq		%[__cy_re],%%rbx		\n\t"\
		"vmovaps	%%ymm0,%%ymm8			\n\t"	/* x copy */\
		"vmovaps	%%ymm1,%%ymm9			\n\t"	/* y copy */\
		/* Inverse negacyclic weight is (wt_re, -wt_im): */\
		"vmulpd		%%ymm10,%%ymm0,%%ymm0	\n\t"	/* wt_re*[x     ] */\
		"vmulpd		%%ymm10,%%ymm1,%%ymm1	\n\t"	/* wt_re*[y     ] */\
	" vfmadd231pd	%%ymm11,%%ymm9,%%ymm0 	\n\t"	/* wt_im*[y copy] ...[a0-3.re] = x*wt_re + y*wt_im */\
	"vfnmadd231pd	%%ymm11,%%ymm8,%%ymm1 	\n\t"	/* wt_im*[x copy] ...[a0-3.im] = y*wt_re - x*wt_im */\
		/* normalize a-pair, compute carryout, compute ROE: */\
		"vaddpd		           (%%rbx),%%ymm0,%%ymm0	\n\t"	/* [a0-3.re] + [cy0-3.re] */\
		"vaddpd		%c[__cy_im](%%rbx),%%ymm1,%%ymm1	\n\t"	/* [a0-3.im] + [cy0-3.im] */\
		"vmovaps	%%ymm0,%%ymm8			\n\t	vmovaps		%%ymm1,%%ymm9		\n\t"	/* copy x|y */\
		"vroundpd	$0,%%ymm0,%%ymm0		\n\t	vroundpd	$0,%%ymm1,%%ymm1	\n\t"	/* temp = DNINT(x|y) */\
		"vmovaps	%c[__odd_radm3](%%rdx),%%ymm12	\n\t"	/* [baseinv0-3] */\
		"vsubpd		%%ymm0 ,%%ymm8 ,%%ymm8 	\n\t	vsubpd		%%ymm1 ,%%ymm9 ,%%ymm9 	\n\t"	/* frac = [x - temp] */\
		"vandpd		%%ymm15,%%ymm8 ,%%ymm8 	\n\t	vandpd		%%ymm15,%%ymm9 ,%%ymm9 	\n\t"	/* frac = fabs(frac) */\
		"vmaxpd		%%ymm13,%%ymm8 ,%%ymm13	\n\t	vmaxpd		%%ymm14,%%ymm9 ,%%ymm14	\n\t"	/* if(frac > maxerr) maxerr=frac */\
		"vmovaps	%%ymm0 ,%%ymm8			\n\t	vmovaps		%%ymm1 ,%%ymm9			\n\t"	/* cpy temp */\
		"vmulpd		%%ymm12,%%ymm8 ,%%ymm8 	\n\t	vmulpd		%%ymm12,%%ymm9 ,%%ymm9 	\n\t"	/* temp*baseinv[0] ... inline the remaining +odd_radix offset in addressing */\
		"vmaxpd		%%ymm13,%%ymm14,%%ymm14	\n\t"	/* merge re|im maxerr vectors */\
		"vroundpd	$0,%%ymm8 ,%%ymm8 		\n\t	vroundpd	$0,%%ymm9,%%ymm9		\n\t"	/* [cy0-3.re] = DNINT(temp*baseinv[0]) */\
		"vmovaps	%c[__odd_radm2](%%rdx),%%ymm12	\n\t"	/* [base0-3] */\
		"vmovaps	%%ymm8,(%%rbx)			\n\t	vmovaps		%%ymm9,%c[__cy_im](%%rbx)\n\t"	/* store [cy0-3.re|im] */\
	"vfnmadd231pd	%%ymm12,%%ymm8,%%ymm0 	\n\t vfnmadd231pd	%%ymm12,%%ymm9,%%ymm1	\n\t"	/* base[0]*[cy0-3.re|im] ... XMM0|1 = [a0-3.re|im] = temp - [cy0-3.re|im]*base[0] */\
		"vmovaps	%%ymm0 ,%%ymm8			\n\t	vmovaps		%%ymm1 ,%%ymm9			\n\t"	/* cpy x|y */\
		/* Forward acyclic-convo weight is (wt_re, +wt_im): */\
		"vmulpd		%%ymm10,%%ymm0,%%ymm0	\n\t"	/* wt_re*[x     ] */\
		"vmulpd		%%ymm10,%%ymm1,%%ymm1	\n\t"	/* wt_re*[y     ] */\
	"vfnmadd231pd	%%ymm11,%%ymm9,%%ymm0 	\n\t"	/* wt_im*[y copy] ... [a0-3.re] = x*wt_re - y*wt_im */\
	" vfmadd231pd	%%ymm11,%%ymm8,%%ymm1 	\n\t"	/* wt_im*[x copy] ... [a0-3.im] = y*wt_re + x*wt_im */\
		/* Up-multiply negacyclic roots stored in ymm10,11 by exp(j*I*Pi/2)/RADIX, for j = 4 */\
		"vmovaps	0x140(%%rax),%%ymm8 	\n\t"	/* x = Re(exp) in ymm10 */\
		"vmovaps	0x160(%%rax),%%ymm9 	\n\t"	/* y = Im(exp) in ymm11 */\
		"vmulpd		%%ymm10,%%ymm9 ,%%ymm12	\n\t"	/* ymm12 = c.y */\
		"vmulpd		%%ymm10,%%ymm8 ,%%ymm10	\n\t"	/* ymm10 = c.x */\
	" vfmadd231pd	%%ymm11,%%ymm8,%%ymm12	\n\t"	/*  ymm8  = s.x ... ymm11 = wt.im = s.x + c.y; ymm8 ,4 free */\
	"vfnmadd231pd	%%ymm11,%%ymm9,%%ymm10	\n\t"	/*  ymm9  = s.y ... ymm10 = wt.re = c.x - s.y; ymm9  free */\
		"vmovaps	%%ymm12,0x20(%%rax)		\n\t"	/* Im part */\
		"vmovaps	%%ymm10,    (%%rax)		\n\t"	/* Store result, overwriting the old base root */\
		/* Get next set of negacyclic roots: */\
		"vmovaps	0x40(%%rax),%%ymm10		\n\t"	/* c = Re part of 2nd base-root quartet */\
		"vmovaps	0x60(%%rax),%%ymm11		\n\t"	/* s = Im part */\
	/* Now do b-quartet: Data in ymm2,ymm3 */\
	"movslq		%[__p1],%%r15	\n\t"	\
	"prefetcht0	(%%rcx,%%r15,8)	\n\t"\
		"movq		%[__half_arr],%%rdx		\n\t"\
		"vmovaps	%%ymm14,%%ymm13			\n\t"	/* maxerr copy */\
		"movslq		%[__jcycleA],%%r15		\n\t"\
		"addq		%%r15,%%rdx				\n\t"	/* jcycle assumed already in left-shifted ptr-byte-offset form */\
		"vmovaps	%%ymm2,%%ymm8			\n\t"	/* x copy */\
		"vmovaps	%%ymm3,%%ymm9			\n\t"	/* y copy */\
		/* Inverse negacyclic weight is (wt_re, -wt_im): */\
		"vmulpd		%%ymm10,%%ymm2,%%ymm2	\n\t"	/* wt_re*[x     ] */\
		"vmulpd		%%ymm10,%%ymm3,%%ymm3	\n\t"	/* wt_re*[y     ] */\
	" vfmadd231pd	%%ymm11,%%ymm9,%%ymm2 	\n\t"	/* wt_im*[y copy] ...[a0-3.re] = x*wt_re + y*wt_im */\
	"vfnmadd231pd	%%ymm11,%%ymm8,%%ymm3 	\n\t"	/* wt_im*[x copy] ...[a0-3.im] = y*wt_re - x*wt_im */\
		/* normalize, compute carryout, compute ROE: */\
		"vaddpd		           (%%rbx),%%ymm2,%%ymm2	\n\t"	/* [a0-3.re] + [cy0-3.re] */\
		"vaddpd		%c[__cy_im](%%rbx),%%ymm3,%%ymm3	\n\t"	/* [a0-3.im] + [cy0-3.im] */\
		"vmovaps	%%ymm2,%%ymm8			\n\t	vmovaps		%%ymm3,%%ymm9		\n\t"	/* copy x|y */\
		"vroundpd	$0,%%ymm2,%%ymm2		\n\t	vroundpd	$0,%%ymm3,%%ymm3	\n\t"	/* temp = DNINT(x|y) */\
		"vmovaps	%c[__odd_radm3](%%rdx),%%ymm12	\n\t"	/* [baseinv0-3] */\
		"vsubpd		%%ymm2 ,%%ymm8 ,%%ymm8 	\n\t	vsubpd		%%ymm3 ,%%ymm9 ,%%ymm9 	\n\t"	/* frac = [x - temp] */\
		"vandpd		%%ymm15,%%ymm8 ,%%ymm8 	\n\t	vandpd		%%ymm15,%%ymm9 ,%%ymm9 	\n\t"	/* frac = fabs(frac) */\
		"vmaxpd		%%ymm13,%%ymm8 ,%%ymm13	\n\t	vmaxpd		%%ymm14,%%ymm9 ,%%ymm14	\n\t"	/* if(frac > maxerr) maxerr=frac */\
		"vmovaps	%%ymm2 ,%%ymm8			\n\t	vmovaps		%%ymm3 ,%%ymm9			\n\t"	/* cpy temp */\
		"vmulpd		%%ymm12,%%ymm8 ,%%ymm8 	\n\t	vmulpd		%%ymm12,%%ymm9 ,%%ymm9 	\n\t"	/* temp*baseinv[0] ... inline the remaining +odd_radix offset in addressing */\
		"vmaxpd		%%ymm13,%%ymm14,%%ymm14	\n\t"	/* merge re|im maxerr vectors */\
		"vroundpd	$0,%%ymm8 ,%%ymm8 		\n\t	vroundpd	$0,%%ymm9,%%ymm9		\n\t"	/* [cy0-3.re] = DNINT(temp*baseinv[0]) */\
		"vmovaps	%c[__odd_radm2](%%rdx),%%ymm12	\n\t"	/* [base0-3] */\
		"vmovaps	%%ymm8,(%%rbx)			\n\t	vmovaps		%%ymm9,%c[__cy_im](%%rbx)\n\t"	/* store [cy0-3.re|im] */\
	"vfnmadd231pd	%%ymm12,%%ymm8,%%ymm2 	\n\t vfnmadd231pd	%%ymm12,%%ymm9,%%ymm3	\n\t"	/* base[0]*[cy0-3.re|im] ... XMM0|1 = [a0-3.re|im] = temp - [cy0-3.re|im]*base[0] */\
		"vmovaps	%%ymm2 ,%%ymm8			\n\t	vmovaps		%%ymm3 ,%%ymm9			\n\t"	/* cpy x|y */\
		/* Forward acyclic-convo weight is (wt_re, +wt_im): */\
		"vmulpd		%%ymm10,%%ymm2,%%ymm2	\n\t"	/* wt_re*[x     ] */\
		"vmulpd		%%ymm10,%%ymm3,%%ymm3	\n\t"	/* wt_re*[y     ] */\
	"vfnmadd231pd	%%ymm11,%%ymm9,%%ymm2 	\n\t"	/* wt_im*[y copy] ... [a0-3.re] = x*wt_re - y*wt_im */\
	" vfmadd231pd	%%ymm11,%%ymm8,%%ymm3 	\n\t"	/* wt_im*[x copy] ... [a0-3.im] = y*wt_re + x*wt_im */\
		/* Up-multiply negacyclic roots stored in ymm10,11 by exp(j*I*Pi/2)/RADIX, for j = 4 */\
		"vmovaps	0x140(%%rax),%%ymm8 	\n\t"	/* x = Re(exp) in ymm10 */\
		"vmovaps	0x160(%%rax),%%ymm9 	\n\t"	/* y = Im(exp) in ymm11 */\
		"vmulpd		%%ymm10,%%ymm9 ,%%ymm12	\n\t"	/* ymm12 = c.y */\
		"vmulpd		%%ymm10,%%ymm8 ,%%ymm10	\n\t"	/* ymm10 = c.x */\
	" vfmadd231pd	%%ymm11,%%ymm8,%%ymm12	\n\t"	/*  ymm8  = s.x ... ymm11 = wt.im = s.x + c.y; ymm8 ,4 free */\
	"vfnmadd231pd	%%ymm11,%%ymm9,%%ymm10	\n\t"	/*  ymm9  = s.y ... ymm10 = wt.re = c.x - s.y; ymm9  free */\
		"vmovaps	%%ymm12,0x60(%%rax)		\n\t"	/* Im part */\
		"vmovaps	%%ymm10,0x40(%%rax)		\n\t"	/* Store result, overwriting the old base root */\
		/* Get next set of negacyclic roots: */\
		"vmovaps	0x80(%%rax),%%ymm10		\n\t"	/* c = Re part of 3rd base-root quartet */\
		"vmovaps	0xa0(%%rax),%%ymm11		\n\t"	/* s = Im part */\
	/* Now do c-quartet: Data in ymm4,ymm5 */\
	"movslq		%[__p2],%%r15	\n\t"	\
	"prefetcht0	(%%rcx,%%r15,8)	\n\t"\
		"movq		%[__half_arr],%%rdx		\n\t"\
		"vmovaps	%%ymm14,%%ymm13			\n\t"	/* maxerr copy */\
		"movslq		%[__kcycleA],%%r15		\n\t"\
		"addq		%%r15,%%rdx				\n\t"	/* kcycle assumed already in left-shifted ptr-byte-offset form */\
		"vmovaps	%%ymm4,%%ymm8			\n\t"	/* x copy */\
		"vmovaps	%%ymm5,%%ymm9			\n\t"	/* y copy */\
		/* Inverse negacyclic weight is (wt_re, -wt_im): */\
		"vmulpd		%%ymm10,%%ymm4,%%ymm4	\n\t"	/* wt_re*[x     ] */\
		"vmulpd		%%ymm10,%%ymm5,%%ymm5	\n\t"	/* wt_re*[y     ] */\
	" vfmadd231pd	%%ymm11,%%ymm9,%%ymm4 	\n\t"	/* wt_im*[y copy] ...[a0-3.re] = x*wt_re + y*wt_im */\
	"vfnmadd231pd	%%ymm11,%%ymm8,%%ymm5 	\n\t"	/* wt_im*[x copy] ...[a0-3.im] = y*wt_re - x*wt_im */\
		/* normalize, compute carryout, compute ROE: */\
		"vaddpd		           (%%rbx),%%ymm4,%%ymm4	\n\t"	/* [a0-3.re] + [cy0-3.re] */\
		"vaddpd		%c[__cy_im](%%rbx),%%ymm5,%%ymm5	\n\t"	/* [a0-3.im] + [cy0-3.im] */\
		"vmovaps	%%ymm4,%%ymm8			\n\t	vmovaps		%%ymm5,%%ymm9		\n\t"	/* copy x|y */\
		"vroundpd	$0,%%ymm4,%%ymm4		\n\t	vroundpd	$0,%%ymm5,%%ymm5	\n\t"	/* temp = DNINT(x|y) */\
		"vmovaps	%c[__odd_radm3](%%rdx),%%ymm12	\n\t"	/* [baseinv0-3] */\
		"vsubpd		%%ymm4 ,%%ymm8 ,%%ymm8 	\n\t	vsubpd		%%ymm5 ,%%ymm9 ,%%ymm9 	\n\t"	/* frac = [x - temp] */\
		"vandpd		%%ymm15,%%ymm8 ,%%ymm8 	\n\t	vandpd		%%ymm15,%%ymm9 ,%%ymm9 	\n\t"	/* frac = fabs(frac) */\
		"vmaxpd		%%ymm13,%%ymm8 ,%%ymm13	\n\t	vmaxpd		%%ymm14,%%ymm9 ,%%ymm14	\n\t"	/* if(frac > maxerr) maxerr=frac */\
		"vmovaps	%%ymm4 ,%%ymm8			\n\t	vmovaps		%%ymm5 ,%%ymm9			\n\t"	/* cpy temp */\
		"vmulpd		%%ymm12,%%ymm8 ,%%ymm8 	\n\t	vmulpd		%%ymm12,%%ymm9 ,%%ymm9 	\n\t"	/* temp*baseinv[0] ... inline the remaining +odd_radix offset in addressing */\
		"vmaxpd		%%ymm13,%%ymm14,%%ymm14	\n\t"	/* merge re|im maxerr vectors */\
		"vroundpd	$0,%%ymm8 ,%%ymm8 		\n\t	vroundpd	$0,%%ymm9,%%ymm9		\n\t"	/* [cy0-3.re] = DNINT(temp*baseinv[0]) */\
		"vmovaps	%c[__odd_radm2](%%rdx),%%ymm12	\n\t"	/* [base0-3] */\
		"vmovaps	%%ymm8,(%%rbx)			\n\t	vmovaps		%%ymm9,%c[__cy_im](%%rbx)\n\t"	/* store [cy0-3.re|im] */\
	"vfnmadd231pd	%%ymm12,%%ymm8,%%ymm4 	\n\t vfnmadd231pd	%%ymm12,%%ymm9,%%ymm5	\n\t"	/* base[0]*[cy0-3.re|im] ... XMM0|1 = [a0-3.re|im] = temp - [cy0-3.re|im]*base[0] */\
		"vmovaps	%%ymm4 ,%%ymm8			\n\t	vmovaps		%%ymm5 ,%%ymm9			\n\t"	/* cpy x|y */\
		/* Forward acyclic-convo weight is (wt_re, +wt_im): */\
		"vmulpd		%%ymm10,%%ymm4,%%ymm4	\n\t"	/* wt_re*[x     ] */\
		"vmulpd		%%ymm10,%%ymm5,%%ymm5	\n\t"	/* wt_re*[y     ] */\
	"vfnmadd231pd	%%ymm11,%%ymm9,%%ymm4 	\n\t"	/* wt_im*[y copy] ... [a0-3.re] = x*wt_re - y*wt_im */\
	" vfmadd231pd	%%ymm11,%%ymm8,%%ymm5 	\n\t"	/* wt_im*[x copy] ... [a0-3.im] = y*wt_re + x*wt_im */\
		/* Up-multiply negacyclic roots stored in ymm10,11 by exp(j*I*Pi/2)/RADIX, for j = 4 */\
		"vmovaps	0x140(%%rax),%%ymm8 	\n\t"	/* x = Re(exp) in ymm10 */\
		"vmovaps	0x160(%%rax),%%ymm9 	\n\t"	/* y = Im(exp) in ymm11 */\
		"vmulpd		%%ymm10,%%ymm9 ,%%ymm12	\n\t"	/* ymm12 = c.y */\
		"vmulpd		%%ymm10,%%ymm8 ,%%ymm10	\n\t"	/* ymm10 = c.x */\
	" vfmadd231pd	%%ymm11,%%ymm8,%%ymm12	\n\t"	/*  ymm8  = s.x ... ymm11 = wt.im = s.x + c.y; ymm8 ,4 free */\
	"vfnmadd231pd	%%ymm11,%%ymm9,%%ymm10	\n\t"	/*  ymm9  = s.y ... ymm10 = wt.re = c.x - s.y; ymm9  free */\
		"vmovaps	%%ymm12,0xa0(%%rax)		\n\t"	/* Im part */\
		"vmovaps	%%ymm10,0x80(%%rax)		\n\t"	/* Store result, overwriting the old base root */\
		/* Get next set of negacyclic roots: */\
		"vmovaps	0xc0(%%rax),%%ymm10		\n\t"	/* c = Re part of 4th base-root quartet */\
		"vmovaps	0xe0(%%rax),%%ymm11		\n\t"	/* s = Im part */\
	/* Now do d-quartet: Data in ymm6,ymm7 */\
	"movslq		%[__p3],%%r15	\n\t"	\
	"prefetcht0	(%%rcx,%%r15,8)	\n\t"\
		"movq		%[__half_arr],%%rdx		\n\t"\
		"vmovaps	%%ymm14,%%ymm13			\n\t"	/* maxerr copy */\
		"movslq		%[__lcycleA],%%r15		\n\t"\
		"addq		%%r15,%%rdx				\n\t"	/* lcycle assumed already in left-shifted ptr-byte-offset form */\
		"vmovaps	%%ymm6,%%ymm8			\n\t"	/* x copy */\
		"vmovaps	%%ymm7,%%ymm9			\n\t"	/* y copy */\
		/* Inverse negacyclic weight is (wt_re, -wt_im): */\
		"vmulpd		%%ymm10,%%ymm6,%%ymm6	\n\t"	/* wt_re*[x     ] */\
		"vmulpd		%%ymm10,%%ymm7,%%ymm7	\n\t"	/* wt_re*[y     ] */\
	" vfmadd231pd	%%ymm11,%%ymm9,%%ymm6 	\n\t"	/* wt_im*[y copy] ...[a0-3.re] = x*wt_re + y*wt_im */\
	"vfnmadd231pd	%%ymm11,%%ymm8,%%ymm7 	\n\t"	/* wt_im*[x copy] ...[a0-3.im] = y*wt_re - x*wt_im */\
		/* normalize, compute carryout, compute ROE: */\
		"vaddpd		           (%%rbx),%%ymm6,%%ymm6	\n\t"	/* [a0-3.re] + [cy0-3.re] */\
		"vaddpd		%c[__cy_im](%%rbx),%%ymm7,%%ymm7	\n\t"	/* [a0-3.im] + [cy0-3.im] */\
		"vmovaps	%%ymm6,%%ymm8			\n\t	vmovaps		%%ymm7,%%ymm9		\n\t"	/* copy x|y */\
		"vroundpd	$0,%%ymm6,%%ymm6		\n\t	vroundpd	$0,%%ymm7,%%ymm7	\n\t"	/* temp = DNINT(x|y) */\
		"vmovaps	%c[__odd_radm3](%%rdx),%%ymm12	\n\t"	/* [baseinv0-3] */\
		"vsubpd		%%ymm6 ,%%ymm8 ,%%ymm8 	\n\t	vsubpd		%%ymm7 ,%%ymm9 ,%%ymm9 	\n\t"	/* frac = [x - temp] */\
		"vandpd		%%ymm15,%%ymm8 ,%%ymm8 	\n\t	vandpd		%%ymm15,%%ymm9 ,%%ymm9 	\n\t"	/* frac = fabs(frac) */\
		"vmaxpd		%%ymm13,%%ymm8 ,%%ymm13	\n\t	vmaxpd		%%ymm14,%%ymm9 ,%%ymm14	\n\t"	/* if(frac > maxerr) maxerr=frac */\
		"vmovaps	%%ymm6 ,%%ymm8			\n\t	vmovaps		%%ymm7 ,%%ymm9			\n\t"	/* cpy temp */\
		"vmulpd		%%ymm12,%%ymm8 ,%%ymm8 	\n\t	vmulpd		%%ymm12,%%ymm9 ,%%ymm9 	\n\t"	/* temp*baseinv[0] ... inline the remaining +odd_radix offset in addressing */\
		"vmaxpd		%%ymm13,%%ymm14,%%ymm14	\n\t"	/* merge re|im maxerr vectors */\
		"vroundpd	$0,%%ymm8 ,%%ymm8 		\n\t	vroundpd	$0,%%ymm9,%%ymm9		\n\t"	/* [cy0-3.re] = DNINT(temp*baseinv[0]) */\
		"vmovaps	%c[__odd_radm2](%%rdx),%%ymm12	\n\t"	/* [base0-3] */\
		"vmovaps	%%ymm8,(%%rbx)			\n\t	vmovaps		%%ymm9,%c[__cy_im](%%rbx)\n\t"	/* store [cy0-3.re|im] */\
	"vfnmadd231pd	%%ymm12,%%ymm8,%%ymm6 	\n\t vfnmadd231pd	%%ymm12,%%ymm9,%%ymm7	\n\t"	/* base[0]*[cy0-3.re|im] ... XMM0|1 = [a0-3.re|im] = temp - [cy0-3.re|im]*base[0] */\
		"vmovaps	%%ymm6 ,%%ymm8			\n\t	vmovaps		%%ymm7 ,%%ymm9			\n\t"	/* cpy x|y */\
		/* Forward acyclic-convo weight is (wt_re, +wt_im): */\
		"vmulpd		%%ymm10,%%ymm6,%%ymm6	\n\t"	/* wt_re*[x     ] */\
		"vmulpd		%%ymm10,%%ymm7,%%ymm7	\n\t"	/* wt_re*[y     ] */\
	"vfnmadd231pd	%%ymm11,%%ymm9,%%ymm6 	\n\t"	/* wt_im*[y copy] ... [a0-3.re] = x*wt_re - y*wt_im */\
	" vfmadd231pd	%%ymm11,%%ymm8,%%ymm7 	\n\t"	/* wt_im*[x copy] ... [a0-3.im] = y*wt_re + x*wt_im */\
		/* Up-multiply negacyclic roots stored in ymm10,11 by exp(j*I*Pi/2)/RADIX, for j = 4 */\
		"vmovaps	0x140(%%rax),%%ymm8 	\n\t"	/* x = Re(exp) in ymm10 */\
		"vmovaps	0x160(%%rax),%%ymm9 	\n\t"	/* y = Im(exp) in ymm11 */\
		"vmulpd		%%ymm10,%%ymm9 ,%%ymm12	\n\t"	/* ymm12 = c.y */\
		"vmulpd		%%ymm10,%%ymm8 ,%%ymm10	\n\t"	/* ymm10 = c.x */\
	" vfmadd231pd	%%ymm11,%%ymm8,%%ymm12	\n\t"	/*  ymm8  = s.x ... ymm11 = wt.im = s.x + c.y; ymm8 ,4 free */\
	"vfnmadd231pd	%%ymm11,%%ymm9,%%ymm10	\n\t"	/*  ymm9  = s.y ... ymm10 = wt.re = c.x - s.y; ymm9  free */\
		"vmovaps	%%ymm12,0xe0(%%rax)		\n\t"	/* Im part */\
		"vmovaps	%%ymm10,0xc0(%%rax)		\n\t"	/* Store result, overwriting the old base root */\
		/* Store maxerr: */\
		"movq		%[__half_arr],%%rdx		\n\t"\
		"vmovaps	%%ymm14,-0x40(%%rdx)	\n\t"\
	/* 4-way transpose of outputs (Re, Im parts separately): Inputs from ymm0-7. Outputs into r0/1,2/3,4/5.6/7: */	\
	/* Because default inputs for our 4 x 4 transpose macro (e.g. the one used at start of this carry macro) */\
	/* are into ymm4/2/8/2, munge inputs into that order, resolving name-conflicts via use of the now-available ymm8-15 for outputs: */\
		"movq		%[__data],%%rax			\n\t"\
		"vshufpd	$15,%%ymm2,%%ymm0,%%ymm10					\n\t		vshufpd	$15,%%ymm3,%%ymm1,%%ymm11						\n\t"\
		"vshufpd	$0 ,%%ymm2,%%ymm0,%%ymm0					\n\t		vshufpd	$0 ,%%ymm3,%%ymm1,%%ymm1						\n\t"\
		"vshufpd	$15,%%ymm6,%%ymm4,%%ymm12					\n\t		vshufpd	$15,%%ymm7,%%ymm5,%%ymm13						\n\t"\
		"vshufpd	$0 ,%%ymm6,%%ymm4,%%ymm4					\n\t		vshufpd	$0 ,%%ymm7,%%ymm5,%%ymm5						\n\t"\
		"vperm2f128 $32,%%ymm12,%%ymm10,%%ymm2 		/* Re B	*/	\n\t		vperm2f128 $32,%%ymm13,%%ymm11,%%ymm3		/* Im B	*/	\n\t"\
		"vperm2f128 $49,%%ymm12,%%ymm10,%%ymm10		/* Re D	*/	\n\t		vperm2f128 $49,%%ymm13,%%ymm11,%%ymm11		/* Im D	*/	\n\t"\
		"vperm2f128 $32,%%ymm4 ,%%ymm0 ,%%ymm12		/* Re A	*/	\n\t		vperm2f128 $32,%%ymm5 ,%%ymm1 ,%%ymm13 		/* Im A	*/	\n\t"\
		"vperm2f128 $49,%%ymm4 ,%%ymm0 ,%%ymm0 		/* Re C	*/	\n\t		vperm2f128 $49,%%ymm5 ,%%ymm1 ,%%ymm1		/* Im C	*/	\n\t"\
		/* Multiply normalized, re-permuted transform outputs by forward IBDWT weights: */\
		"movslq		%[__icycleA],%%r15		\n\t"	\
		"vmovaps	(%%rdx,%%r15),%%ymm4	\n\t"	/* [wt0-3]A */\
		"vmovaps	(%%rdx,%%r9 ),%%ymm5	\n\t"	/* [wt0-3]B */\
		"vmovaps	(%%rdx,%%r8 ),%%ymm6	\n\t"	/* [wt0-3]C */\
		"vmovaps	(%%rdx,%%r10),%%ymm7	\n\t"	/* [wt0-3]D */\
		"vmulpd		%%ymm4,%%ymm12,%%ymm12						\n\t		vmulpd		%%ymm4,%%ymm13,%%ymm13			\n\t"\
		"vmulpd		%%ymm5,%%ymm2 ,%%ymm2 						\n\t		vmulpd		%%ymm5,%%ymm3 ,%%ymm3 			\n\t"\
		"vmulpd		%%ymm6,%%ymm0 ,%%ymm0 						\n\t		vmulpd		%%ymm6,%%ymm1 ,%%ymm1 			\n\t"\
		"vmulpd		%%ymm7,%%ymm10,%%ymm10						\n\t		vmulpd		%%ymm7,%%ymm11,%%ymm11			\n\t"\
		/* And write 'em back to memory: */\
		"vmovaps	%%ymm12,     (%%rax)						\n\t		vmovaps	%%ymm13,0x020(%%rax)				\n\t"\
		"vmovaps	%%ymm2 ,0x040(%%rax)						\n\t		vmovaps	%%ymm3 ,0x060(%%rax)				\n\t"\
		"vmovaps	%%ymm0 ,0x080(%%rax)						\n\t		vmovaps	%%ymm1 ,0x0a0(%%rax)				\n\t"\
		"vmovaps	%%ymm10,0x0c0(%%rax)						\n\t		vmovaps	%%ymm11,0x0e0(%%rax)				\n\t"\
		:						/* outputs: none */\
		:	[__data]		"m" (Xdata)	/* All inputs from memory addresses here */\
		,	[__base_root]	"m" (Xbase_root)\
		,	[__cy_re]		"m" (Xcy_re)\
		,	[__cy_im]		"e" (Xcy_im)	/* Use literal-byte-offset for this ome to save a reg */\
		/* [1,2,3]*odd_radix are the index offsets to the wtinv, base, and base_inv values, resp. - assumed << l2_sz_vd on input: */\
		,	[__odd_radix]   "e" (Xodd_radix)\
		,	[__odd_radm2]   "e" (Xodd_radm2)\
		,	[__odd_radm3]   "e" (Xodd_radm3)\
		,	[__half_arr]	"m" (Xhalf_arr)\
		,	[__sign_mask]	"m" (Xsign_mask)\
		/* Need quartet of ascending [modulo odd_radix] icycle indices for IBDWT weights: */\
		,	[__icycleA]		"m" (XicycleA)\
		,	[__icycleB]		"m" (XicycleB)\
		,	[__icycleC]		"m" (XicycleC)\
		,	[__icycleD]		"m" (XicycleD)\
		/* Need quartet of same-index [i,j,k,l]cycle indices for negacyclic weights and base/baseinv normalizations: */\
		,	[__jcycleA]		"m" (XjcycleA)\
		,	[__kcycleA]		"m" (XkcycleA)\
		,	[__lcycleA]		"m" (XlcycleA)\
		/* Prefetch: base address and 3 index offsets */\
		,	[__add0] "m" (Xadd0)\
		,	[__p1] "m" (Xp1)\
		,	[__p2] "m" (Xp2)\
		,	[__p3] "m" (Xp3)\
		: "cc","memory","rax","rbx","rcx","rdx","r8","r9","r10","r15","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14"   /* Clobbered registers */\
	);\
	}

  #else	// USE_AVX, no AVX2/FMA3 support:

	/* Power-of-2-runlength Fermat-mod acyclic-transform/IBDWT carry macro.

	NOTE: The array indices i/j/k/lcycle declared int in caller but assumed to have been
	byte-shift-converted at time this macro called, thus can use as complex-address-offsets.
	*/
	#define SSE2_fermat_carry_norm_pow2_errcheck_X4(Xdata,Xbase_root,Xcmul_offset,Xcy_re,Xcy_im,Xhalf_arr,Xsign_mask, Xadd0,Xp1,Xp2,Xp3)\
	{\
	__asm__ volatile (\
		"movq		%[__add0],%%r14		\n\t"/* base address for 4 prefetches-from-main-data-array spread through this macro */\
		/* Base negacyclic roots at this address in [0,2,4,6]*0x20 (Re parts), [1,3,5,7]*0x20 (Imag parts) */\
		"movq	%[__base_root] ,%%rax			\n\t"	/* Won't need main-array again until output transpose, so re-use rax for base_root */\
		"movslq	%[__cmul_offset],%%rbx			\n\t"\
		"addq	%%rax,%%rbx	\n\t"	/* Index into complex const multipliers block, each applied to 4 sets of base roots */\
		/* Up-multiply quartet of negacyclic roots used in this macro invocation; store sets 2-4 back into mem, keep set 1 in ymm10,11 [that's why we do sets 1/2 after 3/4] */\
		"vmovaps	    (%%rbx),%%ymm10		\n\t	"	/* Multiply by exp(j*I*Pi/2)/RADIX, for j = 0-3 */\
		"vmovaps	0x20(%%rbx),%%ymm11		\n\t	"	/* c = Re(exp) in ymm0, s = Im(exp) in ymm1 */\
		"/* Sets 3/4: */"\
		"vmovaps	 0x80(%%rax),%%ymm0		\n\t	vmovaps	 0xc0(%%rax),%%ymm2			\n\t"	\
		"vmovaps	 0xa0(%%rax),%%ymm1		\n\t	vmovaps	 0xe0(%%rax),%%ymm3			\n\t"	\
		"vmovaps	%%ymm0,%%ymm4			\n\t	vmovaps		%%ymm2,%%ymm6			\n\t"	\
		"vmovaps	%%ymm1,%%ymm5			\n\t	vmovaps		%%ymm3,%%ymm7			\n\t"	\
		"vmulpd		%%ymm10,%%ymm0,%%ymm0	\n\t	vmulpd		%%ymm10,%%ymm2,%%ymm2	\n\t"	\
		"vmulpd		%%ymm11,%%ymm5,%%ymm5	\n\t	vmulpd		%%ymm11,%%ymm7,%%ymm7	\n\t"	\
		"vmulpd		%%ymm10,%%ymm1,%%ymm1	\n\t	vmulpd		%%ymm10,%%ymm3,%%ymm3	\n\t"	\
		"vmulpd		%%ymm11,%%ymm4,%%ymm4	\n\t	vmulpd		%%ymm11,%%ymm6,%%ymm6	\n\t"	\
		"vsubpd		%%ymm5 ,%%ymm0,%%ymm0	\n\t	vsubpd		%%ymm7 ,%%ymm2,%%ymm2	\n\t"	\
		"vaddpd		%%ymm4 ,%%ymm1,%%ymm1	\n\t	vaddpd		%%ymm6 ,%%ymm3,%%ymm3	\n\t"	\
		"vmovaps	%%ymm0 ,0x80(%%rax)		\n\t	vmovaps		%%ymm2 ,0xc0(%%rax)		\n\t"	/* Store result, overwriting input base root */\
		"vmovaps	%%ymm1 ,0xa0(%%rax)		\n\t	vmovaps		%%ymm3 ,0xe0(%%rax)		\n\t"	/* Im part */\
		"/* Sets 1/2: */"\
		"vmovaps	     (%%rax),%%ymm0		\n\t	vmovaps	 0x40(%%rax),%%ymm2			\n\t"	/* x = Re part of 1st base-root quartet */\
		"vmovaps	 0x20(%%rax),%%ymm1		\n\t	vmovaps	 0x60(%%rax),%%ymm3			\n\t"	/* y = Im part */\
		"vmovaps	%%ymm0,%%ymm4			\n\t	vmovaps		%%ymm2,%%ymm6			\n\t"	/* Copy x */\
		"vmovaps	%%ymm1,%%ymm5			\n\t	vmovaps		%%ymm3,%%ymm7			\n\t"	/* Copy y */\
		"vmulpd		%%ymm10,%%ymm0,%%ymm0	\n\t	vmulpd		%%ymm10,%%ymm2,%%ymm2	\n\t"	/* c.x */\
		"vmulpd		%%ymm11,%%ymm5,%%ymm5	\n\t	vmulpd		%%ymm11,%%ymm7,%%ymm7	\n\t"	/* s.y */\
		"vmulpd		%%ymm10,%%ymm1,%%ymm1	\n\t	vmulpd		%%ymm10,%%ymm3,%%ymm3	\n\t"	/* c.y */\
		"vmulpd		%%ymm11,%%ymm4,%%ymm4	\n\t	vmulpd		%%ymm11,%%ymm6,%%ymm6	\n\t"	/* s.x */\
		"vsubpd		%%ymm5 ,%%ymm0,%%ymm10	\n\t	vsubpd		%%ymm7 ,%%ymm2,%%ymm2	\n\t"	/* Out.re = c.x - s.y */\
		"vaddpd		%%ymm4 ,%%ymm1,%%ymm11	\n\t	vaddpd		%%ymm6 ,%%ymm3,%%ymm3	\n\t"	/* Out.im = c.y + s.x */\
		"											vmovaps		%%ymm2 ,0x40(%%rax)		\n\t"	/* Store result, overwriting input base root */\
		"											vmovaps		%%ymm3 ,0x60(%%rax)		\n\t"	/* Im part */\
		"movq		%[__data],%%rax		\n\t"\
		"movq		%[__half_arr],%%rdx		\n\t"\
		"vmovaps	0x40(%%rdx),%%ymm12	\n\t"	/* XMM12 = scale */\
	/* 4-way transpose of inputs (Re, Im parts separately): Inputs from r0/1,2/3,4/5.6/7. Outputs into ymm0-7: */\
	/* Real parts use ymm0,2,4,6, ymm8 as tmp-reg:					Imag parts use ymm1,3,5,7, ymm9 as tm-reg: */\
		"vmovaps	     (%%rax),%%ymm4						\n\t		vmovaps	0x020(%%rax),%%ymm5							\n\t"\
		"vmovaps	0x040(%%rax),%%ymm2						\n\t		vmovaps	0x060(%%rax),%%ymm3							\n\t"\
		"vshufpd	$15,%%ymm2,%%ymm4,%%ymm6				\n\t		vshufpd	$15,%%ymm3,%%ymm5,%%ymm7					\n\t"\
		"vshufpd	$0 ,%%ymm2,%%ymm4,%%ymm4				\n\t		vshufpd	$0 ,%%ymm3,%%ymm5,%%ymm5					\n\t"\
		"vmovaps	0x080(%%rax),%%ymm8						\n\t		vmovaps	0x0a0(%%rax),%%ymm9							\n\t"\
		"vmovaps	0x0c0(%%rax),%%ymm2						\n\t		vmovaps	0x0e0(%%rax),%%ymm3							\n\t"\
		"vshufpd	$15,%%ymm2,%%ymm8,%%ymm0				\n\t		vshufpd	$15,%%ymm3,%%ymm9,%%ymm1					\n\t"\
		"vshufpd	$0 ,%%ymm2,%%ymm8,%%ymm8				\n\t		vshufpd	$0 ,%%ymm3,%%ymm9,%%ymm9					\n\t"\
		"vperm2f128 $32,%%ymm0,%%ymm6,%%ymm2	/* Re B	*/	\n\t		vperm2f128 $32,%%ymm1,%%ymm7,%%ymm3		/* Im B	*/	\n\t"\
		"vperm2f128 $49,%%ymm0,%%ymm6,%%ymm6	/* Re D	*/	\n\t		vperm2f128 $49,%%ymm1,%%ymm7,%%ymm7		/* Im D	*/	\n\t"\
		"vperm2f128 $32,%%ymm8,%%ymm4,%%ymm0	/* Re A	*/	\n\t		vperm2f128 $32,%%ymm9,%%ymm5,%%ymm1 	/* Im A	*/	\n\t"\
		"vperm2f128 $49,%%ymm8,%%ymm4,%%ymm4	/* Re C	*/	\n\t		vperm2f128 $49,%%ymm9,%%ymm5,%%ymm5 	/* Im C	*/	\n\t"\
	/* Apply inverse-complex-runlength scaling factor to the data: */\
		"vmulpd		%%ymm12,%%ymm2,%%ymm2					\n\t		vmulpd		%%ymm12,%%ymm3,%%ymm3	\n\t"\
		"vmulpd		%%ymm12,%%ymm6,%%ymm6					\n\t		vmulpd		%%ymm12,%%ymm7,%%ymm7	\n\t"\
		"vmulpd		%%ymm12,%%ymm0,%%ymm0					\n\t		vmulpd		%%ymm12,%%ymm1,%%ymm1	\n\t"\
		"vmulpd		%%ymm12,%%ymm4,%%ymm4					\n\t		vmulpd		%%ymm12,%%ymm5,%%ymm5	\n\t"\
		/* Base negacyclic roots at this address in [0,2,4,6]*0x20 (Re parts), [1,3,5,7]*0x20 (Imag parts) */\
		"movq	%[__base_root] ,%%rax		\n\t"	/* Won't need main-array again until output transpose, so re-use rax for base_root */\
		/* half_arr[0,1,2,3] = [base*2, baseinv*2,wt_re*2,wt_im*2] */\
		"vmovaps	-0x40(%%rdx),%%ymm13	\n\t"	/* XMM13 = maxerr */\
		"/*vmovaps	-0x20(%%rdx),%%ymm15	*/\n\t"	/* rnd_const; prefer ROUNDPD in AVX mode, so ymm15 free */\
		"vmovaps	%%ymm13,%%ymm14			\n\t"	/* maxerr copy */\
		"movq		%[__cy_re],%%rbx		\n\t"\
		"movq		%[__cy_im],%%rcx		\n\t"\
	/* Do a-quartet: Data in ymm0,ymm1 */\
	"prefetcht0	(%%r14)	\n\t"\
		/* For a-quartet, needed negacyclic root already in ymm10/11: */\
		/* Data in ymm0,ymm1 */\
		"vmovaps	%%ymm0,%%ymm8			\n\t"	/* x copy */\
		"vmovaps	%%ymm1,%%ymm9			\n\t"	/* y copy */\
		/* Inverse negacyclic weight is (wt_re, -wt_im): */\
		"vmulpd		%%ymm10,%%ymm0,%%ymm0	\n\t"	/* wt_re*[x     ] */\
		"vmulpd		%%ymm11,%%ymm9,%%ymm9	\n\t"	/* wt_im*[y copy] */\
		"vmulpd		%%ymm10,%%ymm1,%%ymm1	\n\t"	/* wt_re*[y     ] */\
		"vmulpd		%%ymm11,%%ymm8,%%ymm8	\n\t"	/* wt_im*[x copy] */\
		"vaddpd		%%ymm9 ,%%ymm0,%%ymm0	\n\t"	/* [a0-3.re] = x*wt_re + y*wt_im */\
		"vsubpd		%%ymm8 ,%%ymm1,%%ymm1	\n\t"	/* [a0-3.im] = y*wt_re - x*wt_im */\
		/* normalize a-pair, compute carryout, compute ROE: */\
		"vaddpd		(%%rbx),%%ymm0,%%ymm0	\n\t"	/* [a0-3.re] + [cy0-3.re] */\
		"vaddpd		(%%rcx),%%ymm1,%%ymm1	\n\t"	/* [a0-3.im] + [cy0-3.im] */\
		"vmovaps	%%ymm0,%%ymm8			\n\t	vmovaps		%%ymm1,%%ymm9		\n\t"	/* copy x|y */\
		"vroundpd	$0,%%ymm0,%%ymm0		\n\t	vroundpd	$0,%%ymm1,%%ymm1	\n\t"	/* temp = DNINT(x|y) */\
		"movq		%[__sign_mask],%%rsi	\n\t"\
		"vmovaps	0x20(%%rdx),%%ymm12	\n\t"	/* [baseinv0-3] */\
		"vsubpd		%%ymm0 ,%%ymm8 ,%%ymm8 	\n\t	vsubpd		%%ymm1 ,%%ymm9 ,%%ymm9 	\n\t"	/* frac = [x - temp] */\
		"vandpd		(%%rsi),%%ymm8 ,%%ymm8 	\n\t	vandpd		(%%rsi),%%ymm9 ,%%ymm9 	\n\t"	/* frac = fabs(frac) */\
		"vmaxpd		%%ymm13,%%ymm8 ,%%ymm13	\n\t	vmaxpd		%%ymm14,%%ymm9 ,%%ymm14	\n\t"	/* if(frac > maxerr) maxerr=frac */\
		"vmovaps	%%ymm0 ,%%ymm8			\n\t	vmovaps		%%ymm1 ,%%ymm9			\n\t"	/* cpy temp */\
		"vmulpd		%%ymm12,%%ymm8 ,%%ymm8 	\n\t	vmulpd		%%ymm12,%%ymm9 ,%%ymm9 	\n\t"	/* temp*baseinv[0] ... inline the remaining +odd_radix offset in addressing */\
		"vmaxpd		%%ymm13,%%ymm14,%%ymm14	\n\t"	/* merge re|im maxerr vectors */\
		"vroundpd	$0,%%ymm8 ,%%ymm8 		\n\t	vroundpd	$0,%%ymm9,%%ymm9		\n\t"	/* [cy0-3.re] = DNINT(temp*baseinv[0]) */\
		"vmovaps	(%%rdx),%%ymm12	\n\t"	/* [base0-3] */\
		"vmovaps	%%ymm8,(%%rbx)			\n\t	vmovaps		%%ymm9,(%%rcx)			\n\t"	/* store [cy0-3.re|im] */\
		"vmulpd		%%ymm12,%%ymm8,%%ymm8	\n\t	vmulpd		%%ymm12,%%ymm9,%%ymm9	\n\t"	/* base[0]*[cy0-3.re|im] */\
		"vsubpd		%%ymm8 ,%%ymm0,%%ymm0	\n\t	vsubpd		%%ymm9 ,%%ymm1,%%ymm1	\n\t"	/* XMM0|1 = [a0-3.re|im] = temp - [cy0-3.re|im]*base[0] */\
		"vmovaps	%%ymm0 ,%%ymm8			\n\t	vmovaps		%%ymm1 ,%%ymm9			\n\t"	/* cpy x|y */\
		/* Forward acyclic-convo weight is (wt_re, +wt_im): */\
		"vmulpd		%%ymm10,%%ymm0,%%ymm0	\n\t"	/* wt_re*[x     ] */\
		"vmulpd		%%ymm11,%%ymm9,%%ymm9	\n\t"	/* wt_im*[y copy] */\
		"vmulpd		%%ymm10,%%ymm1,%%ymm1	\n\t"	/* wt_re*[y     ] */\
		"vmulpd		%%ymm11,%%ymm8,%%ymm8	\n\t"	/* wt_im*[x copy] */\
		"vsubpd		%%ymm9 ,%%ymm0,%%ymm0	\n\t"	/* [a0-3.re] = x*wt_re - y*wt_im */\
		"vaddpd		%%ymm8 ,%%ymm1,%%ymm1	\n\t"	/* [a0-3.im] = y*wt_re + x*wt_im */\
	/* Now do b-quartet: Data in ymm2,ymm3 */\
	"movslq		%[__p1],%%r15	\n\t"	\
	"prefetcht0	(%%r14,%%r15,8)	\n\t"\
		"vmovaps	%%ymm14,%%ymm13			\n\t"	/* maxerr copy */\
		/* Get next set of negacyclic roots: */\
		"vmovaps	0x40(%%rax),%%ymm10		\n\t"	/* c = Re part of 2nd base-root quartet */\
		"vmovaps	0x60(%%rax),%%ymm11		\n\t"	/* s = Im part */\
		/* Data in ymm2,ymm3 */\
		"vmovaps	%%ymm2,%%ymm8			\n\t"	/* x copy */\
		"vmovaps	%%ymm3,%%ymm9			\n\t"	/* y copy */\
		/* Inverse negacyclic weight is (wt_re, -wt_im): */\
		"vmulpd		%%ymm10,%%ymm2,%%ymm2	\n\t"	/* wt_re*[x     ] */\
		"vmulpd		%%ymm11,%%ymm9,%%ymm9	\n\t"	/* wt_im*[y copy] */\
		"vmulpd		%%ymm10,%%ymm3,%%ymm3	\n\t"	/* wt_re*[y     ] */\
		"vmulpd		%%ymm11,%%ymm8,%%ymm8	\n\t"	/* wt_im*[x copy] */\
		"vaddpd		%%ymm9 ,%%ymm2,%%ymm2	\n\t"	/* [a0-3.re] = x*wt_re + y*wt_im */\
		"vsubpd		%%ymm8 ,%%ymm3,%%ymm3	\n\t"	/* [a0-3.im] = y*wt_re - x*wt_im */\
		/* normalize, compute carryout, compute ROE: */\
		"vaddpd		(%%rbx),%%ymm2,%%ymm2	\n\t"	/* [a0-3.re] + [cy0-3.re] */\
		"vaddpd		(%%rcx),%%ymm3,%%ymm3	\n\t"	/* [a0-3.im] + [cy0-3.im] */\
		"vmovaps	%%ymm2,%%ymm8			\n\t	vmovaps		%%ymm3,%%ymm9		\n\t"	/* copy x|y */\
		"vroundpd	$0,%%ymm2,%%ymm2		\n\t	vroundpd	$0,%%ymm3,%%ymm3	\n\t"	/* temp = DNINT(x|y) */\
		"vmovaps	0x20(%%rdx),%%ymm12	\n\t"	/* [baseinv0-3] */\
		"vsubpd		%%ymm2 ,%%ymm8 ,%%ymm8 	\n\t	vsubpd		%%ymm3 ,%%ymm9 ,%%ymm9 	\n\t"	/* frac = [x - temp] */\
		"vandpd		(%%rsi),%%ymm8 ,%%ymm8 	\n\t	vandpd		(%%rsi),%%ymm9 ,%%ymm9 	\n\t"	/* frac = fabs(frac) */\
		"vmaxpd		%%ymm13,%%ymm8 ,%%ymm13	\n\t	vmaxpd		%%ymm14,%%ymm9 ,%%ymm14	\n\t"	/* if(frac > maxerr) maxerr=frac */\
		"vmovaps	%%ymm2 ,%%ymm8			\n\t	vmovaps		%%ymm3 ,%%ymm9			\n\t"	/* cpy temp */\
		"vmulpd		%%ymm12,%%ymm8 ,%%ymm8 	\n\t	vmulpd		%%ymm12,%%ymm9 ,%%ymm9 	\n\t"	/* temp*baseinv[0] ... inline the remaining +odd_radix offset in addressing */\
		"vmaxpd		%%ymm13,%%ymm14,%%ymm14	\n\t"	/* merge re|im maxerr vectors */\
		"vroundpd	$0,%%ymm8 ,%%ymm8 		\n\t	vroundpd	$0,%%ymm9,%%ymm9		\n\t"	/* [cy0-3.re] = DNINT(temp*baseinv[0]) */\
		"vmovaps	(%%rdx),%%ymm12	\n\t"	/* [base0-3] */\
		"vmovaps	%%ymm8,(%%rbx)			\n\t	vmovaps		%%ymm9,(%%rcx)			\n\t"	/* store [cy0-3.re|im] */\
		"vmulpd		%%ymm12,%%ymm8,%%ymm8	\n\t	vmulpd		%%ymm12,%%ymm9,%%ymm9	\n\t"	/* base[0]*[cy0-3.re|im] */\
		"vsubpd		%%ymm8 ,%%ymm2,%%ymm2	\n\t	vsubpd		%%ymm9 ,%%ymm3,%%ymm3	\n\t"	/* XMM0|1 = [a0-3.re|im] = temp - [cy0-3.re|im]*base[0] */\
		"vmovaps	%%ymm2 ,%%ymm8			\n\t	vmovaps		%%ymm3 ,%%ymm9			\n\t"	/* cpy x|y */\
		/* Forward acyclic-convo weight is (wt_re, +wt_im): */\
		"vmulpd		%%ymm10,%%ymm2,%%ymm2	\n\t"	/* wt_re*[x     ] */\
		"vmulpd		%%ymm11,%%ymm9,%%ymm9	\n\t"	/* wt_im*[y copy] */\
		"vmulpd		%%ymm10,%%ymm3,%%ymm3	\n\t"	/* wt_re*[y     ] */\
		"vmulpd		%%ymm11,%%ymm8,%%ymm8	\n\t"	/* wt_im*[x copy] */\
		"vsubpd		%%ymm9 ,%%ymm2,%%ymm2	\n\t"	/* [a0-3.re] = x*wt_re - y*wt_im */\
		"vaddpd		%%ymm8 ,%%ymm3,%%ymm3	\n\t"	/* [a0-3.im] = y*wt_re + x*wt_im */\
	/* Now do c-quartet: Data in ymm4,ymm5 */\
	"movslq		%[__p2],%%r15	\n\t"	\
	"prefetcht0	(%%r14,%%r15,8)	\n\t"\
		"vmovaps	%%ymm14,%%ymm13			\n\t"	/* maxerr copy */\
		/* Get next set of negacyclic roots: */\
		"vmovaps	0x80(%%rax),%%ymm10		\n\t"	/* c = Re part of 3rd base-root quartet */\
		"vmovaps	0xa0(%%rax),%%ymm11		\n\t"	/* s = Im part */\
		/* Data in ymm4,ymm5 */\
		"vmovaps	%%ymm4,%%ymm8			\n\t"	/* x copy */\
		"vmovaps	%%ymm5,%%ymm9			\n\t"	/* y copy */\
		/* Inverse negacyclic weight is (wt_re, -wt_im): */\
		"vmulpd		%%ymm10,%%ymm4,%%ymm4	\n\t"	/* wt_re*[x     ] */\
		"vmulpd		%%ymm11,%%ymm9,%%ymm9	\n\t"	/* wt_im*[y copy] */\
		"vmulpd		%%ymm10,%%ymm5,%%ymm5	\n\t"	/* wt_re*[y     ] */\
		"vmulpd		%%ymm11,%%ymm8,%%ymm8	\n\t"	/* wt_im*[x copy] */\
		"vaddpd		%%ymm9 ,%%ymm4,%%ymm4	\n\t"	/* [a0-3.re] = x*wt_re + y*wt_im */\
		"vsubpd		%%ymm8 ,%%ymm5,%%ymm5	\n\t"	/* [a0-3.im] = y*wt_re - x*wt_im */\
		/* normalize, compute carryout, compute ROE: */\
		"vaddpd		(%%rbx),%%ymm4,%%ymm4	\n\t"	/* [a0-3.re] + [cy0-3.re] */\
		"vaddpd		(%%rcx),%%ymm5,%%ymm5	\n\t"	/* [a0-3.im] + [cy0-3.im] */\
		"vmovaps	%%ymm4,%%ymm8			\n\t	vmovaps		%%ymm5,%%ymm9		\n\t"	/* copy x|y */\
		"vroundpd	$0,%%ymm4,%%ymm4		\n\t	vroundpd	$0,%%ymm5,%%ymm5	\n\t"	/* temp = DNINT(x|y) */\
		"vmovaps	0x20(%%rdx),%%ymm12	\n\t"	/* [baseinv0-3] */\
		"vsubpd		%%ymm4 ,%%ymm8 ,%%ymm8 	\n\t	vsubpd		%%ymm5 ,%%ymm9 ,%%ymm9 	\n\t"	/* frac = [x - temp] */\
		"vandpd		(%%rsi),%%ymm8 ,%%ymm8 	\n\t	vandpd		(%%rsi),%%ymm9 ,%%ymm9 	\n\t"	/* frac = fabs(frac) */\
		"vmaxpd		%%ymm13,%%ymm8 ,%%ymm13	\n\t	vmaxpd		%%ymm14,%%ymm9 ,%%ymm14	\n\t"	/* if(frac > maxerr) maxerr=frac */\
		"vmovaps	%%ymm4 ,%%ymm8			\n\t	vmovaps		%%ymm5 ,%%ymm9			\n\t"	/* cpy temp */\
		"vmulpd		%%ymm12,%%ymm8 ,%%ymm8 	\n\t	vmulpd		%%ymm12,%%ymm9 ,%%ymm9 	\n\t"	/* temp*baseinv[0] ... inline the remaining +odd_radix offset in addressing */\
		"vmaxpd		%%ymm13,%%ymm14,%%ymm14	\n\t"	/* merge re|im maxerr vectors */\
		"vroundpd	$0,%%ymm8 ,%%ymm8 		\n\t	vroundpd	$0,%%ymm9,%%ymm9		\n\t"	/* [cy0-3.re] = DNINT(temp*baseinv[0]) */\
		"vmovaps	(%%rdx),%%ymm12	\n\t"	/* [base0-3] */\
		"vmovaps	%%ymm8,(%%rbx)			\n\t	vmovaps		%%ymm9,(%%rcx)			\n\t"	/* store [cy0-3.re|im] */\
		"vmulpd		%%ymm12,%%ymm8,%%ymm8	\n\t	vmulpd		%%ymm12,%%ymm9,%%ymm9	\n\t"	/* base[0]*[cy0-3.re|im] */\
		"vsubpd		%%ymm8 ,%%ymm4,%%ymm4	\n\t	vsubpd		%%ymm9 ,%%ymm5,%%ymm5	\n\t"	/* XMM0|1 = [a0-3.re|im] = temp - [cy0-3.re|im]*base[0] */\
		"vmovaps	%%ymm4 ,%%ymm8			\n\t	vmovaps		%%ymm5 ,%%ymm9			\n\t"	/* cpy x|y */\
		/* Forward acyclic-convo weight is (wt_re, +wt_im): */\
		"vmulpd		%%ymm10,%%ymm4,%%ymm4	\n\t"	/* wt_re*[x     ] */\
		"vmulpd		%%ymm11,%%ymm9,%%ymm9	\n\t"	/* wt_im*[y copy] */\
		"vmulpd		%%ymm10,%%ymm5,%%ymm5	\n\t"	/* wt_re*[y     ] */\
		"vmulpd		%%ymm11,%%ymm8,%%ymm8	\n\t"	/* wt_im*[x copy] */\
		"vsubpd		%%ymm9 ,%%ymm4,%%ymm4	\n\t"	/* [a0-3.re] = x*wt_re - y*wt_im */\
		"vaddpd		%%ymm8 ,%%ymm5,%%ymm5	\n\t"	/* [a0-3.im] = y*wt_re + x*wt_im */\
	/* Now do d-quartet: Data in ymm6,ymm7 */\
	"movslq		%[__p3],%%r15	\n\t"	\
	"prefetcht0	(%%r14,%%r15,8)	\n\t"\
		"vmovaps	%%ymm14,%%ymm13			\n\t"	/* maxerr copy */\
		/* Get next set of negacyclic roots: */\
		"vmovaps	0xc0(%%rax),%%ymm10		\n\t"	/* c = Re part of 4th base-root quartet */\
		"vmovaps	0xe0(%%rax),%%ymm11		\n\t"	/* s = Im part */\
		/* Data in ymm6,ymm7 */\
		"vmovaps	%%ymm6,%%ymm8			\n\t"	/* x copy */\
		"vmovaps	%%ymm7,%%ymm9			\n\t"	/* y copy */\
		/* Inverse negacyclic weight is (wt_re, -wt_im): */\
		"vmulpd		%%ymm10,%%ymm6,%%ymm6	\n\t"	/* wt_re*[x     ] */\
		"vmulpd		%%ymm11,%%ymm9,%%ymm9	\n\t"	/* wt_im*[y copy] */\
		"vmulpd		%%ymm10,%%ymm7,%%ymm7	\n\t"	/* wt_re*[y     ] */\
		"vmulpd		%%ymm11,%%ymm8,%%ymm8	\n\t"	/* wt_im*[x copy] */\
		"vaddpd		%%ymm9 ,%%ymm6,%%ymm6	\n\t"	/* [a0-3.re] = x*wt_re + y*wt_im */\
		"vsubpd		%%ymm8 ,%%ymm7,%%ymm7	\n\t"	/* [a0-3.im] = y*wt_re - x*wt_im */\
		/* normalize, compute carryout, compute ROE: */\
		"vaddpd		(%%rbx),%%ymm6,%%ymm6	\n\t"	/* [a0-3.re] + [cy0-3.re] */\
		"vaddpd		(%%rcx),%%ymm7,%%ymm7	\n\t"	/* [a0-3.im] + [cy0-3.im] */\
		"vmovaps	%%ymm6,%%ymm8			\n\t	vmovaps		%%ymm7,%%ymm9		\n\t"	/* copy x|y */\
		"vroundpd	$0,%%ymm6,%%ymm6		\n\t	vroundpd	$0,%%ymm7,%%ymm7	\n\t"	/* temp = DNINT(x|y) */\
		"vmovaps	0x20(%%rdx),%%ymm12	\n\t"	/* [baseinv0-3] */\
		"vsubpd		%%ymm6 ,%%ymm8 ,%%ymm8 	\n\t	vsubpd		%%ymm7 ,%%ymm9 ,%%ymm9 	\n\t"	/* frac = [x - temp] */\
		"vandpd		(%%rsi),%%ymm8 ,%%ymm8 	\n\t	vandpd		(%%rsi),%%ymm9 ,%%ymm9 	\n\t"	/* frac = fabs(frac) */\
		"vmaxpd		%%ymm13,%%ymm8 ,%%ymm13	\n\t	vmaxpd		%%ymm14,%%ymm9 ,%%ymm14	\n\t"	/* if(frac > maxerr) maxerr=frac */\
		"vmovaps	%%ymm6 ,%%ymm8			\n\t	vmovaps		%%ymm7 ,%%ymm9			\n\t"	/* cpy temp */\
		"vmulpd		%%ymm12,%%ymm8 ,%%ymm8 	\n\t	vmulpd		%%ymm12,%%ymm9 ,%%ymm9 	\n\t"	/* temp*baseinv[0] ... inline the remaining +odd_radix offset in addressing */\
		"vmaxpd		%%ymm13,%%ymm14,%%ymm14	\n\t"	/* merge re|im maxerr vectors */\
		"vroundpd	$0,%%ymm8 ,%%ymm8 		\n\t	vroundpd	$0,%%ymm9,%%ymm9		\n\t"	/* [cy0-3.re] = DNINT(temp*baseinv[0]) */\
		"vmovaps	(%%rdx),%%ymm12	\n\t"	/* [base0-3] */\
		"vmovaps	%%ymm8,(%%rbx)			\n\t	vmovaps		%%ymm9,(%%rcx)			\n\t"	/* store [cy0-3.re|im] */\
		"vmulpd		%%ymm12,%%ymm8,%%ymm8	\n\t	vmulpd		%%ymm12,%%ymm9,%%ymm9	\n\t"	/* base[0]*[cy0-3.re|im] */\
		"vsubpd		%%ymm8 ,%%ymm6,%%ymm6	\n\t	vsubpd		%%ymm9 ,%%ymm7,%%ymm7	\n\t"	/* XMM0|1 = [a0-3.re|im] = temp - [cy0-3.re|im]*base[0] */\
		"vmovaps	%%ymm6 ,%%ymm8			\n\t	vmovaps		%%ymm7 ,%%ymm9			\n\t"	/* cpy x|y */\
		/* Forward acyclic-convo weight is (wt_re, +wt_im): */\
		"vmulpd		%%ymm10,%%ymm6,%%ymm6	\n\t"	/* wt_re*[x     ] */\
		"vmulpd		%%ymm11,%%ymm9,%%ymm9	\n\t"	/* wt_im*[y copy] */\
		"vmulpd		%%ymm10,%%ymm7,%%ymm7	\n\t"	/* wt_re*[y     ] */\
		"vmulpd		%%ymm11,%%ymm8,%%ymm8	\n\t"	/* wt_im*[x copy] */\
		"vsubpd		%%ymm9 ,%%ymm6,%%ymm6	\n\t"	/* [a0-3.re] = x*wt_re - y*wt_im */\
		"vaddpd		%%ymm8 ,%%ymm7,%%ymm7	\n\t"	/* [a0-3.im] = y*wt_re + x*wt_im */\
		/* Store maxerr: */\
		"movq		%[__half_arr],%%rdx		\n\t"\
		"vmovaps	%%ymm14,-0x40(%%rdx)	\n\t"\
	/* 4-way transpose of outputs (Re, Im parts separately): Inputs from ymm0-7. Outputs into r0/1,2/3,4/5.6/7: */	\
	/* Because default inputs for our 4 x 4 transpose macro (e.g. the one used at start of this carry macro) */\
	/* are into ymm4/2/8/2, munge inputs into that order, resolving name-conflicts via use of the now-available ymm8-15 for outputs: */\
		"movq		%[__data],%%rax			\n\t"\
		"vshufpd	$15,%%ymm2,%%ymm0,%%ymm10					\n\t		vshufpd	$15,%%ymm3,%%ymm1,%%ymm11						\n\t"\
		"vshufpd	$0 ,%%ymm2,%%ymm0,%%ymm0					\n\t		vshufpd	$0 ,%%ymm3,%%ymm1,%%ymm1						\n\t"\
		"vshufpd	$15,%%ymm6,%%ymm4,%%ymm12					\n\t		vshufpd	$15,%%ymm7,%%ymm5,%%ymm13						\n\t"\
		"vshufpd	$0 ,%%ymm6,%%ymm4,%%ymm4					\n\t		vshufpd	$0 ,%%ymm7,%%ymm5,%%ymm5						\n\t"\
		"vperm2f128 $32,%%ymm12,%%ymm10,%%ymm2 		/* Re B	*/	\n\t		vperm2f128 $32,%%ymm13,%%ymm11,%%ymm3		/* Im B	*/	\n\t"\
		"vperm2f128 $49,%%ymm12,%%ymm10,%%ymm10		/* Re D	*/	\n\t		vperm2f128 $49,%%ymm13,%%ymm11,%%ymm11		/* Im D	*/	\n\t"\
		"vperm2f128 $32,%%ymm4 ,%%ymm0 ,%%ymm12		/* Re A	*/	\n\t		vperm2f128 $32,%%ymm5 ,%%ymm1 ,%%ymm13 		/* Im A	*/	\n\t"\
		"vperm2f128 $49,%%ymm4 ,%%ymm0 ,%%ymm0 		/* Re C	*/	\n\t		vperm2f128 $49,%%ymm5 ,%%ymm1 ,%%ymm1		/* Im C	*/	\n\t"\
		/* And write 'em back to memory: */\
		"vmovaps	%%ymm12,     (%%rax)						\n\t		vmovaps	%%ymm13,0x020(%%rax)				\n\t"\
		"vmovaps	%%ymm2 ,0x040(%%rax)						\n\t		vmovaps	%%ymm3 ,0x060(%%rax)				\n\t"\
		"vmovaps	%%ymm0 ,0x080(%%rax)						\n\t		vmovaps	%%ymm1 ,0x0a0(%%rax)				\n\t"\
		"vmovaps	%%ymm10,0x0c0(%%rax)						\n\t		vmovaps	%%ymm11,0x0e0(%%rax)				\n\t"\
		:						/* outputs: none */\
		:	[__data]		"m" (Xdata)	/* All inputs from memory addresses here */\
		,	[__base_root]	"m" (Xbase_root)\
		,	[__cmul_offset] "m" (Xcmul_offset)\
		,	[__cy_re]		"m" (Xcy_re)\
		,	[__cy_im]		"m" (Xcy_im)\
		,	[__half_arr]	"m" (Xhalf_arr)\
		,	[__sign_mask]	"m" (Xsign_mask)\
		/* Prefetch: base address and 3 index offsets */\
		,	[__add0] "m" (Xadd0)\
		,	[__p1] "m" (Xp1)\
		,	[__p2] "m" (Xp2)\
		,	[__p3] "m" (Xp3)\
		: "cc","memory","rax","rbx","rcx","rdx","rsi","r14","r15","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14"	/* Clobbered registers */\
	);\
	}

	/* Non-power-of-2-runlength Fermat-mod acyclic-transform/IBDWT carry macro.

	Key differences vs pow2 version:
	- Use odd_radix as index offset into local storage for IBDWT weights and variable base/baseinv terms;
	- Apply inv/fwd IBDWT weights bookending the negacyclic weights;
	- Value of base/baseinv to be applied to output taken from odd_radix-length array, using same index as for selecting IBDWT weight.

	The array indices i/j/k/lcycle declared int in caller but assumed to have been byte-shift-converted at time this macro called,
	thus can use as complex-address-offsets.  Use bytewise literal offsets to save registers for several args here,as vvv-marked:
												                                           vvvvv The [1,2,3]-multiples of odd_radix assumed << l2_sz_vd on input */
	#define SSE2_fermat_carry_norm_errcheck_X4_hiacc(Xdata,Xbase_root,Xcmul_offset,Xcy_re,Xcy_im,Xodd_radix,Xodd_radm2,Xodd_radm3,Xhalf_arr,Xsign_mask,XicycleA,XicycleB,XicycleC,XicycleD, XjcycleA,XkcycleA,XlcycleA, Xadd0,Xp1,Xp2,Xp3)\
	{\
	__asm__ volatile (\
		"movq		%[__add0],%%rcx		\n\t"/* base address for 4 prefetches-from-main-data-array spread through this macro */\
		/* Base negacyclic roots at this address in [0,2,4,6]*0x20 (Re parts), [1,3,5,7]*0x20 (Imag parts) */\
		"movq	%[__base_root]  ,%%rax			\n\t"	/* Won't need main-array again until output transpose, so re-use rax for base_root */\
		"movslq	%[__cmul_offset],%%rbx			\n\t"\
		"addq	%%rax,%%rbx	\n\t"	/* Index into complex const multipliers block, each applied to 4 sets of base roots */\
		/* Up-multiply quartet of negacyclic roots used in this macro invocation; store sets 2-4 back into mem, keep set 1 in ymm10,11 [that's why we do sets 1/2 after 3/4] */\
		"vmovaps	    (%%rbx),%%ymm10		\n\t	"	/* Multiply by exp(j*I*Pi/2)/RADIX, for j = 0-3 */\
		"vmovaps	0x20(%%rbx),%%ymm11		\n\t	"	/* c = Re(exp) in ymm0, s = Im(exp) in ymm1 */\
		"/* Sets 3/4: */"\
		"vmovaps	 0x80(%%rax),%%ymm0		\n\t	vmovaps	 0xc0(%%rax),%%ymm2			\n\t"	\
		"vmovaps	 0xa0(%%rax),%%ymm1		\n\t	vmovaps	 0xe0(%%rax),%%ymm3			\n\t"	\
		"vmovaps	%%ymm0,%%ymm4			\n\t	vmovaps		%%ymm2,%%ymm6			\n\t"	\
		"vmovaps	%%ymm1,%%ymm5			\n\t	vmovaps		%%ymm3,%%ymm7			\n\t"	\
		"vmulpd		%%ymm10,%%ymm0,%%ymm0	\n\t	vmulpd		%%ymm10,%%ymm2,%%ymm2	\n\t"	\
		"vmulpd		%%ymm11,%%ymm5,%%ymm5	\n\t	vmulpd		%%ymm11,%%ymm7,%%ymm7	\n\t"	\
		"vmulpd		%%ymm10,%%ymm1,%%ymm1	\n\t	vmulpd		%%ymm10,%%ymm3,%%ymm3	\n\t"	\
		"vmulpd		%%ymm11,%%ymm4,%%ymm4	\n\t	vmulpd		%%ymm11,%%ymm6,%%ymm6	\n\t"	\
		"vsubpd		%%ymm5 ,%%ymm0,%%ymm0	\n\t	vsubpd		%%ymm7 ,%%ymm2,%%ymm2	\n\t"	\
		"vaddpd		%%ymm4 ,%%ymm1,%%ymm1	\n\t	vaddpd		%%ymm6 ,%%ymm3,%%ymm3	\n\t"	\
		"vmovaps	%%ymm0 ,0x80(%%rax)		\n\t	vmovaps		%%ymm2 ,0xc0(%%rax)		\n\t"	/* Store result, overwriting input base root */\
		"vmovaps	%%ymm1 ,0xa0(%%rax)		\n\t	vmovaps		%%ymm3 ,0xe0(%%rax)		\n\t"	/* Im part */\
		"/* Sets 1/2: */"\
		"vmovaps	     (%%rax),%%ymm0		\n\t	vmovaps	 0x40(%%rax),%%ymm2			\n\t"	/* x = Re part of 1st base-root quartet */\
		"vmovaps	 0x20(%%rax),%%ymm1		\n\t	vmovaps	 0x60(%%rax),%%ymm3			\n\t"	/* y = Im part */\
		"vmovaps	%%ymm0,%%ymm4			\n\t	vmovaps		%%ymm2,%%ymm6			\n\t"	/* Copy x */\
		"vmovaps	%%ymm1,%%ymm5			\n\t	vmovaps		%%ymm3,%%ymm7			\n\t"	/* Copy y */\
		"vmulpd		%%ymm10,%%ymm0,%%ymm0	\n\t	vmulpd		%%ymm10,%%ymm2,%%ymm2	\n\t"	/* c.x */\
		"vmulpd		%%ymm11,%%ymm5,%%ymm5	\n\t	vmulpd		%%ymm11,%%ymm7,%%ymm7	\n\t"	/* s.y */\
		"vmulpd		%%ymm10,%%ymm1,%%ymm1	\n\t	vmulpd		%%ymm10,%%ymm3,%%ymm3	\n\t"	/* c.y */\
		"vmulpd		%%ymm11,%%ymm4,%%ymm4	\n\t	vmulpd		%%ymm11,%%ymm6,%%ymm6	\n\t"	/* s.x */\
		"vsubpd		%%ymm5 ,%%ymm0,%%ymm10	\n\t	vsubpd		%%ymm7 ,%%ymm2,%%ymm2	\n\t"	/* Out.re = c.x - s.y */\
		"vaddpd		%%ymm4 ,%%ymm1,%%ymm11	\n\t	vaddpd		%%ymm6 ,%%ymm3,%%ymm3	\n\t"	/* Out.im = c.y + s.x */\
		"											vmovaps		%%ymm2 ,0x40(%%rax)		\n\t"	/* Store result, overwriting input base root */\
		"											vmovaps		%%ymm3 ,0x60(%%rax)		\n\t"	/* Im part */\
		"movq		%[__data],%%rax		\n\t"\
		"movq		%[__half_arr],%%rdx		\n\t"\
		"addq		$%c[__odd_radix],%%rdx				\n\t"	/* wt|wtinv|base|baseinv data offset by icycle array slots from resp. base addresses */\
		/* Multiply complex transform outputs [x,y] = [re,im] by inverse IBDWT weights, which include the 2/n scale factor: */\
		"movslq		%[__icycleA],%%rdi		\n\t"	\
		"movslq		%[__icycleB],%%r9 		\n\t"	\
		"movslq		%[__icycleC],%%r8 		\n\t"	\
		"movslq		%[__icycleD],%%r10		\n\t"	\
		"vmovaps	(%%rdx,%%rdi),%%ymm12	\n\t"	/* [wtinv0-3]A */\
		"vmovaps	(%%rdx,%%r9 ),%%ymm13	\n\t"	/* [wtinv0-3]B */\
		"vmovaps	(%%rdx,%%r8 ),%%ymm14	\n\t"	/* [wtinv0-3]C */\
		"vmovaps	(%%rdx,%%r10),%%ymm15	\n\t"	/* [wtinv0-3]D */\
	/* 4-way transpose of inputs (Re, Im parts separately): Inputs from r0/1,2/3,4/5.6/7. Outputs into ymm0-7: */\
	/* Real parts use ymm0,2,4,6, ymm8 as tmp-reg:					Imag parts use ymm1,3,5,7, ymm9 as tm-reg: */\
		"vmovaps	     (%%rax),%%ymm4						\n\t		vmovaps	0x020(%%rax),%%ymm5							\n\t"\
		"vmovaps	0x040(%%rax),%%ymm2						\n\t		vmovaps	0x060(%%rax),%%ymm3							\n\t"\
		"vmulpd		%%ymm12,%%ymm4,%%ymm4					\n\t		vmulpd		%%ymm12,%%ymm5,%%ymm5					\n\t"\
		"vmulpd		%%ymm13,%%ymm2,%%ymm2					\n\t		vmulpd		%%ymm13,%%ymm3,%%ymm3					\n\t"\
		"vshufpd	$15,%%ymm2,%%ymm4,%%ymm6				\n\t		vshufpd	$15,%%ymm3,%%ymm5,%%ymm7					\n\t"\
		"vshufpd	$0 ,%%ymm2,%%ymm4,%%ymm4				\n\t		vshufpd	$0 ,%%ymm3,%%ymm5,%%ymm5					\n\t"\
		"vmovaps	0x080(%%rax),%%ymm8						\n\t		vmovaps	0x0a0(%%rax),%%ymm9							\n\t"\
		"vmovaps	0x0c0(%%rax),%%ymm2						\n\t		vmovaps	0x0e0(%%rax),%%ymm3							\n\t"\
		"vmulpd		%%ymm14,%%ymm8,%%ymm8					\n\t		vmulpd		%%ymm14,%%ymm9,%%ymm9					\n\t"\
		"vmulpd		%%ymm15,%%ymm2,%%ymm2					\n\t		vmulpd		%%ymm15,%%ymm3,%%ymm3					\n\t"\
		"vshufpd	$15,%%ymm2,%%ymm8,%%ymm0				\n\t		vshufpd	$15,%%ymm3,%%ymm9,%%ymm1					\n\t"\
		"vshufpd	$0 ,%%ymm2,%%ymm8,%%ymm8				\n\t		vshufpd	$0 ,%%ymm3,%%ymm9,%%ymm9					\n\t"\
		"vperm2f128 $32,%%ymm0,%%ymm6,%%ymm2	/* Re B	*/	\n\t		vperm2f128 $32,%%ymm1,%%ymm7,%%ymm3		/* Im B	*/	\n\t"\
		"vperm2f128 $49,%%ymm0,%%ymm6,%%ymm6	/* Re D	*/	\n\t		vperm2f128 $49,%%ymm1,%%ymm7,%%ymm7		/* Im D	*/	\n\t"\
		"vperm2f128 $32,%%ymm8,%%ymm4,%%ymm0	/* Re A	*/	\n\t		vperm2f128 $32,%%ymm9,%%ymm5,%%ymm1 	/* Im A	*/	\n\t"\
		"vperm2f128 $49,%%ymm8,%%ymm4,%%ymm4	/* Re C	*/	\n\t		vperm2f128 $49,%%ymm9,%%ymm5,%%ymm5 	/* Im C	*/	\n\t"\
		"subq		$%c[__odd_radix],%%rdx				\n\t"\
		/* Base negacyclic roots at this address in [0,2,4,6]*0x20 (Re parts), [1,3,5,7]*0x20 (Imag parts) */\
		"movq		%[__sign_mask],%%rax	\n\t"\
		"vmovaps	(%%rax),%%ymm15	\n\t"	/* ymm15 free for rest of way; use to store sign_mask needed for floating ABS */\
		"movq	%[__base_root] ,%%rax		\n\t"	/* Won't need main-array again until output transpose, so re-use rax for base_root */\
		/* half_arr[0,1,2,3] = [base*2, baseinv*2,wt_re*2,wt_im*2] */\
		"vmovaps	-0x40(%%rdx),%%ymm13	\n\t"	/* XMM13 = maxerr */\
		"addq		%%rdi,%%rdx				\n\t"	/* icycle assumed already in left-shifted ptr-byte-offset form */\
		"vmovaps	%%ymm13,%%ymm14			\n\t"	/* maxerr copy */\
		"movq		%[__cy_re],%%rbx		\n\t"\
	/* Do a-quartet: Data in ymm0,ymm1 */\
	"prefetcht0	(%%rcx)	\n\t"\
		/* For a-quartet, needed negacyclic root already in ymm10/11: */\
		/* Data in ymm0,ymm1 */\
		"vmovaps	%%ymm0,%%ymm8			\n\t"	/* x copy */\
		"vmovaps	%%ymm1,%%ymm9			\n\t"	/* y copy */\
		/* Inverse negacyclic weight is (wt_re, -wt_im): */\
		"vmulpd		%%ymm10,%%ymm0,%%ymm0	\n\t"	/* wt_re*[x     ] */\
		"vmulpd		%%ymm11,%%ymm9,%%ymm9	\n\t"	/* wt_im*[y copy] */\
		"vmulpd		%%ymm10,%%ymm1,%%ymm1	\n\t"	/* wt_re*[y     ] */\
		"vmulpd		%%ymm11,%%ymm8,%%ymm8	\n\t"	/* wt_im*[x copy] */\
		"vaddpd		%%ymm9 ,%%ymm0,%%ymm0	\n\t"	/* [a0-3.re] = x*wt_re + y*wt_im */\
		"vsubpd		%%ymm8 ,%%ymm1,%%ymm1	\n\t"	/* [a0-3.im] = y*wt_re - x*wt_im */\
		/* normalize a-pair, compute carryout, compute ROE: */\
		"vaddpd		           (%%rbx),%%ymm0,%%ymm0	\n\t"	/* [a0-3.re] + [cy0-3.re] */\
		"vaddpd		%c[__cy_im](%%rbx),%%ymm1,%%ymm1	\n\t"	/* [a0-3.im] + [cy0-3.im] */\
		"vmovaps	%%ymm0,%%ymm8			\n\t	vmovaps		%%ymm1,%%ymm9		\n\t"	/* copy x|y */\
		"vroundpd	$0,%%ymm0,%%ymm0		\n\t	vroundpd	$0,%%ymm1,%%ymm1	\n\t"	/* temp = DNINT(x|y) */\
		"vmovaps	%c[__odd_radm3](%%rdx),%%ymm12	\n\t"	/* [baseinv0-3] */\
		"vsubpd		%%ymm0 ,%%ymm8 ,%%ymm8 	\n\t	vsubpd		%%ymm1 ,%%ymm9 ,%%ymm9 	\n\t"	/* frac = [x - temp] */\
		"vandpd		%%ymm15,%%ymm8 ,%%ymm8 	\n\t	vandpd		%%ymm15,%%ymm9 ,%%ymm9 	\n\t"	/* frac = fabs(frac) */\
		"vmaxpd		%%ymm13,%%ymm8 ,%%ymm13	\n\t	vmaxpd		%%ymm14,%%ymm9 ,%%ymm14	\n\t"	/* if(frac > maxerr) maxerr=frac */\
		"vmovaps	%%ymm0 ,%%ymm8			\n\t	vmovaps		%%ymm1 ,%%ymm9			\n\t"	/* cpy temp */\
		"vmulpd		%%ymm12,%%ymm8 ,%%ymm8 	\n\t	vmulpd		%%ymm12,%%ymm9 ,%%ymm9 	\n\t"	/* temp*baseinv[0] ... inline the remaining +odd_radix offset in addressing */\
		"vmaxpd		%%ymm13,%%ymm14,%%ymm14	\n\t"	/* merge re|im maxerr vectors */\
		"vroundpd	$0,%%ymm8 ,%%ymm8 		\n\t	vroundpd	$0,%%ymm9,%%ymm9		\n\t"	/* [cy0-3.re] = DNINT(temp*baseinv[0]) */\
		"vmovaps	%c[__odd_radm2](%%rdx),%%ymm12	\n\t"	/* [base0-3] */\
		"vmovaps	%%ymm8,(%%rbx)			\n\t	vmovaps		%%ymm9,%c[__cy_im](%%rbx)\n\t"	/* store [cy0-3.re|im] */\
		"vmulpd		%%ymm12,%%ymm8,%%ymm8	\n\t	vmulpd		%%ymm12,%%ymm9,%%ymm9	\n\t"	/* base[0]*[cy0-3.re|im] */\
		"vsubpd		%%ymm8 ,%%ymm0,%%ymm0	\n\t	vsubpd		%%ymm9 ,%%ymm1,%%ymm1	\n\t"	/* XMM0|1 = [a0-3.re|im] = temp - [cy0-3.re|im]*base[0] */\
		"vmovaps	%%ymm0 ,%%ymm8			\n\t	vmovaps		%%ymm1 ,%%ymm9			\n\t"	/* cpy x|y */\
		/* Forward acyclic-convo weight is (wt_re, +wt_im): */\
		"vmulpd		%%ymm10,%%ymm0,%%ymm0	\n\t"	/* wt_re*[x     ] */\
		"vmulpd		%%ymm11,%%ymm9,%%ymm9	\n\t"	/* wt_im*[y copy] */\
		"vmulpd		%%ymm10,%%ymm1,%%ymm1	\n\t"	/* wt_re*[y     ] */\
		"vmulpd		%%ymm11,%%ymm8,%%ymm8	\n\t"	/* wt_im*[x copy] */\
		"vsubpd		%%ymm9 ,%%ymm0,%%ymm0	\n\t"	/* [a0-3.re] = x*wt_re - y*wt_im */\
		"vaddpd		%%ymm8 ,%%ymm1,%%ymm1	\n\t"	/* [a0-3.im] = y*wt_re + x*wt_im */\
	/* Now do b-quartet: Data in ymm2,ymm3 */\
	"movslq		%[__p1],%%r15	\n\t"	\
	"prefetcht0	(%%rcx,%%r15,8)	\n\t"\
		"movq		%[__half_arr],%%rdx		\n\t"\
		"vmovaps	%%ymm14,%%ymm13			\n\t"	/* maxerr copy */\
		"movslq		%[__jcycleA],%%rdi		\n\t"\
		"addq		%%rdi,%%rdx				\n\t"	/* jcycle assumed already in left-shifted ptr-byte-offset form */\
		/* Get next set of negacyclic roots: */\
		"vmovaps	0x40(%%rax),%%ymm10		\n\t"	/* c = Re part of 2nd base-root quartet */\
		"vmovaps	0x60(%%rax),%%ymm11		\n\t"	/* s = Im part */\
		/* Data in ymm2,ymm3 */\
		"vmovaps	%%ymm2,%%ymm8			\n\t"	/* x copy */\
		"vmovaps	%%ymm3,%%ymm9			\n\t"	/* y copy */\
		/* Inverse negacyclic weight is (wt_re, -wt_im): */\
		"vmulpd		%%ymm10,%%ymm2,%%ymm2	\n\t"	/* wt_re*[x     ] */\
		"vmulpd		%%ymm11,%%ymm9,%%ymm9	\n\t"	/* wt_im*[y copy] */\
		"vmulpd		%%ymm10,%%ymm3,%%ymm3	\n\t"	/* wt_re*[y     ] */\
		"vmulpd		%%ymm11,%%ymm8,%%ymm8	\n\t"	/* wt_im*[x copy] */\
		"vaddpd		%%ymm9 ,%%ymm2,%%ymm2	\n\t"	/* [a0-3.re] = x*wt_re + y*wt_im */\
		"vsubpd		%%ymm8 ,%%ymm3,%%ymm3	\n\t"	/* [a0-3.im] = y*wt_re - x*wt_im */\
		/* normalize, compute carryout, compute ROE: */\
		"vaddpd		           (%%rbx),%%ymm2,%%ymm2	\n\t"	/* [a0-3.re] + [cy0-3.re] */\
		"vaddpd		%c[__cy_im](%%rbx),%%ymm3,%%ymm3	\n\t"	/* [a0-3.im] + [cy0-3.im] */\
		"vmovaps	%%ymm2,%%ymm8			\n\t	vmovaps		%%ymm3,%%ymm9		\n\t"	/* copy x|y */\
		"vroundpd	$0,%%ymm2,%%ymm2		\n\t	vroundpd	$0,%%ymm3,%%ymm3	\n\t"	/* temp = DNINT(x|y) */\
		"vmovaps	%c[__odd_radm3](%%rdx),%%ymm12	\n\t"	/* [baseinv0-3] */\
		"vsubpd		%%ymm2 ,%%ymm8 ,%%ymm8 	\n\t	vsubpd		%%ymm3 ,%%ymm9 ,%%ymm9 	\n\t"	/* frac = [x - temp] */\
		"vandpd		%%ymm15,%%ymm8 ,%%ymm8 	\n\t	vandpd		%%ymm15,%%ymm9 ,%%ymm9 	\n\t"	/* frac = fabs(frac) */\
		"vmaxpd		%%ymm13,%%ymm8 ,%%ymm13	\n\t	vmaxpd		%%ymm14,%%ymm9 ,%%ymm14	\n\t"	/* if(frac > maxerr) maxerr=frac */\
		"vmovaps	%%ymm2 ,%%ymm8			\n\t	vmovaps		%%ymm3 ,%%ymm9			\n\t"	/* cpy temp */\
		"vmulpd		%%ymm12,%%ymm8 ,%%ymm8 	\n\t	vmulpd		%%ymm12,%%ymm9 ,%%ymm9 	\n\t"	/* temp*baseinv[0] ... inline the remaining +odd_radix offset in addressing */\
		"vmaxpd		%%ymm13,%%ymm14,%%ymm14	\n\t"	/* merge re|im maxerr vectors */\
		"vroundpd	$0,%%ymm8 ,%%ymm8 		\n\t	vroundpd	$0,%%ymm9,%%ymm9		\n\t"	/* [cy0-3.re] = DNINT(temp*baseinv[0]) */\
		"vmovaps	%c[__odd_radm2](%%rdx),%%ymm12	\n\t"	/* [base0-3] */\
		"vmovaps	%%ymm8,(%%rbx)			\n\t	vmovaps		%%ymm9,%c[__cy_im](%%rbx)\n\t"	/* store [cy0-3.re|im] */\
		"vmulpd		%%ymm12,%%ymm8,%%ymm8	\n\t	vmulpd		%%ymm12,%%ymm9,%%ymm9	\n\t"	/* base[0]*[cy0-3.re|im] */\
		"vsubpd		%%ymm8 ,%%ymm2,%%ymm2	\n\t	vsubpd		%%ymm9 ,%%ymm3,%%ymm3	\n\t"	/* XMM0|1 = [a0-3.re|im] = temp - [cy0-3.re|im]*base[0] */\
		"vmovaps	%%ymm2 ,%%ymm8			\n\t	vmovaps		%%ymm3 ,%%ymm9			\n\t"	/* cpy x|y */\
		/* Forward acyclic-convo weight is (wt_re, +wt_im): */\
		"vmulpd		%%ymm10,%%ymm2,%%ymm2	\n\t"	/* wt_re*[x     ] */\
		"vmulpd		%%ymm11,%%ymm9,%%ymm9	\n\t"	/* wt_im*[y copy] */\
		"vmulpd		%%ymm10,%%ymm3,%%ymm3	\n\t"	/* wt_re*[y     ] */\
		"vmulpd		%%ymm11,%%ymm8,%%ymm8	\n\t"	/* wt_im*[x copy] */\
		"vsubpd		%%ymm9 ,%%ymm2,%%ymm2	\n\t"	/* [a0-3.re] = x*wt_re - y*wt_im */\
		"vaddpd		%%ymm8 ,%%ymm3,%%ymm3	\n\t"	/* [a0-3.im] = y*wt_re + x*wt_im */\
	/* Now do c-quartet: Data in ymm4,ymm5 */\
	"movslq		%[__p2],%%r15	\n\t"	\
	"prefetcht0	(%%rcx,%%r15,8)	\n\t"\
		"movq		%[__half_arr],%%rdx		\n\t"\
		"vmovaps	%%ymm14,%%ymm13			\n\t"	/* maxerr copy */\
		"movslq		%[__kcycleA],%%rdi		\n\t"\
		"addq		%%rdi,%%rdx				\n\t"	/* kcycle assumed already in left-shifted ptr-byte-offset form */\
		/* Get next set of negacyclic roots: */\
		"vmovaps	0x80(%%rax),%%ymm10		\n\t"	/* c = Re part of 3rd base-root quartet */\
		"vmovaps	0xa0(%%rax),%%ymm11		\n\t"	/* s = Im part */\
		/* Data in ymm4,ymm5 */\
		"vmovaps	%%ymm4,%%ymm8			\n\t"	/* x copy */\
		"vmovaps	%%ymm5,%%ymm9			\n\t"	/* y copy */\
		/* Inverse negacyclic weight is (wt_re, -wt_im): */\
		"vmulpd		%%ymm10,%%ymm4,%%ymm4	\n\t"	/* wt_re*[x     ] */\
		"vmulpd		%%ymm11,%%ymm9,%%ymm9	\n\t"	/* wt_im*[y copy] */\
		"vmulpd		%%ymm10,%%ymm5,%%ymm5	\n\t"	/* wt_re*[y     ] */\
		"vmulpd		%%ymm11,%%ymm8,%%ymm8	\n\t"	/* wt_im*[x copy] */\
		"vaddpd		%%ymm9 ,%%ymm4,%%ymm4	\n\t"	/* [a0-3.re] = x*wt_re + y*wt_im */\
		"vsubpd		%%ymm8 ,%%ymm5,%%ymm5	\n\t"	/* [a0-3.im] = y*wt_re - x*wt_im */\
		/* normalize, compute carryout, compute ROE: */\
		"vaddpd		           (%%rbx),%%ymm4,%%ymm4	\n\t"	/* [a0-3.re] + [cy0-3.re] */\
		"vaddpd		%c[__cy_im](%%rbx),%%ymm5,%%ymm5	\n\t"	/* [a0-3.im] + [cy0-3.im] */\
		"vmovaps	%%ymm4,%%ymm8			\n\t	vmovaps		%%ymm5,%%ymm9		\n\t"	/* copy x|y */\
		"vroundpd	$0,%%ymm4,%%ymm4		\n\t	vroundpd	$0,%%ymm5,%%ymm5	\n\t"	/* temp = DNINT(x|y) */\
		"vmovaps	%c[__odd_radm3](%%rdx),%%ymm12	\n\t"	/* [baseinv0-3] */\
		"vsubpd		%%ymm4 ,%%ymm8 ,%%ymm8 	\n\t	vsubpd		%%ymm5 ,%%ymm9 ,%%ymm9 	\n\t"	/* frac = [x - temp] */\
		"vandpd		%%ymm15,%%ymm8 ,%%ymm8 	\n\t	vandpd		%%ymm15,%%ymm9 ,%%ymm9 	\n\t"	/* frac = fabs(frac) */\
		"vmaxpd		%%ymm13,%%ymm8 ,%%ymm13	\n\t	vmaxpd		%%ymm14,%%ymm9 ,%%ymm14	\n\t"	/* if(frac > maxerr) maxerr=frac */\
		"vmovaps	%%ymm4 ,%%ymm8			\n\t	vmovaps		%%ymm5 ,%%ymm9			\n\t"	/* cpy temp */\
		"vmulpd		%%ymm12,%%ymm8 ,%%ymm8 	\n\t	vmulpd		%%ymm12,%%ymm9 ,%%ymm9 	\n\t"	/* temp*baseinv[0] ... inline the remaining +odd_radix offset in addressing */\
		"vmaxpd		%%ymm13,%%ymm14,%%ymm14	\n\t"	/* merge re|im maxerr vectors */\
		"vroundpd	$0,%%ymm8 ,%%ymm8 		\n\t	vroundpd	$0,%%ymm9,%%ymm9		\n\t"	/* [cy0-3.re] = DNINT(temp*baseinv[0]) */\
		"vmovaps	%c[__odd_radm2](%%rdx),%%ymm12	\n\t"	/* [base0-3] */\
		"vmovaps	%%ymm8,(%%rbx)			\n\t	vmovaps		%%ymm9,%c[__cy_im](%%rbx)\n\t"	/* store [cy0-3.re|im] */\
		"vmulpd		%%ymm12,%%ymm8,%%ymm8	\n\t	vmulpd		%%ymm12,%%ymm9,%%ymm9	\n\t"	/* base[0]*[cy0-3.re|im] */\
		"vsubpd		%%ymm8 ,%%ymm4,%%ymm4	\n\t	vsubpd		%%ymm9 ,%%ymm5,%%ymm5	\n\t"	/* XMM0|1 = [a0-3.re|im] = temp - [cy0-3.re|im]*base[0] */\
		"vmovaps	%%ymm4 ,%%ymm8			\n\t	vmovaps		%%ymm5 ,%%ymm9			\n\t"	/* cpy x|y */\
		/* Forward acyclic-convo weight is (wt_re, +wt_im): */\
		"vmulpd		%%ymm10,%%ymm4,%%ymm4	\n\t"	/* wt_re*[x     ] */\
		"vmulpd		%%ymm11,%%ymm9,%%ymm9	\n\t"	/* wt_im*[y copy] */\
		"vmulpd		%%ymm10,%%ymm5,%%ymm5	\n\t"	/* wt_re*[y     ] */\
		"vmulpd		%%ymm11,%%ymm8,%%ymm8	\n\t"	/* wt_im*[x copy] */\
		"vsubpd		%%ymm9 ,%%ymm4,%%ymm4	\n\t"	/* [a0-3.re] = x*wt_re - y*wt_im */\
		"vaddpd		%%ymm8 ,%%ymm5,%%ymm5	\n\t"	/* [a0-3.im] = y*wt_re + x*wt_im */\
	/* Now do d-quartet: Data in ymm6,ymm7 */\
	"movslq		%[__p3],%%r15	\n\t"	\
	"prefetcht0	(%%rcx,%%r15,8)	\n\t"\
		"movq		%[__half_arr],%%rdx		\n\t"\
		"vmovaps	%%ymm14,%%ymm13			\n\t"	/* maxerr copy */\
		"movslq		%[__lcycleA],%%rdi		\n\t"\
		"addq		%%rdi,%%rdx				\n\t"	/* lcycle assumed already in left-shifted ptr-byte-offset form */\
		/* Get next set of negacyclic roots: */\
		"vmovaps	0xc0(%%rax),%%ymm10		\n\t"	/* c = Re part of 4th base-root quartet */\
		"vmovaps	0xe0(%%rax),%%ymm11		\n\t"	/* s = Im part */\
		/* Data in ymm6,ymm7 */\
		"vmovaps	%%ymm6,%%ymm8			\n\t"	/* x copy */\
		"vmovaps	%%ymm7,%%ymm9			\n\t"	/* y copy */\
		/* Inverse negacyclic weight is (wt_re, -wt_im): */\
		"vmulpd		%%ymm10,%%ymm6,%%ymm6	\n\t"	/* wt_re*[x     ] */\
		"vmulpd		%%ymm11,%%ymm9,%%ymm9	\n\t"	/* wt_im*[y copy] */\
		"vmulpd		%%ymm10,%%ymm7,%%ymm7	\n\t"	/* wt_re*[y     ] */\
		"vmulpd		%%ymm11,%%ymm8,%%ymm8	\n\t"	/* wt_im*[x copy] */\
		"vaddpd		%%ymm9 ,%%ymm6,%%ymm6	\n\t"	/* [a0-3.re] = x*wt_re + y*wt_im */\
		"vsubpd		%%ymm8 ,%%ymm7,%%ymm7	\n\t"	/* [a0-3.im] = y*wt_re - x*wt_im */\
		/* normalize, compute carryout, compute ROE: */\
		"vaddpd		           (%%rbx),%%ymm6,%%ymm6	\n\t"	/* [a0-3.re] + [cy0-3.re] */\
		"vaddpd		%c[__cy_im](%%rbx),%%ymm7,%%ymm7	\n\t"	/* [a0-3.im] + [cy0-3.im] */\
		"vmovaps	%%ymm6,%%ymm8			\n\t	vmovaps		%%ymm7,%%ymm9		\n\t"	/* copy x|y */\
		"vroundpd	$0,%%ymm6,%%ymm6		\n\t	vroundpd	$0,%%ymm7,%%ymm7	\n\t"	/* temp = DNINT(x|y) */\
		"vmovaps	%c[__odd_radm3](%%rdx),%%ymm12	\n\t"	/* [baseinv0-3] */\
		"vsubpd		%%ymm6 ,%%ymm8 ,%%ymm8 	\n\t	vsubpd		%%ymm7 ,%%ymm9 ,%%ymm9 	\n\t"	/* frac = [x - temp] */\
		"vandpd		%%ymm15,%%ymm8 ,%%ymm8 	\n\t	vandpd		%%ymm15,%%ymm9 ,%%ymm9 	\n\t"	/* frac = fabs(frac) */\
		"vmaxpd		%%ymm13,%%ymm8 ,%%ymm13	\n\t	vmaxpd		%%ymm14,%%ymm9 ,%%ymm14	\n\t"	/* if(frac > maxerr) maxerr=frac */\
		"vmovaps	%%ymm6 ,%%ymm8			\n\t	vmovaps		%%ymm7 ,%%ymm9			\n\t"	/* cpy temp */\
		"vmulpd		%%ymm12,%%ymm8 ,%%ymm8 	\n\t	vmulpd		%%ymm12,%%ymm9 ,%%ymm9 	\n\t"	/* temp*baseinv[0] ... inline the remaining +odd_radix offset in addressing */\
		"vmaxpd		%%ymm13,%%ymm14,%%ymm14	\n\t"	/* merge re|im maxerr vectors */\
		"vroundpd	$0,%%ymm8 ,%%ymm8 		\n\t	vroundpd	$0,%%ymm9,%%ymm9		\n\t"	/* [cy0-3.re] = DNINT(temp*baseinv[0]) */\
		"vmovaps	%c[__odd_radm2](%%rdx),%%ymm12	\n\t"	/* [base0-3] */\
		"vmovaps	%%ymm8,(%%rbx)			\n\t	vmovaps		%%ymm9,%c[__cy_im](%%rbx)\n\t"	/* store [cy0-3.re|im] */\
		"vmulpd		%%ymm12,%%ymm8,%%ymm8	\n\t	vmulpd		%%ymm12,%%ymm9,%%ymm9	\n\t"	/* base[0]*[cy0-3.re|im] */\
		"vsubpd		%%ymm8 ,%%ymm6,%%ymm6	\n\t	vsubpd		%%ymm9 ,%%ymm7,%%ymm7	\n\t"	/* XMM0|1 = [a0-3.re|im] = temp - [cy0-3.re|im]*base[0] */\
		"vmovaps	%%ymm6 ,%%ymm8			\n\t	vmovaps		%%ymm7 ,%%ymm9			\n\t"	/* cpy x|y */\
		/* Forward acyclic-convo weight is (wt_re, +wt_im): */\
		"vmulpd		%%ymm10,%%ymm6,%%ymm6	\n\t"	/* wt_re*[x     ] */\
		"vmulpd		%%ymm11,%%ymm9,%%ymm9	\n\t"	/* wt_im*[y copy] */\
		"vmulpd		%%ymm10,%%ymm7,%%ymm7	\n\t"	/* wt_re*[y     ] */\
		"vmulpd		%%ymm11,%%ymm8,%%ymm8	\n\t"	/* wt_im*[x copy] */\
		"vsubpd		%%ymm9 ,%%ymm6,%%ymm6	\n\t"	/* [a0-3.re] = x*wt_re - y*wt_im */\
		"vaddpd		%%ymm8 ,%%ymm7,%%ymm7	\n\t"	/* [a0-3.im] = y*wt_re + x*wt_im */\
		/* Store maxerr: */\
		"movq		%[__half_arr],%%rdx		\n\t"\
		"vmovaps	%%ymm14,-0x40(%%rdx)	\n\t"\
	/* 4-way transpose of outputs (Re, Im parts separately): Inputs from ymm0-7. Outputs into r0/1,2/3,4/5.6/7: */	\
	/* Because default inputs for our 4 x 4 transpose macro (e.g. the one used at start of this carry macro) */\
	/* are into ymm4/2/8/2, munge inputs into that order, resolving name-conflicts via use of the now-available ymm8-15 for outputs: */\
		"movq		%[__data],%%rax			\n\t"\
		"vshufpd	$15,%%ymm2,%%ymm0,%%ymm10					\n\t		vshufpd	$15,%%ymm3,%%ymm1,%%ymm11						\n\t"\
		"vshufpd	$0 ,%%ymm2,%%ymm0,%%ymm0					\n\t		vshufpd	$0 ,%%ymm3,%%ymm1,%%ymm1						\n\t"\
		"vshufpd	$15,%%ymm6,%%ymm4,%%ymm12					\n\t		vshufpd	$15,%%ymm7,%%ymm5,%%ymm13						\n\t"\
		"vshufpd	$0 ,%%ymm6,%%ymm4,%%ymm4					\n\t		vshufpd	$0 ,%%ymm7,%%ymm5,%%ymm5						\n\t"\
		"vperm2f128 $32,%%ymm12,%%ymm10,%%ymm2 		/* Re B	*/	\n\t		vperm2f128 $32,%%ymm13,%%ymm11,%%ymm3		/* Im B	*/	\n\t"\
		"vperm2f128 $49,%%ymm12,%%ymm10,%%ymm10		/* Re D	*/	\n\t		vperm2f128 $49,%%ymm13,%%ymm11,%%ymm11		/* Im D	*/	\n\t"\
		"vperm2f128 $32,%%ymm4 ,%%ymm0 ,%%ymm12		/* Re A	*/	\n\t		vperm2f128 $32,%%ymm5 ,%%ymm1 ,%%ymm13 		/* Im A	*/	\n\t"\
		"vperm2f128 $49,%%ymm4 ,%%ymm0 ,%%ymm0 		/* Re C	*/	\n\t		vperm2f128 $49,%%ymm5 ,%%ymm1 ,%%ymm1		/* Im C	*/	\n\t"\
		/* Multiply normalized, re-permuted transform outputs by forward IBDWT weights: */\
		"movslq		%[__icycleA],%%rdi		\n\t"	\
		"vmovaps	(%%rdx,%%rdi),%%ymm4	\n\t"	/* [wt0-3]A */\
		"vmovaps	(%%rdx,%%r9 ),%%ymm5	\n\t"	/* [wt0-3]B */\
		"vmovaps	(%%rdx,%%r8 ),%%ymm6	\n\t"	/* [wt0-3]C */\
		"vmovaps	(%%rdx,%%r10),%%ymm7	\n\t"	/* [wt0-3]D */\
		"vmulpd		%%ymm4,%%ymm12,%%ymm12						\n\t		vmulpd		%%ymm4,%%ymm13,%%ymm13			\n\t"\
		"vmulpd		%%ymm5,%%ymm2 ,%%ymm2 						\n\t		vmulpd		%%ymm5,%%ymm3 ,%%ymm3 			\n\t"\
		"vmulpd		%%ymm6,%%ymm0 ,%%ymm0 						\n\t		vmulpd		%%ymm6,%%ymm1 ,%%ymm1 			\n\t"\
		"vmulpd		%%ymm7,%%ymm10,%%ymm10						\n\t		vmulpd		%%ymm7,%%ymm11,%%ymm11			\n\t"\
		/* And write 'em back to memory: */\
		"vmovaps	%%ymm12,     (%%rax)						\n\t		vmovaps	%%ymm13,0x020(%%rax)				\n\t"\
		"vmovaps	%%ymm2 ,0x040(%%rax)						\n\t		vmovaps	%%ymm3 ,0x060(%%rax)				\n\t"\
		"vmovaps	%%ymm0 ,0x080(%%rax)						\n\t		vmovaps	%%ymm1 ,0x0a0(%%rax)				\n\t"\
		"vmovaps	%%ymm10,0x0c0(%%rax)						\n\t		vmovaps	%%ymm11,0x0e0(%%rax)				\n\t"\
		:						/* outputs: none */\
		:	[__data]		"m" (Xdata)	/* All inputs from memory addresses here */\
		,	[__base_root]	"m" (Xbase_root)\
		,	[__cmul_offset] "m" (Xcmul_offset)\
		,	[__cy_re]		"m" (Xcy_re)\
		,	[__cy_im]	"e" (Xcy_im)	/* Use literal-byte-offset for this ome to save a reg */\
		/* [1,2,3]*odd_radix are the index offsets to the wtinv, base, and base_inv values, resp. - assumed << l2_sz_vd on input: */\
		,	[__odd_radix]   "e" (Xodd_radix)\
		,	[__odd_radm2]   "e" (Xodd_radm2)\
		,	[__odd_radm3]   "e" (Xodd_radm3)\
		,	[__half_arr]	"m" (Xhalf_arr)\
		,	[__sign_mask]	"m" (Xsign_mask)\
		/* Need quartet of ascending [modulo odd_radix] icycle indices for IBDWT weights: */\
		,	[__icycleA]		"m" (XicycleA)\
		,	[__icycleB]		"m" (XicycleB)\
		,	[__icycleC]		"m" (XicycleC)\
		,	[__icycleD]		"m" (XicycleD)\
		/* Need quartet of same-index [i,j,k,l]cycle indices for negacyclic weights and base/baseinv normalizations: */\
		,	[__jcycleA]		"m" (XjcycleA)\
		,	[__kcycleA]		"m" (XkcycleA)\
		,	[__lcycleA]		"m" (XlcycleA)\
		/* Prefetch: base address and 3 index offsets */\
		,	[__add0] "m" (Xadd0)\
		,	[__p1] "m" (Xp1)\
		,	[__p2] "m" (Xp2)\
		,	[__p3] "m" (Xp3)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","r8","r9","r10","r15","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"   /* Clobbered registers */\
	);\
	}

	#define SSE2_fermat_carry_init_loacc(Xbase_root)\
	{\
	__asm__ volatile (\
		"movq		%[__base_root] ,%%rax	\n\t	"	/* Base negacyclic roots at this address +8*0x20 (Re parts), +9*0x20 (Imag parts) */\
		"vmovaps	0x100(%%rax),%%ymm10	\n\t	"	/* Multiply by exp(j*I*Pi/2)/RADIX, for j = 0-3 */\
		"vmovaps	0x120(%%rax),%%ymm11	\n\t	"	/* c = Re(exp) in ymm0, s = Im(exp) in ymm1 */\
		"vmovaps	     (%%rax),%%ymm0		\n\t	vmovaps	 0x40(%%rax),%%ymm2			\n\t"	/* x = Re part of 1st base-root quartet */\
		"vmovaps	 0x20(%%rax),%%ymm1		\n\t	vmovaps	 0x60(%%rax),%%ymm3			\n\t"	/* y = Im part */\
		"vmovaps	%%ymm0,%%ymm4			\n\t	vmovaps		%%ymm2,%%ymm6			\n\t"	/* Copy x */\
		"vmovaps	%%ymm1,%%ymm5			\n\t	vmovaps		%%ymm3,%%ymm7			\n\t"	/* Copy y */\
		"vmulpd		%%ymm10,%%ymm0,%%ymm0	\n\t	vmulpd		%%ymm10,%%ymm2,%%ymm2	\n\t"	/* c.x */\
		"vmulpd		%%ymm11,%%ymm5,%%ymm5	\n\t	vmulpd		%%ymm11,%%ymm7,%%ymm7	\n\t"	/* s.y */\
		"vmulpd		%%ymm10,%%ymm1,%%ymm1	\n\t	vmulpd		%%ymm10,%%ymm3,%%ymm3	\n\t"	/* c.y */\
		"vmulpd		%%ymm11,%%ymm4,%%ymm4	\n\t	vmulpd		%%ymm11,%%ymm6,%%ymm6	\n\t"	/* s.x */\
		"vsubpd		%%ymm5 ,%%ymm0,%%ymm0	\n\t	vsubpd		%%ymm7 ,%%ymm2,%%ymm2	\n\t"	/* Out.re = c.x - s.y */\
		"vaddpd		%%ymm4 ,%%ymm1,%%ymm1	\n\t	vaddpd		%%ymm6 ,%%ymm3,%%ymm3	\n\t"	/* Out.im = c.y + s.x */\
		"vmovaps	%%ymm0 ,    (%%rax)		\n\t	vmovaps		%%ymm2 ,0x40(%%rax)		\n\t"	/* Store result, overwriting input base root */\
		"vmovaps	%%ymm1 ,0x20(%%rax)		\n\t	vmovaps		%%ymm3 ,0x60(%%rax)		\n\t"	/* Im part */\
		"/* Process next 2 base-root quartets: */"\
		"vmovaps	 0x80(%%rax),%%ymm0		\n\t	vmovaps	 0xc0(%%rax),%%ymm2			\n\t"	\
		"vmovaps	 0xa0(%%rax),%%ymm1		\n\t	vmovaps	 0xe0(%%rax),%%ymm3			\n\t"	\
		"vmovaps	%%ymm0,%%ymm4			\n\t	vmovaps		%%ymm2,%%ymm6			\n\t"	\
		"vmovaps	%%ymm1,%%ymm5			\n\t	vmovaps		%%ymm3,%%ymm7			\n\t"	\
		"vmulpd		%%ymm10,%%ymm0,%%ymm0	\n\t	vmulpd		%%ymm10,%%ymm2,%%ymm2	\n\t"	\
		"vmulpd		%%ymm11,%%ymm5,%%ymm5	\n\t	vmulpd		%%ymm11,%%ymm7,%%ymm7	\n\t"	\
		"vmulpd		%%ymm10,%%ymm1,%%ymm1	\n\t	vmulpd		%%ymm10,%%ymm3,%%ymm3	\n\t"	\
		"vmulpd		%%ymm11,%%ymm4,%%ymm4	\n\t	vmulpd		%%ymm11,%%ymm6,%%ymm6	\n\t"	\
		"vsubpd		%%ymm5 ,%%ymm0,%%ymm0	\n\t	vsubpd		%%ymm7 ,%%ymm2,%%ymm2	\n\t"	\
		"vaddpd		%%ymm4 ,%%ymm1,%%ymm1	\n\t	vaddpd		%%ymm6 ,%%ymm3,%%ymm3	\n\t"	\
		"vmovaps	%%ymm0 ,0x80(%%rax)		\n\t	vmovaps		%%ymm2 ,0xc0(%%rax)		\n\t"	/* Store result, overwriting input base root */\
		"vmovaps	%%ymm1 ,0xa0(%%rax)		\n\t	vmovaps		%%ymm3 ,0xe0(%%rax)		\n\t"	/* Im part */\
	:						/* outputs: none */\
	:	[__base_root]	"m" (Xbase_root)	/* All inputs from memory addresses here */\
		: "cc","memory","rax","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm10","xmm11"   /* Clobbered registers */\
	);\
	}

	/* Non-power-of-2-runlength Fermat-mod acyclic-transform/IBDWT carry macro.
	The array indices i/j/k/lcycle declared int in caller but assumed to have been byte-shift-converted at time this macro called,
	thus can use as complex-address-offsets.  Use bytewise literal offsets to save registers for several args here,as vvv-marked:
												                             vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
	#define SSE2_fermat_carry_norm_errcheck_X4_loacc(Xdata,Xbase_root,Xcy_re,Xcy_im,Xodd_radix,Xodd_radm2,Xodd_radm3,Xhalf_arr,Xsign_mask,XicycleA,XicycleB,XicycleC,XicycleD, XjcycleA,XkcycleA,XlcycleA, Xadd0,Xp1,Xp2,Xp3)\
	{\
	__asm__ volatile (\
		"movq		%[__add0],%%rcx		\n\t"/* base address for 4 prefetches-from-main-data-array spread through this macro */\
		"movq		%[__data],%%rax		\n\t"\
		"movq		%[__half_arr],%%rdx		\n\t"\
		"addq		$%c[__odd_radix],%%rdx				\n\t"	/* wt|wtinv|base|baseinv data offset by icycle array slots from resp. base addresses */\
		/* Multiply complex transform outputs [x,y] = [re,im] by inverse IBDWT weights, which include the 2/n scale factor: */\
		"movslq		%[__icycleA],%%r15		\n\t"	\
		"movslq		%[__icycleB],%%r9 		\n\t"	\
		"movslq		%[__icycleC],%%r8 		\n\t"	\
		"movslq		%[__icycleD],%%r10		\n\t"	\
		"vmovaps	(%%rdx,%%r15),%%ymm10	\n\t"	/* [wtinv0-3]A */\
		"vmovaps	(%%rdx,%%r9 ),%%ymm11	\n\t"	/* [wtinv0-3]B */\
		"vmovaps	(%%rdx,%%r8 ),%%ymm12	\n\t"	/* [wtinv0-3]C */\
		"vmovaps	(%%rdx,%%r10),%%ymm13	\n\t"	/* [wtinv0-3]D */\
	/* 4-way transpose of inputs (Re, Im parts separately): Inputs from r0/1,2/3,4/5.6/7. Outputs into ymm0-7: */\
	/* Real parts use ymm0,2,4,6, ymm8 as tmp-reg:					Imag parts use ymm1,3,5,7, ymm9 as tm-reg: */\
		"vmovaps	     (%%rax),%%ymm4						\n\t		vmovaps	0x020(%%rax),%%ymm5							\n\t"\
		"vmovaps	0x040(%%rax),%%ymm2						\n\t		vmovaps	0x060(%%rax),%%ymm3							\n\t"\
		"vmulpd		%%ymm10,%%ymm4,%%ymm4					\n\t		vmulpd		%%ymm10,%%ymm5,%%ymm5					\n\t"\
		"vmulpd		%%ymm11,%%ymm2,%%ymm2					\n\t		vmulpd		%%ymm11,%%ymm3,%%ymm3					\n\t"\
		"vshufpd	$15,%%ymm2,%%ymm4,%%ymm6				\n\t		vshufpd	$15,%%ymm3,%%ymm5,%%ymm7					\n\t"\
		"vshufpd	$0 ,%%ymm2,%%ymm4,%%ymm4				\n\t		vshufpd	$0 ,%%ymm3,%%ymm5,%%ymm5					\n\t"\
		"vmovaps	0x080(%%rax),%%ymm8						\n\t		vmovaps	0x0a0(%%rax),%%ymm9							\n\t"\
		"vmovaps	0x0c0(%%rax),%%ymm2						\n\t		vmovaps	0x0e0(%%rax),%%ymm3							\n\t"\
		"vmulpd		%%ymm12,%%ymm8,%%ymm8					\n\t		vmulpd		%%ymm12,%%ymm9,%%ymm9					\n\t"\
		"vmulpd		%%ymm13,%%ymm2,%%ymm2					\n\t		vmulpd		%%ymm13,%%ymm3,%%ymm3					\n\t"\
		"vshufpd	$15,%%ymm2,%%ymm8,%%ymm0				\n\t		vshufpd	$15,%%ymm3,%%ymm9,%%ymm1					\n\t"\
		"vshufpd	$0 ,%%ymm2,%%ymm8,%%ymm8				\n\t		vshufpd	$0 ,%%ymm3,%%ymm9,%%ymm9					\n\t"\
		"vperm2f128 $32,%%ymm0,%%ymm6,%%ymm2	/* Re B	*/	\n\t		vperm2f128 $32,%%ymm1,%%ymm7,%%ymm3		/* Im B	*/	\n\t"\
		"vperm2f128 $49,%%ymm0,%%ymm6,%%ymm6	/* Re D	*/	\n\t		vperm2f128 $49,%%ymm1,%%ymm7,%%ymm7		/* Im D	*/	\n\t"\
		"vperm2f128 $32,%%ymm8,%%ymm4,%%ymm0	/* Re A	*/	\n\t		vperm2f128 $32,%%ymm9,%%ymm5,%%ymm1 	/* Im A	*/	\n\t"\
		"vperm2f128 $49,%%ymm8,%%ymm4,%%ymm4	/* Re C	*/	\n\t		vperm2f128 $49,%%ymm9,%%ymm5,%%ymm5 	/* Im C	*/	\n\t"\
		"subq		$%c[__odd_radix],%%rdx				\n\t"\
		/* Base negacyclic roots at this address in [0,2,4,6]*0x20 (Re parts), [1,3,5,7]*0x20 (Imag parts) */\
		"movq		%[__sign_mask],%%rax	\n\t"\
		"vmovaps	(%%rax),%%ymm15	\n\t"	/* ymm15 free for rest of way; use to store sign_mask needed for floating ABS */\
		"movq	%[__base_root] ,%%rax		\n\t"	/* Won't need main-array again until output transpose, so re-use rax for base_root */\
		/* half_arr[0,1,2,3] = [base*2, baseinv*2,wt_re*2,wt_im*2] */\
		"vmovaps	-0x40(%%rdx),%%ymm13	\n\t"	/* XMM13 = maxerr */\
		"addq		%%r15,%%rdx				\n\t"	/* icycle assumed already in left-shifted ptr-byte-offset form */\
		/* Get next set of negacyclic roots: */\
		"vmovaps	    (%%rax),%%ymm10		\n\t"	/* c = Re part of 1st base-root quartet */\
		"vmovaps	0x20(%%rax),%%ymm11		\n\t"	/* s = Im part */\
	/* Do a-quartet: Data in ymm0,ymm1 */\
	"prefetcht0	(%%rcx)	\n\t"\
		"vmovaps	%%ymm13,%%ymm14			\n\t"	/* maxerr copy */\
		"movq		%[__cy_re],%%rbx		\n\t"\
		"vmovaps	%%ymm0,%%ymm8			\n\t"	/* x copy */\
		"vmovaps	%%ymm1,%%ymm9			\n\t"	/* y copy */\
		/* Inverse negacyclic weight is (wt_re, -wt_im): */\
		"vmulpd		%%ymm10,%%ymm0,%%ymm0	\n\t"	/* wt_re*[x     ] */\
		"vmulpd		%%ymm11,%%ymm9,%%ymm9	\n\t"	/* wt_im*[y copy] */\
		"vmulpd		%%ymm10,%%ymm1,%%ymm1	\n\t"	/* wt_re*[y     ] */\
		"vmulpd		%%ymm11,%%ymm8,%%ymm8	\n\t"	/* wt_im*[x copy] */\
		"vaddpd		%%ymm9 ,%%ymm0,%%ymm0	\n\t"	/* [a0-3.re] = x*wt_re + y*wt_im */\
		"vsubpd		%%ymm8 ,%%ymm1,%%ymm1	\n\t"	/* [a0-3.im] = y*wt_re - x*wt_im */\
		/* normalize a-pair, compute carryout, compute ROE: */\
		"vaddpd		           (%%rbx),%%ymm0,%%ymm0	\n\t"	/* [a0-3.re] + [cy0-3.re] */\
		"vaddpd		%c[__cy_im](%%rbx),%%ymm1,%%ymm1	\n\t"	/* [a0-3.im] + [cy0-3.im] */\
		"vmovaps	%%ymm0,%%ymm8			\n\t	vmovaps		%%ymm1,%%ymm9		\n\t"	/* copy x|y */\
		"vroundpd	$0,%%ymm0,%%ymm0		\n\t	vroundpd	$0,%%ymm1,%%ymm1	\n\t"	/* temp = DNINT(x|y) */\
		"vmovaps	%c[__odd_radm3](%%rdx),%%ymm12	\n\t"	/* [baseinv0-3] */\
		"vsubpd		%%ymm0 ,%%ymm8 ,%%ymm8 	\n\t	vsubpd		%%ymm1 ,%%ymm9 ,%%ymm9 	\n\t"	/* frac = [x - temp] */\
		"vandpd		%%ymm15,%%ymm8 ,%%ymm8 	\n\t	vandpd		%%ymm15,%%ymm9 ,%%ymm9 	\n\t"	/* frac = fabs(frac) */\
		"vmaxpd		%%ymm13,%%ymm8 ,%%ymm13	\n\t	vmaxpd		%%ymm14,%%ymm9 ,%%ymm14	\n\t"	/* if(frac > maxerr) maxerr=frac */\
		"vmovaps	%%ymm0 ,%%ymm8			\n\t	vmovaps		%%ymm1 ,%%ymm9			\n\t"	/* cpy temp */\
		"vmulpd		%%ymm12,%%ymm8 ,%%ymm8 	\n\t	vmulpd		%%ymm12,%%ymm9 ,%%ymm9 	\n\t"	/* temp*baseinv[0] ... inline the remaining +odd_radix offset in addressing */\
		"vmaxpd		%%ymm13,%%ymm14,%%ymm14	\n\t"	/* merge re|im maxerr vectors */\
		"vroundpd	$0,%%ymm8 ,%%ymm8 		\n\t	vroundpd	$0,%%ymm9,%%ymm9		\n\t"	/* [cy0-3.re] = DNINT(temp*baseinv[0]) */\
		"vmovaps	%c[__odd_radm2](%%rdx),%%ymm12	\n\t"	/* [base0-3] */\
		"vmovaps	%%ymm8,(%%rbx)			\n\t	vmovaps		%%ymm9,%c[__cy_im](%%rbx)\n\t"	/* store [cy0-3.re|im] */\
		"vmulpd		%%ymm12,%%ymm8,%%ymm8	\n\t	vmulpd		%%ymm12,%%ymm9,%%ymm9	\n\t"	/* base[0]*[cy0-3.re|im] */\
		"vsubpd		%%ymm8 ,%%ymm0,%%ymm0	\n\t	vsubpd		%%ymm9 ,%%ymm1,%%ymm1	\n\t"	/* XMM0|1 = [a0-3.re|im] = temp - [cy0-3.re|im]*base[0] */\
		"vmovaps	%%ymm0 ,%%ymm8			\n\t	vmovaps		%%ymm1 ,%%ymm9			\n\t"	/* cpy x|y */\
		/* Forward acyclic-convo weight is (wt_re, +wt_im): */\
		"vmulpd		%%ymm10,%%ymm0,%%ymm0	\n\t"	/* wt_re*[x     ] */\
		"vmulpd		%%ymm11,%%ymm9,%%ymm9	\n\t"	/* wt_im*[y copy] */\
		"vmulpd		%%ymm10,%%ymm1,%%ymm1	\n\t"	/* wt_re*[y     ] */\
		"vmulpd		%%ymm11,%%ymm8,%%ymm8	\n\t"	/* wt_im*[x copy] */\
		"vsubpd		%%ymm9 ,%%ymm0,%%ymm0	\n\t"	/* [a0-3.re] = x*wt_re - y*wt_im */\
		"vaddpd		%%ymm8 ,%%ymm1,%%ymm1	\n\t"	/* [a0-3.im] = y*wt_re + x*wt_im */\
		/* Up-multiply negacyclic roots stored in ymm10,11 by exp(j*I*Pi/2)/RADIX, for j = 4 */\
		"vmovaps	0x140(%%rax),%%ymm8 	\n\t"	/* x = Re(exp) in ymm10 */\
		"vmovaps	0x160(%%rax),%%ymm9 	\n\t"	/* y = Im(exp) in ymm11 */\
		"vmulpd		%%ymm10,%%ymm9 ,%%ymm12	\n\t"	/* ymm12 = c.y */\
		"vmulpd		%%ymm10,%%ymm8 ,%%ymm10	\n\t"	/* ymm10 = c.x */\
		"vmulpd		%%ymm11,%%ymm8 ,%%ymm8 	\n\t"	/* ymm8  = s.x */\
		"vmulpd		%%ymm11,%%ymm9 ,%%ymm9 	\n\t"	/* ymm9  = s.y */\
		"vsubpd		%%ymm9 ,%%ymm10,%%ymm10	\n\t"	/* ymm10 = wt.re = c.x - s.y; ymm9  free */\
		"vaddpd		%%ymm8 ,%%ymm12,%%ymm11	\n\t"	/* ymm11 = wt.im = s.x + c.y; ymm8 ,4 free */\
		"vmovaps	%%ymm10,    (%%rax)		\n\t"	/* Store result, overwriting the old base root */\
		"vmovaps	%%ymm11,0x20(%%rax)		\n\t"	/* Im part */\
		/* Get next set of negacyclic roots: */\
		"vmovaps	0x40(%%rax),%%ymm10		\n\t"	/* c = Re part of 2nd base-root quartet */\
		"vmovaps	0x60(%%rax),%%ymm11		\n\t"	/* s = Im part */\
	/* Now do b-quartet: Data in ymm2,ymm3 */\
	"movslq		%[__p1],%%r15	\n\t"	\
	"prefetcht0	(%%rcx,%%r15,8)	\n\t"\
		"movq		%[__half_arr],%%rdx		\n\t"\
		"vmovaps	%%ymm14,%%ymm13			\n\t"	/* maxerr copy */\
		"movslq		%[__jcycleA],%%r15		\n\t"\
		"addq		%%r15,%%rdx				\n\t"	/* jcycle assumed already in left-shifted ptr-byte-offset form */\
		"vmovaps	%%ymm2,%%ymm8			\n\t"	/* x copy */\
		"vmovaps	%%ymm3,%%ymm9			\n\t"	/* y copy */\
		/* Inverse negacyclic weight is (wt_re, -wt_im): */\
		"vmulpd		%%ymm10,%%ymm2,%%ymm2	\n\t"	/* wt_re*[x     ] */\
		"vmulpd		%%ymm11,%%ymm9,%%ymm9	\n\t"	/* wt_im*[y copy] */\
		"vmulpd		%%ymm10,%%ymm3,%%ymm3	\n\t"	/* wt_re*[y     ] */\
		"vmulpd		%%ymm11,%%ymm8,%%ymm8	\n\t"	/* wt_im*[x copy] */\
		"vaddpd		%%ymm9 ,%%ymm2,%%ymm2	\n\t"	/* [a0-3.re] = x*wt_re + y*wt_im */\
		"vsubpd		%%ymm8 ,%%ymm3,%%ymm3	\n\t"	/* [a0-3.im] = y*wt_re - x*wt_im */\
		/* normalize, compute carryout, compute ROE: */\
		"vaddpd		           (%%rbx),%%ymm2,%%ymm2	\n\t"	/* [a0-3.re] + [cy0-3.re] */\
		"vaddpd		%c[__cy_im](%%rbx),%%ymm3,%%ymm3	\n\t"	/* [a0-3.im] + [cy0-3.im] */\
		"vmovaps	%%ymm2,%%ymm8			\n\t	vmovaps		%%ymm3,%%ymm9		\n\t"	/* copy x|y */\
		"vroundpd	$0,%%ymm2,%%ymm2		\n\t	vroundpd	$0,%%ymm3,%%ymm3	\n\t"	/* temp = DNINT(x|y) */\
		"vmovaps	%c[__odd_radm3](%%rdx),%%ymm12	\n\t"	/* [baseinv0-3] */\
		"vsubpd		%%ymm2 ,%%ymm8 ,%%ymm8 	\n\t	vsubpd		%%ymm3 ,%%ymm9 ,%%ymm9 	\n\t"	/* frac = [x - temp] */\
		"vandpd		%%ymm15,%%ymm8 ,%%ymm8 	\n\t	vandpd		%%ymm15,%%ymm9 ,%%ymm9 	\n\t"	/* frac = fabs(frac) */\
		"vmaxpd		%%ymm13,%%ymm8 ,%%ymm13	\n\t	vmaxpd		%%ymm14,%%ymm9 ,%%ymm14	\n\t"	/* if(frac > maxerr) maxerr=frac */\
		"vmovaps	%%ymm2 ,%%ymm8			\n\t	vmovaps		%%ymm3 ,%%ymm9			\n\t"	/* cpy temp */\
		"vmulpd		%%ymm12,%%ymm8 ,%%ymm8 	\n\t	vmulpd		%%ymm12,%%ymm9 ,%%ymm9 	\n\t"	/* temp*baseinv[0] ... inline the remaining +odd_radix offset in addressing */\
		"vmaxpd		%%ymm13,%%ymm14,%%ymm14	\n\t"	/* merge re|im maxerr vectors */\
		"vroundpd	$0,%%ymm8 ,%%ymm8 		\n\t	vroundpd	$0,%%ymm9,%%ymm9		\n\t"	/* [cy0-3.re] = DNINT(temp*baseinv[0]) */\
		"vmovaps	%c[__odd_radm2](%%rdx),%%ymm12	\n\t"	/* [base0-3] */\
		"vmovaps	%%ymm8,(%%rbx)			\n\t	vmovaps		%%ymm9,%c[__cy_im](%%rbx)\n\t"	/* store [cy0-3.re|im] */\
		"vmulpd		%%ymm12,%%ymm8,%%ymm8	\n\t	vmulpd		%%ymm12,%%ymm9,%%ymm9	\n\t"	/* base[0]*[cy0-3.re|im] */\
		"vsubpd		%%ymm8 ,%%ymm2,%%ymm2	\n\t	vsubpd		%%ymm9 ,%%ymm3,%%ymm3	\n\t"	/* XMM0|1 = [a0-3.re|im] = temp - [cy0-3.re|im]*base[0] */\
		"vmovaps	%%ymm2 ,%%ymm8			\n\t	vmovaps		%%ymm3 ,%%ymm9			\n\t"	/* cpy x|y */\
		/* Forward acyclic-convo weight is (wt_re, +wt_im): */\
		"vmulpd		%%ymm10,%%ymm2,%%ymm2	\n\t"	/* wt_re*[x     ] */\
		"vmulpd		%%ymm11,%%ymm9,%%ymm9	\n\t"	/* wt_im*[y copy] */\
		"vmulpd		%%ymm10,%%ymm3,%%ymm3	\n\t"	/* wt_re*[y     ] */\
		"vmulpd		%%ymm11,%%ymm8,%%ymm8	\n\t"	/* wt_im*[x copy] */\
		"vsubpd		%%ymm9 ,%%ymm2,%%ymm2	\n\t"	/* [a0-3.re] = x*wt_re - y*wt_im */\
		"vaddpd		%%ymm8 ,%%ymm3,%%ymm3	\n\t"	/* [a0-3.im] = y*wt_re + x*wt_im */\
		/* Up-multiply negacyclic roots stored in ymm10,11 by exp(j*I*Pi/2)/RADIX, for j = 4 */\
		"vmovaps	0x140(%%rax),%%ymm8 	\n\t"	/* x = Re(exp) in ymm10 */\
		"vmovaps	0x160(%%rax),%%ymm9 	\n\t"	/* y = Im(exp) in ymm11 */\
		"vmulpd		%%ymm10,%%ymm9 ,%%ymm12	\n\t"	/* ymm12 = c.y */\
		"vmulpd		%%ymm10,%%ymm8 ,%%ymm10	\n\t"	/* ymm10 = c.x */\
		"vmulpd		%%ymm11,%%ymm8 ,%%ymm8 	\n\t"	/* ymm8  = s.x */\
		"vmulpd		%%ymm11,%%ymm9 ,%%ymm9 	\n\t"	/* ymm9  = s.y */\
		"vsubpd		%%ymm9 ,%%ymm10,%%ymm10	\n\t"	/* ymm10 = wt.re = c.x - s.y; ymm9  free */\
		"vaddpd		%%ymm8 ,%%ymm12,%%ymm11	\n\t"	/* ymm11 = wt.im = s.x + c.y; ymm8 ,4 free */\
		"vmovaps	%%ymm10,0x40(%%rax)		\n\t"	/* Store result, overwriting the old base root */\
		"vmovaps	%%ymm11,0x60(%%rax)		\n\t"	/* Im part */\
		/* Get next set of negacyclic roots: */\
		"vmovaps	0x80(%%rax),%%ymm10		\n\t"	/* c = Re part of 3rd base-root quartet */\
		"vmovaps	0xa0(%%rax),%%ymm11		\n\t"	/* s = Im part */\
	/* Now do c-quartet: Data in ymm4,ymm5 */\
	"movslq		%[__p2],%%r15	\n\t"	\
	"prefetcht0	(%%rcx,%%r15,8)	\n\t"\
		"movq		%[__half_arr],%%rdx		\n\t"\
		"vmovaps	%%ymm14,%%ymm13			\n\t"	/* maxerr copy */\
		"movslq		%[__kcycleA],%%r15		\n\t"\
		"addq		%%r15,%%rdx				\n\t"	/* kcycle assumed already in left-shifted ptr-byte-offset form */\
		"vmovaps	%%ymm4,%%ymm8			\n\t"	/* x copy */\
		"vmovaps	%%ymm5,%%ymm9			\n\t"	/* y copy */\
		/* Inverse negacyclic weight is (wt_re, -wt_im): */\
		"vmulpd		%%ymm10,%%ymm4,%%ymm4	\n\t"	/* wt_re*[x     ] */\
		"vmulpd		%%ymm11,%%ymm9,%%ymm9	\n\t"	/* wt_im*[y copy] */\
		"vmulpd		%%ymm10,%%ymm5,%%ymm5	\n\t"	/* wt_re*[y     ] */\
		"vmulpd		%%ymm11,%%ymm8,%%ymm8	\n\t"	/* wt_im*[x copy] */\
		"vaddpd		%%ymm9 ,%%ymm4,%%ymm4	\n\t"	/* [a0-3.re] = x*wt_re + y*wt_im */\
		"vsubpd		%%ymm8 ,%%ymm5,%%ymm5	\n\t"	/* [a0-3.im] = y*wt_re - x*wt_im */\
		/* normalize, compute carryout, compute ROE: */\
		"vaddpd		           (%%rbx),%%ymm4,%%ymm4	\n\t"	/* [a0-3.re] + [cy0-3.re] */\
		"vaddpd		%c[__cy_im](%%rbx),%%ymm5,%%ymm5	\n\t"	/* [a0-3.im] + [cy0-3.im] */\
		"vmovaps	%%ymm4,%%ymm8			\n\t	vmovaps		%%ymm5,%%ymm9		\n\t"	/* copy x|y */\
		"vroundpd	$0,%%ymm4,%%ymm4		\n\t	vroundpd	$0,%%ymm5,%%ymm5	\n\t"	/* temp = DNINT(x|y) */\
		"vmovaps	%c[__odd_radm3](%%rdx),%%ymm12	\n\t"	/* [baseinv0-3] */\
		"vsubpd		%%ymm4 ,%%ymm8 ,%%ymm8 	\n\t	vsubpd		%%ymm5 ,%%ymm9 ,%%ymm9 	\n\t"	/* frac = [x - temp] */\
		"vandpd		%%ymm15,%%ymm8 ,%%ymm8 	\n\t	vandpd		%%ymm15,%%ymm9 ,%%ymm9 	\n\t"	/* frac = fabs(frac) */\
		"vmaxpd		%%ymm13,%%ymm8 ,%%ymm13	\n\t	vmaxpd		%%ymm14,%%ymm9 ,%%ymm14	\n\t"	/* if(frac > maxerr) maxerr=frac */\
		"vmovaps	%%ymm4 ,%%ymm8			\n\t	vmovaps		%%ymm5 ,%%ymm9			\n\t"	/* cpy temp */\
		"vmulpd		%%ymm12,%%ymm8 ,%%ymm8 	\n\t	vmulpd		%%ymm12,%%ymm9 ,%%ymm9 	\n\t"	/* temp*baseinv[0] ... inline the remaining +odd_radix offset in addressing */\
		"vmaxpd		%%ymm13,%%ymm14,%%ymm14	\n\t"	/* merge re|im maxerr vectors */\
		"vroundpd	$0,%%ymm8 ,%%ymm8 		\n\t	vroundpd	$0,%%ymm9,%%ymm9		\n\t"	/* [cy0-3.re] = DNINT(temp*baseinv[0]) */\
		"vmovaps	%c[__odd_radm2](%%rdx),%%ymm12	\n\t"	/* [base0-3] */\
		"vmovaps	%%ymm8,(%%rbx)			\n\t	vmovaps		%%ymm9,%c[__cy_im](%%rbx)\n\t"	/* store [cy0-3.re|im] */\
		"vmulpd		%%ymm12,%%ymm8,%%ymm8	\n\t	vmulpd		%%ymm12,%%ymm9,%%ymm9	\n\t"	/* base[0]*[cy0-3.re|im] */\
		"vsubpd		%%ymm8 ,%%ymm4,%%ymm4	\n\t	vsubpd		%%ymm9 ,%%ymm5,%%ymm5	\n\t"	/* XMM0|1 = [a0-3.re|im] = temp - [cy0-3.re|im]*base[0] */\
		"vmovaps	%%ymm4 ,%%ymm8			\n\t	vmovaps		%%ymm5 ,%%ymm9			\n\t"	/* cpy x|y */\
		/* Forward acyclic-convo weight is (wt_re, +wt_im): */\
		"vmulpd		%%ymm10,%%ymm4,%%ymm4	\n\t"	/* wt_re*[x     ] */\
		"vmulpd		%%ymm11,%%ymm9,%%ymm9	\n\t"	/* wt_im*[y copy] */\
		"vmulpd		%%ymm10,%%ymm5,%%ymm5	\n\t"	/* wt_re*[y     ] */\
		"vmulpd		%%ymm11,%%ymm8,%%ymm8	\n\t"	/* wt_im*[x copy] */\
		"vsubpd		%%ymm9 ,%%ymm4,%%ymm4	\n\t"	/* [a0-3.re] = x*wt_re - y*wt_im */\
		"vaddpd		%%ymm8 ,%%ymm5,%%ymm5	\n\t"	/* [a0-3.im] = y*wt_re + x*wt_im */\
		/* Up-multiply negacyclic roots stored in ymm10,11 by exp(j*I*Pi/2)/RADIX, for j = 4 */\
		"vmovaps	0x140(%%rax),%%ymm8 	\n\t"	/* x = Re(exp) in ymm10 */\
		"vmovaps	0x160(%%rax),%%ymm9 	\n\t"	/* y = Im(exp) in ymm11 */\
		"vmulpd		%%ymm10,%%ymm9 ,%%ymm12	\n\t"	/* ymm12 = c.y */\
		"vmulpd		%%ymm10,%%ymm8 ,%%ymm10	\n\t"	/* ymm10 = c.x */\
		"vmulpd		%%ymm11,%%ymm8 ,%%ymm8 	\n\t"	/* ymm8  = s.x */\
		"vmulpd		%%ymm11,%%ymm9 ,%%ymm9 	\n\t"	/* ymm9  = s.y */\
		"vsubpd		%%ymm9 ,%%ymm10,%%ymm10	\n\t"	/* ymm10 = wt.re = c.x - s.y; ymm9  free */\
		"vaddpd		%%ymm8 ,%%ymm12,%%ymm11	\n\t"	/* ymm11 = wt.im = s.x + c.y; ymm8 ,4 free */\
		"vmovaps	%%ymm10,0x80(%%rax)		\n\t"	/* Store result, overwriting the old base root */\
		"vmovaps	%%ymm11,0xa0(%%rax)		\n\t"	/* Im part */\
		/* Get next set of negacyclic roots: */\
		"vmovaps	0xc0(%%rax),%%ymm10		\n\t"	/* c = Re part of 4th base-root quartet */\
		"vmovaps	0xe0(%%rax),%%ymm11		\n\t"	/* s = Im part */\
	/* Now do d-quartet: Data in ymm6,ymm7 */\
	"movslq		%[__p3],%%r15	\n\t"	\
	"prefetcht0	(%%rcx,%%r15,8)	\n\t"\
		"movq		%[__half_arr],%%rdx		\n\t"\
		"vmovaps	%%ymm14,%%ymm13			\n\t"	/* maxerr copy */\
		"movslq		%[__lcycleA],%%r15		\n\t"\
		"addq		%%r15,%%rdx				\n\t"	/* lcycle assumed already in left-shifted ptr-byte-offset form */\
		"vmovaps	%%ymm6,%%ymm8			\n\t"	/* x copy */\
		"vmovaps	%%ymm7,%%ymm9			\n\t"	/* y copy */\
		/* Inverse negacyclic weight is (wt_re, -wt_im): */\
		"vmulpd		%%ymm10,%%ymm6,%%ymm6	\n\t"	/* wt_re*[x     ] */\
		"vmulpd		%%ymm11,%%ymm9,%%ymm9	\n\t"	/* wt_im*[y copy] */\
		"vmulpd		%%ymm10,%%ymm7,%%ymm7	\n\t"	/* wt_re*[y     ] */\
		"vmulpd		%%ymm11,%%ymm8,%%ymm8	\n\t"	/* wt_im*[x copy] */\
		"vaddpd		%%ymm9 ,%%ymm6,%%ymm6	\n\t"	/* [a0-3.re] = x*wt_re + y*wt_im */\
		"vsubpd		%%ymm8 ,%%ymm7,%%ymm7	\n\t"	/* [a0-3.im] = y*wt_re - x*wt_im */\
		/* normalize, compute carryout, compute ROE: */\
		"vaddpd		           (%%rbx),%%ymm6,%%ymm6	\n\t"	/* [a0-3.re] + [cy0-3.re] */\
		"vaddpd		%c[__cy_im](%%rbx),%%ymm7,%%ymm7	\n\t"	/* [a0-3.im] + [cy0-3.im] */\
		"vmovaps	%%ymm6,%%ymm8			\n\t	vmovaps		%%ymm7,%%ymm9		\n\t"	/* copy x|y */\
		"vroundpd	$0,%%ymm6,%%ymm6		\n\t	vroundpd	$0,%%ymm7,%%ymm7	\n\t"	/* temp = DNINT(x|y) */\
		"vmovaps	%c[__odd_radm3](%%rdx),%%ymm12	\n\t"	/* [baseinv0-3] */\
		"vsubpd		%%ymm6 ,%%ymm8 ,%%ymm8 	\n\t	vsubpd		%%ymm7 ,%%ymm9 ,%%ymm9 	\n\t"	/* frac = [x - temp] */\
		"vandpd		%%ymm15,%%ymm8 ,%%ymm8 	\n\t	vandpd		%%ymm15,%%ymm9 ,%%ymm9 	\n\t"	/* frac = fabs(frac) */\
		"vmaxpd		%%ymm13,%%ymm8 ,%%ymm13	\n\t	vmaxpd		%%ymm14,%%ymm9 ,%%ymm14	\n\t"	/* if(frac > maxerr) maxerr=frac */\
		"vmovaps	%%ymm6 ,%%ymm8			\n\t	vmovaps		%%ymm7 ,%%ymm9			\n\t"	/* cpy temp */\
		"vmulpd		%%ymm12,%%ymm8 ,%%ymm8 	\n\t	vmulpd		%%ymm12,%%ymm9 ,%%ymm9 	\n\t"	/* temp*baseinv[0] ... inline the remaining +odd_radix offset in addressing */\
		"vmaxpd		%%ymm13,%%ymm14,%%ymm14	\n\t"	/* merge re|im maxerr vectors */\
		"vroundpd	$0,%%ymm8 ,%%ymm8 		\n\t	vroundpd	$0,%%ymm9,%%ymm9		\n\t"	/* [cy0-3.re] = DNINT(temp*baseinv[0]) */\
		"vmovaps	%c[__odd_radm2](%%rdx),%%ymm12	\n\t"	/* [base0-3] */\
		"vmovaps	%%ymm8,(%%rbx)			\n\t	vmovaps		%%ymm9,%c[__cy_im](%%rbx)\n\t"	/* store [cy0-3.re|im] */\
		"vmulpd		%%ymm12,%%ymm8,%%ymm8	\n\t	vmulpd		%%ymm12,%%ymm9,%%ymm9	\n\t"	/* base[0]*[cy0-3.re|im] */\
		"vsubpd		%%ymm8 ,%%ymm6,%%ymm6	\n\t	vsubpd		%%ymm9 ,%%ymm7,%%ymm7	\n\t"	/* XMM0|1 = [a0-3.re|im] = temp - [cy0-3.re|im]*base[0] */\
		"vmovaps	%%ymm6 ,%%ymm8			\n\t	vmovaps		%%ymm7 ,%%ymm9			\n\t"	/* cpy x|y */\
		/* Forward acyclic-convo weight is (wt_re, +wt_im): */\
		"vmulpd		%%ymm10,%%ymm6,%%ymm6	\n\t"	/* wt_re*[x     ] */\
		"vmulpd		%%ymm11,%%ymm9,%%ymm9	\n\t"	/* wt_im*[y copy] */\
		"vmulpd		%%ymm10,%%ymm7,%%ymm7	\n\t"	/* wt_re*[y     ] */\
		"vmulpd		%%ymm11,%%ymm8,%%ymm8	\n\t"	/* wt_im*[x copy] */\
		"vsubpd		%%ymm9 ,%%ymm6,%%ymm6	\n\t"	/* [a0-3.re] = x*wt_re - y*wt_im */\
		"vaddpd		%%ymm8 ,%%ymm7,%%ymm7	\n\t"	/* [a0-3.im] = y*wt_re + x*wt_im */\
		/* Up-multiply negacyclic roots stored in ymm10,11 by exp(j*I*Pi/2)/RADIX, for j = 4 */\
		"vmovaps	0x140(%%rax),%%ymm8 	\n\t"	/* x = Re(exp) in ymm10 */\
		"vmovaps	0x160(%%rax),%%ymm9 	\n\t"	/* y = Im(exp) in ymm11 */\
		"vmulpd		%%ymm10,%%ymm9 ,%%ymm12	\n\t"	/* ymm12 = c.y */\
		"vmulpd		%%ymm10,%%ymm8 ,%%ymm10	\n\t"	/* ymm10 = c.x */\
		"vmulpd		%%ymm11,%%ymm8 ,%%ymm8 	\n\t"	/* ymm8  = s.x */\
		"vmulpd		%%ymm11,%%ymm9 ,%%ymm9 	\n\t"	/* ymm9  = s.y */\
		"vsubpd		%%ymm9 ,%%ymm10,%%ymm10	\n\t"	/* ymm10 = wt.re = c.x - s.y; ymm9  free */\
		"vaddpd		%%ymm8 ,%%ymm12,%%ymm11	\n\t"	/* ymm11 = wt.im = s.x + c.y; ymm8 ,4 free */\
		"vmovaps	%%ymm10,0xc0(%%rax)		\n\t"	/* Store result, overwriting the old base root */\
		"vmovaps	%%ymm11,0xe0(%%rax)		\n\t"	/* Im part */\
		/* Store maxerr: */\
		"movq		%[__half_arr],%%rdx		\n\t"\
		"vmovaps	%%ymm14,-0x40(%%rdx)	\n\t"\
	/* 4-way transpose of outputs (Re, Im parts separately): Inputs from ymm0-7. Outputs into r0/1,2/3,4/5.6/7: */	\
	/* Because default inputs for our 4 x 4 transpose macro (e.g. the one used at start of this carry macro) */\
	/* are into ymm4/2/8/2, munge inputs into that order, resolving name-conflicts via use of the now-available ymm8-15 for outputs: */\
		"movq		%[__data],%%rax			\n\t"\
		"vshufpd	$15,%%ymm2,%%ymm0,%%ymm10					\n\t		vshufpd	$15,%%ymm3,%%ymm1,%%ymm11						\n\t"\
		"vshufpd	$0 ,%%ymm2,%%ymm0,%%ymm0					\n\t		vshufpd	$0 ,%%ymm3,%%ymm1,%%ymm1						\n\t"\
		"vshufpd	$15,%%ymm6,%%ymm4,%%ymm12					\n\t		vshufpd	$15,%%ymm7,%%ymm5,%%ymm13						\n\t"\
		"vshufpd	$0 ,%%ymm6,%%ymm4,%%ymm4					\n\t		vshufpd	$0 ,%%ymm7,%%ymm5,%%ymm5						\n\t"\
		"vperm2f128 $32,%%ymm12,%%ymm10,%%ymm2 		/* Re B	*/	\n\t		vperm2f128 $32,%%ymm13,%%ymm11,%%ymm3		/* Im B	*/	\n\t"\
		"vperm2f128 $49,%%ymm12,%%ymm10,%%ymm10		/* Re D	*/	\n\t		vperm2f128 $49,%%ymm13,%%ymm11,%%ymm11		/* Im D	*/	\n\t"\
		"vperm2f128 $32,%%ymm4 ,%%ymm0 ,%%ymm12		/* Re A	*/	\n\t		vperm2f128 $32,%%ymm5 ,%%ymm1 ,%%ymm13 		/* Im A	*/	\n\t"\
		"vperm2f128 $49,%%ymm4 ,%%ymm0 ,%%ymm0 		/* Re C	*/	\n\t		vperm2f128 $49,%%ymm5 ,%%ymm1 ,%%ymm1		/* Im C	*/	\n\t"\
		/* Multiply normalized, re-permuted transform outputs by forward IBDWT weights: */\
		"movslq		%[__icycleA],%%r15		\n\t"	\
		"vmovaps	(%%rdx,%%r15),%%ymm4	\n\t"	/* [wt0-3]A */\
		"vmovaps	(%%rdx,%%r9 ),%%ymm5	\n\t"	/* [wt0-3]B */\
		"vmovaps	(%%rdx,%%r8 ),%%ymm6	\n\t"	/* [wt0-3]C */\
		"vmovaps	(%%rdx,%%r10),%%ymm7	\n\t"	/* [wt0-3]D */\
		"vmulpd		%%ymm4,%%ymm12,%%ymm12						\n\t		vmulpd		%%ymm4,%%ymm13,%%ymm13			\n\t"\
		"vmulpd		%%ymm5,%%ymm2 ,%%ymm2 						\n\t		vmulpd		%%ymm5,%%ymm3 ,%%ymm3 			\n\t"\
		"vmulpd		%%ymm6,%%ymm0 ,%%ymm0 						\n\t		vmulpd		%%ymm6,%%ymm1 ,%%ymm1 			\n\t"\
		"vmulpd		%%ymm7,%%ymm10,%%ymm10						\n\t		vmulpd		%%ymm7,%%ymm11,%%ymm11			\n\t"\
		/* And write 'em back to memory: */\
		"vmovaps	%%ymm12,     (%%rax)						\n\t		vmovaps	%%ymm13,0x020(%%rax)				\n\t"\
		"vmovaps	%%ymm2 ,0x040(%%rax)						\n\t		vmovaps	%%ymm3 ,0x060(%%rax)				\n\t"\
		"vmovaps	%%ymm0 ,0x080(%%rax)						\n\t		vmovaps	%%ymm1 ,0x0a0(%%rax)				\n\t"\
		"vmovaps	%%ymm10,0x0c0(%%rax)						\n\t		vmovaps	%%ymm11,0x0e0(%%rax)				\n\t"\
		:						/* outputs: none */\
		:	[__data]		"m" (Xdata)	/* All inputs from memory addresses here */\
		,	[__base_root]	"m" (Xbase_root)\
		,	[__cy_re]		"m" (Xcy_re)\
		,	[__cy_im]		"e" (Xcy_im)	/* Use literal-byte-offset for this ome to save a reg */\
		/* [1,2,3]*odd_radix are the index offsets to the wtinv, base, and base_inv values, resp. - assumed << l2_sz_vd on input: */\
		,	[__odd_radix]   "e" (Xodd_radix)\
		,	[__odd_radm2]   "e" (Xodd_radm2)\
		,	[__odd_radm3]   "e" (Xodd_radm3)\
		,	[__half_arr]	"m" (Xhalf_arr)\
		,	[__sign_mask]	"m" (Xsign_mask)\
		/* Need quartet of ascending [modulo odd_radix] icycle indices for IBDWT weights: */\
		,	[__icycleA]		"m" (XicycleA)\
		,	[__icycleB]		"m" (XicycleB)\
		,	[__icycleC]		"m" (XicycleC)\
		,	[__icycleD]		"m" (XicycleD)\
		/* Need quartet of same-index [i,j,k,l]cycle indices for negacyclic weights and base/baseinv normalizations: */\
		,	[__jcycleA]		"m" (XjcycleA)\
		,	[__kcycleA]		"m" (XkcycleA)\
		,	[__lcycleA]		"m" (XlcycleA)\
		/* Prefetch: base address and 3 index offsets */\
		,	[__add0] "m" (Xadd0)\
		,	[__p1] "m" (Xp1)\
		,	[__p2] "m" (Xp2)\
		,	[__p3] "m" (Xp3)\
		: "cc","memory","rax","rbx","rcx","rdx","r8","r9","r10","r15","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14"   /* Clobbered registers */\
	);\
	}

  #endif	// AVX2/FMA3?

#else	// SSE2

	/* Power-of-2-runlength Fermat-mod acyclic-transform carry macro. (No IBDWT needed for power-of-2 runlenghts).
	*/
	#define SSE2_fermat_carry_norm_pow2_errcheck(Xdata,Xcy,Xnrt_bits,Xnrtm1,Xidx_offset,Xidx_incr,Xhalf_arr,Xsign_mask,Xadd1,Xadd2, Xadd0)\
	{\
	__asm__ volatile (\
	"movq	%[__add0],%%r14	\n\t"\
	"prefetcht0	(%%r14)		\n\t"\
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
		/* Prefetch address */\
		,	[__add0] "m" (Xadd0)\
		: "cc","memory","rax","rbx","rcx","rdx","rsi","rdi","r14","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"	/* Clobbered registers */\
	);\
	}

	/* Same power-of-2-transform deal as above, but use xmm8-15 to process 2 sets of carries side-by-side.
	Data/Carry #2 assumed offset by +0x20/0x10 from #1 (which are accessed via the [__data/__cy] pointers, resp.):
	*/
	#define SSE2_fermat_carry_norm_pow2_errcheck_X2(Xdata,Xcy,Xnrt_bits,Xnrtm1,Xidx_offset,Xidx_incr,Xhalf_arr,Xsign_mask,Xadd1,Xadd2, Xadd0,Xp1)\
	{\
	__asm__ volatile (\
	"movq	%[__add0],%%r14	\n\t"/* base address for 2 prefetches-from-main-data-array spread through this macro */\
	"prefetcht0	(%%r14)	\n\t"\
		/* lcol -> rcol index analogs: [rsi,rax,rbx] -> [r10,r11,r12], [rcx,rdx,rdi] shared */\
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
	"movslq		%[__p1],%%r15	\n\t"	\
	"prefetcht0	(%%r14,%%r15,8)	\n\t"\
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
		/* Prefetch: base address and 1 index offset */\
		,	[__add0] "m" (Xadd0)\
		,	[__p1] "m" (Xp1)\
		: "cc","memory","rax","rbx","rcx","rdx","rsi","rdi","r10","r11","r12","r14","r15","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14"		/* Clobbered registers */\
	);\
	}

	/* Non-power-of-2-runlength Fermat-mod acyclic-transform/IBDWT carry macro.
	The array indices icycle0/1 declared int in caller but assumed to have been << 4 at time this macro called, thus can use as complex-address-offsets.
	*/
	#define SSE2_fermat_carry_norm_errcheck(Xdata,Xcy,Xnrt_bits,Xnrtm1,Xidx_offset,Xidx_incr,Xodd_radix,Xhalf_arr,Xsign_mask,Xadd1,Xadd2,Xicycle0,Xjcycle0, Xadd0)\
	{\
	__asm__ volatile (\
	"movq	%[__add0],%%r14	\n\t"\
	"prefetcht0	(%%r14)	\n\t"\
		"movslq	%[__idx_offset],%%rsi	\n\t"	/* esi stores [j + idx_offset], idx_offset starts = 0, gets incremented by idx_incr each macro invocation */\
		"movslq %[__odd_radix],%%rdi	\n\t"	/* [1,2,3]*odd_radix are the index offsets to the wtinv, base, and base_inv values, respectively. */\
		"movslq	%[__nrt_bits],%%rcx		\n\t"	\
		"movslq %[__icycle0],%%r10		\n\t"	\
		"movslq	%[__jcycle0],%%r11		\n\t"	\
		"movq		%%rsi,%%rax			\n\t"	/* j + idx_offset */\
		"shrq		$1,%%rax			\n\t"	/* l = ((j + idx_offset) >> 1) */\
		"movq		%%rax,%%rbx			\n\t"	\
		"andq		%[__nrtm1],%%rax	\n\t"	/* k1 = (l & __NRTM1) */\
		"shrq		%%cl,%%rbx			\n\t"	/* k2=(l >> __NRT_BITS) */\
		"shlq		$4,%%rax			\n\t"	/* 16 bytes for array-of-complex */\
		"shlq		$4,%%rbx			\n\t"	/* 16 bytes for array-of-complex */\
		"shlq		$4,%%rdi			\n\t"	/* 16 bytes for array-of-complex */\
		"addq		%[__add1],%%rax		\n\t"	/* rn0[k1] */\
		"addq		%[__add2],%%rbx		\n\t"	/* rn1[k2] */\
		"movaps		(%%rax),%%xmm0		\n\t"	/* [c0,s0] */\
		"movaps		(%%rbx),%%xmm1		\n\t"	/* [x0,y0] */\
		"movq		%%rsi,%%rax			\n\t"	\
		"movaps		%%xmm1,%%xmm2		\n\t"	/* [x0,y0] copy */\
		"shufpd	$1,	%%xmm2,%%xmm2		\n\t"	/* [y0,x0] (swap re <--> im) */\
		"mulpd		%%xmm0,%%xmm1		\n\t"	/* [c0.x0,s0.y0] */\
		"mulpd		%%xmm0,%%xmm2		\n\t"	/* [c0.y0,s0.x0] 1,2 used */\
		/* Get next root for interleaving with the first: */\
		"addq		$2,%%rax			\n\t"	\
		"shrq		$1,%%rax			\n\t"	/* l = ((j + idx_offset) >> 1) */\
		"movq		%%rax,%%rbx			\n\t"	\
		"andq		%[__nrtm1],%%rax	\n\t"	/* k1 = (l & __NRTM1) */\
		"shrq		%%cl,%%rbx			\n\t"	/* k2=(l >> __NRT_BITS) */\
		"shlq		$4,%%rax			\n\t"	/* 16 bytes for array-of-complex */\
		"shlq		$4,%%rbx			\n\t"	/* 16 bytes for array-of-complex */\
		"addq		%[__add1],%%rax		\n\t"	/* rn0[k1] */\
		"addq		%[__add2],%%rbx		\n\t"	/* rn1[k2] */\
		"movaps		(%%rax),%%xmm0		\n\t"	/* [c1,s1] */\
		"movaps		(%%rbx),%%xmm3		\n\t"	/* [x1,y1] 0-3 used*/\
		"movq		%%rsi,%%rax			\n\t"	\
		"movaps		%%xmm3,%%xmm4		\n\t"	/* [x1,y1] copy */\
		"shufpd	$1,	%%xmm4,%%xmm4		\n\t"	/* [y1,x1] (swap re <--> im) */\
		"mulpd		%%xmm0,%%xmm3		\n\t"	/* [c1.x1,s1.y1] */\
		"mulpd		%%xmm0,%%xmm4		\n\t"	/* [c1.y1,s1.x1] 1-4 used */\
		"movaps		%%xmm1,%%xmm0		\n\t"	/* xmm0 <- copy [c0.x0,s0.y0] */\
		"unpcklpd	%%xmm3,%%xmm0		\n\t"	/* [c0.x0,c1.x1] */\
		"unpckhpd	%%xmm3,%%xmm1		\n\t"	/* [s0.y0,s1.y1], 0-2,4 used */\
		"subpd		%%xmm1,%%xmm0		\n\t"	/* XMM0 = [wt_r0,wt_r1] 0,2,4 used */\
		"movaps		%%xmm2,%%xmm1		\n\t"	/* xmm1 <- copy [c0.y0,s0.x0] 0-2,4 used */\
		"unpcklpd	%%xmm4,%%xmm1		\n\t"	/* [c0.y0,c1.y1] */\
		"unpckhpd	%%xmm4,%%xmm2		\n\t"	/* [s0.x0,s1.x1] */\
		"addpd		%%xmm2,%%xmm1		\n\t"	/* XMM1 = [wt_i0,wt_i1] 0-1 used */\
		/* half_arr[0,1,2,3] = [base*2, baseinv*2,wt_re*2,wt_im*2] */\
		"movq		%[__half_arr],%%rcx	\n\t"	/* No longer need __NRT_BITS, so reuse ecx */\
		/* Multiply the complex transform output [x,y] = [re,im] by the inverse IBDWT weight, which includes the scale factor: [x,y] *= wtinv: */\
		"movq		%[__data],%%rdx		\n\t"	\
		"movaps		     (%%rdx),%%xmm4	\n\t"	/* x = [a.re,b.re] */\
		"movaps		 0x10(%%rdx),%%xmm2	\n\t"	/* y = [a.im,b.im] */\
		"addq		%%r10,%%rcx			\n\t"	\
		"movaps	(%%rcx,%%rdi),%%xmm5	\n\t"	/* [wtinv0,wtinv1] */\
		"subq		%%r10,%%rcx			\n\t"	\
		"mulpd		%%xmm5,%%xmm4		\n\t"	\
		"mulpd		%%xmm5,%%xmm2		\n\t"	\
		"movaps		%%xmm4,%%xmm5		\n\t"	/* x copy */\
		"movaps		%%xmm2,%%xmm3		\n\t"	/* y copy */\
		/* Inverse weight is (wt_re, -wt_im): */\
		"mulpd		%%xmm0,%%xmm4		\n\t"	/* [x     ]*wt_re */\
		"mulpd		%%xmm1,%%xmm3		\n\t"	/* [y copy]*wt_im */\
		"mulpd		%%xmm0,%%xmm2		\n\t"	/* [y     ]*wt_re */\
		"mulpd		%%xmm1,%%xmm5		\n\t"	/* [x copy]*wt_im */\
		"addpd		%%xmm3,%%xmm4		\n\t"	/* [a.re,b.re] = x*wt_re + y*wt_im */\
		"subpd		%%xmm5,%%xmm2		\n\t"	/* [a.im,b.im] = y*wt_re - x*wt_im */\
		"movq		%[__cy],%%rbx		\n\t"	\
		"movaps		%%xmm4,%%xmm5		\n\t"	/* [a.re,b.re] copy */\
		"shufpd	$0,	%%xmm2,%%xmm4		\n\t"	/* XMM4 = x = [a.re,a.im] */\
		"shufpd	$3,	%%xmm2,%%xmm5		\n\t"	/* XMM5 = y = [b.re,b.im] 0,1,4,5 uaed */\
		/* normalize a-pair, compute carryout, compute ROE: */\
		"addpd		     (%%rbx),%%xmm4	\n\t"	/* [a.re,a.im] + [cx,cy] */\
		"movaps		-0x20(%%rcx),%%xmm6	\n\t"	/* XMM6 = maxerr */\
		"movaps		-0x10(%%rcx),%%xmm7	\n\t"	/* XMM7 = rnd_const */\
		"addq	   %%r10,%%rcx			\n\t"	\
		"movaps		%%xmm4,%%xmm2		\n\t"	/* copy x */\
		"shlq	   $1,%%rdi				\n\t"	\
		"addpd		%%xmm7,%%xmm4		\n\t"	\
		"subpd		%%xmm7,%%xmm4		\n\t"	/* temp = DNINT(x) */\
		"movq		%[__sign_mask],%%rax\n\t"	\
		"subpd		%%xmm4,%%xmm2		\n\t"	/* frac = [x - temp] */\
		"andpd		(%%rax),%%xmm2		\n\t"	/* frac = fabs(frac) */\
		"maxpd		%%xmm6,%%xmm2		\n\t"	/* if(frac > maxerr) maxerr=frac */\
		"movaps		%%xmm2,%%xmm6		\n\t"	/* Note serialization here! */\
		"addq		%%rdi,%%rcx			\n\t"	\
		"shrq		$1,%%rdi			\n\t"	\
		"movaps		%%xmm4,%%xmm2		\n\t"	/* cpy temp */\
		"mulpd  (%%rcx,%%rdi),%%xmm2	\n\t"	/* temp*baseinv[0] */\
		"addpd		%%xmm7,%%xmm2		\n\t"	\
		"subpd		%%xmm7,%%xmm2		\n\t"	/* [cx,cy] = DNINT(temp*baseinv[0]) */\
		"movaps		%%xmm2,%%xmm3		\n\t"	/* cpy [cx,cy] */\
		"mulpd		(%%rcx),%%xmm3		\n\t"	/* [cx,cy]*base[0] */\
		"subq		%%r10,%%rcx			\n\t"	\
		"subpd		%%xmm3,%%xmm4		\n\t"	/* XMM4 = [a.re,a.im] = temp-[cx,cy]*base[0] */\
		/* Now do b-pair: [b.re,b.im] in xmm5, carry in xmm2, xmm3 free, wt_[re,im] in xmmA,B, xmm6 free, rnd_const in xmm7: */\
		"addpd		%%xmm2,%%xmm5		\n\t"	/* [b.re,b.im] + [cx,cy] */\
		"movaps		%%xmm5,%%xmm2		\n\t"	/* copy y */\
		"addpd		%%xmm7,%%xmm5		\n\t"	\
		"subpd		%%xmm7,%%xmm5		\n\t"	/* temp = DNINT(y) */\
		"subpd		%%xmm5,%%xmm2		\n\t"	/* frac = [y - temp] */\
		"andpd		(%%rax),%%xmm2		\n\t"	/* frac = fabs(frac) */\
		"maxpd		%%xmm6,%%xmm2		\n\t"	/* if(frac > maxerr) maxerr=frac */\
		"movaps		%%xmm2,%%xmm6		\n\t"	/* Note serialization here! */\
		"movaps		%%xmm5,%%xmm2		\n\t"	/* cpy temp */\
		"addq		%%r11,%%rcx			\n\t"	\
		"mulpd  (%%rcx,%%rdi),%%xmm2	\n\t"	/* temp*baseinv[1] */\
		"addpd		%%xmm7,%%xmm2		\n\t"	\
		"subpd		%%xmm7,%%xmm2		\n\t"	/* [cx,cy] = DNINT(temp*baseinv[1]) */\
		"shlq		$1,%%rdi			\n\t"	/* prepare to re-subtract 2*odd_radix from local-store pointer */\
		"movaps		%%xmm2,%%xmm3		\n\t"	/* cpy [cx,cy] */\
		"mulpd		(%%rcx),%%xmm3		\n\t"	/* [cx,cy]*base[1] */\
		"subq		%%r11,%%rcx			\n\t"	\
		"subpd		%%xmm3,%%xmm5		\n\t"	/* XMM5 = [b.re,b.im] = temp-[cx,cy]*base[1] */\
		"movaps		%%xmm2,(%%rbx)		\n\t"	/* store cy_out */\
		"movaps		%%xmm4,%%xmm2		\n\t"	/* [a.re,a.im] copy */\
		"shufpd	$0,	%%xmm5,%%xmm4		\n\t"	/* x = [a.re,b.re] */\
		"shufpd	$3,	%%xmm5,%%xmm2		\n\t"	/* y = [a.im,b.im] */\
		"movaps		%%xmm4,%%xmm5		\n\t"	/* x copy */\
		"movaps		%%xmm2,%%xmm3		\n\t"	/* y copy */\
		/* Forward acyclic-convo weight is (wt_re, +wt_im): */\
		"subq		%%rdi,%%rcx			\n\t"	\
		"mulpd		%%xmm0,%%xmm4		\n\t"	/* [x     ]*wt_re */\
		"mulpd		%%xmm1,%%xmm3		\n\t"	/* [y copy]*wt_im */\
		"movaps	 	%%xmm6,-0x20(%%rcx)	\n\t"	/* Store maxerr */\
		"addq		%%r10,%%rcx			\n\t"	\
		"mulpd		%%xmm0,%%xmm2		\n\t"	/* [y     ]*wt_re */\
		"mulpd		%%xmm1,%%xmm5		\n\t"	/* [x copy]*wt_im */\
		"movaps		(%%rcx),%%xmm0		\n\t"	/* [wt0,wt1] */\
		"subpd		%%xmm3,%%xmm4		\n\t"	/* rt = x*wt_re - y*wt_im */\
		"addpd		%%xmm2,%%xmm5		\n\t"	/* it = x*wt_im + y*wt_re */\
		/* Forward IBDWT weight: */\
		"mulpd		%%xmm0,%%xmm4		\n\t"	\
		"mulpd		%%xmm0,%%xmm5		\n\t"	\
		"movaps		%%xmm4,    (%%rdx)	\n\t"	/* store rt = ~[a.re,b.re] */\
		"movaps		%%xmm5,0x10(%%rdx)	\n\t"	/* store it = ~[a.im,b.im] */\
		/* Prepare for next pair of complex data: */\
		"addq	%[__idx_incr],%%rsi		\n\t"	/* idx_offset += idx_incr */\
		"mov	%%esi, %[__idx_offset]	\n\t"	/* Store incremented idx_offset */\
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
		/* Prefetch address */\
		,	[__add0] "m" (Xadd0)\
		: "cc","memory","rax","rbx","rcx","rdx","rsi","rdi","r10","r11","r14","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"   /* Clobbered registers */\
	);\
	}

	/* Same non-power-of-2-transform deal as above, but use xmm8-15 to process 2 sets of carries side-by-side.
	Data/Carry #2 assumed offset by +0x20/0x10 from #1 (which are accessed via the [__data/__cy] pointers, resp.)
	i/jcycle0 and i/jcycle1 are the address offsets needed for IBDWT array indexing for the 2 resp. carries. */
	#define SSE2_fermat_carry_norm_errcheck_X2(Xdata,Xcy,Xnrt_bits,Xnrtm1,Xidx_offset,Xidx_incr,Xodd_radix,Xhalf_arr,Xsign_mask,Xadd1,Xadd2,Xicycle0,Xjcycle0,Xicycle1,Xjcycle1, Xadd0,Xp1)\
	{\
	__asm__ volatile (\
	"movq	%[__add0],%%r14	\n\t"/* base address for 2 prefetches-from-main-data-array spread through this macro */\
	"prefetcht0	(%%r14)	\n\t"\
		/* lcol -> rcol index analogs: [rsi,rax,rbx] -> [r10,r11,r12], [rcx,rdx,rdi] shared */\
		"movslq	%[__idx_offset],%%rax	\n\t		movslq	%[__idx_incr],%%r10		\n\t"\
		"movslq		%[__nrt_bits],%%rcx	\n\t		addq	%%rax,%%r10				\n\t"\
		"movslq		%[__nrtm1],%%rdi	\n\t		movq		%%r10,%%r11	\n\t"/* r10 contains idx_offset2, i.e. is the rcol-analog of idx_offset1 in lcol: */\
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
		"movslq	%[__idx_offset],%%rax	\n\t		movq		%%r10,%%r11			\n\t"\
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
		"movslq	%[__idx_offset],%%rax	\n\t		movq		%%r10,%%r11			\n\t"\
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
		"movslq	%[__odd_radix],%%rdi	\n\t"\
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
	"movslq		%[__p1],%%r15	\n\t"	\
	"prefetcht0	(%%r14,%%r15,8)	\n\t"\
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
		"movslq	%[__jcycle0],%%r8 		\n\t		movslq	%[__jcycle1],%%r10		\n\t"\
		"addq		%%r8 ,%%rcx			\n\t		addq		%%r10,%%r12			\n\t"\
		"maxpd		%%xmm14,%%xmm6		/* Save larger of maxerr1,2: */	\n\t"\
		"mulpd	(%%rcx,%%rdi),%%xmm2	\n\t		mulpd	(%%r12,%%rdi),%%xmm10	\n\t"\
		"roundpd	$0,%%xmm2,%%xmm2	\n\t		roundpd	$0,%%xmm10,%%xmm10		\n\t"\
		"shlq		$1,%%rdi														\n\t"\
		"movaps		%%xmm2,%%xmm3		\n\t		movaps		%%xmm10,%%xmm11		\n\t"\
		"mulpd		     (%%rcx),%%xmm3	\n\t		mulpd		     (%%r12),%%xmm11	\n\t"\
		"subq		%%r8 ,%%rcx			\n\t		subq		%%r10,%%r12			\n\t"\
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
		"movslq	%[__icycle0],%%r8 		\n\t		movslq	%[__icycle1],%%r10		\n\t"\
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
		,	[__odd_radix]	"m" (Xodd_radix)\
		,	[__half_arr]	"m" (Xhalf_arr)\
		,	[__sign_mask]	"m" (Xsign_mask)\
		,	[__add1]		"m" (Xadd1)\
		,	[__add2]		"m" (Xadd2)\
		,	[__icycle0]		"m" (Xicycle0)\
		,	[__jcycle0]		"m" (Xjcycle0)\
		,	[__icycle1]		"m" (Xicycle1)\
		,	[__jcycle1]		"m" (Xjcycle1)\
		/* Prefetch: base address and 1 index offset */\
		,	[__add0] "m" (Xadd0)\
		,	[__p1] "m" (Xp1)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","r8","r10","r11","r12","r14","r15","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14"		/* Clobbered registers */\
	);\
	}

#endif	// AVX or SSE2?

	/*************************************************************/
	/**************** MERSENNE-MOD CARRY MACROS ******************/
	/*************************************************************/

#ifdef USE_AVX

	/********* Names ending in _X4 denote "Genuine AVX" carry macros, which process 8 AVX-sized vector (= 32 doubles): **********/
	// Note that in this "true AVX" carry, [n_minus_sil,n_minus_silp1,sinwt,sinwtm1] are all pointers to (struct uint32x4) data, rather than 32bit ints as for SSE2-style macros!
	/*
	Aug 2013: The following "true AVX" carry macros - fusing the fancy-indexing footwork of the legacy SSE2 mersenne-mod-DWT
	carry macros and the AVX data-permute aspects of the AVX-based Fermat-mod carry macros - ran incredibly, awfully, unbelievably
	slowly in the initial implementation. I eventually traced that back to the mixing of legacy SSE instructions (using xmm-form registers)
	in the indexing-computation portions of the code with AVX instructions used for weights and carries in the new AVX code. The solution
	was to simply prepend a "v" to the legacy SSE instructions and (for ones where the VEX form of the instruction adds a third operand)
	to duplicate the original SRC+DEST operand (rightmost on this AT&T/GCC-syntax inline ASM) in order to satisfy the 3-operand syntax.

	CF. Intel's own "Mixing SSE and AVX bad, very bad" cautions in the following "Avoiding AVX-SSE Transition Penalties" PDF:

		http://software.intel.com/sites/default/files/m/d/4/1/d/8/11MC12_Avoiding_2BAVX-SSE_2BTransition_2BPenalties_2Brh_2Bfinal.pdf

	Here is the money snippet:

		"When using Intel® AVX instructions, it is important to know that mixing 256-bit Intel® AVX instructions
		 with legacy (non VEX-encoded) Intel® SSE instructions may result in penalties that could impact performance.
		 256-bit Intel® AVX instructions operate on the 256-bit YMM registers which are 256-bit extensions of the
		 existing 128-bit XMM registers. 128-bit Intel® AVX instructions operate on the lower 128 bits of the YMM
		 registers and zero the upper 128 bits. However, legacy Intel® SSE instructions operate on the XMM registers
		 and have no knowledge of the upper 128 bits of the YMM registers. Because of this, the hardware saves the
		 contents of the upper 128 bits of the YMM registers when transitioning from 256-bit Intel® AVX to legacy
		 Intel® SSE, and then restores these values when transitioning back from Intel® SSE to Intel® AVX (256-bit
		 or 128-bit). The save and restore operations both cause a penalty that amounts to several tens of clock
		 cycles for each operation."

	Cf. also Agner Fog's "early in the AVX life cycle" commentary at http://software.intel.com/en-us/forums/topic/301853

	To convey a sense of just how severe the timing penalties resulting from such SSE/AVX instruction mixing can be,
	here are sample timings for my AVX-enabled Mlucas code running at an FFT length of 4096 kdoubles on my 3.4 GHz quad-core
	Haswell system, using 2 threads, with a full-time 4-threaded Mlucas run [ongoing multimonth run of F28 @15360 K] as background load:

	[1] The baseline timing here is set by the well-tested and quite fast pure-AVX Fermat-mod carry macros:
		time ./Mlucas -f26 -fftlen 4096 -iters 100 -radset 0 -nthread 2
		...
		100 iterations of F26 with FFT length 4194304 = 4096 K
		Res64: A42BECD80DAEC4CB. AvgMaxErr = 0.005018834. MaxErr = 0.005859375. Program: E3.0x
		real	0m5.298s
	*	user	0m5.172s

	[2] Now do Mersene-mod run @same FFT length, using AVX-based FFT-pass code but SSE2-based carry macros,
		tweaked to take account of the differing AVX data layout:
		gcc -c -O3 -DUSE_THREADS -DUSE_AVX radix32*cy*c && gcc -o Mlucas *.o -lm -lpthread
		time ./Mlucas -fftlen 4096 -iters 100 -radset 0 -nthread 2
		...
		100 iterations of M77597293 with FFT length 4194304 = 4096 K
		Res64: 8CC30E314BF3E556. AvgMaxErr = 0.293526786. MaxErr = 0.343750000. Program: E3.0x
		real	0m6.673s
	*	user	0m6.760s

	...which is ~35% slower than the Fermat-mod run at the same FFT length - not a complete and total disaster but still
	very bad compared to the expected 10-20% runtime hit here for the Mersenne-mod computation.
	We hope for better using true-AVX mode for the carry step. Alas, our initial tests are "beyond unpromising":

	[3] gcc -c -O3 -DUSE_THREADS -DUSE_AVX -DUSE_AVX_CARRIES radix32*cy*c && gcc -o Mlucas *.o -lm -lpthread
		time ./Mlucas -fftlen 4096 -iters 100 -radset 0 -nthread 2
		...
		real	0m13.521s
	*	user	0m12.857s	<*** Same FFT code + AVX-based carries more than doubles the overall runtime!!! ***

	[4] Somehow I had failed to come across any of the literature discussing the SSE/AVX mixed-code performance
	penalty in my AVX-related reading. Once I realized that the key difference in instruction mix between the Fermat-mod
	and Mersenne-mod "true AVX" carry code was the presence of non-AVX SSE instructions in the latter, I quickly found
	the above docs detailing the performance hit such code entails, fixed the offending macros and immediately saw rather
	more promising timings, to say the least:

		gcc -c -O3 -DUSE_THREADS -DUSE_AVX -DUSE_AVX_CARRIES radix32*cy*c && gcc -o Mlucas *.o -lm -lpthread
		time ./Mlucas -fftlen 4096 -iters 100 -radset 0 -nthread 2
		...
		real	0m6.287s
	*	user	0m5.952s	<*** More than 2x faster than [3]!

	I.e. the performance *penalty* component alone from using the mixed SSE/AVX carry macros here was greater
	than the entire *runtime* needed for either the "pure-AVX-FFT/pure-SSE-carry" hybrid in [2] or the "true AVX 4 all"
	code, which only suffers a ~15% runtime penalty versus the analogous Fermat-mod code.
	*/
	#define AVX_cmplx_carry_norm_pow2_errcheck0_X4(Xdata,XwtA,XwtB,XwtC,Xcy,Xbjmod_0,Xhalf_arr,Xi,Xn_minus_silp1,Xn_minus_sil,Xsign_mask,Xsinwt,Xsinwtm1,Xsse_bw,Xsse_nm1,Xsse_sw, Xadd0,Xp1,Xp2,Xp3)\
	{\
	__asm__ volatile (\
	"movq	%[__add0],%%r14	\n\t"/* base address for 2 prefetches-from-main-data-array spread through this macro */\
		"movq		%[__data],%%rax		\n\t"\
	/* 4-way transpose of inputs (Re, Im parts separately): Inputs from r0/1,2/3,4/5.6/7. Outputs into ymm0-7: */\
	/* Real parts use ymm0,2,4,6, ymm8 as tmp-reg:					Imag parts use ymm1,3,5,7, ymm9 as tm-reg: */\
		"vmovaps	     (%%rax),%%ymm4						\n\t		vmovaps	0x020(%%rax),%%ymm5							\n\t"\
		"vmovaps	0x040(%%rax),%%ymm2						\n\t		vmovaps	0x060(%%rax),%%ymm3							\n\t"\
		"vshufpd	$15,%%ymm2,%%ymm4,%%ymm6				\n\t		vshufpd	$15,%%ymm3,%%ymm5,%%ymm7					\n\t"\
		"vshufpd	$0 ,%%ymm2,%%ymm4,%%ymm4				\n\t		vshufpd	$0 ,%%ymm3,%%ymm5,%%ymm5					\n\t"\
		"vmovaps	0x080(%%rax),%%ymm8						\n\t		vmovaps	0x0a0(%%rax),%%ymm9							\n\t"\
		"vmovaps	0x0c0(%%rax),%%ymm2						\n\t		vmovaps	0x0e0(%%rax),%%ymm3							\n\t"\
		"vshufpd	$15,%%ymm2,%%ymm8,%%ymm0				\n\t		vshufpd	$15,%%ymm3,%%ymm9,%%ymm1					\n\t"\
		"vshufpd	$0 ,%%ymm2,%%ymm8,%%ymm8				\n\t		vshufpd	$0 ,%%ymm3,%%ymm9,%%ymm9					\n\t"\
		"vperm2f128 $32,%%ymm0,%%ymm6,%%ymm2	/* Re B	*/	\n\t		vperm2f128 $32,%%ymm1,%%ymm7,%%ymm3		/* Im B	*/	\n\t"\
		"vperm2f128 $49,%%ymm0,%%ymm6,%%ymm6	/* Re D	*/	\n\t		vperm2f128 $49,%%ymm1,%%ymm7,%%ymm7		/* Im D	*/	\n\t"\
		"vperm2f128 $32,%%ymm8,%%ymm4,%%ymm0	/* Re A	*/	\n\t		vperm2f128 $32,%%ymm9,%%ymm5,%%ymm1 	/* Im A	*/	\n\t"\
		"vperm2f128 $49,%%ymm8,%%ymm4,%%ymm4	/* Re C	*/	\n\t		vperm2f128 $49,%%ymm9,%%ymm5,%%ymm5 	/* Im C	*/	\n\t"\
		/* Won't need main-array again until output transpose, so re-use rax for half_arr */\
		"movq	 %[__half_arr],%%rax	\n\t"\
		/* half_arr + 16*[0,1,2,3] = [wt,wt_inv,base,baseinv] */\
		"movq		%[__cy],%%rbx	\n\t"\
		"vmovaps	(%%rbx),%%ymm14	\n\t"/* ymm14 = cy_in */\
		"vmovaps	-0x40(%%rax),%%ymm15	\n\t"/* ymm15 = maxerr */\
	/**********************************/\
	/* Do A.re-quartet: Data in ymm0: */\
	/**********************************/\
	"prefetcht0	(%%r14)	\n\t"\
		"movq	%[__bjmod_0],%%rsi		\n\t"\
		"vmovaps	(%%rsi),	%%xmm8 		\n\t"/* bjmod[0:3]. PERSISTENT COPY OF BJMOD[0:3] REMAINS IN XMM8. */\
		"movq	%[__sse_sw],%%rsi		\n\t"\
		"vmovaps	(%%rsi),	%%xmm11		\n\t"/* sw[0:3] */\
		"vpsubd		%%xmm8 ,%%xmm11,%%xmm11	\n\t"/* sw[0:3] - bjmod[0:3] */\
		"vmovmskps	%%xmm11,	%%rsi	\n\t"/* Extract sign bits into 4-bit signmask <i3|i2|i1|i0>; idxs into base/inv table */\
		"movslq	%[__i]	,%%rcx			\n\t"/* i0=i for first block: */\
	"andq	$0xfffffffffffffffe,%%rsi	\n\t"/* Mask off lowest bit */\
		"addq	%%rcx	,%%rsi			\n\t"\
		"\n\t"\
		"movq	%[__n_minus_sil],%%rcx	\n\t"\
		"vmovd	0x0(%%rcx)	,%%xmm9 	\n\t"/* n_minus_sil in low 32 bits of xmm9  */\
		"vpshufd	$0,	%%xmm9 ,%%xmm9 		\n\t"/* Broadcast low 32 bits of xmm9  to all 4 slots of xmm9  */\
		"vpsubd		%%xmm8 ,%%xmm9 ,%%xmm9 		\n\t"/* n_minus_sil - bjmod[0:3] */\
		"vmovmskps	%%xmm9 ,%%rcx		\n\t"/* Extract sign bits into 4-bit signmask <m3|m2|m1|m0>; idxs into base/inv tables -> byte[2] of ecx... */\
		"\n\t"\
		"movq	%[__sinwt]	,%%rdx		\n\t"\
		"vmovd	0x0(%%rdx)	,%%xmm10		\n\t"/* sinwt in low 32 bits of xmm10 */\
		"vpshufd	$0,	%%xmm10,%%xmm10		\n\t"/* Broadcast low 32 bits of xmm10 to all 4 slots of xmm10 */\
		"vmovaps		%%xmm8 ,%%xmm11		\n\t"/* xmm11 = bjmod[0:3] copy */\
		"vpsubd		%%xmm10,%%xmm11,%%xmm11		\n\t"/* bjmod[0:3] - sinwt */\
		"vmovmskps	%%xmm11,%%rdx		\n\t"/* Extract sign bits into 4-bit signmask <n3|n2|n1|n0>; idxs into base/inv tables -> byte[1] of edx... */\
		"\n\t"\
	"shlq	$5,%%rsi	\n\t"/* i0123 << 5 (= lshift to give ptr offset for ymm-size data */\
	"shlq	$5,%%rcx	\n\t"/* m0123 << 5 (= lshift to give ptr offset for ymm-size data */\
	"shlq	$5,%%rdx	\n\t"/* n0123 << 5 (= lshift to give ptr offset for ymm-size data */\
		"\n\t"\
		"movq	%[__half_arr],%%rdi	\n\t"\
		"movq	%[__wtA]	,%%rax		\n\t"\
		"movq	%[__wtB]	,%%rbx		\n\t"\
		"vmovaps	     (%%rax),%%ymm12 	\n\t"/* wtA[j  ]; ebx FREE */\
		"vmovups	-0x10(%%rbx),%%ymm13	\n\t"/* wtB[j-1]; load doubles from rcx+[-0x10,-0x08, 0, +0x08] - note this address is ymm (16-byte) but not ymm (32-byte) aligned */\
		/* reverse-running indexing used for inv-wts really means we need to reverse ordering of 4 doubles d0-3 in ymm13 */\
		"vshufpd	$5,%%ymm13,%%ymm13,%%ymm13	\n\t"/* d[0123] -> d[1032] */\
		"vperm2f128 $1,%%ymm13,%%ymm13,%%ymm13	\n\t"/* d[1032] -> d[3210] */\
		"\n\t"\
		"addq	%%rdi,%%rcx		\n\t"\
		"addq	%%rdi,%%rdx		\n\t"\
		/* SSE2 version has double-copies in wtl/wtn ... AVX replaces redundant-data loads with load-with-broadcast: */\
		"vbroadcastsd 0x800(%%rdi),%%ymm9 	\n\t"/* wtl */\
		"vbroadcastsd 0x808(%%rdi),%%ymm10	\n\t"/* wtn */\
		"vmulpd	%%ymm12,%%ymm9 ,%%ymm9 		\n\t"/* wt   =wtA*wtl */\
		"vmulpd	%%ymm13,%%ymm10,%%ymm10		\n\t"/* wtinv=wtB*wtn */\
		"vmulpd	     (%%rcx),%%ymm9 ,%%ymm9 	\n\t"/* wt    *= one_half[m0123] */\
		"vmulpd	0x200(%%rdx),%%ymm10,%%ymm10	\n\t"/* wtinv *= one_half[16+n0123] */\
		"\n\t"\
		"vmulpd		%%ymm10,%%ymm0,%%ymm0		\n\t"/* x *= wtinv; ymm10 FREE */\
		"vaddpd		%%ymm14,%%ymm0,%%ymm0		\n\t"/* x *= wtinv + cy */\
		"vmovaps	%%ymm0,%%ymm10		\n\t"/* temp = x */\
		"vroundpd	$0,%%ymm10,%%ymm10	\n\t"/* temp = DNINT(x) */\
		"\n\t"\
		"movq	%[__sign_mask],%%rbx	\n\t"\
		"vsubpd		%%ymm10,%%ymm0,%%ymm0	\n\t"/* x - temp */\
		"vandpd		(%%rbx),%%ymm0,%%ymm0	\n\t"/* frac = fabs(x-temp) */\
		"vmaxpd		%%ymm15,%%ymm0,%%ymm15	\n\t"/* if(frac > maxerr) maxerr=frac */\
		/* cy   = DNINT(temp*baseinv[i1]): */\
		"addq	%%rdi,%%rsi		\n\t"\
		"vmovaps	%%ymm10	,%%ymm0				\n\t"/* cpy temp */\
		"vmulpd	0x600(%%rsi),%%ymm10,%%ymm10	\n\t"/* temp*baseinv[i0123] */\
		"vroundpd	$0,%%ymm10,%%ymm14			\n\t"/* cy_out */\
		"vmovaps	%%ymm14,%%ymm10				\n\t"/* cy = cpy cy_out */\
		/* x = (temp-cy*base[i1])*wt: */\
		"vmulpd	0x400(%%rsi),%%ymm10,%%ymm10	\n\t"/* cy*base[i0123] */\
		"vsubpd		%%ymm10,%%ymm0 ,%%ymm0 		\n\t"/* x = (temp-cy*base[i0123]) */\
		"vmulpd		%%ymm9 ,%%ymm0 ,%%ymm0 		\n\t"/* x *= wt */\
		/* Get ready for next set [IM0~] : */\
		"movq	%[__sse_bw]	,%%rax		\n\t"\
		"movq	%[__sse_nm1],%%rbx		\n\t"\
		"vpaddd		(%%rax),%%xmm8 ,%%xmm8 	\n\t"/* bjmod[0:3] += bw ; must use packed-INTEGER add [not addpd!] here, severe performance penalty from using addpd. */\
		"vpand		(%%rbx),%%xmm8 ,%%xmm8 	\n\t"/* bjmod[0:3] &= nm1; & doesn't care whether integer [pand] or floating [andpd], but data are int, so use pand for form's sake */\
		"\n\t"\
	/**********************************/\
	/* Do A.im-quartet: Data in ymm1: */\
	/**********************************/\
		"movq	%[__sse_sw],%%rsi		\n\t"\
		"vmovaps	(%%rsi),	%%xmm11		\n\t"/* sw[0:3] */\
		"vpsubd		%%xmm8 ,%%xmm11,%%xmm11	\n\t"/* sw[0:3] - bjmod[0:3] */\
		"vmovmskps	%%xmm11,	%%rsi	\n\t"/* Extract sign bits into 4-bit signmask <i3|i2|i1|i0>; idxs into base/inv table */\
		"\n\t"\
		"movq	%[__n_minus_silp1],%%rcx	\n\t"\
		"vmovd	0x0(%%rcx)	,%%xmm9 	\n\t"/* n_minus_sil in low 32 bits of xmm9  */\
		"vpshufd	$0,	%%xmm9 ,%%xmm9 		\n\t"/* Broadcast low 32 bits of xmm9  to all 4 slots of xmm9  */\
		"vpsubd		%%xmm8 ,%%xmm9 ,%%xmm9 		\n\t"/* n_minus_sil - bjmod[0:3] */\
		"vmovmskps	%%xmm9 ,%%rcx		\n\t"/* Extract sign bits into 4-bit signmask <m3|m2|m1|m0>; idxs into base/inv tables -> byte[2] of ecx... */\
		"\n\t"\
		"movq	%[__sinwtm1]	,%%rdx		\n\t"\
		"vmovd	0x0(%%rdx)	,%%xmm10		\n\t"/* sinwt in low 32 bits of xmm10 */\
		"vpshufd	$0,	%%xmm10,%%xmm10		\n\t"/* Broadcast low 32 bits of xmm10 to all 4 slots of xmm10 */\
		"vmovaps		%%xmm8 ,%%xmm11		\n\t"/* xmm11 = bjmod[0:3] copy */\
		"vpsubd		%%xmm10,%%xmm11,%%xmm11		\n\t"/* bjmod[0:3] - sinwt */\
		"vmovmskps	%%xmm11,%%rdx		\n\t"/* Extract sign bits into 4-bit signmask <n3|n2|n1|n0>; idxs into base/inv tables -> byte[1] of edx... */\
		"\n\t"\
	"shlq	$5,%%rsi	\n\t"/* i0123 << 5 (= lshift to give ptr offset for ymm-size data */\
	"shlq	$5,%%rcx	\n\t"/* m0123 << 5 (= lshift to give ptr offset for ymm-size data */\
	"shlq	$5,%%rdx	\n\t"/* n0123 << 5 (= lshift to give ptr offset for ymm-size data */\
		"\n\t"\
		"movq	%[__wtC]	,%%rbx		\n\t"/* wtA unchanged; wtB == wtC for remaining 7 of 8 sets of carries */\
		"vmovups	-0x10(%%rbx),%%ymm13	\n\t"/* wtC[j-1]; load doubles from rcx+[-0x10,-0x08, 0, +0x08] - note this address is ymm (16-byte) but not ymm (32-byte) aligned */\
		"vshufpd	$5,%%ymm13,%%ymm13,%%ymm13	\n\t"/* d[0123] -> d[1032] */\
		"vperm2f128 $1,%%ymm13,%%ymm13,%%ymm13	\n\t"/* d[1032] -> d[3210] */\
		"\n\t"\
		"addq	%%rdi,%%rcx		\n\t"\
		"addq	%%rdi,%%rdx		\n\t"\
		/* SSE2 version has double-copies in wtl/wtn ... AVX replaces redundant-data loads with load-with-broadcast: */\
		"vbroadcastsd 0x810(%%rdi),%%ymm9 	\n\t"/* wtlp1 */\
		"vbroadcastsd 0x818(%%rdi),%%ymm10	\n\t"/* wtnm1 */\
		"vmulpd	%%ymm12,%%ymm9 ,%%ymm9 		\n\t"/* wt   =wtA*wtlp1 */\
		"vmulpd	%%ymm13,%%ymm10,%%ymm10		\n\t"/* wtinv=wtC*wtnm1 */\
		"vmulpd	     (%%rcx),%%ymm9 ,%%ymm9 	\n\t"/* wt    *= one_half[m0123] */\
		"vmulpd	0x200(%%rdx),%%ymm10,%%ymm10	\n\t"/* wtinv *= one_half[16+n0123] */\
		"\n\t"\
		"vmulpd		%%ymm10,%%ymm1,%%ymm1		\n\t"/* x *= wtinv; ymm10 FREE */\
		"vaddpd		%%ymm14,%%ymm1,%%ymm1		\n\t"/* x *= wtinv + cy */\
		"vmovaps	%%ymm1,%%ymm10		\n\t"/* temp = x */\
		"vroundpd	$0,%%ymm10,%%ymm10	\n\t"/* temp = DNINT(x) */\
		"\n\t"\
		"movq	%[__sign_mask],%%rbx	\n\t"\
		"vsubpd		%%ymm10,%%ymm1,%%ymm1	\n\t"/* x - temp */\
		"vandpd		(%%rbx),%%ymm1,%%ymm1	\n\t"/* frac = fabs(x-temp) */\
		"vmaxpd		%%ymm15,%%ymm1,%%ymm15	\n\t"/* if(frac > maxerr) maxerr=frac */\
		/* cy   = DNINT(temp*baseinv[i1]): */\
		"addq	%%rdi,%%rsi		\n\t"\
		"vmovaps	%%ymm10	,%%ymm1				\n\t"/* cpy temp */\
		"vmulpd	0x600(%%rsi),%%ymm10,%%ymm10	\n\t"/* temp*baseinv[i0123] */\
		"vroundpd	$0,%%ymm10,%%ymm14			\n\t"/* cy_out */\
		"vmovaps	%%ymm14,%%ymm10				\n\t"/* cy = cpy cy_out */\
		/* x = (temp-cy*base[i1])*wt: */\
		"vmulpd	0x400(%%rsi),%%ymm10,%%ymm10	\n\t"/* cy*base[i0123] */\
		"vsubpd		%%ymm10,%%ymm1 ,%%ymm1 		\n\t"/* x = (temp-cy*base[i0123]) */\
		"vmulpd		%%ymm9 ,%%ymm1 ,%%ymm1 		\n\t"/* x *= wt */\
		/* Get ready for next set [RE1~] : */\
		"movq	%[__sse_bw]	,%%rax		\n\t"\
		"movq	%[__sse_nm1],%%rbx		\n\t"\
		"vpaddd		(%%rax),%%xmm8 ,%%xmm8 	\n\t"\
		"vpand		(%%rbx),%%xmm8 ,%%xmm8 	\n\t"\
		"\n\t"\
	/**********************************/\
	/* Do B.re-quartet: Data in ymm2: */\
	/**********************************/\
	"movslq		%[__p1],%%r15	\n\t"	\
	"prefetcht0	(%%r14,%%r15,8)	\n\t"\
		"movq	%[__sse_sw],%%rsi		\n\t"\
		"vmovaps	(%%rsi),	%%xmm11		\n\t"\
		"vpsubd		%%xmm8 ,%%xmm11,%%xmm11	\n\t"\
		"vmovmskps	%%xmm11,	%%rsi	\n\t"\
		"\n\t"\
		"movq	%[__n_minus_sil],%%rcx	\n\t"\
		"vmovd	0x4(%%rcx)	,%%xmm9 	\n\t"/* .d1 term of index quartet */\
		"vpshufd	$0,	%%xmm9 ,%%xmm9 		\n\t"\
		"vpsubd		%%xmm8 ,%%xmm9 ,%%xmm9 		\n\t"\
		"vmovmskps	%%xmm9 ,%%rcx		\n\t"\
		"\n\t"\
		"movq	%[__sinwt]	,%%rdx		\n\t"\
		"vmovd	0x4(%%rdx)	,%%xmm10	\n\t"/* .d1 term of index quartet */\
		"vpshufd	$0,	%%xmm10,%%xmm10		\n\t"\
		"vmovaps		%%xmm8 ,%%xmm11		\n\t"\
		"vpsubd		%%xmm10,%%xmm11,%%xmm11		\n\t"\
		"vmovmskps	%%xmm11,%%rdx		\n\t"\
		"\n\t"\
	"shlq	$5,%%rsi	\n\t"/* i0123 */\
	"shlq	$5,%%rcx	\n\t"/* m0123 */\
	"shlq	$5,%%rdx	\n\t"/* n0123 */\
		"\n\t"\
		"addq	%%rdi,%%rcx		\n\t"\
		"addq	%%rdi,%%rdx		\n\t"\
		"vbroadcastsd 0x820(%%rdi),%%ymm9 	\n\t"/* wtl */\
		"vbroadcastsd 0x828(%%rdi),%%ymm10	\n\t"/* wtn */\
		"vmulpd	%%ymm12,%%ymm9 ,%%ymm9 		\n\t"\
		"vmulpd	%%ymm13,%%ymm10,%%ymm10		\n\t"\
		"vmulpd	     (%%rcx),%%ymm9 ,%%ymm9 \n\t"\
		"vmulpd	0x200(%%rdx),%%ymm10,%%ymm10\n\t"\
		"\n\t"\
		"vmulpd		%%ymm10,%%ymm2,%%ymm2	\n\t"\
		"vaddpd		%%ymm14,%%ymm2,%%ymm2	\n\t"\
		"vmovaps	%%ymm2,%%ymm10			\n\t"\
		"vroundpd	$0,%%ymm10,%%ymm10		\n\t"\
		"\n\t"\
		"movq	%[__sign_mask],%%rbx		\n\t"\
		"vsubpd		%%ymm10,%%ymm2,%%ymm2	\n\t"\
		"vandpd		(%%rbx),%%ymm2,%%ymm2	\n\t"\
		"vmaxpd		%%ymm15,%%ymm2,%%ymm15	\n\t"\
		/* cy   = DNINT(temp*baseinv[i1]): */\
		"addq	%%rdi,%%rsi		\n\t"\
		"vmovaps	%%ymm10	,%%ymm2			\n\t"\
		"vmulpd	0x600(%%rsi),%%ymm10,%%ymm10\n\t"\
		"vroundpd	$0,%%ymm10,%%ymm14		\n\t"\
		"vmovaps	%%ymm14,%%ymm10				\n\t"/* cy = cpy cy_out */\
		/* x = (temp-cy*base[i1])*wt: */\
		"vmulpd	0x400(%%rsi),%%ymm10,%%ymm10	\n\t"/* cy*base[i0123] */\
		"vsubpd		%%ymm10,%%ymm2 ,%%ymm2 	\n\t"\
		"vmulpd		%%ymm9 ,%%ymm2 ,%%ymm2 	\n\t"\
		/* Get ready for next set [IM1~] : */\
		"movq	%[__sse_bw]	,%%rax		\n\t"\
		"movq	%[__sse_nm1],%%rbx		\n\t"\
		"vpaddd		(%%rax),%%xmm8 ,%%xmm8 	\n\t"\
		"vpand		(%%rbx),%%xmm8 ,%%xmm8 	\n\t"\
		"\n\t"\
	/**********************************/\
	/* Do B.im-quartet: Data in ymm3: */\
	/**********************************/\
		"movq	%[__sse_sw],%%rsi		\n\t"\
		"vmovaps	(%%rsi),	%%xmm11		\n\t"\
		"vpsubd		%%xmm8 ,%%xmm11,%%xmm11	\n\t"\
		"vmovmskps	%%xmm11,	%%rsi	\n\t"\
		"\n\t"\
		"movq	%[__n_minus_silp1],%%rcx	\n\t"\
		"vmovd	0x4(%%rcx)	,%%xmm9 	\n\t"/* .d1 term of index quartet */\
		"vpshufd	$0,	%%xmm9 ,%%xmm9 		\n\t"\
		"vpsubd		%%xmm8 ,%%xmm9 ,%%xmm9 		\n\t"\
		"vmovmskps	%%xmm9 ,%%rcx		\n\t"\
		"\n\t"\
		"movq	%[__sinwtm1]	,%%rdx		\n\t"\
		"vmovd	0x4(%%rdx)	,%%xmm10	\n\t"/* .d1 term of index quartet */\
		"vpshufd	$0,	%%xmm10,%%xmm10		\n\t"\
		"vmovaps		%%xmm8 ,%%xmm11		\n\t"\
		"vpsubd		%%xmm10,%%xmm11,%%xmm11		\n\t"\
		"vmovmskps	%%xmm11,%%rdx		\n\t"\
		"\n\t"\
	"shlq	$5,%%rsi	\n\t"/* i0123 */\
	"shlq	$5,%%rcx	\n\t"/* m0123 */\
	"shlq	$5,%%rdx	\n\t"/* n0123 */\
		"\n\t"\
		"addq	%%rdi,%%rcx		\n\t"\
		"addq	%%rdi,%%rdx		\n\t"\
		"vbroadcastsd 0x830(%%rdi),%%ymm9 	\n\t"/* wtlp1 */\
		"vbroadcastsd 0x838(%%rdi),%%ymm10	\n\t"/* wtnm1 */\
		"vmulpd	%%ymm12,%%ymm9 ,%%ymm9 		\n\t"\
		"vmulpd	%%ymm13,%%ymm10,%%ymm10		\n\t"\
		"vmulpd	     (%%rcx),%%ymm9 ,%%ymm9 \n\t"\
		"vmulpd	0x200(%%rdx),%%ymm10,%%ymm10\n\t"\
		"\n\t"\
		"vmulpd		%%ymm10,%%ymm3,%%ymm3	\n\t"\
		"vaddpd		%%ymm14,%%ymm3,%%ymm3	\n\t"\
		"vmovaps	%%ymm3,%%ymm10			\n\t"\
		"vroundpd	$0,%%ymm10,%%ymm10		\n\t"\
		"\n\t"\
		"movq	%[__sign_mask],%%rbx		\n\t"\
		"vsubpd		%%ymm10,%%ymm3,%%ymm3	\n\t"\
		"vandpd		(%%rbx),%%ymm3,%%ymm3	\n\t"\
		"vmaxpd		%%ymm15,%%ymm3,%%ymm15	\n\t"\
		/* cy   = DNINT(temp*baseinv[i1]): */\
		"addq	%%rdi,%%rsi		\n\t"\
		"vmovaps	%%ymm10	,%%ymm3			\n\t"\
		"vmulpd	0x600(%%rsi),%%ymm10,%%ymm10\n\t"\
		"vroundpd	$0,%%ymm10,%%ymm14		\n\t"\
		"vmovaps	%%ymm14,%%ymm10				\n\t"/* cy = cpy cy_out */\
		/* x = (temp-cy*base[i1])*wt: */\
		"vmulpd	0x400(%%rsi),%%ymm10,%%ymm10	\n\t"/* cy*base[i0123] */\
		"vsubpd		%%ymm10,%%ymm3 ,%%ymm3 	\n\t"\
		"vmulpd		%%ymm9 ,%%ymm3 ,%%ymm3 	\n\t"\
		/* Get ready for next set [RE2~] : */\
		"movq	%[__sse_bw]	,%%rax		\n\t"\
		"movq	%[__sse_nm1],%%rbx		\n\t"\
		"vpaddd		(%%rax),%%xmm8 ,%%xmm8 	\n\t"\
		"vpand		(%%rbx),%%xmm8 ,%%xmm8 	\n\t"\
		"\n\t"\
	/**********************************/\
	/* Do C.re-quartet: Data in ymm4: */\
	/**********************************/\
	"movslq		%[__p2],%%r15	\n\t"	\
	"prefetcht0	(%%r14,%%r15,8)	\n\t"\
		"movq	%[__sse_sw],%%rsi		\n\t"\
		"vmovaps	(%%rsi),	%%xmm11		\n\t"\
		"vpsubd		%%xmm8 ,%%xmm11,%%xmm11	\n\t"\
		"vmovmskps	%%xmm11,	%%rsi	\n\t"\
		"\n\t"\
		"movq	%[__n_minus_sil],%%rcx	\n\t"\
		"vmovd	0x8(%%rcx)	,%%xmm9 	\n\t"/* .d2 term of index quartet */\
		"vpshufd	$0,	%%xmm9 ,%%xmm9 		\n\t"\
		"vpsubd		%%xmm8 ,%%xmm9 ,%%xmm9 		\n\t"\
		"vmovmskps	%%xmm9 ,%%rcx		\n\t"\
		"\n\t"\
		"movq	%[__sinwt]	,%%rdx		\n\t"\
		"vmovd	0x8(%%rdx)	,%%xmm10	\n\t"/* .d2 term of index quartet */\
		"vpshufd	$0,	%%xmm10,%%xmm10		\n\t"\
		"vmovaps		%%xmm8 ,%%xmm11		\n\t"\
		"vpsubd		%%xmm10,%%xmm11,%%xmm11		\n\t"\
		"vmovmskps	%%xmm11,%%rdx		\n\t"\
		"\n\t"\
	"shlq	$5,%%rsi	\n\t"/* i0123 */\
	"shlq	$5,%%rcx	\n\t"/* m0123 */\
	"shlq	$5,%%rdx	\n\t"/* n0123 */\
		"\n\t"\
		"addq	%%rdi,%%rcx		\n\t"\
		"addq	%%rdi,%%rdx		\n\t"\
		"vbroadcastsd 0x840(%%rdi),%%ymm9 	\n\t"/* wtl */\
		"vbroadcastsd 0x848(%%rdi),%%ymm10	\n\t"/* wtn */\
		"vmulpd	%%ymm12,%%ymm9 ,%%ymm9 		\n\t"\
		"vmulpd	%%ymm13,%%ymm10,%%ymm10		\n\t"\
		"vmulpd	     (%%rcx),%%ymm9 ,%%ymm9 \n\t"\
		"vmulpd	0x200(%%rdx),%%ymm10,%%ymm10\n\t"\
		"\n\t"\
		"vmulpd		%%ymm10,%%ymm4,%%ymm4	\n\t"\
		"vaddpd		%%ymm14,%%ymm4,%%ymm4	\n\t"\
		"vmovaps	%%ymm4,%%ymm10			\n\t"\
		"vroundpd	$0,%%ymm10,%%ymm10		\n\t"\
		"\n\t"\
		"movq	%[__sign_mask],%%rbx		\n\t"\
		"vsubpd		%%ymm10,%%ymm4,%%ymm4	\n\t"\
		"vandpd		(%%rbx),%%ymm4,%%ymm4	\n\t"\
		"vmaxpd		%%ymm15,%%ymm4,%%ymm15	\n\t"\
		/* cy   = DNINT(temp*baseinv[i1]): */\
		"addq	%%rdi,%%rsi		\n\t"\
		"vmovaps	%%ymm10	,%%ymm4			\n\t"\
		"vmulpd	0x600(%%rsi),%%ymm10,%%ymm10\n\t"\
		"vroundpd	$0,%%ymm10,%%ymm14		\n\t"\
		"vmovaps	%%ymm14,%%ymm10				\n\t"/* cy = cpy cy_out */\
		/* x = (temp-cy*base[i1])*wt: */\
		"vmulpd	0x400(%%rsi),%%ymm10,%%ymm10	\n\t"/* cy*base[i0123] */\
		"vsubpd		%%ymm10,%%ymm4 ,%%ymm4 	\n\t"\
		"vmulpd		%%ymm9 ,%%ymm4 ,%%ymm4 	\n\t"\
		/* Get ready for next set [IM2~] : */\
		"movq	%[__sse_bw]	,%%rax		\n\t"\
		"movq	%[__sse_nm1],%%rbx		\n\t"\
		"vpaddd		(%%rax),%%xmm8 ,%%xmm8 	\n\t"\
		"vpand		(%%rbx),%%xmm8 ,%%xmm8 	\n\t"\
		"\n\t"\
	/**********************************/\
	/* Do C.im-quartet: Data in ymm5: */\
	/**********************************/\
		"movq	%[__sse_sw],%%rsi		\n\t"\
		"vmovaps	(%%rsi),	%%xmm11		\n\t"\
		"vpsubd		%%xmm8 ,%%xmm11,%%xmm11	\n\t"\
		"vmovmskps	%%xmm11,	%%rsi	\n\t"\
		"\n\t"\
		"movq	%[__n_minus_silp1],%%rcx	\n\t"\
		"vmovd	0x8(%%rcx)	,%%xmm9 	\n\t"/* .d2 term of index quartet */\
		"vpshufd	$0,	%%xmm9 ,%%xmm9 		\n\t"\
		"vpsubd		%%xmm8 ,%%xmm9 ,%%xmm9 		\n\t"\
		"vmovmskps	%%xmm9 ,%%rcx		\n\t"\
		"\n\t"\
		"movq	%[__sinwtm1]	,%%rdx		\n\t"\
		"vmovd	0x8(%%rdx)	,%%xmm10	\n\t"/* .d2 term of index quartet */\
		"vpshufd	$0,	%%xmm10,%%xmm10		\n\t"\
		"vmovaps		%%xmm8 ,%%xmm11		\n\t"\
		"vpsubd		%%xmm10,%%xmm11,%%xmm11		\n\t"\
		"vmovmskps	%%xmm11,%%rdx		\n\t"\
		"\n\t"\
	"shlq	$5,%%rsi	\n\t"/* i0123 */\
	"shlq	$5,%%rcx	\n\t"/* m0123 */\
	"shlq	$5,%%rdx	\n\t"/* n0123 */\
		"\n\t"\
		"addq	%%rdi,%%rcx		\n\t"\
		"addq	%%rdi,%%rdx		\n\t"\
		"vbroadcastsd 0x850(%%rdi),%%ymm9 	\n\t"/* wtlp1 */\
		"vbroadcastsd 0x858(%%rdi),%%ymm10	\n\t"/* wtnm1 */\
		"vmulpd	%%ymm12,%%ymm9 ,%%ymm9 		\n\t"\
		"vmulpd	%%ymm13,%%ymm10,%%ymm10		\n\t"\
		"vmulpd	     (%%rcx),%%ymm9 ,%%ymm9 \n\t"\
		"vmulpd	0x200(%%rdx),%%ymm10,%%ymm10\n\t"\
		"\n\t"\
		"vmulpd		%%ymm10,%%ymm5,%%ymm5	\n\t"\
		"vaddpd		%%ymm14,%%ymm5,%%ymm5	\n\t"\
		"vmovaps	%%ymm5,%%ymm10			\n\t"\
		"vroundpd	$0,%%ymm10,%%ymm10		\n\t"\
		"\n\t"\
		"movq	%[__sign_mask],%%rbx		\n\t"\
		"vsubpd		%%ymm10,%%ymm5,%%ymm5	\n\t"\
		"vandpd		(%%rbx),%%ymm5,%%ymm5	\n\t"\
		"vmaxpd		%%ymm15,%%ymm5,%%ymm15	\n\t"\
		/* cy   = DNINT(temp*baseinv[i1]): */\
		"addq	%%rdi,%%rsi		\n\t"\
		"vmovaps	%%ymm10	,%%ymm5			\n\t"\
		"vmulpd	0x600(%%rsi),%%ymm10,%%ymm10\n\t"\
		"vroundpd	$0,%%ymm10,%%ymm14		\n\t"\
		"vmovaps	%%ymm14,%%ymm10				\n\t"/* cy = cpy cy_out */\
		/* x = (temp-cy*base[i1])*wt: */\
		"vmulpd	0x400(%%rsi),%%ymm10,%%ymm10	\n\t"/* cy*base[i0123] */\
		"vsubpd		%%ymm10,%%ymm5 ,%%ymm5 	\n\t"\
		"vmulpd		%%ymm9 ,%%ymm5 ,%%ymm5 	\n\t"\
		/* Get ready for next set [RE3~] : */\
		"movq	%[__sse_bw]	,%%rax		\n\t"\
		"movq	%[__sse_nm1],%%rbx		\n\t"\
		"vpaddd		(%%rax),%%xmm8 ,%%xmm8 	\n\t"\
		"vpand		(%%rbx),%%xmm8 ,%%xmm8 	\n\t"\
		"\n\t"\
	/**********************************/\
	/* Do D.re-quartet: Data in ymm6: */\
	/**********************************/\
	"movslq		%[__p3],%%r15	\n\t"	\
	"prefetcht0	(%%r14,%%r15,8)	\n\t"\
		"movq	%[__sse_sw],%%rsi		\n\t"\
		"vmovaps	(%%rsi),	%%xmm11		\n\t"\
		"vpsubd		%%xmm8 ,%%xmm11,%%xmm11	\n\t"\
		"vmovmskps	%%xmm11,	%%rsi	\n\t"\
		"\n\t"\
		"movq	%[__n_minus_sil],%%rcx	\n\t"\
		"vmovd	0xC(%%rcx)	,%%xmm9 	\n\t"/* .d3 term of index quartet */\
		"vpshufd	$0,	%%xmm9 ,%%xmm9 		\n\t"\
		"vpsubd		%%xmm8 ,%%xmm9 ,%%xmm9 		\n\t"\
		"vmovmskps	%%xmm9 ,%%rcx		\n\t"\
		"\n\t"\
		"movq	%[__sinwt]	,%%rdx		\n\t"\
		"vmovd	0xC(%%rdx)	,%%xmm10	\n\t"/* .d3 term of index quartet */\
		"vpshufd	$0,	%%xmm10,%%xmm10		\n\t"\
		"vmovaps		%%xmm8 ,%%xmm11		\n\t"\
		"vpsubd		%%xmm10,%%xmm11,%%xmm11		\n\t"\
		"vmovmskps	%%xmm11,%%rdx		\n\t"\
		"\n\t"\
	"shlq	$5,%%rsi	\n\t"/* i0123 */\
	"shlq	$5,%%rcx	\n\t"/* m0123 */\
	"shlq	$5,%%rdx	\n\t"/* n0123 */\
		"\n\t"\
		"addq	%%rdi,%%rcx		\n\t"\
		"addq	%%rdi,%%rdx		\n\t"\
		"vbroadcastsd 0x860(%%rdi),%%ymm9 	\n\t"/* wtl */\
		"vbroadcastsd 0x868(%%rdi),%%ymm10	\n\t"/* wtn */\
		"vmulpd	%%ymm12,%%ymm9 ,%%ymm9 		\n\t"\
		"vmulpd	%%ymm13,%%ymm10,%%ymm10		\n\t"\
		"vmulpd	     (%%rcx),%%ymm9 ,%%ymm9 \n\t"\
		"vmulpd	0x200(%%rdx),%%ymm10,%%ymm10\n\t"\
		"\n\t"\
		"vmulpd		%%ymm10,%%ymm6,%%ymm6	\n\t"\
		"vaddpd		%%ymm14,%%ymm6,%%ymm6	\n\t"\
		"vmovaps	%%ymm6,%%ymm10			\n\t"\
		"vroundpd	$0,%%ymm10,%%ymm10		\n\t"\
		"\n\t"\
		"movq	%[__sign_mask],%%rbx		\n\t"\
		"vsubpd		%%ymm10,%%ymm6,%%ymm6	\n\t"\
		"vandpd		(%%rbx),%%ymm6,%%ymm6	\n\t"\
		"vmaxpd		%%ymm15,%%ymm6,%%ymm15	\n\t"\
		/* cy   = DNINT(temp*baseinv[i1]): */\
		"addq	%%rdi,%%rsi		\n\t"\
		"vmovaps	%%ymm10	,%%ymm6			\n\t"\
		"vmulpd	0x600(%%rsi),%%ymm10,%%ymm10\n\t"\
		"vroundpd	$0,%%ymm10,%%ymm14		\n\t"\
		"vmovaps	%%ymm14,%%ymm10				\n\t"/* cy = cpy cy_out */\
		/* x = (temp-cy*base[i1])*wt: */\
		"vmulpd	0x400(%%rsi),%%ymm10,%%ymm10	\n\t"/* cy*base[i0123] */\
		"vsubpd		%%ymm10,%%ymm6 ,%%ymm6 	\n\t"\
		"vmulpd		%%ymm9 ,%%ymm6 ,%%ymm6 	\n\t"\
		/* Get ready for next set [IM3~] : */\
		"movq	%[__sse_bw]	,%%rax		\n\t"\
		"movq	%[__sse_nm1],%%rbx		\n\t"\
		"vpaddd		(%%rax),%%xmm8 ,%%xmm8 	\n\t"\
		"vpand		(%%rbx),%%xmm8 ,%%xmm8 	\n\t"\
		"\n\t"\
	/**********************************/\
	/* Do D.im-quartet: Data in ymm7: */\
	/**********************************/\
		"movq	%[__sse_sw],%%rsi		\n\t"\
		"vmovaps	(%%rsi),	%%xmm11		\n\t"\
		"vpsubd		%%xmm8 ,%%xmm11,%%xmm11	\n\t"\
		"vmovmskps	%%xmm11,	%%rsi	\n\t"\
		"\n\t"\
		"movq	%[__n_minus_silp1],%%rcx	\n\t"\
		"vmovd	0xC(%%rcx)	,%%xmm9 	\n\t"/* .d3 term of index quartet */\
		"vpshufd	$0,	%%xmm9 ,%%xmm9 		\n\t"\
		"vpsubd		%%xmm8 ,%%xmm9 ,%%xmm9 		\n\t"\
		"vmovmskps	%%xmm9 ,%%rcx		\n\t"\
		"\n\t"\
		"movq	%[__sinwtm1]	,%%rdx		\n\t"\
		"vmovd	0xC(%%rdx)	,%%xmm10	\n\t"/* .d3 term of index quartet */\
		"vpshufd	$0,	%%xmm10,%%xmm10		\n\t"\
		"vmovaps		%%xmm8 ,%%xmm11		\n\t"\
		"vpsubd		%%xmm10,%%xmm11,%%xmm11		\n\t"\
		"vmovmskps	%%xmm11,%%rdx		\n\t"\
		"\n\t"\
	"shlq	$5,%%rsi	\n\t"/* i0123 */\
	"shlq	$5,%%rcx	\n\t"/* m0123 */\
	"shlq	$5,%%rdx	\n\t"/* n0123 */\
		"\n\t"\
		"addq	%%rdi,%%rcx		\n\t"\
		"addq	%%rdi,%%rdx		\n\t"\
		"vbroadcastsd 0x870(%%rdi),%%ymm9 	\n\t"/* wtlp1 */\
		"vbroadcastsd 0x878(%%rdi),%%ymm10	\n\t"/* wtnm1 */\
		"vmulpd	%%ymm12,%%ymm9 ,%%ymm9 		\n\t"\
		"vmulpd	%%ymm13,%%ymm10,%%ymm10		\n\t"\
		"vmulpd	     (%%rcx),%%ymm9 ,%%ymm9 \n\t"\
		"vmulpd	0x200(%%rdx),%%ymm10,%%ymm10\n\t"\
		"\n\t"\
		"vmulpd		%%ymm10,%%ymm7,%%ymm7	\n\t"\
		"vaddpd		%%ymm14,%%ymm7,%%ymm7	\n\t"\
		"vmovaps	%%ymm7,%%ymm10			\n\t"\
		"vroundpd	$0,%%ymm10,%%ymm10		\n\t"\
		"\n\t"\
		"movq	%[__sign_mask],%%rbx		\n\t"\
		"vsubpd		%%ymm10,%%ymm7,%%ymm7	\n\t"\
		"vandpd		(%%rbx),%%ymm7,%%ymm7	\n\t"\
		"vmaxpd		%%ymm15,%%ymm7,%%ymm15	\n\t"\
		/* cy   = DNINT(temp*baseinv[i1]): */\
		"addq	%%rdi,%%rsi		\n\t"\
		"vmovaps	%%ymm10	,%%ymm7			\n\t"\
		"vmulpd	0x600(%%rsi),%%ymm10,%%ymm10\n\t"\
		"vroundpd	$0,%%ymm10,%%ymm14		\n\t"\
		"vmovaps	%%ymm14,%%ymm10				\n\t"/* cy = cpy cy_out */\
		/* x = (temp-cy*base[i1])*wt: */\
		"vmulpd	0x400(%%rsi),%%ymm10,%%ymm10	\n\t"/* cy*base[i0123] */\
		"vsubpd		%%ymm10,%%ymm7 ,%%ymm7 	\n\t"\
		"vmulpd		%%ymm9 ,%%ymm7 ,%%ymm7 	\n\t"\
		/* Update wts-array pointers in preparation for next call of the macro: */\
		"movq	%[__wtA]	,%%rax		\n\t"\
		"movq	%[__wtB]	,%%rbx		\n\t"\
		"movq	%[__wtC]	,%%rcx		\n\t"\
		"addq	$0x20	,%%rax			\n\t"/* add0 += 4 */\
		"subq	$0x20	,%%rbx			\n\t"/* add1 -= 4 */\
		"subq	$0x20	,%%rcx			\n\t"/* add2 -= 4 */\
		"movq	%%rax	,%[__wtA]		\n\t"\
		"movq	%%rbx	,%[__wtB]		\n\t"\
		"movq	%%rcx	,%[__wtC]		\n\t"\
		"\n\t"\
		/* Get ready for store of final-updated bjmod[0:3] values: */\
		"movq	%[__sse_bw]	,%%rax		\n\t"\
		"movq	%[__sse_nm1],%%rbx		\n\t"\
		"vpaddd		(%%rax),%%xmm8 ,%%xmm8 	\n\t"\
		"vpand		(%%rbx),%%xmm8 ,%%xmm8 	\n\t"\
		"\n\t"\
		/* Store bjmodn index quartet: */\
		"movq	%[__bjmod_0],%%rcx		\n\t"\
		"vmovaps	%%xmm8,(%%rcx)			\n\t"\
		/* Store cy_out: */\
		"movq		%[__cy],%%rbx	\n\t"\
		"vmovaps	%%ymm14,(%%rbx)	\n\t"/* ymm14 = cy_in */\
		/* Store maxerr: */\
		"movq		%[__half_arr],%%rdx		\n\t"\
		"vmovaps	%%ymm15,-0x40(%%rdx)	\n\t"\
	/* 4-way transpose of outputs (Re, Im parts separately): Inputs from ymm0-7. Outputs into r0/1,2/3,4/5.6/7: */	\
	/* Because default inputs for our 4 x 4 transpose macro (e.g. the one used at start of this carry macro) */\
	/* are into ymm4/2/8/2, munge inputs into that order, resolving name-conflicts via use of the now-available ymm8-15 for outputs: */\
		"movq		%[__data],%%rax			\n\t"\
		"vshufpd	$15,%%ymm2,%%ymm0,%%ymm10					\n\t		vshufpd	$15,%%ymm3,%%ymm1,%%ymm11						\n\t"\
		"vshufpd	$0 ,%%ymm2,%%ymm0,%%ymm0					\n\t		vshufpd	$0 ,%%ymm3,%%ymm1,%%ymm1						\n\t"\
		"vshufpd	$15,%%ymm6,%%ymm4,%%ymm12					\n\t		vshufpd	$15,%%ymm7,%%ymm5,%%ymm13						\n\t"\
		"vshufpd	$0 ,%%ymm6,%%ymm4,%%ymm4					\n\t		vshufpd	$0 ,%%ymm7,%%ymm5,%%ymm5						\n\t"\
		"vperm2f128 $32,%%ymm12,%%ymm10,%%ymm2 		/* Re B	*/	\n\t		vperm2f128 $32,%%ymm13,%%ymm11,%%ymm3		/* Im B	*/	\n\t"\
		"vperm2f128 $49,%%ymm12,%%ymm10,%%ymm10		/* Re D	*/	\n\t		vperm2f128 $49,%%ymm13,%%ymm11,%%ymm11		/* Im D	*/	\n\t"\
		"vperm2f128 $32,%%ymm4 ,%%ymm0 ,%%ymm12		/* Re A	*/	\n\t		vperm2f128 $32,%%ymm5 ,%%ymm1 ,%%ymm13 		/* Im A	*/	\n\t"\
		"vperm2f128 $49,%%ymm4 ,%%ymm0 ,%%ymm0 		/* Re C	*/	\n\t		vperm2f128 $49,%%ymm5 ,%%ymm1 ,%%ymm1		/* Im C	*/	\n\t"\
		/* And write 'em back to memory: */\
		"vmovaps	%%ymm12,     (%%rax)						\n\t		vmovaps	%%ymm13,0x020(%%rax)				\n\t"\
		"vmovaps	%%ymm2 ,0x040(%%rax)						\n\t		vmovaps	%%ymm3 ,0x060(%%rax)				\n\t"\
		"vmovaps	%%ymm0 ,0x080(%%rax)						\n\t		vmovaps	%%ymm1 ,0x0a0(%%rax)				\n\t"\
		"vmovaps	%%ymm10,0x0c0(%%rax)						\n\t		vmovaps	%%ymm11,0x0e0(%%rax)				\n\t"\
		:					/* outputs: none */\
		: [__data]		"m" (Xdata)	/* All inputs from memory addresses here */\
		, [__wtA]		"m" (XwtA)		\
		, [__wtB]		"m" (XwtB)		\
		, [__wtC]		"m" (XwtC)		\
		, [__cy]		"m" (Xcy)		\
		, [__bjmod_0]	"m" (Xbjmod_0)		\
		, [__half_arr]	"m" (Xhalf_arr)		\
		, [__i]			"m" (Xi)			\
		, [__n_minus_silp1] "m" (Xn_minus_silp1)\
		, [__n_minus_sil]	"m" (Xn_minus_sil)	\
		, [__sign_mask]	"m" (Xsign_mask)	\
		, [__sinwt]		"m" (Xsinwt)		\
		, [__sinwtm1]	"m" (Xsinwtm1)		\
		, [__sse_bw]	"m" (Xsse_bw)		\
		, [__sse_nm1]	"m" (Xsse_nm1)		\
		, [__sse_sw]	"m" (Xsse_sw)		\
		/* Prefetch: base address and 3 index offsets */\
		,	[__add0] "m" (Xadd0)\
		,	[__p1] "m" (Xp1)\
		,	[__p2] "m" (Xp2)\
		,	[__p3] "m" (Xp3)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r14","r15","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
	);\
	}

	#define AVX_cmplx_carry_norm_pow2_errcheck1_X4(Xdata,XwtA,XwtB,XwtC,Xcy,Xbjmod_0,Xhalf_arr,   Xn_minus_silp1,Xn_minus_sil,Xsign_mask,Xsinwt,Xsinwtm1,Xsse_bw,Xsse_nm1,Xsse_sw, Xadd0,Xp1,Xp2,Xp3)\
	{\
	__asm__ volatile (\
	"movq	%[__add0],%%r14	\n\t"/* base address for 4 prefetches-from-main-data-array spread through this macro */\
		"movq		%[__data],%%rax		\n\t"\
	/* 4-way transpose of inputs (Re, Im parts separately): Inputs from r0/1,2/3,4/5.6/7. Outputs into ymm0-7: */\
	/* Real parts use ymm0,2,4,6, ymm8 as tmp-reg:					Imag parts use ymm1,3,5,7, ymm9 as tm-reg: */\
		"vmovaps	     (%%rax),%%ymm4						\n\t		vmovaps	0x020(%%rax),%%ymm5							\n\t"\
		"vmovaps	0x040(%%rax),%%ymm2						\n\t		vmovaps	0x060(%%rax),%%ymm3							\n\t"\
		"vshufpd	$15,%%ymm2,%%ymm4,%%ymm6				\n\t		vshufpd	$15,%%ymm3,%%ymm5,%%ymm7					\n\t"\
		"vshufpd	$0 ,%%ymm2,%%ymm4,%%ymm4				\n\t		vshufpd	$0 ,%%ymm3,%%ymm5,%%ymm5					\n\t"\
		"vmovaps	0x080(%%rax),%%ymm8						\n\t		vmovaps	0x0a0(%%rax),%%ymm9							\n\t"\
		"vmovaps	0x0c0(%%rax),%%ymm2						\n\t		vmovaps	0x0e0(%%rax),%%ymm3							\n\t"\
		"vshufpd	$15,%%ymm2,%%ymm8,%%ymm0				\n\t		vshufpd	$15,%%ymm3,%%ymm9,%%ymm1					\n\t"\
		"vshufpd	$0 ,%%ymm2,%%ymm8,%%ymm8				\n\t		vshufpd	$0 ,%%ymm3,%%ymm9,%%ymm9					\n\t"\
		"vperm2f128 $32,%%ymm0,%%ymm6,%%ymm2	/* Re B	*/	\n\t		vperm2f128 $32,%%ymm1,%%ymm7,%%ymm3		/* Im B	*/	\n\t"\
		"vperm2f128 $49,%%ymm0,%%ymm6,%%ymm6	/* Re D	*/	\n\t		vperm2f128 $49,%%ymm1,%%ymm7,%%ymm7		/* Im D	*/	\n\t"\
		"vperm2f128 $32,%%ymm8,%%ymm4,%%ymm0	/* Re A	*/	\n\t		vperm2f128 $32,%%ymm9,%%ymm5,%%ymm1 	/* Im A	*/	\n\t"\
		"vperm2f128 $49,%%ymm8,%%ymm4,%%ymm4	/* Re C	*/	\n\t		vperm2f128 $49,%%ymm9,%%ymm5,%%ymm5 	/* Im C	*/	\n\t"\
		/* Won't need main-array again until output transpose, so re-use rax for half_arr */\
		"movq	 %[__half_arr],%%rax	\n\t"\
		/* half_arr + 16*[0,1,2,3] = [wt,wt_inv,base,baseinv] */\
		"movq		%[__cy],%%rbx	\n\t"\
		"vmovaps	(%%rbx),%%ymm14	\n\t"/* ymm14 = cy_in */\
		"vmovaps	-0x40(%%rax),%%ymm15	\n\t"/* ymm15 = maxerr */\
	/**********************************/\
	/* Do A.re-quartet: Data in ymm0: */\
	/**********************************/\
	"prefetcht0	(%%r14)		\n\t"\
		"movq	%[__bjmod_0],%%rsi		\n\t"\
		"vmovaps	(%%rsi),	%%xmm8 		\n\t"/* bjmod[0:3]. PERSISTENT COPY OF BJMOD[0:3] REMAINS IN XMM8. */\
		"movq	%[__sse_sw],%%rsi		\n\t"\
		"vmovaps	(%%rsi),	%%xmm11		\n\t"/* sw[0:3] */\
		"vpsubd		%%xmm8 ,%%xmm11,%%xmm11	\n\t"/* sw[0:3] - bjmod[0:3] */\
		"vmovmskps	%%xmm11,	%%rsi	\n\t"/* Extract sign bits into 4-bit signmask <i3|i2|i1|i0>; idxs into base/inv table */\
		"\n\t"\
		"movq	%[__n_minus_sil],%%rcx	\n\t"\
		"vmovd	0x0(%%rcx)	,%%xmm9 	\n\t"/* n_minus_sil in low 32 bits of xmm9  */\
		"vpshufd	$0,	%%xmm9 ,%%xmm9 		\n\t"/* Broadcast low 32 bits of xmm9  to all 4 slots of xmm9  */\
		"vpsubd		%%xmm8 ,%%xmm9 ,%%xmm9 		\n\t"/* n_minus_sil - bjmod[0:3] */\
		"vmovmskps	%%xmm9 ,%%rcx		\n\t"/* Extract sign bits into 4-bit signmask <m3|m2|m1|m0>; idxs into base/inv tables -> byte[2] of ecx... */\
		"\n\t"\
		"movq	%[__sinwt]	,%%rdx		\n\t"\
		"vmovd	0x0(%%rdx)	,%%xmm10		\n\t"/* sinwt in low 32 bits of xmm10 */\
		"vpshufd	$0,	%%xmm10,%%xmm10		\n\t"/* Broadcast low 32 bits of xmm10 to all 4 slots of xmm10 */\
		"vmovaps		%%xmm8 ,%%xmm11		\n\t"/* xmm11 = bjmod[0:3] copy */\
		"vpsubd		%%xmm10,%%xmm11,%%xmm11		\n\t"/* bjmod[0:3] - sinwt */\
		"vmovmskps	%%xmm11,%%rdx		\n\t"/* Extract sign bits into 4-bit signmask <n3|n2|n1|n0>; idxs into base/inv tables -> byte[1] of edx... */\
		"\n\t"\
	"shlq	$5,%%rsi	\n\t"/* i0123 << 5 (= lshift to give ptr offset for ymm-size data */\
	"shlq	$5,%%rcx	\n\t"/* m0123 << 5 (= lshift to give ptr offset for ymm-size data */\
	"shlq	$5,%%rdx	\n\t"/* n0123 << 5 (= lshift to give ptr offset for ymm-size data */\
		"\n\t"\
		"movq	%[__half_arr],%%rdi	\n\t"\
		"movq	%[__wtA]	,%%rax		\n\t"\
		"movq	%[__wtB]	,%%rbx		\n\t"\
		"vmovaps	     (%%rax),%%ymm12 	\n\t"/* wtA[j  ]; ebx FREE */\
		"vmovups	-0x10(%%rbx),%%ymm13	\n\t"/* wtB[j-1]; load doubles from rcx+[-0x10,-0x08, 0, +0x08] - note this address is ymm (16-byte) but not ymm (32-byte) aligned */\
		/* reverse-running indexing used for inv-wts really means we need to reverse ordering of 4 doubles d0-3 in ymm13 */\
		"vshufpd	$5,%%ymm13,%%ymm13,%%ymm13	\n\t"/* d[0123] -> d[1032] */\
		"vperm2f128 $1,%%ymm13,%%ymm13,%%ymm13	\n\t"/* d[1032] -> d[3210] */\
		"\n\t"\
		"addq	%%rdi,%%rcx		\n\t"\
		"addq	%%rdi,%%rdx		\n\t"\
		/* SSE2 version has double-copies in wtl/wtn ... AVX replaces redundant-data loads with load-with-broadcast: */\
		"vbroadcastsd 0x800(%%rdi),%%ymm9 	\n\t"/* wtl */\
		"vbroadcastsd 0x808(%%rdi),%%ymm10	\n\t"/* wtn */\
		"vmulpd	%%ymm12,%%ymm9 ,%%ymm9 		\n\t"/* wt   =wtA*wtl */\
		"vmulpd	%%ymm13,%%ymm10,%%ymm10		\n\t"/* wtinv=wtB*wtn */\
		"vmulpd	     (%%rcx),%%ymm9 ,%%ymm9 	\n\t"/* wt    *= one_half[m0123] */\
		"vmulpd	0x200(%%rdx),%%ymm10,%%ymm10	\n\t"/* wtinv *= one_half[16+n0123] */\
		"\n\t"\
		"vmulpd		%%ymm10,%%ymm0,%%ymm0		\n\t"/* x *= wtinv; ymm10 FREE */\
		"vaddpd		%%ymm14,%%ymm0,%%ymm0		\n\t"/* x *= wtinv + cy */\
		"vmovaps	%%ymm0,%%ymm10		\n\t"/* temp = x */\
		"vroundpd	$0,%%ymm10,%%ymm10	\n\t"/* temp = DNINT(x) */\
		"\n\t"\
		"movq	%[__sign_mask],%%rbx	\n\t"\
		"vsubpd		%%ymm10,%%ymm0,%%ymm0	\n\t"/* x - temp */\
		"vandpd		(%%rbx),%%ymm0,%%ymm0	\n\t"/* frac = fabs(x-temp) */\
		"vmaxpd		%%ymm15,%%ymm0,%%ymm15	\n\t"/* if(frac > maxerr) maxerr=frac */\
		/* cy   = DNINT(temp*baseinv[i1]): */\
		"addq	%%rdi,%%rsi		\n\t"\
		"vmovaps	%%ymm10	,%%ymm0				\n\t"/* cpy temp */\
		"vmulpd	0x600(%%rsi),%%ymm10,%%ymm10	\n\t"/* temp*baseinv[i0123] */\
		"vroundpd	$0,%%ymm10,%%ymm14			\n\t"/* cy_out */\
		"vmovaps	%%ymm14,%%ymm10				\n\t"/* cy = cpy cy_out */\
		/* x = (temp-cy*base[i1])*wt: */\
		"vmulpd	0x400(%%rsi),%%ymm10,%%ymm10	\n\t"/* cy*base[i0123] */\
		"vsubpd		%%ymm10,%%ymm0 ,%%ymm0 		\n\t"/* x = (temp-cy*base[i0123]) */\
		"vmulpd		%%ymm9 ,%%ymm0 ,%%ymm0 		\n\t"/* x *= wt */\
		/* Get ready for next set [IM0~] : */\
		"movq	%[__sse_bw]	,%%rax		\n\t"\
		"movq	%[__sse_nm1],%%rbx		\n\t"\
		"vpaddd		(%%rax),%%xmm8 ,%%xmm8 	\n\t"/* bjmod[0:3] += bw ; must use packed-INTEGER add [not addpd!] here, severe performance penalty from using addpd. */\
		"vpand		(%%rbx),%%xmm8 ,%%xmm8 	\n\t"/* bjmod[0:3] &= nm1; & doesn't care whether integer [pand] or floating [andpd], but data are int, so use pand for form's sake */\
		"\n\t"\
	/**********************************/\
	/* Do A.im-quartet: Data in ymm1: */\
	/**********************************/\
		"movq	%[__sse_sw],%%rsi		\n\t"\
		"vmovaps	(%%rsi),	%%xmm11		\n\t"/* sw[0:3] */\
		"vpsubd		%%xmm8 ,%%xmm11,%%xmm11	\n\t"/* sw[0:3] - bjmod[0:3] */\
		"vmovmskps	%%xmm11,	%%rsi	\n\t"/* Extract sign bits into 4-bit signmask <i3|i2|i1|i0>; idxs into base/inv table */\
		"\n\t"\
		"movq	%[__n_minus_silp1],%%rcx	\n\t"\
		"vmovd	0x0(%%rcx)	,%%xmm9 	\n\t"/* n_minus_sil in low 32 bits of xmm9  */\
		"vpshufd	$0,	%%xmm9 ,%%xmm9 		\n\t"/* Broadcast low 32 bits of xmm9  to all 4 slots of xmm9  */\
		"vpsubd		%%xmm8 ,%%xmm9 ,%%xmm9 		\n\t"/* n_minus_sil - bjmod[0:3] */\
		"vmovmskps	%%xmm9 ,%%rcx		\n\t"/* Extract sign bits into 4-bit signmask <m3|m2|m1|m0>; idxs into base/inv tables -> byte[2] of ecx... */\
		"\n\t"\
		"movq	%[__sinwtm1]	,%%rdx		\n\t"\
		"vmovd	0x0(%%rdx)	,%%xmm10		\n\t"/* sinwt in low 32 bits of xmm10 */\
		"vpshufd	$0,	%%xmm10,%%xmm10		\n\t"/* Broadcast low 32 bits of xmm10 to all 4 slots of xmm10 */\
		"vmovaps		%%xmm8 ,%%xmm11		\n\t"/* xmm11 = bjmod[0:3] copy */\
		"vpsubd		%%xmm10,%%xmm11,%%xmm11		\n\t"/* bjmod[0:3] - sinwt */\
		"vmovmskps	%%xmm11,%%rdx		\n\t"/* Extract sign bits into 4-bit signmask <n3|n2|n1|n0>; idxs into base/inv tables -> byte[1] of edx... */\
		"\n\t"\
	"shlq	$5,%%rsi	\n\t"/* i0123 << 5 (= lshift to give ptr offset for ymm-size data */\
	"shlq	$5,%%rcx	\n\t"/* m0123 << 5 (= lshift to give ptr offset for ymm-size data */\
	"shlq	$5,%%rdx	\n\t"/* n0123 << 5 (= lshift to give ptr offset for ymm-size data */\
		"\n\t"\
		"movq	%[__wtC]	,%%rbx		\n\t"/* wtA unchanged; wtB == wtC for remaining 7 of 8 sets of carries */\
		"vmovups	-0x10(%%rbx),%%ymm13	\n\t"/* wtC[j-1]; load doubles from rcx+[-0x10,-0x08, 0, +0x08] - note this address is ymm (16-byte) but not ymm (32-byte) aligned */\
		"vshufpd	$5,%%ymm13,%%ymm13,%%ymm13	\n\t"/* d[0123] -> d[1032] */\
		"vperm2f128 $1,%%ymm13,%%ymm13,%%ymm13	\n\t"/* d[1032] -> d[3210] */\
		"\n\t"\
		"addq	%%rdi,%%rcx		\n\t"\
		"addq	%%rdi,%%rdx		\n\t"\
		/* SSE2 version has double-copies in wtl/wtn ... AVX replaces redundant-data loads with load-with-broadcast: */\
		"vbroadcastsd 0x810(%%rdi),%%ymm9 	\n\t"/* wtlp1 */\
		"vbroadcastsd 0x818(%%rdi),%%ymm10	\n\t"/* wtnm1 */\
		"vmulpd	%%ymm12,%%ymm9 ,%%ymm9 		\n\t"/* wt   =wtA*wtlp1 */\
		"vmulpd	%%ymm13,%%ymm10,%%ymm10		\n\t"/* wtinv=wtC*wtnm1 */\
		"vmulpd	     (%%rcx),%%ymm9 ,%%ymm9 	\n\t"/* wt    *= one_half[m0123] */\
		"vmulpd	0x200(%%rdx),%%ymm10,%%ymm10	\n\t"/* wtinv *= one_half[16+n0123] */\
		"\n\t"\
		"vmulpd		%%ymm10,%%ymm1,%%ymm1		\n\t"/* x *= wtinv; ymm10 FREE */\
		"vaddpd		%%ymm14,%%ymm1,%%ymm1		\n\t"/* x *= wtinv + cy */\
		"vmovaps	%%ymm1,%%ymm10		\n\t"/* temp = x */\
		"vroundpd	$0,%%ymm10,%%ymm10	\n\t"/* temp = DNINT(x) */\
		"\n\t"\
		"movq	%[__sign_mask],%%rbx	\n\t"\
		"vsubpd		%%ymm10,%%ymm1,%%ymm1	\n\t"/* x - temp */\
		"vandpd		(%%rbx),%%ymm1,%%ymm1	\n\t"/* frac = fabs(x-temp) */\
		"vmaxpd		%%ymm15,%%ymm1,%%ymm15	\n\t"/* if(frac > maxerr) maxerr=frac */\
		/* cy   = DNINT(temp*baseinv[i1]): */\
		"addq	%%rdi,%%rsi		\n\t"\
		"vmovaps	%%ymm10	,%%ymm1				\n\t"/* cpy temp */\
		"vmulpd	0x600(%%rsi),%%ymm10,%%ymm10	\n\t"/* temp*baseinv[i0123] */\
		"vroundpd	$0,%%ymm10,%%ymm14			\n\t"/* cy_out */\
		"vmovaps	%%ymm14,%%ymm10				\n\t"/* cy = cpy cy_out */\
		/* x = (temp-cy*base[i1])*wt: */\
		"vmulpd	0x400(%%rsi),%%ymm10,%%ymm10	\n\t"/* cy*base[i0123] */\
		"vsubpd		%%ymm10,%%ymm1 ,%%ymm1 		\n\t"/* x = (temp-cy*base[i0123]) */\
		"vmulpd		%%ymm9 ,%%ymm1 ,%%ymm1 		\n\t"/* x *= wt */\
		/* Get ready for next set [RE1~] : */\
		"movq	%[__sse_bw]	,%%rax		\n\t"\
		"movq	%[__sse_nm1],%%rbx		\n\t"\
		"vpaddd		(%%rax),%%xmm8 ,%%xmm8 	\n\t"\
		"vpand		(%%rbx),%%xmm8 ,%%xmm8 	\n\t"\
		"\n\t"\
	/**********************************/\
	/* Do B.re-quartet: Data in ymm2: */\
	/**********************************/\
	"movslq		%[__p1],%%r15	\n\t"	\
	"prefetcht0	(%%r14,%%r15,8)	\n\t"\
		"movq	%[__sse_sw],%%rsi		\n\t"\
		"vmovaps	(%%rsi),	%%xmm11		\n\t"\
		"vpsubd		%%xmm8 ,%%xmm11,%%xmm11	\n\t"\
		"vmovmskps	%%xmm11,	%%rsi	\n\t"\
		"\n\t"\
		"movq	%[__n_minus_sil],%%rcx	\n\t"\
		"vmovd	0x4(%%rcx)	,%%xmm9 	\n\t"/* .d1 term of index quartet */\
		"vpshufd	$0,	%%xmm9 ,%%xmm9 		\n\t"\
		"vpsubd		%%xmm8 ,%%xmm9 ,%%xmm9 		\n\t"\
		"vmovmskps	%%xmm9 ,%%rcx		\n\t"\
		"\n\t"\
		"movq	%[__sinwt]	,%%rdx		\n\t"\
		"vmovd	0x4(%%rdx)	,%%xmm10	\n\t"/* .d1 term of index quartet */\
		"vpshufd	$0,	%%xmm10,%%xmm10		\n\t"\
		"vmovaps		%%xmm8 ,%%xmm11		\n\t"\
		"vpsubd		%%xmm10,%%xmm11,%%xmm11		\n\t"\
		"vmovmskps	%%xmm11,%%rdx		\n\t"\
		"\n\t"\
	"shlq	$5,%%rsi	\n\t"/* i0123 */\
	"shlq	$5,%%rcx	\n\t"/* m0123 */\
	"shlq	$5,%%rdx	\n\t"/* n0123 */\
		"\n\t"\
		"addq	%%rdi,%%rcx		\n\t"\
		"addq	%%rdi,%%rdx		\n\t"\
		"vbroadcastsd 0x820(%%rdi),%%ymm9 	\n\t"/* wtl */\
		"vbroadcastsd 0x828(%%rdi),%%ymm10	\n\t"/* wtn */\
		"vmulpd	%%ymm12,%%ymm9 ,%%ymm9 		\n\t"\
		"vmulpd	%%ymm13,%%ymm10,%%ymm10		\n\t"\
		"vmulpd	     (%%rcx),%%ymm9 ,%%ymm9 \n\t"\
		"vmulpd	0x200(%%rdx),%%ymm10,%%ymm10\n\t"\
		"\n\t"\
		"vmulpd		%%ymm10,%%ymm2,%%ymm2	\n\t"\
		"vaddpd		%%ymm14,%%ymm2,%%ymm2	\n\t"\
		"vmovaps	%%ymm2,%%ymm10			\n\t"\
		"vroundpd	$0,%%ymm10,%%ymm10		\n\t"\
		"\n\t"\
		"movq	%[__sign_mask],%%rbx		\n\t"\
		"vsubpd		%%ymm10,%%ymm2,%%ymm2	\n\t"\
		"vandpd		(%%rbx),%%ymm2,%%ymm2	\n\t"\
		"vmaxpd		%%ymm15,%%ymm2,%%ymm15	\n\t"\
		/* cy   = DNINT(temp*baseinv[i1]): */\
		"addq	%%rdi,%%rsi		\n\t"\
		"vmovaps	%%ymm10	,%%ymm2			\n\t"\
		"vmulpd	0x600(%%rsi),%%ymm10,%%ymm10\n\t"\
		"vroundpd	$0,%%ymm10,%%ymm14		\n\t"\
		"vmovaps	%%ymm14,%%ymm10				\n\t"/* cy = cpy cy_out */\
		/* x = (temp-cy*base[i1])*wt: */\
		"vmulpd	0x400(%%rsi),%%ymm10,%%ymm10	\n\t"/* cy*base[i0123] */\
		"vsubpd		%%ymm10,%%ymm2 ,%%ymm2 	\n\t"\
		"vmulpd		%%ymm9 ,%%ymm2 ,%%ymm2 	\n\t"\
		/* Get ready for next set [IM1~] : */\
		"movq	%[__sse_bw]	,%%rax		\n\t"\
		"movq	%[__sse_nm1],%%rbx		\n\t"\
		"vpaddd		(%%rax),%%xmm8 ,%%xmm8 	\n\t"\
		"vpand		(%%rbx),%%xmm8 ,%%xmm8 	\n\t"\
		"\n\t"\
	/**********************************/\
	/* Do B.im-quartet: Data in ymm3: */\
	/**********************************/\
		"movq	%[__sse_sw],%%rsi		\n\t"\
		"vmovaps	(%%rsi),	%%xmm11		\n\t"\
		"vpsubd		%%xmm8 ,%%xmm11,%%xmm11	\n\t"\
		"vmovmskps	%%xmm11,	%%rsi	\n\t"\
		"\n\t"\
		"movq	%[__n_minus_silp1],%%rcx	\n\t"\
		"vmovd	0x4(%%rcx)	,%%xmm9 	\n\t"/* .d1 term of index quartet */\
		"vpshufd	$0,	%%xmm9 ,%%xmm9 		\n\t"\
		"vpsubd		%%xmm8 ,%%xmm9 ,%%xmm9 		\n\t"\
		"vmovmskps	%%xmm9 ,%%rcx		\n\t"\
		"\n\t"\
		"movq	%[__sinwtm1]	,%%rdx		\n\t"\
		"vmovd	0x4(%%rdx)	,%%xmm10	\n\t"/* .d1 term of index quartet */\
		"vpshufd	$0,	%%xmm10,%%xmm10		\n\t"\
		"vmovaps		%%xmm8 ,%%xmm11		\n\t"\
		"vpsubd		%%xmm10,%%xmm11,%%xmm11		\n\t"\
		"vmovmskps	%%xmm11,%%rdx		\n\t"\
		"\n\t"\
	"shlq	$5,%%rsi	\n\t"/* i0123 */\
	"shlq	$5,%%rcx	\n\t"/* m0123 */\
	"shlq	$5,%%rdx	\n\t"/* n0123 */\
		"\n\t"\
		"addq	%%rdi,%%rcx		\n\t"\
		"addq	%%rdi,%%rdx		\n\t"\
		"vbroadcastsd 0x830(%%rdi),%%ymm9 	\n\t"/* wtlp1 */\
		"vbroadcastsd 0x838(%%rdi),%%ymm10	\n\t"/* wtnm1 */\
		"vmulpd	%%ymm12,%%ymm9 ,%%ymm9 		\n\t"\
		"vmulpd	%%ymm13,%%ymm10,%%ymm10		\n\t"\
		"vmulpd	     (%%rcx),%%ymm9 ,%%ymm9 \n\t"\
		"vmulpd	0x200(%%rdx),%%ymm10,%%ymm10\n\t"\
		"\n\t"\
		"vmulpd		%%ymm10,%%ymm3,%%ymm3	\n\t"\
		"vaddpd		%%ymm14,%%ymm3,%%ymm3	\n\t"\
		"vmovaps	%%ymm3,%%ymm10			\n\t"\
		"vroundpd	$0,%%ymm10,%%ymm10		\n\t"\
		"\n\t"\
		"movq	%[__sign_mask],%%rbx		\n\t"\
		"vsubpd		%%ymm10,%%ymm3,%%ymm3	\n\t"\
		"vandpd		(%%rbx),%%ymm3,%%ymm3	\n\t"\
		"vmaxpd		%%ymm15,%%ymm3,%%ymm15	\n\t"\
		/* cy   = DNINT(temp*baseinv[i1]): */\
		"addq	%%rdi,%%rsi		\n\t"\
		"vmovaps	%%ymm10	,%%ymm3			\n\t"\
		"vmulpd	0x600(%%rsi),%%ymm10,%%ymm10\n\t"\
		"vroundpd	$0,%%ymm10,%%ymm14		\n\t"\
		"vmovaps	%%ymm14,%%ymm10				\n\t"/* cy = cpy cy_out */\
		/* x = (temp-cy*base[i1])*wt: */\
		"vmulpd	0x400(%%rsi),%%ymm10,%%ymm10	\n\t"/* cy*base[i0123] */\
		"vsubpd		%%ymm10,%%ymm3 ,%%ymm3 	\n\t"\
		"vmulpd		%%ymm9 ,%%ymm3 ,%%ymm3 	\n\t"\
		/* Get ready for next set [RE2~] : */\
		"movq	%[__sse_bw]	,%%rax		\n\t"\
		"movq	%[__sse_nm1],%%rbx		\n\t"\
		"vpaddd		(%%rax),%%xmm8 ,%%xmm8 	\n\t"\
		"vpand		(%%rbx),%%xmm8 ,%%xmm8 	\n\t"\
		"\n\t"\
	/**********************************/\
	/* Do C.re-quartet: Data in ymm4: */\
	/**********************************/\
	"movslq		%[__p2],%%r15	\n\t"	\
	"prefetcht0	(%%r14,%%r15,8)	\n\t"\
		"movq	%[__sse_sw],%%rsi		\n\t"\
		"vmovaps	(%%rsi),	%%xmm11		\n\t"\
		"vpsubd		%%xmm8 ,%%xmm11,%%xmm11	\n\t"\
		"vmovmskps	%%xmm11,	%%rsi	\n\t"\
		"\n\t"\
		"movq	%[__n_minus_sil],%%rcx	\n\t"\
		"vmovd	0x8(%%rcx)	,%%xmm9 	\n\t"/* .d2 term of index quartet */\
		"vpshufd	$0,	%%xmm9 ,%%xmm9 		\n\t"\
		"vpsubd		%%xmm8 ,%%xmm9 ,%%xmm9 		\n\t"\
		"vmovmskps	%%xmm9 ,%%rcx		\n\t"\
		"\n\t"\
		"movq	%[__sinwt]	,%%rdx		\n\t"\
		"vmovd	0x8(%%rdx)	,%%xmm10	\n\t"/* .d2 term of index quartet */\
		"vpshufd	$0,	%%xmm10,%%xmm10		\n\t"\
		"vmovaps		%%xmm8 ,%%xmm11		\n\t"\
		"vpsubd		%%xmm10,%%xmm11,%%xmm11		\n\t"\
		"vmovmskps	%%xmm11,%%rdx		\n\t"\
		"\n\t"\
	"shlq	$5,%%rsi	\n\t"/* i0123 */\
	"shlq	$5,%%rcx	\n\t"/* m0123 */\
	"shlq	$5,%%rdx	\n\t"/* n0123 */\
		"\n\t"\
		"addq	%%rdi,%%rcx		\n\t"\
		"addq	%%rdi,%%rdx		\n\t"\
		"vbroadcastsd 0x840(%%rdi),%%ymm9 	\n\t"/* wtl */\
		"vbroadcastsd 0x848(%%rdi),%%ymm10	\n\t"/* wtn */\
		"vmulpd	%%ymm12,%%ymm9 ,%%ymm9 		\n\t"\
		"vmulpd	%%ymm13,%%ymm10,%%ymm10		\n\t"\
		"vmulpd	     (%%rcx),%%ymm9 ,%%ymm9 \n\t"\
		"vmulpd	0x200(%%rdx),%%ymm10,%%ymm10\n\t"\
		"\n\t"\
		"vmulpd		%%ymm10,%%ymm4,%%ymm4	\n\t"\
		"vaddpd		%%ymm14,%%ymm4,%%ymm4	\n\t"\
		"vmovaps	%%ymm4,%%ymm10			\n\t"\
		"vroundpd	$0,%%ymm10,%%ymm10		\n\t"\
		"\n\t"\
		"movq	%[__sign_mask],%%rbx		\n\t"\
		"vsubpd		%%ymm10,%%ymm4,%%ymm4	\n\t"\
		"vandpd		(%%rbx),%%ymm4,%%ymm4	\n\t"\
		"vmaxpd		%%ymm15,%%ymm4,%%ymm15	\n\t"\
		/* cy   = DNINT(temp*baseinv[i1]): */\
		"addq	%%rdi,%%rsi		\n\t"\
		"vmovaps	%%ymm10	,%%ymm4			\n\t"\
		"vmulpd	0x600(%%rsi),%%ymm10,%%ymm10\n\t"\
		"vroundpd	$0,%%ymm10,%%ymm14		\n\t"\
		"vmovaps	%%ymm14,%%ymm10				\n\t"/* cy = cpy cy_out */\
		/* x = (temp-cy*base[i1])*wt: */\
		"vmulpd	0x400(%%rsi),%%ymm10,%%ymm10	\n\t"/* cy*base[i0123] */\
		"vsubpd		%%ymm10,%%ymm4 ,%%ymm4 	\n\t"\
		"vmulpd		%%ymm9 ,%%ymm4 ,%%ymm4 	\n\t"\
		/* Get ready for next set [IM2~] : */\
		"movq	%[__sse_bw]	,%%rax		\n\t"\
		"movq	%[__sse_nm1],%%rbx		\n\t"\
		"vpaddd		(%%rax),%%xmm8 ,%%xmm8 	\n\t"\
		"vpand		(%%rbx),%%xmm8 ,%%xmm8 	\n\t"\
		"\n\t"\
	/**********************************/\
	/* Do C.im-quartet: Data in ymm5: */\
	/**********************************/\
		"movq	%[__sse_sw],%%rsi		\n\t"\
		"vmovaps	(%%rsi),	%%xmm11		\n\t"\
		"vpsubd		%%xmm8 ,%%xmm11,%%xmm11	\n\t"\
		"vmovmskps	%%xmm11,	%%rsi	\n\t"\
		"\n\t"\
		"movq	%[__n_minus_silp1],%%rcx	\n\t"\
		"vmovd	0x8(%%rcx)	,%%xmm9 	\n\t"/* .d2 term of index quartet */\
		"vpshufd	$0,	%%xmm9 ,%%xmm9 		\n\t"\
		"vpsubd		%%xmm8 ,%%xmm9 ,%%xmm9 		\n\t"\
		"vmovmskps	%%xmm9 ,%%rcx		\n\t"\
		"\n\t"\
		"movq	%[__sinwtm1]	,%%rdx		\n\t"\
		"vmovd	0x8(%%rdx)	,%%xmm10	\n\t"/* .d2 term of index quartet */\
		"vpshufd	$0,	%%xmm10,%%xmm10		\n\t"\
		"vmovaps		%%xmm8 ,%%xmm11		\n\t"\
		"vpsubd		%%xmm10,%%xmm11,%%xmm11		\n\t"\
		"vmovmskps	%%xmm11,%%rdx		\n\t"\
		"\n\t"\
	"shlq	$5,%%rsi	\n\t"/* i0123 */\
	"shlq	$5,%%rcx	\n\t"/* m0123 */\
	"shlq	$5,%%rdx	\n\t"/* n0123 */\
		"\n\t"\
		"addq	%%rdi,%%rcx		\n\t"\
		"addq	%%rdi,%%rdx		\n\t"\
		"vbroadcastsd 0x850(%%rdi),%%ymm9 	\n\t"/* wtlp1 */\
		"vbroadcastsd 0x858(%%rdi),%%ymm10	\n\t"/* wtnm1 */\
		"vmulpd	%%ymm12,%%ymm9 ,%%ymm9 		\n\t"\
		"vmulpd	%%ymm13,%%ymm10,%%ymm10		\n\t"\
		"vmulpd	     (%%rcx),%%ymm9 ,%%ymm9 \n\t"\
		"vmulpd	0x200(%%rdx),%%ymm10,%%ymm10\n\t"\
		"\n\t"\
		"vmulpd		%%ymm10,%%ymm5,%%ymm5	\n\t"\
		"vaddpd		%%ymm14,%%ymm5,%%ymm5	\n\t"\
		"vmovaps	%%ymm5,%%ymm10			\n\t"\
		"vroundpd	$0,%%ymm10,%%ymm10		\n\t"\
		"\n\t"\
		"movq	%[__sign_mask],%%rbx		\n\t"\
		"vsubpd		%%ymm10,%%ymm5,%%ymm5	\n\t"\
		"vandpd		(%%rbx),%%ymm5,%%ymm5	\n\t"\
		"vmaxpd		%%ymm15,%%ymm5,%%ymm15	\n\t"\
		/* cy   = DNINT(temp*baseinv[i1]): */\
		"addq	%%rdi,%%rsi		\n\t"\
		"vmovaps	%%ymm10	,%%ymm5			\n\t"\
		"vmulpd	0x600(%%rsi),%%ymm10,%%ymm10\n\t"\
		"vroundpd	$0,%%ymm10,%%ymm14		\n\t"\
		"vmovaps	%%ymm14,%%ymm10				\n\t"/* cy = cpy cy_out */\
		/* x = (temp-cy*base[i1])*wt: */\
		"vmulpd	0x400(%%rsi),%%ymm10,%%ymm10	\n\t"/* cy*base[i0123] */\
		"vsubpd		%%ymm10,%%ymm5 ,%%ymm5 	\n\t"\
		"vmulpd		%%ymm9 ,%%ymm5 ,%%ymm5 	\n\t"\
		/* Get ready for next set [RE3~] : */\
		"movq	%[__sse_bw]	,%%rax		\n\t"\
		"movq	%[__sse_nm1],%%rbx		\n\t"\
		"vpaddd		(%%rax),%%xmm8 ,%%xmm8 	\n\t"\
		"vpand		(%%rbx),%%xmm8 ,%%xmm8 	\n\t"\
		"\n\t"\
	/**********************************/\
	/* Do D.re-quartet: Data in ymm6: */\
	/**********************************/\
	"movslq		%[__p3],%%r15	\n\t"	\
	"prefetcht0	(%%r14,%%r15,8)	\n\t"\
		"movq	%[__sse_sw],%%rsi		\n\t"\
		"vmovaps	(%%rsi),	%%xmm11		\n\t"\
		"vpsubd		%%xmm8 ,%%xmm11,%%xmm11	\n\t"\
		"vmovmskps	%%xmm11,	%%rsi	\n\t"\
		"\n\t"\
		"movq	%[__n_minus_sil],%%rcx	\n\t"\
		"vmovd	0xC(%%rcx)	,%%xmm9 	\n\t"/* .d3 term of index quartet */\
		"vpshufd	$0,	%%xmm9 ,%%xmm9 		\n\t"\
		"vpsubd		%%xmm8 ,%%xmm9 ,%%xmm9 		\n\t"\
		"vmovmskps	%%xmm9 ,%%rcx		\n\t"\
		"\n\t"\
		"movq	%[__sinwt]	,%%rdx		\n\t"\
		"vmovd	0xC(%%rdx)	,%%xmm10	\n\t"/* .d3 term of index quartet */\
		"vpshufd	$0,	%%xmm10,%%xmm10		\n\t"\
		"vmovaps		%%xmm8 ,%%xmm11		\n\t"\
		"vpsubd		%%xmm10,%%xmm11,%%xmm11		\n\t"\
		"vmovmskps	%%xmm11,%%rdx		\n\t"\
		"\n\t"\
	"shlq	$5,%%rsi	\n\t"/* i0123 */\
	"shlq	$5,%%rcx	\n\t"/* m0123 */\
	"shlq	$5,%%rdx	\n\t"/* n0123 */\
		"\n\t"\
		"addq	%%rdi,%%rcx		\n\t"\
		"addq	%%rdi,%%rdx		\n\t"\
		"vbroadcastsd 0x860(%%rdi),%%ymm9 	\n\t"/* wtl */\
		"vbroadcastsd 0x868(%%rdi),%%ymm10	\n\t"/* wtn */\
		"vmulpd	%%ymm12,%%ymm9 ,%%ymm9 		\n\t"\
		"vmulpd	%%ymm13,%%ymm10,%%ymm10		\n\t"\
		"vmulpd	     (%%rcx),%%ymm9 ,%%ymm9 \n\t"\
		"vmulpd	0x200(%%rdx),%%ymm10,%%ymm10\n\t"\
		"\n\t"\
		"vmulpd		%%ymm10,%%ymm6,%%ymm6	\n\t"\
		"vaddpd		%%ymm14,%%ymm6,%%ymm6	\n\t"\
		"vmovaps	%%ymm6,%%ymm10			\n\t"\
		"vroundpd	$0,%%ymm10,%%ymm10		\n\t"\
		"\n\t"\
		"movq	%[__sign_mask],%%rbx		\n\t"\
		"vsubpd		%%ymm10,%%ymm6,%%ymm6	\n\t"\
		"vandpd		(%%rbx),%%ymm6,%%ymm6	\n\t"\
		"vmaxpd		%%ymm15,%%ymm6,%%ymm15	\n\t"\
		/* cy   = DNINT(temp*baseinv[i1]): */\
		"addq	%%rdi,%%rsi		\n\t"\
		"vmovaps	%%ymm10	,%%ymm6			\n\t"\
		"vmulpd	0x600(%%rsi),%%ymm10,%%ymm10\n\t"\
		"vroundpd	$0,%%ymm10,%%ymm14		\n\t"\
		"vmovaps	%%ymm14,%%ymm10				\n\t"/* cy = cpy cy_out */\
		/* x = (temp-cy*base[i1])*wt: */\
		"vmulpd	0x400(%%rsi),%%ymm10,%%ymm10	\n\t"/* cy*base[i0123] */\
		"vsubpd		%%ymm10,%%ymm6 ,%%ymm6 	\n\t"\
		"vmulpd		%%ymm9 ,%%ymm6 ,%%ymm6 	\n\t"\
		/* Get ready for next set [IM3~] : */\
		"movq	%[__sse_bw]	,%%rax		\n\t"\
		"movq	%[__sse_nm1],%%rbx		\n\t"\
		"vpaddd		(%%rax),%%xmm8 ,%%xmm8 	\n\t"\
		"vpand		(%%rbx),%%xmm8 ,%%xmm8 	\n\t"\
		"\n\t"\
	/**********************************/\
	/* Do D.im-quartet: Data in ymm7: */\
	/**********************************/\
		"movq	%[__sse_sw],%%rsi		\n\t"\
		"vmovaps	(%%rsi),	%%xmm11		\n\t"\
		"vpsubd		%%xmm8 ,%%xmm11,%%xmm11	\n\t"\
		"vmovmskps	%%xmm11,	%%rsi	\n\t"\
		"\n\t"\
		"movq	%[__n_minus_silp1],%%rcx	\n\t"\
		"vmovd	0xC(%%rcx)	,%%xmm9 	\n\t"/* .d3 term of index quartet */\
		"vpshufd	$0,	%%xmm9 ,%%xmm9 		\n\t"\
		"vpsubd		%%xmm8 ,%%xmm9 ,%%xmm9 		\n\t"\
		"vmovmskps	%%xmm9 ,%%rcx		\n\t"\
		"\n\t"\
		"movq	%[__sinwtm1]	,%%rdx		\n\t"\
		"vmovd	0xC(%%rdx)	,%%xmm10	\n\t"/* .d3 term of index quartet */\
		"vpshufd	$0,	%%xmm10,%%xmm10		\n\t"\
		"vmovaps		%%xmm8 ,%%xmm11		\n\t"\
		"vpsubd		%%xmm10,%%xmm11,%%xmm11		\n\t"\
		"vmovmskps	%%xmm11,%%rdx		\n\t"\
		"\n\t"\
	"shlq	$5,%%rsi	\n\t"/* i0123 */\
	"shlq	$5,%%rcx	\n\t"/* m0123 */\
	"shlq	$5,%%rdx	\n\t"/* n0123 */\
		"\n\t"\
		"addq	%%rdi,%%rcx		\n\t"\
		"addq	%%rdi,%%rdx		\n\t"\
		"vbroadcastsd 0x870(%%rdi),%%ymm9 	\n\t"/* wtlp1 */\
		"vbroadcastsd 0x878(%%rdi),%%ymm10	\n\t"/* wtnm1 */\
		"vmulpd	%%ymm12,%%ymm9 ,%%ymm9 		\n\t"\
		"vmulpd	%%ymm13,%%ymm10,%%ymm10		\n\t"\
		"vmulpd	     (%%rcx),%%ymm9 ,%%ymm9 \n\t"\
		"vmulpd	0x200(%%rdx),%%ymm10,%%ymm10\n\t"\
		"\n\t"\
		"vmulpd		%%ymm10,%%ymm7,%%ymm7	\n\t"\
		"vaddpd		%%ymm14,%%ymm7,%%ymm7	\n\t"\
		"vmovaps	%%ymm7,%%ymm10			\n\t"\
		"vroundpd	$0,%%ymm10,%%ymm10		\n\t"\
		"\n\t"\
		"movq	%[__sign_mask],%%rbx		\n\t"\
		"vsubpd		%%ymm10,%%ymm7,%%ymm7	\n\t"\
		"vandpd		(%%rbx),%%ymm7,%%ymm7	\n\t"\
		"vmaxpd		%%ymm15,%%ymm7,%%ymm15	\n\t"\
		/* cy   = DNINT(temp*baseinv[i1]): */\
		"addq	%%rdi,%%rsi		\n\t"\
		"vmovaps	%%ymm10	,%%ymm7			\n\t"\
		"vmulpd	0x600(%%rsi),%%ymm10,%%ymm10\n\t"\
		"vroundpd	$0,%%ymm10,%%ymm14		\n\t"\
		"vmovaps	%%ymm14,%%ymm10				\n\t"/* cy = cpy cy_out */\
		/* x = (temp-cy*base[i1])*wt: */\
		"vmulpd	0x400(%%rsi),%%ymm10,%%ymm10	\n\t"/* cy*base[i0123] */\
		"vsubpd		%%ymm10,%%ymm7 ,%%ymm7 	\n\t"\
		"vmulpd		%%ymm9 ,%%ymm7 ,%%ymm7 	\n\t"\
		/* Update wts-array pointers in preparation for next call of the macro: */\
		"movq	%[__wtA]	,%%rax		\n\t"\
		"movq	%[__wtB]	,%%rbx		\n\t"\
		"movq	%[__wtC]	,%%rcx		\n\t"\
		"addq	$0x20	,%%rax			\n\t"/* add0 += 4 */\
		"subq	$0x20	,%%rbx			\n\t"/* add1 -= 4 */\
		"subq	$0x20	,%%rcx			\n\t"/* add2 -= 4 */\
		"movq	%%rax	,%[__wtA]		\n\t"\
		"movq	%%rbx	,%[__wtB]		\n\t"\
		"movq	%%rcx	,%[__wtC]		\n\t"\
		"\n\t"\
		/* Get ready for store of final-updated bjmod[0:3] values: */\
		"movq	%[__sse_bw]	,%%rax		\n\t"\
		"movq	%[__sse_nm1],%%rbx		\n\t"\
		"vpaddd		(%%rax),%%xmm8 ,%%xmm8 	\n\t"\
		"vpand		(%%rbx),%%xmm8 ,%%xmm8 	\n\t"\
		"\n\t"\
		/* Store bjmodn index quartet: */\
		"movq	%[__bjmod_0],%%rcx		\n\t"\
		"vmovaps	%%xmm8,(%%rcx)			\n\t"\
		/* Store cy_out: */\
		"movq		%[__cy],%%rbx	\n\t"\
		"vmovaps	%%ymm14,(%%rbx)	\n\t"/* ymm14 = cy_in */\
		/* Store maxerr: */\
		"movq		%[__half_arr],%%rdx		\n\t"\
		"vmovaps	%%ymm15,-0x40(%%rdx)	\n\t"\
	/* 4-way transpose of outputs (Re, Im parts separately): Inputs from ymm0-7. Outputs into r0/1,2/3,4/5.6/7: */	\
	/* Because default inputs for our 4 x 4 transpose macro (e.g. the one used at start of this carry macro) */\
	/* are into ymm4/2/8/2, munge inputs into that order, resolving name-conflicts via use of the now-available ymm8-15 for outputs: */\
		"movq		%[__data],%%rax			\n\t"\
		"vshufpd	$15,%%ymm2,%%ymm0,%%ymm10					\n\t		vshufpd	$15,%%ymm3,%%ymm1,%%ymm11						\n\t"\
		"vshufpd	$0 ,%%ymm2,%%ymm0,%%ymm0					\n\t		vshufpd	$0 ,%%ymm3,%%ymm1,%%ymm1						\n\t"\
		"vshufpd	$15,%%ymm6,%%ymm4,%%ymm12					\n\t		vshufpd	$15,%%ymm7,%%ymm5,%%ymm13						\n\t"\
		"vshufpd	$0 ,%%ymm6,%%ymm4,%%ymm4					\n\t		vshufpd	$0 ,%%ymm7,%%ymm5,%%ymm5						\n\t"\
		"vperm2f128 $32,%%ymm12,%%ymm10,%%ymm2 		/* Re B	*/	\n\t		vperm2f128 $32,%%ymm13,%%ymm11,%%ymm3		/* Im B	*/	\n\t"\
		"vperm2f128 $49,%%ymm12,%%ymm10,%%ymm10		/* Re D	*/	\n\t		vperm2f128 $49,%%ymm13,%%ymm11,%%ymm11		/* Im D	*/	\n\t"\
		"vperm2f128 $32,%%ymm4 ,%%ymm0 ,%%ymm12		/* Re A	*/	\n\t		vperm2f128 $32,%%ymm5 ,%%ymm1 ,%%ymm13 		/* Im A	*/	\n\t"\
		"vperm2f128 $49,%%ymm4 ,%%ymm0 ,%%ymm0 		/* Re C	*/	\n\t		vperm2f128 $49,%%ymm5 ,%%ymm1 ,%%ymm1		/* Im C	*/	\n\t"\
		/* And write 'em back to memory: */\
		"vmovaps	%%ymm12,     (%%rax)						\n\t		vmovaps	%%ymm13,0x020(%%rax)				\n\t"\
		"vmovaps	%%ymm2 ,0x040(%%rax)						\n\t		vmovaps	%%ymm3 ,0x060(%%rax)				\n\t"\
		"vmovaps	%%ymm0 ,0x080(%%rax)						\n\t		vmovaps	%%ymm1 ,0x0a0(%%rax)				\n\t"\
		"vmovaps	%%ymm10,0x0c0(%%rax)						\n\t		vmovaps	%%ymm11,0x0e0(%%rax)				\n\t"\
		:					/* outputs: none */\
		: [__data]		"m" (Xdata)	/* All inputs from memory addresses here */\
		, [__wtA]		"m" (XwtA)		\
		, [__wtB]		"m" (XwtB)		\
		, [__wtC]		"m" (XwtC)		\
		, [__cy]		"m" (Xcy)		\
		, [__bjmod_0]	"m" (Xbjmod_0)		\
		, [__half_arr]	"m" (Xhalf_arr)		\
		, [__n_minus_silp1] "m" (Xn_minus_silp1)\
		, [__n_minus_sil]	"m" (Xn_minus_sil)	\
		, [__sign_mask]	"m" (Xsign_mask)	\
		, [__sinwt]		"m" (Xsinwt)		\
		, [__sinwtm1]	"m" (Xsinwtm1)		\
		, [__sse_bw]	"m" (Xsse_bw)		\
		, [__sse_nm1]	"m" (Xsse_nm1)		\
		, [__sse_sw]	"m" (Xsse_sw)		\
		/* Prefetch: base address and 3 index offsets */\
		,	[__add0] "m" (Xadd0)\
		,	[__p1] "m" (Xp1)\
		,	[__p2] "m" (Xp2)\
		,	[__p3] "m" (Xp3)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r14","r15","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
	);\
	}

	/****************************************************************************************************************/
	/*** Non-power-of-2-FFT versions of above AVX_cmplx_carry_norm_pow2_errcheck0|1_X4 Mersenne-mod carry macros: ***/
	/****************************************************************************************************************/

	#define AVX_cmplx_carry_norm_errcheck0_X4(Xdata,XwtA,XwtB,XwtC,Xcy,Xbjmod_0,Xhalf_arr,Xi,Xn_minus_silp1,Xn_minus_sil,Xsign_mask,Xsinwt,Xsinwtm1,Xsse_bw,Xsse_n,Xsse_sw, Xadd0,Xp1,Xp2,Xp3)\
	{\
	__asm__ volatile (\
	"movq	%[__add0],%%r14	\n\t"/* base address for 4 prefetches-from-main-data-array spread through this macro */\
		"movq		%[__data],%%rax		\n\t"\
	/* 4-way transpose of inputs (Re, Im parts separately): Inputs from r0/1,2/3,4/5.6/7. Outputs into ymm0-7: */\
	/* Real parts use ymm0,2,4,6, ymm8 as tmp-reg:					Imag parts use ymm1,3,5,7, ymm9 as tm-reg: */\
		"vmovaps	     (%%rax),%%ymm4						\n\t		vmovaps	0x020(%%rax),%%ymm5							\n\t"\
		"vmovaps	0x040(%%rax),%%ymm2						\n\t		vmovaps	0x060(%%rax),%%ymm3							\n\t"\
		"vshufpd	$15,%%ymm2,%%ymm4,%%ymm6				\n\t		vshufpd	$15,%%ymm3,%%ymm5,%%ymm7					\n\t"\
		"vshufpd	$0 ,%%ymm2,%%ymm4,%%ymm4				\n\t		vshufpd	$0 ,%%ymm3,%%ymm5,%%ymm5					\n\t"\
		"vmovaps	0x080(%%rax),%%ymm8						\n\t		vmovaps	0x0a0(%%rax),%%ymm9							\n\t"\
		"vmovaps	0x0c0(%%rax),%%ymm2						\n\t		vmovaps	0x0e0(%%rax),%%ymm3							\n\t"\
		"vshufpd	$15,%%ymm2,%%ymm8,%%ymm0				\n\t		vshufpd	$15,%%ymm3,%%ymm9,%%ymm1					\n\t"\
		"vshufpd	$0 ,%%ymm2,%%ymm8,%%ymm8				\n\t		vshufpd	$0 ,%%ymm3,%%ymm9,%%ymm9					\n\t"\
		"vperm2f128 $32,%%ymm0,%%ymm6,%%ymm2	/* Re B	*/	\n\t		vperm2f128 $32,%%ymm1,%%ymm7,%%ymm3		/* Im B	*/	\n\t"\
		"vperm2f128 $49,%%ymm0,%%ymm6,%%ymm6	/* Re D	*/	\n\t		vperm2f128 $49,%%ymm1,%%ymm7,%%ymm7		/* Im D	*/	\n\t"\
		"vperm2f128 $32,%%ymm8,%%ymm4,%%ymm0	/* Re A	*/	\n\t		vperm2f128 $32,%%ymm9,%%ymm5,%%ymm1 	/* Im A	*/	\n\t"\
		"vperm2f128 $49,%%ymm8,%%ymm4,%%ymm4	/* Re C	*/	\n\t		vperm2f128 $49,%%ymm9,%%ymm5,%%ymm5 	/* Im C	*/	\n\t"\
		/* Won't need main-array again until output transpose, so re-use rax for half_arr */\
		"movq	 %[__half_arr],%%rax	\n\t"\
		/* half_arr + 16*[0,1,2,3] = [wt,wt_inv,base,baseinv] */\
		"movq		%[__cy],%%rbx	\n\t"\
		"vmovaps	(%%rbx),%%ymm14	\n\t"/* ymm14 = cy_in */\
		"vmovaps	-0x40(%%rax),%%ymm15	\n\t"/* ymm15 = maxerr */\
	/**********************************/\
	/* Do A.re-quartet: Data in ymm0: */\
	/**********************************/\
	"prefetcht0	(%%r14)		\n\t"\
		"movq	%[__bjmod_0],%%rsi		\n\t"\
		"vmovaps	(%%rsi),	%%xmm8 		\n\t"/* bjmod[0:3]. PERSISTENT COPY OF BJMOD[0:3] REMAINS IN XMM8. */\
		"movq	%[__sse_sw],%%rsi		\n\t"\
		"vmovaps	(%%rsi),	%%xmm11		\n\t"/* sw[0:3] */\
		"vpsubd		%%xmm8 ,%%xmm11,%%xmm11	\n\t"/* sw[0:3] - bjmod[0:3] */\
		"vmovmskps	%%xmm11,	%%rsi	\n\t"/* Extract sign bits into 4-bit signmask <i3|i2|i1|i0>; idxs into base/inv table */\
		"movslq	%[__i]	,%%rcx			\n\t"/* i0=i for first block: */\
	"andq	$0xfffffffffffffffe,%%rsi	\n\t"/* Mask off lowest bit */\
		"addq	%%rcx	,%%rsi			\n\t"\
		"\n\t"\
		"movq	%[__n_minus_sil],%%rcx	\n\t"\
		"vmovd	0x0(%%rcx)	,%%xmm9 	\n\t"/* n_minus_sil in low 32 bits of xmm9  */\
		"vpshufd	$0,	%%xmm9 ,%%xmm9 		\n\t"/* Broadcast low 32 bits of xmm9  to all 4 slots of xmm9  */\
		"vpsubd		%%xmm8 ,%%xmm9 ,%%xmm9 		\n\t"/* n_minus_sil - bjmod[0:3] */\
		"vmovmskps	%%xmm9 ,%%rcx		\n\t"/* Extract sign bits into 4-bit signmask <m3|m2|m1|m0>; idxs into base/inv tables -> byte[2] of ecx... */\
		"\n\t"\
		"movq	%[__sinwt]	,%%rdx		\n\t"\
		"vmovd	0x0(%%rdx)	,%%xmm10		\n\t"/* sinwt in low 32 bits of xmm10 */\
		"vpshufd	$0,	%%xmm10,%%xmm10		\n\t"/* Broadcast low 32 bits of xmm10 to all 4 slots of xmm10 */\
		"vmovaps		%%xmm8 ,%%xmm11		\n\t"/* xmm11 = bjmod[0:3] copy */\
		"vpsubd		%%xmm10,%%xmm11,%%xmm11		\n\t"/* bjmod[0:3] - sinwt */\
		"vmovmskps	%%xmm11,%%rdx		\n\t"/* Extract sign bits into 4-bit signmask <n3|n2|n1|n0>; idxs into base/inv tables -> byte[1] of edx... */\
		"\n\t"\
	"shlq	$5,%%rsi	\n\t"/* i0123 << 5 (= lshift to give ptr offset for ymm-size data */\
	"shlq	$5,%%rcx	\n\t"/* m0123 << 5 (= lshift to give ptr offset for ymm-size data */\
	"shlq	$5,%%rdx	\n\t"/* n0123 << 5 (= lshift to give ptr offset for ymm-size data */\
		"\n\t"\
		"movq	%[__half_arr],%%rdi	\n\t"\
		"movq	%[__wtA]	,%%rax		\n\t"\
		"movq	%[__wtB]	,%%rbx		\n\t"\
		"vmovaps	     (%%rax),%%ymm12 	\n\t"/* wtA[j  ]; ebx FREE */\
		"vmovups	-0x10(%%rbx),%%ymm13	\n\t"/* wtB[j-1]; load doubles from rcx+[-0x10,-0x08, 0, +0x08] - note this address is ymm (16-byte) but not ymm (32-byte) aligned */\
		/* reverse-running indexing used for inv-wts really means we need to reverse ordering of 4 doubles d0-3 in ymm13 */\
		"vshufpd	$5,%%ymm13,%%ymm13,%%ymm13	\n\t"/* d[0123] -> d[1032] */\
		"vperm2f128 $1,%%ymm13,%%ymm13,%%ymm13	\n\t"/* d[1032] -> d[3210] */\
		"\n\t"\
		"addq	%%rdi,%%rcx		\n\t"\
		"addq	%%rdi,%%rdx		\n\t"\
		/* SSE2 version has double-copies in wtl/wtn ... AVX replaces redundant-data loads with load-with-broadcast: */\
		"vbroadcastsd 0x800(%%rdi),%%ymm9 	\n\t"/* wtl */\
		"vbroadcastsd 0x808(%%rdi),%%ymm10	\n\t"/* wtn */\
		"vmulpd	%%ymm12,%%ymm9 ,%%ymm9 		\n\t"/* wt   =wtA*wtl */\
		"vmulpd	%%ymm13,%%ymm10,%%ymm10		\n\t"/* wtinv=wtB*wtn */\
		"vmulpd	     (%%rcx),%%ymm9 ,%%ymm9 	\n\t"/* wt    *= one_half[m0123] */\
		"vmulpd	0x200(%%rdx),%%ymm10,%%ymm10	\n\t"/* wtinv *= one_half[16+n0123] */\
		"\n\t"\
		"vmulpd		%%ymm10,%%ymm0,%%ymm0		\n\t"/* x *= wtinv; ymm10 FREE */\
		"vaddpd		%%ymm14,%%ymm0,%%ymm0		\n\t"/* x *= wtinv + cy */\
		"vmovaps	%%ymm0,%%ymm10		\n\t"/* temp = x */\
		"vroundpd	$0,%%ymm10,%%ymm10	\n\t"/* temp = DNINT(x) */\
		"\n\t"\
		"movq	%[__sign_mask],%%rbx	\n\t"\
		"vsubpd		%%ymm10,%%ymm0,%%ymm0	\n\t"/* x - temp */\
		"vandpd		(%%rbx),%%ymm0,%%ymm0	\n\t"/* frac = fabs(x-temp) */\
		"vmaxpd		%%ymm15,%%ymm0,%%ymm15	\n\t"/* if(frac > maxerr) maxerr=frac */\
		/* cy   = DNINT(temp*baseinv[i1]): */\
		"addq	%%rdi,%%rsi		\n\t"\
		"vmovaps	%%ymm10	,%%ymm0				\n\t"/* cpy temp */\
		"vmulpd	0x600(%%rsi),%%ymm10,%%ymm10	\n\t"/* temp*baseinv[i0123] */\
		"vroundpd	$0,%%ymm10,%%ymm14			\n\t"/* cy_out */\
		"vmovaps	%%ymm14,%%ymm10				\n\t"/* cy = cpy cy_out */\
		/* x = (temp-cy*base[i1])*wt: */\
		"vmulpd	0x400(%%rsi),%%ymm10,%%ymm10	\n\t"/* cy*base[i0123] */\
		"vsubpd		%%ymm10,%%ymm0 ,%%ymm0 		\n\t"/* x = (temp-cy*base[i0123]) */\
		"vmulpd		%%ymm9 ,%%ymm0 ,%%ymm0 		\n\t"/* x *= wt */\
		/* Get ready for next set [IM0~] by computing bjmod[0:3] += bw (mod n): */\
		"movq	%[__sse_n]	,%%rbx		\n\t"\
		"vmovaps		(%%rbx)	,%%xmm10	\n\t"\
		"movq	%[__sse_bw]	,%%rax		\n\t"\
		"vpaddd		(%%rax),%%xmm8 ,%%xmm8 	\n\t"\
		"vmovaps	%%xmm8 ,%%xmm9 		\n\t"\
		"vpcmpgtd	%%xmm10,%%xmm9 ,%%xmm9 	\n\t"\
		"vpand		%%xmm10,%%xmm9 ,%%xmm9 	\n\t"\
		"vpsubd		%%xmm9 ,%%xmm8 ,%%xmm8 	\n\t"\
		"\n\t"\
	/**********************************/\
	/* Do A.im-quartet: Data in ymm1: */\
	/**********************************/\
		"movq	%[__sse_sw],%%rsi		\n\t"\
		"vmovaps	(%%rsi),	%%xmm11		\n\t"/* sw[0:3] */\
		"vpsubd		%%xmm8 ,%%xmm11,%%xmm11	\n\t"/* sw[0:3] - bjmod[0:3] */\
		"vmovmskps	%%xmm11,	%%rsi	\n\t"/* Extract sign bits into 4-bit signmask <i3|i2|i1|i0>; idxs into base/inv table */\
		"\n\t"\
		"movq	%[__n_minus_silp1],%%rcx	\n\t"\
		"vmovd	0x0(%%rcx)	,%%xmm9 	\n\t"/* n_minus_sil in low 32 bits of xmm9  */\
		"vpshufd	$0,	%%xmm9 ,%%xmm9 		\n\t"/* Broadcast low 32 bits of xmm9  to all 4 slots of xmm9  */\
		"vpsubd		%%xmm8 ,%%xmm9 ,%%xmm9 		\n\t"/* n_minus_sil - bjmod[0:3] */\
		"vmovmskps	%%xmm9 ,%%rcx		\n\t"/* Extract sign bits into 4-bit signmask <m3|m2|m1|m0>; idxs into base/inv tables -> byte[2] of ecx... */\
		"\n\t"\
		"movq	%[__sinwtm1]	,%%rdx		\n\t"\
		"vmovd	0x0(%%rdx)	,%%xmm10		\n\t"/* sinwt in low 32 bits of xmm10 */\
		"vpshufd	$0,	%%xmm10,%%xmm10		\n\t"/* Broadcast low 32 bits of xmm10 to all 4 slots of xmm10 */\
		"vmovaps		%%xmm8 ,%%xmm11		\n\t"/* xmm11 = bjmod[0:3] copy */\
		"vpsubd		%%xmm10,%%xmm11,%%xmm11		\n\t"/* bjmod[0:3] - sinwt */\
		"vmovmskps	%%xmm11,%%rdx		\n\t"/* Extract sign bits into 4-bit signmask <n3|n2|n1|n0>; idxs into base/inv tables -> byte[1] of edx... */\
		"\n\t"\
	"shlq	$5,%%rsi	\n\t"/* i0123 << 5 (= lshift to give ptr offset for ymm-size data */\
	"shlq	$5,%%rcx	\n\t"/* m0123 << 5 (= lshift to give ptr offset for ymm-size data */\
	"shlq	$5,%%rdx	\n\t"/* n0123 << 5 (= lshift to give ptr offset for ymm-size data */\
		"\n\t"\
		"movq	%[__wtC]	,%%rbx		\n\t"/* wtA unchanged; wtB == wtC for remaining 7 of 8 sets of carries */\
		"vmovups	-0x10(%%rbx),%%ymm13	\n\t"/* wtC[j-1]; load doubles from rcx+[-0x10,-0x08, 0, +0x08] - note this address is ymm (16-byte) but not ymm (32-byte) aligned */\
		"vshufpd	$5,%%ymm13,%%ymm13,%%ymm13	\n\t"/* d[0123] -> d[1032] */\
		"vperm2f128 $1,%%ymm13,%%ymm13,%%ymm13	\n\t"/* d[1032] -> d[3210] */\
		"\n\t"\
		"addq	%%rdi,%%rcx		\n\t"\
		"addq	%%rdi,%%rdx		\n\t"\
		/* SSE2 version has double-copies in wtl/wtn ... AVX replaces redundant-data loads with load-with-broadcast: */\
		"vbroadcastsd 0x810(%%rdi),%%ymm9 	\n\t"/* wtlp1 */\
		"vbroadcastsd 0x818(%%rdi),%%ymm10	\n\t"/* wtnm1 */\
		"vmulpd	%%ymm12,%%ymm9 ,%%ymm9 		\n\t"/* wt   =wtA*wtlp1 */\
		"vmulpd	%%ymm13,%%ymm10,%%ymm10		\n\t"/* wtinv=wtC*wtnm1 */\
		"vmulpd	     (%%rcx),%%ymm9 ,%%ymm9 	\n\t"/* wt    *= one_half[m0123] */\
		"vmulpd	0x200(%%rdx),%%ymm10,%%ymm10	\n\t"/* wtinv *= one_half[16+n0123] */\
		"\n\t"\
		"vmulpd		%%ymm10,%%ymm1,%%ymm1		\n\t"/* x *= wtinv; ymm10 FREE */\
		"vaddpd		%%ymm14,%%ymm1,%%ymm1		\n\t"/* x *= wtinv + cy */\
		"vmovaps	%%ymm1,%%ymm10		\n\t"/* temp = x */\
		"vroundpd	$0,%%ymm10,%%ymm10	\n\t"/* temp = DNINT(x) */\
		"\n\t"\
		"movq	%[__sign_mask],%%rbx	\n\t"\
		"vsubpd		%%ymm10,%%ymm1,%%ymm1	\n\t"/* x - temp */\
		"vandpd		(%%rbx),%%ymm1,%%ymm1	\n\t"/* frac = fabs(x-temp) */\
		"vmaxpd		%%ymm15,%%ymm1,%%ymm15	\n\t"/* if(frac > maxerr) maxerr=frac */\
		/* cy   = DNINT(temp*baseinv[i1]): */\
		"addq	%%rdi,%%rsi		\n\t"\
		"vmovaps	%%ymm10	,%%ymm1				\n\t"/* cpy temp */\
		"vmulpd	0x600(%%rsi),%%ymm10,%%ymm10	\n\t"/* temp*baseinv[i0123] */\
		"vroundpd	$0,%%ymm10,%%ymm14			\n\t"/* cy_out */\
		"vmovaps	%%ymm14,%%ymm10				\n\t"/* cy = cpy cy_out */\
		/* x = (temp-cy*base[i1])*wt: */\
		"vmulpd	0x400(%%rsi),%%ymm10,%%ymm10	\n\t"/* cy*base[i0123] */\
		"vsubpd		%%ymm10,%%ymm1 ,%%ymm1 		\n\t"/* x = (temp-cy*base[i0123]) */\
		"vmulpd		%%ymm9 ,%%ymm1 ,%%ymm1 		\n\t"/* x *= wt */\
		/* Get ready for next set [RE1~] by computing bjmod[0:3] += bw (mod n): */\
		"movq	%[__sse_n]	,%%rbx		\n\t"\
		"vmovaps		(%%rbx)	,%%xmm10	\n\t"\
		"movq	%[__sse_bw]	,%%rax		\n\t"\
		"vpaddd		(%%rax),%%xmm8 ,%%xmm8 	\n\t"\
		"vmovaps	%%xmm8 ,%%xmm9 		\n\t"\
		"vpcmpgtd	%%xmm10,%%xmm9 ,%%xmm9 	\n\t"\
		"vpand		%%xmm10,%%xmm9 ,%%xmm9 	\n\t"\
		"vpsubd		%%xmm9 ,%%xmm8 ,%%xmm8 	\n\t"\
		"\n\t"\
	/**********************************/\
	/* Do B.re-quartet: Data in ymm2: */\
	/**********************************/\
	"movslq		%[__p1],%%r15	\n\t"	\
	"prefetcht0	(%%r14,%%r15,8)	\n\t"\
		"movq	%[__sse_sw],%%rsi		\n\t"\
		"vmovaps	(%%rsi),	%%xmm11		\n\t"\
		"vpsubd		%%xmm8 ,%%xmm11,%%xmm11	\n\t"\
		"vmovmskps	%%xmm11,	%%rsi	\n\t"\
		"\n\t"\
		"movq	%[__n_minus_sil],%%rcx	\n\t"\
		"vmovd	0x4(%%rcx)	,%%xmm9 	\n\t"/* .d1 term of index quartet */\
		"vpshufd	$0,	%%xmm9 ,%%xmm9 		\n\t"\
		"vpsubd		%%xmm8 ,%%xmm9 ,%%xmm9 		\n\t"\
		"vmovmskps	%%xmm9 ,%%rcx		\n\t"\
		"\n\t"\
		"movq	%[__sinwt]	,%%rdx		\n\t"\
		"vmovd	0x4(%%rdx)	,%%xmm10	\n\t"/* .d1 term of index quartet */\
		"vpshufd	$0,	%%xmm10,%%xmm10		\n\t"\
		"vmovaps		%%xmm8 ,%%xmm11		\n\t"\
		"vpsubd		%%xmm10,%%xmm11,%%xmm11		\n\t"\
		"vmovmskps	%%xmm11,%%rdx		\n\t"\
		"\n\t"\
	"shlq	$5,%%rsi	\n\t"/* i0123 */\
	"shlq	$5,%%rcx	\n\t"/* m0123 */\
	"shlq	$5,%%rdx	\n\t"/* n0123 */\
		"\n\t"\
		"addq	%%rdi,%%rcx		\n\t"\
		"addq	%%rdi,%%rdx		\n\t"\
		"vbroadcastsd 0x820(%%rdi),%%ymm9 	\n\t"/* wtl */\
		"vbroadcastsd 0x828(%%rdi),%%ymm10	\n\t"/* wtn */\
		"vmulpd	%%ymm12,%%ymm9 ,%%ymm9 		\n\t"\
		"vmulpd	%%ymm13,%%ymm10,%%ymm10		\n\t"\
		"vmulpd	     (%%rcx),%%ymm9 ,%%ymm9 \n\t"\
		"vmulpd	0x200(%%rdx),%%ymm10,%%ymm10\n\t"\
		"\n\t"\
		"vmulpd		%%ymm10,%%ymm2,%%ymm2	\n\t"\
		"vaddpd		%%ymm14,%%ymm2,%%ymm2	\n\t"\
		"vmovaps	%%ymm2,%%ymm10			\n\t"\
		"vroundpd	$0,%%ymm10,%%ymm10		\n\t"\
		"\n\t"\
		"movq	%[__sign_mask],%%rbx		\n\t"\
		"vsubpd		%%ymm10,%%ymm2,%%ymm2	\n\t"\
		"vandpd		(%%rbx),%%ymm2,%%ymm2	\n\t"\
		"vmaxpd		%%ymm15,%%ymm2,%%ymm15	\n\t"\
		/* cy   = DNINT(temp*baseinv[i1]): */\
		"addq	%%rdi,%%rsi		\n\t"\
		"vmovaps	%%ymm10	,%%ymm2			\n\t"\
		"vmulpd	0x600(%%rsi),%%ymm10,%%ymm10\n\t"\
		"vroundpd	$0,%%ymm10,%%ymm14		\n\t"\
		"vmovaps	%%ymm14,%%ymm10				\n\t"/* cy = cpy cy_out */\
		/* x = (temp-cy*base[i1])*wt: */\
		"vmulpd	0x400(%%rsi),%%ymm10,%%ymm10	\n\t"/* cy*base[i0123] */\
		"vsubpd		%%ymm10,%%ymm2 ,%%ymm2 	\n\t"\
		"vmulpd		%%ymm9 ,%%ymm2 ,%%ymm2 	\n\t"\
		/* Get ready for next set [IM1~] by computing bjmod[0:3] += bw (mod n): */\
		"movq	%[__sse_n]	,%%rbx		\n\t"\
		"vmovaps		(%%rbx)	,%%xmm10	\n\t"\
		"movq	%[__sse_bw]	,%%rax		\n\t"\
		"vpaddd		(%%rax),%%xmm8 ,%%xmm8 	\n\t"\
		"vmovaps	%%xmm8 ,%%xmm9 		\n\t"\
		"vpcmpgtd	%%xmm10,%%xmm9 ,%%xmm9 	\n\t"\
		"vpand		%%xmm10,%%xmm9 ,%%xmm9 	\n\t"\
		"vpsubd		%%xmm9 ,%%xmm8 ,%%xmm8 	\n\t"\
		"\n\t"\
	/**********************************/\
	/* Do B.im-quartet: Data in ymm3: */\
	/**********************************/\
		"movq	%[__sse_sw],%%rsi		\n\t"\
		"vmovaps	(%%rsi),	%%xmm11		\n\t"\
		"vpsubd		%%xmm8 ,%%xmm11,%%xmm11	\n\t"\
		"vmovmskps	%%xmm11,	%%rsi	\n\t"\
		"\n\t"\
		"movq	%[__n_minus_silp1],%%rcx	\n\t"\
		"vmovd	0x4(%%rcx)	,%%xmm9 	\n\t"/* .d1 term of index quartet */\
		"vpshufd	$0,	%%xmm9 ,%%xmm9 		\n\t"\
		"vpsubd		%%xmm8 ,%%xmm9 ,%%xmm9 		\n\t"\
		"vmovmskps	%%xmm9 ,%%rcx		\n\t"\
		"\n\t"\
		"movq	%[__sinwtm1]	,%%rdx		\n\t"\
		"vmovd	0x4(%%rdx)	,%%xmm10	\n\t"/* .d1 term of index quartet */\
		"vpshufd	$0,	%%xmm10,%%xmm10		\n\t"\
		"vmovaps		%%xmm8 ,%%xmm11		\n\t"\
		"vpsubd		%%xmm10,%%xmm11,%%xmm11		\n\t"\
		"vmovmskps	%%xmm11,%%rdx		\n\t"\
		"\n\t"\
	"shlq	$5,%%rsi	\n\t"/* i0123 */\
	"shlq	$5,%%rcx	\n\t"/* m0123 */\
	"shlq	$5,%%rdx	\n\t"/* n0123 */\
		"\n\t"\
		"addq	%%rdi,%%rcx		\n\t"\
		"addq	%%rdi,%%rdx		\n\t"\
		"vbroadcastsd 0x830(%%rdi),%%ymm9 	\n\t"/* wtlp1 */\
		"vbroadcastsd 0x838(%%rdi),%%ymm10	\n\t"/* wtnm1 */\
		"vmulpd	%%ymm12,%%ymm9 ,%%ymm9 		\n\t"\
		"vmulpd	%%ymm13,%%ymm10,%%ymm10		\n\t"\
		"vmulpd	     (%%rcx),%%ymm9 ,%%ymm9 \n\t"\
		"vmulpd	0x200(%%rdx),%%ymm10,%%ymm10\n\t"\
		"\n\t"\
		"vmulpd		%%ymm10,%%ymm3,%%ymm3	\n\t"\
		"vaddpd		%%ymm14,%%ymm3,%%ymm3	\n\t"\
		"vmovaps	%%ymm3,%%ymm10			\n\t"\
		"vroundpd	$0,%%ymm10,%%ymm10		\n\t"\
		"\n\t"\
		"movq	%[__sign_mask],%%rbx		\n\t"\
		"vsubpd		%%ymm10,%%ymm3,%%ymm3	\n\t"\
		"vandpd		(%%rbx),%%ymm3,%%ymm3	\n\t"\
		"vmaxpd		%%ymm15,%%ymm3,%%ymm15	\n\t"\
		/* cy   = DNINT(temp*baseinv[i1]): */\
		"addq	%%rdi,%%rsi		\n\t"\
		"vmovaps	%%ymm10	,%%ymm3			\n\t"\
		"vmulpd	0x600(%%rsi),%%ymm10,%%ymm10\n\t"\
		"vroundpd	$0,%%ymm10,%%ymm14		\n\t"\
		"vmovaps	%%ymm14,%%ymm10				\n\t"/* cy = cpy cy_out */\
		/* x = (temp-cy*base[i1])*wt: */\
		"vmulpd	0x400(%%rsi),%%ymm10,%%ymm10	\n\t"/* cy*base[i0123] */\
		"vsubpd		%%ymm10,%%ymm3 ,%%ymm3 	\n\t"\
		"vmulpd		%%ymm9 ,%%ymm3 ,%%ymm3 	\n\t"\
		/* Get ready for next set [RE2~] by computing bjmod[0:3] += bw (mod n): */\
		"movq	%[__sse_n]	,%%rbx		\n\t"\
		"vmovaps		(%%rbx)	,%%xmm10	\n\t"\
		"movq	%[__sse_bw]	,%%rax		\n\t"\
		"vpaddd		(%%rax),%%xmm8 ,%%xmm8 	\n\t"\
		"vmovaps	%%xmm8 ,%%xmm9 		\n\t"\
		"vpcmpgtd	%%xmm10,%%xmm9 ,%%xmm9 	\n\t"\
		"vpand		%%xmm10,%%xmm9 ,%%xmm9 	\n\t"\
		"vpsubd		%%xmm9 ,%%xmm8 ,%%xmm8 	\n\t"\
		"\n\t"\
	/**********************************/\
	/* Do C.re-quartet: Data in ymm4: */\
	/**********************************/\
	"movslq		%[__p2],%%r15	\n\t"	\
	"prefetcht0	(%%r14,%%r15,8)	\n\t"\
		"movq	%[__sse_sw],%%rsi		\n\t"\
		"vmovaps	(%%rsi),	%%xmm11		\n\t"\
		"vpsubd		%%xmm8 ,%%xmm11,%%xmm11	\n\t"\
		"vmovmskps	%%xmm11,	%%rsi	\n\t"\
		"\n\t"\
		"movq	%[__n_minus_sil],%%rcx	\n\t"\
		"vmovd	0x8(%%rcx)	,%%xmm9 	\n\t"/* .d2 term of index quartet */\
		"vpshufd	$0,	%%xmm9 ,%%xmm9 		\n\t"\
		"vpsubd		%%xmm8 ,%%xmm9 ,%%xmm9 		\n\t"\
		"vmovmskps	%%xmm9 ,%%rcx		\n\t"\
		"\n\t"\
		"movq	%[__sinwt]	,%%rdx		\n\t"\
		"vmovd	0x8(%%rdx)	,%%xmm10	\n\t"/* .d2 term of index quartet */\
		"vpshufd	$0,	%%xmm10,%%xmm10		\n\t"\
		"vmovaps		%%xmm8 ,%%xmm11		\n\t"\
		"vpsubd		%%xmm10,%%xmm11,%%xmm11		\n\t"\
		"vmovmskps	%%xmm11,%%rdx		\n\t"\
		"\n\t"\
	"shlq	$5,%%rsi	\n\t"/* i0123 */\
	"shlq	$5,%%rcx	\n\t"/* m0123 */\
	"shlq	$5,%%rdx	\n\t"/* n0123 */\
		"\n\t"\
		"addq	%%rdi,%%rcx		\n\t"\
		"addq	%%rdi,%%rdx		\n\t"\
		"vbroadcastsd 0x840(%%rdi),%%ymm9 	\n\t"/* wtl */\
		"vbroadcastsd 0x848(%%rdi),%%ymm10	\n\t"/* wtn */\
		"vmulpd	%%ymm12,%%ymm9 ,%%ymm9 		\n\t"\
		"vmulpd	%%ymm13,%%ymm10,%%ymm10		\n\t"\
		"vmulpd	     (%%rcx),%%ymm9 ,%%ymm9 \n\t"\
		"vmulpd	0x200(%%rdx),%%ymm10,%%ymm10\n\t"\
		"\n\t"\
		"vmulpd		%%ymm10,%%ymm4,%%ymm4	\n\t"\
		"vaddpd		%%ymm14,%%ymm4,%%ymm4	\n\t"\
		"vmovaps	%%ymm4,%%ymm10			\n\t"\
		"vroundpd	$0,%%ymm10,%%ymm10		\n\t"\
		"\n\t"\
		"movq	%[__sign_mask],%%rbx		\n\t"\
		"vsubpd		%%ymm10,%%ymm4,%%ymm4	\n\t"\
		"vandpd		(%%rbx),%%ymm4,%%ymm4	\n\t"\
		"vmaxpd		%%ymm15,%%ymm4,%%ymm15	\n\t"\
		/* cy   = DNINT(temp*baseinv[i1]): */\
		"addq	%%rdi,%%rsi		\n\t"\
		"vmovaps	%%ymm10	,%%ymm4			\n\t"\
		"vmulpd	0x600(%%rsi),%%ymm10,%%ymm10\n\t"\
		"vroundpd	$0,%%ymm10,%%ymm14		\n\t"\
		"vmovaps	%%ymm14,%%ymm10				\n\t"/* cy = cpy cy_out */\
		/* x = (temp-cy*base[i1])*wt: */\
		"vmulpd	0x400(%%rsi),%%ymm10,%%ymm10	\n\t"/* cy*base[i0123] */\
		"vsubpd		%%ymm10,%%ymm4 ,%%ymm4 	\n\t"\
		"vmulpd		%%ymm9 ,%%ymm4 ,%%ymm4 	\n\t"\
		/* Get ready for next set [IM2~] by computing bjmod[0:3] += bw (mod n): */\
		"movq	%[__sse_n]	,%%rbx		\n\t"\
		"vmovaps		(%%rbx)	,%%xmm10	\n\t"\
		"movq	%[__sse_bw]	,%%rax		\n\t"\
		"vpaddd		(%%rax),%%xmm8 ,%%xmm8 	\n\t"\
		"vmovaps	%%xmm8 ,%%xmm9 		\n\t"\
		"vpcmpgtd	%%xmm10,%%xmm9 ,%%xmm9 	\n\t"\
		"vpand		%%xmm10,%%xmm9 ,%%xmm9 	\n\t"\
		"vpsubd		%%xmm9 ,%%xmm8 ,%%xmm8 	\n\t"\
		"\n\t"\
	/**********************************/\
	/* Do C.im-quartet: Data in ymm5: */\
	/**********************************/\
		"movq	%[__sse_sw],%%rsi		\n\t"\
		"vmovaps	(%%rsi),	%%xmm11		\n\t"\
		"vpsubd		%%xmm8 ,%%xmm11,%%xmm11	\n\t"\
		"vmovmskps	%%xmm11,	%%rsi	\n\t"\
		"\n\t"\
		"movq	%[__n_minus_silp1],%%rcx	\n\t"\
		"vmovd	0x8(%%rcx)	,%%xmm9 	\n\t"/* .d2 term of index quartet */\
		"vpshufd	$0,	%%xmm9 ,%%xmm9 		\n\t"\
		"vpsubd		%%xmm8 ,%%xmm9 ,%%xmm9 		\n\t"\
		"vmovmskps	%%xmm9 ,%%rcx		\n\t"\
		"\n\t"\
		"movq	%[__sinwtm1]	,%%rdx		\n\t"\
		"vmovd	0x8(%%rdx)	,%%xmm10	\n\t"/* .d2 term of index quartet */\
		"vpshufd	$0,	%%xmm10,%%xmm10		\n\t"\
		"vmovaps		%%xmm8 ,%%xmm11		\n\t"\
		"vpsubd		%%xmm10,%%xmm11,%%xmm11		\n\t"\
		"vmovmskps	%%xmm11,%%rdx		\n\t"\
		"\n\t"\
	"shlq	$5,%%rsi	\n\t"/* i0123 */\
	"shlq	$5,%%rcx	\n\t"/* m0123 */\
	"shlq	$5,%%rdx	\n\t"/* n0123 */\
		"\n\t"\
		"addq	%%rdi,%%rcx		\n\t"\
		"addq	%%rdi,%%rdx		\n\t"\
		"vbroadcastsd 0x850(%%rdi),%%ymm9 	\n\t"/* wtlp1 */\
		"vbroadcastsd 0x858(%%rdi),%%ymm10	\n\t"/* wtnm1 */\
		"vmulpd	%%ymm12,%%ymm9 ,%%ymm9 		\n\t"\
		"vmulpd	%%ymm13,%%ymm10,%%ymm10		\n\t"\
		"vmulpd	     (%%rcx),%%ymm9 ,%%ymm9 \n\t"\
		"vmulpd	0x200(%%rdx),%%ymm10,%%ymm10\n\t"\
		"\n\t"\
		"vmulpd		%%ymm10,%%ymm5,%%ymm5	\n\t"\
		"vaddpd		%%ymm14,%%ymm5,%%ymm5	\n\t"\
		"vmovaps	%%ymm5,%%ymm10			\n\t"\
		"vroundpd	$0,%%ymm10,%%ymm10		\n\t"\
		"\n\t"\
		"movq	%[__sign_mask],%%rbx		\n\t"\
		"vsubpd		%%ymm10,%%ymm5,%%ymm5	\n\t"\
		"vandpd		(%%rbx),%%ymm5,%%ymm5	\n\t"\
		"vmaxpd		%%ymm15,%%ymm5,%%ymm15	\n\t"\
		/* cy   = DNINT(temp*baseinv[i1]): */\
		"addq	%%rdi,%%rsi		\n\t"\
		"vmovaps	%%ymm10	,%%ymm5			\n\t"\
		"vmulpd	0x600(%%rsi),%%ymm10,%%ymm10\n\t"\
		"vroundpd	$0,%%ymm10,%%ymm14		\n\t"\
		"vmovaps	%%ymm14,%%ymm10				\n\t"/* cy = cpy cy_out */\
		/* x = (temp-cy*base[i1])*wt: */\
		"vmulpd	0x400(%%rsi),%%ymm10,%%ymm10	\n\t"/* cy*base[i0123] */\
		"vsubpd		%%ymm10,%%ymm5 ,%%ymm5 	\n\t"\
		"vmulpd		%%ymm9 ,%%ymm5 ,%%ymm5 	\n\t"\
		/* Get ready for next set [RE3~] : by computing bjmod[0:3] += bw (mod n): */\
		"movq	%[__sse_n]	,%%rbx		\n\t"\
		"vmovaps		(%%rbx)	,%%xmm10	\n\t"\
		"movq	%[__sse_bw]	,%%rax		\n\t"\
		"vpaddd		(%%rax),%%xmm8 ,%%xmm8 	\n\t"\
		"vmovaps	%%xmm8 ,%%xmm9 		\n\t"\
		"vpcmpgtd	%%xmm10,%%xmm9 ,%%xmm9 	\n\t"\
		"vpand		%%xmm10,%%xmm9 ,%%xmm9 	\n\t"\
		"vpsubd		%%xmm9 ,%%xmm8 ,%%xmm8 	\n\t"\
		"\n\t"\
	/**********************************/\
	/* Do D.re-quartet: Data in ymm6: */\
	/**********************************/\
	"movslq		%[__p3],%%r15	\n\t"	\
	"prefetcht0	(%%r14,%%r15,8)	\n\t"\
		"movq	%[__sse_sw],%%rsi		\n\t"\
		"vmovaps	(%%rsi),	%%xmm11		\n\t"\
		"vpsubd		%%xmm8 ,%%xmm11,%%xmm11	\n\t"\
		"vmovmskps	%%xmm11,	%%rsi	\n\t"\
		"\n\t"\
		"movq	%[__n_minus_sil],%%rcx	\n\t"\
		"vmovd	0xC(%%rcx)	,%%xmm9 	\n\t"/* .d3 term of index quartet */\
		"vpshufd	$0,	%%xmm9 ,%%xmm9 		\n\t"\
		"vpsubd		%%xmm8 ,%%xmm9 ,%%xmm9 		\n\t"\
		"vmovmskps	%%xmm9 ,%%rcx		\n\t"\
		"\n\t"\
		"movq	%[__sinwt]	,%%rdx		\n\t"\
		"vmovd	0xC(%%rdx)	,%%xmm10	\n\t"/* .d3 term of index quartet */\
		"vpshufd	$0,	%%xmm10,%%xmm10		\n\t"\
		"vmovaps		%%xmm8 ,%%xmm11		\n\t"\
		"vpsubd		%%xmm10,%%xmm11,%%xmm11		\n\t"\
		"vmovmskps	%%xmm11,%%rdx		\n\t"\
		"\n\t"\
	"shlq	$5,%%rsi	\n\t"/* i0123 */\
	"shlq	$5,%%rcx	\n\t"/* m0123 */\
	"shlq	$5,%%rdx	\n\t"/* n0123 */\
		"\n\t"\
		"addq	%%rdi,%%rcx		\n\t"\
		"addq	%%rdi,%%rdx		\n\t"\
		"vbroadcastsd 0x860(%%rdi),%%ymm9 	\n\t"/* wtl */\
		"vbroadcastsd 0x868(%%rdi),%%ymm10	\n\t"/* wtn */\
		"vmulpd	%%ymm12,%%ymm9 ,%%ymm9 		\n\t"\
		"vmulpd	%%ymm13,%%ymm10,%%ymm10		\n\t"\
		"vmulpd	     (%%rcx),%%ymm9 ,%%ymm9 \n\t"\
		"vmulpd	0x200(%%rdx),%%ymm10,%%ymm10\n\t"\
		"\n\t"\
		"vmulpd		%%ymm10,%%ymm6,%%ymm6	\n\t"\
		"vaddpd		%%ymm14,%%ymm6,%%ymm6	\n\t"\
		"vmovaps	%%ymm6,%%ymm10			\n\t"\
		"vroundpd	$0,%%ymm10,%%ymm10		\n\t"\
		"\n\t"\
		"movq	%[__sign_mask],%%rbx		\n\t"\
		"vsubpd		%%ymm10,%%ymm6,%%ymm6	\n\t"\
		"vandpd		(%%rbx),%%ymm6,%%ymm6	\n\t"\
		"vmaxpd		%%ymm15,%%ymm6,%%ymm15	\n\t"\
		/* cy   = DNINT(temp*baseinv[i1]): */\
		"addq	%%rdi,%%rsi		\n\t"\
		"vmovaps	%%ymm10	,%%ymm6			\n\t"\
		"vmulpd	0x600(%%rsi),%%ymm10,%%ymm10\n\t"\
		"vroundpd	$0,%%ymm10,%%ymm14		\n\t"\
		"vmovaps	%%ymm14,%%ymm10				\n\t"/* cy = cpy cy_out */\
		/* x = (temp-cy*base[i1])*wt: */\
		"vmulpd	0x400(%%rsi),%%ymm10,%%ymm10	\n\t"/* cy*base[i0123] */\
		"vsubpd		%%ymm10,%%ymm6 ,%%ymm6 	\n\t"\
		"vmulpd		%%ymm9 ,%%ymm6 ,%%ymm6 	\n\t"\
		/* Get ready for next set [IM3~] : by computing bjmod[0:3] += bw (mod n): */\
		"movq	%[__sse_n]	,%%rbx		\n\t"\
		"vmovaps		(%%rbx)	,%%xmm10	\n\t"\
		"movq	%[__sse_bw]	,%%rax		\n\t"\
		"vpaddd		(%%rax),%%xmm8 ,%%xmm8 	\n\t"\
		"vmovaps	%%xmm8 ,%%xmm9 		\n\t"\
		"vpcmpgtd	%%xmm10,%%xmm9 ,%%xmm9 	\n\t"\
		"vpand		%%xmm10,%%xmm9 ,%%xmm9 	\n\t"\
		"vpsubd		%%xmm9 ,%%xmm8 ,%%xmm8 	\n\t"\
		"\n\t"\
	/**********************************/\
	/* Do D.im-quartet: Data in ymm7: */\
	/**********************************/\
		"movq	%[__sse_sw],%%rsi		\n\t"\
		"vmovaps	(%%rsi),	%%xmm11		\n\t"\
		"vpsubd		%%xmm8 ,%%xmm11,%%xmm11	\n\t"\
		"vmovmskps	%%xmm11,	%%rsi	\n\t"\
		"\n\t"\
		"movq	%[__n_minus_silp1],%%rcx	\n\t"\
		"vmovd	0xC(%%rcx)	,%%xmm9 	\n\t"/* .d3 term of index quartet */\
		"vpshufd	$0,	%%xmm9 ,%%xmm9 		\n\t"\
		"vpsubd		%%xmm8 ,%%xmm9 ,%%xmm9 		\n\t"\
		"vmovmskps	%%xmm9 ,%%rcx		\n\t"\
		"\n\t"\
		"movq	%[__sinwtm1]	,%%rdx		\n\t"\
		"vmovd	0xC(%%rdx)	,%%xmm10	\n\t"/* .d3 term of index quartet */\
		"vpshufd	$0,	%%xmm10,%%xmm10		\n\t"\
		"vmovaps		%%xmm8 ,%%xmm11		\n\t"\
		"vpsubd		%%xmm10,%%xmm11,%%xmm11		\n\t"\
		"vmovmskps	%%xmm11,%%rdx		\n\t"\
		"\n\t"\
	"shlq	$5,%%rsi	\n\t"/* i0123 */\
	"shlq	$5,%%rcx	\n\t"/* m0123 */\
	"shlq	$5,%%rdx	\n\t"/* n0123 */\
		"\n\t"\
		"addq	%%rdi,%%rcx		\n\t"\
		"addq	%%rdi,%%rdx		\n\t"\
		"vbroadcastsd 0x870(%%rdi),%%ymm9 	\n\t"/* wtlp1 */\
		"vbroadcastsd 0x878(%%rdi),%%ymm10	\n\t"/* wtnm1 */\
		"vmulpd	%%ymm12,%%ymm9 ,%%ymm9 		\n\t"\
		"vmulpd	%%ymm13,%%ymm10,%%ymm10		\n\t"\
		"vmulpd	     (%%rcx),%%ymm9 ,%%ymm9 \n\t"\
		"vmulpd	0x200(%%rdx),%%ymm10,%%ymm10\n\t"\
		"\n\t"\
		"vmulpd		%%ymm10,%%ymm7,%%ymm7	\n\t"\
		"vaddpd		%%ymm14,%%ymm7,%%ymm7	\n\t"\
		"vmovaps	%%ymm7,%%ymm10			\n\t"\
		"vroundpd	$0,%%ymm10,%%ymm10		\n\t"\
		"\n\t"\
		"movq	%[__sign_mask],%%rbx		\n\t"\
		"vsubpd		%%ymm10,%%ymm7,%%ymm7	\n\t"\
		"vandpd		(%%rbx),%%ymm7,%%ymm7	\n\t"\
		"vmaxpd		%%ymm15,%%ymm7,%%ymm15	\n\t"\
		/* cy   = DNINT(temp*baseinv[i1]): */\
		"addq	%%rdi,%%rsi		\n\t"\
		"vmovaps	%%ymm10	,%%ymm7			\n\t"\
		"vmulpd	0x600(%%rsi),%%ymm10,%%ymm10\n\t"\
		"vroundpd	$0,%%ymm10,%%ymm14		\n\t"\
		"vmovaps	%%ymm14,%%ymm10				\n\t"/* cy = cpy cy_out */\
		/* x = (temp-cy*base[i1])*wt: */\
		"vmulpd	0x400(%%rsi),%%ymm10,%%ymm10	\n\t"/* cy*base[i0123] */\
		"vsubpd		%%ymm10,%%ymm7 ,%%ymm7 	\n\t"\
		"vmulpd		%%ymm9 ,%%ymm7 ,%%ymm7 	\n\t"\
		/* Update wts-array pointers in preparation for next call of the macro: */\
		"movq	%[__wtA]	,%%rax		\n\t"\
		"movq	%[__wtB]	,%%rbx		\n\t"\
		"movq	%[__wtC]	,%%rcx		\n\t"\
		"addq	$0x20	,%%rax			\n\t"/* add0 += 4 */\
		"subq	$0x20	,%%rbx			\n\t"/* add1 -= 4 */\
		"subq	$0x20	,%%rcx			\n\t"/* add2 -= 4 */\
		"movq	%%rax	,%[__wtA]		\n\t"\
		"movq	%%rbx	,%[__wtB]		\n\t"\
		"movq	%%rcx	,%[__wtC]		\n\t"\
		"\n\t"\
		/* Get ready for store of final-updated bjmod[0:3] values: */\
		"movq	%[__sse_n]	,%%rbx		\n\t"\
		"vmovaps		(%%rbx)	,%%xmm10	\n\t"\
		"movq	%[__sse_bw]	,%%rax		\n\t"\
		"vpaddd		(%%rax),%%xmm8 ,%%xmm8 	\n\t"\
		"vmovaps	%%xmm8 ,%%xmm9 		\n\t"\
		"vpcmpgtd	%%xmm10,%%xmm9 ,%%xmm9 	\n\t"\
		"vpand		%%xmm10,%%xmm9 ,%%xmm9 	\n\t"\
		"vpsubd		%%xmm9 ,%%xmm8 ,%%xmm8 	\n\t"\
		"\n\t"\
		/* Store bjmodn index quartet: */\
		"movq	%[__bjmod_0],%%rcx		\n\t"\
		"vmovaps	%%xmm8,(%%rcx)			\n\t"\
		/* Store cy_out: */\
		"movq		%[__cy],%%rbx	\n\t"\
		"vmovaps	%%ymm14,(%%rbx)	\n\t"/* ymm14 = cy_in */\
		/* Store maxerr: */\
		"movq		%[__half_arr],%%rdx		\n\t"\
		"vmovaps	%%ymm15,-0x40(%%rdx)	\n\t"\
	/* 4-way transpose of outputs (Re, Im parts separately): Inputs from ymm0-7. Outputs into r0/1,2/3,4/5.6/7: */	\
	/* Because default inputs for our 4 x 4 transpose macro (e.g. the one used at start of this carry macro) */\
	/* are into ymm4/2/8/2, munge inputs into that order, resolving name-conflicts via use of the now-available ymm8-15 for outputs: */\
		"movq		%[__data],%%rax			\n\t"\
		"vshufpd	$15,%%ymm2,%%ymm0,%%ymm10					\n\t		vshufpd	$15,%%ymm3,%%ymm1,%%ymm11						\n\t"\
		"vshufpd	$0 ,%%ymm2,%%ymm0,%%ymm0					\n\t		vshufpd	$0 ,%%ymm3,%%ymm1,%%ymm1						\n\t"\
		"vshufpd	$15,%%ymm6,%%ymm4,%%ymm12					\n\t		vshufpd	$15,%%ymm7,%%ymm5,%%ymm13						\n\t"\
		"vshufpd	$0 ,%%ymm6,%%ymm4,%%ymm4					\n\t		vshufpd	$0 ,%%ymm7,%%ymm5,%%ymm5						\n\t"\
		"vperm2f128 $32,%%ymm12,%%ymm10,%%ymm2 		/* Re B	*/	\n\t		vperm2f128 $32,%%ymm13,%%ymm11,%%ymm3		/* Im B	*/	\n\t"\
		"vperm2f128 $49,%%ymm12,%%ymm10,%%ymm10		/* Re D	*/	\n\t		vperm2f128 $49,%%ymm13,%%ymm11,%%ymm11		/* Im D	*/	\n\t"\
		"vperm2f128 $32,%%ymm4 ,%%ymm0 ,%%ymm12		/* Re A	*/	\n\t		vperm2f128 $32,%%ymm5 ,%%ymm1 ,%%ymm13 		/* Im A	*/	\n\t"\
		"vperm2f128 $49,%%ymm4 ,%%ymm0 ,%%ymm0 		/* Re C	*/	\n\t		vperm2f128 $49,%%ymm5 ,%%ymm1 ,%%ymm1		/* Im C	*/	\n\t"\
		/* And write 'em back to memory: */\
		"vmovaps	%%ymm12,     (%%rax)						\n\t		vmovaps	%%ymm13,0x020(%%rax)				\n\t"\
		"vmovaps	%%ymm2 ,0x040(%%rax)						\n\t		vmovaps	%%ymm3 ,0x060(%%rax)				\n\t"\
		"vmovaps	%%ymm0 ,0x080(%%rax)						\n\t		vmovaps	%%ymm1 ,0x0a0(%%rax)				\n\t"\
		"vmovaps	%%ymm10,0x0c0(%%rax)						\n\t		vmovaps	%%ymm11,0x0e0(%%rax)				\n\t"\
		:					/* outputs: none */\
		: [__data]		"m" (Xdata)	/* All inputs from memory addresses here */\
		, [__wtA]		"m" (XwtA)		\
		, [__wtB]		"m" (XwtB)		\
		, [__wtC]		"m" (XwtC)		\
		, [__cy]		"m" (Xcy)		\
		, [__bjmod_0]	"m" (Xbjmod_0)		\
		, [__half_arr]	"m" (Xhalf_arr)		\
		, [__i]			"m" (Xi)			\
		, [__n_minus_silp1] "m" (Xn_minus_silp1)\
		, [__n_minus_sil]	"m" (Xn_minus_sil)	\
		, [__sign_mask]	"m" (Xsign_mask)	\
		, [__sinwt]		"m" (Xsinwt)		\
		, [__sinwtm1]	"m" (Xsinwtm1)		\
		, [__sse_bw]	"m" (Xsse_bw)		\
		, [__sse_n]		"m" (Xsse_n)		\
		, [__sse_sw]	"m" (Xsse_sw)		\
		/* Prefetch: base address and 3 index offsets */\
		,	[__add0] "m" (Xadd0)\
		,	[__p1] "m" (Xp1)\
		,	[__p2] "m" (Xp2)\
		,	[__p3] "m" (Xp3)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r14","r15","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
	);\
	}

	#define AVX_cmplx_carry_norm_errcheck1_X4(Xdata,XwtA,XwtB,XwtC,Xcy,Xbjmod_0,Xhalf_arr,   Xn_minus_silp1,Xn_minus_sil,Xsign_mask,Xsinwt,Xsinwtm1,Xsse_bw,Xsse_n,Xsse_sw, Xadd0,Xp1,Xp2,Xp3)\
	{\
	__asm__ volatile (\
	"movq	%[__add0],%%r14	\n\t"/* base address for 4 prefetches-from-main-data-array spread through this macro */\
		"movq		%[__data],%%rax		\n\t"\
	/* 4-way transpose of inputs (Re, Im parts separately): Inputs from r0/1,2/3,4/5.6/7. Outputs into ymm0-7: */\
	/* Real parts use ymm0,2,4,6, ymm8 as tmp-reg:					Imag parts use ymm1,3,5,7, ymm9 as tm-reg: */\
		"vmovaps	     (%%rax),%%ymm4						\n\t		vmovaps	0x020(%%rax),%%ymm5							\n\t"\
		"vmovaps	0x040(%%rax),%%ymm2						\n\t		vmovaps	0x060(%%rax),%%ymm3							\n\t"\
		"vshufpd	$15,%%ymm2,%%ymm4,%%ymm6				\n\t		vshufpd	$15,%%ymm3,%%ymm5,%%ymm7					\n\t"\
		"vshufpd	$0 ,%%ymm2,%%ymm4,%%ymm4				\n\t		vshufpd	$0 ,%%ymm3,%%ymm5,%%ymm5					\n\t"\
		"vmovaps	0x080(%%rax),%%ymm8						\n\t		vmovaps	0x0a0(%%rax),%%ymm9							\n\t"\
		"vmovaps	0x0c0(%%rax),%%ymm2						\n\t		vmovaps	0x0e0(%%rax),%%ymm3							\n\t"\
		"vshufpd	$15,%%ymm2,%%ymm8,%%ymm0				\n\t		vshufpd	$15,%%ymm3,%%ymm9,%%ymm1					\n\t"\
		"vshufpd	$0 ,%%ymm2,%%ymm8,%%ymm8				\n\t		vshufpd	$0 ,%%ymm3,%%ymm9,%%ymm9					\n\t"\
		"vperm2f128 $32,%%ymm0,%%ymm6,%%ymm2	/* Re B	*/	\n\t		vperm2f128 $32,%%ymm1,%%ymm7,%%ymm3		/* Im B	*/	\n\t"\
		"vperm2f128 $49,%%ymm0,%%ymm6,%%ymm6	/* Re D	*/	\n\t		vperm2f128 $49,%%ymm1,%%ymm7,%%ymm7		/* Im D	*/	\n\t"\
		"vperm2f128 $32,%%ymm8,%%ymm4,%%ymm0	/* Re A	*/	\n\t		vperm2f128 $32,%%ymm9,%%ymm5,%%ymm1 	/* Im A	*/	\n\t"\
		"vperm2f128 $49,%%ymm8,%%ymm4,%%ymm4	/* Re C	*/	\n\t		vperm2f128 $49,%%ymm9,%%ymm5,%%ymm5 	/* Im C	*/	\n\t"\
		/* Won't need main-array again until output transpose, so re-use rax for half_arr */\
		"movq	 %[__half_arr],%%rax	\n\t"\
		/* half_arr + 16*[0,1,2,3] = [wt,wt_inv,base,baseinv] */\
		"movq		%[__cy],%%rbx	\n\t"\
		"vmovaps	(%%rbx),%%ymm14	\n\t"/* ymm14 = cy_in */\
		"vmovaps	-0x40(%%rax),%%ymm15	\n\t"/* ymm15 = maxerr */\
	/**********************************/\
	/* Do A.re-quartet: Data in ymm0: */\
	/**********************************/\
	"prefetcht0	(%%r14)		\n\t"\
		"movq	%[__bjmod_0],%%rsi		\n\t"\
		"vmovaps	(%%rsi),	%%xmm8 		\n\t"/* bjmod[0:3]. PERSISTENT COPY OF BJMOD[0:3] REMAINS IN XMM8. */\
		"movq	%[__sse_sw],%%rsi		\n\t"\
		"vmovaps	(%%rsi),	%%xmm11		\n\t"/* sw[0:3] */\
		"vpsubd		%%xmm8 ,%%xmm11,%%xmm11	\n\t"/* sw[0:3] - bjmod[0:3] */\
		"vmovmskps	%%xmm11,	%%rsi	\n\t"/* Extract sign bits into 4-bit signmask <i3|i2|i1|i0>; idxs into base/inv table */\
		"\n\t"\
		"movq	%[__n_minus_sil],%%rcx	\n\t"\
		"vmovd	0x0(%%rcx)	,%%xmm9 	\n\t"/* n_minus_sil in low 32 bits of xmm9  */\
		"vpshufd	$0,	%%xmm9 ,%%xmm9 		\n\t"/* Broadcast low 32 bits of xmm9  to all 4 slots of xmm9  */\
		"vpsubd		%%xmm8 ,%%xmm9 ,%%xmm9 		\n\t"/* n_minus_sil - bjmod[0:3] */\
		"vmovmskps	%%xmm9 ,%%rcx		\n\t"/* Extract sign bits into 4-bit signmask <m3|m2|m1|m0>; idxs into base/inv tables -> byte[2] of ecx... */\
		"\n\t"\
		"movq	%[__sinwt]	,%%rdx		\n\t"\
		"vmovd	0x0(%%rdx)	,%%xmm10		\n\t"/* sinwt in low 32 bits of xmm10 */\
		"vpshufd	$0,	%%xmm10,%%xmm10		\n\t"/* Broadcast low 32 bits of xmm10 to all 4 slots of xmm10 */\
		"vmovaps		%%xmm8 ,%%xmm11		\n\t"/* xmm11 = bjmod[0:3] copy */\
		"vpsubd		%%xmm10,%%xmm11,%%xmm11		\n\t"/* bjmod[0:3] - sinwt */\
		"vmovmskps	%%xmm11,%%rdx		\n\t"/* Extract sign bits into 4-bit signmask <n3|n2|n1|n0>; idxs into base/inv tables -> byte[1] of edx... */\
		"\n\t"\
	"shlq	$5,%%rsi	\n\t"/* i0123 << 5 (= lshift to give ptr offset for ymm-size data */\
	"shlq	$5,%%rcx	\n\t"/* m0123 << 5 (= lshift to give ptr offset for ymm-size data */\
	"shlq	$5,%%rdx	\n\t"/* n0123 << 5 (= lshift to give ptr offset for ymm-size data */\
		"\n\t"\
		"movq	%[__half_arr],%%rdi	\n\t"\
		"movq	%[__wtA]	,%%rax		\n\t"\
		"movq	%[__wtB]	,%%rbx		\n\t"\
		"vmovaps	     (%%rax),%%ymm12 	\n\t"/* wtA[j  ]; ebx FREE */\
		"vmovups	-0x10(%%rbx),%%ymm13	\n\t"/* wtB[j-1]; load doubles from rcx+[-0x10,-0x08, 0, +0x08] - note this address is ymm (16-byte) but not ymm (32-byte) aligned */\
		/* reverse-running indexing used for inv-wts really means we need to reverse ordering of 4 doubles d0-3 in ymm13 */\
		"vshufpd	$5,%%ymm13,%%ymm13,%%ymm13	\n\t"/* d[0123] -> d[1032] */\
		"vperm2f128 $1,%%ymm13,%%ymm13,%%ymm13	\n\t"/* d[1032] -> d[3210] */\
		"\n\t"\
		"addq	%%rdi,%%rcx		\n\t"\
		"addq	%%rdi,%%rdx		\n\t"\
		/* SSE2 version has double-copies in wtl/wtn ... AVX replaces redundant-data loads with load-with-broadcast: */\
		"vbroadcastsd 0x800(%%rdi),%%ymm9 	\n\t"/* wtl */\
		"vbroadcastsd 0x808(%%rdi),%%ymm10	\n\t"/* wtn */\
		"vmulpd	%%ymm12,%%ymm9 ,%%ymm9 		\n\t"/* wt   =wtA*wtl */\
		"vmulpd	%%ymm13,%%ymm10,%%ymm10		\n\t"/* wtinv=wtB*wtn */\
		"vmulpd	     (%%rcx),%%ymm9 ,%%ymm9 	\n\t"/* wt    *= one_half[m0123] */\
		"vmulpd	0x200(%%rdx),%%ymm10,%%ymm10	\n\t"/* wtinv *= one_half[16+n0123] */\
		"\n\t"\
		"vmulpd		%%ymm10,%%ymm0,%%ymm0		\n\t"/* x *= wtinv; ymm10 FREE */\
		"vaddpd		%%ymm14,%%ymm0,%%ymm0		\n\t"/* x *= wtinv + cy */\
		"vmovaps	%%ymm0,%%ymm10		\n\t"/* temp = x */\
		"vroundpd	$0,%%ymm10,%%ymm10	\n\t"/* temp = DNINT(x) */\
		"\n\t"\
		"movq	%[__sign_mask],%%rbx	\n\t"\
		"vsubpd		%%ymm10,%%ymm0,%%ymm0	\n\t"/* x - temp */\
		"vandpd		(%%rbx),%%ymm0,%%ymm0	\n\t"/* frac = fabs(x-temp) */\
		"vmaxpd		%%ymm15,%%ymm0,%%ymm15	\n\t"/* if(frac > maxerr) maxerr=frac */\
		/* cy   = DNINT(temp*baseinv[i1]): */\
		"addq	%%rdi,%%rsi		\n\t"\
		"vmovaps	%%ymm10	,%%ymm0				\n\t"/* cpy temp */\
		"vmulpd	0x600(%%rsi),%%ymm10,%%ymm10	\n\t"/* temp*baseinv[i0123] */\
		"vroundpd	$0,%%ymm10,%%ymm14			\n\t"/* cy_out */\
		"vmovaps	%%ymm14,%%ymm10				\n\t"/* cy = cpy cy_out */\
		/* x = (temp-cy*base[i1])*wt: */\
		"vmulpd	0x400(%%rsi),%%ymm10,%%ymm10	\n\t"/* cy*base[i0123] */\
		"vsubpd		%%ymm10,%%ymm0 ,%%ymm0 		\n\t"/* x = (temp-cy*base[i0123]) */\
		"vmulpd		%%ymm9 ,%%ymm0 ,%%ymm0 		\n\t"/* x *= wt */\
		/* Get ready for next set [IM0~] by computing bjmod[0:3] += bw (mod n): */\
		"movq	%[__sse_n]	,%%rbx		\n\t"\
		"vmovaps		(%%rbx)	,%%xmm10	\n\t"\
		"movq	%[__sse_bw]	,%%rax		\n\t"\
		"vpaddd		(%%rax),%%xmm8 ,%%xmm8 	\n\t"\
		"vmovaps	%%xmm8 ,%%xmm9 		\n\t"\
		"vpcmpgtd	%%xmm10,%%xmm9 ,%%xmm9 	\n\t"\
		"vpand		%%xmm10,%%xmm9 ,%%xmm9 	\n\t"\
		"vpsubd		%%xmm9 ,%%xmm8 ,%%xmm8 	\n\t"\
		"\n\t"\
	/**********************************/\
	/* Do A.im-quartet: Data in ymm1: */\
	/**********************************/\
		"movq	%[__sse_sw],%%rsi		\n\t"\
		"vmovaps	(%%rsi),	%%xmm11		\n\t"/* sw[0:3] */\
		"vpsubd		%%xmm8 ,%%xmm11,%%xmm11	\n\t"/* sw[0:3] - bjmod[0:3] */\
		"vmovmskps	%%xmm11,	%%rsi	\n\t"/* Extract sign bits into 4-bit signmask <i3|i2|i1|i0>; idxs into base/inv table */\
		"\n\t"\
		"movq	%[__n_minus_silp1],%%rcx	\n\t"\
		"vmovd	0x0(%%rcx)	,%%xmm9 	\n\t"/* n_minus_sil in low 32 bits of xmm9  */\
		"vpshufd	$0,	%%xmm9 ,%%xmm9 		\n\t"/* Broadcast low 32 bits of xmm9  to all 4 slots of xmm9  */\
		"vpsubd		%%xmm8 ,%%xmm9 ,%%xmm9 		\n\t"/* n_minus_sil - bjmod[0:3] */\
		"vmovmskps	%%xmm9 ,%%rcx		\n\t"/* Extract sign bits into 4-bit signmask <m3|m2|m1|m0>; idxs into base/inv tables -> byte[2] of ecx... */\
		"\n\t"\
		"movq	%[__sinwtm1]	,%%rdx		\n\t"\
		"vmovd	0x0(%%rdx)	,%%xmm10		\n\t"/* sinwt in low 32 bits of xmm10 */\
		"vpshufd	$0,	%%xmm10,%%xmm10		\n\t"/* Broadcast low 32 bits of xmm10 to all 4 slots of xmm10 */\
		"vmovaps		%%xmm8 ,%%xmm11		\n\t"/* xmm11 = bjmod[0:3] copy */\
		"vpsubd		%%xmm10,%%xmm11,%%xmm11		\n\t"/* bjmod[0:3] - sinwt */\
		"vmovmskps	%%xmm11,%%rdx		\n\t"/* Extract sign bits into 4-bit signmask <n3|n2|n1|n0>; idxs into base/inv tables -> byte[1] of edx... */\
		"\n\t"\
	"shlq	$5,%%rsi	\n\t"/* i0123 << 5 (= lshift to give ptr offset for ymm-size data */\
	"shlq	$5,%%rcx	\n\t"/* m0123 << 5 (= lshift to give ptr offset for ymm-size data */\
	"shlq	$5,%%rdx	\n\t"/* n0123 << 5 (= lshift to give ptr offset for ymm-size data */\
		"\n\t"\
		"movq	%[__wtC]	,%%rbx		\n\t"/* wtA unchanged; wtB == wtC for remaining 7 of 8 sets of carries */\
		"vmovups	-0x10(%%rbx),%%ymm13	\n\t"/* wtC[j-1]; load doubles from rcx+[-0x10,-0x08, 0, +0x08] - note this address is ymm (16-byte) but not ymm (32-byte) aligned */\
		"vshufpd	$5,%%ymm13,%%ymm13,%%ymm13	\n\t"/* d[0123] -> d[1032] */\
		"vperm2f128 $1,%%ymm13,%%ymm13,%%ymm13	\n\t"/* d[1032] -> d[3210] */\
		"\n\t"\
		"addq	%%rdi,%%rcx		\n\t"\
		"addq	%%rdi,%%rdx		\n\t"\
		/* SSE2 version has double-copies in wtl/wtn ... AVX replaces redundant-data loads with load-with-broadcast: */\
		"vbroadcastsd 0x810(%%rdi),%%ymm9 	\n\t"/* wtlp1 */\
		"vbroadcastsd 0x818(%%rdi),%%ymm10	\n\t"/* wtnm1 */\
		"vmulpd	%%ymm12,%%ymm9 ,%%ymm9 		\n\t"/* wt   =wtA*wtlp1 */\
		"vmulpd	%%ymm13,%%ymm10,%%ymm10		\n\t"/* wtinv=wtC*wtnm1 */\
		"vmulpd	     (%%rcx),%%ymm9 ,%%ymm9 	\n\t"/* wt    *= one_half[m0123] */\
		"vmulpd	0x200(%%rdx),%%ymm10,%%ymm10	\n\t"/* wtinv *= one_half[16+n0123] */\
		"\n\t"\
		"vmulpd		%%ymm10,%%ymm1,%%ymm1		\n\t"/* x *= wtinv; ymm10 FREE */\
		"vaddpd		%%ymm14,%%ymm1,%%ymm1		\n\t"/* x *= wtinv + cy */\
		"vmovaps	%%ymm1,%%ymm10		\n\t"/* temp = x */\
		"vroundpd	$0,%%ymm10,%%ymm10	\n\t"/* temp = DNINT(x) */\
		"\n\t"\
		"movq	%[__sign_mask],%%rbx	\n\t"\
		"vsubpd		%%ymm10,%%ymm1,%%ymm1	\n\t"/* x - temp */\
		"vandpd		(%%rbx),%%ymm1,%%ymm1	\n\t"/* frac = fabs(x-temp) */\
		"vmaxpd		%%ymm15,%%ymm1,%%ymm15	\n\t"/* if(frac > maxerr) maxerr=frac */\
		/* cy   = DNINT(temp*baseinv[i1]): */\
		"addq	%%rdi,%%rsi		\n\t"\
		"vmovaps	%%ymm10	,%%ymm1				\n\t"/* cpy temp */\
		"vmulpd	0x600(%%rsi),%%ymm10,%%ymm10	\n\t"/* temp*baseinv[i0123] */\
		"vroundpd	$0,%%ymm10,%%ymm14			\n\t"/* cy_out */\
		"vmovaps	%%ymm14,%%ymm10				\n\t"/* cy = cpy cy_out */\
		/* x = (temp-cy*base[i1])*wt: */\
		"vmulpd	0x400(%%rsi),%%ymm10,%%ymm10	\n\t"/* cy*base[i0123] */\
		"vsubpd		%%ymm10,%%ymm1 ,%%ymm1 		\n\t"/* x = (temp-cy*base[i0123]) */\
		"vmulpd		%%ymm9 ,%%ymm1 ,%%ymm1 		\n\t"/* x *= wt */\
		/* Get ready for next set [RE1~] by computing bjmod[0:3] += bw (mod n): */\
		"movq	%[__sse_n]	,%%rbx		\n\t"\
		"vmovaps		(%%rbx)	,%%xmm10	\n\t"\
		"movq	%[__sse_bw]	,%%rax		\n\t"\
		"vpaddd		(%%rax),%%xmm8 ,%%xmm8 	\n\t"\
		"vmovaps	%%xmm8 ,%%xmm9 		\n\t"\
		"vpcmpgtd	%%xmm10,%%xmm9 ,%%xmm9 	\n\t"\
		"vpand		%%xmm10,%%xmm9 ,%%xmm9 	\n\t"\
		"vpsubd		%%xmm9 ,%%xmm8 ,%%xmm8 	\n\t"\
		"\n\t"\
	/**********************************/\
	/* Do B.re-quartet: Data in ymm2: */\
	/**********************************/\
	"movslq		%[__p1],%%r15	\n\t"	\
	"prefetcht0	(%%r14,%%r15,8)	\n\t"\
		"movq	%[__sse_sw],%%rsi		\n\t"\
		"vmovaps	(%%rsi),	%%xmm11		\n\t"\
		"vpsubd		%%xmm8 ,%%xmm11,%%xmm11	\n\t"\
		"vmovmskps	%%xmm11,	%%rsi	\n\t"\
		"\n\t"\
		"movq	%[__n_minus_sil],%%rcx	\n\t"\
		"vmovd	0x4(%%rcx)	,%%xmm9 	\n\t"/* .d1 term of index quartet */\
		"vpshufd	$0,	%%xmm9 ,%%xmm9 		\n\t"\
		"vpsubd		%%xmm8 ,%%xmm9 ,%%xmm9 		\n\t"\
		"vmovmskps	%%xmm9 ,%%rcx		\n\t"\
		"\n\t"\
		"movq	%[__sinwt]	,%%rdx		\n\t"\
		"vmovd	0x4(%%rdx)	,%%xmm10	\n\t"/* .d1 term of index quartet */\
		"vpshufd	$0,	%%xmm10,%%xmm10		\n\t"\
		"vmovaps		%%xmm8 ,%%xmm11		\n\t"\
		"vpsubd		%%xmm10,%%xmm11,%%xmm11		\n\t"\
		"vmovmskps	%%xmm11,%%rdx		\n\t"\
		"\n\t"\
	"shlq	$5,%%rsi	\n\t"/* i0123 */\
	"shlq	$5,%%rcx	\n\t"/* m0123 */\
	"shlq	$5,%%rdx	\n\t"/* n0123 */\
		"\n\t"\
		"addq	%%rdi,%%rcx		\n\t"\
		"addq	%%rdi,%%rdx		\n\t"\
		"vbroadcastsd 0x820(%%rdi),%%ymm9 	\n\t"/* wtl */\
		"vbroadcastsd 0x828(%%rdi),%%ymm10	\n\t"/* wtn */\
		"vmulpd	%%ymm12,%%ymm9 ,%%ymm9 		\n\t"\
		"vmulpd	%%ymm13,%%ymm10,%%ymm10		\n\t"\
		"vmulpd	     (%%rcx),%%ymm9 ,%%ymm9 \n\t"\
		"vmulpd	0x200(%%rdx),%%ymm10,%%ymm10\n\t"\
		"\n\t"\
		"vmulpd		%%ymm10,%%ymm2,%%ymm2	\n\t"\
		"vaddpd		%%ymm14,%%ymm2,%%ymm2	\n\t"\
		"vmovaps	%%ymm2,%%ymm10			\n\t"\
		"vroundpd	$0,%%ymm10,%%ymm10		\n\t"\
		"\n\t"\
		"movq	%[__sign_mask],%%rbx		\n\t"\
		"vsubpd		%%ymm10,%%ymm2,%%ymm2	\n\t"\
		"vandpd		(%%rbx),%%ymm2,%%ymm2	\n\t"\
		"vmaxpd		%%ymm15,%%ymm2,%%ymm15	\n\t"\
		/* cy   = DNINT(temp*baseinv[i1]): */\
		"addq	%%rdi,%%rsi		\n\t"\
		"vmovaps	%%ymm10	,%%ymm2			\n\t"\
		"vmulpd	0x600(%%rsi),%%ymm10,%%ymm10\n\t"\
		"vroundpd	$0,%%ymm10,%%ymm14		\n\t"\
		"vmovaps	%%ymm14,%%ymm10				\n\t"/* cy = cpy cy_out */\
		/* x = (temp-cy*base[i1])*wt: */\
		"vmulpd	0x400(%%rsi),%%ymm10,%%ymm10	\n\t"/* cy*base[i0123] */\
		"vsubpd		%%ymm10,%%ymm2 ,%%ymm2 	\n\t"\
		"vmulpd		%%ymm9 ,%%ymm2 ,%%ymm2 	\n\t"\
		/* Get ready for next set [IM1~] by computing bjmod[0:3] += bw (mod n): */\
		"movq	%[__sse_n]	,%%rbx		\n\t"\
		"vmovaps		(%%rbx)	,%%xmm10	\n\t"\
		"movq	%[__sse_bw]	,%%rax		\n\t"\
		"vpaddd		(%%rax),%%xmm8 ,%%xmm8 	\n\t"\
		"vmovaps	%%xmm8 ,%%xmm9 		\n\t"\
		"vpcmpgtd	%%xmm10,%%xmm9 ,%%xmm9 	\n\t"\
		"vpand		%%xmm10,%%xmm9 ,%%xmm9 	\n\t"\
		"vpsubd		%%xmm9 ,%%xmm8 ,%%xmm8 	\n\t"\
		"\n\t"\
	/**********************************/\
	/* Do B.im-quartet: Data in ymm3: */\
	/**********************************/\
		"movq	%[__sse_sw],%%rsi		\n\t"\
		"vmovaps	(%%rsi),	%%xmm11		\n\t"\
		"vpsubd		%%xmm8 ,%%xmm11,%%xmm11	\n\t"\
		"vmovmskps	%%xmm11,	%%rsi	\n\t"\
		"\n\t"\
		"movq	%[__n_minus_silp1],%%rcx	\n\t"\
		"vmovd	0x4(%%rcx)	,%%xmm9 	\n\t"/* .d1 term of index quartet */\
		"vpshufd	$0,	%%xmm9 ,%%xmm9 		\n\t"\
		"vpsubd		%%xmm8 ,%%xmm9 ,%%xmm9 		\n\t"\
		"vmovmskps	%%xmm9 ,%%rcx		\n\t"\
		"\n\t"\
		"movq	%[__sinwtm1]	,%%rdx		\n\t"\
		"vmovd	0x4(%%rdx)	,%%xmm10	\n\t"/* .d1 term of index quartet */\
		"vpshufd	$0,	%%xmm10,%%xmm10		\n\t"\
		"vmovaps		%%xmm8 ,%%xmm11		\n\t"\
		"vpsubd		%%xmm10,%%xmm11,%%xmm11		\n\t"\
		"vmovmskps	%%xmm11,%%rdx		\n\t"\
		"\n\t"\
	"shlq	$5,%%rsi	\n\t"/* i0123 */\
	"shlq	$5,%%rcx	\n\t"/* m0123 */\
	"shlq	$5,%%rdx	\n\t"/* n0123 */\
		"\n\t"\
		"addq	%%rdi,%%rcx		\n\t"\
		"addq	%%rdi,%%rdx		\n\t"\
		"vbroadcastsd 0x830(%%rdi),%%ymm9 	\n\t"/* wtlp1 */\
		"vbroadcastsd 0x838(%%rdi),%%ymm10	\n\t"/* wtnm1 */\
		"vmulpd	%%ymm12,%%ymm9 ,%%ymm9 		\n\t"\
		"vmulpd	%%ymm13,%%ymm10,%%ymm10		\n\t"\
		"vmulpd	     (%%rcx),%%ymm9 ,%%ymm9 \n\t"\
		"vmulpd	0x200(%%rdx),%%ymm10,%%ymm10\n\t"\
		"\n\t"\
		"vmulpd		%%ymm10,%%ymm3,%%ymm3	\n\t"\
		"vaddpd		%%ymm14,%%ymm3,%%ymm3	\n\t"\
		"vmovaps	%%ymm3,%%ymm10			\n\t"\
		"vroundpd	$0,%%ymm10,%%ymm10		\n\t"\
		"\n\t"\
		"movq	%[__sign_mask],%%rbx		\n\t"\
		"vsubpd		%%ymm10,%%ymm3,%%ymm3	\n\t"\
		"vandpd		(%%rbx),%%ymm3,%%ymm3	\n\t"\
		"vmaxpd		%%ymm15,%%ymm3,%%ymm15	\n\t"\
		/* cy   = DNINT(temp*baseinv[i1]): */\
		"addq	%%rdi,%%rsi		\n\t"\
		"vmovaps	%%ymm10	,%%ymm3			\n\t"\
		"vmulpd	0x600(%%rsi),%%ymm10,%%ymm10\n\t"\
		"vroundpd	$0,%%ymm10,%%ymm14		\n\t"\
		"vmovaps	%%ymm14,%%ymm10				\n\t"/* cy = cpy cy_out */\
		/* x = (temp-cy*base[i1])*wt: */\
		"vmulpd	0x400(%%rsi),%%ymm10,%%ymm10	\n\t"/* cy*base[i0123] */\
		"vsubpd		%%ymm10,%%ymm3 ,%%ymm3 	\n\t"\
		"vmulpd		%%ymm9 ,%%ymm3 ,%%ymm3 	\n\t"\
		/* Get ready for next set [RE2~] by computing bjmod[0:3] += bw (mod n): */\
		"movq	%[__sse_n]	,%%rbx		\n\t"\
		"vmovaps		(%%rbx)	,%%xmm10	\n\t"\
		"movq	%[__sse_bw]	,%%rax		\n\t"\
		"vpaddd		(%%rax),%%xmm8 ,%%xmm8 	\n\t"\
		"vmovaps	%%xmm8 ,%%xmm9 		\n\t"\
		"vpcmpgtd	%%xmm10,%%xmm9 ,%%xmm9 	\n\t"\
		"vpand		%%xmm10,%%xmm9 ,%%xmm9 	\n\t"\
		"vpsubd		%%xmm9 ,%%xmm8 ,%%xmm8 	\n\t"\
		"\n\t"\
	/**********************************/\
	/* Do C.re-quartet: Data in ymm4: */\
	/**********************************/\
	"movslq		%[__p2],%%r15	\n\t"	\
	"prefetcht0	(%%r14,%%r15,8)	\n\t"\
		"movq	%[__sse_sw],%%rsi		\n\t"\
		"vmovaps	(%%rsi),	%%xmm11		\n\t"\
		"vpsubd		%%xmm8 ,%%xmm11,%%xmm11	\n\t"\
		"vmovmskps	%%xmm11,	%%rsi	\n\t"\
		"\n\t"\
		"movq	%[__n_minus_sil],%%rcx	\n\t"\
		"vmovd	0x8(%%rcx)	,%%xmm9 	\n\t"/* .d2 term of index quartet */\
		"vpshufd	$0,	%%xmm9 ,%%xmm9 		\n\t"\
		"vpsubd		%%xmm8 ,%%xmm9 ,%%xmm9 		\n\t"\
		"vmovmskps	%%xmm9 ,%%rcx		\n\t"\
		"\n\t"\
		"movq	%[__sinwt]	,%%rdx		\n\t"\
		"vmovd	0x8(%%rdx)	,%%xmm10	\n\t"/* .d2 term of index quartet */\
		"vpshufd	$0,	%%xmm10,%%xmm10		\n\t"\
		"vmovaps		%%xmm8 ,%%xmm11		\n\t"\
		"vpsubd		%%xmm10,%%xmm11,%%xmm11		\n\t"\
		"vmovmskps	%%xmm11,%%rdx		\n\t"\
		"\n\t"\
	"shlq	$5,%%rsi	\n\t"/* i0123 */\
	"shlq	$5,%%rcx	\n\t"/* m0123 */\
	"shlq	$5,%%rdx	\n\t"/* n0123 */\
		"\n\t"\
		"addq	%%rdi,%%rcx		\n\t"\
		"addq	%%rdi,%%rdx		\n\t"\
		"vbroadcastsd 0x840(%%rdi),%%ymm9 	\n\t"/* wtl */\
		"vbroadcastsd 0x848(%%rdi),%%ymm10	\n\t"/* wtn */\
		"vmulpd	%%ymm12,%%ymm9 ,%%ymm9 		\n\t"\
		"vmulpd	%%ymm13,%%ymm10,%%ymm10		\n\t"\
		"vmulpd	     (%%rcx),%%ymm9 ,%%ymm9 \n\t"\
		"vmulpd	0x200(%%rdx),%%ymm10,%%ymm10\n\t"\
		"\n\t"\
		"vmulpd		%%ymm10,%%ymm4,%%ymm4	\n\t"\
		"vaddpd		%%ymm14,%%ymm4,%%ymm4	\n\t"\
		"vmovaps	%%ymm4,%%ymm10			\n\t"\
		"vroundpd	$0,%%ymm10,%%ymm10		\n\t"\
		"\n\t"\
		"movq	%[__sign_mask],%%rbx		\n\t"\
		"vsubpd		%%ymm10,%%ymm4,%%ymm4	\n\t"\
		"vandpd		(%%rbx),%%ymm4,%%ymm4	\n\t"\
		"vmaxpd		%%ymm15,%%ymm4,%%ymm15	\n\t"\
		/* cy   = DNINT(temp*baseinv[i1]): */\
		"addq	%%rdi,%%rsi		\n\t"\
		"vmovaps	%%ymm10	,%%ymm4			\n\t"\
		"vmulpd	0x600(%%rsi),%%ymm10,%%ymm10\n\t"\
		"vroundpd	$0,%%ymm10,%%ymm14		\n\t"\
		"vmovaps	%%ymm14,%%ymm10				\n\t"/* cy = cpy cy_out */\
		/* x = (temp-cy*base[i1])*wt: */\
		"vmulpd	0x400(%%rsi),%%ymm10,%%ymm10	\n\t"/* cy*base[i0123] */\
		"vsubpd		%%ymm10,%%ymm4 ,%%ymm4 	\n\t"\
		"vmulpd		%%ymm9 ,%%ymm4 ,%%ymm4 	\n\t"\
		/* Get ready for next set [IM2~] by computing bjmod[0:3] += bw (mod n): */\
		"movq	%[__sse_n]	,%%rbx		\n\t"\
		"vmovaps		(%%rbx)	,%%xmm10	\n\t"\
		"movq	%[__sse_bw]	,%%rax		\n\t"\
		"vpaddd		(%%rax),%%xmm8 ,%%xmm8 	\n\t"\
		"vmovaps	%%xmm8 ,%%xmm9 		\n\t"\
		"vpcmpgtd	%%xmm10,%%xmm9 ,%%xmm9 	\n\t"\
		"vpand		%%xmm10,%%xmm9 ,%%xmm9 	\n\t"\
		"vpsubd		%%xmm9 ,%%xmm8 ,%%xmm8 	\n\t"\
		"\n\t"\
	/**********************************/\
	/* Do C.im-quartet: Data in ymm5: */\
	/**********************************/\
		"movq	%[__sse_sw],%%rsi		\n\t"\
		"vmovaps	(%%rsi),	%%xmm11		\n\t"\
		"vpsubd		%%xmm8 ,%%xmm11,%%xmm11	\n\t"\
		"vmovmskps	%%xmm11,	%%rsi	\n\t"\
		"\n\t"\
		"movq	%[__n_minus_silp1],%%rcx	\n\t"\
		"vmovd	0x8(%%rcx)	,%%xmm9 	\n\t"/* .d2 term of index quartet */\
		"vpshufd	$0,	%%xmm9 ,%%xmm9 		\n\t"\
		"vpsubd		%%xmm8 ,%%xmm9 ,%%xmm9 		\n\t"\
		"vmovmskps	%%xmm9 ,%%rcx		\n\t"\
		"\n\t"\
		"movq	%[__sinwtm1]	,%%rdx		\n\t"\
		"vmovd	0x8(%%rdx)	,%%xmm10	\n\t"/* .d2 term of index quartet */\
		"vpshufd	$0,	%%xmm10,%%xmm10		\n\t"\
		"vmovaps		%%xmm8 ,%%xmm11		\n\t"\
		"vpsubd		%%xmm10,%%xmm11,%%xmm11		\n\t"\
		"vmovmskps	%%xmm11,%%rdx		\n\t"\
		"\n\t"\
	"shlq	$5,%%rsi	\n\t"/* i0123 */\
	"shlq	$5,%%rcx	\n\t"/* m0123 */\
	"shlq	$5,%%rdx	\n\t"/* n0123 */\
		"\n\t"\
		"addq	%%rdi,%%rcx		\n\t"\
		"addq	%%rdi,%%rdx		\n\t"\
		"vbroadcastsd 0x850(%%rdi),%%ymm9 	\n\t"/* wtlp1 */\
		"vbroadcastsd 0x858(%%rdi),%%ymm10	\n\t"/* wtnm1 */\
		"vmulpd	%%ymm12,%%ymm9 ,%%ymm9 		\n\t"\
		"vmulpd	%%ymm13,%%ymm10,%%ymm10		\n\t"\
		"vmulpd	     (%%rcx),%%ymm9 ,%%ymm9 \n\t"\
		"vmulpd	0x200(%%rdx),%%ymm10,%%ymm10\n\t"\
		"\n\t"\
		"vmulpd		%%ymm10,%%ymm5,%%ymm5	\n\t"\
		"vaddpd		%%ymm14,%%ymm5,%%ymm5	\n\t"\
		"vmovaps	%%ymm5,%%ymm10			\n\t"\
		"vroundpd	$0,%%ymm10,%%ymm10		\n\t"\
		"\n\t"\
		"movq	%[__sign_mask],%%rbx		\n\t"\
		"vsubpd		%%ymm10,%%ymm5,%%ymm5	\n\t"\
		"vandpd		(%%rbx),%%ymm5,%%ymm5	\n\t"\
		"vmaxpd		%%ymm15,%%ymm5,%%ymm15	\n\t"\
		/* cy   = DNINT(temp*baseinv[i1]): */\
		"addq	%%rdi,%%rsi		\n\t"\
		"vmovaps	%%ymm10	,%%ymm5			\n\t"\
		"vmulpd	0x600(%%rsi),%%ymm10,%%ymm10\n\t"\
		"vroundpd	$0,%%ymm10,%%ymm14		\n\t"\
		"vmovaps	%%ymm14,%%ymm10				\n\t"/* cy = cpy cy_out */\
		/* x = (temp-cy*base[i1])*wt: */\
		"vmulpd	0x400(%%rsi),%%ymm10,%%ymm10	\n\t"/* cy*base[i0123] */\
		"vsubpd		%%ymm10,%%ymm5 ,%%ymm5 	\n\t"\
		"vmulpd		%%ymm9 ,%%ymm5 ,%%ymm5 	\n\t"\
		/* Get ready for next set [RE3~] : by computing bjmod[0:3] += bw (mod n): */\
		"movq	%[__sse_n]	,%%rbx		\n\t"\
		"vmovaps		(%%rbx)	,%%xmm10	\n\t"\
		"movq	%[__sse_bw]	,%%rax		\n\t"\
		"vpaddd		(%%rax),%%xmm8 ,%%xmm8 	\n\t"\
		"vmovaps	%%xmm8 ,%%xmm9 		\n\t"\
		"vpcmpgtd	%%xmm10,%%xmm9 ,%%xmm9 	\n\t"\
		"vpand		%%xmm10,%%xmm9 ,%%xmm9 	\n\t"\
		"vpsubd		%%xmm9 ,%%xmm8 ,%%xmm8 	\n\t"\
		"\n\t"\
	/**********************************/\
	/* Do D.re-quartet: Data in ymm6: */\
	/**********************************/\
	"movslq		%[__p3],%%r15	\n\t"	\
	"prefetcht0	(%%r14,%%r15,8)	\n\t"\
		"movq	%[__sse_sw],%%rsi		\n\t"\
		"vmovaps	(%%rsi),	%%xmm11		\n\t"\
		"vpsubd		%%xmm8 ,%%xmm11,%%xmm11	\n\t"\
		"vmovmskps	%%xmm11,	%%rsi	\n\t"\
		"\n\t"\
		"movq	%[__n_minus_sil],%%rcx	\n\t"\
		"vmovd	0xC(%%rcx)	,%%xmm9 	\n\t"/* .d3 term of index quartet */\
		"vpshufd	$0,	%%xmm9 ,%%xmm9 		\n\t"\
		"vpsubd		%%xmm8 ,%%xmm9 ,%%xmm9 		\n\t"\
		"vmovmskps	%%xmm9 ,%%rcx		\n\t"\
		"\n\t"\
		"movq	%[__sinwt]	,%%rdx		\n\t"\
		"vmovd	0xC(%%rdx)	,%%xmm10	\n\t"/* .d3 term of index quartet */\
		"vpshufd	$0,	%%xmm10,%%xmm10		\n\t"\
		"vmovaps		%%xmm8 ,%%xmm11		\n\t"\
		"vpsubd		%%xmm10,%%xmm11,%%xmm11		\n\t"\
		"vmovmskps	%%xmm11,%%rdx		\n\t"\
		"\n\t"\
	"shlq	$5,%%rsi	\n\t"/* i0123 */\
	"shlq	$5,%%rcx	\n\t"/* m0123 */\
	"shlq	$5,%%rdx	\n\t"/* n0123 */\
		"\n\t"\
		"addq	%%rdi,%%rcx		\n\t"\
		"addq	%%rdi,%%rdx		\n\t"\
		"vbroadcastsd 0x860(%%rdi),%%ymm9 	\n\t"/* wtl */\
		"vbroadcastsd 0x868(%%rdi),%%ymm10	\n\t"/* wtn */\
		"vmulpd	%%ymm12,%%ymm9 ,%%ymm9 		\n\t"\
		"vmulpd	%%ymm13,%%ymm10,%%ymm10		\n\t"\
		"vmulpd	     (%%rcx),%%ymm9 ,%%ymm9 \n\t"\
		"vmulpd	0x200(%%rdx),%%ymm10,%%ymm10\n\t"\
		"\n\t"\
		"vmulpd		%%ymm10,%%ymm6,%%ymm6	\n\t"\
		"vaddpd		%%ymm14,%%ymm6,%%ymm6	\n\t"\
		"vmovaps	%%ymm6,%%ymm10			\n\t"\
		"vroundpd	$0,%%ymm10,%%ymm10		\n\t"\
		"\n\t"\
		"movq	%[__sign_mask],%%rbx		\n\t"\
		"vsubpd		%%ymm10,%%ymm6,%%ymm6	\n\t"\
		"vandpd		(%%rbx),%%ymm6,%%ymm6	\n\t"\
		"vmaxpd		%%ymm15,%%ymm6,%%ymm15	\n\t"\
		/* cy   = DNINT(temp*baseinv[i1]): */\
		"addq	%%rdi,%%rsi		\n\t"\
		"vmovaps	%%ymm10	,%%ymm6			\n\t"\
		"vmulpd	0x600(%%rsi),%%ymm10,%%ymm10\n\t"\
		"vroundpd	$0,%%ymm10,%%ymm14		\n\t"\
		"vmovaps	%%ymm14,%%ymm10				\n\t"/* cy = cpy cy_out */\
		/* x = (temp-cy*base[i1])*wt: */\
		"vmulpd	0x400(%%rsi),%%ymm10,%%ymm10	\n\t"/* cy*base[i0123] */\
		"vsubpd		%%ymm10,%%ymm6 ,%%ymm6 	\n\t"\
		"vmulpd		%%ymm9 ,%%ymm6 ,%%ymm6 	\n\t"\
		/* Get ready for next set [IM3~] : by computing bjmod[0:3] += bw (mod n): */\
		"movq	%[__sse_n]	,%%rbx		\n\t"\
		"vmovaps		(%%rbx)	,%%xmm10	\n\t"\
		"movq	%[__sse_bw]	,%%rax		\n\t"\
		"vpaddd		(%%rax),%%xmm8 ,%%xmm8 	\n\t"\
		"vmovaps	%%xmm8 ,%%xmm9 		\n\t"\
		"vpcmpgtd	%%xmm10,%%xmm9 ,%%xmm9 	\n\t"\
		"vpand		%%xmm10,%%xmm9 ,%%xmm9 	\n\t"\
		"vpsubd		%%xmm9 ,%%xmm8 ,%%xmm8 	\n\t"\
		"\n\t"\
	/**********************************/\
	/* Do D.im-quartet: Data in ymm7: */\
	/**********************************/\
		"movq	%[__sse_sw],%%rsi		\n\t"\
		"vmovaps	(%%rsi),	%%xmm11		\n\t"\
		"vpsubd		%%xmm8 ,%%xmm11,%%xmm11	\n\t"\
		"vmovmskps	%%xmm11,	%%rsi	\n\t"\
		"\n\t"\
		"movq	%[__n_minus_silp1],%%rcx	\n\t"\
		"vmovd	0xC(%%rcx)	,%%xmm9 	\n\t"/* .d3 term of index quartet */\
		"vpshufd	$0,	%%xmm9 ,%%xmm9 		\n\t"\
		"vpsubd		%%xmm8 ,%%xmm9 ,%%xmm9 		\n\t"\
		"vmovmskps	%%xmm9 ,%%rcx		\n\t"\
		"\n\t"\
		"movq	%[__sinwtm1]	,%%rdx		\n\t"\
		"vmovd	0xC(%%rdx)	,%%xmm10	\n\t"/* .d3 term of index quartet */\
		"vpshufd	$0,	%%xmm10,%%xmm10		\n\t"\
		"vmovaps		%%xmm8 ,%%xmm11		\n\t"\
		"vpsubd		%%xmm10,%%xmm11,%%xmm11		\n\t"\
		"vmovmskps	%%xmm11,%%rdx		\n\t"\
		"\n\t"\
	"shlq	$5,%%rsi	\n\t"/* i0123 */\
	"shlq	$5,%%rcx	\n\t"/* m0123 */\
	"shlq	$5,%%rdx	\n\t"/* n0123 */\
		"\n\t"\
		"addq	%%rdi,%%rcx		\n\t"\
		"addq	%%rdi,%%rdx		\n\t"\
		"vbroadcastsd 0x870(%%rdi),%%ymm9 	\n\t"/* wtlp1 */\
		"vbroadcastsd 0x878(%%rdi),%%ymm10	\n\t"/* wtnm1 */\
		"vmulpd	%%ymm12,%%ymm9 ,%%ymm9 		\n\t"\
		"vmulpd	%%ymm13,%%ymm10,%%ymm10		\n\t"\
		"vmulpd	     (%%rcx),%%ymm9 ,%%ymm9 \n\t"\
		"vmulpd	0x200(%%rdx),%%ymm10,%%ymm10\n\t"\
		"\n\t"\
		"vmulpd		%%ymm10,%%ymm7,%%ymm7	\n\t"\
		"vaddpd		%%ymm14,%%ymm7,%%ymm7	\n\t"\
		"vmovaps	%%ymm7,%%ymm10			\n\t"\
		"vroundpd	$0,%%ymm10,%%ymm10		\n\t"\
		"\n\t"\
		"movq	%[__sign_mask],%%rbx		\n\t"\
		"vsubpd		%%ymm10,%%ymm7,%%ymm7	\n\t"\
		"vandpd		(%%rbx),%%ymm7,%%ymm7	\n\t"\
		"vmaxpd		%%ymm15,%%ymm7,%%ymm15	\n\t"\
		/* cy   = DNINT(temp*baseinv[i1]): */\
		"addq	%%rdi,%%rsi		\n\t"\
		"vmovaps	%%ymm10	,%%ymm7			\n\t"\
		"vmulpd	0x600(%%rsi),%%ymm10,%%ymm10\n\t"\
		"vroundpd	$0,%%ymm10,%%ymm14		\n\t"\
		"vmovaps	%%ymm14,%%ymm10				\n\t"/* cy = cpy cy_out */\
		/* x = (temp-cy*base[i1])*wt: */\
		"vmulpd	0x400(%%rsi),%%ymm10,%%ymm10	\n\t"/* cy*base[i0123] */\
		"vsubpd		%%ymm10,%%ymm7 ,%%ymm7 	\n\t"\
		"vmulpd		%%ymm9 ,%%ymm7 ,%%ymm7 	\n\t"\
		/* Update wts-array pointers in preparation for next call of the macro: */\
		"movq	%[__wtA]	,%%rax		\n\t"\
		"movq	%[__wtB]	,%%rbx		\n\t"\
		"movq	%[__wtC]	,%%rcx		\n\t"\
		"addq	$0x20	,%%rax			\n\t"/* add0 += 4 */\
		"subq	$0x20	,%%rbx			\n\t"/* add1 -= 4 */\
		"subq	$0x20	,%%rcx			\n\t"/* add2 -= 4 */\
		"movq	%%rax	,%[__wtA]		\n\t"\
		"movq	%%rbx	,%[__wtB]		\n\t"\
		"movq	%%rcx	,%[__wtC]		\n\t"\
		"\n\t"\
		/* Get ready for store of final-updated bjmod[0:3] values: */\
		"movq	%[__sse_n]	,%%rbx		\n\t"\
		"vmovaps		(%%rbx)	,%%xmm10	\n\t"\
		"movq	%[__sse_bw]	,%%rax		\n\t"\
		"vpaddd		(%%rax),%%xmm8 ,%%xmm8 	\n\t"\
		"vmovaps	%%xmm8 ,%%xmm9 		\n\t"\
		"vpcmpgtd	%%xmm10,%%xmm9 ,%%xmm9 	\n\t"\
		"vpand		%%xmm10,%%xmm9 ,%%xmm9 	\n\t"\
		"vpsubd		%%xmm9 ,%%xmm8 ,%%xmm8 	\n\t"\
		"\n\t"\
		/* Store bjmodn index quartet: */\
		"movq	%[__bjmod_0],%%rcx		\n\t"\
		"vmovaps	%%xmm8,(%%rcx)			\n\t"\
		/* Store cy_out: */\
		"movq		%[__cy],%%rbx	\n\t"\
		"vmovaps	%%ymm14,(%%rbx)	\n\t"/* ymm14 = cy_in */\
		/* Store maxerr: */\
		"movq		%[__half_arr],%%rdx		\n\t"\
		"vmovaps	%%ymm15,-0x40(%%rdx)	\n\t"\
	/* 4-way transpose of outputs (Re, Im parts separately): Inputs from ymm0-7. Outputs into r0/1,2/3,4/5.6/7: */	\
	/* Because default inputs for our 4 x 4 transpose macro (e.g. the one used at start of this carry macro) */\
	/* are into ymm4/2/8/2, munge inputs into that order, resolving name-conflicts via use of the now-available ymm8-15 for outputs: */\
		"movq		%[__data],%%rax			\n\t"\
		"vshufpd	$15,%%ymm2,%%ymm0,%%ymm10					\n\t		vshufpd	$15,%%ymm3,%%ymm1,%%ymm11						\n\t"\
		"vshufpd	$0 ,%%ymm2,%%ymm0,%%ymm0					\n\t		vshufpd	$0 ,%%ymm3,%%ymm1,%%ymm1						\n\t"\
		"vshufpd	$15,%%ymm6,%%ymm4,%%ymm12					\n\t		vshufpd	$15,%%ymm7,%%ymm5,%%ymm13						\n\t"\
		"vshufpd	$0 ,%%ymm6,%%ymm4,%%ymm4					\n\t		vshufpd	$0 ,%%ymm7,%%ymm5,%%ymm5						\n\t"\
		"vperm2f128 $32,%%ymm12,%%ymm10,%%ymm2 		/* Re B	*/	\n\t		vperm2f128 $32,%%ymm13,%%ymm11,%%ymm3		/* Im B	*/	\n\t"\
		"vperm2f128 $49,%%ymm12,%%ymm10,%%ymm10		/* Re D	*/	\n\t		vperm2f128 $49,%%ymm13,%%ymm11,%%ymm11		/* Im D	*/	\n\t"\
		"vperm2f128 $32,%%ymm4 ,%%ymm0 ,%%ymm12		/* Re A	*/	\n\t		vperm2f128 $32,%%ymm5 ,%%ymm1 ,%%ymm13 		/* Im A	*/	\n\t"\
		"vperm2f128 $49,%%ymm4 ,%%ymm0 ,%%ymm0 		/* Re C	*/	\n\t		vperm2f128 $49,%%ymm5 ,%%ymm1 ,%%ymm1		/* Im C	*/	\n\t"\
		/* And write 'em back to memory: */\
		"vmovaps	%%ymm12,     (%%rax)						\n\t		vmovaps	%%ymm13,0x020(%%rax)				\n\t"\
		"vmovaps	%%ymm2 ,0x040(%%rax)						\n\t		vmovaps	%%ymm3 ,0x060(%%rax)				\n\t"\
		"vmovaps	%%ymm0 ,0x080(%%rax)						\n\t		vmovaps	%%ymm1 ,0x0a0(%%rax)				\n\t"\
		"vmovaps	%%ymm10,0x0c0(%%rax)						\n\t		vmovaps	%%ymm11,0x0e0(%%rax)				\n\t"\
		:					/* outputs: none */\
		: [__data]		"m" (Xdata)	/* All inputs from memory addresses here */\
		, [__wtA]		"m" (XwtA)		\
		, [__wtB]		"m" (XwtB)		\
		, [__wtC]		"m" (XwtC)		\
		, [__cy]		"m" (Xcy)		\
		, [__bjmod_0]	"m" (Xbjmod_0)		\
		, [__half_arr]	"m" (Xhalf_arr)		\
		, [__n_minus_silp1] "m" (Xn_minus_silp1)\
		, [__n_minus_sil]	"m" (Xn_minus_sil)	\
		, [__sign_mask]	"m" (Xsign_mask)	\
		, [__sinwt]		"m" (Xsinwt)		\
		, [__sinwtm1]	"m" (Xsinwtm1)		\
		, [__sse_bw]	"m" (Xsse_bw)		\
		, [__sse_n]		"m" (Xsse_n)		\
		, [__sse_sw]	"m" (Xsse_sw)		\
		/* Prefetch: base address and 3 index offsets */\
		,	[__add0] "m" (Xadd0)\
		,	[__p1] "m" (Xp1)\
		,	[__p2] "m" (Xp2)\
		,	[__p3] "m" (Xp3)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r14","r15","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
	);\
	}

#else

	/********* Packed 32-bit-int version of SSE2_cmplx_carry_norm_pow2_errcheck0_2x:***********/
	#define SSE2_cmplx_carry_norm_pow2_errcheck0_2B(Xdata,XwtA,XwtB,XwtC,XcyA,XcyB,Xbjmod_0,Xhalf_arr,Xi,Xn_minus_silp1,Xn_minus_sil,Xsign_mask,Xsinwt,Xsinwtm1,Xsse_bw,Xsse_nm1,Xsse_sw, Xadd0,Xp1)\
	{\
	__asm__ volatile (\
	"movq	%[__add0],%%r14	\n\t"/* base address for 2 prefetches-from-main-data-array spread through this macro */\
	"prefetcht0	(%%r14)		\n\t"\
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
		"movaps		-0x10(%%rax),%%xmm4	\n\t"\
		"movaps		     (%%rbx),%%xmm2	\n\t	movaps		 0x10(%%rbx),%%xmm6	\n\t"\
		"movhpd		     (%%rcx),%%xmm3	\n\t	movhpd		-0x10(%%rcx),%%xmm7	\n\t"\
		"movlpd		 0x08(%%rcx),%%xmm3	\n\t	movlpd		-0x08(%%rcx),%%xmm7	\n\t"\
	/* The commented-out-below variant of the above movhpd/movlpd combo is no faster in SSE2 mode, but generalizes better to AVX: */\
	"/*	movaps		     (%%rcx),%%xmm3	\n\t	movaps		-0x10(%%rcx),%%xmm7	\n\t	*/"\
	"/*	shufpd		$1 ,%%xmm3,%%xmm3	\n\t	shufpd		$1 ,%%xmm7,%%xmm7	\n\t	*/"\
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
		"/*roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd		$0,%%xmm7,%%xmm7	***Requires SSE4.1, any tiny speed gain not worth loss of portability***/\n\t"\
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
		"/*roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd		$0,%%xmm7,%%xmm7	*/\n\t"\
		"movaps		%%xmm3	,(%%rcx)	\n\t	movaps		%%xmm7	,(%%rdx)	\n\t"\
		"\n\t"\
		"movq	%[__data]	,%%rax		\n\t"\
		"mulpd	 0x80(%%rdi)	,%%xmm3	\n\t	mulpd	 0x80(%%rbx)	,%%xmm7	\n\t"\
		"subpd		%%xmm3	,	%%xmm1	\n\t	subpd		%%xmm7	,	%%xmm5	\n\t"\
		"mulpd		%%xmm2	,	%%xmm1	\n\t	mulpd		%%xmm6	,	%%xmm5	\n\t"\
		"movaps		%%xmm1	,    (%%rax)\n\t	movaps		%%xmm5	,0x40(%%rax)\n\t"\
		"\n\t"\
		"movq	%[__sse_bw]	,%%rax		\n\t"\
		"movq	%[__sse_nm1]	,%%rbx	\n\t"\
		"paddd		(%%rax)	,%%xmm0		\n\t"\
		"pand		(%%rbx)	,%%xmm0		\n\t"\
	"/**********************************************/\n\t"\
	"/*          Imaginary parts               */\n\t"\
	"/**********************************************/\n\t"\
	"movslq		%[__p1],%%r15	\n\t"	\
	"prefetcht0	(%%r14,%%r15,8)	\n\t"\
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
		"movhpd		     (%%rcx),%%xmm3	\n\t	movhpd		-0x10(%%rcx),%%xmm7	\n\t"\
		"movlpd		 0x08(%%rcx),%%xmm3	\n\t	movlpd		-0x08(%%rcx),%%xmm7	\n\t"\
	"/*	movaps		     (%%rcx),%%xmm3	\n\t	movaps		-0x10(%%rcx),%%xmm7	\n\t	*/"\
	"/*	shufpd		$1 ,%%xmm3,%%xmm3	\n\t	shufpd		$1 ,%%xmm7,%%xmm7	\n\t	*/"\
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
		"addpd		%%xmm4	,%%xmm3		\n\t	addpd		%%xmm4	,%%xmm7		\n\t"\
		"subpd		%%xmm4	,%%xmm3		\n\t	subpd		%%xmm4	,%%xmm7		\n\t"\
		"/*roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd		$0,%%xmm7,%%xmm7	*/\n\t"\
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
		"addpd		%%xmm4	,%%xmm3		\n\t	addpd		%%xmm4	,%%xmm7		\n\t"\
		"subpd		%%xmm4	,%%xmm3		\n\t	subpd		%%xmm4	,%%xmm7		\n\t"\
		"/*roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd		$0,%%xmm7,%%xmm7	*/\n\t"\
		"movaps	%%xmm3	,(%%rcx)		\n\t	movaps	%%xmm7	,(%%rdx)		\n\t"\
		"\n\t"\
		"movq	%[__data]	,%%rax		\n\t"\
		"mulpd	 0x80(%%rdi)	,%%xmm3	\n\t	mulpd	 0x80(%%rbx)	,%%xmm7	\n\t"\
		"subpd	%%xmm3	,%%xmm1			\n\t	subpd	%%xmm7	,%%xmm5			\n\t"\
		"mulpd	%%xmm2	,%%xmm1			\n\t	mulpd	%%xmm6	,%%xmm5			\n\t"\
		"movaps	%%xmm1	, 0x10(%%rax)	\n\t	movaps	%%xmm5	, 0x50(%%rax)	\n\t"\
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
		/* Prefetch: base address and 1 index offset */\
		,	[__add0] "m" (Xadd0)\
		,	[__p1] "m" (Xp1)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r14","r15","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"	/* Clobbered registers */\
	);\
	}

	/********* Packed 32-bit-int version of SSE2_cmplx_carry_norm_pow2_errcheck1_2x:***********/
	#define SSE2_cmplx_carry_norm_pow2_errcheck1_2B(Xdata,XwtA,XwtB,XwtC,XcyA,XcyB,Xbjmod_0,Xhalf_arr,Xn_minus_silp1,Xn_minus_sil,Xsign_mask,Xsinwt,Xsinwtm1,Xsse_bw,Xsse_nm1,Xsse_sw, Xadd0,Xp1)\
	{\
	__asm__ volatile (\
	"movq	%[__add0],%%r14	\n\t"/* base address for 2 prefetches-from-main-data-array spread through this macro */\
	"prefetcht0	(%%r14)		\n\t"\
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
		"/*roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd		$0,%%xmm7,%%xmm7	*/\n\t"\
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
		"addpd		%%xmm4	,%%xmm3		\n\t	addpd		%%xmm4	,%%xmm7		\n\t"\
		"subpd		%%xmm4	,%%xmm3		\n\t	subpd		%%xmm4	,%%xmm7		\n\t"\
		"/*roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd		$0,%%xmm7,%%xmm7	*/\n\t"\
		"movaps		%%xmm3	,(%%rcx)	\n\t	movaps		%%xmm7	,(%%rdx)	\n\t"\
		"\n\t"\
		"movq	%[__data]	,%%rax		\n\t"\
		"mulpd	 0x80(%%rdi)	,%%xmm3	\n\t	mulpd	 0x80(%%rbx)	,%%xmm7	\n\t"\
		"subpd		%%xmm3	,	%%xmm1	\n\t	subpd		%%xmm7	,	%%xmm5	\n\t"\
		"mulpd		%%xmm2	,	%%xmm1	\n\t	mulpd		%%xmm6	,	%%xmm5	\n\t"\
		"movaps		%%xmm1	,    (%%rax)\n\t	movaps		%%xmm5	,0x40(%%rax)\n\t"\
		"\n\t"\
		"movq	%[__sse_bw]	,%%rax		\n\t"\
		"movq	%[__sse_nm1]	,%%rbx	\n\t"\
		"paddd		(%%rax)	,%%xmm0		\n\t"\
		"pand		(%%rbx)	,%%xmm0		\n\t"\
	"/**********************************************/\n\t"\
	"/*          Imaginary parts               */\n\t"\
	"/**********************************************/\n\t"\
	"movslq		%[__p1],%%r15	\n\t"	\
	"prefetcht0	(%%r14,%%r15,8)	\n\t"\
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
		"movhpd		     (%%rcx),%%xmm3	\n\t	movhpd		-0x10(%%rcx),%%xmm7	\n\t"\
		"movlpd		 0x08(%%rcx),%%xmm3	\n\t	movlpd		-0x08(%%rcx),%%xmm7	\n\t"\
	"/*	movaps		     (%%rcx),%%xmm3	\n\t	movaps		-0x10(%%rcx),%%xmm7	\n\t	*/"\
	"/*	shufpd		$1 ,%%xmm3,%%xmm3	\n\t	shufpd		$1 ,%%xmm7,%%xmm7	\n\t	*/"\
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
		"addpd		%%xmm4	,%%xmm3		\n\t	addpd		%%xmm4	,%%xmm7		\n\t"\
		"subpd		%%xmm4	,%%xmm3		\n\t	subpd		%%xmm4	,%%xmm7		\n\t"\
		"/*roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd		$0,%%xmm7,%%xmm7	*/\n\t"\
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
		"addpd		%%xmm4	,%%xmm3		\n\t	addpd		%%xmm4	,%%xmm7		\n\t"\
		"subpd		%%xmm4	,%%xmm3		\n\t	subpd		%%xmm4	,%%xmm7		\n\t"\
		"/*roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd		$0,%%xmm7,%%xmm7	*/\n\t"\
		"movaps	%%xmm3	,(%%rcx)		\n\t	movaps	%%xmm7	,(%%rdx)		\n\t"\
		"\n\t"\
		"movq	%[__data]	,%%rax		\n\t"\
		"mulpd	 0x80(%%rdi)	,%%xmm3	\n\t	mulpd	 0x80(%%rbx)	,%%xmm7	\n\t"\
		"subpd	%%xmm3	,%%xmm1			\n\t	subpd	%%xmm7	,%%xmm5			\n\t"\
		"mulpd	%%xmm2	,%%xmm1			\n\t	mulpd	%%xmm6	,%%xmm5			\n\t"\
		"movaps	%%xmm1	, 0x10(%%rax)	\n\t	movaps	%%xmm5	, 0x50(%%rax)	\n\t"\
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
		/* Prefetch: base address and 1 index offset */\
		,	[__add0] "m" (Xadd0)\
		,	[__p1] "m" (Xp1)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r14","r15","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
	);\
	}

	/********* Packed 32-bit-int version of SSE2_cmplx_carry_norm_pow2_errcheck2_2x:***********/
	#define SSE2_cmplx_carry_norm_pow2_errcheck2_2B(Xdata,XwtA,XwtB,XcyA,XcyB,Xbjmod_0,Xhalf_arr,Xn_minus_silp1,Xn_minus_sil,Xsign_mask,Xsinwt,Xsinwtm1,Xsse_bw,Xsse_nm1,Xsse_sw, Xadd0,Xp2,Xp3)\
	{\
	__asm__ volatile (\
	"movq	%[__add0],%%r14	\n\t"/* base address for 2 prefetches-from-main-data-array spread through this macro */\
	"movslq		%[__p2],%%r15	\n\t"	\
	"prefetcht0	(%%r14,%%r15,8)	\n\t"\
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
		"movaps		-0x10(%%rax),%%xmm4	\n\t"\
		"movaps		     (%%rbx),%%xmm2	\n\t	movaps		 0x10(%%rbx),%%xmm6	\n\t"\
		"movhpd		     (%%rcx),%%xmm3	\n\t	movhpd		-0x10(%%rcx),%%xmm7	\n\t"\
		"movlpd		 0x08(%%rcx),%%xmm3	\n\t	movlpd		-0x08(%%rcx),%%xmm7	\n\t"\
	"/*	movaps		     (%%rcx),%%xmm3	\n\t	movaps		-0x10(%%rcx),%%xmm7	\n\t	*/"\
	"/*	shufpd		$1 ,%%xmm3,%%xmm3	\n\t	shufpd		$1 ,%%xmm7,%%xmm7	\n\t	*/"\
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
		"/*roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd		$0,%%xmm7,%%xmm7	*/\n\t"\
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
		"/*roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd		$0,%%xmm7,%%xmm7	*/\n\t"\
		"movaps		%%xmm3	,(%%rcx)	\n\t	movaps		%%xmm7	,(%%rdx)	\n\t"\
		"\n\t"\
		"movq	%[__data]	,%%rax		\n\t"\
		"mulpd	 0x80(%%rdi)	,%%xmm3	\n\t	mulpd	 0x80(%%rbx)	,%%xmm7	\n\t"\
		"subpd		%%xmm3	,	%%xmm1	\n\t	subpd		%%xmm7	,	%%xmm5	\n\t"\
		"mulpd		%%xmm2	,	%%xmm1	\n\t	mulpd		%%xmm6	,	%%xmm5	\n\t"\
		"movaps		%%xmm1	,0x20(%%rax)\n\t	movaps		%%xmm5	,0x60(%%rax)\n\t"\
		"\n\t"\
		"movq	%[__sse_bw]	,%%rax		\n\t"\
		"movq	%[__sse_nm1]	,%%rbx	\n\t"\
		"paddd		(%%rax)	,%%xmm0		\n\t"\
		"pand		(%%rbx)	,%%xmm0		\n\t"\
	"/**********************************************/\n\t"\
	"/*          Imaginary parts               */\n\t"\
	"/**********************************************/\n\t"\
	"movslq		%[__p3],%%r15	\n\t"	\
	"prefetcht0	(%%r14,%%r15,8)	\n\t"\
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
		"movhpd		     (%%rcx),%%xmm3	\n\t	movhpd		-0x10(%%rcx),%%xmm7	\n\t"\
		"movlpd		 0x08(%%rcx),%%xmm3	\n\t	movlpd		-0x08(%%rcx),%%xmm7	\n\t"\
	"/*	movaps		     (%%rcx),%%xmm3	\n\t	movaps		-0x10(%%rcx),%%xmm7	\n\t	*/"\
	"/*	shufpd		$1 ,%%xmm3,%%xmm3	\n\t	shufpd		$1 ,%%xmm7,%%xmm7	\n\t	*/"\
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
		"addpd		%%xmm4	,%%xmm3		\n\t	addpd		%%xmm4	,%%xmm7		\n\t"\
		"subpd		%%xmm4	,%%xmm3		\n\t	subpd		%%xmm4	,%%xmm7		\n\t"\
		"/*roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd		$0,%%xmm7,%%xmm7	*/\n\t"\
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
		"addpd		%%xmm4	,%%xmm3		\n\t	addpd		%%xmm4	,%%xmm7		\n\t"\
		"subpd		%%xmm4	,%%xmm3		\n\t	subpd		%%xmm4	,%%xmm7		\n\t"\
		"/*roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd		$0,%%xmm7,%%xmm7	*/\n\t"\
		"movaps	%%xmm3	,(%%rcx)		\n\t	movaps	%%xmm7	,(%%rdx)		\n\t"\
		"\n\t"\
		"movq	%[__data]	,%%rax		\n\t"\
		"mulpd	 0x80(%%rdi)	,%%xmm3	\n\t	mulpd	 0x80(%%rbx)	,%%xmm7	\n\t"\
		"subpd	%%xmm3	,%%xmm1			\n\t	subpd	%%xmm7	,%%xmm5			\n\t"\
		"mulpd	%%xmm2	,%%xmm1			\n\t	mulpd	%%xmm6	,%%xmm5			\n\t"\
		"movaps	%%xmm1	, 0x30(%%rax)	\n\t	movaps	%%xmm5	, 0x70(%%rax)	\n\t"\
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
		/* Prefetch: base address and 2 index offsets */\
		,	[__add0] "m" (Xadd0)\
		,	[__p2] "m" (Xp2)\
		,	[__p3] "m" (Xp3)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r14","r15","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
	);\
	}

#endif	// #ifdef USE_AVX

	/******************************************************************************************************************************************************************/
	/********** No-ROE-Check versions of the latter 2 of the above 3 macros - the first is called too infrequently to bother with a special non-ROE version: **********/
	/******************************************************************************************************************************************************************/

	#define SSE2_cmplx_carry_norm_pow2_nocheck1_2B(Xdata,XwtA,XwtB,XwtC,XcyA,XcyB,Xbjmod_0,Xhalf_arr,Xn_minus_silp1,Xn_minus_sil,Xsinwt,Xsinwtm1,Xsse_bw,Xsse_nm1,Xsse_sw, Xadd0,Xp1)\
	{\
	__asm__ volatile (\
	"movq	%[__add0],%%r14	\n\t"/* base address for 2 prefetches-from-main-data-array spread through this macro */\
	"prefetcht0	(%%r14)		\n\t"\
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
		"/*roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd		$0,%%xmm7,%%xmm7	*/\n\t"\
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
		"/*roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd		$0,%%xmm7,%%xmm7	*/\n\t"\
		"movaps		%%xmm3	,(%%rcx)	\n\t	movaps		%%xmm7	,(%%rdx)	\n\t"\
		"\n\t"\
		"movq	%[__data]	,%%rax		\n\t"\
		"mulpd	 0x80(%%rdi)	,%%xmm3	\n\t	mulpd	 0x80(%%rbx)	,%%xmm7	\n\t"\
		"subpd		%%xmm3	,	%%xmm1	\n\t	subpd		%%xmm7	,	%%xmm5	\n\t"\
		"mulpd		%%xmm2	,	%%xmm1	\n\t	mulpd		%%xmm6	,	%%xmm5	\n\t"\
		"movaps		%%xmm1	,    (%%rax)\n\t	movaps		%%xmm5	,0x40(%%rax)\n\t"\
		"\n\t"\
		"movq	%[__sse_bw]	,%%rax		\n\t"\
		"movq	%[__sse_nm1]	,%%rbx	\n\t"\
		"paddd		(%%rax)	,%%xmm0		\n\t"\
		"pand		(%%rbx)	,%%xmm0		\n\t"\
	"/**********************************************/\n\t"\
	"/*          Imaginary parts               */\n\t"\
	"/**********************************************/\n\t"\
	"movslq		%[__p1],%%r15	\n\t"	\
	"prefetcht0	(%%r14,%%r15,8)	\n\t"\
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
		"addpd		%%xmm4	,%%xmm3		\n\t	addpd		%%xmm4	,%%xmm7		\n\t"\
		"subpd		%%xmm4	,%%xmm3		\n\t	subpd		%%xmm4	,%%xmm7		\n\t"\
		"/*roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd		$0,%%xmm7,%%xmm7	*/\n\t"\
		"\n\t"\
		"/* NO ROE HERE */		\n\t"\
		"\n\t"\
		"movq	%%rsi,	%%rdi	\n\t	movq	%%rsi,	%%rbx	\n\t"\
		"shrq	$20,	%%rdi	\n\t	shrq	$22,	%%rbx	\n\t"\
		"andq	$0x0000000000000030	,%%rdi		\n\t	andq	$0x0000000000000030	,%%rbx		\n\t"\
		"addq	%%rax	,%%rdi			\n\t	addq	%%rax	,%%rbx			\n\t"\
		"movaps	%%xmm3	,%%xmm1			\n\t	movaps	%%xmm7	,%%xmm5			\n\t"\
		"mulpd	 0xc0(%%rdi)	,%%xmm3	\n\t	mulpd	 0xc0(%%rbx)	,%%xmm7	\n\t"\
		"addpd		%%xmm4	,%%xmm3		\n\t	addpd		%%xmm4	,%%xmm7		\n\t"\
		"subpd		%%xmm4	,%%xmm3		\n\t	subpd		%%xmm4	,%%xmm7		\n\t"\
		"/*roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd		$0,%%xmm7,%%xmm7	*/\n\t"\
		"movaps	%%xmm3	,(%%rcx)		\n\t	movaps	%%xmm7	,(%%rdx)		\n\t"\
		"\n\t"\
		"movq	%[__data]	,%%rax		\n\t"\
		"mulpd	 0x80(%%rdi)	,%%xmm3	\n\t	mulpd	 0x80(%%rbx)	,%%xmm7	\n\t"\
		"subpd	%%xmm3	,%%xmm1			\n\t	subpd	%%xmm7	,%%xmm5			\n\t"\
		"mulpd	%%xmm2	,%%xmm1			\n\t	mulpd	%%xmm6	,%%xmm5			\n\t"\
		"movaps	%%xmm1	, 0x10(%%rax)	\n\t	movaps	%%xmm5	, 0x50(%%rax)	\n\t"\
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
		/* Prefetch: base address and 1 index offset */\
		,	[__add0] "m" (Xadd0)\
		,	[__p1] "m" (Xp1)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r14","r15","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
	);\
	}

	#define SSE2_cmplx_carry_norm_pow2_nocheck2_2B(Xdata,XwtA,XwtB,XcyA,XcyB,Xbjmod_0,Xhalf_arr,Xn_minus_silp1,Xn_minus_sil,Xsinwt,Xsinwtm1,Xsse_bw,Xsse_nm1,Xsse_sw, Xadd0,Xp2,Xp3)\
	{\
	__asm__ volatile (\
	"movq	%[__add0],%%r14	\n\t"/* base address for 2 prefetches-from-main-data-array spread through this macro */\
	"movslq		%[__p2],%%r15	\n\t"	\
	"prefetcht0	(%%r14,%%r15,8)	\n\t"\
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
		"/*roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd		$0,%%xmm7,%%xmm7	*/\n\t"\
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
		"/*roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd		$0,%%xmm7,%%xmm7	*/\n\t"\
		"movaps		%%xmm3	,(%%rcx)	\n\t	movaps		%%xmm7	,(%%rdx)	\n\t"\
		"\n\t"\
		"movq	%[__data]	,%%rax		\n\t"\
		"mulpd	 0x80(%%rdi)	,%%xmm3	\n\t	mulpd	 0x80(%%rbx)	,%%xmm7	\n\t"\
		"subpd		%%xmm3	,	%%xmm1	\n\t	subpd		%%xmm7	,	%%xmm5	\n\t"\
		"mulpd		%%xmm2	,	%%xmm1	\n\t	mulpd		%%xmm6	,	%%xmm5	\n\t"\
		"movaps		%%xmm1	,0x20(%%rax)\n\t	movaps		%%xmm5	,0x60(%%rax)\n\t"\
		"\n\t"\
		"movq	%[__sse_bw]	,%%rax		\n\t"\
		"movq	%[__sse_nm1]	,%%rbx	\n\t"\
		"paddd		(%%rax)	,%%xmm0		\n\t"\
		"pand		(%%rbx)	,%%xmm0		\n\t"\
	"/**********************************************/\n\t"\
	"/*          Imaginary parts               */\n\t"\
	"/**********************************************/\n\t"\
	"movslq		%[__p3],%%r15	\n\t"	\
	"prefetcht0	(%%r14,%%r15,8)	\n\t"\
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
		"addpd		%%xmm4	,%%xmm3		\n\t	addpd		%%xmm4	,%%xmm7		\n\t"\
		"subpd		%%xmm4	,%%xmm3		\n\t	subpd		%%xmm4	,%%xmm7		\n\t"\
		"/*roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd		$0,%%xmm7,%%xmm7	*/\n\t"\
		"\n\t"\
		"/* NO ROE HERE */		\n\t"\
		"\n\t"\
		"movq	%%rsi,	%%rdi	\n\t	movq	%%rsi,	%%rbx	\n\t"\
		"shrq	$20,	%%rdi	\n\t	shrq	$22,	%%rbx	\n\t"\
		"andq	$0x0000000000000030	,%%rdi		\n\t	andq	$0x0000000000000030	,%%rbx		\n\t"\
		"addq	%%rax	,%%rdi			\n\t	addq	%%rax	,%%rbx			\n\t"\
		"movaps	%%xmm3	,%%xmm1			\n\t	movaps	%%xmm7	,%%xmm5			\n\t"\
		"mulpd	 0xc0(%%rdi)	,%%xmm3	\n\t	mulpd	 0xc0(%%rbx)	,%%xmm7	\n\t"\
		"addpd		%%xmm4	,%%xmm3		\n\t	addpd		%%xmm4	,%%xmm7		\n\t"\
		"subpd		%%xmm4	,%%xmm3		\n\t	subpd		%%xmm4	,%%xmm7		\n\t"\
		"/*roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd		$0,%%xmm7,%%xmm7	*/\n\t"\
		"movaps	%%xmm3	,(%%rcx)		\n\t	movaps	%%xmm7	,(%%rdx)		\n\t"\
		"\n\t"\
		"movq	%[__data]	,%%rax		\n\t"\
		"mulpd	 0x80(%%rdi)	,%%xmm3	\n\t	mulpd	 0x80(%%rbx)	,%%xmm7	\n\t"\
		"subpd	%%xmm3	,%%xmm1			\n\t	subpd	%%xmm7	,%%xmm5			\n\t"\
		"mulpd	%%xmm2	,%%xmm1			\n\t	mulpd	%%xmm6	,%%xmm5			\n\t"\
		"movaps	%%xmm1	, 0x30(%%rax)	\n\t	movaps	%%xmm5	, 0x70(%%rax)	\n\t"\
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
		/* Prefetch: base address and 2 index offsets */\
		,	[__add0] "m" (Xadd0)\
		,	[__p1] "m" (Xp1)\
		,	[__p2] "m" (Xp2)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r14","r15","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
	);\
	}


	/***************************************************************************************************************************************************/
	/********* Non-power-of-2-FFT versions of SSE2_cmplx_carry_norm_pow2_errcheck0_2B,1_2B,2_2B (only give sans-error-check version of latter 2: *******/
	/***************************************************************************************************************************************************/

#ifdef USE_AVX	// Our initial "AVX version" is simply the SSE2-based macros, applied to the AVX data layout
  #if 0
	#define SSE2_cmplx_carry_norm_errcheck0_2B(Xdata,XwtA,XwtB,XwtC,XcyA,XcyB,Xbjmod_0,Xhalf_arr,Xi,Xn_minus_silp1,Xn_minus_sil,Xsign_mask,Xsinwt,Xsinwtm1,Xsse_bw,Xsse_n,Xsse_sw)\
	{\
	__asm__ volatile (\
	"/***************Unpack the data:*************************/\n\t"\
		"movq	%[__data]	,%%rax	\n\t"\
		"movaps		    (%%rax)	,%%xmm1	\n\t	movaps		0x80(%%rax)	,%%xmm5	\n\t"\
		"movaps		    (%%rax)	,%%xmm2	\n\t	movaps		0x80(%%rax)	,%%xmm6	\n\t"\
		"unpcklpd	0x40(%%rax)	,%%xmm1	\n\t	unpcklpd	0xc0(%%rax)	,%%xmm5	\n\t"\
		"unpckhpd	0x40(%%rax)	,%%xmm2	\n\t	unpckhpd	0xc0(%%rax)	,%%xmm6	\n\t"\
		"movaps		%%xmm2, 0x40(%%rax)	\n\t	movaps		%%xmm6, 0xc0(%%rax)	\n\t"\
		"\n\t"\
		"movaps		0x20(%%rax)	,%%xmm2	\n\t	movaps		0xa0(%%rax)	,%%xmm6	\n\t"\
		"movaps		0x20(%%rax)	,%%xmm3	\n\t	movaps		0xa0(%%rax)	,%%xmm7	\n\t"\
		"unpcklpd	0x60(%%rax)	,%%xmm2	\n\t	unpcklpd	0xe0(%%rax)	,%%xmm6	\n\t"\
		"unpckhpd	0x60(%%rax)	,%%xmm3	\n\t	unpckhpd	0xe0(%%rax)	,%%xmm7	\n\t"\
		"movaps		%%xmm2, 0x20(%%rax)	\n\t	movaps		%%xmm6, 0xa0(%%rax)	\n\t"\
		"movaps		%%xmm3, 0x60(%%rax)	\n\t	movaps		%%xmm7, 0xe0(%%rax)	\n\t"\
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
		"maxpd		-0x40(%%rax),%%xmm1	\n\t"\
		"movaps		%%xmm1,-0x40(%%rax)	\n\t"\
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
		"movq	%[__data]	,%%rax		\n\t"\
		"mulpd	 0x80(%%rdi)	,%%xmm3	\n\t	mulpd	 0x80(%%rbx)	,%%xmm7	\n\t"\
		"subpd		%%xmm3	,	%%xmm1	\n\t	subpd		%%xmm7	,	%%xmm5	\n\t"\
		"mulpd		%%xmm2	,	%%xmm1	\n\t	mulpd		%%xmm6	,	%%xmm5	\n\t"\
		"movaps		%%xmm1	,    (%%rax)\n\t	movaps		%%xmm5	,0x80(%%rax)\n\t"\
		"\n\t"\
		"movq	%[__sse_n]	,%%rbx		\n\t"\
		"movaps		(%%rbx)	,%%xmm2		\n\t"\
		"movq	%[__sse_bw]	,%%rax		\n\t"\
		"paddd		(%%rax)	,%%xmm0		\n\t"\
		"movaps		%%xmm0	,%%xmm1		\n\t"\
		"pcmpgtd	%%xmm2	,%%xmm1		\n\t"\
		"pand		%%xmm2	,%%xmm1		\n\t"\
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
		"movaps	 0x20(%%rax)	,%%xmm1	\n\t"\
		"movaps	 0xa0(%%rax)	,%%xmm5	\n\t"\
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
		"maxpd	-0x40(%%rax)	,%%xmm1	\n\t"\
		"movaps	%%xmm1	,-0x40(%%rax)	\n\t"\
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
		"movaps	%%xmm1	, 0x20(%%rax)	\n\t	movaps	%%xmm5	, 0xa0(%%rax)	\n\t"\
		"\n\t"\
		"movq	%[__sse_n]	,%%rbx		\n\t"\
		"movaps		(%%rbx)	,%%xmm2		\n\t"\
		"movq	%[__sse_bw]	,%%rax		\n\t"\
		"paddd		(%%rax)	,%%xmm0		\n\t"\
		"movaps		%%xmm0	,%%xmm1		\n\t"\
		"pcmpgtd	%%xmm2	,%%xmm1		\n\t"\
		"pand		%%xmm2	,%%xmm1		\n\t"\
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
		"movaps		    (%%rax)	,%%xmm1	\n\t	movaps		0x80(%%rax)	,%%xmm5	\n\t"\
		"movaps		    (%%rax)	,%%xmm2	\n\t	movaps		0x80(%%rax)	,%%xmm6	\n\t"\
		"unpcklpd	0x40(%%rax)	,%%xmm1	\n\t	unpcklpd	0xc0(%%rax)	,%%xmm5	\n\t"\
		"unpckhpd	0x40(%%rax)	,%%xmm2	\n\t	unpckhpd	0xc0(%%rax)	,%%xmm6	\n\t"\
		"movaps		%%xmm2, 0x40(%%rax)	\n\t	movaps		%%xmm6, 0xc0(%%rax)	\n\t"\
		"\n\t"\
		"movaps		0x20(%%rax)	,%%xmm2	\n\t	movaps		0xa0(%%rax)	,%%xmm6	\n\t"\
		"movaps		0x20(%%rax)	,%%xmm3	\n\t	movaps		0xa0(%%rax)	,%%xmm7	\n\t"\
		"unpcklpd	0x60(%%rax)	,%%xmm2	\n\t	unpcklpd	0xe0(%%rax)	,%%xmm6	\n\t"\
		"unpckhpd	0x60(%%rax)	,%%xmm3	\n\t	unpckhpd	0xe0(%%rax)	,%%xmm7	\n\t"\
		"movaps		%%xmm2, 0x20(%%rax)	\n\t	movaps		%%xmm6, 0xa0(%%rax)	\n\t"\
		"movaps		%%xmm3, 0x60(%%rax)	\n\t	movaps		%%xmm7, 0xe0(%%rax)	\n\t"\
	"/**********************************************/\n\t"\
	"/*          Real      parts                   */\n\t"\
	"/**********************************************/\n\t"\
		"movq	%[__bjmod_0],	%%rax	\n\t"\
		"movaps		(%%rax)	,	%%xmm0	\n\t"\
		"movq	%[__sse_sw]	,	%%rbx	\n\t"\
		"movaps		(%%rbx)	,	%%xmm7	\n\t"\
		"psubd		%%xmm0	,	%%xmm7	\n\t"\
		"movmskps	%%xmm7	,	%%rsi	\n\t"\
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
		"maxpd	-0x40(%%rax)	,%%xmm1	\n\t"\
		"movaps	%%xmm1	,-0x40(%%rax)	\n\t"\
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
		"movq	%[__data]	,%%rax		\n\t"\
		"mulpd	 0x80(%%rdi)	,%%xmm3	\n\t	mulpd	 0x80(%%rbx)	,%%xmm7	\n\t"\
		"subpd		%%xmm3	,	%%xmm1	\n\t	subpd		%%xmm7	,	%%xmm5	\n\t"\
		"mulpd		%%xmm2	,	%%xmm1	\n\t	mulpd		%%xmm6	,	%%xmm5	\n\t"\
		"movaps		%%xmm1	,    (%%rax)\n\t	movaps		%%xmm5	,0x80(%%rax)\n\t"\
		"\n\t"\
		"movq	%[__sse_n]	,%%rbx		\n\t"\
		"movaps		(%%rbx)	,%%xmm2		\n\t"\
		"movq	%[__sse_bw]	,%%rax		\n\t"\
		"paddd		(%%rax)	,%%xmm0		\n\t"\
		"movaps		%%xmm0	,%%xmm1		\n\t"\
		"pcmpgtd	%%xmm2	,%%xmm1		\n\t"\
		"pand		%%xmm2	,%%xmm1		\n\t"\
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
		"movaps	 0x20(%%rax)	,%%xmm1	\n\t"\
		"movaps	 0xa0(%%rax)	,%%xmm5	\n\t"\
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
		"maxpd	-0x40(%%rax)	,%%xmm1	\n\t"\
		"movaps	%%xmm1	,-0x40(%%rax)	\n\t"\
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
		"movaps	%%xmm1	, 0x20(%%rax)	\n\t	movaps	%%xmm5	, 0xa0(%%rax)	\n\t"\
		"\n\t"\
		"movq	%[__sse_n]	,%%rbx		\n\t"\
		"movaps		(%%rbx)	,%%xmm2		\n\t"\
		"movq	%[__sse_bw]	,%%rax		\n\t"\
		"paddd		(%%rax)	,%%xmm0		\n\t"\
		"movaps		%%xmm0	,%%xmm1		\n\t"\
		"pcmpgtd	%%xmm2	,%%xmm1		\n\t"\
		"pand		%%xmm2	,%%xmm1		\n\t"\
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
		"movaps	 0x40(%%rax),%%xmm1		\n\t	movaps		 0xc0(%%rax),%%xmm5\n\t"\
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
		"maxpd		-0x40(%%rax),%%xmm1	\n\t"\
		"movaps		%%xmm1,-0x40(%%rax)	\n\t"\
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
		"movq	%[__data]	,%%rax		\n\t"\
		"mulpd	 0x80(%%rdi)	,%%xmm3	\n\t	mulpd	 0x80(%%rbx)	,%%xmm7	\n\t"\
		"subpd		%%xmm3	,	%%xmm1	\n\t	subpd		%%xmm7	,	%%xmm5	\n\t"\
		"mulpd		%%xmm2	,	%%xmm1	\n\t	mulpd		%%xmm6	,	%%xmm5	\n\t"\
		"movaps		%%xmm1	,0x40(%%rax)\n\t	movaps		%%xmm5	,0xc0(%%rax)\n\t"\
		"\n\t"\
		"movq	%[__sse_n]	,%%rbx		\n\t"\
		"movaps		(%%rbx)	,%%xmm2		\n\t"\
		"movq	%[__sse_bw]	,%%rax		\n\t"\
		"paddd		(%%rax)	,%%xmm0		\n\t"\
		"movaps		%%xmm0	,%%xmm1		\n\t"\
		"pcmpgtd	%%xmm2	,%%xmm1		\n\t"\
		"pand		%%xmm2	,%%xmm1		\n\t"\
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
		"movaps	 0x60(%%rax)	,%%xmm1	\n\t"\
		"movaps	 0xe0(%%rax)	,%%xmm5	\n\t"\
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
		"maxpd	-0x40(%%rax)	,%%xmm1	\n\t"\
		"movaps	%%xmm1	,-0x40(%%rax)	\n\t"\
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
		"movaps	%%xmm1	, 0x60(%%rax)	\n\t	movaps	%%xmm5	, 0xe0(%%rax)	\n\t"\
		"\n\t"\
		"movq	%[__sse_n]	,%%rbx		\n\t"\
		"movaps		(%%rbx)	,%%xmm2		\n\t"\
		"movq	%[__sse_bw]	,%%rax		\n\t"\
		"paddd		(%%rax)	,%%xmm0		\n\t"\
		"movaps		%%xmm0	,%%xmm1		\n\t"\
		"pcmpgtd	%%xmm2	,%%xmm1		\n\t"\
		"pand		%%xmm2	,%%xmm1		\n\t"\
		"psubd		%%xmm1	,%%xmm0		\n\t"\
		"movq	%[__bjmod_0],%%rcx		\n\t"\
		"movaps	%%xmm0,(%%rcx)			\n\t"\
	"/**********************************************/\n\t"\
	"/*          Repack the data:              */\n\t"\
	"/**********************************************/\n\t"\
		"mov	%[__data],%%rax			\n\t"\
		"movaps		0x20(%%rax)	,%%xmm1	\n\t	movaps		0xa0(%%rax)	,%%xmm5	\n\t"\
		"movaps		    (%%rax)	,%%xmm0	\n\t	movaps		0x80(%%rax)	,%%xmm4	\n\t"\
		"movaps		%%xmm1		,%%xmm3	\n\t	movaps		%%xmm5		,%%xmm7	\n\t"\
		"movaps		%%xmm0		,%%xmm2	\n\t	movaps		%%xmm4		,%%xmm6	\n\t"\
		"unpckhpd	0x60(%%rax)	,%%xmm3	\n\t	unpckhpd	0xe0(%%rax)	,%%xmm7	\n\t"\
		"unpcklpd	0x60(%%rax)	,%%xmm1	\n\t	unpcklpd	0xe0(%%rax)	,%%xmm5	\n\t"\
		"movaps		%%xmm3,0x60(%%rax)	\n\t	movaps		%%xmm7,0xe0(%%rax)	\n\t"\
		"unpckhpd	0x40(%%rax)	,%%xmm2	\n\t	unpckhpd	0xc0(%%rax)	,%%xmm6	\n\t"\
		"unpcklpd	0x40(%%rax)	,%%xmm0	\n\t	unpcklpd	0xc0(%%rax)	,%%xmm4	\n\t"\
		"movaps		%%xmm2,0x40(%%rax)	\n\t	movaps		%%xmm6,0xc0(%%rax)	\n\t"\
		"movaps		%%xmm1,0x20(%%rax)	\n\t	movaps		%%xmm5,0xa0(%%rax)	\n\t"\
		"movaps		%%xmm0,    (%%rax)	\n\t	movaps		%%xmm4,0x80(%%rax)	\n\t"\
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
  #endif	// if 0

#else

	#define SSE2_cmplx_carry_norm_errcheck0_2B(Xdata,XwtA,XwtB,XwtC,XcyA,XcyB,Xbjmod_0,Xhalf_arr,Xi,Xn_minus_silp1,Xn_minus_sil,Xsign_mask,Xsinwt,Xsinwtm1,Xsse_bw,Xsse_n,Xsse_sw, Xadd0,Xp1)\
	{\
	__asm__ volatile (\
	"movq	%[__add0],%%r14	\n\t"/* base address for 2 prefetches-from-main-data-array spread through this macro */\
	"prefetcht0	(%%r14)		\n\t"\
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
		"mulpd	 	0x100(%%rax),%%xmm2	\n\t	mulpd	 	0x100(%%rax),%%xmm6	\n\t"/* wt   =wtA*wtl */\
		"mulpd	 	0x110(%%rax),%%xmm3	\n\t	mulpd	 	0x110(%%rax),%%xmm7	\n\t"/* wtinv=wtB*wtn */\
		"mulpd	 	     (%%rdi),%%xmm2	\n\t	mulpd	 	     (%%rbx),%%xmm6	\n\t"/* wt   =wt   *one_half[m01] */\
		"mulpd	 	0x040(%%rdx),%%xmm3	\n\t	mulpd	 	0x040(%%rcx),%%xmm7	\n\t"/* wtinv=wtinv*one_half[4+m23] */\
		"\n\t"\
		"movq	%[__cyA]	,%%rcx		\n\t	movq	%[__cyB]	,%%rdx		\n\t"\
		"mulpd		%%xmm3	,%%xmm1		\n\t	mulpd		%%xmm7	,%%xmm5		\n\t"\
		"addpd		(%%rcx)	,%%xmm1		\n\t	addpd		(%%rdx)	,%%xmm5		\n\t"\
		"movaps		%%xmm1	,%%xmm3		\n\t	movaps		%%xmm5	,%%xmm7		\n\t"\
		"addpd		%%xmm4	,%%xmm3		\n\t	addpd		%%xmm4	,%%xmm7		\n\t"\
		"subpd		%%xmm4	,%%xmm3		\n\t	subpd		%%xmm4	,%%xmm7		\n\t"\
		"/*roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd		$0,%%xmm7,%%xmm7	*/\n\t"\
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
		"/*roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd		$0,%%xmm7,%%xmm7	*/\n\t"\
		"movaps		%%xmm3	,(%%rcx)	\n\t	movaps		%%xmm7	,(%%rdx)	\n\t"\
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
		"pand		%%xmm2	,%%xmm1		\n\t"\
		"psubd		%%xmm1	,%%xmm0		\n\t"\
	"/**********************************************/\n\t"\
	"/*          Imaginary parts               */\n\t"\
	"/**********************************************/\n\t"\
	"movslq		%[__p1],%%r15	\n\t"	\
	"prefetcht0	(%%r14,%%r15,8)	\n\t"\
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
		"mulpd	 0x120(%%rax)	,%%xmm2	\n\t	mulpd	 0x120(%%rax)	,%%xmm6	\n\t"/* wt   =wtA*wtlp1 */\
		"mulpd	 0x130(%%rax)	,%%xmm3	\n\t	mulpd	 0x130(%%rax)	,%%xmm7	\n\t"/* wtinv=wtC*wtnm1 */\
		"mulpd	      (%%rdi)	,%%xmm2	\n\t	mulpd	      (%%rbx)	,%%xmm6	\n\t"/* wt   =wt   *one_half[m01] */\
		"mulpd	 0x040(%%rdx)	,%%xmm3	\n\t	mulpd	 0x040(%%rcx)	,%%xmm7	\n\t"/* wtinv=wtinv*one_half[4+m23] */\
		"\n\t"\
		"movq	%[__cyA]	,%%rcx		\n\t	movq	%[__cyB]	,%%rdx		\n\t"\
		"mulpd	%%xmm3	,%%xmm1			\n\t	mulpd	%%xmm7	,%%xmm5			\n\t"\
		"addpd	(%%rcx)	,%%xmm1			\n\t	addpd	(%%rdx)	,%%xmm5			\n\t"\
		"movaps	%%xmm1	,%%xmm3			\n\t	movaps	%%xmm5	,%%xmm7			\n\t"\
		"addpd		%%xmm4	,%%xmm3		\n\t	addpd		%%xmm4	,%%xmm7		\n\t"\
		"subpd		%%xmm4	,%%xmm3		\n\t	subpd		%%xmm4	,%%xmm7		\n\t"\
		"/*roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd		$0,%%xmm7,%%xmm7	*/\n\t"\
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
		"addpd		%%xmm4	,%%xmm3		\n\t	addpd		%%xmm4	,%%xmm7		\n\t"\
		"subpd		%%xmm4	,%%xmm3		\n\t	subpd		%%xmm4	,%%xmm7		\n\t"\
		"/*roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd		$0,%%xmm7,%%xmm7	*/\n\t"\
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
		"pand		%%xmm2	,%%xmm1		\n\t"\
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
		/* Prefetch: base address and 1 index offset */\
		,	[__add0] "m" (Xadd0)\
		,	[__p1] "m" (Xp1)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r14","r15","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
	);\
	}

	#define SSE2_cmplx_carry_norm_errcheck1_2B(Xdata,XwtA,XwtB,XwtC,XcyA,XcyB,Xbjmod_0,Xhalf_arr,Xn_minus_silp1,Xn_minus_sil,Xsign_mask,Xsinwt,Xsinwtm1,Xsse_bw,Xsse_n,Xsse_sw, Xadd0,Xp1)\
	{\
	__asm__ volatile (\
	"movq	%[__add0],%%r14	\n\t"/* base address for 2 prefetches-from-main-data-array spread through this macro */\
	"prefetcht0	(%%r14)		\n\t"\
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
		"/*roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd		$0,%%xmm7,%%xmm7	*/\n\t"\
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
		"addpd		%%xmm4	,%%xmm3		\n\t	addpd		%%xmm4	,%%xmm7		\n\t"\
		"subpd		%%xmm4	,%%xmm3		\n\t	subpd		%%xmm4	,%%xmm7		\n\t"\
		"/*roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd		$0,%%xmm7,%%xmm7	*/\n\t"\
		"movaps		%%xmm3	,(%%rcx)	\n\t	movaps		%%xmm7	,(%%rdx)	\n\t"\
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
		"pand		%%xmm2	,%%xmm1		\n\t"\
		"psubd		%%xmm1	,%%xmm0		\n\t"\
	"/**********************************************/\n\t"\
	"/*          Imaginary parts               */\n\t"\
	"/**********************************************/\n\t"\
	"movslq		%[__p1],%%r15	\n\t"	\
	"prefetcht0	(%%r14,%%r15,8)	\n\t"\
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
		"addpd		%%xmm4	,%%xmm3		\n\t	addpd		%%xmm4	,%%xmm7		\n\t"\
		"subpd		%%xmm4	,%%xmm3		\n\t	subpd		%%xmm4	,%%xmm7		\n\t"\
		"/*roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd		$0,%%xmm7,%%xmm7	*/\n\t"\
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
		"addpd		%%xmm4	,%%xmm3		\n\t	addpd		%%xmm4	,%%xmm7		\n\t"\
		"subpd		%%xmm4	,%%xmm3		\n\t	subpd		%%xmm4	,%%xmm7		\n\t"\
		"/*roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd		$0,%%xmm7,%%xmm7	*/\n\t"\
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
		"pand		%%xmm2	,%%xmm1		\n\t"\
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
		/* Prefetch: base address and 1 index offset */\
		,	[__add0] "m" (Xadd0)\
		,	[__p1] "m" (Xp1)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r14","r15","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
	);\
	}

	#define SSE2_cmplx_carry_norm_errcheck2_2B(Xdata,XwtA,XwtB,XcyA,XcyB,Xbjmod_0,Xhalf_arr,Xn_minus_silp1,Xn_minus_sil,Xsign_mask,Xsinwt,Xsinwtm1,Xsse_bw,Xsse_n,Xsse_sw, Xadd0,Xp2,Xp3)\
	{\
	__asm__ volatile (\
	"movq	%[__add0],%%r14	\n\t"/* base address for 2 prefetches-from-main-data-array spread through this macro */\
	"movslq		%[__p2],%%r15	\n\t"	\
	"prefetcht0	(%%r14,%%r15,8)	\n\t"\
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
		"/*roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd		$0,%%xmm7,%%xmm7	*/\n\t"\
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
		"/*roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd		$0,%%xmm7,%%xmm7	*/\n\t"\
		"movaps		%%xmm3	,(%%rcx)	\n\t	movaps		%%xmm7	,(%%rdx)	\n\t"\
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
		"pand		%%xmm2	,%%xmm1		\n\t"\
		"psubd		%%xmm1	,%%xmm0		\n\t"\
	"/**********************************************/\n\t"\
	"/*          Imaginary parts               */\n\t"\
	"/**********************************************/\n\t"\
	"movslq		%[__p3],%%r15	\n\t"	\
	"prefetcht0	(%%r14,%%r15,8)	\n\t"\
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
		"addpd		%%xmm4	,%%xmm3		\n\t	addpd		%%xmm4	,%%xmm7		\n\t"\
		"subpd		%%xmm4	,%%xmm3		\n\t	subpd		%%xmm4	,%%xmm7		\n\t"\
		"/*roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd		$0,%%xmm7,%%xmm7	*/\n\t"\
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
		"addpd		%%xmm4	,%%xmm3		\n\t	addpd		%%xmm4	,%%xmm7		\n\t"\
		"subpd		%%xmm4	,%%xmm3		\n\t	subpd		%%xmm4	,%%xmm7		\n\t"\
		"/*roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd		$0,%%xmm7,%%xmm7	*/\n\t"\
		"movaps	%%xmm3	,(%%rcx)		\n\t	movaps	%%xmm7	,(%%rdx)		\n\t"\
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
		"pand		%%xmm2	,%%xmm1		\n\t"\
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
		/* Prefetch: base address and 2 index offsets */\
		,	[__add0] "m" (Xadd0)\
		,	[__p2] "m" (Xp2)\
		,	[__p3] "m" (Xp3)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r14","r15","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
	);\
	}

#endif	// #ifdef USE_AVX

	/******************************************************************************************************************************************************************/
	/********** No-ROE-Check versions of the latter 2 of the above 3 macros - the first is called too infrequently to bother with a special non-ROE version: **********/
	/******************************************************************************************************************************************************************/

	#define SSE2_cmplx_carry_norm_nocheck1_2B(Xdata,XwtA,XwtB,XwtC,XcyA,XcyB,Xbjmod_0,Xhalf_arr,Xn_minus_silp1,Xn_minus_sil,Xsinwt,Xsinwtm1,Xsse_bw,Xsse_n,Xsse_sw, Xadd0,Xp1)\
	{\
	__asm__ volatile (\
	"movq	%[__add0],%%r14	\n\t"/* base address for 2 prefetches-from-main-data-array spread through this macro */\
	"prefetcht0	(%%r14)		\n\t"\
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
		"/*roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd		$0,%%xmm7,%%xmm7	*/\n\t"\
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
		"/*roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd		$0,%%xmm7,%%xmm7	*/\n\t"\
		"movaps		%%xmm3	,(%%rcx)	\n\t	movaps		%%xmm7	,(%%rdx)	\n\t"\
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
		"pand		%%xmm2	,%%xmm1		\n\t"\
		"psubd		%%xmm1	,%%xmm0		\n\t"\
	"/**********************************************/\n\t"\
	"/*          Imaginary parts               */\n\t"\
	"/**********************************************/\n\t"\
	"movslq		%[__p1],%%r15	\n\t"	\
	"prefetcht0	(%%r14,%%r15,8)	\n\t"\
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
		"addpd		%%xmm4	,%%xmm3		\n\t	addpd		%%xmm4	,%%xmm7		\n\t"\
		"subpd		%%xmm4	,%%xmm3		\n\t	subpd		%%xmm4	,%%xmm7		\n\t"\
		"/*roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd		$0,%%xmm7,%%xmm7	*/\n\t"\
		"\n\t"\
		"/* NO ROE HERE */				\n\t"\
		"\n\t"\
		"movq	%%rsi,	%%rdi	\n\t	movq	%%rsi,	%%rbx	\n\t"\
		"shrq	$20,	%%rdi	\n\t	shrq	$22,	%%rbx	\n\t"\
		"andq	$0x0000000000000030	,%%rdi		\n\t	andq	$0x0000000000000030	,%%rbx		\n\t"\
		"addq	%%rax	,%%rdi			\n\t	addq	%%rax	,%%rbx			\n\t"\
		"movaps	%%xmm3	,%%xmm1			\n\t	movaps	%%xmm7	,%%xmm5			\n\t"\
		"mulpd	 0xc0(%%rdi)	,%%xmm3	\n\t	mulpd	 0xc0(%%rbx)	,%%xmm7	\n\t"\
		"addpd		%%xmm4	,%%xmm3		\n\t	addpd		%%xmm4	,%%xmm7		\n\t"\
		"subpd		%%xmm4	,%%xmm3		\n\t	subpd		%%xmm4	,%%xmm7		\n\t"\
		"/*roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd		$0,%%xmm7,%%xmm7	*/\n\t"\
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
		"pand		%%xmm2	,%%xmm1		\n\t"\
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
		/* Prefetch: base address and 1 index offset */\
		,	[__add0] "m" (Xadd0)\
		,	[__p1] "m" (Xp1)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r14","r15","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
	);\
	}

	#define SSE2_cmplx_carry_norm_nocheck2_2B(Xdata,XwtA,XwtB,XcyA,XcyB,Xbjmod_0,Xhalf_arr,Xn_minus_silp1,Xn_minus_sil,Xsinwt,Xsinwtm1,Xsse_bw,Xsse_n,Xsse_sw, Xadd0,Xp2,Xp3)\
	{\
	__asm__ volatile (\
	"movq	%[__add0],%%r14	\n\t"/* base address for 2 prefetches-from-main-data-array spread through this macro */\
	"movslq		%[__p2],%%r15	\n\t"	\
	"prefetcht0	(%%r14,%%r15,8)	\n\t"\
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
		"/*roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd		$0,%%xmm7,%%xmm7	*/\n\t"\
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
		"/*roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd		$0,%%xmm7,%%xmm7	*/\n\t"\
		"movaps		%%xmm3	,(%%rcx)	\n\t	movaps		%%xmm7	,(%%rdx)	\n\t"\
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
		"pand		%%xmm2	,%%xmm1		\n\t"\
		"psubd		%%xmm1	,%%xmm0		\n\t"\
	"/**********************************************/\n\t"\
	"/*          Imaginary parts               */\n\t"\
	"/**********************************************/\n\t"\
	"movslq		%[__p3],%%r15	\n\t"	\
	"prefetcht0	(%%r14,%%r15,8)	\n\t"\
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
		"addpd		%%xmm4	,%%xmm3		\n\t	addpd		%%xmm4	,%%xmm7		\n\t"\
		"subpd		%%xmm4	,%%xmm3		\n\t	subpd		%%xmm4	,%%xmm7		\n\t"\
		"/*roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd		$0,%%xmm7,%%xmm7	*/\n\t"\
		"\n\t"\
		"/* NO ROE HERE */				\n\t"\
		"\n\t"\
		"movq	%%rsi,	%%rdi	\n\t	movq	%%rsi,	%%rbx	\n\t"\
		"shrq	$20,	%%rdi	\n\t	shrq	$22,	%%rbx	\n\t"\
		"andq	$0x0000000000000030	,%%rdi		\n\t	andq	$0x0000000000000030	,%%rbx		\n\t"\
		"addq	%%rax	,%%rdi			\n\t	addq	%%rax	,%%rbx			\n\t"\
		"movaps	%%xmm3	,%%xmm1			\n\t	movaps	%%xmm7	,%%xmm5			\n\t"\
		"mulpd	 0xc0(%%rdi)	,%%xmm3	\n\t	mulpd	 0xc0(%%rbx)	,%%xmm7	\n\t"\
		"addpd		%%xmm4	,%%xmm3		\n\t	addpd		%%xmm4	,%%xmm7		\n\t"\
		"subpd		%%xmm4	,%%xmm3		\n\t	subpd		%%xmm4	,%%xmm7		\n\t"\
		"/*roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd		$0,%%xmm7,%%xmm7	*/\n\t"\
		"movaps	%%xmm3	,(%%rcx)		\n\t	movaps	%%xmm7	,(%%rdx)		\n\t"\
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
		"pand		%%xmm2	,%%xmm1		\n\t"\
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
		/* Prefetch: base address and 2 index offsets */\
		,	[__add0] "m" (Xadd0)\
		,	[__p1] "m" (Xp1)\
		,	[__p2] "m" (Xp2)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r14","r15","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
	);\
	}

#endif	/* carry_gcc_h_included */

