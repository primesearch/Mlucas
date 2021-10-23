/*******************************************************************************
*                                                                              *
*   (C) 1997-2021 by Ernst W. Mayer.                                           *
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
#ifndef twopmodq100_h_included
#define twopmodq100_h_included

#include "masterdefs.h"
#include "types.h"

#ifdef __cplusplus
extern "C" {
#endif

#ifdef USE_AVX512	// Implies FMA3 support

  #ifndef ALL_FMA	// Default version uses ADDPD,SUBPD:

	#define SSE2_twopmodq100_modmul_q32(Xaq0,Xaqinv0,Xax0,Xbase,Xbinv,Xcrnd,Xpshift,Xj)\
	{\
	__asm__ volatile (\
		"movq	%[__base],%%rax		\n\t	vbroadcastsd	(%%rax),%%zmm31	\n\t"/* BASE */\
		"movq	%[__binv],%%rbx		\n\t	vbroadcastsd	(%%rbx),%%zmm30	\n\t"/* BINV */\
		"movq	%[__crnd],%%rcx		\n\t	vbroadcastsd	(%%rcx),%%zmm29	\n\t"/* CRND */\
		"vmulpd	%%zmm31,%%zmm29,%%zmm29			\n\t"/* CRND*BASE */\
		"movq	%[__ax0],%%rax					\n\t"/* &x */\
		"movq	%[__aq0],%%rbx					\n\t"\
		"movq	%[__aqinv0],%%rdx				\n\t"\
		/*** STREAM 0 ***/							/*** STREAM 1 ***/							/*** STREAM 2 ***/								/*** STREAM 3 ***/					\
		/* Inputs [a|c|e|g]lo0-2 enter in zmm[0-1,4-5,8-9,12-13], resp.: */\
	/************ hi:lo = SQR_LOHI(x) : **************/\
		"vmovaps	(%%rax),%%zmm1				\n\t	vmovaps	0x040(%%rax),%%zmm5				\n\t	vmovaps	0x080(%%rax),%%zmm9				\n\t	vmovaps	0x0c0(%%rax),%%zmm13			\n\t"/* load x.lo */\
		"vaddpd			%%zmm1,%%zmm1,%%zmm2	\n\t	vaddpd			%%zmm5,%%zmm5,%%zmm6	\n\t	vaddpd			%%zmm9,%%zmm9,%%zmm10	\n\t	vaddpd			%%zmm13,%%zmm13,%%zmm14	\n\t"/* 2*lo */\
		/* lo*lo: a = b = x.lo: */\
		"     vmulpd	%%zmm1,%%zmm1,%%zmm3	\n\t	     vmulpd	%%zmm5,%%zmm5,%%zmm7		\n\t	     vmulpd	%%zmm9,%%zmm9,%%zmm11		\n\t	     vmulpd	%%zmm13,%%zmm13,%%zmm15		\n\t"/* hi = fma(a,b, 0  ) - use VMULPD instead of VFMADD231PD, since addend = 0. */\
		"vmovaps	%%zmm3,%%zmm0				\n\t	vmovaps	%%zmm7,%%zmm4					\n\t	vmovaps	%%zmm11,%%zmm8					\n\t	vmovaps	%%zmm15,%%zmm12					\n\t"/* cpy hi into lo-reg */\
		"vfmsub231pd	%%zmm1,%%zmm1,%%zmm0	\n\t	vfmsub231pd	%%zmm5,%%zmm5,%%zmm4		\n\t	vfmsub231pd	%%zmm9,%%zmm9,%%zmm8		\n\t	vfmsub231pd	%%zmm13,%%zmm13,%%zmm12		\n\t"/* lo = fma(a,b, -hi) */\
		"vaddpd			%%zmm29,%%zmm3,%%zmm1	\n\t	vaddpd			%%zmm29,%%zmm7,%%zmm5	\n\t	vaddpd			%%zmm29,%%zmm11,%%zmm9	\n\t	vaddpd			%%zmm29,%%zmm15,%%zmm13	\n\t"\
		"vsubpd			%%zmm29,%%zmm1,%%zmm1	\n\t	vsubpd			%%zmm29,%%zmm5,%%zmm5	\n\t	vsubpd			%%zmm29,%%zmm9,%%zmm9	\n\t	vsubpd			%%zmm29,%%zmm13,%%zmm13	\n\t"/* hh = hi +- crnd50 to round-to-nearest-multiple-of-base */\
		"vsubpd			%%zmm1,%%zmm3,%%zmm3	\n\t	vsubpd			%%zmm5,%%zmm7,%%zmm7	\n\t	vsubpd			%%zmm9,%%zmm11,%%zmm11	\n\t	vsubpd			%%zmm13,%%zmm15,%%zmm15	\n\t"/* hi - hh gives backward-carry... */\
		"vaddpd			%%zmm3,%%zmm0,%%zmm0	\n\t	vaddpd			%%zmm7,%%zmm4,%%zmm4	\n\t	vaddpd			%%zmm11,%%zmm8,%%zmm8	\n\t	vaddpd			%%zmm15,%%zmm12,%%zmm12	\n\t"/* ...which we add to lo (prod[0]). */\
		"vmovaps	%%zmm0,     (%%rax)			\n\t	vmovaps	%%zmm4,0x040(%%rax)				\n\t	vmovaps	%%zmm8,0x080(%%rax)				\n\t	vmovaps	%%zmm12,0x0c0(%%rax)			\n\t"/* write prod[0]; this overwrites x.lo, but only need 2*x.lo for rest of squaring, that is in zmm2 */\
													/*** prod[0] in zmm0, hi*base in zmm1. ***/\
		/* 2*lo*hi: a = 2*x.lo, b = x.hi: */\
		"     vmulpd 0x100(%%rax),%%zmm2,%%zmm3	\n\t	     vmulpd 0x140(%%rax),%%zmm6,%%zmm7	\n\t	     vmulpd 0x180(%%rax),%%zmm10,%%zmm11\n\t	     vmulpd 0x1c0(%%rax),%%zmm14,%%zmm15\n\t"/* hi = fma(a,b, 0  ) - use VMULPD instead of VFMADD231PD, since addend = 0. */\
		"vmovaps	%%zmm3,%%zmm0				\n\t	vmovaps	%%zmm7,%%zmm4					\n\t	vmovaps	%%zmm11,%%zmm8					\n\t	vmovaps	%%zmm15,%%zmm12					\n\t"/* cpy hi into lo-reg */\
		"vfmsub231pd 0x100(%%rax),%%zmm2,%%zmm0	\n\t	vfmsub231pd 0x140(%%rax),%%zmm6,%%zmm4	\n\t	vfmsub231pd 0x180(%%rax),%%zmm10,%%zmm8	\n\t	vfmsub231pd 0x1c0(%%rax),%%zmm14,%%zmm12\n\t"/* lo = fma(a,b, -hi) */\
		"     vmulpd	%%zmm30,%%zmm3,%%zmm2	\n\t	     vmulpd	%%zmm30,%%zmm7,%%zmm6		\n\t	     vmulpd	%%zmm30,%%zmm11,%%zmm10		\n\t	     vmulpd	%%zmm30,%%zmm15,%%zmm14		\n\t"/* hi*binv */\
		"vrndscalepd $1,%%zmm2,%%zmm2			\n\t	vrndscalepd $1,%%zmm6,%%zmm6			\n\t	vrndscalepd $1,%%zmm10,%%zmm10			\n\t	vrndscalepd $1,%%zmm14,%%zmm14			\n\t"/* hh = floor(hi*binv); This part remains in hi */\
		"vfmadd231pd	%%zmm30,%%zmm1,%%zmm0	\n\t	vfmadd231pd	%%zmm30,%%zmm5,%%zmm4		\n\t	vfmadd231pd	%%zmm30,%%zmm9,%%zmm8		\n\t	vfmadd231pd	%%zmm30,%%zmm13,%%zmm12		\n\t"/* Add lo to hi-output of previous lo:hi pair, which also needs *= binv */\
		"vfnmadd231pd	%%zmm31,%%zmm2,%%zmm3	\n\t	vfnmadd231pd	%%zmm31,%%zmm6,%%zmm7	\n\t	vfnmadd231pd	%%zmm31,%%zmm10,%%zmm11	\n\t	vfnmadd231pd	%%zmm31,%%zmm14,%%zmm15	\n\t"/* cy = FMA(hh,-base, hi); backward-carry from hi into lo... */\
		"vaddpd			%%zmm3,%%zmm0,%%zmm1	\n\t	vaddpd			%%zmm7,%%zmm4,%%zmm5	\n\t	vaddpd			%%zmm11,%%zmm8,%%zmm9	\n\t	vaddpd			%%zmm15,%%zmm12,%%zmm13	\n\t"/* ...which we add to lo (prod[1]). */\
	/* MONT_MUL *requires* the h:l = x^2 results to be properly normalized in order to yield correct result;
	at this point prod[1] can be > base, thus must do a normalize-and-carry: */\
		"     vmulpd	%%zmm30,%%zmm1,%%zmm3	\n\t	     vmulpd	%%zmm30,%%zmm5,%%zmm7		\n\t	     vmulpd	%%zmm30,%%zmm9,%%zmm11		\n\t	     vmulpd	%%zmm30,%%zmm13,%%zmm15		\n\t"/* prod[1]*binv */\
		"vrndscalepd $1,%%zmm3,%%zmm3			\n\t	vrndscalepd $1,%%zmm7,%%zmm7			\n\t	vrndscalepd $1,%%zmm11,%%zmm11			\n\t	vrndscalepd $1,%%zmm15,%%zmm15			\n\t"/* cy = floor(prod[1]*binv) */\
		"vfnmadd231pd	%%zmm31,%%zmm3,%%zmm1	\n\t	vfnmadd231pd	%%zmm31,%%zmm7,%%zmm5	\n\t	vfnmadd231pd	%%zmm31,%%zmm11,%%zmm9	\n\t	vfnmadd231pd	%%zmm31,%%zmm15,%%zmm13	\n\t"/* prod[1] -= cy*base */\
		"vaddpd			%%zmm3,%%zmm2,%%zmm2	\n\t	vaddpd			%%zmm7,%%zmm6,%%zmm6	\n\t	vaddpd			%%zmm11,%%zmm10,%%zmm10	\n\t	vaddpd			%%zmm15,%%zmm14,%%zmm14	\n\t"/* prod[2] += cy      */\
		"     vmulpd	%%zmm31,%%zmm2,%%zmm2	\n\t	     vmulpd	%%zmm31,%%zmm6,%%zmm6		\n\t	     vmulpd	%%zmm31,%%zmm10,%%zmm10		\n\t	     vmulpd	%%zmm31,%%zmm14,%%zmm14		\n\t"/* rescale hh *= base */\
		"vmovaps	0x100(%%rax),%%zmm3			\n\t	vmovaps	0x140(%%rax),%%zmm7				\n\t	vmovaps	0x180(%%rax),%%zmm11			\n\t	vmovaps	0x1c0(%%rax),%%zmm15			\n\t"/* reload hi into zmm3... */\
		"vmovaps	%%zmm1,0x100(%%rax)			\n\t	vmovaps	%%zmm5,0x140(%%rax)				\n\t	vmovaps	%%zmm9,0x180(%%rax)				\n\t	vmovaps	%%zmm13,0x1c0(%%rax)			\n\t"/* ...before overwriting with prod[1] */\
													/*** prod[1] in zmm1, hi*base in zmm2, x_in.hi in zmm3. ***/\
		/* hi*hi: a = b = x.hi: */\
		"     vmulpd	%%zmm3,%%zmm3,%%zmm0	\n\t	     vmulpd	%%zmm7,%%zmm7,%%zmm4		\n\t	     vmulpd	%%zmm11,%%zmm11,%%zmm8		\n\t	     vmulpd	%%zmm15,%%zmm15,%%zmm12		\n\t"/* hi = fma(a,b, 0  ) - use VMULPD instead of VFMADD231PD, since addend = 0. */\
		"vmovaps	%%zmm0,%%zmm1				\n\t	vmovaps	%%zmm4,%%zmm5					\n\t	vmovaps	%%zmm8,%%zmm9					\n\t	vmovaps	%%zmm12,%%zmm13					\n\t"/* cpy hi into lo-reg */\
		"vfmsub231pd	%%zmm3,%%zmm3,%%zmm1	\n\t	vfmsub231pd	%%zmm7,%%zmm7,%%zmm5		\n\t	vfmsub231pd	%%zmm11,%%zmm11,%%zmm9		\n\t	vfmsub231pd	%%zmm15,%%zmm15,%%zmm13		\n\t"/* lo = fma(a,b, -hi) */\
		"vaddpd			%%zmm29,%%zmm0,%%zmm3	\n\t	vaddpd			%%zmm29,%%zmm4,%%zmm7	\n\t	vaddpd			%%zmm29,%%zmm8,%%zmm11	\n\t	vaddpd			%%zmm29,%%zmm12,%%zmm15	\n\t"\
		"vsubpd			%%zmm29,%%zmm3,%%zmm3	\n\t	vsubpd			%%zmm29,%%zmm7,%%zmm7	\n\t	vsubpd			%%zmm29,%%zmm11,%%zmm11	\n\t	vsubpd			%%zmm29,%%zmm15,%%zmm15	\n\t"/* hh = hi +- crnd50 to round-to-nearest-multiple-of-base */\
		"vfmadd231pd	%%zmm30,%%zmm2,%%zmm1	\n\t	vfmadd231pd	%%zmm30,%%zmm6,%%zmm5		\n\t	vfmadd231pd	%%zmm30,%%zmm10,%%zmm9		\n\t	vfmadd231pd	%%zmm30,%%zmm14,%%zmm13		\n\t"/* Add lo to hi-output of previous lo:hi pair, which also needs *= binv */\
		"vmovaps	%%zmm0,0x400(%%rax)			\n\t	vmovaps	%%zmm4,0x440(%%rax)				\n\t	vmovaps	%%zmm8,0x480(%%rax)				\n\t	vmovaps	%%zmm12,0x4c0(%%rax)			\n\t"/* write hi-approximant to full 2-word product into scratch slot set aside for it */\
		"vsubpd			%%zmm3,%%zmm0,%%zmm0	\n\t	vsubpd			%%zmm7,%%zmm4,%%zmm4	\n\t	vsubpd			%%zmm11,%%zmm8,%%zmm8	\n\t	vsubpd			%%zmm15,%%zmm12,%%zmm12	\n\t"/* hi - hh gives backward-carry... */\
		"vaddpd			%%zmm0,%%zmm1,%%zmm2	\n\t	vaddpd			%%zmm4,%%zmm5,%%zmm6	\n\t	vaddpd			%%zmm8,%%zmm9,%%zmm10	\n\t	vaddpd			%%zmm12,%%zmm13,%%zmm14	\n\t"/* ...which we add to lo (prod[2]). */\
		"     vmulpd	%%zmm3,%%zmm30,%%zmm3	\n\t	     vmulpd	%%zmm7,%%zmm30,%%zmm7		\n\t	     vmulpd	%%zmm11,%%zmm30,%%zmm11		\n\t	     vmulpd	%%zmm15,%%zmm30,%%zmm15		\n\t"/* prod[3] = hh*binv */\
		"vmovaps	%%zmm2,0x200(%%rax)			\n\t	vmovaps	%%zmm6,0x240(%%rax)				\n\t	vmovaps	%%zmm10,0x280(%%rax)			\n\t	vmovaps	%%zmm14,0x2c0(%%rax)			\n\t"/* write prod[2] */\
		"vmovaps	%%zmm3,0x300(%%rax)			\n\t	vmovaps	%%zmm7,0x340(%%rax)				\n\t	vmovaps	%%zmm11,0x380(%%rax)			\n\t	vmovaps	%%zmm15,0x3c0(%%rax)			\n\t"/* write prod[3] */\
													/*** prod[2,3] in zmm2,3. ***/\
	/********************* lo = MULL(qinv,lo) : ***************************************************************************
		Bitranges (B = base):													Subproducts:
											|----------- 2B-1: 0 -----------|	qinv.lo * lo.lo
							|----------- 3B-1: B -----------|					qinv.hi * lo.lo, qinv.lo * lo.hi
			|----------- 4B-1:2B -----------|									qinv.hi * lo.hi
											|----------- OUTPUT ------------|
		Thus only need full 2B-wide (qinv.lo * lo.lo) subproduct and low halves of the (qinv.hi * lo.lo), (qinv.lo * lo.hi)
		subproducts - would that AVX-512F supported VPMULLQ! [Alas that needs AVX-512DQ, i.e. Cannonlake and beyond.]
	*********************************************************************************************************************/\
		"vmovaps	(%%rax),%%zmm1				\n\t	vmovaps	0x040(%%rax),%%zmm5				\n\t	vmovaps	0x080(%%rax),%%zmm9				\n\t	vmovaps	0x0c0(%%rax),%%zmm13			\n\t"/* load lo.lo and save copy in zmm2 */\
		"vmovaps	%%zmm1, %%zmm2				\n\t	vmovaps	%%zmm5, %%zmm6					\n\t	vmovaps	%%zmm9, %%zmm10					\n\t	vmovaps	%%zmm13, %%zmm14				\n\t"\
		/* qinv.lo*lo.lo: a = qinv.lo; b = lo.lo: */\
		"     vmulpd	(%%rdx),%%zmm1,%%zmm3	\n\t	     vmulpd	0x040(%%rdx),%%zmm5,%%zmm7	\n\t	     vmulpd	0x080(%%rdx),%%zmm9,%%zmm11	\n\t	     vmulpd	0x0c0(%%rdx),%%zmm13,%%zmm15\n\t"/* hi = fma(a,b, 0  ) - use VMULPD instead of VFMADD231PD, since addend = 0. */\
		"vmovaps	%%zmm3,%%zmm0				\n\t	vmovaps	%%zmm7,%%zmm4					\n\t	vmovaps	%%zmm11,%%zmm8					\n\t	vmovaps	%%zmm15,%%zmm12					\n\t"/* cpy hi into lo-reg */\
		"vfmsub231pd	(%%rdx),%%zmm1,%%zmm0	\n\t	vfmsub231pd	0x040(%%rdx),%%zmm5,%%zmm4	\n\t	vfmsub231pd	0x080(%%rdx),%%zmm9,%%zmm8	\n\t	vfmsub231pd	0x0c0(%%rdx),%%zmm13,%%zmm12\n\t"/* lo = fma(a,b, -hi) */\
		"vaddpd			%%zmm29,%%zmm3,%%zmm1	\n\t	vaddpd			%%zmm29,%%zmm7,%%zmm5	\n\t	vaddpd			%%zmm29,%%zmm11,%%zmm9	\n\t	vaddpd			%%zmm29,%%zmm15,%%zmm13	\n\t"\
		"vsubpd			%%zmm29,%%zmm1,%%zmm1	\n\t	vsubpd			%%zmm29,%%zmm5,%%zmm5	\n\t	vsubpd			%%zmm29,%%zmm9,%%zmm9	\n\t	vsubpd			%%zmm29,%%zmm13,%%zmm13	\n\t"/* hh = hi +- crnd50 to round-to-nearest-multiple-of-base */\
		"vsubpd			%%zmm1,%%zmm3,%%zmm3	\n\t	vsubpd			%%zmm5,%%zmm7,%%zmm7	\n\t	vsubpd			%%zmm9,%%zmm11,%%zmm11	\n\t	vsubpd			%%zmm13,%%zmm15,%%zmm15	\n\t"/* hi - hh gives backward-carry... */\
		"vaddpd			%%zmm3,%%zmm0,%%zmm0	\n\t	vaddpd			%%zmm7,%%zmm4,%%zmm4	\n\t	vaddpd			%%zmm11,%%zmm8,%%zmm8	\n\t	vaddpd			%%zmm15,%%zmm12,%%zmm12	\n\t"/* ...which we add to lo (prod[0]). */\
		"vmovaps	%%zmm0,     (%%rax)			\n\t	vmovaps	%%zmm4,0x040(%%rax)				\n\t	vmovaps	%%zmm8,0x080(%%rax)				\n\t	vmovaps	%%zmm12,0x0c0(%%rax)			\n\t"/* write prod[0]; this overwrites lo.lo, but copy of that saved in zmm2 */\
													/*** hi*base in zmm1, lo.lo in zmm2 ***/\
		/* qinv.hi*lo.lo: a = qinv.hi, b = lo.lo: */\
		"     vmulpd 0x100(%%rdx),%%zmm2,%%zmm3	\n\t	     vmulpd 0x140(%%rdx),%%zmm6,%%zmm7	\n\t	     vmulpd 0x180(%%rdx),%%zmm10,%%zmm11\n\t	     vmulpd 0x1c0(%%rdx),%%zmm14,%%zmm15\n\t"/* hi = fma(a,b, 0  ) - use VMULPD instead of VFMADD231PD, since addend = 0. */\
		"vmovaps	%%zmm3,%%zmm0				\n\t	vmovaps	%%zmm7,%%zmm4					\n\t	vmovaps	%%zmm11,%%zmm8					\n\t	vmovaps	%%zmm15,%%zmm12					\n\t"/* cpy hi into lo-reg */\
		"vfmsub231pd 0x100(%%rdx),%%zmm2,%%zmm0	\n\t	vfmsub231pd 0x140(%%rdx),%%zmm6,%%zmm4	\n\t	vfmsub231pd 0x180(%%rdx),%%zmm10,%%zmm8	\n\t	vfmsub231pd 0x1c0(%%rdx),%%zmm14,%%zmm12\n\t"/* lo = fma(a,b, -hi) */\
		/* NOTE: Cannot simply sum raw (unnormalized) hi-outouts here with one from FMA-mul below, because hi-terms generally
		have 53 significant bits but == 0 w.r.to different powers of 2, thus we lose critical bits on the low end in sch an add. */\
		"vaddpd			%%zmm29,%%zmm3,%%zmm2	\n\t	vaddpd			%%zmm29,%%zmm7,%%zmm6	\n\t	vaddpd			%%zmm29,%%zmm11,%%zmm10	\n\t	vaddpd			%%zmm29,%%zmm15,%%zmm14	\n\t"\
		"vsubpd			%%zmm29,%%zmm2,%%zmm2	\n\t	vsubpd			%%zmm29,%%zmm6,%%zmm6	\n\t	vsubpd			%%zmm29,%%zmm10,%%zmm10	\n\t	vsubpd			%%zmm29,%%zmm14,%%zmm14	\n\t"/* hh = hi +- crnd50 to round-to-nearest-multiple-of-base */\
		"vsubpd			%%zmm2,%%zmm3,%%zmm3	\n\t	vsubpd			%%zmm6,%%zmm7,%%zmm7	\n\t	vsubpd			%%zmm10,%%zmm11,%%zmm11	\n\t	vsubpd			%%zmm14,%%zmm15,%%zmm15	\n\t"/* hi - hh gives backward-carry... */\
		"vaddpd			%%zmm3,%%zmm0,%%zmm0	\n\t	vaddpd			%%zmm7,%%zmm4,%%zmm4	\n\t	vaddpd			%%zmm11,%%zmm8,%%zmm8	\n\t	vaddpd			%%zmm15,%%zmm12,%%zmm12	\n\t"/* ...which we add to lo. */\
		"vfmadd231pd	%%zmm30,%%zmm1,%%zmm0	\n\t	vfmadd231pd	%%zmm30,%%zmm5,%%zmm4		\n\t	vfmadd231pd	%%zmm30,%%zmm9,%%zmm8		\n\t	vfmadd231pd	%%zmm30,%%zmm13,%%zmm12		\n\t"/* Add lo to hi-output of previous lo:hi pair, which also needs *= binv */\
													/*** prod[1] (partial sum) in zmm0 ***/\
		/* qinv.lo*lo.hi: a = qinv.lo, b = lo.hi: */\
		"vmovaps	0x100(%%rax),%%zmm1			\n\t	vmovaps	0x140(%%rax),%%zmm5				\n\t	vmovaps	0x180(%%rax),%%zmm9				\n\t	vmovaps	0x1c0(%%rax),%%zmm13			\n\t"/* load lo.hi into zmm1... */\
		"     vmulpd	(%%rdx),%%zmm1,%%zmm3	\n\t	     vmulpd	0x040(%%rdx),%%zmm5,%%zmm7	\n\t	     vmulpd	0x080(%%rdx),%%zmm9,%%zmm11	\n\t	     vmulpd	0x0c0(%%rdx),%%zmm13,%%zmm15\n\t"/* hi = fma(a,b, 0  ) - use VMULPD instead of VFMADD231PD, since addend = 0. */\
		"vmovaps	%%zmm3,%%zmm2				\n\t	vmovaps	%%zmm7,%%zmm6					\n\t	vmovaps	%%zmm11,%%zmm10					\n\t	vmovaps	%%zmm15,%%zmm14					\n\t"/* cpy hi into lo-reg */\
		"vfmsub231pd	(%%rdx),%%zmm1,%%zmm2	\n\t	vfmsub231pd	0x040(%%rdx),%%zmm5,%%zmm6	\n\t	vfmsub231pd	0x080(%%rdx),%%zmm9,%%zmm10	\n\t	vfmsub231pd	0x0c0(%%rdx),%%zmm13,%%zmm14\n\t"/* lo = fma(a,b, -hi) */\
		"vaddpd			%%zmm29,%%zmm3,%%zmm1	\n\t	vaddpd			%%zmm29,%%zmm7,%%zmm5	\n\t	vaddpd			%%zmm29,%%zmm11,%%zmm9	\n\t	vaddpd			%%zmm29,%%zmm15,%%zmm13	\n\t"\
		"vsubpd			%%zmm29,%%zmm1,%%zmm1	\n\t	vsubpd			%%zmm29,%%zmm5,%%zmm5	\n\t	vsubpd			%%zmm29,%%zmm9,%%zmm9	\n\t	vsubpd			%%zmm29,%%zmm13,%%zmm13	\n\t"/* hh = hi +- crnd50 to round-to-nearest-multiple-of-base */\
		"vsubpd			%%zmm1,%%zmm3,%%zmm3	\n\t	vsubpd			%%zmm5,%%zmm7,%%zmm7	\n\t	vsubpd			%%zmm9,%%zmm11,%%zmm11	\n\t	vsubpd			%%zmm13,%%zmm15,%%zmm15	\n\t"/* hi - hh gives backward-carry... */\
		"vaddpd			%%zmm3,%%zmm2,%%zmm2	\n\t	vaddpd			%%zmm7,%%zmm6,%%zmm6	\n\t	vaddpd			%%zmm11,%%zmm10,%%zmm10	\n\t	vaddpd			%%zmm15,%%zmm14,%%zmm14	\n\t"/* ...which we add to lo... */\
		"vaddpd			%%zmm0,%%zmm2,%%zmm0	\n\t	vaddpd			%%zmm4,%%zmm6,%%zmm4	\n\t	vaddpd			%%zmm8,%%zmm10,%%zmm8	\n\t	vaddpd			%%zmm12,%%zmm14,%%zmm12	\n\t"/* ...and add that to zmm0 (prod[1]). */\
		/* Since x = prod[1] was obtained by summing two base-normalized FMA-mul low-half outputs, need to still normalize it: */\
		"     vmulpd	%%zmm30,%%zmm0,%%zmm1	\n\t	     vmulpd	%%zmm30,%%zmm4,%%zmm5		\n\t	     vmulpd	%%zmm30,%%zmm8,%%zmm9		\n\t	     vmulpd	%%zmm30,%%zmm12,%%zmm13		\n\t"/* x*binv */\
		"vrndscalepd $1,%%zmm1,%%zmm1			\n\t	vrndscalepd $1,%%zmm5,%%zmm5			\n\t	vrndscalepd $1,%%zmm9,%%zmm9			\n\t	vrndscalepd $1,%%zmm13,%%zmm13			\n\t"/* cy = floor(x*binv) */\
		"vfnmadd231pd	%%zmm31,%%zmm1,%%zmm0	\n\t	vfnmadd231pd	%%zmm31,%%zmm5,%%zmm4	\n\t	vfnmadd231pd	%%zmm31,%%zmm9,%%zmm8	\n\t	vfnmadd231pd	%%zmm31,%%zmm13,%%zmm12	\n\t"/* x = FMA(cy,-base, x); x in [0,base) */\
		"vmovaps	%%zmm0,0x100(%%rax)			\n\t	vmovaps	%%zmm4,0x140(%%rax)				\n\t	vmovaps	%%zmm8,0x180(%%rax)				\n\t	vmovaps	%%zmm12,0x1c0(%%rax)			\n\t"/* write prod[1] (lo.hi) */\
													/*** lo.lo in (%%rax), lo.hi in zmm0. ***/\
	/********************* lo = UMULH(q,lo) : *****************************************************************************
		Bitranges (B = base):													Subproducts:
											|----------- 2B-1: 0 -----------|	q.lo * lo.lo
							|----------- 3B-1: B -----------|					q.hi * lo.lo, q.lo * lo.hi
			|----------- 4B-1:2B -----------|									q.hi * lo.hi
			|----------- OUTPUT ------------|
		Thus may only need high half of the (q.lo * lo.lo) subproducts - alas, I see few other work-saving possibilities.
		But wait - as long as our base is a few bits less than the 53-bits of an IEEE double mantissa, the hi-part MUL of
		a double-wide FMA-mul with its 53 significant bits will capture a disproportionate share of the full-prodict bits.
		Thus it my be possible to simplify the computation of the lower-half carryin-bits by computing the hi bits of the
		(q.lo * lo.lo) subproduct, feeding that result * binv as the addend to the (q.hi * lo.lo) subproduct (again sans
		a corresponding low-part FMA-mul), and adding that result in turn to the high-part-only (q.lo * lo.hi) subproduct.
		We can then round off the excess bits at the low end of the result, and proceed to compute the exact double-wide
		(q.hi * lo.hi) subproduct.
	*********************************************************************************************************************/\
		"vmulpd		(%%rax),%%zmm30,%%zmm1		\n\t	vmulpd	0x040(%%rax),%%zmm30,%%zmm5		\n\t	vmulpd	0x080(%%rax),%%zmm30,%%zmm9		\n\t	vmulpd	0x0c0(%%rax),%%zmm30,%%zmm13	\n\t"/* Premultiply lo.lo by binv */\
		"vmovaps	(%%rax),%%zmm2				\n\t	vmovaps	0x040(%%rax),%%zmm6				\n\t	vmovaps	0x080(%%rax),%%zmm10			\n\t	vmovaps	0x0c0(%%rax),%%zmm14			\n\t"/* load lo.lo and save copy in zmm2 */\
		"vmovaps	0x100(%%rax),%%zmm3			\n\t	vmovaps	0x140(%%rax),%%zmm7				\n\t	vmovaps	0x180(%%rax),%%zmm11			\n\t	vmovaps	0x1c0(%%rax),%%zmm15			\n\t"/* preload lo.hi */\
		/* q.lo*lo.lo: a = q.lo; b = lo.lo: */\
		"vmulpd		(%%rbx),%%zmm1,%%zmm1		\n\t	vmulpd		0x040(%%rbx),%%zmm5,%%zmm5	\n\t	vmulpd		0x080(%%rbx),%%zmm9,%%zmm9	\n\t	vmulpd		0x0c0(%%rbx),%%zmm13,%%zmm13\n\t"/* [high 53 bits of q.lo * lo.lo]*binv... */\
		"vfmadd231pd 0x100(%%rbx),%%zmm2,%%zmm1	\n\t	vfmadd231pd 0x140(%%rbx),%%zmm6,%%zmm5	\n\t	vfmadd231pd 0x180(%%rbx),%%zmm10,%%zmm9	\n\t	vfmadd231pd 0x1c0(%%rbx),%%zmm14,%%zmm13\n\t"/* ...feed that as addend to (q.hi * lo.lo) subproduct (again high-53-bits only)... */\
		"vfmadd231pd      (%%rbx),%%zmm3,%%zmm1	\n\t	vfmadd231pd 0x040(%%rbx),%%zmm7,%%zmm5	\n\t	vfmadd231pd 0x080(%%rbx),%%zmm11,%%zmm9	\n\t	vfmadd231pd 0x0c0(%%rbx),%%zmm15,%%zmm13\n\t"/* ...add that result in turn to the high-part-only (q.lo * lo.hi) subproduct. */\
		"vmulpd		%%zmm30,%%zmm1,%%zmm0		\n\t	vmulpd		%%zmm30,%%zmm5,%%zmm4		\n\t	vmulpd		%%zmm30,%%zmm9,%%zmm8		\n\t	vmulpd		%%zmm30,%%zmm13,%%zmm12		\n\t"/* tmp *= binv */\
	"vrndscalepd $3,%%zmm0,%%zmm0				\n\t	vrndscalepd $3,%%zmm4,%%zmm4			\n\t	vrndscalepd $3,%%zmm8,%%zmm8			\n\t	vrndscalepd $3,%%zmm12,%%zmm12			\n\t"/* round_toward_0(tmp*binv), this gets added to low half of (q.hi * lo.hi) */\
		/* q.hi*lo.hi: a = q.hi, b = lo.hi: */\
		"     vmulpd 0x100(%%rbx),%%zmm3,%%zmm1	\n\t	     vmulpd 0x140(%%rbx),%%zmm7,%%zmm5	\n\t	     vmulpd 0x180(%%rbx),%%zmm11,%%zmm9	\n\t	     vmulpd 0x1c0(%%rbx),%%zmm15,%%zmm13\n\t"/* hi = fma(a,b, 0  ) - use VMULPD instead of VFMADD231PD, since addend = 0. */\
		"vcmppd	$14,0x400(%%rax),%%zmm1,%%k1	\n\t	vcmppd	$14,0x440(%%rax),%%zmm5,%%k2	\n\t	vcmppd	$14,0x480(%%rax),%%zmm9,%%k3	\n\t	vcmppd	$14,0x4c0(%%rax),%%zmm13,%%k4	\n\t"/* bitmask = (lo53 > hi53)? In AVX512 version, mask-regs k1-4 replace dest-regs m2,6,10,14 */\
		"vmovaps	%%zmm1,%%zmm2				\n\t	vmovaps	%%zmm5,%%zmm6					\n\t	vmovaps	%%zmm9,%%zmm10					\n\t	vmovaps	%%zmm13,%%zmm14					\n\t"/* cpy hi into lo-reg */\
		"vfmsub231pd 0x100(%%rbx),%%zmm3,%%zmm2	\n\t	vfmsub231pd 0x140(%%rbx),%%zmm7,%%zmm6	\n\t	vfmsub231pd 0x180(%%rbx),%%zmm11,%%zmm10\n\t	vfmsub231pd 0x1c0(%%rbx),%%zmm15,%%zmm14\n\t"/* lo = fma(a,b, -hi) */\
		"vaddpd			%%zmm29,%%zmm1,%%zmm3	\n\t	vaddpd			%%zmm29,%%zmm5,%%zmm7	\n\t	vaddpd			%%zmm29,%%zmm9,%%zmm11	\n\t	vaddpd			%%zmm29,%%zmm13,%%zmm15	\n\t"\
		"vsubpd			%%zmm29,%%zmm3,%%zmm3	\n\t	vsubpd			%%zmm29,%%zmm7,%%zmm7	\n\t	vsubpd			%%zmm29,%%zmm11,%%zmm11	\n\t	vsubpd			%%zmm29,%%zmm15,%%zmm15	\n\t"/* hh = hi +- crnd50 to round-to-nearest-multiple-of-base */\
		"vsubpd			%%zmm3,%%zmm1,%%zmm1	\n\t	vsubpd			%%zmm7,%%zmm5,%%zmm5	\n\t	vsubpd			%%zmm11,%%zmm9,%%zmm9	\n\t	vsubpd			%%zmm15,%%zmm13,%%zmm13	\n\t"/* hi - hh gives backward-carry... */\
		"vmulpd		%%zmm30,%%zmm3,%%zmm3		\n\t	vmulpd		%%zmm30,%%zmm7,%%zmm7		\n\t	vmulpd		%%zmm30,%%zmm11,%%zmm11		\n\t	vmulpd		%%zmm30,%%zmm15,%%zmm15		\n\t"/* hh *= binv */\
		"vaddpd			%%zmm1,%%zmm2,%%zmm2	\n\t	vaddpd			%%zmm5,%%zmm6,%%zmm6	\n\t	vaddpd			%%zmm9,%%zmm10,%%zmm10	\n\t	vaddpd			%%zmm13,%%zmm14,%%zmm14	\n\t"/* ...which we add to lo... */\
		"vaddpd			%%zmm0,%%zmm2,%%zmm0	\n\t	vaddpd			%%zmm4,%%zmm6,%%zmm4	\n\t	vaddpd			%%zmm8,%%zmm10,%%zmm8	\n\t	vaddpd			%%zmm12,%%zmm14,%%zmm12	\n\t"/* ...and add that to zmm0. */\
													/*** lo.lo in zmm0, lo.hi in zmm3 ***/\
		/* If h < l (already done above using leading significant 53 bits of both operands, i.e. the unnormalized
		high outputs of the respective 2-word FMA-muls), calculate h-l+q; otherwise h-l. Result is in [0, q). */\
		"vmovaps	0x200(%%rax),%%zmm1			\n\t	vmovaps	0x240(%%rax),%%zmm5				\n\t	vmovaps	0x280(%%rax),%%zmm9				\n\t	vmovaps	0x2c0(%%rax),%%zmm13			\n\t"/* hi.lo */\
		"vmovaps	0x300(%%rax),%%zmm2			\n\t	vmovaps	0x340(%%rax),%%zmm6				\n\t	vmovaps	0x380(%%rax),%%zmm10			\n\t	vmovaps	0x3c0(%%rax),%%zmm14			\n\t"/* hi.hi */\
		"vsubpd	%%zmm0,%%zmm1,%%zmm1			\n\t	vsubpd	%%zmm4,%%zmm5,%%zmm5			\n\t	vsubpd	%%zmm8,%%zmm9,%%zmm9			\n\t	vsubpd	%%zmm12,%%zmm13,%%zmm13			\n\t"/* (hi - lo).lo, zmm0 FREE */\
		"vsubpd	%%zmm3,%%zmm2,%%zmm2			\n\t	vsubpd	%%zmm7,%%zmm6,%%zmm6			\n\t	vsubpd	%%zmm11,%%zmm10,%%zmm10			\n\t	vsubpd	%%zmm15,%%zmm14,%%zmm14			\n\t"/* (hi - lo).hi, zmm3 FREE */\
		"vmovaps	0x100(%%rbx),%%zmm0			\n\t	vmovaps	0x140(%%rbx),%%zmm4				\n\t	vmovaps	0x180(%%rbx),%%zmm8				\n\t	vmovaps	0x1c0(%%rbx),%%zmm12			\n\t"/* q.hi */\
	   "vaddpd (%%rbx),%%zmm1,%%zmm1%{%%k1%}	\n\t vaddpd	0x040(%%rbx),%%zmm5,%%zmm5%{%%k2%}	\n\t vaddpd	0x080(%%rbx),%%zmm9,%%zmm9%{%%k3%}	\n\t vaddpd	0x0c0(%%rbx),%%zmm13,%%zmm13%{%%k4%}\n\t"/* xlo = (hi-lo).lo + (q.lo & bitmask) */\
	   "vaddpd	%%zmm0,%%zmm2,%%zmm2%{%%k1%}	\n\t	vaddpd	%%zmm4,%%zmm6,%%zmm6%{%%k2%}	\n\t	vaddpd	%%zmm8,%%zmm10,%%zmm10%{%%k3%}	\n\t	vaddpd	%%zmm12,%%zmm14,%%zmm14%{%%k4%}	\n\t"/* xhi = (hi-lo).hi + (q.hi & bitmask) */\
													/*** x.lo,hi in zmm1,zmm2; q.hi in zmm0 ***/\
	/* if((pshift >> j) & (uint64)1) { */\
		"movq	%[__pshift],%%rsi				\n\t"\
		"movslq	%[__j],%%rcx					\n\t"\
		"shrq	%%cl,%%rsi						\n\t"\
		"andq	$0x1,%%rsi						\n\t"\
	"je twopmodq100_2wdq32_gcc64				\n\t"\
		"vaddpd	%%zmm1,%%zmm1,%%zmm1			\n\t	vaddpd	%%zmm5,%%zmm5,%%zmm5			\n\t	vaddpd	%%zmm9,%%zmm9,%%zmm9			\n\t	vaddpd	%%zmm13,%%zmm13,%%zmm13			\n\t"/* 2*x.lo */\
		"vaddpd	%%zmm2,%%zmm2,%%zmm2			\n\t	vaddpd	%%zmm6,%%zmm6,%%zmm6			\n\t	vaddpd	%%zmm10,%%zmm10,%%zmm10			\n\t	vaddpd	%%zmm14,%%zmm14,%%zmm14			\n\t"/* 2*x.hi */\
		"vmovaps	%%zmm1,%%zmm3				\n\t	vmovaps	%%zmm5,%%zmm7					\n\t	vmovaps	%%zmm9,%%zmm11					\n\t	vmovaps	%%zmm13,%%zmm15					\n\t"/* cpy of 2*x.lo */\
		"vfmadd231pd %%zmm31,%%zmm2,%%zmm3		\n\t	vfmadd231pd %%zmm31,%%zmm6,%%zmm7		\n\t	vfmadd231pd %%zmm31,%%zmm10,%%zmm11		\n\t	vfmadd231pd %%zmm31,%%zmm14,%%zmm15		\n\t"/* 2*(x.lo + x.hi*base) to get most-sig 53 bits of doubled-x for the following comparison */\
		/* If x > q, subtract q: */\
		"vcmppd $13,0x200(%%rbx),%%zmm3,%%k1	\n\t	vcmppd $13,0x240(%%rbx),%%zmm7,%%k2		\n\t	vcmppd $13,0x280(%%rbx),%%zmm11,%%k3	\n\t	vcmppd $13,0x2c0(%%rbx),%%zmm15,%%k4	\n\t"/* bitmask = (xhi53 >= qhi53)? In AVX512 version, mask-regs k1-4 replace dest-regs m3,7,11,15 */\
	   "vsubpd (%%rbx),%%zmm1,%%zmm1%{%%k1%}	\n\t vsubpd	0x040(%%rbx),%%zmm5,%%zmm5%{%%k2%}	\n\t vsubpd	0x080(%%rbx),%%zmm9,%%zmm9%{%%k3%}	\n\t vsubpd	0x0c0(%%rbx),%%zmm13,%%zmm13%{%%k4%}\n\t"/* xlo = 2*x.lo - (q.lo & bitmask) */\
	   "vsubpd	%%zmm0,%%zmm2,%%zmm2%{%%k1%}	\n\t	vsubpd	%%zmm4,%%zmm6,%%zmm6%{%%k2%}	\n\t	vsubpd	%%zmm8,%%zmm10,%%zmm10%{%%k3%}	\n\t	vsubpd	%%zmm12,%%zmm14,%%zmm14%{%%k4%}	\n\t"/* xhi = 2*x.hi - (q.hi & bitmask) */\
	"twopmodq100_2wdq32_gcc64:					\n\t"\
		"vmovaps	%%zmm1,     (%%rax)			\n\t	vmovaps	%%zmm5,0x040(%%rax)				\n\t	vmovaps	%%zmm9,0x080(%%rax)				\n\t	vmovaps	%%zmm13,0x0c0(%%rax)			\n\t"/* write x.lo */\
		"vmovaps	%%zmm2,0x100(%%rax)			\n\t	vmovaps	%%zmm6,0x140(%%rax)				\n\t	vmovaps	%%zmm10,0x180(%%rax)			\n\t	vmovaps	%%zmm14,0x1c0(%%rax)			\n\t"/* write x.hi */\
	/* } */\
		:					/* outputs: none */\
		: [__aq0]	 "m" (Xaq0)	/* All inputs from memory addresses here */\
		 ,[__aqinv0] "m" (Xaqinv0)	\
		 ,[__ax0]	 "m" (Xax0)		\
		 ,[__base] "m" (Xbase)	\
		 ,[__binv] "m" (Xbinv)	\
		 ,[__crnd] "m" (Xcrnd)\
		 ,[__pshift] "m" (Xpshift)	\
		 ,[__j]		 "m" (Xj)		\
		: "cc","memory","cl","rax","rbx","rcx","rdx","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15", "xmm29","xmm30","xmm31"	/* Clobbered registers */\
	);\
	}

  #else	// ALL_FMA = true; this version uses FMA-with-unity-input in place of ADDPD,SUBPD:
#error *** To-do! propagate more FMAs into ALL_FMA version of macro! ***
	#define SSE2_twopmodq100_modmul_q32(Xaq0,Xaqinv0,Xax0,Xbase,Xbinv,Xcrnd,Xpshift,Xj)\
	{\
	__asm__ volatile (\
		"movq	%[__base],%%rax		\n\t	vbroadcastsd	(%%rax),%%zmm31	\n\t"/* BASE */\
		"movq	%[__binv],%%rbx		\n\t	vbroadcastsd	(%%rbx),%%zmm30	\n\t"/* BINV */\
		"movq	%[__crnd],%%rcx		\n\t	vbroadcastsd	(%%rcx),%%zmm29	\n\t"/* CRND */\
		"vmulpd	%%zmm31,%%zmm29,%%zmm29			\n\t"/* CRND*BASE */\
		"vmulpd	%%zmm31,%%zmm30,%%zmm28			\n\t"/* 1.0 ... no good way to generate 1.0 on-the-fly, so use base * binv */\
		"movq	%[__ax0],%%rax					\n\t"/* &x */\
		"movq	%[__aq0],%%rbx					\n\t"\
		"movq	%[__aqinv0],%%rdx				\n\t"\
		/*** STREAM 0 ***/							/*** STREAM 1 ***/							/*** STREAM 2 ***/								/*** STREAM 3 ***/					\
		/* Inputs [a|c|e|g]lo0-2 enter in zmm[0-1,4-5,8-9,12-13], resp.: */\
	/************ hi:lo = SQR_LOHI(x) : **************/\
		"vmovaps	(%%rax),%%zmm1				\n\t	vmovaps	0x040(%%rax),%%zmm5				\n\t	vmovaps	0x080(%%rax),%%zmm9				\n\t	vmovaps	0x0c0(%%rax),%%zmm13			\n\t"/* load x.lo */\
		"vaddpd			%%zmm1,%%zmm1,%%zmm2	\n\t	vaddpd			%%zmm5,%%zmm5,%%zmm6	\n\t	vaddpd			%%zmm9,%%zmm9,%%zmm10	\n\t	vaddpd			%%zmm13,%%zmm13,%%zmm14	\n\t"/* 2*lo */\
		/* lo*lo: a = b = x.lo: */\
		"     vmulpd	%%zmm1,%%zmm1,%%zmm3	\n\t	     vmulpd	%%zmm5,%%zmm5,%%zmm7		\n\t	     vmulpd	%%zmm9,%%zmm9,%%zmm11		\n\t	     vmulpd	%%zmm13,%%zmm13,%%zmm15		\n\t"/* hi = fma(a,b, 0  ) - use VMULPD instead of VFMADD231PD, since addend = 0. */\
		"vmovaps	%%zmm3,%%zmm0				\n\t	vmovaps	%%zmm7,%%zmm4					\n\t	vmovaps	%%zmm11,%%zmm8					\n\t	vmovaps	%%zmm15,%%zmm12					\n\t"/* cpy hi into lo-reg */\
		"vfmsub231pd	%%zmm1,%%zmm1,%%zmm0	\n\t	vfmsub231pd	%%zmm5,%%zmm5,%%zmm4		\n\t	vfmsub231pd	%%zmm9,%%zmm9,%%zmm8		\n\t	vfmsub231pd	%%zmm13,%%zmm13,%%zmm12		\n\t"/* lo = fma(a,b, -hi) */\
		"vaddpd			%%zmm29,%%zmm3,%%zmm1	\n\t	vaddpd			%%zmm29,%%zmm7,%%zmm5	\n\t	vaddpd			%%zmm29,%%zmm11,%%zmm9	\n\t	vaddpd			%%zmm29,%%zmm15,%%zmm13	\n\t"\
		"vsubpd			%%zmm29,%%zmm1,%%zmm1	\n\t	vsubpd			%%zmm29,%%zmm5,%%zmm5	\n\t	vsubpd			%%zmm29,%%zmm9,%%zmm9	\n\t	vsubpd			%%zmm29,%%zmm13,%%zmm13	\n\t"/* hh = hi +- crnd50 to round-to-nearest-multiple-of-base */\
		"vsubpd			%%zmm1,%%zmm3,%%zmm3	\n\t	vsubpd			%%zmm5,%%zmm7,%%zmm7	\n\t	vsubpd			%%zmm9,%%zmm11,%%zmm11	\n\t	vsubpd			%%zmm13,%%zmm15,%%zmm15	\n\t"/* hi - hh gives backward-carry... */\
		"vfmadd231pd	%%zmm28,%%zmm3,%%zmm0	\n\t	vfmadd231pd	%%zmm28,%%zmm7,%%zmm4		\n\t	vfmadd231pd	%%zmm28,%%zmm11,%%zmm8		\n\t	vfmadd231pd	%%zmm28,%%zmm15,%%zmm12		\n\t"/* ...which we add to lo (prod[0]). */\
		"vmovaps	%%zmm0,     (%%rax)			\n\t	vmovaps	%%zmm4,0x040(%%rax)				\n\t	vmovaps	%%zmm8,0x080(%%rax)				\n\t	vmovaps	%%zmm12,0x0c0(%%rax)			\n\t"/* write prod[0]; this overwrites x.lo, but only need 2*x.lo for rest of squaring, that is in zmm2 */\
													/*** prod[0] in zmm0, hi*base in zmm1. ***/\
		/* 2*lo*hi: a = 2*x.lo, b = x.hi: */\
		"     vmulpd 0x100(%%rax),%%zmm2,%%zmm3	\n\t	     vmulpd 0x140(%%rax),%%zmm6,%%zmm7	\n\t	     vmulpd 0x180(%%rax),%%zmm10,%%zmm11\n\t	     vmulpd 0x1c0(%%rax),%%zmm14,%%zmm15\n\t"/* hi = fma(a,b, 0  ) - use VMULPD instead of VFMADD231PD, since addend = 0. */\
		"vmovaps	%%zmm3,%%zmm0				\n\t	vmovaps	%%zmm7,%%zmm4					\n\t	vmovaps	%%zmm11,%%zmm8					\n\t	vmovaps	%%zmm15,%%zmm12					\n\t"/* cpy hi into lo-reg */\
		"vfmsub231pd 0x100(%%rax),%%zmm2,%%zmm0	\n\t	vfmsub231pd 0x140(%%rax),%%zmm6,%%zmm4	\n\t	vfmsub231pd 0x180(%%rax),%%zmm10,%%zmm8	\n\t	vfmsub231pd 0x1c0(%%rax),%%zmm14,%%zmm12\n\t"/* lo = fma(a,b, -hi) */\
		"     vmulpd	%%zmm30,%%zmm3,%%zmm2	\n\t	     vmulpd	%%zmm30,%%zmm7,%%zmm6		\n\t	     vmulpd	%%zmm30,%%zmm11,%%zmm10		\n\t	     vmulpd	%%zmm30,%%zmm15,%%zmm14		\n\t"/* hi*binv */\
		"vrndscalepd $1,%%zmm2,%%zmm2			\n\t	vrndscalepd $1,%%zmm6,%%zmm6			\n\t	vrndscalepd $1,%%zmm10,%%zmm10			\n\t	vrndscalepd $1,%%zmm14,%%zmm14			\n\t"/* hh = floor(hi*binv); This part remains in hi */\
		"vfmadd231pd	%%zmm30,%%zmm1,%%zmm0	\n\t	vfmadd231pd	%%zmm30,%%zmm5,%%zmm4		\n\t	vfmadd231pd	%%zmm30,%%zmm9,%%zmm8		\n\t	vfmadd231pd	%%zmm30,%%zmm13,%%zmm12		\n\t"/* Add lo to hi-output of previous lo:hi pair, which also needs *= binv */\
		"vfnmadd231pd	%%zmm31,%%zmm2,%%zmm3	\n\t	vfnmadd231pd	%%zmm31,%%zmm6,%%zmm7	\n\t	vfnmadd231pd	%%zmm31,%%zmm10,%%zmm11	\n\t	vfnmadd231pd	%%zmm31,%%zmm14,%%zmm15	\n\t"/* cy = FMA(hh,-base, hi); backward-carry from hi into lo... */\
		"vaddpd			%%zmm3,%%zmm0,%%zmm1	\n\t	vaddpd			%%zmm7,%%zmm4,%%zmm5	\n\t	vaddpd			%%zmm11,%%zmm8,%%zmm9	\n\t	vaddpd			%%zmm15,%%zmm12,%%zmm13	\n\t"/* ...which we add to lo (prod[1]). */\
	/* MONT_MUL *requires* the h:l = x^2 results to be properly normalized in order to yield correct result;
	at this point prod[1] can be > base, thus must do a normalize-and-carry: */\
		"     vmulpd	%%zmm30,%%zmm1,%%zmm3	\n\t	     vmulpd	%%zmm30,%%zmm5,%%zmm7		\n\t	     vmulpd	%%zmm30,%%zmm9,%%zmm11		\n\t	     vmulpd	%%zmm30,%%zmm13,%%zmm15		\n\t"/* prod[1]*binv */\
		"vrndscalepd $1,%%zmm3,%%zmm3			\n\t	vrndscalepd $1,%%zmm7,%%zmm7			\n\t	vrndscalepd $1,%%zmm11,%%zmm11			\n\t	vrndscalepd $1,%%zmm15,%%zmm15			\n\t"/* cy = floor(prod[1]*binv) */\
		"vfnmadd231pd	%%zmm31,%%zmm3,%%zmm1	\n\t	vfnmadd231pd	%%zmm31,%%zmm7,%%zmm5	\n\t	vfnmadd231pd	%%zmm31,%%zmm11,%%zmm9	\n\t	vfnmadd231pd	%%zmm31,%%zmm15,%%zmm13	\n\t"/* prod[1] -= cy*base */\
		"vaddpd			%%zmm3,%%zmm2,%%zmm2	\n\t	vaddpd			%%zmm7,%%zmm6,%%zmm6	\n\t	vaddpd			%%zmm11,%%zmm10,%%zmm10	\n\t	vaddpd			%%zmm15,%%zmm14,%%zmm14	\n\t"/* prod[2] += cy      */\
		"     vmulpd	%%zmm31,%%zmm2,%%zmm2	\n\t	     vmulpd	%%zmm31,%%zmm6,%%zmm6		\n\t	     vmulpd	%%zmm31,%%zmm10,%%zmm10		\n\t	     vmulpd	%%zmm31,%%zmm14,%%zmm14		\n\t"/* rescale hh *= base */\
		"vmovaps	0x100(%%rax),%%zmm3			\n\t	vmovaps	0x140(%%rax),%%zmm7				\n\t	vmovaps	0x180(%%rax),%%zmm11			\n\t	vmovaps	0x1c0(%%rax),%%zmm15			\n\t"/* reload hi into zmm3... */\
		"vmovaps	%%zmm1,0x100(%%rax)			\n\t	vmovaps	%%zmm5,0x140(%%rax)				\n\t	vmovaps	%%zmm9,0x180(%%rax)				\n\t	vmovaps	%%zmm13,0x1c0(%%rax)			\n\t"/* ...before overwriting with prod[1] */\
													/*** prod[1] in zmm1, hi*base in zmm2, x_in.hi in zmm3. ***/\
		/* hi*hi: a = b = x.hi: */\
		"     vmulpd	%%zmm3,%%zmm3,%%zmm0	\n\t	     vmulpd	%%zmm7,%%zmm7,%%zmm4		\n\t	     vmulpd	%%zmm11,%%zmm11,%%zmm8		\n\t	     vmulpd	%%zmm15,%%zmm15,%%zmm12		\n\t"/* hi = fma(a,b, 0  ) - use VMULPD instead of VFMADD231PD, since addend = 0. */\
		"vmovaps	%%zmm0,%%zmm1				\n\t	vmovaps	%%zmm4,%%zmm5					\n\t	vmovaps	%%zmm8,%%zmm9					\n\t	vmovaps	%%zmm12,%%zmm13					\n\t"/* cpy hi into lo-reg */\
		"vfmsub231pd	%%zmm3,%%zmm3,%%zmm1	\n\t	vfmsub231pd	%%zmm7,%%zmm7,%%zmm5		\n\t	vfmsub231pd	%%zmm11,%%zmm11,%%zmm9		\n\t	vfmsub231pd	%%zmm15,%%zmm15,%%zmm13		\n\t"/* lo = fma(a,b, -hi) */\
		"vaddpd			%%zmm29,%%zmm0,%%zmm3	\n\t	vaddpd			%%zmm29,%%zmm4,%%zmm7	\n\t	vaddpd			%%zmm29,%%zmm8,%%zmm11	\n\t	vaddpd			%%zmm29,%%zmm12,%%zmm15	\n\t"\
		"vsubpd			%%zmm29,%%zmm3,%%zmm3	\n\t	vsubpd			%%zmm29,%%zmm7,%%zmm7	\n\t	vsubpd			%%zmm29,%%zmm11,%%zmm11	\n\t	vsubpd			%%zmm29,%%zmm15,%%zmm15	\n\t"/* hh = hi +- crnd50 to round-to-nearest-multiple-of-base */\
		"vfmadd231pd	%%zmm30,%%zmm2,%%zmm1	\n\t	vfmadd231pd	%%zmm30,%%zmm6,%%zmm5		\n\t	vfmadd231pd	%%zmm30,%%zmm10,%%zmm9		\n\t	vfmadd231pd	%%zmm30,%%zmm14,%%zmm13		\n\t"/* Add lo to hi-output of previous lo:hi pair, which also needs *= binv */\
		"vmovaps	%%zmm0,0x400(%%rax)			\n\t	vmovaps	%%zmm4,0x440(%%rax)				\n\t	vmovaps	%%zmm8,0x480(%%rax)				\n\t	vmovaps	%%zmm12,0x4c0(%%rax)			\n\t"/* write hi-approximant to full 2-word product into scratch slot set aside for it */\
		"vsubpd			%%zmm3,%%zmm0,%%zmm0	\n\t	vsubpd			%%zmm7,%%zmm4,%%zmm4	\n\t	vsubpd			%%zmm11,%%zmm8,%%zmm8	\n\t	vsubpd			%%zmm15,%%zmm12,%%zmm12	\n\t"/* hi - hh gives backward-carry... */\
		"vaddpd			%%zmm0,%%zmm1,%%zmm2	\n\t	vaddpd			%%zmm4,%%zmm5,%%zmm6	\n\t	vaddpd			%%zmm8,%%zmm9,%%zmm10	\n\t	vaddpd			%%zmm12,%%zmm13,%%zmm14	\n\t"/* ...which we add to lo (prod[2]). */\
		"     vmulpd	%%zmm3,%%zmm30,%%zmm3	\n\t	     vmulpd	%%zmm7,%%zmm30,%%zmm7		\n\t	     vmulpd	%%zmm11,%%zmm30,%%zmm11		\n\t	     vmulpd	%%zmm15,%%zmm30,%%zmm15		\n\t"/* prod[3] = hh*binv */\
		"vmovaps	%%zmm2,0x200(%%rax)			\n\t	vmovaps	%%zmm6,0x240(%%rax)				\n\t	vmovaps	%%zmm10,0x280(%%rax)			\n\t	vmovaps	%%zmm14,0x2c0(%%rax)			\n\t"/* write prod[2] */\
		"vmovaps	%%zmm3,0x300(%%rax)			\n\t	vmovaps	%%zmm7,0x340(%%rax)				\n\t	vmovaps	%%zmm11,0x380(%%rax)			\n\t	vmovaps	%%zmm15,0x3c0(%%rax)			\n\t"/* write prod[3] */\
													/*** prod[2,3] in zmm2,3. ***/\
	/********************* lo = MULL(qinv,lo) : ***************************************************************************
		Bitranges (B = base):													Subproducts:
											|----------- 2B-1: 0 -----------|	qinv.lo * lo.lo
							|----------- 3B-1: B -----------|					qinv.hi * lo.lo, qinv.lo * lo.hi
			|----------- 4B-1:2B -----------|									qinv.hi * lo.hi
											|----------- OUTPUT ------------|
		Thus only need full 2B-wide (qinv.lo * lo.lo) subproduct and low halves of the (qinv.hi * lo.lo), (qinv.lo * lo.hi)
		subproducts - would that AVX-512F supported VPMULLQ! [Alas that needs AVX-512DQ, i.e. Cannonlake and beyond.]
	*********************************************************************************************************************/\
		"vmovaps	(%%rax),%%zmm1				\n\t	vmovaps	0x040(%%rax),%%zmm5				\n\t	vmovaps	0x080(%%rax),%%zmm9				\n\t	vmovaps	0x0c0(%%rax),%%zmm13			\n\t"/* load lo.lo and save copy in zmm2 */\
		"vmovaps	%%zmm1, %%zmm2				\n\t	vmovaps	%%zmm5, %%zmm6					\n\t	vmovaps	%%zmm9, %%zmm10					\n\t	vmovaps	%%zmm13, %%zmm14				\n\t"\
		/* qinv.lo*lo.lo: a = qinv.lo; b = lo.lo: */\
		"     vmulpd	(%%rdx),%%zmm1,%%zmm3	\n\t	     vmulpd	0x040(%%rdx),%%zmm5,%%zmm7	\n\t	     vmulpd	0x080(%%rdx),%%zmm9,%%zmm11	\n\t	     vmulpd	0x0c0(%%rdx),%%zmm13,%%zmm15\n\t"/* hi = fma(a,b, 0  ) - use VMULPD instead of VFMADD231PD, since addend = 0. */\
		"vmovaps	%%zmm3,%%zmm0				\n\t	vmovaps	%%zmm7,%%zmm4					\n\t	vmovaps	%%zmm11,%%zmm8					\n\t	vmovaps	%%zmm15,%%zmm12					\n\t"/* cpy hi into lo-reg */\
		"vfmsub231pd	(%%rdx),%%zmm1,%%zmm0	\n\t	vfmsub231pd	0x040(%%rdx),%%zmm5,%%zmm4	\n\t	vfmsub231pd	0x080(%%rdx),%%zmm9,%%zmm8	\n\t	vfmsub231pd	0x0c0(%%rdx),%%zmm13,%%zmm12\n\t"/* lo = fma(a,b, -hi) */\
		"vaddpd			%%zmm29,%%zmm3,%%zmm1	\n\t	vaddpd			%%zmm29,%%zmm7,%%zmm5	\n\t	vaddpd			%%zmm29,%%zmm11,%%zmm9	\n\t	vaddpd			%%zmm29,%%zmm15,%%zmm13	\n\t"\
		"vsubpd			%%zmm29,%%zmm1,%%zmm1	\n\t	vsubpd			%%zmm29,%%zmm5,%%zmm5	\n\t	vsubpd			%%zmm29,%%zmm9,%%zmm9	\n\t	vsubpd			%%zmm29,%%zmm13,%%zmm13	\n\t"/* hh = hi +- crnd50 to round-to-nearest-multiple-of-base */\
		"vsubpd			%%zmm1,%%zmm3,%%zmm3	\n\t	vsubpd			%%zmm5,%%zmm7,%%zmm7	\n\t	vsubpd			%%zmm9,%%zmm11,%%zmm11	\n\t	vsubpd			%%zmm13,%%zmm15,%%zmm15	\n\t"/* hi - hh gives backward-carry... */\
		"vaddpd			%%zmm3,%%zmm0,%%zmm0	\n\t	vaddpd			%%zmm7,%%zmm4,%%zmm4	\n\t	vaddpd			%%zmm11,%%zmm8,%%zmm8	\n\t	vaddpd			%%zmm15,%%zmm12,%%zmm12	\n\t"/* ...which we add to lo (prod[0]). */\
		"vmovaps	%%zmm0,     (%%rax)			\n\t	vmovaps	%%zmm4,0x040(%%rax)				\n\t	vmovaps	%%zmm8,0x080(%%rax)				\n\t	vmovaps	%%zmm12,0x0c0(%%rax)			\n\t"/* write prod[0]; this overwrites lo.lo, but copy of that saved in zmm2 */\
													/*** hi*base in zmm1, lo.lo in zmm2 ***/\
		/* qinv.hi*lo.lo: a = qinv.hi, b = lo.lo: */\
		"     vmulpd 0x100(%%rdx),%%zmm2,%%zmm3	\n\t	     vmulpd 0x140(%%rdx),%%zmm6,%%zmm7	\n\t	     vmulpd 0x180(%%rdx),%%zmm10,%%zmm11\n\t	     vmulpd 0x1c0(%%rdx),%%zmm14,%%zmm15\n\t"/* hi = fma(a,b, 0  ) - use VMULPD instead of VFMADD231PD, since addend = 0. */\
		"vmovaps	%%zmm3,%%zmm0				\n\t	vmovaps	%%zmm7,%%zmm4					\n\t	vmovaps	%%zmm11,%%zmm8					\n\t	vmovaps	%%zmm15,%%zmm12					\n\t"/* cpy hi into lo-reg */\
		"vfmsub231pd 0x100(%%rdx),%%zmm2,%%zmm0	\n\t	vfmsub231pd 0x140(%%rdx),%%zmm6,%%zmm4	\n\t	vfmsub231pd 0x180(%%rdx),%%zmm10,%%zmm8	\n\t	vfmsub231pd 0x1c0(%%rdx),%%zmm14,%%zmm12\n\t"/* lo = fma(a,b, -hi) */\
		/* NOTE: Cannot simply sum raw (unnormalized) hi-outouts here with one from FMA-mul below, because hi-terms generally
		have 53 significant bits but == 0 w.r.to different powers of 2, thus we lose critical bits on the low end in sch an add. */\
		"vaddpd			%%zmm29,%%zmm3,%%zmm2	\n\t	vaddpd			%%zmm29,%%zmm7,%%zmm6	\n\t	vaddpd			%%zmm29,%%zmm11,%%zmm10	\n\t	vaddpd			%%zmm29,%%zmm15,%%zmm14	\n\t"\
		"vsubpd			%%zmm29,%%zmm2,%%zmm2	\n\t	vsubpd			%%zmm29,%%zmm6,%%zmm6	\n\t	vsubpd			%%zmm29,%%zmm10,%%zmm10	\n\t	vsubpd			%%zmm29,%%zmm14,%%zmm14	\n\t"/* hh = hi +- crnd50 to round-to-nearest-multiple-of-base */\
		"vsubpd			%%zmm2,%%zmm3,%%zmm3	\n\t	vsubpd			%%zmm6,%%zmm7,%%zmm7	\n\t	vsubpd			%%zmm10,%%zmm11,%%zmm11	\n\t	vsubpd			%%zmm14,%%zmm15,%%zmm15	\n\t"/* hi - hh gives backward-carry... */\
		"vaddpd			%%zmm3,%%zmm0,%%zmm0	\n\t	vaddpd			%%zmm7,%%zmm4,%%zmm4	\n\t	vaddpd			%%zmm11,%%zmm8,%%zmm8	\n\t	vaddpd			%%zmm15,%%zmm12,%%zmm12	\n\t"/* ...which we add to lo. */\
		"vfmadd231pd	%%zmm30,%%zmm1,%%zmm0	\n\t	vfmadd231pd	%%zmm30,%%zmm5,%%zmm4		\n\t	vfmadd231pd	%%zmm30,%%zmm9,%%zmm8		\n\t	vfmadd231pd	%%zmm30,%%zmm13,%%zmm12		\n\t"/* Add lo to hi-output of previous lo:hi pair, which also needs *= binv */\
													/*** prod[1] (partial sum) in zmm0 ***/\
		/* qinv.lo*lo.hi: a = qinv.lo, b = lo.hi: */\
		"vmovaps	0x100(%%rax),%%zmm1			\n\t	vmovaps	0x140(%%rax),%%zmm5				\n\t	vmovaps	0x180(%%rax),%%zmm9				\n\t	vmovaps	0x1c0(%%rax),%%zmm13			\n\t"/* load lo.hi into zmm1... */\
		"     vmulpd	(%%rdx),%%zmm1,%%zmm3	\n\t	     vmulpd	0x040(%%rdx),%%zmm5,%%zmm7	\n\t	     vmulpd	0x080(%%rdx),%%zmm9,%%zmm11	\n\t	     vmulpd	0x0c0(%%rdx),%%zmm13,%%zmm15\n\t"/* hi = fma(a,b, 0  ) - use VMULPD instead of VFMADD231PD, since addend = 0. */\
		"vmovaps	%%zmm3,%%zmm2				\n\t	vmovaps	%%zmm7,%%zmm6					\n\t	vmovaps	%%zmm11,%%zmm10					\n\t	vmovaps	%%zmm15,%%zmm14					\n\t"/* cpy hi into lo-reg */\
		"vfmsub231pd	(%%rdx),%%zmm1,%%zmm2	\n\t	vfmsub231pd	0x040(%%rdx),%%zmm5,%%zmm6	\n\t	vfmsub231pd	0x080(%%rdx),%%zmm9,%%zmm10	\n\t	vfmsub231pd	0x0c0(%%rdx),%%zmm13,%%zmm14\n\t"/* lo = fma(a,b, -hi) */\
		"vaddpd			%%zmm29,%%zmm3,%%zmm1	\n\t	vaddpd			%%zmm29,%%zmm7,%%zmm5	\n\t	vaddpd			%%zmm29,%%zmm11,%%zmm9	\n\t	vaddpd			%%zmm29,%%zmm15,%%zmm13	\n\t"\
		"vsubpd			%%zmm29,%%zmm1,%%zmm1	\n\t	vsubpd			%%zmm29,%%zmm5,%%zmm5	\n\t	vsubpd			%%zmm29,%%zmm9,%%zmm9	\n\t	vsubpd			%%zmm29,%%zmm13,%%zmm13	\n\t"/* hh = hi +- crnd50 to round-to-nearest-multiple-of-base */\
		"vsubpd			%%zmm1,%%zmm3,%%zmm3	\n\t	vsubpd			%%zmm5,%%zmm7,%%zmm7	\n\t	vsubpd			%%zmm9,%%zmm11,%%zmm11	\n\t	vsubpd			%%zmm13,%%zmm15,%%zmm15	\n\t"/* hi - hh gives backward-carry... */\
		"vaddpd			%%zmm3,%%zmm2,%%zmm2	\n\t	vaddpd			%%zmm7,%%zmm6,%%zmm6	\n\t	vaddpd			%%zmm11,%%zmm10,%%zmm10	\n\t	vaddpd			%%zmm15,%%zmm14,%%zmm14	\n\t"/* ...which we add to lo... */\
		"vaddpd			%%zmm0,%%zmm2,%%zmm0	\n\t	vaddpd			%%zmm4,%%zmm6,%%zmm4	\n\t	vaddpd			%%zmm8,%%zmm10,%%zmm8	\n\t	vaddpd			%%zmm12,%%zmm14,%%zmm12	\n\t"/* ...and add that to zmm0 (prod[1]). */\
		/* Since x = prod[1] was obtained by summing two base-normalized FMA-mul low-half outputs, need to still normalize it: */\
		"     vmulpd	%%zmm30,%%zmm0,%%zmm1	\n\t	     vmulpd	%%zmm30,%%zmm4,%%zmm5		\n\t	     vmulpd	%%zmm30,%%zmm8,%%zmm9		\n\t	     vmulpd	%%zmm30,%%zmm12,%%zmm13		\n\t"/* x*binv */\
		"vrndscalepd $1,%%zmm1,%%zmm1			\n\t	vrndscalepd $1,%%zmm5,%%zmm5			\n\t	vrndscalepd $1,%%zmm9,%%zmm9			\n\t	vrndscalepd $1,%%zmm13,%%zmm13			\n\t"/* cy = floor(x*binv) */\
		"vfnmadd231pd	%%zmm31,%%zmm1,%%zmm0	\n\t	vfnmadd231pd	%%zmm31,%%zmm5,%%zmm4	\n\t	vfnmadd231pd	%%zmm31,%%zmm9,%%zmm8	\n\t	vfnmadd231pd	%%zmm31,%%zmm13,%%zmm12	\n\t"/* x = FMA(cy,-base, x); x in [0,base) */\
		"vmovaps	%%zmm0,0x100(%%rax)			\n\t	vmovaps	%%zmm4,0x140(%%rax)				\n\t	vmovaps	%%zmm8,0x180(%%rax)				\n\t	vmovaps	%%zmm12,0x1c0(%%rax)			\n\t"/* write prod[1] (lo.hi) */\
													/*** lo.lo in (%%rax), lo.hi in zmm0. ***/\
	/********************* lo = UMULH(q,lo) : *****************************************************************************
		Bitranges (B = base):													Subproducts:
											|----------- 2B-1: 0 -----------|	q.lo * lo.lo
							|----------- 3B-1: B -----------|					q.hi * lo.lo, q.lo * lo.hi
			|----------- 4B-1:2B -----------|									q.hi * lo.hi
			|----------- OUTPUT ------------|
		Thus may only need high half of the (q.lo * lo.lo) subproducts - alas, I see few other work-saving possibilities.
		But wait - as long as our base is a few bits less than the 53-bits of an IEEE double mantissa, the hi-part MUL of
		a double-wide FMA-mul with its 53 significant bits will capture a disproportionate share of the full-prodict bits.
		Thus it my be possible to simplify the computation of the lower-half carryin-bits by computing the hi bits of the
		(q.lo * lo.lo) subproduct, feeding that result * binv as the addend to the (q.hi * lo.lo) subproduct (again sans
		a corresponding low-part FMA-mul), and adding that result in turn to the high-part-only (q.lo * lo.hi) subproduct.
		We can then round off the excess bits at the low end of the result, and proceed to compute the exact double-wide
		(q.hi * lo.hi) subproduct.
	*********************************************************************************************************************/\
		"vmulpd		(%%rax),%%zmm30,%%zmm1		\n\t	vmulpd	0x040(%%rax),%%zmm30,%%zmm5		\n\t	vmulpd	0x080(%%rax),%%zmm30,%%zmm9		\n\t	vmulpd	0x0c0(%%rax),%%zmm30,%%zmm13	\n\t"/* Premultiply lo.lo by binv */\
		"vmovaps	(%%rax),%%zmm2				\n\t	vmovaps	0x040(%%rax),%%zmm6				\n\t	vmovaps	0x080(%%rax),%%zmm10			\n\t	vmovaps	0x0c0(%%rax),%%zmm14			\n\t"/* load lo.lo and save copy in zmm2 */\
		"vmovaps	0x100(%%rax),%%zmm3			\n\t	vmovaps	0x140(%%rax),%%zmm7				\n\t	vmovaps	0x180(%%rax),%%zmm11			\n\t	vmovaps	0x1c0(%%rax),%%zmm15			\n\t"/* preload lo.hi */\
		/* q.lo*lo.lo: a = q.lo; b = lo.lo: */\
		"vmulpd		(%%rbx),%%zmm1,%%zmm1		\n\t	vmulpd		0x040(%%rbx),%%zmm5,%%zmm5	\n\t	vmulpd		0x080(%%rbx),%%zmm9,%%zmm9	\n\t	vmulpd		0x0c0(%%rbx),%%zmm13,%%zmm13\n\t"/* [high 53 bits of q.lo * lo.lo]*binv... */\
		"vfmadd231pd 0x100(%%rbx),%%zmm2,%%zmm1	\n\t	vfmadd231pd 0x140(%%rbx),%%zmm6,%%zmm5	\n\t	vfmadd231pd 0x180(%%rbx),%%zmm10,%%zmm9	\n\t	vfmadd231pd 0x1c0(%%rbx),%%zmm14,%%zmm13\n\t"/* ...feed that as addend to (q.hi * lo.lo) subproduct (again high-53-bits only)... */\
		"vfmadd231pd      (%%rbx),%%zmm3,%%zmm1	\n\t	vfmadd231pd 0x040(%%rbx),%%zmm7,%%zmm5	\n\t	vfmadd231pd 0x080(%%rbx),%%zmm11,%%zmm9	\n\t	vfmadd231pd 0x0c0(%%rbx),%%zmm15,%%zmm13\n\t"/* ...add that result in turn to the high-part-only (q.lo * lo.hi) subproduct. */\
		"vmulpd		%%zmm30,%%zmm1,%%zmm0		\n\t	vmulpd		%%zmm30,%%zmm5,%%zmm4		\n\t	vmulpd		%%zmm30,%%zmm9,%%zmm8		\n\t	vmulpd		%%zmm30,%%zmm13,%%zmm12		\n\t"/* tmp *= binv */\
	"vrndscalepd $3,%%zmm0,%%zmm0				\n\t	vrndscalepd $3,%%zmm4,%%zmm4			\n\t	vrndscalepd $3,%%zmm8,%%zmm8			\n\t	vrndscalepd $3,%%zmm12,%%zmm12			\n\t"/* round_toward_0(tmp*binv), this gets added to low half of (q.hi * lo.hi) */\
		/* q.hi*lo.hi: a = q.hi, b = lo.hi: */\
		"     vmulpd 0x100(%%rbx),%%zmm3,%%zmm1	\n\t	     vmulpd 0x140(%%rbx),%%zmm7,%%zmm5	\n\t	     vmulpd 0x180(%%rbx),%%zmm11,%%zmm9	\n\t	     vmulpd 0x1c0(%%rbx),%%zmm15,%%zmm13\n\t"/* hi = fma(a,b, 0  ) - use VMULPD instead of VFMADD231PD, since addend = 0. */\
		"vcmppd	$14,0x400(%%rax),%%zmm1,%%k1	\n\t	vcmppd	$14,0x440(%%rax),%%zmm5,%%k2	\n\t	vcmppd	$14,0x480(%%rax),%%zmm9,%%k3	\n\t	vcmppd	$14,0x4c0(%%rax),%%zmm13,%%k4	\n\t"/* bitmask = (lo53 > hi53)? In AVX512 version, mask-regs k1-4 replace dest-regs m2,6,10,14 */\
		"vmovaps	%%zmm1,%%zmm2				\n\t	vmovaps	%%zmm5,%%zmm6					\n\t	vmovaps	%%zmm9,%%zmm10					\n\t	vmovaps	%%zmm13,%%zmm14					\n\t"/* cpy hi into lo-reg */\
		"vfmsub231pd 0x100(%%rbx),%%zmm3,%%zmm2	\n\t	vfmsub231pd 0x140(%%rbx),%%zmm7,%%zmm6	\n\t	vfmsub231pd 0x180(%%rbx),%%zmm11,%%zmm10\n\t	vfmsub231pd 0x1c0(%%rbx),%%zmm15,%%zmm14\n\t"/* lo = fma(a,b, -hi) */\
		"vaddpd			%%zmm29,%%zmm1,%%zmm3	\n\t	vaddpd			%%zmm29,%%zmm5,%%zmm7	\n\t	vaddpd			%%zmm29,%%zmm9,%%zmm11	\n\t	vaddpd			%%zmm29,%%zmm13,%%zmm15	\n\t"\
		"vsubpd			%%zmm29,%%zmm3,%%zmm3	\n\t	vsubpd			%%zmm29,%%zmm7,%%zmm7	\n\t	vsubpd			%%zmm29,%%zmm11,%%zmm11	\n\t	vsubpd			%%zmm29,%%zmm15,%%zmm15	\n\t"/* hh = hi +- crnd50 to round-to-nearest-multiple-of-base */\
		"vsubpd			%%zmm3,%%zmm1,%%zmm1	\n\t	vsubpd			%%zmm7,%%zmm5,%%zmm5	\n\t	vsubpd			%%zmm11,%%zmm9,%%zmm9	\n\t	vsubpd			%%zmm15,%%zmm13,%%zmm13	\n\t"/* hi - hh gives backward-carry... */\
		"vmulpd		%%zmm30,%%zmm3,%%zmm3		\n\t	vmulpd		%%zmm30,%%zmm7,%%zmm7		\n\t	vmulpd		%%zmm30,%%zmm11,%%zmm11		\n\t	vmulpd		%%zmm30,%%zmm15,%%zmm15		\n\t"/* hh *= binv */\
		"vaddpd			%%zmm1,%%zmm2,%%zmm2	\n\t	vaddpd			%%zmm5,%%zmm6,%%zmm6	\n\t	vaddpd			%%zmm9,%%zmm10,%%zmm10	\n\t	vaddpd			%%zmm13,%%zmm14,%%zmm14	\n\t"/* ...which we add to lo... */\
		"vaddpd			%%zmm0,%%zmm2,%%zmm0	\n\t	vaddpd			%%zmm4,%%zmm6,%%zmm4	\n\t	vaddpd			%%zmm8,%%zmm10,%%zmm8	\n\t	vaddpd			%%zmm12,%%zmm14,%%zmm12	\n\t"/* ...and add that to zmm0. */\
													/*** lo.lo in zmm0, lo.hi in zmm3 ***/\
		/* If h < l (already done above using leading significant 53 bits of both operands, i.e. the unnormalized
		high outputs of the respective 2-word FMA-muls), calculate h-l+q; otherwise h-l. Result is in [0, q). */\
		"vmovaps	0x200(%%rax),%%zmm1			\n\t	vmovaps	0x240(%%rax),%%zmm5				\n\t	vmovaps	0x280(%%rax),%%zmm9				\n\t	vmovaps	0x2c0(%%rax),%%zmm13			\n\t"/* hi.lo */\
		"vmovaps	0x300(%%rax),%%zmm2			\n\t	vmovaps	0x340(%%rax),%%zmm6				\n\t	vmovaps	0x380(%%rax),%%zmm10			\n\t	vmovaps	0x3c0(%%rax),%%zmm14			\n\t"/* hi.hi */\
		"vsubpd	%%zmm0,%%zmm1,%%zmm1			\n\t	vsubpd	%%zmm4,%%zmm5,%%zmm5			\n\t	vsubpd	%%zmm8,%%zmm9,%%zmm9			\n\t	vsubpd	%%zmm12,%%zmm13,%%zmm13			\n\t"/* (hi - lo).lo, zmm0 FREE */\
		"vsubpd	%%zmm3,%%zmm2,%%zmm2			\n\t	vsubpd	%%zmm7,%%zmm6,%%zmm6			\n\t	vsubpd	%%zmm11,%%zmm10,%%zmm10			\n\t	vsubpd	%%zmm15,%%zmm14,%%zmm14			\n\t"/* (hi - lo).hi, zmm3 FREE */\
		"vmovaps	0x100(%%rbx),%%zmm0			\n\t	vmovaps	0x140(%%rbx),%%zmm4				\n\t	vmovaps	0x180(%%rbx),%%zmm8				\n\t	vmovaps	0x1c0(%%rbx),%%zmm12			\n\t"/* q.hi */\
	   "vaddpd (%%rbx),%%zmm1,%%zmm1%{%%k1%}	\n\t vaddpd	0x040(%%rbx),%%zmm5,%%zmm5%{%%k2%}	\n\t vaddpd	0x080(%%rbx),%%zmm9,%%zmm9%{%%k3%}	\n\t vaddpd	0x0c0(%%rbx),%%zmm13,%%zmm13%{%%k4%}\n\t"/* xlo = (hi-lo).lo + (q.lo & bitmask) */\
	   "vaddpd	%%zmm0,%%zmm2,%%zmm2%{%%k1%}	\n\t	vaddpd	%%zmm4,%%zmm6,%%zmm6%{%%k2%}	\n\t	vaddpd	%%zmm8,%%zmm10,%%zmm10%{%%k3%}	\n\t	vaddpd	%%zmm12,%%zmm14,%%zmm14%{%%k4%}	\n\t"/* xhi = (hi-lo).hi + (q.hi & bitmask) */\
													/*** x.lo,hi in zmm1,zmm2; q.hi in zmm0 ***/\
	/* if((pshift >> j) & (uint64)1) { */\
		"movq	%[__pshift],%%rsi				\n\t"\
		"movslq	%[__j],%%rcx					\n\t"\
		"shrq	%%cl,%%rsi						\n\t"\
		"andq	$0x1,%%rsi						\n\t"\
	"je twopmodq100_2wdq32_gcc64				\n\t"\
		"vaddpd	%%zmm1,%%zmm1,%%zmm1			\n\t	vaddpd	%%zmm5,%%zmm5,%%zmm5			\n\t	vaddpd	%%zmm9,%%zmm9,%%zmm9			\n\t	vaddpd	%%zmm13,%%zmm13,%%zmm13			\n\t"/* 2*x.lo */\
		"vaddpd	%%zmm2,%%zmm2,%%zmm2			\n\t	vaddpd	%%zmm6,%%zmm6,%%zmm6			\n\t	vaddpd	%%zmm10,%%zmm10,%%zmm10			\n\t	vaddpd	%%zmm14,%%zmm14,%%zmm14			\n\t"/* 2*x.hi */\
		"vmovaps	%%zmm1,%%zmm3				\n\t	vmovaps	%%zmm5,%%zmm7					\n\t	vmovaps	%%zmm9,%%zmm11					\n\t	vmovaps	%%zmm13,%%zmm15					\n\t"/* cpy of 2*x.lo */\
		"vfmadd231pd %%zmm31,%%zmm2,%%zmm3		\n\t	vfmadd231pd %%zmm31,%%zmm6,%%zmm7		\n\t	vfmadd231pd %%zmm31,%%zmm10,%%zmm11		\n\t	vfmadd231pd %%zmm31,%%zmm14,%%zmm15		\n\t"/* 2*(x.lo + x.hi*base) to get most-sig 53 bits of doubled-x for the following comparison */\
		/* If x > q, subtract q: */\
		"vcmppd $13,0x200(%%rbx),%%zmm3,%%k1	\n\t	vcmppd $13,0x240(%%rbx),%%zmm7,%%k2		\n\t	vcmppd $13,0x280(%%rbx),%%zmm11,%%k3	\n\t	vcmppd $13,0x2c0(%%rbx),%%zmm15,%%k4	\n\t"/* bitmask = (xhi53 >= qhi53)? In AVX512 version, mask-regs k1-4 replace dest-regs m3,7,11,15 */\
	   "vsubpd (%%rbx),%%zmm1,%%zmm1%{%%k1%}	\n\t vsubpd	0x040(%%rbx),%%zmm5,%%zmm5%{%%k2%}	\n\t vsubpd	0x080(%%rbx),%%zmm9,%%zmm9%{%%k3%}	\n\t vsubpd	0x0c0(%%rbx),%%zmm13,%%zmm13%{%%k4%}\n\t"/* xlo = 2*x.lo - (q.lo & bitmask) */\
	   "vsubpd	%%zmm0,%%zmm2,%%zmm2%{%%k1%}	\n\t	vsubpd	%%zmm4,%%zmm6,%%zmm6%{%%k2%}	\n\t	vsubpd	%%zmm8,%%zmm10,%%zmm10%{%%k3%}	\n\t	vsubpd	%%zmm12,%%zmm14,%%zmm14%{%%k4%}	\n\t"/* xhi = 2*x.hi - (q.hi & bitmask) */\
	"twopmodq100_2wdq32_gcc64:					\n\t"\
		"vmovaps	%%zmm1,     (%%rax)			\n\t	vmovaps	%%zmm5,0x040(%%rax)				\n\t	vmovaps	%%zmm9,0x080(%%rax)				\n\t	vmovaps	%%zmm13,0x0c0(%%rax)			\n\t"/* write x.lo */\
		"vmovaps	%%zmm2,0x100(%%rax)			\n\t	vmovaps	%%zmm6,0x140(%%rax)				\n\t	vmovaps	%%zmm10,0x180(%%rax)			\n\t	vmovaps	%%zmm14,0x1c0(%%rax)			\n\t"/* write x.hi */\
	/* } */\
		:					/* outputs: none */\
		: [__aq0]	 "m" (Xaq0)	/* All inputs from memory addresses here */\
		 ,[__aqinv0] "m" (Xaqinv0)	\
		 ,[__ax0]	 "m" (Xax0)		\
		 ,[__base] "m" (Xbase)	\
		 ,[__binv] "m" (Xbinv)	\
		 ,[__crnd] "m" (Xcrnd)\
		 ,[__pshift] "m" (Xpshift)	\
		 ,[__j]		 "m" (Xj)		\
		: "cc","memory","cl","rax","rbx","rcx","rdx","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15", "xmm28","xmm29","xmm30","xmm31"	/* Clobbered registers */\
	);\
	}

  #endif	// ALL_FMA ?

#endif	// USE_AVX512 ?

#ifdef __cplusplus
}
#endif

#endif	/* twopmodq100_h_included */

