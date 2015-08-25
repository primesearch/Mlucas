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
#ifndef radix32_dif_dit_pass_gcc_h_included
#define radix32_dif_dit_pass_gcc_h_included

#ifdef USE_AVX2	// FMA-based versions of selected macros in this file for Intel AVX2/FMA3

  #define USE_64BIT_ASM_STYLE	1	// Only support 64-bit-style version using xmm0-15 in AVX2 mode and above.

	// Oct 2014: Aggressively deploy FMA, both to save steps in CMULs and in-place butterflies
	// (e.g. in the frequent non-FMA sequence x = x-y, y *= 2, y += x which produces x +- y in the
	// original y,x-regs, we replace the latter 2 ops with one FMA: y = 2*y + x), and also to
	// replace all the ADD/SUBs, thus trading the lower latency of those (3 cycles vs 5 for FMA)
	// with the higher throughput of FMA (2 per cycle vs just 1 ADD/SUB). In order to distinguish
	// the former FMAs from the latter "trivial" ones (in the sense that one multiplicand is unity),
	// we 'undent' the former one tab leftward relative to the trivial FMAs.
	//
	// For the ADD/SUB -> FMA we use rsi to hold the SIMD-register-width-propagated 1.0 and then replace like so:
	// ADD:
	//		vaddpd		%%ymm2,%%ymm1,%%ymm1		// ymm1 += ymm2
	//	-->	vfmadd132pd	(%%rsi),%%ymm2,%%ymm1		// ymm1 = ymm1*1.0 + ymm2
	// SUB:
	//		vsubpd		%%ymm2,%%ymm1,%%ymm1		// ymm1 -= ymm2
	//	-->	vfmsub132pd	(%%rsi),%%ymm2,%%ymm1		// ymm1 = ymm1*1.0 - ymm2
	//
	// The choice of the 132-variant of Intel FMA3 preserves the name (add or sub) of the replaced ADD/SUB op.
	// If we have a register free (or have > 2 such ops in a row, i.e. can spill/reload a data-reg with less
	// memory traffic than multiple implied-loads-of-1.0 would incur) we stick the 1.0 into a register.
	// The other kind of ADD/SUB we may wish to replace with FMA is one where one addend (subtrahend) is in mem,
	// in which case we can only FMAize if 1.0 is in-reg (call it ymm#), in which case we then replace like so:
	// ADD:
	//		vaddpd		(%%mem),%%ymm1,%%ymm1		// ymm1 += (mem)
	//	-->	vfmadd231pd	(%%mem),%%ymm#,%%ymm1		// ymm1 = +(mem)*1.0 + ymm1
	// SUB:
	//		vsubpd		%%ymm2,%%ymm1,%%ymm1		// ymm1 -= (mem)
	//	-->vfnmadd231pd	(%%mem),%%ymm#,%%ymm1		// ymm1 = -(mem)*1.0 + ymm1
	//
	/*
	For GCC-macro version of this, use that - in terms of vec_dbl pointer-offsets - isrt2 + 1,3,5 = cc0,cc1,cc3,
	isrt2 + 7,15,23,31,39,47,55,63 = c00,04,02,06,01,05,03,07, and isrt2 + 71,72 = one,two
	in order to reduce number of args to <= the GCC-allowed maximum of 30:
	*/
  #ifdef ALL_FMA	// Aggressive-FMA version: Replace most ADD/SUB by FMA-with-one-unity-multiplicand

	// In our FMA-ized version, original mix of [223 ADD, 223 SUB, 326 MUL] ==> [21 ADD, 32 SUB, 385 FMA (224 nontrivial), 102 MUL].
	//
	#define SSE2_RADIX32_DIF_TWIDDLE(Xadd0,Xp01,Xp02,Xp03,Xp04,Xp08,Xp0C,Xp10,Xp18,Xr00)\
	{\
	__asm__ volatile (\
	/*...Block 1: */\
		"movq	%[__add0],%%rax							\n\t"\
		"movslq	%[__p08],%%rbx							\n\t"\
		"movslq	%[__p10],%%rcx							\n\t	movslq	%[__p04],%%r9	\n\t"\
		"movslq	%[__p18],%%rdx							\n\t	movq	%[__r00],%%rsi	\n\t"\
		"shlq	$3,%%rbx								\n\t	shlq	$3,%%r9		\n\t"\
		"shlq	$3,%%rcx								\n\t	leaq	0x1100(%%rsi),%%r8	\n\t"/* two */\
		"shlq	$3,%%rdx								\n\t	leaq	(%%r9,%%rax),%%r10	\n\t"\
		"addq	%%rax,%%rbx								\n\t	leaq	(%%r9,%%rbx),%%r11	\n\t"\
		"addq	%%rax,%%rcx								\n\t	leaq	(%%r9,%%rcx),%%r12	\n\t"\
		"addq	%%rax,%%rdx								\n\t	leaq	(%%r9,%%rdx),%%r13	\n\t"\
		/* SSE2_RADIX4_DIF_4TWIDDLE_A(r00,c10) */				/* SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO_B(r08,c04) */\
		"vmovaps	     (%%rax),%%ymm0					\n\t	vmovaps	     (%%r10),%%ymm8 	\n\t"\
		"vmovaps	     (%%rcx),%%ymm4					\n\t	vmovaps	     (%%r12),%%ymm12	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1					\n\t	vmovaps	0x020(%%r10),%%ymm9 	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5					\n\t	vmovaps	0x020(%%r12),%%ymm13	\n\t"\
		"vmovaps	0x920(%%rsi),%%ymm2	/* c10 */		\n\t	vmovaps	0x9e0(%%rsi),%%ymm14	/* c04 */	\n\t"\
		"vmovaps	0x940(%%rsi),%%ymm3					\n\t	vmovaps	0xa00(%%rsi),%%ymm15	\n\t"\
		"vmovaps	%%ymm4,%%ymm6						\n\t	vmovaps	%%ymm8 ,%%ymm10	\n\t"\
		"vmovaps	%%ymm5,%%ymm7						\n\t	vmovaps	%%ymm9 ,%%ymm11	\n\t"\
		"vmulpd	%%ymm2,%%ymm4,%%ymm4					\n\t	vmulpd	%%ymm14,%%ymm8 ,%%ymm8 	\n\t"\
		"vmulpd	%%ymm2,%%ymm5,%%ymm5					\n\t	vmulpd	%%ymm14,%%ymm9 ,%%ymm9 	\n\t"\
	"vfnmadd231pd	%%ymm3,%%ymm7,%%ymm4			\n\t	vfnmadd231pd	%%ymm15,%%ymm11,%%ymm8 	\n\t"\
	" vfmadd231pd	%%ymm3,%%ymm6,%%ymm5			\n\t	 vfmadd231pd	%%ymm15,%%ymm10,%%ymm9 	\n\t"\
	"vmovaps	-0x20(%%r8),%%ymm6 	\n\t"/* ymm6 free; use for 1.0 */\
		"vmovaps	%%ymm0,%%ymm2						\n\t	vmovaps	%%ymm12,%%ymm14	\n\t"\
		"vmovaps	%%ymm1,%%ymm3						\n\t	vmovaps	%%ymm13,%%ymm15	\n\t"\
		"vfmadd132pd	%%ymm6 ,%%ymm4,%%ymm0			\n\t	vmulpd	0xa20(%%rsi),%%ymm12,%%ymm12	/* c14 */	\n\t"\
		"vfmadd132pd	%%ymm6 ,%%ymm5,%%ymm1			\n\t	vmulpd	0xa20(%%rsi),%%ymm13,%%ymm13	\n\t"\
		"vfmsub132pd	%%ymm6 ,%%ymm4,%%ymm2		\n\t	vfnmadd231pd	0xa40(%%rsi),%%ymm15,%%ymm12	\n\t"\
		"vfmsub132pd	%%ymm6 ,%%ymm5,%%ymm3		\n\t	 vfmadd231pd	0xa40(%%rsi),%%ymm14,%%ymm13	\n\t"\
	"vmovaps	%%ymm6,%%ymm14 	\n\t"/* ymm14 free; use for 1.0 */\
		"vmovaps	     (%%rdx),%%ymm4					\n\t	vmovaps	%%ymm8 ,%%ymm10	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm5					\n\t	vmovaps	%%ymm9 ,%%ymm11	\n\t"\
		"vmovaps	     (%%rdx),%%ymm6					\n\t	vfmadd132pd	%%ymm14,%%ymm12,%%ymm8 	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm7					\n\t	vfmadd132pd	%%ymm14,%%ymm13,%%ymm9 	\n\t"\
		"vmulpd	0x9a0(%%rsi),%%ymm4,%%ymm4	/* c18 */	\n\t	vfmsub132pd	%%ymm14,%%ymm12,%%ymm10	\n\t"\
		"vmulpd	0x9a0(%%rsi),%%ymm5,%%ymm5				\n\t	vfmsub132pd	%%ymm14,%%ymm13,%%ymm11	\n\t"\
	"vfnmadd231pd	0x9c0(%%rsi),%%ymm7,%%ymm4			\n\t	vmovaps	     (%%r13),%%ymm12	\n\t"\
	" vfmadd231pd	0x9c0(%%rsi),%%ymm6,%%ymm5			\n\t	vmovaps	0x020(%%r13),%%ymm13	\n\t"\
		"vmovaps	%%ymm5,0x020(%%rsi)					\n\t	vmovaps	     (%%r13),%%ymm14	\n\t"\
		"vmovaps	%%ymm4,     (%%rsi)	/* tmpstr r00 */\n\t	vmovaps	0x020(%%r13),%%ymm15	\n\t"\
		"vmovaps	     (%%rbx),%%ymm4					\n\t	vmulpd	0xaa0(%%rsi),%%ymm12,%%ymm12	/* c1C */	\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm5					\n\t	vmulpd	0xaa0(%%rsi),%%ymm13,%%ymm13	\n\t"\
		"vmovaps	     (%%rbx),%%ymm6				\n\t	vfnmadd231pd	0xac0(%%rsi),%%ymm15,%%ymm12	\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm7				\n\t	 vfmadd231pd	0xac0(%%rsi),%%ymm14,%%ymm13	\n\t"\
		"vmulpd	0x960(%%rsi),%%ymm4,%%ymm4	/* c08 */	\n\t	vmovaps	%%ymm13,0x120(%%rsi)	\n\t"\
		"vmulpd	0x960(%%rsi),%%ymm5,%%ymm5				\n\t	vmovaps	%%ymm12,0x100(%%rsi)	/* tmpstr r08 */\n\t"\
	"vfnmadd231pd	0x980(%%rsi),%%ymm7,%%ymm4			\n\t	vmovaps	     (%%r11),%%ymm12	\n\t"\
	" vfmadd231pd	0x980(%%rsi),%%ymm6,%%ymm5			\n\t	vmovaps	0x020(%%r11),%%ymm13	\n\t"\
		"vmovaps	%%ymm5,%%ymm7						\n\t	vmovaps	     (%%r11),%%ymm14	\n\t"\
		"vmovaps	%%ymm4,%%ymm6						\n\t	vmovaps	0x020(%%r11),%%ymm15	\n\t"\
		"vsubpd	     (%%rsi),%%ymm4,%%ymm4	/* r00 */	\n\t	vmulpd	0xa60(%%rsi),%%ymm12,%%ymm12	/* c0C */	\n\t"\
		"vsubpd	0x020(%%rsi),%%ymm5,%%ymm5				\n\t	vmulpd	0xa60(%%rsi),%%ymm13,%%ymm13	\n\t"\
		"vaddpd	     (%%rsi),%%ymm6,%%ymm6			\n\t	vfnmadd231pd	0xa80(%%rsi),%%ymm15,%%ymm12	\n\t"\
		"vaddpd	0x020(%%rsi),%%ymm7,%%ymm7			\n\t	 vfmadd231pd	0xa80(%%rsi),%%ymm14,%%ymm13	\n\t"\
	"vmovaps	%%ymm3,(%%rsi)	\n\t	vmovaps	-0x20(%%r8),%%ymm3	\n\t"/* spill to make room for 1.0 */\
		"vfmsub132pd	%%ymm3 ,%%ymm6,%%ymm0			\n\t	vmovaps	%%ymm13,%%ymm15	\n\t"\
		"vfmsub132pd	%%ymm3 ,%%ymm5,%%ymm2			\n\t	vmovaps	%%ymm12,%%ymm14	\n\t"\
		"vfmsub132pd	%%ymm3 ,%%ymm7,%%ymm1			\n\t	vsubpd	0x100(%%rsi),%%ymm12,%%ymm12	\n\t"\
		"vfmsub132pd	(%%rsi),%%ymm4,%%ymm3			\n\t	vsubpd	0x120(%%rsi),%%ymm13,%%ymm13	\n\t"\
		"												\n\t	vaddpd	0x100(%%rsi),%%ymm14,%%ymm14	\n\t"\
		"												\n\t	vaddpd	0x120(%%rsi),%%ymm15,%%ymm15	\n\t"\
		"												\n\t	vsubpd	     %%ymm14,%%ymm8 ,%%ymm8 	\n\t"\
		"												\n\t	vsubpd	     %%ymm15,%%ymm9 ,%%ymm9 	\n\t"\
		"												\n\t	vfmsub132pd	-0x20(%%r8),%%ymm13,%%ymm10	\n\t"\
		"vmovaps	%%ymm8 ,0x180(%%rsi)				\n\t	vfmsub132pd	-0x20(%%r8),%%ymm12,%%ymm11	\n\t"\
	"vmovaps	%%ymm12,(%%rsi) 	\n\t"/* spill ymm12 to make room for 2.0 */"	vmovaps	(%%r8),%%ymm12 	\n\t"/* two */\
	"vfmadd132pd	%%ymm12,%%ymm0,%%ymm6			\n\t	vfmadd132pd	%%ymm12,%%ymm8 ,%%ymm14	\n\t"\
	"vfmadd132pd	%%ymm12,%%ymm2,%%ymm5			\n\t	vfmadd132pd	%%ymm12,%%ymm10,%%ymm13	\n\t"\
	"vfmadd132pd	%%ymm12,%%ymm1,%%ymm7			\n\t	vfmadd132pd	%%ymm12,%%ymm9 ,%%ymm15	\n\t"\
	"vfmadd132pd	%%ymm12,%%ymm3,%%ymm4			\n\t	vfmadd132pd	(%%rsi),%%ymm11,%%ymm12	\n\t"\
		/* SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r00) */\
		"												\n\t	vmovaps	%%ymm14,0x100(%%rsi)	\n\t"\
		"												\n\t	vmovaps	%%ymm15,0x120(%%rsi)	\n\t"\
		"												\n\t	vmovaps	%%ymm10,%%ymm14	\n\t"\
		"												\n\t	vmovaps	%%ymm13,%%ymm15	\n\t"\
		"												\n\t	vmovaps	-0x20(%%r8 ),%%ymm8 	\n\t"/* one */\
		"														vfmsub132pd	%%ymm8 ,%%ymm12,%%ymm10	\n\t"\
		"														vfmsub132pd	%%ymm8 ,%%ymm11,%%ymm13	\n\t"\
		"														vfmadd132pd	%%ymm8 ,%%ymm12,%%ymm14	\n\t"\
		"														vfmadd132pd	%%ymm8 ,%%ymm11,%%ymm15	\n\t"\
		"												\n\t	vmovaps	0x800(%%rsi),%%ymm8 	\n\t"/* isrt2 */\
		"vmovaps	0x100(%%rsi),%%ymm11			\n\t	vfnmadd231pd %%ymm8,%%ymm10,%%ymm2	\n\t"\
		"vmovaps	0x120(%%rsi),%%ymm12			\n\t	vfnmadd231pd %%ymm8,%%ymm13,%%ymm3	\n\t"\
		"vfmsub132pd	-0x20(%%r8),%%ymm11,%%ymm6		\n\t	vfnmadd231pd %%ymm8,%%ymm14,%%ymm4	\n\t"\
		"vfmsub132pd	-0x20(%%r8),%%ymm9 ,%%ymm0		\n\t	vfnmadd231pd %%ymm8,%%ymm15,%%ymm5	\n\t"\
		"vsubpd	0x180(%%rsi),%%ymm1 ,%%ymm1 			\n\t	vfmsub132pd	-0x20(%%r8),%%ymm12 ,%%ymm7	\n\t"\
	"vmovaps	(%%r8),%%ymm8	\n\t"/* rcol: *2 => *sqrt2 due to above FMA leaving ymm10,13-15 unmultiplied-by-isrt2: */\
		"vmovaps	%%ymm6 ,0x100(%%rsi)				\n\t	vmovaps	%%ymm2 ,0x140(%%rsi)	\n\t"\
		"vmovaps	%%ymm0 ,0x080(%%rsi)				\n\t	vmovaps	%%ymm5 ,0x0c0(%%rsi)	\n\t"\
		"vmovaps	%%ymm7 ,0x120(%%rsi)				\n\t	vmovaps	%%ymm4 ,0x160(%%rsi)	\n\t"\
		"vmovaps	%%ymm1 ,0x1a0(%%rsi)				\n\t	vmovaps	%%ymm3 ,0x1e0(%%rsi)	\n\t"\
	"vfmadd132pd	%%ymm8 ,%%ymm6 ,%%ymm11			\n\t	vfmadd132pd	0x20(%%r8),%%ymm2 ,%%ymm10	\n\t"\
	"vfmadd132pd	%%ymm8 ,%%ymm0 ,%%ymm9 			\n\t	vfmadd132pd	0x20(%%r8),%%ymm5 ,%%ymm15	\n\t"\
	"vfmadd132pd	%%ymm8 ,%%ymm7 ,%%ymm12			\n\t	vfmadd132pd	0x20(%%r8),%%ymm4 ,%%ymm14	\n\t"\
	"vfmadd132pd	0x180(%%rsi),%%ymm1 ,%%ymm8 	\n\t	vfmadd132pd	0x20(%%r8),%%ymm3 ,%%ymm13	\n\t"\
		"vmovaps	%%ymm11,     (%%rsi)				\n\t	vmovaps	%%ymm10,0x040(%%rsi)	\n\t"\
		"vmovaps	%%ymm9 ,0x180(%%rsi)				\n\t	vmovaps	%%ymm15,0x1c0(%%rsi)	\n\t"\
		"vmovaps	%%ymm12,0x020(%%rsi)				\n\t	vmovaps	%%ymm14,0x060(%%rsi)	\n\t"\
		"vmovaps	%%ymm8 ,0x0a0(%%rsi)				\n\t	vmovaps	%%ymm13,0x0e0(%%rsi)	\n\t"\
		"\n\t"\
	/***************************************/\
		"\n\t"\
	/*...Block 2: */\
		"movslq	%[__p02],%%rdi							\n\t"\
		"shlq	$3,%%rdi								\n\t"	/* p04<<3 still in r9 */\
		"addq	%%rdi,%%rax	/* &a[j1+p2] */				\n\t	leaq	(%%r9,%%rax),%%r10	\n\t"\
		"addq	%%rdi,%%rbx								\n\t	leaq	(%%r9,%%rbx),%%r11	\n\t"\
		"addq	%%rdi,%%rcx								\n\t	leaq	(%%r9,%%rcx),%%r12	\n\t"\
		"addq	%%rdi,%%rdx								\n\t	leaq	(%%r9,%%rdx),%%r13	\n\t"\
		/* SSE2_RADIX4_DIF_4TWIDDLE_B(r10,c02) */				/* SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO_B(r18,c06) */\
		"addq	$0x200,%%rsi	/* r10 */				\n\t"\
		"vmovaps	     (%%rax),%%ymm0					\n\t	vmovaps	     (%%r10),%%ymm8 	\n\t"\
		"vmovaps	     (%%rcx),%%ymm4					\n\t	vmovaps	     (%%r12),%%ymm12	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1					\n\t	vmovaps	0x020(%%r10),%%ymm9 	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5					\n\t	vmovaps	0x020(%%r12),%%ymm13	\n\t"\
		"vmovaps	0x8e0(%%rsi),%%ymm6	/* c02 */		\n\t	vmovaps	0x9e0(%%rsi),%%ymm14	/* c06 */	\n\t"\
		"vmovaps	0x900(%%rsi),%%ymm7					\n\t	vmovaps	0xa00(%%rsi),%%ymm15	\n\t"\
		"vmovaps	%%ymm0,%%ymm2						\n\t	vmovaps	%%ymm8 ,%%ymm10	\n\t"\
		"vmovaps	%%ymm1,%%ymm3						\n\t	vmovaps	%%ymm9 ,%%ymm11	\n\t"\
		"vmulpd	%%ymm6,%%ymm0,%%ymm0					\n\t	vmulpd	%%ymm14,%%ymm8 ,%%ymm8 	\n\t"\
		"vmulpd	%%ymm6,%%ymm1,%%ymm1					\n\t	vmulpd	%%ymm14,%%ymm9 ,%%ymm9 	\n\t"\
	"vfnmadd231pd	%%ymm7,%%ymm3,%%ymm0			\n\t	vfnmadd231pd	%%ymm15,%%ymm11,%%ymm8 	\n\t"\
	" vfmadd231pd	%%ymm7,%%ymm2,%%ymm1			\n\t	 vfmadd231pd	%%ymm15,%%ymm10,%%ymm9 	\n\t"\
		"vmovaps	%%ymm4,%%ymm6						\n\t	vmovaps	%%ymm12,%%ymm14	\n\t"\
		"vmovaps	%%ymm5,%%ymm7						\n\t	vmovaps	%%ymm13,%%ymm15	\n\t"\
		"vmulpd	0x920(%%rsi),%%ymm4,%%ymm4	/* c12 */	\n\t	vmulpd	0xa20(%%rsi),%%ymm12,%%ymm12	/* c16 */	\n\t"\
		"vmulpd	0x920(%%rsi),%%ymm5,%%ymm5				\n\t	vmulpd	0xa20(%%rsi),%%ymm13,%%ymm13	\n\t"\
	"vfnmadd231pd	0x940(%%rsi),%%ymm7,%%ymm4		\n\t	vfnmadd231pd	0xa40(%%rsi),%%ymm15,%%ymm12	\n\t"\
	" vfmadd231pd	0x940(%%rsi),%%ymm6,%%ymm5		\n\t	 vfmadd231pd	0xa40(%%rsi),%%ymm14,%%ymm13	\n\t"\
		"vmovaps	%%ymm0,%%ymm2						\n\t	vmovaps	%%ymm8 ,%%ymm10	\n\t"\
		"vmovaps	%%ymm1,%%ymm3						\n\t	vmovaps	%%ymm9 ,%%ymm11	\n\t"\
	"vmovaps	-0x20(%%r8),%%ymm14	\n\t"/* ymm14 free; use for 1.0 */\
		"vfmadd132pd	%%ymm14,%%ymm4,%%ymm0			\n\t	vfmadd132pd	%%ymm14,%%ymm12,%%ymm8 	\n\t"\
		"vfmadd132pd	%%ymm14,%%ymm5,%%ymm1			\n\t	vfmadd132pd	%%ymm14,%%ymm13,%%ymm9 	\n\t"\
		"vfmsub132pd	%%ymm14,%%ymm4,%%ymm2			\n\t	vfmsub132pd	%%ymm14,%%ymm12,%%ymm10	\n\t"\
		"vfmsub132pd	%%ymm14,%%ymm5,%%ymm3			\n\t	vfmsub132pd	%%ymm14,%%ymm13,%%ymm11	\n\t"\
		"vmovaps	     (%%rdx),%%ymm4					\n\t	vmovaps	     (%%r13),%%ymm12	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm5					\n\t	vmovaps	0x020(%%r13),%%ymm13	\n\t"\
		"vmovaps	     (%%rdx),%%ymm6					\n\t	vmovaps	     (%%r13),%%ymm14	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm7					\n\t	vmovaps	0x020(%%r13),%%ymm15	\n\t"\
		"vmulpd	0x9a0(%%rsi),%%ymm4,%%ymm4	/* c1A */	\n\t	vmulpd	0xaa0(%%rsi),%%ymm12,%%ymm12	/* c1E */	\n\t"\
		"vmulpd	0x9a0(%%rsi),%%ymm5,%%ymm5				\n\t	vmulpd	0xaa0(%%rsi),%%ymm13,%%ymm13	\n\t"\
	"vfnmadd231pd	0x9c0(%%rsi),%%ymm7,%%ymm4		\n\t	vfnmadd231pd	0xac0(%%rsi),%%ymm15,%%ymm12	\n\t"\
	" vfmadd231pd	0x9c0(%%rsi),%%ymm6,%%ymm5		\n\t	 vfmadd231pd	0xac0(%%rsi),%%ymm14,%%ymm13	\n\t"\
		"vmovaps	%%ymm5,0x020(%%rsi)					\n\t	vmovaps	%%ymm13,0x120(%%rsi)	\n\t"\
		"vmovaps	%%ymm4,     (%%rsi)					\n\t	vmovaps	%%ymm12,0x100(%%rsi)	/* tmpstr r08 */\n\t"\
		"vmovaps	     (%%rbx),%%ymm4					\n\t	vmovaps	     (%%r11),%%ymm12	\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm5					\n\t	vmovaps	0x020(%%r11),%%ymm13	\n\t"\
		"vmovaps	     (%%rbx),%%ymm6					\n\t	vmovaps	     (%%r11),%%ymm14	\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm7					\n\t	vmovaps	0x020(%%r11),%%ymm15	\n\t"\
		"vmulpd	0x960(%%rsi),%%ymm4,%%ymm4	/* c0A */	\n\t	vmulpd	0xa60(%%rsi),%%ymm12,%%ymm12	/* c0E */	\n\t"\
		"vmulpd	0x960(%%rsi),%%ymm5,%%ymm5				\n\t	vmulpd	0xa60(%%rsi),%%ymm13,%%ymm13	\n\t"\
	"vfnmadd231pd	0x980(%%rsi),%%ymm7,%%ymm4		\n\t	vfnmadd231pd	0xa80(%%rsi),%%ymm15,%%ymm12	\n\t"\
	" vfmadd231pd	0x980(%%rsi),%%ymm6,%%ymm5		\n\t	 vfmadd231pd	0xa80(%%rsi),%%ymm14,%%ymm13	\n\t"\
		"vmovaps	%%ymm4,%%ymm6						\n\t	vmovaps	%%ymm12,%%ymm14	\n\t"\
		"vmovaps	%%ymm5,%%ymm7						\n\t	vmovaps	%%ymm13,%%ymm15	\n\t"\
		"vaddpd	     (%%rsi),%%ymm6,%%ymm6				\n\t	vaddpd	0x100(%%rsi),%%ymm14,%%ymm14	\n\t"\
		"vaddpd	0x020(%%rsi),%%ymm7,%%ymm7				\n\t	vaddpd	0x120(%%rsi),%%ymm15,%%ymm15	\n\t"\
		"vsubpd	     (%%rsi),%%ymm4,%%ymm4				\n\t	vsubpd	0x100(%%rsi),%%ymm12,%%ymm12	\n\t"\
		"vsubpd	0x020(%%rsi),%%ymm5,%%ymm5				\n\t	vsubpd	0x120(%%rsi),%%ymm13,%%ymm13	\n\t"\
	"vmovaps	%%ymm11,(%%rsi) 	\n\t"/* spill ymm11 to make room for 1.0 */"	vmovaps	-0x20(%%r8),%%ymm11 	\n\t"/* one */\
		"vfmsub132pd	%%ymm11,%%ymm6,%%ymm0			\n\t	vsubpd	     %%ymm14,%%ymm8 ,%%ymm8 	\n\t"\
		"vfmsub132pd	%%ymm11,%%ymm5,%%ymm2			\n\t	vsubpd	     %%ymm15,%%ymm9 ,%%ymm9 	\n\t	vmovaps	%%ymm8 ,0x180(%%rsi)	\n\t"\
		"vfmsub132pd	%%ymm11,%%ymm7,%%ymm1			\n\t	vfmsub132pd	%%ymm11,%%ymm13,%%ymm10	\n\t"\
		"vfmsub132pd	%%ymm11,%%ymm4,%%ymm3			\n\t	vfmsub132pd	(%%rsi),%%ymm12,%%ymm11	\n\t"\
	"vmovaps	%%ymm12,(%%rsi) 	\n\t"/* spill ymm12 to make room for 2.0 */"	vmovaps	(%%r8),%%ymm12 	\n\t"/* two */\
	"vfmadd132pd	%%ymm12,%%ymm0,%%ymm6			\n\t	vfmadd132pd	%%ymm12,%%ymm8 ,%%ymm14	\n\t"\
	"vfmadd132pd	%%ymm12,%%ymm2,%%ymm5			\n\t	vfmadd132pd	%%ymm12,%%ymm10,%%ymm13	\n\t"\
	"vfmadd132pd	%%ymm12,%%ymm1,%%ymm7			\n\t	vfmadd132pd	%%ymm12,%%ymm9 ,%%ymm15	\n\t"\
	"vfmadd132pd	%%ymm12,%%ymm3,%%ymm4			\n\t	vfmadd132pd	(%%rsi),%%ymm11,%%ymm12	\n\t"\
		/* SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r10) */\
		"												\n\t	vmovaps	%%ymm14,0x100(%%rsi)	\n\t"\
		"												\n\t	vmovaps	%%ymm15,0x120(%%rsi)	\n\t"\
		"												\n\t	vmovaps	%%ymm10,%%ymm14	\n\t"\
		"												\n\t	vmovaps	%%ymm13,%%ymm15	\n\t"\
		"												\n\t	vmovaps	-0x20(%%r8 ),%%ymm8 	\n\t"/* one */\
		"														vfmsub132pd	%%ymm8 ,%%ymm12,%%ymm10	\n\t"\
		"														vfmsub132pd	%%ymm8 ,%%ymm11,%%ymm13	\n\t"\
		"														vfmadd132pd	%%ymm8 ,%%ymm12,%%ymm14	\n\t"\
		"														vfmadd132pd	%%ymm8 ,%%ymm11,%%ymm15	\n\t"\
		"												\n\t	vmovaps	0x600(%%rsi),%%ymm8 	\n\t"/* isrt2 */\
		"vmovaps	0x100(%%rsi),%%ymm11			\n\t	vfnmadd231pd %%ymm8,%%ymm10,%%ymm2	\n\t"\
		"vmovaps	0x120(%%rsi),%%ymm12			\n\t	vfnmadd231pd %%ymm8,%%ymm13,%%ymm3	\n\t"\
		"vfmsub132pd	-0x20(%%r8),%%ymm11,%%ymm6		\n\t	vfnmadd231pd %%ymm8,%%ymm14,%%ymm4	\n\t"\
		"vfmsub132pd	-0x20(%%r8),%%ymm9 ,%%ymm0		\n\t	vfnmadd231pd %%ymm8,%%ymm15,%%ymm5	\n\t"\
		"vsubpd	0x180(%%rsi),%%ymm1 ,%%ymm1 			\n\t	vfmsub132pd	-0x20(%%r8),%%ymm12 ,%%ymm7	\n\t"\
	"vmovaps	(%%r8),%%ymm8	\n\t"/* rcol: *2 => *sqrt2 due to above FMA leaving ymm10,13-15 unmultiplied-by-isrt2: */\
		"vmovaps	%%ymm6 ,0x100(%%rsi)				\n\t	vmovaps	%%ymm2 ,0x140(%%rsi)	\n\t"\
		"vmovaps	%%ymm0 ,0x080(%%rsi)				\n\t	vmovaps	%%ymm5 ,0x0c0(%%rsi)	\n\t"\
		"vmovaps	%%ymm7 ,0x120(%%rsi)				\n\t	vmovaps	%%ymm4 ,0x160(%%rsi)	\n\t"\
		"vmovaps	%%ymm1 ,0x1a0(%%rsi)				\n\t	vmovaps	%%ymm3 ,0x1e0(%%rsi)	\n\t"\
	"vfmadd132pd	%%ymm8 ,%%ymm6 ,%%ymm11			\n\t	vfmadd132pd	0x20(%%r8),%%ymm2 ,%%ymm10	\n\t"\
	"vfmadd132pd	%%ymm8 ,%%ymm0 ,%%ymm9 			\n\t	vfmadd132pd	0x20(%%r8),%%ymm5 ,%%ymm15	\n\t"\
	"vfmadd132pd	%%ymm8 ,%%ymm7 ,%%ymm12			\n\t	vfmadd132pd	0x20(%%r8),%%ymm4 ,%%ymm14	\n\t"\
	"vfmadd132pd	0x180(%%rsi),%%ymm1 ,%%ymm8 	\n\t	vfmadd132pd	0x20(%%r8),%%ymm3 ,%%ymm13	\n\t"\
		"vmovaps	%%ymm11,     (%%rsi)				\n\t	vmovaps	%%ymm10,0x040(%%rsi)	\n\t"\
		"vmovaps	%%ymm9 ,0x180(%%rsi)				\n\t	vmovaps	%%ymm15,0x1c0(%%rsi)	\n\t"\
		"vmovaps	%%ymm12,0x020(%%rsi)				\n\t	vmovaps	%%ymm14,0x060(%%rsi)	\n\t"\
		"vmovaps	%%ymm8 ,0x0a0(%%rsi)				\n\t	vmovaps	%%ymm13,0x0e0(%%rsi)	\n\t"\
		"\n\t"\
	/***************************************/\
		"\n\t"\
	/*...Block 3: */\
		"movq	%[__add0],%%rax	\n\t"\
		"movslq	%[__p01],%%rdi	\n\t"/* Do this way [rather than repeatedly +=p1] since array-padding scheme means p2 == p1+p1 not guaranteed. */\
		"movslq	%[__p08],%%rbx	\n\t"/* Can't get these via simple load-one-and-shift-as-needed due to array padding scheme */\
		"movslq	%[__p10],%%rcx							\n\t"\
		"movslq	%[__p18],%%rdx							\n\t"	/* p04<<3 still in r9 */\
		"leaq	(%%rax,%%rdi,8),%%rax	/* &a[j1+p1] */	\n\t	leaq	(%%r9,%%rax),%%r10	\n\t"\
		"leaq	(%%rax,%%rbx,8),%%rbx					\n\t	leaq	(%%r9,%%rbx),%%r11	\n\t"\
		"leaq	(%%rax,%%rcx,8),%%rcx					\n\t	leaq	(%%r9,%%rcx),%%r12	\n\t"\
		"leaq	(%%rax,%%rdx,8),%%rdx					\n\t	leaq	(%%r9,%%rdx),%%r13	\n\t"\
		/* SSE2_RADIX4_DIF_4TWIDDLE_B(r20,c01) */				/* SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO_B(r28,c05) */\
		"addq	$0x200,%%rsi	/* r20 */				\n\t"\
		"vmovaps	     (%%rax),%%ymm0					\n\t	vmovaps	     (%%r10),%%ymm8 	\n\t"\
		"vmovaps	     (%%rcx),%%ymm4					\n\t	vmovaps	     (%%r12),%%ymm12	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1					\n\t	vmovaps	0x020(%%r10),%%ymm9 	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5					\n\t	vmovaps	0x020(%%r12),%%ymm13	\n\t"\
		"vmovaps	0x8e0(%%rsi),%%ymm6	/* c01 */		\n\t	vmovaps	0x9e0(%%rsi),%%ymm14	/* c05 */	\n\t"\
		"vmovaps	0x900(%%rsi),%%ymm7					\n\t	vmovaps	0xa00(%%rsi),%%ymm15	\n\t"\
		"vmovaps	%%ymm0,%%ymm2						\n\t	vmovaps	%%ymm8 ,%%ymm10	\n\t"\
		"vmovaps	%%ymm1,%%ymm3						\n\t	vmovaps	%%ymm9 ,%%ymm11	\n\t"\
		"vmulpd	%%ymm6,%%ymm0,%%ymm0					\n\t	vmulpd	%%ymm14,%%ymm8 ,%%ymm8 	\n\t"\
		"vmulpd	%%ymm6,%%ymm1,%%ymm1					\n\t	vmulpd	%%ymm14,%%ymm9 ,%%ymm9 	\n\t"\
	"vfnmadd231pd	%%ymm7,%%ymm3,%%ymm0			\n\t	vfnmadd231pd	%%ymm15,%%ymm11,%%ymm8 	\n\t"\
	" vfmadd231pd	%%ymm7,%%ymm2,%%ymm1			\n\t	 vfmadd231pd	%%ymm15,%%ymm10,%%ymm9 	\n\t"\
		"vmovaps	%%ymm4,%%ymm6						\n\t	vmovaps	%%ymm12,%%ymm14	\n\t"\
		"vmovaps	%%ymm5,%%ymm7						\n\t	vmovaps	%%ymm13,%%ymm15	\n\t"\
		"vmulpd	0x920(%%rsi),%%ymm4,%%ymm4	/* c11 */	\n\t	vmulpd	0xa20(%%rsi),%%ymm12,%%ymm12	/* c15 */	\n\t"\
		"vmulpd	0x920(%%rsi),%%ymm5,%%ymm5				\n\t	vmulpd	0xa20(%%rsi),%%ymm13,%%ymm13	\n\t"\
	"vfnmadd231pd	0x940(%%rsi),%%ymm7,%%ymm4		\n\t	vfnmadd231pd	0xa40(%%rsi),%%ymm15,%%ymm12	\n\t"\
	" vfmadd231pd	0x940(%%rsi),%%ymm6,%%ymm5		\n\t	 vfmadd231pd	0xa40(%%rsi),%%ymm14,%%ymm13	\n\t"\
		"vmovaps	%%ymm0,%%ymm2						\n\t	vmovaps	%%ymm8 ,%%ymm10	\n\t"\
		"vmovaps	%%ymm1,%%ymm3						\n\t	vmovaps	%%ymm9 ,%%ymm11	\n\t"\
	"vmovaps	-0x20(%%r8),%%ymm14	\n\t"/* ymm14 free; use for 1.0 */\
		"vfmadd132pd	%%ymm14,%%ymm4,%%ymm0			\n\t	vfmadd132pd	%%ymm14,%%ymm12,%%ymm8 	\n\t"\
		"vfmadd132pd	%%ymm14,%%ymm5,%%ymm1			\n\t	vfmadd132pd	%%ymm14,%%ymm13,%%ymm9 	\n\t"\
		"vfmsub132pd	%%ymm14,%%ymm4,%%ymm2			\n\t	vfmsub132pd	%%ymm14,%%ymm12,%%ymm10	\n\t"\
		"vfmsub132pd	%%ymm14,%%ymm5,%%ymm3			\n\t	vfmsub132pd	%%ymm14,%%ymm13,%%ymm11	\n\t"\
		"vmovaps	     (%%rdx),%%ymm4					\n\t	vmovaps	     (%%r13),%%ymm12	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm5					\n\t	vmovaps	0x020(%%r13),%%ymm13	\n\t"\
		"vmovaps	     (%%rdx),%%ymm6					\n\t	vmovaps	     (%%r13),%%ymm14	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm7					\n\t	vmovaps	0x020(%%r13),%%ymm15	\n\t"\
		"vmulpd	0x9a0(%%rsi),%%ymm4,%%ymm4	/* c19 */	\n\t	vmulpd	0xaa0(%%rsi),%%ymm12,%%ymm12	/* c1D */	\n\t"\
		"vmulpd	0x9a0(%%rsi),%%ymm5,%%ymm5				\n\t	vmulpd	0xaa0(%%rsi),%%ymm13,%%ymm13	\n\t"\
	"vfnmadd231pd	0x9c0(%%rsi),%%ymm7,%%ymm4		\n\t	vfnmadd231pd	0xac0(%%rsi),%%ymm15,%%ymm12	\n\t"\
	" vfmadd231pd	0x9c0(%%rsi),%%ymm6,%%ymm5		\n\t	 vfmadd231pd	0xac0(%%rsi),%%ymm14,%%ymm13	\n\t"\
		"vmovaps	%%ymm5,0x020(%%rsi)					\n\t	vmovaps	%%ymm13,0x120(%%rsi)	\n\t"\
		"vmovaps	%%ymm4,     (%%rsi)					\n\t	vmovaps	%%ymm12,0x100(%%rsi)	/* tmpstr r08 */\n\t"\
		"vmovaps	     (%%rbx),%%ymm4					\n\t	vmovaps	     (%%r11),%%ymm12	\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm5					\n\t	vmovaps	0x020(%%r11),%%ymm13	\n\t"\
		"vmovaps	     (%%rbx),%%ymm6					\n\t	vmovaps	     (%%r11),%%ymm14	\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm7					\n\t	vmovaps	0x020(%%r11),%%ymm15	\n\t"\
		"vmulpd	0x960(%%rsi),%%ymm4,%%ymm4	/* c09 */	\n\t	vmulpd	0xa60(%%rsi),%%ymm12,%%ymm12	/* c0D */	\n\t"\
		"vmulpd	0x960(%%rsi),%%ymm5,%%ymm5				\n\t	vmulpd	0xa60(%%rsi),%%ymm13,%%ymm13	\n\t"\
	"vfnmadd231pd	0x980(%%rsi),%%ymm7,%%ymm4		\n\t	vfnmadd231pd	0xa80(%%rsi),%%ymm15,%%ymm12	\n\t"\
	" vfmadd231pd	0x980(%%rsi),%%ymm6,%%ymm5		\n\t	 vfmadd231pd	0xa80(%%rsi),%%ymm14,%%ymm13	\n\t"\
		"vmovaps	%%ymm4,%%ymm6						\n\t	vmovaps	%%ymm12,%%ymm14	\n\t"\
		"vmovaps	%%ymm5,%%ymm7						\n\t	vmovaps	%%ymm13,%%ymm15	\n\t"\
		"vaddpd	     (%%rsi),%%ymm6,%%ymm6				\n\t	vaddpd	0x100(%%rsi),%%ymm14,%%ymm14	\n\t"\
		"vaddpd	0x020(%%rsi),%%ymm7,%%ymm7				\n\t	vaddpd	0x120(%%rsi),%%ymm15,%%ymm15	\n\t"\
		"vsubpd	     (%%rsi),%%ymm4,%%ymm4				\n\t	vsubpd	0x100(%%rsi),%%ymm12,%%ymm12	\n\t"\
		"vsubpd	0x020(%%rsi),%%ymm5,%%ymm5				\n\t	vsubpd	0x120(%%rsi),%%ymm13,%%ymm13	\n\t"\
	"vmovaps	%%ymm11,(%%rsi) 	\n\t"/* spill ymm11 to make room for 1.0 */"	vmovaps	-0x20(%%r8),%%ymm11 	\n\t"/* one */\
		"vfmsub132pd	%%ymm11,%%ymm6,%%ymm0			\n\t	vsubpd	     %%ymm14,%%ymm8 ,%%ymm8 	\n\t"\
		"vfmsub132pd	%%ymm11,%%ymm5,%%ymm2			\n\t	vsubpd	     %%ymm15,%%ymm9 ,%%ymm9 	\n\t	vmovaps	%%ymm8 ,0x180(%%rsi)	\n\t"\
		"vfmsub132pd	%%ymm11,%%ymm7,%%ymm1			\n\t	vfmsub132pd	%%ymm11,%%ymm13,%%ymm10	\n\t"\
		"vfmsub132pd	%%ymm11,%%ymm4,%%ymm3			\n\t	vfmsub132pd	(%%rsi),%%ymm12,%%ymm11	\n\t"\
	"vmovaps	%%ymm12,(%%rsi) 	\n\t"/* spill ymm12 to make room for 2.0 */"	vmovaps	(%%r8),%%ymm12 	\n\t"/* two */\
	"vfmadd132pd	%%ymm12,%%ymm0,%%ymm6			\n\t	vfmadd132pd	%%ymm12,%%ymm8 ,%%ymm14	\n\t"\
	"vfmadd132pd	%%ymm12,%%ymm2,%%ymm5			\n\t	vfmadd132pd	%%ymm12,%%ymm10,%%ymm13	\n\t"\
	"vfmadd132pd	%%ymm12,%%ymm1,%%ymm7			\n\t	vfmadd132pd	%%ymm12,%%ymm9 ,%%ymm15	\n\t"\
	"vfmadd132pd	%%ymm12,%%ymm3,%%ymm4			\n\t	vfmadd132pd	(%%rsi),%%ymm11,%%ymm12	\n\t"\
		/* SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r20) */\
		"												\n\t	vmovaps	%%ymm14,0x100(%%rsi)	\n\t"\
		"												\n\t	vmovaps	%%ymm15,0x120(%%rsi)	\n\t"\
		"												\n\t	vmovaps	%%ymm10,%%ymm14	\n\t"\
		"												\n\t	vmovaps	%%ymm13,%%ymm15	\n\t"\
		"												\n\t	vmovaps	-0x20(%%r8 ),%%ymm8 	\n\t"/* one */\
		"														vfmsub132pd	%%ymm8 ,%%ymm12,%%ymm10	\n\t"\
		"														vfmsub132pd	%%ymm8 ,%%ymm11,%%ymm13	\n\t"\
		"														vfmadd132pd	%%ymm8 ,%%ymm12,%%ymm14	\n\t"\
		"														vfmadd132pd	%%ymm8 ,%%ymm11,%%ymm15	\n\t"\
		"												\n\t	vmovaps	0x400(%%rsi),%%ymm8 	\n\t"/* isrt2 */\
		"vmovaps	0x100(%%rsi),%%ymm11			\n\t	vfnmadd231pd %%ymm8,%%ymm10,%%ymm2	\n\t"\
		"vmovaps	0x120(%%rsi),%%ymm12			\n\t	vfnmadd231pd %%ymm8,%%ymm13,%%ymm3	\n\t"\
		"vfmsub132pd	-0x20(%%r8),%%ymm11,%%ymm6		\n\t	vfnmadd231pd %%ymm8,%%ymm14,%%ymm4	\n\t"\
		"vfmsub132pd	-0x20(%%r8),%%ymm9 ,%%ymm0		\n\t	vfnmadd231pd %%ymm8,%%ymm15,%%ymm5	\n\t"\
		"vsubpd	0x180(%%rsi),%%ymm1 ,%%ymm1 			\n\t	vfmsub132pd	-0x20(%%r8),%%ymm12 ,%%ymm7	\n\t"\
	"vmovaps	(%%r8),%%ymm8	\n\t"/* rcol: *2 => *sqrt2 due to above FMA leaving ymm10,13-15 unmultiplied-by-isrt2: */\
		"vmovaps	%%ymm6 ,0x100(%%rsi)				\n\t	vmovaps	%%ymm2 ,0x140(%%rsi)	\n\t"\
		"vmovaps	%%ymm0 ,0x080(%%rsi)				\n\t	vmovaps	%%ymm5 ,0x0c0(%%rsi)	\n\t"\
		"vmovaps	%%ymm7 ,0x120(%%rsi)				\n\t	vmovaps	%%ymm4 ,0x160(%%rsi)	\n\t"\
		"vmovaps	%%ymm1 ,0x1a0(%%rsi)				\n\t	vmovaps	%%ymm3 ,0x1e0(%%rsi)	\n\t"\
	"vfmadd132pd	%%ymm8 ,%%ymm6 ,%%ymm11			\n\t	vfmadd132pd	0x20(%%r8),%%ymm2 ,%%ymm10	\n\t"\
	"vfmadd132pd	%%ymm8 ,%%ymm0 ,%%ymm9 			\n\t	vfmadd132pd	0x20(%%r8),%%ymm5 ,%%ymm15	\n\t"\
	"vfmadd132pd	%%ymm8 ,%%ymm7 ,%%ymm12			\n\t	vfmadd132pd	0x20(%%r8),%%ymm4 ,%%ymm14	\n\t"\
	"vfmadd132pd	0x180(%%rsi),%%ymm1 ,%%ymm8 	\n\t	vfmadd132pd	0x20(%%r8),%%ymm3 ,%%ymm13	\n\t"\
		"vmovaps	%%ymm11,     (%%rsi)				\n\t	vmovaps	%%ymm10,0x040(%%rsi)	\n\t"\
		"vmovaps	%%ymm9 ,0x180(%%rsi)				\n\t	vmovaps	%%ymm15,0x1c0(%%rsi)	\n\t"\
		"vmovaps	%%ymm12,0x020(%%rsi)				\n\t	vmovaps	%%ymm14,0x060(%%rsi)	\n\t"\
		"vmovaps	%%ymm8 ,0x0a0(%%rsi)				\n\t	vmovaps	%%ymm13,0x0e0(%%rsi)	\n\t"\
		"\n\t"\
	/***************************************/\
		"\n\t"\
	/*...Block 4: */\
		"movq	%[__add0],%%rax	\n\t"\
		"movslq	%[__p03],%%rdi	\n\t"\
		"movslq	%[__p08],%%rbx	\n\t"\
		"movslq	%[__p10],%%rcx							\n\t"\
		"movslq	%[__p18],%%rdx							\n\t"	/* p04<<3 still in r9 */\
		"leaq	(%%rax,%%rdi,8),%%rax	/* &a[j1+p3] */	\n\t	leaq	(%%r9,%%rax),%%r10	\n\t"\
		"leaq	(%%rax,%%rbx,8),%%rbx					\n\t	leaq	(%%r9,%%rbx),%%r11	\n\t"\
		"leaq	(%%rax,%%rcx,8),%%rcx					\n\t	leaq	(%%r9,%%rcx),%%r12	\n\t"\
		"leaq	(%%rax,%%rdx,8),%%rdx					\n\t	leaq	(%%r9,%%rdx),%%r13	\n\t"\
		/* SSE2_RADIX4_DIF_4TWIDDLE_B(r30,c03) */				/* SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO_B(r38,c07) */\
		"addq	$0x200,%%rsi	/* r30 */				\n\t"\
		"vmovaps	     (%%rax),%%ymm0					\n\t	vmovaps	     (%%r10),%%ymm8 	\n\t"\
		"vmovaps	     (%%rcx),%%ymm4					\n\t	vmovaps	     (%%r12),%%ymm12	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1					\n\t	vmovaps	0x020(%%r10),%%ymm9 	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5					\n\t	vmovaps	0x020(%%r12),%%ymm13	\n\t"\
		"vmovaps	0x8e0(%%rsi),%%ymm6	/* c03 */		\n\t	vmovaps	0x9e0(%%rsi),%%ymm14	/* c07 */	\n\t"\
		"vmovaps	0x900(%%rsi),%%ymm7					\n\t	vmovaps	0xa00(%%rsi),%%ymm15	\n\t"\
		"vmovaps	%%ymm0,%%ymm2						\n\t	vmovaps	%%ymm8 ,%%ymm10	\n\t"\
		"vmovaps	%%ymm1,%%ymm3						\n\t	vmovaps	%%ymm9 ,%%ymm11	\n\t"\
		"vmulpd	%%ymm6,%%ymm0,%%ymm0					\n\t	vmulpd	%%ymm14,%%ymm8 ,%%ymm8 	\n\t"\
		"vmulpd	%%ymm6,%%ymm1,%%ymm1					\n\t	vmulpd	%%ymm14,%%ymm9 ,%%ymm9 	\n\t"\
	"vfnmadd231pd	%%ymm7,%%ymm3,%%ymm0			\n\t	vfnmadd231pd	%%ymm15,%%ymm11,%%ymm8 	\n\t"\
	" vfmadd231pd	%%ymm7,%%ymm2,%%ymm1			\n\t	 vfmadd231pd	%%ymm15,%%ymm10,%%ymm9 	\n\t"\
		"vmovaps	%%ymm4,%%ymm6						\n\t	vmovaps	%%ymm12,%%ymm14	\n\t"\
		"vmovaps	%%ymm5,%%ymm7						\n\t	vmovaps	%%ymm13,%%ymm15	\n\t"\
		"vmulpd	0x920(%%rsi),%%ymm4,%%ymm4	/* c13 */	\n\t	vmulpd	0xa20(%%rsi),%%ymm12,%%ymm12	/* c17 */	\n\t"\
		"vmulpd	0x920(%%rsi),%%ymm5,%%ymm5				\n\t	vmulpd	0xa20(%%rsi),%%ymm13,%%ymm13	\n\t"\
	"vfnmadd231pd	0x940(%%rsi),%%ymm7,%%ymm4		\n\t	vfnmadd231pd	0xa40(%%rsi),%%ymm15,%%ymm12	\n\t"\
	" vfmadd231pd	0x940(%%rsi),%%ymm6,%%ymm5		\n\t	 vfmadd231pd	0xa40(%%rsi),%%ymm14,%%ymm13	\n\t"\
		"vmovaps	%%ymm0,%%ymm2						\n\t	vmovaps	%%ymm8 ,%%ymm10	\n\t"\
		"vmovaps	%%ymm1,%%ymm3						\n\t	vmovaps	%%ymm9 ,%%ymm11	\n\t"\
	"vmovaps	-0x20(%%r8),%%ymm14	\n\t"/* ymm14 free; use for 1.0 */\
		"vfmadd132pd	%%ymm14,%%ymm4,%%ymm0			\n\t	vfmadd132pd	%%ymm14,%%ymm12,%%ymm8 	\n\t"\
		"vfmadd132pd	%%ymm14,%%ymm5,%%ymm1			\n\t	vfmadd132pd	%%ymm14,%%ymm13,%%ymm9 	\n\t"\
		"vfmsub132pd	%%ymm14,%%ymm4,%%ymm2			\n\t	vfmsub132pd	%%ymm14,%%ymm12,%%ymm10	\n\t"\
		"vfmsub132pd	%%ymm14,%%ymm5,%%ymm3			\n\t	vfmsub132pd	%%ymm14,%%ymm13,%%ymm11	\n\t"\
		"vmovaps	     (%%rdx),%%ymm4					\n\t	vmovaps	     (%%r13),%%ymm12	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm5					\n\t	vmovaps	0x020(%%r13),%%ymm13	\n\t"\
		"vmovaps	     (%%rdx),%%ymm6					\n\t	vmovaps	     (%%r13),%%ymm14	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm7					\n\t	vmovaps	0x020(%%r13),%%ymm15	\n\t"\
		"vmulpd	0x9a0(%%rsi),%%ymm4,%%ymm4	/* c1B */	\n\t	vmulpd	0xaa0(%%rsi),%%ymm12,%%ymm12	/* c1F */	\n\t"\
		"vmulpd	0x9a0(%%rsi),%%ymm5,%%ymm5				\n\t	vmulpd	0xaa0(%%rsi),%%ymm13,%%ymm13	\n\t"\
	"vfnmadd231pd	0x9c0(%%rsi),%%ymm7,%%ymm4		\n\t	vfnmadd231pd	0xac0(%%rsi),%%ymm15,%%ymm12	\n\t"\
	" vfmadd231pd	0x9c0(%%rsi),%%ymm6,%%ymm5		\n\t	 vfmadd231pd	0xac0(%%rsi),%%ymm14,%%ymm13	\n\t"\
		"vmovaps	%%ymm5,0x020(%%rsi)					\n\t	vmovaps	%%ymm13,0x120(%%rsi)	\n\t"\
		"vmovaps	%%ymm4,     (%%rsi)					\n\t	vmovaps	%%ymm12,0x100(%%rsi)	/* tmpstr r08 */\n\t"\
		"vmovaps	     (%%rbx),%%ymm4					\n\t	vmovaps	     (%%r11),%%ymm12	\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm5					\n\t	vmovaps	0x020(%%r11),%%ymm13	\n\t"\
		"vmovaps	     (%%rbx),%%ymm6					\n\t	vmovaps	     (%%r11),%%ymm14	\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm7					\n\t	vmovaps	0x020(%%r11),%%ymm15	\n\t"\
		"vmulpd	0x960(%%rsi),%%ymm4,%%ymm4	/* c0B */	\n\t	vmulpd	0xa60(%%rsi),%%ymm12,%%ymm12	/* c0F */	\n\t"\
		"vmulpd	0x960(%%rsi),%%ymm5,%%ymm5				\n\t	vmulpd	0xa60(%%rsi),%%ymm13,%%ymm13	\n\t"\
	"vfnmadd231pd	0x980(%%rsi),%%ymm7,%%ymm4		\n\t	vfnmadd231pd	0xa80(%%rsi),%%ymm15,%%ymm12	\n\t"\
	" vfmadd231pd	0x980(%%rsi),%%ymm6,%%ymm5		\n\t	 vfmadd231pd	0xa80(%%rsi),%%ymm14,%%ymm13	\n\t"\
		"vmovaps	%%ymm4,%%ymm6						\n\t	vmovaps	%%ymm12,%%ymm14	\n\t"\
		"vmovaps	%%ymm5,%%ymm7						\n\t	vmovaps	%%ymm13,%%ymm15	\n\t"\
		"vaddpd	     (%%rsi),%%ymm6,%%ymm6				\n\t	vaddpd	0x100(%%rsi),%%ymm14,%%ymm14	\n\t"\
		"vaddpd	0x020(%%rsi),%%ymm7,%%ymm7				\n\t	vaddpd	0x120(%%rsi),%%ymm15,%%ymm15	\n\t"\
		"vsubpd	     (%%rsi),%%ymm4,%%ymm4				\n\t	vsubpd	0x100(%%rsi),%%ymm12,%%ymm12	\n\t"\
		"vsubpd	0x020(%%rsi),%%ymm5,%%ymm5				\n\t	vsubpd	0x120(%%rsi),%%ymm13,%%ymm13	\n\t"\
	"vmovaps	%%ymm11,(%%rsi) 	\n\t"/* spill ymm11 to make room for 1.0 */"	vmovaps	-0x20(%%r8),%%ymm11 	\n\t"/* one */\
		"vfmsub132pd	%%ymm11,%%ymm6,%%ymm0			\n\t	vsubpd	     %%ymm14,%%ymm8 ,%%ymm8 	\n\t"\
		"vfmsub132pd	%%ymm11,%%ymm5,%%ymm2			\n\t	vsubpd	     %%ymm15,%%ymm9 ,%%ymm9 	\n\t	vmovaps	%%ymm8 ,0x180(%%rsi)	\n\t"\
		"vfmsub132pd	%%ymm11,%%ymm7,%%ymm1			\n\t	vfmsub132pd	%%ymm11,%%ymm13,%%ymm10	\n\t"\
		"vfmsub132pd	%%ymm11,%%ymm4,%%ymm3			\n\t	vfmsub132pd	(%%rsi),%%ymm12,%%ymm11	\n\t"\
	"vmovaps	%%ymm12,(%%rsi) 	\n\t"/* spill ymm12 to make room for 2.0 */"	vmovaps	(%%r8),%%ymm12 	\n\t"/* two */\
	"vfmadd132pd	%%ymm12,%%ymm0,%%ymm6			\n\t	vfmadd132pd	%%ymm12,%%ymm8 ,%%ymm14	\n\t"\
	"vfmadd132pd	%%ymm12,%%ymm2,%%ymm5			\n\t	vfmadd132pd	%%ymm12,%%ymm10,%%ymm13	\n\t"\
	"vfmadd132pd	%%ymm12,%%ymm1,%%ymm7			\n\t	vfmadd132pd	%%ymm12,%%ymm9 ,%%ymm15	\n\t"\
	"vfmadd132pd	%%ymm12,%%ymm3,%%ymm4			\n\t	vfmadd132pd	(%%rsi),%%ymm11,%%ymm12	\n\t"\
		/* SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r30) */\
		"												\n\t	vmovaps	%%ymm14,0x100(%%rsi)	\n\t"\
		"												\n\t	vmovaps	%%ymm15,0x120(%%rsi)	\n\t"\
		"												\n\t	vmovaps	%%ymm10,%%ymm14	\n\t"\
		"												\n\t	vmovaps	%%ymm13,%%ymm15	\n\t"\
		"												\n\t	vmovaps	-0x20(%%r8 ),%%ymm8 	\n\t"/* one */\
		"														vfmsub132pd	%%ymm8 ,%%ymm12,%%ymm10	\n\t"\
		"														vfmsub132pd	%%ymm8 ,%%ymm11,%%ymm13	\n\t"\
		"														vfmadd132pd	%%ymm8 ,%%ymm12,%%ymm14	\n\t"\
		"														vfmadd132pd	%%ymm8 ,%%ymm11,%%ymm15	\n\t"\
		"												\n\t	vmovaps	0x200(%%rsi),%%ymm8 	\n\t"/* isrt2 */\
		"vmovaps	0x100(%%rsi),%%ymm11			\n\t	vfnmadd231pd %%ymm8,%%ymm10,%%ymm2	\n\t"\
		"vmovaps	0x120(%%rsi),%%ymm12			\n\t	vfnmadd231pd %%ymm8,%%ymm13,%%ymm3	\n\t"\
		"vfmsub132pd	-0x20(%%r8),%%ymm11,%%ymm6		\n\t	vfnmadd231pd %%ymm8,%%ymm14,%%ymm4	\n\t"\
		"vfmsub132pd	-0x20(%%r8),%%ymm9 ,%%ymm0		\n\t	vfnmadd231pd %%ymm8,%%ymm15,%%ymm5	\n\t"\
		"vsubpd	0x180(%%rsi),%%ymm1 ,%%ymm1 			\n\t	vfmsub132pd	-0x20(%%r8),%%ymm12 ,%%ymm7	\n\t"\
	"vmovaps	(%%r8),%%ymm8	\n\t"/* rcol: *2 => *sqrt2 due to above FMA leaving ymm10,13-15 unmultiplied-by-isrt2: */\
		"vmovaps	%%ymm6 ,0x100(%%rsi)				\n\t	vmovaps	%%ymm2 ,0x140(%%rsi)	\n\t"\
		"vmovaps	%%ymm0 ,0x080(%%rsi)				\n\t	vmovaps	%%ymm5 ,0x0c0(%%rsi)	\n\t"\
		"vmovaps	%%ymm7 ,0x120(%%rsi)				\n\t	vmovaps	%%ymm4 ,0x160(%%rsi)	\n\t"\
		"vmovaps	%%ymm1 ,0x1a0(%%rsi)				\n\t	vmovaps	%%ymm3 ,0x1e0(%%rsi)	\n\t"\
	"vfmadd132pd	%%ymm8 ,%%ymm6 ,%%ymm11			\n\t	vfmadd132pd	0x20(%%r8),%%ymm2 ,%%ymm10	\n\t"\
	"vfmadd132pd	%%ymm8 ,%%ymm0 ,%%ymm9 			\n\t	vfmadd132pd	0x20(%%r8),%%ymm5 ,%%ymm15	\n\t"\
	"vfmadd132pd	%%ymm8 ,%%ymm7 ,%%ymm12			\n\t	vfmadd132pd	0x20(%%r8),%%ymm4 ,%%ymm14	\n\t"\
	"vfmadd132pd	0x180(%%rsi),%%ymm1 ,%%ymm8 	\n\t	vfmadd132pd	0x20(%%r8),%%ymm3 ,%%ymm13	\n\t"\
		"vmovaps	%%ymm11,     (%%rsi)				\n\t	vmovaps	%%ymm10,0x040(%%rsi)	\n\t"\
		"vmovaps	%%ymm9 ,0x180(%%rsi)				\n\t	vmovaps	%%ymm15,0x1c0(%%rsi)	\n\t"\
		"vmovaps	%%ymm12,0x020(%%rsi)				\n\t	vmovaps	%%ymm14,0x060(%%rsi)	\n\t"\
		"vmovaps	%%ymm8 ,0x0a0(%%rsi)				\n\t	vmovaps	%%ymm13,0x0e0(%%rsi)	\n\t"\
		"\n\t"\
	/**********************************************************************************/\
	/*...and now do eight radix-4 transforms, including the internal twiddle factors: */\
	/**********************************************************************************/\
	/*...Block 1: t00,t10,t20,t30	*/\
		"movq	%[__add0],%%rax	/* &a[j1] */			\n\t"\
		"movslq	%[__p01],%%rbx							\n\t"	/*...Block 5: t08,t18,t28,t38	*/\
		"movslq	%[__p02],%%rcx							\n\t"	/* p04<<3 still in r9 */\
		"movslq	%[__p03],%%rdx							\n\t	leaq	(%%r9,%%rax),%%r10	\n\t"\
		"leaq	(%%rax,%%rbx,8),%%rbx					\n\t	leaq	(%%r9,%%rbx),%%r11	\n\t"\
		"leaq	(%%rax,%%rcx,8),%%rcx					\n\t	leaq	(%%r9,%%rcx),%%r12	\n\t"\
		"leaq	(%%rax,%%rdx,8),%%rdx					\n\t	leaq	(%%r9,%%rdx),%%r13	\n\t"\
		"movq	%[__r00],%%rsi							\n\t	vmovaps	0x800(%%rsi),%%ymm11	\n\t"/* isrt2 */\
		"vmovaps	     (%%rsi),%%ymm0					\n\t	vmovaps	0x500(%%rsi),%%ymm12	\n\t"\
		"vmovaps	0x400(%%rsi),%%ymm4					\n\t	vmovaps	0x520(%%rsi),%%ymm13	\n\t"\
		"vmovaps	0x020(%%rsi),%%ymm1					\n\t	vmovaps	0x700(%%rsi),%%ymm14	\n\t"\
		"vmovaps	0x420(%%rsi),%%ymm5					\n\t	vmovaps	0x720(%%rsi),%%ymm15	\n\t"\
		"vmovaps	0x200(%%rsi),%%ymm2					\n\t	vmulpd	%%ymm11,%%ymm12,%%ymm12	\n\t"\
		"vmovaps	0x600(%%rsi),%%ymm6					\n\t	vmulpd	%%ymm11,%%ymm13,%%ymm13	\n\t"\
		"vmovaps	0x220(%%rsi),%%ymm3					\n\t	vmulpd	%%ymm11,%%ymm14,%%ymm14	\n\t"\
		"vmovaps	0x620(%%rsi),%%ymm7					\n\t	vmulpd	%%ymm11,%%ymm15,%%ymm15	\n\t"\
		"vsubpd	0x200(%%rsi),%%ymm0,%%ymm0				\n\t	vmovaps	0x100(%%rsi),%%ymm8 	\n\t"\
		"vsubpd	0x600(%%rsi),%%ymm4,%%ymm4				\n\t	vmovaps	0x120(%%rsi),%%ymm9 	\n\t"\
		"vsubpd	0x220(%%rsi),%%ymm1,%%ymm1				\n\t	vmovaps	0x300(%%rsi),%%ymm10	\n\t"\
		"vsubpd	0x620(%%rsi),%%ymm5,%%ymm5				\n\t	vmovaps	0x320(%%rsi),%%ymm11	\n\t"\
	"vmovaps	%%ymm1 ,(%%rbx)	\n\t"/* spill ymm1  to make room for 1.0 */"vmovaps	-0x20(%%r8),%%ymm1 	\n\t"/* one */\
		"vaddpd	     (%%rsi),%%ymm2,%%ymm2				\n\t	vfmsub132pd	%%ymm1 ,%%ymm11,%%ymm8 	\n\t"\
		"vaddpd	0x400(%%rsi),%%ymm6,%%ymm6				\n\t	vfmsub132pd	%%ymm1 ,%%ymm13,%%ymm12	\n\t"\
		"vaddpd	0x020(%%rsi),%%ymm3,%%ymm3				\n\t	vfmsub132pd	%%ymm1 ,%%ymm10,%%ymm9 	\n\t"\
		"vaddpd	0x420(%%rsi),%%ymm7,%%ymm7				\n\t	vfmsub132pd	%%ymm1 ,%%ymm14,%%ymm15	\n\t"\
		"vfmsub132pd	%%ymm1 ,%%ymm6,%%ymm2		\n\t	vfmadd132pd	(%%r8) ,%%ymm8 ,%%ymm11	\n\t"\
		"vfmsub132pd	%%ymm1 ,%%ymm5,%%ymm0		\n\t	vfmadd132pd	(%%r8) ,%%ymm12,%%ymm13	\n\t"\
		"vfmsub132pd	%%ymm1 ,%%ymm7,%%ymm3		\n\t	vfmadd132pd	(%%r8) ,%%ymm9 ,%%ymm10	\n\t"\
		"vfmsub132pd	(%%rbx),%%ymm4,%%ymm1		\n\t	vfmadd132pd	(%%r8) ,%%ymm15,%%ymm14	\n\t"\
		"vmovaps	%%ymm2,     (%%rbx)					\n\t	vmovaps	(%%r8),%%ymm2	\n\t"/* 2.0 */\
	"vfmadd213pd	(%%rbx),%%ymm2,%%ymm6				\n\t	vfmsub132pd	-0x20(%%r8),%%ymm14,%%ymm12	\n\t"\
	"vfmadd132pd	%%ymm2 ,%%ymm0,%%ymm5				\n\t	vfmsub132pd	-0x20(%%r8),%%ymm15,%%ymm13	\n\t"\
	"vfmadd132pd	%%ymm2 ,%%ymm3,%%ymm7			\n\t	 vfmadd132pd	%%ymm2 ,%%ymm12,%%ymm14	\n\t"\
	"vfmadd132pd	%%ymm2 ,%%ymm1,%%ymm4			\n\t	 vfmadd132pd	%%ymm2 ,%%ymm13,%%ymm15	\n\t"\
	"vmovaps	%%ymm11,(%%rdx)	\n\t	vmovaps	-0x20(%%r8),%%ymm11	\n\t"/* spill to make room for 1.0 */\
		/* ymm2-datum already spilled to destination */"\n\t	vfmsub132pd	%%ymm11,%%ymm12,%%ymm8 	\n\t"\
		"vmovaps	%%ymm0,     (%%rcx)					\n\t	vfmsub132pd	%%ymm11,%%ymm13,%%ymm10	\n\t"\
		"vmovaps	%%ymm3,0x020(%%rbx)					\n\t	vfmsub132pd	%%ymm11,%%ymm14,%%ymm9 	\n\t"\
		"vmovaps	%%ymm1,0x020(%%rdx)					\n\t	vfmsub132pd	(%%rdx),%%ymm15,%%ymm11	\n\t"\
		"vmovaps	%%ymm6,     (%%rax)				\n\t	vfmadd132pd	%%ymm2 ,%%ymm8 ,%%ymm12	\n\t"\
		"vmovaps	%%ymm5,     (%%rdx)				\n\t	vfmadd132pd	%%ymm2 ,%%ymm11,%%ymm15	\n\t"\
		"vmovaps	%%ymm7,0x020(%%rax)				\n\t	vfmadd132pd	%%ymm2 ,%%ymm10,%%ymm13	\n\t"\
		"vmovaps	%%ymm4,0x020(%%rcx)				\n\t	vfmadd132pd	%%ymm2 ,%%ymm9 ,%%ymm14	\n\t"\
		"												\n\t	vmovaps	%%ymm8 ,     (%%r11)	\n\t"\
		"												\n\t	vmovaps	%%ymm11,     (%%r12)	\n\t"\
		"												\n\t	vmovaps	%%ymm10,0x020(%%r11)	\n\t"\
		"												\n\t	vmovaps	%%ymm9 ,0x020(%%r13)	\n\t"\
		"												\n\t	vmovaps	%%ymm12,     (%%r10)	\n\t"\
		"												\n\t	vmovaps	%%ymm15,     (%%r13)	\n\t"\
		"												\n\t	vmovaps	%%ymm13,0x020(%%r10)	\n\t"\
		"												\n\t	vmovaps	%%ymm14,0x020(%%r12)	\n\t"\
	/*...Block 3: t04,t14,t24,t34	*/\
		"addq	$0x080,%%rsi	/* r04 */				\n\t"\
		"subq	%%rax,%%rbx		\n\t"/* p01 << 3 ... Note the fact that the 3 subtracts leave the pointerfied (left- */\
		"subq	%%rax,%%rcx		\n\t"/* p02 << 3 ... shifted) offset means no ',8' (aka << 3) needed in last 3 LEAs. */\
		"subq	%%rax,%%rdx		\n\t"/* p03 << 3 */				/*...Block 7: t0C,t1C,t2C,t3C	*/\
		"movslq	%[__p08],%%rdi							\n\t"	/* p04<<3 still in r9 */\
		"leaq	(%%rax,%%rdi,8),%%rax	/* &a[j1+p8] */	\n\t	leaq	(%%r9,%%rax),%%r10	\n\t"\
		"leaq	(%%rax,%%rbx),%%rbx						\n\t	leaq	(%%r9,%%rbx),%%r11	\n\t"\
		"leaq	(%%rax,%%rcx),%%rcx						\n\t	leaq	(%%r9,%%rcx),%%r12	\n\t"\
		"leaq	(%%rax,%%rdx),%%rdx						\n\t	leaq	(%%r9,%%rdx),%%r13	\n\t"\
		"vmovaps	0x400(%%rsi),%%ymm4					\n\t	vmovaps	0x500(%%rsi),%%ymm12	\n\t"\
		"vmovaps	0x420(%%rsi),%%ymm5					\n\t	vmovaps	0x520(%%rsi),%%ymm13	\n\t"\
		"vmovaps	0x7a0(%%rsi),%%ymm3	/* cc0 */		\n\t"\
		"vmovaps	0x7c0(%%rsi),%%ymm2					\n\t"\
		"vmovaps	%%ymm4,%%ymm6						\n\t	vmovaps	%%ymm12,%%ymm14	\n\t"\
		"vmovaps	%%ymm5,%%ymm7						\n\t	vmovaps	%%ymm13,%%ymm15	\n\t"\
		"vmulpd	%%ymm3,%%ymm4,%%ymm4					\n\t	vmulpd	%%ymm2 ,%%ymm12,%%ymm12	\n\t"\
		"vmulpd	%%ymm3,%%ymm5,%%ymm5					\n\t	vmulpd	%%ymm2 ,%%ymm13,%%ymm13	\n\t"\
		"vmovaps	0x600(%%rsi),%%ymm0					\n\t	vmovaps	0x700(%%rsi),%%ymm8 	\n\t"\
		"vmovaps	0x620(%%rsi),%%ymm1					\n\t	vmovaps	0x720(%%rsi),%%ymm9 	\n\t"\
	"vfnmadd231pd	%%ymm2,%%ymm7,%%ymm4			\n\t	vfnmadd231pd	%%ymm3 ,%%ymm15,%%ymm12	\n\t"\
	" vfmadd231pd	%%ymm2,%%ymm6,%%ymm5			\n\t	 vfmadd231pd	%%ymm3 ,%%ymm14,%%ymm13	\n\t"\
		"vmovaps	%%ymm0,%%ymm6						\n\t	vmovaps	%%ymm8 ,%%ymm14	\n\t"\
		"vmovaps	%%ymm1,%%ymm7						\n\t	vmovaps	%%ymm9 ,%%ymm15	\n\t"\
		"vmulpd	%%ymm2,%%ymm6,%%ymm6					\n\t	vmulpd	%%ymm3 ,%%ymm14,%%ymm14	\n\t"\
		"vmulpd	%%ymm2,%%ymm7,%%ymm7					\n\t	vmulpd	%%ymm3 ,%%ymm15,%%ymm15	\n\t"\
	"vfnmadd231pd	%%ymm3,%%ymm1,%%ymm6			\n\t	vfnmadd231pd	%%ymm2 ,%%ymm9 ,%%ymm14	\n\t"\
	" vfmadd231pd	%%ymm3,%%ymm0,%%ymm7			\n\t	 vfmadd231pd	%%ymm2 ,%%ymm8 ,%%ymm15	\n\t"\
		"vmovaps	%%ymm4,%%ymm2						\n\t	vmovaps	%%ymm12,%%ymm10	\n\t"\
		"vmovaps	%%ymm5,%%ymm3						\n\t	vmovaps	%%ymm13,(%%rbx)	\n\t	vmovaps	-0x20(%%r8),%%ymm11	\n\t"/* one */\
		"vfmsub132pd	%%ymm11,%%ymm6,%%ymm4			\n\t	vfmsub132pd	%%ymm11,%%ymm14,%%ymm12	\n\t"\
		"vfmsub132pd	%%ymm11,%%ymm7,%%ymm5			\n\t	vfmsub132pd	%%ymm11,%%ymm15,%%ymm13	\n\t"\
		"vfmadd132pd	%%ymm11,%%ymm2,%%ymm6			\n\t	vfmadd132pd	%%ymm11,%%ymm10,%%ymm14	\n\t"\
		"vfmadd132pd	%%ymm11,%%ymm3,%%ymm7			\n\t	vaddpd	(%%rbx),%%ymm15,%%ymm15	\n\t"\
		"vmovaps	0x200(%%rsi),%%ymm2					\n\t	vmovaps	0x300(%%rsi),%%ymm10	\n\t"\
		"vmovaps	0x220(%%rsi),%%ymm3					\n\t	vmovaps	0x320(%%rsi),%%ymm11	\n\t"\
		"vmovaps	0x780(%%rsi),%%ymm1	\n\t"/* isrt2 */\
		"vmovaps	%%ymm2,%%ymm0						\n\t	vmovaps	%%ymm10,%%ymm8 	\n\t"\
	"vmovaps	%%ymm14,(%%rbx)	\n\t	vmovaps	-0x20(%%r8),%%ymm14	\n\t"/* spill to make room for 1.0 */\
		"vfmsub132pd	%%ymm14,%%ymm3,%%ymm2					\n\t	vfmadd132pd	%%ymm14,%%ymm11,%%ymm10	\n\t"\
		"vfmadd132pd	%%ymm14,%%ymm0,%%ymm3					\n\t	vfmsub132pd	%%ymm14,%%ymm8 ,%%ymm11	\n\t"\
		"vmulpd	%%ymm1,%%ymm2,%%ymm2					\n\t	vmulpd	%%ymm1 ,%%ymm10,%%ymm10	\n\t"\
		"vmulpd	%%ymm1,%%ymm3,%%ymm3					\n\t	vmulpd	%%ymm1 ,%%ymm11,%%ymm11	\n\t"\
		"vmovaps	     (%%rsi),%%ymm0					\n\t	vmovaps	0x100(%%rsi),%%ymm8 	\n\t"\
		"vmovaps	0x020(%%rsi),%%ymm1					\n\t	vmovaps	0x120(%%rsi),%%ymm9 	\n\t"\
		"vfmsub132pd	%%ymm14,%%ymm2,%%ymm0			\n\t	vfmsub132pd	%%ymm14,%%ymm10,%%ymm8 	\n\t"\
		"vfmsub132pd	%%ymm14,%%ymm3,%%ymm1			\n\t	vfmsub132pd	%%ymm14,%%ymm11,%%ymm9 	\n\t"\
	"vfmadd132pd	(%%r8),%%ymm0,%%ymm2			\n\t	vfmadd132pd	(%%r8) ,%%ymm8 ,%%ymm10	\n\t"\
	"vfmadd132pd	(%%r8),%%ymm1,%%ymm3			\n\t	vfmadd132pd	(%%r8) ,%%ymm9 ,%%ymm11	\n\t"\
		"vfmsub132pd	%%ymm14,%%ymm6,%%ymm2			\n\t	vfmsub132pd	%%ymm14,%%ymm12,%%ymm8 	\n\t"\
		"vfmsub132pd	%%ymm14,%%ymm5,%%ymm0			\n\t	vfmsub132pd	%%ymm14,%%ymm15,%%ymm10	\n\t"\
		"vfmsub132pd	%%ymm14,%%ymm7,%%ymm3			\n\t	vfmsub132pd	%%ymm14,%%ymm13,%%ymm9 	\n\t"\
		"vfmsub132pd	%%ymm14,%%ymm4,%%ymm1			\n\t	vfmsub213pd	(%%rbx),%%ymm14,%%ymm11	\n\t"\
	"vmovaps	(%%rbx),%%ymm14	\n\t"/* restore spill */\
		"vmovaps	%%ymm2,     (%%rbx)					\n\t	vmovaps	%%ymm8 ,     (%%r11)	\n\t"\
		"vmovaps	%%ymm0,     (%%rcx)					\n\t	vmovaps	%%ymm10,     (%%r12)	\n\t"\
		"vmovaps	%%ymm3,0x020(%%rbx)					\n\t	vmovaps	%%ymm9 ,0x020(%%r11)	\n\t"\
		"vmovaps	%%ymm1,0x020(%%rdx)					\n\t	vmovaps	%%ymm11,0x020(%%r13)	\n\t"\
	"vmovaps	(%%r8),%%ymm2	\n\t"/* 2.0 */\
	"vfmadd213pd	(%%rbx),%%ymm2,%%ymm6			\n\t	vfmadd132pd	%%ymm2 ,%%ymm8 ,%%ymm12	\n\t"\
	"vfmadd132pd	%%ymm2 ,%%ymm0,%%ymm5			\n\t	vfmadd132pd	%%ymm2 ,%%ymm10,%%ymm15	\n\t"\
	"vfmadd132pd	%%ymm2 ,%%ymm3,%%ymm7			\n\t	vfmadd132pd	%%ymm2 ,%%ymm9 ,%%ymm13	\n\t"\
	"vfmadd132pd	%%ymm2 ,%%ymm1,%%ymm4			\n\t	vfmadd132pd	%%ymm2 ,%%ymm11,%%ymm14	\n\t"\
		"vmovaps	%%ymm6,     (%%rax)					\n\t	vmovaps	%%ymm12,     (%%r10)	\n\t"\
		"vmovaps	%%ymm5,     (%%rdx)					\n\t	vmovaps	%%ymm15,     (%%r13)	\n\t"\
		"vmovaps	%%ymm7,0x020(%%rax)					\n\t	vmovaps	%%ymm13,0x020(%%r10)	\n\t"\
		"vmovaps	%%ymm4,0x020(%%rcx)					\n\t	vmovaps	%%ymm14,0x020(%%r12)	\n\t"\
		"\n\t"\
	/*...Block 2: t02,t12,t22,t32	*/\
		"subq	$0x040,%%rsi	/* r02 */				\n\t"\
		"subq	%%rax,%%rbx								\n\t"\
		"subq	%%rax,%%rcx								\n\t"\
		"subq	%%rax,%%rdx								\n\t"\
		"movq	%[__add0],%%rax	/* &a[j1] */			\n\t"	/*...Block 6: t0A,t1A,t2A,t3A	*/\
		"movslq	%[__p10],%%rdi							\n\t"	/* p04<<3 still in r9 */\
		"leaq	(%%rax,%%rdi,8),%%rax	/* &a[j1+p10] */\n\t	leaq	(%%r9,%%rax),%%r10	\n\t"\
		"leaq	(%%rax,%%rbx),%%rbx						\n\t	leaq	(%%r9,%%rbx),%%r11	\n\t"\
		"leaq	(%%rax,%%rcx),%%rcx						\n\t	leaq	(%%r9,%%rcx),%%r12	\n\t"\
		"leaq	(%%rax,%%rdx),%%rdx						\n\t	leaq	(%%r9,%%rdx),%%r13	\n\t"\
		"vmovaps	0x400(%%rsi),%%ymm4					\n\t	vmovaps	0x500(%%rsi),%%ymm12	\n\t"\
		"vmovaps	0x420(%%rsi),%%ymm5					\n\t	vmovaps	0x520(%%rsi),%%ymm13	\n\t"\
		"vmovaps	0x820(%%rsi),%%ymm2	/* cc1 */		\n\t	vmovaps	0x860(%%rsi),%%ymm11	/* cc3 */	\n\t"\
		"vmovaps	0x840(%%rsi),%%ymm3					\n\t	vmovaps	0x880(%%rsi),%%ymm10	\n\t"\
		"vmovaps	%%ymm4,%%ymm6						\n\t	vmovaps	%%ymm12,%%ymm14	\n\t"\
		"vmovaps	%%ymm5,%%ymm7						\n\t	vmovaps	%%ymm13,%%ymm15	\n\t"\
		"vmulpd	%%ymm2,%%ymm4,%%ymm4					\n\t	vmulpd	%%ymm10,%%ymm12,%%ymm12	\n\t"\
		"vmulpd	%%ymm2,%%ymm5,%%ymm5					\n\t	vmulpd	%%ymm10,%%ymm13,%%ymm13	\n\t"\
		"vmovaps	0x600(%%rsi),%%ymm0					\n\t	vmovaps	0x700(%%rsi),%%ymm8 	\n\t"\
		"vmovaps	0x620(%%rsi),%%ymm1					\n\t	vmovaps	0x720(%%rsi),%%ymm9 	\n\t"\
	"vfnmadd231pd	%%ymm3,%%ymm7,%%ymm4			\n\t	vfnmadd231pd	%%ymm11,%%ymm15,%%ymm12	\n\t"\
	" vfmadd231pd	%%ymm3,%%ymm6,%%ymm5			\n\t	 vfmadd231pd	%%ymm11,%%ymm14,%%ymm13	\n\t"\
		"vmovaps	%%ymm0,%%ymm6						\n\t	vmovaps	%%ymm8 ,%%ymm14	\n\t"\
		"vmovaps	%%ymm1,%%ymm7						\n\t	vmovaps	%%ymm9 ,%%ymm15	\n\t"\
		"vmulpd	%%ymm11,%%ymm6,%%ymm6					\n\t	vmulpd	%%ymm2 ,%%ymm14,%%ymm14	\n\t"\
		"vmulpd	%%ymm11,%%ymm7,%%ymm7					\n\t	vmulpd	%%ymm2 ,%%ymm15,%%ymm15	\n\t"\
	"vfnmadd231pd	%%ymm10,%%ymm1,%%ymm6			\n\t	 vfmadd231pd	%%ymm3 ,%%ymm9 ,%%ymm14	\n\t"\
	" vfmadd231pd	%%ymm10,%%ymm0,%%ymm7			\n\t	vfnmadd231pd	%%ymm3 ,%%ymm8 ,%%ymm15	\n\t"\
	"vmovaps	-0x20(%%r8),%%ymm1 	\n\t"/* ymm1 free; use for 1.0 */\
		"vmovaps	%%ymm4,%%ymm2						\n\t	vmovaps	%%ymm12,%%ymm10	\n\t"\
		"vmovaps	%%ymm5,%%ymm3						\n\t	vmovaps	%%ymm13,%%ymm11	\n\t"\
		"vfmsub132pd	%%ymm1 ,%%ymm6,%%ymm4			\n\t	vfmsub132pd	%%ymm1 ,%%ymm14,%%ymm12	\n\t"\
		"vfmsub132pd	%%ymm1 ,%%ymm7,%%ymm5			\n\t	vfmsub132pd	%%ymm1 ,%%ymm15,%%ymm13	\n\t"\
		"vfmadd132pd	%%ymm1 ,%%ymm2,%%ymm6			\n\t	vfmadd132pd	%%ymm1 ,%%ymm10,%%ymm14	\n\t"\
		"vfmadd132pd	%%ymm1 ,%%ymm3,%%ymm7			\n\t	vfmadd132pd	%%ymm1 ,%%ymm11,%%ymm15	\n\t"\
		"vmovaps	0x200(%%rsi),%%ymm1					\n\t	vmovaps	0x300(%%rsi),%%ymm9 	\n\t"\
		"vmovaps	0x220(%%rsi),%%ymm3					\n\t	vmovaps	0x320(%%rsi),%%ymm11	\n\t"\
		"vmovaps	0x800(%%rsi),%%ymm0	/* ss0 */		\n\t	vmovaps	0x7e0(%%rsi),%%ymm8 		/* cc0 */	\n\t"\
		"vmovaps	%%ymm1,%%ymm2						\n\t	vmovaps	%%ymm9 ,%%ymm10	\n\t"\
		"vmulpd	      %%ymm0,%%ymm1,%%ymm1				\n\t	vmulpd	     %%ymm8 ,%%ymm9 ,%%ymm9 	\n\t"\
		"vmulpd	      %%ymm3,%%ymm0,%%ymm0				\n\t	vmulpd	     %%ymm11,%%ymm8 ,%%ymm8 	\n\t"\
	"vfmsub132pd	0x7e0(%%rsi),%%ymm0,%%ymm2		\n\t	vfmadd132pd	0x800(%%rsi),%%ymm8 ,%%ymm10	\n\t"\
	"vfmadd132pd	0x7e0(%%rsi),%%ymm1,%%ymm3		\n\t	vfmsub132pd	0x800(%%rsi),%%ymm9 ,%%ymm11	\n\t"\
	"vmovaps	%%ymm14,(%%rbx)	\n\t	vmovaps	-0x20(%%r8),%%ymm14	\n\t"/* spill to make room for 1.0 */\
		"vmovaps	     (%%rsi),%%ymm0					\n\t	vmovaps	0x100(%%rsi),%%ymm8 	\n\t"\
		"vmovaps	0x020(%%rsi),%%ymm1					\n\t	vmovaps	0x120(%%rsi),%%ymm9 	\n\t"\
		"vfmsub132pd	%%ymm14,%%ymm2,%%ymm0			\n\t	vfmsub132pd	%%ymm14,%%ymm10,%%ymm8 	\n\t"\
		"vfmsub132pd	%%ymm14,%%ymm3,%%ymm1			\n\t	vfmsub132pd	%%ymm14,%%ymm11,%%ymm9 	\n\t"\
		"vfmadd132pd	(%%r8),%%ymm0,%%ymm2			\n\t	vfmadd132pd	(%%r8) ,%%ymm8 ,%%ymm10	\n\t"\
		"vfmadd132pd	(%%r8),%%ymm1,%%ymm3			\n\t	vfmadd132pd	(%%r8) ,%%ymm9 ,%%ymm11	\n\t"\
		"vfmsub132pd	%%ymm14,%%ymm6,%%ymm2			\n\t	vfmsub132pd	%%ymm14,%%ymm12,%%ymm8 	\n\t"\
		"vfmsub132pd	%%ymm14,%%ymm5,%%ymm0			\n\t	vfmsub132pd	%%ymm14,%%ymm15,%%ymm10	\n\t"\
		"vfmsub132pd	%%ymm14,%%ymm7,%%ymm3			\n\t	vfmsub132pd	%%ymm14,%%ymm13,%%ymm9 	\n\t"\
		"vfmsub132pd	%%ymm14,%%ymm4,%%ymm1			\n\t	vfmsub213pd	(%%rbx),%%ymm14,%%ymm11	\n\t"\
		"vmovaps	(%%rbx),%%ymm14	\n\t"/* restore spill */\
		"vmovaps	%%ymm2,     (%%rbx)					\n\t	vmovaps	%%ymm8 ,     (%%r11)	\n\t"\
		"vmovaps	%%ymm0,     (%%rcx)					\n\t	vmovaps	%%ymm10,     (%%r12)	\n\t"\
		"vmovaps	%%ymm3,0x020(%%rbx)					\n\t	vmovaps	%%ymm9 ,0x020(%%r11)	\n\t"\
		"vmovaps	%%ymm1,0x020(%%rdx)					\n\t	vmovaps	%%ymm11,0x020(%%r13)	\n\t"\
	"vmovaps	(%%r8),%%ymm2	\n\t"/* 2.0 */\
	"vfmadd213pd	(%%rbx),%%ymm2,%%ymm6			\n\t	vfmadd132pd	%%ymm2 ,%%ymm8 ,%%ymm12	\n\t"\
	"vfmadd132pd	%%ymm2 ,%%ymm0,%%ymm5			\n\t	vfmadd132pd	%%ymm2 ,%%ymm10,%%ymm15	\n\t"\
	"vfmadd132pd	%%ymm2 ,%%ymm3,%%ymm7			\n\t	vfmadd132pd	%%ymm2 ,%%ymm9 ,%%ymm13	\n\t"\
	"vfmadd132pd	%%ymm2 ,%%ymm1,%%ymm4			\n\t	vfmadd132pd	%%ymm2 ,%%ymm11,%%ymm14	\n\t"\
		"vmovaps	%%ymm6,     (%%rax)					\n\t	vmovaps	%%ymm12,     (%%r10)	\n\t"\
		"vmovaps	%%ymm5,     (%%rdx)					\n\t	vmovaps	%%ymm15,     (%%r13)	\n\t"\
		"vmovaps	%%ymm7,0x020(%%rax)					\n\t	vmovaps	%%ymm13,0x020(%%r10)	\n\t"\
		"vmovaps	%%ymm4,0x020(%%rcx)					\n\t	vmovaps	%%ymm14,0x020(%%r12)	\n\t"\
		"\n\t"\
	/*...Block 4: t06,t16,t26,t36	*/\
		"addq	$0x080,%%rsi	/* r06 */	\n\t"\
		"subq	%%rax,%%rbx								\n\t"\
		"subq	%%rax,%%rcx								\n\t"\
		"subq	%%rax,%%rdx								\n\t"\
		"movq	%[__add0],%%rax	/* &a[j1] */			\n\t"	/*...Block 8: t0E,t1E,t2E,t3E	*/\
		"movslq	%[__p18],%%rdi							\n\t"	/* p04<<3 still in r9 */\
		"leaq	(%%rax,%%rdi,8),%%rax	/* &a[j1+p18] */\n\t	leaq	(%%r9,%%rax),%%r10	\n\t"\
		"leaq	(%%rax,%%rbx),%%rbx						\n\t	leaq	(%%r9,%%rbx),%%r11	\n\t"\
		"leaq	(%%rax,%%rcx),%%rcx						\n\t	leaq	(%%r9,%%rcx),%%r12	\n\t"\
		"leaq	(%%rax,%%rdx),%%rdx						\n\t	leaq	(%%r9,%%rdx),%%r13	\n\t"\
		"vmovaps	0x400(%%rsi),%%ymm4					\n\t	vmovaps	0x500(%%rsi),%%ymm12	\n\t"\
		"vmovaps	0x420(%%rsi),%%ymm5					\n\t	vmovaps	0x520(%%rsi),%%ymm13	\n\t"\
		"vmovaps	0x7e0(%%rsi),%%ymm2	/* cc3 */		\n\t	vmovaps	0x7a0(%%rsi),%%ymm11 /* cc1 */	\n\t"\
		"vmovaps	0x800(%%rsi),%%ymm3					\n\t	vmovaps	0x7c0(%%rsi),%%ymm10 \n\t"\
		"vmovaps	%%ymm4,%%ymm6						\n\t	vmovaps	%%ymm12,%%ymm14	\n\t"\
		"vmovaps	%%ymm5,%%ymm7						\n\t	vmovaps	%%ymm13,%%ymm15	\n\t"\
		"vmulpd	%%ymm2,%%ymm4,%%ymm4					\n\t	vmulpd	%%ymm10,%%ymm12,%%ymm12	\n\t"\
		"vmulpd	%%ymm2,%%ymm5,%%ymm5					\n\t	vmulpd	%%ymm10,%%ymm13,%%ymm13	\n\t"\
	"vfnmadd231pd	%%ymm3,%%ymm7,%%ymm4			\n\t	vfnmadd231pd	%%ymm11,%%ymm15,%%ymm12	\n\t"\
	" vfmadd231pd	%%ymm3,%%ymm6,%%ymm5			\n\t	 vfmadd231pd	%%ymm11,%%ymm14,%%ymm13	\n\t"\
		"vmovaps	0x600(%%rsi),%%ymm0					\n\t	vmovaps	0x700(%%rsi),%%ymm8 	\n\t"\
		"vmovaps	0x620(%%rsi),%%ymm1					\n\t	vmovaps	0x720(%%rsi),%%ymm9 	\n\t"\
		"vmovaps	%%ymm0,%%ymm6						\n\t	vmovaps	%%ymm8 ,%%ymm14	\n\t"\
		"vmovaps	%%ymm1,%%ymm7						\n\t	vmovaps	%%ymm9 ,%%ymm15	\n\t"\
		"vmulpd	%%ymm10,%%ymm6,%%ymm6					\n\t	vmulpd	%%ymm3 ,%%ymm14,%%ymm14	\n\t"\
		"vmulpd	%%ymm10,%%ymm7,%%ymm7					\n\t	vmulpd	%%ymm3 ,%%ymm15,%%ymm15	\n\t"\
	" vfmadd231pd	%%ymm11,%%ymm1,%%ymm6			\n\t	vfnmadd231pd	%%ymm2 ,%%ymm9 ,%%ymm14	\n\t"\
	"vfnmadd231pd	%%ymm11,%%ymm0,%%ymm7			\n\t	 vfmadd231pd	%%ymm2 ,%%ymm8 ,%%ymm15	\n\t"\
	"vmovaps	-0x20(%%r8),%%ymm1 	\n\t"/* ymm1 free; use for 1.0 */\
		"vmovaps	%%ymm4,%%ymm2						\n\t	vmovaps	%%ymm12,%%ymm10	\n\t"\
		"vmovaps	%%ymm5,%%ymm3						\n\t	vmovaps	%%ymm13,%%ymm11	\n\t"\
		"vfmsub132pd	%%ymm1 ,%%ymm6,%%ymm4			\n\t	vfmsub132pd	%%ymm1 ,%%ymm14,%%ymm12	\n\t"\
		"vfmsub132pd	%%ymm1 ,%%ymm7,%%ymm5			\n\t	vfmsub132pd	%%ymm1 ,%%ymm15,%%ymm13	\n\t"\
		"vfmadd132pd	%%ymm1 ,%%ymm2,%%ymm6			\n\t	vfmadd132pd	%%ymm1 ,%%ymm10,%%ymm14	\n\t"\
		"vfmadd132pd	%%ymm1 ,%%ymm3,%%ymm7			\n\t	vfmadd132pd	%%ymm1 ,%%ymm11,%%ymm15	\n\t"\
		"vmovaps	0x200(%%rsi),%%ymm1					\n\t	vmovaps	0x300(%%rsi),%%ymm9 	\n\t"\
		"vmovaps	0x220(%%rsi),%%ymm3					\n\t	vmovaps	0x320(%%rsi),%%ymm11	\n\t"\
		"vmovaps	0x760(%%rsi),%%ymm0		/* cc0 */	\n\t	vmovaps	0x780(%%rsi),%%ymm8 	/* ss0 */	\n\t"\
		"vmovaps	%%ymm1,%%ymm2						\n\t	vmovaps	%%ymm9 ,%%ymm10	\n\t"\
		"vmulpd	      %%ymm0,%%ymm1,%%ymm1				\n\t	vmulpd	     %%ymm8 ,%%ymm9 ,%%ymm9 	\n\t"\
		"vmulpd	      %%ymm3,%%ymm0,%%ymm0				\n\t	vmulpd	     %%ymm11,%%ymm8 ,%%ymm8 	\n\t"\
	"vfmsub132pd	0x780(%%rsi),%%ymm0,%%ymm2		\n\t	vfmadd132pd	0x760(%%rsi),%%ymm8 ,%%ymm10	\n\t"\
	"vfmadd132pd	0x780(%%rsi),%%ymm1,%%ymm3		\n\t	vfmsub132pd	0x760(%%rsi),%%ymm9 ,%%ymm11	\n\t"\
	"vmovaps	%%ymm14,(%%rbx)	\n\t	vmovaps	-0x20(%%r8),%%ymm14	\n\t"/* spill to make room for 1.0 */\
		"vmovaps	     (%%rsi),%%ymm0					\n\t	vmovaps	0x100(%%rsi),%%ymm8 	\n\t"\
		"vmovaps	0x020(%%rsi),%%ymm1					\n\t	vmovaps	0x120(%%rsi),%%ymm9 	\n\t"\
		"vfmsub132pd	%%ymm14,%%ymm2,%%ymm0			\n\t	vfmsub132pd	%%ymm14,%%ymm10,%%ymm8 	\n\t"\
		"vfmsub132pd	%%ymm14,%%ymm3,%%ymm1			\n\t	vfmsub132pd	%%ymm14,%%ymm11,%%ymm9 	\n\t"\
	"vfmadd132pd	(%%r8),%%ymm0,%%ymm2			\n\t	vfmadd132pd	(%%r8) ,%%ymm8 ,%%ymm10	\n\t"\
	"vfmadd132pd	(%%r8),%%ymm1,%%ymm3			\n\t	vfmadd132pd	(%%r8) ,%%ymm9 ,%%ymm11	\n\t"\
		"vfmsub132pd	%%ymm14,%%ymm4,%%ymm2			\n\t	vfmsub132pd	%%ymm14,%%ymm12,%%ymm8 	\n\t"\
		"vfmsub132pd	%%ymm14,%%ymm7,%%ymm0			\n\t	vfmsub132pd	%%ymm14,%%ymm15,%%ymm10	\n\t"\
		"vfmsub132pd	%%ymm14,%%ymm5,%%ymm3			\n\t	vfmsub132pd	%%ymm14,%%ymm13,%%ymm9 	\n\t"\
		"vfmsub132pd	%%ymm14,%%ymm6,%%ymm1			\n\t	vfmsub213pd	(%%rbx),%%ymm14,%%ymm11	\n\t"\
	"vmovaps	(%%rbx),%%ymm14	\n\t"/* restore spill */\
		"vmovaps	%%ymm2,     (%%rbx)					\n\t	vmovaps	%%ymm8 ,     (%%r11)	\n\t"\
		"vmovaps	%%ymm0,     (%%rcx)					\n\t	vmovaps	%%ymm10,     (%%r12)	\n\t"\
		"vmovaps	%%ymm3,0x020(%%rbx)					\n\t	vmovaps	%%ymm9 ,0x020(%%r11)	\n\t"\
		"vmovaps	%%ymm1,0x020(%%rdx)					\n\t	vmovaps	%%ymm11,0x020(%%r13)	\n\t"\
	"vmovaps	(%%r8),%%ymm2	\n\t"/* 2.0 */\
	"vfmadd213pd	(%%rbx),%%ymm2,%%ymm4			\n\t	vfmadd132pd	%%ymm2 ,%%ymm8 ,%%ymm12	\n\t"\
	"vfmadd132pd	%%ymm2 ,%%ymm0,%%ymm7			\n\t	vfmadd132pd	%%ymm2 ,%%ymm10,%%ymm15	\n\t"\
	"vfmadd132pd	%%ymm2 ,%%ymm3,%%ymm5			\n\t	vfmadd132pd	%%ymm2 ,%%ymm9 ,%%ymm13	\n\t"\
	"vfmadd132pd	%%ymm2 ,%%ymm1,%%ymm6			\n\t	vfmadd132pd	%%ymm2 ,%%ymm11,%%ymm14	\n\t"\
		"vmovaps	%%ymm4,     (%%rax)					\n\t	vmovaps	%%ymm12,     (%%r10)	\n\t"\
		"vmovaps	%%ymm7,     (%%rdx)					\n\t	vmovaps	%%ymm15,     (%%r13)	\n\t"\
		"vmovaps	%%ymm5,0x020(%%rax)					\n\t	vmovaps	%%ymm13,0x020(%%r10)	\n\t"\
		"vmovaps	%%ymm6,0x020(%%rcx)					\n\t	vmovaps	%%ymm14,0x020(%%r12)	\n\t"\
		"\n\t"\
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
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r8","r9","r10","r11","r12","r13","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"		/* Clobbered registers */\
	);\
	}

	/*
	For GCC-macro version of this, use that - in terms of vec_dbl pointer-offsets - isrt2 + 1,3,5 = cc0,cc1,cc3,
	isrt2 + 7,15,23,31,39,47,55,63 = c00,04,02,06,01,05,03,07, and isrt2 + 71,72 = one,two
	in order to reduce number of args to <= the GCC-allowed maximum of 30:
	*/
	// In our FMA-ized version, original mix of [220 ADD, 220 SUB, 374 MUL] ==> [0 ADD, 5 SUB, 435 FMA (272 nontrivial), 102 MUL].
	#define SSE2_RADIX32_DIT_TWIDDLE(Xadd0,Xp01,Xp02,Xp03,Xp04,Xp05,Xp06,Xp07,Xp08,Xp10,Xp18,Xr00,Xr10,Xr20,Xr30,Xisrt2)\
	{\
	__asm__ volatile (\
	/*...Block 1: */\
		"movq	%[__r00],%%rsi					\n\t"\
		"movq	%[__add0],%%rax					\n\t"\
		"movslq	%[__p01],%%rbx					\n\t"\
		"movslq	%[__p02],%%rcx					\n\t		movq	%[__isrt2],%%r8	\n\t"\
		"movslq	%[__p03],%%rdx					\n\t		movslq	%[__p04],%%r9	\n\t"\
		"shlq	$3,%%rbx						\n\t		shlq	$3,%%r9		\n\t"\
		"shlq	$3,%%rcx						\n\t		addq	$0x900,%%r8	\n\t"/* two */\
		"shlq	$3,%%rdx						\n\t		leaq	(%%r9,%%rax),%%r10	\n\t"\
		"addq	%%rax,%%rbx						\n\t		leaq	(%%r9,%%rbx),%%r11	\n\t"\
		"addq	%%rax,%%rcx						\n\t		leaq	(%%r9,%%rcx),%%r12	\n\t"\
		"addq	%%rax,%%rdx						\n\t		leaq	(%%r9,%%rdx),%%r13	\n\t"\
		/* SSE2_RADIX4_DIT_0TWIDDLE_B(r00) */			/* SSE2_RADIX4_DIT_0TWIDDLE_2NDOFTWO_B(r08) */\
		"vmovaps	     (%%rax),%%ymm2			\n\t		vmovaps	     (%%r10),%%ymm10		\n\t"\
		"vmovaps	     (%%rcx),%%ymm6			\n\t		vmovaps	     (%%r12),%%ymm14		\n\t"\
		"vmovaps	0x020(%%rax),%%ymm3			\n\t		vmovaps	0x020(%%r10),%%ymm11		\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm7			\n\t"	/*	vmovaps	0x020(%%r12),%%ymm15	*/\
		"vmovaps	     (%%rbx),%%ymm0			\n\t		vmovaps	     (%%r11),%%ymm8 		\n\t"\
		"vmovaps	     (%%rdx),%%ymm4			\n\t		vmovaps	     (%%r13),%%ymm12		\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm1			\n\t		vmovaps	0x020(%%r11),%%ymm9 		\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm5			\n\t		vmovaps	0x020(%%r13),%%ymm13		\n\t"\
	"vmovaps	-0x20(%%r8),%%ymm15	\n\t"/* ymm15 free; use for 1.0 */\
		"vfmsub132pd	%%ymm15,%%ymm0,%%ymm2			\n\t		vfmsub132pd	%%ymm15,%%ymm8 ,%%ymm10			\n\t"\
		"vfmsub132pd	%%ymm15,%%ymm4,%%ymm6			\n\t		vfmsub132pd	%%ymm15,%%ymm12,%%ymm14			\n\t"\
		"vfmsub132pd	%%ymm15,%%ymm1,%%ymm3			\n\t		vfmsub132pd	%%ymm15,%%ymm9 ,%%ymm11			\n\t"\
		"vfmsub132pd	%%ymm15,%%ymm5,%%ymm7			\n\t		vfmsub132pd	0x20(%%r12),%%ymm13,%%ymm15		\n\t"\
	"vmovaps	%%ymm13,(%%rsi) 	\n\t"/* spill ymm13 to make room for 2.0 */"	vmovaps	(%%r8),%%ymm13 	\n\t"/* two */\
	"vfmadd132pd	%%ymm13,%%ymm2,%%ymm0		\n\t	vfmadd132pd	%%ymm13,%%ymm10,%%ymm8 	\n\t"\
	"vfmadd132pd	%%ymm13,%%ymm6,%%ymm4		\n\t	vfmadd132pd	%%ymm13,%%ymm14,%%ymm12	\n\t"\
	"vfmadd132pd	%%ymm13,%%ymm3,%%ymm1		\n\t	vfmadd132pd	%%ymm13,%%ymm11,%%ymm9 	\n\t"\
	"vfmadd132pd	%%ymm13,%%ymm7,%%ymm5		\n\t	vfmadd132pd	(%%rsi),%%ymm15,%%ymm13	\n\t"\
		/* Finish radix-4 butterfly: */\
	"vmovaps	%%ymm14,(%%rsi) 	\n\t"/* spill ymm11 to make room for 1.0 */"vmovaps	-0x20(%%r8),%%ymm14 \n\t"/* one */\
		"vfmsub132pd	%%ymm14,%%ymm4,%%ymm0			\n\t		vfmsub132pd	%%ymm14,%%ymm12,%%ymm8 			\n\t"\
		"vfmsub132pd	%%ymm14,%%ymm7,%%ymm2			\n\t		vfmsub132pd	%%ymm14,%%ymm13,%%ymm9 			\n\t"\
		"vfmsub132pd	%%ymm14,%%ymm5,%%ymm1			\n\t		vfmsub132pd	%%ymm14,%%ymm15,%%ymm10			\n\t"\
		"vfmsub132pd	%%ymm14,%%ymm6,%%ymm3			\n\t		vfmsub213pd	(%%rsi),%%ymm14,%%ymm11			\n\t"\
	"vmovaps	(%%r8),%%ymm14 	\n\t"/* two */\
	"vfmadd132pd	%%ymm14,%%ymm0,%%ymm4		\n\t	vfmadd132pd	%%ymm14,%%ymm8 ,%%ymm12	\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm2,%%ymm7		\n\t	vfmadd132pd	%%ymm14,%%ymm9 ,%%ymm13	\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm1,%%ymm5		\n\t	vfmadd132pd	%%ymm14,%%ymm10,%%ymm15	\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm3,%%ymm6		\n\t	vfmadd132pd	(%%rsi),%%ymm11,%%ymm14	\n\t"\
	"vmovaps	%%ymm8 ,(%%rsi) 	\n\t"/* spill ymm8  to make room for 1.0 */"vmovaps	-0x20(%%r8),%%ymm8 	\n\t"/* one */\
		"vfmsub132pd	%%ymm8 ,%%ymm12,%%ymm4			\n\t		vfmsub132pd	%%ymm8 ,%%ymm15,%%ymm11			\n\t"\
		"vfmsub132pd	%%ymm8 ,%%ymm9 ,%%ymm0			\n\t		vfmsub132pd	%%ymm8 ,%%ymm10,%%ymm14			\n\t"\
		"vfmsub132pd	%%ymm8 ,%%ymm13,%%ymm5			\n\t"\
	"vmovaps	(%%r8),%%ymm8  	\n\t"/* two */			"	\n\t	vfmadd132pd	%%ymm8 ,%%ymm11,%%ymm15			\n\t"\
		"vsubpd	(%%rsi),%%ymm1	,%%ymm1			\n\t	vfmadd132pd	%%ymm8 ,%%ymm14,%%ymm10			\n\t"\
	"vfmadd132pd	%%ymm8 ,%%ymm4 ,%%ymm12		\n\t	vfnmadd231pd	-0x900(%%r8),%%ymm15,%%ymm7		\n\t"/* isrt2 */\
	"vfmadd132pd	%%ymm8 ,%%ymm0 ,%%ymm9 		\n\t	vfnmadd231pd	-0x900(%%r8),%%ymm14,%%ymm2		\n\t"\
	"vfmadd132pd	%%ymm8 ,%%ymm5 ,%%ymm13		\n\t	vfnmadd231pd	-0x900(%%r8),%%ymm11,%%ymm3		\n\t"\
	"vfmadd132pd	(%%rsi),%%ymm1 ,%%ymm8 		\n\t	vfnmadd231pd	-0x900(%%r8),%%ymm10,%%ymm6		\n\t"\
		/* SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS(r00,r02,r04,r06,r08,r0A,r0C,r0E) */\
	"vmovaps	%%ymm12,(%%rsi) 	\n\t"/* store ymm12 to make room for sqrt2 */"	vmovaps	0x20(%%r8),%%ymm12 	\n\t"/* sqrt2 */\
		"vmovaps	%%ymm4 ,0x100(%%rsi)		\n\t	vfmadd132pd	%%ymm12,%%ymm7 ,%%ymm15	\n\t"/* sqrt2, use of FMA above means we must combine 2*isrt2 here */\
		"vmovaps	%%ymm0 ,0x180(%%rsi)		\n\t	vfmadd132pd	%%ymm12,%%ymm2 ,%%ymm14	\n\t"\
		"vmovaps	%%ymm5 ,0x120(%%rsi)		\n\t	vfmadd132pd	%%ymm12,%%ymm3 ,%%ymm11	\n\t"\
		"vmovaps	%%ymm1 ,0x0a0(%%rsi)		\n\t	vfmadd132pd	%%ymm12,%%ymm6 ,%%ymm10	\n\t"\
		/*vmovaps	%%ymm12,     (%%rsi)*/	"	\n\t		vmovaps	%%ymm7 ,0x140(%%rsi)	\n\t"\
		"vmovaps	%%ymm9 ,0x080(%%rsi)		\n\t		vmovaps	%%ymm2 ,0x1c0(%%rsi)	\n\t"\
		"vmovaps	%%ymm13,0x020(%%rsi)		\n\t		vmovaps	%%ymm3 ,0x160(%%rsi)	\n\t"\
		"vmovaps	%%ymm8 ,0x1a0(%%rsi)		\n\t		vmovaps	%%ymm6 ,0x0e0(%%rsi)	\n\t"\
		"													vmovaps	%%ymm15,0x040(%%rsi)	\n\t"\
		"													vmovaps	%%ymm14,0x0c0(%%rsi)	\n\t"\
		"													vmovaps	%%ymm11,0x060(%%rsi)	\n\t"\
		"													vmovaps	%%ymm10,0x1e0(%%rsi)	\n\t"\
		"\n\t"\
	/*...Block 2: */\
		"addq	$0x200,%%rsi	/* r10 */		\n\t"\
		"movslq	%[__p08],%%rdi					\n\t"\
		"shlq	$3,%%rdi						\n\t"\
		"addq	%%rdi,%%rax	/* add0+p08 */		\n\t		leaq	(%%r9,%%rax),%%r10	\n\t"\
		"addq	%%rdi,%%rbx						\n\t		leaq	(%%r9,%%rbx),%%r11	\n\t"\
		"addq	%%rdi,%%rcx						\n\t		leaq	(%%r9,%%rcx),%%r12	\n\t"\
		"addq	%%rdi,%%rdx						\n\t		leaq	(%%r9,%%rdx),%%r13	\n\t"\
		/* SSE2_RADIX4_DIT_0TWIDDLE_B(r10) */			/* SSE2_RADIX4_DIT_0TWIDDLE_2NDOFTWO_B(r18) */\
		"vmovaps	     (%%rax),%%ymm2			\n\t		vmovaps	     (%%r10),%%ymm10		\n\t"\
		"vmovaps	     (%%rcx),%%ymm6			\n\t		vmovaps	     (%%r12),%%ymm14		\n\t"\
		"vmovaps	0x020(%%rax),%%ymm3			\n\t		vmovaps	0x020(%%r10),%%ymm11		\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm7			\n\t"	/*	vmovaps	0x020(%%r12),%%ymm15	*/\
		"vmovaps	     (%%rbx),%%ymm0			\n\t		vmovaps	     (%%r11),%%ymm8 		\n\t"\
		"vmovaps	     (%%rdx),%%ymm4			\n\t		vmovaps	     (%%r13),%%ymm12		\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm1			\n\t		vmovaps	0x020(%%r11),%%ymm9 		\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm5			\n\t		vmovaps	0x020(%%r13),%%ymm13		\n\t"\
	"vmovaps	-0x20(%%r8),%%ymm15	\n\t"/* ymm15 free; use for 1.0 */\
		"vfmsub132pd	%%ymm15,%%ymm0,%%ymm2			\n\t		vfmsub132pd	%%ymm15,%%ymm8 ,%%ymm10			\n\t"\
		"vfmsub132pd	%%ymm15,%%ymm4,%%ymm6			\n\t		vfmsub132pd	%%ymm15,%%ymm12,%%ymm14			\n\t"\
		"vfmsub132pd	%%ymm15,%%ymm1,%%ymm3			\n\t		vfmsub132pd	%%ymm15,%%ymm9 ,%%ymm11			\n\t"\
		"vfmsub132pd	%%ymm15,%%ymm5,%%ymm7			\n\t		vfmsub132pd	0x20(%%r12),%%ymm13,%%ymm15		\n\t"\
	"vmovaps	%%ymm13,(%%rsi) 	\n\t"/* spill ymm13 to make room for 2.0 */"	vmovaps	(%%r8),%%ymm13 	\n\t"/* two */\
	"vfmadd132pd	%%ymm13,%%ymm2,%%ymm0		\n\t	vfmadd132pd	%%ymm13,%%ymm10,%%ymm8 	\n\t"\
	"vfmadd132pd	%%ymm13,%%ymm6,%%ymm4		\n\t	vfmadd132pd	%%ymm13,%%ymm14,%%ymm12	\n\t"\
	"vfmadd132pd	%%ymm13,%%ymm3,%%ymm1		\n\t	vfmadd132pd	%%ymm13,%%ymm11,%%ymm9 	\n\t"\
	"vfmadd132pd	%%ymm13,%%ymm7,%%ymm5		\n\t	vfmadd132pd	(%%rsi),%%ymm15,%%ymm13	\n\t"\
		/* Finish radix-4 butterfly: */\
	"vmovaps	%%ymm14,(%%rsi) 	\n\t"/* spill ymm11 to make room for 1.0 */"vmovaps	-0x20(%%r8),%%ymm14 \n\t"/* one */\
		"vfmsub132pd	%%ymm14,%%ymm4,%%ymm0			\n\t		vfmsub132pd	%%ymm14,%%ymm12,%%ymm8 			\n\t"\
		"vfmsub132pd	%%ymm14,%%ymm7,%%ymm2			\n\t		vfmsub132pd	%%ymm14,%%ymm13,%%ymm9 			\n\t"\
		"vfmsub132pd	%%ymm14,%%ymm5,%%ymm1			\n\t		vfmsub132pd	%%ymm14,%%ymm15,%%ymm10			\n\t"\
		"vfmsub132pd	%%ymm14,%%ymm6,%%ymm3			\n\t		vfmsub213pd	(%%rsi),%%ymm14,%%ymm11			\n\t"\
	"vmovaps	(%%r8),%%ymm14 	\n\t"/* two */\
	"vfmadd132pd	%%ymm14,%%ymm0,%%ymm4		\n\t	vfmadd132pd	%%ymm14,%%ymm8 ,%%ymm12	\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm2,%%ymm7		\n\t	vfmadd132pd	%%ymm14,%%ymm9 ,%%ymm13	\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm1,%%ymm5		\n\t	vfmadd132pd	%%ymm14,%%ymm10,%%ymm15	\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm3,%%ymm6		\n\t	vfmadd132pd	(%%rsi),%%ymm11,%%ymm14	\n\t"\
	"vmovaps	%%ymm8 ,(%%rsi) 	\n\t"/* spill ymm8  to make room for 1.0 */"vmovaps	-0x20(%%r8),%%ymm8 	\n\t"/* one */\
		"vfmsub132pd	%%ymm8 ,%%ymm12,%%ymm4			\n\t		vfmsub132pd	%%ymm8 ,%%ymm15,%%ymm11			\n\t"\
		"vfmsub132pd	%%ymm8 ,%%ymm9 ,%%ymm0			\n\t		vfmsub132pd	%%ymm8 ,%%ymm10,%%ymm14			\n\t"\
		"vfmsub132pd	%%ymm8 ,%%ymm13,%%ymm5			\n\t"\
	"vmovaps	(%%r8),%%ymm8  	\n\t"/* two */			"	\n\t	vfmadd132pd	%%ymm8 ,%%ymm11,%%ymm15			\n\t"\
		"vsubpd	(%%rsi),%%ymm1	,%%ymm1			\n\t	vfmadd132pd	%%ymm8 ,%%ymm14,%%ymm10			\n\t"\
	"vfmadd132pd	%%ymm8 ,%%ymm4 ,%%ymm12		\n\t	vfnmadd231pd	-0x900(%%r8),%%ymm15,%%ymm7		\n\t"/* isrt2 */\
	"vfmadd132pd	%%ymm8 ,%%ymm0 ,%%ymm9 		\n\t	vfnmadd231pd	-0x900(%%r8),%%ymm14,%%ymm2		\n\t"\
	"vfmadd132pd	%%ymm8 ,%%ymm5 ,%%ymm13		\n\t	vfnmadd231pd	-0x900(%%r8),%%ymm11,%%ymm3		\n\t"\
	"vfmadd132pd	(%%rsi),%%ymm1 ,%%ymm8 		\n\t	vfnmadd231pd	-0x900(%%r8),%%ymm10,%%ymm6		\n\t"\
		/* SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS(r10,r12,r14,r16,r18,r1A,r1C,r1E) */\
	"vmovaps	%%ymm12,(%%rsi) 	\n\t"/* store ymm12 to make room for sqrt2 */"	vmovaps	0x20(%%r8),%%ymm12 	\n\t"/* sqrt2 */\
		"vmovaps	%%ymm4 ,0x100(%%rsi)		\n\t	vfmadd132pd	%%ymm12,%%ymm7 ,%%ymm15	\n\t"/* sqrt2, use of FMA above means we must combine 2*isrt2 here */\
		"vmovaps	%%ymm0 ,0x180(%%rsi)		\n\t	vfmadd132pd	%%ymm12,%%ymm2 ,%%ymm14	\n\t"\
		"vmovaps	%%ymm5 ,0x120(%%rsi)		\n\t	vfmadd132pd	%%ymm12,%%ymm3 ,%%ymm11	\n\t"\
		"vmovaps	%%ymm1 ,0x0a0(%%rsi)		\n\t	vfmadd132pd	%%ymm12,%%ymm6 ,%%ymm10	\n\t"\
		/*vmovaps	%%ymm12,     (%%rsi)*/	"	\n\t		vmovaps	%%ymm7 ,0x140(%%rsi)	\n\t"\
		"vmovaps	%%ymm9 ,0x080(%%rsi)		\n\t		vmovaps	%%ymm2 ,0x1c0(%%rsi)	\n\t"\
		"vmovaps	%%ymm13,0x020(%%rsi)		\n\t		vmovaps	%%ymm3 ,0x160(%%rsi)	\n\t"\
		"vmovaps	%%ymm8 ,0x1a0(%%rsi)		\n\t		vmovaps	%%ymm6 ,0x0e0(%%rsi)	\n\t"\
		"													vmovaps	%%ymm15,0x040(%%rsi)	\n\t"\
		"													vmovaps	%%ymm14,0x0c0(%%rsi)	\n\t"\
		"													vmovaps	%%ymm11,0x060(%%rsi)	\n\t"\
		"													vmovaps	%%ymm10,0x1e0(%%rsi)	\n\t"\
		"\n\t"\
	/*...Block 3: */\
		"addq	$0x200,%%rsi	/* r20 */		\n\t"\
		"movslq	%[__p10],%%r14					\n\t"\
		"shlq	$3,%%r14						\n\t"\
		"subq	%%rdi,%%r14	/* p10-p8 */		\n\t"\
		"addq	%%r14,%%rax	/* add0+p10 */		\n\t		leaq	(%%r9,%%rax),%%r10	\n\t"\
		"addq	%%r14,%%rbx						\n\t		leaq	(%%r9,%%rbx),%%r11	\n\t"\
		"addq	%%r14,%%rcx						\n\t		leaq	(%%r9,%%rcx),%%r12	\n\t"\
		"addq	%%r14,%%rdx						\n\t		leaq	(%%r9,%%rdx),%%r13	\n\t"\
		/* SSE2_RADIX4_DIT_0TWIDDLE_B(r20) */			/* SSE2_RADIX4_DIT_0TWIDDLE_2NDOFTWO_B(r28) */\
		"vmovaps	     (%%rax),%%ymm2			\n\t		vmovaps	     (%%r10),%%ymm10		\n\t"\
		"vmovaps	     (%%rcx),%%ymm6			\n\t		vmovaps	     (%%r12),%%ymm14		\n\t"\
		"vmovaps	0x020(%%rax),%%ymm3			\n\t		vmovaps	0x020(%%r10),%%ymm11		\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm7			\n\t"	/*	vmovaps	0x020(%%r12),%%ymm15	*/\
		"vmovaps	     (%%rbx),%%ymm0			\n\t		vmovaps	     (%%r11),%%ymm8 		\n\t"\
		"vmovaps	     (%%rdx),%%ymm4			\n\t		vmovaps	     (%%r13),%%ymm12		\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm1			\n\t		vmovaps	0x020(%%r11),%%ymm9 		\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm5			\n\t		vmovaps	0x020(%%r13),%%ymm13		\n\t"\
	"vmovaps	-0x20(%%r8),%%ymm15	\n\t"/* ymm15 free; use for 1.0 */\
		"vfmsub132pd	%%ymm15,%%ymm0,%%ymm2			\n\t		vfmsub132pd	%%ymm15,%%ymm8 ,%%ymm10			\n\t"\
		"vfmsub132pd	%%ymm15,%%ymm4,%%ymm6			\n\t		vfmsub132pd	%%ymm15,%%ymm12,%%ymm14			\n\t"\
		"vfmsub132pd	%%ymm15,%%ymm1,%%ymm3			\n\t		vfmsub132pd	%%ymm15,%%ymm9 ,%%ymm11			\n\t"\
		"vfmsub132pd	%%ymm15,%%ymm5,%%ymm7			\n\t		vfmsub132pd	0x20(%%r12),%%ymm13,%%ymm15		\n\t"\
	"vmovaps	%%ymm13,(%%rsi) 	\n\t"/* spill ymm13 to make room for 2.0 */"	vmovaps	(%%r8),%%ymm13 	\n\t"/* two */\
	"vfmadd132pd	%%ymm13,%%ymm2,%%ymm0		\n\t	vfmadd132pd	%%ymm13,%%ymm10,%%ymm8 	\n\t"\
	"vfmadd132pd	%%ymm13,%%ymm6,%%ymm4		\n\t	vfmadd132pd	%%ymm13,%%ymm14,%%ymm12	\n\t"\
	"vfmadd132pd	%%ymm13,%%ymm3,%%ymm1		\n\t	vfmadd132pd	%%ymm13,%%ymm11,%%ymm9 	\n\t"\
	"vfmadd132pd	%%ymm13,%%ymm7,%%ymm5		\n\t	vfmadd132pd	(%%rsi),%%ymm15,%%ymm13	\n\t"\
		/* Finish radix-4 butterfly: */\
	"vmovaps	%%ymm14,(%%rsi) 	\n\t"/* spill ymm11 to make room for 1.0 */"vmovaps	-0x20(%%r8),%%ymm14 \n\t"/* one */\
		"vfmsub132pd	%%ymm14,%%ymm4,%%ymm0			\n\t		vfmsub132pd	%%ymm14,%%ymm12,%%ymm8 			\n\t"\
		"vfmsub132pd	%%ymm14,%%ymm7,%%ymm2			\n\t		vfmsub132pd	%%ymm14,%%ymm13,%%ymm9 			\n\t"\
		"vfmsub132pd	%%ymm14,%%ymm5,%%ymm1			\n\t		vfmsub132pd	%%ymm14,%%ymm15,%%ymm10			\n\t"\
		"vfmsub132pd	%%ymm14,%%ymm6,%%ymm3			\n\t		vfmsub213pd	(%%rsi),%%ymm14,%%ymm11			\n\t"\
	"vmovaps	(%%r8),%%ymm14 	\n\t"/* two */\
	"vfmadd132pd	%%ymm14,%%ymm0,%%ymm4		\n\t	vfmadd132pd	%%ymm14,%%ymm8 ,%%ymm12	\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm2,%%ymm7		\n\t	vfmadd132pd	%%ymm14,%%ymm9 ,%%ymm13	\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm1,%%ymm5		\n\t	vfmadd132pd	%%ymm14,%%ymm10,%%ymm15	\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm3,%%ymm6		\n\t	vfmadd132pd	(%%rsi),%%ymm11,%%ymm14	\n\t"\
	"vmovaps	%%ymm8 ,(%%rsi) 	\n\t"/* spill ymm8  to make room for 1.0 */"vmovaps	-0x20(%%r8),%%ymm8 	\n\t"/* one */\
		"vfmsub132pd	%%ymm8 ,%%ymm12,%%ymm4			\n\t		vfmsub132pd	%%ymm8 ,%%ymm15,%%ymm11			\n\t"\
		"vfmsub132pd	%%ymm8 ,%%ymm9 ,%%ymm0			\n\t		vfmsub132pd	%%ymm8 ,%%ymm10,%%ymm14			\n\t"\
		"vfmsub132pd	%%ymm8 ,%%ymm13,%%ymm5			\n\t"\
	"vmovaps	(%%r8),%%ymm8  	\n\t"/* two */			"	\n\t	vfmadd132pd	%%ymm8 ,%%ymm11,%%ymm15			\n\t"\
		"vsubpd	(%%rsi),%%ymm1	,%%ymm1			\n\t	vfmadd132pd	%%ymm8 ,%%ymm14,%%ymm10			\n\t"\
	"vfmadd132pd	%%ymm8 ,%%ymm4 ,%%ymm12		\n\t	vfnmadd231pd	-0x900(%%r8),%%ymm15,%%ymm7		\n\t"/* isrt2 */\
	"vfmadd132pd	%%ymm8 ,%%ymm0 ,%%ymm9 		\n\t	vfnmadd231pd	-0x900(%%r8),%%ymm14,%%ymm2		\n\t"\
	"vfmadd132pd	%%ymm8 ,%%ymm5 ,%%ymm13		\n\t	vfnmadd231pd	-0x900(%%r8),%%ymm11,%%ymm3		\n\t"\
	"vfmadd132pd	(%%rsi),%%ymm1 ,%%ymm8 		\n\t	vfnmadd231pd	-0x900(%%r8),%%ymm10,%%ymm6		\n\t"\
		/* SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS(r20,r22,r24,r26,r28,r2A,r2C,r2E) */\
	"vmovaps	%%ymm12,(%%rsi) 	\n\t"/* store ymm12 to make room for sqrt2 */"	vmovaps	0x20(%%r8),%%ymm12 	\n\t"/* sqrt2 */\
		"vmovaps	%%ymm4 ,0x100(%%rsi)		\n\t	vfmadd132pd	%%ymm12,%%ymm7 ,%%ymm15	\n\t"/* sqrt2, use of FMA above means we must combine 2*isrt2 here */\
		"vmovaps	%%ymm0 ,0x180(%%rsi)		\n\t	vfmadd132pd	%%ymm12,%%ymm2 ,%%ymm14	\n\t"\
		"vmovaps	%%ymm5 ,0x120(%%rsi)		\n\t	vfmadd132pd	%%ymm12,%%ymm3 ,%%ymm11	\n\t"\
		"vmovaps	%%ymm1 ,0x0a0(%%rsi)		\n\t	vfmadd132pd	%%ymm12,%%ymm6 ,%%ymm10	\n\t"\
		/*vmovaps	%%ymm12,     (%%rsi)*/	"	\n\t		vmovaps	%%ymm7 ,0x140(%%rsi)	\n\t"\
		"vmovaps	%%ymm9 ,0x080(%%rsi)		\n\t		vmovaps	%%ymm2 ,0x1c0(%%rsi)	\n\t"\
		"vmovaps	%%ymm13,0x020(%%rsi)		\n\t		vmovaps	%%ymm3 ,0x160(%%rsi)	\n\t"\
		"vmovaps	%%ymm8 ,0x1a0(%%rsi)		\n\t		vmovaps	%%ymm6 ,0x0e0(%%rsi)	\n\t"\
		"													vmovaps	%%ymm15,0x040(%%rsi)	\n\t"\
		"													vmovaps	%%ymm14,0x0c0(%%rsi)	\n\t"\
		"													vmovaps	%%ymm11,0x060(%%rsi)	\n\t"\
		"													vmovaps	%%ymm10,0x1e0(%%rsi)	\n\t"\
		"\n\t"\
	/*...Block 4: */\
		"addq	$0x200,%%rsi	/* r30 */		\n\t"\
		"movslq	%[__p08],%%rdi					\n\t"\
		"shlq	$3,%%rdi						\n\t"\
		"addq	%%rdi,%%rax	/* add0+p18 */		\n\t		leaq	(%%r9,%%rax),%%r10	\n\t"\
		"addq	%%rdi,%%rbx						\n\t		leaq	(%%r9,%%rbx),%%r11	\n\t"\
		"addq	%%rdi,%%rcx						\n\t		leaq	(%%r9,%%rcx),%%r12	\n\t"\
		"addq	%%rdi,%%rdx						\n\t		leaq	(%%r9,%%rdx),%%r13	\n\t"\
		/* SSE2_RADIX4_DIT_0TWIDDLE_B(r30) */		/* SSE2_RADIX4_DIT_0TWIDDLE_2NDOFTWO_B(r38) */\
		"vmovaps	     (%%rax),%%ymm2			\n\t		vmovaps	     (%%r10),%%ymm10		\n\t"\
		"vmovaps	     (%%rcx),%%ymm6			\n\t		vmovaps	     (%%r12),%%ymm14		\n\t"\
		"vmovaps	0x020(%%rax),%%ymm3			\n\t		vmovaps	0x020(%%r10),%%ymm11		\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm7			\n\t"	/*	vmovaps	0x020(%%r12),%%ymm15	*/\
		"vmovaps	     (%%rbx),%%ymm0			\n\t		vmovaps	     (%%r11),%%ymm8 		\n\t"\
		"vmovaps	     (%%rdx),%%ymm4			\n\t		vmovaps	     (%%r13),%%ymm12		\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm1			\n\t		vmovaps	0x020(%%r11),%%ymm9 		\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm5			\n\t		vmovaps	0x020(%%r13),%%ymm13		\n\t"\
	"vmovaps	-0x20(%%r8),%%ymm15	\n\t"/* ymm15 free; use for 1.0 */\
		"vfmsub132pd	%%ymm15,%%ymm0,%%ymm2			\n\t		vfmsub132pd	%%ymm15,%%ymm8 ,%%ymm10			\n\t"\
		"vfmsub132pd	%%ymm15,%%ymm4,%%ymm6			\n\t		vfmsub132pd	%%ymm15,%%ymm12,%%ymm14			\n\t"\
		"vfmsub132pd	%%ymm15,%%ymm1,%%ymm3			\n\t		vfmsub132pd	%%ymm15,%%ymm9 ,%%ymm11			\n\t"\
		"vfmsub132pd	%%ymm15,%%ymm5,%%ymm7			\n\t		vfmsub132pd	0x20(%%r12),%%ymm13,%%ymm15		\n\t"\
	"vmovaps	%%ymm13,(%%rsi) 	\n\t"/* spill ymm13 to make room for 2.0 */"	vmovaps	(%%r8),%%ymm13 	\n\t"/* two */\
	"vfmadd132pd	%%ymm13,%%ymm2,%%ymm0		\n\t	vfmadd132pd	%%ymm13,%%ymm10,%%ymm8 	\n\t"\
	"vfmadd132pd	%%ymm13,%%ymm6,%%ymm4		\n\t	vfmadd132pd	%%ymm13,%%ymm14,%%ymm12	\n\t"\
	"vfmadd132pd	%%ymm13,%%ymm3,%%ymm1		\n\t	vfmadd132pd	%%ymm13,%%ymm11,%%ymm9 	\n\t"\
	"vfmadd132pd	%%ymm13,%%ymm7,%%ymm5		\n\t	vfmadd132pd	(%%rsi),%%ymm15,%%ymm13	\n\t"\
		/* Finish radix-4 butterfly: */\
	"vmovaps	%%ymm14,(%%rsi) 	\n\t"/* spill ymm11 to make room for 1.0 */"vmovaps	-0x20(%%r8),%%ymm14 \n\t"/* one */\
		"vfmsub132pd	%%ymm14,%%ymm4,%%ymm0			\n\t		vfmsub132pd	%%ymm14,%%ymm12,%%ymm8 			\n\t"\
		"vfmsub132pd	%%ymm14,%%ymm7,%%ymm2			\n\t		vfmsub132pd	%%ymm14,%%ymm13,%%ymm9 			\n\t"\
		"vfmsub132pd	%%ymm14,%%ymm5,%%ymm1			\n\t		vfmsub132pd	%%ymm14,%%ymm15,%%ymm10			\n\t"\
		"vfmsub132pd	%%ymm14,%%ymm6,%%ymm3			\n\t		vfmsub213pd	(%%rsi),%%ymm14,%%ymm11			\n\t"\
	"vmovaps	(%%r8),%%ymm14 	\n\t"/* two */\
	"vfmadd132pd	%%ymm14,%%ymm0,%%ymm4		\n\t	vfmadd132pd	%%ymm14,%%ymm8 ,%%ymm12	\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm2,%%ymm7		\n\t	vfmadd132pd	%%ymm14,%%ymm9 ,%%ymm13	\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm1,%%ymm5		\n\t	vfmadd132pd	%%ymm14,%%ymm10,%%ymm15	\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm3,%%ymm6		\n\t	vfmadd132pd	(%%rsi),%%ymm11,%%ymm14	\n\t"\
	"vmovaps	%%ymm8 ,(%%rsi) 	\n\t"/* spill ymm8  to make room for 1.0 */"vmovaps	-0x20(%%r8),%%ymm8 	\n\t"/* one */\
		"vfmsub132pd	%%ymm8 ,%%ymm12,%%ymm4			\n\t		vfmsub132pd	%%ymm8 ,%%ymm15,%%ymm11			\n\t"\
		"vfmsub132pd	%%ymm8 ,%%ymm9 ,%%ymm0			\n\t		vfmsub132pd	%%ymm8 ,%%ymm10,%%ymm14			\n\t"\
		"vfmsub132pd	%%ymm8 ,%%ymm13,%%ymm5			\n\t"\
	"vmovaps	(%%r8),%%ymm8  	\n\t"/* two */			"	\n\t	vfmadd132pd	%%ymm8 ,%%ymm11,%%ymm15			\n\t"\
		"vsubpd	(%%rsi),%%ymm1	,%%ymm1			\n\t	vfmadd132pd	%%ymm8 ,%%ymm14,%%ymm10			\n\t"\
	"vfmadd132pd	%%ymm8 ,%%ymm4 ,%%ymm12		\n\t	vfnmadd231pd	-0x900(%%r8),%%ymm15,%%ymm7		\n\t"/* isrt2 */\
	"vfmadd132pd	%%ymm8 ,%%ymm0 ,%%ymm9 		\n\t	vfnmadd231pd	-0x900(%%r8),%%ymm14,%%ymm2		\n\t"\
	"vfmadd132pd	%%ymm8 ,%%ymm5 ,%%ymm13		\n\t	vfnmadd231pd	-0x900(%%r8),%%ymm11,%%ymm3		\n\t"\
	"vfmadd132pd	(%%rsi),%%ymm1 ,%%ymm8 		\n\t	vfnmadd231pd	-0x900(%%r8),%%ymm10,%%ymm6		\n\t"\
		/* SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS(r30,r32,r34,r36,r38,r3A,r3C,r3E) */\
	"vmovaps	%%ymm12,(%%rsi) 	\n\t"/* store ymm12 to make room for sqrt2 */"	vmovaps	0x20(%%r8),%%ymm12 	\n\t"/* sqrt2 */\
		"vmovaps	%%ymm4 ,0x100(%%rsi)		\n\t	vfmadd132pd	%%ymm12,%%ymm7 ,%%ymm15	\n\t"/* sqrt2, use of FMA above means we must combine 2*isrt2 here */\
		"vmovaps	%%ymm0 ,0x180(%%rsi)		\n\t	vfmadd132pd	%%ymm12,%%ymm2 ,%%ymm14	\n\t"\
		"vmovaps	%%ymm5 ,0x120(%%rsi)		\n\t	vfmadd132pd	%%ymm12,%%ymm3 ,%%ymm11	\n\t"\
		"vmovaps	%%ymm1 ,0x0a0(%%rsi)		\n\t	vfmadd132pd	%%ymm12,%%ymm6 ,%%ymm10	\n\t"\
		/*vmovaps	%%ymm12,     (%%rsi)*/	"	\n\t		vmovaps	%%ymm7 ,0x140(%%rsi)	\n\t"\
		"vmovaps	%%ymm9 ,0x080(%%rsi)		\n\t		vmovaps	%%ymm2 ,0x1c0(%%rsi)	\n\t"\
		"vmovaps	%%ymm13,0x020(%%rsi)		\n\t		vmovaps	%%ymm3 ,0x160(%%rsi)	\n\t"\
		"vmovaps	%%ymm8 ,0x1a0(%%rsi)		\n\t		vmovaps	%%ymm6 ,0x0e0(%%rsi)	\n\t"\
		"													vmovaps	%%ymm15,0x040(%%rsi)	\n\t"\
		"													vmovaps	%%ymm14,0x0c0(%%rsi)	\n\t"\
		"													vmovaps	%%ymm11,0x060(%%rsi)	\n\t"\
		"													vmovaps	%%ymm10,0x1e0(%%rsi)	\n\t"\
		"\n\t"\
	/**********************************************************************************/\
	/*...and now do eight radix-4 transforms, including the internal twiddle factors: */\
	/**********************************************************************************/\
		"\n\t"\
		"movslq	%[__p10],%%rdi	\n\t"/* tdi will store copy of p10 throughout */\
	/*...Block 1: t00,t10,t20,t30	*/							/*...Block 5: t08,t18,t28,t38*/\
		"movq	%[__add0],%%rax							\n\t		movslq	%[__p04],%%rsi			\n\t"\
		"movslq	%[__p08],%%rbx							\n\t		movq	%%rax,%%r10				\n\t"\
		"movq	%[__r00],%%rcx							\n\t		shlq	$3,%%rsi				\n\t"\
		"movq	%[__r10],%%rdx							\n\t		addq	%%rsi	,%%r10			\n\t"/* add0 = &a[j1+p4] */\
		"shlq	$3,%%rbx								\n\t		movq	%[__isrt2],%%rsi		\n\t"\
		"shlq	$3,%%rdi								\n\t		leaq	(%%rbx,%%r10),%%r11		\n\t"/* add1 = add0+p12; Need this rbx-read before add0 gets added to rbx at left */\
		"addq	%%rax,%%rbx	\n\t"/* add1 = add0+p8 */"	\n\t		vmovaps	(%%rsi),%%ymm10			\n\t"/* isrt2 */\
	"vmovaps	-0x20(%%r8),%%ymm11	\n\t"/* one */\
		"vmovaps	     (%%rdx),%%ymm2					\n\t		vmovaps	0x500(%%rcx),%%ymm12		\n\t"\
		"vmovaps	0x400(%%rdx),%%ymm4					\n\t		vmovaps	0x520(%%rcx),%%ymm13		\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm3					\n\t		vmovaps	0x500(%%rdx),%%ymm14		\n\t"\
		"vmovaps	0x420(%%rdx),%%ymm5					\n\t		vmovaps	0x520(%%rdx),%%ymm15		\n\t"\
		"vmovaps	     (%%rcx),%%ymm0					\n\t		vmulpd	%%ymm10,%%ymm13,%%ymm13		\n\t"\
		"vmovaps	0x400(%%rcx),%%ymm6					\n\t		vmulpd	%%ymm10,%%ymm14,%%ymm14		\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm1					\n\t		vmovaps	0x100(%%rcx),%%ymm8 		\n\t"\
		"vmovaps	0x420(%%rcx),%%ymm7					\n\t		vmovaps	0x100(%%rdx),%%ymm9 		\n\t"\
		"vfmsub132pd	%%ymm11,%%ymm2,%%ymm0					\n\t	vfnmadd231pd	%%ymm12,%%ymm10,%%ymm13		\n\t"\
		"vfmsub132pd	%%ymm11,%%ymm4,%%ymm6					\n\t	vfnmadd231pd	%%ymm15,%%ymm10,%%ymm14		\n\t"\
		"vfmsub132pd	%%ymm11,%%ymm3,%%ymm1					\n\t	 vfmadd132pd	0x20(%%r8),%%ymm13,%%ymm12	\n\t"/* sqrt2 */\
		"vfmsub132pd	%%ymm11,%%ymm5,%%ymm7					\n\t	 vfmadd132pd	0x20(%%r8),%%ymm14,%%ymm15	\n\t"\
	"vmovaps	(%%r8) ,%%ymm10		\n\t"/* two */				/*	vmovaps	0x120(%%rdx),%%ymm10*//* Use ymm10 for 2.0 now, so delay this load */\
	"vfmadd132pd	%%ymm10,%%ymm0,%%ymm2				\n\t		vfmsub213pd	0x120(%%rdx),%%ymm11,%%ymm8 		\n\t"\
	"vfmadd132pd	%%ymm10,%%ymm6,%%ymm4				\n\t		vfmsub132pd	     %%ymm11,%%ymm14,%%ymm12		\n\t"\
	"vfmadd132pd	%%ymm10,%%ymm1,%%ymm3				\n\t		vfmsub132pd		 %%ymm11,%%ymm15,%%ymm13		\n\t"\
	"vfmadd132pd	%%ymm10,%%ymm7,%%ymm5				\n\t		vfmsub132pd	0x120(%%rcx),%%ymm9 ,%%ymm11		\n\t"/* lcol: spill ymm6 to make room for 2.0 */\
		/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c00): */"\n\t		vmovaps	0x120(%%rdx),%%ymm10			\n\t"\
	"vmovaps %%ymm6,(%%rax) \n\t vmovaps (%%r8),%%ymm6  \n\t	vfmadd132pd	%%ymm6 ,%%ymm12,%%ymm14		\n\t"\
		"vfmsub132pd	-0x20(%%r8),%%ymm4,%%ymm2		\n\t	vfmadd132pd	%%ymm6 ,%%ymm8 ,%%ymm10		\n\t"\
		"vfmsub132pd	-0x20(%%r8),%%ymm7,%%ymm0		\n\t	vfmadd132pd	%%ymm6 ,%%ymm13,%%ymm15		\n\t"\
		"vfmsub132pd	-0x20(%%r8),%%ymm5,%%ymm3		\n\t	vfmadd132pd	%%ymm6 ,%%ymm11,%%ymm9 		\n\t"\
		"vsubpd	(%%rax),%%ymm1,%%ymm1					\n\t"	/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c04 = c00 + 0x100): */\
	"vfmadd132pd	%%ymm6 ,%%ymm2,%%ymm4				\n\t	vmovaps	%%ymm2,(%%rcx)	\n\t	vmovaps	-0x20(%%r8),%%ymm2 	\n\t"/* one */\
	"vfmadd132pd	%%ymm6 ,%%ymm0,%%ymm7				\n\t		vfmsub132pd	%%ymm2 ,%%ymm12,%%ymm10		\n\t"\
	"vfmadd132pd	%%ymm6 ,%%ymm3,%%ymm5				\n\t		vfmsub132pd	%%ymm2 ,%%ymm15,%%ymm8 		\n\t"\
	"vfmadd132pd	(%%rax),%%ymm1,%%ymm6				\n\t		vfmsub132pd	%%ymm2 ,%%ymm13,%%ymm11		\n\t"\
		"addq	$0x0e0,%%rsi			\n\t"/* c00 */	"			vfmsub132pd	%%ymm2 ,%%ymm14,%%ymm9 		\n\t"\
		"vmovaps	(%%r8),%%ymm2		\n\t"/* two */	"		vfmadd132pd	%%ymm2 ,%%ymm10,%%ymm12		\n\t"\
		"vmovaps	%%ymm0,0x020(%%rdx)					\n\t	vfmadd132pd	%%ymm2 ,%%ymm8 ,%%ymm15		\n\t"\
		"vmovaps	%%ymm3,0x020(%%rcx)					\n\t	vfmadd132pd	%%ymm2 ,%%ymm11,%%ymm13		\n\t"\
		"vmovaps	%%ymm6,     (%%rdx)					\n\t	vfmadd132pd	%%ymm2 ,%%ymm9 ,%%ymm14		\n\t"\
		"vmovaps	%%ymm4,%%ymm2		\n\t	vmovaps	%%ymm10,0x100(%%rcx)	\n\t	vmovaps	%%ymm12,%%ymm10	\n\t"\
		"vmovaps	%%ymm7,%%ymm0		\n\t	vmovaps	%%ymm8 ,0x120(%%rdx)	\n\t	vmovaps	%%ymm15,%%ymm8 	\n\t"\
		"vmovaps	%%ymm5,%%ymm3		\n\t	vmovaps	%%ymm11,0x120(%%rcx)	\n\t	vmovaps	%%ymm13,%%ymm11	\n\t"\
		"vmovaps	%%ymm1,%%ymm6		\n\t	vmovaps	%%ymm14,0x100(%%rdx)	\n\t	vmovaps	%%ymm9 ,%%ymm14	\n\t"\
		"vmulpd		     (%%rsi),%%ymm5,%%ymm5			\n\t		vmulpd		0x100(%%rsi),%%ymm13,%%ymm13	\n\t"\
		"vmulpd		0x040(%%rsi),%%ymm1,%%ymm1			\n\t		vmulpd		0x140(%%rsi),%%ymm9 ,%%ymm9 	\n\t"\
		"vmulpd		     (%%rsi),%%ymm4,%%ymm4			\n\t		vmulpd		0x100(%%rsi),%%ymm12,%%ymm12	\n\t"\
		"vmulpd		0x040(%%rsi),%%ymm7,%%ymm7			\n\t		vmulpd		0x140(%%rsi),%%ymm15,%%ymm15	\n\t"\
	"vfnmadd231pd	0x020(%%rsi),%%ymm2,%%ymm5			\n\t	vfnmadd231pd	0x120(%%rsi),%%ymm10,%%ymm13	\n\t"\
	"vfnmadd231pd	0x060(%%rsi),%%ymm0,%%ymm1			\n\t	vfnmadd231pd	0x160(%%rsi),%%ymm8 ,%%ymm9 	\n\t"\
	" vfmadd231pd	0x020(%%rsi),%%ymm3,%%ymm4			\n\t	 vfmadd231pd	0x120(%%rsi),%%ymm11,%%ymm12	\n\t"\
	" vfmadd231pd	0x060(%%rsi),%%ymm6,%%ymm7			\n\t	 vfmadd231pd	0x160(%%rsi),%%ymm14,%%ymm15	\n\t"\
		"vmovaps	%%ymm5,0x020(%%rax)					\n\t		vmovaps	%%ymm13,0x020(%%r10)	\n\t"\
		"vmovaps	%%ymm1,0x020(%%rbx)					\n\t		vmovaps	%%ymm9 ,0x020(%%r11)	\n\t"\
		"vmovaps	%%ymm4,     (%%rax)					\n\t		vmovaps	%%ymm12,     (%%r10)	\n\t"\
		"vmovaps	%%ymm7,     (%%rbx)					\n\t		vmovaps	%%ymm15,     (%%r11)	\n\t"\
		"addq	%%rdi,%%rax	\n\t	addq	%%rdi,%%rbx	\n\t		addq	%%rdi,%%r10			\n\t"\
		"addq	$0x080,%%rsi			\n\t"/* c10 */	"			addq	%%rdi,%%r11			\n\t"\
		"vmovaps	     (%%rcx),%%ymm4					\n\t		vmovaps	0x100(%%rcx),%%ymm12	\n\t"/* c14 = c10 + 0x100: */\
		"vmovaps	0x020(%%rdx),%%ymm0					\n\t		vmovaps	0x120(%%rdx),%%ymm8 	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5					\n\t		vmovaps	0x120(%%rcx),%%ymm13	\n\t"\
		"vmovaps	     (%%rdx),%%ymm6					\n\t		vmovaps	0x100(%%rdx),%%ymm14	\n\t"\
		"vmovaps	%%ymm4,%%ymm2						\n\t		vmovaps	%%ymm12,%%ymm10		\n\t"\
		"vmovaps	%%ymm0,%%ymm1						\n\t		vmovaps	%%ymm8 ,%%ymm9 		\n\t"\
		"vmovaps	%%ymm5,%%ymm3						\n\t		vmovaps	%%ymm13,%%ymm11		\n\t"\
		"vmovaps	%%ymm6,%%ymm7						\n\t		vmovaps	%%ymm14,%%ymm15		\n\t"\
		"vmulpd		     (%%rsi),%%ymm5,%%ymm5			\n\t		vmulpd		0x100(%%rsi),%%ymm13,%%ymm13	\n\t"\
		"vmulpd		0x040(%%rsi),%%ymm6,%%ymm6			\n\t		vmulpd		0x140(%%rsi),%%ymm14,%%ymm14	\n\t"\
		"vmulpd		     (%%rsi),%%ymm4,%%ymm4			\n\t		vmulpd		0x100(%%rsi),%%ymm12,%%ymm12	\n\t"\
		"vmulpd		0x040(%%rsi),%%ymm0,%%ymm0			\n\t		vmulpd		0x140(%%rsi),%%ymm8 ,%%ymm8 	\n\t"\
	"vfnmadd231pd	0x020(%%rsi),%%ymm2,%%ymm5			\n\t	vfnmadd231pd	0x120(%%rsi),%%ymm10,%%ymm13	\n\t"\
	"vfnmadd231pd	0x060(%%rsi),%%ymm1,%%ymm6			\n\t	vfnmadd231pd	0x160(%%rsi),%%ymm9 ,%%ymm14	\n\t"\
	" vfmadd231pd	0x020(%%rsi),%%ymm3,%%ymm4			\n\t	 vfmadd231pd	0x120(%%rsi),%%ymm11,%%ymm12	\n\t"\
	" vfmadd231pd	0x060(%%rsi),%%ymm7,%%ymm0			\n\t	 vfmadd231pd	0x160(%%rsi),%%ymm15,%%ymm8 	\n\t"\
		"vmovaps	%%ymm5,0x020(%%rax)					\n\t		vmovaps	%%ymm13,0x020(%%r10)	\n\t"\
		"vmovaps	%%ymm6,0x020(%%rbx)					\n\t		vmovaps	%%ymm14,0x020(%%r11)	\n\t"\
		"vmovaps	%%ymm4,     (%%rax)					\n\t		vmovaps	%%ymm12,     (%%r10)	\n\t"\
		"vmovaps	%%ymm0,     (%%rbx)					\n\t		vmovaps	%%ymm8 ,     (%%r11)	\n\t"\
		"\n\t"\
	/*...Block 2: t02,t12,t22,t32	*/							/*...Block 6: t0A,t1A,t2A,t3A*/\
		"movq	%[__add0],%%rax							\n\t		movq	%%rax,%%r10			\n\t"\
		"movslq	%[__p01],%%rsi							\n\t		movslq	%[__p05],%%r9		\n\t"\
		"movslq	%[__p08],%%rbx							\n\t		shlq	$3,%%r9			\n\t"\
		"shlq	$3,%%rsi								\n\t		addq	%%r9 	,%%r10	/* add0 = &a[j1+p5] */\n\t"\
		"shlq	$3,%%rbx								\n\t		movq	%%r10,%%r11			\n\t"\
		"addq	%%rsi	,%%rax	/* add0 = &a[j1+p1] */	\n\t		addq	%%rbx	,%%r11	/* add1 = add0+p13 */\n\t"\
		"addq	%%rax	,%%rbx	/* add1 = add0+p8 */	\n\t"\
		"movq	%[__isrt2],%%rsi						\n\t"\
		"addq	$0x440,%%rcx	/* r22 */				\n\t"\
		"addq	$0x440,%%rdx	/* r32 */				\n\t"\
		"addq	$0x060,%%rsi	/* cc1; cc3 += 0x040: */\n\t"\
		"vmovaps	     (%%rcx),%%ymm4					\n\t		vmovaps	0x100(%%rcx),%%ymm12	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5					\n\t		vmovaps	0x120(%%rcx),%%ymm13	\n\t"\
		"vmovaps	     (%%rsi),%%ymm2		\n\t"/* cc1 */"			vmovaps	0x040(%%rsi),%%ymm11	\n\t"/* cc3 */\
		"vmovaps	0x020(%%rsi),%%ymm3					\n\t		vmovaps	0x060(%%rsi),%%ymm10	\n\t"\
		"vmovaps	%%ymm4,%%ymm6						\n\t		vmovaps	%%ymm12,%%ymm14		\n\t"\
		"vmovaps	%%ymm5,%%ymm7						\n\t		vmovaps	%%ymm13,%%ymm15		\n\t"\
		"vmulpd		%%ymm2,%%ymm4,%%ymm4				\n\t		vmulpd		%%ymm10,%%ymm12,%%ymm12		\n\t"\
		"vmulpd		%%ymm2,%%ymm5,%%ymm5				\n\t		vmulpd		%%ymm10,%%ymm13,%%ymm13		\n\t"\
		"vmovaps	     (%%rdx),%%ymm0					\n\t		vmovaps	0x100(%%rdx),%%ymm8 	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm1					\n\t		vmovaps	0x120(%%rdx),%%ymm9 	/* t3B */\n\t"\
	" vfmadd231pd	%%ymm3,%%ymm7,%%ymm4				\n\t	 vfmadd231pd	%%ymm11,%%ymm15,%%ymm12		\n\t"\
	"vfnmadd231pd	%%ymm3,%%ymm6,%%ymm5				\n\t	vfnmadd231pd	%%ymm11,%%ymm14,%%ymm13		\n\t"\
		"vmovaps	%%ymm0,%%ymm6						\n\t		vmovaps	%%ymm8 ,%%ymm14		\n\t"\
		"vmovaps	%%ymm1,%%ymm7						\n\t		vmovaps	%%ymm9 ,%%ymm15		\n\t"\
		"vmulpd		%%ymm11,%%ymm0,%%ymm0				\n\t		vmulpd		%%ymm2 ,%%ymm8 ,%%ymm8 		\n\t"\
		"vmulpd		%%ymm11,%%ymm1,%%ymm1				\n\t		vmulpd		%%ymm2 ,%%ymm9 ,%%ymm9 		\n\t"\
	" vfmadd231pd	%%ymm10,%%ymm7,%%ymm0				\n\t	vfnmadd231pd	%%ymm3 ,%%ymm15,%%ymm8 		\n\t"\
	"vfnmadd231pd	%%ymm10,%%ymm6,%%ymm1				\n\t	 vfmadd231pd	%%ymm3 ,%%ymm14,%%ymm9 		\n\t"\
		"vmovaps	%%ymm5,%%ymm7						\n\t		vmovaps	%%ymm13,%%ymm15		\n\t"\
		"vmovaps	%%ymm4,%%ymm6						\n\t		vmovaps	%%ymm12,%%ymm14		\n\t"\
	"vmovaps	-0x20(%%r8),%%ymm2 	\n\t"/* ymm2 free; use for 1.0 */\
		"vfmadd132pd	%%ymm2 ,%%ymm0,%%ymm4			\n\t		vfmadd132pd	%%ymm2 ,%%ymm8 ,%%ymm14		\n\t"\
		"vfmadd132pd	%%ymm2 ,%%ymm1,%%ymm5			\n\t		vfmadd132pd	%%ymm2 ,%%ymm9 ,%%ymm15		\n\t"\
		"vfmsub132pd	%%ymm2 ,%%ymm0,%%ymm6			\n\t		vfmsub132pd	%%ymm2 ,%%ymm8 ,%%ymm12		\n\t"\
		"vfmsub132pd	%%ymm2 ,%%ymm1,%%ymm7			\n\t		vfmsub132pd	%%ymm2 ,%%ymm9 ,%%ymm13		\n\t"\
		"subq	$0x400,%%rcx	/* r02 */				\n\t"\
		"subq	$0x400,%%rdx	/* r12 */				\n\t"\
		"subq	$0x040,%%rsi	/* cc0 */				\n\t"\
		"vmovaps	     (%%rdx),%%ymm1					\n\t		vmovaps	0x100(%%rdx),%%ymm8 	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm3					\n\t		vmovaps	0x120(%%rdx),%%ymm10	\n\t"\
		"vmovaps	0x020(%%rsi),%%ymm2					\n\t		vmovaps	0x020(%%rsi),%%ymm9 	\n\t"\
		"vmovaps	%%ymm1,%%ymm0						\n\t		vmovaps	%%ymm8 ,%%ymm11		\n\t"\
		"vmulpd		%%ymm2 ,%%ymm1,%%ymm1				\n\t		vmulpd		%%ymm9 ,%%ymm8 ,%%ymm8 		\n\t"\
		"vmulpd		%%ymm3 ,%%ymm2,%%ymm2				\n\t		vmulpd		%%ymm10,%%ymm9 ,%%ymm9 		\n\t"\
	"vmovaps %%ymm9,(%%rax) \n\t vmovaps (%%rsi),%%ymm9	\n\t"/* isrt2 */\
	" vfmsub132pd	%%ymm9 ,%%ymm1,%%ymm3				\n\t	vfnmadd231pd	%%ymm9 ,%%ymm10,%%ymm8 		\n\t"\
	" vfmadd231pd	%%ymm9 ,%%ymm0,%%ymm2				\n\t	 vfmadd213pd	(%%rax),%%ymm11,%%ymm9 		\n\t"\
		"vmovaps	     (%%rcx),%%ymm0					\n\t		vmovaps	0x100(%%rcx),%%ymm10	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm1					\n\t		vmovaps	0x120(%%rcx),%%ymm11	\n\t"\
	"vmovaps	%%ymm9,(%%rax)	\n\t	vmovaps	-0x20(%%r8),%%ymm9	\n\t"/* spill to make room for 1.0 */\
		"vfmsub132pd	%%ymm9 ,%%ymm2,%%ymm0					\n\t		vfmsub132pd	%%ymm9 ,%%ymm8 ,%%ymm10		\n\t"\
		"vfmsub132pd	%%ymm9 ,%%ymm3,%%ymm1					\n\t		vfmsub213pd	(%%rax),%%ymm9 ,%%ymm11		\n\t"\
	"vmovaps	(%%r8),%%ymm9	\n\t"/* 2.0 */\
	"vfmadd132pd	%%ymm9,%%ymm0,%%ymm2				\n\t	vfmadd132pd	%%ymm9 ,%%ymm10,%%ymm8 		\n\t"\
	"vfmadd132pd	%%ymm9,%%ymm1,%%ymm3				\n\t	vfmadd132pd	(%%rax),%%ymm11,%%ymm9 		\n\t"\
		"\n\t"\
		/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c01): */			/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c05): */\
		"addq	$0x4c0,%%rsi	\n\t"/* c01; c05 += 0x100: */\
	"vmovaps	%%ymm14,(%%rax)	\n\t	vmovaps	-0x20(%%r8),%%ymm14	\n\t"/* spill to make room for 1.0 */\
		"vfmsub132pd	%%ymm14,%%ymm4,%%ymm2					\n\t		vfmsub132pd	%%ymm14,%%ymm12,%%ymm10		\n\t"\
		"vfmsub132pd	%%ymm14,%%ymm7,%%ymm0					\n\t		vfmsub132pd	%%ymm14,%%ymm15,%%ymm8 		\n\t"\
		"vfmsub132pd	%%ymm14,%%ymm5,%%ymm3					\n\t		vfmsub132pd	%%ymm14,%%ymm13,%%ymm11		\n\t"\
		"vfmsub132pd	%%ymm14,%%ymm6,%%ymm1					\n\t		vfmsub213pd	(%%rax),%%ymm14,%%ymm9 		\n\t"\
	"vmovaps (%%r8),%%ymm14\n\t"/* 2.0 */\
	"vfmadd132pd	%%ymm14,%%ymm2,%%ymm4				\n\t	vfmadd132pd	%%ymm14,%%ymm10,%%ymm12		\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm0,%%ymm7				\n\t	vfmadd132pd	%%ymm14,%%ymm8 ,%%ymm15		\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm3,%%ymm5				\n\t	vfmadd132pd	%%ymm14,%%ymm11,%%ymm13		\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm1,%%ymm6				\n\t	vfmadd132pd	(%%rax),%%ymm9 ,%%ymm14		\n\t"\
		"vmovaps	%%ymm2,     (%%rcx)					\n\t		vmovaps	%%ymm10,0x100(%%rcx)	\n\t"\
		"vmovaps	%%ymm0,0x020(%%rdx)					\n\t		vmovaps	%%ymm8 ,0x120(%%rdx)	\n\t"\
		"vmovaps	%%ymm3,0x020(%%rcx)					\n\t		vmovaps	%%ymm11,0x120(%%rcx)	\n\t"\
		"vmovaps	%%ymm6,     (%%rdx)					\n\t		vmovaps	%%ymm14,0x100(%%rdx)	\n\t"\
		"vmovaps	%%ymm4,%%ymm2						\n\t		vmovaps	%%ymm12,%%ymm10		\n\t"\
		"vmovaps	%%ymm7,%%ymm0						\n\t		vmovaps	%%ymm15,%%ymm8 		\n\t"\
		"vmovaps	%%ymm5,%%ymm3						\n\t		vmovaps	%%ymm13,%%ymm11		\n\t"\
		"vmovaps	%%ymm1,%%ymm6						\n\t		vmovaps	%%ymm9 ,%%ymm14		\n\t"\
		"vmulpd		     (%%rsi),%%ymm4,%%ymm4			\n\t		vmulpd		0x100(%%rsi),%%ymm12,%%ymm12	\n\t"\
		"vmulpd		0x040(%%rsi),%%ymm7,%%ymm7			\n\t		vmulpd		0x140(%%rsi),%%ymm15,%%ymm15	\n\t"\
		"vmulpd		     (%%rsi),%%ymm5,%%ymm5			\n\t		vmulpd		0x100(%%rsi),%%ymm13,%%ymm13	\n\t"\
		"vmulpd		0x040(%%rsi),%%ymm1,%%ymm1			\n\t		vmulpd		0x140(%%rsi),%%ymm9 ,%%ymm9 	\n\t"\
	" vfmadd231pd	0x020(%%rsi),%%ymm3,%%ymm4			\n\t	 vfmadd231pd	0x120(%%rsi),%%ymm11,%%ymm12	\n\t"\
	" vfmadd231pd	0x060(%%rsi),%%ymm6,%%ymm7			\n\t	 vfmadd231pd	0x160(%%rsi),%%ymm14,%%ymm15	\n\t"\
	"vfnmadd231pd	0x020(%%rsi),%%ymm2,%%ymm5			\n\t	vfnmadd231pd	0x120(%%rsi),%%ymm10,%%ymm13	\n\t"\
	"vfnmadd231pd	0x060(%%rsi),%%ymm0,%%ymm1			\n\t	vfnmadd231pd	0x160(%%rsi),%%ymm8 ,%%ymm9 	\n\t"\
		"vmovaps	%%ymm5,0x020(%%rax)					\n\t		vmovaps	%%ymm13,0x020(%%r10)	\n\t"\
		"vmovaps	%%ymm1,0x020(%%rbx)					\n\t		vmovaps	%%ymm9 ,0x020(%%r11)	\n\t"\
		"vmovaps	%%ymm4,     (%%rax)					\n\t		vmovaps	%%ymm12,     (%%r10)	\n\t"\
		"vmovaps	%%ymm7,     (%%rbx)					\n\t		vmovaps	%%ymm15,     (%%r11)	\n\t"\
		"addq	%%rdi,%%rax								\n\t		addq	%%rdi,%%r10			\n\t"\
		"addq	%%rdi,%%rbx								\n\t		addq	%%rdi,%%r11			\n\t"\
		"addq	$0x080,%%rsi	\n\t"/* c11; c15 += 0x100: */\
		"vmovaps	     (%%rcx),%%ymm4					\n\t		vmovaps	0x100(%%rcx),%%ymm12	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm0					\n\t		vmovaps	0x120(%%rdx),%%ymm8 	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5					\n\t		vmovaps	0x120(%%rcx),%%ymm13	\n\t"\
		"vmovaps	     (%%rdx),%%ymm6					\n\t		vmovaps	0x100(%%rdx),%%ymm14	\n\t"\
		"vmovaps	%%ymm4,%%ymm2						\n\t		vmovaps	%%ymm12,%%ymm10		\n\t"\
		"vmovaps	%%ymm0,%%ymm1						\n\t		vmovaps	%%ymm8 ,%%ymm9 		\n\t"\
		"vmovaps	%%ymm5,%%ymm3						\n\t		vmovaps	%%ymm13,%%ymm11		\n\t"\
		"vmovaps	%%ymm6,%%ymm7						\n\t		vmovaps	%%ymm14,%%ymm15		\n\t"\
		"vmulpd		     (%%rsi),%%ymm4,%%ymm4			\n\t		vmulpd		0x100(%%rsi),%%ymm12,%%ymm12	\n\t"\
		"vmulpd		0x040(%%rsi),%%ymm0,%%ymm0			\n\t		vmulpd		0x140(%%rsi),%%ymm8 ,%%ymm8 	\n\t"\
		"vmulpd		     (%%rsi),%%ymm5,%%ymm5			\n\t		vmulpd		0x100(%%rsi),%%ymm13,%%ymm13	\n\t"\
		"vmulpd		0x040(%%rsi),%%ymm6,%%ymm6			\n\t		vmulpd		0x140(%%rsi),%%ymm14,%%ymm14	\n\t"\
	" vfmadd231pd	0x020(%%rsi),%%ymm3,%%ymm4			\n\t	 vfmadd231pd	0x120(%%rsi),%%ymm11,%%ymm12	\n\t"\
	" vfmadd231pd	0x060(%%rsi),%%ymm7,%%ymm0			\n\t	 vfmadd231pd	0x160(%%rsi),%%ymm15,%%ymm8 	\n\t"\
	"vfnmadd231pd	0x020(%%rsi),%%ymm2,%%ymm5			\n\t	vfnmadd231pd	0x120(%%rsi),%%ymm10,%%ymm13	\n\t"\
	"vfnmadd231pd	0x060(%%rsi),%%ymm1,%%ymm6			\n\t	vfnmadd231pd	0x160(%%rsi),%%ymm9 ,%%ymm14	\n\t"\
		"vmovaps	%%ymm5,0x020(%%rax)					\n\t		vmovaps	%%ymm13,0x020(%%r10)	\n\t"\
		"vmovaps	%%ymm6,0x020(%%rbx)					\n\t		vmovaps	%%ymm14,0x020(%%r11)	\n\t"\
		"vmovaps	%%ymm4,     (%%rax)					\n\t		vmovaps	%%ymm12,     (%%r10)	\n\t"\
		"vmovaps	%%ymm0,     (%%rbx)					\n\t		vmovaps	%%ymm8 ,     (%%r11)	\n\t"\
		"\n\t"\
	/*...Block 3: t04,t14,t24,t34*/								/*...Block 7: t0C,t1C,t2C,t3C*/\
		"movq	%[__add0],%%rax							\n\t		movq	%%rax,%%r10			\n\t"\
		"movslq	%[__p02],%%rsi							\n\t		movslq	%[__p06],%%r9		\n\t"\
		"movslq	%[__p08],%%rbx							\n\t		shlq	$3,%%r9			\n\t"\
		"shlq	$3,%%rsi								\n\t		addq	%%r9 	,%%r10	/* add0 = &a[j1+p6] */\n\t"\
		"shlq	$3,%%rbx								\n\t		movq	%%r10,%%r11			\n\t"\
		"addq	%%rsi	,%%rax	/* add0 = &a[j1+p2] */	\n\t		addq	%%rbx	,%%r11	/* add1 = add0+p14 */\n\t"\
		"addq	%%rax	,%%rbx	/* add1 = add0+p8 */	\n\t"\
		"movq	%[__isrt2],%%rsi						\n\t"\
		"addq	$0x440,%%rcx	/* r24 */				\n\t"\
		"addq	$0x440,%%rdx	/* r34 */				\n\t"\
		"addq	$0x020,%%rsi	/* cc0 */				\n\t"\
		"vmovaps	     (%%rcx),%%ymm4					\n\t		vmovaps	0x100(%%rcx),%%ymm12	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5					\n\t		vmovaps	0x120(%%rcx),%%ymm13	\n\t"\
		"vmovaps	     (%%rsi),%%ymm2					\n\t"\
		"vmovaps	0x020(%%rsi),%%ymm3					\n\t"\
		"vmovaps	%%ymm4,%%ymm6						\n\t		vmovaps	%%ymm12,%%ymm14		\n\t"\
		"vmovaps	%%ymm5,%%ymm7						\n\t		vmovaps	%%ymm13,%%ymm15		\n\t"\
		"vmulpd	%%ymm2,%%ymm4,%%ymm4					\n\t		vmulpd	%%ymm3,%%ymm12,%%ymm12	\n\t"\
		"vmulpd	%%ymm2,%%ymm5,%%ymm5					\n\t		vmulpd	%%ymm3,%%ymm13,%%ymm13	\n\t"\
		"vmovaps	     (%%rdx),%%ymm0					\n\t		vmovaps	0x100(%%rdx),%%ymm8 	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm1					\n\t		vmovaps	0x120(%%rdx),%%ymm9 	\n\t"\
	" vfmadd231pd	%%ymm3,%%ymm7,%%ymm4				\n\t	 vfmadd231pd	%%ymm2,%%ymm15,%%ymm12		\n\t"\
	"vfnmadd231pd	%%ymm3,%%ymm6,%%ymm5				\n\t	vfnmadd231pd	%%ymm2,%%ymm14,%%ymm13		\n\t"\
		"vmovaps	%%ymm0,%%ymm6						\n\t		vmovaps	%%ymm8,%%ymm14		\n\t"\
		"vmovaps	%%ymm1,%%ymm7						\n\t		vmovaps	%%ymm9,%%ymm15		\n\t"\
		"vmulpd	%%ymm3,%%ymm0,%%ymm0					\n\t		vmulpd	%%ymm2 ,%%ymm8 ,%%ymm8 	\n\t"\
		"vmulpd	%%ymm3,%%ymm1,%%ymm1					\n\t		vmulpd	%%ymm2 ,%%ymm9 ,%%ymm9 	\n\t"\
	" vfmadd231pd	%%ymm2,%%ymm7,%%ymm0				\n\t	 vfmadd231pd	%%ymm3,%%ymm15,%%ymm8 		\n\t"\
	"vfnmadd231pd	%%ymm2,%%ymm6,%%ymm1				\n\t	vfnmadd231pd	%%ymm3,%%ymm14,%%ymm9 		\n\t"\
		"vmovaps	%%ymm5,%%ymm7						\n\t		vmovaps	%%ymm13,%%ymm15		\n\t"\
		"vmovaps	%%ymm4,%%ymm6						\n\t		vmovaps	%%ymm12,%%ymm14		\n\t"\
	"vmovaps	-0x20(%%r8),%%ymm2 	\n\t"/* ymm2 free; use for 1.0 */\
		"vfmadd132pd	%%ymm2 ,%%ymm0,%%ymm4			\n\t		vfmadd132pd	%%ymm2 ,%%ymm8 ,%%ymm14		\n\t"\
		"vfmadd132pd	%%ymm2 ,%%ymm1,%%ymm5			\n\t		vfmadd132pd	%%ymm2 ,%%ymm9 ,%%ymm15		\n\t"\
		"vfmsub132pd	%%ymm2 ,%%ymm0,%%ymm6			\n\t		vfmsub132pd	%%ymm2 ,%%ymm8 ,%%ymm12		\n\t"\
		"vfmsub132pd	%%ymm2 ,%%ymm1,%%ymm7			\n\t		vfmsub132pd	%%ymm2 ,%%ymm9 ,%%ymm13		\n\t"\
		"subq	$0x400,%%rcx							\n\t"\
		"subq	$0x400,%%rdx							\n\t"\
		"subq	$0x020,%%rsi					\n\t"/* isrt2 */\
		"vmovaps	     (%%rdx),%%ymm2					\n\t		vmovaps	0x100(%%rdx),%%ymm8 	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm3					\n\t		vmovaps	0x120(%%rdx),%%ymm9 	\n\t"\
		"vmovaps	(%%rsi),%%ymm1				\n\t"/* isrt2 */\
		"vmovaps	%%ymm3,%%ymm0						\n\t		vmovaps	%%ymm8 ,%%ymm10		\n\t"\
		"vfmsub132pd	-0x20(%%r8),%%ymm2,%%ymm3					\n\t		vfmsub132pd	-0x20(%%r8),%%ymm9 ,%%ymm8 	\n\t"\
		"vfmadd132pd	-0x20(%%r8),%%ymm0,%%ymm2					\n\t		vfmadd132pd	-0x20(%%r8),%%ymm10,%%ymm9 	\n\t"\
		"vmulpd	%%ymm1,%%ymm2,%%ymm2					\n\t		vmulpd	%%ymm1 ,%%ymm8 	,%%ymm8 \n\t"\
		"vmulpd	%%ymm1,%%ymm3,%%ymm3					\n\t		vmulpd	%%ymm1 ,%%ymm9 	,%%ymm9 \n\t"\
		"vmovaps	     (%%rcx),%%ymm0					\n\t		vmovaps	0x100(%%rcx),%%ymm10	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm1					\n\t		vmovaps	0x120(%%rcx),%%ymm11	\n\t"\
	"vmovaps	%%ymm9,(%%rax)	\n\t	vmovaps	-0x20(%%r8),%%ymm9	\n\t"/* spill to make room for 1.0 */\
		"vfmsub132pd	%%ymm9 ,%%ymm2,%%ymm0					\n\t		vfmsub132pd	%%ymm9 ,%%ymm8 ,%%ymm10		\n\t"\
		"vfmsub132pd	%%ymm9 ,%%ymm3,%%ymm1					\n\t		vfmsub213pd	(%%rax),%%ymm9 ,%%ymm11		\n\t"\
	"vmovaps	(%%r8),%%ymm9	\n\t"/* 2.0 */\
	"vfmadd132pd	%%ymm9,%%ymm0,%%ymm2				\n\t	vfmadd132pd	%%ymm9 ,%%ymm10,%%ymm8 		\n\t"\
	"vfmadd132pd	%%ymm9,%%ymm1,%%ymm3				\n\t	vfmadd132pd	(%%rax),%%ymm11,%%ymm9 		\n\t"\
		"\n\t"\
		/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c02): */			/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c06): */\
		"addq	$0x2e0,%%rsi	\n\t"/* c02; c06 += 0x100: */\
	"vmovaps	%%ymm14,(%%rax)	\n\t	vmovaps	-0x20(%%r8),%%ymm14	\n\t"/* spill to make room for 1.0 */\
		"vfmsub132pd	%%ymm14,%%ymm4,%%ymm2					\n\t		vfmsub132pd	%%ymm14,%%ymm12,%%ymm10		\n\t"\
		"vfmsub132pd	%%ymm14,%%ymm7,%%ymm0					\n\t		vfmsub132pd	%%ymm14,%%ymm15,%%ymm8 		\n\t"\
		"vfmsub132pd	%%ymm14,%%ymm5,%%ymm3					\n\t		vfmsub132pd	%%ymm14,%%ymm13,%%ymm11		\n\t"\
		"vfmsub132pd	%%ymm14,%%ymm6,%%ymm1					\n\t		vfmsub213pd	(%%rax),%%ymm14,%%ymm9 		\n\t"\
	"vmovaps (%%r8),%%ymm14\n\t"/* 2.0 */\
	"vfmadd132pd	%%ymm14,%%ymm2,%%ymm4				\n\t	vfmadd132pd	%%ymm14,%%ymm10,%%ymm12		\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm0,%%ymm7				\n\t	vfmadd132pd	%%ymm14,%%ymm8 ,%%ymm15		\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm3,%%ymm5				\n\t	vfmadd132pd	%%ymm14,%%ymm11,%%ymm13		\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm1,%%ymm6				\n\t	vfmadd132pd	(%%rax),%%ymm9 ,%%ymm14		\n\t"\
		"vmovaps	%%ymm2,     (%%rcx)					\n\t		vmovaps	%%ymm10,0x100(%%rcx)	\n\t"\
		"vmovaps	%%ymm0,0x020(%%rdx)					\n\t		vmovaps	%%ymm8 ,0x120(%%rdx)	\n\t"\
		"vmovaps	%%ymm3,0x020(%%rcx)					\n\t		vmovaps	%%ymm11,0x120(%%rcx)	\n\t"\
		"vmovaps	%%ymm6,     (%%rdx)					\n\t		vmovaps	%%ymm14,0x100(%%rdx)	\n\t"\
		"vmovaps	%%ymm4,%%ymm2						\n\t		vmovaps	%%ymm12,%%ymm10		\n\t"\
		"vmovaps	%%ymm7,%%ymm0						\n\t		vmovaps	%%ymm15,%%ymm8 		\n\t"\
		"vmovaps	%%ymm5,%%ymm3						\n\t		vmovaps	%%ymm13,%%ymm11		\n\t"\
		"vmovaps	%%ymm1,%%ymm6						\n\t		vmovaps	%%ymm9 ,%%ymm14		\n\t"\
		"vmulpd		     (%%rsi),%%ymm4,%%ymm4			\n\t		vmulpd		0x100(%%rsi),%%ymm12,%%ymm12	\n\t"\
		"vmulpd		0x040(%%rsi),%%ymm7,%%ymm7			\n\t		vmulpd		0x140(%%rsi),%%ymm15,%%ymm15	\n\t"\
		"vmulpd		     (%%rsi),%%ymm5,%%ymm5			\n\t		vmulpd		0x100(%%rsi),%%ymm13,%%ymm13	\n\t"\
		"vmulpd		0x040(%%rsi),%%ymm1,%%ymm1			\n\t		vmulpd		0x140(%%rsi),%%ymm9 ,%%ymm9 	\n\t"\
	" vfmadd231pd	0x020(%%rsi),%%ymm3,%%ymm4			\n\t	 vfmadd231pd	0x120(%%rsi),%%ymm11,%%ymm12	\n\t"\
	" vfmadd231pd	0x060(%%rsi),%%ymm6,%%ymm7			\n\t	 vfmadd231pd	0x160(%%rsi),%%ymm14,%%ymm15	\n\t"\
	"vfnmadd231pd	0x020(%%rsi),%%ymm2,%%ymm5			\n\t	vfnmadd231pd	0x120(%%rsi),%%ymm10,%%ymm13	\n\t"\
	"vfnmadd231pd	0x060(%%rsi),%%ymm0,%%ymm1			\n\t	vfnmadd231pd	0x160(%%rsi),%%ymm8 ,%%ymm9 	\n\t"\
		"vmovaps	%%ymm5,0x020(%%rax)					\n\t		vmovaps	%%ymm13,0x020(%%r10)	\n\t"\
		"vmovaps	%%ymm1,0x020(%%rbx)					\n\t		vmovaps	%%ymm9 ,0x020(%%r11)	\n\t"\
		"vmovaps	%%ymm4,     (%%rax)					\n\t		vmovaps	%%ymm12,     (%%r10)	\n\t"\
		"vmovaps	%%ymm7,     (%%rbx)					\n\t		vmovaps	%%ymm15,     (%%r11)	\n\t"\
		"addq	%%rdi,%%rax								\n\t		addq	%%rdi,%%r10			\n\t"\
		"addq	%%rdi,%%rbx								\n\t		addq	%%rdi,%%r11			\n\t"\
		"addq	$0x080,%%rsi	\n\t"/* c12; c16 += 0x100: */\
		"vmovaps	     (%%rcx),%%ymm4					\n\t		vmovaps	0x100(%%rcx),%%ymm12	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm0					\n\t		vmovaps	0x120(%%rdx),%%ymm8 	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5					\n\t		vmovaps	0x120(%%rcx),%%ymm13	\n\t"\
		"vmovaps	     (%%rdx),%%ymm6					\n\t		vmovaps	0x100(%%rdx),%%ymm14	\n\t"\
		"vmovaps	%%ymm4,%%ymm2						\n\t		vmovaps	%%ymm12,%%ymm10		\n\t"\
		"vmovaps	%%ymm0,%%ymm1						\n\t		vmovaps	%%ymm8 ,%%ymm9 		\n\t"\
		"vmovaps	%%ymm5,%%ymm3						\n\t		vmovaps	%%ymm13,%%ymm11		\n\t"\
		"vmovaps	%%ymm6,%%ymm7						\n\t		vmovaps	%%ymm14,%%ymm15		\n\t"\
		"vmulpd		     (%%rsi),%%ymm4,%%ymm4			\n\t		vmulpd		0x100(%%rsi),%%ymm12,%%ymm12	\n\t"\
		"vmulpd		0x040(%%rsi),%%ymm0,%%ymm0			\n\t		vmulpd		0x140(%%rsi),%%ymm8 ,%%ymm8 	\n\t"\
		"vmulpd		     (%%rsi),%%ymm5,%%ymm5			\n\t		vmulpd		0x100(%%rsi),%%ymm13,%%ymm13	\n\t"\
		"vmulpd		0x040(%%rsi),%%ymm6,%%ymm6			\n\t		vmulpd		0x140(%%rsi),%%ymm14,%%ymm14	\n\t"\
	" vfmadd231pd	0x020(%%rsi),%%ymm3,%%ymm4			\n\t	 vfmadd231pd	0x120(%%rsi),%%ymm11,%%ymm12	\n\t"\
	" vfmadd231pd	0x060(%%rsi),%%ymm7,%%ymm0			\n\t	 vfmadd231pd	0x160(%%rsi),%%ymm15,%%ymm8 	\n\t"\
	"vfnmadd231pd	0x020(%%rsi),%%ymm2,%%ymm5			\n\t	vfnmadd231pd	0x120(%%rsi),%%ymm10,%%ymm13	\n\t"\
	"vfnmadd231pd	0x060(%%rsi),%%ymm1,%%ymm6			\n\t	vfnmadd231pd	0x160(%%rsi),%%ymm9 ,%%ymm14	\n\t"\
		"vmovaps	%%ymm5,0x020(%%rax)					\n\t		vmovaps	%%ymm13,0x020(%%r10)	\n\t"\
		"vmovaps	%%ymm6,0x020(%%rbx)					\n\t		vmovaps	%%ymm14,0x020(%%r11)	\n\t"\
		"vmovaps	%%ymm4,     (%%rax)					\n\t		vmovaps	%%ymm12,     (%%r10)	\n\t"\
		"vmovaps	%%ymm0,     (%%rbx)					\n\t		vmovaps	%%ymm8 ,     (%%r11)	\n\t"\
		"\n\t"\
	/*...Block 4: t06,t16,t26,t36*/								/*...Block 8: t0E,t1E,t2E,t3E*/\
		"movq	%[__add0],%%rax							\n\t		movq	%%rax,%%r10			\n\t"\
		"movslq	%[__p03],%%rsi							\n\t		movslq	%[__p07],%%r9		\n\t"\
		"movslq	%[__p08],%%rbx							\n\t		shlq	$3,%%r9			\n\t"\
		"shlq	$3,%%rsi								\n\t		addq	%%r9 	,%%r10	/* add0 = &a[j1+p6] */\n\t"\
		"shlq	$3,%%rbx								\n\t		movq	%%r10,%%r11			\n\t"\
		"addq	%%rsi	,%%rax	/* add0 = &a[j1+p3] */	\n\t		addq	%%rbx	,%%r11	/* add1 = add0+p14 */\n\t"\
		"addq	%%rax	,%%rbx	/* add1 = add0+p8 */	\n\t"\
		"movq	%[__isrt2],%%rsi						\n\t"\
		"addq	$0x440,%%rcx	/* r26 */				\n\t"\
		"addq	$0x440,%%rdx	/* r36 */				\n\t"\
		"addq	$0x060,%%rsi	/* cc1; cc3 += 0x040: */\n\t"\
		"vmovaps	     (%%rcx),%%ymm4					\n\t		vmovaps	0x100(%%rcx),%%ymm12	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5					\n\t		vmovaps	0x120(%%rcx),%%ymm13	\n\t"\
		"vmovaps	     (%%rsi),%%ymm2		\n\t"/* cc1 */"			vmovaps	0x040(%%rsi),%%ymm10	\n\t"/* cc3 */\
		"vmovaps	0x020(%%rsi),%%ymm3					\n\t		vmovaps	0x060(%%rsi),%%ymm11	\n\t"\
		"vmovaps	%%ymm4,%%ymm6						\n\t		vmovaps	%%ymm12,%%ymm14		\n\t"\
		"vmovaps	%%ymm5,%%ymm7						\n\t		vmovaps	%%ymm13,%%ymm15		\n\t"\
		"vmulpd		%%ymm10,%%ymm4,%%ymm4				\n\t		vmulpd		%%ymm3,%%ymm12,%%ymm12		\n\t"\
		"vmulpd		%%ymm10,%%ymm5,%%ymm5				\n\t		vmulpd		%%ymm3,%%ymm13,%%ymm13		\n\t"\
		"vmovaps	     (%%rdx),%%ymm0					\n\t		vmovaps	0x100(%%rdx),%%ymm8 	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm1					\n\t		vmovaps	0x120(%%rdx),%%ymm9 	\n\t"\
	" vfmadd231pd	%%ymm11,%%ymm7,%%ymm4				\n\t	 vfmadd231pd	%%ymm2,%%ymm15,%%ymm12		\n\t"\
	"vfnmadd231pd	%%ymm11,%%ymm6,%%ymm5				\n\t	vfnmadd231pd	%%ymm2,%%ymm14,%%ymm13		\n\t"\
		"vmovaps	%%ymm0,%%ymm6						\n\t		vmovaps	%%ymm8 ,%%ymm14		\n\t"\
		"vmovaps	%%ymm1,%%ymm7						\n\t		vmovaps	%%ymm9 ,%%ymm15		\n\t"\
		"vmulpd		%%ymm3,%%ymm0,%%ymm0				\n\t		vmulpd		%%ymm11,%%ymm8 ,%%ymm8 		\n\t"\
		"vmulpd		%%ymm3,%%ymm1,%%ymm1				\n\t		vmulpd		%%ymm11,%%ymm9 ,%%ymm9 		\n\t"\
	"vfnmadd231pd	%%ymm2,%%ymm7,%%ymm0				\n\t	 vfmadd231pd	%%ymm10,%%ymm15,%%ymm8 		\n\t"\
	" vfmadd231pd	%%ymm2,%%ymm6,%%ymm1				\n\t	vfnmadd231pd	%%ymm10,%%ymm14,%%ymm9 		\n\t"\
		"vmovaps	%%ymm5,%%ymm7						\n\t		vmovaps	%%ymm13,%%ymm15		\n\t"\
		"vmovaps	%%ymm4,%%ymm6						\n\t		vmovaps	%%ymm12,%%ymm14		\n\t"\
	"vmovaps	-0x20(%%r8),%%ymm2 	\n\t"/* ymm2 free; use for 1.0 */\
		"vfmadd132pd	%%ymm2 ,%%ymm0,%%ymm6			\n\t		vfmsub132pd	%%ymm2 ,%%ymm8 ,%%ymm12		\n\t"\
		"vfmadd132pd	%%ymm2 ,%%ymm1,%%ymm7			\n\t		vfmsub132pd	%%ymm2 ,%%ymm9 ,%%ymm13		\n\t"\
		"vfmsub132pd	%%ymm2 ,%%ymm0,%%ymm4			\n\t		vfmadd132pd	%%ymm2 ,%%ymm8 ,%%ymm14		\n\t"\
		"vfmsub132pd	%%ymm2 ,%%ymm1,%%ymm5			\n\t		vfmadd132pd	%%ymm2 ,%%ymm9 ,%%ymm15		\n\t"\
		"subq	$0x400,%%rcx			\n\t"\
		"subq	$0x400,%%rdx			\n\t"\
		"subq	$0x040,%%rsi	/* cc0 */	\n\t"\
		"vmovaps	     (%%rdx),%%ymm2					\n\t		vmovaps	0x100(%%rdx),%%ymm11	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm0					\n\t		vmovaps	0x120(%%rdx),%%ymm9 	\n\t"\
		"vmovaps	0x020(%%rsi),%%ymm3					\n\t		vmovaps	0x020(%%rsi),%%ymm10	\n\t"\
		"vmovaps	%%ymm2,%%ymm1						\n\t		vmovaps	%%ymm11,%%ymm8 		\n\t"\
		"vmulpd	%%ymm3 ,%%ymm2,%%ymm2					\n\t		vmulpd	%%ymm10,%%ymm11,%%ymm11		\n\t"\
		"vmulpd	%%ymm0 ,%%ymm3,%%ymm3					\n\t		vmulpd	%%ymm9 ,%%ymm10,%%ymm10		\n\t"\
	"vmovaps %%ymm8,(%%rax) \n\t vmovaps (%%rsi),%%ymm8	\n\t"/* isrt2 */\
	" vfmadd231pd	%%ymm8 ,%%ymm0,%%ymm2				\n\t	 vfmadd132pd	%%ymm8 ,%%ymm11,%%ymm9 		\n\t"\
	"vfnmadd231pd	%%ymm8 ,%%ymm1,%%ymm3				\n\t	 vfmsub132pd	(%%rax),%%ymm10,%%ymm8 		\n\t"\
		"vmovaps	     (%%rcx),%%ymm0					\n\t		vmovaps	0x100(%%rcx),%%ymm10	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm1					\n\t		vmovaps	0x120(%%rcx),%%ymm11	\n\t"\
	"vmovaps	%%ymm9,(%%rax)	\n\t	vmovaps	-0x20(%%r8),%%ymm9	\n\t"/* spill to make room for 1.0 */\
		"vfmsub132pd	%%ymm9 ,%%ymm2,%%ymm0					\n\t		vfmsub132pd	%%ymm9 ,%%ymm8 ,%%ymm10		\n\t"\
		"vfmsub132pd	%%ymm9 ,%%ymm3,%%ymm1					\n\t		vfmsub213pd	(%%rax),%%ymm9 ,%%ymm11		\n\t"\
	"vmovaps	(%%r8),%%ymm9	\n\t"/* 2.0 */\
	"vfmadd132pd	%%ymm9,%%ymm0,%%ymm2				\n\t	vfmadd132pd	%%ymm9 ,%%ymm10,%%ymm8 		\n\t"\
	"vfmadd132pd	%%ymm9,%%ymm1,%%ymm3				\n\t	vfmadd132pd	(%%rax),%%ymm11,%%ymm9 		\n\t"\
		"\n\t"\
		/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c03): */			/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c07): */\
		"addq	$0x6c0,%%rsi	\n\t"/* c03; c07 += 0x100: */\
	"vmovaps	%%ymm14,(%%rax)	\n\t	vmovaps	-0x20(%%r8),%%ymm14	\n\t"/* spill to make room for 1.0 */\
		"vfmsub132pd	%%ymm14,%%ymm4,%%ymm2					\n\t		vfmsub132pd	%%ymm14,%%ymm12,%%ymm10		\n\t"\
		"vfmsub132pd	%%ymm14,%%ymm7,%%ymm0					\n\t		vfmsub132pd	%%ymm14,%%ymm15,%%ymm8 		\n\t"\
		"vfmsub132pd	%%ymm14,%%ymm5,%%ymm3					\n\t		vfmsub132pd	%%ymm14,%%ymm13,%%ymm11		\n\t"\
		"vfmsub132pd	%%ymm14,%%ymm6,%%ymm1					\n\t		vfmsub213pd	(%%rax),%%ymm14,%%ymm9 		\n\t"\
	"vmovaps (%%r8),%%ymm14\n\t"/* 2.0 */\
	"vfmadd132pd	%%ymm14,%%ymm2,%%ymm4				\n\t	vfmadd132pd	%%ymm14,%%ymm10,%%ymm12		\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm0,%%ymm7				\n\t	vfmadd132pd	%%ymm14,%%ymm8 ,%%ymm15		\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm3,%%ymm5				\n\t	vfmadd132pd	%%ymm14,%%ymm11,%%ymm13		\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm1,%%ymm6				\n\t	vfmadd132pd	(%%rax),%%ymm9 ,%%ymm14		\n\t"\
		"vmovaps	%%ymm2,     (%%rcx)					\n\t		vmovaps	%%ymm10,0x100(%%rcx)	\n\t"\
		"vmovaps	%%ymm0,0x020(%%rdx)					\n\t		vmovaps	%%ymm8 ,0x120(%%rdx)	\n\t"\
		"vmovaps	%%ymm3,0x020(%%rcx)					\n\t		vmovaps	%%ymm11,0x120(%%rcx)	\n\t"\
		"vmovaps	%%ymm6,     (%%rdx)					\n\t		vmovaps	%%ymm14,0x100(%%rdx)	\n\t"\
		"vmovaps	%%ymm4,%%ymm2						\n\t		vmovaps	%%ymm12,%%ymm10		\n\t"\
		"vmovaps	%%ymm7,%%ymm0						\n\t		vmovaps	%%ymm15,%%ymm8 		\n\t"\
		"vmovaps	%%ymm5,%%ymm3						\n\t		vmovaps	%%ymm13,%%ymm11		\n\t"\
		"vmovaps	%%ymm1,%%ymm6						\n\t		vmovaps	%%ymm9 ,%%ymm14		\n\t"\
		"vmulpd		     (%%rsi),%%ymm4,%%ymm4			\n\t		vmulpd		0x100(%%rsi),%%ymm12,%%ymm12	\n\t"\
		"vmulpd		0x040(%%rsi),%%ymm7,%%ymm7			\n\t		vmulpd		0x140(%%rsi),%%ymm15,%%ymm15	\n\t"\
		"vmulpd		     (%%rsi),%%ymm5,%%ymm5			\n\t		vmulpd		0x100(%%rsi),%%ymm13,%%ymm13	\n\t"\
		"vmulpd		0x040(%%rsi),%%ymm1,%%ymm1			\n\t		vmulpd		0x140(%%rsi),%%ymm9 ,%%ymm9 	\n\t"\
	" vfmadd231pd	0x020(%%rsi),%%ymm3,%%ymm4			\n\t	 vfmadd231pd	0x120(%%rsi),%%ymm11,%%ymm12	\n\t"\
	" vfmadd231pd	0x060(%%rsi),%%ymm6,%%ymm7			\n\t	 vfmadd231pd	0x160(%%rsi),%%ymm14,%%ymm15	\n\t"\
	"vfnmadd231pd	0x020(%%rsi),%%ymm2,%%ymm5			\n\t	vfnmadd231pd	0x120(%%rsi),%%ymm10,%%ymm13	\n\t"\
	"vfnmadd231pd	0x060(%%rsi),%%ymm0,%%ymm1			\n\t	vfnmadd231pd	0x160(%%rsi),%%ymm8 ,%%ymm9 	\n\t"\
		"vmovaps	%%ymm5,0x020(%%rax)					\n\t		vmovaps	%%ymm13,0x020(%%r10)	\n\t"\
		"vmovaps	%%ymm1,0x020(%%rbx)					\n\t		vmovaps	%%ymm9 ,0x020(%%r11)	\n\t"\
		"vmovaps	%%ymm4,     (%%rax)					\n\t		vmovaps	%%ymm12,     (%%r10)	\n\t"\
		"vmovaps	%%ymm7,     (%%rbx)					\n\t		vmovaps	%%ymm15,     (%%r11)	\n\t"\
		"addq	%%rdi,%%rax								\n\t		addq	%%rdi,%%r10			\n\t"\
		"addq	%%rdi,%%rbx								\n\t		addq	%%rdi,%%r11			\n\t"\
		"addq	$0x080,%%rsi	\n\t"/* c13; c17 += 0x100: */\
		"vmovaps	     (%%rcx),%%ymm4					\n\t		vmovaps	0x100(%%rcx),%%ymm12	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm0					\n\t		vmovaps	0x120(%%rdx),%%ymm8 	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5					\n\t		vmovaps	0x120(%%rcx),%%ymm13	\n\t"\
		"vmovaps	     (%%rdx),%%ymm6					\n\t		vmovaps	0x100(%%rdx),%%ymm14	\n\t"\
		"vmovaps	%%ymm4,%%ymm2						\n\t		vmovaps	%%ymm12,%%ymm10		\n\t"\
		"vmovaps	%%ymm0,%%ymm1						\n\t		vmovaps	%%ymm8 ,%%ymm9 		\n\t"\
		"vmovaps	%%ymm5,%%ymm3						\n\t		vmovaps	%%ymm13,%%ymm11		\n\t"\
		"vmovaps	%%ymm6,%%ymm7						\n\t		vmovaps	%%ymm14,%%ymm15		\n\t"\
		"vmulpd		     (%%rsi),%%ymm4,%%ymm4			\n\t		vmulpd		0x100(%%rsi),%%ymm12,%%ymm12	\n\t"\
		"vmulpd		0x040(%%rsi),%%ymm0,%%ymm0			\n\t		vmulpd		0x140(%%rsi),%%ymm8 ,%%ymm8 	\n\t"\
		"vmulpd		     (%%rsi),%%ymm5,%%ymm5			\n\t		vmulpd		0x100(%%rsi),%%ymm13,%%ymm13	\n\t"\
		"vmulpd		0x040(%%rsi),%%ymm6,%%ymm6			\n\t		vmulpd		0x140(%%rsi),%%ymm14,%%ymm14	\n\t"\
	" vfmadd231pd	0x020(%%rsi),%%ymm3,%%ymm4			\n\t	 vfmadd231pd	0x120(%%rsi),%%ymm11,%%ymm12	\n\t"\
	" vfmadd231pd	0x060(%%rsi),%%ymm7,%%ymm0			\n\t	 vfmadd231pd	0x160(%%rsi),%%ymm15,%%ymm8 	\n\t"\
	"vfnmadd231pd	0x020(%%rsi),%%ymm2,%%ymm5			\n\t	vfnmadd231pd	0x120(%%rsi),%%ymm10,%%ymm13	\n\t"\
	"vfnmadd231pd	0x060(%%rsi),%%ymm1,%%ymm6			\n\t	vfnmadd231pd	0x160(%%rsi),%%ymm9 ,%%ymm14	\n\t"\
		"vmovaps	%%ymm5,0x020(%%rax)					\n\t		vmovaps	%%ymm13,0x020(%%r10)	\n\t"\
		"vmovaps	%%ymm6,0x020(%%rbx)					\n\t		vmovaps	%%ymm14,0x020(%%r11)	\n\t"\
		"vmovaps	%%ymm4,     (%%rax)					\n\t		vmovaps	%%ymm12,     (%%r10)	\n\t"\
		"vmovaps	%%ymm0,     (%%rbx)					\n\t		vmovaps	%%ymm8 ,     (%%r11)	\n\t"\
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
		 ,[__r10] "m" (Xr10)\
		 ,[__r20] "m" (Xr20)\
		 ,[__r30] "m" (Xr30)\
		 ,[__isrt2] "m" (Xisrt2)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r8","r9","r10","r11","r12","r13","r14","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
	);\
	}

  #else	// ALL_FMA = False: All FMAs have non-unity multiplicands:

	// In our FMA-ized version, original mix of [223 ADD, 223 SUB, 326 MUL] ==> [58+156 ADD, 224 FMA, 102 MUL].
	//
	#define SSE2_RADIX32_DIF_TWIDDLE(Xadd0,Xp01,Xp02,Xp03,Xp04,Xp08,Xp0C,Xp10,Xp18,Xr00)\
	{\
	__asm__ volatile (\
	/*...Block 1: */\
		"movq	%[__add0],%%rax							\n\t"\
		"movslq	%[__p08],%%rbx							\n\t"\
		"movslq	%[__p10],%%rcx							\n\t	movslq	%[__p04],%%r9	\n\t"\
		"movslq	%[__p18],%%rdx							\n\t	movq	%[__r00],%%rsi	\n\t"\
		"shlq	$3,%%rbx								\n\t	shlq	$3,%%r9		\n\t"\
		"shlq	$3,%%rcx								\n\t	leaq	0x1100(%%rsi),%%r8	\n\t"/* two */\
		"shlq	$3,%%rdx								\n\t	leaq	(%%r9,%%rax),%%r10	\n\t"\
		"addq	%%rax,%%rbx								\n\t	leaq	(%%r9,%%rbx),%%r11	\n\t"\
		"addq	%%rax,%%rcx								\n\t	leaq	(%%r9,%%rcx),%%r12	\n\t"\
		"addq	%%rax,%%rdx								\n\t	leaq	(%%r9,%%rdx),%%r13	\n\t"\
		/* SSE2_RADIX4_DIF_4TWIDDLE_A(r00,c10) */				/* SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO_B(r08,c04) */\
		"vmovaps	     (%%rax),%%ymm0					\n\t	vmovaps	     (%%r10),%%ymm8 	\n\t"\
		"vmovaps	     (%%rcx),%%ymm4					\n\t	vmovaps	     (%%r12),%%ymm12	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1					\n\t	vmovaps	0x020(%%r10),%%ymm9 	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5					\n\t	vmovaps	0x020(%%r12),%%ymm13	\n\t"\
		"vmovaps	0x920(%%rsi),%%ymm2	/* c10 */		\n\t	vmovaps	0x9e0(%%rsi),%%ymm14	/* c04 */	\n\t"\
		"vmovaps	0x940(%%rsi),%%ymm3					\n\t	vmovaps	0xa00(%%rsi),%%ymm15	\n\t"\
		"vmovaps	%%ymm4,%%ymm6						\n\t	vmovaps	%%ymm8 ,%%ymm10	\n\t"\
		"vmovaps	%%ymm5,%%ymm7						\n\t	vmovaps	%%ymm9 ,%%ymm11	\n\t"\
		"vmulpd	%%ymm2,%%ymm4,%%ymm4					\n\t	vmulpd	%%ymm14,%%ymm8 ,%%ymm8 	\n\t"\
		"vmulpd	%%ymm2,%%ymm5,%%ymm5					\n\t	vmulpd	%%ymm14,%%ymm9 ,%%ymm9 	\n\t"\
	"vfnmadd231pd	%%ymm3,%%ymm7,%%ymm4			\n\t	vfnmadd231pd	%%ymm15,%%ymm11,%%ymm8 	\n\t"\
	" vfmadd231pd	%%ymm3,%%ymm6,%%ymm5			\n\t	 vfmadd231pd	%%ymm15,%%ymm10,%%ymm9 	\n\t"\
		"vmovaps	%%ymm0,%%ymm2						\n\t	vmovaps	%%ymm12,%%ymm14	\n\t"\
		"vmovaps	%%ymm1,%%ymm3						\n\t	vmovaps	%%ymm13,%%ymm15	\n\t"\
		"vaddpd	%%ymm4,%%ymm0,%%ymm0					\n\t	vmulpd	0xa20(%%rsi),%%ymm12,%%ymm12	/* c14 */	\n\t"\
		"vaddpd	%%ymm5,%%ymm1,%%ymm1					\n\t	vmulpd	0xa20(%%rsi),%%ymm13,%%ymm13	\n\t"\
		"vsubpd	%%ymm4,%%ymm2,%%ymm2				\n\t	vfnmadd231pd	0xa40(%%rsi),%%ymm15,%%ymm12	\n\t"\
		"vsubpd	%%ymm5,%%ymm3,%%ymm3				\n\t	 vfmadd231pd	0xa40(%%rsi),%%ymm14,%%ymm13	\n\t"\
		"vmovaps	     (%%rdx),%%ymm4					\n\t	vmovaps	%%ymm8 ,%%ymm10	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm5					\n\t	vmovaps	%%ymm9 ,%%ymm11	\n\t"\
		"vmovaps	     (%%rdx),%%ymm6					\n\t	vaddpd	%%ymm12,%%ymm8 ,%%ymm8 	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm7					\n\t	vaddpd	%%ymm13,%%ymm9 ,%%ymm9 	\n\t"\
		"vmulpd	0x9a0(%%rsi),%%ymm4,%%ymm4	/* c18 */	\n\t	vsubpd	%%ymm12,%%ymm10,%%ymm10	\n\t"\
		"vmulpd	0x9a0(%%rsi),%%ymm5,%%ymm5				\n\t	vsubpd	%%ymm13,%%ymm11,%%ymm11	\n\t"\
	"vfnmadd231pd	0x9c0(%%rsi),%%ymm7,%%ymm4			\n\t	vmovaps	     (%%r13),%%ymm12	\n\t"\
	" vfmadd231pd	0x9c0(%%rsi),%%ymm6,%%ymm5			\n\t	vmovaps	0x020(%%r13),%%ymm13	\n\t"\
		"vmovaps	%%ymm5,0x020(%%rsi)					\n\t	vmovaps	     (%%r13),%%ymm14	\n\t"\
		"vmovaps	%%ymm4,     (%%rsi)	/* tmpstr r00 */\n\t	vmovaps	0x020(%%r13),%%ymm15	\n\t"\
		"vmovaps	     (%%rbx),%%ymm4					\n\t	vmulpd	0xaa0(%%rsi),%%ymm12,%%ymm12	/* c1C */	\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm5					\n\t	vmulpd	0xaa0(%%rsi),%%ymm13,%%ymm13	\n\t"\
		"vmovaps	     (%%rbx),%%ymm6				\n\t	vfnmadd231pd	0xac0(%%rsi),%%ymm15,%%ymm12	\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm7				\n\t	 vfmadd231pd	0xac0(%%rsi),%%ymm14,%%ymm13	\n\t"\
		"vmulpd	0x960(%%rsi),%%ymm4,%%ymm4	/* c08 */	\n\t	vmovaps	%%ymm13,0x120(%%rsi)	\n\t"\
		"vmulpd	0x960(%%rsi),%%ymm5,%%ymm5				\n\t	vmovaps	%%ymm12,0x100(%%rsi)	/* tmpstr r08 */\n\t"\
	"vfnmadd231pd	0x980(%%rsi),%%ymm7,%%ymm4			\n\t	vmovaps	     (%%r11),%%ymm12	\n\t"\
	" vfmadd231pd	0x980(%%rsi),%%ymm6,%%ymm5			\n\t	vmovaps	0x020(%%r11),%%ymm13	\n\t"\
		"vmovaps	%%ymm5,%%ymm7						\n\t	vmovaps	     (%%r11),%%ymm14	\n\t"\
		"vmovaps	%%ymm4,%%ymm6						\n\t	vmovaps	0x020(%%r11),%%ymm15	\n\t"\
		"vsubpd	     (%%rsi),%%ymm4,%%ymm4	/* r00 */	\n\t	vmulpd	0xa60(%%rsi),%%ymm12,%%ymm12	/* c0C */	\n\t"\
		"vsubpd	0x020(%%rsi),%%ymm5,%%ymm5				\n\t	vmulpd	0xa60(%%rsi),%%ymm13,%%ymm13	\n\t"\
		"vaddpd	     (%%rsi),%%ymm6,%%ymm6			\n\t	vfnmadd231pd	0xa80(%%rsi),%%ymm15,%%ymm12	\n\t"\
		"vaddpd	0x020(%%rsi),%%ymm7,%%ymm7			\n\t	 vfmadd231pd	0xa80(%%rsi),%%ymm14,%%ymm13	\n\t"\
		"vsubpd	%%ymm6,%%ymm0,%%ymm0					\n\t	vmovaps	%%ymm13,%%ymm15	\n\t"\
		"vsubpd	%%ymm5,%%ymm2,%%ymm2					\n\t	vmovaps	%%ymm12,%%ymm14	\n\t"\
		"vsubpd	%%ymm7,%%ymm1,%%ymm1					\n\t	vsubpd	0x100(%%rsi),%%ymm12,%%ymm12	\n\t"\
		"vsubpd	%%ymm4,%%ymm3,%%ymm3					\n\t	vsubpd	0x120(%%rsi),%%ymm13,%%ymm13	\n\t"\
		"												\n\t	vaddpd	0x100(%%rsi),%%ymm14,%%ymm14	\n\t"\
		"												\n\t	vaddpd	0x120(%%rsi),%%ymm15,%%ymm15	\n\t"\
		"												\n\t	vsubpd	     %%ymm14,%%ymm8 ,%%ymm8 	\n\t"\
		"												\n\t	vsubpd	     %%ymm15,%%ymm9 ,%%ymm9 	\n\t"\
		"												\n\t	vsubpd	%%ymm13,%%ymm10,%%ymm10	\n\t"\
		"vmovaps	%%ymm8 ,0x180(%%rsi)				\n\t	vsubpd	%%ymm12,%%ymm11,%%ymm11	\n\t"\
	"vmovaps	%%ymm12,(%%rsi) 	\n\t"/* spill ymm12 to make room for 2.0 */"	vmovaps	(%%r8),%%ymm12 	\n\t"/* two */\
	"vfmadd132pd	%%ymm12,%%ymm0,%%ymm6			\n\t	vfmadd132pd	%%ymm12,%%ymm8 ,%%ymm14	\n\t"\
	"vfmadd132pd	%%ymm12,%%ymm2,%%ymm5			\n\t	vfmadd132pd	%%ymm12,%%ymm10,%%ymm13	\n\t"\
	"vfmadd132pd	%%ymm12,%%ymm1,%%ymm7			\n\t	vfmadd132pd	%%ymm12,%%ymm9 ,%%ymm15	\n\t"\
	"vfmadd132pd	%%ymm12,%%ymm3,%%ymm4			\n\t	vfmadd132pd	(%%rsi),%%ymm11,%%ymm12	\n\t"\
		/* SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r00) */\
		"												\n\t	vmovaps	%%ymm14,0x100(%%rsi)	\n\t"\
		"												\n\t	vmovaps	%%ymm15,0x120(%%rsi)	\n\t"\
		"												\n\t	vmovaps	0x800(%%rsi),%%ymm8 	\n\t"/* isrt2 */\
		"												\n\t	vmovaps	%%ymm10,%%ymm14	\n\t"\
		"												\n\t	vmovaps	%%ymm13,%%ymm15	\n\t"\
		"														vsubpd	%%ymm12,%%ymm10,%%ymm10	\n\t"\
		"														vsubpd	%%ymm11,%%ymm13,%%ymm13	\n\t"\
		"														vaddpd	%%ymm12,%%ymm14,%%ymm14	\n\t"\
		"														vaddpd	%%ymm11,%%ymm15,%%ymm15	\n\t"\
		"vmovaps	0x100(%%rsi),%%ymm11			\n\t	vfnmadd231pd %%ymm8,%%ymm10,%%ymm2	\n\t"\
		"vmovaps	0x120(%%rsi),%%ymm12			\n\t	vfnmadd231pd %%ymm8,%%ymm13,%%ymm3	\n\t"\
		"vsubpd	%%ymm11,%%ymm6,%%ymm6				\n\t	vfnmadd231pd %%ymm8,%%ymm14,%%ymm4	\n\t"\
		"vsubpd	%%ymm9 ,%%ymm0,%%ymm0				\n\t	vfnmadd231pd %%ymm8,%%ymm15,%%ymm5	\n\t"\
		"vsubpd	0x180(%%rsi),%%ymm1 ,%%ymm1 			\n\t	vsubpd	%%ymm12 ,%%ymm7,%%ymm7	\n\t"\
	"vmovaps	(%%r8),%%ymm8	\n\t"/* rcol: *2 => *sqrt2 due to above FMA leaving ymm10,13-15 unmultiplied-by-isrt2: */\
		"vmovaps	%%ymm6 ,0x100(%%rsi)				\n\t	vmovaps	%%ymm2 ,0x140(%%rsi)	\n\t"\
		"vmovaps	%%ymm0 ,0x080(%%rsi)				\n\t	vmovaps	%%ymm5 ,0x0c0(%%rsi)	\n\t"\
		"vmovaps	%%ymm7 ,0x120(%%rsi)				\n\t	vmovaps	%%ymm4 ,0x160(%%rsi)	\n\t"\
		"vmovaps	%%ymm1 ,0x1a0(%%rsi)				\n\t	vmovaps	%%ymm3 ,0x1e0(%%rsi)	\n\t"\
	"vfmadd132pd	%%ymm8 ,%%ymm6 ,%%ymm11			\n\t	vfmadd132pd	0x20(%%r8),%%ymm2 ,%%ymm10	\n\t"\
	"vfmadd132pd	%%ymm8 ,%%ymm0 ,%%ymm9 			\n\t	vfmadd132pd	0x20(%%r8),%%ymm5 ,%%ymm15	\n\t"\
	"vfmadd132pd	%%ymm8 ,%%ymm7 ,%%ymm12			\n\t	vfmadd132pd	0x20(%%r8),%%ymm4 ,%%ymm14	\n\t"\
	"vfmadd132pd	0x180(%%rsi),%%ymm1 ,%%ymm8 	\n\t	vfmadd132pd	0x20(%%r8),%%ymm3 ,%%ymm13	\n\t"\
		"vmovaps	%%ymm11,     (%%rsi)				\n\t	vmovaps	%%ymm10,0x040(%%rsi)	\n\t"\
		"vmovaps	%%ymm9 ,0x180(%%rsi)				\n\t	vmovaps	%%ymm15,0x1c0(%%rsi)	\n\t"\
		"vmovaps	%%ymm12,0x020(%%rsi)				\n\t	vmovaps	%%ymm14,0x060(%%rsi)	\n\t"\
		"vmovaps	%%ymm8 ,0x0a0(%%rsi)				\n\t	vmovaps	%%ymm13,0x0e0(%%rsi)	\n\t"\
		"\n\t"\
	/***************************************/\
		"\n\t"\
	/*...Block 2: */\
		"movslq	%[__p02],%%rdi							\n\t"\
		"shlq	$3,%%rdi								\n\t"	/* p04<<3 still in r9 */\
		"addq	%%rdi,%%rax	/* &a[j1+p2] */				\n\t	leaq	(%%r9,%%rax),%%r10	\n\t"\
		"addq	%%rdi,%%rbx								\n\t	leaq	(%%r9,%%rbx),%%r11	\n\t"\
		"addq	%%rdi,%%rcx								\n\t	leaq	(%%r9,%%rcx),%%r12	\n\t"\
		"addq	%%rdi,%%rdx								\n\t	leaq	(%%r9,%%rdx),%%r13	\n\t"\
		/* SSE2_RADIX4_DIF_4TWIDDLE_B(r10,c02) */				/* SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO_B(r18,c06) */\
		"addq	$0x200,%%rsi	/* r10 */				\n\t"\
		"vmovaps	     (%%rax),%%ymm0					\n\t	vmovaps	     (%%r10),%%ymm8 	\n\t"\
		"vmovaps	     (%%rcx),%%ymm4					\n\t	vmovaps	     (%%r12),%%ymm12	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1					\n\t	vmovaps	0x020(%%r10),%%ymm9 	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5					\n\t	vmovaps	0x020(%%r12),%%ymm13	\n\t"\
		"vmovaps	0x8e0(%%rsi),%%ymm6	/* c02 */		\n\t	vmovaps	0x9e0(%%rsi),%%ymm14	/* c06 */	\n\t"\
		"vmovaps	0x900(%%rsi),%%ymm7					\n\t	vmovaps	0xa00(%%rsi),%%ymm15	\n\t"\
		"vmovaps	%%ymm0,%%ymm2						\n\t	vmovaps	%%ymm8 ,%%ymm10	\n\t"\
		"vmovaps	%%ymm1,%%ymm3						\n\t	vmovaps	%%ymm9 ,%%ymm11	\n\t"\
		"vmulpd	%%ymm6,%%ymm0,%%ymm0					\n\t	vmulpd	%%ymm14,%%ymm8 ,%%ymm8 	\n\t"\
		"vmulpd	%%ymm6,%%ymm1,%%ymm1					\n\t	vmulpd	%%ymm14,%%ymm9 ,%%ymm9 	\n\t"\
	"vfnmadd231pd	%%ymm7,%%ymm3,%%ymm0			\n\t	vfnmadd231pd	%%ymm15,%%ymm11,%%ymm8 	\n\t"\
	" vfmadd231pd	%%ymm7,%%ymm2,%%ymm1			\n\t	 vfmadd231pd	%%ymm15,%%ymm10,%%ymm9 	\n\t"\
		"vmovaps	%%ymm4,%%ymm6						\n\t	vmovaps	%%ymm12,%%ymm14	\n\t"\
		"vmovaps	%%ymm5,%%ymm7						\n\t	vmovaps	%%ymm13,%%ymm15	\n\t"\
		"vmulpd	0x920(%%rsi),%%ymm4,%%ymm4	/* c12 */	\n\t	vmulpd	0xa20(%%rsi),%%ymm12,%%ymm12	/* c16 */	\n\t"\
		"vmulpd	0x920(%%rsi),%%ymm5,%%ymm5				\n\t	vmulpd	0xa20(%%rsi),%%ymm13,%%ymm13	\n\t"\
	"vfnmadd231pd	0x940(%%rsi),%%ymm7,%%ymm4		\n\t	vfnmadd231pd	0xa40(%%rsi),%%ymm15,%%ymm12	\n\t"\
	" vfmadd231pd	0x940(%%rsi),%%ymm6,%%ymm5		\n\t	 vfmadd231pd	0xa40(%%rsi),%%ymm14,%%ymm13	\n\t"\
		"vmovaps	%%ymm0,%%ymm2						\n\t	vmovaps	%%ymm8 ,%%ymm10	\n\t"\
		"vmovaps	%%ymm1,%%ymm3						\n\t	vmovaps	%%ymm9 ,%%ymm11	\n\t"\
		"vaddpd	%%ymm4,%%ymm0,%%ymm0					\n\t	vaddpd	%%ymm12,%%ymm8 ,%%ymm8 	\n\t"\
		"vaddpd	%%ymm5,%%ymm1,%%ymm1					\n\t	vaddpd	%%ymm13,%%ymm9 ,%%ymm9 	\n\t"\
		"vsubpd	%%ymm4,%%ymm2,%%ymm2					\n\t	vsubpd	%%ymm12,%%ymm10,%%ymm10	\n\t"\
		"vsubpd	%%ymm5,%%ymm3,%%ymm3					\n\t	vsubpd	%%ymm13,%%ymm11,%%ymm11	\n\t"\
		"vmovaps	     (%%rdx),%%ymm4					\n\t	vmovaps	     (%%r13),%%ymm12	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm5					\n\t	vmovaps	0x020(%%r13),%%ymm13	\n\t"\
		"vmovaps	     (%%rdx),%%ymm6					\n\t	vmovaps	     (%%r13),%%ymm14	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm7					\n\t	vmovaps	0x020(%%r13),%%ymm15	\n\t"\
		"vmulpd	0x9a0(%%rsi),%%ymm4,%%ymm4	/* c1A */	\n\t	vmulpd	0xaa0(%%rsi),%%ymm12,%%ymm12	/* c1E */	\n\t"\
		"vmulpd	0x9a0(%%rsi),%%ymm5,%%ymm5				\n\t	vmulpd	0xaa0(%%rsi),%%ymm13,%%ymm13	\n\t"\
	"vfnmadd231pd	0x9c0(%%rsi),%%ymm7,%%ymm4		\n\t	vfnmadd231pd	0xac0(%%rsi),%%ymm15,%%ymm12	\n\t"\
	" vfmadd231pd	0x9c0(%%rsi),%%ymm6,%%ymm5		\n\t	 vfmadd231pd	0xac0(%%rsi),%%ymm14,%%ymm13	\n\t"\
		"vmovaps	%%ymm5,0x020(%%rsi)					\n\t	vmovaps	%%ymm13,0x120(%%rsi)	\n\t"\
		"vmovaps	%%ymm4,     (%%rsi)					\n\t	vmovaps	%%ymm12,0x100(%%rsi)	/* tmpstr r08 */\n\t"\
		"vmovaps	     (%%rbx),%%ymm4					\n\t	vmovaps	     (%%r11),%%ymm12	\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm5					\n\t	vmovaps	0x020(%%r11),%%ymm13	\n\t"\
		"vmovaps	     (%%rbx),%%ymm6					\n\t	vmovaps	     (%%r11),%%ymm14	\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm7					\n\t	vmovaps	0x020(%%r11),%%ymm15	\n\t"\
		"vmulpd	0x960(%%rsi),%%ymm4,%%ymm4	/* c0A */	\n\t	vmulpd	0xa60(%%rsi),%%ymm12,%%ymm12	/* c0E */	\n\t"\
		"vmulpd	0x960(%%rsi),%%ymm5,%%ymm5				\n\t	vmulpd	0xa60(%%rsi),%%ymm13,%%ymm13	\n\t"\
	"vfnmadd231pd	0x980(%%rsi),%%ymm7,%%ymm4		\n\t	vfnmadd231pd	0xa80(%%rsi),%%ymm15,%%ymm12	\n\t"\
	" vfmadd231pd	0x980(%%rsi),%%ymm6,%%ymm5		\n\t	 vfmadd231pd	0xa80(%%rsi),%%ymm14,%%ymm13	\n\t"\
		"vmovaps	%%ymm5,%%ymm7						\n\t	vmovaps	%%ymm13,%%ymm15	\n\t"\
		"vmovaps	%%ymm4,%%ymm6						\n\t	vmovaps	%%ymm12,%%ymm14	\n\t"\
		"vsubpd	     (%%rsi),%%ymm4,%%ymm4				\n\t	vsubpd	0x100(%%rsi),%%ymm12,%%ymm12	\n\t"\
		"vsubpd	0x020(%%rsi),%%ymm5,%%ymm5				\n\t	vsubpd	0x120(%%rsi),%%ymm13,%%ymm13	\n\t"\
		"vaddpd	     (%%rsi),%%ymm6,%%ymm6				\n\t	vaddpd	0x100(%%rsi),%%ymm14,%%ymm14	\n\t"\
		"vaddpd	0x020(%%rsi),%%ymm7,%%ymm7				\n\t	vaddpd	0x120(%%rsi),%%ymm15,%%ymm15	\n\t"\
		"vsubpd	%%ymm6,%%ymm0,%%ymm0					\n\t	vsubpd	     %%ymm14,%%ymm8 ,%%ymm8 	\n\t"\
		"vsubpd	%%ymm5,%%ymm2,%%ymm2					\n\t	vsubpd	     %%ymm15,%%ymm9 ,%%ymm9 	\n\t	vmovaps	%%ymm8 ,0x180(%%rsi)	\n\t"\
		"vsubpd	%%ymm7,%%ymm1,%%ymm1					\n\t	vsubpd	%%ymm13,%%ymm10,%%ymm10	\n\t"\
		"vsubpd	%%ymm4,%%ymm3,%%ymm3					\n\t	vsubpd	%%ymm12,%%ymm11,%%ymm11	\n\t"\
	"vmovaps	%%ymm12,(%%rsi) 	\n\t"/* spill ymm12 to make room for 2.0 */"	vmovaps	(%%r8),%%ymm12 	\n\t"/* two */\
	"vfmadd132pd	%%ymm12,%%ymm0,%%ymm6			\n\t	vfmadd132pd	%%ymm12,%%ymm8 ,%%ymm14	\n\t"\
	"vfmadd132pd	%%ymm12,%%ymm2,%%ymm5			\n\t	vfmadd132pd	%%ymm12,%%ymm10,%%ymm13	\n\t"\
	"vfmadd132pd	%%ymm12,%%ymm1,%%ymm7			\n\t	vfmadd132pd	%%ymm12,%%ymm9 ,%%ymm15	\n\t"\
	"vfmadd132pd	%%ymm12,%%ymm3,%%ymm4			\n\t	vfmadd132pd	(%%rsi),%%ymm11,%%ymm12	\n\t"\
		/* SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r10) */\
		"												\n\t	vmovaps	%%ymm14,0x100(%%rsi)	\n\t"\
		"												\n\t	vmovaps	%%ymm15,0x120(%%rsi)	\n\t"\
		"												\n\t	vmovaps	0x600(%%rsi),%%ymm8 	\n\t"/* isrt2 */\
		"												\n\t	vmovaps	%%ymm10,%%ymm14	\n\t"\
		"												\n\t	vmovaps	%%ymm13,%%ymm15	\n\t"\
		"														vsubpd	%%ymm12,%%ymm10,%%ymm10	\n\t"\
		"														vsubpd	%%ymm11,%%ymm13,%%ymm13	\n\t"\
		"														vaddpd	%%ymm12,%%ymm14,%%ymm14	\n\t"\
		"														vaddpd	%%ymm11,%%ymm15,%%ymm15	\n\t"\
		"vmovaps	0x100(%%rsi),%%ymm11			\n\t	vfnmadd231pd %%ymm8,%%ymm10,%%ymm2	\n\t"\
		"vmovaps	0x120(%%rsi),%%ymm12			\n\t	vfnmadd231pd %%ymm8,%%ymm13,%%ymm3	\n\t"\
		"vsubpd	%%ymm11,%%ymm6,%%ymm6				\n\t	vfnmadd231pd %%ymm8,%%ymm14,%%ymm4	\n\t"\
		"vsubpd	%%ymm9 ,%%ymm0,%%ymm0				\n\t	vfnmadd231pd %%ymm8,%%ymm15,%%ymm5	\n\t"\
		"vsubpd	0x180(%%rsi),%%ymm1 ,%%ymm1 			\n\t	vsubpd	%%ymm12 ,%%ymm7,%%ymm7	\n\t"\
	"vmovaps	(%%r8),%%ymm8	\n\t"/* rcol: *2 => *sqrt2 due to above FMA leaving ymm10,13-15 unmultiplied-by-isrt2: */\
		"vmovaps	%%ymm6 ,0x100(%%rsi)				\n\t	vmovaps	%%ymm2 ,0x140(%%rsi)	\n\t"\
		"vmovaps	%%ymm0 ,0x080(%%rsi)				\n\t	vmovaps	%%ymm5 ,0x0c0(%%rsi)	\n\t"\
		"vmovaps	%%ymm7 ,0x120(%%rsi)				\n\t	vmovaps	%%ymm4 ,0x160(%%rsi)	\n\t"\
		"vmovaps	%%ymm1 ,0x1a0(%%rsi)				\n\t	vmovaps	%%ymm3 ,0x1e0(%%rsi)	\n\t"\
	"vfmadd132pd	%%ymm8 ,%%ymm6 ,%%ymm11			\n\t	vfmadd132pd	0x20(%%r8),%%ymm2 ,%%ymm10	\n\t"\
	"vfmadd132pd	%%ymm8 ,%%ymm0 ,%%ymm9 			\n\t	vfmadd132pd	0x20(%%r8),%%ymm5 ,%%ymm15	\n\t"\
	"vfmadd132pd	%%ymm8 ,%%ymm7 ,%%ymm12			\n\t	vfmadd132pd	0x20(%%r8),%%ymm4 ,%%ymm14	\n\t"\
	"vfmadd132pd	0x180(%%rsi),%%ymm1 ,%%ymm8 	\n\t	vfmadd132pd	0x20(%%r8),%%ymm3 ,%%ymm13	\n\t"\
		"vmovaps	%%ymm11,     (%%rsi)				\n\t	vmovaps	%%ymm10,0x040(%%rsi)	\n\t"\
		"vmovaps	%%ymm9 ,0x180(%%rsi)				\n\t	vmovaps	%%ymm15,0x1c0(%%rsi)	\n\t"\
		"vmovaps	%%ymm12,0x020(%%rsi)				\n\t	vmovaps	%%ymm14,0x060(%%rsi)	\n\t"\
		"vmovaps	%%ymm8 ,0x0a0(%%rsi)				\n\t	vmovaps	%%ymm13,0x0e0(%%rsi)	\n\t"\
		"\n\t"\
	/***************************************/\
		"\n\t"\
	/*...Block 3: */\
		"movq	%[__add0],%%rax	\n\t"\
		"movslq	%[__p01],%%rdi	\n\t"/* Do this way [rather than repeatedly +=p1] since array-padding scheme means p2 == p1+p1 not guaranteed. */\
		"movslq	%[__p08],%%rbx	\n\t"/* Can't get these via simple load-one-and-shift-as-needed due to array padding scheme */\
		"movslq	%[__p10],%%rcx							\n\t"\
		"movslq	%[__p18],%%rdx							\n\t"	/* p04<<3 still in r9 */\
		"leaq	(%%rax,%%rdi,8),%%rax	/* &a[j1+p1] */	\n\t	leaq	(%%r9,%%rax),%%r10	\n\t"\
		"leaq	(%%rax,%%rbx,8),%%rbx					\n\t	leaq	(%%r9,%%rbx),%%r11	\n\t"\
		"leaq	(%%rax,%%rcx,8),%%rcx					\n\t	leaq	(%%r9,%%rcx),%%r12	\n\t"\
		"leaq	(%%rax,%%rdx,8),%%rdx					\n\t	leaq	(%%r9,%%rdx),%%r13	\n\t"\
		/* SSE2_RADIX4_DIF_4TWIDDLE_B(r20,c01) */				/* SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO_B(r28,c05) */\
		"addq	$0x200,%%rsi	/* r20 */				\n\t"\
		"vmovaps	     (%%rax),%%ymm0					\n\t	vmovaps	     (%%r10),%%ymm8 	\n\t"\
		"vmovaps	     (%%rcx),%%ymm4					\n\t	vmovaps	     (%%r12),%%ymm12	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1					\n\t	vmovaps	0x020(%%r10),%%ymm9 	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5					\n\t	vmovaps	0x020(%%r12),%%ymm13	\n\t"\
		"vmovaps	0x8e0(%%rsi),%%ymm6	/* c01 */		\n\t	vmovaps	0x9e0(%%rsi),%%ymm14	/* c05 */	\n\t"\
		"vmovaps	0x900(%%rsi),%%ymm7					\n\t	vmovaps	0xa00(%%rsi),%%ymm15	\n\t"\
		"vmovaps	%%ymm0,%%ymm2						\n\t	vmovaps	%%ymm8 ,%%ymm10	\n\t"\
		"vmovaps	%%ymm1,%%ymm3						\n\t	vmovaps	%%ymm9 ,%%ymm11	\n\t"\
		"vmulpd	%%ymm6,%%ymm0,%%ymm0					\n\t	vmulpd	%%ymm14,%%ymm8 ,%%ymm8 	\n\t"\
		"vmulpd	%%ymm6,%%ymm1,%%ymm1					\n\t	vmulpd	%%ymm14,%%ymm9 ,%%ymm9 	\n\t"\
	"vfnmadd231pd	%%ymm7,%%ymm3,%%ymm0			\n\t	vfnmadd231pd	%%ymm15,%%ymm11,%%ymm8 	\n\t"\
	" vfmadd231pd	%%ymm7,%%ymm2,%%ymm1			\n\t	 vfmadd231pd	%%ymm15,%%ymm10,%%ymm9 	\n\t"\
		"vmovaps	%%ymm4,%%ymm6						\n\t	vmovaps	%%ymm12,%%ymm14	\n\t"\
		"vmovaps	%%ymm5,%%ymm7						\n\t	vmovaps	%%ymm13,%%ymm15	\n\t"\
		"vmulpd	0x920(%%rsi),%%ymm4,%%ymm4	/* c11 */	\n\t	vmulpd	0xa20(%%rsi),%%ymm12,%%ymm12	/* c15 */	\n\t"\
		"vmulpd	0x920(%%rsi),%%ymm5,%%ymm5				\n\t	vmulpd	0xa20(%%rsi),%%ymm13,%%ymm13	\n\t"\
	"vfnmadd231pd	0x940(%%rsi),%%ymm7,%%ymm4		\n\t	vfnmadd231pd	0xa40(%%rsi),%%ymm15,%%ymm12	\n\t"\
	" vfmadd231pd	0x940(%%rsi),%%ymm6,%%ymm5		\n\t	 vfmadd231pd	0xa40(%%rsi),%%ymm14,%%ymm13	\n\t"\
		"vmovaps	%%ymm0,%%ymm2						\n\t	vmovaps	%%ymm8 ,%%ymm10	\n\t"\
		"vmovaps	%%ymm1,%%ymm3						\n\t	vmovaps	%%ymm9 ,%%ymm11	\n\t"\
		"vaddpd	%%ymm4,%%ymm0,%%ymm0					\n\t	vaddpd	%%ymm12,%%ymm8 ,%%ymm8 	\n\t"\
		"vaddpd	%%ymm5,%%ymm1,%%ymm1					\n\t	vaddpd	%%ymm13,%%ymm9 ,%%ymm9 	\n\t"\
		"vsubpd	%%ymm4,%%ymm2,%%ymm2					\n\t	vsubpd	%%ymm12,%%ymm10,%%ymm10	\n\t"\
		"vsubpd	%%ymm5,%%ymm3,%%ymm3					\n\t	vsubpd	%%ymm13,%%ymm11,%%ymm11	\n\t"\
		"vmovaps	     (%%rdx),%%ymm4					\n\t	vmovaps	     (%%r13),%%ymm12	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm5					\n\t	vmovaps	0x020(%%r13),%%ymm13	\n\t"\
		"vmovaps	     (%%rdx),%%ymm6					\n\t	vmovaps	     (%%r13),%%ymm14	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm7					\n\t	vmovaps	0x020(%%r13),%%ymm15	\n\t"\
		"vmulpd	0x9a0(%%rsi),%%ymm4,%%ymm4	/* c19 */	\n\t	vmulpd	0xaa0(%%rsi),%%ymm12,%%ymm12	/* c1D */	\n\t"\
		"vmulpd	0x9a0(%%rsi),%%ymm5,%%ymm5				\n\t	vmulpd	0xaa0(%%rsi),%%ymm13,%%ymm13	\n\t"\
	"vfnmadd231pd	0x9c0(%%rsi),%%ymm7,%%ymm4		\n\t	vfnmadd231pd	0xac0(%%rsi),%%ymm15,%%ymm12	\n\t"\
	" vfmadd231pd	0x9c0(%%rsi),%%ymm6,%%ymm5		\n\t	 vfmadd231pd	0xac0(%%rsi),%%ymm14,%%ymm13	\n\t"\
		"vmovaps	%%ymm5,0x020(%%rsi)					\n\t	vmovaps	%%ymm13,0x120(%%rsi)	\n\t"\
		"vmovaps	%%ymm4,     (%%rsi)					\n\t	vmovaps	%%ymm12,0x100(%%rsi)	/* tmpstr r08 */\n\t"\
		"vmovaps	     (%%rbx),%%ymm4					\n\t	vmovaps	     (%%r11),%%ymm12	\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm5					\n\t	vmovaps	0x020(%%r11),%%ymm13	\n\t"\
		"vmovaps	     (%%rbx),%%ymm6					\n\t	vmovaps	     (%%r11),%%ymm14	\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm7					\n\t	vmovaps	0x020(%%r11),%%ymm15	\n\t"\
		"vmulpd	0x960(%%rsi),%%ymm4,%%ymm4	/* c09 */	\n\t	vmulpd	0xa60(%%rsi),%%ymm12,%%ymm12	/* c0D */	\n\t"\
		"vmulpd	0x960(%%rsi),%%ymm5,%%ymm5				\n\t	vmulpd	0xa60(%%rsi),%%ymm13,%%ymm13	\n\t"\
	"vfnmadd231pd	0x980(%%rsi),%%ymm7,%%ymm4		\n\t	vfnmadd231pd	0xa80(%%rsi),%%ymm15,%%ymm12	\n\t"\
	" vfmadd231pd	0x980(%%rsi),%%ymm6,%%ymm5		\n\t	 vfmadd231pd	0xa80(%%rsi),%%ymm14,%%ymm13	\n\t"\
		"vmovaps	%%ymm5,%%ymm7						\n\t	vmovaps	%%ymm13,%%ymm15	\n\t"\
		"vmovaps	%%ymm4,%%ymm6						\n\t	vmovaps	%%ymm12,%%ymm14	\n\t"\
		"vsubpd	     (%%rsi),%%ymm4,%%ymm4				\n\t	vsubpd	0x100(%%rsi),%%ymm12,%%ymm12	\n\t"\
		"vsubpd	0x020(%%rsi),%%ymm5,%%ymm5				\n\t	vsubpd	0x120(%%rsi),%%ymm13,%%ymm13	\n\t"\
		"vaddpd	     (%%rsi),%%ymm6,%%ymm6				\n\t	vaddpd	0x100(%%rsi),%%ymm14,%%ymm14	\n\t"\
		"vaddpd	0x020(%%rsi),%%ymm7,%%ymm7				\n\t	vaddpd	0x120(%%rsi),%%ymm15,%%ymm15	\n\t"\
		"vsubpd	%%ymm6,%%ymm0,%%ymm0					\n\t	vsubpd	     %%ymm14,%%ymm8 ,%%ymm8 	\n\t"\
		"vsubpd	%%ymm5,%%ymm2,%%ymm2					\n\t	vsubpd	     %%ymm15,%%ymm9 ,%%ymm9 	\n\t	vmovaps	%%ymm8 ,0x180(%%rsi)	\n\t"\
		"vsubpd	%%ymm7,%%ymm1,%%ymm1					\n\t	vsubpd	%%ymm13,%%ymm10,%%ymm10	\n\t"\
		"vsubpd	%%ymm4,%%ymm3,%%ymm3					\n\t	vsubpd	%%ymm12,%%ymm11,%%ymm11	\n\t"\
	"vmovaps	%%ymm12,(%%rsi) 	\n\t"/* spill ymm12 to make room for 2.0 */"	vmovaps	(%%r8),%%ymm12 	\n\t"/* two */\
	"vfmadd132pd	%%ymm12,%%ymm0,%%ymm6			\n\t	vfmadd132pd	%%ymm12,%%ymm8 ,%%ymm14	\n\t"\
	"vfmadd132pd	%%ymm12,%%ymm2,%%ymm5			\n\t	vfmadd132pd	%%ymm12,%%ymm10,%%ymm13	\n\t"\
	"vfmadd132pd	%%ymm12,%%ymm1,%%ymm7			\n\t	vfmadd132pd	%%ymm12,%%ymm9 ,%%ymm15	\n\t"\
	"vfmadd132pd	%%ymm12,%%ymm3,%%ymm4			\n\t	vfmadd132pd	(%%rsi),%%ymm11,%%ymm12	\n\t"\
		/* SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r20) */\
		"												\n\t	vmovaps	%%ymm14,0x100(%%rsi)	\n\t"\
		"												\n\t	vmovaps	%%ymm15,0x120(%%rsi)	\n\t"\
		"												\n\t	vmovaps	0x400(%%rsi),%%ymm8 	\n\t"/* isrt2 */\
		"												\n\t	vmovaps	%%ymm10,%%ymm14	\n\t"\
		"												\n\t	vmovaps	%%ymm13,%%ymm15	\n\t"\
		"														vsubpd	%%ymm12,%%ymm10,%%ymm10	\n\t"\
		"														vsubpd	%%ymm11,%%ymm13,%%ymm13	\n\t"\
		"														vaddpd	%%ymm12,%%ymm14,%%ymm14	\n\t"\
		"														vaddpd	%%ymm11,%%ymm15,%%ymm15	\n\t"\
		"vmovaps	0x100(%%rsi),%%ymm11			\n\t	vfnmadd231pd %%ymm8,%%ymm10,%%ymm2	\n\t"\
		"vmovaps	0x120(%%rsi),%%ymm12			\n\t	vfnmadd231pd %%ymm8,%%ymm13,%%ymm3	\n\t"\
		"vsubpd	%%ymm11,%%ymm6,%%ymm6				\n\t	vfnmadd231pd %%ymm8,%%ymm14,%%ymm4	\n\t"\
		"vsubpd	%%ymm9 ,%%ymm0,%%ymm0				\n\t	vfnmadd231pd %%ymm8,%%ymm15,%%ymm5	\n\t"\
		"vsubpd	0x180(%%rsi),%%ymm1 ,%%ymm1 			\n\t	vsubpd	%%ymm12 ,%%ymm7,%%ymm7	\n\t"\
	"vmovaps	(%%r8),%%ymm8	\n\t"/* rcol: *2 => *sqrt2 due to above FMA leaving ymm10,13-15 unmultiplied-by-isrt2: */\
		"vmovaps	%%ymm6 ,0x100(%%rsi)				\n\t	vmovaps	%%ymm2 ,0x140(%%rsi)	\n\t"\
		"vmovaps	%%ymm0 ,0x080(%%rsi)				\n\t	vmovaps	%%ymm5 ,0x0c0(%%rsi)	\n\t"\
		"vmovaps	%%ymm7 ,0x120(%%rsi)				\n\t	vmovaps	%%ymm4 ,0x160(%%rsi)	\n\t"\
		"vmovaps	%%ymm1 ,0x1a0(%%rsi)				\n\t	vmovaps	%%ymm3 ,0x1e0(%%rsi)	\n\t"\
	"vfmadd132pd	%%ymm8 ,%%ymm6 ,%%ymm11			\n\t	vfmadd132pd	0x20(%%r8),%%ymm2 ,%%ymm10	\n\t"\
	"vfmadd132pd	%%ymm8 ,%%ymm0 ,%%ymm9 			\n\t	vfmadd132pd	0x20(%%r8),%%ymm5 ,%%ymm15	\n\t"\
	"vfmadd132pd	%%ymm8 ,%%ymm7 ,%%ymm12			\n\t	vfmadd132pd	0x20(%%r8),%%ymm4 ,%%ymm14	\n\t"\
	"vfmadd132pd	0x180(%%rsi),%%ymm1 ,%%ymm8 	\n\t	vfmadd132pd	0x20(%%r8),%%ymm3 ,%%ymm13	\n\t"\
		"vmovaps	%%ymm11,     (%%rsi)				\n\t	vmovaps	%%ymm10,0x040(%%rsi)	\n\t"\
		"vmovaps	%%ymm9 ,0x180(%%rsi)				\n\t	vmovaps	%%ymm15,0x1c0(%%rsi)	\n\t"\
		"vmovaps	%%ymm12,0x020(%%rsi)				\n\t	vmovaps	%%ymm14,0x060(%%rsi)	\n\t"\
		"vmovaps	%%ymm8 ,0x0a0(%%rsi)				\n\t	vmovaps	%%ymm13,0x0e0(%%rsi)	\n\t"\
		"\n\t"\
	/***************************************/\
		"\n\t"\
	/*...Block 4: */\
		"movq	%[__add0],%%rax	\n\t"\
		"movslq	%[__p03],%%rdi	\n\t"\
		"movslq	%[__p08],%%rbx	\n\t"\
		"movslq	%[__p10],%%rcx							\n\t"\
		"movslq	%[__p18],%%rdx							\n\t"	/* p04<<3 still in r9 */\
		"leaq	(%%rax,%%rdi,8),%%rax	/* &a[j1+p3] */	\n\t	leaq	(%%r9,%%rax),%%r10	\n\t"\
		"leaq	(%%rax,%%rbx,8),%%rbx					\n\t	leaq	(%%r9,%%rbx),%%r11	\n\t"\
		"leaq	(%%rax,%%rcx,8),%%rcx					\n\t	leaq	(%%r9,%%rcx),%%r12	\n\t"\
		"leaq	(%%rax,%%rdx,8),%%rdx					\n\t	leaq	(%%r9,%%rdx),%%r13	\n\t"\
		/* SSE2_RADIX4_DIF_4TWIDDLE_B(r30,c03) */				/* SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO_B(r38,c07) */\
		"addq	$0x200,%%rsi	/* r30 */				\n\t"\
		"vmovaps	     (%%rax),%%ymm0					\n\t	vmovaps	     (%%r10),%%ymm8 	\n\t"\
		"vmovaps	     (%%rcx),%%ymm4					\n\t	vmovaps	     (%%r12),%%ymm12	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1					\n\t	vmovaps	0x020(%%r10),%%ymm9 	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5					\n\t	vmovaps	0x020(%%r12),%%ymm13	\n\t"\
		"vmovaps	0x8e0(%%rsi),%%ymm6	/* c03 */		\n\t	vmovaps	0x9e0(%%rsi),%%ymm14	/* c07 */	\n\t"\
		"vmovaps	0x900(%%rsi),%%ymm7					\n\t	vmovaps	0xa00(%%rsi),%%ymm15	\n\t"\
		"vmovaps	%%ymm0,%%ymm2						\n\t	vmovaps	%%ymm8 ,%%ymm10	\n\t"\
		"vmovaps	%%ymm1,%%ymm3						\n\t	vmovaps	%%ymm9 ,%%ymm11	\n\t"\
		"vmulpd	%%ymm6,%%ymm0,%%ymm0					\n\t	vmulpd	%%ymm14,%%ymm8 ,%%ymm8 	\n\t"\
		"vmulpd	%%ymm6,%%ymm1,%%ymm1					\n\t	vmulpd	%%ymm14,%%ymm9 ,%%ymm9 	\n\t"\
	"vfnmadd231pd	%%ymm7,%%ymm3,%%ymm0			\n\t	vfnmadd231pd	%%ymm15,%%ymm11,%%ymm8 	\n\t"\
	" vfmadd231pd	%%ymm7,%%ymm2,%%ymm1			\n\t	 vfmadd231pd	%%ymm15,%%ymm10,%%ymm9 	\n\t"\
		"vmovaps	%%ymm4,%%ymm6						\n\t	vmovaps	%%ymm12,%%ymm14	\n\t"\
		"vmovaps	%%ymm5,%%ymm7						\n\t	vmovaps	%%ymm13,%%ymm15	\n\t"\
		"vmulpd	0x920(%%rsi),%%ymm4,%%ymm4	/* c13 */	\n\t	vmulpd	0xa20(%%rsi),%%ymm12,%%ymm12	/* c17 */	\n\t"\
		"vmulpd	0x920(%%rsi),%%ymm5,%%ymm5				\n\t	vmulpd	0xa20(%%rsi),%%ymm13,%%ymm13	\n\t"\
	"vfnmadd231pd	0x940(%%rsi),%%ymm7,%%ymm4		\n\t	vfnmadd231pd	0xa40(%%rsi),%%ymm15,%%ymm12	\n\t"\
	" vfmadd231pd	0x940(%%rsi),%%ymm6,%%ymm5		\n\t	 vfmadd231pd	0xa40(%%rsi),%%ymm14,%%ymm13	\n\t"\
		"vmovaps	%%ymm0,%%ymm2						\n\t	vmovaps	%%ymm8 ,%%ymm10	\n\t"\
		"vmovaps	%%ymm1,%%ymm3						\n\t	vmovaps	%%ymm9 ,%%ymm11	\n\t"\
		"vaddpd	%%ymm4,%%ymm0,%%ymm0					\n\t	vaddpd	%%ymm12,%%ymm8 ,%%ymm8 	\n\t"\
		"vaddpd	%%ymm5,%%ymm1,%%ymm1					\n\t	vaddpd	%%ymm13,%%ymm9 ,%%ymm9 	\n\t"\
		"vsubpd	%%ymm4,%%ymm2,%%ymm2					\n\t	vsubpd	%%ymm12,%%ymm10,%%ymm10	\n\t"\
		"vsubpd	%%ymm5,%%ymm3,%%ymm3					\n\t	vsubpd	%%ymm13,%%ymm11,%%ymm11	\n\t"\
		"vmovaps	     (%%rdx),%%ymm4					\n\t	vmovaps	     (%%r13),%%ymm12	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm5					\n\t	vmovaps	0x020(%%r13),%%ymm13	\n\t"\
		"vmovaps	     (%%rdx),%%ymm6					\n\t	vmovaps	     (%%r13),%%ymm14	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm7					\n\t	vmovaps	0x020(%%r13),%%ymm15	\n\t"\
		"vmulpd	0x9a0(%%rsi),%%ymm4,%%ymm4	/* c1B */	\n\t	vmulpd	0xaa0(%%rsi),%%ymm12,%%ymm12	/* c1F */	\n\t"\
		"vmulpd	0x9a0(%%rsi),%%ymm5,%%ymm5				\n\t	vmulpd	0xaa0(%%rsi),%%ymm13,%%ymm13	\n\t"\
	"vfnmadd231pd	0x9c0(%%rsi),%%ymm7,%%ymm4		\n\t	vfnmadd231pd	0xac0(%%rsi),%%ymm15,%%ymm12	\n\t"\
	" vfmadd231pd	0x9c0(%%rsi),%%ymm6,%%ymm5		\n\t	 vfmadd231pd	0xac0(%%rsi),%%ymm14,%%ymm13	\n\t"\
		"vmovaps	%%ymm5,0x020(%%rsi)					\n\t	vmovaps	%%ymm13,0x120(%%rsi)	\n\t"\
		"vmovaps	%%ymm4,     (%%rsi)					\n\t	vmovaps	%%ymm12,0x100(%%rsi)	/* tmpstr r08 */\n\t"\
		"vmovaps	     (%%rbx),%%ymm4					\n\t	vmovaps	     (%%r11),%%ymm12	\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm5					\n\t	vmovaps	0x020(%%r11),%%ymm13	\n\t"\
		"vmovaps	     (%%rbx),%%ymm6					\n\t	vmovaps	     (%%r11),%%ymm14	\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm7					\n\t	vmovaps	0x020(%%r11),%%ymm15	\n\t"\
		"vmulpd	0x960(%%rsi),%%ymm4,%%ymm4	/* c0B */	\n\t	vmulpd	0xa60(%%rsi),%%ymm12,%%ymm12	/* c0F */	\n\t"\
		"vmulpd	0x960(%%rsi),%%ymm5,%%ymm5				\n\t	vmulpd	0xa60(%%rsi),%%ymm13,%%ymm13	\n\t"\
	"vfnmadd231pd	0x980(%%rsi),%%ymm7,%%ymm4		\n\t	vfnmadd231pd	0xa80(%%rsi),%%ymm15,%%ymm12	\n\t"\
	" vfmadd231pd	0x980(%%rsi),%%ymm6,%%ymm5		\n\t	 vfmadd231pd	0xa80(%%rsi),%%ymm14,%%ymm13	\n\t"\
		"vmovaps	%%ymm5,%%ymm7						\n\t	vmovaps	%%ymm13,%%ymm15	\n\t"\
		"vmovaps	%%ymm4,%%ymm6						\n\t	vmovaps	%%ymm12,%%ymm14	\n\t"\
		"vsubpd	     (%%rsi),%%ymm4,%%ymm4				\n\t	vsubpd	0x100(%%rsi),%%ymm12,%%ymm12	\n\t"\
		"vsubpd	0x020(%%rsi),%%ymm5,%%ymm5				\n\t	vsubpd	0x120(%%rsi),%%ymm13,%%ymm13	\n\t"\
		"vaddpd	     (%%rsi),%%ymm6,%%ymm6				\n\t	vaddpd	0x100(%%rsi),%%ymm14,%%ymm14	\n\t"\
		"vaddpd	0x020(%%rsi),%%ymm7,%%ymm7				\n\t	vaddpd	0x120(%%rsi),%%ymm15,%%ymm15	\n\t"\
		"vsubpd	%%ymm6,%%ymm0,%%ymm0					\n\t	vsubpd	     %%ymm14,%%ymm8 ,%%ymm8 	\n\t"\
		"vsubpd	%%ymm5,%%ymm2,%%ymm2					\n\t	vsubpd	     %%ymm15,%%ymm9 ,%%ymm9 	\n\t	vmovaps	%%ymm8 ,0x180(%%rsi)	\n\t"\
		"vsubpd	%%ymm7,%%ymm1,%%ymm1					\n\t	vsubpd	%%ymm13,%%ymm10,%%ymm10	\n\t"\
		"vsubpd	%%ymm4,%%ymm3,%%ymm3					\n\t	vsubpd	%%ymm12,%%ymm11,%%ymm11	\n\t"\
	"vmovaps	%%ymm12,(%%rsi) 	\n\t"/* spill ymm12 to make room for 2.0 */"	vmovaps	(%%r8),%%ymm12 	\n\t"/* two */\
	"vfmadd132pd	%%ymm12,%%ymm0,%%ymm6			\n\t	vfmadd132pd	%%ymm12,%%ymm8 ,%%ymm14	\n\t"\
	"vfmadd132pd	%%ymm12,%%ymm2,%%ymm5			\n\t	vfmadd132pd	%%ymm12,%%ymm10,%%ymm13	\n\t"\
	"vfmadd132pd	%%ymm12,%%ymm1,%%ymm7			\n\t	vfmadd132pd	%%ymm12,%%ymm9 ,%%ymm15	\n\t"\
	"vfmadd132pd	%%ymm12,%%ymm3,%%ymm4			\n\t	vfmadd132pd	(%%rsi),%%ymm11,%%ymm12	\n\t"\
		/* SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r30) */\
		"												\n\t	vmovaps	%%ymm14,0x100(%%rsi)	\n\t"\
		"												\n\t	vmovaps	%%ymm15,0x120(%%rsi)	\n\t"\
		"												\n\t	vmovaps	0x200(%%rsi),%%ymm8 	\n\t"/* isrt2 */\
		"												\n\t	vmovaps	%%ymm10,%%ymm14	\n\t"\
		"												\n\t	vmovaps	%%ymm13,%%ymm15	\n\t"\
		"														vsubpd	%%ymm12,%%ymm10,%%ymm10	\n\t"\
		"														vsubpd	%%ymm11,%%ymm13,%%ymm13	\n\t"\
		"														vaddpd	%%ymm12,%%ymm14,%%ymm14	\n\t"\
		"														vaddpd	%%ymm11,%%ymm15,%%ymm15	\n\t"\
		"vmovaps	0x100(%%rsi),%%ymm11			\n\t	vfnmadd231pd %%ymm8,%%ymm10,%%ymm2	\n\t"\
		"vmovaps	0x120(%%rsi),%%ymm12			\n\t	vfnmadd231pd %%ymm8,%%ymm13,%%ymm3	\n\t"\
		"vsubpd	%%ymm11,%%ymm6,%%ymm6				\n\t	vfnmadd231pd %%ymm8,%%ymm14,%%ymm4	\n\t"\
		"vsubpd	%%ymm9 ,%%ymm0,%%ymm0				\n\t	vfnmadd231pd %%ymm8,%%ymm15,%%ymm5	\n\t"\
		"vsubpd	0x180(%%rsi),%%ymm1 ,%%ymm1 			\n\t	vsubpd	%%ymm12 ,%%ymm7,%%ymm7	\n\t"\
	"vmovaps	(%%r8),%%ymm8	\n\t"/* rcol: *2 => *sqrt2 due to above FMA leaving ymm10,13-15 unmultiplied-by-isrt2: */\
		"vmovaps	%%ymm6 ,0x100(%%rsi)				\n\t	vmovaps	%%ymm2 ,0x140(%%rsi)	\n\t"\
		"vmovaps	%%ymm0 ,0x080(%%rsi)				\n\t	vmovaps	%%ymm5 ,0x0c0(%%rsi)	\n\t"\
		"vmovaps	%%ymm7 ,0x120(%%rsi)				\n\t	vmovaps	%%ymm4 ,0x160(%%rsi)	\n\t"\
		"vmovaps	%%ymm1 ,0x1a0(%%rsi)				\n\t	vmovaps	%%ymm3 ,0x1e0(%%rsi)	\n\t"\
	"vfmadd132pd	%%ymm8 ,%%ymm6 ,%%ymm11			\n\t	vfmadd132pd	0x20(%%r8),%%ymm2 ,%%ymm10	\n\t"\
	"vfmadd132pd	%%ymm8 ,%%ymm0 ,%%ymm9 			\n\t	vfmadd132pd	0x20(%%r8),%%ymm5 ,%%ymm15	\n\t"\
	"vfmadd132pd	%%ymm8 ,%%ymm7 ,%%ymm12			\n\t	vfmadd132pd	0x20(%%r8),%%ymm4 ,%%ymm14	\n\t"\
	"vfmadd132pd	0x180(%%rsi),%%ymm1 ,%%ymm8 	\n\t	vfmadd132pd	0x20(%%r8),%%ymm3 ,%%ymm13	\n\t"\
		"vmovaps	%%ymm11,     (%%rsi)				\n\t	vmovaps	%%ymm10,0x040(%%rsi)	\n\t"\
		"vmovaps	%%ymm9 ,0x180(%%rsi)				\n\t	vmovaps	%%ymm15,0x1c0(%%rsi)	\n\t"\
		"vmovaps	%%ymm12,0x020(%%rsi)				\n\t	vmovaps	%%ymm14,0x060(%%rsi)	\n\t"\
		"vmovaps	%%ymm8 ,0x0a0(%%rsi)				\n\t	vmovaps	%%ymm13,0x0e0(%%rsi)	\n\t"\
		"\n\t"\
	/**********************************************************************************/\
	/*...and now do eight radix-4 transforms, including the internal twiddle factors: */\
	/**********************************************************************************/\
	/*...Block 1: t00,t10,t20,t30	*/\
		"movq	%[__add0],%%rax	/* &a[j1] */			\n\t"\
		"movslq	%[__p01],%%rbx							\n\t"	/*...Block 5: t08,t18,t28,t38	*/\
		"movslq	%[__p02],%%rcx							\n\t"	/* p04<<3 still in r9 */\
		"movslq	%[__p03],%%rdx							\n\t	leaq	(%%r9,%%rax),%%r10	\n\t"\
		"leaq	(%%rax,%%rbx,8),%%rbx					\n\t	leaq	(%%r9,%%rbx),%%r11	\n\t"\
		"leaq	(%%rax,%%rcx,8),%%rcx					\n\t	leaq	(%%r9,%%rcx),%%r12	\n\t"\
		"leaq	(%%rax,%%rdx,8),%%rdx					\n\t	leaq	(%%r9,%%rdx),%%r13	\n\t"\
		"movq	%[__r00],%%rsi							\n\t	vmovaps	0x800(%%rsi),%%ymm11	\n\t"/* isrt2 */\
		"vmovaps	     (%%rsi),%%ymm0					\n\t	vmovaps	0x500(%%rsi),%%ymm12	\n\t"\
		"vmovaps	0x400(%%rsi),%%ymm4					\n\t	vmovaps	0x520(%%rsi),%%ymm13	\n\t"\
		"vmovaps	0x020(%%rsi),%%ymm1					\n\t	vmovaps	0x700(%%rsi),%%ymm14	\n\t"\
		"vmovaps	0x420(%%rsi),%%ymm5					\n\t	vmovaps	0x720(%%rsi),%%ymm15	\n\t"\
		"vmovaps	0x200(%%rsi),%%ymm2					\n\t	vmulpd	%%ymm11,%%ymm12,%%ymm12	\n\t"\
		"vmovaps	0x600(%%rsi),%%ymm6					\n\t	vmulpd	%%ymm11,%%ymm13,%%ymm13	\n\t"\
		"vmovaps	0x220(%%rsi),%%ymm3					\n\t	vmulpd	%%ymm11,%%ymm14,%%ymm14	\n\t"\
		"vmovaps	0x620(%%rsi),%%ymm7					\n\t	vmulpd	%%ymm11,%%ymm15,%%ymm15	\n\t"\
		"vsubpd	0x200(%%rsi),%%ymm0,%%ymm0				\n\t	vmovaps	0x100(%%rsi),%%ymm8 	\n\t"\
		"vsubpd	0x600(%%rsi),%%ymm4,%%ymm4				\n\t	vmovaps	0x120(%%rsi),%%ymm9 	\n\t"\
		"vsubpd	0x220(%%rsi),%%ymm1,%%ymm1				\n\t	vmovaps	0x300(%%rsi),%%ymm10	\n\t"\
		"vsubpd	0x620(%%rsi),%%ymm5,%%ymm5				\n\t	vmovaps	0x320(%%rsi),%%ymm11	\n\t"\
		"vaddpd	     (%%rsi),%%ymm2,%%ymm2				\n\t	vsubpd	%%ymm11,%%ymm8 ,%%ymm8 	\n\t"\
		"vaddpd	0x400(%%rsi),%%ymm6,%%ymm6				\n\t	vsubpd	%%ymm13,%%ymm12,%%ymm12	\n\t"\
		"vaddpd	0x020(%%rsi),%%ymm3,%%ymm3				\n\t	vsubpd	%%ymm10,%%ymm9 ,%%ymm9 	\n\t"\
		"vaddpd	0x420(%%rsi),%%ymm7,%%ymm7				\n\t	vsubpd	%%ymm14,%%ymm15,%%ymm15	\n\t"\
		"vsubpd	%%ymm6,%%ymm2,%%ymm2				\n\t	vfmadd132pd	(%%r8) ,%%ymm8 ,%%ymm11	\n\t"\
		"vsubpd	%%ymm5,%%ymm0,%%ymm0				\n\t	vfmadd132pd	(%%r8) ,%%ymm12,%%ymm13	\n\t"\
		"vsubpd	%%ymm7,%%ymm3,%%ymm3				\n\t	vfmadd132pd	(%%r8) ,%%ymm9 ,%%ymm10	\n\t"\
		"vsubpd	%%ymm4,%%ymm1,%%ymm1				\n\t	vfmadd132pd	(%%r8) ,%%ymm15,%%ymm14	\n\t"\
		"vmovaps	%%ymm2,     (%%rbx)					\n\t	vmovaps	(%%r8),%%ymm2	\n\t"/* 2.0 */\
	"vfmadd213pd	(%%rbx),%%ymm2,%%ymm6				\n\t	vsubpd	%%ymm14,%%ymm12,%%ymm12	\n\t"\
	"vfmadd132pd	%%ymm2 ,%%ymm0,%%ymm5				\n\t	vsubpd	%%ymm15,%%ymm13,%%ymm13	\n\t"\
	"vfmadd132pd	%%ymm2 ,%%ymm3,%%ymm7			\n\t	 vfmadd132pd	%%ymm2 ,%%ymm12,%%ymm14	\n\t"\
	"vfmadd132pd	%%ymm2 ,%%ymm1,%%ymm4			\n\t	 vfmadd132pd	%%ymm2 ,%%ymm13,%%ymm15	\n\t"\
		/* ymm2-datum already spilled to destination */"\n\t	vsubpd	%%ymm12,%%ymm8 ,%%ymm8 	\n\t"\
		"vmovaps	%%ymm0,     (%%rcx)					\n\t	vsubpd	%%ymm13,%%ymm10,%%ymm10	\n\t"\
		"vmovaps	%%ymm3,0x020(%%rbx)					\n\t	vsubpd	%%ymm14,%%ymm9 ,%%ymm9 	\n\t"\
		"vmovaps	%%ymm1,0x020(%%rdx)					\n\t	vsubpd	%%ymm15,%%ymm11,%%ymm11	\n\t"\
		"vmovaps	%%ymm6,     (%%rax)				\n\t	vfmadd132pd	%%ymm2 ,%%ymm8 ,%%ymm12	\n\t"\
		"vmovaps	%%ymm5,     (%%rdx)				\n\t	vfmadd132pd	%%ymm2 ,%%ymm11,%%ymm15	\n\t"\
		"vmovaps	%%ymm7,0x020(%%rax)				\n\t	vfmadd132pd	%%ymm2 ,%%ymm10,%%ymm13	\n\t"\
		"vmovaps	%%ymm4,0x020(%%rcx)				\n\t	vfmadd132pd	%%ymm2 ,%%ymm9 ,%%ymm14	\n\t"\
		"												\n\t	vmovaps	%%ymm8 ,     (%%r11)	\n\t"\
		"												\n\t	vmovaps	%%ymm11,     (%%r12)	\n\t"\
		"												\n\t	vmovaps	%%ymm10,0x020(%%r11)	\n\t"\
		"												\n\t	vmovaps	%%ymm9 ,0x020(%%r13)	\n\t"\
		"												\n\t	vmovaps	%%ymm12,     (%%r10)	\n\t"\
		"												\n\t	vmovaps	%%ymm15,     (%%r13)	\n\t"\
		"												\n\t	vmovaps	%%ymm13,0x020(%%r10)	\n\t"\
		"												\n\t	vmovaps	%%ymm14,0x020(%%r12)	\n\t"\
	/*...Block 3: t04,t14,t24,t34	*/\
		"addq	$0x080,%%rsi	/* r04 */				\n\t"\
		"subq	%%rax,%%rbx		\n\t"/* p01 << 3 ... Note the fact that the 3 subtracts leave the pointerfied (left- */\
		"subq	%%rax,%%rcx		\n\t"/* p02 << 3 ... shifted) offset means no ',8' (aka << 3) needed in last 3 LEAs. */\
		"subq	%%rax,%%rdx		\n\t"/* p03 << 3 */				/*...Block 7: t0C,t1C,t2C,t3C	*/\
		"movslq	%[__p08],%%rdi							\n\t"	/* p04<<3 still in r9 */\
		"leaq	(%%rax,%%rdi,8),%%rax	/* &a[j1+p8] */	\n\t	leaq	(%%r9,%%rax),%%r10	\n\t"\
		"leaq	(%%rax,%%rbx),%%rbx						\n\t	leaq	(%%r9,%%rbx),%%r11	\n\t"\
		"leaq	(%%rax,%%rcx),%%rcx						\n\t	leaq	(%%r9,%%rcx),%%r12	\n\t"\
		"leaq	(%%rax,%%rdx),%%rdx						\n\t	leaq	(%%r9,%%rdx),%%r13	\n\t"\
		"vmovaps	0x400(%%rsi),%%ymm4					\n\t	vmovaps	0x500(%%rsi),%%ymm12	\n\t"\
		"vmovaps	0x420(%%rsi),%%ymm5					\n\t	vmovaps	0x520(%%rsi),%%ymm13	\n\t"\
		"vmovaps	0x7a0(%%rsi),%%ymm3	/* cc0 */		\n\t"\
		"vmovaps	0x7c0(%%rsi),%%ymm2					\n\t"\
		"vmovaps	%%ymm4,%%ymm6						\n\t	vmovaps	%%ymm12,%%ymm14	\n\t"\
		"vmovaps	%%ymm5,%%ymm7						\n\t	vmovaps	%%ymm13,%%ymm15	\n\t"\
		"vmulpd	%%ymm3,%%ymm4,%%ymm4					\n\t	vmulpd	%%ymm2 ,%%ymm12,%%ymm12	\n\t"\
		"vmulpd	%%ymm3,%%ymm5,%%ymm5					\n\t	vmulpd	%%ymm2 ,%%ymm13,%%ymm13	\n\t"\
		"vmovaps	0x600(%%rsi),%%ymm0					\n\t	vmovaps	0x700(%%rsi),%%ymm8 	\n\t"\
		"vmovaps	0x620(%%rsi),%%ymm1					\n\t	vmovaps	0x720(%%rsi),%%ymm9 	\n\t"\
	"vfnmadd231pd	%%ymm2,%%ymm7,%%ymm4			\n\t	vfnmadd231pd	%%ymm3 ,%%ymm15,%%ymm12	\n\t"\
	" vfmadd231pd	%%ymm2,%%ymm6,%%ymm5			\n\t	 vfmadd231pd	%%ymm3 ,%%ymm14,%%ymm13	\n\t"\
		"vmovaps	%%ymm0,%%ymm6						\n\t	vmovaps	%%ymm8 ,%%ymm14	\n\t"\
		"vmovaps	%%ymm1,%%ymm7						\n\t	vmovaps	%%ymm9 ,%%ymm15	\n\t"\
		"vmulpd	%%ymm2,%%ymm6,%%ymm6					\n\t	vmulpd	%%ymm3 ,%%ymm14,%%ymm14	\n\t"\
		"vmulpd	%%ymm2,%%ymm7,%%ymm7					\n\t	vmulpd	%%ymm3 ,%%ymm15,%%ymm15	\n\t"\
	"vfnmadd231pd	%%ymm3,%%ymm1,%%ymm6			\n\t	vfnmadd231pd	%%ymm2 ,%%ymm9 ,%%ymm14	\n\t"\
	" vfmadd231pd	%%ymm3,%%ymm0,%%ymm7			\n\t	 vfmadd231pd	%%ymm2 ,%%ymm8 ,%%ymm15	\n\t"\
		"vmovaps	%%ymm4,%%ymm2						\n\t	vmovaps	%%ymm12,%%ymm10	\n\t"\
		"vmovaps	%%ymm5,%%ymm3						\n\t	vmovaps	%%ymm13,%%ymm11	\n\t"\
		"vsubpd	%%ymm6,%%ymm4,%%ymm4					\n\t	vsubpd	%%ymm14,%%ymm12,%%ymm12	\n\t"\
		"vsubpd	%%ymm7,%%ymm5,%%ymm5					\n\t	vsubpd	%%ymm15,%%ymm13,%%ymm13	\n\t"\
		"vaddpd	%%ymm2,%%ymm6,%%ymm6					\n\t	vaddpd	%%ymm10,%%ymm14,%%ymm14	\n\t"\
		"vaddpd	%%ymm3,%%ymm7,%%ymm7					\n\t	vaddpd	%%ymm11,%%ymm15,%%ymm15	\n\t"\
		"vmovaps	0x200(%%rsi),%%ymm2					\n\t	vmovaps	0x300(%%rsi),%%ymm10	\n\t"\
		"vmovaps	0x220(%%rsi),%%ymm3					\n\t	vmovaps	0x320(%%rsi),%%ymm11	\n\t"\
		"vmovaps	0x780(%%rsi),%%ymm1	\n\t"/* isrt2 */\
		"vmovaps	%%ymm2,%%ymm0						\n\t	vmovaps	%%ymm10,%%ymm8 	\n\t"\
		"vsubpd	%%ymm3,%%ymm2,%%ymm2					\n\t	vaddpd	%%ymm11,%%ymm10,%%ymm10	\n\t"\
		"vaddpd	%%ymm0,%%ymm3,%%ymm3					\n\t	vsubpd	%%ymm8 ,%%ymm11,%%ymm11	\n\t"\
		"vmulpd	%%ymm1,%%ymm2,%%ymm2					\n\t	vmulpd	%%ymm1 ,%%ymm10,%%ymm10	\n\t"\
		"vmulpd	%%ymm1,%%ymm3,%%ymm3					\n\t	vmulpd	%%ymm1 ,%%ymm11,%%ymm11	\n\t"\
		"vmovaps	     (%%rsi),%%ymm0					\n\t	vmovaps	0x100(%%rsi),%%ymm8 	\n\t"\
		"vmovaps	0x020(%%rsi),%%ymm1					\n\t	vmovaps	0x120(%%rsi),%%ymm9 	\n\t"\
		"vsubpd	%%ymm2,%%ymm0,%%ymm0					\n\t	vsubpd	%%ymm10,%%ymm8 ,%%ymm8 	\n\t"\
		"vsubpd	%%ymm3,%%ymm1,%%ymm1					\n\t	vsubpd	%%ymm11,%%ymm9 ,%%ymm9 	\n\t"\
	"vfmadd132pd	(%%r8),%%ymm0,%%ymm2			\n\t	vfmadd132pd	(%%r8) ,%%ymm8 ,%%ymm10	\n\t"\
	"vfmadd132pd	(%%r8),%%ymm1,%%ymm3			\n\t	vfmadd132pd	(%%r8) ,%%ymm9 ,%%ymm11	\n\t"\
		"vsubpd	%%ymm6,%%ymm2,%%ymm2					\n\t	vsubpd	%%ymm12,%%ymm8 ,%%ymm8 	\n\t"\
		"vsubpd	%%ymm5,%%ymm0,%%ymm0					\n\t	vsubpd	%%ymm15,%%ymm10,%%ymm10	\n\t"\
		"vsubpd	%%ymm7,%%ymm3,%%ymm3					\n\t	vsubpd	%%ymm13,%%ymm9 ,%%ymm9 	\n\t"\
		"vsubpd	%%ymm4,%%ymm1,%%ymm1					\n\t	vsubpd	%%ymm14,%%ymm11,%%ymm11	\n\t"\
		"vmovaps	%%ymm2,     (%%rbx)					\n\t	vmovaps	%%ymm8 ,     (%%r11)	\n\t"\
		"vmovaps	%%ymm0,     (%%rcx)					\n\t	vmovaps	%%ymm10,     (%%r12)	\n\t"\
		"vmovaps	%%ymm3,0x020(%%rbx)					\n\t	vmovaps	%%ymm9 ,0x020(%%r11)	\n\t"\
		"vmovaps	%%ymm1,0x020(%%rdx)					\n\t	vmovaps	%%ymm11,0x020(%%r13)	\n\t"\
	"vmovaps	(%%r8),%%ymm2	\n\t"/* 2.0 */\
	"vfmadd213pd	(%%rbx),%%ymm2,%%ymm6			\n\t	vfmadd132pd	%%ymm2 ,%%ymm8 ,%%ymm12	\n\t"\
	"vfmadd132pd	%%ymm2 ,%%ymm0,%%ymm5			\n\t	vfmadd132pd	%%ymm2 ,%%ymm10,%%ymm15	\n\t"\
	"vfmadd132pd	%%ymm2 ,%%ymm3,%%ymm7			\n\t	vfmadd132pd	%%ymm2 ,%%ymm9 ,%%ymm13	\n\t"\
	"vfmadd132pd	%%ymm2 ,%%ymm1,%%ymm4			\n\t	vfmadd132pd	%%ymm2 ,%%ymm11,%%ymm14	\n\t"\
		"vmovaps	%%ymm6,     (%%rax)					\n\t	vmovaps	%%ymm12,     (%%r10)	\n\t"\
		"vmovaps	%%ymm5,     (%%rdx)					\n\t	vmovaps	%%ymm15,     (%%r13)	\n\t"\
		"vmovaps	%%ymm7,0x020(%%rax)					\n\t	vmovaps	%%ymm13,0x020(%%r10)	\n\t"\
		"vmovaps	%%ymm4,0x020(%%rcx)					\n\t	vmovaps	%%ymm14,0x020(%%r12)	\n\t"\
		"\n\t"\
	/*...Block 2: t02,t12,t22,t32	*/\
		"subq	$0x040,%%rsi	/* r02 */				\n\t"\
		"subq	%%rax,%%rbx								\n\t"\
		"subq	%%rax,%%rcx								\n\t"\
		"subq	%%rax,%%rdx								\n\t"\
		"movq	%[__add0],%%rax	/* &a[j1] */			\n\t"	/*...Block 6: t0A,t1A,t2A,t3A	*/\
		"movslq	%[__p10],%%rdi							\n\t"	/* p04<<3 still in r9 */\
		"leaq	(%%rax,%%rdi,8),%%rax	/* &a[j1+p10] */\n\t	leaq	(%%r9,%%rax),%%r10	\n\t"\
		"leaq	(%%rax,%%rbx),%%rbx						\n\t	leaq	(%%r9,%%rbx),%%r11	\n\t"\
		"leaq	(%%rax,%%rcx),%%rcx						\n\t	leaq	(%%r9,%%rcx),%%r12	\n\t"\
		"leaq	(%%rax,%%rdx),%%rdx						\n\t	leaq	(%%r9,%%rdx),%%r13	\n\t"\
		"vmovaps	0x400(%%rsi),%%ymm4					\n\t	vmovaps	0x500(%%rsi),%%ymm12	\n\t"\
		"vmovaps	0x420(%%rsi),%%ymm5					\n\t	vmovaps	0x520(%%rsi),%%ymm13	\n\t"\
		"vmovaps	0x820(%%rsi),%%ymm2	/* cc1 */		\n\t	vmovaps	0x860(%%rsi),%%ymm11	/* cc3 */	\n\t"\
		"vmovaps	0x840(%%rsi),%%ymm3					\n\t	vmovaps	0x880(%%rsi),%%ymm10	\n\t"\
		"vmovaps	%%ymm4,%%ymm6						\n\t	vmovaps	%%ymm12,%%ymm14	\n\t"\
		"vmovaps	%%ymm5,%%ymm7						\n\t	vmovaps	%%ymm13,%%ymm15	\n\t"\
		"vmulpd	%%ymm2,%%ymm4,%%ymm4					\n\t	vmulpd	%%ymm10,%%ymm12,%%ymm12	\n\t"\
		"vmulpd	%%ymm2,%%ymm5,%%ymm5					\n\t	vmulpd	%%ymm10,%%ymm13,%%ymm13	\n\t"\
		"vmovaps	0x600(%%rsi),%%ymm0					\n\t	vmovaps	0x700(%%rsi),%%ymm8 	\n\t"\
		"vmovaps	0x620(%%rsi),%%ymm1					\n\t	vmovaps	0x720(%%rsi),%%ymm9 	\n\t"\
	"vfnmadd231pd	%%ymm3,%%ymm7,%%ymm4			\n\t	vfnmadd231pd	%%ymm11,%%ymm15,%%ymm12	\n\t"\
	" vfmadd231pd	%%ymm3,%%ymm6,%%ymm5			\n\t	 vfmadd231pd	%%ymm11,%%ymm14,%%ymm13	\n\t"\
		"vmovaps	%%ymm0,%%ymm6						\n\t	vmovaps	%%ymm8 ,%%ymm14	\n\t"\
		"vmovaps	%%ymm1,%%ymm7						\n\t	vmovaps	%%ymm9 ,%%ymm15	\n\t"\
		"vmulpd	%%ymm11,%%ymm6,%%ymm6					\n\t	vmulpd	%%ymm2 ,%%ymm14,%%ymm14	\n\t"\
		"vmulpd	%%ymm11,%%ymm7,%%ymm7					\n\t	vmulpd	%%ymm2 ,%%ymm15,%%ymm15	\n\t"\
	"vfnmadd231pd	%%ymm10,%%ymm1,%%ymm6			\n\t	 vfmadd231pd	%%ymm3 ,%%ymm9 ,%%ymm14	\n\t"\
	" vfmadd231pd	%%ymm10,%%ymm0,%%ymm7			\n\t	vfnmadd231pd	%%ymm3 ,%%ymm8 ,%%ymm15	\n\t"\
		"vmovaps	%%ymm4,%%ymm2						\n\t	vmovaps	%%ymm12,%%ymm10	\n\t"\
		"vmovaps	%%ymm5,%%ymm3						\n\t	vmovaps	%%ymm13,%%ymm11	\n\t"\
		"vsubpd	%%ymm6,%%ymm4,%%ymm4					\n\t	vsubpd	%%ymm14,%%ymm12,%%ymm12	\n\t"\
		"vsubpd	%%ymm7,%%ymm5,%%ymm5					\n\t	vsubpd	%%ymm15,%%ymm13,%%ymm13	\n\t"\
		"vaddpd	%%ymm2,%%ymm6,%%ymm6					\n\t	vaddpd	%%ymm10,%%ymm14,%%ymm14	\n\t"\
		"vaddpd	%%ymm3,%%ymm7,%%ymm7					\n\t	vaddpd	%%ymm11,%%ymm15,%%ymm15	\n\t"\
		"vmovaps	0x200(%%rsi),%%ymm1					\n\t	vmovaps	0x300(%%rsi),%%ymm9 	\n\t"\
		"vmovaps	0x220(%%rsi),%%ymm3					\n\t	vmovaps	0x320(%%rsi),%%ymm11	\n\t"\
		"vmovaps	0x800(%%rsi),%%ymm0	/* ss0 */		\n\t	vmovaps	0x7e0(%%rsi),%%ymm8 		/* cc0 */	\n\t"\
		"vmovaps	%%ymm1,%%ymm2						\n\t	vmovaps	%%ymm9 ,%%ymm10	\n\t"\
		"vmulpd	      %%ymm0,%%ymm1,%%ymm1				\n\t	vmulpd	     %%ymm8 ,%%ymm9 ,%%ymm9 	\n\t"\
		"vmulpd	      %%ymm3,%%ymm0,%%ymm0				\n\t	vmulpd	     %%ymm11,%%ymm8 ,%%ymm8 	\n\t"\
	"vfmsub132pd	0x7e0(%%rsi),%%ymm0,%%ymm2		\n\t	vfmadd132pd	0x800(%%rsi),%%ymm8 ,%%ymm10	\n\t"\
	"vfmadd132pd	0x7e0(%%rsi),%%ymm1,%%ymm3		\n\t	vfmsub132pd	0x800(%%rsi),%%ymm9 ,%%ymm11	\n\t"\
		"vmovaps	     (%%rsi),%%ymm0					\n\t	vmovaps	0x100(%%rsi),%%ymm8 	\n\t"\
		"vmovaps	0x020(%%rsi),%%ymm1					\n\t	vmovaps	0x120(%%rsi),%%ymm9 	\n\t"\
		"vsubpd	%%ymm2,%%ymm0,%%ymm0					\n\t	vsubpd	%%ymm10,%%ymm8 ,%%ymm8 	\n\t"\
		"vsubpd	%%ymm3,%%ymm1,%%ymm1					\n\t	vsubpd	%%ymm11,%%ymm9 ,%%ymm9 	\n\t"\
	"vfmadd132pd	(%%r8),%%ymm0,%%ymm2			\n\t	vfmadd132pd	(%%r8) ,%%ymm8 ,%%ymm10	\n\t"\
	"vfmadd132pd	(%%r8),%%ymm1,%%ymm3			\n\t	vfmadd132pd	(%%r8) ,%%ymm9 ,%%ymm11	\n\t"\
		"vsubpd	%%ymm6,%%ymm2,%%ymm2					\n\t	vsubpd	%%ymm12,%%ymm8 ,%%ymm8 	\n\t"\
		"vsubpd	%%ymm5,%%ymm0,%%ymm0					\n\t	vsubpd	%%ymm15,%%ymm10,%%ymm10	\n\t"\
		"vsubpd	%%ymm7,%%ymm3,%%ymm3					\n\t	vsubpd	%%ymm13,%%ymm9 ,%%ymm9 	\n\t"\
		"vsubpd	%%ymm4,%%ymm1,%%ymm1					\n\t	vsubpd	%%ymm14,%%ymm11,%%ymm11	\n\t"\
		"vmovaps	%%ymm2,     (%%rbx)					\n\t	vmovaps	%%ymm8 ,     (%%r11)	\n\t"\
		"vmovaps	%%ymm0,     (%%rcx)					\n\t	vmovaps	%%ymm10,     (%%r12)	\n\t"\
		"vmovaps	%%ymm3,0x020(%%rbx)					\n\t	vmovaps	%%ymm9 ,0x020(%%r11)	\n\t"\
		"vmovaps	%%ymm1,0x020(%%rdx)					\n\t	vmovaps	%%ymm11,0x020(%%r13)	\n\t"\
	"vmovaps	(%%r8),%%ymm2	\n\t"/* 2.0 */\
	"vfmadd213pd	(%%rbx),%%ymm2,%%ymm6			\n\t	vfmadd132pd	%%ymm2 ,%%ymm8 ,%%ymm12	\n\t"\
	"vfmadd132pd	%%ymm2 ,%%ymm0,%%ymm5			\n\t	vfmadd132pd	%%ymm2 ,%%ymm10,%%ymm15	\n\t"\
	"vfmadd132pd	%%ymm2 ,%%ymm3,%%ymm7			\n\t	vfmadd132pd	%%ymm2 ,%%ymm9 ,%%ymm13	\n\t"\
	"vfmadd132pd	%%ymm2 ,%%ymm1,%%ymm4			\n\t	vfmadd132pd	%%ymm2 ,%%ymm11,%%ymm14	\n\t"\
		"vmovaps	%%ymm6,     (%%rax)					\n\t	vmovaps	%%ymm12,     (%%r10)	\n\t"\
		"vmovaps	%%ymm5,     (%%rdx)					\n\t	vmovaps	%%ymm15,     (%%r13)	\n\t"\
		"vmovaps	%%ymm7,0x020(%%rax)					\n\t	vmovaps	%%ymm13,0x020(%%r10)	\n\t"\
		"vmovaps	%%ymm4,0x020(%%rcx)					\n\t	vmovaps	%%ymm14,0x020(%%r12)	\n\t"\
		"\n\t"\
	/*...Block 4: t06,t16,t26,t36	*/\
		"addq	$0x080,%%rsi	/* r06 */	\n\t"\
		"subq	%%rax,%%rbx								\n\t"\
		"subq	%%rax,%%rcx								\n\t"\
		"subq	%%rax,%%rdx								\n\t"\
		"movq	%[__add0],%%rax	/* &a[j1] */			\n\t"	/*...Block 8: t0E,t1E,t2E,t3E	*/\
		"movslq	%[__p18],%%rdi							\n\t"	/* p04<<3 still in r9 */\
		"leaq	(%%rax,%%rdi,8),%%rax	/* &a[j1+p18] */\n\t	leaq	(%%r9,%%rax),%%r10	\n\t"\
		"leaq	(%%rax,%%rbx),%%rbx						\n\t	leaq	(%%r9,%%rbx),%%r11	\n\t"\
		"leaq	(%%rax,%%rcx),%%rcx						\n\t	leaq	(%%r9,%%rcx),%%r12	\n\t"\
		"leaq	(%%rax,%%rdx),%%rdx						\n\t	leaq	(%%r9,%%rdx),%%r13	\n\t"\
		"vmovaps	0x400(%%rsi),%%ymm4					\n\t	vmovaps	0x500(%%rsi),%%ymm12	\n\t"\
		"vmovaps	0x420(%%rsi),%%ymm5					\n\t	vmovaps	0x520(%%rsi),%%ymm13	\n\t"\
		"vmovaps	0x7e0(%%rsi),%%ymm2	/* cc3 */		\n\t	vmovaps	0x7a0(%%rsi),%%ymm11 /* cc1 */	\n\t"\
		"vmovaps	0x800(%%rsi),%%ymm3					\n\t	vmovaps	0x7c0(%%rsi),%%ymm10 \n\t"\
		"vmovaps	%%ymm4,%%ymm6						\n\t	vmovaps	%%ymm12,%%ymm14	\n\t"\
		"vmovaps	%%ymm5,%%ymm7						\n\t	vmovaps	%%ymm13,%%ymm15	\n\t"\
		"vmulpd	%%ymm2,%%ymm4,%%ymm4					\n\t	vmulpd	%%ymm10,%%ymm12,%%ymm12	\n\t"\
		"vmulpd	%%ymm2,%%ymm5,%%ymm5					\n\t	vmulpd	%%ymm10,%%ymm13,%%ymm13	\n\t"\
	"vfnmadd231pd	%%ymm3,%%ymm7,%%ymm4			\n\t	vfnmadd231pd	%%ymm11,%%ymm15,%%ymm12	\n\t"\
	" vfmadd231pd	%%ymm3,%%ymm6,%%ymm5			\n\t	 vfmadd231pd	%%ymm11,%%ymm14,%%ymm13	\n\t"\
		"vmovaps	0x600(%%rsi),%%ymm0					\n\t	vmovaps	0x700(%%rsi),%%ymm8 	\n\t"\
		"vmovaps	0x620(%%rsi),%%ymm1					\n\t	vmovaps	0x720(%%rsi),%%ymm9 	\n\t"\
		"vmovaps	%%ymm0,%%ymm6						\n\t	vmovaps	%%ymm8 ,%%ymm14	\n\t"\
		"vmovaps	%%ymm1,%%ymm7						\n\t	vmovaps	%%ymm9 ,%%ymm15	\n\t"\
		"vmulpd	%%ymm10,%%ymm6,%%ymm6					\n\t	vmulpd	%%ymm3 ,%%ymm14,%%ymm14	\n\t"\
		"vmulpd	%%ymm10,%%ymm7,%%ymm7					\n\t	vmulpd	%%ymm3 ,%%ymm15,%%ymm15	\n\t"\
	" vfmadd231pd	%%ymm11,%%ymm1,%%ymm6			\n\t	vfnmadd231pd	%%ymm2 ,%%ymm9 ,%%ymm14	\n\t"\
	"vfnmadd231pd	%%ymm11,%%ymm0,%%ymm7			\n\t	 vfmadd231pd	%%ymm2 ,%%ymm8 ,%%ymm15	\n\t"\
		"vmovaps	%%ymm4,%%ymm2						\n\t	vmovaps	%%ymm12,%%ymm10	\n\t"\
		"vmovaps	%%ymm5,%%ymm3						\n\t	vmovaps	%%ymm13,%%ymm11	\n\t"\
		"vsubpd	%%ymm6,%%ymm4,%%ymm4					\n\t	vsubpd	%%ymm14,%%ymm12,%%ymm12	\n\t"\
		"vsubpd	%%ymm7,%%ymm5,%%ymm5					\n\t	vsubpd	%%ymm15,%%ymm13,%%ymm13	\n\t"\
		"vaddpd	%%ymm2,%%ymm6,%%ymm6					\n\t	vaddpd	%%ymm10,%%ymm14,%%ymm14	\n\t"\
		"vaddpd	%%ymm3,%%ymm7,%%ymm7					\n\t	vaddpd	%%ymm11,%%ymm15,%%ymm15	\n\t"\
		"vmovaps	0x200(%%rsi),%%ymm1					\n\t	vmovaps	0x300(%%rsi),%%ymm9 	\n\t"\
		"vmovaps	0x220(%%rsi),%%ymm3					\n\t	vmovaps	0x320(%%rsi),%%ymm11	\n\t"\
		"vmovaps	0x760(%%rsi),%%ymm0		/* cc0 */	\n\t	vmovaps	0x780(%%rsi),%%ymm8 	/* ss0 */	\n\t"\
		"vmovaps	%%ymm1,%%ymm2						\n\t	vmovaps	%%ymm9 ,%%ymm10	\n\t"\
		"vmulpd	      %%ymm0,%%ymm1,%%ymm1				\n\t	vmulpd	     %%ymm8 ,%%ymm9 ,%%ymm9 	\n\t"\
		"vmulpd	      %%ymm3,%%ymm0,%%ymm0				\n\t	vmulpd	     %%ymm11,%%ymm8 ,%%ymm8 	\n\t"\
	"vfmsub132pd	0x780(%%rsi),%%ymm0,%%ymm2		\n\t	vfmadd132pd	0x760(%%rsi),%%ymm8 ,%%ymm10	\n\t"\
	"vfmadd132pd	0x780(%%rsi),%%ymm1,%%ymm3		\n\t	vfmsub132pd	0x760(%%rsi),%%ymm9 ,%%ymm11	\n\t"\
		"vmovaps	     (%%rsi),%%ymm0					\n\t	vmovaps	0x100(%%rsi),%%ymm8 	\n\t"\
		"vmovaps	0x020(%%rsi),%%ymm1					\n\t	vmovaps	0x120(%%rsi),%%ymm9 	\n\t"\
		"vsubpd	%%ymm2,%%ymm0,%%ymm0					\n\t	vsubpd	%%ymm10,%%ymm8 ,%%ymm8 	\n\t"\
		"vsubpd	%%ymm3,%%ymm1,%%ymm1					\n\t	vsubpd	%%ymm11,%%ymm9 ,%%ymm9 	\n\t"\
	"vfmadd132pd	(%%r8),%%ymm0,%%ymm2			\n\t	vfmadd132pd	(%%r8) ,%%ymm8 ,%%ymm10	\n\t"\
	"vfmadd132pd	(%%r8),%%ymm1,%%ymm3			\n\t	vfmadd132pd	(%%r8) ,%%ymm9 ,%%ymm11	\n\t"\
		"vsubpd	%%ymm4,%%ymm2,%%ymm2					\n\t	vsubpd	%%ymm12,%%ymm8 ,%%ymm8 	\n\t"\
		"vsubpd	%%ymm7,%%ymm0,%%ymm0					\n\t	vsubpd	%%ymm15,%%ymm10,%%ymm10	\n\t"\
		"vsubpd	%%ymm5,%%ymm3,%%ymm3					\n\t	vsubpd	%%ymm13,%%ymm9 ,%%ymm9 	\n\t"\
		"vsubpd	%%ymm6,%%ymm1,%%ymm1					\n\t	vsubpd	%%ymm14,%%ymm11,%%ymm11	\n\t"\
		"vmovaps	%%ymm2,     (%%rbx)					\n\t	vmovaps	%%ymm8 ,     (%%r11)	\n\t"\
		"vmovaps	%%ymm0,     (%%rcx)					\n\t	vmovaps	%%ymm10,     (%%r12)	\n\t"\
		"vmovaps	%%ymm3,0x020(%%rbx)					\n\t	vmovaps	%%ymm9 ,0x020(%%r11)	\n\t"\
		"vmovaps	%%ymm1,0x020(%%rdx)					\n\t	vmovaps	%%ymm11,0x020(%%r13)	\n\t"\
	"vmovaps	(%%r8),%%ymm2	\n\t"/* 2.0 */\
	"vfmadd213pd	(%%rbx),%%ymm2,%%ymm4			\n\t	vfmadd132pd	%%ymm2 ,%%ymm8 ,%%ymm12	\n\t"\
	"vfmadd132pd	%%ymm2 ,%%ymm0,%%ymm7			\n\t	vfmadd132pd	%%ymm2 ,%%ymm10,%%ymm15	\n\t"\
	"vfmadd132pd	%%ymm2 ,%%ymm3,%%ymm5			\n\t	vfmadd132pd	%%ymm2 ,%%ymm9 ,%%ymm13	\n\t"\
	"vfmadd132pd	%%ymm2 ,%%ymm1,%%ymm6			\n\t	vfmadd132pd	%%ymm2 ,%%ymm11,%%ymm14	\n\t"\
		"vmovaps	%%ymm4,     (%%rax)					\n\t	vmovaps	%%ymm12,     (%%r10)	\n\t"\
		"vmovaps	%%ymm7,     (%%rdx)					\n\t	vmovaps	%%ymm15,     (%%r13)	\n\t"\
		"vmovaps	%%ymm5,0x020(%%rax)					\n\t	vmovaps	%%ymm13,0x020(%%r10)	\n\t"\
		"vmovaps	%%ymm6,0x020(%%rcx)					\n\t	vmovaps	%%ymm14,0x020(%%r12)	\n\t"\
		"\n\t"\
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
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r8","r9","r10","r11","r12","r13","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"		/* Clobbered registers */\
	);\
	}

	/*
	For GCC-macro version of this, use that - in terms of vec_dbl pointer-offsets - isrt2 + 1,3,5 = cc0,cc1,cc3,
	isrt2 + 7,15,23,31,39,47,55,63 = c00,04,02,06,01,05,03,07, and isrt2 + 71,72 = one,two
	in order to reduce number of args to <= the GCC-allowed maximum of 30:
	*/
	// In our FMA-ized version, original mix of [220 ADD, 220 SUB, 374 MUL] ==> [14 ADD, 154 SUB, 272 FMA (272 nontrivial), 172 MUL].
	#define SSE2_RADIX32_DIT_TWIDDLE(Xadd0,Xp01,Xp02,Xp03,Xp04,Xp05,Xp06,Xp07,Xp08,Xp10,Xp18,Xr00,Xr10,Xr20,Xr30,Xisrt2)\
	{\
	__asm__ volatile (\
	/*...Block 1: */\
		"movq	%[__r00],%%rsi					\n\t"\
		"movq	%[__add0],%%rax					\n\t"\
		"movslq	%[__p01],%%rbx					\n\t"\
		"movslq	%[__p02],%%rcx					\n\t		movq	%[__isrt2],%%r8	\n\t"\
		"movslq	%[__p03],%%rdx					\n\t		movslq	%[__p04],%%r9	\n\t"\
		"shlq	$3,%%rbx						\n\t		shlq	$3,%%r9		\n\t"\
		"shlq	$3,%%rcx						\n\t		addq	$0x900,%%r8	\n\t"/* two */\
		"shlq	$3,%%rdx						\n\t		leaq	(%%r9,%%rax),%%r10	\n\t"\
		"addq	%%rax,%%rbx						\n\t		leaq	(%%r9,%%rbx),%%r11	\n\t"\
		"addq	%%rax,%%rcx						\n\t		leaq	(%%r9,%%rcx),%%r12	\n\t"\
		"addq	%%rax,%%rdx						\n\t		leaq	(%%r9,%%rdx),%%r13	\n\t"\
		/* SSE2_RADIX4_DIT_0TWIDDLE_B(r00) */			/* SSE2_RADIX4_DIT_0TWIDDLE_2NDOFTWO_B(r08) */\
		"vmovaps	     (%%rax),%%ymm2			\n\t		vmovaps	     (%%r10),%%ymm10		\n\t"\
		"vmovaps	     (%%rcx),%%ymm6			\n\t		vmovaps	     (%%r12),%%ymm14		\n\t"\
		"vmovaps	0x020(%%rax),%%ymm3			\n\t		vmovaps	0x020(%%r10),%%ymm11		\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm7			\n\t		vmovaps	0x020(%%r12),%%ymm15		\n\t"\
		"vmovaps	     (%%rbx),%%ymm0			\n\t		vmovaps	     (%%r11),%%ymm8 		\n\t"\
		"vmovaps	     (%%rdx),%%ymm4			\n\t		vmovaps	     (%%r13),%%ymm12		\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm1			\n\t		vmovaps	0x020(%%r11),%%ymm9 		\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm5			\n\t		vmovaps	0x020(%%r13),%%ymm13		\n\t"\
		"vsubpd	%%ymm0,%%ymm2,%%ymm2			\n\t		vsubpd	%%ymm8 ,%%ymm10,%%ymm10			\n\t"\
		"vsubpd	%%ymm4,%%ymm6,%%ymm6			\n\t		vsubpd	%%ymm12,%%ymm14,%%ymm14			\n\t"\
		"vsubpd	%%ymm1,%%ymm3,%%ymm3			\n\t		vsubpd	%%ymm9 ,%%ymm11,%%ymm11			\n\t"\
		"vsubpd	%%ymm5,%%ymm7,%%ymm7			\n\t		vsubpd	%%ymm13,%%ymm15,%%ymm15			\n\t"\
	"vmovaps	%%ymm13,(%%rsi) 	\n\t"/* spill ymm13 to make room for 2.0 */"	vmovaps	(%%r8),%%ymm13 	\n\t"/* two */\
	"vfmadd132pd	%%ymm13,%%ymm2,%%ymm0		\n\t	vfmadd132pd	%%ymm13,%%ymm10,%%ymm8 	\n\t"\
	"vfmadd132pd	%%ymm13,%%ymm6,%%ymm4		\n\t	vfmadd132pd	%%ymm13,%%ymm14,%%ymm12	\n\t"\
	"vfmadd132pd	%%ymm13,%%ymm3,%%ymm1		\n\t	vfmadd132pd	%%ymm13,%%ymm11,%%ymm9 	\n\t"\
	"vfmadd132pd	%%ymm13,%%ymm7,%%ymm5		\n\t	vfmadd132pd	(%%rsi),%%ymm15,%%ymm13	\n\t"\
		/* Finish radix-4 butterfly: */\
		"vsubpd	%%ymm4,%%ymm0,%%ymm0			\n\t		vsubpd	%%ymm12,%%ymm8 ,%%ymm8 			\n\t"\
		"vsubpd	%%ymm7,%%ymm2,%%ymm2			\n\t		vsubpd	%%ymm13,%%ymm9 ,%%ymm9 			\n\t"\
		"vsubpd	%%ymm5,%%ymm1,%%ymm1			\n\t		vsubpd	%%ymm15,%%ymm10,%%ymm10			\n\t"\
		"vsubpd	%%ymm6,%%ymm3,%%ymm3			\n\t		vsubpd	%%ymm14,%%ymm11,%%ymm11			\n\t"\
	"vmovaps	%%ymm14,(%%rsi) 	\n\t"/* spill ymm14 to make room for 2.0 */"	vmovaps	(%%r8),%%ymm14 	\n\t"/* two */\
	"vfmadd132pd	%%ymm14,%%ymm0,%%ymm4		\n\t	vfmadd132pd	%%ymm14,%%ymm8 ,%%ymm12	\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm2,%%ymm7		\n\t	vfmadd132pd	%%ymm14,%%ymm9 ,%%ymm13	\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm1,%%ymm5		\n\t	vfmadd132pd	%%ymm14,%%ymm10,%%ymm15	\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm3,%%ymm6		\n\t	vfmadd132pd	(%%rsi),%%ymm11,%%ymm14	\n\t"\
		"vsubpd	%%ymm12,%%ymm4	,%%ymm4			\n\t		vsubpd	%%ymm15,%%ymm11,%%ymm11			\n\t"\
		"vsubpd	%%ymm9 ,%%ymm0	,%%ymm0			\n\t		vsubpd	%%ymm10,%%ymm14,%%ymm14			\n\t"\
	"vmovaps	%%ymm8 ,(%%rsi) 	\n\t"/* spill ymm8  to make room for 2.0 */"	vmovaps	(%%r8),%%ymm8  	\n\t"/* two */\
		"vsubpd	%%ymm13,%%ymm5	,%%ymm5			\n\t	vfmadd132pd	%%ymm8 ,%%ymm11,%%ymm15			\n\t"\
		"vsubpd	(%%rsi),%%ymm1	,%%ymm1			\n\t	vfmadd132pd	%%ymm8 ,%%ymm14,%%ymm10			\n\t"\
	"vfmadd132pd	%%ymm8 ,%%ymm4 ,%%ymm12		\n\t	vfnmadd231pd	-0x900(%%r8),%%ymm15,%%ymm7		\n\t"/* isrt2 */\
	"vfmadd132pd	%%ymm8 ,%%ymm0 ,%%ymm9 		\n\t	vfnmadd231pd	-0x900(%%r8),%%ymm14,%%ymm2		\n\t"\
	"vfmadd132pd	%%ymm8 ,%%ymm5 ,%%ymm13		\n\t	vfnmadd231pd	-0x900(%%r8),%%ymm11,%%ymm3		\n\t"\
	"vfmadd132pd	(%%rsi),%%ymm1 ,%%ymm8 		\n\t	vfnmadd231pd	-0x900(%%r8),%%ymm10,%%ymm6		\n\t"\
		/* SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS(r00,r02,r04,r06,r08,r0A,r0C,r0E) */\
	"vmovaps	%%ymm12,(%%rsi) 	\n\t"/* store ymm12 to make room for sqrt2 */"	vmovaps	0x20(%%r8),%%ymm12 	\n\t"/* sqrt2 */\
		"vmovaps	%%ymm4 ,0x100(%%rsi)		\n\t	vfmadd132pd	%%ymm12,%%ymm7 ,%%ymm15	\n\t"/* sqrt2, use of FMA above means we must combine 2*isrt2 here */\
		"vmovaps	%%ymm0 ,0x180(%%rsi)		\n\t	vfmadd132pd	%%ymm12,%%ymm2 ,%%ymm14	\n\t"\
		"vmovaps	%%ymm5 ,0x120(%%rsi)		\n\t	vfmadd132pd	%%ymm12,%%ymm3 ,%%ymm11	\n\t"\
		"vmovaps	%%ymm1 ,0x0a0(%%rsi)		\n\t	vfmadd132pd	%%ymm12,%%ymm6 ,%%ymm10	\n\t"\
		/*vmovaps	%%ymm12,     (%%rsi)*/	"	\n\t		vmovaps	%%ymm7 ,0x140(%%rsi)	\n\t"\
		"vmovaps	%%ymm9 ,0x080(%%rsi)		\n\t		vmovaps	%%ymm2 ,0x1c0(%%rsi)	\n\t"\
		"vmovaps	%%ymm13,0x020(%%rsi)		\n\t		vmovaps	%%ymm3 ,0x160(%%rsi)	\n\t"\
		"vmovaps	%%ymm8 ,0x1a0(%%rsi)		\n\t		vmovaps	%%ymm6 ,0x0e0(%%rsi)	\n\t"\
		"													vmovaps	%%ymm15,0x040(%%rsi)	\n\t"\
		"													vmovaps	%%ymm14,0x0c0(%%rsi)	\n\t"\
		"													vmovaps	%%ymm11,0x060(%%rsi)	\n\t"\
		"													vmovaps	%%ymm10,0x1e0(%%rsi)	\n\t"\
		"\n\t"\
	/*...Block 2: */\
		"addq	$0x200,%%rsi	/* r10 */		\n\t"\
		"movslq	%[__p08],%%rdi					\n\t"\
		"shlq	$3,%%rdi						\n\t"\
		"addq	%%rdi,%%rax	/* add0+p08 */		\n\t		leaq	(%%r9,%%rax),%%r10	\n\t"\
		"addq	%%rdi,%%rbx						\n\t		leaq	(%%r9,%%rbx),%%r11	\n\t"\
		"addq	%%rdi,%%rcx						\n\t		leaq	(%%r9,%%rcx),%%r12	\n\t"\
		"addq	%%rdi,%%rdx						\n\t		leaq	(%%r9,%%rdx),%%r13	\n\t"\
		/* SSE2_RADIX4_DIT_0TWIDDLE_B(r10) */			/* SSE2_RADIX4_DIT_0TWIDDLE_2NDOFTWO_B(r18) */\
		"vmovaps	     (%%rax),%%ymm2			\n\t		vmovaps	     (%%r10),%%ymm10		\n\t"\
		"vmovaps	     (%%rcx),%%ymm6			\n\t		vmovaps	     (%%r12),%%ymm14		\n\t"\
		"vmovaps	0x020(%%rax),%%ymm3			\n\t		vmovaps	0x020(%%r10),%%ymm11		\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm7			\n\t		vmovaps	0x020(%%r12),%%ymm15		\n\t"\
		"vmovaps	     (%%rbx),%%ymm0			\n\t		vmovaps	     (%%r11),%%ymm8 		\n\t"\
		"vmovaps	     (%%rdx),%%ymm4			\n\t		vmovaps	     (%%r13),%%ymm12		\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm1			\n\t		vmovaps	0x020(%%r11),%%ymm9 		\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm5			\n\t		vmovaps	0x020(%%r13),%%ymm13		\n\t"\
		"vsubpd	%%ymm0,%%ymm2,%%ymm2			\n\t		vsubpd	%%ymm8 ,%%ymm10,%%ymm10			\n\t"\
		"vsubpd	%%ymm4,%%ymm6,%%ymm6			\n\t		vsubpd	%%ymm12,%%ymm14,%%ymm14			\n\t"\
		"vsubpd	%%ymm1,%%ymm3,%%ymm3			\n\t		vsubpd	%%ymm9 ,%%ymm11,%%ymm11			\n\t"\
		"vsubpd	%%ymm5,%%ymm7,%%ymm7			\n\t		vsubpd	%%ymm13,%%ymm15,%%ymm15			\n\t"\
	"vmovaps	%%ymm13,(%%rsi) 	\n\t"/* spill ymm13 to make room for 2.0 */"	vmovaps	(%%r8),%%ymm13 	\n\t"/* two */\
	"vfmadd132pd	%%ymm13,%%ymm2,%%ymm0		\n\t	vfmadd132pd	%%ymm13,%%ymm10,%%ymm8 	\n\t"\
	"vfmadd132pd	%%ymm13,%%ymm6,%%ymm4		\n\t	vfmadd132pd	%%ymm13,%%ymm14,%%ymm12	\n\t"\
	"vfmadd132pd	%%ymm13,%%ymm3,%%ymm1		\n\t	vfmadd132pd	%%ymm13,%%ymm11,%%ymm9 	\n\t"\
	"vfmadd132pd	%%ymm13,%%ymm7,%%ymm5		\n\t	vfmadd132pd	(%%rsi),%%ymm15,%%ymm13	\n\t"\
		/* Finish radix-4 butterfly: */\
		"vsubpd	%%ymm4,%%ymm0,%%ymm0			\n\t		vsubpd	%%ymm12,%%ymm8 ,%%ymm8 			\n\t"\
		"vsubpd	%%ymm7,%%ymm2,%%ymm2			\n\t		vsubpd	%%ymm13,%%ymm9 ,%%ymm9 			\n\t"\
		"vsubpd	%%ymm5,%%ymm1,%%ymm1			\n\t		vsubpd	%%ymm15,%%ymm10,%%ymm10			\n\t"\
		"vsubpd	%%ymm6,%%ymm3,%%ymm3			\n\t		vsubpd	%%ymm14,%%ymm11,%%ymm11			\n\t"\
	"vmovaps	%%ymm14,(%%rsi) 	\n\t"/* spill ymm14 to make room for 2.0 */"	vmovaps	(%%r8),%%ymm14 	\n\t"/* two */\
	"vfmadd132pd	%%ymm14,%%ymm0,%%ymm4		\n\t	vfmadd132pd	%%ymm14,%%ymm8 ,%%ymm12	\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm2,%%ymm7		\n\t	vfmadd132pd	%%ymm14,%%ymm9 ,%%ymm13	\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm1,%%ymm5		\n\t	vfmadd132pd	%%ymm14,%%ymm10,%%ymm15	\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm3,%%ymm6		\n\t	vfmadd132pd	(%%rsi),%%ymm11,%%ymm14	\n\t"\
		"vsubpd	%%ymm12,%%ymm4	,%%ymm4			\n\t		vsubpd	%%ymm15,%%ymm11,%%ymm11			\n\t"\
		"vsubpd	%%ymm9 ,%%ymm0	,%%ymm0			\n\t		vsubpd	%%ymm10,%%ymm14,%%ymm14			\n\t"\
	"vmovaps	%%ymm8 ,(%%rsi) 	\n\t"/* spill ymm8  to make room for 2.0 */"	vmovaps	(%%r8),%%ymm8  	\n\t"/* two */\
		"vsubpd	%%ymm13,%%ymm5	,%%ymm5			\n\t	vfmadd132pd	%%ymm8 ,%%ymm11,%%ymm15			\n\t"\
		"vsubpd	(%%rsi),%%ymm1	,%%ymm1			\n\t	vfmadd132pd	%%ymm8 ,%%ymm14,%%ymm10			\n\t"\
	"vfmadd132pd	%%ymm8 ,%%ymm4 ,%%ymm12		\n\t	vfnmadd231pd	-0x900(%%r8),%%ymm15,%%ymm7		\n\t"/* isrt2 */\
	"vfmadd132pd	%%ymm8 ,%%ymm0 ,%%ymm9 		\n\t	vfnmadd231pd	-0x900(%%r8),%%ymm14,%%ymm2		\n\t"\
	"vfmadd132pd	%%ymm8 ,%%ymm5 ,%%ymm13		\n\t	vfnmadd231pd	-0x900(%%r8),%%ymm11,%%ymm3		\n\t"\
	"vfmadd132pd	(%%rsi),%%ymm1 ,%%ymm8 		\n\t	vfnmadd231pd	-0x900(%%r8),%%ymm10,%%ymm6		\n\t"\
		/* SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS(r10,r12,r14,r16,r18,r1A,r1C,r1E) */\
	"vmovaps	%%ymm12,(%%rsi) 	\n\t"/* store ymm12 to make room for sqrt2 */"	vmovaps	0x20(%%r8),%%ymm12 	\n\t"/* sqrt2 */\
		"vmovaps	%%ymm4 ,0x100(%%rsi)		\n\t	vfmadd132pd	%%ymm12,%%ymm7 ,%%ymm15	\n\t"/* sqrt2, use of FMA above means we must combine 2*isrt2 here */\
		"vmovaps	%%ymm0 ,0x180(%%rsi)		\n\t	vfmadd132pd	%%ymm12,%%ymm2 ,%%ymm14	\n\t"\
		"vmovaps	%%ymm5 ,0x120(%%rsi)		\n\t	vfmadd132pd	%%ymm12,%%ymm3 ,%%ymm11	\n\t"\
		"vmovaps	%%ymm1 ,0x0a0(%%rsi)		\n\t	vfmadd132pd	%%ymm12,%%ymm6 ,%%ymm10	\n\t"\
		/*vmovaps	%%ymm12,     (%%rsi)*/	"	\n\t		vmovaps	%%ymm7 ,0x140(%%rsi)	\n\t"\
		"vmovaps	%%ymm9 ,0x080(%%rsi)		\n\t		vmovaps	%%ymm2 ,0x1c0(%%rsi)	\n\t"\
		"vmovaps	%%ymm13,0x020(%%rsi)		\n\t		vmovaps	%%ymm3 ,0x160(%%rsi)	\n\t"\
		"vmovaps	%%ymm8 ,0x1a0(%%rsi)		\n\t		vmovaps	%%ymm6 ,0x0e0(%%rsi)	\n\t"\
		"													vmovaps	%%ymm15,0x040(%%rsi)	\n\t"\
		"													vmovaps	%%ymm14,0x0c0(%%rsi)	\n\t"\
		"													vmovaps	%%ymm11,0x060(%%rsi)	\n\t"\
		"													vmovaps	%%ymm10,0x1e0(%%rsi)	\n\t"\
		"\n\t"\
	/*...Block 3: */\
		"addq	$0x200,%%rsi	/* r20 */		\n\t"\
		"movslq	%[__p10],%%r14					\n\t"\
		"shlq	$3,%%r14						\n\t"\
		"subq	%%rdi,%%r14	/* p10-p8 */		\n\t"\
		"addq	%%r14,%%rax	/* add0+p10 */		\n\t		leaq	(%%r9,%%rax),%%r10	\n\t"\
		"addq	%%r14,%%rbx						\n\t		leaq	(%%r9,%%rbx),%%r11	\n\t"\
		"addq	%%r14,%%rcx						\n\t		leaq	(%%r9,%%rcx),%%r12	\n\t"\
		"addq	%%r14,%%rdx						\n\t		leaq	(%%r9,%%rdx),%%r13	\n\t"\
		/* SSE2_RADIX4_DIT_0TWIDDLE_B(r20) */			/* SSE2_RADIX4_DIT_0TWIDDLE_2NDOFTWO_B(r28) */\
		"vmovaps	     (%%rax),%%ymm2			\n\t		vmovaps	     (%%r10),%%ymm10		\n\t"\
		"vmovaps	     (%%rcx),%%ymm6			\n\t		vmovaps	     (%%r12),%%ymm14		\n\t"\
		"vmovaps	0x020(%%rax),%%ymm3			\n\t		vmovaps	0x020(%%r10),%%ymm11		\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm7			\n\t		vmovaps	0x020(%%r12),%%ymm15		\n\t"\
		"vmovaps	     (%%rbx),%%ymm0			\n\t		vmovaps	     (%%r11),%%ymm8 		\n\t"\
		"vmovaps	     (%%rdx),%%ymm4			\n\t		vmovaps	     (%%r13),%%ymm12		\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm1			\n\t		vmovaps	0x020(%%r11),%%ymm9 		\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm5			\n\t		vmovaps	0x020(%%r13),%%ymm13		\n\t"\
		"vsubpd	%%ymm0,%%ymm2,%%ymm2			\n\t		vsubpd	%%ymm8 ,%%ymm10,%%ymm10			\n\t"\
		"vsubpd	%%ymm4,%%ymm6,%%ymm6			\n\t		vsubpd	%%ymm12,%%ymm14,%%ymm14			\n\t"\
		"vsubpd	%%ymm1,%%ymm3,%%ymm3			\n\t		vsubpd	%%ymm9 ,%%ymm11,%%ymm11			\n\t"\
		"vsubpd	%%ymm5,%%ymm7,%%ymm7			\n\t		vsubpd	%%ymm13,%%ymm15,%%ymm15			\n\t"\
	"vmovaps	%%ymm13,(%%rsi) 	\n\t"/* spill ymm13 to make room for 2.0 */"	vmovaps	(%%r8),%%ymm13 	\n\t"/* two */\
	"vfmadd132pd	%%ymm13,%%ymm2,%%ymm0		\n\t	vfmadd132pd	%%ymm13,%%ymm10,%%ymm8 	\n\t"\
	"vfmadd132pd	%%ymm13,%%ymm6,%%ymm4		\n\t	vfmadd132pd	%%ymm13,%%ymm14,%%ymm12	\n\t"\
	"vfmadd132pd	%%ymm13,%%ymm3,%%ymm1		\n\t	vfmadd132pd	%%ymm13,%%ymm11,%%ymm9 	\n\t"\
	"vfmadd132pd	%%ymm13,%%ymm7,%%ymm5		\n\t	vfmadd132pd	(%%rsi),%%ymm15,%%ymm13	\n\t"\
		/* Finish radix-4 butterfly: */\
		"vsubpd	%%ymm4,%%ymm0,%%ymm0			\n\t		vsubpd	%%ymm12,%%ymm8 ,%%ymm8 			\n\t"\
		"vsubpd	%%ymm7,%%ymm2,%%ymm2			\n\t		vsubpd	%%ymm13,%%ymm9 ,%%ymm9 			\n\t"\
		"vsubpd	%%ymm5,%%ymm1,%%ymm1			\n\t		vsubpd	%%ymm15,%%ymm10,%%ymm10			\n\t"\
		"vsubpd	%%ymm6,%%ymm3,%%ymm3			\n\t		vsubpd	%%ymm14,%%ymm11,%%ymm11			\n\t"\
	"vmovaps	%%ymm14,(%%rsi) 	\n\t"/* spill ymm14 to make room for 2.0 */"	vmovaps	(%%r8),%%ymm14 	\n\t"/* two */\
	"vfmadd132pd	%%ymm14,%%ymm0,%%ymm4		\n\t	vfmadd132pd	%%ymm14,%%ymm8 ,%%ymm12	\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm2,%%ymm7		\n\t	vfmadd132pd	%%ymm14,%%ymm9 ,%%ymm13	\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm1,%%ymm5		\n\t	vfmadd132pd	%%ymm14,%%ymm10,%%ymm15	\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm3,%%ymm6		\n\t	vfmadd132pd	(%%rsi),%%ymm11,%%ymm14	\n\t"\
		"vsubpd	%%ymm12,%%ymm4	,%%ymm4			\n\t		vsubpd	%%ymm15,%%ymm11,%%ymm11			\n\t"\
		"vsubpd	%%ymm9 ,%%ymm0	,%%ymm0			\n\t		vsubpd	%%ymm10,%%ymm14,%%ymm14			\n\t"\
	"vmovaps	%%ymm8 ,(%%rsi) 	\n\t"/* spill ymm8  to make room for 2.0 */"	vmovaps	(%%r8),%%ymm8  	\n\t"/* two */\
		"vsubpd	%%ymm13,%%ymm5	,%%ymm5			\n\t	vfmadd132pd	%%ymm8 ,%%ymm11,%%ymm15			\n\t"\
		"vsubpd	(%%rsi),%%ymm1	,%%ymm1			\n\t	vfmadd132pd	%%ymm8 ,%%ymm14,%%ymm10			\n\t"\
	"vfmadd132pd	%%ymm8 ,%%ymm4 ,%%ymm12		\n\t	vfnmadd231pd	-0x900(%%r8),%%ymm15,%%ymm7		\n\t"/* isrt2 */\
	"vfmadd132pd	%%ymm8 ,%%ymm0 ,%%ymm9 		\n\t	vfnmadd231pd	-0x900(%%r8),%%ymm14,%%ymm2		\n\t"\
	"vfmadd132pd	%%ymm8 ,%%ymm5 ,%%ymm13		\n\t	vfnmadd231pd	-0x900(%%r8),%%ymm11,%%ymm3		\n\t"\
	"vfmadd132pd	(%%rsi),%%ymm1 ,%%ymm8 		\n\t	vfnmadd231pd	-0x900(%%r8),%%ymm10,%%ymm6		\n\t"\
		/* SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS(r20,r22,r24,r26,r28,r2A,r2C,r2E) */\
	"vmovaps	%%ymm12,(%%rsi) 	\n\t"/* store ymm12 to make room for sqrt2 */"	vmovaps	0x20(%%r8),%%ymm12 	\n\t"/* sqrt2 */\
		"vmovaps	%%ymm4 ,0x100(%%rsi)		\n\t	vfmadd132pd	%%ymm12,%%ymm7 ,%%ymm15	\n\t"/* sqrt2, use of FMA above means we must combine 2*isrt2 here */\
		"vmovaps	%%ymm0 ,0x180(%%rsi)		\n\t	vfmadd132pd	%%ymm12,%%ymm2 ,%%ymm14	\n\t"\
		"vmovaps	%%ymm5 ,0x120(%%rsi)		\n\t	vfmadd132pd	%%ymm12,%%ymm3 ,%%ymm11	\n\t"\
		"vmovaps	%%ymm1 ,0x0a0(%%rsi)		\n\t	vfmadd132pd	%%ymm12,%%ymm6 ,%%ymm10	\n\t"\
		/*vmovaps	%%ymm12,     (%%rsi)*/	"	\n\t		vmovaps	%%ymm7 ,0x140(%%rsi)	\n\t"\
		"vmovaps	%%ymm9 ,0x080(%%rsi)		\n\t		vmovaps	%%ymm2 ,0x1c0(%%rsi)	\n\t"\
		"vmovaps	%%ymm13,0x020(%%rsi)		\n\t		vmovaps	%%ymm3 ,0x160(%%rsi)	\n\t"\
		"vmovaps	%%ymm8 ,0x1a0(%%rsi)		\n\t		vmovaps	%%ymm6 ,0x0e0(%%rsi)	\n\t"\
		"													vmovaps	%%ymm15,0x040(%%rsi)	\n\t"\
		"													vmovaps	%%ymm14,0x0c0(%%rsi)	\n\t"\
		"													vmovaps	%%ymm11,0x060(%%rsi)	\n\t"\
		"													vmovaps	%%ymm10,0x1e0(%%rsi)	\n\t"\
		"\n\t"\
	/*...Block 4: */\
		"addq	$0x200,%%rsi	/* r30 */		\n\t"\
		"movslq	%[__p08],%%rdi					\n\t"\
		"shlq	$3,%%rdi						\n\t"\
		"addq	%%rdi,%%rax	/* add0+p18 */		\n\t		leaq	(%%r9,%%rax),%%r10	\n\t"\
		"addq	%%rdi,%%rbx						\n\t		leaq	(%%r9,%%rbx),%%r11	\n\t"\
		"addq	%%rdi,%%rcx						\n\t		leaq	(%%r9,%%rcx),%%r12	\n\t"\
		"addq	%%rdi,%%rdx						\n\t		leaq	(%%r9,%%rdx),%%r13	\n\t"\
		/* SSE2_RADIX4_DIT_0TWIDDLE_B(r30) */		/* SSE2_RADIX4_DIT_0TWIDDLE_2NDOFTWO_B(r38) */\
		"vmovaps	     (%%rax),%%ymm2			\n\t		vmovaps	     (%%r10),%%ymm10		\n\t"\
		"vmovaps	     (%%rcx),%%ymm6			\n\t		vmovaps	     (%%r12),%%ymm14		\n\t"\
		"vmovaps	0x020(%%rax),%%ymm3			\n\t		vmovaps	0x020(%%r10),%%ymm11		\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm7			\n\t		vmovaps	0x020(%%r12),%%ymm15		\n\t"\
		"vmovaps	     (%%rbx),%%ymm0			\n\t		vmovaps	     (%%r11),%%ymm8 		\n\t"\
		"vmovaps	     (%%rdx),%%ymm4			\n\t		vmovaps	     (%%r13),%%ymm12		\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm1			\n\t		vmovaps	0x020(%%r11),%%ymm9 		\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm5			\n\t		vmovaps	0x020(%%r13),%%ymm13		\n\t"\
		"vsubpd	%%ymm0,%%ymm2,%%ymm2			\n\t		vsubpd	%%ymm8 ,%%ymm10,%%ymm10			\n\t"\
		"vsubpd	%%ymm4,%%ymm6,%%ymm6			\n\t		vsubpd	%%ymm12,%%ymm14,%%ymm14			\n\t"\
		"vsubpd	%%ymm1,%%ymm3,%%ymm3			\n\t		vsubpd	%%ymm9 ,%%ymm11,%%ymm11			\n\t"\
		"vsubpd	%%ymm5,%%ymm7,%%ymm7			\n\t		vsubpd	%%ymm13,%%ymm15,%%ymm15			\n\t"\
	"vmovaps	%%ymm13,(%%rsi) 	\n\t"/* spill ymm13 to make room for 2.0 */"	vmovaps	(%%r8),%%ymm13 	\n\t"/* two */\
	"vfmadd132pd	%%ymm13,%%ymm2,%%ymm0		\n\t	vfmadd132pd	%%ymm13,%%ymm10,%%ymm8 	\n\t"\
	"vfmadd132pd	%%ymm13,%%ymm6,%%ymm4		\n\t	vfmadd132pd	%%ymm13,%%ymm14,%%ymm12	\n\t"\
	"vfmadd132pd	%%ymm13,%%ymm3,%%ymm1		\n\t	vfmadd132pd	%%ymm13,%%ymm11,%%ymm9 	\n\t"\
	"vfmadd132pd	%%ymm13,%%ymm7,%%ymm5		\n\t	vfmadd132pd	(%%rsi),%%ymm15,%%ymm13	\n\t"\
		/* Finish radix-4 butterfly: */\
		"vsubpd	%%ymm4,%%ymm0,%%ymm0			\n\t		vsubpd	%%ymm12,%%ymm8 ,%%ymm8 			\n\t"\
		"vsubpd	%%ymm7,%%ymm2,%%ymm2			\n\t		vsubpd	%%ymm13,%%ymm9 ,%%ymm9 			\n\t"\
		"vsubpd	%%ymm5,%%ymm1,%%ymm1			\n\t		vsubpd	%%ymm15,%%ymm10,%%ymm10			\n\t"\
		"vsubpd	%%ymm6,%%ymm3,%%ymm3			\n\t		vsubpd	%%ymm14,%%ymm11,%%ymm11			\n\t"\
	"vmovaps	%%ymm14,(%%rsi) 	\n\t"/* spill ymm14 to make room for 2.0 */"	vmovaps	(%%r8),%%ymm14 	\n\t"/* two */\
	"vfmadd132pd	%%ymm14,%%ymm0,%%ymm4		\n\t	vfmadd132pd	%%ymm14,%%ymm8 ,%%ymm12	\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm2,%%ymm7		\n\t	vfmadd132pd	%%ymm14,%%ymm9 ,%%ymm13	\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm1,%%ymm5		\n\t	vfmadd132pd	%%ymm14,%%ymm10,%%ymm15	\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm3,%%ymm6		\n\t	vfmadd132pd	(%%rsi),%%ymm11,%%ymm14	\n\t"\
		"vsubpd	%%ymm12,%%ymm4	,%%ymm4			\n\t		vsubpd	%%ymm15,%%ymm11,%%ymm11			\n\t"\
		"vsubpd	%%ymm9 ,%%ymm0	,%%ymm0			\n\t		vsubpd	%%ymm10,%%ymm14,%%ymm14			\n\t"\
	"vmovaps	%%ymm8 ,(%%rsi) 	\n\t"/* spill ymm8  to make room for 2.0 */"	vmovaps	(%%r8),%%ymm8  	\n\t"/* two */\
		"vsubpd	%%ymm13,%%ymm5	,%%ymm5			\n\t	vfmadd132pd	%%ymm8 ,%%ymm11,%%ymm15			\n\t"\
		"vsubpd	(%%rsi),%%ymm1	,%%ymm1			\n\t	vfmadd132pd	%%ymm8 ,%%ymm14,%%ymm10			\n\t"\
	"vfmadd132pd	%%ymm8 ,%%ymm4 ,%%ymm12		\n\t	vfnmadd231pd	-0x900(%%r8),%%ymm15,%%ymm7		\n\t"/* isrt2 */\
	"vfmadd132pd	%%ymm8 ,%%ymm0 ,%%ymm9 		\n\t	vfnmadd231pd	-0x900(%%r8),%%ymm14,%%ymm2		\n\t"\
	"vfmadd132pd	%%ymm8 ,%%ymm5 ,%%ymm13		\n\t	vfnmadd231pd	-0x900(%%r8),%%ymm11,%%ymm3		\n\t"\
	"vfmadd132pd	(%%rsi),%%ymm1 ,%%ymm8 		\n\t	vfnmadd231pd	-0x900(%%r8),%%ymm10,%%ymm6		\n\t"\
		/* SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS(r30,r32,r34,r36,r38,r3A,r3C,r3E) */\
	"vmovaps	%%ymm12,(%%rsi) 	\n\t"/* store ymm12 to make room for sqrt2 */"	vmovaps	0x20(%%r8),%%ymm12 	\n\t"/* sqrt2 */\
		"vmovaps	%%ymm4 ,0x100(%%rsi)		\n\t	vfmadd132pd	%%ymm12,%%ymm7 ,%%ymm15	\n\t"/* sqrt2, use of FMA above means we must combine 2*isrt2 here */\
		"vmovaps	%%ymm0 ,0x180(%%rsi)		\n\t	vfmadd132pd	%%ymm12,%%ymm2 ,%%ymm14	\n\t"\
		"vmovaps	%%ymm5 ,0x120(%%rsi)		\n\t	vfmadd132pd	%%ymm12,%%ymm3 ,%%ymm11	\n\t"\
		"vmovaps	%%ymm1 ,0x0a0(%%rsi)		\n\t	vfmadd132pd	%%ymm12,%%ymm6 ,%%ymm10	\n\t"\
		/*vmovaps	%%ymm12,     (%%rsi)*/	"	\n\t		vmovaps	%%ymm7 ,0x140(%%rsi)	\n\t"\
		"vmovaps	%%ymm9 ,0x080(%%rsi)		\n\t		vmovaps	%%ymm2 ,0x1c0(%%rsi)	\n\t"\
		"vmovaps	%%ymm13,0x020(%%rsi)		\n\t		vmovaps	%%ymm3 ,0x160(%%rsi)	\n\t"\
		"vmovaps	%%ymm8 ,0x1a0(%%rsi)		\n\t		vmovaps	%%ymm6 ,0x0e0(%%rsi)	\n\t"\
		"													vmovaps	%%ymm15,0x040(%%rsi)	\n\t"\
		"													vmovaps	%%ymm14,0x0c0(%%rsi)	\n\t"\
		"													vmovaps	%%ymm11,0x060(%%rsi)	\n\t"\
		"													vmovaps	%%ymm10,0x1e0(%%rsi)	\n\t"\
		"\n\t"\
	/**********************************************************************************/\
	/*...and now do eight radix-4 transforms, including the internal twiddle factors: */\
	/**********************************************************************************/\
		"\n\t"\
		"movslq	%[__p10],%%rdi	\n\t"/* tdi will store copy of p10 throughout */\
	/*...Block 1: t00,t10,t20,t30	*/							/*...Block 5: t08,t18,t28,t38*/\
		"movq	%[__add0],%%rax							\n\t		movslq	%[__p04],%%rsi			\n\t"\
		"movslq	%[__p08],%%rbx							\n\t		movq	%%rax,%%r10				\n\t"\
		"movq	%[__r00],%%rcx							\n\t		shlq	$3,%%rsi				\n\t"\
		"movq	%[__r10],%%rdx							\n\t		addq	%%rsi	,%%r10			\n\t"/* add0 = &a[j1+p4] */\
		"shlq	$3,%%rbx								\n\t		movq	%[__isrt2],%%rsi		\n\t"\
		"shlq	$3,%%rdi								\n\t		leaq	(%%rbx,%%r10),%%r11		\n\t"/* add1 = add0+p12; Need this rbx-read before add0 gets added to rbx at left */\
		"addq	%%rax,%%rbx	\n\t"/* add1 = add0+p8 */"	\n\t		vmovaps	(%%rsi),%%ymm10			\n\t"/* isrt2 */\
		"												\n\t		vmovaps	(%%r8) ,%%ymm11			\n\t"/*   two */\
		"vmovaps	     (%%rdx),%%ymm2					\n\t		vmovaps	0x500(%%rcx),%%ymm12		\n\t"\
		"vmovaps	0x400(%%rdx),%%ymm4					\n\t		vmovaps	0x520(%%rcx),%%ymm13		\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm3					\n\t		vmovaps	0x500(%%rdx),%%ymm14		\n\t"\
		"vmovaps	0x420(%%rdx),%%ymm5					\n\t		vmovaps	0x520(%%rdx),%%ymm15		\n\t"\
		"vmovaps	     (%%rcx),%%ymm0					\n\t		vmulpd	%%ymm10,%%ymm13,%%ymm13		\n\t"\
		"vmovaps	0x400(%%rcx),%%ymm6					\n\t		vmulpd	%%ymm10,%%ymm14,%%ymm14		\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm1					\n\t		vmovaps	0x100(%%rcx),%%ymm8 		\n\t"\
		"vmovaps	0x420(%%rcx),%%ymm7					\n\t		vmovaps	0x100(%%rdx),%%ymm9 		\n\t"\
		"vsubpd	%%ymm2,%%ymm0,%%ymm0					\n\t	vfnmadd231pd	%%ymm12,%%ymm10,%%ymm13		\n\t"\
		"vsubpd	%%ymm4,%%ymm6,%%ymm6					\n\t	vfnmadd231pd	%%ymm15,%%ymm10,%%ymm14		\n\t"\
		"vsubpd	%%ymm3,%%ymm1,%%ymm1					\n\t	 vfmadd132pd	0x20(%%r8),%%ymm13,%%ymm12	\n\t"/* sqrt2 */\
		"vsubpd	%%ymm5,%%ymm7,%%ymm7					\n\t	 vfmadd132pd	0x20(%%r8),%%ymm14,%%ymm15	\n\t"\
	"vfmadd132pd	%%ymm11,%%ymm0,%%ymm2				\n\t		vmovaps	0x120(%%rdx),%%ymm10		\n\t"\
	"vfmadd132pd	%%ymm11,%%ymm6,%%ymm4				\n\t		vsubpd	%%ymm10,%%ymm8 ,%%ymm8 		\n\t"\
	"vfmadd132pd	%%ymm11,%%ymm1,%%ymm3				\n\t		vsubpd	%%ymm14,%%ymm12,%%ymm12		\n\t"\
	"vfmadd132pd	%%ymm11,%%ymm7,%%ymm5				\n\t		vmovaps	0x120(%%rcx),%%ymm11		\n\t"\
		/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c00): */"\n\t		vsubpd	%%ymm15,%%ymm13,%%ymm13		\n\t"\
	"vmovaps %%ymm6,(%%rax) \n\t vmovaps (%%r8),%%ymm6  \n\t		vsubpd	%%ymm9 ,%%ymm11,%%ymm11		\n\t"/* lcol: spill ymm6 to make room for 2.0 */\
		"vsubpd	%%ymm4,%%ymm2,%%ymm2					\n\t	vfmadd132pd	%%ymm6 ,%%ymm12,%%ymm14		\n\t"\
		"vsubpd	%%ymm7,%%ymm0,%%ymm0					\n\t	vfmadd132pd	%%ymm6 ,%%ymm8 ,%%ymm10		\n\t"\
		"vsubpd	%%ymm5,%%ymm3,%%ymm3					\n\t	vfmadd132pd	%%ymm6 ,%%ymm13,%%ymm15		\n\t"\
		"vsubpd	(%%rax),%%ymm1,%%ymm1					\n\t	vfmadd132pd	%%ymm6 ,%%ymm11,%%ymm9 		\n\t"\
	"vfmadd132pd	%%ymm6 ,%%ymm2,%%ymm4				\n\t"		/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c04 = c00 + 0x100): */\
	"vfmadd132pd	%%ymm6 ,%%ymm0,%%ymm7				\n\t		vsubpd	%%ymm12,%%ymm10,%%ymm10		\n\t"\
	"vfmadd132pd	%%ymm6 ,%%ymm3,%%ymm5				\n\t		vsubpd	%%ymm15,%%ymm8 ,%%ymm8 		\n\t"\
	"vfmadd132pd	(%%rax),%%ymm1,%%ymm6				\n\t		vsubpd	%%ymm13,%%ymm11,%%ymm11		\n\t"\
		"addq	$0x0e0,%%rsi			\n\t"/* c00 */	"			vsubpd	%%ymm14,%%ymm9 ,%%ymm9 		\n\t"\
		"vmovaps	%%ymm2,     (%%rcx)					\n\t	vfmadd132pd	(%%r8) ,%%ymm10,%%ymm12		\n\t"\
		"vmovaps	%%ymm0,0x020(%%rdx)					\n\t	vfmadd132pd	(%%r8) ,%%ymm8 ,%%ymm15		\n\t"\
		"vmovaps	%%ymm3,0x020(%%rcx)					\n\t	vfmadd132pd	(%%r8) ,%%ymm11,%%ymm13		\n\t"\
		"vmovaps	%%ymm6,     (%%rdx)					\n\t	vfmadd132pd	(%%r8) ,%%ymm9 ,%%ymm14		\n\t"\
		"vmovaps	%%ymm4,%%ymm2		\n\t	vmovaps	%%ymm10,0x100(%%rcx)	\n\t	vmovaps	%%ymm12,%%ymm10	\n\t"\
		"vmovaps	%%ymm7,%%ymm0		\n\t	vmovaps	%%ymm8 ,0x120(%%rdx)	\n\t	vmovaps	%%ymm15,%%ymm8 	\n\t"\
		"vmovaps	%%ymm5,%%ymm3		\n\t	vmovaps	%%ymm11,0x120(%%rcx)	\n\t	vmovaps	%%ymm13,%%ymm11	\n\t"\
		"vmovaps	%%ymm1,%%ymm6		\n\t	vmovaps	%%ymm14,0x100(%%rdx)	\n\t	vmovaps	%%ymm9 ,%%ymm14	\n\t"\
		"vmulpd		     (%%rsi),%%ymm5,%%ymm5			\n\t		vmulpd		0x100(%%rsi),%%ymm13,%%ymm13	\n\t"\
		"vmulpd		0x040(%%rsi),%%ymm1,%%ymm1			\n\t		vmulpd		0x140(%%rsi),%%ymm9 ,%%ymm9 	\n\t"\
		"vmulpd		     (%%rsi),%%ymm4,%%ymm4			\n\t		vmulpd		0x100(%%rsi),%%ymm12,%%ymm12	\n\t"\
		"vmulpd		0x040(%%rsi),%%ymm7,%%ymm7			\n\t		vmulpd		0x140(%%rsi),%%ymm15,%%ymm15	\n\t"\
	"vfnmadd231pd	0x020(%%rsi),%%ymm2,%%ymm5			\n\t	vfnmadd231pd	0x120(%%rsi),%%ymm10,%%ymm13	\n\t"\
	"vfnmadd231pd	0x060(%%rsi),%%ymm0,%%ymm1			\n\t	vfnmadd231pd	0x160(%%rsi),%%ymm8 ,%%ymm9 	\n\t"\
	" vfmadd231pd	0x020(%%rsi),%%ymm3,%%ymm4			\n\t	 vfmadd231pd	0x120(%%rsi),%%ymm11,%%ymm12	\n\t"\
	" vfmadd231pd	0x060(%%rsi),%%ymm6,%%ymm7			\n\t	 vfmadd231pd	0x160(%%rsi),%%ymm14,%%ymm15	\n\t"\
		"vmovaps	%%ymm5,0x020(%%rax)					\n\t		vmovaps	%%ymm13,0x020(%%r10)	\n\t"\
		"vmovaps	%%ymm1,0x020(%%rbx)					\n\t		vmovaps	%%ymm9 ,0x020(%%r11)	\n\t"\
		"vmovaps	%%ymm4,     (%%rax)					\n\t		vmovaps	%%ymm12,     (%%r10)	\n\t"\
		"vmovaps	%%ymm7,     (%%rbx)					\n\t		vmovaps	%%ymm15,     (%%r11)	\n\t"\
		"addq	%%rdi,%%rax	\n\t	addq	%%rdi,%%rbx	\n\t		addq	%%rdi,%%r10			\n\t"\
		"addq	$0x080,%%rsi			\n\t"/* c10 */	"			addq	%%rdi,%%r11			\n\t"\
		"vmovaps	     (%%rcx),%%ymm4					\n\t		vmovaps	0x100(%%rcx),%%ymm12	\n\t"/* c14 = c10 + 0x100: */\
		"vmovaps	0x020(%%rdx),%%ymm0					\n\t		vmovaps	0x120(%%rdx),%%ymm8 	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5					\n\t		vmovaps	0x120(%%rcx),%%ymm13	\n\t"\
		"vmovaps	     (%%rdx),%%ymm6					\n\t		vmovaps	0x100(%%rdx),%%ymm14	\n\t"\
		"vmovaps	%%ymm4,%%ymm2						\n\t		vmovaps	%%ymm12,%%ymm10		\n\t"\
		"vmovaps	%%ymm0,%%ymm1						\n\t		vmovaps	%%ymm8 ,%%ymm9 		\n\t"\
		"vmovaps	%%ymm5,%%ymm3						\n\t		vmovaps	%%ymm13,%%ymm11		\n\t"\
		"vmovaps	%%ymm6,%%ymm7						\n\t		vmovaps	%%ymm14,%%ymm15		\n\t"\
		"vmulpd		     (%%rsi),%%ymm5,%%ymm5			\n\t		vmulpd		0x100(%%rsi),%%ymm13,%%ymm13	\n\t"\
		"vmulpd		0x040(%%rsi),%%ymm6,%%ymm6			\n\t		vmulpd		0x140(%%rsi),%%ymm14,%%ymm14	\n\t"\
		"vmulpd		     (%%rsi),%%ymm4,%%ymm4			\n\t		vmulpd		0x100(%%rsi),%%ymm12,%%ymm12	\n\t"\
		"vmulpd		0x040(%%rsi),%%ymm0,%%ymm0			\n\t		vmulpd		0x140(%%rsi),%%ymm8 ,%%ymm8 	\n\t"\
	"vfnmadd231pd	0x020(%%rsi),%%ymm2,%%ymm5			\n\t	vfnmadd231pd	0x120(%%rsi),%%ymm10,%%ymm13	\n\t"\
	"vfnmadd231pd	0x060(%%rsi),%%ymm1,%%ymm6			\n\t	vfnmadd231pd	0x160(%%rsi),%%ymm9 ,%%ymm14	\n\t"\
	" vfmadd231pd	0x020(%%rsi),%%ymm3,%%ymm4			\n\t	 vfmadd231pd	0x120(%%rsi),%%ymm11,%%ymm12	\n\t"\
	" vfmadd231pd	0x060(%%rsi),%%ymm7,%%ymm0			\n\t	 vfmadd231pd	0x160(%%rsi),%%ymm15,%%ymm8 	\n\t"\
		"vmovaps	%%ymm5,0x020(%%rax)					\n\t		vmovaps	%%ymm13,0x020(%%r10)	\n\t"\
		"vmovaps	%%ymm6,0x020(%%rbx)					\n\t		vmovaps	%%ymm14,0x020(%%r11)	\n\t"\
		"vmovaps	%%ymm4,     (%%rax)					\n\t		vmovaps	%%ymm12,     (%%r10)	\n\t"\
		"vmovaps	%%ymm0,     (%%rbx)					\n\t		vmovaps	%%ymm8 ,     (%%r11)	\n\t"\
		"\n\t"\
	/*...Block 2: t02,t12,t22,t32	*/							/*...Block 6: t0A,t1A,t2A,t3A*/\
		"movq	%[__add0],%%rax							\n\t		movq	%%rax,%%r10			\n\t"\
		"movslq	%[__p01],%%rsi							\n\t		movslq	%[__p05],%%r9		\n\t"\
		"movslq	%[__p08],%%rbx							\n\t		shlq	$3,%%r9			\n\t"\
		"shlq	$3,%%rsi								\n\t		addq	%%r9 	,%%r10	/* add0 = &a[j1+p5] */\n\t"\
		"shlq	$3,%%rbx								\n\t		movq	%%r10,%%r11			\n\t"\
		"addq	%%rsi	,%%rax	/* add0 = &a[j1+p1] */	\n\t		addq	%%rbx	,%%r11	/* add1 = add0+p13 */\n\t"\
		"addq	%%rax	,%%rbx	/* add1 = add0+p8 */	\n\t"\
		"movq	%[__isrt2],%%rsi						\n\t"\
		"addq	$0x440,%%rcx	/* r22 */				\n\t"\
		"addq	$0x440,%%rdx	/* r32 */				\n\t"\
		"addq	$0x060,%%rsi	/* cc1; cc3 += 0x040: */\n\t"\
		"vmovaps	     (%%rcx),%%ymm4					\n\t		vmovaps	0x100(%%rcx),%%ymm12	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5					\n\t		vmovaps	0x120(%%rcx),%%ymm13	\n\t"\
		"vmovaps	     (%%rsi),%%ymm2		\n\t"/* cc1 */"			vmovaps	0x040(%%rsi),%%ymm11	\n\t"/* cc3 */\
		"vmovaps	0x020(%%rsi),%%ymm3					\n\t		vmovaps	0x060(%%rsi),%%ymm10	\n\t"\
		"vmovaps	%%ymm4,%%ymm6						\n\t		vmovaps	%%ymm12,%%ymm14		\n\t"\
		"vmovaps	%%ymm5,%%ymm7						\n\t		vmovaps	%%ymm13,%%ymm15		\n\t"\
		"vmulpd		%%ymm2,%%ymm4,%%ymm4				\n\t		vmulpd		%%ymm10,%%ymm12,%%ymm12		\n\t"\
		"vmulpd		%%ymm2,%%ymm5,%%ymm5				\n\t		vmulpd		%%ymm10,%%ymm13,%%ymm13		\n\t"\
		"vmovaps	     (%%rdx),%%ymm0					\n\t		vmovaps	0x100(%%rdx),%%ymm8 	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm1					\n\t		vmovaps	0x120(%%rdx),%%ymm9 	/* t3B */\n\t"\
	" vfmadd231pd	%%ymm3,%%ymm7,%%ymm4				\n\t	 vfmadd231pd	%%ymm11,%%ymm15,%%ymm12		\n\t"\
	"vfnmadd231pd	%%ymm3,%%ymm6,%%ymm5				\n\t	vfnmadd231pd	%%ymm11,%%ymm14,%%ymm13		\n\t"\
		"vmovaps	%%ymm0,%%ymm6						\n\t		vmovaps	%%ymm8 ,%%ymm14		\n\t"\
		"vmovaps	%%ymm1,%%ymm7						\n\t		vmovaps	%%ymm9 ,%%ymm15		\n\t"\
		"vmulpd		%%ymm11,%%ymm0,%%ymm0				\n\t		vmulpd		%%ymm2 ,%%ymm8 ,%%ymm8 		\n\t"\
		"vmulpd		%%ymm11,%%ymm1,%%ymm1				\n\t		vmulpd		%%ymm2 ,%%ymm9 ,%%ymm9 		\n\t"\
	" vfmadd231pd	%%ymm10,%%ymm7,%%ymm0				\n\t	vfnmadd231pd	%%ymm3 ,%%ymm15,%%ymm8 		\n\t"\
	"vfnmadd231pd	%%ymm10,%%ymm6,%%ymm1				\n\t	 vfmadd231pd	%%ymm3 ,%%ymm14,%%ymm9 		\n\t"\
		"vmovaps	%%ymm5,%%ymm7						\n\t		vmovaps	%%ymm13,%%ymm15		\n\t"\
		"vmovaps	%%ymm4,%%ymm6						\n\t		vmovaps	%%ymm12,%%ymm14		\n\t"\
		"vaddpd	%%ymm0,%%ymm4,%%ymm4					\n\t		vaddpd	%%ymm8 ,%%ymm14,%%ymm14		\n\t"\
		"vaddpd	%%ymm1,%%ymm5,%%ymm5					\n\t		vaddpd	%%ymm9 ,%%ymm15,%%ymm15		\n\t"\
		"vsubpd	%%ymm0,%%ymm6,%%ymm6					\n\t		vsubpd	%%ymm8 ,%%ymm12,%%ymm12		\n\t"\
		"vsubpd	%%ymm1,%%ymm7,%%ymm7					\n\t		vsubpd	%%ymm9 ,%%ymm13,%%ymm13		\n\t"\
		"subq	$0x400,%%rcx	/* r02 */				\n\t"\
		"subq	$0x400,%%rdx	/* r12 */				\n\t"\
		"subq	$0x040,%%rsi	/* cc0 */				\n\t"\
		"vmovaps	     (%%rdx),%%ymm1					\n\t		vmovaps	0x100(%%rdx),%%ymm8 	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm3					\n\t		vmovaps	0x120(%%rdx),%%ymm10	\n\t"\
		"vmovaps	0x020(%%rsi),%%ymm2					\n\t		vmovaps	0x020(%%rsi),%%ymm9 	\n\t"\
		"vmovaps	%%ymm1,%%ymm0						\n\t		vmovaps	%%ymm8 ,%%ymm11		\n\t"\
		"vmulpd		%%ymm2 ,%%ymm1,%%ymm1				\n\t		vmulpd		%%ymm9 ,%%ymm8 ,%%ymm8 		\n\t"\
		"vmulpd		%%ymm3 ,%%ymm2,%%ymm2				\n\t		vmulpd		%%ymm10,%%ymm9 ,%%ymm9 		\n\t"\
	"vmovaps %%ymm9,(%%rax) \n\t vmovaps (%%rsi),%%ymm9	\n\t"/* isrt2 */\
	" vfmsub132pd	%%ymm9 ,%%ymm1,%%ymm3				\n\t	vfnmadd231pd	%%ymm9 ,%%ymm10,%%ymm8 		\n\t"\
	" vfmadd231pd	%%ymm9 ,%%ymm0,%%ymm2				\n\t	 vfmadd213pd	(%%rax),%%ymm11,%%ymm9 		\n\t"\
		"vmovaps	     (%%rcx),%%ymm0					\n\t		vmovaps	0x100(%%rcx),%%ymm10	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm1					\n\t		vmovaps	0x120(%%rcx),%%ymm11	\n\t"\
		"vsubpd	%%ymm2,%%ymm0,%%ymm0					\n\t		vsubpd	%%ymm8 ,%%ymm10,%%ymm10		\n\t"\
		"vsubpd	%%ymm3,%%ymm1,%%ymm1					\n\t		vsubpd	%%ymm9 ,%%ymm11,%%ymm11		\n\t"\
	"vfmadd132pd	(%%r8),%%ymm0,%%ymm2				\n\t	vfmadd132pd	(%%r8) ,%%ymm10,%%ymm8 		\n\t"\
	"vfmadd132pd	(%%r8),%%ymm1,%%ymm3				\n\t	vfmadd132pd	(%%r8) ,%%ymm11,%%ymm9 		\n\t"\
		"\n\t"\
		/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c01): */			/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c05): */\
		"addq	$0x4c0,%%rsi	\n\t"/* c01; c05 += 0x100: */\
		"vsubpd	%%ymm4,%%ymm2,%%ymm2					\n\t		vsubpd	%%ymm12,%%ymm10,%%ymm10		\n\t"\
		"vsubpd	%%ymm7,%%ymm0,%%ymm0					\n\t		vsubpd	%%ymm15,%%ymm8 ,%%ymm8 		\n\t"\
		"vsubpd	%%ymm5,%%ymm3,%%ymm3					\n\t		vsubpd	%%ymm13,%%ymm11,%%ymm11		\n\t"\
		"vsubpd	%%ymm6,%%ymm1,%%ymm1					\n\t		vsubpd	%%ymm14,%%ymm9 ,%%ymm9 		\n\t"\
	"vmovaps %%ymm14,(%%rax) \n\t vmovaps (%%r8),%%ymm14\n\t"/* 2.0 */\
	"vfmadd132pd	%%ymm14,%%ymm2,%%ymm4				\n\t	vfmadd132pd	%%ymm14,%%ymm10,%%ymm12		\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm0,%%ymm7				\n\t	vfmadd132pd	%%ymm14,%%ymm8 ,%%ymm15		\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm3,%%ymm5				\n\t	vfmadd132pd	%%ymm14,%%ymm11,%%ymm13		\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm1,%%ymm6				\n\t	vfmadd132pd	(%%rax),%%ymm9 ,%%ymm14		\n\t"\
		"vmovaps	%%ymm2,     (%%rcx)					\n\t		vmovaps	%%ymm10,0x100(%%rcx)	\n\t"\
		"vmovaps	%%ymm0,0x020(%%rdx)					\n\t		vmovaps	%%ymm8 ,0x120(%%rdx)	\n\t"\
		"vmovaps	%%ymm3,0x020(%%rcx)					\n\t		vmovaps	%%ymm11,0x120(%%rcx)	\n\t"\
		"vmovaps	%%ymm6,     (%%rdx)					\n\t		vmovaps	%%ymm14,0x100(%%rdx)	\n\t"\
		"vmovaps	%%ymm4,%%ymm2						\n\t		vmovaps	%%ymm12,%%ymm10		\n\t"\
		"vmovaps	%%ymm7,%%ymm0						\n\t		vmovaps	%%ymm15,%%ymm8 		\n\t"\
		"vmovaps	%%ymm5,%%ymm3						\n\t		vmovaps	%%ymm13,%%ymm11		\n\t"\
		"vmovaps	%%ymm1,%%ymm6						\n\t		vmovaps	%%ymm9 ,%%ymm14		\n\t"\
		"vmulpd		     (%%rsi),%%ymm4,%%ymm4			\n\t		vmulpd		0x100(%%rsi),%%ymm12,%%ymm12	\n\t"\
		"vmulpd		0x040(%%rsi),%%ymm7,%%ymm7			\n\t		vmulpd		0x140(%%rsi),%%ymm15,%%ymm15	\n\t"\
		"vmulpd		     (%%rsi),%%ymm5,%%ymm5			\n\t		vmulpd		0x100(%%rsi),%%ymm13,%%ymm13	\n\t"\
		"vmulpd		0x040(%%rsi),%%ymm1,%%ymm1			\n\t		vmulpd		0x140(%%rsi),%%ymm9 ,%%ymm9 	\n\t"\
	" vfmadd231pd	0x020(%%rsi),%%ymm3,%%ymm4			\n\t	 vfmadd231pd	0x120(%%rsi),%%ymm11,%%ymm12	\n\t"\
	" vfmadd231pd	0x060(%%rsi),%%ymm6,%%ymm7			\n\t	 vfmadd231pd	0x160(%%rsi),%%ymm14,%%ymm15	\n\t"\
	"vfnmadd231pd	0x020(%%rsi),%%ymm2,%%ymm5			\n\t	vfnmadd231pd	0x120(%%rsi),%%ymm10,%%ymm13	\n\t"\
	"vfnmadd231pd	0x060(%%rsi),%%ymm0,%%ymm1			\n\t	vfnmadd231pd	0x160(%%rsi),%%ymm8 ,%%ymm9 	\n\t"\
		"vmovaps	%%ymm5,0x020(%%rax)					\n\t		vmovaps	%%ymm13,0x020(%%r10)	\n\t"\
		"vmovaps	%%ymm1,0x020(%%rbx)					\n\t		vmovaps	%%ymm9 ,0x020(%%r11)	\n\t"\
		"vmovaps	%%ymm4,     (%%rax)					\n\t		vmovaps	%%ymm12,     (%%r10)	\n\t"\
		"vmovaps	%%ymm7,     (%%rbx)					\n\t		vmovaps	%%ymm15,     (%%r11)	\n\t"\
		"addq	%%rdi,%%rax								\n\t		addq	%%rdi,%%r10			\n\t"\
		"addq	%%rdi,%%rbx								\n\t		addq	%%rdi,%%r11			\n\t"\
		"addq	$0x080,%%rsi	\n\t"/* c11; c15 += 0x100: */\
		"vmovaps	     (%%rcx),%%ymm4					\n\t		vmovaps	0x100(%%rcx),%%ymm12	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm0					\n\t		vmovaps	0x120(%%rdx),%%ymm8 	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5					\n\t		vmovaps	0x120(%%rcx),%%ymm13	\n\t"\
		"vmovaps	     (%%rdx),%%ymm6					\n\t		vmovaps	0x100(%%rdx),%%ymm14	\n\t"\
		"vmovaps	%%ymm4,%%ymm2						\n\t		vmovaps	%%ymm12,%%ymm10		\n\t"\
		"vmovaps	%%ymm0,%%ymm1						\n\t		vmovaps	%%ymm8 ,%%ymm9 		\n\t"\
		"vmovaps	%%ymm5,%%ymm3						\n\t		vmovaps	%%ymm13,%%ymm11		\n\t"\
		"vmovaps	%%ymm6,%%ymm7						\n\t		vmovaps	%%ymm14,%%ymm15		\n\t"\
		"vmulpd		     (%%rsi),%%ymm4,%%ymm4			\n\t		vmulpd		0x100(%%rsi),%%ymm12,%%ymm12	\n\t"\
		"vmulpd		0x040(%%rsi),%%ymm0,%%ymm0			\n\t		vmulpd		0x140(%%rsi),%%ymm8 ,%%ymm8 	\n\t"\
		"vmulpd		     (%%rsi),%%ymm5,%%ymm5			\n\t		vmulpd		0x100(%%rsi),%%ymm13,%%ymm13	\n\t"\
		"vmulpd		0x040(%%rsi),%%ymm6,%%ymm6			\n\t		vmulpd		0x140(%%rsi),%%ymm14,%%ymm14	\n\t"\
	" vfmadd231pd	0x020(%%rsi),%%ymm3,%%ymm4			\n\t	 vfmadd231pd	0x120(%%rsi),%%ymm11,%%ymm12	\n\t"\
	" vfmadd231pd	0x060(%%rsi),%%ymm7,%%ymm0			\n\t	 vfmadd231pd	0x160(%%rsi),%%ymm15,%%ymm8 	\n\t"\
	"vfnmadd231pd	0x020(%%rsi),%%ymm2,%%ymm5			\n\t	vfnmadd231pd	0x120(%%rsi),%%ymm10,%%ymm13	\n\t"\
	"vfnmadd231pd	0x060(%%rsi),%%ymm1,%%ymm6			\n\t	vfnmadd231pd	0x160(%%rsi),%%ymm9 ,%%ymm14	\n\t"\
		"vmovaps	%%ymm5,0x020(%%rax)					\n\t		vmovaps	%%ymm13,0x020(%%r10)	\n\t"\
		"vmovaps	%%ymm6,0x020(%%rbx)					\n\t		vmovaps	%%ymm14,0x020(%%r11)	\n\t"\
		"vmovaps	%%ymm4,     (%%rax)					\n\t		vmovaps	%%ymm12,     (%%r10)	\n\t"\
		"vmovaps	%%ymm0,     (%%rbx)					\n\t		vmovaps	%%ymm8 ,     (%%r11)	\n\t"\
		"\n\t"\
	/*...Block 3: t04,t14,t24,t34*/								/*...Block 7: t0C,t1C,t2C,t3C*/\
		"movq	%[__add0],%%rax							\n\t		movq	%%rax,%%r10			\n\t"\
		"movslq	%[__p02],%%rsi							\n\t		movslq	%[__p06],%%r9		\n\t"\
		"movslq	%[__p08],%%rbx							\n\t		shlq	$3,%%r9			\n\t"\
		"shlq	$3,%%rsi								\n\t		addq	%%r9 	,%%r10	/* add0 = &a[j1+p6] */\n\t"\
		"shlq	$3,%%rbx								\n\t		movq	%%r10,%%r11			\n\t"\
		"addq	%%rsi	,%%rax	/* add0 = &a[j1+p2] */	\n\t		addq	%%rbx	,%%r11	/* add1 = add0+p14 */\n\t"\
		"addq	%%rax	,%%rbx	/* add1 = add0+p8 */	\n\t"\
		"movq	%[__isrt2],%%rsi						\n\t"\
		"addq	$0x440,%%rcx	/* r24 */				\n\t"\
		"addq	$0x440,%%rdx	/* r34 */				\n\t"\
		"addq	$0x020,%%rsi	/* cc0 */				\n\t"\
		"vmovaps	     (%%rcx),%%ymm4					\n\t		vmovaps	0x100(%%rcx),%%ymm12	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5					\n\t		vmovaps	0x120(%%rcx),%%ymm13	\n\t"\
		"vmovaps	     (%%rsi),%%ymm2					\n\t"\
		"vmovaps	0x020(%%rsi),%%ymm3					\n\t"\
		"vmovaps	%%ymm4,%%ymm6						\n\t		vmovaps	%%ymm12,%%ymm14		\n\t"\
		"vmovaps	%%ymm5,%%ymm7						\n\t		vmovaps	%%ymm13,%%ymm15		\n\t"\
		"vmulpd	%%ymm2,%%ymm4,%%ymm4					\n\t		vmulpd	%%ymm3,%%ymm12,%%ymm12	\n\t"\
		"vmulpd	%%ymm2,%%ymm5,%%ymm5					\n\t		vmulpd	%%ymm3,%%ymm13,%%ymm13	\n\t"\
		"vmovaps	     (%%rdx),%%ymm0					\n\t		vmovaps	0x100(%%rdx),%%ymm8 	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm1					\n\t		vmovaps	0x120(%%rdx),%%ymm9 	\n\t"\
	" vfmadd231pd	%%ymm3,%%ymm7,%%ymm4				\n\t	 vfmadd231pd	%%ymm2,%%ymm15,%%ymm12		\n\t"\
	"vfnmadd231pd	%%ymm3,%%ymm6,%%ymm5				\n\t	vfnmadd231pd	%%ymm2,%%ymm14,%%ymm13		\n\t"\
		"vmovaps	%%ymm0,%%ymm6						\n\t		vmovaps	%%ymm8,%%ymm14		\n\t"\
		"vmovaps	%%ymm1,%%ymm7						\n\t		vmovaps	%%ymm9,%%ymm15		\n\t"\
		"vmulpd	%%ymm3,%%ymm0,%%ymm0					\n\t		vmulpd	%%ymm2 ,%%ymm8 ,%%ymm8 	\n\t"\
		"vmulpd	%%ymm3,%%ymm1,%%ymm1					\n\t		vmulpd	%%ymm2 ,%%ymm9 ,%%ymm9 	\n\t"\
	" vfmadd231pd	%%ymm2,%%ymm7,%%ymm0				\n\t	 vfmadd231pd	%%ymm3,%%ymm15,%%ymm8 		\n\t"\
	"vfnmadd231pd	%%ymm2,%%ymm6,%%ymm1				\n\t	vfnmadd231pd	%%ymm3,%%ymm14,%%ymm9 		\n\t"\
		"vmovaps	%%ymm5,%%ymm7						\n\t		vmovaps	%%ymm13,%%ymm15		\n\t"\
		"vmovaps	%%ymm4,%%ymm6						\n\t		vmovaps	%%ymm12,%%ymm14		\n\t"\
		"vaddpd	%%ymm0,%%ymm4,%%ymm4					\n\t		vaddpd	%%ymm8,%%ymm14,%%ymm14		\n\t"\
		"vaddpd	%%ymm1,%%ymm5,%%ymm5					\n\t		vaddpd	%%ymm9,%%ymm15,%%ymm15		\n\t"\
		"vsubpd	%%ymm0,%%ymm6,%%ymm6					\n\t		vsubpd	%%ymm8,%%ymm12,%%ymm12		\n\t"\
		"vsubpd	%%ymm1,%%ymm7,%%ymm7					\n\t		vsubpd	%%ymm9,%%ymm13,%%ymm13		\n\t"\
		"subq	$0x400,%%rcx							\n\t"\
		"subq	$0x400,%%rdx							\n\t"\
		"subq	$0x020,%%rsi					\n\t"/* isrt2 */\
		"vmovaps	     (%%rdx),%%ymm2					\n\t		vmovaps	0x100(%%rdx),%%ymm8 	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm3					\n\t		vmovaps	0x120(%%rdx),%%ymm9 	\n\t"\
		"vmovaps	(%%rsi),%%ymm1				\n\t"/* isrt2 */\
		"vmovaps	%%ymm3,%%ymm0						\n\t		vmovaps	%%ymm8 ,%%ymm10		\n\t"\
		"vsubpd	%%ymm2,%%ymm3,%%ymm3					\n\t		vsubpd	%%ymm9 ,%%ymm8 	,%%ymm8 	\n\t"\
		"vaddpd	%%ymm0,%%ymm2,%%ymm2					\n\t		vaddpd	%%ymm10,%%ymm9 	,%%ymm9 	\n\t"\
		"vmulpd	%%ymm1,%%ymm2,%%ymm2					\n\t		vmulpd	%%ymm1 ,%%ymm8 	,%%ymm8 \n\t"\
		"vmulpd	%%ymm1,%%ymm3,%%ymm3					\n\t		vmulpd	%%ymm1 ,%%ymm9 	,%%ymm9 \n\t"\
		"vmovaps	     (%%rcx),%%ymm0					\n\t		vmovaps	0x100(%%rcx),%%ymm10	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm1					\n\t		vmovaps	0x120(%%rcx),%%ymm11	\n\t"\
		"vsubpd	%%ymm2,%%ymm0,%%ymm0					\n\t		vsubpd	%%ymm8 ,%%ymm10,%%ymm10		\n\t"\
		"vsubpd	%%ymm3,%%ymm1,%%ymm1					\n\t		vsubpd	%%ymm9 ,%%ymm11,%%ymm11		\n\t"\
	"vfmadd132pd	(%%r8),%%ymm0,%%ymm2				\n\t	vfmadd132pd	(%%r8) ,%%ymm10,%%ymm8 		\n\t"\
	"vfmadd132pd	(%%r8),%%ymm1,%%ymm3				\n\t	vfmadd132pd	(%%r8) ,%%ymm11,%%ymm9 		\n\t"\
		"\n\t"\
		/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c02): */			/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c06): */\
		"addq	$0x2e0,%%rsi	\n\t"/* c02; c06 += 0x100: */\
		"vsubpd	%%ymm4,%%ymm2,%%ymm2					\n\t		vsubpd	%%ymm12,%%ymm10,%%ymm10		\n\t"\
		"vsubpd	%%ymm7,%%ymm0,%%ymm0					\n\t		vsubpd	%%ymm15,%%ymm8 ,%%ymm8 		\n\t"\
		"vsubpd	%%ymm5,%%ymm3,%%ymm3					\n\t		vsubpd	%%ymm13,%%ymm11,%%ymm11		\n\t"\
		"vsubpd	%%ymm6,%%ymm1,%%ymm1					\n\t		vsubpd	%%ymm14,%%ymm9 ,%%ymm9 		\n\t"\
	"vmovaps %%ymm14,(%%rax) \n\t vmovaps (%%r8),%%ymm14\n\t"/* 2.0 */\
	"vfmadd132pd	%%ymm14,%%ymm2,%%ymm4				\n\t	vfmadd132pd	%%ymm14,%%ymm10,%%ymm12		\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm0,%%ymm7				\n\t	vfmadd132pd	%%ymm14,%%ymm8 ,%%ymm15		\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm3,%%ymm5				\n\t	vfmadd132pd	%%ymm14,%%ymm11,%%ymm13		\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm1,%%ymm6				\n\t	vfmadd132pd	(%%rax),%%ymm9 ,%%ymm14		\n\t"\
		"vmovaps	%%ymm2,     (%%rcx)					\n\t		vmovaps	%%ymm10,0x100(%%rcx)	\n\t"\
		"vmovaps	%%ymm0,0x020(%%rdx)					\n\t		vmovaps	%%ymm8 ,0x120(%%rdx)	\n\t"\
		"vmovaps	%%ymm3,0x020(%%rcx)					\n\t		vmovaps	%%ymm11,0x120(%%rcx)	\n\t"\
		"vmovaps	%%ymm6,     (%%rdx)					\n\t		vmovaps	%%ymm14,0x100(%%rdx)	\n\t"\
		"vmovaps	%%ymm4,%%ymm2						\n\t		vmovaps	%%ymm12,%%ymm10		\n\t"\
		"vmovaps	%%ymm7,%%ymm0						\n\t		vmovaps	%%ymm15,%%ymm8 		\n\t"\
		"vmovaps	%%ymm5,%%ymm3						\n\t		vmovaps	%%ymm13,%%ymm11		\n\t"\
		"vmovaps	%%ymm1,%%ymm6						\n\t		vmovaps	%%ymm9 ,%%ymm14		\n\t"\
		"vmulpd		     (%%rsi),%%ymm4,%%ymm4			\n\t		vmulpd		0x100(%%rsi),%%ymm12,%%ymm12	\n\t"\
		"vmulpd		0x040(%%rsi),%%ymm7,%%ymm7			\n\t		vmulpd		0x140(%%rsi),%%ymm15,%%ymm15	\n\t"\
		"vmulpd		     (%%rsi),%%ymm5,%%ymm5			\n\t		vmulpd		0x100(%%rsi),%%ymm13,%%ymm13	\n\t"\
		"vmulpd		0x040(%%rsi),%%ymm1,%%ymm1			\n\t		vmulpd		0x140(%%rsi),%%ymm9 ,%%ymm9 	\n\t"\
	" vfmadd231pd	0x020(%%rsi),%%ymm3,%%ymm4			\n\t	 vfmadd231pd	0x120(%%rsi),%%ymm11,%%ymm12	\n\t"\
	" vfmadd231pd	0x060(%%rsi),%%ymm6,%%ymm7			\n\t	 vfmadd231pd	0x160(%%rsi),%%ymm14,%%ymm15	\n\t"\
	"vfnmadd231pd	0x020(%%rsi),%%ymm2,%%ymm5			\n\t	vfnmadd231pd	0x120(%%rsi),%%ymm10,%%ymm13	\n\t"\
	"vfnmadd231pd	0x060(%%rsi),%%ymm0,%%ymm1			\n\t	vfnmadd231pd	0x160(%%rsi),%%ymm8 ,%%ymm9 	\n\t"\
		"vmovaps	%%ymm5,0x020(%%rax)					\n\t		vmovaps	%%ymm13,0x020(%%r10)	\n\t"\
		"vmovaps	%%ymm1,0x020(%%rbx)					\n\t		vmovaps	%%ymm9 ,0x020(%%r11)	\n\t"\
		"vmovaps	%%ymm4,     (%%rax)					\n\t		vmovaps	%%ymm12,     (%%r10)	\n\t"\
		"vmovaps	%%ymm7,     (%%rbx)					\n\t		vmovaps	%%ymm15,     (%%r11)	\n\t"\
		"addq	%%rdi,%%rax								\n\t		addq	%%rdi,%%r10			\n\t"\
		"addq	%%rdi,%%rbx								\n\t		addq	%%rdi,%%r11			\n\t"\
		"addq	$0x080,%%rsi	\n\t"/* c12; c16 += 0x100: */\
		"vmovaps	     (%%rcx),%%ymm4					\n\t		vmovaps	0x100(%%rcx),%%ymm12	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm0					\n\t		vmovaps	0x120(%%rdx),%%ymm8 	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5					\n\t		vmovaps	0x120(%%rcx),%%ymm13	\n\t"\
		"vmovaps	     (%%rdx),%%ymm6					\n\t		vmovaps	0x100(%%rdx),%%ymm14	\n\t"\
		"vmovaps	%%ymm4,%%ymm2						\n\t		vmovaps	%%ymm12,%%ymm10		\n\t"\
		"vmovaps	%%ymm0,%%ymm1						\n\t		vmovaps	%%ymm8 ,%%ymm9 		\n\t"\
		"vmovaps	%%ymm5,%%ymm3						\n\t		vmovaps	%%ymm13,%%ymm11		\n\t"\
		"vmovaps	%%ymm6,%%ymm7						\n\t		vmovaps	%%ymm14,%%ymm15		\n\t"\
		"vmulpd		     (%%rsi),%%ymm4,%%ymm4			\n\t		vmulpd		0x100(%%rsi),%%ymm12,%%ymm12	\n\t"\
		"vmulpd		0x040(%%rsi),%%ymm0,%%ymm0			\n\t		vmulpd		0x140(%%rsi),%%ymm8 ,%%ymm8 	\n\t"\
		"vmulpd		     (%%rsi),%%ymm5,%%ymm5			\n\t		vmulpd		0x100(%%rsi),%%ymm13,%%ymm13	\n\t"\
		"vmulpd		0x040(%%rsi),%%ymm6,%%ymm6			\n\t		vmulpd		0x140(%%rsi),%%ymm14,%%ymm14	\n\t"\
	" vfmadd231pd	0x020(%%rsi),%%ymm3,%%ymm4			\n\t	 vfmadd231pd	0x120(%%rsi),%%ymm11,%%ymm12	\n\t"\
	" vfmadd231pd	0x060(%%rsi),%%ymm7,%%ymm0			\n\t	 vfmadd231pd	0x160(%%rsi),%%ymm15,%%ymm8 	\n\t"\
	"vfnmadd231pd	0x020(%%rsi),%%ymm2,%%ymm5			\n\t	vfnmadd231pd	0x120(%%rsi),%%ymm10,%%ymm13	\n\t"\
	"vfnmadd231pd	0x060(%%rsi),%%ymm1,%%ymm6			\n\t	vfnmadd231pd	0x160(%%rsi),%%ymm9 ,%%ymm14	\n\t"\
		"vmovaps	%%ymm5,0x020(%%rax)					\n\t		vmovaps	%%ymm13,0x020(%%r10)	\n\t"\
		"vmovaps	%%ymm6,0x020(%%rbx)					\n\t		vmovaps	%%ymm14,0x020(%%r11)	\n\t"\
		"vmovaps	%%ymm4,     (%%rax)					\n\t		vmovaps	%%ymm12,     (%%r10)	\n\t"\
		"vmovaps	%%ymm0,     (%%rbx)					\n\t		vmovaps	%%ymm8 ,     (%%r11)	\n\t"\
		"\n\t"\
	/*...Block 4: t06,t16,t26,t36*/								/*...Block 8: t0E,t1E,t2E,t3E*/\
		"movq	%[__add0],%%rax							\n\t		movq	%%rax,%%r10			\n\t"\
		"movslq	%[__p03],%%rsi							\n\t		movslq	%[__p07],%%r9		\n\t"\
		"movslq	%[__p08],%%rbx							\n\t		shlq	$3,%%r9			\n\t"\
		"shlq	$3,%%rsi								\n\t		addq	%%r9 	,%%r10	/* add0 = &a[j1+p6] */\n\t"\
		"shlq	$3,%%rbx								\n\t		movq	%%r10,%%r11			\n\t"\
		"addq	%%rsi	,%%rax	/* add0 = &a[j1+p3] */	\n\t		addq	%%rbx	,%%r11	/* add1 = add0+p14 */\n\t"\
		"addq	%%rax	,%%rbx	/* add1 = add0+p8 */	\n\t"\
		"movq	%[__isrt2],%%rsi						\n\t"\
		"addq	$0x440,%%rcx	/* r26 */				\n\t"\
		"addq	$0x440,%%rdx	/* r36 */				\n\t"\
		"addq	$0x060,%%rsi	/* cc1; cc3 += 0x040: */\n\t"\
		"vmovaps	     (%%rcx),%%ymm4					\n\t		vmovaps	0x100(%%rcx),%%ymm12	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5					\n\t		vmovaps	0x120(%%rcx),%%ymm13	\n\t"\
		"vmovaps	     (%%rsi),%%ymm2		\n\t"/* cc1 */"			vmovaps	0x040(%%rsi),%%ymm10	\n\t"/* cc3 */\
		"vmovaps	0x020(%%rsi),%%ymm3					\n\t		vmovaps	0x060(%%rsi),%%ymm11	\n\t"\
		"vmovaps	%%ymm4,%%ymm6						\n\t		vmovaps	%%ymm12,%%ymm14		\n\t"\
		"vmovaps	%%ymm5,%%ymm7						\n\t		vmovaps	%%ymm13,%%ymm15		\n\t"\
		"vmulpd		%%ymm10,%%ymm4,%%ymm4				\n\t		vmulpd		%%ymm3,%%ymm12,%%ymm12		\n\t"\
		"vmulpd		%%ymm10,%%ymm5,%%ymm5				\n\t		vmulpd		%%ymm3,%%ymm13,%%ymm13		\n\t"\
		"vmovaps	     (%%rdx),%%ymm0					\n\t		vmovaps	0x100(%%rdx),%%ymm8 	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm1					\n\t		vmovaps	0x120(%%rdx),%%ymm9 	\n\t"\
	" vfmadd231pd	%%ymm11,%%ymm7,%%ymm4				\n\t	 vfmadd231pd	%%ymm2,%%ymm15,%%ymm12		\n\t"\
	"vfnmadd231pd	%%ymm11,%%ymm6,%%ymm5				\n\t	vfnmadd231pd	%%ymm2,%%ymm14,%%ymm13		\n\t"\
		"vmovaps	%%ymm0,%%ymm6						\n\t		vmovaps	%%ymm8 ,%%ymm14		\n\t"\
		"vmovaps	%%ymm1,%%ymm7						\n\t		vmovaps	%%ymm9 ,%%ymm15		\n\t"\
		"vmulpd		%%ymm3,%%ymm0,%%ymm0				\n\t		vmulpd		%%ymm11,%%ymm8 ,%%ymm8 		\n\t"\
		"vmulpd		%%ymm3,%%ymm1,%%ymm1				\n\t		vmulpd		%%ymm11,%%ymm9 ,%%ymm9 		\n\t"\
	"vfnmadd231pd	%%ymm2,%%ymm7,%%ymm0				\n\t	 vfmadd231pd	%%ymm10,%%ymm15,%%ymm8 		\n\t"\
	" vfmadd231pd	%%ymm2,%%ymm6,%%ymm1				\n\t	vfnmadd231pd	%%ymm10,%%ymm14,%%ymm9 		\n\t"\
		"vmovaps	%%ymm5,%%ymm7						\n\t		vmovaps	%%ymm13,%%ymm15		\n\t"\
		"vmovaps	%%ymm4,%%ymm6						\n\t		vmovaps	%%ymm12,%%ymm14		\n\t"\
		"vaddpd	%%ymm0,%%ymm6,%%ymm6					\n\t		vsubpd	%%ymm8 ,%%ymm12,%%ymm12		\n\t"\
		"vaddpd	%%ymm1,%%ymm7,%%ymm7					\n\t		vsubpd	%%ymm9 ,%%ymm13,%%ymm13		\n\t"\
		"vsubpd	%%ymm0,%%ymm4,%%ymm4					\n\t		vaddpd	%%ymm8 ,%%ymm14,%%ymm14		\n\t"\
		"vsubpd	%%ymm1,%%ymm5,%%ymm5					\n\t		vaddpd	%%ymm9 ,%%ymm15,%%ymm15		\n\t"\
		"subq	$0x400,%%rcx			\n\t"\
		"subq	$0x400,%%rdx			\n\t"\
		"subq	$0x040,%%rsi	/* cc0 */	\n\t"\
		"vmovaps	     (%%rdx),%%ymm2					\n\t		vmovaps	0x100(%%rdx),%%ymm11	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm0					\n\t		vmovaps	0x120(%%rdx),%%ymm9 	\n\t"\
		"vmovaps	0x020(%%rsi),%%ymm3					\n\t		vmovaps	0x020(%%rsi),%%ymm10	\n\t"\
		"vmovaps	%%ymm2,%%ymm1						\n\t		vmovaps	%%ymm11,%%ymm8 		\n\t"\
		"vmulpd	%%ymm3 ,%%ymm2,%%ymm2					\n\t		vmulpd	%%ymm10,%%ymm11,%%ymm11		\n\t"\
		"vmulpd	%%ymm0 ,%%ymm3,%%ymm3					\n\t		vmulpd	%%ymm9 ,%%ymm10,%%ymm10		\n\t"\
	"vmovaps %%ymm8,(%%rax) \n\t vmovaps (%%rsi),%%ymm8	\n\t"/* isrt2 */\
	" vfmadd231pd	%%ymm8 ,%%ymm0,%%ymm2				\n\t	 vfmadd132pd	%%ymm8 ,%%ymm11,%%ymm9 		\n\t"\
	"vfnmadd231pd	%%ymm8 ,%%ymm1,%%ymm3				\n\t	 vfmsub132pd	(%%rax),%%ymm10,%%ymm8 		\n\t"\
		"vmovaps	     (%%rcx),%%ymm0					\n\t		vmovaps	0x100(%%rcx),%%ymm10	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm1					\n\t		vmovaps	0x120(%%rcx),%%ymm11	\n\t"\
		"vsubpd	%%ymm2,%%ymm0,%%ymm0					\n\t		vsubpd	%%ymm8 ,%%ymm10,%%ymm10		\n\t"\
		"vsubpd	%%ymm3,%%ymm1,%%ymm1					\n\t		vsubpd	%%ymm9 ,%%ymm11,%%ymm11		\n\t"\
	"vfmadd132pd	(%%r8),%%ymm0,%%ymm2				\n\t	vfmadd132pd	(%%r8) ,%%ymm10,%%ymm8 		\n\t"\
	"vfmadd132pd	(%%r8),%%ymm1,%%ymm3				\n\t	vfmadd132pd	(%%r8) ,%%ymm11,%%ymm9 		\n\t"\
		"\n\t"\
		/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c03): */			/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c07): */\
		"addq	$0x6c0,%%rsi	\n\t"/* c03; c07 += 0x100: */\
		"vsubpd	%%ymm4,%%ymm2,%%ymm2					\n\t		vsubpd	%%ymm12,%%ymm10,%%ymm10		\n\t"\
		"vsubpd	%%ymm7,%%ymm0,%%ymm0					\n\t		vsubpd	%%ymm15,%%ymm8 ,%%ymm8 		\n\t"\
		"vsubpd	%%ymm5,%%ymm3,%%ymm3					\n\t		vsubpd	%%ymm13,%%ymm11,%%ymm11		\n\t"\
		"vsubpd	%%ymm6,%%ymm1,%%ymm1					\n\t		vsubpd	%%ymm14,%%ymm9 ,%%ymm9 		\n\t"\
	"vmovaps %%ymm14,(%%rax) \n\t vmovaps (%%r8),%%ymm14\n\t"/* 2.0 */\
	"vfmadd132pd	%%ymm14,%%ymm2,%%ymm4				\n\t	vfmadd132pd	%%ymm14,%%ymm10,%%ymm12		\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm0,%%ymm7				\n\t	vfmadd132pd	%%ymm14,%%ymm8 ,%%ymm15		\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm3,%%ymm5				\n\t	vfmadd132pd	%%ymm14,%%ymm11,%%ymm13		\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm1,%%ymm6				\n\t	vfmadd132pd	(%%rax),%%ymm9 ,%%ymm14		\n\t"\
		"vmovaps	%%ymm2,     (%%rcx)					\n\t		vmovaps	%%ymm10,0x100(%%rcx)	\n\t"\
		"vmovaps	%%ymm0,0x020(%%rdx)					\n\t		vmovaps	%%ymm8 ,0x120(%%rdx)	\n\t"\
		"vmovaps	%%ymm3,0x020(%%rcx)					\n\t		vmovaps	%%ymm11,0x120(%%rcx)	\n\t"\
		"vmovaps	%%ymm6,     (%%rdx)					\n\t		vmovaps	%%ymm14,0x100(%%rdx)	\n\t"\
		"vmovaps	%%ymm4,%%ymm2						\n\t		vmovaps	%%ymm12,%%ymm10		\n\t"\
		"vmovaps	%%ymm7,%%ymm0						\n\t		vmovaps	%%ymm15,%%ymm8 		\n\t"\
		"vmovaps	%%ymm5,%%ymm3						\n\t		vmovaps	%%ymm13,%%ymm11		\n\t"\
		"vmovaps	%%ymm1,%%ymm6						\n\t		vmovaps	%%ymm9 ,%%ymm14		\n\t"\
		"vmulpd		     (%%rsi),%%ymm4,%%ymm4			\n\t		vmulpd		0x100(%%rsi),%%ymm12,%%ymm12	\n\t"\
		"vmulpd		0x040(%%rsi),%%ymm7,%%ymm7			\n\t		vmulpd		0x140(%%rsi),%%ymm15,%%ymm15	\n\t"\
		"vmulpd		     (%%rsi),%%ymm5,%%ymm5			\n\t		vmulpd		0x100(%%rsi),%%ymm13,%%ymm13	\n\t"\
		"vmulpd		0x040(%%rsi),%%ymm1,%%ymm1			\n\t		vmulpd		0x140(%%rsi),%%ymm9 ,%%ymm9 	\n\t"\
	" vfmadd231pd	0x020(%%rsi),%%ymm3,%%ymm4			\n\t	 vfmadd231pd	0x120(%%rsi),%%ymm11,%%ymm12	\n\t"\
	" vfmadd231pd	0x060(%%rsi),%%ymm6,%%ymm7			\n\t	 vfmadd231pd	0x160(%%rsi),%%ymm14,%%ymm15	\n\t"\
	"vfnmadd231pd	0x020(%%rsi),%%ymm2,%%ymm5			\n\t	vfnmadd231pd	0x120(%%rsi),%%ymm10,%%ymm13	\n\t"\
	"vfnmadd231pd	0x060(%%rsi),%%ymm0,%%ymm1			\n\t	vfnmadd231pd	0x160(%%rsi),%%ymm8 ,%%ymm9 	\n\t"\
		"vmovaps	%%ymm5,0x020(%%rax)					\n\t		vmovaps	%%ymm13,0x020(%%r10)	\n\t"\
		"vmovaps	%%ymm1,0x020(%%rbx)					\n\t		vmovaps	%%ymm9 ,0x020(%%r11)	\n\t"\
		"vmovaps	%%ymm4,     (%%rax)					\n\t		vmovaps	%%ymm12,     (%%r10)	\n\t"\
		"vmovaps	%%ymm7,     (%%rbx)					\n\t		vmovaps	%%ymm15,     (%%r11)	\n\t"\
		"addq	%%rdi,%%rax								\n\t		addq	%%rdi,%%r10			\n\t"\
		"addq	%%rdi,%%rbx								\n\t		addq	%%rdi,%%r11			\n\t"\
		"addq	$0x080,%%rsi	\n\t"/* c13; c17 += 0x100: */\
		"vmovaps	     (%%rcx),%%ymm4					\n\t		vmovaps	0x100(%%rcx),%%ymm12	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm0					\n\t		vmovaps	0x120(%%rdx),%%ymm8 	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5					\n\t		vmovaps	0x120(%%rcx),%%ymm13	\n\t"\
		"vmovaps	     (%%rdx),%%ymm6					\n\t		vmovaps	0x100(%%rdx),%%ymm14	\n\t"\
		"vmovaps	%%ymm4,%%ymm2						\n\t		vmovaps	%%ymm12,%%ymm10		\n\t"\
		"vmovaps	%%ymm0,%%ymm1						\n\t		vmovaps	%%ymm8 ,%%ymm9 		\n\t"\
		"vmovaps	%%ymm5,%%ymm3						\n\t		vmovaps	%%ymm13,%%ymm11		\n\t"\
		"vmovaps	%%ymm6,%%ymm7						\n\t		vmovaps	%%ymm14,%%ymm15		\n\t"\
		"vmulpd		     (%%rsi),%%ymm4,%%ymm4			\n\t		vmulpd		0x100(%%rsi),%%ymm12,%%ymm12	\n\t"\
		"vmulpd		0x040(%%rsi),%%ymm0,%%ymm0			\n\t		vmulpd		0x140(%%rsi),%%ymm8 ,%%ymm8 	\n\t"\
		"vmulpd		     (%%rsi),%%ymm5,%%ymm5			\n\t		vmulpd		0x100(%%rsi),%%ymm13,%%ymm13	\n\t"\
		"vmulpd		0x040(%%rsi),%%ymm6,%%ymm6			\n\t		vmulpd		0x140(%%rsi),%%ymm14,%%ymm14	\n\t"\
	" vfmadd231pd	0x020(%%rsi),%%ymm3,%%ymm4			\n\t	 vfmadd231pd	0x120(%%rsi),%%ymm11,%%ymm12	\n\t"\
	" vfmadd231pd	0x060(%%rsi),%%ymm7,%%ymm0			\n\t	 vfmadd231pd	0x160(%%rsi),%%ymm15,%%ymm8 	\n\t"\
	"vfnmadd231pd	0x020(%%rsi),%%ymm2,%%ymm5			\n\t	vfnmadd231pd	0x120(%%rsi),%%ymm10,%%ymm13	\n\t"\
	"vfnmadd231pd	0x060(%%rsi),%%ymm1,%%ymm6			\n\t	vfnmadd231pd	0x160(%%rsi),%%ymm9 ,%%ymm14	\n\t"\
		"vmovaps	%%ymm5,0x020(%%rax)					\n\t		vmovaps	%%ymm13,0x020(%%r10)	\n\t"\
		"vmovaps	%%ymm6,0x020(%%rbx)					\n\t		vmovaps	%%ymm14,0x020(%%r11)	\n\t"\
		"vmovaps	%%ymm4,     (%%rax)					\n\t		vmovaps	%%ymm12,     (%%r10)	\n\t"\
		"vmovaps	%%ymm0,     (%%rbx)					\n\t		vmovaps	%%ymm8 ,     (%%r11)	\n\t"\
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
		 ,[__r10] "m" (Xr10)\
		 ,[__r20] "m" (Xr20)\
		 ,[__r30] "m" (Xr30)\
		 ,[__isrt2] "m" (Xisrt2)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r8","r9","r10","r11","r12","r13","r14","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
	);\
	}

  #endif	// ALL_FMA ?

#elif defined(USE_AVX)	// AVX and AVX2 both use 256-bit registers

  #define USE_64BIT_ASM_STYLE	1	// My x86 timings indicate fancier versions below using xmm0-15 for the radix-32 DFT is faster.

  #if !USE_64BIT_ASM_STYLE	// USE_64BIT_ASM_STYLE = False: Use simple 64-bit-ified version of the analogous 32-bit ASM macros, i.e. using just xmm0-7.

	/*
	For GCC-macro version of this, use that isrt2 + 0x020,0x060,0x0a0 = cc0,cc1,cc3,
	and isrt2 + 0x0e0,0x1e0,0x2e0,0x3e0,0x4e0,0x5e0,0x6e0,0x7e0 = c00,04,02,06,01,05,03,07
	in order to reduce number of args to <= the GCC-allowed maximum of 30:
	*/
	#define SSE2_RADIX32_DIF_TWIDDLE(Xadd0,Xp01,Xp02,Xp03,Xp04,Xp08,Xp0C,Xp10,Xp18,Xr00)\
	{\
	__asm__ volatile (\
		"/*...Block 1: */	\n\t"\
		"movq	%[__add0],%%rax	\n\t"\
		"movslq	%[__p08],%%rbx	\n\t"\
		"movslq	%[__p10],%%rcx	\n\t"\
		"movslq	%[__p18],%%rdx	\n\t"\
		"shlq	$3,%%rbx	\n\t"\
		"shlq	$3,%%rcx	\n\t"\
		"shlq	$3,%%rdx	\n\t"\
		"addq	%%rax,%%rbx	\n\t"\
		"addq	%%rax,%%rcx	\n\t"\
		"addq	%%rax,%%rdx	\n\t"\
		"/* SSE2_RADIX4_DIF_4TWIDDLE_A(r00,c10) */	\n\t"\
		"movq	%[__r00],%%rsi	\n\t"\
		"vmovaps	     (%%rax),%%ymm0	\n\t"\
		"vmovaps	     (%%rcx),%%ymm4	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5	\n\t"\
		"vmovaps	0x920(%%rsi),%%ymm2	/* c10 */	\n\t"\
		"vmovaps	0x940(%%rsi),%%ymm3	\n\t"\
		"vmovaps	%%ymm4,%%ymm6	\n\t"\
		"vmovaps	%%ymm5,%%ymm7	\n\t"\
		"vmulpd	%%ymm2,%%ymm4,%%ymm4	\n\t"\
		"vmulpd	%%ymm2,%%ymm5,%%ymm5	\n\t"\
		"vmulpd	%%ymm3,%%ymm6,%%ymm6	\n\t"\
		"vmovaps	%%ymm0,%%ymm2	\n\t"\
		"vmulpd	%%ymm3,%%ymm7,%%ymm7	\n\t"\
		"vaddpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
		"vmovaps	%%ymm1,%%ymm3	\n\t"\
		"vsubpd	%%ymm7,%%ymm4,%%ymm4	\n\t"\
		"vaddpd	%%ymm4,%%ymm0,%%ymm0	\n\t"\
		"vaddpd	%%ymm5,%%ymm1,%%ymm1	\n\t"\
		"vsubpd	%%ymm4,%%ymm2,%%ymm2	\n\t"\
		"vsubpd	%%ymm5,%%ymm3,%%ymm3	\n\t"\
		"vmovaps	     (%%rdx),%%ymm4	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm5	\n\t"\
		"vmovaps	     (%%rdx),%%ymm6	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm7	\n\t"\
		"vmulpd	0x9a0(%%rsi),%%ymm4,%%ymm4	/* c18 */	\n\t"\
		"vmulpd	0x9a0(%%rsi),%%ymm5,%%ymm5	\n\t"\
		"vmulpd	0x9c0(%%rsi),%%ymm6,%%ymm6	\n\t"\
		"vmulpd	0x9c0(%%rsi),%%ymm7,%%ymm7	\n\t"\
		"vaddpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
		"vsubpd	%%ymm7,%%ymm4,%%ymm4	\n\t"\
		"vmovaps	%%ymm5,0x020(%%rsi)	\n\t"\
		"vmovaps	%%ymm4,     (%%rsi)	/* r00 */	\n\t"\
		"vmovaps	     (%%rbx),%%ymm4	\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm5	\n\t"\
		"vmovaps	     (%%rbx),%%ymm6	\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm7	\n\t"\
		"vmulpd	0x960(%%rsi),%%ymm4,%%ymm4	/* c08 */	\n\t"\
		"vmulpd	0x960(%%rsi),%%ymm5,%%ymm5	\n\t"\
		"vmulpd	0x980(%%rsi),%%ymm6,%%ymm6	\n\t"\
		"vmulpd	0x980(%%rsi),%%ymm7,%%ymm7	\n\t"\
		"vaddpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
		"vsubpd	%%ymm7,%%ymm4,%%ymm4	\n\t"\
		"vmovaps	%%ymm5,%%ymm7	\n\t"\
		"vmovaps	%%ymm4,%%ymm6	\n\t"\
		"vsubpd	     (%%rsi),%%ymm4,%%ymm4	/* r00 */	\n\t"\
		"vsubpd	0x020(%%rsi),%%ymm5,%%ymm5	\n\t"\
		"vaddpd	     (%%rsi),%%ymm6,%%ymm6	\n\t"\
		"vaddpd	0x020(%%rsi),%%ymm7,%%ymm7	\n\t"\
		"vsubpd	%%ymm6,%%ymm0,%%ymm0	\n\t"\
		"vsubpd	%%ymm5,%%ymm2,%%ymm2	\n\t"\
		"vsubpd	%%ymm7,%%ymm1,%%ymm1	\n\t"\
		"vsubpd	%%ymm4,%%ymm3,%%ymm3	\n\t"\
		"vmovaps	%%ymm0,0x080(%%rsi)	\n\t"\
		"vmovaps	%%ymm2,0x040(%%rsi)	\n\t"\
		"vmovaps	%%ymm1,0x0a0(%%rsi)	\n\t"\
		"vmovaps	%%ymm3,0x0e0(%%rsi)	\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6	\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7	\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4	\n\t"\
		"vaddpd	%%ymm0,%%ymm6,%%ymm6	\n\t"\
		"vaddpd	%%ymm2,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm1,%%ymm7,%%ymm7	\n\t"\
		"vaddpd	%%ymm3,%%ymm4,%%ymm4	\n\t"\
		"vmovaps	%%ymm6,     (%%rsi)	\n\t"\
		"vmovaps	%%ymm5,0x0c0(%%rsi)	\n\t"\
		"vmovaps	%%ymm7,0x020(%%rsi)	\n\t"\
		"vmovaps	%%ymm4,0x060(%%rsi)	\n\t"\
		"/***************************************/	\n\t"\
		"movslq	%[__p04],%%rdi	\n\t"\
		"shlq	$3,%%rdi	\n\t"\
		"addq	%%rdi,%%rax	/* &a[j1+p4] */	\n\t"\
		"addq	%%rdi,%%rbx	\n\t"\
		"addq	%%rdi,%%rcx	\n\t"\
		"addq	%%rdi,%%rdx	\n\t"\
		"/* SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO_B(r08,c04) */	\n\t"\
		"vmovaps	     (%%rax),%%ymm0	\n\t"\
		"vmovaps	     (%%rcx),%%ymm4	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5	\n\t"\
		"vmovaps	0x9e0(%%rsi),%%ymm6	/* c04 */	\n\t"\
		"vmovaps	0xa00(%%rsi),%%ymm7	\n\t"\
		"vmovaps	%%ymm0,%%ymm2	\n\t"\
		"vmovaps	%%ymm1,%%ymm3	\n\t"\
		"vmulpd	%%ymm6,%%ymm0,%%ymm0	\n\t"\
		"vmulpd	%%ymm6,%%ymm1,%%ymm1	\n\t"\
		"vmulpd	%%ymm7,%%ymm2,%%ymm2	\n\t"\
		"vmulpd	%%ymm7,%%ymm3,%%ymm3	\n\t"\
		"vmovaps	%%ymm4,%%ymm6	\n\t"\
		"vaddpd	%%ymm2,%%ymm1,%%ymm1	\n\t"\
		"vmovaps	%%ymm5,%%ymm7	\n\t"\
		"vmulpd	0xa20(%%rsi),%%ymm4,%%ymm4	/* c14 */	\n\t"\
		"vsubpd	%%ymm3,%%ymm0,%%ymm0	\n\t"\
		"vmulpd	0xa20(%%rsi),%%ymm5,%%ymm5	\n\t"\
		"vmulpd	0xa40(%%rsi),%%ymm6,%%ymm6	\n\t"\
		"vmovaps	%%ymm0,%%ymm2	\n\t"\
		"vmulpd	0xa40(%%rsi),%%ymm7,%%ymm7	\n\t"\
		"vaddpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
		"vmovaps	%%ymm1,%%ymm3	\n\t"\
		"vsubpd	%%ymm7,%%ymm4,%%ymm4	\n\t"\
		"vaddpd	%%ymm4,%%ymm0,%%ymm0	\n\t"\
		"vaddpd	%%ymm5,%%ymm1,%%ymm1	\n\t"\
		"vsubpd	%%ymm4,%%ymm2,%%ymm2	\n\t"\
		"vsubpd	%%ymm5,%%ymm3,%%ymm3	\n\t"\
		"vmovaps	     (%%rdx),%%ymm4	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm5	\n\t"\
		"vmovaps	     (%%rdx),%%ymm6	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm7	\n\t"\
		"vmulpd	0xaa0(%%rsi),%%ymm4,%%ymm4	/* c1C */	\n\t"\
		"vmulpd	0xaa0(%%rsi),%%ymm5,%%ymm5	\n\t"\
		"vmulpd	0xac0(%%rsi),%%ymm6,%%ymm6	\n\t"\
		"vmulpd	0xac0(%%rsi),%%ymm7,%%ymm7	\n\t"\
		"addq	$0x100,%%rsi	/* r08 */	\n\t"\
		"vaddpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
		"vsubpd	%%ymm7,%%ymm4,%%ymm4	\n\t"\
		"vmovaps	%%ymm5,0x020(%%rsi)	\n\t"\
		"vmovaps	%%ymm4,     (%%rsi)	\n\t"\
		"vmovaps	     (%%rbx),%%ymm4	\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm5	\n\t"\
		"vmovaps	     (%%rbx),%%ymm6	\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm7	\n\t"\
		"vmulpd	0x960(%%rsi),%%ymm4,%%ymm4	/* c0C */	\n\t"\
		"vmulpd	0x960(%%rsi),%%ymm5,%%ymm5	\n\t"\
		"vmulpd	0x980(%%rsi),%%ymm6,%%ymm6	\n\t"\
		"vmulpd	0x980(%%rsi),%%ymm7,%%ymm7	\n\t"\
		"vaddpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
		"vsubpd	%%ymm7,%%ymm4,%%ymm4	\n\t"\
		"vmovaps	%%ymm5,%%ymm7	\n\t"\
		"vmovaps	%%ymm4,%%ymm6	\n\t"\
		"vsubpd	     (%%rsi),%%ymm4,%%ymm4	\n\t"\
		"vsubpd	0x020(%%rsi),%%ymm5,%%ymm5	\n\t"\
		"vaddpd	     (%%rsi),%%ymm6,%%ymm6	\n\t"\
		"vaddpd	0x020(%%rsi),%%ymm7,%%ymm7	\n\t"\
		"vsubpd	%%ymm6,%%ymm0,%%ymm0	\n\t"\
		"vsubpd	%%ymm7,%%ymm1,%%ymm1	\n\t"\
		"vmovaps	%%ymm0,0x080(%%rsi)	\n\t"\
		"vsubpd	%%ymm5,%%ymm2,%%ymm2	\n\t"\
		"vmovaps	%%ymm1,0x0a0(%%rsi)	\n\t"\
		"vsubpd	%%ymm4,%%ymm3,%%ymm3	\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6	\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7	\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4	\n\t"\
		"vaddpd	%%ymm0,%%ymm6,%%ymm6	\n\t"\
		"vaddpd	%%ymm2,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm1,%%ymm7,%%ymm7	\n\t"\
		"vaddpd	%%ymm3,%%ymm4,%%ymm4	\n\t"\
		"vmovaps	%%ymm6,     (%%rsi)	\n\t"\
		"vmovaps	%%ymm7,0x020(%%rsi)	\n\t"\
		"vmovaps	0x700(%%rsi),%%ymm0	/* isrt2 */	\n\t"\
		"vmovaps	%%ymm2,%%ymm6	\n\t"\
		"vmovaps	%%ymm5,%%ymm7	\n\t"\
		"vsubpd	%%ymm4,%%ymm2,%%ymm2	\n\t"\
		"vsubpd	%%ymm3,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm4,%%ymm6,%%ymm6	\n\t"\
		"vaddpd	%%ymm3,%%ymm7,%%ymm7	\n\t"\
		"vmulpd	%%ymm0,%%ymm2,%%ymm2	\n\t"\
		"vmulpd	%%ymm0,%%ymm5,%%ymm5	\n\t"\
		"vmulpd	%%ymm0,%%ymm6,%%ymm6	\n\t"\
		"vmulpd	%%ymm0,%%ymm7,%%ymm7	\n\t"\
		"vmovaps	%%ymm2,0x040(%%rsi)	\n\t"\
		"vmovaps	%%ymm5,0x0c0(%%rsi)	\n\t"\
		"vmovaps	%%ymm6,0x060(%%rsi)	\n\t"\
		"vmovaps	%%ymm7,0x0e0(%%rsi)	\n\t"\
		"/***************************************/	\n\t"\
		"\n\t"\
		"/* SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r00) */	\n\t"\
		"subq	$0x100,%%rsi	/* r00 */	\n\t"\
		"vmovaps	     (%%rsi),%%ymm0	\n\t"\
		"vmovaps	0x080(%%rsi),%%ymm4	\n\t"\
		"vmovaps	0x020(%%rsi),%%ymm1	\n\t"\
		"vmovaps	0x0a0(%%rsi),%%ymm5	\n\t"\
		"vmovaps	0x100(%%rsi),%%ymm2	\n\t"\
		"vmovaps	0x1a0(%%rsi),%%ymm7	\n\t"\
		"vmovaps	0x120(%%rsi),%%ymm3	\n\t"\
		"vmovaps	0x180(%%rsi),%%ymm6	\n\t"\
		"vsubpd   %%ymm2,%%ymm0,%%ymm0	\n\t"\
		"vsubpd   %%ymm7,%%ymm4,%%ymm4	\n\t"\
		"vsubpd   %%ymm3,%%ymm1,%%ymm1	\n\t"\
		"vsubpd   %%ymm6,%%ymm5,%%ymm5	\n\t"\
		"vaddpd   %%ymm2,%%ymm2,%%ymm2	\n\t"\
		"vaddpd   %%ymm7,%%ymm7,%%ymm7	\n\t"\
		"vaddpd   %%ymm3,%%ymm3,%%ymm3	\n\t"\
		"vaddpd   %%ymm6,%%ymm6,%%ymm6	\n\t"\
		"vaddpd   %%ymm0,%%ymm2,%%ymm2	\n\t"\
		"vaddpd   %%ymm4,%%ymm7,%%ymm7	\n\t"\
		"vaddpd   %%ymm1,%%ymm3,%%ymm3	\n\t"\
		"vaddpd   %%ymm5,%%ymm6,%%ymm6	\n\t"\
		"vmovaps	%%ymm0,0x100(%%rsi)	\n\t"\
		"vmovaps	%%ymm4,0x080(%%rsi)	\n\t"\
		"vmovaps	%%ymm1,0x120(%%rsi)	\n\t"\
		"vmovaps	%%ymm5,0x1a0(%%rsi)	\n\t"\
		"vmovaps	%%ymm2,     (%%rsi)	\n\t"\
		"vmovaps	%%ymm7,0x180(%%rsi)	\n\t"\
		"vmovaps	%%ymm3,0x020(%%rsi)	\n\t"\
		"vmovaps	%%ymm6,0x0a0(%%rsi)	\n\t"\
		"vmovaps	0x040(%%rsi),%%ymm0	\n\t"\
		"vmovaps	0x0c0(%%rsi),%%ymm4	\n\t"\
		"vmovaps	0x060(%%rsi),%%ymm1	\n\t"\
		"vmovaps	0x0e0(%%rsi),%%ymm5	\n\t"\
		"vmovaps	0x140(%%rsi),%%ymm2	\n\t"\
		"vmovaps	0x1e0(%%rsi),%%ymm7	\n\t"\
		"vmovaps	0x160(%%rsi),%%ymm3	\n\t"\
		"vmovaps	0x1c0(%%rsi),%%ymm6	\n\t"\
		"vsubpd   %%ymm2,%%ymm0,%%ymm0	\n\t"\
		"vsubpd   %%ymm7,%%ymm4,%%ymm4	\n\t"\
		"vsubpd   %%ymm3,%%ymm1,%%ymm1	\n\t"\
		"vsubpd   %%ymm6,%%ymm5,%%ymm5	\n\t"\
		"vaddpd   %%ymm2,%%ymm2,%%ymm2	\n\t"\
		"vaddpd   %%ymm7,%%ymm7,%%ymm7	\n\t"\
		"vaddpd   %%ymm3,%%ymm3,%%ymm3	\n\t"\
		"vaddpd   %%ymm6,%%ymm6,%%ymm6	\n\t"\
		"vaddpd   %%ymm0,%%ymm2,%%ymm2	\n\t"\
		"vaddpd   %%ymm4,%%ymm7,%%ymm7	\n\t"\
		"vaddpd   %%ymm1,%%ymm3,%%ymm3	\n\t"\
		"vaddpd   %%ymm5,%%ymm6,%%ymm6	\n\t"\
		"vmovaps	%%ymm0,0x140(%%rsi)	\n\t"\
		"vmovaps	%%ymm4,0x0c0(%%rsi)	\n\t"\
		"vmovaps	%%ymm1,0x160(%%rsi)	\n\t"\
		"vmovaps	%%ymm5,0x1e0(%%rsi)	\n\t"\
		"vmovaps	%%ymm2,0x040(%%rsi)	\n\t"\
		"vmovaps	%%ymm7,0x1c0(%%rsi)	\n\t"\
		"vmovaps	%%ymm3,0x060(%%rsi)	\n\t"\
		"vmovaps	%%ymm6,0x0e0(%%rsi)	\n\t"\
		"/***************************************/	\n\t"\
		"\n\t"\
		"/*...Block 2: */	\n\t"\
		"subq	%%rdi,%%rax	/* &a[j1] */	\n\t"\
		"subq	%%rdi,%%rbx	\n\t"\
		"subq	%%rdi,%%rcx	\n\t"\
		"subq	%%rdi,%%rdx	\n\t"\
		"movslq	%[__p02],%%rdi	/* Do things this way [rather than repeatedly adding p1] since array-padding scheme means p2 == p+p1 not guaranteed. */	\n\t"\
		"shlq	$3,%%rdi	\n\t"\
		"addq	%%rdi,%%rax	/* &a[j1+p2] */	\n\t"\
		"addq	%%rdi,%%rbx	\n\t"\
		"addq	%%rdi,%%rcx	\n\t"\
		"addq	%%rdi,%%rdx	\n\t"\
		"/* SSE2_RADIX4_DIF_4TWIDDLE_B(r10,c02) */	\n\t"\
		"addq	$0x200,%%rsi	/* r10 */	\n\t"\
		"vmovaps	     (%%rax),%%ymm0	\n\t"\
		"vmovaps	     (%%rcx),%%ymm4	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5	\n\t"\
		"vmovaps	0x8e0(%%rsi),%%ymm6	/* c02 */	\n\t"\
		"vmovaps	0x900(%%rsi),%%ymm7	\n\t"\
		"vmovaps	%%ymm0,%%ymm2	\n\t"\
		"vmovaps	%%ymm1,%%ymm3	\n\t"\
		"vmulpd   %%ymm6,%%ymm0,%%ymm0	\n\t"\
		"vmulpd   %%ymm6,%%ymm1,%%ymm1	\n\t"\
		"vmulpd   %%ymm7,%%ymm2,%%ymm2	\n\t"\
		"vmulpd   %%ymm7,%%ymm3,%%ymm3	\n\t"\
		"vmovaps	%%ymm4,%%ymm6	\n\t"\
		"vaddpd   %%ymm2,%%ymm1,%%ymm1	\n\t"\
		"vmovaps	%%ymm5,%%ymm7	\n\t"\
		"vmulpd   0x920(%%rsi),%%ymm4,%%ymm4	/* c12 */	\n\t"\
		"vsubpd   %%ymm3,%%ymm0,%%ymm0	\n\t"\
		"vmulpd   0x920(%%rsi),%%ymm5,%%ymm5	\n\t"\
		"vmulpd   0x940(%%rsi),%%ymm6,%%ymm6	\n\t"\
		"vmovaps	%%ymm0,%%ymm2	\n\t"\
		"vmulpd	0x940(%%rsi),%%ymm7,%%ymm7	\n\t"\
		"vaddpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
		"vmovaps	%%ymm1,%%ymm3	\n\t"\
		"vsubpd	%%ymm7,%%ymm4,%%ymm4	\n\t"\
		"vaddpd	%%ymm4,%%ymm0,%%ymm0	\n\t"\
		"vaddpd	%%ymm5,%%ymm1,%%ymm1	\n\t"\
		"vsubpd	%%ymm4,%%ymm2,%%ymm2	\n\t"\
		"vsubpd	%%ymm5,%%ymm3,%%ymm3	\n\t"\
		"vmovaps	     (%%rdx),%%ymm4	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm5	\n\t"\
		"vmovaps	     (%%rdx),%%ymm6	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm7	\n\t"\
		"vmulpd	0x9a0(%%rsi),%%ymm4,%%ymm4	/* c1A */	\n\t"\
		"vmulpd	0x9a0(%%rsi),%%ymm5,%%ymm5	\n\t"\
		"vmulpd	0x9c0(%%rsi),%%ymm6,%%ymm6	\n\t"\
		"vmulpd	0x9c0(%%rsi),%%ymm7,%%ymm7	\n\t"\
		"vaddpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
		"vsubpd	%%ymm7,%%ymm4,%%ymm4	\n\t"\
		"vmovaps	%%ymm5,0x020(%%rsi)	\n\t"\
		"vmovaps	%%ymm4,     (%%rsi)	\n\t"\
		"vmovaps	     (%%rbx),%%ymm4	\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm5	\n\t"\
		"vmovaps	     (%%rbx),%%ymm6	\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm7	\n\t"\
		"vmulpd	0x960(%%rsi),%%ymm4,%%ymm4	/* c0A */	\n\t"\
		"vmulpd	0x960(%%rsi),%%ymm5,%%ymm5	\n\t"\
		"vmulpd	0x980(%%rsi),%%ymm6,%%ymm6	\n\t"\
		"vmulpd	0x980(%%rsi),%%ymm7,%%ymm7	\n\t"\
		"vaddpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
		"vsubpd	%%ymm7,%%ymm4,%%ymm4	\n\t"\
		"vmovaps	%%ymm5,%%ymm7	\n\t"\
		"vmovaps	%%ymm4,%%ymm6	\n\t"\
		"vsubpd	     (%%rsi),%%ymm4,%%ymm4	\n\t"\
		"vsubpd	0x020(%%rsi),%%ymm5,%%ymm5	\n\t"\
		"vaddpd	     (%%rsi),%%ymm6,%%ymm6	\n\t"\
		"vaddpd	0x020(%%rsi),%%ymm7,%%ymm7	\n\t"\
		"vsubpd	%%ymm6,%%ymm0,%%ymm0	\n\t"\
		"vsubpd	%%ymm5,%%ymm2,%%ymm2	\n\t"\
		"vsubpd	%%ymm7,%%ymm1,%%ymm1	\n\t"\
		"vsubpd	%%ymm4,%%ymm3,%%ymm3	\n\t"\
		"vmovaps	%%ymm0,0x080(%%rsi)	\n\t"\
		"vmovaps	%%ymm2,0x040(%%rsi)	\n\t"\
		"vmovaps	%%ymm1,0x0a0(%%rsi)	\n\t"\
		"vmovaps	%%ymm3,0x0e0(%%rsi)	\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6	\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7	\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4	\n\t"\
		"vaddpd	%%ymm0,%%ymm6,%%ymm6	\n\t"\
		"vaddpd	%%ymm2,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm1,%%ymm7,%%ymm7	\n\t"\
		"vaddpd	%%ymm3,%%ymm4,%%ymm4	\n\t"\
		"vmovaps	%%ymm6,     (%%rsi)	\n\t"\
		"vmovaps	%%ymm5,0x0c0(%%rsi)	\n\t"\
		"vmovaps	%%ymm7,0x020(%%rsi)	\n\t"\
		"vmovaps	%%ymm4,0x060(%%rsi)	\n\t"\
		"/***************************************/	\n\t"\
		"movslq	%[__p04],%%rdi	\n\t"\
		"shlq	$3,%%rdi	\n\t"\
		"addq	%%rdi,%%rax	/* &a[j1+p6] */	\n\t"\
		"addq	%%rdi,%%rbx	\n\t"\
		"addq	%%rdi,%%rcx	\n\t"\
		"addq	%%rdi,%%rdx	\n\t"\
		"/* SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO_B(r18,c06) */	\n\t"\
		"/* rsi contains r10 */	\n\t"\
		"vmovaps	     (%%rax),%%ymm0	\n\t"\
		"vmovaps	     (%%rcx),%%ymm4	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5	\n\t"\
		"vmovaps	0x9e0(%%rsi),%%ymm6	/* c06 */	\n\t"\
		"vmovaps	0xa00(%%rsi),%%ymm7	\n\t"\
		"vmovaps	%%ymm0,%%ymm2	\n\t"\
		"vmovaps	%%ymm1,%%ymm3	\n\t"\
		"vmulpd	%%ymm6,%%ymm0,%%ymm0	\n\t"\
		"vmulpd	%%ymm6,%%ymm1,%%ymm1	\n\t"\
		"vmulpd	%%ymm7,%%ymm2,%%ymm2	\n\t"\
		"vmulpd	%%ymm7,%%ymm3,%%ymm3	\n\t"\
		"vmovaps	%%ymm4,%%ymm6	\n\t"\
		"vaddpd	%%ymm2,%%ymm1,%%ymm1	\n\t"\
		"vmovaps	%%ymm5,%%ymm7	\n\t"\
		"vmulpd	0xa20(%%rsi),%%ymm4,%%ymm4	/* c16 */	\n\t"\
		"vsubpd	%%ymm3,%%ymm0,%%ymm0	\n\t"\
		"vmulpd	0xa20(%%rsi),%%ymm5,%%ymm5	\n\t"\
		"vmulpd	0xa40(%%rsi),%%ymm6,%%ymm6	\n\t"\
		"vmovaps	%%ymm0,%%ymm2	\n\t"\
		"vmulpd	0xa40(%%rsi),%%ymm7,%%ymm7	\n\t"\
		"vaddpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
		"vmovaps	%%ymm1,%%ymm3	\n\t"\
		"vsubpd	%%ymm7,%%ymm4,%%ymm4	\n\t"\
		"vaddpd	%%ymm4,%%ymm0,%%ymm0	\n\t"\
		"vaddpd	%%ymm5,%%ymm1,%%ymm1	\n\t"\
		"vsubpd	%%ymm4,%%ymm2,%%ymm2	\n\t"\
		"vsubpd	%%ymm5,%%ymm3,%%ymm3	\n\t"\
		"vmovaps	     (%%rdx),%%ymm4	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm5	\n\t"\
		"vmovaps	     (%%rdx),%%ymm6	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm7	\n\t"\
		"vmulpd	0xaa0(%%rsi),%%ymm4,%%ymm4	/* c1E */	\n\t"\
		"vmulpd	0xaa0(%%rsi),%%ymm5,%%ymm5	\n\t"\
		"vmulpd	0xac0(%%rsi),%%ymm6,%%ymm6	\n\t"\
		"vmulpd	0xac0(%%rsi),%%ymm7,%%ymm7	\n\t"\
		"addq	$0x100,%%rsi	/* r18 */	\n\t"\
		"vaddpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
		"vsubpd	%%ymm7,%%ymm4,%%ymm4	\n\t"\
		"vmovaps	%%ymm5,0x020(%%rsi)	\n\t"\
		"vmovaps	%%ymm4,     (%%rsi)	\n\t"\
		"vmovaps	     (%%rbx),%%ymm4	\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm5	\n\t"\
		"vmovaps	     (%%rbx),%%ymm6	\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm7	\n\t"\
		"vmulpd	0x960(%%rsi),%%ymm4,%%ymm4	/* c0E */	\n\t"\
		"vmulpd	0x960(%%rsi),%%ymm5,%%ymm5	\n\t"\
		"vmulpd	0x980(%%rsi),%%ymm6,%%ymm6	\n\t"\
		"vmulpd	0x980(%%rsi),%%ymm7,%%ymm7	\n\t"\
		"vaddpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
		"vsubpd	%%ymm7,%%ymm4,%%ymm4	\n\t"\
		"vmovaps	%%ymm5,%%ymm7	\n\t"\
		"vmovaps	%%ymm4,%%ymm6	\n\t"\
		"vsubpd	     (%%rsi),%%ymm4,%%ymm4	\n\t"\
		"vsubpd	0x020(%%rsi),%%ymm5,%%ymm5	\n\t"\
		"vaddpd	     (%%rsi),%%ymm6,%%ymm6	\n\t"\
		"vaddpd	0x020(%%rsi),%%ymm7,%%ymm7	\n\t"\
		"vsubpd	%%ymm6,%%ymm0,%%ymm0	\n\t"\
		"vsubpd	%%ymm7,%%ymm1,%%ymm1	\n\t"\
		"vmovaps	%%ymm0,0x080(%%rsi)	\n\t"\
		"vsubpd	%%ymm5,%%ymm2,%%ymm2	\n\t"\
		"vmovaps	%%ymm1,0x0a0(%%rsi)	\n\t"\
		"vsubpd	%%ymm4,%%ymm3,%%ymm3	\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6	\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7	\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4	\n\t"\
		"vaddpd	%%ymm0,%%ymm6,%%ymm6	\n\t"\
		"vaddpd	%%ymm2,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm1,%%ymm7,%%ymm7	\n\t"\
		"vaddpd	%%ymm3,%%ymm4,%%ymm4	\n\t"\
		"vmovaps	%%ymm6,     (%%rsi)	\n\t"\
		"vmovaps	%%ymm7,0x020(%%rsi)	\n\t"\
		"vmovaps	0x500(%%rsi),%%ymm0	/* isrt2 */	\n\t"\
		"vmovaps	%%ymm2,%%ymm6	\n\t"\
		"vmovaps	%%ymm5,%%ymm7	\n\t"\
		"vsubpd	%%ymm4,%%ymm2,%%ymm2	\n\t"\
		"vsubpd	%%ymm3,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm4,%%ymm6,%%ymm6	\n\t"\
		"vaddpd	%%ymm3,%%ymm7,%%ymm7	\n\t"\
		"vmulpd	%%ymm0,%%ymm2,%%ymm2	\n\t"\
		"vmulpd	%%ymm0,%%ymm5,%%ymm5	\n\t"\
		"vmulpd	%%ymm0,%%ymm6,%%ymm6	\n\t"\
		"vmulpd	%%ymm0,%%ymm7,%%ymm7	\n\t"\
		"vmovaps	%%ymm2,0x040(%%rsi)	\n\t"\
		"vmovaps	%%ymm5,0x0c0(%%rsi)	\n\t"\
		"vmovaps	%%ymm6,0x060(%%rsi)	\n\t"\
		"vmovaps	%%ymm7,0x0e0(%%rsi)	\n\t"\
		"/***************************************/	\n\t"\
		"\n\t"\
		"/* SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r10) */	\n\t"\
		"subq		$0x100,%%rsi	/* r10 */	\n\t"\
		"vmovaps		     (%%rsi),%%ymm0	\n\t"\
		"vmovaps		0x080(%%rsi),%%ymm4	\n\t"\
		"vmovaps		0x020(%%rsi),%%ymm1	\n\t"\
		"vmovaps		0x0a0(%%rsi),%%ymm5	\n\t"\
		"vmovaps		0x100(%%rsi),%%ymm2	\n\t"\
		"vmovaps		0x1a0(%%rsi),%%ymm7	\n\t"\
		"vmovaps		0x120(%%rsi),%%ymm3	\n\t"\
		"vmovaps		0x180(%%rsi),%%ymm6	\n\t"\
		"vsubpd	%%ymm2,%%ymm0,%%ymm0	\n\t"\
		"vsubpd	%%ymm7,%%ymm4,%%ymm4	\n\t"\
		"vsubpd	%%ymm3,%%ymm1,%%ymm1	\n\t"\
		"vsubpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm2,%%ymm2,%%ymm2	\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7	\n\t"\
		"vaddpd	%%ymm3,%%ymm3,%%ymm3	\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6	\n\t"\
		"vaddpd	%%ymm0,%%ymm2,%%ymm2	\n\t"\
		"vaddpd	%%ymm4,%%ymm7,%%ymm7	\n\t"\
		"vaddpd	%%ymm1,%%ymm3,%%ymm3	\n\t"\
		"vaddpd	%%ymm5,%%ymm6,%%ymm6	\n\t"\
		"vmovaps		%%ymm0,0x100(%%rsi)	\n\t"\
		"vmovaps		%%ymm4,0x080(%%rsi)	\n\t"\
		"vmovaps		%%ymm1,0x120(%%rsi)	\n\t"\
		"vmovaps		%%ymm5,0x1a0(%%rsi)	\n\t"\
		"vmovaps		%%ymm2,     (%%rsi)	\n\t"\
		"vmovaps		%%ymm7,0x180(%%rsi)	\n\t"\
		"vmovaps		%%ymm3,0x020(%%rsi)	\n\t"\
		"vmovaps		%%ymm6,0x0a0(%%rsi)	\n\t"\
		"vmovaps		0x040(%%rsi),%%ymm0	\n\t"\
		"vmovaps		0x0c0(%%rsi),%%ymm4	\n\t"\
		"vmovaps		0x060(%%rsi),%%ymm1	\n\t"\
		"vmovaps		0x0e0(%%rsi),%%ymm5	\n\t"\
		"vmovaps		0x140(%%rsi),%%ymm2	\n\t"\
		"vmovaps		0x1e0(%%rsi),%%ymm7	\n\t"\
		"vmovaps		0x160(%%rsi),%%ymm3	\n\t"\
		"vmovaps		0x1c0(%%rsi),%%ymm6	\n\t"\
		"vsubpd	%%ymm2,%%ymm0,%%ymm0	\n\t"\
		"vsubpd	%%ymm7,%%ymm4,%%ymm4	\n\t"\
		"vsubpd	%%ymm3,%%ymm1,%%ymm1	\n\t"\
		"vsubpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm2,%%ymm2,%%ymm2	\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7	\n\t"\
		"vaddpd	%%ymm3,%%ymm3,%%ymm3	\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6	\n\t"\
		"vaddpd	%%ymm0,%%ymm2,%%ymm2	\n\t"\
		"vaddpd	%%ymm4,%%ymm7,%%ymm7	\n\t"\
		"vaddpd	%%ymm1,%%ymm3,%%ymm3	\n\t"\
		"vaddpd	%%ymm5,%%ymm6,%%ymm6	\n\t"\
		"vmovaps		%%ymm0,0x140(%%rsi)	\n\t"\
		"vmovaps		%%ymm4,0x0c0(%%rsi)	\n\t"\
		"vmovaps		%%ymm1,0x160(%%rsi)	\n\t"\
		"vmovaps		%%ymm5,0x1e0(%%rsi)	\n\t"\
		"vmovaps		%%ymm2,0x040(%%rsi)	\n\t"\
		"vmovaps		%%ymm7,0x1c0(%%rsi)	\n\t"\
		"vmovaps		%%ymm3,0x060(%%rsi)	\n\t"\
		"vmovaps		%%ymm6,0x0e0(%%rsi)	\n\t"\
		"/***************************************/	\n\t"\
		"\n\t"\
		"/*...Block 3: */	\n\t"\
		"movq	%[__add0],%%rax	\n\t"\
		"movslq	%[__p01],%%rdi	/* Do things this way [rather than repeatedly adding p1] since array-padding scheme means p2 == p+p1 not guaranteed. */	\n\t"\
		"movslq	%[__p08],%%rbx	/* Can't get these via simple load-one-and-shift-as-needed due to array padding scheme */	\n\t"\
		"movslq	%[__p10],%%rcx	\n\t"\
		"movslq	%[__p18],%%rdx	\n\t"\
		"shlq	$3,%%rdi	\n\t"\
		"shlq	$3,%%rbx	\n\t"\
		"shlq	$3,%%rcx	\n\t"\
		"shlq	$3,%%rdx	\n\t"\
		"addq	%%rdi,%%rax	/* &a[j1+p1] */	\n\t"\
		"addq	%%rax,%%rbx	\n\t"\
		"addq	%%rax,%%rcx	\n\t"\
		"addq	%%rax,%%rdx	\n\t"\
		"/* SSE2_RADIX4_DIF_4TWIDDLE_B(r20,c01) */	\n\t"\
		"addq	$0x200,%%rsi	/* r20 */	\n\t"\
		"vmovaps	     (%%rax),%%ymm0	\n\t"\
		"vmovaps	     (%%rcx),%%ymm4	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5	\n\t"\
		"vmovaps	0x8e0(%%rsi),%%ymm6	/* c01 */	\n\t"\
		"vmovaps	0x900(%%rsi),%%ymm7	\n\t"\
		"vmovaps	%%ymm0,%%ymm2	\n\t"\
		"vmovaps	%%ymm1,%%ymm3	\n\t"\
		"vmulpd   %%ymm6,%%ymm0,%%ymm0	\n\t"\
		"vmulpd   %%ymm6,%%ymm1,%%ymm1	\n\t"\
		"vmulpd   %%ymm7,%%ymm2,%%ymm2	\n\t"\
		"vmulpd   %%ymm7,%%ymm3,%%ymm3	\n\t"\
		"vmovaps	%%ymm4,%%ymm6	\n\t"\
		"vaddpd   %%ymm2,%%ymm1,%%ymm1	\n\t"\
		"vmovaps	%%ymm5,%%ymm7	\n\t"\
		"vmulpd   0x920(%%rsi),%%ymm4,%%ymm4	/* c11 */	\n\t"\
		"vsubpd   %%ymm3,%%ymm0,%%ymm0	\n\t"\
		"vmulpd   0x920(%%rsi),%%ymm5,%%ymm5	\n\t"\
		"vmulpd   0x940(%%rsi),%%ymm6,%%ymm6	\n\t"\
		"vmovaps	%%ymm0,%%ymm2	\n\t"\
		"vmulpd	0x940(%%rsi),%%ymm7,%%ymm7	\n\t"\
		"vaddpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
		"vmovaps	%%ymm1,%%ymm3	\n\t"\
		"vsubpd	%%ymm7,%%ymm4,%%ymm4	\n\t"\
		"vaddpd	%%ymm4,%%ymm0,%%ymm0	\n\t"\
		"vaddpd	%%ymm5,%%ymm1,%%ymm1	\n\t"\
		"vsubpd	%%ymm4,%%ymm2,%%ymm2	\n\t"\
		"vsubpd	%%ymm5,%%ymm3,%%ymm3	\n\t"\
		"vmovaps	     (%%rdx),%%ymm4	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm5	\n\t"\
		"vmovaps	     (%%rdx),%%ymm6	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm7	\n\t"\
		"vmulpd	0x9a0(%%rsi),%%ymm4,%%ymm4	/* c19 */	\n\t"\
		"vmulpd	0x9a0(%%rsi),%%ymm5,%%ymm5	\n\t"\
		"vmulpd	0x9c0(%%rsi),%%ymm6,%%ymm6	\n\t"\
		"vmulpd	0x9c0(%%rsi),%%ymm7,%%ymm7	\n\t"\
		"vaddpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
		"vsubpd	%%ymm7,%%ymm4,%%ymm4	\n\t"\
		"vmovaps	%%ymm5,0x020(%%rsi)	\n\t"\
		"vmovaps	%%ymm4,     (%%rsi)	\n\t"\
		"vmovaps	     (%%rbx),%%ymm4	\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm5	\n\t"\
		"vmovaps	     (%%rbx),%%ymm6	\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm7	\n\t"\
		"vmulpd	0x960(%%rsi),%%ymm4,%%ymm4	/* c01 */	\n\t"\
		"vmulpd	0x960(%%rsi),%%ymm5,%%ymm5	\n\t"\
		"vmulpd	0x980(%%rsi),%%ymm6,%%ymm6	\n\t"\
		"vmulpd	0x980(%%rsi),%%ymm7,%%ymm7	\n\t"\
		"vaddpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
		"vsubpd	%%ymm7,%%ymm4,%%ymm4	\n\t"\
		"vmovaps	%%ymm5,%%ymm7	\n\t"\
		"vmovaps	%%ymm4,%%ymm6	\n\t"\
		"vsubpd	     (%%rsi),%%ymm4,%%ymm4	\n\t"\
		"vsubpd	0x020(%%rsi),%%ymm5,%%ymm5	\n\t"\
		"vaddpd	     (%%rsi),%%ymm6,%%ymm6	\n\t"\
		"vaddpd	0x020(%%rsi),%%ymm7,%%ymm7	\n\t"\
		"vsubpd	%%ymm6,%%ymm0,%%ymm0	\n\t"\
		"vsubpd	%%ymm5,%%ymm2,%%ymm2	\n\t"\
		"vsubpd	%%ymm7,%%ymm1,%%ymm1	\n\t"\
		"vsubpd	%%ymm4,%%ymm3,%%ymm3	\n\t"\
		"vmovaps	%%ymm0,0x080(%%rsi)	\n\t"\
		"vmovaps	%%ymm2,0x040(%%rsi)	\n\t"\
		"vmovaps	%%ymm1,0x0a0(%%rsi)	\n\t"\
		"vmovaps	%%ymm3,0x0e0(%%rsi)	\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6	\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7	\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4	\n\t"\
		"vaddpd	%%ymm0,%%ymm6,%%ymm6	\n\t"\
		"vaddpd	%%ymm2,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm1,%%ymm7,%%ymm7	\n\t"\
		"vaddpd	%%ymm3,%%ymm4,%%ymm4	\n\t"\
		"vmovaps	%%ymm6,     (%%rsi)	\n\t"\
		"vmovaps	%%ymm5,0x0c0(%%rsi)	\n\t"\
		"vmovaps	%%ymm7,0x020(%%rsi)	\n\t"\
		"vmovaps	%%ymm4,0x060(%%rsi)	\n\t"\
		"/***************************************/	\n\t"\
		"movslq	%[__p04],%%rdi	\n\t"\
		"shlq	$3,%%rdi	\n\t"\
		"addq	%%rdi,%%rax	/* &a[j1+p5] */	\n\t"\
		"addq	%%rdi,%%rbx	\n\t"\
		"addq	%%rdi,%%rcx	\n\t"\
		"addq	%%rdi,%%rdx	\n\t"\
		"/* SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO_B(r28,c05) */	\n\t"\
		"/* rsi contains r20 */	\n\t"\
		"vmovaps	     (%%rax),%%ymm0	\n\t"\
		"vmovaps	     (%%rcx),%%ymm4	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5	\n\t"\
		"vmovaps	0x9e0(%%rsi),%%ymm6	/* c05 */	\n\t"\
		"vmovaps	0xa00(%%rsi),%%ymm7	\n\t"\
		"vmovaps	%%ymm0,%%ymm2	\n\t"\
		"vmovaps	%%ymm1,%%ymm3	\n\t"\
		"vmulpd	%%ymm6,%%ymm0,%%ymm0	\n\t"\
		"vmulpd	%%ymm6,%%ymm1,%%ymm1	\n\t"\
		"vmulpd	%%ymm7,%%ymm2,%%ymm2	\n\t"\
		"vmulpd	%%ymm7,%%ymm3,%%ymm3	\n\t"\
		"vmovaps	%%ymm4,%%ymm6	\n\t"\
		"vaddpd	%%ymm2,%%ymm1,%%ymm1	\n\t"\
		"vmovaps	%%ymm5,%%ymm7	\n\t"\
		"vmulpd	0xa20(%%rsi),%%ymm4,%%ymm4	/* c15 */	\n\t"\
		"vsubpd	%%ymm3,%%ymm0,%%ymm0	\n\t"\
		"vmulpd	0xa20(%%rsi),%%ymm5,%%ymm5	\n\t"\
		"vmulpd	0xa40(%%rsi),%%ymm6,%%ymm6	\n\t"\
		"vmovaps	%%ymm0,%%ymm2	\n\t"\
		"vmulpd	0xa40(%%rsi),%%ymm7,%%ymm7	\n\t"\
		"vaddpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
		"vmovaps	%%ymm1,%%ymm3	\n\t"\
		"vsubpd	%%ymm7,%%ymm4,%%ymm4	\n\t"\
		"vaddpd	%%ymm4,%%ymm0,%%ymm0	\n\t"\
		"vaddpd	%%ymm5,%%ymm1,%%ymm1	\n\t"\
		"vsubpd	%%ymm4,%%ymm2,%%ymm2	\n\t"\
		"vsubpd	%%ymm5,%%ymm3,%%ymm3	\n\t"\
		"vmovaps	     (%%rdx),%%ymm4	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm5	\n\t"\
		"vmovaps	     (%%rdx),%%ymm6	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm7	\n\t"\
		"vmulpd	0xaa0(%%rsi),%%ymm4,%%ymm4	/* c1D */	\n\t"\
		"vmulpd	0xaa0(%%rsi),%%ymm5,%%ymm5	\n\t"\
		"vmulpd	0xac0(%%rsi),%%ymm6,%%ymm6	\n\t"\
		"vmulpd	0xac0(%%rsi),%%ymm7,%%ymm7	\n\t"\
		"addq	$0x100,%%rsi	/* r28 */	\n\t"\
		"vaddpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
		"vsubpd	%%ymm7,%%ymm4,%%ymm4	\n\t"\
		"vmovaps	%%ymm5,0x020(%%rsi)	\n\t"\
		"vmovaps	%%ymm4,     (%%rsi)	\n\t"\
		"vmovaps	     (%%rbx),%%ymm4	\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm5	\n\t"\
		"vmovaps	     (%%rbx),%%ymm6	\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm7	\n\t"\
		"vmulpd	0x960(%%rsi),%%ymm4,%%ymm4	/* c0D */	\n\t"\
		"vmulpd	0x960(%%rsi),%%ymm5,%%ymm5	\n\t"\
		"vmulpd	0x980(%%rsi),%%ymm6,%%ymm6	\n\t"\
		"vmulpd	0x980(%%rsi),%%ymm7,%%ymm7	\n\t"\
		"vaddpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
		"vsubpd	%%ymm7,%%ymm4,%%ymm4	\n\t"\
		"vmovaps	%%ymm5,%%ymm7	\n\t"\
		"vmovaps	%%ymm4,%%ymm6	\n\t"\
		"vsubpd	     (%%rsi),%%ymm4,%%ymm4	\n\t"\
		"vsubpd	0x020(%%rsi),%%ymm5,%%ymm5	\n\t"\
		"vaddpd	     (%%rsi),%%ymm6,%%ymm6	\n\t"\
		"vaddpd	0x020(%%rsi),%%ymm7,%%ymm7	\n\t"\
		"vsubpd	%%ymm6,%%ymm0,%%ymm0	\n\t"\
		"vsubpd	%%ymm7,%%ymm1,%%ymm1	\n\t"\
		"vmovaps	%%ymm0,0x080(%%rsi)	\n\t"\
		"vsubpd	%%ymm5,%%ymm2,%%ymm2	\n\t"\
		"vmovaps	%%ymm1,0x0a0(%%rsi)	\n\t"\
		"vsubpd	%%ymm4,%%ymm3,%%ymm3	\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6	\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7	\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4	\n\t"\
		"vaddpd	%%ymm0,%%ymm6,%%ymm6	\n\t"\
		"vaddpd	%%ymm2,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm1,%%ymm7,%%ymm7	\n\t"\
		"vaddpd	%%ymm3,%%ymm4,%%ymm4	\n\t"\
		"vmovaps	%%ymm6,     (%%rsi)	\n\t"\
		"vmovaps	%%ymm7,0x020(%%rsi)	\n\t"\
		"vmovaps	0x300(%%rsi),%%ymm0	/* isrt2 */	\n\t"\
		"vmovaps	%%ymm2,%%ymm6	\n\t"\
		"vmovaps	%%ymm5,%%ymm7	\n\t"\
		"vsubpd	%%ymm4,%%ymm2,%%ymm2	\n\t"\
		"vsubpd	%%ymm3,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm4,%%ymm6,%%ymm6	\n\t"\
		"vaddpd	%%ymm3,%%ymm7,%%ymm7	\n\t"\
		"vmulpd	%%ymm0,%%ymm2,%%ymm2	\n\t"\
		"vmulpd	%%ymm0,%%ymm5,%%ymm5	\n\t"\
		"vmulpd	%%ymm0,%%ymm6,%%ymm6	\n\t"\
		"vmulpd	%%ymm0,%%ymm7,%%ymm7	\n\t"\
		"vmovaps	%%ymm2,0x040(%%rsi)	\n\t"\
		"vmovaps	%%ymm5,0x0c0(%%rsi)	\n\t"\
		"vmovaps	%%ymm6,0x060(%%rsi)	\n\t"\
		"vmovaps	%%ymm7,0x0e0(%%rsi)	\n\t"\
		"/***************************************/	\n\t"\
		"\n\t"\
		"/* SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r20) */	\n\t"\
		"subq		$0x100,%%rsi	/* r20 */	\n\t"\
		"vmovaps		     (%%rsi),%%ymm0	\n\t"\
		"vmovaps		0x080(%%rsi),%%ymm4	\n\t"\
		"vmovaps		0x020(%%rsi),%%ymm1	\n\t"\
		"vmovaps		0x0a0(%%rsi),%%ymm5	\n\t"\
		"vmovaps		0x100(%%rsi),%%ymm2	\n\t"\
		"vmovaps		0x1a0(%%rsi),%%ymm7	\n\t"\
		"vmovaps		0x120(%%rsi),%%ymm3	\n\t"\
		"vmovaps		0x180(%%rsi),%%ymm6	\n\t"\
		"vsubpd	%%ymm2,%%ymm0,%%ymm0	\n\t"\
		"vsubpd	%%ymm7,%%ymm4,%%ymm4	\n\t"\
		"vsubpd	%%ymm3,%%ymm1,%%ymm1	\n\t"\
		"vsubpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm2,%%ymm2,%%ymm2	\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7	\n\t"\
		"vaddpd	%%ymm3,%%ymm3,%%ymm3	\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6	\n\t"\
		"vaddpd	%%ymm0,%%ymm2,%%ymm2	\n\t"\
		"vaddpd	%%ymm4,%%ymm7,%%ymm7	\n\t"\
		"vaddpd	%%ymm1,%%ymm3,%%ymm3	\n\t"\
		"vaddpd	%%ymm5,%%ymm6,%%ymm6	\n\t"\
		"vmovaps		%%ymm0,0x100(%%rsi)	\n\t"\
		"vmovaps		%%ymm4,0x080(%%rsi)	\n\t"\
		"vmovaps		%%ymm1,0x120(%%rsi)	\n\t"\
		"vmovaps		%%ymm5,0x1a0(%%rsi)	\n\t"\
		"vmovaps		%%ymm2,     (%%rsi)	\n\t"\
		"vmovaps		%%ymm7,0x180(%%rsi)	\n\t"\
		"vmovaps		%%ymm3,0x020(%%rsi)	\n\t"\
		"vmovaps		%%ymm6,0x0a0(%%rsi)	\n\t"\
		"vmovaps		0x040(%%rsi),%%ymm0	\n\t"\
		"vmovaps		0x0c0(%%rsi),%%ymm4	\n\t"\
		"vmovaps		0x060(%%rsi),%%ymm1	\n\t"\
		"vmovaps		0x0e0(%%rsi),%%ymm5	\n\t"\
		"vmovaps		0x140(%%rsi),%%ymm2	\n\t"\
		"vmovaps		0x1e0(%%rsi),%%ymm7	\n\t"\
		"vmovaps		0x160(%%rsi),%%ymm3	\n\t"\
		"vmovaps		0x1c0(%%rsi),%%ymm6	\n\t"\
		"vsubpd	%%ymm2,%%ymm0,%%ymm0	\n\t"\
		"vsubpd	%%ymm7,%%ymm4,%%ymm4	\n\t"\
		"vsubpd	%%ymm3,%%ymm1,%%ymm1	\n\t"\
		"vsubpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm2,%%ymm2,%%ymm2	\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7	\n\t"\
		"vaddpd	%%ymm3,%%ymm3,%%ymm3	\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6	\n\t"\
		"vaddpd	%%ymm0,%%ymm2,%%ymm2	\n\t"\
		"vaddpd	%%ymm4,%%ymm7,%%ymm7	\n\t"\
		"vaddpd	%%ymm1,%%ymm3,%%ymm3	\n\t"\
		"vaddpd	%%ymm5,%%ymm6,%%ymm6	\n\t"\
		"vmovaps		%%ymm0,0x140(%%rsi)	\n\t"\
		"vmovaps		%%ymm4,0x0c0(%%rsi)	\n\t"\
		"vmovaps		%%ymm1,0x160(%%rsi)	\n\t"\
		"vmovaps		%%ymm5,0x1e0(%%rsi)	\n\t"\
		"vmovaps		%%ymm2,0x040(%%rsi)	\n\t"\
		"vmovaps		%%ymm7,0x1c0(%%rsi)	\n\t"\
		"vmovaps		%%ymm3,0x060(%%rsi)	\n\t"\
		"vmovaps		%%ymm6,0x0e0(%%rsi)	\n\t"\
		"/***************************************/	\n\t"\
		"\n\t"\
		"/*...Block 4: */	\n\t"\
		"movq	%[__add0],%%rax	\n\t"\
		"movslq	%[__p03],%%rdi	/* Do things this way [rather than repeatedly adding p1] since array-padding scheme means p2 == p+p1 not guaranteed. */	\n\t"\
		"movslq	%[__p08],%%rbx	/* Can't get these via simple load-one-and-shift-as-needed due to array padding scheme */	\n\t"\
		"movslq	%[__p10],%%rcx	\n\t"\
		"movslq	%[__p18],%%rdx	\n\t"\
		"shlq	$3,%%rdi	\n\t"\
		"shlq	$3,%%rbx	\n\t"\
		"shlq	$3,%%rcx	\n\t"\
		"shlq	$3,%%rdx	\n\t"\
		"addq	%%rdi,%%rax	/* &a[j1+p3] */	\n\t"\
		"addq	%%rax,%%rbx	\n\t"\
		"addq	%%rax,%%rcx	\n\t"\
		"addq	%%rax,%%rdx	\n\t"\
		"/* SSE2_RADIX4_DIF_4TWIDDLE_B(r30,c03) */	\n\t"\
		"addq	$0x200,%%rsi	/* r30 */	\n\t"\
		"vmovaps	     (%%rax),%%ymm0	\n\t"\
		"vmovaps	     (%%rcx),%%ymm4	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5	\n\t"\
		"vmovaps	0x8e0(%%rsi),%%ymm6	/* c03 */	\n\t"\
		"vmovaps	0x900(%%rsi),%%ymm7	\n\t"\
		"vmovaps	%%ymm0,%%ymm2	\n\t"\
		"vmovaps	%%ymm1,%%ymm3	\n\t"\
		"vmulpd   %%ymm6,%%ymm0,%%ymm0	\n\t"\
		"vmulpd   %%ymm6,%%ymm1,%%ymm1	\n\t"\
		"vmulpd   %%ymm7,%%ymm2,%%ymm2	\n\t"\
		"vmulpd   %%ymm7,%%ymm3,%%ymm3	\n\t"\
		"vmovaps	%%ymm4,%%ymm6	\n\t"\
		"vaddpd   %%ymm2,%%ymm1,%%ymm1	\n\t"\
		"vmovaps	%%ymm5,%%ymm7	\n\t"\
		"vmulpd   0x920(%%rsi),%%ymm4,%%ymm4	/* c13 */	\n\t"\
		"vsubpd   %%ymm3,%%ymm0,%%ymm0	\n\t"\
		"vmulpd   0x920(%%rsi),%%ymm5,%%ymm5	\n\t"\
		"vmulpd   0x940(%%rsi),%%ymm6,%%ymm6	\n\t"\
		"vmovaps	%%ymm0,%%ymm2	\n\t"\
		"vmulpd	0x940(%%rsi),%%ymm7,%%ymm7	\n\t"\
		"vaddpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
		"vmovaps	%%ymm1,%%ymm3	\n\t"\
		"vsubpd	%%ymm7,%%ymm4,%%ymm4	\n\t"\
		"vaddpd	%%ymm4,%%ymm0,%%ymm0	\n\t"\
		"vaddpd	%%ymm5,%%ymm1,%%ymm1	\n\t"\
		"vsubpd	%%ymm4,%%ymm2,%%ymm2	\n\t"\
		"vsubpd	%%ymm5,%%ymm3,%%ymm3	\n\t"\
		"vmovaps	     (%%rdx),%%ymm4	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm5	\n\t"\
		"vmovaps	     (%%rdx),%%ymm6	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm7	\n\t"\
		"vmulpd	0x9a0(%%rsi),%%ymm4,%%ymm4	/* c1B */	\n\t"\
		"vmulpd	0x9a0(%%rsi),%%ymm5,%%ymm5	\n\t"\
		"vmulpd	0x9c0(%%rsi),%%ymm6,%%ymm6	\n\t"\
		"vmulpd	0x9c0(%%rsi),%%ymm7,%%ymm7	\n\t"\
		"vaddpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
		"vsubpd	%%ymm7,%%ymm4,%%ymm4	\n\t"\
		"vmovaps	%%ymm5,0x020(%%rsi)	\n\t"\
		"vmovaps	%%ymm4,     (%%rsi)	\n\t"\
		"vmovaps	     (%%rbx),%%ymm4	\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm5	\n\t"\
		"vmovaps	     (%%rbx),%%ymm6	\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm7	\n\t"\
		"vmulpd	0x960(%%rsi),%%ymm4,%%ymm4	/* c0B */	\n\t"\
		"vmulpd	0x960(%%rsi),%%ymm5,%%ymm5	\n\t"\
		"vmulpd	0x980(%%rsi),%%ymm6,%%ymm6	\n\t"\
		"vmulpd	0x980(%%rsi),%%ymm7,%%ymm7	\n\t"\
		"vaddpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
		"vsubpd	%%ymm7,%%ymm4,%%ymm4	\n\t"\
		"vmovaps	%%ymm5,%%ymm7	\n\t"\
		"vmovaps	%%ymm4,%%ymm6	\n\t"\
		"vsubpd	     (%%rsi),%%ymm4,%%ymm4	\n\t"\
		"vsubpd	0x020(%%rsi),%%ymm5,%%ymm5	\n\t"\
		"vaddpd	     (%%rsi),%%ymm6,%%ymm6	\n\t"\
		"vaddpd	0x020(%%rsi),%%ymm7,%%ymm7	\n\t"\
		"vsubpd	%%ymm6,%%ymm0,%%ymm0	\n\t"\
		"vsubpd	%%ymm5,%%ymm2,%%ymm2	\n\t"\
		"vsubpd	%%ymm7,%%ymm1,%%ymm1	\n\t"\
		"vsubpd	%%ymm4,%%ymm3,%%ymm3	\n\t"\
		"vmovaps	%%ymm0,0x080(%%rsi)	\n\t"\
		"vmovaps	%%ymm2,0x040(%%rsi)	\n\t"\
		"vmovaps	%%ymm1,0x0a0(%%rsi)	\n\t"\
		"vmovaps	%%ymm3,0x0e0(%%rsi)	\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6	\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7	\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4	\n\t"\
		"vaddpd	%%ymm0,%%ymm6,%%ymm6	\n\t"\
		"vaddpd	%%ymm2,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm1,%%ymm7,%%ymm7	\n\t"\
		"vaddpd	%%ymm3,%%ymm4,%%ymm4	\n\t"\
		"vmovaps	%%ymm6,     (%%rsi)	\n\t"\
		"vmovaps	%%ymm5,0x0c0(%%rsi)	\n\t"\
		"vmovaps	%%ymm7,0x020(%%rsi)	\n\t"\
		"vmovaps	%%ymm4,0x060(%%rsi)	\n\t"\
		"/***************************************/	\n\t"\
		"movslq	%[__p04],%%rdi	\n\t"\
		"shlq	$3,%%rdi	\n\t"\
		"addq	%%rdi,%%rax	/* &a[j1+p7] */	\n\t"\
		"addq	%%rdi,%%rbx	\n\t"\
		"addq	%%rdi,%%rcx	\n\t"\
		"addq	%%rdi,%%rdx	\n\t"\
		"/* SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO_B(r38,c07) */	\n\t"\
		"/* rsi contains r30 */	\n\t"\
		"vmovaps	     (%%rax),%%ymm0	\n\t"\
		"vmovaps	     (%%rcx),%%ymm4	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5	\n\t"\
		"vmovaps	0x9e0(%%rsi),%%ymm6	/* c07 */	\n\t"\
		"vmovaps	0xa00(%%rsi),%%ymm7	\n\t"\
		"vmovaps	%%ymm0,%%ymm2	\n\t"\
		"vmovaps	%%ymm1,%%ymm3	\n\t"\
		"vmulpd	%%ymm6,%%ymm0,%%ymm0	\n\t"\
		"vmulpd	%%ymm6,%%ymm1,%%ymm1	\n\t"\
		"vmulpd	%%ymm7,%%ymm2,%%ymm2	\n\t"\
		"vmulpd	%%ymm7,%%ymm3,%%ymm3	\n\t"\
		"vmovaps	%%ymm4,%%ymm6	\n\t"\
		"vaddpd	%%ymm2,%%ymm1,%%ymm1	\n\t"\
		"vmovaps	%%ymm5,%%ymm7	\n\t"\
		"vmulpd	0xa20(%%rsi),%%ymm4,%%ymm4	/* c17 */	\n\t"\
		"vsubpd	%%ymm3,%%ymm0,%%ymm0	\n\t"\
		"vmulpd	0xa20(%%rsi),%%ymm5,%%ymm5	\n\t"\
		"vmulpd	0xa40(%%rsi),%%ymm6,%%ymm6	\n\t"\
		"vmovaps	%%ymm0,%%ymm2	\n\t"\
		"vmulpd	0xa40(%%rsi),%%ymm7,%%ymm7	\n\t"\
		"vaddpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
		"vmovaps	%%ymm1,%%ymm3	\n\t"\
		"vsubpd	%%ymm7,%%ymm4,%%ymm4	\n\t"\
		"vaddpd	%%ymm4,%%ymm0,%%ymm0	\n\t"\
		"vaddpd	%%ymm5,%%ymm1,%%ymm1	\n\t"\
		"vsubpd	%%ymm4,%%ymm2,%%ymm2	\n\t"\
		"vsubpd	%%ymm5,%%ymm3,%%ymm3	\n\t"\
		"vmovaps	     (%%rdx),%%ymm4	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm5	\n\t"\
		"vmovaps	     (%%rdx),%%ymm6	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm7	\n\t"\
		"vmulpd	0xaa0(%%rsi),%%ymm4,%%ymm4	/* c1F */	\n\t"\
		"vmulpd	0xaa0(%%rsi),%%ymm5,%%ymm5	\n\t"\
		"vmulpd	0xac0(%%rsi),%%ymm6,%%ymm6	\n\t"\
		"vmulpd	0xac0(%%rsi),%%ymm7,%%ymm7	\n\t"\
		"addq	$0x100,%%rsi	/* r38 */	\n\t"\
		"vaddpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
		"vsubpd	%%ymm7,%%ymm4,%%ymm4	\n\t"\
		"vmovaps	%%ymm5,0x020(%%rsi)	\n\t"\
		"vmovaps	%%ymm4,     (%%rsi)	\n\t"\
		"vmovaps	     (%%rbx),%%ymm4	\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm5	\n\t"\
		"vmovaps	     (%%rbx),%%ymm6	\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm7	\n\t"\
		"vmulpd	0x960(%%rsi),%%ymm4,%%ymm4	/* c0F */	\n\t"\
		"vmulpd	0x960(%%rsi),%%ymm5,%%ymm5	\n\t"\
		"vmulpd	0x980(%%rsi),%%ymm6,%%ymm6	\n\t"\
		"vmulpd	0x980(%%rsi),%%ymm7,%%ymm7	\n\t"\
		"vaddpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
		"vsubpd	%%ymm7,%%ymm4,%%ymm4	\n\t"\
		"vmovaps	%%ymm5,%%ymm7	\n\t"\
		"vmovaps	%%ymm4,%%ymm6	\n\t"\
		"vsubpd	     (%%rsi),%%ymm4,%%ymm4	\n\t"\
		"vsubpd	0x020(%%rsi),%%ymm5,%%ymm5	\n\t"\
		"vaddpd	     (%%rsi),%%ymm6,%%ymm6	\n\t"\
		"vaddpd	0x020(%%rsi),%%ymm7,%%ymm7	\n\t"\
		"vsubpd	%%ymm6,%%ymm0,%%ymm0	\n\t"\
		"vsubpd	%%ymm7,%%ymm1,%%ymm1	\n\t"\
		"vmovaps	%%ymm0,0x080(%%rsi)	\n\t"\
		"vsubpd	%%ymm5,%%ymm2,%%ymm2	\n\t"\
		"vmovaps	%%ymm1,0x0a0(%%rsi)	\n\t"\
		"vsubpd	%%ymm4,%%ymm3,%%ymm3	\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6	\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7	\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4	\n\t"\
		"vaddpd	%%ymm0,%%ymm6,%%ymm6	\n\t"\
		"vaddpd	%%ymm2,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm1,%%ymm7,%%ymm7	\n\t"\
		"vaddpd	%%ymm3,%%ymm4,%%ymm4	\n\t"\
		"vmovaps	%%ymm6,     (%%rsi)	\n\t"\
		"vmovaps	%%ymm7,0x020(%%rsi)	\n\t"\
		"vmovaps	0x100(%%rsi),%%ymm0	/* isrt2 */	\n\t"\
		"vmovaps	%%ymm2,%%ymm6	\n\t"\
		"vmovaps	%%ymm5,%%ymm7	\n\t"\
		"vsubpd	%%ymm4,%%ymm2,%%ymm2	\n\t"\
		"vsubpd	%%ymm3,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm4,%%ymm6,%%ymm6	\n\t"\
		"vaddpd	%%ymm3,%%ymm7,%%ymm7	\n\t"\
		"vmulpd	%%ymm0,%%ymm2,%%ymm2	\n\t"\
		"vmulpd	%%ymm0,%%ymm5,%%ymm5	\n\t"\
		"vmulpd	%%ymm0,%%ymm6,%%ymm6	\n\t"\
		"vmulpd	%%ymm0,%%ymm7,%%ymm7	\n\t"\
		"vmovaps	%%ymm2,0x040(%%rsi)	\n\t"\
		"vmovaps	%%ymm5,0x0c0(%%rsi)	\n\t"\
		"vmovaps	%%ymm6,0x060(%%rsi)	\n\t"\
		"vmovaps	%%ymm7,0x0e0(%%rsi)	\n\t"\
		"/***************************************/	\n\t"\
		"\n\t"\
		"/* SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r30) */	\n\t"\
		"subq		$0x100,%%rsi	/* r30 */	\n\t"\
		"vmovaps		     (%%rsi),%%ymm0	\n\t"\
		"vmovaps		0x080(%%rsi),%%ymm4	\n\t"\
		"vmovaps		0x020(%%rsi),%%ymm1	\n\t"\
		"vmovaps		0x0a0(%%rsi),%%ymm5	\n\t"\
		"vmovaps		0x100(%%rsi),%%ymm2	\n\t"\
		"vmovaps		0x1a0(%%rsi),%%ymm7	\n\t"\
		"vmovaps		0x120(%%rsi),%%ymm3	\n\t"\
		"vmovaps		0x180(%%rsi),%%ymm6	\n\t"\
		"vsubpd	%%ymm2,%%ymm0,%%ymm0	\n\t"\
		"vsubpd	%%ymm7,%%ymm4,%%ymm4	\n\t"\
		"vsubpd	%%ymm3,%%ymm1,%%ymm1	\n\t"\
		"vsubpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm2,%%ymm2,%%ymm2	\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7	\n\t"\
		"vaddpd	%%ymm3,%%ymm3,%%ymm3	\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6	\n\t"\
		"vaddpd	%%ymm0,%%ymm2,%%ymm2	\n\t"\
		"vaddpd	%%ymm4,%%ymm7,%%ymm7	\n\t"\
		"vaddpd	%%ymm1,%%ymm3,%%ymm3	\n\t"\
		"vaddpd	%%ymm5,%%ymm6,%%ymm6	\n\t"\
		"vmovaps		%%ymm0,0x100(%%rsi)	\n\t"\
		"vmovaps		%%ymm4,0x080(%%rsi)	\n\t"\
		"vmovaps		%%ymm1,0x120(%%rsi)	\n\t"\
		"vmovaps		%%ymm5,0x1a0(%%rsi)	\n\t"\
		"vmovaps		%%ymm2,     (%%rsi)	\n\t"\
		"vmovaps		%%ymm7,0x180(%%rsi)	\n\t"\
		"vmovaps		%%ymm3,0x020(%%rsi)	\n\t"\
		"vmovaps		%%ymm6,0x0a0(%%rsi)	\n\t"\
		"vmovaps		0x040(%%rsi),%%ymm0	\n\t"\
		"vmovaps		0x0c0(%%rsi),%%ymm4	\n\t"\
		"vmovaps		0x060(%%rsi),%%ymm1	\n\t"\
		"vmovaps		0x0e0(%%rsi),%%ymm5	\n\t"\
		"vmovaps		0x140(%%rsi),%%ymm2	\n\t"\
		"vmovaps		0x1e0(%%rsi),%%ymm7	\n\t"\
		"vmovaps		0x160(%%rsi),%%ymm3	\n\t"\
		"vmovaps		0x1c0(%%rsi),%%ymm6	\n\t"\
		"vsubpd	%%ymm2,%%ymm0,%%ymm0	\n\t"\
		"vsubpd	%%ymm7,%%ymm4,%%ymm4	\n\t"\
		"vsubpd	%%ymm3,%%ymm1,%%ymm1	\n\t"\
		"vsubpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm2,%%ymm2,%%ymm2	\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7	\n\t"\
		"vaddpd	%%ymm3,%%ymm3,%%ymm3	\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6	\n\t"\
		"vaddpd	%%ymm0,%%ymm2,%%ymm2	\n\t"\
		"vaddpd	%%ymm4,%%ymm7,%%ymm7	\n\t"\
		"vaddpd	%%ymm1,%%ymm3,%%ymm3	\n\t"\
		"vaddpd	%%ymm5,%%ymm6,%%ymm6	\n\t"\
		"vmovaps		%%ymm0,0x140(%%rsi)	\n\t"\
		"vmovaps		%%ymm4,0x0c0(%%rsi)	\n\t"\
		"vmovaps		%%ymm1,0x160(%%rsi)	\n\t"\
		"vmovaps		%%ymm5,0x1e0(%%rsi)	\n\t"\
		"vmovaps		%%ymm2,0x040(%%rsi)	\n\t"\
		"vmovaps		%%ymm7,0x1c0(%%rsi)	\n\t"\
		"vmovaps		%%ymm3,0x060(%%rsi)	\n\t"\
		"vmovaps		%%ymm6,0x0e0(%%rsi)	\n\t"\
		"\n\t"\
		"/**********************************************************************************/	\n\t"\
		"/*...and now do eight radix-4 transforms, including the internal twiddle factors: */	\n\t"\
		"/**********************************************************************************/	\n\t"\
		"/*...Block 1: t00,t10,t20,t30	*/	\n\t"\
		"movq	%[__add0],%%rax	/* &a[j1] */	\n\t"\
		"movslq	%[__p01],%%rbx	\n\t"\
		"movslq	%[__p02],%%rcx	\n\t"\
		"movslq	%[__p03],%%rdx	\n\t"\
		"shlq	$3,%%rbx	\n\t"\
		"shlq	$3,%%rcx	\n\t"\
		"shlq	$3,%%rdx	\n\t"\
		"addq	%%rax,%%rbx	\n\t"\
		"addq	%%rax,%%rcx	\n\t"\
		"addq	%%rax,%%rdx	\n\t"\
		"movq	%[__r00],%%rsi	\n\t"\
		"vmovaps	     (%%rsi),%%ymm0	\n\t"\
		"vmovaps	0x400(%%rsi),%%ymm4	\n\t"\
		"vmovaps	0x020(%%rsi),%%ymm1	\n\t"\
		"vmovaps	0x420(%%rsi),%%ymm5	\n\t"\
		"vmovaps	0x200(%%rsi),%%ymm2	\n\t"\
		"vmovaps	0x600(%%rsi),%%ymm6	\n\t"\
		"vmovaps	0x220(%%rsi),%%ymm3	\n\t"\
		"vmovaps	0x620(%%rsi),%%ymm7	\n\t"\
		"vsubpd	0x200(%%rsi),%%ymm0,%%ymm0	\n\t"\
		"vsubpd	0x600(%%rsi),%%ymm4,%%ymm4	\n\t"\
		"vsubpd	0x220(%%rsi),%%ymm1,%%ymm1	\n\t"\
		"vsubpd	0x620(%%rsi),%%ymm5,%%ymm5	\n\t"\
		"vaddpd	     (%%rsi),%%ymm2,%%ymm2	\n\t"\
		"vaddpd	0x400(%%rsi),%%ymm6,%%ymm6	\n\t"\
		"vaddpd	0x020(%%rsi),%%ymm3,%%ymm3	\n\t"\
		"vaddpd	0x420(%%rsi),%%ymm7,%%ymm7	\n\t"\
		"vsubpd	%%ymm6,%%ymm2,%%ymm2	\n\t"\
		"vsubpd	%%ymm5,%%ymm0,%%ymm0	\n\t"\
		"vsubpd	%%ymm7,%%ymm3,%%ymm3	\n\t"\
		"vsubpd	%%ymm4,%%ymm1,%%ymm1	\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6	\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7	\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4	\n\t"\
		"vmovaps	%%ymm2,     (%%rbx)	\n\t"\
		"vmovaps	%%ymm0,     (%%rcx)	\n\t"\
		"vmovaps	%%ymm3,0x020(%%rbx)	\n\t"\
		"vmovaps	%%ymm1,0x020(%%rdx)	\n\t"\
		"vaddpd	%%ymm2,%%ymm6,%%ymm6	\n\t"\
		"vaddpd	%%ymm0,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm3,%%ymm7,%%ymm7	\n\t"\
		"vaddpd	%%ymm1,%%ymm4,%%ymm4	\n\t"\
		"vmovaps	%%ymm6,     (%%rax)	\n\t"\
		"vmovaps	%%ymm5,     (%%rdx)	\n\t"\
		"vmovaps	%%ymm7,0x020(%%rax)	\n\t"\
		"vmovaps	%%ymm4,0x020(%%rcx)	\n\t"\
		"\n\t"\
		"/*...Block 5: t08,t18,t28,t38	*/	\n\t"\
		"addq	$0x100,%%rsi	/* r08 */	\n\t"\
		"subq	%%rax,%%rbx	/* p01 << 3 */	\n\t"\
		"subq	%%rax,%%rcx	/* p02 << 3 */	\n\t"\
		"subq	%%rax,%%rdx	/* p03 << 3 */	\n\t"\
		"movslq	%[__p04],%%rdi	\n\t"\
		"shlq	$3,%%rdi	\n\t"\
		"addq	%%rdi,%%rax	/* &a[j1+p04] */	\n\t"\
		"addq	%%rax,%%rbx	\n\t"\
		"addq	%%rax,%%rcx	\n\t"\
		"addq	%%rax,%%rdx	\n\t"\
		"vmovaps	0x700(%%rsi),%%ymm3	/* isrt2 */	\n\t"\
		"vmovaps	0x400(%%rsi),%%ymm4	\n\t"\
		"vmovaps	0x420(%%rsi),%%ymm5	\n\t"\
		"vmovaps	0x600(%%rsi),%%ymm6	\n\t"\
		"vmovaps	0x620(%%rsi),%%ymm7	\n\t"\
		"vmulpd	%%ymm3,%%ymm4,%%ymm4	\n\t"\
		"vmovaps	     (%%rsi),%%ymm0	\n\t"\
		"vmulpd	%%ymm3,%%ymm5,%%ymm5	\n\t"\
		"vmovaps	0x020(%%rsi),%%ymm1	\n\t"\
		"vmulpd	%%ymm3,%%ymm6,%%ymm6	\n\t"\
		"vmovaps	0x200(%%rsi),%%ymm2	\n\t"\
		"vmulpd	%%ymm3,%%ymm7,%%ymm7	\n\t"\
		"vmovaps	0x220(%%rsi),%%ymm3	\n\t"\
		"vsubpd	%%ymm3,%%ymm0,%%ymm0	\n\t"\
		"vsubpd	%%ymm5,%%ymm4,%%ymm4	\n\t"\
		"vsubpd	%%ymm2,%%ymm1,%%ymm1	\n\t"\
		"vsubpd	%%ymm6,%%ymm7,%%ymm7	\n\t"\
		"vaddpd	%%ymm3,%%ymm3,%%ymm3	\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm2,%%ymm2,%%ymm2	\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6	\n\t"\
		"vaddpd	%%ymm0,%%ymm3,%%ymm3	\n\t"\
		"vaddpd	%%ymm4,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm1,%%ymm2,%%ymm2	\n\t"\
		"vaddpd	%%ymm7,%%ymm6,%%ymm6	\n\t"\
		"vsubpd	%%ymm6,%%ymm4,%%ymm4	\n\t"\
		"vsubpd	%%ymm7,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6	\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7	\n\t"\
		"vaddpd	%%ymm4,%%ymm6,%%ymm6	\n\t"\
		"vaddpd	%%ymm5,%%ymm7,%%ymm7	\n\t"\
		"vsubpd	%%ymm4,%%ymm0,%%ymm0	\n\t"\
		"vsubpd	%%ymm7,%%ymm3,%%ymm3	\n\t"\
		"vsubpd	%%ymm5,%%ymm2,%%ymm2	\n\t"\
		"vsubpd	%%ymm6,%%ymm1,%%ymm1	\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4	\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7	\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6	\n\t"\
		"vmovaps	%%ymm0,     (%%rbx)	\n\t"\
		"vmovaps	%%ymm3,     (%%rcx)	\n\t"\
		"vmovaps	%%ymm2,0x020(%%rbx)	\n\t"\
		"vmovaps	%%ymm1,0x020(%%rdx)	\n\t"\
		"vaddpd	%%ymm0,%%ymm4,%%ymm4	\n\t"\
		"vaddpd	%%ymm3,%%ymm7,%%ymm7	\n\t"\
		"vaddpd	%%ymm2,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm1,%%ymm6,%%ymm6	\n\t"\
		"vmovaps	%%ymm4,     (%%rax)	\n\t"\
		"vmovaps	%%ymm7,     (%%rdx)	\n\t"\
		"vmovaps	%%ymm5,0x020(%%rax)	\n\t"\
		"vmovaps	%%ymm6,0x020(%%rcx)	\n\t"\
		"\n\t"\
		"/*...Block 3: t04,t14,t24,t34	*/	\n\t"\
		"subq	$0x80,%%rsi	/* r04 */	\n\t"\
		"subq	%%rax,%%rbx	/* p01 << 3 */	\n\t"\
		"subq	%%rax,%%rcx	/* p02 << 3 */	\n\t"\
		"subq	%%rax,%%rdx	/* p03 << 3 */	\n\t"\
		"subq	%%rdi,%%rax	/* &a[j1] */	\n\t"\
		"movslq	%[__p08],%%rdi	\n\t"\
		"shlq	$3,%%rdi	\n\t"\
		"addq	%%rdi,%%rax	/* &a[j1+p08] */	\n\t"\
		"addq	%%rax,%%rbx	\n\t"\
		"addq	%%rax,%%rcx	\n\t"\
		"addq	%%rax,%%rdx	\n\t"\
		"vmovaps	0x400(%%rsi),%%ymm4	\n\t"\
		"vmovaps	0x420(%%rsi),%%ymm5	\n\t"\
		"vmovaps	0x7a0(%%rsi),%%ymm3	/* cc0 */	\n\t"\
		"vmovaps	0x7c0(%%rsi),%%ymm2	\n\t"\
		"vmovaps	%%ymm4,%%ymm6	\n\t"\
		"vmovaps	%%ymm5,%%ymm7	\n\t"\
		"vmulpd	%%ymm3,%%ymm4,%%ymm4	\n\t"\
		"vmulpd	%%ymm3,%%ymm5,%%ymm5	\n\t"\
		"vmulpd	%%ymm2,%%ymm6,%%ymm6	\n\t"\
		"vmovaps	0x600(%%rsi),%%ymm0	\n\t"\
		"vmulpd	%%ymm2,%%ymm7,%%ymm7	\n\t"\
		"vmovaps	0x620(%%rsi),%%ymm1	\n\t"\
		"vaddpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
		"vmovaps	%%ymm0,%%ymm6	\n\t"\
		"vsubpd	%%ymm7,%%ymm4,%%ymm4	\n\t"\
		"vmovaps	%%ymm1,%%ymm7	\n\t"\
		"vmulpd	%%ymm2,%%ymm6,%%ymm6	\n\t"\
		"vmulpd	%%ymm2,%%ymm7,%%ymm7	\n\t"\
		"vmulpd	%%ymm3,%%ymm0,%%ymm0	\n\t"\
		"vmulpd	%%ymm3,%%ymm1,%%ymm1	\n\t"\
		"vaddpd	%%ymm0,%%ymm7,%%ymm7	\n\t"\
		"vsubpd	%%ymm1,%%ymm6,%%ymm6	\n\t"\
		"vmovaps	%%ymm4,%%ymm2	\n\t"\
		"vmovaps	%%ymm5,%%ymm3	\n\t"\
		"vsubpd	%%ymm6,%%ymm4,%%ymm4	\n\t"\
		"vsubpd	%%ymm7,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm2,%%ymm6,%%ymm6	\n\t"\
		"vaddpd	%%ymm3,%%ymm7,%%ymm7	\n\t"\
		"vmovaps	0x200(%%rsi),%%ymm2	\n\t"\
		"vmovaps	0x220(%%rsi),%%ymm3	\n\t"\
		"vmovaps	0x780(%%rsi),%%ymm1	/* isrt2 */	\n\t"\
		"vmovaps	%%ymm2,%%ymm0	\n\t"\
		"vsubpd	%%ymm3,%%ymm2,%%ymm2	\n\t"\
		"vaddpd	%%ymm0,%%ymm3,%%ymm3	\n\t"\
		"vmulpd	%%ymm1,%%ymm2,%%ymm2	\n\t"\
		"vmulpd	%%ymm1,%%ymm3,%%ymm3	\n\t"\
		"vmovaps	     (%%rsi),%%ymm0	\n\t"\
		"vmovaps	0x020(%%rsi),%%ymm1	\n\t"\
		"vsubpd	%%ymm2,%%ymm0,%%ymm0	\n\t"\
		"vsubpd	%%ymm3,%%ymm1,%%ymm1	\n\t"\
		"vaddpd	%%ymm2,%%ymm2,%%ymm2	\n\t"\
		"vaddpd	%%ymm3,%%ymm3,%%ymm3	\n\t"\
		"vaddpd	%%ymm0,%%ymm2,%%ymm2	\n\t"\
		"vaddpd	%%ymm1,%%ymm3,%%ymm3	\n\t"\
		"vsubpd	%%ymm6,%%ymm2,%%ymm2	\n\t"\
		"vsubpd	%%ymm5,%%ymm0,%%ymm0	\n\t"\
		"vsubpd	%%ymm7,%%ymm3,%%ymm3	\n\t"\
		"vsubpd	%%ymm4,%%ymm1,%%ymm1	\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6	\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7	\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4	\n\t"\
		"vmovaps	%%ymm2,     (%%rbx)	\n\t"\
		"vmovaps	%%ymm0,     (%%rcx)	\n\t"\
		"vmovaps	%%ymm3,0x020(%%rbx)	\n\t"\
		"vmovaps	%%ymm1,0x020(%%rdx)	\n\t"\
		"vaddpd	%%ymm2,%%ymm6,%%ymm6	\n\t"\
		"vaddpd	%%ymm0,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm3,%%ymm7,%%ymm7	\n\t"\
		"vaddpd	%%ymm1,%%ymm4,%%ymm4	\n\t"\
		"vmovaps	%%ymm6,     (%%rax)	\n\t"\
		"vmovaps	%%ymm5,     (%%rdx)	\n\t"\
		"vmovaps	%%ymm7,0x020(%%rax)	\n\t"\
		"vmovaps	%%ymm4,0x020(%%rcx)	\n\t"\
		"\n\t"\
		"/*...Block 7: t0C,t1C,t2C,t3C	*/	\n\t"\
		"addq	$0x100,%%rsi	/* r0C */	\n\t"\
		"subq	%%rax,%%rbx	/* p01 << 3 */	\n\t"\
		"subq	%%rax,%%rcx	/* p02 << 3 */	\n\t"\
		"subq	%%rax,%%rdx	/* p03 << 3 */	\n\t"\
		"subq	%%rdi,%%rax	/* &a[j1] */	\n\t"\
		"movslq	%[__p0C],%%rdi	\n\t"\
		"shlq	$3,%%rdi	\n\t"\
		"addq	%%rdi,%%rax	/* &a[j1+p0C] */	\n\t"\
		"addq	%%rax,%%rbx	\n\t"\
		"addq	%%rax,%%rcx	\n\t"\
		"addq	%%rax,%%rdx	\n\t"\
		"vmovaps	0x400(%%rsi),%%ymm4	\n\t"\
		"vmovaps	0x420(%%rsi),%%ymm5	\n\t"\
		"vmovaps	0x6a0(%%rsi),%%ymm2	/* cc0 */	\n\t"\
		"vmovaps	0x6c0(%%rsi),%%ymm3	\n\t"\
		"vmovaps	%%ymm4,%%ymm6	\n\t"\
		"vmovaps	%%ymm5,%%ymm7	\n\t"\
		"vmulpd	%%ymm3,%%ymm4,%%ymm4	\n\t"\
		"vmulpd	%%ymm3,%%ymm5,%%ymm5	\n\t"\
		"vmulpd	%%ymm2,%%ymm6,%%ymm6	\n\t"\
		"vmovaps	0x600(%%rsi),%%ymm0	\n\t"\
		"vmulpd	%%ymm2,%%ymm7,%%ymm7	\n\t"\
		"vmovaps	0x620(%%rsi),%%ymm1	\n\t"\
		"vaddpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
		"vmovaps	%%ymm0,%%ymm6	\n\t"\
		"vsubpd	%%ymm7,%%ymm4,%%ymm4	\n\t"\
		"vmovaps	%%ymm1,%%ymm7	\n\t"\
		"vmulpd	%%ymm2,%%ymm6,%%ymm6	\n\t"\
		"vmulpd	%%ymm2,%%ymm7,%%ymm7	\n\t"\
		"vmulpd	%%ymm3,%%ymm0,%%ymm0	\n\t"\
		"vmulpd	%%ymm3,%%ymm1,%%ymm1	\n\t"\
		"vaddpd	%%ymm0,%%ymm7,%%ymm7	\n\t"\
		"vsubpd	%%ymm1,%%ymm6,%%ymm6	\n\t"\
		"vmovaps	%%ymm4,%%ymm2	\n\t"\
		"vmovaps	%%ymm5,%%ymm3	\n\t"\
		"vsubpd	%%ymm6,%%ymm4,%%ymm4	\n\t"\
		"vsubpd	%%ymm7,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm2,%%ymm6,%%ymm6	\n\t"\
		"vaddpd	%%ymm3,%%ymm7,%%ymm7	\n\t"\
		"vmovaps	0x200(%%rsi),%%ymm2	\n\t"\
		"vmovaps	0x220(%%rsi),%%ymm3	\n\t"\
		"vmovaps	0x680(%%rsi),%%ymm1	/* isrt2 */	\n\t"\
		"vmovaps	%%ymm2,%%ymm0	\n\t"\
		"vaddpd	%%ymm3,%%ymm2,%%ymm2	\n\t"\
		"vsubpd	%%ymm0,%%ymm3,%%ymm3	\n\t"\
		"vmulpd	%%ymm1,%%ymm2,%%ymm2	\n\t"\
		"vmulpd	%%ymm1,%%ymm3,%%ymm3	\n\t"\
		"vmovaps	     (%%rsi),%%ymm0	\n\t"\
		"vmovaps	0x020(%%rsi),%%ymm1	\n\t"\
		"vsubpd	%%ymm2,%%ymm0,%%ymm0	\n\t"\
		"vsubpd	%%ymm3,%%ymm1,%%ymm1	\n\t"\
		"vaddpd	%%ymm2,%%ymm2,%%ymm2	\n\t"\
		"vaddpd	%%ymm3,%%ymm3,%%ymm3	\n\t"\
		"vaddpd	%%ymm0,%%ymm2,%%ymm2	\n\t"\
		"vaddpd	%%ymm1,%%ymm3,%%ymm3	\n\t"\
		"vsubpd	%%ymm4,%%ymm0,%%ymm0	\n\t"\
		"vsubpd	%%ymm7,%%ymm2,%%ymm2	\n\t"\
		"vsubpd	%%ymm5,%%ymm1,%%ymm1	\n\t"\
		"vsubpd	%%ymm6,%%ymm3,%%ymm3	\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4	\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7	\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6	\n\t"\
		"vmovaps	%%ymm0,     (%%rbx)	\n\t"\
		"vmovaps	%%ymm2,     (%%rcx)	\n\t"\
		"vmovaps	%%ymm1,0x020(%%rbx)	\n\t"\
		"vmovaps	%%ymm3,0x020(%%rdx)	\n\t"\
		"vaddpd	%%ymm0,%%ymm4,%%ymm4	\n\t"\
		"vaddpd	%%ymm2,%%ymm7,%%ymm7	\n\t"\
		"vaddpd	%%ymm1,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm3,%%ymm6,%%ymm6	\n\t"\
		"vmovaps	%%ymm4,     (%%rax)	\n\t"\
		"vmovaps	%%ymm7,     (%%rdx)	\n\t"\
		"vmovaps	%%ymm5,0x020(%%rax)	\n\t"\
		"vmovaps	%%ymm6,0x020(%%rcx)	\n\t"\
		"\n\t"\
		"/*...Block 2: t02,t12,t22,t32	*/	\n\t"\
		"subq	$0x140,%%rsi	/* r02 */	\n\t"\
		"movslq	%[__p10],%%rdi	\n\t"\
		"subq	%%rax,%%rbx	\n\t"\
		"subq	%%rax,%%rcx	\n\t"\
		"subq	%%rax,%%rdx	\n\t"\
		"movq	%[__add0],%%rax	/* &a[j1] */	\n\t"\
		"shlq	$3,%%rdi	\n\t"\
		"addq	%%rdi,%%rax	/* &a[j1+p10) */	\n\t"\
		"addq	%%rax,%%rbx	\n\t"\
		"addq	%%rax,%%rcx	\n\t"\
		"addq	%%rax,%%rdx	\n\t"\
		"vmovaps	0x400(%%rsi),%%ymm4	\n\t"\
		"vmovaps	0x420(%%rsi),%%ymm5	\n\t"\
		"vmovaps	0x820(%%rsi),%%ymm2	/* cc1 */	\n\t"\
		"vmovaps	0x840(%%rsi),%%ymm3	\n\t"\
		"vmovaps	%%ymm4,%%ymm6	\n\t"\
		"vmovaps	%%ymm5,%%ymm7	\n\t"\
		"vmulpd	%%ymm2,%%ymm4,%%ymm4	\n\t"\
		"vmulpd	%%ymm2,%%ymm5,%%ymm5	\n\t"\
		"vmulpd	%%ymm3,%%ymm6,%%ymm6	\n\t"\
		"vmovaps	0x600(%%rsi),%%ymm0	\n\t"\
		"vmulpd	%%ymm3,%%ymm7,%%ymm7	\n\t"\
		"vmovaps	0x620(%%rsi),%%ymm1	\n\t"\
		"vmovaps	0x860(%%rsi),%%ymm2	/* cc3 */	\n\t"\
		"vmovaps	0x880(%%rsi),%%ymm3	\n\t"\
		"vaddpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
		"vmovaps	%%ymm0,%%ymm6	\n\t"\
		"vsubpd	%%ymm7,%%ymm4,%%ymm4	\n\t"\
		"vmovaps	%%ymm1,%%ymm7	\n\t"\
		"vmulpd	%%ymm2,%%ymm6,%%ymm6	\n\t"\
		"vmulpd	%%ymm2,%%ymm7,%%ymm7	\n\t"\
		"vmulpd	%%ymm3,%%ymm0,%%ymm0	\n\t"\
		"vmulpd	%%ymm3,%%ymm1,%%ymm1	\n\t"\
		"vaddpd	%%ymm0,%%ymm7,%%ymm7	\n\t"\
		"vsubpd	%%ymm1,%%ymm6,%%ymm6	\n\t"\
		"vmovaps	%%ymm4,%%ymm2	\n\t"\
		"vmovaps	%%ymm5,%%ymm3	\n\t"\
		"vsubpd	%%ymm6,%%ymm4,%%ymm4	\n\t"\
		"vsubpd	%%ymm7,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm2,%%ymm6,%%ymm6	\n\t"\
		"vaddpd	%%ymm3,%%ymm7,%%ymm7	\n\t"\
		"vmovaps	0x200(%%rsi),%%ymm1	\n\t"\
		"vmovaps	0x220(%%rsi),%%ymm3	\n\t"\
		"vmovaps	0x800(%%rsi),%%ymm0	/* ss1 */	\n\t"\
		"vmovaps	%%ymm1,%%ymm2	\n\t"\
		"vmulpd	%%ymm0,%%ymm1,%%ymm1	\n\t"\
		"vmulpd	%%ymm3,%%ymm0,%%ymm0	\n\t"\
		"vmulpd	0x7e0(%%rsi),%%ymm2,%%ymm2		/* cc1 */	\n\t"\
		"vmulpd	0x7e0(%%rsi),%%ymm3,%%ymm3	\n\t"\
		"vsubpd	%%ymm0,%%ymm2,%%ymm2	\n\t"\
		"vaddpd	%%ymm1,%%ymm3,%%ymm3	\n\t"\
		"vmovaps	     (%%rsi),%%ymm0	\n\t"\
		"vmovaps	0x020(%%rsi),%%ymm1	\n\t"\
		"vsubpd	%%ymm2,%%ymm0,%%ymm0	\n\t"\
		"vsubpd	%%ymm3,%%ymm1,%%ymm1	\n\t"\
		"vaddpd	%%ymm2,%%ymm2,%%ymm2	\n\t"\
		"vaddpd	%%ymm3,%%ymm3,%%ymm3	\n\t"\
		"vaddpd	%%ymm0,%%ymm2,%%ymm2	\n\t"\
		"vaddpd	%%ymm1,%%ymm3,%%ymm3	\n\t"\
		"vsubpd	%%ymm6,%%ymm2,%%ymm2	\n\t"\
		"vsubpd	%%ymm5,%%ymm0,%%ymm0	\n\t"\
		"vsubpd	%%ymm7,%%ymm3,%%ymm3	\n\t"\
		"vsubpd	%%ymm4,%%ymm1,%%ymm1	\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6	\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7	\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4	\n\t"\
		"vmovaps	%%ymm2,     (%%rbx)	\n\t"\
		"vmovaps	%%ymm0,     (%%rcx)	\n\t"\
		"vmovaps	%%ymm3,0x020(%%rbx)	\n\t"\
		"vmovaps	%%ymm1,0x020(%%rdx)	\n\t"\
		"vaddpd	%%ymm2,%%ymm6,%%ymm6	\n\t"\
		"vaddpd	%%ymm0,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm3,%%ymm7,%%ymm7	\n\t"\
		"vaddpd	%%ymm1,%%ymm4,%%ymm4	\n\t"\
		"vmovaps	%%ymm6,     (%%rax)	\n\t"\
		"vmovaps	%%ymm5,     (%%rdx)	\n\t"\
		"vmovaps	%%ymm7,0x020(%%rax)	\n\t"\
		"vmovaps	%%ymm4,0x020(%%rcx)	\n\t"\
		"\n\t"\
		"/*...Block 6: t0A,t1A,t2A,t3A	*/	\n\t"\
		"addq	$0x100,%%rsi	/* r0A */	\n\t"\
		"movslq	%[__p04],%%rdi	\n\t"\
		"shlq	$3,%%rdi	\n\t"\
		"addq	%%rdi,%%rax	/* &a[j1+p14] */	\n\t"\
		"addq	%%rdi,%%rbx	\n\t"\
		"addq	%%rdi,%%rcx	\n\t"\
		"addq	%%rdi,%%rdx	\n\t"\
		"vmovaps	0x400(%%rsi),%%ymm4	\n\t"\
		"vmovaps	0x420(%%rsi),%%ymm5	\n\t"\
		"vmovaps	0x760(%%rsi),%%ymm3	/* cc3 */	\n\t"\
		"vmovaps	0x780(%%rsi),%%ymm2	\n\t"\
		"vmovaps	%%ymm4,%%ymm6	\n\t"\
		"vmovaps	%%ymm5,%%ymm7	\n\t"\
		"vmulpd	%%ymm2,%%ymm4,%%ymm4	\n\t"\
		"vmulpd	%%ymm2,%%ymm5,%%ymm5	\n\t"\
		"vmulpd	%%ymm3,%%ymm6,%%ymm6	\n\t"\
		"vmovaps	0x600(%%rsi),%%ymm0	\n\t"\
		"vmulpd	%%ymm3,%%ymm7,%%ymm7	\n\t"\
		"vmovaps	0x620(%%rsi),%%ymm1	\n\t"\
		"vmovaps	0x720(%%rsi),%%ymm2	/* cc1 */	\n\t"\
		"vmovaps	0x740(%%rsi),%%ymm3	\n\t"\
		"vaddpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
		"vmovaps	%%ymm0,%%ymm6	\n\t"\
		"vsubpd	%%ymm7,%%ymm4,%%ymm4	\n\t"\
		"vmovaps	%%ymm1,%%ymm7	\n\t"\
		"vmulpd	%%ymm2,%%ymm6,%%ymm6	\n\t"\
		"vmulpd	%%ymm2,%%ymm7,%%ymm7	\n\t"\
		"vmulpd	%%ymm3,%%ymm0,%%ymm0	\n\t"\
		"vmulpd	%%ymm3,%%ymm1,%%ymm1	\n\t"\
		"vsubpd	%%ymm0,%%ymm7,%%ymm7	\n\t"\
		"vaddpd	%%ymm1,%%ymm6,%%ymm6	\n\t"\
		"vmovaps	%%ymm4,%%ymm2	\n\t"\
		"vmovaps	%%ymm5,%%ymm3	\n\t"\
		"vsubpd	%%ymm6,%%ymm4,%%ymm4	\n\t"\
		"vsubpd	%%ymm7,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm2,%%ymm6,%%ymm6	\n\t"\
		"vaddpd	%%ymm3,%%ymm7,%%ymm7	\n\t"\
		"vmovaps	0x200(%%rsi),%%ymm1	\n\t"\
		"vmovaps	0x220(%%rsi),%%ymm3	\n\t"\
		"vmovaps	0x6e0(%%rsi),%%ymm0		/* cc0 */	\n\t"\
		"vmovaps	%%ymm1,%%ymm2	\n\t"\
		"vmulpd	%%ymm0,%%ymm1,%%ymm1	\n\t"\
		"vmulpd	%%ymm3,%%ymm0,%%ymm0	\n\t"\
		"vmulpd	0x700(%%rsi),%%ymm2,%%ymm2	/* ss0 */	\n\t"\
		"vmulpd	0x700(%%rsi),%%ymm3,%%ymm3	\n\t"\
		"vaddpd	%%ymm0,%%ymm2,%%ymm2	\n\t"\
		"vsubpd	%%ymm1,%%ymm3,%%ymm3	\n\t"\
		"vmovaps	     (%%rsi),%%ymm0	\n\t"\
		"vmovaps	0x020(%%rsi),%%ymm1	\n\t"\
		"vsubpd	%%ymm2,%%ymm0,%%ymm0	\n\t"\
		"vsubpd	%%ymm3,%%ymm1,%%ymm1	\n\t"\
		"vaddpd	%%ymm2,%%ymm2,%%ymm2	\n\t"\
		"vaddpd	%%ymm3,%%ymm3,%%ymm3	\n\t"\
		"vaddpd	%%ymm0,%%ymm2,%%ymm2	\n\t"\
		"vaddpd	%%ymm1,%%ymm3,%%ymm3	\n\t"\
		"vsubpd	%%ymm4,%%ymm0,%%ymm0	\n\t"\
		"vsubpd	%%ymm7,%%ymm2,%%ymm2	\n\t"\
		"vsubpd	%%ymm5,%%ymm1,%%ymm1	\n\t"\
		"vsubpd	%%ymm6,%%ymm3,%%ymm3	\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4	\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7	\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6	\n\t"\
		"vmovaps	%%ymm0,     (%%rbx)	\n\t"\
		"vmovaps	%%ymm2,     (%%rcx)	\n\t"\
		"vmovaps	%%ymm1,0x020(%%rbx)	\n\t"\
		"vmovaps	%%ymm3,0x020(%%rdx)	\n\t"\
		"vaddpd	%%ymm0,%%ymm4,%%ymm4	\n\t"\
		"vaddpd	%%ymm2,%%ymm7,%%ymm7	\n\t"\
		"vaddpd	%%ymm1,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm3,%%ymm6,%%ymm6	\n\t"\
		"vmovaps	%%ymm4,     (%%rax)	\n\t"\
		"vmovaps	%%ymm7,     (%%rdx)	\n\t"\
		"vmovaps	%%ymm5,0x020(%%rax)	\n\t"\
		"vmovaps	%%ymm6,0x020(%%rcx)	\n\t"\
		"\n\t"\
		"/*...Block 4: t06,t16,t26,t36	*/	\n\t"\
		"subq	$0x80,%%rsi	/* r06 */	\n\t"\
		"movslq	%[__p18],%%rdi	\n\t"\
		"subq	%%rax,%%rbx	\n\t"\
		"subq	%%rax,%%rcx	\n\t"\
		"subq	%%rax,%%rdx	\n\t"\
		"movq	%[__add0],%%rax	/* &a[j1] */	\n\t"\
		"shlq	$3,%%rdi	\n\t"\
		"addq	%%rdi,%%rax	/* &a[j1+p18] */	\n\t"\
		"addq	%%rax,%%rbx	\n\t"\
		"addq	%%rax,%%rcx	\n\t"\
		"addq	%%rax,%%rdx	\n\t"\
		"vmovaps	0x400(%%rsi),%%ymm4	\n\t"\
		"vmovaps	0x420(%%rsi),%%ymm5	\n\t"\
		"vmovaps	0x7e0(%%rsi),%%ymm2	/* cc3 */	\n\t"\
		"vmovaps	0x800(%%rsi),%%ymm3	\n\t"\
		"vmovaps	%%ymm4,%%ymm6	\n\t"\
		"vmovaps	%%ymm5,%%ymm7	\n\t"\
		"vmulpd	%%ymm2,%%ymm4,%%ymm4	\n\t"\
		"vmulpd	%%ymm2,%%ymm5,%%ymm5	\n\t"\
		"vmulpd	%%ymm3,%%ymm6,%%ymm6	\n\t"\
		"vmovaps	0x600(%%rsi),%%ymm0	\n\t"\
		"vmulpd	%%ymm3,%%ymm7,%%ymm7	\n\t"\
		"vmovaps	0x620(%%rsi),%%ymm1	\n\t"\
		"vmovaps	0x7a0(%%rsi),%%ymm3	/* cc1 */	\n\t"\
		"vmovaps	0x7c0(%%rsi),%%ymm2	\n\t"\
		"vaddpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
		"vmovaps	%%ymm0,%%ymm6	\n\t"\
		"vsubpd	%%ymm7,%%ymm4,%%ymm4	\n\t"\
		"vmovaps	%%ymm1,%%ymm7	\n\t"\
		"vmulpd	%%ymm2,%%ymm6,%%ymm6	\n\t"\
		"vmulpd	%%ymm2,%%ymm7,%%ymm7	\n\t"\
		"vmulpd	%%ymm3,%%ymm0,%%ymm0	\n\t"\
		"vmulpd	%%ymm3,%%ymm1,%%ymm1	\n\t"\
		"vsubpd	%%ymm0,%%ymm7,%%ymm7	\n\t"\
		"vaddpd	%%ymm1,%%ymm6,%%ymm6	\n\t"\
		"vmovaps	%%ymm4,%%ymm2	\n\t"\
		"vmovaps	%%ymm5,%%ymm3	\n\t"\
		"vsubpd	%%ymm6,%%ymm4,%%ymm4	\n\t"\
		"vsubpd	%%ymm7,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm2,%%ymm6,%%ymm6	\n\t"\
		"vaddpd	%%ymm3,%%ymm7,%%ymm7	\n\t"\
		"vmovaps	0x200(%%rsi),%%ymm1	\n\t"\
		"vmovaps	0x220(%%rsi),%%ymm3	\n\t"\
		"vmovaps	0x760(%%rsi),%%ymm0		/* cc0 */	\n\t"\
		"vmovaps	%%ymm1,%%ymm2	\n\t"\
		"vmulpd	%%ymm0,%%ymm1,%%ymm1	\n\t"\
		"vmulpd	%%ymm3,%%ymm0,%%ymm0	\n\t"\
		"vmulpd	0x780(%%rsi),%%ymm2,%%ymm2		/* ss0 */	\n\t"\
		"vmulpd	0x780(%%rsi),%%ymm3,%%ymm3	\n\t"\
		"vsubpd	%%ymm0,%%ymm2,%%ymm2	\n\t"\
		"vaddpd	%%ymm1,%%ymm3,%%ymm3	\n\t"\
		"vmovaps	     (%%rsi),%%ymm0	\n\t"\
		"vmovaps	0x020(%%rsi),%%ymm1	\n\t"\
		"vsubpd	%%ymm2,%%ymm0,%%ymm0	\n\t"\
		"vsubpd	%%ymm3,%%ymm1,%%ymm1	\n\t"\
		"vaddpd	%%ymm2,%%ymm2,%%ymm2	\n\t"\
		"vaddpd	%%ymm3,%%ymm3,%%ymm3	\n\t"\
		"vaddpd	%%ymm0,%%ymm2,%%ymm2	\n\t"\
		"vaddpd	%%ymm1,%%ymm3,%%ymm3	\n\t"\
		"vsubpd	%%ymm4,%%ymm2,%%ymm2	\n\t"\
		"vsubpd	%%ymm7,%%ymm0,%%ymm0	\n\t"\
		"vsubpd	%%ymm5,%%ymm3,%%ymm3	\n\t"\
		"vsubpd	%%ymm6,%%ymm1,%%ymm1	\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4	\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7	\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6	\n\t"\
		"vmovaps	%%ymm2,     (%%rbx)	\n\t"\
		"vmovaps	%%ymm0,     (%%rcx)	\n\t"\
		"vmovaps	%%ymm3,0x020(%%rbx)	\n\t"\
		"vmovaps	%%ymm1,0x020(%%rdx)	\n\t"\
		"vaddpd	%%ymm2,%%ymm4,%%ymm4	\n\t"\
		"vaddpd	%%ymm0,%%ymm7,%%ymm7	\n\t"\
		"vaddpd	%%ymm3,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm1,%%ymm6,%%ymm6	\n\t"\
		"vmovaps	%%ymm4,     (%%rax)	\n\t"\
		"vmovaps	%%ymm7,     (%%rdx)	\n\t"\
		"vmovaps	%%ymm5,0x020(%%rax)	\n\t"\
		"vmovaps	%%ymm6,0x020(%%rcx)	\n\t"\
		"\n\t"\
		"/*...Block 8: t0E,t1E,t2E,t3E	*/	\n\t"\
		"addq	$0x100,%%rsi	/* r0E */	\n\t"\
		"subq	%%rax,%%rbx	/* p01 << 3 */	\n\t"\
		"subq	%%rax,%%rcx	/* p02 << 3 */	\n\t"\
		"subq	%%rax,%%rdx	/* p03 << 3 */	\n\t"\
		"movslq	%[__p04],%%rdi	\n\t"\
		"shlq	$3,%%rdi	\n\t"\
		"addq	%%rdi,%%rax	/* &a[j1+p1C] */	\n\t"\
		"addq	%%rax,%%rbx	\n\t"\
		"addq	%%rax,%%rcx	\n\t"\
		"addq	%%rax,%%rdx	\n\t"\
		"vmovaps	0x400(%%rsi),%%ymm4	\n\t"\
		"vmovaps	0x420(%%rsi),%%ymm5	\n\t"\
		"vmovaps	0x6a0(%%rsi),%%ymm3	/* cc1 */	\n\t"\
		"vmovaps	0x6c0(%%rsi),%%ymm2	\n\t"\
		"vmovaps	%%ymm4,%%ymm6	\n\t"\
		"vmovaps	%%ymm5,%%ymm7	\n\t"\
		"vmulpd	%%ymm2,%%ymm4,%%ymm4	\n\t"\
		"vmulpd	%%ymm2,%%ymm5,%%ymm5	\n\t"\
		"vmulpd	%%ymm3,%%ymm6,%%ymm6	\n\t"\
		"vmovaps	0x600(%%rsi),%%ymm0	\n\t"\
		"vmulpd	%%ymm3,%%ymm7,%%ymm7	\n\t"\
		"vmovaps	0x620(%%rsi),%%ymm1	\n\t"\
		"vmovaps	0x6e0(%%rsi),%%ymm3	/* cc3 */	\n\t"\
		"vmovaps	0x700(%%rsi),%%ymm2	\n\t"\
		"vaddpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
		"vmovaps	%%ymm0,%%ymm6	\n\t"\
		"vsubpd	%%ymm7,%%ymm4,%%ymm4	\n\t"\
		"vmovaps	%%ymm1,%%ymm7	\n\t"\
		"vmulpd	%%ymm2,%%ymm6,%%ymm6	\n\t"\
		"vmulpd	%%ymm2,%%ymm7,%%ymm7	\n\t"\
		"vmulpd	%%ymm3,%%ymm0,%%ymm0	\n\t"\
		"vmulpd	%%ymm3,%%ymm1,%%ymm1	\n\t"\
		"vaddpd	%%ymm0,%%ymm7,%%ymm7	\n\t"\
		"vsubpd	%%ymm1,%%ymm6,%%ymm6	\n\t"\
		"vmovaps	%%ymm4,%%ymm2	\n\t"\
		"vmovaps	%%ymm5,%%ymm3	\n\t"\
		"vsubpd	%%ymm6,%%ymm4,%%ymm4	\n\t"\
		"vsubpd	%%ymm7,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm2,%%ymm6,%%ymm6	\n\t"\
		"vaddpd	%%ymm3,%%ymm7,%%ymm7	\n\t"\
		"vmovaps	0x200(%%rsi),%%ymm1	\n\t"\
		"vmovaps	0x220(%%rsi),%%ymm3	\n\t"\
		"vmovaps	0x680(%%rsi),%%ymm0	/* ss0 */	\n\t"\
		"vmovaps	%%ymm1,%%ymm2	\n\t"\
		"vmulpd	%%ymm0,%%ymm1,%%ymm1	\n\t"\
		"vmulpd	%%ymm3,%%ymm0,%%ymm0	\n\t"\
		"vmulpd	0x660(%%rsi),%%ymm2,%%ymm2	/* cc0 */	\n\t"\
		"vmulpd	0x660(%%rsi),%%ymm3,%%ymm3	\n\t"\
		"vaddpd	%%ymm0,%%ymm2,%%ymm2	\n\t"\
		"vsubpd	%%ymm1,%%ymm3,%%ymm3	\n\t"\
		"vmovaps	     (%%rsi),%%ymm0	\n\t"\
		"vmovaps	0x020(%%rsi),%%ymm1	\n\t"\
		"vsubpd	%%ymm2,%%ymm0,%%ymm0	\n\t"\
		"vsubpd	%%ymm3,%%ymm1,%%ymm1	\n\t"\
		"vaddpd	%%ymm2,%%ymm2,%%ymm2	\n\t"\
		"vaddpd	%%ymm3,%%ymm3,%%ymm3	\n\t"\
		"vaddpd	%%ymm0,%%ymm2,%%ymm2	\n\t"\
		"vaddpd	%%ymm1,%%ymm3,%%ymm3	\n\t"\
		"vsubpd	%%ymm4,%%ymm0,%%ymm0	\n\t"\
		"vsubpd	%%ymm7,%%ymm2,%%ymm2	\n\t"\
		"vsubpd	%%ymm5,%%ymm1,%%ymm1	\n\t"\
		"vsubpd	%%ymm6,%%ymm3,%%ymm3	\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4	\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7	\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6	\n\t"\
		"vmovaps	%%ymm0,     (%%rbx)	\n\t"\
		"vmovaps	%%ymm2,     (%%rcx)	\n\t"\
		"vmovaps	%%ymm1,0x020(%%rbx)	\n\t"\
		"vmovaps	%%ymm3,0x020(%%rdx)	\n\t"\
		"vaddpd	%%ymm0,%%ymm4,%%ymm4	\n\t"\
		"vaddpd	%%ymm2,%%ymm7,%%ymm7	\n\t"\
		"vaddpd	%%ymm1,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm3,%%ymm6,%%ymm6	\n\t"\
		"vmovaps	%%ymm4,     (%%rax)	\n\t"\
		"vmovaps	%%ymm7,     (%%rdx)	\n\t"\
		"vmovaps	%%ymm5,0x020(%%rax)	\n\t"\
		"vmovaps	%%ymm6,0x020(%%rcx)	\n\t"\
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
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"	/* Clobbered registers */\
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
		"movq	%[__add0],%%rax\n\t"\
		"movslq	%[__p05],%%rbx\n\t"\
		"movslq	%[__p06],%%rcx\n\t"\
		"movslq	%[__p07],%%rdx\n\t"\
		"movslq	%[__p04],%%rdi	/* rdi will store copy of p4 throughout */\n\t"\
		"shlq	$3,%%rbx\n\t"\
		"shlq	$3,%%rcx\n\t"\
		"shlq	$3,%%rdx\n\t"\
		"shlq	$3,%%rdi\n\t"\
		"addq	%%rax,%%rbx\n\t"\
		"addq	%%rax,%%rcx\n\t"\
		"addq	%%rax,%%rdx\n\t"\
		"addq	%%rdi,%%rax\n\t"\
		"/* SSE2_RADIX8_DIT_0TWIDDLE(r00): */\n\t"\
		"movq	%[__r00]	,%%rsi 		\n\t"\
		"vmovaps	     (%%rax),%%ymm0		\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1		\n\t"\
		"vmovaps	%%ymm0		,%%ymm2		\n\t"\
		"vmovaps	%%ymm1		,%%ymm3		\n\t"\
		"vaddpd	     (%%rbx),%%ymm2,%%ymm2		\n\t"\
		"vaddpd	0x020(%%rbx),%%ymm3,%%ymm3		\n\t"\
		"vsubpd	     (%%rbx),%%ymm0,%%ymm0		\n\t"\
		"vsubpd	0x020(%%rbx),%%ymm1,%%ymm1		\n\t"\
		"vmovaps	     (%%rcx),%%ymm4		\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5		\n\t"\
		"vmovaps	%%ymm4		,%%ymm6		\n\t"\
		"vmovaps	%%ymm5		,%%ymm7		\n\t"\
		"vaddpd	     (%%rdx),%%ymm6,%%ymm6		\n\t"\
		"vaddpd	0x020(%%rdx),%%ymm7,%%ymm7		\n\t"\
		"vsubpd	     (%%rdx),%%ymm4,%%ymm4		\n\t"\
		"vsubpd	0x020(%%rdx),%%ymm5,%%ymm5		\n\t"\
		"vmovaps	%%ymm6		,0x180(%%rsi)\n\t"\
		"vmovaps	%%ymm7		,0x1a0(%%rsi)\n\t"\
		"vaddpd	%%ymm2		,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	%%ymm3		,%%ymm7,%%ymm7		\n\t"\
		"vsubpd	0x180(%%rsi),%%ymm2,%%ymm2		\n\t"\
		"vsubpd	0x1a0(%%rsi),%%ymm3,%%ymm3		\n\t"\
		"vmovaps	%%ymm2		,0x180(%%rsi)\n\t"\
		"vmovaps	%%ymm3		,0x1a0(%%rsi)\n\t"\
		"vmovaps	%%ymm4		,%%ymm2		\n\t"\
		"vmovaps	%%ymm5		,%%ymm3		\n\t"\
		"vaddpd	%%ymm0		,%%ymm5,%%ymm5		\n\t"\
		"vsubpd	%%ymm3		,%%ymm0,%%ymm0		\n\t"\
		"vaddpd	%%ymm1		,%%ymm4,%%ymm4		\n\t"\
		"vsubpd	%%ymm2		,%%ymm1,%%ymm1		\n\t"\
		"movq	%[__isrt2]	,%%rsi		\n\t"\
		"vmovaps	%%ymm5		,%%ymm2		\n\t"\
		"vmovaps	%%ymm1		,%%ymm3		\n\t"\
		"vaddpd	%%ymm1		,%%ymm5,%%ymm5		\n\t"\
		"vmovaps	(%%rsi)		,%%ymm1		\n\t"\
		"vsubpd	%%ymm3		,%%ymm2,%%ymm2		\n\t"\
		"movq	%[__r00]	,%%rsi		\n\t"\
		"vmulpd	%%ymm1		,%%ymm5,%%ymm5		\n\t"\
		"vmulpd	%%ymm1		,%%ymm2,%%ymm2		\n\t"\
		"vmovaps	%%ymm0		,%%ymm3		\n\t"\
		"vmovaps	%%ymm5		,0x140(%%rsi)\n\t"\
		"vmovaps	%%ymm4		,%%ymm5		\n\t"\
		"vaddpd	%%ymm4		,%%ymm0,%%ymm0		\n\t"\
		"vmovaps	%%ymm2		,0x160(%%rsi)\n\t"\
		"vsubpd	%%ymm5		,%%ymm3,%%ymm3		\n\t"\
		"vmulpd	%%ymm1		,%%ymm0,%%ymm0		\n\t"\
		"vmulpd	%%ymm1		,%%ymm3,%%ymm3		\n\t"\
		"vmovaps	%%ymm0		,0x1c0(%%rsi)\n\t"\
		"vmovaps	%%ymm3		,0x1e0(%%rsi)\n\t"\
		"subq	%%rdi		,%%rax		\n\t"\
		"subq	%%rdi		,%%rbx		\n\t"\
		"subq	%%rdi		,%%rcx		\n\t"\
		"subq	%%rdi		,%%rdx		\n\t"\
		"vmovaps	     (%%rax),%%ymm0		\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1		\n\t"\
		"vmovaps	%%ymm0		,%%ymm2		\n\t"\
		"vmovaps	%%ymm1		,%%ymm3		\n\t"\
		"vaddpd	     (%%rbx),%%ymm2,%%ymm2		\n\t"\
		"vaddpd	0x020(%%rbx),%%ymm3,%%ymm3		\n\t"\
		"vsubpd	     (%%rbx),%%ymm0,%%ymm0		\n\t"\
		"vsubpd	0x020(%%rbx),%%ymm1,%%ymm1		\n\t"\
		"vmovaps	%%ymm6		,     (%%rsi)\n\t"\
		"vmovaps	%%ymm7		,0x020(%%rsi)\n\t"\
		"vaddpd	%%ymm6		,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	%%ymm7		,%%ymm7,%%ymm7		\n\t"\
		"vmovaps	%%ymm6		,0x100(%%rsi)\n\t"\
		"vmovaps	%%ymm7		,0x120(%%rsi)\n\t"\
		"vmovaps	     (%%rcx),%%ymm4		\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5		\n\t"\
		"vmovaps	%%ymm4		,%%ymm6		\n\t"\
		"vmovaps	%%ymm5		,%%ymm7		\n\t"\
		"vaddpd	     (%%rdx),%%ymm6,%%ymm6		\n\t"\
		"vaddpd	0x020(%%rdx),%%ymm7,%%ymm7		\n\t"\
		"vsubpd	     (%%rdx),%%ymm4,%%ymm4		\n\t"\
		"vsubpd	0x020(%%rdx),%%ymm5,%%ymm5		\n\t"\
		"vmovaps	%%ymm6		,0x080(%%rsi)\n\t"\
		"vmovaps	%%ymm7		,0x0a0(%%rsi)\n\t"\
		"vaddpd	%%ymm2		,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	%%ymm3		,%%ymm7,%%ymm7		\n\t"\
		"vsubpd	0x080(%%rsi),%%ymm2,%%ymm2		\n\t"\
		"vsubpd	0x0a0(%%rsi),%%ymm3,%%ymm3		\n\t"\
		"vaddpd	     (%%rsi),%%ymm6,%%ymm6		\n\t"\
		"vaddpd	0x020(%%rsi),%%ymm7,%%ymm7		\n\t"\
		"vmovaps	%%ymm6		,     (%%rsi)\n\t"\
		"vmovaps	%%ymm7		,0x020(%%rsi)\n\t"\
		"vsubpd	0x100(%%rsi),%%ymm6,%%ymm6		\n\t"\
		"vsubpd	0x120(%%rsi),%%ymm7,%%ymm7		\n\t"\
		"vmovaps	%%ymm6		,0x100(%%rsi)\n\t"\
		"vmovaps	%%ymm7		,0x120(%%rsi)\n\t"\
		"vmovaps	%%ymm4		,%%ymm6		\n\t"\
		"vmovaps	%%ymm5		,%%ymm7		\n\t"\
		"vaddpd	%%ymm0		,%%ymm5,%%ymm5		\n\t"\
		"vsubpd	%%ymm7		,%%ymm0,%%ymm0		\n\t"\
		"vaddpd	%%ymm1		,%%ymm4,%%ymm4		\n\t"\
		"vsubpd	%%ymm6		,%%ymm1,%%ymm1		\n\t"\
		"vmovaps	%%ymm5		,%%ymm6		\n\t"\
		"vmovaps	%%ymm1		,%%ymm7		\n\t"\
		"vaddpd	0x140(%%rsi),%%ymm5,%%ymm5		\n\t"\
		"vsubpd	0x160(%%rsi),%%ymm1,%%ymm1		\n\t"\
		"vsubpd	0x140(%%rsi),%%ymm6,%%ymm6		\n\t"\
		"vaddpd	0x160(%%rsi),%%ymm7,%%ymm7		\n\t"\
		"vmovaps	%%ymm5		,0x040(%%rsi)\n\t"\
		"vmovaps	%%ymm1		,0x060(%%rsi)\n\t"\
		"vmovaps	%%ymm6		,0x140(%%rsi)\n\t"\
		"vmovaps	%%ymm7		,0x160(%%rsi)\n\t"\
		"vmovaps	%%ymm2		,%%ymm6		\n\t"\
		"vmovaps	%%ymm3		,%%ymm7		\n\t"\
		"vaddpd	0x1a0(%%rsi),%%ymm2,%%ymm2		\n\t"\
		"vsubpd	0x180(%%rsi),%%ymm3,%%ymm3		\n\t"\
		"vsubpd	0x1a0(%%rsi),%%ymm6,%%ymm6		\n\t"\
		"vaddpd	0x180(%%rsi),%%ymm7,%%ymm7		\n\t"\
		"vmovaps	%%ymm2		,0x080(%%rsi)\n\t"\
		"vmovaps	%%ymm3		,0x0a0(%%rsi)\n\t"\
		"vmovaps	%%ymm6		,0x180(%%rsi)\n\t"\
		"vmovaps	%%ymm7		,0x1a0(%%rsi)\n\t"\
		"vmovaps	%%ymm0		,%%ymm6		\n\t"\
		"vmovaps	%%ymm4		,%%ymm7		\n\t"\
		"vsubpd	0x1e0(%%rsi),%%ymm0,%%ymm0		\n\t"\
		"vsubpd	0x1c0(%%rsi),%%ymm4,%%ymm4		\n\t"\
		"vaddpd	0x1e0(%%rsi),%%ymm6,%%ymm6		\n\t"\
		"vaddpd	0x1c0(%%rsi),%%ymm7,%%ymm7		\n\t"\
		"vmovaps	%%ymm0		,0x0c0(%%rsi)\n\t"\
		"vmovaps	%%ymm4		,0x0e0(%%rsi)\n\t"\
		"vmovaps	%%ymm6		,0x1c0(%%rsi)\n\t"\
		"vmovaps	%%ymm7		,0x1e0(%%rsi)\n\t"\
		"\n\t"\
		"movq	%[__add0],%%rax\n\t"\
		"movslq	%[__p08],%%rsi\n\t"\
		"shlq	$3,%%rsi\n\t"\
		"addq	%%rsi,%%rax	/* add1 = add0+p8 */\n\t"\
		"movslq	%[__p05],%%rbx\n\t"\
		"movslq	%[__p06],%%rcx\n\t"\
		"movslq	%[__p07],%%rdx\n\t"\
		"movslq	%[__p04],%%rdi\n\t"\
		"shlq	$3,%%rbx\n\t"\
		"shlq	$3,%%rcx\n\t"\
		"shlq	$3,%%rdx\n\t"\
		"shlq	$3,%%rdi\n\t"\
		"addq	%%rax,%%rbx\n\t"\
		"addq	%%rax,%%rcx\n\t"\
		"addq	%%rax,%%rdx\n\t"\
		"addq	%%rdi,%%rax\n\t"\
		"/* SSE2_RADIX8_DIT_0TWIDDLE(r10): */\n\t"\
		"movq	%[__r10]	,%%rsi 		\n\t"\
		"vmovaps	     (%%rax),%%ymm0		\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1		\n\t"\
		"vmovaps	%%ymm0		,%%ymm2		\n\t"\
		"vmovaps	%%ymm1		,%%ymm3		\n\t"\
		"vaddpd	     (%%rbx),%%ymm2,%%ymm2		\n\t"\
		"vaddpd	0x020(%%rbx),%%ymm3,%%ymm3		\n\t"\
		"vsubpd	     (%%rbx),%%ymm0,%%ymm0		\n\t"\
		"vsubpd	0x020(%%rbx),%%ymm1,%%ymm1		\n\t"\
		"vmovaps	     (%%rcx),%%ymm4		\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5		\n\t"\
		"vmovaps	%%ymm4		,%%ymm6		\n\t"\
		"vmovaps	%%ymm5		,%%ymm7		\n\t"\
		"vaddpd	     (%%rdx),%%ymm6,%%ymm6		\n\t"\
		"vaddpd	0x020(%%rdx),%%ymm7,%%ymm7		\n\t"\
		"vsubpd	     (%%rdx),%%ymm4,%%ymm4		\n\t"\
		"vsubpd	0x020(%%rdx),%%ymm5,%%ymm5		\n\t"\
		"vmovaps	%%ymm6		,0x180(%%rsi)\n\t"\
		"vmovaps	%%ymm7		,0x1a0(%%rsi)\n\t"\
		"vaddpd	%%ymm2		,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	%%ymm3		,%%ymm7,%%ymm7		\n\t"\
		"vsubpd	0x180(%%rsi),%%ymm2,%%ymm2		\n\t"\
		"vsubpd	0x1a0(%%rsi),%%ymm3,%%ymm3		\n\t"\
		"vmovaps	%%ymm2		,0x180(%%rsi)\n\t"\
		"vmovaps	%%ymm3		,0x1a0(%%rsi)\n\t"\
		"vmovaps	%%ymm4		,%%ymm2		\n\t"\
		"vmovaps	%%ymm5		,%%ymm3		\n\t"\
		"vaddpd	%%ymm0		,%%ymm5,%%ymm5		\n\t"\
		"vsubpd	%%ymm3		,%%ymm0,%%ymm0		\n\t"\
		"vaddpd	%%ymm1		,%%ymm4,%%ymm4		\n\t"\
		"vsubpd	%%ymm2		,%%ymm1,%%ymm1		\n\t"\
		"movq	%[__isrt2]	,%%rsi		\n\t"\
		"vmovaps	%%ymm5		,%%ymm2		\n\t"\
		"vmovaps	%%ymm1		,%%ymm3		\n\t"\
		"vaddpd	%%ymm1		,%%ymm5,%%ymm5		\n\t"\
		"vmovaps	(%%rsi)		,%%ymm1		\n\t"\
		"vsubpd	%%ymm3		,%%ymm2,%%ymm2		\n\t"\
		"movq	%[__r10]	,%%rsi		\n\t"\
		"vmulpd	%%ymm1		,%%ymm5,%%ymm5		\n\t"\
		"vmulpd	%%ymm1		,%%ymm2,%%ymm2		\n\t"\
		"vmovaps	%%ymm0		,%%ymm3		\n\t"\
		"vmovaps	%%ymm5		,0x140(%%rsi)\n\t"\
		"vmovaps	%%ymm4		,%%ymm5		\n\t"\
		"vaddpd	%%ymm4		,%%ymm0,%%ymm0		\n\t"\
		"vmovaps	%%ymm2		,0x160(%%rsi)\n\t"\
		"vsubpd	%%ymm5		,%%ymm3,%%ymm3		\n\t"\
		"vmulpd	%%ymm1		,%%ymm0,%%ymm0		\n\t"\
		"vmulpd	%%ymm1		,%%ymm3,%%ymm3		\n\t"\
		"vmovaps	%%ymm0		,0x1c0(%%rsi)\n\t"\
		"vmovaps	%%ymm3		,0x1e0(%%rsi)\n\t"\
		"subq	%%rdi		,%%rax		\n\t"\
		"subq	%%rdi		,%%rbx		\n\t"\
		"subq	%%rdi		,%%rcx		\n\t"\
		"subq	%%rdi		,%%rdx		\n\t"\
		"vmovaps	     (%%rax),%%ymm0		\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1		\n\t"\
		"vmovaps	%%ymm0		,%%ymm2		\n\t"\
		"vmovaps	%%ymm1		,%%ymm3		\n\t"\
		"vaddpd	     (%%rbx),%%ymm2,%%ymm2		\n\t"\
		"vaddpd	0x020(%%rbx),%%ymm3,%%ymm3		\n\t"\
		"vsubpd	     (%%rbx),%%ymm0,%%ymm0		\n\t"\
		"vsubpd	0x020(%%rbx),%%ymm1,%%ymm1		\n\t"\
		"vmovaps	%%ymm6		,     (%%rsi)\n\t"\
		"vmovaps	%%ymm7		,0x020(%%rsi)\n\t"\
		"vaddpd	%%ymm6		,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	%%ymm7		,%%ymm7,%%ymm7		\n\t"\
		"vmovaps	%%ymm6		,0x100(%%rsi)\n\t"\
		"vmovaps	%%ymm7		,0x120(%%rsi)\n\t"\
		"vmovaps	     (%%rcx),%%ymm4		\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5		\n\t"\
		"vmovaps	%%ymm4		,%%ymm6		\n\t"\
		"vmovaps	%%ymm5		,%%ymm7		\n\t"\
		"vaddpd	     (%%rdx),%%ymm6,%%ymm6		\n\t"\
		"vaddpd	0x020(%%rdx),%%ymm7,%%ymm7		\n\t"\
		"vsubpd	     (%%rdx),%%ymm4,%%ymm4		\n\t"\
		"vsubpd	0x020(%%rdx),%%ymm5,%%ymm5		\n\t"\
		"vmovaps	%%ymm6		,0x080(%%rsi)\n\t"\
		"vmovaps	%%ymm7		,0x0a0(%%rsi)\n\t"\
		"vaddpd	%%ymm2		,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	%%ymm3		,%%ymm7,%%ymm7		\n\t"\
		"vsubpd	0x080(%%rsi),%%ymm2,%%ymm2		\n\t"\
		"vsubpd	0x0a0(%%rsi),%%ymm3,%%ymm3		\n\t"\
		"vaddpd	     (%%rsi),%%ymm6,%%ymm6		\n\t"\
		"vaddpd	0x020(%%rsi),%%ymm7,%%ymm7		\n\t"\
		"vmovaps	%%ymm6		,     (%%rsi)\n\t"\
		"vmovaps	%%ymm7		,0x020(%%rsi)\n\t"\
		"vsubpd	0x100(%%rsi),%%ymm6,%%ymm6		\n\t"\
		"vsubpd	0x120(%%rsi),%%ymm7,%%ymm7		\n\t"\
		"vmovaps	%%ymm6		,0x100(%%rsi)\n\t"\
		"vmovaps	%%ymm7		,0x120(%%rsi)\n\t"\
		"vmovaps	%%ymm4		,%%ymm6		\n\t"\
		"vmovaps	%%ymm5		,%%ymm7		\n\t"\
		"vaddpd	%%ymm0		,%%ymm5,%%ymm5		\n\t"\
		"vsubpd	%%ymm7		,%%ymm0,%%ymm0		\n\t"\
		"vaddpd	%%ymm1		,%%ymm4,%%ymm4		\n\t"\
		"vsubpd	%%ymm6		,%%ymm1,%%ymm1		\n\t"\
		"vmovaps	%%ymm5		,%%ymm6		\n\t"\
		"vmovaps	%%ymm1		,%%ymm7		\n\t"\
		"vaddpd	0x140(%%rsi),%%ymm5,%%ymm5		\n\t"\
		"vsubpd	0x160(%%rsi),%%ymm1,%%ymm1		\n\t"\
		"vsubpd	0x140(%%rsi),%%ymm6,%%ymm6		\n\t"\
		"vaddpd	0x160(%%rsi),%%ymm7,%%ymm7		\n\t"\
		"vmovaps	%%ymm5		,0x040(%%rsi)\n\t"\
		"vmovaps	%%ymm1		,0x060(%%rsi)\n\t"\
		"vmovaps	%%ymm6		,0x140(%%rsi)\n\t"\
		"vmovaps	%%ymm7		,0x160(%%rsi)\n\t"\
		"vmovaps	%%ymm2		,%%ymm6		\n\t"\
		"vmovaps	%%ymm3		,%%ymm7		\n\t"\
		"vaddpd	0x1a0(%%rsi),%%ymm2,%%ymm2		\n\t"\
		"vsubpd	0x180(%%rsi),%%ymm3,%%ymm3		\n\t"\
		"vsubpd	0x1a0(%%rsi),%%ymm6,%%ymm6		\n\t"\
		"vaddpd	0x180(%%rsi),%%ymm7,%%ymm7		\n\t"\
		"vmovaps	%%ymm2		,0x080(%%rsi)\n\t"\
		"vmovaps	%%ymm3		,0x0a0(%%rsi)\n\t"\
		"vmovaps	%%ymm6		,0x180(%%rsi)\n\t"\
		"vmovaps	%%ymm7		,0x1a0(%%rsi)\n\t"\
		"vmovaps	%%ymm0		,%%ymm6		\n\t"\
		"vmovaps	%%ymm4		,%%ymm7		\n\t"\
		"vsubpd	0x1e0(%%rsi),%%ymm0,%%ymm0		\n\t"\
		"vsubpd	0x1c0(%%rsi),%%ymm4,%%ymm4		\n\t"\
		"vaddpd	0x1e0(%%rsi),%%ymm6,%%ymm6		\n\t"\
		"vaddpd	0x1c0(%%rsi),%%ymm7,%%ymm7		\n\t"\
		"vmovaps	%%ymm0		,0x0c0(%%rsi)\n\t"\
		"vmovaps	%%ymm4		,0x0e0(%%rsi)\n\t"\
		"vmovaps	%%ymm6		,0x1c0(%%rsi)\n\t"\
		"vmovaps	%%ymm7		,0x1e0(%%rsi)\n\t"\
		"\n\t"\
		"movq	%[__add0],%%rax\n\t"\
		"movslq	%[__p10],%%rsi\n\t"\
		"shlq	$3,%%rsi\n\t"\
		"addq	%%rsi,%%rax	/* add2 = add0+p10*/\n\t"\
		"movslq	%[__p05],%%rbx\n\t"\
		"movslq	%[__p06],%%rcx\n\t"\
		"movslq	%[__p07],%%rdx\n\t"\
		"movslq	%[__p04],%%rdi\n\t"\
		"shlq	$3,%%rbx\n\t"\
		"shlq	$3,%%rcx\n\t"\
		"shlq	$3,%%rdx\n\t"\
		"shlq	$3,%%rdi\n\t"\
		"addq	%%rax,%%rbx\n\t"\
		"addq	%%rax,%%rcx\n\t"\
		"addq	%%rax,%%rdx\n\t"\
		"addq	%%rdi,%%rax\n\t"\
		"/* SSE2_RADIX8_DIT_0TWIDDLE(r20): */\n\t"\
		"movq	%[__r20]	,%%rsi 		\n\t"\
		"vmovaps	     (%%rax),%%ymm0		\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1		\n\t"\
		"vmovaps	%%ymm0		,%%ymm2		\n\t"\
		"vmovaps	%%ymm1		,%%ymm3		\n\t"\
		"vaddpd	     (%%rbx),%%ymm2,%%ymm2		\n\t"\
		"vaddpd	0x020(%%rbx),%%ymm3,%%ymm3		\n\t"\
		"vsubpd	     (%%rbx),%%ymm0,%%ymm0		\n\t"\
		"vsubpd	0x020(%%rbx),%%ymm1,%%ymm1		\n\t"\
		"vmovaps	     (%%rcx),%%ymm4		\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5		\n\t"\
		"vmovaps	%%ymm4		,%%ymm6		\n\t"\
		"vmovaps	%%ymm5		,%%ymm7		\n\t"\
		"vaddpd	     (%%rdx),%%ymm6,%%ymm6		\n\t"\
		"vaddpd	0x020(%%rdx),%%ymm7,%%ymm7		\n\t"\
		"vsubpd	     (%%rdx),%%ymm4,%%ymm4		\n\t"\
		"vsubpd	0x020(%%rdx),%%ymm5,%%ymm5		\n\t"\
		"vmovaps	%%ymm6		,0x180(%%rsi)\n\t"\
		"vmovaps	%%ymm7		,0x1a0(%%rsi)\n\t"\
		"vaddpd	%%ymm2		,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	%%ymm3		,%%ymm7,%%ymm7		\n\t"\
		"vsubpd	0x180(%%rsi),%%ymm2,%%ymm2		\n\t"\
		"vsubpd	0x1a0(%%rsi),%%ymm3,%%ymm3		\n\t"\
		"vmovaps	%%ymm2		,0x180(%%rsi)\n\t"\
		"vmovaps	%%ymm3		,0x1a0(%%rsi)\n\t"\
		"vmovaps	%%ymm4		,%%ymm2		\n\t"\
		"vmovaps	%%ymm5		,%%ymm3		\n\t"\
		"vaddpd	%%ymm0		,%%ymm5,%%ymm5		\n\t"\
		"vsubpd	%%ymm3		,%%ymm0,%%ymm0		\n\t"\
		"vaddpd	%%ymm1		,%%ymm4,%%ymm4		\n\t"\
		"vsubpd	%%ymm2		,%%ymm1,%%ymm1		\n\t"\
		"movq	%[__isrt2]	,%%rsi		\n\t"\
		"vmovaps	%%ymm5		,%%ymm2		\n\t"\
		"vmovaps	%%ymm1		,%%ymm3		\n\t"\
		"vaddpd	%%ymm1		,%%ymm5,%%ymm5		\n\t"\
		"vmovaps	(%%rsi)		,%%ymm1		\n\t"\
		"vsubpd	%%ymm3		,%%ymm2,%%ymm2		\n\t"\
		"movq	%[__r20]	,%%rsi		\n\t"\
		"vmulpd	%%ymm1		,%%ymm5,%%ymm5		\n\t"\
		"vmulpd	%%ymm1		,%%ymm2,%%ymm2		\n\t"\
		"vmovaps	%%ymm0		,%%ymm3		\n\t"\
		"vmovaps	%%ymm5		,0x140(%%rsi)\n\t"\
		"vmovaps	%%ymm4		,%%ymm5		\n\t"\
		"vaddpd	%%ymm4		,%%ymm0,%%ymm0		\n\t"\
		"vmovaps	%%ymm2		,0x160(%%rsi)\n\t"\
		"vsubpd	%%ymm5		,%%ymm3,%%ymm3		\n\t"\
		"vmulpd	%%ymm1		,%%ymm0,%%ymm0		\n\t"\
		"vmulpd	%%ymm1		,%%ymm3,%%ymm3		\n\t"\
		"vmovaps	%%ymm0		,0x1c0(%%rsi)\n\t"\
		"vmovaps	%%ymm3		,0x1e0(%%rsi)\n\t"\
		"subq	%%rdi		,%%rax		\n\t"\
		"subq	%%rdi		,%%rbx		\n\t"\
		"subq	%%rdi		,%%rcx		\n\t"\
		"subq	%%rdi		,%%rdx		\n\t"\
		"vmovaps	     (%%rax),%%ymm0		\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1		\n\t"\
		"vmovaps	%%ymm0		,%%ymm2		\n\t"\
		"vmovaps	%%ymm1		,%%ymm3		\n\t"\
		"vaddpd	     (%%rbx),%%ymm2,%%ymm2		\n\t"\
		"vaddpd	0x020(%%rbx),%%ymm3,%%ymm3		\n\t"\
		"vsubpd	     (%%rbx),%%ymm0,%%ymm0		\n\t"\
		"vsubpd	0x020(%%rbx),%%ymm1,%%ymm1		\n\t"\
		"vmovaps	%%ymm6		,     (%%rsi)\n\t"\
		"vmovaps	%%ymm7		,0x020(%%rsi)\n\t"\
		"vaddpd	%%ymm6		,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	%%ymm7		,%%ymm7,%%ymm7		\n\t"\
		"vmovaps	%%ymm6		,0x100(%%rsi)\n\t"\
		"vmovaps	%%ymm7		,0x120(%%rsi)\n\t"\
		"vmovaps	     (%%rcx),%%ymm4		\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5		\n\t"\
		"vmovaps	%%ymm4		,%%ymm6		\n\t"\
		"vmovaps	%%ymm5		,%%ymm7		\n\t"\
		"vaddpd	     (%%rdx),%%ymm6,%%ymm6		\n\t"\
		"vaddpd	0x020(%%rdx),%%ymm7,%%ymm7		\n\t"\
		"vsubpd	     (%%rdx),%%ymm4,%%ymm4		\n\t"\
		"vsubpd	0x020(%%rdx),%%ymm5,%%ymm5		\n\t"\
		"vmovaps	%%ymm6		,0x080(%%rsi)\n\t"\
		"vmovaps	%%ymm7		,0x0a0(%%rsi)\n\t"\
		"vaddpd	%%ymm2		,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	%%ymm3		,%%ymm7,%%ymm7		\n\t"\
		"vsubpd	0x080(%%rsi),%%ymm2,%%ymm2		\n\t"\
		"vsubpd	0x0a0(%%rsi),%%ymm3,%%ymm3		\n\t"\
		"vaddpd	     (%%rsi),%%ymm6,%%ymm6		\n\t"\
		"vaddpd	0x020(%%rsi),%%ymm7,%%ymm7		\n\t"\
		"vmovaps	%%ymm6		,     (%%rsi)\n\t"\
		"vmovaps	%%ymm7		,0x020(%%rsi)\n\t"\
		"vsubpd	0x100(%%rsi),%%ymm6,%%ymm6		\n\t"\
		"vsubpd	0x120(%%rsi),%%ymm7,%%ymm7		\n\t"\
		"vmovaps	%%ymm6		,0x100(%%rsi)\n\t"\
		"vmovaps	%%ymm7		,0x120(%%rsi)\n\t"\
		"vmovaps	%%ymm4		,%%ymm6		\n\t"\
		"vmovaps	%%ymm5		,%%ymm7		\n\t"\
		"vaddpd	%%ymm0		,%%ymm5,%%ymm5		\n\t"\
		"vsubpd	%%ymm7		,%%ymm0,%%ymm0		\n\t"\
		"vaddpd	%%ymm1		,%%ymm4,%%ymm4		\n\t"\
		"vsubpd	%%ymm6		,%%ymm1,%%ymm1		\n\t"\
		"vmovaps	%%ymm5		,%%ymm6		\n\t"\
		"vmovaps	%%ymm1		,%%ymm7		\n\t"\
		"vaddpd	0x140(%%rsi),%%ymm5,%%ymm5		\n\t"\
		"vsubpd	0x160(%%rsi),%%ymm1,%%ymm1		\n\t"\
		"vsubpd	0x140(%%rsi),%%ymm6,%%ymm6		\n\t"\
		"vaddpd	0x160(%%rsi),%%ymm7,%%ymm7		\n\t"\
		"vmovaps	%%ymm5		,0x040(%%rsi)\n\t"\
		"vmovaps	%%ymm1		,0x060(%%rsi)\n\t"\
		"vmovaps	%%ymm6		,0x140(%%rsi)\n\t"\
		"vmovaps	%%ymm7		,0x160(%%rsi)\n\t"\
		"vmovaps	%%ymm2		,%%ymm6		\n\t"\
		"vmovaps	%%ymm3		,%%ymm7		\n\t"\
		"vaddpd	0x1a0(%%rsi),%%ymm2,%%ymm2		\n\t"\
		"vsubpd	0x180(%%rsi),%%ymm3,%%ymm3		\n\t"\
		"vsubpd	0x1a0(%%rsi),%%ymm6,%%ymm6		\n\t"\
		"vaddpd	0x180(%%rsi),%%ymm7,%%ymm7		\n\t"\
		"vmovaps	%%ymm2		,0x080(%%rsi)\n\t"\
		"vmovaps	%%ymm3		,0x0a0(%%rsi)\n\t"\
		"vmovaps	%%ymm6		,0x180(%%rsi)\n\t"\
		"vmovaps	%%ymm7		,0x1a0(%%rsi)\n\t"\
		"vmovaps	%%ymm0		,%%ymm6		\n\t"\
		"vmovaps	%%ymm4		,%%ymm7		\n\t"\
		"vsubpd	0x1e0(%%rsi),%%ymm0,%%ymm0		\n\t"\
		"vsubpd	0x1c0(%%rsi),%%ymm4,%%ymm4		\n\t"\
		"vaddpd	0x1e0(%%rsi),%%ymm6,%%ymm6		\n\t"\
		"vaddpd	0x1c0(%%rsi),%%ymm7,%%ymm7		\n\t"\
		"vmovaps	%%ymm0		,0x0c0(%%rsi)\n\t"\
		"vmovaps	%%ymm4		,0x0e0(%%rsi)\n\t"\
		"vmovaps	%%ymm6		,0x1c0(%%rsi)\n\t"\
		"vmovaps	%%ymm7		,0x1e0(%%rsi)\n\t"\
		"\n\t"\
		"movq	%[__add0],%%rax\n\t"\
		"movslq	%[__p18],%%rsi\n\t"\
		"shlq	$3,%%rsi\n\t"\
		"addq	%%rsi,%%rax	/* add3 = add0+p18*/\n\t"\
		"movslq	%[__p05],%%rbx\n\t"\
		"movslq	%[__p06],%%rcx\n\t"\
		"movslq	%[__p07],%%rdx\n\t"\
		"movslq	%[__p04],%%rdi\n\t"\
		"shlq	$3,%%rbx\n\t"\
		"shlq	$3,%%rcx\n\t"\
		"shlq	$3,%%rdx\n\t"\
		"shlq	$3,%%rdi\n\t"\
		"addq	%%rax,%%rbx\n\t"\
		"addq	%%rax,%%rcx\n\t"\
		"addq	%%rax,%%rdx\n\t"\
		"addq	%%rdi,%%rax\n\t"\
		"/* SSE2_RADIX8_DIT_0TWIDDLE(r30): */\n\t"\
		"movq	%[__r30]	,%%rsi 		\n\t"\
		"vmovaps	     (%%rax),%%ymm0		\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1		\n\t"\
		"vmovaps	%%ymm0		,%%ymm2		\n\t"\
		"vmovaps	%%ymm1		,%%ymm3		\n\t"\
		"vaddpd	     (%%rbx),%%ymm2,%%ymm2		\n\t"\
		"vaddpd	0x020(%%rbx),%%ymm3,%%ymm3		\n\t"\
		"vsubpd	     (%%rbx),%%ymm0,%%ymm0		\n\t"\
		"vsubpd	0x020(%%rbx),%%ymm1,%%ymm1		\n\t"\
		"vmovaps	     (%%rcx),%%ymm4		\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5		\n\t"\
		"vmovaps	%%ymm4		,%%ymm6		\n\t"\
		"vmovaps	%%ymm5		,%%ymm7		\n\t"\
		"vaddpd	     (%%rdx),%%ymm6,%%ymm6		\n\t"\
		"vaddpd	0x020(%%rdx),%%ymm7,%%ymm7		\n\t"\
		"vsubpd	     (%%rdx),%%ymm4,%%ymm4		\n\t"\
		"vsubpd	0x020(%%rdx),%%ymm5,%%ymm5		\n\t"\
		"vmovaps	%%ymm6		,0x180(%%rsi)\n\t"\
		"vmovaps	%%ymm7		,0x1a0(%%rsi)\n\t"\
		"vaddpd	%%ymm2		,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	%%ymm3		,%%ymm7,%%ymm7		\n\t"\
		"vsubpd	0x180(%%rsi),%%ymm2,%%ymm2		\n\t"\
		"vsubpd	0x1a0(%%rsi),%%ymm3,%%ymm3		\n\t"\
		"vmovaps	%%ymm2		,0x180(%%rsi)\n\t"\
		"vmovaps	%%ymm3		,0x1a0(%%rsi)\n\t"\
		"vmovaps	%%ymm4		,%%ymm2		\n\t"\
		"vmovaps	%%ymm5		,%%ymm3		\n\t"\
		"vaddpd	%%ymm0		,%%ymm5,%%ymm5		\n\t"\
		"vsubpd	%%ymm3		,%%ymm0,%%ymm0		\n\t"\
		"vaddpd	%%ymm1		,%%ymm4,%%ymm4		\n\t"\
		"vsubpd	%%ymm2		,%%ymm1,%%ymm1		\n\t"\
		"movq	%[__isrt2]	,%%rsi		\n\t"\
		"vmovaps	%%ymm5		,%%ymm2		\n\t"\
		"vmovaps	%%ymm1		,%%ymm3		\n\t"\
		"vaddpd	%%ymm1		,%%ymm5,%%ymm5		\n\t"\
		"vmovaps	(%%rsi)		,%%ymm1		\n\t"\
		"vsubpd	%%ymm3		,%%ymm2,%%ymm2		\n\t"\
		"movq	%[__r30]	,%%rsi		\n\t"\
		"vmulpd	%%ymm1		,%%ymm5,%%ymm5		\n\t"\
		"vmulpd	%%ymm1		,%%ymm2,%%ymm2		\n\t"\
		"vmovaps	%%ymm0		,%%ymm3		\n\t"\
		"vmovaps	%%ymm5		,0x140(%%rsi)\n\t"\
		"vmovaps	%%ymm4		,%%ymm5		\n\t"\
		"vaddpd	%%ymm4		,%%ymm0,%%ymm0		\n\t"\
		"vmovaps	%%ymm2		,0x160(%%rsi)\n\t"\
		"vsubpd	%%ymm5		,%%ymm3,%%ymm3		\n\t"\
		"vmulpd	%%ymm1		,%%ymm0,%%ymm0		\n\t"\
		"vmulpd	%%ymm1		,%%ymm3,%%ymm3		\n\t"\
		"vmovaps	%%ymm0		,0x1c0(%%rsi)\n\t"\
		"vmovaps	%%ymm3		,0x1e0(%%rsi)\n\t"\
		"subq	%%rdi		,%%rax		\n\t"\
		"subq	%%rdi		,%%rbx		\n\t"\
		"subq	%%rdi		,%%rcx		\n\t"\
		"subq	%%rdi		,%%rdx		\n\t"\
		"vmovaps	     (%%rax),%%ymm0		\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1		\n\t"\
		"vmovaps	%%ymm0		,%%ymm2		\n\t"\
		"vmovaps	%%ymm1		,%%ymm3		\n\t"\
		"vaddpd	     (%%rbx),%%ymm2,%%ymm2		\n\t"\
		"vaddpd	0x020(%%rbx),%%ymm3,%%ymm3		\n\t"\
		"vsubpd	     (%%rbx),%%ymm0,%%ymm0		\n\t"\
		"vsubpd	0x020(%%rbx),%%ymm1,%%ymm1		\n\t"\
		"vmovaps	%%ymm6		,     (%%rsi)\n\t"\
		"vmovaps	%%ymm7		,0x020(%%rsi)\n\t"\
		"vaddpd	%%ymm6		,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	%%ymm7		,%%ymm7,%%ymm7		\n\t"\
		"vmovaps	%%ymm6		,0x100(%%rsi)\n\t"\
		"vmovaps	%%ymm7		,0x120(%%rsi)\n\t"\
		"vmovaps	     (%%rcx),%%ymm4		\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5		\n\t"\
		"vmovaps	%%ymm4		,%%ymm6		\n\t"\
		"vmovaps	%%ymm5		,%%ymm7		\n\t"\
		"vaddpd	     (%%rdx),%%ymm6,%%ymm6		\n\t"\
		"vaddpd	0x020(%%rdx),%%ymm7,%%ymm7		\n\t"\
		"vsubpd	     (%%rdx),%%ymm4,%%ymm4		\n\t"\
		"vsubpd	0x020(%%rdx),%%ymm5,%%ymm5		\n\t"\
		"vmovaps	%%ymm6		,0x080(%%rsi)\n\t"\
		"vmovaps	%%ymm7		,0x0a0(%%rsi)\n\t"\
		"vaddpd	%%ymm2		,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	%%ymm3		,%%ymm7,%%ymm7		\n\t"\
		"vsubpd	0x080(%%rsi),%%ymm2,%%ymm2		\n\t"\
		"vsubpd	0x0a0(%%rsi),%%ymm3,%%ymm3		\n\t"\
		"vaddpd	     (%%rsi),%%ymm6,%%ymm6		\n\t"\
		"vaddpd	0x020(%%rsi),%%ymm7,%%ymm7		\n\t"\
		"vmovaps	%%ymm6		,     (%%rsi)\n\t"\
		"vmovaps	%%ymm7		,0x020(%%rsi)\n\t"\
		"vsubpd	0x100(%%rsi),%%ymm6,%%ymm6		\n\t"\
		"vsubpd	0x120(%%rsi),%%ymm7,%%ymm7		\n\t"\
		"vmovaps	%%ymm6		,0x100(%%rsi)\n\t"\
		"vmovaps	%%ymm7		,0x120(%%rsi)\n\t"\
		"vmovaps	%%ymm4		,%%ymm6		\n\t"\
		"vmovaps	%%ymm5		,%%ymm7		\n\t"\
		"vaddpd	%%ymm0		,%%ymm5,%%ymm5		\n\t"\
		"vsubpd	%%ymm7		,%%ymm0,%%ymm0		\n\t"\
		"vaddpd	%%ymm1		,%%ymm4,%%ymm4		\n\t"\
		"vsubpd	%%ymm6		,%%ymm1,%%ymm1		\n\t"\
		"vmovaps	%%ymm5		,%%ymm6		\n\t"\
		"vmovaps	%%ymm1		,%%ymm7		\n\t"\
		"vaddpd	0x140(%%rsi),%%ymm5,%%ymm5		\n\t"\
		"vsubpd	0x160(%%rsi),%%ymm1,%%ymm1		\n\t"\
		"vsubpd	0x140(%%rsi),%%ymm6,%%ymm6		\n\t"\
		"vaddpd	0x160(%%rsi),%%ymm7,%%ymm7		\n\t"\
		"vmovaps	%%ymm5		,0x040(%%rsi)\n\t"\
		"vmovaps	%%ymm1		,0x060(%%rsi)\n\t"\
		"vmovaps	%%ymm6		,0x140(%%rsi)\n\t"\
		"vmovaps	%%ymm7		,0x160(%%rsi)\n\t"\
		"vmovaps	%%ymm2		,%%ymm6		\n\t"\
		"vmovaps	%%ymm3		,%%ymm7		\n\t"\
		"vaddpd	0x1a0(%%rsi),%%ymm2,%%ymm2		\n\t"\
		"vsubpd	0x180(%%rsi),%%ymm3,%%ymm3		\n\t"\
		"vsubpd	0x1a0(%%rsi),%%ymm6,%%ymm6		\n\t"\
		"vaddpd	0x180(%%rsi),%%ymm7,%%ymm7		\n\t"\
		"vmovaps	%%ymm2		,0x080(%%rsi)\n\t"\
		"vmovaps	%%ymm3		,0x0a0(%%rsi)\n\t"\
		"vmovaps	%%ymm6		,0x180(%%rsi)\n\t"\
		"vmovaps	%%ymm7		,0x1a0(%%rsi)\n\t"\
		"vmovaps	%%ymm0		,%%ymm6		\n\t"\
		"vmovaps	%%ymm4		,%%ymm7		\n\t"\
		"vsubpd	0x1e0(%%rsi),%%ymm0,%%ymm0		\n\t"\
		"vsubpd	0x1c0(%%rsi),%%ymm4,%%ymm4		\n\t"\
		"vaddpd	0x1e0(%%rsi),%%ymm6,%%ymm6		\n\t"\
		"vaddpd	0x1c0(%%rsi),%%ymm7,%%ymm7		\n\t"\
		"vmovaps	%%ymm0		,0x0c0(%%rsi)\n\t"\
		"vmovaps	%%ymm4		,0x0e0(%%rsi)\n\t"\
		"vmovaps	%%ymm6		,0x1c0(%%rsi)\n\t"\
		"vmovaps	%%ymm7		,0x1e0(%%rsi)\n\t"\
		"\n\t"\
		"/*...and now do eight radix-4 transforms, including the internal twiddle factors: */\n\t"\
		"\n\t"\
		"/*...Block 1: t00,t10,t20,t30	*/\n\t"\
		"movq	%[__add0],%%rax\n\t"\
		"movslq	%[__p08],%%rbx\n\t"\
		"movq	%[__r00],%%rcx\n\t"\
		"movq	%[__r10],%%rdx\n\t"\
		"movslq	%[__p10],%%rdi	/* rdi will store copy of p10 throughout */\n\t"\
		"shlq	$3,%%rbx\n\t"\
		"shlq	$3,%%rdi\n\t"\
		"addq	%%rax	,%%rbx	/* add1 = add0+p8 */\n\t"\
		"\n\t"\
		"vmovaps	     (%%rdx),%%ymm2	/* t10 */\n\t"\
		"vmovaps	0x400(%%rdx),%%ymm4	/* t30 */\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm3	/* t11 */\n\t"\
		"vmovaps	0x420(%%rdx),%%ymm5	/* t31 */\n\t"\
		"vmovaps	     (%%rcx),%%ymm0	/* t00 */\n\t"\
		"vmovaps	0x400(%%rcx),%%ymm6	/* t20 */\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm1	/* t01 */\n\t"\
		"vmovaps	0x420(%%rcx),%%ymm7	/* t21 */\n\t"\
		"\n\t"\
		"vsubpd	%%ymm2,%%ymm0,%%ymm0	/*~t10=t00-t10*/\n\t"\
		"vsubpd	%%ymm4,%%ymm6,%%ymm6	/*~t30=t20-t30*/\n\t"\
		"vsubpd	%%ymm3,%%ymm1,%%ymm1	/*~t11=t01-t11*/\n\t"\
		"vsubpd	%%ymm5,%%ymm7,%%ymm7	/*~t31=t21-t31*/\n\t"\
		"vaddpd	%%ymm2,%%ymm2,%%ymm2	/*       2*t10*/\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4	/*       2*t30*/\n\t"\
		"vaddpd	%%ymm3,%%ymm3,%%ymm3	/*       2*t11*/\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5	/*       2*t31*/\n\t"\
		"vaddpd	%%ymm0,%%ymm2,%%ymm2	/*~t00=t00+t10*/\n\t"\
		"vaddpd	%%ymm6,%%ymm4,%%ymm4	/*~t20=t20+t30*/\n\t"\
		"vaddpd	%%ymm1,%%ymm3,%%ymm3	/*~t01=t01+t11*/\n\t"\
		"vaddpd	%%ymm7,%%ymm5,%%ymm5	/*~t21=t21+t31*/\n\t"\
		"\n\t"\
		"/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c00): */\n\t"\
		"movq	%[__isrt2],%%rsi\n\t"\
		"vsubpd	%%ymm4,%%ymm2,%%ymm2\n\t"\
		"vsubpd	%%ymm7,%%ymm0,%%ymm0\n\t"\
		"vsubpd	%%ymm5,%%ymm3,%%ymm3\n\t"\
		"vsubpd	%%ymm6,%%ymm1,%%ymm1\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6\n\t"\
		"vaddpd	%%ymm2,%%ymm4,%%ymm4\n\t"\
		"vaddpd	%%ymm0,%%ymm7,%%ymm7\n\t"\
		"vaddpd	%%ymm3,%%ymm5,%%ymm5\n\t"\
		"vaddpd	%%ymm1,%%ymm6,%%ymm6\n\t"\
		"addq	$0x0e0,%%rsi	/* c00 */\n\t"\
		"vmovaps	%%ymm2,     (%%rcx)\n\t"\
		"vmovaps	%%ymm0,0x020(%%rdx)\n\t"\
		"vmovaps	%%ymm3,0x020(%%rcx)\n\t"\
		"vmovaps	%%ymm6,     (%%rdx)\n\t"\
		"vmovaps	%%ymm4,%%ymm2\n\t"\
		"vmovaps	%%ymm7,%%ymm0\n\t"\
		"vmovaps	%%ymm5,%%ymm3\n\t"\
		"vmovaps	%%ymm1,%%ymm6\n\t"\
		"vmulpd	     (%%rsi),%%ymm4,%%ymm4\n\t"\
		"vmulpd	0x040(%%rsi),%%ymm7,%%ymm7\n\t"\
		"vmulpd	     (%%rsi),%%ymm5,%%ymm5\n\t"\
		"vmulpd	0x040(%%rsi),%%ymm1,%%ymm1\n\t"\
		"vmulpd	0x020(%%rsi),%%ymm2,%%ymm2\n\t"\
		"vmulpd	0x060(%%rsi),%%ymm0,%%ymm0\n\t"\
		"vmulpd	0x020(%%rsi),%%ymm3,%%ymm3\n\t"\
		"vmulpd	0x060(%%rsi),%%ymm6,%%ymm6\n\t"\
		"vsubpd	%%ymm2,%%ymm5,%%ymm5\n\t"\
		"vsubpd	%%ymm0,%%ymm1,%%ymm1\n\t"\
		"vaddpd	%%ymm3,%%ymm4,%%ymm4\n\t"\
		"vaddpd	%%ymm6,%%ymm7,%%ymm7\n\t"\
		"vmovaps	%%ymm5,0x020(%%rax)\n\t"\
		"vmovaps	%%ymm1,0x020(%%rbx)\n\t"\
		"vmovaps	%%ymm4,     (%%rax)\n\t"\
		"vmovaps	%%ymm7,     (%%rbx)\n\t"\
		"addq	%%rdi,%%rax\n\t"\
		"addq	%%rdi,%%rbx\n\t"\
		"addq	$0x80,%%rsi\n\t"\
		"vmovaps	     (%%rcx),%%ymm4\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm0\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5\n\t"\
		"vmovaps	     (%%rdx),%%ymm6\n\t"\
		"vmovaps	%%ymm4,%%ymm2\n\t"\
		"vmovaps	%%ymm0,%%ymm1\n\t"\
		"vmovaps	%%ymm5,%%ymm3\n\t"\
		"vmovaps	%%ymm6,%%ymm7\n\t"\
		"vmulpd	     (%%rsi),%%ymm4,%%ymm4\n\t"\
		"vmulpd	0x040(%%rsi),%%ymm0,%%ymm0\n\t"\
		"vmulpd	     (%%rsi),%%ymm5,%%ymm5\n\t"\
		"vmulpd	0x040(%%rsi),%%ymm6,%%ymm6\n\t"\
		"vmulpd	0x020(%%rsi),%%ymm2,%%ymm2\n\t"\
		"vmulpd	0x060(%%rsi),%%ymm1,%%ymm1\n\t"\
		"vmulpd	0x020(%%rsi),%%ymm3,%%ymm3\n\t"\
		"vmulpd	0x060(%%rsi),%%ymm7,%%ymm7\n\t"\
		"vsubpd	%%ymm2,%%ymm5,%%ymm5\n\t"\
		"vsubpd	%%ymm1,%%ymm6,%%ymm6\n\t"\
		"vaddpd	%%ymm3,%%ymm4,%%ymm4\n\t"\
		"vaddpd	%%ymm7,%%ymm0,%%ymm0\n\t"\
		"vmovaps	%%ymm5,0x020(%%rax)\n\t"\
		"vmovaps	%%ymm6,0x020(%%rbx)\n\t"\
		"vmovaps	%%ymm4,     (%%rax)\n\t"\
		"vmovaps	%%ymm0,     (%%rbx)\n\t"\
		"\n\t"\
		"/*...Block 2: t02,t12,t22,t32	*/\n\t"\
		"movq	%[__add0],%%rax\n\t"\
		"movslq	%[__p01],%%rsi\n\t"\
		"movslq	%[__p08],%%rbx\n\t"\
		"shlq	$3,%%rsi\n\t"\
		"shlq	$3,%%rbx\n\t"\
		"addq	%%rsi	,%%rax	/* add0 = &a[j1+p1] */\n\t"\
		"addq	%%rax	,%%rbx	/* add1 = add0+p8 */\n\t"\
		"movq	%[__r02],%%rcx\n\t"\
		"movq	%[__r12],%%rdx\n\t"\
		"movq	%[__isrt2],%%rsi\n\t"\
		"addq	$0x400,%%rcx\n\t"\
		"addq	$0x400,%%rdx\n\t"\
		"\n\t"\
		"addq	$0x060,%%rsi	/* cc1 */\n\t"\
		"vmovaps	     (%%rcx),%%ymm4	/* t22 */\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5	/* t23 */\n\t"\
		"vmovaps	     (%%rsi),%%ymm2	/* c32_1 */\n\t"\
		"vmovaps	0x020(%%rsi),%%ymm3	/* s32_1 */\n\t"\
		"vmovaps	%%ymm4,%%ymm6	/* ymm2 <- cpy t22 */\n\t"\
		"vmovaps	%%ymm5,%%ymm7	/* ymm3 <- cpy t23 */\n\t"\
		"\n\t"\
		"vmulpd	%%ymm2,%%ymm4,%%ymm4	/* t22*c32_1 */\n\t"\
		"vmulpd	%%ymm2,%%ymm5,%%ymm5	/* t23*c32_1 */\n\t"\
		"vmulpd	%%ymm3,%%ymm6,%%ymm6	/* t22*s32_1 */\n\t"\
		"vmovaps	     (%%rdx),%%ymm0	/* t32 */\n\t"\
		"vmulpd	%%ymm3,%%ymm7,%%ymm7	/* t23*s32_1 */\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm1	/* t33 */\n\t"\
		"addq	$0x040,%%rsi	/* cc3 */\n\t"\
		"vmovaps	     (%%rsi),%%ymm2	/* c32_3 */\n\t"\
		"vmovaps	0x020(%%rsi),%%ymm3	/* s32_3 */\n\t"\
		"vsubpd	%%ymm6,%%ymm5,%%ymm5	/* ymm5 <-~t23 */\n\t"\
		"vmovaps	%%ymm0,%%ymm6	/* ymm6 <- cpy t32 */\n\t"\
		"vaddpd	%%ymm7,%%ymm4,%%ymm4	/* ymm4 <-~t22 */\n\t"\
		"vmovaps	%%ymm1,%%ymm7	/* ymm7 <- cpy t33 */\n\t"\
		"\n\t"\
		"vmulpd	%%ymm2,%%ymm0,%%ymm0	/* t32*c32_3 */\n\t"\
		"vmulpd	%%ymm2,%%ymm1,%%ymm1	/* t33*c32_3 */\n\t"\
		"vmulpd	%%ymm3,%%ymm6,%%ymm6	/* t32*s32_3 */\n\t"\
		"vmulpd	%%ymm3,%%ymm7,%%ymm7	/* t33*s32_3 */\n\t"\
		"vsubpd	%%ymm6,%%ymm1,%%ymm1	/* ymm1 <- it */\n\t"\
		"vaddpd	%%ymm7,%%ymm0,%%ymm0	/* ymm0 <- rt */\n\t"\
		"\n\t"\
		"vmovaps	%%ymm5,%%ymm7	/* ymm7 <- cpy~t23*/\n\t"\
		"vmovaps	%%ymm4,%%ymm6	/* ymm6 <- cpy~t22*/\n\t"\
		"\n\t"\
		"vaddpd	%%ymm0,%%ymm4,%%ymm4	/* ~t22 <- t22+rt */\n\t"\
		"vaddpd	%%ymm1,%%ymm5,%%ymm5	/* ~t23 <- t23+it */\n\t"\
		"vsubpd	%%ymm0,%%ymm6,%%ymm6	/* ~t32 <- t22-rt */\n\t"\
		"vsubpd	%%ymm1,%%ymm7,%%ymm7	/* ~t33 <- t23-it */\n\t"\
		"\n\t"\
		"subq	$0x400,%%rcx\n\t"\
		"subq	$0x400,%%rdx\n\t"\
		"subq	$0x080,%%rsi	/* cc0 */\n\t"\
		"vmovaps	     (%%rdx),%%ymm1	/* t12 */\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm3	/* t13 */\n\t"\
		"vmovaps	0x020(%%rsi),%%ymm2	/* s */\n\t"\
		"vmovaps	%%ymm1,%%ymm0	/* cpy t12 */\n\t"\
		"vmulpd	%%ymm2,%%ymm1,%%ymm1	/* t12*s */\n\t"\
		"vmulpd	%%ymm3,%%ymm2,%%ymm2	/* t13*s */\n\t"\
		"vmulpd	(%%rsi),%%ymm0,%%ymm0	/* t12*c */\n\t"\
		"vmulpd	(%%rsi),%%ymm3,%%ymm3	/* t13*c */\n\t"\
		"\n\t"\
		"vaddpd	%%ymm0,%%ymm2,%%ymm2	/* rt =t12*c + t13*s */\n\t"\
		"vsubpd	%%ymm1,%%ymm3,%%ymm3	/* it =t13*c - t12*s */\n\t"\
		"\n\t"\
		"vmovaps	     (%%rcx),%%ymm0	/* t02 */\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm1	/* t03 */\n\t"\
		"vsubpd	%%ymm2,%%ymm0,%%ymm0	/*~t12 <- t02- rt */\n\t"\
		"vsubpd	%%ymm3,%%ymm1,%%ymm1	/*~t13 <- t03- it */\n\t"\
		"vaddpd	%%ymm2,%%ymm2,%%ymm2	/*          2* rt */\n\t"\
		"vaddpd	%%ymm3,%%ymm3,%%ymm3	/*          2* it */\n\t"\
		"vaddpd	%%ymm0,%%ymm2,%%ymm2	/*~t02 <- t02+ rt */\n\t"\
		"vaddpd	%%ymm1,%%ymm3,%%ymm3	/*~t03 <- t03+ it */\n\t"\
		"\n\t"\
		"/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c01): */\n\t"\
		"addq	$0x4c0,%%rsi	/* c01 */\n\t"\
		"vsubpd	%%ymm4,%%ymm2,%%ymm2\n\t"\
		"vsubpd	%%ymm7,%%ymm0,%%ymm0\n\t"\
		"vsubpd	%%ymm5,%%ymm3,%%ymm3\n\t"\
		"vsubpd	%%ymm6,%%ymm1,%%ymm1\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6\n\t"\
		"vaddpd	%%ymm2,%%ymm4,%%ymm4\n\t"\
		"vaddpd	%%ymm0,%%ymm7,%%ymm7\n\t"\
		"vaddpd	%%ymm3,%%ymm5,%%ymm5\n\t"\
		"vaddpd	%%ymm1,%%ymm6,%%ymm6\n\t"\
		"vmovaps	%%ymm2,     (%%rcx)\n\t"\
		"vmovaps	%%ymm0,0x020(%%rdx)\n\t"\
		"vmovaps	%%ymm3,0x020(%%rcx)\n\t"\
		"vmovaps	%%ymm6,     (%%rdx)\n\t"\
		"vmovaps	%%ymm4,%%ymm2\n\t"\
		"vmovaps	%%ymm7,%%ymm0\n\t"\
		"vmovaps	%%ymm5,%%ymm3\n\t"\
		"vmovaps	%%ymm1,%%ymm6\n\t"\
		"vmulpd	     (%%rsi),%%ymm4,%%ymm4\n\t"\
		"vmulpd	0x040(%%rsi),%%ymm7,%%ymm7\n\t"\
		"vmulpd	     (%%rsi),%%ymm5,%%ymm5\n\t"\
		"vmulpd	0x040(%%rsi),%%ymm1,%%ymm1\n\t"\
		"vmulpd	0x020(%%rsi),%%ymm2,%%ymm2\n\t"\
		"vmulpd	0x060(%%rsi),%%ymm0,%%ymm0\n\t"\
		"vmulpd	0x020(%%rsi),%%ymm3,%%ymm3\n\t"\
		"vmulpd	0x060(%%rsi),%%ymm6,%%ymm6\n\t"\
		"vsubpd	%%ymm2,%%ymm5,%%ymm5\n\t"\
		"vsubpd	%%ymm0,%%ymm1,%%ymm1\n\t"\
		"vaddpd	%%ymm3,%%ymm4,%%ymm4\n\t"\
		"vaddpd	%%ymm6,%%ymm7,%%ymm7\n\t"\
		"vmovaps	%%ymm5,0x020(%%rax)\n\t"\
		"vmovaps	%%ymm1,0x020(%%rbx)\n\t"\
		"vmovaps	%%ymm4,     (%%rax)\n\t"\
		"vmovaps	%%ymm7,     (%%rbx)\n\t"\
		"addq	%%rdi,%%rax\n\t"\
		"addq	%%rdi,%%rbx\n\t"\
		"addq	$0x80,%%rsi\n\t"\
		"vmovaps	     (%%rcx),%%ymm4\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm0\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5\n\t"\
		"vmovaps	     (%%rdx),%%ymm6\n\t"\
		"vmovaps	%%ymm4,%%ymm2\n\t"\
		"vmovaps	%%ymm0,%%ymm1\n\t"\
		"vmovaps	%%ymm5,%%ymm3\n\t"\
		"vmovaps	%%ymm6,%%ymm7\n\t"\
		"vmulpd	     (%%rsi),%%ymm4,%%ymm4\n\t"\
		"vmulpd	0x040(%%rsi),%%ymm0,%%ymm0\n\t"\
		"vmulpd	     (%%rsi),%%ymm5,%%ymm5\n\t"\
		"vmulpd	0x040(%%rsi),%%ymm6,%%ymm6\n\t"\
		"vmulpd	0x020(%%rsi),%%ymm2,%%ymm2\n\t"\
		"vmulpd	0x060(%%rsi),%%ymm1,%%ymm1\n\t"\
		"vmulpd	0x020(%%rsi),%%ymm3,%%ymm3\n\t"\
		"vmulpd	0x060(%%rsi),%%ymm7,%%ymm7\n\t"\
		"vsubpd	%%ymm2,%%ymm5,%%ymm5\n\t"\
		"vsubpd	%%ymm1,%%ymm6,%%ymm6\n\t"\
		"vaddpd	%%ymm3,%%ymm4,%%ymm4\n\t"\
		"vaddpd	%%ymm7,%%ymm0,%%ymm0\n\t"\
		"vmovaps	%%ymm5,0x020(%%rax)\n\t"\
		"vmovaps	%%ymm6,0x020(%%rbx)\n\t"\
		"vmovaps	%%ymm4,     (%%rax)\n\t"\
		"vmovaps	%%ymm0,     (%%rbx)\n\t"\
		"\n\t"\
		"/*...Block 3: t04,t14,t24,t34	*/\n\t"\
		"movq	%[__add0],%%rax\n\t"\
		"movslq	%[__p02],%%rsi\n\t"\
		"movslq	%[__p08],%%rbx\n\t"\
		"shlq	$3,%%rsi\n\t"\
		"shlq	$3,%%rbx\n\t"\
		"addq	%%rsi	,%%rax	/* add0 = &a[j1+p2] */\n\t"\
		"addq	%%rax	,%%rbx	/* add1 = add0+p8 */\n\t"\
		"movq	%[__r04],%%rcx\n\t"\
		"movq	%[__r14],%%rdx\n\t"\
		"movq	%[__isrt2],%%rsi\n\t"\
		"addq	$0x400,%%rcx	/* t24 */\n\t"\
		"addq	$0x400,%%rdx	/* t25 */\n\t"\
		"addq	$0x020,%%rsi	/* cc0 */\n\t"\
		"\n\t"\
		"vmovaps	     (%%rcx),%%ymm4	/* t24 */\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5	/* t25 */\n\t"\
		"vmovaps	     (%%rsi),%%ymm2	/* c */\n\t"\
		"vmovaps	0x020(%%rsi),%%ymm3	/* s */\n\t"\
		"vmovaps	%%ymm4,%%ymm6	/* ymm2 <- cpy t24 */\n\t"\
		"vmovaps	%%ymm5,%%ymm7	/* ymm3 <- cpy t25 */\n\t"\
		"\n\t"\
		"vmulpd	%%ymm2,%%ymm4,%%ymm4	/* t24*c */\n\t"\
		"vmulpd	%%ymm2,%%ymm5,%%ymm5	/* t25*c */\n\t"\
		"vmulpd	%%ymm3,%%ymm6,%%ymm6	/* t24*s */\n\t"\
		"vmovaps	     (%%rdx),%%ymm0	/* t34 */\n\t"\
		"vmulpd	%%ymm3,%%ymm7,%%ymm7	/* t25*s */\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm1	/* t35 */\n\t"\
		"vsubpd	%%ymm6,%%ymm5,%%ymm5	/* ymm1 <-~t25 */\n\t"\
		"vmovaps	%%ymm0,%%ymm6	/* ymm6 <- cpy t34 */\n\t"\
		"vaddpd	%%ymm7,%%ymm4,%%ymm4	/* ymm0 <-~t24 */\n\t"\
		"vmovaps	%%ymm1,%%ymm7	/* ymm7 <- cpy t35 */\n\t"\
		"\n\t"\
		"vmulpd	%%ymm3,%%ymm0,%%ymm0	/* t34*s */\n\t"\
		"vmulpd	%%ymm3,%%ymm1,%%ymm1	/* t35*s */\n\t"\
		"vmulpd	%%ymm2,%%ymm6,%%ymm6	/* t34*c */\n\t"\
		"vmulpd	%%ymm2,%%ymm7,%%ymm7	/* t35*c */\n\t"\
		"vsubpd	%%ymm6,%%ymm1,%%ymm1	/* ymm5 <- it */\n\t"\
		"vaddpd	%%ymm7,%%ymm0,%%ymm0	/* ymm4 <- rt */\n\t"\
		"\n\t"\
		"vmovaps	%%ymm5,%%ymm7	/* ymm3 <- cpy~t25*/\n\t"\
		"vmovaps	%%ymm4,%%ymm6	/* ymm2 <- cpy~t24*/\n\t"\
		"\n\t"\
		"vaddpd	%%ymm0,%%ymm4,%%ymm4	/* ~t24 <- t24+rt */\n\t"\
		"vaddpd	%%ymm1,%%ymm5,%%ymm5	/* ~t25 <- t25+it */\n\t"\
		"vsubpd	%%ymm0,%%ymm6,%%ymm6	/* ~t34 <- t24-rt */\n\t"\
		"vsubpd	%%ymm1,%%ymm7,%%ymm7	/* ~t35 <- t25-it */\n\t"\
		"\n\t"\
		"subq	$0x400,%%rcx\n\t"\
		"subq	$0x400,%%rdx\n\t"\
		"subq	$0x20,%%rsi	/* isrt2 */\n\t"\
		"vmovaps	     (%%rdx),%%ymm2	/* t14 */\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm3	/* t15 */\n\t"\
		"vmovaps	(%%rsi),%%ymm1	/* isrt2 */\n\t"\
		"vmovaps	%%ymm3,%%ymm0	/* cpy t15 */\n\t"\
		"vsubpd	%%ymm2,%%ymm3,%%ymm3	/*~t15=t15-t14 */\n\t"\
		"vaddpd	%%ymm0,%%ymm2,%%ymm2	/*~t14=t14+t15 */\n\t"\
		"vmulpd	%%ymm1,%%ymm2,%%ymm2	/* rt */\n\t"\
		"vmulpd	%%ymm1,%%ymm3,%%ymm3	/* it */\n\t"\
		"\n\t"\
		"vmovaps	     (%%rcx),%%ymm0	/* t04 */\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm1	/* t05 */\n\t"\
		"vsubpd	%%ymm2,%%ymm0,%%ymm0	/*~t14 <- t04- rt */\n\t"\
		"vsubpd	%%ymm3,%%ymm1,%%ymm1	/*~t15 <- t05- it */\n\t"\
		"vaddpd	%%ymm2,%%ymm2,%%ymm2	/*          2* rt */\n\t"\
		"vaddpd	%%ymm3,%%ymm3,%%ymm3	/*          2* it */\n\t"\
		"vaddpd	%%ymm0,%%ymm2,%%ymm2	/*~t04 <- t04+ rt */\n\t"\
		"vaddpd	%%ymm1,%%ymm3,%%ymm3	/*~t05 <- t05+ it */\n\t"\
		"\n\t"\
		"/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c02): */\n\t"\
		"addq	$0x2e0,%%rsi	/* c02 */\n\t"\
		"vsubpd	%%ymm4,%%ymm2,%%ymm2\n\t"\
		"vsubpd	%%ymm7,%%ymm0,%%ymm0\n\t"\
		"vsubpd	%%ymm5,%%ymm3,%%ymm3\n\t"\
		"vsubpd	%%ymm6,%%ymm1,%%ymm1\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6\n\t"\
		"vaddpd	%%ymm2,%%ymm4,%%ymm4\n\t"\
		"vaddpd	%%ymm0,%%ymm7,%%ymm7\n\t"\
		"vaddpd	%%ymm3,%%ymm5,%%ymm5\n\t"\
		"vaddpd	%%ymm1,%%ymm6,%%ymm6\n\t"\
		"vmovaps	%%ymm2,     (%%rcx)\n\t"\
		"vmovaps	%%ymm0,0x020(%%rdx)\n\t"\
		"vmovaps	%%ymm3,0x020(%%rcx)\n\t"\
		"vmovaps	%%ymm6,     (%%rdx)\n\t"\
		"vmovaps	%%ymm4,%%ymm2\n\t"\
		"vmovaps	%%ymm7,%%ymm0\n\t"\
		"vmovaps	%%ymm5,%%ymm3\n\t"\
		"vmovaps	%%ymm1,%%ymm6\n\t"\
		"vmulpd	     (%%rsi),%%ymm4,%%ymm4\n\t"\
		"vmulpd	0x040(%%rsi),%%ymm7,%%ymm7\n\t"\
		"vmulpd	     (%%rsi),%%ymm5,%%ymm5\n\t"\
		"vmulpd	0x040(%%rsi),%%ymm1,%%ymm1\n\t"\
		"vmulpd	0x020(%%rsi),%%ymm2,%%ymm2\n\t"\
		"vmulpd	0x060(%%rsi),%%ymm0,%%ymm0\n\t"\
		"vmulpd	0x020(%%rsi),%%ymm3,%%ymm3\n\t"\
		"vmulpd	0x060(%%rsi),%%ymm6,%%ymm6\n\t"\
		"vsubpd	%%ymm2,%%ymm5,%%ymm5\n\t"\
		"vsubpd	%%ymm0,%%ymm1,%%ymm1\n\t"\
		"vaddpd	%%ymm3,%%ymm4,%%ymm4\n\t"\
		"vaddpd	%%ymm6,%%ymm7,%%ymm7\n\t"\
		"vmovaps	%%ymm5,0x020(%%rax)\n\t"\
		"vmovaps	%%ymm1,0x020(%%rbx)\n\t"\
		"vmovaps	%%ymm4,     (%%rax)\n\t"\
		"vmovaps	%%ymm7,     (%%rbx)\n\t"\
		"addq	%%rdi,%%rax\n\t"\
		"addq	%%rdi,%%rbx\n\t"\
		"addq	$0x80,%%rsi\n\t"\
		"vmovaps	     (%%rcx),%%ymm4\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm0\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5\n\t"\
		"vmovaps	     (%%rdx),%%ymm6\n\t"\
		"vmovaps	%%ymm4,%%ymm2\n\t"\
		"vmovaps	%%ymm0,%%ymm1\n\t"\
		"vmovaps	%%ymm5,%%ymm3\n\t"\
		"vmovaps	%%ymm6,%%ymm7\n\t"\
		"vmulpd	     (%%rsi),%%ymm4,%%ymm4\n\t"\
		"vmulpd	0x040(%%rsi),%%ymm0,%%ymm0\n\t"\
		"vmulpd	     (%%rsi),%%ymm5,%%ymm5\n\t"\
		"vmulpd	0x040(%%rsi),%%ymm6,%%ymm6\n\t"\
		"vmulpd	0x020(%%rsi),%%ymm2,%%ymm2\n\t"\
		"vmulpd	0x060(%%rsi),%%ymm1,%%ymm1\n\t"\
		"vmulpd	0x020(%%rsi),%%ymm3,%%ymm3\n\t"\
		"vmulpd	0x060(%%rsi),%%ymm7,%%ymm7\n\t"\
		"vsubpd	%%ymm2,%%ymm5,%%ymm5\n\t"\
		"vsubpd	%%ymm1,%%ymm6,%%ymm6\n\t"\
		"vaddpd	%%ymm3,%%ymm4,%%ymm4\n\t"\
		"vaddpd	%%ymm7,%%ymm0,%%ymm0\n\t"\
		"vmovaps	%%ymm5,0x020(%%rax)\n\t"\
		"vmovaps	%%ymm6,0x020(%%rbx)\n\t"\
		"vmovaps	%%ymm4,     (%%rax)\n\t"\
		"vmovaps	%%ymm0,     (%%rbx)\n\t"\
		"\n\t"\
		"/*...Block 4: t06,t16,t26,t36	*/\n\t"\
		"movq	%[__add0],%%rax\n\t"\
		"movslq	%[__p03],%%rsi\n\t"\
		"movslq	%[__p08],%%rbx\n\t"\
		"shlq	$3,%%rsi\n\t"\
		"shlq	$3,%%rbx\n\t"\
		"addq	%%rsi	,%%rax	/* add0 = &a[j1+p3] */\n\t"\
		"addq	%%rax	,%%rbx	/* add1 = add0+p8 */\n\t"\
		"movq	%[__r06],%%rcx\n\t"\
		"movq	%[__r16],%%rdx\n\t"\
		"movq	%[__isrt2],%%rsi\n\t"\
		"addq	$0x400,%%rcx\n\t"\
		"addq	$0x400,%%rdx\n\t"\
		"addq	$0x0a0,%%rsi	/* cc3 */\n\t"\
		"\n\t"\
		"vmovaps	     (%%rcx),%%ymm4	/* t26 */\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5	/* t27 */\n\t"\
		"vmovaps	     (%%rsi),%%ymm2	/* c32_3 */\n\t"\
		"vmovaps	0x020(%%rsi),%%ymm3	/* s32_3 */\n\t"\
		"vmovaps	%%ymm4,%%ymm6	/* ymm2 <- cpy t26 */\n\t"\
		"vmovaps	%%ymm5,%%ymm7	/* ymm3 <- cpy t27 */\n\t"\
		"\n\t"\
		"vmulpd	%%ymm2,%%ymm4,%%ymm4	/* t26*c32_3 */\n\t"\
		"vmulpd	%%ymm2,%%ymm5,%%ymm5	/* t27*c32_3 */\n\t"\
		"vmulpd	%%ymm3,%%ymm6,%%ymm6	/* t26*s32_3 */\n\t"\
		"vmovaps	     (%%rdx),%%ymm0	/* t36 */\n\t"\
		"vmulpd	%%ymm3,%%ymm7,%%ymm7	/* t27*s32_3 */\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm1	/* t37 */\n\t"\
		"subq	$0x40,%%rsi	/* cc1 */\n\t"\
		"vmovaps	     (%%rsi),%%ymm3	/* c32_1 */\n\t"\
		"vmovaps	0x020(%%rsi),%%ymm2	/* s32_1 */\n\t"\
		"vsubpd	%%ymm6,%%ymm5,%%ymm5	/* ymm5 <-~t27 */\n\t"\
		"vmovaps	%%ymm0,%%ymm6	/* ymm6 <- cpy t36 */\n\t"\
		"vaddpd	%%ymm7,%%ymm4,%%ymm4	/* ymm4 <-~t26 */\n\t"\
		"vmovaps	%%ymm1,%%ymm7	/* ymm7 <- cpy t37 */\n\t"\
		"\n\t"\
		"vmulpd	%%ymm2,%%ymm0,%%ymm0	/* t36*s32_1 */\n\t"\
		"vmulpd	%%ymm2,%%ymm1,%%ymm1	/* t37*s32_1 */\n\t"\
		"vmulpd	%%ymm3,%%ymm6,%%ymm6	/* t36*c32_1 */\n\t"\
		"vmulpd	%%ymm3,%%ymm7,%%ymm7	/* t37*c32_1 */\n\t"\
		"vaddpd	%%ymm6,%%ymm1,%%ymm1	/* ymm1 <- it */\n\t"\
		"vsubpd	%%ymm7,%%ymm0,%%ymm0	/* ymm0 <- rt */\n\t"\
		"\n\t"\
		"vmovaps	%%ymm5,%%ymm7	/* ymm7 <- cpy~t27*/\n\t"\
		"vmovaps	%%ymm4,%%ymm6	/* ymm6 <- cpy~t26*/\n\t"\
		"\n\t"\
		"vaddpd	%%ymm0,%%ymm6,%%ymm6	/* ~t36 <- t26+rt */\n\t"\
		"vaddpd	%%ymm1,%%ymm7,%%ymm7	/* ~t37 <- t27+it */\n\t"\
		"vsubpd	%%ymm0,%%ymm4,%%ymm4	/* ~t26 <- t26-rt */\n\t"\
		"vsubpd	%%ymm1,%%ymm5,%%ymm5	/* ~t27 <- t27-it */\n\t"\
		"\n\t"\
		"subq	$0x400,%%rcx\n\t"\
		"subq	$0x400,%%rdx\n\t"\
		"subq	$0x40,%%rsi	/* cc0 */\n\t"\
		"vmovaps	     (%%rdx),%%ymm2	/* t16 */\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm0	/* t17 */\n\t"\
		"vmovaps	0x020(%%rsi),%%ymm3	/* s */\n\t"\
		"vmovaps	%%ymm2,%%ymm1	/* cpy t16 */\n\t"\
		"vmulpd	%%ymm3,%%ymm2,%%ymm2	/* t16*s */\n\t"\
		"vmulpd	%%ymm0,%%ymm3,%%ymm3	/* s*t17 */\n\t"\
		"vmulpd	(%%rsi),%%ymm1,%%ymm1	/* t16*c */\n\t"\
		"vmulpd	(%%rsi),%%ymm0,%%ymm0	/* t17*c */\n\t"\
		"\n\t"\
		"vaddpd	%%ymm0,%%ymm2,%%ymm2	/* rt =t16*s - t17*c */\n\t"\
		"vsubpd	%%ymm1,%%ymm3,%%ymm3	/* it =t17*s + t16*c */\n\t"\
		"\n\t"\
		"vmovaps	     (%%rcx),%%ymm0	/* t06 */\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm1	/* t07 */\n\t"\
		"vsubpd	%%ymm2,%%ymm0,%%ymm0	/*~t16 <- t06- rt */\n\t"\
		"vsubpd	%%ymm3,%%ymm1,%%ymm1	/*~t17 <- t07- it */\n\t"\
		"vaddpd	%%ymm2,%%ymm2,%%ymm2	/*          2* rt */\n\t"\
		"vaddpd	%%ymm3,%%ymm3,%%ymm3	/*          2* it */\n\t"\
		"vaddpd	%%ymm0,%%ymm2,%%ymm2	/*~t06 <- t06+ rt */\n\t"\
		"vaddpd	%%ymm1,%%ymm3,%%ymm3	/*~t07 <- t07+ it */\n\t"\
		"\n\t"\
		"/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c03): */\n\t"\
		"addq	$0x6c0,%%rsi	/* c03 */\n\t"\
		"vsubpd	%%ymm4,%%ymm2,%%ymm2\n\t"\
		"vsubpd	%%ymm7,%%ymm0,%%ymm0\n\t"\
		"vsubpd	%%ymm5,%%ymm3,%%ymm3\n\t"\
		"vsubpd	%%ymm6,%%ymm1,%%ymm1\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6\n\t"\
		"vaddpd	%%ymm2,%%ymm4,%%ymm4\n\t"\
		"vaddpd	%%ymm0,%%ymm7,%%ymm7\n\t"\
		"vaddpd	%%ymm3,%%ymm5,%%ymm5\n\t"\
		"vaddpd	%%ymm1,%%ymm6,%%ymm6\n\t"\
		"vmovaps	%%ymm2,     (%%rcx)\n\t"\
		"vmovaps	%%ymm0,0x020(%%rdx)\n\t"\
		"vmovaps	%%ymm3,0x020(%%rcx)\n\t"\
		"vmovaps	%%ymm6,     (%%rdx)\n\t"\
		"vmovaps	%%ymm4,%%ymm2\n\t"\
		"vmovaps	%%ymm7,%%ymm0\n\t"\
		"vmovaps	%%ymm5,%%ymm3\n\t"\
		"vmovaps	%%ymm1,%%ymm6\n\t"\
		"vmulpd	     (%%rsi),%%ymm4,%%ymm4\n\t"\
		"vmulpd	0x040(%%rsi),%%ymm7,%%ymm7\n\t"\
		"vmulpd	     (%%rsi),%%ymm5,%%ymm5\n\t"\
		"vmulpd	0x040(%%rsi),%%ymm1,%%ymm1\n\t"\
		"vmulpd	0x020(%%rsi),%%ymm2,%%ymm2\n\t"\
		"vmulpd	0x060(%%rsi),%%ymm0,%%ymm0\n\t"\
		"vmulpd	0x020(%%rsi),%%ymm3,%%ymm3\n\t"\
		"vmulpd	0x060(%%rsi),%%ymm6,%%ymm6\n\t"\
		"vsubpd	%%ymm2,%%ymm5,%%ymm5\n\t"\
		"vsubpd	%%ymm0,%%ymm1,%%ymm1\n\t"\
		"vaddpd	%%ymm3,%%ymm4,%%ymm4\n\t"\
		"vaddpd	%%ymm6,%%ymm7,%%ymm7\n\t"\
		"vmovaps	%%ymm5,0x020(%%rax)\n\t"\
		"vmovaps	%%ymm1,0x020(%%rbx)\n\t"\
		"vmovaps	%%ymm4,     (%%rax)\n\t"\
		"vmovaps	%%ymm7,     (%%rbx)\n\t"\
		"addq	%%rdi,%%rax\n\t"\
		"addq	%%rdi,%%rbx\n\t"\
		"addq	$0x80,%%rsi	/* c0B */\n\t"\
		"vmovaps	     (%%rcx),%%ymm4\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm0\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5\n\t"\
		"vmovaps	     (%%rdx),%%ymm6\n\t"\
		"vmovaps	%%ymm4,%%ymm2\n\t"\
		"vmovaps	%%ymm0,%%ymm1\n\t"\
		"vmovaps	%%ymm5,%%ymm3\n\t"\
		"vmovaps	%%ymm6,%%ymm7\n\t"\
		"vmulpd	     (%%rsi),%%ymm4,%%ymm4\n\t"\
		"vmulpd	0x040(%%rsi),%%ymm0,%%ymm0\n\t"\
		"vmulpd	     (%%rsi),%%ymm5,%%ymm5\n\t"\
		"vmulpd	0x040(%%rsi),%%ymm6,%%ymm6\n\t"\
		"vmulpd	0x020(%%rsi),%%ymm2,%%ymm2\n\t"\
		"vmulpd	0x060(%%rsi),%%ymm1,%%ymm1\n\t"\
		"vmulpd	0x020(%%rsi),%%ymm3,%%ymm3\n\t"\
		"vmulpd	0x060(%%rsi),%%ymm7,%%ymm7\n\t"\
		"vsubpd	%%ymm2,%%ymm5,%%ymm5\n\t"\
		"vsubpd	%%ymm1,%%ymm6,%%ymm6\n\t"\
		"vaddpd	%%ymm3,%%ymm4,%%ymm4\n\t"\
		"vaddpd	%%ymm7,%%ymm0,%%ymm0\n\t"\
		"vmovaps	%%ymm5,0x020(%%rax)\n\t"\
		"vmovaps	%%ymm6,0x020(%%rbx)\n\t"\
		"vmovaps	%%ymm4,     (%%rax)\n\t"\
		"vmovaps	%%ymm0,     (%%rbx)\n\t"\
		"\n\t"\
		"/*...Block 5: t08,t18,t28,t38	*/\n\t"\
		"movq	%[__add0],%%rax\n\t"\
		"movslq	%[__p04],%%rsi\n\t"\
		"movslq	%[__p08],%%rbx\n\t"\
		"shlq	$3,%%rsi\n\t"\
		"shlq	$3,%%rbx\n\t"\
		"addq	%%rsi	,%%rax	/* add0 = &a[j1+p4] */\n\t"\
		"addq	%%rax	,%%rbx	/* add1 = add0+p8 */\n\t"\
		"movq	%[__r08],%%rcx\n\t"\
		"movq	%[__r18],%%rdx\n\t"\
		"movq	%[__isrt2],%%rsi\n\t"\
		"addq	$0x400,%%rcx\n\t"\
		"addq	$0x400,%%rdx\n\t"\
		"vmovaps	(%%rsi),%%ymm2	/* isrt2 */\n\t"\
		"\n\t"\
		"vmovaps	     (%%rcx),%%ymm4	/* t28 */\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5	/* t29 */\n\t"\
		"vmovaps	     (%%rdx),%%ymm6	/* t38 */\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm7	/* t39 */\n\t"\
		"subq	$0x400,%%rcx\n\t"\
		"subq	$0x400,%%rdx\n\t"\
		"vmulpd	%%ymm2,%%ymm4,%%ymm4\n\t"\
		"vmulpd	%%ymm2,%%ymm5,%%ymm5\n\t"\
		"vmulpd	%%ymm2,%%ymm6,%%ymm6\n\t"\
		"vmulpd	%%ymm2,%%ymm7,%%ymm7\n\t"\
		"\n\t"\
		"vsubpd	%%ymm4,%%ymm5,%%ymm5	/*~t29=t29-t28*/\n\t"\
		"vmovaps	     (%%rcx),%%ymm0	/* t08 */\n\t"\
		"vsubpd	%%ymm7,%%ymm6,%%ymm6	/* rt =t38-t39*/\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm2	/* t19 */\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4	/*       2*t28*/\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm3	/* t09 */\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7	/*       2*t39*/\n\t"\
		"vmovaps	     (%%rdx),%%ymm1	/* t18 */\n\t"\
		"vaddpd	%%ymm5,%%ymm4,%%ymm4	/*~t28=t28+t29*/\n\t"\
		"vaddpd	%%ymm6,%%ymm7,%%ymm7	/* it =t39+t38*/\n\t"\
		"\n\t"\
		"vsubpd	%%ymm6,%%ymm4,%%ymm4	/*~t28=t28-rt */\n\t"\
		"vsubpd	%%ymm2,%%ymm0,%%ymm0	/*~t18=t08-t19*/\n\t"\
		"vsubpd	%%ymm7,%%ymm5,%%ymm5	/*~t29=t29-it */\n\t"\
		"vsubpd	%%ymm1,%%ymm3,%%ymm3	/*~t09=t09-t18*/\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6	/*       2*rt */\n\t"\
		"vaddpd	%%ymm2,%%ymm2,%%ymm2	/*       2*t08*/\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7	/*       2*it */\n\t"\
		"vaddpd	%%ymm1,%%ymm1,%%ymm1	/*       2*t09*/\n\t"\
		"vaddpd	%%ymm4,%%ymm6,%%ymm6	/*~t38=t28+rt */\n\t"\
		"vaddpd	%%ymm0,%%ymm2,%%ymm2	/*~t08=t19+t08*/\n\t"\
		"vaddpd	%%ymm5,%%ymm7,%%ymm7	/*~t39=t29+it */\n\t"\
		"vaddpd	%%ymm3,%%ymm1,%%ymm1	/*~t19=t18+t09*/\n\t"\
		"\n\t"\
		"/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c04): */\n\t"\
		"addq	$0x1e0,%%rsi	/* c04 */\n\t"\
		"vsubpd	%%ymm4,%%ymm2,%%ymm2\n\t"\
		"vsubpd	%%ymm7,%%ymm0,%%ymm0\n\t"\
		"vsubpd	%%ymm5,%%ymm3,%%ymm3\n\t"\
		"vsubpd	%%ymm6,%%ymm1,%%ymm1\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6\n\t"\
		"vaddpd	%%ymm2,%%ymm4,%%ymm4\n\t"\
		"vaddpd	%%ymm0,%%ymm7,%%ymm7\n\t"\
		"vaddpd	%%ymm3,%%ymm5,%%ymm5\n\t"\
		"vaddpd	%%ymm1,%%ymm6,%%ymm6\n\t"\
		"vmovaps	%%ymm2,     (%%rcx)\n\t"\
		"vmovaps	%%ymm0,0x020(%%rdx)\n\t"\
		"vmovaps	%%ymm3,0x020(%%rcx)\n\t"\
		"vmovaps	%%ymm6,     (%%rdx)\n\t"\
		"vmovaps	%%ymm4,%%ymm2\n\t"\
		"vmovaps	%%ymm7,%%ymm0\n\t"\
		"vmovaps	%%ymm5,%%ymm3\n\t"\
		"vmovaps	%%ymm1,%%ymm6\n\t"\
		"vmulpd	     (%%rsi),%%ymm4,%%ymm4\n\t"\
		"vmulpd	0x040(%%rsi),%%ymm7,%%ymm7\n\t"\
		"vmulpd	     (%%rsi),%%ymm5,%%ymm5\n\t"\
		"vmulpd	0x040(%%rsi),%%ymm1,%%ymm1\n\t"\
		"vmulpd	0x020(%%rsi),%%ymm2,%%ymm2\n\t"\
		"vmulpd	0x060(%%rsi),%%ymm0,%%ymm0\n\t"\
		"vmulpd	0x020(%%rsi),%%ymm3,%%ymm3\n\t"\
		"vmulpd	0x060(%%rsi),%%ymm6,%%ymm6\n\t"\
		"vsubpd	%%ymm2,%%ymm5,%%ymm5\n\t"\
		"vsubpd	%%ymm0,%%ymm1,%%ymm1\n\t"\
		"vaddpd	%%ymm3,%%ymm4,%%ymm4\n\t"\
		"vaddpd	%%ymm6,%%ymm7,%%ymm7\n\t"\
		"vmovaps	%%ymm5,0x020(%%rax)\n\t"\
		"vmovaps	%%ymm1,0x020(%%rbx)\n\t"\
		"vmovaps	%%ymm4,     (%%rax)\n\t"\
		"vmovaps	%%ymm7,     (%%rbx)\n\t"\
		"addq	%%rdi,%%rax\n\t"\
		"addq	%%rdi,%%rbx\n\t"\
		"addq	$0x80,%%rsi	/* c0C */\n\t"\
		"vmovaps	     (%%rcx),%%ymm4\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm0\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5\n\t"\
		"vmovaps	     (%%rdx),%%ymm6\n\t"\
		"vmovaps	%%ymm4,%%ymm2\n\t"\
		"vmovaps	%%ymm0,%%ymm1\n\t"\
		"vmovaps	%%ymm5,%%ymm3\n\t"\
		"vmovaps	%%ymm6,%%ymm7\n\t"\
		"vmulpd	     (%%rsi),%%ymm4,%%ymm4\n\t"\
		"vmulpd	0x040(%%rsi),%%ymm0,%%ymm0\n\t"\
		"vmulpd	     (%%rsi),%%ymm5,%%ymm5\n\t"\
		"vmulpd	0x040(%%rsi),%%ymm6,%%ymm6\n\t"\
		"vmulpd	0x020(%%rsi),%%ymm2,%%ymm2\n\t"\
		"vmulpd	0x060(%%rsi),%%ymm1,%%ymm1\n\t"\
		"vmulpd	0x020(%%rsi),%%ymm3,%%ymm3\n\t"\
		"vmulpd	0x060(%%rsi),%%ymm7,%%ymm7\n\t"\
		"vsubpd	%%ymm2,%%ymm5,%%ymm5\n\t"\
		"vsubpd	%%ymm1,%%ymm6,%%ymm6\n\t"\
		"vaddpd	%%ymm3,%%ymm4,%%ymm4\n\t"\
		"vaddpd	%%ymm7,%%ymm0,%%ymm0\n\t"\
		"vmovaps	%%ymm5,0x020(%%rax)\n\t"\
		"vmovaps	%%ymm6,0x020(%%rbx)\n\t"\
		"vmovaps	%%ymm4,     (%%rax)\n\t"\
		"vmovaps	%%ymm0,     (%%rbx)\n\t"\
		"\n\t"\
		"/*...Block 6: t0A,t1A,t2A,t3A	*/\n\t"\
		"movq	%[__add0],%%rax\n\t"\
		"movslq	%[__p05],%%rsi\n\t"\
		"movslq	%[__p08],%%rbx\n\t"\
		"shlq	$3,%%rsi\n\t"\
		"shlq	$3,%%rbx\n\t"\
		"addq	%%rsi	,%%rax	/* add0 = &a[j1+p5] */\n\t"\
		"addq	%%rax	,%%rbx	/* add1 = add0+p8 */\n\t"\
		"movq	%[__r0A],%%rcx\n\t"\
		"movq	%[__r1A],%%rdx\n\t"\
		"movq	%[__isrt2],%%rsi\n\t"\
		"addq	$0x400,%%rcx\n\t"\
		"addq	$0x400,%%rdx\n\t"\
		"addq	$0x0a0,%%rsi	/* cc3 */\n\t"\
		"\n\t"\
		"vmovaps	     (%%rcx),%%ymm4	/* t2A */\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5	/* t2B */\n\t"\
		"vmovaps	     (%%rsi),%%ymm3	/* c32_3 */\n\t"\
		"vmovaps	0x020(%%rsi),%%ymm2	/* s32_3 */\n\t"\
		"vmovaps	%%ymm4,%%ymm6	/* ymm2 <- cpy t2A */\n\t"\
		"vmovaps	%%ymm5,%%ymm7	/* ymm3 <- cpy t2B */\n\t"\
		"\n\t"\
		"vmulpd	%%ymm2,%%ymm4,%%ymm4	/* t2A*s32_3 */\n\t"\
		"vmulpd	%%ymm2,%%ymm5,%%ymm5	/* t2B*s32_3 */\n\t"\
		"vmulpd	%%ymm3,%%ymm6,%%ymm6	/* t2A*c32_3 */\n\t"\
		"vmovaps	     (%%rdx),%%ymm0	/* t3A */\n\t"\
		"vmulpd	%%ymm3,%%ymm7,%%ymm7	/* t2B*c32_3 */\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm1	/* t3B */\n\t"\
		"subq	$0x40,%%rsi	/* cc1 */\n\t"\
		"vmovaps	     (%%rsi),%%ymm2	/* c32_1 */\n\t"\
		"vmovaps	0x020(%%rsi),%%ymm3	/* s32_1 */\n\t"\
		"vsubpd	%%ymm6,%%ymm5,%%ymm5	/* ymm5 <-~t2B */\n\t"\
		"vmovaps	%%ymm0,%%ymm6	/* ymm6 <- cpy t3A */\n\t"\
		"vaddpd	%%ymm7,%%ymm4,%%ymm4	/* ymm4 <-~t2A */\n\t"\
		"vmovaps	%%ymm1,%%ymm7	/* ymm7 <- cpy t3B */\n\t"\
		"\n\t"\
		"vmulpd	%%ymm2,%%ymm0,%%ymm0	/* t3A*c32_1 */\n\t"\
		"vmulpd	%%ymm2,%%ymm1,%%ymm1	/* t3B*c32_1 */\n\t"\
		"vmulpd	%%ymm3,%%ymm6,%%ymm6	/* t3A*s32_1 */\n\t"\
		"vmulpd	%%ymm3,%%ymm7,%%ymm7	/* t3B*s32_1 */\n\t"\
		"vaddpd	%%ymm6,%%ymm1,%%ymm1	/* ymm1 <- it */\n\t"\
		"vsubpd	%%ymm7,%%ymm0,%%ymm0	/* ymm0 <- rt */\n\t"\
		"\n\t"\
		"vmovaps	%%ymm5,%%ymm7	/* ymm7 <- cpy~t2B*/\n\t"\
		"vmovaps	%%ymm4,%%ymm6	/* ymm6 <- cpy~t2A*/\n\t"\
		"\n\t"\
		"vaddpd	%%ymm0,%%ymm6,%%ymm6	/* ~t3A <- t2A+rt */\n\t"\
		"vaddpd	%%ymm1,%%ymm7,%%ymm7	/* ~t3B <- t2B+it */\n\t"\
		"vsubpd	%%ymm0,%%ymm4,%%ymm4	/* ~t2A <- t2A-rt */\n\t"\
		"vsubpd	%%ymm1,%%ymm5,%%ymm5	/* ~t2B <- t2B-it */\n\t"\
		"\n\t"\
		"subq	$0x400,%%rcx\n\t"\
		"subq	$0x400,%%rdx\n\t"\
		"subq	$0x40,%%rsi	/* cc0 */\n\t"\
		"vmovaps	     (%%rdx),%%ymm0	/* t1A */\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm2	/* t1B */\n\t"\
		"vmovaps	0x020(%%rsi),%%ymm1	/* s */\n\t"\
		"vmovaps	%%ymm0,%%ymm3	/* cpy t1A */\n\t"\
		"vmulpd	%%ymm1,%%ymm0,%%ymm0	/* t1A*s */\n\t"\
		"vmulpd	%%ymm2,%%ymm1,%%ymm1	/* s*t1B */\n\t"\
		"vmulpd	(%%rsi),%%ymm3,%%ymm3	/* t1A*c */\n\t"\
		"vmulpd	(%%rsi),%%ymm2,%%ymm2	/* t1B*c */\n\t"\
		"\n\t"\
		"vsubpd	%%ymm2,%%ymm0,%%ymm0	/* rt =t1A*s - t1B*c */\n\t"\
		"vaddpd	%%ymm3,%%ymm1,%%ymm1	/* it =t1B*s + t1A*c */\n\t"\
		"\n\t"\
		"vmovaps	     (%%rcx),%%ymm2	/* t0A */\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm3	/* t0B */\n\t"\
		"vsubpd	%%ymm0,%%ymm2,%%ymm2	/*~t0A <- t0A- rt */\n\t"\
		"vsubpd	%%ymm1,%%ymm3,%%ymm3	/*~t0B <- t0B- it */\n\t"\
		"vaddpd	%%ymm0,%%ymm0,%%ymm0	/*          2* rt */\n\t"\
		"vaddpd	%%ymm1,%%ymm1,%%ymm1	/*          2* it */\n\t"\
		"vaddpd	%%ymm2,%%ymm0,%%ymm0	/*~t1A <- t0A+ rt */\n\t"\
		"vaddpd	%%ymm3,%%ymm1,%%ymm1	/*~t1B <- t0B+ it */\n\t"\
		"\n\t"\
		"/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c05): */\n\t"\
		"addq	$0x5c0,%%rsi	/* c05 */\n\t"\
		"vsubpd	%%ymm4,%%ymm2,%%ymm2\n\t"\
		"vsubpd	%%ymm7,%%ymm0,%%ymm0\n\t"\
		"vsubpd	%%ymm5,%%ymm3,%%ymm3\n\t"\
		"vsubpd	%%ymm6,%%ymm1,%%ymm1\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6\n\t"\
		"vaddpd	%%ymm2,%%ymm4,%%ymm4\n\t"\
		"vaddpd	%%ymm0,%%ymm7,%%ymm7\n\t"\
		"vaddpd	%%ymm3,%%ymm5,%%ymm5\n\t"\
		"vaddpd	%%ymm1,%%ymm6,%%ymm6\n\t"\
		"vmovaps	%%ymm2,     (%%rcx)\n\t"\
		"vmovaps	%%ymm0,0x020(%%rdx)\n\t"\
		"vmovaps	%%ymm3,0x020(%%rcx)\n\t"\
		"vmovaps	%%ymm6,     (%%rdx)\n\t"\
		"vmovaps	%%ymm4,%%ymm2\n\t"\
		"vmovaps	%%ymm7,%%ymm0\n\t"\
		"vmovaps	%%ymm5,%%ymm3\n\t"\
		"vmovaps	%%ymm1,%%ymm6\n\t"\
		"vmulpd	     (%%rsi),%%ymm4,%%ymm4\n\t"\
		"vmulpd	0x040(%%rsi),%%ymm7,%%ymm7\n\t"\
		"vmulpd	     (%%rsi),%%ymm5,%%ymm5\n\t"\
		"vmulpd	0x040(%%rsi),%%ymm1,%%ymm1\n\t"\
		"vmulpd	0x020(%%rsi),%%ymm2,%%ymm2\n\t"\
		"vmulpd	0x060(%%rsi),%%ymm0,%%ymm0\n\t"\
		"vmulpd	0x020(%%rsi),%%ymm3,%%ymm3\n\t"\
		"vmulpd	0x060(%%rsi),%%ymm6,%%ymm6\n\t"\
		"vsubpd	%%ymm2,%%ymm5,%%ymm5\n\t"\
		"vsubpd	%%ymm0,%%ymm1,%%ymm1\n\t"\
		"vaddpd	%%ymm3,%%ymm4,%%ymm4\n\t"\
		"vaddpd	%%ymm6,%%ymm7,%%ymm7\n\t"\
		"vmovaps	%%ymm5,0x020(%%rax)\n\t"\
		"vmovaps	%%ymm1,0x020(%%rbx)\n\t"\
		"vmovaps	%%ymm4,     (%%rax)\n\t"\
		"vmovaps	%%ymm7,     (%%rbx)\n\t"\
		"addq	%%rdi,%%rax\n\t"\
		"addq	%%rdi,%%rbx\n\t"\
		"addq	$0x80,%%rsi	/* c0D */\n\t"\
		"vmovaps	     (%%rcx),%%ymm4\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm0\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5\n\t"\
		"vmovaps	     (%%rdx),%%ymm6\n\t"\
		"vmovaps	%%ymm4,%%ymm2\n\t"\
		"vmovaps	%%ymm0,%%ymm1\n\t"\
		"vmovaps	%%ymm5,%%ymm3\n\t"\
		"vmovaps	%%ymm6,%%ymm7\n\t"\
		"vmulpd	     (%%rsi),%%ymm4,%%ymm4\n\t"\
		"vmulpd	0x040(%%rsi),%%ymm0,%%ymm0\n\t"\
		"vmulpd	     (%%rsi),%%ymm5,%%ymm5\n\t"\
		"vmulpd	0x040(%%rsi),%%ymm6,%%ymm6\n\t"\
		"vmulpd	0x020(%%rsi),%%ymm2,%%ymm2\n\t"\
		"vmulpd	0x060(%%rsi),%%ymm1,%%ymm1\n\t"\
		"vmulpd	0x020(%%rsi),%%ymm3,%%ymm3\n\t"\
		"vmulpd	0x060(%%rsi),%%ymm7,%%ymm7\n\t"\
		"vsubpd	%%ymm2,%%ymm5,%%ymm5\n\t"\
		"vsubpd	%%ymm1,%%ymm6,%%ymm6\n\t"\
		"vaddpd	%%ymm3,%%ymm4,%%ymm4\n\t"\
		"vaddpd	%%ymm7,%%ymm0,%%ymm0\n\t"\
		"vmovaps	%%ymm5,0x020(%%rax)\n\t"\
		"vmovaps	%%ymm6,0x020(%%rbx)\n\t"\
		"vmovaps	%%ymm4,     (%%rax)\n\t"\
		"vmovaps	%%ymm0,     (%%rbx)\n\t"\
		"\n\t"\
		"/*...Block 7: t0C,t1C,t2C,t3C	*/\n\t"\
		"movq	%[__add0],%%rax\n\t"\
		"movslq	%[__p06],%%rsi\n\t"\
		"movslq	%[__p08],%%rbx\n\t"\
		"shlq	$3,%%rsi\n\t"\
		"shlq	$3,%%rbx\n\t"\
		"addq	%%rsi	,%%rax	/* add0 = &a[j1+p6] */\n\t"\
		"addq	%%rax	,%%rbx	/* add1 = add0+p8 */\n\t"\
		"movq	%[__r0C],%%rcx\n\t"\
		"movq	%[__r1C],%%rdx\n\t"\
		"movq	%[__isrt2],%%rsi\n\t"\
		"addq	$0x400,%%rcx\n\t"\
		"addq	$0x400,%%rdx\n\t"\
		"addq	$0x020,%%rsi	/* cc0 */\n\t"\
		"\n\t"\
		"vmovaps	     (%%rcx),%%ymm4	/* t2C */\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5	/* t2D */\n\t"\
		"vmovaps	     (%%rsi),%%ymm3	/* c */\n\t"\
		"vmovaps	0x020(%%rsi),%%ymm2	/* s */\n\t"\
		"vmovaps	%%ymm4,%%ymm6	/* ymm2 <- cpy t2C */\n\t"\
		"vmovaps	%%ymm5,%%ymm7	/* ymm3 <- cpy t2D */\n\t"\
		"\n\t"\
		"vmulpd	%%ymm2,%%ymm4,%%ymm4	/* t2C*s */\n\t"\
		"vmulpd	%%ymm2,%%ymm5,%%ymm5	/* t2D*s */\n\t"\
		"vmulpd	%%ymm3,%%ymm6,%%ymm6	/* t2C*c */\n\t"\
		"vmovaps	     (%%rdx),%%ymm0	/* t3C */\n\t"\
		"vmulpd	%%ymm3,%%ymm7,%%ymm7	/* t2D*c */\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm1	/* t3D */\n\t"\
		"vsubpd	%%ymm6,%%ymm5,%%ymm5	/* ymm5 <-~t2D */\n\t"\
		"vmovaps	%%ymm0,%%ymm6	/* ymm6 <- cpy t3C */\n\t"\
		"vaddpd	%%ymm7,%%ymm4,%%ymm4	/* ymm4 <-~t2C */\n\t"\
		"vmovaps	%%ymm1,%%ymm7	/* ymm7 <- cpy t3D */\n\t"\
		"\n\t"\
		"vmulpd	%%ymm3,%%ymm0,%%ymm0	/* t3C*c */\n\t"\
		"vmulpd	%%ymm3,%%ymm1,%%ymm1	/* t3D*c */\n\t"\
		"vmulpd	%%ymm2,%%ymm6,%%ymm6	/* t3C*s */\n\t"\
		"vmulpd	%%ymm2,%%ymm7,%%ymm7	/* t3D*s */\n\t"\
		"vsubpd	%%ymm6,%%ymm1,%%ymm1	/* ymm1 <- it */\n\t"\
		"vaddpd	%%ymm7,%%ymm0,%%ymm0	/* ymm0 <- rt */\n\t"\
		"\n\t"\
		"vmovaps	%%ymm5,%%ymm7	/* ymm7 <- cpy~t2D*/\n\t"\
		"vmovaps	%%ymm4,%%ymm6	/* ymm6 <- cpy~t2C*/\n\t"\
		"\n\t"\
		"vaddpd	%%ymm0,%%ymm6,%%ymm6	/* ~t3C <- t2C+rt */\n\t"\
		"vaddpd	%%ymm1,%%ymm7,%%ymm7	/* ~t3D <- t2D+it */\n\t"\
		"vsubpd	%%ymm0,%%ymm4,%%ymm4	/* ~t2C <- t2C-rt */\n\t"\
		"vsubpd	%%ymm1,%%ymm5,%%ymm5	/* ~t2D <- t2D-it */\n\t"\
		"\n\t"\
		"subq	$0x400,%%rcx\n\t"\
		"subq	$0x400,%%rdx\n\t"\
		"subq	$0x20,%%rsi	/* isrt2 */\n\t"\
		"vmovaps	     (%%rdx),%%ymm0	/* t1C */\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm1	/* t1D */\n\t"\
		"vmovaps	(%%rsi),%%ymm3	/* isrt2 */\n\t"\
		"vmovaps	%%ymm0,%%ymm2	/* cpy t1C */\n\t"\
		"vsubpd	%%ymm1,%%ymm0,%%ymm0	/*~t1C=t1C-t1D */\n\t"\
		"vaddpd	%%ymm2,%%ymm1,%%ymm1	/*~t1D=t1D+t1C */\n\t"\
		"vmulpd	%%ymm3,%%ymm0,%%ymm0	/* it */\n\t"\
		"vmulpd	%%ymm3,%%ymm1,%%ymm1	/* rt */\n\t"\
		"\n\t"\
		"vmovaps	     (%%rcx),%%ymm2	/* t0C */\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm3	/* t0D */\n\t"\
		"vsubpd	%%ymm0,%%ymm2,%%ymm2	/*~t0C <- t0C- rt */\n\t"\
		"vsubpd	%%ymm1,%%ymm3,%%ymm3	/*~t0D <- t0D- it */\n\t"\
		"vaddpd	%%ymm0,%%ymm0,%%ymm0	/*          2* rt */\n\t"\
		"vaddpd	%%ymm1,%%ymm1,%%ymm1	/*          2* it */\n\t"\
		"vaddpd	%%ymm2,%%ymm0,%%ymm0	/*~t1C <- t0C+ rt */\n\t"\
		"vaddpd	%%ymm3,%%ymm1,%%ymm1	/*~t1D <- t0D+ it */\n\t"\
		"\n\t"\
		"/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c06): */\n\t"\
		"addq	$0x3e0,%%rsi	/* c06 */\n\t"\
		"vsubpd	%%ymm4,%%ymm2,%%ymm2\n\t"\
		"vsubpd	%%ymm7,%%ymm0,%%ymm0\n\t"\
		"vsubpd	%%ymm5,%%ymm3,%%ymm3\n\t"\
		"vsubpd	%%ymm6,%%ymm1,%%ymm1\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6\n\t"\
		"vaddpd	%%ymm2,%%ymm4,%%ymm4\n\t"\
		"vaddpd	%%ymm0,%%ymm7,%%ymm7\n\t"\
		"vaddpd	%%ymm3,%%ymm5,%%ymm5\n\t"\
		"vaddpd	%%ymm1,%%ymm6,%%ymm6\n\t"\
		"vmovaps	%%ymm2,     (%%rcx)\n\t"\
		"vmovaps	%%ymm0,0x020(%%rdx)\n\t"\
		"vmovaps	%%ymm3,0x020(%%rcx)\n\t"\
		"vmovaps	%%ymm6,     (%%rdx)\n\t"\
		"vmovaps	%%ymm4,%%ymm2\n\t"\
		"vmovaps	%%ymm7,%%ymm0\n\t"\
		"vmovaps	%%ymm5,%%ymm3\n\t"\
		"vmovaps	%%ymm1,%%ymm6\n\t"\
		"vmulpd	     (%%rsi),%%ymm4,%%ymm4\n\t"\
		"vmulpd	0x040(%%rsi),%%ymm7,%%ymm7\n\t"\
		"vmulpd	     (%%rsi),%%ymm5,%%ymm5\n\t"\
		"vmulpd	0x040(%%rsi),%%ymm1,%%ymm1\n\t"\
		"vmulpd	0x020(%%rsi),%%ymm2,%%ymm2\n\t"\
		"vmulpd	0x060(%%rsi),%%ymm0,%%ymm0\n\t"\
		"vmulpd	0x020(%%rsi),%%ymm3,%%ymm3\n\t"\
		"vmulpd	0x060(%%rsi),%%ymm6,%%ymm6\n\t"\
		"vsubpd	%%ymm2,%%ymm5,%%ymm5\n\t"\
		"vsubpd	%%ymm0,%%ymm1,%%ymm1\n\t"\
		"vaddpd	%%ymm3,%%ymm4,%%ymm4\n\t"\
		"vaddpd	%%ymm6,%%ymm7,%%ymm7\n\t"\
		"vmovaps	%%ymm5,0x020(%%rax)\n\t"\
		"vmovaps	%%ymm1,0x020(%%rbx)\n\t"\
		"vmovaps	%%ymm4,     (%%rax)\n\t"\
		"vmovaps	%%ymm7,     (%%rbx)\n\t"\
		"addq	%%rdi,%%rax\n\t"\
		"addq	%%rdi,%%rbx\n\t"\
		"addq	$0x80,%%rsi	/* c0E */\n\t"\
		"vmovaps	     (%%rcx),%%ymm4\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm0\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5\n\t"\
		"vmovaps	     (%%rdx),%%ymm6\n\t"\
		"vmovaps	%%ymm4,%%ymm2\n\t"\
		"vmovaps	%%ymm0,%%ymm1\n\t"\
		"vmovaps	%%ymm5,%%ymm3\n\t"\
		"vmovaps	%%ymm6,%%ymm7\n\t"\
		"vmulpd	     (%%rsi),%%ymm4,%%ymm4\n\t"\
		"vmulpd	0x040(%%rsi),%%ymm0,%%ymm0\n\t"\
		"vmulpd	     (%%rsi),%%ymm5,%%ymm5\n\t"\
		"vmulpd	0x040(%%rsi),%%ymm6,%%ymm6\n\t"\
		"vmulpd	0x020(%%rsi),%%ymm2,%%ymm2\n\t"\
		"vmulpd	0x060(%%rsi),%%ymm1,%%ymm1\n\t"\
		"vmulpd	0x020(%%rsi),%%ymm3,%%ymm3\n\t"\
		"vmulpd	0x060(%%rsi),%%ymm7,%%ymm7\n\t"\
		"vsubpd	%%ymm2,%%ymm5,%%ymm5\n\t"\
		"vsubpd	%%ymm1,%%ymm6,%%ymm6\n\t"\
		"vaddpd	%%ymm3,%%ymm4,%%ymm4\n\t"\
		"vaddpd	%%ymm7,%%ymm0,%%ymm0\n\t"\
		"vmovaps	%%ymm5,0x020(%%rax)\n\t"\
		"vmovaps	%%ymm6,0x020(%%rbx)\n\t"\
		"vmovaps	%%ymm4,     (%%rax)\n\t"\
		"vmovaps	%%ymm0,     (%%rbx)\n\t"\
		"\n\t"\
		"/*...Block 8: t0E,t1E,t2E,t3E	*/\n\t"\
		"movq	%[__add0],%%rax\n\t"\
		"movslq	%[__p07],%%rsi\n\t"\
		"movslq	%[__p08],%%rbx\n\t"\
		"shlq	$3,%%rsi\n\t"\
		"shlq	$3	,%%rbx	/* Pointer offset for floating doubles */\n\t"\
		"addq	%%rsi	,%%rax	/* add0 = &a[j1+p7] */\n\t"\
		"addq	%%rax	,%%rbx	/* add1 = add0+p8 */\n\t"\
		"movq	%[__r0E],%%rcx\n\t"\
		"movq	%[__r1E],%%rdx\n\t"\
		"movq	%[__isrt2],%%rsi\n\t"\
		"addq	$0x400,%%rcx\n\t"\
		"addq	$0x400,%%rdx\n\t"\
		"addq	$0x060,%%rsi	/* cc1 */\n\t"\
		"\n\t"\
		"vmovaps	     (%%rcx),%%ymm4	/* t2E */\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5	/* t2F */\n\t"\
		"vmovaps	     (%%rsi),%%ymm3	/* c32_1 */\n\t"\
		"vmovaps	0x020(%%rsi),%%ymm2	/* s32_1 */\n\t"\
		"vmovaps	%%ymm4,%%ymm6	/* ymm2 <- cpy t2E */\n\t"\
		"vmovaps	%%ymm5,%%ymm7	/* ymm3 <- cpy t2F */\n\t"\
		"\n\t"\
		"vmulpd	%%ymm2,%%ymm4,%%ymm4	/* t2E*s32_1 */\n\t"\
		"vmulpd	%%ymm2,%%ymm5,%%ymm5	/* t2F*s32_1 */\n\t"\
		"vmulpd	%%ymm3,%%ymm6,%%ymm6	/* t2E*c32_1 */\n\t"\
		"vmovaps	     (%%rdx),%%ymm0	/* t3E */\n\t"\
		"vmulpd	%%ymm3,%%ymm7,%%ymm7	/* t2F*c32_1 */\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm1	/* t3F */\n\t"\
		"addq	$0x40,%%rsi	/* cc3 */\n\t"\
		"vmovaps	     (%%rsi),%%ymm3	/* c32_3 */\n\t"\
		"vmovaps	0x020(%%rsi),%%ymm2	/* s32_3 */\n\t"\
		"vsubpd	%%ymm6,%%ymm5,%%ymm5	/* ymm5 <-~t2F */\n\t"\
		"vmovaps	%%ymm0,%%ymm6	/* ymm6 <- cpy t3E */\n\t"\
		"vaddpd	%%ymm7,%%ymm4,%%ymm4	/* ymm4 <-~t2E */\n\t"\
		"vmovaps	%%ymm1,%%ymm7	/* ymm7 <- cpy t3F */\n\t"\
		"\n\t"\
		"vmulpd	%%ymm2,%%ymm0,%%ymm0	/* t3E*s32_3 */\n\t"\
		"vmulpd	%%ymm2,%%ymm1,%%ymm1	/* t3F*s32_3 */\n\t"\
		"vmulpd	%%ymm3,%%ymm6,%%ymm6	/* t3E*c32_3 */\n\t"\
		"vmulpd	%%ymm3,%%ymm7,%%ymm7	/* t3F*c32_3 */\n\t"\
		"vsubpd	%%ymm6,%%ymm1,%%ymm1	/* ymm1 <- it */\n\t"\
		"vaddpd	%%ymm7,%%ymm0,%%ymm0	/* ymm0 <- rt */\n\t"\
		"\n\t"\
		"vmovaps	%%ymm5,%%ymm7	/* ymm7 <- cpy~t2F*/\n\t"\
		"vmovaps	%%ymm4,%%ymm6	/* ymm6 <- cpy~t2E*/\n\t"\
		"\n\t"\
		"vsubpd	%%ymm0,%%ymm4,%%ymm4	/* ~t2E <- t2E-rt */\n\t"\
		"vsubpd	%%ymm1,%%ymm5,%%ymm5	/* ~t2F <- t2F-it */\n\t"\
		"vaddpd	%%ymm0,%%ymm6,%%ymm6	/* ~t3E <- t2E+rt */\n\t"\
		"vaddpd	%%ymm1,%%ymm7,%%ymm7	/* ~t3F <- t2F+it */\n\t"\
		"\n\t"\
		"subq	$0x400,%%rcx\n\t"\
		"subq	$0x400,%%rdx\n\t"\
		"subq	$0x80,%%rsi	 /* cc0 */\n\t"\
		"vmovaps	     (%%rdx),%%ymm3	/* t1E */\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm1	/* t1F */\n\t"\
		"vmovaps	0x020(%%rsi),%%ymm2	/* s */\n\t"\
		"vmovaps	%%ymm3,%%ymm0	/* cpy t1E */\n\t"\
		"vmulpd	%%ymm2,%%ymm3,%%ymm3	/* t1E*s */\n\t"\
		"vmulpd	%%ymm1,%%ymm2,%%ymm2	/* t1F*s */\n\t"\
		"vmulpd	(%%rsi),%%ymm0,%%ymm0	/* t1E*c */\n\t"\
		"vmulpd	(%%rsi),%%ymm1,%%ymm1	/* t1F*c */\n\t"\
		"\n\t"\
		"vsubpd	%%ymm2,%%ymm0,%%ymm0	/* rt =t1E*c - t1F*s */\n\t"\
		"vaddpd	%%ymm3,%%ymm1,%%ymm1	/* it =t1F*c + t1E*s */\n\t"\
		"\n\t"\
		"vmovaps	     (%%rcx),%%ymm2	/* t0E */\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm3	/* t0F */\n\t"\
		"vsubpd	%%ymm0,%%ymm2,%%ymm2	/*~t0E <- t0E- rt */\n\t"\
		"vsubpd	%%ymm1,%%ymm3,%%ymm3	/*~t0F <- t0F- it */\n\t"\
		"vaddpd	%%ymm0,%%ymm0,%%ymm0	/*          2* rt */\n\t"\
		"vaddpd	%%ymm1,%%ymm1,%%ymm1	/*          2* it */\n\t"\
		"vaddpd	%%ymm2,%%ymm0,%%ymm0	/*~t1E <- t0E+ rt */\n\t"\
		"vaddpd	%%ymm3,%%ymm1,%%ymm1	/*~t1F <- t0F+ it */\n\t"\
		"/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c07): */\n\t"\
		"addq	$0x7c0,%%rsi	/* c07 */\n\t"\
		"vsubpd	%%ymm4,%%ymm2,%%ymm2\n\t"\
		"vsubpd	%%ymm7,%%ymm0,%%ymm0\n\t"\
		"vsubpd	%%ymm5,%%ymm3,%%ymm3\n\t"\
		"vsubpd	%%ymm6,%%ymm1,%%ymm1\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6\n\t"\
		"vaddpd	%%ymm2,%%ymm4,%%ymm4\n\t"\
		"vaddpd	%%ymm0,%%ymm7,%%ymm7\n\t"\
		"vaddpd	%%ymm3,%%ymm5,%%ymm5\n\t"\
		"vaddpd	%%ymm1,%%ymm6,%%ymm6\n\t"\
		"vmovaps	%%ymm2,     (%%rcx)\n\t"\
		"vmovaps	%%ymm0,0x020(%%rdx)\n\t"\
		"vmovaps	%%ymm3,0x020(%%rcx)\n\t"\
		"vmovaps	%%ymm6,     (%%rdx)\n\t"\
		"vmovaps	%%ymm4,%%ymm2\n\t"\
		"vmovaps	%%ymm7,%%ymm0\n\t"\
		"vmovaps	%%ymm5,%%ymm3\n\t"\
		"vmovaps	%%ymm1,%%ymm6\n\t"\
		"vmulpd	     (%%rsi),%%ymm4,%%ymm4\n\t"\
		"vmulpd	0x040(%%rsi),%%ymm7,%%ymm7\n\t"\
		"vmulpd	     (%%rsi),%%ymm5,%%ymm5\n\t"\
		"vmulpd	0x040(%%rsi),%%ymm1,%%ymm1\n\t"\
		"vmulpd	0x020(%%rsi),%%ymm2,%%ymm2\n\t"\
		"vmulpd	0x060(%%rsi),%%ymm0,%%ymm0\n\t"\
		"vmulpd	0x020(%%rsi),%%ymm3,%%ymm3\n\t"\
		"vmulpd	0x060(%%rsi),%%ymm6,%%ymm6\n\t"\
		"vsubpd	%%ymm2,%%ymm5,%%ymm5\n\t"\
		"vsubpd	%%ymm0,%%ymm1,%%ymm1\n\t"\
		"vaddpd	%%ymm3,%%ymm4,%%ymm4\n\t"\
		"vaddpd	%%ymm6,%%ymm7,%%ymm7\n\t"\
		"vmovaps	%%ymm5,0x020(%%rax)\n\t"\
		"vmovaps	%%ymm1,0x020(%%rbx)\n\t"\
		"vmovaps	%%ymm4,     (%%rax)\n\t"\
		"vmovaps	%%ymm7,     (%%rbx)\n\t"\
		"addq	%%rdi,%%rax\n\t"\
		"addq	%%rdi,%%rbx\n\t"\
		"addq	$0x80,%%rsi	/* c0F */\n\t"\
		"vmovaps	     (%%rcx),%%ymm4\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm0\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5\n\t"\
		"vmovaps	     (%%rdx),%%ymm6\n\t"\
		"vmovaps	%%ymm4,%%ymm2\n\t"\
		"vmovaps	%%ymm0,%%ymm1\n\t"\
		"vmovaps	%%ymm5,%%ymm3\n\t"\
		"vmovaps	%%ymm6,%%ymm7\n\t"\
		"vmulpd	     (%%rsi),%%ymm4,%%ymm4\n\t"\
		"vmulpd	0x040(%%rsi),%%ymm0,%%ymm0\n\t"\
		"vmulpd	     (%%rsi),%%ymm5,%%ymm5\n\t"\
		"vmulpd	0x040(%%rsi),%%ymm6,%%ymm6\n\t"\
		"vmulpd	0x020(%%rsi),%%ymm2,%%ymm2\n\t"\
		"vmulpd	0x060(%%rsi),%%ymm1,%%ymm1\n\t"\
		"vmulpd	0x020(%%rsi),%%ymm3,%%ymm3\n\t"\
		"vmulpd	0x060(%%rsi),%%ymm7,%%ymm7\n\t"\
		"vsubpd	%%ymm2,%%ymm5,%%ymm5\n\t"\
		"vsubpd	%%ymm1,%%ymm6,%%ymm6\n\t"\
		"vaddpd	%%ymm3,%%ymm4,%%ymm4\n\t"\
		"vaddpd	%%ymm7,%%ymm0,%%ymm0\n\t"\
		"vmovaps	%%ymm5,0x020(%%rax)\n\t"\
		"vmovaps	%%ymm6,0x020(%%rbx)\n\t"\
		"vmovaps	%%ymm4,     (%%rax)\n\t"\
		"vmovaps	%%ymm0,     (%%rbx)\n\t"\
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
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"	/* Clobbered registers */\
	);\
	}

  #else // USE_64BIT_ASM_STYLE = True: Deeper 64-bit-ified version of the above 32-bit ASM macros, using all of ymm0-15

	/*
	For GCC-macro version of this, use that isrt2 + 0x020,0x060,0x0a0 = cc0,cc1,cc3,
	and isrt2 + 0x0e0,0x1e0,0x2e0,0x3e0,0x4e0,0x5e0,0x6e0,0x7e0 = c00,04,02,06,01,05,03,07
	in order to reduce number of args to <= the GCC-allowed maximum of 30:
	*/
	#define SSE2_RADIX32_DIF_TWIDDLE(Xadd0,Xp01,Xp02,Xp03,Xp04,Xp08,Xp0C,Xp10,Xp18,Xr00)\
	{\
	__asm__ volatile (\
		"/*...Block 1: */								\n\t"\
		"movq	%[__add0],%%rax							\n\t"\
		"movslq	%[__p08],%%rbx							\n\t"\
		"movslq	%[__p10],%%rcx							\n\t	movslq	%[__p04],%%r9	\n\t"\
		"movslq	%[__p18],%%rdx							\n\t	movq	%[__r00],%%rsi	\n\t"\
		"shlq	$3,%%rbx								\n\t	shlq	$3,%%r9	\n\t"\
		"shlq	$3,%%rcx								\n\t	movq	%%rsi,%%r8	\n\t"\
		"shlq	$3,%%rdx								\n\t	movq	%%r9,%%r10	\n\t	addq	%%rax,%%r10	\n\t"\
		"addq	%%rax,%%rbx								\n\t	movq	%%r9,%%r11	\n\t	addq	%%rbx,%%r11	\n\t"\
		"addq	%%rax,%%rcx								\n\t	movq	%%r9,%%r12	\n\t	addq	%%rcx,%%r12	\n\t"\
		"addq	%%rax,%%rdx								\n\t	movq	%%r9,%%r13	\n\t	addq	%%rdx,%%r13	\n\t"\
		"/* SSE2_RADIX4_DIF_4TWIDDLE_A(r00,c10) */		\n\t	/* SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO_B(r08,c04) */\n\t"\
		"addq	$0x1100,%%r8	/* two */				\n\t	vmovaps	     (%%r10),%%ymm8 	\n\t"\
		"vmovaps	     (%%rax),%%ymm0					\n\t	vmovaps	     (%%r12),%%ymm12	\n\t"\
		"vmovaps	     (%%rcx),%%ymm4					\n\t	vmovaps	0x020(%%r10),%%ymm9 	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1					\n\t	vmovaps	0x020(%%r12),%%ymm13	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5					\n\t	vmovaps	0x9e0(%%rsi),%%ymm14	/* c04 */	\n\t"\
		"vmovaps	0x920(%%rsi),%%ymm2	/* c10 */		\n\t	vmovaps	0xa00(%%rsi),%%ymm15	\n\t"\
		"vmovaps	0x940(%%rsi),%%ymm3					\n\t	vmovaps	%%ymm8 ,%%ymm10	\n\t"\
		"vmovaps	%%ymm4,%%ymm6						\n\t	vmovaps	%%ymm9 ,%%ymm11	\n\t"\
		"vmovaps	%%ymm5,%%ymm7						\n\t	vmulpd	%%ymm14,%%ymm8 ,%%ymm8 	\n\t"\
		"vmulpd	%%ymm2,%%ymm4,%%ymm4					\n\t	vmulpd	%%ymm14,%%ymm9 ,%%ymm9 	\n\t"\
		"vmulpd	%%ymm2,%%ymm5,%%ymm5					\n\t	vmulpd	%%ymm15,%%ymm10,%%ymm10	\n\t"\
		"vmulpd	%%ymm3,%%ymm6,%%ymm6					\n\t	vmulpd	%%ymm15,%%ymm11,%%ymm11	\n\t"\
		"vmulpd	%%ymm3,%%ymm7,%%ymm7					\n\t	vaddpd	%%ymm10,%%ymm9 ,%%ymm9 	\n\t"\
		"vmovaps	%%ymm0,%%ymm2						\n\t	vmovaps	%%ymm12,%%ymm14	\n\t"\
		"vmovaps	%%ymm1,%%ymm3						\n\t	vmulpd	0xa20(%%rsi),%%ymm12,%%ymm12	/* c14 */	\n\t"\
		"vaddpd	%%ymm6,%%ymm5,%%ymm5					\n\t	vmovaps	%%ymm13,%%ymm15	\n\t"\
		"vsubpd	%%ymm7,%%ymm4,%%ymm4					\n\t	vsubpd	%%ymm11,%%ymm8 ,%%ymm8 	\n\t"\
		"vaddpd	%%ymm4,%%ymm0,%%ymm0					\n\t	vmulpd	0xa20(%%rsi),%%ymm13,%%ymm13	\n\t"\
		"vaddpd	%%ymm5,%%ymm1,%%ymm1					\n\t	vmulpd	0xa40(%%rsi),%%ymm14,%%ymm14	\n\t"\
		"vsubpd	%%ymm4,%%ymm2,%%ymm2					\n\t	vmovaps	%%ymm8 ,%%ymm10	\n\t"\
		"vsubpd	%%ymm5,%%ymm3,%%ymm3					\n\t	vmulpd	0xa40(%%rsi),%%ymm15,%%ymm15	\n\t"\
		"vmovaps	     (%%rdx),%%ymm4					\n\t	vaddpd	%%ymm14,%%ymm13,%%ymm13	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm5					\n\t	vmovaps	%%ymm9 ,%%ymm11	\n\t"\
		"vmovaps	     (%%rdx),%%ymm6					\n\t	vsubpd	%%ymm15,%%ymm12,%%ymm12	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm7					\n\t	vaddpd	%%ymm12,%%ymm8 ,%%ymm8 	\n\t"\
		"vmulpd	0x9a0(%%rsi),%%ymm4,%%ymm4	/* c18 */	\n\t	vaddpd	%%ymm13,%%ymm9 ,%%ymm9 	\n\t"\
		"vmulpd	0x9a0(%%rsi),%%ymm5,%%ymm5				\n\t	vsubpd	%%ymm12,%%ymm10,%%ymm10	\n\t"\
		"vmulpd	0x9c0(%%rsi),%%ymm6,%%ymm6				\n\t	vsubpd	%%ymm13,%%ymm11,%%ymm11	\n\t"\
		"vmulpd	0x9c0(%%rsi),%%ymm7,%%ymm7				\n\t	vmovaps	     (%%r13),%%ymm12	\n\t"\
		"vaddpd	%%ymm6,%%ymm5,%%ymm5					\n\t	vmovaps	0x020(%%r13),%%ymm13	\n\t"\
		"vsubpd	%%ymm7,%%ymm4,%%ymm4					\n\t	vmovaps	     (%%r13),%%ymm14	\n\t"\
		"vmovaps	%%ymm5,0x020(%%rsi)					\n\t	vmovaps	0x020(%%r13),%%ymm15	\n\t"\
		"vmovaps	%%ymm4,     (%%rsi)	/* tmpstr r00 */\n\t	vmulpd	0xaa0(%%rsi),%%ymm12,%%ymm12	/* c1C */	\n\t"\
		"vmovaps	     (%%rbx),%%ymm4					\n\t	vmulpd	0xaa0(%%rsi),%%ymm13,%%ymm13	\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm5					\n\t	vmulpd	0xac0(%%rsi),%%ymm14,%%ymm14	\n\t"\
		"vmovaps	     (%%rbx),%%ymm6					\n\t	vmulpd	0xac0(%%rsi),%%ymm15,%%ymm15	\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm7					\n\t	vaddpd	     %%ymm14,%%ymm13,%%ymm13	\n\t"\
		"vmulpd	0x960(%%rsi),%%ymm4,%%ymm4	/* c08 */	\n\t	vsubpd	     %%ymm15,%%ymm12,%%ymm12	\n\t"\
		"vmulpd	0x960(%%rsi),%%ymm5,%%ymm5				\n\t	vmovaps	%%ymm13,0x120(%%rsi)	\n\t"\
		"vmulpd	0x980(%%rsi),%%ymm6,%%ymm6				\n\t	vmovaps	%%ymm12,0x100(%%rsi)	/* tmpstr r08 */\n\t"\
		"vmulpd	0x980(%%rsi),%%ymm7,%%ymm7				\n\t	vmovaps	     (%%r11),%%ymm12	\n\t"\
		"vaddpd	%%ymm6,%%ymm5,%%ymm5					\n\t	vmovaps	0x020(%%r11),%%ymm13	\n\t"\
		"vsubpd	%%ymm7,%%ymm4,%%ymm4					\n\t	vmovaps	     (%%r11),%%ymm14	\n\t"\
		"vmovaps	%%ymm5,%%ymm7						\n\t	vmovaps	0x020(%%r11),%%ymm15	\n\t"\
		"vmovaps	%%ymm4,%%ymm6						\n\t	vmulpd	0xa60(%%rsi),%%ymm12,%%ymm12	/* c0C */	\n\t"\
		"vsubpd	     (%%rsi),%%ymm4,%%ymm4	/* r00 */	\n\t	vmulpd	0xa60(%%rsi),%%ymm13,%%ymm13	\n\t"\
		"vsubpd	0x020(%%rsi),%%ymm5,%%ymm5				\n\t	vmulpd	0xa80(%%rsi),%%ymm14,%%ymm14	\n\t"\
		"vaddpd	     (%%rsi),%%ymm6,%%ymm6				\n\t	vmulpd	0xa80(%%rsi),%%ymm15,%%ymm15	\n\t"\
		"vaddpd	0x020(%%rsi),%%ymm7,%%ymm7				\n\t	vaddpd	     %%ymm14,%%ymm13,%%ymm13	\n\t"\
		"vsubpd	%%ymm6,%%ymm0,%%ymm0					\n\t	vsubpd	     %%ymm15,%%ymm12,%%ymm12	\n\t"\
		"vsubpd	%%ymm5,%%ymm2,%%ymm2					\n\t	vmovaps	%%ymm13,%%ymm15	\n\t"\
		"vsubpd	%%ymm7,%%ymm1,%%ymm1					\n\t	vmovaps	%%ymm12,%%ymm14	\n\t"\
		"vsubpd	%%ymm4,%%ymm3,%%ymm3					\n\t	vsubpd	0x100(%%rsi),%%ymm12,%%ymm12	\n\t"\
		"/*vmovaps	%%ymm0,0x080(%%rsi)	*/				\n\t	vsubpd	0x120(%%rsi),%%ymm13,%%ymm13	\n\t"\
		"/*vmovaps	%%ymm2,0x040(%%rsi)	*/				\n\t	vaddpd	0x100(%%rsi),%%ymm14,%%ymm14	\n\t"\
		"/*vmovaps	%%ymm1,0x0a0(%%rsi)	*/				\n\t	vaddpd	0x120(%%rsi),%%ymm15,%%ymm15	\n\t"\
		"/*vmovaps	%%ymm3,0x0e0(%%rsi)	*/				\n\t	vsubpd	     %%ymm14,%%ymm8 ,%%ymm8 	\n\t"\
		"vmulpd	(%%r8),%%ymm6,%%ymm6					\n\t	vsubpd	     %%ymm15,%%ymm9 ,%%ymm9 	\n\t"\
		"vmulpd	(%%r8),%%ymm5,%%ymm5					\n\t	vmovaps	%%ymm8 ,0x180(%%rsi)	\n\t"\
		"vmulpd	(%%r8),%%ymm7,%%ymm7					\n\t	vsubpd	%%ymm13,%%ymm10,%%ymm10	\n\t"\
		"vmulpd	(%%r8),%%ymm4,%%ymm4					\n\t/*	vmovaps	%%ymm9 ,0x1a0(%%rsi)	*/\n\t"\
		"vaddpd	%%ymm0,%%ymm6,%%ymm6					\n\t	vsubpd	%%ymm12,%%ymm11,%%ymm11	\n\t"\
		"vaddpd	%%ymm2,%%ymm5,%%ymm5					\n\t	vmulpd	(%%r8) ,%%ymm14,%%ymm14	\n\t"\
		"vaddpd	%%ymm1,%%ymm7,%%ymm7					\n\t	vmulpd	(%%r8) ,%%ymm13,%%ymm13	\n\t"\
		"vaddpd	%%ymm3,%%ymm4,%%ymm4					\n\t	vmulpd	(%%r8) ,%%ymm15,%%ymm15	\n\t"\
		"/*vmovaps	%%ymm6,     (%%rsi)	*/				\n\t	vmulpd	(%%r8) ,%%ymm12,%%ymm12	\n\t"\
		"/*vmovaps	%%ymm5,0x0c0(%%rsi)	*/				\n\t	vaddpd	%%ymm8 ,%%ymm14,%%ymm14	\n\t"\
		"/*vmovaps	%%ymm7,0x020(%%rsi)	*/				\n\t	vaddpd	%%ymm10,%%ymm13,%%ymm13	\n\t"\
		"/*vmovaps	%%ymm4,0x060(%%rsi)	*/				\n\t	vaddpd	%%ymm9 ,%%ymm15,%%ymm15	\n\t"\
		"												\n\t	vaddpd	%%ymm11,%%ymm12,%%ymm12	\n\t"\
		"														vmovaps	%%ymm14,0x100(%%rsi)	\n\t"\
		"														vmovaps	%%ymm15,0x120(%%rsi)	\n\t"\
		"														vmovaps	0x800(%%rsi),%%ymm8 	/* isrt2 */	\n\t"\
		"														vmovaps	%%ymm10,%%ymm14	\n\t"\
		"														vmovaps	%%ymm13,%%ymm15	\n\t"\
		"														vsubpd	%%ymm12,%%ymm10,%%ymm10	\n\t"\
		"														vsubpd	%%ymm11,%%ymm13,%%ymm13	\n\t"\
		"														vaddpd	%%ymm12,%%ymm14,%%ymm14	\n\t"\
		"														vaddpd	%%ymm11,%%ymm15,%%ymm15	\n\t"\
		"														vmulpd	%%ymm8 ,%%ymm10,%%ymm10	\n\t"\
		"														vmulpd	%%ymm8 ,%%ymm13,%%ymm13	\n\t"\
		"														vmulpd	%%ymm8 ,%%ymm14,%%ymm14	\n\t"\
		"vsubpd	%%ymm10,%%ymm2,%%ymm2					\n\t	vmulpd	%%ymm8 ,%%ymm15,%%ymm15	\n\t"\
		"/* Since ymm10 output is first-done of 2nd set, move the computation of ymm2/10 outputs up so can use ymm2 for 2.0 doubling constant below */\n\t"\
		"													/*	vmovaps	%%ymm10,0x140(%%rsi)	*/\n\t"\
		"													/*	vmovaps	%%ymm13,0x1c0(%%rsi)	*/\n\t"\
		"													/*	vmovaps	%%ymm14,0x160(%%rsi)	*/\n\t"\
		"													/*	vmovaps	%%ymm15,0x1e0(%%rsi)	*/\n\t"\
		"/* SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r00) */	\n\t"\
		"vmovaps	0x100(%%rsi),%%ymm11		\n\t"\
		"vmovaps	0x120(%%rsi),%%ymm12		\n\t"\
		"						vmulpd	(%%r8),%%ymm10,%%ymm10	\n\t"\
		"vmovaps	0x180(%%rsi),%%ymm8 		\n\t"\
		"vsubpd	%%ymm11,%%ymm6,%%ymm6					\n\t/*	vsubpd	%%ymm10,%%ymm2,%%ymm2	*/\n\t"\
		"vsubpd	%%ymm9 ,%%ymm0,%%ymm0					\n\t	vsubpd	%%ymm15,%%ymm5,%%ymm5	\n\t"\
		"						vaddpd	%%ymm2 ,%%ymm10,%%ymm10	\n\t"\
		"						vmovaps	%%ymm2 ,0x140(%%rsi)	\n\t"\
		"vsubpd	%%ymm12,%%ymm7,%%ymm7					\n\t	vsubpd	%%ymm14,%%ymm4,%%ymm4	\n\t"\
		"vmovaps	(%%r8),%%ymm2	/* 2.0 ... We could use add-in-place for doubling, but want to load-balance add/mul here. */	\n\t"\
		"vsubpd	%%ymm8 ,%%ymm1 ,%%ymm1 					\n\t	vsubpd	%%ymm13,%%ymm3 ,%%ymm3 	\n\t"\
		"vmulpd	%%ymm2 ,%%ymm11,%%ymm11					\n\t/*	vaddpd	%%ymm10,%%ymm10,%%ymm10	*/\n\t"\
		"vmulpd	%%ymm2 ,%%ymm9 ,%%ymm9 					\n\t	vmulpd	%%ymm2 ,%%ymm15,%%ymm15	\n\t"\
		"vmulpd	%%ymm2 ,%%ymm12,%%ymm12					\n\t	vmulpd	%%ymm2 ,%%ymm14,%%ymm14	\n\t"\
		"vmulpd	%%ymm2 ,%%ymm8 ,%%ymm8 					\n\t	vmulpd	%%ymm2 ,%%ymm13,%%ymm13	\n\t"\
		"vaddpd	%%ymm6 ,%%ymm11,%%ymm11					\n\t/*	vaddpd	%%ymm2 ,%%ymm10,%%ymm10	*/\n\t"\
		"vaddpd	%%ymm0 ,%%ymm9 ,%%ymm9 					\n\t	vaddpd	%%ymm5 ,%%ymm15,%%ymm15	\n\t"\
		"vaddpd	%%ymm7 ,%%ymm12,%%ymm12					\n\t	vaddpd	%%ymm4 ,%%ymm14,%%ymm14	\n\t"\
		"vaddpd	%%ymm1 ,%%ymm8 ,%%ymm8 					\n\t	vaddpd	%%ymm3 ,%%ymm13,%%ymm13	\n\t"\
		"vmovaps	%%ymm6 ,0x100(%%rsi)				\n\t/*	vmovaps	%%ymm2 ,0x140(%%rsi)	*/\n\t"\
		"vmovaps	%%ymm0 ,0x080(%%rsi)				\n\t	vmovaps	%%ymm5 ,0x0c0(%%rsi)	\n\t"\
		"vmovaps	%%ymm7 ,0x120(%%rsi)				\n\t	vmovaps	%%ymm4 ,0x160(%%rsi)	\n\t"\
		"vmovaps	%%ymm1 ,0x1a0(%%rsi)				\n\t	vmovaps	%%ymm3 ,0x1e0(%%rsi)	\n\t"\
		"vmovaps	%%ymm11,     (%%rsi)				\n\t	vmovaps	%%ymm10,0x040(%%rsi)	\n\t"\
		"vmovaps	%%ymm9 ,0x180(%%rsi)				\n\t	vmovaps	%%ymm15,0x1c0(%%rsi)	\n\t"\
		"vmovaps	%%ymm12,0x020(%%rsi)				\n\t	vmovaps	%%ymm14,0x060(%%rsi)	\n\t"\
		"vmovaps	%%ymm8 ,0x0a0(%%rsi)				\n\t	vmovaps	%%ymm13,0x0e0(%%rsi)	\n\t"\
		"\n\t"\
		"/***************************************/\n\t"\
		"\n\t"\
		"/*...Block 2: */\n\t"\
		"movslq	%[__p02],%%rdi							\n\t	movslq	%[__p04],%%r9	\n\t"\
		"shlq	$3,%%rdi								\n\t	shlq	$3,%%r9	\n\t"\
		"addq	%%rdi,%%rax	/* &a[j1+p2] */				\n\t	movq	%%r9,%%r10	\n\t	addq	%%rax,%%r10	\n\t"\
		"addq	%%rdi,%%rbx								\n\t	movq	%%r9,%%r11	\n\t	addq	%%rbx,%%r11	\n\t"\
		"addq	%%rdi,%%rcx								\n\t	movq	%%r9,%%r12	\n\t	addq	%%rcx,%%r12	\n\t"\
		"addq	%%rdi,%%rdx								\n\t	movq	%%r9,%%r13	\n\t	addq	%%rdx,%%r13	\n\t"\
		"/* SSE2_RADIX4_DIF_4TWIDDLE_B(r10,c02) */		\n\t	/* SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO_B(r18,c06) */\n\t"\
		"addq	$0x200,%%rsi	/* r10 */				\n\t	vmovaps	     (%%r10),%%ymm8 	\n\t"\
		"vmovaps	     (%%rax),%%ymm0					\n\t	vmovaps	     (%%r12),%%ymm12	\n\t"\
		"vmovaps	     (%%rcx),%%ymm4					\n\t	vmovaps	0x020(%%r10),%%ymm9 	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1					\n\t	vmovaps	0x020(%%r12),%%ymm13	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5					\n\t	vmovaps	0x9e0(%%rsi),%%ymm14	/* c06 */	\n\t"\
		"vmovaps	0x8e0(%%rsi),%%ymm6	/* c02 */		\n\t	vmovaps	0xa00(%%rsi),%%ymm15	\n\t"\
		"vmovaps	0x900(%%rsi),%%ymm7					\n\t	vmovaps	%%ymm8 ,%%ymm10	\n\t"\
		"vmovaps	%%ymm0,%%ymm2						\n\t	vmovaps	%%ymm9 ,%%ymm11	\n\t"\
		"vmovaps	%%ymm1,%%ymm3						\n\t	vmulpd	%%ymm14,%%ymm8 ,%%ymm8 	\n\t"\
		"vmulpd	%%ymm6,%%ymm0,%%ymm0					\n\t	vmulpd	%%ymm14,%%ymm9 ,%%ymm9 	\n\t"\
		"vmulpd	%%ymm6,%%ymm1,%%ymm1					\n\t	vmulpd	%%ymm15,%%ymm10,%%ymm10	\n\t"\
		"vmulpd	%%ymm7,%%ymm2,%%ymm2					\n\t	vmulpd	%%ymm15,%%ymm11,%%ymm11	\n\t"\
		"vmulpd	%%ymm7,%%ymm3,%%ymm3					\n\t	vmovaps	%%ymm12,%%ymm14	\n\t"\
		"vmovaps	%%ymm4,%%ymm6						\n\t	vaddpd	%%ymm10,%%ymm9,%%ymm9 	\n\t"\
		"vaddpd	%%ymm2,%%ymm1,%%ymm1					\n\t	vmovaps	%%ymm13,%%ymm15	\n\t"\
		"vmovaps	%%ymm5,%%ymm7						\n\t	vmulpd	0xa20(%%rsi),%%ymm12,%%ymm12	/* c16 */	\n\t"\
		"vmulpd	0x920(%%rsi),%%ymm4,%%ymm4	/* c12 */	\n\t	vsubpd	%%ymm11,%%ymm8,%%ymm8 	\n\t"\
		"vsubpd	%%ymm3,%%ymm0,%%ymm0					\n\t	vmulpd	0xa20(%%rsi),%%ymm13,%%ymm13	\n\t"\
		"vmulpd	0x920(%%rsi),%%ymm5,%%ymm5				\n\t	vmulpd	0xa40(%%rsi),%%ymm14,%%ymm14	\n\t"\
		"vmulpd	0x940(%%rsi),%%ymm6,%%ymm6				\n\t	vmovaps	%%ymm8 ,%%ymm10	\n\t"\
		"vmovaps	%%ymm0,%%ymm2						\n\t	vmulpd	0xa40(%%rsi),%%ymm15,%%ymm15	\n\t"\
		"vmulpd	0x940(%%rsi),%%ymm7,%%ymm7				\n\t	vaddpd	%%ymm14,%%ymm13,%%ymm13	\n\t"\
		"vaddpd	%%ymm6,%%ymm5,%%ymm5					\n\t	vmovaps	%%ymm9 ,%%ymm11	\n\t"\
		"vmovaps	%%ymm1,%%ymm3						\n\t	vsubpd	%%ymm15,%%ymm12,%%ymm12	\n\t"\
		"vsubpd	%%ymm7,%%ymm4,%%ymm4					\n\t	vaddpd	%%ymm12,%%ymm8 ,%%ymm8 	\n\t"\
		"vaddpd	%%ymm4,%%ymm0,%%ymm0					\n\t	vaddpd	%%ymm13,%%ymm9 ,%%ymm9 	\n\t"\
		"vaddpd	%%ymm5,%%ymm1,%%ymm1					\n\t	vsubpd	%%ymm12,%%ymm10,%%ymm10	\n\t"\
		"vsubpd	%%ymm4,%%ymm2,%%ymm2					\n\t	vsubpd	%%ymm13,%%ymm11,%%ymm11	\n\t"\
		"vsubpd	%%ymm5,%%ymm3,%%ymm3					\n\t	vmovaps	     (%%r13),%%ymm12	\n\t"\
		"vmovaps	     (%%rdx),%%ymm4					\n\t	vmovaps	0x020(%%r13),%%ymm13	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm5					\n\t	vmovaps	     (%%r13),%%ymm14	\n\t"\
		"vmovaps	     (%%rdx),%%ymm6					\n\t	vmovaps	0x020(%%r13),%%ymm15	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm7					\n\t	vmulpd	0xaa0(%%rsi),%%ymm12,%%ymm12	/* c1E */	\n\t"\
		"vmulpd	0x9a0(%%rsi),%%ymm4,%%ymm4	/* c1A */	\n\t	vmulpd	0xaa0(%%rsi),%%ymm13,%%ymm13	\n\t"\
		"vmulpd	0x9a0(%%rsi),%%ymm5,%%ymm5				\n\t	vmulpd	0xac0(%%rsi),%%ymm14,%%ymm14	\n\t"\
		"vmulpd	0x9c0(%%rsi),%%ymm6,%%ymm6				\n\t	vmulpd	0xac0(%%rsi),%%ymm15,%%ymm15	\n\t"\
		"vmulpd	0x9c0(%%rsi),%%ymm7,%%ymm7				\n\t	vaddpd	     %%ymm14,%%ymm13,%%ymm13	\n\t"\
		"vaddpd	%%ymm6,%%ymm5,%%ymm5					\n\t	vsubpd	     %%ymm15,%%ymm12,%%ymm12	\n\t"\
		"vsubpd	%%ymm7,%%ymm4,%%ymm4					\n\t	vmovaps	%%ymm13,0x120(%%rsi)	\n\t"\
		"vmovaps	%%ymm5,0x020(%%rsi)					\n\t	vmovaps	%%ymm12,0x100(%%rsi)	/* tmpstr r08 */\n\t"\
		"vmovaps	%%ymm4,     (%%rsi)					\n\t	vmovaps	     (%%r11),%%ymm12	\n\t"\
		"vmovaps	     (%%rbx),%%ymm4					\n\t	vmovaps	0x020(%%r11),%%ymm13	\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm5					\n\t	vmovaps	     (%%r11),%%ymm14	\n\t"\
		"vmovaps	     (%%rbx),%%ymm6					\n\t	vmovaps	0x020(%%r11),%%ymm15	\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm7					\n\t	vmulpd	0xa60(%%rsi),%%ymm12,%%ymm12	/* c0E */	\n\t"\
		"vmulpd	0x960(%%rsi),%%ymm4,%%ymm4	/* c0A */	\n\t	vmulpd	0xa60(%%rsi),%%ymm13,%%ymm13	\n\t"\
		"vmulpd	0x960(%%rsi),%%ymm5,%%ymm5				\n\t	vmulpd	0xa80(%%rsi),%%ymm14,%%ymm14	\n\t"\
		"vmulpd	0x980(%%rsi),%%ymm6,%%ymm6				\n\t	vmulpd	0xa80(%%rsi),%%ymm15,%%ymm15	\n\t"\
		"vmulpd	0x980(%%rsi),%%ymm7,%%ymm7				\n\t	vaddpd	     %%ymm14,%%ymm13,%%ymm13	\n\t"\
		"vaddpd	%%ymm6,%%ymm5,%%ymm5					\n\t	vsubpd	     %%ymm15,%%ymm12,%%ymm12	\n\t"\
		"vsubpd	%%ymm7,%%ymm4,%%ymm4					\n\t	vmovaps	%%ymm13,%%ymm15	\n\t"\
		"vmovaps	%%ymm5,%%ymm7						\n\t	vmovaps	%%ymm12,%%ymm14	\n\t"\
		"vmovaps	%%ymm4,%%ymm6						\n\t	vsubpd	0x100(%%rsi),%%ymm12,%%ymm12	\n\t"\
		"vsubpd	     (%%rsi),%%ymm4,%%ymm4				\n\t	vsubpd	0x120(%%rsi),%%ymm13,%%ymm13	\n\t"\
		"vsubpd	0x020(%%rsi),%%ymm5,%%ymm5				\n\t	vaddpd	0x100(%%rsi),%%ymm14,%%ymm14	\n\t"\
		"vaddpd	     (%%rsi),%%ymm6,%%ymm6				\n\t	vaddpd	0x120(%%rsi),%%ymm15,%%ymm15	\n\t"\
		"vaddpd	0x020(%%rsi),%%ymm7,%%ymm7				\n\t	vsubpd	     %%ymm14,%%ymm8 ,%%ymm8 	\n\t"\
		"vsubpd	%%ymm6,%%ymm0,%%ymm0					\n\t	vsubpd	     %%ymm15,%%ymm9 ,%%ymm9 	\n\t"\
		"vsubpd	%%ymm5,%%ymm2,%%ymm2					\n\t	vmovaps	%%ymm8 ,0x180(%%rsi)	\n\t"\
		"vsubpd	%%ymm7,%%ymm1,%%ymm1					\n\t	vsubpd	%%ymm13,%%ymm10,%%ymm10	\n\t"\
		"vsubpd	%%ymm4,%%ymm3,%%ymm3					\n\t/*	vmovaps	%%ymm9 ,0x1a0(%%rsi)	*/\n\t"\
		"/*vmovaps	%%ymm0,0x080(%%rsi)	*/				\n\t	vsubpd	%%ymm12,%%ymm11,%%ymm11	\n\t"\
		"/*vmovaps	%%ymm2,0x040(%%rsi)	*/				\n\t	vmulpd	(%%r8) ,%%ymm14,%%ymm14	\n\t"\
		"/*vmovaps	%%ymm1,0x0a0(%%rsi)	*/				\n\t	vmulpd	(%%r8) ,%%ymm13,%%ymm13	\n\t"\
		"/*vmovaps	%%ymm3,0x0e0(%%rsi)	*/				\n\t	vmulpd	(%%r8) ,%%ymm15,%%ymm15	\n\t"\
		"vmulpd	(%%r8),%%ymm6,%%ymm6					\n\t	vmulpd	(%%r8) ,%%ymm12,%%ymm12	\n\t"\
		"vmulpd	(%%r8),%%ymm5,%%ymm5					\n\t	vaddpd	%%ymm8 ,%%ymm14,%%ymm14	\n\t"\
		"vmulpd	(%%r8),%%ymm7,%%ymm7					\n\t	vaddpd	%%ymm10,%%ymm13,%%ymm13	\n\t"\
		"vmulpd	(%%r8),%%ymm4,%%ymm4					\n\t	vaddpd	%%ymm9 ,%%ymm15,%%ymm15	\n\t"\
		"vaddpd	%%ymm0,%%ymm6,%%ymm6					\n\t	vaddpd	%%ymm11,%%ymm12,%%ymm12	\n\t"\
		"vaddpd	%%ymm2,%%ymm5,%%ymm5					\n\t	vmovaps	%%ymm14,0x100(%%rsi)	\n\t"\
		"vaddpd	%%ymm1,%%ymm7,%%ymm7					\n\t	vmovaps	%%ymm15,0x120(%%rsi)	\n\t"\
		"vaddpd	%%ymm3,%%ymm4,%%ymm4					\n\t	vmovaps	0x600(%%rsi),%%ymm8 	/* isrt2 */	\n\t"\
		"/*vmovaps	%%ymm6,     (%%rsi)	*/				\n\t	vmovaps	%%ymm10,%%ymm14	\n\t"\
		"/*vmovaps	%%ymm5,0x0c0(%%rsi)	*/				\n\t	vmovaps	%%ymm13,%%ymm15	\n\t"\
		"/*vmovaps	%%ymm7,0x020(%%rsi)	*/				\n\t	vsubpd	%%ymm12,%%ymm10,%%ymm10	\n\t"\
		"/*vmovaps	%%ymm4,0x060(%%rsi)	*/				\n\t	vsubpd	%%ymm11,%%ymm13,%%ymm13	\n\t"\
		"														vaddpd	%%ymm12,%%ymm14,%%ymm14	\n\t"\
		"														vaddpd	%%ymm11,%%ymm15,%%ymm15	\n\t"\
		"														vmulpd	%%ymm8 ,%%ymm10,%%ymm10	\n\t"\
		"														vmulpd	%%ymm8 ,%%ymm13,%%ymm13	\n\t"\
		"														vmulpd	%%ymm8 ,%%ymm14,%%ymm14	\n\t"\
		"vsubpd	%%ymm10,%%ymm2,%%ymm2					\n\t	vmulpd	%%ymm8 ,%%ymm15,%%ymm15	\n\t"\
		"													/*	vmovaps	%%ymm10,0x140(%%rsi)	*/\n\t"\
		"													/*	vmovaps	%%ymm13,0x1c0(%%rsi)	*/\n\t"\
		"													/*	vmovaps	%%ymm14,0x160(%%rsi)	*/\n\t"\
		"													/*	vmovaps	%%ymm15,0x1e0(%%rsi)	*/\n\t"\
		"/* SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r10) */	\n\t"\
		"vmovaps	0x100(%%rsi),%%ymm11		\n\t"\
		"vmovaps	0x120(%%rsi),%%ymm12		\n\t"\
		"						vmulpd	(%%r8),%%ymm10,%%ymm10	\n\t"\
		"vmovaps	0x180(%%rsi),%%ymm8 		\n\t"\
		"vsubpd	%%ymm11,%%ymm6,%%ymm6					\n\t/*	vsubpd	%%ymm10,%%ymm2,%%ymm2	*/\n\t"\
		"vsubpd	%%ymm9 ,%%ymm0,%%ymm0					\n\t	vsubpd	%%ymm15,%%ymm5,%%ymm5	\n\t"\
		"						vaddpd	%%ymm2 ,%%ymm10,%%ymm10	\n\t"\
		"						vmovaps	%%ymm2 ,0x140(%%rsi)		\n\t"\
		"vsubpd	%%ymm12,%%ymm7,%%ymm7					\n\t	vsubpd	%%ymm14,%%ymm4,%%ymm4	\n\t"\
		"						vmovaps	(%%r8),%%ymm2	/* 2.0 */	\n\t"\
		"vsubpd	%%ymm8 ,%%ymm1 ,%%ymm1 					\n\t	vsubpd	%%ymm13,%%ymm3,%%ymm3	\n\t"\
		"vmulpd	%%ymm2 ,%%ymm11,%%ymm11					\n\t/*	vaddpd	%%ymm10,%%ymm10,%%ymm10	*/\n\t"\
		"vmulpd	%%ymm2 ,%%ymm9 ,%%ymm9 					\n\t	vmulpd	%%ymm2 ,%%ymm15,%%ymm15	\n\t"\
		"vmulpd	%%ymm2 ,%%ymm12,%%ymm12					\n\t	vmulpd	%%ymm2 ,%%ymm14,%%ymm14	\n\t"\
		"vmulpd	%%ymm2 ,%%ymm8 ,%%ymm8 					\n\t	vmulpd	%%ymm2 ,%%ymm13,%%ymm13	\n\t"\
		"vaddpd	%%ymm6 ,%%ymm11,%%ymm11					\n\t/*	vaddpd	%%ymm2 ,%%ymm10,%%ymm10	*/\n\t"\
		"vaddpd	%%ymm0 ,%%ymm9 ,%%ymm9 					\n\t	vaddpd	%%ymm5 ,%%ymm15,%%ymm15	\n\t"\
		"vaddpd	%%ymm7 ,%%ymm12,%%ymm12					\n\t	vaddpd	%%ymm4 ,%%ymm14,%%ymm14	\n\t"\
		"vaddpd	%%ymm1 ,%%ymm8 ,%%ymm8 					\n\t	vaddpd	%%ymm3 ,%%ymm13,%%ymm13	\n\t"\
		"vmovaps	%%ymm6 ,0x100(%%rsi)				\n\t/*	vmovaps	%%ymm2 ,0x140(%%rsi)	*/\n\t"\
		"vmovaps	%%ymm0 ,0x080(%%rsi)				\n\t	vmovaps	%%ymm5 ,0x0c0(%%rsi)	\n\t"\
		"vmovaps	%%ymm7 ,0x120(%%rsi)				\n\t	vmovaps	%%ymm4 ,0x160(%%rsi)	\n\t"\
		"vmovaps	%%ymm1 ,0x1a0(%%rsi)				\n\t	vmovaps	%%ymm3 ,0x1e0(%%rsi)	\n\t"\
		"vmovaps	%%ymm11,     (%%rsi)				\n\t	vmovaps	%%ymm10,0x040(%%rsi)	\n\t"\
		"vmovaps	%%ymm9 ,0x180(%%rsi)				\n\t	vmovaps	%%ymm15,0x1c0(%%rsi)	\n\t"\
		"vmovaps	%%ymm12,0x020(%%rsi)				\n\t	vmovaps	%%ymm14,0x060(%%rsi)	\n\t"\
		"vmovaps	%%ymm8 ,0x0a0(%%rsi)				\n\t	vmovaps	%%ymm13,0x0e0(%%rsi)	\n\t"\
		"\n\t"\
		"/***************************************/\n\t"\
		"\n\t"\
		"/*...Block 3: */	\n\t"\
		"movq	%[__add0],%%rax	\n\t"\
		"movslq	%[__p01],%%rdi	/* Do things this way [rather than repeatedly adding p1] since array-padding scheme means p2 == p+p1 not guaranteed. */	\n\t"\
		"movslq	%[__p08],%%rbx	/* Can't get these via simple load-one-and-shift-as-needed due to array padding scheme */	\n\t"\
		"movslq	%[__p10],%%rcx	\n\t"\
		"movslq	%[__p18],%%rdx	\n\t"\
		"shlq	$3,%%rdi	\n\t"\
		"shlq	$3,%%rbx	\n\t"\
		"shlq	$3,%%rcx								\n\t	movslq	%[__p04],%%r9	\n\t"\
		"shlq	$3,%%rdx								\n\t	shlq	$3,%%r9	\n\t"\
		"addq	%%rdi,%%rax	/* &a[j1+p1] */				\n\t	movq	%%r9,%%r10	\n\t	addq	%%rax,%%r10	\n\t"\
		"addq	%%rax,%%rbx								\n\t	movq	%%r9,%%r11	\n\t	addq	%%rbx,%%r11	\n\t"\
		"addq	%%rax,%%rcx								\n\t	movq	%%r9,%%r12	\n\t	addq	%%rcx,%%r12	\n\t"\
		"addq	%%rax,%%rdx								\n\t	movq	%%r9,%%r13	\n\t	addq	%%rdx,%%r13	\n\t"\
		"/* SSE2_RADIX4_DIF_4TWIDDLE_B(r20,c01) */		\n\t	/* SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO_B(r28,c05) */\n\t"\
		"addq	$0x200,%%rsi	/* r20 */				\n\t	vmovaps	     (%%r10),%%ymm8 	\n\t"\
		"vmovaps	     (%%rax),%%ymm0					\n\t	vmovaps	     (%%r12),%%ymm12	\n\t"\
		"vmovaps	     (%%rcx),%%ymm4					\n\t	vmovaps	0x020(%%r10),%%ymm9 	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1					\n\t	vmovaps	0x020(%%r12),%%ymm13	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5					\n\t	vmovaps	0x9e0(%%rsi),%%ymm14	/* c05 */	\n\t"\
		"vmovaps	0x8e0(%%rsi),%%ymm6	/* c01 */		\n\t	vmovaps	0xa00(%%rsi),%%ymm15	\n\t"\
		"vmovaps	0x900(%%rsi),%%ymm7					\n\t	vmovaps	%%ymm8 ,%%ymm10	\n\t"\
		"vmovaps	%%ymm0,%%ymm2						\n\t	vmovaps	%%ymm9 ,%%ymm11	\n\t"\
		"vmovaps	%%ymm1,%%ymm3						\n\t	vmulpd	%%ymm14,%%ymm8 ,%%ymm8 	\n\t"\
		"vmulpd	%%ymm6,%%ymm0,%%ymm0					\n\t	vmulpd	%%ymm14,%%ymm9 ,%%ymm9 	\n\t"\
		"vmulpd	%%ymm6,%%ymm1,%%ymm1					\n\t	vmulpd	%%ymm15,%%ymm10,%%ymm10	\n\t"\
		"vmulpd	%%ymm7,%%ymm2,%%ymm2					\n\t	vmulpd	%%ymm15,%%ymm11,%%ymm11	\n\t"\
		"vmulpd	%%ymm7,%%ymm3,%%ymm3					\n\t	vmovaps	%%ymm12,%%ymm14	\n\t"\
		"vmovaps	%%ymm4,%%ymm6						\n\t	vaddpd	%%ymm10,%%ymm9 ,%%ymm9 	\n\t"\
		"vaddpd	%%ymm2,%%ymm1,%%ymm1					\n\t	vmovaps	%%ymm13,%%ymm15	\n\t"\
		"vmovaps	%%ymm5,%%ymm7						\n\t	vmulpd	0xa20(%%rsi),%%ymm12,%%ymm12	/* c15 */	\n\t"\
		"vmulpd	0x920(%%rsi),%%ymm4,%%ymm4	/* c11 */	\n\t	vsubpd	     %%ymm11,%%ymm8 ,%%ymm8 	\n\t"\
		"vsubpd	%%ymm3,%%ymm0,%%ymm0					\n\t	vmulpd	0xa20(%%rsi),%%ymm13,%%ymm13	\n\t"\
		"vmulpd	0x920(%%rsi),%%ymm5,%%ymm5				\n\t	vmulpd	0xa40(%%rsi),%%ymm14,%%ymm14	\n\t"\
		"vmulpd	0x940(%%rsi),%%ymm6,%%ymm6				\n\t	vmovaps	%%ymm8 ,%%ymm10	\n\t"\
		"vmovaps	%%ymm0,%%ymm2						\n\t	vmulpd	0xa40(%%rsi),%%ymm15,%%ymm15	\n\t"\
		"vmulpd	0x940(%%rsi),%%ymm7,%%ymm7				\n\t	vaddpd	%%ymm14,%%ymm13,%%ymm13	\n\t"\
		"vaddpd	%%ymm6,%%ymm5,%%ymm5					\n\t	vmovaps	%%ymm9 ,%%ymm11	\n\t"\
		"vmovaps	%%ymm1,%%ymm3						\n\t	vsubpd	%%ymm15,%%ymm12,%%ymm12	\n\t"\
		"vsubpd	%%ymm7,%%ymm4,%%ymm4					\n\t	vaddpd	%%ymm12,%%ymm8 ,%%ymm8 	\n\t"\
		"vaddpd	%%ymm4,%%ymm0,%%ymm0					\n\t	vaddpd	%%ymm13,%%ymm9 ,%%ymm9 	\n\t"\
		"vaddpd	%%ymm5,%%ymm1,%%ymm1					\n\t	vsubpd	%%ymm12,%%ymm10,%%ymm10	\n\t"\
		"vsubpd	%%ymm4,%%ymm2,%%ymm2					\n\t	vsubpd	%%ymm13,%%ymm11,%%ymm11	\n\t"\
		"vsubpd	%%ymm5,%%ymm3,%%ymm3					\n\t	vmovaps	     (%%r13),%%ymm12	\n\t"\
		"vmovaps	     (%%rdx),%%ymm4					\n\t	vmovaps	0x020(%%r13),%%ymm13	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm5					\n\t	vmovaps	     (%%r13),%%ymm14	\n\t"\
		"vmovaps	     (%%rdx),%%ymm6					\n\t	vmovaps	0x020(%%r13),%%ymm15	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm7					\n\t	vmulpd	0xaa0(%%rsi),%%ymm12,%%ymm12	/* c1D */	\n\t"\
		"vmulpd	0x9a0(%%rsi),%%ymm4,%%ymm4	/* c19 */	\n\t	vmulpd	0xaa0(%%rsi),%%ymm13,%%ymm13	\n\t"\
		"vmulpd	0x9a0(%%rsi),%%ymm5,%%ymm5				\n\t	vmulpd	0xac0(%%rsi),%%ymm14,%%ymm14	\n\t"\
		"vmulpd	0x9c0(%%rsi),%%ymm6,%%ymm6				\n\t	vmulpd	0xac0(%%rsi),%%ymm15,%%ymm15	\n\t"\
		"vmulpd	0x9c0(%%rsi),%%ymm7,%%ymm7				\n\t	vaddpd	     %%ymm14,%%ymm13,%%ymm13	\n\t"\
		"vaddpd	%%ymm6,%%ymm5,%%ymm5					\n\t	vsubpd	     %%ymm15,%%ymm12,%%ymm12	\n\t"\
		"vsubpd	%%ymm7,%%ymm4,%%ymm4					\n\t	vmovaps	%%ymm13,0x120(%%rsi)	\n\t"\
		"vmovaps	%%ymm5,0x020(%%rsi)					\n\t	vmovaps	%%ymm12,0x100(%%rsi)	/* tmpstr r08 */\n\t"\
		"vmovaps	%%ymm4,     (%%rsi)					\n\t	vmovaps	     (%%r11),%%ymm12	\n\t"\
		"vmovaps	     (%%rbx),%%ymm4					\n\t	vmovaps	0x020(%%r11),%%ymm13	\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm5					\n\t	vmovaps	     (%%r11),%%ymm14	\n\t"\
		"vmovaps	     (%%rbx),%%ymm6					\n\t	vmovaps	0x020(%%r11),%%ymm15	\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm7					\n\t	vmulpd	0xa60(%%rsi),%%ymm12,%%ymm12	/* c0D */	\n\t"\
		"vmulpd	0x960(%%rsi),%%ymm4,%%ymm4	/* c01 */	\n\t	vmulpd	0xa60(%%rsi),%%ymm13,%%ymm13	\n\t"\
		"vmulpd	0x960(%%rsi),%%ymm5,%%ymm5				\n\t	vmulpd	0xa80(%%rsi),%%ymm14,%%ymm14	\n\t"\
		"vmulpd	0x980(%%rsi),%%ymm6,%%ymm6				\n\t	vmulpd	0xa80(%%rsi),%%ymm15,%%ymm15	\n\t"\
		"vmulpd	0x980(%%rsi),%%ymm7,%%ymm7				\n\t	vaddpd	     %%ymm14,%%ymm13,%%ymm13	\n\t"\
		"vaddpd	%%ymm6,%%ymm5,%%ymm5					\n\t	vsubpd	     %%ymm15,%%ymm12,%%ymm12	\n\t"\
		"vsubpd	%%ymm7,%%ymm4,%%ymm4					\n\t	vmovaps	%%ymm13,%%ymm15	\n\t"\
		"vmovaps	%%ymm5,%%ymm7						\n\t	vmovaps	%%ymm12,%%ymm14	\n\t"\
		"vmovaps	%%ymm4,%%ymm6						\n\t	vsubpd	0x100(%%rsi),%%ymm12,%%ymm12	\n\t"\
		"vsubpd	     (%%rsi),%%ymm4,%%ymm4				\n\t	vsubpd	0x120(%%rsi),%%ymm13,%%ymm13	\n\t"\
		"vsubpd	0x020(%%rsi),%%ymm5,%%ymm5				\n\t	vaddpd	0x100(%%rsi),%%ymm14,%%ymm14	\n\t"\
		"vaddpd	     (%%rsi),%%ymm6,%%ymm6				\n\t	vaddpd	0x120(%%rsi),%%ymm15,%%ymm15	\n\t"\
		"vaddpd	0x020(%%rsi),%%ymm7,%%ymm7				\n\t	vsubpd	     %%ymm14,%%ymm8 ,%%ymm8 	\n\t"\
		"vsubpd	%%ymm6,%%ymm0,%%ymm0					\n\t	vsubpd	     %%ymm15,%%ymm9 ,%%ymm9 	\n\t"\
		"vsubpd	%%ymm5,%%ymm2,%%ymm2					\n\t	vmovaps	%%ymm8 ,0x180(%%rsi)	\n\t"\
		"vsubpd	%%ymm7,%%ymm1,%%ymm1					\n\t	vsubpd	%%ymm13,%%ymm10,%%ymm10	\n\t"\
		"vsubpd	%%ymm4,%%ymm3,%%ymm3					\n\t/*	vmovaps	%%ymm9 ,0x1a0(%%rsi)	*/\n\t"\
		"/*vmovaps	%%ymm0,0x080(%%rsi)	*/				\n\t	vsubpd	%%ymm12,%%ymm11,%%ymm11	\n\t"\
		"/*vmovaps	%%ymm2,0x040(%%rsi)	*/				\n\t	vmulpd	(%%r8) ,%%ymm14,%%ymm14	\n\t"\
		"/*vmovaps	%%ymm1,0x0a0(%%rsi)	*/				\n\t	vmulpd	(%%r8) ,%%ymm13,%%ymm13	\n\t"\
		"/*vmovaps	%%ymm3,0x0e0(%%rsi)	*/				\n\t	vmulpd	(%%r8) ,%%ymm15,%%ymm15	\n\t"\
		"vmulpd	(%%r8),%%ymm6,%%ymm6					\n\t	vmulpd	(%%r8) ,%%ymm12,%%ymm12	\n\t"\
		"vmulpd	(%%r8),%%ymm5,%%ymm5					\n\t	vaddpd	%%ymm8 ,%%ymm14,%%ymm14	\n\t"\
		"vmulpd	(%%r8),%%ymm7,%%ymm7					\n\t	vaddpd	%%ymm10,%%ymm13,%%ymm13	\n\t"\
		"vmulpd	(%%r8),%%ymm4,%%ymm4					\n\t	vaddpd	%%ymm9 ,%%ymm15,%%ymm15	\n\t"\
		"vaddpd	%%ymm0,%%ymm6,%%ymm6					\n\t	vaddpd	%%ymm11,%%ymm12,%%ymm12	\n\t"\
		"vaddpd	%%ymm2,%%ymm5,%%ymm5					\n\t	vmovaps	%%ymm14,0x100(%%rsi)	\n\t"\
		"vaddpd	%%ymm1,%%ymm7,%%ymm7					\n\t	vmovaps	%%ymm15,0x120(%%rsi)	\n\t"\
		"vaddpd	%%ymm3,%%ymm4,%%ymm4					\n\t	vmovaps	0x400(%%rsi),%%ymm8 	/* isrt2 */	\n\t"\
		"/*vmovaps	%%ymm6,     (%%rsi)	*/				\n\t	vmovaps	%%ymm10,%%ymm14	\n\t"\
		"/*vmovaps	%%ymm5,0x0c0(%%rsi)	*/				\n\t	vmovaps	%%ymm13,%%ymm15	\n\t"\
		"/*vmovaps	%%ymm7,0x020(%%rsi)	*/				\n\t	vsubpd	%%ymm12,%%ymm10,%%ymm10	\n\t"\
		"/*vmovaps	%%ymm4,0x060(%%rsi)	*/				\n\t	vsubpd	%%ymm11,%%ymm13,%%ymm13	\n\t"\
		"														vaddpd	%%ymm12,%%ymm14,%%ymm14	\n\t"\
		"														vaddpd	%%ymm11,%%ymm15,%%ymm15	\n\t"\
		"														vmulpd	%%ymm8 ,%%ymm10,%%ymm10	\n\t"\
		"														vmulpd	%%ymm8 ,%%ymm13,%%ymm13	\n\t"\
		"														vmulpd	%%ymm8 ,%%ymm14,%%ymm14	\n\t"\
		"vsubpd	%%ymm10,%%ymm2,%%ymm2					\n\t	vmulpd	%%ymm8 ,%%ymm15,%%ymm15	\n\t"\
		"													/*	vmovaps	%%ymm10,0x140(%%rsi)	*/\n\t"\
		"													/*	vmovaps	%%ymm13,0x1c0(%%rsi)	*/\n\t"\
		"													/*	vmovaps	%%ymm14,0x160(%%rsi)	*/\n\t"\
		"													/*	vmovaps	%%ymm15,0x1e0(%%rsi)	*/\n\t"\
		"/* SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r20) */	\n\t"\
		"vmovaps	0x100(%%rsi),%%ymm11		\n\t"\
		"vmovaps	0x120(%%rsi),%%ymm12		\n\t"\
		"						vmulpd	(%%r8),%%ymm10,%%ymm10	\n\t"\
		"vmovaps	0x180(%%rsi),%%ymm8 		\n\t"\
		"vsubpd	%%ymm11,%%ymm6,%%ymm6					\n\t/*	vsubpd	%%ymm10,%%ymm2,%%ymm2	*/\n\t"\
		"vsubpd	%%ymm9 ,%%ymm0,%%ymm0					\n\t	vsubpd	%%ymm15,%%ymm5,%%ymm5	\n\t"\
		"						vaddpd	%%ymm2 ,%%ymm10,%%ymm10	\n\t"\
		"						vmovaps	%%ymm2 ,0x140(%%rsi)	\n\t"\
		"vsubpd	%%ymm12,%%ymm7,%%ymm7					\n\t	vsubpd	%%ymm14,%%ymm4,%%ymm4	\n\t"\
		"						vmovaps	(%%r8),%%ymm2	/* 2.0 */	\n\t"\
		"vsubpd	%%ymm8 ,%%ymm1 ,%%ymm1 					\n\t	vsubpd	%%ymm13,%%ymm3,%%ymm3	\n\t"\
		"vmulpd	%%ymm2 ,%%ymm11,%%ymm11					\n\t/*	vaddpd	%%ymm10,%%ymm10,%%ymm10	*/\n\t"\
		"vmulpd	%%ymm2 ,%%ymm9 ,%%ymm9 					\n\t	vmulpd	%%ymm2 ,%%ymm15,%%ymm15	\n\t"\
		"vmulpd	%%ymm2 ,%%ymm12,%%ymm12					\n\t	vmulpd	%%ymm2 ,%%ymm14,%%ymm14	\n\t"\
		"vmulpd	%%ymm2 ,%%ymm8 ,%%ymm8 					\n\t	vmulpd	%%ymm2 ,%%ymm13,%%ymm13	\n\t"\
		"vaddpd	%%ymm6 ,%%ymm11,%%ymm11					\n\t/*	vaddpd	%%ymm2 ,%%ymm10,%%ymm10	*/\n\t"\
		"vaddpd	%%ymm0 ,%%ymm9 ,%%ymm9 					\n\t	vaddpd	%%ymm5 ,%%ymm15,%%ymm15	\n\t"\
		"vaddpd	%%ymm7 ,%%ymm12,%%ymm12					\n\t	vaddpd	%%ymm4 ,%%ymm14,%%ymm14	\n\t"\
		"vaddpd	%%ymm1 ,%%ymm8 ,%%ymm8 					\n\t	vaddpd	%%ymm3 ,%%ymm13,%%ymm13	\n\t"\
		"vmovaps	%%ymm6 ,0x100(%%rsi)				\n\t/*	vmovaps	%%ymm2 ,0x140(%%rsi)	*/\n\t"\
		"vmovaps	%%ymm0 ,0x080(%%rsi)				\n\t	vmovaps	%%ymm5 ,0x0c0(%%rsi)	\n\t"\
		"vmovaps	%%ymm7 ,0x120(%%rsi)				\n\t	vmovaps	%%ymm4 ,0x160(%%rsi)	\n\t"\
		"vmovaps	%%ymm1 ,0x1a0(%%rsi)				\n\t	vmovaps	%%ymm3 ,0x1e0(%%rsi)	\n\t"\
		"vmovaps	%%ymm11,     (%%rsi)				\n\t	vmovaps	%%ymm10,0x040(%%rsi)	\n\t"\
		"vmovaps	%%ymm9 ,0x180(%%rsi)				\n\t	vmovaps	%%ymm15,0x1c0(%%rsi)	\n\t"\
		"vmovaps	%%ymm12,0x020(%%rsi)				\n\t	vmovaps	%%ymm14,0x060(%%rsi)	\n\t"\
		"vmovaps	%%ymm8 ,0x0a0(%%rsi)				\n\t	vmovaps	%%ymm13,0x0e0(%%rsi)	\n\t"\
		"\n\t"\
		"/***************************************/\n\t"\
		"\n\t"\
		"/*...Block 4: */	\n\t"\
		"movq	%[__add0],%%rax	\n\t"\
		"movslq	%[__p03],%%rdi	/* Do things this way [rather than repeatedly adding p1] since array-padding scheme means p2 == p+p1 not guaranteed. */	\n\t"\
		"movslq	%[__p08],%%rbx	/* Can't get these via simple load-one-and-shift-as-needed due to array padding scheme */	\n\t"\
		"movslq	%[__p10],%%rcx	\n\t"\
		"movslq	%[__p18],%%rdx	\n\t"\
		"shlq	$3,%%rdi	\n\t"\
		"shlq	$3,%%rbx	\n\t"\
		"shlq	$3,%%rcx								\n\t	movslq	%[__p04],%%r9	\n\t"\
		"shlq	$3,%%rdx								\n\t	shlq	$3,%%r9	\n\t"\
		"addq	%%rdi,%%rax	/* &a[j1+p3] */				\n\t	movq	%%r9,%%r10	\n\t	addq	%%rax,%%r10	\n\t"\
		"addq	%%rax,%%rbx								\n\t	movq	%%r9,%%r11	\n\t	addq	%%rbx,%%r11	\n\t"\
		"addq	%%rax,%%rcx								\n\t	movq	%%r9,%%r12	\n\t	addq	%%rcx,%%r12	\n\t"\
		"addq	%%rax,%%rdx								\n\t	movq	%%r9,%%r13	\n\t	addq	%%rdx,%%r13	\n\t"\
		"/* SSE2_RADIX4_DIF_4TWIDDLE_B(r30,c03) */		\n\t	/* SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO_B(r38,c07) */\n\t"\
		"addq	$0x200,%%rsi	/* r30 */				\n\t	vmovaps	     (%%r10),%%ymm8 	\n\t"\
		"vmovaps	     (%%rax),%%ymm0					\n\t	vmovaps	     (%%r12),%%ymm12	\n\t"\
		"vmovaps	     (%%rcx),%%ymm4					\n\t	vmovaps	0x020(%%r10),%%ymm9 	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1					\n\t	vmovaps	0x020(%%r12),%%ymm13	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5					\n\t	vmovaps	0x9e0(%%rsi),%%ymm14	/* c07 */	\n\t"\
		"vmovaps	0x8e0(%%rsi),%%ymm6	/* c03 */		\n\t	vmovaps	0xa00(%%rsi),%%ymm15	\n\t"\
		"vmovaps	0x900(%%rsi),%%ymm7					\n\t	vmovaps	%%ymm8 ,%%ymm10	\n\t"\
		"vmovaps	%%ymm0,%%ymm2						\n\t	vmovaps	%%ymm9 ,%%ymm11	\n\t"\
		"vmovaps	%%ymm1,%%ymm3						\n\t	vmulpd	%%ymm14,%%ymm8 ,%%ymm8 	\n\t"\
		"vmulpd	%%ymm6,%%ymm0,%%ymm0					\n\t	vmulpd	%%ymm14,%%ymm9 ,%%ymm9 	\n\t"\
		"vmulpd	%%ymm6,%%ymm1,%%ymm1					\n\t	vmulpd	%%ymm15,%%ymm10,%%ymm10	\n\t"\
		"vmulpd	%%ymm7,%%ymm2,%%ymm2					\n\t	vmulpd	%%ymm15,%%ymm11,%%ymm11	\n\t"\
		"vmulpd	%%ymm7,%%ymm3,%%ymm3					\n\t	vmovaps	%%ymm12,%%ymm14	\n\t"\
		"vmovaps	%%ymm4,%%ymm6						\n\t	vaddpd	%%ymm10,%%ymm9 ,%%ymm9 	\n\t"\
		"vaddpd	%%ymm2,%%ymm1,%%ymm1					\n\t	vmovaps	%%ymm13,%%ymm15	\n\t"\
		"vmovaps	%%ymm5,%%ymm7						\n\t	vmulpd	0xa20(%%rsi),%%ymm12,%%ymm12	/* c17 */	\n\t"\
		"vmulpd	0x920(%%rsi),%%ymm4,%%ymm4	/* c13 */	\n\t	vsubpd	     %%ymm11,%%ymm8 ,%%ymm8 	\n\t"\
		"vsubpd	%%ymm3,%%ymm0,%%ymm0					\n\t	vmulpd	0xa20(%%rsi),%%ymm13,%%ymm13	\n\t"\
		"vmulpd	0x920(%%rsi),%%ymm5,%%ymm5				\n\t	vmulpd	0xa40(%%rsi),%%ymm14,%%ymm14	\n\t"\
		"vmulpd	0x940(%%rsi),%%ymm6,%%ymm6				\n\t	vmovaps	%%ymm8 ,%%ymm10	\n\t"\
		"vmovaps	%%ymm0,%%ymm2						\n\t	vmulpd	0xa40(%%rsi),%%ymm15,%%ymm15	\n\t"\
		"vmulpd	0x940(%%rsi),%%ymm7,%%ymm7				\n\t	vaddpd	%%ymm14,%%ymm13,%%ymm13	\n\t"\
		"vaddpd	%%ymm6,%%ymm5,%%ymm5					\n\t	vmovaps	%%ymm9 ,%%ymm11	\n\t"\
		"vmovaps	%%ymm1,%%ymm3						\n\t	vsubpd	%%ymm15,%%ymm12,%%ymm12	\n\t"\
		"vsubpd	%%ymm7,%%ymm4,%%ymm4					\n\t	vaddpd	%%ymm12,%%ymm8 ,%%ymm8 	\n\t"\
		"vaddpd	%%ymm4,%%ymm0,%%ymm0					\n\t	vaddpd	%%ymm13,%%ymm9 ,%%ymm9 	\n\t"\
		"vaddpd	%%ymm5,%%ymm1,%%ymm1					\n\t	vsubpd	%%ymm12,%%ymm10,%%ymm10	\n\t"\
		"vsubpd	%%ymm4,%%ymm2,%%ymm2					\n\t	vsubpd	%%ymm13,%%ymm11,%%ymm11	\n\t"\
		"vsubpd	%%ymm5,%%ymm3,%%ymm3					\n\t	vmovaps	     (%%r13),%%ymm12	\n\t"\
		"vmovaps	     (%%rdx),%%ymm4					\n\t	vmovaps	0x020(%%r13),%%ymm13	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm5					\n\t	vmovaps	     (%%r13),%%ymm14	\n\t"\
		"vmovaps	     (%%rdx),%%ymm6					\n\t	vmovaps	0x020(%%r13),%%ymm15	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm7					\n\t	vmulpd	0xaa0(%%rsi),%%ymm12,%%ymm12	/* c1F */	\n\t"\
		"vmulpd	0x9a0(%%rsi),%%ymm4,%%ymm4	/* c1B */	\n\t	vmulpd	0xaa0(%%rsi),%%ymm13,%%ymm13	\n\t"\
		"vmulpd	0x9a0(%%rsi),%%ymm5,%%ymm5				\n\t	vmulpd	0xac0(%%rsi),%%ymm14,%%ymm14	\n\t"\
		"vmulpd	0x9c0(%%rsi),%%ymm6,%%ymm6				\n\t	vmulpd	0xac0(%%rsi),%%ymm15,%%ymm15	\n\t"\
		"vmulpd	0x9c0(%%rsi),%%ymm7,%%ymm7				\n\t	vaddpd	     %%ymm14,%%ymm13,%%ymm13	\n\t"\
		"vaddpd	%%ymm6,%%ymm5,%%ymm5					\n\t	vsubpd	     %%ymm15,%%ymm12,%%ymm12	\n\t"\
		"vsubpd	%%ymm7,%%ymm4,%%ymm4					\n\t	vmovaps	%%ymm13,0x120(%%rsi)	\n\t"\
		"vmovaps	%%ymm5,0x020(%%rsi)					\n\t	vmovaps	%%ymm12,0x100(%%rsi)	/* tmpstr r08 */\n\t"\
		"vmovaps	%%ymm4,     (%%rsi)					\n\t	vmovaps	     (%%r11),%%ymm12	\n\t"\
		"vmovaps	     (%%rbx),%%ymm4					\n\t	vmovaps	0x020(%%r11),%%ymm13	\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm5					\n\t	vmovaps	     (%%r11),%%ymm14	\n\t"\
		"vmovaps	     (%%rbx),%%ymm6					\n\t	vmovaps	0x020(%%r11),%%ymm15	\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm7					\n\t	vmulpd	0xa60(%%rsi),%%ymm12,%%ymm12	/* c0F */	\n\t"\
		"vmulpd	0x960(%%rsi),%%ymm4,%%ymm4	/* c0B */	\n\t	vmulpd	0xa60(%%rsi),%%ymm13,%%ymm13	\n\t"\
		"vmulpd	0x960(%%rsi),%%ymm5,%%ymm5				\n\t	vmulpd	0xa80(%%rsi),%%ymm14,%%ymm14	\n\t"\
		"vmulpd	0x980(%%rsi),%%ymm6,%%ymm6				\n\t	vmulpd	0xa80(%%rsi),%%ymm15,%%ymm15	\n\t"\
		"vmulpd	0x980(%%rsi),%%ymm7,%%ymm7				\n\t	vaddpd	     %%ymm14,%%ymm13,%%ymm13	\n\t"\
		"vaddpd	%%ymm6,%%ymm5,%%ymm5					\n\t	vsubpd	     %%ymm15,%%ymm12,%%ymm12	\n\t"\
		"vsubpd	%%ymm7,%%ymm4,%%ymm4					\n\t	vmovaps	%%ymm13,%%ymm15	\n\t"\
		"vmovaps	%%ymm5,%%ymm7						\n\t	vmovaps	%%ymm12,%%ymm14	\n\t"\
		"vmovaps	%%ymm4,%%ymm6						\n\t	vsubpd	0x100(%%rsi),%%ymm12,%%ymm12	\n\t"\
		"vsubpd	     (%%rsi),%%ymm4,%%ymm4				\n\t	vsubpd	0x120(%%rsi),%%ymm13,%%ymm13	\n\t"\
		"vsubpd	0x020(%%rsi),%%ymm5,%%ymm5				\n\t	vaddpd	0x100(%%rsi),%%ymm14,%%ymm14	\n\t"\
		"vaddpd	     (%%rsi),%%ymm6,%%ymm6				\n\t	vaddpd	0x120(%%rsi),%%ymm15,%%ymm15	\n\t"\
		"vaddpd	0x020(%%rsi),%%ymm7,%%ymm7				\n\t	vsubpd	%%ymm14,%%ymm8 ,%%ymm8 	\n\t"\
		"vsubpd	%%ymm6,%%ymm0,%%ymm0					\n\t	vsubpd	%%ymm15,%%ymm9 ,%%ymm9 	\n\t"\
		"vsubpd	%%ymm5,%%ymm2,%%ymm2					\n\t	vmovaps	%%ymm8 ,0x180(%%rsi)	\n\t"\
		"vsubpd	%%ymm7,%%ymm1,%%ymm1					\n\t	vsubpd	%%ymm13,%%ymm10,%%ymm10	\n\t"\
		"vsubpd	%%ymm4,%%ymm3,%%ymm3					\n\t/*	vmovaps	%%ymm9 ,0x1a0(%%rsi)	*/\n\t"\
		"/*vmovaps	%%ymm0,0x080(%%rsi)	*/				\n\t	vsubpd	%%ymm12,%%ymm11,%%ymm11	\n\t"\
		"/*vmovaps	%%ymm2,0x040(%%rsi)	*/				\n\t	vmulpd	(%%r8) ,%%ymm14,%%ymm14	\n\t"\
		"/*vmovaps	%%ymm1,0x0a0(%%rsi)	*/				\n\t	vmulpd	(%%r8) ,%%ymm13,%%ymm13	\n\t"\
		"/*vmovaps	%%ymm3,0x0e0(%%rsi)	*/				\n\t	vmulpd	(%%r8) ,%%ymm15,%%ymm15	\n\t"\
		"vmulpd	(%%r8),%%ymm6,%%ymm6					\n\t	vmulpd	(%%r8) ,%%ymm12,%%ymm12	\n\t"\
		"vmulpd	(%%r8),%%ymm5,%%ymm5					\n\t	vaddpd	%%ymm8 ,%%ymm14,%%ymm14	\n\t"\
		"vmulpd	(%%r8),%%ymm7,%%ymm7					\n\t	vaddpd	%%ymm10,%%ymm13,%%ymm13	\n\t"\
		"vmulpd	(%%r8),%%ymm4,%%ymm4					\n\t	vaddpd	%%ymm9 ,%%ymm15,%%ymm15	\n\t"\
		"vaddpd	%%ymm0,%%ymm6,%%ymm6					\n\t	vaddpd	%%ymm11,%%ymm12,%%ymm12	\n\t"\
		"vaddpd	%%ymm2,%%ymm5,%%ymm5					\n\t	vmovaps	%%ymm14,0x100(%%rsi)	\n\t"\
		"vaddpd	%%ymm1,%%ymm7,%%ymm7					\n\t	vmovaps	%%ymm15,0x120(%%rsi)	\n\t"\
		"vaddpd	%%ymm3,%%ymm4,%%ymm4					\n\t	vmovaps	0x200(%%rsi),%%ymm8 	/* isrt2 */	\n\t"\
		"/*vmovaps	%%ymm6,     (%%rsi)	*/				\n\t	vmovaps	%%ymm10,%%ymm14	\n\t"\
		"/*vmovaps	%%ymm5,0x0c0(%%rsi)	*/				\n\t	vmovaps	%%ymm13,%%ymm15	\n\t"\
		"/*vmovaps	%%ymm7,0x020(%%rsi)	*/				\n\t	vsubpd	%%ymm12,%%ymm10,%%ymm10	\n\t"\
		"/*vmovaps	%%ymm4,0x060(%%rsi)	*/				\n\t	vsubpd	%%ymm11,%%ymm13,%%ymm13	\n\t"\
		"														vaddpd	%%ymm12,%%ymm14,%%ymm14	\n\t"\
		"														vaddpd	%%ymm11,%%ymm15,%%ymm15	\n\t"\
		"														vmulpd	%%ymm8 ,%%ymm10,%%ymm10	\n\t"\
		"														vmulpd	%%ymm8 ,%%ymm13,%%ymm13	\n\t"\
		"														vmulpd	%%ymm8 ,%%ymm14,%%ymm14	\n\t"\
		"vsubpd	%%ymm10,%%ymm2,%%ymm2					\n\t	vmulpd	%%ymm8 ,%%ymm15,%%ymm15	\n\t"\
		"													/*	vmovaps	%%ymm10,0x140(%%rsi)	*/\n\t"\
		"													/*	vmovaps	%%ymm13,0x1c0(%%rsi)	*/\n\t"\
		"													/*	vmovaps	%%ymm14,0x160(%%rsi)	*/\n\t"\
		"													/*	vmovaps	%%ymm15,0x1e0(%%rsi)	*/\n\t"\
		"/* SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r30) */	\n\t"\
		"vmovaps	0x100(%%rsi),%%ymm11		\n\t"\
		"vmovaps	0x120(%%rsi),%%ymm12		\n\t"\
		"						vmulpd	(%%r8),%%ymm10,%%ymm10	\n\t"\
		"vmovaps	0x180(%%rsi),%%ymm8 		\n\t"\
		"vsubpd	%%ymm11,%%ymm6,%%ymm6					\n\t/*	vsubpd	%%ymm10,%%ymm2,%%ymm2	*/\n\t"\
		"vsubpd	%%ymm9 ,%%ymm0,%%ymm0					\n\t	vsubpd	%%ymm15,%%ymm5,%%ymm5	\n\t"\
		"						vaddpd	%%ymm2 ,%%ymm10,%%ymm10	\n\t"\
		"						vmovaps	%%ymm2 ,0x140(%%rsi)	\n\t"\
		"vsubpd	%%ymm12,%%ymm7,%%ymm7					\n\t	vsubpd	%%ymm14,%%ymm4,%%ymm4	\n\t"\
		"						vmovaps	(%%r8),%%ymm2	/* 2.0 */	\n\t"\
		"vsubpd	%%ymm8 ,%%ymm1 ,%%ymm1 					\n\t	vsubpd	%%ymm13,%%ymm3,%%ymm3	\n\t"\
		"vmulpd	%%ymm2 ,%%ymm11,%%ymm11					\n\t/*	vaddpd	%%ymm10,%%ymm10,%%ymm10	*/\n\t"\
		"vmulpd	%%ymm2 ,%%ymm9 ,%%ymm9 					\n\t	vmulpd	%%ymm2 ,%%ymm15,%%ymm15	\n\t"\
		"vmulpd	%%ymm2 ,%%ymm12,%%ymm12					\n\t	vmulpd	%%ymm2 ,%%ymm14,%%ymm14	\n\t"\
		"vmulpd	%%ymm2 ,%%ymm8 ,%%ymm8 					\n\t	vmulpd	%%ymm2 ,%%ymm13,%%ymm13	\n\t"\
		"vaddpd	%%ymm6 ,%%ymm11,%%ymm11					\n\t/*	vaddpd	%%ymm2 ,%%ymm10,%%ymm10	*/\n\t"\
		"vaddpd	%%ymm0 ,%%ymm9 ,%%ymm9 					\n\t	vaddpd	%%ymm5 ,%%ymm15,%%ymm15	\n\t"\
		"vaddpd	%%ymm7 ,%%ymm12,%%ymm12					\n\t	vaddpd	%%ymm4 ,%%ymm14,%%ymm14	\n\t"\
		"vaddpd	%%ymm1 ,%%ymm8 ,%%ymm8 					\n\t	vaddpd	%%ymm3 ,%%ymm13,%%ymm13	\n\t"\
		"vmovaps	%%ymm6 ,0x100(%%rsi)				\n\t/*	vmovaps	%%ymm2 ,0x140(%%rsi)	*/\n\t"\
		"vmovaps	%%ymm0 ,0x080(%%rsi)				\n\t	vmovaps	%%ymm5 ,0x0c0(%%rsi)	\n\t"\
		"vmovaps	%%ymm7 ,0x120(%%rsi)				\n\t	vmovaps	%%ymm4 ,0x160(%%rsi)	\n\t"\
		"vmovaps	%%ymm1 ,0x1a0(%%rsi)				\n\t	vmovaps	%%ymm3 ,0x1e0(%%rsi)	\n\t"\
		"vmovaps	%%ymm11,     (%%rsi)				\n\t	vmovaps	%%ymm10,0x040(%%rsi)	\n\t"\
		"vmovaps	%%ymm9 ,0x180(%%rsi)				\n\t	vmovaps	%%ymm15,0x1c0(%%rsi)	\n\t"\
		"vmovaps	%%ymm12,0x020(%%rsi)				\n\t	vmovaps	%%ymm14,0x060(%%rsi)	\n\t"\
		"vmovaps	%%ymm8 ,0x0a0(%%rsi)				\n\t	vmovaps	%%ymm13,0x0e0(%%rsi)	\n\t"\
		"\n\t"\
		"/**********************************************************************************/\n\t"\
		"/*...and now do eight radix-4 transforms, including the internal twiddle factors: */\n\t"\
		"/**********************************************************************************/\n\t"\
		"/*...Block 1: t00,t10,t20,t30	*/				\n\t"\
		"movq	%[__add0],%%rax	/* &a[j1] */			\n\t"\
		"movslq	%[__p01],%%rbx							\n\t"\
		"movslq	%[__p02],%%rcx							\n\t"\
		"movslq	%[__p03],%%rdx							\n\t		/*...Block 5: t08,t18,t28,t38	*/	\n\t"\
		"shlq	$3,%%rbx								\n\t		movslq	%[__p04],%%r9	\n\t"\
		"shlq	$3,%%rcx								\n\t		shlq	$3,%%r9			\n\t"\
		"shlq	$3,%%rdx								\n\t		movq	%%r9,%%r10	\n\t	addq	%%rax,%%r10	\n\t"\
		"addq	%%rax,%%rbx								\n\t		movq	%%r9,%%r11	\n\t	addq	%%rbx,%%r11	\n\t"\
		"addq	%%rax,%%rcx								\n\t		movq	%%r9,%%r12	\n\t	addq	%%rcx,%%r12	\n\t"\
		"addq	%%rax,%%rdx								\n\t		movq	%%r9,%%r13	\n\t	addq	%%rdx,%%r13	\n\t"\
		"movq	%[__r00],%%rsi							\n\t		vmovaps	0x800(%%rsi),%%ymm11	/* isrt2 */	\n\t"\
		"vmovaps	     (%%rsi),%%ymm0					\n\t		vmovaps	0x500(%%rsi),%%ymm12	\n\t"\
		"vmovaps	0x400(%%rsi),%%ymm4					\n\t		vmovaps	0x520(%%rsi),%%ymm13	\n\t"\
		"vmovaps	0x020(%%rsi),%%ymm1					\n\t		vmovaps	0x700(%%rsi),%%ymm14	\n\t"\
		"vmovaps	0x420(%%rsi),%%ymm5					\n\t		vmovaps	0x720(%%rsi),%%ymm15	\n\t"\
		"vmovaps	0x200(%%rsi),%%ymm2					\n\t		vmulpd	%%ymm11,%%ymm12,%%ymm12	\n\t"\
		"vmovaps	0x600(%%rsi),%%ymm6					\n\t		vmovaps	0x100(%%rsi),%%ymm8 	\n\t"\
		"vmovaps	0x220(%%rsi),%%ymm3					\n\t		vmulpd	%%ymm11,%%ymm13,%%ymm13	\n\t"\
		"vmovaps	0x620(%%rsi),%%ymm7					\n\t		vmovaps	0x120(%%rsi),%%ymm9 	\n\t"\
		"vsubpd	0x200(%%rsi),%%ymm0,%%ymm0				\n\t		vmulpd	%%ymm11,%%ymm14,%%ymm14	\n\t"\
		"vsubpd	0x600(%%rsi),%%ymm4,%%ymm4				\n\t		vmovaps	0x300(%%rsi),%%ymm10	\n\t"\
		"vsubpd	0x220(%%rsi),%%ymm1,%%ymm1				\n\t		vmulpd	%%ymm11,%%ymm15,%%ymm15	\n\t"\
		"vsubpd	0x620(%%rsi),%%ymm5,%%ymm5				\n\t		vmovaps	0x320(%%rsi),%%ymm11	\n\t"\
		"vaddpd	     (%%rsi),%%ymm2,%%ymm2				\n\t		vsubpd	%%ymm11,%%ymm8 ,%%ymm8 	\n\t"\
		"vaddpd	0x400(%%rsi),%%ymm6,%%ymm6				\n\t		vsubpd	%%ymm13,%%ymm12,%%ymm12	\n\t"\
		"vaddpd	0x020(%%rsi),%%ymm3,%%ymm3				\n\t		vsubpd	%%ymm10,%%ymm9 ,%%ymm9 	\n\t"\
		"vaddpd	0x420(%%rsi),%%ymm7,%%ymm7				\n\t		vsubpd	%%ymm14,%%ymm15,%%ymm15	\n\t"\
		"vsubpd	%%ymm6,%%ymm2,%%ymm2					\n\t		vmulpd	(%%r8) ,%%ymm11,%%ymm11	\n\t"\
		"vsubpd	%%ymm5,%%ymm0,%%ymm0					\n\t		vmulpd	(%%r8) ,%%ymm13,%%ymm13	\n\t"\
		"vsubpd	%%ymm7,%%ymm3,%%ymm3					\n\t		vmulpd	(%%r8) ,%%ymm10,%%ymm10	\n\t"\
		"vsubpd	%%ymm4,%%ymm1,%%ymm1					\n\t		vmulpd	(%%r8) ,%%ymm14,%%ymm14	\n\t"\
		"vmovaps	%%ymm2,     (%%rbx)					\n\t		vaddpd	%%ymm8 ,%%ymm11,%%ymm11	\n\t"\
		"vmovaps	(%%r8),%%ymm2	/* 2.0 */			\n\t		vaddpd	%%ymm12,%%ymm13,%%ymm13	\n\t"\
		"vmulpd	%%ymm2,%%ymm6,%%ymm6					\n\t		vaddpd	%%ymm9 ,%%ymm10,%%ymm10	\n\t"\
		"vmulpd	%%ymm2,%%ymm5,%%ymm5					\n\t		vaddpd	%%ymm15,%%ymm14,%%ymm14	\n\t"\
		"vmulpd	%%ymm2,%%ymm7,%%ymm7					\n\t		vsubpd	%%ymm14,%%ymm12,%%ymm12	\n\t"\
		"vmulpd	%%ymm2,%%ymm4,%%ymm4					\n\t		vsubpd	%%ymm15,%%ymm13,%%ymm13	\n\t"\
		"/*vmovaps	%%ymm2,     (%%rbx)	*/				\n\t		vmulpd	%%ymm2 ,%%ymm14,%%ymm14	\n\t"\
		"vmovaps	%%ymm0,     (%%rcx)					\n\t		vmulpd	%%ymm2 ,%%ymm15,%%ymm15	\n\t"\
		"vmovaps	%%ymm3,0x020(%%rbx)					\n\t		vaddpd	%%ymm12,%%ymm14,%%ymm14	\n\t"\
		"vmovaps	%%ymm1,0x020(%%rdx)					\n\t		vaddpd	%%ymm13,%%ymm15,%%ymm15	\n\t"\
		"vaddpd	(%%rbx),%%ymm6,%%ymm6					\n\t		vsubpd	%%ymm12,%%ymm8 ,%%ymm8 	\n\t"\
		"vaddpd	%%ymm0,%%ymm5,%%ymm5					\n\t		vsubpd	%%ymm15,%%ymm11,%%ymm11	\n\t"\
		"vaddpd	%%ymm3,%%ymm7,%%ymm7					\n\t		vsubpd	%%ymm13,%%ymm10,%%ymm10	\n\t"\
		"vaddpd	%%ymm1,%%ymm4,%%ymm4					\n\t		vsubpd	%%ymm14,%%ymm9 ,%%ymm9 	\n\t"\
		"vmovaps	%%ymm6,     (%%rax)					\n\t		vmulpd	%%ymm2 ,%%ymm12,%%ymm12	\n\t"\
		"vmovaps	%%ymm5,     (%%rdx)					\n\t		vmulpd	%%ymm2 ,%%ymm15,%%ymm15	\n\t"\
		"vmovaps	%%ymm7,0x020(%%rax)					\n\t		vmulpd	%%ymm2 ,%%ymm13,%%ymm13	\n\t"\
		"vmovaps	%%ymm4,0x020(%%rcx)					\n\t		vmulpd	%%ymm2 ,%%ymm14,%%ymm14	\n\t"\
		"												\n\t		vmovaps	%%ymm8 ,     (%%r11)	\n\t"\
		"												\n\t		vmovaps	%%ymm11,     (%%r12)	\n\t"\
		"												\n\t		vmovaps	%%ymm10,0x020(%%r11)	\n\t"\
		"												\n\t		vmovaps	%%ymm9 ,0x020(%%r13)	\n\t"\
		"												\n\t		vaddpd	%%ymm8 ,%%ymm12,%%ymm12	\n\t"\
		"												\n\t		vaddpd	%%ymm11,%%ymm15,%%ymm15	\n\t"\
		"												\n\t		vaddpd	%%ymm10,%%ymm13,%%ymm13	\n\t"\
		"												\n\t		vaddpd	%%ymm9 ,%%ymm14,%%ymm14	\n\t"\
		"												\n\t		vmovaps	%%ymm12,     (%%r10)	\n\t"\
		"												\n\t		vmovaps	%%ymm15,     (%%r13)	\n\t"\
		"												\n\t		vmovaps	%%ymm13,0x020(%%r10)	\n\t"\
		"												\n\t		vmovaps	%%ymm14,0x020(%%r12)	\n\t"\
		"\n\t"\
		"/*...Block 3: t04,t14,t24,t34	*/				\n\t"\
		"addq	$0x080,%%rsi	/* r04 */				\n\t"\
		"subq	%%rax,%%rbx	/* p01 << 3 */				\n\t"\
		"subq	%%rax,%%rcx	/* p02 << 3 */				\n\t"\
		"subq	%%rax,%%rdx	/* p03 << 3 */				\n\t		/*...Block 7: t0C,t1C,t2C,t3C	*/\n\t"\
		"movslq	%[__p08],%%rdi							\n\t		movslq	%[__p04],%%r9	\n\t"\
		"shlq	$3,%%rdi								\n\t		shlq	$3,%%r9			\n\t"\
		"addq	%%rdi,%%rax	/* &a[j1+p08] */			\n\t		movq	%%r9,%%r10	\n\t	addq	%%rax,%%r10	\n\t"\
		"addq	%%rax,%%rbx								\n\t		movq	%%r9,%%r11	\n\t	addq	%%rbx,%%r11	\n\t"\
		"addq	%%rax,%%rcx								\n\t		movq	%%r9,%%r12	\n\t	addq	%%rcx,%%r12	\n\t"\
		"addq	%%rax,%%rdx								\n\t		movq	%%r9,%%r13	\n\t	addq	%%rdx,%%r13	\n\t"\
		"vmovaps	0x400(%%rsi),%%ymm4					\n\t		vmovaps	0x500(%%rsi),%%ymm12	\n\t"\
		"vmovaps	0x420(%%rsi),%%ymm5					\n\t		vmovaps	0x520(%%rsi),%%ymm13	\n\t"\
		"vmovaps	0x7a0(%%rsi),%%ymm3	/* cc0 */		\n\t"\
		"vmovaps	0x7c0(%%rsi),%%ymm2					\n\t"\
		"vmovaps	%%ymm4,%%ymm6						\n\t		vmovaps	%%ymm12,%%ymm14	\n\t"\
		"vmovaps	%%ymm5,%%ymm7						\n\t		vmovaps	%%ymm13,%%ymm15	\n\t"\
		"vmulpd	%%ymm3,%%ymm4,%%ymm4					\n\t		vmulpd	%%ymm2 ,%%ymm12,%%ymm12	\n\t"\
		"vmulpd	%%ymm3,%%ymm5,%%ymm5					\n\t		vmulpd	%%ymm2 ,%%ymm13,%%ymm13	\n\t"\
		"vmulpd	%%ymm2,%%ymm6,%%ymm6					\n\t		vmulpd	%%ymm3 ,%%ymm14,%%ymm14	\n\t"\
		"vmulpd	%%ymm2,%%ymm7,%%ymm7					\n\t		vmulpd	%%ymm3 ,%%ymm15,%%ymm15	\n\t"\
		"vmovaps	0x600(%%rsi),%%ymm0					\n\t		vmovaps	0x700(%%rsi),%%ymm8 	\n\t"\
		"vmovaps	0x620(%%rsi),%%ymm1					\n\t		vmovaps	0x720(%%rsi),%%ymm9 	\n\t"\
		"vaddpd	%%ymm6,%%ymm5,%%ymm5					\n\t		vaddpd	%%ymm14,%%ymm13,%%ymm13	\n\t"\
		"vmovaps	%%ymm0,%%ymm6						\n\t		vmovaps	%%ymm8 ,%%ymm14	\n\t"\
		"vsubpd	%%ymm7,%%ymm4,%%ymm4					\n\t		vsubpd	%%ymm15,%%ymm12,%%ymm12	\n\t"\
		"vmovaps	%%ymm1,%%ymm7						\n\t		vmovaps	%%ymm9 ,%%ymm15	\n\t"\
		"vmulpd	%%ymm2,%%ymm6,%%ymm6					\n\t		vmulpd	%%ymm3 ,%%ymm14,%%ymm14	\n\t"\
		"vmulpd	%%ymm2,%%ymm7,%%ymm7					\n\t		vmulpd	%%ymm3 ,%%ymm15,%%ymm15	\n\t"\
		"vmulpd	%%ymm3,%%ymm0,%%ymm0					\n\t		vmulpd	%%ymm2 ,%%ymm8 ,%%ymm8 	\n\t"\
		"vmulpd	%%ymm3,%%ymm1,%%ymm1					\n\t		vmulpd	%%ymm2 ,%%ymm9 ,%%ymm9 	\n\t"\
		"vaddpd	%%ymm0,%%ymm7,%%ymm7					\n\t		vaddpd	%%ymm8 ,%%ymm15,%%ymm15	\n\t"\
		"vsubpd	%%ymm1,%%ymm6,%%ymm6					\n\t		vsubpd	%%ymm9 ,%%ymm14,%%ymm14	\n\t"\
		"vmovaps	%%ymm4,%%ymm2						\n\t		vmovaps	%%ymm12,%%ymm10	\n\t"\
		"vmovaps	%%ymm5,%%ymm3						\n\t		vmovaps	%%ymm13,%%ymm11	\n\t"\
		"vsubpd	%%ymm6,%%ymm4,%%ymm4					\n\t		vsubpd	%%ymm14,%%ymm12,%%ymm12	\n\t"\
		"vsubpd	%%ymm7,%%ymm5,%%ymm5					\n\t		vsubpd	%%ymm15,%%ymm13,%%ymm13	\n\t"\
		"vaddpd	%%ymm2,%%ymm6,%%ymm6					\n\t		vaddpd	%%ymm10,%%ymm14,%%ymm14	\n\t"\
		"vaddpd	%%ymm3,%%ymm7,%%ymm7					\n\t		vaddpd	%%ymm11,%%ymm15,%%ymm15	\n\t"\
		"vmovaps	0x200(%%rsi),%%ymm2					\n\t		vmovaps	0x300(%%rsi),%%ymm10	\n\t"\
		"vmovaps	0x220(%%rsi),%%ymm3					\n\t		vmovaps	0x320(%%rsi),%%ymm11	\n\t"\
		"vmovaps	0x780(%%rsi),%%ymm1	/* isrt2 */		\n\t"\
		"vmovaps	%%ymm2,%%ymm0						\n\t		vmovaps	%%ymm10,%%ymm8 	\n\t"\
		"vsubpd	%%ymm3,%%ymm2,%%ymm2					\n\t		vaddpd	%%ymm11,%%ymm10,%%ymm10	\n\t"\
		"vaddpd	%%ymm0,%%ymm3,%%ymm3					\n\t		vsubpd	%%ymm8 ,%%ymm11,%%ymm11	\n\t"\
		"vmulpd	%%ymm1,%%ymm2,%%ymm2					\n\t		vmulpd	%%ymm1 ,%%ymm10,%%ymm10	\n\t"\
		"vmulpd	%%ymm1,%%ymm3,%%ymm3					\n\t		vmulpd	%%ymm1 ,%%ymm11,%%ymm11	\n\t"\
		"vmovaps	     (%%rsi),%%ymm0					\n\t		vmovaps	0x100(%%rsi),%%ymm8 	\n\t"\
		"vmovaps	0x020(%%rsi),%%ymm1					\n\t		vmovaps	0x120(%%rsi),%%ymm9 	\n\t"\
		"vsubpd	%%ymm2,%%ymm0,%%ymm0					\n\t		vsubpd	%%ymm10,%%ymm8 ,%%ymm8 	\n\t"\
		"vsubpd	%%ymm3,%%ymm1,%%ymm1					\n\t		vsubpd	%%ymm11,%%ymm9 ,%%ymm9 	\n\t"\
		"vmulpd	(%%r8),%%ymm2,%%ymm2					\n\t		vmulpd	(%%r8) ,%%ymm10,%%ymm10	\n\t"\
		"vmulpd	(%%r8),%%ymm3,%%ymm3					\n\t		vmulpd	(%%r8) ,%%ymm11,%%ymm11	\n\t"\
		"vaddpd	%%ymm0,%%ymm2,%%ymm2					\n\t		vaddpd	%%ymm8 ,%%ymm10,%%ymm10	\n\t"\
		"vaddpd	%%ymm1,%%ymm3,%%ymm3					\n\t		vaddpd	%%ymm9 ,%%ymm11,%%ymm11	\n\t"\
		"vsubpd	%%ymm6,%%ymm2,%%ymm2					\n\t		vsubpd	%%ymm12,%%ymm8 ,%%ymm8 	\n\t"\
		"vsubpd	%%ymm5,%%ymm0,%%ymm0					\n\t		vsubpd	%%ymm15,%%ymm10,%%ymm10	\n\t"\
		"vsubpd	%%ymm7,%%ymm3,%%ymm3					\n\t		vsubpd	%%ymm13,%%ymm9 ,%%ymm9 	\n\t"\
		"vsubpd	%%ymm4,%%ymm1,%%ymm1					\n\t		vsubpd	%%ymm14,%%ymm11,%%ymm11	\n\t"\
		"vmulpd	(%%r8),%%ymm6,%%ymm6					\n\t		vmulpd	(%%r8) ,%%ymm12,%%ymm12	\n\t"\
		"vmulpd	(%%r8),%%ymm5,%%ymm5					\n\t		vmulpd	(%%r8) ,%%ymm15,%%ymm15	\n\t"\
		"vmulpd	(%%r8),%%ymm7,%%ymm7					\n\t		vmulpd	(%%r8) ,%%ymm13,%%ymm13	\n\t"\
		"vmulpd	(%%r8),%%ymm4,%%ymm4					\n\t		vmulpd	(%%r8) ,%%ymm14,%%ymm14	\n\t"\
		"vmovaps	%%ymm2,     (%%rbx)					\n\t		vmovaps	%%ymm8 ,     (%%r11)	\n\t"\
		"vmovaps	%%ymm0,     (%%rcx)					\n\t		vmovaps	%%ymm10,     (%%r12)	\n\t"\
		"vmovaps	%%ymm3,0x020(%%rbx)					\n\t		vmovaps	%%ymm9 ,0x020(%%r11)	\n\t"\
		"vmovaps	%%ymm1,0x020(%%rdx)					\n\t		vmovaps	%%ymm11,0x020(%%r13)	\n\t"\
		"vaddpd	%%ymm2,%%ymm6,%%ymm6					\n\t		vaddpd	%%ymm8 ,%%ymm12,%%ymm12		\n\t"\
		"vaddpd	%%ymm0,%%ymm5,%%ymm5					\n\t		vaddpd	%%ymm10,%%ymm15,%%ymm15		\n\t"\
		"vaddpd	%%ymm3,%%ymm7,%%ymm7					\n\t		vaddpd	%%ymm9 ,%%ymm13,%%ymm13		\n\t"\
		"vaddpd	%%ymm1,%%ymm4,%%ymm4					\n\t		vaddpd	%%ymm11,%%ymm14,%%ymm14		\n\t"\
		"vmovaps	%%ymm6,     (%%rax)					\n\t		vmovaps	%%ymm12,     (%%r10)	\n\t"\
		"vmovaps	%%ymm5,     (%%rdx)					\n\t		vmovaps	%%ymm15,     (%%r13)	\n\t"\
		"vmovaps	%%ymm7,0x020(%%rax)					\n\t		vmovaps	%%ymm13,0x020(%%r10)	\n\t"\
		"vmovaps	%%ymm4,0x020(%%rcx)					\n\t		vmovaps	%%ymm14,0x020(%%r12)	\n\t"\
		"\n\t"\
		"/*...Block 2: t02,t12,t22,t32	*/				\n\t"\
		"subq	$0x040,%%rsi	/* r02 */				\n\t"\
		"movslq	%[__p10],%%rdi							\n\t"\
		"subq	%%rax,%%rbx								\n\t"\
		"subq	%%rax,%%rcx								\n\t"\
		"subq	%%rax,%%rdx								\n\t		/*...Block 6: t0A,t1A,t2A,t3A	*/	\n\t"\
		"movq	%[__add0],%%rax	/* &a[j1] */			\n\t		movslq	%[__p04],%%r9	\n\t"\
		"shlq	$3,%%rdi								\n\t		shlq	$3,%%r9			\n\t"\
		"addq	%%rdi,%%rax	/* &a[j1+p10) */			\n\t		movq	%%r9,%%r10	\n\t	addq	%%rax,%%r10	\n\t"\
		"addq	%%rax,%%rbx								\n\t		movq	%%r9,%%r11	\n\t	addq	%%rbx,%%r11	\n\t"\
		"addq	%%rax,%%rcx								\n\t		movq	%%r9,%%r12	\n\t	addq	%%rcx,%%r12	\n\t"\
		"addq	%%rax,%%rdx								\n\t		movq	%%r9,%%r13	\n\t	addq	%%rdx,%%r13	\n\t"\
		"vmovaps	0x400(%%rsi),%%ymm4					\n\t		vmovaps	0x500(%%rsi),%%ymm12	\n\t"\
		"vmovaps	0x420(%%rsi),%%ymm5					\n\t		vmovaps	0x520(%%rsi),%%ymm13	\n\t"\
		"vmovaps	0x820(%%rsi),%%ymm2	/* cc1 */		\n\t		vmovaps	0x860(%%rsi),%%ymm11	/* cc3 */	\n\t"\
		"vmovaps	0x840(%%rsi),%%ymm3					\n\t		vmovaps	0x880(%%rsi),%%ymm10	\n\t"\
		"vmovaps	%%ymm4,%%ymm6						\n\t		vmovaps	%%ymm12,%%ymm14	\n\t"\
		"vmovaps	%%ymm5,%%ymm7						\n\t		vmovaps	%%ymm13,%%ymm15	\n\t"\
		"vmulpd	%%ymm2,%%ymm4,%%ymm4					\n\t		vmulpd	%%ymm10,%%ymm12,%%ymm12	\n\t"\
		"vmulpd	%%ymm2,%%ymm5,%%ymm5					\n\t		vmulpd	%%ymm10,%%ymm13,%%ymm13	\n\t"\
		"vmulpd	%%ymm3,%%ymm6,%%ymm6					\n\t		vmulpd	%%ymm11,%%ymm14,%%ymm14	\n\t"\
		"vmulpd	%%ymm3,%%ymm7,%%ymm7					\n\t		vmulpd	%%ymm11,%%ymm15,%%ymm15	\n\t"\
		"vmovaps	0x600(%%rsi),%%ymm0					\n\t		vmovaps	0x700(%%rsi),%%ymm8 	\n\t"\
		"vmovaps	0x620(%%rsi),%%ymm1					\n\t		vmovaps	0x720(%%rsi),%%ymm9 	\n\t"\
		"vaddpd	%%ymm6,%%ymm5,%%ymm5					\n\t		vaddpd	%%ymm14,%%ymm13,%%ymm13	\n\t"\
		"vmovaps	%%ymm0,%%ymm6						\n\t		vmovaps	%%ymm8 ,%%ymm14	\n\t"\
		"vsubpd	%%ymm7,%%ymm4,%%ymm4					\n\t		vsubpd	%%ymm15,%%ymm12,%%ymm12	\n\t"\
		"vmovaps	%%ymm1,%%ymm7						\n\t		vmovaps	%%ymm9 ,%%ymm15	\n\t"\
		"vmulpd	%%ymm11,%%ymm6,%%ymm6					\n\t		vmulpd	%%ymm2 ,%%ymm14,%%ymm14	\n\t"\
		"vmulpd	%%ymm11,%%ymm7,%%ymm7					\n\t		vmulpd	%%ymm2 ,%%ymm15,%%ymm15	\n\t"\
		"vmulpd	%%ymm10,%%ymm0,%%ymm0					\n\t		vmulpd	%%ymm3 ,%%ymm8 ,%%ymm8 	\n\t"\
		"vmulpd	%%ymm10,%%ymm1,%%ymm1					\n\t		vmulpd	%%ymm3 ,%%ymm9 ,%%ymm9 	\n\t"\
		"vaddpd	%%ymm0 ,%%ymm7,%%ymm7					\n\t		vsubpd	%%ymm8 ,%%ymm15,%%ymm15	\n\t"\
		"vsubpd	%%ymm1 ,%%ymm6,%%ymm6					\n\t		vaddpd	%%ymm9 ,%%ymm14,%%ymm14	\n\t"\
		"vmovaps	%%ymm4,%%ymm2						\n\t		vmovaps	%%ymm12,%%ymm10	\n\t"\
		"vmovaps	%%ymm5,%%ymm3						\n\t		vmovaps	%%ymm13,%%ymm11	\n\t"\
		"vsubpd	%%ymm6,%%ymm4,%%ymm4					\n\t		vsubpd	%%ymm14,%%ymm12,%%ymm12	\n\t"\
		"vsubpd	%%ymm7,%%ymm5,%%ymm5					\n\t		vsubpd	%%ymm15,%%ymm13,%%ymm13	\n\t"\
		"vaddpd	%%ymm2,%%ymm6,%%ymm6					\n\t		vaddpd	%%ymm10,%%ymm14,%%ymm14	\n\t"\
		"vaddpd	%%ymm3,%%ymm7,%%ymm7					\n\t		vaddpd	%%ymm11,%%ymm15,%%ymm15	\n\t"\
		"vmovaps	0x200(%%rsi),%%ymm1					\n\t		vmovaps	0x300(%%rsi),%%ymm9 	\n\t"\
		"vmovaps	0x220(%%rsi),%%ymm3					\n\t		vmovaps	0x320(%%rsi),%%ymm11	\n\t"\
		"vmovaps	0x800(%%rsi),%%ymm0	/* ss0 */		\n\t		vmovaps	0x7e0(%%rsi),%%ymm8 		/* cc0 */	\n\t"\
		"vmovaps	%%ymm1,%%ymm2						\n\t		vmovaps	%%ymm9 ,%%ymm10	\n\t"\
		"vmulpd	      %%ymm0,%%ymm1,%%ymm1				\n\t		vmulpd	     %%ymm8 ,%%ymm9 ,%%ymm9 	\n\t"\
		"vmulpd	      %%ymm3,%%ymm0,%%ymm0				\n\t		vmulpd	     %%ymm11,%%ymm8 ,%%ymm8 	\n\t"\
		"vmulpd	0x7e0(%%rsi),%%ymm2,%%ymm2	/* cc0 */	\n\t		vmulpd	0x800(%%rsi),%%ymm10,%%ymm10	/* ss0 */	\n\t"\
		"vmulpd	0x7e0(%%rsi),%%ymm3,%%ymm3				\n\t		vmulpd	0x800(%%rsi),%%ymm11,%%ymm11	\n\t"\
		"vsubpd	      %%ymm0,%%ymm2,%%ymm2				\n\t		vaddpd	     %%ymm8 ,%%ymm10,%%ymm10	\n\t"\
		"vaddpd	      %%ymm1,%%ymm3,%%ymm3				\n\t		vsubpd	     %%ymm9 ,%%ymm11,%%ymm11	\n\t"\
		"vmovaps	     (%%rsi),%%ymm0					\n\t		vmovaps	0x100(%%rsi),%%ymm8 	\n\t"\
		"vmovaps	0x020(%%rsi),%%ymm1					\n\t		vmovaps	0x120(%%rsi),%%ymm9 	\n\t"\
		"vsubpd	%%ymm2,%%ymm0,%%ymm0					\n\t		vsubpd	%%ymm10,%%ymm8 ,%%ymm8 	\n\t"\
		"vsubpd	%%ymm3,%%ymm1,%%ymm1					\n\t		vsubpd	%%ymm11,%%ymm9 ,%%ymm9 	\n\t"\
		"vmulpd	(%%r8),%%ymm2,%%ymm2					\n\t		vmulpd	(%%r8) ,%%ymm10,%%ymm10	\n\t"\
		"vmulpd	(%%r8),%%ymm3,%%ymm3					\n\t		vmulpd	(%%r8) ,%%ymm11,%%ymm11	\n\t"\
		"vaddpd	%%ymm0,%%ymm2,%%ymm2					\n\t		vaddpd	%%ymm8 ,%%ymm10,%%ymm10	\n\t"\
		"vaddpd	%%ymm1,%%ymm3,%%ymm3					\n\t		vaddpd	%%ymm9 ,%%ymm11,%%ymm11	\n\t"\
		"vsubpd	%%ymm6,%%ymm2,%%ymm2					\n\t		vsubpd	%%ymm12,%%ymm8 ,%%ymm8 	\n\t"\
		"vsubpd	%%ymm5,%%ymm0,%%ymm0					\n\t		vsubpd	%%ymm15,%%ymm10,%%ymm10	\n\t"\
		"vsubpd	%%ymm7,%%ymm3,%%ymm3					\n\t		vsubpd	%%ymm13,%%ymm9 ,%%ymm9 	\n\t"\
		"vsubpd	%%ymm4,%%ymm1,%%ymm1					\n\t		vsubpd	%%ymm14,%%ymm11,%%ymm11	\n\t"\
		"vmulpd	(%%r8),%%ymm6,%%ymm6					\n\t		vmulpd	(%%r8) ,%%ymm12,%%ymm12	\n\t"\
		"vmulpd	(%%r8),%%ymm5,%%ymm5					\n\t		vmulpd	(%%r8) ,%%ymm15,%%ymm15	\n\t"\
		"vmulpd	(%%r8),%%ymm7,%%ymm7					\n\t		vmulpd	(%%r8) ,%%ymm13,%%ymm13	\n\t"\
		"vmulpd	(%%r8),%%ymm4,%%ymm4					\n\t		vmulpd	(%%r8) ,%%ymm14,%%ymm14	\n\t"\
		"vmovaps	%%ymm2,     (%%rbx)					\n\t		vmovaps	%%ymm8 ,     (%%r11)	\n\t"\
		"vmovaps	%%ymm0,     (%%rcx)					\n\t		vmovaps	%%ymm10,     (%%r12)	\n\t"\
		"vmovaps	%%ymm3,0x020(%%rbx)					\n\t		vmovaps	%%ymm9 ,0x020(%%r11)	\n\t"\
		"vmovaps	%%ymm1,0x020(%%rdx)					\n\t		vmovaps	%%ymm11,0x020(%%r13)	\n\t"\
		"vaddpd	%%ymm2,%%ymm6,%%ymm6					\n\t		vaddpd	%%ymm8 ,%%ymm12,%%ymm12	\n\t"\
		"vaddpd	%%ymm0,%%ymm5,%%ymm5					\n\t		vaddpd	%%ymm10,%%ymm15,%%ymm15	\n\t"\
		"vaddpd	%%ymm3,%%ymm7,%%ymm7					\n\t		vaddpd	%%ymm9 ,%%ymm13,%%ymm13	\n\t"\
		"vaddpd	%%ymm1,%%ymm4,%%ymm4					\n\t		vaddpd	%%ymm11,%%ymm14,%%ymm14	\n\t"\
		"vmovaps	%%ymm6,     (%%rax)					\n\t		vmovaps	%%ymm12,     (%%r10)	\n\t"\
		"vmovaps	%%ymm5,     (%%rdx)					\n\t		vmovaps	%%ymm15,     (%%r13)	\n\t"\
		"vmovaps	%%ymm7,0x020(%%rax)					\n\t		vmovaps	%%ymm13,0x020(%%r10)	\n\t"\
		"vmovaps	%%ymm4,0x020(%%rcx)					\n\t		vmovaps	%%ymm14,0x020(%%r12)	\n\t"\
		"\n\t"\
		"/*...Block 4: t06,t16,t26,t36	*/	\n\t"\
		"addq	$0x080,%%rsi	/* r06 */	\n\t"\
		"movslq	%[__p18],%%rdi	\n\t"\
		"subq	%%rax,%%rbx	\n\t"\
		"subq	%%rax,%%rcx	\n\t"\
		"subq	%%rax,%%rdx								\n\t		/*...Block 8: t0E,t1E,t2E,t3E	*/	\n\t"\
		"movq	%[__add0],%%rax	/* &a[j1] */			\n\t		movslq	%[__p04],%%r9	\n\t"\
		"shlq	$3,%%rdi								\n\t		shlq	$3,%%r9			\n\t"\
		"addq	%%rdi,%%rax	/* &a[j1+p18] */			\n\t		movq	%%r9,%%r10	\n\t	addq	%%rax,%%r10	\n\t"\
		"addq	%%rax,%%rbx								\n\t		movq	%%r9,%%r11	\n\t	addq	%%rbx,%%r11	\n\t"\
		"addq	%%rax,%%rcx								\n\t		movq	%%r9,%%r12	\n\t	addq	%%rcx,%%r12	\n\t"\
		"addq	%%rax,%%rdx								\n\t		movq	%%r9,%%r13	\n\t	addq	%%rdx,%%r13	\n\t"\
		"vmovaps	0x400(%%rsi),%%ymm4					\n\t		vmovaps	0x500(%%rsi),%%ymm12	\n\t"\
		"vmovaps	0x420(%%rsi),%%ymm5					\n\t		vmovaps	0x520(%%rsi),%%ymm13	\n\t"\
		"vmovaps	0x7e0(%%rsi),%%ymm2	/* cc3 */		\n\t		vmovaps	0x7a0(%%rsi),%%ymm11 /* cc1 */	\n\t"\
		"vmovaps	0x800(%%rsi),%%ymm3					\n\t		vmovaps	0x7c0(%%rsi),%%ymm10 \n\t"\
		"vmovaps	%%ymm4,%%ymm6						\n\t		vmovaps	%%ymm12,%%ymm14	\n\t"\
		"vmovaps	%%ymm5,%%ymm7						\n\t		vmovaps	%%ymm13,%%ymm15	\n\t"\
		"vmulpd	%%ymm2,%%ymm4,%%ymm4					\n\t		vmulpd	%%ymm10,%%ymm12,%%ymm12	\n\t"\
		"vmulpd	%%ymm2,%%ymm5,%%ymm5					\n\t		vmulpd	%%ymm10,%%ymm13,%%ymm13	\n\t"\
		"vmulpd	%%ymm3,%%ymm6,%%ymm6					\n\t		vmulpd	%%ymm11,%%ymm14,%%ymm14	\n\t"\
		"vmulpd	%%ymm3,%%ymm7,%%ymm7					\n\t		vmulpd	%%ymm11,%%ymm15,%%ymm15	\n\t"\
		"vmovaps	0x600(%%rsi),%%ymm0					\n\t		vmovaps	0x700(%%rsi),%%ymm8 	\n\t"\
		"vmovaps	0x620(%%rsi),%%ymm1					\n\t		vmovaps	0x720(%%rsi),%%ymm9 	\n\t"\
		"vaddpd	%%ymm6,%%ymm5,%%ymm5					\n\t		vaddpd	%%ymm14,%%ymm13,%%ymm13	\n\t"\
		"vmovaps	%%ymm0,%%ymm6						\n\t		vmovaps	%%ymm8 ,%%ymm14	\n\t"\
		"vsubpd	%%ymm7,%%ymm4,%%ymm4					\n\t		vsubpd	%%ymm15,%%ymm12,%%ymm12	\n\t"\
		"vmovaps	%%ymm1,%%ymm7						\n\t		vmovaps	%%ymm9 ,%%ymm15	\n\t"\
		"vmulpd	%%ymm10,%%ymm6,%%ymm6					\n\t		vmulpd	%%ymm3 ,%%ymm14,%%ymm14	\n\t"\
		"vmulpd	%%ymm10,%%ymm7,%%ymm7					\n\t		vmulpd	%%ymm3 ,%%ymm15,%%ymm15	\n\t"\
		"vmulpd	%%ymm11,%%ymm0,%%ymm0					\n\t		vmulpd	%%ymm2 ,%%ymm8 ,%%ymm8 	\n\t"\
		"vmulpd	%%ymm11,%%ymm1,%%ymm1					\n\t		vmulpd	%%ymm2 ,%%ymm9 ,%%ymm9 	\n\t"\
		"vsubpd	%%ymm0 ,%%ymm7,%%ymm7					\n\t		vaddpd	%%ymm8 ,%%ymm15,%%ymm15	\n\t"\
		"vaddpd	%%ymm1 ,%%ymm6,%%ymm6					\n\t		vsubpd	%%ymm9 ,%%ymm14,%%ymm14	\n\t"\
		"vmovaps	%%ymm4,%%ymm2						\n\t		vmovaps	%%ymm12,%%ymm10	\n\t"\
		"vmovaps	%%ymm5,%%ymm3						\n\t		vmovaps	%%ymm13,%%ymm11	\n\t"\
		"vsubpd	%%ymm6,%%ymm4,%%ymm4					\n\t		vsubpd	%%ymm14,%%ymm12,%%ymm12	\n\t"\
		"vsubpd	%%ymm7,%%ymm5,%%ymm5					\n\t		vsubpd	%%ymm15,%%ymm13,%%ymm13	\n\t"\
		"vaddpd	%%ymm2,%%ymm6,%%ymm6					\n\t		vaddpd	%%ymm10,%%ymm14,%%ymm14	\n\t"\
		"vaddpd	%%ymm3,%%ymm7,%%ymm7					\n\t		vaddpd	%%ymm11,%%ymm15,%%ymm15	\n\t"\
		"vmovaps	0x200(%%rsi),%%ymm1					\n\t		vmovaps	0x300(%%rsi),%%ymm9 	\n\t"\
		"vmovaps	0x220(%%rsi),%%ymm3					\n\t		vmovaps	0x320(%%rsi),%%ymm11	\n\t"\
		"vmovaps	0x760(%%rsi),%%ymm0		/* cc0 */	\n\t		vmovaps	0x780(%%rsi),%%ymm8 	/* ss0 */	\n\t"\
		"vmovaps	%%ymm1,%%ymm2						\n\t		vmovaps	%%ymm9 ,%%ymm10	\n\t"\
		"vmulpd	      %%ymm0,%%ymm1,%%ymm1				\n\t		vmulpd	     %%ymm8 ,%%ymm9 ,%%ymm9 	\n\t"\
		"vmulpd	      %%ymm3,%%ymm0,%%ymm0				\n\t		vmulpd	     %%ymm11,%%ymm8 ,%%ymm8 	\n\t"\
		"vmulpd	0x780(%%rsi),%%ymm2,%%ymm2	/* ss0 */	\n\t		vmulpd	0x760(%%rsi),%%ymm10,%%ymm10	/* cc0 */	\n\t"\
		"vmulpd	0x780(%%rsi),%%ymm3,%%ymm3				\n\t		vmulpd	0x760(%%rsi),%%ymm11,%%ymm11	\n\t"\
		"vsubpd	      %%ymm0,%%ymm2,%%ymm2				\n\t		vaddpd	     %%ymm8 ,%%ymm10,%%ymm10	\n\t"\
		"vaddpd	      %%ymm1,%%ymm3,%%ymm3				\n\t		vsubpd	     %%ymm9 ,%%ymm11,%%ymm11	\n\t"\
		"vmovaps	     (%%rsi),%%ymm0					\n\t		vmovaps	0x100(%%rsi),%%ymm8 	\n\t"\
		"vmovaps	0x020(%%rsi),%%ymm1					\n\t		vmovaps	0x120(%%rsi),%%ymm9 	\n\t"\
		"vsubpd	%%ymm2,%%ymm0,%%ymm0					\n\t		vsubpd	%%ymm10,%%ymm8 ,%%ymm8 	\n\t"\
		"vsubpd	%%ymm3,%%ymm1,%%ymm1					\n\t		vsubpd	%%ymm11,%%ymm9 ,%%ymm9 	\n\t"\
		"vmulpd	(%%r8),%%ymm2,%%ymm2					\n\t		vmulpd	(%%r8) ,%%ymm10,%%ymm10	\n\t"\
		"vmulpd	(%%r8),%%ymm3,%%ymm3					\n\t		vmulpd	(%%r8) ,%%ymm11,%%ymm11	\n\t"\
		"vaddpd	%%ymm0,%%ymm2,%%ymm2					\n\t		vaddpd	%%ymm8 ,%%ymm10,%%ymm10	\n\t"\
		"vaddpd	%%ymm1,%%ymm3,%%ymm3					\n\t		vaddpd	%%ymm9 ,%%ymm11,%%ymm11	\n\t"\
		"vsubpd	%%ymm4,%%ymm2,%%ymm2					\n\t		vsubpd	%%ymm12,%%ymm8 ,%%ymm8 	\n\t"\
		"vsubpd	%%ymm7,%%ymm0,%%ymm0					\n\t		vsubpd	%%ymm15,%%ymm10,%%ymm10	\n\t"\
		"vsubpd	%%ymm5,%%ymm3,%%ymm3					\n\t		vsubpd	%%ymm13,%%ymm9 ,%%ymm9 	\n\t"\
		"vsubpd	%%ymm6,%%ymm1,%%ymm1					\n\t		vsubpd	%%ymm14,%%ymm11,%%ymm11	\n\t"\
		"vmulpd	(%%r8),%%ymm4,%%ymm4					\n\t		vmulpd	(%%r8) ,%%ymm12,%%ymm12	\n\t"\
		"vmulpd	(%%r8),%%ymm7,%%ymm7					\n\t		vmulpd	(%%r8) ,%%ymm15,%%ymm15	\n\t"\
		"vmulpd	(%%r8),%%ymm5,%%ymm5					\n\t		vmulpd	(%%r8) ,%%ymm13,%%ymm13	\n\t"\
		"vmulpd	(%%r8),%%ymm6,%%ymm6					\n\t		vmulpd	(%%r8) ,%%ymm14,%%ymm14	\n\t"\
		"vmovaps	%%ymm2,     (%%rbx)					\n\t		vmovaps	%%ymm8 ,     (%%r11)	\n\t"\
		"vmovaps	%%ymm0,     (%%rcx)					\n\t		vmovaps	%%ymm10,     (%%r12)	\n\t"\
		"vmovaps	%%ymm3,0x020(%%rbx)					\n\t		vmovaps	%%ymm9 ,0x020(%%r11)	\n\t"\
		"vmovaps	%%ymm1,0x020(%%rdx)					\n\t		vmovaps	%%ymm11,0x020(%%r13)	\n\t"\
		"vaddpd	%%ymm2,%%ymm4,%%ymm4					\n\t		vaddpd	%%ymm8 ,%%ymm12,%%ymm12	\n\t"\
		"vaddpd	%%ymm0,%%ymm7,%%ymm7					\n\t		vaddpd	%%ymm10,%%ymm15,%%ymm15	\n\t"\
		"vaddpd	%%ymm3,%%ymm5,%%ymm5					\n\t		vaddpd	%%ymm9 ,%%ymm13,%%ymm13	\n\t"\
		"vaddpd	%%ymm1,%%ymm6,%%ymm6					\n\t		vaddpd	%%ymm11,%%ymm14,%%ymm14	\n\t"\
		"vmovaps	%%ymm4,     (%%rax)					\n\t		vmovaps	%%ymm12,     (%%r10)	\n\t"\
		"vmovaps	%%ymm7,     (%%rdx)					\n\t		vmovaps	%%ymm15,     (%%r13)	\n\t"\
		"vmovaps	%%ymm5,0x020(%%rax)					\n\t		vmovaps	%%ymm13,0x020(%%r10)	\n\t"\
		"vmovaps	%%ymm6,0x020(%%rcx)					\n\t		vmovaps	%%ymm14,0x020(%%r12)	\n\t"\
		"\n\t"\
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
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r8","r9","r10","r11","r12","r13","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"		/* Clobbered registers */\
	);\
	}

	/*
	For GCC-macro version of this, use that isrt2 + 0x020,0x060,0x0a0 = cc0,cc1,cc3,
	and isrt2 + 0x0e0,0x1e0,0x2e0,0x3e0,0x4e0,0x5e0,0x6e0,0x7e0 = c00,04,02,06,01,05,03,07
	in order to reduce number of args to <= the GCC-allowed maximum of 30:
	*/
	#define SSE2_RADIX32_DIT_TWIDDLE(Xadd0,Xp01,Xp02,Xp03,Xp04,Xp05,Xp06,Xp07,Xp08,Xp10,Xp18,Xr00,Xr10,Xr20,Xr30,Xisrt2)\
	{\
	__asm__ volatile (\
	"/**************************...Block 1: ****************************/\n\t"\
		"movq	%[__r00],%%rsi					\n\t"\
		"movq	%[__add0],%%rax					\n\t"\
		"movslq	%[__p01],%%rbx					\n\t"\
		"movslq	%[__p02],%%rcx					\n\t"\
		"movslq	%[__p03],%%rdx					\n\t		movq	%[__isrt2],%%r8	\n\t"\
		"shlq	$3,%%rbx						\n\t		movslq	%[__p04],%%r9	\n\t"\
		"shlq	$3,%%rcx						\n\t		shlq	$3,%%r9		\n\t	addq	$0x900,%%r8	/* two */\n\t"\
		"shlq	$3,%%rdx						\n\t		movq	%%r9,%%r10	\n\t	addq	%%rax,%%r10	\n\t"\
		"addq	%%rax,%%rbx						\n\t		movq	%%r9,%%r11	\n\t	addq	%%rbx,%%r11	\n\t"\
		"addq	%%rax,%%rcx						\n\t		movq	%%r9,%%r12	\n\t	addq	%%rcx,%%r12	\n\t"\
		"addq	%%rax,%%rdx						\n\t		movq	%%r9,%%r13	\n\t	addq	%%rdx,%%r13	\n\t"\
		"/* SSE2_RADIX4_DIT_0TWIDDLE_B(r00) */	\n\t	/* SSE2_RADIX4_DIT_0TWIDDLE_2NDOFTWO_B(r08) */\n\t"\
		"vmovaps	     (%%rax),%%ymm2			\n\t		vmovaps	     (%%r10),%%ymm10		\n\t"\
		"vmovaps	     (%%rcx),%%ymm6			\n\t		vmovaps	     (%%r12),%%ymm14		\n\t"\
		"vmovaps	0x020(%%rax),%%ymm3			\n\t		vmovaps	0x020(%%r10),%%ymm11		\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm7			\n\t		vmovaps	0x020(%%r12),%%ymm15		\n\t"\
		"vmovaps	     (%%rbx),%%ymm0			\n\t		vmovaps	     (%%r11),%%ymm8 		\n\t"\
		"vmovaps	     (%%rdx),%%ymm4			\n\t		vmovaps	     (%%r13),%%ymm12		\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm1			\n\t		vmovaps	0x020(%%r11),%%ymm9 		\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm5			\n\t		vmovaps	0x020(%%r13),%%ymm13		\n\t"\
		"vsubpd	%%ymm0,%%ymm2,%%ymm2			\n\t		vsubpd	%%ymm8 ,%%ymm10,%%ymm10			\n\t"\
		"vsubpd	%%ymm4,%%ymm6,%%ymm6			\n\t		vsubpd	%%ymm12,%%ymm14,%%ymm14			\n\t"\
		"vsubpd	%%ymm1,%%ymm3,%%ymm3			\n\t		vsubpd	%%ymm9 ,%%ymm11,%%ymm11			\n\t"\
		"vsubpd	%%ymm5,%%ymm7,%%ymm7			\n\t		vsubpd	%%ymm13,%%ymm15,%%ymm15			\n\t"\
		"vmulpd	(%%r8),%%ymm0,%%ymm0			\n\t		vmulpd	(%%r8) ,%%ymm8 ,%%ymm8 			\n\t"\
		"vmulpd	(%%r8),%%ymm4,%%ymm4			\n\t		vmulpd	(%%r8) ,%%ymm12,%%ymm12			\n\t"\
		"vmulpd	(%%r8),%%ymm1,%%ymm1			\n\t		vmulpd	(%%r8) ,%%ymm9 ,%%ymm9 			\n\t"\
		"vmulpd	(%%r8),%%ymm5,%%ymm5			\n\t		vmulpd	(%%r8) ,%%ymm13,%%ymm13			\n\t"\
		"vaddpd	%%ymm2,%%ymm0,%%ymm0			\n\t		vaddpd	%%ymm10,%%ymm8 ,%%ymm8 			\n\t"\
		"vaddpd	%%ymm6,%%ymm4,%%ymm4			\n\t		vaddpd	%%ymm14,%%ymm12,%%ymm12			\n\t"\
		"vaddpd	%%ymm3,%%ymm1,%%ymm1			\n\t		vaddpd	%%ymm11,%%ymm9 ,%%ymm9 			\n\t"\
		"vaddpd	%%ymm7,%%ymm5,%%ymm5			\n\t		vaddpd	%%ymm15,%%ymm13,%%ymm13			\n\t"\
		"/* Finish radix-4 butterfly: */		\n\t"\
		"vsubpd	%%ymm4,%%ymm0,%%ymm0			\n\t		vsubpd	%%ymm12,%%ymm8 ,%%ymm8 			\n\t"\
		"vsubpd	%%ymm7,%%ymm2,%%ymm2			\n\t		vsubpd	%%ymm13,%%ymm9 ,%%ymm9 			\n\t"\
		"vsubpd	%%ymm5,%%ymm1,%%ymm1			\n\t"\
		"vsubpd	%%ymm6,%%ymm3,%%ymm3			\n\t"\
		"/*vmovaps	%%ymm0,0x080(%%rsi)*/		\n\t		vsubpd	%%ymm15,%%ymm10,%%ymm10			\n\t"\
		"/*vmovaps	%%ymm2,0x0c0(%%rsi)*/		\n\t		vmulpd	(%%r8) ,%%ymm12,%%ymm12			\n\t"\
		"/*vmovaps	%%ymm1,0x0a0(%%rsi)*/		\n\t		vsubpd	%%ymm14,%%ymm11,%%ymm11			\n\t"\
		"/*vmovaps	%%ymm3,0x060(%%rsi)*/		\n\t		vmulpd	(%%r8) ,%%ymm13,%%ymm13			\n\t"\
		"vmulpd	(%%r8),%%ymm4,%%ymm4			\n\t		vmulpd	(%%r8) ,%%ymm15,%%ymm15			\n\t"\
		"vmulpd	(%%r8),%%ymm7,%%ymm7			\n\t		vaddpd	%%ymm8 ,%%ymm12,%%ymm12			\n\t"\
		"vmulpd	(%%r8),%%ymm5,%%ymm5			\n\t		vmulpd	(%%r8) ,%%ymm14,%%ymm14			\n\t"\
		"vmulpd	(%%r8),%%ymm6,%%ymm6			\n\t		vaddpd	%%ymm9 ,%%ymm13,%%ymm13			\n\t"\
		"vaddpd	%%ymm0,%%ymm4,%%ymm4			\n\t		vaddpd	%%ymm10,%%ymm15,%%ymm15			\n\t"\
		"vaddpd	%%ymm2,%%ymm7,%%ymm7			\n\t"\
		"vaddpd	%%ymm1,%%ymm5,%%ymm5			\n\t		vaddpd	%%ymm11,%%ymm14,%%ymm14			\n\t"\
		"vaddpd	%%ymm3,%%ymm6,%%ymm6			\n\t"\
		"/*vmovaps	%%ymm4,     (%%rsi)*/		\n\t		vsubpd	%%ymm15,%%ymm11,%%ymm11			\n\t"\
		"/*vmovaps	%%ymm7,0x040(%%rsi)*/		\n\t		vsubpd	%%ymm10,%%ymm14,%%ymm14			\n\t"\
		"/*vmovaps	%%ymm5,0x020(%%rsi)*/		\n\t		vmulpd	(%%r8) ,%%ymm15,%%ymm15			\n\t"\
		"/*vmovaps	%%ymm6,0x0e0(%%rsi)*/		\n\t		vmulpd	(%%r8) ,%%ymm10,%%ymm10			\n\t"\
		"vsubpd	%%ymm12,%%ymm4	,%%ymm4			\n\t		vaddpd	%%ymm11,%%ymm15,%%ymm15			\n\t"\
		"vsubpd	%%ymm9 ,%%ymm0	,%%ymm0			\n\t		vaddpd	%%ymm14,%%ymm10,%%ymm10			\n\t"\
		"vsubpd	%%ymm13,%%ymm5	,%%ymm5			\n\t		vmulpd	-0x900(%%r8),%%ymm11,%%ymm11	/* isrt2 */\n\t"\
		"vsubpd	%%ymm8 ,%%ymm1	,%%ymm1			\n\t		vmulpd	-0x900(%%r8),%%ymm14,%%ymm14	\n\t"\
		"vmulpd	(%%r8) ,%%ymm12	,%%ymm12		\n\t		vmulpd	-0x900(%%r8),%%ymm15,%%ymm15	\n\t"\
		"vmulpd	(%%r8) ,%%ymm9	,%%ymm9			\n\t		vmulpd	-0x900(%%r8),%%ymm10,%%ymm10	\n\t"\
		"vmulpd	(%%r8) ,%%ymm13	,%%ymm13		\n\t		/*vmovaps	%%ymm11,0x160(%%rsi)*/\n\t"\
		"vmulpd	(%%r8) ,%%ymm8	,%%ymm8			\n\t		/*vmovaps	%%ymm14,0x1e0(%%rsi)*/\n\t"\
		"vaddpd	%%ymm4 ,%%ymm12	,%%ymm12		\n\t		/*vmovaps	%%ymm15,0x140(%%rsi)*/\n\t"\
		"vaddpd	%%ymm0 ,%%ymm9	,%%ymm9			\n\t		/*vmovaps	%%ymm10,0x1c0(%%rsi)*/\n\t"\
		"/* SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS(r00,r02,r04,r06,r08,r0A,r0C,r0E) */\n\t"\
		"vaddpd	%%ymm5 ,%%ymm13	,%%ymm13		\n\t		vsubpd	%%ymm15,%%ymm7	,%%ymm7			\n\t"\
		"vaddpd	%%ymm1 ,%%ymm8	,%%ymm8			\n\t		vsubpd	%%ymm14,%%ymm2	,%%ymm2			\n\t"\
		"vmovaps	%%ymm4 ,0x100(%%rsi)		\n\t		vsubpd	%%ymm11,%%ymm3	,%%ymm3			\n\t"\
		"vmovaps	%%ymm0 ,0x180(%%rsi)		\n\t		vsubpd	%%ymm10,%%ymm6	,%%ymm6			\n\t"\
		"vmovaps	%%ymm5 ,0x120(%%rsi)		\n\t		vmulpd	(%%r8) ,%%ymm15	,%%ymm15		\n\t"\
		"vmovaps	%%ymm1 ,0x0a0(%%rsi)		\n\t		vmulpd	(%%r8) ,%%ymm14	,%%ymm14		\n\t"\
		"vmovaps	%%ymm12,     (%%rsi)		\n\t		vmulpd	(%%r8) ,%%ymm11	,%%ymm11		\n\t"\
		"vmovaps	%%ymm9 ,0x080(%%rsi)		\n\t		vmulpd	(%%r8) ,%%ymm10	,%%ymm10		\n\t"\
		"vmovaps	%%ymm13,0x020(%%rsi)		\n\t		vaddpd	%%ymm7 ,%%ymm15	,%%ymm15		\n\t"\
		"vmovaps	%%ymm8 ,0x1a0(%%rsi)		\n\t		vaddpd	%%ymm2 ,%%ymm14	,%%ymm14		\n\t"\
		"													vaddpd	%%ymm3 ,%%ymm11	,%%ymm11		\n\t"\
		"													vaddpd	%%ymm6 ,%%ymm10	,%%ymm10		\n\t"\
		"													vmovaps	%%ymm7 ,0x140(%%rsi)		\n\t"\
		"													vmovaps	%%ymm2 ,0x1c0(%%rsi)		\n\t"\
		"													vmovaps	%%ymm3 ,0x160(%%rsi)		\n\t"\
		"													vmovaps	%%ymm6 ,0x0e0(%%rsi)		\n\t"\
		"													vmovaps	%%ymm15,0x040(%%rsi)		\n\t"\
		"													vmovaps	%%ymm14,0x0c0(%%rsi)		\n\t"\
		"													vmovaps	%%ymm11,0x060(%%rsi)		\n\t"\
		"													vmovaps	%%ymm10,0x1e0(%%rsi)		\n\t"\
		"\n\t"\
	"/**************************...Block 2: ****************************/\n\t"\
		"addq	$0x200,%%rsi	/* r10 */		\n\t"\
		"movslq	%[__p08],%%rdi					\n\t"\
		"shlq	$3,%%rdi						\n\t"\
		"addq	%%rdi,%%rax	/* add0+p08 */\n\t		movq	%%r9,%%r10	\n\t	addq	%%rax,%%r10	\n\t"\
		"addq	%%rdi,%%rbx				\n\t		movq	%%r9,%%r11	\n\t	addq	%%rbx,%%r11	\n\t"\
		"addq	%%rdi,%%rcx				\n\t		movq	%%r9,%%r12	\n\t	addq	%%rcx,%%r12	\n\t"\
		"addq	%%rdi,%%rdx				\n\t		movq	%%r9,%%r13	\n\t	addq	%%rdx,%%r13	\n\t"\
		"/* SSE2_RADIX4_DIT_0TWIDDLE_B(r10) */	\n\t	/* SSE2_RADIX4_DIT_0TWIDDLE_2NDOFTWO_B(r18) */\n\t"\
		"vmovaps	     (%%rax),%%ymm2			\n\t		vmovaps	     (%%r10),%%ymm10		\n\t"\
		"vmovaps	     (%%rcx),%%ymm6			\n\t		vmovaps	     (%%r12),%%ymm14		\n\t"\
		"vmovaps	0x020(%%rax),%%ymm3			\n\t		vmovaps	0x020(%%r10),%%ymm11		\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm7			\n\t		vmovaps	0x020(%%r12),%%ymm15		\n\t"\
		"vmovaps	     (%%rbx),%%ymm0			\n\t		vmovaps	     (%%r11),%%ymm8 		\n\t"\
		"vmovaps	     (%%rdx),%%ymm4			\n\t		vmovaps	     (%%r13),%%ymm12		\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm1			\n\t		vmovaps	0x020(%%r11),%%ymm9 		\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm5			\n\t		vmovaps	0x020(%%r13),%%ymm13		\n\t"\
		"vsubpd	%%ymm0,%%ymm2,%%ymm2			\n\t		vsubpd	%%ymm8 ,%%ymm10,%%ymm10			\n\t"\
		"vsubpd	%%ymm4,%%ymm6,%%ymm6			\n\t		vsubpd	%%ymm12,%%ymm14,%%ymm14			\n\t"\
		"vsubpd	%%ymm1,%%ymm3,%%ymm3			\n\t		vsubpd	%%ymm9 ,%%ymm11,%%ymm11			\n\t"\
		"vsubpd	%%ymm5,%%ymm7,%%ymm7			\n\t		vsubpd	%%ymm13,%%ymm15,%%ymm15			\n\t"\
		"vmulpd	(%%r8),%%ymm0,%%ymm0			\n\t		vmulpd	(%%r8) ,%%ymm8 ,%%ymm8 			\n\t"\
		"vmulpd	(%%r8),%%ymm4,%%ymm4			\n\t		vmulpd	(%%r8) ,%%ymm12,%%ymm12			\n\t"\
		"vmulpd	(%%r8),%%ymm1,%%ymm1			\n\t		vmulpd	(%%r8) ,%%ymm9 ,%%ymm9 			\n\t"\
		"vmulpd	(%%r8),%%ymm5,%%ymm5			\n\t		vmulpd	(%%r8) ,%%ymm13,%%ymm13			\n\t"\
		"vaddpd	%%ymm2,%%ymm0,%%ymm0			\n\t		vaddpd	%%ymm10,%%ymm8 ,%%ymm8 			\n\t"\
		"vaddpd	%%ymm6,%%ymm4,%%ymm4			\n\t		vaddpd	%%ymm14,%%ymm12,%%ymm12			\n\t"\
		"vaddpd	%%ymm3,%%ymm1,%%ymm1			\n\t		vaddpd	%%ymm11,%%ymm9 ,%%ymm9 			\n\t"\
		"vaddpd	%%ymm7,%%ymm5,%%ymm5			\n\t		vaddpd	%%ymm15,%%ymm13,%%ymm13			\n\t"\
		"/* Finish radix-4 butterfly: */		\n\t"\
		"vsubpd	%%ymm4,%%ymm0,%%ymm0			\n\t		vsubpd	%%ymm12,%%ymm8 ,%%ymm8 			\n\t"\
		"vsubpd	%%ymm7,%%ymm2,%%ymm2			\n\t		vsubpd	%%ymm13,%%ymm9 ,%%ymm9 			\n\t"\
		"vsubpd	%%ymm5,%%ymm1,%%ymm1			\n\t		vsubpd	%%ymm15,%%ymm10,%%ymm10			\n\t"\
		"vsubpd	%%ymm6,%%ymm3,%%ymm3			\n\t		vmulpd	(%%r8) ,%%ymm12,%%ymm12			\n\t"\
		"vmulpd	(%%r8),%%ymm4,%%ymm4			\n\t		vsubpd	%%ymm14,%%ymm11,%%ymm11			\n\t"\
		"vmulpd	(%%r8),%%ymm7,%%ymm7			\n\t		vmulpd	(%%r8) ,%%ymm13,%%ymm13			\n\t"\
		"vmulpd	(%%r8),%%ymm5,%%ymm5			\n\t		vmulpd	(%%r8) ,%%ymm15,%%ymm15			\n\t"\
		"vmulpd	(%%r8),%%ymm6,%%ymm6			\n\t		vaddpd	%%ymm8 ,%%ymm12,%%ymm12			\n\t"\
		"vaddpd	%%ymm0,%%ymm4,%%ymm4			\n\t		vmulpd	(%%r8) ,%%ymm14,%%ymm14			\n\t"\
		"vaddpd	%%ymm2,%%ymm7,%%ymm7			\n\t		vaddpd	%%ymm9 ,%%ymm13,%%ymm13			\n\t"\
		"vaddpd	%%ymm1,%%ymm5,%%ymm5			\n\t		vaddpd	%%ymm10,%%ymm15,%%ymm15			\n\t"\
		"vaddpd	%%ymm3,%%ymm6,%%ymm6			\n\t		vaddpd	%%ymm11,%%ymm14,%%ymm14			\n\t"\
		"/* SSE2_RADIX8_DIT_COMBINE_RAD4_SUB...*/			vsubpd	%%ymm15,%%ymm11,%%ymm11			\n\t"\
		"vsubpd	%%ymm12,%%ymm4	,%%ymm4			\n\t		vsubpd	%%ymm10,%%ymm14,%%ymm14			\n\t"\
		"vsubpd	%%ymm9 ,%%ymm0	,%%ymm0			\n\t		vmulpd	(%%r8) ,%%ymm15,%%ymm15			\n\t"\
		"vsubpd	%%ymm13,%%ymm5	,%%ymm5			\n\t		vmulpd	(%%r8) ,%%ymm10,%%ymm10			\n\t"\
		"vsubpd	%%ymm8 ,%%ymm1	,%%ymm1			\n\t		vaddpd	%%ymm11,%%ymm15,%%ymm15			\n\t"\
		"vmulpd	(%%r8) ,%%ymm12	,%%ymm12		\n\t		vaddpd	%%ymm14,%%ymm10,%%ymm10			\n\t"\
		"vmulpd	(%%r8) ,%%ymm9	,%%ymm9			\n\t		vmulpd	-0x900(%%r8),%%ymm11,%%ymm11			\n\t"\
		"vmulpd	(%%r8) ,%%ymm13	,%%ymm13		\n\t		vmulpd	-0x900(%%r8),%%ymm14,%%ymm14			\n\t"\
		"vmulpd	(%%r8) ,%%ymm8	,%%ymm8			\n\t		vmulpd	-0x900(%%r8),%%ymm15,%%ymm15			\n\t"\
		"vaddpd	%%ymm4 ,%%ymm12	,%%ymm12		\n\t		vmulpd	-0x900(%%r8),%%ymm10,%%ymm10			\n\t"\
		"vaddpd	%%ymm0 ,%%ymm9	,%%ymm9			\n\t	/*...S(r00,r02,r04,r06,r08,r0A,r0C,r0E) */\n\t"\
		"vaddpd	%%ymm5 ,%%ymm13	,%%ymm13		\n\t		vsubpd	%%ymm15,%%ymm7	,%%ymm7			\n\t"\
		"vaddpd	%%ymm1 ,%%ymm8	,%%ymm8			\n\t		vsubpd	%%ymm14,%%ymm2	,%%ymm2			\n\t"\
		"vmovaps	%%ymm4 ,0x100(%%rsi)		\n\t		vsubpd	%%ymm11,%%ymm3	,%%ymm3			\n\t"\
		"vmovaps	%%ymm0 ,0x180(%%rsi)		\n\t		vsubpd	%%ymm10,%%ymm6	,%%ymm6			\n\t"\
		"vmovaps	%%ymm5 ,0x120(%%rsi)		\n\t		vmulpd	(%%r8) ,%%ymm15	,%%ymm15		\n\t"\
		"vmovaps	%%ymm1 ,0x0a0(%%rsi)		\n\t		vmulpd	(%%r8) ,%%ymm14	,%%ymm14		\n\t"\
		"vmovaps	%%ymm12,     (%%rsi)		\n\t		vmulpd	(%%r8) ,%%ymm11	,%%ymm11		\n\t"\
		"vmovaps	%%ymm9 ,0x080(%%rsi)		\n\t		vmulpd	(%%r8) ,%%ymm10	,%%ymm10		\n\t"\
		"vmovaps	%%ymm13,0x020(%%rsi)		\n\t		vaddpd	%%ymm7 ,%%ymm15	,%%ymm15		\n\t"\
		"vmovaps	%%ymm8 ,0x1a0(%%rsi)		\n\t		vaddpd	%%ymm2 ,%%ymm14	,%%ymm14		\n\t"\
		"													vaddpd	%%ymm3 ,%%ymm11	,%%ymm11		\n\t"\
		"													vaddpd	%%ymm6 ,%%ymm10	,%%ymm10		\n\t"\
		"													vmovaps	%%ymm7 ,0x140(%%rsi)		\n\t"\
		"													vmovaps	%%ymm2 ,0x1c0(%%rsi)		\n\t"\
		"													vmovaps	%%ymm3 ,0x160(%%rsi)		\n\t"\
		"													vmovaps	%%ymm6 ,0x0e0(%%rsi)		\n\t"\
		"													vmovaps	%%ymm15,0x040(%%rsi)		\n\t"\
		"													vmovaps	%%ymm14,0x0c0(%%rsi)		\n\t"\
		"													vmovaps	%%ymm11,0x060(%%rsi)		\n\t"\
		"													vmovaps	%%ymm10,0x1e0(%%rsi)		\n\t"\
		"\n\t"\
	"/**************************...Block 3: ****************************/\n\t"\
		"addq	$0x200,%%rsi	/* r20 */	\n\t"\
		"movslq	%[__p10],%%r14			\n\t"\
		"shlq	$3,%%r14				\n\t"\
		"subq	%%rdi,%%r14	/* p10-p8 */\n\t"\
		"addq	%%r14,%%rax	/* add0+p10 */\n\t		movq	%%r9,%%r10	\n\t	addq	%%rax,%%r10	\n\t"\
		"addq	%%r14,%%rbx				\n\t		movq	%%r9,%%r11	\n\t	addq	%%rbx,%%r11	\n\t"\
		"addq	%%r14,%%rcx				\n\t		movq	%%r9,%%r12	\n\t	addq	%%rcx,%%r12	\n\t"\
		"addq	%%r14,%%rdx				\n\t		movq	%%r9,%%r13	\n\t	addq	%%rdx,%%r13	\n\t"\
		"/* SSE2_RADIX4_DIT_0TWIDDLE_B(r20) */	\n\t	/* SSE2_RADIX4_DIT_0TWIDDLE_2NDOFTWO_B(r28) */\n\t"\
		"vmovaps	     (%%rax),%%ymm2			\n\t		vmovaps	     (%%r10),%%ymm10		\n\t"\
		"vmovaps	     (%%rcx),%%ymm6			\n\t		vmovaps	     (%%r12),%%ymm14		\n\t"\
		"vmovaps	0x020(%%rax),%%ymm3			\n\t		vmovaps	0x020(%%r10),%%ymm11		\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm7			\n\t		vmovaps	0x020(%%r12),%%ymm15		\n\t"\
		"vmovaps	     (%%rbx),%%ymm0			\n\t		vmovaps	     (%%r11),%%ymm8 		\n\t"\
		"vmovaps	     (%%rdx),%%ymm4			\n\t		vmovaps	     (%%r13),%%ymm12		\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm1			\n\t		vmovaps	0x020(%%r11),%%ymm9 		\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm5			\n\t		vmovaps	0x020(%%r13),%%ymm13		\n\t"\
		"vsubpd	%%ymm0,%%ymm2,%%ymm2			\n\t		vsubpd	%%ymm8 ,%%ymm10,%%ymm10			\n\t"\
		"vsubpd	%%ymm4,%%ymm6,%%ymm6			\n\t		vsubpd	%%ymm12,%%ymm14,%%ymm14			\n\t"\
		"vsubpd	%%ymm1,%%ymm3,%%ymm3			\n\t		vsubpd	%%ymm9 ,%%ymm11,%%ymm11			\n\t"\
		"vsubpd	%%ymm5,%%ymm7,%%ymm7			\n\t		vsubpd	%%ymm13,%%ymm15,%%ymm15			\n\t"\
		"vmulpd	(%%r8),%%ymm0,%%ymm0			\n\t		vmulpd	(%%r8) ,%%ymm8 ,%%ymm8 			\n\t"\
		"vmulpd	(%%r8),%%ymm4,%%ymm4			\n\t		vmulpd	(%%r8) ,%%ymm12,%%ymm12			\n\t"\
		"vmulpd	(%%r8),%%ymm1,%%ymm1			\n\t		vmulpd	(%%r8) ,%%ymm9 ,%%ymm9 			\n\t"\
		"vmulpd	(%%r8),%%ymm5,%%ymm5			\n\t		vmulpd	(%%r8) ,%%ymm13,%%ymm13			\n\t"\
		"vaddpd	%%ymm2,%%ymm0,%%ymm0			\n\t		vaddpd	%%ymm10,%%ymm8 ,%%ymm8 			\n\t"\
		"vaddpd	%%ymm6,%%ymm4,%%ymm4			\n\t		vaddpd	%%ymm14,%%ymm12,%%ymm12			\n\t"\
		"vaddpd	%%ymm3,%%ymm1,%%ymm1			\n\t		vaddpd	%%ymm11,%%ymm9 ,%%ymm9 			\n\t"\
		"vaddpd	%%ymm7,%%ymm5,%%ymm5			\n\t		vaddpd	%%ymm15,%%ymm13,%%ymm13			\n\t"\
		"/* Finish radix-4 butterfly: */		\n\t"\
		"vsubpd	%%ymm4,%%ymm0,%%ymm0			\n\t		vsubpd	%%ymm12,%%ymm8 ,%%ymm8 			\n\t"\
		"vsubpd	%%ymm7,%%ymm2,%%ymm2			\n\t		vsubpd	%%ymm13,%%ymm9 ,%%ymm9 			\n\t"\
		"vsubpd	%%ymm5,%%ymm1,%%ymm1			\n\t		vsubpd	%%ymm15,%%ymm10,%%ymm10			\n\t"\
		"vsubpd	%%ymm6,%%ymm3,%%ymm3			\n\t		vmulpd	(%%r8) ,%%ymm12,%%ymm12			\n\t"\
		"vmulpd	(%%r8),%%ymm4,%%ymm4			\n\t		vsubpd	%%ymm14,%%ymm11,%%ymm11			\n\t"\
		"vmulpd	(%%r8),%%ymm7,%%ymm7			\n\t		vmulpd	(%%r8) ,%%ymm13,%%ymm13			\n\t"\
		"vmulpd	(%%r8),%%ymm5,%%ymm5			\n\t		vmulpd	(%%r8) ,%%ymm15,%%ymm15			\n\t"\
		"vmulpd	(%%r8),%%ymm6,%%ymm6			\n\t		vaddpd	%%ymm8 ,%%ymm12,%%ymm12			\n\t"\
		"vaddpd	%%ymm0,%%ymm4,%%ymm4			\n\t		vmulpd	(%%r8) ,%%ymm14,%%ymm14			\n\t"\
		"vaddpd	%%ymm2,%%ymm7,%%ymm7			\n\t		vaddpd	%%ymm9 ,%%ymm13,%%ymm13			\n\t"\
		"vaddpd	%%ymm1,%%ymm5,%%ymm5			\n\t		vaddpd	%%ymm10,%%ymm15,%%ymm15			\n\t"\
		"vaddpd	%%ymm3,%%ymm6,%%ymm6			\n\t		vaddpd	%%ymm11,%%ymm14,%%ymm14			\n\t"\
		"/* SSE2_RADIX8_DIT_COMBINE_RAD4_SUB...*/			vsubpd	%%ymm15,%%ymm11,%%ymm11			\n\t"\
		"vsubpd	%%ymm12,%%ymm4	,%%ymm4			\n\t		vsubpd	%%ymm10,%%ymm14,%%ymm14			\n\t"\
		"vsubpd	%%ymm9 ,%%ymm0	,%%ymm0			\n\t		vmulpd	(%%r8) ,%%ymm15,%%ymm15			\n\t"\
		"vsubpd	%%ymm13,%%ymm5	,%%ymm5			\n\t		vmulpd	(%%r8) ,%%ymm10,%%ymm10			\n\t"\
		"vsubpd	%%ymm8 ,%%ymm1	,%%ymm1			\n\t		vaddpd	%%ymm11,%%ymm15,%%ymm15			\n\t"\
		"vmulpd	(%%r8) ,%%ymm12	,%%ymm12		\n\t		vaddpd	%%ymm14,%%ymm10,%%ymm10			\n\t"\
		"vmulpd	(%%r8) ,%%ymm9	,%%ymm9			\n\t		vmulpd	-0x900(%%r8),%%ymm11,%%ymm11			\n\t"\
		"vmulpd	(%%r8) ,%%ymm13	,%%ymm13		\n\t		vmulpd	-0x900(%%r8),%%ymm14,%%ymm14			\n\t"\
		"vmulpd	(%%r8) ,%%ymm8	,%%ymm8			\n\t		vmulpd	-0x900(%%r8),%%ymm15,%%ymm15			\n\t"\
		"vaddpd	%%ymm4 ,%%ymm12	,%%ymm12		\n\t		vmulpd	-0x900(%%r8),%%ymm10,%%ymm10			\n\t"\
		"vaddpd	%%ymm0 ,%%ymm9	,%%ymm9			\n\t	/*...S(r00,r02,r04,r06,r08,r0A,r0C,r0E) */\n\t"\
		"vaddpd	%%ymm5 ,%%ymm13	,%%ymm13		\n\t		vsubpd	%%ymm15,%%ymm7	,%%ymm7			\n\t"\
		"vaddpd	%%ymm1 ,%%ymm8	,%%ymm8			\n\t		vsubpd	%%ymm14,%%ymm2	,%%ymm2			\n\t"\
		"vmovaps	%%ymm4 ,0x100(%%rsi)		\n\t		vsubpd	%%ymm11,%%ymm3	,%%ymm3			\n\t"\
		"vmovaps	%%ymm0 ,0x180(%%rsi)		\n\t		vsubpd	%%ymm10,%%ymm6	,%%ymm6			\n\t"\
		"vmovaps	%%ymm5 ,0x120(%%rsi)		\n\t		vmulpd	(%%r8) ,%%ymm15	,%%ymm15		\n\t"\
		"vmovaps	%%ymm1 ,0x0a0(%%rsi)		\n\t		vmulpd	(%%r8) ,%%ymm14	,%%ymm14		\n\t"\
		"vmovaps	%%ymm12,     (%%rsi)		\n\t		vmulpd	(%%r8) ,%%ymm11	,%%ymm11		\n\t"\
		"vmovaps	%%ymm9 ,0x080(%%rsi)		\n\t		vmulpd	(%%r8) ,%%ymm10	,%%ymm10		\n\t"\
		"vmovaps	%%ymm13,0x020(%%rsi)		\n\t		vaddpd	%%ymm7 ,%%ymm15	,%%ymm15		\n\t"\
		"vmovaps	%%ymm8 ,0x1a0(%%rsi)		\n\t		vaddpd	%%ymm2 ,%%ymm14	,%%ymm14		\n\t"\
		"													vaddpd	%%ymm3 ,%%ymm11	,%%ymm11		\n\t"\
		"													vaddpd	%%ymm6 ,%%ymm10	,%%ymm10		\n\t"\
		"													vmovaps	%%ymm7 ,0x140(%%rsi)		\n\t"\
		"													vmovaps	%%ymm2 ,0x1c0(%%rsi)		\n\t"\
		"													vmovaps	%%ymm3 ,0x160(%%rsi)		\n\t"\
		"													vmovaps	%%ymm6 ,0x0e0(%%rsi)		\n\t"\
		"													vmovaps	%%ymm15,0x040(%%rsi)		\n\t"\
		"													vmovaps	%%ymm14,0x0c0(%%rsi)		\n\t"\
		"													vmovaps	%%ymm11,0x060(%%rsi)		\n\t"\
		"													vmovaps	%%ymm10,0x1e0(%%rsi)		\n\t"\
		"\n\t"\
	"/**************************...Block 4: ****************************/\n\t"\
		"addq	$0x200,%%rsi	/* r30 */	\n\t"\
		"movslq	%[__p08],%%rdi			\n\t"\
		"shlq	$3,%%rdi				\n\t"\
		"addq	%%rdi,%%rax	/* add0+p18 */\n\t		movq	%%r9,%%r10	\n\t	addq	%%rax,%%r10	\n\t"\
		"addq	%%rdi,%%rbx				\n\t		movq	%%r9,%%r11	\n\t	addq	%%rbx,%%r11	\n\t"\
		"addq	%%rdi,%%rcx				\n\t		movq	%%r9,%%r12	\n\t	addq	%%rcx,%%r12	\n\t"\
		"addq	%%rdi,%%rdx				\n\t		movq	%%r9,%%r13	\n\t	addq	%%rdx,%%r13	\n\t"\
		"/* SSE2_RADIX4_DIT_0TWIDDLE_B(r30) */\n\t	/* SSE2_RADIX4_DIT_0TWIDDLE_2NDOFTWO_B(r38) */\n\t"\
		"vmovaps	     (%%rax),%%ymm2			\n\t		vmovaps	     (%%r10),%%ymm10		\n\t"\
		"vmovaps	     (%%rcx),%%ymm6			\n\t		vmovaps	     (%%r12),%%ymm14		\n\t"\
		"vmovaps	0x020(%%rax),%%ymm3			\n\t		vmovaps	0x020(%%r10),%%ymm11		\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm7			\n\t		vmovaps	0x020(%%r12),%%ymm15		\n\t"\
		"vmovaps	     (%%rbx),%%ymm0			\n\t		vmovaps	     (%%r11),%%ymm8 		\n\t"\
		"vmovaps	     (%%rdx),%%ymm4			\n\t		vmovaps	     (%%r13),%%ymm12		\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm1			\n\t		vmovaps	0x020(%%r11),%%ymm9 		\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm5			\n\t		vmovaps	0x020(%%r13),%%ymm13		\n\t"\
		"vsubpd	%%ymm0,%%ymm2,%%ymm2			\n\t		vsubpd	%%ymm8 ,%%ymm10,%%ymm10			\n\t"\
		"vsubpd	%%ymm4,%%ymm6,%%ymm6			\n\t		vsubpd	%%ymm12,%%ymm14,%%ymm14			\n\t"\
		"vsubpd	%%ymm1,%%ymm3,%%ymm3			\n\t		vsubpd	%%ymm9 ,%%ymm11,%%ymm11			\n\t"\
		"vsubpd	%%ymm5,%%ymm7,%%ymm7			\n\t		vsubpd	%%ymm13,%%ymm15,%%ymm15			\n\t"\
		"vmulpd	(%%r8),%%ymm0,%%ymm0			\n\t		vmulpd	(%%r8) ,%%ymm8 ,%%ymm8 			\n\t"\
		"vmulpd	(%%r8),%%ymm4,%%ymm4			\n\t		vmulpd	(%%r8) ,%%ymm12,%%ymm12			\n\t"\
		"vmulpd	(%%r8),%%ymm1,%%ymm1			\n\t		vmulpd	(%%r8) ,%%ymm9 ,%%ymm9 			\n\t"\
		"vmulpd	(%%r8),%%ymm5,%%ymm5			\n\t		vmulpd	(%%r8) ,%%ymm13,%%ymm13			\n\t"\
		"vaddpd	%%ymm2,%%ymm0,%%ymm0			\n\t		vaddpd	%%ymm10,%%ymm8 ,%%ymm8 			\n\t"\
		"vaddpd	%%ymm6,%%ymm4,%%ymm4			\n\t		vaddpd	%%ymm14,%%ymm12,%%ymm12			\n\t"\
		"vaddpd	%%ymm3,%%ymm1,%%ymm1			\n\t		vaddpd	%%ymm11,%%ymm9 ,%%ymm9 			\n\t"\
		"vaddpd	%%ymm7,%%ymm5,%%ymm5			\n\t		vaddpd	%%ymm15,%%ymm13,%%ymm13			\n\t"\
		"/* Finish radix-4 butterfly: */		\n\t"\
		"vsubpd	%%ymm4,%%ymm0,%%ymm0			\n\t		vsubpd	%%ymm12,%%ymm8 ,%%ymm8 			\n\t"\
		"vsubpd	%%ymm7,%%ymm2,%%ymm2			\n\t		vsubpd	%%ymm13,%%ymm9 ,%%ymm9 			\n\t"\
		"vsubpd	%%ymm5,%%ymm1,%%ymm1			\n\t		vsubpd	%%ymm15,%%ymm10,%%ymm10			\n\t"\
		"vsubpd	%%ymm6,%%ymm3,%%ymm3			\n\t		vmulpd	(%%r8) ,%%ymm12,%%ymm12			\n\t"\
		"vmulpd	(%%r8),%%ymm4,%%ymm4			\n\t		vsubpd	%%ymm14,%%ymm11,%%ymm11			\n\t"\
		"vmulpd	(%%r8),%%ymm7,%%ymm7			\n\t		vmulpd	(%%r8) ,%%ymm13,%%ymm13			\n\t"\
		"vmulpd	(%%r8),%%ymm5,%%ymm5			\n\t		vmulpd	(%%r8) ,%%ymm15,%%ymm15			\n\t"\
		"vmulpd	(%%r8),%%ymm6,%%ymm6			\n\t		vaddpd	%%ymm8 ,%%ymm12,%%ymm12			\n\t"\
		"vaddpd	%%ymm0,%%ymm4,%%ymm4			\n\t		vmulpd	(%%r8) ,%%ymm14,%%ymm14			\n\t"\
		"vaddpd	%%ymm2,%%ymm7,%%ymm7			\n\t		vaddpd	%%ymm9 ,%%ymm13,%%ymm13			\n\t"\
		"vaddpd	%%ymm1,%%ymm5,%%ymm5			\n\t		vaddpd	%%ymm10,%%ymm15,%%ymm15			\n\t"\
		"vaddpd	%%ymm3,%%ymm6,%%ymm6			\n\t		vaddpd	%%ymm11,%%ymm14,%%ymm14			\n\t"\
		"/* SSE2_RADIX8_DIT_COMBINE_RAD4_SUB...*/			vsubpd	%%ymm15,%%ymm11,%%ymm11			\n\t"\
		"vsubpd	%%ymm12,%%ymm4	,%%ymm4			\n\t		vsubpd	%%ymm10,%%ymm14,%%ymm14			\n\t"\
		"vsubpd	%%ymm9 ,%%ymm0	,%%ymm0			\n\t		vmulpd	(%%r8) ,%%ymm15,%%ymm15			\n\t"\
		"vsubpd	%%ymm13,%%ymm5	,%%ymm5			\n\t		vmulpd	(%%r8) ,%%ymm10,%%ymm10			\n\t"\
		"vsubpd	%%ymm8 ,%%ymm1	,%%ymm1			\n\t		vaddpd	%%ymm11,%%ymm15,%%ymm15			\n\t"\
		"vmulpd	(%%r8) ,%%ymm12	,%%ymm12		\n\t		vaddpd	%%ymm14,%%ymm10,%%ymm10			\n\t"\
		"vmulpd	(%%r8) ,%%ymm9	,%%ymm9			\n\t		vmulpd	-0x900(%%r8),%%ymm11,%%ymm11			\n\t"\
		"vmulpd	(%%r8) ,%%ymm13	,%%ymm13		\n\t		vmulpd	-0x900(%%r8),%%ymm14,%%ymm14			\n\t"\
		"vmulpd	(%%r8) ,%%ymm8	,%%ymm8			\n\t		vmulpd	-0x900(%%r8),%%ymm15,%%ymm15			\n\t"\
		"vaddpd	%%ymm4 ,%%ymm12	,%%ymm12		\n\t		vmulpd	-0x900(%%r8),%%ymm10,%%ymm10			\n\t"\
		"vaddpd	%%ymm0 ,%%ymm9	,%%ymm9			\n\t	/*...S(r00,r02,r04,r06,r08,r0A,r0C,r0E) */\n\t"\
		"vaddpd	%%ymm5 ,%%ymm13	,%%ymm13		\n\t		vsubpd	%%ymm15,%%ymm7	,%%ymm7			\n\t"\
		"vaddpd	%%ymm1 ,%%ymm8	,%%ymm8			\n\t		vsubpd	%%ymm14,%%ymm2	,%%ymm2			\n\t"\
		"vmovaps	%%ymm4 ,0x100(%%rsi)		\n\t		vsubpd	%%ymm11,%%ymm3	,%%ymm3			\n\t"\
		"vmovaps	%%ymm0 ,0x180(%%rsi)		\n\t		vsubpd	%%ymm10,%%ymm6	,%%ymm6			\n\t"\
		"vmovaps	%%ymm5 ,0x120(%%rsi)		\n\t		vmulpd	(%%r8) ,%%ymm15	,%%ymm15		\n\t"\
		"vmovaps	%%ymm1 ,0x0a0(%%rsi)		\n\t		vmulpd	(%%r8) ,%%ymm14	,%%ymm14		\n\t"\
		"vmovaps	%%ymm12,     (%%rsi)		\n\t		vmulpd	(%%r8) ,%%ymm11	,%%ymm11		\n\t"\
		"vmovaps	%%ymm9 ,0x080(%%rsi)		\n\t		vmulpd	(%%r8) ,%%ymm10	,%%ymm10		\n\t"\
		"vmovaps	%%ymm13,0x020(%%rsi)		\n\t		vaddpd	%%ymm7 ,%%ymm15	,%%ymm15		\n\t"\
		"vmovaps	%%ymm8 ,0x1a0(%%rsi)		\n\t		vaddpd	%%ymm2 ,%%ymm14	,%%ymm14		\n\t"\
		"													vaddpd	%%ymm3 ,%%ymm11	,%%ymm11		\n\t"\
		"													vaddpd	%%ymm6 ,%%ymm10	,%%ymm10		\n\t"\
		"													vmovaps	%%ymm7 ,0x140(%%rsi)		\n\t"\
		"													vmovaps	%%ymm2 ,0x1c0(%%rsi)		\n\t"\
		"													vmovaps	%%ymm3 ,0x160(%%rsi)		\n\t"\
		"													vmovaps	%%ymm6 ,0x0e0(%%rsi)		\n\t"\
		"													vmovaps	%%ymm15,0x040(%%rsi)		\n\t"\
		"													vmovaps	%%ymm14,0x0c0(%%rsi)		\n\t"\
		"													vmovaps	%%ymm11,0x060(%%rsi)		\n\t"\
		"													vmovaps	%%ymm10,0x1e0(%%rsi)		\n\t"\
		"\n\t"\
		"/**********************************************************************************/\n\t"\
		"/*...and now do eight radix-4 transforms, including the internal twiddle factors: */\n\t"\
		"/**********************************************************************************/\n\t"\
		"\n\t"\
		"movslq	%[__p10],%%rdi	/* edi will store copy of p10 throughout */\n\t"\
		"/*...Block 1: t00,t10,t20,t30	*/				/*...Block 5: t08,t18,t28,t38*/\n\t"\
		"movq	%[__add0],%%rax							\n\t		movslq	%[__p04],%%rsi		\n\t"\
		"movslq	%[__p08],%%rbx							\n\t		movq	%%rax,%%r10			\n\t"\
		"movq	%[__r00],%%rcx							\n\t		shlq	$3,%%rsi			\n\t"\
		"movq	%[__r10],%%rdx							\n\t		addq	%%rsi	,%%r10	/* add0 = &a[j1+p4] */\n\t"\
		"shlq	$3,%%rbx								\n\t		movq	%[__isrt2],%%rsi	\n\t"\
		"shlq	$3,%%rdi								\n\t		movq	%%rbx,%%r11		/* Need this register-copy before add0 gets added to rbx at left */	\n\t"\
		"addq	%%rax	,%%rbx	/* add1 = add0+p8 */	\n\t		vmovaps	(%%rsi),%%ymm10	/* isrt2 */\n\t"\
		"												\n\t		addq	%%r10	,%%r11	/* add1 = add0+p12 */\n\t"\
		"vmovaps	     (%%rdx),%%ymm2					\n\t		vmovaps	0x500(%%rcx),%%ymm12	\n\t"\
		"vmovaps	0x400(%%rdx),%%ymm4					\n\t		vmovaps	0x520(%%rcx),%%ymm13	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm3					\n\t		vmovaps	0x500(%%rdx),%%ymm14	\n\t"\
		"vmovaps	0x420(%%rdx),%%ymm5					\n\t		vmovaps	0x520(%%rdx),%%ymm15	\n\t"\
		"vmovaps	     (%%rcx),%%ymm0					\n\t		vmulpd	%%ymm10,%%ymm12,%%ymm12		\n\t"\
		"vmovaps	0x400(%%rcx),%%ymm6					\n\t		vmulpd	%%ymm10,%%ymm13,%%ymm13		\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm1					\n\t		vmulpd	%%ymm10,%%ymm14,%%ymm14		\n\t"\
		"vmovaps	0x420(%%rcx),%%ymm7					\n\t		vmulpd	%%ymm10,%%ymm15,%%ymm15		\n\t"\
		"vsubpd	%%ymm2,%%ymm0,%%ymm0					\n\t		vsubpd	%%ymm12,%%ymm13,%%ymm13		\n\t"\
		"vsubpd	%%ymm4,%%ymm6,%%ymm6					\n\t		vmovaps	0x100(%%rcx),%%ymm8 	\n\t"\
		"vsubpd	%%ymm3,%%ymm1,%%ymm1					\n\t		vsubpd	%%ymm15,%%ymm14,%%ymm14		\n\t"\
		"vsubpd	%%ymm5,%%ymm7,%%ymm7					\n\t		vmovaps	0x120(%%rdx),%%ymm10	\n\t"\
		"vmulpd	(%%r8),%%ymm2,%%ymm2					\n\t		vmulpd	(%%r8) ,%%ymm12,%%ymm12		\n\t"\
		"vmulpd	(%%r8),%%ymm4,%%ymm4					\n\t		vmovaps	0x120(%%rcx),%%ymm11	\n\t"\
		"vmulpd	(%%r8),%%ymm3,%%ymm3					\n\t		vmulpd	(%%r8) ,%%ymm15,%%ymm15		\n\t"\
		"vmulpd	(%%r8),%%ymm5,%%ymm5					\n\t		vmovaps	0x100(%%rdx),%%ymm9 	\n\t"\
		"vaddpd	%%ymm0,%%ymm2,%%ymm2					\n\t		vaddpd	%%ymm13,%%ymm12,%%ymm12		\n\t"\
		"vaddpd	%%ymm6,%%ymm4,%%ymm4					\n\t		vaddpd	%%ymm14,%%ymm15,%%ymm15		\n\t"\
		"vaddpd	%%ymm1,%%ymm3,%%ymm3					\n\t		vsubpd	%%ymm14,%%ymm12,%%ymm12		\n\t"\
		"vaddpd	%%ymm7,%%ymm5,%%ymm5					\n\t		vsubpd	%%ymm10,%%ymm8 ,%%ymm8 		\n\t"\
		"/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c00): */\n\t		vsubpd	%%ymm15,%%ymm13,%%ymm13		\n\t"\
		"vsubpd	%%ymm4,%%ymm2,%%ymm2					\n\t		vsubpd	%%ymm9 ,%%ymm11,%%ymm11		\n\t"\
		"vsubpd	%%ymm7,%%ymm0,%%ymm0					\n\t		vmulpd	(%%r8) ,%%ymm14,%%ymm14		\n\t"\
		"vsubpd	%%ymm5,%%ymm3,%%ymm3					\n\t		vmulpd	(%%r8) ,%%ymm10,%%ymm10		\n\t"\
		"vsubpd	%%ymm6,%%ymm1,%%ymm1					\n\t		vmulpd	(%%r8) ,%%ymm15,%%ymm15		\n\t"\
		"vmulpd	(%%r8),%%ymm4,%%ymm4					\n\t		vmulpd	(%%r8) ,%%ymm9 ,%%ymm9 		\n\t"\
		"vmulpd	(%%r8),%%ymm7,%%ymm7					\n\t		vaddpd	%%ymm12,%%ymm14,%%ymm14		\n\t"\
		"vmulpd	(%%r8),%%ymm5,%%ymm5					\n\t		vaddpd	%%ymm8 ,%%ymm10,%%ymm10		\n\t"\
		"vmulpd	(%%r8),%%ymm6,%%ymm6					\n\t		vaddpd	%%ymm13,%%ymm15,%%ymm15		\n\t"\
		"vaddpd	%%ymm2,%%ymm4,%%ymm4					\n\t		vaddpd	%%ymm11,%%ymm9 ,%%ymm9 		\n\t"\
		"vaddpd	%%ymm0,%%ymm7,%%ymm7								/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c04): */\n\t"\
		"vaddpd	%%ymm3,%%ymm5,%%ymm5								/* c04 = c00 + 0x100: */\n\t"\
		"vaddpd	%%ymm1,%%ymm6,%%ymm6					\n\t		vsubpd	%%ymm12,%%ymm10,%%ymm10		\n\t"\
		"addq	$0x0e0,%%rsi	/* c00 */				\n\t		vsubpd	%%ymm15,%%ymm8 ,%%ymm8 		\n\t"\
		"vmovaps	%%ymm2,     (%%rcx)					\n\t		vsubpd	%%ymm13,%%ymm11,%%ymm11		\n\t"\
		"vmovaps	%%ymm0,0x020(%%rdx)					\n\t		vsubpd	%%ymm14,%%ymm9 ,%%ymm9 		\n\t"\
		"vmovaps	%%ymm3,0x020(%%rcx)					\n\t		vmulpd	(%%r8) ,%%ymm12,%%ymm12		\n\t"\
		"vmovaps	%%ymm6,     (%%rdx)					\n\t		vmulpd	(%%r8) ,%%ymm15,%%ymm15		\n\t"\
		"vmovaps	%%ymm4,%%ymm2						\n\t		vmulpd	(%%r8) ,%%ymm13,%%ymm13		\n\t"\
		"vmovaps	%%ymm7,%%ymm0						\n\t		vmulpd	(%%r8) ,%%ymm14,%%ymm14		\n\t"\
		"vmovaps	%%ymm5,%%ymm3						\n\t		vaddpd	%%ymm10,%%ymm12,%%ymm12		\n\t"\
		"vmovaps	%%ymm1,%%ymm6						\n\t		vaddpd	%%ymm8 ,%%ymm15,%%ymm15		\n\t"\
		"vmulpd	     (%%rsi),%%ymm4,%%ymm4				\n\t		vaddpd	%%ymm11,%%ymm13,%%ymm13		\n\t"\
		"vmulpd	0x040(%%rsi),%%ymm7,%%ymm7				\n\t		vaddpd	%%ymm9 ,%%ymm14,%%ymm14		\n\t"\
		"vmulpd	     (%%rsi),%%ymm5,%%ymm5				\n\t		vmovaps	%%ymm10,0x100(%%rcx)	\n\t"\
		"vmulpd	0x040(%%rsi),%%ymm1,%%ymm1				\n\t		vmovaps	%%ymm8 ,0x120(%%rdx)	\n\t"\
		"vmulpd	0x020(%%rsi),%%ymm2,%%ymm2				\n\t		vmovaps	%%ymm11,0x120(%%rcx)	\n\t"\
		"vmulpd	0x060(%%rsi),%%ymm0,%%ymm0				\n\t		vmovaps	%%ymm14,0x100(%%rdx)	\n\t"\
		"vmulpd	0x020(%%rsi),%%ymm3,%%ymm3				\n\t		vmovaps	%%ymm12,%%ymm10		\n\t"\
		"vmulpd	0x060(%%rsi),%%ymm6,%%ymm6				\n\t		vmovaps	%%ymm15,%%ymm8 		\n\t"\
		"vsubpd	      %%ymm2,%%ymm5,%%ymm5				\n\t		vmovaps	%%ymm13,%%ymm11		\n\t"\
		"vsubpd	      %%ymm0,%%ymm1,%%ymm1				\n\t		vmovaps	%%ymm9 ,%%ymm14		\n\t"\
		"vaddpd	      %%ymm3,%%ymm4,%%ymm4				\n\t		vmulpd	0x100(%%rsi),%%ymm12,%%ymm12	\n\t"\
		"vaddpd	      %%ymm6,%%ymm7,%%ymm7				\n\t		vmulpd	0x140(%%rsi),%%ymm15,%%ymm15	\n\t"\
		"vmovaps	%%ymm5,0x020(%%rax)					\n\t		vmulpd	0x100(%%rsi),%%ymm13,%%ymm13	\n\t"\
		"vmovaps	%%ymm1,0x020(%%rbx)					\n\t		vmulpd	0x140(%%rsi),%%ymm9 ,%%ymm9 	\n\t"\
		"vmovaps	%%ymm4,     (%%rax)					\n\t		vmulpd	0x120(%%rsi),%%ymm10,%%ymm10	\n\t"\
		"vmovaps	%%ymm7,     (%%rbx)					\n\t		vmulpd	0x160(%%rsi),%%ymm8 ,%%ymm8 	\n\t"\
		"addq	%%rdi,%%rax								\n\t		vmulpd	0x120(%%rsi),%%ymm11,%%ymm11	\n\t"\
		"addq	%%rdi,%%rbx								\n\t		vmulpd	0x160(%%rsi),%%ymm14,%%ymm14	\n\t"\
		"addq	$0x080,%%rsi		/* c10 */			/* c14 = c10 + 0x100: */\n\t"\
		"vmovaps	     (%%rcx),%%ymm4					\n\t		vsubpd	%%ymm10,%%ymm13,%%ymm13		\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm0					\n\t		vsubpd	%%ymm8 ,%%ymm9 ,%%ymm9 		\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5					\n\t		vaddpd	%%ymm11,%%ymm12,%%ymm12		\n\t"\
		"vmovaps	     (%%rdx),%%ymm6					\n\t		vaddpd	%%ymm14,%%ymm15,%%ymm15		\n\t"\
		"vmovaps	%%ymm4,%%ymm2						\n\t		vmovaps	%%ymm13,0x020(%%r10)	\n\t"\
		"vmovaps	%%ymm0,%%ymm1						\n\t		vmovaps	%%ymm9 ,0x020(%%r11)	\n\t"\
		"vmovaps	%%ymm5,%%ymm3						\n\t		vmovaps	%%ymm12,     (%%r10)	\n\t"\
		"vmovaps	%%ymm6,%%ymm7						\n\t		vmovaps	%%ymm15,     (%%r11)	\n\t"\
		"vmulpd	     (%%rsi),%%ymm4,%%ymm4				\n\t		addq	%%rdi,%%r10			\n\t"\
		"vmulpd	0x040(%%rsi),%%ymm0,%%ymm0				\n\t		addq	%%rdi,%%r11			\n\t"\
		"vmulpd	     (%%rsi),%%ymm5,%%ymm5				\n\t		vmovaps	0x100(%%rcx),%%ymm12	\n\t"\
		"vmulpd	0x040(%%rsi),%%ymm6,%%ymm6				\n\t		vmovaps	0x120(%%rdx),%%ymm8 	\n\t"\
		"vmulpd	0x020(%%rsi),%%ymm2,%%ymm2				\n\t		vmovaps	0x120(%%rcx),%%ymm13	\n\t"\
		"vmulpd	0x060(%%rsi),%%ymm1,%%ymm1				\n\t		vmovaps	0x100(%%rdx),%%ymm14	\n\t"\
		"vmulpd	0x020(%%rsi),%%ymm3,%%ymm3				\n\t		vmovaps	%%ymm12,%%ymm10		\n\t"\
		"vmulpd	0x060(%%rsi),%%ymm7,%%ymm7				\n\t		vmovaps	%%ymm8 ,%%ymm9 		\n\t"\
		"vsubpd	      %%ymm2,%%ymm5,%%ymm5				\n\t		vmovaps	%%ymm13,%%ymm11		\n\t"\
		"vsubpd	      %%ymm1,%%ymm6,%%ymm6				\n\t		vmovaps	%%ymm14,%%ymm15		\n\t"\
		"vaddpd	      %%ymm3,%%ymm4,%%ymm4				\n\t		vmulpd	0x100(%%rsi),%%ymm12,%%ymm12	\n\t"\
		"vaddpd	      %%ymm7,%%ymm0,%%ymm0				\n\t		vmulpd	0x140(%%rsi),%%ymm8 ,%%ymm8 	\n\t"\
		"vmovaps	%%ymm5,0x020(%%rax)					\n\t		vmulpd	0x100(%%rsi),%%ymm13,%%ymm13	\n\t"\
		"vmovaps	%%ymm6,0x020(%%rbx)					\n\t		vmulpd	0x140(%%rsi),%%ymm14,%%ymm14	\n\t"\
		"vmovaps	%%ymm4,     (%%rax)					\n\t		vmulpd	0x120(%%rsi),%%ymm10,%%ymm10	\n\t"\
		"vmovaps	%%ymm0,     (%%rbx)					\n\t		vmulpd	0x160(%%rsi),%%ymm9 ,%%ymm9 	\n\t"\
		"												\n\t		vmulpd	0x120(%%rsi),%%ymm11,%%ymm11	\n\t"\
		"												\n\t		vmulpd	0x160(%%rsi),%%ymm15,%%ymm15	\n\t"\
		"												\n\t		vsubpd	%%ymm10,%%ymm13,%%ymm13		\n\t"\
		"												\n\t		vsubpd	%%ymm9 ,%%ymm14,%%ymm14		\n\t"\
		"												\n\t		vaddpd	%%ymm11,%%ymm12,%%ymm12		\n\t"\
		"												\n\t		vaddpd	%%ymm15,%%ymm8 ,%%ymm8 		\n\t"\
		"												\n\t		vmovaps	%%ymm13,0x020(%%r10)	\n\t"\
		"												\n\t		vmovaps	%%ymm14,0x020(%%r11)	\n\t"\
		"												\n\t		vmovaps	%%ymm12,     (%%r10)	\n\t"\
		"												\n\t		vmovaps	%%ymm8 ,     (%%r11)	\n\t"\
		"\n\t"\
		"/*...Block 2: t02,t12,t22,t32	*/				/*...Block 6: t0A,t1A,t2A,t3A*/\n\t"\
		"movq	%[__add0],%%rax							\n\t		movq	%%rax,%%r10			\n\t"\
		"movslq	%[__p01],%%rsi							\n\t		movslq	%[__p05],%%r9		\n\t"\
		"movslq	%[__p08],%%rbx							\n\t		shlq	$3,%%r9			\n\t"\
		"shlq	$3,%%rsi								\n\t		addq	%%r9 	,%%r10	/* add0 = &a[j1+p5] */\n\t"\
		"shlq	$3,%%rbx								\n\t		movq	%%r10,%%r11			\n\t"\
		"addq	%%rsi	,%%rax	/* add0 = &a[j1+p1] */	\n\t		addq	%%rbx	,%%r11	/* add1 = add0+p13 */\n\t"\
		"addq	%%rax	,%%rbx	/* add1 = add0+p8 */	\n\t"\
		"movq	%[__isrt2],%%rsi						\n\t"\
		"addq	$0x440,%%rcx	/* r22 */				\n\t"\
		"addq	$0x440,%%rdx	/* r32 */				\n\t"\
		"addq	$0x060,%%rsi	/* cc1 */				/* cc3 = cc1 + 0x040: */\n\t"\
		"vmovaps	     (%%rcx),%%ymm4					\n\t		vmovaps	0x100(%%rcx),%%ymm12	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5					\n\t		vmovaps	0x120(%%rcx),%%ymm13	\n\t"\
		"vmovaps	     (%%rsi),%%ymm2					\n\t		vmovaps	0x040(%%rsi),%%ymm11	\n\t"\
		"vmovaps	0x020(%%rsi),%%ymm3					\n\t		vmovaps	0x060(%%rsi),%%ymm10	\n\t"\
		"vmovaps	%%ymm4,%%ymm6						\n\t		vmovaps	%%ymm12,%%ymm14		\n\t"\
		"vmovaps	%%ymm5,%%ymm7						\n\t		vmovaps	%%ymm13,%%ymm15		\n\t"\
		"\n\t"\
		"vmulpd	%%ymm2,%%ymm4,%%ymm4					\n\t		vmulpd	%%ymm10,%%ymm12,%%ymm12		\n\t"\
		"vmulpd	%%ymm2,%%ymm5,%%ymm5					\n\t		vmulpd	%%ymm10,%%ymm13,%%ymm13		\n\t"\
		"vmulpd	%%ymm3,%%ymm6,%%ymm6					\n\t		vmulpd	%%ymm11,%%ymm14,%%ymm14		\n\t"\
		"vmulpd	%%ymm3,%%ymm7,%%ymm7					\n\t		vmulpd	%%ymm11,%%ymm15,%%ymm15		\n\t"\
		"vmovaps	     (%%rdx),%%ymm0					\n\t		vmovaps	0x100(%%rdx),%%ymm8 	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm1					\n\t		vmovaps	0x120(%%rdx),%%ymm9 	/* t3B */\n\t"\
		"/* cc3 */										/* cc1 */\n\t"\
		"/*vmovaps	0x040(%%rsi),%%ymm2					\n\t		vmovaps	     (%%rsi),%%ymm10	*/\n\t"\
		"/*vmovaps	0x060(%%rsi),%%ymm3					\n\t		vmovaps	0x020(%%rsi),%%ymm11	*/\n\t"\
		"vsubpd	%%ymm6,%%ymm5,%%ymm5					\n\t		vsubpd	%%ymm14,%%ymm13,%%ymm13		\n\t"\
		"vmovaps	%%ymm0,%%ymm6						\n\t		vmovaps	%%ymm8 ,%%ymm14		\n\t"\
		"vaddpd	%%ymm7,%%ymm4,%%ymm4					\n\t		vaddpd	%%ymm15,%%ymm12,%%ymm12		\n\t"\
		"vmovaps	%%ymm1,%%ymm7						\n\t		vmovaps	%%ymm9 ,%%ymm15		\n\t"\
		"\n\t"\
		"vmulpd	%%ymm11,%%ymm0,%%ymm0					\n\t		vmulpd	%%ymm2 ,%%ymm8 ,%%ymm8 		\n\t"\
		"vmulpd	%%ymm11,%%ymm1,%%ymm1					\n\t		vmulpd	%%ymm2 ,%%ymm9 ,%%ymm9 		\n\t"\
		"vmulpd	%%ymm10,%%ymm6,%%ymm6					\n\t		vmulpd	%%ymm3 ,%%ymm14,%%ymm14		\n\t"\
		"vmulpd	%%ymm10,%%ymm7,%%ymm7					\n\t		vmulpd	%%ymm3 ,%%ymm15,%%ymm15		\n\t"\
		"vsubpd	%%ymm6 ,%%ymm1,%%ymm1					\n\t		vaddpd	%%ymm14,%%ymm9 ,%%ymm9 		\n\t"\
		"vaddpd	%%ymm7 ,%%ymm0,%%ymm0					\n\t		vsubpd	%%ymm15,%%ymm8 ,%%ymm8 		\n\t"\
		"\n\t"\
		"vmovaps	%%ymm5,%%ymm7						\n\t		vmovaps	%%ymm13,%%ymm15		\n\t"\
		"vmovaps	%%ymm4,%%ymm6						\n\t		vmovaps	%%ymm12,%%ymm14		\n\t"\
		"\n\t"\
		"vaddpd	%%ymm0,%%ymm4,%%ymm4					\n\t		vaddpd	%%ymm8 ,%%ymm14,%%ymm14		\n\t"\
		"vaddpd	%%ymm1,%%ymm5,%%ymm5					\n\t		vaddpd	%%ymm9 ,%%ymm15,%%ymm15		\n\t"\
		"vsubpd	%%ymm0,%%ymm6,%%ymm6					\n\t		vsubpd	%%ymm8 ,%%ymm12,%%ymm12		\n\t"\
		"vsubpd	%%ymm1,%%ymm7,%%ymm7					\n\t		vsubpd	%%ymm9 ,%%ymm13,%%ymm13		\n\t"\
		"\n\t"\
		"subq	$0x400,%%rcx	/* r02 */				\n\t"\
		"subq	$0x400,%%rdx	/* r12 */				\n\t"\
		"subq	$0x040,%%rsi	/* cc0 */				\n\t"\
		"vmovaps	     (%%rdx),%%ymm1					\n\t		vmovaps	0x100(%%rdx),%%ymm8 	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm3					\n\t		vmovaps	0x120(%%rdx),%%ymm10	\n\t"\
		"vmovaps	0x020(%%rsi),%%ymm2					\n\t		vmovaps	0x020(%%rsi),%%ymm9 	\n\t"\
		"vmovaps	%%ymm1,%%ymm0						\n\t		vmovaps	%%ymm8 ,%%ymm11		\n\t"\
		"vmulpd	%%ymm2 ,%%ymm1,%%ymm1					\n\t		vmulpd	%%ymm9 ,%%ymm8 ,%%ymm8 		\n\t"\
		"vmulpd	%%ymm3 ,%%ymm2,%%ymm2					\n\t		vmulpd	%%ymm10,%%ymm9 ,%%ymm9 		\n\t"\
		"vmulpd	(%%rsi),%%ymm0,%%ymm0					\n\t		vmulpd	(%%rsi),%%ymm11,%%ymm11		\n\t"\
		"vmulpd	(%%rsi),%%ymm3,%%ymm3					\n\t		vmulpd	(%%rsi),%%ymm10,%%ymm10		\n\t"\
		"vaddpd	%%ymm0 ,%%ymm2,%%ymm2					\n\t		vsubpd	%%ymm10,%%ymm8 ,%%ymm8 		\n\t"\
		"vsubpd	%%ymm1 ,%%ymm3,%%ymm3					\n\t		vaddpd	%%ymm11,%%ymm9 ,%%ymm9 		\n\t"\
		"\n\t"\
		"vmovaps	     (%%rcx),%%ymm0					\n\t		vmovaps	0x100(%%rcx),%%ymm10	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm1					\n\t		vmovaps	0x120(%%rcx),%%ymm11	\n\t"\
		"vsubpd	%%ymm2,%%ymm0,%%ymm0					\n\t		vsubpd	%%ymm8 ,%%ymm10,%%ymm10		\n\t"\
		"vsubpd	%%ymm3,%%ymm1,%%ymm1					\n\t		vsubpd	%%ymm9 ,%%ymm11,%%ymm11		\n\t"\
		"vmulpd	(%%r8),%%ymm2,%%ymm2					\n\t		vmulpd	(%%r8) ,%%ymm8 ,%%ymm8 		\n\t"\
		"vmulpd	(%%r8),%%ymm3,%%ymm3					\n\t		vmulpd	(%%r8) ,%%ymm9 ,%%ymm9 		\n\t"\
		"vaddpd	%%ymm0,%%ymm2,%%ymm2					\n\t		vaddpd	%%ymm10,%%ymm8 ,%%ymm8 		\n\t"\
		"vaddpd	%%ymm1,%%ymm3,%%ymm3					\n\t		vaddpd	%%ymm11,%%ymm9 ,%%ymm9 		\n\t"\
		"\n\t"\
		"/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c01): *//* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c05): */\n\t"\
		"addq	$0x4c0,%%rsi	/* c01 */				/* c05 = c01 + 0x100: */\n\t"\
		"vsubpd	%%ymm4,%%ymm2,%%ymm2					\n\t		vsubpd	%%ymm12,%%ymm10,%%ymm10		\n\t"\
		"vsubpd	%%ymm7,%%ymm0,%%ymm0					\n\t		vsubpd	%%ymm15,%%ymm8 ,%%ymm8 		\n\t"\
		"vsubpd	%%ymm5,%%ymm3,%%ymm3					\n\t		vsubpd	%%ymm13,%%ymm11,%%ymm11		\n\t"\
		"vsubpd	%%ymm6,%%ymm1,%%ymm1					\n\t		vsubpd	%%ymm14,%%ymm9 ,%%ymm9 		\n\t"\
		"vmulpd	(%%r8),%%ymm4,%%ymm4					\n\t		vmulpd	(%%r8) ,%%ymm12,%%ymm12		\n\t"\
		"vmulpd	(%%r8),%%ymm7,%%ymm7					\n\t		vmulpd	(%%r8) ,%%ymm15,%%ymm15		\n\t"\
		"vmulpd	(%%r8),%%ymm5,%%ymm5					\n\t		vmulpd	(%%r8) ,%%ymm13,%%ymm13		\n\t"\
		"vmulpd	(%%r8),%%ymm6,%%ymm6					\n\t		vmulpd	(%%r8) ,%%ymm14,%%ymm14		\n\t"\
		"vaddpd	%%ymm2,%%ymm4,%%ymm4					\n\t		vaddpd	%%ymm10,%%ymm12,%%ymm12		\n\t"\
		"vaddpd	%%ymm0,%%ymm7,%%ymm7					\n\t		vaddpd	%%ymm8 ,%%ymm15,%%ymm15		\n\t"\
		"vaddpd	%%ymm3,%%ymm5,%%ymm5					\n\t		vaddpd	%%ymm11,%%ymm13,%%ymm13		\n\t"\
		"vaddpd	%%ymm1,%%ymm6,%%ymm6					\n\t		vaddpd	%%ymm9 ,%%ymm14,%%ymm14		\n\t"\
		"vmovaps	%%ymm2,     (%%rcx)					\n\t		vmovaps	%%ymm10,0x100(%%rcx)	\n\t"\
		"vmovaps	%%ymm0,0x020(%%rdx)					\n\t		vmovaps	%%ymm8 ,0x120(%%rdx)	\n\t"\
		"vmovaps	%%ymm3,0x020(%%rcx)					\n\t		vmovaps	%%ymm11,0x120(%%rcx)	\n\t"\
		"vmovaps	%%ymm6,     (%%rdx)					\n\t		vmovaps	%%ymm14,0x100(%%rdx)	\n\t"\
		"vmovaps	%%ymm4,%%ymm2						\n\t		vmovaps	%%ymm12,%%ymm10		\n\t"\
		"vmovaps	%%ymm7,%%ymm0						\n\t		vmovaps	%%ymm15,%%ymm8 		\n\t"\
		"vmovaps	%%ymm5,%%ymm3						\n\t		vmovaps	%%ymm13,%%ymm11		\n\t"\
		"vmovaps	%%ymm1,%%ymm6						\n\t		vmovaps	%%ymm9 ,%%ymm14		\n\t"\
		"vmulpd	     (%%rsi),%%ymm4,%%ymm4				\n\t		vmulpd	0x100(%%rsi),%%ymm12,%%ymm12	\n\t"\
		"vmulpd	0x040(%%rsi),%%ymm7,%%ymm7				\n\t		vmulpd	0x140(%%rsi),%%ymm15,%%ymm15	\n\t"\
		"vmulpd	     (%%rsi),%%ymm5,%%ymm5				\n\t		vmulpd	0x100(%%rsi),%%ymm13,%%ymm13	\n\t"\
		"vmulpd	0x040(%%rsi),%%ymm1,%%ymm1				\n\t		vmulpd	0x140(%%rsi),%%ymm9 ,%%ymm9 	\n\t"\
		"vmulpd	0x020(%%rsi),%%ymm2,%%ymm2				\n\t		vmulpd	0x120(%%rsi),%%ymm10,%%ymm10	\n\t"\
		"vmulpd	0x060(%%rsi),%%ymm0,%%ymm0				\n\t		vmulpd	0x160(%%rsi),%%ymm8 ,%%ymm8 	\n\t"\
		"vmulpd	0x020(%%rsi),%%ymm3,%%ymm3				\n\t		vmulpd	0x120(%%rsi),%%ymm11,%%ymm11	\n\t"\
		"vmulpd	0x060(%%rsi),%%ymm6,%%ymm6				\n\t		vmulpd	0x160(%%rsi),%%ymm14,%%ymm14	\n\t"\
		"vsubpd	      %%ymm2,%%ymm5,%%ymm5				\n\t		vsubpd	%%ymm10,%%ymm13,%%ymm13		\n\t"\
		"vsubpd	      %%ymm0,%%ymm1,%%ymm1				\n\t		vsubpd	%%ymm8 ,%%ymm9 ,%%ymm9 		\n\t"\
		"vaddpd	      %%ymm3,%%ymm4,%%ymm4				\n\t		vaddpd	%%ymm11,%%ymm12,%%ymm12		\n\t"\
		"vaddpd	      %%ymm6,%%ymm7,%%ymm7				\n\t		vaddpd	%%ymm14,%%ymm15,%%ymm15		\n\t"\
		"vmovaps	%%ymm5,0x020(%%rax)					\n\t		vmovaps	%%ymm13,0x020(%%r10)	\n\t"\
		"vmovaps	%%ymm1,0x020(%%rbx)					\n\t		vmovaps	%%ymm9 ,0x020(%%r11)	\n\t"\
		"vmovaps	%%ymm4,     (%%rax)					\n\t		vmovaps	%%ymm12,     (%%r10)	\n\t"\
		"vmovaps	%%ymm7,     (%%rbx)					\n\t		vmovaps	%%ymm15,     (%%r11)	\n\t"\
		"addq	%%rdi,%%rax								\n\t		addq	%%rdi,%%r10			\n\t"\
		"addq	%%rdi,%%rbx								\n\t		addq	%%rdi,%%r11			\n\t"\
		"addq	$0x080,%%rsi		/* c11 */			/* c15 = c11 + 0x100: */\n\t"\
		"vmovaps	     (%%rcx),%%ymm4					\n\t		vmovaps	0x100(%%rcx),%%ymm12	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm0					\n\t		vmovaps	0x120(%%rdx),%%ymm8 	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5					\n\t		vmovaps	0x120(%%rcx),%%ymm13	\n\t"\
		"vmovaps	     (%%rdx),%%ymm6					\n\t		vmovaps	0x100(%%rdx),%%ymm14	\n\t"\
		"vmovaps	%%ymm4,%%ymm2						\n\t		vmovaps	%%ymm12,%%ymm10		\n\t"\
		"vmovaps	%%ymm0,%%ymm1						\n\t		vmovaps	%%ymm8 ,%%ymm9 		\n\t"\
		"vmovaps	%%ymm5,%%ymm3						\n\t		vmovaps	%%ymm13,%%ymm11		\n\t"\
		"vmovaps	%%ymm6,%%ymm7						\n\t		vmovaps	%%ymm14,%%ymm15		\n\t"\
		"vmulpd	     (%%rsi),%%ymm4,%%ymm4				\n\t		vmulpd	0x100(%%rsi),%%ymm12,%%ymm12	\n\t"\
		"vmulpd	0x040(%%rsi),%%ymm0,%%ymm0				\n\t		vmulpd	0x140(%%rsi),%%ymm8 ,%%ymm8 	\n\t"\
		"vmulpd	     (%%rsi),%%ymm5,%%ymm5				\n\t		vmulpd	0x100(%%rsi),%%ymm13,%%ymm13	\n\t"\
		"vmulpd	0x040(%%rsi),%%ymm6,%%ymm6				\n\t		vmulpd	0x140(%%rsi),%%ymm14,%%ymm14	\n\t"\
		"vmulpd	0x020(%%rsi),%%ymm2,%%ymm2				\n\t		vmulpd	0x120(%%rsi),%%ymm10,%%ymm10	\n\t"\
		"vmulpd	0x060(%%rsi),%%ymm1,%%ymm1				\n\t		vmulpd	0x160(%%rsi),%%ymm9 ,%%ymm9 	\n\t"\
		"vmulpd	0x020(%%rsi),%%ymm3,%%ymm3				\n\t		vmulpd	0x120(%%rsi),%%ymm11,%%ymm11	\n\t"\
		"vmulpd	0x060(%%rsi),%%ymm7,%%ymm7				\n\t		vmulpd	0x160(%%rsi),%%ymm15,%%ymm15	\n\t"\
		"vsubpd	      %%ymm2,%%ymm5,%%ymm5				\n\t		vsubpd	%%ymm10,%%ymm13,%%ymm13		\n\t"\
		"vsubpd	      %%ymm1,%%ymm6,%%ymm6				\n\t		vsubpd	%%ymm9 ,%%ymm14,%%ymm14		\n\t"\
		"vaddpd	      %%ymm3,%%ymm4,%%ymm4				\n\t		vaddpd	%%ymm11,%%ymm12,%%ymm12		\n\t"\
		"vaddpd	      %%ymm7,%%ymm0,%%ymm0				\n\t		vaddpd	%%ymm15,%%ymm8 ,%%ymm8 		\n\t"\
		"vmovaps	%%ymm5,0x020(%%rax)					\n\t		vmovaps	%%ymm13,0x020(%%r10)	\n\t"\
		"vmovaps	%%ymm6,0x020(%%rbx)					\n\t		vmovaps	%%ymm14,0x020(%%r11)	\n\t"\
		"vmovaps	%%ymm4,     (%%rax)					\n\t		vmovaps	%%ymm12,     (%%r10)	\n\t"\
		"vmovaps	%%ymm0,     (%%rbx)					\n\t		vmovaps	%%ymm8 ,     (%%r11)	\n\t"\
		"\n\t"\
		"/*...Block 3: t04,t14,t24,t34*/				/*...Block 7: t0C,t1C,t2C,t3C*/\n\t"\
		"movq	%[__add0],%%rax							\n\t		movq	%%rax,%%r10			\n\t"\
		"movslq	%[__p02],%%rsi							\n\t		movslq	%[__p06],%%r9		\n\t"\
		"movslq	%[__p08],%%rbx							\n\t		shlq	$3,%%r9			\n\t"\
		"shlq	$3,%%rsi								\n\t		addq	%%r9 	,%%r10	/* add0 = &a[j1+p6] */\n\t"\
		"shlq	$3,%%rbx								\n\t		movq	%%r10,%%r11			\n\t"\
		"addq	%%rsi	,%%rax	/* add0 = &a[j1+p2] */	\n\t		addq	%%rbx	,%%r11	/* add1 = add0+p14 */\n\t"\
		"addq	%%rax	,%%rbx	/* add1 = add0+p8 */	\n\t"\
		"movq	%[__isrt2],%%rsi						\n\t"\
		"addq	$0x440,%%rcx	/* r24 */				\n\t"\
		"addq	$0x440,%%rdx	/* r34 */				\n\t"\
		"addq	$0x020,%%rsi	/* cc0 */				\n\t"\
		"vmovaps	     (%%rcx),%%ymm4					\n\t		vmovaps	0x100(%%rcx),%%ymm12	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5					\n\t		vmovaps	0x120(%%rcx),%%ymm13	\n\t"\
		"vmovaps	     (%%rsi),%%ymm2					\n\t"\
		"vmovaps	0x020(%%rsi),%%ymm3					\n\t"\
		"vmovaps	%%ymm4,%%ymm6						\n\t		vmovaps	%%ymm12,%%ymm14		\n\t"\
		"vmovaps	%%ymm5,%%ymm7						\n\t		vmovaps	%%ymm13,%%ymm15		\n\t"\
		"\n\t"\
		"vmulpd	%%ymm2,%%ymm4,%%ymm4					\n\t		vmulpd	%%ymm3,%%ymm12,%%ymm12	\n\t"\
		"vmulpd	%%ymm2,%%ymm5,%%ymm5					\n\t		vmulpd	%%ymm3,%%ymm13,%%ymm13	\n\t"\
		"vmulpd	%%ymm3,%%ymm6,%%ymm6					\n\t		vmulpd	%%ymm2,%%ymm14,%%ymm14	\n\t"\
		"vmulpd	%%ymm3,%%ymm7,%%ymm7					\n\t		vmulpd	%%ymm2,%%ymm15,%%ymm15	\n\t"\
		"vmovaps	     (%%rdx),%%ymm0					\n\t		vmovaps	0x100(%%rdx),%%ymm8 	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm1					\n\t		vmovaps	0x120(%%rdx),%%ymm9 	\n\t"\
		"vsubpd	%%ymm6,%%ymm5,%%ymm5					\n\t		vsubpd	%%ymm14,%%ymm13,%%ymm13		\n\t"\
		"vmovaps	%%ymm0,%%ymm6						\n\t		vmovaps	%%ymm8,%%ymm14		\n\t"\
		"vaddpd	%%ymm7,%%ymm4,%%ymm4					\n\t		vaddpd	%%ymm15,%%ymm12,%%ymm12		\n\t"\
		"vmovaps	%%ymm1,%%ymm7						\n\t		vmovaps	%%ymm9,%%ymm15		\n\t"\
		"\n\t"\
		"vmulpd	%%ymm3,%%ymm0,%%ymm0					\n\t		vmulpd	%%ymm2 ,%%ymm8 ,%%ymm8 	\n\t"\
		"vmulpd	%%ymm3,%%ymm1,%%ymm1					\n\t		vmulpd	%%ymm2 ,%%ymm9 ,%%ymm9 	\n\t"\
		"vmulpd	%%ymm2,%%ymm6,%%ymm6					\n\t		vmulpd	%%ymm3 ,%%ymm14,%%ymm14	\n\t"\
		"vmulpd	%%ymm2,%%ymm7,%%ymm7					\n\t		vmulpd	%%ymm3 ,%%ymm15,%%ymm15	\n\t"\
		"vsubpd	%%ymm6,%%ymm1,%%ymm1					\n\t		vsubpd	%%ymm14,%%ymm9 ,%%ymm9 		\n\t"\
		"vaddpd	%%ymm7,%%ymm0,%%ymm0					\n\t		vaddpd	%%ymm15,%%ymm8 ,%%ymm8 		\n\t"\
		"\n\t"\
		"vmovaps	%%ymm5,%%ymm7						\n\t		vmovaps	%%ymm13,%%ymm15		\n\t"\
		"vmovaps	%%ymm4,%%ymm6						\n\t		vmovaps	%%ymm12,%%ymm14		\n\t"\
		"\n\t"\
		"vaddpd	%%ymm0,%%ymm4,%%ymm4					\n\t		vaddpd	%%ymm8,%%ymm14,%%ymm14		\n\t"\
		"vaddpd	%%ymm1,%%ymm5,%%ymm5					\n\t		vaddpd	%%ymm9,%%ymm15,%%ymm15		\n\t"\
		"vsubpd	%%ymm0,%%ymm6,%%ymm6					\n\t		vsubpd	%%ymm8,%%ymm12,%%ymm12		\n\t"\
		"vsubpd	%%ymm1,%%ymm7,%%ymm7					\n\t		vsubpd	%%ymm9,%%ymm13,%%ymm13		\n\t"\
		"\n\t"\
		"subq	$0x400,%%rcx							\n\t"\
		"subq	$0x400,%%rdx							\n\t"\
		"subq	$0x020,%%rsi	/* isrt2 */				\n\t"\
		"vmovaps	     (%%rdx),%%ymm2					\n\t		vmovaps	0x100(%%rdx),%%ymm8 	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm3					\n\t		vmovaps	0x120(%%rdx),%%ymm9 	\n\t"\
		"vmovaps	(%%rsi),%%ymm1	/* isrt2 */	\n\t"\
		"vmovaps	%%ymm3,%%ymm0						\n\t		vmovaps	%%ymm8 ,%%ymm10		\n\t"\
		"vsubpd	%%ymm2,%%ymm3,%%ymm3					\n\t		vsubpd	%%ymm9 ,%%ymm8 	,%%ymm8 	\n\t"\
		"vaddpd	%%ymm0,%%ymm2,%%ymm2					\n\t		vaddpd	%%ymm10,%%ymm9 	,%%ymm9 	\n\t"\
		"vmulpd	%%ymm1,%%ymm2,%%ymm2					\n\t		vmulpd	%%ymm1 ,%%ymm8 	,%%ymm8 \n\t"\
		"vmulpd	%%ymm1,%%ymm3,%%ymm3					\n\t		vmulpd	%%ymm1 ,%%ymm9 	,%%ymm9 \n\t"\
		"\n\t"\
		"vmovaps	     (%%rcx),%%ymm0					\n\t		vmovaps	0x100(%%rcx),%%ymm10	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm1					\n\t		vmovaps	0x120(%%rcx),%%ymm11	\n\t"\
		"vsubpd	%%ymm2,%%ymm0,%%ymm0					\n\t		vsubpd	%%ymm8 ,%%ymm10,%%ymm10		\n\t"\
		"vsubpd	%%ymm3,%%ymm1,%%ymm1					\n\t		vsubpd	%%ymm9 ,%%ymm11,%%ymm11		\n\t"\
		"vmulpd	(%%r8),%%ymm2,%%ymm2					\n\t		vmulpd	(%%r8) ,%%ymm8 ,%%ymm8 		\n\t"\
		"vmulpd	(%%r8),%%ymm3,%%ymm3					\n\t		vmulpd	(%%r8) ,%%ymm9 ,%%ymm9 		\n\t"\
		"vaddpd	%%ymm0,%%ymm2,%%ymm2					\n\t		vaddpd	%%ymm10,%%ymm8 ,%%ymm8 		\n\t"\
		"vaddpd	%%ymm1,%%ymm3,%%ymm3					\n\t		vaddpd	%%ymm11,%%ymm9 ,%%ymm9 		\n\t"\
		"\n\t"\
		"/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c02): */\n\t		/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c06): */\n\t"\
		"addq	$0x2e0,%%rsi	/* c02 */				/* c06 = c02 + 0x100: */\n\t"\
		"vsubpd	%%ymm4,%%ymm2,%%ymm2					\n\t		vsubpd	%%ymm12,%%ymm10,%%ymm10		\n\t"\
		"vsubpd	%%ymm7,%%ymm0,%%ymm0					\n\t		vsubpd	%%ymm15,%%ymm8 ,%%ymm8 		\n\t"\
		"vsubpd	%%ymm5,%%ymm3,%%ymm3					\n\t		vsubpd	%%ymm13,%%ymm11,%%ymm11		\n\t"\
		"vsubpd	%%ymm6,%%ymm1,%%ymm1					\n\t		vsubpd	%%ymm14,%%ymm9 ,%%ymm9 		\n\t"\
		"vmulpd	(%%r8),%%ymm4,%%ymm4					\n\t		vmulpd	(%%r8) ,%%ymm12,%%ymm12		\n\t"\
		"vmulpd	(%%r8),%%ymm7,%%ymm7					\n\t		vmulpd	(%%r8) ,%%ymm15,%%ymm15		\n\t"\
		"vmulpd	(%%r8),%%ymm5,%%ymm5					\n\t		vmulpd	(%%r8) ,%%ymm13,%%ymm13		\n\t"\
		"vmulpd	(%%r8),%%ymm6,%%ymm6					\n\t		vmulpd	(%%r8) ,%%ymm14,%%ymm14		\n\t"\
		"vaddpd	%%ymm2,%%ymm4,%%ymm4					\n\t		vaddpd	%%ymm10,%%ymm12,%%ymm12		\n\t"\
		"vaddpd	%%ymm0,%%ymm7,%%ymm7					\n\t		vaddpd	%%ymm8 ,%%ymm15,%%ymm15		\n\t"\
		"vaddpd	%%ymm3,%%ymm5,%%ymm5					\n\t		vaddpd	%%ymm11,%%ymm13,%%ymm13		\n\t"\
		"vaddpd	%%ymm1,%%ymm6,%%ymm6					\n\t		vaddpd	%%ymm9 ,%%ymm14,%%ymm14		\n\t"\
		"vmovaps	%%ymm2,     (%%rcx)					\n\t		vmovaps	%%ymm10,0x100(%%rcx)	\n\t"\
		"vmovaps	%%ymm0,0x020(%%rdx)					\n\t		vmovaps	%%ymm8 ,0x120(%%rdx)	\n\t"\
		"vmovaps	%%ymm3,0x020(%%rcx)					\n\t		vmovaps	%%ymm11,0x120(%%rcx)	\n\t"\
		"vmovaps	%%ymm6,     (%%rdx)					\n\t		vmovaps	%%ymm14,0x100(%%rdx)	\n\t"\
		"vmovaps	%%ymm4,%%ymm2						\n\t		vmovaps	%%ymm12,%%ymm10		\n\t"\
		"vmovaps	%%ymm7,%%ymm0						\n\t		vmovaps	%%ymm15,%%ymm8 		\n\t"\
		"vmovaps	%%ymm5,%%ymm3						\n\t		vmovaps	%%ymm13,%%ymm11		\n\t"\
		"vmovaps	%%ymm1,%%ymm6						\n\t		vmovaps	%%ymm9 ,%%ymm14		\n\t"\
		"vmulpd	     (%%rsi),%%ymm4,%%ymm4				\n\t		vmulpd	0x100(%%rsi),%%ymm12,%%ymm12	\n\t"\
		"vmulpd	0x040(%%rsi),%%ymm7,%%ymm7				\n\t		vmulpd	0x140(%%rsi),%%ymm15,%%ymm15	\n\t"\
		"vmulpd	     (%%rsi),%%ymm5,%%ymm5				\n\t		vmulpd	0x100(%%rsi),%%ymm13,%%ymm13	\n\t"\
		"vmulpd	0x040(%%rsi),%%ymm1,%%ymm1				\n\t		vmulpd	0x140(%%rsi),%%ymm9 ,%%ymm9 	\n\t"\
		"vmulpd	0x020(%%rsi),%%ymm2,%%ymm2				\n\t		vmulpd	0x120(%%rsi),%%ymm10,%%ymm10	\n\t"\
		"vmulpd	0x060(%%rsi),%%ymm0,%%ymm0				\n\t		vmulpd	0x160(%%rsi),%%ymm8 ,%%ymm8 	\n\t"\
		"vmulpd	0x020(%%rsi),%%ymm3,%%ymm3				\n\t		vmulpd	0x120(%%rsi),%%ymm11,%%ymm11	\n\t"\
		"vmulpd	0x060(%%rsi),%%ymm6,%%ymm6				\n\t		vmulpd	0x160(%%rsi),%%ymm14,%%ymm14	\n\t"\
		"vsubpd	      %%ymm2,%%ymm5,%%ymm5				\n\t		vsubpd	%%ymm10,%%ymm13,%%ymm13		\n\t"\
		"vsubpd	      %%ymm0,%%ymm1,%%ymm1				\n\t		vsubpd	%%ymm8 ,%%ymm9 ,%%ymm9 		\n\t"\
		"vaddpd	      %%ymm3,%%ymm4,%%ymm4				\n\t		vaddpd	%%ymm11,%%ymm12,%%ymm12		\n\t"\
		"vaddpd	      %%ymm6,%%ymm7,%%ymm7				\n\t		vaddpd	%%ymm14,%%ymm15,%%ymm15		\n\t"\
		"vmovaps	%%ymm5,0x020(%%rax)					\n\t		vmovaps	%%ymm13,0x020(%%r10)	\n\t"\
		"vmovaps	%%ymm1,0x020(%%rbx)					\n\t		vmovaps	%%ymm9 ,0x020(%%r11)	\n\t"\
		"vmovaps	%%ymm4,     (%%rax)					\n\t		vmovaps	%%ymm12,     (%%r10)	\n\t"\
		"vmovaps	%%ymm7,     (%%rbx)					\n\t		vmovaps	%%ymm15,     (%%r11)	\n\t"\
		"addq	%%rdi,%%rax								\n\t		addq	%%rdi,%%r10			\n\t"\
		"addq	%%rdi,%%rbx								\n\t		addq	%%rdi,%%r11			\n\t"\
		"addq	$0x080,%%rsi		/* c12 */			/* c16 = c12 + 0x100: */\n\t"\
		"vmovaps	     (%%rcx),%%ymm4					\n\t		vmovaps	0x100(%%rcx),%%ymm12	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm0					\n\t		vmovaps	0x120(%%rdx),%%ymm8 	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5					\n\t		vmovaps	0x120(%%rcx),%%ymm13	\n\t"\
		"vmovaps	     (%%rdx),%%ymm6					\n\t		vmovaps	0x100(%%rdx),%%ymm14	\n\t"\
		"vmovaps	%%ymm4,%%ymm2						\n\t		vmovaps	%%ymm12,%%ymm10		\n\t"\
		"vmovaps	%%ymm0,%%ymm1						\n\t		vmovaps	%%ymm8 ,%%ymm9 		\n\t"\
		"vmovaps	%%ymm5,%%ymm3						\n\t		vmovaps	%%ymm13,%%ymm11		\n\t"\
		"vmovaps	%%ymm6,%%ymm7						\n\t		vmovaps	%%ymm14,%%ymm15		\n\t"\
		"vmulpd	     (%%rsi),%%ymm4,%%ymm4				\n\t		vmulpd	0x100(%%rsi),%%ymm12,%%ymm12	\n\t"\
		"vmulpd	0x040(%%rsi),%%ymm0,%%ymm0				\n\t		vmulpd	0x140(%%rsi),%%ymm8 ,%%ymm8 	\n\t"\
		"vmulpd	     (%%rsi),%%ymm5,%%ymm5				\n\t		vmulpd	0x100(%%rsi),%%ymm13,%%ymm13	\n\t"\
		"vmulpd	0x040(%%rsi),%%ymm6,%%ymm6				\n\t		vmulpd	0x140(%%rsi),%%ymm14,%%ymm14	\n\t"\
		"vmulpd	0x020(%%rsi),%%ymm2,%%ymm2				\n\t		vmulpd	0x120(%%rsi),%%ymm10,%%ymm10	\n\t"\
		"vmulpd	0x060(%%rsi),%%ymm1,%%ymm1				\n\t		vmulpd	0x160(%%rsi),%%ymm9 ,%%ymm9 	\n\t"\
		"vmulpd	0x020(%%rsi),%%ymm3,%%ymm3				\n\t		vmulpd	0x120(%%rsi),%%ymm11,%%ymm11	\n\t"\
		"vmulpd	0x060(%%rsi),%%ymm7,%%ymm7				\n\t		vmulpd	0x160(%%rsi),%%ymm15,%%ymm15	\n\t"\
		"vsubpd	      %%ymm2,%%ymm5,%%ymm5				\n\t		vsubpd	%%ymm10,%%ymm13,%%ymm13		\n\t"\
		"vsubpd	      %%ymm1,%%ymm6,%%ymm6				\n\t		vsubpd	%%ymm9 ,%%ymm14,%%ymm14		\n\t"\
		"vaddpd	      %%ymm3,%%ymm4,%%ymm4				\n\t		vaddpd	%%ymm11,%%ymm12,%%ymm12		\n\t"\
		"vaddpd	      %%ymm7,%%ymm0,%%ymm0				\n\t		vaddpd	%%ymm15,%%ymm8 ,%%ymm8 		\n\t"\
		"vmovaps	%%ymm5,0x020(%%rax)					\n\t		vmovaps	%%ymm13,0x020(%%r10)	\n\t"\
		"vmovaps	%%ymm6,0x020(%%rbx)					\n\t		vmovaps	%%ymm14,0x020(%%r11)	\n\t"\
		"vmovaps	%%ymm4,     (%%rax)					\n\t		vmovaps	%%ymm12,     (%%r10)	\n\t"\
		"vmovaps	%%ymm0,     (%%rbx)					\n\t		vmovaps	%%ymm8 ,     (%%r11)	\n\t"\
		"\n\t"\
		"/*...Block 4: t06,t16,t26,t36*/				\n\t		/*...Block 8: t0E,t1E,t2E,t3E*/\n\t"\
		"movq	%[__add0],%%rax							\n\t		movq	%%rax,%%r10			\n\t"\
		"movslq	%[__p03],%%rsi							\n\t		movslq	%[__p07],%%r9		\n\t"\
		"movslq	%[__p08],%%rbx							\n\t		shlq	$3,%%r9			\n\t"\
		"shlq	$3,%%rsi								\n\t		addq	%%r9 	,%%r10	/* add0 = &a[j1+p6] */\n\t"\
		"shlq	$3,%%rbx								\n\t		movq	%%r10,%%r11			\n\t"\
		"addq	%%rsi	,%%rax	/* add0 = &a[j1+p3] */	\n\t		addq	%%rbx	,%%r11	/* add1 = add0+p14 */\n\t"\
		"addq	%%rax	,%%rbx	/* add1 = add0+p8 */	\n\t"\
		"movq	%[__isrt2],%%rsi						\n\t"\
		"addq	$0x440,%%rcx	/* r26 */				\n\t"\
		"addq	$0x440,%%rdx	/* r36 */				\n\t"\
		"addq	$0x060,%%rsi	/* cc1 */				\n\t"\
		"\n\t"\
		"vmovaps	     (%%rcx),%%ymm4					\n\t		vmovaps	0x100(%%rcx),%%ymm12	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5					\n\t		vmovaps	0x120(%%rcx),%%ymm13	\n\t"\
		"vmovaps	     (%%rsi),%%ymm2	/* c32_1 */		\n\t		vmovaps	0x040(%%rsi),%%ymm10	/* c32_3 */\n\t"\
		"vmovaps	0x020(%%rsi),%%ymm3	/* s32_1 */		\n\t		vmovaps	0x060(%%rsi),%%ymm11	/* s32_3 */\n\t"\
		"vmovaps	%%ymm4,%%ymm6						\n\t		vmovaps	%%ymm12,%%ymm14		\n\t"\
		"vmovaps	%%ymm5,%%ymm7						\n\t		vmovaps	%%ymm13,%%ymm15		\n\t"\
		"\n\t"\
		"vmulpd	%%ymm10,%%ymm4,%%ymm4					\n\t		vmulpd	%%ymm3,%%ymm12,%%ymm12		\n\t"\
		"vmulpd	%%ymm10,%%ymm5,%%ymm5					\n\t		vmulpd	%%ymm3,%%ymm13,%%ymm13		\n\t"\
		"vmulpd	%%ymm11,%%ymm6,%%ymm6					\n\t		vmulpd	%%ymm2,%%ymm14,%%ymm14		\n\t"\
		"vmulpd	%%ymm11,%%ymm7,%%ymm7					\n\t		vmulpd	%%ymm2,%%ymm15,%%ymm15		\n\t"\
		"vmovaps	     (%%rdx),%%ymm0					\n\t		vmovaps	0x100(%%rdx),%%ymm8 	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm1					\n\t		vmovaps	0x120(%%rdx),%%ymm9 	\n\t"\
		"vsubpd	%%ymm6,%%ymm5,%%ymm5					\n\t		vsubpd	%%ymm14,%%ymm13,%%ymm13		\n\t"\
		"vmovaps	%%ymm0,%%ymm6						\n\t		vmovaps	%%ymm8 ,%%ymm14		\n\t"\
		"vaddpd	%%ymm7,%%ymm4,%%ymm4					\n\t		vaddpd	%%ymm15,%%ymm12,%%ymm12		\n\t"\
		"vmovaps	%%ymm1,%%ymm7						\n\t		vmovaps	%%ymm9 ,%%ymm15		\n\t"\
		"\n\t"\
		"vmulpd	%%ymm3,%%ymm0,%%ymm0					\n\t		vmulpd	%%ymm11,%%ymm8 ,%%ymm8 		\n\t"\
		"vmulpd	%%ymm3,%%ymm1,%%ymm1					\n\t		vmulpd	%%ymm11,%%ymm9 ,%%ymm9 		\n\t"\
		"vmulpd	%%ymm2,%%ymm6,%%ymm6					\n\t		vmulpd	%%ymm10,%%ymm14,%%ymm14		\n\t"\
		"vmulpd	%%ymm2,%%ymm7,%%ymm7					\n\t		vmulpd	%%ymm10,%%ymm15,%%ymm15		\n\t"\
		"vaddpd	%%ymm6,%%ymm1,%%ymm1					\n\t		vsubpd	%%ymm14,%%ymm9 ,%%ymm9 		\n\t"\
		"vsubpd	%%ymm7,%%ymm0,%%ymm0					\n\t		vaddpd	%%ymm15,%%ymm8 ,%%ymm8 		\n\t"\
		"\n\t"\
		"vmovaps	%%ymm5,%%ymm7						\n\t		vmovaps	%%ymm13,%%ymm15		\n\t"\
		"vmovaps	%%ymm4,%%ymm6						\n\t		vmovaps	%%ymm12,%%ymm14		\n\t"\
		"\n\t"\
		"vaddpd	%%ymm0,%%ymm6,%%ymm6					\n\t		vsubpd	%%ymm8 ,%%ymm12,%%ymm12		\n\t"\
		"vaddpd	%%ymm1,%%ymm7,%%ymm7					\n\t		vsubpd	%%ymm9 ,%%ymm13,%%ymm13		\n\t"\
		"vsubpd	%%ymm0,%%ymm4,%%ymm4					\n\t		vaddpd	%%ymm8 ,%%ymm14,%%ymm14		\n\t"\
		"vsubpd	%%ymm1,%%ymm5,%%ymm5					\n\t		vaddpd	%%ymm9 ,%%ymm15,%%ymm15		\n\t"\
		"\n\t"\
		"subq	$0x400,%%rcx			\n\t"\
		"subq	$0x400,%%rdx			\n\t"\
		"subq	$0x040,%%rsi	/* cc0 */	\n\t"\
		"vmovaps	     (%%rdx),%%ymm2					\n\t		vmovaps	0x100(%%rdx),%%ymm11	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm0					\n\t		vmovaps	0x120(%%rdx),%%ymm9 	\n\t"\
		"vmovaps	0x020(%%rsi),%%ymm3					\n\t		vmovaps	0x020(%%rsi),%%ymm10	\n\t"\
		"vmovaps	%%ymm2,%%ymm1						\n\t		vmovaps	%%ymm11,%%ymm8 		\n\t"\
		"vmulpd	%%ymm3 ,%%ymm2,%%ymm2					\n\t		vmulpd	%%ymm10,%%ymm11,%%ymm11		\n\t"\
		"vmulpd	%%ymm0 ,%%ymm3,%%ymm3					\n\t		vmulpd	%%ymm9 ,%%ymm10,%%ymm10		\n\t"\
		"vmulpd	(%%rsi),%%ymm1,%%ymm1					\n\t		vmulpd	(%%rsi),%%ymm8 ,%%ymm8 		\n\t"\
		"vmulpd	(%%rsi),%%ymm0,%%ymm0					\n\t		vmulpd	(%%rsi),%%ymm9 ,%%ymm9 		\n\t"\
		"vaddpd	%%ymm0 ,%%ymm2,%%ymm2					\n\t		vsubpd	%%ymm10,%%ymm8 ,%%ymm8 		\n\t"\
		"vsubpd	%%ymm1 ,%%ymm3,%%ymm3					\n\t		vaddpd	%%ymm11,%%ymm9 ,%%ymm9 		\n\t"\
		"\n\t"\
		"vmovaps	     (%%rcx),%%ymm0					\n\t		vmovaps	0x100(%%rcx),%%ymm10	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm1					\n\t		vmovaps	0x120(%%rcx),%%ymm11	\n\t"\
		"vsubpd	%%ymm2,%%ymm0,%%ymm0					\n\t		vsubpd	%%ymm8 ,%%ymm10,%%ymm10		\n\t"\
		"vsubpd	%%ymm3,%%ymm1,%%ymm1					\n\t		vsubpd	%%ymm9 ,%%ymm11,%%ymm11		\n\t"\
		"vmulpd	(%%r8),%%ymm2,%%ymm2					\n\t		vmulpd	(%%r8) ,%%ymm8 ,%%ymm8 		\n\t"\
		"vmulpd	(%%r8),%%ymm3,%%ymm3					\n\t		vmulpd	(%%r8) ,%%ymm9 ,%%ymm9 		\n\t"\
		"vaddpd	%%ymm0,%%ymm2,%%ymm2					\n\t		vaddpd	%%ymm10,%%ymm8 ,%%ymm8 		\n\t"\
		"vaddpd	%%ymm1,%%ymm3,%%ymm3					\n\t		vaddpd	%%ymm11,%%ymm9 ,%%ymm9 		\n\t"\
		"\n\t"\
		"/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c03): *//* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c07): */\n\t"\
		"addq	$0x6c0,%%rsi	/* c03 */				/* c07 = c03 + 0x100: */\n\t"\
		"vsubpd	%%ymm4,%%ymm2,%%ymm2					\n\t		vsubpd	%%ymm12,%%ymm10,%%ymm10		\n\t"\
		"vsubpd	%%ymm7,%%ymm0,%%ymm0					\n\t		vsubpd	%%ymm15,%%ymm8 ,%%ymm8 		\n\t"\
		"vsubpd	%%ymm5,%%ymm3,%%ymm3					\n\t		vsubpd	%%ymm13,%%ymm11,%%ymm11		\n\t"\
		"vsubpd	%%ymm6,%%ymm1,%%ymm1					\n\t		vsubpd	%%ymm14,%%ymm9 ,%%ymm9 		\n\t"\
		"vmulpd	(%%r8),%%ymm4,%%ymm4					\n\t		vmulpd	(%%r8) ,%%ymm12,%%ymm12		\n\t"\
		"vmulpd	(%%r8),%%ymm7,%%ymm7					\n\t		vmulpd	(%%r8) ,%%ymm15,%%ymm15		\n\t"\
		"vmulpd	(%%r8),%%ymm5,%%ymm5					\n\t		vmulpd	(%%r8) ,%%ymm13,%%ymm13		\n\t"\
		"vmulpd	(%%r8),%%ymm6,%%ymm6					\n\t		vmulpd	(%%r8) ,%%ymm14,%%ymm14		\n\t"\
		"vaddpd	%%ymm2,%%ymm4,%%ymm4					\n\t		vaddpd	%%ymm10,%%ymm12,%%ymm12		\n\t"\
		"vaddpd	%%ymm0,%%ymm7,%%ymm7					\n\t		vaddpd	%%ymm8 ,%%ymm15,%%ymm15		\n\t"\
		"vaddpd	%%ymm3,%%ymm5,%%ymm5					\n\t		vaddpd	%%ymm11,%%ymm13,%%ymm13		\n\t"\
		"vaddpd	%%ymm1,%%ymm6,%%ymm6					\n\t		vaddpd	%%ymm9 ,%%ymm14,%%ymm14		\n\t"\
		"vmovaps	%%ymm2,     (%%rcx)					\n\t		vmovaps	%%ymm10,0x100(%%rcx)	\n\t"\
		"vmovaps	%%ymm0,0x020(%%rdx)					\n\t		vmovaps	%%ymm8 ,0x120(%%rdx)	\n\t"\
		"vmovaps	%%ymm3,0x020(%%rcx)					\n\t		vmovaps	%%ymm11,0x120(%%rcx)	\n\t"\
		"vmovaps	%%ymm6,     (%%rdx)					\n\t		vmovaps	%%ymm14,0x100(%%rdx)	\n\t"\
		"vmovaps	%%ymm4,%%ymm2						\n\t		vmovaps	%%ymm12,%%ymm10		\n\t"\
		"vmovaps	%%ymm7,%%ymm0						\n\t		vmovaps	%%ymm15,%%ymm8 		\n\t"\
		"vmovaps	%%ymm5,%%ymm3						\n\t		vmovaps	%%ymm13,%%ymm11		\n\t"\
		"vmovaps	%%ymm1,%%ymm6						\n\t		vmovaps	%%ymm9 ,%%ymm14		\n\t"\
		"vmulpd	     (%%rsi),%%ymm4,%%ymm4				\n\t		vmulpd	0x100(%%rsi),%%ymm12,%%ymm12	\n\t"\
		"vmulpd	0x040(%%rsi),%%ymm7,%%ymm7				\n\t		vmulpd	0x140(%%rsi),%%ymm15,%%ymm15	\n\t"\
		"vmulpd	     (%%rsi),%%ymm5,%%ymm5				\n\t		vmulpd	0x100(%%rsi),%%ymm13,%%ymm13	\n\t"\
		"vmulpd	0x040(%%rsi),%%ymm1,%%ymm1				\n\t		vmulpd	0x140(%%rsi),%%ymm9 ,%%ymm9 	\n\t"\
		"vmulpd	0x020(%%rsi),%%ymm2,%%ymm2				\n\t		vmulpd	0x120(%%rsi),%%ymm10,%%ymm10	\n\t"\
		"vmulpd	0x060(%%rsi),%%ymm0,%%ymm0				\n\t		vmulpd	0x160(%%rsi),%%ymm8 ,%%ymm8 	\n\t"\
		"vmulpd	0x020(%%rsi),%%ymm3,%%ymm3				\n\t		vmulpd	0x120(%%rsi),%%ymm11,%%ymm11	\n\t"\
		"vmulpd	0x060(%%rsi),%%ymm6,%%ymm6				\n\t		vmulpd	0x160(%%rsi),%%ymm14,%%ymm14	\n\t"\
		"vsubpd	      %%ymm2,%%ymm5,%%ymm5				\n\t		vsubpd	%%ymm10,%%ymm13,%%ymm13		\n\t"\
		"vsubpd	      %%ymm0,%%ymm1,%%ymm1				\n\t		vsubpd	%%ymm8 ,%%ymm9 ,%%ymm9 		\n\t"\
		"vaddpd	      %%ymm3,%%ymm4,%%ymm4				\n\t		vaddpd	%%ymm11,%%ymm12,%%ymm12		\n\t"\
		"vaddpd	      %%ymm6,%%ymm7,%%ymm7				\n\t		vaddpd	%%ymm14,%%ymm15,%%ymm15		\n\t"\
		"vmovaps	%%ymm5,0x020(%%rax)					\n\t		vmovaps	%%ymm13,0x020(%%r10)	\n\t"\
		"vmovaps	%%ymm1,0x020(%%rbx)					\n\t		vmovaps	%%ymm9 ,0x020(%%r11)	\n\t"\
		"vmovaps	%%ymm4,     (%%rax)					\n\t		vmovaps	%%ymm12,     (%%r10)	\n\t"\
		"vmovaps	%%ymm7,     (%%rbx)					\n\t		vmovaps	%%ymm15,     (%%r11)	\n\t"\
		"addq	%%rdi,%%rax								\n\t		addq	%%rdi,%%r10			\n\t"\
		"addq	%%rdi,%%rbx								\n\t		addq	%%rdi,%%r11			\n\t"\
		"addq	$0x080,%%rsi	/* c13 */				\n\t		/* c17 = c13 + 0x100: */\n\t"\
		"vmovaps	     (%%rcx),%%ymm4					\n\t		vmovaps	0x100(%%rcx),%%ymm12	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm0					\n\t		vmovaps	0x120(%%rdx),%%ymm8 	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5					\n\t		vmovaps	0x120(%%rcx),%%ymm13	\n\t"\
		"vmovaps	     (%%rdx),%%ymm6					\n\t		vmovaps	0x100(%%rdx),%%ymm14	\n\t"\
		"vmovaps	%%ymm4,%%ymm2						\n\t		vmovaps	%%ymm12,%%ymm10		\n\t"\
		"vmovaps	%%ymm0,%%ymm1						\n\t		vmovaps	%%ymm8,%%ymm9 		\n\t"\
		"vmovaps	%%ymm5,%%ymm3						\n\t		vmovaps	%%ymm13,%%ymm11		\n\t"\
		"vmovaps	%%ymm6,%%ymm7						\n\t		vmovaps	%%ymm14,%%ymm15		\n\t"\
		"vmulpd	     (%%rsi),%%ymm4,%%ymm4				\n\t		vmulpd	0x100(%%rsi),%%ymm12,%%ymm12	\n\t"\
		"vmulpd	0x040(%%rsi),%%ymm0,%%ymm0				\n\t		vmulpd	0x140(%%rsi),%%ymm8 ,%%ymm8 	\n\t"\
		"vmulpd	     (%%rsi),%%ymm5,%%ymm5				\n\t		vmulpd	0x100(%%rsi),%%ymm13,%%ymm13	\n\t"\
		"vmulpd	0x040(%%rsi),%%ymm6,%%ymm6				\n\t		vmulpd	0x140(%%rsi),%%ymm14,%%ymm14	\n\t"\
		"vmulpd	0x020(%%rsi),%%ymm2,%%ymm2				\n\t		vmulpd	0x120(%%rsi),%%ymm10,%%ymm10	\n\t"\
		"vmulpd	0x060(%%rsi),%%ymm1,%%ymm1				\n\t		vmulpd	0x160(%%rsi),%%ymm9 ,%%ymm9 	\n\t"\
		"vmulpd	0x020(%%rsi),%%ymm3,%%ymm3				\n\t		vmulpd	0x120(%%rsi),%%ymm11,%%ymm11	\n\t"\
		"vmulpd	0x060(%%rsi),%%ymm7,%%ymm7				\n\t		vmulpd	0x160(%%rsi),%%ymm15,%%ymm15	\n\t"\
		"vsubpd	      %%ymm2,%%ymm5,%%ymm5				\n\t		vsubpd	%%ymm10,%%ymm13,%%ymm13		\n\t"\
		"vsubpd	      %%ymm1,%%ymm6,%%ymm6				\n\t		vsubpd	%%ymm9 ,%%ymm14,%%ymm14		\n\t"\
		"vaddpd	      %%ymm3,%%ymm4,%%ymm4				\n\t		vaddpd	%%ymm11,%%ymm12,%%ymm12		\n\t"\
		"vaddpd	      %%ymm7,%%ymm0,%%ymm0				\n\t		vaddpd	%%ymm15,%%ymm8 ,%%ymm8 		\n\t"\
		"vmovaps	%%ymm5,0x020(%%rax)					\n\t		vmovaps	%%ymm13,0x020(%%r10)	\n\t"\
		"vmovaps	%%ymm6,0x020(%%rbx)					\n\t		vmovaps	%%ymm14,0x020(%%r11)	\n\t"\
		"vmovaps	%%ymm4,     (%%rax)					\n\t		vmovaps	%%ymm12,     (%%r10)	\n\t"\
		"vmovaps	%%ymm0,     (%%rbx)					\n\t		vmovaps	%%ymm8 ,     (%%r11)	\n\t"\
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
		 ,[__r10] "m" (Xr10)\
		 ,[__r20] "m" (Xr20)\
		 ,[__r30] "m" (Xr30)\
		 ,[__isrt2] "m" (Xisrt2)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r8","r9","r10","r11","r12","r13","r14","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
	);\
	}

  #endif	// USE_64BIT_ASM_STYLE ?

#elif defined(USE_SSE2)

  #define USE_64BIT_ASM_STYLE	1	// My x86 timings indicate fancier versions below using xmm0-15 for the radix-32 DFT is faster.

  #if !USE_64BIT_ASM_STYLE	// False: Use simple 64-bit-ified version of the analogous 32-bit ASM macros, i.e. using just xmm0-7.

	/*
	For GCC-macro version of this, use that isrt2 + 0x010,0x030,0x050 = cc0,cc1,cc3,
	and isrt2 + 0x070,0x0f0,0x170,0x1f0,0x270,0x2f0,0x370,0x3f0 = c00,04,02,06,01,05,03,07
	in order to reduce number of args to <= the GCC-allowed maximum of 30:
	*/
	#define SSE2_RADIX32_DIF_TWIDDLE(Xadd0,Xp01,Xp02,Xp03,Xp04,Xp08,Xp0C,Xp10,Xp18,Xr00)\
	{\
	__asm__ volatile (\
		"/*...Block 1: */	\n\t"\
		"movq	%[__add0],%%rax	\n\t"\
		"movslq	%[__p08],%%rbx	\n\t"\
		"movslq	%[__p10],%%rcx	\n\t"\
		"movslq	%[__p18],%%rdx	\n\t"\
		"shlq	$3,%%rbx	\n\t"\
		"shlq	$3,%%rcx	\n\t"\
		"shlq	$3,%%rdx	\n\t"\
		"addq	%%rax,%%rbx	\n\t"\
		"addq	%%rax,%%rcx	\n\t"\
		"addq	%%rax,%%rdx	\n\t"\
		"/* SSE2_RADIX4_DIF_4TWIDDLE_A(r00,c10) */	\n\t"\
		"movq	%[__r00],%%rsi	\n\t"\
		"movaps	    (%%rax),%%xmm0	\n\t"\
		"movaps	    (%%rcx),%%xmm4	\n\t"\
		"movaps	0x10(%%rax),%%xmm1	\n\t"\
		"movaps	0x10(%%rcx),%%xmm5	\n\t"\
		"movaps	0x490(%%rsi),%%xmm2	/* c10 */	\n\t"\
		"movaps	0x4a0(%%rsi),%%xmm3	\n\t"\
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
		"movaps	    (%%rdx),%%xmm4	\n\t"\
		"movaps	0x10(%%rdx),%%xmm5	\n\t"\
		"movaps	    (%%rdx),%%xmm6	\n\t"\
		"movaps	0x10(%%rdx),%%xmm7	\n\t"\
		"mulpd	0x4d0(%%rsi),%%xmm4	/* c18 */	\n\t"\
		"mulpd	0x4d0(%%rsi),%%xmm5	\n\t"\
		"mulpd	0x4e0(%%rsi),%%xmm6	\n\t"\
		"mulpd	0x4e0(%%rsi),%%xmm7	\n\t"\
		"addpd	%%xmm6,%%xmm5	\n\t"\
		"subpd	%%xmm7,%%xmm4	\n\t"\
		"movaps	%%xmm5,0x010(%%rsi)	\n\t"\
		"movaps	%%xmm4,     (%%rsi)	/* r00 */	\n\t"\
		"movaps	    (%%rbx),%%xmm4	\n\t"\
		"movaps	0x10(%%rbx),%%xmm5	\n\t"\
		"movaps	    (%%rbx),%%xmm6	\n\t"\
		"movaps	0x10(%%rbx),%%xmm7	\n\t"\
		"mulpd	0x4b0(%%rsi),%%xmm4	/* c08 */	\n\t"\
		"mulpd	0x4b0(%%rsi),%%xmm5	\n\t"\
		"mulpd	0x4c0(%%rsi),%%xmm6	\n\t"\
		"mulpd	0x4c0(%%rsi),%%xmm7	\n\t"\
		"addpd	%%xmm6,%%xmm5	\n\t"\
		"subpd	%%xmm7,%%xmm4	\n\t"\
		"movaps	%%xmm5,%%xmm7	\n\t"\
		"movaps	%%xmm4,%%xmm6	\n\t"\
		"subpd	    (%%rsi),%%xmm4	/* r00 */	\n\t"\
		"subpd	0x10(%%rsi),%%xmm5	\n\t"\
		"addpd	    (%%rsi),%%xmm6	\n\t"\
		"addpd	0x10(%%rsi),%%xmm7	\n\t"\
		"subpd	%%xmm6,%%xmm0	\n\t"\
		"subpd	%%xmm5,%%xmm2	\n\t"\
		"subpd	%%xmm7,%%xmm1	\n\t"\
		"subpd	%%xmm4,%%xmm3	\n\t"\
		"movaps	%%xmm0,0x040(%%rsi)	\n\t"\
		"movaps	%%xmm2,0x020(%%rsi)	\n\t"\
		"movaps	%%xmm1,0x050(%%rsi)	\n\t"\
		"movaps	%%xmm3,0x070(%%rsi)	\n\t"\
		"addpd	%%xmm6,%%xmm6	\n\t"\
		"addpd	%%xmm5,%%xmm5	\n\t"\
		"addpd	%%xmm7,%%xmm7	\n\t"\
		"addpd	%%xmm4,%%xmm4	\n\t"\
		"addpd	%%xmm0,%%xmm6	\n\t"\
		"addpd	%%xmm2,%%xmm5	\n\t"\
		"addpd	%%xmm1,%%xmm7	\n\t"\
		"addpd	%%xmm3,%%xmm4	\n\t"\
		"movaps	%%xmm6,     (%%rsi)	\n\t"\
		"movaps	%%xmm5,0x060(%%rsi)	\n\t"\
		"movaps	%%xmm7,0x010(%%rsi)	\n\t"\
		"movaps	%%xmm4,0x030(%%rsi)	\n\t"\
		"/***************************************/	\n\t"\
		"movslq	%[__p04],%%rdi	\n\t"\
		"shlq	$3,%%rdi	\n\t"\
		"addq	%%rdi,%%rax	/* &a[j1+p4] */	\n\t"\
		"addq	%%rdi,%%rbx	\n\t"\
		"addq	%%rdi,%%rcx	\n\t"\
		"addq	%%rdi,%%rdx	\n\t"\
		"/* SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO_B(r08,c04) */	\n\t"\
		"movaps	    (%%rax),%%xmm0	\n\t"\
		"movaps	    (%%rcx),%%xmm4	\n\t"\
		"movaps	0x10(%%rax),%%xmm1	\n\t"\
		"movaps	0x10(%%rcx),%%xmm5	\n\t"\
		"movaps	0x4f0(%%rsi),%%xmm6	/* c04 */	\n\t"\
		"movaps	0x500(%%rsi),%%xmm7	\n\t"\
		"movaps	%%xmm0,%%xmm2	\n\t"\
		"movaps	%%xmm1,%%xmm3	\n\t"\
		"mulpd	%%xmm6,%%xmm0	\n\t"\
		"mulpd	%%xmm6,%%xmm1	\n\t"\
		"mulpd	%%xmm7,%%xmm2	\n\t"\
		"mulpd	%%xmm7,%%xmm3	\n\t"\
		"movaps	%%xmm4,%%xmm6	\n\t"\
		"addpd	%%xmm2,%%xmm1	\n\t"\
		"movaps	%%xmm5,%%xmm7	\n\t"\
		"mulpd	0x510(%%rsi),%%xmm4	/* c14 */	\n\t"\
		"subpd	%%xmm3,%%xmm0	\n\t"\
		"mulpd	0x510(%%rsi),%%xmm5	\n\t"\
		"mulpd	0x520(%%rsi),%%xmm6	\n\t"\
		"movaps	%%xmm0,%%xmm2	\n\t"\
		"mulpd	0x520(%%rsi),%%xmm7	\n\t"\
		"addpd	%%xmm6,%%xmm5	\n\t"\
		"movaps	%%xmm1,%%xmm3	\n\t"\
		"subpd	%%xmm7,%%xmm4	\n\t"\
		"addpd	%%xmm4,%%xmm0	\n\t"\
		"addpd	%%xmm5,%%xmm1	\n\t"\
		"subpd	%%xmm4,%%xmm2	\n\t"\
		"subpd	%%xmm5,%%xmm3	\n\t"\
		"movaps	    (%%rdx),%%xmm4	\n\t"\
		"movaps	0x10(%%rdx),%%xmm5	\n\t"\
		"movaps	    (%%rdx),%%xmm6	\n\t"\
		"movaps	0x10(%%rdx),%%xmm7	\n\t"\
		"mulpd	0x550(%%rsi),%%xmm4	/* c1C */	\n\t"\
		"mulpd	0x550(%%rsi),%%xmm5	\n\t"\
		"mulpd	0x560(%%rsi),%%xmm6	\n\t"\
		"mulpd	0x560(%%rsi),%%xmm7	\n\t"\
		"addq	$0x80,%%rsi	/* r08 */	\n\t"\
		"addpd	%%xmm6,%%xmm5	\n\t"\
		"subpd	%%xmm7,%%xmm4	\n\t"\
		"movaps	%%xmm5,0x010(%%rsi)	\n\t"\
		"movaps	%%xmm4,     (%%rsi)	\n\t"\
		"movaps	    (%%rbx),%%xmm4	\n\t"\
		"movaps	0x10(%%rbx),%%xmm5	\n\t"\
		"movaps	    (%%rbx),%%xmm6	\n\t"\
		"movaps	0x10(%%rbx),%%xmm7	\n\t"\
		"mulpd	0x4b0(%%rsi),%%xmm4	/* c0C */	\n\t"\
		"mulpd	0x4b0(%%rsi),%%xmm5	\n\t"\
		"mulpd	0x4c0(%%rsi),%%xmm6	\n\t"\
		"mulpd	0x4c0(%%rsi),%%xmm7	\n\t"\
		"addpd	%%xmm6,%%xmm5	\n\t"\
		"subpd	%%xmm7,%%xmm4	\n\t"\
		"movaps	%%xmm5,%%xmm7	\n\t"\
		"movaps	%%xmm4,%%xmm6	\n\t"\
		"subpd	     (%%rsi),%%xmm4	\n\t"\
		"subpd	0x010(%%rsi),%%xmm5	\n\t"\
		"addpd	     (%%rsi),%%xmm6	\n\t"\
		"addpd	0x010(%%rsi),%%xmm7	\n\t"\
		"subpd	%%xmm6,%%xmm0	\n\t"\
		"subpd	%%xmm7,%%xmm1	\n\t"\
		"movaps	%%xmm0,0x040(%%rsi)	\n\t"\
		"subpd	%%xmm5,%%xmm2	\n\t"\
		"movaps	%%xmm1,0x050(%%rsi)	\n\t"\
		"subpd	%%xmm4,%%xmm3	\n\t"\
		"addpd	%%xmm6,%%xmm6	\n\t"\
		"addpd	%%xmm5,%%xmm5	\n\t"\
		"addpd	%%xmm7,%%xmm7	\n\t"\
		"addpd	%%xmm4,%%xmm4	\n\t"\
		"addpd	%%xmm0,%%xmm6	\n\t"\
		"addpd	%%xmm2,%%xmm5	\n\t"\
		"addpd	%%xmm1,%%xmm7	\n\t"\
		"addpd	%%xmm3,%%xmm4	\n\t"\
		"movaps	%%xmm6,     (%%rsi)	\n\t"\
		"movaps	%%xmm7,0x010(%%rsi)	\n\t"\
		"movaps	0x380(%%rsi),%%xmm0	/* isrt2 */	\n\t"\
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
		"movaps	%%xmm2,0x020(%%rsi)	\n\t"\
		"movaps	%%xmm5,0x060(%%rsi)	\n\t"\
		"movaps	%%xmm6,0x030(%%rsi)	\n\t"\
		"movaps	%%xmm7,0x070(%%rsi)	\n\t"\
		"/***************************************/	\n\t"\
		"	\n\t"\
		"/* SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r00) */	\n\t"\
		"subq	$0x80,%%rsi	/* r00 */	\n\t"\
		"movaps	    (%%rsi),%%xmm0	\n\t"\
		"movaps	0x40(%%rsi),%%xmm4	\n\t"\
		"movaps	0x10(%%rsi),%%xmm1	\n\t"\
		"movaps	0x50(%%rsi),%%xmm5	\n\t"\
		"movaps	0x80(%%rsi),%%xmm2	\n\t"\
		"movaps	0xd0(%%rsi),%%xmm7	\n\t"\
		"movaps	0x90(%%rsi),%%xmm3	\n\t"\
		"movaps	0xc0(%%rsi),%%xmm6	\n\t"\
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
		"movaps	%%xmm0,0x80(%%rsi)	\n\t"\
		"movaps	%%xmm4,0x40(%%rsi)	\n\t"\
		"movaps	%%xmm1,0x90(%%rsi)	\n\t"\
		"movaps	%%xmm5,0xd0(%%rsi)	\n\t"\
		"movaps	%%xmm2,    (%%rsi)	\n\t"\
		"movaps	%%xmm7,0xc0(%%rsi)	\n\t"\
		"movaps	%%xmm3,0x10(%%rsi)	\n\t"\
		"movaps	%%xmm6,0x50(%%rsi)	\n\t"\
		"movaps	0x20(%%rsi),%%xmm0	\n\t"\
		"movaps	0x60(%%rsi),%%xmm4	\n\t"\
		"movaps	0x30(%%rsi),%%xmm1	\n\t"\
		"movaps	0x70(%%rsi),%%xmm5	\n\t"\
		"movaps	0xa0(%%rsi),%%xmm2	\n\t"\
		"movaps	0xf0(%%rsi),%%xmm7	\n\t"\
		"movaps	0xb0(%%rsi),%%xmm3	\n\t"\
		"movaps	0xe0(%%rsi),%%xmm6	\n\t"\
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
		"movaps	%%xmm0,0xa0(%%rsi)	\n\t"\
		"movaps	%%xmm4,0x60(%%rsi)	\n\t"\
		"movaps	%%xmm1,0xb0(%%rsi)	\n\t"\
		"movaps	%%xmm5,0xf0(%%rsi)	\n\t"\
		"movaps	%%xmm2,0x20(%%rsi)	\n\t"\
		"movaps	%%xmm7,0xe0(%%rsi)	\n\t"\
		"movaps	%%xmm3,0x30(%%rsi)	\n\t"\
		"movaps	%%xmm6,0x70(%%rsi)	\n\t"\
		"/***************************************/	\n\t"\
		"	\n\t"\
		"/*...Block 2: */	\n\t"\
		"subq	%%rdi,%%rax	/* &a[j1] */	\n\t"\
		"subq	%%rdi,%%rbx	\n\t"\
		"subq	%%rdi,%%rcx	\n\t"\
		"subq	%%rdi,%%rdx	\n\t"\
		"movslq	%[__p02],%%rdi	/* Do things this way [rather than repeatedly adding p1] since array-padding scheme means p2 == p+p1 not guaranteed. */	\n\t"\
		"shlq	$3,%%rdi	\n\t"\
		"addq	%%rdi,%%rax	/* &a[j1+p2] */	\n\t"\
		"addq	%%rdi,%%rbx	\n\t"\
		"addq	%%rdi,%%rcx	\n\t"\
		"addq	%%rdi,%%rdx	\n\t"\
		"/* SSE2_RADIX4_DIF_4TWIDDLE_B(r10,c02) */	\n\t"\
		"addq	$0x100,%%rsi	/* r10 */	\n\t"\
		"movaps	    (%%rax),%%xmm0	\n\t"\
		"movaps	    (%%rcx),%%xmm4	\n\t"\
		"movaps	0x10(%%rax),%%xmm1	\n\t"\
		"movaps	0x10(%%rcx),%%xmm5	\n\t"\
		"movaps	0x470(%%rsi),%%xmm6	/* c02 */	\n\t"\
		"movaps	0x480(%%rsi),%%xmm7	\n\t"\
		"movaps	%%xmm0,%%xmm2	\n\t"\
		"movaps	%%xmm1,%%xmm3	\n\t"\
		"mulpd   %%xmm6,%%xmm0	\n\t"\
		"mulpd   %%xmm6,%%xmm1	\n\t"\
		"mulpd   %%xmm7,%%xmm2	\n\t"\
		"mulpd   %%xmm7,%%xmm3	\n\t"\
		"movaps	%%xmm4,%%xmm6	\n\t"\
		"addpd   %%xmm2,%%xmm1	\n\t"\
		"movaps	%%xmm5,%%xmm7	\n\t"\
		"mulpd   0x490(%%rsi),%%xmm4	/* c12 */	\n\t"\
		"subpd   %%xmm3,%%xmm0	\n\t"\
		"mulpd   0x490(%%rsi),%%xmm5	\n\t"\
		"mulpd   0x4a0(%%rsi),%%xmm6	\n\t"\
		"movaps	%%xmm0,%%xmm2	\n\t"\
		"mulpd	0x4a0(%%rsi),%%xmm7	\n\t"\
		"addpd	%%xmm6,%%xmm5	\n\t"\
		"movaps	%%xmm1,%%xmm3	\n\t"\
		"subpd	%%xmm7,%%xmm4	\n\t"\
		"addpd	%%xmm4,%%xmm0	\n\t"\
		"addpd	%%xmm5,%%xmm1	\n\t"\
		"subpd	%%xmm4,%%xmm2	\n\t"\
		"subpd	%%xmm5,%%xmm3	\n\t"\
		"movaps	    (%%rdx),%%xmm4	\n\t"\
		"movaps	0x10(%%rdx),%%xmm5	\n\t"\
		"movaps	    (%%rdx),%%xmm6	\n\t"\
		"movaps	0x10(%%rdx),%%xmm7	\n\t"\
		"mulpd	0x4d0(%%rsi),%%xmm4	/* c1A */	\n\t"\
		"mulpd	0x4d0(%%rsi),%%xmm5	\n\t"\
		"mulpd	0x4e0(%%rsi),%%xmm6	\n\t"\
		"mulpd	0x4e0(%%rsi),%%xmm7	\n\t"\
		"addpd	%%xmm6,%%xmm5	\n\t"\
		"subpd	%%xmm7,%%xmm4	\n\t"\
		"movaps	%%xmm5,0x010(%%rsi)	\n\t"\
		"movaps	%%xmm4,     (%%rsi)	\n\t"\
		"movaps	    (%%rbx),%%xmm4	\n\t"\
		"movaps	0x10(%%rbx),%%xmm5	\n\t"\
		"movaps	    (%%rbx),%%xmm6	\n\t"\
		"movaps	0x10(%%rbx),%%xmm7	\n\t"\
		"mulpd	0x4b0(%%rsi),%%xmm4	/* c0A */	\n\t"\
		"mulpd	0x4b0(%%rsi),%%xmm5	\n\t"\
		"mulpd	0x4c0(%%rsi),%%xmm6	\n\t"\
		"mulpd	0x4c0(%%rsi),%%xmm7	\n\t"\
		"addpd	%%xmm6,%%xmm5	\n\t"\
		"subpd	%%xmm7,%%xmm4	\n\t"\
		"movaps	%%xmm5,%%xmm7	\n\t"\
		"movaps	%%xmm4,%%xmm6	\n\t"\
		"subpd	    (%%rsi),%%xmm4	\n\t"\
		"subpd	0x10(%%rsi),%%xmm5	\n\t"\
		"addpd	    (%%rsi),%%xmm6	\n\t"\
		"addpd	0x10(%%rsi),%%xmm7	\n\t"\
		"subpd	%%xmm6,%%xmm0	\n\t"\
		"subpd	%%xmm5,%%xmm2	\n\t"\
		"subpd	%%xmm7,%%xmm1	\n\t"\
		"subpd	%%xmm4,%%xmm3	\n\t"\
		"movaps	%%xmm0,0x040(%%rsi)	\n\t"\
		"movaps	%%xmm2,0x020(%%rsi)	\n\t"\
		"movaps	%%xmm1,0x050(%%rsi)	\n\t"\
		"movaps	%%xmm3,0x070(%%rsi)	\n\t"\
		"addpd	%%xmm6,%%xmm6	\n\t"\
		"addpd	%%xmm5,%%xmm5	\n\t"\
		"addpd	%%xmm7,%%xmm7	\n\t"\
		"addpd	%%xmm4,%%xmm4	\n\t"\
		"addpd	%%xmm0,%%xmm6	\n\t"\
		"addpd	%%xmm2,%%xmm5	\n\t"\
		"addpd	%%xmm1,%%xmm7	\n\t"\
		"addpd	%%xmm3,%%xmm4	\n\t"\
		"movaps	%%xmm6,     (%%rsi)	\n\t"\
		"movaps	%%xmm5,0x060(%%rsi)	\n\t"\
		"movaps	%%xmm7,0x010(%%rsi)	\n\t"\
		"movaps	%%xmm4,0x030(%%rsi)	\n\t"\
		"/***************************************/	\n\t"\
		"movslq	%[__p04],%%rdi	\n\t"\
		"shlq	$3,%%rdi	\n\t"\
		"addq	%%rdi,%%rax	/* &a[j1+p6] */	\n\t"\
		"addq	%%rdi,%%rbx	\n\t"\
		"addq	%%rdi,%%rcx	\n\t"\
		"addq	%%rdi,%%rdx	\n\t"\
		"/* SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO_B(r18,c06) */	\n\t"\
		"/* rsi contains r10 */	\n\t"\
		"movaps	    (%%rax),%%xmm0	\n\t"\
		"movaps	    (%%rcx),%%xmm4	\n\t"\
		"movaps	0x10(%%rax),%%xmm1	\n\t"\
		"movaps	0x10(%%rcx),%%xmm5	\n\t"\
		"movaps	0x4f0(%%rsi),%%xmm6	/* c06 */	\n\t"\
		"movaps	0x500(%%rsi),%%xmm7	\n\t"\
		"movaps	%%xmm0,%%xmm2	\n\t"\
		"movaps	%%xmm1,%%xmm3	\n\t"\
		"mulpd	%%xmm6,%%xmm0	\n\t"\
		"mulpd	%%xmm6,%%xmm1	\n\t"\
		"mulpd	%%xmm7,%%xmm2	\n\t"\
		"mulpd	%%xmm7,%%xmm3	\n\t"\
		"movaps	%%xmm4,%%xmm6	\n\t"\
		"addpd	%%xmm2,%%xmm1	\n\t"\
		"movaps	%%xmm5,%%xmm7	\n\t"\
		"mulpd	0x510(%%rsi),%%xmm4	/* c16 */	\n\t"\
		"subpd	%%xmm3,%%xmm0	\n\t"\
		"mulpd	0x510(%%rsi),%%xmm5	\n\t"\
		"mulpd	0x520(%%rsi),%%xmm6	\n\t"\
		"movaps	%%xmm0,%%xmm2	\n\t"\
		"mulpd	0x520(%%rsi),%%xmm7	\n\t"\
		"addpd	%%xmm6,%%xmm5	\n\t"\
		"movaps	%%xmm1,%%xmm3	\n\t"\
		"subpd	%%xmm7,%%xmm4	\n\t"\
		"addpd	%%xmm4,%%xmm0	\n\t"\
		"addpd	%%xmm5,%%xmm1	\n\t"\
		"subpd	%%xmm4,%%xmm2	\n\t"\
		"subpd	%%xmm5,%%xmm3	\n\t"\
		"movaps	    (%%rdx),%%xmm4	\n\t"\
		"movaps	0x10(%%rdx),%%xmm5	\n\t"\
		"movaps	    (%%rdx),%%xmm6	\n\t"\
		"movaps	0x10(%%rdx),%%xmm7	\n\t"\
		"mulpd	0x550(%%rsi),%%xmm4	/* c1E */	\n\t"\
		"mulpd	0x550(%%rsi),%%xmm5	\n\t"\
		"mulpd	0x560(%%rsi),%%xmm6	\n\t"\
		"mulpd	0x560(%%rsi),%%xmm7	\n\t"\
		"addq	$0x80,%%rsi	/* r18 */	\n\t"\
		"addpd	%%xmm6,%%xmm5	\n\t"\
		"subpd	%%xmm7,%%xmm4	\n\t"\
		"movaps	%%xmm5,0x010(%%rsi)	\n\t"\
		"movaps	%%xmm4,     (%%rsi)	\n\t"\
		"movaps	    (%%rbx),%%xmm4	\n\t"\
		"movaps	0x10(%%rbx),%%xmm5	\n\t"\
		"movaps	    (%%rbx),%%xmm6	\n\t"\
		"movaps	0x10(%%rbx),%%xmm7	\n\t"\
		"mulpd	0x4b0(%%rsi),%%xmm4	/* c0E */	\n\t"\
		"mulpd	0x4b0(%%rsi),%%xmm5	\n\t"\
		"mulpd	0x4c0(%%rsi),%%xmm6	\n\t"\
		"mulpd	0x4c0(%%rsi),%%xmm7	\n\t"\
		"addpd	%%xmm6,%%xmm5	\n\t"\
		"subpd	%%xmm7,%%xmm4	\n\t"\
		"movaps	%%xmm5,%%xmm7	\n\t"\
		"movaps	%%xmm4,%%xmm6	\n\t"\
		"subpd	     (%%rsi),%%xmm4	\n\t"\
		"subpd	0x010(%%rsi),%%xmm5	\n\t"\
		"addpd	     (%%rsi),%%xmm6	\n\t"\
		"addpd	0x010(%%rsi),%%xmm7	\n\t"\
		"subpd	%%xmm6,%%xmm0	\n\t"\
		"subpd	%%xmm7,%%xmm1	\n\t"\
		"movaps	%%xmm0,0x040(%%rsi)	\n\t"\
		"subpd	%%xmm5,%%xmm2	\n\t"\
		"movaps	%%xmm1,0x050(%%rsi)	\n\t"\
		"subpd	%%xmm4,%%xmm3	\n\t"\
		"addpd	%%xmm6,%%xmm6	\n\t"\
		"addpd	%%xmm5,%%xmm5	\n\t"\
		"addpd	%%xmm7,%%xmm7	\n\t"\
		"addpd	%%xmm4,%%xmm4	\n\t"\
		"addpd	%%xmm0,%%xmm6	\n\t"\
		"addpd	%%xmm2,%%xmm5	\n\t"\
		"addpd	%%xmm1,%%xmm7	\n\t"\
		"addpd	%%xmm3,%%xmm4	\n\t"\
		"movaps	%%xmm6,     (%%rsi)	\n\t"\
		"movaps	%%xmm7,0x010(%%rsi)	\n\t"\
		"movaps	0x280(%%rsi),%%xmm0	/* isrt2 */	\n\t"\
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
		"movaps	%%xmm2,0x020(%%rsi)	\n\t"\
		"movaps	%%xmm5,0x060(%%rsi)	\n\t"\
		"movaps	%%xmm6,0x030(%%rsi)	\n\t"\
		"movaps	%%xmm7,0x070(%%rsi)	\n\t"\
		"/***************************************/	\n\t"\
		"	\n\t"\
		"/* SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r10) */	\n\t"\
		"subq		$0x80,%%rsi	/* r10 */	\n\t"\
		"movaps		    (%%rsi),%%xmm0	\n\t"\
		"movaps		0x40(%%rsi),%%xmm4	\n\t"\
		"movaps		0x10(%%rsi),%%xmm1	\n\t"\
		"movaps		0x50(%%rsi),%%xmm5	\n\t"\
		"movaps		0x80(%%rsi),%%xmm2	\n\t"\
		"movaps		0xd0(%%rsi),%%xmm7	\n\t"\
		"movaps		0x90(%%rsi),%%xmm3	\n\t"\
		"movaps		0xc0(%%rsi),%%xmm6	\n\t"\
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
		"movaps		%%xmm0,0x80(%%rsi)	\n\t"\
		"movaps		%%xmm4,0x40(%%rsi)	\n\t"\
		"movaps		%%xmm1,0x90(%%rsi)	\n\t"\
		"movaps		%%xmm5,0xd0(%%rsi)	\n\t"\
		"movaps		%%xmm2,    (%%rsi)	\n\t"\
		"movaps		%%xmm7,0xc0(%%rsi)	\n\t"\
		"movaps		%%xmm3,0x10(%%rsi)	\n\t"\
		"movaps		%%xmm6,0x50(%%rsi)	\n\t"\
		"movaps		0x20(%%rsi),%%xmm0	\n\t"\
		"movaps		0x60(%%rsi),%%xmm4	\n\t"\
		"movaps		0x30(%%rsi),%%xmm1	\n\t"\
		"movaps		0x70(%%rsi),%%xmm5	\n\t"\
		"movaps		0xa0(%%rsi),%%xmm2	\n\t"\
		"movaps		0xf0(%%rsi),%%xmm7	\n\t"\
		"movaps		0xb0(%%rsi),%%xmm3	\n\t"\
		"movaps		0xe0(%%rsi),%%xmm6	\n\t"\
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
		"movaps		%%xmm0,0xa0(%%rsi)	\n\t"\
		"movaps		%%xmm4,0x60(%%rsi)	\n\t"\
		"movaps		%%xmm1,0xb0(%%rsi)	\n\t"\
		"movaps		%%xmm5,0xf0(%%rsi)	\n\t"\
		"movaps		%%xmm2,0x20(%%rsi)	\n\t"\
		"movaps		%%xmm7,0xe0(%%rsi)	\n\t"\
		"movaps		%%xmm3,0x30(%%rsi)	\n\t"\
		"movaps		%%xmm6,0x70(%%rsi)	\n\t"\
		"/***************************************/	\n\t"\
		"	\n\t"\
		"/*...Block 3: */	\n\t"\
		"movq	%[__add0],%%rax	\n\t"\
		"movslq	%[__p01],%%rdi	/* Do things this way [rather than repeatedly adding p1] since array-padding scheme means p2 == p+p1 not guaranteed. */	\n\t"\
		"movslq	%[__p08],%%rbx	/* Can't get these via simple load-one-and-shift-as-needed due to array padding scheme */	\n\t"\
		"movslq	%[__p10],%%rcx	\n\t"\
		"movslq	%[__p18],%%rdx	\n\t"\
		"shlq	$3,%%rdi	\n\t"\
		"shlq	$3,%%rbx	\n\t"\
		"shlq	$3,%%rcx	\n\t"\
		"shlq	$3,%%rdx	\n\t"\
		"addq	%%rdi,%%rax	/* &a[j1+p1] */	\n\t"\
		"addq	%%rax,%%rbx	\n\t"\
		"addq	%%rax,%%rcx	\n\t"\
		"addq	%%rax,%%rdx	\n\t"\
		"/* SSE2_RADIX4_DIF_4TWIDDLE_B(r20,c01) */	\n\t"\
		"addq	$0x100,%%rsi	/* r20 */	\n\t"\
		"movaps	    (%%rax),%%xmm0	\n\t"\
		"movaps	    (%%rcx),%%xmm4	\n\t"\
		"movaps	0x10(%%rax),%%xmm1	\n\t"\
		"movaps	0x10(%%rcx),%%xmm5	\n\t"\
		"movaps	0x470(%%rsi),%%xmm6	/* c01 */	\n\t"\
		"movaps	0x480(%%rsi),%%xmm7	\n\t"\
		"movaps	%%xmm0,%%xmm2	\n\t"\
		"movaps	%%xmm1,%%xmm3	\n\t"\
		"mulpd   %%xmm6,%%xmm0	\n\t"\
		"mulpd   %%xmm6,%%xmm1	\n\t"\
		"mulpd   %%xmm7,%%xmm2	\n\t"\
		"mulpd   %%xmm7,%%xmm3	\n\t"\
		"movaps	%%xmm4,%%xmm6	\n\t"\
		"addpd   %%xmm2,%%xmm1	\n\t"\
		"movaps	%%xmm5,%%xmm7	\n\t"\
		"mulpd   0x490(%%rsi),%%xmm4	/* c11 */	\n\t"\
		"subpd   %%xmm3,%%xmm0	\n\t"\
		"mulpd   0x490(%%rsi),%%xmm5	\n\t"\
		"mulpd   0x4a0(%%rsi),%%xmm6	\n\t"\
		"movaps	%%xmm0,%%xmm2	\n\t"\
		"mulpd	0x4a0(%%rsi),%%xmm7	\n\t"\
		"addpd	%%xmm6,%%xmm5	\n\t"\
		"movaps	%%xmm1,%%xmm3	\n\t"\
		"subpd	%%xmm7,%%xmm4	\n\t"\
		"addpd	%%xmm4,%%xmm0	\n\t"\
		"addpd	%%xmm5,%%xmm1	\n\t"\
		"subpd	%%xmm4,%%xmm2	\n\t"\
		"subpd	%%xmm5,%%xmm3	\n\t"\
		"movaps	    (%%rdx),%%xmm4	\n\t"\
		"movaps	0x10(%%rdx),%%xmm5	\n\t"\
		"movaps	    (%%rdx),%%xmm6	\n\t"\
		"movaps	0x10(%%rdx),%%xmm7	\n\t"\
		"mulpd	0x4d0(%%rsi),%%xmm4	/* c19 */	\n\t"\
		"mulpd	0x4d0(%%rsi),%%xmm5	\n\t"\
		"mulpd	0x4e0(%%rsi),%%xmm6	\n\t"\
		"mulpd	0x4e0(%%rsi),%%xmm7	\n\t"\
		"addpd	%%xmm6,%%xmm5	\n\t"\
		"subpd	%%xmm7,%%xmm4	\n\t"\
		"movaps	%%xmm5,0x010(%%rsi)	\n\t"\
		"movaps	%%xmm4,     (%%rsi)	\n\t"\
		"movaps	    (%%rbx),%%xmm4	\n\t"\
		"movaps	0x10(%%rbx),%%xmm5	\n\t"\
		"movaps	    (%%rbx),%%xmm6	\n\t"\
		"movaps	0x10(%%rbx),%%xmm7	\n\t"\
		"mulpd	0x4b0(%%rsi),%%xmm4	/* c01 */	\n\t"\
		"mulpd	0x4b0(%%rsi),%%xmm5	\n\t"\
		"mulpd	0x4c0(%%rsi),%%xmm6	\n\t"\
		"mulpd	0x4c0(%%rsi),%%xmm7	\n\t"\
		"addpd	%%xmm6,%%xmm5	\n\t"\
		"subpd	%%xmm7,%%xmm4	\n\t"\
		"movaps	%%xmm5,%%xmm7	\n\t"\
		"movaps	%%xmm4,%%xmm6	\n\t"\
		"subpd	    (%%rsi),%%xmm4	\n\t"\
		"subpd	0x10(%%rsi),%%xmm5	\n\t"\
		"addpd	    (%%rsi),%%xmm6	\n\t"\
		"addpd	0x10(%%rsi),%%xmm7	\n\t"\
		"subpd	%%xmm6,%%xmm0	\n\t"\
		"subpd	%%xmm5,%%xmm2	\n\t"\
		"subpd	%%xmm7,%%xmm1	\n\t"\
		"subpd	%%xmm4,%%xmm3	\n\t"\
		"movaps	%%xmm0,0x040(%%rsi)	\n\t"\
		"movaps	%%xmm2,0x020(%%rsi)	\n\t"\
		"movaps	%%xmm1,0x050(%%rsi)	\n\t"\
		"movaps	%%xmm3,0x070(%%rsi)	\n\t"\
		"addpd	%%xmm6,%%xmm6	\n\t"\
		"addpd	%%xmm5,%%xmm5	\n\t"\
		"addpd	%%xmm7,%%xmm7	\n\t"\
		"addpd	%%xmm4,%%xmm4	\n\t"\
		"addpd	%%xmm0,%%xmm6	\n\t"\
		"addpd	%%xmm2,%%xmm5	\n\t"\
		"addpd	%%xmm1,%%xmm7	\n\t"\
		"addpd	%%xmm3,%%xmm4	\n\t"\
		"movaps	%%xmm6,     (%%rsi)	\n\t"\
		"movaps	%%xmm5,0x060(%%rsi)	\n\t"\
		"movaps	%%xmm7,0x010(%%rsi)	\n\t"\
		"movaps	%%xmm4,0x030(%%rsi)	\n\t"\
		"/***************************************/	\n\t"\
		"movslq	%[__p04],%%rdi	\n\t"\
		"shlq	$3,%%rdi	\n\t"\
		"addq	%%rdi,%%rax	/* &a[j1+p5] */	\n\t"\
		"addq	%%rdi,%%rbx	\n\t"\
		"addq	%%rdi,%%rcx	\n\t"\
		"addq	%%rdi,%%rdx	\n\t"\
		"/* SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO_B(r28,c05) */	\n\t"\
		"/* rsi contains r20 */	\n\t"\
		"movaps	    (%%rax),%%xmm0	\n\t"\
		"movaps	    (%%rcx),%%xmm4	\n\t"\
		"movaps	0x10(%%rax),%%xmm1	\n\t"\
		"movaps	0x10(%%rcx),%%xmm5	\n\t"\
		"movaps	0x4f0(%%rsi),%%xmm6	/* c05 */	\n\t"\
		"movaps	0x500(%%rsi),%%xmm7	\n\t"\
		"movaps	%%xmm0,%%xmm2	\n\t"\
		"movaps	%%xmm1,%%xmm3	\n\t"\
		"mulpd	%%xmm6,%%xmm0	\n\t"\
		"mulpd	%%xmm6,%%xmm1	\n\t"\
		"mulpd	%%xmm7,%%xmm2	\n\t"\
		"mulpd	%%xmm7,%%xmm3	\n\t"\
		"movaps	%%xmm4,%%xmm6	\n\t"\
		"addpd	%%xmm2,%%xmm1	\n\t"\
		"movaps	%%xmm5,%%xmm7	\n\t"\
		"mulpd	0x510(%%rsi),%%xmm4	/* c15 */	\n\t"\
		"subpd	%%xmm3,%%xmm0	\n\t"\
		"mulpd	0x510(%%rsi),%%xmm5	\n\t"\
		"mulpd	0x520(%%rsi),%%xmm6	\n\t"\
		"movaps	%%xmm0,%%xmm2	\n\t"\
		"mulpd	0x520(%%rsi),%%xmm7	\n\t"\
		"addpd	%%xmm6,%%xmm5	\n\t"\
		"movaps	%%xmm1,%%xmm3	\n\t"\
		"subpd	%%xmm7,%%xmm4	\n\t"\
		"addpd	%%xmm4,%%xmm0	\n\t"\
		"addpd	%%xmm5,%%xmm1	\n\t"\
		"subpd	%%xmm4,%%xmm2	\n\t"\
		"subpd	%%xmm5,%%xmm3	\n\t"\
		"movaps	    (%%rdx),%%xmm4	\n\t"\
		"movaps	0x10(%%rdx),%%xmm5	\n\t"\
		"movaps	    (%%rdx),%%xmm6	\n\t"\
		"movaps	0x10(%%rdx),%%xmm7	\n\t"\
		"mulpd	0x550(%%rsi),%%xmm4	/* c1D */	\n\t"\
		"mulpd	0x550(%%rsi),%%xmm5	\n\t"\
		"mulpd	0x560(%%rsi),%%xmm6	\n\t"\
		"mulpd	0x560(%%rsi),%%xmm7	\n\t"\
		"addq	$0x80,%%rsi	/* r28 */	\n\t"\
		"addpd	%%xmm6,%%xmm5	\n\t"\
		"subpd	%%xmm7,%%xmm4	\n\t"\
		"movaps	%%xmm5,0x010(%%rsi)	\n\t"\
		"movaps	%%xmm4,     (%%rsi)	\n\t"\
		"movaps	    (%%rbx),%%xmm4	\n\t"\
		"movaps	0x10(%%rbx),%%xmm5	\n\t"\
		"movaps	    (%%rbx),%%xmm6	\n\t"\
		"movaps	0x10(%%rbx),%%xmm7	\n\t"\
		"mulpd	0x4b0(%%rsi),%%xmm4	/* c0D */	\n\t"\
		"mulpd	0x4b0(%%rsi),%%xmm5	\n\t"\
		"mulpd	0x4c0(%%rsi),%%xmm6	\n\t"\
		"mulpd	0x4c0(%%rsi),%%xmm7	\n\t"\
		"addpd	%%xmm6,%%xmm5	\n\t"\
		"subpd	%%xmm7,%%xmm4	\n\t"\
		"movaps	%%xmm5,%%xmm7	\n\t"\
		"movaps	%%xmm4,%%xmm6	\n\t"\
		"subpd	     (%%rsi),%%xmm4	\n\t"\
		"subpd	0x010(%%rsi),%%xmm5	\n\t"\
		"addpd	     (%%rsi),%%xmm6	\n\t"\
		"addpd	0x010(%%rsi),%%xmm7	\n\t"\
		"subpd	%%xmm6,%%xmm0	\n\t"\
		"subpd	%%xmm7,%%xmm1	\n\t"\
		"movaps	%%xmm0,0x040(%%rsi)	\n\t"\
		"subpd	%%xmm5,%%xmm2	\n\t"\
		"movaps	%%xmm1,0x050(%%rsi)	\n\t"\
		"subpd	%%xmm4,%%xmm3	\n\t"\
		"addpd	%%xmm6,%%xmm6	\n\t"\
		"addpd	%%xmm5,%%xmm5	\n\t"\
		"addpd	%%xmm7,%%xmm7	\n\t"\
		"addpd	%%xmm4,%%xmm4	\n\t"\
		"addpd	%%xmm0,%%xmm6	\n\t"\
		"addpd	%%xmm2,%%xmm5	\n\t"\
		"addpd	%%xmm1,%%xmm7	\n\t"\
		"addpd	%%xmm3,%%xmm4	\n\t"\
		"movaps	%%xmm6,     (%%rsi)	\n\t"\
		"movaps	%%xmm7,0x010(%%rsi)	\n\t"\
		"movaps	0x180(%%rsi),%%xmm0	/* isrt2 */	\n\t"\
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
		"movaps	%%xmm2,0x020(%%rsi)	\n\t"\
		"movaps	%%xmm5,0x060(%%rsi)	\n\t"\
		"movaps	%%xmm6,0x030(%%rsi)	\n\t"\
		"movaps	%%xmm7,0x070(%%rsi)	\n\t"\
		"/***************************************/	\n\t"\
		"	\n\t"\
		"/* SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r20) */	\n\t"\
		"subq		$0x80,%%rsi	/* r20 */	\n\t"\
		"movaps		    (%%rsi),%%xmm0	\n\t"\
		"movaps		0x40(%%rsi),%%xmm4	\n\t"\
		"movaps		0x10(%%rsi),%%xmm1	\n\t"\
		"movaps		0x50(%%rsi),%%xmm5	\n\t"\
		"movaps		0x80(%%rsi),%%xmm2	\n\t"\
		"movaps		0xd0(%%rsi),%%xmm7	\n\t"\
		"movaps		0x90(%%rsi),%%xmm3	\n\t"\
		"movaps		0xc0(%%rsi),%%xmm6	\n\t"\
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
		"movaps		%%xmm0,0x80(%%rsi)	\n\t"\
		"movaps		%%xmm4,0x40(%%rsi)	\n\t"\
		"movaps		%%xmm1,0x90(%%rsi)	\n\t"\
		"movaps		%%xmm5,0xd0(%%rsi)	\n\t"\
		"movaps		%%xmm2,    (%%rsi)	\n\t"\
		"movaps		%%xmm7,0xc0(%%rsi)	\n\t"\
		"movaps		%%xmm3,0x10(%%rsi)	\n\t"\
		"movaps		%%xmm6,0x50(%%rsi)	\n\t"\
		"movaps		0x20(%%rsi),%%xmm0	\n\t"\
		"movaps		0x60(%%rsi),%%xmm4	\n\t"\
		"movaps		0x30(%%rsi),%%xmm1	\n\t"\
		"movaps		0x70(%%rsi),%%xmm5	\n\t"\
		"movaps		0xa0(%%rsi),%%xmm2	\n\t"\
		"movaps		0xf0(%%rsi),%%xmm7	\n\t"\
		"movaps		0xb0(%%rsi),%%xmm3	\n\t"\
		"movaps		0xe0(%%rsi),%%xmm6	\n\t"\
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
		"movaps		%%xmm0,0xa0(%%rsi)	\n\t"\
		"movaps		%%xmm4,0x60(%%rsi)	\n\t"\
		"movaps		%%xmm1,0xb0(%%rsi)	\n\t"\
		"movaps		%%xmm5,0xf0(%%rsi)	\n\t"\
		"movaps		%%xmm2,0x20(%%rsi)	\n\t"\
		"movaps		%%xmm7,0xe0(%%rsi)	\n\t"\
		"movaps		%%xmm3,0x30(%%rsi)	\n\t"\
		"movaps		%%xmm6,0x70(%%rsi)	\n\t"\
		"/***************************************/	\n\t"\
		"	\n\t"\
		"/*...Block 4: */	\n\t"\
		"movq	%[__add0],%%rax	\n\t"\
		"movslq	%[__p03],%%rdi	/* Do things this way [rather than repeatedly adding p1] since array-padding scheme means p2 == p+p1 not guaranteed. */	\n\t"\
		"movslq	%[__p08],%%rbx	/* Can't get these via simple load-one-and-shift-as-needed due to array padding scheme */	\n\t"\
		"movslq	%[__p10],%%rcx	\n\t"\
		"movslq	%[__p18],%%rdx	\n\t"\
		"shlq	$3,%%rdi	\n\t"\
		"shlq	$3,%%rbx	\n\t"\
		"shlq	$3,%%rcx	\n\t"\
		"shlq	$3,%%rdx	\n\t"\
		"addq	%%rdi,%%rax	/* &a[j1+p3] */	\n\t"\
		"addq	%%rax,%%rbx	\n\t"\
		"addq	%%rax,%%rcx	\n\t"\
		"addq	%%rax,%%rdx	\n\t"\
		"/* SSE2_RADIX4_DIF_4TWIDDLE_B(r30,c03) */	\n\t"\
		"addq	$0x100,%%rsi	/* r30 */	\n\t"\
		"movaps	    (%%rax),%%xmm0	\n\t"\
		"movaps	    (%%rcx),%%xmm4	\n\t"\
		"movaps	0x10(%%rax),%%xmm1	\n\t"\
		"movaps	0x10(%%rcx),%%xmm5	\n\t"\
		"movaps	0x470(%%rsi),%%xmm6	/* c03 */	\n\t"\
		"movaps	0x480(%%rsi),%%xmm7	\n\t"\
		"movaps	%%xmm0,%%xmm2	\n\t"\
		"movaps	%%xmm1,%%xmm3	\n\t"\
		"mulpd   %%xmm6,%%xmm0	\n\t"\
		"mulpd   %%xmm6,%%xmm1	\n\t"\
		"mulpd   %%xmm7,%%xmm2	\n\t"\
		"mulpd   %%xmm7,%%xmm3	\n\t"\
		"movaps	%%xmm4,%%xmm6	\n\t"\
		"addpd   %%xmm2,%%xmm1	\n\t"\
		"movaps	%%xmm5,%%xmm7	\n\t"\
		"mulpd   0x490(%%rsi),%%xmm4	/* c13 */	\n\t"\
		"subpd   %%xmm3,%%xmm0	\n\t"\
		"mulpd   0x490(%%rsi),%%xmm5	\n\t"\
		"mulpd   0x4a0(%%rsi),%%xmm6	\n\t"\
		"movaps	%%xmm0,%%xmm2	\n\t"\
		"mulpd	0x4a0(%%rsi),%%xmm7	\n\t"\
		"addpd	%%xmm6,%%xmm5	\n\t"\
		"movaps	%%xmm1,%%xmm3	\n\t"\
		"subpd	%%xmm7,%%xmm4	\n\t"\
		"addpd	%%xmm4,%%xmm0	\n\t"\
		"addpd	%%xmm5,%%xmm1	\n\t"\
		"subpd	%%xmm4,%%xmm2	\n\t"\
		"subpd	%%xmm5,%%xmm3	\n\t"\
		"movaps	    (%%rdx),%%xmm4	\n\t"\
		"movaps	0x10(%%rdx),%%xmm5	\n\t"\
		"movaps	    (%%rdx),%%xmm6	\n\t"\
		"movaps	0x10(%%rdx),%%xmm7	\n\t"\
		"mulpd	0x4d0(%%rsi),%%xmm4	/* c1B */	\n\t"\
		"mulpd	0x4d0(%%rsi),%%xmm5	\n\t"\
		"mulpd	0x4e0(%%rsi),%%xmm6	\n\t"\
		"mulpd	0x4e0(%%rsi),%%xmm7	\n\t"\
		"addpd	%%xmm6,%%xmm5	\n\t"\
		"subpd	%%xmm7,%%xmm4	\n\t"\
		"movaps	%%xmm5,0x010(%%rsi)	\n\t"\
		"movaps	%%xmm4,     (%%rsi)	\n\t"\
		"movaps	    (%%rbx),%%xmm4	\n\t"\
		"movaps	0x10(%%rbx),%%xmm5	\n\t"\
		"movaps	    (%%rbx),%%xmm6	\n\t"\
		"movaps	0x10(%%rbx),%%xmm7	\n\t"\
		"mulpd	0x4b0(%%rsi),%%xmm4	/* c0B */	\n\t"\
		"mulpd	0x4b0(%%rsi),%%xmm5	\n\t"\
		"mulpd	0x4c0(%%rsi),%%xmm6	\n\t"\
		"mulpd	0x4c0(%%rsi),%%xmm7	\n\t"\
		"addpd	%%xmm6,%%xmm5	\n\t"\
		"subpd	%%xmm7,%%xmm4	\n\t"\
		"movaps	%%xmm5,%%xmm7	\n\t"\
		"movaps	%%xmm4,%%xmm6	\n\t"\
		"subpd	    (%%rsi),%%xmm4	\n\t"\
		"subpd	0x10(%%rsi),%%xmm5	\n\t"\
		"addpd	    (%%rsi),%%xmm6	\n\t"\
		"addpd	0x10(%%rsi),%%xmm7	\n\t"\
		"subpd	%%xmm6,%%xmm0	\n\t"\
		"subpd	%%xmm5,%%xmm2	\n\t"\
		"subpd	%%xmm7,%%xmm1	\n\t"\
		"subpd	%%xmm4,%%xmm3	\n\t"\
		"movaps	%%xmm0,0x040(%%rsi)	\n\t"\
		"movaps	%%xmm2,0x020(%%rsi)	\n\t"\
		"movaps	%%xmm1,0x050(%%rsi)	\n\t"\
		"movaps	%%xmm3,0x070(%%rsi)	\n\t"\
		"addpd	%%xmm6,%%xmm6	\n\t"\
		"addpd	%%xmm5,%%xmm5	\n\t"\
		"addpd	%%xmm7,%%xmm7	\n\t"\
		"addpd	%%xmm4,%%xmm4	\n\t"\
		"addpd	%%xmm0,%%xmm6	\n\t"\
		"addpd	%%xmm2,%%xmm5	\n\t"\
		"addpd	%%xmm1,%%xmm7	\n\t"\
		"addpd	%%xmm3,%%xmm4	\n\t"\
		"movaps	%%xmm6,     (%%rsi)	\n\t"\
		"movaps	%%xmm5,0x060(%%rsi)	\n\t"\
		"movaps	%%xmm7,0x010(%%rsi)	\n\t"\
		"movaps	%%xmm4,0x030(%%rsi)	\n\t"\
		"/***************************************/	\n\t"\
		"movslq	%[__p04],%%rdi	\n\t"\
		"shlq	$3,%%rdi	\n\t"\
		"addq	%%rdi,%%rax	/* &a[j1+p7] */	\n\t"\
		"addq	%%rdi,%%rbx	\n\t"\
		"addq	%%rdi,%%rcx	\n\t"\
		"addq	%%rdi,%%rdx	\n\t"\
		"/* SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO_B(r38,c07) */	\n\t"\
		"/* rsi contains r30 */	\n\t"\
		"movaps	    (%%rax),%%xmm0	\n\t"\
		"movaps	    (%%rcx),%%xmm4	\n\t"\
		"movaps	0x10(%%rax),%%xmm1	\n\t"\
		"movaps	0x10(%%rcx),%%xmm5	\n\t"\
		"movaps	0x4f0(%%rsi),%%xmm6	/* c07 */	\n\t"\
		"movaps	0x500(%%rsi),%%xmm7	\n\t"\
		"movaps	%%xmm0,%%xmm2	\n\t"\
		"movaps	%%xmm1,%%xmm3	\n\t"\
		"mulpd	%%xmm6,%%xmm0	\n\t"\
		"mulpd	%%xmm6,%%xmm1	\n\t"\
		"mulpd	%%xmm7,%%xmm2	\n\t"\
		"mulpd	%%xmm7,%%xmm3	\n\t"\
		"movaps	%%xmm4,%%xmm6	\n\t"\
		"addpd	%%xmm2,%%xmm1	\n\t"\
		"movaps	%%xmm5,%%xmm7	\n\t"\
		"mulpd	0x510(%%rsi),%%xmm4	/* c17 */	\n\t"\
		"subpd	%%xmm3,%%xmm0	\n\t"\
		"mulpd	0x510(%%rsi),%%xmm5	\n\t"\
		"mulpd	0x520(%%rsi),%%xmm6	\n\t"\
		"movaps	%%xmm0,%%xmm2	\n\t"\
		"mulpd	0x520(%%rsi),%%xmm7	\n\t"\
		"addpd	%%xmm6,%%xmm5	\n\t"\
		"movaps	%%xmm1,%%xmm3	\n\t"\
		"subpd	%%xmm7,%%xmm4	\n\t"\
		"addpd	%%xmm4,%%xmm0	\n\t"\
		"addpd	%%xmm5,%%xmm1	\n\t"\
		"subpd	%%xmm4,%%xmm2	\n\t"\
		"subpd	%%xmm5,%%xmm3	\n\t"\
		"movaps	    (%%rdx),%%xmm4	\n\t"\
		"movaps	0x10(%%rdx),%%xmm5	\n\t"\
		"movaps	    (%%rdx),%%xmm6	\n\t"\
		"movaps	0x10(%%rdx),%%xmm7	\n\t"\
		"mulpd	0x550(%%rsi),%%xmm4	/* c1F */	\n\t"\
		"mulpd	0x550(%%rsi),%%xmm5	\n\t"\
		"mulpd	0x560(%%rsi),%%xmm6	\n\t"\
		"mulpd	0x560(%%rsi),%%xmm7	\n\t"\
		"addq	$0x80,%%rsi	/* r38 */	\n\t"\
		"addpd	%%xmm6,%%xmm5	\n\t"\
		"subpd	%%xmm7,%%xmm4	\n\t"\
		"movaps	%%xmm5,0x010(%%rsi)	\n\t"\
		"movaps	%%xmm4,     (%%rsi)	\n\t"\
		"movaps	    (%%rbx),%%xmm4	\n\t"\
		"movaps	0x10(%%rbx),%%xmm5	\n\t"\
		"movaps	    (%%rbx),%%xmm6	\n\t"\
		"movaps	0x10(%%rbx),%%xmm7	\n\t"\
		"mulpd	0x4b0(%%rsi),%%xmm4	/* c0F */	\n\t"\
		"mulpd	0x4b0(%%rsi),%%xmm5	\n\t"\
		"mulpd	0x4c0(%%rsi),%%xmm6	\n\t"\
		"mulpd	0x4c0(%%rsi),%%xmm7	\n\t"\
		"addpd	%%xmm6,%%xmm5	\n\t"\
		"subpd	%%xmm7,%%xmm4	\n\t"\
		"movaps	%%xmm5,%%xmm7	\n\t"\
		"movaps	%%xmm4,%%xmm6	\n\t"\
		"subpd	     (%%rsi),%%xmm4	\n\t"\
		"subpd	0x010(%%rsi),%%xmm5	\n\t"\
		"addpd	     (%%rsi),%%xmm6	\n\t"\
		"addpd	0x010(%%rsi),%%xmm7	\n\t"\
		"subpd	%%xmm6,%%xmm0	\n\t"\
		"subpd	%%xmm7,%%xmm1	\n\t"\
		"movaps	%%xmm0,0x040(%%rsi)	\n\t"\
		"subpd	%%xmm5,%%xmm2	\n\t"\
		"movaps	%%xmm1,0x050(%%rsi)	\n\t"\
		"subpd	%%xmm4,%%xmm3	\n\t"\
		"addpd	%%xmm6,%%xmm6	\n\t"\
		"addpd	%%xmm5,%%xmm5	\n\t"\
		"addpd	%%xmm7,%%xmm7	\n\t"\
		"addpd	%%xmm4,%%xmm4	\n\t"\
		"addpd	%%xmm0,%%xmm6	\n\t"\
		"addpd	%%xmm2,%%xmm5	\n\t"\
		"addpd	%%xmm1,%%xmm7	\n\t"\
		"addpd	%%xmm3,%%xmm4	\n\t"\
		"movaps	%%xmm6,     (%%rsi)	\n\t"\
		"movaps	%%xmm7,0x010(%%rsi)	\n\t"\
		"movaps	0x080(%%rsi),%%xmm0	/* isrt2 */	\n\t"\
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
		"movaps	%%xmm2,0x020(%%rsi)	\n\t"\
		"movaps	%%xmm5,0x060(%%rsi)	\n\t"\
		"movaps	%%xmm6,0x030(%%rsi)	\n\t"\
		"movaps	%%xmm7,0x070(%%rsi)	\n\t"\
		"/***************************************/	\n\t"\
		"	\n\t"\
		"/* SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r30) */	\n\t"\
		"subq		$0x80,%%rsi	/* r30 */	\n\t"\
		"movaps		    (%%rsi),%%xmm0	\n\t"\
		"movaps		0x40(%%rsi),%%xmm4	\n\t"\
		"movaps		0x10(%%rsi),%%xmm1	\n\t"\
		"movaps		0x50(%%rsi),%%xmm5	\n\t"\
		"movaps		0x80(%%rsi),%%xmm2	\n\t"\
		"movaps		0xd0(%%rsi),%%xmm7	\n\t"\
		"movaps		0x90(%%rsi),%%xmm3	\n\t"\
		"movaps		0xc0(%%rsi),%%xmm6	\n\t"\
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
		"movaps		%%xmm0,0x80(%%rsi)	\n\t"\
		"movaps		%%xmm4,0x40(%%rsi)	\n\t"\
		"movaps		%%xmm1,0x90(%%rsi)	\n\t"\
		"movaps		%%xmm5,0xd0(%%rsi)	\n\t"\
		"movaps		%%xmm2,    (%%rsi)	\n\t"\
		"movaps		%%xmm7,0xc0(%%rsi)	\n\t"\
		"movaps		%%xmm3,0x10(%%rsi)	\n\t"\
		"movaps		%%xmm6,0x50(%%rsi)	\n\t"\
		"movaps		0x20(%%rsi),%%xmm0	\n\t"\
		"movaps		0x60(%%rsi),%%xmm4	\n\t"\
		"movaps		0x30(%%rsi),%%xmm1	\n\t"\
		"movaps		0x70(%%rsi),%%xmm5	\n\t"\
		"movaps		0xa0(%%rsi),%%xmm2	\n\t"\
		"movaps		0xf0(%%rsi),%%xmm7	\n\t"\
		"movaps		0xb0(%%rsi),%%xmm3	\n\t"\
		"movaps		0xe0(%%rsi),%%xmm6	\n\t"\
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
		"movaps		%%xmm0,0xa0(%%rsi)	\n\t"\
		"movaps		%%xmm4,0x60(%%rsi)	\n\t"\
		"movaps		%%xmm1,0xb0(%%rsi)	\n\t"\
		"movaps		%%xmm5,0xf0(%%rsi)	\n\t"\
		"movaps		%%xmm2,0x20(%%rsi)	\n\t"\
		"movaps		%%xmm7,0xe0(%%rsi)	\n\t"\
		"movaps		%%xmm3,0x30(%%rsi)	\n\t"\
		"movaps		%%xmm6,0x70(%%rsi)	\n\t"\
		"	\n\t"\
		"/**********************************************************************************/	\n\t"\
		"/*...and now do eight radix-4 transforms, including the internal twiddle factors: */	\n\t"\
		"/**********************************************************************************/	\n\t"\
		"/*...Block 1: t00,t10,t20,t30	*/	\n\t"\
		"movq	%[__add0],%%rax	/* &a[j1] */	\n\t"\
		"movslq	%[__p01],%%rbx	\n\t"\
		"movslq	%[__p02],%%rcx	\n\t"\
		"movslq	%[__p03],%%rdx	\n\t"\
		"shlq	$3,%%rbx	\n\t"\
		"shlq	$3,%%rcx	\n\t"\
		"shlq	$3,%%rdx	\n\t"\
		"addq	%%rax,%%rbx	\n\t"\
		"addq	%%rax,%%rcx	\n\t"\
		"addq	%%rax,%%rdx	\n\t"\
		"movq	%[__r00],%%rsi	\n\t"\
		"movaps	     (%%rsi),%%xmm0	\n\t"\
		"movaps	0x200(%%rsi),%%xmm4	\n\t"\
		"movaps	0x010(%%rsi),%%xmm1	\n\t"\
		"movaps	0x210(%%rsi),%%xmm5	\n\t"\
		"movaps	0x100(%%rsi),%%xmm2	\n\t"\
		"movaps	0x300(%%rsi),%%xmm6	\n\t"\
		"movaps	0x110(%%rsi),%%xmm3	\n\t"\
		"movaps	0x310(%%rsi),%%xmm7	\n\t"\
		"subpd	0x100(%%rsi),%%xmm0	\n\t"\
		"subpd	0x300(%%rsi),%%xmm4	\n\t"\
		"subpd	0x110(%%rsi),%%xmm1	\n\t"\
		"subpd	0x310(%%rsi),%%xmm5	\n\t"\
		"addpd	     (%%rsi),%%xmm2	\n\t"\
		"addpd	0x200(%%rsi),%%xmm6	\n\t"\
		"addpd	0x010(%%rsi),%%xmm3	\n\t"\
		"addpd	0x210(%%rsi),%%xmm7	\n\t"\
		"subpd	%%xmm6,%%xmm2	\n\t"\
		"subpd	%%xmm5,%%xmm0	\n\t"\
		"subpd	%%xmm7,%%xmm3	\n\t"\
		"subpd	%%xmm4,%%xmm1	\n\t"\
		"addpd	%%xmm6,%%xmm6	\n\t"\
		"addpd	%%xmm5,%%xmm5	\n\t"\
		"addpd	%%xmm7,%%xmm7	\n\t"\
		"addpd	%%xmm4,%%xmm4	\n\t"\
		"movaps	%%xmm2,    (%%rbx)	\n\t"\
		"movaps	%%xmm0,    (%%rcx)	\n\t"\
		"movaps	%%xmm3,0x10(%%rbx)	\n\t"\
		"movaps	%%xmm1,0x10(%%rdx)	\n\t"\
		"addpd	%%xmm2,%%xmm6	\n\t"\
		"addpd	%%xmm0,%%xmm5	\n\t"\
		"addpd	%%xmm3,%%xmm7	\n\t"\
		"addpd	%%xmm1,%%xmm4	\n\t"\
		"movaps	%%xmm6,    (%%rax)	\n\t"\
		"movaps	%%xmm5,    (%%rdx)	\n\t"\
		"movaps	%%xmm7,0x10(%%rax)	\n\t"\
		"movaps	%%xmm4,0x10(%%rcx)	\n\t"\
		"	\n\t"\
		"/*...Block 5: t08,t18,t28,t38	*/	\n\t"\
		"addq	$0x80,%%rsi	/* r08 */	\n\t"\
		"subq	%%rax,%%rbx	/* p01 << 3 */	\n\t"\
		"subq	%%rax,%%rcx	/* p02 << 3 */	\n\t"\
		"subq	%%rax,%%rdx	/* p03 << 3 */	\n\t"\
		"movslq	%[__p04],%%rdi	\n\t"\
		"shlq	$3,%%rdi	\n\t"\
		"addq	%%rdi,%%rax	/* &a[j1+p04] */	\n\t"\
		"addq	%%rax,%%rbx	\n\t"\
		"addq	%%rax,%%rcx	\n\t"\
		"addq	%%rax,%%rdx	\n\t"\
		"movaps	0x380(%%rsi),%%xmm3	/* isrt2 */	\n\t"\
		"movaps	0x200(%%rsi),%%xmm4	\n\t"\
		"movaps	0x210(%%rsi),%%xmm5	\n\t"\
		"movaps	0x300(%%rsi),%%xmm6	\n\t"\
		"movaps	0x310(%%rsi),%%xmm7	\n\t"\
		"mulpd	%%xmm3,%%xmm4	\n\t"\
		"movaps	     (%%rsi),%%xmm0	\n\t"\
		"mulpd	%%xmm3,%%xmm5	\n\t"\
		"movaps	0x010(%%rsi),%%xmm1	\n\t"\
		"mulpd	%%xmm3,%%xmm6	\n\t"\
		"movaps	0x100(%%rsi),%%xmm2	\n\t"\
		"mulpd	%%xmm3,%%xmm7	\n\t"\
		"movaps	0x110(%%rsi),%%xmm3	\n\t"\
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
		"movaps	%%xmm0,    (%%rbx)	\n\t"\
		"movaps	%%xmm3,    (%%rcx)	\n\t"\
		"movaps	%%xmm2,0x10(%%rbx)	\n\t"\
		"movaps	%%xmm1,0x10(%%rdx)	\n\t"\
		"addpd	%%xmm0,%%xmm4	\n\t"\
		"addpd	%%xmm3,%%xmm7	\n\t"\
		"addpd	%%xmm2,%%xmm5	\n\t"\
		"addpd	%%xmm1,%%xmm6	\n\t"\
		"movaps	%%xmm4,    (%%rax)	\n\t"\
		"movaps	%%xmm7,    (%%rdx)	\n\t"\
		"movaps	%%xmm5,0x10(%%rax)	\n\t"\
		"movaps	%%xmm6,0x10(%%rcx)	\n\t"\
		"	\n\t"\
		"/*...Block 3: t04,t14,t24,t34	*/	\n\t"\
		"subq	$0x40,%%rsi	/* r04 */	\n\t"\
		"subq	%%rax,%%rbx	/* p01 << 3 */	\n\t"\
		"subq	%%rax,%%rcx	/* p02 << 3 */	\n\t"\
		"subq	%%rax,%%rdx	/* p03 << 3 */	\n\t"\
		"subq	%%rdi,%%rax	/* &a[j1] */	\n\t"\
		"movslq	%[__p08],%%rdi	\n\t"\
		"shlq	$3,%%rdi	\n\t"\
		"addq	%%rdi,%%rax	/* &a[j1+p08] */	\n\t"\
		"addq	%%rax,%%rbx	\n\t"\
		"addq	%%rax,%%rcx	\n\t"\
		"addq	%%rax,%%rdx	\n\t"\
		"movaps	0x200(%%rsi),%%xmm4	\n\t"\
		"movaps	0x210(%%rsi),%%xmm5	\n\t"\
		"movaps	0x3d0(%%rsi),%%xmm3	/* cc0 */	\n\t"\
		"movaps	0x3e0(%%rsi),%%xmm2	\n\t"\
		"movaps	%%xmm4,%%xmm6	\n\t"\
		"movaps	%%xmm5,%%xmm7	\n\t"\
		"mulpd	%%xmm3,%%xmm4	\n\t"\
		"mulpd	%%xmm3,%%xmm5	\n\t"\
		"mulpd	%%xmm2,%%xmm6	\n\t"\
		"movaps	0x300(%%rsi),%%xmm0	\n\t"\
		"mulpd	%%xmm2,%%xmm7	\n\t"\
		"movaps	0x310(%%rsi),%%xmm1	\n\t"\
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
		"movaps	0x100(%%rsi),%%xmm2	\n\t"\
		"movaps	0x110(%%rsi),%%xmm3	\n\t"\
		"movaps	0x3c0(%%rsi),%%xmm1	/* isrt2 */	\n\t"\
		"movaps	%%xmm2,%%xmm0	\n\t"\
		"subpd	%%xmm3,%%xmm2	\n\t"\
		"addpd	%%xmm0,%%xmm3	\n\t"\
		"mulpd	%%xmm1,%%xmm2	\n\t"\
		"mulpd	%%xmm1,%%xmm3	\n\t"\
		"movaps	     (%%rsi),%%xmm0	\n\t"\
		"movaps	0x010(%%rsi),%%xmm1	\n\t"\
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
		"movaps	%%xmm2,    (%%rbx)	\n\t"\
		"movaps	%%xmm0,    (%%rcx)	\n\t"\
		"movaps	%%xmm3,0x10(%%rbx)	\n\t"\
		"movaps	%%xmm1,0x10(%%rdx)	\n\t"\
		"addpd	%%xmm2,%%xmm6	\n\t"\
		"addpd	%%xmm0,%%xmm5	\n\t"\
		"addpd	%%xmm3,%%xmm7	\n\t"\
		"addpd	%%xmm1,%%xmm4	\n\t"\
		"movaps	%%xmm6,    (%%rax)	\n\t"\
		"movaps	%%xmm5,    (%%rdx)	\n\t"\
		"movaps	%%xmm7,0x10(%%rax)	\n\t"\
		"movaps	%%xmm4,0x10(%%rcx)	\n\t"\
		"	\n\t"\
		"/*...Block 7: t0C,t1C,t2C,t3C	*/	\n\t"\
		"addq	$0x80,%%rsi	/* r0C */	\n\t"\
		"subq	%%rax,%%rbx	/* p01 << 3 */	\n\t"\
		"subq	%%rax,%%rcx	/* p02 << 3 */	\n\t"\
		"subq	%%rax,%%rdx	/* p03 << 3 */	\n\t"\
		"subq	%%rdi,%%rax	/* &a[j1] */	\n\t"\
		"movslq	%[__p0C],%%rdi	\n\t"\
		"shlq	$3,%%rdi	\n\t"\
		"addq	%%rdi,%%rax	/* &a[j1+p0C] */	\n\t"\
		"addq	%%rax,%%rbx	\n\t"\
		"addq	%%rax,%%rcx	\n\t"\
		"addq	%%rax,%%rdx	\n\t"\
		"movaps	0x200(%%rsi),%%xmm4	\n\t"\
		"movaps	0x210(%%rsi),%%xmm5	\n\t"\
		"movaps	0x350(%%rsi),%%xmm2	/* cc0 */	\n\t"\
		"movaps	0x360(%%rsi),%%xmm3	\n\t"\
		"movaps	%%xmm4,%%xmm6	\n\t"\
		"movaps	%%xmm5,%%xmm7	\n\t"\
		"mulpd	%%xmm3,%%xmm4	\n\t"\
		"mulpd	%%xmm3,%%xmm5	\n\t"\
		"mulpd	%%xmm2,%%xmm6	\n\t"\
		"movaps	0x300(%%rsi),%%xmm0	\n\t"\
		"mulpd	%%xmm2,%%xmm7	\n\t"\
		"movaps	0x310(%%rsi),%%xmm1	\n\t"\
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
		"movaps	0x100(%%rsi),%%xmm2	\n\t"\
		"movaps	0x110(%%rsi),%%xmm3	\n\t"\
		"movaps	0x340(%%rsi),%%xmm1	/* isrt2 */	\n\t"\
		"movaps	%%xmm2,%%xmm0	\n\t"\
		"addpd	%%xmm3,%%xmm2	\n\t"\
		"subpd	%%xmm0,%%xmm3	\n\t"\
		"mulpd	%%xmm1,%%xmm2	\n\t"\
		"mulpd	%%xmm1,%%xmm3	\n\t"\
		"movaps	     (%%rsi),%%xmm0	\n\t"\
		"movaps	0x010(%%rsi),%%xmm1	\n\t"\
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
		"movaps	%%xmm0,    (%%rbx)	\n\t"\
		"movaps	%%xmm2,    (%%rcx)	\n\t"\
		"movaps	%%xmm1,0x10(%%rbx)	\n\t"\
		"movaps	%%xmm3,0x10(%%rdx)	\n\t"\
		"addpd	%%xmm0,%%xmm4	\n\t"\
		"addpd	%%xmm2,%%xmm7	\n\t"\
		"addpd	%%xmm1,%%xmm5	\n\t"\
		"addpd	%%xmm3,%%xmm6	\n\t"\
		"movaps	%%xmm4,    (%%rax)	\n\t"\
		"movaps	%%xmm7,    (%%rdx)	\n\t"\
		"movaps	%%xmm5,0x10(%%rax)	\n\t"\
		"movaps	%%xmm6,0x10(%%rcx)	\n\t"\
		"	\n\t"\
		"/*...Block 2: t02,t12,t22,t32	*/	\n\t"\
		"subq	$0xa0,%%rsi	/* r02 */	\n\t"\
		"movslq	%[__p10],%%rdi	\n\t"\
		"subq	%%rax,%%rbx	\n\t"\
		"subq	%%rax,%%rcx	\n\t"\
		"subq	%%rax,%%rdx	\n\t"\
		"movq	%[__add0],%%rax	/* &a[j1] */	\n\t"\
		"shlq	$3,%%rdi	\n\t"\
		"addq	%%rdi,%%rax	/* &a[j1+p10) */	\n\t"\
		"addq	%%rax,%%rbx	\n\t"\
		"addq	%%rax,%%rcx	\n\t"\
		"addq	%%rax,%%rdx	\n\t"\
		"movaps	0x200(%%rsi),%%xmm4	\n\t"\
		"movaps	0x210(%%rsi),%%xmm5	\n\t"\
		"movaps	0x410(%%rsi),%%xmm2	/* cc1 */	\n\t"\
		"movaps	0x420(%%rsi),%%xmm3	\n\t"\
		"movaps	%%xmm4,%%xmm6	\n\t"\
		"movaps	%%xmm5,%%xmm7	\n\t"\
		"mulpd	%%xmm2,%%xmm4	\n\t"\
		"mulpd	%%xmm2,%%xmm5	\n\t"\
		"mulpd	%%xmm3,%%xmm6	\n\t"\
		"movaps	0x300(%%rsi),%%xmm0	\n\t"\
		"mulpd	%%xmm3,%%xmm7	\n\t"\
		"movaps	0x310(%%rsi),%%xmm1	\n\t"\
		"movaps	0x430(%%rsi),%%xmm2	/* cc3 */	\n\t"\
		"movaps	0x440(%%rsi),%%xmm3	\n\t"\
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
		"movaps	0x100(%%rsi),%%xmm1	\n\t"\
		"movaps	0x110(%%rsi),%%xmm3	\n\t"\
		"movaps	0x400(%%rsi),%%xmm0	/* ss1 */	\n\t"\
		"movaps	%%xmm1,%%xmm2	\n\t"\
		"mulpd	%%xmm0,%%xmm1	\n\t"\
		"mulpd	%%xmm3,%%xmm0	\n\t"\
		"mulpd	0x3f0(%%rsi),%%xmm2		/* cc1 */	\n\t"\
		"mulpd	0x3f0(%%rsi),%%xmm3	\n\t"\
		"subpd	%%xmm0,%%xmm2	\n\t"\
		"addpd	%%xmm1,%%xmm3	\n\t"\
		"movaps	     (%%rsi),%%xmm0	\n\t"\
		"movaps	0x010(%%rsi),%%xmm1	\n\t"\
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
		"movaps	%%xmm2,    (%%rbx)	\n\t"\
		"movaps	%%xmm0,    (%%rcx)	\n\t"\
		"movaps	%%xmm3,0x10(%%rbx)	\n\t"\
		"movaps	%%xmm1,0x10(%%rdx)	\n\t"\
		"addpd	%%xmm2,%%xmm6	\n\t"\
		"addpd	%%xmm0,%%xmm5	\n\t"\
		"addpd	%%xmm3,%%xmm7	\n\t"\
		"addpd	%%xmm1,%%xmm4	\n\t"\
		"movaps	%%xmm6,    (%%rax)	\n\t"\
		"movaps	%%xmm5,    (%%rdx)	\n\t"\
		"movaps	%%xmm7,0x10(%%rax)	\n\t"\
		"movaps	%%xmm4,0x10(%%rcx)	\n\t"\
		"	\n\t"\
		"/*...Block 6: t0A,t1A,t2A,t3A	*/	\n\t"\
		"addq	$0x80,%%rsi	/* r0A */	\n\t"\
		"movslq	%[__p04],%%rdi	\n\t"\
		"shlq	$3,%%rdi	\n\t"\
		"addq	%%rdi,%%rax	/* &a[j1+p14] */	\n\t"\
		"addq	%%rdi,%%rbx	\n\t"\
		"addq	%%rdi,%%rcx	\n\t"\
		"addq	%%rdi,%%rdx	\n\t"\
		"movaps	0x200(%%rsi),%%xmm4	\n\t"\
		"movaps	0x210(%%rsi),%%xmm5	\n\t"\
		"movaps	0x3b0(%%rsi),%%xmm3	/* cc3 */	\n\t"\
		"movaps	0x3c0(%%rsi),%%xmm2	\n\t"\
		"movaps	%%xmm4,%%xmm6	\n\t"\
		"movaps	%%xmm5,%%xmm7	\n\t"\
		"mulpd	%%xmm2,%%xmm4	\n\t"\
		"mulpd	%%xmm2,%%xmm5	\n\t"\
		"mulpd	%%xmm3,%%xmm6	\n\t"\
		"movaps	0x300(%%rsi),%%xmm0	\n\t"\
		"mulpd	%%xmm3,%%xmm7	\n\t"\
		"movaps	0x310(%%rsi),%%xmm1	\n\t"\
		"movaps	0x390(%%rsi),%%xmm2	/* cc1 */	\n\t"\
		"movaps	0x3a0(%%rsi),%%xmm3	\n\t"\
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
		"movaps	0x100(%%rsi),%%xmm1	\n\t"\
		"movaps	0x110(%%rsi),%%xmm3	\n\t"\
		"movaps	0x370(%%rsi),%%xmm0		/* cc0 */	\n\t"\
		"movaps	%%xmm1,%%xmm2	\n\t"\
		"mulpd	%%xmm0,%%xmm1	\n\t"\
		"mulpd	%%xmm3,%%xmm0	\n\t"\
		"mulpd	0x380(%%rsi),%%xmm2	/* ss0 */	\n\t"\
		"mulpd	0x380(%%rsi),%%xmm3	\n\t"\
		"addpd	%%xmm0,%%xmm2	\n\t"\
		"subpd	%%xmm1,%%xmm3	\n\t"\
		"movaps	     (%%rsi),%%xmm0	\n\t"\
		"movaps	0x010(%%rsi),%%xmm1	\n\t"\
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
		"movaps	%%xmm0,    (%%rbx)	\n\t"\
		"movaps	%%xmm2,    (%%rcx)	\n\t"\
		"movaps	%%xmm1,0x10(%%rbx)	\n\t"\
		"movaps	%%xmm3,0x10(%%rdx)	\n\t"\
		"addpd	%%xmm0,	%%xmm4	\n\t"\
		"addpd	%%xmm2,	%%xmm7	\n\t"\
		"addpd	%%xmm1,	%%xmm5	\n\t"\
		"addpd	%%xmm3,	%%xmm6	\n\t"\
		"movaps	%%xmm4,    (%%rax)	\n\t"\
		"movaps	%%xmm7,    (%%rdx)	\n\t"\
		"movaps	%%xmm5,0x10(%%rax)	\n\t"\
		"movaps	%%xmm6,0x10(%%rcx)	\n\t"\
		"	\n\t"\
		"/*...Block 4: t06,t16,t26,t36	*/	\n\t"\
		"subq	$0x40,%%rsi	/* r06 */	\n\t"\
		"movslq	%[__p18],%%rdi	\n\t"\
		"subq	%%rax,%%rbx	\n\t"\
		"subq	%%rax,%%rcx	\n\t"\
		"subq	%%rax,%%rdx	\n\t"\
		"movq	%[__add0],%%rax	/* &a[j1] */	\n\t"\
		"shlq	$3,%%rdi	\n\t"\
		"addq	%%rdi,%%rax	/* &a[j1+p18] */	\n\t"\
		"addq	%%rax,%%rbx	\n\t"\
		"addq	%%rax,%%rcx	\n\t"\
		"addq	%%rax,%%rdx	\n\t"\
		"movaps	0x200(%%rsi),%%xmm4	\n\t"\
		"movaps	0x210(%%rsi),%%xmm5	\n\t"\
		"movaps	0x3f0(%%rsi),%%xmm2	/* cc3 */	\n\t"\
		"movaps	0x400(%%rsi),%%xmm3	\n\t"\
		"movaps	%%xmm4,%%xmm6	\n\t"\
		"movaps	%%xmm5,%%xmm7	\n\t"\
		"mulpd	%%xmm2,%%xmm4	\n\t"\
		"mulpd	%%xmm2,%%xmm5	\n\t"\
		"mulpd	%%xmm3,%%xmm6	\n\t"\
		"movaps	0x300(%%rsi),%%xmm0	\n\t"\
		"mulpd	%%xmm3,%%xmm7	\n\t"\
		"movaps	0x310(%%rsi),%%xmm1	\n\t"\
		"movaps	0x3d0(%%rsi),%%xmm3	/* cc1 */	\n\t"\
		"movaps	0x3e0(%%rsi),%%xmm2	\n\t"\
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
		"movaps	0x100(%%rsi),%%xmm1	\n\t"\
		"movaps	0x110(%%rsi),%%xmm3	\n\t"\
		"movaps	0x3b0(%%rsi),%%xmm0		/* cc0 */	\n\t"\
		"movaps	%%xmm1,%%xmm2	\n\t"\
		"mulpd	%%xmm0,%%xmm1	\n\t"\
		"mulpd	%%xmm3,%%xmm0	\n\t"\
		"mulpd	0x3c0(%%rsi),%%xmm2		/* ss0 */	\n\t"\
		"mulpd	0x3c0(%%rsi),%%xmm3	\n\t"\
		"subpd	%%xmm0,%%xmm2	\n\t"\
		"addpd	%%xmm1,%%xmm3	\n\t"\
		"movaps	     (%%rsi),%%xmm0	\n\t"\
		"movaps	0x010(%%rsi),%%xmm1	\n\t"\
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
		"movaps	%%xmm2,    (%%rbx)	\n\t"\
		"movaps	%%xmm0,    (%%rcx)	\n\t"\
		"movaps	%%xmm3,0x10(%%rbx)	\n\t"\
		"movaps	%%xmm1,0x10(%%rdx)	\n\t"\
		"addpd	%%xmm2,	%%xmm4	\n\t"\
		"addpd	%%xmm0,	%%xmm7	\n\t"\
		"addpd	%%xmm3,	%%xmm5	\n\t"\
		"addpd	%%xmm1,	%%xmm6	\n\t"\
		"movaps	%%xmm4,    (%%rax)	\n\t"\
		"movaps	%%xmm7,    (%%rdx)	\n\t"\
		"movaps	%%xmm5,0x10(%%rax)	\n\t"\
		"movaps	%%xmm6,0x10(%%rcx)	\n\t"\
		"	\n\t"\
		"/*...Block 8: t0E,t1E,t2E,t3E	*/	\n\t"\
		"addq	$0x80,%%rsi	/* r0E */	\n\t"\
		"subq	%%rax,%%rbx	/* p01 << 3 */	\n\t"\
		"subq	%%rax,%%rcx	/* p02 << 3 */	\n\t"\
		"subq	%%rax,%%rdx	/* p03 << 3 */	\n\t"\
		"movslq	%[__p04],%%rdi	\n\t"\
		"shlq	$3,%%rdi	\n\t"\
		"addq	%%rdi,%%rax	/* &a[j1+p1C] */	\n\t"\
		"addq	%%rax,%%rbx	\n\t"\
		"addq	%%rax,%%rcx	\n\t"\
		"addq	%%rax,%%rdx	\n\t"\
		"movaps	0x200(%%rsi),%%xmm4	\n\t"\
		"movaps	0x210(%%rsi),%%xmm5	\n\t"\
		"movaps	0x350(%%rsi),%%xmm3	/* cc1 */	\n\t"\
		"movaps	0x360(%%rsi),%%xmm2	\n\t"\
		"movaps	%%xmm4,%%xmm6	\n\t"\
		"movaps	%%xmm5,%%xmm7	\n\t"\
		"mulpd	%%xmm2,%%xmm4	\n\t"\
		"mulpd	%%xmm2,%%xmm5	\n\t"\
		"mulpd	%%xmm3,%%xmm6	\n\t"\
		"movaps	0x300(%%rsi),%%xmm0	\n\t"\
		"mulpd	%%xmm3,%%xmm7	\n\t"\
		"movaps	0x310(%%rsi),%%xmm1	\n\t"\
		"movaps	0x370(%%rsi),%%xmm3	/* cc3 */	\n\t"\
		"movaps	0x380(%%rsi),%%xmm2	\n\t"\
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
		"movaps	0x100(%%rsi),%%xmm1	\n\t"\
		"movaps	0x110(%%rsi),%%xmm3	\n\t"\
		"movaps	0x340(%%rsi),%%xmm0	/* ss0 */	\n\t"\
		"movaps	%%xmm1,%%xmm2	\n\t"\
		"mulpd	%%xmm0,%%xmm1	\n\t"\
		"mulpd	%%xmm3,%%xmm0	\n\t"\
		"mulpd	0x330(%%rsi),%%xmm2	/* cc0 */	\n\t"\
		"mulpd	0x330(%%rsi),%%xmm3	\n\t"\
		"addpd	%%xmm0,%%xmm2	\n\t"\
		"subpd	%%xmm1,%%xmm3	\n\t"\
		"movaps	     (%%rsi),%%xmm0	\n\t"\
		"movaps	0x010(%%rsi),%%xmm1	\n\t"\
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
		"movaps	%%xmm0,    (%%rbx)	\n\t"\
		"movaps	%%xmm2,    (%%rcx)	\n\t"\
		"movaps	%%xmm1,0x10(%%rbx)	\n\t"\
		"movaps	%%xmm3,0x10(%%rdx)	\n\t"\
		"addpd	%%xmm0,	%%xmm4	\n\t"\
		"addpd	%%xmm2,	%%xmm7	\n\t"\
		"addpd	%%xmm1,	%%xmm5	\n\t"\
		"addpd	%%xmm3,	%%xmm6	\n\t"\
		"movaps	%%xmm4,    (%%rax)	\n\t"\
		"movaps	%%xmm7,    (%%rdx)	\n\t"\
		"movaps	%%xmm5,0x10(%%rax)	\n\t"\
		"movaps	%%xmm6,0x10(%%rcx)	\n\t"\
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
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"	/* Clobbered registers */\
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
		"movq	%[__add0],%%rax\n\t"\
		"movslq	%[__p05],%%rbx\n\t"\
		"movslq	%[__p06],%%rcx\n\t"\
		"movslq	%[__p07],%%rdx\n\t"\
		"movslq	%[__p04],%%rdi	/* rdi will store copy of p4 throughout */\n\t"\
		"shlq	$3,%%rbx\n\t"\
		"shlq	$3,%%rcx\n\t"\
		"shlq	$3,%%rdx\n\t"\
		"shlq	$3,%%rdi\n\t"\
		"addq	%%rax,%%rbx\n\t"\
		"addq	%%rax,%%rcx\n\t"\
		"addq	%%rax,%%rdx\n\t"\
		"addq	%%rdi,%%rax\n\t"\
		"/* SSE2_RADIX8_DIT_0TWIDDLE(r00): */\n\t"\
		"movq	%[__r00]	,%%rsi 		\n\t"\
		"movaps	    (%%rax)	,%%xmm0		\n\t"\
		"movaps	0x10(%%rax)	,%%xmm1		\n\t"\
		"movaps	%%xmm0		,%%xmm2		\n\t"\
		"movaps	%%xmm1		,%%xmm3		\n\t"\
		"addpd	    (%%rbx)	,%%xmm2		\n\t"\
		"addpd	0x10(%%rbx)	,%%xmm3		\n\t"\
		"subpd	    (%%rbx)	,%%xmm0		\n\t"\
		"subpd	0x10(%%rbx)	,%%xmm1		\n\t"\
		"movaps	    (%%rcx)	,%%xmm4		\n\t"\
		"movaps	0x10(%%rcx)	,%%xmm5		\n\t"\
		"movaps	%%xmm4		,%%xmm6		\n\t"\
		"movaps	%%xmm5		,%%xmm7		\n\t"\
		"addpd	    (%%rdx)	,%%xmm6		\n\t"\
		"addpd	0x10(%%rdx)	,%%xmm7		\n\t"\
		"subpd	    (%%rdx)	,%%xmm4		\n\t"\
		"subpd	0x10(%%rdx)	,%%xmm5		\n\t"\
		"movaps	%%xmm6		,0xc0(%%rsi)\n\t"\
		"movaps	%%xmm7		,0xd0(%%rsi)\n\t"\
		"addpd	%%xmm2		,%%xmm6		\n\t"\
		"addpd	%%xmm3		,%%xmm7		\n\t"\
		"subpd	0xc0(%%rsi)	,%%xmm2		\n\t"\
		"subpd	0xd0(%%rsi)	,%%xmm3		\n\t"\
		"movaps	%%xmm2		,0xc0(%%rsi)\n\t"\
		"movaps	%%xmm3		,0xd0(%%rsi)\n\t"\
		"movaps	%%xmm4		,%%xmm2		\n\t"\
		"movaps	%%xmm5		,%%xmm3		\n\t"\
		"addpd	%%xmm0		,%%xmm5		\n\t"\
		"subpd	%%xmm3		,%%xmm0		\n\t"\
		"addpd	%%xmm1		,%%xmm4		\n\t"\
		"subpd	%%xmm2		,%%xmm1		\n\t"\
		"movq	%[__isrt2]	,%%rsi		\n\t"\
		"movaps	%%xmm5		,%%xmm2		\n\t"\
		"movaps	%%xmm1		,%%xmm3		\n\t"\
		"addpd	%%xmm1		,%%xmm5		\n\t"\
		"movaps	(%%rsi)		,%%xmm1		\n\t"\
		"subpd	%%xmm3		,%%xmm2		\n\t"\
		"movq	%[__r00]	,%%rsi		\n\t"\
		"mulpd	%%xmm1		,%%xmm5		\n\t"\
		"mulpd	%%xmm1		,%%xmm2		\n\t"\
		"movaps	%%xmm0		,%%xmm3		\n\t"\
		"movaps	%%xmm5		,0xa0(%%rsi)\n\t"\
		"movaps	%%xmm4		,%%xmm5		\n\t"\
		"addpd	%%xmm4		,%%xmm0		\n\t"\
		"movaps	%%xmm2		,0xb0(%%rsi)\n\t"\
		"subpd	%%xmm5		,%%xmm3		\n\t"\
		"mulpd	%%xmm1		,%%xmm0		\n\t"\
		"mulpd	%%xmm1		,%%xmm3		\n\t"\
		"movaps	%%xmm0		,0xe0(%%rsi)\n\t"\
		"movaps	%%xmm3		,0xf0(%%rsi)\n\t"\
		"subq	%%rdi		,%%rax		\n\t"\
		"subq	%%rdi		,%%rbx		\n\t"\
		"subq	%%rdi		,%%rcx		\n\t"\
		"subq	%%rdi		,%%rdx		\n\t"\
		"movaps	    (%%rax)	,%%xmm0		\n\t"\
		"movaps	0x10(%%rax)	,%%xmm1		\n\t"\
		"movaps	%%xmm0		,%%xmm2		\n\t"\
		"movaps	%%xmm1		,%%xmm3		\n\t"\
		"addpd	    (%%rbx)	,%%xmm2		\n\t"\
		"addpd	0x10(%%rbx)	,%%xmm3		\n\t"\
		"subpd	    (%%rbx)	,%%xmm0		\n\t"\
		"subpd	0x10(%%rbx)	,%%xmm1		\n\t"\
		"movaps	%%xmm6		,    (%%rsi)\n\t"\
		"movaps	%%xmm7		,0x10(%%rsi)\n\t"\
		"addpd	%%xmm6		,%%xmm6		\n\t"\
		"addpd	%%xmm7		,%%xmm7		\n\t"\
		"movaps	%%xmm6		,0x80(%%rsi)\n\t"\
		"movaps	%%xmm7		,0x90(%%rsi)\n\t"\
		"movaps	    (%%rcx)	,%%xmm4		\n\t"\
		"movaps	0x10(%%rcx)	,%%xmm5		\n\t"\
		"movaps	%%xmm4		,%%xmm6		\n\t"\
		"movaps	%%xmm5		,%%xmm7		\n\t"\
		"addpd	    (%%rdx)	,%%xmm6		\n\t"\
		"addpd	0x10(%%rdx)	,%%xmm7		\n\t"\
		"subpd	    (%%rdx)	,%%xmm4		\n\t"\
		"subpd	0x10(%%rdx)	,%%xmm5		\n\t"\
		"movaps	%%xmm6		,0x40(%%rsi)\n\t"\
		"movaps	%%xmm7		,0x50(%%rsi)\n\t"\
		"addpd	%%xmm2		,%%xmm6		\n\t"\
		"addpd	%%xmm3		,%%xmm7		\n\t"\
		"subpd	0x40(%%rsi)	,%%xmm2		\n\t"\
		"subpd	0x50(%%rsi)	,%%xmm3		\n\t"\
		"addpd	    (%%rsi)	,%%xmm6		\n\t"\
		"addpd	0x10(%%rsi)	,%%xmm7		\n\t"\
		"movaps	%%xmm6		,    (%%rsi)\n\t"\
		"movaps	%%xmm7		,0x10(%%rsi)\n\t"\
		"subpd	0x80(%%rsi)	,%%xmm6		\n\t"\
		"subpd	0x90(%%rsi)	,%%xmm7		\n\t"\
		"movaps	%%xmm6		,0x80(%%rsi)\n\t"\
		"movaps	%%xmm7		,0x90(%%rsi)\n\t"\
		"movaps	%%xmm4		,%%xmm6		\n\t"\
		"movaps	%%xmm5		,%%xmm7		\n\t"\
		"addpd	%%xmm0		,%%xmm5		\n\t"\
		"subpd	%%xmm7		,%%xmm0		\n\t"\
		"addpd	%%xmm1		,%%xmm4		\n\t"\
		"subpd	%%xmm6		,%%xmm1		\n\t"\
		"movaps	%%xmm5		,%%xmm6		\n\t"\
		"movaps	%%xmm1		,%%xmm7		\n\t"\
		"addpd	0xa0(%%rsi)	,%%xmm5		\n\t"\
		"subpd	0xb0(%%rsi)	,%%xmm1		\n\t"\
		"subpd	0xa0(%%rsi)	,%%xmm6		\n\t"\
		"addpd	0xb0(%%rsi)	,%%xmm7		\n\t"\
		"movaps	%%xmm5		,0x20(%%rsi)\n\t"\
		"movaps	%%xmm1		,0x30(%%rsi)\n\t"\
		"movaps	%%xmm6		,0xa0(%%rsi)\n\t"\
		"movaps	%%xmm7		,0xb0(%%rsi)\n\t"\
		"movaps	%%xmm2		,%%xmm6		\n\t"\
		"movaps	%%xmm3		,%%xmm7		\n\t"\
		"addpd	0xd0(%%rsi)	,%%xmm2		\n\t"\
		"subpd	0xc0(%%rsi)	,%%xmm3		\n\t"\
		"subpd	0xd0(%%rsi)	,%%xmm6		\n\t"\
		"addpd	0xc0(%%rsi)	,%%xmm7		\n\t"\
		"movaps	%%xmm2		,0x40(%%rsi)\n\t"\
		"movaps	%%xmm3		,0x50(%%rsi)\n\t"\
		"movaps	%%xmm6		,0xc0(%%rsi)\n\t"\
		"movaps	%%xmm7		,0xd0(%%rsi)\n\t"\
		"movaps	%%xmm0		,%%xmm6		\n\t"\
		"movaps	%%xmm4		,%%xmm7		\n\t"\
		"subpd	0xf0(%%rsi)	,%%xmm0		\n\t"\
		"subpd	0xe0(%%rsi)	,%%xmm4		\n\t"\
		"addpd	0xf0(%%rsi)	,%%xmm6		\n\t"\
		"addpd	0xe0(%%rsi)	,%%xmm7		\n\t"\
		"movaps	%%xmm0		,0x60(%%rsi)\n\t"\
		"movaps	%%xmm4		,0x70(%%rsi)\n\t"\
		"movaps	%%xmm6		,0xe0(%%rsi)\n\t"\
		"movaps	%%xmm7		,0xf0(%%rsi)\n\t"\
		"\n\t"\
		"movq	%[__add0],%%rax\n\t"\
		"movslq	%[__p08],%%rsi\n\t"\
		"shlq	$3,%%rsi\n\t"\
		"addq	%%rsi,%%rax	/* add1 = add0+p8 */\n\t"\
		"movslq	%[__p05],%%rbx\n\t"\
		"movslq	%[__p06],%%rcx\n\t"\
		"movslq	%[__p07],%%rdx\n\t"\
		"movslq	%[__p04],%%rdi\n\t"\
		"shlq	$3,%%rbx\n\t"\
		"shlq	$3,%%rcx\n\t"\
		"shlq	$3,%%rdx\n\t"\
		"shlq	$3,%%rdi\n\t"\
		"addq	%%rax,%%rbx\n\t"\
		"addq	%%rax,%%rcx\n\t"\
		"addq	%%rax,%%rdx\n\t"\
		"addq	%%rdi,%%rax\n\t"\
		"/* SSE2_RADIX8_DIT_0TWIDDLE(r10): */\n\t"\
		"movq	%[__r10]	,%%rsi 		\n\t"\
		"movaps	    (%%rax)	,%%xmm0		\n\t"\
		"movaps	0x10(%%rax)	,%%xmm1		\n\t"\
		"movaps	%%xmm0		,%%xmm2		\n\t"\
		"movaps	%%xmm1		,%%xmm3		\n\t"\
		"addpd	    (%%rbx)	,%%xmm2		\n\t"\
		"addpd	0x10(%%rbx)	,%%xmm3		\n\t"\
		"subpd	    (%%rbx)	,%%xmm0		\n\t"\
		"subpd	0x10(%%rbx)	,%%xmm1		\n\t"\
		"movaps	    (%%rcx)	,%%xmm4		\n\t"\
		"movaps	0x10(%%rcx)	,%%xmm5		\n\t"\
		"movaps	%%xmm4		,%%xmm6		\n\t"\
		"movaps	%%xmm5		,%%xmm7		\n\t"\
		"addpd	    (%%rdx)	,%%xmm6		\n\t"\
		"addpd	0x10(%%rdx)	,%%xmm7		\n\t"\
		"subpd	    (%%rdx)	,%%xmm4		\n\t"\
		"subpd	0x10(%%rdx)	,%%xmm5		\n\t"\
		"movaps	%%xmm6		,0xc0(%%rsi)\n\t"\
		"movaps	%%xmm7		,0xd0(%%rsi)\n\t"\
		"addpd	%%xmm2		,%%xmm6		\n\t"\
		"addpd	%%xmm3		,%%xmm7		\n\t"\
		"subpd	0xc0(%%rsi)	,%%xmm2		\n\t"\
		"subpd	0xd0(%%rsi)	,%%xmm3		\n\t"\
		"movaps	%%xmm2		,0xc0(%%rsi)\n\t"\
		"movaps	%%xmm3		,0xd0(%%rsi)\n\t"\
		"movaps	%%xmm4		,%%xmm2		\n\t"\
		"movaps	%%xmm5		,%%xmm3		\n\t"\
		"addpd	%%xmm0		,%%xmm5		\n\t"\
		"subpd	%%xmm3		,%%xmm0		\n\t"\
		"addpd	%%xmm1		,%%xmm4		\n\t"\
		"subpd	%%xmm2		,%%xmm1		\n\t"\
		"movq	%[__isrt2]	,%%rsi		\n\t"\
		"movaps	%%xmm5		,%%xmm2		\n\t"\
		"movaps	%%xmm1		,%%xmm3		\n\t"\
		"addpd	%%xmm1		,%%xmm5		\n\t"\
		"movaps	(%%rsi)		,%%xmm1		\n\t"\
		"subpd	%%xmm3		,%%xmm2		\n\t"\
		"movq	%[__r10]	,%%rsi		\n\t"\
		"mulpd	%%xmm1		,%%xmm5		\n\t"\
		"mulpd	%%xmm1		,%%xmm2		\n\t"\
		"movaps	%%xmm0		,%%xmm3		\n\t"\
		"movaps	%%xmm5		,0xa0(%%rsi)\n\t"\
		"movaps	%%xmm4		,%%xmm5		\n\t"\
		"addpd	%%xmm4		,%%xmm0		\n\t"\
		"movaps	%%xmm2		,0xb0(%%rsi)\n\t"\
		"subpd	%%xmm5		,%%xmm3		\n\t"\
		"mulpd	%%xmm1		,%%xmm0		\n\t"\
		"mulpd	%%xmm1		,%%xmm3		\n\t"\
		"movaps	%%xmm0		,0xe0(%%rsi)\n\t"\
		"movaps	%%xmm3		,0xf0(%%rsi)\n\t"\
		"subq	%%rdi		,%%rax		\n\t"\
		"subq	%%rdi		,%%rbx		\n\t"\
		"subq	%%rdi		,%%rcx		\n\t"\
		"subq	%%rdi		,%%rdx		\n\t"\
		"movaps	    (%%rax)	,%%xmm0		\n\t"\
		"movaps	0x10(%%rax)	,%%xmm1		\n\t"\
		"movaps	%%xmm0		,%%xmm2		\n\t"\
		"movaps	%%xmm1		,%%xmm3		\n\t"\
		"addpd	    (%%rbx)	,%%xmm2		\n\t"\
		"addpd	0x10(%%rbx)	,%%xmm3		\n\t"\
		"subpd	    (%%rbx)	,%%xmm0		\n\t"\
		"subpd	0x10(%%rbx)	,%%xmm1		\n\t"\
		"movaps	%%xmm6		,    (%%rsi)\n\t"\
		"movaps	%%xmm7		,0x10(%%rsi)\n\t"\
		"addpd	%%xmm6		,%%xmm6		\n\t"\
		"addpd	%%xmm7		,%%xmm7		\n\t"\
		"movaps	%%xmm6		,0x80(%%rsi)\n\t"\
		"movaps	%%xmm7		,0x90(%%rsi)\n\t"\
		"movaps	    (%%rcx)	,%%xmm4		\n\t"\
		"movaps	0x10(%%rcx)	,%%xmm5		\n\t"\
		"movaps	%%xmm4		,%%xmm6		\n\t"\
		"movaps	%%xmm5		,%%xmm7		\n\t"\
		"addpd	    (%%rdx)	,%%xmm6		\n\t"\
		"addpd	0x10(%%rdx)	,%%xmm7		\n\t"\
		"subpd	    (%%rdx)	,%%xmm4		\n\t"\
		"subpd	0x10(%%rdx)	,%%xmm5		\n\t"\
		"movaps	%%xmm6		,0x40(%%rsi)\n\t"\
		"movaps	%%xmm7		,0x50(%%rsi)\n\t"\
		"addpd	%%xmm2		,%%xmm6		\n\t"\
		"addpd	%%xmm3		,%%xmm7		\n\t"\
		"subpd	0x40(%%rsi)	,%%xmm2		\n\t"\
		"subpd	0x50(%%rsi)	,%%xmm3		\n\t"\
		"addpd	    (%%rsi)	,%%xmm6		\n\t"\
		"addpd	0x10(%%rsi)	,%%xmm7		\n\t"\
		"movaps	%%xmm6		,    (%%rsi)\n\t"\
		"movaps	%%xmm7		,0x10(%%rsi)\n\t"\
		"subpd	0x80(%%rsi)	,%%xmm6		\n\t"\
		"subpd	0x90(%%rsi)	,%%xmm7		\n\t"\
		"movaps	%%xmm6		,0x80(%%rsi)\n\t"\
		"movaps	%%xmm7		,0x90(%%rsi)\n\t"\
		"movaps	%%xmm4		,%%xmm6		\n\t"\
		"movaps	%%xmm5		,%%xmm7		\n\t"\
		"addpd	%%xmm0		,%%xmm5		\n\t"\
		"subpd	%%xmm7		,%%xmm0		\n\t"\
		"addpd	%%xmm1		,%%xmm4		\n\t"\
		"subpd	%%xmm6		,%%xmm1		\n\t"\
		"movaps	%%xmm5		,%%xmm6		\n\t"\
		"movaps	%%xmm1		,%%xmm7		\n\t"\
		"addpd	0xa0(%%rsi)	,%%xmm5		\n\t"\
		"subpd	0xb0(%%rsi)	,%%xmm1		\n\t"\
		"subpd	0xa0(%%rsi)	,%%xmm6		\n\t"\
		"addpd	0xb0(%%rsi)	,%%xmm7		\n\t"\
		"movaps	%%xmm5		,0x20(%%rsi)\n\t"\
		"movaps	%%xmm1		,0x30(%%rsi)\n\t"\
		"movaps	%%xmm6		,0xa0(%%rsi)\n\t"\
		"movaps	%%xmm7		,0xb0(%%rsi)\n\t"\
		"movaps	%%xmm2		,%%xmm6		\n\t"\
		"movaps	%%xmm3		,%%xmm7		\n\t"\
		"addpd	0xd0(%%rsi)	,%%xmm2		\n\t"\
		"subpd	0xc0(%%rsi)	,%%xmm3		\n\t"\
		"subpd	0xd0(%%rsi)	,%%xmm6		\n\t"\
		"addpd	0xc0(%%rsi)	,%%xmm7		\n\t"\
		"movaps	%%xmm2		,0x40(%%rsi)\n\t"\
		"movaps	%%xmm3		,0x50(%%rsi)\n\t"\
		"movaps	%%xmm6		,0xc0(%%rsi)\n\t"\
		"movaps	%%xmm7		,0xd0(%%rsi)\n\t"\
		"movaps	%%xmm0		,%%xmm6		\n\t"\
		"movaps	%%xmm4		,%%xmm7		\n\t"\
		"subpd	0xf0(%%rsi)	,%%xmm0		\n\t"\
		"subpd	0xe0(%%rsi)	,%%xmm4		\n\t"\
		"addpd	0xf0(%%rsi)	,%%xmm6		\n\t"\
		"addpd	0xe0(%%rsi)	,%%xmm7		\n\t"\
		"movaps	%%xmm0		,0x60(%%rsi)\n\t"\
		"movaps	%%xmm4		,0x70(%%rsi)\n\t"\
		"movaps	%%xmm6		,0xe0(%%rsi)\n\t"\
		"movaps	%%xmm7		,0xf0(%%rsi)\n\t"\
		"\n\t"\
		"movq	%[__add0],%%rax\n\t"\
		"movslq	%[__p10],%%rsi\n\t"\
		"shlq	$3,%%rsi\n\t"\
		"addq	%%rsi,%%rax	/* add2 = add0+p10*/\n\t"\
		"movslq	%[__p05],%%rbx\n\t"\
		"movslq	%[__p06],%%rcx\n\t"\
		"movslq	%[__p07],%%rdx\n\t"\
		"movslq	%[__p04],%%rdi\n\t"\
		"shlq	$3,%%rbx\n\t"\
		"shlq	$3,%%rcx\n\t"\
		"shlq	$3,%%rdx\n\t"\
		"shlq	$3,%%rdi\n\t"\
		"addq	%%rax,%%rbx\n\t"\
		"addq	%%rax,%%rcx\n\t"\
		"addq	%%rax,%%rdx\n\t"\
		"addq	%%rdi,%%rax\n\t"\
		"/* SSE2_RADIX8_DIT_0TWIDDLE(r20): */\n\t"\
		"movq	%[__r20]	,%%rsi 		\n\t"\
		"movaps	    (%%rax)	,%%xmm0		\n\t"\
		"movaps	0x10(%%rax)	,%%xmm1		\n\t"\
		"movaps	%%xmm0		,%%xmm2		\n\t"\
		"movaps	%%xmm1		,%%xmm3		\n\t"\
		"addpd	    (%%rbx)	,%%xmm2		\n\t"\
		"addpd	0x10(%%rbx)	,%%xmm3		\n\t"\
		"subpd	    (%%rbx)	,%%xmm0		\n\t"\
		"subpd	0x10(%%rbx)	,%%xmm1		\n\t"\
		"movaps	    (%%rcx)	,%%xmm4		\n\t"\
		"movaps	0x10(%%rcx)	,%%xmm5		\n\t"\
		"movaps	%%xmm4		,%%xmm6		\n\t"\
		"movaps	%%xmm5		,%%xmm7		\n\t"\
		"addpd	    (%%rdx)	,%%xmm6		\n\t"\
		"addpd	0x10(%%rdx)	,%%xmm7		\n\t"\
		"subpd	    (%%rdx)	,%%xmm4		\n\t"\
		"subpd	0x10(%%rdx)	,%%xmm5		\n\t"\
		"movaps	%%xmm6		,0xc0(%%rsi)\n\t"\
		"movaps	%%xmm7		,0xd0(%%rsi)\n\t"\
		"addpd	%%xmm2		,%%xmm6		\n\t"\
		"addpd	%%xmm3		,%%xmm7		\n\t"\
		"subpd	0xc0(%%rsi)	,%%xmm2		\n\t"\
		"subpd	0xd0(%%rsi)	,%%xmm3		\n\t"\
		"movaps	%%xmm2		,0xc0(%%rsi)\n\t"\
		"movaps	%%xmm3		,0xd0(%%rsi)\n\t"\
		"movaps	%%xmm4		,%%xmm2		\n\t"\
		"movaps	%%xmm5		,%%xmm3		\n\t"\
		"addpd	%%xmm0		,%%xmm5		\n\t"\
		"subpd	%%xmm3		,%%xmm0		\n\t"\
		"addpd	%%xmm1		,%%xmm4		\n\t"\
		"subpd	%%xmm2		,%%xmm1		\n\t"\
		"movq	%[__isrt2]	,%%rsi		\n\t"\
		"movaps	%%xmm5		,%%xmm2		\n\t"\
		"movaps	%%xmm1		,%%xmm3		\n\t"\
		"addpd	%%xmm1		,%%xmm5		\n\t"\
		"movaps	(%%rsi)		,%%xmm1		\n\t"\
		"subpd	%%xmm3		,%%xmm2		\n\t"\
		"movq	%[__r20]	,%%rsi		\n\t"\
		"mulpd	%%xmm1		,%%xmm5		\n\t"\
		"mulpd	%%xmm1		,%%xmm2		\n\t"\
		"movaps	%%xmm0		,%%xmm3		\n\t"\
		"movaps	%%xmm5		,0xa0(%%rsi)\n\t"\
		"movaps	%%xmm4		,%%xmm5		\n\t"\
		"addpd	%%xmm4		,%%xmm0		\n\t"\
		"movaps	%%xmm2		,0xb0(%%rsi)\n\t"\
		"subpd	%%xmm5		,%%xmm3		\n\t"\
		"mulpd	%%xmm1		,%%xmm0		\n\t"\
		"mulpd	%%xmm1		,%%xmm3		\n\t"\
		"movaps	%%xmm0		,0xe0(%%rsi)\n\t"\
		"movaps	%%xmm3		,0xf0(%%rsi)\n\t"\
		"subq	%%rdi		,%%rax		\n\t"\
		"subq	%%rdi		,%%rbx		\n\t"\
		"subq	%%rdi		,%%rcx		\n\t"\
		"subq	%%rdi		,%%rdx		\n\t"\
		"movaps	    (%%rax)	,%%xmm0		\n\t"\
		"movaps	0x10(%%rax)	,%%xmm1		\n\t"\
		"movaps	%%xmm0		,%%xmm2		\n\t"\
		"movaps	%%xmm1		,%%xmm3		\n\t"\
		"addpd	    (%%rbx)	,%%xmm2		\n\t"\
		"addpd	0x10(%%rbx)	,%%xmm3		\n\t"\
		"subpd	    (%%rbx)	,%%xmm0		\n\t"\
		"subpd	0x10(%%rbx)	,%%xmm1		\n\t"\
		"movaps	%%xmm6		,    (%%rsi)\n\t"\
		"movaps	%%xmm7		,0x10(%%rsi)\n\t"\
		"addpd	%%xmm6		,%%xmm6		\n\t"\
		"addpd	%%xmm7		,%%xmm7		\n\t"\
		"movaps	%%xmm6		,0x80(%%rsi)\n\t"\
		"movaps	%%xmm7		,0x90(%%rsi)\n\t"\
		"movaps	    (%%rcx)	,%%xmm4		\n\t"\
		"movaps	0x10(%%rcx)	,%%xmm5		\n\t"\
		"movaps	%%xmm4		,%%xmm6		\n\t"\
		"movaps	%%xmm5		,%%xmm7		\n\t"\
		"addpd	    (%%rdx)	,%%xmm6		\n\t"\
		"addpd	0x10(%%rdx)	,%%xmm7		\n\t"\
		"subpd	    (%%rdx)	,%%xmm4		\n\t"\
		"subpd	0x10(%%rdx)	,%%xmm5		\n\t"\
		"movaps	%%xmm6		,0x40(%%rsi)\n\t"\
		"movaps	%%xmm7		,0x50(%%rsi)\n\t"\
		"addpd	%%xmm2		,%%xmm6		\n\t"\
		"addpd	%%xmm3		,%%xmm7		\n\t"\
		"subpd	0x40(%%rsi)	,%%xmm2		\n\t"\
		"subpd	0x50(%%rsi)	,%%xmm3		\n\t"\
		"addpd	    (%%rsi)	,%%xmm6		\n\t"\
		"addpd	0x10(%%rsi)	,%%xmm7		\n\t"\
		"movaps	%%xmm6		,    (%%rsi)\n\t"\
		"movaps	%%xmm7		,0x10(%%rsi)\n\t"\
		"subpd	0x80(%%rsi)	,%%xmm6		\n\t"\
		"subpd	0x90(%%rsi)	,%%xmm7		\n\t"\
		"movaps	%%xmm6		,0x80(%%rsi)\n\t"\
		"movaps	%%xmm7		,0x90(%%rsi)\n\t"\
		"movaps	%%xmm4		,%%xmm6		\n\t"\
		"movaps	%%xmm5		,%%xmm7		\n\t"\
		"addpd	%%xmm0		,%%xmm5		\n\t"\
		"subpd	%%xmm7		,%%xmm0		\n\t"\
		"addpd	%%xmm1		,%%xmm4		\n\t"\
		"subpd	%%xmm6		,%%xmm1		\n\t"\
		"movaps	%%xmm5		,%%xmm6		\n\t"\
		"movaps	%%xmm1		,%%xmm7		\n\t"\
		"addpd	0xa0(%%rsi)	,%%xmm5		\n\t"\
		"subpd	0xb0(%%rsi)	,%%xmm1		\n\t"\
		"subpd	0xa0(%%rsi)	,%%xmm6		\n\t"\
		"addpd	0xb0(%%rsi)	,%%xmm7		\n\t"\
		"movaps	%%xmm5		,0x20(%%rsi)\n\t"\
		"movaps	%%xmm1		,0x30(%%rsi)\n\t"\
		"movaps	%%xmm6		,0xa0(%%rsi)\n\t"\
		"movaps	%%xmm7		,0xb0(%%rsi)\n\t"\
		"movaps	%%xmm2		,%%xmm6		\n\t"\
		"movaps	%%xmm3		,%%xmm7		\n\t"\
		"addpd	0xd0(%%rsi)	,%%xmm2		\n\t"\
		"subpd	0xc0(%%rsi)	,%%xmm3		\n\t"\
		"subpd	0xd0(%%rsi)	,%%xmm6		\n\t"\
		"addpd	0xc0(%%rsi)	,%%xmm7		\n\t"\
		"movaps	%%xmm2		,0x40(%%rsi)\n\t"\
		"movaps	%%xmm3		,0x50(%%rsi)\n\t"\
		"movaps	%%xmm6		,0xc0(%%rsi)\n\t"\
		"movaps	%%xmm7		,0xd0(%%rsi)\n\t"\
		"movaps	%%xmm0		,%%xmm6		\n\t"\
		"movaps	%%xmm4		,%%xmm7		\n\t"\
		"subpd	0xf0(%%rsi)	,%%xmm0		\n\t"\
		"subpd	0xe0(%%rsi)	,%%xmm4		\n\t"\
		"addpd	0xf0(%%rsi)	,%%xmm6		\n\t"\
		"addpd	0xe0(%%rsi)	,%%xmm7		\n\t"\
		"movaps	%%xmm0		,0x60(%%rsi)\n\t"\
		"movaps	%%xmm4		,0x70(%%rsi)\n\t"\
		"movaps	%%xmm6		,0xe0(%%rsi)\n\t"\
		"movaps	%%xmm7		,0xf0(%%rsi)\n\t"\
		"\n\t"\
		"movq	%[__add0],%%rax\n\t"\
		"movslq	%[__p18],%%rsi\n\t"\
		"shlq	$3,%%rsi\n\t"\
		"addq	%%rsi,%%rax	/* add3 = add0+p18*/\n\t"\
		"movslq	%[__p05],%%rbx\n\t"\
		"movslq	%[__p06],%%rcx\n\t"\
		"movslq	%[__p07],%%rdx\n\t"\
		"movslq	%[__p04],%%rdi\n\t"\
		"shlq	$3,%%rbx\n\t"\
		"shlq	$3,%%rcx\n\t"\
		"shlq	$3,%%rdx\n\t"\
		"shlq	$3,%%rdi\n\t"\
		"addq	%%rax,%%rbx\n\t"\
		"addq	%%rax,%%rcx\n\t"\
		"addq	%%rax,%%rdx\n\t"\
		"addq	%%rdi,%%rax\n\t"\
		"/* SSE2_RADIX8_DIT_0TWIDDLE(r30): */\n\t"\
		"movq	%[__r30]	,%%rsi 		\n\t"\
		"movaps	    (%%rax)	,%%xmm0		\n\t"\
		"movaps	0x10(%%rax)	,%%xmm1		\n\t"\
		"movaps	%%xmm0		,%%xmm2		\n\t"\
		"movaps	%%xmm1		,%%xmm3		\n\t"\
		"addpd	    (%%rbx)	,%%xmm2		\n\t"\
		"addpd	0x10(%%rbx)	,%%xmm3		\n\t"\
		"subpd	    (%%rbx)	,%%xmm0		\n\t"\
		"subpd	0x10(%%rbx)	,%%xmm1		\n\t"\
		"movaps	    (%%rcx)	,%%xmm4		\n\t"\
		"movaps	0x10(%%rcx)	,%%xmm5		\n\t"\
		"movaps	%%xmm4		,%%xmm6		\n\t"\
		"movaps	%%xmm5		,%%xmm7		\n\t"\
		"addpd	    (%%rdx)	,%%xmm6		\n\t"\
		"addpd	0x10(%%rdx)	,%%xmm7		\n\t"\
		"subpd	    (%%rdx)	,%%xmm4		\n\t"\
		"subpd	0x10(%%rdx)	,%%xmm5		\n\t"\
		"movaps	%%xmm6		,0xc0(%%rsi)\n\t"\
		"movaps	%%xmm7		,0xd0(%%rsi)\n\t"\
		"addpd	%%xmm2		,%%xmm6		\n\t"\
		"addpd	%%xmm3		,%%xmm7		\n\t"\
		"subpd	0xc0(%%rsi)	,%%xmm2		\n\t"\
		"subpd	0xd0(%%rsi)	,%%xmm3		\n\t"\
		"movaps	%%xmm2		,0xc0(%%rsi)\n\t"\
		"movaps	%%xmm3		,0xd0(%%rsi)\n\t"\
		"movaps	%%xmm4		,%%xmm2		\n\t"\
		"movaps	%%xmm5		,%%xmm3		\n\t"\
		"addpd	%%xmm0		,%%xmm5		\n\t"\
		"subpd	%%xmm3		,%%xmm0		\n\t"\
		"addpd	%%xmm1		,%%xmm4		\n\t"\
		"subpd	%%xmm2		,%%xmm1		\n\t"\
		"movq	%[__isrt2]	,%%rsi		\n\t"\
		"movaps	%%xmm5		,%%xmm2		\n\t"\
		"movaps	%%xmm1		,%%xmm3		\n\t"\
		"addpd	%%xmm1		,%%xmm5		\n\t"\
		"movaps	(%%rsi)		,%%xmm1		\n\t"\
		"subpd	%%xmm3		,%%xmm2		\n\t"\
		"movq	%[__r30]	,%%rsi		\n\t"\
		"mulpd	%%xmm1		,%%xmm5		\n\t"\
		"mulpd	%%xmm1		,%%xmm2		\n\t"\
		"movaps	%%xmm0		,%%xmm3		\n\t"\
		"movaps	%%xmm5		,0xa0(%%rsi)\n\t"\
		"movaps	%%xmm4		,%%xmm5		\n\t"\
		"addpd	%%xmm4		,%%xmm0		\n\t"\
		"movaps	%%xmm2		,0xb0(%%rsi)\n\t"\
		"subpd	%%xmm5		,%%xmm3		\n\t"\
		"mulpd	%%xmm1		,%%xmm0		\n\t"\
		"mulpd	%%xmm1		,%%xmm3		\n\t"\
		"movaps	%%xmm0		,0xe0(%%rsi)\n\t"\
		"movaps	%%xmm3		,0xf0(%%rsi)\n\t"\
		"subq	%%rdi		,%%rax		\n\t"\
		"subq	%%rdi		,%%rbx		\n\t"\
		"subq	%%rdi		,%%rcx		\n\t"\
		"subq	%%rdi		,%%rdx		\n\t"\
		"movaps	    (%%rax)	,%%xmm0		\n\t"\
		"movaps	0x10(%%rax)	,%%xmm1		\n\t"\
		"movaps	%%xmm0		,%%xmm2		\n\t"\
		"movaps	%%xmm1		,%%xmm3		\n\t"\
		"addpd	    (%%rbx)	,%%xmm2		\n\t"\
		"addpd	0x10(%%rbx)	,%%xmm3		\n\t"\
		"subpd	    (%%rbx)	,%%xmm0		\n\t"\
		"subpd	0x10(%%rbx)	,%%xmm1		\n\t"\
		"movaps	%%xmm6		,    (%%rsi)\n\t"\
		"movaps	%%xmm7		,0x10(%%rsi)\n\t"\
		"addpd	%%xmm6		,%%xmm6		\n\t"\
		"addpd	%%xmm7		,%%xmm7		\n\t"\
		"movaps	%%xmm6		,0x80(%%rsi)\n\t"\
		"movaps	%%xmm7		,0x90(%%rsi)\n\t"\
		"movaps	    (%%rcx)	,%%xmm4		\n\t"\
		"movaps	0x10(%%rcx)	,%%xmm5		\n\t"\
		"movaps	%%xmm4		,%%xmm6		\n\t"\
		"movaps	%%xmm5		,%%xmm7		\n\t"\
		"addpd	    (%%rdx)	,%%xmm6		\n\t"\
		"addpd	0x10(%%rdx)	,%%xmm7		\n\t"\
		"subpd	    (%%rdx)	,%%xmm4		\n\t"\
		"subpd	0x10(%%rdx)	,%%xmm5		\n\t"\
		"movaps	%%xmm6		,0x40(%%rsi)\n\t"\
		"movaps	%%xmm7		,0x50(%%rsi)\n\t"\
		"addpd	%%xmm2		,%%xmm6		\n\t"\
		"addpd	%%xmm3		,%%xmm7		\n\t"\
		"subpd	0x40(%%rsi)	,%%xmm2		\n\t"\
		"subpd	0x50(%%rsi)	,%%xmm3		\n\t"\
		"addpd	    (%%rsi)	,%%xmm6		\n\t"\
		"addpd	0x10(%%rsi)	,%%xmm7		\n\t"\
		"movaps	%%xmm6		,    (%%rsi)\n\t"\
		"movaps	%%xmm7		,0x10(%%rsi)\n\t"\
		"subpd	0x80(%%rsi)	,%%xmm6		\n\t"\
		"subpd	0x90(%%rsi)	,%%xmm7		\n\t"\
		"movaps	%%xmm6		,0x80(%%rsi)\n\t"\
		"movaps	%%xmm7		,0x90(%%rsi)\n\t"\
		"movaps	%%xmm4		,%%xmm6		\n\t"\
		"movaps	%%xmm5		,%%xmm7		\n\t"\
		"addpd	%%xmm0		,%%xmm5		\n\t"\
		"subpd	%%xmm7		,%%xmm0		\n\t"\
		"addpd	%%xmm1		,%%xmm4		\n\t"\
		"subpd	%%xmm6		,%%xmm1		\n\t"\
		"movaps	%%xmm5		,%%xmm6		\n\t"\
		"movaps	%%xmm1		,%%xmm7		\n\t"\
		"addpd	0xa0(%%rsi)	,%%xmm5		\n\t"\
		"subpd	0xb0(%%rsi)	,%%xmm1		\n\t"\
		"subpd	0xa0(%%rsi)	,%%xmm6		\n\t"\
		"addpd	0xb0(%%rsi)	,%%xmm7		\n\t"\
		"movaps	%%xmm5		,0x20(%%rsi)\n\t"\
		"movaps	%%xmm1		,0x30(%%rsi)\n\t"\
		"movaps	%%xmm6		,0xa0(%%rsi)\n\t"\
		"movaps	%%xmm7		,0xb0(%%rsi)\n\t"\
		"movaps	%%xmm2		,%%xmm6		\n\t"\
		"movaps	%%xmm3		,%%xmm7		\n\t"\
		"addpd	0xd0(%%rsi)	,%%xmm2		\n\t"\
		"subpd	0xc0(%%rsi)	,%%xmm3		\n\t"\
		"subpd	0xd0(%%rsi)	,%%xmm6		\n\t"\
		"addpd	0xc0(%%rsi)	,%%xmm7		\n\t"\
		"movaps	%%xmm2		,0x40(%%rsi)\n\t"\
		"movaps	%%xmm3		,0x50(%%rsi)\n\t"\
		"movaps	%%xmm6		,0xc0(%%rsi)\n\t"\
		"movaps	%%xmm7		,0xd0(%%rsi)\n\t"\
		"movaps	%%xmm0		,%%xmm6		\n\t"\
		"movaps	%%xmm4		,%%xmm7		\n\t"\
		"subpd	0xf0(%%rsi)	,%%xmm0		\n\t"\
		"subpd	0xe0(%%rsi)	,%%xmm4		\n\t"\
		"addpd	0xf0(%%rsi)	,%%xmm6		\n\t"\
		"addpd	0xe0(%%rsi)	,%%xmm7		\n\t"\
		"movaps	%%xmm0		,0x60(%%rsi)\n\t"\
		"movaps	%%xmm4		,0x70(%%rsi)\n\t"\
		"movaps	%%xmm6		,0xe0(%%rsi)\n\t"\
		"movaps	%%xmm7		,0xf0(%%rsi)\n\t"\
		"\n\t"\
		"/*...and now do eight radix-4 transforms, including the internal twiddle factors: */\n\t"\
		"\n\t"\
		"/*...Block 1: t00,t10,t20,t30	*/\n\t"\
		"movq	%[__add0],%%rax\n\t"\
		"movslq	%[__p08],%%rbx\n\t"\
		"movq	%[__r00],%%rcx\n\t"\
		"movq	%[__r10],%%rdx\n\t"\
		"movslq	%[__p10],%%rdi	/* rdi will store copy of p10 throughout */\n\t"\
		"shlq	$3,%%rbx\n\t"\
		"shlq	$3,%%rdi\n\t"\
		"addq	%%rax	,%%rbx	/* add1 = add0+p8 */\n\t"\
		"\n\t"\
		"movaps	     (%%rdx),%%xmm2	/* t10 */\n\t"\
		"movaps	0x200(%%rdx),%%xmm4	/* t30 */\n\t"\
		"movaps	0x010(%%rdx),%%xmm3	/* t11 */\n\t"\
		"movaps	0x210(%%rdx),%%xmm5	/* t31 */\n\t"\
		"movaps	     (%%rcx),%%xmm0	/* t00 */\n\t"\
		"movaps	0x200(%%rcx),%%xmm6	/* t20 */\n\t"\
		"movaps	0x010(%%rcx),%%xmm1	/* t01 */\n\t"\
		"movaps	0x210(%%rcx),%%xmm7	/* t21 */\n\t"\
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
		"movq	%[__isrt2],%%rsi\n\t"\
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
		"addq	$0x070,%%rsi	/* c00 */\n\t"\
		"movaps	%%xmm2,     (%%rcx)\n\t"\
		"movaps	%%xmm0,0x010(%%rdx)\n\t"\
		"movaps	%%xmm3,0x010(%%rcx)\n\t"\
		"movaps	%%xmm6,     (%%rdx)\n\t"\
		"movaps	%%xmm4,%%xmm2\n\t"\
		"movaps	%%xmm7,%%xmm0\n\t"\
		"movaps	%%xmm5,%%xmm3\n\t"\
		"movaps	%%xmm1,%%xmm6\n\t"\
		"mulpd	     (%%rsi),%%xmm4\n\t"\
		"mulpd	0x020(%%rsi),%%xmm7\n\t"\
		"mulpd	     (%%rsi),%%xmm5\n\t"\
		"mulpd	0x020(%%rsi),%%xmm1\n\t"\
		"mulpd	0x010(%%rsi),%%xmm2\n\t"\
		"mulpd	0x030(%%rsi),%%xmm0\n\t"\
		"mulpd	0x010(%%rsi),%%xmm3\n\t"\
		"mulpd	0x030(%%rsi),%%xmm6\n\t"\
		"subpd	%%xmm2,%%xmm5\n\t"\
		"subpd	%%xmm0,%%xmm1\n\t"\
		"addpd	%%xmm3,%%xmm4\n\t"\
		"addpd	%%xmm6,%%xmm7\n\t"\
		"movaps	%%xmm5,0x10(%%rax)\n\t"\
		"movaps	%%xmm1,0x10(%%rbx)\n\t"\
		"movaps	%%xmm4,    (%%rax)\n\t"\
		"movaps	%%xmm7,    (%%rbx)\n\t"\
		"addq	%%rdi,%%rax\n\t"\
		"addq	%%rdi,%%rbx\n\t"\
		"addq	$0x40,%%rsi\n\t"\
		"movaps	     (%%rcx),%%xmm4\n\t"\
		"movaps	0x010(%%rdx),%%xmm0\n\t"\
		"movaps	0x010(%%rcx),%%xmm5\n\t"\
		"movaps	     (%%rdx),%%xmm6\n\t"\
		"movaps	%%xmm4,%%xmm2\n\t"\
		"movaps	%%xmm0,%%xmm1\n\t"\
		"movaps	%%xmm5,%%xmm3\n\t"\
		"movaps	%%xmm6,%%xmm7\n\t"\
		"mulpd	     (%%rsi),%%xmm4\n\t"\
		"mulpd	0x020(%%rsi),%%xmm0\n\t"\
		"mulpd	     (%%rsi),%%xmm5\n\t"\
		"mulpd	0x020(%%rsi),%%xmm6\n\t"\
		"mulpd	0x010(%%rsi),%%xmm2\n\t"\
		"mulpd	0x030(%%rsi),%%xmm1\n\t"\
		"mulpd	0x010(%%rsi),%%xmm3\n\t"\
		"mulpd	0x030(%%rsi),%%xmm7\n\t"\
		"subpd	%%xmm2,%%xmm5\n\t"\
		"subpd	%%xmm1,%%xmm6\n\t"\
		"addpd	%%xmm3,%%xmm4\n\t"\
		"addpd	%%xmm7,%%xmm0\n\t"\
		"movaps	%%xmm5,0x10(%%rax)\n\t"\
		"movaps	%%xmm6,0x10(%%rbx)\n\t"\
		"movaps	%%xmm4,    (%%rax)\n\t"\
		"movaps	%%xmm0,    (%%rbx)\n\t"\
		"\n\t"\
		"/*...Block 2: t02,t12,t22,t32	*/\n\t"\
		"movq	%[__add0],%%rax\n\t"\
		"movslq	%[__p01],%%rsi\n\t"\
		"movslq	%[__p08],%%rbx\n\t"\
		"shlq	$3,%%rsi\n\t"\
		"shlq	$3,%%rbx\n\t"\
		"addq	%%rsi	,%%rax	/* add0 = &a[j1+p1] */\n\t"\
		"addq	%%rax	,%%rbx	/* add1 = add0+p8 */\n\t"\
		"movq	%[__r02],%%rcx\n\t"\
		"movq	%[__r12],%%rdx\n\t"\
		"movq	%[__isrt2],%%rsi\n\t"\
		"addq	$0x200,%%rcx\n\t"\
		"addq	$0x200,%%rdx\n\t"\
		"\n\t"\
		"addq	$0x030,%%rsi	/* cc1 */\n\t"\
		"movaps	     (%%rcx),%%xmm4	/* t22 */\n\t"\
		"movaps	0x010(%%rcx),%%xmm5	/* t23 */\n\t"\
		"movaps	     (%%rsi),%%xmm2	/* c32_1 */\n\t"\
		"movaps	0x010(%%rsi),%%xmm3	/* s32_1 */\n\t"\
		"movaps	%%xmm4,%%xmm6	/* xmm2 <- cpy t22 */\n\t"\
		"movaps	%%xmm5,%%xmm7	/* xmm3 <- cpy t23 */\n\t"\
		"\n\t"\
		"mulpd	%%xmm2,%%xmm4	/* t22*c32_1 */\n\t"\
		"mulpd	%%xmm2,%%xmm5	/* t23*c32_1 */\n\t"\
		"mulpd	%%xmm3,%%xmm6	/* t22*s32_1 */\n\t"\
		"movaps	     (%%rdx),%%xmm0	/* t32 */\n\t"\
		"mulpd	%%xmm3,%%xmm7	/* t23*s32_1 */\n\t"\
		"movaps	0x010(%%rdx),%%xmm1	/* t33 */\n\t"\
		"addq	$0x020,%%rsi	/* cc3 */\n\t"\
		"movaps	     (%%rsi),%%xmm2	/* c32_3 */\n\t"\
		"movaps	0x010(%%rsi),%%xmm3	/* s32_3 */\n\t"\
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
		"subq	$0x200,%%rcx\n\t"\
		"subq	$0x200,%%rdx\n\t"\
		"subq	$0x040,%%rsi	/* cc0 */\n\t"\
		"movaps	     (%%rdx),%%xmm1	/* t12 */\n\t"\
		"movaps	0x010(%%rdx),%%xmm3	/* t13 */\n\t"\
		"movaps	0x010(%%rsi),%%xmm2	/* s */\n\t"\
		"movaps	%%xmm1,%%xmm0	/* cpy t12 */\n\t"\
		"mulpd	%%xmm2,%%xmm1	/* t12*s */\n\t"\
		"mulpd	%%xmm3,%%xmm2	/* t13*s */\n\t"\
		"mulpd	(%%rsi),%%xmm0	/* t12*c */\n\t"\
		"mulpd	(%%rsi),%%xmm3	/* t13*c */\n\t"\
		"\n\t"\
		"addpd	%%xmm0,%%xmm2	/* rt =t12*c + t13*s */\n\t"\
		"subpd	%%xmm1,%%xmm3	/* it =t13*c - t12*s */\n\t"\
		"\n\t"\
		"movaps	     (%%rcx),%%xmm0	/* t02 */\n\t"\
		"movaps	0x010(%%rcx),%%xmm1	/* t03 */\n\t"\
		"subpd	%%xmm2,%%xmm0	/*~t12 <- t02- rt */\n\t"\
		"subpd	%%xmm3,%%xmm1	/*~t13 <- t03- it */\n\t"\
		"addpd	%%xmm2,%%xmm2	/*          2* rt */\n\t"\
		"addpd	%%xmm3,%%xmm3	/*          2* it */\n\t"\
		"addpd	%%xmm0,%%xmm2	/*~t02 <- t02+ rt */\n\t"\
		"addpd	%%xmm1,%%xmm3	/*~t03 <- t03+ it */\n\t"\
		"\n\t"\
		"/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c01): */\n\t"\
		"addq	$0x260,%%rsi	/* c01 */\n\t"\
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
		"movaps	%%xmm2,     (%%rcx)\n\t"\
		"movaps	%%xmm0,0x010(%%rdx)\n\t"\
		"movaps	%%xmm3,0x010(%%rcx)\n\t"\
		"movaps	%%xmm6,     (%%rdx)\n\t"\
		"movaps	%%xmm4,%%xmm2\n\t"\
		"movaps	%%xmm7,%%xmm0\n\t"\
		"movaps	%%xmm5,%%xmm3\n\t"\
		"movaps	%%xmm1,%%xmm6\n\t"\
		"mulpd	     (%%rsi),%%xmm4\n\t"\
		"mulpd	0x020(%%rsi),%%xmm7\n\t"\
		"mulpd	     (%%rsi),%%xmm5\n\t"\
		"mulpd	0x020(%%rsi),%%xmm1\n\t"\
		"mulpd	0x010(%%rsi),%%xmm2\n\t"\
		"mulpd	0x030(%%rsi),%%xmm0\n\t"\
		"mulpd	0x010(%%rsi),%%xmm3\n\t"\
		"mulpd	0x030(%%rsi),%%xmm6\n\t"\
		"subpd	%%xmm2,%%xmm5\n\t"\
		"subpd	%%xmm0,%%xmm1\n\t"\
		"addpd	%%xmm3,%%xmm4\n\t"\
		"addpd	%%xmm6,%%xmm7\n\t"\
		"movaps	%%xmm5,0x10(%%rax)\n\t"\
		"movaps	%%xmm1,0x10(%%rbx)\n\t"\
		"movaps	%%xmm4,    (%%rax)\n\t"\
		"movaps	%%xmm7,    (%%rbx)\n\t"\
		"addq	%%rdi,%%rax\n\t"\
		"addq	%%rdi,%%rbx\n\t"\
		"addq	$0x40,%%rsi\n\t"\
		"movaps	     (%%rcx),%%xmm4\n\t"\
		"movaps	0x010(%%rdx),%%xmm0\n\t"\
		"movaps	0x010(%%rcx),%%xmm5\n\t"\
		"movaps	     (%%rdx),%%xmm6\n\t"\
		"movaps	%%xmm4,%%xmm2\n\t"\
		"movaps	%%xmm0,%%xmm1\n\t"\
		"movaps	%%xmm5,%%xmm3\n\t"\
		"movaps	%%xmm6,%%xmm7\n\t"\
		"mulpd	     (%%rsi),%%xmm4\n\t"\
		"mulpd	0x020(%%rsi),%%xmm0\n\t"\
		"mulpd	     (%%rsi),%%xmm5\n\t"\
		"mulpd	0x020(%%rsi),%%xmm6\n\t"\
		"mulpd	0x010(%%rsi),%%xmm2\n\t"\
		"mulpd	0x030(%%rsi),%%xmm1\n\t"\
		"mulpd	0x010(%%rsi),%%xmm3\n\t"\
		"mulpd	0x030(%%rsi),%%xmm7\n\t"\
		"subpd	%%xmm2,%%xmm5\n\t"\
		"subpd	%%xmm1,%%xmm6\n\t"\
		"addpd	%%xmm3,%%xmm4\n\t"\
		"addpd	%%xmm7,%%xmm0\n\t"\
		"movaps	%%xmm5,0x10(%%rax)\n\t"\
		"movaps	%%xmm6,0x10(%%rbx)\n\t"\
		"movaps	%%xmm4,    (%%rax)\n\t"\
		"movaps	%%xmm0,    (%%rbx)\n\t"\
		"\n\t"\
		"/*...Block 3: t04,t14,t24,t34	*/\n\t"\
		"movq	%[__add0],%%rax\n\t"\
		"movslq	%[__p02],%%rsi\n\t"\
		"movslq	%[__p08],%%rbx\n\t"\
		"shlq	$3,%%rsi\n\t"\
		"shlq	$3,%%rbx\n\t"\
		"addq	%%rsi	,%%rax	/* add0 = &a[j1+p2] */\n\t"\
		"addq	%%rax	,%%rbx	/* add1 = add0+p8 */\n\t"\
		"movq	%[__r04],%%rcx\n\t"\
		"movq	%[__r14],%%rdx\n\t"\
		"movq	%[__isrt2],%%rsi\n\t"\
		"addq	$0x200,%%rcx	/* t24 */\n\t"\
		"addq	$0x200,%%rdx	/* t25 */\n\t"\
		"addq	$0x010,%%rsi	/* cc0 */\n\t"\
		"\n\t"\
		"movaps	     (%%rcx),%%xmm4	/* t24 */\n\t"\
		"movaps	0x010(%%rcx),%%xmm5	/* t25 */\n\t"\
		"movaps	     (%%rsi),%%xmm2	/* c */\n\t"\
		"movaps	0x010(%%rsi),%%xmm3	/* s */\n\t"\
		"movaps	%%xmm4,%%xmm6	/* xmm2 <- cpy t24 */\n\t"\
		"movaps	%%xmm5,%%xmm7	/* xmm3 <- cpy t25 */\n\t"\
		"\n\t"\
		"mulpd	%%xmm2,%%xmm4	/* t24*c */\n\t"\
		"mulpd	%%xmm2,%%xmm5	/* t25*c */\n\t"\
		"mulpd	%%xmm3,%%xmm6	/* t24*s */\n\t"\
		"movaps	     (%%rdx),%%xmm0	/* t34 */\n\t"\
		"mulpd	%%xmm3,%%xmm7	/* t25*s */\n\t"\
		"movaps	0x010(%%rdx),%%xmm1	/* t35 */\n\t"\
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
		"subq	$0x200,%%rcx\n\t"\
		"subq	$0x200,%%rdx\n\t"\
		"subq	$0x10,%%rsi	/* isrt2 */\n\t"\
		"movaps	     (%%rdx),%%xmm2	/* t14 */\n\t"\
		"movaps	0x010(%%rdx),%%xmm3	/* t15 */\n\t"\
		"movaps	(%%rsi),%%xmm1	/* isrt2 */\n\t"\
		"movaps	%%xmm3,%%xmm0	/* cpy t15 */\n\t"\
		"subpd	%%xmm2,%%xmm3	/*~t15=t15-t14 */\n\t"\
		"addpd	%%xmm0,%%xmm2	/*~t14=t14+t15 */\n\t"\
		"mulpd	%%xmm1,%%xmm2	/* rt */\n\t"\
		"mulpd	%%xmm1,%%xmm3	/* it */\n\t"\
		"\n\t"\
		"movaps	     (%%rcx),%%xmm0	/* t04 */\n\t"\
		"movaps	0x010(%%rcx),%%xmm1	/* t05 */\n\t"\
		"subpd	%%xmm2,%%xmm0	/*~t14 <- t04- rt */\n\t"\
		"subpd	%%xmm3,%%xmm1	/*~t15 <- t05- it */\n\t"\
		"addpd	%%xmm2,%%xmm2	/*          2* rt */\n\t"\
		"addpd	%%xmm3,%%xmm3	/*          2* it */\n\t"\
		"addpd	%%xmm0,%%xmm2	/*~t04 <- t04+ rt */\n\t"\
		"addpd	%%xmm1,%%xmm3	/*~t05 <- t05+ it */\n\t"\
		"\n\t"\
		"/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c02): */\n\t"\
		"addq	$0x170,%%rsi	/* c02 */\n\t"\
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
		"movaps	%%xmm2,     (%%rcx)\n\t"\
		"movaps	%%xmm0,0x010(%%rdx)\n\t"\
		"movaps	%%xmm3,0x010(%%rcx)\n\t"\
		"movaps	%%xmm6,     (%%rdx)\n\t"\
		"movaps	%%xmm4,%%xmm2\n\t"\
		"movaps	%%xmm7,%%xmm0\n\t"\
		"movaps	%%xmm5,%%xmm3\n\t"\
		"movaps	%%xmm1,%%xmm6\n\t"\
		"mulpd	     (%%rsi),%%xmm4\n\t"\
		"mulpd	0x020(%%rsi),%%xmm7\n\t"\
		"mulpd	     (%%rsi),%%xmm5\n\t"\
		"mulpd	0x020(%%rsi),%%xmm1\n\t"\
		"mulpd	0x010(%%rsi),%%xmm2\n\t"\
		"mulpd	0x030(%%rsi),%%xmm0\n\t"\
		"mulpd	0x010(%%rsi),%%xmm3\n\t"\
		"mulpd	0x030(%%rsi),%%xmm6\n\t"\
		"subpd	%%xmm2,%%xmm5\n\t"\
		"subpd	%%xmm0,%%xmm1\n\t"\
		"addpd	%%xmm3,%%xmm4\n\t"\
		"addpd	%%xmm6,%%xmm7\n\t"\
		"movaps	%%xmm5,0x10(%%rax)\n\t"\
		"movaps	%%xmm1,0x10(%%rbx)\n\t"\
		"movaps	%%xmm4,    (%%rax)\n\t"\
		"movaps	%%xmm7,    (%%rbx)\n\t"\
		"addq	%%rdi,%%rax\n\t"\
		"addq	%%rdi,%%rbx\n\t"\
		"addq	$0x40,%%rsi\n\t"\
		"movaps	     (%%rcx),%%xmm4\n\t"\
		"movaps	0x010(%%rdx),%%xmm0\n\t"\
		"movaps	0x010(%%rcx),%%xmm5\n\t"\
		"movaps	     (%%rdx),%%xmm6\n\t"\
		"movaps	%%xmm4,%%xmm2\n\t"\
		"movaps	%%xmm0,%%xmm1\n\t"\
		"movaps	%%xmm5,%%xmm3\n\t"\
		"movaps	%%xmm6,%%xmm7\n\t"\
		"mulpd	     (%%rsi),%%xmm4\n\t"\
		"mulpd	0x020(%%rsi),%%xmm0\n\t"\
		"mulpd	     (%%rsi),%%xmm5\n\t"\
		"mulpd	0x020(%%rsi),%%xmm6\n\t"\
		"mulpd	0x010(%%rsi),%%xmm2\n\t"\
		"mulpd	0x030(%%rsi),%%xmm1\n\t"\
		"mulpd	0x010(%%rsi),%%xmm3\n\t"\
		"mulpd	0x030(%%rsi),%%xmm7\n\t"\
		"subpd	%%xmm2,%%xmm5\n\t"\
		"subpd	%%xmm1,%%xmm6\n\t"\
		"addpd	%%xmm3,%%xmm4\n\t"\
		"addpd	%%xmm7,%%xmm0\n\t"\
		"movaps	%%xmm5,0x10(%%rax)\n\t"\
		"movaps	%%xmm6,0x10(%%rbx)\n\t"\
		"movaps	%%xmm4,    (%%rax)\n\t"\
		"movaps	%%xmm0,    (%%rbx)\n\t"\
		"\n\t"\
		"/*...Block 4: t06,t16,t26,t36	*/\n\t"\
		"movq	%[__add0],%%rax\n\t"\
		"movslq	%[__p03],%%rsi\n\t"\
		"movslq	%[__p08],%%rbx\n\t"\
		"shlq	$3,%%rsi\n\t"\
		"shlq	$3,%%rbx\n\t"\
		"addq	%%rsi	,%%rax	/* add0 = &a[j1+p3] */\n\t"\
		"addq	%%rax	,%%rbx	/* add1 = add0+p8 */\n\t"\
		"movq	%[__r06],%%rcx\n\t"\
		"movq	%[__r16],%%rdx\n\t"\
		"movq	%[__isrt2],%%rsi\n\t"\
		"addq	$0x200,%%rcx\n\t"\
		"addq	$0x200,%%rdx\n\t"\
		"addq	$0x050,%%rsi	/* cc3 */\n\t"\
		"\n\t"\
		"movaps	     (%%rcx),%%xmm4	/* t26 */\n\t"\
		"movaps	0x010(%%rcx),%%xmm5	/* t27 */\n\t"\
		"movaps	     (%%rsi),%%xmm2	/* c32_3 */\n\t"\
		"movaps	0x010(%%rsi),%%xmm3	/* s32_3 */\n\t"\
		"movaps	%%xmm4,%%xmm6	/* xmm2 <- cpy t26 */\n\t"\
		"movaps	%%xmm5,%%xmm7	/* xmm3 <- cpy t27 */\n\t"\
		"\n\t"\
		"mulpd	%%xmm2,%%xmm4	/* t26*c32_3 */\n\t"\
		"mulpd	%%xmm2,%%xmm5	/* t27*c32_3 */\n\t"\
		"mulpd	%%xmm3,%%xmm6	/* t26*s32_3 */\n\t"\
		"movaps	     (%%rdx),%%xmm0	/* t36 */\n\t"\
		"mulpd	%%xmm3,%%xmm7	/* t27*s32_3 */\n\t"\
		"movaps	0x010(%%rdx),%%xmm1	/* t37 */\n\t"\
		"subq	$0x20,%%rsi	/* cc1 */\n\t"\
		"movaps	     (%%rsi),%%xmm3	/* c32_1 */\n\t"\
		"movaps	0x010(%%rsi),%%xmm2	/* s32_1 */\n\t"\
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
		"subq	$0x200,%%rcx\n\t"\
		"subq	$0x200,%%rdx\n\t"\
		"subq	$0x20,%%rsi	/* cc0 */\n\t"\
		"movaps	     (%%rdx),%%xmm2	/* t16 */\n\t"\
		"movaps	0x010(%%rdx),%%xmm0	/* t17 */\n\t"\
		"movaps	0x010(%%rsi),%%xmm3	/* s */\n\t"\
		"movaps	%%xmm2,%%xmm1	/* cpy t16 */\n\t"\
		"mulpd	%%xmm3,%%xmm2	/* t16*s */\n\t"\
		"mulpd	%%xmm0,%%xmm3	/* s*t17 */\n\t"\
		"mulpd	(%%rsi),%%xmm1	/* t16*c */\n\t"\
		"mulpd	(%%rsi),%%xmm0	/* t17*c */\n\t"\
		"\n\t"\
		"addpd	%%xmm0,%%xmm2	/* rt =t16*s - t17*c */\n\t"\
		"subpd	%%xmm1,%%xmm3	/* it =t17*s + t16*c */\n\t"\
		"\n\t"\
		"movaps	     (%%rcx),%%xmm0	/* t06 */\n\t"\
		"movaps	0x010(%%rcx),%%xmm1	/* t07 */\n\t"\
		"subpd	%%xmm2,%%xmm0	/*~t16 <- t06- rt */\n\t"\
		"subpd	%%xmm3,%%xmm1	/*~t17 <- t07- it */\n\t"\
		"addpd	%%xmm2,%%xmm2	/*          2* rt */\n\t"\
		"addpd	%%xmm3,%%xmm3	/*          2* it */\n\t"\
		"addpd	%%xmm0,%%xmm2	/*~t06 <- t06+ rt */\n\t"\
		"addpd	%%xmm1,%%xmm3	/*~t07 <- t07+ it */\n\t"\
		"\n\t"\
		"/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c03): */\n\t"\
		"addq	$0x360,%%rsi	/* c03 */\n\t"\
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
		"movaps	%%xmm2,     (%%rcx)\n\t"\
		"movaps	%%xmm0,0x010(%%rdx)\n\t"\
		"movaps	%%xmm3,0x010(%%rcx)\n\t"\
		"movaps	%%xmm6,     (%%rdx)\n\t"\
		"movaps	%%xmm4,%%xmm2\n\t"\
		"movaps	%%xmm7,%%xmm0\n\t"\
		"movaps	%%xmm5,%%xmm3\n\t"\
		"movaps	%%xmm1,%%xmm6\n\t"\
		"mulpd	     (%%rsi),%%xmm4\n\t"\
		"mulpd	0x020(%%rsi),%%xmm7\n\t"\
		"mulpd	     (%%rsi),%%xmm5\n\t"\
		"mulpd	0x020(%%rsi),%%xmm1\n\t"\
		"mulpd	0x010(%%rsi),%%xmm2\n\t"\
		"mulpd	0x030(%%rsi),%%xmm0\n\t"\
		"mulpd	0x010(%%rsi),%%xmm3\n\t"\
		"mulpd	0x030(%%rsi),%%xmm6\n\t"\
		"subpd	%%xmm2,%%xmm5\n\t"\
		"subpd	%%xmm0,%%xmm1\n\t"\
		"addpd	%%xmm3,%%xmm4\n\t"\
		"addpd	%%xmm6,%%xmm7\n\t"\
		"movaps	%%xmm5,0x10(%%rax)\n\t"\
		"movaps	%%xmm1,0x10(%%rbx)\n\t"\
		"movaps	%%xmm4,    (%%rax)\n\t"\
		"movaps	%%xmm7,    (%%rbx)\n\t"\
		"addq	%%rdi,%%rax\n\t"\
		"addq	%%rdi,%%rbx\n\t"\
		"addq	$0x40,%%rsi	/* c0B */\n\t"\
		"movaps	     (%%rcx),%%xmm4\n\t"\
		"movaps	0x010(%%rdx),%%xmm0\n\t"\
		"movaps	0x010(%%rcx),%%xmm5\n\t"\
		"movaps	     (%%rdx),%%xmm6\n\t"\
		"movaps	%%xmm4,%%xmm2\n\t"\
		"movaps	%%xmm0,%%xmm1\n\t"\
		"movaps	%%xmm5,%%xmm3\n\t"\
		"movaps	%%xmm6,%%xmm7\n\t"\
		"mulpd	     (%%rsi),%%xmm4\n\t"\
		"mulpd	0x020(%%rsi),%%xmm0\n\t"\
		"mulpd	     (%%rsi),%%xmm5\n\t"\
		"mulpd	0x020(%%rsi),%%xmm6\n\t"\
		"mulpd	0x010(%%rsi),%%xmm2\n\t"\
		"mulpd	0x030(%%rsi),%%xmm1\n\t"\
		"mulpd	0x010(%%rsi),%%xmm3\n\t"\
		"mulpd	0x030(%%rsi),%%xmm7\n\t"\
		"subpd	%%xmm2,%%xmm5\n\t"\
		"subpd	%%xmm1,%%xmm6\n\t"\
		"addpd	%%xmm3,%%xmm4\n\t"\
		"addpd	%%xmm7,%%xmm0\n\t"\
		"movaps	%%xmm5,0x10(%%rax)\n\t"\
		"movaps	%%xmm6,0x10(%%rbx)\n\t"\
		"movaps	%%xmm4,    (%%rax)\n\t"\
		"movaps	%%xmm0,    (%%rbx)\n\t"\
		"\n\t"\
		"/*...Block 5: t08,t18,t28,t38	*/\n\t"\
		"movq	%[__add0],%%rax\n\t"\
		"movslq	%[__p04],%%rsi\n\t"\
		"movslq	%[__p08],%%rbx\n\t"\
		"shlq	$3,%%rsi\n\t"\
		"shlq	$3,%%rbx\n\t"\
		"addq	%%rsi	,%%rax	/* add0 = &a[j1+p4] */\n\t"\
		"addq	%%rax	,%%rbx	/* add1 = add0+p8 */\n\t"\
		"movq	%[__r08],%%rcx\n\t"\
		"movq	%[__r18],%%rdx\n\t"\
		"movq	%[__isrt2],%%rsi\n\t"\
		"addq	$0x200,%%rcx\n\t"\
		"addq	$0x200,%%rdx\n\t"\
		"movaps	(%%rsi),%%xmm2	/* isrt2 */\n\t"\
		"\n\t"\
		"movaps	     (%%rcx),%%xmm4	/* t28 */\n\t"\
		"movaps	0x010(%%rcx),%%xmm5	/* t29 */\n\t"\
		"movaps	     (%%rdx),%%xmm6	/* t38 */\n\t"\
		"movaps	0x010(%%rdx),%%xmm7	/* t39 */\n\t"\
		"subq	$0x200,%%rcx\n\t"\
		"subq	$0x200,%%rdx\n\t"\
		"mulpd	%%xmm2,%%xmm4\n\t"\
		"mulpd	%%xmm2,%%xmm5\n\t"\
		"mulpd	%%xmm2,%%xmm6\n\t"\
		"mulpd	%%xmm2,%%xmm7\n\t"\
		"\n\t"\
		"subpd	%%xmm4,%%xmm5	/*~t29=t29-t28*/\n\t"\
		"movaps	     (%%rcx),%%xmm0	/* t08 */\n\t"\
		"subpd	%%xmm7,%%xmm6	/* rt =t38-t39*/\n\t"\
		"movaps	0x010(%%rdx),%%xmm2	/* t19 */\n\t"\
		"addpd	%%xmm4,%%xmm4	/*       2*t28*/\n\t"\
		"movaps	0x010(%%rcx),%%xmm3	/* t09 */\n\t"\
		"addpd	%%xmm7,%%xmm7	/*       2*t39*/\n\t"\
		"movaps	     (%%rdx),%%xmm1	/* t18 */\n\t"\
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
		"addq	$0x0f0,%%rsi	/* c04 */\n\t"\
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
		"movaps	%%xmm2,     (%%rcx)\n\t"\
		"movaps	%%xmm0,0x010(%%rdx)\n\t"\
		"movaps	%%xmm3,0x010(%%rcx)\n\t"\
		"movaps	%%xmm6,     (%%rdx)\n\t"\
		"movaps	%%xmm4,%%xmm2\n\t"\
		"movaps	%%xmm7,%%xmm0\n\t"\
		"movaps	%%xmm5,%%xmm3\n\t"\
		"movaps	%%xmm1,%%xmm6\n\t"\
		"mulpd	     (%%rsi),%%xmm4\n\t"\
		"mulpd	0x020(%%rsi),%%xmm7\n\t"\
		"mulpd	     (%%rsi),%%xmm5\n\t"\
		"mulpd	0x020(%%rsi),%%xmm1\n\t"\
		"mulpd	0x010(%%rsi),%%xmm2\n\t"\
		"mulpd	0x030(%%rsi),%%xmm0\n\t"\
		"mulpd	0x010(%%rsi),%%xmm3\n\t"\
		"mulpd	0x030(%%rsi),%%xmm6\n\t"\
		"subpd	%%xmm2,%%xmm5\n\t"\
		"subpd	%%xmm0,%%xmm1\n\t"\
		"addpd	%%xmm3,%%xmm4\n\t"\
		"addpd	%%xmm6,%%xmm7\n\t"\
		"movaps	%%xmm5,0x10(%%rax)\n\t"\
		"movaps	%%xmm1,0x10(%%rbx)\n\t"\
		"movaps	%%xmm4,    (%%rax)\n\t"\
		"movaps	%%xmm7,    (%%rbx)\n\t"\
		"addq	%%rdi,%%rax\n\t"\
		"addq	%%rdi,%%rbx\n\t"\
		"addq	$0x40,%%rsi	/* c0C */\n\t"\
		"movaps	     (%%rcx),%%xmm4\n\t"\
		"movaps	0x010(%%rdx),%%xmm0\n\t"\
		"movaps	0x010(%%rcx),%%xmm5\n\t"\
		"movaps	     (%%rdx),%%xmm6\n\t"\
		"movaps	%%xmm4,%%xmm2\n\t"\
		"movaps	%%xmm0,%%xmm1\n\t"\
		"movaps	%%xmm5,%%xmm3\n\t"\
		"movaps	%%xmm6,%%xmm7\n\t"\
		"mulpd	     (%%rsi),%%xmm4\n\t"\
		"mulpd	0x020(%%rsi),%%xmm0\n\t"\
		"mulpd	     (%%rsi),%%xmm5\n\t"\
		"mulpd	0x020(%%rsi),%%xmm6\n\t"\
		"mulpd	0x010(%%rsi),%%xmm2\n\t"\
		"mulpd	0x030(%%rsi),%%xmm1\n\t"\
		"mulpd	0x010(%%rsi),%%xmm3\n\t"\
		"mulpd	0x030(%%rsi),%%xmm7\n\t"\
		"subpd	%%xmm2,%%xmm5\n\t"\
		"subpd	%%xmm1,%%xmm6\n\t"\
		"addpd	%%xmm3,%%xmm4\n\t"\
		"addpd	%%xmm7,%%xmm0\n\t"\
		"movaps	%%xmm5,0x10(%%rax)\n\t"\
		"movaps	%%xmm6,0x10(%%rbx)\n\t"\
		"movaps	%%xmm4,    (%%rax)\n\t"\
		"movaps	%%xmm0,    (%%rbx)\n\t"\
		"\n\t"\
		"/*...Block 6: t0A,t1A,t2A,t3A	*/\n\t"\
		"movq	%[__add0],%%rax\n\t"\
		"movslq	%[__p05],%%rsi\n\t"\
		"movslq	%[__p08],%%rbx\n\t"\
		"shlq	$3,%%rsi\n\t"\
		"shlq	$3,%%rbx\n\t"\
		"addq	%%rsi	,%%rax	/* add0 = &a[j1+p5] */\n\t"\
		"addq	%%rax	,%%rbx	/* add1 = add0+p8 */\n\t"\
		"movq	%[__r0A],%%rcx\n\t"\
		"movq	%[__r1A],%%rdx\n\t"\
		"movq	%[__isrt2],%%rsi\n\t"\
		"addq	$0x200,%%rcx\n\t"\
		"addq	$0x200,%%rdx\n\t"\
		"addq	$0x050,%%rsi	/* cc3 */\n\t"\
		"\n\t"\
		"movaps	     (%%rcx),%%xmm4	/* t2A */\n\t"\
		"movaps	0x010(%%rcx),%%xmm5	/* t2B */\n\t"\
		"movaps	     (%%rsi),%%xmm3	/* c32_3 */\n\t"\
		"movaps	0x010(%%rsi),%%xmm2	/* s32_3 */\n\t"\
		"movaps	%%xmm4,%%xmm6	/* xmm2 <- cpy t2A */\n\t"\
		"movaps	%%xmm5,%%xmm7	/* xmm3 <- cpy t2B */\n\t"\
		"\n\t"\
		"mulpd	%%xmm2,%%xmm4	/* t2A*s32_3 */\n\t"\
		"mulpd	%%xmm2,%%xmm5	/* t2B*s32_3 */\n\t"\
		"mulpd	%%xmm3,%%xmm6	/* t2A*c32_3 */\n\t"\
		"movaps	     (%%rdx),%%xmm0	/* t3A */\n\t"\
		"mulpd	%%xmm3,%%xmm7	/* t2B*c32_3 */\n\t"\
		"movaps	0x010(%%rdx),%%xmm1	/* t3B */\n\t"\
		"subq	$0x20,%%rsi	/* cc1 */\n\t"\
		"movaps	     (%%rsi),%%xmm2	/* c32_1 */\n\t"\
		"movaps	0x010(%%rsi),%%xmm3	/* s32_1 */\n\t"\
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
		"subq	$0x200,%%rcx\n\t"\
		"subq	$0x200,%%rdx\n\t"\
		"subq	$0x20,%%rsi	/* cc0 */\n\t"\
		"movaps	     (%%rdx),%%xmm0	/* t1A */\n\t"\
		"movaps	0x010(%%rdx),%%xmm2	/* t1B */\n\t"\
		"movaps	0x010(%%rsi),%%xmm1	/* s */\n\t"\
		"movaps	%%xmm0,%%xmm3	/* cpy t1A */\n\t"\
		"mulpd	%%xmm1,%%xmm0	/* t1A*s */\n\t"\
		"mulpd	%%xmm2,%%xmm1	/* s*t1B */\n\t"\
		"mulpd	(%%rsi),%%xmm3	/* t1A*c */\n\t"\
		"mulpd	(%%rsi),%%xmm2	/* t1B*c */\n\t"\
		"\n\t"\
		"subpd	%%xmm2,%%xmm0	/* rt =t1A*s - t1B*c */\n\t"\
		"addpd	%%xmm3,%%xmm1	/* it =t1B*s + t1A*c */\n\t"\
		"\n\t"\
		"movaps	     (%%rcx),%%xmm2	/* t0A */\n\t"\
		"movaps	0x010(%%rcx),%%xmm3	/* t0B */\n\t"\
		"subpd	%%xmm0,%%xmm2	/*~t0A <- t0A- rt */\n\t"\
		"subpd	%%xmm1,%%xmm3	/*~t0B <- t0B- it */\n\t"\
		"addpd	%%xmm0,%%xmm0	/*          2* rt */\n\t"\
		"addpd	%%xmm1,%%xmm1	/*          2* it */\n\t"\
		"addpd	%%xmm2,%%xmm0	/*~t1A <- t0A+ rt */\n\t"\
		"addpd	%%xmm3,%%xmm1	/*~t1B <- t0B+ it */\n\t"\
		"\n\t"\
		"/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c05): */\n\t"\
		"addq	$0x2e0,%%rsi	/* c05 */\n\t"\
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
		"movaps	%%xmm2,     (%%rcx)\n\t"\
		"movaps	%%xmm0,0x010(%%rdx)\n\t"\
		"movaps	%%xmm3,0x010(%%rcx)\n\t"\
		"movaps	%%xmm6,     (%%rdx)\n\t"\
		"movaps	%%xmm4,%%xmm2\n\t"\
		"movaps	%%xmm7,%%xmm0\n\t"\
		"movaps	%%xmm5,%%xmm3\n\t"\
		"movaps	%%xmm1,%%xmm6\n\t"\
		"mulpd	     (%%rsi),%%xmm4\n\t"\
		"mulpd	0x020(%%rsi),%%xmm7\n\t"\
		"mulpd	     (%%rsi),%%xmm5\n\t"\
		"mulpd	0x020(%%rsi),%%xmm1\n\t"\
		"mulpd	0x010(%%rsi),%%xmm2\n\t"\
		"mulpd	0x030(%%rsi),%%xmm0\n\t"\
		"mulpd	0x010(%%rsi),%%xmm3\n\t"\
		"mulpd	0x030(%%rsi),%%xmm6\n\t"\
		"subpd	%%xmm2,%%xmm5\n\t"\
		"subpd	%%xmm0,%%xmm1\n\t"\
		"addpd	%%xmm3,%%xmm4\n\t"\
		"addpd	%%xmm6,%%xmm7\n\t"\
		"movaps	%%xmm5,0x10(%%rax)\n\t"\
		"movaps	%%xmm1,0x10(%%rbx)\n\t"\
		"movaps	%%xmm4,    (%%rax)\n\t"\
		"movaps	%%xmm7,    (%%rbx)\n\t"\
		"addq	%%rdi,%%rax\n\t"\
		"addq	%%rdi,%%rbx\n\t"\
		"addq	$0x40,%%rsi	/* c0D */\n\t"\
		"movaps	     (%%rcx),%%xmm4\n\t"\
		"movaps	0x010(%%rdx),%%xmm0\n\t"\
		"movaps	0x010(%%rcx),%%xmm5\n\t"\
		"movaps	     (%%rdx),%%xmm6\n\t"\
		"movaps	%%xmm4,%%xmm2\n\t"\
		"movaps	%%xmm0,%%xmm1\n\t"\
		"movaps	%%xmm5,%%xmm3\n\t"\
		"movaps	%%xmm6,%%xmm7\n\t"\
		"mulpd	     (%%rsi),%%xmm4\n\t"\
		"mulpd	0x020(%%rsi),%%xmm0\n\t"\
		"mulpd	     (%%rsi),%%xmm5\n\t"\
		"mulpd	0x020(%%rsi),%%xmm6\n\t"\
		"mulpd	0x010(%%rsi),%%xmm2\n\t"\
		"mulpd	0x030(%%rsi),%%xmm1\n\t"\
		"mulpd	0x010(%%rsi),%%xmm3\n\t"\
		"mulpd	0x030(%%rsi),%%xmm7\n\t"\
		"subpd	%%xmm2,%%xmm5\n\t"\
		"subpd	%%xmm1,%%xmm6\n\t"\
		"addpd	%%xmm3,%%xmm4\n\t"\
		"addpd	%%xmm7,%%xmm0\n\t"\
		"movaps	%%xmm5,0x10(%%rax)\n\t"\
		"movaps	%%xmm6,0x10(%%rbx)\n\t"\
		"movaps	%%xmm4,    (%%rax)\n\t"\
		"movaps	%%xmm0,    (%%rbx)\n\t"\
		"\n\t"\
		"/*...Block 7: t0C,t1C,t2C,t3C	*/\n\t"\
		"movq	%[__add0],%%rax\n\t"\
		"movslq	%[__p06],%%rsi\n\t"\
		"movslq	%[__p08],%%rbx\n\t"\
		"shlq	$3,%%rsi\n\t"\
		"shlq	$3,%%rbx\n\t"\
		"addq	%%rsi	,%%rax	/* add0 = &a[j1+p6] */\n\t"\
		"addq	%%rax	,%%rbx	/* add1 = add0+p8 */\n\t"\
		"movq	%[__r0C],%%rcx\n\t"\
		"movq	%[__r1C],%%rdx\n\t"\
		"movq	%[__isrt2],%%rsi\n\t"\
		"addq	$0x200,%%rcx\n\t"\
		"addq	$0x200,%%rdx\n\t"\
		"addq	$0x010,%%rsi	/* cc0 */\n\t"\
		"\n\t"\
		"movaps	     (%%rcx),%%xmm4	/* t2C */\n\t"\
		"movaps	0x010(%%rcx),%%xmm5	/* t2D */\n\t"\
		"movaps	     (%%rsi),%%xmm3	/* c */\n\t"\
		"movaps	0x010(%%rsi),%%xmm2	/* s */\n\t"\
		"movaps	%%xmm4,%%xmm6	/* xmm2 <- cpy t2C */\n\t"\
		"movaps	%%xmm5,%%xmm7	/* xmm3 <- cpy t2D */\n\t"\
		"\n\t"\
		"mulpd	%%xmm2,%%xmm4	/* t2C*s */\n\t"\
		"mulpd	%%xmm2,%%xmm5	/* t2D*s */\n\t"\
		"mulpd	%%xmm3,%%xmm6	/* t2C*c */\n\t"\
		"movaps	     (%%rdx),%%xmm0	/* t3C */\n\t"\
		"mulpd	%%xmm3,%%xmm7	/* t2D*c */\n\t"\
		"movaps	0x010(%%rdx),%%xmm1	/* t3D */\n\t"\
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
		"subq	$0x200,%%rcx\n\t"\
		"subq	$0x200,%%rdx\n\t"\
		"subq	$0x10,%%rsi	/* isrt2 */\n\t"\
		"movaps	     (%%rdx),%%xmm0	/* t1C */\n\t"\
		"movaps	0x010(%%rdx),%%xmm1	/* t1D */\n\t"\
		"movaps	(%%rsi),%%xmm3	/* isrt2 */\n\t"\
		"movaps	%%xmm0,%%xmm2	/* cpy t1C */\n\t"\
		"subpd	%%xmm1,%%xmm0	/*~t1C=t1C-t1D */\n\t"\
		"addpd	%%xmm2,%%xmm1	/*~t1D=t1D+t1C */\n\t"\
		"mulpd	%%xmm3,%%xmm0	/* it */\n\t"\
		"mulpd	%%xmm3,%%xmm1	/* rt */\n\t"\
		"\n\t"\
		"movaps	     (%%rcx),%%xmm2	/* t0C */\n\t"\
		"movaps	0x010(%%rcx),%%xmm3	/* t0D */\n\t"\
		"subpd	%%xmm0,%%xmm2	/*~t0C <- t0C- rt */\n\t"\
		"subpd	%%xmm1,%%xmm3	/*~t0D <- t0D- it */\n\t"\
		"addpd	%%xmm0,%%xmm0	/*          2* rt */\n\t"\
		"addpd	%%xmm1,%%xmm1	/*          2* it */\n\t"\
		"addpd	%%xmm2,%%xmm0	/*~t1C <- t0C+ rt */\n\t"\
		"addpd	%%xmm3,%%xmm1	/*~t1D <- t0D+ it */\n\t"\
		"\n\t"\
		"/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c06): */\n\t"\
		"addq	$0x1f0,%%rsi	/* c06 */\n\t"\
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
		"movaps	%%xmm2,     (%%rcx)\n\t"\
		"movaps	%%xmm0,0x010(%%rdx)\n\t"\
		"movaps	%%xmm3,0x010(%%rcx)\n\t"\
		"movaps	%%xmm6,     (%%rdx)\n\t"\
		"movaps	%%xmm4,%%xmm2\n\t"\
		"movaps	%%xmm7,%%xmm0\n\t"\
		"movaps	%%xmm5,%%xmm3\n\t"\
		"movaps	%%xmm1,%%xmm6\n\t"\
		"mulpd	     (%%rsi),%%xmm4\n\t"\
		"mulpd	0x020(%%rsi),%%xmm7\n\t"\
		"mulpd	     (%%rsi),%%xmm5\n\t"\
		"mulpd	0x020(%%rsi),%%xmm1\n\t"\
		"mulpd	0x010(%%rsi),%%xmm2\n\t"\
		"mulpd	0x030(%%rsi),%%xmm0\n\t"\
		"mulpd	0x010(%%rsi),%%xmm3\n\t"\
		"mulpd	0x030(%%rsi),%%xmm6\n\t"\
		"subpd	%%xmm2,%%xmm5\n\t"\
		"subpd	%%xmm0,%%xmm1\n\t"\
		"addpd	%%xmm3,%%xmm4\n\t"\
		"addpd	%%xmm6,%%xmm7\n\t"\
		"movaps	%%xmm5,0x10(%%rax)\n\t"\
		"movaps	%%xmm1,0x10(%%rbx)\n\t"\
		"movaps	%%xmm4,    (%%rax)\n\t"\
		"movaps	%%xmm7,    (%%rbx)\n\t"\
		"addq	%%rdi,%%rax\n\t"\
		"addq	%%rdi,%%rbx\n\t"\
		"addq	$0x40,%%rsi	/* c0E */\n\t"\
		"movaps	     (%%rcx),%%xmm4\n\t"\
		"movaps	0x010(%%rdx),%%xmm0\n\t"\
		"movaps	0x010(%%rcx),%%xmm5\n\t"\
		"movaps	     (%%rdx),%%xmm6\n\t"\
		"movaps	%%xmm4,%%xmm2\n\t"\
		"movaps	%%xmm0,%%xmm1\n\t"\
		"movaps	%%xmm5,%%xmm3\n\t"\
		"movaps	%%xmm6,%%xmm7\n\t"\
		"mulpd	     (%%rsi),%%xmm4\n\t"\
		"mulpd	0x020(%%rsi),%%xmm0\n\t"\
		"mulpd	     (%%rsi),%%xmm5\n\t"\
		"mulpd	0x020(%%rsi),%%xmm6\n\t"\
		"mulpd	0x010(%%rsi),%%xmm2\n\t"\
		"mulpd	0x030(%%rsi),%%xmm1\n\t"\
		"mulpd	0x010(%%rsi),%%xmm3\n\t"\
		"mulpd	0x030(%%rsi),%%xmm7\n\t"\
		"subpd	%%xmm2,%%xmm5\n\t"\
		"subpd	%%xmm1,%%xmm6\n\t"\
		"addpd	%%xmm3,%%xmm4\n\t"\
		"addpd	%%xmm7,%%xmm0\n\t"\
		"movaps	%%xmm5,0x10(%%rax)\n\t"\
		"movaps	%%xmm6,0x10(%%rbx)\n\t"\
		"movaps	%%xmm4,    (%%rax)\n\t"\
		"movaps	%%xmm0,    (%%rbx)\n\t"\
		"\n\t"\
		"/*...Block 8: t0E,t1E,t2E,t3E	*/\n\t"\
		"movq	%[__add0],%%rax\n\t"\
		"movslq	%[__p07],%%rsi\n\t"\
		"movslq	%[__p08],%%rbx\n\t"\
		"shlq	$3,%%rsi\n\t"\
		"shlq	$3	,%%rbx	/* Pointer offset for floating doubles */\n\t"\
		"addq	%%rsi	,%%rax	/* add0 = &a[j1+p7] */\n\t"\
		"addq	%%rax	,%%rbx	/* add1 = add0+p8 */\n\t"\
		"movq	%[__r0E],%%rcx\n\t"\
		"movq	%[__r1E],%%rdx\n\t"\
		"movq	%[__isrt2],%%rsi\n\t"\
		"addq	$0x200,%%rcx\n\t"\
		"addq	$0x200,%%rdx\n\t"\
		"addq	$0x030,%%rsi	/* cc1 */\n\t"\
		"\n\t"\
		"movaps	     (%%rcx),%%xmm4	/* t2E */\n\t"\
		"movaps	0x010(%%rcx),%%xmm5	/* t2F */\n\t"\
		"movaps	     (%%rsi),%%xmm3	/* c32_1 */\n\t"\
		"movaps	0x010(%%rsi),%%xmm2	/* s32_1 */\n\t"\
		"movaps	%%xmm4,%%xmm6	/* xmm2 <- cpy t2E */\n\t"\
		"movaps	%%xmm5,%%xmm7	/* xmm3 <- cpy t2F */\n\t"\
		"\n\t"\
		"mulpd	%%xmm2,%%xmm4	/* t2E*s32_1 */\n\t"\
		"mulpd	%%xmm2,%%xmm5	/* t2F*s32_1 */\n\t"\
		"mulpd	%%xmm3,%%xmm6	/* t2E*c32_1 */\n\t"\
		"movaps	     (%%rdx),%%xmm0	/* t3E */\n\t"\
		"mulpd	%%xmm3,%%xmm7	/* t2F*c32_1 */\n\t"\
		"movaps	0x010(%%rdx),%%xmm1	/* t3F */\n\t"\
		"addq	$0x20,%%rsi	/* cc3 */\n\t"\
		"movaps	     (%%rsi),%%xmm3	/* c32_3 */\n\t"\
		"movaps	0x010(%%rsi),%%xmm2	/* s32_3 */\n\t"\
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
		"subq	$0x200,%%rcx\n\t"\
		"subq	$0x200,%%rdx\n\t"\
		"subq	$0x40,%%rsi	 /* cc0 */\n\t"\
		"movaps	     (%%rdx),%%xmm3	/* t1E */\n\t"\
		"movaps	0x010(%%rdx),%%xmm1	/* t1F */\n\t"\
		"movaps	0x010(%%rsi),%%xmm2	/* s */\n\t"\
		"movaps	%%xmm3,%%xmm0	/* cpy t1E */\n\t"\
		"mulpd	%%xmm2,%%xmm3	/* t1E*s */\n\t"\
		"mulpd	%%xmm1,%%xmm2	/* t1F*s */\n\t"\
		"mulpd	(%%rsi),%%xmm0	/* t1E*c */\n\t"\
		"mulpd	(%%rsi),%%xmm1	/* t1F*c */\n\t"\
		"\n\t"\
		"subpd	%%xmm2,%%xmm0	/* rt =t1E*c - t1F*s */\n\t"\
		"addpd	%%xmm3,%%xmm1	/* it =t1F*c + t1E*s */\n\t"\
		"\n\t"\
		"movaps	     (%%rcx),%%xmm2	/* t0E */\n\t"\
		"movaps	0x010(%%rcx),%%xmm3	/* t0F */\n\t"\
		"subpd	%%xmm0,%%xmm2	/*~t0E <- t0E- rt */\n\t"\
		"subpd	%%xmm1,%%xmm3	/*~t0F <- t0F- it */\n\t"\
		"addpd	%%xmm0,%%xmm0	/*          2* rt */\n\t"\
		"addpd	%%xmm1,%%xmm1	/*          2* it */\n\t"\
		"addpd	%%xmm2,%%xmm0	/*~t1E <- t0E+ rt */\n\t"\
		"addpd	%%xmm3,%%xmm1	/*~t1F <- t0F+ it */\n\t"\
		"/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c07): */\n\t"\
		"addq	$0x3e0,%%rsi	/* c07 */\n\t"\
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
		"movaps	%%xmm2,     (%%rcx)\n\t"\
		"movaps	%%xmm0,0x010(%%rdx)\n\t"\
		"movaps	%%xmm3,0x010(%%rcx)\n\t"\
		"movaps	%%xmm6,     (%%rdx)\n\t"\
		"movaps	%%xmm4,%%xmm2\n\t"\
		"movaps	%%xmm7,%%xmm0\n\t"\
		"movaps	%%xmm5,%%xmm3\n\t"\
		"movaps	%%xmm1,%%xmm6\n\t"\
		"mulpd	     (%%rsi),%%xmm4\n\t"\
		"mulpd	0x020(%%rsi),%%xmm7\n\t"\
		"mulpd	     (%%rsi),%%xmm5\n\t"\
		"mulpd	0x020(%%rsi),%%xmm1\n\t"\
		"mulpd	0x010(%%rsi),%%xmm2\n\t"\
		"mulpd	0x030(%%rsi),%%xmm0\n\t"\
		"mulpd	0x010(%%rsi),%%xmm3\n\t"\
		"mulpd	0x030(%%rsi),%%xmm6\n\t"\
		"subpd	%%xmm2,%%xmm5\n\t"\
		"subpd	%%xmm0,%%xmm1\n\t"\
		"addpd	%%xmm3,%%xmm4\n\t"\
		"addpd	%%xmm6,%%xmm7\n\t"\
		"movaps	%%xmm5,0x10(%%rax)\n\t"\
		"movaps	%%xmm1,0x10(%%rbx)\n\t"\
		"movaps	%%xmm4,    (%%rax)\n\t"\
		"movaps	%%xmm7,    (%%rbx)\n\t"\
		"addq	%%rdi,%%rax\n\t"\
		"addq	%%rdi,%%rbx\n\t"\
		"addq	$0x40,%%rsi	/* c0F */\n\t"\
		"movaps	     (%%rcx),%%xmm4\n\t"\
		"movaps	0x010(%%rdx),%%xmm0\n\t"\
		"movaps	0x010(%%rcx),%%xmm5\n\t"\
		"movaps	     (%%rdx),%%xmm6\n\t"\
		"movaps	%%xmm4,%%xmm2\n\t"\
		"movaps	%%xmm0,%%xmm1\n\t"\
		"movaps	%%xmm5,%%xmm3\n\t"\
		"movaps	%%xmm6,%%xmm7\n\t"\
		"mulpd	     (%%rsi),%%xmm4\n\t"\
		"mulpd	0x020(%%rsi),%%xmm0\n\t"\
		"mulpd	     (%%rsi),%%xmm5\n\t"\
		"mulpd	0x020(%%rsi),%%xmm6\n\t"\
		"mulpd	0x010(%%rsi),%%xmm2\n\t"\
		"mulpd	0x030(%%rsi),%%xmm1\n\t"\
		"mulpd	0x010(%%rsi),%%xmm3\n\t"\
		"mulpd	0x030(%%rsi),%%xmm7\n\t"\
		"subpd	%%xmm2,%%xmm5\n\t"\
		"subpd	%%xmm1,%%xmm6\n\t"\
		"addpd	%%xmm3,%%xmm4\n\t"\
		"addpd	%%xmm7,%%xmm0\n\t"\
		"movaps	%%xmm5,0x10(%%rax)\n\t"\
		"movaps	%%xmm6,0x10(%%rbx)\n\t"\
		"movaps	%%xmm4,    (%%rax)\n\t"\
		"movaps	%%xmm0,    (%%rbx)\n\t"\
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
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"	/* Clobbered registers */\
	);\
	}

  #else // USE_64BIT_ASM_STYLE = True: Deeper 64-bit-ified version of the above 32-bit ASM macros, using all of xmm0-15

	/*
	For GCC-macro version of this, use that isrt2 + 0x010,0x030,0x050 = cc0,cc1,cc3,
	and isrt2 + 0x070,0x0f0,0x170,0x1f0,0x270,0x2f0,0x370,0x3f0 = c00,04,02,06,01,05,03,07
	in order to reduce number of args to <= the GCC-allowed maximum of 30:
	*/
	#define SSE2_RADIX32_DIF_TWIDDLE(Xadd0,Xp01,Xp02,Xp03,Xp04,Xp08,Xp0C,Xp10,Xp18,Xr00)\
	{\
	__asm__ volatile (\
		"/*...Block 1: */							\n\t"\
		"movq	%[__add0],%%rax						\n\t"\
		"movslq	%[__p08],%%rbx						\n\t"\
		"movslq	%[__p10],%%rcx						\n\t	movslq	%[__p04],%%r9	\n\t"\
		"movslq	%[__p18],%%rdx						\n\t	movq	%[__r00],%%rsi	\n\t"\
		"shlq	$3,%%rbx							\n\t	shlq	$3,%%r9	\n\t"\
		"shlq	$3,%%rcx							\n\t	movq	%%rsi,%%r8	\n\t"\
		"shlq	$3,%%rdx							\n\t	movq	%%r9,%%r10	\n\t	addq	%%rax,%%r10	\n\t"\
		"addq	%%rax,%%rbx							\n\t	movq	%%r9,%%r11	\n\t	addq	%%rbx,%%r11	\n\t"\
		"addq	%%rax,%%rcx							\n\t	movq	%%r9,%%r12	\n\t	addq	%%rcx,%%r12	\n\t"\
		"addq	%%rax,%%rdx							\n\t	movq	%%r9,%%r13	\n\t	addq	%%rdx,%%r13	\n\t"\
		"/* SSE2_RADIX4_DIF_4TWIDDLE_A(r00,c10) */	\n\t	/* SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO_B(r08,c04) */\n\t"\
		"addq	$0x880,%%r8	/* two */				\n\t	movaps	    (%%r10),%%xmm8	\n\t"\
		"movaps	    (%%rax),%%xmm0					\n\t	movaps	    (%%r12),%%xmm12	\n\t"\
		"movaps	    (%%rcx),%%xmm4					\n\t	movaps	0x10(%%r10),%%xmm9	\n\t"\
		"movaps	0x10(%%rax),%%xmm1					\n\t	movaps	0x10(%%r12),%%xmm13	\n\t"\
		"movaps	0x10(%%rcx),%%xmm5					\n\t	movaps	0x4f0(%%rsi),%%xmm14	/* c04 */	\n\t"\
		"movaps	0x490(%%rsi),%%xmm2	/* c10 */		\n\t	movaps	0x500(%%rsi),%%xmm15	\n\t"\
		"movaps	0x4a0(%%rsi),%%xmm3					\n\t	movaps	%%xmm8 ,%%xmm10	\n\t"\
		"movaps	%%xmm4,%%xmm6						\n\t	movaps	%%xmm9 ,%%xmm11	\n\t"\
		"movaps	%%xmm5,%%xmm7						\n\t	mulpd	%%xmm14,%%xmm8	\n\t"\
		"mulpd	%%xmm2,%%xmm4						\n\t	mulpd	%%xmm14,%%xmm9	\n\t"\
		"mulpd	%%xmm2,%%xmm5						\n\t	mulpd	%%xmm15,%%xmm10	\n\t"\
		"mulpd	%%xmm3,%%xmm6						\n\t	mulpd	%%xmm15,%%xmm11	\n\t"\
		"movaps	%%xmm0,%%xmm2						\n\t	movaps	%%xmm12,%%xmm14	\n\t"\
		"mulpd	%%xmm3,%%xmm7						\n\t	addpd	%%xmm10,%%xmm9	\n\t"\
		"addpd	%%xmm6,%%xmm5						\n\t	movaps	%%xmm13,%%xmm15	\n\t"\
		"movaps	%%xmm1,%%xmm3						\n\t	mulpd	0x510(%%rsi),%%xmm12	/* c14 */	\n\t"\
		"subpd	%%xmm7,%%xmm4						\n\t	subpd	%%xmm11,%%xmm8	\n\t"\
		"addpd	%%xmm4,%%xmm0						\n\t	mulpd	0x510(%%rsi),%%xmm13	\n\t"\
		"addpd	%%xmm5,%%xmm1						\n\t	mulpd	0x520(%%rsi),%%xmm14	\n\t"\
		"subpd	%%xmm4,%%xmm2						\n\t	movaps	%%xmm8 ,%%xmm10	\n\t"\
		"subpd	%%xmm5,%%xmm3						\n\t	mulpd	0x520(%%rsi),%%xmm15	\n\t"\
		"movaps	    (%%rdx),%%xmm4					\n\t	addpd	%%xmm14,%%xmm13	\n\t"\
		"movaps	0x10(%%rdx),%%xmm5					\n\t	movaps	%%xmm9 ,%%xmm11	\n\t"\
		"movaps	    (%%rdx),%%xmm6					\n\t	subpd	%%xmm15,%%xmm12	\n\t"\
		"movaps	0x10(%%rdx),%%xmm7					\n\t	addpd	%%xmm12,%%xmm8	\n\t"\
		"mulpd	0x4d0(%%rsi),%%xmm4	/* c18 */		\n\t	addpd	%%xmm13,%%xmm9	\n\t"\
		"mulpd	0x4d0(%%rsi),%%xmm5					\n\t	subpd	%%xmm12,%%xmm10	\n\t"\
		"mulpd	0x4e0(%%rsi),%%xmm6					\n\t	subpd	%%xmm13,%%xmm11	\n\t"\
		"mulpd	0x4e0(%%rsi),%%xmm7					\n\t	movaps	    (%%r13),%%xmm12	\n\t"\
		"addpd	%%xmm6,%%xmm5						\n\t	movaps	0x10(%%r13),%%xmm13	\n\t"\
		"subpd	%%xmm7,%%xmm4						\n\t	movaps	    (%%r13),%%xmm14	\n\t"\
		"movaps	%%xmm5,0x10(%%rsi)					\n\t	movaps	0x10(%%r13),%%xmm15	\n\t"\
		"movaps	%%xmm4,    (%%rsi)	/* tmpstr r00 */\n\t	mulpd	0x550(%%rsi),%%xmm12	/* c1C */	\n\t"\
		"movaps	    (%%rbx),%%xmm4					\n\t	mulpd	0x550(%%rsi),%%xmm13	\n\t"\
		"movaps	0x10(%%rbx),%%xmm5					\n\t	mulpd	0x560(%%rsi),%%xmm14	\n\t"\
		"movaps	    (%%rbx),%%xmm6					\n\t	mulpd	0x560(%%rsi),%%xmm15	\n\t"\
		"movaps	0x10(%%rbx),%%xmm7					\n\t	addpd	%%xmm14,%%xmm13	\n\t"\
		"mulpd	0x4b0(%%rsi),%%xmm4	/* c08 */		\n\t	subpd	%%xmm15,%%xmm12	\n\t"\
		"mulpd	0x4b0(%%rsi),%%xmm5					\n\t	movaps	%%xmm13,0x90(%%rsi)	\n\t"\
		"mulpd	0x4c0(%%rsi),%%xmm6					\n\t	movaps	%%xmm12,0x80(%%rsi)	/* tmpstr r08 */\n\t"\
		"mulpd	0x4c0(%%rsi),%%xmm7					\n\t	movaps	    (%%r11),%%xmm12	\n\t"\
		"addpd	%%xmm6,%%xmm5						\n\t	movaps	0x10(%%r11),%%xmm13	\n\t"\
		"subpd	%%xmm7,%%xmm4						\n\t	movaps	    (%%r11),%%xmm14	\n\t"\
		"movaps	%%xmm5,%%xmm7						\n\t	movaps	0x10(%%r11),%%xmm15	\n\t"\
		"movaps	%%xmm4,%%xmm6						\n\t	mulpd	0x530(%%rsi),%%xmm12	/* c0C */	\n\t"\
		"subpd	    (%%rsi),%%xmm4	/* r00 */		\n\t	mulpd	0x530(%%rsi),%%xmm13	\n\t"\
		"subpd	0x10(%%rsi),%%xmm5					\n\t	mulpd	0x540(%%rsi),%%xmm14	\n\t"\
		"addpd	    (%%rsi),%%xmm6					\n\t	mulpd	0x540(%%rsi),%%xmm15	\n\t"\
		"addpd	0x10(%%rsi),%%xmm7					\n\t	addpd	%%xmm14,%%xmm13	\n\t"\
		"subpd	%%xmm6,%%xmm0						\n\t	subpd	%%xmm15,%%xmm12	\n\t"\
		"subpd	%%xmm5,%%xmm2						\n\t	movaps	%%xmm13,%%xmm15	\n\t"\
		"subpd	%%xmm7,%%xmm1						\n\t	movaps	%%xmm12,%%xmm14	\n\t"\
		"subpd	%%xmm4,%%xmm3						\n\t	subpd	0x80(%%rsi),%%xmm12	\n\t"\
		"/*movaps	%%xmm0,0x40(%%rsi)	*/			\n\t	subpd	0x90(%%rsi),%%xmm13	\n\t"\
		"/*movaps	%%xmm2,0x20(%%rsi)	*/			\n\t	addpd	0x80(%%rsi),%%xmm14	\n\t"\
		"/*movaps	%%xmm1,0x50(%%rsi)	*/			\n\t	addpd	0x90(%%rsi),%%xmm15	\n\t"\
		"/*movaps	%%xmm3,0x70(%%rsi)	*/			\n\t	subpd	%%xmm14,%%xmm8	\n\t"\
		"mulpd	(%%r8),%%xmm6						\n\t	subpd	%%xmm15,%%xmm9	\n\t"\
		"mulpd	(%%r8),%%xmm5						\n\t	movaps	%%xmm8 ,0xc0(%%rsi)	\n\t"\
		"mulpd	(%%r8),%%xmm7						\n\t	subpd	%%xmm13,%%xmm10	\n\t"\
		"mulpd	(%%r8),%%xmm4						\n\t/*	movaps	%%xmm9 ,0xd0(%%rsi)	*/\n\t"\
		"addpd	%%xmm0,%%xmm6						\n\t	subpd	%%xmm12,%%xmm11	\n\t"\
		"addpd	%%xmm2,%%xmm5						\n\t	mulpd	(%%r8),%%xmm14	\n\t"\
		"addpd	%%xmm1,%%xmm7						\n\t	mulpd	(%%r8),%%xmm13	\n\t"\
		"addpd	%%xmm3,%%xmm4						\n\t	mulpd	(%%r8),%%xmm15	\n\t"\
		"/*movaps	%%xmm6,    (%%rsi)	*/			\n\t	mulpd	(%%r8),%%xmm12	\n\t"\
		"/*movaps	%%xmm5,0x60(%%rsi)	*/			\n\t	addpd	%%xmm8 ,%%xmm14	\n\t"\
		"/*movaps	%%xmm7,0x10(%%rsi)	*/			\n\t	addpd	%%xmm10,%%xmm13	\n\t"\
		"/*movaps	%%xmm4,0x30(%%rsi)	*/			\n\t	addpd	%%xmm9 ,%%xmm15	\n\t"\
		"											\n\t	addpd	%%xmm11,%%xmm12	\n\t"\
		"													movaps	%%xmm14,0x80(%%rsi)	\n\t"\
		"													movaps	%%xmm15,0x90(%%rsi)	\n\t"\
		"													movaps	0x400(%%rsi),%%xmm8	/* isrt2 */	\n\t"\
		"													movaps	%%xmm10,%%xmm14	\n\t"\
		"													movaps	%%xmm13,%%xmm15	\n\t"\
		"													subpd	%%xmm12,%%xmm10	\n\t"\
		"													subpd	%%xmm11,%%xmm13	\n\t"\
		"													addpd	%%xmm12,%%xmm14	\n\t"\
		"													addpd	%%xmm11,%%xmm15	\n\t"\
		"													mulpd	%%xmm8 ,%%xmm10	\n\t"\
		"													mulpd	%%xmm8 ,%%xmm13	\n\t"\
		"													mulpd	%%xmm8 ,%%xmm14	\n\t"\
		"subpd	%%xmm10,%%xmm2	/* Since xmm10 output is first-done of 2nd set, move the computation of xmm2/10 outputs up so can use xmm2 for 2.0 doubling constant below */\n\t"\
		"													mulpd	%%xmm8 ,%%xmm15	\n\t"\
		"												/*	movaps	%%xmm10,0xa0(%%rsi)	*/\n\t"\
		"												/*	movaps	%%xmm13,0xe0(%%rsi)	*/\n\t"\
		"												/*	movaps	%%xmm14,0xb0(%%rsi)	*/\n\t"\
		"												/*	movaps	%%xmm15,0xf0(%%rsi)	*/\n\t"\
		"/* SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r00) */	\n\t"\
		"movaps	0x80(%%rsi),%%xmm11		\n\t"\
		"movaps	0x90(%%rsi),%%xmm12		\n\t"\
		"						mulpd	(%%r8),%%xmm10	\n\t"\
		"movaps	0xc0(%%rsi),%%xmm8		\n\t"\
		"subpd	%%xmm11,%%xmm6						\n\t/*	subpd	%%xmm10,%%xmm2	*/\n\t"\
		"subpd	%%xmm9 ,%%xmm0						\n\t	subpd	%%xmm15,%%xmm5	\n\t"\
		"						addpd	%%xmm2 ,%%xmm10	\n\t"\
		"						movaps	%%xmm2 ,0xa0(%%rsi)	\n\t"\
		"subpd	%%xmm12,%%xmm7						\n\t	subpd	%%xmm14,%%xmm4	\n\t"\
		"movaps	(%%r8),%%xmm2	/* 2.0 ... We could use add-in-place for doubling, but want to load-balance add/mul here. */	\n\t"\
		"subpd	%%xmm8 ,%%xmm1						\n\t	subpd	%%xmm13,%%xmm3	\n\t"\
		"mulpd	%%xmm2 ,%%xmm11						\n\t/*	addpd	%%xmm10,%%xmm10	*/\n\t"\
		"mulpd	%%xmm2 ,%%xmm9						\n\t	mulpd	%%xmm2 ,%%xmm15	\n\t"\
		"mulpd	%%xmm2 ,%%xmm12						\n\t	mulpd	%%xmm2 ,%%xmm14	\n\t"\
		"mulpd	%%xmm2 ,%%xmm8						\n\t	mulpd	%%xmm2 ,%%xmm13	\n\t"\
		"addpd	%%xmm6 ,%%xmm11						\n\t/*	addpd	%%xmm2 ,%%xmm10	*/\n\t"\
		"addpd	%%xmm0 ,%%xmm9						\n\t	addpd	%%xmm5 ,%%xmm15	\n\t"\
		"addpd	%%xmm7 ,%%xmm12						\n\t	addpd	%%xmm4 ,%%xmm14	\n\t"\
		"addpd	%%xmm1 ,%%xmm8						\n\t	addpd	%%xmm3 ,%%xmm13	\n\t"\
		"movaps	%%xmm6 ,0x80(%%rsi)					\n\t/*	movaps	%%xmm2 ,0xa0(%%rsi)	*/\n\t"\
		"movaps	%%xmm0 ,0x40(%%rsi)					\n\t	movaps	%%xmm5 ,0x60(%%rsi)	\n\t"\
		"movaps	%%xmm7 ,0x90(%%rsi)					\n\t	movaps	%%xmm4 ,0xb0(%%rsi)	\n\t"\
		"movaps	%%xmm1 ,0xd0(%%rsi)					\n\t	movaps	%%xmm3 ,0xf0(%%rsi)	\n\t"\
		"movaps	%%xmm11,    (%%rsi)					\n\t	movaps	%%xmm10,0x20(%%rsi)	\n\t"\
		"movaps	%%xmm9 ,0xc0(%%rsi)					\n\t	movaps	%%xmm15,0xe0(%%rsi)	\n\t"\
		"movaps	%%xmm12,0x10(%%rsi)					\n\t	movaps	%%xmm14,0x30(%%rsi)	\n\t"\
		"movaps	%%xmm8 ,0x50(%%rsi)					\n\t	movaps	%%xmm13,0x70(%%rsi)	\n\t"\
		"\n\t"\
		"/***************************************/\n\t"\
		"\n\t"\
		"/*...Block 2: */\n\t"\
		"movslq	%[__p02],%%rdi						\n\t	movslq	%[__p04],%%r9	\n\t"\
		"shlq	$3,%%rdi							\n\t	shlq	$3,%%r9	\n\t"\
		"addq	%%rdi,%%rax	/* &a[j1+p2] */			\n\t	movq	%%r9,%%r10	\n\t	addq	%%rax,%%r10	\n\t"\
		"addq	%%rdi,%%rbx							\n\t	movq	%%r9,%%r11	\n\t	addq	%%rbx,%%r11	\n\t"\
		"addq	%%rdi,%%rcx							\n\t	movq	%%r9,%%r12	\n\t	addq	%%rcx,%%r12	\n\t"\
		"addq	%%rdi,%%rdx							\n\t	movq	%%r9,%%r13	\n\t	addq	%%rdx,%%r13	\n\t"\
		"/* SSE2_RADIX4_DIF_4TWIDDLE_B(r10,c02) */	\n\t	/* SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO_B(r18,c06) */\n\t"\
		"addq	$0x100,%%rsi	/* r10 */			\n\t	movaps	    (%%r10),%%xmm8	\n\t"\
		"movaps	    (%%rax),%%xmm0					\n\t	movaps	    (%%r12),%%xmm12	\n\t"\
		"movaps	    (%%rcx),%%xmm4					\n\t	movaps	0x10(%%r10),%%xmm9	\n\t"\
		"movaps	0x10(%%rax),%%xmm1					\n\t	movaps	0x10(%%r12),%%xmm13	\n\t"\
		"movaps	0x10(%%rcx),%%xmm5					\n\t	movaps	0x4f0(%%rsi),%%xmm14	/* c06 */	\n\t"\
		"movaps	0x470(%%rsi),%%xmm6	/* c02 */		\n\t	movaps	0x500(%%rsi),%%xmm15	\n\t"\
		"movaps	0x480(%%rsi),%%xmm7					\n\t	movaps	%%xmm8 ,%%xmm10	\n\t"\
		"movaps	%%xmm0,%%xmm2						\n\t	movaps	%%xmm9 ,%%xmm11	\n\t"\
		"movaps	%%xmm1,%%xmm3						\n\t	mulpd	%%xmm14,%%xmm8	\n\t"\
		"mulpd   %%xmm6,%%xmm0						\n\t	mulpd	%%xmm14,%%xmm9	\n\t"\
		"mulpd   %%xmm6,%%xmm1						\n\t	mulpd	%%xmm15,%%xmm10	\n\t"\
		"mulpd   %%xmm7,%%xmm2						\n\t	mulpd	%%xmm15,%%xmm11	\n\t"\
		"mulpd   %%xmm7,%%xmm3						\n\t	movaps	%%xmm12,%%xmm14	\n\t"\
		"movaps	%%xmm4,%%xmm6						\n\t	addpd	%%xmm10,%%xmm9	\n\t"\
		"addpd   %%xmm2,%%xmm1						\n\t	movaps	%%xmm13,%%xmm15	\n\t"\
		"movaps	%%xmm5,%%xmm7						\n\t	mulpd	0x510(%%rsi),%%xmm12	/* c16 */	\n\t"\
		"mulpd   0x490(%%rsi),%%xmm4	/* c12 */	\n\t	subpd	%%xmm11,%%xmm8	\n\t"\
		"subpd   %%xmm3,%%xmm0						\n\t	mulpd	0x510(%%rsi),%%xmm13	\n\t"\
		"mulpd   0x490(%%rsi),%%xmm5				\n\t	mulpd	0x520(%%rsi),%%xmm14	\n\t"\
		"mulpd   0x4a0(%%rsi),%%xmm6				\n\t	movaps	%%xmm8 ,%%xmm10	\n\t"\
		"movaps	%%xmm0,%%xmm2						\n\t	mulpd	0x520(%%rsi),%%xmm15	\n\t"\
		"mulpd	0x4a0(%%rsi),%%xmm7					\n\t	addpd	%%xmm14,%%xmm13	\n\t"\
		"addpd	%%xmm6,%%xmm5						\n\t	movaps	%%xmm9 ,%%xmm11	\n\t"\
		"movaps	%%xmm1,%%xmm3						\n\t	subpd	%%xmm15,%%xmm12	\n\t"\
		"subpd	%%xmm7,%%xmm4						\n\t	addpd	%%xmm12,%%xmm8	\n\t"\
		"addpd	%%xmm4,%%xmm0						\n\t	addpd	%%xmm13,%%xmm9	\n\t"\
		"addpd	%%xmm5,%%xmm1						\n\t	subpd	%%xmm12,%%xmm10	\n\t"\
		"subpd	%%xmm4,%%xmm2						\n\t	subpd	%%xmm13,%%xmm11	\n\t"\
		"subpd	%%xmm5,%%xmm3						\n\t	movaps	    (%%r13),%%xmm12	\n\t"\
		"movaps	    (%%rdx),%%xmm4					\n\t	movaps	0x10(%%r13),%%xmm13	\n\t"\
		"movaps	0x10(%%rdx),%%xmm5					\n\t	movaps	    (%%r13),%%xmm14	\n\t"\
		"movaps	    (%%rdx),%%xmm6					\n\t	movaps	0x10(%%r13),%%xmm15	\n\t"\
		"movaps	0x10(%%rdx),%%xmm7					\n\t	mulpd	0x550(%%rsi),%%xmm12	/* c1E */	\n\t"\
		"mulpd	0x4d0(%%rsi),%%xmm4	/* c1A */		\n\t	mulpd	0x550(%%rsi),%%xmm13	\n\t"\
		"mulpd	0x4d0(%%rsi),%%xmm5					\n\t	mulpd	0x560(%%rsi),%%xmm14	\n\t"\
		"mulpd	0x4e0(%%rsi),%%xmm6					\n\t	mulpd	0x560(%%rsi),%%xmm15	\n\t"\
		"mulpd	0x4e0(%%rsi),%%xmm7					\n\t	addpd	%%xmm14,%%xmm13	\n\t"\
		"addpd	%%xmm6,%%xmm5						\n\t	subpd	%%xmm15,%%xmm12	\n\t"\
		"subpd	%%xmm7,%%xmm4						\n\t	movaps	%%xmm13,0x90(%%rsi)	\n\t"\
		"movaps	%%xmm5,0x010(%%rsi)					\n\t	movaps	%%xmm12,0x80(%%rsi)	/* tmpstr r08 */\n\t"\
		"movaps	%%xmm4,     (%%rsi)					\n\t	movaps	    (%%r11),%%xmm12	\n\t"\
		"movaps	    (%%rbx),%%xmm4					\n\t	movaps	0x10(%%r11),%%xmm13	\n\t"\
		"movaps	0x10(%%rbx),%%xmm5					\n\t	movaps	    (%%r11),%%xmm14	\n\t"\
		"movaps	    (%%rbx),%%xmm6					\n\t	movaps	0x10(%%r11),%%xmm15	\n\t"\
		"movaps	0x10(%%rbx),%%xmm7					\n\t	mulpd	0x530(%%rsi),%%xmm12	/* c0E */	\n\t"\
		"mulpd	0x4b0(%%rsi),%%xmm4	/* c0A */		\n\t	mulpd	0x530(%%rsi),%%xmm13	\n\t"\
		"mulpd	0x4b0(%%rsi),%%xmm5					\n\t	mulpd	0x540(%%rsi),%%xmm14	\n\t"\
		"mulpd	0x4c0(%%rsi),%%xmm6					\n\t	mulpd	0x540(%%rsi),%%xmm15	\n\t"\
		"mulpd	0x4c0(%%rsi),%%xmm7					\n\t	addpd	%%xmm14,%%xmm13	\n\t"\
		"addpd	%%xmm6,%%xmm5						\n\t	subpd	%%xmm15,%%xmm12	\n\t"\
		"subpd	%%xmm7,%%xmm4						\n\t	movaps	%%xmm13,%%xmm15	\n\t"\
		"movaps	%%xmm5,%%xmm7						\n\t	movaps	%%xmm12,%%xmm14	\n\t"\
		"movaps	%%xmm4,%%xmm6						\n\t	subpd	0x80(%%rsi),%%xmm12	\n\t"\
		"subpd	    (%%rsi),%%xmm4					\n\t	subpd	0x90(%%rsi),%%xmm13	\n\t"\
		"subpd	0x10(%%rsi),%%xmm5					\n\t	addpd	0x80(%%rsi),%%xmm14	\n\t"\
		"addpd	    (%%rsi),%%xmm6					\n\t	addpd	0x90(%%rsi),%%xmm15	\n\t"\
		"addpd	0x10(%%rsi),%%xmm7					\n\t	subpd	%%xmm14,%%xmm8	\n\t"\
		"subpd	%%xmm6,%%xmm0						\n\t	subpd	%%xmm15,%%xmm9	\n\t"\
		"subpd	%%xmm5,%%xmm2						\n\t	movaps	%%xmm8 ,0xc0(%%rsi)	\n\t"\
		"subpd	%%xmm7,%%xmm1						\n\t	subpd	%%xmm13,%%xmm10	\n\t"\
		"subpd	%%xmm4,%%xmm3						\n\t/*	movaps	%%xmm9 ,0xd0(%%rsi)	*/\n\t"\
		"movaps	%%xmm0,0x040(%%rsi)					\n\t	subpd	%%xmm12,%%xmm11	\n\t"\
		"movaps	%%xmm2,0x020(%%rsi)					\n\t	mulpd	(%%r8),%%xmm14	\n\t"\
		"movaps	%%xmm1,0x050(%%rsi)					\n\t	mulpd	(%%r8),%%xmm13	\n\t"\
		"movaps	%%xmm3,0x070(%%rsi)					\n\t	mulpd	(%%r8),%%xmm15	\n\t"\
		"mulpd	(%%r8),%%xmm6						\n\t	mulpd	(%%r8),%%xmm12	\n\t"\
		"mulpd	(%%r8),%%xmm5						\n\t	addpd	%%xmm8 ,%%xmm14	\n\t"\
		"mulpd	(%%r8),%%xmm7						\n\t	addpd	%%xmm10,%%xmm13	\n\t"\
		"mulpd	(%%r8),%%xmm4						\n\t	addpd	%%xmm9 ,%%xmm15	\n\t"\
		"addpd	%%xmm0,%%xmm6						\n\t	addpd	%%xmm11,%%xmm12	\n\t"\
		"addpd	%%xmm2,%%xmm5						\n\t	movaps	%%xmm14,0x80(%%rsi)	\n\t"\
		"addpd	%%xmm1,%%xmm7						\n\t	movaps	%%xmm15,0x90(%%rsi)	\n\t"\
		"addpd	%%xmm3,%%xmm4						\n\t	movaps	0x300(%%rsi),%%xmm8	/* isrt2 */	\n\t"\
		"movaps	%%xmm6,     (%%rsi)					\n\t	movaps	%%xmm10,%%xmm14	\n\t"\
		"movaps	%%xmm5,0x060(%%rsi)					\n\t	movaps	%%xmm13,%%xmm15	\n\t"\
		"movaps	%%xmm7,0x010(%%rsi)					\n\t	subpd	%%xmm12,%%xmm10	\n\t"\
		"movaps	%%xmm4,0x030(%%rsi)					\n\t	subpd	%%xmm11,%%xmm13	\n\t"\
		"													addpd	%%xmm12,%%xmm14	\n\t"\
		"													addpd	%%xmm11,%%xmm15	\n\t"\
		"													mulpd	%%xmm8 ,%%xmm10	\n\t"\
		"													mulpd	%%xmm8 ,%%xmm13	\n\t"\
		"													mulpd	%%xmm8 ,%%xmm14	\n\t"\
		"						subpd	%%xmm10,%%xmm2	\n\t"\
		"													mulpd	%%xmm8 ,%%xmm15	\n\t"\
		"												/*	movaps	%%xmm10,0xa0(%%rsi)	*/\n\t"\
		"												/*	movaps	%%xmm13,0xe0(%%rsi)	*/\n\t"\
		"												/*	movaps	%%xmm14,0xb0(%%rsi)	*/\n\t"\
		"												/*	movaps	%%xmm15,0xf0(%%rsi)	*/\n\t"\
		"/* SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r10) */	\n\t"\
		"movaps	0x80(%%rsi),%%xmm11		\n\t"\
		"movaps	0x90(%%rsi),%%xmm12		\n\t"\
		"						mulpd	(%%r8),%%xmm10	\n\t"\
		"movaps	0xc0(%%rsi),%%xmm8		\n\t"\
		"subpd	%%xmm11,%%xmm6						\n\t/*	subpd	%%xmm10,%%xmm2	*/\n\t"\
		"subpd	%%xmm9 ,%%xmm0						\n\t	subpd	%%xmm15,%%xmm5	\n\t"\
		"						addpd	%%xmm2 ,%%xmm10	\n\t"\
		"						movaps	%%xmm2 ,0xa0(%%rsi)	\n\t"\
		"subpd	%%xmm12,%%xmm7						\n\t	subpd	%%xmm14,%%xmm4	\n\t"\
		"						movaps	(%%r8),%%xmm2	/* 2.0 */	\n\t"\
		"subpd	%%xmm8 ,%%xmm1						\n\t	subpd	%%xmm13,%%xmm3	\n\t"\
		"mulpd	%%xmm2 ,%%xmm11						\n\t/*	addpd	%%xmm10,%%xmm10	*/\n\t"\
		"mulpd	%%xmm2 ,%%xmm9						\n\t	mulpd	%%xmm2 ,%%xmm15	\n\t"\
		"mulpd	%%xmm2 ,%%xmm12						\n\t	mulpd	%%xmm2 ,%%xmm14	\n\t"\
		"mulpd	%%xmm2 ,%%xmm8						\n\t	mulpd	%%xmm2 ,%%xmm13	\n\t"\
		"addpd	%%xmm6 ,%%xmm11						\n\t/*	addpd	%%xmm2 ,%%xmm10	*/\n\t"\
		"addpd	%%xmm0 ,%%xmm9						\n\t	addpd	%%xmm5 ,%%xmm15	\n\t"\
		"addpd	%%xmm7 ,%%xmm12						\n\t	addpd	%%xmm4 ,%%xmm14	\n\t"\
		"addpd	%%xmm1 ,%%xmm8						\n\t	addpd	%%xmm3 ,%%xmm13	\n\t"\
		"movaps	%%xmm6 ,0x80(%%rsi)					\n\t/*	movaps	%%xmm2 ,0xa0(%%rsi)	*/\n\t"\
		"movaps	%%xmm0 ,0x40(%%rsi)					\n\t	movaps	%%xmm5 ,0x60(%%rsi)	\n\t"\
		"movaps	%%xmm7 ,0x90(%%rsi)					\n\t	movaps	%%xmm4 ,0xb0(%%rsi)	\n\t"\
		"movaps	%%xmm1 ,0xd0(%%rsi)					\n\t	movaps	%%xmm3 ,0xf0(%%rsi)	\n\t"\
		"movaps	%%xmm11,    (%%rsi)					\n\t	movaps	%%xmm10,0x20(%%rsi)	\n\t"\
		"movaps	%%xmm9 ,0xc0(%%rsi)					\n\t	movaps	%%xmm15,0xe0(%%rsi)	\n\t"\
		"movaps	%%xmm12,0x10(%%rsi)					\n\t	movaps	%%xmm14,0x30(%%rsi)	\n\t"\
		"movaps	%%xmm8 ,0x50(%%rsi)					\n\t	movaps	%%xmm13,0x70(%%rsi)	\n\t"\
		"\n\t"\
		"/***************************************/\n\t"\
		"\n\t"\
		"/*...Block 3: */	\n\t"\
		"movq	%[__add0],%%rax	\n\t"\
		"movslq	%[__p01],%%rdi	/* Do things this way [rather than repeatedly adding p1] since array-padding scheme means p2 == p+p1 not guaranteed. */	\n\t"\
		"movslq	%[__p08],%%rbx	/* Can't get these via simple load-one-and-shift-as-needed due to array padding scheme */	\n\t"\
		"movslq	%[__p10],%%rcx	\n\t"\
		"movslq	%[__p18],%%rdx	\n\t"\
		"shlq	$3,%%rdi	\n\t"\
		"shlq	$3,%%rbx	\n\t"\
		"shlq	$3,%%rcx							\n\t	movslq	%[__p04],%%r9	\n\t"\
		"shlq	$3,%%rdx							\n\t	shlq	$3,%%r9	\n\t"\
		"addq	%%rdi,%%rax	/* &a[j1+p1] */			\n\t	movq	%%r9,%%r10	\n\t	addq	%%rax,%%r10	\n\t"\
		"addq	%%rax,%%rbx							\n\t	movq	%%r9,%%r11	\n\t	addq	%%rbx,%%r11	\n\t"\
		"addq	%%rax,%%rcx							\n\t	movq	%%r9,%%r12	\n\t	addq	%%rcx,%%r12	\n\t"\
		"addq	%%rax,%%rdx							\n\t	movq	%%r9,%%r13	\n\t	addq	%%rdx,%%r13	\n\t"\
		"/* SSE2_RADIX4_DIF_4TWIDDLE_B(r20,c01) */	\n\t	/* SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO_B(r28,c05) */\n\t"\
		"addq	$0x100,%%rsi	/* r20 */			\n\t	movaps	    (%%r10),%%xmm8	\n\t"\
		"movaps	    (%%rax),%%xmm0					\n\t	movaps	    (%%r12),%%xmm12	\n\t"\
		"movaps	    (%%rcx),%%xmm4					\n\t	movaps	0x10(%%r10),%%xmm9	\n\t"\
		"movaps	0x10(%%rax),%%xmm1					\n\t	movaps	0x10(%%r12),%%xmm13	\n\t"\
		"movaps	0x10(%%rcx),%%xmm5					\n\t	movaps	0x4f0(%%rsi),%%xmm14	/* c05 */	\n\t"\
		"movaps	0x470(%%rsi),%%xmm6	/* c01 */		\n\t	movaps	0x500(%%rsi),%%xmm15	\n\t"\
		"movaps	0x480(%%rsi),%%xmm7					\n\t	movaps	%%xmm8 ,%%xmm10	\n\t"\
		"movaps	%%xmm0,%%xmm2						\n\t	movaps	%%xmm9 ,%%xmm11	\n\t"\
		"movaps	%%xmm1,%%xmm3						\n\t	mulpd	%%xmm14,%%xmm8	\n\t"\
		"mulpd   %%xmm6,%%xmm0						\n\t	mulpd	%%xmm14,%%xmm9	\n\t"\
		"mulpd   %%xmm6,%%xmm1						\n\t	mulpd	%%xmm15,%%xmm10	\n\t"\
		"mulpd   %%xmm7,%%xmm2						\n\t	mulpd	%%xmm15,%%xmm11	\n\t"\
		"mulpd   %%xmm7,%%xmm3						\n\t	movaps	%%xmm12,%%xmm14	\n\t"\
		"movaps	%%xmm4,%%xmm6						\n\t	addpd	%%xmm10,%%xmm9	\n\t"\
		"addpd   %%xmm2,%%xmm1						\n\t	movaps	%%xmm13,%%xmm15	\n\t"\
		"movaps	%%xmm5,%%xmm7						\n\t	mulpd	0x510(%%rsi),%%xmm12	/* c15 */	\n\t"\
		"mulpd   0x490(%%rsi),%%xmm4	/* c11 */	\n\t	subpd	%%xmm11,%%xmm8	\n\t"\
		"subpd   %%xmm3,%%xmm0						\n\t	mulpd	0x510(%%rsi),%%xmm13	\n\t"\
		"mulpd   0x490(%%rsi),%%xmm5				\n\t	mulpd	0x520(%%rsi),%%xmm14	\n\t"\
		"mulpd   0x4a0(%%rsi),%%xmm6				\n\t	movaps	%%xmm8 ,%%xmm10	\n\t"\
		"movaps	%%xmm0,%%xmm2						\n\t	mulpd	0x520(%%rsi),%%xmm15	\n\t"\
		"mulpd	0x4a0(%%rsi),%%xmm7					\n\t	addpd	%%xmm14,%%xmm13	\n\t"\
		"addpd	%%xmm6,%%xmm5						\n\t	movaps	%%xmm9 ,%%xmm11	\n\t"\
		"movaps	%%xmm1,%%xmm3						\n\t	subpd	%%xmm15,%%xmm12	\n\t"\
		"subpd	%%xmm7,%%xmm4						\n\t	addpd	%%xmm12,%%xmm8	\n\t"\
		"addpd	%%xmm4,%%xmm0						\n\t	addpd	%%xmm13,%%xmm9	\n\t"\
		"addpd	%%xmm5,%%xmm1						\n\t	subpd	%%xmm12,%%xmm10	\n\t"\
		"subpd	%%xmm4,%%xmm2						\n\t	subpd	%%xmm13,%%xmm11	\n\t"\
		"subpd	%%xmm5,%%xmm3						\n\t	movaps	    (%%r13),%%xmm12	\n\t"\
		"movaps	    (%%rdx),%%xmm4					\n\t	movaps	0x10(%%r13),%%xmm13	\n\t"\
		"movaps	0x10(%%rdx),%%xmm5					\n\t	movaps	    (%%r13),%%xmm14	\n\t"\
		"movaps	    (%%rdx),%%xmm6					\n\t	movaps	0x10(%%r13),%%xmm15	\n\t"\
		"movaps	0x10(%%rdx),%%xmm7					\n\t	mulpd	0x550(%%rsi),%%xmm12	/* c1D */	\n\t"\
		"mulpd	0x4d0(%%rsi),%%xmm4	/* c19 */		\n\t	mulpd	0x550(%%rsi),%%xmm13	\n\t"\
		"mulpd	0x4d0(%%rsi),%%xmm5					\n\t	mulpd	0x560(%%rsi),%%xmm14	\n\t"\
		"mulpd	0x4e0(%%rsi),%%xmm6					\n\t	mulpd	0x560(%%rsi),%%xmm15	\n\t"\
		"mulpd	0x4e0(%%rsi),%%xmm7					\n\t	addpd	%%xmm14,%%xmm13	\n\t"\
		"addpd	%%xmm6,%%xmm5						\n\t	subpd	%%xmm15,%%xmm12	\n\t"\
		"subpd	%%xmm7,%%xmm4						\n\t	movaps	%%xmm13,0x90(%%rsi)	\n\t"\
		"movaps	%%xmm5,0x010(%%rsi)					\n\t	movaps	%%xmm12,0x80(%%rsi)	/* tmpstr r08 */\n\t"\
		"movaps	%%xmm4,     (%%rsi)					\n\t	movaps	    (%%r11),%%xmm12	\n\t"\
		"movaps	    (%%rbx),%%xmm4					\n\t	movaps	0x10(%%r11),%%xmm13	\n\t"\
		"movaps	0x10(%%rbx),%%xmm5					\n\t	movaps	    (%%r11),%%xmm14	\n\t"\
		"movaps	    (%%rbx),%%xmm6					\n\t	movaps	0x10(%%r11),%%xmm15	\n\t"\
		"movaps	0x10(%%rbx),%%xmm7					\n\t	mulpd	0x530(%%rsi),%%xmm12	/* c0D */	\n\t"\
		"mulpd	0x4b0(%%rsi),%%xmm4	/* c01 */		\n\t	mulpd	0x530(%%rsi),%%xmm13	\n\t"\
		"mulpd	0x4b0(%%rsi),%%xmm5					\n\t	mulpd	0x540(%%rsi),%%xmm14	\n\t"\
		"mulpd	0x4c0(%%rsi),%%xmm6					\n\t	mulpd	0x540(%%rsi),%%xmm15	\n\t"\
		"mulpd	0x4c0(%%rsi),%%xmm7					\n\t	addpd	%%xmm14,%%xmm13	\n\t"\
		"addpd	%%xmm6,%%xmm5						\n\t	subpd	%%xmm15,%%xmm12	\n\t"\
		"subpd	%%xmm7,%%xmm4						\n\t	movaps	%%xmm13,%%xmm15	\n\t"\
		"movaps	%%xmm5,%%xmm7						\n\t	movaps	%%xmm12,%%xmm14	\n\t"\
		"movaps	%%xmm4,%%xmm6						\n\t	subpd	0x80(%%rsi),%%xmm12	\n\t"\
		"subpd	    (%%rsi),%%xmm4					\n\t	subpd	0x90(%%rsi),%%xmm13	\n\t"\
		"subpd	0x10(%%rsi),%%xmm5					\n\t	addpd	0x80(%%rsi),%%xmm14	\n\t"\
		"addpd	    (%%rsi),%%xmm6					\n\t	addpd	0x90(%%rsi),%%xmm15	\n\t"\
		"addpd	0x10(%%rsi),%%xmm7					\n\t	subpd	%%xmm14,%%xmm8	\n\t"\
		"subpd	%%xmm6,%%xmm0						\n\t	subpd	%%xmm15,%%xmm9	\n\t"\
		"subpd	%%xmm5,%%xmm2						\n\t	movaps	%%xmm8 ,0xc0(%%rsi)	\n\t"\
		"subpd	%%xmm7,%%xmm1						\n\t	subpd	%%xmm13,%%xmm10	\n\t"\
		"subpd	%%xmm4,%%xmm3						\n\t/*	movaps	%%xmm9 ,0xd0(%%rsi)	*/\n\t"\
		"movaps	%%xmm0,0x040(%%rsi)					\n\t	subpd	%%xmm12,%%xmm11	\n\t"\
		"movaps	%%xmm2,0x020(%%rsi)					\n\t	mulpd	(%%r8),%%xmm14	\n\t"\
		"movaps	%%xmm1,0x050(%%rsi)					\n\t	mulpd	(%%r8),%%xmm13	\n\t"\
		"movaps	%%xmm3,0x070(%%rsi)					\n\t	mulpd	(%%r8),%%xmm15	\n\t"\
		"mulpd	(%%r8),%%xmm6						\n\t	mulpd	(%%r8),%%xmm12	\n\t"\
		"mulpd	(%%r8),%%xmm5						\n\t	addpd	%%xmm8 ,%%xmm14	\n\t"\
		"mulpd	(%%r8),%%xmm7						\n\t	addpd	%%xmm10,%%xmm13	\n\t"\
		"mulpd	(%%r8),%%xmm4						\n\t	addpd	%%xmm9 ,%%xmm15	\n\t"\
		"addpd	%%xmm0,%%xmm6						\n\t	addpd	%%xmm11,%%xmm12	\n\t"\
		"addpd	%%xmm2,%%xmm5						\n\t	movaps	%%xmm14,0x80(%%rsi)	\n\t"\
		"addpd	%%xmm1,%%xmm7						\n\t	movaps	%%xmm15,0x90(%%rsi)	\n\t"\
		"addpd	%%xmm3,%%xmm4						\n\t	movaps	0x200(%%rsi),%%xmm8	/* isrt2 */	\n\t"\
		"movaps	%%xmm6,     (%%rsi)					\n\t	movaps	%%xmm10,%%xmm14	\n\t"\
		"movaps	%%xmm5,0x060(%%rsi)					\n\t	movaps	%%xmm13,%%xmm15	\n\t"\
		"movaps	%%xmm7,0x010(%%rsi)					\n\t	subpd	%%xmm12,%%xmm10	\n\t"\
		"movaps	%%xmm4,0x030(%%rsi)					\n\t	subpd	%%xmm11,%%xmm13	\n\t"\
		"													addpd	%%xmm12,%%xmm14	\n\t"\
		"													addpd	%%xmm11,%%xmm15	\n\t"\
		"													mulpd	%%xmm8 ,%%xmm10	\n\t"\
		"													mulpd	%%xmm8 ,%%xmm13	\n\t"\
		"													mulpd	%%xmm8 ,%%xmm14	\n\t"\
		"						subpd	%%xmm10,%%xmm2	\n\t"\
		"													mulpd	%%xmm8 ,%%xmm15	\n\t"\
		"												/*	movaps	%%xmm10,0xa0(%%rsi)	*/\n\t"\
		"												/*	movaps	%%xmm13,0xe0(%%rsi)	*/\n\t"\
		"												/*	movaps	%%xmm14,0xb0(%%rsi)	*/\n\t"\
		"												/*	movaps	%%xmm15,0xf0(%%rsi)	*/\n\t"\
		"/* SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r20) */	\n\t"\
		"movaps	0x80(%%rsi),%%xmm11		\n\t"\
		"movaps	0x90(%%rsi),%%xmm12		\n\t"\
		"						mulpd	(%%r8),%%xmm10	\n\t"\
		"movaps	0xc0(%%rsi),%%xmm8		\n\t"\
		"subpd	%%xmm11,%%xmm6						\n\t/*	subpd	%%xmm10,%%xmm2	*/\n\t"\
		"subpd	%%xmm9 ,%%xmm0						\n\t	subpd	%%xmm15,%%xmm5	\n\t"\
		"						addpd	%%xmm2 ,%%xmm10	\n\t"\
		"						movaps	%%xmm2 ,0xa0(%%rsi)	\n\t"\
		"subpd	%%xmm12,%%xmm7						\n\t	subpd	%%xmm14,%%xmm4	\n\t"\
		"						movaps	(%%r8),%%xmm2	/* 2.0 */	\n\t"\
		"subpd	%%xmm8 ,%%xmm1						\n\t	subpd	%%xmm13,%%xmm3	\n\t"\
		"mulpd	%%xmm2 ,%%xmm11						\n\t/*	addpd	%%xmm10,%%xmm10	*/\n\t"\
		"mulpd	%%xmm2 ,%%xmm9						\n\t	mulpd	%%xmm2 ,%%xmm15	\n\t"\
		"mulpd	%%xmm2 ,%%xmm12						\n\t	mulpd	%%xmm2 ,%%xmm14	\n\t"\
		"mulpd	%%xmm2 ,%%xmm8						\n\t	mulpd	%%xmm2 ,%%xmm13	\n\t"\
		"addpd	%%xmm6 ,%%xmm11						\n\t/*	addpd	%%xmm2 ,%%xmm10	*/\n\t"\
		"addpd	%%xmm0 ,%%xmm9						\n\t	addpd	%%xmm5 ,%%xmm15	\n\t"\
		"addpd	%%xmm7 ,%%xmm12						\n\t	addpd	%%xmm4 ,%%xmm14	\n\t"\
		"addpd	%%xmm1 ,%%xmm8						\n\t	addpd	%%xmm3 ,%%xmm13	\n\t"\
		"movaps	%%xmm6 ,0x80(%%rsi)					\n\t/*	movaps	%%xmm2 ,0xa0(%%rsi)	*/\n\t"\
		"movaps	%%xmm0 ,0x40(%%rsi)					\n\t	movaps	%%xmm5 ,0x60(%%rsi)	\n\t"\
		"movaps	%%xmm7 ,0x90(%%rsi)					\n\t	movaps	%%xmm4 ,0xb0(%%rsi)	\n\t"\
		"movaps	%%xmm1 ,0xd0(%%rsi)					\n\t	movaps	%%xmm3 ,0xf0(%%rsi)	\n\t"\
		"movaps	%%xmm11,    (%%rsi)					\n\t	movaps	%%xmm10,0x20(%%rsi)	\n\t"\
		"movaps	%%xmm9 ,0xc0(%%rsi)					\n\t	movaps	%%xmm15,0xe0(%%rsi)	\n\t"\
		"movaps	%%xmm12,0x10(%%rsi)					\n\t	movaps	%%xmm14,0x30(%%rsi)	\n\t"\
		"movaps	%%xmm8 ,0x50(%%rsi)					\n\t	movaps	%%xmm13,0x70(%%rsi)	\n\t"\
		"\n\t"\
		"/***************************************/\n\t"\
		"\n\t"\
		"/*...Block 4: */	\n\t"\
		"movq	%[__add0],%%rax	\n\t"\
		"movslq	%[__p03],%%rdi	/* Do things this way [rather than repeatedly adding p1] since array-padding scheme means p2 == p+p1 not guaranteed. */	\n\t"\
		"movslq	%[__p08],%%rbx	/* Can't get these via simple load-one-and-shift-as-needed due to array padding scheme */	\n\t"\
		"movslq	%[__p10],%%rcx	\n\t"\
		"movslq	%[__p18],%%rdx	\n\t"\
		"shlq	$3,%%rdi	\n\t"\
		"shlq	$3,%%rbx	\n\t"\
		"shlq	$3,%%rcx							\n\t	movslq	%[__p04],%%r9	\n\t"\
		"shlq	$3,%%rdx							\n\t	shlq	$3,%%r9	\n\t"\
		"addq	%%rdi,%%rax	/* &a[j1+p3] */			\n\t	movq	%%r9,%%r10	\n\t	addq	%%rax,%%r10	\n\t"\
		"addq	%%rax,%%rbx							\n\t	movq	%%r9,%%r11	\n\t	addq	%%rbx,%%r11	\n\t"\
		"addq	%%rax,%%rcx							\n\t	movq	%%r9,%%r12	\n\t	addq	%%rcx,%%r12	\n\t"\
		"addq	%%rax,%%rdx							\n\t	movq	%%r9,%%r13	\n\t	addq	%%rdx,%%r13	\n\t"\
		"/* SSE2_RADIX4_DIF_4TWIDDLE_B(r30,c03) */	\n\t	/* SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO_B(r38,c07) */\n\t"\
		"addq	$0x100,%%rsi	/* r30 */			\n\t	movaps	    (%%r10),%%xmm8	\n\t"\
		"movaps	    (%%rax),%%xmm0					\n\t	movaps	    (%%r12),%%xmm12	\n\t"\
		"movaps	    (%%rcx),%%xmm4					\n\t	movaps	0x10(%%r10),%%xmm9	\n\t"\
		"movaps	0x10(%%rax),%%xmm1					\n\t	movaps	0x10(%%r12),%%xmm13	\n\t"\
		"movaps	0x10(%%rcx),%%xmm5					\n\t	movaps	0x4f0(%%rsi),%%xmm14	/* c07 */	\n\t"\
		"movaps	0x470(%%rsi),%%xmm6	/* c03 */		\n\t	movaps	0x500(%%rsi),%%xmm15	\n\t"\
		"movaps	0x480(%%rsi),%%xmm7					\n\t	movaps	%%xmm8 ,%%xmm10	\n\t"\
		"movaps	%%xmm0,%%xmm2						\n\t	movaps	%%xmm9 ,%%xmm11	\n\t"\
		"movaps	%%xmm1,%%xmm3						\n\t	mulpd	%%xmm14,%%xmm8	\n\t"\
		"mulpd   %%xmm6,%%xmm0						\n\t	mulpd	%%xmm14,%%xmm9	\n\t"\
		"mulpd   %%xmm6,%%xmm1						\n\t	mulpd	%%xmm15,%%xmm10	\n\t"\
		"mulpd   %%xmm7,%%xmm2						\n\t	mulpd	%%xmm15,%%xmm11	\n\t"\
		"mulpd   %%xmm7,%%xmm3						\n\t	movaps	%%xmm12,%%xmm14	\n\t"\
		"movaps	%%xmm4,%%xmm6						\n\t	addpd	%%xmm10,%%xmm9	\n\t"\
		"addpd   %%xmm2,%%xmm1						\n\t	movaps	%%xmm13,%%xmm15	\n\t"\
		"movaps	%%xmm5,%%xmm7						\n\t	mulpd	0x510(%%rsi),%%xmm12	/* c17 */	\n\t"\
		"mulpd   0x490(%%rsi),%%xmm4	/* c13 */	\n\t	subpd	%%xmm11,%%xmm8	\n\t"\
		"subpd   %%xmm3,%%xmm0						\n\t	mulpd	0x510(%%rsi),%%xmm13	\n\t"\
		"mulpd   0x490(%%rsi),%%xmm5				\n\t	mulpd	0x520(%%rsi),%%xmm14	\n\t"\
		"mulpd   0x4a0(%%rsi),%%xmm6				\n\t	movaps	%%xmm8 ,%%xmm10	\n\t"\
		"movaps	%%xmm0,%%xmm2						\n\t	mulpd	0x520(%%rsi),%%xmm15	\n\t"\
		"mulpd	0x4a0(%%rsi),%%xmm7					\n\t	addpd	%%xmm14,%%xmm13	\n\t"\
		"addpd	%%xmm6,%%xmm5						\n\t	movaps	%%xmm9 ,%%xmm11	\n\t"\
		"movaps	%%xmm1,%%xmm3						\n\t	subpd	%%xmm15,%%xmm12	\n\t"\
		"subpd	%%xmm7,%%xmm4						\n\t	addpd	%%xmm12,%%xmm8	\n\t"\
		"addpd	%%xmm4,%%xmm0						\n\t	addpd	%%xmm13,%%xmm9	\n\t"\
		"addpd	%%xmm5,%%xmm1						\n\t	subpd	%%xmm12,%%xmm10	\n\t"\
		"subpd	%%xmm4,%%xmm2						\n\t	subpd	%%xmm13,%%xmm11	\n\t"\
		"subpd	%%xmm5,%%xmm3						\n\t	movaps	    (%%r13),%%xmm12	\n\t"\
		"movaps	    (%%rdx),%%xmm4					\n\t	movaps	0x10(%%r13),%%xmm13	\n\t"\
		"movaps	0x10(%%rdx),%%xmm5					\n\t	movaps	    (%%r13),%%xmm14	\n\t"\
		"movaps	    (%%rdx),%%xmm6					\n\t	movaps	0x10(%%r13),%%xmm15	\n\t"\
		"movaps	0x10(%%rdx),%%xmm7					\n\t	mulpd	0x550(%%rsi),%%xmm12	/* c1F */	\n\t"\
		"mulpd	0x4d0(%%rsi),%%xmm4	/* c1B */		\n\t	mulpd	0x550(%%rsi),%%xmm13	\n\t"\
		"mulpd	0x4d0(%%rsi),%%xmm5					\n\t	mulpd	0x560(%%rsi),%%xmm14	\n\t"\
		"mulpd	0x4e0(%%rsi),%%xmm6					\n\t	mulpd	0x560(%%rsi),%%xmm15	\n\t"\
		"mulpd	0x4e0(%%rsi),%%xmm7					\n\t	addpd	%%xmm14,%%xmm13	\n\t"\
		"addpd	%%xmm6,%%xmm5						\n\t	subpd	%%xmm15,%%xmm12	\n\t"\
		"subpd	%%xmm7,%%xmm4						\n\t	movaps	%%xmm13,0x90(%%rsi)	\n\t"\
		"movaps	%%xmm5,0x010(%%rsi)					\n\t	movaps	%%xmm12,0x80(%%rsi)	/* tmpstr r08 */\n\t"\
		"movaps	%%xmm4,     (%%rsi)					\n\t	movaps	    (%%r11),%%xmm12	\n\t"\
		"movaps	    (%%rbx),%%xmm4					\n\t	movaps	0x10(%%r11),%%xmm13	\n\t"\
		"movaps	0x10(%%rbx),%%xmm5					\n\t	movaps	    (%%r11),%%xmm14	\n\t"\
		"movaps	    (%%rbx),%%xmm6					\n\t	movaps	0x10(%%r11),%%xmm15	\n\t"\
		"movaps	0x10(%%rbx),%%xmm7					\n\t	mulpd	0x530(%%rsi),%%xmm12	/* c0F */	\n\t"\
		"mulpd	0x4b0(%%rsi),%%xmm4	/* c0B */		\n\t	mulpd	0x530(%%rsi),%%xmm13	\n\t"\
		"mulpd	0x4b0(%%rsi),%%xmm5					\n\t	mulpd	0x540(%%rsi),%%xmm14	\n\t"\
		"mulpd	0x4c0(%%rsi),%%xmm6					\n\t	mulpd	0x540(%%rsi),%%xmm15	\n\t"\
		"mulpd	0x4c0(%%rsi),%%xmm7					\n\t	addpd	%%xmm14,%%xmm13	\n\t"\
		"addpd	%%xmm6,%%xmm5						\n\t	subpd	%%xmm15,%%xmm12	\n\t"\
		"subpd	%%xmm7,%%xmm4						\n\t	movaps	%%xmm13,%%xmm15	\n\t"\
		"movaps	%%xmm5,%%xmm7						\n\t	movaps	%%xmm12,%%xmm14	\n\t"\
		"movaps	%%xmm4,%%xmm6						\n\t	subpd	0x80(%%rsi),%%xmm12	\n\t"\
		"subpd	    (%%rsi),%%xmm4					\n\t	subpd	0x90(%%rsi),%%xmm13	\n\t"\
		"subpd	0x10(%%rsi),%%xmm5					\n\t	addpd	0x80(%%rsi),%%xmm14	\n\t"\
		"addpd	    (%%rsi),%%xmm6					\n\t	addpd	0x90(%%rsi),%%xmm15	\n\t"\
		"addpd	0x10(%%rsi),%%xmm7					\n\t	subpd	%%xmm14,%%xmm8	\n\t"\
		"subpd	%%xmm6,%%xmm0						\n\t	subpd	%%xmm15,%%xmm9	\n\t"\
		"subpd	%%xmm5,%%xmm2						\n\t	movaps	%%xmm8 ,0xc0(%%rsi)	\n\t"\
		"subpd	%%xmm7,%%xmm1						\n\t	subpd	%%xmm13,%%xmm10	\n\t"\
		"subpd	%%xmm4,%%xmm3						\n\t/*	movaps	%%xmm9 ,0xd0(%%rsi)	*/\n\t"\
		"movaps	%%xmm0,0x040(%%rsi)					\n\t	subpd	%%xmm12,%%xmm11	\n\t"\
		"movaps	%%xmm2,0x020(%%rsi)					\n\t	mulpd	(%%r8),%%xmm14	\n\t"\
		"movaps	%%xmm1,0x050(%%rsi)					\n\t	mulpd	(%%r8),%%xmm13	\n\t"\
		"movaps	%%xmm3,0x070(%%rsi)					\n\t	mulpd	(%%r8),%%xmm15	\n\t"\
		"mulpd	(%%r8),%%xmm6						\n\t	mulpd	(%%r8),%%xmm12	\n\t"\
		"mulpd	(%%r8),%%xmm5						\n\t	addpd	%%xmm8 ,%%xmm14	\n\t"\
		"mulpd	(%%r8),%%xmm7						\n\t	addpd	%%xmm10,%%xmm13	\n\t"\
		"mulpd	(%%r8),%%xmm4						\n\t	addpd	%%xmm9 ,%%xmm15	\n\t"\
		"addpd	%%xmm0,%%xmm6						\n\t	addpd	%%xmm11,%%xmm12	\n\t"\
		"addpd	%%xmm2,%%xmm5						\n\t	movaps	%%xmm14,0x80(%%rsi)	\n\t"\
		"addpd	%%xmm1,%%xmm7						\n\t	movaps	%%xmm15,0x90(%%rsi)	\n\t"\
		"addpd	%%xmm3,%%xmm4						\n\t	movaps	0x100(%%rsi),%%xmm8	/* isrt2 */	\n\t"\
		"movaps	%%xmm6,     (%%rsi)					\n\t	movaps	%%xmm10,%%xmm14	\n\t"\
		"movaps	%%xmm5,0x060(%%rsi)					\n\t	movaps	%%xmm13,%%xmm15	\n\t"\
		"movaps	%%xmm7,0x010(%%rsi)					\n\t	subpd	%%xmm12,%%xmm10	\n\t"\
		"movaps	%%xmm4,0x030(%%rsi)					\n\t	subpd	%%xmm11,%%xmm13	\n\t"\
		"													addpd	%%xmm12,%%xmm14	\n\t"\
		"													addpd	%%xmm11,%%xmm15	\n\t"\
		"													mulpd	%%xmm8 ,%%xmm10	\n\t"\
		"													mulpd	%%xmm8 ,%%xmm13	\n\t"\
		"													mulpd	%%xmm8 ,%%xmm14	\n\t"\
		"						subpd	%%xmm10,%%xmm2	\n\t"\
		"													mulpd	%%xmm8 ,%%xmm15	\n\t"\
		"												/*	movaps	%%xmm10,0xa0(%%rsi)	*/\n\t"\
		"												/*	movaps	%%xmm13,0xe0(%%rsi)	*/\n\t"\
		"												/*	movaps	%%xmm14,0xb0(%%rsi)	*/\n\t"\
		"												/*	movaps	%%xmm15,0xf0(%%rsi)	*/\n\t"\
		"/* SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r30) */	\n\t"\
		"movaps	0x80(%%rsi),%%xmm11		\n\t"\
		"movaps	0x90(%%rsi),%%xmm12		\n\t"\
		"						mulpd	(%%r8),%%xmm10	\n\t"\
		"movaps	0xc0(%%rsi),%%xmm8		\n\t"\
		"subpd	%%xmm11,%%xmm6						\n\t/*	subpd	%%xmm10,%%xmm2	*/\n\t"\
		"subpd	%%xmm9 ,%%xmm0						\n\t	subpd	%%xmm15,%%xmm5	\n\t"\
		"						addpd	%%xmm2 ,%%xmm10	\n\t"\
		"						movaps	%%xmm2 ,0xa0(%%rsi)	\n\t"\
		"subpd	%%xmm12,%%xmm7						\n\t	subpd	%%xmm14,%%xmm4	\n\t"\
		"						movaps	(%%r8),%%xmm2	/* 2.0 */	\n\t"\
		"subpd	%%xmm8 ,%%xmm1						\n\t	subpd	%%xmm13,%%xmm3	\n\t"\
		"mulpd	%%xmm2 ,%%xmm11						\n\t/*	addpd	%%xmm10,%%xmm10	*/\n\t"\
		"mulpd	%%xmm2 ,%%xmm9						\n\t	mulpd	%%xmm2 ,%%xmm15	\n\t"\
		"mulpd	%%xmm2 ,%%xmm12						\n\t	mulpd	%%xmm2 ,%%xmm14	\n\t"\
		"mulpd	%%xmm2 ,%%xmm8						\n\t	mulpd	%%xmm2 ,%%xmm13	\n\t"\
		"addpd	%%xmm6 ,%%xmm11						\n\t/*\n\t	addpd	%%xmm2 ,%%xmm10	*/\n\t"\
		"addpd	%%xmm0 ,%%xmm9						\n\t	addpd	%%xmm5 ,%%xmm15	\n\t"\
		"addpd	%%xmm7 ,%%xmm12						\n\t	addpd	%%xmm4 ,%%xmm14	\n\t"\
		"addpd	%%xmm1 ,%%xmm8						\n\t	addpd	%%xmm3 ,%%xmm13	\n\t"\
		"movaps	%%xmm6 ,0x80(%%rsi)					\n\t/*	movaps	%%xmm2 ,0xa0(%%rsi)	*/\n\t"\
		"movaps	%%xmm0 ,0x40(%%rsi)					\n\t	movaps	%%xmm5 ,0x60(%%rsi)	\n\t"\
		"movaps	%%xmm7 ,0x90(%%rsi)					\n\t	movaps	%%xmm4 ,0xb0(%%rsi)	\n\t"\
		"movaps	%%xmm1 ,0xd0(%%rsi)					\n\t	movaps	%%xmm3 ,0xf0(%%rsi)	\n\t"\
		"movaps	%%xmm11,    (%%rsi)					\n\t	movaps	%%xmm10,0x20(%%rsi)	\n\t"\
		"movaps	%%xmm9 ,0xc0(%%rsi)					\n\t	movaps	%%xmm15,0xe0(%%rsi)	\n\t"\
		"movaps	%%xmm12,0x10(%%rsi)					\n\t	movaps	%%xmm14,0x30(%%rsi)	\n\t"\
		"movaps	%%xmm8 ,0x50(%%rsi)					\n\t	movaps	%%xmm13,0x70(%%rsi)	\n\t"\
		"\n\t"\
		"/**********************************************************************************/\n\t"\
		"/*...and now do eight radix-4 transforms, including the internal twiddle factors: */\n\t"\
		"/**********************************************************************************/\n\t"\
		"/*...Block 1: t00,t10,t20,t30	*/			\n\t"\
		"movq	%[__add0],%%rax	/* &a[j1] */		\n\t"\
		"movslq	%[__p01],%%rbx						\n\t"\
		"movslq	%[__p02],%%rcx						\n\t"\
		"movslq	%[__p03],%%rdx						\n\t		/*...Block 5: t08,t18,t28,t38	*/	\n\t"\
		"shlq	$3,%%rbx							\n\t		movslq	%[__p04],%%r9	\n\t"\
		"shlq	$3,%%rcx							\n\t		shlq	$3,%%r9			\n\t"\
		"shlq	$3,%%rdx							\n\t		movq	%%r9,%%r10	\n\t	addq	%%rax,%%r10	\n\t"\
		"addq	%%rax,%%rbx							\n\t		movq	%%r9,%%r11	\n\t	addq	%%rbx,%%r11	\n\t"\
		"addq	%%rax,%%rcx							\n\t		movq	%%r9,%%r12	\n\t	addq	%%rcx,%%r12	\n\t"\
		"addq	%%rax,%%rdx							\n\t		movq	%%r9,%%r13	\n\t	addq	%%rdx,%%r13	\n\t"\
		"movq	%[__r00],%%rsi						\n\t		movaps	0x400(%%rsi),%%xmm11	/* isrt2 */	\n\t"\
		"movaps	     (%%rsi),%%xmm0					\n\t		movaps	0x280(%%rsi),%%xmm12	\n\t"\
		"movaps	0x200(%%rsi),%%xmm4					\n\t		movaps	0x290(%%rsi),%%xmm13	\n\t"\
		"movaps	0x010(%%rsi),%%xmm1					\n\t		movaps	0x380(%%rsi),%%xmm14	\n\t"\
		"movaps	0x210(%%rsi),%%xmm5					\n\t		movaps	0x390(%%rsi),%%xmm15	\n\t"\
		"movaps	0x100(%%rsi),%%xmm2					\n\t		mulpd	%%xmm11,%%xmm12	\n\t"\
		"movaps	0x300(%%rsi),%%xmm6					\n\t		movaps	0x080(%%rsi),%%xmm8	\n\t"\
		"movaps	0x110(%%rsi),%%xmm3					\n\t		mulpd	%%xmm11,%%xmm13	\n\t"\
		"movaps	0x310(%%rsi),%%xmm7					\n\t		movaps	0x090(%%rsi),%%xmm9	\n\t"\
		"subpd	0x100(%%rsi),%%xmm0					\n\t		mulpd	%%xmm11,%%xmm14	\n\t"\
		"subpd	0x300(%%rsi),%%xmm4					\n\t		movaps	0x180(%%rsi),%%xmm10	\n\t"\
		"subpd	0x110(%%rsi),%%xmm1					\n\t		mulpd	%%xmm11,%%xmm15	\n\t"\
		"subpd	0x310(%%rsi),%%xmm5					\n\t		movaps	0x190(%%rsi),%%xmm11	\n\t"\
		"addpd	     (%%rsi),%%xmm2					\n\t		subpd	%%xmm11,%%xmm8	\n\t"\
		"addpd	0x200(%%rsi),%%xmm6					\n\t		subpd	%%xmm13,%%xmm12	\n\t"\
		"addpd	0x010(%%rsi),%%xmm3					\n\t		subpd	%%xmm10,%%xmm9	\n\t"\
		"addpd	0x210(%%rsi),%%xmm7					\n\t		subpd	%%xmm14,%%xmm15	\n\t"\
		"subpd	%%xmm6,%%xmm2						\n\t		mulpd	(%%r8) ,%%xmm11	\n\t"\
		"subpd	%%xmm5,%%xmm0						\n\t		mulpd	(%%r8) ,%%xmm13	\n\t"\
		"subpd	%%xmm7,%%xmm3						\n\t		mulpd	(%%r8) ,%%xmm10	\n\t"\
		"subpd	%%xmm4,%%xmm1						\n\t		mulpd	(%%r8) ,%%xmm14	\n\t"\
		"movaps	%%xmm2,    (%%rbx)					\n\t		addpd	%%xmm8 ,%%xmm11	\n\t"\
		"movaps	(%%r8),%%xmm2	/* 2.0 */			\n\t		addpd	%%xmm12,%%xmm13	\n\t"\
		"mulpd	%%xmm2,%%xmm6						\n\t		addpd	%%xmm9 ,%%xmm10	\n\t"\
		"mulpd	%%xmm2,%%xmm5						\n\t		addpd	%%xmm15,%%xmm14	\n\t"\
		"mulpd	%%xmm2,%%xmm7						\n\t		subpd	%%xmm14,%%xmm12	\n\t"\
		"mulpd	%%xmm2,%%xmm4						\n\t		subpd	%%xmm15,%%xmm13	\n\t"\
		"/*movaps	%%xmm2,    (%%rbx)	*/			\n\t		mulpd	%%xmm2,%%xmm14	\n\t"\
		"movaps	%%xmm0,    (%%rcx)					\n\t		mulpd	%%xmm2,%%xmm15	\n\t"\
		"movaps	%%xmm3,0x10(%%rbx)					\n\t		addpd	%%xmm12,%%xmm14	\n\t"\
		"movaps	%%xmm1,0x10(%%rdx)					\n\t		addpd	%%xmm13,%%xmm15	\n\t"\
		"addpd	(%%rbx),%%xmm6						\n\t		subpd	%%xmm12,%%xmm8	\n\t"\
		"addpd	%%xmm0,%%xmm5						\n\t		subpd	%%xmm15,%%xmm11	\n\t"\
		"addpd	%%xmm3,%%xmm7						\n\t		subpd	%%xmm13,%%xmm10	\n\t"\
		"addpd	%%xmm1,%%xmm4						\n\t		subpd	%%xmm14,%%xmm9	\n\t"\
		"movaps	%%xmm6,    (%%rax)					\n\t		mulpd	%%xmm2,%%xmm12	\n\t"\
		"movaps	%%xmm5,    (%%rdx)					\n\t		mulpd	%%xmm2,%%xmm15	\n\t"\
		"movaps	%%xmm7,0x10(%%rax)					\n\t		mulpd	%%xmm2,%%xmm13	\n\t"\
		"movaps	%%xmm4,0x10(%%rcx)					\n\t		mulpd	%%xmm2,%%xmm14	\n\t"\
		"											\n\t		movaps	%%xmm8 ,    (%%r11)	\n\t"\
		"											\n\t		movaps	%%xmm11,    (%%r12)	\n\t"\
		"											\n\t		movaps	%%xmm10,0x10(%%r11)	\n\t"\
		"											\n\t		movaps	%%xmm9 ,0x10(%%r13)	\n\t"\
		"											\n\t		addpd	%%xmm8 ,%%xmm12		\n\t"\
		"											\n\t		addpd	%%xmm11,%%xmm15		\n\t"\
		"											\n\t		addpd	%%xmm10,%%xmm13		\n\t"\
		"											\n\t		addpd	%%xmm9 ,%%xmm14		\n\t"\
		"											\n\t		movaps	%%xmm12,    (%%r10)	\n\t"\
		"											\n\t		movaps	%%xmm15,    (%%r13)	\n\t"\
		"											\n\t		movaps	%%xmm13,0x10(%%r10)	\n\t"\
		"											\n\t		movaps	%%xmm14,0x10(%%r12)	\n\t"\
		"\n\t"\
		"/*...Block 3: t04,t14,t24,t34	*/			\n\t"\
		"addq	$0x40,%%rsi	/* r04 */				\n\t"\
		"subq	%%rax,%%rbx	/* p01 << 3 */			\n\t"\
		"subq	%%rax,%%rcx	/* p02 << 3 */			\n\t"\
		"subq	%%rax,%%rdx	/* p03 << 3 */			\n\t		/*...Block 7: t0C,t1C,t2C,t3C	*/\n\t"\
		"movslq	%[__p08],%%rdi						\n\t		movslq	%[__p04],%%r9	\n\t"\
		"shlq	$3,%%rdi							\n\t		shlq	$3,%%r9			\n\t"\
		"addq	%%rdi,%%rax	/* &a[j1+p08] */		\n\t		movq	%%r9,%%r10	\n\t	addq	%%rax,%%r10	\n\t"\
		"addq	%%rax,%%rbx							\n\t		movq	%%r9,%%r11	\n\t	addq	%%rbx,%%r11	\n\t"\
		"addq	%%rax,%%rcx							\n\t		movq	%%r9,%%r12	\n\t	addq	%%rcx,%%r12	\n\t"\
		"addq	%%rax,%%rdx							\n\t		movq	%%r9,%%r13	\n\t	addq	%%rdx,%%r13	\n\t"\
		"movaps	0x200(%%rsi),%%xmm4					\n\t		movaps	0x280(%%rsi),%%xmm12	\n\t"\
		"movaps	0x210(%%rsi),%%xmm5					\n\t		movaps	0x290(%%rsi),%%xmm13	\n\t"\
		"movaps	0x3d0(%%rsi),%%xmm3	/* cc0 */		\n\t"\
		"movaps	0x3e0(%%rsi),%%xmm2					\n\t"\
		"movaps	%%xmm4,%%xmm6						\n\t		movaps	%%xmm12,%%xmm14	\n\t"\
		"movaps	%%xmm5,%%xmm7						\n\t		movaps	%%xmm13,%%xmm15	\n\t"\
		"mulpd	%%xmm3,%%xmm4						\n\t		mulpd	%%xmm2 ,%%xmm12	\n\t"\
		"mulpd	%%xmm3,%%xmm5						\n\t		mulpd	%%xmm2 ,%%xmm13	\n\t"\
		"mulpd	%%xmm2,%%xmm6						\n\t		mulpd	%%xmm3 ,%%xmm14	\n\t"\
		"movaps	0x300(%%rsi),%%xmm0					\n\t		movaps	0x380(%%rsi),%%xmm8	\n\t"\
		"mulpd	%%xmm2,%%xmm7						\n\t		mulpd	%%xmm3 ,%%xmm15	\n\t"\
		"movaps	0x310(%%rsi),%%xmm1					\n\t		movaps	0x390(%%rsi),%%xmm9	\n\t"\
		"addpd	%%xmm6,%%xmm5						\n\t		addpd	%%xmm14,%%xmm13	\n\t"\
		"movaps	%%xmm0,%%xmm6						\n\t		movaps	%%xmm8 ,%%xmm14	\n\t"\
		"subpd	%%xmm7,%%xmm4						\n\t		subpd	%%xmm15,%%xmm12	\n\t"\
		"movaps	%%xmm1,%%xmm7						\n\t		movaps	%%xmm9 ,%%xmm15	\n\t"\
		"mulpd	%%xmm2,%%xmm6						\n\t		mulpd	%%xmm3 ,%%xmm14	\n\t"\
		"mulpd	%%xmm2,%%xmm7						\n\t		mulpd	%%xmm3 ,%%xmm15	\n\t"\
		"mulpd	%%xmm3,%%xmm0						\n\t		mulpd	%%xmm2 ,%%xmm8	\n\t"\
		"mulpd	%%xmm3,%%xmm1						\n\t		mulpd	%%xmm2 ,%%xmm9	\n\t"\
		"addpd	%%xmm0,%%xmm7						\n\t		addpd	%%xmm8 ,%%xmm15	\n\t"\
		"subpd	%%xmm1,%%xmm6						\n\t		subpd	%%xmm9 ,%%xmm14	\n\t"\
		"movaps	%%xmm4,%%xmm2						\n\t		movaps	%%xmm12,%%xmm10	\n\t"\
		"movaps	%%xmm5,%%xmm3						\n\t		movaps	%%xmm13,%%xmm11	\n\t"\
		"subpd	%%xmm6,%%xmm4						\n\t		subpd	%%xmm14,%%xmm12	\n\t"\
		"subpd	%%xmm7,%%xmm5						\n\t		subpd	%%xmm15,%%xmm13	\n\t"\
		"addpd	%%xmm2,%%xmm6						\n\t		addpd	%%xmm10,%%xmm14	\n\t"\
		"addpd	%%xmm3,%%xmm7						\n\t		addpd	%%xmm11,%%xmm15	\n\t"\
		"movaps	0x100(%%rsi),%%xmm2					\n\t		movaps	0x180(%%rsi),%%xmm10	\n\t"\
		"movaps	0x110(%%rsi),%%xmm3					\n\t		movaps	0x190(%%rsi),%%xmm11	\n\t"\
		"movaps	0x3c0(%%rsi),%%xmm1	/* isrt2 */		\n\t"\
		"movaps	%%xmm2,%%xmm0						\n\t		movaps	%%xmm10,%%xmm8	\n\t"\
		"subpd	%%xmm3,%%xmm2						\n\t		addpd	%%xmm11,%%xmm10	\n\t"\
		"addpd	%%xmm0,%%xmm3						\n\t		subpd	%%xmm8 ,%%xmm11	\n\t"\
		"mulpd	%%xmm1,%%xmm2						\n\t		mulpd	%%xmm1 ,%%xmm10	\n\t"\
		"mulpd	%%xmm1,%%xmm3						\n\t		mulpd	%%xmm1 ,%%xmm11	\n\t"\
		"movaps	     (%%rsi),%%xmm0					\n\t		movaps	0x080(%%rsi),%%xmm8	\n\t"\
		"movaps	0x010(%%rsi),%%xmm1					\n\t		movaps	0x090(%%rsi),%%xmm9	\n\t"\
		"subpd	%%xmm2,%%xmm0						\n\t		subpd	%%xmm10,%%xmm8	\n\t"\
		"subpd	%%xmm3,%%xmm1						\n\t		subpd	%%xmm11,%%xmm9	\n\t"\
		"mulpd	(%%r8),%%xmm2						\n\t		mulpd	(%%r8),%%xmm10	\n\t"\
		"mulpd	(%%r8),%%xmm3						\n\t		mulpd	(%%r8),%%xmm11	\n\t"\
		"addpd	%%xmm0,%%xmm2						\n\t		addpd	%%xmm8 ,%%xmm10	\n\t"\
		"addpd	%%xmm1,%%xmm3						\n\t		addpd	%%xmm9 ,%%xmm11	\n\t"\
		"subpd	%%xmm6,%%xmm2						\n\t		subpd	%%xmm12,%%xmm8	\n\t"\
		"subpd	%%xmm5,%%xmm0						\n\t		subpd	%%xmm15,%%xmm10	\n\t"\
		"subpd	%%xmm7,%%xmm3						\n\t		subpd	%%xmm13,%%xmm9	\n\t"\
		"subpd	%%xmm4,%%xmm1						\n\t		subpd	%%xmm14,%%xmm11	\n\t"\
		"mulpd	(%%r8),%%xmm6						\n\t		mulpd	(%%r8),%%xmm12	\n\t"\
		"mulpd	(%%r8),%%xmm5						\n\t		mulpd	(%%r8),%%xmm15	\n\t"\
		"mulpd	(%%r8),%%xmm7						\n\t		mulpd	(%%r8),%%xmm13	\n\t"\
		"mulpd	(%%r8),%%xmm4						\n\t		mulpd	(%%r8),%%xmm14	\n\t"\
		"movaps	%%xmm2,    (%%rbx)					\n\t		movaps	%%xmm8 ,    (%%r11)	\n\t"\
		"movaps	%%xmm0,    (%%rcx)					\n\t		movaps	%%xmm10,    (%%r12)	\n\t"\
		"movaps	%%xmm3,0x10(%%rbx)					\n\t		movaps	%%xmm9 ,0x10(%%r11)	\n\t"\
		"movaps	%%xmm1,0x10(%%rdx)					\n\t		movaps	%%xmm11,0x10(%%r13)	\n\t"\
		"addpd	%%xmm2,%%xmm6						\n\t		addpd	%%xmm8 ,%%xmm12		\n\t"\
		"addpd	%%xmm0,%%xmm5						\n\t		addpd	%%xmm10,%%xmm15		\n\t"\
		"addpd	%%xmm3,%%xmm7						\n\t		addpd	%%xmm9 ,%%xmm13		\n\t"\
		"addpd	%%xmm1,%%xmm4						\n\t		addpd	%%xmm11,%%xmm14		\n\t"\
		"movaps	%%xmm6,    (%%rax)					\n\t		movaps	%%xmm12,    (%%r10)	\n\t"\
		"movaps	%%xmm5,    (%%rdx)					\n\t		movaps	%%xmm15,    (%%r13)	\n\t"\
		"movaps	%%xmm7,0x10(%%rax)					\n\t		movaps	%%xmm13,0x10(%%r10)	\n\t"\
		"movaps	%%xmm4,0x10(%%rcx)					\n\t		movaps	%%xmm14,0x10(%%r12)	\n\t"\
		"\n\t"\
		"/*...Block 2: t02,t12,t22,t32	*/			\n\t"\
		"subq	$0x20,%%rsi	/* r02 */				\n\t"\
		"movslq	%[__p10],%%rdi						\n\t"\
		"subq	%%rax,%%rbx							\n\t"\
		"subq	%%rax,%%rcx							\n\t"\
		"subq	%%rax,%%rdx							\n\t		/*...Block 6: t0A,t1A,t2A,t3A	*/	\n\t"\
		"movq	%[__add0],%%rax	/* &a[j1] */		\n\t		movslq	%[__p04],%%r9	\n\t"\
		"shlq	$3,%%rdi							\n\t		shlq	$3,%%r9			\n\t"\
		"addq	%%rdi,%%rax	/* &a[j1+p10) */		\n\t		movq	%%r9,%%r10	\n\t	addq	%%rax,%%r10	\n\t"\
		"addq	%%rax,%%rbx							\n\t		movq	%%r9,%%r11	\n\t	addq	%%rbx,%%r11	\n\t"\
		"addq	%%rax,%%rcx							\n\t		movq	%%r9,%%r12	\n\t	addq	%%rcx,%%r12	\n\t"\
		"addq	%%rax,%%rdx							\n\t		movq	%%r9,%%r13	\n\t	addq	%%rdx,%%r13	\n\t"\
		"movaps	0x200(%%rsi),%%xmm4					\n\t		movaps	0x280(%%rsi),%%xmm12	\n\t"\
		"movaps	0x210(%%rsi),%%xmm5					\n\t		movaps	0x290(%%rsi),%%xmm13	\n\t"\
		"movaps	0x410(%%rsi),%%xmm2	/* cc1 */		\n\t		movaps	0x430(%%rsi),%%xmm11	/* cc3 */	\n\t"\
		"movaps	0x420(%%rsi),%%xmm3					\n\t		movaps	0x440(%%rsi),%%xmm10	\n\t"\
		"movaps	%%xmm4,%%xmm6						\n\t		movaps	%%xmm12,%%xmm14	\n\t"\
		"movaps	%%xmm5,%%xmm7						\n\t		movaps	%%xmm13,%%xmm15	\n\t"\
		"mulpd	%%xmm2,%%xmm4						\n\t		mulpd	%%xmm10,%%xmm12	\n\t"\
		"mulpd	%%xmm2,%%xmm5						\n\t		mulpd	%%xmm10,%%xmm13	\n\t"\
		"mulpd	%%xmm3,%%xmm6						\n\t		mulpd	%%xmm11,%%xmm14	\n\t"\
		"movaps	0x300(%%rsi),%%xmm0					\n\t		movaps	0x380(%%rsi),%%xmm8	\n\t"\
		"mulpd	%%xmm3,%%xmm7						\n\t		mulpd	%%xmm11,%%xmm15	\n\t"\
		"movaps	0x310(%%rsi),%%xmm1					\n\t		movaps	0x390(%%rsi),%%xmm9	\n\t"\
		"addpd	%%xmm6,%%xmm5						\n\t		addpd	%%xmm14,%%xmm13	\n\t"\
		"movaps	%%xmm0,%%xmm6						\n\t		movaps	%%xmm8 ,%%xmm14	\n\t"\
		"subpd	%%xmm7,%%xmm4						\n\t		subpd	%%xmm15,%%xmm12	\n\t"\
		"movaps	%%xmm1,%%xmm7						\n\t		movaps	%%xmm9 ,%%xmm15	\n\t"\
		"mulpd	%%xmm11,%%xmm6						\n\t		mulpd	%%xmm2 ,%%xmm14	\n\t"\
		"mulpd	%%xmm11,%%xmm7						\n\t		mulpd	%%xmm2 ,%%xmm15	\n\t"\
		"mulpd	%%xmm10,%%xmm0						\n\t		mulpd	%%xmm3 ,%%xmm8	\n\t"\
		"mulpd	%%xmm10,%%xmm1						\n\t		mulpd	%%xmm3 ,%%xmm9	\n\t"\
		"addpd	%%xmm0,%%xmm7						\n\t		subpd	%%xmm8 ,%%xmm15	\n\t"\
		"subpd	%%xmm1,%%xmm6						\n\t		addpd	%%xmm9 ,%%xmm14	\n\t"\
		"movaps	%%xmm4,%%xmm2						\n\t		movaps	%%xmm12,%%xmm10	\n\t"\
		"movaps	%%xmm5,%%xmm3						\n\t		movaps	%%xmm13,%%xmm11	\n\t"\
		"subpd	%%xmm6,%%xmm4						\n\t		subpd	%%xmm14,%%xmm12	\n\t"\
		"subpd	%%xmm7,%%xmm5						\n\t		subpd	%%xmm15,%%xmm13	\n\t"\
		"addpd	%%xmm2,%%xmm6						\n\t		addpd	%%xmm10,%%xmm14	\n\t"\
		"addpd	%%xmm3,%%xmm7						\n\t		addpd	%%xmm11,%%xmm15	\n\t"\
		"movaps	0x100(%%rsi),%%xmm1					\n\t		movaps	0x180(%%rsi),%%xmm9	\n\t"\
		"movaps	0x110(%%rsi),%%xmm3					\n\t		movaps	0x190(%%rsi),%%xmm11	\n\t"\
		"movaps	0x400(%%rsi),%%xmm0	/* ss0 */		\n\t		movaps	0x3f0(%%rsi),%%xmm8		/* cc0 */	\n\t"\
		"movaps	%%xmm1,%%xmm2						\n\t		movaps	%%xmm9 ,%%xmm10	\n\t"\
		"mulpd	%%xmm0,%%xmm1						\n\t		mulpd	%%xmm8 ,%%xmm9	\n\t"\
		"mulpd	%%xmm3,%%xmm0						\n\t		mulpd	%%xmm11,%%xmm8	\n\t"\
		"mulpd	0x3f0(%%rsi),%%xmm2		/* cc0 */	\n\t		mulpd	0x400(%%rsi),%%xmm10	/* ss0 */	\n\t"\
		"mulpd	0x3f0(%%rsi),%%xmm3					\n\t		mulpd	0x400(%%rsi),%%xmm11	\n\t"\
		"subpd	%%xmm0,%%xmm2						\n\t		addpd	%%xmm8,%%xmm10	\n\t"\
		"addpd	%%xmm1,%%xmm3						\n\t		subpd	%%xmm9,%%xmm11	\n\t"\
		"movaps	     (%%rsi),%%xmm0					\n\t		movaps	0x080(%%rsi),%%xmm8	\n\t"\
		"movaps	0x010(%%rsi),%%xmm1					\n\t		movaps	0x090(%%rsi),%%xmm9	\n\t"\
		"subpd	%%xmm2,%%xmm0						\n\t		subpd	%%xmm10,%%xmm8	\n\t"\
		"subpd	%%xmm3,%%xmm1						\n\t		subpd	%%xmm11,%%xmm9	\n\t"\
		"mulpd	(%%r8),%%xmm2						\n\t		mulpd	(%%r8),%%xmm10	\n\t"\
		"mulpd	(%%r8),%%xmm3						\n\t		mulpd	(%%r8),%%xmm11	\n\t"\
		"addpd	%%xmm0,%%xmm2						\n\t		addpd	%%xmm8 ,%%xmm10	\n\t"\
		"addpd	%%xmm1,%%xmm3						\n\t		addpd	%%xmm9 ,%%xmm11	\n\t"\
		"subpd	%%xmm6,%%xmm2						\n\t		subpd	%%xmm12,%%xmm8	\n\t"\
		"subpd	%%xmm5,%%xmm0						\n\t		subpd	%%xmm15,%%xmm10	\n\t"\
		"subpd	%%xmm7,%%xmm3						\n\t		subpd	%%xmm13,%%xmm9	\n\t"\
		"subpd	%%xmm4,%%xmm1						\n\t		subpd	%%xmm14,%%xmm11	\n\t"\
		"mulpd	(%%r8),%%xmm6						\n\t		mulpd	(%%r8),%%xmm12	\n\t"\
		"mulpd	(%%r8),%%xmm5						\n\t		mulpd	(%%r8),%%xmm15	\n\t"\
		"mulpd	(%%r8),%%xmm7						\n\t		mulpd	(%%r8),%%xmm13	\n\t"\
		"mulpd	(%%r8),%%xmm4						\n\t		mulpd	(%%r8),%%xmm14	\n\t"\
		"movaps	%%xmm2,    (%%rbx)					\n\t		movaps	%%xmm8 ,    (%%r11)	\n\t"\
		"movaps	%%xmm0,    (%%rcx)					\n\t		movaps	%%xmm10,    (%%r12)	\n\t"\
		"movaps	%%xmm3,0x10(%%rbx)					\n\t		movaps	%%xmm9 ,0x10(%%r11)	\n\t"\
		"movaps	%%xmm1,0x10(%%rdx)					\n\t		movaps	%%xmm11,0x10(%%r13)	\n\t"\
		"addpd	%%xmm2,%%xmm6						\n\t		addpd	%%xmm8 ,	%%xmm12	\n\t"\
		"addpd	%%xmm0,%%xmm5						\n\t		addpd	%%xmm10,	%%xmm15	\n\t"\
		"addpd	%%xmm3,%%xmm7						\n\t		addpd	%%xmm9 ,	%%xmm13	\n\t"\
		"addpd	%%xmm1,%%xmm4						\n\t		addpd	%%xmm11,	%%xmm14	\n\t"\
		"movaps	%%xmm6,    (%%rax)					\n\t		movaps	%%xmm12,    (%%r10)	\n\t"\
		"movaps	%%xmm5,    (%%rdx)					\n\t		movaps	%%xmm15,    (%%r13)	\n\t"\
		"movaps	%%xmm7,0x10(%%rax)					\n\t		movaps	%%xmm13,0x10(%%r10)	\n\t"\
		"movaps	%%xmm4,0x10(%%rcx)					\n\t		movaps	%%xmm14,0x10(%%r12)	\n\t"\
		"\n\t"\
		"/*...Block 4: t06,t16,t26,t36	*/	\n\t"\
		"addq	$0x40,%%rsi	/* r06 */	\n\t"\
		"movslq	%[__p18],%%rdi	\n\t"\
		"subq	%%rax,%%rbx	\n\t"\
		"subq	%%rax,%%rcx	\n\t"\
		"subq	%%rax,%%rdx							\n\t		/*...Block 8: t0E,t1E,t2E,t3E	*/	\n\t"\
		"movq	%[__add0],%%rax	/* &a[j1] */		\n\t		movslq	%[__p04],%%r9	\n\t"\
		"shlq	$3,%%rdi							\n\t		shlq	$3,%%r9			\n\t"\
		"addq	%%rdi,%%rax	/* &a[j1+p18] */		\n\t		movq	%%r9,%%r10	\n\t	addq	%%rax,%%r10	\n\t"\
		"addq	%%rax,%%rbx							\n\t		movq	%%r9,%%r11	\n\t	addq	%%rbx,%%r11	\n\t"\
		"addq	%%rax,%%rcx							\n\t		movq	%%r9,%%r12	\n\t	addq	%%rcx,%%r12	\n\t"\
		"addq	%%rax,%%rdx							\n\t		movq	%%r9,%%r13	\n\t	addq	%%rdx,%%r13	\n\t"\
		"movaps	0x200(%%rsi),%%xmm4					\n\t		movaps	0x280(%%rsi),%%xmm12	\n\t"\
		"movaps	0x210(%%rsi),%%xmm5					\n\t		movaps	0x290(%%rsi),%%xmm13	\n\t"\
		"movaps	0x3f0(%%rsi),%%xmm2	/* cc3 */		\n\t		movaps	0x3d0(%%rsi),%%xmm11 /* cc1 */	\n\t"\
		"movaps	0x400(%%rsi),%%xmm3					\n\t		movaps	0x3e0(%%rsi),%%xmm10 \n\t"\
		"movaps	%%xmm4,%%xmm6						\n\t		movaps	%%xmm12,%%xmm14	\n\t"\
		"movaps	%%xmm5,%%xmm7						\n\t		movaps	%%xmm13,%%xmm15	\n\t"\
		"mulpd	%%xmm2,%%xmm4						\n\t		mulpd	%%xmm10,%%xmm12	\n\t"\
		"mulpd	%%xmm2,%%xmm5						\n\t		mulpd	%%xmm10,%%xmm13	\n\t"\
		"mulpd	%%xmm3,%%xmm6						\n\t		mulpd	%%xmm11,%%xmm14	\n\t"\
		"movaps	0x300(%%rsi),%%xmm0					\n\t		movaps	0x380(%%rsi),%%xmm8	\n\t"\
		"mulpd	%%xmm3,%%xmm7						\n\t		mulpd	%%xmm11,%%xmm15	\n\t"\
		"movaps	0x310(%%rsi),%%xmm1					\n\t		movaps	0x390(%%rsi),%%xmm9	\n\t"\
		"addpd	%%xmm6,%%xmm5						\n\t		addpd	%%xmm14,%%xmm13	\n\t"\
		"movaps	%%xmm0,%%xmm6						\n\t		movaps	%%xmm8,%%xmm14	\n\t"\
		"subpd	%%xmm7,%%xmm4						\n\t		subpd	%%xmm15,%%xmm12	\n\t"\
		"movaps	%%xmm1,%%xmm7						\n\t		movaps	%%xmm9,%%xmm15	\n\t"\
		"mulpd	%%xmm10,%%xmm6						\n\t		mulpd	%%xmm3,%%xmm14	\n\t"\
		"mulpd	%%xmm10,%%xmm7						\n\t		mulpd	%%xmm3,%%xmm15	\n\t"\
		"mulpd	%%xmm11,%%xmm0						\n\t		mulpd	%%xmm2,%%xmm8	\n\t"\
		"mulpd	%%xmm11,%%xmm1						\n\t		mulpd	%%xmm2,%%xmm9	\n\t"\
		"subpd	%%xmm0,%%xmm7						\n\t		addpd	%%xmm8,%%xmm15	\n\t"\
		"addpd	%%xmm1,%%xmm6						\n\t		subpd	%%xmm9,%%xmm14	\n\t"\
		"movaps	%%xmm4,%%xmm2						\n\t		movaps	%%xmm12,%%xmm10	\n\t"\
		"movaps	%%xmm5,%%xmm3						\n\t		movaps	%%xmm13,%%xmm11	\n\t"\
		"subpd	%%xmm6,%%xmm4						\n\t		subpd	%%xmm14,%%xmm12	\n\t"\
		"subpd	%%xmm7,%%xmm5						\n\t		subpd	%%xmm15,%%xmm13	\n\t"\
		"addpd	%%xmm2,%%xmm6						\n\t		addpd	%%xmm10,%%xmm14	\n\t"\
		"addpd	%%xmm3,%%xmm7						\n\t		addpd	%%xmm11,%%xmm15	\n\t"\
		"movaps	0x100(%%rsi),%%xmm1					\n\t		movaps	0x180(%%rsi),%%xmm9	\n\t"\
		"movaps	0x110(%%rsi),%%xmm3					\n\t		movaps	0x190(%%rsi),%%xmm11	\n\t"\
		"movaps	0x3b0(%%rsi),%%xmm0		/* cc0 */	\n\t		movaps	0x3c0(%%rsi),%%xmm8	/* ss0 */	\n\t"\
		"movaps	%%xmm1,%%xmm2						\n\t		movaps	%%xmm9,%%xmm10	\n\t"\
		"mulpd	%%xmm0,%%xmm1						\n\t		mulpd	%%xmm8,%%xmm9	\n\t"\
		"mulpd	%%xmm3,%%xmm0						\n\t		mulpd	%%xmm11,%%xmm8	\n\t"\
		"mulpd	0x3c0(%%rsi),%%xmm2		/* ss0 */	\n\t		mulpd	0x3b0(%%rsi),%%xmm10	/* cc0 */	\n\t"\
		"mulpd	0x3c0(%%rsi),%%xmm3					\n\t		mulpd	0x3b0(%%rsi),%%xmm11	\n\t"\
		"subpd	%%xmm0,%%xmm2						\n\t		addpd	%%xmm8,%%xmm10	\n\t"\
		"addpd	%%xmm1,%%xmm3						\n\t		subpd	%%xmm9,%%xmm11	\n\t"\
		"movaps	     (%%rsi),%%xmm0					\n\t		movaps	0x080(%%rsi),%%xmm8	\n\t"\
		"movaps	0x010(%%rsi),%%xmm1					\n\t		movaps	0x090(%%rsi),%%xmm9	\n\t"\
		"subpd	%%xmm2,%%xmm0						\n\t		subpd	%%xmm10,%%xmm8	\n\t"\
		"subpd	%%xmm3,%%xmm1						\n\t		subpd	%%xmm11,%%xmm9	\n\t"\
		"mulpd	(%%r8),%%xmm2						\n\t		mulpd	(%%r8),%%xmm10	\n\t"\
		"mulpd	(%%r8),%%xmm3						\n\t		mulpd	(%%r8),%%xmm11	\n\t"\
		"addpd	%%xmm0,%%xmm2						\n\t		addpd	%%xmm8,%%xmm10	\n\t"\
		"addpd	%%xmm1,%%xmm3						\n\t		addpd	%%xmm9,%%xmm11	\n\t"\
		"subpd	%%xmm4,%%xmm2						\n\t		subpd	%%xmm12,%%xmm8	\n\t"\
		"subpd	%%xmm7,%%xmm0						\n\t		subpd	%%xmm15,%%xmm10	\n\t"\
		"subpd	%%xmm5,%%xmm3						\n\t		subpd	%%xmm13,%%xmm9	\n\t"\
		"subpd	%%xmm6,%%xmm1						\n\t		subpd	%%xmm14,%%xmm11	\n\t"\
		"mulpd	(%%r8),%%xmm4						\n\t		mulpd	(%%r8),%%xmm12	\n\t"\
		"mulpd	(%%r8),%%xmm7						\n\t		mulpd	(%%r8),%%xmm15	\n\t"\
		"mulpd	(%%r8),%%xmm5						\n\t		mulpd	(%%r8),%%xmm13	\n\t"\
		"mulpd	(%%r8),%%xmm6						\n\t		mulpd	(%%r8),%%xmm14	\n\t"\
		"movaps	%%xmm2,    (%%rbx)					\n\t		movaps	%%xmm8 ,    (%%r11)	\n\t"\
		"movaps	%%xmm0,    (%%rcx)					\n\t		movaps	%%xmm10,    (%%r12)	\n\t"\
		"movaps	%%xmm3,0x10(%%rbx)					\n\t		movaps	%%xmm9 ,0x10(%%r11)	\n\t"\
		"movaps	%%xmm1,0x10(%%rdx)					\n\t		movaps	%%xmm11,0x10(%%r13)	\n\t"\
		"addpd	%%xmm2,	%%xmm4						\n\t		addpd	%%xmm8 ,	%%xmm12	\n\t"\
		"addpd	%%xmm0,	%%xmm7						\n\t		addpd	%%xmm10,	%%xmm15	\n\t"\
		"addpd	%%xmm3,	%%xmm5						\n\t		addpd	%%xmm9 ,	%%xmm13	\n\t"\
		"addpd	%%xmm1,	%%xmm6						\n\t		addpd	%%xmm11,	%%xmm14	\n\t"\
		"movaps	%%xmm4,    (%%rax)					\n\t		movaps	%%xmm12,    (%%r10)	\n\t"\
		"movaps	%%xmm7,    (%%rdx)					\n\t		movaps	%%xmm15,    (%%r13)	\n\t"\
		"movaps	%%xmm5,0x10(%%rax)					\n\t		movaps	%%xmm13,0x10(%%r10)	\n\t"\
		"movaps	%%xmm6,0x10(%%rcx)					\n\t		movaps	%%xmm14,0x10(%%r12)	\n\t"\
		"\n\t"\
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
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r8","r9","r10","r11","r12","r13","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"		/* Clobbered registers */\
	);\
	}

	/*
	For GCC-macro version of this, use that isrt2 + 0x010,0x030,0x050 = cc0,cc1,cc3,
	and isrt2 + 0x070,0x0f0,0x170,0x1f0,0x270,0x2f0,0x370,0x3f0 = c00,04,02,06,01,05,03,07
	in order to reduce number of args to <= the GCC-allowed maximum of 30:
	*/
	#define SSE2_RADIX32_DIT_TWIDDLE(Xadd0,Xp01,Xp02,Xp03,Xp04,Xp05,Xp06,Xp07,Xp08,Xp10,Xp18,Xr00,Xr10,Xr20,Xr30,Xisrt2)\
	{\
	__asm__ volatile (\
	"/**************************...Block 1: ****************************/\n\t"\
		"movq	%[__r00],%%rsi			\n\t"\
		"movq	%[__add0],%%rax			\n\t"\
		"movslq	%[__p01],%%rbx			\n\t"\
		"movslq	%[__p02],%%rcx			\n\t"\
		"movslq	%[__p03],%%rdx			\n\t		movq	%[__isrt2],%%r8	\n\t"\
		"shlq	$3,%%rbx				\n\t		movslq	%[__p04],%%r9	\n\t"\
		"shlq	$3,%%rcx				\n\t		shlq	$3,%%r9		\n\t	addq	$0x480,%%r8	/* two */\n\t"\
		"shlq	$3,%%rdx				\n\t		movq	%%r9,%%r10	\n\t	addq	%%rax,%%r10	\n\t"\
		"addq	%%rax,%%rbx				\n\t		movq	%%r9,%%r11	\n\t	addq	%%rbx,%%r11	\n\t"\
		"addq	%%rax,%%rcx				\n\t		movq	%%r9,%%r12	\n\t	addq	%%rcx,%%r12	\n\t"\
		"addq	%%rax,%%rdx				\n\t		movq	%%r9,%%r13	\n\t	addq	%%rdx,%%r13	\n\t"\
		"/* SSE2_RADIX4_DIT_0TWIDDLE_B(r00) */\n\t	/* SSE2_RADIX4_DIT_0TWIDDLE_2NDOFTWO_B(r08) */\n\t"\
		"movaps	    (%%rax),%%xmm2		\n\t		movaps	    (%%r10),%%xmm10		\n\t"\
		"movaps	    (%%rcx),%%xmm6		\n\t		movaps	    (%%r12),%%xmm14		\n\t"\
		"movaps	0x10(%%rax),%%xmm3		\n\t		movaps	0x10(%%r10),%%xmm11		\n\t"\
		"movaps	0x10(%%rcx),%%xmm7		\n\t		movaps	0x10(%%r12),%%xmm15		\n\t"\
		"movaps	    (%%rbx),%%xmm0		\n\t		movaps	    (%%r11),%%xmm8		\n\t"\
		"movaps	    (%%rdx),%%xmm4		\n\t		movaps	    (%%r13),%%xmm12		\n\t"\
		"movaps	0x10(%%rbx),%%xmm1		\n\t		movaps	0x10(%%r11),%%xmm9		\n\t"\
		"movaps	0x10(%%rdx),%%xmm5		\n\t		movaps	0x10(%%r13),%%xmm13		\n\t"\
		"subpd	%%xmm0,%%xmm2			\n\t		subpd	%%xmm8 ,%%xmm10			\n\t"\
		"subpd	%%xmm4,%%xmm6			\n\t		subpd	%%xmm12,%%xmm14			\n\t"\
		"subpd	%%xmm1,%%xmm3			\n\t		subpd	%%xmm9 ,%%xmm11			\n\t"\
		"subpd	%%xmm5,%%xmm7			\n\t		subpd	%%xmm13,%%xmm15			\n\t"\
		"mulpd	(%%r8),%%xmm0			\n\t		mulpd	(%%r8) ,%%xmm8			\n\t"\
		"mulpd	(%%r8),%%xmm4			\n\t		mulpd	(%%r8) ,%%xmm12			\n\t"\
		"mulpd	(%%r8),%%xmm1			\n\t		mulpd	(%%r8) ,%%xmm9			\n\t"\
		"mulpd	(%%r8),%%xmm5			\n\t		mulpd	(%%r8) ,%%xmm13			\n\t"\
		"addpd	%%xmm2,%%xmm0			\n\t		addpd	%%xmm10,%%xmm8			\n\t"\
		"addpd	%%xmm6,%%xmm4			\n\t		addpd	%%xmm14,%%xmm12			\n\t"\
		"addpd	%%xmm3,%%xmm1			\n\t		addpd	%%xmm11,%%xmm9			\n\t"\
		"addpd	%%xmm7,%%xmm5			\n\t		addpd	%%xmm15,%%xmm13			\n\t"\
		"/* Finish radix-4 butterfly: */\n\t		/* Finish radix-4 butterfly: */	\n\t"\
		"subpd	%%xmm4,%%xmm0			\n\t		subpd	%%xmm12,%%xmm8			\n\t"\
		"subpd	%%xmm7,%%xmm2			\n\t		subpd	%%xmm13,%%xmm9			\n\t"\
		"subpd	%%xmm5,%%xmm1			\n\t		/*movaps	%%xmm8,0xc0(%%rsi)*/\n\t"\
		"subpd	%%xmm6,%%xmm3			\n\t		/*movaps	%%xmm9,0xd0(%%rsi)*/\n\t"\
		"/*movaps	%%xmm0,0x40(%%rsi)*/\n\t		subpd	%%xmm15,%%xmm10			\n\t"\
		"/*movaps	%%xmm2,0x60(%%rsi)*/\n\t		mulpd	(%%r8) ,%%xmm12			\n\t"\
		"/*movaps	%%xmm1,0x50(%%rsi)*/\n\t		subpd	%%xmm14,%%xmm11			\n\t"\
		"/*movaps	%%xmm3,0x30(%%rsi)*/\n\t		mulpd	(%%r8) ,%%xmm13			\n\t"\
		"mulpd	(%%r8),%%xmm4			\n\t		mulpd	(%%r8) ,%%xmm15			\n\t"\
		"mulpd	(%%r8),%%xmm7			\n\t		addpd	%%xmm8 ,%%xmm12			\n\t"\
		"mulpd	(%%r8),%%xmm5			\n\t		mulpd	(%%r8) ,%%xmm14			\n\t"\
		"mulpd	(%%r8),%%xmm6			\n\t		addpd	%%xmm9 ,%%xmm13			\n\t"\
		"addpd	%%xmm0,%%xmm4			\n\t		addpd	%%xmm10,%%xmm15			\n\t"\
		"addpd	%%xmm2,%%xmm7			\n\t		/*movaps	%%xmm12,0x80(%%rsi)*/\n\t"\
		"addpd	%%xmm1,%%xmm5			\n\t		addpd	%%xmm11,%%xmm14			\n\t"\
		"addpd	%%xmm3,%%xmm6			\n\t		/*movaps	%%xmm13,0x90(%%rsi)*/\n\t"\
		"/*movaps	%%xmm4,    (%%rsi)*/\n\t		subpd	%%xmm15,%%xmm11			\n\t"\
		"/*movaps	%%xmm7,0x20(%%rsi)*/\n\t		subpd	%%xmm10,%%xmm14			\n\t"\
		"/*movaps	%%xmm5,0x10(%%rsi)*/\n\t		mulpd	(%%r8) ,%%xmm15			\n\t"\
		"/*movaps	%%xmm6,0x70(%%rsi)*/\n\t		mulpd	(%%r8) ,%%xmm10			\n\t"\
		"subpd	%%xmm12,%%xmm4			\n\t		addpd	%%xmm11,%%xmm15			\n\t"\
		"subpd	%%xmm9 ,%%xmm0			\n\t		addpd	%%xmm14,%%xmm10			\n\t"\
		"subpd	%%xmm13,%%xmm5			\n\t		mulpd	-0x480(%%r8),%%xmm11	/* isrt2 */\n\t"\
		"subpd	%%xmm8 ,%%xmm1			\n\t		mulpd	-0x480(%%r8),%%xmm14	\n\t"\
		"mulpd	(%%r8),%%xmm12			\n\t		mulpd	-0x480(%%r8),%%xmm15	\n\t"\
		"mulpd	(%%r8),%%xmm9			\n\t		mulpd	-0x480(%%r8),%%xmm10	\n\t"\
		"mulpd	(%%r8),%%xmm13			\n\t		/*movaps	%%xmm11,0xb0(%%rsi)*/\n\t"\
		"mulpd	(%%r8),%%xmm8			\n\t		/*movaps	%%xmm14,0xf0(%%rsi)*/\n\t"\
		"addpd	%%xmm4 ,%%xmm12			\n\t		/*movaps	%%xmm15,0xa0(%%rsi)*/\n\t"\
		"addpd	%%xmm0 ,%%xmm9			\n\t		/*movaps	%%xmm10,0xe0(%%rsi)*/\n\t"\
		"/* SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS(r00,r02,r04,r06,r08,r0A,r0C,r0E) */\n\t"\
		"addpd	%%xmm5 ,%%xmm13			\n\t		subpd	%%xmm15,%%xmm7			\n\t"\
		"addpd	%%xmm1 ,%%xmm8			\n\t		subpd	%%xmm14,%%xmm2			\n\t"\
		"movaps	%%xmm4 ,0x80(%%rsi)		\n\t		subpd	%%xmm11,%%xmm3			\n\t"\
		"movaps	%%xmm0 ,0xc0(%%rsi)		\n\t		subpd	%%xmm10,%%xmm6			\n\t"\
		"movaps	%%xmm5 ,0x90(%%rsi)		\n\t		mulpd	(%%r8) ,%%xmm15			\n\t"\
		"movaps	%%xmm1 ,0x50(%%rsi)		\n\t		mulpd	(%%r8) ,%%xmm14			\n\t"\
		"movaps	%%xmm12,    (%%rsi)		\n\t		mulpd	(%%r8) ,%%xmm11			\n\t"\
		"movaps	%%xmm9 ,0x40(%%rsi)		\n\t		mulpd	(%%r8) ,%%xmm10			\n\t"\
		"movaps	%%xmm13,0x10(%%rsi)		\n\t		addpd	%%xmm7 ,%%xmm15			\n\t"\
		"movaps	%%xmm8 ,0xd0(%%rsi)		\n\t		addpd	%%xmm2 ,%%xmm14			\n\t"\
		"											addpd	%%xmm3 ,%%xmm11			\n\t"\
		"											addpd	%%xmm6 ,%%xmm10			\n\t"\
		"											movaps	%%xmm7 ,0xa0(%%rsi)		\n\t"\
		"											movaps	%%xmm2 ,0xe0(%%rsi)		\n\t"\
		"											movaps	%%xmm3 ,0xb0(%%rsi)		\n\t"\
		"											movaps	%%xmm6 ,0x70(%%rsi)		\n\t"\
		"											movaps	%%xmm15,0x20(%%rsi)		\n\t"\
		"											movaps	%%xmm14,0x60(%%rsi)		\n\t"\
		"											movaps	%%xmm11,0x30(%%rsi)		\n\t"\
		"											movaps	%%xmm10,0xf0(%%rsi)		\n\t"\
		"\n\t"\
	"/**************************...Block 2: ****************************/\n\t"\
		"addq	$0x100,%%rsi	/* r10 */	\n\t"\
		"movslq	%[__p08],%%rdi			\n\t"\
		"shlq	$3,%%rdi				\n\t"\
		"addq	%%rdi,%%rax	/* add0+p08 */\n\t		movq	%%r9,%%r10	\n\t	addq	%%rax,%%r10	\n\t"\
		"addq	%%rdi,%%rbx				\n\t		movq	%%r9,%%r11	\n\t	addq	%%rbx,%%r11	\n\t"\
		"addq	%%rdi,%%rcx				\n\t		movq	%%r9,%%r12	\n\t	addq	%%rcx,%%r12	\n\t"\
		"addq	%%rdi,%%rdx				\n\t		movq	%%r9,%%r13	\n\t	addq	%%rdx,%%r13	\n\t"\
		"/* SSE2_RADIX4_DIT_0TWIDDLE_B(r10) */\n\t	/* SSE2_RADIX4_DIT_0TWIDDLE_2NDOFTWO_B(r18) */\n\t"\
		"movaps	    (%%rax),%%xmm2		\n\t		movaps	    (%%r10),%%xmm10		\n\t"\
		"movaps	    (%%rcx),%%xmm6		\n\t		movaps	    (%%r12),%%xmm14		\n\t"\
		"movaps	0x10(%%rax),%%xmm3		\n\t		movaps	0x10(%%r10),%%xmm11		\n\t"\
		"movaps	0x10(%%rcx),%%xmm7		\n\t		movaps	0x10(%%r12),%%xmm15		\n\t"\
		"movaps	    (%%rbx),%%xmm0		\n\t		movaps	    (%%r11),%%xmm8		\n\t"\
		"movaps	    (%%rdx),%%xmm4		\n\t		movaps	    (%%r13),%%xmm12		\n\t"\
		"movaps	0x10(%%rbx),%%xmm1		\n\t		movaps	0x10(%%r11),%%xmm9		\n\t"\
		"movaps	0x10(%%rdx),%%xmm5		\n\t		movaps	0x10(%%r13),%%xmm13		\n\t"\
		"subpd	%%xmm0,%%xmm2			\n\t		subpd	%%xmm8 ,%%xmm10			\n\t"\
		"subpd	%%xmm4,%%xmm6			\n\t		subpd	%%xmm12,%%xmm14			\n\t"\
		"subpd	%%xmm1,%%xmm3			\n\t		subpd	%%xmm9 ,%%xmm11			\n\t"\
		"subpd	%%xmm5,%%xmm7			\n\t		subpd	%%xmm13,%%xmm15			\n\t"\
		"mulpd	(%%r8),%%xmm0			\n\t		mulpd	(%%r8) ,%%xmm8			\n\t"\
		"mulpd	(%%r8),%%xmm4			\n\t		mulpd	(%%r8) ,%%xmm12			\n\t"\
		"mulpd	(%%r8),%%xmm1			\n\t		mulpd	(%%r8) ,%%xmm9			\n\t"\
		"mulpd	(%%r8),%%xmm5			\n\t		mulpd	(%%r8) ,%%xmm13			\n\t"\
		"addpd	%%xmm2,%%xmm0			\n\t		addpd	%%xmm10,%%xmm8			\n\t"\
		"addpd	%%xmm6,%%xmm4			\n\t		addpd	%%xmm14,%%xmm12			\n\t"\
		"addpd	%%xmm3,%%xmm1			\n\t		addpd	%%xmm11,%%xmm9			\n\t"\
		"addpd	%%xmm7,%%xmm5			\n\t		addpd	%%xmm15,%%xmm13			\n\t"\
		"/* Finish radix-4 butterfly: */\n\t		/* Finish radix-4 butterfly: */	\n\t"\
		"subpd	%%xmm4,%%xmm0			\n\t		subpd	%%xmm12,%%xmm8			\n\t"\
		"subpd	%%xmm7,%%xmm2			\n\t		subpd	%%xmm13,%%xmm9			\n\t"\
		"subpd	%%xmm5,%%xmm1			\n\t		subpd	%%xmm15,%%xmm10			\n\t"\
		"subpd	%%xmm6,%%xmm3			\n\t		mulpd	(%%r8) ,%%xmm12			\n\t"\
		"mulpd	(%%r8),%%xmm4			\n\t		subpd	%%xmm14,%%xmm11			\n\t"\
		"mulpd	(%%r8),%%xmm7			\n\t		mulpd	(%%r8) ,%%xmm13			\n\t"\
		"mulpd	(%%r8),%%xmm5			\n\t		mulpd	(%%r8) ,%%xmm15			\n\t"\
		"mulpd	(%%r8),%%xmm6			\n\t		addpd	%%xmm8 ,%%xmm12			\n\t"\
		"addpd	%%xmm0,%%xmm4			\n\t		mulpd	(%%r8) ,%%xmm14			\n\t"\
		"addpd	%%xmm2,%%xmm7			\n\t		addpd	%%xmm9 ,%%xmm13			\n\t"\
		"addpd	%%xmm1,%%xmm5			\n\t		addpd	%%xmm10,%%xmm15			\n\t"\
		"addpd	%%xmm3,%%xmm6			\n\t		addpd	%%xmm11,%%xmm14			\n\t"\
		"/* SSE2_RADIX8_DIT_COMBINE_RAD4_SUB...*/	subpd	%%xmm15,%%xmm11			\n\t"\
		"subpd	%%xmm12,%%xmm4			\n\t		subpd	%%xmm10,%%xmm14			\n\t"\
		"subpd	%%xmm9 ,%%xmm0			\n\t		mulpd	(%%r8) ,%%xmm15			\n\t"\
		"subpd	%%xmm13,%%xmm5			\n\t		mulpd	(%%r8) ,%%xmm10			\n\t"\
		"subpd	%%xmm8 ,%%xmm1			\n\t		addpd	%%xmm11,%%xmm15			\n\t"\
		"mulpd	(%%r8),%%xmm12			\n\t		addpd	%%xmm14,%%xmm10			\n\t"\
		"mulpd	(%%r8),%%xmm9			\n\t		mulpd	-0x480(%%r8),%%xmm11			\n\t"\
		"mulpd	(%%r8),%%xmm13			\n\t		mulpd	-0x480(%%r8),%%xmm14			\n\t"\
		"mulpd	(%%r8),%%xmm8			\n\t		mulpd	-0x480(%%r8),%%xmm15			\n\t"\
		"addpd	%%xmm4 ,%%xmm12			\n\t		mulpd	-0x480(%%r8),%%xmm10			\n\t"\
		"addpd	%%xmm0 ,%%xmm9			\n\t	/*...S(r00,r02,r04,r06,r08,r0A,r0C,r0E) */\n\t"\
		"addpd	%%xmm5 ,%%xmm13			\n\t		subpd	%%xmm15,%%xmm7			\n\t"\
		"addpd	%%xmm1 ,%%xmm8			\n\t		subpd	%%xmm14,%%xmm2			\n\t"\
		"movaps	%%xmm4 ,0x80(%%rsi)		\n\t		subpd	%%xmm11,%%xmm3			\n\t"\
		"movaps	%%xmm0 ,0xc0(%%rsi)		\n\t		subpd	%%xmm10,%%xmm6			\n\t"\
		"movaps	%%xmm5 ,0x90(%%rsi)		\n\t		mulpd	(%%r8) ,%%xmm15			\n\t"\
		"movaps	%%xmm1 ,0x50(%%rsi)		\n\t		mulpd	(%%r8) ,%%xmm14			\n\t"\
		"movaps	%%xmm12,    (%%rsi)		\n\t		mulpd	(%%r8) ,%%xmm11			\n\t"\
		"movaps	%%xmm9 ,0x40(%%rsi)		\n\t		mulpd	(%%r8) ,%%xmm10			\n\t"\
		"movaps	%%xmm13,0x10(%%rsi)		\n\t		addpd	%%xmm7 ,%%xmm15			\n\t"\
		"movaps	%%xmm8 ,0xd0(%%rsi)		\n\t		addpd	%%xmm2 ,%%xmm14			\n\t"\
		"											addpd	%%xmm3 ,%%xmm11			\n\t"\
		"											addpd	%%xmm6 ,%%xmm10			\n\t"\
		"											movaps	%%xmm7 ,0xa0(%%rsi)		\n\t"\
		"											movaps	%%xmm2 ,0xe0(%%rsi)		\n\t"\
		"											movaps	%%xmm3 ,0xb0(%%rsi)		\n\t"\
		"											movaps	%%xmm6 ,0x70(%%rsi)		\n\t"\
		"											movaps	%%xmm15,0x20(%%rsi)		\n\t"\
		"											movaps	%%xmm14,0x60(%%rsi)		\n\t"\
		"											movaps	%%xmm11,0x30(%%rsi)		\n\t"\
		"											movaps	%%xmm10,0xf0(%%rsi)		\n\t"\
		"\n\t"\
	"/**************************...Block 3: ****************************/\n\t"\
		"addq	$0x100,%%rsi	/* r20 */	\n\t"\
		"movslq	%[__p10],%%r14			\n\t"\
		"shlq	$3,%%r14				\n\t"\
		"subq	%%rdi,%%r14	/* p10-p8 */\n\t"\
		"addq	%%r14,%%rax	/* add0+p10 */\n\t		movq	%%r9,%%r10	\n\t	addq	%%rax,%%r10	\n\t"\
		"addq	%%r14,%%rbx				\n\t		movq	%%r9,%%r11	\n\t	addq	%%rbx,%%r11	\n\t"\
		"addq	%%r14,%%rcx				\n\t		movq	%%r9,%%r12	\n\t	addq	%%rcx,%%r12	\n\t"\
		"addq	%%r14,%%rdx				\n\t		movq	%%r9,%%r13	\n\t	addq	%%rdx,%%r13	\n\t"\
		"/* SSE2_RADIX4_DIT_0TWIDDLE_B(r20) */\n\t	/* SSE2_RADIX4_DIT_0TWIDDLE_2NDOFTWO_B(r28) */\n\t"\
		"movaps	    (%%rax),%%xmm2		\n\t		movaps	    (%%r10),%%xmm10		\n\t"\
		"movaps	    (%%rcx),%%xmm6		\n\t		movaps	    (%%r12),%%xmm14		\n\t"\
		"movaps	0x10(%%rax),%%xmm3		\n\t		movaps	0x10(%%r10),%%xmm11		\n\t"\
		"movaps	0x10(%%rcx),%%xmm7		\n\t		movaps	0x10(%%r12),%%xmm15		\n\t"\
		"movaps	    (%%rbx),%%xmm0		\n\t		movaps	    (%%r11),%%xmm8		\n\t"\
		"movaps	    (%%rdx),%%xmm4		\n\t		movaps	    (%%r13),%%xmm12		\n\t"\
		"movaps	0x10(%%rbx),%%xmm1		\n\t		movaps	0x10(%%r11),%%xmm9		\n\t"\
		"movaps	0x10(%%rdx),%%xmm5		\n\t		movaps	0x10(%%r13),%%xmm13		\n\t"\
		"subpd	%%xmm0,%%xmm2			\n\t		subpd	%%xmm8 ,%%xmm10			\n\t"\
		"subpd	%%xmm4,%%xmm6			\n\t		subpd	%%xmm12,%%xmm14			\n\t"\
		"subpd	%%xmm1,%%xmm3			\n\t		subpd	%%xmm9 ,%%xmm11			\n\t"\
		"subpd	%%xmm5,%%xmm7			\n\t		subpd	%%xmm13,%%xmm15			\n\t"\
		"mulpd	(%%r8),%%xmm0			\n\t		mulpd	(%%r8) ,%%xmm8			\n\t"\
		"mulpd	(%%r8),%%xmm4			\n\t		mulpd	(%%r8) ,%%xmm12			\n\t"\
		"mulpd	(%%r8),%%xmm1			\n\t		mulpd	(%%r8) ,%%xmm9			\n\t"\
		"mulpd	(%%r8),%%xmm5			\n\t		mulpd	(%%r8) ,%%xmm13			\n\t"\
		"addpd	%%xmm2,%%xmm0			\n\t		addpd	%%xmm10,%%xmm8			\n\t"\
		"addpd	%%xmm6,%%xmm4			\n\t		addpd	%%xmm14,%%xmm12			\n\t"\
		"addpd	%%xmm3,%%xmm1			\n\t		addpd	%%xmm11,%%xmm9			\n\t"\
		"addpd	%%xmm7,%%xmm5			\n\t		addpd	%%xmm15,%%xmm13			\n\t"\
		"/* Finish radix-4 butterfly: */\n\t		/* Finish radix-4 butterfly: */	\n\t"\
		"subpd	%%xmm4,%%xmm0			\n\t		subpd	%%xmm12,%%xmm8			\n\t"\
		"subpd	%%xmm7,%%xmm2			\n\t		subpd	%%xmm13,%%xmm9			\n\t"\
		"subpd	%%xmm5,%%xmm1			\n\t		subpd	%%xmm15,%%xmm10			\n\t"\
		"subpd	%%xmm6,%%xmm3			\n\t		mulpd	(%%r8) ,%%xmm12			\n\t"\
		"mulpd	(%%r8),%%xmm4			\n\t		subpd	%%xmm14,%%xmm11			\n\t"\
		"mulpd	(%%r8),%%xmm7			\n\t		mulpd	(%%r8) ,%%xmm13			\n\t"\
		"mulpd	(%%r8),%%xmm5			\n\t		mulpd	(%%r8) ,%%xmm15			\n\t"\
		"mulpd	(%%r8),%%xmm6			\n\t		addpd	%%xmm8 ,%%xmm12			\n\t"\
		"addpd	%%xmm0,%%xmm4			\n\t		mulpd	(%%r8) ,%%xmm14			\n\t"\
		"addpd	%%xmm2,%%xmm7			\n\t		addpd	%%xmm9 ,%%xmm13			\n\t"\
		"addpd	%%xmm1,%%xmm5			\n\t		addpd	%%xmm10,%%xmm15			\n\t"\
		"addpd	%%xmm3,%%xmm6			\n\t		addpd	%%xmm11,%%xmm14			\n\t"\
		"/* SSE2_RADIX8_DIT_COMBINE_RAD4_SUB...*/	subpd	%%xmm15,%%xmm11			\n\t"\
		"subpd	%%xmm12,%%xmm4			\n\t		subpd	%%xmm10,%%xmm14			\n\t"\
		"subpd	%%xmm9 ,%%xmm0			\n\t		mulpd	(%%r8) ,%%xmm15			\n\t"\
		"subpd	%%xmm13,%%xmm5			\n\t		mulpd	(%%r8) ,%%xmm10			\n\t"\
		"subpd	%%xmm8 ,%%xmm1			\n\t		addpd	%%xmm11,%%xmm15			\n\t"\
		"mulpd	(%%r8),%%xmm12			\n\t		addpd	%%xmm14,%%xmm10			\n\t"\
		"mulpd	(%%r8),%%xmm9			\n\t		mulpd	-0x480(%%r8),%%xmm11			\n\t"\
		"mulpd	(%%r8),%%xmm13			\n\t		mulpd	-0x480(%%r8),%%xmm14			\n\t"\
		"mulpd	(%%r8),%%xmm8			\n\t		mulpd	-0x480(%%r8),%%xmm15			\n\t"\
		"addpd	%%xmm4 ,%%xmm12			\n\t		mulpd	-0x480(%%r8),%%xmm10			\n\t"\
		"addpd	%%xmm0 ,%%xmm9			\n\t	/*...S(r00,r02,r04,r06,r08,r0A,r0C,r0E) */\n\t"\
		"addpd	%%xmm5 ,%%xmm13			\n\t		subpd	%%xmm15,%%xmm7			\n\t"\
		"addpd	%%xmm1 ,%%xmm8			\n\t		subpd	%%xmm14,%%xmm2			\n\t"\
		"movaps	%%xmm4 ,0x80(%%rsi)		\n\t		subpd	%%xmm11,%%xmm3			\n\t"\
		"movaps	%%xmm0 ,0xc0(%%rsi)		\n\t		subpd	%%xmm10,%%xmm6			\n\t"\
		"movaps	%%xmm5 ,0x90(%%rsi)		\n\t		mulpd	(%%r8) ,%%xmm15			\n\t"\
		"movaps	%%xmm1 ,0x50(%%rsi)		\n\t		mulpd	(%%r8) ,%%xmm14			\n\t"\
		"movaps	%%xmm12,    (%%rsi)		\n\t		mulpd	(%%r8) ,%%xmm11			\n\t"\
		"movaps	%%xmm9 ,0x40(%%rsi)		\n\t		mulpd	(%%r8) ,%%xmm10			\n\t"\
		"movaps	%%xmm13,0x10(%%rsi)		\n\t		addpd	%%xmm7 ,%%xmm15			\n\t"\
		"movaps	%%xmm8 ,0xd0(%%rsi)		\n\t		addpd	%%xmm2 ,%%xmm14			\n\t"\
		"											addpd	%%xmm3 ,%%xmm11			\n\t"\
		"											addpd	%%xmm6 ,%%xmm10			\n\t"\
		"											movaps	%%xmm7 ,0xa0(%%rsi)		\n\t"\
		"											movaps	%%xmm2 ,0xe0(%%rsi)		\n\t"\
		"											movaps	%%xmm3 ,0xb0(%%rsi)		\n\t"\
		"											movaps	%%xmm6 ,0x70(%%rsi)		\n\t"\
		"											movaps	%%xmm15,0x20(%%rsi)		\n\t"\
		"											movaps	%%xmm14,0x60(%%rsi)		\n\t"\
		"											movaps	%%xmm11,0x30(%%rsi)		\n\t"\
		"											movaps	%%xmm10,0xf0(%%rsi)		\n\t"\
		"\n\t"\
	"/**************************...Block 4: ****************************/\n\t"\
		"addq	$0x100,%%rsi	/* r30 */	\n\t"\
		"movslq	%[__p08],%%rdi			\n\t"\
		"shlq	$3,%%rdi				\n\t"\
		"addq	%%rdi,%%rax	/* add0+p18 */\n\t		movq	%%r9,%%r10	\n\t	addq	%%rax,%%r10	\n\t"\
		"addq	%%rdi,%%rbx				\n\t		movq	%%r9,%%r11	\n\t	addq	%%rbx,%%r11	\n\t"\
		"addq	%%rdi,%%rcx				\n\t		movq	%%r9,%%r12	\n\t	addq	%%rcx,%%r12	\n\t"\
		"addq	%%rdi,%%rdx				\n\t		movq	%%r9,%%r13	\n\t	addq	%%rdx,%%r13	\n\t"\
		"/* SSE2_RADIX4_DIT_0TWIDDLE_B(r30) */\n\t	/* SSE2_RADIX4_DIT_0TWIDDLE_2NDOFTWO_B(r38) */\n\t"\
		"movaps	    (%%rax),%%xmm2		\n\t		movaps	    (%%r10),%%xmm10		\n\t"\
		"movaps	    (%%rcx),%%xmm6		\n\t		movaps	    (%%r12),%%xmm14		\n\t"\
		"movaps	0x10(%%rax),%%xmm3		\n\t		movaps	0x10(%%r10),%%xmm11		\n\t"\
		"movaps	0x10(%%rcx),%%xmm7		\n\t		movaps	0x10(%%r12),%%xmm15		\n\t"\
		"movaps	    (%%rbx),%%xmm0		\n\t		movaps	    (%%r11),%%xmm8		\n\t"\
		"movaps	    (%%rdx),%%xmm4		\n\t		movaps	    (%%r13),%%xmm12		\n\t"\
		"movaps	0x10(%%rbx),%%xmm1		\n\t		movaps	0x10(%%r11),%%xmm9		\n\t"\
		"movaps	0x10(%%rdx),%%xmm5		\n\t		movaps	0x10(%%r13),%%xmm13		\n\t"\
		"subpd	%%xmm0,%%xmm2			\n\t		subpd	%%xmm8 ,%%xmm10			\n\t"\
		"subpd	%%xmm4,%%xmm6			\n\t		subpd	%%xmm12,%%xmm14			\n\t"\
		"subpd	%%xmm1,%%xmm3			\n\t		subpd	%%xmm9 ,%%xmm11			\n\t"\
		"subpd	%%xmm5,%%xmm7			\n\t		subpd	%%xmm13,%%xmm15			\n\t"\
		"mulpd	(%%r8),%%xmm0			\n\t		mulpd	(%%r8) ,%%xmm8			\n\t"\
		"mulpd	(%%r8),%%xmm4			\n\t		mulpd	(%%r8) ,%%xmm12			\n\t"\
		"mulpd	(%%r8),%%xmm1			\n\t		mulpd	(%%r8) ,%%xmm9			\n\t"\
		"mulpd	(%%r8),%%xmm5			\n\t		mulpd	(%%r8) ,%%xmm13			\n\t"\
		"addpd	%%xmm2,%%xmm0			\n\t		addpd	%%xmm10,%%xmm8			\n\t"\
		"addpd	%%xmm6,%%xmm4			\n\t		addpd	%%xmm14,%%xmm12			\n\t"\
		"addpd	%%xmm3,%%xmm1			\n\t		addpd	%%xmm11,%%xmm9			\n\t"\
		"addpd	%%xmm7,%%xmm5			\n\t		addpd	%%xmm15,%%xmm13			\n\t"\
		"/* Finish radix-4 butterfly: */\n\t		/* Finish radix-4 butterfly: */	\n\t"\
		"subpd	%%xmm4,%%xmm0			\n\t		subpd	%%xmm12,%%xmm8			\n\t"\
		"subpd	%%xmm7,%%xmm2			\n\t		subpd	%%xmm13,%%xmm9			\n\t"\
		"subpd	%%xmm5,%%xmm1			\n\t		subpd	%%xmm15,%%xmm10			\n\t"\
		"subpd	%%xmm6,%%xmm3			\n\t		mulpd	(%%r8) ,%%xmm12			\n\t"\
		"mulpd	(%%r8),%%xmm4			\n\t		subpd	%%xmm14,%%xmm11			\n\t"\
		"mulpd	(%%r8),%%xmm7			\n\t		mulpd	(%%r8) ,%%xmm13			\n\t"\
		"mulpd	(%%r8),%%xmm5			\n\t		mulpd	(%%r8) ,%%xmm15			\n\t"\
		"mulpd	(%%r8),%%xmm6			\n\t		addpd	%%xmm8 ,%%xmm12			\n\t"\
		"addpd	%%xmm0,%%xmm4			\n\t		mulpd	(%%r8) ,%%xmm14			\n\t"\
		"addpd	%%xmm2,%%xmm7			\n\t		addpd	%%xmm9 ,%%xmm13			\n\t"\
		"addpd	%%xmm1,%%xmm5			\n\t		addpd	%%xmm10,%%xmm15			\n\t"\
		"addpd	%%xmm3,%%xmm6			\n\t		addpd	%%xmm11,%%xmm14			\n\t"\
		"/* SSE2_RADIX8_DIT_COMBINE_RAD4_SUB...*/	subpd	%%xmm15,%%xmm11			\n\t"\
		"subpd	%%xmm12,%%xmm4			\n\t		subpd	%%xmm10,%%xmm14			\n\t"\
		"subpd	%%xmm9 ,%%xmm0			\n\t		mulpd	(%%r8) ,%%xmm15			\n\t"\
		"subpd	%%xmm13,%%xmm5			\n\t		mulpd	(%%r8) ,%%xmm10			\n\t"\
		"subpd	%%xmm8 ,%%xmm1			\n\t		addpd	%%xmm11,%%xmm15			\n\t"\
		"mulpd	(%%r8),%%xmm12			\n\t		addpd	%%xmm14,%%xmm10			\n\t"\
		"mulpd	(%%r8),%%xmm9			\n\t		mulpd	-0x480(%%r8),%%xmm11			\n\t"\
		"mulpd	(%%r8),%%xmm13			\n\t		mulpd	-0x480(%%r8),%%xmm14			\n\t"\
		"mulpd	(%%r8),%%xmm8			\n\t		mulpd	-0x480(%%r8),%%xmm15			\n\t"\
		"addpd	%%xmm4 ,%%xmm12			\n\t		mulpd	-0x480(%%r8),%%xmm10			\n\t"\
		"addpd	%%xmm0 ,%%xmm9			\n\t	/*...S(r00,r02,r04,r06,r08,r0A,r0C,r0E) */\n\t"\
		"addpd	%%xmm5 ,%%xmm13			\n\t		subpd	%%xmm15,%%xmm7			\n\t"\
		"addpd	%%xmm1 ,%%xmm8			\n\t		subpd	%%xmm14,%%xmm2			\n\t"\
		"movaps	%%xmm4 ,0x80(%%rsi)		\n\t		subpd	%%xmm11,%%xmm3			\n\t"\
		"movaps	%%xmm0 ,0xc0(%%rsi)		\n\t		subpd	%%xmm10,%%xmm6			\n\t"\
		"movaps	%%xmm5 ,0x90(%%rsi)		\n\t		mulpd	(%%r8) ,%%xmm15			\n\t"\
		"movaps	%%xmm1 ,0x50(%%rsi)		\n\t		mulpd	(%%r8) ,%%xmm14			\n\t"\
		"movaps	%%xmm12,    (%%rsi)		\n\t		mulpd	(%%r8) ,%%xmm11			\n\t"\
		"movaps	%%xmm9 ,0x40(%%rsi)		\n\t		mulpd	(%%r8) ,%%xmm10			\n\t"\
		"movaps	%%xmm13,0x10(%%rsi)		\n\t		addpd	%%xmm7 ,%%xmm15			\n\t"\
		"movaps	%%xmm8 ,0xd0(%%rsi)		\n\t		addpd	%%xmm2 ,%%xmm14			\n\t"\
		"											addpd	%%xmm3 ,%%xmm11			\n\t"\
		"											addpd	%%xmm6 ,%%xmm10			\n\t"\
		"											movaps	%%xmm7 ,0xa0(%%rsi)		\n\t"\
		"											movaps	%%xmm2 ,0xe0(%%rsi)		\n\t"\
		"											movaps	%%xmm3 ,0xb0(%%rsi)		\n\t"\
		"											movaps	%%xmm6 ,0x70(%%rsi)		\n\t"\
		"											movaps	%%xmm15,0x20(%%rsi)		\n\t"\
		"											movaps	%%xmm14,0x60(%%rsi)		\n\t"\
		"											movaps	%%xmm11,0x30(%%rsi)		\n\t"\
		"											movaps	%%xmm10,0xf0(%%rsi)		\n\t"\
		"\n\t"\
		"/**********************************************************************************/\n\t"\
		"/*...and now do eight radix-4 transforms, including the internal twiddle factors: */\n\t"\
		"/**********************************************************************************/\n\t"\
		"\n\t"\
		"movslq	%[__p10],%%rdi	/* edi will store copy of p10 throughout */\n\t"\
		"/*...Block 1: t00,t10,t20,t30	*/				/*...Block 5: t08,t18,t28,t38*/\n\t"\
		"movq	%[__add0],%%rax							\n\t		movslq	%[__p04],%%rsi		\n\t"\
		"movslq	%[__p08],%%rbx							\n\t		movq	%%rax,%%r10			\n\t"\
		"movq	%[__r00],%%rcx							\n\t		shlq	$3,%%rsi			\n\t"\
		"movq	%[__r10],%%rdx							\n\t		addq	%%rsi	,%%r10	/* add0 = &a[j1+p4] */\n\t"\
		"shlq	$3,%%rbx								\n\t		movq	%[__isrt2],%%rsi	\n\t"\
		"shlq	$3,%%rdi								\n\t		movq	%%rbx,%%r11		/* Need this register-copy before add0 gets added to rbx at left */	\n\t"\
		"addq	%%rax	,%%rbx	/* add1 = add0+p8 */	\n\t		movaps	(%%rsi),%%xmm10	/* isrt2 */\n\t"\
		"												\n\t		addq	%%r10	,%%r11	/* add1 = add0+p12 */\n\t"\
		"movaps	     (%%rdx),%%xmm2						\n\t		movaps	0x280(%%rcx),%%xmm12	\n\t"\
		"movaps	0x200(%%rdx),%%xmm4						\n\t		movaps	0x290(%%rcx),%%xmm13	\n\t"\
		"movaps	0x010(%%rdx),%%xmm3						\n\t		movaps	0x280(%%rdx),%%xmm14	\n\t"\
		"movaps	0x210(%%rdx),%%xmm5						\n\t		movaps	0x290(%%rdx),%%xmm15	\n\t"\
		"movaps	     (%%rcx),%%xmm0						\n\t		mulpd	%%xmm10,%%xmm12		\n\t"\
		"movaps	0x200(%%rcx),%%xmm6						\n\t		mulpd	%%xmm10,%%xmm13		\n\t"\
		"movaps	0x010(%%rcx),%%xmm1						\n\t		mulpd	%%xmm10,%%xmm14		\n\t"\
		"movaps	0x210(%%rcx),%%xmm7						\n\t		mulpd	%%xmm10,%%xmm15		\n\t"\
		"\n\t"\
		"subpd	%%xmm2,%%xmm0							\n\t		subpd	%%xmm12,%%xmm13		\n\t"\
		"subpd	%%xmm4,%%xmm6							\n\t		movaps	0x080(%%rcx),%%xmm8	\n\t"\
		"subpd	%%xmm3,%%xmm1							\n\t		subpd	%%xmm15,%%xmm14		\n\t"\
		"subpd	%%xmm5,%%xmm7							\n\t		movaps	0x090(%%rdx),%%xmm10	\n\t"\
		"mulpd	(%%r8),%%xmm2							\n\t		mulpd	(%%r8),%%xmm12		\n\t"\
		"mulpd	(%%r8),%%xmm4							\n\t		movaps	0x090(%%rcx),%%xmm11	\n\t"\
		"mulpd	(%%r8),%%xmm3							\n\t		mulpd	(%%r8),%%xmm15		\n\t"\
		"mulpd	(%%r8),%%xmm5							\n\t		movaps	0x080(%%rdx),%%xmm9	\n\t"\
		"addpd	%%xmm0,%%xmm2							\n\t		addpd	%%xmm13,%%xmm12		\n\t"\
		"addpd	%%xmm6,%%xmm4							\n\t		addpd	%%xmm14,%%xmm15		\n\t"\
		"addpd	%%xmm1,%%xmm3							\n\t"\
		"addpd	%%xmm7,%%xmm5							\n\t		subpd	%%xmm14,%%xmm12		\n\t"\
		"												\n\t		subpd	%%xmm10,%%xmm8		\n\t"\
		"/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c00): */\n\t		subpd	%%xmm15,%%xmm13		\n\t"\
		"subpd	%%xmm4,%%xmm2							\n\t		subpd	%%xmm9 ,%%xmm11		\n\t"\
		"subpd	%%xmm7,%%xmm0							\n\t		mulpd	(%%r8),%%xmm14		\n\t"\
		"subpd	%%xmm5,%%xmm3							\n\t		mulpd	(%%r8),%%xmm10		\n\t"\
		"subpd	%%xmm6,%%xmm1							\n\t		mulpd	(%%r8),%%xmm15		\n\t"\
		"mulpd	(%%r8),%%xmm4							\n\t		mulpd	(%%r8),%%xmm9		\n\t"\
		"mulpd	(%%r8),%%xmm7							\n\t		addpd	%%xmm12,%%xmm14		\n\t"\
		"mulpd	(%%r8),%%xmm5							\n\t		addpd	%%xmm8 ,%%xmm10		\n\t"\
		"mulpd	(%%r8),%%xmm6							\n\t		addpd	%%xmm13,%%xmm15		\n\t"\
		"addpd	%%xmm2,%%xmm4							\n\t		addpd	%%xmm11,%%xmm9		\n\t"\
		"addpd	%%xmm0,%%xmm7							/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c04): */\n\t"\
		"addpd	%%xmm3,%%xmm5							/* c04 = c00 + 0x80: */\n\t"\
		"addpd	%%xmm1,%%xmm6							\n\t		subpd	%%xmm12,%%xmm10		\n\t"\
		"addq	$0x070,%%rsi	/* c00 */				\n\t		subpd	%%xmm15,%%xmm8		\n\t"\
		"movaps	%%xmm2,     (%%rcx)						\n\t		subpd	%%xmm13,%%xmm11		\n\t"\
		"movaps	%%xmm0,0x010(%%rdx)						\n\t		subpd	%%xmm14,%%xmm9		\n\t"\
		"movaps	%%xmm3,0x010(%%rcx)						\n\t		mulpd	(%%r8),%%xmm12		\n\t"\
		"movaps	%%xmm6,     (%%rdx)						\n\t		mulpd	(%%r8),%%xmm15		\n\t"\
		"movaps	%%xmm4,%%xmm2							\n\t		mulpd	(%%r8),%%xmm13		\n\t"\
		"movaps	%%xmm7,%%xmm0							\n\t		mulpd	(%%r8),%%xmm14		\n\t"\
		"movaps	%%xmm5,%%xmm3							\n\t		addpd	%%xmm10,%%xmm12		\n\t"\
		"movaps	%%xmm1,%%xmm6							\n\t		addpd	%%xmm8 ,%%xmm15		\n\t"\
		"mulpd	     (%%rsi),%%xmm4						\n\t		addpd	%%xmm11,%%xmm13		\n\t"\
		"mulpd	0x020(%%rsi),%%xmm7						\n\t		addpd	%%xmm9 ,%%xmm14		\n\t"\
		"mulpd	     (%%rsi),%%xmm5						\n\t		movaps	%%xmm10,0x080(%%rcx)	\n\t"\
		"mulpd	0x020(%%rsi),%%xmm1						\n\t		movaps	%%xmm8 ,0x090(%%rdx)	\n\t"\
		"mulpd	0x010(%%rsi),%%xmm2						\n\t		movaps	%%xmm11,0x090(%%rcx)	\n\t"\
		"mulpd	0x030(%%rsi),%%xmm0						\n\t		movaps	%%xmm14,0x080(%%rdx)	\n\t"\
		"mulpd	0x010(%%rsi),%%xmm3						\n\t		movaps	%%xmm12,%%xmm10		\n\t"\
		"mulpd	0x030(%%rsi),%%xmm6						\n\t		movaps	%%xmm15,%%xmm8		\n\t"\
		"subpd	%%xmm2,%%xmm5							\n\t		movaps	%%xmm13,%%xmm11		\n\t"\
		"subpd	%%xmm0,%%xmm1							\n\t		movaps	%%xmm9 ,%%xmm14		\n\t"\
		"addpd	%%xmm3,%%xmm4							\n\t		mulpd	0x080(%%rsi),%%xmm12	\n\t"\
		"addpd	%%xmm6,%%xmm7							\n\t		mulpd	0x0a0(%%rsi),%%xmm15	\n\t"\
		"movaps	%%xmm5,0x10(%%rax)						\n\t		mulpd	0x080(%%rsi),%%xmm13	\n\t"\
		"movaps	%%xmm1,0x10(%%rbx)						\n\t		mulpd	0x0a0(%%rsi),%%xmm9	\n\t"\
		"movaps	%%xmm4,    (%%rax)						\n\t		mulpd	0x090(%%rsi),%%xmm10	\n\t"\
		"movaps	%%xmm7,    (%%rbx)						\n\t		mulpd	0x0b0(%%rsi),%%xmm8	\n\t"\
		"addq	%%rdi,%%rax								\n\t		mulpd	0x090(%%rsi),%%xmm11	\n\t"\
		"addq	%%rdi,%%rbx								\n\t		mulpd	0x0b0(%%rsi),%%xmm14	\n\t"\
		"addq	$0x40,%%rsi		/* c10 */				/* c14 = c10 + 0x80: */\n\t"\
		"movaps	     (%%rcx),%%xmm4						\n\t		subpd	%%xmm10,%%xmm13		\n\t"\
		"movaps	0x010(%%rdx),%%xmm0						\n\t		subpd	%%xmm8 ,%%xmm9		\n\t"\
		"movaps	0x010(%%rcx),%%xmm5						\n\t		addpd	%%xmm11,%%xmm12		\n\t"\
		"movaps	     (%%rdx),%%xmm6						\n\t		addpd	%%xmm14,%%xmm15		\n\t"\
		"movaps	%%xmm4,%%xmm2							\n\t		movaps	%%xmm13,0x10(%%r10)	\n\t"\
		"movaps	%%xmm0,%%xmm1							\n\t		movaps	%%xmm9 ,0x10(%%r11)	\n\t"\
		"movaps	%%xmm5,%%xmm3							\n\t		movaps	%%xmm12,    (%%r10)	\n\t"\
		"movaps	%%xmm6,%%xmm7							\n\t		movaps	%%xmm15,    (%%r11)	\n\t"\
		"mulpd	     (%%rsi),%%xmm4						\n\t		addq	%%rdi,%%r10			\n\t"\
		"mulpd	0x020(%%rsi),%%xmm0						\n\t		addq	%%rdi,%%r11			\n\t"\
		"mulpd	     (%%rsi),%%xmm5						\n\t		movaps	0x080(%%rcx),%%xmm12	\n\t"\
		"mulpd	0x020(%%rsi),%%xmm6						\n\t		movaps	0x090(%%rdx),%%xmm8	\n\t"\
		"mulpd	0x010(%%rsi),%%xmm2						\n\t		movaps	0x090(%%rcx),%%xmm13	\n\t"\
		"mulpd	0x030(%%rsi),%%xmm1						\n\t		movaps	0x080(%%rdx),%%xmm14	\n\t"\
		"mulpd	0x010(%%rsi),%%xmm3						\n\t		movaps	%%xmm12,%%xmm10		\n\t"\
		"mulpd	0x030(%%rsi),%%xmm7						\n\t		movaps	%%xmm8 ,%%xmm9		\n\t"\
		"subpd	%%xmm2,%%xmm5							\n\t		movaps	%%xmm13,%%xmm11		\n\t"\
		"subpd	%%xmm1,%%xmm6							\n\t		movaps	%%xmm14,%%xmm15		\n\t"\
		"addpd	%%xmm3,%%xmm4							\n\t		mulpd	0x080(%%rsi),%%xmm12	\n\t"\
		"addpd	%%xmm7,%%xmm0							\n\t		mulpd	0x0a0(%%rsi),%%xmm8	\n\t"\
		"movaps	%%xmm5,0x10(%%rax)						\n\t		mulpd	0x080(%%rsi),%%xmm13	\n\t"\
		"movaps	%%xmm6,0x10(%%rbx)						\n\t		mulpd	0x0a0(%%rsi),%%xmm14	\n\t"\
		"movaps	%%xmm4,    (%%rax)						\n\t		mulpd	0x090(%%rsi),%%xmm10	\n\t"\
		"movaps	%%xmm0,    (%%rbx)						\n\t		mulpd	0x0b0(%%rsi),%%xmm9	\n\t"\
		"												\n\t		mulpd	0x090(%%rsi),%%xmm11	\n\t"\
		"												\n\t		mulpd	0x0b0(%%rsi),%%xmm15	\n\t"\
		"												\n\t		subpd	%%xmm10,%%xmm13		\n\t"\
		"												\n\t		subpd	%%xmm9 ,%%xmm14		\n\t"\
		"												\n\t		addpd	%%xmm11,%%xmm12		\n\t"\
		"												\n\t		addpd	%%xmm15,%%xmm8		\n\t"\
		"												\n\t		movaps	%%xmm13,0x10(%%r10)	\n\t"\
		"												\n\t		movaps	%%xmm14,0x10(%%r11)	\n\t"\
		"												\n\t		movaps	%%xmm12,    (%%r10)	\n\t"\
		"												\n\t		movaps	%%xmm8 ,    (%%r11)	\n\t"\
		"\n\t"\
		"/*...Block 2: t02,t12,t22,t32	*/				/*...Block 6: t0A,t1A,t2A,t3A*/\n\t"\
		"movq	%[__add0],%%rax							\n\t		movq	%%rax,%%r10			\n\t"\
		"movslq	%[__p01],%%rsi							\n\t		movslq	%[__p05],%%r9		\n\t"\
		"movslq	%[__p08],%%rbx							\n\t		shlq	$3,%%r9			\n\t"\
		"shlq	$3,%%rsi								\n\t		addq	%%r9 	,%%r10	/* add0 = &a[j1+p5] */\n\t"\
		"shlq	$3,%%rbx								\n\t		movq	%%r10,%%r11			\n\t"\
		"addq	%%rsi	,%%rax	/* add0 = &a[j1+p1] */	\n\t		addq	%%rbx	,%%r11	/* add1 = add0+p13 */\n\t"\
		"addq	%%rax	,%%rbx	/* add1 = add0+p8 */	\n\t"\
		"movq	%[__isrt2],%%rsi						\n\t"\
		"addq	$0x220,%%rcx	/* r22 */				\n\t"\
		"addq	$0x220,%%rdx	/* r32 */				\n\t"\
		"\n\t"\
		"addq	$0x030,%%rsi	/* cc1 */				/* cc3 = cc1 + 0x020: */\n\t"\
		"movaps	     (%%rcx),%%xmm4						\n\t		movaps	0x080(%%rcx),%%xmm12	\n\t"\
		"movaps	0x010(%%rcx),%%xmm5						\n\t		movaps	0x090(%%rcx),%%xmm13	\n\t"\
		"movaps	     (%%rsi),%%xmm2						\n\t		movaps	0x020(%%rsi),%%xmm11	\n\t"\
		"movaps	0x010(%%rsi),%%xmm3						\n\t		movaps	0x030(%%rsi),%%xmm10	\n\t"\
		"movaps	%%xmm4,%%xmm6							\n\t		movaps	%%xmm12,%%xmm14		\n\t"\
		"movaps	%%xmm5,%%xmm7							\n\t		movaps	%%xmm13,%%xmm15		\n\t"\
		"\n\t"\
		"mulpd	%%xmm2,%%xmm4							\n\t		mulpd	%%xmm10,%%xmm12		\n\t"\
		"mulpd	%%xmm2,%%xmm5							\n\t		mulpd	%%xmm10,%%xmm13		\n\t"\
		"mulpd	%%xmm3,%%xmm6							\n\t		mulpd	%%xmm11,%%xmm14		\n\t"\
		"movaps	     (%%rdx),%%xmm0						\n\t		movaps	0x080(%%rdx),%%xmm8	\n\t"\
		"mulpd	%%xmm3,%%xmm7							\n\t		mulpd	%%xmm11,%%xmm15		\n\t"\
		"movaps	0x010(%%rdx),%%xmm1						\n\t		movaps	0x090(%%rdx),%%xmm9	/* t3B */\n\t"\
		"/* cc3 */										/* cc1 */\n\t"\
		"/*movaps	0x020(%%rsi),%%xmm2					\n\t		movaps	     (%%rsi),%%xmm10	*/\n\t"\
		"/*movaps	0x030(%%rsi),%%xmm3					\n\t		movaps	0x010(%%rsi),%%xmm11	*/\n\t"\
		"subpd	%%xmm6,%%xmm5							\n\t		subpd	%%xmm14,%%xmm13		\n\t"\
		"movaps	%%xmm0,%%xmm6							\n\t		movaps	%%xmm8 ,%%xmm14		\n\t"\
		"addpd	%%xmm7,%%xmm4							\n\t		addpd	%%xmm15,%%xmm12		\n\t"\
		"movaps	%%xmm1,%%xmm7							\n\t		movaps	%%xmm9 ,%%xmm15		\n\t"\
		"\n\t"\
		"mulpd	%%xmm11,%%xmm0							\n\t		mulpd	%%xmm2 ,%%xmm8		\n\t"\
		"mulpd	%%xmm11,%%xmm1							\n\t		mulpd	%%xmm2 ,%%xmm9		\n\t"\
		"mulpd	%%xmm10,%%xmm6							\n\t		mulpd	%%xmm3 ,%%xmm14		\n\t"\
		"mulpd	%%xmm10,%%xmm7							\n\t		mulpd	%%xmm3 ,%%xmm15		\n\t"\
		"subpd	%%xmm6,%%xmm1							\n\t		addpd	%%xmm14,%%xmm9		\n\t"\
		"addpd	%%xmm7,%%xmm0							\n\t		subpd	%%xmm15,%%xmm8		\n\t"\
		"\n\t"\
		"movaps	%%xmm5,%%xmm7							\n\t		movaps	%%xmm13,%%xmm15		\n\t"\
		"movaps	%%xmm4,%%xmm6							\n\t		movaps	%%xmm12,%%xmm14		\n\t"\
		"\n\t"\
		"addpd	%%xmm0,%%xmm4							\n\t		addpd	%%xmm8 ,%%xmm14		\n\t"\
		"addpd	%%xmm1,%%xmm5							\n\t		addpd	%%xmm9 ,%%xmm15		\n\t"\
		"subpd	%%xmm0,%%xmm6							\n\t		subpd	%%xmm8 ,%%xmm12		\n\t"\
		"subpd	%%xmm1,%%xmm7							\n\t		subpd	%%xmm9 ,%%xmm13		\n\t"\
		"\n\t"\
		"subq	$0x200,%%rcx	/* r02 */				\n\t"\
		"subq	$0x200,%%rdx	/* r12 */				\n\t"\
		"subq	$0x020,%%rsi	/* cc0 */				\n\t"\
		"movaps	     (%%rdx),%%xmm1						\n\t		movaps	0x080(%%rdx),%%xmm8	\n\t"\
		"movaps	0x010(%%rdx),%%xmm3						\n\t		movaps	0x090(%%rdx),%%xmm10	\n\t"\
		"movaps	0x010(%%rsi),%%xmm2						\n\t		movaps	0x010(%%rsi),%%xmm9	\n\t"\
		"movaps	%%xmm1,%%xmm0							\n\t		movaps	%%xmm8 ,%%xmm11		\n\t"\
		"mulpd	%%xmm2,%%xmm1							\n\t		mulpd	%%xmm9 ,%%xmm8		\n\t"\
		"mulpd	%%xmm3,%%xmm2							\n\t		mulpd	%%xmm10,%%xmm9		\n\t"\
		"mulpd	(%%rsi),%%xmm0							\n\t		mulpd	(%%rsi),%%xmm11		\n\t"\
		"mulpd	(%%rsi),%%xmm3							\n\t		mulpd	(%%rsi),%%xmm10		\n\t"\
		"\n\t"\
		"addpd	%%xmm0,%%xmm2							\n\t		subpd	%%xmm10,%%xmm8		\n\t"\
		"subpd	%%xmm1,%%xmm3							\n\t		addpd	%%xmm11,%%xmm9		\n\t"\
		"\n\t"\
		"movaps	     (%%rcx),%%xmm0						\n\t		movaps	0x080(%%rcx),%%xmm10	\n\t"\
		"movaps	0x010(%%rcx),%%xmm1						\n\t		movaps	0x090(%%rcx),%%xmm11	\n\t"\
		"subpd	%%xmm2,%%xmm0							\n\t		subpd	%%xmm8 ,%%xmm10		\n\t"\
		"subpd	%%xmm3,%%xmm1							\n\t		subpd	%%xmm9 ,%%xmm11		\n\t"\
		"mulpd	(%%r8),%%xmm2							\n\t		mulpd	(%%r8),%%xmm8		\n\t"\
		"mulpd	(%%r8),%%xmm3							\n\t		mulpd	(%%r8),%%xmm9		\n\t"\
		"addpd	%%xmm0,%%xmm2							\n\t		addpd	%%xmm10,%%xmm8		\n\t"\
		"addpd	%%xmm1,%%xmm3							\n\t		addpd	%%xmm11,%%xmm9		\n\t"\
		"\n\t"\
		"/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c01): *//* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c05): */\n\t"\
		"addq	$0x260,%%rsi	/* c01 */				/* c05 = c01 + 0x80: */\n\t"\
		"subpd	%%xmm4,%%xmm2							\n\t		subpd	%%xmm12,%%xmm10		\n\t"\
		"subpd	%%xmm7,%%xmm0							\n\t		subpd	%%xmm15,%%xmm8		\n\t"\
		"subpd	%%xmm5,%%xmm3							\n\t		subpd	%%xmm13,%%xmm11		\n\t"\
		"subpd	%%xmm6,%%xmm1							\n\t		subpd	%%xmm14,%%xmm9		\n\t"\
		"mulpd	(%%r8),%%xmm4							\n\t		mulpd	(%%r8),%%xmm12		\n\t"\
		"mulpd	(%%r8),%%xmm7							\n\t		mulpd	(%%r8),%%xmm15		\n\t"\
		"mulpd	(%%r8),%%xmm5							\n\t		mulpd	(%%r8),%%xmm13		\n\t"\
		"mulpd	(%%r8),%%xmm6							\n\t		mulpd	(%%r8),%%xmm14		\n\t"\
		"addpd	%%xmm2,%%xmm4							\n\t		addpd	%%xmm10,%%xmm12		\n\t"\
		"addpd	%%xmm0,%%xmm7							\n\t		addpd	%%xmm8 ,%%xmm15		\n\t"\
		"addpd	%%xmm3,%%xmm5							\n\t		addpd	%%xmm11,%%xmm13		\n\t"\
		"addpd	%%xmm1,%%xmm6							\n\t		addpd	%%xmm9 ,%%xmm14		\n\t"\
		"movaps	%%xmm2,     (%%rcx)						\n\t		movaps	%%xmm10,0x080(%%rcx)	\n\t"\
		"movaps	%%xmm0,0x010(%%rdx)						\n\t		movaps	%%xmm8 ,0x090(%%rdx)	\n\t"\
		"movaps	%%xmm3,0x010(%%rcx)						\n\t		movaps	%%xmm11,0x090(%%rcx)	\n\t"\
		"movaps	%%xmm6,     (%%rdx)						\n\t		movaps	%%xmm14,0x080(%%rdx)	\n\t"\
		"movaps	%%xmm4,%%xmm2							\n\t		movaps	%%xmm12,%%xmm10		\n\t"\
		"movaps	%%xmm7,%%xmm0							\n\t		movaps	%%xmm15,%%xmm8		\n\t"\
		"movaps	%%xmm5,%%xmm3							\n\t		movaps	%%xmm13,%%xmm11		\n\t"\
		"movaps	%%xmm1,%%xmm6							\n\t		movaps	%%xmm9 ,%%xmm14		\n\t"\
		"mulpd	     (%%rsi),%%xmm4						\n\t		mulpd	0x080(%%rsi),%%xmm12	\n\t"\
		"mulpd	0x020(%%rsi),%%xmm7						\n\t		mulpd	0x0a0(%%rsi),%%xmm15	\n\t"\
		"mulpd	     (%%rsi),%%xmm5						\n\t		mulpd	0x080(%%rsi),%%xmm13	\n\t"\
		"mulpd	0x020(%%rsi),%%xmm1						\n\t		mulpd	0x0a0(%%rsi),%%xmm9	\n\t"\
		"mulpd	0x010(%%rsi),%%xmm2						\n\t		mulpd	0x090(%%rsi),%%xmm10	\n\t"\
		"mulpd	0x030(%%rsi),%%xmm0						\n\t		mulpd	0x0b0(%%rsi),%%xmm8	\n\t"\
		"mulpd	0x010(%%rsi),%%xmm3						\n\t		mulpd	0x090(%%rsi),%%xmm11	\n\t"\
		"mulpd	0x030(%%rsi),%%xmm6						\n\t		mulpd	0x0b0(%%rsi),%%xmm14	\n\t"\
		"subpd	%%xmm2,%%xmm5							\n\t		subpd	%%xmm10,%%xmm13		\n\t"\
		"subpd	%%xmm0,%%xmm1							\n\t		subpd	%%xmm8 ,%%xmm9		\n\t"\
		"addpd	%%xmm3,%%xmm4							\n\t		addpd	%%xmm11,%%xmm12		\n\t"\
		"addpd	%%xmm6,%%xmm7							\n\t		addpd	%%xmm14,%%xmm15		\n\t"\
		"movaps	%%xmm5,0x10(%%rax)						\n\t		movaps	%%xmm13,0x10(%%r10)	\n\t"\
		"movaps	%%xmm1,0x10(%%rbx)						\n\t		movaps	%%xmm9 ,0x10(%%r11)	\n\t"\
		"movaps	%%xmm4,    (%%rax)						\n\t		movaps	%%xmm12,    (%%r10)	\n\t"\
		"movaps	%%xmm7,    (%%rbx)						\n\t		movaps	%%xmm15,    (%%r11)	\n\t"\
		"addq	%%rdi,%%rax								\n\t		addq	%%rdi,%%r10			\n\t"\
		"addq	%%rdi,%%rbx								\n\t		addq	%%rdi,%%r11			\n\t"\
		"addq	$0x40,%%rsi		/* c11 */				/* c15 = c11 + 0x80: */\n\t"\
		"movaps	     (%%rcx),%%xmm4						\n\t		movaps	0x080(%%rcx),%%xmm12	\n\t"\
		"movaps	0x010(%%rdx),%%xmm0						\n\t		movaps	0x090(%%rdx),%%xmm8	\n\t"\
		"movaps	0x010(%%rcx),%%xmm5						\n\t		movaps	0x090(%%rcx),%%xmm13	\n\t"\
		"movaps	     (%%rdx),%%xmm6						\n\t		movaps	0x080(%%rdx),%%xmm14	\n\t"\
		"movaps	%%xmm4,%%xmm2							\n\t		movaps	%%xmm12,%%xmm10		\n\t"\
		"movaps	%%xmm0,%%xmm1							\n\t		movaps	%%xmm8 ,%%xmm9		\n\t"\
		"movaps	%%xmm5,%%xmm3							\n\t		movaps	%%xmm13,%%xmm11		\n\t"\
		"movaps	%%xmm6,%%xmm7							\n\t		movaps	%%xmm14,%%xmm15		\n\t"\
		"mulpd	     (%%rsi),%%xmm4						\n\t		mulpd	0x080(%%rsi),%%xmm12	\n\t"\
		"mulpd	0x020(%%rsi),%%xmm0						\n\t		mulpd	0x0a0(%%rsi),%%xmm8	\n\t"\
		"mulpd	     (%%rsi),%%xmm5						\n\t		mulpd	0x080(%%rsi),%%xmm13	\n\t"\
		"mulpd	0x020(%%rsi),%%xmm6						\n\t		mulpd	0x0a0(%%rsi),%%xmm14	\n\t"\
		"mulpd	0x010(%%rsi),%%xmm2						\n\t		mulpd	0x090(%%rsi),%%xmm10	\n\t"\
		"mulpd	0x030(%%rsi),%%xmm1						\n\t		mulpd	0x0b0(%%rsi),%%xmm9	\n\t"\
		"mulpd	0x010(%%rsi),%%xmm3						\n\t		mulpd	0x090(%%rsi),%%xmm11	\n\t"\
		"mulpd	0x030(%%rsi),%%xmm7						\n\t		mulpd	0x0b0(%%rsi),%%xmm15	\n\t"\
		"subpd	%%xmm2,%%xmm5							\n\t		subpd	%%xmm10,%%xmm13		\n\t"\
		"subpd	%%xmm1,%%xmm6							\n\t		subpd	%%xmm9 ,%%xmm14		\n\t"\
		"addpd	%%xmm3,%%xmm4							\n\t		addpd	%%xmm11,%%xmm12		\n\t"\
		"addpd	%%xmm7,%%xmm0							\n\t		addpd	%%xmm15,%%xmm8		\n\t"\
		"movaps	%%xmm5,0x10(%%rax)						\n\t		movaps	%%xmm13,0x10(%%r10)	\n\t"\
		"movaps	%%xmm6,0x10(%%rbx)						\n\t		movaps	%%xmm14,0x10(%%r11)	\n\t"\
		"movaps	%%xmm4,    (%%rax)						\n\t		movaps	%%xmm12,    (%%r10)	\n\t"\
		"movaps	%%xmm0,    (%%rbx)						\n\t		movaps	%%xmm8 ,    (%%r11)	\n\t"\
		"\n\t"\
		"/*...Block 3: t04,t14,t24,t34*/				/*...Block 7: t0C,t1C,t2C,t3C*/\n\t"\
		"movq	%[__add0],%%rax							\n\t		movq	%%rax,%%r10			\n\t"\
		"movslq	%[__p02],%%rsi							\n\t		movslq	%[__p06],%%r9		\n\t"\
		"movslq	%[__p08],%%rbx							\n\t		shlq	$3,%%r9			\n\t"\
		"shlq	$3,%%rsi								\n\t		addq	%%r9 	,%%r10	/* add0 = &a[j1+p6] */\n\t"\
		"shlq	$3,%%rbx								\n\t		movq	%%r10,%%r11			\n\t"\
		"addq	%%rsi	,%%rax	/* add0 = &a[j1+p2] */	\n\t		addq	%%rbx	,%%r11	/* add1 = add0+p14 */\n\t"\
		"addq	%%rax	,%%rbx	/* add1 = add0+p8 */	\n\t"\
		"movq	%[__isrt2],%%rsi						\n\t"\
		"addq	$0x220,%%rcx	/* r24 */				\n\t"\
		"addq	$0x220,%%rdx	/* r34 */				\n\t"\
		"\n\t"\
		"addq	$0x010,%%rsi	/* cc0 */\n\t"\
		"movaps	     (%%rcx),%%xmm4						\n\t		movaps	0x080(%%rcx),%%xmm12	\n\t"\
		"movaps	0x010(%%rcx),%%xmm5						\n\t		movaps	0x090(%%rcx),%%xmm13	\n\t"\
		"movaps	     (%%rsi),%%xmm2		\n\t"\
		"movaps	0x010(%%rsi),%%xmm3		\n\t"\
		"movaps	%%xmm4,%%xmm6							\n\t		movaps	%%xmm12,%%xmm14		\n\t"\
		"movaps	%%xmm5,%%xmm7							\n\t		movaps	%%xmm13,%%xmm15		\n\t"\
		"\n\t"\
		"mulpd	%%xmm2,%%xmm4							\n\t		mulpd	%%xmm3,%%xmm12	\n\t"\
		"mulpd	%%xmm2,%%xmm5							\n\t		mulpd	%%xmm3,%%xmm13	\n\t"\
		"mulpd	%%xmm3,%%xmm6							\n\t		mulpd	%%xmm2,%%xmm14	\n\t"\
		"movaps	     (%%rdx),%%xmm0						\n\t		movaps	0x080(%%rdx),%%xmm8	\n\t"\
		"mulpd	%%xmm3,%%xmm7							\n\t		mulpd	%%xmm2,%%xmm15	\n\t"\
		"movaps	0x010(%%rdx),%%xmm1						\n\t		movaps	0x090(%%rdx),%%xmm9	\n\t"\
		"subpd	%%xmm6,%%xmm5							\n\t		subpd	%%xmm14,%%xmm13		\n\t"\
		"movaps	%%xmm0,%%xmm6							\n\t		movaps	%%xmm8,%%xmm14		\n\t"\
		"addpd	%%xmm7,%%xmm4							\n\t		addpd	%%xmm15,%%xmm12		\n\t"\
		"movaps	%%xmm1,%%xmm7							\n\t		movaps	%%xmm9,%%xmm15		\n\t"\
		"\n\t"\
		"mulpd	%%xmm3,%%xmm0							\n\t		mulpd	%%xmm2,%%xmm8	\n\t"\
		"mulpd	%%xmm3,%%xmm1							\n\t		mulpd	%%xmm2,%%xmm9	\n\t"\
		"mulpd	%%xmm2,%%xmm6							\n\t		mulpd	%%xmm3,%%xmm14	\n\t"\
		"mulpd	%%xmm2,%%xmm7							\n\t		mulpd	%%xmm3,%%xmm15	\n\t"\
		"subpd	%%xmm6,%%xmm1							\n\t		subpd	%%xmm14,%%xmm9		\n\t"\
		"addpd	%%xmm7,%%xmm0							\n\t		addpd	%%xmm15,%%xmm8		\n\t"\
		"\n\t"\
		"movaps	%%xmm5,%%xmm7							\n\t		movaps	%%xmm13,%%xmm15		\n\t"\
		"movaps	%%xmm4,%%xmm6							\n\t		movaps	%%xmm12,%%xmm14		\n\t"\
		"\n\t"\
		"addpd	%%xmm0,%%xmm4							\n\t		addpd	%%xmm8,%%xmm14		\n\t"\
		"addpd	%%xmm1,%%xmm5							\n\t		addpd	%%xmm9,%%xmm15		\n\t"\
		"subpd	%%xmm0,%%xmm6							\n\t		subpd	%%xmm8,%%xmm12		\n\t"\
		"subpd	%%xmm1,%%xmm7							\n\t		subpd	%%xmm9,%%xmm13		\n\t"\
		"\n\t"\
		"subq	$0x200,%%rcx				\n\t"\
		"subq	$0x200,%%rdx				\n\t"\
		"subq	$0x10,%%rsi	/* isrt2 */		\n\t"\
		"movaps	     (%%rdx),%%xmm2						\n\t		movaps	0x080(%%rdx),%%xmm8	\n\t"\
		"movaps	0x010(%%rdx),%%xmm3						\n\t		movaps	0x090(%%rdx),%%xmm9	\n\t"\
		"movaps	(%%rsi),%%xmm1	/* isrt2 */	\n\t"\
		"movaps	%%xmm3,%%xmm0							\n\t		movaps	%%xmm8 ,%%xmm10		\n\t"\
		"subpd	%%xmm2,%%xmm3							\n\t		subpd	%%xmm9 ,%%xmm8		\n\t"\
		"addpd	%%xmm0,%%xmm2							\n\t		addpd	%%xmm10,%%xmm9		\n\t"\
		"mulpd	%%xmm1,%%xmm2							\n\t		mulpd	%%xmm1,%%xmm8	\n\t"\
		"mulpd	%%xmm1,%%xmm3							\n\t		mulpd	%%xmm1,%%xmm9	\n\t"\
		"\n\t"\
		"movaps	     (%%rcx),%%xmm0						\n\t		movaps	0x080(%%rcx),%%xmm10	\n\t"\
		"movaps	0x010(%%rcx),%%xmm1						\n\t		movaps	0x090(%%rcx),%%xmm11	\n\t"\
		"subpd	%%xmm2,%%xmm0							\n\t		subpd	%%xmm8 ,%%xmm10		\n\t"\
		"subpd	%%xmm3,%%xmm1							\n\t		subpd	%%xmm9 ,%%xmm11		\n\t"\
		"mulpd	(%%r8),%%xmm2							\n\t		mulpd	(%%r8),%%xmm8		\n\t"\
		"mulpd	(%%r8),%%xmm3							\n\t		mulpd	(%%r8),%%xmm9		\n\t"\
		"addpd	%%xmm0,%%xmm2							\n\t		addpd	%%xmm10,%%xmm8		\n\t"\
		"addpd	%%xmm1,%%xmm3							\n\t		addpd	%%xmm11,%%xmm9		\n\t"\
		"\n\t"\
		"/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c02): */\n\t		/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c06): */\n\t"\
		"addq	$0x170,%%rsi	/* c02 */				/* c06 = c02 + 0x80: */\n\t"\
		"subpd	%%xmm4,%%xmm2							\n\t		subpd	%%xmm12,%%xmm10		\n\t"\
		"subpd	%%xmm7,%%xmm0							\n\t		subpd	%%xmm15,%%xmm8		\n\t"\
		"subpd	%%xmm5,%%xmm3							\n\t		subpd	%%xmm13,%%xmm11		\n\t"\
		"subpd	%%xmm6,%%xmm1							\n\t		subpd	%%xmm14,%%xmm9		\n\t"\
		"mulpd	(%%r8),%%xmm4							\n\t		mulpd	(%%r8),%%xmm12		\n\t"\
		"mulpd	(%%r8),%%xmm7							\n\t		mulpd	(%%r8),%%xmm15		\n\t"\
		"mulpd	(%%r8),%%xmm5							\n\t		mulpd	(%%r8),%%xmm13		\n\t"\
		"mulpd	(%%r8),%%xmm6							\n\t		mulpd	(%%r8),%%xmm14		\n\t"\
		"addpd	%%xmm2,%%xmm4							\n\t		addpd	%%xmm10,%%xmm12		\n\t"\
		"addpd	%%xmm0,%%xmm7							\n\t		addpd	%%xmm8 ,%%xmm15		\n\t"\
		"addpd	%%xmm3,%%xmm5							\n\t		addpd	%%xmm11,%%xmm13		\n\t"\
		"addpd	%%xmm1,%%xmm6							\n\t		addpd	%%xmm9 ,%%xmm14		\n\t"\
		"movaps	%%xmm2,     (%%rcx)						\n\t		movaps	%%xmm10,0x080(%%rcx)	\n\t"\
		"movaps	%%xmm0,0x010(%%rdx)						\n\t		movaps	%%xmm8 ,0x090(%%rdx)	\n\t"\
		"movaps	%%xmm3,0x010(%%rcx)						\n\t		movaps	%%xmm11,0x090(%%rcx)	\n\t"\
		"movaps	%%xmm6,     (%%rdx)						\n\t		movaps	%%xmm14,0x080(%%rdx)	\n\t"\
		"movaps	%%xmm4,%%xmm2							\n\t		movaps	%%xmm12,%%xmm10		\n\t"\
		"movaps	%%xmm7,%%xmm0							\n\t		movaps	%%xmm15,%%xmm8		\n\t"\
		"movaps	%%xmm5,%%xmm3							\n\t		movaps	%%xmm13,%%xmm11		\n\t"\
		"movaps	%%xmm1,%%xmm6							\n\t		movaps	%%xmm9 ,%%xmm14		\n\t"\
		"mulpd	     (%%rsi),%%xmm4						\n\t		mulpd	0x080(%%rsi),%%xmm12	\n\t"\
		"mulpd	0x020(%%rsi),%%xmm7						\n\t		mulpd	0x0a0(%%rsi),%%xmm15	\n\t"\
		"mulpd	     (%%rsi),%%xmm5						\n\t		mulpd	0x080(%%rsi),%%xmm13	\n\t"\
		"mulpd	0x020(%%rsi),%%xmm1						\n\t		mulpd	0x0a0(%%rsi),%%xmm9	\n\t"\
		"mulpd	0x010(%%rsi),%%xmm2						\n\t		mulpd	0x090(%%rsi),%%xmm10	\n\t"\
		"mulpd	0x030(%%rsi),%%xmm0						\n\t		mulpd	0x0b0(%%rsi),%%xmm8	\n\t"\
		"mulpd	0x010(%%rsi),%%xmm3						\n\t		mulpd	0x090(%%rsi),%%xmm11	\n\t"\
		"mulpd	0x030(%%rsi),%%xmm6						\n\t		mulpd	0x0b0(%%rsi),%%xmm14	\n\t"\
		"subpd	%%xmm2,%%xmm5							\n\t		subpd	%%xmm10,%%xmm13		\n\t"\
		"subpd	%%xmm0,%%xmm1							\n\t		subpd	%%xmm8 ,%%xmm9		\n\t"\
		"addpd	%%xmm3,%%xmm4							\n\t		addpd	%%xmm11,%%xmm12		\n\t"\
		"addpd	%%xmm6,%%xmm7							\n\t		addpd	%%xmm14,%%xmm15		\n\t"\
		"movaps	%%xmm5,0x10(%%rax)						\n\t		movaps	%%xmm13,0x10(%%r10)	\n\t"\
		"movaps	%%xmm1,0x10(%%rbx)						\n\t		movaps	%%xmm9 ,0x10(%%r11)	\n\t"\
		"movaps	%%xmm4,    (%%rax)						\n\t		movaps	%%xmm12,    (%%r10)	\n\t"\
		"movaps	%%xmm7,    (%%rbx)						\n\t		movaps	%%xmm15,    (%%r11)	\n\t"\
		"addq	%%rdi,%%rax								\n\t		addq	%%rdi,%%r10			\n\t"\
		"addq	%%rdi,%%rbx								\n\t		addq	%%rdi,%%r11			\n\t"\
		"addq	$0x40,%%rsi		/* c12 */				/* c16 = c12 + 0x80: */\n\t"\
		"movaps	     (%%rcx),%%xmm4						\n\t		movaps	0x080(%%rcx),%%xmm12	\n\t"\
		"movaps	0x010(%%rdx),%%xmm0						\n\t		movaps	0x090(%%rdx),%%xmm8	\n\t"\
		"movaps	0x010(%%rcx),%%xmm5						\n\t		movaps	0x090(%%rcx),%%xmm13	\n\t"\
		"movaps	     (%%rdx),%%xmm6						\n\t		movaps	0x080(%%rdx),%%xmm14	\n\t"\
		"movaps	%%xmm4,%%xmm2							\n\t		movaps	%%xmm12,%%xmm10		\n\t"\
		"movaps	%%xmm0,%%xmm1							\n\t		movaps	%%xmm8 ,%%xmm9		\n\t"\
		"movaps	%%xmm5,%%xmm3							\n\t		movaps	%%xmm13,%%xmm11		\n\t"\
		"movaps	%%xmm6,%%xmm7							\n\t		movaps	%%xmm14,%%xmm15		\n\t"\
		"mulpd	     (%%rsi),%%xmm4						\n\t		mulpd	0x080(%%rsi),%%xmm12	\n\t"\
		"mulpd	0x020(%%rsi),%%xmm0						\n\t		mulpd	0x0a0(%%rsi),%%xmm8	\n\t"\
		"mulpd	     (%%rsi),%%xmm5						\n\t		mulpd	0x080(%%rsi),%%xmm13	\n\t"\
		"mulpd	0x020(%%rsi),%%xmm6						\n\t		mulpd	0x0a0(%%rsi),%%xmm14	\n\t"\
		"mulpd	0x010(%%rsi),%%xmm2						\n\t		mulpd	0x090(%%rsi),%%xmm10	\n\t"\
		"mulpd	0x030(%%rsi),%%xmm1						\n\t		mulpd	0x0b0(%%rsi),%%xmm9	\n\t"\
		"mulpd	0x010(%%rsi),%%xmm3						\n\t		mulpd	0x090(%%rsi),%%xmm11	\n\t"\
		"mulpd	0x030(%%rsi),%%xmm7						\n\t		mulpd	0x0b0(%%rsi),%%xmm15	\n\t"\
		"subpd	%%xmm2,%%xmm5							\n\t		subpd	%%xmm10,%%xmm13		\n\t"\
		"subpd	%%xmm1,%%xmm6							\n\t		subpd	%%xmm9 ,%%xmm14		\n\t"\
		"addpd	%%xmm3,%%xmm4							\n\t		addpd	%%xmm11,%%xmm12		\n\t"\
		"addpd	%%xmm7,%%xmm0							\n\t		addpd	%%xmm15,%%xmm8		\n\t"\
		"movaps	%%xmm5,0x10(%%rax)						\n\t		movaps	%%xmm13,0x10(%%r10)	\n\t"\
		"movaps	%%xmm6,0x10(%%rbx)						\n\t		movaps	%%xmm14,0x10(%%r11)	\n\t"\
		"movaps	%%xmm4,    (%%rax)						\n\t		movaps	%%xmm12,    (%%r10)	\n\t"\
		"movaps	%%xmm0,    (%%rbx)						\n\t		movaps	%%xmm8 ,    (%%r11)	\n\t"\
		"\n\t"\
		"/*...Block 4: t06,t16,t26,t36*/				\n\t		/*...Block 8: t0E,t1E,t2E,t3E*/\n\t"\
		"movq	%[__add0],%%rax							\n\t		movq	%%rax,%%r10			\n\t"\
		"movslq	%[__p03],%%rsi							\n\t		movslq	%[__p07],%%r9		\n\t"\
		"movslq	%[__p08],%%rbx							\n\t		shlq	$3,%%r9			\n\t"\
		"shlq	$3,%%rsi								\n\t		addq	%%r9 	,%%r10	/* add0 = &a[j1+p6] */\n\t"\
		"shlq	$3,%%rbx								\n\t		movq	%%r10,%%r11			\n\t"\
		"addq	%%rsi	,%%rax	/* add0 = &a[j1+p3] */	\n\t		addq	%%rbx	,%%r11	/* add1 = add0+p14 */\n\t"\
		"addq	%%rax	,%%rbx	/* add1 = add0+p8 */	\n\t"\
		"movq	%[__isrt2],%%rsi						\n\t"\
		"addq	$0x220,%%rcx	/* r26 */				\n\t"\
		"addq	$0x220,%%rdx	/* r36 */				\n\t"\
		"addq	$0x030,%%rsi	/* cc1 */				\n\t"\
		"\n\t"\
		"movaps	     (%%rcx),%%xmm4						\n\t		movaps	0x080(%%rcx),%%xmm12	\n\t"\
		"movaps	0x010(%%rcx),%%xmm5						\n\t		movaps	0x090(%%rcx),%%xmm13	\n\t"\
		"movaps	     (%%rsi),%%xmm2	/* c32_1 */			\n\t		movaps	0x020(%%rsi),%%xmm10	/* c32_3 */\n\t"\
		"movaps	0x010(%%rsi),%%xmm3	/* s32_1 */			\n\t		movaps	0x030(%%rsi),%%xmm11	/* s32_3 */\n\t"\
		"movaps	%%xmm4,%%xmm6							\n\t		movaps	%%xmm12,%%xmm14		\n\t"\
		"movaps	%%xmm5,%%xmm7							\n\t		movaps	%%xmm13,%%xmm15		\n\t"\
		"\n\t"\
		"mulpd	%%xmm10,%%xmm4							\n\t		mulpd	%%xmm3,%%xmm12		\n\t"\
		"mulpd	%%xmm10,%%xmm5							\n\t		mulpd	%%xmm3,%%xmm13		\n\t"\
		"mulpd	%%xmm11,%%xmm6							\n\t		mulpd	%%xmm2,%%xmm14		\n\t"\
		"movaps	     (%%rdx),%%xmm0						\n\t		movaps	0x080(%%rdx),%%xmm8	\n\t"\
		"mulpd	%%xmm11,%%xmm7							\n\t		mulpd	%%xmm2,%%xmm15		\n\t"\
		"movaps	0x010(%%rdx),%%xmm1						\n\t		movaps	0x090(%%rdx),%%xmm9	\n\t"\
		"subpd	%%xmm6,%%xmm5							\n\t		subpd	%%xmm14,%%xmm13		\n\t"\
		"movaps	%%xmm0,%%xmm6							\n\t		movaps	%%xmm8 ,%%xmm14		\n\t"\
		"addpd	%%xmm7,%%xmm4							\n\t		addpd	%%xmm15,%%xmm12		\n\t"\
		"movaps	%%xmm1,%%xmm7							\n\t		movaps	%%xmm9 ,%%xmm15		\n\t"\
		"\n\t"\
		"mulpd	%%xmm3,%%xmm0							\n\t		mulpd	%%xmm11,%%xmm8		\n\t"\
		"mulpd	%%xmm3,%%xmm1							\n\t		mulpd	%%xmm11,%%xmm9		\n\t"\
		"mulpd	%%xmm2,%%xmm6							\n\t		mulpd	%%xmm10,%%xmm14		\n\t"\
		"mulpd	%%xmm2,%%xmm7							\n\t		mulpd	%%xmm10,%%xmm15		\n\t"\
		"addpd	%%xmm6,%%xmm1							\n\t		subpd	%%xmm14,%%xmm9		\n\t"\
		"subpd	%%xmm7,%%xmm0							\n\t		addpd	%%xmm15,%%xmm8		\n\t"\
		"\n\t"\
		"movaps	%%xmm5,%%xmm7							\n\t		movaps	%%xmm13,%%xmm15		\n\t"\
		"movaps	%%xmm4,%%xmm6							\n\t		movaps	%%xmm12,%%xmm14		\n\t"\
		"\n\t"\
		"addpd	%%xmm0,%%xmm6							\n\t		subpd	%%xmm8 ,%%xmm12		\n\t"\
		"addpd	%%xmm1,%%xmm7							\n\t		subpd	%%xmm9 ,%%xmm13		\n\t"\
		"subpd	%%xmm0,%%xmm4							\n\t		addpd	%%xmm8 ,%%xmm14		\n\t"\
		"subpd	%%xmm1,%%xmm5							\n\t		addpd	%%xmm9 ,%%xmm15		\n\t"\
		"\n\t"\
		"subq	$0x200,%%rcx			\n\t"\
		"subq	$0x200,%%rdx			\n\t"\
		"subq	$0x20,%%rsi	/* cc0 */	\n\t"\
		"movaps	     (%%rdx),%%xmm2						\n\t		movaps	0x080(%%rdx),%%xmm11	\n\t"\
		"movaps	0x010(%%rdx),%%xmm0						\n\t		movaps	0x090(%%rdx),%%xmm9	\n\t"\
		"movaps	0x010(%%rsi),%%xmm3						\n\t		movaps	0x010(%%rsi),%%xmm10	\n\t"\
		"movaps	%%xmm2,%%xmm1							\n\t		movaps	%%xmm11,%%xmm8		\n\t"\
		"mulpd	%%xmm3,%%xmm2							\n\t		mulpd	%%xmm10,%%xmm11		\n\t"\
		"mulpd	%%xmm0,%%xmm3							\n\t		mulpd	%%xmm9 ,%%xmm10		\n\t"\
		"mulpd	(%%rsi),%%xmm1							\n\t		mulpd	(%%rsi),%%xmm8		\n\t"\
		"mulpd	(%%rsi),%%xmm0							\n\t		mulpd	(%%rsi),%%xmm9		\n\t"\
		"\n\t"\
		"addpd	%%xmm0,%%xmm2							\n\t		subpd	%%xmm10,%%xmm8		\n\t"\
		"subpd	%%xmm1,%%xmm3							\n\t		addpd	%%xmm11,%%xmm9		\n\t"\
		"\n\t"\
		"movaps	     (%%rcx),%%xmm0						\n\t		movaps	0x080(%%rcx),%%xmm10	\n\t"\
		"movaps	0x010(%%rcx),%%xmm1						\n\t		movaps	0x090(%%rcx),%%xmm11	\n\t"\
		"subpd	%%xmm2,%%xmm0							\n\t		subpd	%%xmm8 ,%%xmm10		\n\t"\
		"subpd	%%xmm3,%%xmm1							\n\t		subpd	%%xmm9 ,%%xmm11		\n\t"\
		"mulpd	(%%r8),%%xmm2							\n\t		mulpd	(%%r8),%%xmm8		\n\t"\
		"mulpd	(%%r8),%%xmm3							\n\t		mulpd	(%%r8),%%xmm9		\n\t"\
		"addpd	%%xmm0,%%xmm2							\n\t		addpd	%%xmm10,%%xmm8		\n\t"\
		"addpd	%%xmm1,%%xmm3							\n\t		addpd	%%xmm11,%%xmm9		\n\t"\
		"\n\t"\
		"/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c03): *//* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c07): */\n\t"\
		"addq	$0x360,%%rsi	/* c03 */				/* c07 = c03 + 0x80: */\n\t"\
		"subpd	%%xmm4,%%xmm2							\n\t		subpd	%%xmm12,%%xmm10		\n\t"\
		"subpd	%%xmm7,%%xmm0							\n\t		subpd	%%xmm15,%%xmm8		\n\t"\
		"subpd	%%xmm5,%%xmm3							\n\t		subpd	%%xmm13,%%xmm11		\n\t"\
		"subpd	%%xmm6,%%xmm1							\n\t		subpd	%%xmm14,%%xmm9		\n\t"\
		"mulpd	(%%r8),%%xmm4							\n\t		mulpd	(%%r8),%%xmm12		\n\t"\
		"mulpd	(%%r8),%%xmm7							\n\t		mulpd	(%%r8),%%xmm15		\n\t"\
		"mulpd	(%%r8),%%xmm5							\n\t		mulpd	(%%r8),%%xmm13		\n\t"\
		"mulpd	(%%r8),%%xmm6							\n\t		mulpd	(%%r8),%%xmm14		\n\t"\
		"addpd	%%xmm2,%%xmm4							\n\t		addpd	%%xmm10,%%xmm12		\n\t"\
		"addpd	%%xmm0,%%xmm7							\n\t		addpd	%%xmm8 ,%%xmm15		\n\t"\
		"addpd	%%xmm3,%%xmm5							\n\t		addpd	%%xmm11,%%xmm13		\n\t"\
		"addpd	%%xmm1,%%xmm6							\n\t		addpd	%%xmm9 ,%%xmm14		\n\t"\
		"movaps	%%xmm2,     (%%rcx)						\n\t		movaps	%%xmm10,0x080(%%rcx)	\n\t"\
		"movaps	%%xmm0,0x010(%%rdx)						\n\t		movaps	%%xmm8 ,0x090(%%rdx)	\n\t"\
		"movaps	%%xmm3,0x010(%%rcx)						\n\t		movaps	%%xmm11,0x090(%%rcx)	\n\t"\
		"movaps	%%xmm6,     (%%rdx)						\n\t		movaps	%%xmm14,0x080(%%rdx)	\n\t"\
		"movaps	%%xmm4,%%xmm2							\n\t		movaps	%%xmm12,%%xmm10		\n\t"\
		"movaps	%%xmm7,%%xmm0							\n\t		movaps	%%xmm15,%%xmm8		\n\t"\
		"movaps	%%xmm5,%%xmm3							\n\t		movaps	%%xmm13,%%xmm11		\n\t"\
		"movaps	%%xmm1,%%xmm6							\n\t		movaps	%%xmm9 ,%%xmm14		\n\t"\
		"mulpd	     (%%rsi),%%xmm4						\n\t		mulpd	0x080(%%rsi),%%xmm12	\n\t"\
		"mulpd	0x020(%%rsi),%%xmm7						\n\t		mulpd	0x0a0(%%rsi),%%xmm15	\n\t"\
		"mulpd	     (%%rsi),%%xmm5						\n\t		mulpd	0x080(%%rsi),%%xmm13	\n\t"\
		"mulpd	0x020(%%rsi),%%xmm1						\n\t		mulpd	0x0a0(%%rsi),%%xmm9	\n\t"\
		"mulpd	0x010(%%rsi),%%xmm2						\n\t		mulpd	0x090(%%rsi),%%xmm10	\n\t"\
		"mulpd	0x030(%%rsi),%%xmm0						\n\t		mulpd	0x0b0(%%rsi),%%xmm8	\n\t"\
		"mulpd	0x010(%%rsi),%%xmm3						\n\t		mulpd	0x090(%%rsi),%%xmm11	\n\t"\
		"mulpd	0x030(%%rsi),%%xmm6						\n\t		mulpd	0x0b0(%%rsi),%%xmm14	\n\t"\
		"subpd	%%xmm2,%%xmm5							\n\t		subpd	%%xmm10,%%xmm13		\n\t"\
		"subpd	%%xmm0,%%xmm1							\n\t		subpd	%%xmm8 ,%%xmm9		\n\t"\
		"addpd	%%xmm3,%%xmm4							\n\t		addpd	%%xmm11,%%xmm12		\n\t"\
		"addpd	%%xmm6,%%xmm7							\n\t		addpd	%%xmm14,%%xmm15		\n\t"\
		"movaps	%%xmm5,0x10(%%rax)						\n\t		movaps	%%xmm13,0x10(%%r10)	\n\t"\
		"movaps	%%xmm1,0x10(%%rbx)						\n\t		movaps	%%xmm9 ,0x10(%%r11)	\n\t"\
		"movaps	%%xmm4,    (%%rax)						\n\t		movaps	%%xmm12,    (%%r10)	\n\t"\
		"movaps	%%xmm7,    (%%rbx)						\n\t		movaps	%%xmm15,    (%%r11)	\n\t"\
		"addq	%%rdi,%%rax								\n\t		addq	%%rdi,%%r10			\n\t"\
		"addq	%%rdi,%%rbx								\n\t		addq	%%rdi,%%r11			\n\t"\
		"addq	$0x40,%%rsi	/* c13 */					\n\t		/* c17 = c13 + 0x80: */\n\t"\
		"movaps	     (%%rcx),%%xmm4						\n\t		movaps	0x080(%%rcx),%%xmm12	\n\t"\
		"movaps	0x010(%%rdx),%%xmm0						\n\t		movaps	0x090(%%rdx),%%xmm8	\n\t"\
		"movaps	0x010(%%rcx),%%xmm5						\n\t		movaps	0x090(%%rcx),%%xmm13	\n\t"\
		"movaps	     (%%rdx),%%xmm6						\n\t		movaps	0x080(%%rdx),%%xmm14	\n\t"\
		"movaps	%%xmm4,%%xmm2							\n\t		movaps	%%xmm12,%%xmm10		\n\t"\
		"movaps	%%xmm0,%%xmm1							\n\t		movaps	%%xmm8,%%xmm9		\n\t"\
		"movaps	%%xmm5,%%xmm3							\n\t		movaps	%%xmm13,%%xmm11		\n\t"\
		"movaps	%%xmm6,%%xmm7							\n\t		movaps	%%xmm14,%%xmm15		\n\t"\
		"mulpd	     (%%rsi),%%xmm4						\n\t		mulpd	0x080(%%rsi),%%xmm12	\n\t"\
		"mulpd	0x020(%%rsi),%%xmm0						\n\t		mulpd	0x0a0(%%rsi),%%xmm8	\n\t"\
		"mulpd	     (%%rsi),%%xmm5						\n\t		mulpd	0x080(%%rsi),%%xmm13	\n\t"\
		"mulpd	0x020(%%rsi),%%xmm6						\n\t		mulpd	0x0a0(%%rsi),%%xmm14	\n\t"\
		"mulpd	0x010(%%rsi),%%xmm2						\n\t		mulpd	0x090(%%rsi),%%xmm10	\n\t"\
		"mulpd	0x030(%%rsi),%%xmm1						\n\t		mulpd	0x0b0(%%rsi),%%xmm9	\n\t"\
		"mulpd	0x010(%%rsi),%%xmm3						\n\t		mulpd	0x090(%%rsi),%%xmm11	\n\t"\
		"mulpd	0x030(%%rsi),%%xmm7						\n\t		mulpd	0x0b0(%%rsi),%%xmm15	\n\t"\
		"subpd	%%xmm2,%%xmm5							\n\t		subpd	%%xmm10,%%xmm13		\n\t"\
		"subpd	%%xmm1,%%xmm6							\n\t		subpd	%%xmm9 ,%%xmm14		\n\t"\
		"addpd	%%xmm3,%%xmm4							\n\t		addpd	%%xmm11,%%xmm12		\n\t"\
		"addpd	%%xmm7,%%xmm0							\n\t		addpd	%%xmm15,%%xmm8		\n\t"\
		"movaps	%%xmm5,0x10(%%rax)						\n\t		movaps	%%xmm13,0x10(%%r10)	\n\t"\
		"movaps	%%xmm6,0x10(%%rbx)						\n\t		movaps	%%xmm14,0x10(%%r11)	\n\t"\
		"movaps	%%xmm4,    (%%rax)						\n\t		movaps	%%xmm12,    (%%r10)	\n\t"\
		"movaps	%%xmm0,    (%%rbx)						\n\t		movaps	%%xmm8 ,    (%%r11)	\n\t"\
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
		 ,[__r10] "m" (Xr10)\
		 ,[__r20] "m" (Xr20)\
		 ,[__r30] "m" (Xr30)\
		 ,[__isrt2] "m" (Xisrt2)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r8","r9","r10","r11","r12","r13","r14","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
	);\
	}

  #endif // USE_64BIT_ASM_STYLE ?

#endif	// AVX / SSE2 toggle

#endif	/* radix32_dif_dit_pass_gcc_h_included */


	// If debugging register contents from source file, stick this code into the requisite inline ASM:
	#define FOO(Xr00)\
	{\
	__asm__ volatile (\
	"/******* DEBUG DUMP OF REGISTER CONTENTS ********/\n\t"\
		"movq	%[__r00],%%rsi	\n\t"\
		"vmovaps	%%ymm0 ,     (%%rsi)	\n\t"\
		"vmovaps	%%ymm1 ,0x020(%%rsi)	\n\t"\
		"vmovaps	%%ymm2 ,0x040(%%rsi)	\n\t"\
		"vmovaps	%%ymm3 ,0x060(%%rsi)	\n\t"\
		"vmovaps	%%ymm4 ,0x080(%%rsi)	\n\t"\
		"vmovaps	%%ymm5 ,0x0a0(%%rsi)	\n\t"\
		"vmovaps	%%ymm6 ,0x0c0(%%rsi)	\n\t"\
		"vmovaps	%%ymm7 ,0x0e0(%%rsi)	\n\t"\
		"vmovaps	%%ymm8, 0x100(%%rsi)	\n\t"\
		"vmovaps	%%ymm9, 0x120(%%rsi)	\n\t"\
		"vmovaps	%%ymm10,0x140(%%rsi)	\n\t"\
		"vmovaps	%%ymm11,0x160(%%rsi)	\n\t"\
		"vmovaps	%%ymm12,0x180(%%rsi)	\n\t"\
		"vmovaps	%%ymm13,0x1a0(%%rsi)	\n\t"\
		"vmovaps	%%ymm14,0x1c0(%%rsi)	\n\t"\
		"vmovaps	%%ymm15,0x1e0(%%rsi)	\n\t"\
	"/******* DEBUG DUMP OF REGISTER CONTENTS ********/\n\t"\
		"movq	%[__r00],%%rsi	\n\t"\
		"movaps	%%xmm0 ,     (%%rsi)	\n\t"\
		"movaps	%%xmm1 ,0x010(%%rsi)	\n\t"\
		"movaps	%%xmm2 ,0x020(%%rsi)	\n\t"\
		"movaps	%%xmm3 ,0x030(%%rsi)	\n\t"\
		"movaps	%%xmm4 ,0x040(%%rsi)	\n\t"\
		"movaps	%%xmm5 ,0x050(%%rsi)	\n\t"\
		"movaps	%%xmm6 ,0x060(%%rsi)	\n\t"\
		"movaps	%%xmm7 ,0x070(%%rsi)	\n\t"\
		"movaps	%%xmm8, 0x080(%%rsi)	\n\t"\
		"movaps	%%xmm9, 0x090(%%rsi)	\n\t"\
		"movaps	%%xmm10,0x0a0(%%rsi)	\n\t"\
		"movaps	%%xmm11,0x0b0(%%rsi)	\n\t"\
		"movaps	%%xmm12,0x0c0(%%rsi)	\n\t"\
		"movaps	%%xmm13,0x0d0(%%rsi)	\n\t"\
		"movaps	%%xmm14,0x0e0(%%rsi)	\n\t"\
		"movaps	%%xmm15,0x0f0(%%rsi)	\n\t"\
		:					/* outputs: none */\
		: [__r00] "m" (Xr00)	/* All inputs from memory addresses here */\
		: "cc","memory","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
	);\
	}

