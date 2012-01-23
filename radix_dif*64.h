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
#ifndef radix32_dif_dit_pass_gcc_h_included
#define radix32_dif_dit_pass_gcc_h_included

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
		"movslq	%[__p10],%%rcx						\n\t"\
		"movslq	%[__p18],%%rdx						\n\t"\
		"shlq	$3,%%rbx							\n\t	movslq	%[__p04],%%r8	\n\t"\
		"shlq	$3,%%rcx							\n\t	shlq	$3,%%r8	\n\t"\
		"shlq	$3,%%rdx							\n\t	movq	%%r8,%%r10	\n\t	addq	%%rax,%%r10	\n\t"\
		"addq	%%rax,%%rbx							\n\t	movq	%%r8,%%r11	\n\t	addq	%%rbx,%%r11	\n\t"\
		"addq	%%rax,%%rcx							\n\t	movq	%%r8,%%r12	\n\t	addq	%%rcx,%%r12	\n\t"\
		"addq	%%rax,%%rdx							\n\t	movq	%%r8,%%r13	\n\t	addq	%%rdx,%%r13	\n\t"\
		"/* SSE2_RADIX4_DIF_4TWIDDLE_A(r00,c10) */	\n\t	/* SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO_B(r08,c04) */\n\t"\
		"movq	%[__r00],%%rsi						\n\t	movaps	    (%%r10),%%xmm8	\n\t"\
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
		"addpd	%%xmm6,%%xmm6						\n\t	subpd	%%xmm15,%%xmm9	\n\t"\
		"addpd	%%xmm5,%%xmm5						\n\t	movaps	%%xmm8 ,0xc0(%%rsi)	\n\t"\
		"addpd	%%xmm7,%%xmm7						\n\t	subpd	%%xmm13,%%xmm10	\n\t"\
		"addpd	%%xmm4,%%xmm4						\n\t/*	movaps	%%xmm9 ,0xd0(%%rsi)	*/\n\t"\
		"addpd	%%xmm0,%%xmm6						\n\t	subpd	%%xmm12,%%xmm11	\n\t"\
		"addpd	%%xmm2,%%xmm5						\n\t	addpd	%%xmm14,%%xmm14	\n\t"\
		"addpd	%%xmm1,%%xmm7						\n\t	addpd	%%xmm13,%%xmm13	\n\t"\
		"addpd	%%xmm3,%%xmm4						\n\t	addpd	%%xmm15,%%xmm15	\n\t"\
		"/*movaps	%%xmm6,    (%%rsi)	*/			\n\t	addpd	%%xmm12,%%xmm12	\n\t"\
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
		"													mulpd	%%xmm8 ,%%xmm15	\n\t"\
		"												/*	movaps	%%xmm10,0xa0(%%rsi)	*/\n\t"\
		"												/*	movaps	%%xmm13,0xe0(%%rsi)	*/\n\t"\
		"												/*	movaps	%%xmm14,0xb0(%%rsi)	*/\n\t"\
		"												/*	movaps	%%xmm15,0xf0(%%rsi)	*/\n\t"\
		"/* SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r00) */	\n\t"\
		"movaps	0x80(%%rsi),%%xmm11		\n\t"\
		"movaps	0x90(%%rsi),%%xmm12		\n\t"\
		"movaps	0xc0(%%rsi),%%xmm8		\n\t"\
		"subpd	%%xmm11,%%xmm6						\n\t	subpd	%%xmm10,%%xmm2	\n\t"\
		"subpd	%%xmm9 ,%%xmm0						\n\t	subpd	%%xmm15,%%xmm5	\n\t"\
		"subpd	%%xmm12,%%xmm7						\n\t	subpd	%%xmm14,%%xmm4	\n\t"\
		"subpd	%%xmm8 ,%%xmm1						\n\t	subpd	%%xmm13,%%xmm3	\n\t"\
		"addpd	%%xmm11,%%xmm11						\n\t	addpd	%%xmm10,%%xmm10	\n\t"\
		"addpd	%%xmm9 ,%%xmm9						\n\t	addpd	%%xmm15,%%xmm15	\n\t"\
		"addpd	%%xmm12,%%xmm12						\n\t	addpd	%%xmm14,%%xmm14	\n\t"\
		"addpd	%%xmm8 ,%%xmm8						\n\t	addpd	%%xmm13,%%xmm13	\n\t"\
		"addpd	%%xmm6 ,%%xmm11						\n\t	addpd	%%xmm2 ,%%xmm10	\n\t"\
		"addpd	%%xmm0 ,%%xmm9						\n\t	addpd	%%xmm5 ,%%xmm15	\n\t"\
		"addpd	%%xmm7 ,%%xmm12						\n\t	addpd	%%xmm4 ,%%xmm14	\n\t"\
		"addpd	%%xmm1 ,%%xmm8						\n\t	addpd	%%xmm3 ,%%xmm13	\n\t"\
		"movaps	%%xmm6 ,0x80(%%rsi)					\n\t	movaps	%%xmm2 ,0xa0(%%rsi)	\n\t"\
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
		"movslq	%[__p02],%%rdi						\n\t	movslq	%[__p04],%%r8	\n\t"\
		"shlq	$3,%%rdi							\n\t	shlq	$3,%%r8	\n\t"\
		"addq	%%rdi,%%rax	/* &a[j1+p2] */			\n\t	movq	%%r8,%%r10	\n\t	addq	%%rax,%%r10	\n\t"\
		"addq	%%rdi,%%rbx							\n\t	movq	%%r8,%%r11	\n\t	addq	%%rbx,%%r11	\n\t"\
		"addq	%%rdi,%%rcx							\n\t	movq	%%r8,%%r12	\n\t	addq	%%rcx,%%r12	\n\t"\
		"addq	%%rdi,%%rdx							\n\t	movq	%%r8,%%r13	\n\t	addq	%%rdx,%%r13	\n\t"\
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
		"movaps	%%xmm2,0x020(%%rsi)					\n\t	addpd	%%xmm14,%%xmm14	\n\t"\
		"movaps	%%xmm1,0x050(%%rsi)					\n\t	addpd	%%xmm13,%%xmm13	\n\t"\
		"movaps	%%xmm3,0x070(%%rsi)					\n\t	addpd	%%xmm15,%%xmm15	\n\t"\
		"addpd	%%xmm6,%%xmm6						\n\t	addpd	%%xmm12,%%xmm12	\n\t"\
		"addpd	%%xmm5,%%xmm5						\n\t	addpd	%%xmm8 ,%%xmm14	\n\t"\
		"addpd	%%xmm7,%%xmm7						\n\t	addpd	%%xmm10,%%xmm13	\n\t"\
		"addpd	%%xmm4,%%xmm4						\n\t	addpd	%%xmm9 ,%%xmm15	\n\t"\
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
		"													mulpd	%%xmm8 ,%%xmm15	\n\t"\
		"												/*	movaps	%%xmm10,0xa0(%%rsi)	*/\n\t"\
		"												/*	movaps	%%xmm13,0xe0(%%rsi)	*/\n\t"\
		"												/*	movaps	%%xmm14,0xb0(%%rsi)	*/\n\t"\
		"												/*	movaps	%%xmm15,0xf0(%%rsi)	*/\n\t"\
		"/* SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r10) */	\n\t"\
		"movaps	0x80(%%rsi),%%xmm11		\n\t"\
		"movaps	0x90(%%rsi),%%xmm12		\n\t"\
		"movaps	0xc0(%%rsi),%%xmm8		\n\t"\
		"subpd	%%xmm11,%%xmm6						\n\t	subpd	%%xmm10,%%xmm2	\n\t"\
		"subpd	%%xmm9 ,%%xmm0						\n\t	subpd	%%xmm15,%%xmm5	\n\t"\
		"subpd	%%xmm12,%%xmm7						\n\t	subpd	%%xmm14,%%xmm4	\n\t"\
		"subpd	%%xmm8 ,%%xmm1						\n\t	subpd	%%xmm13,%%xmm3	\n\t"\
		"addpd	%%xmm11,%%xmm11						\n\t	addpd	%%xmm10,%%xmm10	\n\t"\
		"addpd	%%xmm9 ,%%xmm9						\n\t	addpd	%%xmm15,%%xmm15	\n\t"\
		"addpd	%%xmm12,%%xmm12						\n\t	addpd	%%xmm14,%%xmm14	\n\t"\
		"addpd	%%xmm8 ,%%xmm8						\n\t	addpd	%%xmm13,%%xmm13	\n\t"\
		"addpd	%%xmm6 ,%%xmm11						\n\t	addpd	%%xmm2 ,%%xmm10	\n\t"\
		"addpd	%%xmm0 ,%%xmm9						\n\t	addpd	%%xmm5 ,%%xmm15	\n\t"\
		"addpd	%%xmm7 ,%%xmm12						\n\t	addpd	%%xmm4 ,%%xmm14	\n\t"\
		"addpd	%%xmm1 ,%%xmm8						\n\t	addpd	%%xmm3 ,%%xmm13	\n\t"\
		"movaps	%%xmm6 ,0x80(%%rsi)					\n\t	movaps	%%xmm2 ,0xa0(%%rsi)	\n\t"\
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
		"shlq	$3,%%rcx							\n\t	movslq	%[__p04],%%r8	\n\t"\
		"shlq	$3,%%rdx							\n\t	shlq	$3,%%r8	\n\t"\
		"addq	%%rdi,%%rax	/* &a[j1+p1] */			\n\t	movq	%%r8,%%r10	\n\t	addq	%%rax,%%r10	\n\t"\
		"addq	%%rax,%%rbx							\n\t	movq	%%r8,%%r11	\n\t	addq	%%rbx,%%r11	\n\t"\
		"addq	%%rax,%%rcx							\n\t	movq	%%r8,%%r12	\n\t	addq	%%rcx,%%r12	\n\t"\
		"addq	%%rax,%%rdx							\n\t	movq	%%r8,%%r13	\n\t	addq	%%rdx,%%r13	\n\t"\
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
		"movaps	%%xmm2,0x020(%%rsi)					\n\t	addpd	%%xmm14,%%xmm14	\n\t"\
		"movaps	%%xmm1,0x050(%%rsi)					\n\t	addpd	%%xmm13,%%xmm13	\n\t"\
		"movaps	%%xmm3,0x070(%%rsi)					\n\t	addpd	%%xmm15,%%xmm15	\n\t"\
		"addpd	%%xmm6,%%xmm6						\n\t	addpd	%%xmm12,%%xmm12	\n\t"\
		"addpd	%%xmm5,%%xmm5						\n\t	addpd	%%xmm8 ,%%xmm14	\n\t"\
		"addpd	%%xmm7,%%xmm7						\n\t	addpd	%%xmm10,%%xmm13	\n\t"\
		"addpd	%%xmm4,%%xmm4						\n\t	addpd	%%xmm9 ,%%xmm15	\n\t"\
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
		"													mulpd	%%xmm8 ,%%xmm15	\n\t"\
		"												/*	movaps	%%xmm10,0xa0(%%rsi)	*/\n\t"\
		"												/*	movaps	%%xmm13,0xe0(%%rsi)	*/\n\t"\
		"												/*	movaps	%%xmm14,0xb0(%%rsi)	*/\n\t"\
		"												/*	movaps	%%xmm15,0xf0(%%rsi)	*/\n\t"\
		"/* SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r20) */	\n\t"\
		"movaps	0x80(%%rsi),%%xmm11		\n\t"\
		"movaps	0x90(%%rsi),%%xmm12		\n\t"\
		"movaps	0xc0(%%rsi),%%xmm8		\n\t"\
		"subpd	%%xmm11,%%xmm6						\n\t	subpd	%%xmm10,%%xmm2	\n\t"\
		"subpd	%%xmm9 ,%%xmm0						\n\t	subpd	%%xmm15,%%xmm5	\n\t"\
		"subpd	%%xmm12,%%xmm7						\n\t	subpd	%%xmm14,%%xmm4	\n\t"\
		"subpd	%%xmm8 ,%%xmm1						\n\t	subpd	%%xmm13,%%xmm3	\n\t"\
		"addpd	%%xmm11,%%xmm11						\n\t	addpd	%%xmm10,%%xmm10	\n\t"\
		"addpd	%%xmm9 ,%%xmm9						\n\t	addpd	%%xmm15,%%xmm15	\n\t"\
		"addpd	%%xmm12,%%xmm12						\n\t	addpd	%%xmm14,%%xmm14	\n\t"\
		"addpd	%%xmm8 ,%%xmm8						\n\t	addpd	%%xmm13,%%xmm13	\n\t"\
		"addpd	%%xmm6 ,%%xmm11						\n\t	addpd	%%xmm2 ,%%xmm10	\n\t"\
		"addpd	%%xmm0 ,%%xmm9						\n\t	addpd	%%xmm5 ,%%xmm15	\n\t"\
		"addpd	%%xmm7 ,%%xmm12						\n\t	addpd	%%xmm4 ,%%xmm14	\n\t"\
		"addpd	%%xmm1 ,%%xmm8						\n\t	addpd	%%xmm3 ,%%xmm13	\n\t"\
		"movaps	%%xmm6 ,0x80(%%rsi)					\n\t	movaps	%%xmm2 ,0xa0(%%rsi)	\n\t"\
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
		"shlq	$3,%%rcx							\n\t	movslq	%[__p04],%%r8	\n\t"\
		"shlq	$3,%%rdx							\n\t	shlq	$3,%%r8	\n\t"\
		"addq	%%rdi,%%rax	/* &a[j1+p3] */			\n\t	movq	%%r8,%%r10	\n\t	addq	%%rax,%%r10	\n\t"\
		"addq	%%rax,%%rbx							\n\t	movq	%%r8,%%r11	\n\t	addq	%%rbx,%%r11	\n\t"\
		"addq	%%rax,%%rcx							\n\t	movq	%%r8,%%r12	\n\t	addq	%%rcx,%%r12	\n\t"\
		"addq	%%rax,%%rdx							\n\t	movq	%%r8,%%r13	\n\t	addq	%%rdx,%%r13	\n\t"\
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
		"movaps	%%xmm2,0x020(%%rsi)					\n\t	addpd	%%xmm14,%%xmm14	\n\t"\
		"movaps	%%xmm1,0x050(%%rsi)					\n\t	addpd	%%xmm13,%%xmm13	\n\t"\
		"movaps	%%xmm3,0x070(%%rsi)					\n\t	addpd	%%xmm15,%%xmm15	\n\t"\
		"addpd	%%xmm6,%%xmm6						\n\t	addpd	%%xmm12,%%xmm12	\n\t"\
		"addpd	%%xmm5,%%xmm5						\n\t	addpd	%%xmm8 ,%%xmm14	\n\t"\
		"addpd	%%xmm7,%%xmm7						\n\t	addpd	%%xmm10,%%xmm13	\n\t"\
		"addpd	%%xmm4,%%xmm4						\n\t	addpd	%%xmm9 ,%%xmm15	\n\t"\
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
		"													mulpd	%%xmm8 ,%%xmm15	\n\t"\
		"												/*	movaps	%%xmm10,0xa0(%%rsi)	*/\n\t"\
		"												/*	movaps	%%xmm13,0xe0(%%rsi)	*/\n\t"\
		"												/*	movaps	%%xmm14,0xb0(%%rsi)	*/\n\t"\
		"												/*	movaps	%%xmm15,0xf0(%%rsi)	*/\n\t"\
		"/* SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r30) */	\n\t"\
		"movaps	0x80(%%rsi),%%xmm11		\n\t"\
		"movaps	0x90(%%rsi),%%xmm12		\n\t"\
		"movaps	0xc0(%%rsi),%%xmm8		\n\t"\
		"subpd	%%xmm11,%%xmm6						\n\t	subpd	%%xmm10,%%xmm2	\n\t"\
		"subpd	%%xmm9 ,%%xmm0						\n\t	subpd	%%xmm15,%%xmm5	\n\t"\
		"subpd	%%xmm12,%%xmm7						\n\t	subpd	%%xmm14,%%xmm4	\n\t"\
		"subpd	%%xmm8 ,%%xmm1						\n\t	subpd	%%xmm13,%%xmm3	\n\t"\
		"addpd	%%xmm11,%%xmm11						\n\t	addpd	%%xmm10,%%xmm10	\n\t"\
		"addpd	%%xmm9 ,%%xmm9						\n\t	addpd	%%xmm15,%%xmm15	\n\t"\
		"addpd	%%xmm12,%%xmm12						\n\t	addpd	%%xmm14,%%xmm14	\n\t"\
		"addpd	%%xmm8 ,%%xmm8						\n\t	addpd	%%xmm13,%%xmm13	\n\t"\
		"addpd	%%xmm6 ,%%xmm11						\n\t	addpd	%%xmm2 ,%%xmm10	\n\t"\
		"addpd	%%xmm0 ,%%xmm9						\n\t	addpd	%%xmm5 ,%%xmm15	\n\t"\
		"addpd	%%xmm7 ,%%xmm12						\n\t	addpd	%%xmm4 ,%%xmm14	\n\t"\
		"addpd	%%xmm1 ,%%xmm8						\n\t	addpd	%%xmm3 ,%%xmm13	\n\t"\
		"movaps	%%xmm6 ,0x80(%%rsi)					\n\t	movaps	%%xmm2 ,0xa0(%%rsi)	\n\t"\
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
		"/*...Block 1: t00,t10,t20,t30	*/	\n\t"\
		"movq	%[__add0],%%rax	/* &a[j1] */\n\t"\
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
		: "rax","rbx","rcx","rdx","rdi","rsi","r8","r10","r11","r12","r13","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"		/* Clobbered registers */\
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
		"movq	%[__add0],%%rax		\n\t"\
		"movslq	%[__p05],%%rbx		\n\t"\
		"movslq	%[__p06],%%rcx		\n\t"\
		"movslq	%[__p07],%%rdx		\n\t"\
		"movslq	%[__p04],%%rdi	/* edi will store copy of p4 throughout */\n\t"\
		"shlq	$3,%%rbx			\n\t"\
		"shlq	$3,%%rcx			\n\t"\
		"shlq	$3,%%rdx			\n\t"\
		"shlq	$3,%%rdi			\n\t"\
		"addq	%%rax,%%rbx			\n\t"\
		"addq	%%rax,%%rcx			\n\t"\
		"addq	%%rax,%%rdx			\n\t"\
		"addq	%%rdi,%%rax			\n\t"\
		"/* SSE2_RADIX8_DIT_0TWIDDLE(r00): */\n\t"\
		"movq	%[__r00]	,%%rsi 			\n\t"\
		"movaps	    (%%rax)	,%%xmm0			\n\t"\
		"movaps	0x10(%%rax)	,%%xmm1			\n\t"\
		"movaps	%%xmm0		,%%xmm2			\n\t"\
		"movaps	%%xmm1		,%%xmm3			\n\t"\
		"addpd	    (%%rbx)	,%%xmm2			\n\t"\
		"addpd	0x10(%%rbx)	,%%xmm3			\n\t"\
		"subpd	    (%%rbx)	,%%xmm0			\n\t"\
		"subpd	0x10(%%rbx)	,%%xmm1			\n\t"\
		"movaps	    (%%rcx)	,%%xmm4			\n\t"\
		"movaps	0x10(%%rcx)	,%%xmm5			\n\t"\
		"movaps	%%xmm4		,%%xmm6			\n\t"\
		"movaps	%%xmm5		,%%xmm7			\n\t"\
		"addpd	    (%%rdx)	,%%xmm6			\n\t"\
		"addpd	0x10(%%rdx)	,%%xmm7			\n\t"\
		"subpd	    (%%rdx)	,%%xmm4			\n\t"\
		"subpd	0x10(%%rdx)	,%%xmm5			\n\t"\
		"movaps	%%xmm6		,0xc0(%%rsi)	\n\t"\
		"movaps	%%xmm7		,0xd0(%%rsi)	\n\t"\
		"addpd	%%xmm2		,%%xmm6			\n\t"\
		"addpd	%%xmm3		,%%xmm7			\n\t"\
		"subpd	0xc0(%%rsi)	,%%xmm2			\n\t"\
		"subpd	0xd0(%%rsi)	,%%xmm3			\n\t"\
		"movaps	%%xmm2		,0xc0(%%rsi)	\n\t"\
		"movaps	%%xmm3		,0xd0(%%rsi)	\n\t"\
		"movaps	%%xmm4		,%%xmm2			\n\t"\
		"movaps	%%xmm5		,%%xmm3			\n\t"\
		"addpd	%%xmm0		,%%xmm5			\n\t"\
		"subpd	%%xmm3		,%%xmm0			\n\t"\
		"addpd	%%xmm1		,%%xmm4			\n\t"\
		"subpd	%%xmm2		,%%xmm1			\n\t"\
		"movq	%[__isrt2]	,%%rsi			\n\t"\
		"movaps	%%xmm5		,%%xmm2			\n\t"\
		"movaps	%%xmm1		,%%xmm3			\n\t"\
		"addpd	%%xmm1		,%%xmm5			\n\t"\
		"movaps	(%%rsi)		,%%xmm1			\n\t"\
		"subpd	%%xmm3		,%%xmm2			\n\t"\
		"movq	%[__r00]	,%%rsi			\n\t"\
		"mulpd	%%xmm1		,%%xmm5			\n\t"\
		"mulpd	%%xmm1		,%%xmm2			\n\t"\
		"movaps	%%xmm0		,%%xmm3			\n\t"\
		"movaps	%%xmm5		,0xa0(%%rsi)	\n\t"\
		"movaps	%%xmm4		,%%xmm5			\n\t"\
		"addpd	%%xmm4		,%%xmm0			\n\t"\
		"movaps	%%xmm2		,0xb0(%%rsi)	\n\t"\
		"subpd	%%xmm5		,%%xmm3			\n\t"\
		"mulpd	%%xmm1		,%%xmm0			\n\t"\
		"mulpd	%%xmm1		,%%xmm3			\n\t"\
		"movaps	%%xmm0		,0xe0(%%rsi)	\n\t"\
		"movaps	%%xmm3		,0xf0(%%rsi)	\n\t"\
		"subq	%%rdi		,%%rax			\n\t"\
		"subq	%%rdi		,%%rbx			\n\t"\
		"subq	%%rdi		,%%rcx			\n\t"\
		"subq	%%rdi		,%%rdx			\n\t"\
		"movaps	    (%%rax)	,%%xmm0			\n\t"\
		"movaps	0x10(%%rax)	,%%xmm1			\n\t"\
		"movaps	%%xmm0		,%%xmm2			\n\t"\
		"movaps	%%xmm1		,%%xmm3			\n\t"\
		"addpd	    (%%rbx)	,%%xmm2			\n\t"\
		"addpd	0x10(%%rbx)	,%%xmm3			\n\t"\
		"subpd	    (%%rbx)	,%%xmm0			\n\t"\
		"subpd	0x10(%%rbx)	,%%xmm1			\n\t"\
		"movaps	%%xmm6		,    (%%rsi)	\n\t"\
		"movaps	%%xmm7		,0x10(%%rsi)	\n\t"\
		"addpd	%%xmm6		,%%xmm6			\n\t"\
		"addpd	%%xmm7		,%%xmm7			\n\t"\
		"movaps	%%xmm6		,0x80(%%rsi)	\n\t"\
		"movaps	%%xmm7		,0x90(%%rsi)	\n\t"\
		"movaps	    (%%rcx)	,%%xmm4			\n\t"\
		"movaps	0x10(%%rcx)	,%%xmm5			\n\t"\
		"movaps	%%xmm4		,%%xmm6			\n\t"\
		"movaps	%%xmm5		,%%xmm7			\n\t"\
		"addpd	    (%%rdx)	,%%xmm6			\n\t"\
		"addpd	0x10(%%rdx)	,%%xmm7			\n\t"\
		"subpd	    (%%rdx)	,%%xmm4			\n\t"\
		"subpd	0x10(%%rdx)	,%%xmm5			\n\t"\
		"movaps	%%xmm6		,0x40(%%rsi)	\n\t"\
		"movaps	%%xmm7		,0x50(%%rsi)	\n\t"\
		"addpd	%%xmm2		,%%xmm6			\n\t"\
		"addpd	%%xmm3		,%%xmm7			\n\t"\
		"subpd	0x40(%%rsi)	,%%xmm2			\n\t"\
		"subpd	0x50(%%rsi)	,%%xmm3			\n\t"\
		"addpd	    (%%rsi)	,%%xmm6			\n\t"\
		"addpd	0x10(%%rsi)	,%%xmm7			\n\t"\
		"movaps	%%xmm6		,    (%%rsi)	\n\t"\
		"movaps	%%xmm7		,0x10(%%rsi)	\n\t"\
		"subpd	0x80(%%rsi)	,%%xmm6			\n\t"\
		"subpd	0x90(%%rsi)	,%%xmm7			\n\t"\
		"movaps	%%xmm6		,0x80(%%rsi)	\n\t"\
		"movaps	%%xmm7		,0x90(%%rsi)	\n\t"\
		"movaps	%%xmm4		,%%xmm6			\n\t"\
		"movaps	%%xmm5		,%%xmm7			\n\t"\
		"addpd	%%xmm0		,%%xmm5			\n\t"\
		"subpd	%%xmm7		,%%xmm0			\n\t"\
		"addpd	%%xmm1		,%%xmm4			\n\t"\
		"subpd	%%xmm6		,%%xmm1			\n\t"\
		"movaps	%%xmm5		,%%xmm6			\n\t"\
		"movaps	%%xmm1		,%%xmm7			\n\t"\
		"addpd	0xa0(%%rsi)	,%%xmm5			\n\t"\
		"subpd	0xb0(%%rsi)	,%%xmm1			\n\t"\
		"subpd	0xa0(%%rsi)	,%%xmm6			\n\t"\
		"addpd	0xb0(%%rsi)	,%%xmm7			\n\t"\
		"movaps	%%xmm5		,0x20(%%rsi)	\n\t"\
		"movaps	%%xmm1		,0x30(%%rsi)	\n\t"\
		"movaps	%%xmm6		,0xa0(%%rsi)	\n\t"\
		"movaps	%%xmm7		,0xb0(%%rsi)	\n\t"\
		"movaps	%%xmm2		,%%xmm6			\n\t"\
		"movaps	%%xmm3		,%%xmm7			\n\t"\
		"addpd	0xd0(%%rsi)	,%%xmm2			\n\t"\
		"subpd	0xc0(%%rsi)	,%%xmm3			\n\t"\
		"subpd	0xd0(%%rsi)	,%%xmm6			\n\t"\
		"addpd	0xc0(%%rsi)	,%%xmm7			\n\t"\
		"movaps	%%xmm2		,0x40(%%rsi)	\n\t"\
		"movaps	%%xmm3		,0x50(%%rsi)	\n\t"\
		"movaps	%%xmm6		,0xc0(%%rsi)	\n\t"\
		"movaps	%%xmm7		,0xd0(%%rsi)	\n\t"\
		"movaps	%%xmm0		,%%xmm6			\n\t"\
		"movaps	%%xmm4		,%%xmm7			\n\t"\
		"subpd	0xf0(%%rsi)	,%%xmm0			\n\t"\
		"subpd	0xe0(%%rsi)	,%%xmm4			\n\t"\
		"addpd	0xf0(%%rsi)	,%%xmm6			\n\t"\
		"addpd	0xe0(%%rsi)	,%%xmm7			\n\t"\
		"movaps	%%xmm0		,0x60(%%rsi)	\n\t"\
		"movaps	%%xmm4		,0x70(%%rsi)	\n\t"\
		"movaps	%%xmm6		,0xe0(%%rsi)	\n\t"\
		"movaps	%%xmm7		,0xf0(%%rsi)	\n\t"\
		"\n\t"\
		"movq	%[__add0],%%rax		\n\t"\
		"movslq	%[__p08],%%rsi		\n\t"\
		"shlq	$3,%%rsi			\n\t"\
		"addq	%%rsi,%%rax	/* add1 = add0+p8 */\n\t"\
		"movslq	%[__p05],%%rbx		\n\t"\
		"movslq	%[__p06],%%rcx		\n\t"\
		"movslq	%[__p07],%%rdx		\n\t"\
		"movslq	%[__p04],%%rdi		\n\t"\
		"shlq	$3,%%rbx			\n\t"\
		"shlq	$3,%%rcx			\n\t"\
		"shlq	$3,%%rdx			\n\t"\
		"shlq	$3,%%rdi			\n\t"\
		"addq	%%rax,%%rbx			\n\t"\
		"addq	%%rax,%%rcx			\n\t"\
		"addq	%%rax,%%rdx			\n\t"\
		"addq	%%rdi,%%rax			\n\t"\
		"/* SSE2_RADIX8_DIT_0TWIDDLE(r10): */\n\t"\
		"movq	%[__r10]	,%%rsi 			\n\t"\
		"movaps	    (%%rax)	,%%xmm0			\n\t"\
		"movaps	0x10(%%rax)	,%%xmm1			\n\t"\
		"movaps	%%xmm0		,%%xmm2			\n\t"\
		"movaps	%%xmm1		,%%xmm3			\n\t"\
		"addpd	    (%%rbx)	,%%xmm2			\n\t"\
		"addpd	0x10(%%rbx)	,%%xmm3			\n\t"\
		"subpd	    (%%rbx)	,%%xmm0			\n\t"\
		"subpd	0x10(%%rbx)	,%%xmm1			\n\t"\
		"movaps	    (%%rcx)	,%%xmm4			\n\t"\
		"movaps	0x10(%%rcx)	,%%xmm5			\n\t"\
		"movaps	%%xmm4		,%%xmm6			\n\t"\
		"movaps	%%xmm5		,%%xmm7			\n\t"\
		"addpd	    (%%rdx)	,%%xmm6			\n\t"\
		"addpd	0x10(%%rdx)	,%%xmm7			\n\t"\
		"subpd	    (%%rdx)	,%%xmm4			\n\t"\
		"subpd	0x10(%%rdx)	,%%xmm5			\n\t"\
		"movaps	%%xmm6		,0xc0(%%rsi)	\n\t"\
		"movaps	%%xmm7		,0xd0(%%rsi)	\n\t"\
		"addpd	%%xmm2		,%%xmm6			\n\t"\
		"addpd	%%xmm3		,%%xmm7			\n\t"\
		"subpd	0xc0(%%rsi)	,%%xmm2			\n\t"\
		"subpd	0xd0(%%rsi)	,%%xmm3			\n\t"\
		"movaps	%%xmm2		,0xc0(%%rsi)	\n\t"\
		"movaps	%%xmm3		,0xd0(%%rsi)	\n\t"\
		"movaps	%%xmm4		,%%xmm2			\n\t"\
		"movaps	%%xmm5		,%%xmm3			\n\t"\
		"addpd	%%xmm0		,%%xmm5			\n\t"\
		"subpd	%%xmm3		,%%xmm0			\n\t"\
		"addpd	%%xmm1		,%%xmm4			\n\t"\
		"subpd	%%xmm2		,%%xmm1			\n\t"\
		"movq	%[__isrt2]	,%%rsi			\n\t"\
		"movaps	%%xmm5		,%%xmm2			\n\t"\
		"movaps	%%xmm1		,%%xmm3			\n\t"\
		"addpd	%%xmm1		,%%xmm5			\n\t"\
		"movaps	(%%rsi)		,%%xmm1			\n\t"\
		"subpd	%%xmm3		,%%xmm2			\n\t"\
		"movq	%[__r10]	,%%rsi			\n\t"\
		"mulpd	%%xmm1		,%%xmm5			\n\t"\
		"mulpd	%%xmm1		,%%xmm2			\n\t"\
		"movaps	%%xmm0		,%%xmm3			\n\t"\
		"movaps	%%xmm5		,0xa0(%%rsi)	\n\t"\
		"movaps	%%xmm4		,%%xmm5			\n\t"\
		"addpd	%%xmm4		,%%xmm0			\n\t"\
		"movaps	%%xmm2		,0xb0(%%rsi)	\n\t"\
		"subpd	%%xmm5		,%%xmm3			\n\t"\
		"mulpd	%%xmm1		,%%xmm0			\n\t"\
		"mulpd	%%xmm1		,%%xmm3			\n\t"\
		"movaps	%%xmm0		,0xe0(%%rsi)	\n\t"\
		"movaps	%%xmm3		,0xf0(%%rsi)	\n\t"\
		"subq	%%rdi		,%%rax			\n\t"\
		"subq	%%rdi		,%%rbx			\n\t"\
		"subq	%%rdi		,%%rcx			\n\t"\
		"subq	%%rdi		,%%rdx			\n\t"\
		"movaps	    (%%rax)	,%%xmm0			\n\t"\
		"movaps	0x10(%%rax)	,%%xmm1			\n\t"\
		"movaps	%%xmm0		,%%xmm2			\n\t"\
		"movaps	%%xmm1		,%%xmm3			\n\t"\
		"addpd	    (%%rbx)	,%%xmm2			\n\t"\
		"addpd	0x10(%%rbx)	,%%xmm3			\n\t"\
		"subpd	    (%%rbx)	,%%xmm0			\n\t"\
		"subpd	0x10(%%rbx)	,%%xmm1			\n\t"\
		"movaps	%%xmm6		,    (%%rsi)	\n\t"\
		"movaps	%%xmm7		,0x10(%%rsi)	\n\t"\
		"addpd	%%xmm6		,%%xmm6			\n\t"\
		"addpd	%%xmm7		,%%xmm7			\n\t"\
		"movaps	%%xmm6		,0x80(%%rsi)	\n\t"\
		"movaps	%%xmm7		,0x90(%%rsi)	\n\t"\
		"movaps	    (%%rcx)	,%%xmm4			\n\t"\
		"movaps	0x10(%%rcx)	,%%xmm5			\n\t"\
		"movaps	%%xmm4		,%%xmm6			\n\t"\
		"movaps	%%xmm5		,%%xmm7			\n\t"\
		"addpd	    (%%rdx)	,%%xmm6			\n\t"\
		"addpd	0x10(%%rdx)	,%%xmm7			\n\t"\
		"subpd	    (%%rdx)	,%%xmm4			\n\t"\
		"subpd	0x10(%%rdx)	,%%xmm5			\n\t"\
		"movaps	%%xmm6		,0x40(%%rsi)	\n\t"\
		"movaps	%%xmm7		,0x50(%%rsi)	\n\t"\
		"addpd	%%xmm2		,%%xmm6			\n\t"\
		"addpd	%%xmm3		,%%xmm7			\n\t"\
		"subpd	0x40(%%rsi)	,%%xmm2			\n\t"\
		"subpd	0x50(%%rsi)	,%%xmm3			\n\t"\
		"addpd	    (%%rsi)	,%%xmm6			\n\t"\
		"addpd	0x10(%%rsi)	,%%xmm7			\n\t"\
		"movaps	%%xmm6		,    (%%rsi)	\n\t"\
		"movaps	%%xmm7		,0x10(%%rsi)	\n\t"\
		"subpd	0x80(%%rsi)	,%%xmm6			\n\t"\
		"subpd	0x90(%%rsi)	,%%xmm7			\n\t"\
		"movaps	%%xmm6		,0x80(%%rsi)	\n\t"\
		"movaps	%%xmm7		,0x90(%%rsi)	\n\t"\
		"movaps	%%xmm4		,%%xmm6			\n\t"\
		"movaps	%%xmm5		,%%xmm7			\n\t"\
		"addpd	%%xmm0		,%%xmm5			\n\t"\
		"subpd	%%xmm7		,%%xmm0			\n\t"\
		"addpd	%%xmm1		,%%xmm4			\n\t"\
		"subpd	%%xmm6		,%%xmm1			\n\t"\
		"movaps	%%xmm5		,%%xmm6			\n\t"\
		"movaps	%%xmm1		,%%xmm7			\n\t"\
		"addpd	0xa0(%%rsi)	,%%xmm5			\n\t"\
		"subpd	0xb0(%%rsi)	,%%xmm1			\n\t"\
		"subpd	0xa0(%%rsi)	,%%xmm6			\n\t"\
		"addpd	0xb0(%%rsi)	,%%xmm7			\n\t"\
		"movaps	%%xmm5		,0x20(%%rsi)	\n\t"\
		"movaps	%%xmm1		,0x30(%%rsi)	\n\t"\
		"movaps	%%xmm6		,0xa0(%%rsi)	\n\t"\
		"movaps	%%xmm7		,0xb0(%%rsi)	\n\t"\
		"movaps	%%xmm2		,%%xmm6			\n\t"\
		"movaps	%%xmm3		,%%xmm7			\n\t"\
		"addpd	0xd0(%%rsi)	,%%xmm2			\n\t"\
		"subpd	0xc0(%%rsi)	,%%xmm3			\n\t"\
		"subpd	0xd0(%%rsi)	,%%xmm6			\n\t"\
		"addpd	0xc0(%%rsi)	,%%xmm7			\n\t"\
		"movaps	%%xmm2		,0x40(%%rsi)	\n\t"\
		"movaps	%%xmm3		,0x50(%%rsi)	\n\t"\
		"movaps	%%xmm6		,0xc0(%%rsi)	\n\t"\
		"movaps	%%xmm7		,0xd0(%%rsi)	\n\t"\
		"movaps	%%xmm0		,%%xmm6			\n\t"\
		"movaps	%%xmm4		,%%xmm7			\n\t"\
		"subpd	0xf0(%%rsi)	,%%xmm0			\n\t"\
		"subpd	0xe0(%%rsi)	,%%xmm4			\n\t"\
		"addpd	0xf0(%%rsi)	,%%xmm6			\n\t"\
		"addpd	0xe0(%%rsi)	,%%xmm7			\n\t"\
		"movaps	%%xmm0		,0x60(%%rsi)	\n\t"\
		"movaps	%%xmm4		,0x70(%%rsi)	\n\t"\
		"movaps	%%xmm6		,0xe0(%%rsi)	\n\t"\
		"movaps	%%xmm7		,0xf0(%%rsi)	\n\t"\
		"\n\t"\
		"movq	%[__add0],%%rax		\n\t"\
		"movslq	%[__p10],%%rsi		\n\t"\
		"shlq	$3,%%rsi			\n\t"\
		"addq	%%rsi,%%rax	/* add2 = add0+p10*/\n\t"\
		"movslq	%[__p05],%%rbx		\n\t"\
		"movslq	%[__p06],%%rcx		\n\t"\
		"movslq	%[__p07],%%rdx		\n\t"\
		"movslq	%[__p04],%%rdi		\n\t"\
		"shlq	$3,%%rbx			\n\t"\
		"shlq	$3,%%rcx			\n\t"\
		"shlq	$3,%%rdx			\n\t"\
		"shlq	$3,%%rdi			\n\t"\
		"addq	%%rax,%%rbx			\n\t"\
		"addq	%%rax,%%rcx			\n\t"\
		"addq	%%rax,%%rdx			\n\t"\
		"addq	%%rdi,%%rax			\n\t"\
		"/* SSE2_RADIX8_DIT_0TWIDDLE(r20): */\n\t"\
		"movq	%[__r20]	,%%rsi 			\n\t"\
		"movaps	    (%%rax)	,%%xmm0			\n\t"\
		"movaps	0x10(%%rax)	,%%xmm1			\n\t"\
		"movaps	%%xmm0		,%%xmm2			\n\t"\
		"movaps	%%xmm1		,%%xmm3			\n\t"\
		"addpd	    (%%rbx)	,%%xmm2			\n\t"\
		"addpd	0x10(%%rbx)	,%%xmm3			\n\t"\
		"subpd	    (%%rbx)	,%%xmm0			\n\t"\
		"subpd	0x10(%%rbx)	,%%xmm1			\n\t"\
		"movaps	    (%%rcx)	,%%xmm4			\n\t"\
		"movaps	0x10(%%rcx)	,%%xmm5			\n\t"\
		"movaps	%%xmm4		,%%xmm6			\n\t"\
		"movaps	%%xmm5		,%%xmm7			\n\t"\
		"addpd	    (%%rdx)	,%%xmm6			\n\t"\
		"addpd	0x10(%%rdx)	,%%xmm7			\n\t"\
		"subpd	    (%%rdx)	,%%xmm4			\n\t"\
		"subpd	0x10(%%rdx)	,%%xmm5			\n\t"\
		"movaps	%%xmm6		,0xc0(%%rsi)	\n\t"\
		"movaps	%%xmm7		,0xd0(%%rsi)	\n\t"\
		"addpd	%%xmm2		,%%xmm6			\n\t"\
		"addpd	%%xmm3		,%%xmm7			\n\t"\
		"subpd	0xc0(%%rsi)	,%%xmm2			\n\t"\
		"subpd	0xd0(%%rsi)	,%%xmm3			\n\t"\
		"movaps	%%xmm2		,0xc0(%%rsi)	\n\t"\
		"movaps	%%xmm3		,0xd0(%%rsi)	\n\t"\
		"movaps	%%xmm4		,%%xmm2			\n\t"\
		"movaps	%%xmm5		,%%xmm3			\n\t"\
		"addpd	%%xmm0		,%%xmm5			\n\t"\
		"subpd	%%xmm3		,%%xmm0			\n\t"\
		"addpd	%%xmm1		,%%xmm4			\n\t"\
		"subpd	%%xmm2		,%%xmm1			\n\t"\
		"movq	%[__isrt2]	,%%rsi			\n\t"\
		"movaps	%%xmm5		,%%xmm2			\n\t"\
		"movaps	%%xmm1		,%%xmm3			\n\t"\
		"addpd	%%xmm1		,%%xmm5			\n\t"\
		"movaps	(%%rsi)		,%%xmm1			\n\t"\
		"subpd	%%xmm3		,%%xmm2			\n\t"\
		"movq	%[__r20]	,%%rsi			\n\t"\
		"mulpd	%%xmm1		,%%xmm5			\n\t"\
		"mulpd	%%xmm1		,%%xmm2			\n\t"\
		"movaps	%%xmm0		,%%xmm3			\n\t"\
		"movaps	%%xmm5		,0xa0(%%rsi)	\n\t"\
		"movaps	%%xmm4		,%%xmm5			\n\t"\
		"addpd	%%xmm4		,%%xmm0			\n\t"\
		"movaps	%%xmm2		,0xb0(%%rsi)	\n\t"\
		"subpd	%%xmm5		,%%xmm3			\n\t"\
		"mulpd	%%xmm1		,%%xmm0			\n\t"\
		"mulpd	%%xmm1		,%%xmm3			\n\t"\
		"movaps	%%xmm0		,0xe0(%%rsi)	\n\t"\
		"movaps	%%xmm3		,0xf0(%%rsi)	\n\t"\
		"subq	%%rdi		,%%rax			\n\t"\
		"subq	%%rdi		,%%rbx			\n\t"\
		"subq	%%rdi		,%%rcx			\n\t"\
		"subq	%%rdi		,%%rdx			\n\t"\
		"movaps	    (%%rax)	,%%xmm0			\n\t"\
		"movaps	0x10(%%rax)	,%%xmm1			\n\t"\
		"movaps	%%xmm0		,%%xmm2			\n\t"\
		"movaps	%%xmm1		,%%xmm3			\n\t"\
		"addpd	    (%%rbx)	,%%xmm2			\n\t"\
		"addpd	0x10(%%rbx)	,%%xmm3			\n\t"\
		"subpd	    (%%rbx)	,%%xmm0			\n\t"\
		"subpd	0x10(%%rbx)	,%%xmm1			\n\t"\
		"movaps	%%xmm6		,    (%%rsi)	\n\t"\
		"movaps	%%xmm7		,0x10(%%rsi)	\n\t"\
		"addpd	%%xmm6		,%%xmm6			\n\t"\
		"addpd	%%xmm7		,%%xmm7			\n\t"\
		"movaps	%%xmm6		,0x80(%%rsi)	\n\t"\
		"movaps	%%xmm7		,0x90(%%rsi)	\n\t"\
		"movaps	    (%%rcx)	,%%xmm4			\n\t"\
		"movaps	0x10(%%rcx)	,%%xmm5			\n\t"\
		"movaps	%%xmm4		,%%xmm6			\n\t"\
		"movaps	%%xmm5		,%%xmm7			\n\t"\
		"addpd	    (%%rdx)	,%%xmm6			\n\t"\
		"addpd	0x10(%%rdx)	,%%xmm7			\n\t"\
		"subpd	    (%%rdx)	,%%xmm4			\n\t"\
		"subpd	0x10(%%rdx)	,%%xmm5			\n\t"\
		"movaps	%%xmm6		,0x40(%%rsi)	\n\t"\
		"movaps	%%xmm7		,0x50(%%rsi)	\n\t"\
		"addpd	%%xmm2		,%%xmm6			\n\t"\
		"addpd	%%xmm3		,%%xmm7			\n\t"\
		"subpd	0x40(%%rsi)	,%%xmm2			\n\t"\
		"subpd	0x50(%%rsi)	,%%xmm3			\n\t"\
		"addpd	    (%%rsi)	,%%xmm6			\n\t"\
		"addpd	0x10(%%rsi)	,%%xmm7			\n\t"\
		"movaps	%%xmm6		,    (%%rsi)	\n\t"\
		"movaps	%%xmm7		,0x10(%%rsi)	\n\t"\
		"subpd	0x80(%%rsi)	,%%xmm6			\n\t"\
		"subpd	0x90(%%rsi)	,%%xmm7			\n\t"\
		"movaps	%%xmm6		,0x80(%%rsi)	\n\t"\
		"movaps	%%xmm7		,0x90(%%rsi)	\n\t"\
		"movaps	%%xmm4		,%%xmm6			\n\t"\
		"movaps	%%xmm5		,%%xmm7			\n\t"\
		"addpd	%%xmm0		,%%xmm5			\n\t"\
		"subpd	%%xmm7		,%%xmm0			\n\t"\
		"addpd	%%xmm1		,%%xmm4			\n\t"\
		"subpd	%%xmm6		,%%xmm1			\n\t"\
		"movaps	%%xmm5		,%%xmm6			\n\t"\
		"movaps	%%xmm1		,%%xmm7			\n\t"\
		"addpd	0xa0(%%rsi)	,%%xmm5			\n\t"\
		"subpd	0xb0(%%rsi)	,%%xmm1			\n\t"\
		"subpd	0xa0(%%rsi)	,%%xmm6			\n\t"\
		"addpd	0xb0(%%rsi)	,%%xmm7			\n\t"\
		"movaps	%%xmm5		,0x20(%%rsi)	\n\t"\
		"movaps	%%xmm1		,0x30(%%rsi)	\n\t"\
		"movaps	%%xmm6		,0xa0(%%rsi)	\n\t"\
		"movaps	%%xmm7		,0xb0(%%rsi)	\n\t"\
		"movaps	%%xmm2		,%%xmm6			\n\t"\
		"movaps	%%xmm3		,%%xmm7			\n\t"\
		"addpd	0xd0(%%rsi)	,%%xmm2			\n\t"\
		"subpd	0xc0(%%rsi)	,%%xmm3			\n\t"\
		"subpd	0xd0(%%rsi)	,%%xmm6			\n\t"\
		"addpd	0xc0(%%rsi)	,%%xmm7			\n\t"\
		"movaps	%%xmm2		,0x40(%%rsi)	\n\t"\
		"movaps	%%xmm3		,0x50(%%rsi)	\n\t"\
		"movaps	%%xmm6		,0xc0(%%rsi)	\n\t"\
		"movaps	%%xmm7		,0xd0(%%rsi)	\n\t"\
		"movaps	%%xmm0		,%%xmm6			\n\t"\
		"movaps	%%xmm4		,%%xmm7			\n\t"\
		"subpd	0xf0(%%rsi)	,%%xmm0			\n\t"\
		"subpd	0xe0(%%rsi)	,%%xmm4			\n\t"\
		"addpd	0xf0(%%rsi)	,%%xmm6			\n\t"\
		"addpd	0xe0(%%rsi)	,%%xmm7			\n\t"\
		"movaps	%%xmm0		,0x60(%%rsi)	\n\t"\
		"movaps	%%xmm4		,0x70(%%rsi)	\n\t"\
		"movaps	%%xmm6		,0xe0(%%rsi)	\n\t"\
		"movaps	%%xmm7		,0xf0(%%rsi)	\n\t"\
		"\n\t"\
		"movq	%[__add0],%%rax			\n\t"\
		"movslq	%[__p18],%%rsi			\n\t"\
		"shlq	$3,%%rsi			\n\t"\
		"addq	%%rsi,%%rax	/* add3 = add0+p18*/\n\t"\
		"movslq	%[__p05],%%rbx			\n\t"\
		"movslq	%[__p06],%%rcx			\n\t"\
		"movslq	%[__p07],%%rdx			\n\t"\
		"movslq	%[__p04],%%rdi			\n\t"\
		"shlq	$3,%%rbx			\n\t"\
		"shlq	$3,%%rcx			\n\t"\
		"shlq	$3,%%rdx			\n\t"\
		"shlq	$3,%%rdi			\n\t"\
		"addq	%%rax,%%rbx			\n\t"\
		"addq	%%rax,%%rcx			\n\t"\
		"addq	%%rax,%%rdx			\n\t"\
		"addq	%%rdi,%%rax			\n\t"\
		"/* SSE2_RADIX8_DIT_0TWIDDLE(r30): */\n\t"\
		"movq	%[__r30]	,%%rsi 			\n\t"\
		"movaps	    (%%rax)	,%%xmm0			\n\t"\
		"movaps	0x10(%%rax)	,%%xmm1			\n\t"\
		"movaps	%%xmm0		,%%xmm2			\n\t"\
		"movaps	%%xmm1		,%%xmm3			\n\t"\
		"addpd	    (%%rbx)	,%%xmm2			\n\t"\
		"addpd	0x10(%%rbx)	,%%xmm3			\n\t"\
		"subpd	    (%%rbx)	,%%xmm0			\n\t"\
		"subpd	0x10(%%rbx)	,%%xmm1			\n\t"\
		"movaps	    (%%rcx)	,%%xmm4			\n\t"\
		"movaps	0x10(%%rcx)	,%%xmm5			\n\t"\
		"movaps	%%xmm4		,%%xmm6			\n\t"\
		"movaps	%%xmm5		,%%xmm7			\n\t"\
		"addpd	    (%%rdx)	,%%xmm6			\n\t"\
		"addpd	0x10(%%rdx)	,%%xmm7			\n\t"\
		"subpd	    (%%rdx)	,%%xmm4			\n\t"\
		"subpd	0x10(%%rdx)	,%%xmm5			\n\t"\
		"movaps	%%xmm6		,0xc0(%%rsi)	\n\t"\
		"movaps	%%xmm7		,0xd0(%%rsi)	\n\t"\
		"addpd	%%xmm2		,%%xmm6			\n\t"\
		"addpd	%%xmm3		,%%xmm7			\n\t"\
		"subpd	0xc0(%%rsi)	,%%xmm2			\n\t"\
		"subpd	0xd0(%%rsi)	,%%xmm3			\n\t"\
		"movaps	%%xmm2		,0xc0(%%rsi)	\n\t"\
		"movaps	%%xmm3		,0xd0(%%rsi)	\n\t"\
		"movaps	%%xmm4		,%%xmm2			\n\t"\
		"movaps	%%xmm5		,%%xmm3			\n\t"\
		"addpd	%%xmm0		,%%xmm5			\n\t"\
		"subpd	%%xmm3		,%%xmm0			\n\t"\
		"addpd	%%xmm1		,%%xmm4			\n\t"\
		"subpd	%%xmm2		,%%xmm1			\n\t"\
		"movq	%[__isrt2]	,%%rsi			\n\t"\
		"movaps	%%xmm5		,%%xmm2			\n\t"\
		"movaps	%%xmm1		,%%xmm3			\n\t"\
		"addpd	%%xmm1		,%%xmm5			\n\t"\
		"movaps	(%%rsi)		,%%xmm1			\n\t"\
		"subpd	%%xmm3		,%%xmm2			\n\t"\
		"movq	%[__r30]	,%%rsi			\n\t"\
		"mulpd	%%xmm1		,%%xmm5			\n\t"\
		"mulpd	%%xmm1		,%%xmm2			\n\t"\
		"movaps	%%xmm0		,%%xmm3			\n\t"\
		"movaps	%%xmm5		,0xa0(%%rsi)	\n\t"\
		"movaps	%%xmm4		,%%xmm5			\n\t"\
		"addpd	%%xmm4		,%%xmm0			\n\t"\
		"movaps	%%xmm2		,0xb0(%%rsi)	\n\t"\
		"subpd	%%xmm5		,%%xmm3			\n\t"\
		"mulpd	%%xmm1		,%%xmm0			\n\t"\
		"mulpd	%%xmm1		,%%xmm3			\n\t"\
		"movaps	%%xmm0		,0xe0(%%rsi)	\n\t"\
		"movaps	%%xmm3		,0xf0(%%rsi)	\n\t"\
		"subq	%%rdi		,%%rax			\n\t"\
		"subq	%%rdi		,%%rbx			\n\t"\
		"subq	%%rdi		,%%rcx			\n\t"\
		"subq	%%rdi		,%%rdx			\n\t"\
		"movaps	    (%%rax)	,%%xmm0			\n\t"\
		"movaps	0x10(%%rax)	,%%xmm1			\n\t"\
		"movaps	%%xmm0		,%%xmm2			\n\t"\
		"movaps	%%xmm1		,%%xmm3			\n\t"\
		"addpd	    (%%rbx)	,%%xmm2			\n\t"\
		"addpd	0x10(%%rbx)	,%%xmm3			\n\t"\
		"subpd	    (%%rbx)	,%%xmm0			\n\t"\
		"subpd	0x10(%%rbx)	,%%xmm1			\n\t"\
		"movaps	%%xmm6		,    (%%rsi)	\n\t"\
		"movaps	%%xmm7		,0x10(%%rsi)	\n\t"\
		"addpd	%%xmm6		,%%xmm6			\n\t"\
		"addpd	%%xmm7		,%%xmm7			\n\t"\
		"movaps	%%xmm6		,0x80(%%rsi)	\n\t"\
		"movaps	%%xmm7		,0x90(%%rsi)	\n\t"\
		"movaps	    (%%rcx)	,%%xmm4			\n\t"\
		"movaps	0x10(%%rcx)	,%%xmm5			\n\t"\
		"movaps	%%xmm4		,%%xmm6			\n\t"\
		"movaps	%%xmm5		,%%xmm7			\n\t"\
		"addpd	    (%%rdx)	,%%xmm6			\n\t"\
		"addpd	0x10(%%rdx)	,%%xmm7			\n\t"\
		"subpd	    (%%rdx)	,%%xmm4			\n\t"\
		"subpd	0x10(%%rdx)	,%%xmm5			\n\t"\
		"movaps	%%xmm6		,0x40(%%rsi)	\n\t"\
		"movaps	%%xmm7		,0x50(%%rsi)	\n\t"\
		"addpd	%%xmm2		,%%xmm6			\n\t"\
		"addpd	%%xmm3		,%%xmm7			\n\t"\
		"subpd	0x40(%%rsi)	,%%xmm2			\n\t"\
		"subpd	0x50(%%rsi)	,%%xmm3			\n\t"\
		"addpd	    (%%rsi)	,%%xmm6			\n\t"\
		"addpd	0x10(%%rsi)	,%%xmm7			\n\t"\
		"movaps	%%xmm6		,    (%%rsi)	\n\t"\
		"movaps	%%xmm7		,0x10(%%rsi)	\n\t"\
		"subpd	0x80(%%rsi)	,%%xmm6			\n\t"\
		"subpd	0x90(%%rsi)	,%%xmm7			\n\t"\
		"movaps	%%xmm6		,0x80(%%rsi)	\n\t"\
		"movaps	%%xmm7		,0x90(%%rsi)	\n\t"\
		"movaps	%%xmm4		,%%xmm6			\n\t"\
		"movaps	%%xmm5		,%%xmm7			\n\t"\
		"addpd	%%xmm0		,%%xmm5			\n\t"\
		"subpd	%%xmm7		,%%xmm0			\n\t"\
		"addpd	%%xmm1		,%%xmm4			\n\t"\
		"subpd	%%xmm6		,%%xmm1			\n\t"\
		"movaps	%%xmm5		,%%xmm6			\n\t"\
		"movaps	%%xmm1		,%%xmm7			\n\t"\
		"addpd	0xa0(%%rsi)	,%%xmm5			\n\t"\
		"subpd	0xb0(%%rsi)	,%%xmm1			\n\t"\
		"subpd	0xa0(%%rsi)	,%%xmm6			\n\t"\
		"addpd	0xb0(%%rsi)	,%%xmm7			\n\t"\
		"movaps	%%xmm5		,0x20(%%rsi)	\n\t"\
		"movaps	%%xmm1		,0x30(%%rsi)	\n\t"\
		"movaps	%%xmm6		,0xa0(%%rsi)	\n\t"\
		"movaps	%%xmm7		,0xb0(%%rsi)	\n\t"\
		"movaps	%%xmm2		,%%xmm6			\n\t"\
		"movaps	%%xmm3		,%%xmm7			\n\t"\
		"addpd	0xd0(%%rsi)	,%%xmm2			\n\t"\
		"subpd	0xc0(%%rsi)	,%%xmm3			\n\t"\
		"subpd	0xd0(%%rsi)	,%%xmm6			\n\t"\
		"addpd	0xc0(%%rsi)	,%%xmm7			\n\t"\
		"movaps	%%xmm2		,0x40(%%rsi)	\n\t"\
		"movaps	%%xmm3		,0x50(%%rsi)	\n\t"\
		"movaps	%%xmm6		,0xc0(%%rsi)	\n\t"\
		"movaps	%%xmm7		,0xd0(%%rsi)	\n\t"\
		"movaps	%%xmm0		,%%xmm6			\n\t"\
		"movaps	%%xmm4		,%%xmm7			\n\t"\
		"subpd	0xf0(%%rsi)	,%%xmm0			\n\t"\
		"subpd	0xe0(%%rsi)	,%%xmm4			\n\t"\
		"addpd	0xf0(%%rsi)	,%%xmm6			\n\t"\
		"addpd	0xe0(%%rsi)	,%%xmm7			\n\t"\
		"movaps	%%xmm0		,0x60(%%rsi)	\n\t"\
		"movaps	%%xmm4		,0x70(%%rsi)	\n\t"\
		"movaps	%%xmm6		,0xe0(%%rsi)	\n\t"\
		"movaps	%%xmm7		,0xf0(%%rsi)	\n\t"\
		"\n\t"\
		"/*...and now do eight radix-4 transforms, including the internal twiddle factors: */\n\t"\
		"\n\t"\
		"/*...Block 1: t00,t10,t20,t30	*/	\n\t"\
		"movq	%[__add0],%%rax		\n\t"\
		"movslq	%[__p08],%%rbx		\n\t"\
		"movslq	%[__p10],%%rdi	/* edi will store copy of p10 throughout */\n\t"\
		"movq	%[__r00],%%rcx		\n\t"\
		"movq	%[__r10],%%rdx		\n\t"\
		"shlq	$3,%%rbx			\n\t"\
		"shlq	$3,%%rdi			\n\t"\
		"addq	%%rax	,%%rbx	/* add1 = add0+p8 */\n\t"\
		"\n\t"\
		"movaps	     (%%rdx),%%xmm2	\n\t"\
		"movaps	0x200(%%rdx),%%xmm4	\n\t"\
		"movaps	0x010(%%rdx),%%xmm3	\n\t"\
		"movaps	0x210(%%rdx),%%xmm5	\n\t"\
		"movaps	     (%%rcx),%%xmm0	\n\t"\
		"movaps	0x200(%%rcx),%%xmm6	\n\t"\
		"movaps	0x010(%%rcx),%%xmm1	\n\t"\
		"movaps	0x210(%%rcx),%%xmm7	\n\t"\
		"\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm4,%%xmm6		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"subpd	%%xmm5,%%xmm7		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm6,%%xmm4		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"addpd	%%xmm7,%%xmm5		\n\t"\
		"\n\t"\
		"/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c00): */\n\t"\
		"movq	%[__isrt2],%%rsi	\n\t"\
		"subpd	%%xmm4,%%xmm2		\n\t"\
		"subpd	%%xmm7,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm3		\n\t"\
		"subpd	%%xmm6,%%xmm1		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm2,%%xmm4		\n\t"\
		"addpd	%%xmm0,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm5		\n\t"\
		"addpd	%%xmm1,%%xmm6		\n\t"\
		"addq	$0x070,%%rsi	/* c00 */\n\t"\
		"movaps	%%xmm2,     (%%rcx)	\n\t"\
		"movaps	%%xmm0,0x010(%%rdx)	\n\t"\
		"movaps	%%xmm3,0x010(%%rcx)	\n\t"\
		"movaps	%%xmm6,     (%%rdx)	\n\t"\
		"movaps	%%xmm4,%%xmm2		\n\t"\
		"movaps	%%xmm7,%%xmm0		\n\t"\
		"movaps	%%xmm5,%%xmm3		\n\t"\
		"movaps	%%xmm1,%%xmm6		\n\t"\
		"mulpd	     (%%rsi),%%xmm4	\n\t"\
		"mulpd	0x020(%%rsi),%%xmm7	\n\t"\
		"mulpd	     (%%rsi),%%xmm5	\n\t"\
		"mulpd	0x020(%%rsi),%%xmm1	\n\t"\
		"mulpd	0x010(%%rsi),%%xmm2	\n\t"\
		"mulpd	0x030(%%rsi),%%xmm0	\n\t"\
		"mulpd	0x010(%%rsi),%%xmm3	\n\t"\
		"mulpd	0x030(%%rsi),%%xmm6	\n\t"\
		"subpd	%%xmm2,%%xmm5		\n\t"\
		"subpd	%%xmm0,%%xmm1		\n\t"\
		"addpd	%%xmm3,%%xmm4		\n\t"\
		"addpd	%%xmm6,%%xmm7		\n\t"\
		"movaps	%%xmm5,0x10(%%rax)	\n\t"\
		"movaps	%%xmm1,0x10(%%rbx)	\n\t"\
		"movaps	%%xmm4,    (%%rax)	\n\t"\
		"movaps	%%xmm7,    (%%rbx)	\n\t"\
		"addq	%%rdi,%%rax			\n\t"\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	$0x40,%%rsi			\n\t"\
		"movaps	     (%%rcx),%%xmm4	\n\t"\
		"movaps	0x010(%%rdx),%%xmm0	\n\t"\
		"movaps	0x010(%%rcx),%%xmm5	\n\t"\
		"movaps	     (%%rdx),%%xmm6	\n\t"\
		"movaps	%%xmm4,%%xmm2		\n\t"\
		"movaps	%%xmm0,%%xmm1		\n\t"\
		"movaps	%%xmm5,%%xmm3		\n\t"\
		"movaps	%%xmm6,%%xmm7		\n\t"\
		"mulpd	     (%%rsi),%%xmm4	\n\t"\
		"mulpd	0x020(%%rsi),%%xmm0	\n\t"\
		"mulpd	     (%%rsi),%%xmm5	\n\t"\
		"mulpd	0x020(%%rsi),%%xmm6	\n\t"\
		"mulpd	0x010(%%rsi),%%xmm2	\n\t"\
		"mulpd	0x030(%%rsi),%%xmm1	\n\t"\
		"mulpd	0x010(%%rsi),%%xmm3	\n\t"\
		"mulpd	0x030(%%rsi),%%xmm7	\n\t"\
		"subpd	%%xmm2,%%xmm5		\n\t"\
		"subpd	%%xmm1,%%xmm6		\n\t"\
		"addpd	%%xmm3,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm0		\n\t"\
		"movaps	%%xmm5,0x10(%%rax)	\n\t"\
		"movaps	%%xmm6,0x10(%%rbx)	\n\t"\
		"movaps	%%xmm4,    (%%rax)	\n\t"\
		"movaps	%%xmm0,    (%%rbx)	\n\t"\
		"\n\t"\
		"/*...Block 2: t02,t12,t22,t32	*/\n\t"\
		"movq	%[__add0],%%rax		\n\t"\
		"movslq	%[__p01],%%rsi		\n\t"\
		"movslq	%[__p08],%%rbx		\n\t"\
		"shlq	$3,%%rsi			\n\t"\
		"shlq	$3,%%rbx			\n\t"\
		"addq	%%rsi	,%%rax	/* add0 = &a[j1+p1] */\n\t"\
		"addq	%%rax	,%%rbx	/* add1 = add0+p8 */\n\t"\
		"movq	%[__isrt2],%%rsi	\n\t"\
		"addq	$0x220,%%rcx	/* r22 */\n\t"\
		"addq	$0x220,%%rdx	/* r32 */\n\t"\
		"\n\t"\
		"addq	$0x030,%%rsi	/* cc1 */\n\t"\
		"movaps	     (%%rcx),%%xmm4	\n\t"\
		"movaps	0x010(%%rcx),%%xmm5	\n\t"\
		"movaps	     (%%rsi),%%xmm2	\n\t"\
		"movaps	0x010(%%rsi),%%xmm3	\n\t"\
		"movaps	%%xmm4,%%xmm6		\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t"\
		"\n\t"\
		"mulpd	%%xmm2,%%xmm4		\n\t"\
		"mulpd	%%xmm2,%%xmm5		\n\t"\
		"mulpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	     (%%rdx),%%xmm0	\n\t"\
		"mulpd	%%xmm3,%%xmm7		\n\t"\
		"movaps	0x010(%%rdx),%%xmm1	\n\t"\
		"addq	$0x020,%%rsi	/* cc3 */\n\t"\
		"movaps	     (%%rsi),%%xmm2	\n\t"\
		"movaps	0x010(%%rsi),%%xmm3	\n\t"\
		"subpd	%%xmm6,%%xmm5		\n\t"\
		"movaps	%%xmm0,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm4		\n\t"\
		"movaps	%%xmm1,%%xmm7		\n\t"\
		"\n\t"\
		"mulpd	%%xmm2,%%xmm0		\n\t"\
		"mulpd	%%xmm2,%%xmm1		\n\t"\
		"mulpd	%%xmm3,%%xmm6		\n\t"\
		"mulpd	%%xmm3,%%xmm7		\n\t"\
		"subpd	%%xmm6,%%xmm1		\n\t"\
		"addpd	%%xmm7,%%xmm0		\n\t"\
		"\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t"\
		"movaps	%%xmm4,%%xmm6		\n\t"\
		"\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"subpd	%%xmm0,%%xmm6		\n\t"\
		"subpd	%%xmm1,%%xmm7		\n\t"\
		"\n\t"\
		"subq	$0x200,%%rcx	/* r02 */\n\t"\
		"subq	$0x200,%%rdx	/* r12 */\n\t"\
		"subq	$0x040,%%rsi	/* cc0 */\n\t"\
		"movaps	     (%%rdx),%%xmm1	\n\t"\
		"movaps	0x010(%%rdx),%%xmm3	\n\t"\
		"movaps	0x010(%%rsi),%%xmm2	\n\t"\
		"movaps	%%xmm1,%%xmm0		\n\t"\
		"mulpd	%%xmm2,%%xmm1		\n\t"\
		"mulpd	%%xmm3,%%xmm2		\n\t"\
		"mulpd	(%%rsi),%%xmm0		\n\t"\
		"mulpd	(%%rsi),%%xmm3		\n\t"\
		"\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"subpd	%%xmm1,%%xmm3		\n\t"\
		"\n\t"\
		"movaps	     (%%rcx),%%xmm0	\n\t"\
		"movaps	0x010(%%rcx),%%xmm1	\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"\n\t"\
		"/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c01): */\n\t"\
		"addq	$0x260,%%rsi	/* c01 */\n\t"\
		"subpd	%%xmm4,%%xmm2		\n\t"\
		"subpd	%%xmm7,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm3		\n\t"\
		"subpd	%%xmm6,%%xmm1		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm2,%%xmm4		\n\t"\
		"addpd	%%xmm0,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm5		\n\t"\
		"addpd	%%xmm1,%%xmm6		\n\t"\
		"movaps	%%xmm2,     (%%rcx)	\n\t"\
		"movaps	%%xmm0,0x010(%%rdx)	\n\t"\
		"movaps	%%xmm3,0x010(%%rcx)	\n\t"\
		"movaps	%%xmm6,     (%%rdx)	\n\t"\
		"movaps	%%xmm4,%%xmm2		\n\t"\
		"movaps	%%xmm7,%%xmm0		\n\t"\
		"movaps	%%xmm5,%%xmm3		\n\t"\
		"movaps	%%xmm1,%%xmm6		\n\t"\
		"mulpd	     (%%rsi),%%xmm4	\n\t"\
		"mulpd	0x020(%%rsi),%%xmm7	\n\t"\
		"mulpd	     (%%rsi),%%xmm5	\n\t"\
		"mulpd	0x020(%%rsi),%%xmm1	\n\t"\
		"mulpd	0x010(%%rsi),%%xmm2	\n\t"\
		"mulpd	0x030(%%rsi),%%xmm0	\n\t"\
		"mulpd	0x010(%%rsi),%%xmm3	\n\t"\
		"mulpd	0x030(%%rsi),%%xmm6	\n\t"\
		"subpd	%%xmm2,%%xmm5		\n\t"\
		"subpd	%%xmm0,%%xmm1		\n\t"\
		"addpd	%%xmm3,%%xmm4		\n\t"\
		"addpd	%%xmm6,%%xmm7		\n\t"\
		"movaps	%%xmm5,0x10(%%rax)	\n\t"\
		"movaps	%%xmm1,0x10(%%rbx)	\n\t"\
		"movaps	%%xmm4,    (%%rax)	\n\t"\
		"movaps	%%xmm7,    (%%rbx)	\n\t"\
		"addq	%%rdi,%%rax			\n\t"\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	$0x40,%%rsi			\n\t"\
		"movaps	     (%%rcx),%%xmm4	\n\t"\
		"movaps	0x010(%%rdx),%%xmm0	\n\t"\
		"movaps	0x010(%%rcx),%%xmm5	\n\t"\
		"movaps	     (%%rdx),%%xmm6	\n\t"\
		"movaps	%%xmm4,%%xmm2		\n\t"\
		"movaps	%%xmm0,%%xmm1		\n\t"\
		"movaps	%%xmm5,%%xmm3		\n\t"\
		"movaps	%%xmm6,%%xmm7		\n\t"\
		"mulpd	     (%%rsi),%%xmm4	\n\t"\
		"mulpd	0x020(%%rsi),%%xmm0	\n\t"\
		"mulpd	     (%%rsi),%%xmm5	\n\t"\
		"mulpd	0x020(%%rsi),%%xmm6	\n\t"\
		"mulpd	0x010(%%rsi),%%xmm2	\n\t"\
		"mulpd	0x030(%%rsi),%%xmm1	\n\t"\
		"mulpd	0x010(%%rsi),%%xmm3	\n\t"\
		"mulpd	0x030(%%rsi),%%xmm7	\n\t"\
		"subpd	%%xmm2,%%xmm5		\n\t"\
		"subpd	%%xmm1,%%xmm6		\n\t"\
		"addpd	%%xmm3,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm0		\n\t"\
		"movaps	%%xmm5,0x10(%%rax)	\n\t"\
		"movaps	%%xmm6,0x10(%%rbx)	\n\t"\
		"movaps	%%xmm4,    (%%rax)	\n\t"\
		"movaps	%%xmm0,    (%%rbx)	\n\t"\
		"\n\t"\
		"/*...Block 3: t04,t14,t24,t34*/\n\t"\
		"movq	%[__add0],%%rax		\n\t"\
		"movslq	%[__p02],%%rsi		\n\t"\
		"movslq	%[__p08],%%rbx		\n\t"\
		"shlq	$3,%%rsi			\n\t"\
		"shlq	$3,%%rbx			\n\t"\
		"addq	%%rsi	,%%rax	/* add0 = &a[j1+p2] */\n\t"\
		"addq	%%rax	,%%rbx	/* add1 = add0+p8 */\n\t"\
		"movq	%[__isrt2],%%rsi	\n\t"\
		"addq	$0x220,%%rcx	/* r24 */\n\t"\
		"addq	$0x220,%%rdx	/* r34 */\n\t"\
		"addq	$0x010,%%rsi	/* cc0 */\n\t"\
		"\n\t"\
		"movaps	     (%%rcx),%%xmm4	\n\t"\
		"movaps	0x010(%%rcx),%%xmm5	\n\t"\
		"movaps	     (%%rsi),%%xmm2	\n\t"\
		"movaps	0x010(%%rsi),%%xmm3	\n\t"\
		"movaps	%%xmm4,%%xmm6		\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t"\
		"\n\t"\
		"mulpd	%%xmm2,%%xmm4		\n\t"\
		"mulpd	%%xmm2,%%xmm5		\n\t"\
		"mulpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	     (%%rdx),%%xmm0	\n\t"\
		"mulpd	%%xmm3,%%xmm7		\n\t"\
		"movaps	0x010(%%rdx),%%xmm1	\n\t"\
		"subpd	%%xmm6,%%xmm5		\n\t"\
		"movaps	%%xmm0,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm4		\n\t"\
		"movaps	%%xmm1,%%xmm7		\n\t"\
		"\n\t"\
		"mulpd	%%xmm3,%%xmm0		\n\t"\
		"mulpd	%%xmm3,%%xmm1		\n\t"\
		"mulpd	%%xmm2,%%xmm6		\n\t"\
		"mulpd	%%xmm2,%%xmm7		\n\t"\
		"subpd	%%xmm6,%%xmm1		\n\t"\
		"addpd	%%xmm7,%%xmm0		\n\t"\
		"\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t"\
		"movaps	%%xmm4,%%xmm6		\n\t"\
		"\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"subpd	%%xmm0,%%xmm6		\n\t"\
		"subpd	%%xmm1,%%xmm7		\n\t"\
		"\n\t"\
		"subq	$0x200,%%rcx		\n\t"\
		"subq	$0x200,%%rdx		\n\t"\
		"subq	$0x10,%%rsi	/* isrt2 */\n\t"\
		"movaps	     (%%rdx),%%xmm2	\n\t"\
		"movaps	0x010(%%rdx),%%xmm3	\n\t"\
		"movaps	(%%rsi),%%xmm1	/* isrt2 */\n\t"\
		"movaps	%%xmm3,%%xmm0		\n\t"\
		"subpd	%%xmm2,%%xmm3		\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"mulpd	%%xmm1,%%xmm2		\n\t"\
		"mulpd	%%xmm1,%%xmm3		\n\t"\
		"\n\t"\
		"movaps	     (%%rcx),%%xmm0	\n\t"\
		"movaps	0x010(%%rcx),%%xmm1	\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"\n\t"\
		"/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c02): */\n\t"\
		"addq	$0x170,%%rsi	/* c02 */\n\t"\
		"subpd	%%xmm4,%%xmm2		\n\t"\
		"subpd	%%xmm7,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm3		\n\t"\
		"subpd	%%xmm6,%%xmm1		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm2,%%xmm4		\n\t"\
		"addpd	%%xmm0,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm5		\n\t"\
		"addpd	%%xmm1,%%xmm6		\n\t"\
		"movaps	%%xmm2,     (%%rcx)	\n\t"\
		"movaps	%%xmm0,0x010(%%rdx)	\n\t"\
		"movaps	%%xmm3,0x010(%%rcx)	\n\t"\
		"movaps	%%xmm6,     (%%rdx)	\n\t"\
		"movaps	%%xmm4,%%xmm2		\n\t"\
		"movaps	%%xmm7,%%xmm0		\n\t"\
		"movaps	%%xmm5,%%xmm3		\n\t"\
		"movaps	%%xmm1,%%xmm6		\n\t"\
		"mulpd	     (%%rsi),%%xmm4	\n\t"\
		"mulpd	0x020(%%rsi),%%xmm7	\n\t"\
		"mulpd	     (%%rsi),%%xmm5	\n\t"\
		"mulpd	0x020(%%rsi),%%xmm1	\n\t"\
		"mulpd	0x010(%%rsi),%%xmm2	\n\t"\
		"mulpd	0x030(%%rsi),%%xmm0	\n\t"\
		"mulpd	0x010(%%rsi),%%xmm3	\n\t"\
		"mulpd	0x030(%%rsi),%%xmm6	\n\t"\
		"subpd	%%xmm2,%%xmm5		\n\t"\
		"subpd	%%xmm0,%%xmm1		\n\t"\
		"addpd	%%xmm3,%%xmm4		\n\t"\
		"addpd	%%xmm6,%%xmm7		\n\t"\
		"movaps	%%xmm5,0x10(%%rax)	\n\t"\
		"movaps	%%xmm1,0x10(%%rbx)	\n\t"\
		"movaps	%%xmm4,    (%%rax)	\n\t"\
		"movaps	%%xmm7,    (%%rbx)	\n\t"\
		"addq	%%rdi,%%rax			\n\t"\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	$0x40,%%rsi			\n\t"\
		"movaps	     (%%rcx),%%xmm4	\n\t"\
		"movaps	0x010(%%rdx),%%xmm0	\n\t"\
		"movaps	0x010(%%rcx),%%xmm5	\n\t"\
		"movaps	     (%%rdx),%%xmm6	\n\t"\
		"movaps	%%xmm4,%%xmm2		\n\t"\
		"movaps	%%xmm0,%%xmm1		\n\t"\
		"movaps	%%xmm5,%%xmm3		\n\t"\
		"movaps	%%xmm6,%%xmm7		\n\t"\
		"mulpd	     (%%rsi),%%xmm4	\n\t"\
		"mulpd	0x020(%%rsi),%%xmm0	\n\t"\
		"mulpd	     (%%rsi),%%xmm5	\n\t"\
		"mulpd	0x020(%%rsi),%%xmm6	\n\t"\
		"mulpd	0x010(%%rsi),%%xmm2	\n\t"\
		"mulpd	0x030(%%rsi),%%xmm1	\n\t"\
		"mulpd	0x010(%%rsi),%%xmm3	\n\t"\
		"mulpd	0x030(%%rsi),%%xmm7	\n\t"\
		"subpd	%%xmm2,%%xmm5		\n\t"\
		"subpd	%%xmm1,%%xmm6		\n\t"\
		"addpd	%%xmm3,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm0		\n\t"\
		"movaps	%%xmm5,0x10(%%rax)	\n\t"\
		"movaps	%%xmm6,0x10(%%rbx)	\n\t"\
		"movaps	%%xmm4,    (%%rax)	\n\t"\
		"movaps	%%xmm0,    (%%rbx)	\n\t"\
		"\n\t"\
		"/*...Block 4: t06,t16,t26,t36*/\n\t"\
		"movq	%[__add0],%%rax		\n\t"\
		"movslq	%[__p03],%%rsi		\n\t"\
		"movslq	%[__p08],%%rbx		\n\t"\
		"shlq	$3,%%rsi			\n\t"\
		"shlq	$3,%%rbx			\n\t"\
		"addq	%%rsi	,%%rax	/* add0 = &a[j1+p3] */\n\t"\
		"addq	%%rax	,%%rbx	/* add1 = add0+p8 */\n\t"\
		"movq	%[__isrt2],%%rsi	\n\t"\
		"addq	$0x220,%%rcx	/* r26 */\n\t"\
		"addq	$0x220,%%rdx	/* r36 */\n\t"\
		"addq	$0x050,%%rsi	/* cc3 */\n\t"\
		"\n\t"\
		"movaps	     (%%rcx),%%xmm4	\n\t"\
		"movaps	0x010(%%rcx),%%xmm5	\n\t"\
		"movaps	     (%%rsi),%%xmm2	\n\t"\
		"movaps	0x010(%%rsi),%%xmm3	\n\t"\
		"movaps	%%xmm4,%%xmm6		\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t"\
		"\n\t"\
		"mulpd	%%xmm2,%%xmm4		\n\t"\
		"mulpd	%%xmm2,%%xmm5		\n\t"\
		"mulpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	     (%%rdx),%%xmm0	\n\t"\
		"mulpd	%%xmm3,%%xmm7		\n\t"\
		"movaps	0x010(%%rdx),%%xmm1	\n\t"\
		"subq	$0x20,%%rsi			/* cc1 */\n\t"\
		"movaps	     (%%rsi),%%xmm3	/* c32_1 */\n\t"\
		"movaps	0x010(%%rsi),%%xmm2	/* s32_1 */\n\t"\
		"subpd	%%xmm6,%%xmm5		\n\t"\
		"movaps	%%xmm0,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm4		\n\t"\
		"movaps	%%xmm1,%%xmm7		\n\t"\
		"\n\t"\
		"mulpd	%%xmm2,%%xmm0		\n\t"\
		"mulpd	%%xmm2,%%xmm1		\n\t"\
		"mulpd	%%xmm3,%%xmm6		\n\t"\
		"mulpd	%%xmm3,%%xmm7		\n\t"\
		"addpd	%%xmm6,%%xmm1		\n\t"\
		"subpd	%%xmm7,%%xmm0		\n\t"\
		"\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t"\
		"movaps	%%xmm4,%%xmm6		\n\t"\
		"\n\t"\
		"addpd	%%xmm0,%%xmm6		\n\t"\
		"addpd	%%xmm1,%%xmm7		\n\t"\
		"subpd	%%xmm0,%%xmm4		\n\t"\
		"subpd	%%xmm1,%%xmm5		\n\t"\
		"\n\t"\
		"subq	$0x200,%%rcx		\n\t"\
		"subq	$0x200,%%rdx		\n\t"\
		"subq	$0x20,%%rsi	/* cc0 */\n\t"\
		"movaps	     (%%rdx),%%xmm2	\n\t"\
		"movaps	0x010(%%rdx),%%xmm0	\n\t"\
		"movaps	0x010(%%rsi),%%xmm3	\n\t"\
		"movaps	%%xmm2,%%xmm1		\n\t"\
		"mulpd	%%xmm3,%%xmm2		\n\t"\
		"mulpd	%%xmm0,%%xmm3		\n\t"\
		"mulpd	(%%rsi),%%xmm1		\n\t"\
		"mulpd	(%%rsi),%%xmm0		\n\t"\
		"\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"subpd	%%xmm1,%%xmm3		\n\t"\
		"\n\t"\
		"movaps	     (%%rcx),%%xmm0	\n\t"\
		"movaps	0x010(%%rcx),%%xmm1	\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"\n\t"\
		"/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c03): */\n\t"\
		"addq	$0x360,%%rsi	/* c03 */\n\t"\
		"subpd	%%xmm4,%%xmm2		\n\t"\
		"subpd	%%xmm7,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm3		\n\t"\
		"subpd	%%xmm6,%%xmm1		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm2,%%xmm4		\n\t"\
		"addpd	%%xmm0,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm5		\n\t"\
		"addpd	%%xmm1,%%xmm6		\n\t"\
		"movaps	%%xmm2,     (%%rcx)	\n\t"\
		"movaps	%%xmm0,0x010(%%rdx)	\n\t"\
		"movaps	%%xmm3,0x010(%%rcx)	\n\t"\
		"movaps	%%xmm6,     (%%rdx)	\n\t"\
		"movaps	%%xmm4,%%xmm2		\n\t"\
		"movaps	%%xmm7,%%xmm0		\n\t"\
		"movaps	%%xmm5,%%xmm3		\n\t"\
		"movaps	%%xmm1,%%xmm6		\n\t"\
		"mulpd	     (%%rsi),%%xmm4	\n\t"\
		"mulpd	0x020(%%rsi),%%xmm7	\n\t"\
		"mulpd	     (%%rsi),%%xmm5	\n\t"\
		"mulpd	0x020(%%rsi),%%xmm1	\n\t"\
		"mulpd	0x010(%%rsi),%%xmm2	\n\t"\
		"mulpd	0x030(%%rsi),%%xmm0	\n\t"\
		"mulpd	0x010(%%rsi),%%xmm3	\n\t"\
		"mulpd	0x030(%%rsi),%%xmm6	\n\t"\
		"subpd	%%xmm2,%%xmm5		\n\t"\
		"subpd	%%xmm0,%%xmm1		\n\t"\
		"addpd	%%xmm3,%%xmm4		\n\t"\
		"addpd	%%xmm6,%%xmm7		\n\t"\
		"movaps	%%xmm5,0x10(%%rax)	\n\t"\
		"movaps	%%xmm1,0x10(%%rbx)	\n\t"\
		"movaps	%%xmm4,    (%%rax)	\n\t"\
		"movaps	%%xmm7,    (%%rbx)	\n\t"\
		"addq	%%rdi,%%rax			\n\t"\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	$0x40,%%rsi	/* c0B */\n\t"\
		"movaps	     (%%rcx),%%xmm4	\n\t"\
		"movaps	0x010(%%rdx),%%xmm0	\n\t"\
		"movaps	0x010(%%rcx),%%xmm5	\n\t"\
		"movaps	     (%%rdx),%%xmm6	\n\t"\
		"movaps	%%xmm4,%%xmm2		\n\t"\
		"movaps	%%xmm0,%%xmm1		\n\t"\
		"movaps	%%xmm5,%%xmm3		\n\t"\
		"movaps	%%xmm6,%%xmm7		\n\t"\
		"mulpd	     (%%rsi),%%xmm4	\n\t"\
		"mulpd	0x020(%%rsi),%%xmm0	\n\t"\
		"mulpd	     (%%rsi),%%xmm5	\n\t"\
		"mulpd	0x020(%%rsi),%%xmm6	\n\t"\
		"mulpd	0x010(%%rsi),%%xmm2	\n\t"\
		"mulpd	0x030(%%rsi),%%xmm1	\n\t"\
		"mulpd	0x010(%%rsi),%%xmm3	\n\t"\
		"mulpd	0x030(%%rsi),%%xmm7	\n\t"\
		"subpd	%%xmm2,%%xmm5		\n\t"\
		"subpd	%%xmm1,%%xmm6		\n\t"\
		"addpd	%%xmm3,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm0		\n\t"\
		"movaps	%%xmm5,0x10(%%rax)	\n\t"\
		"movaps	%%xmm6,0x10(%%rbx)	\n\t"\
		"movaps	%%xmm4,    (%%rax)	\n\t"\
		"movaps	%%xmm0,    (%%rbx)	\n\t"\
		"\n\t"\
		"/*...Block 5: t08,t18,t28,t38*/\n\t"\
		"movq	%[__add0],%%rax		\n\t"\
		"movslq	%[__p04],%%rsi		\n\t"\
		"movslq	%[__p08],%%rbx		\n\t"\
		"shlq	$3,%%rsi			\n\t"\
		"shlq	$3,%%rbx			\n\t"\
		"addq	%%rsi	,%%rax	/* add0 = &a[j1+p4] */\n\t"\
		"addq	%%rax	,%%rbx	/* add1 = add0+p8 */\n\t"\
		"movq	%[__isrt2],%%rsi	\n\t"\
		"addq	$0x220,%%rcx	/* r28 */\n\t"\
		"addq	$0x220,%%rdx	/* r38 */\n\t"\
		"movaps	(%%rsi),%%xmm2	/* isrt2 */\n\t"\
		"\n\t"\
		"movaps	     (%%rcx),%%xmm4	\n\t"\
		"movaps	0x010(%%rcx),%%xmm5	\n\t"\
		"movaps	     (%%rdx),%%xmm6	\n\t"\
		"movaps	0x010(%%rdx),%%xmm7	\n\t"\
		"subq	$0x200,%%rcx		\n\t"\
		"subq	$0x200,%%rdx		\n\t"\
		"mulpd	%%xmm2,%%xmm4		\n\t"\
		"mulpd	%%xmm2,%%xmm5		\n\t"\
		"mulpd	%%xmm2,%%xmm6		\n\t"\
		"mulpd	%%xmm2,%%xmm7		\n\t"\
		"\n\t"\
		"subpd	%%xmm4,%%xmm5		\n\t"\
		"movaps	     (%%rcx),%%xmm0	\n\t"\
		"subpd	%%xmm7,%%xmm6		\n\t"\
		"movaps	0x010(%%rdx),%%xmm2	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"movaps	0x010(%%rcx),%%xmm3	\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"movaps	     (%%rdx),%%xmm1	\n\t"\
		"addpd	%%xmm5,%%xmm4		\n\t"\
		"addpd	%%xmm6,%%xmm7		\n\t"\
		"\n\t"\
		"subpd	%%xmm6,%%xmm4		\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm5		\n\t"\
		"subpd	%%xmm1,%%xmm3		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm1,%%xmm1		\n\t"\
		"addpd	%%xmm4,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm5,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm1		\n\t"\
		"\n\t"\
		"/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c04): */\n\t"\
		"addq	$0x0f0,%%rsi	/* c04 */\n\t"\
		"subpd	%%xmm4,%%xmm2		\n\t"\
		"subpd	%%xmm7,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm3		\n\t"\
		"subpd	%%xmm6,%%xmm1		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm2,%%xmm4		\n\t"\
		"addpd	%%xmm0,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm5		\n\t"\
		"addpd	%%xmm1,%%xmm6		\n\t"\
		"movaps	%%xmm2,     (%%rcx)	\n\t"\
		"movaps	%%xmm0,0x010(%%rdx)	\n\t"\
		"movaps	%%xmm3,0x010(%%rcx)	\n\t"\
		"movaps	%%xmm6,     (%%rdx)	\n\t"\
		"movaps	%%xmm4,%%xmm2		\n\t"\
		"movaps	%%xmm7,%%xmm0		\n\t"\
		"movaps	%%xmm5,%%xmm3		\n\t"\
		"movaps	%%xmm1,%%xmm6		\n\t"\
		"mulpd	     (%%rsi),%%xmm4	\n\t"\
		"mulpd	0x020(%%rsi),%%xmm7	\n\t"\
		"mulpd	     (%%rsi),%%xmm5	\n\t"\
		"mulpd	0x020(%%rsi),%%xmm1	\n\t"\
		"mulpd	0x010(%%rsi),%%xmm2	\n\t"\
		"mulpd	0x030(%%rsi),%%xmm0	\n\t"\
		"mulpd	0x010(%%rsi),%%xmm3	\n\t"\
		"mulpd	0x030(%%rsi),%%xmm6	\n\t"\
		"subpd	%%xmm2,%%xmm5		\n\t"\
		"subpd	%%xmm0,%%xmm1		\n\t"\
		"addpd	%%xmm3,%%xmm4		\n\t"\
		"addpd	%%xmm6,%%xmm7		\n\t"\
		"movaps	%%xmm5,0x10(%%rax)	\n\t"\
		"movaps	%%xmm1,0x10(%%rbx)	\n\t"\
		"movaps	%%xmm4,    (%%rax)	\n\t"\
		"movaps	%%xmm7,    (%%rbx)	\n\t"\
		"addq	%%rdi,%%rax			\n\t"\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	$0x40,%%rsi		/* c0C */\n\t"\
		"movaps	     (%%rcx),%%xmm4	\n\t"\
		"movaps	0x010(%%rdx),%%xmm0	\n\t"\
		"movaps	0x010(%%rcx),%%xmm5	\n\t"\
		"movaps	     (%%rdx),%%xmm6	\n\t"\
		"movaps	%%xmm4,%%xmm2		\n\t"\
		"movaps	%%xmm0,%%xmm1		\n\t"\
		"movaps	%%xmm5,%%xmm3		\n\t"\
		"movaps	%%xmm6,%%xmm7		\n\t"\
		"mulpd	     (%%rsi),%%xmm4	\n\t"\
		"mulpd	0x020(%%rsi),%%xmm0	\n\t"\
		"mulpd	     (%%rsi),%%xmm5	\n\t"\
		"mulpd	0x020(%%rsi),%%xmm6	\n\t"\
		"mulpd	0x010(%%rsi),%%xmm2	\n\t"\
		"mulpd	0x030(%%rsi),%%xmm1	\n\t"\
		"mulpd	0x010(%%rsi),%%xmm3	\n\t"\
		"mulpd	0x030(%%rsi),%%xmm7	\n\t"\
		"subpd	%%xmm2,%%xmm5		\n\t"\
		"subpd	%%xmm1,%%xmm6		\n\t"\
		"addpd	%%xmm3,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm0		\n\t"\
		"movaps	%%xmm5,0x10(%%rax)	\n\t"\
		"movaps	%%xmm6,0x10(%%rbx)	\n\t"\
		"movaps	%%xmm4,    (%%rax)	\n\t"\
		"movaps	%%xmm0,    (%%rbx)	\n\t"\
		"\n\t"\
		"/*...Block 6: t0A,t1A,t2A,t3A*/\n\t"\
		"movq	%[__add0],%%rax		\n\t"\
		"movslq	%[__p05],%%rsi		\n\t"\
		"movslq	%[__p08],%%rbx		\n\t"\
		"shlq	$3,%%rsi			\n\t"\
		"shlq	$3,%%rbx			\n\t"\
		"addq	%%rsi	,%%rax	/* add0 = &a[j1+p5] */\n\t"\
		"addq	%%rax	,%%rbx	/* add1 = add0+p8 */\n\t"\
		"movq	%[__isrt2],%%rsi	\n\t"\
		"addq	$0x220,%%rcx	/* r2A */\n\t"\
		"addq	$0x220,%%rdx	/* r3A */\n\t"\
		"addq	$0x050,%%rsi	/* cc3 */\n\t"\
		"\n\t"\
		"movaps	     (%%rcx),%%xmm4	\n\t"\
		"movaps	0x010(%%rcx),%%xmm5	\n\t"\
		"movaps	     (%%rsi),%%xmm3	\n\t"\
		"movaps	0x010(%%rsi),%%xmm2	\n\t"\
		"movaps	%%xmm4,%%xmm6		\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t"\
		"\n\t"\
		"mulpd	%%xmm2,%%xmm4		\n\t"\
		"mulpd	%%xmm2,%%xmm5		\n\t"\
		"mulpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	     (%%rdx),%%xmm0	\n\t"\
		"mulpd	%%xmm3,%%xmm7		\n\t"\
		"movaps	0x010(%%rdx),%%xmm1	/* t3B */\n\t"\
		"subq	$0x20,%%rsi		/* cc1 */\n\t"\
		"movaps	     (%%rsi),%%xmm2	\n\t"\
		"movaps	0x010(%%rsi),%%xmm3	\n\t"\
		"subpd	%%xmm6,%%xmm5		\n\t"\
		"movaps	%%xmm0,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm4		\n\t"\
		"movaps	%%xmm1,%%xmm7		\n\t"\
		"\n\t"\
		"mulpd	%%xmm2,%%xmm0		\n\t"\
		"mulpd	%%xmm2,%%xmm1		\n\t"\
		"mulpd	%%xmm3,%%xmm6		\n\t"\
		"mulpd	%%xmm3,%%xmm7		\n\t"\
		"addpd	%%xmm6,%%xmm1		\n\t"\
		"subpd	%%xmm7,%%xmm0		\n\t"\
		"\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t"\
		"movaps	%%xmm4,%%xmm6		\n\t"\
		"\n\t"\
		"addpd	%%xmm0,%%xmm6		\n\t"\
		"addpd	%%xmm1,%%xmm7		\n\t"\
		"subpd	%%xmm0,%%xmm4		\n\t"\
		"subpd	%%xmm1,%%xmm5		\n\t"\
		"\n\t"\
		"subq	$0x200,%%rcx		\n\t"\
		"subq	$0x200,%%rdx		\n\t"\
		"subq	$0x20,%%rsi		/* cc0 */\n\t"\
		"movaps	     (%%rdx),%%xmm0	\n\t"\
		"movaps	0x010(%%rdx),%%xmm2	\n\t"\
		"movaps	0x010(%%rsi),%%xmm1	\n\t"\
		"movaps	%%xmm0,%%xmm3		\n\t"\
		"mulpd	%%xmm1,%%xmm0		\n\t"\
		"mulpd	%%xmm2,%%xmm1		\n\t"\
		"mulpd	(%%rsi),%%xmm3		\n\t"\
		"mulpd	(%%rsi),%%xmm2		\n\t"\
		"\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"addpd	%%xmm3,%%xmm1		\n\t"\
		"\n\t"\
		"movaps	     (%%rcx),%%xmm2	\n\t"\
		"movaps	0x010(%%rcx),%%xmm3	\n\t"\
		"subpd	%%xmm0,%%xmm2		\n\t"\
		"subpd	%%xmm1,%%xmm3		\n\t"\
		"addpd	%%xmm0,%%xmm0		\n\t"\
		"addpd	%%xmm1,%%xmm1		\n\t"\
		"addpd	%%xmm2,%%xmm0		\n\t"\
		"addpd	%%xmm3,%%xmm1		\n\t"\
		"\n\t"\
		"/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c05): */\n\t"\
		"addq	$0x2e0,%%rsi	/* c05 */\n\t"\
		"subpd	%%xmm4,%%xmm2		\n\t"\
		"subpd	%%xmm7,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm3		\n\t"\
		"subpd	%%xmm6,%%xmm1		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm2,%%xmm4		\n\t"\
		"addpd	%%xmm0,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm5		\n\t"\
		"addpd	%%xmm1,%%xmm6		\n\t"\
		"movaps	%%xmm2,     (%%rcx)	\n\t"\
		"movaps	%%xmm0,0x010(%%rdx)	\n\t"\
		"movaps	%%xmm3,0x010(%%rcx)	\n\t"\
		"movaps	%%xmm6,     (%%rdx)	\n\t"\
		"movaps	%%xmm4,%%xmm2		\n\t"\
		"movaps	%%xmm7,%%xmm0		\n\t"\
		"movaps	%%xmm5,%%xmm3		\n\t"\
		"movaps	%%xmm1,%%xmm6		\n\t"\
		"mulpd	     (%%rsi),%%xmm4	\n\t"\
		"mulpd	0x020(%%rsi),%%xmm7	\n\t"\
		"mulpd	     (%%rsi),%%xmm5	\n\t"\
		"mulpd	0x020(%%rsi),%%xmm1	\n\t"\
		"mulpd	0x010(%%rsi),%%xmm2	\n\t"\
		"mulpd	0x030(%%rsi),%%xmm0	\n\t"\
		"mulpd	0x010(%%rsi),%%xmm3	\n\t"\
		"mulpd	0x030(%%rsi),%%xmm6	\n\t"\
		"subpd	%%xmm2,%%xmm5		\n\t"\
		"subpd	%%xmm0,%%xmm1		\n\t"\
		"addpd	%%xmm3,%%xmm4		\n\t"\
		"addpd	%%xmm6,%%xmm7		\n\t"\
		"movaps	%%xmm5,0x10(%%rax)	\n\t"\
		"movaps	%%xmm1,0x10(%%rbx)	\n\t"\
		"movaps	%%xmm4,    (%%rax)	\n\t"\
		"movaps	%%xmm7,    (%%rbx)	\n\t"\
		"addq	%%rdi,%%rax			\n\t"\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	$0x40,%%rsi		/* c0D */\n\t"\
		"movaps	     (%%rcx),%%xmm4	\n\t"\
		"movaps	0x010(%%rdx),%%xmm0	\n\t"\
		"movaps	0x010(%%rcx),%%xmm5	\n\t"\
		"movaps	     (%%rdx),%%xmm6	\n\t"\
		"movaps	%%xmm4,%%xmm2		\n\t"\
		"movaps	%%xmm0,%%xmm1		\n\t"\
		"movaps	%%xmm5,%%xmm3		\n\t"\
		"movaps	%%xmm6,%%xmm7		\n\t"\
		"mulpd	     (%%rsi),%%xmm4	\n\t"\
		"mulpd	0x020(%%rsi),%%xmm0	\n\t"\
		"mulpd	     (%%rsi),%%xmm5	\n\t"\
		"mulpd	0x020(%%rsi),%%xmm6	\n\t"\
		"mulpd	0x010(%%rsi),%%xmm2	\n\t"\
		"mulpd	0x030(%%rsi),%%xmm1	\n\t"\
		"mulpd	0x010(%%rsi),%%xmm3	\n\t"\
		"mulpd	0x030(%%rsi),%%xmm7	\n\t"\
		"subpd	%%xmm2,%%xmm5		\n\t"\
		"subpd	%%xmm1,%%xmm6		\n\t"\
		"addpd	%%xmm3,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm0		\n\t"\
		"movaps	%%xmm5,0x10(%%rax)	\n\t"\
		"movaps	%%xmm6,0x10(%%rbx)	\n\t"\
		"movaps	%%xmm4,    (%%rax)	\n\t"\
		"movaps	%%xmm0,    (%%rbx)	\n\t"\
		"\n\t"\
		"/*...Block 7: t0C,t1C,t2C,t3C*/\n\t"\
		"movq	%[__add0],%%rax		\n\t"\
		"movslq	%[__p06],%%rsi		\n\t"\
		"movslq	%[__p08],%%rbx		\n\t"\
		"shlq	$3,%%rsi			\n\t"\
		"shlq	$3,%%rbx			\n\t"\
		"addq	%%rsi	,%%rax	/* add0 = &a[j1+p6] */\n\t"\
		"addq	%%rax	,%%rbx	/* add1 = add0+p8 */\n\t"\
		"movq	%[__isrt2],%%rsi	\n\t"\
		"addq	$0x220,%%rcx	/* r2C */\n\t"\
		"addq	$0x220,%%rdx	/* r3C */\n\t"\
		"addq	$0x010,%%rsi	/* cc0 */\n\t"\
		"\n\t"\
		"movaps	     (%%rcx),%%xmm4	\n\t"\
		"movaps	0x010(%%rcx),%%xmm5	\n\t"\
		"movaps	     (%%rsi),%%xmm3	\n\t"\
		"movaps	0x010(%%rsi),%%xmm2	\n\t"\
		"movaps	%%xmm4,%%xmm6		\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t"\
		"\n\t"\
		"mulpd	%%xmm2,%%xmm4	\n\t"\
		"mulpd	%%xmm2,%%xmm5	\n\t"\
		"mulpd	%%xmm3,%%xmm6	\n\t"\
		"movaps	     (%%rdx),%%xmm0	\n\t"\
		"mulpd	%%xmm3,%%xmm7	\n\t"\
		"movaps	0x010(%%rdx),%%xmm1	\n\t"\
		"subpd	%%xmm6,%%xmm5		\n\t"\
		"movaps	%%xmm0,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm4		\n\t"\
		"movaps	%%xmm1,%%xmm7		\n\t"\
		"\n\t"\
		"mulpd	%%xmm3,%%xmm0	\n\t"\
		"mulpd	%%xmm3,%%xmm1	\n\t"\
		"mulpd	%%xmm2,%%xmm6	\n\t"\
		"mulpd	%%xmm2,%%xmm7	\n\t"\
		"subpd	%%xmm6,%%xmm1		\n\t"\
		"addpd	%%xmm7,%%xmm0		\n\t"\
		"\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t"\
		"movaps	%%xmm4,%%xmm6		\n\t"\
		"\n\t"\
		"addpd	%%xmm0,%%xmm6		\n\t"\
		"addpd	%%xmm1,%%xmm7		\n\t"\
		"subpd	%%xmm0,%%xmm4		\n\t"\
		"subpd	%%xmm1,%%xmm5		\n\t"\
		"\n\t"\
		"subq	$0x200,%%rcx		\n\t"\
		"subq	$0x200,%%rdx		\n\t"\
		"subq	$0x10,%%rsi	/* isrt2 */\n\t"\
		"movaps	     (%%rdx),%%xmm0	\n\t"\
		"movaps	0x010(%%rdx),%%xmm1	\n\t"\
		"movaps	(%%rsi),%%xmm3	/* isrt2 */\n\t"\
		"movaps	%%xmm0,%%xmm2		\n\t"\
		"subpd	%%xmm1,%%xmm0		\n\t"\
		"addpd	%%xmm2,%%xmm1		\n\t"\
		"mulpd	%%xmm3,%%xmm0	\n\t"\
		"mulpd	%%xmm3,%%xmm1	\n\t"\
		"\n\t"\
		"movaps	     (%%rcx),%%xmm2	\n\t"\
		"movaps	0x010(%%rcx),%%xmm3	\n\t"\
		"subpd	%%xmm0,%%xmm2		\n\t"\
		"subpd	%%xmm1,%%xmm3		\n\t"\
		"addpd	%%xmm0,%%xmm0		\n\t"\
		"addpd	%%xmm1,%%xmm1		\n\t"\
		"addpd	%%xmm2,%%xmm0		\n\t"\
		"addpd	%%xmm3,%%xmm1		\n\t"\
		"\n\t"\
		"/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c06): */\n\t"\
		"addq	$0x1f0,%%rsi	/* c06 */\n\t"\
		"subpd	%%xmm4,%%xmm2		\n\t"\
		"subpd	%%xmm7,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm3		\n\t"\
		"subpd	%%xmm6,%%xmm1		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm2,%%xmm4		\n\t"\
		"addpd	%%xmm0,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm5		\n\t"\
		"addpd	%%xmm1,%%xmm6		\n\t"\
		"movaps	%%xmm2,     (%%rcx)	\n\t"\
		"movaps	%%xmm0,0x010(%%rdx)	\n\t"\
		"movaps	%%xmm3,0x010(%%rcx)	\n\t"\
		"movaps	%%xmm6,     (%%rdx)	\n\t"\
		"movaps	%%xmm4,%%xmm2		\n\t"\
		"movaps	%%xmm7,%%xmm0		\n\t"\
		"movaps	%%xmm5,%%xmm3		\n\t"\
		"movaps	%%xmm1,%%xmm6		\n\t"\
		"mulpd	     (%%rsi),%%xmm4	\n\t"\
		"mulpd	0x020(%%rsi),%%xmm7	\n\t"\
		"mulpd	     (%%rsi),%%xmm5	\n\t"\
		"mulpd	0x020(%%rsi),%%xmm1	\n\t"\
		"mulpd	0x010(%%rsi),%%xmm2	\n\t"\
		"mulpd	0x030(%%rsi),%%xmm0	\n\t"\
		"mulpd	0x010(%%rsi),%%xmm3	\n\t"\
		"mulpd	0x030(%%rsi),%%xmm6	\n\t"\
		"subpd	%%xmm2,%%xmm5		\n\t"\
		"subpd	%%xmm0,%%xmm1		\n\t"\
		"addpd	%%xmm3,%%xmm4		\n\t"\
		"addpd	%%xmm6,%%xmm7		\n\t"\
		"movaps	%%xmm5,0x10(%%rax)	\n\t"\
		"movaps	%%xmm1,0x10(%%rbx)	\n\t"\
		"movaps	%%xmm4,    (%%rax)	\n\t"\
		"movaps	%%xmm7,    (%%rbx)	\n\t"\
		"addq	%%rdi,%%rax			\n\t"\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	$0x40,%%rsi		/* c0E */\n\t"\
		"movaps	     (%%rcx),%%xmm4	\n\t"\
		"movaps	0x010(%%rdx),%%xmm0	\n\t"\
		"movaps	0x010(%%rcx),%%xmm5	\n\t"\
		"movaps	     (%%rdx),%%xmm6	\n\t"\
		"movaps	%%xmm4,%%xmm2		\n\t"\
		"movaps	%%xmm0,%%xmm1		\n\t"\
		"movaps	%%xmm5,%%xmm3		\n\t"\
		"movaps	%%xmm6,%%xmm7		\n\t"\
		"mulpd	     (%%rsi),%%xmm4	\n\t"\
		"mulpd	0x020(%%rsi),%%xmm0	\n\t"\
		"mulpd	     (%%rsi),%%xmm5	\n\t"\
		"mulpd	0x020(%%rsi),%%xmm6	\n\t"\
		"mulpd	0x010(%%rsi),%%xmm2	\n\t"\
		"mulpd	0x030(%%rsi),%%xmm1	\n\t"\
		"mulpd	0x010(%%rsi),%%xmm3	\n\t"\
		"mulpd	0x030(%%rsi),%%xmm7	\n\t"\
		"subpd	%%xmm2,%%xmm5		\n\t"\
		"subpd	%%xmm1,%%xmm6		\n\t"\
		"addpd	%%xmm3,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm0		\n\t"\
		"movaps	%%xmm5,0x10(%%rax)	\n\t"\
		"movaps	%%xmm6,0x10(%%rbx)	\n\t"\
		"movaps	%%xmm4,    (%%rax)	\n\t"\
		"movaps	%%xmm0,    (%%rbx)	\n\t"\
		"\n\t"\
		"/*...Block 8: t0E,t1E,t2E,t3E*/\n\t"\
		"movq	%[__add0],%%rax		\n\t"\
		"movslq	%[__p07],%%rsi		\n\t"\
		"movslq	%[__p08],%%rbx		\n\t"\
		"shlq	$3,%%rsi			\n\t"\
		"shlq	$3,%%rbx			\n\t"\
		"addq	%%rsi	,%%rax	/* add0 = &a[j1+p7] */\n\t"\
		"addq	%%rax	,%%rbx	/* add1 = add0+p8 */\n\t"\
		"movq	%[__isrt2],%%rsi	\n\t"\
		"addq	$0x220,%%rcx	/* r2E */\n\t"\
		"addq	$0x220,%%rdx	/* r3E */\n\t"\
		"addq	$0x030,%%rsi	/* cc1 */\n\t"\
		"\n\t"\
		"movaps	     (%%rcx),%%xmm4	\n\t"\
		"movaps	0x010(%%rcx),%%xmm5	\n\t"\
		"movaps	     (%%rsi),%%xmm3	\n\t"\
		"movaps	0x010(%%rsi),%%xmm2	\n\t"\
		"movaps	%%xmm4,%%xmm6		\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t"\
		"\n\t"\
		"mulpd	%%xmm2,%%xmm4		\n\t"\
		"mulpd	%%xmm2,%%xmm5		\n\t"\
		"mulpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	     (%%rdx),%%xmm0	\n\t"\
		"mulpd	%%xmm3,%%xmm7		\n\t"\
		"movaps	0x010(%%rdx),%%xmm1	\n\t"\
		"addq	$0x20,%%rsi		/* cc3 */\n\t"\
		"movaps	     (%%rsi),%%xmm3	/* c32_3 */\n\t"\
		"movaps	0x010(%%rsi),%%xmm2	/* s32_3 */\n\t"\
		"subpd	%%xmm6,%%xmm5		\n\t"\
		"movaps	%%xmm0,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm4		\n\t"\
		"movaps	%%xmm1,%%xmm7		\n\t"\
		"\n\t"\
		"mulpd	%%xmm2,%%xmm0		\n\t"\
		"mulpd	%%xmm2,%%xmm1		\n\t"\
		"mulpd	%%xmm3,%%xmm6		\n\t"\
		"mulpd	%%xmm3,%%xmm7		\n\t"\
		"subpd	%%xmm6,%%xmm1		\n\t"\
		"addpd	%%xmm7,%%xmm0		\n\t"\
		"\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t"\
		"movaps	%%xmm4,%%xmm6		\n\t"\
		"\n\t"\
		"subpd	%%xmm0,%%xmm4		\n\t"\
		"subpd	%%xmm1,%%xmm5		\n\t"\
		"addpd	%%xmm0,%%xmm6		\n\t"\
		"addpd	%%xmm1,%%xmm7		\n\t"\
		"\n\t"\
		"subq	$0x200,%%rcx		\n\t"\
		"subq	$0x200,%%rdx		\n\t"\
		"subq	$0x40,%%rsi	 	/* cc0 */\n\t"\
		"movaps	     (%%rdx),%%xmm3	\n\t"\
		"movaps	0x010(%%rdx),%%xmm1	\n\t"\
		"movaps	0x010(%%rsi),%%xmm2	\n\t"\
		"movaps	%%xmm3,%%xmm0		\n\t"\
		"mulpd	%%xmm2,%%xmm3		\n\t"\
		"mulpd	%%xmm1,%%xmm2		\n\t"\
		"mulpd	(%%rsi),%%xmm0		\n\t"\
		"mulpd	(%%rsi),%%xmm1		\n\t"\
		"\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"addpd	%%xmm3,%%xmm1		\n\t"\
		"\n\t"\
		"movaps	     (%%rcx),%%xmm2	\n\t"\
		"movaps	0x010(%%rcx),%%xmm3	\n\t"\
		"subpd	%%xmm0,%%xmm2		\n\t"\
		"subpd	%%xmm1,%%xmm3		\n\t"\
		"addpd	%%xmm0,%%xmm0		\n\t"\
		"addpd	%%xmm1,%%xmm1		\n\t"\
		"addpd	%%xmm2,%%xmm0		\n\t"\
		"addpd	%%xmm3,%%xmm1		\n\t"\
		"/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c07): */\n\t"\
		"addq	$0x3e0,%%rsi	/* c07 */\n\t"\
		"subpd	%%xmm4,%%xmm2		\n\t"\
		"subpd	%%xmm7,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm3		\n\t"\
		"subpd	%%xmm6,%%xmm1		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm2,%%xmm4		\n\t"\
		"addpd	%%xmm0,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm5		\n\t"\
		"addpd	%%xmm1,%%xmm6		\n\t"\
		"movaps	%%xmm2,     (%%rcx)	\n\t"\
		"movaps	%%xmm0,0x010(%%rdx)	\n\t"\
		"movaps	%%xmm3,0x010(%%rcx)	\n\t"\
		"movaps	%%xmm6,     (%%rdx)	\n\t"\
		"movaps	%%xmm4,%%xmm2		\n\t"\
		"movaps	%%xmm7,%%xmm0		\n\t"\
		"movaps	%%xmm5,%%xmm3		\n\t"\
		"movaps	%%xmm1,%%xmm6		\n\t"\
		"mulpd	     (%%rsi),%%xmm4	\n\t"\
		"mulpd	0x020(%%rsi),%%xmm7	\n\t"\
		"mulpd	     (%%rsi),%%xmm5	\n\t"\
		"mulpd	0x020(%%rsi),%%xmm1	\n\t"\
		"mulpd	0x010(%%rsi),%%xmm2	\n\t"\
		"mulpd	0x030(%%rsi),%%xmm0	\n\t"\
		"mulpd	0x010(%%rsi),%%xmm3	\n\t"\
		"mulpd	0x030(%%rsi),%%xmm6	\n\t"\
		"subpd	%%xmm2,%%xmm5		\n\t"\
		"subpd	%%xmm0,%%xmm1		\n\t"\
		"addpd	%%xmm3,%%xmm4		\n\t"\
		"addpd	%%xmm6,%%xmm7		\n\t"\
		"movaps	%%xmm5,0x10(%%rax)	\n\t"\
		"movaps	%%xmm1,0x10(%%rbx)	\n\t"\
		"movaps	%%xmm4,    (%%rax)	\n\t"\
		"movaps	%%xmm7,    (%%rbx)	\n\t"\
		"addq	%%rdi,%%rax			\n\t"\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	$0x40,%%rsi		/* c0F */\n\t"\
		"movaps	     (%%rcx),%%xmm4	\n\t"\
		"movaps	0x010(%%rdx),%%xmm0	\n\t"\
		"movaps	0x010(%%rcx),%%xmm5	\n\t"\
		"movaps	     (%%rdx),%%xmm6	\n\t"\
		"movaps	%%xmm4,%%xmm2		\n\t"\
		"movaps	%%xmm0,%%xmm1		\n\t"\
		"movaps	%%xmm5,%%xmm3		\n\t"\
		"movaps	%%xmm6,%%xmm7		\n\t"\
		"mulpd	     (%%rsi),%%xmm4	\n\t"\
		"mulpd	0x020(%%rsi),%%xmm0	\n\t"\
		"mulpd	     (%%rsi),%%xmm5	\n\t"\
		"mulpd	0x020(%%rsi),%%xmm6	\n\t"\
		"mulpd	0x010(%%rsi),%%xmm2	\n\t"\
		"mulpd	0x030(%%rsi),%%xmm1	\n\t"\
		"mulpd	0x010(%%rsi),%%xmm3	\n\t"\
		"mulpd	0x030(%%rsi),%%xmm7	\n\t"\
		"subpd	%%xmm2,%%xmm5		\n\t"\
		"subpd	%%xmm1,%%xmm6		\n\t"\
		"addpd	%%xmm3,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm0		\n\t"\
		"movaps	%%xmm5,0x10(%%rax)	\n\t"\
		"movaps	%%xmm6,0x10(%%rbx)	\n\t"\
		"movaps	%%xmm4,    (%%rax)	\n\t"\
		"movaps	%%xmm0,    (%%rbx)	\n\t"\
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
		: "rax","rbx","rcx","rdx","rdi","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
	);\
	}

#endif	/* radix32_dif_dit_pass_gcc_h_included */

