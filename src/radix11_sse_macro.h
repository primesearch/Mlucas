/*******************************************************************************
*                                                                              *
*   (C) 1997-2019 by Ernst W. Mayer.                                           *
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
#ifndef radix11_sse_macro_h_included
#define radix11_sse_macro_h_included

#include "sse2_macro_gcc64.h"

#ifdef USE_ARM_V8_SIMD

	// Cost: [30 ADD, 120 MUL/FMA, 57 LDP/STP], better than general target of 1 vector-memref per vec_dbl arithmetic op.
	#define SSE2_RADIX_11_DFT(Xin0,Xi_off, Xcc, Xout0,Xo_off)\
	{\
	__asm__ volatile (		/*** Rcol does Imaginary Parts: ***/\
		/* v31 = two, v0 = one, but discard the latter: */\
		"ldr	x11,%[__cc]		\n\t	ldp	q31,q0,[x11,#-0x20]	\n\t"\
		/* Compute &i[0-A] from 32-bit offsets: */\
		"ldr	x15,%[in0]		\n\t	ldr	x16,%[i_off]	\n\t"\
		"ldr	w0,[x16]		\n\t	add	x0,x15,x0	\n\t"/* I0 */\
		"ldp	w1 ,w2 ,[x16,#0x04]	\n\t"\
		"ldp	w3 ,w4 ,[x16,#0x0c]	\n\t"\
		"ldp	w5 ,w6 ,[x16,#0x14]	\n\t"\
		"ldp	w7 ,w8 ,[x16,#0x1c]	\n\t"\
		"ldp	w9 ,w10,[x16,#0x24]	\n\t"\
	/* ...Do the + parts of the opening wave of radix-2 butterflies: */\
		"add	x1 ,x15,x1			\n\t	ldp	q1 ,q6 ,[x1]	\n\t"/* x1 */\
		"add	x10,x15,x10			\n\t	ldp	q15,q20,[x10]	\n\t"/* xA */\
		"fadd	v1.2d,v1.2d,v15.2d	\n\t	fadd	v6.2d ,v6.2d ,v20.2d	\n\t"/* t1 = A1 + Aa */\
		"add	x2 ,x15,x2			\n\t	ldp	q2 ,q7 ,[x2]	\n\t"/* x2 */\
		"add	x9 ,x15,x9			\n\t	ldp	q14,q19,[x9]	\n\t"/* x9 */\
		"fadd	v2.2d,v2.2d,v14.2d	\n\t	fadd	v7.2d ,v7.2d ,v19.2d	\n\t"/* t2 = A2 + A9 */\
		"add	x3 ,x15,x3			\n\t	ldp	q3 ,q8 ,[x3]	\n\t"/* x3 */\
		"add	x8 ,x15,x8			\n\t	ldp	q13,q18,[x8]	\n\t"/* x8 */\
		"fadd	v3.2d,v3.2d,v13.2d	\n\t	fadd	v8.2d ,v8.2d ,v18.2d	\n\t"/* t3 = A3 + A8 */\
		"add	x4 ,x15,x4			\n\t	ldp	q4 ,q9 ,[x4]	\n\t"/* x4 */\
		"add	x7 ,x15,x7			\n\t	ldp	q12,q17,[x7]	\n\t"/* x7 */\
		"fadd	v4.2d,v4.2d,v12.2d	\n\t	fadd	v9.2d ,v9.2d ,v17.2d	\n\t"/* t4 = A4 + A7 */\
		"add	x5 ,x15,x5			\n\t	ldp	q5 ,q10,[x5]	\n\t"/* x5 */\
		"add	x6 ,x15,x6			\n\t	ldp	q11,q16,[x6]	\n\t"/* x6 */\
		"fadd	v5.2d,v5.2d,v11.2d	\n\t	fadd	v10.2d,v10.2d,v16.2d	\n\t"/* t5 = A5 + A6 */\
	/* Write resulting lower-half outputs t1-5 back to memory before ensuing FMA-block to compute - parts overwrites each + result: */\
		"stp	q1 ,q6 ,[x1]	\n\t	fmls	v1.2d,v31.2d,v15.2d	\n\t	fmls	v6.2d ,v31.2d,v20.2d	\n\t"/* ta = A1 - Aa */\
		"stp	q2 ,q7 ,[x2]	\n\t	fmls	v2.2d,v31.2d,v14.2d	\n\t	fmls	v7.2d ,v31.2d,v19.2d	\n\t"/* t9 = A2 - A9 */\
		"stp	q3 ,q8 ,[x3]	\n\t	fmls	v3.2d,v31.2d,v13.2d	\n\t	fmls	v8.2d ,v31.2d,v18.2d	\n\t"/* t8 = A3 - A8 */\
		"stp	q4 ,q9 ,[x4]	\n\t	fmls	v4.2d,v31.2d,v12.2d	\n\t	fmls	v9.2d ,v31.2d,v17.2d	\n\t"/* t7 = A4 - A7 */\
		"stp	q5 ,q10,[x5]	\n\t	fmls	v5.2d,v31.2d,v11.2d	\n\t	fmls	v10.2d,v31.2d,v16.2d	\n\t"/* t6 = A5 - A6 */\
	/* ...And write the upper-half outputs back to memory to free up registers for the 5x5 sine-term-subconvo computation: */\
	/*	"stp	q1 ,q6 ,[x10]	\n\t"// tar,i */\
		"stp	q2 ,q7 ,[x9 ]	\n\t"/* t9r,i */\
		"stp	q3 ,q8 ,[x8 ]	\n\t"/* t8r,i */\
		"stp	q4 ,q9 ,[x7 ]	\n\t"/* t7r,i */\
		"stp	q5 ,q10,[x6 ]	\n\t"/* t6r,i */\
/*
	Sr1 =      ss1*tar+ss2*t9r+ss3*t8r+ss4*t7r+ss5*t6r;		Si1 =      ss1*tai+ss2*t9i+ss3*t8i+ss4*t7i+ss5*t6i;	// S1
	Sr2 =      ss2*tar+ss4*t9r-ss5*t8r-ss3*t7r-ss1*t6r;		Si2 =      ss2*tai+ss4*t9i-ss5*t8i-ss3*t7i-ss1*t6i;	// S2
	Sr3 =      ss3*tar-ss5*t9r-ss2*t8r+ss1*t7r+ss4*t6r;		Si3 =      ss3*tai-ss5*t9i-ss2*t8i+ss1*t7i+ss4*t6i;	// S3
	Sr4 =      ss4*tar-ss3*t9r+ss1*t8r+ss5*t7r-ss2*t6r;		Si4 =      ss4*tai-ss3*t9i+ss1*t8i+ss5*t7i-ss2*t6i;	// S4
	Sr5 =      ss5*tar-ss1*t9r+ss4*t8r-ss2*t7r+ss3*t6r;		Si5 =      ss5*tai-ss1*t9i+ss4*t8i-ss2*t7i+ss3*t6i;	// S5

	In a 16-reg FMA model, can keep the 10 S-terms and all 5 shared sincos data in registers, but that
	requires us to reload one of the 2 repeated t*-mults (e.g. tai in the 1st block) 5 times per block.
	OTOH if we keep both t*r,i-mults and 4 of the 5 ss-mults in-reg, we just need 2 loads-from-mem per block:
*/\
	"add	x11,x11,#0x50	\n\t"/* Incr trig ptr to point to ss1 */\
		"ld1r	{v0.2d},[x11]		\n\t"/* ss1 */\
		"ldp	q11,q12,[x11,#0x10]	\n\t"/* ss2,ss3 */\
		"ldp	q13,q14,[x11,#0x30]	\n\t"/* ss4,ss5 */\
		/* S1 = ta already in v1,6; save on reg-copies by taking advantage of opening 3-operand MULs */\
		"fmul	v2.2d,v1.2d,v11.2d	\n\t	fmul	v7.2d ,v6.2d ,v11.2d	\n\t"/* S2  = ss2*ta; */\
		"fmul	v3.2d,v1.2d,v12.2d	\n\t	fmul	v8.2d ,v6.2d ,v12.2d	\n\t"/* S3  = ss3*ta; */\
		"fmul	v4.2d,v1.2d,v13.2d	\n\t	fmul	v9.2d ,v6.2d ,v13.2d	\n\t"/* S4  = ss4*ta; */\
		"fmul	v5.2d,v1.2d,v14.2d	\n\t	fmul	v10.2d,v6.2d ,v14.2d	\n\t"/* S5  = ss5*ta; */\
		"fmul	v1.2d,v1.2d,v0.2d	\n\t	fmul	v6.2d ,v6.2d ,v0.2d		\n\t"/* S1  = ss1*ta; Moved this one to last MUL slot because unmultiplied v1,6 needed by above MULs */\
		"ldp	q16,q15,[x9]	\n\t"/* load t9 */\
		"fmla	v1.2d,v11.2d,v16.2d	\n\t	fmla	v6.2d ,v11.2d,v15.2d	\n\t"/* S1 += ss2*t9; */\
		"fmla	v2.2d,v13.2d,v16.2d	\n\t	fmla	v7.2d ,v13.2d,v15.2d	\n\t"/* S2 += ss4*t9; */\
		"fmls	v3.2d,v14.2d,v16.2d	\n\t	fmls	v8.2d ,v14.2d,v15.2d	\n\t"/* S3 -= ss5*t9; */\
		"fmls	v4.2d,v12.2d,v16.2d	\n\t	fmls	v9.2d ,v12.2d,v15.2d	\n\t"/* S4 -= ss3*t9; */\
		"fmls	v5.2d,v0.2d ,v16.2d	\n\t	fmls	v10.2d,v0.2d ,v15.2d	\n\t"/* S5 -= ss1*t9; */\
		"ldp	q16,q15,[x8]	\n\t"/* load t8 */\
		"fmla	v1.2d,v12.2d,v16.2d	\n\t	fmla	v6.2d ,v12.2d,v15.2d	\n\t"/* S1 += ss3*t8; */\
		"fmls	v2.2d,v14.2d,v16.2d	\n\t	fmls	v7.2d ,v14.2d,v15.2d	\n\t"/* S2 -= ss5*t8; */\
		"fmls	v3.2d,v11.2d,v16.2d	\n\t	fmls	v8.2d ,v11.2d,v15.2d	\n\t"/* S3 -= ss2*t8; */\
		"fmla	v4.2d,v0.2d ,v16.2d	\n\t	fmla	v9.2d ,v0.2d ,v15.2d	\n\t"/* S4 += ss1*t8; */\
		"fmla	v5.2d,v13.2d,v16.2d	\n\t	fmla	v10.2d,v13.2d,v15.2d	\n\t"/* S5 += ss4*t8; */\
		"ldp	q16,q15,[x7]	\n\t"/* load t7 */\
		"fmla	v1.2d,v13.2d,v16.2d	\n\t	fmla	v6.2d ,v13.2d,v15.2d	\n\t"/* S1 += ss4*t7; */\
		"fmls	v2.2d,v12.2d,v16.2d	\n\t	fmls	v7.2d ,v12.2d,v15.2d	\n\t"/* S2 -= ss3*t7; */\
		"fmla	v3.2d,v0.2d ,v16.2d	\n\t	fmla	v8.2d ,v0.2d ,v15.2d	\n\t"/* S3 += ss1*t7; */\
		"fmla	v4.2d,v14.2d,v16.2d	\n\t	fmla	v9.2d ,v14.2d,v15.2d	\n\t"/* S4 += ss5*t7; */\
		"fmls	v5.2d,v11.2d,v16.2d	\n\t	fmls	v10.2d,v11.2d,v15.2d	\n\t"/* S5 -= ss2*t7; */\
		"ldp	q16,q15,[x6]	\n\t"/* load t6 */\
		"fmla	v1.2d,v14.2d,v16.2d	\n\t	fmla	v6.2d ,v14.2d,v15.2d	\n\t"/* S1 += ss5*t6; */\
		"fmls	v2.2d,v0.2d ,v16.2d	\n\t	fmls	v7.2d ,v0.2d ,v15.2d	\n\t"/* S2 -= ss1*t6; */\
		"fmla	v3.2d,v13.2d,v16.2d	\n\t	fmla	v8.2d ,v13.2d,v15.2d	\n\t"/* S3 += ss4*t6; */\
		"fmls	v4.2d,v11.2d,v16.2d	\n\t	fmls	v9.2d ,v11.2d,v15.2d	\n\t"/* S4 -= ss2*t6; */\
		"fmla	v5.2d,v12.2d,v16.2d	\n\t	fmla	v10.2d,v12.2d,v15.2d	\n\t"/* S5 += ss3*t6; */\
		/* Compute &o[6-A] from 32-bit offsets: */\
		"ldr	x15,%[out0]		\n\t	ldr	x16,%[o_off]	\n\t"\
		"ldr	w10,  [x16,#0x28]	\n\t	add	x10,x15,x10	\n\t"\
		"ldp	w8,w9,[x16,#0x20]	\n\t	add	x8,x15,x8	\n\t	add	x9,x15,x9	\n\t"\
		"ldp	w6,w7,[x16,#0x18]	\n\t	add	x6,x15,x6	\n\t	add	x7,x15,x7	\n\t"\
	/* Write sine-term-subconvo outputs back to memory to free up registers for the 5x5 cosine-term-subconvo computation: */\
		"stp	q1 ,q6 ,[x10]	\n\t"/* s1r,i */\
		"stp	q2 ,q7 ,[x9 ]	\n\t"/* s2r,i */\
		"stp	q3 ,q8 ,[x8 ]	\n\t"/* s3r,i */\
		"stp	q4 ,q9 ,[x7 ]	\n\t"/* s4r,i */\
		"stp	q5 ,q10,[x6 ]	\n\t"/* s5r,i */\
/*
	Cr1 = t0r+cc1*t1r+cc2*t2r+cc3*t3r+cc4*t4r+cc5*t5r;		Ci1 = t0i+cc1*t1i+cc2*t2i+cc3*t3i+cc4*t4i+cc5*t5i;	// C1
	Cr2 = t0r+cc2*t1r+cc4*t2r+cc5*t3r+cc3*t4r+cc1*t5r;		Ci2 = t0i+cc2*t1i+cc4*t2i+cc5*t3i+cc3*t4i+cc1*t5i;	// C2
	Cr3 = t0r+cc3*t1r+cc5*t2r+cc2*t3r+cc1*t4r+cc4*t5r;		Ci3 = t0i+cc3*t1i+cc5*t2i+cc2*t3i+cc1*t4i+cc4*t5i;	// C3
	Cr4 = t0r+cc4*t1r+cc3*t2r+cc1*t3r+cc5*t4r+cc2*t5r;		Ci4 = t0i+cc4*t1i+cc3*t2i+cc1*t3i+cc5*t4i+cc2*t5i;	// C4
	Cr5 = t0r+cc5*t1r+cc1*t2r+cc4*t3r+cc2*t4r+cc3*t5r;		Ci5 = t0i+cc5*t1i+cc1*t2i+cc4*t3i+cc2*t4i+cc3*t5i;	// C5
	B0r = t0r+    t1r+    t2r+    t3r+    t4r+    t5r;		B0i = t0i+    t1i+    t2i+    t3i+    t4i+    t5i;	// X0

	In a 16-reg FMA model, makes sense to init the c-terms to the t0-summand,
	and keep the 10 c-terms and 4 of the 5 shared sincos data in registers. That uses 14 regs,
	leaving 2 for the 2 t*[r,i] data shared by each such block. Those need to be in-reg (at least
	in the cosine-mults section) because of the need to accumulate the DC terms: E.g. at end of the
	first block below (after doing the 10 c-term-update FMAs) we add the current values of B0r,i
	(which we init to t0r,i in the ASM version) to t1r,i, then write the result to memory and
	load t2r,i into the same 2 regs in preparation for the next such block.
*/\
	"sub	x11,x11,#0x50	\n\t"/* Decr trig ptr to point to cc1 */\
		"ld1r	{v0.2d},[x11]		\n\t"/* cc1 */\
		"ldp	q11,q12,[x11,#0x10]	\n\t"/* cc2,cc3 */\
		"ldp	q13,q14,[x11,#0x30]	\n\t"/* cc4,cc5 */\
		"ldp	q1,q6,[x0]	\n\t"	/* c1 = A0 */\
		"ldr	w0,[x16]		\n\t	add	x0,x15,x0	\n\t"/* Overwrite I0 with O0 */\
		"mov	v2.16b,v1.16b	\n\t	mov	v7.16b,v6.16b	\n\t"/* c2 = A0 */\
		"mov	v3.16b,v1.16b	\n\t	mov	v8.16b,v6.16b	\n\t"/* c3 = A0 */\
		"mov	v4.16b,v1.16b	\n\t	mov	v9.16b,v6.16b	\n\t"/* c4 = A0 */\
		"mov	v5.16b,v1.16b	\n\t	mov	v10.16b,v6.16b	\n\t"/* c5 = A0 */\
		"mov	v20.16b,v1.16b	\n\t	mov	v21.16b,v6.16b	\n\t"/* B0 = A0; Init registers used to accumulate DC-component */\
		"ldp	q16,q15,[x1]	\n\t"/* load t1 */\
		"fmla	v1.2d,v0.2d ,v16.2d	\n\t	fmla	v6.2d ,v0.2d ,v15.2d	\n\t"/* C1 += cc1*ta; */\
		"fmla	v2.2d,v11.2d,v16.2d	\n\t	fmla	v7.2d ,v11.2d,v15.2d	\n\t"/* C2 += cc2*ta; */\
		"fmla	v3.2d,v12.2d,v16.2d	\n\t	fmla	v8.2d ,v12.2d,v15.2d	\n\t"/* C3 += cc3*ta; */\
		"fmla	v4.2d,v13.2d,v16.2d	\n\t	fmla	v9.2d ,v13.2d,v15.2d	\n\t"/* C4 += cc4*ta; */\
		"fmla	v5.2d,v14.2d,v16.2d	\n\t	fmla	v10.2d,v14.2d,v15.2d	\n\t"/* C5 += cc5*ta; */\
		"fadd	v20.2d,v20.2d,v16.2d	\n\t	fadd	v21.2d,v21.2d,v15.2d	\n\t"/* B0 += t1 */\
		"ldp	q16,q15,[x2]	\n\t"/* load t2 */\
		"fmla	v1.2d,v11.2d,v16.2d	\n\t	fmla	v6.2d ,v11.2d,v15.2d	\n\t"/* C1 += cc2*t9; */\
		"fmla	v2.2d,v13.2d,v16.2d	\n\t	fmla	v7.2d ,v13.2d,v15.2d	\n\t"/* C2 += cc4*t9; */\
		"fmla	v3.2d,v14.2d,v16.2d	\n\t	fmla	v8.2d ,v14.2d,v15.2d	\n\t"/* C3 += cc5*t9; */\
		"fmla	v4.2d,v12.2d,v16.2d	\n\t	fmla	v9.2d ,v12.2d,v15.2d	\n\t"/* C4 += cc3*t9; */\
		"fmla	v5.2d,v0.2d ,v16.2d	\n\t	fmla	v10.2d,v0.2d ,v15.2d	\n\t"/* C5 += cc1*t9; */\
		"fadd	v20.2d,v20.2d,v16.2d	\n\t	fadd	v21.2d,v21.2d,v15.2d	\n\t"/* B0 += t2 */\
		"ldp	q16,q15,[x3]	\n\t"/* load t3 */\
		"fmla	v1.2d,v12.2d,v16.2d	\n\t	fmla	v6.2d ,v12.2d,v15.2d	\n\t"/* C1 += cc3*t8; */\
		"fmla	v2.2d,v14.2d,v16.2d	\n\t	fmla	v7.2d ,v14.2d,v15.2d	\n\t"/* C2 += cc5*t8; */\
		"fmla	v3.2d,v11.2d,v16.2d	\n\t	fmla	v8.2d ,v11.2d,v15.2d	\n\t"/* C3 += cc2*t8; */\
		"fmla	v4.2d,v0.2d ,v16.2d	\n\t	fmla	v9.2d ,v0.2d ,v15.2d	\n\t"/* C4 += cc1*t8; */\
		"fmla	v5.2d,v13.2d,v16.2d	\n\t	fmla	v10.2d,v13.2d,v15.2d	\n\t"/* C5 += cc4*t8; */\
		"fadd	v20.2d,v20.2d,v16.2d	\n\t	fadd	v21.2d,v21.2d,v15.2d	\n\t"/* B0 += t3 */\
		"ldp	q16,q15,[x4]	\n\t"/* load t4 */\
		"fmla	v1.2d,v13.2d,v16.2d	\n\t	fmla	v6.2d ,v13.2d,v15.2d	\n\t"/* C1 += cc4*t7; */\
		"fmla	v2.2d,v12.2d,v16.2d	\n\t	fmla	v7.2d ,v12.2d,v15.2d	\n\t"/* C2 += cc3*t7; */\
		"fmla	v3.2d,v0.2d ,v16.2d	\n\t	fmla	v8.2d ,v0.2d ,v15.2d	\n\t"/* C3 += cc1*t7; */\
		"fmla	v4.2d,v14.2d,v16.2d	\n\t	fmla	v9.2d ,v14.2d,v15.2d	\n\t"/* C4 += cc5*t7; */\
		"fmla	v5.2d,v11.2d,v16.2d	\n\t	fmla	v10.2d,v11.2d,v15.2d	\n\t"/* C5 += cc2*t7; */\
		"fadd	v20.2d,v20.2d,v16.2d	\n\t	fadd	v21.2d,v21.2d,v15.2d	\n\t"/* B0 += t4 */\
		"ldp	q16,q15,[x5]	\n\t"/* load t5 */\
		"fmla	v1.2d,v14.2d,v16.2d	\n\t	fmla	v6.2d ,v14.2d,v15.2d	\n\t"/* C1 += cc5*t6; */\
		"fmla	v2.2d,v0.2d ,v16.2d	\n\t	fmla	v7.2d ,v0.2d ,v15.2d	\n\t"/* C2 += cc1*t6; */\
		"fmla	v3.2d,v13.2d,v16.2d	\n\t	fmla	v8.2d ,v13.2d,v15.2d	\n\t"/* C3 += cc4*t6; */\
		"fmla	v4.2d,v11.2d,v16.2d	\n\t	fmla	v9.2d ,v11.2d,v15.2d	\n\t"/* C4 += cc2*t6; */\
		"fmla	v5.2d,v12.2d,v16.2d	\n\t	fmla	v10.2d,v12.2d,v15.2d	\n\t"/* C5 += cc3*t6; */\
		"fadd	v20.2d,v20.2d,v16.2d	\n\t	fadd	v21.2d,v21.2d,v15.2d	\n\t"/* B0 += t5 */\
		"stp	q20,q21,[x0]	\n\t"/* Store B0 */\
		/* Compute &o[1-5] from 32-bit offsets: */\
		"ldr	w1,   [x16,#0x04]	\n\t	add	x1,x15,x1	\n\t"\
		"ldp	w2,w3,[x16,#0x08]	\n\t	add	x2,x15,x2	\n\t	add	x3,x15,x3	\n\t"\
		"ldp	w4,w5,[x16,#0x10]	\n\t	add	x4,x15,x4	\n\t	add	x5,x15,x5	\n\t"\
		/* Reload S-terms: */\
		"ldp	q11,q16,[x10]	\n\t"/* s1r,i */\
		"ldp	q12,q17,[x9 ]	\n\t"/* s2r,i */\
		"ldp	q13,q18,[x8 ]	\n\t"/* s3r,i */\
		"ldp	q14,q19,[x7 ]	\n\t"/* s4r,i */\
		"ldp	q15,q20,[x6 ]	\n\t"/* s5r,i */\
	/* ...Do the lower-half-index portion of the closing wave of radix-2 butterflies: */\
		"fsub	v1.2d,v1.2d,v16.2d	\n\t	fadd	v6.2d ,v6.2d ,v11.2d	\n\t"/* B1 = C1 + I*S1 */\
		"fsub	v2.2d,v2.2d,v17.2d	\n\t	fadd	v7.2d ,v7.2d ,v12.2d	\n\t"/* B2 = C2 + I*S2 */\
		"fsub	v3.2d,v3.2d,v18.2d	\n\t	fadd	v8.2d ,v8.2d ,v13.2d	\n\t"/* B3 = C3 + I*S3 */\
		"fsub	v4.2d,v4.2d,v19.2d	\n\t	fadd	v9.2d ,v9.2d ,v14.2d	\n\t"/* B4 = C4 + I*S4 */\
		"fsub	v5.2d,v5.2d,v20.2d	\n\t	fadd	v10.2d,v10.2d,v15.2d	\n\t"/* B5 = C5 + I*S5 */\
	/* Write resulting lower-half outputs B1-5 back to memory before ensuing FMA-block to compute - parts overwrites each + result: */\
		"stp	q1 ,q6 ,[x1]	\n\t	fmla	v1.2d,v31.2d,v16.2d	\n\t	fmls	v6.2d ,v31.2d,v11.2d	\n\t"/* B10= C1 - I*S1 */\
		"stp	q2 ,q7 ,[x2]	\n\t	fmla	v2.2d,v31.2d,v17.2d	\n\t	fmls	v7.2d ,v31.2d,v12.2d	\n\t"/* B9 = C2 - I*S2 */\
		"stp	q3 ,q8 ,[x3]	\n\t	fmla	v3.2d,v31.2d,v18.2d	\n\t	fmls	v8.2d ,v31.2d,v13.2d	\n\t"/* B8 = C3 - I*S3 */\
		"stp	q4 ,q9 ,[x4]	\n\t	fmla	v4.2d,v31.2d,v19.2d	\n\t	fmls	v9.2d ,v31.2d,v14.2d	\n\t"/* B7 = C4 - I*S4 */\
		"stp	q5 ,q10,[x5]	\n\t	fmla	v5.2d,v31.2d,v20.2d	\n\t	fmls	v10.2d,v31.2d,v15.2d	\n\t"/* B6 = C5 - I*S5 */\
	/* ...And write the upper-half outputs back to memory. */\
		"stp	q1 ,q6 ,[x10]	\n\t"/* Bar,i */\
		"stp	q2 ,q7 ,[x9 ]	\n\t"/* B9r,i */\
		"stp	q3 ,q8 ,[x8 ]	\n\t"/* B8r,i */\
		"stp	q4 ,q9 ,[x7 ]	\n\t"/* B7r,i */\
		"stp	q5 ,q10,[x6 ]	\n\t"/* B6r,i */\
		:					/* outputs: none */\
		: [in0] "m" (Xin0)	/* All inputs from memory addresses here */\
		 ,[i_off] "m" (Xi_off)	/* I-address pointer-strides */\
		 ,[__cc] "m" (Xcc)\
		 ,[out0] "m" (Xout0)\
		 ,[o_off] "m" (Xo_off)	/* O-address pointer-strides */\
		: "cc","memory","x0","x1","x2","x3","x4","x5","x6","x7","x8","x9","x10","x11","x15","x16",\
		"v0","v1","v2","v3","v4","v5","v6","v7","v8","v9","v10","v11","v12","v13","v14","v15","v16","v17","v18","v19","v20","v21","v31"	/* Clobbered registers */\
		);\
	}

#elif defined(USE_AVX512)	// AVX512 version, based on RADIX_11_DFT_FMA in dft_macro.h

   #if DFT_11_FMA	// [1] Naive but good-ROE and FMA-friendly impl based on RADIX_11_DFT_BASIC in dft_macro.h:

	// FMAs used for all arithmetic, including 'trivial' ones (one mult = 1.0) to replace ADD/SUB:
	//
	// Arithmetic opcount: [10 ADD, 140 FMA (20 trivial, incl 10 MUL), 166 memref], close to general target of 1 memref per vec_dbl arithmetic op.
	// potentially much lower cycle count (at a max theoretical rate of 2 FMA/cycle) than non-FMA version.
	//
	// Compare to non-FMA: [182 ADD, 44 MUL), 154 memref], less accurate and bottlenecked by 1 ADD/cycle max issue rate.
	// And to DFT_11_FMA = 0 version below: [134 ADD, 54 MUL/FMA, 172 memref].
	//
	#define SSE2_RADIX_11_DFT(Xin0,Xi_off, Xcc, Xout0,Xo_off)\
	{\
	__asm__ volatile (\
	"movq	%[in0],%%rax		\n\t	movq	%[out0],%%rcx	\n\t"/* rax,rcx hold to I/O base-pointers */\
	"movq	%[i_off],%%rbx		\n\t	movq	%[o_off],%%rdx	\n\t"/* rbx,rdx hold pointers to uint32[] I/O pointer-offset-stride arrays */\
		/********************************************/\
		/*       Here are the 5 cosine terms:       */\
		/********************************************/\
		/* Compute &i[1-A] from 32-bit offsets: */\
		"movslq	0x04(%%rbx),%%rsi	\n\t	addq	%%rax,%%rsi		\n\t"\
		"movslq	0x08(%%rbx),%%rdi	\n\t	addq	%%rax,%%rdi		\n\t"\
		"movslq	0x0c(%%rbx),%%r8 	\n\t	addq	%%rax,%%r8 		\n\t"\
		"movslq	0x10(%%rbx),%%r9 	\n\t	addq	%%rax,%%r9 		\n\t"\
		"movslq	0x14(%%rbx),%%r10	\n\t	addq	%%rax,%%r10		\n\t"\
		"movslq	0x18(%%rbx),%%r11	\n\t	addq	%%rax,%%r11		\n\t"\
		"movslq	0x1c(%%rbx),%%r12	\n\t	addq	%%rax,%%r12		\n\t"\
		"movslq	0x20(%%rbx),%%r13	\n\t	addq	%%rax,%%r13		\n\t"\
		"movslq	0x24(%%rbx),%%r14	\n\t	addq	%%rax,%%r14		\n\t"\
		"movslq	0x28(%%rbx),%%r15	\n\t	addq	%%rax,%%r15		\n\t"\
		"vmovaps	(%%rsi),%%zmm1		\n\t"/* A1r */"		vmovaps	0x40(%%rsi),%%zmm6	\n\t"/* A1i */\
		"vmovaps	(%%rdi),%%zmm2		\n\t"/* A2r */"		vmovaps	0x40(%%rdi),%%zmm7	\n\t"/* A2i */\
		"vmovaps	(%%r8 ),%%zmm3		\n\t"/* A3r */"		vmovaps	0x40(%%r8 ),%%zmm8	\n\t"/* A3i */\
		"vmovaps	(%%r9 ),%%zmm4		\n\t"/* A4r */"		vmovaps	0x40(%%r9 ),%%zmm9	\n\t"/* A4i */\
		"vmovaps	(%%r10),%%zmm5		\n\t"/* A5r */"		vmovaps	0x40(%%r10),%%zmm10	\n\t"/* A5i */\
		/* Have enough registers to also hold the real-upper-half inputs: */\
		"vmovaps	(%%r11),%%zmm11		\n\t"/* A6r */\
		"vmovaps	(%%r12),%%zmm12		\n\t"/* A7r */\
		"vmovaps	(%%r13),%%zmm13		\n\t"/* A8r */\
		"vmovaps	(%%r14),%%zmm14		\n\t"/* A9r */\
		"vmovaps	(%%r15),%%zmm15		\n\t"/* AAr */\
	/* ...Do the + parts of the opening wave of radix-2 butterflies: */\
		"vaddpd %%zmm15,%%zmm1,%%zmm1	\n\t	vaddpd 0x40(%%r15),%%zmm6 ,%%zmm6 	\n\t"/* t1 = A1 + Aa */\
		"vaddpd %%zmm14,%%zmm2,%%zmm2	\n\t	vaddpd 0x40(%%r14),%%zmm7 ,%%zmm7 	\n\t"/* t2 = A2 + A9 */\
		"vaddpd %%zmm13,%%zmm3,%%zmm3	\n\t	vaddpd 0x40(%%r13),%%zmm8 ,%%zmm8 	\n\t"/* t3 = A3 + A8 */\
		"vaddpd %%zmm12,%%zmm4,%%zmm4	\n\t	vaddpd 0x40(%%r12),%%zmm9 ,%%zmm9 	\n\t"/* t4 = A4 + A7 */\
		"vaddpd %%zmm11,%%zmm5,%%zmm5	\n\t	vaddpd 0x40(%%r11),%%zmm10,%%zmm10	\n\t"/* t5 = A5 + A6 */\
	/* Write lower-half outputs back to memory... */\
		"vmovaps	%%zmm1,(%%rsi)		\n\t"/* t1r */"		vmovaps	%%zmm6 ,0x40(%%rsi)	\n\t"/* t1i */\
		"vmovaps	%%zmm2,(%%rdi)		\n\t"/* t2r */"		vmovaps	%%zmm7 ,0x40(%%rdi)	\n\t"/* t2i */\
		"vmovaps	%%zmm3,(%%r8 )		\n\t"/* t3r */"		vmovaps	%%zmm8 ,0x40(%%r8 )	\n\t"/* t3i */\
		"movq $0x4000000000000000,%%r8	\n\t	vpbroadcastq  %%r8,%%zmm0	\n\t"/* 2.0 */\
		"vmovaps	%%zmm4,(%%r9 )		\n\t"/* t4r */"		vmovaps	%%zmm9 ,0x40(%%r9 )	\n\t"/* t4i */\
		"vmovaps	%%zmm5,(%%r10)		\n\t"/* t5r */"		vmovaps	%%zmm10,0x40(%%r10)	\n\t"/* t5i */\
	/* ...Do the - parts of the radix-2 butterflies: */\
	"vfnmadd231pd %%zmm15,%%zmm0,%%zmm1		\n\t	vfnmadd231pd 0x40(%%r15),%%zmm0,%%zmm6	\n\t"/* ta = A1 - Aa */\
	"vfnmadd231pd %%zmm14,%%zmm0,%%zmm2		\n\t	vfnmadd231pd 0x40(%%r14),%%zmm0,%%zmm7	\n\t"/* t9 = A2 - A9 */\
	"vfnmadd231pd %%zmm13,%%zmm0,%%zmm3		\n\t	vfnmadd231pd 0x40(%%r13),%%zmm0,%%zmm8	\n\t"/* t8 = A3 - A8 */\
	"vfnmadd231pd %%zmm12,%%zmm0,%%zmm4		\n\t	vfnmadd231pd 0x40(%%r12),%%zmm0,%%zmm9	\n\t"/* t7 = A4 - A7 */\
	"vfnmadd231pd %%zmm11,%%zmm0,%%zmm5		\n\t	vfnmadd231pd 0x40(%%r11),%%zmm0,%%zmm10	\n\t"/* t6 = A5 - A6 */\
	/* ...And write upper-half outputs 6-9 (A stays in-reg) back to memory to free up registers for the 5x5 sine-term-subconvo computation: */\
		"vmovaps	%%zmm2,(%%r14)		\n\t"/* t9r */"		vmovaps	%%zmm7 ,0x40(%%r14)	\n\t"/* t9i */\
		"vmovaps	%%zmm3,(%%r13)		\n\t"/* t8r */"		vmovaps	%%zmm8 ,0x40(%%r13)	\n\t"/* t8i */\
		"vmovaps	%%zmm4,(%%r12)		\n\t"/* t7r */"		vmovaps	%%zmm9 ,0x40(%%r12)	\n\t"/* t7i */\
		"vmovaps	%%zmm5,(%%r11)		\n\t"/* t6r */"		vmovaps	%%zmm10,0x40(%%r11)	\n\t"/* t6i */\
	/*
		Sr1 =      ss1*tar+ss2*t9r+ss3*t8r+ss4*t7r+ss5*t6r;		Si1 =      ss1*tai+ss2*t9i+ss3*t8i+ss4*t7i+ss5*t6i;	// S1
		Sr2 =      ss2*tar+ss4*t9r-ss5*t8r-ss3*t7r-ss1*t6r;		Si2 =      ss2*tai+ss4*t9i-ss5*t8i-ss3*t7i-ss1*t6i;	// S2
		Sr3 =      ss3*tar-ss5*t9r-ss2*t8r+ss1*t7r+ss4*t6r;		Si3 =      ss3*tai-ss5*t9i-ss2*t8i+ss1*t7i+ss4*t6i;	// S3
		Sr4 =      ss4*tar-ss3*t9r+ss1*t8r+ss5*t7r-ss2*t6r;		Si4 =      ss4*tai-ss3*t9i+ss1*t8i+ss5*t7i-ss2*t6i;	// S4
		Sr5 =      ss5*tar-ss1*t9r+ss4*t8r-ss2*t7r+ss3*t6r;		Si5 =      ss5*tai-ss1*t9i+ss4*t8i-ss2*t7i+ss3*t6i;	// S5

		In a 16-reg FMA model, can keep the 10 S-terms and all 5 shared sincos data in registers, but that
		requires us to reload one of the 2 repeated t*-mults (e.g. tai in the 1st block) 5 times per block.
		OTOH if we keep both t*r,i-mults and 4 of the 5 ss-mults in-reg, we just need 2 loads-from-mem per block:
	*/\
	"movq	%[__cc],%%r9	\n\t"/* cc1 */\
		"vmovaps	%%zmm1,%%zmm2		\n\t	vmovaps	%%zmm6 ,%%zmm7		\n\t"/* S2 = ta */\
		"vmovaps	%%zmm1,%%zmm3		\n\t	vmovaps	%%zmm6 ,%%zmm8		\n\t"/* S3 = ta */\
		"vmovaps	%%zmm1,%%zmm4		\n\t	vmovaps	%%zmm6 ,%%zmm9		\n\t"/* S4 = ta */\
		"vmovaps	%%zmm1,%%zmm5		\n\t	vmovaps	%%zmm6 ,%%zmm10		\n\t"/* S5 = ta */\
	"addq	$0x140,%%r9	\n\t"/* Incr trig ptr to point to ss1 */\
		"vmovaps	0x040(%%r9),%%zmm11	\n\t	vmovaps	0x080(%%r9),%%zmm12		\n\t"/* ss2,ss3 */\
		"vmovaps	0x0c0(%%r9),%%zmm13	\n\t	vmovaps	0x100(%%r9),%%zmm14		\n\t"/* ss4,ss5 */\
\
	"	vmulpd	(%%r9),%%zmm1,%%zmm1	\n\t	vmulpd	(%%r9),%%zmm6 ,%%zmm6	\n\t"/* S1  = ss1*ta; */\
	"	vmulpd	%%zmm11,%%zmm2,%%zmm2	\n\t	vmulpd	%%zmm11,%%zmm7 ,%%zmm7	\n\t"/* S2  = ss2*ta; */\
	"	vmulpd	%%zmm12,%%zmm3,%%zmm3	\n\t	vmulpd	%%zmm12,%%zmm8 ,%%zmm8	\n\t"/* S3  = ss3*ta; */\
	"	vmulpd	%%zmm13,%%zmm4,%%zmm4	\n\t	vmulpd	%%zmm13,%%zmm9 ,%%zmm9	\n\t"/* S4  = ss4*ta; */\
	"	vmulpd	%%zmm14,%%zmm5,%%zmm5	\n\t	vmulpd	%%zmm14,%%zmm10,%%zmm10	\n\t"/* S5  = ss5*ta; */\
\
		"vmovaps	    (%%r14),%%zmm0	\n\t	vmovaps	0x40(%%r14),%%zmm15	\n\t"/* t9 */\
	" vfmadd231pd %%zmm11,%%zmm0,%%zmm1	\n\t  vfmadd231pd %%zmm11,%%zmm15,%%zmm6	\n\t"/* S1 += ss2*t9; */\
	" vfmadd231pd %%zmm13,%%zmm0,%%zmm2	\n\t  vfmadd231pd %%zmm13,%%zmm15,%%zmm7	\n\t"/* S2 += ss4*t9; */\
	"vfnmadd231pd %%zmm14,%%zmm0,%%zmm3	\n\t vfnmadd231pd %%zmm14,%%zmm15,%%zmm8	\n\t"/* S3 -= ss5*t9; */\
	"vfnmadd231pd %%zmm12,%%zmm0,%%zmm4	\n\t vfnmadd231pd %%zmm12,%%zmm15,%%zmm9	\n\t"/* S4 -= ss3*t9; */\
	"vfnmadd231pd (%%r9 ),%%zmm0,%%zmm5	\n\t vfnmadd231pd (%%r9 ),%%zmm15,%%zmm10	\n\t"/* S5 -= ss1*t9; */\
\
		"vmovaps	    (%%r13),%%zmm0	\n\t	vmovaps	0x40(%%r13),%%zmm15	\n\t"/* t8 */\
	" vfmadd231pd %%zmm12,%%zmm0,%%zmm1	\n\t  vfmadd231pd %%zmm12,%%zmm15,%%zmm6	\n\t"/* S1 += ss3*t8; */\
	"vfnmadd231pd %%zmm14,%%zmm0,%%zmm2	\n\t vfnmadd231pd %%zmm14,%%zmm15,%%zmm7	\n\t"/* S2 -= ss5*t8; */\
	"vfnmadd231pd %%zmm11,%%zmm0,%%zmm3	\n\t vfnmadd231pd %%zmm11,%%zmm15,%%zmm8	\n\t"/* S3 -= ss2*t8; */\
	" vfmadd231pd (%%r9 ),%%zmm0,%%zmm4	\n\t  vfmadd231pd (%%r9 ),%%zmm15,%%zmm9	\n\t"/* S4 += ss1*t8; */\
	" vfmadd231pd %%zmm13,%%zmm0,%%zmm5	\n\t  vfmadd231pd %%zmm13,%%zmm15,%%zmm10	\n\t"/* S5 += ss4*t8; */\
\
		"vmovaps	    (%%r12),%%zmm0	\n\t	vmovaps	0x40(%%r12),%%zmm15	\n\t"/* t7 */\
	" vfmadd231pd %%zmm13,%%zmm0,%%zmm1	\n\t  vfmadd231pd %%zmm13,%%zmm15,%%zmm6	\n\t"/* S1 += ss4*t7; */\
	"vfnmadd231pd %%zmm12,%%zmm0,%%zmm2	\n\t vfnmadd231pd %%zmm12,%%zmm15,%%zmm7	\n\t"/* S2 -= ss3*t7; */\
	" vfmadd231pd (%%r9 ),%%zmm0,%%zmm3	\n\t  vfmadd231pd (%%r9 ),%%zmm15,%%zmm8	\n\t"/* S3 += ss1*t7; */\
	" vfmadd231pd %%zmm14,%%zmm0,%%zmm4	\n\t  vfmadd231pd %%zmm14,%%zmm15,%%zmm9	\n\t"/* S4 += ss5*t7; */\
	"vfnmadd231pd %%zmm11,%%zmm0,%%zmm5	\n\t vfnmadd231pd %%zmm11,%%zmm15,%%zmm10	\n\t"/* S5 -= ss2*t7; */\
\
		"vmovaps	    (%%r11),%%zmm0	\n\t	vmovaps	0x40(%%r11),%%zmm15	\n\t"/* t6 */\
	" vfmadd231pd %%zmm14,%%zmm0,%%zmm1	\n\t  vfmadd231pd %%zmm14,%%zmm15,%%zmm6	\n\t"/* S1 += ss5*t6; */\
	"vfnmadd231pd (%%r9 ),%%zmm0,%%zmm2	\n\t vfnmadd231pd (%%r9 ),%%zmm15,%%zmm7	\n\t"/* S2 -= ss1*t6; */\
	" vfmadd231pd %%zmm13,%%zmm0,%%zmm3	\n\t  vfmadd231pd %%zmm13,%%zmm15,%%zmm8	\n\t"/* S3 += ss4*t6; */\
	"vfnmadd231pd %%zmm11,%%zmm0,%%zmm4	\n\t vfnmadd231pd %%zmm11,%%zmm15,%%zmm9	\n\t"/* S4 -= ss2*t6; */\
	" vfmadd231pd %%zmm12,%%zmm0,%%zmm5	\n\t  vfmadd231pd %%zmm12,%%zmm15,%%zmm10	\n\t"/* S5 += ss3*t6; */\
\
	/* Write sine-term-subconvo outputs back to memory to free up registers for the 5x5 cosine-term-subconvo computation: */\
	/* These i-ptrs no longer needed, so load o-ptrs in their place: */\
		"movslq	0x28(%%rdx),%%r15	\n\t	addq	%%rcx,%%r15		\n\t"/* &BA */\
		"movslq	0x24(%%rdx),%%r14	\n\t	addq	%%rcx,%%r14		\n\t"/* &B9 */\
		"movslq	0x20(%%rdx),%%r13	\n\t	addq	%%rcx,%%r13		\n\t"/* &B8 */\
		"movslq	0x1c(%%rdx),%%r12	\n\t	addq	%%rcx,%%r12		\n\t"/* &B7 */\
		"movslq	0x18(%%rdx),%%r11	\n\t	addq	%%rcx,%%r11		\n\t"/* &B6 */\
		"vmovaps	%%zmm1,(%%r15)		\n\t"/* s1r */"		vmovaps	%%zmm6 ,0x40(%%r15)	\n\t"/* s1i */\
		"vmovaps	%%zmm2,(%%r14)		\n\t"/* s2r */"		vmovaps	%%zmm7 ,0x40(%%r14)	\n\t"/* s2i */\
		"vmovaps	%%zmm3,(%%r13)		\n\t"/* s3r */"		vmovaps	%%zmm8 ,0x40(%%r13)	\n\t"/* s3i */\
		"vmovaps	%%zmm4,(%%r12)		\n\t"/* s4r */"		vmovaps	%%zmm9 ,0x40(%%r12)	\n\t"/* s4i */\
		"vmovaps	%%zmm5,(%%r11)		\n\t"/* s5r */"		vmovaps	%%zmm10,0x40(%%r11)	\n\t"/* s5i */\
	/*
		Cr1 = t0r+cc1*t1r+cc2*t2r+cc3*t3r+cc4*t4r+cc5*t5r;		Ci1 = t0i+cc1*t1i+cc2*t2i+cc3*t3i+cc4*t4i+cc5*t5i;	// C1
		Cr2 = t0r+cc2*t1r+cc4*t2r+cc5*t3r+cc3*t4r+cc1*t5r;		Ci2 = t0i+cc2*t1i+cc4*t2i+cc5*t3i+cc3*t4i+cc1*t5i;	// C2
		Cr3 = t0r+cc3*t1r+cc5*t2r+cc2*t3r+cc1*t4r+cc4*t5r;		Ci3 = t0i+cc3*t1i+cc5*t2i+cc2*t3i+cc1*t4i+cc4*t5i;	// C3
		Cr4 = t0r+cc4*t1r+cc3*t2r+cc1*t3r+cc5*t4r+cc2*t5r;		Ci4 = t0i+cc4*t1i+cc3*t2i+cc1*t3i+cc5*t4i+cc2*t5i;	// C4
		Cr5 = t0r+cc5*t1r+cc1*t2r+cc4*t3r+cc2*t4r+cc3*t5r;		Ci5 = t0i+cc5*t1i+cc1*t2i+cc4*t3i+cc2*t4i+cc3*t5i;	// C5
		B0r = t0r+    t1r+    t2r+    t3r+    t4r+    t5r;		B0i = t0i+    t1i+    t2i+    t3i+    t4i+    t5i;	// X0

		In a 16-reg FMA model, makes sense to init the c-terms to the t0-summand,
		and keep the 10 c-terms and 4 of the 5 shared sincos data in registers. That uses 14 regs,
		leaving 2 for the 2 t*[r,i] data shared by each such block. Those need to be in-reg (at least
		in the cosine-mults section) because of the need to accumulate the DC terms: E.g. at end of the
		first block below (after doing the 10 c-term-update FMAs) we add the current values of B0r,i
		(which we init to t0r,i in the ASM version) to t1r,i, then write the result to memory and
		load t2r,i into the same 2 regs in preparation for the next such block.
	*/\
	/*
	Current GPR assignments:
		rax	rbx	rcx	rdx	rsi	rdi	r8	r9	r10	r11	r12	r13	r14	r15
		in	iff	out	off	i1	i2	2.0	ss1	i5	o6	o7	o8	o9	oA
	*/\
	"subq	$0x140,%%r9	\n\t"/* Decr trig ptr to point to cc1 */\
	"movslq	(%%rbx),%%r8	\n\t	addq	%%rax,%%r8	\n\t"/* i0-from-offset */\
		"vmovaps	(%%r8 ),%%zmm1		\n\t	vmovaps	0x40(%%r8 ),%%zmm6	\n\t"/* c1 = A0 */\
		"vmovaps	%%zmm1,%%zmm2		\n\t	vmovaps	%%zmm6 ,%%zmm7		\n\t"/* c2 = A0 */\
		"vmovaps	%%zmm1,%%zmm3		\n\t	vmovaps	%%zmm6 ,%%zmm8		\n\t"/* c3 = A0 */\
		"vmovaps	%%zmm1,%%zmm4		\n\t	vmovaps	%%zmm6 ,%%zmm9		\n\t"/* c4 = A0 */\
		"vmovaps	%%zmm1,%%zmm5		\n\t	vmovaps	%%zmm6 ,%%zmm10		\n\t"/* c5 = A0 */\
\
		"vmovaps	0x040(%%r9),%%zmm11	\n\t	vmovaps	0x080(%%r9),%%zmm12		\n\t"/* cc2,cc3 */\
		"vmovaps	0x0c0(%%r9),%%zmm13	\n\t	vmovaps	0x100(%%r9),%%zmm14		\n\t"/* cc4,cc5 */\
\
		"vmovaps	    (%%rsi),%%zmm0	\n\t	vmovaps	0x40(%%rsi),%%zmm15	\n\t"/* t1 */\
	"vfmadd231pd (%%r9 ),%%zmm0,%%zmm1	\n\t vfmadd231pd (%%r9 ),%%zmm15,%%zmm6	\n\t"/* c1 += cc1*t1; */\
	"vfmadd231pd %%zmm11,%%zmm0,%%zmm2	\n\t vfmadd231pd %%zmm11,%%zmm15,%%zmm7	\n\t"/* c2 += cc2*t1; */\
	"vfmadd231pd %%zmm12,%%zmm0,%%zmm3	\n\t vfmadd231pd %%zmm12,%%zmm15,%%zmm8	\n\t"/* c3 += cc3*t1; */\
	"vfmadd231pd %%zmm13,%%zmm0,%%zmm4	\n\t vfmadd231pd %%zmm13,%%zmm15,%%zmm9	\n\t"/* c4 += cc4*t1; */\
	"vfmadd231pd %%zmm14,%%zmm0,%%zmm5	\n\t vfmadd231pd %%zmm14,%%zmm15,%%zmm10	\n\t"/* c5 += cc5*t1; */\
		"vaddpd (%%r8 ),%%zmm0,%%zmm0	\n\t	vaddpd 0x40(%%r8 ),%%zmm15,%%zmm15	\n\t"/* B0 += t1; */\
		"vmovaps	%%zmm0,(%%r8 )		\n\t	vmovaps	%%zmm15,0x40(%%r8 )			\n\t"/* Store B0 */\
\
		"vmovaps	    (%%rdi),%%zmm0	\n\t	vmovaps	0x40(%%rdi),%%zmm15	\n\t"/* t2 */\
		"movslq	0x0c(%%rbx),%%rsi	\n\t	addq	%%rax,%%rsi		\n\t"/* Compute &i3,4 from 32-bit offsets */\
		"movslq	0x10(%%rbx),%%rdi	\n\t	addq	%%rax,%%rdi		\n\t"\
	"vfmadd231pd %%zmm11,%%zmm0,%%zmm1	\n\t vfmadd231pd %%zmm11,%%zmm15,%%zmm6	\n\t"/* c1 += cc2*t2; */\
	"vfmadd231pd %%zmm13,%%zmm0,%%zmm2	\n\t vfmadd231pd %%zmm13,%%zmm15,%%zmm7	\n\t"/* c2 += cc4*t2; */\
	"vfmadd231pd %%zmm14,%%zmm0,%%zmm3	\n\t vfmadd231pd %%zmm14,%%zmm15,%%zmm8	\n\t"/* c3 += cc5*t2; */\
	"vfmadd231pd %%zmm12,%%zmm0,%%zmm4	\n\t vfmadd231pd %%zmm12,%%zmm15,%%zmm9	\n\t"/* c4 += cc3*t2; */\
	"vfmadd231pd (%%r9 ),%%zmm0,%%zmm5	\n\t vfmadd231pd (%%r9 ),%%zmm15,%%zmm10	\n\t"/* c5 += cc1*t2; */\
		"vaddpd (%%r8 ),%%zmm0,%%zmm0	\n\t	vaddpd 0x40(%%r8 ),%%zmm15,%%zmm15	\n\t"/* B0 += t2; */\
		"vmovaps	%%zmm0,(%%r8 )		\n\t	vmovaps	%%zmm15,0x40(%%r8 )			\n\t"/* Store B0 */\
\
		"vmovaps	    (%%rsi),%%zmm0	\n\t	vmovaps	0x40(%%rsi),%%zmm15	\n\t"/* t3 */\
	"vfmadd231pd %%zmm12,%%zmm0,%%zmm1	\n\t vfmadd231pd %%zmm12,%%zmm15,%%zmm6	\n\t"/* c1 += cc3*t3; */\
	"vfmadd231pd %%zmm14,%%zmm0,%%zmm2	\n\t vfmadd231pd %%zmm14,%%zmm15,%%zmm7	\n\t"/* c2 += cc5*t3; */\
	"vfmadd231pd %%zmm11,%%zmm0,%%zmm3	\n\t vfmadd231pd %%zmm11,%%zmm15,%%zmm8	\n\t"/* c3 += cc2*t3; */\
	"vfmadd231pd (%%r9 ),%%zmm0,%%zmm4	\n\t vfmadd231pd (%%r9 ),%%zmm15,%%zmm9	\n\t"/* c4 += cc1*t3; */\
	"vfmadd231pd %%zmm13,%%zmm0,%%zmm5	\n\t vfmadd231pd %%zmm13,%%zmm15,%%zmm10	\n\t"/* c5 += cc4*t3; */\
		"vaddpd (%%r8 ),%%zmm0,%%zmm0	\n\t	vaddpd 0x40(%%r8 ),%%zmm15,%%zmm15	\n\t"/* B0 += t3; */\
		"vmovaps	%%zmm0,(%%r8 )		\n\t	vmovaps	%%zmm15,0x40(%%r8 )			\n\t"/* Store B0 */\
\
		"vmovaps	    (%%rdi),%%zmm0	\n\t	vmovaps	0x40(%%rdi),%%zmm15	\n\t"/* t4 */\
	"vfmadd231pd %%zmm13,%%zmm0,%%zmm1	\n\t vfmadd231pd %%zmm13,%%zmm15,%%zmm6	\n\t"/* c1 += cc4*t4; */\
	"vfmadd231pd %%zmm12,%%zmm0,%%zmm2	\n\t vfmadd231pd %%zmm12,%%zmm15,%%zmm7	\n\t"/* c2 += cc3*t4; */\
	"vfmadd231pd (%%r9 ),%%zmm0,%%zmm3	\n\t vfmadd231pd (%%r9 ),%%zmm15,%%zmm8	\n\t"/* c3 += cc1*t4; */\
	"vfmadd231pd %%zmm14,%%zmm0,%%zmm4	\n\t vfmadd231pd %%zmm14,%%zmm15,%%zmm9	\n\t"/* c4 += cc5*t4; */\
	"vfmadd231pd %%zmm11,%%zmm0,%%zmm5	\n\t vfmadd231pd %%zmm11,%%zmm15,%%zmm10	\n\t"/* c5 += cc2*t4; */\
		"vaddpd (%%r8 ),%%zmm0,%%zmm0	\n\t	vaddpd 0x40(%%r8 ),%%zmm15,%%zmm15	\n\t"/* B0 += t4; */\
		"vmovaps	%%zmm0,(%%r8 )		\n\t	vmovaps	%%zmm15,0x40(%%r8 )			\n\t"/* Store B0 */\
\
		"vmovaps	    (%%r10),%%zmm0	\n\t	vmovaps	0x40(%%r10),%%zmm15	\n\t"/* t5 */\
	"vfmadd231pd %%zmm14,%%zmm0,%%zmm1	\n\t vfmadd231pd %%zmm14,%%zmm15,%%zmm6	\n\t"/* c1 += cc5*t5; */\
	"vfmadd231pd (%%r9 ),%%zmm0,%%zmm2	\n\t vfmadd231pd (%%r9 ),%%zmm15,%%zmm7	\n\t"/* c2 += cc1*t5; */\
	"vfmadd231pd %%zmm13,%%zmm0,%%zmm3	\n\t vfmadd231pd %%zmm13,%%zmm15,%%zmm8	\n\t"/* c3 += cc4*t5; */\
	"vfmadd231pd %%zmm11,%%zmm0,%%zmm4	\n\t vfmadd231pd %%zmm11,%%zmm15,%%zmm9	\n\t"/* c4 += cc2*t5; */\
	"vfmadd231pd %%zmm12,%%zmm0,%%zmm5	\n\t vfmadd231pd %%zmm12,%%zmm15,%%zmm10	\n\t"/* c5 += cc3*t5; */\
		"vaddpd (%%r8 ),%%zmm0,%%zmm0	\n\t	vaddpd 0x40(%%r8 ),%%zmm15,%%zmm15	\n\t"/* B0 += t5; */\
	/*
	Current GPR assignments:
		rax	rbx	rcx	rdx	rsi	rdi	r8	r9	r10	r11	r12	r13	r14	r15
		in	iff	out	off	i3	i4	i0	cc1	i5	o6	o7	o8	o9	oA
	*/\
		"movslq	(%%rdx),%%r8	\n\t	addq	%%rcx,%%r8	\n\t"/* o0-from-offset */\
		"vmovaps	%%zmm0,(%%r8 )		\n\t	vmovaps	%%zmm15,0x40(%%r8 )			\n\t"/* Store B0, now into output slot */\
		/* Have enough registers to also hold the real parts of S-terms: */\
		"vmovaps	(%%r15),%%zmm11		\n\t"/* s1r */\
		"vmovaps	(%%r14),%%zmm12		\n\t"/* s2r */\
		"vmovaps	(%%r13),%%zmm13		\n\t"/* s3r */\
		"vmovaps	(%%r12),%%zmm14		\n\t"/* s4r */\
		"vmovaps	(%%r11),%%zmm15		\n\t"/* s5r */\
		/* Need rcx,rdx for O-address computations, but done with I-addesses so rax,rbx free: */\
		"movslq	0x04(%%rdx),%%rsi	\n\t	addq	%%rcx,%%rsi		\n\t"/* &B1 */\
		"movslq	0x08(%%rdx),%%rdi	\n\t	addq	%%rcx,%%rdi		\n\t"/* &B2 */\
		"movslq	0x0c(%%rdx),%%rax	\n\t	addq	%%rcx,%%rax		\n\t"/* &B3 */\
		"movslq	0x10(%%rdx),%%rbx	\n\t	addq	%%rcx,%%rbx		\n\t"/* &B4 */\
		"movslq	0x14(%%rdx),%%r10	\n\t	addq	%%rcx,%%r10		\n\t"/* &B5 */\
	/* ...Do the lower-half-index portion of the closing wave of radix-2 butterflies: */\
		"vsubpd 0x40(%%r15),%%zmm1,%%zmm1		\n\t	 vaddpd %%zmm11,%%zmm6 ,%%zmm6	\n\t"/* B1 = C1 + I*S1 */\
		"vsubpd 0x40(%%r14),%%zmm2,%%zmm2		\n\t	 vaddpd %%zmm12,%%zmm7 ,%%zmm7	\n\t"/* B2 = C2 + I*S2 */\
		"vsubpd 0x40(%%r13),%%zmm3,%%zmm3		\n\t	 vaddpd %%zmm13,%%zmm8 ,%%zmm8	\n\t"/* B3 = C3 + I*S3 */\
		"vsubpd 0x40(%%r12),%%zmm4,%%zmm4		\n\t	 vaddpd %%zmm14,%%zmm9 ,%%zmm9	\n\t"/* B4 = C4 + I*S4 */\
		"vsubpd 0x40(%%r11),%%zmm5,%%zmm5		\n\t	 vaddpd %%zmm15,%%zmm10,%%zmm10	\n\t"/* B5 = C5 + I*S5 */\
		"movq $0x4000000000000000,%%r8	\n\t	vpbroadcastq  %%r8,%%zmm0	\n\t"/* 2.0 */\
	/* Write lower-half outputs back to memory... */\
		"vmovaps	%%zmm1,(%%rsi)		\n\t"/* B1r */"		vmovaps	%%zmm6 ,0x40(%%rsi)	\n\t"/* B1i */\
		"vmovaps	%%zmm2,(%%rdi)		\n\t"/* B2r */"		vmovaps	%%zmm7 ,0x40(%%rdi)	\n\t"/* B2i */\
		"vmovaps	%%zmm3,(%%rax)		\n\t"/* B3r */"		vmovaps	%%zmm8 ,0x40(%%rax)	\n\t"/* B3i */\
		"vmovaps	%%zmm4,(%%rbx)		\n\t"/* B4r */"		vmovaps	%%zmm9 ,0x40(%%rbx)	\n\t"/* B4i */\
		"vmovaps	%%zmm5,(%%r10)		\n\t"/* B5r */"		vmovaps	%%zmm10,0x40(%%r10)	\n\t"/* B5i */\
	/* ...Do the upper-half-index portion of the radix-2 butterflies: */\
		" vfmadd231pd 0x40(%%r15),%%zmm0,%%zmm1		\n\t	vfnmadd231pd %%zmm11,%%zmm0,%%zmm6	\n\t"/* B10= C1 - I*S1 */\
		" vfmadd231pd 0x40(%%r14),%%zmm0,%%zmm2		\n\t	vfnmadd231pd %%zmm12,%%zmm0,%%zmm7	\n\t"/* B9 = C2 - I*S2 */\
		" vfmadd231pd 0x40(%%r13),%%zmm0,%%zmm3		\n\t	vfnmadd231pd %%zmm13,%%zmm0,%%zmm8	\n\t"/* B8 = C3 - I*S3 */\
		" vfmadd231pd 0x40(%%r12),%%zmm0,%%zmm4		\n\t	vfnmadd231pd %%zmm14,%%zmm0,%%zmm9	\n\t"/* B7 = C4 - I*S4 */\
		" vfmadd231pd 0x40(%%r11),%%zmm0,%%zmm5		\n\t	vfnmadd231pd %%zmm15,%%zmm0,%%zmm10	\n\t"/* B6 = C5 - I*S5 */\
	/* ...And write the upper-half outputs back to memory. */\
		"vmovaps	%%zmm1,(%%r15)		\n\t"/* Bar */"		vmovaps	%%zmm6 ,0x40(%%r15)	\n\t"/* Bai */\
		"vmovaps	%%zmm2,(%%r14)		\n\t"/* B9r */"		vmovaps	%%zmm7 ,0x40(%%r14)	\n\t"/* B9i */\
		"vmovaps	%%zmm3,(%%r13)		\n\t"/* B8r */"		vmovaps	%%zmm8 ,0x40(%%r13)	\n\t"/* B8i */\
		"vmovaps	%%zmm4,(%%r12)		\n\t"/* B7r */"		vmovaps	%%zmm9 ,0x40(%%r12)	\n\t"/* B7i */\
		"vmovaps	%%zmm5,(%%r11)		\n\t"/* B6r */"		vmovaps	%%zmm10,0x40(%%r11)	\n\t"/* B6i */\
		:					/* outputs: none */\
		: [in0] "m" (Xin0)	/* All inputs from memory addresses here */\
		 ,[i_off] "m" (Xi_off)	/* I-address pointer-strides */\
		 ,[__cc] "m" (Xcc)\
		 ,[out0] "m" (Xout0)\
		 ,[o_off] "m" (Xo_off)	/* O-address pointer-strides */\
		: "cc","memory","rax","rbx","rcx","rdx","rsi","rdi","r8","r9","r10","r11","r12","r13","r14","r15",\
			"xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
		);\
	}

   #else		// [2] [LOACC, no-flags default] based on Radix-11 DFT impl based on RADIX_11_DFT in dft_macro.h,
				// using length-5 cyclic convolution scheme for the 5 x 5 matrix submultiplies.
				// Totals: [134 ADD, 54 MUL/FMA, 153 memref] in this FMA-lite version (Use extra AVX512 vregs to save 19 memref vs AVX2 version).
				//
	#define SSE2_RADIX_11_DFT(XI0,Xi1,Xi2,Xi3,Xi4,Xi5,Xi6,Xi7,Xi8,Xi9,XiA, Xcc, XO0,Xo1,Xo2,Xo3,Xo4,Xo5,Xo6,Xo7,Xo8,Xo9,XoA)\
	{\
	__asm__ volatile (\
		/********************************************/\
		/*       Here are the 5 cosine terms:       */\
		/********************************************/\
		"movq	%[__I0],%%rax			\n\t			movq	%[__O0],%%rcx		\n\t"/* rax/rcx point to Re parts of I/Os */\
		"leaq	0x40(%%rax),%%rbx		\n\t			leaq	0x40(%%rcx),%%rdx	\n\t"/* rbx/rdx point to Im parts of I/Os */\
		"movq	%[__cc],%%rsi			\n\t	vmovaps -0x80(%%rsi),%%zmm31	\n\t"/* cc1 ... two in -0x80(%%rsi) */\
		/*****************************/				/******************************/\
		/*        Real Parts:        */				/*     Imaginary Parts:       */\
		/*****************************/				/******************************/\
	"movq	%[__i1],%%r8	\n\t	vmovaps (%%r8 ),%%zmm1	\n\t	vmovaps	0x40(%%r8 ),%%zmm9	\n\t"\
	"movq	%[__iA],%%r9	\n\t	vmovaps (%%r9 ),%%zmm5	\n\t	vmovaps	0x40(%%r9 ),%%zmm13	\n\t"\
	"movq	%[__i2],%%r10	\n\t	vmovaps (%%r10),%%zmm2	\n\t	vmovaps	0x40(%%r10),%%zmm10	\n\t"\
	"movq	%[__i9],%%r11	\n\t	vmovaps (%%r11),%%zmm6	\n\t	vmovaps	0x40(%%r11),%%zmm14	\n\t"\
	"movq	%[__i3],%%r12	\n\t	vmovaps (%%r12),%%zmm3	\n\t	vmovaps	0x40(%%r12),%%zmm11	\n\t"\
	"movq	%[__i8],%%r13	\n\t	vmovaps (%%r13),%%zmm7	\n\t	vmovaps	0x40(%%r13),%%zmm15	\n\t"\
	"movq	%[__i4],%%r14	\n\t	vmovaps (%%r14),%%zmm4	\n\t	vmovaps	0x40(%%r14),%%zmm12	\n\t"\
	"movq	%[__i7],%%r15	\n\t	vmovaps (%%r15),%%zmm0	\n\t	vmovaps	0x40(%%r15),%%zmm8	\n\t"\
		"vsubpd	%%zmm5,%%zmm1,%%zmm1	\n\t			vsubpd	%%zmm13,%%zmm9 ,%%zmm9			\n\t"/* t21 = x1-xA */\
		"vsubpd	%%zmm6,%%zmm2,%%zmm2	\n\t			vsubpd	%%zmm14,%%zmm10,%%zmm10			\n\t"/* t19 = x2-x9 */\
		"vsubpd	%%zmm7,%%zmm3,%%zmm3	\n\t			vsubpd	%%zmm15,%%zmm11,%%zmm11			\n\t"/* t17 = x3-x8 */\
		"vsubpd	%%zmm0,%%zmm4,%%zmm4	\n\t			vsubpd	%%zmm8 ,%%zmm12,%%zmm12			\n\t"/* t15 = x4-x7 */\
	"movq	%[__oA],%%r8	\n\t	vmovaps %%zmm1,0x40(%%r8 )	\n\t	vmovaps %%zmm9 ,(%%r8 )	\n\t"/* ...and store the above 4 pairs */\
	"movq	%[__o9],%%r9	\n\t	vmovaps %%zmm2,0x40(%%r9 )	\n\t	vmovaps %%zmm10,(%%r9 )	\n\t"/* of re/im-part diffs back to mem */\
	"movq	%[__o8],%%r10	\n\t	vmovaps %%zmm3,0x40(%%r10)	\n\t	vmovaps %%zmm11,(%%r10)	\n\t"/* slots oA,8,8,7, swapping the re/im. */\
	"movq	%[__o7],%%r11	\n\t	vmovaps %%zmm4,0x40(%%r11)	\n\t	vmovaps %%zmm12,(%%r11)	\n\t"\
	"vfmadd231pd 	%%zmm31,%%zmm5,%%zmm1	\n\t	vfmadd231pd 	%%zmm31,%%zmm13,%%zmm9	\n\t"/* t3  = x1+xA */\
	"vfmadd231pd 	%%zmm31,%%zmm6,%%zmm2	\n\t	vfmadd231pd 	%%zmm31,%%zmm14,%%zmm10	\n\t"/* t5  = x2+x9 */\
	"movq	%[__i5],%%r9	\n\t	vmovaps (%%r9 ),%%zmm5	\n\t	vmovaps 0x40(%%r9 ),%%zmm13	\n\t"/* x5 */\
	"movq	%[__i6],%%r11	\n\t	vmovaps (%%r11),%%zmm6	\n\t	vmovaps 0x40(%%r11),%%zmm14	\n\t"/* x6 */\
	"vfmadd231pd 	%%zmm31,%%zmm7,%%zmm3	\n\t	vfmadd231pd 	%%zmm31,%%zmm15,%%zmm11	\n\t"/* t7  = x3+x8 */\
	"vfmadd231pd 	%%zmm31,%%zmm0,%%zmm4	\n\t	vfmadd231pd 	%%zmm31,%%zmm8 ,%%zmm12	\n\t"/* t9  = x4+x7 */\
		"vsubpd	%%zmm6,%%zmm5,%%zmm5	\n\t			vsubpd	%%zmm14,%%zmm13,%%zmm13			\n\t"/* t13 = x5-x6 */\
	"movq	%[__o6],%%r8	\n\t	vmovaps %%zmm5,0x40(%%r8 )	\n\t	vmovaps %%zmm13,(%%r8 )	\n\t"/* store t13 = (x5-x6) in o6 */\
		"vmovaps	(%%rax),%%zmm0		\n\t			vmovaps	(%%rbx),%%zmm8					\n\t"/* x0 */\
	"vfmadd231pd 	%%zmm31,%%zmm6,%%zmm5	\n\t	vfmadd231pd 	%%zmm31,%%zmm14,%%zmm13	\n\t"/* t11 = x5+x6 */\
	"movq	%[__o1],%%r11	\n\t"\
	"movq	%[__o2],%%r12	\n\t"\
	"movq	%[__o3],%%r13	\n\t"\
	"movq	%[__o4],%%r14	\n\t"\
	"movq	%[__o5],%%r15	\n\t"\
		"vsubpd	%%zmm2,%%zmm1,%%zmm1			\n\t	vsubpd	%%zmm10,%%zmm9 ,%%zmm9			\n\t"/* c1 = t3 -t5 */\
		"vsubpd	%%zmm2,%%zmm5,%%zmm5			\n\t	vsubpd	%%zmm10,%%zmm13,%%zmm13			\n\t"/* c2 = t11-t5 */\
		"vmovaps	%%zmm1,%%zmm6				\n\t	vmovaps	%%zmm9 ,%%zmm14					\n\t"/* cpy c1 */\
		"vsubpd	%%zmm2,%%zmm3,%%zmm3			\n\t	vsubpd	%%zmm10,%%zmm11,%%zmm11			\n\t"/* c4 = t7 -t5 */\
		"vmovaps	%%zmm1,%%zmm7				\n\t	vmovaps	%%zmm9 ,%%zmm15					\n\t"/* cpy c1 */\
		"vmulpd	      (%%rsi),%%zmm1,%%zmm1		\n\t	vmulpd	      (%%rsi),%%zmm9 ,%%zmm9	\n\t"/* a0.c1 */\
		"vmovaps	%%zmm1,(%%r11)				\n\t	vmovaps	%%zmm9 ,0x40(%%r11)				\n\t"/* store copy of a0.c1 in mem */\
		"vsubpd	%%zmm2,%%zmm4,%%zmm4			\n\t	vsubpd	%%zmm10,%%zmm12,%%zmm12			\n\t"/* c5 = t9-t5 */\
		"vaddpd	%%zmm5,%%zmm6,%%zmm6			\n\t	vaddpd	%%zmm13,%%zmm14,%%zmm14			\n\t"/* c3 = c1+c2 */\
		"vsubpd	%%zmm3,%%zmm7,%%zmm7			\n\t	vsubpd	%%zmm11,%%zmm15,%%zmm15			\n\t"/* c7 = c1-c4 */\
		"vmulpd	 0x180(%%rsi),%%zmm7,%%zmm7		\n\t	vmulpd	 0x180(%%rsi),%%zmm15,%%zmm15	\n\t"/* a6.c7 */\
		"vmovaps	%%zmm7,(%%r13)				\n\t	vmovaps	%%zmm15,0x40(%%r13)				\n\t"/* store copy of a6.c7 in mem */\
		"vmovaps	%%zmm3,%%zmm1				\n\t	vmovaps	%%zmm11,%%zmm9					\n\t"/* cpy c4 */\
		"vaddpd	%%zmm4,%%zmm3,%%zmm3			\n\t	vaddpd	%%zmm12,%%zmm11,%%zmm11			\n\t"/* c6 = c4+c5 */\
		"vmulpd	 0x0c0(%%rsi),%%zmm1,%%zmm1		\n\t	vmulpd	 0x0c0(%%rsi),%%zmm9 ,%%zmm9	\n\t"/* a3.c4 */\
		"vmovaps	%%zmm1,(%%r12)				\n\t	vmovaps	%%zmm9 ,0x40(%%r12)				\n\t"/* store copy of a3.c4 in mem */\
		"vmovaps	%%zmm5,%%zmm7				\n\t	vmovaps	%%zmm13,%%zmm15					\n\t"/* cpy c2 */\
		"vsubpd	%%zmm4,%%zmm5,%%zmm5			\n\t	vsubpd	%%zmm12,%%zmm13,%%zmm13			\n\t"/* c8 = c2-c5 */\
	"vfmadd132pd -0x40(%%rsi),%%zmm3,%%zmm2		\n\t vfmadd132pd -0x40(%%rsi),%%zmm11,%%zmm10	\n\t"/* 5.t5+c6 */\
		"vmovaps	%%zmm6,%%zmm1				\n\t	vmovaps	%%zmm14,%%zmm9					\n\t"/* cpy c3 */\
		"vsubpd	%%zmm3,%%zmm6,%%zmm6			\n\t	vsubpd	%%zmm11,%%zmm14,%%zmm14			\n\t"/* c9 = c3-c6 */\
		"vmulpd	 0x200(%%rsi),%%zmm6,%%zmm6		\n\t	vmulpd	 0x200(%%rsi),%%zmm14,%%zmm14	\n\t"/* c9 = a8.c9 */\
		"vmulpd	 0x140(%%rsi),%%zmm3,%%zmm3		\n\t	vmulpd	 0x140(%%rsi),%%zmm11,%%zmm11	\n\t"/* c6 = a5.c6 */\
		"vaddpd	%%zmm1,%%zmm2,%%zmm2			\n\t	vaddpd	%%zmm9 ,%%zmm10,%%zmm10			\n\t"/* c10= c3 +c6+5*t5 */\
		"vmulpd	 0x080(%%rsi),%%zmm1,%%zmm1		\n\t	vmulpd	 0x080(%%rsi),%%zmm9 ,%%zmm9	\n\t"/* c3 = a2.c3 */\
		"vaddpd	%%zmm2,%%zmm0,%%zmm0			\n\t	vaddpd	%%zmm10,%%zmm8 ,%%zmm8			\n\t"/* B0 = A0+c10 */\
		"vmovaps	%%zmm0,(%%rcx)				\n\t	vmovaps	%%zmm8 ,(%%rdx)					\n\t"/* store B0 */\
	"vfmadd132pd 0x040(%%rsi),%%zmm1,%%zmm7		\n\t vfmadd132pd 0x040(%%rsi),%%zmm9 ,%%zmm15	\n\t"/* c2 = a1.c2+c3 */\
	"vfmadd132pd 0x100(%%rsi),%%zmm3,%%zmm4		\n\t vfmadd132pd 0x100(%%rsi),%%zmm11,%%zmm12	\n\t"/* c5 = a4.c5+c6 */\
	"vfmadd132pd 0x1c0(%%rsi),%%zmm6,%%zmm5		\n\t vfmadd132pd 0x1c0(%%rsi),%%zmm14,%%zmm13	\n\t"/* c8 = a7.c8+c9 */\
	"vfmadd132pd 0x240(%%rsi),%%zmm0,%%zmm2		\n\t vfmadd132pd 0x240(%%rsi),%%zmm8 ,%%zmm10	\n\t"/* c10= a9.c10+B0 */\
		"vaddpd	(%%r11),%%zmm1,%%zmm1			\n\t	vaddpd	0x40(%%r11),%%zmm9 ,%%zmm9		\n\t"/* c1 = a0.c1+c3 */\
		"vaddpd	(%%r12),%%zmm3,%%zmm3			\n\t	vaddpd	0x40(%%r12),%%zmm11,%%zmm11		\n\t"/* c4 = a3.c4+c6 */\
		"vaddpd	(%%r13),%%zmm6,%%zmm6			\n\t	vaddpd	0x40(%%r13),%%zmm14,%%zmm14		\n\t"/* c7 = a6.c7+c9 */\
		"vmovaps	%%zmm2,%%zmm0				\n\t	vmovaps	%%zmm10,%%zmm8					\n\t"/* cpy c10 */\
		"vsubpd	%%zmm1,%%zmm2,%%zmm2			\n\t	vsubpd	%%zmm9 ,%%zmm10,%%zmm10			\n\t"/* c10-c1 */\
		"vaddpd	%%zmm0,%%zmm1,%%zmm1			\n\t	vaddpd	%%zmm8 ,%%zmm9 ,%%zmm9			\n\t"/* c10+c1 */\
		"vsubpd	%%zmm7,%%zmm2,%%zmm2			\n\t	vsubpd	%%zmm15,%%zmm10,%%zmm10			\n\t"/* c10-c1-c2 */\
		"vaddpd	%%zmm0,%%zmm7,%%zmm7			\n\t	vaddpd	%%zmm8 ,%%zmm15,%%zmm15			\n\t"/* c10+c2 */\
		"vsubpd	%%zmm6,%%zmm1,%%zmm1			\n\t	vsubpd	%%zmm14,%%zmm9 ,%%zmm9			\n\t"/* cr1 = c10+c1-c7 */\
		"vsubpd	%%zmm3,%%zmm2,%%zmm2			\n\t	vsubpd	%%zmm11,%%zmm10,%%zmm10			\n\t"/* c10-c1-c2-c4 */\
		"vaddpd	%%zmm0,%%zmm3,%%zmm3			\n\t	vaddpd	%%zmm8 ,%%zmm11,%%zmm11			\n\t"/* c10+c4 */\
		"vsubpd	%%zmm5,%%zmm7,%%zmm7			\n\t	vsubpd	%%zmm13,%%zmm15,%%zmm15			\n\t"/* cr5 = c10+c2-c8 */\
		"vsubpd	%%zmm4,%%zmm2,%%zmm2			\n\t	vsubpd	%%zmm12,%%zmm10,%%zmm10			\n\t"/* cr2 = c10-c1-c2-c4-c5 */\
		"vaddpd	%%zmm0,%%zmm4,%%zmm4			\n\t	vaddpd	%%zmm8 ,%%zmm12,%%zmm12			\n\t"/* c10+c5 */\
		"vaddpd	%%zmm6,%%zmm3,%%zmm3			\n\t	vaddpd	%%zmm14,%%zmm11,%%zmm11			\n\t"/* cr3 = c10+c4+c7 */\
		"vaddpd	%%zmm5,%%zmm4,%%zmm4			\n\t	vaddpd	%%zmm13,%%zmm12,%%zmm12			\n\t"/* cr4 = c10+c5+c8 */\
		"vmovaps	%%zmm1,(%%r11)				\n\t	vmovaps	%%zmm9 ,0x40(%%r11)				\n\t"/* write cr1 */\
		"vmovaps	%%zmm7,(%%r15)				\n\t	vmovaps	%%zmm15,0x40(%%r15)				\n\t"/* write cr5 */\
		"vmovaps	%%zmm2,(%%r12)				\n\t	vmovaps	%%zmm10,0x40(%%r12)				\n\t"/* write cr2 */\
		"vmovaps	%%zmm3,(%%r13)				\n\t	vmovaps	%%zmm11,0x40(%%r13)				\n\t"/* write cr3 */\
		"vmovaps	%%zmm4,(%%r14)				\n\t	vmovaps	%%zmm12,0x40(%%r14)				\n\t"/* write cr4 */\
	/***********************************************************************************************************/\
	/* Here are 5 sine terms: Similar to cosines, but t21,19,17,15,13 replace t3,5,7,9,11 and some sign flips: */\
	/***********************************************************************************************************/\
		"addq	$0x280,%%rsi				\n\t"/* b0 */\
		/*****************************/				/******************************/\
		/*     Imaginary Parts:      */				/*        Real Parts:         *//*** Note re/im swap ***/\
		/*****************************/				/******************************/\
	"movq	%[__oA],%%r11	\n\t	vmovaps (%%r11),%%zmm1	\n\t	vmovaps 0x40(%%r11),%%zmm9	\n\t"/* t21 */\
	"movq	%[__o9],%%r12	\n\t	vmovaps (%%r12),%%zmm2	\n\t	vmovaps 0x40(%%r12),%%zmm10	\n\t"/* t19 */\
	"movq	%[__o8],%%r13	\n\t	vmovaps (%%r13),%%zmm3	\n\t	vmovaps 0x40(%%r13),%%zmm11	\n\t"/* t17 */\
	"movq	%[__o7],%%r14	\n\t	vmovaps (%%r14),%%zmm4	\n\t	vmovaps 0x40(%%r14),%%zmm12	\n\t"/* t15 */\
	"movq	%[__o6],%%r15	\n\t	vmovaps (%%r15),%%zmm5	\n\t	vmovaps 0x40(%%r15),%%zmm13	\n\t"/* t13 */\
		"vaddpd	%%zmm2,%%zmm1,%%zmm1			\n\t	vaddpd	%%zmm10,%%zmm9 ,%%zmm9			\n\t"/* c1 = t21+t19 */\
		"vaddpd	%%zmm2,%%zmm5,%%zmm5			\n\t	vaddpd	%%zmm10,%%zmm13,%%zmm13			\n\t"/* c2 = t13+t19 */\
		"vmovaps	%%zmm1,%%zmm6				\n\t	vmovaps	%%zmm9 ,%%zmm14					\n\t"/* cpy c1 */\
		"vaddpd	%%zmm2,%%zmm3,%%zmm3			\n\t	vaddpd	%%zmm10,%%zmm11,%%zmm11			\n\t"/* c4 = t17+t19 */\
		"vmovaps	%%zmm1,%%zmm7				\n\t	vmovaps	%%zmm9 ,%%zmm15					\n\t"/* cpy c1 */\
		"vmulpd	      (%%rsi),%%zmm1,%%zmm1		\n\t	vmulpd	      (%%rsi),%%zmm9 ,%%zmm9	\n\t"/* b0.c1 */\
		"vmovaps	%%zmm1,(%%r15)				\n\t	vmovaps	%%zmm9 ,0x40(%%r15)				\n\t"/* store b0.c1 */\
		"vaddpd	%%zmm2,%%zmm4,%%zmm4			\n\t	vaddpd	%%zmm10,%%zmm12,%%zmm12			\n\t"/* c5 = t15+t19 */\
		"vaddpd	%%zmm5,%%zmm6,%%zmm6			\n\t	vaddpd	%%zmm13,%%zmm14,%%zmm14			\n\t"/* c3 = c1+c2 */\
		"vsubpd	%%zmm3,%%zmm7,%%zmm7			\n\t	vsubpd	%%zmm11,%%zmm15,%%zmm15			\n\t"/* c7 = c1-c4 */\
		"vmulpd	 0x180(%%rsi),%%zmm7,%%zmm7		\n\t	vmulpd	 0x180(%%rsi),%%zmm15,%%zmm15	\n\t"/* b6.c7 */\
		"vmovaps	%%zmm7,(%%r13)				\n\t	vmovaps	%%zmm15,0x40(%%r13)				\n\t"/* store b6.c7 */\
		"vmovaps	%%zmm3,%%zmm1				\n\t	vmovaps	%%zmm11,%%zmm9					\n\t"/* cpy c4 */\
		"vaddpd	%%zmm4,%%zmm3,%%zmm3			\n\t	vaddpd	%%zmm12,%%zmm11,%%zmm11			\n\t"/* c6 = c4+c5 */\
		"vmulpd	 0x0c0(%%rsi),%%zmm1,%%zmm1		\n\t	vmulpd	 0x0c0(%%rsi),%%zmm9 ,%%zmm9	\n\t"/* b3.c4 */\
		"vmovaps	%%zmm1,(%%r14)				\n\t	vmovaps	%%zmm9 ,0x40(%%r14)				\n\t"/* store b3.c4 */\
		"vmovaps	%%zmm5,%%zmm7				\n\t	vmovaps	%%zmm13,%%zmm15					\n\t"/* cpy c2 */\
		"vsubpd	%%zmm4,%%zmm5,%%zmm5			\n\t	vsubpd	%%zmm12,%%zmm13,%%zmm13			\n\t"/* c8 = c2-c5 */\
	"vfmsub132pd -0x2c0(%%rsi),%%zmm3,%%zmm2	\n\t vfmsub132pd -0x2c0(%%rsi),%%zmm11,%%zmm10	\n\t"/* 5.t19-c6 */\
		"vmovaps	%%zmm6,%%zmm1				\n\t	vmovaps	%%zmm14,%%zmm9					\n\t"/* cpy c3 */\
		"vsubpd	%%zmm3,%%zmm6,%%zmm6			\n\t	vsubpd	%%zmm11,%%zmm14,%%zmm14			\n\t"/* c9 = c3-c6 */\
		"vmulpd	 0x200(%%rsi),%%zmm6,%%zmm6		\n\t	vmulpd	 0x200(%%rsi),%%zmm14,%%zmm14	\n\t"/* c9 = b8.c9 */\
		"vmulpd	 0x140(%%rsi),%%zmm3,%%zmm3		\n\t	vmulpd	 0x140(%%rsi),%%zmm11,%%zmm11	\n\t"/* c6 = b5.c6 */\
		"vsubpd	%%zmm1,%%zmm2,%%zmm2			\n\t	vsubpd	%%zmm9 ,%%zmm10,%%zmm10			\n\t"/* -c10 = 5.t19-c6-c3 (c6 prior to above mul-by-b5) */\
		"vmulpd	 0x080(%%rsi),%%zmm1,%%zmm1		\n\t	vmulpd	 0x080(%%rsi),%%zmm9 ,%%zmm9	\n\t"/* c3 = b2.c3 */\
		"vmulpd	 0x240(%%rsi),%%zmm2,%%zmm2		\n\t	vmulpd	 0x240(%%rsi),%%zmm10,%%zmm10	\n\t"/* c10 = -b9.-c10 */\
	"vfmadd132pd 0x1c0(%%rsi),%%zmm6,%%zmm5		\n\t vfmadd132pd 0x1c0(%%rsi),%%zmm14,%%zmm13	\n\t"/* c6 = b7.c8+c9 */\
		"vaddpd	(%%r13),%%zmm6,%%zmm6			\n\t	vaddpd	0x40(%%r13),%%zmm14,%%zmm14		\n\t"/* c5 = b6.c7+c9 */\
	"vfmadd132pd 0x100(%%rsi),%%zmm3,%%zmm4		\n\t vfmadd132pd 0x100(%%rsi),%%zmm11,%%zmm12	\n\t"/* c4 = b4.c5+c6 */\
		"vaddpd	(%%r14),%%zmm3,%%zmm3			\n\t	vaddpd	0x40(%%r14),%%zmm11,%%zmm11		\n\t"/* c3 = b3.c4+c6 */\
	"vfmadd132pd 0x040(%%rsi),%%zmm1,%%zmm7		\n\t vfmadd132pd 0x040(%%rsi),%%zmm9 ,%%zmm15	\n\t"/* c2 = b1.c2+c3 */\
		"vaddpd	(%%r15),%%zmm1,%%zmm1			\n\t	vaddpd	0x40(%%r15),%%zmm9 ,%%zmm9		\n\t"/* c1 = b0.c1+c3 */\
		"vpxorq	%%zmm0,%%zmm0,%%zmm0			\n\t	vpxorq	%%zmm8 ,%%zmm8 ,%%zmm8			\n\t"/* 0 */\
		"vsubpd	%%zmm2,%%zmm0,%%zmm0			\n\t	vsubpd	%%zmm10,%%zmm8 ,%%zmm8			\n\t"/* -c10 */\
		"vaddpd	%%zmm1,%%zmm0,%%zmm0			\n\t	vaddpd	%%zmm9 ,%%zmm8 ,%%zmm8			\n\t"/* c1-c10 */\
		"vaddpd	%%zmm2,%%zmm1,%%zmm1			\n\t	vaddpd	%%zmm10,%%zmm9 ,%%zmm9			\n\t"/* c10+c1 */\
		"vaddpd	%%zmm7,%%zmm0,%%zmm0			\n\t	vaddpd	%%zmm15,%%zmm8 ,%%zmm8			\n\t"/* c1+c2-c10 */\
		"vaddpd	%%zmm2,%%zmm7,%%zmm7			\n\t	vaddpd	%%zmm10,%%zmm15,%%zmm15			\n\t"/* c10+c2 */\
		"vsubpd	%%zmm6,%%zmm1,%%zmm1			\n\t	vsubpd	%%zmm14,%%zmm9 ,%%zmm9			\n\t"/* si1 = c10+c1-c5 */\
		"vaddpd	%%zmm3,%%zmm0,%%zmm0			\n\t	vaddpd	%%zmm11,%%zmm8 ,%%zmm8			\n\t"/* c1+c2+c3-c10 */\
		"vaddpd	%%zmm2,%%zmm3,%%zmm3			\n\t	vaddpd	%%zmm10,%%zmm11,%%zmm11			\n\t"/* c10+c3 */\
		"vsubpd	%%zmm5,%%zmm7,%%zmm7			\n\t	vsubpd	%%zmm13,%%zmm15,%%zmm15			\n\t"/* si5 = c10+c2-c6 */\
		"vaddpd	%%zmm4,%%zmm0,%%zmm0			\n\t	vaddpd	%%zmm12,%%zmm8 ,%%zmm8			\n\t"/* si2 = c1+c2+c3+c4-c10 */\
		"vaddpd	%%zmm2,%%zmm4,%%zmm4			\n\t	vaddpd	%%zmm10,%%zmm12,%%zmm12			\n\t"/* c10+c4 */\
		"vaddpd	%%zmm6,%%zmm3,%%zmm3			\n\t	vaddpd	%%zmm14,%%zmm11,%%zmm11			\n\t"/* si3 = c10+c3+c5 */\
		"vaddpd	%%zmm5,%%zmm4,%%zmm4			\n\t	vaddpd	%%zmm13,%%zmm12,%%zmm12			\n\t"/* si4 = c10+c4+c6 */\
	/* si[1-5] in mm1,0,3,4,7; sr[1-5] in mm9,8,11,12,15: */\
	"movq	%[__o1],%%rax					\n\t"\
	"movq	%[__o2],%%rbx					\n\t"\
	"movq	%[__o3],%%r8					\n\t"\
	"movq	%[__o4],%%r9					\n\t"\
	"movq	%[__o5],%%r10					\n\t"\
		"vmovaps		(%%rax),%%zmm5 		\n\t"/* cr1 */\
		"vmovaps		(%%rbx),%%zmm6 		\n\t"/* cr2 */\
		"vmovaps		(%%r8 ),%%zmm10		\n\t"/* cr3 */\
		"vmovaps		(%%r9 ),%%zmm13		\n\t"/* cr4 */\
		"vmovaps		(%%r10),%%zmm14		\n\t"/* cr5 */\
		"vsubpd	%%zmm1 ,%%zmm5 ,%%zmm16			\n\t	vaddpd	%%zmm5 ,%%zmm1,%%zmm1	\n\t"/* B1r = cr1-si1;  Bar = cr1+si1 */\
		"vsubpd	%%zmm0 ,%%zmm6 ,%%zmm17			\n\t	vaddpd	%%zmm6 ,%%zmm0,%%zmm0	\n\t"/* B2r = cr2-si2;  B9r = cr2+si2 */\
		"vsubpd	%%zmm3 ,%%zmm10,%%zmm18			\n\t	vaddpd	%%zmm10,%%zmm3,%%zmm3	\n\t"/* B3r = cr3-si3;  B8r = cr3+si3 */\
		"vsubpd	%%zmm4 ,%%zmm13,%%zmm19			\n\t	vaddpd	%%zmm13,%%zmm4,%%zmm4	\n\t"/* B4r = cr4-si4;  B7r = cr4+si4 */\
		"vsubpd	%%zmm7 ,%%zmm14,%%zmm20			\n\t	vaddpd	%%zmm14,%%zmm7,%%zmm7	\n\t"/* B5r = cr5-si5;  B6r = cr5+si5 */\
		"vmovaps	%%zmm16,    (%%rax)			\n\t	vmovaps	%%zmm1,    (%%r11)			\n\t"/* write B1r,Bar */\
		"vmovaps	%%zmm17,    (%%rbx)			\n\t	vmovaps	%%zmm0,    (%%r12)			\n\t"/* write B2r,B9r */\
		"vmovaps	%%zmm18,    (%%r8 )			\n\t	vmovaps	%%zmm3,    (%%r13)			\n\t"/* write B3r,B8r */\
		"vmovaps	%%zmm19,    (%%r9 )			\n\t	vmovaps	%%zmm4,    (%%r14)			\n\t"/* write B4r,B7r */\
		"vmovaps	%%zmm20,    (%%r10)			\n\t	vmovaps	%%zmm7,    (%%r15)			\n\t"/* write B5r,B6r */\
		"vmovaps	0x40(%%rax),%%zmm5		\n\t"/* ci1 */\
		"vmovaps	0x40(%%rbx),%%zmm6		\n\t"/* ci2 */\
		"vmovaps	0x40(%%r8 ),%%zmm10		\n\t"/* ci3 */\
		"vmovaps	0x40(%%r9 ),%%zmm13		\n\t"/* ci4 */\
		"vmovaps	0x40(%%r10),%%zmm14		\n\t"/* ci5 */\
		"vsubpd	%%zmm9 ,%%zmm5 ,%%zmm16			\n\t	vaddpd	%%zmm5 ,%%zmm9 ,%%zmm9 	\n\t"/* Bai = ci1-si1;  B1i = ci1+si1 */\
		"vsubpd	%%zmm8 ,%%zmm6 ,%%zmm17			\n\t	vaddpd	%%zmm6 ,%%zmm8 ,%%zmm8 	\n\t"/* B9i = ci2-si2;  B2i = ci2+si2 */\
		"vsubpd	%%zmm11,%%zmm10,%%zmm18			\n\t	vaddpd	%%zmm10,%%zmm11,%%zmm11	\n\t"/* B8i = ci3-si3;  B3i = ci3+si3 */\
		"vsubpd	%%zmm12,%%zmm13,%%zmm19			\n\t	vaddpd	%%zmm13,%%zmm12,%%zmm12	\n\t"/* B7i = ci4-si4;  B4i = ci4+si4 */\
		"vsubpd	%%zmm15,%%zmm14,%%zmm20			\n\t	vaddpd	%%zmm14,%%zmm15,%%zmm15	\n\t"/* B6i = ci5-si5;  B5i = ci5+si5 */\
		"vmovaps	%%zmm16,0x40(%%r11)			\n\t	vmovaps	%%zmm9 ,0x40(%%rax)			\n\t"/* write Bai,B1i */\
		"vmovaps	%%zmm17,0x40(%%r12)			\n\t	vmovaps	%%zmm8 ,0x40(%%rbx)			\n\t"/* write B9i,B2i */\
		"vmovaps	%%zmm18,0x40(%%r13)			\n\t	vmovaps	%%zmm11,0x40(%%r8 )			\n\t"/* write B8i,B3i */\
		"vmovaps	%%zmm19,0x40(%%r14)			\n\t	vmovaps	%%zmm12,0x40(%%r9 )			\n\t"/* write B7i,B4i */\
		"vmovaps	%%zmm20,0x40(%%r15)			\n\t	vmovaps	%%zmm15,0x40(%%r10)			\n\t"/* write B6i,B5i */\
		:					/* outputs: none */\
		: [__I0] "m" (XI0)	/* All inputs from memory addresses here */\
		 ,[__i1] "m" (Xi1)\
		 ,[__i2] "m" (Xi2)\
		 ,[__i3] "m" (Xi3)\
		 ,[__i4] "m" (Xi4)\
		 ,[__i5] "m" (Xi5)\
		 ,[__i6] "m" (Xi6)\
		 ,[__i7] "m" (Xi7)\
		 ,[__i8] "m" (Xi8)\
		 ,[__i9] "m" (Xi9)\
		 ,[__iA] "m" (XiA)\
		 ,[__cc] "m" (Xcc)\
		 ,[__O0] "m" (XO0)\
		 ,[__o1] "m" (Xo1)\
		 ,[__o2] "m" (Xo2)\
		 ,[__o3] "m" (Xo3)\
		 ,[__o4] "m" (Xo4)\
		 ,[__o5] "m" (Xo5)\
		 ,[__o6] "m" (Xo6)\
		 ,[__o7] "m" (Xo7)\
		 ,[__o8] "m" (Xo8)\
		 ,[__o9] "m" (Xo9)\
		 ,[__oA] "m" (XoA)\
		: "cc","memory","rax","rbx","rcx","rdx","rsi","r8","r9","r10","r11","r12","r13","r14","r15",\
			"xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15","xmm16","xmm17","xmm18","xmm19","xmm20","xmm31"	/* Clobbered registers */\
		);\
	}

   #endif	// ? DFT_11_FMA

#elif defined(USE_AVX2)	// AVX+FMA version, based on RADIX_11_DFT_FMA in dft_macro.h

   #if DFT_11_FMA	// [1] Naive but good-ROE and FMA-friendly impl based on RADIX_11_DFT_BASIC in dft_macro.h:

	// FMAs used for all arithmetic, including 'trivial' ones (one mult = 1.0) to replace ADD/SUB:
	//
	// Arithmetic opcount: [20 ADD, 130 FMA (incl 10 MUL), 166 memref], close to general target of 1 memref per vec_dbl arithmetic op.
	// potentially much lower cycle count (at a max theoretical rate of 2 FMA/cycle) than non-FMA version.
	//
	// Compare to non-FMA: [182 ADD, 44 MUL), 154 memref], less accurate and bottlenecked by 1 ADD/cycle max issue rate.
	// Compare to DFT_11_FMA = 0 version below: [134 ADD, 54 MUL/FMA, 172 memref].
	//
	#define SSE2_RADIX_11_DFT(Xin0,Xi_off, Xcc, Xout0,Xo_off)\
	{\
	__asm__ volatile (\
		/* Jan 2021: In reduced-#args revamp of this macro, Due to GPR-usage constraints, instead of
		keeping a ptr to vec_dbl*two in GPR r15 throughout, create 2.0 in low double, of an XMM, copy
		to an MMX via MOVDQ2Q, and copy-back-via-MOVQ2DQ/shuffle-propagate-to-rest-of-XMM as needed.
		Note that in AVX-512 we can create consts directly into GPRs and use VPBROADCAST to broadcast
		to a vector register, but in AVX2 and below the register-source version of this instruction is
		not available.
		***UPDATE***: Alas, the use-MMX trick is HUGELY SLOW on my Haswell CPU, restore load-from-memory.
		*/\
		"movq	%[__cc],%%r9	\n\t"/* cc1 */\
		"vmovaps	-0x40(%%r9),%%ymm0		\n\t"/* 2.0 */\
	"movq	%[in0],%%rax		\n\t	movq	%[out0],%%rcx	\n\t"/* rax,rcx hold to I/O base-pointers */\
	"movq	%[i_off],%%rbx		\n\t	movq	%[o_off],%%rdx	\n\t"/* rbx,rdx hold pointers to uint32[] I/O pointer-offset-stride arrays */\
		/********************************************/\
		/*       Here are the 5 cosine terms:       */\
		/********************************************/\
		/* Compute &i[1-A] from 32-bit offsets: */\
		"movslq	0x04(%%rbx),%%rsi	\n\t	addq	%%rax,%%rsi		\n\t"\
		"movslq	0x08(%%rbx),%%rdi	\n\t	addq	%%rax,%%rdi		\n\t"\
		"movslq	0x0c(%%rbx),%%r8 	\n\t	addq	%%rax,%%r8 		\n\t"\
		"movslq	0x10(%%rbx),%%r9 	\n\t	addq	%%rax,%%r9 		\n\t"\
		"movslq	0x14(%%rbx),%%r10	\n\t	addq	%%rax,%%r10		\n\t"\
		"movslq	0x18(%%rbx),%%r11	\n\t	addq	%%rax,%%r11		\n\t"\
		"movslq	0x1c(%%rbx),%%r12	\n\t	addq	%%rax,%%r12		\n\t"\
		"movslq	0x20(%%rbx),%%r13	\n\t	addq	%%rax,%%r13		\n\t"\
		"movslq	0x24(%%rbx),%%r14	\n\t	addq	%%rax,%%r14		\n\t"\
		"movslq	0x28(%%rbx),%%r15	\n\t	addq	%%rax,%%r15		\n\t"\
		"vmovaps	(%%rsi),%%ymm1		\n\t"/* A1r */"		vmovaps	0x20(%%rsi),%%ymm6	\n\t"/* A1i */\
		"vmovaps	(%%rdi),%%ymm2		\n\t"/* A2r */"		vmovaps	0x20(%%rdi),%%ymm7	\n\t"/* A2i */\
		"vmovaps	(%%r8 ),%%ymm3		\n\t"/* A3r */"		vmovaps	0x20(%%r8 ),%%ymm8	\n\t"/* A3i */\
		"vmovaps	(%%r9 ),%%ymm4		\n\t"/* A4r */"		vmovaps	0x20(%%r9 ),%%ymm9	\n\t"/* A4i */\
		"vmovaps	(%%r10),%%ymm5		\n\t"/* A5r */"		vmovaps	0x20(%%r10),%%ymm10	\n\t"/* A5i */\
		/* Have enough registers to also hold the real-upper-half inputs: */\
		"vmovaps	(%%r11),%%ymm11		\n\t"/* A6r */\
		"vmovaps	(%%r12),%%ymm12		\n\t"/* A7r */\
		"vmovaps	(%%r13),%%ymm13		\n\t"/* A8r */\
		"vmovaps	(%%r14),%%ymm14		\n\t"/* A9r */\
		"vmovaps	(%%r15),%%ymm15		\n\t"/* AAr */\
	/* ...Do the + parts of the opening wave of radix-2 butterflies: */\
		"vaddpd %%ymm15,%%ymm1,%%ymm1	\n\t	vaddpd 0x20(%%r15),%%ymm6 ,%%ymm6 	\n\t"/* t1 = A1 + Aa */\
		"vaddpd %%ymm14,%%ymm2,%%ymm2	\n\t	vaddpd 0x20(%%r14),%%ymm7 ,%%ymm7 	\n\t"/* t2 = A2 + A9 */\
		"vaddpd %%ymm13,%%ymm3,%%ymm3	\n\t	vaddpd 0x20(%%r13),%%ymm8 ,%%ymm8 	\n\t"/* t3 = A3 + A8 */\
		"vaddpd %%ymm12,%%ymm4,%%ymm4	\n\t	vaddpd 0x20(%%r12),%%ymm9 ,%%ymm9 	\n\t"/* t4 = A4 + A7 */\
		"vaddpd %%ymm11,%%ymm5,%%ymm5	\n\t	vaddpd 0x20(%%r11),%%ymm10,%%ymm10	\n\t"/* t5 = A5 + A6 */\
	/* Write lower-half outputs back to memory... */\
		"vmovaps	%%ymm1,(%%rsi)		\n\t"/* t1r */"		vmovaps	%%ymm6 ,0x20(%%rsi)	\n\t"/* t1i */\
		"vmovaps	%%ymm2,(%%rdi)		\n\t"/* t2r */"		vmovaps	%%ymm7 ,0x20(%%rdi)	\n\t"/* t2i */\
		"vmovaps	%%ymm3,(%%r8 )		\n\t"/* t3r */"		vmovaps	%%ymm8 ,0x20(%%r8 )	\n\t"/* t3i */\
		"vmovaps	%%ymm4,(%%r9 )		\n\t"/* t4r */"		vmovaps	%%ymm9 ,0x20(%%r9 )	\n\t"/* t4i */\
		"vmovaps	%%ymm5,(%%r10)		\n\t"/* t5r */"		vmovaps	%%ymm10,0x20(%%r10)	\n\t"/* t5i */\
	/* ...Do the - parts of the radix-2 butterflies: */\
	"vfnmadd231pd %%ymm15,%%ymm0,%%ymm1		\n\t	vfnmadd231pd 0x20(%%r15),%%ymm0,%%ymm6	\n\t"/* ta = A1 - Aa */\
	"vfnmadd231pd %%ymm14,%%ymm0,%%ymm2		\n\t	vfnmadd231pd 0x20(%%r14),%%ymm0,%%ymm7	\n\t"/* t9 = A2 - A9 */\
	"vfnmadd231pd %%ymm13,%%ymm0,%%ymm3		\n\t	vfnmadd231pd 0x20(%%r13),%%ymm0,%%ymm8	\n\t"/* t8 = A3 - A8 */\
	"vfnmadd231pd %%ymm12,%%ymm0,%%ymm4		\n\t	vfnmadd231pd 0x20(%%r12),%%ymm0,%%ymm9	\n\t"/* t7 = A4 - A7 */\
	"vfnmadd231pd %%ymm11,%%ymm0,%%ymm5		\n\t	vfnmadd231pd 0x20(%%r11),%%ymm0,%%ymm10	\n\t"/* t6 = A5 - A6 */\
	/* ...And write upper-half outputs 6-9 (A stays in-reg) back to memory to free up registers for the 5x5 sine-term-subconvo computation: */\
		"vmovaps	%%ymm2,(%%r14)		\n\t"/* t9r */"		vmovaps	%%ymm7 ,0x20(%%r14)	\n\t"/* t9i */\
		"vmovaps	%%ymm3,(%%r13)		\n\t"/* t8r */"		vmovaps	%%ymm8 ,0x20(%%r13)	\n\t"/* t8i */\
		"vmovaps	%%ymm4,(%%r12)		\n\t"/* t7r */"		vmovaps	%%ymm9 ,0x20(%%r12)	\n\t"/* t7i */\
		"vmovaps	%%ymm5,(%%r11)		\n\t"/* t6r */"		vmovaps	%%ymm10,0x20(%%r11)	\n\t"/* t6i */\
	/*
		Sr1 =      ss1*tar+ss2*t9r+ss3*t8r+ss4*t7r+ss5*t6r;		Si1 =      ss1*tai+ss2*t9i+ss3*t8i+ss4*t7i+ss5*t6i;	// S1
		Sr2 =      ss2*tar+ss4*t9r-ss5*t8r-ss3*t7r-ss1*t6r;		Si2 =      ss2*tai+ss4*t9i-ss5*t8i-ss3*t7i-ss1*t6i;	// S2
		Sr3 =      ss3*tar-ss5*t9r-ss2*t8r+ss1*t7r+ss4*t6r;		Si3 =      ss3*tai-ss5*t9i-ss2*t8i+ss1*t7i+ss4*t6i;	// S3
		Sr4 =      ss4*tar-ss3*t9r+ss1*t8r+ss5*t7r-ss2*t6r;		Si4 =      ss4*tai-ss3*t9i+ss1*t8i+ss5*t7i-ss2*t6i;	// S4
		Sr5 =      ss5*tar-ss1*t9r+ss4*t8r-ss2*t7r+ss3*t6r;		Si5 =      ss5*tai-ss1*t9i+ss4*t8i-ss2*t7i+ss3*t6i;	// S5

		In a 16-reg FMA model, can keep the 10 S-terms and all 5 shared sincos data in registers, but that
		requires us to reload one of the 2 repeated t*-mults (e.g. tai in the 1st block) 5 times per block.
		OTOH if we keep both t*r,i-mults and 4 of the 5 ss-mults in-reg, we just need 2 loads-from-mem per block:
	*/\
	"movq	%[__cc],%%r9	\n\t"/* cc1 */\
		"vmovaps	%%ymm1,%%ymm2		\n\t	vmovaps	%%ymm6 ,%%ymm7		\n\t"/* S2 = ta */\
		"vmovaps	%%ymm1,%%ymm3		\n\t	vmovaps	%%ymm6 ,%%ymm8		\n\t"/* S3 = ta */\
		"vmovaps	%%ymm1,%%ymm4		\n\t	vmovaps	%%ymm6 ,%%ymm9		\n\t"/* S4 = ta */\
		"vmovaps	%%ymm1,%%ymm5		\n\t	vmovaps	%%ymm6 ,%%ymm10		\n\t"/* S5 = ta */\
	"addq	$0xa0,%%r9	\n\t"/* Incr trig ptr to point to ss1 */\
		"vmovaps	0x20(%%r9 ),%%ymm11	\n\t	vmovaps	0x40(%%r9 ),%%ymm12		\n\t"/* ss2,ss3 */\
		"vmovaps	0x60(%%r9 ),%%ymm13	\n\t	vmovaps	0x80(%%r9 ),%%ymm14		\n\t"/* ss4,ss5 */\
\
	"	vmulpd	(%%r9),%%ymm1,%%ymm1	\n\t	vmulpd	(%%r9),%%ymm6 ,%%ymm6	\n\t"/* S1  = ss1*ta; */\
	"	vmulpd	%%ymm11,%%ymm2,%%ymm2	\n\t	vmulpd	%%ymm11,%%ymm7 ,%%ymm7	\n\t"/* S2  = ss2*ta; */\
	"	vmulpd	%%ymm12,%%ymm3,%%ymm3	\n\t	vmulpd	%%ymm12,%%ymm8 ,%%ymm8	\n\t"/* S3  = ss3*ta; */\
	"	vmulpd	%%ymm13,%%ymm4,%%ymm4	\n\t	vmulpd	%%ymm13,%%ymm9 ,%%ymm9	\n\t"/* S4  = ss4*ta; */\
	"	vmulpd	%%ymm14,%%ymm5,%%ymm5	\n\t	vmulpd	%%ymm14,%%ymm10,%%ymm10	\n\t"/* S5  = ss5*ta; */\
\
		"vmovaps	    (%%r14),%%ymm0	\n\t	vmovaps	0x20(%%r14),%%ymm15	\n\t"/* t9 */\
	" vfmadd231pd %%ymm11,%%ymm0,%%ymm1	\n\t  vfmadd231pd %%ymm11,%%ymm15,%%ymm6	\n\t"/* S1 += ss2*t9; */\
	" vfmadd231pd %%ymm13,%%ymm0,%%ymm2	\n\t  vfmadd231pd %%ymm13,%%ymm15,%%ymm7	\n\t"/* S2 += ss4*t9; */\
	"vfnmadd231pd %%ymm14,%%ymm0,%%ymm3	\n\t vfnmadd231pd %%ymm14,%%ymm15,%%ymm8	\n\t"/* S3 -= ss5*t9; */\
	"vfnmadd231pd %%ymm12,%%ymm0,%%ymm4	\n\t vfnmadd231pd %%ymm12,%%ymm15,%%ymm9	\n\t"/* S4 -= ss3*t9; */\
	"vfnmadd231pd (%%r9 ),%%ymm0,%%ymm5	\n\t vfnmadd231pd (%%r9 ),%%ymm15,%%ymm10	\n\t"/* S5 -= ss1*t9; */\
\
		"vmovaps	    (%%r13),%%ymm0	\n\t	vmovaps	0x20(%%r13),%%ymm15	\n\t"/* t8 */\
	" vfmadd231pd %%ymm12,%%ymm0,%%ymm1	\n\t  vfmadd231pd %%ymm12,%%ymm15,%%ymm6	\n\t"/* S1 += ss3*t8; */\
	"vfnmadd231pd %%ymm14,%%ymm0,%%ymm2	\n\t vfnmadd231pd %%ymm14,%%ymm15,%%ymm7	\n\t"/* S2 -= ss5*t8; */\
	"vfnmadd231pd %%ymm11,%%ymm0,%%ymm3	\n\t vfnmadd231pd %%ymm11,%%ymm15,%%ymm8	\n\t"/* S3 -= ss2*t8; */\
	" vfmadd231pd (%%r9 ),%%ymm0,%%ymm4	\n\t  vfmadd231pd (%%r9 ),%%ymm15,%%ymm9	\n\t"/* S4 += ss1*t8; */\
	" vfmadd231pd %%ymm13,%%ymm0,%%ymm5	\n\t  vfmadd231pd %%ymm13,%%ymm15,%%ymm10	\n\t"/* S5 += ss4*t8; */\
\
		"vmovaps	    (%%r12),%%ymm0	\n\t	vmovaps	0x20(%%r12),%%ymm15	\n\t"/* t7 */\
	" vfmadd231pd %%ymm13,%%ymm0,%%ymm1	\n\t  vfmadd231pd %%ymm13,%%ymm15,%%ymm6	\n\t"/* S1 += ss4*t7; */\
	"vfnmadd231pd %%ymm12,%%ymm0,%%ymm2	\n\t vfnmadd231pd %%ymm12,%%ymm15,%%ymm7	\n\t"/* S2 -= ss3*t7; */\
	" vfmadd231pd (%%r9 ),%%ymm0,%%ymm3	\n\t  vfmadd231pd (%%r9 ),%%ymm15,%%ymm8	\n\t"/* S3 += ss1*t7; */\
	" vfmadd231pd %%ymm14,%%ymm0,%%ymm4	\n\t  vfmadd231pd %%ymm14,%%ymm15,%%ymm9	\n\t"/* S4 += ss5*t7; */\
	"vfnmadd231pd %%ymm11,%%ymm0,%%ymm5	\n\t vfnmadd231pd %%ymm11,%%ymm15,%%ymm10	\n\t"/* S5 -= ss2*t7; */\
\
		"vmovaps	    (%%r11),%%ymm0	\n\t	vmovaps	0x20(%%r11),%%ymm15	\n\t"/* t6 */\
	" vfmadd231pd %%ymm14,%%ymm0,%%ymm1	\n\t  vfmadd231pd %%ymm14,%%ymm15,%%ymm6	\n\t"/* S1 += ss5*t6; */\
	"vfnmadd231pd (%%r9 ),%%ymm0,%%ymm2	\n\t vfnmadd231pd (%%r9 ),%%ymm15,%%ymm7	\n\t"/* S2 -= ss1*t6; */\
	" vfmadd231pd %%ymm13,%%ymm0,%%ymm3	\n\t  vfmadd231pd %%ymm13,%%ymm15,%%ymm8	\n\t"/* S3 += ss4*t6; */\
	"vfnmadd231pd %%ymm11,%%ymm0,%%ymm4	\n\t vfnmadd231pd %%ymm11,%%ymm15,%%ymm9	\n\t"/* S4 -= ss2*t6; */\
	" vfmadd231pd %%ymm12,%%ymm0,%%ymm5	\n\t  vfmadd231pd %%ymm12,%%ymm15,%%ymm10	\n\t"/* S5 += ss3*t6; */\
\
	/* Write sine-term-subconvo outputs back to memory to free up registers for the 5x5 cosine-term-subconvo computation: */\
	/* These i-ptrs no longer needed, so load o-ptrs in their place: */\
		"movslq	0x28(%%rdx),%%r15	\n\t	addq	%%rcx,%%r15		\n\t"/* &BA */\
		"movslq	0x24(%%rdx),%%r14	\n\t	addq	%%rcx,%%r14		\n\t"/* &B9 */\
		"movslq	0x20(%%rdx),%%r13	\n\t	addq	%%rcx,%%r13		\n\t"/* &B8 */\
		"movslq	0x1c(%%rdx),%%r12	\n\t	addq	%%rcx,%%r12		\n\t"/* &B7 */\
		"movslq	0x18(%%rdx),%%r11	\n\t	addq	%%rcx,%%r11		\n\t"/* &B6 */\
		"vmovaps	%%ymm1,(%%r15)		\n\t"/* s1r */"		vmovaps	%%ymm6 ,0x20(%%r15)	\n\t"/* s1i */\
		"vmovaps	%%ymm2,(%%r14)		\n\t"/* s2r */"		vmovaps	%%ymm7 ,0x20(%%r14)	\n\t"/* s2i */\
		"vmovaps	%%ymm3,(%%r13)		\n\t"/* s3r */"		vmovaps	%%ymm8 ,0x20(%%r13)	\n\t"/* s3i */\
		"vmovaps	%%ymm4,(%%r12)		\n\t"/* s4r */"		vmovaps	%%ymm9 ,0x20(%%r12)	\n\t"/* s4i */\
		"vmovaps	%%ymm5,(%%r11)		\n\t"/* s5r */"		vmovaps	%%ymm10,0x20(%%r11)	\n\t"/* s5i */\
	/*
		Cr1 = t0r+cc1*t1r+cc2*t2r+cc3*t3r+cc4*t4r+cc5*t5r;		Ci1 = t0i+cc1*t1i+cc2*t2i+cc3*t3i+cc4*t4i+cc5*t5i;	// C1
		Cr2 = t0r+cc2*t1r+cc4*t2r+cc5*t3r+cc3*t4r+cc1*t5r;		Ci2 = t0i+cc2*t1i+cc4*t2i+cc5*t3i+cc3*t4i+cc1*t5i;	// C2
		Cr3 = t0r+cc3*t1r+cc5*t2r+cc2*t3r+cc1*t4r+cc4*t5r;		Ci3 = t0i+cc3*t1i+cc5*t2i+cc2*t3i+cc1*t4i+cc4*t5i;	// C3
		Cr4 = t0r+cc4*t1r+cc3*t2r+cc1*t3r+cc5*t4r+cc2*t5r;		Ci4 = t0i+cc4*t1i+cc3*t2i+cc1*t3i+cc5*t4i+cc2*t5i;	// C4
		Cr5 = t0r+cc5*t1r+cc1*t2r+cc4*t3r+cc2*t4r+cc3*t5r;		Ci5 = t0i+cc5*t1i+cc1*t2i+cc4*t3i+cc2*t4i+cc3*t5i;	// C5
		B0r = t0r+    t1r+    t2r+    t3r+    t4r+    t5r;		B0i = t0i+    t1i+    t2i+    t3i+    t4i+    t5i;	// X0

		In a 16-reg FMA model, makes sense to init the c-terms to the t0-summand,
		and keep the 10 c-terms and 4 of the 5 shared sincos data in registers. That uses 14 regs,
		leaving 2 for the 2 t*[r,i] data shared by each such block. Those need to be in-reg (at least
		in the cosine-mults section) because of the need to accumulate the DC terms: E.g. at end of the
		first block below (after doing the 10 c-term-update FMAs) we add the current values of B0r,i
		(which we init to t0r,i in the ASM version) to t1r,i, then write the result to memory and
		load t2r,i into the same 2 regs in preparation for the next such block.
	*/\
	/*
	Current GPR assignments:
		rax	rbx	rcx	rdx	rsi	rdi	r8	r9	r10	r11	r12	r13	r14	r15
		in	iff	out	off	i1	i2	i3	ss1	i5	o6	o7	o8	o9	oA
	*/\
	"subq	$0xa0,%%r9	\n\t"/* Decr trig ptr to point to cc1 */\
	"movslq	(%%rbx),%%r8	\n\t	addq	%%rax,%%r8	\n\t"/* i0-from-offset */\
		"vmovaps	(%%r8 ),%%ymm1		\n\t	vmovaps	0x20(%%r8 ),%%ymm6	\n\t"/* c1 = A0 */\
		"vmovaps	%%ymm1,%%ymm2		\n\t	vmovaps	%%ymm6 ,%%ymm7		\n\t"/* c2 = A0 */\
		"vmovaps	%%ymm1,%%ymm3		\n\t	vmovaps	%%ymm6 ,%%ymm8		\n\t"/* c3 = A0 */\
		"vmovaps	%%ymm1,%%ymm4		\n\t	vmovaps	%%ymm6 ,%%ymm9		\n\t"/* c4 = A0 */\
		"vmovaps	%%ymm1,%%ymm5		\n\t	vmovaps	%%ymm6 ,%%ymm10		\n\t"/* c5 = A0 */\
\
		"vmovaps	0x20(%%r9 ),%%ymm11	\n\t	vmovaps	0x40(%%r9 ),%%ymm12		\n\t"/* cc2,cc3 */\
		"vmovaps	0x60(%%r9 ),%%ymm13	\n\t	vmovaps	0x80(%%r9 ),%%ymm14		\n\t"/* cc4,cc5 */\
\
		"vmovaps	    (%%rsi),%%ymm0	\n\t	vmovaps	0x20(%%rsi),%%ymm15	\n\t"/* t1 */\
	"vfmadd231pd (%%r9 ),%%ymm0,%%ymm1	\n\t vfmadd231pd (%%r9 ),%%ymm15,%%ymm6	\n\t"/* c1 += cc1*t1; */\
	"vfmadd231pd %%ymm11,%%ymm0,%%ymm2	\n\t vfmadd231pd %%ymm11,%%ymm15,%%ymm7	\n\t"/* c2 += cc2*t1; */\
	"vfmadd231pd %%ymm12,%%ymm0,%%ymm3	\n\t vfmadd231pd %%ymm12,%%ymm15,%%ymm8	\n\t"/* c3 += cc3*t1; */\
	"vfmadd231pd %%ymm13,%%ymm0,%%ymm4	\n\t vfmadd231pd %%ymm13,%%ymm15,%%ymm9	\n\t"/* c4 += cc4*t1; */\
	"vfmadd231pd %%ymm14,%%ymm0,%%ymm5	\n\t vfmadd231pd %%ymm14,%%ymm15,%%ymm10	\n\t"/* c5 += cc5*t1; */\
		"vaddpd (%%r8 ),%%ymm0,%%ymm0	\n\t	vaddpd 0x20(%%r8 ),%%ymm15,%%ymm15	\n\t"/* B0 += t1; */\
		"vmovaps	%%ymm0,(%%r8 )		\n\t	vmovaps	%%ymm15,0x20(%%r8 )			\n\t"/* Store B0 */\
\
		"vmovaps	    (%%rdi),%%ymm0	\n\t	vmovaps	0x20(%%rdi),%%ymm15	\n\t"/* t2 */\
		"movslq	0x0c(%%rbx),%%rsi	\n\t	addq	%%rax,%%rsi		\n\t"/* Compute &i3,4 from 32-bit offsets */\
		"movslq	0x10(%%rbx),%%rdi	\n\t	addq	%%rax,%%rdi		\n\t"\
	"vfmadd231pd %%ymm11,%%ymm0,%%ymm1	\n\t vfmadd231pd %%ymm11,%%ymm15,%%ymm6	\n\t"/* c1 += cc2*t2; */\
	"vfmadd231pd %%ymm13,%%ymm0,%%ymm2	\n\t vfmadd231pd %%ymm13,%%ymm15,%%ymm7	\n\t"/* c2 += cc4*t2; */\
	"vfmadd231pd %%ymm14,%%ymm0,%%ymm3	\n\t vfmadd231pd %%ymm14,%%ymm15,%%ymm8	\n\t"/* c3 += cc5*t2; */\
	"vfmadd231pd %%ymm12,%%ymm0,%%ymm4	\n\t vfmadd231pd %%ymm12,%%ymm15,%%ymm9	\n\t"/* c4 += cc3*t2; */\
	"vfmadd231pd (%%r9 ),%%ymm0,%%ymm5	\n\t vfmadd231pd (%%r9 ),%%ymm15,%%ymm10	\n\t"/* c5 += cc1*t2; */\
		"vaddpd (%%r8 ),%%ymm0,%%ymm0	\n\t	vaddpd 0x20(%%r8 ),%%ymm15,%%ymm15	\n\t"/* B0 += t2; */\
		"vmovaps	%%ymm0,(%%r8 )		\n\t	vmovaps	%%ymm15,0x20(%%r8 )			\n\t"/* Store B0 */\
\
		"vmovaps	    (%%rsi),%%ymm0	\n\t	vmovaps	0x20(%%rsi),%%ymm15	\n\t"/* t3 */\
	"vfmadd231pd %%ymm12,%%ymm0,%%ymm1	\n\t vfmadd231pd %%ymm12,%%ymm15,%%ymm6	\n\t"/* c1 += cc3*t3; */\
	"vfmadd231pd %%ymm14,%%ymm0,%%ymm2	\n\t vfmadd231pd %%ymm14,%%ymm15,%%ymm7	\n\t"/* c2 += cc5*t3; */\
	"vfmadd231pd %%ymm11,%%ymm0,%%ymm3	\n\t vfmadd231pd %%ymm11,%%ymm15,%%ymm8	\n\t"/* c3 += cc2*t3; */\
	"vfmadd231pd (%%r9 ),%%ymm0,%%ymm4	\n\t vfmadd231pd (%%r9 ),%%ymm15,%%ymm9	\n\t"/* c4 += cc1*t3; */\
	"vfmadd231pd %%ymm13,%%ymm0,%%ymm5	\n\t vfmadd231pd %%ymm13,%%ymm15,%%ymm10	\n\t"/* c5 += cc4*t3; */\
		"vaddpd (%%r8 ),%%ymm0,%%ymm0	\n\t	vaddpd 0x20(%%r8 ),%%ymm15,%%ymm15	\n\t"/* B0 += t3; */\
		"vmovaps	%%ymm0,(%%r8 )		\n\t	vmovaps	%%ymm15,0x20(%%r8 )			\n\t"/* Store B0 */\
\
		"vmovaps	    (%%rdi),%%ymm0	\n\t	vmovaps	0x20(%%rdi),%%ymm15	\n\t"/* t4 */\
	"vfmadd231pd %%ymm13,%%ymm0,%%ymm1	\n\t vfmadd231pd %%ymm13,%%ymm15,%%ymm6	\n\t"/* c1 += cc4*t4; */\
	"vfmadd231pd %%ymm12,%%ymm0,%%ymm2	\n\t vfmadd231pd %%ymm12,%%ymm15,%%ymm7	\n\t"/* c2 += cc3*t4; */\
	"vfmadd231pd (%%r9 ),%%ymm0,%%ymm3	\n\t vfmadd231pd (%%r9 ),%%ymm15,%%ymm8	\n\t"/* c3 += cc1*t4; */\
	"vfmadd231pd %%ymm14,%%ymm0,%%ymm4	\n\t vfmadd231pd %%ymm14,%%ymm15,%%ymm9	\n\t"/* c4 += cc5*t4; */\
	"vfmadd231pd %%ymm11,%%ymm0,%%ymm5	\n\t vfmadd231pd %%ymm11,%%ymm15,%%ymm10	\n\t"/* c5 += cc2*t4; */\
		"vaddpd (%%r8 ),%%ymm0,%%ymm0	\n\t	vaddpd 0x20(%%r8 ),%%ymm15,%%ymm15	\n\t"/* B0 += t4; */\
		"vmovaps	%%ymm0,(%%r8 )		\n\t	vmovaps	%%ymm15,0x20(%%r8 )			\n\t"/* Store B0 */\
\
		"vmovaps	    (%%r10),%%ymm0	\n\t	vmovaps	0x20(%%r10),%%ymm15	\n\t"/* t5 */\
	"vfmadd231pd %%ymm14,%%ymm0,%%ymm1	\n\t vfmadd231pd %%ymm14,%%ymm15,%%ymm6	\n\t"/* c1 += cc5*t5; */\
	"vfmadd231pd (%%r9 ),%%ymm0,%%ymm2	\n\t vfmadd231pd (%%r9 ),%%ymm15,%%ymm7	\n\t"/* c2 += cc1*t5; */\
	"vfmadd231pd %%ymm13,%%ymm0,%%ymm3	\n\t vfmadd231pd %%ymm13,%%ymm15,%%ymm8	\n\t"/* c3 += cc4*t5; */\
	"vfmadd231pd %%ymm11,%%ymm0,%%ymm4	\n\t vfmadd231pd %%ymm11,%%ymm15,%%ymm9	\n\t"/* c4 += cc2*t5; */\
	"vfmadd231pd %%ymm12,%%ymm0,%%ymm5	\n\t vfmadd231pd %%ymm12,%%ymm15,%%ymm10	\n\t"/* c5 += cc3*t5; */\
		"vaddpd (%%r8 ),%%ymm0,%%ymm0	\n\t	vaddpd 0x20(%%r8 ),%%ymm15,%%ymm15	\n\t"/* B0 += t5; */\
	/*
	Current GPR assignments:
		rax	rbx	rcx	rdx	rsi	rdi	r8	r9	r10	r11	r12	r13	r14	r15
		in	iff	out	off	i3	i4	i0	cc1	i5	o6	o7	o8	o9	oA
	*/\
		"movslq	(%%rdx),%%r8	\n\t	addq	%%rcx,%%r8	\n\t"/* o0-from-offset */\
		"vmovaps	%%ymm0,(%%r8 )		\n\t	vmovaps	%%ymm15,0x20(%%r8 )			\n\t"/* Store B0, now into output slot */\
		/* Have enough registers to also hold the real parts of S-terms: */\
		"vmovaps	(%%r15),%%ymm11		\n\t"/* s1r */\
		"vmovaps	(%%r14),%%ymm12		\n\t"/* s2r */\
		"vmovaps	(%%r13),%%ymm13		\n\t"/* s3r */\
		"vmovaps	(%%r12),%%ymm14		\n\t"/* s4r */\
		"vmovaps	(%%r11),%%ymm15		\n\t"/* s5r */\
		/* Need rcx,rdx for O-address computations, but done with I-addesses so rax,rbx free: */\
		"movslq	0x04(%%rdx),%%rsi	\n\t	addq	%%rcx,%%rsi		\n\t"/* &B1 */\
		"movslq	0x08(%%rdx),%%rdi	\n\t	addq	%%rcx,%%rdi		\n\t"/* &B2 */\
		"movslq	0x0c(%%rdx),%%rax	\n\t	addq	%%rcx,%%rax		\n\t"/* &B3 */\
		"movslq	0x10(%%rdx),%%rbx	\n\t	addq	%%rcx,%%rbx		\n\t"/* &B4 */\
		"movslq	0x14(%%rdx),%%r10	\n\t	addq	%%rcx,%%r10		\n\t"/* &B5 */\
	/* ...Do the lower-half-index portion of the closing wave of radix-2 butterflies: */\
		"vsubpd 0x20(%%r15),%%ymm1,%%ymm1		\n\t	 vaddpd %%ymm11,%%ymm6 ,%%ymm6	\n\t"/* B1 = C1 + I*S1 */\
		"vsubpd 0x20(%%r14),%%ymm2,%%ymm2		\n\t	 vaddpd %%ymm12,%%ymm7 ,%%ymm7	\n\t"/* B2 = C2 + I*S2 */\
		"vsubpd 0x20(%%r13),%%ymm3,%%ymm3		\n\t	 vaddpd %%ymm13,%%ymm8 ,%%ymm8	\n\t"/* B3 = C3 + I*S3 */\
		"vsubpd 0x20(%%r12),%%ymm4,%%ymm4		\n\t	 vaddpd %%ymm14,%%ymm9 ,%%ymm9	\n\t"/* B4 = C4 + I*S4 */\
		"vsubpd 0x20(%%r11),%%ymm5,%%ymm5		\n\t	 vaddpd %%ymm15,%%ymm10,%%ymm10	\n\t"/* B5 = C5 + I*S5 */\
	"vmovaps	-0x40(%%r9),%%ymm0		\n\t"/* 2.0 */\
	/* Write lower-half outputs back to memory... */\
		"vmovaps	%%ymm1,(%%rsi)		\n\t"/* B1r */"		vmovaps	%%ymm6 ,0x20(%%rsi)	\n\t"/* B1i */\
		"vmovaps	%%ymm2,(%%rdi)		\n\t"/* B2r */"		vmovaps	%%ymm7 ,0x20(%%rdi)	\n\t"/* B2i */\
		"vmovaps	%%ymm3,(%%rax)		\n\t"/* B3r */"		vmovaps	%%ymm8 ,0x20(%%rax)	\n\t"/* B3i */\
		"vmovaps	%%ymm4,(%%rbx)		\n\t"/* B4r */"		vmovaps	%%ymm9 ,0x20(%%rbx)	\n\t"/* B4i */\
		"vmovaps	%%ymm5,(%%r10)		\n\t"/* B5r */"		vmovaps	%%ymm10,0x20(%%r10)	\n\t"/* B5i */\
	/* ...Do the upper-half-index portion of the radix-2 butterflies: */\
		" vfmadd231pd 0x20(%%r15),%%ymm0,%%ymm1		\n\t	vfnmadd231pd %%ymm11,%%ymm0,%%ymm6	\n\t"/* B10= C1 - I*S1 */\
		" vfmadd231pd 0x20(%%r14),%%ymm0,%%ymm2		\n\t	vfnmadd231pd %%ymm12,%%ymm0,%%ymm7	\n\t"/* B9 = C2 - I*S2 */\
		" vfmadd231pd 0x20(%%r13),%%ymm0,%%ymm3		\n\t	vfnmadd231pd %%ymm13,%%ymm0,%%ymm8	\n\t"/* B8 = C3 - I*S3 */\
		" vfmadd231pd 0x20(%%r12),%%ymm0,%%ymm4		\n\t	vfnmadd231pd %%ymm14,%%ymm0,%%ymm9	\n\t"/* B7 = C4 - I*S4 */\
		" vfmadd231pd 0x20(%%r11),%%ymm0,%%ymm5		\n\t	vfnmadd231pd %%ymm15,%%ymm0,%%ymm10	\n\t"/* B6 = C5 - I*S5 */\
	/* ...And write the upper-half outputs back to memory. */\
		"vmovaps	%%ymm1,(%%r15)		\n\t"/* Bar */"		vmovaps	%%ymm6 ,0x20(%%r15)	\n\t"/* Bai */\
		"vmovaps	%%ymm2,(%%r14)		\n\t"/* B9r */"		vmovaps	%%ymm7 ,0x20(%%r14)	\n\t"/* B9i */\
		"vmovaps	%%ymm3,(%%r13)		\n\t"/* B8r */"		vmovaps	%%ymm8 ,0x20(%%r13)	\n\t"/* B8i */\
		"vmovaps	%%ymm4,(%%r12)		\n\t"/* B7r */"		vmovaps	%%ymm9 ,0x20(%%r12)	\n\t"/* B7i */\
		"vmovaps	%%ymm5,(%%r11)		\n\t"/* B6r */"		vmovaps	%%ymm10,0x20(%%r11)	\n\t"/* B6i */\
		:					/* outputs: none */\
		: [in0] "m" (Xin0)	/* All inputs from memory addresses here */\
		 ,[i_off] "m" (Xi_off)	/* I-address pointer-strides */\
		 ,[__cc] "m" (Xcc)\
		 ,[out0] "m" (Xout0)\
		 ,[o_off] "m" (Xo_off)	/* O-address pointer-strides */\
		: "cc","memory","rax","rbx","rcx","rdx","rsi","rdi","r8","r9","r10","r11","r12","r13","r14","r15",\
			"xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
		);\
	}

   #else		// [2] [LOACC, no-flags default] based on Radix-11 DFT impl based on RADIX_11_DFT in dft_macro.h,
				// using length-5 cyclic convolution scheme for the 5 x 5 matrix submultiplies.
				// Totals: [134 ADD, 54 MUL/FMA, 172 memref] in this FMA-lite version.
				//
	#define SSE2_RADIX_11_DFT(XI0,Xi1,Xi2,Xi3,Xi4,Xi5,Xi6,Xi7,Xi8,Xi9,XiA, Xcc, XO0,Xo1,Xo2,Xo3,Xo4,Xo5,Xo6,Xo7,Xo8,Xo9,XoA)\
	{\
	__asm__ volatile (\
		/********************************************/\
		/*       Here are the 5 cosine terms:       */\
		/********************************************/\
		"movq	%[__I0],%%rax			\n\t			movq	%[__O0],%%rcx		\n\t"/* rax/rcx point to Re parts of I/Os */\
		"leaq	0x20(%%rax),%%rbx		\n\t			leaq	0x20(%%rcx),%%rdx	\n\t"/* rbx/rdx point to Im parts of I/Os */\
		"movq	%[__cc],%%rsi			\n\t"/* cc1 ... two in -0x40(%%rsi) */\
		/*****************************/				/******************************/\
		/*        Real Parts:        */				/*     Imaginary Parts:       */\
		/*****************************/				/******************************/\
	"movq	%[__i1],%%r8	\n\t	vmovaps (%%r8 ),%%ymm1	\n\t	vmovaps	0x20(%%r8 ),%%ymm9	\n\t"\
	"movq	%[__iA],%%r9	\n\t	vmovaps (%%r9 ),%%ymm5	\n\t	vmovaps	0x20(%%r9 ),%%ymm13	\n\t"\
	"movq	%[__i2],%%r10	\n\t	vmovaps (%%r10),%%ymm2	\n\t	vmovaps	0x20(%%r10),%%ymm10	\n\t"\
	"movq	%[__i9],%%r11	\n\t	vmovaps (%%r11),%%ymm6	\n\t	vmovaps	0x20(%%r11),%%ymm14	\n\t"\
	"movq	%[__i3],%%r12	\n\t	vmovaps (%%r12),%%ymm3	\n\t	vmovaps	0x20(%%r12),%%ymm11	\n\t"\
	"movq	%[__i8],%%r13	\n\t	vmovaps (%%r13),%%ymm7	\n\t	vmovaps	0x20(%%r13),%%ymm15	\n\t"\
	"movq	%[__i4],%%r14	\n\t	vmovaps (%%r14),%%ymm4	\n\t	vmovaps	0x20(%%r14),%%ymm12	\n\t"\
	"movq	%[__i7],%%r15	\n\t	vmovaps (%%r15),%%ymm0	\n\t	vmovaps	0x20(%%r15),%%ymm8	\n\t"\
		"vsubpd	%%ymm5,%%ymm1,%%ymm1	\n\t			vsubpd	%%ymm13,%%ymm9 ,%%ymm9			\n\t"/* t21 = x1-xA */\
		"vsubpd	%%ymm6,%%ymm2,%%ymm2	\n\t			vsubpd	%%ymm14,%%ymm10,%%ymm10			\n\t"/* t19 = x2-x9 */\
		"vsubpd	%%ymm7,%%ymm3,%%ymm3	\n\t			vsubpd	%%ymm15,%%ymm11,%%ymm11			\n\t"/* t17 = x3-x8 */\
		"vsubpd	%%ymm0,%%ymm4,%%ymm4	\n\t			vsubpd	%%ymm8 ,%%ymm12,%%ymm12			\n\t"/* t15 = x4-x7 */\
	"movq	%[__oA],%%r8	\n\t	vmovaps %%ymm1,0x20(%%r8 )	\n\t	vmovaps %%ymm9 ,(%%r8 )	\n\t"/* ...and store the above 4 pairs */\
	"movq	%[__o9],%%r9	\n\t	vmovaps %%ymm2,0x20(%%r9 )	\n\t	vmovaps %%ymm10,(%%r9 )	\n\t"/* of re/im-part diffs back to mem */\
	"movq	%[__o8],%%r10	\n\t	vmovaps %%ymm3,0x20(%%r10)	\n\t	vmovaps %%ymm11,(%%r10)	\n\t"/* slots oA,8,8,7, swapping the re/im. */\
	"movq	%[__o7],%%r11	\n\t	vmovaps %%ymm4,0x20(%%r11)	\n\t	vmovaps %%ymm12,(%%r11)	\n\t"\
	"vfmadd231pd -0x40(%%rsi),%%ymm5,%%ymm1	\n\t	vfmadd231pd -0x40(%%rsi),%%ymm13,%%ymm9	\n\t"/* t3  = x1+xA */\
	"vfmadd231pd -0x40(%%rsi),%%ymm6,%%ymm2	\n\t	vfmadd231pd -0x40(%%rsi),%%ymm14,%%ymm10	\n\t"/* t5  = x2+x9 */\
	"movq	%[__i5],%%r9	\n\t	vmovaps (%%r9 ),%%ymm5	\n\t	vmovaps 0x20(%%r9 ),%%ymm13	\n\t"/* x5 */\
	"movq	%[__i6],%%r11	\n\t	vmovaps (%%r11),%%ymm6	\n\t	vmovaps 0x20(%%r11),%%ymm14	\n\t"/* x6 */\
	"vfmadd231pd -0x40(%%rsi),%%ymm7,%%ymm3	\n\t	vfmadd231pd -0x40(%%rsi),%%ymm15,%%ymm11	\n\t"/* t7  = x3+x8 */\
	"vfmadd231pd -0x40(%%rsi),%%ymm0,%%ymm4	\n\t	vfmadd231pd -0x40(%%rsi),%%ymm8 ,%%ymm12	\n\t"/* t9  = x4+x7 */\
		"vsubpd	%%ymm6,%%ymm5,%%ymm5	\n\t			vsubpd	%%ymm14,%%ymm13,%%ymm13			\n\t"/* t13 = x5-x6 */\
	"movq	%[__o6],%%r8	\n\t	vmovaps %%ymm5,0x20(%%r8 )	\n\t	vmovaps %%ymm13,(%%r8 )	\n\t"/* store t13 = (x5-x6) in o6 */\
		"vmovaps	(%%rax),%%ymm0		\n\t			vmovaps	(%%rbx),%%ymm8					\n\t"/* x0 */\
	"vfmadd231pd -0x40(%%rsi),%%ymm6,%%ymm5	\n\t	vfmadd231pd -0x40(%%rsi),%%ymm14,%%ymm13	\n\t"/* t11 = x5+x6 */\
	"movq	%[__o1],%%r11	\n\t"\
	"movq	%[__o2],%%r12	\n\t"\
	"movq	%[__o3],%%r13	\n\t"\
	"movq	%[__o4],%%r14	\n\t"\
	"movq	%[__o5],%%r15	\n\t"\
		"vsubpd	%%ymm2,%%ymm1,%%ymm1			\n\t	vsubpd	%%ymm10,%%ymm9 ,%%ymm9			\n\t"/* c1 = t3 -t5 */\
		"vsubpd	%%ymm2,%%ymm5,%%ymm5			\n\t	vsubpd	%%ymm10,%%ymm13,%%ymm13			\n\t"/* c2 = t11-t5 */\
		"vmovaps	%%ymm1,%%ymm6				\n\t	vmovaps	%%ymm9 ,%%ymm14					\n\t"/* cpy c1 */\
		"vsubpd	%%ymm2,%%ymm3,%%ymm3			\n\t	vsubpd	%%ymm10,%%ymm11,%%ymm11			\n\t"/* c4 = t7 -t5 */\
		"vmovaps	%%ymm1,%%ymm7				\n\t	vmovaps	%%ymm9 ,%%ymm15					\n\t"/* cpy c1 */\
	/*	"vmulpd	      (%%rsi),%%ymm1,%%ymm1		\n\t	vmulpd	      (%%rsi),%%ymm9 ,%%ymm9	\n\t"// a0.c1 */\
		"vmovaps	%%ymm1,(%%r11)				\n\t	vmovaps	%%ymm9 ,0x20(%%r11)				\n\t"/* store copy of a0.c1 in mem */\
		"vsubpd	%%ymm2,%%ymm4,%%ymm4			\n\t	vsubpd	%%ymm10,%%ymm12,%%ymm12			\n\t"/* c5 = t9-t5 */\
		"vaddpd	%%ymm5,%%ymm6,%%ymm6			\n\t	vaddpd	%%ymm13,%%ymm14,%%ymm14			\n\t"/* c3 = c1+c2 */\
		"vsubpd	%%ymm3,%%ymm7,%%ymm7			\n\t	vsubpd	%%ymm11,%%ymm15,%%ymm15			\n\t"/* c7 = c1-c4 */\
	/*	"vmulpd	 0x0c0(%%rsi),%%ymm7,%%ymm7		\n\t	vmulpd	 0x0c0(%%rsi),%%ymm15,%%ymm15	\n\t"// a6.c7 */\
		"vmovaps	%%ymm7,(%%r13)				\n\t	vmovaps	%%ymm15,0x20(%%r13)				\n\t"/* store copy of a6.c7 in mem */\
		"vmovaps	%%ymm3,%%ymm1				\n\t	vmovaps	%%ymm11,%%ymm9					\n\t"/* cpy c4 */\
		"vaddpd	%%ymm4,%%ymm3,%%ymm3			\n\t	vaddpd	%%ymm12,%%ymm11,%%ymm11			\n\t"/* c6 = c4+c5 */\
	/*	"vmulpd	 0x060(%%rsi),%%ymm1,%%ymm1		\n\t	vmulpd	 0x060(%%rsi),%%ymm9 ,%%ymm9	\n\t"// a3.c4 */\
		"vmovaps	%%ymm1,(%%r12)				\n\t	vmovaps	%%ymm9 ,0x20(%%r12)				\n\t"/* store copy of a3.c4 in mem */\
		"vmovaps	%%ymm5,%%ymm7				\n\t	vmovaps	%%ymm13,%%ymm15					\n\t"/* cpy c2 */\
		"vsubpd	%%ymm4,%%ymm5,%%ymm5			\n\t	vsubpd	%%ymm12,%%ymm13,%%ymm13			\n\t"/* c8 = c2-c5 */\
	"vfmadd132pd -0x20(%%rsi),%%ymm3,%%ymm2		\n\t vfmadd132pd -0x20(%%rsi),%%ymm11,%%ymm10	\n\t"/* 5.t5+c6 */\
		"vmovaps	%%ymm6,%%ymm1				\n\t	vmovaps	%%ymm14,%%ymm9					\n\t"/* cpy c3 */\
		"vsubpd	%%ymm3,%%ymm6,%%ymm6			\n\t	vsubpd	%%ymm11,%%ymm14,%%ymm14			\n\t"/* c9 = c3-c6 */\
		"vmulpd	 0x100(%%rsi),%%ymm6,%%ymm6		\n\t	vmulpd	 0x100(%%rsi),%%ymm14,%%ymm14	\n\t"/* c9 = a8.c9 */\
		"vmulpd	 0x0a0(%%rsi),%%ymm3,%%ymm3		\n\t	vmulpd	 0x0a0(%%rsi),%%ymm11,%%ymm11	\n\t"/* c6 = a5.c6 */\
		"vaddpd	%%ymm1,%%ymm2,%%ymm2			\n\t	vaddpd	%%ymm9 ,%%ymm10,%%ymm10			\n\t"/* c10= c3 +c6+5*t5 */\
		"vmulpd	 0x040(%%rsi),%%ymm1,%%ymm1		\n\t	vmulpd	 0x040(%%rsi),%%ymm9 ,%%ymm9	\n\t"/* c3 = a2.c3 */\
		"vaddpd	%%ymm2,%%ymm0,%%ymm0			\n\t	vaddpd	%%ymm10,%%ymm8 ,%%ymm8			\n\t"/* B0 = A0+c10 */\
		"vmovaps	%%ymm0,(%%rcx)				\n\t	vmovaps	%%ymm8 ,(%%rdx)					\n\t"/* store B0 */\
	"vfmadd132pd 0x020(%%rsi),%%ymm1,%%ymm7		\n\t vfmadd132pd 0x020(%%rsi),%%ymm9 ,%%ymm15	\n\t"/* c2 = a1.c2+c3 */\
	"vfmadd132pd 0x080(%%rsi),%%ymm3,%%ymm4		\n\t vfmadd132pd 0x080(%%rsi),%%ymm11,%%ymm12	\n\t"/* c5 = a4.c5+c6 */\
	"vfmadd132pd 0x0e0(%%rsi),%%ymm6,%%ymm5		\n\t vfmadd132pd 0x0e0(%%rsi),%%ymm14,%%ymm13	\n\t"/* c8 = a7.c8+c9 */\
	"vfmadd132pd 0x120(%%rsi),%%ymm0,%%ymm2		\n\t vfmadd132pd 0x120(%%rsi),%%ymm8 ,%%ymm10	\n\t"/* c10= a9.c10+B0 */\
/* ymm0,8 free - use to reload c1,4,7 prior to FMA: */\
		"vmovaps 	(%%r11),%%ymm0				\n\t	vmovaps	0x20(%%r11),%%ymm8	\n\t"\
	"vfmadd231pd      (%%rsi),%%ymm0,%%ymm1		\n\t vfmadd231pd      (%%r11),%%ymm8 ,%%ymm9 	\n\t"/* c1 = a0.c1+c3 */\
		"vmovaps 	(%%r12),%%ymm0				\n\t	vmovaps	0x20(%%r12),%%ymm8	\n\t"\
	"vfmadd231pd 0x060(%%rsi),%%ymm0,%%ymm3		\n\t vfmadd231pd 0x060(%%r12),%%ymm8 ,%%ymm11	\n\t"/* c4 = a3.c4+c6 */\
		"vmovaps	 (%%r13),%%ymm0				\n\t	vmovaps	0x20(%%r13),%%ymm8	\n\t"\
	"vfmadd231pd 0x0c0(%%rsi),%%ymm0,%%ymm6		\n\t vfmadd231pd 0x0c0(%%r13),%%ymm8 ,%%ymm14	\n\t"/* c7 = a6.c7+c9 */\
		"vmovaps	%%ymm2,%%ymm0				\n\t	vmovaps	%%ymm10,%%ymm8					\n\t"/* cpy c10 */\
		"vsubpd	%%ymm1,%%ymm2,%%ymm2			\n\t	vsubpd	%%ymm9 ,%%ymm10,%%ymm10			\n\t"/* c10-c1 */\
		"vaddpd	%%ymm0,%%ymm1,%%ymm1			\n\t	vaddpd	%%ymm8 ,%%ymm9 ,%%ymm9			\n\t"/* c10+c1 */\
		"vsubpd	%%ymm7,%%ymm2,%%ymm2			\n\t	vsubpd	%%ymm15,%%ymm10,%%ymm10			\n\t"/* c10-c1-c2 */\
		"vaddpd	%%ymm0,%%ymm7,%%ymm7			\n\t	vaddpd	%%ymm8 ,%%ymm15,%%ymm15			\n\t"/* c10+c2 */\
		"vsubpd	%%ymm6,%%ymm1,%%ymm1			\n\t	vsubpd	%%ymm14,%%ymm9 ,%%ymm9			\n\t"/* cr1 = c10+c1-c7 */\
		"vsubpd	%%ymm3,%%ymm2,%%ymm2			\n\t	vsubpd	%%ymm11,%%ymm10,%%ymm10			\n\t"/* c10-c1-c2-c4 */\
		"vaddpd	%%ymm0,%%ymm3,%%ymm3			\n\t	vaddpd	%%ymm8 ,%%ymm11,%%ymm11			\n\t"/* c10+c4 */\
		"vsubpd	%%ymm5,%%ymm7,%%ymm7			\n\t	vsubpd	%%ymm13,%%ymm15,%%ymm15			\n\t"/* cr5 = c10+c2-c8 */\
		"vsubpd	%%ymm4,%%ymm2,%%ymm2			\n\t	vsubpd	%%ymm12,%%ymm10,%%ymm10			\n\t"/* cr2 = c10-c1-c2-c4-c5 */\
		"vaddpd	%%ymm0,%%ymm4,%%ymm4			\n\t	vaddpd	%%ymm8 ,%%ymm12,%%ymm12			\n\t"/* c10+c5 */\
		"vaddpd	%%ymm6,%%ymm3,%%ymm3			\n\t	vaddpd	%%ymm14,%%ymm11,%%ymm11			\n\t"/* cr3 = c10+c4+c7 */\
		"vaddpd	%%ymm5,%%ymm4,%%ymm4			\n\t	vaddpd	%%ymm13,%%ymm12,%%ymm12			\n\t"/* cr4 = c10+c5+c8 */\
		"vmovaps	%%ymm1,(%%r11)				\n\t	vmovaps	%%ymm9 ,0x20(%%r11)				\n\t"/* write cr1 */\
		"vmovaps	%%ymm7,(%%r15)				\n\t	vmovaps	%%ymm15,0x20(%%r15)				\n\t"/* write cr5 */\
		"vmovaps	%%ymm2,(%%r12)				\n\t	vmovaps	%%ymm10,0x20(%%r12)				\n\t"/* write cr2 */\
		"vmovaps	%%ymm3,(%%r13)				\n\t	vmovaps	%%ymm11,0x20(%%r13)				\n\t"/* write cr3 */\
		"vmovaps	%%ymm4,(%%r14)				\n\t	vmovaps	%%ymm12,0x20(%%r14)				\n\t"/* write cr4 */\
	/***********************************************************************************************************/\
	/* Here are 5 sine terms: Similar to cosines, but t21,19,17,15,13 replace t3,5,7,9,11 and some sign flips: */\
	/***********************************************************************************************************/\
		"addq	$0x140,%%rsi				\n\t"/* b0 */\
		/*****************************/				/******************************/\
		/*     Imaginary Parts:      */				/*        Real Parts:         *//*** Note re/im swap ***/\
		/*****************************/				/******************************/\
	"movq	%[__oA],%%r11	\n\t	vmovaps (%%r11),%%ymm1	\n\t	vmovaps 0x20(%%r11),%%ymm9	\n\t"/* t21 */\
	"movq	%[__o9],%%r12	\n\t	vmovaps (%%r12),%%ymm2	\n\t	vmovaps 0x20(%%r12),%%ymm10	\n\t"/* t19 */\
	"movq	%[__o8],%%r13	\n\t	vmovaps (%%r13),%%ymm3	\n\t	vmovaps 0x20(%%r13),%%ymm11	\n\t"/* t17 */\
	"movq	%[__o7],%%r14	\n\t	vmovaps (%%r14),%%ymm4	\n\t	vmovaps 0x20(%%r14),%%ymm12	\n\t"/* t15 */\
	"movq	%[__o6],%%r15	\n\t	vmovaps (%%r15),%%ymm5	\n\t	vmovaps 0x20(%%r15),%%ymm13	\n\t"/* t13 */\
		"vaddpd	%%ymm2,%%ymm1,%%ymm1			\n\t	vaddpd	%%ymm10,%%ymm9 ,%%ymm9			\n\t"/* c1 = t21+t19 */\
		"vaddpd	%%ymm2,%%ymm5,%%ymm5			\n\t	vaddpd	%%ymm10,%%ymm13,%%ymm13			\n\t"/* c2 = t13+t19 */\
		"vmovaps	%%ymm1,%%ymm6				\n\t	vmovaps	%%ymm9 ,%%ymm14					\n\t"/* cpy c1 */\
		"vaddpd	%%ymm2,%%ymm3,%%ymm3			\n\t	vaddpd	%%ymm10,%%ymm11,%%ymm11			\n\t"/* c4 = t17+t19 */\
		"vmovaps	%%ymm1,%%ymm7				\n\t	vmovaps	%%ymm9 ,%%ymm15					\n\t"/* cpy c1 */\
		"vmulpd	      (%%rsi),%%ymm1,%%ymm1		\n\t	vmulpd	      (%%rsi),%%ymm9 ,%%ymm9	\n\t"/* b0.c1 */\
		"vmovaps	%%ymm1,(%%r15)				\n\t	vmovaps	%%ymm9 ,0x20(%%r15)				\n\t"/* store b0.c1 */\
		"vaddpd	%%ymm2,%%ymm4,%%ymm4			\n\t	vaddpd	%%ymm10,%%ymm12,%%ymm12			\n\t"/* c5 = t15+t19 */\
		"vaddpd	%%ymm5,%%ymm6,%%ymm6			\n\t	vaddpd	%%ymm13,%%ymm14,%%ymm14			\n\t"/* c3 = c1+c2 */\
		"vsubpd	%%ymm3,%%ymm7,%%ymm7			\n\t	vsubpd	%%ymm11,%%ymm15,%%ymm15			\n\t"/* c7 = c1-c4 */\
		"vmulpd	 0x0c0(%%rsi),%%ymm7,%%ymm7		\n\t	vmulpd	 0x0c0(%%rsi),%%ymm15,%%ymm15	\n\t"/* b6.c7 */\
		"vmovaps	%%ymm7,(%%r13)				\n\t	vmovaps	%%ymm15,0x20(%%r13)				\n\t"/* store b6.c7 */\
		"vmovaps	%%ymm3,%%ymm1				\n\t	vmovaps	%%ymm11,%%ymm9					\n\t"/* cpy c4 */\
		"vaddpd	%%ymm4,%%ymm3,%%ymm3			\n\t	vaddpd	%%ymm12,%%ymm11,%%ymm11			\n\t"/* c6 = c4+c5 */\
		"vmulpd	 0x060(%%rsi),%%ymm1,%%ymm1		\n\t	vmulpd	 0x060(%%rsi),%%ymm9 ,%%ymm9	\n\t"/* b3.c4 */\
		"vmovaps	%%ymm1,(%%r14)				\n\t	vmovaps	%%ymm9 ,0x20(%%r14)				\n\t"/* store b3.c4 */\
		"vmovaps	%%ymm5,%%ymm7				\n\t	vmovaps	%%ymm13,%%ymm15					\n\t"/* cpy c2 */\
		"vsubpd	%%ymm4,%%ymm5,%%ymm5			\n\t	vsubpd	%%ymm12,%%ymm13,%%ymm13			\n\t"/* c8 = c2-c5 */\
	"vfmsub132pd -0x160(%%rsi),%%ymm3,%%ymm2	\n\t vfmsub132pd -0x160(%%rsi),%%ymm11,%%ymm10	\n\t"/* 5.t19-c6 */\
		"vmovaps	%%ymm6,%%ymm1				\n\t	vmovaps	%%ymm14,%%ymm9					\n\t"/* cpy c3 */\
		"vsubpd	%%ymm3,%%ymm6,%%ymm6			\n\t	vsubpd	%%ymm11,%%ymm14,%%ymm14			\n\t"/* c9 = c3-c6 */\
		"vmulpd	 0x100(%%rsi),%%ymm6,%%ymm6		\n\t	vmulpd	 0x100(%%rsi),%%ymm14,%%ymm14	\n\t"/* c9 = b8.c9 */\
		"vmulpd	 0x0a0(%%rsi),%%ymm3,%%ymm3		\n\t	vmulpd	 0x0a0(%%rsi),%%ymm11,%%ymm11	\n\t"/* c6 = b5.c6 */\
		"vsubpd	%%ymm1,%%ymm2,%%ymm2			\n\t	vsubpd	%%ymm9 ,%%ymm10,%%ymm10			\n\t"/* -c10 = 5.t19-c6-c3 (c6 prior to above mul-by-b5) */\
		"vmulpd	 0x040(%%rsi),%%ymm1,%%ymm1		\n\t	vmulpd	 0x040(%%rsi),%%ymm9 ,%%ymm9	\n\t"/* c3 = b2.c3 */\
		"vmulpd	 0x120(%%rsi),%%ymm2,%%ymm2		\n\t	vmulpd	 0x120(%%rsi),%%ymm10,%%ymm10	\n\t"/* c10 = -b9.-c10 */\
	"vfmadd132pd 0x0e0(%%rsi),%%ymm6,%%ymm5		\n\t vfmadd132pd 0x0e0(%%rsi),%%ymm14,%%ymm13	\n\t"/* c6 = b7.c8+c9 */\
		"vaddpd	(%%r13),%%ymm6,%%ymm6			\n\t	vaddpd	0x20(%%r13),%%ymm14,%%ymm14		\n\t"/* c5 = b6.c7+c9 */\
	"vfmadd132pd 0x080(%%rsi),%%ymm3,%%ymm4		\n\t vfmadd132pd 0x080(%%rsi),%%ymm11,%%ymm12	\n\t"/* c4 = b4.c5+c6 */\
		"vaddpd	(%%r14),%%ymm3,%%ymm3			\n\t	vaddpd	0x20(%%r14),%%ymm11,%%ymm11		\n\t"/* c3 = b3.c4+c6 */\
	"vfmadd132pd 0x020(%%rsi),%%ymm1,%%ymm7		\n\t vfmadd132pd 0x020(%%rsi),%%ymm9 ,%%ymm15	\n\t"/* c2 = b1.c2+c3 */\
		"vaddpd	(%%r15),%%ymm1,%%ymm1			\n\t	vaddpd	0x20(%%r15),%%ymm9 ,%%ymm9		\n\t"/* c1 = b0.c1+c3 */\
		"vxorpd	%%ymm0,%%ymm0,%%ymm0			\n\t	vxorpd	%%ymm8 ,%%ymm8 ,%%ymm8			\n\t"/* 0 */\
		"vsubpd	%%ymm2,%%ymm0,%%ymm0			\n\t	vsubpd	%%ymm10,%%ymm8 ,%%ymm8			\n\t"/* -c10 */\
		"vaddpd	%%ymm1,%%ymm0,%%ymm0			\n\t	vaddpd	%%ymm9 ,%%ymm8 ,%%ymm8			\n\t"/* c1-c10 */\
		"vaddpd	%%ymm2,%%ymm1,%%ymm1			\n\t	vaddpd	%%ymm10,%%ymm9 ,%%ymm9			\n\t"/* c10+c1 */\
		"vaddpd	%%ymm7,%%ymm0,%%ymm0			\n\t	vaddpd	%%ymm15,%%ymm8 ,%%ymm8			\n\t"/* c1+c2-c10 */\
		"vaddpd	%%ymm2,%%ymm7,%%ymm7			\n\t	vaddpd	%%ymm10,%%ymm15,%%ymm15			\n\t"/* c10+c2 */\
		"vsubpd	%%ymm6,%%ymm1,%%ymm1			\n\t	vsubpd	%%ymm14,%%ymm9 ,%%ymm9			\n\t"/* si1 = c10+c1-c5 */\
		"vaddpd	%%ymm3,%%ymm0,%%ymm0			\n\t	vaddpd	%%ymm11,%%ymm8 ,%%ymm8			\n\t"/* c1+c2+c3-c10 */\
		"vaddpd	%%ymm2,%%ymm3,%%ymm3			\n\t	vaddpd	%%ymm10,%%ymm11,%%ymm11			\n\t"/* c10+c3 */\
		"vsubpd	%%ymm5,%%ymm7,%%ymm7			\n\t	vsubpd	%%ymm13,%%ymm15,%%ymm15			\n\t"/* si5 = c10+c2-c6 */\
		"vaddpd	%%ymm4,%%ymm0,%%ymm0			\n\t	vaddpd	%%ymm12,%%ymm8 ,%%ymm8			\n\t"/* si2 = c1+c2+c3+c4-c10 */\
		"vaddpd	%%ymm2,%%ymm4,%%ymm4			\n\t	vaddpd	%%ymm10,%%ymm12,%%ymm12			\n\t"/* c10+c4 */\
		"vaddpd	%%ymm6,%%ymm3,%%ymm3			\n\t	vaddpd	%%ymm14,%%ymm11,%%ymm11			\n\t"/* si3 = c10+c3+c5 */\
		"vaddpd	%%ymm5,%%ymm4,%%ymm4			\n\t	vaddpd	%%ymm13,%%ymm12,%%ymm12			\n\t"/* si4 = c10+c4+c6 */\
	/* si[1-5] in mm1,0,3,4,7; sr[1-5] in mm9,8,11,12,15: */\
	"movq	%[__o1],%%rax					\n\t"\
	"movq	%[__o2],%%rbx					\n\t"\
	"movq	%[__o3],%%r8					\n\t"\
	"movq	%[__o4],%%r9					\n\t"\
	"movq	%[__o5],%%r10					\n\t"\
		"vmovaps		(%%rax),%%ymm5		\n\t"/* cr1 */\
		"vmovaps		(%%rbx),%%ymm6		\n\t"/* cr2 */\
		"vmovaps		(%%r8 ),%%ymm10		\n\t"/* cr3 */\
		"vmovaps		(%%r9 ),%%ymm13		\n\t"/* cr4 */\
		"vmovaps		(%%r10),%%ymm14		\n\t"/* cr5 */\
		"vsubpd	%%ymm1 ,%%ymm5 ,%%ymm5			\n\t	vaddpd	    (%%rax),%%ymm1,%%ymm1	\n\t"/* B1r = cr1-si1;  Bar = cr1+si1 */\
		"vsubpd	%%ymm0 ,%%ymm6 ,%%ymm6			\n\t	vaddpd	    (%%rbx),%%ymm0,%%ymm0	\n\t"/* B2r = cr2-si2;  B9r = cr2+si2 */\
		"vsubpd	%%ymm3 ,%%ymm10,%%ymm10			\n\t	vaddpd	    (%%r8 ),%%ymm3,%%ymm3	\n\t"/* B3r = cr3-si3;  B8r = cr3+si3 */\
		"vsubpd	%%ymm4 ,%%ymm13,%%ymm13			\n\t	vaddpd	    (%%r9 ),%%ymm4,%%ymm4	\n\t"/* B4r = cr4-si4;  B7r = cr4+si4 */\
		"vsubpd	%%ymm7 ,%%ymm14,%%ymm14			\n\t	vaddpd	    (%%r10),%%ymm7,%%ymm7	\n\t"/* B5r = cr5-si5;  B6r = cr5+si5 */\
		"vmovaps	%%ymm5 ,    (%%rax)			\n\t	vmovaps	%%ymm1,    (%%r11)			\n\t"/* write B1r,Bar */\
		"vmovaps	%%ymm6 ,    (%%rbx)			\n\t	vmovaps	%%ymm0,    (%%r12)			\n\t"/* write B2r,B9r */\
		"vmovaps	%%ymm10,    (%%r8 )			\n\t	vmovaps	%%ymm3,    (%%r13)			\n\t"/* write B3r,B8r */\
		"vmovaps	%%ymm13,    (%%r9 )			\n\t	vmovaps	%%ymm4,    (%%r14)			\n\t"/* write B4r,B7r */\
		"vmovaps	%%ymm14,    (%%r10)			\n\t	vmovaps	%%ymm7,    (%%r15)			\n\t"/* write B5r,B6r */\
		"vmovaps	0x20(%%rax),%%ymm5		\n\t"/* ci1 */\
		"vmovaps	0x20(%%rbx),%%ymm6		\n\t"/* ci2 */\
		"vmovaps	0x20(%%r8 ),%%ymm10		\n\t"/* ci3 */\
		"vmovaps	0x20(%%r9 ),%%ymm13		\n\t"/* ci4 */\
		"vmovaps	0x20(%%r10),%%ymm14		\n\t"/* ci5 */\
		"vsubpd	%%ymm9 ,%%ymm5 ,%%ymm5			\n\t	vaddpd	0x20(%%rax),%%ymm9 ,%%ymm9 	\n\t"/* Bai = ci1-si1;  B1i = ci1+si1 */\
		"vsubpd	%%ymm8 ,%%ymm6 ,%%ymm6			\n\t	vaddpd	0x20(%%rbx),%%ymm8 ,%%ymm8 	\n\t"/* B9i = ci2-si2;  B2i = ci2+si2 */\
		"vsubpd	%%ymm11,%%ymm10,%%ymm10			\n\t	vaddpd	0x20(%%r8 ),%%ymm11,%%ymm11	\n\t"/* B8i = ci3-si3;  B3i = ci3+si3 */\
		"vsubpd	%%ymm12,%%ymm13,%%ymm13			\n\t	vaddpd	0x20(%%r9 ),%%ymm12,%%ymm12	\n\t"/* B7i = ci4-si4;  B4i = ci4+si4 */\
		"vsubpd	%%ymm15,%%ymm14,%%ymm14			\n\t	vaddpd	0x20(%%r10),%%ymm15,%%ymm15	\n\t"/* B6i = ci5-si5;  B5i = ci5+si5 */\
		"vmovaps	%%ymm5 ,0x20(%%r11)			\n\t	vmovaps	%%ymm9 ,0x20(%%rax)			\n\t"/* write Bai,B1i */\
		"vmovaps	%%ymm6 ,0x20(%%r12)			\n\t	vmovaps	%%ymm8 ,0x20(%%rbx)			\n\t"/* write B9i,B2i */\
		"vmovaps	%%ymm10,0x20(%%r13)			\n\t	vmovaps	%%ymm11,0x20(%%r8 )			\n\t"/* write B8i,B3i */\
		"vmovaps	%%ymm13,0x20(%%r14)			\n\t	vmovaps	%%ymm12,0x20(%%r9 )			\n\t"/* write B7i,B4i */\
		"vmovaps	%%ymm14,0x20(%%r15)			\n\t	vmovaps	%%ymm15,0x20(%%r10)			\n\t"/* write B6i,B5i */\
		:					/* outputs: none */\
		: [in0] "m" (Xin0)	/* All inputs from memory addresses here */\
		 ,[i_off] "m" (Xi_off)	/* I-address pointer-strides */\
		 ,[__cc] "m" (Xcc)\
		 ,[out0] "m" (Xout0)\
		 ,[o_off] "m" (Xo_off)	/* O-address pointer-strides */\
		: "cc","memory","rax","rbx","rcx","rdx","rsi","r8","r9","r10","r11","r12","r13","r14","r15",\
			"xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"		/* Clobbered registers */\
		);\
	}

   #endif	// ? DFT_11_FMA

#elif defined(USE_AVX)
	//
	// NB: GCC builds gave "Error: invalid character '(' in mnemonic" for this macro ... I'd forgotten
	// to restore whites to several ...pd (... pairings after having removed it for purposes of col-
	// alignment during AVX-ification of what started as the SSE2 version of the macro:
	//
	// ******> A good way to track down such non-obvious syntax errors output by the assembler (typically w/o
	// helpful localization information such as line numbers) is to insert a single asm line with a different
	// kind of asm error - e.g. a vxorpd with just 2 register args in this case - and move 1 or 2 such around
	// in order to ever-more-closely bracket the first error of the type whose origin is being sought.
	//
	// Arithmetic opcount: [182 ADD, 44 MUL), 154 memref], under our general target of 1 memref per vec_dbl arithmetic op.
	//
	#define SSE2_RADIX_11_DFT(Xin0,Xi_off, Xcc, Xout0,Xo_off)\
	{\
	__asm__ volatile (\
	"movq	%[in0],%%rax		\n\t	movq	%[out0],%%rcx	\n\t"/* rax,rcx hold to I/O base-pointers */\
	"movq	%[i_off],%%rbx		\n\t	movq	%[o_off],%%rdx	\n\t"/* rbx,rdx hold pointers to uint32[] I/O pointer-offset-stride arrays */\
		/********************************************/\
		/*       Here are the 5 cosine terms:       */\
		/********************************************/\
		"movq	%[__cc],%%rsi			\n\t"/* cc1 ... two in -0x40(%%rsi) */\
		/*****************************/				/******************************/\
		/*        Real Parts:        */				/*     Imaginary Parts:       */\
		/*****************************/				/******************************/\
	"movslq	0x04(%%rbx),%%r8	\n\t"/* i1 offset */\
	"movslq	0x28(%%rbx),%%r9	\n\t"/* iA offset */\
	"movslq	0x08(%%rbx),%%r10	\n\t"/* i2 offset */\
	"movslq	0x24(%%rbx),%%r11	\n\t"/* i9 offset */\
	"movslq	0x0c(%%rbx),%%r12	\n\t"/* i3 offset */\
	"movslq	0x20(%%rbx),%%r13	\n\t"/* i8 offset */\
	"movslq	0x10(%%rbx),%%r14	\n\t"/* i4 offset */\
	"movslq	0x1c(%%rbx),%%r15	\n\t"/* i7 offset */\
	"addq	%%rax,%%r8		\n\t	vmovaps (%%r8 ),%%ymm1	\n\t	vmovaps	0x20(%%r8 ),%%ymm9	\n\t"/* x1 */\
	"addq	%%rax,%%r9		\n\t	vmovaps (%%r9 ),%%ymm5	\n\t	vmovaps	0x20(%%r9 ),%%ymm13	\n\t"/* xA */\
	"addq	%%rax,%%r10		\n\t	vmovaps (%%r10),%%ymm2	\n\t	vmovaps	0x20(%%r10),%%ymm10	\n\t"/* x2 */\
	"addq	%%rax,%%r11		\n\t	vmovaps (%%r11),%%ymm6	\n\t	vmovaps	0x20(%%r11),%%ymm14	\n\t"/* x9 */\
	"addq	%%rax,%%r12		\n\t	vmovaps (%%r12),%%ymm3	\n\t	vmovaps	0x20(%%r12),%%ymm11	\n\t"/* x3 */\
	"addq	%%rax,%%r13		\n\t	vmovaps (%%r13),%%ymm7	\n\t	vmovaps	0x20(%%r13),%%ymm15	\n\t"/* x8 */\
	"addq	%%rax,%%r14		\n\t	vmovaps (%%r14),%%ymm4	\n\t	vmovaps	0x20(%%r14),%%ymm12	\n\t"/* x4 */\
	"addq	%%rax,%%r15		\n\t	vmovaps (%%r15),%%ymm0	\n\t	vmovaps	0x20(%%r15),%%ymm8	\n\t"/* x7 */\
		"vsubpd	%%ymm5,%%ymm1,%%ymm1	\n\t			vsubpd	%%ymm13,%%ymm9 ,%%ymm9			\n\t"/* t21 = x1-xA */\
		"vsubpd	%%ymm6,%%ymm2,%%ymm2	\n\t			vsubpd	%%ymm14,%%ymm10,%%ymm10			\n\t"/* t19 = x2-x9 */\
		"vsubpd	%%ymm7,%%ymm3,%%ymm3	\n\t			vsubpd	%%ymm15,%%ymm11,%%ymm11			\n\t"/* t17 = x3-x8 */\
		"vsubpd	%%ymm0,%%ymm4,%%ymm4	\n\t			vsubpd	%%ymm8 ,%%ymm12,%%ymm12			\n\t"/* t15 = x4-x7 */\
	"movslq	0x28(%%rdx),%%r8	\n\t"/* oA offset */\
	"movslq	0x24(%%rdx),%%r9	\n\t"/* o9 offset */\
	"movslq	0x20(%%rdx),%%r10	\n\t"/* o8 offset */\
	"movslq	0x1c(%%rdx),%%r11	\n\t"/* o7 offset */\
	"addq	%%rcx,%%r8		\n\t	vmovaps %%ymm1,0x20(%%r8 )	\n\t	vmovaps %%ymm9 ,(%%r8 )	\n\t"/* ...and store the above 4 pairs */\
	"addq	%%rcx,%%r9		\n\t	vmovaps %%ymm2,0x20(%%r9 )	\n\t	vmovaps %%ymm10,(%%r9 )	\n\t"/* of re/im-part diffs back to mem */\
	"addq	%%rcx,%%r10		\n\t	vmovaps %%ymm3,0x20(%%r10)	\n\t	vmovaps %%ymm11,(%%r10)	\n\t"/* slots oA,8,8,7, swapping the re/im. */\
	"addq	%%rcx,%%r11		\n\t	vmovaps %%ymm4,0x20(%%r11)	\n\t	vmovaps %%ymm12,(%%r11)	\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5	\n\t			vaddpd	%%ymm13,%%ymm13,%%ymm13			\n\t"/* 2.xA */\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6	\n\t			vaddpd	%%ymm14,%%ymm14,%%ymm14			\n\t"/* 2.x9 */\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7	\n\t			vaddpd	%%ymm15,%%ymm15,%%ymm15			\n\t"/* 2.x8 */\
		"vaddpd	%%ymm0,%%ymm0,%%ymm0	\n\t			vaddpd	%%ymm8 ,%%ymm8 ,%%ymm8			\n\t"/* 2.x7 */\
		"vaddpd	%%ymm5,%%ymm1,%%ymm1	\n\t			vaddpd	%%ymm13,%%ymm9 ,%%ymm9			\n\t"/* t3 = x1+xA */\
		"vaddpd	%%ymm6,%%ymm2,%%ymm2	\n\t			vaddpd	%%ymm14,%%ymm10,%%ymm10			\n\t"/* t5 = x2+x9 */\
	"movslq	0x14(%%rbx),%%r9	\n\t"/* i5 offset */\
	"movslq	0x18(%%rbx),%%r11	\n\t"/* i6 offset */\
	"addq	%%rax,%%r9		\n\t	vmovaps (%%r9 ),%%ymm5	\n\t	vmovaps 0x20(%%r9 ),%%ymm13	\n\t"/* x5 */\
	"addq	%%rax,%%r11		\n\t	vmovaps (%%r11),%%ymm6	\n\t	vmovaps 0x20(%%r11),%%ymm14	\n\t"/* x6 */\
		"vaddpd	%%ymm7,%%ymm3,%%ymm3	\n\t			vaddpd	%%ymm15,%%ymm11,%%ymm11			\n\t"/* t7 = x3+x8 */\
		"vaddpd	%%ymm0,%%ymm4,%%ymm4	\n\t			vaddpd	%%ymm8 ,%%ymm12,%%ymm12			\n\t"/* t9 = x4+x7 */\
		"vsubpd	%%ymm6,%%ymm5,%%ymm5	\n\t			vsubpd	%%ymm14,%%ymm13,%%ymm13			\n\t"/* t13 = x5-x6 */\
	"movslq	0x18(%%rdx),%%r8	\n\t"/* o6 offset */\
	"addq	%%rcx,%%r8	\n\t	vmovaps %%ymm5,0x20(%%r8 )	\n\t	vmovaps %%ymm13,(%%r8 )	\n\t"/* store t13 */\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6	\n\t			vaddpd	%%ymm14,%%ymm14,%%ymm14			\n\t"/* 2.x6 */\
		"movslq	(%%rbx),%%r8	\n\t	movslq	(%%rdx),%%r9	\n\t"/* i0,o0 offsets */\
		"addq	%%rax,%%r8		\n\t	addq	%%rcx,%%r9		\n\t"\
		"vmovaps	(%%r8),%%ymm0			\n\t		vmovaps	0x20(%%r8),%%ymm8		\n\t"/* Load 0-input: */\
		"vaddpd	%%ymm6,%%ymm5,%%ymm5	\n\t			vaddpd	%%ymm14,%%ymm13,%%ymm13			\n\t"/* t11 = x5+x6 */\
	"movslq	0x04(%%rdx),%%r11	\n\t	addq	%%rcx,%%r11	\n\t"/* o1 */\
	"movslq	0x08(%%rdx),%%r12	\n\t	addq	%%rcx,%%r12	\n\t"/* o2 */\
	"movslq	0x0c(%%rdx),%%r13	\n\t	addq	%%rcx,%%r13	\n\t"/* o3 */\
	"movslq	0x10(%%rdx),%%r14	\n\t	addq	%%rcx,%%r14	\n\t"/* o4 */\
	"movslq	0x14(%%rdx),%%r15	\n\t	addq	%%rcx,%%r15	\n\t"/* o5 */\
		"vsubpd	%%ymm2,%%ymm1,%%ymm1			\n\t	vsubpd	%%ymm10,%%ymm9 ,%%ymm9			\n\t"/* c1 = t3 -t5 */\
		"vsubpd	%%ymm2,%%ymm5,%%ymm5			\n\t	vsubpd	%%ymm10,%%ymm13,%%ymm13			\n\t"/* c2 = t11-t5 */\
		"vmovaps	%%ymm1,%%ymm6				\n\t	vmovaps	%%ymm9 ,%%ymm14					\n\t"/* cpy c1 */\
		"vsubpd	%%ymm2,%%ymm3,%%ymm3			\n\t	vsubpd	%%ymm10,%%ymm11,%%ymm11			\n\t"/* c4 = t7 -t5 */\
		"vmovaps	%%ymm1,%%ymm7				\n\t	vmovaps	%%ymm9 ,%%ymm15					\n\t"/* cpy c1 */\
		"vmulpd	      (%%rsi),%%ymm1,%%ymm1		\n\t	vmulpd	      (%%rsi),%%ymm9 ,%%ymm9	\n\t"/* a0.c1 */\
		"vmovaps	%%ymm1,(%%r11)				\n\t	vmovaps	%%ymm9 ,0x20(%%r11)				\n\t"/* store copy of a0.c1 in mem */\
		"vsubpd	%%ymm2,%%ymm4,%%ymm4			\n\t	vsubpd	%%ymm10,%%ymm12,%%ymm12			\n\t"/* c5 = t9-t5 */\
		"vmulpd	-0x020(%%rsi),%%ymm2,%%ymm2		\n\t	vmulpd	-0x020(%%rsi),%%ymm10,%%ymm10	\n\t"/* 5.t5 */\
		"vaddpd	%%ymm5,%%ymm6,%%ymm6			\n\t	vaddpd	%%ymm13,%%ymm14,%%ymm14			\n\t"/* c3 = c1+c2 */\
		"vsubpd	%%ymm3,%%ymm7,%%ymm7			\n\t	vsubpd	%%ymm11,%%ymm15,%%ymm15			\n\t"/* c7 = c1-c4 */\
		"vmulpd	 0x0c0(%%rsi),%%ymm7,%%ymm7		\n\t	vmulpd	 0x0c0(%%rsi),%%ymm15,%%ymm15	\n\t"/* a6.c7 */\
		"vmovaps	%%ymm7,(%%r13)				\n\t	vmovaps	%%ymm15,0x20(%%r13)				\n\t"/* store copy of a6.c7 in mem */\
		"vmovaps	%%ymm3,%%ymm1				\n\t	vmovaps	%%ymm11,%%ymm9					\n\t"/* cpy c4 */\
		"vaddpd	%%ymm4,%%ymm3,%%ymm3			\n\t	vaddpd	%%ymm12,%%ymm11,%%ymm11			\n\t"/* c6 = c4+c5 */\
		"vmulpd	 0x060(%%rsi),%%ymm1,%%ymm1		\n\t	vmulpd	 0x060(%%rsi),%%ymm9 ,%%ymm9	\n\t"/* a3.c4 */\
		"vmovaps	%%ymm1,(%%r12)				\n\t	vmovaps	%%ymm9 ,0x20(%%r12)				\n\t"/* store copy of a3.c4 in mem */\
		"vmovaps	%%ymm5,%%ymm7				\n\t	vmovaps	%%ymm13,%%ymm15					\n\t"/* cpy c2 */\
		"vsubpd	%%ymm4,%%ymm5,%%ymm5			\n\t	vsubpd	%%ymm12,%%ymm13,%%ymm13			\n\t"/* c8 = c2-c5 */\
		"vmulpd	 0x080(%%rsi),%%ymm4,%%ymm4		\n\t	vmulpd	 0x080(%%rsi),%%ymm12,%%ymm12	\n\t"/* a4.c5 */\
		"vmulpd	 0x0e0(%%rsi),%%ymm5,%%ymm5		\n\t	vmulpd	 0x0e0(%%rsi),%%ymm13,%%ymm13	\n\t"/* a7.c8 */\
		"vmulpd	 0x020(%%rsi),%%ymm7,%%ymm7		\n\t	vmulpd	 0x020(%%rsi),%%ymm15,%%ymm15	\n\t"/* a1.c2 */\
		"vaddpd	%%ymm3,%%ymm2,%%ymm2			\n\t	vaddpd	%%ymm11,%%ymm10,%%ymm10			\n\t"/* 5.t5+c6 */\
		"vmovaps	%%ymm6,%%ymm1				\n\t	vmovaps	%%ymm14,%%ymm9					\n\t"/* cpy c3 */\
		"vsubpd	%%ymm3,%%ymm6,%%ymm6			\n\t	vsubpd	%%ymm11,%%ymm14,%%ymm14			\n\t"/* c9 = c3-c6 */\
		"vmulpd	 0x100(%%rsi),%%ymm6,%%ymm6		\n\t	vmulpd	 0x100(%%rsi),%%ymm14,%%ymm14	\n\t"/* c9 = a8.c9 */\
		"vmulpd	 0x0a0(%%rsi),%%ymm3,%%ymm3		\n\t	vmulpd	 0x0a0(%%rsi),%%ymm11,%%ymm11	\n\t"/* c6 = a5.c6 */\
		"vaddpd	%%ymm1,%%ymm2,%%ymm2			\n\t	vaddpd	%%ymm9 ,%%ymm10,%%ymm10			\n\t"/* c10= c3 +c6+5.t5 */\
		"vmulpd	 0x040(%%rsi),%%ymm1,%%ymm1		\n\t	vmulpd	 0x040(%%rsi),%%ymm9 ,%%ymm9	\n\t"/* c3 = a2.c3 */\
		"vaddpd	%%ymm2,%%ymm0,%%ymm0			\n\t	vaddpd	%%ymm10,%%ymm8 ,%%ymm8			\n\t"/* B0 = A0+c10 */\
		"vmulpd	 0x120(%%rsi),%%ymm2,%%ymm2		\n\t	vmulpd	 0x120(%%rsi),%%ymm10,%%ymm10	\n\t"/* c10=a9.c10 */\
		"vmovaps	%%ymm0,(%%r9)				\n\t	vmovaps	%%ymm8 ,0x20(%%r9)			\n\t"/* store B0 */\
		"vaddpd	%%ymm0,%%ymm2,%%ymm2			\n\t	vaddpd	%%ymm8 ,%%ymm10,%%ymm10			\n\t"/* c10= a9.c10+B0 */\
		"vaddpd	%%ymm1,%%ymm7,%%ymm7			\n\t	vaddpd	%%ymm9 ,%%ymm15,%%ymm15			\n\t"/* c2 = a1.c2+c3 */\
		"vaddpd	(%%r11),%%ymm1,%%ymm1			\n\t	vaddpd	0x20(%%r11),%%ymm9 ,%%ymm9		\n\t"/* c1 = a0.c1+c3 */\
		"vaddpd	%%ymm3,%%ymm4,%%ymm4			\n\t	vaddpd	%%ymm11,%%ymm12,%%ymm12			\n\t"/* c5 = a4.c5+c6 */\
		"vaddpd	(%%r12),%%ymm3,%%ymm3			\n\t	vaddpd	0x20(%%r12),%%ymm11,%%ymm11		\n\t"/* c4 = a3.c4+c6 */\
		"vaddpd	%%ymm6,%%ymm5,%%ymm5			\n\t	vaddpd	%%ymm14,%%ymm13,%%ymm13			\n\t"/* c8 = a7.c8+c9 */\
		"vaddpd	(%%r13),%%ymm6,%%ymm6			\n\t	vaddpd	0x20(%%r13),%%ymm14,%%ymm14		\n\t"/* c7 = a6.c7+c9 */\
		"vmovaps	%%ymm2,%%ymm0				\n\t	vmovaps	%%ymm10,%%ymm8					\n\t"/* cpy c10 */\
		"vsubpd	%%ymm1,%%ymm2,%%ymm2			\n\t	vsubpd	%%ymm9 ,%%ymm10,%%ymm10			\n\t"/* c10-c1 */\
		"vaddpd	%%ymm0,%%ymm1,%%ymm1			\n\t	vaddpd	%%ymm8 ,%%ymm9 ,%%ymm9			\n\t"/* c10+c1 */\
		"vsubpd	%%ymm7,%%ymm2,%%ymm2			\n\t	vsubpd	%%ymm15,%%ymm10,%%ymm10			\n\t"/* c10-c1-c2 */\
		"vaddpd	%%ymm0,%%ymm7,%%ymm7			\n\t	vaddpd	%%ymm8 ,%%ymm15,%%ymm15			\n\t"/* c10+c2 */\
		"vsubpd	%%ymm6,%%ymm1,%%ymm1			\n\t	vsubpd	%%ymm14,%%ymm9 ,%%ymm9			\n\t"/* cr1 = c10+c1-c7 */\
		"vsubpd	%%ymm3,%%ymm2,%%ymm2			\n\t	vsubpd	%%ymm11,%%ymm10,%%ymm10			\n\t"/* c10-c1-c2-c4 */\
		"vmovaps	%%ymm1,(%%r11)				\n\t	vmovaps	%%ymm9 ,0x20(%%r11)				\n\t"/* write cr1 */\
		"vaddpd	%%ymm0,%%ymm3,%%ymm3			\n\t	vaddpd	%%ymm8 ,%%ymm11,%%ymm11			\n\t"/* c10+c4 */\
		"vsubpd	%%ymm5,%%ymm7,%%ymm7			\n\t	vsubpd	%%ymm13,%%ymm15,%%ymm15			\n\t"/* cr5 = c10+c2-c8 */\
		"vsubpd	%%ymm4,%%ymm2,%%ymm2			\n\t	vsubpd	%%ymm12,%%ymm10,%%ymm10			\n\t"/* cr2 = c10-c1-c2-c4-c5 */\
		"vmovaps	%%ymm7,(%%r15)				\n\t	vmovaps	%%ymm15,0x20(%%r15)				\n\t"/* write cr5 */\
		"vmovaps	%%ymm2,(%%r12)				\n\t	vmovaps	%%ymm10,0x20(%%r12)				\n\t"/* write cr2 */\
		"vaddpd	%%ymm0,%%ymm4,%%ymm4			\n\t	vaddpd	%%ymm8 ,%%ymm12,%%ymm12			\n\t"/* c10+c5 */\
		"vaddpd	%%ymm6,%%ymm3,%%ymm3			\n\t	vaddpd	%%ymm14,%%ymm11,%%ymm11			\n\t"/* cr3 = c10+c4+c7 */\
		"vaddpd	%%ymm5,%%ymm4,%%ymm4			\n\t	vaddpd	%%ymm13,%%ymm12,%%ymm12			\n\t"/* cr4 = c10+c5+c8 */\
		"vmovaps	%%ymm3,(%%r13)				\n\t	vmovaps	%%ymm11,0x20(%%r13)				\n\t"/* write cr3 */\
		"vmovaps	%%ymm4,(%%r14)				\n\t	vmovaps	%%ymm12,0x20(%%r14)				\n\t"/* write cr4 */\
	/***********************************************************************************************************/\
	/* Here are 5 sine terms: Similar to cosines, but t21,19,17,15,13 replace t3,5,7,9,11 and some sign flips: */\
	/***********************************************************************************************************/\
		"addq	$0x140,%%rsi				\n\t"/* b0 */\
		/*****************************/				/******************************/\
		/*     Imaginary Parts:      */				/*        Real Parts:         *//*** Note re/im swap ***/\
		/*****************************/				/******************************/\
	"movslq	0x28(%%rdx),%%r11	\n\t"/* oA */\
	"movslq	0x24(%%rdx),%%r12	\n\t"/* o9 */\
	"movslq	0x20(%%rdx),%%r13	\n\t"/* o8 */\
	"movslq	0x1c(%%rdx),%%r14	\n\t"/* o7 */\
	"movslq	0x18(%%rdx),%%r15	\n\t"/* o6 */\
	"addq	%%rcx,%%r11		\n\t	vmovaps (%%r11),%%ymm1	\n\t	vmovaps 0x20(%%r11),%%ymm9	\n\t"/* t21 */\
	"addq	%%rcx,%%r12		\n\t	vmovaps (%%r12),%%ymm2	\n\t	vmovaps 0x20(%%r12),%%ymm10	\n\t"/* t19 */\
	"addq	%%rcx,%%r13		\n\t	vmovaps (%%r13),%%ymm3	\n\t	vmovaps 0x20(%%r13),%%ymm11	\n\t"/* t17 */\
	"addq	%%rcx,%%r14		\n\t	vmovaps (%%r14),%%ymm4	\n\t	vmovaps 0x20(%%r14),%%ymm12	\n\t"/* t15 */\
	"addq	%%rcx,%%r15		\n\t	vmovaps (%%r15),%%ymm5	\n\t	vmovaps 0x20(%%r15),%%ymm13	\n\t"/* t13 */\
		"vaddpd	%%ymm2,%%ymm1,%%ymm1			\n\t	vaddpd	%%ymm10,%%ymm9 ,%%ymm9			\n\t"/* c1 = t21+t19 */\
		"vaddpd	%%ymm2,%%ymm5,%%ymm5			\n\t	vaddpd	%%ymm10,%%ymm13,%%ymm13			\n\t"/* c2 = t13+t19 */\
		"vmovaps	%%ymm1,%%ymm6				\n\t	vmovaps	%%ymm9 ,%%ymm14					\n\t"/* cpy c1 */\
		"vaddpd	%%ymm2,%%ymm3,%%ymm3			\n\t	vaddpd	%%ymm10,%%ymm11,%%ymm11			\n\t"/* c4 = t17+t19 */\
		"vmovaps	%%ymm1,%%ymm7				\n\t	vmovaps	%%ymm9 ,%%ymm15					\n\t"/* cpy c1 */\
		"vmulpd	      (%%rsi),%%ymm1,%%ymm1		\n\t	vmulpd	      (%%rsi),%%ymm9 ,%%ymm9	\n\t"/* b0.c1 */\
		"vmovaps	%%ymm1,(%%r15)				\n\t	vmovaps	%%ymm9 ,0x20(%%r15)				\n\t"/* store b0.c1 */\
		"vaddpd	%%ymm2,%%ymm4,%%ymm4			\n\t	vaddpd	%%ymm10,%%ymm12,%%ymm12			\n\t"/* c5 = t15+t19 */\
		"vmulpd	-0x160(%%rsi),%%ymm2,%%ymm2		\n\t	vmulpd	-0x160(%%rsi),%%ymm10,%%ymm10	\n\t"/* 5.t19 */\
		"vaddpd	%%ymm5,%%ymm6,%%ymm6			\n\t	vaddpd	%%ymm13,%%ymm14,%%ymm14			\n\t"/* c3 = c1+c2 */\
		"vsubpd	%%ymm3,%%ymm7,%%ymm7			\n\t	vsubpd	%%ymm11,%%ymm15,%%ymm15			\n\t"/* c7 = c1-c4 */\
		"vmulpd	 0x0c0(%%rsi),%%ymm7,%%ymm7		\n\t	vmulpd	 0x0c0(%%rsi),%%ymm15,%%ymm15	\n\t"/* b6.c7 */\
		"vmovaps	%%ymm7,(%%r13)				\n\t	vmovaps	%%ymm15,0x20(%%r13)				\n\t"/* store b6.c7 */\
		"vmovaps	%%ymm3,%%ymm1				\n\t	vmovaps	%%ymm11,%%ymm9					\n\t"/* cpy c4 */\
		"vaddpd	%%ymm4,%%ymm3,%%ymm3			\n\t	vaddpd	%%ymm12,%%ymm11,%%ymm11			\n\t"/* c6 = c4+c5 */\
		"vmulpd	 0x060(%%rsi),%%ymm1,%%ymm1		\n\t	vmulpd	 0x060(%%rsi),%%ymm9 ,%%ymm9	\n\t"/* b3.c4 */\
		"vmovaps	%%ymm1,(%%r14)				\n\t	vmovaps	%%ymm9 ,0x20(%%r14)				\n\t"/* store b3.c4 */\
		"vmovaps	%%ymm5,%%ymm7				\n\t	vmovaps	%%ymm13,%%ymm15					\n\t"/* cpy c2 */\
		"vsubpd	%%ymm4,%%ymm5,%%ymm5			\n\t	vsubpd	%%ymm12,%%ymm13,%%ymm13			\n\t"/* c8 = c2-c5 */\
		"vmulpd	 0x080(%%rsi),%%ymm4,%%ymm4		\n\t	vmulpd	 0x080(%%rsi),%%ymm12,%%ymm12	\n\t"/* b4.c5 */\
		"vmulpd	 0x0e0(%%rsi),%%ymm5,%%ymm5		\n\t	vmulpd	 0x0e0(%%rsi),%%ymm13,%%ymm13	\n\t"/* b7.c8 */\
		"vmulpd	 0x020(%%rsi),%%ymm7,%%ymm7		\n\t	vmulpd	 0x020(%%rsi),%%ymm15,%%ymm15	\n\t"/* b1.c2 */\
		"vsubpd	%%ymm3,%%ymm2,%%ymm2			\n\t	vsubpd	%%ymm11,%%ymm10,%%ymm10			\n\t"/* 5.t19-c6 */\
		"vmovaps	%%ymm6,%%ymm1				\n\t	vmovaps	%%ymm14,%%ymm9					\n\t"/* cpy c3 */\
		"vsubpd	%%ymm3,%%ymm6,%%ymm6			\n\t	vsubpd	%%ymm11,%%ymm14,%%ymm14			\n\t"/* c9 = c3-c6 */\
		"vmulpd	 0x100(%%rsi),%%ymm6,%%ymm6		\n\t	vmulpd	 0x100(%%rsi),%%ymm14,%%ymm14	\n\t"/* c9 = b8.c9 */\
		"vmulpd	 0x0a0(%%rsi),%%ymm3,%%ymm3		\n\t	vmulpd	 0x0a0(%%rsi),%%ymm11,%%ymm11	\n\t"/* c6 = b5.c6 */\
		"vsubpd	%%ymm1,%%ymm2,%%ymm2			\n\t	vsubpd	%%ymm9 ,%%ymm10,%%ymm10			\n\t"/* -c10 = 5.t19-c6-c3 */\
		"vmulpd	 0x040(%%rsi),%%ymm1,%%ymm1		\n\t	vmulpd	 0x040(%%rsi),%%ymm9 ,%%ymm9	\n\t"/* c3 = b2.c3 */\
		"vmulpd	 0x120(%%rsi),%%ymm2,%%ymm2		\n\t	vmulpd	 0x120(%%rsi),%%ymm10,%%ymm10	\n\t"/* c10 = -b9.-c10 */\
		"vaddpd	%%ymm1,%%ymm7,%%ymm7			\n\t	vaddpd	%%ymm9 ,%%ymm15,%%ymm15			\n\t"/* c2 = b1.c2+c3 */\
		"vaddpd	(%%r15),%%ymm1,%%ymm1			\n\t	vaddpd	0x20(%%r15),%%ymm9 ,%%ymm9		\n\t"/* c1 = b0.c1+c3 */\
		"vaddpd	%%ymm3,%%ymm4,%%ymm4			\n\t	vaddpd	%%ymm11,%%ymm12,%%ymm12			\n\t"/* c4 = b4.c5+c6 */\
		"vaddpd	(%%r14),%%ymm3,%%ymm3			\n\t	vaddpd	0x20(%%r14),%%ymm11,%%ymm11		\n\t"/* c3 = b3.c4+c6 */\
		"vaddpd	%%ymm6,%%ymm5,%%ymm5			\n\t	vaddpd	%%ymm14,%%ymm13,%%ymm13			\n\t"/* c6 = b7.c8+c9 */\
		"vaddpd	(%%r13),%%ymm6,%%ymm6			\n\t	vaddpd	0x20(%%r13),%%ymm14,%%ymm14		\n\t"/* c5 = b6.c7+c9 */\
		"vxorpd	%%ymm0,%%ymm0,%%ymm0			\n\t	vxorpd	%%ymm8 ,%%ymm8 ,%%ymm8			\n\t"/* @ */\
		"vsubpd	%%ymm2,%%ymm0,%%ymm0			\n\t	vsubpd	%%ymm10,%%ymm8 ,%%ymm8			\n\t"/* @ */\
		"vaddpd	%%ymm1,%%ymm0,%%ymm0			\n\t	vaddpd	%%ymm9 ,%%ymm8 ,%%ymm8			\n\t"/* @ */\
		"vaddpd	%%ymm2,%%ymm1,%%ymm1			\n\t	vaddpd	%%ymm10,%%ymm9 ,%%ymm9			\n\t"/* @ */\
		"vaddpd	%%ymm7,%%ymm0,%%ymm0			\n\t	vaddpd	%%ymm15,%%ymm8 ,%%ymm8			\n\t"/* @ */\
		"vaddpd	%%ymm2,%%ymm7,%%ymm7			\n\t	vaddpd	%%ymm10,%%ymm15,%%ymm15			\n\t"/* @ */\
		"vsubpd	%%ymm6,%%ymm1,%%ymm1			\n\t	vsubpd	%%ymm14,%%ymm9 ,%%ymm9			\n\t"/* @ */\
		"vaddpd	%%ymm3,%%ymm0,%%ymm0			\n\t	vaddpd	%%ymm11,%%ymm8 ,%%ymm8			\n\t"/* @ */\
		"vaddpd	%%ymm2,%%ymm3,%%ymm3			\n\t	vaddpd	%%ymm10,%%ymm11,%%ymm11			\n\t"/* @ */\
		"vsubpd	%%ymm5,%%ymm7,%%ymm7			\n\t	vsubpd	%%ymm13,%%ymm15,%%ymm15			\n\t"/* @ */\
		"vaddpd	%%ymm4,%%ymm0,%%ymm0			\n\t	vaddpd	%%ymm12,%%ymm8 ,%%ymm8			\n\t"/* @ */\
		"vaddpd	%%ymm2,%%ymm4,%%ymm4			\n\t	vaddpd	%%ymm10,%%ymm12,%%ymm12			\n\t"/* @ */\
		"vaddpd	%%ymm6,%%ymm3,%%ymm3			\n\t	vaddpd	%%ymm14,%%ymm11,%%ymm11			\n\t"/* @ */\
		"vaddpd	%%ymm5,%%ymm4,%%ymm4			\n\t	vaddpd	%%ymm13,%%ymm12,%%ymm12			\n\t"/* @ */\
		"vmovaps	%%ymm1,%%ymm2				\n\t	vmovaps	%%ymm9 ,%%ymm10					\n\t"/* @ */\
	"movslq	0x04(%%rdx),%%rax	\n\t	addq	%%rcx,%%rax	\n\t"/* o1 */\
	"movslq	0x08(%%rdx),%%rbx	\n\t	addq	%%rcx,%%rbx	\n\t"/* o2 */\
	"movslq	0x0c(%%rdx),%%r8	\n\t	addq	%%rcx,%%r8	\n\t"/* o3 */\
	"movslq	0x10(%%rdx),%%r9	\n\t	addq	%%rcx,%%r9	\n\t"/* o4 */\
	"movslq	0x14(%%rdx),%%r10	\n\t	addq	%%rcx,%%r10	\n\t"/* o5 */\
	/* o1,A: */\
		"vaddpd	(%%rax),%%ymm1,%%ymm1			\n\t	vaddpd	0x20(%%rax),%%ymm9 ,%%ymm9		\n\t"/* @ */\
		"vaddpd	%%ymm2,%%ymm2,%%ymm2			\n\t	vaddpd	%%ymm10,%%ymm10,%%ymm10			\n\t"/* @ */\
		"vmovaps	%%ymm1,(%%r11)				\n\t	vmovaps	%%ymm9 ,0x20(%%rax)				\n\t"/* @ */\
		"vsubpd	%%ymm2,%%ymm1,%%ymm1			\n\t	vsubpd	%%ymm10,%%ymm9 ,%%ymm9			\n\t"/* @ */\
		"vmovaps	%%ymm1,(%%rax)				\n\t	vmovaps	%%ymm9 ,0x20(%%r11)				\n\t"/* @ */\
		"vmovaps	%%ymm0,%%ymm5				\n\t	vmovaps	%%ymm8 ,%%ymm13					\n\t"/* @ */\
	/* o2,9: */\
		"vaddpd	(%%rbx),%%ymm0,%%ymm0			\n\t	vaddpd	0x20(%%rbx),%%ymm8 ,%%ymm8		\n\t"/* @ */\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5			\n\t	vaddpd	%%ymm13,%%ymm13,%%ymm13			\n\t"/* @ */\
		"vmovaps	%%ymm0,(%%r12)				\n\t	vmovaps	%%ymm8 ,0x20(%%rbx)				\n\t"/* @ */\
		"vsubpd	%%ymm5,%%ymm0,%%ymm0			\n\t	vsubpd	%%ymm13,%%ymm8 ,%%ymm8			\n\t"/* @ */\
		"vmovaps	%%ymm0,(%%rbx)				\n\t	vmovaps	%%ymm8 ,0x20(%%r12)				\n\t"/* @ */\
		"vmovaps	%%ymm3,%%ymm6				\n\t	vmovaps	%%ymm11,%%ymm14					\n\t"/* @ */\
	/* o3,8: */\
		"vaddpd	(%%r8 ),%%ymm3,%%ymm3			\n\t	vaddpd	0x20(%%r8 ),%%ymm11,%%ymm11		\n\t"/* @ */\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6			\n\t	vaddpd	%%ymm14,%%ymm14,%%ymm14			\n\t"/* @ */\
		"vmovaps	%%ymm3,(%%r13)				\n\t	vmovaps	%%ymm11,0x20(%%r8 )				\n\t"/* @ */\
		"vsubpd	%%ymm6,%%ymm3,%%ymm3			\n\t	vsubpd	%%ymm14,%%ymm11,%%ymm11			\n\t"/* @ */\
		"vmovaps	%%ymm3,(%%r8 )				\n\t	vmovaps	%%ymm11,0x20(%%r13)				\n\t"/* @ */\
		"vmovaps	%%ymm4,%%ymm2				\n\t	vmovaps	%%ymm12,%%ymm10					\n\t"/* @ */\
	/* o4,7: */\
		"vaddpd	(%%r9 ),%%ymm4,%%ymm4			\n\t	vaddpd	0x20(%%r9 ),%%ymm12,%%ymm12		\n\t"/* @ */\
		"vaddpd	%%ymm2,%%ymm2,%%ymm2			\n\t	vaddpd	%%ymm10,%%ymm10,%%ymm10			\n\t"/* @ */\
		"vmovaps	%%ymm4,(%%r14)				\n\t	vmovaps	%%ymm12,0x20(%%r9 )				\n\t"/* @ */\
		"vsubpd	%%ymm2,%%ymm4,%%ymm4			\n\t	vsubpd	%%ymm10,%%ymm12,%%ymm12			\n\t"/* @ */\
		"vmovaps	%%ymm4,(%%r9 )				\n\t	vmovaps	%%ymm12,0x20(%%r14)				\n\t"/* @ */\
		"vmovaps	%%ymm7,%%ymm5				\n\t	vmovaps	%%ymm15,%%ymm13					\n\t"/* @ */\
	/* o5,6: */\
		"vaddpd	(%%r10),%%ymm7,%%ymm7			\n\t	vaddpd	0x20(%%r10),%%ymm15,%%ymm15		\n\t"/* @ */\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5			\n\t	vaddpd	%%ymm13,%%ymm13,%%ymm13			\n\t"/* @ */\
		"vmovaps	%%ymm7,(%%r15)				\n\t	vmovaps	%%ymm15,0x20(%%r10)				\n\t"/* @ */\
		"vsubpd	%%ymm5,%%ymm7,%%ymm7			\n\t	vsubpd	%%ymm13,%%ymm15,%%ymm15			\n\t"/* @ */\
		"vmovaps	%%ymm7,(%%r10)				\n\t	vmovaps	%%ymm15,0x20(%%r15)				\n\t"/* @ */\
		:					/* outputs: none */\
		: [in0] "m" (Xin0)	/* All inputs from memory addresses here */\
		 ,[i_off] "m" (Xi_off)	/* I-address pointer-strides */\
		 ,[__cc] "m" (Xcc)\
		 ,[out0] "m" (Xout0)\
		 ,[o_off] "m" (Xo_off)	/* O-address pointer-strides */\
		: "cc","memory","rax","rbx","rcx","rdx","rsi","r8","r9","r10","r11","r12","r13","r14","r15","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"		/* Clobbered registers */\
		);\
	}

#elif defined(USE_SSE2) && (OS_BITS == 64)

	// 2-Instruction-stream-overlapped 64-bit-ified version of the 32-bit-style ASM macro, using all of xmm0-15
	#define SSE2_RADIX_11_DFT(Xin0,Xi_off, Xcc, Xout0,Xo_off)\
	{\
	__asm__ volatile (\
	"movq	%[in0],%%rax		\n\t	movq	%[out0],%%rcx	\n\t"/* rax,rcx hold to I/O base-pointers */\
	"movq	%[i_off],%%rbx		\n\t	movq	%[o_off],%%rdx	\n\t"/* rbx,rdx hold pointers to uint32[] I/O pointer-offset-stride arrays */\
		/********************************************/\
		/*       Here are the 5 cosine terms:       */\
		/********************************************/\
		"movq	%[__cc],%%rsi			\n\t"\
		/*****************************/				/******************************/\
		/*        Real Parts:        */				/*     Imaginary Parts:       */\
		/*****************************/				/******************************/\
	"movslq	0x04(%%rbx),%%r8	\n\t"/* i1 offset */\
	"movslq	0x28(%%rbx),%%r9	\n\t"/* iA offset */\
	"movslq	0x08(%%rbx),%%r10	\n\t"/* i2 offset */\
	"movslq	0x24(%%rbx),%%r11	\n\t"/* i9 offset */\
	"movslq	0x0c(%%rbx),%%r12	\n\t"/* i3 offset */\
	"movslq	0x20(%%rbx),%%r13	\n\t"/* i8 offset */\
	"movslq	0x10(%%rbx),%%r14	\n\t"/* i4 offset */\
	"movslq	0x1c(%%rbx),%%r15	\n\t"/* i7 offset */\
	"addq	%%rax,%%r8	\n\t	movaps (%%r8 ),%%xmm1	\n\t	movaps	0x10(%%r8 ),%%xmm9	\n\t"\
	"addq	%%rax,%%r9	\n\t	movaps (%%r9 ),%%xmm5	\n\t	movaps	0x10(%%r9 ),%%xmm13	\n\t"\
	"addq	%%rax,%%r10	\n\t	movaps (%%r10),%%xmm2	\n\t	movaps	0x10(%%r10),%%xmm10	\n\t"\
	"addq	%%rax,%%r11	\n\t	movaps (%%r11),%%xmm6	\n\t	movaps	0x10(%%r11),%%xmm14	\n\t"\
	"addq	%%rax,%%r12	\n\t	movaps (%%r12),%%xmm3	\n\t	movaps	0x10(%%r12),%%xmm11	\n\t"\
	"addq	%%rax,%%r13	\n\t	movaps (%%r13),%%xmm7	\n\t	movaps	0x10(%%r13),%%xmm15	\n\t"\
	"addq	%%rax,%%r14	\n\t	movaps (%%r14),%%xmm4	\n\t	movaps	0x10(%%r14),%%xmm12	\n\t"\
	"addq	%%rax,%%r15	\n\t	movaps (%%r15),%%xmm0	\n\t	movaps	0x10(%%r15),%%xmm8	\n\t"\
		"subpd	%%xmm5,%%xmm1			\n\t			subpd	%%xmm13,%%xmm9			\n\t"\
		"subpd	%%xmm6,%%xmm2			\n\t			subpd	%%xmm14,%%xmm10			\n\t"\
		"subpd	%%xmm7,%%xmm3			\n\t			subpd	%%xmm15,%%xmm11			\n\t"\
		"subpd	%%xmm0,%%xmm4			\n\t			subpd	%%xmm8 ,%%xmm12			\n\t"\
	"movslq	0x28(%%rdx),%%r8	\n\t"/* oA offset */\
	"movslq	0x24(%%rdx),%%r9	\n\t"/* o9 offset */\
	"movslq	0x20(%%rdx),%%r10	\n\t"/* o8 offset */\
	"movslq	0x1c(%%rdx),%%r11	\n\t"/* o7 offset */\
	"addq	%%rcx,%%r8	\n\t	movaps %%xmm1,0x10(%%r8 )	\n\t	movaps %%xmm9 ,(%%r8 )	\n\t"\
	"addq	%%rcx,%%r9	\n\t	movaps %%xmm2,0x10(%%r9 )	\n\t	movaps %%xmm10,(%%r9 )	\n\t"\
	"addq	%%rcx,%%r10	\n\t	movaps %%xmm3,0x10(%%r10)	\n\t	movaps %%xmm11,(%%r10)	\n\t"\
	"addq	%%rcx,%%r11	\n\t	movaps %%xmm4,0x10(%%r11)	\n\t	movaps %%xmm12,(%%r11)	\n\t"\
		"addpd	%%xmm5,%%xmm5			\n\t			addpd	%%xmm13,%%xmm13			\n\t"\
		"addpd	%%xmm6,%%xmm6			\n\t			addpd	%%xmm14,%%xmm14			\n\t"\
		"addpd	%%xmm7,%%xmm7			\n\t			addpd	%%xmm15,%%xmm15			\n\t"\
		"addpd	%%xmm0,%%xmm0			\n\t			addpd	%%xmm8 ,%%xmm8			\n\t"\
		"addpd	%%xmm5,%%xmm1			\n\t			addpd	%%xmm13,%%xmm9			\n\t"\
		"addpd	%%xmm6,%%xmm2			\n\t			addpd	%%xmm14,%%xmm10			\n\t"\
	"movslq	0x14(%%rbx),%%r9	\n\t"/* i5 offset */\
	"movslq	0x18(%%rbx),%%r11	\n\t"/* i6 offset */\
	"addq	%%rax,%%r9	\n\t	movaps (%%r9 ),%%xmm5	\n\t	movaps 0x10(%%r9 ),%%xmm13	\n\t"\
	"addq	%%rax,%%r11	\n\t	movaps (%%r11),%%xmm6	\n\t	movaps 0x10(%%r11),%%xmm14	\n\t"\
		"addpd	%%xmm7,%%xmm3			\n\t			addpd	%%xmm15,%%xmm11			\n\t"\
		"addpd	%%xmm0,%%xmm4			\n\t			addpd	%%xmm8 ,%%xmm12			\n\t"\
		"subpd	%%xmm6,%%xmm5			\n\t			subpd	%%xmm14,%%xmm13			\n\t"\
	"movslq	0x18(%%rdx),%%r8	\n\t"/* o6 offset */\
	"addq	%%rcx,%%r8	\n\t	movaps %%xmm5,0x10(%%r8 )	\n\t	movaps %%xmm13,(%%r8 )	\n\t"\
		"addpd	%%xmm6,%%xmm6			\n\t			addpd	%%xmm14,%%xmm14			\n\t"\
		"movslq	(%%rbx),%%r8	\n\t	movslq	(%%rdx),%%r9	\n\t"/* i0,o0 offsets */\
		"addq	%%rax,%%r8		\n\t	addq	%%rcx,%%r9		\n\t"\
		"movaps	(%%r8),%%xmm0			\n\t			movaps	0x10(%%r8),%%xmm8		\n\t"/* Load 0-input: */\
		"addpd	%%xmm6,%%xmm5			\n\t			addpd	%%xmm14,%%xmm13			\n\t"\
	"movslq	0x04(%%rdx),%%r11	\n\t	addq	%%rcx,%%r11	\n\t"/* o1 */\
	"movslq	0x08(%%rdx),%%r12	\n\t	addq	%%rcx,%%r12	\n\t"/* o2 */\
	"movslq	0x0c(%%rdx),%%r13	\n\t	addq	%%rcx,%%r13	\n\t"/* o3 */\
	"movslq	0x10(%%rdx),%%r14	\n\t	addq	%%rcx,%%r14	\n\t"/* o4 */\
	"movslq	0x14(%%rdx),%%r15	\n\t	addq	%%rcx,%%r15	\n\t"/* o5 */\
		"subpd	%%xmm2,%%xmm1			\n\t			subpd	%%xmm10,%%xmm9			\n\t"\
		"subpd	%%xmm2,%%xmm5			\n\t			subpd	%%xmm10,%%xmm13			\n\t"\
		"movaps	%%xmm1,%%xmm6			\n\t			movaps	%%xmm9 ,%%xmm14			\n\t"\
		"subpd	%%xmm2,%%xmm3			\n\t			subpd	%%xmm10,%%xmm11			\n\t"\
		"movaps	%%xmm1,%%xmm7			\n\t			movaps	%%xmm9 ,%%xmm15			\n\t"\
		"mulpd	     (%%rsi),%%xmm1		\n\t			mulpd	     (%%rsi),%%xmm9	\n\t"\
		"movaps	%%xmm1,(%%r11)			\n\t			movaps	%%xmm9 ,0x10(%%r11)		\n\t"\
		"subpd	%%xmm2,%%xmm4			\n\t			subpd	%%xmm10,%%xmm12			\n\t"\
		"mulpd	-0x10(%%rsi),%%xmm2		\n\t			mulpd	-0x10(%%rsi),%%xmm10	\n\t"\
		"addpd	%%xmm5,%%xmm6			\n\t			addpd	%%xmm13,%%xmm14			\n\t"\
		"subpd	%%xmm3,%%xmm7			\n\t			subpd	%%xmm11,%%xmm15			\n\t"\
		"mulpd	 0x60(%%rsi),%%xmm7		\n\t			mulpd	 0x60(%%rsi),%%xmm15	\n\t"\
		"movaps	%%xmm7,(%%r13)			\n\t			movaps	%%xmm15,0x10(%%r13)		\n\t"\
		"movaps	%%xmm3,%%xmm1			\n\t			movaps	%%xmm11,%%xmm9			\n\t"\
		"addpd	%%xmm4,%%xmm3			\n\t			addpd	%%xmm12,%%xmm11			\n\t"\
		"mulpd	 0x30(%%rsi),%%xmm1		\n\t			mulpd	 0x30(%%rsi),%%xmm9	\n\t"\
		"movaps	%%xmm1,(%%r12)			\n\t			movaps	%%xmm9 ,0x10(%%r12)		\n\t"\
		"movaps	%%xmm5,%%xmm7			\n\t			movaps	%%xmm13,%%xmm15			\n\t"\
		"subpd	%%xmm4,%%xmm5			\n\t			subpd	%%xmm12,%%xmm13			\n\t"\
		"mulpd	 0x40(%%rsi),%%xmm4		\n\t			mulpd	 0x40(%%rsi),%%xmm12	\n\t"\
		"mulpd	 0x70(%%rsi),%%xmm5		\n\t			mulpd	 0x70(%%rsi),%%xmm13	\n\t"\
		"mulpd	 0x10(%%rsi),%%xmm7		\n\t			mulpd	 0x10(%%rsi),%%xmm15	\n\t"\
		"addpd	%%xmm3,%%xmm2			\n\t			addpd	%%xmm11,%%xmm10			\n\t"\
		"movaps	%%xmm6,%%xmm1			\n\t			movaps	%%xmm14,%%xmm9			\n\t"\
		"subpd	%%xmm3,%%xmm6			\n\t			subpd	%%xmm11,%%xmm14			\n\t"\
		"mulpd	 0x80(%%rsi),%%xmm6		\n\t			mulpd	 0x80(%%rsi),%%xmm14	\n\t"\
		"mulpd	 0x50(%%rsi),%%xmm3		\n\t			mulpd	 0x50(%%rsi),%%xmm11	\n\t"\
		"addpd	%%xmm1,%%xmm2			\n\t			addpd	%%xmm9 ,%%xmm10			\n\t"\
		"mulpd	 0x20(%%rsi),%%xmm1		\n\t			mulpd	 0x20(%%rsi),%%xmm9	\n\t"\
		"addpd	%%xmm2,%%xmm0			\n\t			addpd	%%xmm10,%%xmm8			\n\t"\
		"mulpd	 0x90(%%rsi),%%xmm2		\n\t			mulpd	 0x90(%%rsi),%%xmm10	\n\t"\
		"movaps	%%xmm0,(%%r9)			\n\t			movaps	%%xmm8 ,0x10(%%r9)		\n\t"/* Write 0-output */\
		"addpd	%%xmm0,%%xmm2			\n\t			addpd	%%xmm8 ,%%xmm10			\n\t"\
		"addpd	%%xmm1,%%xmm7			\n\t			addpd	%%xmm9 ,%%xmm15			\n\t"\
		"addpd	(%%r11),%%xmm1			\n\t			addpd	0x10(%%r11),%%xmm9		\n\t"\
		"addpd	%%xmm3,%%xmm4			\n\t			addpd	%%xmm11,%%xmm12			\n\t"\
		"addpd	(%%r12),%%xmm3			\n\t			addpd	0x10(%%r12),%%xmm11		\n\t"\
		"addpd	%%xmm6,%%xmm5			\n\t			addpd	%%xmm14,%%xmm13			\n\t"\
		"addpd	(%%r13),%%xmm6			\n\t			addpd	0x10(%%r13),%%xmm14		\n\t"\
		"movaps	%%xmm2,%%xmm0			\n\t			movaps	%%xmm10,%%xmm8			\n\t"\
		"subpd	%%xmm1,%%xmm2			\n\t			subpd	%%xmm9 ,%%xmm10			\n\t"\
		"addpd	%%xmm0,%%xmm1			\n\t			addpd	%%xmm8 ,%%xmm9			\n\t"\
		"subpd	%%xmm7,%%xmm2			\n\t			subpd	%%xmm15,%%xmm10			\n\t"\
		"addpd	%%xmm0,%%xmm7			\n\t			addpd	%%xmm8 ,%%xmm15			\n\t"\
		"subpd	%%xmm6,%%xmm1			\n\t			subpd	%%xmm14,%%xmm9			\n\t"\
		"subpd	%%xmm3,%%xmm2			\n\t			subpd	%%xmm11,%%xmm10			\n\t"\
		"movaps	%%xmm1,(%%r11)			\n\t			movaps	%%xmm9 ,0x10(%%r11)		\n\t"\
		"addpd	%%xmm0,%%xmm3			\n\t			addpd	%%xmm8 ,%%xmm11			\n\t"\
		"subpd	%%xmm5,%%xmm7			\n\t			subpd	%%xmm13,%%xmm15			\n\t"\
		"subpd	%%xmm4,%%xmm2			\n\t			subpd	%%xmm12,%%xmm10			\n\t"\
		"movaps	%%xmm7,(%%r15)			\n\t			movaps	%%xmm15,0x10(%%r15)		\n\t"\
		"movaps	%%xmm2,(%%r12)			\n\t			movaps	%%xmm10,0x10(%%r12)		\n\t"\
		"addpd	%%xmm0,%%xmm4			\n\t			addpd	%%xmm8 ,%%xmm12			\n\t"\
		"addpd	%%xmm6,%%xmm3			\n\t			addpd	%%xmm14,%%xmm11			\n\t"\
		"addpd	%%xmm5,%%xmm4			\n\t			addpd	%%xmm13,%%xmm12			\n\t"\
		"movaps	%%xmm3,(%%r13)			\n\t			movaps	%%xmm11,0x10(%%r13)		\n\t"\
		"movaps	%%xmm4,(%%r14)			\n\t			movaps	%%xmm12,0x10(%%r14)		\n\t"\
	/***********************************************************************************************************/\
	/* Here are 5 sine terms: Similar to cosines, but t21,19,17,15,13 replace t3,5,7,9,11 and some sign flips: */\
	/***********************************************************************************************************/\
		"addq	$0xa0,%%rsi				\n\t"\
		/*****************************/				/******************************/\
		/*        Real Parts:        */				/*     Imaginary Parts:       */\
		/*****************************/				/******************************/\
	"movslq	0x28(%%rdx),%%r11	\n\t"/* oA */\
	"movslq	0x24(%%rdx),%%r12	\n\t"/* o9 */\
	"movslq	0x20(%%rdx),%%r13	\n\t"/* o8 */\
	"movslq	0x1c(%%rdx),%%r14	\n\t"/* o7 */\
	"movslq	0x18(%%rdx),%%r15	\n\t"/* o6 */\
	"addq	%%rcx,%%r11	\n\t	movaps (%%r11),%%xmm1	\n\t	movaps 0x10(%%r11),%%xmm9	\n\t"\
	"addq	%%rcx,%%r12	\n\t	movaps (%%r12),%%xmm2	\n\t	movaps 0x10(%%r12),%%xmm10	\n\t"\
	"addq	%%rcx,%%r13	\n\t	movaps (%%r13),%%xmm3	\n\t	movaps 0x10(%%r13),%%xmm11	\n\t"\
	"addq	%%rcx,%%r14	\n\t	movaps (%%r14),%%xmm4	\n\t	movaps 0x10(%%r14),%%xmm12	\n\t"\
	"addq	%%rcx,%%r15	\n\t	movaps (%%r15),%%xmm5	\n\t	movaps 0x10(%%r15),%%xmm13	\n\t"\
		"addpd	%%xmm2,%%xmm1			\n\t			addpd	%%xmm10,%%xmm9			\n\t"\
		"addpd	%%xmm2,%%xmm5			\n\t			addpd	%%xmm10,%%xmm13			\n\t"\
		"movaps	%%xmm1,%%xmm6			\n\t			movaps	%%xmm9 ,%%xmm14			\n\t"\
		"addpd	%%xmm2,%%xmm3			\n\t			addpd	%%xmm10,%%xmm11			\n\t"\
		"movaps	%%xmm1,%%xmm7			\n\t			movaps	%%xmm9 ,%%xmm15			\n\t"\
		"mulpd	     (%%rsi),%%xmm1		\n\t			mulpd	     (%%rsi),%%xmm9	\n\t"\
		"movaps	%%xmm1,(%%r15)			\n\t			movaps	%%xmm9 ,0x10(%%r15)		\n\t"\
		"addpd	%%xmm2,%%xmm4			\n\t			addpd	%%xmm10,%%xmm12			\n\t"\
		"mulpd	-0xb0(%%rsi),%%xmm2		\n\t			mulpd	-0xb0(%%rsi),%%xmm10	\n\t"\
		"addpd	%%xmm5,%%xmm6			\n\t			addpd	%%xmm13,%%xmm14			\n\t"\
		"subpd	%%xmm3,%%xmm7			\n\t			subpd	%%xmm11,%%xmm15			\n\t"\
		"mulpd	 0x60(%%rsi),%%xmm7		\n\t			mulpd	 0x60(%%rsi),%%xmm15	\n\t"\
		"movaps	%%xmm7,(%%r13)			\n\t			movaps	%%xmm15,0x10(%%r13)		\n\t"\
		"movaps	%%xmm3,%%xmm1			\n\t			movaps	%%xmm11,%%xmm9			\n\t"\
		"addpd	%%xmm4,%%xmm3			\n\t			addpd	%%xmm12,%%xmm11			\n\t"\
		"mulpd	 0x30(%%rsi),%%xmm1		\n\t			mulpd	 0x30(%%rsi),%%xmm9	\n\t"\
		"movaps	%%xmm1,(%%r14)			\n\t			movaps	%%xmm9 ,0x10(%%r14)		\n\t"\
		"movaps	%%xmm5,%%xmm7			\n\t			movaps	%%xmm13,%%xmm15			\n\t"\
		"subpd	%%xmm4,%%xmm5			\n\t			subpd	%%xmm12,%%xmm13			\n\t"\
		"mulpd	 0x40(%%rsi),%%xmm4		\n\t			mulpd	 0x40(%%rsi),%%xmm12	\n\t"\
		"mulpd	 0x70(%%rsi),%%xmm5		\n\t			mulpd	 0x70(%%rsi),%%xmm13	\n\t"\
		"mulpd	 0x10(%%rsi),%%xmm7		\n\t			mulpd	 0x10(%%rsi),%%xmm15	\n\t"\
		"subpd	%%xmm3,%%xmm2			\n\t			subpd	%%xmm11,%%xmm10			\n\t"\
		"movaps	%%xmm6,%%xmm1			\n\t			movaps	%%xmm14,%%xmm9			\n\t"\
		"subpd	%%xmm3,%%xmm6			\n\t			subpd	%%xmm11,%%xmm14			\n\t"\
		"mulpd	 0x80(%%rsi),%%xmm6		\n\t			mulpd	 0x80(%%rsi),%%xmm14	\n\t"\
		"mulpd	 0x50(%%rsi),%%xmm3		\n\t			mulpd	 0x50(%%rsi),%%xmm11	\n\t"\
		"subpd	%%xmm1,%%xmm2			\n\t			subpd	%%xmm9 ,%%xmm10			\n\t"\
		"mulpd	 0x20(%%rsi),%%xmm1		\n\t			mulpd	 0x20(%%rsi),%%xmm9	\n\t"\
		"mulpd	 0x90(%%rsi),%%xmm2		\n\t			mulpd	 0x90(%%rsi),%%xmm10	\n\t"\
		"addpd	%%xmm1,%%xmm7			\n\t			addpd	%%xmm9 ,%%xmm15			\n\t"\
		"addpd	(%%r15),%%xmm1			\n\t			addpd	0x10(%%r15),%%xmm9		\n\t"\
		"addpd	%%xmm3,%%xmm4			\n\t			addpd	%%xmm11,%%xmm12			\n\t"\
		"addpd	(%%r14),%%xmm3			\n\t			addpd	0x10(%%r14),%%xmm11		\n\t"\
		"addpd	%%xmm6,%%xmm5			\n\t			addpd	%%xmm14,%%xmm13			\n\t"\
		"addpd	(%%r13),%%xmm6			\n\t			addpd	0x10(%%r13),%%xmm14		\n\t"\
		"xorpd	%%xmm0,%%xmm0			\n\t			xorpd	%%xmm8 ,%%xmm8			\n\t"\
		"subpd	%%xmm2,%%xmm0			\n\t			subpd	%%xmm10,%%xmm8			\n\t"\
		"addpd	%%xmm1,%%xmm0			\n\t			addpd	%%xmm9 ,%%xmm8			\n\t"\
		"addpd	%%xmm2,%%xmm1			\n\t			addpd	%%xmm10,%%xmm9			\n\t"\
		"addpd	%%xmm7,%%xmm0			\n\t			addpd	%%xmm15,%%xmm8			\n\t"\
		"addpd	%%xmm2,%%xmm7			\n\t			addpd	%%xmm10,%%xmm15			\n\t"\
		"subpd	%%xmm6,%%xmm1			\n\t			subpd	%%xmm14,%%xmm9			\n\t"\
		"addpd	%%xmm3,%%xmm0			\n\t			addpd	%%xmm11,%%xmm8			\n\t"\
		"addpd	%%xmm2,%%xmm3			\n\t			addpd	%%xmm10,%%xmm11			\n\t"\
		"subpd	%%xmm5,%%xmm7			\n\t			subpd	%%xmm13,%%xmm15			\n\t"\
		"addpd	%%xmm4,%%xmm0			\n\t			addpd	%%xmm12,%%xmm8			\n\t"\
		"addpd	%%xmm2,%%xmm4			\n\t			addpd	%%xmm10,%%xmm12			\n\t"\
		"addpd	%%xmm6,%%xmm3			\n\t			addpd	%%xmm14,%%xmm11			\n\t"\
		"addpd	%%xmm5,%%xmm4			\n\t			addpd	%%xmm13,%%xmm12			\n\t"\
		"movaps	%%xmm1,%%xmm2			\n\t			movaps	%%xmm9 ,%%xmm10			\n\t"\
	"movslq	0x04(%%rdx),%%rax	\n\t	addq	%%rcx,%%rax	\n\t"/* o1 */\
	"movslq	0x08(%%rdx),%%rbx	\n\t	addq	%%rcx,%%rbx	\n\t"/* o2 */\
	"movslq	0x0c(%%rdx),%%r8	\n\t	addq	%%rcx,%%r8	\n\t"/* o3 */\
	"movslq	0x10(%%rdx),%%r9	\n\t	addq	%%rcx,%%r9	\n\t"/* o4 */\
	"movslq	0x14(%%rdx),%%r10	\n\t	addq	%%rcx,%%r10	\n\t"/* o5 */\
	/* o1: */\
		"addpd	(%%rax),%%xmm1			\n\t			addpd	0x10(%%rax),%%xmm9		\n\t"\
		"addpd	%%xmm2,%%xmm2			\n\t			addpd	%%xmm10,%%xmm10			\n\t"\
		"movaps	%%xmm1,(%%r11)			\n\t			movaps	%%xmm9 ,0x10(%%rax)		\n\t"\
		"subpd	%%xmm2,%%xmm1			\n\t			subpd	%%xmm10,%%xmm9			\n\t"\
		"movaps	%%xmm1,(%%rax)			\n\t			movaps	%%xmm9 ,0x10(%%r11)		\n\t"\
		"movaps	%%xmm0,%%xmm5			\n\t			movaps	%%xmm8 ,%%xmm13			\n\t"\
	/* o2: */\
		"addpd	(%%rbx),%%xmm0			\n\t			addpd	0x10(%%rbx),%%xmm8		\n\t"\
		"addpd	%%xmm5,%%xmm5			\n\t			addpd	%%xmm13,%%xmm13			\n\t"\
		"movaps	%%xmm0,(%%r12)			\n\t			movaps	%%xmm8 ,0x10(%%rbx)		\n\t"\
		"subpd	%%xmm5,%%xmm0			\n\t			subpd	%%xmm13,%%xmm8			\n\t"\
		"movaps	%%xmm0,(%%rbx)			\n\t			movaps	%%xmm8 ,0x10(%%r12)		\n\t"\
		"movaps	%%xmm3,%%xmm6			\n\t			movaps	%%xmm11,%%xmm14			\n\t"\
	/* o3: */\
		"addpd	(%%r8 ),%%xmm3			\n\t			addpd	0x10(%%r8 ),%%xmm11		\n\t"\
		"addpd	%%xmm6,%%xmm6			\n\t			addpd	%%xmm14,%%xmm14			\n\t"\
		"movaps	%%xmm3,(%%r13)			\n\t			movaps	%%xmm11,0x10(%%r8 )		\n\t"\
		"subpd	%%xmm6,%%xmm3			\n\t			subpd	%%xmm14,%%xmm11			\n\t"\
		"movaps	%%xmm3,(%%r8 )			\n\t			movaps	%%xmm11,0x10(%%r13)		\n\t"\
		"movaps	%%xmm4,%%xmm2			\n\t			movaps	%%xmm12,%%xmm10			\n\t"\
	/* o4: */\
		"addpd	(%%r9 ),%%xmm4			\n\t			addpd	0x10(%%r9 ),%%xmm12		\n\t"\
		"addpd	%%xmm2,%%xmm2			\n\t			addpd	%%xmm10,%%xmm10			\n\t"\
		"movaps	%%xmm4,(%%r14)			\n\t			movaps	%%xmm12,0x10(%%r9 )		\n\t"\
		"subpd	%%xmm2,%%xmm4			\n\t			subpd	%%xmm10,%%xmm12			\n\t"\
		"movaps	%%xmm4,(%%r9 )			\n\t			movaps	%%xmm12,0x10(%%r14)		\n\t"\
		"movaps	%%xmm7,%%xmm5			\n\t			movaps	%%xmm15,%%xmm13			\n\t"\
	/* o5: */\
		"addpd	(%%r10),%%xmm7			\n\t			addpd	0x10(%%r10),%%xmm15		\n\t"\
		"addpd	%%xmm5,%%xmm5			\n\t			addpd	%%xmm13,%%xmm13			\n\t"\
		"movaps	%%xmm7,(%%r15)			\n\t			movaps	%%xmm15,0x10(%%r10)		\n\t"\
		"subpd	%%xmm5,%%xmm7			\n\t			subpd	%%xmm13,%%xmm15			\n\t"\
		"movaps	%%xmm7,(%%r10)			\n\t			movaps	%%xmm15,0x10(%%r15)		\n\t"\
		:					/* outputs: none */\
		: [in0] "m" (Xin0)	/* All inputs from memory addresses here */\
		 ,[i_off] "m" (Xi_off)	/* I-address pointer-strides */\
		 ,[__cc] "m" (Xcc)\
		 ,[out0] "m" (Xout0)\
		 ,[o_off] "m" (Xo_off)	/* O-address pointer-strides */\
		: "cc","memory","rax","rbx","rcx","rdx","rsi","r8","r9","r10","r11","r12","r13","r14","r15","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"		/* Clobbered registers */\
		);\
	}

#elif defined(USE_SSE2) && (OS_BITS == 32)

	#error 32-bit OSes no longer supported for SIMD builds!

#else

	#error Unhandled combination of preprocessr flags!

#endif	// x86 simd version ?

#endif	/* radix11_sse_macro_h_included */

