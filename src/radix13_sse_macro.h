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
#ifndef radix13_sse_macro_h_included
#define radix13_sse_macro_h_included

#include "sse2_macro_gcc64.h"

#ifdef USE_ARM_V8_SIMD

	// Cost: [36 ADD, 168 MUL/FMA, 68 LDP/STP], better than general target of 1 vector-memref per vec_dbl arithmetic op.
	#define SSE2_RADIX_13_DFT(Xin0,Xi_off, Xcc, Xout0,Xo_off)\
	{\
	__asm__ volatile (		/*** Rcol does Imaginary Parts: ***/\
		/* Compute &i[0-C] from 32-bit offsets: */\
		"ldr	x15,%[in0]		\n\t	ldr	x16,%[i_off]	\n\t"\
		"ldr	w0,[x16]		\n\t	add	x0,x15,x0	\n\t"/* I0 */\
		"ldp	w1 ,w2 ,[x16,#0x04]	\n\t"\
		"ldp	w3 ,w4 ,[x16,#0x0c]	\n\t"\
		"ldp	w5 ,w6 ,[x16,#0x14]	\n\t"\
		"ldp	w7 ,w8 ,[x16,#0x1c]	\n\t"\
		"ldp	w9 ,w10,[x16,#0x24]	\n\t"\
		"ldp	w11,w12,[x16,#0x2c]	\n\t"\
		"ldr	x13,%[__cc]	\n\t	ldp	q31,q0,[x13,#-0x10]	\n\t"/* v31 = two, v0 = one, but discard the latter */\
	/* ...Do the + parts of the opening wave of radix-2 butterflies: */\
		"add	x1 ,x15,x1			\n\t	ldp	q1 ,q7 ,[x1]	\n\t"/* x1 */\
		"add	x12,x15,x12			\n\t	ldp	q18,q24,[x12]	\n\t"/* xC */\
		"fadd	v1.2d,v1.2d,v18.2d	\n\t	fadd	v7.2d ,v7.2d ,v24.2d 	\n\t"/* t1 = A1 + Ac */\
		"add	x2 ,x15,x2			\n\t	ldp	q2 ,q8 ,[x2]	\n\t"/* x2 */\
		"add	x11,x15,x11			\n\t	ldp	q17,q23,[x11]	\n\t"/* xB */\
		"fadd	v2.2d,v2.2d,v17.2d	\n\t	fadd	v8.2d ,v8.2d ,v23.2d 	\n\t"/* t2 = A2 + Ab */\
		"add	x3 ,x15,x3			\n\t	ldp	q3 ,q9 ,[x3]	\n\t"/* x3 */\
		"add	x10,x15,x10			\n\t	ldp	q16,q22,[x10]	\n\t"/* xA */\
		"fadd	v3.2d,v3.2d,v16.2d	\n\t	fadd	v9.2d ,v9.2d ,v22.2d 	\n\t"/* t3 = A3 + Aa */\
		"add	x4 ,x15,x4			\n\t	ldp	q4 ,q10,[x4]	\n\t"/* x4 */\
		"add	x9 ,x15,x9			\n\t	ldp	q15,q21,[x9]	\n\t"/* x9 */\
		"fadd	v4.2d,v4.2d,v15.2d	\n\t	fadd	v10.2d,v10.2d,v21.2d 	\n\t"/* t4 = A4 + A9 */\
		"add	x5 ,x15,x5			\n\t	ldp	q5 ,q11,[x5]	\n\t"/* x5 */\
		"add	x8 ,x15,x8			\n\t	ldp	q14,q20,[x8]	\n\t"/* x8 */\
		"fadd	v5.2d,v5.2d,v14.2d	\n\t	fadd	v11.2d,v11.2d,v20.2d 	\n\t"/* t5 = A5 + A8 */\
		"add	x6 ,x15,x6			\n\t	ldp	q6 ,q12,[x6]	\n\t"/* x6 */\
		"add	x7 ,x15,x7			\n\t	ldp	q13,q19,[x7]	\n\t"/* x7 */\
		"fadd	v6.2d,v6.2d,v13.2d	\n\t	fadd	v12.2d,v12.2d,v19.2d 	\n\t"/* t6 = A6 + A7 */\
	/* Write resulting lower-half outputs t1-6 back to memory before ensuing FMA-block to compute - parts overwrites each + result: */\
		"stp	q1 ,q7 ,[x1]	\n\t	fmls	v1.2d,v31.2d,v18.2d	\n\t	fmls	v7.2d ,v31.2d,v24.2d 	\n\t"/* tc = A1 - Ac */\
		"stp	q2 ,q8 ,[x2]	\n\t	fmls	v2.2d,v31.2d,v17.2d	\n\t	fmls	v8.2d ,v31.2d,v23.2d 	\n\t"/* tb = A2 - Ab */\
		"stp	q3 ,q9 ,[x3]	\n\t	fmls	v3.2d,v31.2d,v16.2d	\n\t	fmls	v9.2d ,v31.2d,v22.2d 	\n\t"/* ta = A3 - Aa */\
		"stp	q4 ,q10,[x4]	\n\t	fmls	v4.2d,v31.2d,v15.2d	\n\t	fmls	v10.2d,v31.2d,v21.2d 	\n\t"/* t9 = A4 - A9 */\
		"stp	q5 ,q11,[x5]	\n\t	fmls	v5.2d,v31.2d,v14.2d	\n\t	fmls	v11.2d,v31.2d,v20.2d 	\n\t"/* t8 = A5 - A8 */\
		"stp	q6 ,q12,[x6]	\n\t	fmls	v6.2d,v31.2d,v13.2d	\n\t	fmls	v12.2d,v31.2d,v19.2d 	\n\t"/* t7 = A6 - A7 */\
	/* ...And write the upper-half outputs back to memory to free up registers for the 6x6 sine-term-subconvo computation: */\
	/*	"stp	q1 ,q7 ,[x12]	\n\t"// tcr,i */\
		"stp	q2 ,q8 ,[x11]	\n\t"/* tbr,i */\
		"stp	q3 ,q9 ,[x10]	\n\t"/* tar,i */\
		"stp	q4 ,q10,[x9 ]	\n\t"/* t9r,i */\
		"stp	q5 ,q11,[x8 ]	\n\t"/* t8r,i */\
		"stp	q6 ,q12,[x7 ]	\n\t"/* t7r,i */\
/*
	S1 =      ss1*tc + ss2*tb + ss3*ta + ss4*t9 + ss5*t8 + ss6*t7
	S2 =      ss2*tc + ss4*tb + ss6*ta - ss5*t9 - ss3*t8 - ss1*t7
	S3 =      ss3*tc + ss6*tb - ss4*ta - ss1*t9 + ss2*t8 + ss5*t7
	S4 =      ss4*tc - ss5*tb - ss1*ta + ss3*t9 - ss6*t8 - ss2*t7
	S5 =      ss5*tc - ss3*tb + ss2*ta - ss6*t9 - ss1*t8 + ss4*t7
	S6 =      ss6*tc - ss1*tb + ss5*ta - ss2*t9 + ss4*t8 - ss3*t7
*/\
	"add	x13,x13,#0x60	\n\t"/* Incr trig ptr to point to ss1 */\
		"ldp	q13,q14,[x13      ]	\n\t"/* ss1,2 */\
		"ldp	q15,q16,[x13,#0x20]	\n\t"/* ss3,4 */\
		"ldp	q17,q18,[x13,#0x40]	\n\t"/* ss5,6 */\
		/* S1 = tc already in v1,7; save on reg-copies by taking advantage of opening 3-operand MULs: */\
		"fmul	v2.2d,v1.2d,v14.2d 	\n\t	fmul	v8.2d ,v7.2d,v14.2d 	\n\t"/* S2  = ss2*tc; */\
		"fmul	v3.2d,v1.2d,v15.2d 	\n\t	fmul	v9.2d ,v7.2d,v15.2d 	\n\t"/* S3  = ss3*tc; */\
		"fmul	v4.2d,v1.2d,v16.2d 	\n\t	fmul	v10.2d,v7.2d,v16.2d 	\n\t"/* S4  = ss4*tc; */\
		"fmul	v5.2d,v1.2d,v17.2d 	\n\t	fmul	v11.2d,v7.2d,v17.2d 	\n\t"/* S5  = ss5*tc; */\
		"fmul	v6.2d,v1.2d,v18.2d 	\n\t	fmul	v12.2d,v7.2d,v18.2d 	\n\t"/* S6  = ss6*tc; */\
		"fmul	v1.2d,v1.2d,v13.2d 	\n\t	fmul	v7.2d ,v7.2d,v13.2d 	\n\t"/* S1  = ss1*tc; Moved this one to last MUL slot because unmultiplied v1,7 needed by above MULs */\
		"ldp	q19,q20,[x11] 	\n\t"/* load tB */\
		"fmla	v1.2d,v14.2d,v19.2d 	\n\t	fmla	v7.2d ,v14.2d,v20.2d 	\n\t"/* S1 += ss2*tb; */\
		"fmla	v2.2d,v16.2d,v19.2d 	\n\t	fmla	v8.2d ,v16.2d,v20.2d 	\n\t"/* S2 += ss4*tb; */\
		"fmla	v3.2d,v18.2d,v19.2d 	\n\t	fmla	v9.2d ,v18.2d,v20.2d 	\n\t"/* S3 += ss6*tb; */\
		"fmls	v4.2d,v17.2d,v19.2d 	\n\t	fmls	v10.2d,v17.2d,v20.2d 	\n\t"/* S4 -= ss5*tb; */\
		"fmls	v5.2d,v15.2d,v19.2d 	\n\t	fmls	v11.2d,v15.2d,v20.2d 	\n\t"/* S5 -= ss3*tb; */\
		"fmls	v6.2d,v13.2d,v19.2d 	\n\t	fmls	v12.2d,v13.2d,v20.2d 	\n\t"/* S6 -= ss1*tb; */\
		"ldp	q19,q20,[x10] 	\n\t"/* load tA */\
		"fmla	v1.2d,v15.2d,v19.2d 	\n\t	fmla	v7.2d ,v15.2d,v20.2d 	\n\t"/* S1 += ss3*ta; */\
		"fmla	v2.2d,v18.2d,v19.2d 	\n\t	fmla	v8.2d ,v18.2d,v20.2d 	\n\t"/* S2 += ss6*ta; */\
		"fmls	v3.2d,v16.2d,v19.2d 	\n\t	fmls	v9.2d ,v16.2d,v20.2d 	\n\t"/* S3 -= ss4*ta; */\
		"fmls	v4.2d,v13.2d,v19.2d 	\n\t	fmls	v10.2d,v13.2d,v20.2d 	\n\t"/* S4 -= ss1*ta; */\
		"fmla	v5.2d,v14.2d,v19.2d 	\n\t	fmla	v11.2d,v14.2d,v20.2d 	\n\t"/* S5 += ss2*ta; */\
		"fmla	v6.2d,v17.2d,v19.2d 	\n\t	fmla	v12.2d,v17.2d,v20.2d 	\n\t"/* S6 += ss5*ta; */\
		"ldp	q19,q20,[x9] 	\n\t"/* load t9 */\
		"fmla	v1.2d,v16.2d,v19.2d 	\n\t	fmla	v7.2d ,v16.2d,v20.2d 	\n\t"/* S1 += ss4*t9; */\
		"fmls	v2.2d,v17.2d,v19.2d 	\n\t	fmls	v8.2d ,v17.2d,v20.2d 	\n\t"/* S2 -= ss5*t9; */\
		"fmls	v3.2d,v13.2d,v19.2d 	\n\t	fmls	v9.2d ,v13.2d,v20.2d 	\n\t"/* S3 -= ss1*t9; */\
		"fmla	v4.2d,v15.2d,v19.2d 	\n\t	fmla	v10.2d,v15.2d,v20.2d 	\n\t"/* S4 += ss3*t9; */\
		"fmls	v5.2d,v18.2d,v19.2d 	\n\t	fmls	v11.2d,v18.2d,v20.2d 	\n\t"/* S5 -= ss6*t9; */\
		"fmls	v6.2d,v14.2d,v19.2d 	\n\t	fmls	v12.2d,v14.2d,v20.2d 	\n\t"/* S6 -= ss2*t9; */\
		"ldp	q19,q20,[x8] 	\n\t"/* load t8 */\
		"fmla	v1.2d,v17.2d,v19.2d 	\n\t	fmla	v7.2d ,v17.2d,v20.2d 	\n\t"/* S1 += ss5*t8; */\
		"fmls	v2.2d,v15.2d,v19.2d 	\n\t	fmls	v8.2d ,v15.2d,v20.2d 	\n\t"/* S2 -= ss3*t8; */\
		"fmla	v3.2d,v14.2d,v19.2d 	\n\t	fmla	v9.2d ,v14.2d,v20.2d 	\n\t"/* S3 += ss2*t8; */\
		"fmls	v4.2d,v18.2d,v19.2d 	\n\t	fmls	v10.2d,v18.2d,v20.2d 	\n\t"/* S4 -= ss6*t8; */\
		"fmls	v5.2d,v13.2d,v19.2d 	\n\t	fmls	v11.2d,v13.2d,v20.2d 	\n\t"/* S5 -= ss1*t8; */\
		"fmla	v6.2d,v16.2d,v19.2d 	\n\t	fmla	v12.2d,v16.2d,v20.2d 	\n\t"/* S6 += ss4*t8; */\
		"ldp	q19,q20,[x7] 	\n\t"/* load t7 */\
		"fmla	v1.2d,v18.2d,v19.2d 	\n\t	fmla	v7.2d ,v18.2d,v20.2d 	\n\t"/* S1 += ss6*t7; */\
		"fmls	v2.2d,v13.2d,v19.2d 	\n\t	fmls	v8.2d ,v13.2d,v20.2d 	\n\t"/* S2 -= ss1*t7; */\
		"fmla	v3.2d,v17.2d,v19.2d 	\n\t	fmla	v9.2d ,v17.2d,v20.2d 	\n\t"/* S3 += ss5*t7; */\
		"fmls	v4.2d,v14.2d,v19.2d 	\n\t	fmls	v10.2d,v14.2d,v20.2d 	\n\t"/* S4 -= ss2*t7; */\
		"fmla	v5.2d,v16.2d,v19.2d 	\n\t	fmla	v11.2d,v16.2d,v20.2d 	\n\t"/* S5 += ss4*t7; */\
		"fmls	v6.2d,v15.2d,v19.2d 	\n\t	fmls	v12.2d,v15.2d,v20.2d 	\n\t"/* S6 -= ss3*t7; */\
		/* Compute &o[7-C] from 32-bit offsets: */\
		"ldr	x15,%[out0]		\n\t	ldr	x16,%[o_off]	\n\t"\
		"ldp	w11,w12,[x16,#0x2c]	\n\t	add	x11,x15,x11	\n\t	add	x12,x15,x12	\n\t"\
		"ldp	w9 ,w10,[x16,#0x24]	\n\t	add	x9 ,x15,x9	\n\t	add	x10,x15,x10	\n\t"\
		"ldp	w7 ,w8 ,[x16,#0x1c]	\n\t	add	x7 ,x15,x7	\n\t	add	x8 ,x15,x8	\n\t"\
	/* Write sine-term-subconvo outputs back to memory to free up registers for the 6x6 cosine-term-subconvo computation: */\
		"stp	q1 ,q7 ,[x12]	\n\t"/* s1r,i */\
		"stp	q2 ,q8 ,[x11]	\n\t"/* s2r,i */\
		"stp	q3 ,q9 ,[x10]	\n\t"/* s3r,i */\
		"stp	q4 ,q10,[x9 ]	\n\t"/* s4r,i */\
		"stp	q5 ,q11,[x8 ]	\n\t"/* s5r,i */\
		"stp	q6 ,q12,[x7 ]	\n\t"/* s6r,i */\
/*
	C1 = t0 + cc1*t1 + cc2*t2 + cc3*t3 + cc4*t4 + cc5*t5 + cc6*t6
	C2 = t0 + cc2*t1 + cc4*t2 + cc6*t3 + cc5*t4 + cc3*t5 + cc1*t6
	C3 = t0 + cc3*t1 + cc6*t2 + cc4*t3 + cc1*t4 + cc2*t5 + cc5*t6
	C4 = t0 + cc4*t1 + cc5*t2 + cc1*t3 + cc3*t4 + cc6*t5 + cc2*t6
	C5 = t0 + cc5*t1 + cc3*t2 + cc2*t3 + cc6*t4 + cc1*t5 + cc4*t6
	C6 = t0 + cc6*t1 + cc1*t2 + cc5*t3 + cc2*t4 + cc4*t5 + cc3*t6
	B0 = t0 +     t1 +     t2 +     t3 +     t4 +     t5 +     t6	// X0
*/\
	"sub	x13,x13,#0x60	\n\t"/* Decr trig ptr to point to cc1 */\
		"ldp	q13,q14,[x13      ]	\n\t"/* cc1,2 */\
		"ldp	q15,q16,[x13,#0x20]	\n\t"/* cc3,4 */\
		"ldp	q17,q18,[x13,#0x40]	\n\t"/* cc5,6 */\
		"ldp	q1,q7,[x0] 	\n\t"	/* c1 = A0 */\
		"ldr	w0,[x16]		\n\t	add	x0,x15,x0	\n\t"/* Overwrite I0 with O0 */\
		"mov	v2.16b,v1.16b	\n\t	mov	v8.16b ,v7.16b	\n\t"/* c2 = A0 */\
		"mov	v3.16b,v1.16b	\n\t	mov	v9.16b ,v7.16b	\n\t"/* c3 = A0 */\
		"mov	v4.16b,v1.16b	\n\t	mov	v10.16b,v7.16b	\n\t"/* c4 = A0 */\
		"mov	v5.16b,v1.16b	\n\t	mov	v11.16b,v7.16b	\n\t"/* c5 = A0 */\
		"mov	v6.16b,v1.16b	\n\t	mov	v12.16b,v7.16b	\n\t"/* c6 = A0 */\
		"mov	v21.16b,v1.16b	\n\t	mov	v22.16b,v7.16b	\n\t"/* B0 = A0; Init registers used to accumulate DC-component */\
		"ldp	q19,q20,[x1] 	\n\t"/* load t1 */\
		"fmla	v1.2d,v13.2d,v19.2d 	\n\t	fmla	v7.2d ,v13.2d,v20.2d 	\n\t"/* C1 += cc1*tc; */\
		"fmla	v2.2d,v14.2d,v19.2d 	\n\t	fmla	v8.2d ,v14.2d,v20.2d 	\n\t"/* C2 += cc2*tc; */\
		"fmla	v3.2d,v15.2d,v19.2d 	\n\t	fmla	v9.2d ,v15.2d,v20.2d 	\n\t"/* C3 += cc3*tc; */\
		"fmla	v4.2d,v16.2d,v19.2d 	\n\t	fmla	v10.2d,v16.2d,v20.2d 	\n\t"/* C4 += cc4*tc; */\
		"fmla	v5.2d,v17.2d,v19.2d 	\n\t	fmla	v11.2d,v17.2d,v20.2d 	\n\t"/* C5 += cc5*tc; */\
		"fmla	v6.2d,v18.2d,v19.2d 	\n\t	fmla	v12.2d,v18.2d,v20.2d 	\n\t"/* C6 += cc6*tc; */\
		"fadd	v21.2d,v21.2d,v19.2d \n\t	fadd	v22.2d,v22.2d,v20.2d 	\n\t"/* B0 += t1 */\
		"ldp	q19,q20,[x2] 	\n\t"/* load t2 */\
		"fmla	v1.2d,v14.2d,v19.2d 	\n\t	fmla	v7.2d ,v14.2d,v20.2d 	\n\t"/* C1 += cc2*tb; */\
		"fmla	v2.2d,v16.2d,v19.2d 	\n\t	fmla	v8.2d ,v16.2d,v20.2d 	\n\t"/* C2 += cc4*tb; */\
		"fmla	v3.2d,v18.2d,v19.2d 	\n\t	fmla	v9.2d ,v18.2d,v20.2d 	\n\t"/* C3 += cc6*tb; */\
		"fmla	v4.2d,v17.2d,v19.2d 	\n\t	fmla	v10.2d,v17.2d,v20.2d 	\n\t"/* C4 += cc5*tb; */\
		"fmla	v5.2d,v15.2d,v19.2d 	\n\t	fmla	v11.2d,v15.2d,v20.2d 	\n\t"/* C5 += cc3*tb; */\
		"fmla	v6.2d,v13.2d,v19.2d 	\n\t	fmla	v12.2d,v13.2d,v20.2d 	\n\t"/* C6 += cc1*tb; */\
		"fadd	v21.2d,v21.2d,v19.2d \n\t	fadd	v22.2d,v22.2d,v20.2d 	\n\t"/* B0 += t2 */\
		"ldp	q19,q20,[x3] 	\n\t"/* load t3 */\
		"fmla	v1.2d,v15.2d,v19.2d 	\n\t	fmla	v7.2d ,v15.2d,v20.2d 	\n\t"/* C1 += cc3*ta; */\
		"fmla	v2.2d,v18.2d,v19.2d 	\n\t	fmla	v8.2d ,v18.2d,v20.2d 	\n\t"/* C2 += cc6*ta; */\
		"fmla	v3.2d,v16.2d,v19.2d 	\n\t	fmla	v9.2d ,v16.2d,v20.2d 	\n\t"/* C3 += cc4*ta; */\
		"fmla	v4.2d,v13.2d,v19.2d 	\n\t	fmla	v10.2d,v13.2d,v20.2d 	\n\t"/* C4 += cc1*ta; */\
		"fmla	v5.2d,v14.2d,v19.2d 	\n\t	fmla	v11.2d,v14.2d,v20.2d 	\n\t"/* C5 += cc2*ta; */\
		"fmla	v6.2d,v17.2d,v19.2d 	\n\t	fmla	v12.2d,v17.2d,v20.2d 	\n\t"/* C6 += cc5*ta; */\
		"fadd	v21.2d,v21.2d,v19.2d \n\t	fadd	v22.2d,v22.2d,v20.2d 	\n\t"/* B0 += t3 */\
		"ldp	q19,q20,[x4] 	\n\t"/* load t4 */\
		"fmla	v1.2d,v16.2d,v19.2d 	\n\t	fmla	v7.2d ,v16.2d,v20.2d 	\n\t"/* C1 += cc4*t9; */\
		"fmla	v2.2d,v17.2d,v19.2d 	\n\t	fmla	v8.2d ,v17.2d,v20.2d 	\n\t"/* C2 += cc5*t9; */\
		"fmla	v3.2d,v13.2d,v19.2d 	\n\t	fmla	v9.2d ,v13.2d,v20.2d 	\n\t"/* C3 += cc1*t9; */\
		"fmla	v4.2d,v15.2d,v19.2d 	\n\t	fmla	v10.2d,v15.2d,v20.2d 	\n\t"/* C4 += cc3*t9; */\
		"fmla	v5.2d,v18.2d,v19.2d 	\n\t	fmla	v11.2d,v18.2d,v20.2d 	\n\t"/* C5 += cc6*t9; */\
		"fmla	v6.2d,v14.2d,v19.2d 	\n\t	fmla	v12.2d,v14.2d,v20.2d 	\n\t"/* C6 += cc2*t9; */\
		"fadd	v21.2d,v21.2d,v19.2d \n\t	fadd	v22.2d,v22.2d,v20.2d 	\n\t"/* B0 += t4 */\
		"ldp	q19,q20,[x5] 	\n\t"/* load t5 */\
		"fmla	v1.2d,v17.2d,v19.2d 	\n\t	fmla	v7.2d ,v17.2d,v20.2d 	\n\t"/* C1 += cc5*t8; */\
		"fmla	v2.2d,v15.2d,v19.2d 	\n\t	fmla	v8.2d ,v15.2d,v20.2d 	\n\t"/* C2 += cc3*t8; */\
		"fmla	v3.2d,v14.2d,v19.2d 	\n\t	fmla	v9.2d ,v14.2d,v20.2d 	\n\t"/* C3 += cc2*t8; */\
		"fmla	v4.2d,v18.2d,v19.2d 	\n\t	fmla	v10.2d,v18.2d,v20.2d 	\n\t"/* C4 += cc6*t8; */\
		"fmla	v5.2d,v13.2d,v19.2d 	\n\t	fmla	v11.2d,v13.2d,v20.2d 	\n\t"/* C5 += cc1*t8; */\
		"fmla	v6.2d,v16.2d,v19.2d 	\n\t	fmla	v12.2d,v16.2d,v20.2d 	\n\t"/* C6 += cc4*t8; */\
		"fadd	v21.2d,v21.2d,v19.2d \n\t	fadd	v22.2d,v22.2d,v20.2d 	\n\t"/* B0 += t5 */\
		"ldp	q19,q20,[x6] 	\n\t"/* load t6 */\
		"fmla	v1.2d,v18.2d,v19.2d 	\n\t	fmla	v7.2d ,v18.2d,v20.2d 	\n\t"/* C1 += cc6*t7; */\
		"fmla	v2.2d,v13.2d,v19.2d 	\n\t	fmla	v8.2d ,v13.2d,v20.2d 	\n\t"/* C2 += cc1*t7; */\
		"fmla	v3.2d,v17.2d,v19.2d 	\n\t	fmla	v9.2d ,v17.2d,v20.2d 	\n\t"/* C3 += cc5*t7; */\
		"fmla	v4.2d,v14.2d,v19.2d 	\n\t	fmla	v10.2d,v14.2d,v20.2d 	\n\t"/* C4 += cc2*t7; */\
		"fmla	v5.2d,v16.2d,v19.2d 	\n\t	fmla	v11.2d,v16.2d,v20.2d 	\n\t"/* C5 += cc4*t7; */\
		"fmla	v6.2d,v15.2d,v19.2d 	\n\t	fmla	v12.2d,v15.2d,v20.2d 	\n\t"/* C6 += cc3*t7; */\
		"fadd	v21.2d,v21.2d,v19.2d \n\t	fadd	v22.2d,v22.2d,v20.2d 	\n\t"/* B0 += t6 */\
		"stp	q21,q22,[x0]	\n\t"/* Store B0 */\
		/* Compute &o[1-6] from 32-bit offsets: */\
		"ldp	w1,w2,[x16,#0x04]	\n\t	add	x1,x15,x1	\n\t	add	x2,x15,x2	\n\t"\
		"ldp	w3,w4,[x16,#0x0c]	\n\t	add	x3,x15,x3	\n\t	add	x4,x15,x4	\n\t"\
		"ldp	w5,w6,[x16,#0x14]	\n\t	add	x5,x15,x5	\n\t	add	x6,x15,x6	\n\t"\
		/* Reload S-terms: */\
		"ldp	q13,q19,[x12]	\n\t"/* s1r,i */\
		"ldp	q14,q20,[x11]	\n\t"/* s2r,i */\
		"ldp	q15,q21,[x10]	\n\t"/* s3r,i */\
		"ldp	q16,q22,[x9 ]	\n\t"/* s4r,i */\
		"ldp	q17,q23,[x8 ]	\n\t"/* s5r,i */\
		"ldp	q18,q24,[x7 ]	\n\t"/* s6r,i */\
	/* ...Do the lower-half-index portion of the closing wave of radix-2 butterflies: */\
		"fsub	v1.2d,v1.2d,v19.2d 	\n\t	fadd	v7.2d ,v7.2d ,v13.2d 	\n\t"/* B1 = C1 + I*S1 */\
		"fsub	v2.2d,v2.2d,v20.2d 	\n\t	fadd	v8.2d ,v8.2d ,v14.2d 	\n\t"/* B2 = C2 + I*S2 */\
		"fsub	v3.2d,v3.2d,v21.2d 	\n\t	fadd	v9.2d ,v9.2d ,v15.2d 	\n\t"/* B3 = C3 + I*S3 */\
		"fsub	v4.2d,v4.2d,v22.2d 	\n\t	fadd	v10.2d,v10.2d,v16.2d 	\n\t"/* B4 = C4 + I*S4 */\
		"fsub	v5.2d,v5.2d,v23.2d 	\n\t	fadd	v11.2d,v11.2d,v17.2d 	\n\t"/* B5 = C5 + I*S5 */\
		"fsub	v6.2d,v6.2d,v24.2d 	\n\t	fadd	v12.2d,v12.2d,v18.2d 	\n\t"/* B6 = C6 + I*S6 */\
	/* Write resulting lower-half outputs B1-6 back to memory before ensuing FMA-block to compute - parts overwrites each + result: */\
		"stp	q1 ,q7 ,[x1]	\n\t	fmla	v1.2d,v31.2d,v19.2d 	\n\t	fmls	v7.2d ,v31.2d,v13.2d 	\n\t"/* BC = C1 - I*S1 */\
		"stp	q2 ,q8 ,[x2]	\n\t	fmla	v2.2d,v31.2d,v20.2d 	\n\t	fmls	v8.2d ,v31.2d,v14.2d 	\n\t"/* BB = C2 - I*S2 */\
		"stp	q3 ,q9 ,[x3]	\n\t	fmla	v3.2d,v31.2d,v21.2d 	\n\t	fmls	v9.2d ,v31.2d,v15.2d 	\n\t"/* BA = C3 - I*S3 */\
		"stp	q4 ,q10,[x4]	\n\t	fmla	v4.2d,v31.2d,v22.2d 	\n\t	fmls	v10.2d,v31.2d,v16.2d 	\n\t"/* B9 = C4 - I*S4 */\
		"stp	q5 ,q11,[x5]	\n\t	fmla	v5.2d,v31.2d,v23.2d 	\n\t	fmls	v11.2d,v31.2d,v17.2d 	\n\t"/* B8 = C5 - I*S5 */\
		"stp	q6 ,q12,[x6]	\n\t	fmla	v6.2d,v31.2d,v24.2d 	\n\t	fmls	v12.2d,v31.2d,v18.2d 	\n\t"/* B7 = C6 - I*S6 */\
	/* ...And write the upper-half outputs back to memory. */\
		"stp	q1 ,q7 ,[x12]	\n\t"/* BCr,i */\
		"stp	q2 ,q8 ,[x11]	\n\t"/* BBr,i */\
		"stp	q3 ,q9 ,[x10]	\n\t"/* BAr,i */\
		"stp	q4 ,q10,[x9 ]	\n\t"/* B9r,i */\
		"stp	q5 ,q11,[x8 ]	\n\t"/* B8r,i */\
		"stp	q6 ,q12,[x7 ]	\n\t"/* B7r,i */\
		:					/* outputs: none */\
		: [in0] "m" (Xin0)	/* All inputs from memory addresses here */\
		 ,[i_off] "m" (Xi_off)	/* I-address pointer-strides */\
		 ,[__cc] "m" (Xcc)\
		 ,[out0] "m" (Xout0)\
		 ,[o_off] "m" (Xo_off)	/* O-address pointer-strides */\
		: "cc","memory","x0","x1","x2","x3","x4","x5","x6","x7","x8","x9","x10","x11","x12","x13", "x15","x16",\
			"v0","v1","v2","v3","v4","v5","v6","v7","v8","v9","v10","v11","v12","v13","v14","v15",\
			"v16","v17","v18","v19","v20","v21","v22","v23","v24","v31"	/* Clobbered registers - x18,x29 are reserved on Apple platforms */\
		);\
	}

#elif defined(USE_AVX512)	// AVX512 version; cf. comments in AVX2 block below for algo derivations here:

   #if DFT_13_FMA	// [1] DFT_13_FMA = True: Naive but good-ROE and FMA-friendly impl based on RADIX_13_DFT_BASIC in dft_macro.h:

	// 32-vector-regs allows us to replace the sincos-memrefs with reg-data and this save ~100 memrefs.
	// Arithmetic opcount: [36 ADD, 168 FMA (incl 12 MUL), 189 memref], under our general target of 1 memref per vec_dbl arithmetic op.
	#define SSE2_RADIX_13_DFT(Xin0,Xi_off, Xcc, Xout0,Xo_off)\
	{\
	__asm__ volatile (\
	/* Jan 2021: In reduced-#args revamp of this macro, Due to GPR-usage constraints, instead of using a GPR to
	hold vec_dbl*cc (const-offset from two) throughout, broadcast the int64 bitfield for 2.0 across a ZMM, copy
	low qword of result to an MMX via MOVDQ2Q, copy-back-via-MOVQ2DQ/shuffle-propagate-to-rest-of-ZMM as needed.
	***NOTE***: Alas, the use-MMX trick is HUGELY SLOW on my Haswell CPU, looks like Intel significantly sped
	it up on AVX-512-supporting CPUs.
	*/\
	"movq $0x4000000000000000,%%r8	\n\t	vpbroadcastq  %%r8,%%zmm0	\n\t"/* (vec_dbl)2.0 */\
	"movdq2q	%%xmm0,%%mm2		\n\t"/* mm2 holds (double)2, get it? :) */\
	"movq	%[in0],%%r14		\n\t"/* r14 holds Input base-pointer */\
	"movq	%[i_off],%%r15		\n\t"/* r15 holds pointer to uint32[] Outnput pointer-offset-stride array */\
	/* Compute i1-6 from offsets: */\
	"movslq	    (%%r15),%%r8 	\n\t"\
	"movslq	0x04(%%r15),%%rax	\n\t"\
	"movslq	0x08(%%r15),%%rbx	\n\t"\
	"movslq	0x0c(%%r15),%%rcx	\n\t"\
	"movslq	0x10(%%r15),%%rdx	\n\t"\
	"movslq	0x14(%%r15),%%rsi	\n\t"\
	"movslq	0x18(%%r15),%%rdi	\n\t"\
	"addq	%%r14,%%r8 	\n\t	movq	%%r8,%%mm0	\n\t"/* Stash the I0-ptr in MMX for later reload */\
	"addq	%%r14,%%rax	\n\t	vmovaps	(%%rax),%%zmm1	\n\t	vmovaps	0x40(%%rax),%%zmm7 	\n\t"/* A1r,i */\
	"addq	%%r14,%%rbx	\n\t	vmovaps	(%%rbx),%%zmm2	\n\t	vmovaps	0x40(%%rbx),%%zmm8 	\n\t"/* A2r,i */\
	"addq	%%r14,%%rcx	\n\t	vmovaps	(%%rcx),%%zmm3	\n\t	vmovaps	0x40(%%rcx),%%zmm9 	\n\t"/* A3r,i */\
	"addq	%%r14,%%rdx	\n\t	vmovaps	(%%rdx),%%zmm4	\n\t	vmovaps	0x40(%%rdx),%%zmm10	\n\t"/* A4r,i */\
	"addq	%%r14,%%rsi	\n\t	vmovaps	(%%rsi),%%zmm5	\n\t	vmovaps	0x40(%%rsi),%%zmm11	\n\t"/* A5r,i */\
	"addq	%%r14,%%rdi	\n\t	vmovaps	(%%rdi),%%zmm6	\n\t	vmovaps	0x40(%%rdi),%%zmm12	\n\t"/* A6r,i */\
	/* Compute iC-7 from offsets: */\
	"movslq	0x30(%%r15),%%r13	\n\t"\
	"movslq	0x2c(%%r15),%%r12	\n\t"\
	"movslq	0x28(%%r15),%%r11	\n\t"\
	"movslq	0x24(%%r15),%%r10	\n\t"\
	"movslq	0x20(%%r15),%%r9 	\n\t"\
	"movslq	0x1c(%%r15),%%r8 	\n\t"\
	/* ...Do the + parts of the opening wave of radix-2 butterflies: */\
	"addq	%%r14,%%r13	\n\t	vaddpd (%%r13),%%zmm1,%%zmm1	\n\t	vaddpd 0x40(%%r13),%%zmm7 ,%%zmm7 	\n\t"/* t1 = A1 + Ac */\
	"addq	%%r14,%%r12	\n\t	vaddpd (%%r12),%%zmm2,%%zmm2	\n\t	vaddpd 0x40(%%r12),%%zmm8 ,%%zmm8 	\n\t"/* t2 = A2 + Ab */\
	"addq	%%r14,%%r11	\n\t	vaddpd (%%r11),%%zmm3,%%zmm3	\n\t	vaddpd 0x40(%%r11),%%zmm9 ,%%zmm9 	\n\t"/* t3 = A3 + Aa */\
	"addq	%%r14,%%r10	\n\t	vaddpd (%%r10),%%zmm4,%%zmm4	\n\t	vaddpd 0x40(%%r10),%%zmm10,%%zmm10	\n\t"/* t4 = A4 + A9 */\
	"addq	%%r14,%%r9 	\n\t	vaddpd (%%r9 ),%%zmm5,%%zmm5	\n\t	vaddpd 0x40(%%r9 ),%%zmm11,%%zmm11	\n\t"/* t5 = A5 + A8 */\
	"addq	%%r14,%%r8 	\n\t	vaddpd (%%r8 ),%%zmm6,%%zmm6	\n\t	vaddpd 0x40(%%r8 ),%%zmm12,%%zmm12	\n\t"/* t6 = A6 + A7 */\
	/* Write lower-half outputs back to memory... */\
		"vmovaps	%%zmm1,(%%rax)		\n\t"/* t1r */"		vmovaps	%%zmm7 ,0x40(%%rax)	\n\t"/* t1i */\
		"vmovaps	%%zmm2,(%%rbx)		\n\t"/* t2r */"		vmovaps	%%zmm8 ,0x40(%%rbx)	\n\t"/* t2i */\
		"vmovaps	%%zmm3,(%%rcx)		\n\t"/* t3r */"		vmovaps	%%zmm9 ,0x40(%%rcx)	\n\t"/* t3i */\
		"vmovaps	%%zmm4,(%%rdx)		\n\t"/* t4r */"		vmovaps	%%zmm10,0x40(%%rdx)	\n\t"/* t4i */\
		"vmovaps	%%zmm5,(%%rsi)		\n\t"/* t5r */"		vmovaps	%%zmm11,0x40(%%rsi)	\n\t"/* t5i */\
		"vmovaps	%%zmm6,(%%rdi)		\n\t"/* t6r */"		vmovaps	%%zmm12,0x40(%%rdi)	\n\t"/* t6i */\
	/* ...Do the - parts of the radix-2 butterflies: */\
		"vfnmadd231pd (%%r13),%%zmm0,%%zmm1		\n\t	vfnmadd231pd 0x40(%%r13),%%zmm0,%%zmm7 	\n\t"/* ta = A1 - Ac */\
		"vfnmadd231pd (%%r12),%%zmm0,%%zmm2		\n\t	vfnmadd231pd 0x40(%%r12),%%zmm0,%%zmm8 	\n\t"/* t9 = A2 - Ab */\
		"vfnmadd231pd (%%r11),%%zmm0,%%zmm3		\n\t	vfnmadd231pd 0x40(%%r11),%%zmm0,%%zmm9 	\n\t"/* t8 = A3 - Aa */\
		"vfnmadd231pd (%%r10),%%zmm0,%%zmm4		\n\t	vfnmadd231pd 0x40(%%r10),%%zmm0,%%zmm10	\n\t"/* t7 = A4 - A9 */\
		"vfnmadd231pd (%%r9 ),%%zmm0,%%zmm5		\n\t	vfnmadd231pd 0x40(%%r9 ),%%zmm0,%%zmm11	\n\t"/* t6 = A5 - A8 */\
		"vfnmadd231pd (%%r8 ),%%zmm0,%%zmm6		\n\t	vfnmadd231pd 0x40(%%r8 ),%%zmm0,%%zmm12	\n\t"/* t6 = A6 - A7 */\
	/* ...And write the upper-half outputs back to memory to free up registers for the 5x5 sine-term-subconvo computation: */\
		"vmovaps	%%zmm1,(%%r13)		\n\t"/* tar */"		vmovaps	%%zmm7 ,0x40(%%r13)	\n\t"/* tai */\
		"vmovaps	%%zmm2,(%%r12)		\n\t"/* t9r */"		vmovaps	%%zmm8 ,0x40(%%r12)	\n\t"/* t9i */\
		"vmovaps	%%zmm3,(%%r11)		\n\t"/* t8r */"		vmovaps	%%zmm9 ,0x40(%%r11)	\n\t"/* t8i */\
		"vmovaps	%%zmm4,(%%r10)		\n\t"/* t7r */"		vmovaps	%%zmm10,0x40(%%r10)	\n\t"/* t7i */\
		"vmovaps	%%zmm5,(%%r9 )		\n\t"/* t6r */"		vmovaps	%%zmm11,0x40(%%r9 )	\n\t"/* t6i */\
		"vmovaps	%%zmm6,(%%r8 )		\n\t"/* t6r */"		vmovaps	%%zmm12,0x40(%%r8 )	\n\t"/* t6i */\
/*
	S1 =      ss1*tc + ss2*tb + ss3*ta + ss4*t9 + ss5*t8 + ss6*t7
	S2 =      ss2*tc + ss4*tb + ss6*ta - ss5*t9 - ss3*t8 - ss1*t7
	S3 =      ss3*tc + ss6*tb - ss4*ta - ss1*t9 + ss2*t8 + ss5*t7
	S4 =      ss4*tc - ss5*tb - ss1*ta + ss3*t9 - ss6*t8 - ss2*t7
	S5 =      ss5*tc - ss3*tb + ss2*ta - ss6*t9 - ss1*t8 + ss4*t7
	S6 =      ss6*tc - ss1*tb + ss5*ta - ss2*t9 + ss4*t8 - ss3*t7

	In a 16-reg FMA model, can keep the 12 S-terms, both t*r,i-mults and 2 of the 6 ss-mults in-reg, thus
	need 8 loads-from-mem per block.
*/\
	/* Done with rax,rbx,rcx,rdx,rsi,rdi for now, use rsi to hold trig pointer: */\
	"movq	%[__cc],%%rsi	\n\t"/* cc1 */\
		"vmovaps	(%%r13),%%zmm1		\n\t	vmovaps	0x40(%%r13),%%zmm7	\n\t"/* S1 = tc */\
		"vmovaps	%%zmm1,%%zmm2		\n\t	vmovaps	%%zmm7 ,%%zmm8 		\n\t"/* S2 = tc */\
		"vmovaps	%%zmm1,%%zmm3		\n\t	vmovaps	%%zmm7 ,%%zmm9 		\n\t"/* S3 = tc */\
		"vmovaps	%%zmm1,%%zmm4		\n\t	vmovaps	%%zmm7 ,%%zmm10		\n\t"/* S4 = tc */\
		"vmovaps	%%zmm1,%%zmm5		\n\t	vmovaps	%%zmm7 ,%%zmm11		\n\t"/* S5 = tc */\
		"vmovaps	%%zmm1,%%zmm6		\n\t	vmovaps	%%zmm7 ,%%zmm12		\n\t"/* S6 = tc */\
	"addq	$0x180,%%rsi	\n\t"/* Incr trig ptr to point to ss1 */\
	/* AVX-512 gives us 16 more vector-regs, so keep all 6 mults in-registers: */\
		"vmovaps	    (%%rsi),%%zmm13	\n\t	vmovaps	 0x40(%%rsi),%%zmm14	\n\t"/* ss1,ss2 */\
		"vmovaps	0x80(%%rsi),%%zmm16	\n\t	vmovaps	 0xc0(%%rsi),%%zmm17	\n\t"/* ss3,ss4 */\
		"vmovaps   0x100(%%rsi),%%zmm18	\n\t	vmovaps	0x140(%%rsi),%%zmm19	\n\t"/* ss5,ss6 */\
\
	"vmulpd		%%zmm13,%%zmm1,%%zmm1	\n\t	vmulpd		%%zmm13,%%zmm7 ,%%zmm7 	\n\t"/* S1  = ss1*tc */\
	"vmulpd		%%zmm14,%%zmm2,%%zmm2	\n\t	vmulpd		%%zmm14,%%zmm8 ,%%zmm8 	\n\t"/* S2  = ss2*tc */\
	"vmulpd		%%zmm16,%%zmm3,%%zmm3	\n\t	vmulpd		%%zmm16,%%zmm9 ,%%zmm9 	\n\t"/* S3  = ss3*tc */\
	"vmulpd		%%zmm17,%%zmm4,%%zmm4	\n\t	vmulpd		%%zmm17,%%zmm10,%%zmm10	\n\t"/* S4  = ss4*tc */\
	"vmulpd		%%zmm18,%%zmm5,%%zmm5	\n\t	vmulpd		%%zmm18,%%zmm11,%%zmm11	\n\t"/* S5  = ss5*tc */\
	"vmulpd		%%zmm19,%%zmm6,%%zmm6	\n\t	vmulpd		%%zmm19,%%zmm12,%%zmm12	\n\t"/* S6  = ss6*tc */\
\
		"vmovaps	    (%%r12),%%zmm0	\n\t	vmovaps	0x40(%%r12),%%zmm15	\n\t"/* tb */\
	" vfmadd231pd 	%%zmm14,%%zmm0,%%zmm1	\n\t  vfmadd231pd 	%%zmm14,%%zmm15,%%zmm7 	\n\t"/* S1 += ss2*tb */\
	" vfmadd231pd 	%%zmm17,%%zmm0,%%zmm2	\n\t  vfmadd231pd 	%%zmm17,%%zmm15,%%zmm8 	\n\t"/* S2 += ss4*tb */\
	" vfmadd231pd 	%%zmm19,%%zmm0,%%zmm3	\n\t  vfmadd231pd 	%%zmm19,%%zmm15,%%zmm9 	\n\t"/* S3 += ss6*tb */\
	"vfnmadd231pd 	%%zmm18,%%zmm0,%%zmm4	\n\t vfnmadd231pd 	%%zmm18,%%zmm15,%%zmm10	\n\t"/* S4 -= ss5*tb */\
	"vfnmadd231pd 	%%zmm16,%%zmm0,%%zmm5	\n\t vfnmadd231pd 	%%zmm16,%%zmm15,%%zmm11	\n\t"/* S5 -= ss3*tb */\
	"vfnmadd231pd 	%%zmm13,%%zmm0,%%zmm6	\n\t vfnmadd231pd 	%%zmm13,%%zmm15,%%zmm12	\n\t"/* S6 -= ss1*tb */\
\
		"vmovaps	    (%%r11),%%zmm0	\n\t	vmovaps	0x40(%%r11),%%zmm15	\n\t"/* ta */\
	" vfmadd231pd 	%%zmm16,%%zmm0,%%zmm1	\n\t  vfmadd231pd 	%%zmm16,%%zmm15,%%zmm7 	\n\t"/* S1 += ss3*ta */\
	" vfmadd231pd 	%%zmm19,%%zmm0,%%zmm2	\n\t  vfmadd231pd 	%%zmm19,%%zmm15,%%zmm8 	\n\t"/* S2 += ss6*ta */\
	"vfnmadd231pd 	%%zmm17,%%zmm0,%%zmm3	\n\t vfnmadd231pd 	%%zmm17,%%zmm15,%%zmm9 	\n\t"/* S3 -= ss4*ta */\
	"vfnmadd231pd 	%%zmm13,%%zmm0,%%zmm4	\n\t vfnmadd231pd 	%%zmm13,%%zmm15,%%zmm10	\n\t"/* S4 -= ss1*ta */\
	" vfmadd231pd 	%%zmm14,%%zmm0,%%zmm5	\n\t  vfmadd231pd 	%%zmm14,%%zmm15,%%zmm11	\n\t"/* S5 += ss2*ta */\
	" vfmadd231pd 	%%zmm18,%%zmm0,%%zmm6	\n\t  vfmadd231pd 	%%zmm18,%%zmm15,%%zmm12	\n\t"/* S6 += ss5*ta */\
\
		"vmovaps	    (%%r10),%%zmm0	\n\t	vmovaps	0x40(%%r10),%%zmm15	\n\t"/* t9 */\
	" vfmadd231pd 	%%zmm17,%%zmm0,%%zmm1	\n\t  vfmadd231pd 	%%zmm17,%%zmm15,%%zmm7 	\n\t"/* S1 += ss4*t9 */\
	"vfnmadd231pd 	%%zmm18,%%zmm0,%%zmm2	\n\t vfnmadd231pd 	%%zmm18,%%zmm15,%%zmm8 	\n\t"/* S2 -= ss5*t9 */\
	"vfnmadd231pd 	%%zmm13,%%zmm0,%%zmm3	\n\t vfnmadd231pd 	%%zmm13,%%zmm15,%%zmm9 	\n\t"/* S3 -= ss1*t9 */\
	" vfmadd231pd 	%%zmm16,%%zmm0,%%zmm4	\n\t  vfmadd231pd 	%%zmm16,%%zmm15,%%zmm10	\n\t"/* S4 += ss3*t9 */\
	"vfnmadd231pd 	%%zmm19,%%zmm0,%%zmm5	\n\t vfnmadd231pd 	%%zmm19,%%zmm15,%%zmm11	\n\t"/* S5 -= ss6*t9 */\
	"vfnmadd231pd 	%%zmm14,%%zmm0,%%zmm6	\n\t vfnmadd231pd 	%%zmm14,%%zmm15,%%zmm12	\n\t"/* S6 -= ss2*t9 */\
\
		"vmovaps	    (%%r9 ),%%zmm0	\n\t	vmovaps	0x40(%%r9 ),%%zmm15	\n\t"/* t8 */\
	" vfmadd231pd 	%%zmm18,%%zmm0,%%zmm1	\n\t  vfmadd231pd 	%%zmm18,%%zmm15,%%zmm7 	\n\t"/* S1 += ss5*t8 */\
	"vfnmadd231pd 	%%zmm16,%%zmm0,%%zmm2	\n\t vfnmadd231pd 	%%zmm16,%%zmm15,%%zmm8 	\n\t"/* S2 -= ss3*t8 */\
	" vfmadd231pd 	%%zmm14,%%zmm0,%%zmm3	\n\t  vfmadd231pd 	%%zmm14,%%zmm15,%%zmm9 	\n\t"/* S3 += ss2*t8 */\
	"vfnmadd231pd 	%%zmm19,%%zmm0,%%zmm4	\n\t vfnmadd231pd 	%%zmm19,%%zmm15,%%zmm10	\n\t"/* S4 -= ss6*t8 */\
	"vfnmadd231pd 	%%zmm13,%%zmm0,%%zmm5	\n\t vfnmadd231pd 	%%zmm13,%%zmm15,%%zmm11	\n\t"/* S5 -= ss1*t8 */\
	" vfmadd231pd 	%%zmm17,%%zmm0,%%zmm6	\n\t  vfmadd231pd 	%%zmm17,%%zmm15,%%zmm12	\n\t"/* S6 += ss4*t8 */\
\
		"vmovaps	    (%%r8 ),%%zmm0	\n\t	vmovaps	0x40(%%r8 ),%%zmm15	\n\t"/* t7 */\
	" vfmadd231pd 	%%zmm19,%%zmm0,%%zmm1	\n\t  vfmadd231pd 	%%zmm19,%%zmm15,%%zmm7 	\n\t"/* S1 += ss6*t7 */\
	"vfnmadd231pd 	%%zmm13,%%zmm0,%%zmm2	\n\t vfnmadd231pd 	%%zmm13,%%zmm15,%%zmm8 	\n\t"/* S2 -= ss1*t7 */\
	" vfmadd231pd 	%%zmm18,%%zmm0,%%zmm3	\n\t  vfmadd231pd 	%%zmm18,%%zmm15,%%zmm9 	\n\t"/* S3 += ss5*t7 */\
	"vfnmadd231pd 	%%zmm14,%%zmm0,%%zmm4	\n\t vfnmadd231pd 	%%zmm14,%%zmm15,%%zmm10	\n\t"/* S4 -= ss2*t7 */\
	" vfmadd231pd 	%%zmm17,%%zmm0,%%zmm5	\n\t  vfmadd231pd 	%%zmm17,%%zmm15,%%zmm11	\n\t"/* S5 += ss4*t7 */\
	"vfnmadd231pd 	%%zmm16,%%zmm0,%%zmm6	\n\t vfnmadd231pd 	%%zmm16,%%zmm15,%%zmm12	\n\t"/* S6 -= ss3*t7 */\
\
	/* Before overwriting r14,15 with O-side base|offset-ptrs, restore i5 and load A0 into zmm0:*/\
	"movslq	0x14(%%r15),%%rsi	\n\t	addq	%%r14,%%rsi	\n\t"\
	/* Write sine-term-subconvo outputs back to memory to free up registers for the 5x5 cosine-term-subconvo computation: */\
	/* i7-C ptrs no longer needed, so load o-ptrs in their place: */\
	"movq	%[out0],%%r14		\n\t"/* r14 holds Output base-pointer */\
	"movq	%[o_off],%%r15		\n\t"/* r15 holds pointer to uint32[] Outnput pointer-offset-stride array */\
	"movslq	    (%%r15),%%r8 	\n\t"\
	"addq	%%r14,%%r8 	\n\t	movq	%%r8,%%mm1	\n\t"/* Stash the O0-ptr in MMX for later reload */\
	/* Compute oC-7 from offsets: */\
	"movslq	0x30(%%r15),%%r13	\n\t"\
	"movslq	0x2c(%%r15),%%r12	\n\t"\
	"movslq	0x28(%%r15),%%r11	\n\t"\
	"movslq	0x24(%%r15),%%r10	\n\t"\
	"movslq	0x20(%%r15),%%r9 	\n\t"\
	"movslq	0x1c(%%r15),%%r8 	\n\t"\
	"addq	%%r14,%%r13	\n\t	vmovaps	%%zmm1,(%%r13)	\n\t	vmovaps	%%zmm7 ,0x40(%%r13)	\n\t"/* S1 */\
	"addq	%%r14,%%r12	\n\t	vmovaps	%%zmm2,(%%r12)	\n\t	vmovaps	%%zmm8 ,0x40(%%r12)	\n\t"/* S2 */\
	"addq	%%r14,%%r11	\n\t	vmovaps	%%zmm3,(%%r11)	\n\t	vmovaps	%%zmm9 ,0x40(%%r11)	\n\t"/* S3 */\
	"addq	%%r14,%%r10	\n\t	vmovaps	%%zmm4,(%%r10)	\n\t	vmovaps	%%zmm10,0x40(%%r10)	\n\t"/* S4 */\
	"addq	%%r14,%%r9 	\n\t	vmovaps	%%zmm5,(%%r9 )	\n\t	vmovaps	%%zmm11,0x40(%%r9 )	\n\t"/* S5 */\
	"addq	%%r14,%%r8 	\n\t	vmovaps	%%zmm6,(%%r8 )	\n\t	vmovaps	%%zmm12,0x40(%%r8 )	\n\t"/* S6 */\
/*
	C1 = t0 + cc1*t1 + cc2*t2 + cc3*t3 + cc4*t4 + cc5*t5 + cc6*t6
	C2 = t0 + cc2*t1 + cc4*t2 + cc6*t3 + cc5*t4 + cc3*t5 + cc1*t6
	C3 = t0 + cc3*t1 + cc6*t2 + cc4*t3 + cc1*t4 + cc2*t5 + cc5*t6
	C4 = t0 + cc4*t1 + cc5*t2 + cc1*t3 + cc3*t4 + cc6*t5 + cc2*t6
	C5 = t0 + cc5*t1 + cc3*t2 + cc2*t3 + cc6*t4 + cc1*t5 + cc4*t6
	C6 = t0 + cc6*t1 + cc1*t2 + cc5*t3 + cc2*t4 + cc4*t5 + cc3*t6
	B0 = t0 +     t1 +     t2 +     t3 +     t4 +     t5 +     t6	// X0

	In a 16-reg FMA model, makes sense to init the C-terms to the t0-summand,
	and keep the 12 C-terms and 2 of the 6 shared sincos data in registers. That uses 14 regs,
	leaving 2 for the 2 t*[r,i] data shared by each such block. Those need to be in-reg (at least
	in the cosine-mults section) because of the need to accumulate the DC terms: E.g. at end of the
	first block below (after doing the 12 C-term-update FMAs) we add the current values of B0r,i
	(which we init to t0r,i in the ASM version) to t1r,i, then write the result to memory and
	load t2r,i into the same 2 regs in preparation for the next such block.
*/\
/*** rax,rbx,rcx,rdx,rsi,rdi hold i1-6, mm0-2 hold I0,O0,2.0 ***/\
	"movq	%[__cc],%%r15	\n\t"/* cc1 */\
		"movq	%%mm0,%%r8 		\n\t"/* reload I0 into r8 */\
		"vmovaps	(%%r8 ),%%zmm1		\n\t	vmovaps	0x40(%%r8 ),%%zmm7	\n\t"/* c1 = A0 */\
		"vmovaps	%%zmm1,%%zmm2		\n\t	vmovaps	%%zmm7 ,%%zmm8 		\n\t"/* c2 = A0 */\
		"vmovaps	%%zmm1,%%zmm3		\n\t	vmovaps	%%zmm7 ,%%zmm9 		\n\t"/* c3 = A0 */\
		"vmovaps	%%zmm1,%%zmm4		\n\t	vmovaps	%%zmm7 ,%%zmm10		\n\t"/* c4 = A0 */\
		"vmovaps	%%zmm1,%%zmm5		\n\t	vmovaps	%%zmm7 ,%%zmm11		\n\t"/* c5 = A0 */\
		"vmovaps	%%zmm1,%%zmm6		\n\t	vmovaps	%%zmm7 ,%%zmm12		\n\t"/* c6 = A0 */\
\
	/* AVX-512 gives us 16 more vector-regs, so keep all 6 mults in-registers: */\
		"vmovaps	    (%%r15),%%zmm13	\n\t	vmovaps	 0x40(%%r15),%%zmm14	\n\t"/* cc1,cc2 */\
		"vmovaps	0x80(%%r15),%%zmm16	\n\t	vmovaps	 0xc0(%%r15),%%zmm17	\n\t"/* cc3,cc4 */\
		"vmovaps   0x100(%%r15),%%zmm18	\n\t	vmovaps	0x140(%%r15),%%zmm19	\n\t"/* cc5,cc6 */\
\
		"vmovaps	    (%%rax),%%zmm0	\n\t	vmovaps	0x40(%%rax),%%zmm15	\n\t"/* t1 */\
	"vfmadd231pd 	%%zmm13,%%zmm0,%%zmm1	\n\t vfmadd231pd 	%%zmm13,%%zmm15,%%zmm7 	\n\t"/* C1 += cc1*t1 */\
	"vfmadd231pd 	%%zmm14,%%zmm0,%%zmm2	\n\t vfmadd231pd 	%%zmm14,%%zmm15,%%zmm8 	\n\t"/* C2 += cc2*t1 */\
	"vfmadd231pd 	%%zmm16,%%zmm0,%%zmm3	\n\t vfmadd231pd 	%%zmm16,%%zmm15,%%zmm9 	\n\t"/* C3 += cc3*t1 */\
	"vfmadd231pd 	%%zmm17,%%zmm0,%%zmm4	\n\t vfmadd231pd 	%%zmm17,%%zmm15,%%zmm10	\n\t"/* C4 += cc4*t1 */\
	"vfmadd231pd 	%%zmm18,%%zmm0,%%zmm5	\n\t vfmadd231pd 	%%zmm18,%%zmm15,%%zmm11	\n\t"/* C5 += cc5*t1 */\
	"vfmadd231pd 	%%zmm19,%%zmm0,%%zmm6	\n\t vfmadd231pd 	%%zmm19,%%zmm15,%%zmm12	\n\t"/* C6 += cc6*t1 */\
		"vaddpd (%%r8 ),%%zmm0,%%zmm0	\n\t	vaddpd 0x40(%%r8 ),%%zmm15,%%zmm15	\n\t"/* B0 += t1; */\
		"vmovaps	%%zmm0,(%%r8 )		\n\t	vmovaps	%%zmm15,0x40(%%r8 )			\n\t"/* Store B0 */\
\
		"vmovaps	    (%%rbx),%%zmm0	\n\t	vmovaps	0x40(%%rbx),%%zmm15	\n\t"/* t2 */\
	"vfmadd231pd 	%%zmm14,%%zmm0,%%zmm1	\n\t vfmadd231pd 	%%zmm14,%%zmm15,%%zmm7 	\n\t"/* C1 += cc2*t2 */\
	"vfmadd231pd 	%%zmm17,%%zmm0,%%zmm2	\n\t vfmadd231pd 	%%zmm17,%%zmm15,%%zmm8 	\n\t"/* C2 += cc4*t2 */\
	"vfmadd231pd 	%%zmm19,%%zmm0,%%zmm3	\n\t vfmadd231pd 	%%zmm19,%%zmm15,%%zmm9 	\n\t"/* C3 += cc6*t2 */\
	"vfmadd231pd 	%%zmm18,%%zmm0,%%zmm4	\n\t vfmadd231pd 	%%zmm18,%%zmm15,%%zmm10	\n\t"/* C4 += cc5*t2 */\
	"vfmadd231pd 	%%zmm16,%%zmm0,%%zmm5	\n\t vfmadd231pd 	%%zmm16,%%zmm15,%%zmm11	\n\t"/* C5 += cc3*t2 */\
	"vfmadd231pd 	%%zmm13,%%zmm0,%%zmm6	\n\t vfmadd231pd 	%%zmm13,%%zmm15,%%zmm12	\n\t"/* C6 += cc1*t2 */\
		"vaddpd (%%r8 ),%%zmm0,%%zmm0	\n\t	vaddpd 0x40(%%r8 ),%%zmm15,%%zmm15	\n\t"/* B0 += t2; */\
		"vmovaps	%%zmm0,(%%r8 )		\n\t	vmovaps	%%zmm15,0x40(%%r8 )			\n\t"/* Store B0 */\
\
		"vmovaps	    (%%rcx),%%zmm0	\n\t	vmovaps	0x40(%%rcx),%%zmm15	\n\t"/* t3 */\
	"vfmadd231pd 	%%zmm16,%%zmm0,%%zmm1	\n\t vfmadd231pd 	%%zmm16,%%zmm15,%%zmm7 	\n\t"/* C1 += cc3*t3 */\
	"vfmadd231pd 	%%zmm19,%%zmm0,%%zmm2	\n\t vfmadd231pd 	%%zmm19,%%zmm15,%%zmm8 	\n\t"/* C2 += cc6*t3 */\
	"vfmadd231pd 	%%zmm17,%%zmm0,%%zmm3	\n\t vfmadd231pd 	%%zmm17,%%zmm15,%%zmm9 	\n\t"/* C3 += cc4*t3 */\
	"vfmadd231pd 	%%zmm13,%%zmm0,%%zmm4	\n\t vfmadd231pd 	%%zmm13,%%zmm15,%%zmm10	\n\t"/* C4 += cc1*t3 */\
	"vfmadd231pd 	%%zmm14,%%zmm0,%%zmm5	\n\t vfmadd231pd 	%%zmm14,%%zmm15,%%zmm11	\n\t"/* C5 += cc2*t3 */\
	"vfmadd231pd 	%%zmm18,%%zmm0,%%zmm6	\n\t vfmadd231pd 	%%zmm18,%%zmm15,%%zmm12	\n\t"/* C6 += cc5*t3 */\
		"vaddpd (%%r8 ),%%zmm0,%%zmm0	\n\t	vaddpd 0x40(%%r8 ),%%zmm15,%%zmm15	\n\t"/* B0 += t3; */\
		"vmovaps	%%zmm0,(%%r8 )		\n\t	vmovaps	%%zmm15,0x40(%%r8 )			\n\t"/* Store B0 */\
\
		"vmovaps	    (%%rdx),%%zmm0	\n\t	vmovaps	0x40(%%rdx),%%zmm15	\n\t"/* t4 */\
	"vfmadd231pd 	%%zmm17,%%zmm0,%%zmm1	\n\t vfmadd231pd 	%%zmm17,%%zmm15,%%zmm7 	\n\t"/* C1 += cc4*t4 */\
	"vfmadd231pd 	%%zmm18,%%zmm0,%%zmm2	\n\t vfmadd231pd 	%%zmm18,%%zmm15,%%zmm8 	\n\t"/* C2 += cc5*t4 */\
	"vfmadd231pd 	%%zmm13,%%zmm0,%%zmm3	\n\t vfmadd231pd 	%%zmm13,%%zmm15,%%zmm9 	\n\t"/* C3 += cc1*t4 */\
	"vfmadd231pd 	%%zmm16,%%zmm0,%%zmm4	\n\t vfmadd231pd 	%%zmm16,%%zmm15,%%zmm10	\n\t"/* C4 += cc3*t4 */\
	"vfmadd231pd 	%%zmm19,%%zmm0,%%zmm5	\n\t vfmadd231pd 	%%zmm19,%%zmm15,%%zmm11	\n\t"/* C5 += cc6*t4 */\
	"vfmadd231pd 	%%zmm14,%%zmm0,%%zmm6	\n\t vfmadd231pd 	%%zmm14,%%zmm15,%%zmm12	\n\t"/* C6 += cc2*t4 */\
		"vaddpd (%%r8 ),%%zmm0,%%zmm0	\n\t	vaddpd 0x40(%%r8 ),%%zmm15,%%zmm15	\n\t"/* B0 += t4; */\
		"vmovaps	%%zmm0,(%%r8 )		\n\t	vmovaps	%%zmm15,0x40(%%r8 )			\n\t"/* Store B0 */\
\
		"vmovaps	    (%%rsi),%%zmm0	\n\t	vmovaps	0x40(%%rsi),%%zmm15	\n\t"/* t5 */\
	"vfmadd231pd 	%%zmm18,%%zmm0,%%zmm1	\n\t vfmadd231pd 	%%zmm18,%%zmm15,%%zmm7 	\n\t"/* C1 += cc5*t5 */\
	"vfmadd231pd 	%%zmm16,%%zmm0,%%zmm2	\n\t vfmadd231pd 	%%zmm16,%%zmm15,%%zmm8 	\n\t"/* C2 += cc3*t5 */\
	"vfmadd231pd 	%%zmm14,%%zmm0,%%zmm3	\n\t vfmadd231pd 	%%zmm14,%%zmm15,%%zmm9 	\n\t"/* C3 += cc2*t5 */\
	"vfmadd231pd 	%%zmm19,%%zmm0,%%zmm4	\n\t vfmadd231pd 	%%zmm19,%%zmm15,%%zmm10	\n\t"/* C4 += cc6*t5 */\
	"vfmadd231pd 	%%zmm13,%%zmm0,%%zmm5	\n\t vfmadd231pd 	%%zmm13,%%zmm15,%%zmm11	\n\t"/* C5 += cc1*t5 */\
	"vfmadd231pd 	%%zmm17,%%zmm0,%%zmm6	\n\t vfmadd231pd 	%%zmm17,%%zmm15,%%zmm12	\n\t"/* C6 += cc4*t5 */\
		"vaddpd (%%r8 ),%%zmm0,%%zmm0	\n\t	vaddpd 0x40(%%r8 ),%%zmm15,%%zmm15	\n\t"/* B0 += t5; */\
		"vmovaps	%%zmm0,(%%r8 )		\n\t	vmovaps	%%zmm15,0x40(%%r8 )			\n\t"/* Store B0 */\
\
		"vmovaps	    (%%rdi),%%zmm0	\n\t	vmovaps	0x40(%%rdi),%%zmm15	\n\t"/* t6 */\
	"vfmadd231pd 	%%zmm19,%%zmm0,%%zmm1	\n\t vfmadd231pd 	%%zmm19,%%zmm15,%%zmm7 	\n\t"/* C1 += cc6*t6 */\
	"vfmadd231pd 	%%zmm13,%%zmm0,%%zmm2	\n\t vfmadd231pd 	%%zmm13,%%zmm15,%%zmm8 	\n\t"/* C2 += cc1*t6 */\
	"vfmadd231pd 	%%zmm18,%%zmm0,%%zmm3	\n\t vfmadd231pd 	%%zmm18,%%zmm15,%%zmm9 	\n\t"/* C3 += cc5*t6 */\
	"vfmadd231pd 	%%zmm14,%%zmm0,%%zmm4	\n\t vfmadd231pd 	%%zmm14,%%zmm15,%%zmm10	\n\t"/* C4 += cc2*t6 */\
	"vfmadd231pd 	%%zmm17,%%zmm0,%%zmm5	\n\t vfmadd231pd 	%%zmm17,%%zmm15,%%zmm11	\n\t"/* C5 += cc4*t6 */\
	"vfmadd231pd 	%%zmm16,%%zmm0,%%zmm6	\n\t vfmadd231pd 	%%zmm16,%%zmm15,%%zmm12	\n\t"/* C6 += cc3*t6 */\
		"vaddpd (%%r8 ),%%zmm0,%%zmm0	\n\t	vaddpd 0x40(%%r8 ),%%zmm15,%%zmm15	\n\t"/* B0 += t6; */\
		"movq	%%mm1,%%r8 		\n\t"/* reload O0 into r8 */\
		"vmovaps	%%zmm0,(%%r8 )		\n\t	vmovaps	%%zmm15,0x40(%%r8 )			\n\t"/* Store B0, now into output slot */\
\
	"movq	%[o_off],%%r15		\n\t"/* r15 holds pointer to uint32[] Outnput pointer-offset-stride array */\
	"movslq	    (%%r15),%%r8 	\n\t"\
	/* Compute o1-6 from offsets and recompute o7; o8-c still in r9-13): */\
	"movslq	0x04(%%r15),%%rax	\n\t"\
	"movslq	0x08(%%r15),%%rbx	\n\t"\
	"movslq	0x0c(%%r15),%%rcx	\n\t"\
	"movslq	0x10(%%r15),%%rdx	\n\t"\
	"movslq	0x14(%%r15),%%rsi	\n\t"\
	"movslq	0x18(%%r15),%%rdi	\n\t"\
	"movslq	0x1c(%%r15),%%r8 	\n\t"\
	"addq	%%r14,%%r8 	\n\t"\
	"addq	%%r14,%%rax	\n\t"\
	"addq	%%r14,%%rbx	\n\t"\
	"addq	%%r14,%%rcx	\n\t"	/* Have enough registers to also hold 3 of the 6 real parts of S-terms: */\
	"addq	%%r14,%%rdx	\n\t	vmovaps	(%%r10),%%zmm13		\n\t"/* s4r */\
	"addq	%%r14,%%rsi	\n\t	vmovaps	(%%r9 ),%%zmm14		\n\t"/* s5r */\
	"addq	%%r14,%%rdi	\n\t	vmovaps	(%%r8 ),%%zmm15		\n\t"/* s5r */\
	/* ...Do the lower-half-index portion of the closing wave of radix-2 butterflies: */\
		"vsubpd 0x40(%%r13),%%zmm1,%%zmm1		\n\t	 vaddpd (%%r13),%%zmm7 ,%%zmm7 	\n\t"/* B1 = C1 + I*S1 */\
		"vsubpd 0x40(%%r12),%%zmm2,%%zmm2		\n\t	 vaddpd (%%r12),%%zmm8 ,%%zmm8 	\n\t"/* B2 = C2 + I*S2 */\
		"vsubpd 0x40(%%r11),%%zmm3,%%zmm3		\n\t	 vaddpd (%%r11),%%zmm9 ,%%zmm9 	\n\t"/* B3 = C3 + I*S3 */\
		"vsubpd 0x40(%%r10),%%zmm4,%%zmm4		\n\t	 vaddpd %%zmm13,%%zmm10,%%zmm10	\n\t"/* B4 = C4 + I*S4 */\
		"vsubpd 0x40(%%r9 ),%%zmm5,%%zmm5		\n\t	 vaddpd %%zmm14,%%zmm11,%%zmm11	\n\t"/* B5 = C5 + I*S5 */\
		"vsubpd 0x40(%%r8 ),%%zmm6,%%zmm6		\n\t	 vaddpd %%zmm15,%%zmm12,%%zmm12	\n\t"/* B6 = C6 + I*S6 */\
	/* Write lower-half outputs back to memory... */\
		"movq2dq	%%mm2,%%xmm0	\n\t	vbroadcastsd	%%xmm0,%%zmm0	\n\t"/* reload 2.0 from MMX0 */\
		"vmovaps	%%zmm1,(%%rax)		\n\t"/* B1r */"		vmovaps	%%zmm7 ,0x40(%%rax)	\n\t"/* B1i */\
		"vmovaps	%%zmm2,(%%rbx)		\n\t"/* B2r */"		vmovaps	%%zmm8 ,0x40(%%rbx)	\n\t"/* B2i */\
		"vmovaps	%%zmm3,(%%rcx)		\n\t"/* B3r */"		vmovaps	%%zmm9 ,0x40(%%rcx)	\n\t"/* B3i */\
		"vmovaps	%%zmm4,(%%rdx)		\n\t"/* B4r */"		vmovaps	%%zmm10,0x40(%%rdx)	\n\t"/* B4i */\
		"vmovaps	%%zmm5,(%%rsi)		\n\t"/* B5r */"		vmovaps	%%zmm11,0x40(%%rsi)	\n\t"/* B5i */\
		"vmovaps	%%zmm6,(%%rdi)		\n\t"/* B6r */"		vmovaps	%%zmm12,0x40(%%rdi)	\n\t"/* B6i */\
	/* ...Do the upper-half-index portion of the radix-2 butterflies: */\
		" vfmadd231pd 0x40(%%r13),%%zmm0,%%zmm1		\n\t	vfnmadd231pd (%%r13),%%zmm0,%%zmm7 	\n\t"/* Bc = C1 - I*S1 */\
		" vfmadd231pd 0x40(%%r12),%%zmm0,%%zmm2		\n\t	vfnmadd231pd (%%r12),%%zmm0,%%zmm8 	\n\t"/* Bb = C2 - I*S2 */\
		" vfmadd231pd 0x40(%%r11),%%zmm0,%%zmm3		\n\t	vfnmadd231pd (%%r11),%%zmm0,%%zmm9 	\n\t"/* Ba = C3 - I*S3 */\
		" vfmadd231pd 0x40(%%r10),%%zmm0,%%zmm4		\n\t	vfnmadd231pd %%zmm13,%%zmm0,%%zmm10	\n\t"/* B9 = C4 - I*S4 */\
		" vfmadd231pd 0x40(%%r9 ),%%zmm0,%%zmm5		\n\t	vfnmadd231pd %%zmm14,%%zmm0,%%zmm11	\n\t"/* B8 = C5 - I*S5 */\
		" vfmadd231pd 0x40(%%r8 ),%%zmm0,%%zmm6		\n\t	vfnmadd231pd %%zmm15,%%zmm0,%%zmm12	\n\t"/* B7 = C6 - I*S6 */\
	/* ...And write the upper-half outputs back to memory. */\
		"vmovaps	%%zmm1,(%%r13)		\n\t"/* Bcr */"		vmovaps	%%zmm7 ,0x40(%%r13)	\n\t"/* Bci */\
		"vmovaps	%%zmm2,(%%r12)		\n\t"/* Bbr */"		vmovaps	%%zmm8 ,0x40(%%r12)	\n\t"/* Bbi */\
		"vmovaps	%%zmm3,(%%r11)		\n\t"/* Bar */"		vmovaps	%%zmm9 ,0x40(%%r11)	\n\t"/* Bai */\
		"vmovaps	%%zmm4,(%%r10)		\n\t"/* B9r */"		vmovaps	%%zmm10,0x40(%%r10)	\n\t"/* B9i */\
		"vmovaps	%%zmm5,(%%r9 )		\n\t"/* B8r */"		vmovaps	%%zmm11,0x40(%%r9 )	\n\t"/* B8i */\
		"vmovaps	%%zmm6,(%%r8 )		\n\t"/* B7r */"		vmovaps	%%zmm12,0x40(%%r8 )	\n\t"/* B7i */\
		:					/* outputs: none */\
		: [in0] "m" (Xin0)	/* All inputs from memory addresses here */\
		 ,[i_off] "m" (Xi_off)	/* I-address pointer-strides */\
		 ,[__cc] "m" (Xcc)\
		 ,[out0] "m" (Xout0)\
		 ,[o_off] "m" (Xo_off)	/* O-address pointer-strides */\
		: "cc","memory","mm0","mm1","mm2",\
			"rax","rbx","rcx","rdx","rsi","rdi","r8","r9","r10","r11","r12","r13","r14","r15",\
			"xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15",\
			"xmm16","xmm17","xmm18","xmm19"	/* Clobbered registers */\
		);\
	}

   #else	// [2] DFT_13_FMA = False: [LOACC, no-flags default] based on FMAized version of same van Buskirk tan-DFT used in no-FMA mode.
	// Arithmetic opcount: [120 ADD, 78 FMA (incl 34 MUL), 168 memref], under our general target of 1 memref per vec_dbl arithmetic op.

	#define SSE2_RADIX_13_DFT(Xin0,Xi_off, Xcc, Xout0,Xo_off)\
	{\
	__asm__ volatile (\
		"movq	%[__cc],%%rsi	\n\t"\
	/* rax,rcx hold [In|Out]put base-pointers: */\
	"movq	%[in0]  ,%%rax		\n\t	movq	%[out0]  ,%%rcx		\n\t"\
	/* rbx,rdx hold pointers to uint32[] [In|Out]put pointer-offset-stride arrays: */\
	"movq	%[i_off],%%rbx		\n\t	movq	%[o_off],%%rdx		\n\t"\
	/* Compute i1-6 from offsets:		I0,O0 get stashed in r8,r9: */\
	"movslq	0x04(%%rbx),%%r15	\n\t	movslq	(%%rbx),%%r8	\n\t"\
	"movslq	0x08(%%rbx),%%r14	\n\t	movslq	(%%rdx),%%r9	\n\t"\
	"movslq	0x0c(%%rbx),%%r13	\n\t	addq	%%rax,%%r8		\n\t"/* I0 */\
	"movslq	0x10(%%rbx),%%r12	\n\t	addq	%%rcx,%%r9		\n\t"/* O0 */\
	"movslq	0x14(%%rbx),%%r11	\n\t"\
	"movslq	0x18(%%rbx),%%r10	\n\t"\
	/* xr-terms:                                                      yi-terms: */\
		"addq	%%rax,%%r15	\n\t	vmovaps	(%%r15),%%zmm1			\n\t	vmovaps	0x40(%%r15),%%zmm9 	\n\t"\
		"addq	%%rax,%%r14	\n\t	vmovaps	(%%r14),%%zmm3			\n\t	vmovaps	0x40(%%r14),%%zmm10	\n\t"\
		"addq	%%rax,%%r13	\n\t	vmovaps	(%%r13),%%zmm6			\n\t	vmovaps	0x40(%%r13),%%zmm13	\n\t"\
		"addq	%%rax,%%r12	\n\t	vmovaps	(%%r12),%%zmm4			\n\t	vmovaps	0x40(%%r12),%%zmm11	\n\t"\
		"addq	%%rax,%%r11	\n\t	vmovaps	(%%r11),%%zmm5			\n\t	vmovaps	0x40(%%r11),%%zmm12	\n\t"\
		"addq	%%rax,%%r10	\n\t	vmovaps	(%%r10),%%zmm7			\n\t	vmovaps	0x40(%%r10),%%zmm14	\n\t"\
	/* Compute iC-7 from offsets: */\
	"movslq	0x30(%%rbx),%%r15	\n\t"\
	"movslq	0x2c(%%rbx),%%r14	\n\t"\
	"movslq	0x28(%%rbx),%%r13	\n\t"\
	"movslq	0x24(%%rbx),%%r12	\n\t"\
	"movslq	0x20(%%rbx),%%r11	\n\t"\
	"movslq	0x1c(%%rbx),%%r10	\n\t"\
		"addq	%%rax,%%r15	\n\t	vaddpd	(%%r15),%%zmm1,%%zmm1	\n\t	vsubpd	0x40(%%r15),%%zmm9 ,%%zmm9 	\n\t"\
		"addq	%%rax,%%r14	\n\t	vaddpd	(%%r14),%%zmm3,%%zmm3	\n\t	vsubpd	0x40(%%r14),%%zmm10,%%zmm10	\n\t"\
		"addq	%%rax,%%r13	\n\t	vaddpd	(%%r13),%%zmm6,%%zmm6	\n\t	vsubpd	0x40(%%r13),%%zmm13,%%zmm13	\n\t"\
		"addq	%%rax,%%r12	\n\t	vaddpd	(%%r12),%%zmm4,%%zmm4	\n\t	vsubpd	0x40(%%r12),%%zmm11,%%zmm11	\n\t"\
		"addq	%%rax,%%r11	\n\t	vaddpd	(%%r11),%%zmm5,%%zmm5	\n\t	vsubpd	0x40(%%r11),%%zmm12,%%zmm12	\n\t"\
		"addq	%%rax,%%r10	\n\t	vaddpd	(%%r10),%%zmm7,%%zmm7	\n\t	vsubpd	0x40(%%r10),%%zmm14,%%zmm14	\n\t"\
	/*** Done with i1-C for now ***/\
		/* Next part (up to 'lcol section' comment) identical for Re,Im-part sections: */\
	"vmovaps	-0x40(%%rsi),%%zmm0		\n\t"/* two */\
		"														vmovaps	%%zmm9 ,%%zmm15	\n\t	vmovaps	%%zmm10,%%zmm8	\n\t"\
		"vsubpd	%%zmm5,%%zmm1,%%zmm1					\n\t			vsubpd	%%zmm11,%%zmm15,%%zmm15			\n\t"\
		"vsubpd	%%zmm6,%%zmm3,%%zmm3					\n\t			vaddpd	%%zmm12,%%zmm8 ,%%zmm8 			\n\t"\
		"vsubpd	%%zmm7,%%zmm4,%%zmm4					\n\t			vaddpd	%%zmm11,%%zmm9 ,%%zmm9 			\n\t"\
	"vfmadd132pd	%%zmm0,%%zmm1,%%zmm5				\n\t			vaddpd	%%zmm13,%%zmm15,%%zmm15			\n\t"\
	"vfmadd132pd	%%zmm0,%%zmm3,%%zmm6				\n\t			vaddpd	%%zmm11,%%zmm13,%%zmm13			\n\t"\
	"vfmadd132pd	%%zmm0,%%zmm4,%%zmm7				\n\t			vaddpd	%%zmm14,%%zmm8 ,%%zmm8 			\n\t"\
		"vmovaps	%%zmm5,%%zmm2						\n\t			vsubpd	%%zmm12,%%zmm10,%%zmm10			\n\t"\
		"vaddpd	%%zmm6,%%zmm2,%%zmm2					\n\t			vsubpd	%%zmm12,%%zmm14,%%zmm14			\n\t"\
		"vaddpd	%%zmm7,%%zmm2,%%zmm2					\n\t			vmovaps	%%zmm8 ,%%zmm12					\n\t"\
/*I0:*/	"vaddpd	(%%r8 ),%%zmm2,%%zmm0					\n\t			vmovaps	%%zmm15,%%zmm11					\n\t"\
/*O0:*/	"vmovaps	%%zmm0,(%%r9 )						\n\t			vmulpd	0x100(%%rsi),%%zmm15,%%zmm15	\n\t"\
	"vfmadd231pd	(%%rsi),%%zmm2,%%zmm0				\n\t			vmulpd	0x100(%%rsi),%%zmm8 ,%%zmm8 	\n\t"\
		"												\n\t		vfmsub132pd	0xc0(%%rsi),%%zmm15,%%zmm12	\n\t"\
		"vmovaps	0x40(%%rsi),%%zmm2					\n\t		vfmadd132pd	0xc0(%%rsi),%%zmm8 ,%%zmm11	\n\t"\
	/* lcol|rcol o-offsets 1-3|B-C - compute in order used below, store into r10-12|13-14: */\
	"movslq	0x30(%%rdx),%%r15	\n\t	addq	%%rcx,%%r15	\n\t"/* oC */\
	"movslq	0x04(%%rdx),%%r10	\n\t	addq	%%rcx,%%r10	\n\t"/* o1 */\
	"movslq	0x08(%%rdx),%%r11	\n\t	addq	%%rcx,%%r11	\n\t"/* o2 */\
	"movslq	0x0c(%%rdx),%%r12	\n\t	addq	%%rcx,%%r12	\n\t"/* o3 */\
	"movslq	0x2c(%%rdx),%%r14	\n\t	addq	%%rcx,%%r14	\n\t"/* oB */\
		"vmulpd	%%zmm2,%%zmm6,%%zmm6					\n\t			vmovaps	%%zmm12,(%%r15)			\n\t"\
		"vmulpd	%%zmm2,%%zmm5,%%zmm5					\n\t			vmovaps	%%zmm9 ,%%zmm8 			\n\t"\
		"vmulpd	%%zmm2,%%zmm7,%%zmm7					\n\t			vmovaps	%%zmm13,%%zmm12			\n\t"\
		"vmovaps %%zmm6,(%%r10)							\n\t			vaddpd	%%zmm10,%%zmm8 ,%%zmm8 			\n\t"\
		"vmovaps %%zmm5,(%%r11)							\n\t			vaddpd	%%zmm14,%%zmm12,%%zmm12			\n\t"\
		"vmovaps %%zmm7,(%%r12)							\n\t"\
		"vmovaps 0x200(%%rsi),%%zmm2					\n\t			vmovaps	%%zmm12,%%zmm15			\n\t"\
	"vfmadd132pd	%%zmm2,%%zmm0,%%zmm5				\n\t		vfmadd132pd 0x2c0(%%rsi),%%zmm8 ,%%zmm12	\n\t"\
	"vfmadd132pd	%%zmm2,%%zmm0,%%zmm7				\n\t			vmulpd 0x380(%%rsi),%%zmm8 ,%%zmm8 		\n\t"\
	"vfmadd132pd	%%zmm2,%%zmm0,%%zmm6				\n\t"\
		"												\n\t			vaddpd	%%zmm15,%%zmm8 ,%%zmm8 			\n\t"\
		"												\n\t			vmulpd	0x140(%%rsi),%%zmm12,%%zmm12	\n\t"\
		"												\n\t			vmulpd	0x140(%%rsi),%%zmm8 ,%%zmm8 		\n\t"\
		"vaddpd	(%%r10),%%zmm5,%%zmm5					\n\t			vmovaps	%%zmm10,(%%r14)			\n\t"\
		"vaddpd	(%%r11),%%zmm7,%%zmm7					\n\t			vmovaps	%%zmm14,%%zmm15			\n\t"\
		"vaddpd	(%%r12),%%zmm6,%%zmm6					\n\t			vmulpd 0x300(%%rsi),%%zmm14,%%zmm14	\n\t"\
		"vmovaps 0x280(%%rsi),%%zmm2					\n\t			vmulpd 0x3c0(%%rsi),%%zmm10,%%zmm10	\n\t"\
		"vmovaps 0x240(%%rsi),%%zmm0					\n\t			vaddpd	(%%r14),%%zmm14,%%zmm14			\n\t"\
		"vmovaps %%zmm4,(%%r10)							\n\t			vaddpd	%%zmm15,%%zmm10,%%zmm10			\n\t"\
		"vmovaps %%zmm1,(%%r11)							\n\t			vmulpd	0x180(%%rsi),%%zmm14,%%zmm14	\n\t"\
		"vmovaps %%zmm3,(%%r12)							\n\t			vmulpd	0x180(%%rsi),%%zmm10,%%zmm10	\n\t"\
		"												\n\t			vaddpd	%%zmm12,%%zmm14,%%zmm14			\n\t"\
		"												\n\t			vaddpd	%%zmm8 ,%%zmm10,%%zmm10			\n\t"\
		"												\n\t			vmovaps	%%zmm9 ,(%%r14)			\n\t"\
	"vfmadd213pd	(%%r12),%%zmm2,%%zmm4				\n\t			vmovaps	%%zmm13,%%zmm15			\n\t"\
	"vfmsub213pd	(%%r10),%%zmm2,%%zmm1				\n\t			vmulpd 0x340(%%rsi),%%zmm13,%%zmm13	\n\t"\
	"vfmadd213pd	(%%r11),%%zmm2,%%zmm3				\n\t			vmulpd 0x400(%%rsi),%%zmm9 ,%%zmm9 		\n\t"\
		"												\n\t			vaddpd	(%%r14),%%zmm13,%%zmm13			\n\t"\
		"												\n\t			vaddpd	%%zmm15,%%zmm9 ,%%zmm9 			\n\t"\
		"												\n\t			vmulpd	0x1c0(%%rsi),%%zmm13,%%zmm13	\n\t"\
		"vmovaps 0x80(%%rsi),%%zmm2						\n\t			vmulpd	0x1c0(%%rsi),%%zmm9 ,%%zmm9 		\n\t"\
	"vfmadd213pd	(%%r11),%%zmm0,%%zmm4				\n\t			vaddpd	%%zmm12,%%zmm13,%%zmm13			\n\t"\
	"vfmsub213pd	(%%r12),%%zmm0,%%zmm1				\n\t			vaddpd	%%zmm8 ,%%zmm9 ,%%zmm9 			\n\t"\
	"vfmsub213pd	(%%r10),%%zmm0,%%zmm3				\n\t			vmovaps	(%%r15),%%zmm12			\n\t"\
	"vmulpd	-0x40(%%rsi),%%zmm2,%%zmm0					\n\t			vmovaps	%%zmm11,%%zmm15			\n\t"\
	"vfnmadd231pd	%%zmm4,%%zmm2,%%zmm6				\n\t			vmovaps	%%zmm12,%%zmm8 			\n\t"\
	"vfnmadd231pd	%%zmm1,%%zmm2,%%zmm7				\n\t			vaddpd	%%zmm14,%%zmm15,%%zmm15			\n\t"\
	"vfnmadd231pd	%%zmm3,%%zmm2,%%zmm5				\n\t			vaddpd	%%zmm13,%%zmm8 ,%%zmm8 			\n\t"\
		"												\n\t			vsubpd	%%zmm11,%%zmm14,%%zmm14			\n\t"\
		"												\n\t			vsubpd	%%zmm12,%%zmm13,%%zmm13			\n\t"\
		"												\n\t			vaddpd	%%zmm10,%%zmm14,%%zmm14			\n\t"\
	"vfmadd213pd	%%zmm6,%%zmm0,%%zmm4				\n\t			vaddpd	%%zmm9 ,%%zmm13,%%zmm13			\n\t"\
	"vfmadd213pd	%%zmm7,%%zmm0,%%zmm1				\n\t			vaddpd	%%zmm11,%%zmm10,%%zmm10			\n\t"\
	"vfmadd213pd	%%zmm5,%%zmm0,%%zmm3				\n\t			vaddpd	%%zmm9 ,%%zmm12,%%zmm12			\n\t"\
	/* yi-data in zmm8,10,12,13,14,15; zmm0,2,9,11 free: */\
/*** o1-3|B,C already in r10-12|14,15: ***/\
	/*** o3 already in r12; clobber rsi for now, leaving rsi,rdi,r13 for o6,7,A: ***/\
	"movslq	0x18(%%rdx),%%rsi	\n\t	addq	%%rcx,%%rsi	\n\t"/* o6 */\
	"movslq	0x1c(%%rdx),%%rdi	\n\t	addq	%%rcx,%%rdi	\n\t"/* o7 */\
	"movslq	0x28(%%rdx),%%r13	\n\t	addq	%%rcx,%%r13	\n\t"/* oA */\
		"vmovaps	%%zmm4 ,%%zmm0						\n\t			vmovaps	%%zmm7 ,%%zmm9 			\n\t"\
		"vsubpd	%%zmm15,%%zmm4,%%zmm4					\n\t			vsubpd	%%zmm8 ,%%zmm7,%%zmm7	\n\t"\
		"vaddpd	%%zmm15,%%zmm0,%%zmm0					\n\t			vaddpd	%%zmm8 ,%%zmm9,%%zmm9	\n\t"\
		"vmovaps	%%zmm4 ,(%%rsi)						\n\t			vmovaps	%%zmm7 ,(%%r12)	\n\t"\
		"vmovaps	%%zmm0 ,(%%rdi)						\n\t			vmovaps	%%zmm9 ,(%%r13)	\n\t"\
	/* Compute o8,5,4,9, store into rsi,rdi,r12-13: */\
	"movslq	0x20(%%rdx),%%rsi	\n\t	addq	%%rcx,%%rsi	\n\t"/* o8 */\
	"movslq	0x14(%%rdx),%%rdi	\n\t	addq	%%rcx,%%rdi	\n\t"/* o5 */\
	"movslq	0x10(%%rdx),%%r12	\n\t	addq	%%rcx,%%r12	\n\t"/* o4 */\
	"movslq	0x24(%%rdx),%%r13	\n\t	addq	%%rcx,%%r13	\n\t"/* o9 */\
		"vmovaps	%%zmm5 ,%%zmm0						\n\t			vmovaps	%%zmm6 ,%%zmm9 			\n\t"\
		"vsubpd	%%zmm14,%%zmm5,%%zmm5					\n\t			vsubpd	%%zmm13,%%zmm6,%%zmm6	\n\t"\
		"vaddpd	%%zmm14,%%zmm0,%%zmm0					\n\t			vaddpd	%%zmm13,%%zmm9,%%zmm9	\n\t"\
		"vmovaps	%%zmm5 ,(%%rsi)						\n\t			vmovaps	%%zmm6 ,(%%r12)	\n\t"\
		"vmovaps	%%zmm0 ,(%%rdi)						\n\t			vmovaps	%%zmm9 ,(%%r13)	\n\t"\
	/*** o2,B,1,C already in r11,r14,r10,r15: ***/\
		"vmovaps	%%zmm1 ,%%zmm0						\n\t			vmovaps	%%zmm3 ,%%zmm9 			\n\t"\
		"vsubpd	%%zmm10,%%zmm1,%%zmm1					\n\t			vsubpd	%%zmm12,%%zmm3,%%zmm3	\n\t"\
		"vaddpd	%%zmm10,%%zmm0,%%zmm0					\n\t			vaddpd	%%zmm12,%%zmm9,%%zmm9	\n\t"\
		"vmovaps	%%zmm1 ,(%%r11)						\n\t			vmovaps	%%zmm3 ,(%%r10)	\n\t"\
		"vmovaps	%%zmm0 ,(%%r14)						\n\t			vmovaps	%%zmm9 ,(%%r15)	\n\t"\
		/****************************************************************/\
		/*                      IMAG PARTS:                             */\
		/****************************************************************/\
		"movq	%[__cc],%%rsi	\n\t"/* Reload cc-ptr, which we clobbered in above closing-butterfly section */\
		"addq	$0x40,%%r8		\n\t"/* I0.imag */\
		"addq	$0x40,%%r9		\n\t"/* O0.imag */\
	/* Compute i1-6 from offsets: */\
	"movslq	0x04(%%rbx),%%r15	\n\t"\
	"movslq	0x08(%%rbx),%%r14	\n\t"\
	"movslq	0x0c(%%rbx),%%r13	\n\t"\
	"movslq	0x10(%%rbx),%%r12	\n\t"\
	"movslq	0x14(%%rbx),%%r11	\n\t"\
	"movslq	0x18(%%rbx),%%r10	\n\t"\
	/* xi-terms:                                                      yr-terms: */\
		"addq	%%rax,%%r15	\n\t	vmovaps	0x40(%%r15),%%zmm1			\n\t	vmovaps	(%%r15),%%zmm9 	\n\t"\
		"addq	%%rax,%%r14	\n\t	vmovaps	0x40(%%r14),%%zmm3			\n\t	vmovaps	(%%r14),%%zmm10	\n\t"\
		"addq	%%rax,%%r13	\n\t	vmovaps	0x40(%%r13),%%zmm6			\n\t	vmovaps	(%%r13),%%zmm13	\n\t"\
		"addq	%%rax,%%r12	\n\t	vmovaps	0x40(%%r12),%%zmm4			\n\t	vmovaps	(%%r12),%%zmm11	\n\t"\
		"addq	%%rax,%%r11	\n\t	vmovaps	0x40(%%r11),%%zmm5			\n\t	vmovaps	(%%r11),%%zmm12	\n\t"\
		"addq	%%rax,%%r10	\n\t	vmovaps	0x40(%%r10),%%zmm7			\n\t	vmovaps	(%%r10),%%zmm14	\n\t"\
	/* Compute iC-7 from offsets: */\
	"movslq	0x30(%%rbx),%%r15	\n\t"\
	"movslq	0x2c(%%rbx),%%r14	\n\t"\
	"movslq	0x28(%%rbx),%%r13	\n\t"\
	"movslq	0x24(%%rbx),%%r12	\n\t"\
	"movslq	0x20(%%rbx),%%r11	\n\t"\
	"movslq	0x1c(%%rbx),%%r10	\n\t"\
		"addq	%%rax,%%r15	\n\t	vaddpd	0x40(%%r15),%%zmm1,%%zmm1	\n\t	vsubpd	(%%r15),%%zmm9 ,%%zmm9 	\n\t"\
		"addq	%%rax,%%r14	\n\t	vaddpd	0x40(%%r14),%%zmm3,%%zmm3	\n\t	vsubpd	(%%r14),%%zmm10,%%zmm10	\n\t"\
		"addq	%%rax,%%r13	\n\t	vaddpd	0x40(%%r13),%%zmm6,%%zmm6	\n\t	vsubpd	(%%r13),%%zmm13,%%zmm13	\n\t"\
		"addq	%%rax,%%r12	\n\t	vaddpd	0x40(%%r12),%%zmm4,%%zmm4	\n\t	vsubpd	(%%r12),%%zmm11,%%zmm11	\n\t"\
		"addq	%%rax,%%r11	\n\t	vaddpd	0x40(%%r11),%%zmm5,%%zmm5	\n\t	vsubpd	(%%r11),%%zmm12,%%zmm12	\n\t"\
		"addq	%%rax,%%r10	\n\t	vaddpd	0x40(%%r10),%%zmm7,%%zmm7	\n\t	vsubpd	(%%r10),%%zmm14,%%zmm14	\n\t"\
/*** Done w/i1-C - Still need 5 GPRs [rcx,rdx,rsi,r8,r9], to hold out0,o_off,__c,I0,O0 ***/\
		/* Next part (up to 'lcol section' comment) identical for Re,Im-part sections: */\
	"vmovaps	-0x40(%%rsi),%%zmm0		\n\t"/* two */\
		"														vmovaps	%%zmm9 ,%%zmm15	\n\t	vmovaps	%%zmm10,%%zmm8	\n\t"\
		"vsubpd	%%zmm5,%%zmm1,%%zmm1					\n\t			vsubpd	%%zmm11,%%zmm15,%%zmm15			\n\t"\
		"vsubpd	%%zmm6,%%zmm3,%%zmm3					\n\t			vaddpd	%%zmm12,%%zmm8 ,%%zmm8 			\n\t"\
		"vsubpd	%%zmm7,%%zmm4,%%zmm4					\n\t			vaddpd	%%zmm11,%%zmm9 ,%%zmm9 			\n\t"\
	"vfmadd132pd	%%zmm0,%%zmm1,%%zmm5				\n\t			vaddpd	%%zmm13,%%zmm15,%%zmm15			\n\t"\
	"vfmadd132pd	%%zmm0,%%zmm3,%%zmm6				\n\t			vaddpd	%%zmm11,%%zmm13,%%zmm13			\n\t"\
	"vfmadd132pd	%%zmm0,%%zmm4,%%zmm7				\n\t			vaddpd	%%zmm14,%%zmm8 ,%%zmm8 			\n\t"\
		"vmovaps	%%zmm5,%%zmm2						\n\t			vsubpd	%%zmm12,%%zmm10,%%zmm10			\n\t"\
		"vaddpd	%%zmm6,%%zmm2,%%zmm2					\n\t			vsubpd	%%zmm12,%%zmm14,%%zmm14			\n\t"\
		"vaddpd	%%zmm7,%%zmm2,%%zmm2					\n\t			vmovaps	%%zmm8 ,%%zmm12					\n\t"\
/*I0:*/	"vaddpd	(%%r8 ),%%zmm2,%%zmm0					\n\t			vmovaps	%%zmm15,%%zmm11					\n\t"\
/*O0:*/	"vmovaps	%%zmm0,(%%r9 )						\n\t			vmulpd	0x100(%%rsi),%%zmm15,%%zmm15	\n\t"\
/*** Done w/I0,O0 - At this point could use any 12 GPRs of 13 available [14 minus rsi] to hold o1-C, if desired ***/\
	"vfmadd231pd	(%%rsi),%%zmm2,%%zmm0				\n\t			vmulpd	0x100(%%rsi),%%zmm8 ,%%zmm8 	\n\t"\
		"												\n\t		vfmsub132pd	0xc0(%%rsi),%%zmm15,%%zmm12	\n\t"\
		"vmovaps	0x40(%%rsi),%%zmm2					\n\t		vfmadd132pd	0xc0(%%rsi),%%zmm8 ,%%zmm11	\n\t"\
	/* lcol|rcol o-offsets A-C|1-2 - compute in order used below, store into r13-15|10-11: */\
	"movslq	0x04(%%rdx),%%r10	\n\t	addq	%%rcx,%%r10	\n\t"/* o1 */\
	"movslq	0x30(%%rdx),%%r15	\n\t	addq	%%rcx,%%r15	\n\t"/* oC */\
	"movslq	0x2c(%%rdx),%%r14	\n\t	addq	%%rcx,%%r14	\n\t"/* oB */\
	"movslq	0x28(%%rdx),%%r13	\n\t	addq	%%rcx,%%r13	\n\t"/* oA */\
	"movslq	0x08(%%rdx),%%r11	\n\t	addq	%%rcx,%%r11	\n\t"/* o2 */\
		/* Next part (up to 'yi-data' comment) identical for Re,Im, except for replacements (%%r10-15) ==> 0x40(%%r15-10): */\
		"vmulpd	%%zmm2,%%zmm6,%%zmm6					\n\t			vmovaps	%%zmm12,0x40(%%r10)			\n\t"\
		"vmulpd	%%zmm2,%%zmm5,%%zmm5					\n\t			vmovaps	%%zmm9 ,%%zmm8 			\n\t"\
		"vmulpd	%%zmm2,%%zmm7,%%zmm7					\n\t			vmovaps	%%zmm13,%%zmm12			\n\t"\
		"vmovaps %%zmm6,0x40(%%r15)						\n\t			vaddpd	%%zmm10,%%zmm8 ,%%zmm8 			\n\t"\
		"vmovaps %%zmm5,0x40(%%r14)						\n\t			vaddpd	%%zmm14,%%zmm12,%%zmm12			\n\t"\
		"vmovaps %%zmm7,0x40(%%r13)						\n\t"\
		"vmovaps 0x200(%%rsi),%%zmm2					\n\t			vmovaps	%%zmm12,%%zmm15			\n\t"\
	"vfmadd132pd	%%zmm2,%%zmm0,%%zmm5				\n\t		vfmadd132pd 0x2c0(%%rsi),%%zmm8 ,%%zmm12	\n\t"\
	"vfmadd132pd	%%zmm2,%%zmm0,%%zmm7				\n\t			vmulpd 0x380(%%rsi),%%zmm8 ,%%zmm8 		\n\t"\
	"vfmadd132pd	%%zmm2,%%zmm0,%%zmm6				\n\t"\
		"												\n\t			vaddpd	%%zmm15,%%zmm8 ,%%zmm8 			\n\t"\
		"												\n\t			vmulpd	0x140(%%rsi),%%zmm12,%%zmm12	\n\t"\
		"												\n\t			vmulpd	0x140(%%rsi),%%zmm8 ,%%zmm8 		\n\t"\
		"vaddpd	0x40(%%r15),%%zmm5,%%zmm5				\n\t			vmovaps	%%zmm10,0x40(%%r11)			\n\t"\
		"vaddpd	0x40(%%r14),%%zmm7,%%zmm7				\n\t			vmovaps	%%zmm14,%%zmm15			\n\t"\
		"vaddpd	0x40(%%r13),%%zmm6,%%zmm6				\n\t			vmulpd 0x300(%%rsi),%%zmm14,%%zmm14	\n\t"\
		"vmovaps 0x280(%%rsi),%%zmm2					\n\t			vmulpd 0x3c0(%%rsi),%%zmm10,%%zmm10	\n\t"\
		"vmovaps 0x240(%%rsi),%%zmm0					\n\t			vaddpd	0x40(%%r11),%%zmm14,%%zmm14			\n\t"\
		"vmovaps %%zmm4,0x40(%%r15)						\n\t			vaddpd	%%zmm15,%%zmm10,%%zmm10			\n\t"\
		"vmovaps %%zmm1,0x40(%%r14)						\n\t			vmulpd	0x180(%%rsi),%%zmm14,%%zmm14	\n\t"\
		"vmovaps %%zmm3,0x40(%%r13)						\n\t			vmulpd	0x180(%%rsi),%%zmm10,%%zmm10	\n\t"\
		"												\n\t			vaddpd	%%zmm12,%%zmm14,%%zmm14			\n\t"\
		"												\n\t			vaddpd	%%zmm8 ,%%zmm10,%%zmm10			\n\t"\
		"												\n\t			vmovaps	%%zmm9 ,0x40(%%r11)			\n\t"\
	"vfmadd213pd	0x40(%%r13),%%zmm2,%%zmm4			\n\t			vmovaps	%%zmm13,%%zmm15			\n\t"\
	"vfmsub213pd	0x40(%%r15),%%zmm2,%%zmm1			\n\t			vmulpd 0x340(%%rsi),%%zmm13,%%zmm13	\n\t"\
	"vfmadd213pd	0x40(%%r14),%%zmm2,%%zmm3			\n\t			vmulpd 0x400(%%rsi),%%zmm9 ,%%zmm9 		\n\t"\
		"												\n\t			vaddpd	0x40(%%r11),%%zmm13,%%zmm13			\n\t"\
		"												\n\t			vaddpd	%%zmm15,%%zmm9 ,%%zmm9 			\n\t"\
		"												\n\t			vmulpd	0x1c0(%%rsi),%%zmm13,%%zmm13	\n\t"\
		"vmovaps 0x80(%%rsi),%%zmm2						\n\t			vmulpd	0x1c0(%%rsi),%%zmm9 ,%%zmm9 		\n\t"\
	"vfmadd213pd	0x40(%%r14),%%zmm0,%%zmm4			\n\t			vaddpd	%%zmm12,%%zmm13,%%zmm13			\n\t"\
	"vfmsub213pd	0x40(%%r13),%%zmm0,%%zmm1			\n\t			vaddpd	%%zmm8 ,%%zmm9 ,%%zmm9 			\n\t"\
	"vfmsub213pd	0x40(%%r15),%%zmm0,%%zmm3			\n\t			vmovaps	0x40(%%r10),%%zmm12			\n\t"\
	"vmulpd	-0x40(%%rsi),%%zmm2,%%zmm0					\n\t			vmovaps	%%zmm11,%%zmm15			\n\t"\
	"vfnmadd231pd	%%zmm4,%%zmm2,%%zmm6				\n\t			vmovaps	%%zmm12,%%zmm8 			\n\t"\
	"vfnmadd231pd	%%zmm1,%%zmm2,%%zmm7				\n\t			vaddpd	%%zmm14,%%zmm15,%%zmm15			\n\t"\
	"vfnmadd231pd	%%zmm3,%%zmm2,%%zmm5				\n\t			vaddpd	%%zmm13,%%zmm8 ,%%zmm8 			\n\t"\
		"												\n\t			vsubpd	%%zmm11,%%zmm14,%%zmm14			\n\t"\
		"												\n\t			vsubpd	%%zmm12,%%zmm13,%%zmm13			\n\t"\
		"												\n\t			vaddpd	%%zmm10,%%zmm14,%%zmm14			\n\t"\
	"vfmadd213pd	%%zmm6,%%zmm0,%%zmm4				\n\t			vaddpd	%%zmm9 ,%%zmm13,%%zmm13			\n\t"\
	"vfmadd213pd	%%zmm7,%%zmm0,%%zmm1				\n\t			vaddpd	%%zmm11,%%zmm10,%%zmm10			\n\t"\
	"vfmadd213pd	%%zmm5,%%zmm0,%%zmm3				\n\t			vaddpd	%%zmm9 ,%%zmm12,%%zmm12			\n\t"\
	/* yi-data in zmm8,10,12,13,14,15; zmm0,2,9,11 free: */\
/*** o1,2|A-C already in r10,11|13-15: ***/\
	/*** oA already in r13; Compute o3,6,7 into r12,8,9: ***/\
	"movslq	0x0c(%%rdx),%%r12	\n\t	addq	%%rcx,%%r12	\n\t"/* o3 */\
	"movslq	0x18(%%rdx),%%r8	\n\t	addq	%%rcx,%%r8	\n\t"/* o6 */\
	"movslq	0x1c(%%rdx),%%r9	\n\t	addq	%%rcx,%%r9	\n\t"/* o7 */\
		"vmovaps	%%zmm4 ,%%zmm0						\n\t			vmovaps	%%zmm7 ,%%zmm9		\n\t"\
		"vsubpd	%%zmm15,%%zmm4,%%zmm4					\n\t			vsubpd	%%zmm8 ,%%zmm7,%%zmm7		\n\t"\
		"vaddpd	%%zmm15,%%zmm0,%%zmm0					\n\t			vaddpd	%%zmm8 ,%%zmm9,%%zmm9		\n\t"\
		"vmovaps	%%zmm4 ,0x40(%%r9 )					\n\t			vmovaps	%%zmm7 ,0x40(%%r13)	\n\t"\
		"vmovaps	%%zmm0 ,0x40(%%r8 )					\n\t			vmovaps	%%zmm9 ,0x40(%%r12)	\n\t"\
	/* Compute o5,8,9,4, store into rsi,rdi,r12-13: */\
	"movslq	0x14(%%rdx),%%rsi	\n\t	addq	%%rcx,%%rsi	\n\t"/* o5 */\
	"movslq	0x20(%%rdx),%%rdi	\n\t	addq	%%rcx,%%rdi	\n\t"/* o8 */\
	"movslq	0x24(%%rdx),%%r13	\n\t	addq	%%rcx,%%r13	\n\t"/* o9 */\
	"movslq	0x10(%%rdx),%%r12	\n\t	addq	%%rcx,%%r12	\n\t"/* o4 */\
		"vmovaps	%%zmm5 ,%%zmm0						\n\t			vmovaps	%%zmm6 ,%%zmm9		\n\t"\
		"vsubpd	%%zmm14,%%zmm5,%%zmm5					\n\t			vsubpd	%%zmm13,%%zmm6,%%zmm6		\n\t"\
		"vaddpd	%%zmm14,%%zmm0,%%zmm0					\n\t			vaddpd	%%zmm13,%%zmm9,%%zmm9		\n\t"\
		"vmovaps	%%zmm5 ,0x40(%%rsi)					\n\t			vmovaps	%%zmm6 ,0x40(%%r13)	\n\t"\
		"vmovaps	%%zmm0 ,0x40(%%rdi)					\n\t			vmovaps	%%zmm9 ,0x40(%%r12)	\n\t"\
	/*** o1,2,B,C already in r10,11,14,15: ***/\
		"vmovaps	%%zmm1 ,%%zmm0						\n\t			vmovaps	%%zmm3 ,%%zmm9		\n\t"\
		"vsubpd	%%zmm10,%%zmm1,%%zmm1					\n\t			vsubpd	%%zmm12,%%zmm3,%%zmm3		\n\t"\
		"vaddpd	%%zmm10,%%zmm0,%%zmm0					\n\t			vaddpd	%%zmm12,%%zmm9,%%zmm9		\n\t"\
		"vmovaps	%%zmm1 ,0x40(%%r14)					\n\t			vmovaps	%%zmm3 ,0x40(%%r15)	\n\t"\
		"vmovaps	%%zmm0 ,0x40(%%r11)					\n\t			vmovaps	%%zmm9 ,0x40(%%r10)	\n\t"\
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

   #endif	// ? DFT_13_FMA

#elif defined(USE_AVX2)	// AVX+FMA versions: there are 2 of these:
	// [1] Naive but good-ROE and FMA-friendly impl based on RADIX_13_DFT_BASIC in dft_macro.h;
	// [2] Based on FMAized version of same van Buskirk tangent-DFT used in no-FMA mode.
	// Sample all-4-Haswell-cores timings/ROEs of 10kiter at FFT length 3328K [Res64: 69BB5D04E781424E]:
	//					DFT_13_FMA											Tangent-DFT
	//	AvgMaxErr = 0.240814071. MaxErr = 0.328125000		AvgMaxErr = 0.251917766, MaxErr = 0.375000000
	//	8.92 ms/iter										8.75 ms/iter
	/*
	Jan 2021: slightly lower self-test expos in v19.1, repeat above experiment on Broadwell NUC, 2-core/4-thread:
			10000 iterations of M62705077 with FFT length 3407872 = 3328 K, Res64: AF81565673A4D22B
		Radset = 208,16,16,32:
			AvgMaxErr = 0.216766192. MaxErr = 0.312500000		AvgMaxErr = 0.225895141. MaxErr = 0.343750000
			27.6 ms/iter										26.1 ms/iter, 5% faster for whole-iteration
		Radset = 52,32,32,32:
			AvgMaxErr = 0.181361694. MaxErr = 0.250000000		AvgMaxErr = 0.189262038. MaxErr = 0.281250000
			26.4 ms/iter										25.9 ms/iter
	NUC8I3CYSM, avx-512 build, 2-core/4-thread:
		Radset = 208,16,16,32:
			11.4 ms/iter										11.1 ms/iter
		Radset = 52,32,32,32:
			13.5 ms/iter										13.3 ms/iter
	KNL, avx-512 build, -cpu 64:67 (with F30 job running on cores 0-63):
		Radset = 208,16,16,32:
			26.2 ms/iter										25.6 ms/iter
		Radset = 52,32,32,32:
			32.5 ms/iter										30.1 ms/iter
	*/
   #if DFT_13_FMA	// [1] Naive but good-ROE and FMA-friendly impl based on RADIX_13_DFT_BASIC in dft_macro.h:

	// Arithmetic opcount: [36 ADD, 168 FMA (incl 12 MUL), 268 memref], well above general target of 1 memref per vec_dbl arithmetic op.
	// Potentially lower cycle count (at a max theoretical rate of 2 FMA/cycle) than non-FMA version, but quadratic nature of core convo
	// computation really starting to bite here, less obviously a win than for radix-11.
	//
	// Compare to van-Buskirk + FMA: [120 ADD, 78 FMA (incl 34 MUL), 168 memref], only slightly less accurate, appreciably faster.
	//
	#define SSE2_RADIX_13_DFT(Xin0,Xi_off, Xcc, Xout0,Xo_off)\
	{\
	__asm__ volatile (\
		/* Jan 2021: In reduced-#args revamp of this macro, Due to GPR-usage constraints, instead of
		keeping a ptr to vec_dbl*two in GPR r15 throughout, create 2.0 in low double, of an XMM, copy
		to an MMX via MOVDQ2Q, and copy-back-via-MOVQ2DQ/shuffle-propagate-to-rest-of-XMM as needed.
		Note that in AVX-512 we can create consts directly into GPRs and use VPBROADCAST to broadcast
		to a vector register, but in AVX2 and below the register-source version of this instruction is
		not available.
		***UPDATE***: Alas, the use-MMX trick is HUGELY SLOW on my Haswell CPU, so instead use the fact
		that in DFT_13_FMA mode, high 5 slots of precomputed-trig-consts block, __c + [12-16], unused.
		*/\
		"movq	%[__cc],%%rbx				\n\t"\
		"vmovaps	-0x20(%%rbx),%%ymm0		\n\t"/* two */\
	"movq	%[in0],%%r14		\n\t"/* r14 holds Input base-pointer */\
	"movq	%[i_off],%%r15		\n\t"/* r15 holds pointer to uint32[] Outnput pointer-offset-stride array */\
	/* Stash I0-ptr in qword 0 of unused (vec_dbl*)__cc + 16 slot for later reload: */\
	"movslq	    (%%r15),%%r8 	\n\t"\
	"addq	%%r14,%%r8 	\n\t	movq	%%r8,0x200(%%rbx)	\n\t"\
	/* Compute i1-6 from offsets: */\
	"movslq	0x04(%%r15),%%rax	\n\t"\
	"movslq	0x08(%%r15),%%rbx	\n\t"\
	"movslq	0x0c(%%r15),%%rcx	\n\t"\
	"movslq	0x10(%%r15),%%rdx	\n\t"\
	"movslq	0x14(%%r15),%%rsi	\n\t"\
	"movslq	0x18(%%r15),%%rdi	\n\t"\
	"addq	%%r14,%%rax	\n\t	vmovaps	(%%rax),%%ymm1	\n\t	vmovaps	0x20(%%rax),%%ymm7 	\n\t"/* A1r,i */\
	"addq	%%r14,%%rbx	\n\t	vmovaps	(%%rbx),%%ymm2	\n\t	vmovaps	0x20(%%rbx),%%ymm8 	\n\t"/* A2r,i */\
	"addq	%%r14,%%rcx	\n\t	vmovaps	(%%rcx),%%ymm3	\n\t	vmovaps	0x20(%%rcx),%%ymm9 	\n\t"/* A3r,i */\
	"addq	%%r14,%%rdx	\n\t	vmovaps	(%%rdx),%%ymm4	\n\t	vmovaps	0x20(%%rdx),%%ymm10	\n\t"/* A4r,i */\
	"addq	%%r14,%%rsi	\n\t	vmovaps	(%%rsi),%%ymm5	\n\t	vmovaps	0x20(%%rsi),%%ymm11	\n\t"/* A5r,i */\
	"addq	%%r14,%%rdi	\n\t	vmovaps	(%%rdi),%%ymm6	\n\t	vmovaps	0x20(%%rdi),%%ymm12	\n\t"/* A6r,i */\
	/* Compute iC-7 from offsets: */\
	"movslq	0x30(%%r15),%%r13	\n\t"\
	"movslq	0x2c(%%r15),%%r12	\n\t"\
	"movslq	0x28(%%r15),%%r11	\n\t"\
	"movslq	0x24(%%r15),%%r10	\n\t"\
	"movslq	0x20(%%r15),%%r9 	\n\t"\
	"movslq	0x1c(%%r15),%%r8 	\n\t"\
	/* ...Do the + parts of the opening wave of radix-2 butterflies: */\
	"addq	%%r14,%%r13	\n\t	vaddpd (%%r13),%%ymm1,%%ymm1	\n\t	vaddpd 0x20(%%r13),%%ymm7 ,%%ymm7 	\n\t"/* t1 = A1 + Ac */\
	"addq	%%r14,%%r12	\n\t	vaddpd (%%r12),%%ymm2,%%ymm2	\n\t	vaddpd 0x20(%%r12),%%ymm8 ,%%ymm8 	\n\t"/* t2 = A2 + Ab */\
	"addq	%%r14,%%r11	\n\t	vaddpd (%%r11),%%ymm3,%%ymm3	\n\t	vaddpd 0x20(%%r11),%%ymm9 ,%%ymm9 	\n\t"/* t3 = A3 + Aa */\
	"addq	%%r14,%%r10	\n\t	vaddpd (%%r10),%%ymm4,%%ymm4	\n\t	vaddpd 0x20(%%r10),%%ymm10,%%ymm10	\n\t"/* t4 = A4 + A9 */\
	"addq	%%r14,%%r9 	\n\t	vaddpd (%%r9 ),%%ymm5,%%ymm5	\n\t	vaddpd 0x20(%%r9 ),%%ymm11,%%ymm11	\n\t"/* t5 = A5 + A8 */\
	"addq	%%r14,%%r8 	\n\t	vaddpd (%%r8 ),%%ymm6,%%ymm6	\n\t	vaddpd 0x20(%%r8 ),%%ymm12,%%ymm12	\n\t"/* t6 = A6 + A7 */\
	/* Write lower-half outputs back to memory... */\
		"vmovaps	%%ymm1,(%%rax)		\n\t"/* t1r */"		vmovaps	%%ymm7 ,0x20(%%rax)	\n\t"/* t1i */\
		"vmovaps	%%ymm2,(%%rbx)		\n\t"/* t2r */"		vmovaps	%%ymm8 ,0x20(%%rbx)	\n\t"/* t2i */\
		"vmovaps	%%ymm3,(%%rcx)		\n\t"/* t3r */"		vmovaps	%%ymm9 ,0x20(%%rcx)	\n\t"/* t3i */\
		"vmovaps	%%ymm4,(%%rdx)		\n\t"/* t4r */"		vmovaps	%%ymm10,0x20(%%rdx)	\n\t"/* t4i */\
		"vmovaps	%%ymm5,(%%rsi)		\n\t"/* t5r */"		vmovaps	%%ymm11,0x20(%%rsi)	\n\t"/* t5i */\
		"vmovaps	%%ymm6,(%%rdi)		\n\t"/* t6r */"		vmovaps	%%ymm12,0x20(%%rdi)	\n\t"/* t6i */\
	/* ...Do the - parts of the radix-2 butterflies: */\
		"vfnmadd231pd (%%r13),%%ymm0,%%ymm1		\n\t	vfnmadd231pd 0x20(%%r13),%%ymm0,%%ymm7 	\n\t"/* ta = A1 - Ac */\
		"vfnmadd231pd (%%r12),%%ymm0,%%ymm2		\n\t	vfnmadd231pd 0x20(%%r12),%%ymm0,%%ymm8 	\n\t"/* t9 = A2 - Ab */\
		"vfnmadd231pd (%%r11),%%ymm0,%%ymm3		\n\t	vfnmadd231pd 0x20(%%r11),%%ymm0,%%ymm9 	\n\t"/* t8 = A3 - Aa */\
		"vfnmadd231pd (%%r10),%%ymm0,%%ymm4		\n\t	vfnmadd231pd 0x20(%%r10),%%ymm0,%%ymm10	\n\t"/* t7 = A4 - A9 */\
		"vfnmadd231pd (%%r9 ),%%ymm0,%%ymm5		\n\t	vfnmadd231pd 0x20(%%r9 ),%%ymm0,%%ymm11	\n\t"/* t6 = A5 - A8 */\
		"vfnmadd231pd (%%r8 ),%%ymm0,%%ymm6		\n\t	vfnmadd231pd 0x20(%%r8 ),%%ymm0,%%ymm12	\n\t"/* t6 = A6 - A7 */\
	/* ...And write the upper-half outputs back to memory to free up registers for the 5x5 sine-term-subconvo computation: */\
		"vmovaps	%%ymm1,(%%r13)		\n\t"/* tar */"		vmovaps	%%ymm7 ,0x20(%%r13)	\n\t"/* tai */\
		"vmovaps	%%ymm2,(%%r12)		\n\t"/* t9r */"		vmovaps	%%ymm8 ,0x20(%%r12)	\n\t"/* t9i */\
		"vmovaps	%%ymm3,(%%r11)		\n\t"/* t8r */"		vmovaps	%%ymm9 ,0x20(%%r11)	\n\t"/* t8i */\
		"vmovaps	%%ymm4,(%%r10)		\n\t"/* t7r */"		vmovaps	%%ymm10,0x20(%%r10)	\n\t"/* t7i */\
		"vmovaps	%%ymm5,(%%r9 )		\n\t"/* t6r */"		vmovaps	%%ymm11,0x20(%%r9 )	\n\t"/* t6i */\
		"vmovaps	%%ymm6,(%%r8 )		\n\t"/* t6r */"		vmovaps	%%ymm12,0x20(%%r8 )	\n\t"/* t6i */\
/*
	S1 =      ss1*tc + ss2*tb + ss3*ta + ss4*t9 + ss5*t8 + ss6*t7
	S2 =      ss2*tc + ss4*tb + ss6*ta - ss5*t9 - ss3*t8 - ss1*t7
	S3 =      ss3*tc + ss6*tb - ss4*ta - ss1*t9 + ss2*t8 + ss5*t7
	S4 =      ss4*tc - ss5*tb - ss1*ta + ss3*t9 - ss6*t8 - ss2*t7
	S5 =      ss5*tc - ss3*tb + ss2*ta - ss6*t9 - ss1*t8 + ss4*t7
	S6 =      ss6*tc - ss1*tb + ss5*ta - ss2*t9 + ss4*t8 - ss3*t7

	In a 16-reg FMA model, can keep the 12 S-terms, both t*r,i-mults and 2 of the 6 ss-mults in-reg, thus need 8 loads-from-mem per block:
*/\
	/* Done with rax,rbx,rcx,rdx,rsi,rdi for now, use rsi to hold trig pointer: */\
	"movq	%[__cc],%%rsi	\n\t"/* cc1 */\
		"vmovaps	(%%r13),%%ymm1		\n\t	vmovaps	0x20(%%r13),%%ymm7	\n\t"/* S1 = tc */\
		"vmovaps	%%ymm1,%%ymm2		\n\t	vmovaps	%%ymm7 ,%%ymm8 		\n\t"/* S2 = tc */\
		"vmovaps	%%ymm1,%%ymm3		\n\t	vmovaps	%%ymm7 ,%%ymm9 		\n\t"/* S3 = tc */\
		"vmovaps	%%ymm1,%%ymm4		\n\t	vmovaps	%%ymm7 ,%%ymm10		\n\t"/* S4 = tc */\
		"vmovaps	%%ymm1,%%ymm5		\n\t	vmovaps	%%ymm7 ,%%ymm11		\n\t"/* S5 = tc */\
		"vmovaps	%%ymm1,%%ymm6		\n\t	vmovaps	%%ymm7 ,%%ymm12		\n\t"/* S6 = tc */\
	"addq	$0xc0,%%rsi	\n\t"/* Incr trig ptr to point to ss1 */\
		"vmovaps	    (%%rsi),%%ymm13	\n\t	vmovaps	0x20(%%rsi),%%ymm14		\n\t"/* ss1,ss2 */\
\
	"vmulpd	    %%ymm13,%%ymm1,%%ymm1	\n\t	vmulpd	    %%ymm13,%%ymm7 ,%%ymm7 	\n\t"/* S1  = ss1*tc */\
	"vmulpd	    %%ymm14,%%ymm2,%%ymm2	\n\t	vmulpd	    %%ymm14,%%ymm8 ,%%ymm8 	\n\t"/* S2  = ss2*tc */\
	"vmulpd	0x40(%%rsi),%%ymm3,%%ymm3	\n\t	vmulpd	0x40(%%rsi),%%ymm9 ,%%ymm9 	\n\t"/* S3  = ss3*tc */\
	"vmulpd	0x60(%%rsi),%%ymm4,%%ymm4	\n\t	vmulpd	0x60(%%rsi),%%ymm10,%%ymm10	\n\t"/* S4  = ss4*tc */\
	"vmulpd	0x80(%%rsi),%%ymm5,%%ymm5	\n\t	vmulpd	0x80(%%rsi),%%ymm11,%%ymm11	\n\t"/* S5  = ss5*tc */\
	"vmulpd	0xa0(%%rsi),%%ymm6,%%ymm6	\n\t	vmulpd	0xa0(%%rsi),%%ymm12,%%ymm12	\n\t"/* S6  = ss6*tc */\
\
		"vmovaps	    (%%r12),%%ymm0	\n\t	vmovaps	0x20(%%r12),%%ymm15	\n\t"/* tb */\
	" vfmadd231pd     %%ymm14,%%ymm0,%%ymm1	\n\t  vfmadd231pd     %%ymm14,%%ymm15,%%ymm7 	\n\t"/* S1 += ss2*tb */\
	" vfmadd231pd 0x60(%%rsi),%%ymm0,%%ymm2	\n\t  vfmadd231pd 0x60(%%rsi),%%ymm15,%%ymm8 	\n\t"/* S2 += ss4*tb */\
	" vfmadd231pd 0xa0(%%rsi),%%ymm0,%%ymm3	\n\t  vfmadd231pd 0xa0(%%rsi),%%ymm15,%%ymm9 	\n\t"/* S3 += ss6*tb */\
	"vfnmadd231pd 0x80(%%rsi),%%ymm0,%%ymm4	\n\t vfnmadd231pd 0x80(%%rsi),%%ymm15,%%ymm10	\n\t"/* S4 -= ss5*tb */\
	"vfnmadd231pd 0x40(%%rsi),%%ymm0,%%ymm5	\n\t vfnmadd231pd 0x40(%%rsi),%%ymm15,%%ymm11	\n\t"/* S5 -= ss3*tb */\
	"vfnmadd231pd     %%ymm13,%%ymm0,%%ymm6	\n\t vfnmadd231pd     %%ymm13,%%ymm15,%%ymm12	\n\t"/* S6 -= ss1*tb */\
\
		"vmovaps	    (%%r11),%%ymm0	\n\t	vmovaps	0x20(%%r11),%%ymm15	\n\t"/* ta */\
	" vfmadd231pd 0x40(%%rsi),%%ymm0,%%ymm1	\n\t  vfmadd231pd 0x40(%%rsi),%%ymm15,%%ymm7 	\n\t"/* S1 += ss3*ta */\
	" vfmadd231pd 0xa0(%%rsi),%%ymm0,%%ymm2	\n\t  vfmadd231pd 0xa0(%%rsi),%%ymm15,%%ymm8 	\n\t"/* S2 += ss6*ta */\
	"vfnmadd231pd 0x60(%%rsi),%%ymm0,%%ymm3	\n\t vfnmadd231pd 0x60(%%rsi),%%ymm15,%%ymm9 	\n\t"/* S3 -= ss4*ta */\
	"vfnmadd231pd     %%ymm13,%%ymm0,%%ymm4	\n\t vfnmadd231pd     %%ymm13,%%ymm15,%%ymm10	\n\t"/* S4 -= ss1*ta */\
	" vfmadd231pd     %%ymm14,%%ymm0,%%ymm5	\n\t  vfmadd231pd     %%ymm14,%%ymm15,%%ymm11	\n\t"/* S5 += ss2*ta */\
	" vfmadd231pd 0x80(%%rsi),%%ymm0,%%ymm6	\n\t  vfmadd231pd 0x80(%%rsi),%%ymm15,%%ymm12	\n\t"/* S6 += ss5*ta */\
\
		"vmovaps	    (%%r10),%%ymm0	\n\t	vmovaps	0x20(%%r10),%%ymm15	\n\t"/* t9 */\
	" vfmadd231pd 0x60(%%rsi),%%ymm0,%%ymm1	\n\t  vfmadd231pd 0x60(%%rsi),%%ymm15,%%ymm7 	\n\t"/* S1 += ss4*t9 */\
	"vfnmadd231pd 0x80(%%rsi),%%ymm0,%%ymm2	\n\t vfnmadd231pd 0x80(%%rsi),%%ymm15,%%ymm8 	\n\t"/* S2 -= ss5*t9 */\
	"vfnmadd231pd     %%ymm13,%%ymm0,%%ymm3	\n\t vfnmadd231pd     %%ymm13,%%ymm15,%%ymm9 	\n\t"/* S3 -= ss1*t9 */\
	" vfmadd231pd 0x40(%%rsi),%%ymm0,%%ymm4	\n\t  vfmadd231pd 0x40(%%rsi),%%ymm15,%%ymm10	\n\t"/* S4 += ss3*t9 */\
	"vfnmadd231pd 0xa0(%%rsi),%%ymm0,%%ymm5	\n\t vfnmadd231pd 0xa0(%%rsi),%%ymm15,%%ymm11	\n\t"/* S5 -= ss6*t9 */\
	"vfnmadd231pd     %%ymm14,%%ymm0,%%ymm6	\n\t vfnmadd231pd     %%ymm14,%%ymm15,%%ymm12	\n\t"/* S6 -= ss2*t9 */\
\
		"vmovaps	    (%%r9 ),%%ymm0	\n\t	vmovaps	0x20(%%r9 ),%%ymm15	\n\t"/* t8 */\
	" vfmadd231pd 0x80(%%rsi),%%ymm0,%%ymm1	\n\t  vfmadd231pd 0x80(%%rsi),%%ymm15,%%ymm7 	\n\t"/* S1 += ss5*t8 */\
	"vfnmadd231pd 0x40(%%rsi),%%ymm0,%%ymm2	\n\t vfnmadd231pd 0x40(%%rsi),%%ymm15,%%ymm8 	\n\t"/* S2 -= ss3*t8 */\
	" vfmadd231pd     %%ymm14,%%ymm0,%%ymm3	\n\t  vfmadd231pd     %%ymm14,%%ymm15,%%ymm9 	\n\t"/* S3 += ss2*t8 */\
	"vfnmadd231pd 0xa0(%%rsi),%%ymm0,%%ymm4	\n\t vfnmadd231pd 0xa0(%%rsi),%%ymm15,%%ymm10	\n\t"/* S4 -= ss6*t8 */\
	"vfnmadd231pd     %%ymm13,%%ymm0,%%ymm5	\n\t vfnmadd231pd     %%ymm13,%%ymm15,%%ymm11	\n\t"/* S5 -= ss1*t8 */\
	" vfmadd231pd 0x60(%%rsi),%%ymm0,%%ymm6	\n\t  vfmadd231pd 0x60(%%rsi),%%ymm15,%%ymm12	\n\t"/* S6 += ss4*t8 */\
\
		"vmovaps	    (%%r8 ),%%ymm0	\n\t	vmovaps	0x20(%%r8 ),%%ymm15	\n\t"/* t7 */\
	" vfmadd231pd 0xa0(%%rsi),%%ymm0,%%ymm1	\n\t  vfmadd231pd 0xa0(%%rsi),%%ymm15,%%ymm7 	\n\t"/* S1 += ss6*t7 */\
	"vfnmadd231pd     %%ymm13,%%ymm0,%%ymm2	\n\t vfnmadd231pd     %%ymm13,%%ymm15,%%ymm8 	\n\t"/* S2 -= ss1*t7 */\
	" vfmadd231pd 0x80(%%rsi),%%ymm0,%%ymm3	\n\t  vfmadd231pd 0x80(%%rsi),%%ymm15,%%ymm9 	\n\t"/* S3 += ss5*t7 */\
	"vfnmadd231pd     %%ymm14,%%ymm0,%%ymm4	\n\t vfnmadd231pd     %%ymm14,%%ymm15,%%ymm10	\n\t"/* S4 -= ss2*t7 */\
	" vfmadd231pd 0x60(%%rsi),%%ymm0,%%ymm5	\n\t  vfmadd231pd 0x60(%%rsi),%%ymm15,%%ymm11	\n\t"/* S5 += ss4*t7 */\
	"vfnmadd231pd 0x40(%%rsi),%%ymm0,%%ymm6	\n\t vfnmadd231pd 0x40(%%rsi),%%ymm15,%%ymm12	\n\t"/* S6 -= ss3*t7 */\
\
	/* Before overwriting r14,15 with O-side base|offset-ptrs, restore i5 and load A0 into ymm0:*/\
	"movslq	0x14(%%r15),%%rsi	\n\t	addq	%%r14,%%rsi	\n\t"\
	/* Write sine-term-subconvo outputs back to memory to free up registers for the 5x5 cosine-term-subconvo computation: */\
	/* i7-C ptrs no longer needed, so load o-ptrs in their place: */\
	"movq	%[out0],%%r14		\n\t"/* r14 holds Output base-pointer */\
	"movq	%[o_off],%%r15		\n\t"/* r15 holds pointer to uint32[] Outnput pointer-offset-stride array */\
	/* Compute oC-7 from offsets: */\
	"movslq	0x30(%%r15),%%r13	\n\t	addq	%%r14,%%r13	\n\t"\
	"movslq	0x2c(%%r15),%%r12	\n\t	addq	%%r14,%%r12	\n\t"\
	"movslq	0x28(%%r15),%%r11	\n\t	addq	%%r14,%%r11	\n\t"\
	"movslq	0x24(%%r15),%%r10	\n\t	addq	%%r14,%%r10	\n\t"\
	"movslq	0x20(%%r15),%%r9 	\n\t	addq	%%r14,%%r9 	\n\t"\
	"movslq	0x1c(%%r15),%%r8 	\n\t	addq	%%r14,%%r8 	\n\t"\
	/* Want to stash O0-ptr in qword 1 of unused (vec_dbl*)__cc + 16 slot for later reload,
	but have no GPR to reload __cc into, so need some register-juggling to acomplish that: */\
	"movslq	(%%r15),%%r15 	\n\t"/* Can't use ADDQ (r15),r14 here since only want low 32 bits of (r15), so instead dereference r15 into-itself */\
	"addq	%%r15,%%r14		\n\t"/* r14 has O0, now load __cc into r15 */\
	"movq	%[__cc],%%r15	\n\t"/* cc1 */\
	"movq	%%r14,0x208(%%r15)	\n\t"/* ...and then stash O0 into qword 1 of that (vec_dbl)-width memloc. */\
	"vmovaps	%%ymm1,(%%r13)	\n\t	vmovaps	%%ymm7 ,0x20(%%r13)	\n\t"/* S1 */\
	"vmovaps	%%ymm2,(%%r12)	\n\t	vmovaps	%%ymm8 ,0x20(%%r12)	\n\t"/* S2 */\
	"vmovaps	%%ymm3,(%%r11)	\n\t	vmovaps	%%ymm9 ,0x20(%%r11)	\n\t"/* S3 */\
	"vmovaps	%%ymm4,(%%r10)	\n\t	vmovaps	%%ymm10,0x20(%%r10)	\n\t"/* S4 */\
	"vmovaps	%%ymm5,(%%r9 )	\n\t	vmovaps	%%ymm11,0x20(%%r9 )	\n\t"/* S5 */\
	"vmovaps	%%ymm6,(%%r8 )	\n\t	vmovaps	%%ymm12,0x20(%%r8 )	\n\t"/* S6 */\
/*
	C1 = t0 + cc1*t1 + cc2*t2 + cc3*t3 + cc4*t4 + cc5*t5 + cc6*t6
	C2 = t0 + cc2*t1 + cc4*t2 + cc6*t3 + cc5*t4 + cc3*t5 + cc1*t6
	C3 = t0 + cc3*t1 + cc6*t2 + cc4*t3 + cc1*t4 + cc2*t5 + cc5*t6
	C4 = t0 + cc4*t1 + cc5*t2 + cc1*t3 + cc3*t4 + cc6*t5 + cc2*t6
	C5 = t0 + cc5*t1 + cc3*t2 + cc2*t3 + cc6*t4 + cc1*t5 + cc4*t6
	C6 = t0 + cc6*t1 + cc1*t2 + cc5*t3 + cc2*t4 + cc4*t5 + cc3*t6
	B0 = t0 +     t1 +     t2 +     t3 +     t4 +     t5 +     t6	// X0

	In a 16-reg FMA model, makes sense to init the C-terms to the t0-summand,
	and keep the 12 C-terms and 2 of the 6 shared sincos data in registers. That uses 14 regs,
	leaving 2 for the 2 t*[r,i] data shared by each such block. Those need to be in-reg (at least
	in the cosine-mults section) because of the need to accumulate the DC terms: E.g. at end of the
	first block below (after doing the 12 C-term-update FMAs) we add the current values of B0r,i
	(which we init to t0r,i in the ASM version) to t1r,i, then write the result to memory and
	load t2r,i into the same 2 regs in preparation for the next such block.
*/\
/*** rax,rbx,rcx,rdx,rsi,rdi hold i1-6, mm0-2 hold I0,O0,2.0 ***/\
	"movq	0x200(%%r15),%%r8 		\n\t"/* reload I0 into r8 */\
		"vmovaps	(%%r8 ),%%ymm1		\n\t	vmovaps	0x20(%%r8 ),%%ymm7	\n\t"/* c1 = A0 */\
		"vmovaps	%%ymm1,%%ymm2		\n\t	vmovaps	%%ymm7 ,%%ymm8 		\n\t"/* c2 = A0 */\
		"vmovaps	%%ymm1,%%ymm3		\n\t	vmovaps	%%ymm7 ,%%ymm9 		\n\t"/* c3 = A0 */\
		"vmovaps	%%ymm1,%%ymm4		\n\t	vmovaps	%%ymm7 ,%%ymm10		\n\t"/* c4 = A0 */\
		"vmovaps	%%ymm1,%%ymm5		\n\t	vmovaps	%%ymm7 ,%%ymm11		\n\t"/* c5 = A0 */\
		"vmovaps	%%ymm1,%%ymm6		\n\t	vmovaps	%%ymm7 ,%%ymm12		\n\t"/* c6 = A0 */\
\
		"vmovaps	    (%%r15),%%ymm13	\n\t	vmovaps	0x20(%%r15),%%ymm14		\n\t"/* cc1,cc2 */\
\
		"vmovaps	    (%%rax),%%ymm0	\n\t	vmovaps	0x20(%%rax),%%ymm15	\n\t"/* t1 */\
	"vfmadd231pd     %%ymm13,%%ymm0,%%ymm1	\n\t vfmadd231pd     %%ymm13,%%ymm15,%%ymm7 	\n\t"/* C1 += cc1*t1 */\
	"vfmadd231pd     %%ymm14,%%ymm0,%%ymm2	\n\t vfmadd231pd     %%ymm14,%%ymm15,%%ymm8 	\n\t"/* C2 += cc2*t1 */\
	"vfmadd231pd 0x40(%%r15),%%ymm0,%%ymm3	\n\t vfmadd231pd 0x40(%%r15),%%ymm15,%%ymm9 	\n\t"/* C3 += cc3*t1 */\
	"vfmadd231pd 0x60(%%r15),%%ymm0,%%ymm4	\n\t vfmadd231pd 0x60(%%r15),%%ymm15,%%ymm10	\n\t"/* C4 += cc4*t1 */\
	"vfmadd231pd 0x80(%%r15),%%ymm0,%%ymm5	\n\t vfmadd231pd 0x80(%%r15),%%ymm15,%%ymm11	\n\t"/* C5 += cc5*t1 */\
	"vfmadd231pd 0xa0(%%r15),%%ymm0,%%ymm6	\n\t vfmadd231pd 0xa0(%%r15),%%ymm15,%%ymm12	\n\t"/* C6 += cc6*t1 */\
		"vaddpd (%%r8 ),%%ymm0,%%ymm0	\n\t	vaddpd 0x20(%%r8 ),%%ymm15,%%ymm15	\n\t"/* B0 += t1; */\
		"vmovaps	%%ymm0,(%%r8 )		\n\t	vmovaps	%%ymm15,0x20(%%r8 )			\n\t"/* Store B0 */\
\
		"vmovaps	    (%%rbx),%%ymm0	\n\t	vmovaps	0x20(%%rbx),%%ymm15	\n\t"/* t2 */\
	"vfmadd231pd     %%ymm14,%%ymm0,%%ymm1	\n\t vfmadd231pd     %%ymm14,%%ymm15,%%ymm7 	\n\t"/* C1 += cc2*t2 */\
	"vfmadd231pd 0x60(%%r15),%%ymm0,%%ymm2	\n\t vfmadd231pd 0x60(%%r15),%%ymm15,%%ymm8 	\n\t"/* C2 += cc4*t2 */\
	"vfmadd231pd 0xa0(%%r15),%%ymm0,%%ymm3	\n\t vfmadd231pd 0xa0(%%r15),%%ymm15,%%ymm9 	\n\t"/* C3 += cc6*t2 */\
	"vfmadd231pd 0x80(%%r15),%%ymm0,%%ymm4	\n\t vfmadd231pd 0x80(%%r15),%%ymm15,%%ymm10	\n\t"/* C4 += cc5*t2 */\
	"vfmadd231pd 0x40(%%r15),%%ymm0,%%ymm5	\n\t vfmadd231pd 0x40(%%r15),%%ymm15,%%ymm11	\n\t"/* C5 += cc3*t2 */\
	"vfmadd231pd     %%ymm13,%%ymm0,%%ymm6	\n\t vfmadd231pd     %%ymm13,%%ymm15,%%ymm12	\n\t"/* C6 += cc1*t2 */\
		"vaddpd (%%r8 ),%%ymm0,%%ymm0	\n\t	vaddpd 0x20(%%r8 ),%%ymm15,%%ymm15	\n\t"/* B0 += t2; */\
		"vmovaps	%%ymm0,(%%r8 )		\n\t	vmovaps	%%ymm15,0x20(%%r8 )			\n\t"/* Store B0 */\
\
		"vmovaps	    (%%rcx),%%ymm0	\n\t	vmovaps	0x20(%%rcx),%%ymm15	\n\t"/* t3 */\
	"vfmadd231pd 0x40(%%r15),%%ymm0,%%ymm1	\n\t vfmadd231pd 0x40(%%r15),%%ymm15,%%ymm7 	\n\t"/* C1 += cc3*t3 */\
	"vfmadd231pd 0xa0(%%r15),%%ymm0,%%ymm2	\n\t vfmadd231pd 0xa0(%%r15),%%ymm15,%%ymm8 	\n\t"/* C2 += cc6*t3 */\
	"vfmadd231pd 0x60(%%r15),%%ymm0,%%ymm3	\n\t vfmadd231pd 0x60(%%r15),%%ymm15,%%ymm9 	\n\t"/* C3 += cc4*t3 */\
	"vfmadd231pd     %%ymm13,%%ymm0,%%ymm4	\n\t vfmadd231pd     %%ymm13,%%ymm15,%%ymm10	\n\t"/* C4 += cc1*t3 */\
	"vfmadd231pd     %%ymm14,%%ymm0,%%ymm5	\n\t vfmadd231pd     %%ymm14,%%ymm15,%%ymm11	\n\t"/* C5 += cc2*t3 */\
	"vfmadd231pd 0x80(%%r15),%%ymm0,%%ymm6	\n\t vfmadd231pd 0x80(%%r15),%%ymm15,%%ymm12	\n\t"/* C6 += cc5*t3 */\
		"vaddpd (%%r8 ),%%ymm0,%%ymm0	\n\t	vaddpd 0x20(%%r8 ),%%ymm15,%%ymm15	\n\t"/* B0 += t3; */\
		"vmovaps	%%ymm0,(%%r8 )		\n\t	vmovaps	%%ymm15,0x20(%%r8 )			\n\t"/* Store B0 */\
\
		"vmovaps	    (%%rdx),%%ymm0	\n\t	vmovaps	0x20(%%rdx),%%ymm15	\n\t"/* t4 */\
	"vfmadd231pd 0x60(%%r15),%%ymm0,%%ymm1	\n\t vfmadd231pd 0x60(%%r15),%%ymm15,%%ymm7 	\n\t"/* C1 += cc4*t4 */\
	"vfmadd231pd 0x80(%%r15),%%ymm0,%%ymm2	\n\t vfmadd231pd 0x80(%%r15),%%ymm15,%%ymm8 	\n\t"/* C2 += cc5*t4 */\
	"vfmadd231pd     %%ymm13,%%ymm0,%%ymm3	\n\t vfmadd231pd     %%ymm13,%%ymm15,%%ymm9 	\n\t"/* C3 += cc1*t4 */\
	"vfmadd231pd 0x40(%%r15),%%ymm0,%%ymm4	\n\t vfmadd231pd 0x40(%%r15),%%ymm15,%%ymm10	\n\t"/* C4 += cc3*t4 */\
	"vfmadd231pd 0xa0(%%r15),%%ymm0,%%ymm5	\n\t vfmadd231pd 0xa0(%%r15),%%ymm15,%%ymm11	\n\t"/* C5 += cc6*t4 */\
	"vfmadd231pd     %%ymm14,%%ymm0,%%ymm6	\n\t vfmadd231pd     %%ymm14,%%ymm15,%%ymm12	\n\t"/* C6 += cc2*t4 */\
		"vaddpd (%%r8 ),%%ymm0,%%ymm0	\n\t	vaddpd 0x20(%%r8 ),%%ymm15,%%ymm15	\n\t"/* B0 += t4; */\
		"vmovaps	%%ymm0,(%%r8 )		\n\t	vmovaps	%%ymm15,0x20(%%r8 )			\n\t"/* Store B0 */\
\
		"vmovaps	    (%%rsi),%%ymm0	\n\t	vmovaps	0x20(%%rsi),%%ymm15	\n\t"/* t5 */\
	"vfmadd231pd 0x80(%%r15),%%ymm0,%%ymm1	\n\t vfmadd231pd 0x80(%%r15),%%ymm15,%%ymm7 	\n\t"/* C1 += cc5*t5 */\
	"vfmadd231pd 0x40(%%r15),%%ymm0,%%ymm2	\n\t vfmadd231pd 0x40(%%r15),%%ymm15,%%ymm8 	\n\t"/* C2 += cc3*t5 */\
	"vfmadd231pd     %%ymm14,%%ymm0,%%ymm3	\n\t vfmadd231pd     %%ymm14,%%ymm15,%%ymm9 	\n\t"/* C3 += cc2*t5 */\
	"vfmadd231pd 0xa0(%%r15),%%ymm0,%%ymm4	\n\t vfmadd231pd 0xa0(%%r15),%%ymm15,%%ymm10	\n\t"/* C4 += cc6*t5 */\
	"vfmadd231pd     %%ymm13,%%ymm0,%%ymm5	\n\t vfmadd231pd     %%ymm13,%%ymm15,%%ymm11	\n\t"/* C5 += cc1*t5 */\
	"vfmadd231pd 0x60(%%r15),%%ymm0,%%ymm6	\n\t vfmadd231pd 0x60(%%r15),%%ymm15,%%ymm12	\n\t"/* C6 += cc4*t5 */\
		"vaddpd (%%r8 ),%%ymm0,%%ymm0	\n\t	vaddpd 0x20(%%r8 ),%%ymm15,%%ymm15	\n\t"/* B0 += t5; */\
		"vmovaps	%%ymm0,(%%r8 )		\n\t	vmovaps	%%ymm15,0x20(%%r8 )			\n\t"/* Store B0 */\
\
		"vmovaps	    (%%rdi),%%ymm0	\n\t	vmovaps	0x20(%%rdi),%%ymm15	\n\t"/* t6 */\
	"vfmadd231pd 0xa0(%%r15),%%ymm0,%%ymm1	\n\t vfmadd231pd 0xa0(%%r15),%%ymm15,%%ymm7 	\n\t"/* C1 += cc6*t6 */\
	"vfmadd231pd     %%ymm13,%%ymm0,%%ymm2	\n\t vfmadd231pd     %%ymm13,%%ymm15,%%ymm8 	\n\t"/* C2 += cc1*t6 */\
	"vfmadd231pd 0x80(%%r15),%%ymm0,%%ymm3	\n\t vfmadd231pd 0x80(%%r15),%%ymm15,%%ymm9 	\n\t"/* C3 += cc5*t6 */\
	"vfmadd231pd     %%ymm14,%%ymm0,%%ymm4	\n\t vfmadd231pd     %%ymm14,%%ymm15,%%ymm10	\n\t"/* C4 += cc2*t6 */\
	"vfmadd231pd 0x60(%%r15),%%ymm0,%%ymm5	\n\t vfmadd231pd 0x60(%%r15),%%ymm15,%%ymm11	\n\t"/* C5 += cc4*t6 */\
	"vfmadd231pd 0x40(%%r15),%%ymm0,%%ymm6	\n\t vfmadd231pd 0x40(%%r15),%%ymm15,%%ymm12	\n\t"/* C6 += cc3*t6 */\
		"vaddpd (%%r8 ),%%ymm0,%%ymm0	\n\t	vaddpd 0x20(%%r8 ),%%ymm15,%%ymm15	\n\t"/* B0 += t6; */\
	"movq	0x208(%%r15),%%r8 		\n\t"/* reload O0 into r8 */\
		"vmovaps	%%ymm0,(%%r8 )		\n\t	vmovaps	%%ymm15,0x20(%%r8 )			\n\t"/* Store B0, now into output slot */\
\
	"movq	%[out0],%%r14		\n\t"/* r14 holds Output base-pointer */\
	"movq	%[o_off],%%r15		\n\t"/* r15 holds pointer to uint32[] Outnput pointer-offset-stride array */\
	"movslq	    (%%r15),%%r8 	\n\t"\
	/* Compute o1-6 from offsets and recompute o7; o8-c still in r9-13): */\
	"movslq	0x04(%%r15),%%rax	\n\t"\
	"movslq	0x08(%%r15),%%rbx	\n\t"\
	"movslq	0x0c(%%r15),%%rcx	\n\t"\
	"movslq	0x10(%%r15),%%rdx	\n\t"\
	"movslq	0x14(%%r15),%%rsi	\n\t"\
	"movslq	0x18(%%r15),%%rdi	\n\t"\
	"movslq	0x1c(%%r15),%%r8 	\n\t"\
	"addq	%%r14,%%r8 	\n\t"\
	"addq	%%r14,%%rax	\n\t"\
	"addq	%%r14,%%rbx	\n\t"\
	"addq	%%r14,%%rcx	\n\t"	/* Have enough registers to also hold 3 of the 6 real parts of S-terms: */\
	"addq	%%r14,%%rdx	\n\t	vmovaps	(%%r10),%%ymm13		\n\t"/* s4r */\
	"addq	%%r14,%%rsi	\n\t	vmovaps	(%%r9 ),%%ymm14		\n\t"/* s5r */\
	"addq	%%r14,%%rdi	\n\t	vmovaps	(%%r8 ),%%ymm15		\n\t"/* s5r */\
	/* ...Do the lower-half-index portion of the closing wave of radix-2 butterflies: */\
		"vsubpd 0x20(%%r13),%%ymm1,%%ymm1		\n\t	 vaddpd (%%r13),%%ymm7 ,%%ymm7 	\n\t"/* B1 = C1 + I*S1 */\
		"vsubpd 0x20(%%r12),%%ymm2,%%ymm2		\n\t	 vaddpd (%%r12),%%ymm8 ,%%ymm8 	\n\t"/* B2 = C2 + I*S2 */\
		"vsubpd 0x20(%%r11),%%ymm3,%%ymm3		\n\t	 vaddpd (%%r11),%%ymm9 ,%%ymm9 	\n\t"/* B3 = C3 + I*S3 */\
		"vsubpd 0x20(%%r10),%%ymm4,%%ymm4		\n\t	 vaddpd %%ymm13,%%ymm10,%%ymm10	\n\t"/* B4 = C4 + I*S4 */\
		"vsubpd 0x20(%%r9 ),%%ymm5,%%ymm5		\n\t	 vaddpd %%ymm14,%%ymm11,%%ymm11	\n\t"/* B5 = C5 + I*S5 */\
		"vsubpd 0x20(%%r8 ),%%ymm6,%%ymm6		\n\t	 vaddpd %%ymm15,%%ymm12,%%ymm12	\n\t"/* B6 = C6 + I*S6 */\
	/* Write lower-half outputs back to memory... */\
	"movq	%[__cc],%%r15				\n\t	vmovaps	-0x20(%%r15),%%ymm0		\n\t"/* two */\
		"vmovaps	%%ymm1,(%%rax)		\n\t"/* B1r */"		vmovaps	%%ymm7 ,0x20(%%rax)	\n\t"/* B1i */\
		"vmovaps	%%ymm2,(%%rbx)		\n\t"/* B2r */"		vmovaps	%%ymm8 ,0x20(%%rbx)	\n\t"/* B2i */\
		"vmovaps	%%ymm3,(%%rcx)		\n\t"/* B3r */"		vmovaps	%%ymm9 ,0x20(%%rcx)	\n\t"/* B3i */\
		"vmovaps	%%ymm4,(%%rdx)		\n\t"/* B4r */"		vmovaps	%%ymm10,0x20(%%rdx)	\n\t"/* B4i */\
		"vmovaps	%%ymm5,(%%rsi)		\n\t"/* B5r */"		vmovaps	%%ymm11,0x20(%%rsi)	\n\t"/* B5i */\
		"vmovaps	%%ymm6,(%%rdi)		\n\t"/* B6r */"		vmovaps	%%ymm12,0x20(%%rdi)	\n\t"/* B6i */\
	/* ...Do the upper-half-index portion of the radix-2 butterflies: */\
		" vfmadd231pd 0x20(%%r13),%%ymm0,%%ymm1		\n\t	vfnmadd231pd (%%r13),%%ymm0,%%ymm7 	\n\t"/* Bc = C1 - I*S1 */\
		" vfmadd231pd 0x20(%%r12),%%ymm0,%%ymm2		\n\t	vfnmadd231pd (%%r12),%%ymm0,%%ymm8 	\n\t"/* Bb = C2 - I*S2 */\
		" vfmadd231pd 0x20(%%r11),%%ymm0,%%ymm3		\n\t	vfnmadd231pd (%%r11),%%ymm0,%%ymm9 	\n\t"/* Ba = C3 - I*S3 */\
		" vfmadd231pd 0x20(%%r10),%%ymm0,%%ymm4		\n\t	vfnmadd231pd %%ymm13,%%ymm0,%%ymm10	\n\t"/* B9 = C4 - I*S4 */\
		" vfmadd231pd 0x20(%%r9 ),%%ymm0,%%ymm5		\n\t	vfnmadd231pd %%ymm14,%%ymm0,%%ymm11	\n\t"/* B8 = C5 - I*S5 */\
		" vfmadd231pd 0x20(%%r8 ),%%ymm0,%%ymm6		\n\t	vfnmadd231pd %%ymm15,%%ymm0,%%ymm12	\n\t"/* B7 = C6 - I*S6 */\
	/* ...And write the upper-half outputs back to memory. */\
		"vmovaps	%%ymm1,(%%r13)		\n\t"/* Bcr */"		vmovaps	%%ymm7 ,0x20(%%r13)	\n\t"/* Bci */\
		"vmovaps	%%ymm2,(%%r12)		\n\t"/* Bbr */"		vmovaps	%%ymm8 ,0x20(%%r12)	\n\t"/* Bbi */\
		"vmovaps	%%ymm3,(%%r11)		\n\t"/* Bar */"		vmovaps	%%ymm9 ,0x20(%%r11)	\n\t"/* Bai */\
		"vmovaps	%%ymm4,(%%r10)		\n\t"/* B9r */"		vmovaps	%%ymm10,0x20(%%r10)	\n\t"/* B9i */\
		"vmovaps	%%ymm5,(%%r9 )		\n\t"/* B8r */"		vmovaps	%%ymm11,0x20(%%r9 )	\n\t"/* B8i */\
		"vmovaps	%%ymm6,(%%r8 )		\n\t"/* B7r */"		vmovaps	%%ymm12,0x20(%%r8 )	\n\t"/* B7i */\
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

   #else		// [2] [LOACC, no-flags default] based on FMAized version of same van Buskirk tan-DFT used in no-FMA mode.

	// 4 Dec 2014: First-look FMAization changes opcount from [164 ADD, 74 MUL, 170 memref] ==> [120 ADD, 78 FMA (incl 34 MUL), 168 memref]
	#define SSE2_RADIX_13_DFT(Xin0,Xi_off, Xcc, Xout0,Xo_off)\
	{\
	__asm__ volatile (\
		"movq	%[__cc],%%rsi	\n\t"\
	/* rax,rcx hold [In|Out]put base-pointers: */\
	"movq	%[in0]  ,%%rax		\n\t	movq	%[out0]  ,%%rcx		\n\t"\
	/* rbx,rdx hold pointers to uint32[] [In|Out]put pointer-offset-stride arrays: */\
	"movq	%[i_off],%%rbx		\n\t	movq	%[o_off],%%rdx		\n\t"\
	/* Compute i1-6 from offsets:		I0,O0 get stashed in r8,r9: */\
	"movslq	0x04(%%rbx),%%r15	\n\t	movslq	(%%rbx),%%r8	\n\t"\
	"movslq	0x08(%%rbx),%%r14	\n\t	movslq	(%%rdx),%%r9	\n\t"\
	"movslq	0x0c(%%rbx),%%r13	\n\t	addq	%%rax,%%r8		\n\t"/* I0 */\
	"movslq	0x10(%%rbx),%%r12	\n\t	addq	%%rcx,%%r9		\n\t"/* O0 */\
	"movslq	0x14(%%rbx),%%r11	\n\t"\
	"movslq	0x18(%%rbx),%%r10	\n\t"\
	/* xr-terms:                                                      yi-terms: */\
		"addq	%%rax,%%r15	\n\t	vmovaps	(%%r15),%%ymm1			\n\t	vmovaps	0x20(%%r15),%%ymm9 	\n\t"\
		"addq	%%rax,%%r14	\n\t	vmovaps	(%%r14),%%ymm3			\n\t	vmovaps	0x20(%%r14),%%ymm10	\n\t"\
		"addq	%%rax,%%r13	\n\t	vmovaps	(%%r13),%%ymm6			\n\t	vmovaps	0x20(%%r13),%%ymm13	\n\t"\
		"addq	%%rax,%%r12	\n\t	vmovaps	(%%r12),%%ymm4			\n\t	vmovaps	0x20(%%r12),%%ymm11	\n\t"\
		"addq	%%rax,%%r11	\n\t	vmovaps	(%%r11),%%ymm5			\n\t	vmovaps	0x20(%%r11),%%ymm12	\n\t"\
		"addq	%%rax,%%r10	\n\t	vmovaps	(%%r10),%%ymm7			\n\t	vmovaps	0x20(%%r10),%%ymm14	\n\t"\
	/* Compute iC-7 from offsets: */\
	"movslq	0x30(%%rbx),%%r15	\n\t"\
	"movslq	0x2c(%%rbx),%%r14	\n\t"\
	"movslq	0x28(%%rbx),%%r13	\n\t"\
	"movslq	0x24(%%rbx),%%r12	\n\t"\
	"movslq	0x20(%%rbx),%%r11	\n\t"\
	"movslq	0x1c(%%rbx),%%r10	\n\t"\
		"addq	%%rax,%%r15	\n\t	vaddpd	(%%r15),%%ymm1,%%ymm1	\n\t	vsubpd	0x20(%%r15),%%ymm9 ,%%ymm9 	\n\t"\
		"addq	%%rax,%%r14	\n\t	vaddpd	(%%r14),%%ymm3,%%ymm3	\n\t	vsubpd	0x20(%%r14),%%ymm10,%%ymm10	\n\t"\
		"addq	%%rax,%%r13	\n\t	vaddpd	(%%r13),%%ymm6,%%ymm6	\n\t	vsubpd	0x20(%%r13),%%ymm13,%%ymm13	\n\t"\
		"addq	%%rax,%%r12	\n\t	vaddpd	(%%r12),%%ymm4,%%ymm4	\n\t	vsubpd	0x20(%%r12),%%ymm11,%%ymm11	\n\t"\
		"addq	%%rax,%%r11	\n\t	vaddpd	(%%r11),%%ymm5,%%ymm5	\n\t	vsubpd	0x20(%%r11),%%ymm12,%%ymm12	\n\t"\
		"addq	%%rax,%%r10	\n\t	vaddpd	(%%r10),%%ymm7,%%ymm7	\n\t	vsubpd	0x20(%%r10),%%ymm14,%%ymm14	\n\t"\
	/*** Done with i1-C for now ***/\
		/* Next part (up to 'lcol section' comment) identical for Re,Im-part sections: */\
	"vmovaps	-0x20(%%rsi),%%ymm0		\n\t"/* two */\
		"														vmovaps	%%ymm9 ,%%ymm15	\n\t	vmovaps	%%ymm10,%%ymm8	\n\t"\
		"vsubpd	%%ymm5,%%ymm1,%%ymm1					\n\t			vsubpd	%%ymm11,%%ymm15,%%ymm15			\n\t"\
		"vsubpd	%%ymm6,%%ymm3,%%ymm3					\n\t			vaddpd	%%ymm12,%%ymm8 ,%%ymm8 			\n\t"\
		"vsubpd	%%ymm7,%%ymm4,%%ymm4					\n\t			vaddpd	%%ymm11,%%ymm9 ,%%ymm9 			\n\t"\
	"vfmadd132pd	%%ymm0,%%ymm1,%%ymm5				\n\t			vaddpd	%%ymm13,%%ymm15,%%ymm15			\n\t"\
	"vfmadd132pd	%%ymm0,%%ymm3,%%ymm6				\n\t			vaddpd	%%ymm11,%%ymm13,%%ymm13			\n\t"\
	"vfmadd132pd	%%ymm0,%%ymm4,%%ymm7				\n\t			vaddpd	%%ymm14,%%ymm8 ,%%ymm8 			\n\t"\
		"vmovaps	%%ymm5,%%ymm2						\n\t			vsubpd	%%ymm12,%%ymm10,%%ymm10			\n\t"\
		"vaddpd	%%ymm6,%%ymm2,%%ymm2					\n\t			vsubpd	%%ymm12,%%ymm14,%%ymm14			\n\t"\
		"vaddpd	%%ymm7,%%ymm2,%%ymm2					\n\t			vmovaps	%%ymm8 ,%%ymm12					\n\t"\
/*I0:*/	"vaddpd	(%%r8 ),%%ymm2,%%ymm0					\n\t			vmovaps	%%ymm15,%%ymm11					\n\t"\
/*O0:*/	"vmovaps	%%ymm0,(%%r9 )						\n\t			vmulpd	0x80(%%rsi),%%ymm15,%%ymm15	\n\t"\
	"vfmadd231pd	(%%rsi),%%ymm2,%%ymm0				\n\t			vmulpd	0x80(%%rsi),%%ymm8 ,%%ymm8 	\n\t"\
		"												\n\t		vfmsub132pd	0x60(%%rsi),%%ymm15,%%ymm12	\n\t"\
		"vmovaps	0x20(%%rsi),%%ymm2					\n\t		vfmadd132pd	0x60(%%rsi),%%ymm8 ,%%ymm11	\n\t"\
	/* lcol|rcol o-offsets 1-3|B-C - compute in order used below, store into r10-12|13-14: */\
	"movslq	0x30(%%rdx),%%r15	\n\t	addq	%%rcx,%%r15	\n\t"/* oC */\
	"movslq	0x04(%%rdx),%%r10	\n\t	addq	%%rcx,%%r10	\n\t"/* o1 */\
	"movslq	0x08(%%rdx),%%r11	\n\t	addq	%%rcx,%%r11	\n\t"/* o2 */\
	"movslq	0x0c(%%rdx),%%r12	\n\t	addq	%%rcx,%%r12	\n\t"/* o3 */\
	"movslq	0x2c(%%rdx),%%r14	\n\t	addq	%%rcx,%%r14	\n\t"/* oB */\
		"vmulpd	%%ymm2,%%ymm6,%%ymm6					\n\t			vmovaps	%%ymm12,(%%r15)			\n\t"\
		"vmulpd	%%ymm2,%%ymm5,%%ymm5					\n\t			vmovaps	%%ymm9 ,%%ymm8 			\n\t"\
		"vmulpd	%%ymm2,%%ymm7,%%ymm7					\n\t			vmovaps	%%ymm13,%%ymm12			\n\t"\
		"vmovaps %%ymm6,(%%r10)							\n\t			vaddpd	%%ymm10,%%ymm8 ,%%ymm8 			\n\t"\
		"vmovaps %%ymm5,(%%r11)							\n\t			vaddpd	%%ymm14,%%ymm12,%%ymm12			\n\t"\
		"vmovaps %%ymm7,(%%r12)							\n\t"\
		"vmovaps 0x100(%%rsi),%%ymm2					\n\t			vmovaps	%%ymm12,%%ymm15			\n\t"\
	"vfmadd132pd	%%ymm2,%%ymm0,%%ymm5				\n\t		vfmadd132pd 0x160(%%rsi),%%ymm8 ,%%ymm12	\n\t"\
	"vfmadd132pd	%%ymm2,%%ymm0,%%ymm7				\n\t			vmulpd 0x1c0(%%rsi),%%ymm8 ,%%ymm8 		\n\t"\
	"vfmadd132pd	%%ymm2,%%ymm0,%%ymm6				\n\t"\
		"												\n\t			vaddpd	%%ymm15,%%ymm8 ,%%ymm8 			\n\t"\
		"												\n\t			vmulpd	0xa0(%%rsi),%%ymm12,%%ymm12	\n\t"\
		"												\n\t			vmulpd	0xa0(%%rsi),%%ymm8 ,%%ymm8 		\n\t"\
		"vaddpd	(%%r10),%%ymm5,%%ymm5					\n\t			vmovaps	%%ymm10,(%%r14)			\n\t"\
		"vaddpd	(%%r11),%%ymm7,%%ymm7					\n\t			vmovaps	%%ymm14,%%ymm15			\n\t"\
		"vaddpd	(%%r12),%%ymm6,%%ymm6					\n\t			vmulpd 0x180(%%rsi),%%ymm14,%%ymm14	\n\t"\
		"vmovaps 0x140(%%rsi),%%ymm2					\n\t			vmulpd 0x1e0(%%rsi),%%ymm10,%%ymm10	\n\t"\
		"vmovaps 0x120(%%rsi),%%ymm0					\n\t			vaddpd	(%%r14),%%ymm14,%%ymm14			\n\t"\
		"vmovaps %%ymm4,(%%r10)							\n\t			vaddpd	%%ymm15,%%ymm10,%%ymm10			\n\t"\
		"vmovaps %%ymm1,(%%r11)							\n\t			vmulpd	0xc0(%%rsi),%%ymm14,%%ymm14	\n\t"\
		"vmovaps %%ymm3,(%%r12)							\n\t			vmulpd	0xc0(%%rsi),%%ymm10,%%ymm10	\n\t"\
		"												\n\t			vaddpd	%%ymm12,%%ymm14,%%ymm14			\n\t"\
		"												\n\t			vaddpd	%%ymm8 ,%%ymm10,%%ymm10			\n\t"\
		"												\n\t			vmovaps	%%ymm9 ,(%%r14)			\n\t"\
	"vfmadd213pd	(%%r12),%%ymm2,%%ymm4				\n\t			vmovaps	%%ymm13,%%ymm15			\n\t"\
	"vfmsub213pd	(%%r10),%%ymm2,%%ymm1				\n\t			vmulpd 0x1a0(%%rsi),%%ymm13,%%ymm13	\n\t"\
	"vfmadd213pd	(%%r11),%%ymm2,%%ymm3				\n\t			vmulpd 0x200(%%rsi),%%ymm9 ,%%ymm9 		\n\t"\
		"												\n\t			vaddpd	(%%r14),%%ymm13,%%ymm13			\n\t"\
		"												\n\t			vaddpd	%%ymm15,%%ymm9 ,%%ymm9 			\n\t"\
		"												\n\t			vmulpd	0xe0(%%rsi),%%ymm13,%%ymm13	\n\t"\
		"vmovaps 0x40(%%rsi),%%ymm2						\n\t			vmulpd	0xe0(%%rsi),%%ymm9 ,%%ymm9 		\n\t"\
	"vfmadd213pd	(%%r11),%%ymm0,%%ymm4				\n\t			vaddpd	%%ymm12,%%ymm13,%%ymm13			\n\t"\
	"vfmsub213pd	(%%r12),%%ymm0,%%ymm1				\n\t			vaddpd	%%ymm8 ,%%ymm9 ,%%ymm9 			\n\t"\
	"vfmsub213pd	(%%r10),%%ymm0,%%ymm3				\n\t			vmovaps	(%%r15),%%ymm12			\n\t"\
	"vmulpd	-0x20(%%rsi),%%ymm2,%%ymm0					\n\t			vmovaps	%%ymm11,%%ymm15			\n\t"\
	"vfnmadd231pd	%%ymm4,%%ymm2,%%ymm6				\n\t			vmovaps	%%ymm12,%%ymm8 			\n\t"\
	"vfnmadd231pd	%%ymm1,%%ymm2,%%ymm7				\n\t			vaddpd	%%ymm14,%%ymm15,%%ymm15			\n\t"\
	"vfnmadd231pd	%%ymm3,%%ymm2,%%ymm5				\n\t			vaddpd	%%ymm13,%%ymm8 ,%%ymm8 			\n\t"\
		"												\n\t			vsubpd	%%ymm11,%%ymm14,%%ymm14			\n\t"\
		"												\n\t			vsubpd	%%ymm12,%%ymm13,%%ymm13			\n\t"\
		"												\n\t			vaddpd	%%ymm10,%%ymm14,%%ymm14			\n\t"\
	"vfmadd213pd	%%ymm6,%%ymm0,%%ymm4				\n\t			vaddpd	%%ymm9 ,%%ymm13,%%ymm13			\n\t"\
	"vfmadd213pd	%%ymm7,%%ymm0,%%ymm1				\n\t			vaddpd	%%ymm11,%%ymm10,%%ymm10			\n\t"\
	"vfmadd213pd	%%ymm5,%%ymm0,%%ymm3				\n\t			vaddpd	%%ymm9 ,%%ymm12,%%ymm12			\n\t"\
	/* yi-data in ymm8,10,12,13,14,15; ymm0,2,9,11 free: */\
/*** o1-3|B,C already in r10-12|14,15: ***/\
	/*** o3 already in r12; clobber rsi for now, leaving rsi,rdi,r13 for o6,7,A: ***/\
	"movslq	0x18(%%rdx),%%rsi	\n\t	addq	%%rcx,%%rsi	\n\t"/* o6 */\
	"movslq	0x1c(%%rdx),%%rdi	\n\t	addq	%%rcx,%%rdi	\n\t"/* o7 */\
	"movslq	0x28(%%rdx),%%r13	\n\t	addq	%%rcx,%%r13	\n\t"/* oA */\
		"vmovaps	%%ymm4 ,%%ymm0						\n\t			vmovaps	%%ymm7 ,%%ymm9 			\n\t"\
		"vsubpd	%%ymm15,%%ymm4,%%ymm4					\n\t			vsubpd	%%ymm8 ,%%ymm7,%%ymm7	\n\t"\
		"vaddpd	%%ymm15,%%ymm0,%%ymm0					\n\t			vaddpd	%%ymm8 ,%%ymm9,%%ymm9	\n\t"\
		"vmovaps	%%ymm4 ,(%%rsi)						\n\t			vmovaps	%%ymm7 ,(%%r12)	\n\t"\
		"vmovaps	%%ymm0 ,(%%rdi)						\n\t			vmovaps	%%ymm9 ,(%%r13)	\n\t"\
	/* Compute o8,5,4,9, store into rsi,rdi,r12-13: */\
	"movslq	0x20(%%rdx),%%rsi	\n\t	addq	%%rcx,%%rsi	\n\t"/* o8 */\
	"movslq	0x14(%%rdx),%%rdi	\n\t	addq	%%rcx,%%rdi	\n\t"/* o5 */\
	"movslq	0x10(%%rdx),%%r12	\n\t	addq	%%rcx,%%r12	\n\t"/* o4 */\
	"movslq	0x24(%%rdx),%%r13	\n\t	addq	%%rcx,%%r13	\n\t"/* o9 */\
		"vmovaps	%%ymm5 ,%%ymm0						\n\t			vmovaps	%%ymm6 ,%%ymm9 			\n\t"\
		"vsubpd	%%ymm14,%%ymm5,%%ymm5					\n\t			vsubpd	%%ymm13,%%ymm6,%%ymm6	\n\t"\
		"vaddpd	%%ymm14,%%ymm0,%%ymm0					\n\t			vaddpd	%%ymm13,%%ymm9,%%ymm9	\n\t"\
		"vmovaps	%%ymm5 ,(%%rsi)						\n\t			vmovaps	%%ymm6 ,(%%r12)	\n\t"\
		"vmovaps	%%ymm0 ,(%%rdi)						\n\t			vmovaps	%%ymm9 ,(%%r13)	\n\t"\
	/*** o2,B,1,C already in r11,r14,r10,r15: ***/\
		"vmovaps	%%ymm1 ,%%ymm0						\n\t			vmovaps	%%ymm3 ,%%ymm9 			\n\t"\
		"vsubpd	%%ymm10,%%ymm1,%%ymm1					\n\t			vsubpd	%%ymm12,%%ymm3,%%ymm3	\n\t"\
		"vaddpd	%%ymm10,%%ymm0,%%ymm0					\n\t			vaddpd	%%ymm12,%%ymm9,%%ymm9	\n\t"\
		"vmovaps	%%ymm1 ,(%%r11)						\n\t			vmovaps	%%ymm3 ,(%%r10)	\n\t"\
		"vmovaps	%%ymm0 ,(%%r14)						\n\t			vmovaps	%%ymm9 ,(%%r15)	\n\t"\
		/****************************************************************/\
		/*                      IMAG PARTS:                             */\
		/****************************************************************/\
		"movq	%[__cc],%%rsi	\n\t"/* Reload cc-ptr, which we clobbered in above closing-butterfly section */\
		"addq	$0x20,%%r8		\n\t"/* I0.imag */\
		"addq	$0x20,%%r9		\n\t"/* O0.imag */\
	/* Compute i1-6 from offsets: */\
	"movslq	0x04(%%rbx),%%r15	\n\t"\
	"movslq	0x08(%%rbx),%%r14	\n\t"\
	"movslq	0x0c(%%rbx),%%r13	\n\t"\
	"movslq	0x10(%%rbx),%%r12	\n\t"\
	"movslq	0x14(%%rbx),%%r11	\n\t"\
	"movslq	0x18(%%rbx),%%r10	\n\t"\
	/* xi-terms:                                                      yr-terms: */\
		"addq	%%rax,%%r15	\n\t	vmovaps	0x20(%%r15),%%ymm1			\n\t	vmovaps	(%%r15),%%ymm9 	\n\t"\
		"addq	%%rax,%%r14	\n\t	vmovaps	0x20(%%r14),%%ymm3			\n\t	vmovaps	(%%r14),%%ymm10	\n\t"\
		"addq	%%rax,%%r13	\n\t	vmovaps	0x20(%%r13),%%ymm6			\n\t	vmovaps	(%%r13),%%ymm13	\n\t"\
		"addq	%%rax,%%r12	\n\t	vmovaps	0x20(%%r12),%%ymm4			\n\t	vmovaps	(%%r12),%%ymm11	\n\t"\
		"addq	%%rax,%%r11	\n\t	vmovaps	0x20(%%r11),%%ymm5			\n\t	vmovaps	(%%r11),%%ymm12	\n\t"\
		"addq	%%rax,%%r10	\n\t	vmovaps	0x20(%%r10),%%ymm7			\n\t	vmovaps	(%%r10),%%ymm14	\n\t"\
	/* Compute iC-7 from offsets: */\
	"movslq	0x30(%%rbx),%%r15	\n\t"\
	"movslq	0x2c(%%rbx),%%r14	\n\t"\
	"movslq	0x28(%%rbx),%%r13	\n\t"\
	"movslq	0x24(%%rbx),%%r12	\n\t"\
	"movslq	0x20(%%rbx),%%r11	\n\t"\
	"movslq	0x1c(%%rbx),%%r10	\n\t"\
		"addq	%%rax,%%r15	\n\t	vaddpd	0x20(%%r15),%%ymm1,%%ymm1	\n\t	vsubpd	(%%r15),%%ymm9 ,%%ymm9 	\n\t"\
		"addq	%%rax,%%r14	\n\t	vaddpd	0x20(%%r14),%%ymm3,%%ymm3	\n\t	vsubpd	(%%r14),%%ymm10,%%ymm10	\n\t"\
		"addq	%%rax,%%r13	\n\t	vaddpd	0x20(%%r13),%%ymm6,%%ymm6	\n\t	vsubpd	(%%r13),%%ymm13,%%ymm13	\n\t"\
		"addq	%%rax,%%r12	\n\t	vaddpd	0x20(%%r12),%%ymm4,%%ymm4	\n\t	vsubpd	(%%r12),%%ymm11,%%ymm11	\n\t"\
		"addq	%%rax,%%r11	\n\t	vaddpd	0x20(%%r11),%%ymm5,%%ymm5	\n\t	vsubpd	(%%r11),%%ymm12,%%ymm12	\n\t"\
		"addq	%%rax,%%r10	\n\t	vaddpd	0x20(%%r10),%%ymm7,%%ymm7	\n\t	vsubpd	(%%r10),%%ymm14,%%ymm14	\n\t"\
/*** Done w/i1-C - Still need 5 GPRs [rcx,rdx,rsi,r8,r9], to hold out0,o_off,__c,I0,O0 ***/\
		/* Next part (up to 'lcol section' comment) identical for Re,Im-part sections: */\
	"vmovaps	-0x20(%%rsi),%%ymm0		\n\t"/* two */\
		"														vmovaps	%%ymm9 ,%%ymm15	\n\t	vmovaps	%%ymm10,%%ymm8	\n\t"\
		"vsubpd	%%ymm5,%%ymm1,%%ymm1					\n\t			vsubpd	%%ymm11,%%ymm15,%%ymm15			\n\t"\
		"vsubpd	%%ymm6,%%ymm3,%%ymm3					\n\t			vaddpd	%%ymm12,%%ymm8 ,%%ymm8 			\n\t"\
		"vsubpd	%%ymm7,%%ymm4,%%ymm4					\n\t			vaddpd	%%ymm11,%%ymm9 ,%%ymm9 			\n\t"\
	"vfmadd132pd	%%ymm0,%%ymm1,%%ymm5				\n\t			vaddpd	%%ymm13,%%ymm15,%%ymm15			\n\t"\
	"vfmadd132pd	%%ymm0,%%ymm3,%%ymm6				\n\t			vaddpd	%%ymm11,%%ymm13,%%ymm13			\n\t"\
	"vfmadd132pd	%%ymm0,%%ymm4,%%ymm7				\n\t			vaddpd	%%ymm14,%%ymm8 ,%%ymm8 			\n\t"\
		"vmovaps	%%ymm5,%%ymm2						\n\t			vsubpd	%%ymm12,%%ymm10,%%ymm10			\n\t"\
		"vaddpd	%%ymm6,%%ymm2,%%ymm2					\n\t			vsubpd	%%ymm12,%%ymm14,%%ymm14			\n\t"\
		"vaddpd	%%ymm7,%%ymm2,%%ymm2					\n\t			vmovaps	%%ymm8 ,%%ymm12					\n\t"\
/*I0:*/	"vaddpd	(%%r8 ),%%ymm2,%%ymm0					\n\t			vmovaps	%%ymm15,%%ymm11					\n\t"\
/*O0:*/	"vmovaps	%%ymm0,(%%r9 )						\n\t			vmulpd	0x80(%%rsi),%%ymm15,%%ymm15	\n\t"\
/*** Done w/I0,O0 - At this point could use any 12 GPRs of 13 available [14 minus rsi] to hold o1-C, if desired ***/\
	"vfmadd231pd	(%%rsi),%%ymm2,%%ymm0				\n\t			vmulpd	0x80(%%rsi),%%ymm8 ,%%ymm8 	\n\t"\
		"												\n\t		vfmsub132pd	0x60(%%rsi),%%ymm15,%%ymm12	\n\t"\
		"vmovaps	0x20(%%rsi),%%ymm2					\n\t		vfmadd132pd	0x60(%%rsi),%%ymm8 ,%%ymm11	\n\t"\
	/* lcol|rcol o-offsets A-C|1-2 - compute in order used below, store into r13-15|10-11: */\
	"movslq	0x04(%%rdx),%%r10	\n\t	addq	%%rcx,%%r10	\n\t"/* o1 */\
	"movslq	0x30(%%rdx),%%r15	\n\t	addq	%%rcx,%%r15	\n\t"/* oC */\
	"movslq	0x2c(%%rdx),%%r14	\n\t	addq	%%rcx,%%r14	\n\t"/* oB */\
	"movslq	0x28(%%rdx),%%r13	\n\t	addq	%%rcx,%%r13	\n\t"/* oA */\
	"movslq	0x08(%%rdx),%%r11	\n\t	addq	%%rcx,%%r11	\n\t"/* o2 */\
		/* Next part (up to 'yi-data' comment) identical for Re,Im, except for replacements (%%r10-15) ==> 0x20(%%r15-10): */\
		"vmulpd	%%ymm2,%%ymm6,%%ymm6					\n\t			vmovaps	%%ymm12,0x20(%%r10)			\n\t"\
		"vmulpd	%%ymm2,%%ymm5,%%ymm5					\n\t			vmovaps	%%ymm9 ,%%ymm8 			\n\t"\
		"vmulpd	%%ymm2,%%ymm7,%%ymm7					\n\t			vmovaps	%%ymm13,%%ymm12			\n\t"\
		"vmovaps %%ymm6,0x20(%%r15)						\n\t			vaddpd	%%ymm10,%%ymm8 ,%%ymm8 			\n\t"\
		"vmovaps %%ymm5,0x20(%%r14)						\n\t			vaddpd	%%ymm14,%%ymm12,%%ymm12			\n\t"\
		"vmovaps %%ymm7,0x20(%%r13)						\n\t"\
		"vmovaps 0x100(%%rsi),%%ymm2					\n\t			vmovaps	%%ymm12,%%ymm15			\n\t"\
	"vfmadd132pd	%%ymm2,%%ymm0,%%ymm5				\n\t		vfmadd132pd 0x160(%%rsi),%%ymm8 ,%%ymm12	\n\t"\
	"vfmadd132pd	%%ymm2,%%ymm0,%%ymm7				\n\t			vmulpd 0x1c0(%%rsi),%%ymm8 ,%%ymm8 		\n\t"\
	"vfmadd132pd	%%ymm2,%%ymm0,%%ymm6				\n\t"\
		"												\n\t			vaddpd	%%ymm15,%%ymm8 ,%%ymm8 			\n\t"\
		"												\n\t			vmulpd	0xa0(%%rsi),%%ymm12,%%ymm12	\n\t"\
		"												\n\t			vmulpd	0xa0(%%rsi),%%ymm8 ,%%ymm8 		\n\t"\
		"vaddpd	0x20(%%r15),%%ymm5,%%ymm5				\n\t			vmovaps	%%ymm10,0x20(%%r11)			\n\t"\
		"vaddpd	0x20(%%r14),%%ymm7,%%ymm7				\n\t			vmovaps	%%ymm14,%%ymm15			\n\t"\
		"vaddpd	0x20(%%r13),%%ymm6,%%ymm6				\n\t			vmulpd 0x180(%%rsi),%%ymm14,%%ymm14	\n\t"\
		"vmovaps 0x140(%%rsi),%%ymm2					\n\t			vmulpd 0x1e0(%%rsi),%%ymm10,%%ymm10	\n\t"\
		"vmovaps 0x120(%%rsi),%%ymm0					\n\t			vaddpd	0x20(%%r11),%%ymm14,%%ymm14			\n\t"\
		"vmovaps %%ymm4,0x20(%%r15)						\n\t			vaddpd	%%ymm15,%%ymm10,%%ymm10			\n\t"\
		"vmovaps %%ymm1,0x20(%%r14)						\n\t			vmulpd	0xc0(%%rsi),%%ymm14,%%ymm14	\n\t"\
		"vmovaps %%ymm3,0x20(%%r13)						\n\t			vmulpd	0xc0(%%rsi),%%ymm10,%%ymm10	\n\t"\
		"												\n\t			vaddpd	%%ymm12,%%ymm14,%%ymm14			\n\t"\
		"												\n\t			vaddpd	%%ymm8 ,%%ymm10,%%ymm10			\n\t"\
		"												\n\t			vmovaps	%%ymm9 ,0x20(%%r11)			\n\t"\
	"vfmadd213pd	0x20(%%r13),%%ymm2,%%ymm4			\n\t			vmovaps	%%ymm13,%%ymm15			\n\t"\
	"vfmsub213pd	0x20(%%r15),%%ymm2,%%ymm1			\n\t			vmulpd 0x1a0(%%rsi),%%ymm13,%%ymm13	\n\t"\
	"vfmadd213pd	0x20(%%r14),%%ymm2,%%ymm3			\n\t			vmulpd 0x200(%%rsi),%%ymm9 ,%%ymm9 		\n\t"\
		"												\n\t			vaddpd	0x20(%%r11),%%ymm13,%%ymm13			\n\t"\
		"												\n\t			vaddpd	%%ymm15,%%ymm9 ,%%ymm9 			\n\t"\
		"												\n\t			vmulpd	0xe0(%%rsi),%%ymm13,%%ymm13	\n\t"\
		"vmovaps 0x40(%%rsi),%%ymm2						\n\t			vmulpd	0xe0(%%rsi),%%ymm9 ,%%ymm9 		\n\t"\
	"vfmadd213pd	0x20(%%r14),%%ymm0,%%ymm4			\n\t			vaddpd	%%ymm12,%%ymm13,%%ymm13			\n\t"\
	"vfmsub213pd	0x20(%%r13),%%ymm0,%%ymm1			\n\t			vaddpd	%%ymm8 ,%%ymm9 ,%%ymm9 			\n\t"\
	"vfmsub213pd	0x20(%%r15),%%ymm0,%%ymm3			\n\t			vmovaps	0x20(%%r10),%%ymm12			\n\t"\
	"vmulpd	-0x20(%%rsi),%%ymm2,%%ymm0					\n\t			vmovaps	%%ymm11,%%ymm15			\n\t"\
	"vfnmadd231pd	%%ymm4,%%ymm2,%%ymm6				\n\t			vmovaps	%%ymm12,%%ymm8 			\n\t"\
	"vfnmadd231pd	%%ymm1,%%ymm2,%%ymm7				\n\t			vaddpd	%%ymm14,%%ymm15,%%ymm15			\n\t"\
	"vfnmadd231pd	%%ymm3,%%ymm2,%%ymm5				\n\t			vaddpd	%%ymm13,%%ymm8 ,%%ymm8 			\n\t"\
		"												\n\t			vsubpd	%%ymm11,%%ymm14,%%ymm14			\n\t"\
		"												\n\t			vsubpd	%%ymm12,%%ymm13,%%ymm13			\n\t"\
		"												\n\t			vaddpd	%%ymm10,%%ymm14,%%ymm14			\n\t"\
	"vfmadd213pd	%%ymm6,%%ymm0,%%ymm4				\n\t			vaddpd	%%ymm9 ,%%ymm13,%%ymm13			\n\t"\
	"vfmadd213pd	%%ymm7,%%ymm0,%%ymm1				\n\t			vaddpd	%%ymm11,%%ymm10,%%ymm10			\n\t"\
	"vfmadd213pd	%%ymm5,%%ymm0,%%ymm3				\n\t			vaddpd	%%ymm9 ,%%ymm12,%%ymm12			\n\t"\
	/* yi-data in ymm8,10,12,13,14,15; ymm0,2,9,11 free: */\
/*** o1,2|A-C already in r10,11|13-15: ***/\
	/*** oA already in r13; Compute o3,6,7 into r12,8,9: ***/\
	"movslq	0x0c(%%rdx),%%r12	\n\t	addq	%%rcx,%%r12	\n\t"/* o3 */\
	"movslq	0x18(%%rdx),%%r8	\n\t	addq	%%rcx,%%r8	\n\t"/* o6 */\
	"movslq	0x1c(%%rdx),%%r9	\n\t	addq	%%rcx,%%r9	\n\t"/* o7 */\
		"vmovaps	%%ymm4 ,%%ymm0						\n\t			vmovaps	%%ymm7 ,%%ymm9		\n\t"\
		"vsubpd	%%ymm15,%%ymm4,%%ymm4					\n\t			vsubpd	%%ymm8 ,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	%%ymm15,%%ymm0,%%ymm0					\n\t			vaddpd	%%ymm8 ,%%ymm9,%%ymm9		\n\t"\
		"vmovaps	%%ymm4 ,0x20(%%r9 )					\n\t			vmovaps	%%ymm7 ,0x20(%%r13)	\n\t"\
		"vmovaps	%%ymm0 ,0x20(%%r8 )					\n\t			vmovaps	%%ymm9 ,0x20(%%r12)	\n\t"\
	/* Compute o5,8,9,4, store into rsi,rdi,r12-13: */\
	"movslq	0x14(%%rdx),%%rsi	\n\t	addq	%%rcx,%%rsi	\n\t"/* o5 */\
	"movslq	0x20(%%rdx),%%rdi	\n\t	addq	%%rcx,%%rdi	\n\t"/* o8 */\
	"movslq	0x24(%%rdx),%%r13	\n\t	addq	%%rcx,%%r13	\n\t"/* o9 */\
	"movslq	0x10(%%rdx),%%r12	\n\t	addq	%%rcx,%%r12	\n\t"/* o4 */\
		"vmovaps	%%ymm5 ,%%ymm0						\n\t			vmovaps	%%ymm6 ,%%ymm9		\n\t"\
		"vsubpd	%%ymm14,%%ymm5,%%ymm5					\n\t			vsubpd	%%ymm13,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	%%ymm14,%%ymm0,%%ymm0					\n\t			vaddpd	%%ymm13,%%ymm9,%%ymm9		\n\t"\
		"vmovaps	%%ymm5 ,0x20(%%rsi)					\n\t			vmovaps	%%ymm6 ,0x20(%%r13)	\n\t"\
		"vmovaps	%%ymm0 ,0x20(%%rdi)					\n\t			vmovaps	%%ymm9 ,0x20(%%r12)	\n\t"\
	/*** o1,2,B,C already in r10,11,14,15: ***/\
		"vmovaps	%%ymm1 ,%%ymm0						\n\t			vmovaps	%%ymm3 ,%%ymm9		\n\t"\
		"vsubpd	%%ymm10,%%ymm1,%%ymm1					\n\t			vsubpd	%%ymm12,%%ymm3,%%ymm3		\n\t"\
		"vaddpd	%%ymm10,%%ymm0,%%ymm0					\n\t			vaddpd	%%ymm12,%%ymm9,%%ymm9		\n\t"\
		"vmovaps	%%ymm1 ,0x20(%%r14)					\n\t			vmovaps	%%ymm3 ,0x20(%%r15)	\n\t"\
		"vmovaps	%%ymm0 ,0x20(%%r11)					\n\t			vmovaps	%%ymm9 ,0x20(%%r10)	\n\t"\
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

   #endif	// ? DFT_13_FMA

#elif defined(USE_AVX)

	// Opcount: [164 ADD, 74 MUL, 170 memref]
	//
	#define SSE2_RADIX_13_DFT(Xin0,Xi_off, Xcc, Xout0,Xo_off)\
	{\
	__asm__ volatile (\
	"movq	%[in0],%%rax		\n\t	movq	%[out0],%%rcx	\n\t"/* rax,rcx hold to I/O base-pointers */\
	"movq	%[i_off],%%rbx		\n\t	movq	%[o_off],%%rdx	\n\t"/* rbx,rdx hold pointers to uint32[] I/O pointer-offset-stride arrays */\
	"movq	%[__cc],%%rsi			\n\t"\
	"movslq	(%%rbx),%%r8	\n\t	addq	%%rax,%%r8	\n\t"/* I0 */\
	"movslq	(%%rdx),%%r9	\n\t	addq	%%rcx,%%r9	\n\t"/* O0 */\
	/* xr-terms:                                                      yi-terms: */\
	"movslq	0x04(%%rbx),%%r15	\n\t"/* Compute i1-6 from offsets: */\
	"movslq	0x08(%%rbx),%%r14	\n\t"\
	"movslq	0x0c(%%rbx),%%r13	\n\t"\
	"movslq	0x10(%%rbx),%%r12	\n\t"\
	"movslq	0x14(%%rbx),%%r11	\n\t"\
	"movslq	0x18(%%rbx),%%r10	\n\t"\
		"addq	%%rax,%%r15	\n\t	vmovaps (%%r15),%%ymm1	\n\t	vmovaps	0x20(%%r15),%%ymm9	\n\t"/* load A1-6 */\
		"addq	%%rax,%%r14	\n\t	vmovaps (%%r14),%%ymm3	\n\t	vmovaps	0x20(%%r14),%%ymm10	\n\t"\
		"addq	%%rax,%%r13	\n\t	vmovaps (%%r13),%%ymm6	\n\t	vmovaps	0x20(%%r13),%%ymm13	\n\t"\
		"addq	%%rax,%%r12	\n\t	vmovaps (%%r12),%%ymm4	\n\t	vmovaps	0x20(%%r12),%%ymm11	\n\t"\
		"addq	%%rax,%%r11	\n\t	vmovaps (%%r11),%%ymm5	\n\t	vmovaps	0x20(%%r11),%%ymm12	\n\t"\
		"addq	%%rax,%%r10	\n\t	vmovaps (%%r10),%%ymm7	\n\t	vmovaps	0x20(%%r10),%%ymm14	\n\t"\
		"vmovaps	-0x020(%%rsi),%%ymm0	\n\t"/* 2.0 */\
	"movslq	0x30(%%rbx),%%r15	\n\t"/* Compute iC-7 from offsets: */\
	"movslq	0x2c(%%rbx),%%r14	\n\t"\
	"movslq	0x28(%%rbx),%%r13	\n\t"\
	"movslq	0x24(%%rbx),%%r12	\n\t"\
	"movslq	0x20(%%rbx),%%r11	\n\t"\
	"movslq	0x1c(%%rbx),%%r10	\n\t"\
		"addq	%%rax,%%r15	\n\t	vaddpd	(%%r15),%%ymm1,%%ymm1	\n\t	vsubpd	0x20(%%r15),%%ymm9 ,%%ymm9 	\n\t"\
		"addq	%%rax,%%r14	\n\t	vaddpd	(%%r14),%%ymm3,%%ymm3	\n\t	vsubpd	0x20(%%r14),%%ymm10,%%ymm10	\n\t"\
		"addq	%%rax,%%r13	\n\t	vaddpd	(%%r13),%%ymm6,%%ymm6	\n\t	vsubpd	0x20(%%r13),%%ymm13,%%ymm13	\n\t"\
		"addq	%%rax,%%r12	\n\t	vaddpd	(%%r12),%%ymm4,%%ymm4	\n\t	vsubpd	0x20(%%r12),%%ymm11,%%ymm11	\n\t"\
		"addq	%%rax,%%r11	\n\t	vaddpd	(%%r11),%%ymm5,%%ymm5	\n\t	vsubpd	0x20(%%r11),%%ymm12,%%ymm12	\n\t"\
		"addq	%%rax,%%r10	\n\t	vaddpd	(%%r10),%%ymm7,%%ymm7	\n\t	vsubpd	0x20(%%r10),%%ymm14,%%ymm14	\n\t"\
		"vsubpd	%%ymm5,%%ymm1,%%ymm1					\n\t			vmovaps	%%ymm9 ,%%ymm15					\n\t"\
		"vmulpd	%%ymm0,%%ymm5,%%ymm5					\n\t			vmovaps	%%ymm10,%%ymm8 		 			\n\t"\
		"vsubpd	%%ymm6,%%ymm3,%%ymm3					\n\t			vsubpd	%%ymm11,%%ymm15,%%ymm15			\n\t"\
		"vmulpd	%%ymm0,%%ymm6,%%ymm6					\n\t			vaddpd	%%ymm12,%%ymm8 ,%%ymm8 			\n\t"\
		"vsubpd	%%ymm7,%%ymm4,%%ymm4					\n\t			vaddpd	%%ymm13,%%ymm15,%%ymm15			\n\t"\
		"vmulpd	%%ymm0,%%ymm7,%%ymm7					\n\t			vaddpd	%%ymm14,%%ymm8 ,%%ymm8 			\n\t"\
		"vaddpd	%%ymm1,%%ymm5,%%ymm5					\n\t			vaddpd	%%ymm11,%%ymm9 ,%%ymm9 			\n\t"\
		"vaddpd	%%ymm3,%%ymm6,%%ymm6					\n\t			vaddpd	%%ymm11,%%ymm13,%%ymm13			\n\t"\
		"vaddpd	%%ymm4,%%ymm7,%%ymm7					\n\t			vsubpd	%%ymm12,%%ymm10,%%ymm10			\n\t"\
		"vmovaps	%%ymm5,%%ymm2						\n\t			vsubpd	%%ymm12,%%ymm14,%%ymm14			\n\t"\
		"vaddpd	%%ymm6,%%ymm2,%%ymm2					\n\t			vmovaps	%%ymm8 ,%%ymm12					\n\t"\
		"vaddpd	%%ymm7,%%ymm2,%%ymm2					\n\t			vmovaps	%%ymm15,%%ymm11					\n\t"\
		"vmovaps	(%%r8 ),%%ymm0						\n\t			vmulpd	0x80(%%rsi),%%ymm8 ,%%ymm8 	\n\t"\
		"vaddpd	%%ymm2,%%ymm0,%%ymm0					\n\t			vmulpd	0x80(%%rsi),%%ymm15,%%ymm15	\n\t"\
		"vmulpd (%%rsi),%%ymm2,%%ymm2					\n\t			vmulpd	0x60(%%rsi),%%ymm12,%%ymm12	\n\t"\
		"vmovaps	%%ymm0,(%%r9 )						\n\t			vmulpd	0x60(%%rsi),%%ymm11,%%ymm11	\n\t"\
		"vaddpd	%%ymm2,%%ymm0,%%ymm0					\n\t			vsubpd	%%ymm15,%%ymm12,%%ymm12			\n\t"\
		"vmovaps	0x20(%%rsi),%%ymm2					\n\t			vaddpd	%%ymm8 ,%%ymm11,%%ymm11			\n\t"\
/* lcol|rcol section uses o-offsets 1,2,3|B,C - prestore those into r10-12|14,15: */\
	"movslq	0x04(%%rdx),%%r10	\n\t	addq	%%rcx,%%r10	\n\t"/* o1 */\
	"movslq	0x08(%%rdx),%%r11	\n\t	addq	%%rcx,%%r11	\n\t"/* o2 */\
	"movslq	0x0c(%%rdx),%%r12	\n\t	addq	%%rcx,%%r12	\n\t"/* o3 */\
	"movslq	0x2c(%%rdx),%%r14	\n\t	addq	%%rcx,%%r14	\n\t"/* oB */\
	"movslq	0x30(%%rdx),%%r15	\n\t	addq	%%rcx,%%r15	\n\t"/* oC */\
		"vmulpd	%%ymm2,%%ymm6,%%ymm6					\n\t			vmovaps	%%ymm12,(%%r15)			\n\t"\
		"vmulpd	%%ymm2,%%ymm5,%%ymm5					\n\t			vmovaps	%%ymm9 ,%%ymm8 			\n\t"\
		"vmulpd	%%ymm2,%%ymm7,%%ymm7					\n\t			vmovaps	%%ymm13,%%ymm12			\n\t"\
		"vmovaps %%ymm6,(%%r10)							\n\t			vaddpd	%%ymm10,%%ymm8 ,%%ymm8 			\n\t"\
		"vmovaps %%ymm5,(%%r11)							\n\t			vaddpd	%%ymm14,%%ymm12,%%ymm12			\n\t"\
		"vmovaps %%ymm7,(%%r12)							\n\t			vmovaps	%%ymm8 ,(%%r14)			\n\t"\
		"vmovaps 0x100(%%rsi),%%ymm2					\n\t			vmovaps	%%ymm12,%%ymm15			\n\t"\
		"vmulpd	%%ymm2,%%ymm5,%%ymm5					\n\t			vmulpd 0x160(%%rsi),%%ymm12,%%ymm12	\n\t"\
		"vmulpd	%%ymm2,%%ymm7,%%ymm7					\n\t			vmulpd 0x1c0(%%rsi),%%ymm8 ,%%ymm8 		\n\t"\
		"vmulpd	%%ymm2,%%ymm6,%%ymm6					\n\t			vaddpd	(%%r14),%%ymm12,%%ymm12			\n\t"\
		"vaddpd	%%ymm0,%%ymm5,%%ymm5					\n\t			vaddpd	%%ymm15,%%ymm8 ,%%ymm8 			\n\t"\
		"vaddpd	%%ymm0,%%ymm7,%%ymm7					\n\t			vmulpd	0xa0(%%rsi),%%ymm12,%%ymm12	\n\t"\
		"vaddpd	%%ymm0,%%ymm6,%%ymm6					\n\t			vmulpd	0xa0(%%rsi),%%ymm8 ,%%ymm8 		\n\t"\
		"vaddpd (%%r10),%%ymm5,%%ymm5					\n\t			vmovaps	%%ymm10,(%%r14)			\n\t"\
		"vaddpd (%%r11),%%ymm7,%%ymm7					\n\t			vmovaps	%%ymm14,%%ymm15			\n\t"\
		"vaddpd (%%r12),%%ymm6,%%ymm6					\n\t			vmulpd 0x180(%%rsi),%%ymm14,%%ymm14	\n\t"\
		"vmovaps 0x140(%%rsi),%%ymm2					\n\t			vmulpd 0x1e0(%%rsi),%%ymm10,%%ymm10	\n\t"\
		"vmovaps 0x120(%%rsi),%%ymm0					\n\t			vaddpd	(%%r14),%%ymm14,%%ymm14			\n\t"\
		"vmovaps %%ymm4,(%%r10)							\n\t			vaddpd	%%ymm15,%%ymm10,%%ymm10			\n\t"\
		"vmovaps %%ymm1,(%%r11)							\n\t			vmulpd	0xc0(%%rsi),%%ymm14,%%ymm14	\n\t"\
		"vmovaps %%ymm3,(%%r12)							\n\t			vmulpd	0xc0(%%rsi),%%ymm10,%%ymm10	\n\t"\
		"vmulpd	%%ymm2,%%ymm4,%%ymm4					\n\t			vaddpd	%%ymm12,%%ymm14,%%ymm14			\n\t"\
		"vmulpd	%%ymm2,%%ymm1,%%ymm1					\n\t			vaddpd	%%ymm8 ,%%ymm10,%%ymm10			\n\t"\
		"vmulpd	%%ymm2,%%ymm3,%%ymm3					\n\t			vmovaps	%%ymm9 ,(%%r14)			\n\t"\
		"vaddpd (%%r12),%%ymm4,%%ymm4					\n\t			vmovaps	%%ymm13,%%ymm15			\n\t"\
		"vsubpd (%%r10),%%ymm1,%%ymm1					\n\t			vmulpd 0x1a0(%%rsi),%%ymm13,%%ymm13	\n\t"\
		"vaddpd (%%r11),%%ymm3,%%ymm3					\n\t			vmulpd 0x200(%%rsi),%%ymm9 ,%%ymm9 		\n\t"\
		"vmulpd	%%ymm0,%%ymm4,%%ymm4					\n\t			vaddpd	(%%r14),%%ymm13,%%ymm13			\n\t"\
		"vmulpd	%%ymm0,%%ymm1,%%ymm1					\n\t			vaddpd	%%ymm15,%%ymm9 ,%%ymm9 			\n\t"\
		"vmulpd	%%ymm0,%%ymm3,%%ymm3					\n\t			vmulpd	0xe0(%%rsi),%%ymm13,%%ymm13	\n\t"\
		"vmovaps 0x40(%%rsi),%%ymm2						\n\t			vmulpd	0xe0(%%rsi),%%ymm9 ,%%ymm9 		\n\t"\
		"vaddpd (%%r11),%%ymm4,%%ymm4					\n\t			vaddpd	%%ymm12,%%ymm13,%%ymm13			\n\t"\
		"vsubpd (%%r12),%%ymm1,%%ymm1					\n\t			vaddpd	%%ymm8 ,%%ymm9 ,%%ymm9 			\n\t"\
		"vsubpd (%%r10),%%ymm3,%%ymm3					\n\t			vmovaps	(%%r15),%%ymm12			\n\t"\
		"vmulpd	%%ymm2,%%ymm4,%%ymm4					\n\t			vmovaps	%%ymm11,%%ymm15			\n\t"\
		"vmulpd	%%ymm2,%%ymm1,%%ymm1					\n\t			vmovaps	%%ymm12,%%ymm8 			\n\t"\
		"vmulpd	%%ymm2,%%ymm3,%%ymm3					\n\t			vaddpd	%%ymm14,%%ymm15,%%ymm15			\n\t"\
		"vmovaps -0x20(%%rsi),%%ymm0					\n\t			vaddpd	%%ymm13,%%ymm8 ,%%ymm8 			\n\t"\
		"vsubpd	%%ymm4,%%ymm6,%%ymm6					\n\t			vsubpd	%%ymm11,%%ymm14,%%ymm14			\n\t"\
		"vmulpd	%%ymm0,%%ymm4,%%ymm4					\n\t			vsubpd	%%ymm12,%%ymm13,%%ymm13			\n\t"\
		"vsubpd	%%ymm1,%%ymm7,%%ymm7					\n\t			vaddpd	%%ymm10,%%ymm14,%%ymm14			\n\t"\
		"vmulpd	%%ymm0,%%ymm1,%%ymm1					\n\t			vaddpd	%%ymm9 ,%%ymm13,%%ymm13			\n\t"\
		"vsubpd	%%ymm3,%%ymm5,%%ymm5					\n\t			vaddpd	%%ymm11,%%ymm10,%%ymm10			\n\t"\
		"vmulpd	%%ymm0,%%ymm3,%%ymm3					\n\t			vaddpd	%%ymm9 ,%%ymm12,%%ymm12			\n\t"\
		"vaddpd	%%ymm6,%%ymm4,%%ymm4					\n\t"\
		"vaddpd	%%ymm7,%%ymm1,%%ymm1					\n\t"\
		"vaddpd	%%ymm5,%%ymm3,%%ymm3					\n\t"\
	/* yi-data in ymm8,10,12,13,14,15; ymm0,2,9,11 free: */\
	"movslq	0x18(%%rdx),%%r10	\n\t	addq	%%rcx,%%r10	\n\t"/* o6 */\
	"movslq	0x1c(%%rdx),%%r11	\n\t	addq	%%rcx,%%r11	\n\t"/* o7 */\
	"movslq	0x0c(%%rdx),%%r12	\n\t	addq	%%rcx,%%r12	\n\t"/* o3 */\
	"movslq	0x28(%%rdx),%%r13	\n\t	addq	%%rcx,%%r13	\n\t"/* oA */\
		"vmovaps	%%ymm4 ,%%ymm0						\n\t			vmovaps	%%ymm7 ,%%ymm9 			\n\t"\
		"vsubpd	%%ymm15,%%ymm4,%%ymm4					\n\t			vsubpd	%%ymm8 ,%%ymm7,%%ymm7	\n\t"\
		"vaddpd	%%ymm15,%%ymm0,%%ymm0					\n\t			vaddpd	%%ymm8 ,%%ymm9,%%ymm9	\n\t"\
		"vmovaps	%%ymm4 ,(%%r10)						\n\t			vmovaps	%%ymm7 ,(%%r12)	\n\t"\
		"vmovaps	%%ymm0 ,(%%r11)						\n\t			vmovaps	%%ymm9 ,(%%r13)	\n\t"\
	"movslq	0x20(%%rdx),%%r10	\n\t	addq	%%rcx,%%r10	\n\t"/* o8 */\
	"movslq	0x14(%%rdx),%%r11	\n\t	addq	%%rcx,%%r11	\n\t"/* o5 */\
	"movslq	0x10(%%rdx),%%r12	\n\t	addq	%%rcx,%%r12	\n\t"/* o4 */\
	"movslq	0x24(%%rdx),%%r13	\n\t	addq	%%rcx,%%r13	\n\t"/* o9 */\
		"vmovaps	%%ymm5 ,%%ymm0						\n\t			vmovaps	%%ymm6 ,%%ymm9 			\n\t"\
		"vsubpd	%%ymm14,%%ymm5,%%ymm5					\n\t			vsubpd	%%ymm13,%%ymm6,%%ymm6	\n\t"\
		"vaddpd	%%ymm14,%%ymm0,%%ymm0					\n\t			vaddpd	%%ymm13,%%ymm9,%%ymm9	\n\t"\
		"vmovaps	%%ymm5 ,(%%r10)						\n\t			vmovaps	%%ymm6 ,(%%r12)	\n\t"\
		"vmovaps	%%ymm0 ,(%%r11)						\n\t			vmovaps	%%ymm9 ,(%%r13)	\n\t"\
	"movslq	0x08(%%rdx),%%r10	\n\t	addq	%%rcx,%%r10	\n\t"/* o2 */\
	"movslq	0x2c(%%rdx),%%r11	\n\t	addq	%%rcx,%%r11	\n\t"/* oB */\
	"movslq	0x04(%%rdx),%%r12	\n\t	addq	%%rcx,%%r12	\n\t"/* o1 */\
	"movslq	0x30(%%rdx),%%r13	\n\t	addq	%%rcx,%%r13	\n\t"/* oC */\
		"vmovaps	%%ymm1 ,%%ymm0						\n\t			vmovaps	%%ymm3 ,%%ymm9 			\n\t"\
		"vsubpd	%%ymm10,%%ymm1,%%ymm1					\n\t			vsubpd	%%ymm12,%%ymm3,%%ymm3	\n\t"\
		"vaddpd	%%ymm10,%%ymm0,%%ymm0					\n\t			vaddpd	%%ymm12,%%ymm9,%%ymm9	\n\t"\
		"vmovaps	%%ymm1 ,(%%r10)						\n\t			vmovaps	%%ymm3 ,(%%r12)	\n\t"\
		"vmovaps	%%ymm0 ,(%%r11)						\n\t			vmovaps	%%ymm9 ,(%%r13)	\n\t"\
		/****************************************************************/\
		/*                      IMAG PARTS:                             */\
		/****************************************************************/\
		"addq	$0x20,%%r8 		 		\n\t"/* I0.imag */\
				"addq	$0x20,%%r9 		\n\t"/* O0.imag */\
	/* xi-terms:                                                      yr-terms: */\
	"movslq	0x04(%%rbx),%%r15	\n\t"/* Compute i1-6 from offsets: */\
	"movslq	0x08(%%rbx),%%r14	\n\t"\
	"movslq	0x0c(%%rbx),%%r13	\n\t"\
	"movslq	0x10(%%rbx),%%r12	\n\t"\
	"movslq	0x14(%%rbx),%%r11	\n\t"\
	"movslq	0x18(%%rbx),%%r10	\n\t"\
		"addq	%%rax,%%r15	\n\t	vmovaps 0x20(%%r15),%%ymm1	\n\t	vmovaps	(%%r15),%%ymm9	\n\t"/* load A1-6 */\
		"addq	%%rax,%%r14	\n\t	vmovaps 0x20(%%r14),%%ymm3	\n\t	vmovaps	(%%r14),%%ymm10	\n\t"\
		"addq	%%rax,%%r13	\n\t	vmovaps 0x20(%%r13),%%ymm6	\n\t	vmovaps	(%%r13),%%ymm13	\n\t"\
		"addq	%%rax,%%r12	\n\t	vmovaps 0x20(%%r12),%%ymm4	\n\t	vmovaps	(%%r12),%%ymm11	\n\t"\
		"addq	%%rax,%%r11	\n\t	vmovaps 0x20(%%r11),%%ymm5	\n\t	vmovaps	(%%r11),%%ymm12	\n\t"\
		"addq	%%rax,%%r10	\n\t	vmovaps 0x20(%%r10),%%ymm7	\n\t	vmovaps	(%%r10),%%ymm14	\n\t"\
		"vmovaps	-0x020(%%rsi),%%ymm0	\n\t"/* 2.0 */\
	"movslq	0x30(%%rbx),%%r15	\n\t"/* Compute iC-7 from offsets: */\
	"movslq	0x2c(%%rbx),%%r14	\n\t"\
	"movslq	0x28(%%rbx),%%r13	\n\t"\
	"movslq	0x24(%%rbx),%%r12	\n\t"\
	"movslq	0x20(%%rbx),%%r11	\n\t"\
	"movslq	0x1c(%%rbx),%%r10	\n\t"\
		"vmovaps -0x20(%%rsi),%%ymm0					\n\t"\
		"addq	%%rax,%%r15	\n\t	vaddpd	0x20(%%r15),%%ymm1,%%ymm1	\n\t	vsubpd	(%%r15),%%ymm9 ,%%ymm9 	\n\t"\
		"addq	%%rax,%%r14	\n\t	vaddpd	0x20(%%r14),%%ymm3,%%ymm3	\n\t	vsubpd	(%%r14),%%ymm10,%%ymm10	\n\t"\
		"addq	%%rax,%%r13	\n\t	vaddpd	0x20(%%r13),%%ymm6,%%ymm6	\n\t	vsubpd	(%%r13),%%ymm13,%%ymm13	\n\t"\
		"addq	%%rax,%%r12	\n\t	vaddpd	0x20(%%r12),%%ymm4,%%ymm4	\n\t	vsubpd	(%%r12),%%ymm11,%%ymm11	\n\t"\
		"addq	%%rax,%%r11	\n\t	vaddpd	0x20(%%r11),%%ymm5,%%ymm5	\n\t	vsubpd	(%%r11),%%ymm12,%%ymm12	\n\t"\
		"addq	%%rax,%%r10	\n\t	vaddpd	0x20(%%r10),%%ymm7,%%ymm7	\n\t	vsubpd	(%%r10),%%ymm14,%%ymm14	\n\t"\
		"vsubpd	%%ymm5,%%ymm1,%%ymm1					\n\t			vmovaps	%%ymm9 ,%%ymm15			\n\t"\
		"vmulpd	%%ymm0,%%ymm5,%%ymm5					\n\t			vmovaps	%%ymm10,%%ymm8 			\n\t"\
		"vsubpd	%%ymm6,%%ymm3,%%ymm3					\n\t			vsubpd	%%ymm11,%%ymm15,%%ymm15			\n\t"\
		"vmulpd	%%ymm0,%%ymm6,%%ymm6					\n\t			vaddpd	%%ymm12,%%ymm8 ,%%ymm8 			\n\t"\
		"vsubpd	%%ymm7,%%ymm4,%%ymm4					\n\t			vaddpd	%%ymm13,%%ymm15,%%ymm15			\n\t"\
		"vmulpd	%%ymm0,%%ymm7,%%ymm7					\n\t			vaddpd	%%ymm14,%%ymm8 ,%%ymm8 			\n\t"\
		"vaddpd	%%ymm1,%%ymm5,%%ymm5					\n\t			vaddpd	%%ymm11,%%ymm9 ,%%ymm9 			\n\t"\
		"vaddpd	%%ymm3,%%ymm6,%%ymm6					\n\t			vaddpd	%%ymm11,%%ymm13,%%ymm13			\n\t"\
		"vaddpd	%%ymm4,%%ymm7,%%ymm7					\n\t			vsubpd	%%ymm12,%%ymm10,%%ymm10			\n\t"\
		"vmovaps	%%ymm5,%%ymm2						\n\t			vsubpd	%%ymm12,%%ymm14,%%ymm14			\n\t"\
		"vaddpd	%%ymm6,%%ymm2,%%ymm2					\n\t			vmovaps	%%ymm8 ,%%ymm12			\n\t"\
		"vaddpd	%%ymm7,%%ymm2,%%ymm2					\n\t			vmovaps	%%ymm15,%%ymm11			\n\t"\
		"vmovaps	(%%r8 ),%%ymm0						\n\t			vmulpd	0x80(%%rsi),%%ymm8 ,%%ymm8 	\n\t"\
		"vaddpd	%%ymm2,%%ymm0,%%ymm0					\n\t			vmulpd	0x80(%%rsi),%%ymm15,%%ymm15	\n\t"\
		"vmulpd	(%%rsi),%%ymm2,%%ymm2					\n\t			vmulpd	0x60(%%rsi),%%ymm12,%%ymm12	\n\t"\
		"vmovaps	%%ymm0,(%%r9 )						\n\t			vmulpd	0x60(%%rsi),%%ymm11,%%ymm11	\n\t"\
		"vaddpd	%%ymm2,%%ymm0,%%ymm0					\n\t			vsubpd	%%ymm15,%%ymm12,%%ymm12			\n\t"\
		"vmovaps	0x20(%%rsi),%%ymm2					\n\t			vaddpd	%%ymm8 ,%%ymm11,%%ymm11			\n\t"\
/* rcol|lcol (note LR switcherooski) section uses o-offsets 1,2|A,B,C - prestore those into r10,11|13-15: */\
	"movslq	0x04(%%rdx),%%r10	\n\t	addq	%%rcx,%%r10	\n\t"/* o1 */\
	"movslq	0x08(%%rdx),%%r11	\n\t	addq	%%rcx,%%r11	\n\t"/* o2 */\
	"movslq	0x28(%%rdx),%%r13	\n\t	addq	%%rcx,%%r13	\n\t"/* oA */\
	"movslq	0x2c(%%rdx),%%r14	\n\t	addq	%%rcx,%%r14	\n\t"/* oB */\
	"movslq	0x30(%%rdx),%%r15	\n\t	addq	%%rcx,%%r15	\n\t"/* oC */\
		"vmulpd	%%ymm2,%%ymm6,%%ymm6					\n\t			vmovaps	%%ymm12,0x20(%%r10)	\n\t"\
		"vmulpd	%%ymm2,%%ymm5,%%ymm5					\n\t			vmovaps	%%ymm9 ,%%ymm8 			\n\t"\
		"vmulpd	%%ymm2,%%ymm7,%%ymm7					\n\t			vmovaps	%%ymm13,%%ymm12			\n\t"\
		"vmovaps %%ymm6,0x20(%%r15)						\n\t			vaddpd	%%ymm10,%%ymm8 ,%%ymm8 			\n\t"\
		"vmovaps %%ymm5,0x20(%%r14)						\n\t			vaddpd	%%ymm14,%%ymm12,%%ymm12			\n\t"\
		"vmovaps %%ymm7,0x20(%%r13)						\n\t			vmovaps	%%ymm8 ,0x20(%%r11)	\n\t"\
		"vmovaps 0x100(%%rsi),%%ymm2					\n\t			vmovaps	%%ymm12,%%ymm15			\n\t"\
		"vmulpd	%%ymm2,%%ymm5,%%ymm5					\n\t			vmulpd 0x160(%%rsi),%%ymm12,%%ymm12	\n\t"\
		"vmulpd	%%ymm2,%%ymm7,%%ymm7					\n\t			vmulpd 0x1c0(%%rsi),%%ymm8 ,%%ymm8 		\n\t"\
		"vmulpd	%%ymm2,%%ymm6,%%ymm6					\n\t			vaddpd	0x20(%%r11),%%ymm12,%%ymm12	\n\t"\
		"vaddpd	%%ymm0,%%ymm5,%%ymm5					\n\t			vaddpd	%%ymm15,%%ymm8 ,%%ymm8 			\n\t"\
		"vaddpd	%%ymm0,%%ymm7,%%ymm7					\n\t			vmulpd	0xa0(%%rsi),%%ymm12,%%ymm12	\n\t"\
		"vaddpd	%%ymm0,%%ymm6,%%ymm6					\n\t			vmulpd	0xa0(%%rsi),%%ymm8 ,%%ymm8 		\n\t"\
		"vaddpd	0x20(%%r15),%%ymm5,%%ymm5				\n\t			vmovaps	%%ymm10,0x20(%%r11)	\n\t"\
		"vaddpd	0x20(%%r14),%%ymm7,%%ymm7				\n\t			vmovaps	%%ymm14,%%ymm15			\n\t"\
		"vaddpd	0x20(%%r13),%%ymm6,%%ymm6				\n\t			vmulpd 0x180(%%rsi),%%ymm14,%%ymm14	\n\t"\
		"vmovaps 0x140(%%rsi),%%ymm2					\n\t			vmulpd 0x1e0(%%rsi),%%ymm10,%%ymm10	\n\t"\
		"vmovaps 0x120(%%rsi),%%ymm0					\n\t			vaddpd	0x20(%%r11),%%ymm14,%%ymm14	\n\t"\
		"vmovaps %%ymm4,0x20(%%r15)						\n\t			vaddpd	%%ymm15,%%ymm10,%%ymm10			\n\t"\
		"vmovaps %%ymm1,0x20(%%r14)						\n\t			vmulpd	0xc0(%%rsi),%%ymm14,%%ymm14	\n\t"\
		"vmovaps %%ymm3,0x20(%%r13)						\n\t			vmulpd	0xc0(%%rsi),%%ymm10,%%ymm10	\n\t"\
		"vmulpd	%%ymm2,%%ymm4,%%ymm4					\n\t			vaddpd	%%ymm12,%%ymm14,%%ymm14			\n\t"\
		"vmulpd	%%ymm2,%%ymm1,%%ymm1					\n\t			vaddpd	%%ymm8 ,%%ymm10,%%ymm10			\n\t"\
		"vmulpd	%%ymm2,%%ymm3,%%ymm3					\n\t			vmovaps	%%ymm9 ,0x20(%%r11)	\n\t"\
		"vaddpd	0x20(%%r13),%%ymm4,%%ymm4				\n\t			vmovaps	%%ymm13,%%ymm15			\n\t"\
		"vsubpd	0x20(%%r15),%%ymm1,%%ymm1				\n\t			vmulpd 0x1a0(%%rsi),%%ymm13,%%ymm13	\n\t"\
		"vaddpd	0x20(%%r14),%%ymm3,%%ymm3				\n\t			vmulpd 0x200(%%rsi),%%ymm9 ,%%ymm9 		\n\t"\
		"vmulpd	%%ymm0,%%ymm4,%%ymm4					\n\t			vaddpd	0x20(%%r11),%%ymm13,%%ymm13	\n\t"\
		"vmulpd	%%ymm0,%%ymm1,%%ymm1					\n\t			vaddpd	%%ymm15,%%ymm9 ,%%ymm9 			\n\t"\
		"vmulpd	%%ymm0,%%ymm3,%%ymm3					\n\t			vmulpd	0xe0(%%rsi),%%ymm13,%%ymm13	\n\t"\
		"vmovaps 0x40(%%rsi),%%ymm2						\n\t			vmulpd	0xe0(%%rsi),%%ymm9 ,%%ymm9 		\n\t"\
		"vaddpd	0x20(%%r14),%%ymm4,%%ymm4				\n\t			vaddpd	%%ymm12,%%ymm13,%%ymm13			\n\t"\
		"vsubpd	0x20(%%r13),%%ymm1,%%ymm1				\n\t			vaddpd	%%ymm8 ,%%ymm9 ,%%ymm9 			\n\t"\
		"vsubpd	0x20(%%r15),%%ymm3,%%ymm3				\n\t			vmovaps	0x20(%%r10),%%ymm12	\n\t"\
		"vmulpd	%%ymm2,%%ymm4,%%ymm4					\n\t			vmovaps	%%ymm11,%%ymm15			\n\t"\
		"vmulpd	%%ymm2,%%ymm1,%%ymm1					\n\t			vmovaps	%%ymm12,%%ymm8 			\n\t"\
		"vmulpd	%%ymm2,%%ymm3,%%ymm3					\n\t			vaddpd	%%ymm14,%%ymm15,%%ymm15			\n\t"\
		"vmovaps -0x20(%%rsi),%%ymm0					\n\t			vaddpd	%%ymm13,%%ymm8 ,%%ymm8 			\n\t"\
		"vsubpd	%%ymm4,%%ymm6,%%ymm6					\n\t			vsubpd	%%ymm11,%%ymm14,%%ymm14			\n\t"\
		"vmulpd	%%ymm0,%%ymm4,%%ymm4					\n\t			vsubpd	%%ymm12,%%ymm13,%%ymm13			\n\t"\
		"vsubpd	%%ymm1,%%ymm7,%%ymm7					\n\t			vaddpd	%%ymm10,%%ymm14,%%ymm14			\n\t"\
		"vmulpd	%%ymm0,%%ymm1,%%ymm1					\n\t			vaddpd	%%ymm9 ,%%ymm13,%%ymm13			\n\t"\
		"vsubpd	%%ymm3,%%ymm5,%%ymm5					\n\t			vaddpd	%%ymm11,%%ymm10,%%ymm10			\n\t"\
		"vmulpd	%%ymm0,%%ymm3,%%ymm3					\n\t			vaddpd	%%ymm9 ,%%ymm12,%%ymm12			\n\t"\
		"vaddpd	%%ymm6,%%ymm4,%%ymm4					\n\t"\
		"vaddpd	%%ymm7,%%ymm1,%%ymm1					\n\t"\
		"vaddpd	%%ymm5,%%ymm3,%%ymm3					\n\t"\
	/* yi-data in ymm8,10,12,13,14,15; ymm0,2,9,11 free: */\
	"movslq	0x18(%%rdx),%%r10	\n\t	addq	%%rcx,%%r10	\n\t"/* o6 */\
	"movslq	0x1c(%%rdx),%%r11	\n\t	addq	%%rcx,%%r11	\n\t"/* o7 */\
	"movslq	0x0c(%%rdx),%%r12	\n\t	addq	%%rcx,%%r12	\n\t"/* o3 */\
	"movslq	0x28(%%rdx),%%r13	\n\t	addq	%%rcx,%%r13	\n\t"/* oA */\
		"vmovaps	%%ymm4 ,%%ymm0						\n\t			vmovaps	%%ymm7 ,%%ymm9		\n\t"\
		"vsubpd	%%ymm15,%%ymm4,%%ymm4					\n\t			vsubpd	%%ymm8 ,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	%%ymm15,%%ymm0,%%ymm0					\n\t			vaddpd	%%ymm8 ,%%ymm9,%%ymm9		\n\t"\
		"vmovaps	%%ymm4 ,0x20(%%r11)					\n\t			vmovaps	%%ymm7 ,0x20(%%r13)	\n\t"\
		"vmovaps	%%ymm0 ,0x20(%%r10)					\n\t			vmovaps	%%ymm9 ,0x20(%%r12)	\n\t"\
	"movslq	0x20(%%rdx),%%r10	\n\t	addq	%%rcx,%%r10	\n\t"/* o8 */\
	"movslq	0x14(%%rdx),%%r11	\n\t	addq	%%rcx,%%r11	\n\t"/* o5 */\
	"movslq	0x10(%%rdx),%%r12	\n\t	addq	%%rcx,%%r12	\n\t"/* o4 */\
	"movslq	0x24(%%rdx),%%r13	\n\t	addq	%%rcx,%%r13	\n\t"/* o9 */\
		"vmovaps	%%ymm5 ,%%ymm0						\n\t			vmovaps	%%ymm6 ,%%ymm9		\n\t"\
		"vsubpd	%%ymm14,%%ymm5,%%ymm5					\n\t			vsubpd	%%ymm13,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	%%ymm14,%%ymm0,%%ymm0					\n\t			vaddpd	%%ymm13,%%ymm9,%%ymm9		\n\t"\
		"vmovaps	%%ymm5 ,0x20(%%r11)					\n\t			vmovaps	%%ymm6 ,0x20(%%r13)	\n\t"\
		"vmovaps	%%ymm0 ,0x20(%%r10)					\n\t			vmovaps	%%ymm9 ,0x20(%%r12)	\n\t"\
	"movslq	0x08(%%rdx),%%r10	\n\t	addq	%%rcx,%%r10	\n\t"/* o2 */\
	"movslq	0x2c(%%rdx),%%r11	\n\t	addq	%%rcx,%%r11	\n\t"/* oB */\
	"movslq	0x04(%%rdx),%%r12	\n\t	addq	%%rcx,%%r12	\n\t"/* o1 */\
	"movslq	0x30(%%rdx),%%r13	\n\t	addq	%%rcx,%%r13	\n\t"/* oC */\
		"vmovaps	%%ymm1 ,%%ymm0						\n\t			vmovaps	%%ymm3 ,%%ymm9		\n\t"\
		"vsubpd	%%ymm10,%%ymm1,%%ymm1					\n\t			vsubpd	%%ymm12,%%ymm3,%%ymm3		\n\t"\
		"vaddpd	%%ymm10,%%ymm0,%%ymm0					\n\t			vaddpd	%%ymm12,%%ymm9,%%ymm9		\n\t"\
		"vmovaps	%%ymm1 ,0x20(%%r11)					\n\t			vmovaps	%%ymm3 ,0x20(%%r13)	\n\t"\
		"vmovaps	%%ymm0 ,0x20(%%r10)					\n\t			vmovaps	%%ymm9 ,0x20(%%r12)	\n\t"\
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

	#define SSE2_RADIX_13_DFT(Xin0,Xi_off, Xcc, Xout0,Xo_off)\
	{\
	__asm__ volatile (\
	"movq	%[in0],%%rax		\n\t	movq	%[out0],%%rcx	\n\t"/* rax,rcx hold to I/O base-pointers */\
	"movq	%[i_off],%%rbx		\n\t	movq	%[o_off],%%rdx	\n\t"/* rbx,rdx hold pointers to uint32[] I/O pointer-offset-stride arrays */\
	"movq	%[__cc],%%rsi			\n\t"\
	"movslq	(%%rbx),%%r8	\n\t	addq	%%rax,%%r8	\n\t"/* I0 */\
	"movslq	(%%rdx),%%r9	\n\t	addq	%%rcx,%%r9	\n\t"/* O0 */\
	/* xr-terms:                                                      yi-terms: */\
	"movslq	0x04(%%rbx),%%r15	\n\t"/* Compute i1-6 from offsets: */\
	"movslq	0x08(%%rbx),%%r14	\n\t"\
	"movslq	0x0c(%%rbx),%%r13	\n\t"\
	"movslq	0x10(%%rbx),%%r12	\n\t"\
	"movslq	0x14(%%rbx),%%r11	\n\t"\
	"movslq	0x18(%%rbx),%%r10	\n\t"\
		"addq	%%rax,%%r15	\n\t	movaps (%%r15),%%xmm1	\n\t	movaps	0x10(%%r15),%%xmm9	\n\t"/* load A1-6 */\
		"addq	%%rax,%%r14	\n\t	movaps (%%r14),%%xmm3	\n\t	movaps	0x10(%%r14),%%xmm10	\n\t"\
		"addq	%%rax,%%r13	\n\t	movaps (%%r13),%%xmm6	\n\t	movaps	0x10(%%r13),%%xmm13	\n\t"\
		"addq	%%rax,%%r12	\n\t	movaps (%%r12),%%xmm4	\n\t	movaps	0x10(%%r12),%%xmm11	\n\t"\
		"addq	%%rax,%%r11	\n\t	movaps (%%r11),%%xmm5	\n\t	movaps	0x10(%%r11),%%xmm12	\n\t"\
		"addq	%%rax,%%r10	\n\t	movaps (%%r10),%%xmm7	\n\t	movaps	0x10(%%r10),%%xmm14	\n\t"\
		"movaps	-0x010(%%rsi),%%xmm0	\n\t"/* 2.0 */\
	"movslq	0x30(%%rbx),%%r15	\n\t"/* Compute iC-7 from offsets: */\
	"movslq	0x2c(%%rbx),%%r14	\n\t"\
	"movslq	0x28(%%rbx),%%r13	\n\t"\
	"movslq	0x24(%%rbx),%%r12	\n\t"\
	"movslq	0x20(%%rbx),%%r11	\n\t"\
	"movslq	0x1c(%%rbx),%%r10	\n\t"\
		"addq	%%rax,%%r15	\n\t	addpd (%%r15),%%xmm1	\n\t	subpd	0x10(%%r15),%%xmm9	\n\t"/* A1-6r += AC-7r; A1-6i -= AC-7i */\
		"addq	%%rax,%%r14	\n\t	addpd (%%r14),%%xmm3	\n\t	subpd	0x10(%%r14),%%xmm10	\n\t"\
		"addq	%%rax,%%r13	\n\t	addpd (%%r13),%%xmm6	\n\t	subpd	0x10(%%r13),%%xmm13	\n\t"\
		"addq	%%rax,%%r12	\n\t	addpd (%%r12),%%xmm4	\n\t	subpd	0x10(%%r12),%%xmm11	\n\t"\
		"addq	%%rax,%%r11	\n\t	addpd (%%r11),%%xmm5	\n\t	subpd	0x10(%%r11),%%xmm12	\n\t"\
		"addq	%%rax,%%r10	\n\t	addpd (%%r10),%%xmm7	\n\t	subpd	0x10(%%r10),%%xmm14	\n\t"\
		"subpd	%%xmm5,%%xmm1							\n\t			movaps	%%xmm9 ,%%xmm15			\n\t"\
		"mulpd	%%xmm0,%%xmm5							\n\t			movaps	%%xmm10,%%xmm8			\n\t"\
		"subpd	%%xmm6,%%xmm3							\n\t			subpd	%%xmm11,%%xmm15			\n\t"\
		"mulpd	%%xmm0,%%xmm6							\n\t			addpd	%%xmm12,%%xmm8			\n\t"\
		"subpd	%%xmm7,%%xmm4							\n\t			addpd	%%xmm13,%%xmm15			\n\t"\
		"mulpd	%%xmm0,%%xmm7							\n\t			addpd	%%xmm14,%%xmm8			\n\t"\
		"addpd	%%xmm1,%%xmm5							\n\t			addpd	%%xmm11,%%xmm9			\n\t"\
		"addpd	%%xmm3,%%xmm6							\n\t			addpd	%%xmm11,%%xmm13			\n\t"\
		"addpd	%%xmm4,%%xmm7							\n\t			subpd	%%xmm12,%%xmm10			\n\t"\
		"movaps	%%xmm5,%%xmm2							\n\t			subpd	%%xmm12,%%xmm14			\n\t"\
		"addpd	%%xmm6,%%xmm2							\n\t			movaps	%%xmm8 ,%%xmm12			\n\t"\
		"addpd	%%xmm7,%%xmm2							\n\t			movaps	%%xmm15,%%xmm11			\n\t"\
		"movaps	(%%r8 ),%%xmm0							\n\t			mulpd	0x040(%%rsi),%%xmm8		\n\t"\
		"addpd	%%xmm2,%%xmm0							\n\t			mulpd	0x040(%%rsi),%%xmm15	\n\t"\
		"mulpd	(%%rsi),%%xmm2							\n\t			mulpd	0x030(%%rsi),%%xmm12	\n\t"\
		"movaps	%%xmm0,(%%r9 )							\n\t			mulpd	0x030(%%rsi),%%xmm11	\n\t"\
		"addpd	%%xmm2,%%xmm0							\n\t			subpd	%%xmm15,%%xmm12			\n\t"\
		"movaps	0x010(%%rsi),%%xmm2						\n\t			addpd	%%xmm8 ,%%xmm11			\n\t"\
/* lcol|rcol section uses o-offsets 1,2,3|B,C - prestore those into r10-12|14,15: */\
	"movslq	0x04(%%rdx),%%r10	\n\t	addq	%%rcx,%%r10	\n\t"/* o1 */\
	"movslq	0x08(%%rdx),%%r11	\n\t	addq	%%rcx,%%r11	\n\t"/* o2 */\
	"movslq	0x0c(%%rdx),%%r12	\n\t	addq	%%rcx,%%r12	\n\t"/* o3 */\
	"movslq	0x2c(%%rdx),%%r14	\n\t	addq	%%rcx,%%r14	\n\t"/* oB */\
	"movslq	0x30(%%rdx),%%r15	\n\t	addq	%%rcx,%%r15	\n\t"/* oC */\
		"mulpd	%%xmm2,%%xmm6							\n\t			movaps	%%xmm12,(%%r15)			\n\t"\
		"mulpd	%%xmm2,%%xmm5							\n\t			movaps	%%xmm9 ,%%xmm8			\n\t"\
		"mulpd	%%xmm2,%%xmm7							\n\t			movaps	%%xmm13,%%xmm12			\n\t"\
		"movaps %%xmm6,(%%r10)							\n\t			addpd	%%xmm10,%%xmm8			\n\t"\
		"movaps %%xmm5,(%%r11)							\n\t			addpd	%%xmm14,%%xmm12			\n\t"\
		"movaps %%xmm7,(%%r12)							\n\t			movaps	%%xmm8 ,(%%r14)			\n\t"\
		"movaps 0x080(%%rsi),%%xmm2						\n\t			movaps	%%xmm12,%%xmm15			\n\t"\
		"mulpd	%%xmm2,%%xmm5							\n\t			mulpd	0x0b0(%%rsi),%%xmm12	\n\t"\
		"mulpd	%%xmm2,%%xmm7							\n\t			mulpd	0x0e0(%%rsi),%%xmm8		\n\t"\
		"mulpd	%%xmm2,%%xmm6							\n\t			addpd	(%%r14),%%xmm12			\n\t"\
		"addpd	%%xmm0,%%xmm5							\n\t			addpd	%%xmm15,%%xmm8			\n\t"\
		"addpd	%%xmm0,%%xmm7							\n\t			mulpd	0x050(%%rsi),%%xmm12	\n\t"\
		"addpd	%%xmm0,%%xmm6							\n\t			mulpd	0x050(%%rsi),%%xmm8		\n\t"\
		"addpd	(%%r10),%%xmm5							\n\t			movaps	%%xmm10,(%%r14)			\n\t"\
		"addpd	(%%r11),%%xmm7							\n\t			movaps	%%xmm14,%%xmm15			\n\t"\
		"addpd	(%%r12),%%xmm6							\n\t			mulpd	0x0c0(%%rsi),%%xmm14	\n\t"\
		"movaps 0x0a0(%%rsi),%%xmm2						\n\t			mulpd	0x0f0(%%rsi),%%xmm10	\n\t"\
		"movaps 0x090(%%rsi),%%xmm0						\n\t			addpd	(%%r14),%%xmm14			\n\t"\
		"movaps %%xmm4,(%%r10)							\n\t			addpd	%%xmm15,%%xmm10			\n\t"\
		"movaps %%xmm1,(%%r11)							\n\t			mulpd	0x060(%%rsi),%%xmm14	\n\t"\
		"movaps %%xmm3,(%%r12)							\n\t			mulpd	0x060(%%rsi),%%xmm10	\n\t"\
		"mulpd	%%xmm2,%%xmm4							\n\t			addpd	%%xmm12,%%xmm14			\n\t"\
		"mulpd	%%xmm2,%%xmm1							\n\t			addpd	%%xmm8 ,%%xmm10			\n\t"\
		"mulpd	%%xmm2,%%xmm3							\n\t			movaps	%%xmm9 ,(%%r14)			\n\t"\
		"addpd	(%%r12),%%xmm4							\n\t			movaps	%%xmm13,%%xmm15			\n\t"\
		"subpd	(%%r10),%%xmm1							\n\t			mulpd	0x0d0(%%rsi),%%xmm13	\n\t"\
		"addpd	(%%r11),%%xmm3							\n\t			mulpd	0x100(%%rsi),%%xmm9		\n\t"\
		"mulpd	%%xmm0,%%xmm4							\n\t			addpd	(%%r14),%%xmm13			\n\t"\
		"mulpd	%%xmm0,%%xmm1							\n\t			addpd	%%xmm15,%%xmm9			\n\t"\
		"mulpd	%%xmm0,%%xmm3							\n\t			mulpd	0x070(%%rsi),%%xmm13	\n\t"\
		"movaps 0x020(%%rsi),%%xmm2						\n\t			mulpd	0x070(%%rsi),%%xmm9		\n\t"\
		"addpd	(%%r11),%%xmm4							\n\t			addpd	%%xmm12,%%xmm13			\n\t"\
		"subpd	(%%r12),%%xmm1							\n\t			addpd	%%xmm8 ,%%xmm9			\n\t"\
		"subpd	(%%r10),%%xmm3							\n\t			movaps	(%%r15),%%xmm12			\n\t"\
		"mulpd	%%xmm2,%%xmm4							\n\t			movaps	%%xmm11,%%xmm15			\n\t"\
		"mulpd	%%xmm2,%%xmm1							\n\t			movaps	%%xmm12,%%xmm8			\n\t"\
		"mulpd	%%xmm2,%%xmm3							\n\t			addpd	%%xmm14,%%xmm15			\n\t"\
		"movaps -0x010(%%rsi),%%xmm0					\n\t			addpd	%%xmm13,%%xmm8			\n\t"\
		"subpd	%%xmm4,%%xmm6							\n\t			subpd	%%xmm11,%%xmm14			\n\t"\
		"mulpd	%%xmm0,%%xmm4							\n\t			subpd	%%xmm12,%%xmm13			\n\t"\
		"subpd	%%xmm1,%%xmm7							\n\t			addpd	%%xmm10,%%xmm14			\n\t"\
		"mulpd	%%xmm0,%%xmm1							\n\t			addpd	%%xmm9 ,%%xmm13			\n\t"\
		"subpd	%%xmm3,%%xmm5							\n\t			addpd	%%xmm11,%%xmm10			\n\t"\
		"mulpd	%%xmm0,%%xmm3							\n\t			addpd	%%xmm9 ,%%xmm12			\n\t"\
		"addpd	%%xmm6,%%xmm4							\n\t"\
		"addpd	%%xmm7,%%xmm1							\n\t"\
		"addpd	%%xmm5,%%xmm3							\n\t"\
  /* yi-data in xmm8,10,12,13,14,15; xmm0,2,9,11 free: */\
	"movslq	0x18(%%rdx),%%r10	\n\t	addq	%%rcx,%%r10	\n\t"/* o6 */\
	"movslq	0x1c(%%rdx),%%r11	\n\t	addq	%%rcx,%%r11	\n\t"/* o7 */\
	"movslq	0x0c(%%rdx),%%r12	\n\t	addq	%%rcx,%%r12	\n\t"/* o3 */\
	"movslq	0x28(%%rdx),%%r13	\n\t	addq	%%rcx,%%r13	\n\t"/* oA */\
		"movaps	%%xmm4 ,%%xmm0					\n\t			movaps	%%xmm7 ,%%xmm9			\n\t"\
		"subpd	%%xmm15,%%xmm4					\n\t			subpd	%%xmm8 ,%%xmm7			\n\t"\
		"addpd	%%xmm15,%%xmm0					\n\t			addpd	%%xmm8 ,%%xmm9			\n\t"\
		"movaps	%%xmm4 ,(%%r10)					\n\t			movaps	%%xmm7 ,(%%r12)	\n\t"\
		"movaps	%%xmm0 ,(%%r11)					\n\t			movaps	%%xmm9 ,(%%r13)	\n\t"\
	"movslq	0x20(%%rdx),%%r10	\n\t	addq	%%rcx,%%r10	\n\t"/* o8 */\
	"movslq	0x14(%%rdx),%%r11	\n\t	addq	%%rcx,%%r11	\n\t"/* o5 */\
	"movslq	0x10(%%rdx),%%r12	\n\t	addq	%%rcx,%%r12	\n\t"/* o4 */\
	"movslq	0x24(%%rdx),%%r13	\n\t	addq	%%rcx,%%r13	\n\t"/* o9 */\
		"movaps	%%xmm5 ,%%xmm0					\n\t			movaps	%%xmm6 ,%%xmm9			\n\t"\
		"subpd	%%xmm14,%%xmm5					\n\t			subpd	%%xmm13,%%xmm6			\n\t"\
		"addpd	%%xmm14,%%xmm0					\n\t			addpd	%%xmm13,%%xmm9			\n\t"\
		"movaps	%%xmm5 ,(%%r10)					\n\t			movaps	%%xmm6 ,(%%r12)	\n\t"\
		"movaps	%%xmm0 ,(%%r11)					\n\t			movaps	%%xmm9 ,(%%r13)	\n\t"\
	"movslq	0x08(%%rdx),%%r10	\n\t	addq	%%rcx,%%r10	\n\t"/* o2 */\
	"movslq	0x2c(%%rdx),%%r11	\n\t	addq	%%rcx,%%r11	\n\t"/* oB */\
	"movslq	0x04(%%rdx),%%r12	\n\t	addq	%%rcx,%%r12	\n\t"/* o1 */\
	"movslq	0x30(%%rdx),%%r13	\n\t	addq	%%rcx,%%r13	\n\t"/* oC */\
		"movaps	%%xmm1 ,%%xmm0					\n\t			movaps	%%xmm3 ,%%xmm9			\n\t"\
		"subpd	%%xmm10,%%xmm1					\n\t			subpd	%%xmm12,%%xmm3			\n\t"\
		"addpd	%%xmm10,%%xmm0					\n\t			addpd	%%xmm12,%%xmm9			\n\t"\
		"movaps	%%xmm1 ,(%%r10)					\n\t			movaps	%%xmm3 ,(%%r12)	\n\t"\
		"movaps	%%xmm0 ,(%%r11)					\n\t			movaps	%%xmm9 ,(%%r13)	\n\t"\
		/****************************************************************/\
		/*                      IMAG PARTS:                             */\
		/****************************************************************/\
		"addq	$0x10,%%r8 		 		\n\t"/* I0.imag */\
				"addq	$0x10,%%r9 		\n\t"/* O0.imag */\
	/* xi-terms:                                                      yr-terms: */\
	"movslq	0x04(%%rbx),%%r15	\n\t"/* Compute i1-6 from offsets: */\
	"movslq	0x08(%%rbx),%%r14	\n\t"\
	"movslq	0x0c(%%rbx),%%r13	\n\t"\
	"movslq	0x10(%%rbx),%%r12	\n\t"\
	"movslq	0x14(%%rbx),%%r11	\n\t"\
	"movslq	0x18(%%rbx),%%r10	\n\t"\
		"addq	%%rax,%%r15	\n\t	movaps 0x10(%%r15),%%xmm1	\n\t	movaps	(%%r15),%%xmm9	\n\t"/* load A1-6 */\
		"addq	%%rax,%%r14	\n\t	movaps 0x10(%%r14),%%xmm3	\n\t	movaps	(%%r14),%%xmm10	\n\t"\
		"addq	%%rax,%%r13	\n\t	movaps 0x10(%%r13),%%xmm6	\n\t	movaps	(%%r13),%%xmm13	\n\t"\
		"addq	%%rax,%%r12	\n\t	movaps 0x10(%%r12),%%xmm4	\n\t	movaps	(%%r12),%%xmm11	\n\t"\
		"addq	%%rax,%%r11	\n\t	movaps 0x10(%%r11),%%xmm5	\n\t	movaps	(%%r11),%%xmm12	\n\t"\
		"addq	%%rax,%%r10	\n\t	movaps 0x10(%%r10),%%xmm7	\n\t	movaps	(%%r10),%%xmm14	\n\t"\
		"movaps	-0x010(%%rsi),%%xmm0	\n\t"/* 2.0 */\
	"movslq	0x30(%%rbx),%%r15	\n\t"/* Compute iC-7 from offsets: */\
	"movslq	0x2c(%%rbx),%%r14	\n\t"\
	"movslq	0x28(%%rbx),%%r13	\n\t"\
	"movslq	0x24(%%rbx),%%r12	\n\t"\
	"movslq	0x20(%%rbx),%%r11	\n\t"\
	"movslq	0x1c(%%rbx),%%r10	\n\t"\
		"addq	%%rax,%%r15	\n\t	addpd  0x10(%%r15),%%xmm1	\n\t	subpd	(%%r15),%%xmm9	\n\t"/* A1-6i += AC-7i; A1-6r -= AC-7r */\
		"addq	%%rax,%%r14	\n\t	addpd  0x10(%%r14),%%xmm3	\n\t	subpd	(%%r14),%%xmm10	\n\t"\
		"addq	%%rax,%%r13	\n\t	addpd  0x10(%%r13),%%xmm6	\n\t	subpd	(%%r13),%%xmm13	\n\t"\
		"addq	%%rax,%%r12	\n\t	addpd  0x10(%%r12),%%xmm4	\n\t	subpd	(%%r12),%%xmm11	\n\t"\
		"addq	%%rax,%%r11	\n\t	addpd  0x10(%%r11),%%xmm5	\n\t	subpd	(%%r11),%%xmm12	\n\t"\
		"addq	%%rax,%%r10	\n\t	addpd  0x10(%%r10),%%xmm7	\n\t	subpd	(%%r10),%%xmm14	\n\t"\
		"subpd	%%xmm5,%%xmm1							\n\t			movaps	%%xmm9 ,%%xmm15			\n\t"\
		"mulpd	%%xmm0,%%xmm5							\n\t			movaps	%%xmm10,%%xmm8			\n\t"\
		"subpd	%%xmm6,%%xmm3							\n\t			subpd	%%xmm11,%%xmm15			\n\t"\
		"mulpd	%%xmm0,%%xmm6							\n\t			addpd	%%xmm12,%%xmm8			\n\t"\
		"subpd	%%xmm7,%%xmm4							\n\t			addpd	%%xmm13,%%xmm15			\n\t"\
		"mulpd	%%xmm0,%%xmm7							\n\t			addpd	%%xmm14,%%xmm8			\n\t"\
		"addpd	%%xmm1,%%xmm5							\n\t			addpd	%%xmm11,%%xmm9			\n\t"\
		"addpd	%%xmm3,%%xmm6							\n\t			addpd	%%xmm11,%%xmm13			\n\t"\
		"addpd	%%xmm4,%%xmm7							\n\t			subpd	%%xmm12,%%xmm10			\n\t"\
		"movaps	%%xmm5,%%xmm2							\n\t			subpd	%%xmm12,%%xmm14			\n\t"\
		"addpd	%%xmm6,%%xmm2							\n\t			movaps	%%xmm8 ,%%xmm12			\n\t"\
		"addpd	%%xmm7,%%xmm2							\n\t			movaps	%%xmm15,%%xmm11			\n\t"\
		"movaps	(%%r8 ),%%xmm0							\n\t			mulpd	0x040(%%rsi),%%xmm8		\n\t"\
		"addpd	%%xmm2,%%xmm0							\n\t			mulpd	0x040(%%rsi),%%xmm15	\n\t"\
		"mulpd	(%%rsi),%%xmm2							\n\t			mulpd	0x030(%%rsi),%%xmm12	\n\t"\
		"movaps	%%xmm0,(%%r9 )							\n\t			mulpd	0x030(%%rsi),%%xmm11	\n\t"\
		"addpd	%%xmm2,%%xmm0							\n\t			subpd	%%xmm15,%%xmm12			\n\t"\
		"movaps	0x010(%%rsi),%%xmm2						\n\t			addpd	%%xmm8 ,%%xmm11			\n\t"\
/* rcol|lcol (note LR switcherooski) section uses o-offsets 1,2|A,B,C - prestore those into r10,11|13-15: */\
	"movslq	0x04(%%rdx),%%r10	\n\t	addq	%%rcx,%%r10	\n\t"/* o1 */\
	"movslq	0x08(%%rdx),%%r11	\n\t	addq	%%rcx,%%r11	\n\t"/* o2 */\
	"movslq	0x28(%%rdx),%%r13	\n\t	addq	%%rcx,%%r13	\n\t"/* oA */\
	"movslq	0x2c(%%rdx),%%r14	\n\t	addq	%%rcx,%%r14	\n\t"/* oB */\
	"movslq	0x30(%%rdx),%%r15	\n\t	addq	%%rcx,%%r15	\n\t"/* oC */\
		"mulpd	%%xmm2,%%xmm6							\n\t			movaps	%%xmm12,0x10(%%r10)	\n\t"\
		"mulpd	%%xmm2,%%xmm5							\n\t			movaps	%%xmm9 ,%%xmm8			\n\t"\
		"mulpd	%%xmm2,%%xmm7							\n\t			movaps	%%xmm13,%%xmm12			\n\t"\
		"movaps %%xmm6,0x10(%%r15)						\n\t			addpd	%%xmm10,%%xmm8			\n\t"\
		"movaps %%xmm5,0x10(%%r14)						\n\t			addpd	%%xmm14,%%xmm12			\n\t"\
		"movaps %%xmm7,0x10(%%r13)						\n\t			movaps	%%xmm8 ,0x10(%%r11)	\n\t"\
		"movaps 0x080(%%rsi),%%xmm2						\n\t			movaps	%%xmm12,%%xmm15			\n\t"\
		"mulpd	%%xmm2,%%xmm5							\n\t			mulpd	0x0b0(%%rsi),%%xmm12	\n\t"\
		"mulpd	%%xmm2,%%xmm7							\n\t			mulpd	0x0e0(%%rsi),%%xmm8		\n\t"\
		"mulpd	%%xmm2,%%xmm6							\n\t			addpd	0x10(%%r11),%%xmm12	\n\t"\
		"addpd	%%xmm0,%%xmm5							\n\t			addpd	%%xmm15,%%xmm8			\n\t"\
		"addpd	%%xmm0,%%xmm7							\n\t			mulpd	0x050(%%rsi),%%xmm12	\n\t"\
		"addpd	%%xmm0,%%xmm6							\n\t			mulpd	0x050(%%rsi),%%xmm8		\n\t"\
		"addpd	0x10(%%r15),%%xmm5						\n\t			movaps	%%xmm10,0x10(%%r11)	\n\t"\
		"addpd	0x10(%%r14),%%xmm7						\n\t			movaps	%%xmm14,%%xmm15			\n\t"\
		"addpd	0x10(%%r13),%%xmm6						\n\t			mulpd	0x0c0(%%rsi),%%xmm14	\n\t"\
		"movaps 0x0a0(%%rsi),%%xmm2						\n\t			mulpd	0x0f0(%%rsi),%%xmm10	\n\t"\
		"movaps 0x090(%%rsi),%%xmm0						\n\t			addpd	0x10(%%r11),%%xmm14	\n\t"\
		"movaps %%xmm4,0x10(%%r15)						\n\t			addpd	%%xmm15,%%xmm10			\n\t"\
		"movaps %%xmm1,0x10(%%r14)						\n\t			mulpd	0x060(%%rsi),%%xmm14	\n\t"\
		"movaps %%xmm3,0x10(%%r13)						\n\t			mulpd	0x060(%%rsi),%%xmm10	\n\t"\
		"mulpd	%%xmm2,%%xmm4							\n\t			addpd	%%xmm12,%%xmm14			\n\t"\
		"mulpd	%%xmm2,%%xmm1							\n\t			addpd	%%xmm8 ,%%xmm10			\n\t"\
		"mulpd	%%xmm2,%%xmm3							\n\t			movaps	%%xmm9 ,0x10(%%r11)	\n\t"\
		"addpd	0x10(%%r13),%%xmm4						\n\t			movaps	%%xmm13,%%xmm15			\n\t"\
		"subpd	0x10(%%r15),%%xmm1						\n\t			mulpd	0x0d0(%%rsi),%%xmm13	\n\t"\
		"addpd	0x10(%%r14),%%xmm3						\n\t			mulpd	0x100(%%rsi),%%xmm9		\n\t"\
		"mulpd	%%xmm0,%%xmm4							\n\t			addpd	0x10(%%r11),%%xmm13	\n\t"\
		"mulpd	%%xmm0,%%xmm1							\n\t			addpd	%%xmm15,%%xmm9			\n\t"\
		"mulpd	%%xmm0,%%xmm3							\n\t			mulpd	0x070(%%rsi),%%xmm13	\n\t"\
		"movaps 0x020(%%rsi),%%xmm2						\n\t			mulpd	0x070(%%rsi),%%xmm9		\n\t"\
		"addpd	0x10(%%r14),%%xmm4						\n\t			addpd	%%xmm12,%%xmm13			\n\t"\
		"subpd	0x10(%%r13),%%xmm1						\n\t			addpd	%%xmm8 ,%%xmm9			\n\t"\
		"subpd	0x10(%%r15),%%xmm3						\n\t			movaps	0x10(%%r10),%%xmm12	\n\t"\
		"mulpd	%%xmm2,%%xmm4							\n\t			movaps	%%xmm11,%%xmm15			\n\t"\
		"mulpd	%%xmm2,%%xmm1							\n\t			movaps	%%xmm12,%%xmm8			\n\t"\
		"mulpd	%%xmm2,%%xmm3							\n\t			addpd	%%xmm14,%%xmm15			\n\t"\
		"movaps -0x010(%%rsi),%%xmm0					\n\t			addpd	%%xmm13,%%xmm8			\n\t"\
		"subpd	%%xmm4,%%xmm6							\n\t			subpd	%%xmm11,%%xmm14			\n\t"\
		"mulpd	%%xmm0,%%xmm4							\n\t			subpd	%%xmm12,%%xmm13			\n\t"\
		"subpd	%%xmm1,%%xmm7							\n\t			addpd	%%xmm10,%%xmm14			\n\t"\
		"mulpd	%%xmm0,%%xmm1							\n\t			addpd	%%xmm9 ,%%xmm13			\n\t"\
		"subpd	%%xmm3,%%xmm5							\n\t			addpd	%%xmm11,%%xmm10			\n\t"\
		"mulpd	%%xmm0,%%xmm3							\n\t			addpd	%%xmm9 ,%%xmm12			\n\t"\
		"addpd	%%xmm6,%%xmm4							\n\t"\
		"addpd	%%xmm7,%%xmm1							\n\t"\
		"addpd	%%xmm5,%%xmm3							\n\t"\
  /* yi-data in xmm8,10,12,13,14,15; xmm0,2,9,11 free: */\
	"movslq	0x18(%%rdx),%%r10	\n\t	addq	%%rcx,%%r10	\n\t"/* o6 */\
	"movslq	0x1c(%%rdx),%%r11	\n\t	addq	%%rcx,%%r11	\n\t"/* o7 */\
	"movslq	0x0c(%%rdx),%%r12	\n\t	addq	%%rcx,%%r12	\n\t"/* o3 */\
	"movslq	0x28(%%rdx),%%r13	\n\t	addq	%%rcx,%%r13	\n\t"/* oA */\
		"movaps	%%xmm4 ,%%xmm0					\n\t			movaps	%%xmm7 ,%%xmm9		\n\t"\
		"subpd	%%xmm15,%%xmm4					\n\t			subpd	%%xmm8 ,%%xmm7		\n\t"\
		"addpd	%%xmm15,%%xmm0					\n\t			addpd	%%xmm8 ,%%xmm9		\n\t"\
		"movaps	%%xmm4 ,0x10(%%r11)				\n\t			movaps	%%xmm7 ,0x10(%%r13)	\n\t"\
		"movaps	%%xmm0 ,0x10(%%r10)				\n\t			movaps	%%xmm9 ,0x10(%%r12)	\n\t"\
	"movslq	0x20(%%rdx),%%r10	\n\t	addq	%%rcx,%%r10	\n\t"/* o8 */\
	"movslq	0x14(%%rdx),%%r11	\n\t	addq	%%rcx,%%r11	\n\t"/* o5 */\
	"movslq	0x10(%%rdx),%%r12	\n\t	addq	%%rcx,%%r12	\n\t"/* o4 */\
	"movslq	0x24(%%rdx),%%r13	\n\t	addq	%%rcx,%%r13	\n\t"/* o9 */\
		"movaps	%%xmm5 ,%%xmm0					\n\t			movaps	%%xmm6 ,%%xmm9		\n\t"\
		"subpd	%%xmm14,%%xmm5					\n\t			subpd	%%xmm13,%%xmm6		\n\t"\
		"addpd	%%xmm14,%%xmm0					\n\t			addpd	%%xmm13,%%xmm9		\n\t"\
		"movaps	%%xmm5 ,0x10(%%r11)				\n\t			movaps	%%xmm6 ,0x10(%%r13)	\n\t"\
		"movaps	%%xmm0 ,0x10(%%r10)				\n\t			movaps	%%xmm9 ,0x10(%%r12)	\n\t"\
	"movslq	0x08(%%rdx),%%r10	\n\t	addq	%%rcx,%%r10	\n\t"/* o2 */\
	"movslq	0x2c(%%rdx),%%r11	\n\t	addq	%%rcx,%%r11	\n\t"/* oB */\
	"movslq	0x04(%%rdx),%%r12	\n\t	addq	%%rcx,%%r12	\n\t"/* o1 */\
	"movslq	0x30(%%rdx),%%r13	\n\t	addq	%%rcx,%%r13	\n\t"/* oC */\
		"movaps	%%xmm1 ,%%xmm0					\n\t			movaps	%%xmm3 ,%%xmm9		\n\t"\
		"subpd	%%xmm10,%%xmm1					\n\t			subpd	%%xmm12,%%xmm3		\n\t"\
		"addpd	%%xmm10,%%xmm0					\n\t			addpd	%%xmm12,%%xmm9		\n\t"\
		"movaps	%%xmm1 ,0x10(%%r11)				\n\t			movaps	%%xmm3 ,0x10(%%r13)	\n\t"\
		"movaps	%%xmm0 ,0x10(%%r10)				\n\t			movaps	%%xmm9 ,0x10(%%r12)	\n\t"\
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

	#error Unhandled combination of preprocessor flags!

#endif	// x86 simd version ?

#endif	/* radix13_sse_macro_h_included */

