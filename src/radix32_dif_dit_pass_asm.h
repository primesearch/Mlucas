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
#ifndef radix32_dif_dit_pass_asm_h_included
#define radix32_dif_dit_pass_asm_h_included

#ifdef USE_ARM_V8_SIMD
	#if 0	// twiddles mem-map:
					cc[00 01 02 03 04 05 06 07 08 09 0A 0B 0C 0D 0E 0F 10 11 12 13 14 15 16 17 18 19 1A 1B 1C 1D 1E 1F]
		(cc0,ss0) + 0x[06,26,16,36|0e,2e,1e,3e|0a,2a,1a,3a|12,32,22,42|08,28,18,38|10,30,20,40|0c,2c,1c,3c|14,34,24,44].
	#endif
	// Vector-opcount: 95 LDP, 64 STP, 280 ADD, 176 MUL; compare to 16-DIF opcount: 50 LDP, 32 STP, 128 ADD, 96 MUL.
	#define SSE2_RADIX32_DIF_TWIDDLE(Xadd0,Xp01,Xp02,Xp03,Xp04,Xp08,Xp0C,Xp10,Xp18,Xr00)\
	{\
	__asm__ volatile (\
		"ldr	x14,%[__r00]		\n\t	add	x10, x14,#0x400			\n\t"/* x10 = &isrt2 */\
/* ~0x10, not -0x10: */"ldr w11,=0x10 \n\t	ld1r	{v29.2d},[x10],x11	\n\t"/* x10 has isrt2-ptr for LD1R, gets post-incr'ed by 0x10 to cc0 */\
		"ldr	x0,%[__add0]		\n\t	ldr	w9,%[__p02]				\n\t"\
		"ldr	w2,%[__p08]			\n\t	add	x2, x0,x2,lsl #3		\n\t"/* add0 + p8 */\
		"ldr	w4,%[__p10]			\n\t	add	x4, x0,x4,lsl #3		\n\t"/* add0 + p16 */\
		"ldr	w6,%[__p18]			\n\t	add	x6, x0,x6,lsl #3		\n\t"/* add0 + p24 */\
		"ldr	w8,%[__p04]			\n\t	add	x1, x0,x8,lsl #3		\n\t"/* add0 + p4 */\
		"ldr	w11,%[__p01]		\n\t	add	x3, x2,x8,lsl #3		\n\t"/* add0 + p12 */\
		"ldr	w12,%[__p02]		\n\t	add	x5, x4,x8,lsl #3		\n\t"/* add0 + p20 */\
		"ldr	w13,%[__p03]		\n\t	add	x7, x6,x8,lsl #3		\n\t"/* add0 + p28 */\
		/*...Block 1:
		o Inputs from add0 + p[0,4,8,12,16,20,24,28] (with resp. 4-DIFs using the [0,8,16,24],[4,12,20,28] subsets of same);
		o Outputs into r00 + [0-15];
		o Twiddles [NB: leftmost Lcol twiddle = c00 = 1.0, i.e. no-op]
			Lcol: c00,10,08,18 = cc0 + 0x[06,08,0a,0c]0
			Rcol: c04,14,0C,1C = cc0 + 0x[0e,10,12,14]0 */\
		"									ldp	q12,q13,[x1]				\n\t"\
		"ldp	q0,q1,[x0]			\n\t	ldp	q14,q15,[x5]				\n\t"\
		"ldp	q2,q3,[x4]			\n\t	ldp	q16,q17,[x10,#0x0e0]		\n\t"/* c04*/\
/* c10*/"ldp	q8,q9,[x10,#0x080]	\n\t	fmul	v18.2d,v12.2d,v16.2d	\n\t"\
		"fmul	v4.2d,v2.2d,v8.2d	\n\t	fmul	v19.2d,v13.2d,v16.2d	\n\t"\
		"fmul	v5.2d,v3.2d,v8.2d	\n\t	fmls	v18.2d,v13.2d,v17.2d	\n\t"\
		"fmls	v4.2d,v3.2d,v9.2d	\n\t	fmla	v19.2d,v12.2d,v17.2d	\n\t"\
		"fmla	v5.2d,v2.2d,v9.2d	\n\t	ldp	q16,q17,[x10,#0x100]		\n\t"/* c14*/\
		"fsub	v2.2d,v0.2d,v4.2d	\n\t	fmul	v12.2d,v14.2d,v16.2d	\n\t"\
		"fsub	v3.2d,v1.2d,v5.2d	\n\t	fmul	v13.2d,v15.2d,v16.2d	\n\t"\
		"fadd	v0.2d,v0.2d,v4.2d	\n\t	fmls	v12.2d,v15.2d,v17.2d	\n\t"\
		"fadd	v1.2d,v1.2d,v5.2d	\n\t	fmla	v13.2d,v14.2d,v17.2d	\n\t"\
		"ldp	q4,q5,[x2]			\n\t	fsub	v14.2d,v18.2d,v12.2d	\n\t"\
		"ldp	q6,q7,[x6]			\n\t	fsub	v15.2d,v19.2d,v13.2d	\n\t"\
/* c08*/"ldp	q8,q9,[x10,#0x0a0]	\n\t	fadd	v12.2d,v18.2d,v12.2d	\n\t"\
		"fmul	v10.2d,v4.2d,v8.2d	\n\t	fadd	v13.2d,v19.2d,v13.2d	\n\t"\
		"fmul	v11.2d,v5.2d,v8.2d	\n\t	ldp	q16,q17,[x3]				\n\t"\
		"fmls	v10.2d,v5.2d,v9.2d	\n\t	ldp	q18,q19,[x7]				\n\t"\
		"fmla	v11.2d,v4.2d,v9.2d	\n\t	ldp	q20,q21,[x10,#0x120]		\n\t"/* c0C*/\
/* c18*/"ldp	q8,q9,[x10,#0x0c0]	\n\t	fmul	v22.2d,v16.2d,v20.2d	\n\t"\
		"fmul	v4.2d,v6.2d,v8.2d	\n\t	fmul	v23.2d,v17.2d,v20.2d	\n\t"\
		"fmul	v5.2d,v7.2d,v8.2d	\n\t	fmls	v22.2d,v17.2d,v21.2d	\n\t"\
		"fmls	v4.2d,v7.2d,v9.2d	\n\t	fmla	v23.2d,v16.2d,v21.2d	\n\t"\
		"fmla	v5.2d,v6.2d,v9.2d	\n\t	ldp	q20,q21,[x10,#0x140]		\n\t"/* c1C*/\
		"fsub	v6.2d,v10.2d,v4.2d	\n\t	fmul	v16.2d,v18.2d,v20.2d	\n\t"\
		"fsub	v7.2d,v11.2d,v5.2d	\n\t	fmul	v17.2d,v19.2d,v20.2d	\n\t"\
		"fadd	v4.2d,v10.2d,v4.2d	\n\t	fmls	v16.2d,v19.2d,v21.2d	\n\t"\
		"fadd	v5.2d,v11.2d,v5.2d	\n\t	fmla	v17.2d,v18.2d,v21.2d	\n\t"\
										"	fsub	v18.2d,v22.2d,v16.2d	\n\t"\
										"	fsub	v19.2d,v23.2d,v17.2d	\n\t"\
										"	fadd	v16.2d,v22.2d,v16.2d	\n\t"\
										"	fadd	v17.2d,v23.2d,v17.2d	\n\t"\
		/* combine to get 2 length-4 output subtransforms... */\
		"fsub	v8.2d,v0.2d,v4.2d	\n\t	fsub	v20.2d,v12.2d,v16.2d	\n\t"\
		"fsub	v9.2d,v1.2d,v5.2d	\n\t	fsub	v21.2d,v13.2d,v17.2d	\n\t"\
		"fadd	v0.2d,v0.2d,v4.2d	\n\t	fadd	v12.2d,v12.2d,v16.2d	\n\t"\
		"fadd	v1.2d,v1.2d,v5.2d	\n\t	fadd	v13.2d,v13.2d,v17.2d	\n\t"\
		"fadd	v4.2d,v2.2d,v7.2d	\n\t	fadd	v17.2d,v14.2d,v19.2d	\n\t"\
		"fsub	v5.2d,v3.2d,v6.2d	\n\t	fsub	v16.2d,v15.2d,v18.2d	\n\t"\
		"fsub	v2.2d,v2.2d,v7.2d	\n\t	fsub	v14.2d,v14.2d,v19.2d	\n\t"\
		"fadd	v3.2d,v3.2d,v6.2d	\n\t	fadd	v15.2d,v15.2d,v18.2d	\n\t"\
	/* v6,7,18,19 free */\
		"fsub	v6.2d,v0.2d,v12.2d		\n\t"/* 1r: xmm0 = 0-12 */"fsub	v18.2d,v14.2d,v15.2d	\n\t"/* 18 =14-15 */\
		"fsub	v7.2d,v1.2d,v13.2d		\n\t"/* 1i: xmm1 = 1-13 */"fadd	v19.2d,v14.2d,v15.2d	\n\t"/* 19 =14+15 */\
		"fadd	v0.2d,v0.2d,v12.2d		\n\t"/* 0r: xmm2 = 0+12 */"fmul	v18.2d,v18.2d,v29.2d	\n\t"/* 18 *=isrt2 */\
		"fadd	v1.2d,v1.2d,v13.2d		\n\t"/* 0i: xmm3 = 1+13 */"fmul	v19.2d,v19.2d,v29.2d	\n\t"/* 19 *=isrt2 */\
		"fsub	v12.2d,v8.2d,v21.2d		\n\t"/* 2r: xmm4 = 8-21 */"fadd	v14.2d,v16.2d,v17.2d	\n\t"/* 14=16+17 */\
		"fsub	v13.2d,v9.2d,v20.2d		\n\t"/* 3i: xmm5 = 9-20 */"fsub	v15.2d,v16.2d,v17.2d	\n\t"/* 15=16-17 */\
		"fadd	v8.2d,v8.2d,v21.2d		\n\t"/* 3r: xmm7 = 8+21 */"fmul	v14.2d,v14.2d,v29.2d	\n\t"/* 14 *=isrt2 */\
		"fadd	v9.2d,v9.2d,v20.2d		\n\t"/* 2i: xmm6 = 9+20 */"fmul	v15.2d,v15.2d,v29.2d	\n\t"/* 15 *=isrt2 */\
		"fsub	v16.2d,v2.2d ,v18.2d	\n\t"/* 5r: xmm8 = 2-18 */"stp	q0 ,q1 ,[x14      ]		\n\t"\
		"fsub	v17.2d,v3.2d ,v19.2d	\n\t"/* 5i: xmm9 = 3-19 */"stp	q6 ,q7 ,[x14,#0x80]		\n\t"\
		"fadd	v2.2d ,v2.2d ,v18.2d	\n\t"/* 4r: xmm10= 2+18 */"stp	q12,q9 ,[x14,#0x40]		\n\t"\
		"fadd	v3.2d ,v3.2d ,v19.2d	\n\t"/* 4i: xmm13= 3+19 */"stp	q8 ,q13,[x14,#0xc0]		\n\t"\
		"fsub	v18.2d,v4.2d ,v14.2d	\n\t"/* 6r: xmm14= 4-14 */"stp	q2, q3 ,[x14,#0x20]		\n\t"\
		"fsub	v19.2d,v5.2d ,v15.2d	\n\t"/* 6i: xmm11= 5-15 */"stp	q16,q17,[x14,#0xa0]		\n\t"\
		"fadd	v4.2d ,v4.2d ,v14.2d	\n\t"/* 7r: xmm12= 4+14 */"stp	q18,q19,[x14,#0x60]		\n\t"\
		"fadd	v5.2d ,v5.2d ,v15.2d	\n\t"/* 7i: xmm15= 5+15 */"stp	q4 ,q5 ,[x14,#0xe0]		\n\t"\
	/***************************************/\
		/*...Block 2:
		o Inputs from add0 + p02 + p[0,4,8,12,16,20,24,28];
		o Outputs into r00 + 16 + [0-15];
		o Twiddles:
			Lcol: c02,12,0A,1A = cc0 + 0x[06,08,0a,0c]0 + 0x100
			Rcol: c06,16,0E,1E = cc0 + 0x[0e,10,12,14]0 + 0x100 */\
		"add	x0, x0,x9,lsl #3	\n\t	add	x1, x1,x9,lsl #3		\n\t"\
		"add	x2, x2,x9,lsl #3	\n\t	add	x3, x3,x9,lsl #3		\n\t"\
		"add	x4, x4,x9,lsl #3	\n\t	add	x5, x5,x9,lsl #3		\n\t"\
		"add	x6, x6,x9,lsl #3	\n\t	add	x7, x7,x9,lsl #3		\n\t"\
		"ldp	q0,q1,[x0]			\n\t	add	x10,x10,#0x100				\n\t"/* Local-array pointer += 0x100 */\
		"ldp	q2,q3,[x4]			\n\t	add	x14,x14,#0x100				\n\t"/* Twiddle-pointer += 0x100: */\
		"ldp	q8,q9,[x10,#0x060]	\n\t	ldp	q12,q13,[x1]				\n\t"\
		"fmul	v6.2d,v0.2d,v8.2d	\n\t	ldp	q14,q15,[x5]				\n\t"\
		"fmul	v7.2d,v1.2d,v8.2d	\n\t	ldp	q16,q17,[x10,#0x0e0]		\n\t"\
		"fmls	v6.2d,v1.2d,v9.2d	\n\t	fmul	v18.2d,v12.2d,v16.2d	\n\t"\
		"fmla	v7.2d,v0.2d,v9.2d	\n\t	fmul	v19.2d,v13.2d,v16.2d	\n\t"\
		"ldp	q8,q9,[x10,#0x080]	\n\t	fmls	v18.2d,v13.2d,v17.2d	\n\t"\
		"fmul	v4.2d,v2.2d,v8.2d	\n\t	fmla	v19.2d,v12.2d,v17.2d	\n\t"\
		"fmul	v5.2d,v3.2d,v8.2d	\n\t"\
		"fmls	v4.2d,v3.2d,v9.2d	\n\t"\
		"fmla	v5.2d,v2.2d,v9.2d	\n\t	ldp	q16,q17,[x10,#0x100]		\n\t"\
		"fsub	v2.2d,v6.2d,v4.2d	\n\t	fmul	v12.2d,v14.2d,v16.2d	\n\t"\
		"fsub	v3.2d,v7.2d,v5.2d	\n\t	fmul	v13.2d,v15.2d,v16.2d	\n\t"\
		"fadd	v0.2d,v6.2d,v4.2d	\n\t	fmls	v12.2d,v15.2d,v17.2d	\n\t"\
		"fadd	v1.2d,v7.2d,v5.2d	\n\t	fmla	v13.2d,v14.2d,v17.2d	\n\t"\
		"ldp	q4,q5,[x2]			\n\t	fsub	v14.2d,v18.2d,v12.2d	\n\t"\
		"ldp	q6,q7,[x6]			\n\t	fsub	v15.2d,v19.2d,v13.2d	\n\t"\
		"ldp	q8,q9,[x10,#0x0a0]	\n\t	fadd	v12.2d,v18.2d,v12.2d	\n\t"\
		"fmul	v10.2d,v4.2d,v8.2d	\n\t	fadd	v13.2d,v19.2d,v13.2d	\n\t"\
		"fmul	v11.2d,v5.2d,v8.2d	\n\t	ldp	q16,q17,[x3]				\n\t"\
		"fmls	v10.2d,v5.2d,v9.2d	\n\t	ldp	q18,q19,[x7]				\n\t"\
		"fmla	v11.2d,v4.2d,v9.2d	\n\t	ldp	q20,q21,[x10,#0x120]		\n\t"\
		"ldp	q8,q9,[x10,#0x0c0]	\n\t	fmul	v22.2d,v16.2d,v20.2d	\n\t"\
		"fmul	v4.2d,v6.2d,v8.2d	\n\t	fmul	v23.2d,v17.2d,v20.2d	\n\t"\
		"fmul	v5.2d,v7.2d,v8.2d	\n\t	fmls	v22.2d,v17.2d,v21.2d	\n\t"\
		"fmls	v4.2d,v7.2d,v9.2d	\n\t	fmla	v23.2d,v16.2d,v21.2d	\n\t"\
		"fmla	v5.2d,v6.2d,v9.2d	\n\t	ldp	q20,q21,[x10,#0x140]		\n\t"\
		"fsub	v6.2d,v10.2d,v4.2d	\n\t	fmul	v16.2d,v18.2d,v20.2d	\n\t"\
		"fsub	v7.2d,v11.2d,v5.2d	\n\t	fmul	v17.2d,v19.2d,v20.2d	\n\t"\
		"fadd	v4.2d,v10.2d,v4.2d	\n\t	fmls	v16.2d,v19.2d,v21.2d	\n\t"\
		"fadd	v5.2d,v11.2d,v5.2d	\n\t	fmla	v17.2d,v18.2d,v21.2d	\n\t"\
										"	fsub	v18.2d,v22.2d,v16.2d	\n\t"\
										"	fsub	v19.2d,v23.2d,v17.2d	\n\t"\
										"	fadd	v16.2d,v22.2d,v16.2d	\n\t"\
										"	fadd	v17.2d,v23.2d,v17.2d	\n\t"\
		/* combine to get 2 length-4 output subtransforms... */\
		"fsub	v8.2d,v0.2d,v4.2d	\n\t	fsub	v20.2d,v12.2d,v16.2d	\n\t"\
		"fsub	v9.2d,v1.2d,v5.2d	\n\t	fsub	v21.2d,v13.2d,v17.2d	\n\t"\
		"fadd	v0.2d,v0.2d,v4.2d	\n\t	fadd	v12.2d,v12.2d,v16.2d	\n\t"\
		"fadd	v1.2d,v1.2d,v5.2d	\n\t	fadd	v13.2d,v13.2d,v17.2d	\n\t"\
		"fadd	v4.2d,v2.2d,v7.2d	\n\t	fadd	v17.2d,v14.2d,v19.2d	\n\t"\
		"fsub	v5.2d,v3.2d,v6.2d	\n\t	fsub	v16.2d,v15.2d,v18.2d	\n\t"\
		"fsub	v2.2d,v2.2d,v7.2d	\n\t	fsub	v14.2d,v14.2d,v19.2d	\n\t"\
		"fadd	v3.2d,v3.2d,v6.2d	\n\t	fadd	v15.2d,v15.2d,v18.2d	\n\t"\
	/* v6,7,18,19 free */\
		"fsub	v6.2d,v0.2d,v12.2d		\n\t	fsub v18.2d,v14.2d,v15.2d	\n\t"\
		"fsub	v7.2d,v1.2d,v13.2d		\n\t	fadd v19.2d,v14.2d,v15.2d	\n\t"\
		"fadd	v0.2d,v0.2d,v12.2d		\n\t	fmul v18.2d,v18.2d,v29.2d	\n\t"\
		"fadd	v1.2d,v1.2d,v13.2d		\n\t	fmul v19.2d,v19.2d,v29.2d	\n\t"\
		"fsub	v12.2d,v8.2d,v21.2d		\n\t	fadd v14.2d,v16.2d,v17.2d	\n\t"\
		"fsub	v13.2d,v9.2d,v20.2d		\n\t	fsub v15.2d,v16.2d,v17.2d	\n\t"\
		"fadd	v8.2d,v8.2d,v21.2d		\n\t	fmul v14.2d,v14.2d,v29.2d	\n\t"\
		"fadd	v9.2d,v9.2d,v20.2d		\n\t	fmul v15.2d,v15.2d,v29.2d	\n\t"\
		"fsub	v16.2d,v2.2d ,v18.2d	\n\t	stp	q0 ,q1 ,[x14      ]		\n\t"\
		"fsub	v17.2d,v3.2d ,v19.2d	\n\t	stp	q6 ,q7 ,[x14,#0x80]		\n\t"\
		"fadd	v2.2d ,v2.2d ,v18.2d	\n\t	stp	q12,q9 ,[x14,#0x40]		\n\t"\
		"fadd	v3.2d ,v3.2d ,v19.2d	\n\t	stp	q8 ,q13,[x14,#0xc0]		\n\t"\
		"fsub	v18.2d,v4.2d ,v14.2d	\n\t	stp	q2, q3 ,[x14,#0x20]		\n\t"\
		"fsub	v19.2d,v5.2d ,v15.2d	\n\t	stp	q16,q17,[x14,#0xa0]		\n\t"\
		"fadd	v4.2d ,v4.2d ,v14.2d	\n\t	stp	q18,q19,[x14,#0x60]		\n\t"\
		"fadd	v5.2d ,v5.2d ,v15.2d	\n\t	stp	q4 ,q5 ,[x14,#0xe0]		\n\t"\
	/***************************************/\
		/*...Block 3:
		o Inputs from add0 + p01 + p[0,4,8,12,16,20,24,28];
		o Outputs into r00 + 32 + [0-15];
		o Twiddles:
			Lcol: c01,11,09,19 = cc0 + 0x[06,08,0a,0c]0 + 0x200
			Rcol: c05,15,0D,1D = cc0 + 0x[0e,10,12,14]0 + 0x200 */\
		"sub	w15,w12,w11			\n\t"/* x15 = p02-p01; due to index-padding this diff is not nec. == p01,
										i.e. to get add0+p01 from add0+p02, need to -= (p02-p01), not -= p01. */\
		"sub	x0, x0,x15,lsl #3	\n\t	sub	x1, x1,x15,lsl #3		\n\t"\
		"sub	x2, x2,x15,lsl #3	\n\t	sub	x3, x3,x15,lsl #3		\n\t"\
		"sub	x4, x4,x15,lsl #3	\n\t	sub	x5, x5,x15,lsl #3		\n\t"\
		"sub	x6, x6,x15,lsl #3	\n\t	sub	x7, x7,x15,lsl #3		\n\t"\
		"ldp	q0,q1,[x0]			\n\t	add	x10,x10,#0x100				\n\t"/* Local-array pointer += 0x100 */\
		"ldp	q2,q3,[x4]			\n\t	add	x14,x14,#0x100				\n\t"/* Twiddle-pointer += 0x100: */\
		"ldp	q8,q9,[x10,#0x060]	\n\t	ldp	q12,q13,[x1]				\n\t"\
		"fmul	v6.2d,v0.2d,v8.2d	\n\t	ldp	q14,q15,[x5]				\n\t"\
		"fmul	v7.2d,v1.2d,v8.2d	\n\t	ldp	q16,q17,[x10,#0x0e0]		\n\t"\
		"fmls	v6.2d,v1.2d,v9.2d	\n\t	fmul	v18.2d,v12.2d,v16.2d	\n\t"\
		"fmla	v7.2d,v0.2d,v9.2d	\n\t	fmul	v19.2d,v13.2d,v16.2d	\n\t"\
		"ldp	q8,q9,[x10,#0x080]	\n\t	fmls	v18.2d,v13.2d,v17.2d	\n\t"\
		"fmul	v4.2d,v2.2d,v8.2d	\n\t	fmla	v19.2d,v12.2d,v17.2d	\n\t"\
		"fmul	v5.2d,v3.2d,v8.2d	\n\t"\
		"fmls	v4.2d,v3.2d,v9.2d	\n\t"\
		"fmla	v5.2d,v2.2d,v9.2d	\n\t	ldp	q16,q17,[x10,#0x100]		\n\t"\
		"fsub	v2.2d,v6.2d,v4.2d	\n\t	fmul	v12.2d,v14.2d,v16.2d	\n\t"\
		"fsub	v3.2d,v7.2d,v5.2d	\n\t	fmul	v13.2d,v15.2d,v16.2d	\n\t"\
		"fadd	v0.2d,v6.2d,v4.2d	\n\t	fmls	v12.2d,v15.2d,v17.2d	\n\t"\
		"fadd	v1.2d,v7.2d,v5.2d	\n\t	fmla	v13.2d,v14.2d,v17.2d	\n\t"\
		"ldp	q4,q5,[x2]			\n\t	fsub	v14.2d,v18.2d,v12.2d	\n\t"\
		"ldp	q6,q7,[x6]			\n\t	fsub	v15.2d,v19.2d,v13.2d	\n\t"\
		"ldp	q8,q9,[x10,#0x0a0]	\n\t	fadd	v12.2d,v18.2d,v12.2d	\n\t"\
		"fmul	v10.2d,v4.2d,v8.2d	\n\t	fadd	v13.2d,v19.2d,v13.2d	\n\t"\
		"fmul	v11.2d,v5.2d,v8.2d	\n\t	ldp	q16,q17,[x3]				\n\t"\
		"fmls	v10.2d,v5.2d,v9.2d	\n\t	ldp	q18,q19,[x7]				\n\t"\
		"fmla	v11.2d,v4.2d,v9.2d	\n\t	ldp	q20,q21,[x10,#0x120]		\n\t"\
		"ldp	q8,q9,[x10,#0x0c0]	\n\t	fmul	v22.2d,v16.2d,v20.2d	\n\t"\
		"fmul	v4.2d,v6.2d,v8.2d	\n\t	fmul	v23.2d,v17.2d,v20.2d	\n\t"\
		"fmul	v5.2d,v7.2d,v8.2d	\n\t	fmls	v22.2d,v17.2d,v21.2d	\n\t"\
		"fmls	v4.2d,v7.2d,v9.2d	\n\t	fmla	v23.2d,v16.2d,v21.2d	\n\t"\
		"fmla	v5.2d,v6.2d,v9.2d	\n\t	ldp	q20,q21,[x10,#0x140]		\n\t"\
		"fsub	v6.2d,v10.2d,v4.2d	\n\t	fmul	v16.2d,v18.2d,v20.2d	\n\t"\
		"fsub	v7.2d,v11.2d,v5.2d	\n\t	fmul	v17.2d,v19.2d,v20.2d	\n\t"\
		"fadd	v4.2d,v10.2d,v4.2d	\n\t	fmls	v16.2d,v19.2d,v21.2d	\n\t"\
		"fadd	v5.2d,v11.2d,v5.2d	\n\t	fmla	v17.2d,v18.2d,v21.2d	\n\t"\
										"	fsub	v18.2d,v22.2d,v16.2d	\n\t"\
										"	fsub	v19.2d,v23.2d,v17.2d	\n\t"\
										"	fadd	v16.2d,v22.2d,v16.2d	\n\t"\
										"	fadd	v17.2d,v23.2d,v17.2d	\n\t"\
		/* combine to get 2 length-4 output subtransforms... */\
		"fsub	v8.2d,v0.2d,v4.2d	\n\t	fsub	v20.2d,v12.2d,v16.2d	\n\t"\
		"fsub	v9.2d,v1.2d,v5.2d	\n\t	fsub	v21.2d,v13.2d,v17.2d	\n\t"\
		"fadd	v0.2d,v0.2d,v4.2d	\n\t	fadd	v12.2d,v12.2d,v16.2d	\n\t"\
		"fadd	v1.2d,v1.2d,v5.2d	\n\t	fadd	v13.2d,v13.2d,v17.2d	\n\t"\
		"fadd	v4.2d,v2.2d,v7.2d	\n\t	fadd	v17.2d,v14.2d,v19.2d	\n\t"\
		"fsub	v5.2d,v3.2d,v6.2d	\n\t	fsub	v16.2d,v15.2d,v18.2d	\n\t"\
		"fsub	v2.2d,v2.2d,v7.2d	\n\t	fsub	v14.2d,v14.2d,v19.2d	\n\t"\
		"fadd	v3.2d,v3.2d,v6.2d	\n\t	fadd	v15.2d,v15.2d,v18.2d	\n\t"\
	/* v6,7,18,19 free */\
		"fsub	v6.2d,v0.2d,v12.2d		\n\t	fsub v18.2d,v14.2d,v15.2d	\n\t"\
		"fsub	v7.2d,v1.2d,v13.2d		\n\t	fadd v19.2d,v14.2d,v15.2d	\n\t"\
		"fadd	v0.2d,v0.2d,v12.2d		\n\t	fmul v18.2d,v18.2d,v29.2d	\n\t"\
		"fadd	v1.2d,v1.2d,v13.2d		\n\t	fmul v19.2d,v19.2d,v29.2d	\n\t"\
		"fsub	v12.2d,v8.2d,v21.2d		\n\t	fadd v14.2d,v16.2d,v17.2d	\n\t"\
		"fsub	v13.2d,v9.2d,v20.2d		\n\t	fsub v15.2d,v16.2d,v17.2d	\n\t"\
		"fadd	v8.2d,v8.2d,v21.2d		\n\t	fmul v14.2d,v14.2d,v29.2d	\n\t"\
		"fadd	v9.2d,v9.2d,v20.2d		\n\t	fmul v15.2d,v15.2d,v29.2d	\n\t"\
		"fsub	v16.2d,v2.2d ,v18.2d	\n\t	stp	q0 ,q1 ,[x14      ]		\n\t"\
		"fsub	v17.2d,v3.2d ,v19.2d	\n\t	stp	q6 ,q7 ,[x14,#0x80]		\n\t"\
		"fadd	v2.2d ,v2.2d ,v18.2d	\n\t	stp	q12,q9 ,[x14,#0x40]		\n\t"\
		"fadd	v3.2d ,v3.2d ,v19.2d	\n\t	stp	q8 ,q13,[x14,#0xc0]		\n\t"\
		"fsub	v18.2d,v4.2d ,v14.2d	\n\t	stp	q2, q3 ,[x14,#0x20]		\n\t"\
		"fsub	v19.2d,v5.2d ,v15.2d	\n\t	stp	q16,q17,[x14,#0xa0]		\n\t"\
		"fadd	v4.2d ,v4.2d ,v14.2d	\n\t	stp	q18,q19,[x14,#0x60]		\n\t"\
		"fadd	v5.2d ,v5.2d ,v15.2d	\n\t	stp	q4 ,q5 ,[x14,#0xe0]		\n\t"\
	/***************************************/\
		/*...Block 4:
		o Inputs from add0 + p03 + p[0,4,8,12,16,20,24,28];
		o Outputs into r00 + 48 + [0-15];
		o Twiddles:
			Lcol: c03,13,0B,1B = cc0 + 0x[06,08,0a,0c]0 + 0x300
			Rcol: c07,17,0F,1F = cc0 + 0x[0e,10,12,14]0 + 0x300 */\
		/* Index-padding scheme guarantees that p01+p02 == p03, so to get add0+p03 from add0+p01, simply += p02: */\
		"add	x0, x0,x9,lsl #3	\n\t	add	x1, x1,x9,lsl #3		\n\t"\
		"add	x2, x2,x9,lsl #3	\n\t	add	x3, x3,x9,lsl #3		\n\t"\
		"add	x4, x4,x9,lsl #3	\n\t	add	x5, x5,x9,lsl #3		\n\t"\
		"add	x6, x6,x9,lsl #3	\n\t	add	x7, x7,x9,lsl #3		\n\t"\
		"ldp	q0,q1,[x0]			\n\t	add	x10,x10,#0x100				\n\t"/* Local-array pointer += 0x100 */\
		"ldp	q2,q3,[x4]			\n\t	add	x14,x14,#0x100				\n\t"/* Twiddle-pointer += 0x100: */\
		"ldp	q8,q9,[x10,#0x060]	\n\t	ldp	q12,q13,[x1]				\n\t"\
		"fmul	v6.2d,v0.2d,v8.2d	\n\t	ldp	q14,q15,[x5]				\n\t"\
		"fmul	v7.2d,v1.2d,v8.2d	\n\t	ldp	q16,q17,[x10,#0x0e0]		\n\t"\
		"fmls	v6.2d,v1.2d,v9.2d	\n\t	fmul	v18.2d,v12.2d,v16.2d	\n\t"\
		"fmla	v7.2d,v0.2d,v9.2d	\n\t	fmul	v19.2d,v13.2d,v16.2d	\n\t"\
		"ldp	q8,q9,[x10,#0x080]	\n\t	fmls	v18.2d,v13.2d,v17.2d	\n\t"\
		"fmul	v4.2d,v2.2d,v8.2d	\n\t	fmla	v19.2d,v12.2d,v17.2d	\n\t"\
		"fmul	v5.2d,v3.2d,v8.2d	\n\t"\
		"fmls	v4.2d,v3.2d,v9.2d	\n\t"\
		"fmla	v5.2d,v2.2d,v9.2d	\n\t	ldp	q16,q17,[x10,#0x100]		\n\t"\
		"fsub	v2.2d,v6.2d,v4.2d	\n\t	fmul	v12.2d,v14.2d,v16.2d	\n\t"\
		"fsub	v3.2d,v7.2d,v5.2d	\n\t	fmul	v13.2d,v15.2d,v16.2d	\n\t"\
		"fadd	v0.2d,v6.2d,v4.2d	\n\t	fmls	v12.2d,v15.2d,v17.2d	\n\t"\
		"fadd	v1.2d,v7.2d,v5.2d	\n\t	fmla	v13.2d,v14.2d,v17.2d	\n\t"\
		"ldp	q4,q5,[x2]			\n\t	fsub	v14.2d,v18.2d,v12.2d	\n\t"\
		"ldp	q6,q7,[x6]			\n\t	fsub	v15.2d,v19.2d,v13.2d	\n\t"\
		"ldp	q8,q9,[x10,#0x0a0]	\n\t	fadd	v12.2d,v18.2d,v12.2d	\n\t"\
		"fmul	v10.2d,v4.2d,v8.2d	\n\t	fadd	v13.2d,v19.2d,v13.2d	\n\t"\
		"fmul	v11.2d,v5.2d,v8.2d	\n\t	ldp	q16,q17,[x3]				\n\t"\
		"fmls	v10.2d,v5.2d,v9.2d	\n\t	ldp	q18,q19,[x7]				\n\t"\
		"fmla	v11.2d,v4.2d,v9.2d	\n\t	ldp	q20,q21,[x10,#0x120]		\n\t"\
		"ldp	q8,q9,[x10,#0x0c0]	\n\t	fmul	v22.2d,v16.2d,v20.2d	\n\t"\
		"fmul	v4.2d,v6.2d,v8.2d	\n\t	fmul	v23.2d,v17.2d,v20.2d	\n\t"\
		"fmul	v5.2d,v7.2d,v8.2d	\n\t	fmls	v22.2d,v17.2d,v21.2d	\n\t"\
		"fmls	v4.2d,v7.2d,v9.2d	\n\t	fmla	v23.2d,v16.2d,v21.2d	\n\t"\
		"fmla	v5.2d,v6.2d,v9.2d	\n\t	ldp	q20,q21,[x10,#0x140]		\n\t"\
		"fsub	v6.2d,v10.2d,v4.2d	\n\t	fmul	v16.2d,v18.2d,v20.2d	\n\t"\
		"fsub	v7.2d,v11.2d,v5.2d	\n\t	fmul	v17.2d,v19.2d,v20.2d	\n\t"\
		"fadd	v4.2d,v10.2d,v4.2d	\n\t	fmls	v16.2d,v19.2d,v21.2d	\n\t"\
		"fadd	v5.2d,v11.2d,v5.2d	\n\t	fmla	v17.2d,v18.2d,v21.2d	\n\t"\
										"	fsub	v18.2d,v22.2d,v16.2d	\n\t"\
										"	fsub	v19.2d,v23.2d,v17.2d	\n\t"\
										"	fadd	v16.2d,v22.2d,v16.2d	\n\t"\
										"	fadd	v17.2d,v23.2d,v17.2d	\n\t"\
		/* combine to get 2 length-4 output subtransforms... */\
		"fsub	v8.2d,v0.2d,v4.2d	\n\t	fsub	v20.2d,v12.2d,v16.2d	\n\t"\
		"fsub	v9.2d,v1.2d,v5.2d	\n\t	fsub	v21.2d,v13.2d,v17.2d	\n\t"\
		"fadd	v0.2d,v0.2d,v4.2d	\n\t	fadd	v12.2d,v12.2d,v16.2d	\n\t"\
		"fadd	v1.2d,v1.2d,v5.2d	\n\t	fadd	v13.2d,v13.2d,v17.2d	\n\t"\
		"fadd	v4.2d,v2.2d,v7.2d	\n\t	fadd	v17.2d,v14.2d,v19.2d	\n\t"\
		"fsub	v5.2d,v3.2d,v6.2d	\n\t	fsub	v16.2d,v15.2d,v18.2d	\n\t"\
		"fsub	v2.2d,v2.2d,v7.2d	\n\t	fsub	v14.2d,v14.2d,v19.2d	\n\t"\
		"fadd	v3.2d,v3.2d,v6.2d	\n\t	fadd	v15.2d,v15.2d,v18.2d	\n\t"\
	/* v6,7,18,19 free */\
		"fsub	v6.2d,v0.2d,v12.2d		\n\t	fsub v18.2d,v14.2d,v15.2d	\n\t"\
		"fsub	v7.2d,v1.2d,v13.2d		\n\t	fadd v19.2d,v14.2d,v15.2d	\n\t"\
		"fadd	v0.2d,v0.2d,v12.2d		\n\t	fmul v18.2d,v18.2d,v29.2d	\n\t"\
		"fadd	v1.2d,v1.2d,v13.2d		\n\t	fmul v19.2d,v19.2d,v29.2d	\n\t"\
		"fsub	v12.2d,v8.2d,v21.2d		\n\t	fadd v14.2d,v16.2d,v17.2d	\n\t"\
		"fsub	v13.2d,v9.2d,v20.2d		\n\t	fsub v15.2d,v16.2d,v17.2d	\n\t"\
		"fadd	v8.2d,v8.2d,v21.2d		\n\t	fmul v14.2d,v14.2d,v29.2d	\n\t"\
		"fadd	v9.2d,v9.2d,v20.2d		\n\t	fmul v15.2d,v15.2d,v29.2d	\n\t"\
		"fsub	v16.2d,v2.2d ,v18.2d	\n\t	stp	q0 ,q1 ,[x14      ]		\n\t"\
		"fsub	v17.2d,v3.2d ,v19.2d	\n\t	stp	q6 ,q7 ,[x14,#0x80]		\n\t"\
		"fadd	v2.2d ,v2.2d ,v18.2d	\n\t	stp	q12,q9 ,[x14,#0x40]		\n\t"\
		"fadd	v3.2d ,v3.2d ,v19.2d	\n\t	stp	q8 ,q13,[x14,#0xc0]		\n\t"\
		"fsub	v18.2d,v4.2d ,v14.2d	\n\t	stp	q2, q3 ,[x14,#0x20]		\n\t"\
		"fsub	v19.2d,v5.2d ,v15.2d	\n\t	stp	q16,q17,[x14,#0xa0]		\n\t"\
		"fadd	v4.2d ,v4.2d ,v14.2d	\n\t	stp	q18,q19,[x14,#0x60]		\n\t"\
		"fadd	v5.2d ,v5.2d ,v15.2d	\n\t	stp	q4 ,q5 ,[x14,#0xe0]		\n\t"\
		"sub	x10,x10,#0x300			\n\t"/* Twiddle-pointer -= 0x300 to get it to point to cc0: */\
	/**********************************************************************************/\
	/*...and now do eight radix-4 transforms, including the internal twiddle factors: */\
	/**********************************************************************************/\
	/*...Block 1: t00,t10,t20,t30	*/		/*...Block 5: t08,t18,t28,t38	*/\
		"ldr	x14,%[__r00]			\n\t	ldr	w9,%[__p08]	\n\t"/* p04 still in w8 */\
		"ldr	x0,%[__add0]			\n\t	ldp	q8 ,q9 ,[x14,#0x080]	\n\t"/* lcol: add0 + p0 */\
		"add	x1,x0,x11,lsl #3		\n\t	ldp	q10,q11,[x14,#0x180]	\n\t"/* lcol: add0 + p1 */\
		"add	x2,x0,x12,lsl #3		\n\t	ldp	q12,q13,[x14,#0x280]	\n\t"/* lcol: add0 + p2 */\
		"add	x3,x0,x13,lsl #3		\n\t	ldp	q14,q15,[x14,#0x380]	\n\t"/* lcol: add0 + p3 */\
		"ldp	q0,q1,[x14       ]		\n\t	add	x4,x0,x8,lsl #3			\n\t"/* rcol: add0 + p4 */\
		"ldp	q2,q3,[x14,#0x100]		\n\t	add	x5,x1,x8,lsl #3			\n\t"/* rcol: add0 + p5 */\
		"ldp	q4,q5,[x14,#0x200]		\n\t	add	x6,x2,x8,lsl #3			\n\t"/* rcol: add0 + p6 */\
		"ldp	q6,q7,[x14,#0x300]		\n\t	add	x7,x3,x8,lsl #3			\n\t"/* rcol: add0 + p7 */\
/* 2 */	"fadd	v16.2d,v0.2d,v2.2d		\n\t	fmul	v12.2d,v12.2d,v29.2d	\n\t"\
/* 3 */	"fadd	v17.2d,v1.2d,v3.2d		\n\t	fmul	v13.2d,v13.2d,v29.2d	\n\t"\
/* 6 */	"fadd	v18.2d,v4.2d,v6.2d		\n\t	fmul	v14.2d,v14.2d,v29.2d	\n\t"\
/* 7 */	"fadd	v19.2d,v5.2d,v7.2d		\n\t	fmul	v15.2d,v15.2d,v29.2d	\n\t"\
		"fsub	v0.2d ,v0.2d,v2.2d		\n\t	fadd	v20.2d,v8.2d ,v11.2d	\n\t"/* b=8+b */\
		"fsub	v1.2d ,v1.2d,v3.2d		\n\t	fadd	v21.2d,v12.2d,v13.2d	\n\t"/* d=c+d */\
		"fsub	v4.2d ,v4.2d,v6.2d		\n\t	fadd	v22.2d,v9.2d ,v10.2d	\n\t"/* a=9+a */\
		"fsub	v5.2d ,v5.2d,v7.2d		\n\t	fadd	v23.2d,v15.2d,v14.2d	\n\t"/* e=f+e */\
		"fsub	v2.2d,v16.2d,v18.2d		\n\t	fsub	v8.2d ,v8.2d ,v11.2d	\n\t"/* 2=2-6 ; 8=8-b */\
		"fsub	v3.2d,v17.2d,v19.2d		\n\t	fsub	v12.2d,v12.2d,v13.2d	\n\t"/* 3=3-7 ; c=c-d */\
		"fadd	v6.2d,v16.2d,v18.2d		\n\t	fsub	v9.2d ,v9.2d ,v10.2d	\n\t"/* 6=2+6 ; 9=9-a */\
		"fadd	v7.2d,v17.2d,v19.2d		\n\t	fsub	v15.2d,v15.2d,v14.2d	\n\t"/* 7=3+7 ; f=f-e */\
/* 0 */	"fsub	v16.2d,v0.2d,v5.2d		\n\t	fadd	v14.2d,v12.2d,v23.2d	\n\t"/* e=c+e */\
/* 1 */	"fsub	v17.2d,v1.2d,v4.2d		\n\t	fsub	v12.2d,v12.2d,v23.2d	\n\t"/* c=c-e */\
		"fadd	v5.2d ,v0.2d,v5.2d		\n\t	fsub	v13.2d,v21.2d,v15.2d	\n\t"/* d=d-f */\
		"fadd	v4.2d ,v1.2d,v4.2d		\n\t	fadd	v15.2d,v21.2d,v15.2d	\n\t"/* f=d+f */\
		/* Lcol Output order is 6,7,2,3,0,4,5,1: */\
		"stp	q6 ,q7 ,[x0]			\n\t	fsub	v21.2d,v8.2d ,v12.2d	\n\t"/* 8=8-c */\
		"stp	q2, q3 ,[x1]			\n\t	fsub	v10.2d,v22.2d,v13.2d	\n\t"/* a=a-d */\
		"stp	q16,q4 ,[x2]			\n\t	fadd	v12.2d,v8.2d ,v12.2d	\n\t"/* c=8+c */\
		"stp	q5 ,q17,[x3]			\n\t	fadd	v13.2d,v22.2d,v13.2d	\n\t"/* d=a+d */\
		/* Rcol Output order is c,d,8,a,b,e,f,9 - fiddle actual store order to follow compute-order mandated by register-data dependencies: */\
		"add	x0,x0,x9,lsl #3			\n\t	fsub	v11.2d,v20.2d,v15.2d	\n\t"/* b=b-f */"	stp	q12,q13,[x4]	\n\t"\
		"add	x1,x1,x9,lsl #3			\n\t	fsub	v22.2d,v9.2d ,v14.2d	\n\t"/* 9=9-e */"	stp	q21,q10,[x5]	\n\t"\
		"add	x2,x2,x9,lsl #3			\n\t	fadd	v15.2d,v20.2d,v15.2d	\n\t"/* f=b+f */"	stp	q15,q22,[x7]	\n\t"\
		"add	x3,x3,x9,lsl #3			\n\t	fadd	v14.2d,v9.2d ,v14.2d	\n\t"/* e=9+e */"	stp	q11,q14,[x6]	\n\t"\
	/*...Block 3: t04,t14,t24,t34; outputs into add0 + p8-11: */\
		"add	x14,x14,#0x40	\n\t"/* r04 */"	ldp	q16,q17,[x10]			\n\t"/* cc0 */\
											/*...Block 7: t0C,t1C,t2C,t3C; outputs into add0 + p12-15: */\
		"ldp	q0 ,q1 ,[x14       ]	\n\t	add	x4,x4,x9,lsl #3			\n\t"/* rcol: add0 + p12 */\
		"ldp	q2 ,q3 ,[x14,#0x100]	\n\t	add	x5,x5,x9,lsl #3			\n\t"/* rcol: add0 + p13 */\
		"ldp	q4 ,q5 ,[x14,#0x200]	\n\t	add	x6,x6,x9,lsl #3			\n\t"/* rcol: add0 + p14 */\
		"ldp	q20,q21,[x14,#0x300]	\n\t	add	x7,x7,x9,lsl #3			\n\t"/* rcol: add0 + p15 */\
		"fmul	v18.2d,v4.2d,v16.2d		\n\t	ldp	q8 ,q9 ,[x14,#0x080]	\n\t"\
		"fmul	v19.2d,v5.2d,v16.2d		\n\t	ldp	q10,q11,[x14,#0x180]	\n\t"\
		"fmls	v18.2d,v5.2d,v17.2d		\n\t	ldp	q12,q13,[x14,#0x280]	\n\t"\
		"fmla	v19.2d,v4.2d,v17.2d		\n\t	ldp	q24,q25,[x14,#0x380]	\n\t"\
		"fmul	v6.2d,v20.2d,v17.2d		\n\t	fmul	v22.2d,v12.2d,v17.2d	\n\t"/* 12s */\
		"fmul	v7.2d,v21.2d,v17.2d		\n\t	fmul	v23.2d,v13.2d,v17.2d	\n\t"/* 13s */\
		"fmls	v6.2d,v21.2d,v16.2d		\n\t	fmls	v22.2d,v13.2d,v16.2d	\n\t"/* 12s-13c */\
		"fmla	v7.2d,v20.2d,v16.2d		\n\t	fmla	v23.2d,v12.2d,v16.2d	\n\t"/* 13s+12c */\
		"fsub	v4.2d,v18.2d,v6.2d		\n\t	fmul	v14.2d,v24.2d,v16.2d	\n\t"/* 14c */\
		"fsub	v5.2d,v19.2d,v7.2d		\n\t	fmul	v15.2d,v25.2d,v16.2d	\n\t"/* 15c */\
		"fadd	v6.2d,v18.2d,v6.2d		\n\t	fmls	v14.2d,v25.2d,v17.2d	\n\t"/* 14=14c-15s */\
		"fadd	v7.2d,v19.2d,v7.2d		\n\t	fmla	v15.2d,v24.2d,v17.2d	\n\t"/* 15=15c+14s */\
		"fsub	v20.2d,v2.2d,v3.2d		\n\t	fsub	v12.2d,v22.2d,v14.2d	\n\t"/* 12=12-14 */\
		"fadd	v21.2d,v2.2d,v3.2d		\n\t	fsub	v13.2d,v23.2d,v15.2d	\n\t"/* 13=13-15 */\
		"fmul	v20.2d,v20.2d,v29.2d	\n\t	fadd	v14.2d,v22.2d,v14.2d	\n\t"/* 14=12+14 */\
		"fmul	v21.2d,v21.2d,v29.2d	\n\t	fadd	v15.2d,v23.2d,v15.2d	\n\t"/* 15=13+15 */\
		"fadd	v2.2d,v0.2d,v20.2d		\n\t	fadd	v22.2d,v10.2d,v11.2d	\n\t"/* 10=10+11 */\
		"fadd	v3.2d,v1.2d,v21.2d		\n\t	fsub	v23.2d,v11.2d,v10.2d	\n\t"/* 11=11-10 */\
		"fsub	v0.2d,v0.2d,v20.2d		\n\t	fmul	v22.2d,v22.2d,v29.2d	\n\t"/* 10 *= isrt2 */\
		"fsub	v1.2d,v1.2d,v21.2d		\n\t	fmul	v23.2d,v23.2d,v29.2d	\n\t"/* 11 *= isrt2 */\
/* 6=2+6*/"fadd	v18.2d,v2.2d,v6.2d		\n\t	fadd	v10.2d,v8.2d,v22.2d		\n\t"/* 10=8+10 */\
/* 7=3+7*/"fadd	v19.2d,v3.2d,v7.2d		\n\t	fadd	v11.2d,v9.2d,v23.2d		\n\t"/* 11=9+11 */\
/* 2=2-6*/"fsub	v2.2d ,v2.2d,v6.2d		\n\t	fsub	v8.2d ,v8.2d,v22.2d		\n\t"/*  8=8-10 */\
/* 3=3-7*/"fsub	v3.2d ,v3.2d,v7.2d		\n\t	fsub	v9.2d ,v9.2d,v23.2d		\n\t"/*  9=9-11 */\
/* 0=0-5*/"fsub	v20.2d,v0.2d,v5.2d		\n\t	fadd	v25.2d,v10.2d,v15.2d	\n\t"/* 15=10+15 */\
/* 4=1+4*/"fadd	v21.2d,v1.2d,v4.2d		\n\t	fadd	v24.2d,v11.2d,v14.2d	\n\t"/* 14=11+14 */\
/* 5=0+5*/"fadd	v5.2d ,v0.2d,v5.2d		\n\t	fsub	v10.2d,v10.2d,v15.2d	\n\t"/* 10=10-15 */\
/* 1=1-4*/"fsub	v1.2d ,v1.2d,v4.2d		\n\t	fsub	v11.2d,v11.2d,v14.2d	\n\t"/* 11=11-14 */\
		/* Lcol Output order: 6,7,2,3,0,4,5,1: */\
		"stp	q18,q19,[x0]			\n\t	fsub	v22.2d,v8.2d ,v12.2d	\n\t"/*  8= 8-12 */\
		"stp	q2 ,q3 ,[x1]			\n\t	fsub	v23.2d,v9.2d ,v13.2d	\n\t"/*  9= 9-13 */\
		"stp	q20,q21,[x2]			\n\t	fadd	v12.2d,v8.2d ,v12.2d	\n\t"/* 12= 8+12 */\
		"stp	q5 ,q1 ,[x3]			\n\t	fadd	v13.2d,v9.2d ,v13.2d	\n\t"/* 13= 9+13 */\
		"ldr	x0,%[__add0]			\n\t	ldr	w9,%[__p10]	\n\t"/* Due to maths of padded-indexing, can't just incr. by another += p08 */\
		/* Rcol Output order: 12,13,8,9,10,14,15,11; Fiddle order of Rcol stores to reflect order in which outputs computed: */\
		"add	x0,x0,x9 ,lsl #3		\n\t	stp	q10,q24,[x6]		\n\t"/* lcol: add0 + p16 */\
		"add	x1,x0,x11,lsl #3		\n\t	stp	q25,q11,[x7]		\n\t"/* lcol: add0 + p17 */\
		"add	x2,x0,x12,lsl #3		\n\t	stp	q22,q23,[x5]		\n\t"/* lcol: add0 + p18 */\
		"add	x3,x0,x13,lsl #3		\n\t	stp	q12,q13,[x4]		\n\t"/* lcol: add0 + p19 */\
	/*...Block 2: t02,t12,t22,t32: */		/*...Block 6: t0A,t1A,t2A,t3A; outputs into add0 + p20-23: */\
		"sub	x14,x14,#0x20			\n\t"/* r02 */\
		"ldp	q16,q17,[x10,#0x20]		\n\t	ldp	q14,q15,[x10,#0x40]		\n\t"/* cc1,cc3 */\
												/* p04 still in w8: */\
		"ldp	q0 ,q1 ,[x14       ]	\n\t	add	x4,x0,x8,lsl #3			\n\t"/* rcol: add0 + p20 */\
		"ldp	q2 ,q3 ,[x14,#0x100]	\n\t	add	x5,x1,x8,lsl #3			\n\t"/* rcol: add0 + p21 */\
		"ldp	q4 ,q5 ,[x14,#0x200]	\n\t	add	x6,x2,x8,lsl #3			\n\t"/* rcol: add0 + p22 */\
		"ldp	q20,q21,[x14,#0x300]	\n\t	add	x7,x3,x8,lsl #3			\n\t"/* rcol: add0 + p23 */\
		"fmul	v18.2d,v4.2d,v16.2d		\n\t	ldp	q8 ,q9 ,[x14,#0x080]	\n\t"\
		"fmul	v19.2d,v5.2d,v16.2d		\n\t	ldp	q10,q11,[x14,#0x180]	\n\t"\
		"fmls	v18.2d,v5.2d,v17.2d		\n\t	ldp	q12,q13,[x14,#0x280]	\n\t"\
		"fmla	v19.2d,v4.2d,v17.2d		\n\t	ldp	q24,q25,[x14,#0x380]	\n\t"\
		/* Next set of tewiddles swap [c,s] pairs between the 2 columns: */\
		"fmul	v6.2d,v20.2d,v14.2d		\n\t	fmul	v22.2d,v12.2d,v15.2d	\n\t"\
		"fmul	v7.2d,v21.2d,v14.2d		\n\t	fmul	v23.2d,v13.2d,v15.2d	\n\t"\
		"fmls	v6.2d,v21.2d,v15.2d		\n\t	fmls	v22.2d,v13.2d,v14.2d	\n\t"\
		"fmla	v7.2d,v20.2d,v15.2d		\n\t	fmla	v23.2d,v12.2d,v14.2d	\n\t"\
		"fsub	v4.2d,v18.2d,v6.2d		\n\t	fmul	v14.2d,v24.2d,v16.2d	\n\t"\
		"fsub	v5.2d,v19.2d,v7.2d		\n\t	fmul	v15.2d,v25.2d,v16.2d	\n\t"\
		"fadd	v6.2d,v18.2d,v6.2d		\n\t	fmla	v14.2d,v25.2d,v17.2d	\n\t"\
		"fadd	v7.2d,v19.2d,v7.2d		\n\t	fmls	v15.2d,v24.2d,v17.2d	\n\t"\
		"ldp	q16,q17,[x10     ]		\n\t"/* cc0 */\
		"fmul	v18.2d,v2.2d,v16.2d		\n\t	fsub	v12.2d,v22.2d,v14.2d	\n\t"\
		"fmul	v19.2d,v3.2d,v16.2d		\n\t	fsub	v13.2d,v23.2d,v15.2d	\n\t"\
		"fmls	v18.2d,v3.2d,v17.2d		\n\t	fadd	v14.2d,v22.2d,v14.2d	\n\t"\
		"fmla	v19.2d,v2.2d,v17.2d		\n\t	fadd	v15.2d,v23.2d,v15.2d	\n\t"\
		/* Post-twiddle sub-DFT sequence same as pvs blockpair, just note 2,3 in v18,v19: */\
		"fadd	v2.2d,v0.2d,v18.2d		\n\t	fmul	v22.2d,v10.2d,v17.2d	\n\t"\
		"fadd	v3.2d,v1.2d,v19.2d		\n\t	fmul	v23.2d,v11.2d,v17.2d	\n\t"\
		"fsub	v0.2d,v0.2d,v18.2d		\n\t	fmla	v22.2d,v11.2d,v16.2d	\n\t"\
		"fsub	v1.2d,v1.2d,v19.2d		\n\t	fmls	v23.2d,v10.2d,v16.2d	\n\t"\
		"fadd	v18.2d,v2.2d,v6.2d		\n\t	fadd	v10.2d,v8.2d,v22.2d		\n\t"\
		"fadd	v19.2d,v3.2d,v7.2d		\n\t	fadd	v11.2d,v9.2d,v23.2d		\n\t"\
		"fsub	v2.2d ,v2.2d,v6.2d		\n\t	fsub	v8.2d ,v8.2d,v22.2d		\n\t"\
		"fsub	v3.2d ,v3.2d,v7.2d		\n\t	fsub	v9.2d ,v9.2d,v23.2d		\n\t"\
		"fsub	v20.2d,v0.2d,v5.2d		\n\t	fadd	v25.2d,v10.2d,v15.2d	\n\t"\
		"fadd	v21.2d,v1.2d,v4.2d		\n\t	fadd	v24.2d,v11.2d,v14.2d	\n\t"\
		"fadd	v5.2d ,v0.2d,v5.2d		\n\t	fsub	v10.2d,v10.2d,v15.2d	\n\t"\
		"fsub	v1.2d ,v1.2d,v4.2d		\n\t	fsub	v11.2d,v11.2d,v14.2d	\n\t"\
		/* Lcol Output order: 6,7,2,3,0,4,5,1: */\
		"stp	q18,q19,[x0]			\n\t	fsub	v22.2d,v8.2d ,v12.2d	\n\t"\
		"stp	q2 ,q3 ,[x1]			\n\t	fsub	v23.2d,v9.2d ,v13.2d	\n\t"\
		"stp	q20,q21,[x2]			\n\t	fadd	v12.2d,v8.2d ,v12.2d	\n\t"\
		"stp	q5 ,q1 ,[x3]			\n\t	fadd	v13.2d,v9.2d ,v13.2d	\n\t"\
		"ldr	w9,%[__p08]		\n\t"\
		/* Rcol Output order: 12,13,8,9,10,14,15,11; Fiddle order of Rcol stores to reflect order in which outputs computed: */\
		"add	x0,x0,x9 ,lsl #3		\n\t	stp	q10,q24,[x6]		\n\t"/* lcol: add0 + p24 */\
		"add	x1,x1,x9 ,lsl #3		\n\t	stp	q25,q11,[x7]		\n\t"/* lcol: add0 + p25 */\
		"add	x2,x2,x9 ,lsl #3		\n\t	stp	q22,q23,[x5]		\n\t"/* lcol: add0 + p26 */\
		"add	x3,x3,x9 ,lsl #3		\n\t	stp	q12,q13,[x4]		\n\t"/* lcol: add0 + p27 */\
	/*...Block 4: t06,t16,t26,t36: */		/*...Block 8: t0E,t1E,t2E,t3E: */\
		"add	x14,x14,#0x40			\n\t"/* r06 */\
		"ldp	q16,q17,[x10,#0x40]		\n\t	ldp	q14,q15,[x10,#0x20]		\n\t"/* cc3,cc1 */\
												/* p04 still in w8: */\
		"ldp	q0 ,q1 ,[x14       ]	\n\t	add	x4,x4,x9,lsl #3			\n\t"/* rcol: add0 + p28 */\
		"ldp	q2 ,q3 ,[x14,#0x100]	\n\t	add	x5,x5,x9,lsl #3			\n\t"/* rcol: add0 + p29 */\
		"ldp	q4 ,q5 ,[x14,#0x200]	\n\t	add	x6,x6,x9,lsl #3			\n\t"/* rcol: add0 + p30 */\
		"ldp	q20,q21,[x14,#0x300]	\n\t	add	x7,x7,x9,lsl #3			\n\t"/* rcol: add0 + p31 */\
		"fmul	v18.2d,v4.2d,v16.2d		\n\t	ldp	q8 ,q9 ,[x14,#0x080]	\n\t"\
		"fmul	v19.2d,v5.2d,v16.2d		\n\t	ldp	q10,q11,[x14,#0x180]	\n\t"\
		"fmls	v18.2d,v5.2d,v17.2d		\n\t	ldp	q12,q13,[x14,#0x280]	\n\t"\
		"fmla	v19.2d,v4.2d,v17.2d		\n\t	ldp	q24,q25,[x14,#0x380]	\n\t"\
		/* Next set of tewiddles swap [c,s] pairs between the 2 columns: */\
		"fmul	v6.2d,v20.2d,v15.2d		\n\t	fmul	v22.2d,v12.2d,v15.2d	\n\t"\
		"fmul	v7.2d,v21.2d,v15.2d		\n\t	fmul	v23.2d,v13.2d,v15.2d	\n\t"\
		"fmla	v6.2d,v21.2d,v14.2d		\n\t	fmls	v22.2d,v13.2d,v14.2d	\n\t"\
		"fmls	v7.2d,v20.2d,v14.2d		\n\t	fmla	v23.2d,v12.2d,v14.2d	\n\t"\
		"fsub	v4.2d,v18.2d,v6.2d		\n\t	fmul	v14.2d,v24.2d,v17.2d	\n\t"\
		"fsub	v5.2d,v19.2d,v7.2d		\n\t	fmul	v15.2d,v25.2d,v17.2d	\n\t"\
		"fadd	v6.2d,v18.2d,v6.2d		\n\t	fmls	v14.2d,v25.2d,v16.2d	\n\t"\
		"fadd	v7.2d,v19.2d,v7.2d		\n\t	fmla	v15.2d,v24.2d,v16.2d	\n\t"\
		"ldp	q16,q17,[x10     ]		\n\t"/* cc0 */\
		"fmul	v18.2d,v2.2d,v17.2d		\n\t	fsub	v12.2d,v22.2d,v14.2d	\n\t"\
		"fmul	v19.2d,v3.2d,v17.2d		\n\t	fsub	v13.2d,v23.2d,v15.2d	\n\t"\
		"fmls	v18.2d,v3.2d,v16.2d		\n\t	fadd	v14.2d,v22.2d,v14.2d	\n\t"\
		"fmla	v19.2d,v2.2d,v16.2d		\n\t	fadd	v15.2d,v23.2d,v15.2d	\n\t"\
		/* Post-twiddle sub-DFT sequence *almost* [see next comment] same as pvs blockpair, just note 2,3 in v18,v19: */\
		"fadd	v2.2d,v0.2d,v18.2d		\n\t	fmul	v22.2d,v10.2d,v16.2d	\n\t"\
		"fadd	v3.2d,v1.2d,v19.2d		\n\t	fmul	v23.2d,v11.2d,v16.2d	\n\t"\
		"fsub	v0.2d,v0.2d,v18.2d		\n\t	fmla	v22.2d,v11.2d,v17.2d	\n\t"\
		"fsub	v1.2d,v1.2d,v19.2d		\n\t	fmls	v23.2d,v10.2d,v17.2d	\n\t"\
	/***** Key difference vs pvs 2 blockpairs - In lcol we swap reg-apirs [4,5] <--> [6,7]: *****/\
		"fadd	v18.2d,v2.2d,v4.2d		\n\t	fadd	v10.2d,v8.2d,v22.2d		\n\t"\
		"fadd	v19.2d,v3.2d,v5.2d		\n\t	fadd	v11.2d,v9.2d,v23.2d		\n\t"\
		"fsub	v2.2d ,v2.2d,v4.2d		\n\t	fsub	v8.2d ,v8.2d,v22.2d		\n\t"\
		"fsub	v3.2d ,v3.2d,v5.2d		\n\t	fsub	v9.2d ,v9.2d,v23.2d		\n\t"\
		"fsub	v20.2d,v0.2d,v7.2d		\n\t	fadd	v25.2d,v10.2d,v15.2d	\n\t"\
		"fadd	v21.2d,v1.2d,v6.2d		\n\t	fadd	v24.2d,v11.2d,v14.2d	\n\t"\
		"fadd	v5.2d ,v0.2d,v7.2d		\n\t	fsub	v10.2d,v10.2d,v15.2d	\n\t"\
		"fsub	v1.2d ,v1.2d,v6.2d		\n\t	fsub	v11.2d,v11.2d,v14.2d	\n\t"\
		/* Lcol Output order: 6,7,2,3,0,4,5,1: */\
		"stp	q18,q19,[x0]			\n\t	fsub	v22.2d,v8.2d ,v12.2d	\n\t"\
		"stp	q2 ,q3 ,[x1]			\n\t	fsub	v23.2d,v9.2d ,v13.2d	\n\t"\
		"stp	q20,q21,[x2]			\n\t	fadd	v12.2d,v8.2d ,v12.2d	\n\t"\
		"stp	q5 ,q1 ,[x3]			\n\t	fadd	v13.2d,v9.2d ,v13.2d	\n\t"\
		"ldr	x0,%[__add0]			\n\t	ldr	w9,%[__p08]	\n\t"\
		/* Rcol Output order: 12,13,8,9,10,14,15,11; Fiddle order of Rcol stores to reflect order in which outputs computed: */\
		"										stp	q10,q24,[x6]		\n\t"/* lcol: add0 + p24 */\
		"										stp	q25,q11,[x7]		\n\t"/* lcol: add0 + p25 */\
		"										stp	q22,q23,[x5]		\n\t"/* lcol: add0 + p26 */\
		"										stp	q12,q13,[x4]		\n\t"/* lcol: add0 + p27 */\
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
		: "cc","memory","x0","x1","x2","x3","x4","x5","x6","x7","x8","x9","x10","x11","x12","x13","x14","x15",\
		"v0","v1","v2","v3","v4","v5","v6","v7","v8","v9","v10","v11","v12","v13","v14","v15",\
		"v16","v17","v18","v19","v20","v21","v22","v23","v24","v25", "v29"	/* Clobbered registers */\
	);\
	}

	#if 0	// twiddles mem-map:
					cc[00 01 02 03 04 05 06 07 08 09 0A 0B 0C 0D 0E 0F 10 11 12 13 14 15 16 17 18 19 1A 1B 1C 1D 1E 1F]
		(cc0,ss0) + 0x[06,26,16,36|0e,2e,1e,3e|08,28,18,38|10,30,20,40|0a,2a,1a,3a|12,32,22,42|0c,2c,1c,3c|14,34,24,44].
	*** NOTE: This is the same pattern as for DIF, but with the middle 2 roots octets [08-40] and [0a-42] swapped ***
	#endif
	// Vector-opcount: 86 LDP, 48 STP, 236 ADD, 148 MUL; compare to above 32-DIF opcount 95 LDP, 64 STP, 280 ADD, 176 MUL:
	#define SSE2_RADIX32_DIT_TWIDDLE(Xadd0,Xp01,Xp02,Xp03,Xp04,Xp05,Xp06,Xp07,Xp08,Xp10,Xr00,Xisrt2)\
	{\
	__asm__ volatile (\
		"ldr	x14,%[__r00]		\n\t	add	x10, x14,#0x400			\n\t"/* x10 = &isrt2 */\
		"ldr	w11,=0x10			\n\t	ld1r	{v29.2d},[x10],x11	\n\t"/* x10 has isrt2-ptr for LD1R, gets post-incr'ed by 0x10 to cc0 */\
		"ldr	x0,%[__add0]		\n\t	ldr	w8,%[__p04]				\n\t"\
		"ldr	w11,%[__p01]		\n\t	add	x1, x0,x11,lsl #3		\n\t"/* add0 + p1 */\
		"ldr	w12,%[__p02]		\n\t	add	x2, x0,x12,lsl #3		\n\t"/* add0 + p2 */\
		"ldr	w13,%[__p03]		\n\t	add	x3, x0,x13,lsl #3		\n\t"/* add0 + p3 */\
		"ldr	w9 ,%[__p08]		\n\t	add	x4, x0,x8 ,lsl #3		\n\t"/* add0 + p4 */\
		"									add	x5, x1,x8 ,lsl #3		\n\t"/* add0 + p5 */\
		"									add	x6, x2,x8 ,lsl #3		\n\t"/* add0 + p6 */\
		"									add	x7, x3,x8 ,lsl #3		\n\t"/* add0 + p7 */\
	/*...Block 1: Inputs from add0 + p[0-7] (with resp. 4-DITs using the [0-3],[4-7] subsets of same), Outputs into r00 + [0-15]: */\
		"ldp	q0,q1,[x0]			\n\t	ldp	q12,q13,[x4]			\n\t"\
		"ldp	q8,q9,[x1]			\n\t	ldp	q20,q21,[x5]			\n\t"\
		"fsub	v2.2d,v0.2d,v8.2d	\n\t	fsub	v14.2d,v12.2d,v20.2d	\n\t"\
		"fsub	v3.2d,v1.2d,v9.2d	\n\t	fsub	v15.2d,v13.2d,v21.2d	\n\t"\
		"fadd	v0.2d,v0.2d,v8.2d	\n\t	fadd	v12.2d,v12.2d,v20.2d	\n\t"\
		"fadd	v1.2d,v1.2d,v9.2d	\n\t	fadd	v13.2d,v13.2d,v21.2d	\n\t"\
		"ldp	q4,q5,[x2]			\n\t	ldp	q16,q17,[x6]			\n\t"\
		"ldp	q8,q9,[x3]			\n\t	ldp	q20,q21,[x7]			\n\t"\
		"fsub	v6.2d,v4.2d,v8.2d	\n\t	fsub	v18.2d,v16.2d,v20.2d	\n\t"\
		"fsub	v7.2d,v5.2d,v9.2d	\n\t	fsub	v19.2d,v17.2d,v21.2d	\n\t"\
		"fadd	v4.2d,v4.2d,v8.2d	\n\t	fadd	v16.2d,v16.2d,v20.2d	\n\t"\
		"fadd	v5.2d,v5.2d,v9.2d	\n\t	fadd	v17.2d,v17.2d,v21.2d	\n\t"\
		/* combine to get the 2 length-4 transforms: */\
		"fsub	v8.2d,v0.2d,v4.2d	\n\t	fsub	v20.2d,v12.2d,v16.2d	\n\t"/* t5,13 */ \
		"fsub	v9.2d,v1.2d,v5.2d	\n\t	fsub	v21.2d,v13.2d,v17.2d	\n\t"/* t6,14 */ \
		"fadd	v0.2d,v0.2d,v4.2d	\n\t	fadd	v12.2d,v12.2d,v16.2d	\n\t"/* t1, 9 */ \
		"fadd	v1.2d,v1.2d,v5.2d	\n\t	fadd	v13.2d,v13.2d,v17.2d	\n\t"/* t2,10 */ \
		"fadd	v4.2d,v2.2d,v7.2d	\n\t	fadd	v16.2d,v14.2d,v19.2d	\n\t"/* t3,11 */\
		"fsub	v5.2d,v3.2d,v6.2d	\n\t	fsub	v17.2d,v15.2d,v18.2d	\n\t"/* t4,12 */\
		"fsub	v2.2d,v2.2d,v7.2d	\n\t	fsub	v14.2d,v14.2d,v19.2d	\n\t"/* t7,15 */\
		"fadd	v3.2d,v3.2d,v6.2d	\n\t	fadd	v15.2d,v15.2d,v18.2d	\n\t"/* t8,16 */\
	/* now combine the two half-transforms: */\
	  /* [0]	a[j1   ]=t1+t9;				a[j2   ]=t2+t10;
				t1      =t1-t9;				t2      =t2-t10; */\
		"fadd	v6.2d,v0.2d,v12.2d	\n\t	fadd	v18.2d,v1.2d,v13.2d	\n\t"\
		"fsub	v7.2d,v0.2d,v12.2d	\n\t	fsub	v19.2d,v1.2d,v13.2d	\n\t"\
		"stp	q6,q18,[x14]		\n\t	stp	q7,q19,[x14,#0x80]		\n\t"\
	  /* [1]	rt=(t11+t12)*ISRT2;			it=(t11-t12)*ISRT2;
				t11     =t3+rt;				t12       =t4-it;
				t3      =t3-rt;				t4        =t4+it; */\
		"fadd	v6.2d,v16.2d,v17.2d	\n\t	fsub	v7.2d,v16.2d,v17.2d	\n\t"/* t3,4 = v4,5; t11,12 = v16,17 */\
		"fmul	v16.2d,v29.2d,v6.2d	\n\t	fmul	v17.2d,v29.2d,v7.2d	\n\t"/* rt,it; */\
		"fadd	v6.2d,v4.2d,v16.2d	\n\t	fsub	v18.2d,v5.2d,v17.2d	\n\t"/* t11,12 = v6,18 */\
		"fsub	v7.2d,v4.2d,v16.2d	\n\t	fadd	v19.2d,v5.2d,v17.2d	\n\t"/* t3 ,4  = v7,19 */\
		"stp	q6,q18,[x14,#0x20]	\n\t	stp	q7,q19,[x14,#0xa0]		\n\t"\
	  /*rt      =t5+t14;			it        =t6-t13;
		t5      =t5-t14;			t6        =t6+t13;	ARM: On input, t5,6 = v8,9; t13,14 = v20,21: */\
		"fadd	v6.2d,v8.2d,v21.2d	\n\t	fsub	v18.2d,v9.2d,v20.2d	\n\t"/* rt,it = v6,18 */\
		"fsub	v7.2d,v8.2d,v21.2d	\n\t	fadd	v19.2d,v9.2d,v20.2d	\n\t"/* t3,4  = v7,19 */\
		"stp	q6,q18,[x14,#0x40]	\n\t	stp	q7,q19,[x14,#0xc0]		\n\t"\
	  /*rt=(t15-t16)*ISRT2;			it=(t15+t16)*ISRT2;
		t15     =t7-rt;				t16       =t8-it;
		t7      =t7+rt;				t8        =t8+it; */\
		"fsub	v6.2d,v14.2d,v15.2d	\n\t	fadd	v7.2d,v14.2d,v15.2d	\n\t"/* t7,8 = v2,3; t15,16 = v14,15 */\
		"fmul	v14.2d,v29.2d,v6.2d	\n\t	fmul	v15.2d,v29.2d,v7.2d	\n\t"/* rt,it; */\
		"fsub	v6.2d,v2.2d,v14.2d	\n\t	fsub	v18.2d,v3.2d,v15.2d	\n\t"/* t15,16 = v6,18 */\
		"fadd	v7.2d,v2.2d,v14.2d	\n\t	fadd	v19.2d,v3.2d,v15.2d	\n\t"/* t7 ,8  = v7,19 */\
		"stp	q6,q18,[x14,#0x60]	\n\t	stp	q7,q19,[x14,#0xe0]		\n\t"\
		/* */\
	/*...Block 2: Inputs from add0 + p[0-7] + p8, Outputs into r00 + [0-15] + 0x100: */\
		/* */\
		"add	x14,x14,#0x100		\n\t"/* r10 */\
		"add	x0, x0,x9,lsl #3	\n\t"\
		"add	x1, x1,x9,lsl #3	\n\t"\
		"add	x2, x2,x9,lsl #3	\n\t"\
		"add	x3, x3,x9,lsl #3	\n\t"\
		"add	x4, x4,x9,lsl #3	\n\t"\
		"add	x5, x5,x9,lsl #3	\n\t"\
		"add	x6, x6,x9,lsl #3	\n\t"\
		"add	x7, x7,x9,lsl #3	\n\t"\
		"ldp	q0,q1,[x0]			\n\t	ldp	q12,q13,[x4]			\n\t"\
		"ldp	q8,q9,[x1]			\n\t	ldp	q20,q21,[x5]			\n\t"\
		"fsub	v2.2d,v0.2d,v8.2d	\n\t	fsub	v14.2d,v12.2d,v20.2d	\n\t"\
		"fsub	v3.2d,v1.2d,v9.2d	\n\t	fsub	v15.2d,v13.2d,v21.2d	\n\t"\
		"fadd	v0.2d,v0.2d,v8.2d	\n\t	fadd	v12.2d,v12.2d,v20.2d	\n\t"\
		"fadd	v1.2d,v1.2d,v9.2d	\n\t	fadd	v13.2d,v13.2d,v21.2d	\n\t"\
		"ldp	q4,q5,[x2]			\n\t	ldp	q16,q17,[x6]			\n\t"\
		"ldp	q8,q9,[x3]			\n\t	ldp	q20,q21,[x7]			\n\t"\
		"fsub	v6.2d,v4.2d,v8.2d	\n\t	fsub	v18.2d,v16.2d,v20.2d	\n\t"\
		"fsub	v7.2d,v5.2d,v9.2d	\n\t	fsub	v19.2d,v17.2d,v21.2d	\n\t"\
		"fadd	v4.2d,v4.2d,v8.2d	\n\t	fadd	v16.2d,v16.2d,v20.2d	\n\t"\
		"fadd	v5.2d,v5.2d,v9.2d	\n\t	fadd	v17.2d,v17.2d,v21.2d	\n\t"\
		/* combine to get the 2 length-4 transforms: */\
		"fsub	v8.2d,v0.2d,v4.2d	\n\t	fsub	v20.2d,v12.2d,v16.2d	\n\t"\
		"fsub	v9.2d,v1.2d,v5.2d	\n\t	fsub	v21.2d,v13.2d,v17.2d	\n\t"\
		"fadd	v0.2d,v0.2d,v4.2d	\n\t	fadd	v12.2d,v12.2d,v16.2d	\n\t"\
		"fadd	v1.2d,v1.2d,v5.2d	\n\t	fadd	v13.2d,v13.2d,v17.2d	\n\t"\
		"fadd	v4.2d,v2.2d,v7.2d	\n\t	fadd	v16.2d,v14.2d,v19.2d	\n\t"\
		"fsub	v5.2d,v3.2d,v6.2d	\n\t	fsub	v17.2d,v15.2d,v18.2d	\n\t"\
		"fsub	v2.2d,v2.2d,v7.2d	\n\t	fsub	v14.2d,v14.2d,v19.2d	\n\t"\
		"fadd	v3.2d,v3.2d,v6.2d	\n\t	fadd	v15.2d,v15.2d,v18.2d	\n\t"\
		/* now combine the two half-transforms: */\
		"fadd	v6.2d,v0.2d,v12.2d	\n\t	fadd	v18.2d,v1.2d,v13.2d	\n\t"\
		"fsub	v7.2d,v0.2d,v12.2d	\n\t	fsub	v19.2d,v1.2d,v13.2d	\n\t"\
		"stp	q6,q18,[x14]		\n\t	stp	q7,q19,[x14,#0x80]		\n\t"\
		"fadd	v6.2d,v16.2d,v17.2d	\n\t	fsub	v7.2d,v16.2d,v17.2d	\n\t"\
		"fmul	v16.2d,v29.2d,v6.2d	\n\t	fmul	v17.2d,v29.2d,v7.2d	\n\t"\
		"fadd	v6.2d,v4.2d,v16.2d	\n\t	fsub	v18.2d,v5.2d,v17.2d	\n\t"\
		"fsub	v7.2d,v4.2d,v16.2d	\n\t	fadd	v19.2d,v5.2d,v17.2d	\n\t"\
		"stp	q6,q18,[x14,#0x20]	\n\t	stp	q7,q19,[x14,#0xa0]		\n\t"\
		"fadd	v6.2d,v8.2d,v21.2d	\n\t	fsub	v18.2d,v9.2d,v20.2d	\n\t"\
		"fsub	v7.2d,v8.2d,v21.2d	\n\t	fadd	v19.2d,v9.2d,v20.2d	\n\t"\
		"stp	q6,q18,[x14,#0x40]	\n\t	stp	q7,q19,[x14,#0xc0]		\n\t"\
		"fsub	v6.2d,v14.2d,v15.2d	\n\t	fadd	v7.2d,v14.2d,v15.2d	\n\t"\
		"fmul	v14.2d,v29.2d,v6.2d	\n\t	fmul	v15.2d,v29.2d,v7.2d	\n\t"\
		"fsub	v6.2d,v2.2d,v14.2d	\n\t	fsub	v18.2d,v3.2d,v15.2d	\n\t"\
		"fadd	v7.2d,v2.2d,v14.2d	\n\t	fadd	v19.2d,v3.2d,v15.2d	\n\t"\
		"stp	q6,q18,[x14,#0x60]	\n\t	stp	q7,q19,[x14,#0xe0]		\n\t"\
		/* */\
	/*...Block 3: Inputs from add0 + p[0-7] + p16, Outputs into r00 + [0-15] + 0x200: */\
		/* */\
		"add	x14,x14,#0x100		\n\t"/* r20 */\
		"ldr	x0,%[__add0]		\n\t	ldr	w9,%[__p10]	\n\t"/* Due to maths of padded-indexing, can't just incr. by another += p08 */\
		"add	x0,x0,x9 ,lsl #3		\n\t"/* add0 + p16 */\
		"add	x1,x0,x11,lsl #3		\n\t"/* add0 + p17 */\
		"add	x2,x0,x12,lsl #3		\n\t"/* add0 + p18 */\
		"add	x3,x0,x13,lsl #3		\n\t"/* add0 + p19 */\
		"add	x4,x0,x8 ,lsl #3		\n\t"/* add0 + p20 */\
		"add	x5,x1,x8 ,lsl #3		\n\t"/* add0 + p21 */\
		"add	x6,x2,x8 ,lsl #3		\n\t"/* add0 + p22 */\
		"add	x7,x3,x8 ,lsl #3		\n\t"/* add0 + p23 */\
		"ldp	q0,q1,[x0]			\n\t	ldp	q12,q13,[x4]			\n\t"\
		"ldp	q8,q9,[x1]			\n\t	ldp	q20,q21,[x5]			\n\t"\
		"fsub	v2.2d,v0.2d,v8.2d	\n\t	fsub	v14.2d,v12.2d,v20.2d	\n\t"\
		"fsub	v3.2d,v1.2d,v9.2d	\n\t	fsub	v15.2d,v13.2d,v21.2d	\n\t"\
		"fadd	v0.2d,v0.2d,v8.2d	\n\t	fadd	v12.2d,v12.2d,v20.2d	\n\t"\
		"fadd	v1.2d,v1.2d,v9.2d	\n\t	fadd	v13.2d,v13.2d,v21.2d	\n\t"\
		"ldp	q4,q5,[x2]			\n\t	ldp	q16,q17,[x6]			\n\t"\
		"ldp	q8,q9,[x3]			\n\t	ldp	q20,q21,[x7]			\n\t"\
		"fsub	v6.2d,v4.2d,v8.2d	\n\t	fsub	v18.2d,v16.2d,v20.2d	\n\t"\
		"fsub	v7.2d,v5.2d,v9.2d	\n\t	fsub	v19.2d,v17.2d,v21.2d	\n\t"\
		"fadd	v4.2d,v4.2d,v8.2d	\n\t	fadd	v16.2d,v16.2d,v20.2d	\n\t"\
		"fadd	v5.2d,v5.2d,v9.2d	\n\t	fadd	v17.2d,v17.2d,v21.2d	\n\t"\
		/* combine to get the 2 length-4 transforms: */\
		"fsub	v8.2d,v0.2d,v4.2d	\n\t	fsub	v20.2d,v12.2d,v16.2d	\n\t"\
		"fsub	v9.2d,v1.2d,v5.2d	\n\t	fsub	v21.2d,v13.2d,v17.2d	\n\t"\
		"fadd	v0.2d,v0.2d,v4.2d	\n\t	fadd	v12.2d,v12.2d,v16.2d	\n\t"\
		"fadd	v1.2d,v1.2d,v5.2d	\n\t	fadd	v13.2d,v13.2d,v17.2d	\n\t"\
		"fadd	v4.2d,v2.2d,v7.2d	\n\t	fadd	v16.2d,v14.2d,v19.2d	\n\t"\
		"fsub	v5.2d,v3.2d,v6.2d	\n\t	fsub	v17.2d,v15.2d,v18.2d	\n\t"\
		"fsub	v2.2d,v2.2d,v7.2d	\n\t	fsub	v14.2d,v14.2d,v19.2d	\n\t"\
		"fadd	v3.2d,v3.2d,v6.2d	\n\t	fadd	v15.2d,v15.2d,v18.2d	\n\t"\
		/* now combine the two half-transforms: */\
		"fadd	v6.2d,v0.2d,v12.2d	\n\t	fadd	v18.2d,v1.2d,v13.2d	\n\t"\
		"fsub	v7.2d,v0.2d,v12.2d	\n\t	fsub	v19.2d,v1.2d,v13.2d	\n\t"\
		"stp	q6,q18,[x14]		\n\t	stp	q7,q19,[x14,#0x80]		\n\t"\
		"fadd	v6.2d,v16.2d,v17.2d	\n\t	fsub	v7.2d,v16.2d,v17.2d	\n\t"\
		"fmul	v16.2d,v29.2d,v6.2d	\n\t	fmul	v17.2d,v29.2d,v7.2d	\n\t"\
		"fadd	v6.2d,v4.2d,v16.2d	\n\t	fsub	v18.2d,v5.2d,v17.2d	\n\t"\
		"fsub	v7.2d,v4.2d,v16.2d	\n\t	fadd	v19.2d,v5.2d,v17.2d	\n\t"\
		"stp	q6,q18,[x14,#0x20]	\n\t	stp	q7,q19,[x14,#0xa0]		\n\t"\
		"fadd	v6.2d,v8.2d,v21.2d	\n\t	fsub	v18.2d,v9.2d,v20.2d	\n\t"\
		"fsub	v7.2d,v8.2d,v21.2d	\n\t	fadd	v19.2d,v9.2d,v20.2d	\n\t"\
		"stp	q6,q18,[x14,#0x40]	\n\t	stp	q7,q19,[x14,#0xc0]		\n\t"\
		"fsub	v6.2d,v14.2d,v15.2d	\n\t	fadd	v7.2d,v14.2d,v15.2d	\n\t"\
		"fmul	v14.2d,v29.2d,v6.2d	\n\t	fmul	v15.2d,v29.2d,v7.2d	\n\t"\
		"fsub	v6.2d,v2.2d,v14.2d	\n\t	fsub	v18.2d,v3.2d,v15.2d	\n\t"\
		"fadd	v7.2d,v2.2d,v14.2d	\n\t	fadd	v19.2d,v3.2d,v15.2d	\n\t"\
		"stp	q6,q18,[x14,#0x60]	\n\t	stp	q7,q19,[x14,#0xe0]		\n\t"\
		/* */\
	/*...Block 4: Inputs from add0 + p[0-7] + p24, Outputs into r00 + [0-15] + 0x300: */\
		/* */\
		"add	x14,x14,#0x100		\n\t"/* r30 */\
		"ldr	w9,%[__p08]			\n\t"\
		"add	x0, x0,x9,lsl #3	\n\t"\
		"add	x1, x1,x9,lsl #3	\n\t"\
		"add	x2, x2,x9,lsl #3	\n\t"\
		"add	x3, x3,x9,lsl #3	\n\t"\
		"add	x4, x4,x9,lsl #3	\n\t"\
		"add	x5, x5,x9,lsl #3	\n\t"\
		"add	x6, x6,x9,lsl #3	\n\t"\
		"add	x7, x7,x9,lsl #3	\n\t"\
		"ldp	q0,q1,[x0]			\n\t	ldp	q12,q13,[x4]			\n\t"\
		"ldp	q8,q9,[x1]			\n\t	ldp	q20,q21,[x5]			\n\t"\
		"fsub	v2.2d,v0.2d,v8.2d	\n\t	fsub	v14.2d,v12.2d,v20.2d	\n\t"\
		"fsub	v3.2d,v1.2d,v9.2d	\n\t	fsub	v15.2d,v13.2d,v21.2d	\n\t"\
		"fadd	v0.2d,v0.2d,v8.2d	\n\t	fadd	v12.2d,v12.2d,v20.2d	\n\t"\
		"fadd	v1.2d,v1.2d,v9.2d	\n\t	fadd	v13.2d,v13.2d,v21.2d	\n\t"\
		"ldp	q4,q5,[x2]			\n\t	ldp	q16,q17,[x6]			\n\t"\
		"ldp	q8,q9,[x3]			\n\t	ldp	q20,q21,[x7]			\n\t"\
		"fsub	v6.2d,v4.2d,v8.2d	\n\t	fsub	v18.2d,v16.2d,v20.2d	\n\t"\
		"fsub	v7.2d,v5.2d,v9.2d	\n\t	fsub	v19.2d,v17.2d,v21.2d	\n\t"\
		"fadd	v4.2d,v4.2d,v8.2d	\n\t	fadd	v16.2d,v16.2d,v20.2d	\n\t"\
		"fadd	v5.2d,v5.2d,v9.2d	\n\t	fadd	v17.2d,v17.2d,v21.2d	\n\t"\
		/* combine to get the 2 length-4 transforms: */\
		"fsub	v8.2d,v0.2d,v4.2d	\n\t	fsub	v20.2d,v12.2d,v16.2d	\n\t"\
		"fsub	v9.2d,v1.2d,v5.2d	\n\t	fsub	v21.2d,v13.2d,v17.2d	\n\t"\
		"fadd	v0.2d,v0.2d,v4.2d	\n\t	fadd	v12.2d,v12.2d,v16.2d	\n\t"\
		"fadd	v1.2d,v1.2d,v5.2d	\n\t	fadd	v13.2d,v13.2d,v17.2d	\n\t"\
		"fadd	v4.2d,v2.2d,v7.2d	\n\t	fadd	v16.2d,v14.2d,v19.2d	\n\t"\
		"fsub	v5.2d,v3.2d,v6.2d	\n\t	fsub	v17.2d,v15.2d,v18.2d	\n\t"\
		"fsub	v2.2d,v2.2d,v7.2d	\n\t	fsub	v14.2d,v14.2d,v19.2d	\n\t"\
		"fadd	v3.2d,v3.2d,v6.2d	\n\t	fadd	v15.2d,v15.2d,v18.2d	\n\t"\
		/* now combine the two half-transforms: */\
		"fadd	v6.2d,v0.2d,v12.2d	\n\t	fadd	v18.2d,v1.2d,v13.2d	\n\t"\
		"fsub	v7.2d,v0.2d,v12.2d	\n\t	fsub	v19.2d,v1.2d,v13.2d	\n\t"\
		"stp	q6,q18,[x14]		\n\t	stp	q7,q19,[x14,#0x80]		\n\t"\
		"fadd	v6.2d,v16.2d,v17.2d	\n\t	fsub	v7.2d,v16.2d,v17.2d	\n\t"\
		"fmul	v16.2d,v29.2d,v6.2d	\n\t	fmul	v17.2d,v29.2d,v7.2d	\n\t"\
		"fadd	v6.2d,v4.2d,v16.2d	\n\t	fsub	v18.2d,v5.2d,v17.2d	\n\t"\
		"fsub	v7.2d,v4.2d,v16.2d	\n\t	fadd	v19.2d,v5.2d,v17.2d	\n\t"\
		"stp	q6,q18,[x14,#0x20]	\n\t	stp	q7,q19,[x14,#0xa0]		\n\t"\
		"fadd	v6.2d,v8.2d,v21.2d	\n\t	fsub	v18.2d,v9.2d,v20.2d	\n\t"\
		"fsub	v7.2d,v8.2d,v21.2d	\n\t	fadd	v19.2d,v9.2d,v20.2d	\n\t"\
		"stp	q6,q18,[x14,#0x40]	\n\t	stp	q7,q19,[x14,#0xc0]		\n\t"\
		"fsub	v6.2d,v14.2d,v15.2d	\n\t	fadd	v7.2d,v14.2d,v15.2d	\n\t"\
		"fmul	v14.2d,v29.2d,v6.2d	\n\t	fmul	v15.2d,v29.2d,v7.2d	\n\t"\
		"fsub	v6.2d,v2.2d,v14.2d	\n\t	fsub	v18.2d,v3.2d,v15.2d	\n\t"\
		"fadd	v7.2d,v2.2d,v14.2d	\n\t	fadd	v19.2d,v3.2d,v15.2d	\n\t"\
		"stp	q6,q18,[x14,#0x60]	\n\t	stp	q7,q19,[x14,#0xe0]		\n\t"\
	/**********************************************************************************/\
	/*...and now do eight radix-4 transforms, including the internal twiddle factors: */\
	/**********************************************************************************/\
		/*...Block 1: t00,t10,t20,t30	*/		/*...Block 5: t08,t18,t28,t38*/\
		"sub	x14,x14,#0x300		\n\t"/* r00; p04 still in w8 */\
		"ldr	w11,%[__p08]		\n\t"\
		"ldr	w12,%[__p10]		\n\t"\
		"add	w13,w12,w11		\n\t"/* __p18 */\
		"ldr	x0,%[__add0]			\n\t	ldp	q8 ,q9 ,[x14,#0x080]	\n\t"/* lcol: add0 + p0  */\
		"add	x1,x0,x11,lsl #3		\n\t	ldp	q10,q11,[x14,#0x180]	\n\t"/* lcol: add0 + p8  */\
		"add	x2,x0,x12,lsl #3		\n\t	ldp	q12,q20,[x14,#0x280]	\n\t"/* lcol: add0 + p16 */\
		"add	x3,x0,x13,lsl #3		\n\t	ldp	q21,q15,[x14,#0x380]	\n\t"/* lcol: add0 + p24 */\
		"ldp	q0,q1,[x14       ]		\n\t	add	x4,x0,x8,lsl #3			\n\t"/* rcol: add0 + p4  */\
		"ldp	q2,q3,[x14,#0x100]		\n\t	add	x5,x1,x8,lsl #3			\n\t"/* rcol: add0 + p12 */\
		"ldp	q4,q5,[x14,#0x200]		\n\t	add	x6,x2,x8,lsl #3			\n\t"/* rcol: add0 + p20 */\
		"ldp	q6,q7,[x14,#0x300]		\n\t	add	x7,x3,x8,lsl #3			\n\t"/* rcol: add0 + p28 */\
/* 2 */	"fadd	v16.2d,v0.2d,v2.2d		\n\t	fmul	v12.2d,v12.2d,v29.2d	\n\t"/* c *= isrt2 */\
/* 3 */	"fadd	v17.2d,v1.2d,v3.2d		\n\t	fmul	v20.2d,v20.2d,v29.2d	\n\t"/* d *= isrt2 */\
/* 6 */	"fadd	v18.2d,v4.2d,v6.2d		\n\t	fmul	v21.2d,v21.2d,v29.2d	\n\t"/* e *= isrt2 */\
/* 7 */	"fadd	v19.2d,v5.2d,v7.2d		\n\t	fmul	v15.2d,v15.2d,v29.2d	\n\t"/* f *= isrt2 */\
		"fsub	v0.2d ,v0.2d,v2.2d		\n\t	fsub	v13.2d,v20.2d,v12.2d	\n\t"/* d=d-c */\
		"fsub	v1.2d ,v1.2d,v3.2d		\n\t	fsub	v14.2d,v21.2d,v15.2d	\n\t"/* e=e-f */\
		"fsub	v4.2d ,v4.2d,v6.2d		\n\t	fadd	v12.2d,v20.2d,v12.2d	\n\t"/* c=d+c */\
		"fsub	v5.2d ,v5.2d,v7.2d		\n\t	fadd	v15.2d,v21.2d,v15.2d	\n\t"/* f=e+f */\
/*2=2-6*/"fsub	v2.2d,v16.2d,v18.2d		\n\t	fsub	v20.2d,v8.2d ,v11.2d	\n\t"/* 8=8-b */\
/*3=3-7*/"fsub	v3.2d,v17.2d,v19.2d		\n\t	fsub	v21.2d,v9.2d ,v10.2d	\n\t"/* 9=9-a */\
/*6=2+6*/"fadd	v6.2d,v16.2d,v18.2d		\n\t	fsub	v22.2d,v12.2d,v14.2d	\n\t"/* c=c-e */\
/*7=3+7*/"fadd	v7.2d,v17.2d,v19.2d		\n\t	fsub	v23.2d,v13.2d,v15.2d	\n\t"/* d=d-f */\
		"										fadd	v11.2d,v8.2d ,v11.2d	\n\t"/* b=8+b */\
		"										fadd	v10.2d,v9.2d ,v10.2d	\n\t"/* a=9+a */\
		"										fadd	v14.2d,v12.2d,v14.2d	\n\t"/* e=c+e */\
/* c08 */"ldp	q18,q19,[x10,#0x80]		\n\t	fadd	v15.2d,v13.2d,v15.2d	\n\t"/* f=d+f */\
/*0=0-5*/"fsub	v16.2d,v0.2d,v5.2d		\n\t	fsub	v8.2d ,v20.2d,v15.2d	\n\t"/* 8=8-f */\
/*1=1-4*/"fsub	v17.2d,v1.2d,v4.2d		\n\t	fsub	v9.2d ,v21.2d,v23.2d	\n\t"/* 9=9-d */\
/*5=0+5*/"fadd	v5.2d ,v0.2d,v5.2d		\n\t	fsub	v12.2d,v10.2d,v14.2d	\n\t"/* a=a-e */\
/*4=1+4*/"fadd	v4.2d ,v1.2d,v4.2d		\n\t	fsub	v13.2d,v11.2d,v22.2d	\n\t"/* b=b-c */\
		"										ldp	q24,q25,[x10,#0x0e0]	\n\t"/* c04 = [a,b] */\
		"										ldp	q26,q27,[x10,#0x100]	\n\t"/* c0C = [c,s] */\
		 "fmul	v0.2d,v17.2d,v18.2d		\n\t	fadd	v15.2d,v20.2d,v15.2d	\n\t"/* f=8+f */\
		 "fmul	v1.2d,v5.2d ,v18.2d		\n\t	fadd	v23.2d,v21.2d,v23.2d	\n\t"/* d=9+d */\
/*1c-5s*/"fmls	v0.2d,v5.2d ,v19.2d		\n\t	fadd	v14.2d,v10.2d,v14.2d	\n\t"/* e=a+e */\
/*1s+5c*/"fmla	v1.2d,v17.2d,v19.2d		\n\t	fadd	v22.2d,v11.2d,v22.2d	\n\t"/* c=b+c ... switch a-f to decimal 10-15 in ensuing twiddle-muls to avoid confusion with sincos data */\
/* c10=[a,b] */"ldp q5 ,q17,[x10,#0xa0]	\n\t	fmul	v10.2d,v23.2d,v24.2d	\n\t"\
/* c18=[c,s] */"ldp q18,q19,[x10,#0xc0]	\n\t	fmul	v11.2d,v22.2d,v24.2d	\n\t"\
		"stp	q6 ,q7 ,[x0]			\n\t	fmls	v10.2d,v22.2d,v25.2d	\n\t"/* 13a-12b */\
		"stp	q1 ,q0 ,[x1]			\n\t	fmla	v11.2d,v23.2d,v25.2d	\n\t"/* 13b+12a */\
		"fmul	v0.2d,v3.2d ,v5.2d		\n\t	fmul	v20.2d,v12.2d,v26.2d	\n\t"\
		"fmul	v1.2d,v2.2d ,v5.2d		\n\t	fmul	v21.2d,v15.2d,v26.2d	\n\t"\
/*3a-2b*/"fmls	v0.2d,v2.2d ,v17.2d		\n\t	fmls	v20.2d,v15.2d,v27.2d	\n\t"/* 10c-15s */\
/*3b+2a*/"fmla	v1.2d,v3.2d ,v17.2d		\n\t	fmla	v21.2d,v12.2d,v27.2d	\n\t"/* 10s+15c */\
		"fmul	v6.2d,v4.2d ,v18.2d		\n\t	ldp	q24,q25,[x10,#0x120]	\n\t"/* c14 = [a,b] */\
		"fmul	v7.2d,v16.2d,v18.2d		\n\t	ldp	q26,q27,[x10,#0x140]	\n\t"/* c1C = [c,s] */\
/*4c-0s*/"fmls	v6.2d,v16.2d,v19.2d		\n\t	stp	q11,q10,[x4]			\n\t"\
/*4s+0c*/"fmla	v7.2d,v4.2d ,v19.2d		\n\t	stp	q21,q20,[x5]			\n\t"\
		"ldr	w11,%[__p01]		\n\t"\
		"										fmul	v10.2d,v9.2d ,v24.2d	\n\t"\
		"										fmul	v11.2d,v13.2d,v24.2d	\n\t"\
		"stp	q1 ,q0 ,[x2]			\n\t	fmls	v10.2d,v13.2d,v25.2d	\n\t"/*  9a-11b */\
		"stp	q7 ,q6 ,[x3]			\n\t	fmla	v11.2d,v9.2d ,v25.2d	\n\t"/*  9b+11a */\
		"add	x0,x0,x11,lsl #3		\n\t	fmul	v20.2d,v14.2d,v26.2d	\n\t"/* lcol: add0 + p1  */\
		"add	x1,x1,x11,lsl #3		\n\t	fmul	v21.2d,v8.2d ,v26.2d	\n\t"/* lcol: add0 + p9  */\
		"										fmls	v20.2d,v8.2d ,v27.2d	\n\t"/* 14c- 8s */\
		"										fmla	v21.2d,v14.2d,v27.2d	\n\t"/* 14s+ 8c */\
		"add	x2,x2,x11,lsl #3		\n\t	stp	q11,q10,[x6]				\n\t"/* lcol: add0 + p17 */\
		"add	x3,x3,x11,lsl #3		\n\t	stp	q21,q20,[x7]				\n\t"/* lcol: add0 + p25 */\
		/*...Block 2: t02,t12,t22,t32	*/		/*...Block 6: t0A,t1A,t2A,t3A*/\
/*r02:*/"add	x14,x14,#0x20			\n\t"/* In comments, let cc1=[c,s], cc3=[a,b]: */\
		"ldp	q0,q1,[x14,#0x200]		\n\t	add	x4,x4,x11,lsl #3			\n\t"/* rcol: add0 + p5  */\
		"ldp	q2,q3,[x14,#0x300]		\n\t	add	x5,x5,x11,lsl #3			\n\t"/* rcol: add0 + p13 */\
/*cc1:*/"ldp	q18,q19,[x10,#0x20]		\n\t	add	x6,x6,x11,lsl #3			\n\t"/* rcol: add0 + p21 */\
/*cc3:*/"ldp	q20,q21,[x10,#0x40]		\n\t	add	x7,x7,x11,lsl #3			\n\t"/* rcol: add0 + p29 */\
		"fmul	v4.2d,v0.2d,v18.2d		\n\t	ldp	q8 ,q9 ,[x14,#0x280]		\n\t"\
		"fmul	v5.2d,v1.2d,v18.2d		\n\t	ldp	q10,q11,[x14,#0x380]		\n\t"\
/*4c+5s*/"fmla	v4.2d,v1.2d,v19.2d		\n\t	fmul	v12.2d,v8.2d,v21.2d		\n\t"\
/*5c-4s*/"fmls	v5.2d,v0.2d,v19.2d		\n\t	fmul	v13.2d,v9.2d,v21.2d		\n\t"\
		"fmul	v16.2d,v2.2d,v20.2d		\n\t	fmla	v12.2d,v9.2d,v20.2d		\n\t"/*12b+13a*/\
		"fmul	v17.2d,v3.2d,v20.2d		\n\t	fmls	v13.2d,v8.2d,v20.2d		\n\t"/*13b-12a*/\
/*6a+7b*/"fmla	v16.2d,v3.2d,v21.2d		\n\t	fmul	v24.2d,v10.2d,v18.2d	\n\t"\
/*7a-6b*/"fmls	v17.2d,v2.2d,v21.2d		\n\t	fmul	v25.2d,v11.2d,v18.2d	\n\t"\
		"fsub	v6.2d,v4.2d,v16.2d		\n\t	fmls	v24.2d,v11.2d,v19.2d	\n\t"/*14c-15s*/\
		"fsub	v7.2d,v5.2d,v17.2d		\n\t	fmla	v25.2d,v10.2d,v19.2d	\n\t"/*15c+14s*/\
		"fadd	v4.2d,v4.2d,v16.2d		\n\t	fadd	v14.2d,v12.2d,v24.2d	\n\t"\
		"fadd	v5.2d,v5.2d,v17.2d		\n\t	fadd	v15.2d,v13.2d,v25.2d	\n\t"\
		"ldp	q0,q1,[x14       ]		\n\t	fsub	v12.2d,v12.2d,v24.2d	\n\t"\
		"ldp	q2,q3,[x14,#0x100]		\n\t	fsub	v13.2d,v13.2d,v25.2d	\n\t"\
		/* cc0 := [c,s]: */\
		"ldp	q20,q21,[x10]			\n\t"\
		"fmul	v16.2d,v2.2d,v20.2d		\n\t	ldp	q8 ,q9 ,[x14,#0x080]		\n\t"\
		"fmul	v17.2d,v3.2d,v20.2d		\n\t	ldp	q10,q11,[x14,#0x180]		\n\t"\
/*2c+3s*/"fmla	v16.2d,v3.2d,v21.2d		\n\t	fmul	v18.2d,v10.2d,v21.2d	\n\t"\
/*3c-2s*/"fmls	v17.2d,v2.2d,v21.2d		\n\t	fmul	v19.2d,v11.2d,v21.2d	\n\t"\
		"fadd	v2.2d,v0.2d,v16.2d		\n\t	fmls	v18.2d,v11.2d,v20.2d	\n\t"/*10s-11c*/\
		"fadd	v3.2d,v1.2d,v17.2d		\n\t	fmla	v19.2d,v10.2d,v20.2d	\n\t"/*11s+10c*/\
		"fsub	v0.2d,v0.2d,v16.2d		\n\t	fsub	v10.2d,v8.2d,v18.2d		\n\t"\
		"fsub	v1.2d,v1.2d,v17.2d		\n\t	fsub	v11.2d,v9.2d,v19.2d		\n\t"\
/* 0 */	"fsub	v16.2d,v0.2d,v7.2d		\n\t	fadd	v8.2d ,v8.2d,v18.2d		\n\t"\
/* 1 */	"fsub	v17.2d,v1.2d,v6.2d		\n\t	fadd	v9.2d ,v9.2d,v19.2d		\n\t"\
/* 2 */	"fsub	v18.2d,v2.2d,v4.2d		\n\t	fsub	v20.2d,v8.2d ,v15.2d	\n\t"/*  8 */\
/* 3 */	"fsub	v19.2d,v3.2d,v5.2d		\n\t	fsub	v21.2d,v9.2d ,v14.2d	\n\t"/*  9 */\
		"fadd	v7.2d,v0.2d,v7.2d		\n\t	fsub	v22.2d,v10.2d,v12.2d	\n\t"/* 10 */\
		"fadd	v6.2d,v1.2d,v6.2d		\n\t	fsub	v23.2d,v11.2d,v13.2d	\n\t"/* 11 */\
		"fadd	v4.2d,v2.2d,v4.2d		\n\t	fadd	v15.2d,v8.2d ,v15.2d	\n\t"\
		"fadd	v5.2d,v3.2d,v5.2d		\n\t	fadd	v14.2d,v9.2d ,v14.2d	\n\t"\
		/* Let c01 = [c,s], c09 = [a,b]: */"	fadd	v12.2d,v10.2d,v12.2d	\n\t"\
		"ldp	q0,q1,[x10,#0x260]		\n\t	fadd	v13.2d,v11.2d,v13.2d	\n\t"\
		"ldp	q2,q3,[x10,#0x280]		\n\t"	/* Let c05 = [c,s], c0D = [a,b]: */\
		"fmul	v24.2d,v4.2d,v0.2d		\n\t	ldp	q8 ,q9 ,[x10,#0x2e0]		\n\t"\
		"fmul	v25.2d,v5.2d,v0.2d		\n\t	ldp	q10,q11,[x10,#0x300]		\n\t"\
/*4c+5s*/"fmla	v24.2d,v5.2d,v1.2d		\n\t	fmul	v26.2d,v12.2d,v8.2d		\n\t"\
/*5c-4s*/"fmls	v25.2d,v4.2d,v1.2d		\n\t	fmul	v27.2d,v13.2d,v8.2d		\n\t"\
		"fmul	v4.2d,v17.2d,v2.2d		\n\t	fmla	v26.2d,v13.2d,v9.2d		\n\t"/*12c+13s*/\
		"fmul	v5.2d,v7.2d ,v2.2d		\n\t	fmls	v27.2d,v12.2d,v9.2d		\n\t"/*13c-12s*/\
/*1a-7b*/"fmls	v4.2d,v7.2d ,v3.2d		\n\t	fmul	v12.2d,v21.2d,v10.2d	\n\t"\
/*7a+1b*/"fmla	v5.2d,v17.2d,v3.2d		\n\t	fmul	v13.2d,v15.2d,v10.2d	\n\t"\
		"stp	q24,q25,[x0]			\n\t	fmls	v12.2d,v15.2d,v11.2d	\n\t"/* 9a-15b*/\
		"stp	q5 ,q4 ,[x1]			\n\t	fmla	v13.2d,v21.2d,v11.2d	\n\t"/*15a+ 9b*/\
		/* Let c11 = [c,s], c19 = [a,b]: */		/* Let c15 = [c,s], c1D = [a,b]: */\
		"ldp	q0,q1,[x10,#0x2a0]		\n\t	ldr	w12,%[__p02]		\n\t"/* x15 = p02-p01; due to index-padding this diff is not nec. == p01, */\
		"ldp	q2,q3,[x10,#0x2c0]		\n\t	sub	w12,w12,w11			\n\t"/* i.e. to get add0+p02 from add0+p01, need += (p02-p01), not += p01. */\
		"fmul	v24.2d,v18.2d,v0.2d		\n\t	stp	q26,q27,[x4]			\n\t"\
		"fmul	v25.2d,v19.2d,v0.2d		\n\t	stp	q13,q12,[x5]			\n\t"\
/*2c+3s*/"fmla	v24.2d,v19.2d,v1.2d		\n\t	ldp	q8 ,q9 ,[x10,#0x320]		\n\t"\
/*3c-2s*/"fmls	v25.2d,v18.2d,v1.2d		\n\t	ldp	q10,q11,[x10,#0x340]		\n\t"\
		"fmul	v4.2d,v6.2d ,v2.2d		\n\t	fmul	v26.2d,v22.2d,v8.2d		\n\t"\
		"fmul	v5.2d,v16.2d,v2.2d		\n\t	fmul	v27.2d,v23.2d,v8.2d		\n\t"\
/*6a-0b*/"fmls	v4.2d,v16.2d,v3.2d		\n\t	fmla	v26.2d,v23.2d,v9.2d		\n\t"/*10c+11s*/\
/*0a+6b*/"fmla	v5.2d,v6.2d ,v3.2d		\n\t	fmls	v27.2d,v22.2d,v9.2d		\n\t"/*11c-10s*/\
		"stp	q24,q25,[x2]			\n\t	fmul	v12.2d,v14.2d,v10.2d	\n\t"\
		"stp	q5 ,q4 ,[x3]			\n\t	fmul	v13.2d,v20.2d,v10.2d	\n\t"\
		"add	x0,x0,x12,lsl #3		\n\t	fmls	v12.2d,v20.2d,v11.2d	\n\t"/*14a-8b*/\
		"add	x1,x1,x12,lsl #3		\n\t	fmla	v13.2d,v14.2d,v11.2d	\n\t"/*8a+14b*/\
		"add	x2,x2,x12,lsl #3		\n\t	stp	q26,q27,[x6]			\n\t"\
		"add	x3,x3,x12,lsl #3		\n\t	stp	q13,q12,[x7]			\n\t"\
		/*** Subsequent blocks have different internal-twiddles arithmetic, but external-
			twiddles final stages all identical save for offsets into the twiddles-array: ***/\
		/*...Block 3: t04,t14,t24,t34*/			/*...Block 7: t0C,t1C,t2C,t3C*/\
		"add	x14,x14,#0x20			\n\t"/* r04 */\
		"ldp	q0,q1,[x14,#0x200]		\n\t	add	x4,x4,x12,lsl #3			\n\t"/* rcol: add0 + p6  */\
		"ldp	q2,q3,[x14,#0x300]		\n\t	add	x5,x5,x12,lsl #3			\n\t"/* rcol: add0 + p14 */\
/*cc0:*/"ldp	q18,q19,[x10]			\n\t	add	x6,x6,x12,lsl #3			\n\t"/* rcol: add0 + p22 */\
		"										add	x7,x7,x12,lsl #3			\n\t"/* rcol: add0 + p30 */\
		"fmul	v4.2d,v0.2d,v18.2d		\n\t	ldp	q8 ,q9 ,[x14,#0x280]		\n\t"\
		"fmul	v5.2d,v1.2d,v18.2d		\n\t	ldp	q10,q11,[x14,#0x380]		\n\t"\
/*4c+5s*/"fmla	v4.2d,v1.2d,v19.2d		\n\t	fmul	v12.2d,v8.2d,v19.2d		\n\t"\
/*5c-4s*/"fmls	v5.2d,v0.2d,v19.2d		\n\t	fmul	v13.2d,v9.2d,v19.2d		\n\t"\
		"fmul	v16.2d,v2.2d,v19.2d		\n\t	fmla	v12.2d,v9.2d,v18.2d		\n\t"/*12s+13c*/\
		"fmul	v17.2d,v3.2d,v19.2d		\n\t	fmls	v13.2d,v8.2d,v18.2d		\n\t"/*13s-12c*/\
/*6s+7c*/"fmla	v16.2d,v3.2d,v18.2d		\n\t	fmul	v24.2d,v10.2d,v18.2d	\n\t"\
/*7s-6c*/"fmls	v17.2d,v2.2d,v18.2d		\n\t	fmul	v25.2d,v11.2d,v18.2d	\n\t"\
		"fsub	v6.2d,v4.2d,v16.2d		\n\t	fmla	v24.2d,v11.2d,v19.2d	\n\t"/*14c+15s*/\
		"fsub	v7.2d,v5.2d,v17.2d		\n\t	fmls	v25.2d,v10.2d,v19.2d	\n\t"/*15c-14s*/\
		"fadd	v4.2d,v4.2d,v16.2d		\n\t	fadd	v14.2d,v12.2d,v24.2d	\n\t"\
		"fadd	v5.2d,v5.2d,v17.2d		\n\t	fadd	v15.2d,v13.2d,v25.2d	\n\t"\
		"ldp	q0,q1,[x14       ]		\n\t	fsub	v12.2d,v12.2d,v24.2d	\n\t"\
		"ldp	q2,q3,[x14,#0x100]		\n\t	fsub	v13.2d,v13.2d,v25.2d	\n\t"\
		"fadd	v16.2d,v3.2d,v2.2d		\n\t	ldp	q8 ,q9 ,[x14,#0x080]		\n\t"\
		"fsub	v17.2d,v3.2d,v2.2d		\n\t	ldp	q10,q11,[x14,#0x180]		\n\t"\
		"fmul	v16.2d,v16.2d,v29.2d	\n\t	fsub	v18.2d,v10.2d,v11.2d	\n\t"\
		"fmul	v17.2d,v17.2d,v29.2d	\n\t	fadd	v19.2d,v10.2d,v11.2d	\n\t"\
		"fadd	v2.2d,v0.2d,v16.2d		\n\t	fmul	v18.2d,v18.2d,v29.2d	\n\t"\
		"fadd	v3.2d,v1.2d,v17.2d		\n\t	fmul	v19.2d,v19.2d,v29.2d	\n\t"\
		"fsub	v0.2d,v0.2d,v16.2d		\n\t	fsub	v10.2d,v8.2d,v18.2d		\n\t"\
		"fsub	v1.2d,v1.2d,v17.2d		\n\t	fsub	v11.2d,v9.2d,v19.2d		\n\t"\
/* 0 */	"fsub	v16.2d,v0.2d,v7.2d		\n\t	fadd	v8.2d ,v8.2d,v18.2d		\n\t"\
/* 1 */	"fsub	v17.2d,v1.2d,v6.2d		\n\t	fadd	v9.2d ,v9.2d,v19.2d		\n\t"\
/* 2 */	"fsub	v18.2d,v2.2d,v4.2d		\n\t	fsub	v20.2d,v8.2d ,v15.2d	\n\t"/*  8 */\
/* 3 */	"fsub	v19.2d,v3.2d,v5.2d		\n\t	fsub	v21.2d,v9.2d ,v14.2d	\n\t"/*  9 */\
		"fadd	v7.2d,v0.2d,v7.2d		\n\t	fsub	v22.2d,v10.2d,v12.2d	\n\t"/* 10 */\
		"fadd	v6.2d,v1.2d,v6.2d		\n\t	fsub	v23.2d,v11.2d,v13.2d	\n\t"/* 11 */\
		"fadd	v4.2d,v2.2d,v4.2d		\n\t	fadd	v15.2d,v8.2d ,v15.2d	\n\t"\
		"fadd	v5.2d,v3.2d,v5.2d		\n\t	fadd	v14.2d,v9.2d ,v14.2d	\n\t"\
		/* Let c02 = [c,s], c0A = [a,b]: */"	fadd	v12.2d,v10.2d,v12.2d	\n\t"\
		"ldp	q0,q1,[x10,#0x160]		\n\t	fadd	v13.2d,v11.2d,v13.2d	\n\t"\
		"ldp	q2,q3,[x10,#0x180]		\n\t"	/* Let c06 = [c,s], c0E = [a,b]: */\
		"fmul	v24.2d,v4.2d,v0.2d		\n\t	ldp	q8 ,q9 ,[x10,#0x1e0]		\n\t"\
		"fmul	v25.2d,v5.2d,v0.2d		\n\t	ldp	q10,q11,[x10,#0x200]		\n\t"\
/*4c+5s*/"fmla	v24.2d,v5.2d,v1.2d		\n\t	fmul	v26.2d,v12.2d,v8.2d		\n\t"\
/*5c-4s*/"fmls	v25.2d,v4.2d,v1.2d		\n\t	fmul	v27.2d,v13.2d,v8.2d		\n\t"\
		"fmul	v4.2d,v17.2d,v2.2d		\n\t	fmla	v26.2d,v13.2d,v9.2d		\n\t"/*12c+13s*/\
		"fmul	v5.2d,v7.2d ,v2.2d		\n\t	fmls	v27.2d,v12.2d,v9.2d		\n\t"/*13c-12s*/\
/*1a-7b*/"fmls	v4.2d,v7.2d ,v3.2d		\n\t	fmul	v12.2d,v21.2d,v10.2d	\n\t"\
/*7a+1b*/"fmla	v5.2d,v17.2d,v3.2d		\n\t	fmul	v13.2d,v15.2d,v10.2d	\n\t"\
		"stp	q24,q25,[x0]			\n\t	fmls	v12.2d,v15.2d,v11.2d	\n\t"/* 9a-15b*/\
		"stp	q5 ,q4 ,[x1]			\n\t	fmla	v13.2d,v21.2d,v11.2d	\n\t"/*15a+ 9b*/\
		/* Let c12 = [c,s], c1A = [a,b]: */		/* Let c16 = [c,s], c1E = [a,b]: */\
		"ldp	q0,q1,[x10,#0x1a0]		\n\t"\
		"ldp	q2,q3,[x10,#0x1c0]		\n\t"\
		"fmul	v24.2d,v18.2d,v0.2d		\n\t	stp	q26,q27,[x4]			\n\t"\
		"fmul	v25.2d,v19.2d,v0.2d		\n\t	stp	q13,q12,[x5]			\n\t"\
/*2c+3s*/"fmla	v24.2d,v19.2d,v1.2d		\n\t	ldp	q8 ,q9 ,[x10,#0x220]		\n\t"\
/*3c-2s*/"fmls	v25.2d,v18.2d,v1.2d		\n\t	ldp	q10,q11,[x10,#0x240]		\n\t"\
		"fmul	v4.2d,v6.2d ,v2.2d		\n\t	fmul	v26.2d,v22.2d,v8.2d		\n\t"\
		"fmul	v5.2d,v16.2d,v2.2d		\n\t	fmul	v27.2d,v23.2d,v8.2d		\n\t"\
/*6a-0b*/"fmls	v4.2d,v16.2d,v3.2d		\n\t	fmla	v26.2d,v23.2d,v9.2d		\n\t"/*10c+11s*/\
/*0a+6b*/"fmla	v5.2d,v6.2d ,v3.2d		\n\t	fmls	v27.2d,v22.2d,v9.2d		\n\t"/*11c-10s*/\
		"stp	q24,q25,[x2]			\n\t	fmul	v12.2d,v14.2d,v10.2d	\n\t"\
		"stp	q5 ,q4 ,[x3]			\n\t	fmul	v13.2d,v20.2d,v10.2d	\n\t"\
		"add	x0,x0,x11,lsl #3		\n\t	fmls	v12.2d,v20.2d,v11.2d	\n\t"/*14a-8b*/\
		"add	x1,x1,x11,lsl #3		\n\t	fmla	v13.2d,v14.2d,v11.2d	\n\t"/*8a+14b*/\
		"add	x2,x2,x11,lsl #3		\n\t	stp	q26,q27,[x6]			\n\t"\
		"add	x3,x3,x11,lsl #3		\n\t	stp	q13,q12,[x7]			\n\t"\
		/*...Block 4: t06,t16,t26,t36*/			/*...Block 8: t0E,t1E,t2E,t3E*/\
		"add	x14,x14,#0x20			\n\t"/* r06 */\
		"ldp	q0,q1,[x14,#0x200]		\n\t	add	x4,x4,x11,lsl #3			\n\t"/* rcol: add0 + p7  */\
		"ldp	q2,q3,[x14,#0x300]		\n\t	add	x5,x5,x11,lsl #3			\n\t"/* rcol: add0 + p15 */\
/*cc1:=[c,s]*/"ldp	q18,q19,[x10,#0x20]	\n\t	add	x6,x6,x11,lsl #3			\n\t"/* rcol: add0 + p23 */\
/*cc3:=[a,b]*/"ldp	q20,q21,[x10,#0x40]	\n\t	add	x7,x7,x11,lsl #3			\n\t"/* rcol: add0 + p31 */\
		"fmul	v4.2d,v0.2d,v20.2d		\n\t	ldp	q8 ,q9 ,[x14,#0x280]		\n\t"\
		"fmul	v5.2d,v1.2d,v20.2d		\n\t	ldp	q10,q11,[x14,#0x380]		\n\t"\
/*4a+5b*/"fmla	v4.2d,v1.2d,v21.2d		\n\t	fmul	v12.2d,v8.2d,v19.2d		\n\t"\
/*5a-4b*/"fmls	v5.2d,v0.2d,v21.2d		\n\t	fmul	v13.2d,v9.2d,v19.2d		\n\t"\
		"fmul	v16.2d,v2.2d,v19.2d		\n\t	fmla	v12.2d,v9.2d,v18.2d		\n\t"/*12s+13c*/\
		"fmul	v17.2d,v3.2d,v19.2d		\n\t	fmls	v13.2d,v8.2d,v18.2d		\n\t"/*13s-12c*/\
/*6s-7c*/"fmls	v16.2d,v3.2d,v18.2d		\n\t	fmul	v24.2d,v10.2d,v21.2d	\n\t"\
/*7s+6c*/"fmla	v17.2d,v2.2d,v18.2d		\n\t	fmul	v25.2d,v11.2d,v21.2d	\n\t"\
		"fadd	v6.2d,v4.2d,v16.2d		\n\t	fmla	v24.2d,v11.2d,v20.2d	\n\t"/*14b+15a*/\
		"fadd	v7.2d,v5.2d,v17.2d		\n\t	fmls	v25.2d,v10.2d,v20.2d	\n\t"/*15b-14a*/\
		"fsub	v4.2d,v4.2d,v16.2d		\n\t	fadd	v14.2d,v12.2d,v24.2d	\n\t"\
		"fsub	v5.2d,v5.2d,v17.2d		\n\t	fadd	v15.2d,v13.2d,v25.2d	\n\t"\
		"ldp	q0,q1,[x14       ]		\n\t	fsub	v12.2d,v12.2d,v24.2d	\n\t"\
		"ldp	q2,q3,[x14,#0x100]		\n\t	fsub	v13.2d,v13.2d,v25.2d	\n\t"\
		/* cc0 := [c,s]: */\
		"ldp	q20,q21,[x10]			\n\t"\
		"fmul	v16.2d,v2.2d,v21.2d		\n\t	ldp	q8 ,q9 ,[x14,#0x080]		\n\t"\
		"fmul	v17.2d,v3.2d,v21.2d		\n\t	ldp	q10,q11,[x14,#0x180]		\n\t"\
/*2s+3c*/"fmla	v16.2d,v3.2d,v20.2d		\n\t	fmul	v18.2d,v10.2d,v20.2d	\n\t"\
/*3s-2c*/"fmls	v17.2d,v2.2d,v20.2d		\n\t	fmul	v19.2d,v11.2d,v20.2d	\n\t"\
		"fadd	v2.2d,v0.2d,v16.2d		\n\t	fmls	v18.2d,v11.2d,v21.2d	\n\t"/*10c-11s*/\
		"fadd	v3.2d,v1.2d,v17.2d		\n\t	fmla	v19.2d,v10.2d,v21.2d	\n\t"/*11c+10s*/\
		"fsub	v0.2d,v0.2d,v16.2d		\n\t	fsub	v10.2d,v8.2d,v18.2d		\n\t"\
		"fsub	v1.2d,v1.2d,v17.2d		\n\t	fsub	v11.2d,v9.2d,v19.2d		\n\t"\
/* 0 */	"fsub	v16.2d,v0.2d,v7.2d		\n\t	fadd	v8.2d ,v8.2d,v18.2d		\n\t"\
/* 1 */	"fsub	v17.2d,v1.2d,v6.2d		\n\t	fadd	v9.2d ,v9.2d,v19.2d		\n\t"\
/* 2 */	"fsub	v18.2d,v2.2d,v4.2d		\n\t	fsub	v20.2d,v8.2d ,v15.2d	\n\t"/*  8 */\
/* 3 */	"fsub	v19.2d,v3.2d,v5.2d		\n\t	fsub	v21.2d,v9.2d ,v14.2d	\n\t"/*  9 */\
		"fadd	v7.2d,v0.2d,v7.2d		\n\t	fsub	v22.2d,v10.2d,v12.2d	\n\t"/* 10 */\
		"fadd	v6.2d,v1.2d,v6.2d		\n\t	fsub	v23.2d,v11.2d,v13.2d	\n\t"/* 11 */\
		"fadd	v4.2d,v2.2d,v4.2d		\n\t	fadd	v15.2d,v8.2d ,v15.2d	\n\t"\
		"fadd	v5.2d,v3.2d,v5.2d		\n\t	fadd	v14.2d,v9.2d ,v14.2d	\n\t"\
		/* LDP byte-offsets must be in [-1024,1008], so incr x10 by 0x360 for final set of twiddles: */\
		"add	x10,x10,#0x360	\n\t"\
		/* Let c03 = [c,s], c0B = [a,b]: */"	fadd	v12.2d,v10.2d,v12.2d	\n\t"\
		"ldp	q0,q1,[x10      ]		\n\t	fadd	v13.2d,v11.2d,v13.2d	\n\t"\
		"ldp	q2,q3,[x10,#0x20]		\n\t"	/* Let c05 = [c,s], c0F = [a,b]: */\
		"fmul	v24.2d,v4.2d,v0.2d		\n\t	ldp	q8 ,q9 ,[x10,#0x80]		\n\t"\
		"fmul	v25.2d,v5.2d,v0.2d		\n\t	ldp	q10,q11,[x10,#0xa0]		\n\t"\
/*4c+5s*/"fmla	v24.2d,v5.2d,v1.2d		\n\t	fmul	v26.2d,v12.2d,v8.2d		\n\t"\
/*5c-4s*/"fmls	v25.2d,v4.2d,v1.2d		\n\t	fmul	v27.2d,v13.2d,v8.2d		\n\t"\
		"fmul	v4.2d,v17.2d,v2.2d		\n\t	fmla	v26.2d,v13.2d,v9.2d		\n\t"/*12c+13s*/\
		"fmul	v5.2d,v7.2d ,v2.2d		\n\t	fmls	v27.2d,v12.2d,v9.2d		\n\t"/*13c-12s*/\
/*1a-7b*/"fmls	v4.2d,v7.2d ,v3.2d		\n\t	fmul	v12.2d,v21.2d,v10.2d	\n\t"\
/*7a+1b*/"fmla	v5.2d,v17.2d,v3.2d		\n\t	fmul	v13.2d,v15.2d,v10.2d	\n\t"\
		"stp	q24,q25,[x0]			\n\t	fmls	v12.2d,v15.2d,v11.2d	\n\t"/* 9a-15b*/\
		"stp	q5 ,q4 ,[x1]			\n\t	fmla	v13.2d,v21.2d,v11.2d	\n\t"/*15a+ 9b*/\
		/* Let c13 = [c,s], c1B = [a,b]: */		/* Let c17 = [c,s], c1F = [a,b]: */\
		"ldp	q0,q1,[x10,#0x40]		\n\t"\
		"ldp	q2,q3,[x10,#0x60]		\n\t"\
		"fmul	v24.2d,v18.2d,v0.2d		\n\t	stp	q26,q27,[x4]			\n\t"\
		"fmul	v25.2d,v19.2d,v0.2d		\n\t	stp	q13,q12,[x5]			\n\t"\
/*2c+3s*/"fmla	v24.2d,v19.2d,v1.2d		\n\t	ldp	q8 ,q9 ,[x10,#0xc0]		\n\t"\
/*3c-2s*/"fmls	v25.2d,v18.2d,v1.2d		\n\t	ldp	q10,q11,[x10,#0xe0]		\n\t"\
		"fmul	v4.2d,v6.2d ,v2.2d		\n\t	fmul	v26.2d,v22.2d,v8.2d		\n\t"\
		"fmul	v5.2d,v16.2d,v2.2d		\n\t	fmul	v27.2d,v23.2d,v8.2d		\n\t"\
/*6a-0b*/"fmls	v4.2d,v16.2d,v3.2d		\n\t	fmla	v26.2d,v23.2d,v9.2d		\n\t"/*10c+11s*/\
/*0a+6b*/"fmla	v5.2d,v6.2d ,v3.2d		\n\t	fmls	v27.2d,v22.2d,v9.2d		\n\t"/*11c-10s*/\
		"stp	q24,q25,[x2]			\n\t	fmul	v12.2d,v14.2d,v10.2d	\n\t"\
		"stp	q5 ,q4 ,[x3]			\n\t	fmul	v13.2d,v20.2d,v10.2d	\n\t"\
		"										fmls	v12.2d,v20.2d,v11.2d	\n\t"/*14a-8b*/\
		"										fmla	v13.2d,v14.2d,v11.2d	\n\t"/*8a+14b*/\
		"										stp	q26,q27,[x6]			\n\t"\
		"										stp	q13,q12,[x7]			\n\t"\
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
		 ,[__r00] "m" (Xr00)\
		 ,[__isrt2] "m" (Xisrt2)\
		: "cc","memory","x0","x1","x2","x3","x4","x5","x6","x7","x8","x9","x10","x11","x12","x13","x14","x15",\
		"v0","v1","v2","v3","v4","v5","v6","v7","v8","v9","v10","v11","v12","v13","v14","v15",\
		"v16","v17","v18","v19","v20","v21","v22","v23","v24","v25","v26","v27", "v29"	/* Clobbered registers */\
	);\
	}

#elif defined(USE_AVX512)	// See AVX2 commentary for mem-layout and algo details...
	// Default here is to adapt the ALL_FMA versions of the AVX2 macros for AVX512 use,
	// and use the extra vector registers to store e.g. 1.0.

	// Cost [vector-ops only]: 386 MEM, 102 MUL, 368 FMA, 70 ADD
	#define SSE2_RADIX32_DIF_TWIDDLE(Xadd0,Xp01,Xp02,Xp03,Xp04,Xp08,Xp0C,Xp10,Xp18,Xr00)\
	{\
	__asm__ volatile (\
		/*...Block 1: */\
		"movq	%[__add0],%%rax							\n\t"\
		"movslq	%[__p08],%%rbx							\n\t"\
		"movslq	%[__p10],%%rcx							\n\t	movslq	%[__p04],%%r9		\n\t"\
		"movslq	%[__p18],%%rdx							\n\t	movq	%[__r00],%%rsi		\n\t"\
		"movq	%%rsi,%%r8								\n\t	leaq	0x2200(%%rsi),%%r8	\n\t"/* &two */\
		"leaq	(%%rax,%%rbx,8),%%rbx					\n\t vmovaps 0x1000(%%rsi), %%zmm29	\n\t"/* isrt2 */\
		"leaq	(%%rax,%%rcx,8),%%rcx					\n\t	vmovaps	-0x40(%%r8),%%zmm30	\n\t"/* 1.0 */\
		"leaq	(%%rax,%%rdx,8),%%rdx					\n\t	vmovaps	     (%%r8),%%zmm31	\n\t"/* 2.0 */\
		"														vmovaps	 0x40(%%r8),%%zmm28	\n\t"/* sqrt2 */\
		/* SSE2_RADIX4_DIF_4TWIDDLE_A(r00,c10) */				/* SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO_B(r08,c04) */\
		"vmovaps	     (%%rcx),%%zmm4					\n\t	vmovaps	     (%%rax,%%r9,8),%%zmm8 	\n\t"\
		"vmovaps	0x040(%%rcx),%%zmm5					\n\t	vmovaps	0x040(%%rax,%%r9,8),%%zmm9 	\n\t"\
		"vmovaps	0x1240(%%rsi),%%zmm2	/* c10 */	\n\t	vmovaps	0x13c0(%%rsi),%%zmm14	/* c04 */	\n\t"\
		"vmovaps	0x1280(%%rsi),%%zmm3				\n\t	vmovaps	0x1400(%%rsi),%%zmm15	\n\t"\
		"vmovaps	%%zmm4,%%zmm6						\n\t	vmovaps	%%zmm8 ,%%zmm10	\n\t"\
		"vmovaps	%%zmm5,%%zmm7						\n\t	vmovaps	%%zmm9 ,%%zmm11	\n\t"\
		"vmulpd	%%zmm2,%%zmm4,%%zmm4					\n\t	vmulpd	%%zmm14,%%zmm8 ,%%zmm8 	\n\t"\
		"vmulpd	%%zmm2,%%zmm5,%%zmm5					\n\t	vmulpd	%%zmm14,%%zmm9 ,%%zmm9 	\n\t"\
		"vmovaps	     (%%rax),%%zmm0					\n\t	vmovaps	     (%%rcx,%%r9,8),%%zmm12	\n\t"\
		"vmovaps	0x040(%%rax),%%zmm1					\n\t	vmovaps	0x040(%%rcx,%%r9,8),%%zmm13	\n\t"\
		"vfnmadd231pd	%%zmm3,%%zmm7,%%zmm4			\n\t	vfnmadd231pd	%%zmm15,%%zmm11,%%zmm8 	\n\t"\
		" vfmadd231pd	%%zmm3,%%zmm6,%%zmm5			\n\t	 vfmadd231pd	%%zmm15,%%zmm10,%%zmm9 	\n\t"\
		"vmovaps	0x1300(%%rsi),%%zmm16		\n\t	vmovaps		0x1440(%%rsi),%%zmm20		\n\t"\
		"vmovaps	0x1340(%%rsi),%%zmm17		\n\t	vmovaps		0x1480(%%rsi),%%zmm21		\n\t"\
		"vmovaps	0x1380(%%rsi),%%zmm18		\n\t	vmovaps		0x1540(%%rsi),%%zmm22		\n\t"\
		"vmovaps	0x12c0(%%rsi),%%zmm19		\n\t	vmovaps		0x1580(%%rsi),%%zmm23		\n\t"\
		"vmovaps	%%zmm0,%%zmm2						\n\t	vmovaps	%%zmm12,%%zmm14	\n\t"\
		"vmovaps	%%zmm1,%%zmm3						\n\t	vmovaps	%%zmm13,%%zmm15	\n\t"\
	"vaddpd	%%zmm4,%%zmm0,%%zmm0				\n\t		vmulpd		%%zmm20,%%zmm12,%%zmm12	\n\t"/* c14 */\
	"vaddpd	%%zmm5,%%zmm1,%%zmm1				\n\t		vmulpd		%%zmm20,%%zmm13,%%zmm13	\n\t"\
	"vsubpd	%%zmm4,%%zmm2,%%zmm2				\n\t	vfnmadd231pd	%%zmm21,%%zmm15,%%zmm12	\n\t"\
	"vsubpd	%%zmm5,%%zmm3,%%zmm3				\n\t	 vfmadd231pd	%%zmm21,%%zmm14,%%zmm13	\n\t"\
		"vmovaps	     (%%rdx),%%zmm4					\n\t	vmovaps	%%zmm8 ,%%zmm10	\n\t"\
		"vmovaps	0x040(%%rdx),%%zmm5					\n\t	vmovaps	%%zmm9 ,%%zmm11	\n\t"\
		"vmovaps	     (%%rdx),%%zmm6					\n\t	vaddpd	%%zmm12,%%zmm8 ,%%zmm8 	\n\t"\
		"vmovaps	0x040(%%rdx),%%zmm7					\n\t	vaddpd	%%zmm13,%%zmm9 ,%%zmm9 	\n\t"\
/* c18 */"vmulpd	%%zmm17,%%zmm4,%%zmm4				\n\t	vsubpd	%%zmm12,%%zmm10,%%zmm10	\n\t"\
		"vmulpd		%%zmm17,%%zmm5,%%zmm5				\n\t	vsubpd	%%zmm13,%%zmm11,%%zmm11	\n\t"\
	"vfnmadd231pd	%%zmm18,%%zmm7,%%zmm4				\n\t	vmovaps	     (%%rdx,%%r9,8),%%zmm12	\n\t"\
	" vfmadd231pd	%%zmm18,%%zmm6,%%zmm5				\n\t	vmovaps	0x040(%%rdx,%%r9,8),%%zmm13	\n\t"\
		"vmovaps	%%zmm5,0x040(%%rsi)					\n\t	vmovaps	     (%%rdx,%%r9,8),%%zmm14	\n\t"\
		"vmovaps	%%zmm4,     (%%rsi)	/* tmpstr r00 */\n\t	vmovaps	0x040(%%rdx,%%r9,8),%%zmm15	\n\t"\
		"vmovaps	     (%%rbx),%%zmm4					\n\t		vmulpd		%%zmm22,%%zmm12,%%zmm12	\n\t"/* c1C */\
		"vmovaps	0x040(%%rbx),%%zmm5					\n\t		vmulpd		%%zmm22,%%zmm13,%%zmm13	\n\t"\
		"vmovaps	     (%%rbx),%%zmm6					\n\t	vfnmadd231pd	%%zmm23,%%zmm15,%%zmm12	\n\t"\
		"vmovaps	0x040(%%rbx),%%zmm7					\n\t	 vfmadd231pd	%%zmm23,%%zmm14,%%zmm13	\n\t"\
/* c08 */"vmulpd	%%zmm19,%%zmm4,%%zmm4				\n\t	vmovaps	%%zmm13,%%zmm23	\n\t"\
		"vmulpd		%%zmm19,%%zmm5,%%zmm5				\n\t	vmovaps	%%zmm12,%%zmm22	/* tmpstr r08 */\n\t"\
	"vfnmadd231pd	%%zmm16,%%zmm7,%%zmm4				\n\t	vmovaps	     (%%rbx,%%r9,8),%%zmm12	\n\t"\
	" vfmadd231pd	%%zmm16,%%zmm6,%%zmm5				\n\t	vmovaps	0x040(%%rbx,%%r9,8),%%zmm13	\n\t"\
		"vmovaps	%%zmm5,%%zmm7						\n\t	vmovaps		0x14c0(%%rsi),%%zmm20		\n\t"\
		"vmovaps	%%zmm4,%%zmm6						\n\t	vmovaps		0x1500(%%rsi),%%zmm21		\n\t"\
		"vsubpd	     (%%rsi),%%zmm4,%%zmm4	/* r00 */	\n\t		vmulpd		%%zmm20,%%zmm12,%%zmm14	\n\t"/* c0C */\
		"vsubpd	0x040(%%rsi),%%zmm5,%%zmm5				\n\t		vmulpd		%%zmm20,%%zmm13,%%zmm15	\n\t"\
		"vaddpd	     (%%rsi),%%zmm6,%%zmm6				\n\t	vfnmadd231pd	%%zmm21,%%zmm13,%%zmm14	\n\t"\
		"vaddpd	0x040(%%rsi),%%zmm7,%%zmm7				\n\t	 vfmadd231pd	%%zmm21,%%zmm12,%%zmm15	\n\t"\
		"														vmovaps	%%zmm15,%%zmm13	\n\t"\
		"														vmovaps	%%zmm14,%%zmm12	\n\t"\
		"vsubpd	%%zmm4,%%zmm3,%%zmm3			\n\t	vsubpd	%%zmm22,%%zmm12,%%zmm12	\n\t"\
		"vsubpd	%%zmm5,%%zmm2,%%zmm2			\n\t	vsubpd	%%zmm23,%%zmm13,%%zmm13	\n\t"\
		"vsubpd	%%zmm6,%%zmm0,%%zmm0			\n\t	vaddpd	%%zmm22,%%zmm14,%%zmm14	\n\t"\
		"vsubpd	%%zmm7,%%zmm1,%%zmm1			\n\t	vaddpd	%%zmm23,%%zmm15,%%zmm15	\n\t"\
		"												\n\t	vsubpd	%%zmm14,%%zmm8 ,%%zmm8 	\n\t"\
		"												\n\t	vsubpd	%%zmm15,%%zmm9 ,%%zmm9 	\n\t"\
		"												\n\t	vsubpd	%%zmm13,%%zmm10,%%zmm10	\n\t"\
		"														vsubpd	%%zmm12,%%zmm11,%%zmm11	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm3,%%zmm4			\n\t	vfmadd132pd	%%zmm31,%%zmm8 ,%%zmm14	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm2,%%zmm5			\n\t	vfmadd132pd	%%zmm31,%%zmm10,%%zmm13	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm0,%%zmm6			\n\t	vfmadd132pd	%%zmm31,%%zmm9 ,%%zmm15	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm1,%%zmm7			\n\t	vfmadd132pd	%%zmm31,%%zmm11,%%zmm12	\n\t"\
		/* SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r00) */\
		"														vsubpd	%%zmm12,%%zmm10,%%zmm16	\n\t"\
		"														vsubpd	%%zmm11,%%zmm13,%%zmm17	\n\t"\
		"														vaddpd	%%zmm12,%%zmm10,%%zmm10	\n\t"\
		"														vaddpd	%%zmm11,%%zmm13,%%zmm13	\n\t"\
		"														vfnmadd231pd %%zmm29,%%zmm16,%%zmm2	\n\t"\
		"														vfnmadd231pd %%zmm29,%%zmm17,%%zmm3	\n\t"\
		"vsubpd	%%zmm14,%%zmm6,%%zmm6			\n\t	vfnmadd231pd %%zmm29,%%zmm10,%%zmm4	\n\t"\
		"vsubpd	%%zmm9 ,%%zmm0,%%zmm0			\n\t	vfnmadd231pd %%zmm29,%%zmm13,%%zmm5	\n\t"\
		"vsubpd	%%zmm8 ,%%zmm1 ,%%zmm1 					\n\t	vsubpd	%%zmm15,%%zmm7 ,%%zmm7 	\n\t"\
		"vmovaps	%%zmm6 ,0x200(%%rsi)				\n\t	vmovaps	%%zmm2 ,0x280(%%rsi)	\n\t"\
		"vmovaps	%%zmm0 ,0x100(%%rsi)				\n\t	vmovaps	%%zmm5 ,0x180(%%rsi)	\n\t"\
		"vmovaps	%%zmm7 ,0x240(%%rsi)				\n\t	vmovaps	%%zmm4 ,0x2c0(%%rsi)	\n\t"\
		"vmovaps	%%zmm1 ,0x340(%%rsi)				\n\t	vmovaps	%%zmm3 ,0x3c0(%%rsi)	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm6 ,%%zmm14			\n\t	vfmadd132pd	%%zmm28,%%zmm2 ,%%zmm16	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm0 ,%%zmm9 			\n\t	vfmadd132pd	%%zmm28,%%zmm5 ,%%zmm13	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm7 ,%%zmm15			\n\t	vfmadd132pd	%%zmm28,%%zmm4 ,%%zmm10	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm1 ,%%zmm8 			\n\t	vfmadd132pd	%%zmm28,%%zmm3 ,%%zmm17	\n\t"\
		"vmovaps	%%zmm14,     (%%rsi)				\n\t	vmovaps	%%zmm16,0x080(%%rsi)	\n\t"\
		"vmovaps	%%zmm9 ,0x300(%%rsi)				\n\t	vmovaps	%%zmm13,0x380(%%rsi)	\n\t"\
		"vmovaps	%%zmm15,0x040(%%rsi)				\n\t	vmovaps	%%zmm10,0x0c0(%%rsi)	\n\t"\
		"vmovaps	%%zmm8 ,0x140(%%rsi)				\n\t	vmovaps	%%zmm17,0x1c0(%%rsi)	\n\t"\
		"\n\t"\
	/***************************************/\
		"\n\t"\
		/*...Block 2: */\
		"movslq	%[__p02],%%rdi							\n\t	addq	%%rdi,%%r9	\n\t"\
		/* SSE2_RADIX4_DIF_4TWIDDLE_B(r10,c02) */				/* SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO_B(r18,c06) */\
		"addq	$0x400,%%rsi	/* r10 */				\n\t"\
		"vmovaps	     (%%rax,%%rdi,8),%%zmm0			\n\t	vmovaps	     (%%rax,%%r9,8),%%zmm8 	\n\t"\
		"vmovaps	     (%%rcx,%%rdi,8),%%zmm4			\n\t	vmovaps	     (%%rcx,%%r9,8),%%zmm12	\n\t"\
		"vmovaps	0x040(%%rax,%%rdi,8),%%zmm1			\n\t	vmovaps	0x040(%%rax,%%r9,8),%%zmm9 	\n\t"\
		"vmovaps	0x040(%%rcx,%%rdi,8),%%zmm5			\n\t	vmovaps	0x040(%%rcx,%%r9,8),%%zmm13	\n\t"\
		"vmovaps	0x11c0(%%rsi),%%zmm6				\n\t	vmovaps	0x13c0(%%rsi),%%zmm14	\n\t"\
		"vmovaps	0x1200(%%rsi),%%zmm7				\n\t	vmovaps	0x1400(%%rsi),%%zmm15	\n\t"\
		"vmovaps	%%zmm0,%%zmm2						\n\t	vmovaps	%%zmm8 ,%%zmm10	\n\t"\
		"vmovaps	%%zmm1,%%zmm3						\n\t	vmovaps	%%zmm9 ,%%zmm11	\n\t"\
		"vmovaps	0x1240(%%rsi),%%zmm16		\n\t	vmovaps		0x1440(%%rsi),%%zmm20		\n\t"\
		"vmovaps	0x1280(%%rsi),%%zmm17		\n\t	vmovaps		0x1480(%%rsi),%%zmm21		\n\t"\
		"vmovaps	0x1340(%%rsi),%%zmm18		\n\t	vmovaps		0x1540(%%rsi),%%zmm22		\n\t"\
		"vmovaps	0x1380(%%rsi),%%zmm19		\n\t	vmovaps		0x1580(%%rsi),%%zmm23		\n\t"\
/* c02 */"	vmulpd		%%zmm6,%%zmm0,%%zmm0		\n\t		vmulpd		%%zmm14,%%zmm8 ,%%zmm8 	\n\t"/* c06 */\
		"	vmulpd		%%zmm6,%%zmm1,%%zmm1		\n\t		vmulpd		%%zmm14,%%zmm9 ,%%zmm9 	\n\t"\
		"vfnmadd231pd	%%zmm7,%%zmm3,%%zmm0		\n\t	vfnmadd231pd	%%zmm15,%%zmm11,%%zmm8 	\n\t"\
		" vfmadd231pd	%%zmm7,%%zmm2,%%zmm1		\n\t	 vfmadd231pd	%%zmm15,%%zmm10,%%zmm9 	\n\t"\
		"vmovaps	%%zmm4,%%zmm6						\n\t	vmovaps	%%zmm12,%%zmm14	\n\t"\
		"vmovaps	%%zmm5,%%zmm7						\n\t	vmovaps	%%zmm13,%%zmm15	\n\t"\
/* c12 */"	vmulpd		%%zmm16,%%zmm4,%%zmm4		\n\t		vmulpd		%%zmm20,%%zmm12,%%zmm12	\n\t"/* c16 */\
		"	vmulpd		%%zmm16,%%zmm5,%%zmm5		\n\t		vmulpd		%%zmm20,%%zmm13,%%zmm13	\n\t"\
		"vfnmadd231pd	%%zmm17,%%zmm7,%%zmm4		\n\t	vfnmadd231pd	%%zmm21,%%zmm15,%%zmm12	\n\t"\
		" vfmadd231pd	%%zmm17,%%zmm6,%%zmm5		\n\t	 vfmadd231pd	%%zmm21,%%zmm14,%%zmm13	\n\t"\
		"vmovaps	%%zmm0,%%zmm2						\n\t	vmovaps	%%zmm8 ,%%zmm10	\n\t"\
		"vmovaps	%%zmm1,%%zmm3						\n\t	vmovaps	%%zmm9 ,%%zmm11	\n\t"\
		"vaddpd	%%zmm4,%%zmm0,%%zmm0			\n\t	vaddpd	%%zmm12,%%zmm8 ,%%zmm8 	\n\t"\
		"vaddpd	%%zmm5,%%zmm1,%%zmm1			\n\t	vaddpd	%%zmm13,%%zmm9 ,%%zmm9 	\n\t"\
		"vsubpd	%%zmm4,%%zmm2,%%zmm2			\n\t	vsubpd	%%zmm12,%%zmm10,%%zmm10	\n\t"\
		"vsubpd	%%zmm5,%%zmm3,%%zmm3			\n\t	vsubpd	%%zmm13,%%zmm11,%%zmm11	\n\t"\
		"vmovaps	     (%%rdx,%%rdi,8),%%zmm6			\n\t	vmovaps	     (%%rdx,%%r9,8),%%zmm14	\n\t"\
		"vmovaps	0x040(%%rdx,%%rdi,8),%%zmm7			\n\t	vmovaps	0x040(%%rdx,%%r9,8),%%zmm15	\n\t"\
/* c1A */"	vmulpd		%%zmm18,%%zmm6,%%zmm4		\n\t		vmulpd		%%zmm22,%%zmm14,%%zmm12	\n\t"/* c1E */\
		"	vmulpd		%%zmm18,%%zmm7,%%zmm5		\n\t		vmulpd		%%zmm22,%%zmm15,%%zmm13	\n\t"\
		"vfnmadd231pd	%%zmm19,%%zmm7,%%zmm4		\n\t	vfnmadd231pd	%%zmm23,%%zmm15,%%zmm12	\n\t"\
		" vfmadd231pd	%%zmm19,%%zmm6,%%zmm5		\n\t	 vfmadd231pd	%%zmm23,%%zmm14,%%zmm13	\n\t"\
		"vmovaps	0x12c0(%%rsi),%%zmm16		\n\t	vmovaps		0x14c0(%%rsi),%%zmm20		\n\t"\
		"vmovaps	0x1300(%%rsi),%%zmm17		\n\t	vmovaps		0x1500(%%rsi),%%zmm21		\n\t"\
		"vmovaps	%%zmm5,0x040(%%rsi)					\n\t	vmovaps	%%zmm13,0x240(%%rsi)	\n\t"\
		"vmovaps	%%zmm4,     (%%rsi)					\n\t	vmovaps	%%zmm12,0x200(%%rsi)	/* tmpstr r08 */\n\t"\
		"vmovaps	     (%%rbx,%%rdi,8),%%zmm6			\n\t	vmovaps	     (%%rbx,%%r9,8),%%zmm14	\n\t"\
		"vmovaps	0x040(%%rbx,%%rdi,8),%%zmm7			\n\t	vmovaps	0x040(%%rbx,%%r9,8),%%zmm15	\n\t"\
/* c0A */"	vmulpd		%%zmm16,%%zmm6,%%zmm4		\n\t		vmulpd		%%zmm20,%%zmm14,%%zmm12	\n\t"/* c0E */\
		"	vmulpd		%%zmm16,%%zmm7,%%zmm5		\n\t		vmulpd		%%zmm20,%%zmm15,%%zmm13	\n\t"\
		"vfnmadd231pd	%%zmm17,%%zmm7,%%zmm4		\n\t	vfnmadd231pd	%%zmm21,%%zmm15,%%zmm12	\n\t"\
		" vfmadd231pd	%%zmm17,%%zmm6,%%zmm5		\n\t	 vfmadd231pd	%%zmm21,%%zmm14,%%zmm13	\n\t"\
		"vmovaps	     (%%rsi),%%zmm16			\n\t	vmovaps	0x200(%%rsi),%%zmm18	\n\t"\
		"vmovaps	0x040(%%rsi),%%zmm17			\n\t	vmovaps	0x240(%%rsi),%%zmm19	\n\t"\
		"vaddpd		%%zmm16,%%zmm4,%%zmm6				\n\t	vaddpd		%%zmm18,%%zmm12,%%zmm14	\n\t"\
		"vaddpd		%%zmm17,%%zmm5,%%zmm7				\n\t	vaddpd		%%zmm19,%%zmm13,%%zmm15	\n\t"\
		"vsubpd		%%zmm16,%%zmm4,%%zmm4				\n\t	vsubpd		%%zmm18,%%zmm12,%%zmm12	\n\t"\
		"vsubpd		%%zmm17,%%zmm5,%%zmm5				\n\t	vsubpd		%%zmm19,%%zmm13,%%zmm13	\n\t"\
		"vsubpd	%%zmm6,%%zmm0,%%zmm0			\n\t	vsubpd	    %%zmm14,%%zmm8 ,%%zmm8 	\n\t"\
		"vsubpd	%%zmm5,%%zmm2,%%zmm2			\n\t	vsubpd	    %%zmm15,%%zmm9 ,%%zmm9 	\n\t"\
		"vsubpd	%%zmm7,%%zmm1,%%zmm1			\n\t	vsubpd	%%zmm13,%%zmm10,%%zmm10	\n\t"\
		"vsubpd	%%zmm4,%%zmm3,%%zmm3			\n\t	vsubpd	%%zmm12,%%zmm11,%%zmm11	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm0,%%zmm6			\n\t	vfmadd132pd	%%zmm31,%%zmm8 ,%%zmm14	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm2,%%zmm5			\n\t	vfmadd132pd	%%zmm31,%%zmm10,%%zmm13	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm1,%%zmm7			\n\t	vfmadd132pd	%%zmm31,%%zmm9 ,%%zmm15	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm3,%%zmm4			\n\t	vfmadd132pd	%%zmm31,%%zmm11,%%zmm12	\n\t"\
		/* SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r10) */\
		"														vsubpd	%%zmm12,%%zmm10,%%zmm16	\n\t"\
		"														vsubpd	%%zmm11,%%zmm13,%%zmm17	\n\t"\
		"														vaddpd	%%zmm12,%%zmm10,%%zmm10	\n\t"\
		"														vaddpd	%%zmm11,%%zmm13,%%zmm13	\n\t"\
		"														vfnmadd231pd %%zmm29,%%zmm16,%%zmm2	\n\t"\
		"														vfnmadd231pd %%zmm29,%%zmm17,%%zmm3	\n\t"\
		"vsubpd	%%zmm14,%%zmm6,%%zmm6			\n\t	vfnmadd231pd %%zmm29,%%zmm10,%%zmm4	\n\t"\
		"vsubpd	%%zmm9 ,%%zmm0,%%zmm0			\n\t	vfnmadd231pd %%zmm29,%%zmm13,%%zmm5	\n\t"\
		"vsubpd	%%zmm8 ,%%zmm1 ,%%zmm1 					\n\t	vsubpd	%%zmm15,%%zmm7 ,%%zmm7 	\n\t"\
		"vmovaps	%%zmm6 ,0x200(%%rsi)				\n\t	vmovaps	%%zmm2 ,0x280(%%rsi)	\n\t"\
		"vmovaps	%%zmm0 ,0x100(%%rsi)				\n\t	vmovaps	%%zmm5 ,0x180(%%rsi)	\n\t"\
		"vmovaps	%%zmm7 ,0x240(%%rsi)				\n\t	vmovaps	%%zmm4 ,0x2c0(%%rsi)	\n\t"\
		"vmovaps	%%zmm1 ,0x340(%%rsi)				\n\t	vmovaps	%%zmm3 ,0x3c0(%%rsi)	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm6 ,%%zmm14			\n\t	vfmadd132pd	%%zmm28,%%zmm2 ,%%zmm16	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm0 ,%%zmm9 			\n\t	vfmadd132pd	%%zmm28,%%zmm5 ,%%zmm13	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm7 ,%%zmm15			\n\t	vfmadd132pd	%%zmm28,%%zmm4 ,%%zmm10	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm1 ,%%zmm8 			\n\t	vfmadd132pd	%%zmm28,%%zmm3 ,%%zmm17	\n\t"\
		"vmovaps	%%zmm14,     (%%rsi)				\n\t	vmovaps	%%zmm16,0x080(%%rsi)	\n\t"\
		"vmovaps	%%zmm9 ,0x300(%%rsi)				\n\t	vmovaps	%%zmm13,0x380(%%rsi)	\n\t"\
		"vmovaps	%%zmm15,0x040(%%rsi)				\n\t	vmovaps	%%zmm10,0x0c0(%%rsi)	\n\t"\
		"vmovaps	%%zmm8 ,0x140(%%rsi)				\n\t	vmovaps	%%zmm17,0x1c0(%%rsi)	\n\t"\
		"\n\t"\
	/***************************************/\
		"\n\t"\
		/*...Block 3: */\
														"		subq	%%rdi,%%r9	\n\t"\
		"movslq	%[__p01],%%rdi							\n\t	addq	%%rdi,%%r9	\n\t"\
		/* SSE2_RADIX4_DIF_4TWIDDLE_B(r20,c01) */				/* SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO_B(r28,c05) */\
		"addq	$0x400,%%rsi	/* r20 */				\n\t"\
		"vmovaps	     (%%rax,%%rdi,8),%%zmm0			\n\t	vmovaps	     (%%rax,%%r9,8),%%zmm8 	\n\t"\
		"vmovaps	     (%%rcx,%%rdi,8),%%zmm4			\n\t	vmovaps	     (%%rcx,%%r9,8),%%zmm12	\n\t"\
		"vmovaps	0x040(%%rax,%%rdi,8),%%zmm1			\n\t	vmovaps	0x040(%%rax,%%r9,8),%%zmm9 	\n\t"\
		"vmovaps	0x040(%%rcx,%%rdi,8),%%zmm5			\n\t	vmovaps	0x040(%%rcx,%%r9,8),%%zmm13	\n\t"\
		"vmovaps	0x11c0(%%rsi),%%zmm6				\n\t	vmovaps	0x13c0(%%rsi),%%zmm14	\n\t"\
		"vmovaps	0x1200(%%rsi),%%zmm7				\n\t	vmovaps	0x1400(%%rsi),%%zmm15	\n\t"\
		"vmovaps	%%zmm0,%%zmm2						\n\t	vmovaps	%%zmm8 ,%%zmm10	\n\t"\
		"vmovaps	%%zmm1,%%zmm3						\n\t	vmovaps	%%zmm9 ,%%zmm11	\n\t"\
		"vmovaps	0x1240(%%rsi),%%zmm16		\n\t	vmovaps		0x1440(%%rsi),%%zmm20		\n\t"\
		"vmovaps	0x1280(%%rsi),%%zmm17		\n\t	vmovaps		0x1480(%%rsi),%%zmm21		\n\t"\
		"vmovaps	0x1340(%%rsi),%%zmm18		\n\t	vmovaps		0x1540(%%rsi),%%zmm22		\n\t"\
		"vmovaps	0x1380(%%rsi),%%zmm19		\n\t	vmovaps		0x1580(%%rsi),%%zmm23		\n\t"\
/* c01 */"	vmulpd		%%zmm6,%%zmm0,%%zmm0		\n\t		vmulpd		%%zmm14,%%zmm8 ,%%zmm8 	\n\t"/* c05 */\
		"	vmulpd		%%zmm6,%%zmm1,%%zmm1		\n\t		vmulpd		%%zmm14,%%zmm9 ,%%zmm9 	\n\t"\
		"vfnmadd231pd	%%zmm7,%%zmm3,%%zmm0		\n\t	vfnmadd231pd	%%zmm15,%%zmm11,%%zmm8 	\n\t"\
		" vfmadd231pd	%%zmm7,%%zmm2,%%zmm1		\n\t	 vfmadd231pd	%%zmm15,%%zmm10,%%zmm9 	\n\t"\
		"vmovaps	%%zmm4,%%zmm6						\n\t	vmovaps	%%zmm12,%%zmm14	\n\t"\
		"vmovaps	%%zmm5,%%zmm7						\n\t	vmovaps	%%zmm13,%%zmm15	\n\t"\
/* c11 */"	vmulpd		%%zmm16,%%zmm4,%%zmm4		\n\t		vmulpd		%%zmm20,%%zmm12,%%zmm12	\n\t"/* c15 */\
		"	vmulpd		%%zmm16,%%zmm5,%%zmm5		\n\t		vmulpd		%%zmm20,%%zmm13,%%zmm13	\n\t"\
		"vfnmadd231pd	%%zmm17,%%zmm7,%%zmm4		\n\t	vfnmadd231pd	%%zmm21,%%zmm15,%%zmm12	\n\t"\
		" vfmadd231pd	%%zmm17,%%zmm6,%%zmm5		\n\t	 vfmadd231pd	%%zmm21,%%zmm14,%%zmm13	\n\t"\
		"vmovaps	%%zmm0,%%zmm2						\n\t	vmovaps	%%zmm8 ,%%zmm10	\n\t"\
		"vmovaps	%%zmm1,%%zmm3						\n\t	vmovaps	%%zmm9 ,%%zmm11	\n\t"\
		"vaddpd	%%zmm4,%%zmm0,%%zmm0			\n\t	vaddpd	%%zmm12,%%zmm8 ,%%zmm8 	\n\t"\
		"vaddpd	%%zmm5,%%zmm1,%%zmm1			\n\t	vaddpd	%%zmm13,%%zmm9 ,%%zmm9 	\n\t"\
		"vsubpd	%%zmm4,%%zmm2,%%zmm2			\n\t	vsubpd	%%zmm12,%%zmm10,%%zmm10	\n\t"\
		"vsubpd	%%zmm5,%%zmm3,%%zmm3			\n\t	vsubpd	%%zmm13,%%zmm11,%%zmm11	\n\t"\
		"vmovaps	     (%%rdx,%%rdi,8),%%zmm6			\n\t	vmovaps	     (%%rdx,%%r9,8),%%zmm14	\n\t"\
		"vmovaps	0x040(%%rdx,%%rdi,8),%%zmm7			\n\t	vmovaps	0x040(%%rdx,%%r9,8),%%zmm15	\n\t"\
/* c19 */"	vmulpd		%%zmm18,%%zmm6,%%zmm4		\n\t		vmulpd		%%zmm22,%%zmm14,%%zmm12	\n\t"/* c1D */\
		"	vmulpd		%%zmm18,%%zmm7,%%zmm5		\n\t		vmulpd		%%zmm22,%%zmm15,%%zmm13	\n\t"\
		"vfnmadd231pd	%%zmm19,%%zmm7,%%zmm4		\n\t	vfnmadd231pd	%%zmm23,%%zmm15,%%zmm12	\n\t"\
		" vfmadd231pd	%%zmm19,%%zmm6,%%zmm5		\n\t	 vfmadd231pd	%%zmm23,%%zmm14,%%zmm13	\n\t"\
		"vmovaps	0x12c0(%%rsi),%%zmm16		\n\t	vmovaps		0x14c0(%%rsi),%%zmm20		\n\t"\
		"vmovaps	0x1300(%%rsi),%%zmm17		\n\t	vmovaps		0x1500(%%rsi),%%zmm21		\n\t"\
		"vmovaps	%%zmm5,0x040(%%rsi)					\n\t	vmovaps	%%zmm13,0x240(%%rsi)	\n\t"\
		"vmovaps	%%zmm4,     (%%rsi)					\n\t	vmovaps	%%zmm12,0x200(%%rsi)	/* tmpstr r08 */\n\t"\
		"vmovaps	     (%%rbx,%%rdi,8),%%zmm6			\n\t	vmovaps	     (%%rbx,%%r9,8),%%zmm14	\n\t"\
		"vmovaps	0x040(%%rbx,%%rdi,8),%%zmm7			\n\t	vmovaps	0x040(%%rbx,%%r9,8),%%zmm15	\n\t"\
/* c09 */"	vmulpd		%%zmm16,%%zmm6,%%zmm4		\n\t		vmulpd		%%zmm20,%%zmm14,%%zmm12	\n\t"/* c0D */\
		"	vmulpd		%%zmm16,%%zmm7,%%zmm5		\n\t		vmulpd		%%zmm20,%%zmm15,%%zmm13	\n\t"\
		"vfnmadd231pd	%%zmm17,%%zmm7,%%zmm4		\n\t	vfnmadd231pd	%%zmm21,%%zmm15,%%zmm12	\n\t"\
		" vfmadd231pd	%%zmm17,%%zmm6,%%zmm5		\n\t	 vfmadd231pd	%%zmm21,%%zmm14,%%zmm13	\n\t"\
		"vmovaps	     (%%rsi),%%zmm16			\n\t	vmovaps	0x200(%%rsi),%%zmm18	\n\t"\
		"vmovaps	0x040(%%rsi),%%zmm17			\n\t	vmovaps	0x240(%%rsi),%%zmm19	\n\t"\
		"vaddpd		%%zmm16,%%zmm4,%%zmm6				\n\t	vaddpd		%%zmm18,%%zmm12,%%zmm14	\n\t"\
		"vaddpd		%%zmm17,%%zmm5,%%zmm7				\n\t	vaddpd		%%zmm19,%%zmm13,%%zmm15	\n\t"\
		"vsubpd		%%zmm16,%%zmm4,%%zmm4				\n\t	vsubpd		%%zmm18,%%zmm12,%%zmm12	\n\t"\
		"vsubpd		%%zmm17,%%zmm5,%%zmm5				\n\t	vsubpd		%%zmm19,%%zmm13,%%zmm13	\n\t"\
		"vsubpd	%%zmm6,%%zmm0,%%zmm0			\n\t	vsubpd	    %%zmm14,%%zmm8 ,%%zmm8 	\n\t"\
		"vsubpd	%%zmm5,%%zmm2,%%zmm2			\n\t	vsubpd	    %%zmm15,%%zmm9 ,%%zmm9 	\n\t"\
		"vsubpd	%%zmm7,%%zmm1,%%zmm1			\n\t	vsubpd	%%zmm13,%%zmm10,%%zmm10	\n\t"\
		"vsubpd	%%zmm4,%%zmm3,%%zmm3			\n\t	vsubpd	%%zmm12,%%zmm11,%%zmm11	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm0,%%zmm6			\n\t	vfmadd132pd	%%zmm31,%%zmm8 ,%%zmm14	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm2,%%zmm5			\n\t	vfmadd132pd	%%zmm31,%%zmm10,%%zmm13	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm1,%%zmm7			\n\t	vfmadd132pd	%%zmm31,%%zmm9 ,%%zmm15	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm3,%%zmm4			\n\t	vfmadd132pd	%%zmm31,%%zmm11,%%zmm12	\n\t"\
		/* SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r20) */\
		"														vsubpd	%%zmm12,%%zmm10,%%zmm16	\n\t"\
		"														vsubpd	%%zmm11,%%zmm13,%%zmm17	\n\t"\
		"														vaddpd	%%zmm12,%%zmm10,%%zmm10	\n\t"\
		"														vaddpd	%%zmm11,%%zmm13,%%zmm13	\n\t"\
		"														vfnmadd231pd %%zmm29,%%zmm16,%%zmm2	\n\t"\
		"														vfnmadd231pd %%zmm29,%%zmm17,%%zmm3	\n\t"\
		"vsubpd	%%zmm14,%%zmm6,%%zmm6			\n\t	vfnmadd231pd %%zmm29,%%zmm10,%%zmm4	\n\t"\
		"vsubpd	%%zmm9 ,%%zmm0,%%zmm0			\n\t	vfnmadd231pd %%zmm29,%%zmm13,%%zmm5	\n\t"\
		"vsubpd	%%zmm8 ,%%zmm1 ,%%zmm1 					\n\t	vsubpd	%%zmm15,%%zmm7 ,%%zmm7 	\n\t"\
		"vmovaps	%%zmm6 ,0x200(%%rsi)				\n\t	vmovaps	%%zmm2 ,0x280(%%rsi)	\n\t"\
		"vmovaps	%%zmm0 ,0x100(%%rsi)				\n\t	vmovaps	%%zmm5 ,0x180(%%rsi)	\n\t"\
		"vmovaps	%%zmm7 ,0x240(%%rsi)				\n\t	vmovaps	%%zmm4 ,0x2c0(%%rsi)	\n\t"\
		"vmovaps	%%zmm1 ,0x340(%%rsi)				\n\t	vmovaps	%%zmm3 ,0x3c0(%%rsi)	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm6 ,%%zmm14			\n\t	vfmadd132pd	%%zmm28,%%zmm2 ,%%zmm16	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm0 ,%%zmm9 			\n\t	vfmadd132pd	%%zmm28,%%zmm5 ,%%zmm13	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm7 ,%%zmm15			\n\t	vfmadd132pd	%%zmm28,%%zmm4 ,%%zmm10	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm1 ,%%zmm8 			\n\t	vfmadd132pd	%%zmm28,%%zmm3 ,%%zmm17	\n\t"\
		"vmovaps	%%zmm14,     (%%rsi)				\n\t	vmovaps	%%zmm16,0x080(%%rsi)	\n\t"\
		"vmovaps	%%zmm9 ,0x300(%%rsi)				\n\t	vmovaps	%%zmm13,0x380(%%rsi)	\n\t"\
		"vmovaps	%%zmm15,0x040(%%rsi)				\n\t	vmovaps	%%zmm10,0x0c0(%%rsi)	\n\t"\
		"vmovaps	%%zmm8 ,0x140(%%rsi)				\n\t	vmovaps	%%zmm17,0x1c0(%%rsi)	\n\t"\
		"\n\t"\
	/***************************************/\
		"\n\t"\
		/*...Block 4: */\
														"		subq	%%rdi,%%r9	\n\t"\
		"movslq	%[__p03],%%rdi							\n\t	addq	%%rdi,%%r9	\n\t"\
		/* SSE2_RADIX4_DIF_4TWIDDLE_B(r30,c03) */				/* SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO_B(r38,c07) */\
		"addq	$0x400,%%rsi	/* r30 */				\n\t"\
		"vmovaps	     (%%rax,%%rdi,8),%%zmm0			\n\t	vmovaps	     (%%rax,%%r9,8),%%zmm8 	\n\t"\
		"vmovaps	     (%%rcx,%%rdi,8),%%zmm4			\n\t	vmovaps	     (%%rcx,%%r9,8),%%zmm12	\n\t"\
		"vmovaps	0x040(%%rax,%%rdi,8),%%zmm1			\n\t	vmovaps	0x040(%%rax,%%r9,8),%%zmm9 	\n\t"\
		"vmovaps	0x040(%%rcx,%%rdi,8),%%zmm5			\n\t	vmovaps	0x040(%%rcx,%%r9,8),%%zmm13	\n\t"\
		"vmovaps	0x11c0(%%rsi),%%zmm6				\n\t	vmovaps	0x13c0(%%rsi),%%zmm14	\n\t"\
		"vmovaps	0x1200(%%rsi),%%zmm7				\n\t	vmovaps	0x1400(%%rsi),%%zmm15	\n\t"\
		"vmovaps	%%zmm0,%%zmm2						\n\t	vmovaps	%%zmm8 ,%%zmm10	\n\t"\
		"vmovaps	%%zmm1,%%zmm3						\n\t	vmovaps	%%zmm9 ,%%zmm11	\n\t"\
		"vmovaps	0x1240(%%rsi),%%zmm16		\n\t	vmovaps		0x1440(%%rsi),%%zmm20		\n\t"\
		"vmovaps	0x1280(%%rsi),%%zmm17		\n\t	vmovaps		0x1480(%%rsi),%%zmm21		\n\t"\
		"vmovaps	0x1340(%%rsi),%%zmm18		\n\t	vmovaps		0x1540(%%rsi),%%zmm22		\n\t"\
		"vmovaps	0x1380(%%rsi),%%zmm19		\n\t	vmovaps		0x1580(%%rsi),%%zmm23		\n\t"\
/* c03 */"	vmulpd		%%zmm6,%%zmm0,%%zmm0		\n\t		vmulpd		%%zmm14,%%zmm8 ,%%zmm8 	\n\t"/* c07 */\
		"	vmulpd		%%zmm6,%%zmm1,%%zmm1		\n\t		vmulpd		%%zmm14,%%zmm9 ,%%zmm9 	\n\t"\
		"vfnmadd231pd	%%zmm7,%%zmm3,%%zmm0		\n\t	vfnmadd231pd	%%zmm15,%%zmm11,%%zmm8 	\n\t"\
		" vfmadd231pd	%%zmm7,%%zmm2,%%zmm1		\n\t	 vfmadd231pd	%%zmm15,%%zmm10,%%zmm9 	\n\t"\
		"vmovaps	%%zmm4,%%zmm6						\n\t	vmovaps	%%zmm12,%%zmm14	\n\t"\
		"vmovaps	%%zmm5,%%zmm7						\n\t	vmovaps	%%zmm13,%%zmm15	\n\t"\
/* c13 */"	vmulpd		%%zmm16,%%zmm4,%%zmm4		\n\t		vmulpd		%%zmm20,%%zmm12,%%zmm12	\n\t"/* c17 */\
		"	vmulpd		%%zmm16,%%zmm5,%%zmm5		\n\t		vmulpd		%%zmm20,%%zmm13,%%zmm13	\n\t"\
		"vfnmadd231pd	%%zmm17,%%zmm7,%%zmm4		\n\t	vfnmadd231pd	%%zmm21,%%zmm15,%%zmm12	\n\t"\
		" vfmadd231pd	%%zmm17,%%zmm6,%%zmm5		\n\t	 vfmadd231pd	%%zmm21,%%zmm14,%%zmm13	\n\t"\
		"vmovaps	%%zmm0,%%zmm2						\n\t	vmovaps	%%zmm8 ,%%zmm10	\n\t"\
		"vmovaps	%%zmm1,%%zmm3						\n\t	vmovaps	%%zmm9 ,%%zmm11	\n\t"\
		"vaddpd	%%zmm4,%%zmm0,%%zmm0			\n\t	vaddpd	%%zmm12,%%zmm8 ,%%zmm8 	\n\t"\
		"vaddpd	%%zmm5,%%zmm1,%%zmm1			\n\t	vaddpd	%%zmm13,%%zmm9 ,%%zmm9 	\n\t"\
		"vsubpd	%%zmm4,%%zmm2,%%zmm2			\n\t	vsubpd	%%zmm12,%%zmm10,%%zmm10	\n\t"\
		"vsubpd	%%zmm5,%%zmm3,%%zmm3			\n\t	vsubpd	%%zmm13,%%zmm11,%%zmm11	\n\t"\
		"vmovaps	     (%%rdx,%%rdi,8),%%zmm6			\n\t	vmovaps	     (%%rdx,%%r9,8),%%zmm14	\n\t"\
		"vmovaps	0x040(%%rdx,%%rdi,8),%%zmm7			\n\t	vmovaps	0x040(%%rdx,%%r9,8),%%zmm15	\n\t"\
/* c1B */"	vmulpd		%%zmm18,%%zmm6,%%zmm4		\n\t		vmulpd		%%zmm22,%%zmm14,%%zmm12	\n\t"/* c1F */\
		"	vmulpd		%%zmm18,%%zmm7,%%zmm5		\n\t		vmulpd		%%zmm22,%%zmm15,%%zmm13	\n\t"\
		"vfnmadd231pd	%%zmm19,%%zmm7,%%zmm4		\n\t	vfnmadd231pd	%%zmm23,%%zmm15,%%zmm12	\n\t"\
		" vfmadd231pd	%%zmm19,%%zmm6,%%zmm5		\n\t	 vfmadd231pd	%%zmm23,%%zmm14,%%zmm13	\n\t"\
		"vmovaps	0x12c0(%%rsi),%%zmm16		\n\t	vmovaps		0x14c0(%%rsi),%%zmm20		\n\t"\
		"vmovaps	0x1300(%%rsi),%%zmm17		\n\t	vmovaps		0x1500(%%rsi),%%zmm21		\n\t"\
		"vmovaps	%%zmm5,0x040(%%rsi)					\n\t	vmovaps	%%zmm13,0x240(%%rsi)	\n\t"\
		"vmovaps	%%zmm4,     (%%rsi)					\n\t	vmovaps	%%zmm12,0x200(%%rsi)	/* tmpstr r08 */\n\t"\
		"vmovaps	     (%%rbx,%%rdi,8),%%zmm6			\n\t	vmovaps	     (%%rbx,%%r9,8),%%zmm14	\n\t"\
		"vmovaps	0x040(%%rbx,%%rdi,8),%%zmm7			\n\t	vmovaps	0x040(%%rbx,%%r9,8),%%zmm15	\n\t"\
/* c0B */"	vmulpd		%%zmm16,%%zmm6,%%zmm4		\n\t		vmulpd		%%zmm20,%%zmm14,%%zmm12	\n\t"/* c0F */\
		"	vmulpd		%%zmm16,%%zmm7,%%zmm5		\n\t		vmulpd		%%zmm20,%%zmm15,%%zmm13	\n\t"\
		"vfnmadd231pd	%%zmm17,%%zmm7,%%zmm4		\n\t	vfnmadd231pd	%%zmm21,%%zmm15,%%zmm12	\n\t"\
		" vfmadd231pd	%%zmm17,%%zmm6,%%zmm5		\n\t	 vfmadd231pd	%%zmm21,%%zmm14,%%zmm13	\n\t"\
		"vmovaps	     (%%rsi),%%zmm16			\n\t	vmovaps	0x200(%%rsi),%%zmm18	\n\t"\
		"vmovaps	0x040(%%rsi),%%zmm17			\n\t	vmovaps	0x240(%%rsi),%%zmm19	\n\t"\
		"vaddpd		%%zmm16,%%zmm4,%%zmm6				\n\t	vaddpd		%%zmm18,%%zmm12,%%zmm14	\n\t"\
		"vaddpd		%%zmm17,%%zmm5,%%zmm7				\n\t	vaddpd		%%zmm19,%%zmm13,%%zmm15	\n\t"\
		"vsubpd		%%zmm16,%%zmm4,%%zmm4				\n\t	vsubpd		%%zmm18,%%zmm12,%%zmm12	\n\t"\
		"vsubpd		%%zmm17,%%zmm5,%%zmm5				\n\t	vsubpd		%%zmm19,%%zmm13,%%zmm13	\n\t"\
		"vsubpd	%%zmm6,%%zmm0,%%zmm0			\n\t	vsubpd	    %%zmm14,%%zmm8 ,%%zmm8 	\n\t"\
		"vsubpd	%%zmm5,%%zmm2,%%zmm2			\n\t	vsubpd	    %%zmm15,%%zmm9 ,%%zmm9 	\n\t"\
		"vsubpd	%%zmm7,%%zmm1,%%zmm1			\n\t	vsubpd	%%zmm13,%%zmm10,%%zmm10	\n\t"\
		"vsubpd	%%zmm4,%%zmm3,%%zmm3			\n\t	vsubpd	%%zmm12,%%zmm11,%%zmm11	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm0,%%zmm6			\n\t	vfmadd132pd	%%zmm31,%%zmm8 ,%%zmm14	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm2,%%zmm5			\n\t	vfmadd132pd	%%zmm31,%%zmm10,%%zmm13	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm1,%%zmm7			\n\t	vfmadd132pd	%%zmm31,%%zmm9 ,%%zmm15	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm3,%%zmm4			\n\t	vfmadd132pd	%%zmm31,%%zmm11,%%zmm12	\n\t"\
		/* SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r30) */\
		"														vsubpd	%%zmm12,%%zmm10,%%zmm16	\n\t"\
		"														vsubpd	%%zmm11,%%zmm13,%%zmm17	\n\t"\
		"														vaddpd	%%zmm12,%%zmm10,%%zmm10	\n\t"\
		"														vaddpd	%%zmm11,%%zmm13,%%zmm13	\n\t"\
		"														vfnmadd231pd %%zmm29,%%zmm16,%%zmm2	\n\t"\
		"														vfnmadd231pd %%zmm29,%%zmm17,%%zmm3	\n\t"\
		"vsubpd	%%zmm14,%%zmm6,%%zmm6			\n\t	vfnmadd231pd %%zmm29,%%zmm10,%%zmm4	\n\t"\
		"vsubpd	%%zmm9 ,%%zmm0,%%zmm0			\n\t	vfnmadd231pd %%zmm29,%%zmm13,%%zmm5	\n\t"\
		"vsubpd	%%zmm8 ,%%zmm1 ,%%zmm1 					\n\t	vsubpd	%%zmm15,%%zmm7 ,%%zmm7 	\n\t"\
		"vmovaps	%%zmm6 ,0x200(%%rsi)				\n\t	vmovaps	%%zmm2 ,0x280(%%rsi)	\n\t"\
		"vmovaps	%%zmm0 ,0x100(%%rsi)				\n\t	vmovaps	%%zmm5 ,0x180(%%rsi)	\n\t"\
		"vmovaps	%%zmm7 ,0x240(%%rsi)				\n\t	vmovaps	%%zmm4 ,0x2c0(%%rsi)	\n\t"\
		"vmovaps	%%zmm1 ,0x340(%%rsi)				\n\t	vmovaps	%%zmm3 ,0x3c0(%%rsi)	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm6 ,%%zmm14			\n\t	vfmadd132pd	%%zmm28,%%zmm2 ,%%zmm16	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm0 ,%%zmm9 			\n\t	vfmadd132pd	%%zmm28,%%zmm5 ,%%zmm13	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm7 ,%%zmm15			\n\t	vfmadd132pd	%%zmm28,%%zmm4 ,%%zmm10	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm1 ,%%zmm8 			\n\t	vfmadd132pd	%%zmm28,%%zmm3 ,%%zmm17	\n\t"\
		"vmovaps	%%zmm14,     (%%rsi)				\n\t	vmovaps	%%zmm16,0x080(%%rsi)	\n\t"\
		"vmovaps	%%zmm9 ,0x300(%%rsi)				\n\t	vmovaps	%%zmm13,0x380(%%rsi)	\n\t"\
		"vmovaps	%%zmm15,0x040(%%rsi)				\n\t	vmovaps	%%zmm10,0x0c0(%%rsi)	\n\t"\
		"vmovaps	%%zmm8 ,0x140(%%rsi)				\n\t	vmovaps	%%zmm17,0x1c0(%%rsi)	\n\t"\
		"\n\t"\
	/**********************************************************************************/\
	/*...and now do eight radix-4 transforms, including the internal twiddle factors: */\
	/**********************************************************************************/\
		/*...Block 1: t00,t10,t20,t30	*/\
		"movq	%[__add0],%%rax	/* &a[j1] */			\n\t"\
		"movslq	%[__p01],%%rbx							\n\t"\
		"movslq	%[__p02],%%rcx							\n\t"\
		"movslq	%[__p03],%%rdx							\n\t"		/*...Block 5: t08,t18,t28,t38	*/\
		"leaq	(%%rax,%%rbx,8),%%rbx					\n\t		subq	%%rdi,%%r9	\n\t"\
		"leaq	(%%rax,%%rcx,8),%%rcx					\n\t"\
		"leaq	(%%rax,%%rdx,8),%%rdx					\n\t"\
		"movq	%[__r00],%%rsi							\n\t"\
		"vmovaps	     (%%rsi),%%zmm0					\n\t	vmovaps	0xa00(%%rsi),%%zmm12	\n\t"\
		"vmovaps	0x800(%%rsi),%%zmm4					\n\t	vmovaps	0xa40(%%rsi),%%zmm13	\n\t"\
		"vmovaps	0x040(%%rsi),%%zmm1					\n\t	vmovaps	0xe00(%%rsi),%%zmm14	\n\t"\
		"vmovaps	0x840(%%rsi),%%zmm5					\n\t	vmovaps	0xe40(%%rsi),%%zmm15	\n\t"\
		"vmovaps	0x400(%%rsi),%%zmm2					\n\t	vmulpd	%%zmm29,%%zmm12,%%zmm12	\n\t"/* *isrt2 */\
		"vmovaps	0xc00(%%rsi),%%zmm6					\n\t	vmulpd	%%zmm29,%%zmm13,%%zmm13	\n\t"\
		"vmovaps	0x440(%%rsi),%%zmm3					\n\t	vmulpd	%%zmm29,%%zmm14,%%zmm14	\n\t"\
		"vmovaps	0xc40(%%rsi),%%zmm7					\n\t	vmulpd	%%zmm29,%%zmm15,%%zmm15	\n\t"\
		"vsubpd	0x400(%%rsi),%%zmm0,%%zmm0				\n\t	vmovaps	0x200(%%rsi),%%zmm8 	\n\t"\
		"vsubpd	0xc00(%%rsi),%%zmm4,%%zmm4				\n\t	vmovaps	0x240(%%rsi),%%zmm9 	\n\t"\
		"vsubpd	0x440(%%rsi),%%zmm1,%%zmm1				\n\t	vmovaps	0x600(%%rsi),%%zmm10	\n\t"\
		"vsubpd	0xc40(%%rsi),%%zmm5,%%zmm5				\n\t	vmovaps	0x640(%%rsi),%%zmm11	\n\t"\
		"vaddpd	     (%%rsi),%%zmm2,%%zmm2				\n\t	vsubpd	%%zmm11,%%zmm8 ,%%zmm8 	\n\t"\
		"vaddpd	0x800(%%rsi),%%zmm6,%%zmm6				\n\t	vsubpd	%%zmm13,%%zmm12,%%zmm12	\n\t"\
		"vaddpd	0x040(%%rsi),%%zmm3,%%zmm3				\n\t	vsubpd	%%zmm10,%%zmm9 ,%%zmm9 	\n\t"\
		"vaddpd	0x840(%%rsi),%%zmm7,%%zmm7				\n\t	vsubpd	%%zmm14,%%zmm15,%%zmm15	\n\t"\
		"vsubpd	%%zmm6,%%zmm2,%%zmm2			\n\t	vfmadd132pd	%%zmm31,%%zmm8 ,%%zmm11	\n\t"\
		"vsubpd	%%zmm5,%%zmm0,%%zmm0			\n\t	vfmadd132pd	%%zmm31,%%zmm12,%%zmm13	\n\t"\
		"vsubpd	%%zmm7,%%zmm3,%%zmm3			\n\t	vfmadd132pd	%%zmm31,%%zmm9 ,%%zmm10	\n\t"\
		"vsubpd	%%zmm4,%%zmm1,%%zmm1			\n\t	vfmadd132pd	%%zmm31,%%zmm15,%%zmm14	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm2,%%zmm6			\n\t	vsubpd	%%zmm14,%%zmm12,%%zmm12	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm0,%%zmm5			\n\t	vsubpd	%%zmm15,%%zmm13,%%zmm13	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm3,%%zmm7			\n\t	 vfmadd132pd	%%zmm31,%%zmm12,%%zmm14	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm1,%%zmm4			\n\t	 vfmadd132pd	%%zmm31,%%zmm13,%%zmm15	\n\t"\
		"vmovaps	%%zmm2,     (%%rbx)					\n\t	vsubpd	%%zmm12,%%zmm8 ,%%zmm8 	\n\t"\
		"vmovaps	%%zmm0,     (%%rcx)					\n\t	vsubpd	%%zmm13,%%zmm10,%%zmm10	\n\t"\
		"vmovaps	%%zmm3,0x040(%%rbx)					\n\t	vsubpd	%%zmm14,%%zmm9 ,%%zmm9 	\n\t"\
		"vmovaps	%%zmm1,0x040(%%rdx)					\n\t	vsubpd	%%zmm15,%%zmm11,%%zmm11	\n\t"\
		"vmovaps	%%zmm6,     (%%rax)					\n\t	vfmadd132pd	%%zmm31,%%zmm8 ,%%zmm12	\n\t"\
		"vmovaps	%%zmm5,     (%%rdx)					\n\t	vfmadd132pd	%%zmm31,%%zmm11,%%zmm15	\n\t"\
		"vmovaps	%%zmm7,0x040(%%rax)					\n\t	vfmadd132pd	%%zmm31,%%zmm10,%%zmm13	\n\t"\
		"vmovaps	%%zmm4,0x040(%%rcx)					\n\t	vfmadd132pd	%%zmm31,%%zmm9 ,%%zmm14	\n\t"\
		"												\n\t	vmovaps	%%zmm8 ,     (%%rbx,%%r9,8)	\n\t"\
		"												\n\t	vmovaps	%%zmm11,     (%%rcx,%%r9,8)	\n\t"\
		"												\n\t	vmovaps	%%zmm10,0x040(%%rbx,%%r9,8)	\n\t"\
		"												\n\t	vmovaps	%%zmm9 ,0x040(%%rdx,%%r9,8)	\n\t"\
		"												\n\t	vmovaps	%%zmm12,     (%%rax,%%r9,8)	\n\t"\
		"												\n\t	vmovaps	%%zmm15,     (%%rdx,%%r9,8)	\n\t"\
		"												\n\t	vmovaps	%%zmm13,0x040(%%rax,%%r9,8)	\n\t"\
		"												\n\t	vmovaps	%%zmm14,0x040(%%rcx,%%r9,8)	\n\t"\
	/*...Block 3: t04,t14,t24,t34	*/\
		"addq	$0x100,%%rsi	/* r04 */				\n\t"\
		"movslq	%[__p08],%%rdi							\n\t		addq	%%rdi,%%r9	\n\t"\
		"vmovaps	0x800(%%rsi),%%zmm4					\n\t	vmovaps	0xa00(%%rsi),%%zmm12	\n\t"\
		"vmovaps	0x840(%%rsi),%%zmm5					\n\t	vmovaps	0xa40(%%rsi),%%zmm13	\n\t"\
		"vmovaps	0xf40(%%rsi),%%zmm3	/* cc0 */		\n\t"\
		"vmovaps	0xf80(%%rsi),%%zmm2					\n\t"\
		"vmovaps	%%zmm4,%%zmm6						\n\t	vmovaps	%%zmm12,%%zmm14	\n\t"\
		"vmovaps	%%zmm5,%%zmm7						\n\t	vmovaps	%%zmm13,%%zmm15	\n\t"\
		"vmulpd	%%zmm3,%%zmm4,%%zmm4					\n\t	vmulpd	%%zmm2 ,%%zmm12,%%zmm12	\n\t"\
		"vmulpd	%%zmm3,%%zmm5,%%zmm5					\n\t	vmulpd	%%zmm2 ,%%zmm13,%%zmm13	\n\t"\
		"vmovaps	0xc00(%%rsi),%%zmm0					\n\t	vmovaps	0xe00(%%rsi),%%zmm8 	\n\t"\
		"vmovaps	0xc40(%%rsi),%%zmm1					\n\t	vmovaps	0xe40(%%rsi),%%zmm9 	\n\t"\
		"vfnmadd231pd	%%zmm2,%%zmm7,%%zmm4			\n\t	vfnmadd231pd	%%zmm3 ,%%zmm15,%%zmm12	\n\t"\
		" vfmadd231pd	%%zmm2,%%zmm6,%%zmm5			\n\t	 vfmadd231pd	%%zmm3 ,%%zmm14,%%zmm13	\n\t"\
		"vmovaps	%%zmm0,%%zmm6						\n\t	vmovaps	%%zmm8 ,%%zmm14	\n\t"\
		"vmovaps	%%zmm1,%%zmm7						\n\t	vmovaps	%%zmm9 ,%%zmm15	\n\t"\
		"vmulpd	%%zmm2,%%zmm6,%%zmm6					\n\t	vmulpd	%%zmm3 ,%%zmm14,%%zmm14	\n\t"\
		"vmulpd	%%zmm2,%%zmm7,%%zmm7					\n\t	vmulpd	%%zmm3 ,%%zmm15,%%zmm15	\n\t"\
		"vfnmadd231pd	%%zmm3,%%zmm1,%%zmm6			\n\t	vfnmadd231pd	%%zmm2 ,%%zmm9 ,%%zmm14	\n\t"\
		" vfmadd231pd	%%zmm3,%%zmm0,%%zmm7			\n\t	 vfmadd231pd	%%zmm2 ,%%zmm8 ,%%zmm15	\n\t"\
		"vmovaps	%%zmm4,%%zmm2						\n\t	vmovaps	%%zmm12,%%zmm10	\n\t"\
		"vmovaps	%%zmm5,%%zmm3						\n\t	vmovaps	%%zmm13,%%zmm11	\n\t"\
		"vsubpd	%%zmm6,%%zmm4,%%zmm4			\n\t	vsubpd	%%zmm14,%%zmm12,%%zmm12	\n\t"\
		"vsubpd	%%zmm7,%%zmm5,%%zmm5			\n\t	vsubpd	%%zmm15,%%zmm13,%%zmm13	\n\t"\
		"vaddpd	%%zmm2,%%zmm6,%%zmm6			\n\t	vaddpd	%%zmm10,%%zmm14,%%zmm14	\n\t"\
		"vaddpd	%%zmm3,%%zmm7,%%zmm7			\n\t	vaddpd	%%zmm11,%%zmm15,%%zmm15	\n\t"\
		"vmovaps	0x400(%%rsi),%%zmm2					\n\t	vmovaps	0x600(%%rsi),%%zmm10	\n\t"\
		"vmovaps	0x440(%%rsi),%%zmm3					\n\t	vmovaps	0x640(%%rsi),%%zmm11	\n\t"\
		"vmovaps	%%zmm2,%%zmm0						\n\t	vmovaps	%%zmm10,%%zmm8 	\n\t"\
		"vsubpd	%%zmm3,%%zmm2,%%zmm2			\n\t	vaddpd	%%zmm11,%%zmm10,%%zmm10	\n\t"\
		"vaddpd	%%zmm0,%%zmm3,%%zmm3			\n\t	vsubpd	%%zmm8 ,%%zmm11,%%zmm11	\n\t"\
		"vmulpd	%%zmm29,%%zmm2,%%zmm2					\n\t	vmulpd	%%zmm29,%%zmm10,%%zmm10	\n\t"\
		"vmulpd	%%zmm29,%%zmm3,%%zmm3					\n\t	vmulpd	%%zmm29,%%zmm11,%%zmm11	\n\t"\
		"vmovaps	     (%%rsi),%%zmm0					\n\t	vmovaps	0x200(%%rsi),%%zmm8 	\n\t"\
		"vmovaps	0x040(%%rsi),%%zmm1					\n\t	vmovaps	0x240(%%rsi),%%zmm9 	\n\t"\
		"vsubpd	%%zmm2,%%zmm0,%%zmm0			\n\t	vsubpd	%%zmm10,%%zmm8,%%zmm8  	\n\t"\
		"vsubpd	%%zmm3,%%zmm1,%%zmm1			\n\t	vsubpd	%%zmm11,%%zmm9,%%zmm9  	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm0,%%zmm2			\n\t	vfmadd132pd	%%zmm31,%%zmm8 ,%%zmm10	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm1,%%zmm3			\n\t	vfmadd132pd	%%zmm31,%%zmm9 ,%%zmm11	\n\t"\
		"vsubpd	%%zmm6,%%zmm2,%%zmm2			\n\t	vsubpd	%%zmm12,%%zmm8 ,%%zmm8 	\n\t"\
		"vsubpd	%%zmm5,%%zmm0,%%zmm0			\n\t	vsubpd	%%zmm15,%%zmm10,%%zmm10	\n\t"\
		"vsubpd	%%zmm7,%%zmm3,%%zmm3			\n\t	vsubpd	%%zmm13,%%zmm9 ,%%zmm9 	\n\t"\
		"vsubpd	%%zmm4,%%zmm1,%%zmm1			\n\t	vsubpd	%%zmm14,%%zmm11,%%zmm11	\n\t"\
		"vmovaps	%%zmm2,     (%%rbx,%%rdi,8)			\n\t	vmovaps	%%zmm8 ,     (%%rbx,%%r9,8)	\n\t"\
		"vmovaps	%%zmm0,     (%%rcx,%%rdi,8)			\n\t	vmovaps	%%zmm10,     (%%rcx,%%r9,8)	\n\t"\
		"vmovaps	%%zmm3,0x040(%%rbx,%%rdi,8)			\n\t	vmovaps	%%zmm9 ,0x040(%%rbx,%%r9,8)	\n\t"\
		"vmovaps	%%zmm1,0x040(%%rdx,%%rdi,8)			\n\t	vmovaps	%%zmm11,0x040(%%rdx,%%r9,8)	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm2,%%zmm6			\n\t	vfmadd132pd	%%zmm31,%%zmm8 ,%%zmm12	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm0,%%zmm5			\n\t	vfmadd132pd	%%zmm31,%%zmm10,%%zmm15	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm3,%%zmm7			\n\t	vfmadd132pd	%%zmm31,%%zmm9 ,%%zmm13	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm1,%%zmm4			\n\t	vfmadd132pd	%%zmm31,%%zmm11,%%zmm14	\n\t"\
		"vmovaps	%%zmm6,     (%%rax,%%rdi,8)			\n\t	vmovaps	%%zmm12,     (%%rax,%%r9,8)	\n\t"\
		"vmovaps	%%zmm5,     (%%rdx,%%rdi,8)			\n\t	vmovaps	%%zmm15,     (%%rdx,%%r9,8)	\n\t"\
		"vmovaps	%%zmm7,0x040(%%rax,%%rdi,8)			\n\t	vmovaps	%%zmm13,0x040(%%rax,%%r9,8)	\n\t"\
		"vmovaps	%%zmm4,0x040(%%rcx,%%rdi,8)			\n\t	vmovaps	%%zmm14,0x040(%%rcx,%%r9,8)	\n\t"\
		"\n\t"\
	/*...Block 2: t02,t12,t22,t32	*/\
		"subq	$0x80,%%rsi	/* r02 */					\n\t		subq	%%rdi,%%r9	\n\t"\
		"movslq	%[__p10],%%rdi							\n\t		addq	%%rdi,%%r9	\n\t"\
		"vmovaps	0x800(%%rsi),%%zmm4					\n\t	vmovaps	0xa00(%%rsi),%%zmm12	\n\t"\
		"vmovaps	0x840(%%rsi),%%zmm5					\n\t	vmovaps	0xa40(%%rsi),%%zmm13	\n\t"\
		"vmovaps	0x1040(%%rsi),%%zmm2	/* cc1 */	\n\t	vmovaps	0x10c0(%%rsi),%%zmm11	/* cc3 */	\n\t"\
		"vmovaps	0x1080(%%rsi),%%zmm3				\n\t	vmovaps	0x1100(%%rsi),%%zmm10	\n\t"\
		"vmovaps	%%zmm4,%%zmm6						\n\t	vmovaps	%%zmm12,%%zmm14	\n\t"\
		"vmovaps	%%zmm5,%%zmm7						\n\t	vmovaps	%%zmm13,%%zmm15	\n\t"\
		"vmulpd	%%zmm2,%%zmm4,%%zmm4					\n\t	vmulpd	%%zmm10,%%zmm12,%%zmm12	\n\t"\
		"vmulpd	%%zmm2,%%zmm5,%%zmm5					\n\t	vmulpd	%%zmm10,%%zmm13,%%zmm13	\n\t"\
		"vmovaps	0xc00(%%rsi),%%zmm0					\n\t	vmovaps	0xe00(%%rsi),%%zmm8 	\n\t"\
		"vmovaps	0xc40(%%rsi),%%zmm1					\n\t	vmovaps	0xe40(%%rsi),%%zmm9 	\n\t"\
		"vfnmadd231pd	%%zmm3,%%zmm7,%%zmm4			\n\t	vfnmadd231pd	%%zmm11,%%zmm15,%%zmm12	\n\t"\
		" vfmadd231pd	%%zmm3,%%zmm6,%%zmm5			\n\t	 vfmadd231pd	%%zmm11,%%zmm14,%%zmm13	\n\t"\
		"vmovaps	%%zmm0,%%zmm6						\n\t	vmovaps	%%zmm8 ,%%zmm14	\n\t"\
		"vmovaps	%%zmm1,%%zmm7						\n\t	vmovaps	%%zmm9 ,%%zmm15	\n\t"\
		"vmulpd	%%zmm11,%%zmm6,%%zmm6					\n\t	vmulpd	%%zmm2 ,%%zmm14,%%zmm14	\n\t"\
		"vmulpd	%%zmm11,%%zmm7,%%zmm7					\n\t	vmulpd	%%zmm2 ,%%zmm15,%%zmm15	\n\t"\
		"vfnmadd231pd	%%zmm10,%%zmm1,%%zmm6			\n\t	 vfmadd231pd	%%zmm3 ,%%zmm9 ,%%zmm14	\n\t"\
		" vfmadd231pd	%%zmm10,%%zmm0,%%zmm7			\n\t	vfnmadd231pd	%%zmm3 ,%%zmm8 ,%%zmm15	\n\t"\
		"vmovaps	%%zmm4,%%zmm2						\n\t	vmovaps	%%zmm12,%%zmm10	\n\t"\
		"vmovaps	%%zmm5,%%zmm3						\n\t	vmovaps	%%zmm13,%%zmm11	\n\t"\
		"vsubpd	%%zmm6,%%zmm4,%%zmm4			\n\t	vsubpd	%%zmm14,%%zmm12,%%zmm12	\n\t"\
		"vsubpd	%%zmm7,%%zmm5,%%zmm5			\n\t	vsubpd	%%zmm15,%%zmm13,%%zmm13	\n\t"\
		"vaddpd	%%zmm2,%%zmm6,%%zmm6			\n\t	vaddpd	%%zmm10,%%zmm14,%%zmm14	\n\t"\
		"vaddpd	%%zmm3,%%zmm7,%%zmm7			\n\t	vaddpd	%%zmm11,%%zmm15,%%zmm15	\n\t"\
		"vmovaps	0x400(%%rsi),%%zmm1					\n\t	vmovaps	0x600(%%rsi),%%zmm9 	\n\t"\
		"vmovaps	0x440(%%rsi),%%zmm3					\n\t	vmovaps	0x640(%%rsi),%%zmm11	\n\t"\
		"vmovaps	0x1000(%%rsi),%%zmm0	/* ss0 */	\n\t	vmovaps	0xfc0(%%rsi),%%zmm8 		/* cc0 */	\n\t"\
		"vmovaps	%%zmm1,%%zmm2						\n\t	vmovaps	%%zmm9 ,%%zmm10	\n\t"\
		"vmulpd	      %%zmm0,%%zmm1,%%zmm1				\n\t	vmulpd	     %%zmm8 ,%%zmm9 ,%%zmm9 	\n\t"\
		"vmulpd	      %%zmm3,%%zmm0,%%zmm0				\n\t	vmulpd	     %%zmm11,%%zmm8 ,%%zmm8 	\n\t"\
		"vfmsub132pd	0xfc0(%%rsi),%%zmm0,%%zmm2		\n\t	vfmadd132pd	0x1000(%%rsi),%%zmm8 ,%%zmm10	\n\t"\
		"vfmadd132pd	0xfc0(%%rsi),%%zmm1,%%zmm3		\n\t	vfmsub132pd	0x1000(%%rsi),%%zmm9 ,%%zmm11	\n\t"\
		"vmovaps	     (%%rsi),%%zmm0					\n\t	vmovaps	0x200(%%rsi),%%zmm8 	\n\t"\
		"vmovaps	0x040(%%rsi),%%zmm1					\n\t	vmovaps	0x240(%%rsi),%%zmm9 	\n\t"\
		"vsubpd	%%zmm2,%%zmm0,%%zmm0			\n\t	vsubpd	%%zmm10,%%zmm8 ,%%zmm8  	\n\t"\
		"vsubpd	%%zmm3,%%zmm1,%%zmm1			\n\t	vsubpd	%%zmm11,%%zmm9 ,%%zmm9  	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm0,%%zmm2			\n\t	vfmadd132pd	%%zmm31,%%zmm8 ,%%zmm10	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm1,%%zmm3			\n\t	vfmadd132pd	%%zmm31,%%zmm9 ,%%zmm11	\n\t"\
		"vsubpd	%%zmm6,%%zmm2,%%zmm2			\n\t	vsubpd	%%zmm12,%%zmm8 ,%%zmm8 	\n\t"\
		"vsubpd	%%zmm5,%%zmm0,%%zmm0			\n\t	vsubpd	%%zmm15,%%zmm10,%%zmm10	\n\t"\
		"vsubpd	%%zmm7,%%zmm3,%%zmm3			\n\t	vsubpd	%%zmm13,%%zmm9 ,%%zmm9 	\n\t"\
		"vsubpd	%%zmm4,%%zmm1,%%zmm1			\n\t	vsubpd	%%zmm14,%%zmm11,%%zmm11	\n\t"\
		"vmovaps	%%zmm2,     (%%rbx,%%rdi,8)			\n\t	vmovaps	%%zmm8 ,     (%%rbx,%%r9,8)	\n\t"\
		"vmovaps	%%zmm0,     (%%rcx,%%rdi,8)			\n\t	vmovaps	%%zmm10,     (%%rcx,%%r9,8)	\n\t"\
		"vmovaps	%%zmm3,0x040(%%rbx,%%rdi,8)			\n\t	vmovaps	%%zmm9 ,0x040(%%rbx,%%r9,8)	\n\t"\
		"vmovaps	%%zmm1,0x040(%%rdx,%%rdi,8)			\n\t	vmovaps	%%zmm11,0x040(%%rdx,%%r9,8)	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm2,%%zmm6			\n\t	vfmadd132pd	%%zmm31,%%zmm8 ,%%zmm12	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm0,%%zmm5			\n\t	vfmadd132pd	%%zmm31,%%zmm10,%%zmm15	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm3,%%zmm7			\n\t	vfmadd132pd	%%zmm31,%%zmm9 ,%%zmm13	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm1,%%zmm4			\n\t	vfmadd132pd	%%zmm31,%%zmm11,%%zmm14	\n\t"\
		"vmovaps	%%zmm6,     (%%rax,%%rdi,8)			\n\t	vmovaps	%%zmm12,     (%%rax,%%r9,8)	\n\t"\
		"vmovaps	%%zmm5,     (%%rdx,%%rdi,8)			\n\t	vmovaps	%%zmm15,     (%%rdx,%%r9,8)	\n\t"\
		"vmovaps	%%zmm7,0x040(%%rax,%%rdi,8)			\n\t	vmovaps	%%zmm13,0x040(%%rax,%%r9,8)	\n\t"\
		"vmovaps	%%zmm4,0x040(%%rcx,%%rdi,8)			\n\t	vmovaps	%%zmm14,0x040(%%rcx,%%r9,8)	\n\t"\
		"\n\t"\
	/*...Block 4: t06,t16,t26,t36	*/\
		"addq	$0x100,%%rsi	/* r06 */				\n\t		subq	%%rdi,%%r9	\n\t"\
		"movslq	%[__p18],%%rdi							\n\t		addq	%%rdi,%%r9	\n\t"\
		"vmovaps	0x800(%%rsi),%%zmm4					\n\t	vmovaps	0xa00(%%rsi),%%zmm12	\n\t"\
		"vmovaps	0x840(%%rsi),%%zmm5					\n\t	vmovaps	0xa40(%%rsi),%%zmm13	\n\t"\
		"vmovaps	0xfc0(%%rsi),%%zmm2	/* cc3 */		\n\t	vmovaps	0xf40(%%rsi),%%zmm11 /* cc1 */	\n\t"\
		"vmovaps	0x1000(%%rsi),%%zmm3				\n\t	vmovaps	0xf80(%%rsi),%%zmm10 \n\t"\
		"vmovaps	%%zmm4,%%zmm6						\n\t	vmovaps	%%zmm12,%%zmm14	\n\t"\
		"vmovaps	%%zmm5,%%zmm7						\n\t	vmovaps	%%zmm13,%%zmm15	\n\t"\
		"vmulpd	%%zmm2,%%zmm4,%%zmm4					\n\t	vmulpd	%%zmm10,%%zmm12,%%zmm12	\n\t"\
		"vmulpd	%%zmm2,%%zmm5,%%zmm5					\n\t	vmulpd	%%zmm10,%%zmm13,%%zmm13	\n\t"\
		"vfnmadd231pd	%%zmm3,%%zmm7,%%zmm4			\n\t	vfnmadd231pd	%%zmm11,%%zmm15,%%zmm12	\n\t"\
		" vfmadd231pd	%%zmm3,%%zmm6,%%zmm5			\n\t	 vfmadd231pd	%%zmm11,%%zmm14,%%zmm13	\n\t"\
		"vmovaps	0xc00(%%rsi),%%zmm0					\n\t	vmovaps	0xe00(%%rsi),%%zmm8 	\n\t"\
		"vmovaps	0xc40(%%rsi),%%zmm1					\n\t	vmovaps	0xe40(%%rsi),%%zmm9 	\n\t"\
		"vmovaps	%%zmm0,%%zmm6						\n\t	vmovaps	%%zmm8 ,%%zmm14	\n\t"\
		"vmovaps	%%zmm1,%%zmm7						\n\t	vmovaps	%%zmm9 ,%%zmm15	\n\t"\
		"vmulpd	%%zmm10,%%zmm6,%%zmm6					\n\t	vmulpd	%%zmm3 ,%%zmm14,%%zmm14	\n\t"\
		"vmulpd	%%zmm10,%%zmm7,%%zmm7					\n\t	vmulpd	%%zmm3 ,%%zmm15,%%zmm15	\n\t"\
		" vfmadd231pd	%%zmm11,%%zmm1,%%zmm6			\n\t	vfnmadd231pd	%%zmm2 ,%%zmm9 ,%%zmm14	\n\t"\
		"vfnmadd231pd	%%zmm11,%%zmm0,%%zmm7			\n\t	 vfmadd231pd	%%zmm2 ,%%zmm8 ,%%zmm15	\n\t"\
		"vmovaps	%%zmm4,%%zmm2						\n\t	vmovaps	%%zmm12,%%zmm10	\n\t"\
		"vmovaps	%%zmm5,%%zmm3						\n\t	vmovaps	%%zmm13,%%zmm11	\n\t"\
		"vsubpd	%%zmm6,%%zmm4,%%zmm4			\n\t	vsubpd	%%zmm14,%%zmm12,%%zmm12	\n\t"\
		"vsubpd	%%zmm7,%%zmm5,%%zmm5			\n\t	vsubpd	%%zmm15,%%zmm13,%%zmm13	\n\t"\
		"vaddpd	%%zmm2,%%zmm6,%%zmm6			\n\t	vaddpd	%%zmm10,%%zmm14,%%zmm14	\n\t"\
		"vaddpd	%%zmm3,%%zmm7,%%zmm7			\n\t	vaddpd	%%zmm11,%%zmm15,%%zmm15	\n\t"\
		"vmovaps	0x400(%%rsi),%%zmm1					\n\t	vmovaps	0x600(%%rsi),%%zmm9 	\n\t"\
		"vmovaps	0x440(%%rsi),%%zmm3					\n\t	vmovaps	0x640(%%rsi),%%zmm11	\n\t"\
		"vmovaps	0xec0(%%rsi),%%zmm0		/* cc0 */	\n\t	vmovaps	0xf00(%%rsi),%%zmm8 	/* ss0 */	\n\t"\
		"vmovaps	%%zmm1,%%zmm2						\n\t	vmovaps	%%zmm9 ,%%zmm10	\n\t"\
		"vmulpd	      %%zmm0,%%zmm1,%%zmm1				\n\t	vmulpd	     %%zmm8 ,%%zmm9 ,%%zmm9 	\n\t"\
		"vmulpd	      %%zmm3,%%zmm0,%%zmm0				\n\t	vmulpd	     %%zmm11,%%zmm8 ,%%zmm8 	\n\t"\
		"vfmsub132pd	0xf00(%%rsi),%%zmm0,%%zmm2		\n\t	vfmadd132pd	0xec0(%%rsi),%%zmm8 ,%%zmm10	\n\t"\
		"vfmadd132pd	0xf00(%%rsi),%%zmm1,%%zmm3		\n\t	vfmsub132pd	0xec0(%%rsi),%%zmm9 ,%%zmm11	\n\t"\
		"vmovaps	     (%%rsi),%%zmm0					\n\t	vmovaps	0x200(%%rsi),%%zmm8 	\n\t"\
		"vmovaps	0x040(%%rsi),%%zmm1					\n\t	vmovaps	0x240(%%rsi),%%zmm9 	\n\t"\
		"vsubpd	%%zmm2,%%zmm0,%%zmm0			\n\t	vsubpd	%%zmm10,%%zmm8 ,%%zmm8  	\n\t"\
		"vsubpd	%%zmm3,%%zmm1,%%zmm1			\n\t	vsubpd	%%zmm11,%%zmm9 ,%%zmm9  	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm0,%%zmm2			\n\t	vfmadd132pd	%%zmm31,%%zmm8 ,%%zmm10	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm1,%%zmm3			\n\t	vfmadd132pd	%%zmm31,%%zmm9 ,%%zmm11	\n\t"\
		"vsubpd	%%zmm4,%%zmm2,%%zmm2			\n\t	vsubpd	%%zmm12,%%zmm8 ,%%zmm8 	\n\t"\
		"vsubpd	%%zmm7,%%zmm0,%%zmm0			\n\t	vsubpd	%%zmm15,%%zmm10,%%zmm10	\n\t"\
		"vsubpd	%%zmm5,%%zmm3,%%zmm3			\n\t	vsubpd	%%zmm13,%%zmm9 ,%%zmm9 	\n\t"\
		"vsubpd	%%zmm6,%%zmm1,%%zmm1			\n\t	vsubpd	%%zmm14,%%zmm11,%%zmm11	\n\t"\
		"vmovaps	%%zmm2,     (%%rbx,%%rdi,8)			\n\t	vmovaps	%%zmm8 ,     (%%rbx,%%r9,8)	\n\t"\
		"vmovaps	%%zmm0,     (%%rcx,%%rdi,8)			\n\t	vmovaps	%%zmm10,     (%%rcx,%%r9,8)	\n\t"\
		"vmovaps	%%zmm3,0x040(%%rbx,%%rdi,8)			\n\t	vmovaps	%%zmm9 ,0x040(%%rbx,%%r9,8)	\n\t"\
		"vmovaps	%%zmm1,0x040(%%rdx,%%rdi,8)			\n\t	vmovaps	%%zmm11,0x040(%%rdx,%%r9,8)	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm2,%%zmm4			\n\t	vfmadd132pd	%%zmm31,%%zmm8 ,%%zmm12	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm0,%%zmm7			\n\t	vfmadd132pd	%%zmm31,%%zmm10,%%zmm15	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm3,%%zmm5			\n\t	vfmadd132pd	%%zmm31,%%zmm9 ,%%zmm13	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm1,%%zmm6			\n\t	vfmadd132pd	%%zmm31,%%zmm11,%%zmm14	\n\t"\
		"vmovaps	%%zmm4,     (%%rax,%%rdi,8)			\n\t	vmovaps	%%zmm12,     (%%rax,%%r9,8)	\n\t"\
		"vmovaps	%%zmm7,     (%%rdx,%%rdi,8)			\n\t	vmovaps	%%zmm15,     (%%rdx,%%r9,8)	\n\t"\
		"vmovaps	%%zmm5,0x040(%%rax,%%rdi,8)			\n\t	vmovaps	%%zmm13,0x040(%%rax,%%r9,8)	\n\t"\
		"vmovaps	%%zmm6,0x040(%%rcx,%%rdi,8)			\n\t	vmovaps	%%zmm14,0x040(%%rcx,%%r9,8)	\n\t"\
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
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r8","r9","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15","xmm16","xmm17","xmm18","xmm19","xmm20","xmm21","xmm22","xmm23","xmm24","xmm25","xmm26","xmm27","xmm28","xmm29","xmm30","xmm31"	/* Clobbered registers */\
	);\
	}

  #ifndef ALL_FMA
	#define ALL_FMA	0	// Jan 2021: See no discernable speed difference between the 2 variants,
  #endif				// so we hope that using the lower-FMA one slightly reduces power consumption

  #if ALL_FMA

	#warning Using ALL_FMA version of AVX-512 SSE2_RADIX32_DIT_TWIDDLE macro.
	// Cost [vector-ops only]: 374 MEM, 102 MUL, 407 FMA, 33 ADD
	#define SSE2_RADIX32_DIT_TWIDDLE(Xadd0,Xp01,Xp02,Xp03,Xp04,Xp05,Xp06,Xp07,Xp08,Xp10,Xr00,Xisrt2)\
	{\
	__asm__ volatile (\
		"movq	%[__add0],%%rax							\n\t	movq	%[__isrt2],%%r8		\n\t"\
		"movslq	%[__p01],%%rbx							\n\t	movslq	%[__p04],%%r9		\n\t		movq	%[__r00],%%rsi	\n\t"\
		"movslq	%[__p02],%%rcx							\n\t	addq	$0x1200,%%r8		\n\t"/* two */\
		"movslq	%[__p03],%%rdx							\n\t  vmovaps 0x1000(%%rsi),%%zmm29	\n\t"/* isrt2 */\
		"leaq	(%%rax,%%rbx,8),%%rbx					\n\t	vmovaps	-0x40(%%r8),%%zmm30	\n\t"/* 1.0 */\
		"leaq	(%%rax,%%rcx,8),%%rcx					\n\t	vmovaps	     (%%r8),%%zmm31	\n\t"/* 2.0 */\
		"leaq	(%%rax,%%rdx,8),%%rdx					\n\t	vmovaps	 0x40(%%r8),%%zmm28	\n\t"/* sqrt2 */\
		/* SSE2_RADIX4_DIT_0TWIDDLE_B(r00) */				/* SSE2_RADIX4_DIT_0TWIDDLE_2NDOFTWO_B(r08) */\
		"vmovaps	     (%%rax),%%zmm2					\n\t	vmovaps	     (%%rax,%%r9,8),%%zmm10		\n\t"\
		"vmovaps	     (%%rbx),%%zmm0					\n\t	vmovaps	     (%%rbx,%%r9,8),%%zmm8 		\n\t"\
		"vmovaps	     (%%rcx),%%zmm6					\n\t	vmovaps	     (%%rcx,%%r9,8),%%zmm14		\n\t"\
		"vmovaps	     (%%rdx),%%zmm4					\n\t	vmovaps	     (%%rdx,%%r9,8),%%zmm12		\n\t"\
		"vmovaps	0x040(%%rax),%%zmm3					\n\t	vmovaps	0x040(%%rax,%%r9,8),%%zmm11		\n\t"\
		"vmovaps	0x040(%%rbx),%%zmm1					\n\t	vmovaps	0x040(%%rbx,%%r9,8),%%zmm9 		\n\t"\
		"vmovaps	0x040(%%rcx),%%zmm7					\n\t	vmovaps	0x040(%%rcx,%%r9,8),%%zmm15		\n\t"\
		"vmovaps	0x040(%%rdx),%%zmm5					\n\t	vmovaps	0x040(%%rdx,%%r9,8),%%zmm13		\n\t"\
		"vfmsub132pd	%%zmm30,%%zmm0,%%zmm2			\n\t	 vfmsub132pd	%%zmm30,%%zmm8 ,%%zmm10	\n\t"\
		"vfmsub132pd	%%zmm30,%%zmm4,%%zmm6			\n\t	 vfmsub132pd	%%zmm30,%%zmm12,%%zmm14	\n\t"\
		"vfmsub132pd	%%zmm30,%%zmm1,%%zmm3			\n\t	 vfmsub132pd	%%zmm30,%%zmm9 ,%%zmm11	\n\t"\
		"vfmsub132pd	%%zmm30,%%zmm5,%%zmm7			\n\t	 vfmsub132pd	%%zmm30,%%zmm13,%%zmm15	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm2,%%zmm0			\n\t	 vfmadd132pd	%%zmm31,%%zmm10,%%zmm8 	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm6,%%zmm4			\n\t	 vfmadd132pd	%%zmm31,%%zmm14,%%zmm12	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm3,%%zmm1			\n\t	 vfmadd132pd	%%zmm31,%%zmm11,%%zmm9 	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm7,%%zmm5			\n\t	 vfmadd132pd	%%zmm31,%%zmm15,%%zmm13	\n\t"\
		/* Finish radix-4 butterfly: */\
		"vfmsub132pd	%%zmm30,%%zmm4,%%zmm0			\n\t	 vfmsub132pd	%%zmm30,%%zmm12,%%zmm8 	\n\t"\
		"vfmsub132pd	%%zmm30,%%zmm7,%%zmm2			\n\t	 vfmsub132pd	%%zmm30,%%zmm13,%%zmm9 	\n\t"\
		"vfmsub132pd	%%zmm30,%%zmm5,%%zmm1			\n\t	 vfmsub132pd	%%zmm30,%%zmm15,%%zmm10	\n\t"\
		"vfmsub132pd	%%zmm30,%%zmm6,%%zmm3			\n\t	 vfmsub132pd	%%zmm30,%%zmm14,%%zmm11	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm0,%%zmm4			\n\t	 vfmadd132pd	%%zmm31,%%zmm8 ,%%zmm12	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm2,%%zmm7			\n\t	 vfmadd132pd	%%zmm31,%%zmm9 ,%%zmm13	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm1,%%zmm5			\n\t	 vfmadd132pd	%%zmm31,%%zmm10,%%zmm15	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm3,%%zmm6			\n\t	 vfmadd132pd	%%zmm31,%%zmm11,%%zmm14	\n\t"\
		"vfmsub132pd	%%zmm30,%%zmm12,%%zmm4			\n\t	 vsubpd			%%zmm15,%%zmm11,%%zmm11	\n\t"\
		"vfmsub132pd	%%zmm30,%%zmm9 ,%%zmm0			\n\t	 vsubpd			%%zmm10,%%zmm14,%%zmm14	\n\t"\
		"vsubpd			%%zmm13,%%zmm5 ,%%zmm5			\n\t	 vfmadd132pd	%%zmm31,%%zmm11,%%zmm15	\n\t"\
		"vsubpd			%%zmm8 ,%%zmm1 ,%%zmm1			\n\t	 vfmadd132pd	%%zmm31,%%zmm14,%%zmm10	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm4 ,%%zmm12			\n\t	vfnmadd231pd	%%zmm29,%%zmm11,%%zmm3	\n\t"/* isrt2 */\
		"vfmadd132pd	%%zmm31,%%zmm0 ,%%zmm9 			\n\t	vfnmadd231pd	%%zmm29,%%zmm14,%%zmm2	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm5 ,%%zmm13			\n\t	vfnmadd231pd	%%zmm29,%%zmm15,%%zmm7	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm1 ,%%zmm8 			\n\t	vfnmadd231pd	%%zmm29,%%zmm10,%%zmm6	\n\t"\
		/* SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS(r00,r02,r04,r06,r08,r0A,r0C,r0E) */\
		"vmovaps	%%zmm4 ,0x200(%%rsi)				\n\t	vfmadd132pd	%%zmm28,%%zmm3 ,%%zmm11	\n\t"/* sqrt2, use of FMA above means we must combine 2*isrt2 here */\
		"vmovaps	%%zmm0 ,0x300(%%rsi)				\n\t	vfmadd132pd	%%zmm28,%%zmm2 ,%%zmm14	\n\t"\
		"vmovaps	%%zmm5 ,0x240(%%rsi)				\n\t	vfmadd132pd	%%zmm28,%%zmm7 ,%%zmm15	\n\t"\
		"vmovaps	%%zmm1 ,0x140(%%rsi)				\n\t	vfmadd132pd	%%zmm28,%%zmm6 ,%%zmm10	\n\t"\
		"vmovaps	%%zmm12,     (%%rsi)				\n\t	vmovaps	%%zmm7 ,0x280(%%rsi)	\n\t"\
		"vmovaps	%%zmm9 ,0x100(%%rsi)				\n\t	vmovaps	%%zmm2 ,0x380(%%rsi)	\n\t"\
		"vmovaps	%%zmm13,0x040(%%rsi)				\n\t	vmovaps	%%zmm3 ,0x2c0(%%rsi)	\n\t"\
		"vmovaps	%%zmm8 ,0x340(%%rsi)				\n\t	vmovaps	%%zmm6 ,0x1c0(%%rsi)	\n\t"\
		"														vmovaps	%%zmm15,0x080(%%rsi)	\n\t"\
		"														vmovaps	%%zmm14,0x180(%%rsi)	\n\t"\
		"														vmovaps	%%zmm11,0x0c0(%%rsi)	\n\t"\
		"														vmovaps	%%zmm10,0x3c0(%%rsi)	\n\t"\
	/*...Block 2: */\
		"addq	$0x400,%%rsi							\n\t"/* r10 */\
		"movslq	%[__p08],%%rdi							\n\t	addq	%%rdi,%%r9	\n\t"\
		/* SSE2_RADIX4_DIT_0TWIDDLE_B(r10) */			/* SSE2_RADIX4_DIT_0TWIDDLE_2NDOFTWO_B(r18) */\
		"vmovaps	     (%%rax,%%rdi,8),%%zmm2			\n\t	vmovaps	     (%%rax,%%r9,8),%%zmm10		\n\t"\
		"vmovaps	     (%%rbx,%%rdi,8),%%zmm0			\n\t	vmovaps	     (%%rbx,%%r9,8),%%zmm8 		\n\t"\
		"vmovaps	     (%%rcx,%%rdi,8),%%zmm6			\n\t	vmovaps	     (%%rcx,%%r9,8),%%zmm14		\n\t"\
		"vmovaps	     (%%rdx,%%rdi,8),%%zmm4			\n\t	vmovaps	     (%%rdx,%%r9,8),%%zmm12		\n\t"\
		"vmovaps	0x040(%%rax,%%rdi,8),%%zmm3			\n\t	vmovaps	0x040(%%rax,%%r9,8),%%zmm11		\n\t"\
		"vmovaps	0x040(%%rbx,%%rdi,8),%%zmm1			\n\t	vmovaps	0x040(%%rbx,%%r9,8),%%zmm9 		\n\t"\
		"vmovaps	0x040(%%rcx,%%rdi,8),%%zmm7			\n\t	vmovaps	0x040(%%rcx,%%r9,8),%%zmm15		\n\t"\
		"vmovaps	0x040(%%rdx,%%rdi,8),%%zmm5			\n\t	vmovaps	0x040(%%rdx,%%r9,8),%%zmm13		\n\t"\
		"vfmsub132pd	%%zmm30,%%zmm0,%%zmm2			\n\t	 vfmsub132pd	%%zmm30,%%zmm8 ,%%zmm10	\n\t"\
		"vfmsub132pd	%%zmm30,%%zmm4,%%zmm6			\n\t	 vfmsub132pd	%%zmm30,%%zmm12,%%zmm14	\n\t"\
		"vfmsub132pd	%%zmm30,%%zmm1,%%zmm3			\n\t	 vfmsub132pd	%%zmm30,%%zmm9 ,%%zmm11	\n\t"\
		"vfmsub132pd	%%zmm30,%%zmm5,%%zmm7			\n\t	 vfmsub132pd	%%zmm30,%%zmm13,%%zmm15	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm2,%%zmm0			\n\t	 vfmadd132pd	%%zmm31,%%zmm10,%%zmm8 	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm6,%%zmm4			\n\t	 vfmadd132pd	%%zmm31,%%zmm14,%%zmm12	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm3,%%zmm1			\n\t	 vfmadd132pd	%%zmm31,%%zmm11,%%zmm9 	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm7,%%zmm5			\n\t	 vfmadd132pd	%%zmm31,%%zmm15,%%zmm13	\n\t"\
		/* Finish radix-4 butterfly: */\
		"vfmsub132pd	%%zmm30,%%zmm4,%%zmm0			\n\t	 vfmsub132pd	%%zmm30,%%zmm12,%%zmm8 	\n\t"\
		"vfmsub132pd	%%zmm30,%%zmm7,%%zmm2			\n\t	 vfmsub132pd	%%zmm30,%%zmm13,%%zmm9 	\n\t"\
		"vfmsub132pd	%%zmm30,%%zmm5,%%zmm1			\n\t	 vfmsub132pd	%%zmm30,%%zmm15,%%zmm10	\n\t"\
		"vfmsub132pd	%%zmm30,%%zmm6,%%zmm3			\n\t	 vfmsub132pd	%%zmm30,%%zmm14,%%zmm11	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm0,%%zmm4			\n\t	 vfmadd132pd	%%zmm31,%%zmm8 ,%%zmm12	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm2,%%zmm7			\n\t	 vfmadd132pd	%%zmm31,%%zmm9 ,%%zmm13	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm1,%%zmm5			\n\t	 vfmadd132pd	%%zmm31,%%zmm10,%%zmm15	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm3,%%zmm6			\n\t	 vfmadd132pd	%%zmm31,%%zmm11,%%zmm14	\n\t"\
		"vfmsub132pd	%%zmm30,%%zmm12,%%zmm4			\n\t	 vsubpd			%%zmm15,%%zmm11,%%zmm11	\n\t"\
		"vfmsub132pd	%%zmm30,%%zmm9 ,%%zmm0			\n\t	 vsubpd			%%zmm10,%%zmm14,%%zmm14	\n\t"\
		"vsubpd			%%zmm13,%%zmm5 ,%%zmm5			\n\t	 vfmadd132pd	%%zmm31,%%zmm11,%%zmm15	\n\t"\
		"vsubpd			%%zmm8 ,%%zmm1 ,%%zmm1			\n\t	 vfmadd132pd	%%zmm31,%%zmm14,%%zmm10	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm4 ,%%zmm12			\n\t	vfnmadd231pd	%%zmm29,%%zmm11,%%zmm3	\n\t"/* isrt2 */\
		"vfmadd132pd	%%zmm31,%%zmm0 ,%%zmm9 			\n\t	vfnmadd231pd	%%zmm29,%%zmm14,%%zmm2	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm5 ,%%zmm13			\n\t	vfnmadd231pd	%%zmm29,%%zmm15,%%zmm7	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm1 ,%%zmm8 			\n\t	vfnmadd231pd	%%zmm29,%%zmm10,%%zmm6	\n\t"\
		/* SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS(r10,r12,r14,r16,r18,r1A,r1C,r1E) */\
		"vmovaps	%%zmm4 ,0x200(%%rsi)				\n\t	vfmadd132pd	%%zmm28,%%zmm3 ,%%zmm11	\n\t"/* sqrt2, use of FMA above means we must combine 2*isrt2 here */\
		"vmovaps	%%zmm0 ,0x300(%%rsi)				\n\t	vfmadd132pd	%%zmm28,%%zmm2 ,%%zmm14	\n\t"\
		"vmovaps	%%zmm5 ,0x240(%%rsi)				\n\t	vfmadd132pd	%%zmm28,%%zmm7 ,%%zmm15	\n\t"\
		"vmovaps	%%zmm1 ,0x140(%%rsi)				\n\t	vfmadd132pd	%%zmm28,%%zmm6 ,%%zmm10	\n\t"\
		"vmovaps	%%zmm12,     (%%rsi)				\n\t	vmovaps	%%zmm7 ,0x280(%%rsi)	\n\t"\
		"vmovaps	%%zmm9 ,0x100(%%rsi)				\n\t	vmovaps	%%zmm2 ,0x380(%%rsi)	\n\t"\
		"vmovaps	%%zmm13,0x040(%%rsi)				\n\t	vmovaps	%%zmm3 ,0x2c0(%%rsi)	\n\t"\
		"vmovaps	%%zmm8 ,0x340(%%rsi)				\n\t	vmovaps	%%zmm6 ,0x1c0(%%rsi)	\n\t"\
		"														vmovaps	%%zmm15,0x080(%%rsi)	\n\t"\
		"														vmovaps	%%zmm14,0x180(%%rsi)	\n\t"\
		"														vmovaps	%%zmm11,0x0c0(%%rsi)	\n\t"\
		"														vmovaps	%%zmm10,0x3c0(%%rsi)	\n\t"\
		"\n\t"\
	/*...Block 3: */\
		"addq	$0x400,%%rsi /* r20 */					\n\t	subq	%%rdi,%%r9	\n\t"\
		"movslq	%[__p10],%%rdi							\n\t	addq	%%rdi,%%r9	\n\t"\
		/* SSE2_RADIX4_DIT_0TWIDDLE_B(r20) */			/* SSE2_RADIX4_DIT_0TWIDDLE_2NDOFTWO_B(r28) */\
		"vmovaps	     (%%rax,%%rdi,8),%%zmm2			\n\t	vmovaps	     (%%rax,%%r9,8),%%zmm10		\n\t"\
		"vmovaps	     (%%rbx,%%rdi,8),%%zmm0			\n\t	vmovaps	     (%%rbx,%%r9,8),%%zmm8 		\n\t"\
		"vmovaps	     (%%rcx,%%rdi,8),%%zmm6			\n\t	vmovaps	     (%%rcx,%%r9,8),%%zmm14		\n\t"\
		"vmovaps	     (%%rdx,%%rdi,8),%%zmm4			\n\t	vmovaps	     (%%rdx,%%r9,8),%%zmm12		\n\t"\
		"vmovaps	0x040(%%rax,%%rdi,8),%%zmm3			\n\t	vmovaps	0x040(%%rax,%%r9,8),%%zmm11		\n\t"\
		"vmovaps	0x040(%%rbx,%%rdi,8),%%zmm1			\n\t	vmovaps	0x040(%%rbx,%%r9,8),%%zmm9 		\n\t"\
		"vmovaps	0x040(%%rcx,%%rdi,8),%%zmm7			\n\t	vmovaps	0x040(%%rcx,%%r9,8),%%zmm15		\n\t"\
		"vmovaps	0x040(%%rdx,%%rdi,8),%%zmm5			\n\t	vmovaps	0x040(%%rdx,%%r9,8),%%zmm13		\n\t"\
		"vfmsub132pd	%%zmm30,%%zmm0,%%zmm2			\n\t	 vfmsub132pd	%%zmm30,%%zmm8 ,%%zmm10	\n\t"\
		"vfmsub132pd	%%zmm30,%%zmm4,%%zmm6			\n\t	 vfmsub132pd	%%zmm30,%%zmm12,%%zmm14	\n\t"\
		"vfmsub132pd	%%zmm30,%%zmm1,%%zmm3			\n\t	 vfmsub132pd	%%zmm30,%%zmm9 ,%%zmm11	\n\t"\
		"vfmsub132pd	%%zmm30,%%zmm5,%%zmm7			\n\t	 vfmsub132pd	%%zmm30,%%zmm13,%%zmm15	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm2,%%zmm0			\n\t	 vfmadd132pd	%%zmm31,%%zmm10,%%zmm8 	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm6,%%zmm4			\n\t	 vfmadd132pd	%%zmm31,%%zmm14,%%zmm12	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm3,%%zmm1			\n\t	 vfmadd132pd	%%zmm31,%%zmm11,%%zmm9 	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm7,%%zmm5			\n\t	 vfmadd132pd	%%zmm31,%%zmm15,%%zmm13	\n\t"\
		/* Finish radix-4 butterfly: */\
		"vfmsub132pd	%%zmm30,%%zmm4,%%zmm0			\n\t	 vfmsub132pd	%%zmm30,%%zmm12,%%zmm8 	\n\t"\
		"vfmsub132pd	%%zmm30,%%zmm7,%%zmm2			\n\t	 vfmsub132pd	%%zmm30,%%zmm13,%%zmm9 	\n\t"\
		"vfmsub132pd	%%zmm30,%%zmm5,%%zmm1			\n\t	 vfmsub132pd	%%zmm30,%%zmm15,%%zmm10	\n\t"\
		"vfmsub132pd	%%zmm30,%%zmm6,%%zmm3			\n\t	 vfmsub132pd	%%zmm30,%%zmm14,%%zmm11	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm0,%%zmm4			\n\t	 vfmadd132pd	%%zmm31,%%zmm8 ,%%zmm12	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm2,%%zmm7			\n\t	 vfmadd132pd	%%zmm31,%%zmm9 ,%%zmm13	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm1,%%zmm5			\n\t	 vfmadd132pd	%%zmm31,%%zmm10,%%zmm15	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm3,%%zmm6			\n\t	 vfmadd132pd	%%zmm31,%%zmm11,%%zmm14	\n\t"\
		"vfmsub132pd	%%zmm30,%%zmm12,%%zmm4			\n\t	 vsubpd			%%zmm15,%%zmm11,%%zmm11	\n\t"\
		"vfmsub132pd	%%zmm30,%%zmm9 ,%%zmm0			\n\t	 vsubpd			%%zmm10,%%zmm14,%%zmm14	\n\t"\
		"vsubpd			%%zmm13,%%zmm5 ,%%zmm5			\n\t	 vfmadd132pd	%%zmm31,%%zmm11,%%zmm15	\n\t"\
		"vsubpd			%%zmm8 ,%%zmm1 ,%%zmm1			\n\t	 vfmadd132pd	%%zmm31,%%zmm14,%%zmm10	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm4 ,%%zmm12			\n\t	vfnmadd231pd	%%zmm29,%%zmm11,%%zmm3	\n\t"/* isrt2 */\
		"vfmadd132pd	%%zmm31,%%zmm0 ,%%zmm9 			\n\t	vfnmadd231pd	%%zmm29,%%zmm14,%%zmm2	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm5 ,%%zmm13			\n\t	vfnmadd231pd	%%zmm29,%%zmm15,%%zmm7	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm1 ,%%zmm8 			\n\t	vfnmadd231pd	%%zmm29,%%zmm10,%%zmm6	\n\t"\
		/* SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS(r20,r22,r24,r26,r28,r2A,r2C,r2E) */\
		"vmovaps	%%zmm4 ,0x200(%%rsi)				\n\t	vfmadd132pd	%%zmm28,%%zmm3 ,%%zmm11	\n\t"/* sqrt2, use of FMA above means we must combine 2*isrt2 here */\
		"vmovaps	%%zmm0 ,0x300(%%rsi)				\n\t	vfmadd132pd	%%zmm28,%%zmm2 ,%%zmm14	\n\t"\
		"vmovaps	%%zmm5 ,0x240(%%rsi)				\n\t	vfmadd132pd	%%zmm28,%%zmm7 ,%%zmm15	\n\t"\
		"vmovaps	%%zmm1 ,0x140(%%rsi)				\n\t	vfmadd132pd	%%zmm28,%%zmm6 ,%%zmm10	\n\t"\
		"vmovaps	%%zmm12,     (%%rsi)				\n\t	vmovaps	%%zmm7 ,0x280(%%rsi)	\n\t"\
		"vmovaps	%%zmm9 ,0x100(%%rsi)				\n\t	vmovaps	%%zmm2 ,0x380(%%rsi)	\n\t"\
		"vmovaps	%%zmm13,0x040(%%rsi)				\n\t	vmovaps	%%zmm3 ,0x2c0(%%rsi)	\n\t"\
		"vmovaps	%%zmm8 ,0x340(%%rsi)				\n\t	vmovaps	%%zmm6 ,0x1c0(%%rsi)	\n\t"\
		"														vmovaps	%%zmm15,0x080(%%rsi)	\n\t"\
		"														vmovaps	%%zmm14,0x180(%%rsi)	\n\t"\
		"														vmovaps	%%zmm11,0x0c0(%%rsi)	\n\t"\
		"														vmovaps	%%zmm10,0x3c0(%%rsi)	\n\t"\
		"\n\t"\
	/*...Block 4: */\
		"addq	$0x400,%%rsi /* r30 */					\n\t	subq	%%rdi,%%r9	\n\t"\
		"movslq	%[__p08],%%r14					\n\t"\
		"movslq	%[__p10],%%rdi					\n\t"\
		"addq	%%r14,%%rdi		/* p18 */				\n\t	addq	%%rdi,%%r9	\n\t"\
		/* SSE2_RADIX4_DIT_0TWIDDLE_B(r30) */		/* SSE2_RADIX4_DIT_0TWIDDLE_2NDOFTWO_B(r38) */\
		"vmovaps	     (%%rax,%%rdi,8),%%zmm2			\n\t	vmovaps	     (%%rax,%%r9,8),%%zmm10		\n\t"\
		"vmovaps	     (%%rbx,%%rdi,8),%%zmm0			\n\t	vmovaps	     (%%rbx,%%r9,8),%%zmm8 		\n\t"\
		"vmovaps	     (%%rcx,%%rdi,8),%%zmm6			\n\t	vmovaps	     (%%rcx,%%r9,8),%%zmm14		\n\t"\
		"vmovaps	     (%%rdx,%%rdi,8),%%zmm4			\n\t	vmovaps	     (%%rdx,%%r9,8),%%zmm12		\n\t"\
		"vmovaps	0x040(%%rax,%%rdi,8),%%zmm3			\n\t	vmovaps	0x040(%%rax,%%r9,8),%%zmm11		\n\t"\
		"vmovaps	0x040(%%rbx,%%rdi,8),%%zmm1			\n\t	vmovaps	0x040(%%rbx,%%r9,8),%%zmm9 		\n\t"\
		"vmovaps	0x040(%%rcx,%%rdi,8),%%zmm7			\n\t	vmovaps	0x040(%%rcx,%%r9,8),%%zmm15		\n\t"\
		"vmovaps	0x040(%%rdx,%%rdi,8),%%zmm5			\n\t	vmovaps	0x040(%%rdx,%%r9,8),%%zmm13		\n\t"\
		"vfmsub132pd	%%zmm30,%%zmm0,%%zmm2			\n\t	 vfmsub132pd	%%zmm30,%%zmm8 ,%%zmm10	\n\t"\
		"vfmsub132pd	%%zmm30,%%zmm4,%%zmm6			\n\t	 vfmsub132pd	%%zmm30,%%zmm12,%%zmm14	\n\t"\
		"vfmsub132pd	%%zmm30,%%zmm1,%%zmm3			\n\t	 vfmsub132pd	%%zmm30,%%zmm9 ,%%zmm11	\n\t"\
		"vfmsub132pd	%%zmm30,%%zmm5,%%zmm7			\n\t	 vfmsub132pd	%%zmm30,%%zmm13,%%zmm15	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm2,%%zmm0			\n\t	 vfmadd132pd	%%zmm31,%%zmm10,%%zmm8 	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm6,%%zmm4			\n\t	 vfmadd132pd	%%zmm31,%%zmm14,%%zmm12	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm3,%%zmm1			\n\t	 vfmadd132pd	%%zmm31,%%zmm11,%%zmm9 	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm7,%%zmm5			\n\t	 vfmadd132pd	%%zmm31,%%zmm15,%%zmm13	\n\t"\
		/* Finish radix-4 butterfly: */\
		"vfmsub132pd	%%zmm30,%%zmm4,%%zmm0			\n\t	 vfmsub132pd	%%zmm30,%%zmm12,%%zmm8 	\n\t"\
		"vfmsub132pd	%%zmm30,%%zmm7,%%zmm2			\n\t	 vfmsub132pd	%%zmm30,%%zmm13,%%zmm9 	\n\t"\
		"vfmsub132pd	%%zmm30,%%zmm5,%%zmm1			\n\t	 vfmsub132pd	%%zmm30,%%zmm15,%%zmm10	\n\t"\
		"vfmsub132pd	%%zmm30,%%zmm6,%%zmm3			\n\t	 vfmsub132pd	%%zmm30,%%zmm14,%%zmm11	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm0,%%zmm4			\n\t	 vfmadd132pd	%%zmm31,%%zmm8 ,%%zmm12	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm2,%%zmm7			\n\t	 vfmadd132pd	%%zmm31,%%zmm9 ,%%zmm13	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm1,%%zmm5			\n\t	 vfmadd132pd	%%zmm31,%%zmm10,%%zmm15	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm3,%%zmm6			\n\t	 vfmadd132pd	%%zmm31,%%zmm11,%%zmm14	\n\t"\
		"vfmsub132pd	%%zmm30,%%zmm12,%%zmm4			\n\t	 vsubpd			%%zmm15,%%zmm11,%%zmm11	\n\t"\
		"vfmsub132pd	%%zmm30,%%zmm9 ,%%zmm0			\n\t	 vsubpd			%%zmm10,%%zmm14,%%zmm14	\n\t"\
		"vsubpd			%%zmm13,%%zmm5 ,%%zmm5			\n\t	 vfmadd132pd	%%zmm31,%%zmm11,%%zmm15	\n\t"\
		"vsubpd			%%zmm8 ,%%zmm1 ,%%zmm1			\n\t	 vfmadd132pd	%%zmm31,%%zmm14,%%zmm10	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm4 ,%%zmm12			\n\t	vfnmadd231pd	%%zmm29,%%zmm11,%%zmm3	\n\t"/* isrt2 */\
		"vfmadd132pd	%%zmm31,%%zmm0 ,%%zmm9 			\n\t	vfnmadd231pd	%%zmm29,%%zmm14,%%zmm2	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm5 ,%%zmm13			\n\t	vfnmadd231pd	%%zmm29,%%zmm15,%%zmm7	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm1 ,%%zmm8 			\n\t	vfnmadd231pd	%%zmm29,%%zmm10,%%zmm6	\n\t"\
		/* SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS(r30,r32,r34,r36,r38,r3A,r3C,r3E) */\
		"vmovaps	%%zmm4 ,0x200(%%rsi)				\n\t	vfmadd132pd	%%zmm28,%%zmm3 ,%%zmm11	\n\t"/* sqrt2, use of FMA above means we must combine 2*isrt2 here */\
		"vmovaps	%%zmm0 ,0x300(%%rsi)				\n\t	vfmadd132pd	%%zmm28,%%zmm2 ,%%zmm14	\n\t"\
		"vmovaps	%%zmm5 ,0x240(%%rsi)				\n\t	vfmadd132pd	%%zmm28,%%zmm7 ,%%zmm15	\n\t"\
		"vmovaps	%%zmm1 ,0x140(%%rsi)				\n\t	vfmadd132pd	%%zmm28,%%zmm6 ,%%zmm10	\n\t"\
		"vmovaps	%%zmm12,     (%%rsi)				\n\t	vmovaps	%%zmm7 ,0x280(%%rsi)	\n\t"\
		"vmovaps	%%zmm9 ,0x100(%%rsi)				\n\t	vmovaps	%%zmm2 ,0x380(%%rsi)	\n\t"\
		"vmovaps	%%zmm13,0x040(%%rsi)				\n\t	vmovaps	%%zmm3 ,0x2c0(%%rsi)	\n\t"\
		"vmovaps	%%zmm8 ,0x340(%%rsi)				\n\t	vmovaps	%%zmm6 ,0x1c0(%%rsi)	\n\t"\
		"														vmovaps	%%zmm15,0x080(%%rsi)	\n\t"\
		"														vmovaps	%%zmm14,0x180(%%rsi)	\n\t"\
		"														vmovaps	%%zmm11,0x0c0(%%rsi)	\n\t"\
		"														vmovaps	%%zmm10,0x3c0(%%rsi)	\n\t"\
		"\n\t"\
	/**********************************************************************************/\
	/*...and now do eight radix-4 transforms, including the internal twiddle factors: */\
	/**********************************************************************************/\
		"\n\t"\
	/*...Block 1: t00,t10,t20,t30	*/							/*...Block 5: t08,t18,t28,t38*/\
		"movq	%[__r00],%%r8	\n\t"/* r10,20,30 += 1,2,3*0x400 ... inputs to our octet of radix-4 pass 2 DFTs from this local-store */\
		"movq	%[__add0],%%rax					\n\t	movslq	%[__p04],%%r9				\n\t"\
		"movslq	%[__p08],%%rbx					\n\t	movq	%[__isrt2],%%rsi			\n\t"\
		"movslq	%[__p10],%%rcx					\n\t	vmovaps	0xa40(%%r8),%%zmm13			\n\t"\
		"leaq	(%%rbx,%%rcx),%%rdx	/* p18 */	\n\t	vmovaps	0xe00(%%r8),%%zmm14			\n\t"\
		"leaq	(%%rax,%%rbx,8),%%rbx			\n\t"/* add0 + p08 */\
		"leaq	(%%rax,%%rcx,8),%%rcx			\n\t"/* add0 + p10 */\
		"leaq	(%%rax,%%rdx,8),%%rdx			\n\t"/* add0 + p18 */\
		"vmovaps	0x400(%%r8),%%zmm2			\n\t	vmulpd	%%zmm29,%%zmm13,%%zmm13		\n\t"\
		"vmovaps	     (%%r8),%%zmm0			\n\t	vmulpd	%%zmm29,%%zmm14,%%zmm14		\n\t"\
		"vmovaps	0xc00(%%r8),%%zmm4			\n\t	vmovaps	0xa00(%%r8),%%zmm12			\n\t"\
		"vmovaps	0x800(%%r8),%%zmm6			\n\t	vmovaps	0xe40(%%r8),%%zmm15			\n\t"\
		"vmovaps	0x440(%%r8),%%zmm3			\n\t	vmovaps	0x200(%%r8),%%zmm8 			\n\t"\
		"vmovaps	0x040(%%r8),%%zmm1			\n\t	vmovaps	0x640(%%r8),%%zmm10			\n\t"\
		"vmovaps	0xc40(%%r8),%%zmm5			\n\t	vmovaps	0x240(%%r8),%%zmm11			\n\t"\
		"vmovaps	0x840(%%r8),%%zmm7			\n\t	vmovaps	0x600(%%r8),%%zmm9 			\n\t"\
	"vfmsub132pd	%%zmm30,%%zmm2,%%zmm0		\n\t	vfnmadd231pd	%%zmm12,%%zmm29,%%zmm13	\n\t"\
	"vfmsub132pd	%%zmm30,%%zmm4,%%zmm6		\n\t	vfnmadd231pd	%%zmm15,%%zmm29,%%zmm14	\n\t"\
	"vfmsub132pd	%%zmm30,%%zmm3,%%zmm1		\n\t	 vfmadd132pd	%%zmm28,%%zmm13,%%zmm12	\n\t"/* sqrt2 */\
	"vfmsub132pd	%%zmm30,%%zmm5,%%zmm7		\n\t	 vfmadd132pd	%%zmm28,%%zmm14,%%zmm15	\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm0,%%zmm2		\n\t	vfmsub132pd		%%zmm30,%%zmm10,%%zmm8 	\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm6,%%zmm4		\n\t	vfmsub132pd		%%zmm30,%%zmm9 ,%%zmm11	\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm1,%%zmm3		\n\t	vfmsub132pd		%%zmm30,%%zmm14,%%zmm12	\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm7,%%zmm5		\n\t	vfmsub132pd		%%zmm30,%%zmm15,%%zmm13	\n\t"\
		/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c00): */\
	"vfmsub132pd	%%zmm30,%%zmm4,%%zmm2		\n\t	vfmadd132pd	%%zmm31,%%zmm8 ,%%zmm10		\n\t"\
	"vfmsub132pd	%%zmm30,%%zmm7,%%zmm0		\n\t	vfmadd132pd	%%zmm31,%%zmm11,%%zmm9 		\n\t"\
	"vfmsub132pd	%%zmm30,%%zmm5,%%zmm3		\n\t	vfmadd132pd	%%zmm31,%%zmm12,%%zmm14		\n\t"\
	"vsubpd			%%zmm6 ,%%zmm1,%%zmm1		\n\t	vfmadd132pd	%%zmm31,%%zmm13,%%zmm15		\n\t"\
		"addq	$0x1c0,%%rsi	\n\t"/* c00 */		/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c04 = c00 + 0x200): */\
	"vfmadd132pd	%%zmm31,%%zmm2,%%zmm4		\n\t	vfmsub132pd	%%zmm30,%%zmm13,%%zmm11		\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm0,%%zmm7		\n\t	vfmsub132pd	%%zmm30,%%zmm12,%%zmm10		\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm3,%%zmm5		\n\t	vfmsub132pd	%%zmm30,%%zmm14,%%zmm9 		\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm1,%%zmm6		\n\t	vfmsub132pd	%%zmm30,%%zmm15,%%zmm8 		\n\t"\
		"vmovaps	%%zmm2,     (%%r8)			\n\t	vfmadd132pd	%%zmm31,%%zmm11,%%zmm13		\n\t"\
		"vmovaps	%%zmm0,0x440(%%r8)			\n\t	vfmadd132pd	%%zmm31,%%zmm10,%%zmm12		\n\t"\
		"vmovaps	%%zmm3,0x040(%%r8)			\n\t	vfmadd132pd	%%zmm31,%%zmm9 ,%%zmm14		\n\t"\
		"vmovaps	%%zmm6,0x400(%%r8)			\n\t	vfmadd132pd	%%zmm31,%%zmm8 ,%%zmm15		\n\t"\
		"vmovaps		  (%%rsi),%%zmm16		\n\t	vmovaps		0x200(%%rsi),%%zmm20		\n\t"\
		"vmovaps	 0x080(%%rsi),%%zmm17		\n\t	vmovaps		0x280(%%rsi),%%zmm21		\n\t"\
		"vmovaps	 0x040(%%rsi),%%zmm18		\n\t	vmovaps		0x240(%%rsi),%%zmm22		\n\t"\
		"vmovaps	 0x0c0(%%rsi),%%zmm19		\n\t	vmovaps		0x2c0(%%rsi),%%zmm23		\n\t"\
		"vmulpd		%%zmm16,%%zmm5,%%zmm3		\n\t		vmulpd		%%zmm20,%%zmm13,%%zmm24	\n\t"\
		"vmulpd		%%zmm17,%%zmm1,%%zmm6		\n\t		vmulpd		%%zmm21,%%zmm9 ,%%zmm25	\n\t"\
		"vmulpd		%%zmm16,%%zmm4,%%zmm2		\n\t		vmulpd		%%zmm20,%%zmm12,%%zmm26	\n\t"\
		"vmulpd		%%zmm17,%%zmm7,%%zmm0		\n\t		vmulpd		%%zmm21,%%zmm15,%%zmm27	\n\t"\
	"vfnmadd231pd	%%zmm18,%%zmm4,%%zmm3		\n\t	vfnmadd231pd	%%zmm22,%%zmm12,%%zmm24	\n\t"\
	"vfnmadd231pd	%%zmm19,%%zmm7,%%zmm6		\n\t	vfnmadd231pd	%%zmm23,%%zmm15,%%zmm25	\n\t"\
	" vfmadd231pd	%%zmm18,%%zmm5,%%zmm2		\n\t	 vfmadd231pd	%%zmm22,%%zmm13,%%zmm26	\n\t"\
	" vfmadd231pd	%%zmm19,%%zmm1,%%zmm0		\n\t	 vfmadd231pd	%%zmm23,%%zmm9 ,%%zmm27	\n\t"\
		"vmovaps	%%zmm3,0x040(%%rax)			\n\t	vmovaps	%%zmm24,0x040(%%rax,%%r9,8)		\n\t"\
		"vmovaps	%%zmm6,0x040(%%rbx)			\n\t	vmovaps	%%zmm25,0x040(%%rbx,%%r9,8)		\n\t"\
		"vmovaps	%%zmm2,     (%%rax)			\n\t	vmovaps	%%zmm26,     (%%rax,%%r9,8)		\n\t"\
		"vmovaps	%%zmm0,     (%%rbx)			\n\t	vmovaps	%%zmm27,     (%%rbx,%%r9,8)		\n\t"\
		"vmovaps	0x040(%%r8),%%zmm5			\n\t"\
		"vmovaps	0x400(%%r8),%%zmm6			\n\t"\
		"vmovaps	     (%%r8),%%zmm4			\n\t"\
		"vmovaps	0x440(%%r8),%%zmm0			\n\t"\
		"addq	$0x100,%%rsi	\n\t"/* c10 */			/* c14 = c10 + 0x200: */\
		"vmovaps	     (%%rsi),%%zmm16		\n\t	vmovaps		0x200(%%rsi),%%zmm20		\n\t"\
		"vmovaps	0x080(%%rsi),%%zmm17		\n\t	vmovaps		0x280(%%rsi),%%zmm21		\n\t"\
		"vmovaps	0x040(%%rsi),%%zmm18		\n\t	vmovaps		0x240(%%rsi),%%zmm22		\n\t"\
		"vmovaps	0x0c0(%%rsi),%%zmm19		\n\t	vmovaps		0x2c0(%%rsi),%%zmm23		\n\t"\
		"vmulpd			%%zmm16,%%zmm5,%%zmm3	\n\t	vmulpd			%%zmm20,%%zmm11,%%zmm13	\n\t"\
		"vmulpd			%%zmm17,%%zmm6,%%zmm7	\n\t	vmulpd			%%zmm21,%%zmm14,%%zmm15	\n\t"\
		"vmulpd			%%zmm16,%%zmm4,%%zmm2	\n\t	vmulpd			%%zmm20,%%zmm10,%%zmm12	\n\t"\
		"vmulpd			%%zmm17,%%zmm0,%%zmm1	\n\t	vmulpd			%%zmm21,%%zmm8 ,%%zmm9 	\n\t"\
	"vfnmadd231pd		%%zmm18,%%zmm4,%%zmm3	\n\t	vfnmadd231pd	%%zmm22,%%zmm10,%%zmm13	\n\t"\
	"vfnmadd231pd		%%zmm19,%%zmm0,%%zmm7	\n\t	vfnmadd231pd	%%zmm23,%%zmm8 ,%%zmm15	\n\t"\
	" vfmadd231pd		%%zmm18,%%zmm5,%%zmm2	\n\t	 vfmadd231pd	%%zmm22,%%zmm11,%%zmm12	\n\t"\
	" vfmadd231pd		%%zmm19,%%zmm6,%%zmm1	\n\t	 vfmadd231pd	%%zmm23,%%zmm14,%%zmm9 	\n\t"\
		"vmovaps	%%zmm3,0x040(%%rcx)			\n\t	vmovaps	%%zmm13,0x040(%%rcx,%%r9,8)		\n\t"\
		"vmovaps	%%zmm7,0x040(%%rdx)			\n\t	vmovaps	%%zmm15,0x040(%%rdx,%%r9,8)		\n\t"\
		"vmovaps	%%zmm2,     (%%rcx)			\n\t	vmovaps	%%zmm12,     (%%rcx,%%r9,8)		\n\t"\
		"vmovaps	%%zmm1,     (%%rdx)			\n\t	vmovaps	%%zmm9 ,     (%%rdx,%%r9,8)		\n\t"\
		"\n\t"\
	/*...Block 2: t02,t12,t22,t32	*/					/*...Block 6: t0A,t1A,t2A,t3A*/\
		"movslq	%[__p01],%%rdi					\n\t	addq	%%rdi,%%r9		\n\t"\
		"movq	%[__isrt2],%%rsi				\n\t"\
		"addq	$0x880,%%r8						\n\t"/* r22 */\
		"addq	$0x0c0,%%rsi					\n\t"/* cc1; cc3 += 0x080: */\
		"vmovaps	     (%%r8),%%zmm4			\n\t	vmovaps	0x200(%%r8),%%zmm12				\n\t"\
		"vmovaps	0x040(%%r8),%%zmm5			\n\t	vmovaps	0x240(%%r8),%%zmm13				\n\t"\
		"vmovaps	     (%%rsi),%%zmm2			\n\t	vmovaps	0x080(%%rsi),%%zmm11			\n\t"/* cc1,cc3 */\
		"vmovaps	0x040(%%rsi),%%zmm3			\n\t	vmovaps	0x0c0(%%rsi),%%zmm10			\n\t"\
		"vmovaps	%%zmm4,%%zmm6				\n\t	vmovaps	%%zmm12,%%zmm14					\n\t"\
		"vmovaps	%%zmm5,%%zmm7				\n\t	vmovaps	%%zmm13,%%zmm15					\n\t"\
		"vmulpd		%%zmm2,%%zmm4,%%zmm4		\n\t	vmulpd		%%zmm10,%%zmm12,%%zmm12		\n\t"\
		"vmulpd		%%zmm2,%%zmm5,%%zmm5		\n\t	vmulpd		%%zmm10,%%zmm13,%%zmm13		\n\t"\
		"vmovaps	0x400(%%r8),%%zmm0			\n\t	vmovaps	0x600(%%r8),%%zmm8 				\n\t"\
		"vmovaps	0x440(%%r8),%%zmm1			\n\t	vmovaps	0x640(%%r8),%%zmm9 	/* t3B */	\n\t"\
		" vfmadd231pd	%%zmm3,%%zmm7,%%zmm4	\n\t	 vfmadd231pd	%%zmm11,%%zmm15,%%zmm12	\n\t"\
		"vfnmadd231pd	%%zmm3,%%zmm6,%%zmm5	\n\t	vfnmadd231pd	%%zmm11,%%zmm14,%%zmm13	\n\t"\
		"vmovaps	%%zmm0,%%zmm6				\n\t	vmovaps	%%zmm8 ,%%zmm14					\n\t"\
		"vmovaps	%%zmm1,%%zmm7				\n\t	vmovaps	%%zmm9 ,%%zmm15					\n\t"\
		"vmulpd		%%zmm11,%%zmm0,%%zmm0		\n\t	vmulpd		%%zmm2 ,%%zmm8 ,%%zmm8 		\n\t"\
		"vmulpd		%%zmm11,%%zmm1,%%zmm1		\n\t	vmulpd		%%zmm2 ,%%zmm9 ,%%zmm9 		\n\t"\
		" vfmadd231pd	%%zmm10,%%zmm7,%%zmm0	\n\t	vfnmadd231pd	%%zmm3 ,%%zmm15,%%zmm8 	\n\t"\
		"vfnmadd231pd	%%zmm10,%%zmm6,%%zmm1	\n\t	 vfmadd231pd	%%zmm3 ,%%zmm14,%%zmm9 	\n\t"\
		"vmovaps	%%zmm5,%%zmm7				\n\t	vmovaps	%%zmm13,%%zmm15					\n\t"\
		"vmovaps	%%zmm4,%%zmm6				\n\t	vmovaps	%%zmm12,%%zmm14					\n\t"\
		"vfmadd132pd	%%zmm30,%%zmm0,%%zmm4	\n\t	vfmadd132pd	%%zmm30,%%zmm8 ,%%zmm14		\n\t"\
		"vfmadd132pd	%%zmm30,%%zmm1,%%zmm5	\n\t	vfmadd132pd	%%zmm30,%%zmm9 ,%%zmm15		\n\t"\
		"vfmsub132pd	%%zmm30,%%zmm0,%%zmm6	\n\t	vfmsub132pd	%%zmm30,%%zmm8 ,%%zmm12		\n\t"\
		"vfmsub132pd	%%zmm30,%%zmm1,%%zmm7	\n\t	vfmsub132pd	%%zmm30,%%zmm9 ,%%zmm13		\n\t"\
		"subq	$0x800,%%r8		/* r02 */		\n\t"\
		"subq	$0x080,%%rsi	/* cc0 */		\n\t	vmovaps	     (%%rsi),%%zmm28			\n\t"/* cc0 */\
		"vmovaps	0x400(%%r8),%%zmm1			\n\t	vmovaps	0x600(%%r8),%%zmm8 				\n\t"\
		"vmovaps	0x440(%%r8),%%zmm3			\n\t	vmovaps	0x640(%%r8),%%zmm10				\n\t"\
		"vmovaps	0x040(%%rsi),%%zmm2			\n\t	vmovaps	0x040(%%rsi),%%zmm9 			\n\t"\
		"vmovaps	%%zmm1,%%zmm0				\n\t	vmovaps	%%zmm8 ,%%zmm11					\n\t"\
		"vmulpd		%%zmm2 ,%%zmm1,%%zmm1		\n\t	vmulpd		%%zmm9 ,%%zmm8 ,%%zmm8 		\n\t"\
		"vmulpd		%%zmm3 ,%%zmm2,%%zmm2		\n\t	vmulpd		%%zmm10,%%zmm9 ,%%zmm9 		\n\t"\
		" vfmsub132pd	%%zmm28,%%zmm1,%%zmm3	\n\t	vfnmadd231pd	%%zmm28,%%zmm10,%%zmm8 	\n\t"/* cc0 */\
		" vfmadd231pd	%%zmm28,%%zmm0,%%zmm2	\n\t	 vfmadd231pd	%%zmm28,%%zmm11,%%zmm9 	\n\t"\
		"vmovaps	     (%%r8),%%zmm0			\n\t	vmovaps	0x200(%%r8),%%zmm10				\n\t"\
		"vmovaps	0x040(%%r8),%%zmm1			\n\t	vmovaps	0x240(%%r8),%%zmm11				\n\t"\
		"vsubpd	%%zmm2,%%zmm0,%%zmm0			\n\t	vsubpd	%%zmm8 ,%%zmm10,%%zmm10			\n\t"\
		"vsubpd	%%zmm3,%%zmm1,%%zmm1			\n\t	vsubpd	%%zmm9 ,%%zmm11,%%zmm11			\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm0,%%zmm2	\n\t	vfmadd132pd	%%zmm31,%%zmm10,%%zmm8 		\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm1,%%zmm3	\n\t	vfmadd132pd	%%zmm31,%%zmm11,%%zmm9 		\n\t"\
		"\n\t"\
		/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c01): */			/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c05): */\
		"addq	$0x980,%%rsi	\n\t"/* c01; c05 += 0x200: */\
		"vfmsub132pd	%%zmm30,%%zmm4,%%zmm2	\n\t	vfmsub132pd	%%zmm30,%%zmm12,%%zmm10		\n\t"\
		"vfmsub132pd	%%zmm30,%%zmm7,%%zmm0	\n\t	vfmsub132pd	%%zmm30,%%zmm15,%%zmm8 		\n\t"\
		"vfmsub132pd	%%zmm30,%%zmm5,%%zmm3	\n\t	vfmsub132pd	%%zmm30,%%zmm13,%%zmm11		\n\t"\
		"vfmsub132pd	%%zmm30,%%zmm6,%%zmm1	\n\t	vfmsub132pd	%%zmm30,%%zmm14,%%zmm9 		\n\t"\
		"vmovaps		  (%%rsi),%%zmm16		\n\t		vmovaps		0x200(%%rsi),%%zmm20	\n\t"\
		"vmovaps	 0x080(%%rsi),%%zmm17		\n\t		vmovaps		0x280(%%rsi),%%zmm21	\n\t"\
		"vmovaps	 0x040(%%rsi),%%zmm18		\n\t		vmovaps		0x240(%%rsi),%%zmm22	\n\t"\
		"vmovaps	 0x0c0(%%rsi),%%zmm19		\n\t		vmovaps		0x2c0(%%rsi),%%zmm23	\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm2,%%zmm4		\n\t	vfmadd132pd	%%zmm31,%%zmm10,%%zmm12		\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm0,%%zmm7		\n\t	vfmadd132pd	%%zmm31,%%zmm8 ,%%zmm15		\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm3,%%zmm5		\n\t	vfmadd132pd	%%zmm31,%%zmm11,%%zmm13		\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm1,%%zmm6		\n\t	vfmadd132pd	%%zmm31,%%zmm9 ,%%zmm14		\n\t"\
		"vmovaps	%%zmm2,     (%%r8)			\n\t"\
		"vmovaps	%%zmm0,0x440(%%r8)			\n\t"\
		"vmovaps	%%zmm3,0x040(%%r8)			\n\t"\
		"vmovaps	%%zmm6,0x400(%%r8)			\n\t"\
		"vmulpd		%%zmm16,%%zmm4,%%zmm2		\n\t		vmulpd		%%zmm20,%%zmm12,%%zmm24	\n\t"\
		"vmulpd		%%zmm17,%%zmm7,%%zmm0		\n\t		vmulpd		%%zmm21,%%zmm15,%%zmm25	\n\t"\
		"vmulpd		%%zmm16,%%zmm5,%%zmm3		\n\t		vmulpd		%%zmm20,%%zmm13,%%zmm26	\n\t"\
		"vmulpd		%%zmm17,%%zmm1,%%zmm6		\n\t		vmulpd		%%zmm21,%%zmm9 ,%%zmm27	\n\t"\
	" vfmadd231pd	%%zmm18,%%zmm5,%%zmm2		\n\t	 vfmadd231pd	%%zmm22,%%zmm13,%%zmm24	\n\t"\
	" vfmadd231pd	%%zmm19,%%zmm1,%%zmm0		\n\t	 vfmadd231pd	%%zmm23,%%zmm9 ,%%zmm25	\n\t"\
	"vfnmadd231pd	%%zmm18,%%zmm4,%%zmm3		\n\t	vfnmadd231pd	%%zmm22,%%zmm12,%%zmm26	\n\t"\
	"vfnmadd231pd	%%zmm19,%%zmm7,%%zmm6		\n\t	vfnmadd231pd	%%zmm23,%%zmm15,%%zmm27	\n\t"\
		"vmovaps	%%zmm2,     (%%rax,%%rdi,8)	\n\t	vmovaps	%%zmm24,     (%%rax,%%r9,8)		\n\t"\
		"vmovaps	%%zmm0,     (%%rbx,%%rdi,8)	\n\t	vmovaps	%%zmm25,     (%%rbx,%%r9,8)		\n\t"\
		"vmovaps	%%zmm3,0x040(%%rax,%%rdi,8)	\n\t	vmovaps	%%zmm26,0x040(%%rax,%%r9,8)		\n\t"\
		"vmovaps	%%zmm6,0x040(%%rbx,%%rdi,8)	\n\t	vmovaps	%%zmm27,0x040(%%rbx,%%r9,8)		\n\t"\
		"addq	$0x100,%%rsi	\n\t"/* c11; c15+= 0x200: */\
		"vmovaps		  (%%rsi),%%zmm16		\n\t	vmovaps		0x200(%%rsi),%%zmm20		\n\t"\
		"vmovaps	 0x080(%%rsi),%%zmm17		\n\t	vmovaps		0x280(%%rsi),%%zmm21		\n\t"\
		"vmovaps	 0x040(%%rsi),%%zmm18		\n\t	vmovaps		0x240(%%rsi),%%zmm22		\n\t"\
		"vmovaps	 0x0c0(%%rsi),%%zmm19		\n\t	vmovaps		0x2c0(%%rsi),%%zmm23		\n\t"\
		"vmovaps	     (%%r8),%%zmm4			\n\t"\
		"vmovaps	0x440(%%r8),%%zmm0			\n\t"\
		"vmovaps	0x040(%%r8),%%zmm5			\n\t"\
		"vmovaps	0x400(%%r8),%%zmm6			\n\t"\
		"vmulpd		%%zmm16,%%zmm4,%%zmm2		\n\t		vmulpd		%%zmm20,%%zmm10,%%zmm12	\n\t"\
		"vmulpd		%%zmm17,%%zmm0,%%zmm1		\n\t		vmulpd		%%zmm21,%%zmm8 ,%%zmm9 	\n\t"\
		"vmulpd		%%zmm16,%%zmm5,%%zmm3		\n\t		vmulpd		%%zmm20,%%zmm11,%%zmm13	\n\t"\
		"vmulpd		%%zmm17,%%zmm6,%%zmm7		\n\t		vmulpd		%%zmm21,%%zmm14,%%zmm15	\n\t"\
	" vfmadd231pd	%%zmm18,%%zmm5,%%zmm2		\n\t	 vfmadd231pd	%%zmm22,%%zmm11,%%zmm12	\n\t"\
	" vfmadd231pd	%%zmm19,%%zmm6,%%zmm1		\n\t	 vfmadd231pd	%%zmm23,%%zmm14,%%zmm9 	\n\t"\
	"vfnmadd231pd	%%zmm18,%%zmm4,%%zmm3		\n\t	vfnmadd231pd	%%zmm22,%%zmm10,%%zmm13	\n\t"\
	"vfnmadd231pd	%%zmm19,%%zmm0,%%zmm7		\n\t	vfnmadd231pd	%%zmm23,%%zmm8 ,%%zmm15	\n\t"\
		"vmovaps	%%zmm2,     (%%rcx,%%rdi,8)	\n\t	vmovaps	%%zmm12,     (%%rcx,%%r9,8)		\n\t"\
		"vmovaps	%%zmm1,     (%%rdx,%%rdi,8)	\n\t	vmovaps	%%zmm9 ,     (%%rdx,%%r9,8)		\n\t"\
		"vmovaps	%%zmm3,0x040(%%rcx,%%rdi,8)	\n\t	vmovaps	%%zmm13,0x040(%%rcx,%%r9,8)		\n\t"\
		"vmovaps	%%zmm7,0x040(%%rdx,%%rdi,8)	\n\t	vmovaps	%%zmm15,0x040(%%rdx,%%r9,8)		\n\t"\
		"\n\t"\
	/*...Block 3: t04,t14,t24,t34*/					/*...Block 7: t0C,t1C,t2C,t3C*/\
		"movq	%[__isrt2],%%rsi				\n\t	subq	%%rdi,%%r9		\n\t"\
		"movslq	%[__p02],%%rdi					\n\t	addq	%%rdi,%%r9		\n\t"\
		"addq	$0x880,%%r8		/* r24 */		\n\t"\
		"addq	$0x040,%%rsi	/* cc0 */		\n\t"\
		"vmovaps	     (%%r8),%%zmm4			\n\t	vmovaps	0x200(%%r8),%%zmm12				\n\t"\
		"vmovaps	0x040(%%r8),%%zmm5			\n\t	vmovaps	0x240(%%r8),%%zmm13				\n\t"\
		"vmovaps	     (%%rsi),%%zmm2			\n\t"\
		"vmovaps	0x040(%%rsi),%%zmm3			\n\t"\
		"vmovaps	%%zmm4,%%zmm6				\n\t	vmovaps	%%zmm12,%%zmm14					\n\t"\
		"vmovaps	%%zmm5,%%zmm7				\n\t	vmovaps	%%zmm13,%%zmm15					\n\t"\
		"vmulpd	%%zmm2,%%zmm4,%%zmm4			\n\t	vmulpd	%%zmm3,%%zmm12,%%zmm12			\n\t"\
		"vmulpd	%%zmm2,%%zmm5,%%zmm5			\n\t	vmulpd	%%zmm3,%%zmm13,%%zmm13			\n\t"\
		"vmovaps	0x400(%%r8),%%zmm0			\n\t	vmovaps	0x600(%%r8),%%zmm8 				\n\t"\
		"vmovaps	0x440(%%r8),%%zmm1			\n\t	vmovaps	0x640(%%r8),%%zmm9 				\n\t"\
		" vfmadd231pd	%%zmm3,%%zmm7,%%zmm4	\n\t	 vfmadd231pd	%%zmm2,%%zmm15,%%zmm12	\n\t"\
		"vfnmadd231pd	%%zmm3,%%zmm6,%%zmm5	\n\t	vfnmadd231pd	%%zmm2,%%zmm14,%%zmm13	\n\t"\
		"vmovaps	%%zmm0,%%zmm6				\n\t	vmovaps	%%zmm8,%%zmm14					\n\t"\
		"vmovaps	%%zmm1,%%zmm7				\n\t	vmovaps	%%zmm9,%%zmm15					\n\t"\
		"vmulpd	%%zmm3,%%zmm0,%%zmm0			\n\t	vmulpd	%%zmm2 ,%%zmm8 ,%%zmm8 			\n\t"\
		"vmulpd	%%zmm3,%%zmm1,%%zmm1			\n\t	vmulpd	%%zmm2 ,%%zmm9 ,%%zmm9 			\n\t"\
		" vfmadd231pd	%%zmm2,%%zmm7,%%zmm0	\n\t	 vfmadd231pd	%%zmm3,%%zmm15,%%zmm8 	\n\t"\
		"vfnmadd231pd	%%zmm2,%%zmm6,%%zmm1	\n\t	vfnmadd231pd	%%zmm3,%%zmm14,%%zmm9 	\n\t"\
		"vmovaps	%%zmm5,%%zmm7				\n\t	vmovaps	%%zmm13,%%zmm15					\n\t"\
		"vmovaps	%%zmm4,%%zmm6				\n\t	vmovaps	%%zmm12,%%zmm14					\n\t"\
		"vfmadd132pd	%%zmm30,%%zmm0,%%zmm4	\n\t	vfmadd132pd	%%zmm30,%%zmm8 ,%%zmm14		\n\t"\
		"vfmadd132pd	%%zmm30,%%zmm1,%%zmm5	\n\t	vfmadd132pd	%%zmm30,%%zmm9 ,%%zmm15		\n\t"\
		"vfmsub132pd	%%zmm30,%%zmm0,%%zmm6	\n\t	vfmsub132pd	%%zmm30,%%zmm8 ,%%zmm12		\n\t"\
		"vfmsub132pd	%%zmm30,%%zmm1,%%zmm7	\n\t	vfmsub132pd	%%zmm30,%%zmm9 ,%%zmm13		\n\t"\
		"subq	$0x800,%%r8						\n\t"\
		"subq	$0x040,%%rsi					\n\t"/* isrt2 */\
		"vmovaps	0x400(%%r8),%%zmm2			\n\t	vmovaps	0x600(%%r8),%%zmm8 				\n\t"\
		"vmovaps	0x440(%%r8),%%zmm3			\n\t	vmovaps	0x640(%%r8),%%zmm9 				\n\t"\
		"vmovaps	%%zmm3,%%zmm0				\n\t	vmovaps	%%zmm8 ,%%zmm10					\n\t"\
		"vsubpd	%%zmm2,%%zmm3,%%zmm3			\n\t	vsubpd	%%zmm9 ,%%zmm8,%%zmm8 			\n\t"\
		"vaddpd	%%zmm0,%%zmm2,%%zmm2			\n\t	vaddpd	%%zmm10,%%zmm9,%%zmm9 			\n\t"\
		"vmulpd	%%zmm29,%%zmm2,%%zmm2			\n\t	vmulpd	%%zmm29 ,%%zmm8 	,%%zmm8 	\n\t"/* isrt2 */\
		"vmulpd	%%zmm29,%%zmm3,%%zmm3			\n\t	vmulpd	%%zmm29 ,%%zmm9 	,%%zmm9 	\n\t"\
		"vmovaps	     (%%r8),%%zmm0			\n\t	vmovaps	0x200(%%r8),%%zmm10				\n\t"\
		"vmovaps	0x040(%%r8),%%zmm1			\n\t	vmovaps	0x240(%%r8),%%zmm11				\n\t"\
		"vsubpd	%%zmm2,%%zmm0,%%zmm0			\n\t	vsubpd	%%zmm8,%%zmm10,%%zmm10			\n\t"\
		"vsubpd	%%zmm3,%%zmm1,%%zmm1			\n\t	vsubpd	%%zmm9,%%zmm11,%%zmm11			\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm0,%%zmm2	\n\t	vfmadd132pd	%%zmm31,%%zmm10,%%zmm8 		\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm1,%%zmm3	\n\t	vfmadd132pd	%%zmm31,%%zmm11,%%zmm9 		\n\t"\
		"\n\t"\
		/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c02): */		/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c06): */\
		"addq	$0x5c0,%%rsi	\n\t"/* c02; c06 += 0x200: */\
	"vfmsub132pd	%%zmm30,%%zmm4,%%zmm2		\n\t	vfmsub132pd	%%zmm30,%%zmm12,%%zmm10		\n\t"\
	"vfmsub132pd	%%zmm30,%%zmm7,%%zmm0		\n\t	vfmsub132pd	%%zmm30,%%zmm15,%%zmm8 		\n\t"\
	"vfmsub132pd	%%zmm30,%%zmm5,%%zmm3		\n\t	vfmsub132pd	%%zmm30,%%zmm13,%%zmm11		\n\t"\
	"vfmsub132pd	%%zmm30,%%zmm6,%%zmm1		\n\t	vfmsub132pd	%%zmm30,%%zmm14,%%zmm9 		\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm2,%%zmm4		\n\t	vfmadd132pd	%%zmm31,%%zmm10,%%zmm12		\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm0,%%zmm7		\n\t	vfmadd132pd	%%zmm31,%%zmm8 ,%%zmm15		\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm3,%%zmm5		\n\t	vfmadd132pd	%%zmm31,%%zmm11,%%zmm13		\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm1,%%zmm6		\n\t	vfmadd132pd	%%zmm31,%%zmm9 ,%%zmm14		\n\t"\
		"vmovaps		  (%%rsi),%%zmm16		\n\t		vmovaps		0x200(%%rsi),%%zmm20	\n\t"\
		"vmovaps	 0x080(%%rsi),%%zmm17		\n\t		vmovaps		0x280(%%rsi),%%zmm21	\n\t"\
		"vmovaps	 0x040(%%rsi),%%zmm18		\n\t		vmovaps		0x240(%%rsi),%%zmm22	\n\t"\
		"vmovaps	 0x0c0(%%rsi),%%zmm19		\n\t		vmovaps		0x2c0(%%rsi),%%zmm23	\n\t"\
		"vmovaps	%%zmm2,     (%%r8)			\n\t"\
		"vmovaps	%%zmm0,0x440(%%r8)			\n\t"\
		"vmovaps	%%zmm3,0x040(%%r8)			\n\t"\
		"vmovaps	%%zmm6,0x400(%%r8)			\n\t"\
		"vmulpd		%%zmm16,%%zmm4,%%zmm2		\n\t		vmulpd		%%zmm20,%%zmm12,%%zmm24	\n\t"\
		"vmulpd		%%zmm17,%%zmm7,%%zmm0		\n\t		vmulpd		%%zmm21,%%zmm15,%%zmm25	\n\t"\
		"vmulpd		%%zmm16,%%zmm5,%%zmm3		\n\t		vmulpd		%%zmm20,%%zmm13,%%zmm26	\n\t"\
		"vmulpd		%%zmm17,%%zmm1,%%zmm6		\n\t		vmulpd		%%zmm21,%%zmm9 ,%%zmm27	\n\t"\
	" vfmadd231pd	%%zmm18,%%zmm5,%%zmm2		\n\t	 vfmadd231pd	%%zmm22,%%zmm13,%%zmm24	\n\t"\
	" vfmadd231pd	%%zmm19,%%zmm1,%%zmm0		\n\t	 vfmadd231pd	%%zmm23,%%zmm9 ,%%zmm25	\n\t"\
	"vfnmadd231pd	%%zmm18,%%zmm4,%%zmm3		\n\t	vfnmadd231pd	%%zmm22,%%zmm12,%%zmm26	\n\t"\
	"vfnmadd231pd	%%zmm19,%%zmm7,%%zmm6		\n\t	vfnmadd231pd	%%zmm23,%%zmm15,%%zmm27	\n\t"\
		"vmovaps	%%zmm2,     (%%rax,%%rdi,8)	\n\t	vmovaps	%%zmm24,     (%%rax,%%r9,8)		\n\t"\
		"vmovaps	%%zmm0,     (%%rbx,%%rdi,8)	\n\t	vmovaps	%%zmm25,     (%%rbx,%%r9,8)		\n\t"\
		"vmovaps	%%zmm3,0x040(%%rax,%%rdi,8)	\n\t	vmovaps	%%zmm26,0x040(%%rax,%%r9,8)		\n\t"\
		"vmovaps	%%zmm6,0x040(%%rbx,%%rdi,8)	\n\t	vmovaps	%%zmm27,0x040(%%rbx,%%r9,8)		\n\t"\
		"addq	$0x100,%%rsi	\n\t"/* c12; c16+= 0x200: */\
		"vmovaps		  (%%rsi),%%zmm16		\n\t	vmovaps		0x200(%%rsi),%%zmm20		\n\t"\
		"vmovaps	 0x080(%%rsi),%%zmm17		\n\t	vmovaps		0x280(%%rsi),%%zmm21		\n\t"\
		"vmovaps	 0x040(%%rsi),%%zmm18		\n\t	vmovaps		0x240(%%rsi),%%zmm22		\n\t"\
		"vmovaps	 0x0c0(%%rsi),%%zmm19		\n\t	vmovaps		0x2c0(%%rsi),%%zmm23		\n\t"\
		"vmovaps	     (%%r8),%%zmm4			\n\t"\
		"vmovaps	0x440(%%r8),%%zmm0			\n\t"\
		"vmovaps	0x040(%%r8),%%zmm5			\n\t"\
		"vmovaps	0x400(%%r8),%%zmm6			\n\t"\
		"vmulpd		%%zmm16,%%zmm4,%%zmm2		\n\t		vmulpd		%%zmm20,%%zmm10,%%zmm12	\n\t"\
		"vmulpd		%%zmm17,%%zmm0,%%zmm1		\n\t		vmulpd		%%zmm21,%%zmm8 ,%%zmm9 	\n\t"\
		"vmulpd		%%zmm16,%%zmm5,%%zmm3		\n\t		vmulpd		%%zmm20,%%zmm11,%%zmm13	\n\t"\
		"vmulpd		%%zmm17,%%zmm6,%%zmm7		\n\t		vmulpd		%%zmm21,%%zmm14,%%zmm15	\n\t"\
	" vfmadd231pd	%%zmm18,%%zmm5,%%zmm2		\n\t	 vfmadd231pd	%%zmm22,%%zmm11,%%zmm12	\n\t"\
	" vfmadd231pd	%%zmm19,%%zmm6,%%zmm1		\n\t	 vfmadd231pd	%%zmm23,%%zmm14,%%zmm9 	\n\t"\
	"vfnmadd231pd	%%zmm18,%%zmm4,%%zmm3		\n\t	vfnmadd231pd	%%zmm22,%%zmm10,%%zmm13	\n\t"\
	"vfnmadd231pd	%%zmm19,%%zmm0,%%zmm7		\n\t	vfnmadd231pd	%%zmm23,%%zmm8 ,%%zmm15	\n\t"\
		"vmovaps	%%zmm2,     (%%rcx,%%rdi,8)	\n\t	vmovaps	%%zmm12,     (%%rcx,%%r9,8)		\n\t"\
		"vmovaps	%%zmm1,     (%%rdx,%%rdi,8)	\n\t	vmovaps	%%zmm9 ,     (%%rdx,%%r9,8)		\n\t"\
		"vmovaps	%%zmm3,0x040(%%rcx,%%rdi,8)	\n\t	vmovaps	%%zmm13,0x040(%%rcx,%%r9,8)		\n\t"\
		"vmovaps	%%zmm7,0x040(%%rdx,%%rdi,8)	\n\t	vmovaps	%%zmm15,0x040(%%rdx,%%r9,8)		\n\t"\
		"\n\t"\
	/*...Block 4: t06,t16,t26,t36*/					/*...Block 8: t0E,t1E,t2E,t3E*/\
		"movq	%[__isrt2],%%rsi				\n\t	subq	%%rdi,%%r9	\n\t"\
		"movslq	%[__p03],%%rdi					\n\t	addq	%%rdi,%%r9	\n\t"\
		"addq	$0x880,%%r8		/* r26 */		\n\t"\
		"addq	$0x0c0,%%rsi					\n\t"/* cc1; cc3 += 0x080: */\
		"vmovaps	     (%%r8),%%zmm4			\n\t	vmovaps	0x200(%%r8),%%zmm12				\n\t"\
		"vmovaps	0x040(%%r8),%%zmm5			\n\t	vmovaps	0x240(%%r8),%%zmm13				\n\t"\
		"vmovaps	     (%%rsi),%%zmm2			\n\t	vmovaps	0x080(%%rsi),%%zmm10			\n\t"/* cc1,cc3 */\
		"vmovaps	0x040(%%rsi),%%zmm3			\n\t	vmovaps	0x0c0(%%rsi),%%zmm11			\n\t"\
		"vmovaps	%%zmm4,%%zmm6				\n\t	vmovaps	%%zmm12,%%zmm14					\n\t"\
		"vmovaps	%%zmm5,%%zmm7				\n\t	vmovaps	%%zmm13,%%zmm15					\n\t"\
		"vmulpd		%%zmm10,%%zmm4,%%zmm4		\n\t	vmulpd		%%zmm3,%%zmm12,%%zmm12		\n\t"\
		"vmulpd		%%zmm10,%%zmm5,%%zmm5		\n\t	vmulpd		%%zmm3,%%zmm13,%%zmm13		\n\t"\
		"vmovaps	0x400(%%r8),%%zmm0			\n\t	vmovaps	0x600(%%r8),%%zmm8 				\n\t"\
		"vmovaps	0x440(%%r8),%%zmm1			\n\t	vmovaps	0x640(%%r8),%%zmm9 				\n\t"\
		" vfmadd231pd	%%zmm11,%%zmm7,%%zmm4	\n\t	 vfmadd231pd	%%zmm2,%%zmm15,%%zmm12	\n\t"\
		"vfnmadd231pd	%%zmm11,%%zmm6,%%zmm5	\n\t	vfnmadd231pd	%%zmm2,%%zmm14,%%zmm13	\n\t"\
		"vmovaps	%%zmm0,%%zmm6				\n\t	vmovaps	%%zmm8 ,%%zmm14					\n\t"\
		"vmovaps	%%zmm1,%%zmm7				\n\t	vmovaps	%%zmm9 ,%%zmm15					\n\t"\
		"vmulpd		%%zmm3,%%zmm0,%%zmm0		\n\t	vmulpd		%%zmm11,%%zmm8 ,%%zmm8 		\n\t"\
		"vmulpd		%%zmm3,%%zmm1,%%zmm1		\n\t	vmulpd		%%zmm11,%%zmm9 ,%%zmm9 		\n\t"\
		"vfnmadd231pd	%%zmm2,%%zmm7,%%zmm0	\n\t	 vfmadd231pd	%%zmm10,%%zmm15,%%zmm8 	\n\t"\
		" vfmadd231pd	%%zmm2,%%zmm6,%%zmm1	\n\t	vfnmadd231pd	%%zmm10,%%zmm14,%%zmm9 	\n\t"\
		"vmovaps	%%zmm5,%%zmm7				\n\t	vmovaps	%%zmm13,%%zmm15					\n\t"\
		"vmovaps	%%zmm4,%%zmm6				\n\t	vmovaps	%%zmm12,%%zmm14					\n\t"\
		"vfmadd132pd	%%zmm30,%%zmm0,%%zmm6	\n\t	vfmsub132pd	%%zmm30,%%zmm8 ,%%zmm12		\n\t"\
		"vfmadd132pd	%%zmm30,%%zmm1,%%zmm7	\n\t	vfmsub132pd	%%zmm30,%%zmm9 ,%%zmm13		\n\t"\
		"vfmsub132pd	%%zmm30,%%zmm0,%%zmm4	\n\t	vfmadd132pd	%%zmm30,%%zmm8 ,%%zmm14		\n\t"\
		"vfmsub132pd	%%zmm30,%%zmm1,%%zmm5	\n\t	vfmadd132pd	%%zmm30,%%zmm9 ,%%zmm15		\n\t"\
		"subq	$0x800,%%r8						\n\t"\
		"subq	$0x080,%%rsi					\n\t"/* cc0 */\
		"vmovaps	0x400(%%r8),%%zmm2			\n\t	vmovaps	0x600(%%r8),%%zmm11				\n\t"\
		"vmovaps	0x440(%%r8),%%zmm0			\n\t	vmovaps	0x640(%%r8),%%zmm9 				\n\t"\
		"vmovaps	0x040(%%rsi),%%zmm3			\n\t	vmovaps	0x040(%%rsi),%%zmm10			\n\t"\
		"vmovaps	%%zmm2,%%zmm1				\n\t	vmovaps	%%zmm11,%%zmm8 					\n\t"\
		"vmulpd	%%zmm3 ,%%zmm2,%%zmm2			\n\t	vmulpd	%%zmm10,%%zmm11,%%zmm11			\n\t"\
		"vmulpd	%%zmm0 ,%%zmm3,%%zmm3			\n\t	vmulpd	%%zmm9 ,%%zmm10,%%zmm10			\n\t"\
		" vfmadd231pd	%%zmm28,%%zmm0,%%zmm2	\n\t	vfmadd132pd	%%zmm28,%%zmm11,%%zmm9 		\n\t"/* cc0 */\
		"vfnmadd231pd	%%zmm28,%%zmm1,%%zmm3	\n\t	vfmsub132pd	%%zmm28,%%zmm10,%%zmm8 		\n\t"\
		"vmovaps	     (%%r8),%%zmm0			\n\t	vmovaps	0x200(%%r8),%%zmm10				\n\t"\
		"vmovaps	0x040(%%r8),%%zmm1			\n\t	vmovaps	0x240(%%r8),%%zmm11				\n\t"\
		"vsubpd	%%zmm2,%%zmm0,%%zmm0			\n\t	vsubpd	%%zmm8,%%zmm10,%%zmm10			\n\t"\
		"vsubpd	%%zmm3,%%zmm1,%%zmm1			\n\t	vsubpd	%%zmm9,%%zmm11,%%zmm11			\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm0,%%zmm2	\n\t	vfmadd132pd	%%zmm31,%%zmm10,%%zmm8 		\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm1,%%zmm3	\n\t	vfmadd132pd	%%zmm31,%%zmm11,%%zmm9 		\n\t"\
		"\n\t"\
		/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c03): */	/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c07): */\
		"addq	$0xd80,%%rsi	\n\t"/* c03; c07 += 0x200: */\
	"vfmsub132pd	%%zmm30,%%zmm4,%%zmm2		\n\t	vfmsub132pd	%%zmm30,%%zmm12,%%zmm10		\n\t"\
	"vfmsub132pd	%%zmm30,%%zmm7,%%zmm0		\n\t	vfmsub132pd	%%zmm30,%%zmm15,%%zmm8 		\n\t"\
	"vfmsub132pd	%%zmm30,%%zmm5,%%zmm3		\n\t	vfmsub132pd	%%zmm30,%%zmm13,%%zmm11		\n\t"\
	"vfmsub132pd	%%zmm30,%%zmm6,%%zmm1		\n\t	vfmsub132pd	%%zmm30,%%zmm14,%%zmm9 		\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm2,%%zmm4		\n\t	vfmadd132pd	%%zmm31,%%zmm10,%%zmm12		\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm0,%%zmm7		\n\t	vfmadd132pd	%%zmm31,%%zmm8 ,%%zmm15		\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm3,%%zmm5		\n\t	vfmadd132pd	%%zmm31,%%zmm11,%%zmm13		\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm1,%%zmm6		\n\t	vfmadd132pd	%%zmm31,%%zmm9 ,%%zmm14		\n\t"\
		"vmovaps		  (%%rsi),%%zmm16		\n\t		vmovaps		0x200(%%rsi),%%zmm20	\n\t"\
		"vmovaps	 0x080(%%rsi),%%zmm17		\n\t		vmovaps		0x280(%%rsi),%%zmm21	\n\t"\
		"vmovaps	 0x040(%%rsi),%%zmm18		\n\t		vmovaps		0x240(%%rsi),%%zmm22	\n\t"\
		"vmovaps	 0x0c0(%%rsi),%%zmm19		\n\t		vmovaps		0x2c0(%%rsi),%%zmm23	\n\t"\
		"vmovaps	%%zmm2,     (%%r8)			\n\t"\
		"vmovaps	%%zmm0,0x440(%%r8)			\n\t"\
		"vmovaps	%%zmm3,0x040(%%r8)			\n\t"\
		"vmovaps	%%zmm6,0x400(%%r8)			\n\t"\
		"vmulpd		%%zmm16,%%zmm4,%%zmm2		\n\t		vmulpd		%%zmm20,%%zmm12,%%zmm24	\n\t"\
		"vmulpd		%%zmm17,%%zmm7,%%zmm0		\n\t		vmulpd		%%zmm21,%%zmm15,%%zmm25	\n\t"\
		"vmulpd		%%zmm16,%%zmm5,%%zmm3		\n\t		vmulpd		%%zmm20,%%zmm13,%%zmm26	\n\t"\
		"vmulpd		%%zmm17,%%zmm1,%%zmm6		\n\t		vmulpd		%%zmm21,%%zmm9 ,%%zmm27	\n\t"\
	" vfmadd231pd	%%zmm18,%%zmm5,%%zmm2		\n\t	 vfmadd231pd	%%zmm22,%%zmm13,%%zmm24	\n\t"\
	" vfmadd231pd	%%zmm19,%%zmm1,%%zmm0		\n\t	 vfmadd231pd	%%zmm23,%%zmm9 ,%%zmm25	\n\t"\
	"vfnmadd231pd	%%zmm18,%%zmm4,%%zmm3		\n\t	vfnmadd231pd	%%zmm22,%%zmm12,%%zmm26	\n\t"\
	"vfnmadd231pd	%%zmm19,%%zmm7,%%zmm6		\n\t	vfnmadd231pd	%%zmm23,%%zmm15,%%zmm27	\n\t"\
		"vmovaps	%%zmm2,     (%%rax,%%rdi,8)	\n\t	vmovaps	%%zmm24,     (%%rax,%%r9,8)		\n\t"\
		"vmovaps	%%zmm0,     (%%rbx,%%rdi,8)	\n\t	vmovaps	%%zmm25,     (%%rbx,%%r9,8)		\n\t"\
		"vmovaps	%%zmm3,0x040(%%rax,%%rdi,8)	\n\t	vmovaps	%%zmm26,0x040(%%rax,%%r9,8)		\n\t"\
		"vmovaps	%%zmm6,0x040(%%rbx,%%rdi,8)	\n\t	vmovaps	%%zmm27,0x040(%%rbx,%%r9,8)		\n\t"\
		"addq	$0x100,%%rsi	\n\t"/* c13; c17 += 0x200: */\
		"vmovaps		  (%%rsi),%%zmm16		\n\t		vmovaps		0x200(%%rsi),%%zmm20	\n\t"\
		"vmovaps	 0x080(%%rsi),%%zmm17		\n\t		vmovaps		0x280(%%rsi),%%zmm21	\n\t"\
		"vmovaps	 0x040(%%rsi),%%zmm18		\n\t		vmovaps		0x240(%%rsi),%%zmm22	\n\t"\
		"vmovaps	 0x0c0(%%rsi),%%zmm19		\n\t		vmovaps		0x2c0(%%rsi),%%zmm23	\n\t"\
		"vmovaps	     (%%r8),%%zmm4			\n\t"\
		"vmovaps	0x440(%%r8),%%zmm0			\n\t"\
		"vmovaps	0x040(%%r8),%%zmm5			\n\t"\
		"vmovaps	0x400(%%r8),%%zmm6			\n\t"\
		"vmulpd		%%zmm16,%%zmm4,%%zmm2		\n\t		vmulpd		%%zmm20,%%zmm10,%%zmm12	\n\t"\
		"vmulpd		%%zmm17,%%zmm0,%%zmm1		\n\t		vmulpd		%%zmm21,%%zmm8 ,%%zmm9 	\n\t"\
		"vmulpd		%%zmm16,%%zmm5,%%zmm3		\n\t		vmulpd		%%zmm20,%%zmm11,%%zmm13	\n\t"\
		"vmulpd		%%zmm17,%%zmm6,%%zmm7		\n\t		vmulpd		%%zmm21,%%zmm14,%%zmm15	\n\t"\
	" vfmadd231pd	%%zmm18,%%zmm5,%%zmm2		\n\t	 vfmadd231pd	%%zmm22,%%zmm11,%%zmm12	\n\t"\
	" vfmadd231pd	%%zmm19,%%zmm6,%%zmm1		\n\t	 vfmadd231pd	%%zmm23,%%zmm14,%%zmm9 	\n\t"\
	"vfnmadd231pd	%%zmm18,%%zmm4,%%zmm3		\n\t	vfnmadd231pd	%%zmm22,%%zmm10,%%zmm13	\n\t"\
	"vfnmadd231pd	%%zmm19,%%zmm0,%%zmm7		\n\t	vfnmadd231pd	%%zmm23,%%zmm8 ,%%zmm15	\n\t"\
		"vmovaps	%%zmm2,     (%%rcx,%%rdi,8)	\n\t	vmovaps	%%zmm12,     (%%rcx,%%r9,8)		\n\t"\
		"vmovaps	%%zmm1,     (%%rdx,%%rdi,8)	\n\t	vmovaps	%%zmm9 ,     (%%rdx,%%r9,8)		\n\t"\
		"vmovaps	%%zmm3,0x040(%%rcx,%%rdi,8)	\n\t	vmovaps	%%zmm13,0x040(%%rcx,%%r9,8)		\n\t"\
		"vmovaps	%%zmm7,0x040(%%rdx,%%rdi,8)	\n\t	vmovaps	%%zmm15,0x040(%%rdx,%%r9,8)		\n\t"\
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
		 ,[__r00] "m" (Xr00)\
		 ,[__isrt2] "m" (Xisrt2)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r8","r9","r14","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15","xmm16","xmm17","xmm18","xmm19","xmm20","xmm21","xmm22","xmm23","xmm24","xmm25","xmm26","xmm27","xmm28","xmm29","xmm30","xmm31"	/* Clobbered registers */\
	);\
	}

  #else	// FMA-using (but not ALL_FMA) version of AVX2 32-DIT macro, adapted to AVX-512

	#warning Using No-unity-multiplicand lower-FMA version of AVX-512 SSE2_RADIX32_DIT_TWIDDLE macro.
	// Cost [vector-ops only]: 374 MEM, 102 MUL, 272 FMA, 168 ADD
	#define SSE2_RADIX32_DIT_TWIDDLE(Xadd0,Xp01,Xp02,Xp03,Xp04,Xp05,Xp06,Xp07,Xp08,Xp10,Xr00,Xisrt2)\
	{\
	__asm__ volatile (\
		"movq	%[__add0],%%rax							\n\t	movq	%[__isrt2],%%r8		\n\t		movq	%[__r00],%%rsi	\n\t"\
		"movslq	%[__p01],%%rbx							\n\t	movslq	%[__p04],%%r9		\n\t		shlq	$3,%%r9	\n\t"\
		"movslq	%[__p02],%%rcx							\n\t	vmovaps		(%%r8),%%zmm29	\n\t"/* isrt2 */\
		"movslq	%[__p03],%%rdx							\n\t	addq	$0x1200,%%r8		\n\t"/* &two */\
		"leaq	(%%rax,%%rbx,8),%%rbx					\n\t	vmovaps	-0x40(%%r8),%%zmm30	\n\t"/* 1.0 */\
		"leaq	(%%rax,%%rcx,8),%%rcx					\n\t	vmovaps	     (%%r8),%%zmm31	\n\t"/* 2.0 */\
		"leaq	(%%rax,%%rdx,8),%%rdx					\n\t	vmovaps	 0x40(%%r8),%%zmm28	\n\t"/* sqrt2 */\
		/* SSE2_RADIX4_DIT_0TWIDDLE_B(r00) */				/* SSE2_RADIX4_DIT_0TWIDDLE_2NDOFTWO_B(r08) */\
		"vmovaps	     (%%rax),%%zmm2					\n\t	vmovaps	     (%%rax,%%r9),%%zmm10		\n\t"\
		"vmovaps	     (%%rbx),%%zmm0					\n\t	vmovaps	     (%%rbx,%%r9),%%zmm8 		\n\t"\
		"vmovaps	     (%%rcx),%%zmm6					\n\t	vmovaps	     (%%rcx,%%r9),%%zmm14		\n\t"\
		"vmovaps	     (%%rdx),%%zmm4					\n\t	vmovaps	     (%%rdx,%%r9),%%zmm12		\n\t"\
		"vmovaps	0x040(%%rax),%%zmm3					\n\t	vmovaps	0x040(%%rax,%%r9),%%zmm11		\n\t"\
		"vmovaps	0x040(%%rbx),%%zmm1					\n\t	vmovaps	0x040(%%rbx,%%r9),%%zmm9 		\n\t"\
		"vmovaps	0x040(%%rcx),%%zmm7					\n\t	vmovaps	0x040(%%rcx,%%r9),%%zmm15		\n\t"\
		"vmovaps	0x040(%%rdx),%%zmm5					\n\t	vmovaps	0x040(%%rdx,%%r9),%%zmm13		\n\t"\
		" vsubpd		%%zmm0,%%zmm2,%%zmm2			\n\t	 vsubpd		%%zmm8 ,%%zmm10,%%zmm10	\n\t"\
		" vsubpd		%%zmm4,%%zmm6,%%zmm6			\n\t	 vsubpd		%%zmm12,%%zmm14,%%zmm14	\n\t"\
		" vsubpd		%%zmm1,%%zmm3,%%zmm3			\n\t	 vsubpd		%%zmm9 ,%%zmm11,%%zmm11	\n\t"\
		" vsubpd		%%zmm5,%%zmm7,%%zmm7			\n\t	 vsubpd		%%zmm13,%%zmm15,%%zmm15	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm2,%%zmm0			\n\t	 vfmadd132pd	%%zmm31,%%zmm10,%%zmm8 	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm6,%%zmm4			\n\t	 vfmadd132pd	%%zmm31,%%zmm14,%%zmm12	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm3,%%zmm1			\n\t	 vfmadd132pd	%%zmm31,%%zmm11,%%zmm9 	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm7,%%zmm5			\n\t	 vfmadd132pd	%%zmm31,%%zmm15,%%zmm13	\n\t"\
		/* Finish radix-4 butterfly: */\
		" vsubpd		%%zmm4,%%zmm0,%%zmm0			\n\t	  vsubpd		%%zmm12,%%zmm8 ,%%zmm8 	\n\t"\
		" vsubpd		%%zmm7,%%zmm2,%%zmm2			\n\t	  vsubpd		%%zmm13,%%zmm9 ,%%zmm9 	\n\t"\
		" vsubpd		%%zmm5,%%zmm1,%%zmm1			\n\t	  vsubpd		%%zmm15,%%zmm10,%%zmm10	\n\t"\
		" vsubpd		%%zmm6,%%zmm3,%%zmm3			\n\t	  vsubpd		%%zmm14,%%zmm11,%%zmm11	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm0,%%zmm4			\n\t	 vfmadd132pd	%%zmm31,%%zmm8 ,%%zmm12	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm2,%%zmm7			\n\t	 vfmadd132pd	%%zmm31,%%zmm9 ,%%zmm13	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm1,%%zmm5			\n\t	 vfmadd132pd	%%zmm31,%%zmm10,%%zmm15	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm3,%%zmm6			\n\t	 vfmadd132pd	%%zmm31,%%zmm11,%%zmm14	\n\t"\
		" vsubpd		%%zmm12,%%zmm4,%%zmm4			\n\t	 vsubpd			%%zmm15,%%zmm11,%%zmm11	\n\t"\
		" vsubpd		%%zmm9 ,%%zmm0,%%zmm0			\n\t	 vsubpd			%%zmm10,%%zmm14,%%zmm14	\n\t"\
		"vsubpd			%%zmm13,%%zmm5 ,%%zmm5			\n\t	 vfmadd132pd	%%zmm31,%%zmm11,%%zmm15	\n\t"\
		"vsubpd			%%zmm8 ,%%zmm1 ,%%zmm1			\n\t	 vfmadd132pd	%%zmm31,%%zmm14,%%zmm10	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm4 ,%%zmm12			\n\t	vfnmadd231pd	%%zmm29,%%zmm11,%%zmm3	\n\t"/* isrt2 */\
		"vfmadd132pd	%%zmm31,%%zmm0 ,%%zmm9 			\n\t	vfnmadd231pd	%%zmm29,%%zmm14,%%zmm2	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm5 ,%%zmm13			\n\t	vfnmadd231pd	%%zmm29,%%zmm15,%%zmm7	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm1 ,%%zmm8 			\n\t	vfnmadd231pd	%%zmm29,%%zmm10,%%zmm6	\n\t"\
		/* SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS(r00,r02,r04,r06,r08,r0A,r0C,r0E) */\
		"vmovaps	%%zmm4 ,0x200(%%rsi)				\n\t	vfmadd132pd	%%zmm28,%%zmm3 ,%%zmm11	\n\t"/* sqrt2, use of FMA above means we must combine 2*isrt2 here */\
		"vmovaps	%%zmm0 ,0x300(%%rsi)				\n\t	vfmadd132pd	%%zmm28,%%zmm2 ,%%zmm14	\n\t"\
		"vmovaps	%%zmm5 ,0x240(%%rsi)				\n\t	vfmadd132pd	%%zmm28,%%zmm7 ,%%zmm15	\n\t"\
		"vmovaps	%%zmm1 ,0x140(%%rsi)				\n\t	vfmadd132pd	%%zmm28,%%zmm6 ,%%zmm10	\n\t"\
		"vmovaps	%%zmm12,     (%%rsi)				\n\t	vmovaps	%%zmm7 ,0x280(%%rsi)	\n\t"\
		"vmovaps	%%zmm9 ,0x100(%%rsi)				\n\t	vmovaps	%%zmm2 ,0x380(%%rsi)	\n\t"\
		"vmovaps	%%zmm13,0x040(%%rsi)				\n\t	vmovaps	%%zmm3 ,0x2c0(%%rsi)	\n\t"\
		"vmovaps	%%zmm8 ,0x340(%%rsi)				\n\t	vmovaps	%%zmm6 ,0x1c0(%%rsi)	\n\t"\
		"														vmovaps	%%zmm15,0x080(%%rsi)	\n\t"\
		"														vmovaps	%%zmm14,0x180(%%rsi)	\n\t"\
		"														vmovaps	%%zmm11,0x0c0(%%rsi)	\n\t"\
		"														vmovaps	%%zmm10,0x3c0(%%rsi)	\n\t"\
	/*...Block 2: */\
		"addq	$0x400,%%rsi							\n\t"/* r10 */\
		"movslq	%[__p08],%%rdi	\n\t	shlq $3,%%rdi	\n\t	addq	%%rdi,%%r9	\n\t"\
		/* SSE2_RADIX4_DIT_0TWIDDLE_B(r10) */			/* SSE2_RADIX4_DIT_0TWIDDLE_2NDOFTWO_B(r18) */\
		"vmovaps	     (%%rax,%%rdi),%%zmm2			\n\t	vmovaps	     (%%rax,%%r9),%%zmm10		\n\t"\
		"vmovaps	     (%%rbx,%%rdi),%%zmm0			\n\t	vmovaps	     (%%rbx,%%r9),%%zmm8 		\n\t"\
		"vmovaps	     (%%rcx,%%rdi),%%zmm6			\n\t	vmovaps	     (%%rcx,%%r9),%%zmm14		\n\t"\
		"vmovaps	     (%%rdx,%%rdi),%%zmm4			\n\t	vmovaps	     (%%rdx,%%r9),%%zmm12		\n\t"\
		"vmovaps	0x040(%%rax,%%rdi),%%zmm3			\n\t	vmovaps	0x040(%%rax,%%r9),%%zmm11		\n\t"\
		"vmovaps	0x040(%%rbx,%%rdi),%%zmm1			\n\t	vmovaps	0x040(%%rbx,%%r9),%%zmm9 		\n\t"\
		"vmovaps	0x040(%%rcx,%%rdi),%%zmm7			\n\t	vmovaps	0x040(%%rcx,%%r9),%%zmm15		\n\t"\
		"vmovaps	0x040(%%rdx,%%rdi),%%zmm5			\n\t	vmovaps	0x040(%%rdx,%%r9),%%zmm13		\n\t"\
		" vsubpd		%%zmm0,%%zmm2,%%zmm2			\n\t	  vsubpd		%%zmm8 ,%%zmm10,%%zmm10	\n\t"\
		" vsubpd		%%zmm4,%%zmm6,%%zmm6			\n\t	  vsubpd		%%zmm12,%%zmm14,%%zmm14	\n\t"\
		" vsubpd		%%zmm1,%%zmm3,%%zmm3			\n\t	  vsubpd		%%zmm9 ,%%zmm11,%%zmm11	\n\t"\
		" vsubpd		%%zmm5,%%zmm7,%%zmm7			\n\t	  vsubpd		%%zmm13,%%zmm15,%%zmm15	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm2,%%zmm0			\n\t	 vfmadd132pd	%%zmm31,%%zmm10,%%zmm8 	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm6,%%zmm4			\n\t	 vfmadd132pd	%%zmm31,%%zmm14,%%zmm12	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm3,%%zmm1			\n\t	 vfmadd132pd	%%zmm31,%%zmm11,%%zmm9 	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm7,%%zmm5			\n\t	 vfmadd132pd	%%zmm31,%%zmm15,%%zmm13	\n\t"\
		/* Finish radix-4 butterfly: */\
		" vsubpd		%%zmm4,%%zmm0,%%zmm0			\n\t	  vsubpd		%%zmm12,%%zmm8 ,%%zmm8 	\n\t"\
		" vsubpd		%%zmm7,%%zmm2,%%zmm2			\n\t	  vsubpd		%%zmm13,%%zmm9 ,%%zmm9 	\n\t"\
		" vsubpd		%%zmm5,%%zmm1,%%zmm1			\n\t	  vsubpd		%%zmm15,%%zmm10,%%zmm10	\n\t"\
		" vsubpd		%%zmm6,%%zmm3,%%zmm3			\n\t	  vsubpd		%%zmm14,%%zmm11,%%zmm11	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm0,%%zmm4			\n\t	 vfmadd132pd	%%zmm31,%%zmm8 ,%%zmm12	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm2,%%zmm7			\n\t	 vfmadd132pd	%%zmm31,%%zmm9 ,%%zmm13	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm1,%%zmm5			\n\t	 vfmadd132pd	%%zmm31,%%zmm10,%%zmm15	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm3,%%zmm6			\n\t	 vfmadd132pd	%%zmm31,%%zmm11,%%zmm14	\n\t"\
		" vsubpd		%%zmm12,%%zmm4,%%zmm4			\n\t	 vsubpd			%%zmm15,%%zmm11,%%zmm11	\n\t"\
		" vsubpd		%%zmm9 ,%%zmm0,%%zmm0			\n\t	 vsubpd			%%zmm10,%%zmm14,%%zmm14	\n\t"\
		"vsubpd			%%zmm13,%%zmm5 ,%%zmm5			\n\t	 vfmadd132pd	%%zmm31,%%zmm11,%%zmm15	\n\t"\
		"vsubpd			%%zmm8 ,%%zmm1 ,%%zmm1			\n\t	 vfmadd132pd	%%zmm31,%%zmm14,%%zmm10	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm4 ,%%zmm12			\n\t	vfnmadd231pd	%%zmm29,%%zmm11,%%zmm3	\n\t"/* isrt2 */\
		"vfmadd132pd	%%zmm31,%%zmm0 ,%%zmm9 			\n\t	vfnmadd231pd	%%zmm29,%%zmm14,%%zmm2	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm5 ,%%zmm13			\n\t	vfnmadd231pd	%%zmm29,%%zmm15,%%zmm7	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm1 ,%%zmm8 			\n\t	vfnmadd231pd	%%zmm29,%%zmm10,%%zmm6	\n\t"\
		/* SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS(r10,r12,r14,r16,r18,r1A,r1C,r1E) */\
		"vmovaps	%%zmm4 ,0x200(%%rsi)				\n\t	vfmadd132pd	%%zmm28,%%zmm3 ,%%zmm11	\n\t"/* sqrt2, use of FMA above means we must combine 2*isrt2 here */\
		"vmovaps	%%zmm0 ,0x300(%%rsi)				\n\t	vfmadd132pd	%%zmm28,%%zmm2 ,%%zmm14	\n\t"\
		"vmovaps	%%zmm5 ,0x240(%%rsi)				\n\t	vfmadd132pd	%%zmm28,%%zmm7 ,%%zmm15	\n\t"\
		"vmovaps	%%zmm1 ,0x140(%%rsi)				\n\t	vfmadd132pd	%%zmm28,%%zmm6 ,%%zmm10	\n\t"\
		"vmovaps	%%zmm12,     (%%rsi)				\n\t	vmovaps	%%zmm7 ,0x280(%%rsi)	\n\t"\
		"vmovaps	%%zmm9 ,0x100(%%rsi)				\n\t	vmovaps	%%zmm2 ,0x380(%%rsi)	\n\t"\
		"vmovaps	%%zmm13,0x040(%%rsi)				\n\t	vmovaps	%%zmm3 ,0x2c0(%%rsi)	\n\t"\
		"vmovaps	%%zmm8 ,0x340(%%rsi)				\n\t	vmovaps	%%zmm6 ,0x1c0(%%rsi)	\n\t"\
		"														vmovaps	%%zmm15,0x080(%%rsi)	\n\t"\
		"														vmovaps	%%zmm14,0x180(%%rsi)	\n\t"\
		"														vmovaps	%%zmm11,0x0c0(%%rsi)	\n\t"\
		"														vmovaps	%%zmm10,0x3c0(%%rsi)	\n\t"\
		"\n\t"\
	/*...Block 3: */\
		"addq	$0x400,%%rsi /* r20 */					\n\t	subq	%%rdi,%%r9	\n\t"\
		"movslq	%[__p10],%%rdi	\n\t	shlq $3,%%rdi	\n\t	addq	%%rdi,%%r9	\n\t"\
		/* SSE2_RADIX4_DIT_0TWIDDLE_B(r20) */			/* SSE2_RADIX4_DIT_0TWIDDLE_2NDOFTWO_B(r28) */\
		"vmovaps	     (%%rax,%%rdi),%%zmm2			\n\t	vmovaps	     (%%rax,%%r9),%%zmm10		\n\t"\
		"vmovaps	     (%%rbx,%%rdi),%%zmm0			\n\t	vmovaps	     (%%rbx,%%r9),%%zmm8 		\n\t"\
		"vmovaps	     (%%rcx,%%rdi),%%zmm6			\n\t	vmovaps	     (%%rcx,%%r9),%%zmm14		\n\t"\
		"vmovaps	     (%%rdx,%%rdi),%%zmm4			\n\t	vmovaps	     (%%rdx,%%r9),%%zmm12		\n\t"\
		"vmovaps	0x040(%%rax,%%rdi),%%zmm3			\n\t	vmovaps	0x040(%%rax,%%r9),%%zmm11		\n\t"\
		"vmovaps	0x040(%%rbx,%%rdi),%%zmm1			\n\t	vmovaps	0x040(%%rbx,%%r9),%%zmm9 		\n\t"\
		"vmovaps	0x040(%%rcx,%%rdi),%%zmm7			\n\t	vmovaps	0x040(%%rcx,%%r9),%%zmm15		\n\t"\
		"vmovaps	0x040(%%rdx,%%rdi),%%zmm5			\n\t	vmovaps	0x040(%%rdx,%%r9),%%zmm13		\n\t"\
		" vsubpd		%%zmm0,%%zmm2,%%zmm2			\n\t	  vsubpd		%%zmm8 ,%%zmm10,%%zmm10	\n\t"\
		" vsubpd		%%zmm4,%%zmm6,%%zmm6			\n\t	  vsubpd		%%zmm12,%%zmm14,%%zmm14	\n\t"\
		" vsubpd		%%zmm1,%%zmm3,%%zmm3			\n\t	  vsubpd		%%zmm9 ,%%zmm11,%%zmm11	\n\t"\
		" vsubpd		%%zmm5,%%zmm7,%%zmm7			\n\t	  vsubpd		%%zmm13,%%zmm15,%%zmm15	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm2,%%zmm0			\n\t	 vfmadd132pd	%%zmm31,%%zmm10,%%zmm8 	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm6,%%zmm4			\n\t	 vfmadd132pd	%%zmm31,%%zmm14,%%zmm12	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm3,%%zmm1			\n\t	 vfmadd132pd	%%zmm31,%%zmm11,%%zmm9 	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm7,%%zmm5			\n\t	 vfmadd132pd	%%zmm31,%%zmm15,%%zmm13	\n\t"\
		/* Finish radix-4 butterfly: */\
		" vsubpd		%%zmm4,%%zmm0,%%zmm0			\n\t	  vsubpd		%%zmm12,%%zmm8 ,%%zmm8 	\n\t"\
		" vsubpd		%%zmm7,%%zmm2,%%zmm2			\n\t	  vsubpd		%%zmm13,%%zmm9 ,%%zmm9 	\n\t"\
		" vsubpd		%%zmm5,%%zmm1,%%zmm1			\n\t	  vsubpd		%%zmm15,%%zmm10,%%zmm10	\n\t"\
		" vsubpd		%%zmm6,%%zmm3,%%zmm3			\n\t	  vsubpd		%%zmm14,%%zmm11,%%zmm11	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm0,%%zmm4			\n\t	 vfmadd132pd	%%zmm31,%%zmm8 ,%%zmm12	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm2,%%zmm7			\n\t	 vfmadd132pd	%%zmm31,%%zmm9 ,%%zmm13	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm1,%%zmm5			\n\t	 vfmadd132pd	%%zmm31,%%zmm10,%%zmm15	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm3,%%zmm6			\n\t	 vfmadd132pd	%%zmm31,%%zmm11,%%zmm14	\n\t"\
		" vsubpd		%%zmm12,%%zmm4,%%zmm4			\n\t	 vsubpd			%%zmm15,%%zmm11,%%zmm11	\n\t"\
		" vsubpd		%%zmm9 ,%%zmm0,%%zmm0			\n\t	 vsubpd			%%zmm10,%%zmm14,%%zmm14	\n\t"\
		"vsubpd			%%zmm13,%%zmm5 ,%%zmm5			\n\t	 vfmadd132pd	%%zmm31,%%zmm11,%%zmm15	\n\t"\
		"vsubpd			%%zmm8 ,%%zmm1 ,%%zmm1			\n\t	 vfmadd132pd	%%zmm31,%%zmm14,%%zmm10	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm4 ,%%zmm12			\n\t	vfnmadd231pd	%%zmm29,%%zmm11,%%zmm3	\n\t"/* isrt2 */\
		"vfmadd132pd	%%zmm31,%%zmm0 ,%%zmm9 			\n\t	vfnmadd231pd	%%zmm29,%%zmm14,%%zmm2	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm5 ,%%zmm13			\n\t	vfnmadd231pd	%%zmm29,%%zmm15,%%zmm7	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm1 ,%%zmm8 			\n\t	vfnmadd231pd	%%zmm29,%%zmm10,%%zmm6	\n\t"\
		/* SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS(r20,r22,r24,r26,r28,r2A,r2C,r2E) */\
		"vmovaps	%%zmm4 ,0x200(%%rsi)				\n\t	vfmadd132pd	%%zmm28,%%zmm3 ,%%zmm11	\n\t"/* sqrt2, use of FMA above means we must combine 2*isrt2 here */\
		"vmovaps	%%zmm0 ,0x300(%%rsi)				\n\t	vfmadd132pd	%%zmm28,%%zmm2 ,%%zmm14	\n\t"\
		"vmovaps	%%zmm5 ,0x240(%%rsi)				\n\t	vfmadd132pd	%%zmm28,%%zmm7 ,%%zmm15	\n\t"\
		"vmovaps	%%zmm1 ,0x140(%%rsi)				\n\t	vfmadd132pd	%%zmm28,%%zmm6 ,%%zmm10	\n\t"\
		"vmovaps	%%zmm12,     (%%rsi)				\n\t	vmovaps	%%zmm7 ,0x280(%%rsi)	\n\t"\
		"vmovaps	%%zmm9 ,0x100(%%rsi)				\n\t	vmovaps	%%zmm2 ,0x380(%%rsi)	\n\t"\
		"vmovaps	%%zmm13,0x040(%%rsi)				\n\t	vmovaps	%%zmm3 ,0x2c0(%%rsi)	\n\t"\
		"vmovaps	%%zmm8 ,0x340(%%rsi)				\n\t	vmovaps	%%zmm6 ,0x1c0(%%rsi)	\n\t"\
		"														vmovaps	%%zmm15,0x080(%%rsi)	\n\t"\
		"														vmovaps	%%zmm14,0x180(%%rsi)	\n\t"\
		"														vmovaps	%%zmm11,0x0c0(%%rsi)	\n\t"\
		"														vmovaps	%%zmm10,0x3c0(%%rsi)	\n\t"\
		"\n\t"\
	/*...Block 4: */\
		"addq	$0x400,%%rsi /* r30 */					\n\t	subq	%%rdi,%%r9	\n\t"\
		"movslq	%[__p08],%%r14					\n\t"\
		"movslq	%[__p10],%%rdi					\n\t"\
/* p18 */"addq	%%r14,%%rdi		\n\t	shlq $3,%%rdi	\n\t	addq	%%rdi,%%r9	\n\t"\
		/* SSE2_RADIX4_DIT_0TWIDDLE_B(r30) */		/* SSE2_RADIX4_DIT_0TWIDDLE_2NDOFTWO_B(r38) */\
		"vmovaps	     (%%rax,%%rdi),%%zmm2			\n\t	vmovaps	     (%%rax,%%r9),%%zmm10		\n\t"\
		"vmovaps	     (%%rbx,%%rdi),%%zmm0			\n\t	vmovaps	     (%%rbx,%%r9),%%zmm8 		\n\t"\
		"vmovaps	     (%%rcx,%%rdi),%%zmm6			\n\t	vmovaps	     (%%rcx,%%r9),%%zmm14		\n\t"\
		"vmovaps	     (%%rdx,%%rdi),%%zmm4			\n\t	vmovaps	     (%%rdx,%%r9),%%zmm12		\n\t"\
		"vmovaps	0x040(%%rax,%%rdi),%%zmm3			\n\t	vmovaps	0x040(%%rax,%%r9),%%zmm11		\n\t"\
		"vmovaps	0x040(%%rbx,%%rdi),%%zmm1			\n\t	vmovaps	0x040(%%rbx,%%r9),%%zmm9 		\n\t"\
		"vmovaps	0x040(%%rcx,%%rdi),%%zmm7			\n\t	vmovaps	0x040(%%rcx,%%r9),%%zmm15		\n\t"\
		"vmovaps	0x040(%%rdx,%%rdi),%%zmm5			\n\t	vmovaps	0x040(%%rdx,%%r9),%%zmm13		\n\t"\
		" vsubpd		%%zmm0,%%zmm2,%%zmm2			\n\t	  vsubpd		%%zmm8 ,%%zmm10,%%zmm10	\n\t"\
		" vsubpd		%%zmm4,%%zmm6,%%zmm6			\n\t	  vsubpd		%%zmm12,%%zmm14,%%zmm14	\n\t"\
		" vsubpd		%%zmm1,%%zmm3,%%zmm3			\n\t	  vsubpd		%%zmm9 ,%%zmm11,%%zmm11	\n\t"\
		" vsubpd		%%zmm5,%%zmm7,%%zmm7			\n\t	  vsubpd		%%zmm13,%%zmm15,%%zmm15	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm2,%%zmm0			\n\t	 vfmadd132pd	%%zmm31,%%zmm10,%%zmm8 	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm6,%%zmm4			\n\t	 vfmadd132pd	%%zmm31,%%zmm14,%%zmm12	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm3,%%zmm1			\n\t	 vfmadd132pd	%%zmm31,%%zmm11,%%zmm9 	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm7,%%zmm5			\n\t	 vfmadd132pd	%%zmm31,%%zmm15,%%zmm13	\n\t"\
		/* Finish radix-4 butterfly: */\
		" vsubpd		%%zmm4,%%zmm0,%%zmm0			\n\t	  vsubpd		%%zmm12,%%zmm8 ,%%zmm8 	\n\t"\
		" vsubpd		%%zmm7,%%zmm2,%%zmm2			\n\t	  vsubpd		%%zmm13,%%zmm9 ,%%zmm9 	\n\t"\
		" vsubpd		%%zmm5,%%zmm1,%%zmm1			\n\t	  vsubpd		%%zmm15,%%zmm10,%%zmm10	\n\t"\
		" vsubpd		%%zmm6,%%zmm3,%%zmm3			\n\t	  vsubpd		%%zmm14,%%zmm11,%%zmm11	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm0,%%zmm4			\n\t	 vfmadd132pd	%%zmm31,%%zmm8 ,%%zmm12	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm2,%%zmm7			\n\t	 vfmadd132pd	%%zmm31,%%zmm9 ,%%zmm13	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm1,%%zmm5			\n\t	 vfmadd132pd	%%zmm31,%%zmm10,%%zmm15	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm3,%%zmm6			\n\t	 vfmadd132pd	%%zmm31,%%zmm11,%%zmm14	\n\t"\
		" vsubpd		%%zmm12,%%zmm4,%%zmm4			\n\t	 vsubpd			%%zmm15,%%zmm11,%%zmm11	\n\t"\
		" vsubpd		%%zmm9 ,%%zmm0,%%zmm0			\n\t	 vsubpd			%%zmm10,%%zmm14,%%zmm14	\n\t"\
		"vsubpd			%%zmm13,%%zmm5 ,%%zmm5			\n\t	 vfmadd132pd	%%zmm31,%%zmm11,%%zmm15	\n\t"\
		"vsubpd			%%zmm8 ,%%zmm1 ,%%zmm1			\n\t	 vfmadd132pd	%%zmm31,%%zmm14,%%zmm10	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm4 ,%%zmm12			\n\t	vfnmadd231pd	%%zmm29,%%zmm11,%%zmm3	\n\t"/* isrt2 */\
		"vfmadd132pd	%%zmm31,%%zmm0 ,%%zmm9 			\n\t	vfnmadd231pd	%%zmm29,%%zmm14,%%zmm2	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm5 ,%%zmm13			\n\t	vfnmadd231pd	%%zmm29,%%zmm15,%%zmm7	\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm1 ,%%zmm8 			\n\t	vfnmadd231pd	%%zmm29,%%zmm10,%%zmm6	\n\t"\
		/* SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS(r30,r32,r34,r36,r38,r3A,r3C,r3E) */\
		"vmovaps	%%zmm4 ,0x200(%%rsi)				\n\t	vfmadd132pd	%%zmm28,%%zmm3 ,%%zmm11	\n\t"/* sqrt2, use of FMA above means we must combine 2*isrt2 here */\
		"vmovaps	%%zmm0 ,0x300(%%rsi)				\n\t	vfmadd132pd	%%zmm28,%%zmm2 ,%%zmm14	\n\t"\
		"vmovaps	%%zmm5 ,0x240(%%rsi)				\n\t	vfmadd132pd	%%zmm28,%%zmm7 ,%%zmm15	\n\t"\
		"vmovaps	%%zmm1 ,0x140(%%rsi)				\n\t	vfmadd132pd	%%zmm28,%%zmm6 ,%%zmm10	\n\t"\
		"vmovaps	%%zmm12,     (%%rsi)				\n\t	vmovaps	%%zmm7 ,0x280(%%rsi)	\n\t"\
		"vmovaps	%%zmm9 ,0x100(%%rsi)				\n\t	vmovaps	%%zmm2 ,0x380(%%rsi)	\n\t"\
		"vmovaps	%%zmm13,0x040(%%rsi)				\n\t	vmovaps	%%zmm3 ,0x2c0(%%rsi)	\n\t"\
		"vmovaps	%%zmm8 ,0x340(%%rsi)				\n\t	vmovaps	%%zmm6 ,0x1c0(%%rsi)	\n\t"\
		"														vmovaps	%%zmm15,0x080(%%rsi)	\n\t"\
		"														vmovaps	%%zmm14,0x180(%%rsi)	\n\t"\
		"														vmovaps	%%zmm11,0x0c0(%%rsi)	\n\t"\
		"														vmovaps	%%zmm10,0x3c0(%%rsi)	\n\t"\
		"\n\t"\
	/**********************************************************************************/\
	/*...and now do eight radix-4 transforms, including the internal twiddle factors: */\
	/**********************************************************************************/\
		"\n\t"\
	/*...Block 1: t00,t10,t20,t30	*/							/*...Block 5: t08,t18,t28,t38*/\
		"movq	%[__r00],%%r8	\n\t"/* r10,20,30 += 1,2,3*0x400 ... inputs to our octet of radix-4 pass 2 DFTs from this local-store */\
		"movq	%[__add0],%%rax					\n\t	movslq	%[__p04],%%r9				\n\t"\
		"movslq	%[__p08],%%rbx					\n\t	movq	%[__isrt2],%%rsi			\n\t"\
		"movslq	%[__p10],%%rcx					\n\t	vmovaps	0xa40(%%r8),%%zmm13			\n\t"\
		"leaq	(%%rbx,%%rcx),%%rdx	/* p18 */	\n\t	vmovaps	0xe00(%%r8),%%zmm14			\n\t"\
		"leaq	(%%rax,%%rbx,8),%%rbx			\n\t"/* add0 + p08 */\
		"leaq	(%%rax,%%rcx,8),%%rcx			\n\t"/* add0 + p10 */\
		"leaq	(%%rax,%%rdx,8),%%rdx			\n\t"/* add0 + p18 */\
		"vmovaps	0x400(%%r8),%%zmm2			\n\t	vmulpd	%%zmm29,%%zmm13,%%zmm13		\n\t"\
		"vmovaps	     (%%r8),%%zmm0			\n\t	vmulpd	%%zmm29,%%zmm14,%%zmm14		\n\t"\
		"vmovaps	0xc00(%%r8),%%zmm4			\n\t	vmovaps	0xa00(%%r8),%%zmm12			\n\t"\
		"vmovaps	0x800(%%r8),%%zmm6			\n\t	vmovaps	0xe40(%%r8),%%zmm15			\n\t"\
		"vmovaps	0x440(%%r8),%%zmm3			\n\t	vmovaps	0x200(%%r8),%%zmm8 			\n\t"\
		"vmovaps	0x040(%%r8),%%zmm1			\n\t	vmovaps	0x640(%%r8),%%zmm10			\n\t"\
		"vmovaps	0xc40(%%r8),%%zmm5			\n\t	vmovaps	0x240(%%r8),%%zmm11			\n\t"\
		"vmovaps	0x840(%%r8),%%zmm7			\n\t	vmovaps	0x600(%%r8),%%zmm9 			\n\t"\
	" vsubpd		%%zmm2,%%zmm0,%%zmm0		\n\t	vfnmadd231pd	%%zmm12,%%zmm29,%%zmm13	\n\t"\
	" vsubpd		%%zmm4,%%zmm6,%%zmm6		\n\t	vfnmadd231pd	%%zmm15,%%zmm29,%%zmm14	\n\t"\
	" vsubpd		%%zmm3,%%zmm1,%%zmm1		\n\t	 vfmadd132pd	%%zmm28,%%zmm13,%%zmm12	\n\t"/* sqrt2 */\
	" vsubpd		%%zmm5,%%zmm7,%%zmm7		\n\t	 vfmadd132pd	%%zmm28,%%zmm14,%%zmm15	\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm0,%%zmm2		\n\t	 vsubpd		%%zmm10,%%zmm8 ,%%zmm8 	\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm6,%%zmm4		\n\t	 vsubpd		%%zmm9 ,%%zmm11,%%zmm11	\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm1,%%zmm3		\n\t	 vsubpd		%%zmm14,%%zmm12,%%zmm12	\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm7,%%zmm5		\n\t	 vsubpd		%%zmm15,%%zmm13,%%zmm13	\n\t"\
		/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c00): */\
	" vsubpd		%%zmm4,%%zmm2,%%zmm2		\n\t	vfmadd132pd	%%zmm31,%%zmm8 ,%%zmm10		\n\t"\
	" vsubpd		%%zmm7,%%zmm0,%%zmm0		\n\t	vfmadd132pd	%%zmm31,%%zmm11,%%zmm9 		\n\t"\
	" vsubpd		%%zmm5,%%zmm3,%%zmm3		\n\t	vfmadd132pd	%%zmm31,%%zmm12,%%zmm14		\n\t"\
	" vsubpd		%%zmm6,%%zmm1,%%zmm1		\n\t	vfmadd132pd	%%zmm31,%%zmm13,%%zmm15		\n\t"\
		"addq	$0x1c0,%%rsi	\n\t"/* c00 */		/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c04 = c00 + 0x200): */\
	"vfmadd132pd	%%zmm31,%%zmm2,%%zmm4		\n\t	 vsubpd		%%zmm13,%%zmm11,%%zmm11		\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm0,%%zmm7		\n\t	 vsubpd		%%zmm12,%%zmm10,%%zmm10		\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm3,%%zmm5		\n\t	 vsubpd		%%zmm14,%%zmm9 ,%%zmm9 		\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm1,%%zmm6		\n\t	 vsubpd		%%zmm15,%%zmm8 ,%%zmm8 		\n\t"\
		"vmovaps	%%zmm2,     (%%r8)			\n\t	vfmadd132pd	%%zmm31,%%zmm11,%%zmm13		\n\t"\
		"vmovaps	%%zmm0,0x440(%%r8)			\n\t	vfmadd132pd	%%zmm31,%%zmm10,%%zmm12		\n\t"\
		"vmovaps	%%zmm3,0x040(%%r8)			\n\t	vfmadd132pd	%%zmm31,%%zmm9 ,%%zmm14		\n\t"\
		"vmovaps	%%zmm6,0x400(%%r8)			\n\t	vfmadd132pd	%%zmm31,%%zmm8 ,%%zmm15		\n\t"\
		"vmovaps		  (%%rsi),%%zmm16		\n\t	vmovaps		0x200(%%rsi),%%zmm20		\n\t"\
		"vmovaps	 0x080(%%rsi),%%zmm17		\n\t	vmovaps		0x280(%%rsi),%%zmm21		\n\t"\
		"vmovaps	 0x040(%%rsi),%%zmm18		\n\t	vmovaps		0x240(%%rsi),%%zmm22		\n\t"\
		"vmovaps	 0x0c0(%%rsi),%%zmm19		\n\t	vmovaps		0x2c0(%%rsi),%%zmm23		\n\t"\
		"vmulpd		%%zmm16,%%zmm5,%%zmm3		\n\t		vmulpd		%%zmm20,%%zmm13,%%zmm24	\n\t"\
		"vmulpd		%%zmm17,%%zmm1,%%zmm6		\n\t		vmulpd		%%zmm21,%%zmm9 ,%%zmm25	\n\t"\
		"vmulpd		%%zmm16,%%zmm4,%%zmm2		\n\t		vmulpd		%%zmm20,%%zmm12,%%zmm26	\n\t"\
		"vmulpd		%%zmm17,%%zmm7,%%zmm0		\n\t		vmulpd		%%zmm21,%%zmm15,%%zmm27	\n\t"\
	"vfnmadd231pd	%%zmm18,%%zmm4,%%zmm3		\n\t	vfnmadd231pd	%%zmm22,%%zmm12,%%zmm24	\n\t"\
	"vfnmadd231pd	%%zmm19,%%zmm7,%%zmm6		\n\t	vfnmadd231pd	%%zmm23,%%zmm15,%%zmm25	\n\t"\
	" vfmadd231pd	%%zmm18,%%zmm5,%%zmm2		\n\t	 vfmadd231pd	%%zmm22,%%zmm13,%%zmm26	\n\t"\
	" vfmadd231pd	%%zmm19,%%zmm1,%%zmm0		\n\t	 vfmadd231pd	%%zmm23,%%zmm9 ,%%zmm27	\n\t"\
		"vmovaps	%%zmm3,0x040(%%rax)			\n\t	vmovaps	%%zmm24,0x040(%%rax,%%r9,8)		\n\t"\
		"vmovaps	%%zmm6,0x040(%%rbx)			\n\t	vmovaps	%%zmm25,0x040(%%rbx,%%r9,8)		\n\t"\
		"vmovaps	%%zmm2,     (%%rax)			\n\t	vmovaps	%%zmm26,     (%%rax,%%r9,8)		\n\t"\
		"vmovaps	%%zmm0,     (%%rbx)			\n\t	vmovaps	%%zmm27,     (%%rbx,%%r9,8)		\n\t"\
		"vmovaps	0x040(%%r8),%%zmm5			\n\t"\
		"vmovaps	0x400(%%r8),%%zmm6			\n\t"\
		"vmovaps	     (%%r8),%%zmm4			\n\t"\
		"vmovaps	0x440(%%r8),%%zmm0			\n\t"\
		"addq	$0x100,%%rsi	\n\t"/* c10 */			/* c14 = c10 + 0x200: */\
		"vmovaps	     (%%rsi),%%zmm16		\n\t	vmovaps		0x200(%%rsi),%%zmm20		\n\t"\
		"vmovaps	0x080(%%rsi),%%zmm17		\n\t	vmovaps		0x280(%%rsi),%%zmm21		\n\t"\
		"vmovaps	0x040(%%rsi),%%zmm18		\n\t	vmovaps		0x240(%%rsi),%%zmm22		\n\t"\
		"vmovaps	0x0c0(%%rsi),%%zmm19		\n\t	vmovaps		0x2c0(%%rsi),%%zmm23		\n\t"\
		"vmulpd			%%zmm16,%%zmm5,%%zmm3	\n\t	vmulpd			%%zmm20,%%zmm11,%%zmm13	\n\t"\
		"vmulpd			%%zmm17,%%zmm6,%%zmm7	\n\t	vmulpd			%%zmm21,%%zmm14,%%zmm15	\n\t"\
		"vmulpd			%%zmm16,%%zmm4,%%zmm2	\n\t	vmulpd			%%zmm20,%%zmm10,%%zmm12	\n\t"\
		"vmulpd			%%zmm17,%%zmm0,%%zmm1	\n\t	vmulpd			%%zmm21,%%zmm8 ,%%zmm9 	\n\t"\
	"vfnmadd231pd		%%zmm18,%%zmm4,%%zmm3	\n\t	vfnmadd231pd	%%zmm22,%%zmm10,%%zmm13	\n\t"\
	"vfnmadd231pd		%%zmm19,%%zmm0,%%zmm7	\n\t	vfnmadd231pd	%%zmm23,%%zmm8 ,%%zmm15	\n\t"\
	" vfmadd231pd		%%zmm18,%%zmm5,%%zmm2	\n\t	 vfmadd231pd	%%zmm22,%%zmm11,%%zmm12	\n\t"\
	" vfmadd231pd		%%zmm19,%%zmm6,%%zmm1	\n\t	 vfmadd231pd	%%zmm23,%%zmm14,%%zmm9 	\n\t"\
		"vmovaps	%%zmm3,0x040(%%rcx)			\n\t	vmovaps	%%zmm13,0x040(%%rcx,%%r9,8)		\n\t"\
		"vmovaps	%%zmm7,0x040(%%rdx)			\n\t	vmovaps	%%zmm15,0x040(%%rdx,%%r9,8)		\n\t"\
		"vmovaps	%%zmm2,     (%%rcx)			\n\t	vmovaps	%%zmm12,     (%%rcx,%%r9,8)		\n\t"\
		"vmovaps	%%zmm1,     (%%rdx)			\n\t	vmovaps	%%zmm9 ,     (%%rdx,%%r9,8)		\n\t"\
		"\n\t"\
	/*...Block 2: t02,t12,t22,t32	*/					/*...Block 6: t0A,t1A,t2A,t3A*/\
		"movslq	%[__p01],%%rdi					\n\t	addq	%%rdi,%%r9		\n\t"\
		"movq	%[__isrt2],%%rsi				\n\t"\
		"addq	$0x880,%%r8						\n\t"/* r22 */\
		"addq	$0x0c0,%%rsi					\n\t"/* cc1; cc3 += 0x080: */\
		"vmovaps	     (%%r8),%%zmm4			\n\t	vmovaps	0x200(%%r8),%%zmm12				\n\t"\
		"vmovaps	0x040(%%r8),%%zmm5			\n\t	vmovaps	0x240(%%r8),%%zmm13				\n\t"\
		"vmovaps	     (%%rsi),%%zmm2			\n\t	vmovaps	0x080(%%rsi),%%zmm11			\n\t"/* cc1,cc3 */\
		"vmovaps	0x040(%%rsi),%%zmm3			\n\t	vmovaps	0x0c0(%%rsi),%%zmm10			\n\t"\
		"vmovaps	%%zmm4,%%zmm6				\n\t	vmovaps	%%zmm12,%%zmm14					\n\t"\
		"vmovaps	%%zmm5,%%zmm7				\n\t	vmovaps	%%zmm13,%%zmm15					\n\t"\
		"vmulpd		%%zmm2,%%zmm4,%%zmm4		\n\t	vmulpd		%%zmm10,%%zmm12,%%zmm12		\n\t"\
		"vmulpd		%%zmm2,%%zmm5,%%zmm5		\n\t	vmulpd		%%zmm10,%%zmm13,%%zmm13		\n\t"\
		"vmovaps	0x400(%%r8),%%zmm0			\n\t	vmovaps	0x600(%%r8),%%zmm8 				\n\t"\
		"vmovaps	0x440(%%r8),%%zmm1			\n\t	vmovaps	0x640(%%r8),%%zmm9 	/* t3B */	\n\t"\
		" vfmadd231pd	%%zmm3,%%zmm7,%%zmm4	\n\t	 vfmadd231pd	%%zmm11,%%zmm15,%%zmm12	\n\t"\
		"vfnmadd231pd	%%zmm3,%%zmm6,%%zmm5	\n\t	vfnmadd231pd	%%zmm11,%%zmm14,%%zmm13	\n\t"\
		"vmovaps	%%zmm0,%%zmm6				\n\t	vmovaps	%%zmm8 ,%%zmm14					\n\t"\
		"vmovaps	%%zmm1,%%zmm7				\n\t	vmovaps	%%zmm9 ,%%zmm15					\n\t"\
		"vmulpd		%%zmm11,%%zmm0,%%zmm0		\n\t	vmulpd		%%zmm2 ,%%zmm8 ,%%zmm8 		\n\t"\
		"vmulpd		%%zmm11,%%zmm1,%%zmm1		\n\t	vmulpd		%%zmm2 ,%%zmm9 ,%%zmm9 		\n\t"\
		" vfmadd231pd	%%zmm10,%%zmm7,%%zmm0	\n\t	vfnmadd231pd	%%zmm3 ,%%zmm15,%%zmm8 	\n\t"\
		"vfnmadd231pd	%%zmm10,%%zmm6,%%zmm1	\n\t	 vfmadd231pd	%%zmm3 ,%%zmm14,%%zmm9 	\n\t"\
		"vmovaps	%%zmm5,%%zmm7				\n\t	vmovaps	%%zmm13,%%zmm15					\n\t"\
		"vmovaps	%%zmm4,%%zmm6				\n\t	vmovaps	%%zmm12,%%zmm14					\n\t"\
		" vaddpd		%%zmm0,%%zmm4,%%zmm4	\n\t	 vaddpd		%%zmm8 ,%%zmm14,%%zmm14		\n\t"\
		" vaddpd		%%zmm1,%%zmm5,%%zmm5	\n\t	 vaddpd		%%zmm9 ,%%zmm15,%%zmm15		\n\t"\
		" vsubpd		%%zmm0,%%zmm6,%%zmm6	\n\t	 vsubpd		%%zmm8 ,%%zmm12,%%zmm12		\n\t"\
		" vsubpd		%%zmm1,%%zmm7,%%zmm7	\n\t	 vsubpd		%%zmm9 ,%%zmm13,%%zmm13		\n\t"\
		"subq	$0x800,%%r8		/* r02 */		\n\t"\
		"subq	$0x080,%%rsi	/* cc0 */		\n\t	vmovaps	     (%%rsi),%%zmm28			\n\t"/* cc0 */\
		"vmovaps	0x400(%%r8),%%zmm1			\n\t	vmovaps	0x600(%%r8),%%zmm8 				\n\t"\
		"vmovaps	0x440(%%r8),%%zmm3			\n\t	vmovaps	0x640(%%r8),%%zmm10				\n\t"\
		"vmovaps	0x040(%%rsi),%%zmm2			\n\t	vmovaps	0x040(%%rsi),%%zmm9 			\n\t"\
		"vmovaps	%%zmm1,%%zmm0				\n\t	vmovaps	%%zmm8 ,%%zmm11					\n\t"\
		"vmulpd		%%zmm2 ,%%zmm1,%%zmm1		\n\t	vmulpd		%%zmm9 ,%%zmm8 ,%%zmm8 		\n\t"\
		"vmulpd		%%zmm3 ,%%zmm2,%%zmm2		\n\t	vmulpd		%%zmm10,%%zmm9 ,%%zmm9 		\n\t"\
		" vfmsub132pd	%%zmm28,%%zmm1,%%zmm3	\n\t	vfnmadd231pd	%%zmm28,%%zmm10,%%zmm8 	\n\t"/* cc0 */\
		" vfmadd231pd	%%zmm28,%%zmm0,%%zmm2	\n\t	 vfmadd231pd	%%zmm28,%%zmm11,%%zmm9 	\n\t"\
		"vmovaps	     (%%r8),%%zmm0			\n\t	vmovaps	0x200(%%r8),%%zmm10				\n\t"\
		"vmovaps	0x040(%%r8),%%zmm1			\n\t	vmovaps	0x240(%%r8),%%zmm11				\n\t"\
		"vsubpd	%%zmm2,%%zmm0,%%zmm0			\n\t	vsubpd	%%zmm8 ,%%zmm10,%%zmm10			\n\t"\
		"vsubpd	%%zmm3,%%zmm1,%%zmm1			\n\t	vsubpd	%%zmm9 ,%%zmm11,%%zmm11			\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm0,%%zmm2	\n\t	vfmadd132pd	%%zmm31,%%zmm10,%%zmm8 		\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm1,%%zmm3	\n\t	vfmadd132pd	%%zmm31,%%zmm11,%%zmm9 		\n\t"\
		"\n\t"\
		/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c01): */			/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c05): */\
		"addq	$0x980,%%rsi	\n\t"/* c01; c05 += 0x200: */\
		" vsubpd		%%zmm4,%%zmm2,%%zmm2	\n\t	 vsubpd		%%zmm12,%%zmm10,%%zmm10		\n\t"\
		" vsubpd		%%zmm7,%%zmm0,%%zmm0	\n\t	 vsubpd		%%zmm15,%%zmm8 ,%%zmm8 		\n\t"\
		" vsubpd		%%zmm5,%%zmm3,%%zmm3	\n\t	 vsubpd		%%zmm13,%%zmm11,%%zmm11		\n\t"\
		" vsubpd		%%zmm6,%%zmm1,%%zmm1	\n\t	 vsubpd		%%zmm14,%%zmm9 ,%%zmm9 		\n\t"\
		"vmovaps		  (%%rsi),%%zmm16		\n\t		vmovaps		0x200(%%rsi),%%zmm20	\n\t"\
		"vmovaps	 0x080(%%rsi),%%zmm17		\n\t		vmovaps		0x280(%%rsi),%%zmm21	\n\t"\
		"vmovaps	 0x040(%%rsi),%%zmm18		\n\t		vmovaps		0x240(%%rsi),%%zmm22	\n\t"\
		"vmovaps	 0x0c0(%%rsi),%%zmm19		\n\t		vmovaps		0x2c0(%%rsi),%%zmm23	\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm2,%%zmm4		\n\t	vfmadd132pd	%%zmm31,%%zmm10,%%zmm12		\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm0,%%zmm7		\n\t	vfmadd132pd	%%zmm31,%%zmm8 ,%%zmm15		\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm3,%%zmm5		\n\t	vfmadd132pd	%%zmm31,%%zmm11,%%zmm13		\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm1,%%zmm6		\n\t	vfmadd132pd	%%zmm31,%%zmm9 ,%%zmm14		\n\t"\
		"vmovaps	%%zmm2,     (%%r8)			\n\t"\
		"vmovaps	%%zmm0,0x440(%%r8)			\n\t"\
		"vmovaps	%%zmm3,0x040(%%r8)			\n\t"\
		"vmovaps	%%zmm6,0x400(%%r8)			\n\t"\
		"vmulpd		%%zmm16,%%zmm4,%%zmm2		\n\t		vmulpd		%%zmm20,%%zmm12,%%zmm24	\n\t"\
		"vmulpd		%%zmm17,%%zmm7,%%zmm0		\n\t		vmulpd		%%zmm21,%%zmm15,%%zmm25	\n\t"\
		"vmulpd		%%zmm16,%%zmm5,%%zmm3		\n\t		vmulpd		%%zmm20,%%zmm13,%%zmm26	\n\t"\
		"vmulpd		%%zmm17,%%zmm1,%%zmm6		\n\t		vmulpd		%%zmm21,%%zmm9 ,%%zmm27	\n\t"\
	" vfmadd231pd	%%zmm18,%%zmm5,%%zmm2		\n\t	 vfmadd231pd	%%zmm22,%%zmm13,%%zmm24	\n\t"\
	" vfmadd231pd	%%zmm19,%%zmm1,%%zmm0		\n\t	 vfmadd231pd	%%zmm23,%%zmm9 ,%%zmm25	\n\t"\
	"vfnmadd231pd	%%zmm18,%%zmm4,%%zmm3		\n\t	vfnmadd231pd	%%zmm22,%%zmm12,%%zmm26	\n\t"\
	"vfnmadd231pd	%%zmm19,%%zmm7,%%zmm6		\n\t	vfnmadd231pd	%%zmm23,%%zmm15,%%zmm27	\n\t"\
		"vmovaps	%%zmm2,     (%%rax,%%rdi,8)	\n\t	vmovaps	%%zmm24,     (%%rax,%%r9,8)		\n\t"\
		"vmovaps	%%zmm0,     (%%rbx,%%rdi,8)	\n\t	vmovaps	%%zmm25,     (%%rbx,%%r9,8)		\n\t"\
		"vmovaps	%%zmm3,0x040(%%rax,%%rdi,8)	\n\t	vmovaps	%%zmm26,0x040(%%rax,%%r9,8)		\n\t"\
		"vmovaps	%%zmm6,0x040(%%rbx,%%rdi,8)	\n\t	vmovaps	%%zmm27,0x040(%%rbx,%%r9,8)		\n\t"\
		"addq	$0x100,%%rsi	\n\t"/* c11; c15+= 0x200: */\
		"vmovaps		  (%%rsi),%%zmm16		\n\t	vmovaps		0x200(%%rsi),%%zmm20		\n\t"\
		"vmovaps	 0x080(%%rsi),%%zmm17		\n\t	vmovaps		0x280(%%rsi),%%zmm21		\n\t"\
		"vmovaps	 0x040(%%rsi),%%zmm18		\n\t	vmovaps		0x240(%%rsi),%%zmm22		\n\t"\
		"vmovaps	 0x0c0(%%rsi),%%zmm19		\n\t	vmovaps		0x2c0(%%rsi),%%zmm23		\n\t"\
		"vmovaps	     (%%r8),%%zmm4			\n\t"\
		"vmovaps	0x440(%%r8),%%zmm0			\n\t"\
		"vmovaps	0x040(%%r8),%%zmm5			\n\t"\
		"vmovaps	0x400(%%r8),%%zmm6			\n\t"\
		"vmulpd		%%zmm16,%%zmm4,%%zmm2		\n\t		vmulpd		%%zmm20,%%zmm10,%%zmm12	\n\t"\
		"vmulpd		%%zmm17,%%zmm0,%%zmm1		\n\t		vmulpd		%%zmm21,%%zmm8 ,%%zmm9 	\n\t"\
		"vmulpd		%%zmm16,%%zmm5,%%zmm3		\n\t		vmulpd		%%zmm20,%%zmm11,%%zmm13	\n\t"\
		"vmulpd		%%zmm17,%%zmm6,%%zmm7		\n\t		vmulpd		%%zmm21,%%zmm14,%%zmm15	\n\t"\
	" vfmadd231pd	%%zmm18,%%zmm5,%%zmm2		\n\t	 vfmadd231pd	%%zmm22,%%zmm11,%%zmm12	\n\t"\
	" vfmadd231pd	%%zmm19,%%zmm6,%%zmm1		\n\t	 vfmadd231pd	%%zmm23,%%zmm14,%%zmm9 	\n\t"\
	"vfnmadd231pd	%%zmm18,%%zmm4,%%zmm3		\n\t	vfnmadd231pd	%%zmm22,%%zmm10,%%zmm13	\n\t"\
	"vfnmadd231pd	%%zmm19,%%zmm0,%%zmm7		\n\t	vfnmadd231pd	%%zmm23,%%zmm8 ,%%zmm15	\n\t"\
		"vmovaps	%%zmm2,     (%%rcx,%%rdi,8)	\n\t	vmovaps	%%zmm12,     (%%rcx,%%r9,8)		\n\t"\
		"vmovaps	%%zmm1,     (%%rdx,%%rdi,8)	\n\t	vmovaps	%%zmm9 ,     (%%rdx,%%r9,8)		\n\t"\
		"vmovaps	%%zmm3,0x040(%%rcx,%%rdi,8)	\n\t	vmovaps	%%zmm13,0x040(%%rcx,%%r9,8)		\n\t"\
		"vmovaps	%%zmm7,0x040(%%rdx,%%rdi,8)	\n\t	vmovaps	%%zmm15,0x040(%%rdx,%%r9,8)		\n\t"\
		"\n\t"\
	/*...Block 3: t04,t14,t24,t34*/					/*...Block 7: t0C,t1C,t2C,t3C*/\
		"movq	%[__isrt2],%%rsi				\n\t	subq	%%rdi,%%r9		\n\t"\
		"movslq	%[__p02],%%rdi					\n\t	addq	%%rdi,%%r9		\n\t"\
		"addq	$0x880,%%r8		/* r24 */		\n\t"\
		"addq	$0x040,%%rsi	/* cc0 */		\n\t"\
		"vmovaps	     (%%r8),%%zmm4			\n\t	vmovaps	0x200(%%r8),%%zmm12				\n\t"\
		"vmovaps	0x040(%%r8),%%zmm5			\n\t	vmovaps	0x240(%%r8),%%zmm13				\n\t"\
		"vmovaps	     (%%rsi),%%zmm2			\n\t"\
		"vmovaps	0x040(%%rsi),%%zmm3			\n\t"\
		"vmovaps	%%zmm4,%%zmm6				\n\t	vmovaps	%%zmm12,%%zmm14					\n\t"\
		"vmovaps	%%zmm5,%%zmm7				\n\t	vmovaps	%%zmm13,%%zmm15					\n\t"\
		"vmulpd	%%zmm2,%%zmm4,%%zmm4			\n\t	vmulpd	%%zmm3,%%zmm12,%%zmm12			\n\t"\
		"vmulpd	%%zmm2,%%zmm5,%%zmm5			\n\t	vmulpd	%%zmm3,%%zmm13,%%zmm13			\n\t"\
		"vmovaps	0x400(%%r8),%%zmm0			\n\t	vmovaps	0x600(%%r8),%%zmm8 				\n\t"\
		"vmovaps	0x440(%%r8),%%zmm1			\n\t	vmovaps	0x640(%%r8),%%zmm9 				\n\t"\
		" vfmadd231pd	%%zmm3,%%zmm7,%%zmm4	\n\t	 vfmadd231pd	%%zmm2,%%zmm15,%%zmm12	\n\t"\
		"vfnmadd231pd	%%zmm3,%%zmm6,%%zmm5	\n\t	vfnmadd231pd	%%zmm2,%%zmm14,%%zmm13	\n\t"\
		"vmovaps	%%zmm0,%%zmm6				\n\t	vmovaps	%%zmm8,%%zmm14					\n\t"\
		"vmovaps	%%zmm1,%%zmm7				\n\t	vmovaps	%%zmm9,%%zmm15					\n\t"\
		"vmulpd	%%zmm3,%%zmm0,%%zmm0			\n\t	vmulpd	%%zmm2 ,%%zmm8 ,%%zmm8 			\n\t"\
		"vmulpd	%%zmm3,%%zmm1,%%zmm1			\n\t	vmulpd	%%zmm2 ,%%zmm9 ,%%zmm9 			\n\t"\
		" vfmadd231pd	%%zmm2,%%zmm7,%%zmm0	\n\t	 vfmadd231pd	%%zmm3,%%zmm15,%%zmm8 	\n\t"\
		"vfnmadd231pd	%%zmm2,%%zmm6,%%zmm1	\n\t	vfnmadd231pd	%%zmm3,%%zmm14,%%zmm9 	\n\t"\
		"vmovaps	%%zmm5,%%zmm7				\n\t	vmovaps	%%zmm13,%%zmm15					\n\t"\
		"vmovaps	%%zmm4,%%zmm6				\n\t	vmovaps	%%zmm12,%%zmm14					\n\t"\
		" vaddpd		%%zmm0,%%zmm4,%%zmm4	\n\t	 vaddpd		%%zmm8 ,%%zmm14,%%zmm14		\n\t"\
		" vaddpd		%%zmm1,%%zmm5,%%zmm5	\n\t	 vaddpd		%%zmm9 ,%%zmm15,%%zmm15		\n\t"\
		" vsubpd		%%zmm0,%%zmm6,%%zmm6	\n\t	 vsubpd		%%zmm8 ,%%zmm12,%%zmm12		\n\t"\
		" vsubpd		%%zmm1,%%zmm7,%%zmm7	\n\t	 vsubpd		%%zmm9 ,%%zmm13,%%zmm13		\n\t"\
		"subq	$0x800,%%r8						\n\t"\
		"subq	$0x040,%%rsi					\n\t"/* isrt2 */\
		"vmovaps	0x400(%%r8),%%zmm2			\n\t	vmovaps	0x600(%%r8),%%zmm8 				\n\t"\
		"vmovaps	0x440(%%r8),%%zmm3			\n\t	vmovaps	0x640(%%r8),%%zmm9 				\n\t"\
		"vmovaps	%%zmm3,%%zmm0				\n\t	vmovaps	%%zmm8 ,%%zmm10					\n\t"\
		"vsubpd	%%zmm2,%%zmm3,%%zmm3			\n\t	vsubpd	%%zmm9 ,%%zmm8,%%zmm8 			\n\t"\
		"vaddpd	%%zmm0,%%zmm2,%%zmm2			\n\t	vaddpd	%%zmm10,%%zmm9,%%zmm9 			\n\t"\
		"vmulpd	%%zmm29,%%zmm2,%%zmm2			\n\t	vmulpd	%%zmm29 ,%%zmm8 	,%%zmm8 	\n\t"/* isrt2 */\
		"vmulpd	%%zmm29,%%zmm3,%%zmm3			\n\t	vmulpd	%%zmm29 ,%%zmm9 	,%%zmm9 	\n\t"\
		"vmovaps	     (%%r8),%%zmm0			\n\t	vmovaps	0x200(%%r8),%%zmm10				\n\t"\
		"vmovaps	0x040(%%r8),%%zmm1			\n\t	vmovaps	0x240(%%r8),%%zmm11				\n\t"\
		"vsubpd	%%zmm2,%%zmm0,%%zmm0			\n\t	vsubpd	%%zmm8,%%zmm10,%%zmm10			\n\t"\
		"vsubpd	%%zmm3,%%zmm1,%%zmm1			\n\t	vsubpd	%%zmm9,%%zmm11,%%zmm11			\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm0,%%zmm2	\n\t	vfmadd132pd	%%zmm31,%%zmm10,%%zmm8 		\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm1,%%zmm3	\n\t	vfmadd132pd	%%zmm31,%%zmm11,%%zmm9 		\n\t"\
		"\n\t"\
		/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c02): */		/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c06): */\
		"addq	$0x5c0,%%rsi	\n\t"/* c02; c06 += 0x200: */\
		" vsubpd		%%zmm4,%%zmm2,%%zmm2	\n\t	 vsubpd		%%zmm12,%%zmm10,%%zmm10		\n\t"\
		" vsubpd		%%zmm7,%%zmm0,%%zmm0	\n\t	 vsubpd		%%zmm15,%%zmm8 ,%%zmm8 		\n\t"\
		" vsubpd		%%zmm5,%%zmm3,%%zmm3	\n\t	 vsubpd		%%zmm13,%%zmm11,%%zmm11		\n\t"\
		" vsubpd		%%zmm6,%%zmm1,%%zmm1	\n\t	 vsubpd		%%zmm14,%%zmm9 ,%%zmm9 		\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm2,%%zmm4		\n\t	vfmadd132pd	%%zmm31,%%zmm10,%%zmm12		\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm0,%%zmm7		\n\t	vfmadd132pd	%%zmm31,%%zmm8 ,%%zmm15		\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm3,%%zmm5		\n\t	vfmadd132pd	%%zmm31,%%zmm11,%%zmm13		\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm1,%%zmm6		\n\t	vfmadd132pd	%%zmm31,%%zmm9 ,%%zmm14		\n\t"\
		"vmovaps		  (%%rsi),%%zmm16		\n\t		vmovaps		0x200(%%rsi),%%zmm20	\n\t"\
		"vmovaps	 0x080(%%rsi),%%zmm17		\n\t		vmovaps		0x280(%%rsi),%%zmm21	\n\t"\
		"vmovaps	 0x040(%%rsi),%%zmm18		\n\t		vmovaps		0x240(%%rsi),%%zmm22	\n\t"\
		"vmovaps	 0x0c0(%%rsi),%%zmm19		\n\t		vmovaps		0x2c0(%%rsi),%%zmm23	\n\t"\
		"vmovaps	%%zmm2,     (%%r8)			\n\t"\
		"vmovaps	%%zmm0,0x440(%%r8)			\n\t"\
		"vmovaps	%%zmm3,0x040(%%r8)			\n\t"\
		"vmovaps	%%zmm6,0x400(%%r8)			\n\t"\
		"vmulpd		%%zmm16,%%zmm4,%%zmm2		\n\t		vmulpd		%%zmm20,%%zmm12,%%zmm24	\n\t"\
		"vmulpd		%%zmm17,%%zmm7,%%zmm0		\n\t		vmulpd		%%zmm21,%%zmm15,%%zmm25	\n\t"\
		"vmulpd		%%zmm16,%%zmm5,%%zmm3		\n\t		vmulpd		%%zmm20,%%zmm13,%%zmm26	\n\t"\
		"vmulpd		%%zmm17,%%zmm1,%%zmm6		\n\t		vmulpd		%%zmm21,%%zmm9 ,%%zmm27	\n\t"\
	" vfmadd231pd	%%zmm18,%%zmm5,%%zmm2		\n\t	 vfmadd231pd	%%zmm22,%%zmm13,%%zmm24	\n\t"\
	" vfmadd231pd	%%zmm19,%%zmm1,%%zmm0		\n\t	 vfmadd231pd	%%zmm23,%%zmm9 ,%%zmm25	\n\t"\
	"vfnmadd231pd	%%zmm18,%%zmm4,%%zmm3		\n\t	vfnmadd231pd	%%zmm22,%%zmm12,%%zmm26	\n\t"\
	"vfnmadd231pd	%%zmm19,%%zmm7,%%zmm6		\n\t	vfnmadd231pd	%%zmm23,%%zmm15,%%zmm27	\n\t"\
		"vmovaps	%%zmm2,     (%%rax,%%rdi,8)	\n\t	vmovaps	%%zmm24,     (%%rax,%%r9,8)		\n\t"\
		"vmovaps	%%zmm0,     (%%rbx,%%rdi,8)	\n\t	vmovaps	%%zmm25,     (%%rbx,%%r9,8)		\n\t"\
		"vmovaps	%%zmm3,0x040(%%rax,%%rdi,8)	\n\t	vmovaps	%%zmm26,0x040(%%rax,%%r9,8)		\n\t"\
		"vmovaps	%%zmm6,0x040(%%rbx,%%rdi,8)	\n\t	vmovaps	%%zmm27,0x040(%%rbx,%%r9,8)		\n\t"\
		"addq	$0x100,%%rsi	\n\t"/* c12; c16+= 0x200: */\
		"vmovaps		  (%%rsi),%%zmm16		\n\t	vmovaps		0x200(%%rsi),%%zmm20		\n\t"\
		"vmovaps	 0x080(%%rsi),%%zmm17		\n\t	vmovaps		0x280(%%rsi),%%zmm21		\n\t"\
		"vmovaps	 0x040(%%rsi),%%zmm18		\n\t	vmovaps		0x240(%%rsi),%%zmm22		\n\t"\
		"vmovaps	 0x0c0(%%rsi),%%zmm19		\n\t	vmovaps		0x2c0(%%rsi),%%zmm23		\n\t"\
		"vmovaps	     (%%r8),%%zmm4			\n\t"\
		"vmovaps	0x440(%%r8),%%zmm0			\n\t"\
		"vmovaps	0x040(%%r8),%%zmm5			\n\t"\
		"vmovaps	0x400(%%r8),%%zmm6			\n\t"\
		"vmulpd		%%zmm16,%%zmm4,%%zmm2		\n\t		vmulpd		%%zmm20,%%zmm10,%%zmm12	\n\t"\
		"vmulpd		%%zmm17,%%zmm0,%%zmm1		\n\t		vmulpd		%%zmm21,%%zmm8 ,%%zmm9 	\n\t"\
		"vmulpd		%%zmm16,%%zmm5,%%zmm3		\n\t		vmulpd		%%zmm20,%%zmm11,%%zmm13	\n\t"\
		"vmulpd		%%zmm17,%%zmm6,%%zmm7		\n\t		vmulpd		%%zmm21,%%zmm14,%%zmm15	\n\t"\
	" vfmadd231pd	%%zmm18,%%zmm5,%%zmm2		\n\t	 vfmadd231pd	%%zmm22,%%zmm11,%%zmm12	\n\t"\
	" vfmadd231pd	%%zmm19,%%zmm6,%%zmm1		\n\t	 vfmadd231pd	%%zmm23,%%zmm14,%%zmm9 	\n\t"\
	"vfnmadd231pd	%%zmm18,%%zmm4,%%zmm3		\n\t	vfnmadd231pd	%%zmm22,%%zmm10,%%zmm13	\n\t"\
	"vfnmadd231pd	%%zmm19,%%zmm0,%%zmm7		\n\t	vfnmadd231pd	%%zmm23,%%zmm8 ,%%zmm15	\n\t"\
		"vmovaps	%%zmm2,     (%%rcx,%%rdi,8)	\n\t	vmovaps	%%zmm12,     (%%rcx,%%r9,8)		\n\t"\
		"vmovaps	%%zmm1,     (%%rdx,%%rdi,8)	\n\t	vmovaps	%%zmm9 ,     (%%rdx,%%r9,8)		\n\t"\
		"vmovaps	%%zmm3,0x040(%%rcx,%%rdi,8)	\n\t	vmovaps	%%zmm13,0x040(%%rcx,%%r9,8)		\n\t"\
		"vmovaps	%%zmm7,0x040(%%rdx,%%rdi,8)	\n\t	vmovaps	%%zmm15,0x040(%%rdx,%%r9,8)		\n\t"\
		"\n\t"\
	/*...Block 4: t06,t16,t26,t36*/					/*...Block 8: t0E,t1E,t2E,t3E*/\
		"movq	%[__isrt2],%%rsi				\n\t	subq	%%rdi,%%r9	\n\t"\
		"movslq	%[__p03],%%rdi					\n\t	addq	%%rdi,%%r9	\n\t"\
		"addq	$0x880,%%r8		/* r26 */		\n\t"\
		"addq	$0x0c0,%%rsi					\n\t"/* cc1; cc3 += 0x080: */\
		"vmovaps	     (%%r8),%%zmm4			\n\t	vmovaps	0x200(%%r8),%%zmm12				\n\t"\
		"vmovaps	0x040(%%r8),%%zmm5			\n\t	vmovaps	0x240(%%r8),%%zmm13				\n\t"\
		"vmovaps	     (%%rsi),%%zmm2			\n\t	vmovaps	0x080(%%rsi),%%zmm10			\n\t"/* cc1,cc3 */\
		"vmovaps	0x040(%%rsi),%%zmm3			\n\t	vmovaps	0x0c0(%%rsi),%%zmm11			\n\t"\
		"vmovaps	%%zmm4,%%zmm6				\n\t	vmovaps	%%zmm12,%%zmm14					\n\t"\
		"vmovaps	%%zmm5,%%zmm7				\n\t	vmovaps	%%zmm13,%%zmm15					\n\t"\
		"vmulpd		%%zmm10,%%zmm4,%%zmm4		\n\t	vmulpd		%%zmm3,%%zmm12,%%zmm12		\n\t"\
		"vmulpd		%%zmm10,%%zmm5,%%zmm5		\n\t	vmulpd		%%zmm3,%%zmm13,%%zmm13		\n\t"\
		"vmovaps	0x400(%%r8),%%zmm0			\n\t	vmovaps	0x600(%%r8),%%zmm8 				\n\t"\
		"vmovaps	0x440(%%r8),%%zmm1			\n\t	vmovaps	0x640(%%r8),%%zmm9 				\n\t"\
		" vfmadd231pd	%%zmm11,%%zmm7,%%zmm4	\n\t	 vfmadd231pd	%%zmm2,%%zmm15,%%zmm12	\n\t"\
		"vfnmadd231pd	%%zmm11,%%zmm6,%%zmm5	\n\t	vfnmadd231pd	%%zmm2,%%zmm14,%%zmm13	\n\t"\
		"vmovaps	%%zmm0,%%zmm6				\n\t	vmovaps	%%zmm8 ,%%zmm14					\n\t"\
		"vmovaps	%%zmm1,%%zmm7				\n\t	vmovaps	%%zmm9 ,%%zmm15					\n\t"\
		"vmulpd		%%zmm3,%%zmm0,%%zmm0		\n\t	vmulpd		%%zmm11,%%zmm8 ,%%zmm8 		\n\t"\
		"vmulpd		%%zmm3,%%zmm1,%%zmm1		\n\t	vmulpd		%%zmm11,%%zmm9 ,%%zmm9 		\n\t"\
		"vfnmadd231pd	%%zmm2,%%zmm7,%%zmm0	\n\t	 vfmadd231pd	%%zmm10,%%zmm15,%%zmm8 	\n\t"\
		" vfmadd231pd	%%zmm2,%%zmm6,%%zmm1	\n\t	vfnmadd231pd	%%zmm10,%%zmm14,%%zmm9 	\n\t"\
		"vmovaps	%%zmm5,%%zmm7				\n\t	vmovaps	%%zmm13,%%zmm15					\n\t"\
		"vmovaps	%%zmm4,%%zmm6				\n\t	vmovaps	%%zmm12,%%zmm14					\n\t"\
		" vaddpd		%%zmm0,%%zmm6,%%zmm6	\n\t	 vsubpd		%%zmm8 ,%%zmm12,%%zmm12		\n\t"\
		" vaddpd		%%zmm1,%%zmm7,%%zmm7	\n\t	 vsubpd		%%zmm9 ,%%zmm13,%%zmm13		\n\t"\
		" vsubpd		%%zmm0,%%zmm4,%%zmm4	\n\t	 vaddpd		%%zmm8 ,%%zmm14,%%zmm14		\n\t"\
		" vsubpd		%%zmm1,%%zmm5,%%zmm5	\n\t	 vaddpd		%%zmm9 ,%%zmm15,%%zmm15		\n\t"\
		"subq	$0x800,%%r8						\n\t"\
		"subq	$0x080,%%rsi					\n\t"/* cc0 */\
		"vmovaps	0x400(%%r8),%%zmm2			\n\t	vmovaps	0x600(%%r8),%%zmm11				\n\t"\
		"vmovaps	0x440(%%r8),%%zmm0			\n\t	vmovaps	0x640(%%r8),%%zmm9 				\n\t"\
		"vmovaps	0x040(%%rsi),%%zmm3			\n\t	vmovaps	0x040(%%rsi),%%zmm10			\n\t"\
		"vmovaps	%%zmm2,%%zmm1				\n\t	vmovaps	%%zmm11,%%zmm8 					\n\t"\
		"vmulpd	%%zmm3 ,%%zmm2,%%zmm2			\n\t	vmulpd	%%zmm10,%%zmm11,%%zmm11			\n\t"\
		"vmulpd	%%zmm0 ,%%zmm3,%%zmm3			\n\t	vmulpd	%%zmm9 ,%%zmm10,%%zmm10			\n\t"\
		" vfmadd231pd	%%zmm28,%%zmm0,%%zmm2	\n\t	vfmadd132pd	%%zmm28,%%zmm11,%%zmm9 		\n\t"/* cc0 */\
		"vfnmadd231pd	%%zmm28,%%zmm1,%%zmm3	\n\t	vfmsub132pd	%%zmm28,%%zmm10,%%zmm8 		\n\t"\
		"vmovaps	     (%%r8),%%zmm0			\n\t	vmovaps	0x200(%%r8),%%zmm10				\n\t"\
		"vmovaps	0x040(%%r8),%%zmm1			\n\t	vmovaps	0x240(%%r8),%%zmm11				\n\t"\
		"vsubpd	%%zmm2,%%zmm0,%%zmm0			\n\t	vsubpd	%%zmm8,%%zmm10,%%zmm10			\n\t"\
		"vsubpd	%%zmm3,%%zmm1,%%zmm1			\n\t	vsubpd	%%zmm9,%%zmm11,%%zmm11			\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm0,%%zmm2	\n\t	vfmadd132pd	%%zmm31,%%zmm10,%%zmm8 		\n\t"\
		"vfmadd132pd	%%zmm31,%%zmm1,%%zmm3	\n\t	vfmadd132pd	%%zmm31,%%zmm11,%%zmm9 		\n\t"\
		"\n\t"\
		/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c03): */	/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c07): */\
		"addq	$0xd80,%%rsi	\n\t"/* c03; c07 += 0x200: */\
		" vsubpd		%%zmm4,%%zmm2,%%zmm2	\n\t	 vsubpd		%%zmm12,%%zmm10,%%zmm10		\n\t"\
		" vsubpd		%%zmm7,%%zmm0,%%zmm0	\n\t	 vsubpd		%%zmm15,%%zmm8 ,%%zmm8 		\n\t"\
		" vsubpd		%%zmm5,%%zmm3,%%zmm3	\n\t	 vsubpd		%%zmm13,%%zmm11,%%zmm11		\n\t"\
		" vsubpd		%%zmm6,%%zmm1,%%zmm1	\n\t	 vsubpd		%%zmm14,%%zmm9 ,%%zmm9 		\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm2,%%zmm4		\n\t	vfmadd132pd	%%zmm31,%%zmm10,%%zmm12		\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm0,%%zmm7		\n\t	vfmadd132pd	%%zmm31,%%zmm8 ,%%zmm15		\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm3,%%zmm5		\n\t	vfmadd132pd	%%zmm31,%%zmm11,%%zmm13		\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm1,%%zmm6		\n\t	vfmadd132pd	%%zmm31,%%zmm9 ,%%zmm14		\n\t"\
		"vmovaps		  (%%rsi),%%zmm16		\n\t		vmovaps		0x200(%%rsi),%%zmm20	\n\t"\
		"vmovaps	 0x080(%%rsi),%%zmm17		\n\t		vmovaps		0x280(%%rsi),%%zmm21	\n\t"\
		"vmovaps	 0x040(%%rsi),%%zmm18		\n\t		vmovaps		0x240(%%rsi),%%zmm22	\n\t"\
		"vmovaps	 0x0c0(%%rsi),%%zmm19		\n\t		vmovaps		0x2c0(%%rsi),%%zmm23	\n\t"\
		"vmovaps	%%zmm2,     (%%r8)			\n\t"\
		"vmovaps	%%zmm0,0x440(%%r8)			\n\t"\
		"vmovaps	%%zmm3,0x040(%%r8)			\n\t"\
		"vmovaps	%%zmm6,0x400(%%r8)			\n\t"\
		"vmulpd		%%zmm16,%%zmm4,%%zmm2		\n\t		vmulpd		%%zmm20,%%zmm12,%%zmm24	\n\t"\
		"vmulpd		%%zmm17,%%zmm7,%%zmm0		\n\t		vmulpd		%%zmm21,%%zmm15,%%zmm25	\n\t"\
		"vmulpd		%%zmm16,%%zmm5,%%zmm3		\n\t		vmulpd		%%zmm20,%%zmm13,%%zmm26	\n\t"\
		"vmulpd		%%zmm17,%%zmm1,%%zmm6		\n\t		vmulpd		%%zmm21,%%zmm9 ,%%zmm27	\n\t"\
	" vfmadd231pd	%%zmm18,%%zmm5,%%zmm2		\n\t	 vfmadd231pd	%%zmm22,%%zmm13,%%zmm24	\n\t"\
	" vfmadd231pd	%%zmm19,%%zmm1,%%zmm0		\n\t	 vfmadd231pd	%%zmm23,%%zmm9 ,%%zmm25	\n\t"\
	"vfnmadd231pd	%%zmm18,%%zmm4,%%zmm3		\n\t	vfnmadd231pd	%%zmm22,%%zmm12,%%zmm26	\n\t"\
	"vfnmadd231pd	%%zmm19,%%zmm7,%%zmm6		\n\t	vfnmadd231pd	%%zmm23,%%zmm15,%%zmm27	\n\t"\
		"vmovaps	%%zmm2,     (%%rax,%%rdi,8)	\n\t	vmovaps	%%zmm24,     (%%rax,%%r9,8)		\n\t"\
		"vmovaps	%%zmm0,     (%%rbx,%%rdi,8)	\n\t	vmovaps	%%zmm25,     (%%rbx,%%r9,8)		\n\t"\
		"vmovaps	%%zmm3,0x040(%%rax,%%rdi,8)	\n\t	vmovaps	%%zmm26,0x040(%%rax,%%r9,8)		\n\t"\
		"vmovaps	%%zmm6,0x040(%%rbx,%%rdi,8)	\n\t	vmovaps	%%zmm27,0x040(%%rbx,%%r9,8)		\n\t"\
		"addq	$0x100,%%rsi	\n\t"/* c13; c17 += 0x200: */\
		"vmovaps		  (%%rsi),%%zmm16		\n\t		vmovaps		0x200(%%rsi),%%zmm20	\n\t"\
		"vmovaps	 0x080(%%rsi),%%zmm17		\n\t		vmovaps		0x280(%%rsi),%%zmm21	\n\t"\
		"vmovaps	 0x040(%%rsi),%%zmm18		\n\t		vmovaps		0x240(%%rsi),%%zmm22	\n\t"\
		"vmovaps	 0x0c0(%%rsi),%%zmm19		\n\t		vmovaps		0x2c0(%%rsi),%%zmm23	\n\t"\
		"vmovaps	     (%%r8),%%zmm4			\n\t"\
		"vmovaps	0x440(%%r8),%%zmm0			\n\t"\
		"vmovaps	0x040(%%r8),%%zmm5			\n\t"\
		"vmovaps	0x400(%%r8),%%zmm6			\n\t"\
		"vmulpd		%%zmm16,%%zmm4,%%zmm2		\n\t		vmulpd		%%zmm20,%%zmm10,%%zmm12	\n\t"\
		"vmulpd		%%zmm17,%%zmm0,%%zmm1		\n\t		vmulpd		%%zmm21,%%zmm8 ,%%zmm9 	\n\t"\
		"vmulpd		%%zmm16,%%zmm5,%%zmm3		\n\t		vmulpd		%%zmm20,%%zmm11,%%zmm13	\n\t"\
		"vmulpd		%%zmm17,%%zmm6,%%zmm7		\n\t		vmulpd		%%zmm21,%%zmm14,%%zmm15	\n\t"\
	" vfmadd231pd	%%zmm18,%%zmm5,%%zmm2		\n\t	 vfmadd231pd	%%zmm22,%%zmm11,%%zmm12	\n\t"\
	" vfmadd231pd	%%zmm19,%%zmm6,%%zmm1		\n\t	 vfmadd231pd	%%zmm23,%%zmm14,%%zmm9 	\n\t"\
	"vfnmadd231pd	%%zmm18,%%zmm4,%%zmm3		\n\t	vfnmadd231pd	%%zmm22,%%zmm10,%%zmm13	\n\t"\
	"vfnmadd231pd	%%zmm19,%%zmm0,%%zmm7		\n\t	vfnmadd231pd	%%zmm23,%%zmm8 ,%%zmm15	\n\t"\
		"vmovaps	%%zmm2,     (%%rcx,%%rdi,8)	\n\t	vmovaps	%%zmm12,     (%%rcx,%%r9,8)		\n\t"\
		"vmovaps	%%zmm1,     (%%rdx,%%rdi,8)	\n\t	vmovaps	%%zmm9 ,     (%%rdx,%%r9,8)		\n\t"\
		"vmovaps	%%zmm3,0x040(%%rcx,%%rdi,8)	\n\t	vmovaps	%%zmm13,0x040(%%rcx,%%r9,8)		\n\t"\
		"vmovaps	%%zmm7,0x040(%%rdx,%%rdi,8)	\n\t	vmovaps	%%zmm15,0x040(%%rdx,%%r9,8)		\n\t"\
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
		 ,[__r00] "m" (Xr00)\
		 ,[__isrt2] "m" (Xisrt2)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r8","r9","r14",\
		"xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15",\
		"xmm16","xmm17","xmm18","xmm19","xmm20","xmm21","xmm22","xmm23","xmm24","xmm25","xmm26","xmm27","xmm28","xmm29","xmm30","xmm31"	/* Clobbered registers */\
	);\
	}

  #endif	// ALL_FMA ?

#elif defined(USE_AVX2)	// FMA-based versions of selected macros in this file for Intel AVX2/FMA3

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
  #ifndef ALL_FMA
	#define ALL_FMA	0	// Jan 2021: See no discernable speed difference between the 2 variants,
  #endif				// so we hope that using the lower-FMA one slightly reduces power consumption

  #if ALL_FMA	// Aggressive-FMA version: Replace most ADD/SUB by FMA-with-one-unity-multiplicand
	#warning Using ALL_FMA version of AVX2 SSE2_RADIX32_DI[F,T}_TWIDDLE macros.

	// In our FMA-ized version, original mix of [488 MOV, 224 FMA, 102 MUL, 214 ADD/SUB]
	//										==> [521 MOV, 385 FMA, 102 MUL,  53 ADD/SUB].
	//
	#define SSE2_RADIX32_DIF_TWIDDLE(Xadd0,Xp01,Xp02,Xp03,Xp04,Xp08,Xp0C,Xp10,Xp18,Xr00)\
	{\
	__asm__ volatile (\
		"movq	%[__add0],%%rax						\n\t"\
	/*...Block 1: */\
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
	#define SSE2_RADIX32_DIT_TWIDDLE(Xadd0,Xp01,Xp02,Xp03,Xp04,Xp05,Xp06,Xp07,Xp08,Xp10,Xr00,Xisrt2)\
	{\
	__asm__ volatile (\
		"movq	%[__add0],%%rax						\n\t"\
	/*...Block 1: */\
		"movq	%[__r00],%%rsi					\n\t"\
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
		"movslq	%[__p10],%%r13					\n\t"\
		"shlq	$3,%%r13						\n\t"\
		"subq	%%rdi,%%r13	/* p10-p8 */		\n\t"\
		"addq	%%r13,%%rax	/* add0+p10 */		\n\t		leaq	(%%r9,%%rax),%%r10	\n\t"\
		"addq	%%r13,%%rbx						\n\t		leaq	(%%r9,%%rbx),%%r11	\n\t"\
		"addq	%%r13,%%rcx						\n\t		leaq	(%%r9,%%rcx),%%r12	\n\t"\
		"addq	%%r13,%%rdx						\n\t		leaq	(%%r9,%%rdx),%%r13	\n\t"\
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
		"leaq	0x200(%%rcx),%%rdx	/* __r10 */			\n\t		addq	%%rsi	,%%r10			\n\t"/* add0 = &a[j1+p4] */\
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
		 ,[__r00] "m" (Xr00)\
		 ,[__isrt2] "m" (Xisrt2)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r8","r9","r10","r11","r12","r13","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
	);\
	}

  #else	// ALL_FMA = False: All FMAs have non-unity multiplicands:
	#warning Using No-unity-multiplicand lower-FMA version of AVX2 SSE2_RADIX32_DI[F|T]_TWIDDLE macros.

	// In our FMA-ized version, original mix of [223 ADD, 223 SUB, 326 MUL] ==> [58+156 ADD, 224 FMA, 102 MUL].
	//
	#define SSE2_RADIX32_DIF_TWIDDLE(Xadd0,Xp01,Xp02,Xp03,Xp04,Xp08,Xp0C,Xp10,Xp18,Xr00)\
	{\
	__asm__ volatile (\
		"movq	%[__add0],%%rax						\n\t"\
	/*...Block 1: */\
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
	#define SSE2_RADIX32_DIT_TWIDDLE(Xadd0,Xp01,Xp02,Xp03,Xp04,Xp05,Xp06,Xp07,Xp08,Xp10,Xr00,Xisrt2)\
	{\
	__asm__ volatile (\
		"movq	%[__add0],%%rax						\n\t"\
	/*...Block 1: */\
		"movq	%[__r00],%%rsi					\n\t"\
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
		"movslq	%[__p10],%%r13					\n\t"\
		"shlq	$3,%%r13						\n\t"\
		"subq	%%rdi,%%r13	/* p10-p8 */		\n\t"\
		"addq	%%r13,%%rax	/* add0+p10 */		\n\t		leaq	(%%r9,%%rax),%%r10	\n\t"\
		"addq	%%r13,%%rbx						\n\t		leaq	(%%r9,%%rbx),%%r11	\n\t"\
		"addq	%%r13,%%rcx						\n\t		leaq	(%%r9,%%rcx),%%r12	\n\t"\
		"addq	%%r13,%%rdx						\n\t		leaq	(%%r9,%%rdx),%%r13	\n\t"\
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
		"leaq	0x200(%%rcx),%%rdx	/* __r10 */			\n\t		addq	%%rsi	,%%r10			\n\t"/* add0 = &a[j1+p4] */\
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
		 ,[__r00] "m" (Xr00)\
		 ,[__isrt2] "m" (Xisrt2)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r8","r9","r10","r11","r12","r13","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
	);\
	}

  #endif	// ALL_FMA ?

#elif defined(USE_AVX)	// AVX and AVX2 both use 256-bit registers

	/*
	For GCC-macro version of this, use that isrt2 + 0x020,0x060,0x0a0 = cc0,cc1,cc3,
	and isrt2 + 0x0e0,0x1e0,0x2e0,0x3e0,0x4e0,0x5e0,0x6e0,0x7e0 = c00,04,02,06,01,05,03,07
	in order to reduce number of args to <= the GCC-allowed maximum of 30:
	*/
	#define SSE2_RADIX32_DIF_TWIDDLE(Xadd0,Xp01,Xp02,Xp03,Xp04,Xp08,Xp0C,Xp10,Xp18,Xr00)\
	{\
	__asm__ volatile (\
		"movq	%[__add0],%%rax						\n\t"\
		/*...Block 1: */\
		"movslq	%[__p08],%%rbx							\n\t"\
		"movslq	%[__p10],%%rcx							\n\t	movslq	%[__p04],%%r9	\n\t"\
		"movslq	%[__p18],%%rdx							\n\t	movq	%[__r00],%%rsi	\n\t"\
		"shlq	$3,%%rbx								\n\t	shlq	$3,%%r9	\n\t"\
		"shlq	$3,%%rcx								\n\t	movq	%%rsi,%%r8	\n\t"\
		"shlq	$3,%%rdx								\n\t	movq	%%r9,%%r10	\n\t	addq	%%rax,%%r10	\n\t"\
		"addq	%%rax,%%rbx								\n\t	movq	%%r9,%%r11	\n\t	addq	%%rbx,%%r11	\n\t"\
		"addq	%%rax,%%rcx								\n\t	movq	%%r9,%%r12	\n\t	addq	%%rcx,%%r12	\n\t"\
		"addq	%%rax,%%rdx								\n\t	movq	%%r9,%%r13	\n\t	addq	%%rdx,%%r13	\n\t"\
		/* SSE2_RADIX4_DIF_4TWIDDLE_A(r00,c10) */				/* SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO_B(r08,c04) */\
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
		/*vmovaps	%%ymm0,0x080(%%rsi)	*/					"	vsubpd	0x120(%%rsi),%%ymm13,%%ymm13	\n\t"\
		/*vmovaps	%%ymm2,0x040(%%rsi)	*/					"	vaddpd	0x100(%%rsi),%%ymm14,%%ymm14	\n\t"\
		/*vmovaps	%%ymm1,0x0a0(%%rsi)	*/					"	vaddpd	0x120(%%rsi),%%ymm15,%%ymm15	\n\t"\
		/*vmovaps	%%ymm3,0x0e0(%%rsi)	*/					"	vsubpd	     %%ymm14,%%ymm8 ,%%ymm8 	\n\t"\
		"vmulpd	(%%r8),%%ymm6,%%ymm6					\n\t	vsubpd	     %%ymm15,%%ymm9 ,%%ymm9 	\n\t"\
		"vmulpd	(%%r8),%%ymm5,%%ymm5					\n\t	vmovaps	%%ymm8 ,0x180(%%rsi)	\n\t"\
		"vmulpd	(%%r8),%%ymm7,%%ymm7					\n\t	vsubpd	%%ymm13,%%ymm10,%%ymm10	\n\t"\
		"vmulpd	(%%r8),%%ymm4,%%ymm4					\n\t"/*	vmovaps	%%ymm9 ,0x1a0(%%rsi)	*/\
		"vaddpd	%%ymm0,%%ymm6,%%ymm6					\n\t	vsubpd	%%ymm12,%%ymm11,%%ymm11	\n\t"\
		"vaddpd	%%ymm2,%%ymm5,%%ymm5					\n\t	vmulpd	(%%r8) ,%%ymm14,%%ymm14	\n\t"\
		"vaddpd	%%ymm1,%%ymm7,%%ymm7					\n\t	vmulpd	(%%r8) ,%%ymm13,%%ymm13	\n\t"\
		"vaddpd	%%ymm3,%%ymm4,%%ymm4					\n\t	vmulpd	(%%r8) ,%%ymm15,%%ymm15	\n\t"\
		/*vmovaps	%%ymm6,     (%%rsi)	*/					"	vmulpd	(%%r8) ,%%ymm12,%%ymm12	\n\t"\
		/*vmovaps	%%ymm5,0x0c0(%%rsi)	*/					"	vaddpd	%%ymm8 ,%%ymm14,%%ymm14	\n\t"\
		/*vmovaps	%%ymm7,0x020(%%rsi)	*/					"	vaddpd	%%ymm10,%%ymm13,%%ymm13	\n\t"\
		/*vmovaps	%%ymm4,0x060(%%rsi)	*/					"	vaddpd	%%ymm9 ,%%ymm15,%%ymm15	\n\t"\
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
		/* Since ymm10 output is first-done of 2nd set, move the computation of ymm2/10 outputs up so can use ymm2 for 2.0 doubling constant below */\
															/*	vmovaps	%%ymm10,0x140(%%rsi)	*/\
															/*	vmovaps	%%ymm13,0x1c0(%%rsi)	*/\
															/*	vmovaps	%%ymm14,0x160(%%rsi)	*/\
															/*	vmovaps	%%ymm15,0x1e0(%%rsi)	*/\
		/* SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r00) */\
		"vmovaps	0x100(%%rsi),%%ymm11		\n\t"\
		"vmovaps	0x120(%%rsi),%%ymm12		\n\t"\
		"						vmulpd	(%%r8),%%ymm10,%%ymm10	\n\t"\
		"vmovaps	0x180(%%rsi),%%ymm8 		\n\t"\
		"vsubpd	%%ymm11,%%ymm6,%%ymm6					\n\t"/*	vsubpd	%%ymm10,%%ymm2,%%ymm2	*/\
		"vsubpd	%%ymm9 ,%%ymm0,%%ymm0					\n\t	vsubpd	%%ymm15,%%ymm5,%%ymm5	\n\t"\
		"						vaddpd	%%ymm2 ,%%ymm10,%%ymm10	\n\t"\
		"						vmovaps	%%ymm2 ,0x140(%%rsi)	\n\t"\
		"vsubpd	%%ymm12,%%ymm7,%%ymm7					\n\t	vsubpd	%%ymm14,%%ymm4,%%ymm4	\n\t"\
		"vmovaps	(%%r8),%%ymm2	\n\t"/* 2.0 ... We could use add-in-place for doubling, but want to load-balance add/mul here. */\
		"vsubpd	%%ymm8 ,%%ymm1 ,%%ymm1 					\n\t	vsubpd	%%ymm13,%%ymm3 ,%%ymm3 	\n\t"\
		"vmulpd	%%ymm2 ,%%ymm11,%%ymm11					\n\t"/*	vaddpd	%%ymm10,%%ymm10,%%ymm10	*/\
		"vmulpd	%%ymm2 ,%%ymm9 ,%%ymm9 					\n\t	vmulpd	%%ymm2 ,%%ymm15,%%ymm15	\n\t"\
		"vmulpd	%%ymm2 ,%%ymm12,%%ymm12					\n\t	vmulpd	%%ymm2 ,%%ymm14,%%ymm14	\n\t"\
		"vmulpd	%%ymm2 ,%%ymm8 ,%%ymm8 					\n\t	vmulpd	%%ymm2 ,%%ymm13,%%ymm13	\n\t"\
		"vaddpd	%%ymm6 ,%%ymm11,%%ymm11					\n\t"/*	vaddpd	%%ymm2 ,%%ymm10,%%ymm10	*/\
		"vaddpd	%%ymm0 ,%%ymm9 ,%%ymm9 					\n\t	vaddpd	%%ymm5 ,%%ymm15,%%ymm15	\n\t"\
		"vaddpd	%%ymm7 ,%%ymm12,%%ymm12					\n\t	vaddpd	%%ymm4 ,%%ymm14,%%ymm14	\n\t"\
		"vaddpd	%%ymm1 ,%%ymm8 ,%%ymm8 					\n\t	vaddpd	%%ymm3 ,%%ymm13,%%ymm13	\n\t"\
		"vmovaps	%%ymm6 ,0x100(%%rsi)				\n\t"/*	vmovaps	%%ymm2 ,0x140(%%rsi)	*/\
		"vmovaps	%%ymm0 ,0x080(%%rsi)				\n\t	vmovaps	%%ymm5 ,0x0c0(%%rsi)	\n\t"\
		"vmovaps	%%ymm7 ,0x120(%%rsi)				\n\t	vmovaps	%%ymm4 ,0x160(%%rsi)	\n\t"\
		"vmovaps	%%ymm1 ,0x1a0(%%rsi)				\n\t	vmovaps	%%ymm3 ,0x1e0(%%rsi)	\n\t"\
		"vmovaps	%%ymm11,     (%%rsi)				\n\t	vmovaps	%%ymm10,0x040(%%rsi)	\n\t"\
		"vmovaps	%%ymm9 ,0x180(%%rsi)				\n\t	vmovaps	%%ymm15,0x1c0(%%rsi)	\n\t"\
		"vmovaps	%%ymm12,0x020(%%rsi)				\n\t	vmovaps	%%ymm14,0x060(%%rsi)	\n\t"\
		"vmovaps	%%ymm8 ,0x0a0(%%rsi)				\n\t	vmovaps	%%ymm13,0x0e0(%%rsi)	\n\t"\
		"\n\t"\
		/***************************************/\
		"\n\t"\
		/*...Block 2: */\
		"movslq	%[__p02],%%rdi							\n\t	movslq	%[__p04],%%r9	\n\t"\
		"shlq	$3,%%rdi								\n\t	shlq	$3,%%r9	\n\t"\
		"addq	%%rdi,%%rax	/* &a[j1+p2] */				\n\t	movq	%%r9,%%r10	\n\t	addq	%%rax,%%r10	\n\t"\
		"addq	%%rdi,%%rbx								\n\t	movq	%%r9,%%r11	\n\t	addq	%%rbx,%%r11	\n\t"\
		"addq	%%rdi,%%rcx								\n\t	movq	%%r9,%%r12	\n\t	addq	%%rcx,%%r12	\n\t"\
		"addq	%%rdi,%%rdx								\n\t	movq	%%r9,%%r13	\n\t	addq	%%rdx,%%r13	\n\t"\
		/* SSE2_RADIX4_DIF_4TWIDDLE_B(r10,c02) */				/* SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO_B(r18,c06) */\
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
		"vsubpd	%%ymm4,%%ymm3,%%ymm3					\n\t"/*	vmovaps	%%ymm9 ,0x1a0(%%rsi)	*/\
		/*vmovaps	%%ymm0,0x080(%%rsi)	*/					"	vsubpd	%%ymm12,%%ymm11,%%ymm11	\n\t"\
		/*vmovaps	%%ymm2,0x040(%%rsi)	*/					"	vmulpd	(%%r8) ,%%ymm14,%%ymm14	\n\t"\
		/*vmovaps	%%ymm1,0x0a0(%%rsi)	*/					"	vmulpd	(%%r8) ,%%ymm13,%%ymm13	\n\t"\
		/*vmovaps	%%ymm3,0x0e0(%%rsi)	*/					"	vmulpd	(%%r8) ,%%ymm15,%%ymm15	\n\t"\
		"vmulpd	(%%r8),%%ymm6,%%ymm6					\n\t	vmulpd	(%%r8) ,%%ymm12,%%ymm12	\n\t"\
		"vmulpd	(%%r8),%%ymm5,%%ymm5					\n\t	vaddpd	%%ymm8 ,%%ymm14,%%ymm14	\n\t"\
		"vmulpd	(%%r8),%%ymm7,%%ymm7					\n\t	vaddpd	%%ymm10,%%ymm13,%%ymm13	\n\t"\
		"vmulpd	(%%r8),%%ymm4,%%ymm4					\n\t	vaddpd	%%ymm9 ,%%ymm15,%%ymm15	\n\t"\
		"vaddpd	%%ymm0,%%ymm6,%%ymm6					\n\t	vaddpd	%%ymm11,%%ymm12,%%ymm12	\n\t"\
		"vaddpd	%%ymm2,%%ymm5,%%ymm5					\n\t	vmovaps	%%ymm14,0x100(%%rsi)	\n\t"\
		"vaddpd	%%ymm1,%%ymm7,%%ymm7					\n\t	vmovaps	%%ymm15,0x120(%%rsi)	\n\t"\
		"vaddpd	%%ymm3,%%ymm4,%%ymm4					\n\t	vmovaps	0x600(%%rsi),%%ymm8 	/* isrt2 */	\n\t"\
		/*vmovaps	%%ymm6,     (%%rsi)	*/					"	vmovaps	%%ymm10,%%ymm14	\n\t"\
		/*vmovaps	%%ymm5,0x0c0(%%rsi)	*/					"	vmovaps	%%ymm13,%%ymm15	\n\t"\
		/*vmovaps	%%ymm7,0x020(%%rsi)	*/					"	vsubpd	%%ymm12,%%ymm10,%%ymm10	\n\t"\
		/*vmovaps	%%ymm4,0x060(%%rsi)	*/					"	vsubpd	%%ymm11,%%ymm13,%%ymm13	\n\t"\
		"														vaddpd	%%ymm12,%%ymm14,%%ymm14	\n\t"\
		"														vaddpd	%%ymm11,%%ymm15,%%ymm15	\n\t"\
		"														vmulpd	%%ymm8 ,%%ymm10,%%ymm10	\n\t"\
		"														vmulpd	%%ymm8 ,%%ymm13,%%ymm13	\n\t"\
		"														vmulpd	%%ymm8 ,%%ymm14,%%ymm14	\n\t"\
		"vsubpd	%%ymm10,%%ymm2,%%ymm2					\n\t	vmulpd	%%ymm8 ,%%ymm15,%%ymm15	\n\t"\
															/*	vmovaps	%%ymm10,0x140(%%rsi)	*/\
															/*	vmovaps	%%ymm13,0x1c0(%%rsi)	*/\
															/*	vmovaps	%%ymm14,0x160(%%rsi)	*/\
															/*	vmovaps	%%ymm15,0x1e0(%%rsi)	*/\
		/* SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r10) */\
		"vmovaps	0x100(%%rsi),%%ymm11		\n\t"\
		"vmovaps	0x120(%%rsi),%%ymm12		\n\t"\
		"						vmulpd	(%%r8),%%ymm10,%%ymm10	\n\t"\
		"vmovaps	0x180(%%rsi),%%ymm8 		\n\t"\
		"vsubpd	%%ymm11,%%ymm6,%%ymm6					\n\t"/*	vsubpd	%%ymm10,%%ymm2,%%ymm2	*/\
		"vsubpd	%%ymm9 ,%%ymm0,%%ymm0					\n\t	vsubpd	%%ymm15,%%ymm5,%%ymm5	\n\t"\
		"						vaddpd	%%ymm2 ,%%ymm10,%%ymm10	\n\t"\
		"						vmovaps	%%ymm2 ,0x140(%%rsi)		\n\t"\
		"vsubpd	%%ymm12,%%ymm7,%%ymm7					\n\t	vsubpd	%%ymm14,%%ymm4,%%ymm4	\n\t"\
		"						vmovaps	(%%r8),%%ymm2	\n\t"/* 2.0 */\
		"vsubpd	%%ymm8 ,%%ymm1 ,%%ymm1 					\n\t	vsubpd	%%ymm13,%%ymm3,%%ymm3	\n\t"\
		"vmulpd	%%ymm2 ,%%ymm11,%%ymm11					\n\t"/*	vaddpd	%%ymm10,%%ymm10,%%ymm10	*/\
		"vmulpd	%%ymm2 ,%%ymm9 ,%%ymm9 					\n\t	vmulpd	%%ymm2 ,%%ymm15,%%ymm15	\n\t"\
		"vmulpd	%%ymm2 ,%%ymm12,%%ymm12					\n\t	vmulpd	%%ymm2 ,%%ymm14,%%ymm14	\n\t"\
		"vmulpd	%%ymm2 ,%%ymm8 ,%%ymm8 					\n\t	vmulpd	%%ymm2 ,%%ymm13,%%ymm13	\n\t"\
		"vaddpd	%%ymm6 ,%%ymm11,%%ymm11					\n\t"/*	vaddpd	%%ymm2 ,%%ymm10,%%ymm10	*/\
		"vaddpd	%%ymm0 ,%%ymm9 ,%%ymm9 					\n\t	vaddpd	%%ymm5 ,%%ymm15,%%ymm15	\n\t"\
		"vaddpd	%%ymm7 ,%%ymm12,%%ymm12					\n\t	vaddpd	%%ymm4 ,%%ymm14,%%ymm14	\n\t"\
		"vaddpd	%%ymm1 ,%%ymm8 ,%%ymm8 					\n\t	vaddpd	%%ymm3 ,%%ymm13,%%ymm13	\n\t"\
		"vmovaps	%%ymm6 ,0x100(%%rsi)				\n\t"/*	vmovaps	%%ymm2 ,0x140(%%rsi)	*/\
		"vmovaps	%%ymm0 ,0x080(%%rsi)				\n\t	vmovaps	%%ymm5 ,0x0c0(%%rsi)	\n\t"\
		"vmovaps	%%ymm7 ,0x120(%%rsi)				\n\t	vmovaps	%%ymm4 ,0x160(%%rsi)	\n\t"\
		"vmovaps	%%ymm1 ,0x1a0(%%rsi)				\n\t	vmovaps	%%ymm3 ,0x1e0(%%rsi)	\n\t"\
		"vmovaps	%%ymm11,     (%%rsi)				\n\t	vmovaps	%%ymm10,0x040(%%rsi)	\n\t"\
		"vmovaps	%%ymm9 ,0x180(%%rsi)				\n\t	vmovaps	%%ymm15,0x1c0(%%rsi)	\n\t"\
		"vmovaps	%%ymm12,0x020(%%rsi)				\n\t	vmovaps	%%ymm14,0x060(%%rsi)	\n\t"\
		"vmovaps	%%ymm8 ,0x0a0(%%rsi)				\n\t	vmovaps	%%ymm13,0x0e0(%%rsi)	\n\t"\
		"\n\t"\
		/***************************************/\
		"\n\t"\
		/*...Block 3: */\
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
		/* SSE2_RADIX4_DIF_4TWIDDLE_B(r20,c01) */				/* SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO_B(r28,c05) */\
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
		"vsubpd	%%ymm4,%%ymm3,%%ymm3					\n\t"/*	vmovaps	%%ymm9 ,0x1a0(%%rsi)	*/\
		/*vmovaps	%%ymm0,0x080(%%rsi)	*/					"	vsubpd	%%ymm12,%%ymm11,%%ymm11	\n\t"\
		/*vmovaps	%%ymm2,0x040(%%rsi)	*/					"	vmulpd	(%%r8) ,%%ymm14,%%ymm14	\n\t"\
		/*vmovaps	%%ymm1,0x0a0(%%rsi)	*/					"	vmulpd	(%%r8) ,%%ymm13,%%ymm13	\n\t"\
		/*vmovaps	%%ymm3,0x0e0(%%rsi)	*/					"	vmulpd	(%%r8) ,%%ymm15,%%ymm15	\n\t"\
		"vmulpd	(%%r8),%%ymm6,%%ymm6					\n\t	vmulpd	(%%r8) ,%%ymm12,%%ymm12	\n\t"\
		"vmulpd	(%%r8),%%ymm5,%%ymm5					\n\t	vaddpd	%%ymm8 ,%%ymm14,%%ymm14	\n\t"\
		"vmulpd	(%%r8),%%ymm7,%%ymm7					\n\t	vaddpd	%%ymm10,%%ymm13,%%ymm13	\n\t"\
		"vmulpd	(%%r8),%%ymm4,%%ymm4					\n\t	vaddpd	%%ymm9 ,%%ymm15,%%ymm15	\n\t"\
		"vaddpd	%%ymm0,%%ymm6,%%ymm6					\n\t	vaddpd	%%ymm11,%%ymm12,%%ymm12	\n\t"\
		"vaddpd	%%ymm2,%%ymm5,%%ymm5					\n\t	vmovaps	%%ymm14,0x100(%%rsi)	\n\t"\
		"vaddpd	%%ymm1,%%ymm7,%%ymm7					\n\t	vmovaps	%%ymm15,0x120(%%rsi)	\n\t"\
		"vaddpd	%%ymm3,%%ymm4,%%ymm4					\n\t	vmovaps	0x400(%%rsi),%%ymm8 	/* isrt2 */	\n\t"\
		/*vmovaps	%%ymm6,     (%%rsi)	*/					"	vmovaps	%%ymm10,%%ymm14	\n\t"\
		/*vmovaps	%%ymm5,0x0c0(%%rsi)	*/					"	vmovaps	%%ymm13,%%ymm15	\n\t"\
		/*vmovaps	%%ymm7,0x020(%%rsi)	*/					"	vsubpd	%%ymm12,%%ymm10,%%ymm10	\n\t"\
		/*vmovaps	%%ymm4,0x060(%%rsi)	*/					"	vsubpd	%%ymm11,%%ymm13,%%ymm13	\n\t"\
		"														vaddpd	%%ymm12,%%ymm14,%%ymm14	\n\t"\
		"														vaddpd	%%ymm11,%%ymm15,%%ymm15	\n\t"\
		"														vmulpd	%%ymm8 ,%%ymm10,%%ymm10	\n\t"\
		"														vmulpd	%%ymm8 ,%%ymm13,%%ymm13	\n\t"\
		"														vmulpd	%%ymm8 ,%%ymm14,%%ymm14	\n\t"\
		"vsubpd	%%ymm10,%%ymm2,%%ymm2					\n\t	vmulpd	%%ymm8 ,%%ymm15,%%ymm15	\n\t"\
															/*	vmovaps	%%ymm10,0x140(%%rsi)	*/\
															/*	vmovaps	%%ymm13,0x1c0(%%rsi)	*/\
															/*	vmovaps	%%ymm14,0x160(%%rsi)	*/\
															/*	vmovaps	%%ymm15,0x1e0(%%rsi)	*/\
		/* SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r20) */\
		"vmovaps	0x100(%%rsi),%%ymm11		\n\t"\
		"vmovaps	0x120(%%rsi),%%ymm12		\n\t"\
		"						vmulpd	(%%r8),%%ymm10,%%ymm10	\n\t"\
		"vmovaps	0x180(%%rsi),%%ymm8 		\n\t"\
		"vsubpd	%%ymm11,%%ymm6,%%ymm6					\n\t"/*	vsubpd	%%ymm10,%%ymm2,%%ymm2	*/\
		"vsubpd	%%ymm9 ,%%ymm0,%%ymm0					\n\t	vsubpd	%%ymm15,%%ymm5,%%ymm5	\n\t"\
		"						vaddpd	%%ymm2 ,%%ymm10,%%ymm10	\n\t"\
		"						vmovaps	%%ymm2 ,0x140(%%rsi)	\n\t"\
		"vsubpd	%%ymm12,%%ymm7,%%ymm7					\n\t	vsubpd	%%ymm14,%%ymm4,%%ymm4	\n\t"\
		"						vmovaps	(%%r8),%%ymm2	\n\t"/* 2.0 */\
		"vsubpd	%%ymm8 ,%%ymm1 ,%%ymm1 					\n\t	vsubpd	%%ymm13,%%ymm3,%%ymm3	\n\t"\
		"vmulpd	%%ymm2 ,%%ymm11,%%ymm11					\n\t"/*	vaddpd	%%ymm10,%%ymm10,%%ymm10	*/\
		"vmulpd	%%ymm2 ,%%ymm9 ,%%ymm9 					\n\t	vmulpd	%%ymm2 ,%%ymm15,%%ymm15	\n\t"\
		"vmulpd	%%ymm2 ,%%ymm12,%%ymm12					\n\t	vmulpd	%%ymm2 ,%%ymm14,%%ymm14	\n\t"\
		"vmulpd	%%ymm2 ,%%ymm8 ,%%ymm8 					\n\t	vmulpd	%%ymm2 ,%%ymm13,%%ymm13	\n\t"\
		"vaddpd	%%ymm6 ,%%ymm11,%%ymm11					\n\t"/*	vaddpd	%%ymm2 ,%%ymm10,%%ymm10	*/\
		"vaddpd	%%ymm0 ,%%ymm9 ,%%ymm9 					\n\t	vaddpd	%%ymm5 ,%%ymm15,%%ymm15	\n\t"\
		"vaddpd	%%ymm7 ,%%ymm12,%%ymm12					\n\t	vaddpd	%%ymm4 ,%%ymm14,%%ymm14	\n\t"\
		"vaddpd	%%ymm1 ,%%ymm8 ,%%ymm8 					\n\t	vaddpd	%%ymm3 ,%%ymm13,%%ymm13	\n\t"\
		"vmovaps	%%ymm6 ,0x100(%%rsi)				\n\t"/*	vmovaps	%%ymm2 ,0x140(%%rsi)	*/\
		"vmovaps	%%ymm0 ,0x080(%%rsi)				\n\t	vmovaps	%%ymm5 ,0x0c0(%%rsi)	\n\t"\
		"vmovaps	%%ymm7 ,0x120(%%rsi)				\n\t	vmovaps	%%ymm4 ,0x160(%%rsi)	\n\t"\
		"vmovaps	%%ymm1 ,0x1a0(%%rsi)				\n\t	vmovaps	%%ymm3 ,0x1e0(%%rsi)	\n\t"\
		"vmovaps	%%ymm11,     (%%rsi)				\n\t	vmovaps	%%ymm10,0x040(%%rsi)	\n\t"\
		"vmovaps	%%ymm9 ,0x180(%%rsi)				\n\t	vmovaps	%%ymm15,0x1c0(%%rsi)	\n\t"\
		"vmovaps	%%ymm12,0x020(%%rsi)				\n\t	vmovaps	%%ymm14,0x060(%%rsi)	\n\t"\
		"vmovaps	%%ymm8 ,0x0a0(%%rsi)				\n\t	vmovaps	%%ymm13,0x0e0(%%rsi)	\n\t"\
		"\n\t"\
		/***************************************/\
		"\n\t"\
		/*...Block 4: */\
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
		/* SSE2_RADIX4_DIF_4TWIDDLE_B(r30,c03) */				/* SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO_B(r38,c07) */\
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
		"vsubpd	%%ymm4,%%ymm3,%%ymm3					\n\t"/*	vmovaps	%%ymm9 ,0x1a0(%%rsi)	*/\
		/*vmovaps	%%ymm0,0x080(%%rsi)	*/					"	vsubpd	%%ymm12,%%ymm11,%%ymm11	\n\t"\
		/*vmovaps	%%ymm2,0x040(%%rsi)	*/					"	vmulpd	(%%r8) ,%%ymm14,%%ymm14	\n\t"\
		/*vmovaps	%%ymm1,0x0a0(%%rsi)	*/					"	vmulpd	(%%r8) ,%%ymm13,%%ymm13	\n\t"\
		/*vmovaps	%%ymm3,0x0e0(%%rsi)	*/					"	vmulpd	(%%r8) ,%%ymm15,%%ymm15	\n\t"\
		"vmulpd	(%%r8),%%ymm6,%%ymm6					\n\t	vmulpd	(%%r8) ,%%ymm12,%%ymm12	\n\t"\
		"vmulpd	(%%r8),%%ymm5,%%ymm5					\n\t	vaddpd	%%ymm8 ,%%ymm14,%%ymm14	\n\t"\
		"vmulpd	(%%r8),%%ymm7,%%ymm7					\n\t	vaddpd	%%ymm10,%%ymm13,%%ymm13	\n\t"\
		"vmulpd	(%%r8),%%ymm4,%%ymm4					\n\t	vaddpd	%%ymm9 ,%%ymm15,%%ymm15	\n\t"\
		"vaddpd	%%ymm0,%%ymm6,%%ymm6					\n\t	vaddpd	%%ymm11,%%ymm12,%%ymm12	\n\t"\
		"vaddpd	%%ymm2,%%ymm5,%%ymm5					\n\t	vmovaps	%%ymm14,0x100(%%rsi)	\n\t"\
		"vaddpd	%%ymm1,%%ymm7,%%ymm7					\n\t	vmovaps	%%ymm15,0x120(%%rsi)	\n\t"\
		"vaddpd	%%ymm3,%%ymm4,%%ymm4					\n\t	vmovaps	0x200(%%rsi),%%ymm8 	/* isrt2 */	\n\t"\
		/*vmovaps	%%ymm6,     (%%rsi)	*/					"	vmovaps	%%ymm10,%%ymm14	\n\t"\
		/*vmovaps	%%ymm5,0x0c0(%%rsi)	*/					"	vmovaps	%%ymm13,%%ymm15	\n\t"\
		/*vmovaps	%%ymm7,0x020(%%rsi)	*/					"	vsubpd	%%ymm12,%%ymm10,%%ymm10	\n\t"\
		/*vmovaps	%%ymm4,0x060(%%rsi)	*/					"	vsubpd	%%ymm11,%%ymm13,%%ymm13	\n\t"\
		"														vaddpd	%%ymm12,%%ymm14,%%ymm14	\n\t"\
		"														vaddpd	%%ymm11,%%ymm15,%%ymm15	\n\t"\
		"														vmulpd	%%ymm8 ,%%ymm10,%%ymm10	\n\t"\
		"														vmulpd	%%ymm8 ,%%ymm13,%%ymm13	\n\t"\
		"														vmulpd	%%ymm8 ,%%ymm14,%%ymm14	\n\t"\
		"vsubpd	%%ymm10,%%ymm2,%%ymm2					\n\t	vmulpd	%%ymm8 ,%%ymm15,%%ymm15	\n\t"\
															/*	vmovaps	%%ymm10,0x140(%%rsi)	*/\
															/*	vmovaps	%%ymm13,0x1c0(%%rsi)	*/\
															/*	vmovaps	%%ymm14,0x160(%%rsi)	*/\
															/*	vmovaps	%%ymm15,0x1e0(%%rsi)	*/\
		/* SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r30) */\
		"vmovaps	0x100(%%rsi),%%ymm11		\n\t"\
		"vmovaps	0x120(%%rsi),%%ymm12		\n\t"\
		"						vmulpd	(%%r8),%%ymm10,%%ymm10	\n\t"\
		"vmovaps	0x180(%%rsi),%%ymm8 		\n\t"\
		"vsubpd	%%ymm11,%%ymm6,%%ymm6					\n\t"/*	vsubpd	%%ymm10,%%ymm2,%%ymm2	*/\
		"vsubpd	%%ymm9 ,%%ymm0,%%ymm0					\n\t	vsubpd	%%ymm15,%%ymm5,%%ymm5	\n\t"\
		"						vaddpd	%%ymm2 ,%%ymm10,%%ymm10	\n\t"\
		"						vmovaps	%%ymm2 ,0x140(%%rsi)	\n\t"\
		"vsubpd	%%ymm12,%%ymm7,%%ymm7					\n\t	vsubpd	%%ymm14,%%ymm4,%%ymm4	\n\t"\
		"						vmovaps	(%%r8),%%ymm2	\n\t"/* 2.0 */\
		"vsubpd	%%ymm8 ,%%ymm1 ,%%ymm1 					\n\t	vsubpd	%%ymm13,%%ymm3,%%ymm3	\n\t"\
		"vmulpd	%%ymm2 ,%%ymm11,%%ymm11					\n\t"/*	vaddpd	%%ymm10,%%ymm10,%%ymm10	*/\
		"vmulpd	%%ymm2 ,%%ymm9 ,%%ymm9 					\n\t	vmulpd	%%ymm2 ,%%ymm15,%%ymm15	\n\t"\
		"vmulpd	%%ymm2 ,%%ymm12,%%ymm12					\n\t	vmulpd	%%ymm2 ,%%ymm14,%%ymm14	\n\t"\
		"vmulpd	%%ymm2 ,%%ymm8 ,%%ymm8 					\n\t	vmulpd	%%ymm2 ,%%ymm13,%%ymm13	\n\t"\
		"vaddpd	%%ymm6 ,%%ymm11,%%ymm11					\n\t"/*	vaddpd	%%ymm2 ,%%ymm10,%%ymm10	*/\
		"vaddpd	%%ymm0 ,%%ymm9 ,%%ymm9 					\n\t	vaddpd	%%ymm5 ,%%ymm15,%%ymm15	\n\t"\
		"vaddpd	%%ymm7 ,%%ymm12,%%ymm12					\n\t	vaddpd	%%ymm4 ,%%ymm14,%%ymm14	\n\t"\
		"vaddpd	%%ymm1 ,%%ymm8 ,%%ymm8 					\n\t	vaddpd	%%ymm3 ,%%ymm13,%%ymm13	\n\t"\
		"vmovaps	%%ymm6 ,0x100(%%rsi)				\n\t"/*	vmovaps	%%ymm2 ,0x140(%%rsi)	*/\
		"vmovaps	%%ymm0 ,0x080(%%rsi)				\n\t	vmovaps	%%ymm5 ,0x0c0(%%rsi)	\n\t"\
		"vmovaps	%%ymm7 ,0x120(%%rsi)				\n\t	vmovaps	%%ymm4 ,0x160(%%rsi)	\n\t"\
		"vmovaps	%%ymm1 ,0x1a0(%%rsi)				\n\t	vmovaps	%%ymm3 ,0x1e0(%%rsi)	\n\t"\
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
		"movslq	%[__p01],%%rbx							\n\t"\
		"movslq	%[__p02],%%rcx							\n\t"\
		"movslq	%[__p03],%%rdx							\n\t"		/*...Block 5: t08,t18,t28,t38	*/\
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
		/*vmovaps	%%ymm2,     (%%rbx)	*/					"		vmulpd	%%ymm2 ,%%ymm14,%%ymm14	\n\t"\
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
		/*...Block 3: t04,t14,t24,t34	*/\
		"addq	$0x080,%%rsi	/* r04 */				\n\t"\
		"subq	%%rax,%%rbx	/* p01 << 3 */				\n\t"\
		"subq	%%rax,%%rcx	/* p02 << 3 */				\n\t"\
		"subq	%%rax,%%rdx	/* p03 << 3 */				\n\t"		/*...Block 7: t0C,t1C,t2C,t3C	*/\
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
		/*...Block 2: t02,t12,t22,t32	*/\
		"subq	$0x040,%%rsi	/* r02 */				\n\t"\
		"movslq	%[__p10],%%rdi							\n\t"\
		"subq	%%rax,%%rbx								\n\t"\
		"subq	%%rax,%%rcx								\n\t"\
		"subq	%%rax,%%rdx								\n\t"		/*...Block 6: t0A,t1A,t2A,t3A	*/\
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
		/*...Block 4: t06,t16,t26,t36	*/\
		"addq	$0x080,%%rsi	/* r06 */	\n\t"\
		"movslq	%[__p18],%%rdi	\n\t"\
		"subq	%%rax,%%rbx	\n\t"\
		"subq	%%rax,%%rcx	\n\t"\
		"subq	%%rax,%%rdx								\n\t"		/*...Block 8: t0E,t1E,t2E,t3E	*/\
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
	#define SSE2_RADIX32_DIT_TWIDDLE(Xadd0,Xp01,Xp02,Xp03,Xp04,Xp05,Xp06,Xp07,Xp08,Xp10,Xr00,Xisrt2)\
	{\
	__asm__ volatile (\
		"movq	%[__add0],%%rax						\n\t"\
	/**************************...Block 1: ****************************/\
		"movq	%[__r00],%%rsi					\n\t"\
		"movslq	%[__p01],%%rbx					\n\t"\
		"movslq	%[__p02],%%rcx					\n\t"\
		"movslq	%[__p03],%%rdx					\n\t		movq	%[__isrt2],%%r8	\n\t"\
		"shlq	$3,%%rbx						\n\t		movslq	%[__p04],%%r9	\n\t"\
		"shlq	$3,%%rcx						\n\t		shlq	$3,%%r9		\n\t	addq	$0x900,%%r8	\n\t"/* two */\
		"shlq	$3,%%rdx						\n\t		movq	%%r9,%%r10	\n\t	addq	%%rax,%%r10	\n\t"\
		"addq	%%rax,%%rbx						\n\t		movq	%%r9,%%r11	\n\t	addq	%%rbx,%%r11	\n\t"\
		"addq	%%rax,%%rcx						\n\t		movq	%%r9,%%r12	\n\t	addq	%%rcx,%%r12	\n\t"\
		"addq	%%rax,%%rdx						\n\t		movq	%%r9,%%r13	\n\t	addq	%%rdx,%%r13	\n\t"\
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
		"vmulpd	(%%r8),%%ymm0,%%ymm0			\n\t		vmulpd	(%%r8) ,%%ymm8 ,%%ymm8 			\n\t"\
		"vmulpd	(%%r8),%%ymm4,%%ymm4			\n\t		vmulpd	(%%r8) ,%%ymm12,%%ymm12			\n\t"\
		"vmulpd	(%%r8),%%ymm1,%%ymm1			\n\t		vmulpd	(%%r8) ,%%ymm9 ,%%ymm9 			\n\t"\
		"vmulpd	(%%r8),%%ymm5,%%ymm5			\n\t		vmulpd	(%%r8) ,%%ymm13,%%ymm13			\n\t"\
		"vaddpd	%%ymm2,%%ymm0,%%ymm0			\n\t		vaddpd	%%ymm10,%%ymm8 ,%%ymm8 			\n\t"\
		"vaddpd	%%ymm6,%%ymm4,%%ymm4			\n\t		vaddpd	%%ymm14,%%ymm12,%%ymm12			\n\t"\
		"vaddpd	%%ymm3,%%ymm1,%%ymm1			\n\t		vaddpd	%%ymm11,%%ymm9 ,%%ymm9 			\n\t"\
		"vaddpd	%%ymm7,%%ymm5,%%ymm5			\n\t		vaddpd	%%ymm15,%%ymm13,%%ymm13			\n\t"\
		/* Finish radix-4 butterfly: */\
		"vsubpd	%%ymm4,%%ymm0,%%ymm0			\n\t		vsubpd	%%ymm12,%%ymm8 ,%%ymm8 			\n\t"\
		"vsubpd	%%ymm7,%%ymm2,%%ymm2			\n\t		vsubpd	%%ymm13,%%ymm9 ,%%ymm9 			\n\t"\
		"vsubpd	%%ymm5,%%ymm1,%%ymm1			\n\t"\
		"vsubpd	%%ymm6,%%ymm3,%%ymm3			\n\t"\
		/*vmovaps	%%ymm0,0x080(%%rsi)*/			"		vsubpd	%%ymm15,%%ymm10,%%ymm10			\n\t"\
		/*vmovaps	%%ymm2,0x0c0(%%rsi)*/			"		vmulpd	(%%r8) ,%%ymm12,%%ymm12			\n\t"\
		/*vmovaps	%%ymm1,0x0a0(%%rsi)*/			"		vsubpd	%%ymm14,%%ymm11,%%ymm11			\n\t"\
		/*vmovaps	%%ymm3,0x060(%%rsi)*/			"		vmulpd	(%%r8) ,%%ymm13,%%ymm13			\n\t"\
		"vmulpd	(%%r8),%%ymm4,%%ymm4			\n\t		vmulpd	(%%r8) ,%%ymm15,%%ymm15			\n\t"\
		"vmulpd	(%%r8),%%ymm7,%%ymm7			\n\t		vaddpd	%%ymm8 ,%%ymm12,%%ymm12			\n\t"\
		"vmulpd	(%%r8),%%ymm5,%%ymm5			\n\t		vmulpd	(%%r8) ,%%ymm14,%%ymm14			\n\t"\
		"vmulpd	(%%r8),%%ymm6,%%ymm6			\n\t		vaddpd	%%ymm9 ,%%ymm13,%%ymm13			\n\t"\
		"vaddpd	%%ymm0,%%ymm4,%%ymm4			\n\t		vaddpd	%%ymm10,%%ymm15,%%ymm15			\n\t"\
		"vaddpd	%%ymm2,%%ymm7,%%ymm7			\n\t"\
		"vaddpd	%%ymm1,%%ymm5,%%ymm5			\n\t		vaddpd	%%ymm11,%%ymm14,%%ymm14			\n\t"\
		"vaddpd	%%ymm3,%%ymm6,%%ymm6			\n\t"\
		/*vmovaps	%%ymm4,     (%%rsi)*/			"		vsubpd	%%ymm15,%%ymm11,%%ymm11			\n\t"\
		/*vmovaps	%%ymm7,0x040(%%rsi)*/			"		vsubpd	%%ymm10,%%ymm14,%%ymm14			\n\t"\
		/*vmovaps	%%ymm5,0x020(%%rsi)*/			"		vmulpd	(%%r8) ,%%ymm15,%%ymm15			\n\t"\
		/*vmovaps	%%ymm6,0x0e0(%%rsi)*/			"		vmulpd	(%%r8) ,%%ymm10,%%ymm10			\n\t"\
		"vsubpd	%%ymm12,%%ymm4	,%%ymm4			\n\t		vaddpd	%%ymm11,%%ymm15,%%ymm15			\n\t"\
		"vsubpd	%%ymm9 ,%%ymm0	,%%ymm0			\n\t		vaddpd	%%ymm14,%%ymm10,%%ymm10			\n\t"\
		"vsubpd	%%ymm13,%%ymm5	,%%ymm5			\n\t		vmulpd	-0x900(%%r8),%%ymm11,%%ymm11	\n\t"/* isrt2 */\
		"vsubpd	%%ymm8 ,%%ymm1	,%%ymm1			\n\t		vmulpd	-0x900(%%r8),%%ymm14,%%ymm14	\n\t"\
		"vmulpd	(%%r8) ,%%ymm12	,%%ymm12		\n\t		vmulpd	-0x900(%%r8),%%ymm15,%%ymm15	\n\t"\
		"vmulpd	(%%r8) ,%%ymm9	,%%ymm9			\n\t		vmulpd	-0x900(%%r8),%%ymm10,%%ymm10	\n\t"\
		"vmulpd	(%%r8) ,%%ymm13	,%%ymm13		\n\t"		/*vmovaps	%%ymm11,0x160(%%rsi)*/\
		"vmulpd	(%%r8) ,%%ymm8	,%%ymm8			\n\t"		/*vmovaps	%%ymm14,0x1e0(%%rsi)*/\
		"vaddpd	%%ymm4 ,%%ymm12	,%%ymm12		\n\t"		/*vmovaps	%%ymm15,0x140(%%rsi)*/\
		"vaddpd	%%ymm0 ,%%ymm9	,%%ymm9			\n\t"		/*vmovaps	%%ymm10,0x1c0(%%rsi)*/\
		/* SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS(r00,r02,r04,r06,r08,r0A,r0C,r0E) */\
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
	/**************************...Block 2: ****************************/\
		"addq	$0x200,%%rsi	/* r10 */		\n\t"\
		"movslq	%[__p08],%%rdi					\n\t"\
		"shlq	$3,%%rdi						\n\t"\
		"addq	%%rdi,%%rax	/* add0+p08 */\n\t		movq	%%r9,%%r10	\n\t	addq	%%rax,%%r10	\n\t"\
		"addq	%%rdi,%%rbx				\n\t		movq	%%r9,%%r11	\n\t	addq	%%rbx,%%r11	\n\t"\
		"addq	%%rdi,%%rcx				\n\t		movq	%%r9,%%r12	\n\t	addq	%%rcx,%%r12	\n\t"\
		"addq	%%rdi,%%rdx				\n\t		movq	%%r9,%%r13	\n\t	addq	%%rdx,%%r13	\n\t"\
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
		"vmulpd	(%%r8),%%ymm0,%%ymm0			\n\t		vmulpd	(%%r8) ,%%ymm8 ,%%ymm8 			\n\t"\
		"vmulpd	(%%r8),%%ymm4,%%ymm4			\n\t		vmulpd	(%%r8) ,%%ymm12,%%ymm12			\n\t"\
		"vmulpd	(%%r8),%%ymm1,%%ymm1			\n\t		vmulpd	(%%r8) ,%%ymm9 ,%%ymm9 			\n\t"\
		"vmulpd	(%%r8),%%ymm5,%%ymm5			\n\t		vmulpd	(%%r8) ,%%ymm13,%%ymm13			\n\t"\
		"vaddpd	%%ymm2,%%ymm0,%%ymm0			\n\t		vaddpd	%%ymm10,%%ymm8 ,%%ymm8 			\n\t"\
		"vaddpd	%%ymm6,%%ymm4,%%ymm4			\n\t		vaddpd	%%ymm14,%%ymm12,%%ymm12			\n\t"\
		"vaddpd	%%ymm3,%%ymm1,%%ymm1			\n\t		vaddpd	%%ymm11,%%ymm9 ,%%ymm9 			\n\t"\
		"vaddpd	%%ymm7,%%ymm5,%%ymm5			\n\t		vaddpd	%%ymm15,%%ymm13,%%ymm13			\n\t"\
		/* Finish radix-4 butterfly: */\
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
		/* SSE2_RADIX8_DIT_COMBINE_RAD4_SUB...*/	"		vsubpd	%%ymm15,%%ymm11,%%ymm11			\n\t"\
		"vsubpd	%%ymm12,%%ymm4	,%%ymm4			\n\t		vsubpd	%%ymm10,%%ymm14,%%ymm14			\n\t"\
		"vsubpd	%%ymm9 ,%%ymm0	,%%ymm0			\n\t		vmulpd	(%%r8) ,%%ymm15,%%ymm15			\n\t"\
		"vsubpd	%%ymm13,%%ymm5	,%%ymm5			\n\t		vmulpd	(%%r8) ,%%ymm10,%%ymm10			\n\t"\
		"vsubpd	%%ymm8 ,%%ymm1	,%%ymm1			\n\t		vaddpd	%%ymm11,%%ymm15,%%ymm15			\n\t"\
		"vmulpd	(%%r8) ,%%ymm12	,%%ymm12		\n\t		vaddpd	%%ymm14,%%ymm10,%%ymm10			\n\t"\
		"vmulpd	(%%r8) ,%%ymm9	,%%ymm9			\n\t		vmulpd	-0x900(%%r8),%%ymm11,%%ymm11			\n\t"\
		"vmulpd	(%%r8) ,%%ymm13	,%%ymm13		\n\t		vmulpd	-0x900(%%r8),%%ymm14,%%ymm14			\n\t"\
		"vmulpd	(%%r8) ,%%ymm8	,%%ymm8			\n\t		vmulpd	-0x900(%%r8),%%ymm15,%%ymm15			\n\t"\
		"vaddpd	%%ymm4 ,%%ymm12	,%%ymm12		\n\t		vmulpd	-0x900(%%r8),%%ymm10,%%ymm10			\n\t"\
		"vaddpd	%%ymm0 ,%%ymm9	,%%ymm9			\n\t"	/*...S(r00,r02,r04,r06,r08,r0A,r0C,r0E) */\
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
	/**************************...Block 3: ****************************/\
		"addq	$0x200,%%rsi	/* r20 */	\n\t"\
		"movslq	%[__p10],%%r13					\n\t"\
		"shlq	$3,%%r13						\n\t"\
		"subq	%%rdi,%%r13	/* p10-p8 */		\n\t"\
		"addq	%%r13,%%rax	/* add0+p10 */		\n\t		leaq	(%%r9,%%rax),%%r10	\n\t"\
		"addq	%%r13,%%rbx						\n\t		leaq	(%%r9,%%rbx),%%r11	\n\t"\
		"addq	%%r13,%%rcx						\n\t		leaq	(%%r9,%%rcx),%%r12	\n\t"\
		"addq	%%r13,%%rdx						\n\t		leaq	(%%r9,%%rdx),%%r13	\n\t"\
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
		"vmulpd	(%%r8),%%ymm0,%%ymm0			\n\t		vmulpd	(%%r8) ,%%ymm8 ,%%ymm8 			\n\t"\
		"vmulpd	(%%r8),%%ymm4,%%ymm4			\n\t		vmulpd	(%%r8) ,%%ymm12,%%ymm12			\n\t"\
		"vmulpd	(%%r8),%%ymm1,%%ymm1			\n\t		vmulpd	(%%r8) ,%%ymm9 ,%%ymm9 			\n\t"\
		"vmulpd	(%%r8),%%ymm5,%%ymm5			\n\t		vmulpd	(%%r8) ,%%ymm13,%%ymm13			\n\t"\
		"vaddpd	%%ymm2,%%ymm0,%%ymm0			\n\t		vaddpd	%%ymm10,%%ymm8 ,%%ymm8 			\n\t"\
		"vaddpd	%%ymm6,%%ymm4,%%ymm4			\n\t		vaddpd	%%ymm14,%%ymm12,%%ymm12			\n\t"\
		"vaddpd	%%ymm3,%%ymm1,%%ymm1			\n\t		vaddpd	%%ymm11,%%ymm9 ,%%ymm9 			\n\t"\
		"vaddpd	%%ymm7,%%ymm5,%%ymm5			\n\t		vaddpd	%%ymm15,%%ymm13,%%ymm13			\n\t"\
		/* Finish radix-4 butterfly: */\
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
		/* SSE2_RADIX8_DIT_COMBINE_RAD4_SUB...*/	"		vsubpd	%%ymm15,%%ymm11,%%ymm11			\n\t"\
		"vsubpd	%%ymm12,%%ymm4	,%%ymm4			\n\t		vsubpd	%%ymm10,%%ymm14,%%ymm14			\n\t"\
		"vsubpd	%%ymm9 ,%%ymm0	,%%ymm0			\n\t		vmulpd	(%%r8) ,%%ymm15,%%ymm15			\n\t"\
		"vsubpd	%%ymm13,%%ymm5	,%%ymm5			\n\t		vmulpd	(%%r8) ,%%ymm10,%%ymm10			\n\t"\
		"vsubpd	%%ymm8 ,%%ymm1	,%%ymm1			\n\t		vaddpd	%%ymm11,%%ymm15,%%ymm15			\n\t"\
		"vmulpd	(%%r8) ,%%ymm12	,%%ymm12		\n\t		vaddpd	%%ymm14,%%ymm10,%%ymm10			\n\t"\
		"vmulpd	(%%r8) ,%%ymm9	,%%ymm9			\n\t		vmulpd	-0x900(%%r8),%%ymm11,%%ymm11			\n\t"\
		"vmulpd	(%%r8) ,%%ymm13	,%%ymm13		\n\t		vmulpd	-0x900(%%r8),%%ymm14,%%ymm14			\n\t"\
		"vmulpd	(%%r8) ,%%ymm8	,%%ymm8			\n\t		vmulpd	-0x900(%%r8),%%ymm15,%%ymm15			\n\t"\
		"vaddpd	%%ymm4 ,%%ymm12	,%%ymm12		\n\t		vmulpd	-0x900(%%r8),%%ymm10,%%ymm10			\n\t"\
		"vaddpd	%%ymm0 ,%%ymm9	,%%ymm9			\n\t"	/*...S(r00,r02,r04,r06,r08,r0A,r0C,r0E) */\
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
	/**************************...Block 4: ****************************/\
		"addq	$0x200,%%rsi	/* r30 */	\n\t"\
		"movslq	%[__p08],%%rdi			\n\t"\
		"shlq	$3,%%rdi				\n\t"\
		"addq	%%rdi,%%rax	/* add0+p18 */\n\t		movq	%%r9,%%r10	\n\t	addq	%%rax,%%r10	\n\t"\
		"addq	%%rdi,%%rbx				\n\t		movq	%%r9,%%r11	\n\t	addq	%%rbx,%%r11	\n\t"\
		"addq	%%rdi,%%rcx				\n\t		movq	%%r9,%%r12	\n\t	addq	%%rcx,%%r12	\n\t"\
		"addq	%%rdi,%%rdx				\n\t		movq	%%r9,%%r13	\n\t	addq	%%rdx,%%r13	\n\t"\
		/* SSE2_RADIX4_DIT_0TWIDDLE_B(r30) */			/* SSE2_RADIX4_DIT_0TWIDDLE_2NDOFTWO_B(r38) */\
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
		/* Finish radix-4 butterfly: */\
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
		/* SSE2_RADIX8_DIT_COMBINE_RAD4_SUB...*/	"		vsubpd	%%ymm15,%%ymm11,%%ymm11			\n\t"\
		"vsubpd	%%ymm12,%%ymm4	,%%ymm4			\n\t		vsubpd	%%ymm10,%%ymm14,%%ymm14			\n\t"\
		"vsubpd	%%ymm9 ,%%ymm0	,%%ymm0			\n\t		vmulpd	(%%r8) ,%%ymm15,%%ymm15			\n\t"\
		"vsubpd	%%ymm13,%%ymm5	,%%ymm5			\n\t		vmulpd	(%%r8) ,%%ymm10,%%ymm10			\n\t"\
		"vsubpd	%%ymm8 ,%%ymm1	,%%ymm1			\n\t		vaddpd	%%ymm11,%%ymm15,%%ymm15			\n\t"\
		"vmulpd	(%%r8) ,%%ymm12	,%%ymm12		\n\t		vaddpd	%%ymm14,%%ymm10,%%ymm10			\n\t"\
		"vmulpd	(%%r8) ,%%ymm9	,%%ymm9			\n\t		vmulpd	-0x900(%%r8),%%ymm11,%%ymm11			\n\t"\
		"vmulpd	(%%r8) ,%%ymm13	,%%ymm13		\n\t		vmulpd	-0x900(%%r8),%%ymm14,%%ymm14			\n\t"\
		"vmulpd	(%%r8) ,%%ymm8	,%%ymm8			\n\t		vmulpd	-0x900(%%r8),%%ymm15,%%ymm15			\n\t"\
		"vaddpd	%%ymm4 ,%%ymm12	,%%ymm12		\n\t		vmulpd	-0x900(%%r8),%%ymm10,%%ymm10			\n\t"\
		"vaddpd	%%ymm0 ,%%ymm9	,%%ymm9			\n\t"	/*...S(r00,r02,r04,r06,r08,r0A,r0C,r0E) */\
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
	/**********************************************************************************/\
	/*...and now do eight radix-4 transforms, including the internal twiddle factors: */\
	/**********************************************************************************/\
		"\n\t"\
		"movslq	%[__p10],%%rdi	\n\t"/* rdi will store copy of p10 throughout */\
		/*...Block 1: t00,t10,t20,t30	*/				/*...Block 5: t08,t18,t28,t38 */\
		"movq	%[__add0],%%rax							\n\t		movslq	%[__p04],%%rsi		\n\t"\
		"movslq	%[__p08],%%rbx							\n\t		movq	%%rax,%%r10			\n\t"\
		"movq	%[__r00],%%rcx							\n\t		shlq	$3,%%rsi			\n\t"\
		"leaq	0x200(%%rcx),%%rdx	/* __r10 */			\n\t		addq	%%rsi	,%%r10	/* add0 = &a[j1+p4] */\n\t"\
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
		/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c00): */		"		vsubpd	%%ymm15,%%ymm13,%%ymm13		\n\t"\
		"vsubpd	%%ymm4,%%ymm2,%%ymm2					\n\t		vsubpd	%%ymm9 ,%%ymm11,%%ymm11		\n\t"\
		"vsubpd	%%ymm7,%%ymm0,%%ymm0					\n\t		vmulpd	(%%r8) ,%%ymm14,%%ymm14		\n\t"\
		"vsubpd	%%ymm5,%%ymm3,%%ymm3					\n\t		vmulpd	(%%r8) ,%%ymm10,%%ymm10		\n\t"\
		"vsubpd	%%ymm6,%%ymm1,%%ymm1					\n\t		vmulpd	(%%r8) ,%%ymm15,%%ymm15		\n\t"\
		"vmulpd	(%%r8),%%ymm4,%%ymm4					\n\t		vmulpd	(%%r8) ,%%ymm9 ,%%ymm9 		\n\t"\
		"vmulpd	(%%r8),%%ymm7,%%ymm7					\n\t		vaddpd	%%ymm12,%%ymm14,%%ymm14		\n\t"\
		"vmulpd	(%%r8),%%ymm5,%%ymm5					\n\t		vaddpd	%%ymm8 ,%%ymm10,%%ymm10		\n\t"\
		"vmulpd	(%%r8),%%ymm6,%%ymm6					\n\t		vaddpd	%%ymm13,%%ymm15,%%ymm15		\n\t"\
		"vaddpd	%%ymm2,%%ymm4,%%ymm4					\n\t		vaddpd	%%ymm11,%%ymm9 ,%%ymm9 		\n\t"\
		"vaddpd	%%ymm0,%%ymm7,%%ymm7					\n\t"		/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c04): */\
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
		/*...Block 2: t02,t12,t22,t32	*/						/*...Block 6: t0A,t1A,t2A,t3A */\
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
		/* cc3 */													/* cc1 */\
		/*vmovaps	0x040(%%rsi),%%ymm2					\n\t		vmovaps	     (%%rsi),%%ymm10	*/\
		/*vmovaps	0x060(%%rsi),%%ymm3					\n\t		vmovaps	0x020(%%rsi),%%ymm11	*/\
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
		/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c01): */			/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c05): */\
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
		/*...Block 3: t04,t14,t24,t34*/							/*...Block 7: t0C,t1C,t2C,t3C*/\
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
		/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c02): */			/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c06): */\
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
		/*...Block 4: t06,t16,t26,t36*/							/*...Block 8: t0E,t1E,t2E,t3E */\
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
		/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c03): */			/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c07): */\
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
		 ,[__r00] "m" (Xr00)\
		 ,[__isrt2] "m" (Xisrt2)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r8","r9","r10","r11","r12","r13","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
	);\
	}

#elif defined(USE_SSE2) && (OS_BITS == 64)

	/*
	[Cf. prefetch.txt re. strategy here] Spread 3 prefetches thru macro, incr. prefetch_addr by 64 bytes each time.
	For GCC-macro version of this, use that isrt2 + 0x010,0x030,0x050 = cc0,cc1,cc3,
	and isrt2 + 0x070,0x0f0,0x170,0x1f0,0x270,0x2f0,0x370,0x3f0 = c00,04,02,06,01,05,03,07
	in order to reduce number of args to <= the GCC-allowed maximum of 30:
	*/
	#define SSE2_RADIX32_DIF_TWIDDLE(Xadd0,Xp01,Xp02,Xp03,Xp04,Xp08,Xp0C,Xp10,Xp18,Xr00)\
	{\
	__asm__ volatile (\
		/*...Block 1: */\
		"movq	%[__add0],%%rax						\n\t"\
		"movslq	%[__p08],%%rbx						\n\t"\
		"movslq	%[__p10],%%rcx						\n\t	movslq	%[__p04],%%r9	\n\t"\
		"movslq	%[__p18],%%rdx						\n\t	movq	%[__r00],%%rsi	\n\t"\
		"movq	%%rsi,%%r8							\n\t"\
		"leaq	(%%rax,%%rbx,8),%%rbx				\n\t"\
		"leaq	(%%rax,%%rcx,8),%%rcx				\n\t"\
		"leaq	(%%rax,%%rdx,8),%%rdx				\n\t"\
		/* SSE2_RADIX4_DIF_4TWIDDLE_A(r00,c10) */			/* SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO_B(r08,c04) */\
		"addq	$0x880,%%r8	/* two */				\n\t	movaps	    (%%rax,%%r9,8),%%xmm8	\n\t"\
		"movaps	    (%%rax),%%xmm0					\n\t	movaps	    (%%rcx,%%r9,8),%%xmm12	\n\t"\
		"movaps	    (%%rcx),%%xmm4					\n\t	movaps	0x10(%%rax,%%r9,8),%%xmm9	\n\t"\
		"movaps	0x10(%%rax),%%xmm1					\n\t	movaps	0x10(%%rcx,%%r9,8),%%xmm13	\n\t"\
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
		"mulpd	0x4e0(%%rsi),%%xmm7					\n\t	movaps	    (%%rdx,%%r9,8),%%xmm12	\n\t"\
		"addpd	%%xmm6,%%xmm5						\n\t	movaps	0x10(%%rdx,%%r9,8),%%xmm13	\n\t"\
		"subpd	%%xmm7,%%xmm4						\n\t	movaps	    (%%rdx,%%r9,8),%%xmm14	\n\t"\
		"movaps	%%xmm5,0x10(%%rsi)					\n\t	movaps	0x10(%%rdx,%%r9,8),%%xmm15	\n\t"\
		"movaps	%%xmm4,    (%%rsi)	/* tmpstr r00 */\n\t	mulpd	0x550(%%rsi),%%xmm12	/* c1C */	\n\t"\
		"movaps	    (%%rbx),%%xmm4					\n\t	mulpd	0x550(%%rsi),%%xmm13	\n\t"\
		"movaps	0x10(%%rbx),%%xmm5					\n\t	mulpd	0x560(%%rsi),%%xmm14	\n\t"\
		"movaps	    (%%rbx),%%xmm6					\n\t	mulpd	0x560(%%rsi),%%xmm15	\n\t"\
		"movaps	0x10(%%rbx),%%xmm7					\n\t	addpd	%%xmm14,%%xmm13	\n\t"\
		"mulpd	0x4b0(%%rsi),%%xmm4	/* c08 */		\n\t	subpd	%%xmm15,%%xmm12	\n\t"\
		"mulpd	0x4b0(%%rsi),%%xmm5					\n\t	movaps	%%xmm13,0x90(%%rsi)	\n\t"\
		"mulpd	0x4c0(%%rsi),%%xmm6					\n\t	movaps	%%xmm12,0x80(%%rsi)	/* tmpstr r08 */\n\t"\
		"mulpd	0x4c0(%%rsi),%%xmm7					\n\t	movaps	    (%%rbx,%%r9,8),%%xmm12	\n\t"\
		"addpd	%%xmm6,%%xmm5						\n\t	movaps	0x10(%%rbx,%%r9,8),%%xmm13	\n\t"\
		"subpd	%%xmm7,%%xmm4						\n\t	movaps	    (%%rbx,%%r9,8),%%xmm14	\n\t"\
		"movaps	%%xmm5,%%xmm7						\n\t	movaps	0x10(%%rbx,%%r9,8),%%xmm15	\n\t"\
		"movaps	%%xmm4,%%xmm6						\n\t	mulpd	0x530(%%rsi),%%xmm12	/* c0C */	\n\t"\
		"subpd	    (%%rsi),%%xmm4	/* r00 */		\n\t	mulpd	0x530(%%rsi),%%xmm13	\n\t"\
		"subpd	0x10(%%rsi),%%xmm5					\n\t	mulpd	0x540(%%rsi),%%xmm14	\n\t"\
		"addpd	    (%%rsi),%%xmm6					\n\t	mulpd	0x540(%%rsi),%%xmm15	\n\t"\
		"addpd	0x10(%%rsi),%%xmm7					\n\t	addpd	%%xmm14,%%xmm13	\n\t"\
		"subpd	%%xmm6,%%xmm0						\n\t	subpd	%%xmm15,%%xmm12	\n\t"\
		"subpd	%%xmm5,%%xmm2						\n\t	movaps	%%xmm13,%%xmm15	\n\t"\
		"subpd	%%xmm7,%%xmm1						\n\t	movaps	%%xmm12,%%xmm14	\n\t"\
		"subpd	%%xmm4,%%xmm3						\n\t	subpd	0x80(%%rsi),%%xmm12	\n\t"\
		/*movaps	%%xmm0,0x40(%%rsi)	*/				"	subpd	0x90(%%rsi),%%xmm13	\n\t"\
		/*movaps	%%xmm2,0x20(%%rsi)	*/				"	addpd	0x80(%%rsi),%%xmm14	\n\t"\
		/*movaps	%%xmm1,0x50(%%rsi)	*/				"	addpd	0x90(%%rsi),%%xmm15	\n\t"\
		/*movaps	%%xmm3,0x70(%%rsi)	*/				"	subpd	%%xmm14,%%xmm8	\n\t"\
		"mulpd	(%%r8),%%xmm6						\n\t	subpd	%%xmm15,%%xmm9	\n\t"\
		"mulpd	(%%r8),%%xmm5						\n\t	movaps	%%xmm8 ,0xc0(%%rsi)	\n\t"\
		"mulpd	(%%r8),%%xmm7						\n\t	subpd	%%xmm13,%%xmm10	\n\t"\
		"mulpd	(%%r8),%%xmm4						\n\t"/*	movaps	%%xmm9 ,0xd0(%%rsi)	*/\
		"addpd	%%xmm0,%%xmm6						\n\t	subpd	%%xmm12,%%xmm11	\n\t"\
		"addpd	%%xmm2,%%xmm5						\n\t	mulpd	(%%r8),%%xmm14	\n\t"\
		"addpd	%%xmm1,%%xmm7						\n\t	mulpd	(%%r8),%%xmm13	\n\t"\
		"addpd	%%xmm3,%%xmm4						\n\t	mulpd	(%%r8),%%xmm15	\n\t"\
		/*movaps	%%xmm6,    (%%rsi)	*/				"	mulpd	(%%r8),%%xmm12	\n\t"\
		/*movaps	%%xmm5,0x60(%%rsi)	*/				"	addpd	%%xmm8 ,%%xmm14	\n\t"\
		/*movaps	%%xmm7,0x10(%%rsi)	*/				"	addpd	%%xmm10,%%xmm13	\n\t"\
		/*movaps	%%xmm4,0x30(%%rsi)	*/				"	addpd	%%xmm9 ,%%xmm15	\n\t"\
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
														/*	movaps	%%xmm10,0xa0(%%rsi)	*/\
														/*	movaps	%%xmm13,0xe0(%%rsi)	*/\
														/*	movaps	%%xmm14,0xb0(%%rsi)	*/\
														/*	movaps	%%xmm15,0xf0(%%rsi)	*/\
		/* SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r00) */\
		"movaps	0x80(%%rsi),%%xmm11		\n\t"\
		"movaps	0x90(%%rsi),%%xmm12		\n\t"\
		"						mulpd	(%%r8),%%xmm10	\n\t"\
		"movaps	0xc0(%%rsi),%%xmm8		\n\t"\
		"subpd	%%xmm11,%%xmm6						\n\t/*	subpd	%%xmm10,%%xmm2	*/\n\t"\
		"subpd	%%xmm9 ,%%xmm0						\n\t	subpd	%%xmm15,%%xmm5	\n\t"\
		"						addpd	%%xmm2 ,%%xmm10	\n\t"\
		"						movaps	%%xmm2 ,0xa0(%%rsi)	\n\t"\
		"subpd	%%xmm12,%%xmm7						\n\t	subpd	%%xmm14,%%xmm4	\n\t"\
		"movaps	(%%r8),%%xmm2	\n\t"/* 2.0 ... We could use add-in-place for doubling, but want to load-balance add/mul here. */\
		"subpd	%%xmm8 ,%%xmm1						\n\t	subpd	%%xmm13,%%xmm3	\n\t"\
		"mulpd	%%xmm2 ,%%xmm11						\n\t"/*	addpd	%%xmm10,%%xmm10	*/\
		"mulpd	%%xmm2 ,%%xmm9						\n\t	mulpd	%%xmm2 ,%%xmm15	\n\t"\
		"mulpd	%%xmm2 ,%%xmm12						\n\t	mulpd	%%xmm2 ,%%xmm14	\n\t"\
		"mulpd	%%xmm2 ,%%xmm8						\n\t	mulpd	%%xmm2 ,%%xmm13	\n\t"\
		"addpd	%%xmm6 ,%%xmm11						\n\t"/*	addpd	%%xmm2 ,%%xmm10	*/\
		"addpd	%%xmm0 ,%%xmm9						\n\t	addpd	%%xmm5 ,%%xmm15	\n\t"\
		"addpd	%%xmm7 ,%%xmm12						\n\t	addpd	%%xmm4 ,%%xmm14	\n\t"\
		"addpd	%%xmm1 ,%%xmm8						\n\t	addpd	%%xmm3 ,%%xmm13	\n\t"\
		"movaps	%%xmm6 ,0x80(%%rsi)					\n\t"/*	movaps	%%xmm2 ,0xa0(%%rsi)	*/\
		"movaps	%%xmm0 ,0x40(%%rsi)					\n\t	movaps	%%xmm5 ,0x60(%%rsi)	\n\t"\
		"movaps	%%xmm7 ,0x90(%%rsi)					\n\t	movaps	%%xmm4 ,0xb0(%%rsi)	\n\t"\
		"movaps	%%xmm1 ,0xd0(%%rsi)					\n\t	movaps	%%xmm3 ,0xf0(%%rsi)	\n\t"\
		"movaps	%%xmm11,    (%%rsi)					\n\t	movaps	%%xmm10,0x20(%%rsi)	\n\t"\
		"movaps	%%xmm9 ,0xc0(%%rsi)					\n\t	movaps	%%xmm15,0xe0(%%rsi)	\n\t"\
		"movaps	%%xmm12,0x10(%%rsi)					\n\t	movaps	%%xmm14,0x30(%%rsi)	\n\t"\
		"movaps	%%xmm8 ,0x50(%%rsi)					\n\t	movaps	%%xmm13,0x70(%%rsi)	\n\t"\
		"\n\t"\
		/***************************************/\
		"\n\t"\
		/*...Block 2: */\
		"movslq	%[__p02],%%rdi						\n\t	addq	%%rdi,%%r9	\n\t"\
		/* SSE2_RADIX4_DIF_4TWIDDLE_B(r10,c02) */			/* SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO_B(r18,c06) */\
		"addq	$0x100,%%rsi	/* r10 */			\n\t	movaps	    (%%rax,%%r9,8),%%xmm8	\n\t"\
		"movaps	    (%%rax,%%rdi,8),%%xmm0			\n\t	movaps	    (%%rcx,%%r9,8),%%xmm12	\n\t"\
		"movaps	    (%%rcx,%%rdi,8),%%xmm4			\n\t	movaps	0x10(%%rax,%%r9,8),%%xmm9	\n\t"\
		"movaps	0x10(%%rax,%%rdi,8),%%xmm1			\n\t	movaps	0x10(%%rcx,%%r9,8),%%xmm13	\n\t"\
		"movaps	0x10(%%rcx,%%rdi,8),%%xmm5			\n\t	movaps	0x4f0(%%rsi),%%xmm14	/* c06 */	\n\t"\
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
		"subpd	%%xmm5,%%xmm3						\n\t	movaps	    (%%rdx,%%r9,8),%%xmm12	\n\t"\
		"movaps	    (%%rdx,%%rdi,8),%%xmm4			\n\t	movaps	0x10(%%rdx,%%r9,8),%%xmm13	\n\t"\
		"movaps	0x10(%%rdx,%%rdi,8),%%xmm5			\n\t	movaps	    (%%rdx,%%r9,8),%%xmm14	\n\t"\
		"movaps	    (%%rdx,%%rdi,8),%%xmm6			\n\t	movaps	0x10(%%rdx,%%r9,8),%%xmm15	\n\t"\
		"movaps	0x10(%%rdx,%%rdi,8),%%xmm7			\n\t	mulpd	0x550(%%rsi),%%xmm12	/* c1E */	\n\t"\
		"mulpd	0x4d0(%%rsi),%%xmm4	/* c1A */		\n\t	mulpd	0x550(%%rsi),%%xmm13	\n\t"\
		"mulpd	0x4d0(%%rsi),%%xmm5					\n\t	mulpd	0x560(%%rsi),%%xmm14	\n\t"\
		"mulpd	0x4e0(%%rsi),%%xmm6					\n\t	mulpd	0x560(%%rsi),%%xmm15	\n\t"\
		"mulpd	0x4e0(%%rsi),%%xmm7					\n\t	addpd	%%xmm14,%%xmm13	\n\t"\
		"addpd	%%xmm6,%%xmm5						\n\t	subpd	%%xmm15,%%xmm12	\n\t"\
		"subpd	%%xmm7,%%xmm4						\n\t	movaps	%%xmm13,0x90(%%rsi)	\n\t"\
		"movaps	%%xmm5,0x010(%%rsi)					\n\t	movaps	%%xmm12,0x80(%%rsi)	/* tmpstr r08 */\n\t"\
		"movaps	%%xmm4,     (%%rsi)					\n\t	movaps	    (%%rbx,%%r9,8),%%xmm12	\n\t"\
		"movaps	    (%%rbx,%%rdi,8),%%xmm4			\n\t	movaps	0x10(%%rbx,%%r9,8),%%xmm13	\n\t"\
		"movaps	0x10(%%rbx,%%rdi,8),%%xmm5			\n\t	movaps	    (%%rbx,%%r9,8),%%xmm14	\n\t"\
		"movaps	    (%%rbx,%%rdi,8),%%xmm6			\n\t	movaps	0x10(%%rbx,%%r9,8),%%xmm15	\n\t"\
		"movaps	0x10(%%rbx,%%rdi,8),%%xmm7			\n\t	mulpd	0x530(%%rsi),%%xmm12	/* c0E */	\n\t"\
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
		"subpd	%%xmm4,%%xmm3						\n\t"/*	movaps	%%xmm9 ,0xd0(%%rsi)	*/\
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
														/*	movaps	%%xmm10,0xa0(%%rsi)	*/\
														/*	movaps	%%xmm13,0xe0(%%rsi)	*/\
														/*	movaps	%%xmm14,0xb0(%%rsi)	*/\
														/*	movaps	%%xmm15,0xf0(%%rsi)	*/\
		/* SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r10) */\
		"movaps	0x80(%%rsi),%%xmm11		\n\t"\
		"movaps	0x90(%%rsi),%%xmm12		\n\t"\
		"						mulpd	(%%r8),%%xmm10	\n\t"\
		"movaps	0xc0(%%rsi),%%xmm8		\n\t"\
		"subpd	%%xmm11,%%xmm6						\n\t"/*	subpd	%%xmm10,%%xmm2	*/\
		"subpd	%%xmm9 ,%%xmm0						\n\t	subpd	%%xmm15,%%xmm5	\n\t"\
		"						addpd	%%xmm2 ,%%xmm10	\n\t"\
		"						movaps	%%xmm2 ,0xa0(%%rsi)	\n\t"\
		"subpd	%%xmm12,%%xmm7						\n\t	subpd	%%xmm14,%%xmm4	\n\t"\
		"						movaps	(%%r8),%%xmm2	\n\t"/* 2.0 */\
		"subpd	%%xmm8 ,%%xmm1						\n\t	subpd	%%xmm13,%%xmm3	\n\t"\
		"mulpd	%%xmm2 ,%%xmm11						\n\t"/*	addpd	%%xmm10,%%xmm10	*/\
		"mulpd	%%xmm2 ,%%xmm9						\n\t	mulpd	%%xmm2 ,%%xmm15	\n\t"\
		"mulpd	%%xmm2 ,%%xmm12						\n\t	mulpd	%%xmm2 ,%%xmm14	\n\t"\
		"mulpd	%%xmm2 ,%%xmm8						\n\t	mulpd	%%xmm2 ,%%xmm13	\n\t"\
		"addpd	%%xmm6 ,%%xmm11						\n\t"/*	addpd	%%xmm2 ,%%xmm10	*/\
		"addpd	%%xmm0 ,%%xmm9						\n\t	addpd	%%xmm5 ,%%xmm15	\n\t"\
		"addpd	%%xmm7 ,%%xmm12						\n\t	addpd	%%xmm4 ,%%xmm14	\n\t"\
		"addpd	%%xmm1 ,%%xmm8						\n\t	addpd	%%xmm3 ,%%xmm13	\n\t"\
		"movaps	%%xmm6 ,0x80(%%rsi)					\n\t"/*	movaps	%%xmm2 ,0xa0(%%rsi)	*/\
		"movaps	%%xmm0 ,0x40(%%rsi)					\n\t	movaps	%%xmm5 ,0x60(%%rsi)	\n\t"\
		"movaps	%%xmm7 ,0x90(%%rsi)					\n\t	movaps	%%xmm4 ,0xb0(%%rsi)	\n\t"\
		"movaps	%%xmm1 ,0xd0(%%rsi)					\n\t	movaps	%%xmm3 ,0xf0(%%rsi)	\n\t"\
		"movaps	%%xmm11,    (%%rsi)					\n\t	movaps	%%xmm10,0x20(%%rsi)	\n\t"\
		"movaps	%%xmm9 ,0xc0(%%rsi)					\n\t	movaps	%%xmm15,0xe0(%%rsi)	\n\t"\
		"movaps	%%xmm12,0x10(%%rsi)					\n\t	movaps	%%xmm14,0x30(%%rsi)	\n\t"\
		"movaps	%%xmm8 ,0x50(%%rsi)					\n\t	movaps	%%xmm13,0x70(%%rsi)	\n\t"\
		"\n\t"\
		/***************************************/\
		"\n\t"\
		/*...Block 3: */\
														"	subq	%%rdi,%%r9	\n\t"\
		"movslq	%[__p01],%%rdi						\n\t	addq	%%rdi,%%r9	\n\t"\
		/* SSE2_RADIX4_DIF_4TWIDDLE_B(r20,c01) */			/* SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO_B(r28,c05) */\
		"addq	$0x100,%%rsi	/* r20 */			\n\t	movaps	    (%%rax,%%r9,8),%%xmm8	\n\t"\
		"movaps	    (%%rax,%%rdi,8),%%xmm0			\n\t	movaps	    (%%rcx,%%r9,8),%%xmm12	\n\t"\
		"movaps	    (%%rcx,%%rdi,8),%%xmm4			\n\t	movaps	0x10(%%rax,%%r9,8),%%xmm9	\n\t"\
		"movaps	0x10(%%rax,%%rdi,8),%%xmm1			\n\t	movaps	0x10(%%rcx,%%r9,8),%%xmm13	\n\t"\
		"movaps	0x10(%%rcx,%%rdi,8),%%xmm5			\n\t	movaps	0x4f0(%%rsi),%%xmm14	/* c05 */	\n\t"\
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
		"subpd	%%xmm5,%%xmm3						\n\t	movaps	    (%%rdx,%%r9,8),%%xmm12	\n\t"\
		"movaps	    (%%rdx,%%rdi,8),%%xmm4			\n\t	movaps	0x10(%%rdx,%%r9,8),%%xmm13	\n\t"\
		"movaps	0x10(%%rdx,%%rdi,8),%%xmm5			\n\t	movaps	    (%%rdx,%%r9,8),%%xmm14	\n\t"\
		"movaps	    (%%rdx,%%rdi,8),%%xmm6			\n\t	movaps	0x10(%%rdx,%%r9,8),%%xmm15	\n\t"\
		"movaps	0x10(%%rdx,%%rdi,8),%%xmm7			\n\t	mulpd	0x550(%%rsi),%%xmm12	/* c1D */	\n\t"\
		"mulpd	0x4d0(%%rsi),%%xmm4	/* c19 */		\n\t	mulpd	0x550(%%rsi),%%xmm13	\n\t"\
		"mulpd	0x4d0(%%rsi),%%xmm5					\n\t	mulpd	0x560(%%rsi),%%xmm14	\n\t"\
		"mulpd	0x4e0(%%rsi),%%xmm6					\n\t	mulpd	0x560(%%rsi),%%xmm15	\n\t"\
		"mulpd	0x4e0(%%rsi),%%xmm7					\n\t	addpd	%%xmm14,%%xmm13	\n\t"\
		"addpd	%%xmm6,%%xmm5						\n\t	subpd	%%xmm15,%%xmm12	\n\t"\
		"subpd	%%xmm7,%%xmm4						\n\t	movaps	%%xmm13,0x90(%%rsi)	\n\t"\
		"movaps	%%xmm5,0x010(%%rsi)					\n\t	movaps	%%xmm12,0x80(%%rsi)	/* tmpstr r08 */\n\t"\
		"movaps	%%xmm4,     (%%rsi)					\n\t	movaps	    (%%rbx,%%r9,8),%%xmm12	\n\t"\
		"movaps	    (%%rbx,%%rdi,8),%%xmm4			\n\t	movaps	0x10(%%rbx,%%r9,8),%%xmm13	\n\t"\
		"movaps	0x10(%%rbx,%%rdi,8),%%xmm5			\n\t	movaps	    (%%rbx,%%r9,8),%%xmm14	\n\t"\
		"movaps	    (%%rbx,%%rdi,8),%%xmm6			\n\t	movaps	0x10(%%rbx,%%r9,8),%%xmm15	\n\t"\
		"movaps	0x10(%%rbx,%%rdi,8),%%xmm7			\n\t	mulpd	0x530(%%rsi),%%xmm12	/* c0D */	\n\t"\
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
		"subpd	%%xmm4,%%xmm3						\n\t"/*	movaps	%%xmm9 ,0xd0(%%rsi)	*/\
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
														/*	movaps	%%xmm10,0xa0(%%rsi)	*/\
														/*	movaps	%%xmm13,0xe0(%%rsi)	*/\
														/*	movaps	%%xmm14,0xb0(%%rsi)	*/\
														/*	movaps	%%xmm15,0xf0(%%rsi)	*/\
		/* SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r20) */\
		"movaps	0x80(%%rsi),%%xmm11		\n\t"\
		"movaps	0x90(%%rsi),%%xmm12		\n\t"\
		"						mulpd	(%%r8),%%xmm10	\n\t"\
		"movaps	0xc0(%%rsi),%%xmm8		\n\t"\
		"subpd	%%xmm11,%%xmm6						\n\t/*	subpd	%%xmm10,%%xmm2	*/\n\t"\
		"subpd	%%xmm9 ,%%xmm0						\n\t	subpd	%%xmm15,%%xmm5	\n\t"\
		"						addpd	%%xmm2 ,%%xmm10	\n\t"\
		"						movaps	%%xmm2 ,0xa0(%%rsi)	\n\t"\
		"subpd	%%xmm12,%%xmm7						\n\t	subpd	%%xmm14,%%xmm4	\n\t"\
		"						movaps	(%%r8),%%xmm2	\n\t"/* 2.0 */\
		"subpd	%%xmm8 ,%%xmm1						\n\t	subpd	%%xmm13,%%xmm3	\n\t"\
		"mulpd	%%xmm2 ,%%xmm11						\n\t"/*	addpd	%%xmm10,%%xmm10	*/\
		"mulpd	%%xmm2 ,%%xmm9						\n\t	mulpd	%%xmm2 ,%%xmm15	\n\t"\
		"mulpd	%%xmm2 ,%%xmm12						\n\t	mulpd	%%xmm2 ,%%xmm14	\n\t"\
		"mulpd	%%xmm2 ,%%xmm8						\n\t	mulpd	%%xmm2 ,%%xmm13	\n\t"\
		"addpd	%%xmm6 ,%%xmm11						\n\t"/*	addpd	%%xmm2 ,%%xmm10	*/\
		"addpd	%%xmm0 ,%%xmm9						\n\t	addpd	%%xmm5 ,%%xmm15	\n\t"\
		"addpd	%%xmm7 ,%%xmm12						\n\t	addpd	%%xmm4 ,%%xmm14	\n\t"\
		"addpd	%%xmm1 ,%%xmm8						\n\t	addpd	%%xmm3 ,%%xmm13	\n\t"\
		"movaps	%%xmm6 ,0x80(%%rsi)					\n\t"/*	movaps	%%xmm2 ,0xa0(%%rsi)	*/\
		"movaps	%%xmm0 ,0x40(%%rsi)					\n\t	movaps	%%xmm5 ,0x60(%%rsi)	\n\t"\
		"movaps	%%xmm7 ,0x90(%%rsi)					\n\t	movaps	%%xmm4 ,0xb0(%%rsi)	\n\t"\
		"movaps	%%xmm1 ,0xd0(%%rsi)					\n\t	movaps	%%xmm3 ,0xf0(%%rsi)	\n\t"\
		"movaps	%%xmm11,    (%%rsi)					\n\t	movaps	%%xmm10,0x20(%%rsi)	\n\t"\
		"movaps	%%xmm9 ,0xc0(%%rsi)					\n\t	movaps	%%xmm15,0xe0(%%rsi)	\n\t"\
		"movaps	%%xmm12,0x10(%%rsi)					\n\t	movaps	%%xmm14,0x30(%%rsi)	\n\t"\
		"movaps	%%xmm8 ,0x50(%%rsi)					\n\t	movaps	%%xmm13,0x70(%%rsi)	\n\t"\
		"\n\t"\
		/***************************************/\
		"\n\t"\
		/*...Block 4: */\
														"	subq	%%rdi,%%r9	\n\t"\
		"movslq	%[__p03],%%rdi						\n\t	addq	%%rdi,%%r9	\n\t"\
		/* SSE2_RADIX4_DIF_4TWIDDLE_B(r30,c03) */			/* SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO_B(r38,c07) */\
		"addq	$0x100,%%rsi	/* r30 */			\n\t	movaps	    (%%rax,%%r9,8),%%xmm8	\n\t"\
		"movaps	    (%%rax,%%rdi,8),%%xmm0			\n\t	movaps	    (%%rcx,%%r9,8),%%xmm12	\n\t"\
		"movaps	    (%%rcx,%%rdi,8),%%xmm4			\n\t	movaps	0x10(%%rax,%%r9,8),%%xmm9	\n\t"\
		"movaps	0x10(%%rax,%%rdi,8),%%xmm1			\n\t	movaps	0x10(%%rcx,%%r9,8),%%xmm13	\n\t"\
		"movaps	0x10(%%rcx,%%rdi,8),%%xmm5			\n\t	movaps	0x4f0(%%rsi),%%xmm14	/* c07 */	\n\t"\
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
		"subpd	%%xmm5,%%xmm3						\n\t	movaps	    (%%rdx,%%r9,8),%%xmm12	\n\t"\
		"movaps	    (%%rdx,%%rdi,8),%%xmm4			\n\t	movaps	0x10(%%rdx,%%r9,8),%%xmm13	\n\t"\
		"movaps	0x10(%%rdx,%%rdi,8),%%xmm5			\n\t	movaps	    (%%rdx,%%r9,8),%%xmm14	\n\t"\
		"movaps	    (%%rdx,%%rdi,8),%%xmm6			\n\t	movaps	0x10(%%rdx,%%r9,8),%%xmm15	\n\t"\
		"movaps	0x10(%%rdx,%%rdi,8),%%xmm7			\n\t	mulpd	0x550(%%rsi),%%xmm12	/* c1F */	\n\t"\
		"mulpd	0x4d0(%%rsi),%%xmm4	/* c1B */		\n\t	mulpd	0x550(%%rsi),%%xmm13	\n\t"\
		"mulpd	0x4d0(%%rsi),%%xmm5					\n\t	mulpd	0x560(%%rsi),%%xmm14	\n\t"\
		"mulpd	0x4e0(%%rsi),%%xmm6					\n\t	mulpd	0x560(%%rsi),%%xmm15	\n\t"\
		"mulpd	0x4e0(%%rsi),%%xmm7					\n\t	addpd	%%xmm14,%%xmm13	\n\t"\
		"addpd	%%xmm6,%%xmm5						\n\t	subpd	%%xmm15,%%xmm12	\n\t"\
		"subpd	%%xmm7,%%xmm4						\n\t	movaps	%%xmm13,0x90(%%rsi)	\n\t"\
		"movaps	%%xmm5,0x010(%%rsi)					\n\t	movaps	%%xmm12,0x80(%%rsi)	/* tmpstr r08 */\n\t"\
		"movaps	%%xmm4,     (%%rsi)					\n\t	movaps	    (%%rbx,%%r9,8),%%xmm12	\n\t"\
		"movaps	    (%%rbx,%%rdi,8),%%xmm4			\n\t	movaps	0x10(%%rbx,%%r9,8),%%xmm13	\n\t"\
		"movaps	0x10(%%rbx,%%rdi,8),%%xmm5			\n\t	movaps	    (%%rbx,%%r9,8),%%xmm14	\n\t"\
		"movaps	    (%%rbx,%%rdi,8),%%xmm6			\n\t	movaps	0x10(%%rbx,%%r9,8),%%xmm15	\n\t"\
		"movaps	0x10(%%rbx,%%rdi,8),%%xmm7			\n\t	mulpd	0x530(%%rsi),%%xmm12	/* c0F */	\n\t"\
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
		"subpd	%%xmm4,%%xmm3						\n\t"/*	movaps	%%xmm9 ,0xd0(%%rsi)	*/\
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
														/*	movaps	%%xmm10,0xa0(%%rsi)	*/\
														/*	movaps	%%xmm13,0xe0(%%rsi)	*/\
														/*	movaps	%%xmm14,0xb0(%%rsi)	*/\
														/*	movaps	%%xmm15,0xf0(%%rsi)	*/\
		/* SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r30) */\
		"movaps	0x80(%%rsi),%%xmm11		\n\t"\
		"movaps	0x90(%%rsi),%%xmm12		\n\t"\
		"						mulpd	(%%r8),%%xmm10	\n\t"\
		"movaps	0xc0(%%rsi),%%xmm8		\n\t"\
		"subpd	%%xmm11,%%xmm6						\n\t"/*	subpd	%%xmm10,%%xmm2	*/\
		"subpd	%%xmm9 ,%%xmm0						\n\t	subpd	%%xmm15,%%xmm5	\n\t"\
		"						addpd	%%xmm2 ,%%xmm10	\n\t"\
		"						movaps	%%xmm2 ,0xa0(%%rsi)	\n\t"\
		"subpd	%%xmm12,%%xmm7						\n\t	subpd	%%xmm14,%%xmm4	\n\t"\
		"						movaps	(%%r8),%%xmm2	\n\t"/* 2.0 */\
		"subpd	%%xmm8 ,%%xmm1						\n\t	subpd	%%xmm13,%%xmm3	\n\t"\
		"mulpd	%%xmm2 ,%%xmm11						\n\t"/*	addpd	%%xmm10,%%xmm10	*/\
		"mulpd	%%xmm2 ,%%xmm9						\n\t	mulpd	%%xmm2 ,%%xmm15	\n\t"\
		"mulpd	%%xmm2 ,%%xmm12						\n\t	mulpd	%%xmm2 ,%%xmm14	\n\t"\
		"mulpd	%%xmm2 ,%%xmm8						\n\t	mulpd	%%xmm2 ,%%xmm13	\n\t"\
		"addpd	%%xmm6 ,%%xmm11						\n\t"/*\n\t	addpd	%%xmm2 ,%%xmm10	*/\
		"addpd	%%xmm0 ,%%xmm9						\n\t	addpd	%%xmm5 ,%%xmm15	\n\t"\
		"addpd	%%xmm7 ,%%xmm12						\n\t	addpd	%%xmm4 ,%%xmm14	\n\t"\
		"addpd	%%xmm1 ,%%xmm8						\n\t	addpd	%%xmm3 ,%%xmm13	\n\t"\
		"movaps	%%xmm6 ,0x80(%%rsi)					\n\t"/*	movaps	%%xmm2 ,0xa0(%%rsi)	*/\
		"movaps	%%xmm0 ,0x40(%%rsi)					\n\t	movaps	%%xmm5 ,0x60(%%rsi)	\n\t"\
		"movaps	%%xmm7 ,0x90(%%rsi)					\n\t	movaps	%%xmm4 ,0xb0(%%rsi)	\n\t"\
		"movaps	%%xmm1 ,0xd0(%%rsi)					\n\t	movaps	%%xmm3 ,0xf0(%%rsi)	\n\t"\
		"movaps	%%xmm11,    (%%rsi)					\n\t	movaps	%%xmm10,0x20(%%rsi)	\n\t"\
		"movaps	%%xmm9 ,0xc0(%%rsi)					\n\t	movaps	%%xmm15,0xe0(%%rsi)	\n\t"\
		"movaps	%%xmm12,0x10(%%rsi)					\n\t	movaps	%%xmm14,0x30(%%rsi)	\n\t"\
		"movaps	%%xmm8 ,0x50(%%rsi)					\n\t	movaps	%%xmm13,0x70(%%rsi)	\n\t"\
		"\n\t"\
	/**********************************************************************************/\
	/*...and now do eight radix-4 transforms, including the internal twiddle factors: */\
	/**********************************************************************************/\
		/*...Block 1: t00,t10,t20,t30	*/\
		"movq	%[__add0],%%rax	/* &a[j1] */		\n\t"\
		"movslq	%[__p01],%%rbx						\n\t"\
		"movslq	%[__p02],%%rcx						\n\t"\
		"movslq	%[__p03],%%rdx						\n\t"		/*...Block 5: t08,t18,t28,t38	*/\
		"leaq	(%%rax,%%rbx,8),%%rbx				\n\t		subq	%%rdi,%%r9	\n\t"\
		"leaq	(%%rax,%%rcx,8),%%rcx				\n\t"\
		"leaq	(%%rax,%%rdx,8),%%rdx				\n\t"\
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
		/*movaps	%%xmm2,    (%%rbx)	*/				"		mulpd	%%xmm2,%%xmm14	\n\t"\
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
		"											\n\t		movaps	%%xmm8 ,    (%%rbx,%%r9,8)	\n\t"\
		"											\n\t		movaps	%%xmm11,    (%%rcx,%%r9,8)	\n\t"\
		"											\n\t		movaps	%%xmm10,0x10(%%rbx,%%r9,8)	\n\t"\
		"											\n\t		movaps	%%xmm9 ,0x10(%%rdx,%%r9,8)	\n\t"\
		"											\n\t		addpd	%%xmm8 ,%%xmm12		\n\t"\
		"											\n\t		addpd	%%xmm11,%%xmm15		\n\t"\
		"											\n\t		addpd	%%xmm10,%%xmm13		\n\t"\
		"											\n\t		addpd	%%xmm9 ,%%xmm14		\n\t"\
		"											\n\t		movaps	%%xmm12,    (%%rax,%%r9,8)	\n\t"\
		"											\n\t		movaps	%%xmm15,    (%%rdx,%%r9,8)	\n\t"\
		"											\n\t		movaps	%%xmm13,0x10(%%rax,%%r9,8)	\n\t"\
		"											\n\t		movaps	%%xmm14,0x10(%%rcx,%%r9,8)	\n\t"\
		"\n\t"\
		/*...Block 3: t04,t14,t24,t34	*/\
		"addq	$0x40,%%rsi	/* r04 */				\n\t"		/*...Block 7: t0C,t1C,t2C,t3C	*/\
		"movslq	%[__p08],%%rdi						\n\t		addq	%%rdi,%%r9	\n\t"\
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
		"movaps	%%xmm2,    (%%rbx,%%rdi,8)			\n\t		movaps	%%xmm8 ,    (%%rbx,%%r9,8)	\n\t"\
		"movaps	%%xmm0,    (%%rcx,%%rdi,8)			\n\t		movaps	%%xmm10,    (%%rcx,%%r9,8)	\n\t"\
		"movaps	%%xmm3,0x10(%%rbx,%%rdi,8)			\n\t		movaps	%%xmm9 ,0x10(%%rbx,%%r9,8)	\n\t"\
		"movaps	%%xmm1,0x10(%%rdx,%%rdi,8)			\n\t		movaps	%%xmm11,0x10(%%rdx,%%r9,8)	\n\t"\
		"addpd	%%xmm2,%%xmm6						\n\t		addpd	%%xmm8 ,%%xmm12		\n\t"\
		"addpd	%%xmm0,%%xmm5						\n\t		addpd	%%xmm10,%%xmm15		\n\t"\
		"addpd	%%xmm3,%%xmm7						\n\t		addpd	%%xmm9 ,%%xmm13		\n\t"\
		"addpd	%%xmm1,%%xmm4						\n\t		addpd	%%xmm11,%%xmm14		\n\t"\
		"movaps	%%xmm6,    (%%rax,%%rdi,8)			\n\t		movaps	%%xmm12,    (%%rax,%%r9,8)	\n\t"\
		"movaps	%%xmm5,    (%%rdx,%%rdi,8)			\n\t		movaps	%%xmm15,    (%%rdx,%%r9,8)	\n\t"\
		"movaps	%%xmm7,0x10(%%rax,%%rdi,8)			\n\t		movaps	%%xmm13,0x10(%%rax,%%r9,8)	\n\t"\
		"movaps	%%xmm4,0x10(%%rcx,%%rdi,8)			\n\t		movaps	%%xmm14,0x10(%%rcx,%%r9,8)	\n\t"\
		"\n\t"\
		/*...Block 2: t02,t12,t22,t32	*/\
		"subq	$0x20,%%rsi	/* r02 */				\n\t		subq	%%rdi,%%r9	\n\t"\
		"movslq	%[__p10],%%rdi						\n\t		addq	%%rdi,%%r9	\n\t"\
																/*...Block 6: t0A,t1A,t2A,t3A	*/\
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
		"movaps	%%xmm2,    (%%rbx,%%rdi,8)			\n\t		movaps	%%xmm8 ,    (%%rbx,%%r9,8)	\n\t"\
		"movaps	%%xmm0,    (%%rcx,%%rdi,8)			\n\t		movaps	%%xmm10,    (%%rcx,%%r9,8)	\n\t"\
		"movaps	%%xmm3,0x10(%%rbx,%%rdi,8)			\n\t		movaps	%%xmm9 ,0x10(%%rbx,%%r9,8)	\n\t"\
		"movaps	%%xmm1,0x10(%%rdx,%%rdi,8)			\n\t		movaps	%%xmm11,0x10(%%rdx,%%r9,8)	\n\t"\
		"addpd	%%xmm2,%%xmm6						\n\t		addpd	%%xmm8 ,%%xmm12		\n\t"\
		"addpd	%%xmm0,%%xmm5						\n\t		addpd	%%xmm10,%%xmm15		\n\t"\
		"addpd	%%xmm3,%%xmm7						\n\t		addpd	%%xmm9 ,%%xmm13		\n\t"\
		"addpd	%%xmm1,%%xmm4						\n\t		addpd	%%xmm11,%%xmm14		\n\t"\
		"movaps	%%xmm6,    (%%rax,%%rdi,8)			\n\t		movaps	%%xmm12,    (%%rax,%%r9,8)	\n\t"\
		"movaps	%%xmm5,    (%%rdx,%%rdi,8)			\n\t		movaps	%%xmm15,    (%%rdx,%%r9,8)	\n\t"\
		"movaps	%%xmm7,0x10(%%rax,%%rdi,8)			\n\t		movaps	%%xmm13,0x10(%%rax,%%r9,8)	\n\t"\
		"movaps	%%xmm4,0x10(%%rcx,%%rdi,8)			\n\t		movaps	%%xmm14,0x10(%%rcx,%%r9,8)	\n\t"\
		"\n\t"\
		/*...Block 4: t06,t16,t26,t36	*/\
		"addq	$0x40,%%rsi	/* r06 */				\n\t		subq	%%rdi,%%r9	\n\t"\
		"movslq	%[__p18],%%rdi						\n\t		addq	%%rdi,%%r9	\n\t"\
																/*...Block 8: t0E,t1E,t2E,t3E	*/\
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
		"movaps	%%xmm2,    (%%rbx,%%rdi,8)			\n\t		movaps	%%xmm8 ,    (%%rbx,%%r9,8)	\n\t"\
		"movaps	%%xmm0,    (%%rcx,%%rdi,8)			\n\t		movaps	%%xmm10,    (%%rcx,%%r9,8)	\n\t"\
		"movaps	%%xmm3,0x10(%%rbx,%%rdi,8)			\n\t		movaps	%%xmm9 ,0x10(%%rbx,%%r9,8)	\n\t"\
		"movaps	%%xmm1,0x10(%%rdx,%%rdi,8)			\n\t		movaps	%%xmm11,0x10(%%rdx,%%r9,8)	\n\t"\
		"addpd	%%xmm2,	%%xmm4						\n\t		addpd	%%xmm8 ,%%xmm12		\n\t"\
		"addpd	%%xmm0,	%%xmm7						\n\t		addpd	%%xmm10,%%xmm15		\n\t"\
		"addpd	%%xmm3,	%%xmm5						\n\t		addpd	%%xmm9 ,%%xmm13		\n\t"\
		"addpd	%%xmm1,	%%xmm6						\n\t		addpd	%%xmm11,%%xmm14		\n\t"\
		"movaps	%%xmm4,    (%%rax,%%rdi,8)			\n\t		movaps	%%xmm12,    (%%rax,%%r9,8)	\n\t"\
		"movaps	%%xmm7,    (%%rdx,%%rdi,8)			\n\t		movaps	%%xmm15,    (%%rdx,%%r9,8)	\n\t"\
		"movaps	%%xmm5,0x10(%%rax,%%rdi,8)			\n\t		movaps	%%xmm13,0x10(%%rax,%%r9,8)	\n\t"\
		"movaps	%%xmm6,0x10(%%rcx,%%rdi,8)			\n\t		movaps	%%xmm14,0x10(%%rcx,%%r9,8)	\n\t"\
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
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r8","r9","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"		/* Clobbered registers */\
	);\
	}

	/*
	[Cf. prefetch.txt re. strategy here] Spread 3 prefetches thru macro, incr. prefetch_addr by 64 bytes each time.
	For GCC-macro version of this, use that isrt2 + 0x010,0x030,0x050 = cc0,cc1,cc3,
	and isrt2 + 0x070,0x0f0,0x170,0x1f0,0x270,0x2f0,0x370,0x3f0 = c00,04,02,06,01,05,03,07
	in order to reduce number of args to <= the GCC-allowed maximum of 30:
	*/
	#define SSE2_RADIX32_DIT_TWIDDLE(Xadd0,Xp01,Xp02,Xp03,Xp04,Xp05,Xp06,Xp07,Xp08,Xp10,Xr00,Xisrt2)\
	{\
	__asm__ volatile (\
	/**************************...Block 1: ****************************/\
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
	/**************************...Block 2: ****************************/\
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
	/**************************...Block 3: ****************************/\
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
	/**************************...Block 4: ****************************/\
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
	/**********************************************************************************/\
	/*...and now do eight radix-4 transforms, including the internal twiddle factors: */\
	/**********************************************************************************/\
		"movslq	%[__p10],%%rdi	/* edi will store copy of p10 throughout */\n\t"\
		/*...Block 1: t00,t10,t20,t30	*/				/*...Block 5: t08,t18,t28,t38*/\
		"movq	%[__add0],%%rax							\n\t		movslq	%[__p04],%%rsi		\n\t"\
		"movslq	%[__p08],%%rbx							\n\t		movq	%%rax,%%r10			\n\t"\
		"movq	%[__r00],%%rcx							\n\t		shlq	$3,%%rsi			\n\t"\
		"leaq	0x100(%%rcx),%%rdx	/* __r10 */			\n\t		addq	%%rsi	,%%r10	/* add0 = &a[j1+p4] */\n\t"\
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
		/*...Block 2: t02,t12,t22,t32	*/				/*...Block 6: t0A,t1A,t2A,t3A*/\
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
		/*...Block 3: t04,t14,t24,t34*/				/*...Block 7: t0C,t1C,t2C,t3C*/\
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
		/*...Block 4: t06,t16,t26,t36*/						/*...Block 8: t0E,t1E,t2E,t3E*/\
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
		 ,[__r00] "m" (Xr00)\
		 ,[__isrt2] "m" (Xisrt2)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r8","r9","r10","r11","r12","r13","r14","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
	);\
	}

#elif defined(USE_SSE2) && (OS_BITS == 32)

	#error 32-bit OSes no longer supported for SIMD builds!

#endif	// AVX / SSE2 toggle

#endif	/* radix32_dif_dit_pass_asm_h_included */

