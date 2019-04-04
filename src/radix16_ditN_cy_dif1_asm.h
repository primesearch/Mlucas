/*******************************************************************************
*                                                                              *
*   (C) 1997-2018 by Ernst W. Mayer.                                           *
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
#ifndef radix16_ditN_cy_dif1_gcc_h_included
#define radix16_ditN_cy_dif1_gcc_h_included

#ifdef USE_ARM_V8_SIMD

	// Vector-opcount: 32 LDP, 32 STP, 68 ADD, 24 MUL, which is [24,60,72] fewer [LDP,ADD,MUL] than with-twiddles version.
	#define SSE2_RADIX16_DIT_NOTWIDDLE(Xadd0,Xp1,Xp2,Xp3,Xp4,Xp8,Xr1,Xr3,Xr5,Xr7,Xr9,Xr11,Xr17,Xr25,Xisrt2,Xcc0)\
	{\
	__asm__ volatile (\
		"ldr	x0,%[__add0]		\n\t	ldr	x5,%[__r1]			\n\t"\
		"ldr	w1,%[__p1]			\n\t	ldr	w6,%[__p8]			\n\t"\
		"ldr	w2,%[__p2]			\n\t	add	x1, x0,x1,lsl #3	\n\t"\
		"ldr	w3,%[__p3]			\n\t	add	x2, x0,x2,lsl #3	\n\t"\
		"ldr	w4,%[__p4]			\n\t	add	x3, x0,x3,lsl #3	\n\t"\
		"ldr	x7,%[__isrt2]		\n\t	ld1r	{v29.2d},[x7]	\n\t"/* isrt2 not needed until pass 2 */\
		/* SSE2_RADIX4_DIT_0TWIDDLE_B(r1 ):	SSE2_RADIX4_DIT_0TWIDDLE_B(r17): */\
		"ldp	q8 ,q9 ,[x0]		\n\t	add	x8 , x0,x6,lsl #3	\n\t"\
		"ldp	q0 ,q1 ,[x1]		\n\t	add	x9 , x1,x6,lsl #3	\n\t"\
		"ldp	q10,q11,[x2]		\n\t	add	x10, x2,x6,lsl #3	\n\t"\
		"ldp	q4 ,q5 ,[x3]		\n\t	add	x11, x3,x6,lsl #3	\n\t"\
		"fsub	v2.2d ,v8.2d ,v0.2d	\n\t	ldp	q20,q21,[x8 ]		\n\t"\
		"fsub	v3.2d ,v9.2d ,v1.2d	\n\t	ldp	q12,q13,[x9 ]		\n\t"\
		"fadd	v0.2d ,v8.2d ,v0.2d	\n\t	ldp	q22,q23,[x10]		\n\t"\
		"fadd	v1.2d ,v9.2d ,v1.2d	\n\t	ldp	q16,q17,[x11]		\n\t"\
		"fsub	v6.2d ,v10.2d,v4.2d	\n\t	fsub	v14.2d,v20.2d,v12.2d	\n\t"\
		"fsub	v7.2d ,v11.2d,v5.2d	\n\t	fsub	v15.2d,v21.2d,v13.2d	\n\t"\
		"fadd	v4.2d ,v10.2d,v4.2d	\n\t	fadd	v12.2d,v20.2d,v12.2d	\n\t"\
		"fadd	v5.2d ,v11.2d,v5.2d	\n\t	fadd	v13.2d,v21.2d,v13.2d	\n\t"\
		"fsub	v8.2d ,v0.2d,v4.2d	\n\t	fsub	v18.2d,v22.2d,v16.2d	\n\t"\
		"fsub	v9.2d ,v1.2d,v5.2d	\n\t	fsub	v19.2d,v23.2d,v17.2d	\n\t"\
		"fadd	v4.2d ,v0.2d,v4.2d	\n\t	fadd	v16.2d,v22.2d,v16.2d	\n\t"\
		"fadd	v5.2d ,v1.2d,v5.2d	\n\t	fadd	v17.2d,v23.2d,v17.2d	\n\t"\
		"fsub	v10.2d,v2.2d,v7.2d	\n\t	fsub	v20.2d,v12.2d,v16.2d	\n\t"\
		"fsub	v11.2d,v3.2d,v6.2d	\n\t	fsub	v21.2d,v13.2d,v17.2d	\n\t"\
		"fadd	v7.2d ,v2.2d,v7.2d	\n\t	fadd	v16.2d,v12.2d,v16.2d	\n\t"\
		"fadd	v6.2d ,v3.2d,v6.2d	\n\t	fadd	v17.2d,v13.2d,v17.2d	\n\t"\
		"stp	q4 ,q5 ,[x5      ]	\n\t	fsub	v22.2d,v14.2d,v19.2d	\n\t"\
		"stp	q7 ,q11,[x5,#0x20]	\n\t	fsub	v23.2d,v15.2d,v18.2d	\n\t"\
		"stp	q8 ,q9 ,[x5,#0x40]	\n\t	fadd	v19.2d,v14.2d,v19.2d	\n\t"\
		"stp	q10,q6 ,[x5,#0x60]	\n\t	fadd	v18.2d,v15.2d,v18.2d	\n\t"\
		"add	x0, x0,x4,lsl #3	\n\t	stp	q16,q17,[x5,#0x100]	\n\t"\
		"add	x1, x1,x4,lsl #3	\n\t	stp	q19,q23,[x5,#0x120]	\n\t"\
		"add	x2, x2,x4,lsl #3	\n\t	stp	q20,q21,[x5,#0x140]	\n\t"\
		"add	x3, x3,x4,lsl #3	\n\t	stp	q22,q18,[x5,#0x160]	\n\t"\
		/* SSE2_RADIX4_DIT_0TWIDDLE_B(r9 ):	SSE2_RADIX4_DIT_0TWIDDLE_B(r25): */\
		"ldp	q8 ,q9 ,[x0]		\n\t	add	x8 , x8 ,x4,lsl #3	\n\t"\
		"ldp	q0 ,q1 ,[x1]		\n\t	add	x9 , x9 ,x4,lsl #3	\n\t"\
		"ldp	q10,q11,[x2]		\n\t	add	x10, x10,x4,lsl #3	\n\t"\
		"ldp	q4 ,q5 ,[x3]		\n\t	add	x11, x11,x4,lsl #3	\n\t"\
		"fsub	v2.2d ,v8.2d ,v0.2d	\n\t	ldp	q20,q21,[x8 ]		\n\t"\
		"fsub	v3.2d ,v9.2d ,v1.2d	\n\t	ldp	q12,q13,[x9 ]		\n\t"\
		"fadd	v0.2d ,v8.2d ,v0.2d	\n\t	ldp	q22,q23,[x10]		\n\t"\
		"fadd	v1.2d ,v9.2d ,v1.2d	\n\t	ldp	q16,q17,[x11]		\n\t"\
		"fsub	v6.2d ,v10.2d,v4.2d	\n\t	fsub	v14.2d,v20.2d,v12.2d	\n\t"\
		"fsub	v7.2d ,v11.2d,v5.2d	\n\t	fsub	v15.2d,v21.2d,v13.2d	\n\t"\
		"fadd	v4.2d ,v10.2d,v4.2d	\n\t	fadd	v12.2d,v20.2d,v12.2d	\n\t"\
		"fadd	v5.2d ,v11.2d,v5.2d	\n\t	fadd	v13.2d,v21.2d,v13.2d	\n\t"\
		"fsub	v8.2d ,v0.2d,v4.2d	\n\t	fsub	v18.2d,v22.2d,v16.2d	\n\t"\
		"fsub	v9.2d ,v1.2d,v5.2d	\n\t	fsub	v19.2d,v23.2d,v17.2d	\n\t"\
		"fadd	v4.2d ,v0.2d,v4.2d	\n\t	fadd	v16.2d,v22.2d,v16.2d	\n\t"\
		"fadd	v5.2d ,v1.2d,v5.2d	\n\t	fadd	v17.2d,v23.2d,v17.2d	\n\t"\
		"fsub	v10.2d,v2.2d,v7.2d	\n\t	fsub	v20.2d,v12.2d,v16.2d	\n\t"\
		"fsub	v11.2d,v3.2d,v6.2d	\n\t	fsub	v21.2d,v13.2d,v17.2d	\n\t"\
		"fadd	v7.2d ,v2.2d,v7.2d	\n\t	fadd	v16.2d,v12.2d,v16.2d	\n\t"\
		"fadd	v6.2d ,v3.2d,v6.2d	\n\t	fadd	v17.2d,v13.2d,v17.2d	\n\t"\
		"stp	q4 ,q5 ,[x5,#0x80]	\n\t	fsub	v22.2d,v14.2d,v19.2d	\n\t"\
		"stp	q7 ,q11,[x5,#0xa0]	\n\t	fsub	v23.2d,v15.2d,v18.2d	\n\t"\
		"stp	q8 ,q9 ,[x5,#0xc0]	\n\t	fadd	v19.2d,v14.2d,v19.2d	\n\t"\
		"stp	q10,q6 ,[x5,#0xe0]	\n\t	fadd	v18.2d,v15.2d,v18.2d	\n\t"\
		/* v20,21 used in next rcol 4-DFT, instead of storing, do reg-copy 20->16,21->30,
		But must carefully arrange order of STP, LDP and MOV here to ensure each lcol datum
		is safely stored before register overwritten by rcol LDP or reg-copy MOV: */\
		"stp	q16,q17,[x5,#0x180]	\n\t	mov	v30.16b,v21.16b			\n\t"/* 16,30 have 6,7 */\
		"mov	v16.16b,v20.16b		\n\t	ldp	q20,q21,[x5,#0x040]		\n\t"/* 20,21 have 0,1 */\
		"stp	q19,q23,[x5,#0x1a0]	\n\t	ldp	q12,q13,[x5,#0x0c0]		\n\t"/* 12,13 have 2,3 */\
		"stp	q22,q18,[x5,#0x1e0]	\n\t	ldp	q22,q23,[x5,#0x140]		\n\t"/* 22,23 have 4,5 */\
	/*** Now do four pass-2 4-DFTs, inputs from local store, outputs back to same: ***/\
		/* Block 0:							Block 2 (loads for same above): */\
		"ldp	q8 ,q9 ,[x5       ]	\n\t	fadd	v18.2d,v23.2d,v22.2d	\n\t"\
		"ldp	q0 ,q1 ,[x5,#0x080]	\n\t	fsub	v19.2d,v23.2d,v22.2d	\n\t"\
		"ldp	q10,q11,[x5,#0x100]	\n\t	fadd	v17.2d,v16.2d,v30.2d	\n\t"\
		"ldp	q4 ,q5 ,[x5,#0x180]	\n\t	fsub	v16.2d,v16.2d,v30.2d	\n\t"\
		"fsub	v2.2d ,v8.2d ,v0.2d	\n\t	fmul	v18.2d,v29.2d,v18.2d	\n\t"\
		"fsub	v3.2d ,v9.2d ,v1.2d	\n\t	fmul	v19.2d,v29.2d,v19.2d	\n\t"\
		"fadd	v0.2d ,v8.2d ,v0.2d	\n\t	fmul	v16.2d,v29.2d,v16.2d	\n\t"\
		"fadd	v1.2d ,v9.2d ,v1.2d	\n\t	fmul	v17.2d,v29.2d,v17.2d	\n\t"\
		"fsub	v6.2d ,v10.2d,v4.2d	\n\t	fsub	v22.2d,v18.2d,v16.2d	\n\t"\
		"fsub	v7.2d ,v11.2d,v5.2d	\n\t	fadd	v23.2d,v18.2d,v16.2d	\n\t"\
		"fadd	v4.2d ,v10.2d,v4.2d	\n\t	fsub	v16.2d,v19.2d,v17.2d	\n\t"\
		"fadd	v5.2d ,v11.2d,v5.2d	\n\t	fadd	v17.2d,v19.2d,v17.2d	\n\t"\
		"fsub	v8.2d ,v0.2d,v4.2d	\n\t	fsub	v18.2d,v20.2d,v13.2d	\n\t"\
		"fsub	v9.2d ,v1.2d,v5.2d	\n\t	fadd	v19.2d,v20.2d,v13.2d	\n\t"\
		"fadd	v4.2d ,v0.2d,v4.2d	\n\t	fsub	v20.2d,v21.2d,v12.2d	\n\t"\
		"fadd	v5.2d ,v1.2d,v5.2d	\n\t	fadd	v21.2d,v21.2d,v12.2d	\n\t"\
		"fsub	v10.2d,v2.2d,v7.2d	\n\t	fsub	v12.2d,v18.2d,v17.2d	\n\t"\
		"fsub	v11.2d,v3.2d,v6.2d	\n\t	fadd	v13.2d,v18.2d,v17.2d	\n\t"\
		"fadd	v7.2d ,v2.2d,v7.2d	\n\t	fsub	v18.2d,v19.2d,v22.2d	\n\t"\
		"fadd	v6.2d ,v3.2d,v6.2d	\n\t	fadd	v19.2d,v19.2d,v22.2d	\n\t"\
		"stp	q4 ,q5 ,[x5       ]	\n\t	fsub	v14.2d,v20.2d,v16.2d	\n\t"\
		"stp	q7 ,q11,[x5,#0x080]	\n\t	fadd	v15.2d,v20.2d,v16.2d	\n\t"\
		"stp	q8 ,q9 ,[x5,#0x100]	\n\t	fsub	v20.2d,v21.2d,v23.2d	\n\t"\
		"stp	q10,q6 ,[x5,#0x180]	\n\t	fadd	v21.2d,v21.2d,v23.2d	\n\t"\
		"ldp	q30,q31,[x7,#0x10]	\n\t"/* cc0,ss0; also do Block 1 loads in lcol: */\
		"ldp	q10,q11,[x5,#0x020]	\n\t	stp	q19,q15,[x5,#0x040]	\n\t"\
		"ldp	q2 ,q3 ,[x5,#0x0a0]	\n\t	stp	q13,q20,[x5,#0x0c0]	\n\t"\
		"ldp	q6 ,q7 ,[x5,#0x120]	\n\t	stp	q18,q14,[x5,#0x140]	\n\t"\
		"ldp	q8 ,q9 ,[x5,#0x1a0]	\n\t	stp	q12,q21,[x5,#0x1c0]	\n\t"\
		/* Block 1:							Block 3: */\
		"fmul	v0.2d,v8.2d,v31.2d	\n\t	ldp	q22,q23,[x5,#0x060]	\n\t"\
		"fmul	v1.2d,v9.2d,v31.2d	\n\t	ldp	q14,q15,[x5,#0x0e0]	\n\t"\
		"fmla	v0.2d,v9.2d,v30.2d	\n\t	ldp	q18,q19,[x5,#0x160]	\n\t"\
		"fmls	v1.2d,v8.2d,v30.2d	\n\t	ldp	q20,q21,[x5,#0x1e0]	\n\t"\
		"fmul	v8.2d,v6.2d,v30.2d	\n\t	fmul	v12.2d,v20.2d,v30.2d	\n\t"\
		"fmul	v9.2d,v7.2d,v30.2d	\n\t	fmul	v13.2d,v21.2d,v30.2d	\n\t"\
		"fmla	v8.2d,v7.2d,v31.2d	\n\t	fmla	v12.2d,v21.2d,v31.2d	\n\t"\
		"fmls	v9.2d,v6.2d,v31.2d	\n\t	fmls	v13.2d,v20.2d,v31.2d	\n\t"\
		"fadd	v4.2d,v8.2d,v0.2d	\n\t	fmul	v20.2d,v18.2d,v31.2d	\n\t"\
		"fadd	v5.2d,v9.2d,v1.2d	\n\t	fmul	v21.2d,v19.2d,v31.2d	\n\t"\
		"fsub	v6.2d,v8.2d,v0.2d	\n\t	fmla	v20.2d,v19.2d,v30.2d	\n\t"\
		"fsub	v7.2d,v9.2d,v1.2d	\n\t	fmls	v21.2d,v18.2d,v30.2d	\n\t"\
		"fadd	v8.2d,v2.2d,v3.2d	\n\t	fadd	v16.2d,v20.2d,v12.2d	\n\t"\
		"fsub	v9.2d,v3.2d,v2.2d	\n\t	fadd	v17.2d,v21.2d,v13.2d	\n\t"\
		"fmul	v8.2d,v8.2d,v29.2d	\n\t	fsub	v18.2d,v20.2d,v12.2d	\n\t"\
		"fmul	v9.2d,v9.2d,v29.2d	\n\t	fsub	v19.2d,v21.2d,v13.2d	\n\t"\
		"fsub	v0.2d,v10.2d,v8.2d	\n\t	fsub	v20.2d,v14.2d,v15.2d	\n\t"\
		"fsub	v1.2d,v11.2d,v9.2d	\n\t	fadd	v21.2d,v15.2d,v14.2d	\n\t"\
		"fadd	v2.2d,v10.2d,v8.2d	\n\t	fmul	v20.2d,v20.2d,v29.2d	\n\t"\
		"fadd	v3.2d,v11.2d,v9.2d	\n\t	fmul	v21.2d,v21.2d,v29.2d	\n\t"\
		"fadd	v8.2d,v2.2d,v4.2d	\n\t	fsub	v12.2d,v22.2d,v20.2d	\n\t"\
		"fadd	v9.2d,v3.2d,v5.2d	\n\t	fsub	v13.2d,v23.2d,v21.2d	\n\t"\
		"fsub	v2.2d,v2.2d,v4.2d	\n\t	fadd	v14.2d,v22.2d,v20.2d	\n\t"\
		"fsub	v3.2d,v3.2d,v5.2d	\n\t	fadd	v15.2d,v23.2d,v21.2d	\n\t"\
		"fadd	v4.2d,v0.2d,v7.2d	\n\t	fadd	v20.2d,v12.2d,v18.2d	\n\t"\
		"fadd	v5.2d,v1.2d,v6.2d	\n\t	fadd	v21.2d,v13.2d,v19.2d	\n\t"\
		"fsub	v0.2d,v0.2d,v7.2d	\n\t	fsub	v12.2d,v12.2d,v18.2d	\n\t"\
		"fsub	v1.2d,v1.2d,v6.2d	\n\t	fsub	v13.2d,v13.2d,v19.2d	\n\t"\
		"stp	q8 ,q9 ,[x5,#0x020]	\n\t	fadd	v18.2d,v14.2d,v17.2d	\n\t"\
		"stp	q2 ,q3 ,[x5,#0x120]	\n\t	fadd	v19.2d,v15.2d,v16.2d	\n\t"\
		"stp	q4 ,q1 ,[x5,#0x0a0]	\n\t	fsub	v14.2d,v14.2d,v17.2d	\n\t"\
		"stp	q0 ,q5 ,[x5,#0x1a0]	\n\t	fsub	v15.2d,v15.2d,v16.2d	\n\t"\
		"									stp	q20,q21,[x5,#0x060]	\n\t"\
		"									stp	q12,q13,[x5,#0x160]	\n\t"\
		"									stp	q18,q15,[x5,#0x0e0]	\n\t"\
		"									stp	q14,q19,[x5,#0x1e0]	\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__p1] "m" (Xp1)\
		 ,[__p2] "m" (Xp2)\
		 ,[__p3] "m" (Xp3)\
		 ,[__p4] "m" (Xp4)\
		 ,[__p8] "m" (Xp8)\
		 ,[__r1] "m" (Xr1)\
		 ,[__r3] "m" (Xr3)\
		 ,[__r5] "m" (Xr5)\
		 ,[__r7] "m" (Xr7)\
		 ,[__r9] "m" (Xr9)\
		 ,[__r17] "m" (Xr17)\
		 ,[__r25] "m" (Xr25)\
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__cc0] "m" (Xcc0)\
		: "cc","memory","x0","x1","x2","x3","x4","x5","x6","x7","v0","v1","v2","v3","v4","v5","v6","v7","v8","v9","v10","v11",\
			"x8","x9","x10","x11","v12","v13","v14","v15","v16","v17","v18","v19","v20","v21","v22","v23", "v29","v30","v31"	/* Clobbered registers */\
	);\
	}

	#define SSE2_RADIX16_DIF_NOTWIDDLE(Xadd0,Xp1,Xp2,Xp3,Xp4,Xp8,Xr1,Xr3,Xr5,Xr7,Xr9,Xr11,Xr13,Xr15,Xr17,Xr19,Xr21,Xr23,Xr25,Xr27,Xr29,Xr31,Xisrt2,Xcc0)\
	{\
	__asm__ volatile (\
		"ldr	x0,%[__add0]		\n\t	ldr	x5,%[__r1]			\n\t"\
		"ldr	w1,%[__p1]			\n\t	ldr	w6,%[__p8]			\n\t"\
		"ldr	w2,%[__p2]			\n\t	add	x1, x0,x1,lsl #3	\n\t"\
		"ldr	w3,%[__p3]			\n\t	add	x2, x0,x2,lsl #3	\n\t"\
		"ldr	w4,%[__p4]			\n\t	add	x3, x0,x3,lsl #3	\n\t"\
		/* isrt2 and add0-offset addresses not needed until pass 2, but no gain
		from delaying their computation, so just re-use the DIT-macro code: */\
		"ldr	x7,%[__isrt2]		\n\t	ld1r	{v29.2d},[x7]	\n\t"\
		/* SSE2_RADIX4_DIF_IN_PLACE(r1,17,9,25)		(r5,21,13,29): */\
		"ldp	q8 ,q9 ,[x5       ]	\n\t	add	x8 , x0,x4,lsl #3	\n\t"\
		"ldp	q0 ,q1 ,[x5,#0x100]	\n\t	add	x9 , x1,x4,lsl #3	\n\t"\
		"ldp	q10,q11,[x5,#0x080]	\n\t	add	x10, x2,x4,lsl #3	\n\t"\
		"ldp	q4 ,q5 ,[x5,#0x180]	\n\t	add	x11, x3,x4,lsl #3	\n\t"\
		"fsub	v2.2d ,v8.2d ,v0.2d	\n\t	ldp	q20,q21,[x5,#0x040]	\n\t"\
		"fsub	v3.2d ,v9.2d ,v1.2d	\n\t	ldp	q12,q13,[x5,#0x140]	\n\t"\
		"fadd	v0.2d ,v8.2d ,v0.2d	\n\t	ldp	q22,q23,[x5,#0x0c0]	\n\t"\
		"fadd	v1.2d ,v9.2d ,v1.2d	\n\t	ldp	q16,q17,[x5,#0x1c0]	\n\t"\
		"fsub	v6.2d ,v10.2d,v4.2d	\n\t	fsub	v14.2d,v20.2d,v12.2d	\n\t"\
		"fsub	v7.2d ,v11.2d,v5.2d	\n\t	fsub	v15.2d,v21.2d,v13.2d	\n\t"\
		"fadd	v4.2d ,v10.2d,v4.2d	\n\t	fadd	v12.2d,v20.2d,v12.2d	\n\t"\
		"fadd	v5.2d ,v11.2d,v5.2d	\n\t	fadd	v13.2d,v21.2d,v13.2d	\n\t"\
		"fsub	v8.2d ,v0.2d,v4.2d	\n\t	fsub	v18.2d,v22.2d,v16.2d	\n\t"\
		"fsub	v9.2d ,v1.2d,v5.2d	\n\t	fsub	v19.2d,v23.2d,v17.2d	\n\t"\
		"fadd	v4.2d ,v0.2d,v4.2d	\n\t	fadd	v16.2d,v22.2d,v16.2d	\n\t"\
		"fadd	v5.2d ,v1.2d,v5.2d	\n\t	fadd	v17.2d,v23.2d,v17.2d	\n\t"\
		"fsub	v10.2d,v2.2d,v7.2d	\n\t	fsub	v20.2d,v12.2d,v16.2d	\n\t"\
		"fsub	v11.2d,v3.2d,v6.2d	\n\t	fsub	v21.2d,v13.2d,v17.2d	\n\t"\
		"fadd	v7.2d ,v2.2d,v7.2d	\n\t	fadd	v16.2d,v12.2d,v16.2d	\n\t"\
		"fadd	v6.2d ,v3.2d,v6.2d	\n\t	fadd	v17.2d,v13.2d,v17.2d	\n\t"\
		"stp	q4 ,q5 ,[x5       ]	\n\t	fsub	v22.2d,v14.2d,v19.2d	\n\t"\
		"stp	q10,q6 ,[x5,#0x080]	\n\t	fsub	v23.2d,v15.2d,v18.2d	\n\t"\
		"stp	q8 ,q9 ,[x5,#0x100]	\n\t	fadd	v19.2d,v14.2d,v19.2d	\n\t"\
		"stp	q7 ,q11,[x5,#0x180]	\n\t	fadd	v18.2d,v15.2d,v18.2d	\n\t"\
		"ldp	q8 ,q9 ,[x5,#0x020]	\n\t	stp	q16,q17,[x5,#0x040]	\n\t"\
		"ldp	q0 ,q1 ,[x5,#0x120]	\n\t	stp	q22,q18,[x5,#0x0c0]	\n\t"\
		"ldp	q10,q11,[x5,#0x0a0]	\n\t	stp	q20,q21,[x5,#0x140]	\n\t"\
		"ldp	q4 ,q5 ,[x5,#0x1a0]	\n\t	stp	q19,q23,[x5,#0x1c0]	\n\t"\
		/* SSE2_RADIX4_DIF_IN_PLACE(r3,19,11,27)		(r7,23,15,31): */\
		"fsub	v2.2d ,v8.2d ,v0.2d	\n\t	ldp	q20,q21,[x5,#0x060]	\n\t"\
		"fsub	v3.2d ,v9.2d ,v1.2d	\n\t	ldp	q12,q13,[x5,#0x160]	\n\t"\
		"fadd	v0.2d ,v8.2d ,v0.2d	\n\t	ldp	q22,q23,[x5,#0x0e0]	\n\t"\
		"fadd	v1.2d ,v9.2d ,v1.2d	\n\t	ldp	q16,q17,[x5,#0x1e0]	\n\t"\
		"fsub	v6.2d ,v10.2d,v4.2d	\n\t	fsub	v14.2d,v20.2d,v12.2d	\n\t"\
		"fsub	v7.2d ,v11.2d,v5.2d	\n\t	fsub	v15.2d,v21.2d,v13.2d	\n\t"\
		"fadd	v4.2d ,v10.2d,v4.2d	\n\t	fadd	v12.2d,v20.2d,v12.2d	\n\t"\
		"fadd	v5.2d ,v11.2d,v5.2d	\n\t	fadd	v13.2d,v21.2d,v13.2d	\n\t"\
		"fsub	v8.2d ,v0.2d,v4.2d	\n\t	fsub	v18.2d,v22.2d,v16.2d	\n\t"\
		"fsub	v9.2d ,v1.2d,v5.2d	\n\t	fsub	v19.2d,v23.2d,v17.2d	\n\t"\
		"fadd	v4.2d ,v0.2d,v4.2d	\n\t	fadd	v16.2d,v22.2d,v16.2d	\n\t"\
		"fadd	v5.2d ,v1.2d,v5.2d	\n\t	fadd	v17.2d,v23.2d,v17.2d	\n\t"\
		"fsub	v10.2d,v2.2d,v7.2d	\n\t	fsub	v20.2d,v12.2d,v16.2d	\n\t"\
		"fsub	v11.2d,v3.2d,v6.2d	\n\t	fsub	v21.2d,v13.2d,v17.2d	\n\t"\
		"fadd	v7.2d ,v2.2d,v7.2d	\n\t	fadd	v16.2d,v12.2d,v16.2d	\n\t"\
		"fadd	v6.2d ,v3.2d,v6.2d	\n\t	fadd	v17.2d,v13.2d,v17.2d	\n\t"\
		"stp	q4 ,q5 ,[x5,#0x020]	\n\t	fsub	v22.2d,v14.2d,v19.2d	\n\t"\
		"stp	q10,q6 ,[x5,#0x0a0]	\n\t	fsub	v23.2d,v15.2d,v18.2d	\n\t"\
		"stp	q8 ,q9 ,[x5,#0x120]	\n\t	fadd	v19.2d,v14.2d,v19.2d	\n\t"\
		"stp	q7 ,q11,[x5,#0x1a0]	\n\t	fadd	v18.2d,v15.2d,v18.2d	\n\t"\
		/* v20,21 used in next rcol 4-DFT, instead of storing, do reg-copy 20->14,21->15,
		But must carefully arrange order of STP, LDP and MOV here to ensure each lcol datum
		is safely stored before register overwritten by rcol LDP or reg-copy MOV: */\
		"stp	q16,q17,[x5,#0x060]	\n\t	ldp	q24,q25,[x5,#0x140]	\n\t"/* 2,3 */\
		"stp	q22,q18,[x5,#0x0e0]	\n\t	ldp	q12,q13,[x5,#0x120]	\n\t"/* 4,5 */\
		"mov	v14.16b,v20.16b		\n\t	mov	v15.16b,v21.16b		\n\t"/* 6,7 */\
		"stp	q19,q23,[x5,#0x1e0]	\n\t	ldp	q22,q23,[x5,#0x100]	\n\t"/* 0,1 - load this pair last because need to wait for lcol writes of q22,23 to issue */\
	/*** Now do four pass-2 4-DFTs, inputs from local store, outputs back to main-array: ***/\
		/* Block 0:							Block 2 (loads for same above):		rcol ops: */\
		"ldp	q8 ,q9 ,[x5       ]	\n\t	fsub	v16.2d,v12.2d,v13.2d	\n\t"/* 4 = 4-5 */\
		"ldp	q2 ,q3 ,[x5,#0x040]	\n\t	fadd	v17.2d,v12.2d,v13.2d	\n\t"/* 5 = 4+5 */\
		"ldp	q10,q11,[x5,#0x020]	\n\t	fadd	v18.2d,v15.2d,v14.2d	\n\t"/* 6 = 7+6 */\
		"ldp	q6 ,q7 ,[x5,#0x060]	\n\t	fsub	v19.2d,v15.2d,v14.2d	\n\t"/* 7 = 7-6 */\
		"fsub	v0.2d ,v8.2d ,v2.2d	\n\t	fmul	v20.2d,v29.2d,v16.2d	\n\t"/* 4*isrt2 */\
		"fsub	v1.2d ,v9.2d ,v3.2d	\n\t	fmul	v21.2d,v29.2d,v17.2d	\n\t"/* 5*isrt2 */\
		"fadd	v2.2d ,v8.2d ,v2.2d	\n\t	fmul	v18.2d,v29.2d,v18.2d	\n\t"/* 6*isrt2 */\
		"fadd	v3.2d ,v9.2d ,v3.2d	\n\t	fmul	v19.2d,v29.2d,v19.2d	\n\t"/* 7*isrt2 */\
		"fsub	v4.2d ,v10.2d,v6.2d	\n\t	fsub	v12.2d,v22.2d,v25.2d	\n\t"/* 0 = 0-3 */\
		"fsub	v5.2d ,v11.2d,v7.2d	\n\t	fsub	v13.2d,v23.2d,v24.2d	\n\t"/* 1 = 1-2 */\
		"fadd	v6.2d ,v10.2d,v6.2d	\n\t	fadd	v14.2d,v23.2d,v24.2d	\n\t"/* 2 = 1+2 */\
		"fadd	v7.2d ,v11.2d,v7.2d	\n\t	fadd	v15.2d,v22.2d,v25.2d	\n\t"/* 3 = 0+3 */\
		"fsub	v8.2d ,v2.2d,v6.2d	\n\t	fsub	v16.2d,v20.2d,v18.2d	\n\t"/* 4 = 4-6 */\
		"fsub	v9.2d ,v3.2d,v7.2d	\n\t	fsub	v17.2d,v21.2d,v19.2d	\n\t"/* 5 = 5-7 */\
		"fadd	v6.2d ,v2.2d,v6.2d	\n\t	fadd	v18.2d,v20.2d,v18.2d	\n\t"/* 6 = 4+6 */\
		"fadd	v7.2d ,v3.2d,v7.2d	\n\t	fadd	v19.2d,v21.2d,v19.2d	\n\t"/* 7 = 5+7 */\
		"fsub	v10.2d,v0.2d,v5.2d	\n\t	fsub	v20.2d,v12.2d,v16.2d	\n\t"/* 0 = 0-4 */\
		"fsub	v11.2d,v1.2d,v4.2d	\n\t	fadd	v16.2d,v12.2d,v16.2d	\n\t"/* 4 = 0+4 */\
		"fadd	v5.2d ,v0.2d,v5.2d	\n\t	fsub	v21.2d,v14.2d,v17.2d	\n\t"/* 2 = 2-5 */\
		"fadd	v4.2d ,v1.2d,v4.2d	\n\t	fadd	v17.2d,v14.2d,v17.2d	\n\t"/* 5 = 2+5 */\
		"stp	q6 ,q7 ,[x0]		\n\t	fsub	v22.2d,v13.2d,v18.2d	\n\t"/* 1 = 1-6 */\
		"stp	q8 ,q9 ,[x1]		\n\t	fadd	v18.2d,v13.2d,v18.2d	\n\t"/* 6 = 1+6 */\
		"stp	q10,q4 ,[x2]		\n\t	fsub	v23.2d,v15.2d,v19.2d	\n\t"/* 3 = 3-7 */\
		"stp	q5 ,q11,[x3]		\n\t	fadd	v19.2d,v15.2d,v19.2d	\n\t"/* 7 = 3+7 */\
		"ldp	q30,q31,[x7,#0x10]	\n\t"/* cc0,ss0 */\
		"add	x0, x0,x6,lsl #3	\n\t	stp	q16,q17,[x8 ]	\n\t"\
		"add	x1, x1,x6,lsl #3	\n\t	stp	q20,q21,[x9 ]	\n\t"\
		"add	x2, x2,x6,lsl #3	\n\t	stp	q23,q18,[x10]	\n\t"\
		"add	x3, x3,x6,lsl #3	\n\t	stp	q19,q22,[x11]	\n\t"\
		/* Block 1:							Block 3: */\
		"ldp	q0 ,q1 ,[x5,#0x080]	\n\t	add	x8 , x8 ,x6,lsl #3	\n\t"\
		"ldp	q2 ,q3 ,[x5,#0x0c0]	\n\t	add	x9 , x9 ,x6,lsl #3	\n\t"\
		"ldp	q4 ,q5 ,[x5,#0x0a0]	\n\t	add	x10, x10,x6,lsl #3	\n\t"\
		"ldp	q10,q11,[x5,#0x0e0]	\n\t	add	x11, x11,x6,lsl #3	\n\t"\
		"fmul	v8.2d,v4.2d,v30.2d	\n\t	ldp	q12,q13,[x5,#0x180]		\n\t"\
		"fmul	v9.2d,v5.2d,v30.2d	\n\t	ldp	q14,q15,[x5,#0x1c0]		\n\t"/* lcol ops: */\
		"fmls	v8.2d,v5.2d,v31.2d	\n\t	ldp	q16,q17,[x5,#0x1a0]		\n\t"/* 4 = 4c-5s */\
		"fmla	v9.2d,v4.2d,v31.2d	\n\t	ldp	q18,q19,[x5,#0x1e0]		\n\t"/* 5 = 4s+5c */\
		"fmul	v6.2d,v10.2d,v31.2d	\n\t	fmul	v20.2d,v16.2d,v31.2d \n\t"\
		"fmul	v7.2d,v11.2d,v31.2d	\n\t	fmul	v21.2d,v17.2d,v31.2d \n\t"/*			rcol ops: */\
		"fmls	v6.2d,v11.2d,v30.2d	\n\t	fmls	v20.2d,v17.2d,v30.2d \n\t"/* 6 = 6s-7c	4 = 4s-5c */\
		"fmla	v7.2d,v10.2d,v30.2d	\n\t	fmla	v21.2d,v16.2d,v30.2d \n\t"/* 7 = 7s+6c	5 = 4c+5s */\
		"fsub	v4.2d,v8.2d,v6.2d	\n\t	fmul	v22.2d,v18.2d,v30.2d \n\t"/* 4 = 4-6 */\
		"fsub	v5.2d,v9.2d,v7.2d	\n\t	fmul	v23.2d,v19.2d,v30.2d \n\t"/* 5 = 5-7 */\
		"fadd	v6.2d,v8.2d,v6.2d	\n\t	fmls	v22.2d,v19.2d,v31.2d \n\t"/* 6 = 4+6	6 = 6c-7s */\
		"fadd	v7.2d,v9.2d,v7.2d	\n\t	fmla	v23.2d,v18.2d,v31.2d \n\t"/* 7 = 5+7	7 = 6s+7c */\
		"fsub	v8.2d,v2.2d,v3.2d	\n\t	fsub	v16.2d,v20.2d,v22.2d \n\t"/* 2 = 2-3	4 = 4-6 */\
		"fadd	v9.2d,v2.2d,v3.2d	\n\t	fsub	v17.2d,v21.2d,v23.2d \n\t"/* 3 = 2+3	5 = 5-7 */\
		"fmul	v8.2d,v8.2d,v29.2d	\n\t	fadd	v18.2d,v20.2d,v22.2d \n\t"/* 2*isrt2	6 = 4+6 */\
		"fmul	v9.2d,v9.2d,v29.2d	\n\t	fadd	v19.2d,v21.2d,v23.2d \n\t"/* 3*isrt2	7 = 5+7 */\
		"fadd	v2.2d,v0.2d,v8.2d	\n\t	fadd	v20.2d,v15.2d,v14.2d \n\t"/* 2 = 0+2	2 = 3+2 */\
		"fadd	v3.2d,v1.2d,v9.2d	\n\t	fsub	v21.2d,v15.2d,v14.2d \n\t"/* 3 = 1+3	3 = 3-2 */\
		"fsub	v0.2d,v0.2d,v8.2d	\n\t	fmul	v20.2d,v20.2d,v29.2d \n\t"/* 0 = 0-2	2*isrt2 */\
		"fsub	v1.2d,v1.2d,v9.2d	\n\t	fmul	v21.2d,v21.2d,v29.2d \n\t"/* 1 = 1-3	3*isrt2 */\
		"fadd	v8.2d,v2.2d,v6.2d	\n\t	fadd	v14.2d,v12.2d,v20.2d \n\t"/* 6 = 2+6	2 = 0+2 */\
		"fadd	v9.2d,v3.2d,v7.2d	\n\t	fadd	v15.2d,v13.2d,v21.2d \n\t"/* 7 = 3+7	3 = 1+3 */\
		"fsub	v2.2d,v2.2d,v6.2d	\n\t	fsub	v12.2d,v12.2d,v20.2d \n\t"/* 2 = 2-6	0 = 0-2 */\
		"fsub	v3.2d,v3.2d,v7.2d	\n\t	fsub	v13.2d,v13.2d,v21.2d \n\t"/* 3 = 3-7	1 = 1-3 */\
		"fadd	v10.2d,v0.2d,v5.2d	\n\t	fadd	v20.2d,v12.2d,v16.2d \n\t"/* 5 = 0+5	4 = 0+4 */\
		"fadd	v11.2d,v1.2d,v4.2d	\n\t	fadd	v21.2d,v13.2d,v17.2d \n\t"/* 4 = 1+4	5 = 1+5 */\
		"fsub	v0.2d,v0.2d,v5.2d	\n\t	fsub	v12.2d,v12.2d,v16.2d \n\t"/* 0 = 0-5	0 = 0-4 */\
		"fsub	v1.2d,v1.2d,v4.2d	\n\t	fsub	v13.2d,v13.2d,v17.2d \n\t"/* 1 = 1-4	1 = 1-5 */\
		"stp	q8 ,q9 ,[x0]		\n\t	fadd	v22.2d,v14.2d,v19.2d \n\t"/* 6,7		7 = 2+7 */\
		"stp	q2 ,q3 ,[x1]		\n\t	fadd	v23.2d,v15.2d,v18.2d \n\t"/* 2,3		6 = 3+6 */\
		"stp	q0 ,q11,[x2]		\n\t	fsub	v14.2d,v14.2d,v19.2d \n\t"/* 0,4		2 = 2-7 */\
		"stp	q10,q1 ,[x3]		\n\t	fsub	v15.2d,v15.2d,v18.2d \n\t"/* 5,1		3 = 3-6 */\
		"									stp	q20,q21,[x8 ]	\n\t"\
		"									stp	q12,q13,[x9 ]	\n\t"\
		"									stp	q14,q23,[x10]	\n\t"\
		"									stp	q22,q15,[x11]	\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__p1] "m" (Xp1)\
		 ,[__p2] "m" (Xp2)\
		 ,[__p3] "m" (Xp3)\
		 ,[__p4] "m" (Xp4)\
		 ,[__p8] "m" (Xp8)\
		 ,[__r1] "m" (Xr1)\
		 ,[__r3] "m" (Xr3)\
		 ,[__r5] "m" (Xr5)\
		 ,[__r7] "m" (Xr7)\
		 ,[__r9] "m" (Xr9)\
		 ,[__r11] "m" (Xr11)\
		 ,[__r13] "m" (Xr13)\
		 ,[__r15] "m" (Xr15)\
		 ,[__r17] "m" (Xr17)\
		 ,[__r19] "m" (Xr19)\
		 ,[__r21] "m" (Xr21)\
		 ,[__r23] "m" (Xr23)\
		 ,[__r25] "m" (Xr25)\
		 ,[__r27] "m" (Xr27)\
		 ,[__r29] "m" (Xr29)\
		 ,[__r31] "m" (Xr31)\
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__cc0] "m" (Xcc0)\
		: "cc","memory","x0","x1","x2","x3","x4","x5","x6","x7","v0","v1","v2","v3","v4","v5","v6","v7","v8","v9","v10","v11",\
		"x8","x9","x10","x11","v12","v13","v14","v15","v16","v17","v18","v19","v20","v21","v22","v23","v24","v25", "v29","v30","v31"	/* Clobbered registers */\
	);\
	}

#elif defined(USE_AVX512)

	#define SSE2_RADIX16_DIT_NOTWIDDLE(Xadd0,Xp1,Xp2,Xp3,Xp4,Xp8,Xr1,Xr3,Xr5,Xr7,Xr9,Xr11,Xr17,Xr25,Xisrt2,Xcc0)\
	{\
	__asm__ volatile (\
		"movq	%[__add0],%%rax		\n\t"\
		"movslq	%[__p1],%%rbx		\n\t"\
		"movslq	%[__p2],%%rcx		\n\t"\
		"movslq	%[__p3],%%rdx		\n\t"\
		"movslq	%[__p4],%%rdi		\n\t"\
		"shlq	$3,%%rbx			\n\t"\
		"shlq	$3,%%rcx			\n\t"\
		"shlq	$3,%%rdx			\n\t"\
		"shlq	$3,%%rdi			\n\t"\
		"addq	%%rax,%%rbx			\n\t"\
		"addq	%%rax,%%rcx			\n\t"\
		"addq	%%rax,%%rdx			\n\t"\
		/* SSE2_RADIX4_DIT_0TWIDDLE_B(r1 ): */\
		"vmovaps	     (%%rax),%%zmm2	\n\t"\
		"vmovaps	     (%%rcx),%%zmm6	\n\t"\
		"vmovaps	0x040(%%rax),%%zmm3	\n\t"\
		"vmovaps	0x040(%%rcx),%%zmm7	\n\t"\
		"vmovaps	     (%%rbx),%%zmm0	\n\t"\
		"vmovaps	     (%%rdx),%%zmm4	\n\t"\
		"vmovaps	0x040(%%rbx),%%zmm1	\n\t"\
		"vmovaps	0x040(%%rdx),%%zmm5	\n\t"\
		"vsubpd	%%zmm0,%%zmm2,%%zmm2		\n\t"\
		"vsubpd	%%zmm4,%%zmm6,%%zmm6		\n\t"\
		"vsubpd	%%zmm1,%%zmm3,%%zmm3		\n\t"\
		"vsubpd	%%zmm5,%%zmm7,%%zmm7		\n\t"\
		"vaddpd	%%zmm0,%%zmm0,%%zmm0		\n\t"\
		"vaddpd	%%zmm4,%%zmm4,%%zmm4		\n\t"\
		"vaddpd	%%zmm1,%%zmm1,%%zmm1		\n\t"\
		"vaddpd	%%zmm5,%%zmm5,%%zmm5		\n\t"\
		"vaddpd	%%zmm2,%%zmm0,%%zmm0		\n\t"\
		"vaddpd	%%zmm6,%%zmm4,%%zmm4		\n\t"\
		"vaddpd	%%zmm3,%%zmm1,%%zmm1		\n\t"\
		"vaddpd	%%zmm7,%%zmm5,%%zmm5		\n\t"\
		"vsubpd	%%zmm4,%%zmm0,%%zmm0		\n\t"\
		"vsubpd	%%zmm7,%%zmm2,%%zmm2		\n\t"\
		"vsubpd	%%zmm5,%%zmm1,%%zmm1		\n\t"\
		"vsubpd	%%zmm6,%%zmm3,%%zmm3		\n\t"\
		"movq	%[__r1], %%rsi		\n\t"\
		"vmovaps	%%zmm0,0x100(%%rsi)	\n\t"\
		"vmovaps	%%zmm2,0x180(%%rsi)	\n\t"\
		"vmovaps	%%zmm1,0x140(%%rsi)	\n\t"\
		"vmovaps	%%zmm3,0x0c0(%%rsi)	\n\t"\
		"vaddpd	%%zmm4,%%zmm4,%%zmm4		\n\t"\
		"vaddpd	%%zmm7,%%zmm7,%%zmm7		\n\t"\
		"vaddpd	%%zmm5,%%zmm5,%%zmm5		\n\t"\
		"vaddpd	%%zmm6,%%zmm6,%%zmm6		\n\t"\
		"vaddpd	%%zmm0,%%zmm4,%%zmm4		\n\t"\
		"vaddpd	%%zmm2,%%zmm7,%%zmm7		\n\t"\
		"vaddpd	%%zmm1,%%zmm5,%%zmm5		\n\t"\
		"vaddpd	%%zmm3,%%zmm6,%%zmm6		\n\t"\
		"vmovaps	%%zmm4,     (%%rsi)	\n\t"\
		"vmovaps	%%zmm7,0x080(%%rsi)	\n\t"\
		"vmovaps	%%zmm5,0x040(%%rsi)	\n\t"\
		"vmovaps	%%zmm6,0x1c0(%%rsi)	\n\t"\
		"addq	%%rdi,%%rax			\n\t"\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	%%rdi,%%rcx			\n\t"\
		"addq	%%rdi,%%rdx			\n\t"\
		/* SSE2_RADIX4_DIT_0TWIDDLE_B(r9 ): */\
		"vmovaps	     (%%rax),%%zmm2	\n\t"\
		"vmovaps	     (%%rcx),%%zmm6	\n\t"\
		"vmovaps	0x040(%%rax),%%zmm3	\n\t"\
		"vmovaps	0x040(%%rcx),%%zmm7	\n\t"\
		"vmovaps	     (%%rbx),%%zmm0	\n\t"\
		"vmovaps	     (%%rdx),%%zmm4	\n\t"\
		"vmovaps	0x040(%%rbx),%%zmm1	\n\t"\
		"vmovaps	0x040(%%rdx),%%zmm5	\n\t"\
		"vsubpd	%%zmm0,%%zmm2,%%zmm2		\n\t"\
		"vsubpd	%%zmm4,%%zmm6,%%zmm6		\n\t"\
		"vsubpd	%%zmm1,%%zmm3,%%zmm3		\n\t"\
		"vsubpd	%%zmm5,%%zmm7,%%zmm7		\n\t"\
		"vaddpd	%%zmm0,%%zmm0,%%zmm0		\n\t"\
		"vaddpd	%%zmm4,%%zmm4,%%zmm4		\n\t"\
		"vaddpd	%%zmm1,%%zmm1,%%zmm1		\n\t"\
		"vaddpd	%%zmm5,%%zmm5,%%zmm5		\n\t"\
		"vaddpd	%%zmm2,%%zmm0,%%zmm0		\n\t"\
		"vaddpd	%%zmm6,%%zmm4,%%zmm4		\n\t"\
		"vaddpd	%%zmm3,%%zmm1,%%zmm1		\n\t"\
		"vaddpd	%%zmm7,%%zmm5,%%zmm5		\n\t"\
		"vsubpd	%%zmm4,%%zmm0,%%zmm0		\n\t"\
		"vsubpd	%%zmm7,%%zmm2,%%zmm2		\n\t"\
		"vsubpd	%%zmm5,%%zmm1,%%zmm1		\n\t"\
		"vsubpd	%%zmm6,%%zmm3,%%zmm3		\n\t"\
		"movq	%[__r9], %%rsi		\n\t"\
		"vmovaps	%%zmm0,0x100(%%rsi)	\n\t"\
		"vmovaps	%%zmm2,0x180(%%rsi)	\n\t"\
		"vmovaps	%%zmm1,0x140(%%rsi)	\n\t"\
		"vmovaps	%%zmm3,0x0c0(%%rsi)	\n\t"\
		"vaddpd	%%zmm4,%%zmm4,%%zmm4		\n\t"\
		"vaddpd	%%zmm7,%%zmm7,%%zmm7		\n\t"\
		"vaddpd	%%zmm5,%%zmm5,%%zmm5		\n\t"\
		"vaddpd	%%zmm6,%%zmm6,%%zmm6		\n\t"\
		"vaddpd	%%zmm0,%%zmm4,%%zmm4		\n\t"\
		"vaddpd	%%zmm2,%%zmm7,%%zmm7		\n\t"\
		"vaddpd	%%zmm1,%%zmm5,%%zmm5		\n\t"\
		"vaddpd	%%zmm3,%%zmm6,%%zmm6		\n\t"\
		"vmovaps	%%zmm4,     (%%rsi)	\n\t"\
		"vmovaps	%%zmm7,0x080(%%rsi)	\n\t"\
		"vmovaps	%%zmm5,0x040(%%rsi)	\n\t"\
		"vmovaps	%%zmm6,0x1c0(%%rsi)	\n\t"\
		"addq	%%rdi,%%rax			\n\t"\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	%%rdi,%%rcx			\n\t"\
		"addq	%%rdi,%%rdx			\n\t"\
		/* SSE2_RADIX4_DIT_0TWIDDLE_B(r17): */\
		"vmovaps	     (%%rax),%%zmm2	\n\t"\
		"vmovaps	     (%%rcx),%%zmm6	\n\t"\
		"vmovaps	0x040(%%rax),%%zmm3	\n\t"\
		"vmovaps	0x040(%%rcx),%%zmm7	\n\t"\
		"vmovaps	     (%%rbx),%%zmm0	\n\t"\
		"vmovaps	     (%%rdx),%%zmm4	\n\t"\
		"vmovaps	0x040(%%rbx),%%zmm1	\n\t"\
		"vmovaps	0x040(%%rdx),%%zmm5	\n\t"\
		"vsubpd	%%zmm0,%%zmm2,%%zmm2		\n\t"\
		"vsubpd	%%zmm4,%%zmm6,%%zmm6		\n\t"\
		"vsubpd	%%zmm1,%%zmm3,%%zmm3		\n\t"\
		"vsubpd	%%zmm5,%%zmm7,%%zmm7		\n\t"\
		"vaddpd	%%zmm0,%%zmm0,%%zmm0		\n\t"\
		"vaddpd	%%zmm4,%%zmm4,%%zmm4		\n\t"\
		"vaddpd	%%zmm1,%%zmm1,%%zmm1		\n\t"\
		"vaddpd	%%zmm5,%%zmm5,%%zmm5		\n\t"\
		"vaddpd	%%zmm2,%%zmm0,%%zmm0		\n\t"\
		"vaddpd	%%zmm6,%%zmm4,%%zmm4		\n\t"\
		"vaddpd	%%zmm3,%%zmm1,%%zmm1		\n\t"\
		"vaddpd	%%zmm7,%%zmm5,%%zmm5		\n\t"\
		"vsubpd	%%zmm4,%%zmm0,%%zmm0		\n\t"\
		"vsubpd	%%zmm7,%%zmm2,%%zmm2		\n\t"\
		"vsubpd	%%zmm5,%%zmm1,%%zmm1		\n\t"\
		"vsubpd	%%zmm6,%%zmm3,%%zmm3		\n\t"\
		"movq	%[__r17], %%rsi		\n\t"\
		"vmovaps	%%zmm0,0x100(%%rsi)	\n\t"\
		"vmovaps	%%zmm2,0x180(%%rsi)	\n\t"\
		"vmovaps	%%zmm1,0x140(%%rsi)	\n\t"\
		"vmovaps	%%zmm3,0x0c0(%%rsi)	\n\t"\
		"vaddpd	%%zmm4,%%zmm4,%%zmm4		\n\t"\
		"vaddpd	%%zmm7,%%zmm7,%%zmm7		\n\t"\
		"vaddpd	%%zmm5,%%zmm5,%%zmm5		\n\t"\
		"vaddpd	%%zmm6,%%zmm6,%%zmm6		\n\t"\
		"vaddpd	%%zmm0,%%zmm4,%%zmm4		\n\t"\
		"vaddpd	%%zmm2,%%zmm7,%%zmm7		\n\t"\
		"vaddpd	%%zmm1,%%zmm5,%%zmm5		\n\t"\
		"vaddpd	%%zmm3,%%zmm6,%%zmm6		\n\t"\
		"vmovaps	%%zmm4,     (%%rsi)	\n\t"\
		"vmovaps	%%zmm7,0x080(%%rsi)	\n\t"\
		"vmovaps	%%zmm5,0x040(%%rsi)	\n\t"\
		"vmovaps	%%zmm6,0x1c0(%%rsi)	\n\t"\
		"addq	%%rdi,%%rax			\n\t"\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	%%rdi,%%rcx			\n\t"\
		"addq	%%rdi,%%rdx			\n\t"\
		/* SSE2_RADIX4_DIT_0TWIDDLE_B(r25): */\
		"vmovaps	     (%%rax),%%zmm2	\n\t"\
		"vmovaps	     (%%rcx),%%zmm6	\n\t"\
		"vmovaps	0x040(%%rax),%%zmm3	\n\t"\
		"vmovaps	0x040(%%rcx),%%zmm7	\n\t"\
		"vmovaps	     (%%rbx),%%zmm0	\n\t"\
		"vmovaps	     (%%rdx),%%zmm4	\n\t"\
		"vmovaps	0x040(%%rbx),%%zmm1	\n\t"\
		"vmovaps	0x040(%%rdx),%%zmm5	\n\t"\
		"vsubpd	%%zmm0,%%zmm2,%%zmm2		\n\t"\
		"vsubpd	%%zmm4,%%zmm6,%%zmm6		\n\t"\
		"vsubpd	%%zmm1,%%zmm3,%%zmm3		\n\t"\
		"vsubpd	%%zmm5,%%zmm7,%%zmm7		\n\t"\
		"vaddpd	%%zmm0,%%zmm0,%%zmm0		\n\t"\
		"vaddpd	%%zmm4,%%zmm4,%%zmm4		\n\t"\
		"vaddpd	%%zmm1,%%zmm1,%%zmm1		\n\t"\
		"vaddpd	%%zmm5,%%zmm5,%%zmm5		\n\t"\
		"vaddpd	%%zmm2,%%zmm0,%%zmm0		\n\t"\
		"vaddpd	%%zmm6,%%zmm4,%%zmm4		\n\t"\
		"vaddpd	%%zmm3,%%zmm1,%%zmm1		\n\t"\
		"vaddpd	%%zmm7,%%zmm5,%%zmm5		\n\t"\
		"vsubpd	%%zmm4,%%zmm0,%%zmm0		\n\t"\
		"vsubpd	%%zmm7,%%zmm2,%%zmm2		\n\t"\
		"vsubpd	%%zmm5,%%zmm1,%%zmm1		\n\t"\
		"vsubpd	%%zmm6,%%zmm3,%%zmm3		\n\t"\
		"movq	%[__r25], %%rsi		\n\t"\
		"vmovaps	%%zmm0,0x100(%%rsi)	\n\t"\
		"vmovaps	%%zmm2,0x180(%%rsi)	\n\t"\
		"vmovaps	%%zmm1,0x140(%%rsi)	\n\t"\
		"vmovaps	%%zmm3,0x0c0(%%rsi)	\n\t"\
		"vaddpd	%%zmm4,%%zmm4,%%zmm4		\n\t"\
		"vaddpd	%%zmm7,%%zmm7,%%zmm7		\n\t"\
		"vaddpd	%%zmm5,%%zmm5,%%zmm5		\n\t"\
		"vaddpd	%%zmm6,%%zmm6,%%zmm6		\n\t"\
		"vaddpd	%%zmm0,%%zmm4,%%zmm4		\n\t"\
		"vaddpd	%%zmm2,%%zmm7,%%zmm7		\n\t"\
		"vaddpd	%%zmm1,%%zmm5,%%zmm5		\n\t"\
		"vaddpd	%%zmm3,%%zmm6,%%zmm6		\n\t"\
		"vmovaps	%%zmm4,     (%%rsi)	\n\t"\
		"vmovaps	%%zmm7,0x080(%%rsi)	\n\t"\
		"vmovaps	%%zmm5,0x040(%%rsi)	\n\t"\
		"vmovaps	%%zmm6,0x1c0(%%rsi)	\n\t"\
	/*** Now do four pass-2 4-DFTs, inputs from local store, outputs back to same: ***/\
		"movq	%[__r1],%%rax		\n\t"\
		"vmovaps	     (%%rax),%%zmm0	\n\t"\
		"vmovaps	0x040(%%rax),%%zmm1	\n\t"\
		"vmovaps	0x200(%%rax),%%zmm2	\n\t"\
		"vmovaps	0x240(%%rax),%%zmm3	\n\t"\
		"vsubpd	0x200(%%rax),%%zmm0,%%zmm0	\n\t"\
		"vsubpd	0x240(%%rax),%%zmm1,%%zmm1	\n\t"\
		"vaddpd	     (%%rax),%%zmm2,%%zmm2	\n\t"\
		"vaddpd	0x040(%%rax),%%zmm3,%%zmm3	\n\t"\
		"vmovaps	0x400(%%rax),%%zmm4	\n\t"\
		"vmovaps	0x440(%%rax),%%zmm5	\n\t"\
		"vmovaps	0x600(%%rax),%%zmm6	\n\t"\
		"vmovaps	0x640(%%rax),%%zmm7	\n\t"\
		"vsubpd	0x600(%%rax),%%zmm4,%%zmm4	\n\t"\
		"vsubpd	0x640(%%rax),%%zmm5,%%zmm5	\n\t"\
		"vaddpd	0x400(%%rax),%%zmm6,%%zmm6	\n\t"\
		"vaddpd	0x440(%%rax),%%zmm7,%%zmm7	\n\t"\
		"vsubpd	      %%zmm6,%%zmm2,%%zmm2		\n\t"\
		"vsubpd	      %%zmm7,%%zmm3,%%zmm3		\n\t"\
		"vmovaps	%%zmm2,0x400(%%rax)	\n\t"\
		"vmovaps	%%zmm3,0x440(%%rax)	\n\t"\
		"vaddpd	      %%zmm6,%%zmm6,%%zmm6		\n\t"\
		"vaddpd	      %%zmm7,%%zmm7,%%zmm7		\n\t"\
		"vaddpd	      %%zmm2,%%zmm6,%%zmm6		\n\t"\
		"vaddpd	      %%zmm3,%%zmm7,%%zmm7		\n\t"\
		"vmovaps	%%zmm6,     (%%rax)	\n\t"\
		"vmovaps	%%zmm7,0x040(%%rax)	\n\t"\
		"vsubpd	      %%zmm5,%%zmm0,%%zmm0		\n\t"\
		"vsubpd	      %%zmm4,%%zmm1,%%zmm1		\n\t"\
		"vmovaps	%%zmm0,0x600(%%rax)	\n\t"\
		"vmovaps	%%zmm1,0x240(%%rax)	\n\t"\
		"vaddpd	      %%zmm5,%%zmm5,%%zmm5		\n\t"\
		"vaddpd	      %%zmm4,%%zmm4,%%zmm4		\n\t"\
		"vaddpd	      %%zmm0,%%zmm5,%%zmm5		\n\t"\
		"vaddpd	      %%zmm1,%%zmm4,%%zmm4		\n\t"\
		"vmovaps	%%zmm5,0x200(%%rax)	\n\t"\
		"vmovaps	%%zmm4,0x640(%%rax)	\n\t"\
		"movq	%[__r5],%%rax		\n\t"\
		"movq	%[__isrt2],%%rbx	\n\t"\
		"vmovaps	(%%rbx),%%zmm2		\n\t"\
		"vmovaps	0x400(%%rax),%%zmm4	\n\t"\
		"vmovaps	0x440(%%rax),%%zmm5	\n\t"\
		"vmovaps	0x600(%%rax),%%zmm0	\n\t"\
		"vmovaps	0x640(%%rax),%%zmm1	\n\t"\
		"vaddpd	0x440(%%rax),%%zmm4,%%zmm4	\n\t"\
		"vsubpd	0x400(%%rax),%%zmm5,%%zmm5	\n\t"\
		"vsubpd	0x640(%%rax),%%zmm0,%%zmm0	\n\t"\
		"vaddpd	0x600(%%rax),%%zmm1,%%zmm1	\n\t"\
		"vmulpd	      %%zmm2,%%zmm4,%%zmm4		\n\t"\
		"vmulpd	      %%zmm2,%%zmm5,%%zmm5		\n\t"\
		"vmulpd	      %%zmm2,%%zmm0,%%zmm0		\n\t"\
		"vmulpd	      %%zmm2,%%zmm1,%%zmm1		\n\t"\
		"vmovaps	%%zmm4,%%zmm6		\n\t"\
		"vmovaps	%%zmm5,%%zmm7		\n\t"\
		"vsubpd	      %%zmm0,%%zmm4,%%zmm4		\n\t"\
		"vsubpd	      %%zmm1,%%zmm5,%%zmm5		\n\t"\
		"vaddpd	      %%zmm0,%%zmm6,%%zmm6		\n\t"\
		"vaddpd	      %%zmm1,%%zmm7,%%zmm7		\n\t"\
		"vmovaps	     (%%rax),%%zmm0	\n\t"\
		"vmovaps	0x040(%%rax),%%zmm1	\n\t"\
		"vmovaps	0x200(%%rax),%%zmm2	\n\t"\
		"vmovaps	0x240(%%rax),%%zmm3	\n\t"\
		"vsubpd	0x240(%%rax),%%zmm0,%%zmm0	\n\t"\
		"vsubpd	0x200(%%rax),%%zmm1,%%zmm1	\n\t"\
		"vaddpd	     (%%rax),%%zmm3,%%zmm3	\n\t"\
		"vaddpd	0x040(%%rax),%%zmm2,%%zmm2	\n\t"\
		"vsubpd	      %%zmm4,%%zmm3,%%zmm3		\n\t"\
		"vsubpd	      %%zmm5,%%zmm1,%%zmm1		\n\t"\
		"vmovaps	%%zmm3,0x400(%%rax)	\n\t"\
		"vmovaps	%%zmm1,0x440(%%rax)	\n\t"\
		"vaddpd	      %%zmm4,%%zmm4,%%zmm4		\n\t"\
		"vaddpd	      %%zmm5,%%zmm5,%%zmm5		\n\t"\
		"vaddpd	      %%zmm3,%%zmm4,%%zmm4		\n\t"\
		"vaddpd	      %%zmm1,%%zmm5,%%zmm5		\n\t"\
		"vmovaps	%%zmm4,     (%%rax)	\n\t"\
		"vmovaps	%%zmm5,0x040(%%rax)	\n\t"\
		"vsubpd	      %%zmm7,%%zmm0,%%zmm0		\n\t"\
		"vsubpd	      %%zmm6,%%zmm2,%%zmm2		\n\t"\
		"vmovaps	%%zmm0,0x600(%%rax)	\n\t"\
		"vmovaps	%%zmm2,0x240(%%rax)	\n\t"\
		"vaddpd	      %%zmm7,%%zmm7,%%zmm7		\n\t"\
		"vaddpd	      %%zmm6,%%zmm6,%%zmm6		\n\t"\
		"vaddpd	      %%zmm0,%%zmm7,%%zmm7		\n\t"\
		"vaddpd	      %%zmm2,%%zmm6,%%zmm6		\n\t"\
		"vmovaps	%%zmm7,0x200(%%rax)	\n\t"\
		"vmovaps	%%zmm6,0x640(%%rax)	\n\t"\
		"movq	%[__r3],%%rax		\n\t"\
		"movq	%[__isrt2],%%rbx	\n\t"\
		"movq	%[__cc0],%%rcx		\n\t"\
		"vmovaps	0x600(%%rax),%%zmm0	\n\t"\
		"vmovaps	0x640(%%rax),%%zmm1	\n\t"\
		"vmovaps	0x600(%%rax),%%zmm2	\n\t"\
		"vmovaps	0x640(%%rax),%%zmm3	\n\t"\
		"vmulpd	0x040(%%rcx),%%zmm0,%%zmm0	\n\t"\
		"vmulpd	0x040(%%rcx),%%zmm1,%%zmm1	\n\t"\
		"vmulpd	     (%%rcx),%%zmm2,%%zmm2	\n\t"\
		"vmulpd	     (%%rcx),%%zmm3,%%zmm3	\n\t"\
		"vsubpd	      %%zmm2,%%zmm1,%%zmm1		\n\t"\
		"vaddpd	      %%zmm3,%%zmm0,%%zmm0		\n\t"\
		"vmovaps	0x400(%%rax),%%zmm4	\n\t"\
		"vmovaps	0x440(%%rax),%%zmm5	\n\t"\
		"vmovaps	0x400(%%rax),%%zmm6	\n\t"\
		"vmovaps	0x440(%%rax),%%zmm7	\n\t"\
		"vmulpd	     (%%rcx),%%zmm4,%%zmm4	\n\t"\
		"vmulpd	     (%%rcx),%%zmm5,%%zmm5	\n\t"\
		"vmulpd	0x040(%%rcx),%%zmm6,%%zmm6	\n\t"\
		"vmulpd	0x040(%%rcx),%%zmm7,%%zmm7	\n\t"\
		"vsubpd	      %%zmm6,%%zmm5,%%zmm5		\n\t"\
		"vaddpd	      %%zmm7,%%zmm4,%%zmm4		\n\t"\
		"vmovaps	%%zmm5,%%zmm7		\n\t"\
		"vmovaps	%%zmm4,%%zmm6		\n\t"\
		"vaddpd	      %%zmm0,%%zmm4,%%zmm4		\n\t"\
		"vaddpd	      %%zmm1,%%zmm5,%%zmm5		\n\t"\
		"vsubpd	      %%zmm0,%%zmm6,%%zmm6		\n\t"\
		"vsubpd	      %%zmm1,%%zmm7,%%zmm7		\n\t"\
		"vmovaps	0x200(%%rax),%%zmm2	\n\t"\
		"vmovaps	0x240(%%rax),%%zmm3	\n\t"\
		"vmovaps	     (%%rax),%%zmm0	\n\t"\
		"vmovaps	0x040(%%rax),%%zmm1	\n\t"\
		"vaddpd	0x240(%%rax),%%zmm2,%%zmm2	\n\t"\
		"vsubpd	0x200(%%rax),%%zmm3,%%zmm3	\n\t"\
		"vmulpd	     (%%rbx),%%zmm2,%%zmm2	\n\t"\
		"vmulpd	     (%%rbx),%%zmm3,%%zmm3	\n\t"\
		"vsubpd	      %%zmm2,%%zmm0,%%zmm0		\n\t"\
		"vsubpd	      %%zmm3,%%zmm1,%%zmm1		\n\t"\
		"vaddpd	      %%zmm2,%%zmm2,%%zmm2		\n\t"\
		"vaddpd	      %%zmm3,%%zmm3,%%zmm3		\n\t"\
		"vaddpd	      %%zmm0,%%zmm2,%%zmm2		\n\t"\
		"vaddpd	      %%zmm1,%%zmm3,%%zmm3		\n\t"\
		"vsubpd	      %%zmm4,%%zmm2,%%zmm2		\n\t"\
		"vsubpd	      %%zmm5,%%zmm3,%%zmm3		\n\t"\
		"vmovaps	%%zmm2,0x400(%%rax)	\n\t"\
		"vmovaps	%%zmm3,0x440(%%rax)	\n\t"\
		"vaddpd	      %%zmm4,%%zmm4,%%zmm4		\n\t"\
		"vaddpd	      %%zmm5,%%zmm5,%%zmm5		\n\t"\
		"vaddpd	      %%zmm2,%%zmm4,%%zmm4		\n\t"\
		"vaddpd	      %%zmm3,%%zmm5,%%zmm5		\n\t"\
		"vmovaps	%%zmm4,     (%%rax)	\n\t"\
		"vmovaps	%%zmm5,0x040(%%rax)	\n\t"\
		"vsubpd	      %%zmm7,%%zmm0,%%zmm0		\n\t"\
		"vsubpd	      %%zmm6,%%zmm1,%%zmm1		\n\t"\
		"vmovaps	%%zmm0,0x600(%%rax)	\n\t"\
		"vmovaps	%%zmm1,0x240(%%rax)	\n\t"\
		"vaddpd	      %%zmm7,%%zmm7,%%zmm7		\n\t"\
		"vaddpd	      %%zmm6,%%zmm6,%%zmm6		\n\t"\
		"vaddpd	      %%zmm0,%%zmm7,%%zmm7		\n\t"\
		"vaddpd	      %%zmm1,%%zmm6,%%zmm6		\n\t"\
		"vmovaps	%%zmm7,0x200(%%rax)	\n\t"\
		"vmovaps	%%zmm6,0x640(%%rax)	\n\t"\
		"movq	%[__r7],%%rax		\n\t"\
		"movq	%[__isrt2],%%rbx	\n\t"\
		"movq	%[__cc0],%%rcx		\n\t"\
		"vmovaps	0x600(%%rax),%%zmm0	\n\t"\
		"vmovaps	0x640(%%rax),%%zmm1	\n\t"\
		"vmovaps	0x600(%%rax),%%zmm2	\n\t"\
		"vmovaps	0x640(%%rax),%%zmm3	\n\t"\
		"vmulpd	     (%%rcx),%%zmm0,%%zmm0	\n\t"\
		"vmulpd	     (%%rcx),%%zmm1,%%zmm1	\n\t"\
		"vmulpd	0x040(%%rcx),%%zmm2,%%zmm2	\n\t"\
		"vmulpd	0x040(%%rcx),%%zmm3,%%zmm3	\n\t"\
		"vsubpd	      %%zmm2,%%zmm1,%%zmm1		\n\t"\
		"vaddpd	      %%zmm3,%%zmm0,%%zmm0		\n\t"\
		"vmovaps	0x400(%%rax),%%zmm4	\n\t"\
		"vmovaps	0x440(%%rax),%%zmm5	\n\t"\
		"vmovaps	0x400(%%rax),%%zmm6	\n\t"\
		"vmovaps	0x440(%%rax),%%zmm7	\n\t"\
		"vmulpd	0x040(%%rcx),%%zmm4,%%zmm4	\n\t"\
		"vmulpd	0x040(%%rcx),%%zmm5,%%zmm5	\n\t"\
		"vmulpd	     (%%rcx),%%zmm6,%%zmm6	\n\t"\
		"vmulpd	     (%%rcx),%%zmm7,%%zmm7	\n\t"\
		"vsubpd	      %%zmm6,%%zmm5,%%zmm5		\n\t"\
		"vaddpd	      %%zmm7,%%zmm4,%%zmm4		\n\t"\
		"vmovaps	%%zmm5,%%zmm7		\n\t"\
		"vmovaps	%%zmm4,%%zmm6		\n\t"\
		"vaddpd	      %%zmm0,%%zmm4,%%zmm4		\n\t"\
		"vaddpd	      %%zmm1,%%zmm5,%%zmm5		\n\t"\
		"vsubpd	      %%zmm0,%%zmm6,%%zmm6		\n\t"\
		"vsubpd	      %%zmm1,%%zmm7,%%zmm7		\n\t"\
		"vmovaps	0x200(%%rax),%%zmm2	\n\t"\
		"vmovaps	0x240(%%rax),%%zmm3	\n\t"\
		"vmovaps	     (%%rax),%%zmm0	\n\t"\
		"vmovaps	0x040(%%rax),%%zmm1	\n\t"\
		"vsubpd	0x240(%%rax),%%zmm2,%%zmm2	\n\t"\
		"vaddpd	0x200(%%rax),%%zmm3,%%zmm3	\n\t"\
		"vmulpd	     (%%rbx),%%zmm2,%%zmm2	\n\t"\
		"vmulpd	     (%%rbx),%%zmm3,%%zmm3	\n\t"\
		"vsubpd	      %%zmm2,%%zmm0,%%zmm0		\n\t"\
		"vsubpd	      %%zmm3,%%zmm1,%%zmm1		\n\t"\
		"vaddpd	      %%zmm2,%%zmm2,%%zmm2		\n\t"\
		"vaddpd	      %%zmm3,%%zmm3,%%zmm3		\n\t"\
		"vaddpd	      %%zmm0,%%zmm2,%%zmm2		\n\t"\
		"vaddpd	      %%zmm1,%%zmm3,%%zmm3		\n\t"\
		"vsubpd	      %%zmm6,%%zmm0,%%zmm0		\n\t"\
		"vsubpd	      %%zmm7,%%zmm1,%%zmm1		\n\t"\
		"vmovaps	%%zmm0,0x400(%%rax)	\n\t"\
		"vmovaps	%%zmm1,0x440(%%rax)	\n\t"\
		"vaddpd	      %%zmm6,%%zmm6,%%zmm6		\n\t"\
		"vaddpd	      %%zmm7,%%zmm7,%%zmm7		\n\t"\
		"vaddpd	      %%zmm0,%%zmm6,%%zmm6		\n\t"\
		"vaddpd	      %%zmm1,%%zmm7,%%zmm7		\n\t"\
		"vmovaps	%%zmm6,     (%%rax)	\n\t"\
		"vmovaps	%%zmm7,0x040(%%rax)	\n\t"\
		"vsubpd	      %%zmm5,%%zmm2,%%zmm2		\n\t"\
		"vsubpd	      %%zmm4,%%zmm3,%%zmm3		\n\t"\
		"vmovaps	%%zmm2,0x600(%%rax)	\n\t"\
		"vmovaps	%%zmm3,0x240(%%rax)	\n\t"\
		"vaddpd	      %%zmm5,%%zmm5,%%zmm5		\n\t"\
		"vaddpd	      %%zmm4,%%zmm4,%%zmm4		\n\t"\
		"vaddpd	      %%zmm2,%%zmm5,%%zmm5		\n\t"\
		"vaddpd	      %%zmm3,%%zmm4,%%zmm4		\n\t"\
		"vmovaps	%%zmm5,0x200(%%rax)	\n\t"\
		"vmovaps	%%zmm4,0x640(%%rax)	\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__p1] "m" (Xp1)\
		 ,[__p2] "m" (Xp2)\
		 ,[__p3] "m" (Xp3)\
		 ,[__p4] "m" (Xp4)\
		 ,[__p8] "m" (Xp8)\
		 ,[__r1] "m" (Xr1)\
		 ,[__r3] "m" (Xr3)\
		 ,[__r5] "m" (Xr5)\
		 ,[__r7] "m" (Xr7)\
		 ,[__r9] "m" (Xr9)\
		 ,[__r17] "m" (Xr17)\
		 ,[__r25] "m" (Xr25)\
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__cc0] "m" (Xcc0)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
	);\
	}

	#define SSE2_RADIX16_DIF_NOTWIDDLE(Xadd0,Xp1,Xp2,Xp3,Xp4,Xp8,Xr1,Xr3,Xr5,Xr7,Xr9,Xr11,Xr13,Xr15,Xr17,Xr19,Xr21,Xr23,Xr25,Xr27,Xr29,Xr31,Xisrt2,Xcc0)\
	{\
	__asm__ volatile (\
		/* SSE2_RADIX4_DIF_IN_PLACE(r1 , r17, r9 , r25): */\
		"movq	%[__r1] ,%%rax		\n\t"\
		"movq	%[__r17],%%rbx		\n\t"\
		"movq	%[__r9] ,%%rcx		\n\t"\
		"movq	%[__r25],%%rdx		\n\t"\
		"vmovaps	     (%%rax),%%zmm0	\n\t"\
		"vmovaps	0x040(%%rax),%%zmm1	\n\t"\
		"vmovaps	     (%%rax),%%zmm2	\n\t"\
		"vmovaps	0x040(%%rax),%%zmm3	\n\t"\
		"vaddpd	     (%%rbx),%%zmm0,%%zmm0	\n\t"\
		"vaddpd	0x040(%%rbx),%%zmm1,%%zmm1	\n\t"\
		"vsubpd	     (%%rbx),%%zmm2,%%zmm2	\n\t"\
		"vsubpd	0x040(%%rbx),%%zmm3,%%zmm3	\n\t"\
		"vmovaps	     (%%rcx),%%zmm4	\n\t"\
		"vmovaps	0x040(%%rcx),%%zmm5	\n\t"\
		"vmovaps	     (%%rcx),%%zmm6	\n\t"\
		"vmovaps	0x040(%%rcx),%%zmm7	\n\t"\
		"vaddpd	     (%%rdx),%%zmm4,%%zmm4	\n\t"\
		"vaddpd	0x040(%%rdx),%%zmm5,%%zmm5	\n\t"\
		"vsubpd	     (%%rdx),%%zmm6,%%zmm6	\n\t"\
		"vsubpd	0x040(%%rdx),%%zmm7,%%zmm7	\n\t"\
		"vsubpd	      %%zmm4,%%zmm0,%%zmm0		\n\t"\
		"vsubpd	      %%zmm5,%%zmm1,%%zmm1		\n\t"\
		"vmovaps	%%zmm0,     (%%rbx)	\n\t"\
		"vmovaps	%%zmm1,0x040(%%rbx)	\n\t"\
		"vaddpd	      %%zmm4,%%zmm4,%%zmm4		\n\t"\
		"vaddpd	      %%zmm5,%%zmm5,%%zmm5		\n\t"\
		"vaddpd	      %%zmm0,%%zmm4,%%zmm4		\n\t"\
		"vaddpd	      %%zmm1,%%zmm5,%%zmm5		\n\t"\
		"vmovaps	%%zmm4,     (%%rax)	\n\t"\
		"vmovaps	%%zmm5,0x040(%%rax)	\n\t"\
		"vsubpd	      %%zmm7,%%zmm2,%%zmm2		\n\t"\
		"vsubpd	      %%zmm6,%%zmm3,%%zmm3		\n\t"\
		"vmovaps	%%zmm2,     (%%rcx)	\n\t"\
		"vmovaps	%%zmm3,0x040(%%rdx)	\n\t"\
		"vaddpd	      %%zmm7,%%zmm7,%%zmm7		\n\t"\
		"vaddpd	      %%zmm6,%%zmm6,%%zmm6		\n\t"\
		"vaddpd	      %%zmm2,%%zmm7,%%zmm7		\n\t"\
		"vaddpd	      %%zmm3,%%zmm6,%%zmm6		\n\t"\
		"vmovaps	%%zmm7,     (%%rdx)	\n\t"\
		"vmovaps	%%zmm6,0x040(%%rcx)	\n\t"\
		/* SSE2_RADIX4_DIF_IN_PLACE(r5 , r21, r13, r29): */\
		"movq	%[__r5] ,%%rax		\n\t"\
		"movq	%[__r21],%%rbx		\n\t"\
		"movq	%[__r13],%%rcx		\n\t"\
		"movq	%[__r29],%%rdx		\n\t"\
		"vmovaps	     (%%rax),%%zmm0	\n\t"\
		"vmovaps	0x040(%%rax),%%zmm1	\n\t"\
		"vmovaps	     (%%rax),%%zmm2	\n\t"\
		"vmovaps	0x040(%%rax),%%zmm3	\n\t"\
		"vaddpd	     (%%rbx),%%zmm0,%%zmm0	\n\t"\
		"vaddpd	0x040(%%rbx),%%zmm1,%%zmm1	\n\t"\
		"vsubpd	     (%%rbx),%%zmm2,%%zmm2	\n\t"\
		"vsubpd	0x040(%%rbx),%%zmm3,%%zmm3	\n\t"\
		"vmovaps	     (%%rcx),%%zmm4	\n\t"\
		"vmovaps	0x040(%%rcx),%%zmm5	\n\t"\
		"vmovaps	     (%%rcx),%%zmm6	\n\t"\
		"vmovaps	0x040(%%rcx),%%zmm7	\n\t"\
		"vaddpd	     (%%rdx),%%zmm4,%%zmm4	\n\t"\
		"vaddpd	0x040(%%rdx),%%zmm5,%%zmm5	\n\t"\
		"vsubpd	     (%%rdx),%%zmm6,%%zmm6	\n\t"\
		"vsubpd	0x040(%%rdx),%%zmm7,%%zmm7	\n\t"\
		"vsubpd	      %%zmm4,%%zmm0,%%zmm0		\n\t"\
		"vsubpd	      %%zmm5,%%zmm1,%%zmm1		\n\t"\
		"vmovaps	%%zmm0,     (%%rbx)	\n\t"\
		"vmovaps	%%zmm1,0x040(%%rbx)	\n\t"\
		"vaddpd	      %%zmm4,%%zmm4,%%zmm4		\n\t"\
		"vaddpd	      %%zmm5,%%zmm5,%%zmm5		\n\t"\
		"vaddpd	      %%zmm0,%%zmm4,%%zmm4		\n\t"\
		"vaddpd	      %%zmm1,%%zmm5,%%zmm5		\n\t"\
		"vmovaps	%%zmm4,     (%%rax)	\n\t"\
		"vmovaps	%%zmm5,0x040(%%rax)	\n\t"\
		"vsubpd	      %%zmm7,%%zmm2,%%zmm2		\n\t"\
		"vsubpd	      %%zmm6,%%zmm3,%%zmm3		\n\t"\
		"vmovaps	%%zmm2,     (%%rcx)	\n\t"\
		"vmovaps	%%zmm3,0x040(%%rdx)	\n\t"\
		"vaddpd	      %%zmm7,%%zmm7,%%zmm7		\n\t"\
		"vaddpd	      %%zmm6,%%zmm6,%%zmm6		\n\t"\
		"vaddpd	      %%zmm2,%%zmm7,%%zmm7		\n\t"\
		"vaddpd	      %%zmm3,%%zmm6,%%zmm6		\n\t"\
		"vmovaps	%%zmm7,     (%%rdx)	\n\t"\
		"vmovaps	%%zmm6,0x040(%%rcx)	\n\t"\
		/* SSE2_RADIX4_DIF_IN_PLACE(r3 , r19, r11, r27): */\
		"movq	%[__r3] ,%%rax		\n\t"\
		"movq	%[__r19],%%rbx		\n\t"\
		"movq	%[__r11],%%rcx		\n\t"\
		"movq	%[__r27],%%rdx		\n\t"\
		"vmovaps	     (%%rax),%%zmm0	\n\t"\
		"vmovaps	0x040(%%rax),%%zmm1	\n\t"\
		"vmovaps	     (%%rax),%%zmm2	\n\t"\
		"vmovaps	0x040(%%rax),%%zmm3	\n\t"\
		"vaddpd	     (%%rbx),%%zmm0,%%zmm0	\n\t"\
		"vaddpd	0x040(%%rbx),%%zmm1,%%zmm1	\n\t"\
		"vsubpd	     (%%rbx),%%zmm2,%%zmm2	\n\t"\
		"vsubpd	0x040(%%rbx),%%zmm3,%%zmm3	\n\t"\
		"vmovaps	     (%%rcx),%%zmm4	\n\t"\
		"vmovaps	0x040(%%rcx),%%zmm5	\n\t"\
		"vmovaps	     (%%rcx),%%zmm6	\n\t"\
		"vmovaps	0x040(%%rcx),%%zmm7	\n\t"\
		"vaddpd	     (%%rdx),%%zmm4,%%zmm4	\n\t"\
		"vaddpd	0x040(%%rdx),%%zmm5,%%zmm5	\n\t"\
		"vsubpd	     (%%rdx),%%zmm6,%%zmm6	\n\t"\
		"vsubpd	0x040(%%rdx),%%zmm7,%%zmm7	\n\t"\
		"vsubpd	      %%zmm4,%%zmm0,%%zmm0		\n\t"\
		"vsubpd	      %%zmm5,%%zmm1,%%zmm1		\n\t"\
		"vmovaps	%%zmm0,     (%%rbx)	\n\t"\
		"vmovaps	%%zmm1,0x040(%%rbx)	\n\t"\
		"vaddpd	      %%zmm4,%%zmm4,%%zmm4		\n\t"\
		"vaddpd	      %%zmm5,%%zmm5,%%zmm5		\n\t"\
		"vaddpd	      %%zmm0,%%zmm4,%%zmm4		\n\t"\
		"vaddpd	      %%zmm1,%%zmm5,%%zmm5		\n\t"\
		"vmovaps	%%zmm4,     (%%rax)	\n\t"\
		"vmovaps	%%zmm5,0x040(%%rax)	\n\t"\
		"vsubpd	      %%zmm7,%%zmm2,%%zmm2		\n\t"\
		"vsubpd	      %%zmm6,%%zmm3,%%zmm3		\n\t"\
		"vmovaps	%%zmm2,     (%%rcx)	\n\t"\
		"vmovaps	%%zmm3,0x040(%%rdx)	\n\t"\
		"vaddpd	      %%zmm7,%%zmm7,%%zmm7		\n\t"\
		"vaddpd	      %%zmm6,%%zmm6,%%zmm6		\n\t"\
		"vaddpd	      %%zmm2,%%zmm7,%%zmm7		\n\t"\
		"vaddpd	      %%zmm3,%%zmm6,%%zmm6		\n\t"\
		"vmovaps	%%zmm7,     (%%rdx)	\n\t"\
		"vmovaps	%%zmm6,0x040(%%rcx)	\n\t"\
		/* SSE2_RADIX4_DIF_IN_PLACE(r7 , r23, r15, r31): */\
		"movq	%[__r7] ,%%rax		\n\t"\
		"movq	%[__r23],%%rbx		\n\t"\
		"movq	%[__r15],%%rcx		\n\t"\
		"movq	%[__r31],%%rdx		\n\t"\
		"vmovaps	     (%%rax),%%zmm0	\n\t"\
		"vmovaps	0x040(%%rax),%%zmm1	\n\t"\
		"vmovaps	     (%%rax),%%zmm2	\n\t"\
		"vmovaps	0x040(%%rax),%%zmm3	\n\t"\
		"vaddpd	     (%%rbx),%%zmm0,%%zmm0	\n\t"\
		"vaddpd	0x040(%%rbx),%%zmm1,%%zmm1	\n\t"\
		"vsubpd	     (%%rbx),%%zmm2,%%zmm2	\n\t"\
		"vsubpd	0x040(%%rbx),%%zmm3,%%zmm3	\n\t"\
		"vmovaps	     (%%rcx),%%zmm4	\n\t"\
		"vmovaps	0x040(%%rcx),%%zmm5	\n\t"\
		"vmovaps	     (%%rcx),%%zmm6	\n\t"\
		"vmovaps	0x040(%%rcx),%%zmm7	\n\t"\
		"vaddpd	     (%%rdx),%%zmm4,%%zmm4	\n\t"\
		"vaddpd	0x040(%%rdx),%%zmm5,%%zmm5	\n\t"\
		"vsubpd	     (%%rdx),%%zmm6,%%zmm6	\n\t"\
		"vsubpd	0x040(%%rdx),%%zmm7,%%zmm7	\n\t"\
		"vsubpd	      %%zmm4,%%zmm0,%%zmm0		\n\t"\
		"vsubpd	      %%zmm5,%%zmm1,%%zmm1		\n\t"\
		"vmovaps	%%zmm0,     (%%rbx)	\n\t"\
		"vmovaps	%%zmm1,0x040(%%rbx)	\n\t"\
		"vaddpd	      %%zmm4,%%zmm4,%%zmm4		\n\t"\
		"vaddpd	      %%zmm5,%%zmm5,%%zmm5		\n\t"\
		"vaddpd	      %%zmm0,%%zmm4,%%zmm4		\n\t"\
		"vaddpd	      %%zmm1,%%zmm5,%%zmm5		\n\t"\
		"vmovaps	%%zmm4,     (%%rax)	\n\t"\
		"vmovaps	%%zmm5,0x040(%%rax)	\n\t"\
		"vsubpd	      %%zmm7,%%zmm2,%%zmm2		\n\t"\
		"vsubpd	      %%zmm6,%%zmm3,%%zmm3		\n\t"\
		"vmovaps	%%zmm2,     (%%rcx)	\n\t"\
		"vmovaps	%%zmm3,0x040(%%rdx)	\n\t"\
		"vaddpd	      %%zmm7,%%zmm7,%%zmm7		\n\t"\
		"vaddpd	      %%zmm6,%%zmm6,%%zmm6		\n\t"\
		"vaddpd	      %%zmm2,%%zmm7,%%zmm7		\n\t"\
		"vaddpd	      %%zmm3,%%zmm6,%%zmm6		\n\t"\
		"vmovaps	%%zmm7,     (%%rdx)	\n\t"\
		"vmovaps	%%zmm6,0x040(%%rcx)	\n\t"\
	/*** Now do four pass-2 4-DFTs, inputs from local store, outputs back to main array: ***/\
		"movq	%[__add0],%%rax		\n\t"\
		"movslq	%[__p1],%%rbx		\n\t"\
		"movslq	%[__p2],%%rcx		\n\t"\
		"movslq	%[__p3],%%rdx		\n\t"\
		"movslq	%[__p4],%%rdi		\n\t"\
		"shlq	$3,%%rbx			\n\t"\
		"shlq	$3,%%rcx			\n\t"\
		"shlq	$3,%%rdx			\n\t"\
		"shlq	$3,%%rdi			\n\t"\
		"addq	%%rax,%%rbx			\n\t"\
		"addq	%%rax,%%rcx			\n\t"\
		"addq	%%rax,%%rdx			\n\t"\
		"movq	%[__r1],%%rsi		\n\t"\
		"vmovaps	     (%%rsi),%%zmm0	\n\t"\
		"vmovaps	0x040(%%rsi),%%zmm1	\n\t"\
		"vmovaps	0x100(%%rsi),%%zmm2	\n\t"\
		"vmovaps	0x140(%%rsi),%%zmm3	\n\t"\
		"vsubpd	0x100(%%rsi),%%zmm0,%%zmm0	\n\t"\
		"vsubpd	0x140(%%rsi),%%zmm1,%%zmm1	\n\t"\
		"vaddpd	     (%%rsi),%%zmm2,%%zmm2	\n\t"\
		"vaddpd	0x040(%%rsi),%%zmm3,%%zmm3	\n\t"\
		"vmovaps	0x080(%%rsi),%%zmm4	\n\t"\
		"vmovaps	0x0c0(%%rsi),%%zmm5	\n\t"\
		"vmovaps	0x180(%%rsi),%%zmm6	\n\t"\
		"vmovaps	0x1c0(%%rsi),%%zmm7	\n\t"\
		"vsubpd	0x180(%%rsi),%%zmm4,%%zmm4	\n\t"\
		"vsubpd	0x1c0(%%rsi),%%zmm5,%%zmm5	\n\t"\
		"vaddpd	0x080(%%rsi),%%zmm6,%%zmm6	\n\t"\
		"vaddpd	0x0c0(%%rsi),%%zmm7,%%zmm7	\n\t"\
		"vsubpd	      %%zmm6,%%zmm2,%%zmm2		\n\t"\
		"vsubpd	      %%zmm7,%%zmm3,%%zmm3		\n\t"\
		"vaddpd	      %%zmm6,%%zmm6,%%zmm6		\n\t"\
		"vaddpd	      %%zmm7,%%zmm7,%%zmm7		\n\t"\
		"vmovaps	%%zmm2,    (%%rbx)	\n\t"\
		"vmovaps	%%zmm3,0x040(%%rbx)	\n\t"\
		"vaddpd	      %%zmm2,%%zmm6,%%zmm6	\n\t"\
		"vaddpd	      %%zmm3,%%zmm7,%%zmm7	\n\t"\
		"vmovaps	%%zmm6,     (%%rax)	\n\t"\
		"vmovaps	%%zmm7,0x040(%%rax)	\n\t"\
		"vsubpd	      %%zmm5,%%zmm0,%%zmm0	\n\t"\
		"vsubpd	      %%zmm4,%%zmm1,%%zmm1	\n\t"\
		"vaddpd	      %%zmm5,%%zmm5,%%zmm5	\n\t"\
		"vaddpd	      %%zmm4,%%zmm4,%%zmm4	\n\t"\
		"vmovaps	%%zmm0,     (%%rcx)	\n\t"\
		"vmovaps	%%zmm1,0x040(%%rdx)	\n\t"\
		"vaddpd	      %%zmm0,%%zmm5,%%zmm5	\n\t"\
		"vaddpd	      %%zmm1,%%zmm4,%%zmm4	\n\t"\
		"vmovaps	%%zmm5,     (%%rdx)	\n\t"\
		"vmovaps	%%zmm4,0x040(%%rcx)	\n\t"\
		"addq	%%rdi,%%rax			\n\t"\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	%%rdi,%%rcx			\n\t"\
		"addq	%%rdi,%%rdx			\n\t"\
		"movq	%[__r17],%%rsi		\n\t"\
		"vmovaps	     (%%rsi),%%zmm0	\n\t"\
		"vmovaps	0x040(%%rsi),%%zmm1	\n\t"\
		"vmovaps	0x100(%%rsi),%%zmm2	\n\t"\
		"vmovaps	0x140(%%rsi),%%zmm3	\n\t"\
		"vsubpd	0x140(%%rsi),%%zmm0,%%zmm0	\n\t"\
		"vsubpd	0x100(%%rsi),%%zmm1,%%zmm1	\n\t"\
		"vaddpd	0x040(%%rsi),%%zmm2,%%zmm2	\n\t"\
		"vaddpd	     (%%rsi),%%zmm3,%%zmm3	\n\t"\
		"vmovaps	0x080(%%rsi),%%zmm4	\n\t"\
		"vmovaps	0x0c0(%%rsi),%%zmm5	\n\t"\
		"vmovaps	0x180(%%rsi),%%zmm6	\n\t"\
		"vmovaps	0x1c0(%%rsi),%%zmm7	\n\t"\
		"vsubpd	0x0c0(%%rsi),%%zmm4,%%zmm4	\n\t"\
		"vaddpd	0x080(%%rsi),%%zmm5,%%zmm5	\n\t"\
		"vaddpd	0x1c0(%%rsi),%%zmm6,%%zmm6	\n\t"\
		"vsubpd	0x180(%%rsi),%%zmm7,%%zmm7	\n\t"\
		"movq	%[__isrt2],%%rsi	\n\t"\
		"vmulpd	     (%%rsi),%%zmm4,%%zmm4		\n\t"\
		"vmulpd	     (%%rsi),%%zmm5,%%zmm5		\n\t"\
		"vmulpd	     (%%rsi),%%zmm6,%%zmm6		\n\t"\
		"vmulpd	     (%%rsi),%%zmm7,%%zmm7		\n\t"\
		"vsubpd	      %%zmm6,%%zmm4,%%zmm4		\n\t"\
		"vsubpd	      %%zmm7,%%zmm5,%%zmm5		\n\t"\
		"vaddpd	      %%zmm6,%%zmm6,%%zmm6		\n\t"\
		"vaddpd	      %%zmm7,%%zmm7,%%zmm7		\n\t"\
		"vaddpd	      %%zmm4,%%zmm6,%%zmm6		\n\t"\
		"vaddpd	      %%zmm5,%%zmm7,%%zmm7		\n\t"\
		"vsubpd	      %%zmm4,%%zmm0,%%zmm0		\n\t"\
		"vsubpd	      %%zmm5,%%zmm2,%%zmm2		\n\t"\
		"vaddpd	      %%zmm4,%%zmm4,%%zmm4		\n\t"\
		"vaddpd	      %%zmm5,%%zmm5,%%zmm5		\n\t"\
		"vmovaps	%%zmm0,    (%%rbx)	\n\t"\
		"vmovaps	%%zmm2,0x040(%%rbx)	\n\t"\
		"vaddpd	      %%zmm0,%%zmm4,%%zmm4	\n\t"\
		"vaddpd	      %%zmm2,%%zmm5,%%zmm5	\n\t"\
		"vmovaps	%%zmm4,     (%%rax)	\n\t"\
		"vmovaps	%%zmm5,0x040(%%rax)	\n\t"\
		"vsubpd	      %%zmm7,%%zmm3,%%zmm3	\n\t"\
		"vsubpd	      %%zmm6,%%zmm1,%%zmm1	\n\t"\
		"vaddpd	      %%zmm7,%%zmm7,%%zmm7	\n\t"\
		"vaddpd	      %%zmm6,%%zmm6,%%zmm6	\n\t"\
		"vmovaps	%%zmm3,     (%%rcx)	\n\t"\
		"vmovaps	%%zmm1,0x040(%%rdx)	\n\t"\
		"vaddpd	      %%zmm3,%%zmm7,%%zmm7	\n\t"\
		"vaddpd	      %%zmm1,%%zmm6,%%zmm6	\n\t"\
		"vmovaps	%%zmm7,     (%%rdx)	\n\t"\
		"vmovaps	%%zmm6,0x040(%%rcx)	\n\t"\
		"addq	%%rdi,%%rax			\n\t"\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	%%rdi,%%rcx			\n\t"\
		"addq	%%rdi,%%rdx			\n\t"\
		"movq	%[__r9] ,%%rsi		\n\t"\
		"movq	%[__cc0],%%rdi		\n\t"\
		"vmovaps	0x080(%%rsi),%%zmm4	\n\t"\
		"vmovaps	0x180(%%rsi),%%zmm6	\n\t"\
		"vmovaps	0x0c0(%%rsi),%%zmm5	\n\t"\
		"vmovaps	0x1c0(%%rsi),%%zmm7	\n\t"\
		"vmovaps	%%zmm4,%%zmm0		\n\t"\
		"vmovaps	%%zmm6,%%zmm2		\n\t"\
		"vmovaps	%%zmm5,%%zmm1		\n\t"\
		"vmovaps	%%zmm7,%%zmm3		\n\t"\
		"vmulpd	     (%%rdi),%%zmm4,%%zmm4	\n\t"\
		"vmulpd	0x040(%%rdi),%%zmm6,%%zmm6	\n\t"\
		"vmulpd	0x040(%%rdi),%%zmm1,%%zmm1	\n\t"\
		"vmulpd	     (%%rdi),%%zmm3,%%zmm3	\n\t"\
		"vmulpd	     (%%rdi),%%zmm5,%%zmm5	\n\t"\
		"vmulpd	0x040(%%rdi),%%zmm7,%%zmm7	\n\t"\
		"vmulpd	0x040(%%rdi),%%zmm0,%%zmm0	\n\t"\
		"vmulpd	     (%%rdi),%%zmm2,%%zmm2	\n\t"\
		"vsubpd	      %%zmm1,%%zmm4,%%zmm4		\n\t"\
		"vsubpd	      %%zmm3,%%zmm6,%%zmm6		\n\t"\
		"vaddpd	      %%zmm0,%%zmm5,%%zmm5		\n\t"\
		"vaddpd	      %%zmm2,%%zmm7,%%zmm7		\n\t"\
		"vsubpd	      %%zmm6,%%zmm4,%%zmm4		\n\t"\
		"vsubpd	      %%zmm7,%%zmm5,%%zmm5		\n\t"\
		"vaddpd	      %%zmm6,%%zmm6,%%zmm6		\n\t"\
		"vaddpd	      %%zmm7,%%zmm7,%%zmm7		\n\t"\
		"vaddpd	      %%zmm4,%%zmm6,%%zmm6		\n\t"\
		"vaddpd	      %%zmm5,%%zmm7,%%zmm7		\n\t"\
		"movq	%[__isrt2],%%rdi	\n\t"\
		"vmovaps	0x100(%%rsi),%%zmm2	\n\t"\
		"vmovaps	0x140(%%rsi),%%zmm3	\n\t"\
		"vsubpd	0x140(%%rsi),%%zmm2,%%zmm2	\n\t"\
		"vaddpd	0x100(%%rsi),%%zmm3,%%zmm3	\n\t"\
		"vmulpd	     (%%rdi),%%zmm2,%%zmm2		\n\t"\
		"vmulpd	     (%%rdi),%%zmm3,%%zmm3		\n\t"\
		"vmovaps	     (%%rsi),%%zmm0	\n\t"\
		"vmovaps	0x040(%%rsi),%%zmm1	\n\t"\
		"movslq	%[__p4],%%rdi		\n\t"\
		"shlq	$3,%%rdi			\n\t"\
		"vsubpd	      %%zmm2,%%zmm0,%%zmm0		\n\t"\
		"vsubpd	      %%zmm3,%%zmm1,%%zmm1		\n\t"\
		"vaddpd	     (%%rsi),%%zmm2,%%zmm2	\n\t"\
		"vaddpd	0x040(%%rsi),%%zmm3,%%zmm3	\n\t"\
		"vsubpd	      %%zmm6,%%zmm2,%%zmm2		\n\t"\
		"vsubpd	      %%zmm7,%%zmm3,%%zmm3		\n\t"\
		"vaddpd	      %%zmm6,%%zmm6,%%zmm6		\n\t"\
		"vaddpd	      %%zmm7,%%zmm7,%%zmm7		\n\t"\
		"vmovaps	%%zmm2,    (%%rbx)	\n\t"\
		"vmovaps	%%zmm3,0x040(%%rbx)	\n\t"\
		"vaddpd	      %%zmm2,%%zmm6,%%zmm6	\n\t"\
		"vaddpd	      %%zmm3,%%zmm7,%%zmm7	\n\t"\
		"vmovaps	%%zmm6,     (%%rax)	\n\t"\
		"vmovaps	%%zmm7,0x040(%%rax)	\n\t"\
		"vsubpd	      %%zmm5,%%zmm0,%%zmm0	\n\t"\
		"vsubpd	      %%zmm4,%%zmm1,%%zmm1	\n\t"\
		"vaddpd	      %%zmm5,%%zmm5,%%zmm5	\n\t"\
		"vaddpd	      %%zmm4,%%zmm4,%%zmm4	\n\t"\
		"vmovaps	%%zmm0,     (%%rcx)	\n\t"\
		"vmovaps	%%zmm1,0x040(%%rdx)	\n\t"\
		"vaddpd	      %%zmm0,%%zmm5,%%zmm5	\n\t"\
		"vaddpd	      %%zmm1,%%zmm4,%%zmm4	\n\t"\
		"vmovaps	%%zmm5,     (%%rdx)	\n\t"\
		"vmovaps	%%zmm4,0x040(%%rcx)	\n\t"\
		"addq	%%rdi,%%rax			\n\t"\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	%%rdi,%%rcx			\n\t"\
		"addq	%%rdi,%%rdx			\n\t"\
		"movq	%[__r25],%%rsi		\n\t"\
		"movq	%[__cc0],%%rdi		\n\t"\
		"vmovaps	0x080(%%rsi),%%zmm4	\n\t"\
		"vmovaps	0x180(%%rsi),%%zmm6	\n\t"\
		"vmovaps	0x0c0(%%rsi),%%zmm5	\n\t"\
		"vmovaps	0x1c0(%%rsi),%%zmm7	\n\t"\
		"vmovaps	%%zmm4,%%zmm0		\n\t"\
		"vmovaps	%%zmm6,%%zmm2		\n\t"\
		"vmovaps	%%zmm5,%%zmm1		\n\t"\
		"vmovaps	%%zmm7,%%zmm3		\n\t"\
		"vmulpd	0x040(%%rdi),%%zmm4,%%zmm4	\n\t"\
		"vmulpd	     (%%rdi),%%zmm6,%%zmm6	\n\t"\
		"vmulpd	     (%%rdi),%%zmm1,%%zmm1	\n\t"\
		"vmulpd	0x040(%%rdi),%%zmm3,%%zmm3	\n\t"\
		"vmulpd	0x040(%%rdi),%%zmm5,%%zmm5	\n\t"\
		"vmulpd	     (%%rdi),%%zmm7,%%zmm7	\n\t"\
		"vmulpd	     (%%rdi),%%zmm0,%%zmm0	\n\t"\
		"vmulpd	0x040(%%rdi),%%zmm2,%%zmm2	\n\t"\
		"vsubpd	      %%zmm1,%%zmm4,%%zmm4		\n\t"\
		"vsubpd	      %%zmm3,%%zmm6,%%zmm6		\n\t"\
		"vaddpd	      %%zmm0,%%zmm5,%%zmm5		\n\t"\
		"vaddpd	      %%zmm2,%%zmm7,%%zmm7		\n\t"\
		"vsubpd	      %%zmm6,%%zmm4,%%zmm4		\n\t"\
		"vsubpd	      %%zmm7,%%zmm5,%%zmm5		\n\t"\
		"vaddpd	      %%zmm6,%%zmm6,%%zmm6		\n\t"\
		"vaddpd	      %%zmm7,%%zmm7,%%zmm7		\n\t"\
		"vaddpd	      %%zmm4,%%zmm6,%%zmm6		\n\t"\
		"vaddpd	      %%zmm5,%%zmm7,%%zmm7		\n\t"\
		"movq	%[__isrt2],%%rdi	\n\t"\
		"vmovaps	0x100(%%rsi),%%zmm2	\n\t"\
		"vmovaps	0x140(%%rsi),%%zmm3	\n\t"\
		"vaddpd	0x140(%%rsi),%%zmm2,%%zmm2	\n\t"\
		"vsubpd	0x100(%%rsi),%%zmm3,%%zmm3	\n\t"\
		"vmulpd	     (%%rdi),%%zmm2,%%zmm2		\n\t"\
		"vmulpd	     (%%rdi),%%zmm3,%%zmm3		\n\t"\
		"vmovaps	     (%%rsi),%%zmm0	\n\t"\
		"vmovaps	0x040(%%rsi),%%zmm1	\n\t"\
		"vsubpd	      %%zmm2,%%zmm0,%%zmm0		\n\t"\
		"vsubpd	      %%zmm3,%%zmm1,%%zmm1		\n\t"\
		"vaddpd	     (%%rsi),%%zmm2,%%zmm2	\n\t"\
		"vaddpd	0x040(%%rsi),%%zmm3,%%zmm3	\n\t"\
		"vsubpd	      %%zmm4,%%zmm0,%%zmm0		\n\t"\
		"vsubpd	      %%zmm5,%%zmm1,%%zmm1		\n\t"\
		"vaddpd	      %%zmm4,%%zmm4,%%zmm4		\n\t"\
		"vaddpd	      %%zmm5,%%zmm5,%%zmm5		\n\t"\
		"vmovaps	%%zmm0,     (%%rbx)	\n\t"\
		"vmovaps	%%zmm1,0x040(%%rbx)	\n\t"\
		"vaddpd	      %%zmm0,%%zmm4,%%zmm4	\n\t"\
		"vaddpd	      %%zmm1,%%zmm5,%%zmm5	\n\t"\
		"vmovaps	%%zmm4,     (%%rax)	\n\t"\
		"vmovaps	%%zmm5,0x040(%%rax)	\n\t"\
		"vsubpd	      %%zmm7,%%zmm2,%%zmm2	\n\t"\
		"vsubpd	      %%zmm6,%%zmm3,%%zmm3	\n\t"\
		"vaddpd	      %%zmm7,%%zmm7,%%zmm7	\n\t"\
		"vaddpd	      %%zmm6,%%zmm6,%%zmm6	\n\t"\
		"vmovaps	%%zmm2,     (%%rcx)	\n\t"\
		"vmovaps	%%zmm3,0x040(%%rdx)	\n\t"\
		"vaddpd	      %%zmm2,%%zmm7,%%zmm7	\n\t"\
		"vaddpd	      %%zmm3,%%zmm6,%%zmm6	\n\t"\
		"vmovaps	%%zmm7,     (%%rdx)	\n\t"\
		"vmovaps	%%zmm6,0x040(%%rcx)	\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__p1] "m" (Xp1)\
		 ,[__p2] "m" (Xp2)\
		 ,[__p3] "m" (Xp3)\
		 ,[__p4] "m" (Xp4)\
		 ,[__p8] "m" (Xp8)\
		 ,[__r1] "m" (Xr1)\
		 ,[__r3] "m" (Xr3)\
		 ,[__r5] "m" (Xr5)\
		 ,[__r7] "m" (Xr7)\
		 ,[__r9] "m" (Xr9)\
		 ,[__r11] "m" (Xr11)\
		 ,[__r13] "m" (Xr13)\
		 ,[__r15] "m" (Xr15)\
		 ,[__r17] "m" (Xr17)\
		 ,[__r19] "m" (Xr19)\
		 ,[__r21] "m" (Xr21)\
		 ,[__r23] "m" (Xr23)\
		 ,[__r25] "m" (Xr25)\
		 ,[__r27] "m" (Xr27)\
		 ,[__r29] "m" (Xr29)\
		 ,[__r31] "m" (Xr31)\
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__cc0] "m" (Xcc0)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
	);\
	}

#elif defined(USE_AVX)	// AVX and AVX2 both use 256-bit registers

	#define SSE2_RADIX16_DIT_NOTWIDDLE(Xadd0,Xp1,Xp2,Xp3,Xp4,Xp8,Xr1,Xr3,Xr5,Xr7,Xr9,Xr11,Xr17,Xr25,Xisrt2,Xcc0)\
	{\
	__asm__ volatile (\
		"movq	%[__add0],%%rax		\n\t"\
		"movslq	%[__p1],%%rbx		\n\t"\
		"movslq	%[__p2],%%rcx		\n\t"\
		"movslq	%[__p3],%%rdx		\n\t"\
		"movslq	%[__p4],%%rdi		\n\t"\
		"shlq	$3,%%rbx			\n\t"\
		"shlq	$3,%%rcx			\n\t"\
		"shlq	$3,%%rdx			\n\t"\
		"shlq	$3,%%rdi			\n\t"\
		"addq	%%rax,%%rbx			\n\t"\
		"addq	%%rax,%%rcx			\n\t"\
		"addq	%%rax,%%rdx			\n\t"\
		"/* SSE2_RADIX4_DIT_0TWIDDLE_B(r1 ): */\n\t"\
		"vmovaps	     (%%rax),%%ymm2	\n\t"\
		"vmovaps	     (%%rcx),%%ymm6	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm3	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm7	\n\t"\
		"vmovaps	     (%%rbx),%%ymm0	\n\t"\
		"vmovaps	     (%%rdx),%%ymm4	\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm1	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm5	\n\t"\
		"vsubpd	%%ymm0,%%ymm2,%%ymm2		\n\t"\
		"vsubpd	%%ymm4,%%ymm6,%%ymm6		\n\t"\
		"vsubpd	%%ymm1,%%ymm3,%%ymm3		\n\t"\
		"vsubpd	%%ymm5,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	%%ymm0,%%ymm0,%%ymm0		\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	%%ymm1,%%ymm1,%%ymm1		\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	%%ymm2,%%ymm0,%%ymm0		\n\t"\
		"vaddpd	%%ymm6,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	%%ymm3,%%ymm1,%%ymm1		\n\t"\
		"vaddpd	%%ymm7,%%ymm5,%%ymm5		\n\t"\
		"vsubpd	%%ymm4,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	%%ymm7,%%ymm2,%%ymm2		\n\t"\
		"vsubpd	%%ymm5,%%ymm1,%%ymm1		\n\t"\
		"vsubpd	%%ymm6,%%ymm3,%%ymm3		\n\t"\
		"movq	%[__r1], %%rsi		\n\t"\
		"vmovaps	%%ymm0,0x080(%%rsi)	\n\t"\
		"vmovaps	%%ymm2,0x0c0(%%rsi)	\n\t"\
		"vmovaps	%%ymm1,0x0a0(%%rsi)	\n\t"\
		"vmovaps	%%ymm3,0x060(%%rsi)	\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	%%ymm0,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	%%ymm2,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	%%ymm1,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	%%ymm3,%%ymm6,%%ymm6		\n\t"\
		"vmovaps	%%ymm4,     (%%rsi)	\n\t"\
		"vmovaps	%%ymm7,0x040(%%rsi)	\n\t"\
		"vmovaps	%%ymm5,0x020(%%rsi)	\n\t"\
		"vmovaps	%%ymm6,0x0e0(%%rsi)	\n\t"\
		"addq	%%rdi,%%rax			\n\t"\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	%%rdi,%%rcx			\n\t"\
		"addq	%%rdi,%%rdx			\n\t"\
		"/* SSE2_RADIX4_DIT_0TWIDDLE_B(r9 ): */\n\t"\
		"vmovaps	     (%%rax),%%ymm2	\n\t"\
		"vmovaps	     (%%rcx),%%ymm6	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm3	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm7	\n\t"\
		"vmovaps	     (%%rbx),%%ymm0	\n\t"\
		"vmovaps	     (%%rdx),%%ymm4	\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm1	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm5	\n\t"\
		"vsubpd	%%ymm0,%%ymm2,%%ymm2		\n\t"\
		"vsubpd	%%ymm4,%%ymm6,%%ymm6		\n\t"\
		"vsubpd	%%ymm1,%%ymm3,%%ymm3		\n\t"\
		"vsubpd	%%ymm5,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	%%ymm0,%%ymm0,%%ymm0		\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	%%ymm1,%%ymm1,%%ymm1		\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	%%ymm2,%%ymm0,%%ymm0		\n\t"\
		"vaddpd	%%ymm6,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	%%ymm3,%%ymm1,%%ymm1		\n\t"\
		"vaddpd	%%ymm7,%%ymm5,%%ymm5		\n\t"\
		"vsubpd	%%ymm4,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	%%ymm7,%%ymm2,%%ymm2		\n\t"\
		"vsubpd	%%ymm5,%%ymm1,%%ymm1		\n\t"\
		"vsubpd	%%ymm6,%%ymm3,%%ymm3		\n\t"\
		"movq	%[__r9], %%rsi		\n\t"\
		"vmovaps	%%ymm0,0x080(%%rsi)	\n\t"\
		"vmovaps	%%ymm2,0x0c0(%%rsi)	\n\t"\
		"vmovaps	%%ymm1,0x0a0(%%rsi)	\n\t"\
		"vmovaps	%%ymm3,0x060(%%rsi)	\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	%%ymm0,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	%%ymm2,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	%%ymm1,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	%%ymm3,%%ymm6,%%ymm6		\n\t"\
		"vmovaps	%%ymm4,     (%%rsi)	\n\t"\
		"vmovaps	%%ymm7,0x040(%%rsi)	\n\t"\
		"vmovaps	%%ymm5,0x020(%%rsi)	\n\t"\
		"vmovaps	%%ymm6,0x0e0(%%rsi)	\n\t"\
		"addq	%%rdi,%%rax			\n\t"\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	%%rdi,%%rcx			\n\t"\
		"addq	%%rdi,%%rdx			\n\t"\
		"/* SSE2_RADIX4_DIT_0TWIDDLE_B(r17): */\n\t"\
		"vmovaps	     (%%rax),%%ymm2	\n\t"\
		"vmovaps	     (%%rcx),%%ymm6	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm3	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm7	\n\t"\
		"vmovaps	     (%%rbx),%%ymm0	\n\t"\
		"vmovaps	     (%%rdx),%%ymm4	\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm1	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm5	\n\t"\
		"vsubpd	%%ymm0,%%ymm2,%%ymm2		\n\t"\
		"vsubpd	%%ymm4,%%ymm6,%%ymm6		\n\t"\
		"vsubpd	%%ymm1,%%ymm3,%%ymm3		\n\t"\
		"vsubpd	%%ymm5,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	%%ymm0,%%ymm0,%%ymm0		\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	%%ymm1,%%ymm1,%%ymm1		\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	%%ymm2,%%ymm0,%%ymm0		\n\t"\
		"vaddpd	%%ymm6,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	%%ymm3,%%ymm1,%%ymm1		\n\t"\
		"vaddpd	%%ymm7,%%ymm5,%%ymm5		\n\t"\
		"vsubpd	%%ymm4,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	%%ymm7,%%ymm2,%%ymm2		\n\t"\
		"vsubpd	%%ymm5,%%ymm1,%%ymm1		\n\t"\
		"vsubpd	%%ymm6,%%ymm3,%%ymm3		\n\t"\
		"movq	%[__r17], %%rsi		\n\t"\
		"vmovaps	%%ymm0,0x080(%%rsi)	\n\t"\
		"vmovaps	%%ymm2,0x0c0(%%rsi)	\n\t"\
		"vmovaps	%%ymm1,0x0a0(%%rsi)	\n\t"\
		"vmovaps	%%ymm3,0x060(%%rsi)	\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	%%ymm0,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	%%ymm2,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	%%ymm1,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	%%ymm3,%%ymm6,%%ymm6		\n\t"\
		"vmovaps	%%ymm4,     (%%rsi)	\n\t"\
		"vmovaps	%%ymm7,0x040(%%rsi)	\n\t"\
		"vmovaps	%%ymm5,0x020(%%rsi)	\n\t"\
		"vmovaps	%%ymm6,0x0e0(%%rsi)	\n\t"\
		"addq	%%rdi,%%rax			\n\t"\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	%%rdi,%%rcx			\n\t"\
		"addq	%%rdi,%%rdx			\n\t"\
		"/* SSE2_RADIX4_DIT_0TWIDDLE_B(r25): */\n\t"\
		"vmovaps	     (%%rax),%%ymm2	\n\t"\
		"vmovaps	     (%%rcx),%%ymm6	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm3	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm7	\n\t"\
		"vmovaps	     (%%rbx),%%ymm0	\n\t"\
		"vmovaps	     (%%rdx),%%ymm4	\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm1	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm5	\n\t"\
		"vsubpd	%%ymm0,%%ymm2,%%ymm2		\n\t"\
		"vsubpd	%%ymm4,%%ymm6,%%ymm6		\n\t"\
		"vsubpd	%%ymm1,%%ymm3,%%ymm3		\n\t"\
		"vsubpd	%%ymm5,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	%%ymm0,%%ymm0,%%ymm0		\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	%%ymm1,%%ymm1,%%ymm1		\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	%%ymm2,%%ymm0,%%ymm0		\n\t"\
		"vaddpd	%%ymm6,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	%%ymm3,%%ymm1,%%ymm1		\n\t"\
		"vaddpd	%%ymm7,%%ymm5,%%ymm5		\n\t"\
		"vsubpd	%%ymm4,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	%%ymm7,%%ymm2,%%ymm2		\n\t"\
		"vsubpd	%%ymm5,%%ymm1,%%ymm1		\n\t"\
		"vsubpd	%%ymm6,%%ymm3,%%ymm3		\n\t"\
		"movq	%[__r25], %%rsi		\n\t"\
		"vmovaps	%%ymm0,0x080(%%rsi)	\n\t"\
		"vmovaps	%%ymm2,0x0c0(%%rsi)	\n\t"\
		"vmovaps	%%ymm1,0x0a0(%%rsi)	\n\t"\
		"vmovaps	%%ymm3,0x060(%%rsi)	\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	%%ymm0,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	%%ymm2,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	%%ymm1,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	%%ymm3,%%ymm6,%%ymm6		\n\t"\
		"vmovaps	%%ymm4,     (%%rsi)	\n\t"\
		"vmovaps	%%ymm7,0x040(%%rsi)	\n\t"\
		"vmovaps	%%ymm5,0x020(%%rsi)	\n\t"\
		"vmovaps	%%ymm6,0x0e0(%%rsi)	\n\t"\
		"movq	%[__r1],%%rax		\n\t"\
		"vmovaps	     (%%rax),%%ymm0	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1	\n\t"\
		"vmovaps	0x100(%%rax),%%ymm2	\n\t"\
		"vmovaps	0x120(%%rax),%%ymm3	\n\t"\
		"vsubpd	0x100(%%rax),%%ymm0,%%ymm0	\n\t"\
		"vsubpd	0x120(%%rax),%%ymm1,%%ymm1	\n\t"\
		"vaddpd	     (%%rax),%%ymm2,%%ymm2	\n\t"\
		"vaddpd	0x020(%%rax),%%ymm3,%%ymm3	\n\t"\
		"vmovaps	0x200(%%rax),%%ymm4	\n\t"\
		"vmovaps	0x220(%%rax),%%ymm5	\n\t"\
		"vmovaps	0x300(%%rax),%%ymm6	\n\t"\
		"vmovaps	0x320(%%rax),%%ymm7	\n\t"\
		"vsubpd	0x300(%%rax),%%ymm4,%%ymm4	\n\t"\
		"vsubpd	0x320(%%rax),%%ymm5,%%ymm5	\n\t"\
		"vaddpd	0x200(%%rax),%%ymm6,%%ymm6	\n\t"\
		"vaddpd	0x220(%%rax),%%ymm7,%%ymm7	\n\t"\
		"vsubpd	      %%ymm6,%%ymm2,%%ymm2		\n\t"\
		"vsubpd	      %%ymm7,%%ymm3,%%ymm3		\n\t"\
		"vmovaps	%%ymm2,0x200(%%rax)	\n\t"\
		"vmovaps	%%ymm3,0x220(%%rax)	\n\t"\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm2,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm3,%%ymm7,%%ymm7		\n\t"\
		"vmovaps	%%ymm6,     (%%rax)	\n\t"\
		"vmovaps	%%ymm7,0x020(%%rax)	\n\t"\
		"vsubpd	      %%ymm5,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	      %%ymm4,%%ymm1,%%ymm1		\n\t"\
		"vmovaps	%%ymm0,0x300(%%rax)	\n\t"\
		"vmovaps	%%ymm1,0x120(%%rax)	\n\t"\
		"vaddpd	      %%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	      %%ymm0,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm1,%%ymm4,%%ymm4		\n\t"\
		"vmovaps	%%ymm5,0x100(%%rax)	\n\t"\
		"vmovaps	%%ymm4,0x320(%%rax)	\n\t"\
		"movq	%[__r5],%%rax		\n\t"\
		"movq	%[__isrt2],%%rbx	\n\t"\
		"vmovaps	(%%rbx),%%ymm2		\n\t"\
		"vmovaps	0x200(%%rax),%%ymm4	\n\t"\
		"vmovaps	0x220(%%rax),%%ymm5	\n\t"\
		"vmovaps	0x300(%%rax),%%ymm0	\n\t"\
		"vmovaps	0x320(%%rax),%%ymm1	\n\t"\
		"vaddpd	0x220(%%rax),%%ymm4,%%ymm4	\n\t"\
		"vsubpd	0x200(%%rax),%%ymm5,%%ymm5	\n\t"\
		"vsubpd	0x320(%%rax),%%ymm0,%%ymm0	\n\t"\
		"vaddpd	0x300(%%rax),%%ymm1,%%ymm1	\n\t"\
		"vmulpd	      %%ymm2,%%ymm4,%%ymm4		\n\t"\
		"vmulpd	      %%ymm2,%%ymm5,%%ymm5		\n\t"\
		"vmulpd	      %%ymm2,%%ymm0,%%ymm0		\n\t"\
		"vmulpd	      %%ymm2,%%ymm1,%%ymm1		\n\t"\
		"vmovaps	%%ymm4,%%ymm6		\n\t"\
		"vmovaps	%%ymm5,%%ymm7		\n\t"\
		"vsubpd	      %%ymm0,%%ymm4,%%ymm4		\n\t"\
		"vsubpd	      %%ymm1,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm0,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm1,%%ymm7,%%ymm7		\n\t"\
		"vmovaps	     (%%rax),%%ymm0	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1	\n\t"\
		"vmovaps	0x100(%%rax),%%ymm2	\n\t"\
		"vmovaps	0x120(%%rax),%%ymm3	\n\t"\
		"vsubpd	0x120(%%rax),%%ymm0,%%ymm0	\n\t"\
		"vsubpd	0x100(%%rax),%%ymm1,%%ymm1	\n\t"\
		"vaddpd	     (%%rax),%%ymm3,%%ymm3	\n\t"\
		"vaddpd	0x020(%%rax),%%ymm2,%%ymm2	\n\t"\
		"vsubpd	      %%ymm4,%%ymm3,%%ymm3		\n\t"\
		"vsubpd	      %%ymm5,%%ymm1,%%ymm1		\n\t"\
		"vmovaps	%%ymm3,0x200(%%rax)	\n\t"\
		"vmovaps	%%ymm1,0x220(%%rax)	\n\t"\
		"vaddpd	      %%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	      %%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm3,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	      %%ymm1,%%ymm5,%%ymm5		\n\t"\
		"vmovaps	%%ymm4,     (%%rax)	\n\t"\
		"vmovaps	%%ymm5,0x020(%%rax)	\n\t"\
		"vsubpd	      %%ymm7,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	      %%ymm6,%%ymm2,%%ymm2		\n\t"\
		"vmovaps	%%ymm0,0x300(%%rax)	\n\t"\
		"vmovaps	%%ymm2,0x120(%%rax)	\n\t"\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm0,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm2,%%ymm6,%%ymm6		\n\t"\
		"vmovaps	%%ymm7,0x100(%%rax)	\n\t"\
		"vmovaps	%%ymm6,0x320(%%rax)	\n\t"\
		"movq	%[__r3],%%rax		\n\t"\
		"movq	%[__isrt2],%%rbx	\n\t"\
		"movq	%[__cc0],%%rcx		\n\t"\
		"vmovaps	0x300(%%rax),%%ymm0	\n\t"\
		"vmovaps	0x320(%%rax),%%ymm1	\n\t"\
		"vmovaps	0x300(%%rax),%%ymm2	\n\t"\
		"vmovaps	0x320(%%rax),%%ymm3	\n\t"\
		"vmulpd	0x020(%%rcx),%%ymm0,%%ymm0	\n\t"\
		"vmulpd	0x020(%%rcx),%%ymm1,%%ymm1	\n\t"\
		"vmulpd	     (%%rcx),%%ymm2,%%ymm2	\n\t"\
		"vmulpd	     (%%rcx),%%ymm3,%%ymm3	\n\t"\
		"vsubpd	      %%ymm2,%%ymm1,%%ymm1		\n\t"\
		"vaddpd	      %%ymm3,%%ymm0,%%ymm0		\n\t"\
		"vmovaps	0x200(%%rax),%%ymm4	\n\t"\
		"vmovaps	0x220(%%rax),%%ymm5	\n\t"\
		"vmovaps	0x200(%%rax),%%ymm6	\n\t"\
		"vmovaps	0x220(%%rax),%%ymm7	\n\t"\
		"vmulpd	     (%%rcx),%%ymm4,%%ymm4	\n\t"\
		"vmulpd	     (%%rcx),%%ymm5,%%ymm5	\n\t"\
		"vmulpd	0x020(%%rcx),%%ymm6,%%ymm6	\n\t"\
		"vmulpd	0x020(%%rcx),%%ymm7,%%ymm7	\n\t"\
		"vsubpd	      %%ymm6,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm7,%%ymm4,%%ymm4		\n\t"\
		"vmovaps	%%ymm5,%%ymm7		\n\t"\
		"vmovaps	%%ymm4,%%ymm6		\n\t"\
		"vaddpd	      %%ymm0,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	      %%ymm1,%%ymm5,%%ymm5		\n\t"\
		"vsubpd	      %%ymm0,%%ymm6,%%ymm6		\n\t"\
		"vsubpd	      %%ymm1,%%ymm7,%%ymm7		\n\t"\
		"vmovaps	0x100(%%rax),%%ymm2	\n\t"\
		"vmovaps	0x120(%%rax),%%ymm3	\n\t"\
		"vmovaps	     (%%rax),%%ymm0	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1	\n\t"\
		"vaddpd	0x120(%%rax),%%ymm2,%%ymm2	\n\t"\
		"vsubpd	0x100(%%rax),%%ymm3,%%ymm3	\n\t"\
		"vmulpd	     (%%rbx),%%ymm2,%%ymm2	\n\t"\
		"vmulpd	     (%%rbx),%%ymm3,%%ymm3	\n\t"\
		"vsubpd	      %%ymm2,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	      %%ymm3,%%ymm1,%%ymm1		\n\t"\
		"vaddpd	      %%ymm2,%%ymm2,%%ymm2		\n\t"\
		"vaddpd	      %%ymm3,%%ymm3,%%ymm3		\n\t"\
		"vaddpd	      %%ymm0,%%ymm2,%%ymm2		\n\t"\
		"vaddpd	      %%ymm1,%%ymm3,%%ymm3		\n\t"\
		"vsubpd	      %%ymm4,%%ymm2,%%ymm2		\n\t"\
		"vsubpd	      %%ymm5,%%ymm3,%%ymm3		\n\t"\
		"vmovaps	%%ymm2,0x200(%%rax)	\n\t"\
		"vmovaps	%%ymm3,0x220(%%rax)	\n\t"\
		"vaddpd	      %%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	      %%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm2,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	      %%ymm3,%%ymm5,%%ymm5		\n\t"\
		"vmovaps	%%ymm4,     (%%rax)	\n\t"\
		"vmovaps	%%ymm5,0x020(%%rax)	\n\t"\
		"vsubpd	      %%ymm7,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	      %%ymm6,%%ymm1,%%ymm1		\n\t"\
		"vmovaps	%%ymm0,0x300(%%rax)	\n\t"\
		"vmovaps	%%ymm1,0x120(%%rax)	\n\t"\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm0,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm1,%%ymm6,%%ymm6		\n\t"\
		"vmovaps	%%ymm7,0x100(%%rax)	\n\t"\
		"vmovaps	%%ymm6,0x320(%%rax)	\n\t"\
		"movq	%[__r7],%%rax		\n\t"\
		"movq	%[__isrt2],%%rbx	\n\t"\
		"movq	%[__cc0],%%rcx		\n\t"\
		"vmovaps	0x300(%%rax),%%ymm0	\n\t"\
		"vmovaps	0x320(%%rax),%%ymm1	\n\t"\
		"vmovaps	0x300(%%rax),%%ymm2	\n\t"\
		"vmovaps	0x320(%%rax),%%ymm3	\n\t"\
		"vmulpd	     (%%rcx),%%ymm0,%%ymm0	\n\t"\
		"vmulpd	     (%%rcx),%%ymm1,%%ymm1	\n\t"\
		"vmulpd	0x020(%%rcx),%%ymm2,%%ymm2	\n\t"\
		"vmulpd	0x020(%%rcx),%%ymm3,%%ymm3	\n\t"\
		"vsubpd	      %%ymm2,%%ymm1,%%ymm1		\n\t"\
		"vaddpd	      %%ymm3,%%ymm0,%%ymm0		\n\t"\
		"vmovaps	0x200(%%rax),%%ymm4	\n\t"\
		"vmovaps	0x220(%%rax),%%ymm5	\n\t"\
		"vmovaps	0x200(%%rax),%%ymm6	\n\t"\
		"vmovaps	0x220(%%rax),%%ymm7	\n\t"\
		"vmulpd	0x020(%%rcx),%%ymm4,%%ymm4	\n\t"\
		"vmulpd	0x020(%%rcx),%%ymm5,%%ymm5	\n\t"\
		"vmulpd	     (%%rcx),%%ymm6,%%ymm6	\n\t"\
		"vmulpd	     (%%rcx),%%ymm7,%%ymm7	\n\t"\
		"vsubpd	      %%ymm6,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm7,%%ymm4,%%ymm4		\n\t"\
		"vmovaps	%%ymm5,%%ymm7		\n\t"\
		"vmovaps	%%ymm4,%%ymm6		\n\t"\
		"vaddpd	      %%ymm0,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	      %%ymm1,%%ymm5,%%ymm5		\n\t"\
		"vsubpd	      %%ymm0,%%ymm6,%%ymm6		\n\t"\
		"vsubpd	      %%ymm1,%%ymm7,%%ymm7		\n\t"\
		"vmovaps	0x100(%%rax),%%ymm2	\n\t"\
		"vmovaps	0x120(%%rax),%%ymm3	\n\t"\
		"vmovaps	     (%%rax),%%ymm0	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1	\n\t"\
		"vsubpd	0x120(%%rax),%%ymm2,%%ymm2	\n\t"\
		"vaddpd	0x100(%%rax),%%ymm3,%%ymm3	\n\t"\
		"vmulpd	     (%%rbx),%%ymm2,%%ymm2	\n\t"\
		"vmulpd	     (%%rbx),%%ymm3,%%ymm3	\n\t"\
		"vsubpd	      %%ymm2,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	      %%ymm3,%%ymm1,%%ymm1		\n\t"\
		"vaddpd	      %%ymm2,%%ymm2,%%ymm2		\n\t"\
		"vaddpd	      %%ymm3,%%ymm3,%%ymm3		\n\t"\
		"vaddpd	      %%ymm0,%%ymm2,%%ymm2		\n\t"\
		"vaddpd	      %%ymm1,%%ymm3,%%ymm3		\n\t"\
		"vsubpd	      %%ymm6,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	      %%ymm7,%%ymm1,%%ymm1		\n\t"\
		"vmovaps	%%ymm0,0x200(%%rax)	\n\t"\
		"vmovaps	%%ymm1,0x220(%%rax)	\n\t"\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm0,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm1,%%ymm7,%%ymm7		\n\t"\
		"vmovaps	%%ymm6,     (%%rax)	\n\t"\
		"vmovaps	%%ymm7,0x020(%%rax)	\n\t"\
		"vsubpd	      %%ymm5,%%ymm2,%%ymm2		\n\t"\
		"vsubpd	      %%ymm4,%%ymm3,%%ymm3		\n\t"\
		"vmovaps	%%ymm2,0x300(%%rax)	\n\t"\
		"vmovaps	%%ymm3,0x120(%%rax)	\n\t"\
		"vaddpd	      %%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	      %%ymm2,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm3,%%ymm4,%%ymm4		\n\t"\
		"vmovaps	%%ymm5,0x100(%%rax)	\n\t"\
		"vmovaps	%%ymm4,0x320(%%rax)	\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__p1] "m" (Xp1)\
		 ,[__p2] "m" (Xp2)\
		 ,[__p3] "m" (Xp3)\
		 ,[__p4] "m" (Xp4)\
		 ,[__p8] "m" (Xp8)\
		 ,[__r1] "m" (Xr1)\
		 ,[__r3] "m" (Xr3)\
		 ,[__r5] "m" (Xr5)\
		 ,[__r7] "m" (Xr7)\
		 ,[__r9] "m" (Xr9)\
		 ,[__r17] "m" (Xr17)\
		 ,[__r25] "m" (Xr25)\
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__cc0] "m" (Xcc0)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
	);\
	}

	#define SSE2_RADIX16_DIF_NOTWIDDLE(Xadd0,Xp1,Xp2,Xp3,Xp4,Xp8,Xr1,Xr3,Xr5,Xr7,Xr9,Xr11,Xr13,Xr15,Xr17,Xr19,Xr21,Xr23,Xr25,Xr27,Xr29,Xr31,Xisrt2,Xcc0)\
	{\
	__asm__ volatile (\
		"/* SSE2_RADIX4_DIF_IN_PLACE(r1 , r17, r9 , r25): */\n\t"\
		"movq	%[__r1] ,%%rax		\n\t"\
		"movq	%[__r17],%%rbx		\n\t"\
		"movq	%[__r9] ,%%rcx		\n\t"\
		"movq	%[__r25],%%rdx		\n\t"\
		"vmovaps	     (%%rax),%%ymm0	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1	\n\t"\
		"vmovaps	     (%%rax),%%ymm2	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm3	\n\t"\
		"vaddpd	     (%%rbx),%%ymm0,%%ymm0	\n\t"\
		"vaddpd	0x020(%%rbx),%%ymm1,%%ymm1	\n\t"\
		"vsubpd	     (%%rbx),%%ymm2,%%ymm2	\n\t"\
		"vsubpd	0x020(%%rbx),%%ymm3,%%ymm3	\n\t"\
		"vmovaps	     (%%rcx),%%ymm4	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5	\n\t"\
		"vmovaps	     (%%rcx),%%ymm6	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm7	\n\t"\
		"vaddpd	     (%%rdx),%%ymm4,%%ymm4	\n\t"\
		"vaddpd	0x020(%%rdx),%%ymm5,%%ymm5	\n\t"\
		"vsubpd	     (%%rdx),%%ymm6,%%ymm6	\n\t"\
		"vsubpd	0x020(%%rdx),%%ymm7,%%ymm7	\n\t"\
		"vsubpd	      %%ymm4,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	      %%ymm5,%%ymm1,%%ymm1		\n\t"\
		"vmovaps	%%ymm0,     (%%rbx)	\n\t"\
		"vmovaps	%%ymm1,0x020(%%rbx)	\n\t"\
		"vaddpd	      %%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	      %%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm0,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	      %%ymm1,%%ymm5,%%ymm5		\n\t"\
		"vmovaps	%%ymm4,     (%%rax)	\n\t"\
		"vmovaps	%%ymm5,0x020(%%rax)	\n\t"\
		"vsubpd	      %%ymm7,%%ymm2,%%ymm2		\n\t"\
		"vsubpd	      %%ymm6,%%ymm3,%%ymm3		\n\t"\
		"vmovaps	%%ymm2,     (%%rcx)	\n\t"\
		"vmovaps	%%ymm3,0x020(%%rdx)	\n\t"\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm2,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm3,%%ymm6,%%ymm6		\n\t"\
		"vmovaps	%%ymm7,     (%%rdx)	\n\t"\
		"vmovaps	%%ymm6,0x020(%%rcx)	\n\t"\
		"/* SSE2_RADIX4_DIF_IN_PLACE(r5 , r21, r13, r29): */\n\t"\
		"movq	%[__r5] ,%%rax		\n\t"\
		"movq	%[__r21],%%rbx		\n\t"\
		"movq	%[__r13],%%rcx		\n\t"\
		"movq	%[__r29],%%rdx		\n\t"\
		"vmovaps	     (%%rax),%%ymm0	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1	\n\t"\
		"vmovaps	     (%%rax),%%ymm2	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm3	\n\t"\
		"vaddpd	     (%%rbx),%%ymm0,%%ymm0	\n\t"\
		"vaddpd	0x020(%%rbx),%%ymm1,%%ymm1	\n\t"\
		"vsubpd	     (%%rbx),%%ymm2,%%ymm2	\n\t"\
		"vsubpd	0x020(%%rbx),%%ymm3,%%ymm3	\n\t"\
		"vmovaps	     (%%rcx),%%ymm4	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5	\n\t"\
		"vmovaps	     (%%rcx),%%ymm6	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm7	\n\t"\
		"vaddpd	     (%%rdx),%%ymm4,%%ymm4	\n\t"\
		"vaddpd	0x020(%%rdx),%%ymm5,%%ymm5	\n\t"\
		"vsubpd	     (%%rdx),%%ymm6,%%ymm6	\n\t"\
		"vsubpd	0x020(%%rdx),%%ymm7,%%ymm7	\n\t"\
		"vsubpd	      %%ymm4,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	      %%ymm5,%%ymm1,%%ymm1		\n\t"\
		"vmovaps	%%ymm0,     (%%rbx)	\n\t"\
		"vmovaps	%%ymm1,0x020(%%rbx)	\n\t"\
		"vaddpd	      %%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	      %%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm0,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	      %%ymm1,%%ymm5,%%ymm5		\n\t"\
		"vmovaps	%%ymm4,     (%%rax)	\n\t"\
		"vmovaps	%%ymm5,0x020(%%rax)	\n\t"\
		"vsubpd	      %%ymm7,%%ymm2,%%ymm2		\n\t"\
		"vsubpd	      %%ymm6,%%ymm3,%%ymm3		\n\t"\
		"vmovaps	%%ymm2,     (%%rcx)	\n\t"\
		"vmovaps	%%ymm3,0x020(%%rdx)	\n\t"\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm2,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm3,%%ymm6,%%ymm6		\n\t"\
		"vmovaps	%%ymm7,     (%%rdx)	\n\t"\
		"vmovaps	%%ymm6,0x020(%%rcx)	\n\t"\
		"/* SSE2_RADIX4_DIF_IN_PLACE(r3 , r19, r11, r27): */\n\t"\
		"movq	%[__r3] ,%%rax		\n\t"\
		"movq	%[__r19],%%rbx		\n\t"\
		"movq	%[__r11],%%rcx		\n\t"\
		"movq	%[__r27],%%rdx		\n\t"\
		"vmovaps	     (%%rax),%%ymm0	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1	\n\t"\
		"vmovaps	     (%%rax),%%ymm2	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm3	\n\t"\
		"vaddpd	     (%%rbx),%%ymm0,%%ymm0	\n\t"\
		"vaddpd	0x020(%%rbx),%%ymm1,%%ymm1	\n\t"\
		"vsubpd	     (%%rbx),%%ymm2,%%ymm2	\n\t"\
		"vsubpd	0x020(%%rbx),%%ymm3,%%ymm3	\n\t"\
		"vmovaps	     (%%rcx),%%ymm4	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5	\n\t"\
		"vmovaps	     (%%rcx),%%ymm6	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm7	\n\t"\
		"vaddpd	     (%%rdx),%%ymm4,%%ymm4	\n\t"\
		"vaddpd	0x020(%%rdx),%%ymm5,%%ymm5	\n\t"\
		"vsubpd	     (%%rdx),%%ymm6,%%ymm6	\n\t"\
		"vsubpd	0x020(%%rdx),%%ymm7,%%ymm7	\n\t"\
		"vsubpd	      %%ymm4,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	      %%ymm5,%%ymm1,%%ymm1		\n\t"\
		"vmovaps	%%ymm0,     (%%rbx)	\n\t"\
		"vmovaps	%%ymm1,0x020(%%rbx)	\n\t"\
		"vaddpd	      %%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	      %%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm0,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	      %%ymm1,%%ymm5,%%ymm5		\n\t"\
		"vmovaps	%%ymm4,     (%%rax)	\n\t"\
		"vmovaps	%%ymm5,0x020(%%rax)	\n\t"\
		"vsubpd	      %%ymm7,%%ymm2,%%ymm2		\n\t"\
		"vsubpd	      %%ymm6,%%ymm3,%%ymm3		\n\t"\
		"vmovaps	%%ymm2,     (%%rcx)	\n\t"\
		"vmovaps	%%ymm3,0x020(%%rdx)	\n\t"\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm2,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm3,%%ymm6,%%ymm6		\n\t"\
		"vmovaps	%%ymm7,     (%%rdx)	\n\t"\
		"vmovaps	%%ymm6,0x020(%%rcx)	\n\t"\
		"/* SSE2_RADIX4_DIF_IN_PLACE(r7 , r23, r15, r31): */\n\t"\
		"movq	%[__r7] ,%%rax		\n\t"\
		"movq	%[__r23],%%rbx		\n\t"\
		"movq	%[__r15],%%rcx		\n\t"\
		"movq	%[__r31],%%rdx		\n\t"\
		"vmovaps	     (%%rax),%%ymm0	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1	\n\t"\
		"vmovaps	     (%%rax),%%ymm2	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm3	\n\t"\
		"vaddpd	     (%%rbx),%%ymm0,%%ymm0	\n\t"\
		"vaddpd	0x020(%%rbx),%%ymm1,%%ymm1	\n\t"\
		"vsubpd	     (%%rbx),%%ymm2,%%ymm2	\n\t"\
		"vsubpd	0x020(%%rbx),%%ymm3,%%ymm3	\n\t"\
		"vmovaps	     (%%rcx),%%ymm4	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5	\n\t"\
		"vmovaps	     (%%rcx),%%ymm6	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm7	\n\t"\
		"vaddpd	     (%%rdx),%%ymm4,%%ymm4	\n\t"\
		"vaddpd	0x020(%%rdx),%%ymm5,%%ymm5	\n\t"\
		"vsubpd	     (%%rdx),%%ymm6,%%ymm6	\n\t"\
		"vsubpd	0x020(%%rdx),%%ymm7,%%ymm7	\n\t"\
		"vsubpd	      %%ymm4,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	      %%ymm5,%%ymm1,%%ymm1		\n\t"\
		"vmovaps	%%ymm0,     (%%rbx)	\n\t"\
		"vmovaps	%%ymm1,0x020(%%rbx)	\n\t"\
		"vaddpd	      %%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	      %%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm0,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	      %%ymm1,%%ymm5,%%ymm5		\n\t"\
		"vmovaps	%%ymm4,     (%%rax)	\n\t"\
		"vmovaps	%%ymm5,0x020(%%rax)	\n\t"\
		"vsubpd	      %%ymm7,%%ymm2,%%ymm2		\n\t"\
		"vsubpd	      %%ymm6,%%ymm3,%%ymm3		\n\t"\
		"vmovaps	%%ymm2,     (%%rcx)	\n\t"\
		"vmovaps	%%ymm3,0x020(%%rdx)	\n\t"\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm2,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm3,%%ymm6,%%ymm6		\n\t"\
		"vmovaps	%%ymm7,     (%%rdx)	\n\t"\
		"vmovaps	%%ymm6,0x020(%%rcx)	\n\t"\
		"movq	%[__add0],%%rax		\n\t"\
		"movslq	%[__p1],%%rbx		\n\t"\
		"movslq	%[__p2],%%rcx		\n\t"\
		"movslq	%[__p3],%%rdx		\n\t"\
		"movslq	%[__p4],%%rdi		\n\t"\
		"shlq	$3,%%rbx			\n\t"\
		"shlq	$3,%%rcx			\n\t"\
		"shlq	$3,%%rdx			\n\t"\
		"shlq	$3,%%rdi			\n\t"\
		"addq	%%rax,%%rbx			\n\t"\
		"addq	%%rax,%%rcx			\n\t"\
		"addq	%%rax,%%rdx			\n\t"\
		"movq	%[__r1],%%rsi		\n\t"\
		"vmovaps	     (%%rsi),%%ymm0	\n\t"\
		"vmovaps	0x020(%%rsi),%%ymm1	\n\t"\
		"vmovaps	0x080(%%rsi),%%ymm2	\n\t"\
		"vmovaps	0x0a0(%%rsi),%%ymm3	\n\t"\
		"vsubpd	0x080(%%rsi),%%ymm0,%%ymm0	\n\t"\
		"vsubpd	0x0a0(%%rsi),%%ymm1,%%ymm1	\n\t"\
		"vaddpd	     (%%rsi),%%ymm2,%%ymm2	\n\t"\
		"vaddpd	0x020(%%rsi),%%ymm3,%%ymm3	\n\t"\
		"vmovaps	0x040(%%rsi),%%ymm4	\n\t"\
		"vmovaps	0x060(%%rsi),%%ymm5	\n\t"\
		"vmovaps	0x0c0(%%rsi),%%ymm6	\n\t"\
		"vmovaps	0x0e0(%%rsi),%%ymm7	\n\t"\
		"vsubpd	0x0c0(%%rsi),%%ymm4,%%ymm4	\n\t"\
		"vsubpd	0x0e0(%%rsi),%%ymm5,%%ymm5	\n\t"\
		"vaddpd	0x040(%%rsi),%%ymm6,%%ymm6	\n\t"\
		"vaddpd	0x060(%%rsi),%%ymm7,%%ymm7	\n\t"\
		"vsubpd	      %%ymm6,%%ymm2,%%ymm2		\n\t"\
		"vsubpd	      %%ymm7,%%ymm3,%%ymm3		\n\t"\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vmovaps	%%ymm2,    (%%rbx)	\n\t"\
		"vmovaps	%%ymm3,0x020(%%rbx)	\n\t"\
		"vaddpd	      %%ymm2,%%ymm6,%%ymm6	\n\t"\
		"vaddpd	      %%ymm3,%%ymm7,%%ymm7	\n\t"\
		"vmovaps	%%ymm6,     (%%rax)	\n\t"\
		"vmovaps	%%ymm7,0x020(%%rax)	\n\t"\
		"vsubpd	      %%ymm5,%%ymm0,%%ymm0	\n\t"\
		"vsubpd	      %%ymm4,%%ymm1,%%ymm1	\n\t"\
		"vaddpd	      %%ymm5,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	      %%ymm4,%%ymm4,%%ymm4	\n\t"\
		"vmovaps	%%ymm0,     (%%rcx)	\n\t"\
		"vmovaps	%%ymm1,0x020(%%rdx)	\n\t"\
		"vaddpd	      %%ymm0,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	      %%ymm1,%%ymm4,%%ymm4	\n\t"\
		"vmovaps	%%ymm5,     (%%rdx)	\n\t"\
		"vmovaps	%%ymm4,0x020(%%rcx)	\n\t"\
		"addq	%%rdi,%%rax			\n\t"\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	%%rdi,%%rcx			\n\t"\
		"addq	%%rdi,%%rdx			\n\t"\
		"movq	%[__r17],%%rsi		\n\t"\
		"vmovaps	     (%%rsi),%%ymm0	\n\t"\
		"vmovaps	0x020(%%rsi),%%ymm1	\n\t"\
		"vmovaps	0x080(%%rsi),%%ymm2	\n\t"\
		"vmovaps	0x0a0(%%rsi),%%ymm3	\n\t"\
		"vsubpd	0x0a0(%%rsi),%%ymm0,%%ymm0	\n\t"\
		"vsubpd	0x080(%%rsi),%%ymm1,%%ymm1	\n\t"\
		"vaddpd	0x020(%%rsi),%%ymm2,%%ymm2	\n\t"\
		"vaddpd	     (%%rsi),%%ymm3,%%ymm3	\n\t"\
		"vmovaps	0x040(%%rsi),%%ymm4	\n\t"\
		"vmovaps	0x060(%%rsi),%%ymm5	\n\t"\
		"vmovaps	0x0c0(%%rsi),%%ymm6	\n\t"\
		"vmovaps	0x0e0(%%rsi),%%ymm7	\n\t"\
		"vsubpd	0x060(%%rsi),%%ymm4,%%ymm4	\n\t"\
		"vaddpd	0x040(%%rsi),%%ymm5,%%ymm5	\n\t"\
		"vaddpd	0x0e0(%%rsi),%%ymm6,%%ymm6	\n\t"\
		"vsubpd	0x0c0(%%rsi),%%ymm7,%%ymm7	\n\t"\
		"movq	%[__isrt2],%%rsi	\n\t"\
		"vmulpd	     (%%rsi),%%ymm4,%%ymm4		\n\t"\
		"vmulpd	     (%%rsi),%%ymm5,%%ymm5		\n\t"\
		"vmulpd	     (%%rsi),%%ymm6,%%ymm6		\n\t"\
		"vmulpd	     (%%rsi),%%ymm7,%%ymm7		\n\t"\
		"vsubpd	      %%ymm6,%%ymm4,%%ymm4		\n\t"\
		"vsubpd	      %%ymm7,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm4,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm5,%%ymm7,%%ymm7		\n\t"\
		"vsubpd	      %%ymm4,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	      %%ymm5,%%ymm2,%%ymm2		\n\t"\
		"vaddpd	      %%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	      %%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vmovaps	%%ymm0,    (%%rbx)	\n\t"\
		"vmovaps	%%ymm2,0x020(%%rbx)	\n\t"\
		"vaddpd	      %%ymm0,%%ymm4,%%ymm4	\n\t"\
		"vaddpd	      %%ymm2,%%ymm5,%%ymm5	\n\t"\
		"vmovaps	%%ymm4,     (%%rax)	\n\t"\
		"vmovaps	%%ymm5,0x020(%%rax)	\n\t"\
		"vsubpd	      %%ymm7,%%ymm3,%%ymm3	\n\t"\
		"vsubpd	      %%ymm6,%%ymm1,%%ymm1	\n\t"\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7	\n\t"\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6	\n\t"\
		"vmovaps	%%ymm3,     (%%rcx)	\n\t"\
		"vmovaps	%%ymm1,0x020(%%rdx)	\n\t"\
		"vaddpd	      %%ymm3,%%ymm7,%%ymm7	\n\t"\
		"vaddpd	      %%ymm1,%%ymm6,%%ymm6	\n\t"\
		"vmovaps	%%ymm7,     (%%rdx)	\n\t"\
		"vmovaps	%%ymm6,0x020(%%rcx)	\n\t"\
		"addq	%%rdi,%%rax			\n\t"\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	%%rdi,%%rcx			\n\t"\
		"addq	%%rdi,%%rdx			\n\t"\
		"movq	%[__r9] ,%%rsi		\n\t"\
		"movq	%[__cc0],%%rdi		\n\t"\
		"vmovaps	0x040(%%rsi),%%ymm4	\n\t"\
		"vmovaps	0x0c0(%%rsi),%%ymm6	\n\t"\
		"vmovaps	0x060(%%rsi),%%ymm5	\n\t"\
		"vmovaps	0x0e0(%%rsi),%%ymm7	\n\t"\
		"vmovaps	%%ymm4,%%ymm0		\n\t"\
		"vmovaps	%%ymm6,%%ymm2		\n\t"\
		"vmovaps	%%ymm5,%%ymm1		\n\t"\
		"vmovaps	%%ymm7,%%ymm3		\n\t"\
		"vmulpd	     (%%rdi),%%ymm4,%%ymm4	\n\t"\
		"vmulpd	0x020(%%rdi),%%ymm6,%%ymm6	\n\t"\
		"vmulpd	0x020(%%rdi),%%ymm1,%%ymm1	\n\t"\
		"vmulpd	     (%%rdi),%%ymm3,%%ymm3	\n\t"\
		"vmulpd	     (%%rdi),%%ymm5,%%ymm5	\n\t"\
		"vmulpd	0x020(%%rdi),%%ymm7,%%ymm7	\n\t"\
		"vmulpd	0x020(%%rdi),%%ymm0,%%ymm0	\n\t"\
		"vmulpd	     (%%rdi),%%ymm2,%%ymm2	\n\t"\
		"vsubpd	      %%ymm1,%%ymm4,%%ymm4		\n\t"\
		"vsubpd	      %%ymm3,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm0,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm2,%%ymm7,%%ymm7		\n\t"\
		"vsubpd	      %%ymm6,%%ymm4,%%ymm4		\n\t"\
		"vsubpd	      %%ymm7,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm4,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm5,%%ymm7,%%ymm7		\n\t"\
		"movq	%[__isrt2],%%rdi	\n\t"\
		"vmovaps	0x080(%%rsi),%%ymm2	\n\t"\
		"vmovaps	0x0a0(%%rsi),%%ymm3	\n\t"\
		"vsubpd	0x0a0(%%rsi),%%ymm2,%%ymm2	\n\t"\
		"vaddpd	0x080(%%rsi),%%ymm3,%%ymm3	\n\t"\
		"vmulpd	     (%%rdi),%%ymm2,%%ymm2		\n\t"\
		"vmulpd	     (%%rdi),%%ymm3,%%ymm3		\n\t"\
		"vmovaps	     (%%rsi),%%ymm0	\n\t"\
		"vmovaps	0x020(%%rsi),%%ymm1	\n\t"\
		"movslq	%[__p4],%%rdi		\n\t"\
		"shlq	$3,%%rdi			\n\t"\
		"vsubpd	      %%ymm2,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	      %%ymm3,%%ymm1,%%ymm1		\n\t"\
		"vaddpd	     (%%rsi),%%ymm2,%%ymm2	\n\t"\
		"vaddpd	0x020(%%rsi),%%ymm3,%%ymm3	\n\t"\
		"vsubpd	      %%ymm6,%%ymm2,%%ymm2		\n\t"\
		"vsubpd	      %%ymm7,%%ymm3,%%ymm3		\n\t"\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vmovaps	%%ymm2,    (%%rbx)	\n\t"\
		"vmovaps	%%ymm3,0x020(%%rbx)	\n\t"\
		"vaddpd	      %%ymm2,%%ymm6,%%ymm6	\n\t"\
		"vaddpd	      %%ymm3,%%ymm7,%%ymm7	\n\t"\
		"vmovaps	%%ymm6,     (%%rax)	\n\t"\
		"vmovaps	%%ymm7,0x020(%%rax)	\n\t"\
		"vsubpd	      %%ymm5,%%ymm0,%%ymm0	\n\t"\
		"vsubpd	      %%ymm4,%%ymm1,%%ymm1	\n\t"\
		"vaddpd	      %%ymm5,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	      %%ymm4,%%ymm4,%%ymm4	\n\t"\
		"vmovaps	%%ymm0,     (%%rcx)	\n\t"\
		"vmovaps	%%ymm1,0x020(%%rdx)	\n\t"\
		"vaddpd	      %%ymm0,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	      %%ymm1,%%ymm4,%%ymm4	\n\t"\
		"vmovaps	%%ymm5,     (%%rdx)	\n\t"\
		"vmovaps	%%ymm4,0x020(%%rcx)	\n\t"\
		"addq	%%rdi,%%rax			\n\t"\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	%%rdi,%%rcx			\n\t"\
		"addq	%%rdi,%%rdx			\n\t"\
		"movq	%[__r25],%%rsi		\n\t"\
		"movq	%[__cc0],%%rdi		\n\t"\
		"vmovaps	0x040(%%rsi),%%ymm4	\n\t"\
		"vmovaps	0x0c0(%%rsi),%%ymm6	\n\t"\
		"vmovaps	0x060(%%rsi),%%ymm5	\n\t"\
		"vmovaps	0x0e0(%%rsi),%%ymm7	\n\t"\
		"vmovaps	%%ymm4,%%ymm0		\n\t"\
		"vmovaps	%%ymm6,%%ymm2		\n\t"\
		"vmovaps	%%ymm5,%%ymm1		\n\t"\
		"vmovaps	%%ymm7,%%ymm3		\n\t"\
		"vmulpd	0x020(%%rdi),%%ymm4,%%ymm4	\n\t"\
		"vmulpd	     (%%rdi),%%ymm6,%%ymm6	\n\t"\
		"vmulpd	     (%%rdi),%%ymm1,%%ymm1	\n\t"\
		"vmulpd	0x020(%%rdi),%%ymm3,%%ymm3	\n\t"\
		"vmulpd	0x020(%%rdi),%%ymm5,%%ymm5	\n\t"\
		"vmulpd	     (%%rdi),%%ymm7,%%ymm7	\n\t"\
		"vmulpd	     (%%rdi),%%ymm0,%%ymm0	\n\t"\
		"vmulpd	0x020(%%rdi),%%ymm2,%%ymm2	\n\t"\
		"vsubpd	      %%ymm1,%%ymm4,%%ymm4		\n\t"\
		"vsubpd	      %%ymm3,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm0,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm2,%%ymm7,%%ymm7		\n\t"\
		"vsubpd	      %%ymm6,%%ymm4,%%ymm4		\n\t"\
		"vsubpd	      %%ymm7,%%ymm5,%%ymm5		\n\t"\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7		\n\t"\
		"vaddpd	      %%ymm4,%%ymm6,%%ymm6		\n\t"\
		"vaddpd	      %%ymm5,%%ymm7,%%ymm7		\n\t"\
		"movq	%[__isrt2],%%rdi	\n\t"\
		"vmovaps	0x080(%%rsi),%%ymm2	\n\t"\
		"vmovaps	0x0a0(%%rsi),%%ymm3	\n\t"\
		"vaddpd	0x0a0(%%rsi),%%ymm2,%%ymm2	\n\t"\
		"vsubpd	0x080(%%rsi),%%ymm3,%%ymm3	\n\t"\
		"vmulpd	     (%%rdi),%%ymm2,%%ymm2		\n\t"\
		"vmulpd	     (%%rdi),%%ymm3,%%ymm3		\n\t"\
		"vmovaps	     (%%rsi),%%ymm0	\n\t"\
		"vmovaps	0x020(%%rsi),%%ymm1	\n\t"\
		"vsubpd	      %%ymm2,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	      %%ymm3,%%ymm1,%%ymm1		\n\t"\
		"vaddpd	     (%%rsi),%%ymm2,%%ymm2	\n\t"\
		"vaddpd	0x020(%%rsi),%%ymm3,%%ymm3	\n\t"\
		"vsubpd	      %%ymm4,%%ymm0,%%ymm0		\n\t"\
		"vsubpd	      %%ymm5,%%ymm1,%%ymm1		\n\t"\
		"vaddpd	      %%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vaddpd	      %%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vmovaps	%%ymm0,     (%%rbx)	\n\t"\
		"vmovaps	%%ymm1,0x020(%%rbx)	\n\t"\
		"vaddpd	      %%ymm0,%%ymm4,%%ymm4	\n\t"\
		"vaddpd	      %%ymm1,%%ymm5,%%ymm5	\n\t"\
		"vmovaps	%%ymm4,     (%%rax)	\n\t"\
		"vmovaps	%%ymm5,0x020(%%rax)	\n\t"\
		"vsubpd	      %%ymm7,%%ymm2,%%ymm2	\n\t"\
		"vsubpd	      %%ymm6,%%ymm3,%%ymm3	\n\t"\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7	\n\t"\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6	\n\t"\
		"vmovaps	%%ymm2,     (%%rcx)	\n\t"\
		"vmovaps	%%ymm3,0x020(%%rdx)	\n\t"\
		"vaddpd	      %%ymm2,%%ymm7,%%ymm7	\n\t"\
		"vaddpd	      %%ymm3,%%ymm6,%%ymm6	\n\t"\
		"vmovaps	%%ymm7,     (%%rdx)	\n\t"\
		"vmovaps	%%ymm6,0x020(%%rcx)	\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__p1] "m" (Xp1)\
		 ,[__p2] "m" (Xp2)\
		 ,[__p3] "m" (Xp3)\
		 ,[__p4] "m" (Xp4)\
		 ,[__p8] "m" (Xp8)\
		 ,[__r1] "m" (Xr1)\
		 ,[__r3] "m" (Xr3)\
		 ,[__r5] "m" (Xr5)\
		 ,[__r7] "m" (Xr7)\
		 ,[__r9] "m" (Xr9)\
		 ,[__r11] "m" (Xr11)\
		 ,[__r13] "m" (Xr13)\
		 ,[__r15] "m" (Xr15)\
		 ,[__r17] "m" (Xr17)\
		 ,[__r19] "m" (Xr19)\
		 ,[__r21] "m" (Xr21)\
		 ,[__r23] "m" (Xr23)\
		 ,[__r25] "m" (Xr25)\
		 ,[__r27] "m" (Xr27)\
		 ,[__r29] "m" (Xr29)\
		 ,[__r31] "m" (Xr31)\
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__cc0] "m" (Xcc0)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
	);\
	}

#elif OS_BITS == 64	// 64-bit SSE2

	#define SSE2_RADIX16_DIT_NOTWIDDLE(Xadd0,Xp1,Xp2,Xp3,Xp4,Xp8,Xr1,Xr3,Xr5,Xr7,Xr9,Xr11,Xr17,Xr25,Xisrt2,Xcc0)\
	{\
	__asm__ volatile (\
		"movq	%[__add0],%%rax		\n\t"\
		"movslq	%[__p1],%%rbx		\n\t"\
		"movslq	%[__p2],%%rcx		\n\t"\
		"movslq	%[__p3],%%rdx		\n\t"\
		"movslq	%[__p4],%%rdi		\n\t"\
		"shlq	$3,%%rbx			\n\t"\
		"shlq	$3,%%rcx			\n\t"\
		"shlq	$3,%%rdx			\n\t"\
		"shlq	$3,%%rdi			\n\t"\
		"addq	%%rax,%%rbx			\n\t"\
		"addq	%%rax,%%rcx			\n\t"\
		"addq	%%rax,%%rdx			\n\t"\
		"/* SSE2_RADIX4_DIT_0TWIDDLE_B(r1 ): */\n\t"\
		"movaps	    (%%rax),%%xmm2	\n\t"\
		"movaps	    (%%rcx),%%xmm6	\n\t"\
		"movaps	0x10(%%rax),%%xmm3	\n\t"\
		"movaps	0x10(%%rcx),%%xmm7	\n\t"\
		"movaps	    (%%rbx),%%xmm0	\n\t"\
		"movaps	    (%%rdx),%%xmm4	\n\t"\
		"movaps	0x10(%%rbx),%%xmm1	\n\t"\
		"movaps	0x10(%%rdx),%%xmm5	\n\t"\
		"subpd	%%xmm0,%%xmm2		\n\t"\
		"subpd	%%xmm4,%%xmm6		\n\t"\
		"subpd	%%xmm1,%%xmm3		\n\t"\
		"subpd	%%xmm5,%%xmm7		\n\t"\
		"addpd	%%xmm0,%%xmm0		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm1,%%xmm1		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm2,%%xmm0		\n\t"\
		"addpd	%%xmm6,%%xmm4		\n\t"\
		"addpd	%%xmm3,%%xmm1		\n\t"\
		"addpd	%%xmm7,%%xmm5		\n\t"\
		"movq	%[__r1], %%rsi		\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm2		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"subpd	%%xmm6,%%xmm3		\n\t"\
		"movaps	%%xmm0,0x040(%%rsi)	\n\t"\
		"movaps	%%xmm2,0x060(%%rsi)	\n\t"\
		"movaps	%%xmm1,0x050(%%rsi)	\n\t"\
		"movaps	%%xmm3,0x030(%%rsi)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"addpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	%%xmm4,     (%%rsi)	\n\t"\
		"movaps	%%xmm7,0x020(%%rsi)	\n\t"\
		"movaps	%%xmm5,0x010(%%rsi)	\n\t"\
		"movaps	%%xmm6,0x070(%%rsi)	\n\t"\
		"addq	%%rdi,%%rax			\n\t"\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	%%rdi,%%rcx			\n\t"\
		"addq	%%rdi,%%rdx			\n\t"\
		"/* SSE2_RADIX4_DIT_0TWIDDLE_B(r9 ): */\n\t"\
		"movaps	    (%%rax),%%xmm2	\n\t"\
		"movaps	    (%%rcx),%%xmm6	\n\t"\
		"movaps	0x10(%%rax),%%xmm3	\n\t"\
		"movaps	0x10(%%rcx),%%xmm7	\n\t"\
		"movaps	    (%%rbx),%%xmm0	\n\t"\
		"movaps	    (%%rdx),%%xmm4	\n\t"\
		"movaps	0x10(%%rbx),%%xmm1	\n\t"\
		"movaps	0x10(%%rdx),%%xmm5	\n\t"\
		"subpd	%%xmm0,%%xmm2		\n\t"\
		"subpd	%%xmm4,%%xmm6		\n\t"\
		"subpd	%%xmm1,%%xmm3		\n\t"\
		"subpd	%%xmm5,%%xmm7		\n\t"\
		"addpd	%%xmm0,%%xmm0		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm1,%%xmm1		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm2,%%xmm0		\n\t"\
		"addpd	%%xmm6,%%xmm4		\n\t"\
		"addpd	%%xmm3,%%xmm1		\n\t"\
		"addpd	%%xmm7,%%xmm5		\n\t"\
		"movq	%[__r9], %%rsi		\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm2		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"subpd	%%xmm6,%%xmm3		\n\t"\
		"movaps	%%xmm0,0x040(%%rsi)	\n\t"\
		"movaps	%%xmm2,0x060(%%rsi)	\n\t"\
		"movaps	%%xmm1,0x050(%%rsi)	\n\t"\
		"movaps	%%xmm3,0x030(%%rsi)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"addpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	%%xmm4,     (%%rsi)	\n\t"\
		"movaps	%%xmm7,0x020(%%rsi)	\n\t"\
		"movaps	%%xmm5,0x010(%%rsi)	\n\t"\
		"movaps	%%xmm6,0x070(%%rsi)	\n\t"\
		"addq	%%rdi,%%rax			\n\t"\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	%%rdi,%%rcx			\n\t"\
		"addq	%%rdi,%%rdx			\n\t"\
		"/* SSE2_RADIX4_DIT_0TWIDDLE_B(r17): */\n\t"\
		"movaps	    (%%rax),%%xmm2	\n\t"\
		"movaps	    (%%rcx),%%xmm6	\n\t"\
		"movaps	0x10(%%rax),%%xmm3	\n\t"\
		"movaps	0x10(%%rcx),%%xmm7	\n\t"\
		"movaps	    (%%rbx),%%xmm0	\n\t"\
		"movaps	    (%%rdx),%%xmm4	\n\t"\
		"movaps	0x10(%%rbx),%%xmm1	\n\t"\
		"movaps	0x10(%%rdx),%%xmm5	\n\t"\
		"subpd	%%xmm0,%%xmm2		\n\t"\
		"subpd	%%xmm4,%%xmm6		\n\t"\
		"subpd	%%xmm1,%%xmm3		\n\t"\
		"subpd	%%xmm5,%%xmm7		\n\t"\
		"addpd	%%xmm0,%%xmm0		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm1,%%xmm1		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm2,%%xmm0		\n\t"\
		"addpd	%%xmm6,%%xmm4		\n\t"\
		"addpd	%%xmm3,%%xmm1		\n\t"\
		"addpd	%%xmm7,%%xmm5		\n\t"\
		"movq	%[__r17], %%rsi		\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm2		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"subpd	%%xmm6,%%xmm3		\n\t"\
		"movaps	%%xmm0,0x040(%%rsi)	\n\t"\
		"movaps	%%xmm2,0x060(%%rsi)	\n\t"\
		"movaps	%%xmm1,0x050(%%rsi)	\n\t"\
		"movaps	%%xmm3,0x030(%%rsi)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"addpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	%%xmm4,     (%%rsi)	\n\t"\
		"movaps	%%xmm7,0x020(%%rsi)	\n\t"\
		"movaps	%%xmm5,0x010(%%rsi)	\n\t"\
		"movaps	%%xmm6,0x070(%%rsi)	\n\t"\
		"addq	%%rdi,%%rax			\n\t"\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	%%rdi,%%rcx			\n\t"\
		"addq	%%rdi,%%rdx			\n\t"\
		"/* SSE2_RADIX4_DIT_0TWIDDLE_B(r25): */\n\t"\
		"movaps	    (%%rax),%%xmm2	\n\t"\
		"movaps	    (%%rcx),%%xmm6	\n\t"\
		"movaps	0x10(%%rax),%%xmm3	\n\t"\
		"movaps	0x10(%%rcx),%%xmm7	\n\t"\
		"movaps	    (%%rbx),%%xmm0	\n\t"\
		"movaps	    (%%rdx),%%xmm4	\n\t"\
		"movaps	0x10(%%rbx),%%xmm1	\n\t"\
		"movaps	0x10(%%rdx),%%xmm5	\n\t"\
		"subpd	%%xmm0,%%xmm2		\n\t"\
		"subpd	%%xmm4,%%xmm6		\n\t"\
		"subpd	%%xmm1,%%xmm3		\n\t"\
		"subpd	%%xmm5,%%xmm7		\n\t"\
		"addpd	%%xmm0,%%xmm0		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm1,%%xmm1		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm2,%%xmm0		\n\t"\
		"addpd	%%xmm6,%%xmm4		\n\t"\
		"addpd	%%xmm3,%%xmm1		\n\t"\
		"addpd	%%xmm7,%%xmm5		\n\t"\
		"movq	%[__r25], %%rsi		\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm2		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"subpd	%%xmm6,%%xmm3		\n\t"\
		"movaps	%%xmm0,0x040(%%rsi)	\n\t"\
		"movaps	%%xmm2,0x060(%%rsi)	\n\t"\
		"movaps	%%xmm1,0x050(%%rsi)	\n\t"\
		"movaps	%%xmm3,0x030(%%rsi)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"addpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	%%xmm4,     (%%rsi)	\n\t"\
		"movaps	%%xmm7,0x020(%%rsi)	\n\t"\
		"movaps	%%xmm5,0x010(%%rsi)	\n\t"\
		"movaps	%%xmm6,0x070(%%rsi)	\n\t"\
	/*** Now do four pass-2 4-DFTs, inputs from local store, outputs back to same: ***/\
		"movq	%[__r1],%%rax		\n\t"\
		"movaps	     (%%rax),%%xmm0	\n\t"\
		"movaps	0x010(%%rax),%%xmm1	\n\t"\
		"movaps	0x080(%%rax),%%xmm2	\n\t"\
		"movaps	0x090(%%rax),%%xmm3	\n\t"\
		"subpd	0x080(%%rax),%%xmm0	\n\t"\
		"subpd	0x090(%%rax),%%xmm1	\n\t"\
		"addpd	     (%%rax),%%xmm2	\n\t"\
		"addpd	0x010(%%rax),%%xmm3	\n\t"\
		"movaps	0x100(%%rax),%%xmm4	\n\t"\
		"movaps	0x110(%%rax),%%xmm5	\n\t"\
		"movaps	0x180(%%rax),%%xmm6	\n\t"\
		"movaps	0x190(%%rax),%%xmm7	\n\t"\
		"subpd	0x180(%%rax),%%xmm4	\n\t"\
		"subpd	0x190(%%rax),%%xmm5	\n\t"\
		"addpd	0x100(%%rax),%%xmm6	\n\t"\
		"addpd	0x110(%%rax),%%xmm7	\n\t"\
		"subpd	%%xmm6,%%xmm2		\n\t"\
		"subpd	%%xmm7,%%xmm3		\n\t"\
		"movaps	%%xmm2,0x100(%%rax)	\n\t"\
		"movaps	%%xmm3,0x110(%%rax)	\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm2,%%xmm6		\n\t"\
		"addpd	%%xmm3,%%xmm7		\n\t"\
		"movaps	%%xmm6,     (%%rax)	\n\t"\
		"movaps	%%xmm7,0x010(%%rax)	\n\t"\
		"subpd	%%xmm5,%%xmm0		\n\t"\
		"subpd	%%xmm4,%%xmm1		\n\t"\
		"movaps	%%xmm0,0x180(%%rax)	\n\t"\
		"movaps	%%xmm1,0x090(%%rax)	\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm0,%%xmm5		\n\t"\
		"addpd	%%xmm1,%%xmm4		\n\t"\
		"movaps	%%xmm5,0x080(%%rax)	\n\t"\
		"movaps	%%xmm4,0x190(%%rax)	\n\t"\
		"movq	%[__r5],%%rax		\n\t"\
		"movq	%[__isrt2],%%rbx	\n\t"\
		"movaps	(%%rbx),%%xmm2		\n\t"\
		"movaps	0x100(%%rax),%%xmm4	\n\t"\
		"movaps	0x110(%%rax),%%xmm5	\n\t"\
		"movaps	0x180(%%rax),%%xmm0	\n\t"\
		"movaps	0x190(%%rax),%%xmm1	\n\t"\
		"addpd	0x110(%%rax),%%xmm4	\n\t"\
		"subpd	0x100(%%rax),%%xmm5	\n\t"\
		"subpd	0x190(%%rax),%%xmm0	\n\t"\
		"addpd	0x180(%%rax),%%xmm1	\n\t"\
		"mulpd	%%xmm2,%%xmm4		\n\t"\
		"mulpd	%%xmm2,%%xmm5		\n\t"\
		"mulpd	%%xmm2,%%xmm0		\n\t"\
		"mulpd	%%xmm2,%%xmm1		\n\t"\
		"movaps	%%xmm4,%%xmm6		\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t"\
		"subpd	%%xmm0,%%xmm4		\n\t"\
		"subpd	%%xmm1,%%xmm5		\n\t"\
		"addpd	%%xmm0,%%xmm6		\n\t"\
		"addpd	%%xmm1,%%xmm7		\n\t"\
		"movaps	     (%%rax),%%xmm0	\n\t"\
		"movaps	0x010(%%rax),%%xmm1	\n\t"\
		"movaps	0x080(%%rax),%%xmm2	\n\t"\
		"movaps	0x090(%%rax),%%xmm3	\n\t"\
		"subpd	0x090(%%rax),%%xmm0	\n\t"\
		"subpd	0x080(%%rax),%%xmm1	\n\t"\
		"addpd	     (%%rax),%%xmm3	\n\t"\
		"addpd	0x010(%%rax),%%xmm2	\n\t"\
		"subpd	%%xmm4,%%xmm3		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"movaps	%%xmm3,0x100(%%rax)	\n\t"\
		"movaps	%%xmm1,0x110(%%rax)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm3,%%xmm4		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"movaps	%%xmm4,     (%%rax)	\n\t"\
		"movaps	%%xmm5,0x010(%%rax)	\n\t"\
		"subpd	%%xmm7,%%xmm0		\n\t"\
		"subpd	%%xmm6,%%xmm2		\n\t"\
		"movaps	%%xmm0,0x180(%%rax)	\n\t"\
		"movaps	%%xmm2,0x090(%%rax)	\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm7		\n\t"\
		"addpd	%%xmm2,%%xmm6		\n\t"\
		"movaps	%%xmm7,0x080(%%rax)	\n\t"\
		"movaps	%%xmm6,0x190(%%rax)	\n\t"\
		"movq	%[__r3],%%rax		\n\t"\
		"movq	%[__isrt2],%%rbx	\n\t"\
		"movq	%[__cc0],%%rcx		\n\t"\
		"movaps	0x180(%%rax),%%xmm0	\n\t"\
		"movaps	0x190(%%rax),%%xmm1	\n\t"\
		"movaps	0x180(%%rax),%%xmm2	\n\t"\
		"movaps	0x190(%%rax),%%xmm3	\n\t"\
		"mulpd	0x010(%%rcx),%%xmm0	\n\t"\
		"mulpd	0x010(%%rcx),%%xmm1	\n\t"\
		"mulpd	     (%%rcx),%%xmm2	\n\t"\
		"mulpd	     (%%rcx),%%xmm3	\n\t"\
		"subpd	%%xmm2,%%xmm1		\n\t"\
		"addpd	%%xmm3,%%xmm0		\n\t"\
		"movaps	0x100(%%rax),%%xmm4	\n\t"\
		"movaps	0x110(%%rax),%%xmm5	\n\t"\
		"movaps	0x100(%%rax),%%xmm6	\n\t"\
		"movaps	0x110(%%rax),%%xmm7	\n\t"\
		"mulpd	     (%%rcx),%%xmm4	\n\t"\
		"mulpd	     (%%rcx),%%xmm5	\n\t"\
		"mulpd	0x010(%%rcx),%%xmm6	\n\t"\
		"mulpd	0x010(%%rcx),%%xmm7	\n\t"\
		"subpd	%%xmm6,%%xmm5		\n\t"\
		"addpd	%%xmm7,%%xmm4		\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t"\
		"movaps	%%xmm4,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"subpd	%%xmm0,%%xmm6		\n\t"\
		"subpd	%%xmm1,%%xmm7		\n\t"\
		"movaps	0x080(%%rax),%%xmm2	\n\t"\
		"movaps	0x090(%%rax),%%xmm3	\n\t"\
		"movaps	     (%%rax),%%xmm0	\n\t"\
		"movaps	0x010(%%rax),%%xmm1	\n\t"\
		"addpd	0x090(%%rax),%%xmm2	\n\t"\
		"subpd	0x080(%%rax),%%xmm3	\n\t"\
		"mulpd	     (%%rbx),%%xmm2	\n\t"\
		"mulpd	     (%%rbx),%%xmm3	\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"subpd	%%xmm4,%%xmm2		\n\t"\
		"subpd	%%xmm5,%%xmm3		\n\t"\
		"movaps	%%xmm2,0x100(%%rax)	\n\t"\
		"movaps	%%xmm3,0x110(%%rax)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm2,%%xmm4		\n\t"\
		"addpd	%%xmm3,%%xmm5		\n\t"\
		"movaps	%%xmm4,     (%%rax)	\n\t"\
		"movaps	%%xmm5,0x010(%%rax)	\n\t"\
		"subpd	%%xmm7,%%xmm0		\n\t"\
		"subpd	%%xmm6,%%xmm1		\n\t"\
		"movaps	%%xmm0,0x180(%%rax)	\n\t"\
		"movaps	%%xmm1,0x090(%%rax)	\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm7		\n\t"\
		"addpd	%%xmm1,%%xmm6		\n\t"\
		"movaps	%%xmm7,0x080(%%rax)	\n\t"\
		"movaps	%%xmm6,0x190(%%rax)	\n\t"\
		"movq	%[__r7],%%rax		\n\t"\
		"movq	%[__isrt2],%%rbx	\n\t"\
		"movq	%[__cc0],%%rcx		\n\t"\
		"movaps	0x180(%%rax),%%xmm0	\n\t"\
		"movaps	0x190(%%rax),%%xmm1	\n\t"\
		"movaps	0x180(%%rax),%%xmm2	\n\t"\
		"movaps	0x190(%%rax),%%xmm3	\n\t"\
		"mulpd	     (%%rcx),%%xmm0	\n\t"\
		"mulpd	     (%%rcx),%%xmm1	\n\t"\
		"mulpd	0x010(%%rcx),%%xmm2	\n\t"\
		"mulpd	0x010(%%rcx),%%xmm3	\n\t"\
		"subpd	%%xmm2,%%xmm1		\n\t"\
		"addpd	%%xmm3,%%xmm0		\n\t"\
		"movaps	0x100(%%rax),%%xmm4	\n\t"\
		"movaps	0x110(%%rax),%%xmm5	\n\t"\
		"movaps	0x100(%%rax),%%xmm6	\n\t"\
		"movaps	0x110(%%rax),%%xmm7	\n\t"\
		"mulpd	0x010(%%rcx),%%xmm4	\n\t"\
		"mulpd	0x010(%%rcx),%%xmm5	\n\t"\
		"mulpd	     (%%rcx),%%xmm6	\n\t"\
		"mulpd	     (%%rcx),%%xmm7	\n\t"\
		"subpd	%%xmm6,%%xmm5		\n\t"\
		"addpd	%%xmm7,%%xmm4		\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t"\
		"movaps	%%xmm4,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"subpd	%%xmm0,%%xmm6		\n\t"\
		"subpd	%%xmm1,%%xmm7		\n\t"\
		"movaps	0x080(%%rax),%%xmm2	\n\t"\
		"movaps	0x090(%%rax),%%xmm3	\n\t"\
		"movaps	     (%%rax),%%xmm0	\n\t"\
		"movaps	0x010(%%rax),%%xmm1	\n\t"\
		"subpd	0x090(%%rax),%%xmm2	\n\t"\
		"addpd	0x080(%%rax),%%xmm3	\n\t"\
		"mulpd	     (%%rbx),%%xmm2	\n\t"\
		"mulpd	     (%%rbx),%%xmm3	\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"subpd	%%xmm6,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm1		\n\t"\
		"movaps	%%xmm0,0x100(%%rax)	\n\t"\
		"movaps	%%xmm1,0x110(%%rax)	\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm0,%%xmm6		\n\t"\
		"addpd	%%xmm1,%%xmm7		\n\t"\
		"movaps	%%xmm6,     (%%rax)	\n\t"\
		"movaps	%%xmm7,0x010(%%rax)	\n\t"\
		"subpd	%%xmm5,%%xmm2		\n\t"\
		"subpd	%%xmm4,%%xmm3		\n\t"\
		"movaps	%%xmm2,0x180(%%rax)	\n\t"\
		"movaps	%%xmm3,0x090(%%rax)	\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm2,%%xmm5		\n\t"\
		"addpd	%%xmm3,%%xmm4		\n\t"\
		"movaps	%%xmm5,0x080(%%rax)	\n\t"\
		"movaps	%%xmm4,0x190(%%rax)	\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__p1] "m" (Xp1)\
		 ,[__p2] "m" (Xp2)\
		 ,[__p3] "m" (Xp3)\
		 ,[__p4] "m" (Xp4)\
		 ,[__p8] "m" (Xp8)\
		 ,[__r1] "m" (Xr1)\
		 ,[__r3] "m" (Xr3)\
		 ,[__r5] "m" (Xr5)\
		 ,[__r7] "m" (Xr7)\
		 ,[__r9] "m" (Xr9)\
		 ,[__r17] "m" (Xr17)\
		 ,[__r25] "m" (Xr25)\
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__cc0] "m" (Xcc0)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
	);\
	}

	#define SSE2_RADIX16_DIF_NOTWIDDLE(Xadd0,Xp1,Xp2,Xp3,Xp4,Xp8,Xr1,Xr3,Xr5,Xr7,Xr9,Xr11,Xr13,Xr15,Xr17,Xr19,Xr21,Xr23,Xr25,Xr27,Xr29,Xr31,Xisrt2,Xcc0)\
	{\
	__asm__ volatile (\
		"/* SSE2_RADIX4_DIF_IN_PLACE(r1 , r17, r9 , r25): */\n\t"\
		"movq	%[__r1] ,%%rax		\n\t"\
		"movq	%[__r17],%%rbx		\n\t"\
		"movq	%[__r9] ,%%rcx		\n\t"\
		"movq	%[__r25],%%rdx		\n\t"\
		"movaps	    (%%rax),%%xmm0	\n\t"\
		"movaps	0x10(%%rax),%%xmm1	\n\t"\
		"movaps	    (%%rax),%%xmm2	\n\t"\
		"movaps	0x10(%%rax),%%xmm3	\n\t"\
		"addpd	    (%%rbx),%%xmm0	\n\t"\
		"addpd	0x10(%%rbx),%%xmm1	\n\t"\
		"subpd	    (%%rbx),%%xmm2	\n\t"\
		"subpd	0x10(%%rbx),%%xmm3	\n\t"\
		"movaps	    (%%rcx),%%xmm4	\n\t"\
		"movaps	0x10(%%rcx),%%xmm5	\n\t"\
		"movaps	    (%%rcx),%%xmm6	\n\t"\
		"movaps	0x10(%%rcx),%%xmm7	\n\t"\
		"addpd	    (%%rdx),%%xmm4	\n\t"\
		"addpd	0x10(%%rdx),%%xmm5	\n\t"\
		"subpd	    (%%rdx),%%xmm6	\n\t"\
		"subpd	0x10(%%rdx),%%xmm7	\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"movaps	%%xmm0,     (%%rbx)	\n\t"\
		"movaps	%%xmm1,0x010(%%rbx)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"movaps	%%xmm4,     (%%rax)	\n\t"\
		"movaps	%%xmm5,0x010(%%rax)	\n\t"\
		"subpd	%%xmm7,%%xmm2		\n\t"\
		"subpd	%%xmm6,%%xmm3		\n\t"\
		"movaps	%%xmm2,     (%%rcx)	\n\t"\
		"movaps	%%xmm3,0x010(%%rdx)	\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	%%xmm7,     (%%rdx)	\n\t"\
		"movaps	%%xmm6,0x010(%%rcx)	\n\t"\
		"/* SSE2_RADIX4_DIF_IN_PLACE(r5 , r21, r13, r29): */\n\t"\
		"movq	%[__r5] ,%%rax		\n\t"\
		"movq	%[__r21],%%rbx		\n\t"\
		"movq	%[__r13],%%rcx		\n\t"\
		"movq	%[__r29],%%rdx		\n\t"\
		"movaps	    (%%rax),%%xmm0	\n\t"\
		"movaps	0x10(%%rax),%%xmm1	\n\t"\
		"movaps	    (%%rax),%%xmm2	\n\t"\
		"movaps	0x10(%%rax),%%xmm3	\n\t"\
		"addpd	    (%%rbx),%%xmm0	\n\t"\
		"addpd	0x10(%%rbx),%%xmm1	\n\t"\
		"subpd	    (%%rbx),%%xmm2	\n\t"\
		"subpd	0x10(%%rbx),%%xmm3	\n\t"\
		"movaps	    (%%rcx),%%xmm4	\n\t"\
		"movaps	0x10(%%rcx),%%xmm5	\n\t"\
		"movaps	    (%%rcx),%%xmm6	\n\t"\
		"movaps	0x10(%%rcx),%%xmm7	\n\t"\
		"addpd	    (%%rdx),%%xmm4	\n\t"\
		"addpd	0x10(%%rdx),%%xmm5	\n\t"\
		"subpd	    (%%rdx),%%xmm6	\n\t"\
		"subpd	0x10(%%rdx),%%xmm7	\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"movaps	%%xmm0,     (%%rbx)	\n\t"\
		"movaps	%%xmm1,0x010(%%rbx)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"movaps	%%xmm4,     (%%rax)	\n\t"\
		"movaps	%%xmm5,0x010(%%rax)	\n\t"\
		"subpd	%%xmm7,%%xmm2		\n\t"\
		"subpd	%%xmm6,%%xmm3		\n\t"\
		"movaps	%%xmm2,     (%%rcx)	\n\t"\
		"movaps	%%xmm3,0x010(%%rdx)	\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	%%xmm7,     (%%rdx)	\n\t"\
		"movaps	%%xmm6,0x010(%%rcx)	\n\t"\
		"/* SSE2_RADIX4_DIF_IN_PLACE(r3 , r19, r11, r27): */\n\t"\
		"movq	%[__r3] ,%%rax		\n\t"\
		"movq	%[__r19],%%rbx		\n\t"\
		"movq	%[__r11],%%rcx		\n\t"\
		"movq	%[__r27],%%rdx		\n\t"\
		"movaps	    (%%rax),%%xmm0	\n\t"\
		"movaps	0x10(%%rax),%%xmm1	\n\t"\
		"movaps	    (%%rax),%%xmm2	\n\t"\
		"movaps	0x10(%%rax),%%xmm3	\n\t"\
		"addpd	    (%%rbx),%%xmm0	\n\t"\
		"addpd	0x10(%%rbx),%%xmm1	\n\t"\
		"subpd	    (%%rbx),%%xmm2	\n\t"\
		"subpd	0x10(%%rbx),%%xmm3	\n\t"\
		"movaps	    (%%rcx),%%xmm4	\n\t"\
		"movaps	0x10(%%rcx),%%xmm5	\n\t"\
		"movaps	    (%%rcx),%%xmm6	\n\t"\
		"movaps	0x10(%%rcx),%%xmm7	\n\t"\
		"addpd	    (%%rdx),%%xmm4	\n\t"\
		"addpd	0x10(%%rdx),%%xmm5	\n\t"\
		"subpd	    (%%rdx),%%xmm6	\n\t"\
		"subpd	0x10(%%rdx),%%xmm7	\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"movaps	%%xmm0,     (%%rbx)	\n\t"\
		"movaps	%%xmm1,0x010(%%rbx)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"movaps	%%xmm4,     (%%rax)	\n\t"\
		"movaps	%%xmm5,0x010(%%rax)	\n\t"\
		"subpd	%%xmm7,%%xmm2		\n\t"\
		"subpd	%%xmm6,%%xmm3		\n\t"\
		"movaps	%%xmm2,     (%%rcx)	\n\t"\
		"movaps	%%xmm3,0x010(%%rdx)	\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	%%xmm7,     (%%rdx)	\n\t"\
		"movaps	%%xmm6,0x010(%%rcx)	\n\t"\
		"/* SSE2_RADIX4_DIF_IN_PLACE(r7 , r23, r15, r31): */\n\t"\
		"movq	%[__r7] ,%%rax		\n\t"\
		"movq	%[__r23],%%rbx		\n\t"\
		"movq	%[__r15],%%rcx		\n\t"\
		"movq	%[__r31],%%rdx		\n\t"\
		"movaps	    (%%rax),%%xmm0	\n\t"\
		"movaps	0x10(%%rax),%%xmm1	\n\t"\
		"movaps	    (%%rax),%%xmm2	\n\t"\
		"movaps	0x10(%%rax),%%xmm3	\n\t"\
		"addpd	    (%%rbx),%%xmm0	\n\t"\
		"addpd	0x10(%%rbx),%%xmm1	\n\t"\
		"subpd	    (%%rbx),%%xmm2	\n\t"\
		"subpd	0x10(%%rbx),%%xmm3	\n\t"\
		"movaps	    (%%rcx),%%xmm4	\n\t"\
		"movaps	0x10(%%rcx),%%xmm5	\n\t"\
		"movaps	    (%%rcx),%%xmm6	\n\t"\
		"movaps	0x10(%%rcx),%%xmm7	\n\t"\
		"addpd	    (%%rdx),%%xmm4	\n\t"\
		"addpd	0x10(%%rdx),%%xmm5	\n\t"\
		"subpd	    (%%rdx),%%xmm6	\n\t"\
		"subpd	0x10(%%rdx),%%xmm7	\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"movaps	%%xmm0,     (%%rbx)	\n\t"\
		"movaps	%%xmm1,0x010(%%rbx)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"movaps	%%xmm4,     (%%rax)	\n\t"\
		"movaps	%%xmm5,0x010(%%rax)	\n\t"\
		"subpd	%%xmm7,%%xmm2		\n\t"\
		"subpd	%%xmm6,%%xmm3		\n\t"\
		"movaps	%%xmm2,     (%%rcx)	\n\t"\
		"movaps	%%xmm3,0x010(%%rdx)	\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	%%xmm7,     (%%rdx)	\n\t"\
		"movaps	%%xmm6,0x010(%%rcx)	\n\t"\
	/*** Now do four pass-2 4-DFTs, inputs from local store, outputs back to main-array: ***/\
		"movq	%[__add0],%%rax		\n\t"\
		"movslq	%[__p1],%%rbx		\n\t"\
		"movslq	%[__p2],%%rcx		\n\t"\
		"movslq	%[__p3],%%rdx		\n\t"\
		"movslq	%[__p4],%%rdi		\n\t"\
		"shlq	$3,%%rbx			\n\t"\
		"shlq	$3,%%rcx			\n\t"\
		"shlq	$3,%%rdx			\n\t"\
		"shlq	$3,%%rdi			\n\t"\
		"addq	%%rax,%%rbx			\n\t"\
		"addq	%%rax,%%rcx			\n\t"\
		"addq	%%rax,%%rdx			\n\t"\
		"movq	%[__r1],%%rsi		\n\t"\
		"movaps	     (%%rsi),%%xmm0	\n\t"\
		"movaps	0x010(%%rsi),%%xmm1	\n\t"\
		"movaps	0x040(%%rsi),%%xmm2	\n\t"\
		"movaps	0x050(%%rsi),%%xmm3	\n\t"\
		"subpd	0x040(%%rsi),%%xmm0	\n\t"\
		"subpd	0x050(%%rsi),%%xmm1	\n\t"\
		"addpd	     (%%rsi),%%xmm2	\n\t"\
		"addpd	0x010(%%rsi),%%xmm3	\n\t"\
		"movaps	0x020(%%rsi),%%xmm4	\n\t"\
		"movaps	0x030(%%rsi),%%xmm5	\n\t"\
		"movaps	0x060(%%rsi),%%xmm6	\n\t"\
		"movaps	0x070(%%rsi),%%xmm7	\n\t"\
		"subpd	0x060(%%rsi),%%xmm4	\n\t"\
		"subpd	0x070(%%rsi),%%xmm5	\n\t"\
		"addpd	0x020(%%rsi),%%xmm6	\n\t"\
		"addpd	0x030(%%rsi),%%xmm7	\n\t"\
		"subpd	%%xmm6,%%xmm2		\n\t"\
		"subpd	%%xmm7,%%xmm3		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"movaps	%%xmm2,    (%%rbx)	\n\t"\
		"movaps	%%xmm3,0x10(%%rbx)	\n\t"\
		"addpd	%%xmm2,		%%xmm6	\n\t"\
		"addpd	%%xmm3,		%%xmm7	\n\t"\
		"movaps	%%xmm6,    (%%rax)	\n\t"\
		"movaps	%%xmm7,0x10(%%rax)	\n\t"\
		"subpd	%%xmm5,		%%xmm0	\n\t"\
		"subpd	%%xmm4,		%%xmm1	\n\t"\
		"addpd	%%xmm5,		%%xmm5	\n\t"\
		"addpd	%%xmm4,		%%xmm4	\n\t"\
		"movaps	%%xmm0,    (%%rcx)	\n\t"\
		"movaps	%%xmm1,0x10(%%rdx)	\n\t"\
		"addpd	%%xmm0,		%%xmm5	\n\t"\
		"addpd	%%xmm1,		%%xmm4	\n\t"\
		"movaps	%%xmm5,    (%%rdx)	\n\t"\
		"movaps	%%xmm4,0x10(%%rcx)	\n\t"\
		"addq	%%rdi,%%rax			\n\t"\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	%%rdi,%%rcx			\n\t"\
		"addq	%%rdi,%%rdx			\n\t"\
		"movq	%[__r17],%%rsi		\n\t"\
		"movaps	     (%%rsi),%%xmm0	\n\t"\
		"movaps	0x010(%%rsi),%%xmm1	\n\t"\
		"movaps	0x040(%%rsi),%%xmm2	\n\t"\
		"movaps	0x050(%%rsi),%%xmm3	\n\t"\
		"subpd	0x050(%%rsi),%%xmm0	\n\t"\
		"subpd	0x040(%%rsi),%%xmm1	\n\t"\
		"addpd	0x010(%%rsi),%%xmm2	\n\t"\
		"addpd	     (%%rsi),%%xmm3	\n\t"\
		"movaps	0x020(%%rsi),%%xmm4	\n\t"\
		"movaps	0x030(%%rsi),%%xmm5	\n\t"\
		"movaps	0x060(%%rsi),%%xmm6	\n\t"\
		"movaps	0x070(%%rsi),%%xmm7	\n\t"\
		"subpd	0x030(%%rsi),%%xmm4	\n\t"\
		"addpd	0x020(%%rsi),%%xmm5	\n\t"\
		"addpd	0x070(%%rsi),%%xmm6	\n\t"\
		"subpd	0x060(%%rsi),%%xmm7	\n\t"\
		"movq	%[__isrt2],%%rsi	\n\t"\
		"mulpd	(%%rsi),%%xmm4		\n\t"\
		"mulpd	(%%rsi),%%xmm5		\n\t"\
		"mulpd	(%%rsi),%%xmm6		\n\t"\
		"mulpd	(%%rsi),%%xmm7		\n\t"\
		"subpd	%%xmm6,%%xmm4		\n\t"\
		"subpd	%%xmm7,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm6		\n\t"\
		"addpd	%%xmm5,%%xmm7		\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm2		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"movaps	%%xmm0,    (%%rbx)	\n\t"\
		"movaps	%%xmm2,0x10(%%rbx)	\n\t"\
		"addpd	%%xmm0,		%%xmm4	\n\t"\
		"addpd	%%xmm2,		%%xmm5	\n\t"\
		"movaps	%%xmm4,    (%%rax)	\n\t"\
		"movaps	%%xmm5,0x10(%%rax)	\n\t"\
		"subpd	%%xmm7,		%%xmm3	\n\t"\
		"subpd	%%xmm6,		%%xmm1	\n\t"\
		"addpd	%%xmm7,		%%xmm7	\n\t"\
		"addpd	%%xmm6,		%%xmm6	\n\t"\
		"movaps	%%xmm3,    (%%rcx)	\n\t"\
		"movaps	%%xmm1,0x10(%%rdx)	\n\t"\
		"addpd	%%xmm3,		%%xmm7	\n\t"\
		"addpd	%%xmm1,		%%xmm6	\n\t"\
		"movaps	%%xmm7,    (%%rdx)	\n\t"\
		"movaps	%%xmm6,0x10(%%rcx)	\n\t"\
		"addq	%%rdi,%%rax			\n\t"\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	%%rdi,%%rcx			\n\t"\
		"addq	%%rdi,%%rdx			\n\t"\
		"movq	%[__r9] ,%%rsi		\n\t"\
		"movq	%[__cc0],%%rdi		\n\t"\
		"movaps	0x020(%%rsi),%%xmm4	\n\t"\
		"movaps	0x060(%%rsi),%%xmm6	\n\t"\
		"movaps	0x030(%%rsi),%%xmm5	\n\t"\
		"movaps	0x070(%%rsi),%%xmm7	\n\t"\
		"movaps	%%xmm4,%%xmm0		\n\t"\
		"movaps	%%xmm6,%%xmm2		\n\t"\
		"movaps	%%xmm5,%%xmm1		\n\t"\
		"movaps	%%xmm7,%%xmm3		\n\t"\
		"mulpd	     (%%rdi),%%xmm4	\n\t"\
		"mulpd	0x010(%%rdi),%%xmm6	\n\t"\
		"mulpd	0x010(%%rdi),%%xmm1	\n\t"\
		"mulpd	     (%%rdi),%%xmm3	\n\t"\
		"mulpd	     (%%rdi),%%xmm5	\n\t"\
		"mulpd	0x010(%%rdi),%%xmm7	\n\t"\
		"mulpd	0x010(%%rdi),%%xmm0	\n\t"\
		"mulpd	     (%%rdi),%%xmm2	\n\t"\
		"subpd	%%xmm1,%%xmm4		\n\t"\
		"subpd	%%xmm3,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm5		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"subpd	%%xmm6,%%xmm4		\n\t"\
		"subpd	%%xmm7,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm6		\n\t"\
		"addpd	%%xmm5,%%xmm7		\n\t"\
		"movq	%[__isrt2],%%rdi	\n\t"\
		"movaps	0x040(%%rsi),%%xmm2	\n\t"\
		"movaps	0x050(%%rsi),%%xmm3	\n\t"\
		"subpd	0x050(%%rsi),%%xmm2	\n\t"\
		"addpd	0x040(%%rsi),%%xmm3	\n\t"\
		"mulpd	(%%rdi),%%xmm2		\n\t"\
		"mulpd	(%%rdi),%%xmm3		\n\t"\
		"movaps	     (%%rsi),%%xmm0	\n\t"\
		"movaps	0x010(%%rsi),%%xmm1	\n\t"\
		"movslq	%[__p4],%%rdi		\n\t"\
		"shlq	$3,%%rdi			\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"addpd	     (%%rsi),%%xmm2	\n\t"\
		"addpd	0x010(%%rsi),%%xmm3	\n\t"\
		"subpd	%%xmm6,%%xmm2		\n\t"\
		"subpd	%%xmm7,%%xmm3		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"movaps	%%xmm2,    (%%rbx)	\n\t"\
		"movaps	%%xmm3,0x10(%%rbx)	\n\t"\
		"addpd	%%xmm2,		%%xmm6	\n\t"\
		"addpd	%%xmm3,		%%xmm7	\n\t"\
		"movaps	%%xmm6,    (%%rax)	\n\t"\
		"movaps	%%xmm7,0x10(%%rax)	\n\t"\
		"subpd	%%xmm5,		%%xmm0	\n\t"\
		"subpd	%%xmm4,		%%xmm1	\n\t"\
		"addpd	%%xmm5,		%%xmm5	\n\t"\
		"addpd	%%xmm4,		%%xmm4	\n\t"\
		"movaps	%%xmm0,    (%%rcx)	\n\t"\
		"movaps	%%xmm1,0x10(%%rdx)	\n\t"\
		"addpd	%%xmm0,		%%xmm5	\n\t"\
		"addpd	%%xmm1,		%%xmm4	\n\t"\
		"movaps	%%xmm5,    (%%rdx)	\n\t"\
		"movaps	%%xmm4,0x10(%%rcx)	\n\t"\
		"addq	%%rdi,%%rax			\n\t"\
		"addq	%%rdi,%%rbx			\n\t"\
		"addq	%%rdi,%%rcx			\n\t"\
		"addq	%%rdi,%%rdx			\n\t"\
		"movq	%[__r25],%%rsi		\n\t"\
		"movq	%[__cc0],%%rdi		\n\t"\
		"movaps	0x020(%%rsi),%%xmm4	\n\t"\
		"movaps	0x060(%%rsi),%%xmm6	\n\t"\
		"movaps	0x030(%%rsi),%%xmm5	\n\t"\
		"movaps	0x070(%%rsi),%%xmm7	\n\t"\
		"movaps	%%xmm4,%%xmm0		\n\t"\
		"movaps	%%xmm6,%%xmm2		\n\t"\
		"movaps	%%xmm5,%%xmm1		\n\t"\
		"movaps	%%xmm7,%%xmm3		\n\t"\
		"mulpd	0x010(%%rdi),%%xmm4	\n\t"\
		"mulpd	     (%%rdi),%%xmm6	\n\t"\
		"mulpd	     (%%rdi),%%xmm1	\n\t"\
		"mulpd	0x010(%%rdi),%%xmm3	\n\t"\
		"mulpd	0x010(%%rdi),%%xmm5	\n\t"\
		"mulpd	     (%%rdi),%%xmm7	\n\t"\
		"mulpd	     (%%rdi),%%xmm0	\n\t"\
		"mulpd	0x010(%%rdi),%%xmm2	\n\t"\
		"subpd	%%xmm1,%%xmm4		\n\t"\
		"subpd	%%xmm3,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm5		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"subpd	%%xmm6,%%xmm4		\n\t"\
		"subpd	%%xmm7,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm6		\n\t"\
		"addpd	%%xmm5,%%xmm7		\n\t"\
		"movq	%[__isrt2],%%rdi	\n\t"\
		"movaps	0x040(%%rsi),%%xmm2	\n\t"\
		"movaps	0x050(%%rsi),%%xmm3	\n\t"\
		"addpd	0x050(%%rsi),%%xmm2	\n\t"\
		"subpd	0x040(%%rsi),%%xmm3	\n\t"\
		"mulpd	(%%rdi),%%xmm2		\n\t"\
		"mulpd	(%%rdi),%%xmm3		\n\t"\
		"movaps	     (%%rsi),%%xmm0	\n\t"\
		"movaps	0x010(%%rsi),%%xmm1	\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"addpd	     (%%rsi),%%xmm2	\n\t"\
		"addpd	0x010(%%rsi),%%xmm3	\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"movaps	%%xmm0,    (%%rbx)	\n\t"\
		"movaps	%%xmm1,0x10(%%rbx)	\n\t"\
		"addpd	%%xmm0,		%%xmm4	\n\t"\
		"addpd	%%xmm1,		%%xmm5	\n\t"\
		"movaps	%%xmm4,    (%%rax)	\n\t"\
		"movaps	%%xmm5,0x10(%%rax)	\n\t"\
		"subpd	%%xmm7,		%%xmm2	\n\t"\
		"subpd	%%xmm6,		%%xmm3	\n\t"\
		"addpd	%%xmm7,		%%xmm7	\n\t"\
		"addpd	%%xmm6,		%%xmm6	\n\t"\
		"movaps	%%xmm2,    (%%rcx)	\n\t"\
		"movaps	%%xmm3,0x10(%%rdx)	\n\t"\
		"addpd	%%xmm2,		%%xmm7	\n\t"\
		"addpd	%%xmm3,		%%xmm6	\n\t"\
		"movaps	%%xmm7,    (%%rdx)	\n\t"\
		"movaps	%%xmm6,0x10(%%rcx)	\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__p1] "m" (Xp1)\
		 ,[__p2] "m" (Xp2)\
		 ,[__p3] "m" (Xp3)\
		 ,[__p4] "m" (Xp4)\
		 ,[__p8] "m" (Xp8)\
		 ,[__r1] "m" (Xr1)\
		 ,[__r3] "m" (Xr3)\
		 ,[__r5] "m" (Xr5)\
		 ,[__r7] "m" (Xr7)\
		 ,[__r9] "m" (Xr9)\
		 ,[__r11] "m" (Xr11)\
		 ,[__r13] "m" (Xr13)\
		 ,[__r15] "m" (Xr15)\
		 ,[__r17] "m" (Xr17)\
		 ,[__r19] "m" (Xr19)\
		 ,[__r21] "m" (Xr21)\
		 ,[__r23] "m" (Xr23)\
		 ,[__r25] "m" (Xr25)\
		 ,[__r27] "m" (Xr27)\
		 ,[__r29] "m" (Xr29)\
		 ,[__r31] "m" (Xr31)\
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__cc0] "m" (Xcc0)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
	);\
	}

#elif OS_BITS == 32	// 32-bit SSE2

	#define SSE2_RADIX16_DIT_NOTWIDDLE(Xadd0,Xp1,Xp2,Xp3,Xp4,Xp8,Xr1,Xr3,Xr5,Xr7,Xr9,Xr11,Xr17,Xr25,Xisrt2,Xcc0)\
	{\
	__asm__ volatile (\
	"pushl %%ebx	\n\t"/* Explicit save/restore of PIC register */\
		"movl	%[__add0],%%eax		\n\t"\
		"movl	%[__p1],%%ebx		\n\t"\
		"movl	%[__p2],%%ecx		\n\t"\
		"movl	%[__p3],%%edx		\n\t"\
		"movl	%[__p4],%%edi		\n\t"\
		"shll	$3,%%ebx			\n\t"\
		"shll	$3,%%ecx			\n\t"\
		"shll	$3,%%edx			\n\t"\
		"shll	$3,%%edi			\n\t"\
		"addl	%%eax,%%ebx			\n\t"\
		"addl	%%eax,%%ecx			\n\t"\
		"addl	%%eax,%%edx			\n\t"\
		"/* SSE2_RADIX4_DIT_0TWIDDLE_B(r1 ): */\n\t"\
		"movaps	    (%%eax),%%xmm2	\n\t"\
		"movaps	    (%%ecx),%%xmm6	\n\t"\
		"movaps	0x10(%%eax),%%xmm3	\n\t"\
		"movaps	0x10(%%ecx),%%xmm7	\n\t"\
		"movaps	    (%%ebx),%%xmm0	\n\t"\
		"movaps	    (%%edx),%%xmm4	\n\t"\
		"movaps	0x10(%%ebx),%%xmm1	\n\t"\
		"movaps	0x10(%%edx),%%xmm5	\n\t"\
		"subpd	%%xmm0,%%xmm2		\n\t"\
		"subpd	%%xmm4,%%xmm6		\n\t"\
		"subpd	%%xmm1,%%xmm3		\n\t"\
		"subpd	%%xmm5,%%xmm7		\n\t"\
		"addpd	%%xmm0,%%xmm0		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm1,%%xmm1		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm2,%%xmm0		\n\t"\
		"addpd	%%xmm6,%%xmm4		\n\t"\
		"addpd	%%xmm3,%%xmm1		\n\t"\
		"addpd	%%xmm7,%%xmm5		\n\t"\
		"movl	%[__r1], %%esi		\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm2		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"subpd	%%xmm6,%%xmm3		\n\t"\
		"movaps	%%xmm0,0x040(%%esi)	\n\t"\
		"movaps	%%xmm2,0x060(%%esi)	\n\t"\
		"movaps	%%xmm1,0x050(%%esi)	\n\t"\
		"movaps	%%xmm3,0x030(%%esi)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"addpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	%%xmm4,     (%%esi)	\n\t"\
		"movaps	%%xmm7,0x020(%%esi)	\n\t"\
		"movaps	%%xmm5,0x010(%%esi)	\n\t"\
		"movaps	%%xmm6,0x070(%%esi)	\n\t"\
		"addl	%%edi,%%eax			\n\t"\
		"addl	%%edi,%%ebx			\n\t"\
		"addl	%%edi,%%ecx			\n\t"\
		"addl	%%edi,%%edx			\n\t"\
		"/* SSE2_RADIX4_DIT_0TWIDDLE_B(r9 ): */\n\t"\
		"movaps	    (%%eax),%%xmm2	\n\t"\
		"movaps	    (%%ecx),%%xmm6	\n\t"\
		"movaps	0x10(%%eax),%%xmm3	\n\t"\
		"movaps	0x10(%%ecx),%%xmm7	\n\t"\
		"movaps	    (%%ebx),%%xmm0	\n\t"\
		"movaps	    (%%edx),%%xmm4	\n\t"\
		"movaps	0x10(%%ebx),%%xmm1	\n\t"\
		"movaps	0x10(%%edx),%%xmm5	\n\t"\
		"subpd	%%xmm0,%%xmm2		\n\t"\
		"subpd	%%xmm4,%%xmm6		\n\t"\
		"subpd	%%xmm1,%%xmm3		\n\t"\
		"subpd	%%xmm5,%%xmm7		\n\t"\
		"addpd	%%xmm0,%%xmm0		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm1,%%xmm1		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm2,%%xmm0		\n\t"\
		"addpd	%%xmm6,%%xmm4		\n\t"\
		"addpd	%%xmm3,%%xmm1		\n\t"\
		"addpd	%%xmm7,%%xmm5		\n\t"\
		"movl	%[__r9], %%esi		\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm2		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"subpd	%%xmm6,%%xmm3		\n\t"\
		"movaps	%%xmm0,0x040(%%esi)	\n\t"\
		"movaps	%%xmm2,0x060(%%esi)	\n\t"\
		"movaps	%%xmm1,0x050(%%esi)	\n\t"\
		"movaps	%%xmm3,0x030(%%esi)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"addpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	%%xmm4,     (%%esi)	\n\t"\
		"movaps	%%xmm7,0x020(%%esi)	\n\t"\
		"movaps	%%xmm5,0x010(%%esi)	\n\t"\
		"movaps	%%xmm6,0x070(%%esi)	\n\t"\
		"addl	%%edi,%%eax			\n\t"\
		"addl	%%edi,%%ebx			\n\t"\
		"addl	%%edi,%%ecx			\n\t"\
		"addl	%%edi,%%edx			\n\t"\
		"/* SSE2_RADIX4_DIT_0TWIDDLE_B(r17): */\n\t"\
		"movaps	    (%%eax),%%xmm2	\n\t"\
		"movaps	    (%%ecx),%%xmm6	\n\t"\
		"movaps	0x10(%%eax),%%xmm3	\n\t"\
		"movaps	0x10(%%ecx),%%xmm7	\n\t"\
		"movaps	    (%%ebx),%%xmm0	\n\t"\
		"movaps	    (%%edx),%%xmm4	\n\t"\
		"movaps	0x10(%%ebx),%%xmm1	\n\t"\
		"movaps	0x10(%%edx),%%xmm5	\n\t"\
		"subpd	%%xmm0,%%xmm2		\n\t"\
		"subpd	%%xmm4,%%xmm6		\n\t"\
		"subpd	%%xmm1,%%xmm3		\n\t"\
		"subpd	%%xmm5,%%xmm7		\n\t"\
		"addpd	%%xmm0,%%xmm0		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm1,%%xmm1		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm2,%%xmm0		\n\t"\
		"addpd	%%xmm6,%%xmm4		\n\t"\
		"addpd	%%xmm3,%%xmm1		\n\t"\
		"addpd	%%xmm7,%%xmm5		\n\t"\
		"movl	%[__r17], %%esi		\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm2		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"subpd	%%xmm6,%%xmm3		\n\t"\
		"movaps	%%xmm0,0x040(%%esi)	\n\t"\
		"movaps	%%xmm2,0x060(%%esi)	\n\t"\
		"movaps	%%xmm1,0x050(%%esi)	\n\t"\
		"movaps	%%xmm3,0x030(%%esi)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"addpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	%%xmm4,     (%%esi)	\n\t"\
		"movaps	%%xmm7,0x020(%%esi)	\n\t"\
		"movaps	%%xmm5,0x010(%%esi)	\n\t"\
		"movaps	%%xmm6,0x070(%%esi)	\n\t"\
		"addl	%%edi,%%eax			\n\t"\
		"addl	%%edi,%%ebx			\n\t"\
		"addl	%%edi,%%ecx			\n\t"\
		"addl	%%edi,%%edx			\n\t"\
		"/* SSE2_RADIX4_DIT_0TWIDDLE_B(r25): */\n\t"\
		"movaps	    (%%eax),%%xmm2	\n\t"\
		"movaps	    (%%ecx),%%xmm6	\n\t"\
		"movaps	0x10(%%eax),%%xmm3	\n\t"\
		"movaps	0x10(%%ecx),%%xmm7	\n\t"\
		"movaps	    (%%ebx),%%xmm0	\n\t"\
		"movaps	    (%%edx),%%xmm4	\n\t"\
		"movaps	0x10(%%ebx),%%xmm1	\n\t"\
		"movaps	0x10(%%edx),%%xmm5	\n\t"\
		"subpd	%%xmm0,%%xmm2		\n\t"\
		"subpd	%%xmm4,%%xmm6		\n\t"\
		"subpd	%%xmm1,%%xmm3		\n\t"\
		"subpd	%%xmm5,%%xmm7		\n\t"\
		"addpd	%%xmm0,%%xmm0		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm1,%%xmm1		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm2,%%xmm0		\n\t"\
		"addpd	%%xmm6,%%xmm4		\n\t"\
		"addpd	%%xmm3,%%xmm1		\n\t"\
		"addpd	%%xmm7,%%xmm5		\n\t"\
		"movl	%[__r25], %%esi		\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm2		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"subpd	%%xmm6,%%xmm3		\n\t"\
		"movaps	%%xmm0,0x040(%%esi)	\n\t"\
		"movaps	%%xmm2,0x060(%%esi)	\n\t"\
		"movaps	%%xmm1,0x050(%%esi)	\n\t"\
		"movaps	%%xmm3,0x030(%%esi)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"addpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	%%xmm4,     (%%esi)	\n\t"\
		"movaps	%%xmm7,0x020(%%esi)	\n\t"\
		"movaps	%%xmm5,0x010(%%esi)	\n\t"\
		"movaps	%%xmm6,0x070(%%esi)	\n\t"\
		"movl	%[__r1],%%eax		\n\t"\
		"movaps	     (%%eax),%%xmm0	\n\t"\
		"movaps	0x010(%%eax),%%xmm1	\n\t"\
		"movaps	0x080(%%eax),%%xmm2	\n\t"\
		"movaps	0x090(%%eax),%%xmm3	\n\t"\
		"subpd	0x080(%%eax),%%xmm0	\n\t"\
		"subpd	0x090(%%eax),%%xmm1	\n\t"\
		"addpd	     (%%eax),%%xmm2	\n\t"\
		"addpd	0x010(%%eax),%%xmm3	\n\t"\
		"movaps	0x100(%%eax),%%xmm4	\n\t"\
		"movaps	0x110(%%eax),%%xmm5	\n\t"\
		"movaps	0x180(%%eax),%%xmm6	\n\t"\
		"movaps	0x190(%%eax),%%xmm7	\n\t"\
		"subpd	0x180(%%eax),%%xmm4	\n\t"\
		"subpd	0x190(%%eax),%%xmm5	\n\t"\
		"addpd	0x100(%%eax),%%xmm6	\n\t"\
		"addpd	0x110(%%eax),%%xmm7	\n\t"\
		"subpd	%%xmm6,%%xmm2		\n\t"\
		"subpd	%%xmm7,%%xmm3		\n\t"\
		"movaps	%%xmm2,0x100(%%eax)	\n\t"\
		"movaps	%%xmm3,0x110(%%eax)	\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm2,%%xmm6		\n\t"\
		"addpd	%%xmm3,%%xmm7		\n\t"\
		"movaps	%%xmm6,     (%%eax)	\n\t"\
		"movaps	%%xmm7,0x010(%%eax)	\n\t"\
		"subpd	%%xmm5,%%xmm0		\n\t"\
		"subpd	%%xmm4,%%xmm1		\n\t"\
		"movaps	%%xmm0,0x180(%%eax)	\n\t"\
		"movaps	%%xmm1,0x090(%%eax)	\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm0,%%xmm5		\n\t"\
		"addpd	%%xmm1,%%xmm4		\n\t"\
		"movaps	%%xmm5,0x080(%%eax)	\n\t"\
		"movaps	%%xmm4,0x190(%%eax)	\n\t"\
		"movl	%[__r5],%%eax		\n\t"\
		"movl	%[__isrt2],%%ebx	\n\t"\
		"movaps	(%%ebx),%%xmm2		\n\t"\
		"movaps	0x100(%%eax),%%xmm4	\n\t"\
		"movaps	0x110(%%eax),%%xmm5	\n\t"\
		"movaps	0x180(%%eax),%%xmm0	\n\t"\
		"movaps	0x190(%%eax),%%xmm1	\n\t"\
		"addpd	0x110(%%eax),%%xmm4	\n\t"\
		"subpd	0x100(%%eax),%%xmm5	\n\t"\
		"subpd	0x190(%%eax),%%xmm0	\n\t"\
		"addpd	0x180(%%eax),%%xmm1	\n\t"\
		"mulpd	%%xmm2,%%xmm4		\n\t"\
		"mulpd	%%xmm2,%%xmm5		\n\t"\
		"mulpd	%%xmm2,%%xmm0		\n\t"\
		"mulpd	%%xmm2,%%xmm1		\n\t"\
		"movaps	%%xmm4,%%xmm6		\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t"\
		"subpd	%%xmm0,%%xmm4		\n\t"\
		"subpd	%%xmm1,%%xmm5		\n\t"\
		"addpd	%%xmm0,%%xmm6		\n\t"\
		"addpd	%%xmm1,%%xmm7		\n\t"\
		"movaps	     (%%eax),%%xmm0	\n\t"\
		"movaps	0x010(%%eax),%%xmm1	\n\t"\
		"movaps	0x080(%%eax),%%xmm2	\n\t"\
		"movaps	0x090(%%eax),%%xmm3	\n\t"\
		"subpd	0x090(%%eax),%%xmm0	\n\t"\
		"subpd	0x080(%%eax),%%xmm1	\n\t"\
		"addpd	     (%%eax),%%xmm3	\n\t"\
		"addpd	0x010(%%eax),%%xmm2	\n\t"\
		"subpd	%%xmm4,%%xmm3		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"movaps	%%xmm3,0x100(%%eax)	\n\t"\
		"movaps	%%xmm1,0x110(%%eax)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm3,%%xmm4		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"movaps	%%xmm4,     (%%eax)	\n\t"\
		"movaps	%%xmm5,0x010(%%eax)	\n\t"\
		"subpd	%%xmm7,%%xmm0		\n\t"\
		"subpd	%%xmm6,%%xmm2		\n\t"\
		"movaps	%%xmm0,0x180(%%eax)	\n\t"\
		"movaps	%%xmm2,0x090(%%eax)	\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm7		\n\t"\
		"addpd	%%xmm2,%%xmm6		\n\t"\
		"movaps	%%xmm7,0x080(%%eax)	\n\t"\
		"movaps	%%xmm6,0x190(%%eax)	\n\t"\
		"movl	%[__r3],%%eax		\n\t"\
		"movl	%[__isrt2],%%ebx	\n\t"\
		"movl	%[__cc0],%%ecx		\n\t"\
		"movaps	0x180(%%eax),%%xmm0	\n\t"\
		"movaps	0x190(%%eax),%%xmm1	\n\t"\
		"movaps	0x180(%%eax),%%xmm2	\n\t"\
		"movaps	0x190(%%eax),%%xmm3	\n\t"\
		"mulpd	0x010(%%ecx),%%xmm0	\n\t"\
		"mulpd	0x010(%%ecx),%%xmm1	\n\t"\
		"mulpd	     (%%ecx),%%xmm2	\n\t"\
		"mulpd	     (%%ecx),%%xmm3	\n\t"\
		"subpd	%%xmm2,%%xmm1		\n\t"\
		"addpd	%%xmm3,%%xmm0		\n\t"\
		"movaps	0x100(%%eax),%%xmm4	\n\t"\
		"movaps	0x110(%%eax),%%xmm5	\n\t"\
		"movaps	0x100(%%eax),%%xmm6	\n\t"\
		"movaps	0x110(%%eax),%%xmm7	\n\t"\
		"mulpd	     (%%ecx),%%xmm4	\n\t"\
		"mulpd	     (%%ecx),%%xmm5	\n\t"\
		"mulpd	0x010(%%ecx),%%xmm6	\n\t"\
		"mulpd	0x010(%%ecx),%%xmm7	\n\t"\
		"subpd	%%xmm6,%%xmm5		\n\t"\
		"addpd	%%xmm7,%%xmm4		\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t"\
		"movaps	%%xmm4,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"subpd	%%xmm0,%%xmm6		\n\t"\
		"subpd	%%xmm1,%%xmm7		\n\t"\
		"movaps	0x080(%%eax),%%xmm2	\n\t"\
		"movaps	0x090(%%eax),%%xmm3	\n\t"\
		"movaps	     (%%eax),%%xmm0	\n\t"\
		"movaps	0x010(%%eax),%%xmm1	\n\t"\
		"addpd	0x090(%%eax),%%xmm2	\n\t"\
		"subpd	0x080(%%eax),%%xmm3	\n\t"\
		"mulpd	     (%%ebx),%%xmm2	\n\t"\
		"mulpd	     (%%ebx),%%xmm3	\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"subpd	%%xmm4,%%xmm2		\n\t"\
		"subpd	%%xmm5,%%xmm3		\n\t"\
		"movaps	%%xmm2,0x100(%%eax)	\n\t"\
		"movaps	%%xmm3,0x110(%%eax)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm2,%%xmm4		\n\t"\
		"addpd	%%xmm3,%%xmm5		\n\t"\
		"movaps	%%xmm4,     (%%eax)	\n\t"\
		"movaps	%%xmm5,0x010(%%eax)	\n\t"\
		"subpd	%%xmm7,%%xmm0		\n\t"\
		"subpd	%%xmm6,%%xmm1		\n\t"\
		"movaps	%%xmm0,0x180(%%eax)	\n\t"\
		"movaps	%%xmm1,0x090(%%eax)	\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm7		\n\t"\
		"addpd	%%xmm1,%%xmm6		\n\t"\
		"movaps	%%xmm7,0x080(%%eax)	\n\t"\
		"movaps	%%xmm6,0x190(%%eax)	\n\t"\
		"movl	%[__r7],%%eax		\n\t"\
		"movl	%[__isrt2],%%ebx	\n\t"\
		"movl	%[__cc0],%%ecx		\n\t"\
		"movaps	0x180(%%eax),%%xmm0	\n\t"\
		"movaps	0x190(%%eax),%%xmm1	\n\t"\
		"movaps	0x180(%%eax),%%xmm2	\n\t"\
		"movaps	0x190(%%eax),%%xmm3	\n\t"\
		"mulpd	     (%%ecx),%%xmm0	\n\t"\
		"mulpd	     (%%ecx),%%xmm1	\n\t"\
		"mulpd	0x010(%%ecx),%%xmm2	\n\t"\
		"mulpd	0x010(%%ecx),%%xmm3	\n\t"\
		"subpd	%%xmm2,%%xmm1		\n\t"\
		"addpd	%%xmm3,%%xmm0		\n\t"\
		"movaps	0x100(%%eax),%%xmm4	\n\t"\
		"movaps	0x110(%%eax),%%xmm5	\n\t"\
		"movaps	0x100(%%eax),%%xmm6	\n\t"\
		"movaps	0x110(%%eax),%%xmm7	\n\t"\
		"mulpd	0x010(%%ecx),%%xmm4	\n\t"\
		"mulpd	0x010(%%ecx),%%xmm5	\n\t"\
		"mulpd	     (%%ecx),%%xmm6	\n\t"\
		"mulpd	     (%%ecx),%%xmm7	\n\t"\
		"subpd	%%xmm6,%%xmm5		\n\t"\
		"addpd	%%xmm7,%%xmm4		\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t"\
		"movaps	%%xmm4,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"subpd	%%xmm0,%%xmm6		\n\t"\
		"subpd	%%xmm1,%%xmm7		\n\t"\
		"movaps	0x080(%%eax),%%xmm2	\n\t"\
		"movaps	0x090(%%eax),%%xmm3	\n\t"\
		"movaps	     (%%eax),%%xmm0	\n\t"\
		"movaps	0x010(%%eax),%%xmm1	\n\t"\
		"subpd	0x090(%%eax),%%xmm2	\n\t"\
		"addpd	0x080(%%eax),%%xmm3	\n\t"\
		"mulpd	     (%%ebx),%%xmm2	\n\t"\
		"mulpd	     (%%ebx),%%xmm3	\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"subpd	%%xmm6,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm1		\n\t"\
		"movaps	%%xmm0,0x100(%%eax)	\n\t"\
		"movaps	%%xmm1,0x110(%%eax)	\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm0,%%xmm6		\n\t"\
		"addpd	%%xmm1,%%xmm7		\n\t"\
		"movaps	%%xmm6,     (%%eax)	\n\t"\
		"movaps	%%xmm7,0x010(%%eax)	\n\t"\
		"subpd	%%xmm5,%%xmm2		\n\t"\
		"subpd	%%xmm4,%%xmm3		\n\t"\
		"movaps	%%xmm2,0x180(%%eax)	\n\t"\
		"movaps	%%xmm3,0x090(%%eax)	\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm2,%%xmm5		\n\t"\
		"addpd	%%xmm3,%%xmm4		\n\t"\
		"movaps	%%xmm5,0x080(%%eax)	\n\t"\
		"movaps	%%xmm4,0x190(%%eax)	\n\t"\
	"popl %%ebx	\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__p1] "m" (Xp1)\
		 ,[__p2] "m" (Xp2)\
		 ,[__p3] "m" (Xp3)\
		 ,[__p4] "m" (Xp4)\
		 ,[__p8] "m" (Xp8)\
		 ,[__r1] "m" (Xr1)\
		 ,[__r3] "m" (Xr3)\
		 ,[__r5] "m" (Xr5)\
		 ,[__r7] "m" (Xr7)\
		 ,[__r9] "m" (Xr9)\
		 ,[__r17] "m" (Xr17)\
		 ,[__r25] "m" (Xr25)\
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__cc0] "m" (Xcc0)\
		: "cc","memory","eax",/*"ebx",*/"ecx","edx","edi","esi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"	/* Clobbered registers */\
	);\
	}

	#define SSE2_RADIX16_DIF_NOTWIDDLE(Xadd0,Xp1,Xp2,Xp3,Xp4,Xp8,Xr1,Xr3,Xr5,Xr7,Xr9,Xr11,Xr13,Xr15,Xr17,Xr19,Xr21,Xr23,Xr25,Xr27,Xr29,Xr31,Xisrt2,Xcc0)\
	{\
	__asm__ volatile (\
	"pushl %%ebx	\n\t"/* Explicit save/restore of PIC register */\
		"/* SSE2_RADIX4_DIF_IN_PLACE(r1 , r17, r9 , r25): */\n\t"\
		"movl	%[__r1] ,%%eax		\n\t"\
		"movl	%[__r17],%%ebx		\n\t"\
		"movl	%[__r9] ,%%ecx		\n\t"\
		"movl	%[__r25],%%edx		\n\t"\
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"movaps	    (%%eax),%%xmm2	\n\t"\
		"movaps	0x10(%%eax),%%xmm3	\n\t"\
		"addpd	    (%%ebx),%%xmm0	\n\t"\
		"addpd	0x10(%%ebx),%%xmm1	\n\t"\
		"subpd	    (%%ebx),%%xmm2	\n\t"\
		"subpd	0x10(%%ebx),%%xmm3	\n\t"\
		"movaps	    (%%ecx),%%xmm4	\n\t"\
		"movaps	0x10(%%ecx),%%xmm5	\n\t"\
		"movaps	    (%%ecx),%%xmm6	\n\t"\
		"movaps	0x10(%%ecx),%%xmm7	\n\t"\
		"addpd	    (%%edx),%%xmm4	\n\t"\
		"addpd	0x10(%%edx),%%xmm5	\n\t"\
		"subpd	    (%%edx),%%xmm6	\n\t"\
		"subpd	0x10(%%edx),%%xmm7	\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"movaps	%%xmm0,     (%%ebx)	\n\t"\
		"movaps	%%xmm1,0x010(%%ebx)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"movaps	%%xmm4,     (%%eax)	\n\t"\
		"movaps	%%xmm5,0x010(%%eax)	\n\t"\
		"subpd	%%xmm7,%%xmm2		\n\t"\
		"subpd	%%xmm6,%%xmm3		\n\t"\
		"movaps	%%xmm2,     (%%ecx)	\n\t"\
		"movaps	%%xmm3,0x010(%%edx)	\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	%%xmm7,     (%%edx)	\n\t"\
		"movaps	%%xmm6,0x010(%%ecx)	\n\t"\
		"/* SSE2_RADIX4_DIF_IN_PLACE(r5 , r21, r13, r29): */\n\t"\
		"movl	%[__r5] ,%%eax		\n\t"\
		"movl	%[__r21],%%ebx		\n\t"\
		"movl	%[__r13],%%ecx		\n\t"\
		"movl	%[__r29],%%edx		\n\t"\
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"movaps	    (%%eax),%%xmm2	\n\t"\
		"movaps	0x10(%%eax),%%xmm3	\n\t"\
		"addpd	    (%%ebx),%%xmm0	\n\t"\
		"addpd	0x10(%%ebx),%%xmm1	\n\t"\
		"subpd	    (%%ebx),%%xmm2	\n\t"\
		"subpd	0x10(%%ebx),%%xmm3	\n\t"\
		"movaps	    (%%ecx),%%xmm4	\n\t"\
		"movaps	0x10(%%ecx),%%xmm5	\n\t"\
		"movaps	    (%%ecx),%%xmm6	\n\t"\
		"movaps	0x10(%%ecx),%%xmm7	\n\t"\
		"addpd	    (%%edx),%%xmm4	\n\t"\
		"addpd	0x10(%%edx),%%xmm5	\n\t"\
		"subpd	    (%%edx),%%xmm6	\n\t"\
		"subpd	0x10(%%edx),%%xmm7	\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"movaps	%%xmm0,     (%%ebx)	\n\t"\
		"movaps	%%xmm1,0x010(%%ebx)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"movaps	%%xmm4,     (%%eax)	\n\t"\
		"movaps	%%xmm5,0x010(%%eax)	\n\t"\
		"subpd	%%xmm7,%%xmm2		\n\t"\
		"subpd	%%xmm6,%%xmm3		\n\t"\
		"movaps	%%xmm2,     (%%ecx)	\n\t"\
		"movaps	%%xmm3,0x010(%%edx)	\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	%%xmm7,     (%%edx)	\n\t"\
		"movaps	%%xmm6,0x010(%%ecx)	\n\t"\
		"/* SSE2_RADIX4_DIF_IN_PLACE(r3 , r19, r11, r27): */\n\t"\
		"movl	%[__r3] ,%%eax		\n\t"\
		"movl	%[__r19],%%ebx		\n\t"\
		"movl	%[__r11],%%ecx		\n\t"\
		"movl	%[__r27],%%edx		\n\t"\
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"movaps	    (%%eax),%%xmm2	\n\t"\
		"movaps	0x10(%%eax),%%xmm3	\n\t"\
		"addpd	    (%%ebx),%%xmm0	\n\t"\
		"addpd	0x10(%%ebx),%%xmm1	\n\t"\
		"subpd	    (%%ebx),%%xmm2	\n\t"\
		"subpd	0x10(%%ebx),%%xmm3	\n\t"\
		"movaps	    (%%ecx),%%xmm4	\n\t"\
		"movaps	0x10(%%ecx),%%xmm5	\n\t"\
		"movaps	    (%%ecx),%%xmm6	\n\t"\
		"movaps	0x10(%%ecx),%%xmm7	\n\t"\
		"addpd	    (%%edx),%%xmm4	\n\t"\
		"addpd	0x10(%%edx),%%xmm5	\n\t"\
		"subpd	    (%%edx),%%xmm6	\n\t"\
		"subpd	0x10(%%edx),%%xmm7	\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"movaps	%%xmm0,     (%%ebx)	\n\t"\
		"movaps	%%xmm1,0x010(%%ebx)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"movaps	%%xmm4,     (%%eax)	\n\t"\
		"movaps	%%xmm5,0x010(%%eax)	\n\t"\
		"subpd	%%xmm7,%%xmm2		\n\t"\
		"subpd	%%xmm6,%%xmm3		\n\t"\
		"movaps	%%xmm2,     (%%ecx)	\n\t"\
		"movaps	%%xmm3,0x010(%%edx)	\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	%%xmm7,     (%%edx)	\n\t"\
		"movaps	%%xmm6,0x010(%%ecx)	\n\t"\
		"/* SSE2_RADIX4_DIF_IN_PLACE(r7 , r23, r15, r31): */\n\t"\
		"movl	%[__r7] ,%%eax		\n\t"\
		"movl	%[__r23],%%ebx		\n\t"\
		"movl	%[__r15],%%ecx		\n\t"\
		"movl	%[__r31],%%edx		\n\t"\
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"movaps	    (%%eax),%%xmm2	\n\t"\
		"movaps	0x10(%%eax),%%xmm3	\n\t"\
		"addpd	    (%%ebx),%%xmm0	\n\t"\
		"addpd	0x10(%%ebx),%%xmm1	\n\t"\
		"subpd	    (%%ebx),%%xmm2	\n\t"\
		"subpd	0x10(%%ebx),%%xmm3	\n\t"\
		"movaps	    (%%ecx),%%xmm4	\n\t"\
		"movaps	0x10(%%ecx),%%xmm5	\n\t"\
		"movaps	    (%%ecx),%%xmm6	\n\t"\
		"movaps	0x10(%%ecx),%%xmm7	\n\t"\
		"addpd	    (%%edx),%%xmm4	\n\t"\
		"addpd	0x10(%%edx),%%xmm5	\n\t"\
		"subpd	    (%%edx),%%xmm6	\n\t"\
		"subpd	0x10(%%edx),%%xmm7	\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"movaps	%%xmm0,     (%%ebx)	\n\t"\
		"movaps	%%xmm1,0x010(%%ebx)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"movaps	%%xmm4,     (%%eax)	\n\t"\
		"movaps	%%xmm5,0x010(%%eax)	\n\t"\
		"subpd	%%xmm7,%%xmm2		\n\t"\
		"subpd	%%xmm6,%%xmm3		\n\t"\
		"movaps	%%xmm2,     (%%ecx)	\n\t"\
		"movaps	%%xmm3,0x010(%%edx)	\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	%%xmm7,     (%%edx)	\n\t"\
		"movaps	%%xmm6,0x010(%%ecx)	\n\t"\
		"movl	%[__add0],%%eax		\n\t"\
		"movl	%[__p1],%%ebx		\n\t"\
		"movl	%[__p2],%%ecx		\n\t"\
		"movl	%[__p3],%%edx		\n\t"\
		"movl	%[__p4],%%edi		\n\t"\
		"shll	$3,%%ebx			\n\t"\
		"shll	$3,%%ecx			\n\t"\
		"shll	$3,%%edx			\n\t"\
		"shll	$3,%%edi			\n\t"\
		"addl	%%eax,%%ebx			\n\t"\
		"addl	%%eax,%%ecx			\n\t"\
		"addl	%%eax,%%edx			\n\t"\
		"movl	%[__r1],%%esi		\n\t"\
		"movaps	     (%%esi),%%xmm0	\n\t"\
		"movaps	0x010(%%esi),%%xmm1	\n\t"\
		"movaps	0x040(%%esi),%%xmm2	\n\t"\
		"movaps	0x050(%%esi),%%xmm3	\n\t"\
		"subpd	0x040(%%esi),%%xmm0	\n\t"\
		"subpd	0x050(%%esi),%%xmm1	\n\t"\
		"addpd	     (%%esi),%%xmm2	\n\t"\
		"addpd	0x010(%%esi),%%xmm3	\n\t"\
		"movaps	0x020(%%esi),%%xmm4	\n\t"\
		"movaps	0x030(%%esi),%%xmm5	\n\t"\
		"movaps	0x060(%%esi),%%xmm6	\n\t"\
		"movaps	0x070(%%esi),%%xmm7	\n\t"\
		"subpd	0x060(%%esi),%%xmm4	\n\t"\
		"subpd	0x070(%%esi),%%xmm5	\n\t"\
		"addpd	0x020(%%esi),%%xmm6	\n\t"\
		"addpd	0x030(%%esi),%%xmm7	\n\t"\
		"subpd	%%xmm6,%%xmm2		\n\t"\
		"subpd	%%xmm7,%%xmm3		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"movaps	%%xmm2,    (%%ebx)	\n\t"\
		"movaps	%%xmm3,0x10(%%ebx)	\n\t"\
		"addpd	%%xmm2,		%%xmm6	\n\t"\
		"addpd	%%xmm3,		%%xmm7	\n\t"\
		"movaps	%%xmm6,    (%%eax)	\n\t"\
		"movaps	%%xmm7,0x10(%%eax)	\n\t"\
		"subpd	%%xmm5,		%%xmm0	\n\t"\
		"subpd	%%xmm4,		%%xmm1	\n\t"\
		"addpd	%%xmm5,		%%xmm5	\n\t"\
		"addpd	%%xmm4,		%%xmm4	\n\t"\
		"movaps	%%xmm0,    (%%ecx)	\n\t"\
		"movaps	%%xmm1,0x10(%%edx)	\n\t"\
		"addpd	%%xmm0,		%%xmm5	\n\t"\
		"addpd	%%xmm1,		%%xmm4	\n\t"\
		"movaps	%%xmm5,    (%%edx)	\n\t"\
		"movaps	%%xmm4,0x10(%%ecx)	\n\t"\
		"addl	%%edi,%%eax			\n\t"\
		"addl	%%edi,%%ebx			\n\t"\
		"addl	%%edi,%%ecx			\n\t"\
		"addl	%%edi,%%edx			\n\t"\
		"movl	%[__r17],%%esi		\n\t"\
		"movaps	     (%%esi),%%xmm0	\n\t"\
		"movaps	0x010(%%esi),%%xmm1	\n\t"\
		"movaps	0x040(%%esi),%%xmm2	\n\t"\
		"movaps	0x050(%%esi),%%xmm3	\n\t"\
		"subpd	0x050(%%esi),%%xmm0	\n\t"\
		"subpd	0x040(%%esi),%%xmm1	\n\t"\
		"addpd	0x010(%%esi),%%xmm2	\n\t"\
		"addpd	     (%%esi),%%xmm3	\n\t"\
		"movaps	0x020(%%esi),%%xmm4	\n\t"\
		"movaps	0x030(%%esi),%%xmm5	\n\t"\
		"movaps	0x060(%%esi),%%xmm6	\n\t"\
		"movaps	0x070(%%esi),%%xmm7	\n\t"\
		"subpd	0x030(%%esi),%%xmm4	\n\t"\
		"addpd	0x020(%%esi),%%xmm5	\n\t"\
		"addpd	0x070(%%esi),%%xmm6	\n\t"\
		"subpd	0x060(%%esi),%%xmm7	\n\t"\
		"movl	%[__isrt2],%%esi	\n\t"\
		"mulpd	(%%esi),%%xmm4		\n\t"\
		"mulpd	(%%esi),%%xmm5		\n\t"\
		"mulpd	(%%esi),%%xmm6		\n\t"\
		"mulpd	(%%esi),%%xmm7		\n\t"\
		"subpd	%%xmm6,%%xmm4		\n\t"\
		"subpd	%%xmm7,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm6		\n\t"\
		"addpd	%%xmm5,%%xmm7		\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm2		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"movaps	%%xmm0,    (%%ebx)	\n\t"\
		"movaps	%%xmm2,0x10(%%ebx)	\n\t"\
		"addpd	%%xmm0,		%%xmm4	\n\t"\
		"addpd	%%xmm2,		%%xmm5	\n\t"\
		"movaps	%%xmm4,    (%%eax)	\n\t"\
		"movaps	%%xmm5,0x10(%%eax)	\n\t"\
		"subpd	%%xmm7,		%%xmm3	\n\t"\
		"subpd	%%xmm6,		%%xmm1	\n\t"\
		"addpd	%%xmm7,		%%xmm7	\n\t"\
		"addpd	%%xmm6,		%%xmm6	\n\t"\
		"movaps	%%xmm3,    (%%ecx)	\n\t"\
		"movaps	%%xmm1,0x10(%%edx)	\n\t"\
		"addpd	%%xmm3,		%%xmm7	\n\t"\
		"addpd	%%xmm1,		%%xmm6	\n\t"\
		"movaps	%%xmm7,    (%%edx)	\n\t"\
		"movaps	%%xmm6,0x10(%%ecx)	\n\t"\
		"addl	%%edi,%%eax			\n\t"\
		"addl	%%edi,%%ebx			\n\t"\
		"addl	%%edi,%%ecx			\n\t"\
		"addl	%%edi,%%edx			\n\t"\
		"movl	%[__r9] ,%%esi		\n\t"\
		"movl	%[__cc0],%%edi		\n\t"\
		"movaps	0x020(%%esi),%%xmm4	\n\t"\
		"movaps	0x060(%%esi),%%xmm6	\n\t"\
		"movaps	0x030(%%esi),%%xmm5	\n\t"\
		"movaps	0x070(%%esi),%%xmm7	\n\t"\
		"movaps	%%xmm4,%%xmm0		\n\t"\
		"movaps	%%xmm6,%%xmm2		\n\t"\
		"movaps	%%xmm5,%%xmm1		\n\t"\
		"movaps	%%xmm7,%%xmm3		\n\t"\
		"mulpd	     (%%edi),%%xmm4	\n\t"\
		"mulpd	0x010(%%edi),%%xmm6	\n\t"\
		"mulpd	0x010(%%edi),%%xmm1	\n\t"\
		"mulpd	     (%%edi),%%xmm3	\n\t"\
		"mulpd	     (%%edi),%%xmm5	\n\t"\
		"mulpd	0x010(%%edi),%%xmm7	\n\t"\
		"mulpd	0x010(%%edi),%%xmm0	\n\t"\
		"mulpd	     (%%edi),%%xmm2	\n\t"\
		"subpd	%%xmm1,%%xmm4		\n\t"\
		"subpd	%%xmm3,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm5		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"subpd	%%xmm6,%%xmm4		\n\t"\
		"subpd	%%xmm7,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm6		\n\t"\
		"addpd	%%xmm5,%%xmm7		\n\t"\
		"movl	%[__isrt2],%%edi	\n\t"\
		"movaps	0x040(%%esi),%%xmm2	\n\t"\
		"movaps	0x050(%%esi),%%xmm3	\n\t"\
		"subpd	0x050(%%esi),%%xmm2	\n\t"\
		"addpd	0x040(%%esi),%%xmm3	\n\t"\
		"mulpd	(%%edi),%%xmm2		\n\t"\
		"mulpd	(%%edi),%%xmm3		\n\t"\
		"movaps	     (%%esi),%%xmm0	\n\t"\
		"movaps	0x010(%%esi),%%xmm1	\n\t"\
		"movl	%[__p4],%%edi		\n\t"\
		"shll	$3,%%edi			\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"addpd	     (%%esi),%%xmm2	\n\t"\
		"addpd	0x010(%%esi),%%xmm3	\n\t"\
		"subpd	%%xmm6,%%xmm2		\n\t"\
		"subpd	%%xmm7,%%xmm3		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"movaps	%%xmm2,    (%%ebx)	\n\t"\
		"movaps	%%xmm3,0x10(%%ebx)	\n\t"\
		"addpd	%%xmm2,		%%xmm6	\n\t"\
		"addpd	%%xmm3,		%%xmm7	\n\t"\
		"movaps	%%xmm6,    (%%eax)	\n\t"\
		"movaps	%%xmm7,0x10(%%eax)	\n\t"\
		"subpd	%%xmm5,		%%xmm0	\n\t"\
		"subpd	%%xmm4,		%%xmm1	\n\t"\
		"addpd	%%xmm5,		%%xmm5	\n\t"\
		"addpd	%%xmm4,		%%xmm4	\n\t"\
		"movaps	%%xmm0,    (%%ecx)	\n\t"\
		"movaps	%%xmm1,0x10(%%edx)	\n\t"\
		"addpd	%%xmm0,		%%xmm5	\n\t"\
		"addpd	%%xmm1,		%%xmm4	\n\t"\
		"movaps	%%xmm5,    (%%edx)	\n\t"\
		"movaps	%%xmm4,0x10(%%ecx)	\n\t"\
		"addl	%%edi,%%eax			\n\t"\
		"addl	%%edi,%%ebx			\n\t"\
		"addl	%%edi,%%ecx			\n\t"\
		"addl	%%edi,%%edx			\n\t"\
		"movl	%[__r25],%%esi		\n\t"\
		"movl	%[__cc0],%%edi		\n\t"\
		"movaps	0x020(%%esi),%%xmm4	\n\t"\
		"movaps	0x060(%%esi),%%xmm6	\n\t"\
		"movaps	0x030(%%esi),%%xmm5	\n\t"\
		"movaps	0x070(%%esi),%%xmm7	\n\t"\
		"movaps	%%xmm4,%%xmm0		\n\t"\
		"movaps	%%xmm6,%%xmm2		\n\t"\
		"movaps	%%xmm5,%%xmm1		\n\t"\
		"movaps	%%xmm7,%%xmm3		\n\t"\
		"mulpd	0x010(%%edi),%%xmm4	\n\t"\
		"mulpd	     (%%edi),%%xmm6	\n\t"\
		"mulpd	     (%%edi),%%xmm1	\n\t"\
		"mulpd	0x010(%%edi),%%xmm3	\n\t"\
		"mulpd	0x010(%%edi),%%xmm5	\n\t"\
		"mulpd	     (%%edi),%%xmm7	\n\t"\
		"mulpd	     (%%edi),%%xmm0	\n\t"\
		"mulpd	0x010(%%edi),%%xmm2	\n\t"\
		"subpd	%%xmm1,%%xmm4		\n\t"\
		"subpd	%%xmm3,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm5		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"subpd	%%xmm6,%%xmm4		\n\t"\
		"subpd	%%xmm7,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm6		\n\t"\
		"addpd	%%xmm5,%%xmm7		\n\t"\
		"movl	%[__isrt2],%%edi	\n\t"\
		"movaps	0x040(%%esi),%%xmm2	\n\t"\
		"movaps	0x050(%%esi),%%xmm3	\n\t"\
		"addpd	0x050(%%esi),%%xmm2	\n\t"\
		"subpd	0x040(%%esi),%%xmm3	\n\t"\
		"mulpd	(%%edi),%%xmm2		\n\t"\
		"mulpd	(%%edi),%%xmm3		\n\t"\
		"movaps	     (%%esi),%%xmm0	\n\t"\
		"movaps	0x010(%%esi),%%xmm1	\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"addpd	     (%%esi),%%xmm2	\n\t"\
		"addpd	0x010(%%esi),%%xmm3	\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"movaps	%%xmm0,    (%%ebx)	\n\t"\
		"movaps	%%xmm1,0x10(%%ebx)	\n\t"\
		"addpd	%%xmm0,		%%xmm4	\n\t"\
		"addpd	%%xmm1,		%%xmm5	\n\t"\
		"movaps	%%xmm4,    (%%eax)	\n\t"\
		"movaps	%%xmm5,0x10(%%eax)	\n\t"\
		"subpd	%%xmm7,		%%xmm2	\n\t"\
		"subpd	%%xmm6,		%%xmm3	\n\t"\
		"addpd	%%xmm7,		%%xmm7	\n\t"\
		"addpd	%%xmm6,		%%xmm6	\n\t"\
		"movaps	%%xmm2,    (%%ecx)	\n\t"\
		"movaps	%%xmm3,0x10(%%edx)	\n\t"\
		"addpd	%%xmm2,		%%xmm7	\n\t"\
		"addpd	%%xmm3,		%%xmm6	\n\t"\
		"movaps	%%xmm7,    (%%edx)	\n\t"\
		"movaps	%%xmm6,0x10(%%ecx)	\n\t"\
	"popl %%ebx	\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__p1] "m" (Xp1)\
		 ,[__p2] "m" (Xp2)\
		 ,[__p3] "m" (Xp3)\
		 ,[__p4] "m" (Xp4)\
		 ,[__p8] "m" (Xp8)\
		 ,[__r1] "m" (Xr1)\
		 ,[__r3] "m" (Xr3)\
		 ,[__r5] "m" (Xr5)\
		 ,[__r7] "m" (Xr7)\
		 ,[__r9] "m" (Xr9)\
		 ,[__r11] "m" (Xr11)\
		 ,[__r13] "m" (Xr13)\
		 ,[__r15] "m" (Xr15)\
		 ,[__r17] "m" (Xr17)\
		 ,[__r19] "m" (Xr19)\
		 ,[__r21] "m" (Xr21)\
		 ,[__r23] "m" (Xr23)\
		 ,[__r25] "m" (Xr25)\
		 ,[__r27] "m" (Xr27)\
		 ,[__r29] "m" (Xr29)\
		 ,[__r31] "m" (Xr31)\
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__cc0] "m" (Xcc0)\
		: "cc","memory","eax",/*"ebx",*/"ecx","edx","edi","esi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"	/* Clobbered registers */\
	);\
	}

#else

	#error Unhandled combination of #defs!

#endif	// SSE2 or AVX?

#endif	/* radix16_ditN_cy_dif1_gcc_h_included */

