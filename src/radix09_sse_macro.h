/*******************************************************************************
*                                                                              *
*   (C) 1997-2020 by Ernst W. Mayer.                                           *
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
#ifndef radix09_sse_macro_h_included
#define radix09_sse_macro_h_included

#include "sse2_macro_gcc64.h"

#ifdef USE_ARM_V8_SIMD

	/*...Radix-9 DIF: Ins in memory locations __i0-8.\
		Outs at memory locations __o0-8, assumed disjoint with inputs:\
	*/\
  #ifdef DFT9_V1	// Basic version:
	#define SSE2_RADIX_09_DIF(Xi0,Xi1,Xi2,Xi3,Xi4,Xi5,Xi6,Xi7,Xi8, Xc1, Xo0,Xo1,Xo2,Xo3,Xo4,Xo5,Xo6,Xo7,Xo8)\
	{\
	__asm__ volatile (\
		"ldr	x3,%[__c1]		\n\t	ldp	q30,q31,[x3,#0x40]	\n\t"/* cc1 */\
	/* Block 1: */\
		"ldr	x1,%[__i3]			\n\t"\
		"ldr	x2,%[__i6]			\n\t	ldp	q4,q5,[x1]		\n\t"\
		"ldr	x0,%[__i0]			\n\t	ldp	q6,q7,[x2]		\n\t"\
		"fadd	v8.2d,v4.2d,v6.2d	\n\t"\
		"fadd	v9.2d,v5.2d,v7.2d	\n\t	ldp	q0,q1,[x0]		\n\t"\
		"fsub	v4.2d,v4.2d,v6.2d	\n\t"\
		"fsub	v5.2d,v5.2d,v7.2d	\n\t	ldr	x0,%[__o0]		\n\t"\
		"fadd	v0.2d,v0.2d,v8.2d	\n\t"\
		"fadd	v1.2d,v1.2d,v9.2d	\n\t"\
		"stp	q0,q1,[x0]		\n\t"\
		"mov	v2.16b,v0.16b	\n\t"\
		"mov	v3.16b,v1.16b	\n\t"\
		"fmla	v2.2d,v8.2d,v30.2d	\n\t"\
		"fmla	v3.2d,v9.2d,v30.2d	\n\t"\
		"mov	v0.16b,v2.16b	\n\t"\
		"mov	v1.16b,v3.16b	\n\t"\
		"fmls	v0.2d,v31.2d,v5.2d	\n\t	ldr	x2,%[__o2]		\n\t"\
		"fmla	v1.2d,v31.2d,v4.2d	\n\t"\
		"fmla	v2.2d,v31.2d,v5.2d	\n\t	ldr	x1,%[__o1]		\n\t"\
		"fmls	v3.2d,v31.2d,v4.2d	\n\t"\
		"stp	q0,q1,[x1]			\n\t"\
		"stp	q2,q3,[x2]			\n\t"\
	/* Block 2: */\
		"ldr	x1,%[__i4]			\n\t"\
		"ldr	x2,%[__i7]			\n\t	ldp	q4,q5,[x1]		\n\t"\
		"ldr	x0,%[__i1]			\n\t	ldp	q6,q7,[x2]		\n\t"\
		"fadd	v8.2d,v4.2d,v6.2d	\n\t"\
		"fadd	v9.2d,v5.2d,v7.2d	\n\t	ldp	q0,q1,[x0]		\n\t"\
		"fsub	v4.2d,v4.2d,v6.2d	\n\t"\
		"fsub	v5.2d,v5.2d,v7.2d	\n\t	ldr	x0,%[__o3]		\n\t"\
		"fadd	v0.2d,v0.2d,v8.2d	\n\t"\
		"fadd	v1.2d,v1.2d,v9.2d	\n\t"\
		"stp	q0,q1,[x0]		\n\t"\
		"mov	v2.16b,v0.16b	\n\t"\
		"mov	v3.16b,v1.16b	\n\t"\
		"fmla	v2.2d,v8.2d,v30.2d	\n\t"\
		"fmla	v3.2d,v9.2d,v30.2d	\n\t"\
		"mov	v0.16b,v2.16b	\n\t"\
		"mov	v1.16b,v3.16b	\n\t"\
		"fmls	v0.2d,v31.2d,v5.2d	\n\t	ldr	x2,%[__o5]		\n\t"\
		"fmla	v1.2d,v31.2d,v4.2d	\n\t"\
		"fmla	v2.2d,v31.2d,v5.2d	\n\t	ldr	x1,%[__o4]		\n\t"\
		"fmls	v3.2d,v31.2d,v4.2d	\n\t"\
		"stp	q0,q1,[x1]			\n\t"\
		"stp	q2,q3,[x2]			\n\t"\
	/* Block 3: */\
		"ldr	x1,%[__i5]			\n\t"\
		"ldr	x2,%[__i8]			\n\t	ldp	q4,q5,[x1]		\n\t"\
		"ldr	x0,%[__i2]			\n\t	ldp	q6,q7,[x2]		\n\t"\
		"fadd	v8.2d,v4.2d,v6.2d	\n\t"\
		"fadd	v9.2d,v5.2d,v7.2d	\n\t	ldp	q0,q1,[x0]		\n\t"\
		"fsub	v4.2d,v4.2d,v6.2d	\n\t"\
		"fsub	v5.2d,v5.2d,v7.2d	\n\t	ldr	x0,%[__o6]		\n\t"\
		"fadd	v0.2d,v0.2d,v8.2d	\n\t"\
		"fadd	v1.2d,v1.2d,v9.2d	\n\t"\
		"stp	q0,q1,[x0]		\n\t"\
		"mov	v2.16b,v0.16b	\n\t"\
		"mov	v3.16b,v1.16b	\n\t"\
		"fmla	v2.2d,v8.2d,v30.2d	\n\t"\
		"fmla	v3.2d,v9.2d,v30.2d	\n\t"\
		"mov	v0.16b,v2.16b	\n\t"\
		"mov	v1.16b,v3.16b	\n\t"\
		"fmls	v0.2d,v31.2d,v5.2d	\n\t	ldr	x2,%[__o8]		\n\t"\
		"fmla	v1.2d,v31.2d,v4.2d	\n\t"\
		"fmla	v2.2d,v31.2d,v5.2d	\n\t	ldr	x1,%[__o7]		\n\t"\
		"fmls	v3.2d,v31.2d,v4.2d	\n\t"\
		"stp	q0,q1,[x1]			\n\t"\
		"stp	q2,q3,[x2]			\n\t"\
	/**********************************************************************************************/\
	/* now do three more radix-3 transforms, including 2 internal twiddles in each of blocks 2,3: */\
	/**********************************************************************************************/\
	/* Block 1: */\
		"ldr	x1,%[__o3]			\n\t"\
		"ldr	x2,%[__o6]			\n\t	ldp	q4,q5,[x1]		\n\t"\
		"ldr	x0,%[__o0]			\n\t	ldp	q6,q7,[x2]		\n\t"\
		"fadd	v8.2d,v4.2d,v6.2d	\n\t"\
		"fadd	v9.2d,v5.2d,v7.2d	\n\t	ldp	q0,q1,[x0]		\n\t"\
		"fsub	v4.2d,v4.2d,v6.2d	\n\t"\
		"fsub	v5.2d,v5.2d,v7.2d	\n\t"\
		"fadd	v0.2d,v0.2d,v8.2d	\n\t"\
		"fadd	v1.2d,v1.2d,v9.2d	\n\t"\
		"stp	q0,q1,[x0]		\n\t"\
		"mov	v2.16b,v0.16b	\n\t"\
		"mov	v3.16b,v1.16b	\n\t"\
		"fmla	v2.2d,v8.2d,v30.2d	\n\t"\
		"fmla	v3.2d,v9.2d,v30.2d	\n\t"\
		"mov	v0.16b,v2.16b	\n\t"\
		"mov	v1.16b,v3.16b	\n\t"\
		"fmls	v0.2d,v31.2d,v5.2d	\n\t"\
		"fmla	v1.2d,v31.2d,v4.2d	\n\t"\
		"fmla	v2.2d,v31.2d,v5.2d	\n\t"\
		"fmls	v3.2d,v31.2d,v4.2d	\n\t"\
		"stp	q0,q1,[x1]			\n\t"\
		"stp	q2,q3,[x2]			\n\t"\
	/* Block 2: */\
		"ldr	x1,%[__o4]			\n\t ldp q2,q3,[x3      ]	\n\t"/* c1 */\
		"ldr	x2,%[__o7]			\n\t	ldp	q8,q9,[x1]		\n\t"\
		"ldr	x0,%[__o1]			\n\t	ldp	q0,q1,[x2]		\n\t"\
		"fmul	v4.2d,v8.2d,v2.2d	\n\t"\
		"fmul	v5.2d,v9.2d,v2.2d	\n\t"\
		"fmls	v4.2d,v9.2d,v3.2d	\n\t"\
		"fmla	v5.2d,v8.2d,v3.2d	\n\t ldp q2,q3,[x3,#0x20]	\n\t"/* c2 */\
		"fmul	v6.2d,v0.2d,v2.2d	\n\t"\
		"fmul	v7.2d,v1.2d,v2.2d	\n\t"\
		"fmls	v6.2d,v1.2d,v3.2d	\n\t"\
		"fmla	v7.2d,v0.2d,v3.2d	\n\t"\
		"fadd	v8.2d,v4.2d,v6.2d	\n\t"\
		"fadd	v9.2d,v5.2d,v7.2d	\n\t	ldp	q0,q1,[x0]		\n\t"\
		"fsub	v4.2d,v4.2d,v6.2d	\n\t"\
		"fsub	v5.2d,v5.2d,v7.2d	\n\t"\
		"fadd	v0.2d,v0.2d,v8.2d	\n\t"\
		"fadd	v1.2d,v1.2d,v9.2d	\n\t"\
		"stp	q0,q1,[x0]		\n\t"\
		"mov	v2.16b,v0.16b	\n\t"\
		"mov	v3.16b,v1.16b	\n\t"\
		"fmla	v2.2d,v8.2d,v30.2d	\n\t"\
		"fmla	v3.2d,v9.2d,v30.2d	\n\t"\
		"mov	v0.16b,v2.16b	\n\t"\
		"mov	v1.16b,v3.16b	\n\t"\
		"fmls	v0.2d,v31.2d,v5.2d	\n\t"\
		"fmla	v1.2d,v31.2d,v4.2d	\n\t"\
		"fmla	v2.2d,v31.2d,v5.2d	\n\t"\
		"fmls	v3.2d,v31.2d,v4.2d	\n\t"\
		"stp	q0,q1,[x1]			\n\t"\
		"stp	q2,q3,[x2]			\n\t"\
	/* Block 3: */\
		"ldr	x1,%[__o5]			\n\t ldp q2,q3,[x3,#0x20]	\n\t"/* c2 */\
		"ldr	x2,%[__o8]			\n\t	ldp	q8,q9,[x1]		\n\t"\
		"ldr	x0,%[__o2]			\n\t	ldp	q0,q1,[x2]		\n\t"\
		"fmul	v4.2d,v8.2d,v2.2d	\n\t"\
		"fmul	v5.2d,v9.2d,v2.2d	\n\t"\
		"fmls	v4.2d,v9.2d,v3.2d	\n\t"\
		"fmla	v5.2d,v8.2d,v3.2d	\n\t ldp q2,q3,[x3,#0x60]	\n\t"/* c4 */\
		"fmul	v6.2d,v0.2d,v2.2d	\n\t"\
		"fmul	v7.2d,v1.2d,v2.2d	\n\t"\
		"fmls	v6.2d,v1.2d,v3.2d	\n\t"\
		"fmla	v7.2d,v0.2d,v3.2d	\n\t"\
		"fadd	v8.2d,v4.2d,v6.2d	\n\t"\
		"fadd	v9.2d,v5.2d,v7.2d	\n\t	ldp	q0,q1,[x0]		\n\t"\
		"fsub	v4.2d,v4.2d,v6.2d	\n\t"\
		"fsub	v5.2d,v5.2d,v7.2d	\n\t"\
		"fadd	v0.2d,v0.2d,v8.2d	\n\t"\
		"fadd	v1.2d,v1.2d,v9.2d	\n\t"\
		"stp	q0,q1,[x0]		\n\t"\
		"mov	v2.16b,v0.16b	\n\t"\
		"mov	v3.16b,v1.16b	\n\t"\
		"fmla	v2.2d,v8.2d,v30.2d	\n\t"\
		"fmla	v3.2d,v9.2d,v30.2d	\n\t"\
		"mov	v0.16b,v2.16b	\n\t"\
		"mov	v1.16b,v3.16b	\n\t"\
		"fmls	v0.2d,v31.2d,v5.2d	\n\t"\
		"fmla	v1.2d,v31.2d,v4.2d	\n\t"\
		"fmla	v2.2d,v31.2d,v5.2d	\n\t"\
		"fmls	v3.2d,v31.2d,v4.2d	\n\t"\
		"stp	q0,q1,[x1]			\n\t"\
		"stp	q2,q3,[x2]			\n\t"\
		:					/* outputs: none */\
		: [__i0] "m" (Xi0)	/* All inputs from memory addresses here */\
		 ,[__i1] "m" (Xi1)\
		 ,[__i2] "m" (Xi2)\
		 ,[__i3] "m" (Xi3)\
		 ,[__i4] "m" (Xi4)\
		 ,[__i5] "m" (Xi5)\
		 ,[__i6] "m" (Xi6)\
		 ,[__i7] "m" (Xi7)\
		 ,[__i8] "m" (Xi8)\
		 ,[__c1] "m" (Xc1)\
		 ,[__o0] "m" (Xo0)\
		 ,[__o1] "m" (Xo1)\
		 ,[__o2] "m" (Xo2)\
		 ,[__o3] "m" (Xo3)\
		 ,[__o4] "m" (Xo4)\
		 ,[__o5] "m" (Xo5)\
		 ,[__o6] "m" (Xo6)\
		 ,[__o7] "m" (Xo7)\
		 ,[__o8] "m" (Xo8)\
		: "cc","memory","x0","x1","x2","x3","v0","v1","v2","v3","v4","v5","v6","v7","v8","v9", "v30","v31"	/* Clobbered registers */\
		);\
	}
  #else	// Version 2 folds 3 3-DFTs in each of pass1 and 2 into side-by-side 3-column form.
		// This runs insignificantly faster on my Odroid C2, so don;t bother doign similar with the 9-DIT macro.
	#define SSE2_RADIX_09_DIF(Xi0,Xi1,Xi2,Xi3,Xi4,Xi5,Xi6,Xi7,Xi8, Xc1, Xo0,Xo1,Xo2,Xo3,Xo4,Xo5,Xo6,Xo7,Xo8)\
	{\
	__asm__ volatile (\
		"ldr	x9,%[__c1]		\n\t	ldp	q30,q31,[x9,#0x40]	\n\t"/* cc1 */\
		/* Block 1: */						/* Block 2: */							/* Block 3: */\
		"ldr	x1,%[__i3]			\n\t	ldr	x4,%[__i4]					\n\t	ldr	x7,%[__i5]					\n\t"\
		"ldr	x2,%[__i6]			\n\t	ldr	x5,%[__i7]					\n\t	ldr	x8,%[__i8]					\n\t"\
		"ldp	q4,q5,[x1]			\n\t	ldp	q14,q15,[x4]				\n\t	ldp	q24,q25,[x7]				\n\t"\
		"ldr	x0,%[__i0]			\n\t	ldr	x3,%[__i1]					\n\t	ldr	x6,%[__i2]					\n\t"\
		"ldp	q6,q7,[x2]			\n\t	ldp	q16,q17,[x5]				\n\t	ldp	q26,q27,[x8]				\n\t"\
		"fadd	v8.2d,v4.2d,v6.2d	\n\t	fadd	v18.2d,v14.2d,v16.2d	\n\t	fadd	v28.2d,v24.2d,v26.2d	\n\t"\
		"fadd	v9.2d,v5.2d,v7.2d	\n\t	fadd	v19.2d,v15.2d,v17.2d	\n\t	fadd	v29.2d,v25.2d,v27.2d	\n\t"\
		"ldp	q0,q1,[x0]			\n\t	ldp	q10,q11,[x3]				\n\t	ldp	q20,q21,[x6]				\n\t"\
		"fsub	v4.2d,v4.2d,v6.2d	\n\t	fsub	v14.2d,v14.2d,v16.2d	\n\t	fsub	v24.2d,v24.2d,v26.2d	\n\t"\
		"fsub	v5.2d,v5.2d,v7.2d	\n\t	fsub	v15.2d,v15.2d,v17.2d	\n\t	fsub	v25.2d,v25.2d,v27.2d	\n\t"\
		"ldr	x0,%[__o0]			\n\t	ldr	x3,%[__o3]					\n\t	ldr	x6,%[__o6]					\n\t"\
		"fadd	v0.2d,v0.2d,v8.2d	\n\t	fadd	v10.2d,v10.2d,v18.2d	\n\t	fadd	v20.2d,v20.2d,v28.2d	\n\t"\
		"fadd	v1.2d,v1.2d,v9.2d	\n\t	fadd	v11.2d,v11.2d,v19.2d	\n\t	fadd	v21.2d,v21.2d,v29.2d	\n\t"\
		"stp	q0,q1,[x0]			\n\t	stp	q10,q11,[x3]				\n\t	stp	q20,q21,[x6]				\n\t"\
		"mov	v2.16b,v0.16b		\n\t	mov	v12.16b,v10.16b				\n\t	mov	v22.16b,v20.16b				\n\t"\
		"mov	v3.16b,v1.16b		\n\t	mov	v13.16b,v11.16b				\n\t	mov	v23.16b,v21.16b				\n\t"\
		"fmla	v2.2d,v8.2d,v30.2d	\n\t	fmla	v12.2d,v18.2d,v30.2d	\n\t	fmla	v22.2d,v28.2d,v30.2d	\n\t"\
		"fmla	v3.2d,v9.2d,v30.2d	\n\t	fmla	v13.2d,v19.2d,v30.2d	\n\t	fmla	v23.2d,v29.2d,v30.2d	\n\t"\
		"mov	v0.16b,v2.16b		\n\t	mov	v10.16b,v12.16b				\n\t	mov	v20.16b,v22.16b				\n\t"\
		"mov	v1.16b,v3.16b		\n\t	mov	v11.16b,v13.16b				\n\t	mov	v21.16b,v23.16b				\n\t"\
		"fmls	v0.2d,v31.2d,v5.2d	\n\t	fmls	v10.2d,v31.2d,v15.2d	\n\t	fmls	v20.2d,v31.2d,v25.2d	\n\t"\
		"ldr	x2,%[__o2]			\n\t	ldr	x5,%[__o5]					\n\t	ldr	x8,%[__o8]					\n\t"\
		"fmla	v1.2d,v31.2d,v4.2d	\n\t	fmla	v11.2d,v31.2d,v14.2d	\n\t	fmla	v21.2d,v31.2d,v24.2d	\n\t"\
		"fmla	v2.2d,v31.2d,v5.2d	\n\t	fmla	v12.2d,v31.2d,v15.2d	\n\t	fmla	v22.2d,v31.2d,v25.2d	\n\t"\
		"ldr	x1,%[__o1]			\n\t	ldr	x4,%[__o4]					\n\t	ldr	x7,%[__o7]					\n\t"\
		"fmls	v3.2d,v31.2d,v4.2d	\n\t	fmls	v13.2d,v31.2d,v14.2d	\n\t	fmls	v23.2d,v31.2d,v24.2d	\n\t"\
		"stp	q0,q1,[x1]			\n\t	stp	q10,q11,[x4]				\n\t	stp	q20,q21,[x7]				\n\t"\
		"stp	q2,q3,[x2]			\n\t	stp	q12,q13,[x5]				\n\t	stp	q22,q23,[x8]				\n\t"\
	/**********************************************************************************************/\
	/* now do three more radix-3 transforms, including 2 internal twiddles in each of blocks 2,3: */\
	/**********************************************************************************************/\
		/* Block 1: */						/* Block 2: */							/* Block 3: */\
		"ldr	x1,%[__o3]			\n\t	ldr	x4,%[__o4]					\n\t	ldr	x7,%[__o5]					\n\t"\
		"ldr	x2,%[__o6]			\n\t	ldp q12,q13,[x9      ]	/* c1 */\n\t	ldp q22,q23,[x9,#0x20]	/* c2 */\n\t"\
		"ldp	q4,q5,[x1]			\n\t	ldr	x5,%[__o7]					\n\t	ldr	x8,%[__o8]					\n\t"\
		"ldr	x0,%[__o0]			\n\t	ldp	q18,q19,[x4]				\n\t	ldp	q28,q29,[x7]				\n\t"\
		"ldp	q6,q7,[x2]			\n\t	ldr	x3,%[__o1]					\n\t	ldr	x6,%[__o2]					\n\t"\
		"fadd	v8.2d,v4.2d,v6.2d	\n\t	ldp	q10,q11,[x5]				\n\t	ldp	q20,q21,[x8]				\n\t"\
		"fadd	v9.2d,v5.2d,v7.2d	\n\t	fmul	v14.2d,v18.2d,v12.2d	\n\t	fmul	v24.2d,v28.2d,v22.2d	\n\t"\
		"ldp	q0,q1,[x0]			\n\t	fmul	v15.2d,v19.2d,v12.2d	\n\t	fmul	v25.2d,v29.2d,v22.2d	\n\t"\
		"fsub	v4.2d,v4.2d,v6.2d	\n\t	fmls	v14.2d,v19.2d,v13.2d	\n\t	fmls	v24.2d,v29.2d,v23.2d	\n\t"\
		"fsub	v5.2d,v5.2d,v7.2d	\n\t	fmla	v15.2d,v18.2d,v13.2d	\n\t	fmla	v25.2d,v28.2d,v23.2d	\n\t"\
		"fadd	v0.2d,v0.2d,v8.2d	\n\t	ldp q12,q13,[x9,#0x20]	/* c2 */\n\t	ldp q22,q23,[x9,#0x60]	/* c4 */\n\t"\
		"fadd	v1.2d,v1.2d,v9.2d	\n\t	fmul	v16.2d,v10.2d,v12.2d	\n\t	fmul	v26.2d,v20.2d,v22.2d	\n\t"\
		"stp	q0,q1,[x0]			\n\t	fmul	v17.2d,v11.2d,v12.2d	\n\t	fmul	v27.2d,v21.2d,v22.2d	\n\t"\
		"mov	v2.16b,v0.16b		\n\t	fmls	v16.2d,v11.2d,v13.2d	\n\t	fmls	v26.2d,v21.2d,v23.2d	\n\t"\
		"mov	v3.16b,v1.16b		\n\t	fmla	v17.2d,v10.2d,v13.2d	\n\t	fmla	v27.2d,v20.2d,v23.2d	\n\t"\
		"fmla	v2.2d,v8.2d,v30.2d	\n\t	fadd	v18.2d,v14.2d,v16.2d	\n\t	fadd	v28.2d,v24.2d,v26.2d	\n\t"\
		"fmla	v3.2d,v9.2d,v30.2d	\n\t	fadd	v19.2d,v15.2d,v17.2d	\n\t	fadd	v29.2d,v25.2d,v27.2d	\n\t"\
		"mov	v0.16b,v2.16b		\n\t	ldp	q10,q11,[x3]				\n\t	ldp	q20,q21,[x6]				\n\t"\
		"mov	v1.16b,v3.16b		\n\t	fsub	v14.2d,v14.2d,v16.2d	\n\t	fsub	v24.2d,v24.2d,v26.2d	\n\t"\
		"fmls	v0.2d,v31.2d,v5.2d	\n\t	fsub	v15.2d,v15.2d,v17.2d	\n\t	fsub	v25.2d,v25.2d,v27.2d	\n\t"\
		"fmla	v1.2d,v31.2d,v4.2d	\n\t	fadd	v10.2d,v10.2d,v18.2d	\n\t	fadd	v20.2d,v20.2d,v28.2d	\n\t"\
		"fmla	v2.2d,v31.2d,v5.2d	\n\t	fadd	v11.2d,v11.2d,v19.2d	\n\t	fadd	v21.2d,v21.2d,v29.2d	\n\t"\
		"fmls	v3.2d,v31.2d,v4.2d	\n\t	stp	q10,q11,[x3]				\n\t	stp	q20,q21,[x6]				\n\t"\
		"stp	q0,q1,[x1]			\n\t	mov	v12.16b,v10.16b				\n\t	mov	v22.16b,v20.16b				\n\t"\
		"stp	q2,q3,[x2]			\n\t	mov	v13.16b,v11.16b				\n\t	mov	v23.16b,v21.16b				\n\t"\
		"									fmla	v12.2d,v18.2d,v30.2d	\n\t	fmla	v22.2d,v28.2d,v30.2d	\n\t"\
		"									fmla	v13.2d,v19.2d,v30.2d	\n\t	fmla	v23.2d,v29.2d,v30.2d	\n\t"\
		"									mov	v10.16b,v12.16b				\n\t	mov	v20.16b,v22.16b				\n\t"\
		"									mov	v11.16b,v13.16b				\n\t	mov	v21.16b,v23.16b				\n\t"\
		"									fmls	v10.2d,v31.2d,v15.2d	\n\t	fmls	v20.2d,v31.2d,v25.2d	\n\t"\
		"									fmla	v11.2d,v31.2d,v14.2d	\n\t	fmla	v21.2d,v31.2d,v24.2d	\n\t"\
		"									fmla	v12.2d,v31.2d,v15.2d	\n\t	fmla	v22.2d,v31.2d,v25.2d	\n\t"\
		"									fmls	v13.2d,v31.2d,v14.2d	\n\t	fmls	v23.2d,v31.2d,v24.2d	\n\t"\
		"									stp	q10,q11,[x4]				\n\t	stp	q20,q21,[x7]				\n\t"\
		"									stp	q12,q13,[x5]				\n\t	stp	q22,q23,[x8]				\n\t"\
		:					/* outputs: none */\
		: [__i0] "m" (Xi0)	/* All inputs from memory addresses here */\
		 ,[__i1] "m" (Xi1)\
		 ,[__i2] "m" (Xi2)\
		 ,[__i3] "m" (Xi3)\
		 ,[__i4] "m" (Xi4)\
		 ,[__i5] "m" (Xi5)\
		 ,[__i6] "m" (Xi6)\
		 ,[__i7] "m" (Xi7)\
		 ,[__i8] "m" (Xi8)\
		 ,[__c1] "m" (Xc1)\
		 ,[__o0] "m" (Xo0)\
		 ,[__o1] "m" (Xo1)\
		 ,[__o2] "m" (Xo2)\
		 ,[__o3] "m" (Xo3)\
		 ,[__o4] "m" (Xo4)\
		 ,[__o5] "m" (Xo5)\
		 ,[__o6] "m" (Xo6)\
		 ,[__o7] "m" (Xo7)\
		 ,[__o8] "m" (Xo8)\
		: "cc","memory","x0","x1","x2","x3","x4","x5","x6","x7","x8","x9",\
						"v0","v1","v2","v3","v4","v5","v6","v7","v8","v9",\
						"v10","v11","v12","v13","v14","v15","v16","v17","v18","v19",\
						"v20","v21","v22","v23","v24","v25","v26","v27","v28","v29","v30","v31"	/* Clobbered registers */\
		);\
	}
  #endif

	/*...Radix-9 DIT: Ins in memory locations __i0-8.\
		Outs at memory locations __o0-8, assumed disjoint with inputs:\
	*/\
	#define SSE2_RADIX_09_DIT(Xi0,Xi1,Xi2,Xi3,Xi4,Xi5,Xi6,Xi7,Xi8, Xc1, Xo0,Xo1,Xo2,Xo3,Xo4,Xo5,Xo6,Xo7,Xo8)\
	{\
	__asm__ volatile (\
		"ldr	x3,%[__c1]		\n\t	ldp	q30,q31,[x3,#0x40]	\n\t"/* cc1 */\
	/* Block 1: */\
		"ldr	x1,%[__i1]			\n\t"\
		"ldr	x2,%[__i2]			\n\t	ldp	q4,q5,[x1]		\n\t"\
		"ldr	x0,%[__i0]			\n\t	ldp	q6,q7,[x2]		\n\t"\
		"fadd	v8.2d,v4.2d,v6.2d	\n\t"\
		"fadd	v9.2d,v5.2d,v7.2d	\n\t	ldp	q0,q1,[x0]		\n\t"\
		"fsub	v4.2d,v4.2d,v6.2d	\n\t"\
		"fsub	v5.2d,v5.2d,v7.2d	\n\t"\
		"fadd	v0.2d,v0.2d,v8.2d	\n\t"\
		"fadd	v1.2d,v1.2d,v9.2d	\n\t"\
		"stp	q0,q1,[x0]		\n\t"\
		"mov	v2.16b,v0.16b	\n\t"\
		"mov	v3.16b,v1.16b	\n\t"\
		"fmla	v2.2d,v8.2d,v30.2d	\n\t"\
		"fmla	v3.2d,v9.2d,v30.2d	\n\t"\
		"mov	v0.16b,v2.16b	\n\t"\
		"mov	v1.16b,v3.16b	\n\t"\
		"fmla	v0.2d,v31.2d,v5.2d	\n\t"\
		"fmls	v1.2d,v31.2d,v4.2d	\n\t"\
		"fmls	v2.2d,v31.2d,v5.2d	\n\t"\
		"fmla	v3.2d,v31.2d,v4.2d	\n\t"\
		"stp	q0,q1,[x1]			\n\t"\
		"stp	q2,q3,[x2]			\n\t"\
	/* Block 2: */\
		"ldr	x1,%[__i4]			\n\t"\
		"ldr	x2,%[__i5]			\n\t	ldp	q4,q5,[x1]		\n\t"\
		"ldr	x0,%[__i3]			\n\t	ldp	q6,q7,[x2]		\n\t"\
		"fadd	v8.2d,v4.2d,v6.2d	\n\t"\
		"fadd	v9.2d,v5.2d,v7.2d	\n\t	ldp	q0,q1,[x0]		\n\t"\
		"fsub	v4.2d,v4.2d,v6.2d	\n\t"\
		"fsub	v5.2d,v5.2d,v7.2d	\n\t"\
		"fadd	v0.2d,v0.2d,v8.2d	\n\t"\
		"fadd	v1.2d,v1.2d,v9.2d	\n\t"\
		"stp	q0,q1,[x0]		\n\t"\
		"mov	v2.16b,v0.16b	\n\t"\
		"mov	v3.16b,v1.16b	\n\t"\
		"fmla	v2.2d,v8.2d,v30.2d	\n\t"\
		"fmla	v3.2d,v9.2d,v30.2d	\n\t"\
		"mov	v0.16b,v2.16b	\n\t"\
		"mov	v1.16b,v3.16b	\n\t"\
		"fmla	v0.2d,v31.2d,v5.2d	\n\t"\
		"fmls	v1.2d,v31.2d,v4.2d	\n\t"\
		"fmls	v2.2d,v31.2d,v5.2d	\n\t"\
		"fmla	v3.2d,v31.2d,v4.2d	\n\t"\
		"stp	q0,q1,[x1]			\n\t"\
		"stp	q2,q3,[x2]			\n\t"\
	/* Block 3: */\
		"ldr	x1,%[__i7]			\n\t"\
		"ldr	x2,%[__i8]			\n\t	ldp	q4,q5,[x1]		\n\t"\
		"ldr	x0,%[__i6]			\n\t	ldp	q6,q7,[x2]		\n\t"\
		"fadd	v8.2d,v4.2d,v6.2d	\n\t"\
		"fadd	v9.2d,v5.2d,v7.2d	\n\t	ldp	q0,q1,[x0]		\n\t"\
		"fsub	v4.2d,v4.2d,v6.2d	\n\t"\
		"fsub	v5.2d,v5.2d,v7.2d	\n\t"\
		"fadd	v0.2d,v0.2d,v8.2d	\n\t"\
		"fadd	v1.2d,v1.2d,v9.2d	\n\t"\
		"stp	q0,q1,[x0]		\n\t"\
		"mov	v2.16b,v0.16b	\n\t"\
		"mov	v3.16b,v1.16b	\n\t"\
		"fmla	v2.2d,v8.2d,v30.2d	\n\t"\
		"fmla	v3.2d,v9.2d,v30.2d	\n\t"\
		"mov	v0.16b,v2.16b	\n\t"\
		"mov	v1.16b,v3.16b	\n\t"\
		"fmla	v0.2d,v31.2d,v5.2d	\n\t"\
		"fmls	v1.2d,v31.2d,v4.2d	\n\t"\
		"fmls	v2.2d,v31.2d,v5.2d	\n\t"\
		"fmla	v3.2d,v31.2d,v4.2d	\n\t"\
		"stp	q0,q1,[x1]			\n\t"\
		"stp	q2,q3,[x2]			\n\t"\
	/**********************************************************************************************/\
	/* now do three more radix-3 transforms, including 2 internal twiddles in each of blocks 2,3: */\
	/**********************************************************************************************/\
	/* Block 1: */\
		"ldr	x1,%[__i3]			\n\t"\
		"ldr	x2,%[__i6]			\n\t	ldp	q4,q5,[x1]		\n\t"\
		"ldr	x0,%[__i0]			\n\t	ldp	q6,q7,[x2]		\n\t"\
		"fadd	v8.2d,v4.2d,v6.2d	\n\t"\
		"fadd	v9.2d,v5.2d,v7.2d	\n\t	ldp	q0,q1,[x0]		\n\t"\
		"fsub	v4.2d,v4.2d,v6.2d	\n\t"\
		"fsub	v5.2d,v5.2d,v7.2d	\n\t	ldr	x0,%[__o0]		\n\t"\
		"fadd	v0.2d,v0.2d,v8.2d	\n\t"\
		"fadd	v1.2d,v1.2d,v9.2d	\n\t"\
		"stp	q0,q1,[x0]		\n\t"\
		"mov	v2.16b,v0.16b	\n\t"\
		"mov	v3.16b,v1.16b	\n\t"\
		"fmla	v2.2d,v8.2d,v30.2d	\n\t"\
		"fmla	v3.2d,v9.2d,v30.2d	\n\t"\
		"mov	v0.16b,v2.16b	\n\t"\
		"mov	v1.16b,v3.16b	\n\t"\
		"fmla	v0.2d,v31.2d,v5.2d	\n\t	ldr	x2,%[__o6]		\n\t"\
		"fmls	v1.2d,v31.2d,v4.2d	\n\t"\
		"fmls	v2.2d,v31.2d,v5.2d	\n\t	ldr	x1,%[__o3]		\n\t"\
		"fmla	v3.2d,v31.2d,v4.2d	\n\t"\
		"stp	q0,q1,[x1]			\n\t"\
		"stp	q2,q3,[x2]			\n\t"\
	/* Block 2: */\
		"ldr	x1,%[__i4]			\n\t ldp q28,q29,[x3      ]	\n\t"/* c1 */\
		"ldr	x2,%[__i7]			\n\t	ldp	q8,q9,[x1]		\n\t"\
		"ldr	x0,%[__i1]			\n\t	ldp	q0,q1,[x2]		\n\t"\
		"fmul	v4.2d,v8.2d,v28.2d	\n\t"\
		"fmul	v5.2d,v9.2d,v28.2d	\n\t"\
		"fmla	v4.2d,v9.2d,v29.2d	\n\t"\
		"fmls	v5.2d,v8.2d,v29.2d	\n\t ldp q28,q29,[x3,#0x20]	\n\t"/* c2 */\
		"fmul	v6.2d,v0.2d,v28.2d	\n\t"\
		"fmul	v7.2d,v1.2d,v28.2d	\n\t"\
		"fmla	v6.2d,v1.2d,v29.2d	\n\t"\
		"fmls	v7.2d,v0.2d,v29.2d	\n\t"\
		"fadd	v8.2d,v4.2d,v6.2d	\n\t"\
		"fadd	v9.2d,v5.2d,v7.2d	\n\t	ldp	q0,q1,[x0]		\n\t"\
		"fsub	v4.2d,v4.2d,v6.2d	\n\t"\
		"fsub	v5.2d,v5.2d,v7.2d	\n\t	ldr	x0,%[__o7]		\n\t"\
		"fadd	v0.2d,v0.2d,v8.2d	\n\t"\
		"fadd	v1.2d,v1.2d,v9.2d	\n\t"\
		"stp	q0,q1,[x0]		\n\t"\
		"mov	v2.16b,v0.16b	\n\t"\
		"mov	v3.16b,v1.16b	\n\t"\
		"fmla	v2.2d,v8.2d,v30.2d	\n\t"\
		"fmla	v3.2d,v9.2d,v30.2d	\n\t"\
		"mov	v0.16b,v2.16b	\n\t"\
		"mov	v1.16b,v3.16b	\n\t"\
		"fmla	v0.2d,v31.2d,v5.2d	\n\t	ldr	x2,%[__o4]		\n\t"\
		"fmls	v1.2d,v31.2d,v4.2d	\n\t"\
		"fmls	v2.2d,v31.2d,v5.2d	\n\t	ldr	x1,%[__o1]		\n\t"\
		"fmla	v3.2d,v31.2d,v4.2d	\n\t"\
		"stp	q0,q1,[x1]			\n\t"\
		"stp	q2,q3,[x2]			\n\t"\
	/* Block 3: */\
		"ldr	x1,%[__i5]			\n\t ldp q28,q29,[x3,#0x20]	\n\t"/* c2 */\
		"ldr	x2,%[__i8]			\n\t	ldp	q8,q9,[x1]		\n\t"\
		"ldr	x0,%[__i2]			\n\t	ldp	q0,q1,[x2]		\n\t"\
		"fmul	v4.2d,v8.2d,v28.2d	\n\t"\
		"fmul	v5.2d,v9.2d,v28.2d	\n\t"\
		"fmla	v4.2d,v9.2d,v29.2d	\n\t"\
		"fmls	v5.2d,v8.2d,v29.2d	\n\t ldp q28,q29,[x3,#0x60]	\n\t"/* c4 */\
		"fmul	v6.2d,v0.2d,v28.2d	\n\t"\
		"fmul	v7.2d,v1.2d,v28.2d	\n\t"\
		"fmla	v6.2d,v1.2d,v29.2d	\n\t"\
		"fmls	v7.2d,v0.2d,v29.2d	\n\t"\
		"fadd	v8.2d,v4.2d,v6.2d	\n\t"\
		"fadd	v9.2d,v5.2d,v7.2d	\n\t	ldp	q0,q1,[x0]		\n\t"\
		"fsub	v4.2d,v4.2d,v6.2d	\n\t"\
		"fsub	v5.2d,v5.2d,v7.2d	\n\t	ldr	x0,%[__o5]		\n\t"\
		"fadd	v0.2d,v0.2d,v8.2d	\n\t"\
		"fadd	v1.2d,v1.2d,v9.2d	\n\t"\
		"stp	q0,q1,[x0]		\n\t"\
		"mov	v2.16b,v0.16b	\n\t"\
		"mov	v3.16b,v1.16b	\n\t"\
		"fmla	v2.2d,v8.2d,v30.2d	\n\t"\
		"fmla	v3.2d,v9.2d,v30.2d	\n\t"\
		"mov	v0.16b,v2.16b	\n\t"\
		"mov	v1.16b,v3.16b	\n\t"\
		"fmla	v0.2d,v31.2d,v5.2d	\n\t	ldr	x2,%[__o2]		\n\t"\
		"fmls	v1.2d,v31.2d,v4.2d	\n\t"\
		"fmls	v2.2d,v31.2d,v5.2d	\n\t	ldr	x1,%[__o8]		\n\t"\
		"fmla	v3.2d,v31.2d,v4.2d	\n\t"\
		"stp	q0,q1,[x1]			\n\t"\
		"stp	q2,q3,[x2]			\n\t"\
		:					/* outputs: none */\
		: [__i0] "m" (Xi0)	/* All inputs from memory addresses here */\
		 ,[__i1] "m" (Xi1)\
		 ,[__i2] "m" (Xi2)\
		 ,[__i3] "m" (Xi3)\
		 ,[__i4] "m" (Xi4)\
		 ,[__i5] "m" (Xi5)\
		 ,[__i6] "m" (Xi6)\
		 ,[__i7] "m" (Xi7)\
		 ,[__i8] "m" (Xi8)\
		 ,[__c1] "m" (Xc1)\
		 ,[__o0] "m" (Xo0)\
		 ,[__o1] "m" (Xo1)\
		 ,[__o2] "m" (Xo2)\
		 ,[__o3] "m" (Xo3)\
		 ,[__o4] "m" (Xo4)\
		 ,[__o5] "m" (Xo5)\
		 ,[__o6] "m" (Xo6)\
		 ,[__o7] "m" (Xo7)\
		 ,[__o8] "m" (Xo8)\
		: "cc","memory","x0","x1","x2","x3","v0","v1","v2","v3","v4","v5","v6","v7","v8","v9", "v28","v29","v30","v31"	/* Clobbered registers */\
		);\
	}

#elif defined(USE_AVX512)	// Use AVX2/FMA3 macros as starting point for these

	//...Radix-9 DIF: First set of Ins in memory locations __i0-8; Outs at memory locations __o0-8,
	// assumed disjoint with inputs.
	// Due to GCC macro argc limit of 30, to enable 16-register data-doubled version
	// of the radix-9 macros, use 2 length-9 ptr arrays for second set of IOs:
	//
	#define SSE2_RADIX_09_DIF_X2(Xi0,Xi1,Xi2,Xi3,Xi4,Xi5,Xi6,Xi7,Xi8, Xc1,Xtwo, Xo0,Xo1,Xo2,Xo3,Xo4,Xo5,Xo6,Xo7,Xo8, Xiptr,Xoptr)\
	{\
	__asm__ volatile (\
	"movq	%[__two],%%rcx	\n\t vmovaps	(%%rcx),%%zmm31	\n\t"/* 2.0 */\
																/*********** 2nd set of data: ********/\
	"movq	%[__i0]	,%%rax	\n\t"/* __i0; r[abc]x store set1 I-addrs */"movq	%[__iptr],%%r8		\n\t"/* Pointer to  input-ptr array */\
	"movq	%[__o0]	,%%rsi	\n\t"/* __o0; rsi,rdi store set1 O-addrs */"movq	%[__optr],%%r9		\n\t"/* Pointer to output-ptr array */\
	"movq	%[__c1]	,%%rdx	\n\t"/* rdx has trig addr for both sides */"movq	(%%r8),%%r10 		\n\t"/* __i0; r10-12 store  input addresses */\
	"																	movq	(%%r9),%%r13 		\n\t"/* __o0; r13,14 store output addresses */\
	/* Block 1: */\
		"movq	%%rax		,%%rbx							\n\t		movq	%%r10		,%%r11			\n\t"\
		"movq	%%rax		,%%rcx							\n\t		movq	%%r10		,%%r12			\n\t"\
		"movq	%[__i3]	,%%rbx 								\n\t		movq	0x18(%%r8),%%r11 			\n\t"/* i3 */\
		"movq	%[__i6]	,%%rcx 								\n\t		movq	0x30(%%r8),%%r12 			\n\t"/* i6 */\
		"vmovaps	    (%%rbx)	,%%zmm4						\n\t		vmovaps	    (%%r11)	,%%zmm12		\n\t"\
		"vmovaps	0x40(%%rbx)	,%%zmm5						\n\t		vmovaps	0x40(%%r11)	,%%zmm13		\n\t"\
		"vmovaps	    (%%rax)	,%%zmm0						\n\t		vmovaps	    (%%r10)	,%%zmm8 		\n\t"\
		"vmovaps	0x40(%%rax)	,%%zmm1						\n\t		vmovaps	0x40(%%r10)	,%%zmm9 		\n\t"\
		"vmovaps	    (%%rcx)	,%%zmm6						\n\t		vmovaps	    (%%r12)	,%%zmm14		\n\t"\
		"vmovaps	0x40(%%rcx)	,%%zmm7						\n\t		vmovaps	0x40(%%r12)	,%%zmm15		\n\t"\
		"vaddpd	%%zmm6,%%zmm4,%%zmm2						\n\t		vaddpd	%%zmm14,%%zmm12,%%zmm10		\n\t"\
		"vaddpd	%%zmm7,%%zmm5,%%zmm3						\n\t		vaddpd	%%zmm15,%%zmm13,%%zmm11		\n\t"\
		"vsubpd	%%zmm6,%%zmm4,%%zmm4						\n\t		vsubpd	%%zmm14,%%zmm12,%%zmm12		\n\t"\
		"vsubpd	%%zmm7,%%zmm5,%%zmm5						\n\t		vsubpd	%%zmm15,%%zmm13,%%zmm13		\n\t"\
		"vaddpd	%%zmm2,%%zmm0,%%zmm0						\n\t		vaddpd	%%zmm10,%%zmm8 ,%%zmm8 		\n\t"\
		"vaddpd	%%zmm3,%%zmm1,%%zmm1						\n\t		vaddpd	%%zmm11,%%zmm9 ,%%zmm9 		\n\t"\
		"vmovaps	0x100(%%rdx)	,%%zmm6			\n\t"/* c3m1 */\
		"vmovaps	0x140(%%rdx)	,%%zmm7			\n\t"\
		"vmovaps	%%zmm0		,    (%%rsi)	/* <- t00 */\n\t		vmovaps	%%zmm8 		,    (%%r13)	\n\t"/* <- t00 */\
		"vmovaps	%%zmm1		,0x40(%%rsi)	/* <- t01 */\n\t		vmovaps	%%zmm9 		,0x40(%%r13)	\n\t"/* <- t01 */\
	"movq	%[__o1]	,%%rdi									\n\t	movq	0x08(%%r9),%%r14				\n\t"/* o1 */\
	"movq	%[__o2]	,%%rsi									\n\t	movq	0x10(%%r9),%%r13				\n\t"/* o2 */\
	" vfmadd132pd	%%zmm6,%%zmm0,%%zmm2 					\n\t	 vfmadd132pd	%%zmm6,%%zmm8,%%zmm10	\n\t"\
	" vfmadd132pd	%%zmm6,%%zmm1,%%zmm3 					\n\t	 vfmadd132pd	%%zmm6,%%zmm9,%%zmm11	\n\t"\
		"vmovaps	%%zmm2,%%zmm0							\n\t		vmovaps	%%zmm10,%%zmm8 				\n\t"\
		"vmovaps	%%zmm3,%%zmm1							\n\t		vmovaps	%%zmm11,%%zmm9 				\n\t"\
	" vfmadd231pd	%%zmm7,%%zmm5,%%zmm0 					\n\t	 vfmadd231pd	%%zmm7,%%zmm13,%%zmm8 	\n\t"\
	"vfnmadd231pd	%%zmm7,%%zmm4,%%zmm1 					\n\t	vfnmadd231pd	%%zmm7,%%zmm12,%%zmm9 	\n\t"\
	"vfnmadd231pd	%%zmm7,%%zmm5,%%zmm2 					\n\t	vfnmadd231pd	%%zmm7,%%zmm13,%%zmm10	\n\t"\
	" vfmadd231pd	%%zmm7,%%zmm4,%%zmm3 					\n\t	 vfmadd231pd	%%zmm7,%%zmm12,%%zmm11	\n\t"\
		"vmovaps	%%zmm0		,    (%%rsi)	/* <- t04 */\n\t		vmovaps	%%zmm8 		,    (%%r13)	\n\t"/* <- t04 */\
		"vmovaps	%%zmm1		,0x40(%%rsi)	/* <- t05 */\n\t		vmovaps	%%zmm9 		,0x40(%%r13)	\n\t"/* <- t05 */\
		"vmovaps	%%zmm2		,    (%%rdi)	/* <- t02 */\n\t		vmovaps	%%zmm10		,    (%%r14)	\n\t"/* <- t02 */\
		"vmovaps	%%zmm3		,0x40(%%rdi)	/* <- t03 */\n\t		vmovaps	%%zmm11		,0x40(%%r14)	\n\t"/* <- t03 */\
	/* Block 2: */\
		"movq	%[__i1]	,%%rax 								\n\t		movq	0x08(%%r8),%%r10 			\n\t"/* i1 */\
		"movq	%[__i4]	,%%rbx 								\n\t		movq	0x20(%%r8),%%r11 			\n\t"/* i4 */\
		"movq	%[__i7]	,%%rcx 								\n\t		movq	0x38(%%r8),%%r12 			\n\t"/* i7 */\
		"vmovaps	    (%%rbx)	,%%zmm4						\n\t		vmovaps	    (%%r11)	,%%zmm12		\n\t"\
		"vmovaps	0x40(%%rbx)	,%%zmm5						\n\t		vmovaps	0x40(%%r11)	,%%zmm13		\n\t"\
		"vmovaps	    (%%rax)	,%%zmm0						\n\t		vmovaps	    (%%r10)	,%%zmm8 		\n\t"\
		"vmovaps	0x40(%%rax)	,%%zmm1						\n\t		vmovaps	0x40(%%r10)	,%%zmm9 		\n\t"\
		"vmovaps	    (%%rcx)	,%%zmm2						\n\t		vmovaps	    (%%r12)	,%%zmm10		\n\t"\
		"vmovaps	0x40(%%rcx)	,%%zmm3						\n\t		vmovaps	0x40(%%r12)	,%%zmm11		\n\t"\
	"movq	%[__o3]	,%%rdi									\n\t	movq	0x18(%%r9),%%r14				\n\t"/* o3 */\
		"vsubpd	%%zmm2,%%zmm4,%%zmm4						\n\t		vsubpd	%%zmm10,%%zmm12,%%zmm12		\n\t"\
		"vsubpd	%%zmm3,%%zmm5,%%zmm5						\n\t		vsubpd	%%zmm11,%%zmm13,%%zmm13		\n\t"\
	" vfmadd132pd	%%zmm31,%%zmm4,%%zmm2 					\n\t	 vfmadd132pd	%%zmm31,%%zmm12,%%zmm10	\n\t"\
	" vfmadd132pd	%%zmm31,%%zmm5,%%zmm3 					\n\t	 vfmadd132pd	%%zmm31,%%zmm13,%%zmm11	\n\t"\
		"vaddpd	%%zmm2,%%zmm0,%%zmm0						\n\t		vaddpd	%%zmm10,%%zmm8 ,%%zmm8 		\n\t"\
		"vaddpd	%%zmm3,%%zmm1,%%zmm1						\n\t		vaddpd	%%zmm11,%%zmm9 ,%%zmm9 		\n\t"\
		"vmovaps	%%zmm0		,    (%%rdi)	/* <- t06 */\n\t		vmovaps	%%zmm8 		,    (%%r14)	\n\t"/* <- t06 */\
		"vmovaps	%%zmm1		,0x40(%%rdi)	/* <- t07 */\n\t		vmovaps	%%zmm9 		,0x40(%%r14)	\n\t"/* <- t07 */\
	"movq	%[__o4]	,%%rsi									\n\t	movq	0x20(%%r9),%%r13				\n\t"/* o4 */\
	"movq	%[__o5]	,%%rdi									\n\t	movq	0x28(%%r9),%%r14				\n\t"/* o5 */\
	" vfmadd132pd	%%zmm6,%%zmm0,%%zmm2 					\n\t	 vfmadd132pd	%%zmm6,%%zmm8,%%zmm10	\n\t"\
	" vfmadd132pd	%%zmm6,%%zmm1,%%zmm3 					\n\t	 vfmadd132pd	%%zmm6,%%zmm9,%%zmm11	\n\t"\
		"vmovaps	%%zmm2,%%zmm0							\n\t		vmovaps	%%zmm10,%%zmm8 				\n\t"\
		"vmovaps	%%zmm3,%%zmm1							\n\t		vmovaps	%%zmm11,%%zmm9 				\n\t"\
	" vfmadd231pd	%%zmm7,%%zmm5,%%zmm0 					\n\t	 vfmadd231pd	%%zmm7,%%zmm13,%%zmm8 	\n\t"\
	"vfnmadd231pd	%%zmm7,%%zmm4,%%zmm1 					\n\t	vfnmadd231pd	%%zmm7,%%zmm12,%%zmm9 	\n\t"\
	"vfnmadd231pd	%%zmm7,%%zmm5,%%zmm2 					\n\t	vfnmadd231pd	%%zmm7,%%zmm13,%%zmm10	\n\t"\
	" vfmadd231pd	%%zmm7,%%zmm4,%%zmm3 					\n\t	 vfmadd231pd	%%zmm7,%%zmm12,%%zmm11	\n\t"\
		"vmovaps	%%zmm0		,    (%%rdi)	/* <- t0a */\n\t		vmovaps	%%zmm8 		,    (%%r14)	\n\t"/* <- t0a */\
		"vmovaps	%%zmm1		,0x40(%%rdi)	/* <- t0b */\n\t		vmovaps	%%zmm9 		,0x40(%%r14)	\n\t"/* <- t0b */\
		"vmovaps	%%zmm2		,    (%%rsi)	/* <- t08 */\n\t		vmovaps	%%zmm10		,    (%%r13)	\n\t"/* <- t08 */\
		"vmovaps	%%zmm3		,0x40(%%rsi)	/* <- t09 */\n\t		vmovaps	%%zmm11		,0x40(%%r13)	\n\t"/* <- t09 */\
	/* Block 3: */\
		"movq	%[__i2]	,%%rax 								\n\t		movq	0x10(%%r8),%%r10 			\n\t"/* i2 */\
		"movq	%[__i5]	,%%rbx 								\n\t		movq	0x28(%%r8),%%r11 			\n\t"/* i5 */\
		"movq	%[__i8]	,%%rcx 								\n\t		movq	0x40(%%r8),%%r12 			\n\t"/* i8 */\
		"vmovaps	    (%%rbx)	,%%zmm4						\n\t		vmovaps	    (%%r11)	,%%zmm12		\n\t"\
		"vmovaps	0x40(%%rbx)	,%%zmm5						\n\t		vmovaps	0x40(%%r11)	,%%zmm13		\n\t"\
		"vmovaps	    (%%rax)	,%%zmm0						\n\t		vmovaps	    (%%r10)	,%%zmm8 		\n\t"\
		"vmovaps	0x40(%%rax)	,%%zmm1						\n\t		vmovaps	0x40(%%r10)	,%%zmm9 		\n\t"\
		"vmovaps	    (%%rcx)	,%%zmm2						\n\t		vmovaps	    (%%r12)	,%%zmm10		\n\t"\
		"vmovaps	0x40(%%rcx)	,%%zmm3						\n\t		vmovaps	0x40(%%r12)	,%%zmm11		\n\t"\
	"movq	%[__o6]	,%%rsi									\n\t	movq	0x30(%%r9),%%r13				\n\t"/* o6 */\
		"vsubpd	%%zmm2,%%zmm4,%%zmm4						\n\t		vsubpd	%%zmm10,%%zmm12,%%zmm12		\n\t"\
		"vsubpd	%%zmm3,%%zmm5,%%zmm5						\n\t		vsubpd	%%zmm11,%%zmm13,%%zmm13		\n\t"\
	" vfmadd132pd	%%zmm31,%%zmm4,%%zmm2 					\n\t	 vfmadd132pd	%%zmm31,%%zmm12,%%zmm10	\n\t"\
	" vfmadd132pd	%%zmm31,%%zmm5,%%zmm3 					\n\t	 vfmadd132pd	%%zmm31,%%zmm13,%%zmm11	\n\t"\
		"vaddpd	%%zmm2,%%zmm0,%%zmm0						\n\t		vaddpd	%%zmm10,%%zmm8 ,%%zmm8 		\n\t"\
		"vaddpd	%%zmm3,%%zmm1,%%zmm1						\n\t		vaddpd	%%zmm11,%%zmm9 ,%%zmm9 		\n\t"\
		"vmovaps	%%zmm0		,    (%%rsi)	/* <- t0c */\n\t		vmovaps	%%zmm8 		,    (%%r13)	\n\t"/* <- t0c */\
		"vmovaps	%%zmm1		,0x40(%%rsi)	/* <- t0d */\n\t		vmovaps	%%zmm9 		,0x40(%%r13)	\n\t"/* <- t0d */\
	"movq	%[__o7]	,%%rdi									\n\t	movq	0x38(%%r9),%%r14				\n\t"/* o7 */\
	"movq	%[__o8]	,%%rsi									\n\t	movq	0x40(%%r9),%%r13				\n\t"/* o8 */\
	" vfmadd132pd	%%zmm6,%%zmm0,%%zmm2 					\n\t	 vfmadd132pd	%%zmm6,%%zmm8,%%zmm10	\n\t"\
	" vfmadd132pd	%%zmm6,%%zmm1,%%zmm3 					\n\t	 vfmadd132pd	%%zmm6,%%zmm9,%%zmm11	\n\t"\
		"vmovaps	%%zmm2,%%zmm0							\n\t		vmovaps	%%zmm10,%%zmm8 				\n\t"\
		"vmovaps	%%zmm3,%%zmm1							\n\t		vmovaps	%%zmm11,%%zmm9 				\n\t"\
	" vfmadd231pd	%%zmm7,%%zmm5,%%zmm0 					\n\t	 vfmadd231pd	%%zmm7,%%zmm13,%%zmm8 	\n\t"\
	"vfnmadd231pd	%%zmm7,%%zmm4,%%zmm1 					\n\t	vfnmadd231pd	%%zmm7,%%zmm12,%%zmm9 	\n\t"\
	"vfnmadd231pd	%%zmm7,%%zmm5,%%zmm2 					\n\t	vfnmadd231pd	%%zmm7,%%zmm13,%%zmm10	\n\t"\
	" vfmadd231pd	%%zmm7,%%zmm4,%%zmm3 					\n\t	 vfmadd231pd	%%zmm7,%%zmm12,%%zmm11	\n\t"\
		"vmovaps	%%zmm0		,    (%%rsi)	/* <- t0g */\n\t		vmovaps	%%zmm8 		,    (%%r13)	\n\t"/* <- t0g */\
		"vmovaps	%%zmm1		,0x40(%%rsi)	/* <- t0h */\n\t		vmovaps	%%zmm9 		,0x40(%%r13)	\n\t"/* <- t0h */\
		"vmovaps	%%zmm2		,    (%%rdi)	/* <- t0e */\n\t		vmovaps	%%zmm10		,    (%%r14)	\n\t"/* <- t0e */\
		"vmovaps	%%zmm3		,0x40(%%rdi)	/* <- t0f */\n\t		vmovaps	%%zmm11		,0x40(%%r14)	\n\t"/* <- t0f */\
	/******************************************************************************************/\
	/*                        now do three more radix-3 transforms:                           */\
	/******************************************************************************************/\
	/* Block 1: */													/* Block 1: */\
	"movq	%[__o0]	,%%rax									\n\t	movq	    (%%r9),%%r10				\n\t"/* o0 */\
	"movq	%[__o3]	,%%rbx									\n\t	movq	0x18(%%r9),%%r11				\n\t"/* o3 */\
	"movq	%[__o6]	,%%rcx									\n\t	movq	0x30(%%r9),%%r12				\n\t"/* o6 */\
		"vmovaps	    (%%rbx)	,%%zmm4						\n\t		vmovaps	    (%%r11)	,%%zmm12		\n\t"\
		"vmovaps	0x40(%%rbx)	,%%zmm5						\n\t		vmovaps	0x40(%%r11)	,%%zmm13		\n\t"\
		"vmovaps	    (%%rcx)	,%%zmm2						\n\t		vmovaps	    (%%r12)	,%%zmm10		\n\t"\
		"vmovaps	0x40(%%rcx)	,%%zmm3						\n\t		vmovaps	0x40(%%r12)	,%%zmm11		\n\t"\
		"vmovaps	    (%%rax)	,%%zmm0						\n\t		vmovaps	    (%%r10)	,%%zmm8 		\n\t"\
		"vmovaps	0x40(%%rax)	,%%zmm1						\n\t		vmovaps	0x40(%%r10)	,%%zmm9 		\n\t"\
		"vsubpd	%%zmm2,%%zmm4,%%zmm4						\n\t		vsubpd	%%zmm10,%%zmm12,%%zmm12		\n\t"\
		"vsubpd	%%zmm3,%%zmm5,%%zmm5						\n\t		vsubpd	%%zmm11,%%zmm13,%%zmm13		\n\t"\
	" vfmadd132pd	%%zmm31,%%zmm4,%%zmm2 					\n\t	 vfmadd132pd	%%zmm31,%%zmm12,%%zmm10	\n\t"\
	" vfmadd132pd	%%zmm31,%%zmm5,%%zmm3 					\n\t	 vfmadd132pd	%%zmm31,%%zmm13,%%zmm11	\n\t"\
		"vaddpd	%%zmm2,%%zmm0,%%zmm0						\n\t		vaddpd	%%zmm10,%%zmm8 ,%%zmm8 		\n\t"\
		"vaddpd	%%zmm3,%%zmm1,%%zmm1						\n\t		vaddpd	%%zmm11,%%zmm9 ,%%zmm9 		\n\t"\
		"vmovaps	%%zmm0		,    (%%rax)				\n\t		vmovaps	%%zmm8 		,    (%%r10)	\n\t"\
		"vmovaps	%%zmm1		,0x40(%%rax)				\n\t		vmovaps	%%zmm9 		,0x40(%%r10)	\n\t"\
	" vfmadd132pd	%%zmm6,%%zmm0,%%zmm2 					\n\t	 vfmadd132pd	%%zmm6,%%zmm8,%%zmm10	\n\t"\
	" vfmadd132pd	%%zmm6,%%zmm1,%%zmm3 					\n\t	 vfmadd132pd	%%zmm6,%%zmm9,%%zmm11	\n\t"\
		"vmovaps	%%zmm2,%%zmm0							\n\t		vmovaps	%%zmm10,%%zmm8 				\n\t"\
		"vmovaps	%%zmm3,%%zmm1							\n\t		vmovaps	%%zmm11,%%zmm9 				\n\t"\
	" vfmadd231pd	%%zmm7,%%zmm5,%%zmm0 					\n\t	 vfmadd231pd	%%zmm7,%%zmm13,%%zmm8 	\n\t"\
	"vfnmadd231pd	%%zmm7,%%zmm4,%%zmm1 					\n\t	vfnmadd231pd	%%zmm7,%%zmm12,%%zmm9 	\n\t"\
	"vfnmadd231pd	%%zmm7,%%zmm5,%%zmm2 					\n\t	vfnmadd231pd	%%zmm7,%%zmm13,%%zmm10	\n\t"\
	" vfmadd231pd	%%zmm7,%%zmm4,%%zmm3 					\n\t	 vfmadd231pd	%%zmm7,%%zmm12,%%zmm11	\n\t"\
		"vmovaps	%%zmm0		,    (%%rcx)				\n\t		vmovaps	%%zmm8 		,    (%%r12)	\n\t"\
		"vmovaps	%%zmm1		,0x40(%%rcx)				\n\t		vmovaps	%%zmm9 		,0x40(%%r12)	\n\t"\
		"vmovaps	%%zmm2		,    (%%rbx)				\n\t		vmovaps	%%zmm10		,    (%%r11)	\n\t"\
		"vmovaps	%%zmm3		,0x40(%%rbx)				\n\t		vmovaps	%%zmm11		,0x40(%%r11)	\n\t"\
	/* Block 2: */\
	"movq	%[__o1]	,%%rax									\n\t	movq	0x08(%%r9),%%r10				\n\t"/* o1 */\
	"movq	%[__o4]	,%%rbx									\n\t	movq	0x20(%%r9),%%r11				\n\t"/* o4 */\
	"movq	%[__o7]	,%%rcx									\n\t	movq	0x38(%%r9),%%r12				\n\t"/* o7 */\
		"vmovaps	    (%%rbx)	,%%zmm2						\n\t		vmovaps	    (%%r11)	,%%zmm10		\n\t"\
		"vmovaps	0x40(%%rbx)	,%%zmm3						\n\t		vmovaps	0x40(%%r11)	,%%zmm11		\n\t"\
		"vmovaps	    (%%rdx)	,%%zmm6				\n\t"/* c1 */\
		"vmovaps	0x40(%%rdx)	,%%zmm7				\n\t"\
		"vmovaps	%%zmm2		,%%zmm0						\n\t		vmovaps	%%zmm10		,%%zmm8 		\n\t"\
		"vmovaps	%%zmm3		,%%zmm1						\n\t		vmovaps	%%zmm11		,%%zmm9 		\n\t"\
		"vmulpd		%%zmm6,%%zmm2,%%zmm2					\n\t		vmulpd		%%zmm6,%%zmm10,%%zmm10	\n\t"\
		"vmulpd		%%zmm6,%%zmm3,%%zmm3					\n\t		vmulpd		%%zmm6,%%zmm11,%%zmm11	\n\t"\
	"vfnmadd231pd	%%zmm7,%%zmm1,%%zmm2 					\n\t	vfnmadd231pd	%%zmm7,%%zmm9 ,%%zmm10	\n\t"\
	" vfmadd231pd	%%zmm7,%%zmm0,%%zmm3 					\n\t	 vfmadd231pd	%%zmm7,%%zmm8 ,%%zmm11	\n\t"\
		"vmovaps	    (%%rcx)	,%%zmm4						\n\t		vmovaps	    (%%r12)	,%%zmm12		\n\t"\
		"vmovaps	0x40(%%rcx)	,%%zmm5						\n\t		vmovaps	0x40(%%r12)	,%%zmm13		\n\t"\
		"vmovaps	0x80(%%rdx)	,%%zmm6				\n\t"/* c2 */\
		"vmovaps	0xc0(%%rdx)	,%%zmm7				\n\t"\
		"vmovaps	%%zmm4		,%%zmm0						\n\t		vmovaps	%%zmm12		,%%zmm8 		\n\t"\
		"vmovaps	%%zmm5		,%%zmm1						\n\t		vmovaps	%%zmm13		,%%zmm9 		\n\t"\
		"vmulpd		%%zmm6,%%zmm4,%%zmm4					\n\t		vmulpd		%%zmm6,%%zmm12,%%zmm12	\n\t"\
		"vmulpd		%%zmm6,%%zmm5,%%zmm5					\n\t		vmulpd		%%zmm6,%%zmm13,%%zmm13	\n\t"\
	"vfnmadd231pd	%%zmm7,%%zmm1,%%zmm4 					\n\t	vfnmadd231pd	%%zmm7,%%zmm9 ,%%zmm12	\n\t"\
	" vfmadd231pd	%%zmm7,%%zmm0,%%zmm5 					\n\t	 vfmadd231pd	%%zmm7,%%zmm8 ,%%zmm13	\n\t"\
		"vmovaps	0x100(%%rdx)	,%%zmm6				\n\t"/* c3m1 */\
		"vmovaps	0x140(%%rdx)	,%%zmm7				\n\t"\
		"vmovaps	    (%%rax)	,%%zmm0						\n\t		vmovaps	    (%%r10)	,%%zmm8 		\n\t"\
		"vmovaps	0x40(%%rax)	,%%zmm1						\n\t		vmovaps	0x40(%%r10)	,%%zmm9 		\n\t"\
		"vsubpd	%%zmm4,%%zmm2,%%zmm2						\n\t		vsubpd	%%zmm12,%%zmm10,%%zmm10		\n\t"\
		"vsubpd	%%zmm5,%%zmm3,%%zmm3						\n\t		vsubpd	%%zmm13,%%zmm11,%%zmm11		\n\t"\
	" vfmadd132pd	%%zmm31,%%zmm2,%%zmm4 					\n\t	 vfmadd132pd	%%zmm31,%%zmm10,%%zmm12	\n\t"\
	" vfmadd132pd	%%zmm31,%%zmm3,%%zmm5 					\n\t	 vfmadd132pd	%%zmm31,%%zmm11,%%zmm13	\n\t"\
		"vaddpd	%%zmm4,%%zmm0,%%zmm0						\n\t		vaddpd	%%zmm12,%%zmm8 ,%%zmm8 		\n\t"\
		"vaddpd	%%zmm5,%%zmm1,%%zmm1						\n\t		vaddpd	%%zmm13,%%zmm9 ,%%zmm9 		\n\t"\
		"vmovaps	%%zmm0		,    (%%rax)				\n\t		vmovaps	%%zmm8 		,    (%%r10)	\n\t"\
		"vmovaps	%%zmm1		,0x40(%%rax)				\n\t		vmovaps	%%zmm9 		,0x40(%%r10)	\n\t"\
	" vfmadd132pd	%%zmm6,%%zmm0,%%zmm4 					\n\t	 vfmadd132pd	%%zmm6,%%zmm8,%%zmm12	\n\t"\
	" vfmadd132pd	%%zmm6,%%zmm1,%%zmm5 					\n\t	 vfmadd132pd	%%zmm6,%%zmm9,%%zmm13	\n\t"\
		"vmovaps	%%zmm4,%%zmm0							\n\t		vmovaps	%%zmm12,%%zmm8 				\n\t"\
		"vmovaps	%%zmm5,%%zmm1							\n\t		vmovaps	%%zmm13,%%zmm9 				\n\t"\
	" vfmadd231pd	%%zmm7,%%zmm3,%%zmm0 					\n\t	 vfmadd231pd	%%zmm7,%%zmm11,%%zmm8 	\n\t"\
	"vfnmadd231pd	%%zmm7,%%zmm2,%%zmm1 					\n\t	vfnmadd231pd	%%zmm7,%%zmm10,%%zmm9 	\n\t"\
	"vfnmadd231pd	%%zmm7,%%zmm3,%%zmm4 					\n\t	vfnmadd231pd	%%zmm7,%%zmm11,%%zmm12	\n\t"\
	" vfmadd231pd	%%zmm7,%%zmm2,%%zmm5 					\n\t	 vfmadd231pd	%%zmm7,%%zmm10,%%zmm13	\n\t"\
		"vmovaps	%%zmm0		,    (%%rcx)				\n\t		vmovaps	%%zmm8 		,    (%%r12)	\n\t"\
		"vmovaps	%%zmm1		,0x40(%%rcx)				\n\t		vmovaps	%%zmm9 		,0x40(%%r12)	\n\t"\
		"vmovaps	%%zmm4		,    (%%rbx)				\n\t		vmovaps	%%zmm12		,    (%%r11)	\n\t"\
		"vmovaps	%%zmm5		,0x40(%%rbx)				\n\t		vmovaps	%%zmm13		,0x40(%%r11)	\n\t"\
	/* Block 3: */\
	"movq	%[__o2]	,%%rax									\n\t	movq	0x10(%%r9),%%r10				\n\t"/* o2 */\
	"movq	%[__o5]	,%%rbx									\n\t	movq	0x28(%%r9),%%r11				\n\t"/* o5 */\
	"movq	%[__o8]	,%%rcx									\n\t	movq	0x40(%%r9),%%r12				\n\t"/* o8 */\
		"vmovaps	    (%%rbx)	,%%zmm2						\n\t		vmovaps	    (%%r11)	,%%zmm10		\n\t"\
		"vmovaps	0x40(%%rbx)	,%%zmm3						\n\t		vmovaps	0x40(%%r11)	,%%zmm11		\n\t"\
		"vmovaps	0x80(%%rdx)	,%%zmm6				\n\t"/* c2 */\
		"vmovaps	0xc0(%%rdx)	,%%zmm7				\n\t"\
		"vmovaps	%%zmm2		,%%zmm0						\n\t		vmovaps	%%zmm10		,%%zmm8 		\n\t"\
		"vmovaps	%%zmm3		,%%zmm1						\n\t		vmovaps	%%zmm11		,%%zmm9 		\n\t"\
		"vmulpd		%%zmm6,%%zmm2,%%zmm2					\n\t		vmulpd		%%zmm6,%%zmm10,%%zmm10	\n\t"\
		"vmulpd		%%zmm6,%%zmm3,%%zmm3					\n\t		vmulpd		%%zmm6,%%zmm11,%%zmm11	\n\t"\
	"vfnmadd231pd	%%zmm7,%%zmm1,%%zmm2 					\n\t	vfnmadd231pd	%%zmm7,%%zmm9 ,%%zmm10	\n\t"\
	" vfmadd231pd	%%zmm7,%%zmm0,%%zmm3 					\n\t	 vfmadd231pd	%%zmm7,%%zmm8 ,%%zmm11	\n\t"\
		"vmovaps	    (%%rcx)	,%%zmm4						\n\t		vmovaps	    (%%r12)	,%%zmm12		\n\t"\
		"vmovaps	0x40(%%rcx)	,%%zmm5						\n\t		vmovaps	0x40(%%r12)	,%%zmm13		\n\t"\
		"vmovaps	0x180(%%rdx)	,%%zmm6				\n\t"/* c4 */\
		"vmovaps	0x1c0(%%rdx)	,%%zmm7				\n\t"\
		"vmovaps	%%zmm4		,%%zmm0						\n\t		vmovaps	%%zmm12		,%%zmm8 		\n\t"\
		"vmovaps	%%zmm5		,%%zmm1						\n\t		vmovaps	%%zmm13		,%%zmm9 		\n\t"\
		"vmulpd		%%zmm6,%%zmm4,%%zmm4					\n\t		vmulpd		%%zmm6,%%zmm12,%%zmm12	\n\t"\
		"vmulpd		%%zmm6,%%zmm5,%%zmm5					\n\t		vmulpd		%%zmm6,%%zmm13,%%zmm13	\n\t"\
	"vfnmadd231pd	%%zmm7,%%zmm1,%%zmm4 					\n\t	vfnmadd231pd	%%zmm7,%%zmm9 ,%%zmm12	\n\t"\
	" vfmadd231pd	%%zmm7,%%zmm0,%%zmm5 					\n\t	 vfmadd231pd	%%zmm7,%%zmm8 ,%%zmm13	\n\t"\
		"vmovaps	0x100(%%rdx)	,%%zmm6				\n\t"/* c3m1 */\
		"vmovaps	0x140(%%rdx)	,%%zmm7				\n\t"\
		"vmovaps	    (%%rax)	,%%zmm0						\n\t		vmovaps	    (%%r10)	,%%zmm8 		\n\t"\
		"vmovaps	0x40(%%rax)	,%%zmm1						\n\t		vmovaps	0x40(%%r10)	,%%zmm9 		\n\t"\
		"vsubpd	%%zmm4,%%zmm2,%%zmm2						\n\t		vsubpd	%%zmm12,%%zmm10,%%zmm10		\n\t"\
		"vsubpd	%%zmm5,%%zmm3,%%zmm3						\n\t		vsubpd	%%zmm13,%%zmm11,%%zmm11		\n\t"\
	" vfmadd132pd	%%zmm31,%%zmm2,%%zmm4 					\n\t	 vfmadd132pd	%%zmm31,%%zmm10,%%zmm12	\n\t"\
	" vfmadd132pd	%%zmm31,%%zmm3,%%zmm5 					\n\t	 vfmadd132pd	%%zmm31,%%zmm11,%%zmm13	\n\t"\
		"vaddpd	%%zmm4,%%zmm0,%%zmm0						\n\t		vaddpd	%%zmm12,%%zmm8 ,%%zmm8 		\n\t"\
		"vaddpd	%%zmm5,%%zmm1,%%zmm1						\n\t		vaddpd	%%zmm13,%%zmm9 ,%%zmm9 		\n\t"\
		"vmovaps	%%zmm0		,    (%%rax)				\n\t		vmovaps	%%zmm8 		,    (%%r10)	\n\t"\
		"vmovaps	%%zmm1		,0x40(%%rax)				\n\t		vmovaps	%%zmm9 		,0x40(%%r10)	\n\t"\
	" vfmadd132pd	%%zmm6,%%zmm0,%%zmm4 					\n\t	 vfmadd132pd	%%zmm6,%%zmm8,%%zmm12	\n\t"\
	" vfmadd132pd	%%zmm6,%%zmm1,%%zmm5 					\n\t	 vfmadd132pd	%%zmm6,%%zmm9,%%zmm13	\n\t"\
		"vmovaps	%%zmm4,%%zmm0							\n\t		vmovaps	%%zmm12,%%zmm8 				\n\t"\
		"vmovaps	%%zmm5,%%zmm1							\n\t		vmovaps	%%zmm13,%%zmm9 				\n\t"\
	" vfmadd231pd	%%zmm7,%%zmm3,%%zmm0 					\n\t	 vfmadd231pd	%%zmm7,%%zmm11,%%zmm8 	\n\t"\
	"vfnmadd231pd	%%zmm7,%%zmm2,%%zmm1 					\n\t	vfnmadd231pd	%%zmm7,%%zmm10,%%zmm9 	\n\t"\
	"vfnmadd231pd	%%zmm7,%%zmm3,%%zmm4 					\n\t	vfnmadd231pd	%%zmm7,%%zmm11,%%zmm12	\n\t"\
	" vfmadd231pd	%%zmm7,%%zmm2,%%zmm5 					\n\t	 vfmadd231pd	%%zmm7,%%zmm10,%%zmm13	\n\t"\
		"vmovaps	%%zmm0		,    (%%rcx)				\n\t		vmovaps	%%zmm8 		,    (%%r12)	\n\t"\
		"vmovaps	%%zmm1		,0x40(%%rcx)				\n\t		vmovaps	%%zmm9 		,0x40(%%r12)	\n\t"\
		"vmovaps	%%zmm4		,    (%%rbx)				\n\t		vmovaps	%%zmm12		,    (%%r11)	\n\t"\
		"vmovaps	%%zmm5		,0x40(%%rbx)				\n\t		vmovaps	%%zmm13		,0x40(%%r11)	\n\t"\
		:					/* outputs: none */\
		: [__i0] "m" (Xi0)	/* All inputs from memory addresses here */\
		 ,[__i1] "m" (Xi1)\
		 ,[__i2] "m" (Xi2)\
		 ,[__i3] "m" (Xi3)\
		 ,[__i4] "m" (Xi4)\
		 ,[__i5] "m" (Xi5)\
		 ,[__i6] "m" (Xi6)\
		 ,[__i7] "m" (Xi7)\
		 ,[__i8] "m" (Xi8)\
		 ,[__c1] "m" (Xc1)\
		 ,[__two] "m" (Xtwo)\
		 ,[__o0] "m" (Xo0)\
		 ,[__o1] "m" (Xo1)\
		 ,[__o2] "m" (Xo2)\
		 ,[__o3] "m" (Xo3)\
		 ,[__o4] "m" (Xo4)\
		 ,[__o5] "m" (Xo5)\
		 ,[__o6] "m" (Xo6)\
		 ,[__o7] "m" (Xo7)\
		 ,[__o8] "m" (Xo8)\
		 ,[__iptr] "m" (Xiptr)	/* Due to GCC macro argc limit of 30, to enable 16-register data-doubled version */\
		 ,[__optr] "m" (Xoptr)	/* of the radix-9 macros use 2 length-9 ptr arrays for second set of IOs. */\
		: "cc","memory","rax","rbx","rcx","rdx","rsi","rdi","r8","r9","r10","r11","r12","r13","r14","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15","xmm31"	/* Clobbered registers */\
		);\
	}

	/*...Radix-9 DIT: Ins in memory locations __i0-8.\
		Outs at memory locations __o0-8, assumed disjoint with inputs:\
	*/\
	#define SSE2_RADIX_09_DIT_X2(Xi0,Xi1,Xi2,Xi3,Xi4,Xi5,Xi6,Xi7,Xi8, Xc1,Xtwo, Xo0,Xo1,Xo2,Xo3,Xo4,Xo5,Xo6,Xo7,Xo8, Xiptr,Xoptr)\
	{\
	__asm__ volatile (\
	"movq	%[__two],%%rcx	\n\t vmovaps	(%%rcx),%%zmm31	\n\t"/* 2.0 */\
																		/*********** 2nd set of data: ********/\
	"																	movq	%[__iptr],%%r8		\n\t"/* Pointer to  input-ptr array */\
	"																	movq	%[__optr],%%r9		\n\t"/* Pointer to output-ptr array */\
	"																	movq	(%%r8),%%r10 		\n\t"/* __i0; r10-12 store  input addresses */\
	"movq	%[__c1]	,%%rdx	\n\t"/* rdx has trig addr for both sides */"movq	(%%r9),%%r13 		\n\t"/* __o0; r13,14 store output addresses */\
	"movq	%[__i0]	,%%rax	\n\t"/* __i0; r[abc]x store set1 I-addrs */"movq	    (%%r8),%%r10				\n\t"/* i0 */\
	"movq	%[__i1]	,%%rbx	\n\t										movq	0x08(%%r8),%%r11				\n\t"/* i1 */\
	"movq	%[__i2]	,%%rcx	\n\t										movq	0x10(%%r8),%%r12				\n\t"/* i2 */\
	/* Block 1: */\
		"vmovaps	    (%%rbx)	,%%zmm4						\n\t		vmovaps	    (%%r11)	,%%zmm12		\n\t"\
		"vmovaps	0x40(%%rbx)	,%%zmm5						\n\t		vmovaps	0x40(%%r11)	,%%zmm13		\n\t"\
		"vmovaps	    (%%rax)	,%%zmm0						\n\t		vmovaps	    (%%r10)	,%%zmm8 		\n\t"\
		"vmovaps	0x40(%%rax)	,%%zmm1						\n\t		vmovaps	0x40(%%r10)	,%%zmm9 		\n\t"\
		"vmovaps	    (%%rcx)	,%%zmm6						\n\t		vmovaps	    (%%r12)	,%%zmm14		\n\t"\
		"vmovaps	0x40(%%rcx)	,%%zmm7						\n\t		vmovaps	0x40(%%r12)	,%%zmm15		\n\t"\
		"vaddpd	%%zmm6,%%zmm4,%%zmm2						\n\t		vaddpd	%%zmm14,%%zmm12,%%zmm10		\n\t"\
		"vaddpd	%%zmm7,%%zmm5,%%zmm3						\n\t		vaddpd	%%zmm15,%%zmm13,%%zmm11		\n\t"\
		"vsubpd	%%zmm6,%%zmm4,%%zmm4						\n\t		vsubpd	%%zmm14,%%zmm12,%%zmm12		\n\t"\
		"vsubpd	%%zmm7,%%zmm5,%%zmm5						\n\t		vsubpd	%%zmm15,%%zmm13,%%zmm13		\n\t"\
		"vaddpd	%%zmm2,%%zmm0,%%zmm0						\n\t		vaddpd	%%zmm10,%%zmm8 ,%%zmm8 		\n\t"\
		"vaddpd	%%zmm3,%%zmm1,%%zmm1						\n\t		vaddpd	%%zmm11,%%zmm9 ,%%zmm9 		\n\t"\
		"vmovaps	0x100(%%rdx)	,%%zmm6						\n\t"/* c3m1 */\
		"vmovaps	0x140(%%rdx)	,%%zmm7						\n\t"\
		"vmovaps	%%zmm0		,    (%%rax)				\n\t		vmovaps	%%zmm8 		,    (%%r10)	\n\t"\
		"vmovaps	%%zmm1		,0x40(%%rax)				\n\t		vmovaps	%%zmm9 		,0x40(%%r10)	\n\t"\
	" vfmadd132pd	%%zmm6,%%zmm0,%%zmm2 					\n\t	 vfmadd132pd	%%zmm6,%%zmm8,%%zmm10	\n\t"\
	" vfmadd132pd	%%zmm6,%%zmm1,%%zmm3 					\n\t	 vfmadd132pd	%%zmm6,%%zmm9,%%zmm11	\n\t"\
		"vmovaps	%%zmm2,%%zmm0							\n\t		vmovaps	%%zmm10,%%zmm8 				\n\t"\
		"vmovaps	%%zmm3,%%zmm1							\n\t		vmovaps	%%zmm11,%%zmm9 				\n\t"\
	"vfnmadd231pd	%%zmm7,%%zmm5,%%zmm0 					\n\t	vfnmadd231pd	%%zmm7,%%zmm13,%%zmm8 	\n\t"\
	" vfmadd231pd	%%zmm7,%%zmm4,%%zmm1 					\n\t	 vfmadd231pd	%%zmm7,%%zmm12,%%zmm9 	\n\t"\
	" vfmadd231pd	%%zmm7,%%zmm5,%%zmm2 					\n\t	 vfmadd231pd	%%zmm7,%%zmm13,%%zmm10	\n\t"\
	"vfnmadd231pd	%%zmm7,%%zmm4,%%zmm3 					\n\t	vfnmadd231pd	%%zmm7,%%zmm12,%%zmm11	\n\t"\
		"vmovaps	%%zmm0		,    (%%rcx)				\n\t		vmovaps	%%zmm8 		,    (%%r12)	\n\t"\
		"vmovaps	%%zmm1		,0x40(%%rcx)				\n\t		vmovaps	%%zmm9 		,0x40(%%r12)	\n\t"\
		"vmovaps	%%zmm2		,    (%%rbx)				\n\t		vmovaps	%%zmm10		,    (%%r11)	\n\t"\
		"vmovaps	%%zmm3		,0x40(%%rbx)				\n\t		vmovaps	%%zmm11		,0x40(%%r11)	\n\t"\
	/* Block 2: */\
	"movq	%[__i3]	,%%rax									\n\t	movq	0x18(%%r8),%%r10				\n\t"/* i3 */\
	"movq	%[__i4]	,%%rbx									\n\t	movq	0x20(%%r8),%%r11				\n\t"/* i4 */\
	"movq	%[__i5]	,%%rcx									\n\t	movq	0x28(%%r8),%%r12				\n\t"/* i5 */\
		"vmovaps	    (%%rbx)	,%%zmm4						\n\t		vmovaps	    (%%r11)	,%%zmm12		\n\t"\
		"vmovaps	0x40(%%rbx)	,%%zmm5						\n\t		vmovaps	0x40(%%r11)	,%%zmm13		\n\t"\
		"vmovaps	    (%%rcx)	,%%zmm2						\n\t		vmovaps	    (%%r12)	,%%zmm10		\n\t"\
		"vmovaps	0x40(%%rcx)	,%%zmm3						\n\t		vmovaps	0x40(%%r12)	,%%zmm11		\n\t"\
		"vmovaps	    (%%rax)	,%%zmm0						\n\t		vmovaps	    (%%r10)	,%%zmm8 		\n\t"\
		"vmovaps	0x40(%%rax)	,%%zmm1						\n\t		vmovaps	0x40(%%r10)	,%%zmm9 		\n\t"\
		"vsubpd	%%zmm2,%%zmm4,%%zmm4						\n\t		vsubpd	%%zmm10,%%zmm12,%%zmm12		\n\t"\
		"vsubpd	%%zmm3,%%zmm5,%%zmm5						\n\t		vsubpd	%%zmm11,%%zmm13,%%zmm13		\n\t"\
	" vfmadd132pd	%%zmm31,%%zmm4,%%zmm2 					\n\t	 vfmadd132pd	%%zmm31,%%zmm12,%%zmm10	\n\t"\
	" vfmadd132pd	%%zmm31,%%zmm5,%%zmm3 					\n\t	 vfmadd132pd	%%zmm31,%%zmm13,%%zmm11	\n\t"\
		"vaddpd	%%zmm2,%%zmm0,%%zmm0						\n\t		vaddpd	%%zmm10,%%zmm8 ,%%zmm8 		\n\t"\
		"vaddpd	%%zmm3,%%zmm1,%%zmm1						\n\t		vaddpd	%%zmm11,%%zmm9 ,%%zmm9 		\n\t"\
		"vmovaps	%%zmm0		,    (%%rax)				\n\t		vmovaps	%%zmm8 		,    (%%r10)	\n\t"\
		"vmovaps	%%zmm1		,0x40(%%rax)				\n\t		vmovaps	%%zmm9 		,0x40(%%r10)	\n\t"\
	" vfmadd132pd	%%zmm6,%%zmm0,%%zmm2 					\n\t	 vfmadd132pd	%%zmm6,%%zmm8,%%zmm10	\n\t"\
	" vfmadd132pd	%%zmm6,%%zmm1,%%zmm3 					\n\t	 vfmadd132pd	%%zmm6,%%zmm9,%%zmm11	\n\t"\
		"vmovaps	%%zmm2,%%zmm0							\n\t		vmovaps	%%zmm10,%%zmm8 				\n\t"\
		"vmovaps	%%zmm3,%%zmm1							\n\t		vmovaps	%%zmm11,%%zmm9 				\n\t"\
	"vfnmadd231pd	%%zmm7,%%zmm5,%%zmm0 					\n\t	vfnmadd231pd	%%zmm7,%%zmm13,%%zmm8 	\n\t"\
	" vfmadd231pd	%%zmm7,%%zmm4,%%zmm1 					\n\t	 vfmadd231pd	%%zmm7,%%zmm12,%%zmm9 	\n\t"\
	" vfmadd231pd	%%zmm7,%%zmm5,%%zmm2 					\n\t	 vfmadd231pd	%%zmm7,%%zmm13,%%zmm10	\n\t"\
	"vfnmadd231pd	%%zmm7,%%zmm4,%%zmm3 					\n\t	vfnmadd231pd	%%zmm7,%%zmm12,%%zmm11	\n\t"\
		"vmovaps	%%zmm0		,    (%%rcx)				\n\t		vmovaps	%%zmm8 		,    (%%r12)	\n\t"\
		"vmovaps	%%zmm1		,0x40(%%rcx)				\n\t		vmovaps	%%zmm9 		,0x40(%%r12)	\n\t"\
		"vmovaps	%%zmm2		,    (%%rbx)				\n\t		vmovaps	%%zmm10		,    (%%r11)	\n\t"\
		"vmovaps	%%zmm3		,0x40(%%rbx)				\n\t		vmovaps	%%zmm11		,0x40(%%r11)	\n\t"\
	/* Block 3: */\
	"movq	%[__i6]	,%%rax									\n\t	movq	0x30(%%r8),%%r10				\n\t"/* i6 */\
	"movq	%[__i7]	,%%rbx									\n\t	movq	0x38(%%r8),%%r11				\n\t"/* i7 */\
	"movq	%[__i8]	,%%rcx									\n\t	movq	0x40(%%r8),%%r12				\n\t"/* i8 */\
		"vmovaps	    (%%rbx)	,%%zmm4						\n\t		vmovaps	    (%%r11)	,%%zmm12		\n\t"\
		"vmovaps	0x40(%%rbx)	,%%zmm5						\n\t		vmovaps	0x40(%%r11)	,%%zmm13		\n\t"\
		"vmovaps	    (%%rcx)	,%%zmm2						\n\t		vmovaps	    (%%r12)	,%%zmm10		\n\t"\
		"vmovaps	0x40(%%rcx)	,%%zmm3						\n\t		vmovaps	0x40(%%r12)	,%%zmm11		\n\t"\
		"vmovaps	    (%%rax)	,%%zmm0						\n\t		vmovaps	    (%%r10)	,%%zmm8 		\n\t"\
		"vmovaps	0x40(%%rax)	,%%zmm1						\n\t		vmovaps	0x40(%%r10)	,%%zmm9 		\n\t"\
		"vsubpd	%%zmm2,%%zmm4,%%zmm4						\n\t		vsubpd	%%zmm10,%%zmm12,%%zmm12		\n\t"\
		"vsubpd	%%zmm3,%%zmm5,%%zmm5						\n\t		vsubpd	%%zmm11,%%zmm13,%%zmm13		\n\t"\
	" vfmadd132pd	%%zmm31,%%zmm4,%%zmm2 					\n\t	 vfmadd132pd	%%zmm31,%%zmm12,%%zmm10	\n\t"\
	" vfmadd132pd	%%zmm31,%%zmm5,%%zmm3 					\n\t	 vfmadd132pd	%%zmm31,%%zmm13,%%zmm11	\n\t"\
		"vaddpd	%%zmm2,%%zmm0,%%zmm0						\n\t		vaddpd	%%zmm10,%%zmm8 ,%%zmm8 		\n\t"\
		"vaddpd	%%zmm3,%%zmm1,%%zmm1						\n\t		vaddpd	%%zmm11,%%zmm9 ,%%zmm9 		\n\t"\
		"vmovaps	%%zmm0		,    (%%rax)				\n\t		vmovaps	%%zmm8 		,    (%%r10)	\n\t"\
		"vmovaps	%%zmm1		,0x40(%%rax)				\n\t		vmovaps	%%zmm9 		,0x40(%%r10)	\n\t"\
	" vfmadd132pd	%%zmm6,%%zmm0,%%zmm2 					\n\t	 vfmadd132pd	%%zmm6,%%zmm8,%%zmm10	\n\t"\
	" vfmadd132pd	%%zmm6,%%zmm1,%%zmm3 					\n\t	 vfmadd132pd	%%zmm6,%%zmm9,%%zmm11	\n\t"\
		"vmovaps	%%zmm2,%%zmm0							\n\t		vmovaps	%%zmm10,%%zmm8 				\n\t"\
		"vmovaps	%%zmm3,%%zmm1							\n\t		vmovaps	%%zmm11,%%zmm9 				\n\t"\
	"vfnmadd231pd	%%zmm7,%%zmm5,%%zmm0 					\n\t	vfnmadd231pd	%%zmm7,%%zmm13,%%zmm8 	\n\t"\
	" vfmadd231pd	%%zmm7,%%zmm4,%%zmm1 					\n\t	 vfmadd231pd	%%zmm7,%%zmm12,%%zmm9 	\n\t"\
	" vfmadd231pd	%%zmm7,%%zmm5,%%zmm2 					\n\t	 vfmadd231pd	%%zmm7,%%zmm13,%%zmm10	\n\t"\
	"vfnmadd231pd	%%zmm7,%%zmm4,%%zmm3 					\n\t	vfnmadd231pd	%%zmm7,%%zmm12,%%zmm11	\n\t"\
		"vmovaps	%%zmm0		,    (%%rcx)				\n\t		vmovaps	%%zmm8 		,    (%%r12)	\n\t"\
		"vmovaps	%%zmm1		,0x40(%%rcx)				\n\t		vmovaps	%%zmm9 		,0x40(%%r12)	\n\t"\
		"vmovaps	%%zmm2		,    (%%rbx)				\n\t		vmovaps	%%zmm10		,    (%%r11)	\n\t"\
		"vmovaps	%%zmm3		,0x40(%%rbx)				\n\t		vmovaps	%%zmm11		,0x40(%%r11)	\n\t"\
	/*****************************************/\
	/* now do three more radix-3 transforms: */\
	/*****************************************/\
	/* Block 1: */\
	"movq	%[__o0]	,%%rsi 									\n\t	movq	    (%%r9),%%r13 				\n\t"/* o0 */\
	"movq	%[__i0]	,%%rax									\n\t	movq	    (%%r8),%%r10				\n\t"/* i0 */\
	"movq	%[__i3]	,%%rbx									\n\t	movq	0x18(%%r8),%%r11				\n\t"/* i3 */\
	"movq	%[__i6]	,%%rcx									\n\t	movq	0x30(%%r8),%%r12				\n\t"/* i6 */\
		"vmovaps	    (%%rbx)	,%%zmm4						\n\t		vmovaps	    (%%r11)	,%%zmm12		\n\t"\
		"vmovaps	0x40(%%rbx)	,%%zmm5						\n\t		vmovaps	0x40(%%r11)	,%%zmm13		\n\t"\
		"vmovaps	    (%%rcx)	,%%zmm2						\n\t		vmovaps	    (%%r12)	,%%zmm10		\n\t"\
		"vmovaps	0x40(%%rcx)	,%%zmm3						\n\t		vmovaps	0x40(%%r12)	,%%zmm11		\n\t"\
		"vmovaps	    (%%rax)	,%%zmm0						\n\t		vmovaps	    (%%r10)	,%%zmm8 		\n\t"\
		"vmovaps	0x40(%%rax)	,%%zmm1						\n\t		vmovaps	0x40(%%r10)	,%%zmm9 		\n\t"\
		"vsubpd	%%zmm2,%%zmm4,%%zmm4						\n\t		vsubpd	%%zmm10,%%zmm12,%%zmm12		\n\t"\
		"vsubpd	%%zmm3,%%zmm5,%%zmm5						\n\t		vsubpd	%%zmm11,%%zmm13,%%zmm13		\n\t"\
	" vfmadd132pd	%%zmm31,%%zmm4,%%zmm2 					\n\t	 vfmadd132pd	%%zmm31,%%zmm12,%%zmm10	\n\t"\
	" vfmadd132pd	%%zmm31,%%zmm5,%%zmm3 					\n\t	 vfmadd132pd	%%zmm31,%%zmm13,%%zmm11	\n\t"\
		"vaddpd	%%zmm2,%%zmm0,%%zmm0						\n\t		vaddpd	%%zmm10,%%zmm8 ,%%zmm8 		\n\t"\
		"vaddpd	%%zmm3,%%zmm1,%%zmm1						\n\t		vaddpd	%%zmm11,%%zmm9 ,%%zmm9 		\n\t"\
		"vmovaps	%%zmm0		,    (%%rsi)				\n\t		vmovaps	%%zmm8 		,    (%%r13)	\n\t"\
		"vmovaps	%%zmm1		,0x40(%%rsi)				\n\t		vmovaps	%%zmm9 		,0x40(%%r13)	\n\t"\
	" vfmadd132pd	%%zmm6,%%zmm0,%%zmm2 					\n\t	 vfmadd132pd	%%zmm6,%%zmm8,%%zmm10	\n\t"\
	" vfmadd132pd	%%zmm6,%%zmm1,%%zmm3 					\n\t	 vfmadd132pd	%%zmm6,%%zmm9,%%zmm11	\n\t"\
		"vmovaps	%%zmm2,%%zmm0							\n\t		vmovaps	%%zmm10,%%zmm8 				\n\t"\
		"vmovaps	%%zmm3,%%zmm1							\n\t		vmovaps	%%zmm11,%%zmm9 				\n\t"\
	"movq	%[__o3]	,%%rdi 									\n\t	movq	0x18(%%r9),%%r14				\n\t"/* o3 */\
	"movq	%[__o6]	,%%rsi 									\n\t	movq	0x30(%%r9),%%r13				\n\t"/* o6 */\
	"vfnmadd231pd	%%zmm7,%%zmm5,%%zmm0 					\n\t	vfnmadd231pd	%%zmm7,%%zmm13,%%zmm8 	\n\t"\
	" vfmadd231pd	%%zmm7,%%zmm4,%%zmm1 					\n\t	 vfmadd231pd	%%zmm7,%%zmm12,%%zmm9 	\n\t"\
	" vfmadd231pd	%%zmm7,%%zmm5,%%zmm2 					\n\t	 vfmadd231pd	%%zmm7,%%zmm13,%%zmm10	\n\t"\
	"vfnmadd231pd	%%zmm7,%%zmm4,%%zmm3 					\n\t	vfnmadd231pd	%%zmm7,%%zmm12,%%zmm11	\n\t"\
		"vmovaps	%%zmm0		,    (%%rsi)				\n\t		vmovaps	%%zmm8 		,    (%%r13)	\n\t"\
		"vmovaps	%%zmm1		,0x40(%%rsi)				\n\t		vmovaps	%%zmm9 		,0x40(%%r13)	\n\t"\
		"vmovaps	%%zmm2		,    (%%rdi)				\n\t		vmovaps	%%zmm10		,    (%%r14)	\n\t"\
		"vmovaps	%%zmm3		,0x40(%%rdi)				\n\t		vmovaps	%%zmm11		,0x40(%%r14)	\n\t"\
	/* Block 2: */\
	"movq	%[__i1]	,%%rax									\n\t	movq	0x08(%%r8),%%r10				\n\t"/* i1 */\
	"movq	%[__i4]	,%%rbx									\n\t	movq	0x20(%%r8),%%r11				\n\t"/* i4 */\
	"movq	%[__i7]	,%%rcx									\n\t	movq	0x38(%%r8),%%r12				\n\t"/* i7 */\
		"vmovaps	    (%%rbx)	,%%zmm2						\n\t		vmovaps	    (%%r11)	,%%zmm10		\n\t"\
		"vmovaps	0x40(%%rbx)	,%%zmm3						\n\t		vmovaps	0x40(%%r11)	,%%zmm11		\n\t"\
		"vmovaps	    (%%rdx)	,%%zmm6						\n\t"/* c1 */\
		"vmovaps	0x40(%%rdx)	,%%zmm7						\n\t"\
		"vmovaps	%%zmm2		,%%zmm0						\n\t		vmovaps	%%zmm10		,%%zmm8 		\n\t"\
		"vmovaps	%%zmm3		,%%zmm1						\n\t		vmovaps	%%zmm11		,%%zmm9 		\n\t"\
		"vmulpd		%%zmm6,%%zmm2,%%zmm2					\n\t		vmulpd		%%zmm6,%%zmm10,%%zmm10	\n\t"\
		"vmulpd		%%zmm6,%%zmm3,%%zmm3					\n\t		vmulpd		%%zmm6,%%zmm11,%%zmm11	\n\t"\
	" vfmadd231pd	%%zmm7,%%zmm1,%%zmm2 					\n\t	 vfmadd231pd	%%zmm7,%%zmm9 ,%%zmm10	\n\t"\
	"vfnmadd231pd	%%zmm7,%%zmm0,%%zmm3 					\n\t	vfnmadd231pd	%%zmm7,%%zmm8 ,%%zmm11	\n\t"\
		"vmovaps	    (%%rcx)	,%%zmm4						\n\t		vmovaps	    (%%r12)	,%%zmm12		\n\t"\
		"vmovaps	0x40(%%rcx)	,%%zmm5						\n\t		vmovaps	0x40(%%r12)	,%%zmm13		\n\t"\
		"vmovaps	0x80(%%rdx)	,%%zmm6						\n\t"/* c2 */\
		"vmovaps	0xc0(%%rdx)	,%%zmm7						\n\t"\
		"vmovaps	%%zmm4		,%%zmm0						\n\t		vmovaps	%%zmm12		,%%zmm8 		\n\t"\
		"vmovaps	%%zmm5		,%%zmm1						\n\t		vmovaps	%%zmm13		,%%zmm9 		\n\t"\
		"vmulpd		%%zmm6,%%zmm4,%%zmm4					\n\t		vmulpd		%%zmm6,%%zmm12,%%zmm12	\n\t"\
		"vmulpd		%%zmm6,%%zmm5,%%zmm5					\n\t		vmulpd		%%zmm6,%%zmm13,%%zmm13	\n\t"\
	" vfmadd231pd	%%zmm7,%%zmm1,%%zmm4 					\n\t	 vfmadd231pd	%%zmm7,%%zmm9 ,%%zmm12	\n\t"\
	"vfnmadd231pd	%%zmm7,%%zmm0,%%zmm5 					\n\t	vfnmadd231pd	%%zmm7,%%zmm8 ,%%zmm13	\n\t"\
		"vmovaps	0x100(%%rdx)	,%%zmm6						\n\t"/* c3m1 */\
		"vmovaps	0x140(%%rdx)	,%%zmm7						\n\t"\
		"vmovaps	    (%%rax)	,%%zmm0						\n\t		vmovaps	    (%%r10)	,%%zmm8 		\n\t"\
		"vmovaps	0x40(%%rax)	,%%zmm1						\n\t		vmovaps	0x40(%%r10)	,%%zmm9 		\n\t"\
	"movq	%[__o7]	,%%rsi 									\n\t	movq	0x38(%%r9),%%r13				\n\t"/* o7 */\
		"vsubpd	%%zmm4,%%zmm2,%%zmm2						\n\t		vsubpd	%%zmm12,%%zmm10,%%zmm10		\n\t"\
		"vsubpd	%%zmm5,%%zmm3,%%zmm3						\n\t		vsubpd	%%zmm13,%%zmm11,%%zmm11		\n\t"\
	" vfmadd132pd	%%zmm31,%%zmm2,%%zmm4 					\n\t	 vfmadd132pd	%%zmm31,%%zmm10,%%zmm12	\n\t"\
	" vfmadd132pd	%%zmm31,%%zmm3,%%zmm5 					\n\t	 vfmadd132pd	%%zmm31,%%zmm11,%%zmm13	\n\t"\
		"vaddpd	%%zmm4,%%zmm0,%%zmm0						\n\t		vaddpd	%%zmm12,%%zmm8 ,%%zmm8 		\n\t"\
		"vaddpd	%%zmm5,%%zmm1,%%zmm1						\n\t		vaddpd	%%zmm13,%%zmm9 ,%%zmm9 		\n\t"\
		"vmovaps	%%zmm0		,    (%%rsi)				\n\t		vmovaps	%%zmm8 		,    (%%r13)	\n\t"\
		"vmovaps	%%zmm1		,0x40(%%rsi)				\n\t		vmovaps	%%zmm9 		,0x40(%%r13)	\n\t"\
	" vfmadd132pd	%%zmm6,%%zmm0,%%zmm4 					\n\t	 vfmadd132pd	%%zmm6,%%zmm8,%%zmm12	\n\t"\
	" vfmadd132pd	%%zmm6,%%zmm1,%%zmm5 					\n\t	 vfmadd132pd	%%zmm6,%%zmm9,%%zmm13	\n\t"\
		"vmovaps	%%zmm4,%%zmm0							\n\t		vmovaps	%%zmm12,%%zmm8 				\n\t"\
		"vmovaps	%%zmm5,%%zmm1							\n\t		vmovaps	%%zmm13,%%zmm9 				\n\t"\
	"movq	%[__o1]	,%%rdi 									\n\t	movq	0x08(%%r9),%%r14				\n\t"/* o1 */\
	"movq	%[__o4]	,%%rsi 									\n\t	movq	0x20(%%r9),%%r13				\n\t"/* o4 */\
	"vfnmadd231pd	%%zmm7,%%zmm3,%%zmm0 					\n\t	vfnmadd231pd	%%zmm7,%%zmm11,%%zmm8 	\n\t"\
	" vfmadd231pd	%%zmm7,%%zmm2,%%zmm1 					\n\t	 vfmadd231pd	%%zmm7,%%zmm10,%%zmm9 	\n\t"\
	" vfmadd231pd	%%zmm7,%%zmm3,%%zmm4 					\n\t	 vfmadd231pd	%%zmm7,%%zmm11,%%zmm12	\n\t"\
	"vfnmadd231pd	%%zmm7,%%zmm2,%%zmm5 					\n\t	vfnmadd231pd	%%zmm7,%%zmm10,%%zmm13	\n\t"\
		"vmovaps	%%zmm0		,    (%%rsi)				\n\t		vmovaps	%%zmm8 		,    (%%r13)	\n\t"\
		"vmovaps	%%zmm1		,0x40(%%rsi)				\n\t		vmovaps	%%zmm9 		,0x40(%%r13)	\n\t"\
		"vmovaps	%%zmm4		,    (%%rdi)				\n\t		vmovaps	%%zmm12		,    (%%r14)	\n\t"\
		"vmovaps	%%zmm5		,0x40(%%rdi)				\n\t		vmovaps	%%zmm13		,0x40(%%r14)	\n\t"\
	/* Block 3: */\
	"movq	%[__i2]	,%%rax									\n\t	movq	0x10(%%r8),%%r10				\n\t"/* i2 */\
	"movq	%[__i5]	,%%rbx									\n\t	movq	0x28(%%r8),%%r11				\n\t"/* i5 */\
	"movq	%[__i8]	,%%rcx									\n\t	movq	0x40(%%r8),%%r12				\n\t"/* i8 */\
		"vmovaps	    (%%rbx)	,%%zmm2						\n\t		vmovaps	    (%%r11)	,%%zmm10		\n\t"\
		"vmovaps	0x40(%%rbx)	,%%zmm3						\n\t		vmovaps	0x40(%%r11)	,%%zmm11		\n\t"\
		"vmovaps	0x80(%%rdx)	,%%zmm6						\n\t"/* c2 */\
		"vmovaps	0xc0(%%rdx)	,%%zmm7						\n\t"\
		"vmovaps	%%zmm2		,%%zmm0						\n\t		vmovaps	%%zmm10		,%%zmm8 		\n\t"\
		"vmovaps	%%zmm3		,%%zmm1						\n\t		vmovaps	%%zmm11		,%%zmm9 		\n\t"\
		"vmulpd		%%zmm6,%%zmm2,%%zmm2					\n\t		vmulpd		%%zmm6,%%zmm10,%%zmm10	\n\t"\
		"vmulpd		%%zmm6,%%zmm3,%%zmm3					\n\t		vmulpd		%%zmm6,%%zmm11,%%zmm11	\n\t"\
	" vfmadd231pd	%%zmm7,%%zmm1,%%zmm2 					\n\t	 vfmadd231pd	%%zmm7,%%zmm9 ,%%zmm10	\n\t"\
	"vfnmadd231pd	%%zmm7,%%zmm0,%%zmm3 					\n\t	vfnmadd231pd	%%zmm7,%%zmm8 ,%%zmm11	\n\t"\
		"vmovaps	    (%%rcx)	,%%zmm4						\n\t		vmovaps	    (%%r12)	,%%zmm12		\n\t"\
		"vmovaps	0x40(%%rcx)	,%%zmm5						\n\t		vmovaps	0x40(%%r12)	,%%zmm13		\n\t"\
		"vmovaps	0x180(%%rdx)	,%%zmm6						\n\t"/* c4 */\
		"vmovaps	0x1c0(%%rdx)	,%%zmm7						\n\t"\
		"vmovaps	%%zmm4		,%%zmm0						\n\t		vmovaps	%%zmm12		,%%zmm8 		\n\t"\
		"vmovaps	%%zmm5		,%%zmm1						\n\t		vmovaps	%%zmm13		,%%zmm9 		\n\t"\
		"vmulpd		%%zmm6,%%zmm4,%%zmm4					\n\t		vmulpd		%%zmm6,%%zmm12,%%zmm12	\n\t"\
		"vmulpd		%%zmm6,%%zmm5,%%zmm5					\n\t		vmulpd		%%zmm6,%%zmm13,%%zmm13	\n\t"\
	" vfmadd231pd	%%zmm7,%%zmm1,%%zmm4 					\n\t	 vfmadd231pd	%%zmm7,%%zmm9 ,%%zmm12	\n\t"\
	"vfnmadd231pd	%%zmm7,%%zmm0,%%zmm5 					\n\t	vfnmadd231pd	%%zmm7,%%zmm8 ,%%zmm13	\n\t"\
		"vmovaps	0x100(%%rdx)	,%%zmm6						\n\t"/* c3m1 */\
		"vmovaps	0x140(%%rdx)	,%%zmm7						\n\t"\
		"vmovaps	    (%%rax)	,%%zmm0						\n\t		vmovaps	    (%%r10)	,%%zmm8 		\n\t"\
		"vmovaps	0x40(%%rax)	,%%zmm1						\n\t		vmovaps	0x40(%%r10)	,%%zmm9 		\n\t"\
	"movq	%[__o5]	,%%rsi 									\n\t	movq	0x28(%%r9),%%r13				\n\t"/* o5 */\
		"vsubpd	%%zmm4,%%zmm2,%%zmm2						\n\t		vsubpd	%%zmm12,%%zmm10,%%zmm10		\n\t"\
		"vsubpd	%%zmm5,%%zmm3,%%zmm3						\n\t		vsubpd	%%zmm13,%%zmm11,%%zmm11		\n\t"\
	" vfmadd132pd	%%zmm31,%%zmm2,%%zmm4 					\n\t	 vfmadd132pd	%%zmm31,%%zmm10,%%zmm12	\n\t"\
	" vfmadd132pd	%%zmm31,%%zmm3,%%zmm5 					\n\t	 vfmadd132pd	%%zmm31,%%zmm11,%%zmm13	\n\t"\
		"vaddpd	%%zmm4,%%zmm0,%%zmm0						\n\t		vaddpd	%%zmm12,%%zmm8 ,%%zmm8 		\n\t"\
		"vaddpd	%%zmm5,%%zmm1,%%zmm1						\n\t		vaddpd	%%zmm13,%%zmm9 ,%%zmm9 		\n\t"\
		"vmovaps	%%zmm0		,    (%%rsi)				\n\t		vmovaps	%%zmm8 		,    (%%r13)	\n\t"\
		"vmovaps	%%zmm1		,0x40(%%rsi)				\n\t		vmovaps	%%zmm9 		,0x40(%%r13)	\n\t"\
	" vfmadd132pd	%%zmm6,%%zmm0,%%zmm4 					\n\t	 vfmadd132pd	%%zmm6,%%zmm8,%%zmm12	\n\t"\
	" vfmadd132pd	%%zmm6,%%zmm1,%%zmm5 					\n\t	 vfmadd132pd	%%zmm6,%%zmm9,%%zmm13	\n\t"\
		"vmovaps	%%zmm4,%%zmm0							\n\t		vmovaps	%%zmm12,%%zmm8 				\n\t"\
		"vmovaps	%%zmm5,%%zmm1							\n\t		vmovaps	%%zmm13,%%zmm9 				\n\t"\
	"movq	%[__o8]	,%%rdi 									\n\t	movq	0x40(%%r9),%%r14				\n\t"/* o8 */\
	"movq	%[__o2]	,%%rsi 									\n\t	movq	0x10(%%r9),%%r13				\n\t"/* o2 */\
	"vfnmadd231pd	%%zmm7,%%zmm3,%%zmm0 					\n\t	vfnmadd231pd	%%zmm7,%%zmm11,%%zmm8 	\n\t"\
	" vfmadd231pd	%%zmm7,%%zmm2,%%zmm1 					\n\t	 vfmadd231pd	%%zmm7,%%zmm10,%%zmm9 	\n\t"\
	" vfmadd231pd	%%zmm7,%%zmm3,%%zmm4 					\n\t	 vfmadd231pd	%%zmm7,%%zmm11,%%zmm12	\n\t"\
	"vfnmadd231pd	%%zmm7,%%zmm2,%%zmm5 					\n\t	vfnmadd231pd	%%zmm7,%%zmm10,%%zmm13	\n\t"\
		"vmovaps	%%zmm0		,    (%%rsi)				\n\t		vmovaps	%%zmm8 		,    (%%r13)	\n\t"\
		"vmovaps	%%zmm1		,0x40(%%rsi)				\n\t		vmovaps	%%zmm9 		,0x40(%%r13)	\n\t"\
		"vmovaps	%%zmm4		,    (%%rdi)				\n\t		vmovaps	%%zmm12		,    (%%r14)	\n\t"\
		"vmovaps	%%zmm5		,0x40(%%rdi)				\n\t		vmovaps	%%zmm13		,0x40(%%r14)	\n\t"\
		:					/* outputs: none */\
		: [__i0] "m" (Xi0)	/* All inputs from memory addresses here */\
		 ,[__i1] "m" (Xi1)\
		 ,[__i2] "m" (Xi2)\
		 ,[__i3] "m" (Xi3)\
		 ,[__i4] "m" (Xi4)\
		 ,[__i5] "m" (Xi5)\
		 ,[__i6] "m" (Xi6)\
		 ,[__i7] "m" (Xi7)\
		 ,[__i8] "m" (Xi8)\
		 ,[__c1] "m" (Xc1)\
		 ,[__two] "m" (Xtwo)\
		 ,[__o0] "m" (Xo0)\
		 ,[__o1] "m" (Xo1)\
		 ,[__o2] "m" (Xo2)\
		 ,[__o3] "m" (Xo3)\
		 ,[__o4] "m" (Xo4)\
		 ,[__o5] "m" (Xo5)\
		 ,[__o6] "m" (Xo6)\
		 ,[__o7] "m" (Xo7)\
		 ,[__o8] "m" (Xo8)\
		 ,[__iptr] "m" (Xiptr)	/* Due to GCC macro argc limit of 30, to enable 16-register data-doubled version */\
		 ,[__optr] "m" (Xoptr)	/* of the radix-9 macros use 2 length-9 ptr arrays for second set of IOs. */\
		: "cc","memory","rax","rbx","rcx","rdx","rsi","rdi","r8","r9","r10","r11","r12","r13","r14","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15","xmm31"	/* Clobbered registers */\
		);\
	}

	/*...Radix-9 DIF: Ins in memory locations __i0-8.\
		Outs at memory locations __o0-8, assumed disjoint with inputs:\
	*/\
	#define SSE2_RADIX_09_DIF(Xi0,Xi1,Xi2,Xi3,Xi4,Xi5,Xi6,Xi7,Xi8, Xc1, Xo0,Xo1,Xo2,Xo3,Xo4,Xo5,Xo6,Xo7,Xo8)\
	{\
	__asm__ volatile (\
	"movq	%[__i0]	,%%rax 	\n\t"/* __i0-8; e[abc]x store input addresses */\
	"movq	%[__o0]	,%%rsi	\n\t"\
	"movq	%[__c1]	,%%rdx 	\n\t"/* edx stores trig addresses throughout */\
	/* Block 1: */\
		"movq	%%rax		,%%rbx\n\t"\
		"movq	%%rax		,%%rcx\n\t"\
		"movq	%[__i3]	,%%rbx 	\n\t"\
		"movq	%[__i6]	,%%rcx 	\n\t"\
		"vmovaps	    (%%rbx)	,%%zmm2\n\t"\
		"vmovaps	0x40(%%rbx)	,%%zmm3\n\t"\
		"vmovaps	    (%%rax)	,%%zmm0\n\t"\
		"vmovaps	0x40(%%rax)	,%%zmm1\n\t"\
		"vmovaps	    (%%rcx)	,%%zmm6\n\t"\
		"vmovaps	0x40(%%rcx)	,%%zmm7\n\t"\
		"vmovaps	%%zmm2		,%%zmm4\n\t"\
		"vmovaps	%%zmm3		,%%zmm5\n\t"\
		"vaddpd	%%zmm6,%%zmm2,%%zmm2\n\t"\
		"vaddpd	%%zmm7,%%zmm3,%%zmm3\n\t"\
		"vsubpd	%%zmm6,%%zmm4,%%zmm4\n\t"\
		"vsubpd	%%zmm7,%%zmm5,%%zmm5\n\t"\
		"vaddpd	%%zmm2,%%zmm0,%%zmm0\n\t"\
		"vaddpd	%%zmm3,%%zmm1,%%zmm1\n\t"\
		"vmovaps	0x100(%%rdx)	,%%zmm6\n\t"/* c3m1 */\
		"vmovaps	0x140(%%rdx)	,%%zmm7\n\t"\
		"vmovaps	%%zmm0		,    (%%rsi)	/* <- t00 */\n\t"\
		"vmovaps	%%zmm1		,0x40(%%rsi)	/* <- t01 */\n\t"\
	"movq	%[__o1]	,%%rdi	\n\t"\
	"movq	%[__o2]	,%%rsi	\n\t"\
		"vmulpd	%%zmm6,%%zmm2,%%zmm2\n\t"\
		"vmulpd	%%zmm6,%%zmm3,%%zmm3\n\t"\
		"vmulpd	%%zmm7,%%zmm4,%%zmm4\n\t"\
		"vmulpd	%%zmm7,%%zmm5,%%zmm5\n\t"\
		"vaddpd	%%zmm0,%%zmm2,%%zmm2\n\t"\
		"vaddpd	%%zmm1,%%zmm3,%%zmm3\n\t"\
		"vmovaps	%%zmm2		,%%zmm0\n\t"\
		"vmovaps	%%zmm3		,%%zmm1\n\t"\
		"vsubpd	%%zmm5,%%zmm2,%%zmm2\n\t"\
		"vaddpd	%%zmm4,%%zmm3,%%zmm3\n\t"\
		"vaddpd	%%zmm5,%%zmm0,%%zmm0\n\t"\
		"vsubpd	%%zmm4,%%zmm1,%%zmm1\n\t"\
		"vmovaps	%%zmm2		,    (%%rdi)	/* <- t02 */\n\t"\
		"vmovaps	%%zmm3		,0x40(%%rdi)	/* <- t03 */\n\t"\
		"vmovaps	%%zmm0		,    (%%rsi)	/* <- t04 */\n\t"\
		"vmovaps	%%zmm1		,0x40(%%rsi)	/* <- t05 */\n\t"\
	/* Block 2: */\
		"movq	%[__i1]	,%%rax 	\n\t"\
		"movq	%[__i4]	,%%rbx 	\n\t"\
		"movq	%[__i7]	,%%rcx 	\n\t"\
		"vmovaps	    (%%rbx)	,%%zmm4\n\t"\
		"vmovaps	0x40(%%rbx)	,%%zmm5\n\t"\
		"vmovaps	    (%%rax)	,%%zmm0\n\t"\
		"vmovaps	0x40(%%rax)	,%%zmm1\n\t"\
		"vmovaps	    (%%rcx)	,%%zmm2\n\t"\
		"vmovaps	0x40(%%rcx)	,%%zmm3\n\t"\
	"movq	%[__o3]	,%%rdi	\n\t"\
		"vsubpd	%%zmm2,%%zmm4,%%zmm4\n\t"\
		"vsubpd	%%zmm3,%%zmm5,%%zmm5\n\t"\
		"vaddpd	%%zmm2,%%zmm2,%%zmm2\n\t"\
		"vaddpd	%%zmm3,%%zmm3,%%zmm3\n\t"\
		"vaddpd	%%zmm4,%%zmm2,%%zmm2\n\t"\
		"vaddpd	%%zmm5,%%zmm3,%%zmm3\n\t"\
		"vaddpd	%%zmm2,%%zmm0,%%zmm0\n\t"\
		"vaddpd	%%zmm3,%%zmm1,%%zmm1\n\t"\
		"vmovaps	%%zmm0		,    (%%rdi)	/* <- t06 */\n\t"\
		"vmovaps	%%zmm1		,0x40(%%rdi)	/* <- t07 */\n\t"\
	"movq	%[__o4]	,%%rsi	\n\t"\
	"movq	%[__o5]	,%%rdi	\n\t"\
		"vmulpd	%%zmm6,%%zmm2,%%zmm2\n\t"\
		"vmulpd	%%zmm6,%%zmm3,%%zmm3\n\t"\
		"vmulpd	%%zmm7,%%zmm4,%%zmm4\n\t"\
		"vmulpd	%%zmm7,%%zmm5,%%zmm5\n\t"\
		"vaddpd	%%zmm0,%%zmm2,%%zmm2\n\t"\
		"vaddpd	%%zmm1,%%zmm3,%%zmm3\n\t"\
		"vmovaps	%%zmm2		,%%zmm0\n\t"\
		"vmovaps	%%zmm3		,%%zmm1\n\t"\
		"vsubpd	%%zmm5,%%zmm2,%%zmm2\n\t"\
		"vaddpd	%%zmm4,%%zmm3,%%zmm3\n\t"\
		"vaddpd	%%zmm5,%%zmm0,%%zmm0\n\t"\
		"vsubpd	%%zmm4,%%zmm1,%%zmm1\n\t"\
		"vmovaps	%%zmm2		,    (%%rsi)	/* <- t08 */\n\t"\
		"vmovaps	%%zmm3		,0x40(%%rsi)	/* <- t09 */\n\t"\
		"vmovaps	%%zmm0		,    (%%rdi)	/* <- t0a */\n\t"\
		"vmovaps	%%zmm1		,0x40(%%rdi)	/* <- t0b */\n\t"\
	/* Block 3: */\
		"movq	%[__i2]	,%%rax 	\n\t"\
		"movq	%[__i5]	,%%rbx 	\n\t"\
		"movq	%[__i8]	,%%rcx 	\n\t"\
		"vmovaps	    (%%rbx)	,%%zmm4\n\t"\
		"vmovaps	0x40(%%rbx)	,%%zmm5\n\t"\
		"vmovaps	    (%%rax)	,%%zmm0\n\t"\
		"vmovaps	0x40(%%rax)	,%%zmm1\n\t"\
		"vmovaps	    (%%rcx)	,%%zmm2\n\t"\
		"vmovaps	0x40(%%rcx)	,%%zmm3\n\t"\
	"movq	%[__o6]	,%%rsi	\n\t"\
		"vsubpd	%%zmm2,%%zmm4,%%zmm4\n\t"\
		"vsubpd	%%zmm3,%%zmm5,%%zmm5\n\t"\
		"vaddpd	%%zmm2,%%zmm2,%%zmm2\n\t"\
		"vaddpd	%%zmm3,%%zmm3,%%zmm3\n\t"\
		"vaddpd	%%zmm4,%%zmm2,%%zmm2\n\t"\
		"vaddpd	%%zmm5,%%zmm3,%%zmm3\n\t"\
		"vaddpd	%%zmm2,%%zmm0,%%zmm0\n\t"\
		"vaddpd	%%zmm3,%%zmm1,%%zmm1\n\t"\
		"vmovaps	%%zmm0		,    (%%rsi)	/* <- t0c */\n\t"\
		"vmovaps	%%zmm1		,0x40(%%rsi)	/* <- t0d */\n\t"\
	"movq	%[__o7]	,%%rdi	\n\t"\
	"movq	%[__o8]	,%%rsi	\n\t"\
		"vmulpd	%%zmm6,%%zmm2,%%zmm2\n\t"\
		"vmulpd	%%zmm6,%%zmm3,%%zmm3\n\t"\
		"vmulpd	%%zmm7,%%zmm4,%%zmm4\n\t"\
		"vmulpd	%%zmm7,%%zmm5,%%zmm5\n\t"\
		"vaddpd	%%zmm0,%%zmm2,%%zmm2\n\t"\
		"vaddpd	%%zmm1,%%zmm3,%%zmm3\n\t"\
		"vmovaps	%%zmm2		,%%zmm0\n\t"\
		"vmovaps	%%zmm3		,%%zmm1\n\t"\
		"vsubpd	%%zmm5,%%zmm2,%%zmm2\n\t"\
		"vaddpd	%%zmm4,%%zmm3,%%zmm3\n\t"\
		"vaddpd	%%zmm5,%%zmm0,%%zmm0\n\t"\
		"vsubpd	%%zmm4,%%zmm1,%%zmm1\n\t"\
		"vmovaps	%%zmm2		,    (%%rdi)	/* <- t0e */\n\t"\
		"vmovaps	%%zmm3		,0x40(%%rdi)	/* <- t0f */\n\t"\
		"vmovaps	%%zmm0		,    (%%rsi)	/* <- t0g */\n\t"\
		"vmovaps	%%zmm1		,0x40(%%rsi)	/* <- t0h */\n\t"\
	/*****************************************/\
	/* now do three more radix-3 transforms: */\
	/*****************************************/\
	/* Block 1: */\
	"movq	%[__o0]	,%%rax	\n\t"\
	"movq	%[__o3]	,%%rbx	\n\t"\
	"movq	%[__o6]	,%%rcx	\n\t"\
		"vmovaps	    (%%rbx)	,%%zmm4\n\t"\
		"vmovaps	0x40(%%rbx)	,%%zmm5\n\t"\
		"vmovaps	    (%%rcx)	,%%zmm2\n\t"\
		"vmovaps	0x40(%%rcx)	,%%zmm3\n\t"\
		"vmovaps	    (%%rax)	,%%zmm0\n\t"\
		"vmovaps	0x40(%%rax)	,%%zmm1\n\t"\
		"vsubpd	%%zmm2,%%zmm4,%%zmm4\n\t"\
		"vsubpd	%%zmm3,%%zmm5,%%zmm5\n\t"\
		"vaddpd	%%zmm2,%%zmm2,%%zmm2\n\t"\
		"vaddpd	%%zmm3,%%zmm3,%%zmm3\n\t"\
		"vaddpd	%%zmm4,%%zmm2,%%zmm2\n\t"\
		"vaddpd	%%zmm5,%%zmm3,%%zmm3\n\t"\
		"vaddpd	%%zmm2,%%zmm0,%%zmm0\n\t"\
		"vaddpd	%%zmm3,%%zmm1,%%zmm1\n\t"\
		"vmovaps	%%zmm0		,    (%%rax)\n\t"\
		"vmovaps	%%zmm1		,0x40(%%rax)\n\t"\
		"vmulpd	%%zmm6,%%zmm2,%%zmm2\n\t"\
		"vmulpd	%%zmm6,%%zmm3,%%zmm3\n\t"\
		"vmulpd	%%zmm7,%%zmm4,%%zmm4\n\t"\
		"vmulpd	%%zmm7,%%zmm5,%%zmm5\n\t"\
		"vaddpd	%%zmm0,%%zmm2,%%zmm2\n\t"\
		"vaddpd	%%zmm1,%%zmm3,%%zmm3\n\t"\
		"vmovaps	%%zmm2		,%%zmm0\n\t"\
		"vmovaps	%%zmm3		,%%zmm1\n\t"\
		"vsubpd	%%zmm5,%%zmm2,%%zmm2\n\t"\
		"vaddpd	%%zmm4,%%zmm3,%%zmm3\n\t"\
		"vaddpd	%%zmm5,%%zmm0,%%zmm0\n\t"\
		"vsubpd	%%zmm4,%%zmm1,%%zmm1\n\t"\
		"vmovaps	%%zmm2		,    (%%rbx)\n\t"\
		"vmovaps	%%zmm3		,0x40(%%rbx)\n\t"\
		"vmovaps	%%zmm0		,    (%%rcx)\n\t"\
		"vmovaps	%%zmm1		,0x40(%%rcx)\n\t"\
	/* Block 2: */\
	"movq	%[__o1]	,%%rax	\n\t"\
	"movq	%[__o4]	,%%rbx	\n\t"\
	"movq	%[__o7]	,%%rcx	\n\t"\
		"vmovaps	    (%%rbx)	,%%zmm2\n\t"\
		"vmovaps	0x40(%%rbx)	,%%zmm3\n\t"\
		"vmovaps	    (%%rdx)	,%%zmm6\n\t"/* c1 */\
		"vmovaps	0x40(%%rdx)	,%%zmm7\n\t"\
		"vmovaps	%%zmm2		,%%zmm0\n\t"\
		"vmovaps	%%zmm3		,%%zmm1\n\t"\
		"vmulpd	%%zmm6,%%zmm2,%%zmm2\n\t"\
		"vmulpd	%%zmm6,%%zmm3,%%zmm3\n\t"\
		"vmulpd	%%zmm7,%%zmm0,%%zmm0\n\t"\
		"vmulpd	%%zmm7,%%zmm1,%%zmm1\n\t"\
		"vsubpd	%%zmm1,%%zmm2,%%zmm2\n\t"\
		"vaddpd	%%zmm0,%%zmm3,%%zmm3\n\t"\
		"vmovaps	    (%%rcx)	,%%zmm4\n\t"\
		"vmovaps	0x40(%%rcx)	,%%zmm5\n\t"\
		"vmovaps	0x80(%%rdx)	,%%zmm6\n\t"/* c2 */\
		"vmovaps	0xc0(%%rdx)	,%%zmm7\n\t"\
		"vmovaps	%%zmm4		,%%zmm0\n\t"\
		"vmovaps	%%zmm5		,%%zmm1\n\t"\
		"vmulpd	%%zmm6,%%zmm4,%%zmm4\n\t"\
		"vmulpd	%%zmm6,%%zmm5,%%zmm5\n\t"\
		"vmulpd	%%zmm7,%%zmm0,%%zmm0\n\t"\
		"vmulpd	%%zmm7,%%zmm1,%%zmm1\n\t"\
		"vsubpd	%%zmm1,%%zmm4,%%zmm4\n\t"\
		"vaddpd	%%zmm0,%%zmm5,%%zmm5\n\t"\
		"vmovaps	0x100(%%rdx)	,%%zmm6\n\t"/* c3m1 */\
		"vmovaps	0x140(%%rdx)	,%%zmm7\n\t"\
		"vmovaps	    (%%rax)	,%%zmm0\n\t"\
		"vmovaps	0x40(%%rax)	,%%zmm1\n\t"\
		"vsubpd	%%zmm4,%%zmm2,%%zmm2\n\t"\
		"vsubpd	%%zmm5,%%zmm3,%%zmm3\n\t"\
		"vaddpd	%%zmm4,%%zmm4,%%zmm4\n\t"\
		"vaddpd	%%zmm5,%%zmm5,%%zmm5\n\t"\
		"vaddpd	%%zmm2,%%zmm4,%%zmm4\n\t"\
		"vaddpd	%%zmm3,%%zmm5,%%zmm5\n\t"\
		"vaddpd	%%zmm4,%%zmm0,%%zmm0\n\t"\
		"vaddpd	%%zmm5,%%zmm1,%%zmm1\n\t"\
		"vmovaps	%%zmm0		,    (%%rax)\n\t"\
		"vmovaps	%%zmm1		,0x40(%%rax)\n\t"\
		"vmulpd	%%zmm6,%%zmm4,%%zmm4\n\t"\
		"vmulpd	%%zmm6,%%zmm5,%%zmm5\n\t"\
		"vmulpd	%%zmm7,%%zmm2,%%zmm2\n\t"\
		"vmulpd	%%zmm7,%%zmm3,%%zmm3\n\t"\
		"vaddpd	%%zmm0,%%zmm4,%%zmm4\n\t"\
		"vaddpd	%%zmm1,%%zmm5,%%zmm5\n\t"\
		"vmovaps	%%zmm4		,%%zmm0\n\t"\
		"vmovaps	%%zmm5		,%%zmm1\n\t"\
		"vsubpd	%%zmm3,%%zmm4,%%zmm4\n\t"\
		"vaddpd	%%zmm2,%%zmm5,%%zmm5\n\t"\
		"vaddpd	%%zmm3,%%zmm0,%%zmm0\n\t"\
		"vsubpd	%%zmm2,%%zmm1,%%zmm1\n\t"\
		"vmovaps	%%zmm4		,    (%%rbx)\n\t"\
		"vmovaps	%%zmm5		,0x40(%%rbx)\n\t"\
		"vmovaps	%%zmm0		,    (%%rcx)\n\t"\
		"vmovaps	%%zmm1		,0x40(%%rcx)\n\t"\
	/* Block 3: */\
	"movq	%[__o2]	,%%rax	\n\t"\
	"movq	%[__o5]	,%%rbx	\n\t"\
	"movq	%[__o8]	,%%rcx	\n\t"\
		"vmovaps	    (%%rbx)	,%%zmm2\n\t"\
		"vmovaps	0x40(%%rbx)	,%%zmm3\n\t"\
		"vmovaps	0x80(%%rdx)	,%%zmm6\n\t"/* c2 */\
		"vmovaps	0xc0(%%rdx)	,%%zmm7\n\t"\
		"vmovaps	%%zmm2		,%%zmm0\n\t"\
		"vmovaps	%%zmm3		,%%zmm1\n\t"\
		"vmulpd	%%zmm6,%%zmm2,%%zmm2\n\t"\
		"vmulpd	%%zmm6,%%zmm3,%%zmm3\n\t"\
		"vmulpd	%%zmm7,%%zmm0,%%zmm0\n\t"\
		"vmulpd	%%zmm7,%%zmm1,%%zmm1\n\t"\
		"vsubpd	%%zmm1,%%zmm2,%%zmm2\n\t"\
		"vaddpd	%%zmm0,%%zmm3,%%zmm3\n\t"\
		"vmovaps	    (%%rcx)	,%%zmm4\n\t"\
		"vmovaps	0x40(%%rcx)	,%%zmm5\n\t"\
		"vmovaps	0x180(%%rdx)	,%%zmm6\n\t"/* c4 */\
		"vmovaps	0x1c0(%%rdx)	,%%zmm7\n\t"\
		"vmovaps	%%zmm4		,%%zmm0\n\t"\
		"vmovaps	%%zmm5		,%%zmm1\n\t"\
		"vmulpd	%%zmm6,%%zmm4,%%zmm4\n\t"\
		"vmulpd	%%zmm6,%%zmm5,%%zmm5\n\t"\
		"vmulpd	%%zmm7,%%zmm0,%%zmm0\n\t"\
		"vmulpd	%%zmm7,%%zmm1,%%zmm1\n\t"\
		"vsubpd	%%zmm1,%%zmm4,%%zmm4\n\t"\
		"vaddpd	%%zmm0,%%zmm5,%%zmm5\n\t"\
		"vmovaps	0x100(%%rdx)	,%%zmm6\n\t"/* c3m1 */\
		"vmovaps	0x140(%%rdx)	,%%zmm7\n\t"\
		"vmovaps	    (%%rax)	,%%zmm0\n\t"\
		"vmovaps	0x40(%%rax)	,%%zmm1\n\t"\
		"\n\t"\
		"vsubpd	%%zmm4,%%zmm2,%%zmm2\n\t"\
		"vsubpd	%%zmm5,%%zmm3,%%zmm3\n\t"\
		"vaddpd	%%zmm4,%%zmm4,%%zmm4\n\t"\
		"vaddpd	%%zmm5,%%zmm5,%%zmm5\n\t"\
		"vaddpd	%%zmm2,%%zmm4,%%zmm4\n\t"\
		"vaddpd	%%zmm3,%%zmm5,%%zmm5\n\t"\
		"vaddpd	%%zmm4,%%zmm0,%%zmm0\n\t"\
		"vaddpd	%%zmm5,%%zmm1,%%zmm1\n\t"\
		"vmovaps	%%zmm0		,    (%%rax)\n\t"\
		"vmovaps	%%zmm1		,0x40(%%rax)\n\t"\
		"vmulpd	%%zmm6,%%zmm4,%%zmm4\n\t"\
		"vmulpd	%%zmm6,%%zmm5,%%zmm5\n\t"\
		"vmulpd	%%zmm7,%%zmm2,%%zmm2\n\t"\
		"vmulpd	%%zmm7,%%zmm3,%%zmm3\n\t"\
		"vaddpd	%%zmm0,%%zmm4,%%zmm4\n\t"\
		"vaddpd	%%zmm1,%%zmm5,%%zmm5\n\t"\
		"\n\t"\
		"vmovaps	%%zmm4		,%%zmm0\n\t"\
		"vmovaps	%%zmm5		,%%zmm1\n\t"\
		"vsubpd	%%zmm3,%%zmm4,%%zmm4\n\t"\
		"vaddpd	%%zmm2,%%zmm5,%%zmm5\n\t"\
		"vaddpd	%%zmm3,%%zmm0,%%zmm0\n\t"\
		"vsubpd	%%zmm2,%%zmm1,%%zmm1\n\t"\
		"vmovaps	%%zmm4		,    (%%rbx)\n\t"\
		"vmovaps	%%zmm5		,0x40(%%rbx)\n\t"\
		"vmovaps	%%zmm0		,    (%%rcx)\n\t"\
		"vmovaps	%%zmm1		,0x40(%%rcx)\n\t"\
		:					/* outputs: none */\
		: [__i0] "m" (Xi0)	/* All inputs from memory addresses here */\
		 ,[__i1] "m" (Xi1)\
		 ,[__i2] "m" (Xi2)\
		 ,[__i3] "m" (Xi3)\
		 ,[__i4] "m" (Xi4)\
		 ,[__i5] "m" (Xi5)\
		 ,[__i6] "m" (Xi6)\
		 ,[__i7] "m" (Xi7)\
		 ,[__i8] "m" (Xi8)\
		 ,[__c1] "m" (Xc1)\
		 ,[__o0] "m" (Xo0)\
		 ,[__o1] "m" (Xo1)\
		 ,[__o2] "m" (Xo2)\
		 ,[__o3] "m" (Xo3)\
		 ,[__o4] "m" (Xo4)\
		 ,[__o5] "m" (Xo5)\
		 ,[__o6] "m" (Xo6)\
		 ,[__o7] "m" (Xo7)\
		 ,[__o8] "m" (Xo8)\
		: "cc","memory","rax","rbx","rcx","rdx","rsi","rdi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
		);\
	}

	/*...Radix-9 DIT: Ins in memory locations __i0-8.\
		Outs at memory locations __o0-8, assumed disjoint with inputs:\
	*/\
	#define SSE2_RADIX_09_DIT(Xi0,Xi1,Xi2,Xi3,Xi4,Xi5,Xi6,Xi7,Xi8, Xc1, Xo0,Xo1,Xo2,Xo3,Xo4,Xo5,Xo6,Xo7,Xo8)\
	{\
	__asm__ volatile (\
	"movq	%[__i0]	,%%rax\n\t"\
	"movq	%[__i1]	,%%rbx\n\t"\
	"movq	%[__i2]	,%%rcx\n\t"\
	"movq	%[__c1]	,%%rdx 	\n\t"/* edx stores trig addresses throughout */\
	/* Block 1: */\
		"vmovaps	    (%%rbx)	,%%zmm2\n\t"\
		"vmovaps	0x40(%%rbx)	,%%zmm3\n\t"\
		"vmovaps	    (%%rax)	,%%zmm0\n\t"\
		"vmovaps	0x40(%%rax)	,%%zmm1\n\t"\
		"vmovaps	    (%%rcx)	,%%zmm6\n\t"\
		"vmovaps	0x40(%%rcx)	,%%zmm7\n\t"\
		"vmovaps	%%zmm2		,%%zmm4\n\t"\
		"vmovaps	%%zmm3		,%%zmm5\n\t"\
		"vaddpd	%%zmm6,%%zmm2,%%zmm2\n\t"\
		"vaddpd	%%zmm7,%%zmm3,%%zmm3\n\t"\
		"vsubpd	%%zmm6,%%zmm4,%%zmm4\n\t"\
		"vsubpd	%%zmm7,%%zmm5,%%zmm5\n\t"\
		"vaddpd	%%zmm2,%%zmm0,%%zmm0\n\t"\
		"vaddpd	%%zmm3,%%zmm1,%%zmm1\n\t"\
		"vmovaps	0x100(%%rdx)	,%%zmm6\n\t"/* c3m1 */\
		"vmovaps	0x140(%%rdx)	,%%zmm7\n\t"\
		"vmovaps	%%zmm0		,    (%%rax)\n\t"\
		"vmovaps	%%zmm1		,0x40(%%rax)\n\t"\
		"vmulpd	%%zmm6,%%zmm2,%%zmm2\n\t"\
		"vmulpd	%%zmm6,%%zmm3,%%zmm3\n\t"\
		"vmulpd	%%zmm7,%%zmm4,%%zmm4\n\t"\
		"vmulpd	%%zmm7,%%zmm5,%%zmm5\n\t"\
		"vaddpd	%%zmm0,%%zmm2,%%zmm2\n\t"\
		"vaddpd	%%zmm1,%%zmm3,%%zmm3\n\t"\
		"vmovaps	%%zmm2		,%%zmm0\n\t"\
		"vmovaps	%%zmm3		,%%zmm1\n\t"\
		"vaddpd	%%zmm5,%%zmm2,%%zmm2\n\t"\
		"vsubpd	%%zmm4,%%zmm3,%%zmm3\n\t"\
		"vsubpd	%%zmm5,%%zmm0,%%zmm0\n\t"\
		"vaddpd	%%zmm4,%%zmm1,%%zmm1\n\t"\
		"vmovaps	%%zmm2		,    (%%rbx)\n\t"\
		"vmovaps	%%zmm3		,0x40(%%rbx)\n\t"\
		"vmovaps	%%zmm0		,    (%%rcx)\n\t"\
		"vmovaps	%%zmm1		,0x40(%%rcx)\n\t"\
	/* Block 2: */\
	"movq	%[__i3]	,%%rax\n\t"\
	"movq	%[__i4]	,%%rbx\n\t"\
	"movq	%[__i5]	,%%rcx\n\t"\
		"vmovaps	    (%%rbx)	,%%zmm4\n\t"\
		"vmovaps	0x40(%%rbx)	,%%zmm5\n\t"\
		"vmovaps	    (%%rcx)	,%%zmm2\n\t"\
		"vmovaps	0x40(%%rcx)	,%%zmm3\n\t"\
		"vmovaps	    (%%rax)	,%%zmm0\n\t"\
		"vmovaps	0x40(%%rax)	,%%zmm1\n\t"\
		"vsubpd	%%zmm2,%%zmm4,%%zmm4\n\t"\
		"vsubpd	%%zmm3,%%zmm5,%%zmm5\n\t"\
		"vaddpd	%%zmm2,%%zmm2,%%zmm2\n\t"\
		"vaddpd	%%zmm3,%%zmm3,%%zmm3\n\t"\
		"vaddpd	%%zmm4,%%zmm2,%%zmm2\n\t"\
		"vaddpd	%%zmm5,%%zmm3,%%zmm3\n\t"\
		"vaddpd	%%zmm2,%%zmm0,%%zmm0\n\t"\
		"vaddpd	%%zmm3,%%zmm1,%%zmm1\n\t"\
		"vmovaps	%%zmm0		,    (%%rax)\n\t"\
		"vmovaps	%%zmm1		,0x40(%%rax)\n\t"\
		"vmulpd	%%zmm6,%%zmm2,%%zmm2\n\t"\
		"vmulpd	%%zmm6,%%zmm3,%%zmm3\n\t"\
		"vmulpd	%%zmm7,%%zmm4,%%zmm4\n\t"\
		"vmulpd	%%zmm7,%%zmm5,%%zmm5\n\t"\
		"vaddpd	%%zmm0,%%zmm2,%%zmm2\n\t"\
		"vaddpd	%%zmm1,%%zmm3,%%zmm3\n\t"\
		"vmovaps	%%zmm2		,%%zmm0\n\t"\
		"vmovaps	%%zmm3		,%%zmm1\n\t"\
		"vaddpd	%%zmm5,%%zmm2,%%zmm2\n\t"\
		"vsubpd	%%zmm4,%%zmm3,%%zmm3\n\t"\
		"vsubpd	%%zmm5,%%zmm0,%%zmm0\n\t"\
		"vaddpd	%%zmm4,%%zmm1,%%zmm1\n\t"\
		"vmovaps	%%zmm2		,    (%%rbx)\n\t"\
		"vmovaps	%%zmm3		,0x40(%%rbx)\n\t"\
		"vmovaps	%%zmm0		,    (%%rcx)\n\t"\
		"vmovaps	%%zmm1		,0x40(%%rcx)\n\t"\
	/* Block 3: */\
	"movq	%[__i6]	,%%rax\n\t"\
	"movq	%[__i7]	,%%rbx\n\t"\
	"movq	%[__i8]	,%%rcx\n\t"\
		"vmovaps	    (%%rbx)	,%%zmm4\n\t"\
		"vmovaps	0x40(%%rbx)	,%%zmm5\n\t"\
		"vmovaps	    (%%rcx)	,%%zmm2\n\t"\
		"vmovaps	0x40(%%rcx)	,%%zmm3\n\t"\
		"vmovaps	    (%%rax)	,%%zmm0\n\t"\
		"vmovaps	0x40(%%rax)	,%%zmm1\n\t"\
		"vsubpd	%%zmm2,%%zmm4,%%zmm4\n\t"\
		"vsubpd	%%zmm3,%%zmm5,%%zmm5\n\t"\
		"vaddpd	%%zmm2,%%zmm2,%%zmm2\n\t"\
		"vaddpd	%%zmm3,%%zmm3,%%zmm3\n\t"\
		"vaddpd	%%zmm4,%%zmm2,%%zmm2\n\t"\
		"vaddpd	%%zmm5,%%zmm3,%%zmm3\n\t"\
		"vaddpd	%%zmm2,%%zmm0,%%zmm0\n\t"\
		"vaddpd	%%zmm3,%%zmm1,%%zmm1\n\t"\
		"vmovaps	%%zmm0		,    (%%rax)\n\t"\
		"vmovaps	%%zmm1		,0x40(%%rax)\n\t"\
		"vmulpd	%%zmm6,%%zmm2,%%zmm2\n\t"\
		"vmulpd	%%zmm6,%%zmm3,%%zmm3\n\t"\
		"vmulpd	%%zmm7,%%zmm4,%%zmm4\n\t"\
		"vmulpd	%%zmm7,%%zmm5,%%zmm5\n\t"\
		"vaddpd	%%zmm0,%%zmm2,%%zmm2\n\t"\
		"vaddpd	%%zmm1,%%zmm3,%%zmm3\n\t"\
		"vmovaps	%%zmm2		,%%zmm0\n\t"\
		"vmovaps	%%zmm3		,%%zmm1\n\t"\
		"vaddpd	%%zmm5,%%zmm2,%%zmm2\n\t"\
		"vsubpd	%%zmm4,%%zmm3,%%zmm3\n\t"\
		"vsubpd	%%zmm5,%%zmm0,%%zmm0\n\t"\
		"vaddpd	%%zmm4,%%zmm1,%%zmm1\n\t"\
		"vmovaps	%%zmm2		,    (%%rbx)\n\t"\
		"vmovaps	%%zmm3		,0x40(%%rbx)\n\t"\
		"vmovaps	%%zmm0		,    (%%rcx)\n\t"\
		"vmovaps	%%zmm1		,0x40(%%rcx)\n\t"\
	/*****************************************/\
	/* now do three more radix-3 transforms: */\
	/*****************************************/\
	/* Block 1: */\
	"movq	%[__o0]	,%%rsi 	\n\t"/* __o0-8: esi,edi store output addresses throughout */\
	"movq	%[__i0]	,%%rax\n\t"\
	"movq	%[__i3]	,%%rbx\n\t"\
	"movq	%[__i6]	,%%rcx\n\t"\
		"vmovaps	    (%%rbx)	,%%zmm4\n\t"\
		"vmovaps	0x40(%%rbx)	,%%zmm5\n\t"\
		"vmovaps	    (%%rcx)	,%%zmm2\n\t"\
		"vmovaps	0x40(%%rcx)	,%%zmm3\n\t"\
		"vmovaps	    (%%rax)	,%%zmm0\n\t"\
		"vmovaps	0x40(%%rax)	,%%zmm1\n\t"\
		"vsubpd	%%zmm2,%%zmm4,%%zmm4\n\t"\
		"vsubpd	%%zmm3,%%zmm5,%%zmm5\n\t"\
		"vaddpd	%%zmm2,%%zmm2,%%zmm2\n\t"\
		"vaddpd	%%zmm3,%%zmm3,%%zmm3\n\t"\
		"vaddpd	%%zmm4,%%zmm2,%%zmm2\n\t"\
		"vaddpd	%%zmm5,%%zmm3,%%zmm3\n\t"\
		"vaddpd	%%zmm2,%%zmm0,%%zmm0\n\t"\
		"vaddpd	%%zmm3,%%zmm1,%%zmm1\n\t"\
		"vmovaps	%%zmm0		,    (%%rsi)\n\t"\
		"vmovaps	%%zmm1		,0x40(%%rsi)\n\t"\
		"vmulpd	%%zmm6,%%zmm2,%%zmm2\n\t"\
		"vmulpd	%%zmm6,%%zmm3,%%zmm3\n\t"\
		"vmulpd	%%zmm7,%%zmm4,%%zmm4\n\t"\
		"vmulpd	%%zmm7,%%zmm5,%%zmm5\n\t"\
		"vaddpd	%%zmm0,%%zmm2,%%zmm2\n\t"\
		"vaddpd	%%zmm1,%%zmm3,%%zmm3\n\t"\
		"vmovaps	%%zmm2		,%%zmm0\n\t"\
		"vmovaps	%%zmm3		,%%zmm1\n\t"\
	"movq	%[__o3]	,%%rdi 	\n\t"\
	"movq	%[__o6]	,%%rsi 	\n\t"\
		"vaddpd	%%zmm5,%%zmm2,%%zmm2\n\t"\
		"vsubpd	%%zmm4,%%zmm3,%%zmm3\n\t"\
		"vsubpd	%%zmm5,%%zmm0,%%zmm0\n\t"\
		"vaddpd	%%zmm4,%%zmm1,%%zmm1\n\t"\
		"vmovaps	%%zmm2		,    (%%rdi)\n\t"\
		"vmovaps	%%zmm3		,0x40(%%rdi)\n\t"\
		"vmovaps	%%zmm0		,    (%%rsi)\n\t"\
		"vmovaps	%%zmm1		,0x40(%%rsi)\n\t"\
	/* Block 2: */\
	"movq	%[__i1]	,%%rax\n\t"\
	"movq	%[__i4]	,%%rbx\n\t"\
	"movq	%[__i7]	,%%rcx\n\t"\
		"vmovaps	    (%%rbx)	,%%zmm2\n\t"\
		"vmovaps	0x40(%%rbx)	,%%zmm3\n\t"\
		"vmovaps	    (%%rdx)	,%%zmm6\n\t"/* c1 */\
		"vmovaps	0x40(%%rdx)	,%%zmm7\n\t"\
		"vmovaps	%%zmm2		,%%zmm0\n\t"\
		"vmovaps	%%zmm3		,%%zmm1\n\t"\
		"vmulpd	%%zmm6,%%zmm2,%%zmm2\n\t"\
		"vmulpd	%%zmm6,%%zmm3,%%zmm3\n\t"\
		"vmulpd	%%zmm7,%%zmm0,%%zmm0\n\t"\
		"vmulpd	%%zmm7,%%zmm1,%%zmm1\n\t"\
		"vaddpd	%%zmm1,%%zmm2,%%zmm2\n\t"\
		"vsubpd	%%zmm0,%%zmm3,%%zmm3\n\t"\
		"vmovaps	    (%%rcx)	,%%zmm4\n\t"\
		"vmovaps	0x40(%%rcx)	,%%zmm5\n\t"\
		"vmovaps	0x80(%%rdx)	,%%zmm6\n\t"/* c2 */\
		"vmovaps	0xc0(%%rdx)	,%%zmm7\n\t"\
		"vmovaps	%%zmm4		,%%zmm0\n\t"\
		"vmovaps	%%zmm5		,%%zmm1\n\t"\
		"vmulpd	%%zmm6,%%zmm4,%%zmm4\n\t"\
		"vmulpd	%%zmm6,%%zmm5,%%zmm5\n\t"\
		"vmulpd	%%zmm7,%%zmm0,%%zmm0\n\t"\
		"vmulpd	%%zmm7,%%zmm1,%%zmm1\n\t"\
		"vaddpd	%%zmm1,%%zmm4,%%zmm4\n\t"\
		"vsubpd	%%zmm0,%%zmm5,%%zmm5\n\t"\
		"vmovaps	0x100(%%rdx)	,%%zmm6\n\t"/* c3m1 */\
		"vmovaps	0x140(%%rdx)	,%%zmm7\n\t"\
		"vmovaps	    (%%rax)	,%%zmm0\n\t"\
		"vmovaps	0x40(%%rax)	,%%zmm1\n\t"\
	"movq	%[__o7]	,%%rsi 	\n\t"\
		"vsubpd	%%zmm4,%%zmm2,%%zmm2\n\t"\
		"vsubpd	%%zmm5,%%zmm3,%%zmm3\n\t"\
		"vaddpd	%%zmm4,%%zmm4,%%zmm4\n\t"\
		"vaddpd	%%zmm5,%%zmm5,%%zmm5\n\t"\
		"vaddpd	%%zmm2,%%zmm4,%%zmm4\n\t"\
		"vaddpd	%%zmm3,%%zmm5,%%zmm5\n\t"\
		"vaddpd	%%zmm4,%%zmm0,%%zmm0\n\t"\
		"vaddpd	%%zmm5,%%zmm1,%%zmm1\n\t"\
		"vmovaps	%%zmm0		,    (%%rsi)\n\t"\
		"vmovaps	%%zmm1		,0x40(%%rsi)\n\t"\
		"vmulpd	%%zmm6,%%zmm4,%%zmm4\n\t"\
		"vmulpd	%%zmm6,%%zmm5,%%zmm5\n\t"\
		"vmulpd	%%zmm7,%%zmm2,%%zmm2\n\t"\
		"vmulpd	%%zmm7,%%zmm3,%%zmm3\n\t"\
		"vaddpd	%%zmm0,%%zmm4,%%zmm4\n\t"\
		"vaddpd	%%zmm1,%%zmm5,%%zmm5\n\t"\
	"movq	%[__o1]	,%%rdi 	\n\t"\
	"movq	%[__o4]	,%%rsi 	\n\t"\
		"vmovaps	%%zmm4		,%%zmm0\n\t"\
		"vmovaps	%%zmm5		,%%zmm1\n\t"\
		"vaddpd	%%zmm3,%%zmm4,%%zmm4\n\t"\
		"vsubpd	%%zmm2,%%zmm5,%%zmm5\n\t"\
		"vsubpd	%%zmm3,%%zmm0,%%zmm0\n\t"\
		"vaddpd	%%zmm2,%%zmm1,%%zmm1\n\t"\
		"vmovaps	%%zmm4		,    (%%rdi)\n\t"\
		"vmovaps	%%zmm5		,0x40(%%rdi)\n\t"\
		"vmovaps	%%zmm0		,    (%%rsi)\n\t"\
		"vmovaps	%%zmm1		,0x40(%%rsi)\n\t"\
	/* Block 3: */\
	"movq	%[__i2]	,%%rax\n\t"\
	"movq	%[__i5]	,%%rbx\n\t"\
	"movq	%[__i8]	,%%rcx\n\t"\
		"vmovaps	    (%%rbx)	,%%zmm2\n\t"\
		"vmovaps	0x40(%%rbx)	,%%zmm3\n\t"\
		"vmovaps	0x80(%%rdx)	,%%zmm6\n\t"/* c2 */\
		"vmovaps	0xc0(%%rdx)	,%%zmm7\n\t"\
		"vmovaps	%%zmm2		,%%zmm0\n\t"\
		"vmovaps	%%zmm3		,%%zmm1\n\t"\
		"vmulpd	%%zmm6,%%zmm2,%%zmm2\n\t"\
		"vmulpd	%%zmm6,%%zmm3,%%zmm3\n\t"\
		"vmulpd	%%zmm7,%%zmm0,%%zmm0\n\t"\
		"vmulpd	%%zmm7,%%zmm1,%%zmm1\n\t"\
		"vaddpd	%%zmm1,%%zmm2,%%zmm2\n\t"\
		"vsubpd	%%zmm0,%%zmm3,%%zmm3\n\t"\
		"vmovaps	    (%%rcx)	,%%zmm4\n\t"\
		"vmovaps	0x40(%%rcx)	,%%zmm5\n\t"\
		"vmovaps	0x180(%%rdx)	,%%zmm6\n\t"/* c4 */\
		"vmovaps	0x1c0(%%rdx)	,%%zmm7\n\t"\
		"vmovaps	%%zmm4		,%%zmm0\n\t"\
		"vmovaps	%%zmm5		,%%zmm1\n\t"\
		"vmulpd	%%zmm6,%%zmm4,%%zmm4\n\t"\
		"vmulpd	%%zmm6,%%zmm5,%%zmm5\n\t"\
		"vmulpd	%%zmm7,%%zmm0,%%zmm0\n\t"\
		"vmulpd	%%zmm7,%%zmm1,%%zmm1\n\t"\
		"vaddpd	%%zmm1,%%zmm4,%%zmm4\n\t"\
		"vsubpd	%%zmm0,%%zmm5,%%zmm5\n\t"\
		"vmovaps	0x100(%%rdx)	,%%zmm6\n\t"\
		"vmovaps	0x140(%%rdx)	,%%zmm7\n\t"\
		"vmovaps	    (%%rax)	,%%zmm0\n\t"\
		"vmovaps	0x40(%%rax)	,%%zmm1\n\t"\
	"movq	%[__o5]	,%%rsi 	\n\t"\
		"vsubpd	%%zmm4,%%zmm2,%%zmm2\n\t"\
		"vsubpd	%%zmm5,%%zmm3,%%zmm3\n\t"\
		"vaddpd	%%zmm4,%%zmm4,%%zmm4\n\t"\
		"vaddpd	%%zmm5,%%zmm5,%%zmm5\n\t"\
		"vaddpd	%%zmm2,%%zmm4,%%zmm4\n\t"\
		"vaddpd	%%zmm3,%%zmm5,%%zmm5\n\t"\
		"vaddpd	%%zmm4,%%zmm0,%%zmm0\n\t"\
		"vaddpd	%%zmm5,%%zmm1,%%zmm1\n\t"\
		"vmovaps	%%zmm0		,    (%%rsi)\n\t"\
		"vmovaps	%%zmm1		,0x40(%%rsi)\n\t"\
		"vmulpd	%%zmm6,%%zmm4,%%zmm4\n\t"\
		"vmulpd	%%zmm6,%%zmm5,%%zmm5\n\t"\
		"vmulpd	%%zmm7,%%zmm2,%%zmm2\n\t"\
		"vmulpd	%%zmm7,%%zmm3,%%zmm3\n\t"\
		"vaddpd	%%zmm0,%%zmm4,%%zmm4\n\t"\
		"vaddpd	%%zmm1,%%zmm5,%%zmm5\n\t"\
	"movq	%[__o8]	,%%rdi 	\n\t"\
	"movq	%[__o2]	,%%rsi 	\n\t"\
		"vmovaps	%%zmm4		,%%zmm0\n\t"\
		"vmovaps	%%zmm5		,%%zmm1\n\t"\
		"vaddpd	%%zmm3,%%zmm4,%%zmm4\n\t"\
		"vsubpd	%%zmm2,%%zmm5,%%zmm5\n\t"\
		"vsubpd	%%zmm3,%%zmm0,%%zmm0\n\t"\
		"vaddpd	%%zmm2,%%zmm1,%%zmm1\n\t"\
		"vmovaps	%%zmm4		,    (%%rdi)\n\t"\
		"vmovaps	%%zmm5		,0x40(%%rdi)\n\t"\
		"vmovaps	%%zmm0		,    (%%rsi)\n\t"\
		"vmovaps	%%zmm1		,0x40(%%rsi)\n\t"\
		:					/* outputs: none */\
		: [__i0] "m" (Xi0)	/* All inputs from memory addresses here */\
		 ,[__i1] "m" (Xi1)\
		 ,[__i2] "m" (Xi2)\
		 ,[__i3] "m" (Xi3)\
		 ,[__i4] "m" (Xi4)\
		 ,[__i5] "m" (Xi5)\
		 ,[__i6] "m" (Xi6)\
		 ,[__i7] "m" (Xi7)\
		 ,[__i8] "m" (Xi8)\
		 ,[__c1] "m" (Xc1)\
		 ,[__o0] "m" (Xo0)\
		 ,[__o1] "m" (Xo1)\
		 ,[__o2] "m" (Xo2)\
		 ,[__o3] "m" (Xo3)\
		 ,[__o4] "m" (Xo4)\
		 ,[__o5] "m" (Xo5)\
		 ,[__o6] "m" (Xo6)\
		 ,[__o7] "m" (Xo7)\
		 ,[__o8] "m" (Xo8)\
		: "cc","memory","rax","rbx","rcx","rdx","rsi","rdi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
		);\
	}

#elif defined(USE_AVX2)

	// Nov 2014: Doubled-data versions of the AVX radix-9 DFT macros are no faster in AVX mode, but will help
	// mitigate FMA latency in AVX2 mode. Opcount: [92 ADD, 88 MUL/FMA].
	//...Radix-9 DIF: First set of Ins in memory locations __i0-8; Outs at memory locations __o0-8, assumed
	// disjoint with inputs.
	// Due to GCC macro argc limit of 30, to enable 16-register data-doubled version
	// of the radix-9 macros use 2 length-9 ptr arrays for second set of IOs:
	//
	#define SSE2_RADIX_09_DIF_X2(Xi0,Xi1,Xi2,Xi3,Xi4,Xi5,Xi6,Xi7,Xi8, Xc1,Xtwo, Xo0,Xo1,Xo2,Xo3,Xo4,Xo5,Xo6,Xo7,Xo8, Xiptr,Xoptr)\
	{\
	__asm__ volatile (												/*********** 2nd set of data: ********/\
	"movq	%[__i0]	,%%rax	\n\t"/* __i0; r[abc]x store set1 I-addrs */"movq	%[__iptr],%%r8		\n\t"/* Pointer to  input-ptr array */\
	"movq	%[__o0]	,%%rsi	\n\t"/* __o0; rsi,rdi store set1 O-addrs */"movq	%[__optr],%%r9		\n\t"/* Pointer to output-ptr array */\
	"movq	%[__c1]	,%%rdx	\n\t"/* rdx has trig addr for both sides */"movq	(%%r8),%%r10 		\n\t"/* __i0; r10-12 store  input addresses */\
	"																	movq	(%%r9),%%r13 		\n\t"/* __o0; r13,14 store output addresses */\
	/* Block 1: */\
		"movq	%%rax		,%%rbx							\n\t		movq	%%r10		,%%r11			\n\t"\
		"movq	%%rax		,%%rcx							\n\t		movq	%%r10		,%%r12			\n\t"\
		"movq	%[__i3]	,%%rbx 								\n\t		movq	0x18(%%r8),%%r11 			\n\t"/* i3 */\
		"movq	%[__i6]	,%%rcx 								\n\t		movq	0x30(%%r8),%%r12 			\n\t"/* i6 */\
		"vmovaps	    (%%rbx)	,%%ymm2						\n\t		vmovaps	    (%%r11)	,%%ymm10		\n\t"\
		"vmovaps	0x20(%%rbx)	,%%ymm3						\n\t		vmovaps	0x20(%%r11)	,%%ymm11		\n\t"\
		"vmovaps	    (%%rax)	,%%ymm0						\n\t		vmovaps	    (%%r10)	,%%ymm8 		\n\t"\
		"vmovaps	0x20(%%rax)	,%%ymm1						\n\t		vmovaps	0x20(%%r10)	,%%ymm9 		\n\t"\
		"vmovaps	    (%%rcx)	,%%ymm6						\n\t		vmovaps	    (%%r12)	,%%ymm14		\n\t"\
		"vmovaps	0x20(%%rcx)	,%%ymm7						\n\t		vmovaps	0x20(%%r12)	,%%ymm15		\n\t"\
		"vmovaps	%%ymm2		,%%ymm4						\n\t		vmovaps	%%ymm10		,%%ymm12		\n\t"\
		"vmovaps	%%ymm3		,%%ymm5						\n\t		vmovaps	%%ymm11		,%%ymm13		\n\t"\
		"vaddpd	%%ymm6,%%ymm2,%%ymm2						\n\t		vaddpd	%%ymm14,%%ymm10,%%ymm10		\n\t"\
		"vaddpd	%%ymm7,%%ymm3,%%ymm3						\n\t		vaddpd	%%ymm15,%%ymm11,%%ymm11		\n\t"\
		"vsubpd	%%ymm6,%%ymm4,%%ymm4						\n\t		vsubpd	%%ymm14,%%ymm12,%%ymm12		\n\t"\
		"vsubpd	%%ymm7,%%ymm5,%%ymm5						\n\t		vsubpd	%%ymm15,%%ymm13,%%ymm13		\n\t"\
		"vaddpd	%%ymm2,%%ymm0,%%ymm0						\n\t		vaddpd	%%ymm10,%%ymm8 ,%%ymm8 		\n\t"\
		"vaddpd	%%ymm3,%%ymm1,%%ymm1						\n\t		vaddpd	%%ymm11,%%ymm9 ,%%ymm9 		\n\t"\
		"vmovaps	0x80(%%rdx)	,%%ymm6			\n\t"/* c3m1 */\
		"vmovaps	0xa0(%%rdx)	,%%ymm7			\n\t"\
		"vmovaps	%%ymm0		,    (%%rsi)	/* <- t00 */\n\t		vmovaps	%%ymm8 		,    (%%r13)	\n\t"/* <- t00 */\
		"vmovaps	%%ymm1		,0x20(%%rsi)	/* <- t01 */\n\t		vmovaps	%%ymm9 		,0x20(%%r13)	\n\t"/* <- t01 */\
	"movq	%[__o1]	,%%rdi									\n\t	movq	0x08(%%r9),%%r14				\n\t"/* o1 */\
	"movq	%[__o2]	,%%rsi									\n\t	movq	0x10(%%r9),%%r13				\n\t"/* o2 */\
	" vfmadd132pd	%%ymm6,%%ymm0,%%ymm2 					\n\t	 vfmadd132pd	%%ymm6,%%ymm8,%%ymm10	\n\t"\
	" vfmadd132pd	%%ymm6,%%ymm1,%%ymm3 					\n\t	 vfmadd132pd	%%ymm6,%%ymm9,%%ymm11	\n\t"\
		"vmovaps	%%ymm2,%%ymm0							\n\t		vmovaps	%%ymm10,%%ymm8 				\n\t"\
		"vmovaps	%%ymm3,%%ymm1							\n\t		vmovaps	%%ymm11,%%ymm9 				\n\t"\
	" vfmadd231pd	%%ymm7,%%ymm5,%%ymm0 					\n\t	 vfmadd231pd	%%ymm7,%%ymm13,%%ymm8 	\n\t"\
	"vfnmadd231pd	%%ymm7,%%ymm4,%%ymm1 					\n\t	vfnmadd231pd	%%ymm7,%%ymm12,%%ymm9 	\n\t"\
	"vfnmadd231pd	%%ymm7,%%ymm5,%%ymm2 					\n\t	vfnmadd231pd	%%ymm7,%%ymm13,%%ymm10	\n\t"\
	" vfmadd231pd	%%ymm7,%%ymm4,%%ymm3 					\n\t	 vfmadd231pd	%%ymm7,%%ymm12,%%ymm11	\n\t"\
		"vmovaps	%%ymm0		,    (%%rsi)	/* <- t04 */\n\t		vmovaps	%%ymm8 		,    (%%r13)	\n\t"/* <- t04 */\
		"vmovaps	%%ymm1		,0x20(%%rsi)	/* <- t05 */\n\t		vmovaps	%%ymm9 		,0x20(%%r13)	\n\t"/* <- t05 */\
		"vmovaps	%%ymm2		,    (%%rdi)	/* <- t02 */\n\t		vmovaps	%%ymm10		,    (%%r14)	\n\t"/* <- t02 */\
		"vmovaps	%%ymm3		,0x20(%%rdi)	/* <- t03 */\n\t		vmovaps	%%ymm11		,0x20(%%r14)	\n\t"/* <- t03 */\
	/* Block 2: */\
		"movq	%[__i1]	,%%rax 								\n\t		movq	0x08(%%r8),%%r10 			\n\t"/* i1 */\
		"movq	%[__i4]	,%%rbx 								\n\t		movq	0x20(%%r8),%%r11 			\n\t"/* i4 */\
		"movq	%[__i7]	,%%rcx 								\n\t		movq	0x38(%%r8),%%r12 			\n\t"/* i7 */\
		"vmovaps	    (%%rbx)	,%%ymm4						\n\t		vmovaps	    (%%r11)	,%%ymm12		\n\t"\
		"vmovaps	0x20(%%rbx)	,%%ymm5						\n\t		vmovaps	0x20(%%r11)	,%%ymm13		\n\t"\
		"vmovaps	    (%%rax)	,%%ymm0						\n\t		vmovaps	    (%%r10)	,%%ymm8 		\n\t"\
		"vmovaps	0x20(%%rax)	,%%ymm1						\n\t		vmovaps	0x20(%%r10)	,%%ymm9 		\n\t"\
		"vmovaps	    (%%rcx)	,%%ymm2						\n\t		vmovaps	    (%%r12)	,%%ymm10		\n\t"\
		"vmovaps	0x20(%%rcx)	,%%ymm3						\n\t		vmovaps	0x20(%%r12)	,%%ymm11		\n\t"\
	"movq	%[__o3]	,%%rdi									\n\t	movq	0x18(%%r9),%%r14				\n\t"/* o3 */\
		"vsubpd	%%ymm2,%%ymm4,%%ymm4						\n\t		vsubpd	%%ymm10,%%ymm12,%%ymm12		\n\t"\
		"vsubpd	%%ymm3,%%ymm5,%%ymm5						\n\t		vsubpd	%%ymm11,%%ymm13,%%ymm13		\n\t"\
		"vaddpd	%%ymm2,%%ymm2,%%ymm2						\n\t		vaddpd	%%ymm10,%%ymm10,%%ymm10		\n\t"\
		"vaddpd	%%ymm3,%%ymm3,%%ymm3						\n\t		vaddpd	%%ymm11,%%ymm11,%%ymm11		\n\t"\
		"vaddpd	%%ymm4,%%ymm2,%%ymm2						\n\t		vaddpd	%%ymm12,%%ymm10,%%ymm10		\n\t"\
		"vaddpd	%%ymm5,%%ymm3,%%ymm3						\n\t		vaddpd	%%ymm13,%%ymm11,%%ymm11		\n\t"\
		"vaddpd	%%ymm2,%%ymm0,%%ymm0						\n\t		vaddpd	%%ymm10,%%ymm8 ,%%ymm8 		\n\t"\
		"vaddpd	%%ymm3,%%ymm1,%%ymm1						\n\t		vaddpd	%%ymm11,%%ymm9 ,%%ymm9 		\n\t"\
		"vmovaps	%%ymm0		,    (%%rdi)	/* <- t06 */\n\t		vmovaps	%%ymm8 		,    (%%r14)	\n\t"/* <- t06 */\
		"vmovaps	%%ymm1		,0x20(%%rdi)	/* <- t07 */\n\t		vmovaps	%%ymm9 		,0x20(%%r14)	\n\t"/* <- t07 */\
	"movq	%[__o4]	,%%rsi									\n\t	movq	0x20(%%r9),%%r13				\n\t"/* o4 */\
	"movq	%[__o5]	,%%rdi									\n\t	movq	0x28(%%r9),%%r14				\n\t"/* o5 */\
	" vfmadd132pd	%%ymm6,%%ymm0,%%ymm2 					\n\t	 vfmadd132pd	%%ymm6,%%ymm8,%%ymm10	\n\t"\
	" vfmadd132pd	%%ymm6,%%ymm1,%%ymm3 					\n\t	 vfmadd132pd	%%ymm6,%%ymm9,%%ymm11	\n\t"\
		"vmovaps	%%ymm2,%%ymm0							\n\t		vmovaps	%%ymm10,%%ymm8 				\n\t"\
		"vmovaps	%%ymm3,%%ymm1							\n\t		vmovaps	%%ymm11,%%ymm9 				\n\t"\
	" vfmadd231pd	%%ymm7,%%ymm5,%%ymm0 					\n\t	 vfmadd231pd	%%ymm7,%%ymm13,%%ymm8 	\n\t"\
	"vfnmadd231pd	%%ymm7,%%ymm4,%%ymm1 					\n\t	vfnmadd231pd	%%ymm7,%%ymm12,%%ymm9 	\n\t"\
	"vfnmadd231pd	%%ymm7,%%ymm5,%%ymm2 					\n\t	vfnmadd231pd	%%ymm7,%%ymm13,%%ymm10	\n\t"\
	" vfmadd231pd	%%ymm7,%%ymm4,%%ymm3 					\n\t	 vfmadd231pd	%%ymm7,%%ymm12,%%ymm11	\n\t"\
		"vmovaps	%%ymm0		,    (%%rdi)	/* <- t0a */\n\t		vmovaps	%%ymm8 		,    (%%r14)	\n\t"/* <- t0a */\
		"vmovaps	%%ymm1		,0x20(%%rdi)	/* <- t0b */\n\t		vmovaps	%%ymm9 		,0x20(%%r14)	\n\t"/* <- t0b */\
		"vmovaps	%%ymm2		,    (%%rsi)	/* <- t08 */\n\t		vmovaps	%%ymm10		,    (%%r13)	\n\t"/* <- t08 */\
		"vmovaps	%%ymm3		,0x20(%%rsi)	/* <- t09 */\n\t		vmovaps	%%ymm11		,0x20(%%r13)	\n\t"/* <- t09 */\
	/* Block 3: */\
		"movq	%[__i2]	,%%rax 								\n\t		movq	0x10(%%r8),%%r10 			\n\t"/* i2 */\
		"movq	%[__i5]	,%%rbx 								\n\t		movq	0x28(%%r8),%%r11 			\n\t"/* i5 */\
		"movq	%[__i8]	,%%rcx 								\n\t		movq	0x40(%%r8),%%r12 			\n\t"/* i8 */\
		"vmovaps	    (%%rbx)	,%%ymm4						\n\t		vmovaps	    (%%r11)	,%%ymm12		\n\t"\
		"vmovaps	0x20(%%rbx)	,%%ymm5						\n\t		vmovaps	0x20(%%r11)	,%%ymm13		\n\t"\
		"vmovaps	    (%%rax)	,%%ymm0						\n\t		vmovaps	    (%%r10)	,%%ymm8 		\n\t"\
		"vmovaps	0x20(%%rax)	,%%ymm1						\n\t		vmovaps	0x20(%%r10)	,%%ymm9 		\n\t"\
		"vmovaps	    (%%rcx)	,%%ymm2						\n\t		vmovaps	    (%%r12)	,%%ymm10		\n\t"\
		"vmovaps	0x20(%%rcx)	,%%ymm3						\n\t		vmovaps	0x20(%%r12)	,%%ymm11		\n\t"\
	"movq	%[__o6]	,%%rsi									\n\t	movq	0x30(%%r9),%%r13				\n\t"/* o6 */\
		"vsubpd	%%ymm2,%%ymm4,%%ymm4						\n\t		vsubpd	%%ymm10,%%ymm12,%%ymm12		\n\t"\
		"vsubpd	%%ymm3,%%ymm5,%%ymm5						\n\t		vsubpd	%%ymm11,%%ymm13,%%ymm13		\n\t"\
		"vaddpd	%%ymm2,%%ymm2,%%ymm2						\n\t		vaddpd	%%ymm10,%%ymm10,%%ymm10		\n\t"\
		"vaddpd	%%ymm3,%%ymm3,%%ymm3						\n\t		vaddpd	%%ymm11,%%ymm11,%%ymm11		\n\t"\
		"vaddpd	%%ymm4,%%ymm2,%%ymm2						\n\t		vaddpd	%%ymm12,%%ymm10,%%ymm10		\n\t"\
		"vaddpd	%%ymm5,%%ymm3,%%ymm3						\n\t		vaddpd	%%ymm13,%%ymm11,%%ymm11		\n\t"\
		"vaddpd	%%ymm2,%%ymm0,%%ymm0						\n\t		vaddpd	%%ymm10,%%ymm8 ,%%ymm8 		\n\t"\
		"vaddpd	%%ymm3,%%ymm1,%%ymm1						\n\t		vaddpd	%%ymm11,%%ymm9 ,%%ymm9 		\n\t"\
		"vmovaps	%%ymm0		,    (%%rsi)	/* <- t0c */\n\t		vmovaps	%%ymm8 		,    (%%r13)	\n\t"/* <- t0c */\
		"vmovaps	%%ymm1		,0x20(%%rsi)	/* <- t0d */\n\t		vmovaps	%%ymm9 		,0x20(%%r13)	\n\t"/* <- t0d */\
	"movq	%[__o7]	,%%rdi									\n\t	movq	0x38(%%r9),%%r14				\n\t"/* o7 */\
	"movq	%[__o8]	,%%rsi									\n\t	movq	0x40(%%r9),%%r13				\n\t"/* o8 */\
	" vfmadd132pd	%%ymm6,%%ymm0,%%ymm2 					\n\t	 vfmadd132pd	%%ymm6,%%ymm8,%%ymm10	\n\t"\
	" vfmadd132pd	%%ymm6,%%ymm1,%%ymm3 					\n\t	 vfmadd132pd	%%ymm6,%%ymm9,%%ymm11	\n\t"\
		"vmovaps	%%ymm2,%%ymm0							\n\t		vmovaps	%%ymm10,%%ymm8 				\n\t"\
		"vmovaps	%%ymm3,%%ymm1							\n\t		vmovaps	%%ymm11,%%ymm9 				\n\t"\
	" vfmadd231pd	%%ymm7,%%ymm5,%%ymm0 					\n\t	 vfmadd231pd	%%ymm7,%%ymm13,%%ymm8 	\n\t"\
	"vfnmadd231pd	%%ymm7,%%ymm4,%%ymm1 					\n\t	vfnmadd231pd	%%ymm7,%%ymm12,%%ymm9 	\n\t"\
	"vfnmadd231pd	%%ymm7,%%ymm5,%%ymm2 					\n\t	vfnmadd231pd	%%ymm7,%%ymm13,%%ymm10	\n\t"\
	" vfmadd231pd	%%ymm7,%%ymm4,%%ymm3 					\n\t	 vfmadd231pd	%%ymm7,%%ymm12,%%ymm11	\n\t"\
		"vmovaps	%%ymm0		,    (%%rsi)	/* <- t0g */\n\t		vmovaps	%%ymm8 		,    (%%r13)	\n\t"/* <- t0g */\
		"vmovaps	%%ymm1		,0x20(%%rsi)	/* <- t0h */\n\t		vmovaps	%%ymm9 		,0x20(%%r13)	\n\t"/* <- t0h */\
		"vmovaps	%%ymm2		,    (%%rdi)	/* <- t0e */\n\t		vmovaps	%%ymm10		,    (%%r14)	\n\t"/* <- t0e */\
		"vmovaps	%%ymm3		,0x20(%%rdi)	/* <- t0f */\n\t		vmovaps	%%ymm11		,0x20(%%r14)	\n\t"/* <- t0f */\
	/******************************************************************************************/\
	/*                        now do three more radix-3 transforms:                           */\
	/******************************************************************************************/\
	/* Block 1: */													/* Block 1: */\
	"movq	%[__o0]	,%%rax									\n\t	movq	    (%%r9),%%r10				\n\t"/* o0 */\
	"movq	%[__o3]	,%%rbx									\n\t	movq	0x18(%%r9),%%r11				\n\t"/* o3 */\
	"movq	%[__o6]	,%%rcx									\n\t	movq	0x30(%%r9),%%r12				\n\t"/* o6 */\
		"vmovaps	    (%%rbx)	,%%ymm4						\n\t		vmovaps	    (%%r11)	,%%ymm12		\n\t"\
		"vmovaps	0x20(%%rbx)	,%%ymm5						\n\t		vmovaps	0x20(%%r11)	,%%ymm13		\n\t"\
		"vmovaps	    (%%rcx)	,%%ymm2						\n\t		vmovaps	    (%%r12)	,%%ymm10		\n\t"\
		"vmovaps	0x20(%%rcx)	,%%ymm3						\n\t		vmovaps	0x20(%%r12)	,%%ymm11		\n\t"\
		"vmovaps	    (%%rax)	,%%ymm0						\n\t		vmovaps	    (%%r10)	,%%ymm8 		\n\t"\
		"vmovaps	0x20(%%rax)	,%%ymm1						\n\t		vmovaps	0x20(%%r10)	,%%ymm9 		\n\t"\
		"vsubpd	%%ymm2,%%ymm4,%%ymm4						\n\t		vsubpd	%%ymm10,%%ymm12,%%ymm12		\n\t"\
		"vsubpd	%%ymm3,%%ymm5,%%ymm5						\n\t		vsubpd	%%ymm11,%%ymm13,%%ymm13		\n\t"\
		"vaddpd	%%ymm2,%%ymm2,%%ymm2						\n\t		vaddpd	%%ymm10,%%ymm10,%%ymm10		\n\t"\
		"vaddpd	%%ymm3,%%ymm3,%%ymm3						\n\t		vaddpd	%%ymm11,%%ymm11,%%ymm11		\n\t"\
		"vaddpd	%%ymm4,%%ymm2,%%ymm2						\n\t		vaddpd	%%ymm12,%%ymm10,%%ymm10		\n\t"\
		"vaddpd	%%ymm5,%%ymm3,%%ymm3						\n\t		vaddpd	%%ymm13,%%ymm11,%%ymm11		\n\t"\
		"vaddpd	%%ymm2,%%ymm0,%%ymm0						\n\t		vaddpd	%%ymm10,%%ymm8 ,%%ymm8 		\n\t"\
		"vaddpd	%%ymm3,%%ymm1,%%ymm1						\n\t		vaddpd	%%ymm11,%%ymm9 ,%%ymm9 		\n\t"\
		"vmovaps	%%ymm0		,    (%%rax)				\n\t		vmovaps	%%ymm8 		,    (%%r10)	\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rax)				\n\t		vmovaps	%%ymm9 		,0x20(%%r10)	\n\t"\
	" vfmadd132pd	%%ymm6,%%ymm0,%%ymm2 					\n\t	 vfmadd132pd	%%ymm6,%%ymm8,%%ymm10	\n\t"\
	" vfmadd132pd	%%ymm6,%%ymm1,%%ymm3 					\n\t	 vfmadd132pd	%%ymm6,%%ymm9,%%ymm11	\n\t"\
		"vmovaps	%%ymm2,%%ymm0							\n\t		vmovaps	%%ymm10,%%ymm8 				\n\t"\
		"vmovaps	%%ymm3,%%ymm1							\n\t		vmovaps	%%ymm11,%%ymm9 				\n\t"\
	" vfmadd231pd	%%ymm7,%%ymm5,%%ymm0 					\n\t	 vfmadd231pd	%%ymm7,%%ymm13,%%ymm8 	\n\t"\
	"vfnmadd231pd	%%ymm7,%%ymm4,%%ymm1 					\n\t	vfnmadd231pd	%%ymm7,%%ymm12,%%ymm9 	\n\t"\
	"vfnmadd231pd	%%ymm7,%%ymm5,%%ymm2 					\n\t	vfnmadd231pd	%%ymm7,%%ymm13,%%ymm10	\n\t"\
	" vfmadd231pd	%%ymm7,%%ymm4,%%ymm3 					\n\t	 vfmadd231pd	%%ymm7,%%ymm12,%%ymm11	\n\t"\
		"vmovaps	%%ymm0		,    (%%rcx)				\n\t		vmovaps	%%ymm8 		,    (%%r12)	\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rcx)				\n\t		vmovaps	%%ymm9 		,0x20(%%r12)	\n\t"\
		"vmovaps	%%ymm2		,    (%%rbx)				\n\t		vmovaps	%%ymm10		,    (%%r11)	\n\t"\
		"vmovaps	%%ymm3		,0x20(%%rbx)				\n\t		vmovaps	%%ymm11		,0x20(%%r11)	\n\t"\
	/* Block 2: */\
	"movq	%[__o1]	,%%rax									\n\t	movq	0x08(%%r9),%%r10				\n\t"/* o1 */\
	"movq	%[__o4]	,%%rbx									\n\t	movq	0x20(%%r9),%%r11				\n\t"/* o4 */\
	"movq	%[__o7]	,%%rcx									\n\t	movq	0x38(%%r9),%%r12				\n\t"/* o7 */\
		"vmovaps	    (%%rbx)	,%%ymm2						\n\t		vmovaps	    (%%r11)	,%%ymm10		\n\t"\
		"vmovaps	0x20(%%rbx)	,%%ymm3						\n\t		vmovaps	0x20(%%r11)	,%%ymm11		\n\t"\
		"vmovaps	    (%%rdx)	,%%ymm6				\n\t"/* c1 */\
		"vmovaps	0x20(%%rdx)	,%%ymm7				\n\t"\
		"vmovaps	%%ymm2		,%%ymm0						\n\t		vmovaps	%%ymm10		,%%ymm8 		\n\t"\
		"vmovaps	%%ymm3		,%%ymm1						\n\t		vmovaps	%%ymm11		,%%ymm9 		\n\t"\
		"vmulpd		%%ymm6,%%ymm2,%%ymm2					\n\t		vmulpd		%%ymm6,%%ymm10,%%ymm10	\n\t"\
		"vmulpd		%%ymm6,%%ymm3,%%ymm3					\n\t		vmulpd		%%ymm6,%%ymm11,%%ymm11	\n\t"\
	"vfnmadd231pd	%%ymm7,%%ymm1,%%ymm2 					\n\t	vfnmadd231pd	%%ymm7,%%ymm9 ,%%ymm10	\n\t"\
	" vfmadd231pd	%%ymm7,%%ymm0,%%ymm3 					\n\t	 vfmadd231pd	%%ymm7,%%ymm8 ,%%ymm11	\n\t"\
		"vmovaps	    (%%rcx)	,%%ymm4						\n\t		vmovaps	    (%%r12)	,%%ymm12		\n\t"\
		"vmovaps	0x20(%%rcx)	,%%ymm5						\n\t		vmovaps	0x20(%%r12)	,%%ymm13		\n\t"\
		"vmovaps	0x40(%%rdx)	,%%ymm6				\n\t"/* c2 */\
		"vmovaps	0x60(%%rdx)	,%%ymm7				\n\t"\
		"vmovaps	%%ymm4		,%%ymm0						\n\t		vmovaps	%%ymm12		,%%ymm8 		\n\t"\
		"vmovaps	%%ymm5		,%%ymm1						\n\t		vmovaps	%%ymm13		,%%ymm9 		\n\t"\
		"vmulpd		%%ymm6,%%ymm4,%%ymm4					\n\t		vmulpd		%%ymm6,%%ymm12,%%ymm12	\n\t"\
		"vmulpd		%%ymm6,%%ymm5,%%ymm5					\n\t		vmulpd		%%ymm6,%%ymm13,%%ymm13	\n\t"\
	"vfnmadd231pd	%%ymm7,%%ymm1,%%ymm4 					\n\t	vfnmadd231pd	%%ymm7,%%ymm9 ,%%ymm12	\n\t"\
	" vfmadd231pd	%%ymm7,%%ymm0,%%ymm5 					\n\t	 vfmadd231pd	%%ymm7,%%ymm8 ,%%ymm13	\n\t"\
		"vmovaps	0x80(%%rdx)	,%%ymm6				\n\t"/* c3m1 */\
		"vmovaps	0xa0(%%rdx)	,%%ymm7				\n\t"\
		"vmovaps	    (%%rax)	,%%ymm0						\n\t		vmovaps	    (%%r10)	,%%ymm8 		\n\t"\
		"vmovaps	0x20(%%rax)	,%%ymm1						\n\t		vmovaps	0x20(%%r10)	,%%ymm9 		\n\t"\
		"vsubpd	%%ymm4,%%ymm2,%%ymm2						\n\t		vsubpd	%%ymm12,%%ymm10,%%ymm10		\n\t"\
		"vsubpd	%%ymm5,%%ymm3,%%ymm3						\n\t		vsubpd	%%ymm13,%%ymm11,%%ymm11		\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4						\n\t		vaddpd	%%ymm12,%%ymm12,%%ymm12		\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5						\n\t		vaddpd	%%ymm13,%%ymm13,%%ymm13		\n\t"\
		"vaddpd	%%ymm2,%%ymm4,%%ymm4						\n\t		vaddpd	%%ymm10,%%ymm12,%%ymm12		\n\t"\
		"vaddpd	%%ymm3,%%ymm5,%%ymm5						\n\t		vaddpd	%%ymm11,%%ymm13,%%ymm13		\n\t"\
		"vaddpd	%%ymm4,%%ymm0,%%ymm0						\n\t		vaddpd	%%ymm12,%%ymm8 ,%%ymm8 		\n\t"\
		"vaddpd	%%ymm5,%%ymm1,%%ymm1						\n\t		vaddpd	%%ymm13,%%ymm9 ,%%ymm9 		\n\t"\
		"vmovaps	%%ymm0		,    (%%rax)				\n\t		vmovaps	%%ymm8 		,    (%%r10)	\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rax)				\n\t		vmovaps	%%ymm9 		,0x20(%%r10)	\n\t"\
	" vfmadd132pd	%%ymm6,%%ymm0,%%ymm4 					\n\t	 vfmadd132pd	%%ymm6,%%ymm8,%%ymm12	\n\t"\
	" vfmadd132pd	%%ymm6,%%ymm1,%%ymm5 					\n\t	 vfmadd132pd	%%ymm6,%%ymm9,%%ymm13	\n\t"\
		"vmovaps	%%ymm4,%%ymm0							\n\t		vmovaps	%%ymm12,%%ymm8 				\n\t"\
		"vmovaps	%%ymm5,%%ymm1							\n\t		vmovaps	%%ymm13,%%ymm9 				\n\t"\
	" vfmadd231pd	%%ymm7,%%ymm3,%%ymm0 					\n\t	 vfmadd231pd	%%ymm7,%%ymm11,%%ymm8 	\n\t"\
	"vfnmadd231pd	%%ymm7,%%ymm2,%%ymm1 					\n\t	vfnmadd231pd	%%ymm7,%%ymm10,%%ymm9 	\n\t"\
	"vfnmadd231pd	%%ymm7,%%ymm3,%%ymm4 					\n\t	vfnmadd231pd	%%ymm7,%%ymm11,%%ymm12	\n\t"\
	" vfmadd231pd	%%ymm7,%%ymm2,%%ymm5 					\n\t	 vfmadd231pd	%%ymm7,%%ymm10,%%ymm13	\n\t"\
		"vmovaps	%%ymm0		,    (%%rcx)				\n\t		vmovaps	%%ymm8 		,    (%%r12)	\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rcx)				\n\t		vmovaps	%%ymm9 		,0x20(%%r12)	\n\t"\
		"vmovaps	%%ymm4		,    (%%rbx)				\n\t		vmovaps	%%ymm12		,    (%%r11)	\n\t"\
		"vmovaps	%%ymm5		,0x20(%%rbx)				\n\t		vmovaps	%%ymm13		,0x20(%%r11)	\n\t"\
	/* Block 3: */\
	"movq	%[__o2]	,%%rax									\n\t	movq	0x10(%%r9),%%r10				\n\t"/* o2 */\
	"movq	%[__o5]	,%%rbx									\n\t	movq	0x28(%%r9),%%r11				\n\t"/* o5 */\
	"movq	%[__o8]	,%%rcx									\n\t	movq	0x40(%%r9),%%r12				\n\t"/* o8 */\
		"vmovaps	    (%%rbx)	,%%ymm2						\n\t		vmovaps	    (%%r11)	,%%ymm10		\n\t"\
		"vmovaps	0x20(%%rbx)	,%%ymm3						\n\t		vmovaps	0x20(%%r11)	,%%ymm11		\n\t"\
		"vmovaps	0x40(%%rdx)	,%%ymm6				\n\t"/* c2 */\
		"vmovaps	0x60(%%rdx)	,%%ymm7				\n\t"\
		"vmovaps	%%ymm2		,%%ymm0						\n\t		vmovaps	%%ymm10		,%%ymm8 		\n\t"\
		"vmovaps	%%ymm3		,%%ymm1						\n\t		vmovaps	%%ymm11		,%%ymm9 		\n\t"\
		"vmulpd		%%ymm6,%%ymm2,%%ymm2					\n\t		vmulpd		%%ymm6,%%ymm10,%%ymm10	\n\t"\
		"vmulpd		%%ymm6,%%ymm3,%%ymm3					\n\t		vmulpd		%%ymm6,%%ymm11,%%ymm11	\n\t"\
	"vfnmadd231pd	%%ymm7,%%ymm1,%%ymm2 					\n\t	vfnmadd231pd	%%ymm7,%%ymm9 ,%%ymm10	\n\t"\
	" vfmadd231pd	%%ymm7,%%ymm0,%%ymm3 					\n\t	 vfmadd231pd	%%ymm7,%%ymm8 ,%%ymm11	\n\t"\
		"vmovaps	    (%%rcx)	,%%ymm4						\n\t		vmovaps	    (%%r12)	,%%ymm12		\n\t"\
		"vmovaps	0x20(%%rcx)	,%%ymm5						\n\t		vmovaps	0x20(%%r12)	,%%ymm13		\n\t"\
		"vmovaps	0xc0(%%rdx)	,%%ymm6				\n\t"/* c4 */\
		"vmovaps	0xe0(%%rdx)	,%%ymm7				\n\t"\
		"vmovaps	%%ymm4		,%%ymm0						\n\t		vmovaps	%%ymm12		,%%ymm8 		\n\t"\
		"vmovaps	%%ymm5		,%%ymm1						\n\t		vmovaps	%%ymm13		,%%ymm9 		\n\t"\
		"vmulpd		%%ymm6,%%ymm4,%%ymm4					\n\t		vmulpd		%%ymm6,%%ymm12,%%ymm12	\n\t"\
		"vmulpd		%%ymm6,%%ymm5,%%ymm5					\n\t		vmulpd		%%ymm6,%%ymm13,%%ymm13	\n\t"\
	"vfnmadd231pd	%%ymm7,%%ymm1,%%ymm4 					\n\t	vfnmadd231pd	%%ymm7,%%ymm9 ,%%ymm12	\n\t"\
	" vfmadd231pd	%%ymm7,%%ymm0,%%ymm5 					\n\t	 vfmadd231pd	%%ymm7,%%ymm8 ,%%ymm13	\n\t"\
		"vmovaps	0x80(%%rdx)	,%%ymm6				\n\t"/* c3m1 */\
		"vmovaps	0xa0(%%rdx)	,%%ymm7				\n\t"\
		"vmovaps	    (%%rax)	,%%ymm0						\n\t		vmovaps	    (%%r10)	,%%ymm8 		\n\t"\
		"vmovaps	0x20(%%rax)	,%%ymm1						\n\t		vmovaps	0x20(%%r10)	,%%ymm9 		\n\t"\
		"vsubpd	%%ymm4,%%ymm2,%%ymm2						\n\t		vsubpd	%%ymm12,%%ymm10,%%ymm10		\n\t"\
		"vsubpd	%%ymm5,%%ymm3,%%ymm3						\n\t		vsubpd	%%ymm13,%%ymm11,%%ymm11		\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4						\n\t		vaddpd	%%ymm12,%%ymm12,%%ymm12		\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5						\n\t		vaddpd	%%ymm13,%%ymm13,%%ymm13		\n\t"\
		"vaddpd	%%ymm2,%%ymm4,%%ymm4						\n\t		vaddpd	%%ymm10,%%ymm12,%%ymm12		\n\t"\
		"vaddpd	%%ymm3,%%ymm5,%%ymm5						\n\t		vaddpd	%%ymm11,%%ymm13,%%ymm13		\n\t"\
		"vaddpd	%%ymm4,%%ymm0,%%ymm0						\n\t		vaddpd	%%ymm12,%%ymm8 ,%%ymm8 		\n\t"\
		"vaddpd	%%ymm5,%%ymm1,%%ymm1						\n\t		vaddpd	%%ymm13,%%ymm9 ,%%ymm9 		\n\t"\
		"vmovaps	%%ymm0		,    (%%rax)				\n\t		vmovaps	%%ymm8 		,    (%%r10)	\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rax)				\n\t		vmovaps	%%ymm9 		,0x20(%%r10)	\n\t"\
	" vfmadd132pd	%%ymm6,%%ymm0,%%ymm4 					\n\t	 vfmadd132pd	%%ymm6,%%ymm8,%%ymm12	\n\t"\
	" vfmadd132pd	%%ymm6,%%ymm1,%%ymm5 					\n\t	 vfmadd132pd	%%ymm6,%%ymm9,%%ymm13	\n\t"\
		"vmovaps	%%ymm4,%%ymm0							\n\t		vmovaps	%%ymm12,%%ymm8 				\n\t"\
		"vmovaps	%%ymm5,%%ymm1							\n\t		vmovaps	%%ymm13,%%ymm9 				\n\t"\
	" vfmadd231pd	%%ymm7,%%ymm3,%%ymm0 					\n\t	 vfmadd231pd	%%ymm7,%%ymm11,%%ymm8 	\n\t"\
	"vfnmadd231pd	%%ymm7,%%ymm2,%%ymm1 					\n\t	vfnmadd231pd	%%ymm7,%%ymm10,%%ymm9 	\n\t"\
	"vfnmadd231pd	%%ymm7,%%ymm3,%%ymm4 					\n\t	vfnmadd231pd	%%ymm7,%%ymm11,%%ymm12	\n\t"\
	" vfmadd231pd	%%ymm7,%%ymm2,%%ymm5 					\n\t	 vfmadd231pd	%%ymm7,%%ymm10,%%ymm13	\n\t"\
		"vmovaps	%%ymm0		,    (%%rcx)				\n\t		vmovaps	%%ymm8 		,    (%%r12)	\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rcx)				\n\t		vmovaps	%%ymm9 		,0x20(%%r12)	\n\t"\
		"vmovaps	%%ymm4		,    (%%rbx)				\n\t		vmovaps	%%ymm12		,    (%%r11)	\n\t"\
		"vmovaps	%%ymm5		,0x20(%%rbx)				\n\t		vmovaps	%%ymm13		,0x20(%%r11)	\n\t"\
		:					/* outputs: none */\
		: [__i0] "m" (Xi0)	/* All inputs from memory addresses here */\
		 ,[__i1] "m" (Xi1)\
		 ,[__i2] "m" (Xi2)\
		 ,[__i3] "m" (Xi3)\
		 ,[__i4] "m" (Xi4)\
		 ,[__i5] "m" (Xi5)\
		 ,[__i6] "m" (Xi6)\
		 ,[__i7] "m" (Xi7)\
		 ,[__i8] "m" (Xi8)\
		 ,[__c1] "m" (Xc1)\
		 ,[__two] "m" (Xtwo)\
		 ,[__o0] "m" (Xo0)\
		 ,[__o1] "m" (Xo1)\
		 ,[__o2] "m" (Xo2)\
		 ,[__o3] "m" (Xo3)\
		 ,[__o4] "m" (Xo4)\
		 ,[__o5] "m" (Xo5)\
		 ,[__o6] "m" (Xo6)\
		 ,[__o7] "m" (Xo7)\
		 ,[__o8] "m" (Xo8)\
		 ,[__iptr] "m" (Xiptr)	/* Due to GCC macro argc limit of 30, to enable 16-register data-doubled version */\
		 ,[__optr] "m" (Xoptr)	/* of the radix-9 macros use 2 length-9 ptr arrays for second set of IOs. */\
		: "cc","memory","rax","rbx","rcx","rdx","rsi","rdi","r8","r9","r10","r11","r12","r13","r14","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
		);\
	}

	/*...Radix-9 DIT: Ins in memory locations __i0-8.\
		Outs at memory locations __o0-8, assumed disjoint with inputs:\
		Opcount: [92 ADD, 88 MUL/FMA].
	*/\
	#define SSE2_RADIX_09_DIT_X2(Xi0,Xi1,Xi2,Xi3,Xi4,Xi5,Xi6,Xi7,Xi8, Xc1,Xtwo, Xo0,Xo1,Xo2,Xo3,Xo4,Xo5,Xo6,Xo7,Xo8, Xiptr,Xoptr)\
	{\
	__asm__ volatile (\
																		/*********** 2nd set of data: ********/\
	"																	movq	%[__iptr],%%r8		\n\t"/* Pointer to  input-ptr array */\
	"																	movq	%[__optr],%%r9		\n\t"/* Pointer to output-ptr array */\
	"																	movq	(%%r8),%%r10 		\n\t"/* __i0; r10-12 store  input addresses */\
	"movq	%[__c1]	,%%rdx	\n\t"/* rdx has trig addr for both sides */"movq	(%%r9),%%r13 		\n\t"/* __o0; r13,14 store output addresses */\
	"movq	%[__i0]	,%%rax	\n\t"/* __i0; r[abc]x store set1 I-addrs */"movq	    (%%r8),%%r10				\n\t"/* i0 */\
	"movq	%[__i1]	,%%rbx	\n\t										movq	0x08(%%r8),%%r11				\n\t"/* i1 */\
	"movq	%[__i2]	,%%rcx	\n\t										movq	0x10(%%r8),%%r12				\n\t"/* i2 */\
	/* Block 1: */\
		"vmovaps	    (%%rbx)	,%%ymm2						\n\t		vmovaps	    (%%r11)	,%%ymm10		\n\t"\
		"vmovaps	0x20(%%rbx)	,%%ymm3						\n\t		vmovaps	0x20(%%r11)	,%%ymm11		\n\t"\
		"vmovaps	    (%%rax)	,%%ymm0						\n\t		vmovaps	    (%%r10)	,%%ymm8 		\n\t"\
		"vmovaps	0x20(%%rax)	,%%ymm1						\n\t		vmovaps	0x20(%%r10)	,%%ymm9 		\n\t"\
		"vmovaps	    (%%rcx)	,%%ymm6						\n\t		vmovaps	    (%%r12)	,%%ymm14		\n\t"\
		"vmovaps	0x20(%%rcx)	,%%ymm7						\n\t		vmovaps	0x20(%%r12)	,%%ymm15		\n\t"\
		"vmovaps	%%ymm2		,%%ymm4						\n\t		vmovaps	%%ymm10		,%%ymm12		\n\t"\
		"vmovaps	%%ymm3		,%%ymm5						\n\t		vmovaps	%%ymm11		,%%ymm13		\n\t"\
		"vaddpd	%%ymm6,%%ymm2,%%ymm2						\n\t		vaddpd	%%ymm14,%%ymm10,%%ymm10		\n\t"\
		"vaddpd	%%ymm7,%%ymm3,%%ymm3						\n\t		vaddpd	%%ymm15,%%ymm11,%%ymm11		\n\t"\
		"vsubpd	%%ymm6,%%ymm4,%%ymm4						\n\t		vsubpd	%%ymm14,%%ymm12,%%ymm12		\n\t"\
		"vsubpd	%%ymm7,%%ymm5,%%ymm5						\n\t		vsubpd	%%ymm15,%%ymm13,%%ymm13		\n\t"\
		"vaddpd	%%ymm2,%%ymm0,%%ymm0						\n\t		vaddpd	%%ymm10,%%ymm8 ,%%ymm8 		\n\t"\
		"vaddpd	%%ymm3,%%ymm1,%%ymm1						\n\t		vaddpd	%%ymm11,%%ymm9 ,%%ymm9 		\n\t"\
		"vmovaps	0x80(%%rdx)	,%%ymm6						\n\t"/* c3m1 */\
		"vmovaps	0xa0(%%rdx)	,%%ymm7						\n\t"\
		"vmovaps	%%ymm0		,    (%%rax)				\n\t		vmovaps	%%ymm8 		,    (%%r10)	\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rax)				\n\t		vmovaps	%%ymm9 		,0x20(%%r10)	\n\t"\
	" vfmadd132pd	%%ymm6,%%ymm0,%%ymm2 					\n\t	 vfmadd132pd	%%ymm6,%%ymm8,%%ymm10	\n\t"\
	" vfmadd132pd	%%ymm6,%%ymm1,%%ymm3 					\n\t	 vfmadd132pd	%%ymm6,%%ymm9,%%ymm11	\n\t"\
		"vmovaps	%%ymm2,%%ymm0							\n\t		vmovaps	%%ymm10,%%ymm8 				\n\t"\
		"vmovaps	%%ymm3,%%ymm1							\n\t		vmovaps	%%ymm11,%%ymm9 				\n\t"\
	"vfnmadd231pd	%%ymm7,%%ymm5,%%ymm0 					\n\t	vfnmadd231pd	%%ymm7,%%ymm13,%%ymm8 	\n\t"\
	" vfmadd231pd	%%ymm7,%%ymm4,%%ymm1 					\n\t	 vfmadd231pd	%%ymm7,%%ymm12,%%ymm9 	\n\t"\
	" vfmadd231pd	%%ymm7,%%ymm5,%%ymm2 					\n\t	 vfmadd231pd	%%ymm7,%%ymm13,%%ymm10	\n\t"\
	"vfnmadd231pd	%%ymm7,%%ymm4,%%ymm3 					\n\t	vfnmadd231pd	%%ymm7,%%ymm12,%%ymm11	\n\t"\
		"vmovaps	%%ymm0		,    (%%rcx)				\n\t		vmovaps	%%ymm8 		,    (%%r12)	\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rcx)				\n\t		vmovaps	%%ymm9 		,0x20(%%r12)	\n\t"\
		"vmovaps	%%ymm2		,    (%%rbx)				\n\t		vmovaps	%%ymm10		,    (%%r11)	\n\t"\
		"vmovaps	%%ymm3		,0x20(%%rbx)				\n\t		vmovaps	%%ymm11		,0x20(%%r11)	\n\t"\
	/* Block 2: */\
	"movq	%[__i3]	,%%rax									\n\t	movq	0x18(%%r8),%%r10				\n\t"/* i3 */\
	"movq	%[__i4]	,%%rbx									\n\t	movq	0x20(%%r8),%%r11				\n\t"/* i4 */\
	"movq	%[__i5]	,%%rcx									\n\t	movq	0x28(%%r8),%%r12				\n\t"/* i5 */\
		"vmovaps	    (%%rbx)	,%%ymm4						\n\t		vmovaps	    (%%r11)	,%%ymm12		\n\t"\
		"vmovaps	0x20(%%rbx)	,%%ymm5						\n\t		vmovaps	0x20(%%r11)	,%%ymm13		\n\t"\
		"vmovaps	    (%%rcx)	,%%ymm2						\n\t		vmovaps	    (%%r12)	,%%ymm10		\n\t"\
		"vmovaps	0x20(%%rcx)	,%%ymm3						\n\t		vmovaps	0x20(%%r12)	,%%ymm11		\n\t"\
		"vmovaps	    (%%rax)	,%%ymm0						\n\t		vmovaps	    (%%r10)	,%%ymm8 		\n\t"\
		"vmovaps	0x20(%%rax)	,%%ymm1						\n\t		vmovaps	0x20(%%r10)	,%%ymm9 		\n\t"\
		"vsubpd	%%ymm2,%%ymm4,%%ymm4						\n\t		vsubpd	%%ymm10,%%ymm12,%%ymm12		\n\t"\
		"vsubpd	%%ymm3,%%ymm5,%%ymm5						\n\t		vsubpd	%%ymm11,%%ymm13,%%ymm13		\n\t"\
		"vaddpd	%%ymm2,%%ymm2,%%ymm2						\n\t		vaddpd	%%ymm10,%%ymm10,%%ymm10		\n\t"\
		"vaddpd	%%ymm3,%%ymm3,%%ymm3						\n\t		vaddpd	%%ymm11,%%ymm11,%%ymm11		\n\t"\
		"vaddpd	%%ymm4,%%ymm2,%%ymm2						\n\t		vaddpd	%%ymm12,%%ymm10,%%ymm10		\n\t"\
		"vaddpd	%%ymm5,%%ymm3,%%ymm3						\n\t		vaddpd	%%ymm13,%%ymm11,%%ymm11		\n\t"\
		"vaddpd	%%ymm2,%%ymm0,%%ymm0						\n\t		vaddpd	%%ymm10,%%ymm8 ,%%ymm8 		\n\t"\
		"vaddpd	%%ymm3,%%ymm1,%%ymm1						\n\t		vaddpd	%%ymm11,%%ymm9 ,%%ymm9 		\n\t"\
		"vmovaps	%%ymm0		,    (%%rax)				\n\t		vmovaps	%%ymm8 		,    (%%r10)	\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rax)				\n\t		vmovaps	%%ymm9 		,0x20(%%r10)	\n\t"\
	" vfmadd132pd	%%ymm6,%%ymm0,%%ymm2 					\n\t	 vfmadd132pd	%%ymm6,%%ymm8,%%ymm10	\n\t"\
	" vfmadd132pd	%%ymm6,%%ymm1,%%ymm3 					\n\t	 vfmadd132pd	%%ymm6,%%ymm9,%%ymm11	\n\t"\
		"vmovaps	%%ymm2,%%ymm0							\n\t		vmovaps	%%ymm10,%%ymm8 				\n\t"\
		"vmovaps	%%ymm3,%%ymm1							\n\t		vmovaps	%%ymm11,%%ymm9 				\n\t"\
	"vfnmadd231pd	%%ymm7,%%ymm5,%%ymm0 					\n\t	vfnmadd231pd	%%ymm7,%%ymm13,%%ymm8 	\n\t"\
	" vfmadd231pd	%%ymm7,%%ymm4,%%ymm1 					\n\t	 vfmadd231pd	%%ymm7,%%ymm12,%%ymm9 	\n\t"\
	" vfmadd231pd	%%ymm7,%%ymm5,%%ymm2 					\n\t	 vfmadd231pd	%%ymm7,%%ymm13,%%ymm10	\n\t"\
	"vfnmadd231pd	%%ymm7,%%ymm4,%%ymm3 					\n\t	vfnmadd231pd	%%ymm7,%%ymm12,%%ymm11	\n\t"\
		"vmovaps	%%ymm0		,    (%%rcx)				\n\t		vmovaps	%%ymm8 		,    (%%r12)	\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rcx)				\n\t		vmovaps	%%ymm9 		,0x20(%%r12)	\n\t"\
		"vmovaps	%%ymm2		,    (%%rbx)				\n\t		vmovaps	%%ymm10		,    (%%r11)	\n\t"\
		"vmovaps	%%ymm3		,0x20(%%rbx)				\n\t		vmovaps	%%ymm11		,0x20(%%r11)	\n\t"\
	/* Block 3: */\
	"movq	%[__i6]	,%%rax									\n\t	movq	0x30(%%r8),%%r10				\n\t"/* i6 */\
	"movq	%[__i7]	,%%rbx									\n\t	movq	0x38(%%r8),%%r11				\n\t"/* i7 */\
	"movq	%[__i8]	,%%rcx									\n\t	movq	0x40(%%r8),%%r12				\n\t"/* i8 */\
		"vmovaps	    (%%rbx)	,%%ymm4						\n\t		vmovaps	    (%%r11)	,%%ymm12		\n\t"\
		"vmovaps	0x20(%%rbx)	,%%ymm5						\n\t		vmovaps	0x20(%%r11)	,%%ymm13		\n\t"\
		"vmovaps	    (%%rcx)	,%%ymm2						\n\t		vmovaps	    (%%r12)	,%%ymm10		\n\t"\
		"vmovaps	0x20(%%rcx)	,%%ymm3						\n\t		vmovaps	0x20(%%r12)	,%%ymm11		\n\t"\
		"vmovaps	    (%%rax)	,%%ymm0						\n\t		vmovaps	    (%%r10)	,%%ymm8 		\n\t"\
		"vmovaps	0x20(%%rax)	,%%ymm1						\n\t		vmovaps	0x20(%%r10)	,%%ymm9 		\n\t"\
		"vsubpd	%%ymm2,%%ymm4,%%ymm4						\n\t		vsubpd	%%ymm10,%%ymm12,%%ymm12		\n\t"\
		"vsubpd	%%ymm3,%%ymm5,%%ymm5						\n\t		vsubpd	%%ymm11,%%ymm13,%%ymm13		\n\t"\
		"vaddpd	%%ymm2,%%ymm2,%%ymm2						\n\t		vaddpd	%%ymm10,%%ymm10,%%ymm10		\n\t"\
		"vaddpd	%%ymm3,%%ymm3,%%ymm3						\n\t		vaddpd	%%ymm11,%%ymm11,%%ymm11		\n\t"\
		"vaddpd	%%ymm4,%%ymm2,%%ymm2						\n\t		vaddpd	%%ymm12,%%ymm10,%%ymm10		\n\t"\
		"vaddpd	%%ymm5,%%ymm3,%%ymm3						\n\t		vaddpd	%%ymm13,%%ymm11,%%ymm11		\n\t"\
		"vaddpd	%%ymm2,%%ymm0,%%ymm0						\n\t		vaddpd	%%ymm10,%%ymm8 ,%%ymm8 		\n\t"\
		"vaddpd	%%ymm3,%%ymm1,%%ymm1						\n\t		vaddpd	%%ymm11,%%ymm9 ,%%ymm9 		\n\t"\
		"vmovaps	%%ymm0		,    (%%rax)				\n\t		vmovaps	%%ymm8 		,    (%%r10)	\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rax)				\n\t		vmovaps	%%ymm9 		,0x20(%%r10)	\n\t"\
	" vfmadd132pd	%%ymm6,%%ymm0,%%ymm2 					\n\t	 vfmadd132pd	%%ymm6,%%ymm8,%%ymm10	\n\t"\
	" vfmadd132pd	%%ymm6,%%ymm1,%%ymm3 					\n\t	 vfmadd132pd	%%ymm6,%%ymm9,%%ymm11	\n\t"\
		"vmovaps	%%ymm2,%%ymm0							\n\t		vmovaps	%%ymm10,%%ymm8 				\n\t"\
		"vmovaps	%%ymm3,%%ymm1							\n\t		vmovaps	%%ymm11,%%ymm9 				\n\t"\
	"vfnmadd231pd	%%ymm7,%%ymm5,%%ymm0 					\n\t	vfnmadd231pd	%%ymm7,%%ymm13,%%ymm8 	\n\t"\
	" vfmadd231pd	%%ymm7,%%ymm4,%%ymm1 					\n\t	 vfmadd231pd	%%ymm7,%%ymm12,%%ymm9 	\n\t"\
	" vfmadd231pd	%%ymm7,%%ymm5,%%ymm2 					\n\t	 vfmadd231pd	%%ymm7,%%ymm13,%%ymm10	\n\t"\
	"vfnmadd231pd	%%ymm7,%%ymm4,%%ymm3 					\n\t	vfnmadd231pd	%%ymm7,%%ymm12,%%ymm11	\n\t"\
		"vmovaps	%%ymm0		,    (%%rcx)				\n\t		vmovaps	%%ymm8 		,    (%%r12)	\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rcx)				\n\t		vmovaps	%%ymm9 		,0x20(%%r12)	\n\t"\
		"vmovaps	%%ymm2		,    (%%rbx)				\n\t		vmovaps	%%ymm10		,    (%%r11)	\n\t"\
		"vmovaps	%%ymm3		,0x20(%%rbx)				\n\t		vmovaps	%%ymm11		,0x20(%%r11)	\n\t"\
	/*****************************************/\
	/* now do three more radix-3 transforms: */\
	/*****************************************/\
	/* Block 1: */\
	"movq	%[__o0]	,%%rsi 									\n\t	movq	    (%%r9),%%r13 				\n\t"/* o0 */\
	"movq	%[__i0]	,%%rax									\n\t	movq	    (%%r8),%%r10				\n\t"/* i0 */\
	"movq	%[__i3]	,%%rbx									\n\t	movq	0x18(%%r8),%%r11				\n\t"/* i3 */\
	"movq	%[__i6]	,%%rcx									\n\t	movq	0x30(%%r8),%%r12				\n\t"/* i6 */\
		"vmovaps	    (%%rbx)	,%%ymm4						\n\t		vmovaps	    (%%r11)	,%%ymm12		\n\t"\
		"vmovaps	0x20(%%rbx)	,%%ymm5						\n\t		vmovaps	0x20(%%r11)	,%%ymm13		\n\t"\
		"vmovaps	    (%%rcx)	,%%ymm2						\n\t		vmovaps	    (%%r12)	,%%ymm10		\n\t"\
		"vmovaps	0x20(%%rcx)	,%%ymm3						\n\t		vmovaps	0x20(%%r12)	,%%ymm11		\n\t"\
		"vmovaps	    (%%rax)	,%%ymm0						\n\t		vmovaps	    (%%r10)	,%%ymm8 		\n\t"\
		"vmovaps	0x20(%%rax)	,%%ymm1						\n\t		vmovaps	0x20(%%r10)	,%%ymm9 		\n\t"\
		"vsubpd	%%ymm2,%%ymm4,%%ymm4						\n\t		vsubpd	%%ymm10,%%ymm12,%%ymm12		\n\t"\
		"vsubpd	%%ymm3,%%ymm5,%%ymm5						\n\t		vsubpd	%%ymm11,%%ymm13,%%ymm13		\n\t"\
		"vaddpd	%%ymm2,%%ymm2,%%ymm2						\n\t		vaddpd	%%ymm10,%%ymm10,%%ymm10		\n\t"\
		"vaddpd	%%ymm3,%%ymm3,%%ymm3						\n\t		vaddpd	%%ymm11,%%ymm11,%%ymm11		\n\t"\
		"vaddpd	%%ymm4,%%ymm2,%%ymm2						\n\t		vaddpd	%%ymm12,%%ymm10,%%ymm10		\n\t"\
		"vaddpd	%%ymm5,%%ymm3,%%ymm3						\n\t		vaddpd	%%ymm13,%%ymm11,%%ymm11		\n\t"\
		"vaddpd	%%ymm2,%%ymm0,%%ymm0						\n\t		vaddpd	%%ymm10,%%ymm8 ,%%ymm8 		\n\t"\
		"vaddpd	%%ymm3,%%ymm1,%%ymm1						\n\t		vaddpd	%%ymm11,%%ymm9 ,%%ymm9 		\n\t"\
		"vmovaps	%%ymm0		,    (%%rsi)				\n\t		vmovaps	%%ymm8 		,    (%%r13)	\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rsi)				\n\t		vmovaps	%%ymm9 		,0x20(%%r13)	\n\t"\
	" vfmadd132pd	%%ymm6,%%ymm0,%%ymm2 					\n\t	 vfmadd132pd	%%ymm6,%%ymm8,%%ymm10	\n\t"\
	" vfmadd132pd	%%ymm6,%%ymm1,%%ymm3 					\n\t	 vfmadd132pd	%%ymm6,%%ymm9,%%ymm11	\n\t"\
		"vmovaps	%%ymm2,%%ymm0							\n\t		vmovaps	%%ymm10,%%ymm8 				\n\t"\
		"vmovaps	%%ymm3,%%ymm1							\n\t		vmovaps	%%ymm11,%%ymm9 				\n\t"\
	"movq	%[__o3]	,%%rdi 									\n\t	movq	0x18(%%r9),%%r14				\n\t"/* o3 */\
	"movq	%[__o6]	,%%rsi 									\n\t	movq	0x30(%%r9),%%r13				\n\t"/* o6 */\
	"vfnmadd231pd	%%ymm7,%%ymm5,%%ymm0 					\n\t	vfnmadd231pd	%%ymm7,%%ymm13,%%ymm8 	\n\t"\
	" vfmadd231pd	%%ymm7,%%ymm4,%%ymm1 					\n\t	 vfmadd231pd	%%ymm7,%%ymm12,%%ymm9 	\n\t"\
	" vfmadd231pd	%%ymm7,%%ymm5,%%ymm2 					\n\t	 vfmadd231pd	%%ymm7,%%ymm13,%%ymm10	\n\t"\
	"vfnmadd231pd	%%ymm7,%%ymm4,%%ymm3 					\n\t	vfnmadd231pd	%%ymm7,%%ymm12,%%ymm11	\n\t"\
		"vmovaps	%%ymm0		,    (%%rsi)				\n\t		vmovaps	%%ymm8 		,    (%%r13)	\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rsi)				\n\t		vmovaps	%%ymm9 		,0x20(%%r13)	\n\t"\
		"vmovaps	%%ymm2		,    (%%rdi)				\n\t		vmovaps	%%ymm10		,    (%%r14)	\n\t"\
		"vmovaps	%%ymm3		,0x20(%%rdi)				\n\t		vmovaps	%%ymm11		,0x20(%%r14)	\n\t"\
	/* Block 2: */\
	"movq	%[__i1]	,%%rax									\n\t	movq	0x08(%%r8),%%r10				\n\t"/* i1 */\
	"movq	%[__i4]	,%%rbx									\n\t	movq	0x20(%%r8),%%r11				\n\t"/* i4 */\
	"movq	%[__i7]	,%%rcx									\n\t	movq	0x38(%%r8),%%r12				\n\t"/* i7 */\
		"vmovaps	    (%%rbx)	,%%ymm2						\n\t		vmovaps	    (%%r11)	,%%ymm10		\n\t"\
		"vmovaps	0x20(%%rbx)	,%%ymm3						\n\t		vmovaps	0x20(%%r11)	,%%ymm11		\n\t"\
		"vmovaps	    (%%rdx)	,%%ymm6						\n\t"/* c1 */\
		"vmovaps	0x20(%%rdx)	,%%ymm7						\n\t"\
		"vmovaps	%%ymm2		,%%ymm0						\n\t		vmovaps	%%ymm10		,%%ymm8 		\n\t"\
		"vmovaps	%%ymm3		,%%ymm1						\n\t		vmovaps	%%ymm11		,%%ymm9 		\n\t"\
		"vmulpd		%%ymm6,%%ymm2,%%ymm2					\n\t		vmulpd		%%ymm6,%%ymm10,%%ymm10	\n\t"\
		"vmulpd		%%ymm6,%%ymm3,%%ymm3					\n\t		vmulpd		%%ymm6,%%ymm11,%%ymm11	\n\t"\
	" vfmadd231pd	%%ymm7,%%ymm1,%%ymm2 					\n\t	 vfmadd231pd	%%ymm7,%%ymm9 ,%%ymm10	\n\t"\
	"vfnmadd231pd	%%ymm7,%%ymm0,%%ymm3 					\n\t	vfnmadd231pd	%%ymm7,%%ymm8 ,%%ymm11	\n\t"\
		"vmovaps	    (%%rcx)	,%%ymm4						\n\t		vmovaps	    (%%r12)	,%%ymm12		\n\t"\
		"vmovaps	0x20(%%rcx)	,%%ymm5						\n\t		vmovaps	0x20(%%r12)	,%%ymm13		\n\t"\
		"vmovaps	0x40(%%rdx)	,%%ymm6						\n\t"/* c2 */\
		"vmovaps	0x60(%%rdx)	,%%ymm7						\n\t"\
		"vmovaps	%%ymm4		,%%ymm0						\n\t		vmovaps	%%ymm12		,%%ymm8 		\n\t"\
		"vmovaps	%%ymm5		,%%ymm1						\n\t		vmovaps	%%ymm13		,%%ymm9 		\n\t"\
		"vmulpd		%%ymm6,%%ymm4,%%ymm4					\n\t		vmulpd		%%ymm6,%%ymm12,%%ymm12	\n\t"\
		"vmulpd		%%ymm6,%%ymm5,%%ymm5					\n\t		vmulpd		%%ymm6,%%ymm13,%%ymm13	\n\t"\
	" vfmadd231pd	%%ymm7,%%ymm1,%%ymm4 					\n\t	 vfmadd231pd	%%ymm7,%%ymm9 ,%%ymm12	\n\t"\
	"vfnmadd231pd	%%ymm7,%%ymm0,%%ymm5 					\n\t	vfnmadd231pd	%%ymm7,%%ymm8 ,%%ymm13	\n\t"\
		"vmovaps	0x80(%%rdx)	,%%ymm6						\n\t"/* c3m1 */\
		"vmovaps	0xa0(%%rdx)	,%%ymm7						\n\t"\
		"vmovaps	    (%%rax)	,%%ymm0						\n\t		vmovaps	    (%%r10)	,%%ymm8 		\n\t"\
		"vmovaps	0x20(%%rax)	,%%ymm1						\n\t		vmovaps	0x20(%%r10)	,%%ymm9 		\n\t"\
	"movq	%[__o7]	,%%rsi 									\n\t	movq	0x38(%%r9),%%r13				\n\t"/* o7 */\
		"vsubpd	%%ymm4,%%ymm2,%%ymm2						\n\t		vsubpd	%%ymm12,%%ymm10,%%ymm10		\n\t"\
		"vsubpd	%%ymm5,%%ymm3,%%ymm3						\n\t		vsubpd	%%ymm13,%%ymm11,%%ymm11		\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4						\n\t		vaddpd	%%ymm12,%%ymm12,%%ymm12		\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5						\n\t		vaddpd	%%ymm13,%%ymm13,%%ymm13		\n\t"\
		"vaddpd	%%ymm2,%%ymm4,%%ymm4						\n\t		vaddpd	%%ymm10,%%ymm12,%%ymm12		\n\t"\
		"vaddpd	%%ymm3,%%ymm5,%%ymm5						\n\t		vaddpd	%%ymm11,%%ymm13,%%ymm13		\n\t"\
		"vaddpd	%%ymm4,%%ymm0,%%ymm0						\n\t		vaddpd	%%ymm12,%%ymm8 ,%%ymm8 		\n\t"\
		"vaddpd	%%ymm5,%%ymm1,%%ymm1						\n\t		vaddpd	%%ymm13,%%ymm9 ,%%ymm9 		\n\t"\
		"vmovaps	%%ymm0		,    (%%rsi)				\n\t		vmovaps	%%ymm8 		,    (%%r13)	\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rsi)				\n\t		vmovaps	%%ymm9 		,0x20(%%r13)	\n\t"\
	" vfmadd132pd	%%ymm6,%%ymm0,%%ymm4 					\n\t	 vfmadd132pd	%%ymm6,%%ymm8,%%ymm12	\n\t"\
	" vfmadd132pd	%%ymm6,%%ymm1,%%ymm5 					\n\t	 vfmadd132pd	%%ymm6,%%ymm9,%%ymm13	\n\t"\
		"vmovaps	%%ymm4,%%ymm0							\n\t		vmovaps	%%ymm12,%%ymm8 				\n\t"\
		"vmovaps	%%ymm5,%%ymm1							\n\t		vmovaps	%%ymm13,%%ymm9 				\n\t"\
	"movq	%[__o1]	,%%rdi 									\n\t	movq	0x08(%%r9),%%r14				\n\t"/* o1 */\
	"movq	%[__o4]	,%%rsi 									\n\t	movq	0x20(%%r9),%%r13				\n\t"/* o4 */\
	"vfnmadd231pd	%%ymm7,%%ymm3,%%ymm0 					\n\t	vfnmadd231pd	%%ymm7,%%ymm11,%%ymm8 	\n\t"\
	" vfmadd231pd	%%ymm7,%%ymm2,%%ymm1 					\n\t	 vfmadd231pd	%%ymm7,%%ymm10,%%ymm9 	\n\t"\
	" vfmadd231pd	%%ymm7,%%ymm3,%%ymm4 					\n\t	 vfmadd231pd	%%ymm7,%%ymm11,%%ymm12	\n\t"\
	"vfnmadd231pd	%%ymm7,%%ymm2,%%ymm5 					\n\t	vfnmadd231pd	%%ymm7,%%ymm10,%%ymm13	\n\t"\
		"vmovaps	%%ymm0		,    (%%rsi)				\n\t		vmovaps	%%ymm8 		,    (%%r13)	\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rsi)				\n\t		vmovaps	%%ymm9 		,0x20(%%r13)	\n\t"\
		"vmovaps	%%ymm4		,    (%%rdi)				\n\t		vmovaps	%%ymm12		,    (%%r14)	\n\t"\
		"vmovaps	%%ymm5		,0x20(%%rdi)				\n\t		vmovaps	%%ymm13		,0x20(%%r14)	\n\t"\
	/* Block 3: */\
	"movq	%[__i2]	,%%rax									\n\t	movq	0x10(%%r8),%%r10				\n\t"/* i2 */\
	"movq	%[__i5]	,%%rbx									\n\t	movq	0x28(%%r8),%%r11				\n\t"/* i5 */\
	"movq	%[__i8]	,%%rcx									\n\t	movq	0x40(%%r8),%%r12				\n\t"/* i8 */\
		"vmovaps	    (%%rbx)	,%%ymm2						\n\t		vmovaps	    (%%r11)	,%%ymm10		\n\t"\
		"vmovaps	0x20(%%rbx)	,%%ymm3						\n\t		vmovaps	0x20(%%r11)	,%%ymm11		\n\t"\
		"vmovaps	0x40(%%rdx)	,%%ymm6						\n\t"/* c2 */\
		"vmovaps	0x60(%%rdx)	,%%ymm7						\n\t"\
		"vmovaps	%%ymm2		,%%ymm0						\n\t		vmovaps	%%ymm10		,%%ymm8 		\n\t"\
		"vmovaps	%%ymm3		,%%ymm1						\n\t		vmovaps	%%ymm11		,%%ymm9 		\n\t"\
		"vmulpd		%%ymm6,%%ymm2,%%ymm2					\n\t		vmulpd		%%ymm6,%%ymm10,%%ymm10	\n\t"\
		"vmulpd		%%ymm6,%%ymm3,%%ymm3					\n\t		vmulpd		%%ymm6,%%ymm11,%%ymm11	\n\t"\
	" vfmadd231pd	%%ymm7,%%ymm1,%%ymm2 					\n\t	 vfmadd231pd	%%ymm7,%%ymm9 ,%%ymm10	\n\t"\
	"vfnmadd231pd	%%ymm7,%%ymm0,%%ymm3 					\n\t	vfnmadd231pd	%%ymm7,%%ymm8 ,%%ymm11	\n\t"\
		"vmovaps	    (%%rcx)	,%%ymm4						\n\t		vmovaps	    (%%r12)	,%%ymm12		\n\t"\
		"vmovaps	0x20(%%rcx)	,%%ymm5						\n\t		vmovaps	0x20(%%r12)	,%%ymm13		\n\t"\
		"vmovaps	0xc0(%%rdx)	,%%ymm6						\n\t"/* c4 */\
		"vmovaps	0xe0(%%rdx)	,%%ymm7						\n\t"\
		"vmovaps	%%ymm4		,%%ymm0						\n\t		vmovaps	%%ymm12		,%%ymm8 		\n\t"\
		"vmovaps	%%ymm5		,%%ymm1						\n\t		vmovaps	%%ymm13		,%%ymm9 		\n\t"\
		"vmulpd		%%ymm6,%%ymm4,%%ymm4					\n\t		vmulpd		%%ymm6,%%ymm12,%%ymm12	\n\t"\
		"vmulpd		%%ymm6,%%ymm5,%%ymm5					\n\t		vmulpd		%%ymm6,%%ymm13,%%ymm13	\n\t"\
	" vfmadd231pd	%%ymm7,%%ymm1,%%ymm4 					\n\t	 vfmadd231pd	%%ymm7,%%ymm9 ,%%ymm12	\n\t"\
	"vfnmadd231pd	%%ymm7,%%ymm0,%%ymm5 					\n\t	vfnmadd231pd	%%ymm7,%%ymm8 ,%%ymm13	\n\t"\
		"vmovaps	0x80(%%rdx)	,%%ymm6						\n\t"/* c3m1 */\
		"vmovaps	0xa0(%%rdx)	,%%ymm7						\n\t"\
		"vmovaps	    (%%rax)	,%%ymm0						\n\t		vmovaps	    (%%r10)	,%%ymm8 		\n\t"\
		"vmovaps	0x20(%%rax)	,%%ymm1						\n\t		vmovaps	0x20(%%r10)	,%%ymm9 		\n\t"\
	"movq	%[__o5]	,%%rsi 									\n\t	movq	0x28(%%r9),%%r13				\n\t"/* o5 */\
		"vsubpd	%%ymm4,%%ymm2,%%ymm2						\n\t		vsubpd	%%ymm12,%%ymm10,%%ymm10		\n\t"\
		"vsubpd	%%ymm5,%%ymm3,%%ymm3						\n\t		vsubpd	%%ymm13,%%ymm11,%%ymm11		\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4						\n\t		vaddpd	%%ymm12,%%ymm12,%%ymm12		\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5						\n\t		vaddpd	%%ymm13,%%ymm13,%%ymm13		\n\t"\
		"vaddpd	%%ymm2,%%ymm4,%%ymm4						\n\t		vaddpd	%%ymm10,%%ymm12,%%ymm12		\n\t"\
		"vaddpd	%%ymm3,%%ymm5,%%ymm5						\n\t		vaddpd	%%ymm11,%%ymm13,%%ymm13		\n\t"\
		"vaddpd	%%ymm4,%%ymm0,%%ymm0						\n\t		vaddpd	%%ymm12,%%ymm8 ,%%ymm8 		\n\t"\
		"vaddpd	%%ymm5,%%ymm1,%%ymm1						\n\t		vaddpd	%%ymm13,%%ymm9 ,%%ymm9 		\n\t"\
		"vmovaps	%%ymm0		,    (%%rsi)				\n\t		vmovaps	%%ymm8 		,    (%%r13)	\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rsi)				\n\t		vmovaps	%%ymm9 		,0x20(%%r13)	\n\t"\
	" vfmadd132pd	%%ymm6,%%ymm0,%%ymm4 					\n\t	 vfmadd132pd	%%ymm6,%%ymm8,%%ymm12	\n\t"\
	" vfmadd132pd	%%ymm6,%%ymm1,%%ymm5 					\n\t	 vfmadd132pd	%%ymm6,%%ymm9,%%ymm13	\n\t"\
		"vmovaps	%%ymm4,%%ymm0							\n\t		vmovaps	%%ymm12,%%ymm8 				\n\t"\
		"vmovaps	%%ymm5,%%ymm1							\n\t		vmovaps	%%ymm13,%%ymm9 				\n\t"\
	"movq	%[__o8]	,%%rdi 									\n\t	movq	0x40(%%r9),%%r14				\n\t"/* o8 */\
	"movq	%[__o2]	,%%rsi 									\n\t	movq	0x10(%%r9),%%r13				\n\t"/* o2 */\
	"vfnmadd231pd	%%ymm7,%%ymm3,%%ymm0 					\n\t	vfnmadd231pd	%%ymm7,%%ymm11,%%ymm8 	\n\t"\
	" vfmadd231pd	%%ymm7,%%ymm2,%%ymm1 					\n\t	 vfmadd231pd	%%ymm7,%%ymm10,%%ymm9 	\n\t"\
	" vfmadd231pd	%%ymm7,%%ymm3,%%ymm4 					\n\t	 vfmadd231pd	%%ymm7,%%ymm11,%%ymm12	\n\t"\
	"vfnmadd231pd	%%ymm7,%%ymm2,%%ymm5 					\n\t	vfnmadd231pd	%%ymm7,%%ymm10,%%ymm13	\n\t"\
		"vmovaps	%%ymm0		,    (%%rsi)				\n\t		vmovaps	%%ymm8 		,    (%%r13)	\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rsi)				\n\t		vmovaps	%%ymm9 		,0x20(%%r13)	\n\t"\
		"vmovaps	%%ymm4		,    (%%rdi)				\n\t		vmovaps	%%ymm12		,    (%%r14)	\n\t"\
		"vmovaps	%%ymm5		,0x20(%%rdi)				\n\t		vmovaps	%%ymm13		,0x20(%%r14)	\n\t"\
		:					/* outputs: none */\
		: [__i0] "m" (Xi0)	/* All inputs from memory addresses here */\
		 ,[__i1] "m" (Xi1)\
		 ,[__i2] "m" (Xi2)\
		 ,[__i3] "m" (Xi3)\
		 ,[__i4] "m" (Xi4)\
		 ,[__i5] "m" (Xi5)\
		 ,[__i6] "m" (Xi6)\
		 ,[__i7] "m" (Xi7)\
		 ,[__i8] "m" (Xi8)\
		 ,[__c1] "m" (Xc1)\
		 ,[__two] "m" (Xtwo)\
		 ,[__o0] "m" (Xo0)\
		 ,[__o1] "m" (Xo1)\
		 ,[__o2] "m" (Xo2)\
		 ,[__o3] "m" (Xo3)\
		 ,[__o4] "m" (Xo4)\
		 ,[__o5] "m" (Xo5)\
		 ,[__o6] "m" (Xo6)\
		 ,[__o7] "m" (Xo7)\
		 ,[__o8] "m" (Xo8)\
		 ,[__iptr] "m" (Xiptr)	/* Due to GCC macro argc limit of 30, to enable 16-register data-doubled version */\
		 ,[__optr] "m" (Xoptr)	/* of the radix-9 macros use 2 length-9 ptr arrays for second set of IOs. */\
		: "cc","memory","rax","rbx","rcx","rdx","rsi","rdi","r8","r9","r10","r11","r12","r13","r14","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
		);\
	}

#endif	// AVX512/AVX2 ?

#if defined(USE_AVX) && !defined(USE_AVX512)	// Need to def these for both AVX and AVX2 builds:

	// Simple 64-bit-ified version of the same-named 32-bit ASM macros, using just xmm0-7:

	/*...Radix-9 DIF: Ins in memory locations __i0-8.\
		Outs at memory locations __o0-8, assumed disjoint with inputs:\
	*/\
	#define SSE2_RADIX_09_DIF(Xi0,Xi1,Xi2,Xi3,Xi4,Xi5,Xi6,Xi7,Xi8, Xc1, Xo0,Xo1,Xo2,Xo3,Xo4,Xo5,Xo6,Xo7,Xo8)\
	{\
	__asm__ volatile (\
	"movq	%[__i0]	,%%rax 	\n\t"/* __i0-8; e[abc]x store input addresses */\
	"movq	%[__o0]	,%%rsi	\n\t"\
	"movq	%[__c1]	,%%rdx 	\n\t"/* edx stores trig addresses throughout */\
	/* Block 1: */\
		"movq	%%rax		,%%rbx\n\t"\
		"movq	%%rax		,%%rcx\n\t"\
		"movq	%[__i3]	,%%rbx 	\n\t"\
		"movq	%[__i6]	,%%rcx 	\n\t"\
		"vmovaps	    (%%rbx)	,%%ymm2\n\t"\
		"vmovaps	0x20(%%rbx)	,%%ymm3\n\t"\
		"vmovaps	    (%%rax)	,%%ymm0\n\t"\
		"vmovaps	0x20(%%rax)	,%%ymm1\n\t"\
		"vmovaps	    (%%rcx)	,%%ymm6\n\t"\
		"vmovaps	0x20(%%rcx)	,%%ymm7\n\t"\
		"vmovaps	%%ymm2		,%%ymm4\n\t"\
		"vmovaps	%%ymm3		,%%ymm5\n\t"\
		"vaddpd	%%ymm6,%%ymm2,%%ymm2\n\t"\
		"vaddpd	%%ymm7,%%ymm3,%%ymm3\n\t"\
		"vsubpd	%%ymm6,%%ymm4,%%ymm4\n\t"\
		"vsubpd	%%ymm7,%%ymm5,%%ymm5\n\t"\
		"vaddpd	%%ymm2,%%ymm0,%%ymm0\n\t"\
		"vaddpd	%%ymm3,%%ymm1,%%ymm1\n\t"\
		"vmovaps	0x80(%%rdx)	,%%ymm6\n\t"/* c3m1 */\
		"vmovaps	0xa0(%%rdx)	,%%ymm7\n\t"\
		"vmovaps	%%ymm0		,    (%%rsi)	/* <- t00 */\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rsi)	/* <- t01 */\n\t"\
	"movq	%[__o1]	,%%rdi	\n\t"\
	"movq	%[__o2]	,%%rsi	\n\t"\
		"vmulpd	%%ymm6,%%ymm2,%%ymm2\n\t"\
		"vmulpd	%%ymm6,%%ymm3,%%ymm3\n\t"\
		"vmulpd	%%ymm7,%%ymm4,%%ymm4\n\t"\
		"vmulpd	%%ymm7,%%ymm5,%%ymm5\n\t"\
		"vaddpd	%%ymm0,%%ymm2,%%ymm2\n\t"\
		"vaddpd	%%ymm1,%%ymm3,%%ymm3\n\t"\
		"vmovaps	%%ymm2		,%%ymm0\n\t"\
		"vmovaps	%%ymm3		,%%ymm1\n\t"\
		"vsubpd	%%ymm5,%%ymm2,%%ymm2\n\t"\
		"vaddpd	%%ymm4,%%ymm3,%%ymm3\n\t"\
		"vaddpd	%%ymm5,%%ymm0,%%ymm0\n\t"\
		"vsubpd	%%ymm4,%%ymm1,%%ymm1\n\t"\
		"vmovaps	%%ymm2		,    (%%rdi)	/* <- t02 */\n\t"\
		"vmovaps	%%ymm3		,0x20(%%rdi)	/* <- t03 */\n\t"\
		"vmovaps	%%ymm0		,    (%%rsi)	/* <- t04 */\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rsi)	/* <- t05 */\n\t"\
	/* Block 2: */\
		"movq	%[__i1]	,%%rax 	\n\t"\
		"movq	%[__i4]	,%%rbx 	\n\t"\
		"movq	%[__i7]	,%%rcx 	\n\t"\
		"vmovaps	    (%%rbx)	,%%ymm4\n\t"\
		"vmovaps	0x20(%%rbx)	,%%ymm5\n\t"\
		"vmovaps	    (%%rax)	,%%ymm0\n\t"\
		"vmovaps	0x20(%%rax)	,%%ymm1\n\t"\
		"vmovaps	    (%%rcx)	,%%ymm2\n\t"\
		"vmovaps	0x20(%%rcx)	,%%ymm3\n\t"\
	"movq	%[__o3]	,%%rdi	\n\t"\
		"vsubpd	%%ymm2,%%ymm4,%%ymm4\n\t"\
		"vsubpd	%%ymm3,%%ymm5,%%ymm5\n\t"\
		"vaddpd	%%ymm2,%%ymm2,%%ymm2\n\t"\
		"vaddpd	%%ymm3,%%ymm3,%%ymm3\n\t"\
		"vaddpd	%%ymm4,%%ymm2,%%ymm2\n\t"\
		"vaddpd	%%ymm5,%%ymm3,%%ymm3\n\t"\
		"vaddpd	%%ymm2,%%ymm0,%%ymm0\n\t"\
		"vaddpd	%%ymm3,%%ymm1,%%ymm1\n\t"\
		"vmovaps	%%ymm0		,    (%%rdi)	/* <- t06 */\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rdi)	/* <- t07 */\n\t"\
	"movq	%[__o4]	,%%rsi	\n\t"\
	"movq	%[__o5]	,%%rdi	\n\t"\
		"vmulpd	%%ymm6,%%ymm2,%%ymm2\n\t"\
		"vmulpd	%%ymm6,%%ymm3,%%ymm3\n\t"\
		"vmulpd	%%ymm7,%%ymm4,%%ymm4\n\t"\
		"vmulpd	%%ymm7,%%ymm5,%%ymm5\n\t"\
		"vaddpd	%%ymm0,%%ymm2,%%ymm2\n\t"\
		"vaddpd	%%ymm1,%%ymm3,%%ymm3\n\t"\
		"vmovaps	%%ymm2		,%%ymm0\n\t"\
		"vmovaps	%%ymm3		,%%ymm1\n\t"\
		"vsubpd	%%ymm5,%%ymm2,%%ymm2\n\t"\
		"vaddpd	%%ymm4,%%ymm3,%%ymm3\n\t"\
		"vaddpd	%%ymm5,%%ymm0,%%ymm0\n\t"\
		"vsubpd	%%ymm4,%%ymm1,%%ymm1\n\t"\
		"vmovaps	%%ymm2		,    (%%rsi)	/* <- t08 */\n\t"\
		"vmovaps	%%ymm3		,0x20(%%rsi)	/* <- t09 */\n\t"\
		"vmovaps	%%ymm0		,    (%%rdi)	/* <- t0a */\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rdi)	/* <- t0b */\n\t"\
	/* Block 3: */\
		"movq	%[__i2]	,%%rax 	\n\t"\
		"movq	%[__i5]	,%%rbx 	\n\t"\
		"movq	%[__i8]	,%%rcx 	\n\t"\
		"vmovaps	    (%%rbx)	,%%ymm4\n\t"\
		"vmovaps	0x20(%%rbx)	,%%ymm5\n\t"\
		"vmovaps	    (%%rax)	,%%ymm0\n\t"\
		"vmovaps	0x20(%%rax)	,%%ymm1\n\t"\
		"vmovaps	    (%%rcx)	,%%ymm2\n\t"\
		"vmovaps	0x20(%%rcx)	,%%ymm3\n\t"\
	"movq	%[__o6]	,%%rsi	\n\t"\
		"vsubpd	%%ymm2,%%ymm4,%%ymm4\n\t"\
		"vsubpd	%%ymm3,%%ymm5,%%ymm5\n\t"\
		"vaddpd	%%ymm2,%%ymm2,%%ymm2\n\t"\
		"vaddpd	%%ymm3,%%ymm3,%%ymm3\n\t"\
		"vaddpd	%%ymm4,%%ymm2,%%ymm2\n\t"\
		"vaddpd	%%ymm5,%%ymm3,%%ymm3\n\t"\
		"vaddpd	%%ymm2,%%ymm0,%%ymm0\n\t"\
		"vaddpd	%%ymm3,%%ymm1,%%ymm1\n\t"\
		"vmovaps	%%ymm0		,    (%%rsi)	/* <- t0c */\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rsi)	/* <- t0d */\n\t"\
	"movq	%[__o7]	,%%rdi	\n\t"\
	"movq	%[__o8]	,%%rsi	\n\t"\
		"vmulpd	%%ymm6,%%ymm2,%%ymm2\n\t"\
		"vmulpd	%%ymm6,%%ymm3,%%ymm3\n\t"\
		"vmulpd	%%ymm7,%%ymm4,%%ymm4\n\t"\
		"vmulpd	%%ymm7,%%ymm5,%%ymm5\n\t"\
		"vaddpd	%%ymm0,%%ymm2,%%ymm2\n\t"\
		"vaddpd	%%ymm1,%%ymm3,%%ymm3\n\t"\
		"vmovaps	%%ymm2		,%%ymm0\n\t"\
		"vmovaps	%%ymm3		,%%ymm1\n\t"\
		"vsubpd	%%ymm5,%%ymm2,%%ymm2\n\t"\
		"vaddpd	%%ymm4,%%ymm3,%%ymm3\n\t"\
		"vaddpd	%%ymm5,%%ymm0,%%ymm0\n\t"\
		"vsubpd	%%ymm4,%%ymm1,%%ymm1\n\t"\
		"vmovaps	%%ymm2		,    (%%rdi)	/* <- t0e */\n\t"\
		"vmovaps	%%ymm3		,0x20(%%rdi)	/* <- t0f */\n\t"\
		"vmovaps	%%ymm0		,    (%%rsi)	/* <- t0g */\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rsi)	/* <- t0h */\n\t"\
	/*****************************************/\
	/* now do three more radix-3 transforms: */\
	/*****************************************/\
	/* Block 1: */\
	"movq	%[__o0]	,%%rax	\n\t"\
	"movq	%[__o3]	,%%rbx	\n\t"\
	"movq	%[__o6]	,%%rcx	\n\t"\
		"vmovaps	    (%%rbx)	,%%ymm4\n\t"\
		"vmovaps	0x20(%%rbx)	,%%ymm5\n\t"\
		"vmovaps	    (%%rcx)	,%%ymm2\n\t"\
		"vmovaps	0x20(%%rcx)	,%%ymm3\n\t"\
		"vmovaps	    (%%rax)	,%%ymm0\n\t"\
		"vmovaps	0x20(%%rax)	,%%ymm1\n\t"\
		"vsubpd	%%ymm2,%%ymm4,%%ymm4\n\t"\
		"vsubpd	%%ymm3,%%ymm5,%%ymm5\n\t"\
		"vaddpd	%%ymm2,%%ymm2,%%ymm2\n\t"\
		"vaddpd	%%ymm3,%%ymm3,%%ymm3\n\t"\
		"vaddpd	%%ymm4,%%ymm2,%%ymm2\n\t"\
		"vaddpd	%%ymm5,%%ymm3,%%ymm3\n\t"\
		"vaddpd	%%ymm2,%%ymm0,%%ymm0\n\t"\
		"vaddpd	%%ymm3,%%ymm1,%%ymm1\n\t"\
		"vmovaps	%%ymm0		,    (%%rax)\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rax)\n\t"\
		"vmulpd	%%ymm6,%%ymm2,%%ymm2\n\t"\
		"vmulpd	%%ymm6,%%ymm3,%%ymm3\n\t"\
		"vmulpd	%%ymm7,%%ymm4,%%ymm4\n\t"\
		"vmulpd	%%ymm7,%%ymm5,%%ymm5\n\t"\
		"vaddpd	%%ymm0,%%ymm2,%%ymm2\n\t"\
		"vaddpd	%%ymm1,%%ymm3,%%ymm3\n\t"\
		"vmovaps	%%ymm2		,%%ymm0\n\t"\
		"vmovaps	%%ymm3		,%%ymm1\n\t"\
		"vsubpd	%%ymm5,%%ymm2,%%ymm2\n\t"\
		"vaddpd	%%ymm4,%%ymm3,%%ymm3\n\t"\
		"vaddpd	%%ymm5,%%ymm0,%%ymm0\n\t"\
		"vsubpd	%%ymm4,%%ymm1,%%ymm1\n\t"\
		"vmovaps	%%ymm2		,    (%%rbx)\n\t"\
		"vmovaps	%%ymm3		,0x20(%%rbx)\n\t"\
		"vmovaps	%%ymm0		,    (%%rcx)\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rcx)\n\t"\
	/* Block 2: */\
	"movq	%[__o1]	,%%rax	\n\t"\
	"movq	%[__o4]	,%%rbx	\n\t"\
	"movq	%[__o7]	,%%rcx	\n\t"\
		"vmovaps	    (%%rbx)	,%%ymm2\n\t"\
		"vmovaps	0x20(%%rbx)	,%%ymm3\n\t"\
		"vmovaps	    (%%rdx)	,%%ymm6\n\t"/* c1 */\
		"vmovaps	0x20(%%rdx)	,%%ymm7\n\t"\
		"vmovaps	%%ymm2		,%%ymm0\n\t"\
		"vmovaps	%%ymm3		,%%ymm1\n\t"\
		"vmulpd	%%ymm6,%%ymm2,%%ymm2\n\t"\
		"vmulpd	%%ymm6,%%ymm3,%%ymm3\n\t"\
		"vmulpd	%%ymm7,%%ymm0,%%ymm0\n\t"\
		"vmulpd	%%ymm7,%%ymm1,%%ymm1\n\t"\
		"vsubpd	%%ymm1,%%ymm2,%%ymm2\n\t"\
		"vaddpd	%%ymm0,%%ymm3,%%ymm3\n\t"\
		"vmovaps	    (%%rcx)	,%%ymm4\n\t"\
		"vmovaps	0x20(%%rcx)	,%%ymm5\n\t"\
		"vmovaps	0x40(%%rdx)	,%%ymm6\n\t"/* c2 */\
		"vmovaps	0x60(%%rdx)	,%%ymm7\n\t"\
		"vmovaps	%%ymm4		,%%ymm0\n\t"\
		"vmovaps	%%ymm5		,%%ymm1\n\t"\
		"vmulpd	%%ymm6,%%ymm4,%%ymm4\n\t"\
		"vmulpd	%%ymm6,%%ymm5,%%ymm5\n\t"\
		"vmulpd	%%ymm7,%%ymm0,%%ymm0\n\t"\
		"vmulpd	%%ymm7,%%ymm1,%%ymm1\n\t"\
		"vsubpd	%%ymm1,%%ymm4,%%ymm4\n\t"\
		"vaddpd	%%ymm0,%%ymm5,%%ymm5\n\t"\
		"vmovaps	0x80(%%rdx)	,%%ymm6\n\t"/* c3m1 */\
		"vmovaps	0xa0(%%rdx)	,%%ymm7\n\t"\
		"vmovaps	    (%%rax)	,%%ymm0\n\t"\
		"vmovaps	0x20(%%rax)	,%%ymm1\n\t"\
		"vsubpd	%%ymm4,%%ymm2,%%ymm2\n\t"\
		"vsubpd	%%ymm5,%%ymm3,%%ymm3\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5\n\t"\
		"vaddpd	%%ymm2,%%ymm4,%%ymm4\n\t"\
		"vaddpd	%%ymm3,%%ymm5,%%ymm5\n\t"\
		"vaddpd	%%ymm4,%%ymm0,%%ymm0\n\t"\
		"vaddpd	%%ymm5,%%ymm1,%%ymm1\n\t"\
		"vmovaps	%%ymm0		,    (%%rax)\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rax)\n\t"\
		"vmulpd	%%ymm6,%%ymm4,%%ymm4\n\t"\
		"vmulpd	%%ymm6,%%ymm5,%%ymm5\n\t"\
		"vmulpd	%%ymm7,%%ymm2,%%ymm2\n\t"\
		"vmulpd	%%ymm7,%%ymm3,%%ymm3\n\t"\
		"vaddpd	%%ymm0,%%ymm4,%%ymm4\n\t"\
		"vaddpd	%%ymm1,%%ymm5,%%ymm5\n\t"\
		"vmovaps	%%ymm4		,%%ymm0\n\t"\
		"vmovaps	%%ymm5		,%%ymm1\n\t"\
		"vsubpd	%%ymm3,%%ymm4,%%ymm4\n\t"\
		"vaddpd	%%ymm2,%%ymm5,%%ymm5\n\t"\
		"vaddpd	%%ymm3,%%ymm0,%%ymm0\n\t"\
		"vsubpd	%%ymm2,%%ymm1,%%ymm1\n\t"\
		"vmovaps	%%ymm4		,    (%%rbx)\n\t"\
		"vmovaps	%%ymm5		,0x20(%%rbx)\n\t"\
		"vmovaps	%%ymm0		,    (%%rcx)\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rcx)\n\t"\
	/* Block 3: */\
	"movq	%[__o2]	,%%rax	\n\t"\
	"movq	%[__o5]	,%%rbx	\n\t"\
	"movq	%[__o8]	,%%rcx	\n\t"\
		"vmovaps	    (%%rbx)	,%%ymm2\n\t"\
		"vmovaps	0x20(%%rbx)	,%%ymm3\n\t"\
		"vmovaps	0x40(%%rdx)	,%%ymm6\n\t"/* c2 */\
		"vmovaps	0x60(%%rdx)	,%%ymm7\n\t"\
		"vmovaps	%%ymm2		,%%ymm0\n\t"\
		"vmovaps	%%ymm3		,%%ymm1\n\t"\
		"vmulpd	%%ymm6,%%ymm2,%%ymm2\n\t"\
		"vmulpd	%%ymm6,%%ymm3,%%ymm3\n\t"\
		"vmulpd	%%ymm7,%%ymm0,%%ymm0\n\t"\
		"vmulpd	%%ymm7,%%ymm1,%%ymm1\n\t"\
		"vsubpd	%%ymm1,%%ymm2,%%ymm2\n\t"\
		"vaddpd	%%ymm0,%%ymm3,%%ymm3\n\t"\
		"vmovaps	    (%%rcx)	,%%ymm4\n\t"\
		"vmovaps	0x20(%%rcx)	,%%ymm5\n\t"\
		"vmovaps	0xc0(%%rdx)	,%%ymm6\n\t"/* c4 */\
		"vmovaps	0xe0(%%rdx)	,%%ymm7\n\t"\
		"vmovaps	%%ymm4		,%%ymm0\n\t"\
		"vmovaps	%%ymm5		,%%ymm1\n\t"\
		"vmulpd	%%ymm6,%%ymm4,%%ymm4\n\t"\
		"vmulpd	%%ymm6,%%ymm5,%%ymm5\n\t"\
		"vmulpd	%%ymm7,%%ymm0,%%ymm0\n\t"\
		"vmulpd	%%ymm7,%%ymm1,%%ymm1\n\t"\
		"vsubpd	%%ymm1,%%ymm4,%%ymm4\n\t"\
		"vaddpd	%%ymm0,%%ymm5,%%ymm5\n\t"\
		"vmovaps	0x80(%%rdx)	,%%ymm6\n\t"/* c3m1 */\
		"vmovaps	0xa0(%%rdx)	,%%ymm7\n\t"\
		"vmovaps	    (%%rax)	,%%ymm0\n\t"\
		"vmovaps	0x20(%%rax)	,%%ymm1\n\t"\
		"\n\t"\
		"vsubpd	%%ymm4,%%ymm2,%%ymm2\n\t"\
		"vsubpd	%%ymm5,%%ymm3,%%ymm3\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5\n\t"\
		"vaddpd	%%ymm2,%%ymm4,%%ymm4\n\t"\
		"vaddpd	%%ymm3,%%ymm5,%%ymm5\n\t"\
		"vaddpd	%%ymm4,%%ymm0,%%ymm0\n\t"\
		"vaddpd	%%ymm5,%%ymm1,%%ymm1\n\t"\
		"vmovaps	%%ymm0		,    (%%rax)\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rax)\n\t"\
		"vmulpd	%%ymm6,%%ymm4,%%ymm4\n\t"\
		"vmulpd	%%ymm6,%%ymm5,%%ymm5\n\t"\
		"vmulpd	%%ymm7,%%ymm2,%%ymm2\n\t"\
		"vmulpd	%%ymm7,%%ymm3,%%ymm3\n\t"\
		"vaddpd	%%ymm0,%%ymm4,%%ymm4\n\t"\
		"vaddpd	%%ymm1,%%ymm5,%%ymm5\n\t"\
		"\n\t"\
		"vmovaps	%%ymm4		,%%ymm0\n\t"\
		"vmovaps	%%ymm5		,%%ymm1\n\t"\
		"vsubpd	%%ymm3,%%ymm4,%%ymm4\n\t"\
		"vaddpd	%%ymm2,%%ymm5,%%ymm5\n\t"\
		"vaddpd	%%ymm3,%%ymm0,%%ymm0\n\t"\
		"vsubpd	%%ymm2,%%ymm1,%%ymm1\n\t"\
		"vmovaps	%%ymm4		,    (%%rbx)\n\t"\
		"vmovaps	%%ymm5		,0x20(%%rbx)\n\t"\
		"vmovaps	%%ymm0		,    (%%rcx)\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rcx)\n\t"\
		:					/* outputs: none */\
		: [__i0] "m" (Xi0)	/* All inputs from memory addresses here */\
		 ,[__i1] "m" (Xi1)\
		 ,[__i2] "m" (Xi2)\
		 ,[__i3] "m" (Xi3)\
		 ,[__i4] "m" (Xi4)\
		 ,[__i5] "m" (Xi5)\
		 ,[__i6] "m" (Xi6)\
		 ,[__i7] "m" (Xi7)\
		 ,[__i8] "m" (Xi8)\
		 ,[__c1] "m" (Xc1)\
		 ,[__o0] "m" (Xo0)\
		 ,[__o1] "m" (Xo1)\
		 ,[__o2] "m" (Xo2)\
		 ,[__o3] "m" (Xo3)\
		 ,[__o4] "m" (Xo4)\
		 ,[__o5] "m" (Xo5)\
		 ,[__o6] "m" (Xo6)\
		 ,[__o7] "m" (Xo7)\
		 ,[__o8] "m" (Xo8)\
		: "cc","memory","rax","rbx","rcx","rdx","rsi","rdi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
		);\
	}

	/*...Radix-9 DIT: Ins in memory locations __i0-8.\
		Outs at memory locations __o0-8, assumed disjoint with inputs:\
	*/\
	#define SSE2_RADIX_09_DIT(Xi0,Xi1,Xi2,Xi3,Xi4,Xi5,Xi6,Xi7,Xi8, Xc1, Xo0,Xo1,Xo2,Xo3,Xo4,Xo5,Xo6,Xo7,Xo8)\
	{\
	__asm__ volatile (\
	"movq	%[__i0]	,%%rax\n\t"\
	"movq	%[__i1]	,%%rbx\n\t"\
	"movq	%[__i2]	,%%rcx\n\t"\
	"movq	%[__c1]	,%%rdx 	\n\t"/* edx stores trig addresses throughout */\
	/* Block 1: */\
		"vmovaps	    (%%rbx)	,%%ymm2\n\t"\
		"vmovaps	0x20(%%rbx)	,%%ymm3\n\t"\
		"vmovaps	    (%%rax)	,%%ymm0\n\t"\
		"vmovaps	0x20(%%rax)	,%%ymm1\n\t"\
		"vmovaps	    (%%rcx)	,%%ymm6\n\t"\
		"vmovaps	0x20(%%rcx)	,%%ymm7\n\t"\
		"vmovaps	%%ymm2		,%%ymm4\n\t"\
		"vmovaps	%%ymm3		,%%ymm5\n\t"\
		"vaddpd	%%ymm6,%%ymm2,%%ymm2\n\t"\
		"vaddpd	%%ymm7,%%ymm3,%%ymm3\n\t"\
		"vsubpd	%%ymm6,%%ymm4,%%ymm4\n\t"\
		"vsubpd	%%ymm7,%%ymm5,%%ymm5\n\t"\
		"vaddpd	%%ymm2,%%ymm0,%%ymm0\n\t"\
		"vaddpd	%%ymm3,%%ymm1,%%ymm1\n\t"\
		"vmovaps	0x80(%%rdx)	,%%ymm6\n\t"/* c3m1 */\
		"vmovaps	0xa0(%%rdx)	,%%ymm7\n\t"\
		"vmovaps	%%ymm0		,    (%%rax)\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rax)\n\t"\
		"vmulpd	%%ymm6,%%ymm2,%%ymm2\n\t"\
		"vmulpd	%%ymm6,%%ymm3,%%ymm3\n\t"\
		"vmulpd	%%ymm7,%%ymm4,%%ymm4\n\t"\
		"vmulpd	%%ymm7,%%ymm5,%%ymm5\n\t"\
		"vaddpd	%%ymm0,%%ymm2,%%ymm2\n\t"\
		"vaddpd	%%ymm1,%%ymm3,%%ymm3\n\t"\
		"vmovaps	%%ymm2		,%%ymm0\n\t"\
		"vmovaps	%%ymm3		,%%ymm1\n\t"\
		"vaddpd	%%ymm5,%%ymm2,%%ymm2\n\t"\
		"vsubpd	%%ymm4,%%ymm3,%%ymm3\n\t"\
		"vsubpd	%%ymm5,%%ymm0,%%ymm0\n\t"\
		"vaddpd	%%ymm4,%%ymm1,%%ymm1\n\t"\
		"vmovaps	%%ymm2		,    (%%rbx)\n\t"\
		"vmovaps	%%ymm3		,0x20(%%rbx)\n\t"\
		"vmovaps	%%ymm0		,    (%%rcx)\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rcx)\n\t"\
	/* Block 2: */\
	"movq	%[__i3]	,%%rax\n\t"\
	"movq	%[__i4]	,%%rbx\n\t"\
	"movq	%[__i5]	,%%rcx\n\t"\
		"vmovaps	    (%%rbx)	,%%ymm4\n\t"\
		"vmovaps	0x20(%%rbx)	,%%ymm5\n\t"\
		"vmovaps	    (%%rcx)	,%%ymm2\n\t"\
		"vmovaps	0x20(%%rcx)	,%%ymm3\n\t"\
		"vmovaps	    (%%rax)	,%%ymm0\n\t"\
		"vmovaps	0x20(%%rax)	,%%ymm1\n\t"\
		"vsubpd	%%ymm2,%%ymm4,%%ymm4\n\t"\
		"vsubpd	%%ymm3,%%ymm5,%%ymm5\n\t"\
		"vaddpd	%%ymm2,%%ymm2,%%ymm2\n\t"\
		"vaddpd	%%ymm3,%%ymm3,%%ymm3\n\t"\
		"vaddpd	%%ymm4,%%ymm2,%%ymm2\n\t"\
		"vaddpd	%%ymm5,%%ymm3,%%ymm3\n\t"\
		"vaddpd	%%ymm2,%%ymm0,%%ymm0\n\t"\
		"vaddpd	%%ymm3,%%ymm1,%%ymm1\n\t"\
		"vmovaps	%%ymm0		,    (%%rax)\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rax)\n\t"\
		"vmulpd	%%ymm6,%%ymm2,%%ymm2\n\t"\
		"vmulpd	%%ymm6,%%ymm3,%%ymm3\n\t"\
		"vmulpd	%%ymm7,%%ymm4,%%ymm4\n\t"\
		"vmulpd	%%ymm7,%%ymm5,%%ymm5\n\t"\
		"vaddpd	%%ymm0,%%ymm2,%%ymm2\n\t"\
		"vaddpd	%%ymm1,%%ymm3,%%ymm3\n\t"\
		"vmovaps	%%ymm2		,%%ymm0\n\t"\
		"vmovaps	%%ymm3		,%%ymm1\n\t"\
		"vaddpd	%%ymm5,%%ymm2,%%ymm2\n\t"\
		"vsubpd	%%ymm4,%%ymm3,%%ymm3\n\t"\
		"vsubpd	%%ymm5,%%ymm0,%%ymm0\n\t"\
		"vaddpd	%%ymm4,%%ymm1,%%ymm1\n\t"\
		"vmovaps	%%ymm2		,    (%%rbx)\n\t"\
		"vmovaps	%%ymm3		,0x20(%%rbx)\n\t"\
		"vmovaps	%%ymm0		,    (%%rcx)\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rcx)\n\t"\
	/* Block 3: */\
	"movq	%[__i6]	,%%rax\n\t"\
	"movq	%[__i7]	,%%rbx\n\t"\
	"movq	%[__i8]	,%%rcx\n\t"\
		"vmovaps	    (%%rbx)	,%%ymm4\n\t"\
		"vmovaps	0x20(%%rbx)	,%%ymm5\n\t"\
		"vmovaps	    (%%rcx)	,%%ymm2\n\t"\
		"vmovaps	0x20(%%rcx)	,%%ymm3\n\t"\
		"vmovaps	    (%%rax)	,%%ymm0\n\t"\
		"vmovaps	0x20(%%rax)	,%%ymm1\n\t"\
		"vsubpd	%%ymm2,%%ymm4,%%ymm4\n\t"\
		"vsubpd	%%ymm3,%%ymm5,%%ymm5\n\t"\
		"vaddpd	%%ymm2,%%ymm2,%%ymm2\n\t"\
		"vaddpd	%%ymm3,%%ymm3,%%ymm3\n\t"\
		"vaddpd	%%ymm4,%%ymm2,%%ymm2\n\t"\
		"vaddpd	%%ymm5,%%ymm3,%%ymm3\n\t"\
		"vaddpd	%%ymm2,%%ymm0,%%ymm0\n\t"\
		"vaddpd	%%ymm3,%%ymm1,%%ymm1\n\t"\
		"vmovaps	%%ymm0		,    (%%rax)\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rax)\n\t"\
		"vmulpd	%%ymm6,%%ymm2,%%ymm2\n\t"\
		"vmulpd	%%ymm6,%%ymm3,%%ymm3\n\t"\
		"vmulpd	%%ymm7,%%ymm4,%%ymm4\n\t"\
		"vmulpd	%%ymm7,%%ymm5,%%ymm5\n\t"\
		"vaddpd	%%ymm0,%%ymm2,%%ymm2\n\t"\
		"vaddpd	%%ymm1,%%ymm3,%%ymm3\n\t"\
		"vmovaps	%%ymm2		,%%ymm0\n\t"\
		"vmovaps	%%ymm3		,%%ymm1\n\t"\
		"vaddpd	%%ymm5,%%ymm2,%%ymm2\n\t"\
		"vsubpd	%%ymm4,%%ymm3,%%ymm3\n\t"\
		"vsubpd	%%ymm5,%%ymm0,%%ymm0\n\t"\
		"vaddpd	%%ymm4,%%ymm1,%%ymm1\n\t"\
		"vmovaps	%%ymm2		,    (%%rbx)\n\t"\
		"vmovaps	%%ymm3		,0x20(%%rbx)\n\t"\
		"vmovaps	%%ymm0		,    (%%rcx)\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rcx)\n\t"\
	/*****************************************/\
	/* now do three more radix-3 transforms: */\
	/*****************************************/\
	/* Block 1: */\
	"movq	%[__o0]	,%%rsi 	\n\t"/* __o0-8: esi,edi store output addresses throughout */\
	"movq	%[__i0]	,%%rax\n\t"\
	"movq	%[__i3]	,%%rbx\n\t"\
	"movq	%[__i6]	,%%rcx\n\t"\
		"vmovaps	    (%%rbx)	,%%ymm4\n\t"\
		"vmovaps	0x20(%%rbx)	,%%ymm5\n\t"\
		"vmovaps	    (%%rcx)	,%%ymm2\n\t"\
		"vmovaps	0x20(%%rcx)	,%%ymm3\n\t"\
		"vmovaps	    (%%rax)	,%%ymm0\n\t"\
		"vmovaps	0x20(%%rax)	,%%ymm1\n\t"\
		"vsubpd	%%ymm2,%%ymm4,%%ymm4\n\t"\
		"vsubpd	%%ymm3,%%ymm5,%%ymm5\n\t"\
		"vaddpd	%%ymm2,%%ymm2,%%ymm2\n\t"\
		"vaddpd	%%ymm3,%%ymm3,%%ymm3\n\t"\
		"vaddpd	%%ymm4,%%ymm2,%%ymm2\n\t"\
		"vaddpd	%%ymm5,%%ymm3,%%ymm3\n\t"\
		"vaddpd	%%ymm2,%%ymm0,%%ymm0\n\t"\
		"vaddpd	%%ymm3,%%ymm1,%%ymm1\n\t"\
		"vmovaps	%%ymm0		,    (%%rsi)\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rsi)\n\t"\
		"vmulpd	%%ymm6,%%ymm2,%%ymm2\n\t"\
		"vmulpd	%%ymm6,%%ymm3,%%ymm3\n\t"\
		"vmulpd	%%ymm7,%%ymm4,%%ymm4\n\t"\
		"vmulpd	%%ymm7,%%ymm5,%%ymm5\n\t"\
		"vaddpd	%%ymm0,%%ymm2,%%ymm2\n\t"\
		"vaddpd	%%ymm1,%%ymm3,%%ymm3\n\t"\
		"vmovaps	%%ymm2		,%%ymm0\n\t"\
		"vmovaps	%%ymm3		,%%ymm1\n\t"\
	"movq	%[__o3]	,%%rdi 	\n\t"\
	"movq	%[__o6]	,%%rsi 	\n\t"\
		"vaddpd	%%ymm5,%%ymm2,%%ymm2\n\t"\
		"vsubpd	%%ymm4,%%ymm3,%%ymm3\n\t"\
		"vsubpd	%%ymm5,%%ymm0,%%ymm0\n\t"\
		"vaddpd	%%ymm4,%%ymm1,%%ymm1\n\t"\
		"vmovaps	%%ymm2		,    (%%rdi)\n\t"\
		"vmovaps	%%ymm3		,0x20(%%rdi)\n\t"\
		"vmovaps	%%ymm0		,    (%%rsi)\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rsi)\n\t"\
	/* Block 2: */\
	"movq	%[__i1]	,%%rax\n\t"\
	"movq	%[__i4]	,%%rbx\n\t"\
	"movq	%[__i7]	,%%rcx\n\t"\
		"vmovaps	    (%%rbx)	,%%ymm2\n\t"\
		"vmovaps	0x20(%%rbx)	,%%ymm3\n\t"\
		"vmovaps	    (%%rdx)	,%%ymm6\n\t"/* c1 */\
		"vmovaps	0x20(%%rdx)	,%%ymm7\n\t"\
		"vmovaps	%%ymm2		,%%ymm0\n\t"\
		"vmovaps	%%ymm3		,%%ymm1\n\t"\
		"vmulpd	%%ymm6,%%ymm2,%%ymm2\n\t"\
		"vmulpd	%%ymm6,%%ymm3,%%ymm3\n\t"\
		"vmulpd	%%ymm7,%%ymm0,%%ymm0\n\t"\
		"vmulpd	%%ymm7,%%ymm1,%%ymm1\n\t"\
		"vaddpd	%%ymm1,%%ymm2,%%ymm2\n\t"\
		"vsubpd	%%ymm0,%%ymm3,%%ymm3\n\t"\
		"vmovaps	    (%%rcx)	,%%ymm4\n\t"\
		"vmovaps	0x20(%%rcx)	,%%ymm5\n\t"\
		"vmovaps	0x40(%%rdx)	,%%ymm6\n\t"/* c2 */\
		"vmovaps	0x60(%%rdx)	,%%ymm7\n\t"\
		"vmovaps	%%ymm4		,%%ymm0\n\t"\
		"vmovaps	%%ymm5		,%%ymm1\n\t"\
		"vmulpd	%%ymm6,%%ymm4,%%ymm4\n\t"\
		"vmulpd	%%ymm6,%%ymm5,%%ymm5\n\t"\
		"vmulpd	%%ymm7,%%ymm0,%%ymm0\n\t"\
		"vmulpd	%%ymm7,%%ymm1,%%ymm1\n\t"\
		"vaddpd	%%ymm1,%%ymm4,%%ymm4\n\t"\
		"vsubpd	%%ymm0,%%ymm5,%%ymm5\n\t"\
		"vmovaps	0x80(%%rdx)	,%%ymm6\n\t"/* c3m1 */\
		"vmovaps	0xa0(%%rdx)	,%%ymm7\n\t"\
		"vmovaps	    (%%rax)	,%%ymm0\n\t"\
		"vmovaps	0x20(%%rax)	,%%ymm1\n\t"\
	"movq	%[__o7]	,%%rsi 	\n\t"\
		"vsubpd	%%ymm4,%%ymm2,%%ymm2\n\t"\
		"vsubpd	%%ymm5,%%ymm3,%%ymm3\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5\n\t"\
		"vaddpd	%%ymm2,%%ymm4,%%ymm4\n\t"\
		"vaddpd	%%ymm3,%%ymm5,%%ymm5\n\t"\
		"vaddpd	%%ymm4,%%ymm0,%%ymm0\n\t"\
		"vaddpd	%%ymm5,%%ymm1,%%ymm1\n\t"\
		"vmovaps	%%ymm0		,    (%%rsi)\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rsi)\n\t"\
		"vmulpd	%%ymm6,%%ymm4,%%ymm4\n\t"\
		"vmulpd	%%ymm6,%%ymm5,%%ymm5\n\t"\
		"vmulpd	%%ymm7,%%ymm2,%%ymm2\n\t"\
		"vmulpd	%%ymm7,%%ymm3,%%ymm3\n\t"\
		"vaddpd	%%ymm0,%%ymm4,%%ymm4\n\t"\
		"vaddpd	%%ymm1,%%ymm5,%%ymm5\n\t"\
	"movq	%[__o1]	,%%rdi 	\n\t"\
	"movq	%[__o4]	,%%rsi 	\n\t"\
		"vmovaps	%%ymm4		,%%ymm0\n\t"\
		"vmovaps	%%ymm5		,%%ymm1\n\t"\
		"vaddpd	%%ymm3,%%ymm4,%%ymm4\n\t"\
		"vsubpd	%%ymm2,%%ymm5,%%ymm5\n\t"\
		"vsubpd	%%ymm3,%%ymm0,%%ymm0\n\t"\
		"vaddpd	%%ymm2,%%ymm1,%%ymm1\n\t"\
		"vmovaps	%%ymm4		,    (%%rdi)\n\t"\
		"vmovaps	%%ymm5		,0x20(%%rdi)\n\t"\
		"vmovaps	%%ymm0		,    (%%rsi)\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rsi)\n\t"\
	/* Block 3: */\
	"movq	%[__i2]	,%%rax\n\t"\
	"movq	%[__i5]	,%%rbx\n\t"\
	"movq	%[__i8]	,%%rcx\n\t"\
		"vmovaps	    (%%rbx)	,%%ymm2\n\t"\
		"vmovaps	0x20(%%rbx)	,%%ymm3\n\t"\
		"vmovaps	0x40(%%rdx)	,%%ymm6\n\t"/* c2 */\
		"vmovaps	0x60(%%rdx)	,%%ymm7\n\t"\
		"vmovaps	%%ymm2		,%%ymm0\n\t"\
		"vmovaps	%%ymm3		,%%ymm1\n\t"\
		"vmulpd	%%ymm6,%%ymm2,%%ymm2\n\t"\
		"vmulpd	%%ymm6,%%ymm3,%%ymm3\n\t"\
		"vmulpd	%%ymm7,%%ymm0,%%ymm0\n\t"\
		"vmulpd	%%ymm7,%%ymm1,%%ymm1\n\t"\
		"vaddpd	%%ymm1,%%ymm2,%%ymm2\n\t"\
		"vsubpd	%%ymm0,%%ymm3,%%ymm3\n\t"\
		"vmovaps	    (%%rcx)	,%%ymm4\n\t"\
		"vmovaps	0x20(%%rcx)	,%%ymm5\n\t"\
		"vmovaps	0xc0(%%rdx)	,%%ymm6\n\t"/* c4 */\
		"vmovaps	0xe0(%%rdx)	,%%ymm7\n\t"\
		"vmovaps	%%ymm4		,%%ymm0\n\t"\
		"vmovaps	%%ymm5		,%%ymm1\n\t"\
		"vmulpd	%%ymm6,%%ymm4,%%ymm4\n\t"\
		"vmulpd	%%ymm6,%%ymm5,%%ymm5\n\t"\
		"vmulpd	%%ymm7,%%ymm0,%%ymm0\n\t"\
		"vmulpd	%%ymm7,%%ymm1,%%ymm1\n\t"\
		"vaddpd	%%ymm1,%%ymm4,%%ymm4\n\t"\
		"vsubpd	%%ymm0,%%ymm5,%%ymm5\n\t"\
		"vmovaps	0x80(%%rdx)	,%%ymm6\n\t"\
		"vmovaps	0xa0(%%rdx)	,%%ymm7\n\t"\
		"vmovaps	    (%%rax)	,%%ymm0\n\t"\
		"vmovaps	0x20(%%rax)	,%%ymm1\n\t"\
	"movq	%[__o5]	,%%rsi 	\n\t"\
		"vsubpd	%%ymm4,%%ymm2,%%ymm2\n\t"\
		"vsubpd	%%ymm5,%%ymm3,%%ymm3\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5\n\t"\
		"vaddpd	%%ymm2,%%ymm4,%%ymm4\n\t"\
		"vaddpd	%%ymm3,%%ymm5,%%ymm5\n\t"\
		"vaddpd	%%ymm4,%%ymm0,%%ymm0\n\t"\
		"vaddpd	%%ymm5,%%ymm1,%%ymm1\n\t"\
		"vmovaps	%%ymm0		,    (%%rsi)\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rsi)\n\t"\
		"vmulpd	%%ymm6,%%ymm4,%%ymm4\n\t"\
		"vmulpd	%%ymm6,%%ymm5,%%ymm5\n\t"\
		"vmulpd	%%ymm7,%%ymm2,%%ymm2\n\t"\
		"vmulpd	%%ymm7,%%ymm3,%%ymm3\n\t"\
		"vaddpd	%%ymm0,%%ymm4,%%ymm4\n\t"\
		"vaddpd	%%ymm1,%%ymm5,%%ymm5\n\t"\
	"movq	%[__o8]	,%%rdi 	\n\t"\
	"movq	%[__o2]	,%%rsi 	\n\t"\
		"vmovaps	%%ymm4		,%%ymm0\n\t"\
		"vmovaps	%%ymm5		,%%ymm1\n\t"\
		"vaddpd	%%ymm3,%%ymm4,%%ymm4\n\t"\
		"vsubpd	%%ymm2,%%ymm5,%%ymm5\n\t"\
		"vsubpd	%%ymm3,%%ymm0,%%ymm0\n\t"\
		"vaddpd	%%ymm2,%%ymm1,%%ymm1\n\t"\
		"vmovaps	%%ymm4		,    (%%rdi)\n\t"\
		"vmovaps	%%ymm5		,0x20(%%rdi)\n\t"\
		"vmovaps	%%ymm0		,    (%%rsi)\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rsi)\n\t"\
		:					/* outputs: none */\
		: [__i0] "m" (Xi0)	/* All inputs from memory addresses here */\
		 ,[__i1] "m" (Xi1)\
		 ,[__i2] "m" (Xi2)\
		 ,[__i3] "m" (Xi3)\
		 ,[__i4] "m" (Xi4)\
		 ,[__i5] "m" (Xi5)\
		 ,[__i6] "m" (Xi6)\
		 ,[__i7] "m" (Xi7)\
		 ,[__i8] "m" (Xi8)\
		 ,[__c1] "m" (Xc1)\
		 ,[__o0] "m" (Xo0)\
		 ,[__o1] "m" (Xo1)\
		 ,[__o2] "m" (Xo2)\
		 ,[__o3] "m" (Xo3)\
		 ,[__o4] "m" (Xo4)\
		 ,[__o5] "m" (Xo5)\
		 ,[__o6] "m" (Xo6)\
		 ,[__o7] "m" (Xo7)\
		 ,[__o8] "m" (Xo8)\
		: "cc","memory","rax","rbx","rcx","rdx","rsi","rdi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
		);\
	}

#endif	// ifdef (USE_AVX && USE_AVX2)

#if defined(USE_SSE2) && !defined(USE_AVX) && !defined(USE_ARM_V8_SIMD)	// SSE2-only, 32 and 64-bit modes:

  #if OS_BITS == 64

	// Simple 64-bit-ified version of the same-named 32-bit ASM macros, using just xmm0-7:

	/*...Radix-9 DIF: Ins in memory locations __i0-8.\
		Outs at memory locations __o0-8, assumed disjoint with inputs:\
	*/\
	#define SSE2_RADIX_09_DIF(Xi0,Xi1,Xi2,Xi3,Xi4,Xi5,Xi6,Xi7,Xi8, Xc1, Xo0,Xo1,Xo2,Xo3,Xo4,Xo5,Xo6,Xo7,Xo8)\
	{\
	__asm__ volatile (\
	"movq	%[__i0]	,%%rax 	\n\t"/* __i0-8; e[abc]x store input addresses */\
	"movq	%[__o0]	,%%rsi	\n\t"\
	"movq	%[__c1]	,%%rdx 	\n\t"/* edx stores trig addresses throughout */\
	/* Block 1: */\
		"movq	%%rax		,%%rbx\n\t"\
		"movq	%%rax		,%%rcx\n\t"\
		"movq	%[__i3]	,%%rbx 	\n\t"\
		"movq	%[__i6]	,%%rcx 	\n\t"\
		"movaps	    (%%rbx)	,%%xmm2\n\t"\
		"movaps	0x10(%%rbx)	,%%xmm3\n\t"\
		"movaps	    (%%rax)	,%%xmm0\n\t"\
		"movaps	0x10(%%rax)	,%%xmm1\n\t"\
		"movaps	    (%%rcx)	,%%xmm6\n\t"\
		"movaps	0x10(%%rcx)	,%%xmm7\n\t"\
		"movaps	%%xmm2		,%%xmm4\n\t"\
		"movaps	%%xmm3		,%%xmm5\n\t"\
		"addpd	%%xmm6		,%%xmm2\n\t"\
		"addpd	%%xmm7		,%%xmm3\n\t"\
		"subpd	%%xmm6		,%%xmm4\n\t"\
		"subpd	%%xmm7		,%%xmm5\n\t"\
		"addpd	%%xmm2		,%%xmm0\n\t"\
		"addpd	%%xmm3		,%%xmm1\n\t"\
		"movaps	0x40(%%rdx)	,%%xmm6\n\t"/* c3m1 */\
		"movaps	0x50(%%rdx)	,%%xmm7\n\t"\
		"movaps	%%xmm0		,    (%%rsi)	/* <- t00 */\n\t"\
		"movaps	%%xmm1		,0x10(%%rsi)	/* <- t01 */\n\t"\
	"movq	%[__o1]	,%%rdi	\n\t"\
	"movq	%[__o2]	,%%rsi	\n\t"\
		"mulpd	%%xmm6		,%%xmm2\n\t"\
		"mulpd	%%xmm6		,%%xmm3\n\t"\
		"mulpd	%%xmm7		,%%xmm4\n\t"\
		"mulpd	%%xmm7		,%%xmm5\n\t"\
		"addpd	%%xmm0		,%%xmm2\n\t"\
		"addpd	%%xmm1		,%%xmm3\n\t"\
		"movaps	%%xmm2		,%%xmm0\n\t"\
		"movaps	%%xmm3		,%%xmm1\n\t"\
		"subpd	%%xmm5		,%%xmm2\n\t"\
		"addpd	%%xmm4		,%%xmm3\n\t"\
		"addpd	%%xmm5		,%%xmm0\n\t"\
		"subpd	%%xmm4		,%%xmm1\n\t"\
		"movaps	%%xmm2		,    (%%rdi)	/* <- t02 */\n\t"\
		"movaps	%%xmm3		,0x10(%%rdi)	/* <- t03 */\n\t"\
		"movaps	%%xmm0		,    (%%rsi)	/* <- t04 */\n\t"\
		"movaps	%%xmm1		,0x10(%%rsi)	/* <- t05 */\n\t"\
	/* Block 2: */\
		"movq	%[__i1]	,%%rax 	\n\t"\
		"movq	%[__i4]	,%%rbx 	\n\t"\
		"movq	%[__i7]	,%%rcx 	\n\t"\
		"movaps	    (%%rbx)	,%%xmm4\n\t"\
		"movaps	0x10(%%rbx)	,%%xmm5\n\t"\
		"movaps	    (%%rax)	,%%xmm0\n\t"\
		"movaps	0x10(%%rax)	,%%xmm1\n\t"\
		"movaps	    (%%rcx)	,%%xmm2\n\t"\
		"movaps	0x10(%%rcx)	,%%xmm3\n\t"\
	"movq	%[__o3]	,%%rdi	\n\t"\
		"subpd	%%xmm2		,%%xmm4\n\t"\
		"subpd	%%xmm3		,%%xmm5\n\t"\
		"addpd	%%xmm2		,%%xmm2\n\t"\
		"addpd	%%xmm3		,%%xmm3\n\t"\
		"addpd	%%xmm4		,%%xmm2\n\t"\
		"addpd	%%xmm5		,%%xmm3\n\t"\
		"addpd	%%xmm2		,%%xmm0\n\t"\
		"addpd	%%xmm3		,%%xmm1\n\t"\
		"movaps	%%xmm0		,    (%%rdi)	/* <- t06 */\n\t"\
		"movaps	%%xmm1		,0x10(%%rdi)	/* <- t07 */\n\t"\
	"movq	%[__o4]	,%%rsi	\n\t"\
	"movq	%[__o5]	,%%rdi	\n\t"\
		"mulpd	%%xmm6		,%%xmm2\n\t"\
		"mulpd	%%xmm6		,%%xmm3\n\t"\
		"mulpd	%%xmm7		,%%xmm4\n\t"\
		"mulpd	%%xmm7		,%%xmm5\n\t"\
		"addpd	%%xmm0		,%%xmm2\n\t"\
		"addpd	%%xmm1		,%%xmm3\n\t"\
		"movaps	%%xmm2		,%%xmm0\n\t"\
		"movaps	%%xmm3		,%%xmm1\n\t"\
		"subpd	%%xmm5		,%%xmm2\n\t"\
		"addpd	%%xmm4		,%%xmm3\n\t"\
		"addpd	%%xmm5		,%%xmm0\n\t"\
		"subpd	%%xmm4		,%%xmm1\n\t"\
		"movaps	%%xmm2		,    (%%rsi)	/* <- t08 */\n\t"\
		"movaps	%%xmm3		,0x10(%%rsi)	/* <- t09 */\n\t"\
		"movaps	%%xmm0		,    (%%rdi)	/* <- t0a */\n\t"\
		"movaps	%%xmm1		,0x10(%%rdi)	/* <- t0b */\n\t"\
	/* Block 3: */\
		"movq	%[__i2]	,%%rax 	\n\t"\
		"movq	%[__i5]	,%%rbx 	\n\t"\
		"movq	%[__i8]	,%%rcx 	\n\t"\
		"movaps	    (%%rbx)	,%%xmm4\n\t"\
		"movaps	0x10(%%rbx)	,%%xmm5\n\t"\
		"movaps	    (%%rax)	,%%xmm0\n\t"\
		"movaps	0x10(%%rax)	,%%xmm1\n\t"\
		"movaps	    (%%rcx)	,%%xmm2\n\t"\
		"movaps	0x10(%%rcx)	,%%xmm3\n\t"\
	"movq	%[__o6]	,%%rsi	\n\t"\
		"subpd	%%xmm2		,%%xmm4\n\t"\
		"subpd	%%xmm3		,%%xmm5\n\t"\
		"addpd	%%xmm2		,%%xmm2\n\t"\
		"addpd	%%xmm3		,%%xmm3\n\t"\
		"addpd	%%xmm4		,%%xmm2\n\t"\
		"addpd	%%xmm5		,%%xmm3\n\t"\
		"addpd	%%xmm2		,%%xmm0\n\t"\
		"addpd	%%xmm3		,%%xmm1\n\t"\
		"movaps	%%xmm0		,    (%%rsi)	/* <- t0c */\n\t"\
		"movaps	%%xmm1		,0x10(%%rsi)	/* <- t0d */\n\t"\
	"movq	%[__o7]	,%%rdi	\n\t"\
	"movq	%[__o8]	,%%rsi	\n\t"\
		"mulpd	%%xmm6		,%%xmm2\n\t"\
		"mulpd	%%xmm6		,%%xmm3\n\t"\
		"mulpd	%%xmm7		,%%xmm4\n\t"\
		"mulpd	%%xmm7		,%%xmm5\n\t"\
		"addpd	%%xmm0		,%%xmm2\n\t"\
		"addpd	%%xmm1		,%%xmm3\n\t"\
		"movaps	%%xmm2		,%%xmm0\n\t"\
		"movaps	%%xmm3		,%%xmm1\n\t"\
		"subpd	%%xmm5		,%%xmm2\n\t"\
		"addpd	%%xmm4		,%%xmm3\n\t"\
		"addpd	%%xmm5		,%%xmm0\n\t"\
		"subpd	%%xmm4		,%%xmm1\n\t"\
		"movaps	%%xmm2		,    (%%rdi)	/* <- t0e */\n\t"\
		"movaps	%%xmm3		,0x10(%%rdi)	/* <- t0f */\n\t"\
		"movaps	%%xmm0		,    (%%rsi)	/* <- t0g */\n\t"\
		"movaps	%%xmm1		,0x10(%%rsi)	/* <- t0h */\n\t"\
	/*****************************************/\
	/* now do three more radix-3 transforms: */\
	/*****************************************/\
	/* Block 1: */\
	"movq	%[__o0]	,%%rax	\n\t"\
	"movq	%[__o3]	,%%rbx	\n\t"\
	"movq	%[__o6]	,%%rcx	\n\t"\
		"movaps	    (%%rbx)	,%%xmm4\n\t"\
		"movaps	0x10(%%rbx)	,%%xmm5\n\t"\
		"movaps	    (%%rcx)	,%%xmm2\n\t"\
		"movaps	0x10(%%rcx)	,%%xmm3\n\t"\
		"movaps	    (%%rax)	,%%xmm0\n\t"\
		"movaps	0x10(%%rax)	,%%xmm1\n\t"\
		"subpd	%%xmm2		,%%xmm4\n\t"\
		"subpd	%%xmm3		,%%xmm5\n\t"\
		"addpd	%%xmm2		,%%xmm2\n\t"\
		"addpd	%%xmm3		,%%xmm3\n\t"\
		"addpd	%%xmm4		,%%xmm2\n\t"\
		"addpd	%%xmm5		,%%xmm3\n\t"\
		"addpd	%%xmm2		,%%xmm0\n\t"\
		"addpd	%%xmm3		,%%xmm1\n\t"\
		"movaps	%%xmm0		,    (%%rax)\n\t"\
		"movaps	%%xmm1		,0x10(%%rax)\n\t"\
		"mulpd	%%xmm6		,%%xmm2\n\t"\
		"mulpd	%%xmm6		,%%xmm3\n\t"\
		"mulpd	%%xmm7		,%%xmm4\n\t"\
		"mulpd	%%xmm7		,%%xmm5\n\t"\
		"addpd	%%xmm0		,%%xmm2\n\t"\
		"addpd	%%xmm1		,%%xmm3\n\t"\
		"movaps	%%xmm2		,%%xmm0\n\t"\
		"movaps	%%xmm3		,%%xmm1\n\t"\
		"subpd	%%xmm5		,%%xmm2\n\t"\
		"addpd	%%xmm4		,%%xmm3\n\t"\
		"addpd	%%xmm5		,%%xmm0\n\t"\
		"subpd	%%xmm4		,%%xmm1\n\t"\
		"movaps	%%xmm2		,    (%%rbx)\n\t"\
		"movaps	%%xmm3		,0x10(%%rbx)\n\t"\
		"movaps	%%xmm0		,    (%%rcx)\n\t"\
		"movaps	%%xmm1		,0x10(%%rcx)\n\t"\
	/* Block 2: */\
	"movq	%[__o1]	,%%rax	\n\t"\
	"movq	%[__o4]	,%%rbx	\n\t"\
	"movq	%[__o7]	,%%rcx	\n\t"\
		"movaps	    (%%rbx)	,%%xmm2\n\t"\
		"movaps	0x10(%%rbx)	,%%xmm3\n\t"\
		"movaps	    (%%rdx)	,%%xmm6\n\t"/* c1 */\
		"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
		"movaps	%%xmm2		,%%xmm0\n\t"\
		"movaps	%%xmm3		,%%xmm1\n\t"\
		"mulpd	%%xmm6		,%%xmm2\n\t"\
		"mulpd	%%xmm6		,%%xmm3\n\t"\
		"mulpd	%%xmm7		,%%xmm0\n\t"\
		"mulpd	%%xmm7		,%%xmm1\n\t"\
		"subpd	%%xmm1		,%%xmm2\n\t"\
		"addpd	%%xmm0		,%%xmm3\n\t"\
		"movaps	    (%%rcx)	,%%xmm4\n\t"\
		"movaps	0x10(%%rcx)	,%%xmm5\n\t"\
		"movaps	0x20(%%rdx)	,%%xmm6\n\t"/* c2 */\
		"movaps	0x30(%%rdx)	,%%xmm7\n\t"\
		"movaps	%%xmm4		,%%xmm0\n\t"\
		"movaps	%%xmm5		,%%xmm1\n\t"\
		"mulpd	%%xmm6		,%%xmm4\n\t"\
		"mulpd	%%xmm6		,%%xmm5\n\t"\
		"mulpd	%%xmm7		,%%xmm0\n\t"\
		"mulpd	%%xmm7		,%%xmm1\n\t"\
		"subpd	%%xmm1		,%%xmm4\n\t"\
		"addpd	%%xmm0		,%%xmm5\n\t"\
		"movaps	0x40(%%rdx)	,%%xmm6\n\t"/* c3m1 */\
		"movaps	0x50(%%rdx)	,%%xmm7\n\t"\
		"movaps	    (%%rax)	,%%xmm0\n\t"\
		"movaps	0x10(%%rax)	,%%xmm1\n\t"\
		"subpd	%%xmm4		,%%xmm2\n\t"\
		"subpd	%%xmm5		,%%xmm3\n\t"\
		"addpd	%%xmm4		,%%xmm4\n\t"\
		"addpd	%%xmm5		,%%xmm5\n\t"\
		"addpd	%%xmm2		,%%xmm4\n\t"\
		"addpd	%%xmm3		,%%xmm5\n\t"\
		"addpd	%%xmm4		,%%xmm0\n\t"\
		"addpd	%%xmm5		,%%xmm1\n\t"\
		"movaps	%%xmm0		,    (%%rax)\n\t"\
		"movaps	%%xmm1		,0x10(%%rax)\n\t"\
		"mulpd	%%xmm6		,%%xmm4\n\t"\
		"mulpd	%%xmm6		,%%xmm5\n\t"\
		"mulpd	%%xmm7		,%%xmm2\n\t"\
		"mulpd	%%xmm7		,%%xmm3\n\t"\
		"addpd	%%xmm0		,%%xmm4\n\t"\
		"addpd	%%xmm1		,%%xmm5\n\t"\
		"movaps	%%xmm4		,%%xmm0\n\t"\
		"movaps	%%xmm5		,%%xmm1\n\t"\
		"subpd	%%xmm3		,%%xmm4\n\t"\
		"addpd	%%xmm2		,%%xmm5\n\t"\
		"addpd	%%xmm3		,%%xmm0\n\t"\
		"subpd	%%xmm2		,%%xmm1\n\t"\
		"movaps	%%xmm4		,    (%%rbx)\n\t"\
		"movaps	%%xmm5		,0x10(%%rbx)\n\t"\
		"movaps	%%xmm0		,    (%%rcx)\n\t"\
		"movaps	%%xmm1		,0x10(%%rcx)\n\t"\
	/* Block 3: */\
	"movq	%[__o2]	,%%rax	\n\t"\
	"movq	%[__o5]	,%%rbx	\n\t"\
	"movq	%[__o8]	,%%rcx	\n\t"\
		"movaps	    (%%rbx)	,%%xmm2\n\t"\
		"movaps	0x10(%%rbx)	,%%xmm3\n\t"\
		"movaps	0x20(%%rdx)	,%%xmm6\n\t"/* c2 */\
		"movaps	0x30(%%rdx)	,%%xmm7\n\t"\
		"movaps	%%xmm2		,%%xmm0\n\t"\
		"movaps	%%xmm3		,%%xmm1\n\t"\
		"mulpd	%%xmm6		,%%xmm2\n\t"\
		"mulpd	%%xmm6		,%%xmm3\n\t"\
		"mulpd	%%xmm7		,%%xmm0\n\t"\
		"mulpd	%%xmm7		,%%xmm1\n\t"\
		"subpd	%%xmm1		,%%xmm2\n\t"\
		"addpd	%%xmm0		,%%xmm3\n\t"\
		"movaps	    (%%rcx)	,%%xmm4\n\t"\
		"movaps	0x10(%%rcx)	,%%xmm5\n\t"\
		"movaps	0x60(%%rdx)	,%%xmm6\n\t"/* c4 */\
		"movaps	0x70(%%rdx)	,%%xmm7\n\t"\
		"movaps	%%xmm4		,%%xmm0\n\t"\
		"movaps	%%xmm5		,%%xmm1\n\t"\
		"mulpd	%%xmm6		,%%xmm4\n\t"\
		"mulpd	%%xmm6		,%%xmm5\n\t"\
		"mulpd	%%xmm7		,%%xmm0\n\t"\
		"mulpd	%%xmm7		,%%xmm1\n\t"\
		"subpd	%%xmm1		,%%xmm4\n\t"\
		"addpd	%%xmm0		,%%xmm5\n\t"\
		"movaps	0x40(%%rdx)	,%%xmm6\n\t"/* c3m1 */\
		"movaps	0x50(%%rdx)	,%%xmm7\n\t"\
		"movaps	    (%%rax)	,%%xmm0\n\t"\
		"movaps	0x10(%%rax)	,%%xmm1\n\t"\
		"\n\t"\
		"subpd	%%xmm4		,%%xmm2\n\t"\
		"subpd	%%xmm5		,%%xmm3\n\t"\
		"addpd	%%xmm4		,%%xmm4\n\t"\
		"addpd	%%xmm5		,%%xmm5\n\t"\
		"addpd	%%xmm2		,%%xmm4\n\t"\
		"addpd	%%xmm3		,%%xmm5\n\t"\
		"addpd	%%xmm4		,%%xmm0\n\t"\
		"addpd	%%xmm5		,%%xmm1\n\t"\
		"movaps	%%xmm0		,    (%%rax)\n\t"\
		"movaps	%%xmm1		,0x10(%%rax)\n\t"\
		"mulpd	%%xmm6		,%%xmm4\n\t"\
		"mulpd	%%xmm6		,%%xmm5\n\t"\
		"mulpd	%%xmm7		,%%xmm2\n\t"\
		"mulpd	%%xmm7		,%%xmm3\n\t"\
		"addpd	%%xmm0		,%%xmm4\n\t"\
		"addpd	%%xmm1		,%%xmm5\n\t"\
		"\n\t"\
		"movaps	%%xmm4		,%%xmm0\n\t"\
		"movaps	%%xmm5		,%%xmm1\n\t"\
		"subpd	%%xmm3		,%%xmm4\n\t"\
		"addpd	%%xmm2		,%%xmm5\n\t"\
		"addpd	%%xmm3		,%%xmm0\n\t"\
		"subpd	%%xmm2		,%%xmm1\n\t"\
		"movaps	%%xmm4		,    (%%rbx)\n\t"\
		"movaps	%%xmm5		,0x10(%%rbx)\n\t"\
		"movaps	%%xmm0		,    (%%rcx)\n\t"\
		"movaps	%%xmm1		,0x10(%%rcx)\n\t"\
		:					/* outputs: none */\
		: [__i0] "m" (Xi0)	/* All inputs from memory addresses here */\
		 ,[__i1] "m" (Xi1)\
		 ,[__i2] "m" (Xi2)\
		 ,[__i3] "m" (Xi3)\
		 ,[__i4] "m" (Xi4)\
		 ,[__i5] "m" (Xi5)\
		 ,[__i6] "m" (Xi6)\
		 ,[__i7] "m" (Xi7)\
		 ,[__i8] "m" (Xi8)\
		 ,[__c1] "m" (Xc1)\
		 ,[__o0] "m" (Xo0)\
		 ,[__o1] "m" (Xo1)\
		 ,[__o2] "m" (Xo2)\
		 ,[__o3] "m" (Xo3)\
		 ,[__o4] "m" (Xo4)\
		 ,[__o5] "m" (Xo5)\
		 ,[__o6] "m" (Xo6)\
		 ,[__o7] "m" (Xo7)\
		 ,[__o8] "m" (Xo8)\
		: "cc","memory","rax","rbx","rcx","rdx","rsi","rdi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
		);\
	}

	/*...Radix-9 DIT: Ins in memory locations __i0-8.\
		Outs at memory locations __o0-8, assumed disjoint with inputs:\
	*/\
	#define SSE2_RADIX_09_DIT(Xi0,Xi1,Xi2,Xi3,Xi4,Xi5,Xi6,Xi7,Xi8, Xc1, Xo0,Xo1,Xo2,Xo3,Xo4,Xo5,Xo6,Xo7,Xo8)\
	{\
	__asm__ volatile (\
	"movq	%[__i0]	,%%rax\n\t"\
	"movq	%[__i1]	,%%rbx\n\t"\
	"movq	%[__i2]	,%%rcx\n\t"\
	"movq	%[__c1]	,%%rdx 	\n\t"/* edx stores trig addresses throughout */\
	/* Block 1: */\
		"movaps	    (%%rbx)	,%%xmm2\n\t"\
		"movaps	0x10(%%rbx)	,%%xmm3\n\t"\
		"movaps	    (%%rax)	,%%xmm0\n\t"\
		"movaps	0x10(%%rax)	,%%xmm1\n\t"\
		"movaps	    (%%rcx)	,%%xmm6\n\t"\
		"movaps	0x10(%%rcx)	,%%xmm7\n\t"\
		"movaps	%%xmm2		,%%xmm4\n\t"\
		"movaps	%%xmm3		,%%xmm5\n\t"\
		"addpd	%%xmm6		,%%xmm2\n\t"\
		"addpd	%%xmm7		,%%xmm3\n\t"\
		"subpd	%%xmm6		,%%xmm4\n\t"\
		"subpd	%%xmm7		,%%xmm5\n\t"\
		"addpd	%%xmm2		,%%xmm0\n\t"\
		"addpd	%%xmm3		,%%xmm1\n\t"\
		"movaps	0x40(%%rdx)	,%%xmm6\n\t"/* c3m1 */\
		"movaps	0x50(%%rdx)	,%%xmm7\n\t"\
		"movaps	%%xmm0		,    (%%rax)\n\t"\
		"movaps	%%xmm1		,0x10(%%rax)\n\t"\
		"mulpd	%%xmm6		,%%xmm2\n\t"\
		"mulpd	%%xmm6		,%%xmm3\n\t"\
		"mulpd	%%xmm7		,%%xmm4\n\t"\
		"mulpd	%%xmm7		,%%xmm5\n\t"\
		"addpd	%%xmm0		,%%xmm2\n\t"\
		"addpd	%%xmm1		,%%xmm3\n\t"\
		"movaps	%%xmm2		,%%xmm0\n\t"\
		"movaps	%%xmm3		,%%xmm1\n\t"\
		"addpd	%%xmm5		,%%xmm2\n\t"\
		"subpd	%%xmm4		,%%xmm3\n\t"\
		"subpd	%%xmm5		,%%xmm0\n\t"\
		"addpd	%%xmm4		,%%xmm1\n\t"\
		"movaps	%%xmm2		,    (%%rbx)\n\t"\
		"movaps	%%xmm3		,0x10(%%rbx)\n\t"\
		"movaps	%%xmm0		,    (%%rcx)\n\t"\
		"movaps	%%xmm1		,0x10(%%rcx)\n\t"\
	/* Block 2: */\
	"movq	%[__i3]	,%%rax\n\t"\
	"movq	%[__i4]	,%%rbx\n\t"\
	"movq	%[__i5]	,%%rcx\n\t"\
		"movaps	    (%%rbx)	,%%xmm4\n\t"\
		"movaps	0x10(%%rbx)	,%%xmm5\n\t"\
		"movaps	    (%%rcx)	,%%xmm2\n\t"\
		"movaps	0x10(%%rcx)	,%%xmm3\n\t"\
		"movaps	    (%%rax)	,%%xmm0\n\t"\
		"movaps	0x10(%%rax)	,%%xmm1\n\t"\
		"subpd	%%xmm2		,%%xmm4\n\t"\
		"subpd	%%xmm3		,%%xmm5\n\t"\
		"addpd	%%xmm2		,%%xmm2\n\t"\
		"addpd	%%xmm3		,%%xmm3\n\t"\
		"addpd	%%xmm4		,%%xmm2\n\t"\
		"addpd	%%xmm5		,%%xmm3\n\t"\
		"addpd	%%xmm2		,%%xmm0\n\t"\
		"addpd	%%xmm3		,%%xmm1\n\t"\
		"movaps	%%xmm0		,    (%%rax)\n\t"\
		"movaps	%%xmm1		,0x10(%%rax)\n\t"\
		"mulpd	%%xmm6		,%%xmm2\n\t"\
		"mulpd	%%xmm6		,%%xmm3\n\t"\
		"mulpd	%%xmm7		,%%xmm4\n\t"\
		"mulpd	%%xmm7		,%%xmm5\n\t"\
		"addpd	%%xmm0		,%%xmm2\n\t"\
		"addpd	%%xmm1		,%%xmm3\n\t"\
		"movaps	%%xmm2		,%%xmm0\n\t"\
		"movaps	%%xmm3		,%%xmm1\n\t"\
		"addpd	%%xmm5		,%%xmm2\n\t"\
		"subpd	%%xmm4		,%%xmm3\n\t"\
		"subpd	%%xmm5		,%%xmm0\n\t"\
		"addpd	%%xmm4		,%%xmm1\n\t"\
		"movaps	%%xmm2		,    (%%rbx)\n\t"\
		"movaps	%%xmm3		,0x10(%%rbx)\n\t"\
		"movaps	%%xmm0		,    (%%rcx)\n\t"\
		"movaps	%%xmm1		,0x10(%%rcx)\n\t"\
	/* Block 3: */\
	"movq	%[__i6]	,%%rax\n\t"\
	"movq	%[__i7]	,%%rbx\n\t"\
	"movq	%[__i8]	,%%rcx\n\t"\
		"movaps	    (%%rbx)	,%%xmm4\n\t"\
		"movaps	0x10(%%rbx)	,%%xmm5\n\t"\
		"movaps	    (%%rcx)	,%%xmm2\n\t"\
		"movaps	0x10(%%rcx)	,%%xmm3\n\t"\
		"movaps	    (%%rax)	,%%xmm0\n\t"\
		"movaps	0x10(%%rax)	,%%xmm1\n\t"\
		"subpd	%%xmm2		,%%xmm4\n\t"\
		"subpd	%%xmm3		,%%xmm5\n\t"\
		"addpd	%%xmm2		,%%xmm2\n\t"\
		"addpd	%%xmm3		,%%xmm3\n\t"\
		"addpd	%%xmm4		,%%xmm2\n\t"\
		"addpd	%%xmm5		,%%xmm3\n\t"\
		"addpd	%%xmm2		,%%xmm0\n\t"\
		"addpd	%%xmm3		,%%xmm1\n\t"\
		"movaps	%%xmm0		,    (%%rax)\n\t"\
		"movaps	%%xmm1		,0x10(%%rax)\n\t"\
		"mulpd	%%xmm6		,%%xmm2\n\t"\
		"mulpd	%%xmm6		,%%xmm3\n\t"\
		"mulpd	%%xmm7		,%%xmm4\n\t"\
		"mulpd	%%xmm7		,%%xmm5\n\t"\
		"addpd	%%xmm0		,%%xmm2\n\t"\
		"addpd	%%xmm1		,%%xmm3\n\t"\
		"movaps	%%xmm2		,%%xmm0\n\t"\
		"movaps	%%xmm3		,%%xmm1\n\t"\
		"addpd	%%xmm5		,%%xmm2\n\t"\
		"subpd	%%xmm4		,%%xmm3\n\t"\
		"subpd	%%xmm5		,%%xmm0\n\t"\
		"addpd	%%xmm4		,%%xmm1\n\t"\
		"movaps	%%xmm2		,    (%%rbx)\n\t"\
		"movaps	%%xmm3		,0x10(%%rbx)\n\t"\
		"movaps	%%xmm0		,    (%%rcx)\n\t"\
		"movaps	%%xmm1		,0x10(%%rcx)\n\t"\
	/*****************************************/\
	/* now do three more radix-3 transforms: */\
	/*****************************************/\
	/* Block 1: */\
	"movq	%[__o0]	,%%rsi 	\n\t"/* __o0-8: esi,edi store output addresses throughout */\
	"movq	%[__i0]	,%%rax\n\t"\
	"movq	%[__i3]	,%%rbx\n\t"\
	"movq	%[__i6]	,%%rcx\n\t"\
		"movaps	    (%%rbx)	,%%xmm4\n\t"\
		"movaps	0x10(%%rbx)	,%%xmm5\n\t"\
		"movaps	    (%%rcx)	,%%xmm2\n\t"\
		"movaps	0x10(%%rcx)	,%%xmm3\n\t"\
		"movaps	    (%%rax)	,%%xmm0\n\t"\
		"movaps	0x10(%%rax)	,%%xmm1\n\t"\
		"subpd	%%xmm2		,%%xmm4\n\t"\
		"subpd	%%xmm3		,%%xmm5\n\t"\
		"addpd	%%xmm2		,%%xmm2\n\t"\
		"addpd	%%xmm3		,%%xmm3\n\t"\
		"addpd	%%xmm4		,%%xmm2\n\t"\
		"addpd	%%xmm5		,%%xmm3\n\t"\
		"addpd	%%xmm2		,%%xmm0\n\t"\
		"addpd	%%xmm3		,%%xmm1\n\t"\
		"movaps	%%xmm0		,    (%%rsi)\n\t"\
		"movaps	%%xmm1		,0x10(%%rsi)\n\t"\
		"mulpd	%%xmm6		,%%xmm2\n\t"\
		"mulpd	%%xmm6		,%%xmm3\n\t"\
		"mulpd	%%xmm7		,%%xmm4\n\t"\
		"mulpd	%%xmm7		,%%xmm5\n\t"\
		"addpd	%%xmm0		,%%xmm2\n\t"\
		"addpd	%%xmm1		,%%xmm3\n\t"\
		"movaps	%%xmm2		,%%xmm0\n\t"\
		"movaps	%%xmm3		,%%xmm1\n\t"\
	"movq	%[__o3]	,%%rdi 	\n\t"\
	"movq	%[__o6]	,%%rsi 	\n\t"\
		"addpd	%%xmm5		,%%xmm2\n\t"\
		"subpd	%%xmm4		,%%xmm3\n\t"\
		"subpd	%%xmm5		,%%xmm0\n\t"\
		"addpd	%%xmm4		,%%xmm1\n\t"\
		"movaps	%%xmm2		,    (%%rdi)\n\t"\
		"movaps	%%xmm3		,0x10(%%rdi)\n\t"\
		"movaps	%%xmm0		,    (%%rsi)\n\t"\
		"movaps	%%xmm1		,0x10(%%rsi)\n\t"\
	/* Block 2: */\
	"movq	%[__i1]	,%%rax\n\t"\
	"movq	%[__i4]	,%%rbx\n\t"\
	"movq	%[__i7]	,%%rcx\n\t"\
		"movaps	    (%%rbx)	,%%xmm2\n\t"\
		"movaps	0x10(%%rbx)	,%%xmm3\n\t"\
		"movaps	    (%%rdx)	,%%xmm6\n\t"/* c1 */\
		"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
		"movaps	%%xmm2		,%%xmm0\n\t"\
		"movaps	%%xmm3		,%%xmm1\n\t"\
		"mulpd	%%xmm6		,%%xmm2\n\t"\
		"mulpd	%%xmm6		,%%xmm3\n\t"\
		"mulpd	%%xmm7		,%%xmm0\n\t"\
		"mulpd	%%xmm7		,%%xmm1\n\t"\
		"addpd	%%xmm1		,%%xmm2\n\t"\
		"subpd	%%xmm0		,%%xmm3\n\t"\
		"movaps	    (%%rcx)	,%%xmm4\n\t"\
		"movaps	0x10(%%rcx)	,%%xmm5\n\t"\
		"movaps	0x20(%%rdx)	,%%xmm6\n\t"/* c2 */\
		"movaps	0x30(%%rdx)	,%%xmm7\n\t"\
		"movaps	%%xmm4		,%%xmm0\n\t"\
		"movaps	%%xmm5		,%%xmm1\n\t"\
		"mulpd	%%xmm6		,%%xmm4\n\t"\
		"mulpd	%%xmm6		,%%xmm5\n\t"\
		"mulpd	%%xmm7		,%%xmm0\n\t"\
		"mulpd	%%xmm7		,%%xmm1\n\t"\
		"addpd	%%xmm1		,%%xmm4\n\t"\
		"subpd	%%xmm0		,%%xmm5\n\t"\
		"movaps	0x40(%%rdx)	,%%xmm6\n\t"/* c3m1 */\
		"movaps	0x50(%%rdx)	,%%xmm7\n\t"\
		"movaps	    (%%rax)	,%%xmm0\n\t"\
		"movaps	0x10(%%rax)	,%%xmm1\n\t"\
	"movq	%[__o7]	,%%rsi 	\n\t"\
		"subpd	%%xmm4		,%%xmm2\n\t"\
		"subpd	%%xmm5		,%%xmm3\n\t"\
		"addpd	%%xmm4		,%%xmm4\n\t"\
		"addpd	%%xmm5		,%%xmm5\n\t"\
		"addpd	%%xmm2		,%%xmm4\n\t"\
		"addpd	%%xmm3		,%%xmm5\n\t"\
		"addpd	%%xmm4		,%%xmm0\n\t"\
		"addpd	%%xmm5		,%%xmm1\n\t"\
		"movaps	%%xmm0		,    (%%rsi)\n\t"\
		"movaps	%%xmm1		,0x10(%%rsi)\n\t"\
		"mulpd	%%xmm6		,%%xmm4\n\t"\
		"mulpd	%%xmm6		,%%xmm5\n\t"\
		"mulpd	%%xmm7		,%%xmm2\n\t"\
		"mulpd	%%xmm7		,%%xmm3\n\t"\
		"addpd	%%xmm0		,%%xmm4\n\t"\
		"addpd	%%xmm1		,%%xmm5\n\t"\
	"movq	%[__o1]	,%%rdi 	\n\t"\
	"movq	%[__o4]	,%%rsi 	\n\t"\
		"movaps	%%xmm4		,%%xmm0\n\t"\
		"movaps	%%xmm5		,%%xmm1\n\t"\
		"addpd	%%xmm3		,%%xmm4\n\t"\
		"subpd	%%xmm2		,%%xmm5\n\t"\
		"subpd	%%xmm3		,%%xmm0\n\t"\
		"addpd	%%xmm2		,%%xmm1\n\t"\
		"movaps	%%xmm4		,    (%%rdi)\n\t"\
		"movaps	%%xmm5		,0x10(%%rdi)\n\t"\
		"movaps	%%xmm0		,    (%%rsi)\n\t"\
		"movaps	%%xmm1		,0x10(%%rsi)\n\t"\
	/* Block 3: */\
	"movq	%[__i2]	,%%rax\n\t"\
	"movq	%[__i5]	,%%rbx\n\t"\
	"movq	%[__i8]	,%%rcx\n\t"\
		"movaps	    (%%rbx)	,%%xmm2\n\t"\
		"movaps	0x10(%%rbx)	,%%xmm3\n\t"\
		"movaps	0x20(%%rdx)	,%%xmm6\n\t"/* c2 */\
		"movaps	0x30(%%rdx)	,%%xmm7\n\t"\
		"movaps	%%xmm2		,%%xmm0\n\t"\
		"movaps	%%xmm3		,%%xmm1\n\t"\
		"mulpd	%%xmm6		,%%xmm2\n\t"\
		"mulpd	%%xmm6		,%%xmm3\n\t"\
		"mulpd	%%xmm7		,%%xmm0\n\t"\
		"mulpd	%%xmm7		,%%xmm1\n\t"\
		"addpd	%%xmm1		,%%xmm2\n\t"\
		"subpd	%%xmm0		,%%xmm3\n\t"\
		"movaps	    (%%rcx)	,%%xmm4\n\t"\
		"movaps	0x10(%%rcx)	,%%xmm5\n\t"\
		"movaps	0x60(%%rdx)	,%%xmm6\n\t"/* c4 */\
		"movaps	0x70(%%rdx)	,%%xmm7\n\t"\
		"movaps	%%xmm4		,%%xmm0\n\t"\
		"movaps	%%xmm5		,%%xmm1\n\t"\
		"mulpd	%%xmm6		,%%xmm4\n\t"\
		"mulpd	%%xmm6		,%%xmm5\n\t"\
		"mulpd	%%xmm7		,%%xmm0\n\t"\
		"mulpd	%%xmm7		,%%xmm1\n\t"\
		"addpd	%%xmm1		,%%xmm4\n\t"\
		"subpd	%%xmm0		,%%xmm5\n\t"\
		"movaps	0x40(%%rdx)	,%%xmm6\n\t"\
		"movaps	0x50(%%rdx)	,%%xmm7\n\t"\
		"movaps	    (%%rax)	,%%xmm0\n\t"\
		"movaps	0x10(%%rax)	,%%xmm1\n\t"\
	"movq	%[__o5]	,%%rsi 	\n\t"\
		"subpd	%%xmm4		,%%xmm2\n\t"\
		"subpd	%%xmm5		,%%xmm3\n\t"\
		"addpd	%%xmm4		,%%xmm4\n\t"\
		"addpd	%%xmm5		,%%xmm5\n\t"\
		"addpd	%%xmm2		,%%xmm4\n\t"\
		"addpd	%%xmm3		,%%xmm5\n\t"\
		"addpd	%%xmm4		,%%xmm0\n\t"\
		"addpd	%%xmm5		,%%xmm1\n\t"\
		"movaps	%%xmm0		,    (%%rsi)\n\t"\
		"movaps	%%xmm1		,0x10(%%rsi)\n\t"\
		"mulpd	%%xmm6		,%%xmm4\n\t"\
		"mulpd	%%xmm6		,%%xmm5\n\t"\
		"mulpd	%%xmm7		,%%xmm2\n\t"\
		"mulpd	%%xmm7		,%%xmm3\n\t"\
		"addpd	%%xmm0		,%%xmm4\n\t"\
		"addpd	%%xmm1		,%%xmm5\n\t"\
	"movq	%[__o8]	,%%rdi 	\n\t"\
	"movq	%[__o2]	,%%rsi 	\n\t"\
		"movaps	%%xmm4		,%%xmm0\n\t"\
		"movaps	%%xmm5		,%%xmm1\n\t"\
		"addpd	%%xmm3		,%%xmm4\n\t"\
		"subpd	%%xmm2		,%%xmm5\n\t"\
		"subpd	%%xmm3		,%%xmm0\n\t"\
		"addpd	%%xmm2		,%%xmm1\n\t"\
		"movaps	%%xmm4		,    (%%rdi)\n\t"\
		"movaps	%%xmm5		,0x10(%%rdi)\n\t"\
		"movaps	%%xmm0		,    (%%rsi)\n\t"\
		"movaps	%%xmm1		,0x10(%%rsi)\n\t"\
		:					/* outputs: none */\
		: [__i0] "m" (Xi0)	/* All inputs from memory addresses here */\
		 ,[__i1] "m" (Xi1)\
		 ,[__i2] "m" (Xi2)\
		 ,[__i3] "m" (Xi3)\
		 ,[__i4] "m" (Xi4)\
		 ,[__i5] "m" (Xi5)\
		 ,[__i6] "m" (Xi6)\
		 ,[__i7] "m" (Xi7)\
		 ,[__i8] "m" (Xi8)\
		 ,[__c1] "m" (Xc1)\
		 ,[__o0] "m" (Xo0)\
		 ,[__o1] "m" (Xo1)\
		 ,[__o2] "m" (Xo2)\
		 ,[__o3] "m" (Xo3)\
		 ,[__o4] "m" (Xo4)\
		 ,[__o5] "m" (Xo5)\
		 ,[__o6] "m" (Xo6)\
		 ,[__o7] "m" (Xo7)\
		 ,[__o8] "m" (Xo8)\
		: "cc","memory","rax","rbx","rcx","rdx","rsi","rdi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
		);\
	}

  #elif OS_BITS == 32

	#error 32-bit OSes no longer supported for SIMD builds!

  #endif	// SSE2-only, 32 and 64-bit modes

#endif	// x86 simd version ?

#endif	/* radix09_sse_macro_h_included */

