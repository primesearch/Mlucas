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
#ifndef radix16_wrapper_square_gcc_h_included
#define radix16_wrapper_square_gcc_h_included

// See the corresponding radix16_dyadic_square header file for notes re. the AVX-mode data layout
#ifdef USE_ARM_V8_SIMD

	// Cost: 57 LDP, 40 STP, 136 ADD, 84 MUL/FMA, 32 permute (TRN1,2). Compare the 210-line macro here with the 340-line hot mess
	// that is my old 64-bit SSE2 version (which, to be fair, lacks FMA and pairwise vector-load and has just 16 vector registers).
	#define SSE2_RADIX16_WRAPPER_DIF(Xadd0,Xadd1,Xr1,Xr9,Xr17,Xr25,Xisrt2,Xcc0,Xc1,Xc2,Xc3,Xc4,Xc5,Xc6,Xc7,Xc8,Xc9,Xc10,Xc11,Xc12,Xc13,Xc14,Xc15)\
	{\
	__asm__ volatile (\
	/*************************************************************/\
	/*                  1st set of inputs:                       */\
	/*************************************************************/\
		/*...Block 1: */						/*...Block 2: */\
		"ldr	x0,%[__add0]		\n\t"\
		"ldr	x1,%[__add1]		\n\t"\
		/* vector-data-to-be-interleaved not contiguous, but can load each of the 2 vector [Re,Im]-pairs via LDP and
		TRN1 (ARM analog of x86 UNPCKLPD) and TRN2 (ARM analog of x86 UNPCKHPD) to do the needed 2x2 transposes: */\
		"ldr	x2,%[__r1]			\n\t	ldp	q14,q15,[x0,#0x20]	\n\t"\
		"ldr	x3,%[__c4]			\n\t	ldp	q16,q17,[x1,#0x20]	\n\t"\
		"ldp	q2,q3,[x0,#0x40]	\n\t	trn1	v18.2d,v14.2d,v16.2d	\n\t"\
		"ldp	q4,q5,[x1,#0x40]	\n\t	trn1	v19.2d,v15.2d,v17.2d	\n\t"\
		"trn1	v6.2d,v2.2d,v4.2d	\n\t	trn2	v12.2d,v14.2d,v16.2d	\n\t"\
		"trn1	v7.2d,v3.2d,v5.2d	\n\t	trn2	v13.2d,v15.2d,v17.2d	\n\t"\
		"trn2	v0.2d,v2.2d,v4.2d	\n\t	stp	q12,q13,[x2,#0x180]	\n\t"\
		"trn2	v1.2d,v3.2d,v5.2d	\n\t	ldp	q14,q15,[x3,#0x40]	\n\t"/* c2 */\
		"stp	q0,q1,[x2,#0x140]	\n\t	fmul	v12.2d,v18.2d,v14.2d	\n\t"\
		"ldp	q2,q3,[x3]			\n\t	fmul	v13.2d,v19.2d,v14.2d	\n\t"\
		"fmul	v0.2d,v6.2d,v2.2d	\n\t	fmls	v12.2d,v19.2d,v15.2d	\n\t"\
		"fmul	v1.2d,v7.2d,v2.2d	\n\t	fmla	v13.2d,v18.2d,v15.2d	\n\t"\
		"fmls	v0.2d,v7.2d,v3.2d	\n\t	ldp	q14,q15,[x0,#0xa0]	\n\t"\
		"fmla	v1.2d,v6.2d,v3.2d	\n\t	ldp	q16,q17,[x1,#0xa0]	\n\t"\
		"ldr	x3,%[__c12]			\n\t	trn1	v18.2d,v14.2d,v16.2d	\n\t"\
		"ldp	q2,q3,[x0,#0xc0]	\n\t	trn1	v19.2d,v15.2d,v17.2d	\n\t"\
		"ldp	q4,q5,[x1,#0xc0]	\n\t	trn2	v16.2d,v14.2d,v16.2d	\n\t"\
		"trn1	v6.2d,v2.2d,v4.2d	\n\t	trn2	v17.2d,v15.2d,v17.2d	\n\t"\
		"trn1	v7.2d,v3.2d,v5.2d	\n\t	stp	q16,q17,[x2,#0x1a0]	\n\t"\
		"trn2	v4.2d,v2.2d,v4.2d	\n\t	ldp	q14,q15,[x3,#0x40]	\n\t"/* c10 */\
		"trn2	v5.2d,v3.2d,v5.2d	\n\t	fmul	v16.2d,v18.2d,v14.2d	\n\t"\
		"stp	q4,q5,[x2,#0x160]	\n\t	fmul	v17.2d,v19.2d,v14.2d	\n\t"\
		"ldp	q2,q3,[x3]			\n\t	fmls	v16.2d,v19.2d,v15.2d	\n\t"\
		"fmul	v4.2d,v6.2d,v2.2d	\n\t	fmla	v17.2d,v18.2d,v15.2d	\n\t"\
		"fmul	v5.2d,v7.2d,v2.2d	\n\t	fsub	v14.2d,v12.2d,v16.2d	\n\t"\
		"fmls	v4.2d,v7.2d,v3.2d	\n\t	fsub	v15.2d,v13.2d,v17.2d	\n\t"\
		"fmla	v5.2d,v6.2d,v3.2d	\n\t	fadd	v12.2d,v12.2d,v16.2d	\n\t"\
		"fsub	v2.2d,v0.2d,v4.2d	\n\t	fadd	v13.2d,v13.2d,v17.2d	\n\t"\
		"fsub	v3.2d,v1.2d,v5.2d	\n\t	ldp	q18,q19,[x0,#0xe0]	\n\t"\
		"fadd	v0.2d,v0.2d,v4.2d	\n\t	ldp	q20,q21,[x1,#0xe0]	\n\t"\
		"fadd	v1.2d,v1.2d,v5.2d	\n\t	trn1	v16.2d,v18.2d,v20.2d	\n\t"\
		"ldr	x3,%[__c8]			\n\t	trn1	v17.2d,v19.2d,v21.2d	\n\t"\
		"ldp	q6,q7,[x0,#0x80]	\n\t	trn2	v20.2d,v18.2d,v20.2d	\n\t"\
		"ldp	q8,q9,[x1,#0x80]	\n\t	trn2	v21.2d,v19.2d,v21.2d	\n\t"\
		"trn1	v4.2d,v6.2d,v8.2d	\n\t	stp	q20,q21,[x2,#0x1e0]	\n\t"\
		"trn1	v5.2d,v7.2d,v9.2d	\n\t	ldp	q18,q19,[x3,#0xc0]	\n\t"/* c14 */\
		"trn2	v8.2d,v6.2d,v8.2d	\n\t	fmul	v20.2d,v16.2d,v18.2d	\n\t"\
		"trn2	v9.2d,v7.2d,v9.2d	\n\t	fmul	v21.2d,v17.2d,v18.2d	\n\t"\
		"stp	q8,q9,[x2,#0x120]	\n\t	fmls	v20.2d,v17.2d,v19.2d	\n\t"\
		"ldp	q6,q7,[x3]			\n\t	fmla	v21.2d,v16.2d,v19.2d	\n\t"\
		"fmul	v8.2d,v4.2d,v6.2d	\n\t	ldp	q22,q23,[x0,#0x60]	\n\t"\
		"fmul	v9.2d,v5.2d,v6.2d	\n\t	ldp	q18,q19,[x1,#0x60]	\n\t"\
		"fmls	v8.2d,v5.2d,v7.2d	\n\t	trn1	v16.2d,v22.2d,v18.2d	\n\t"/* lcol: v8,9 in xmm4,5 in SSE2 version */\
		"fmla	v9.2d,v4.2d,v7.2d	\n\t	trn1	v17.2d,v23.2d,v19.2d	\n\t"\
		"ldp	q10,q11,[x0]		\n\t	trn2	v18.2d,v22.2d,v18.2d	\n\t"\
		"ldp	q4 ,q5 ,[x1]		\n\t	trn2	v19.2d,v23.2d,v19.2d	\n\t"\
		"trn1	v6.2d,v10.2d,v4.2d	\n\t	stp	q18,q19,[x2,#0x1c0]	\n\t"\
		"trn1	v7.2d,v11.2d,v5.2d	\n\t	ldp	q22,q23,[x3,#0xa0]	\n\t"/* c6 */\
		"trn2	v4.2d,v10.2d,v4.2d	\n\t	fmul	v18.2d,v16.2d,v22.2d	\n\t"\
		"trn2	v5.2d,v11.2d,v5.2d	\n\t	fmul	v19.2d,v17.2d,v22.2d	\n\t"\
		"stp	q4,q5,[x2,#0x100]	\n\t	fmls	v18.2d,v17.2d,v23.2d	\n\t"\
		"fadd	v4.2d,v6.2d,v8.2d	\n\t	fmla	v19.2d,v16.2d,v23.2d	\n\t"\
		"fadd	v5.2d,v7.2d,v9.2d	\n\t	fadd	v16.2d,v18.2d,v20.2d	\n\t"/* [16,17]|[18,19] swapped */\
		"fsub	v6.2d,v6.2d,v8.2d	\n\t	fadd	v17.2d,v19.2d,v21.2d	\n\t"/* vs respective lcol pairs */\
		"fsub	v7.2d,v7.2d,v9.2d	\n\t	fsub	v18.2d,v18.2d,v20.2d	\n\t"\
		/* Finish radix-4 butterfly: */	"	fsub	v19.2d,v19.2d,v21.2d	\n\t"\
		"fadd	v8.2d,v4.2d,v0.2d	\n\t"\
		"fadd	v9.2d,v5.2d,v1.2d	\n\t"\
		"fsub	v0.2d,v4.2d,v0.2d	\n\t	fadd	v20.2d,v12.2d,v16.2d	\n\t"\
		"fsub	v1.2d,v5.2d,v1.2d	\n\t	fadd	v21.2d,v13.2d,v17.2d	\n\t"\
		"stp	q8,q9,[x2      ]	\n\t	fsub	v12.2d,v12.2d,v16.2d	\n\t"\
		"stp	q0,q1,[x2,#0x40]	\n\t	fsub	v13.2d,v13.2d,v17.2d	\n\t"\
		"fadd	v8.2d,v6.2d,v3.2d	\n\t	stp	q20,q21,[x2,#0x80]	\n\t"\
		"fadd	v9.2d,v7.2d,v2.2d	\n\t	stp	q12,q13,[x2,#0xc0]	\n\t"\
		"fsub	v3.2d,v6.2d,v3.2d	\n\t	fadd	v20.2d,v15.2d,v18.2d	\n\t"\
		"fsub	v2.2d,v7.2d,v2.2d	\n\t	fadd	v21.2d,v14.2d,v19.2d	\n\t"\
		"stp	q3,q9,[x2,#0x20]	\n\t	fsub	v15.2d,v15.2d,v18.2d	\n\t"\
		"stp	q8,q2,[x2,#0x60]	\n\t	fsub	v14.2d,v14.2d,v19.2d	\n\t"\
		"									stp	q14,q20,[x2,#0xa0]	\n\t"\
		"									stp	q21,q15,[x2,#0xe0]	\n\t"\
	/****************************************************************************************************/\
	/* Each of the above 2 radix-4 blocks is unique, but the ensuing 2 blocks [operating on odd-indexed */\
	/* elements from the unpck*pd commands which were stored to temporaries] can use a common macro:    */\
	/****************************************************************************************************/\
		"ldr	x3,%[__cc0]	\n\t"/* Use cc0 as common base-addr for both lcol&rcol roots in this section */\
		/*...Block 3: */						/*...Block 4: */\
		"									ldp q18,q19,[x2,#0x180] \n\t"/* r25 */\
		"									ldp q14,q15,[x3,#0x1a0] \n\t"/* c3  */\
		"ldp q6,q7,[x2,#0x100] \n\t"/*r17*/"fmul	v12.2d,v18.2d,v14.2d	\n\t"\
		"ldp q2,q3,[x3,#0x120] \n\t"/*c1 */"fmul	v13.2d,v19.2d,v14.2d	\n\t"\
		"fmul	v0.2d,v6.2d,v2.2d	\n\t	fmls	v12.2d,v19.2d,v15.2d	\n\t"\
		"fmul	v1.2d,v7.2d,v2.2d	\n\t	fmla	v13.2d,v18.2d,v15.2d	\n\t"\
		"fmls	v0.2d,v7.2d,v3.2d	\n\t	ldp q18,q19,[x2,#0x1a0] \n\t"/* r27 */\
		"fmla	v1.2d,v6.2d,v3.2d	\n\t	ldp q14,q15,[x3,#0x1c0] \n\t"/* c11 */\
		"ldp q6,q7,[x2,#0x120] \n\t"/*r19*/"fmul	v16.2d,v18.2d,v14.2d	\n\t"\
		"ldp q2,q3,[x3,#0x140] \n\t"/*c9 */"fmul	v17.2d,v19.2d,v14.2d	\n\t"\
		"fmul	v4.2d,v6.2d,v2.2d	\n\t	fmls	v16.2d,v19.2d,v15.2d	\n\t"\
		"fmul	v5.2d,v7.2d,v2.2d	\n\t	fmla	v17.2d,v18.2d,v15.2d	\n\t"\
		"fmls	v4.2d,v7.2d,v3.2d	\n\t	fsub	v14.2d,v12.2d,v16.2d	\n\t"\
		"fmla	v5.2d,v6.2d,v3.2d	\n\t	fsub	v15.2d,v13.2d,v17.2d	\n\t"\
		"fsub	v2.2d,v0.2d,v4.2d	\n\t	fadd	v12.2d,v12.2d,v16.2d	\n\t"\
		"fsub	v3.2d,v1.2d,v5.2d	\n\t	fadd	v13.2d,v13.2d,v17.2d	\n\t"\
		"fadd	v0.2d,v0.2d,v4.2d	\n\t	ldp q16,q17,[x2,#0x1c0] \n\t"/* r29 */\
		"fadd	v1.2d,v1.2d,v5.2d	\n\t	ldp q18,q19,[x3,#0x1e0] \n\t"/* c7  */\
		"ldp q4,q5,[x2,#0x140] \n\t"/*r21*/"fmul	v20.2d,v16.2d,v18.2d	\n\t"\
		"ldp q6,q7,[x3,#0x160] \n\t"/*c5 */"fmul	v21.2d,v17.2d,v18.2d	\n\t"\
		"fmul	v8.2d,v4.2d,v6.2d	\n\t	fmls	v20.2d,v17.2d,v19.2d	\n\t"\
		"fmul	v9.2d,v5.2d,v6.2d	\n\t	fmla	v21.2d,v16.2d,v19.2d	\n\t"\
		"fmls	v8.2d,v5.2d,v7.2d	\n\t	ldp q16,q17,[x2,#0x1e0] \n\t"/* r31 */\
		"fmla	v9.2d,v4.2d,v7.2d	\n\t	ldp q22,q23,[x3,#0x200] \n\t"/* c15 */\
		"ldp q4 ,q5 ,[x2,#0x160]\n\t"/*r23*/"fmul	v18.2d,v16.2d,v22.2d	\n\t"\
		"ldp q10,q11,[x3,#0x180]\n\t"/*c13*/"fmul	v19.2d,v17.2d,v22.2d	\n\t"\
		"fmul	v6.2d,v4.2d,v10.2d	\n\t	fmls	v18.2d,v17.2d,v23.2d	\n\t"\
		"fmul	v7.2d,v5.2d,v10.2d	\n\t	fmla	v19.2d,v16.2d,v23.2d	\n\t"\
		"fmls	v6.2d,v5.2d,v11.2d	\n\t	fsub	v16.2d,v20.2d,v18.2d	\n\t"\
		"fmla	v7.2d,v4.2d,v11.2d	\n\t	fsub	v17.2d,v21.2d,v19.2d	\n\t"\
		"fsub	v4.2d,v8.2d,v6.2d	\n\t	fadd	v18.2d,v20.2d,v18.2d	\n\t"\
		"fsub	v5.2d,v9.2d,v7.2d	\n\t	fadd	v19.2d,v21.2d,v19.2d	\n\t"\
		"fadd	v6.2d,v8.2d,v6.2d	\n\t"\
		"fadd	v7.2d,v9.2d,v7.2d	\n\t"\
	/* Finish radix-4 fly; data in v0-7,*/"	fsub	v20.2d ,v12.2d,v18.2d	\n\t"\
	/* 12-19 corr. to SSE2 xmm0-7,8-15: */"	fsub	v21.2d ,v13.2d,v19.2d	\n\t"\
		"fsub	v8.2d ,v0.2d,v6.2d	\n\t	fadd	v12.2d,v12.2d,v18.2d	\n\t"\
		"fsub	v9.2d ,v1.2d,v7.2d	\n\t	fadd	v13.2d,v13.2d,v19.2d	\n\t"\
		"fadd	v0.2d,v0.2d,v6.2d	\n\t	stp	q20,q21,[x2,#0x1c0]		\n\t"\
		"fadd	v1.2d,v1.2d,v7.2d	\n\t	stp	q12,q13,[x2,#0x180]		\n\t"\
		"stp	q8 ,q9 ,[x2,#0x140]	\n\t	fsub	v22.2d,v14.2d,v17.2d	\n\t"\
		"stp	q0 ,q1 ,[x2,#0x100]	\n\t	fsub	v23.2d,v15.2d,v16.2d	\n\t"\
		"fsub	v10.2d,v2.2d,v5.2d	\n\t	fadd	v14.2d,v14.2d,v17.2d	\n\t"\
		"fsub	v11.2d,v3.2d,v4.2d	\n\t	fadd	v15.2d,v15.2d,v16.2d	\n\t"\
		"fadd	v2.2d,v2.2d,v5.2d	\n\t	stp	q22,q15,[x2,#0x1a0]		\n\t"\
		"fadd	v3.2d,v3.2d,v4.2d	\n\t	stp	q14,q23,[x2,#0x1e0]		\n\t"\
		"stp	q10,q3 ,[x2,#0x120]	\n\t"\
		"stp	q2 ,q11,[x2,#0x160]	\n\t"\
	/**************************************************************************************/\
	/*...and now do four more radix-4 transforms, including the internal twiddle factors: */\
	/**************************************************************************************/\
		"\n\t"\
		"ldp q29,q0,[x3,#-0x10]		\n\t"/* isrt2 into v29, will discard upper-register load data */\
		/*...Block 1: t1,9,17,25 */			/*...Block 3: t5,13,21,10: All x2-offsets incr +0x40 in rcol w.r.to lcol: */\
		"									ldp q12,q13,[x2,#0x040]		\n\t"\
		"									ldp q20,q21,[x2,#0x0c0]		\n\t"\
		"ldp q0,q1,[x2       ]		\n\t	fadd	v15.2d,v12.2d,v21.2d	\n\t"\
		"ldp q8,q9,[x2,#0x080]		\n\t	fadd	v14.2d,v13.2d,v20.2d	\n\t"\
		"fadd	v2.2d,v0.2d,v8.2d	\n\t	fsub	v12.2d,v12.2d,v21.2d	\n\t"\
		"fadd	v3.2d,v1.2d,v9.2d	\n\t	fsub	v13.2d,v13.2d,v20.2d	\n\t"\
		"fsub	v0.2d,v0.2d,v8.2d	\n\t	ldp q22,q23,[x2,#0x140]		\n\t"\
		"fsub	v1.2d,v1.2d,v9.2d	\n\t	ldp q20,q21,[x2,#0x1c0]		\n\t"\
		"ldp q4,q5,[x2,#0x100]		\n\t	fsub	v16.2d,v22.2d,v23.2d	\n\t"\
		"ldp q8,q9,[x2,#0x180]		\n\t	fadd	v17.2d,v23.2d,v22.2d	\n\t"\
		"fadd	v6.2d,v4.2d,v8.2d	\n\t	fadd	v18.2d,v20.2d,v21.2d	\n\t"\
		"fadd	v7.2d,v5.2d,v9.2d	\n\t	fsub	v19.2d,v21.2d,v20.2d	\n\t"\
		"fsub	v4.2d,v4.2d,v8.2d	\n\t	fmul	v20.2d,v16.2d,v29.2d	\n\t"\
		"fsub	v5.2d,v5.2d,v9.2d	\n\t	fmul	v21.2d,v17.2d,v29.2d	\n\t"\
		"fsub	v8.2d,v2.2d,v6.2d	\n\t	fmul	v18.2d,v18.2d,v29.2d	\n\t"\
		"fsub	v9.2d,v3.2d,v7.2d	\n\t	fmul	v19.2d,v19.2d,v29.2d	\n\t"\
		"fadd	v6.2d,v2.2d,v6.2d	\n\t	fsub	v16.2d,v20.2d,v18.2d	\n\t"\
		"fadd	v7.2d,v3.2d,v7.2d	\n\t	fsub	v17.2d,v21.2d,v19.2d	\n\t"\
		"stp q8,q9,[x2,#0x100]		\n\t	fadd	v18.2d,v20.2d,v18.2d	\n\t"\
		"stp q6,q7,[x2       ]		\n\t	fadd	v19.2d,v21.2d,v19.2d	\n\t"\
		"fsub	v8.2d,v0.2d,v5.2d	\n\t	fsub	v20.2d,v12.2d,v16.2d	\n\t"\
		"fsub	v9.2d,v1.2d,v4.2d	\n\t	fsub	v21.2d,v14.2d,v17.2d	\n\t"\
		"fadd	v5.2d,v0.2d,v5.2d	\n\t	fadd	v16.2d,v12.2d,v16.2d	\n\t"\
		"fadd	v4.2d,v1.2d,v4.2d	\n\t	fadd	v17.2d,v14.2d,v17.2d	\n\t"\
		"stp q5,q9,[x2,#0x180]		\n\t	stp q20,q21,[x2,#0x140]		\n\t"\
		"stp q8,q4,[x2,#0x080]		\n\t	stp q16,q17,[x2,#0x040]		\n\t"\
		"									fsub	v20.2d,v15.2d,v19.2d	\n\t"\
		"ldp q30,q31,[x3]	\n\t"/* c,s */"	fsub	v21.2d,v13.2d,v18.2d	\n\t"\
		"ldp q8 ,q9 ,[x2,#0x120]	\n\t	fadd	v19.2d,v15.2d,v19.2d	\n\t"/* lcol preloads first 2 pairs of Block 2 data */\
		"ldp q10,q11,[x2,#0x1a0]	\n\t	fadd	v18.2d,v13.2d,v18.2d	\n\t"\
		"									stp q19,q21,[x2,#0x1c0]		\n\t"\
		"									stp q20,q18,[x2,#0x0c0]		\n\t"\
		/*...Block 2: t3,11,19,27 */		/*...Block 4: t7,15,23,31 */\
		"fmul	v4.2d,v8.2d ,v30.2d	\n\t	ldp q20,q21,[x2,#0x160]		\n\t"\
		"fmul	v6.2d,v10.2d,v31.2d	\n\t	ldp q22,q23,[x2,#0x1e0]		\n\t"\
		"fmul	v5.2d,v9.2d ,v30.2d	\n\t	fmul	v16.2d,v20.2d,v31.2d	\n\t"\
		"fmul	v7.2d,v11.2d,v31.2d	\n\t	fmul	v18.2d,v22.2d,v30.2d	\n\t"\
		"fmls	v4.2d,v9.2d ,v31.2d	\n\t	fmul	v17.2d,v21.2d,v31.2d	\n\t"\
		"fmls	v6.2d,v11.2d,v30.2d	\n\t	fmul	v19.2d,v23.2d,v30.2d	\n\t"\
		"fmla	v5.2d,v8.2d ,v31.2d	\n\t	fmls	v16.2d,v21.2d,v30.2d	\n\t"\
		"fmla	v7.2d,v10.2d,v30.2d	\n\t	fmls	v18.2d,v23.2d,v31.2d	\n\t"\
		"fsub	v10.2d,v4.2d,v6.2d	\n\t	fmla	v17.2d,v20.2d,v30.2d	\n\t"\
		"fsub	v11.2d,v5.2d,v7.2d	\n\t	fmla	v19.2d,v22.2d,v31.2d	\n\t"\
		"fadd	v6.2d,v4.2d,v6.2d	\n\t	fsub	v22.2d,v16.2d,v18.2d	\n\t"\
		"fadd	v7.2d,v5.2d,v7.2d	\n\t	fsub	v23.2d,v17.2d,v19.2d	\n\t"\
		"ldp q8,q9,[x2,#0x0a0]		\n\t	fadd	v18.2d,v16.2d,v18.2d	\n\t"\
		"ldp q0,q1,[x2,#0x020]		\n\t	fadd	v19.2d,v17.2d,v19.2d	\n\t"\
		"fsub	v2.2d,v8.2d,v9.2d	\n\t	ldp q20,q21,[x2,#0x0e0]		\n\t"\
		"fadd	v3.2d,v9.2d,v8.2d	\n\t	ldp q12,q13,[x2,#0x060]		\n\t"\
		"fmul	v8.2d,v2.2d,v29.2d	\n\t	fadd	v14.2d,v20.2d,v21.2d	\n\t"\
		"fmul	v9.2d,v3.2d,v29.2d	\n\t	fsub	v15.2d,v21.2d,v20.2d	\n\t"\
		"fadd	v2.2d,v0.2d,v8.2d	\n\t	fmul	v20.2d,v14.2d,v29.2d	\n\t"\
		"fadd	v3.2d,v1.2d,v9.2d	\n\t	fmul	v21.2d,v15.2d,v29.2d	\n\t"\
		"fsub	v0.2d,v0.2d,v8.2d	\n\t	fadd	v14.2d,v12.2d,v20.2d	\n\t"\
		"fsub	v1.2d,v1.2d,v9.2d	\n\t	fadd	v15.2d,v13.2d,v21.2d	\n\t"\
		"fsub	v8.2d,v2.2d,v6.2d	\n\t	fsub	v12.2d,v12.2d,v20.2d	\n\t"\
		"fsub	v9.2d,v3.2d,v7.2d	\n\t	fsub	v13.2d,v13.2d,v21.2d	\n\t"\
		"fadd	v6.2d,v2.2d,v6.2d	\n\t	fsub	v16.2d,v12.2d,v22.2d	\n\t"\
		"fadd	v7.2d,v3.2d,v7.2d	\n\t	fsub	v17.2d,v13.2d,v23.2d	\n\t"\
		"stp q8,q9,[x2,#0x120]		\n\t	fadd	v12.2d,v12.2d,v22.2d	\n\t"\
		"stp q6,q7,[x2,#0x020]		\n\t	fadd	v13.2d,v13.2d,v23.2d	\n\t"\
		"fadd	v5.2d,v0.2d,v11.2d	\n\t	stp q16,q17,[x2,#0x160]		\n\t"\
		"fadd	v4.2d,v1.2d,v10.2d	\n\t	stp q12,q13,[x2,#0x060]		\n\t"\
		"fsub	v0.2d,v0.2d,v11.2d	\n\t	fadd	v20.2d,v14.2d,v19.2d	\n\t"\
		"fsub	v1.2d,v1.2d,v10.2d	\n\t	fadd	v21.2d,v15.2d,v18.2d	\n\t"\
		"stp q0,q4,[x2,#0x0a0]		\n\t	fsub	v14.2d,v14.2d,v19.2d	\n\t"\
		"stp q5,q1,[x2,#0x1a0]		\n\t	fsub	v15.2d,v15.2d,v18.2d	\n\t"\
		"									stp q14,q21,[x2,#0x0e0]		\n\t"\
		"									stp q20,q15,[x2,#0x1e0]		\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__add1] "m" (Xadd1)\
		 ,[__r1] "m" (Xr1)\
		 ,[__r9] "m" (Xr9)\
		 ,[__r17] "m" (Xr17)\
		 ,[__r25] "m" (Xr25)\
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__cc0] "m" (Xcc0)\
		 ,[__c1] "m" (Xc1)\
		 ,[__c2] "m" (Xc2)\
		 ,[__c3] "m" (Xc3)\
		 ,[__c4] "m" (Xc4)\
		 ,[__c5] "m" (Xc5)\
		 ,[__c6] "m" (Xc6)\
		 ,[__c7] "m" (Xc7)\
		 ,[__c8] "m" (Xc8)\
		 ,[__c9] "m" (Xc9)\
		 ,[__c10] "m" (Xc10)\
		 ,[__c11] "m" (Xc11)\
		 ,[__c12] "m" (Xc12)\
		 ,[__c13] "m" (Xc13)\
		 ,[__c14] "m" (Xc14)\
		 ,[__c15] "m" (Xc15)\
		: "cc","memory","x0","x1","x2","x3","v0","v1","v2","v3","v4","v5","v6","v7","v8","v9","v10","v11",\
				"v12","v13","v14","v15","v16","v17","v18","v19","v20","v21","v22","v23", "v29","v30","v31"	/* Clobbered registers */\
	);\
	}

	#define SSE2_RADIX16_WRAPPER_DIT(Xadd0,Xadd1,Xr1,Xr9,Xr17,Xr25,Xisrt2,Xcc0,Xc1,Xc2,Xc3,Xc4,Xc5,Xc6,Xc7,Xc8,Xc9,Xc10,Xc11,Xc12,Xc13,Xc14,Xc15)\
	{\
	__asm__ volatile (\
		"ldr	x2,%[__r1]			\n\t"\
		/*...Block 1: r1,9,17,25 */			/*...Block 3: r3,11,19,27 */\
		"ldr	x0,%[__add0]		\n\t	ldp q12,q13,[x2,#0x020]			\n\t"\
		"ldr	x1,%[__add1]		\n\t	ldp q20,q21,[x2,#0x120]			\n\t"\
		"ldp q0,q1,[x2       ]		\n\t	fsub	v14.2d,v12.2d,v20.2d	\n\t"\
		"ldp q8,q9,[x2,#0x100]		\n\t	fsub	v15.2d,v13.2d,v21.2d	\n\t"\
		"fsub	v2.2d,v0.2d,v8.2d	\n\t	fadd	v12.2d,v12.2d,v20.2d	\n\t"\
		"fsub	v3.2d,v1.2d,v9.2d	\n\t	fadd	v13.2d,v13.2d,v21.2d	\n\t"\
		"fadd	v0.2d,v0.2d,v8.2d	\n\t	ldp q16,q17,[x2,#0x0a0]			\n\t"\
		"fadd	v1.2d,v1.2d,v9.2d	\n\t	ldp q20,q21,[x2,#0x1a0]			\n\t"\
		"ldp q4,q5,[x2,#0x080]		\n\t	fsub	v18.2d,v16.2d,v20.2d	\n\t"\
		"ldp q8,q9,[x2,#0x180]		\n\t	fsub	v19.2d,v17.2d,v21.2d	\n\t"\
		"fsub	v6.2d,v4.2d,v8.2d	\n\t	fadd	v16.2d,v16.2d,v20.2d	\n\t"\
		"fsub	v7.2d,v5.2d,v9.2d	\n\t	fadd	v17.2d,v17.2d,v21.2d	\n\t"\
		"fadd	v4.2d,v4.2d,v8.2d	\n\t	fsub	v20.2d,v12.2d,v16.2d	\n\t"\
		"fadd	v5.2d,v5.2d,v9.2d	\n\t	fsub	v21.2d,v13.2d,v17.2d	\n\t"\
		"fsub	v8.2d,v0.2d,v4.2d	\n\t	fadd	v16.2d,v12.2d,v16.2d	\n\t"\
		"fsub	v9.2d,v1.2d,v5.2d	\n\t	fadd	v17.2d,v13.2d,v17.2d	\n\t"\
		"fadd	v4.2d,v0.2d,v4.2d	\n\t	stp q20,q21,[x2,#0x120]			\n\t"\
		"fadd	v5.2d,v1.2d,v5.2d	\n\t	stp q16,q17,[x2,#0x020]			\n\t"\
		"stp q8,q9,[x2,#0x100]		\n\t	fsub	v22.2d,v14.2d,v19.2d	\n\t"\
		"stp q4,q5,[x2       ]		\n\t	fsub	v23.2d,v15.2d,v18.2d	\n\t"\
		"fsub	v10.2d,v2.2d,v7.2d	\n\t	fadd	v19.2d ,v14.2d,v19.2d	\n\t"\
		"fsub	v11.2d,v3.2d,v6.2d	\n\t	fadd	v18.2d ,v15.2d,v18.2d	\n\t"\
		"fadd	v7.2d ,v2.2d,v7.2d	\n\t	stp q19 ,q23,[x2,#0x0a0]		\n\t"\
		"fadd	v6.2d ,v3.2d,v6.2d	\n\t	stp q22,q18 ,[x2,#0x1a0]		\n\t"\
		"stp q7 ,q11,[x2,#0x080]	\n\t"\
		"stp q10,q6 ,[x2,#0x180]	\n\t"\
		/*...Block 1: r1,9,17,25 */			/*...Block 3: r3,11,19,27 */\
		"									ldp q12,q13,[x2,#0x060]			\n\t"\
		"									ldp q20,q21,[x2,#0x160]			\n\t"\
		"ldp q0,q1,[x2,#0x040]		\n\t	fsub	v14.2d,v12.2d,v20.2d	\n\t"\
		"ldp q8,q9,[x2,#0x140]		\n\t	fsub	v15.2d,v13.2d,v21.2d	\n\t"\
		"fsub	v2.2d,v0.2d,v8.2d	\n\t	fadd	v12.2d,v12.2d,v20.2d	\n\t"\
		"fsub	v3.2d,v1.2d,v9.2d	\n\t	fadd	v13.2d,v13.2d,v21.2d	\n\t"\
		"fadd	v0.2d,v0.2d,v8.2d	\n\t	ldp q16,q17,[x2,#0x0e0]			\n\t"\
		"fadd	v1.2d,v1.2d,v9.2d	\n\t	ldp q20,q21,[x2,#0x1e0]			\n\t"\
		"ldp q4,q5,[x2,#0x0c0]		\n\t	fsub	v18.2d,v16.2d,v20.2d	\n\t"\
		"ldp q8,q9,[x2,#0x1c0]		\n\t	fsub	v19.2d,v17.2d,v21.2d	\n\t"\
		"fsub	v6.2d,v4.2d,v8.2d	\n\t	fadd	v16.2d,v16.2d,v20.2d	\n\t"\
		"fsub	v7.2d,v5.2d,v9.2d	\n\t	fadd	v17.2d,v17.2d,v21.2d	\n\t"\
		"fadd	v4.2d,v4.2d,v8.2d	\n\t	fsub	v20.2d,v12.2d,v16.2d	\n\t"\
		"fadd	v5.2d,v5.2d,v9.2d	\n\t	fsub	v21.2d,v13.2d,v17.2d	\n\t"\
		"fsub	v8.2d,v0.2d,v4.2d	\n\t	fadd	v16.2d,v12.2d,v16.2d	\n\t"\
		"fsub	v9.2d,v1.2d,v5.2d	\n\t	fadd	v17.2d,v13.2d,v17.2d	\n\t"\
		"fadd	v4.2d,v0.2d,v4.2d	\n\t	stp q20,q21,[x2,#0x160]			\n\t"\
		"fadd	v5.2d,v1.2d,v5.2d	\n\t	stp q16,q17,[x2,#0x060]			\n\t"\
		"stp q8,q9,[x2,#0x140]		\n\t	fsub	v22.2d,v14.2d,v19.2d	\n\t"\
		"stp q4,q5,[x2,#0x040]		\n\t	fsub	v23.2d,v15.2d,v18.2d	\n\t"\
		"fsub	v10.2d,v2.2d,v7.2d	\n\t	fadd	v19.2d ,v14.2d,v19.2d	\n\t"\
		"fsub	v11.2d,v3.2d,v6.2d	\n\t	fadd	v18.2d ,v15.2d,v18.2d	\n\t"\
		"fadd	v7.2d ,v2.2d,v7.2d	\n\t	stp q19 ,q23,[x2,#0x0e0]		\n\t"\
		"fadd	v6.2d ,v3.2d,v6.2d	\n\t	stp q22,q18 ,[x2,#0x1e0]		\n\t"\
		"stp q7 ,q11,[x2,#0x0c0]	\n\t"\
		"stp q10,q6 ,[x2,#0x1c0]	\n\t"\
	/***************************************************************************************************/\
	/*..and now do four more radix-4 transforms, including the internal and external twiddle factors.  */\
	/*  Write even-index 16-byte output pairs to a(j1), odd-index to a(j2), unpack same as on inputs.  */\
	/*  We do the last 2 radix-4 blocks first, to make the unpack-interleaving of outputs smoother.    */\
	/***************************************************************************************************/\
	/* Main-array addresses add0,1 in x0,1;
	all rax-offsets incr +0x100 in rcol w.r.to lcol;
	all twiddle addr-offsets +0x80 in rcol w.r.to lcol, e.g. Block uses c1,9,5,13, Block uses c3,11,7,15: */\
		"ldr	x3,%[__cc0]			\n\t	ldp q30,q31,[x3]	\n\t"/* c,s */\
		"ldp q29,q0,[x3,#-0x10]		\n\t"/* isrt2 into v29, will discard upper-register load data */\
	/*...Block 3: t3,11,19,27->r9,13,11,15	Block 4: t7,15,23,31->r25,29,27,31: */\
		"ldp q6,q7,[x2,#0x0a0]		\n\t"\
		"ldp q2,q3,[x2,#0x0e0]		\n\t"\
		"fmul	v4.2d,v6.2d,v30.2d	\n\t	ldp q18,q19,[x2,#0x1a0]			\n\t"\
		"fmul	v5.2d,v7.2d,v30.2d	\n\t	ldp q14,q15,[x2,#0x1e0]			\n\t"\
		"fmul	v0.2d,v2.2d,v31.2d	\n\t	fmul	v16.2d,v18.2d,v31.2d	\n\t"\
		"fmul	v1.2d,v3.2d,v31.2d	\n\t	fmul	v17.2d,v19.2d,v31.2d	\n\t"\
		"fmla	v4.2d,v7.2d,v31.2d	\n\t	fmul	v12.2d,v14.2d,v30.2d	\n\t"\
		"fmls	v5.2d,v6.2d,v31.2d	\n\t	fmul	v13.2d,v15.2d,v30.2d	\n\t"\
		"fmla	v0.2d,v3.2d,v30.2d	\n\t	fmla	v16.2d,v19.2d,v30.2d	\n\t"\
		"fmls	v1.2d,v2.2d,v30.2d	\n\t	fmls	v17.2d,v18.2d,v30.2d	\n\t"\
		"fsub	v6.2d,v4.2d,v0.2d	\n\t	fmla	v12.2d,v15.2d,v31.2d	\n\t"\
		"fsub	v7.2d,v5.2d,v1.2d	\n\t	fmls	v13.2d,v14.2d,v31.2d	\n\t"\
		"fadd	v4.2d,v4.2d,v0.2d	\n\t	fsub	v18.2d,v16.2d,v12.2d	\n\t"\
		"fadd	v5.2d,v5.2d,v1.2d	\n\t	fsub	v19.2d,v17.2d,v13.2d	\n\t"\
		"ldp q0,q1,[x2,#0x080]		\n\t	fadd	v16.2d,v16.2d,v12.2d	\n\t"\
		"ldp q2,q3,[x2,#0x0c0]		\n\t	fadd	v17.2d,v17.2d,v13.2d	\n\t"\
		"fadd	v8.2d,v2.2d,v3.2d	\n\t	ldp q12,q13,[x2,#0x180]			\n\t"\
		"fsub	v9.2d,v3.2d,v2.2d	\n\t	ldp q14,q15,[x2,#0x1c0]			\n\t"\
		"fmul	v8.2d,v8.2d,v29.2d	\n\t	fsub	v20.2d,v14.2d,v15.2d	\n\t"\
		"fmul	v9.2d,v9.2d,v29.2d	\n\t	fadd	v21.2d,v15.2d,v14.2d	\n\t"\
		"fadd	v2.2d,v0.2d,v8.2d	\n\t	fmul	v20.2d,v20.2d,v29.2d	\n\t"\
		"fadd	v3.2d,v1.2d,v9.2d	\n\t	fmul	v21.2d,v21.2d,v29.2d	\n\t"\
		"fsub	v0.2d,v0.2d,v8.2d	\n\t	fadd	v14.2d,v12.2d,v20.2d	\n\t"\
		"fsub	v1.2d,v1.2d,v9.2d	\n\t	fadd	v15.2d,v13.2d,v21.2d	\n\t"\
		"ldp q28,q29,[x3,#0x120]	\n\t	fsub	v12.2d,v12.2d,v20.2d	\n\t"/* lcol: c1 */\
		"fadd	v8.2d,v2.2d,v4.2d	\n\t	fsub	v13.2d,v13.2d,v21.2d	\n\t"\
		"fadd	v9.2d,v3.2d,v5.2d	\n\t	ldp q30,q31,[x3,#0x1a0]			\n\t"/* rcol: c3 */\
		"fsub	v2.2d,v2.2d,v4.2d	\n\t	fadd	v20.2d,v12.2d,v18.2d	\n\t"\
		"fsub	v3.2d,v3.2d,v5.2d	\n\t	fadd	v21.2d,v13.2d,v19.2d	\n\t"\
		"fmul	v4.2d,v8.2d,v28.2d	\n\t	fsub	v12.2d,v12.2d,v18.2d	\n\t"\
		"fmul	v5.2d,v9.2d,v28.2d	\n\t	fsub	v13.2d,v13.2d,v19.2d	\n\t"\
		"fmla	v4.2d,v9.2d,v29.2d	\n\t	fmul	v18.2d,v20.2d,v30.2d	\n\t"\
		"fmls	v5.2d,v8.2d,v29.2d	\n\t	fmul	v19.2d,v21.2d,v30.2d	\n\t"\
		"ldp q28,q29,[x3,#0x140]	\n\t	fmla	v18.2d,v21.2d,v31.2d	\n\t"/* lcol: c9 */\
		"stp	q4,q5,[x1      ]	\n\t	fmls	v19.2d,v20.2d,v31.2d	\n\t"\
/* Tmps for later unpack into add1+0,1,4,5,8,9,c,d: */	/* Tmps for later unpack into add1+2,3,6,7,a,b,e,f: */\
		"fmul	v4.2d,v2.2d,v28.2d	\n\t	ldp q30,q31,[x3,#0x1c0]			\n\t"/* rcol: c11 */\
		"fmul	v5.2d,v3.2d,v28.2d	\n\t	stp	q18,q19,[x1,#0x20]			\n\t"\
		"fmla	v4.2d,v3.2d,v29.2d	\n\t	fmul	v18.2d,v12.2d,v30.2d	\n\t"\
		"fmls	v5.2d,v2.2d,v29.2d	\n\t	fmul	v19.2d,v13.2d,v30.2d	\n\t"\
		"stp	q4,q5,[x1,#0x80]	\n\t	fmla	v18.2d,v13.2d,v31.2d	\n\t"\
		"ldp q28,q29,[x3,#0x160]	\n\t	fmls	v19.2d,v12.2d,v31.2d	\n\t"/* lcol: c5 */\
		"fadd	v8.2d,v0.2d,v7.2d	\n\t	stp	q18,q19,[x1,#0xa0]			\n\t"\
		"fadd	v9.2d,v1.2d,v6.2d	\n\t	ldp q30,q31,[x3,#0x1e0]			\n\t"/* rcol: c7 */\
		"fsub	v0.2d,v0.2d,v7.2d	\n\t	fadd	v20.2d,v14.2d,v17.2d	\n\t"\
		"fsub	v1.2d,v1.2d,v6.2d	\n\t	fadd	v21.2d,v15.2d,v16.2d	\n\t"\
		"fmul	v2.2d,v8.2d,v28.2d	\n\t	fsub	v14.2d,v14.2d,v17.2d	\n\t"\
		"fmul	v3.2d,v1.2d,v28.2d	\n\t	fsub	v15.2d,v15.2d,v16.2d	\n\t"\
		"fmla	v2.2d,v1.2d,v29.2d	\n\t	fmul	v12.2d,v20.2d,v30.2d	\n\t"\
		"fmls	v3.2d,v8.2d,v29.2d	\n\t	fmul	v13.2d,v15.2d,v30.2d	\n\t"\
		"stp	q2,q3,[x1,#0x40]	\n\t	fmla	v12.2d,v15.2d,v31.2d	\n\t"\
		"ldp q28,q29,[x3,#0x180]	\n\t	fmls	v13.2d,v20.2d,v31.2d	\n\t"/* lcol: c13 */\
		"fmul	v2.2d,v0.2d,v28.2d	\n\t	stp	q12,q13,[x1,#0x60]			\n\t"\
		"fmul	v3.2d,v9.2d,v28.2d	\n\t	ldp q30,q31,[x3,#0x200]			\n\t"/* rcol: c15 */\
		"fmla	v2.2d,v9.2d,v29.2d	\n\t	fmul	v16.2d,v14.2d,v30.2d	\n\t"\
		"fmls	v3.2d,v0.2d,v29.2d	\n\t	fmul	v17.2d,v21.2d,v30.2d	\n\t"\
		"stp	q2,q3,[x1,#0xc0]	\n\t	fmla	v16.2d,v21.2d,v31.2d	\n\t"\
	"ldp q29,q0,[x3,#-0x10]	\n\t"/* isrt2*/"fmls	v17.2d,v14.2d,v31.2d	\n\t"\
		"									stp	q16,q17,[x1,#0xe0]			\n\t"\
	/*...Block 1: t1,9,17,25 -> r1,5,3,7	Block 3: t5,13,21,29 -> r17,21,19,23: Note rcol vreg indices are +4 w.r.to SSE2 code, to reflect */\
		"									ldp q18,q19,[x2,#0x120]	\n\t"/* greater #regs available in ARM Neon and allow for full utilization of */\
		"									ldp q20,q21,[x2,#0x160]	\n\t"/* 3-operand instructions format in order to avoid reg-copying & spills */\
		"ldp q8,q9,[x2       ]		\n\t	fsub	v12.2d,v20.2d,v21.2d	\n\t"\
		"ldp q2,q3,[x2,#0x040]		\n\t	fadd	v13.2d,v21.2d,v20.2d	\n\t"\
		"fsub	v0.2d,v8.2d,v2.2d	\n\t	fadd	v16.2d,v18.2d,v19.2d	\n\t"\
		"fsub	v1.2d,v9.2d,v3.2d	\n\t	fsub	v17.2d,v19.2d,v18.2d	\n\t"\
		"fadd	v2.2d,v8.2d,v2.2d	\n\t	fmul	v18.2d,v12.2d,v29.2d	\n\t"\
		"fadd	v3.2d,v9.2d,v3.2d	\n\t	fmul	v19.2d,v13.2d,v29.2d	\n\t"\
		"ldp q10,q11,[x3,#0x40]\n\t"/* c8*/"fmul	v20.2d,v16.2d,v29.2d	\n\t"\
		"ldp q30,q31,[x3,#0xa0]\n\t"/* c2*/"fmul	v21.2d,v17.2d,v29.2d	\n\t"/* c2 is for rcol */\
		"ldp q8,q9,[x2,#0x020]		\n\t	fsub	v16.2d,v20.2d,v18.2d	\n\t"\
		"ldp q6,q7,[x2,#0x060]		\n\t	fsub	v17.2d,v21.2d,v19.2d	\n\t"\
		"fsub	v4.2d,v8.2d,v6.2d	\n\t	fadd	v18.2d,v20.2d,v18.2d	\n\t"/* v18,19 = xmm14,15, need to be preserved until later use */\
		"fsub	v5.2d,v9.2d,v7.2d	\n\t	fadd	v19.2d,v21.2d,v19.2d	\n\t"\
		"fadd	v6.2d,v8.2d,v6.2d	\n\t	ldp q20,q21,[x2,#0x100]		\n\t"\
		"fadd	v7.2d,v9.2d,v7.2d	\n\t	ldp q14,q15,[x2,#0x140]		\n\t"\
		"fadd	v8.2d,v2.2d,v6.2d	\n\t	fsub	v12.2d,v20.2d,v15.2d	\n\t"\
		"fadd	v9.2d,v3.2d,v7.2d	\n\t	fsub	v13.2d,v21.2d,v14.2d	\n\t"\
		"fsub	v2.2d,v2.2d,v6.2d	\n\t	fadd	v14.2d,v21.2d,v14.2d	\n\t"/* v12,14 = xmm8,10, need to be preserved until later use */\
		"fsub	v3.2d,v3.2d,v7.2d	\n\t	fadd	v15.2d,v20.2d,v15.2d	\n\t"\
"stp	q8,q9,[x0      ]			\n\t	fsub	v20.2d,v15.2d,v16.2d	\n\t"\
		"fmul	v8.2d,v2.2d,v10.2d	\n\t	fsub	v21.2d,v13.2d,v17.2d	\n\t"\
		"fmul	v9.2d,v3.2d,v10.2d	\n\t	fadd	v16.2d,v15.2d,v16.2d	\n\t"\
		"fmla	v8.2d,v3.2d,v11.2d	\n\t	fadd	v17.2d,v13.2d,v17.2d	\n\t"\
		"fmls	v9.2d,v2.2d,v11.2d	\n\t	stp	q20,q21,[x2,#0x100]	\n\t"/* v8,9 = xmm2,3*/\
	/*********************************************************************************************/\
	/************ ALL STORES FROM HERE ON ARE FINAL **********************************************/\
	/*********************************************************************************************/\
					/* rcol: Use v22,23 in place of 18,19 because 18,19 data must be preserved for c6 multiply sction */\
		"ldp q10,q11,[x1,#0x80]		\n\t	ldp q20,q21,[x1,#0x20]	\n\t"\
		"trn1	v2.2d,v8.2d,v10.2d	\n\t	fmul	v22.2d,v16.2d,v30.2d	\n\t"\
		"trn1	v3.2d,v9.2d,v11.2d	\n\t	fmul	v23.2d,v17.2d,v30.2d	\n\t"\
		"trn2	v6.2d,v8.2d,v10.2d	\n\t	fmla	v22.2d,v17.2d,v31.2d	\n\t"/* v22,23 = xmm16,17 */\
		"trn2	v7.2d,v9.2d,v11.2d	\n\t	fmls	v23.2d,v16.2d,v31.2d	\n\t"\
		"stp	q2,q3,[x0,#0x80]	\n\t	trn1	v16.2d,v22.2d,v20.2d	\n\t"\
		"stp	q6,q7,[x1,#0x80]	\n\t	trn1	v17.2d,v23.2d,v21.2d	\n\t"\
"ldp	q8,q9,[x0      ]			\n\t	trn2	v22.2d,v22.2d,v20.2d	\n\t"\
		"ldp q10,q11,[x1      ]		\n\t	trn2	v23.2d,v23.2d,v21.2d	\n\t"\
		"trn1	v2.2d,v8.2d,v10.2d	\n\t	stp	q16,q17,[x0,#0x20]		\n\t"\
		"trn1	v3.2d,v9.2d,v11.2d	\n\t	stp	q22,q23,[x1,#0x20]		\n\t"\
		"trn2	v6.2d,v8.2d,v10.2d	\n\t	ldp q30,q31,[x3,#0xc0]		\n\t"/* rcol: c10 */\
		"trn2	v7.2d,v9.2d,v11.2d	\n\t	ldp q22,q23,[x2,#0x100]		\n\t"\
		"stp	q2,q3,[x0      ]	\n\t	fmul	v16.2d,v22.2d,v30.2d	\n\t"\
		"stp	q6,q7,[x1      ]	\n\t	fmul	v17.2d,v23.2d,v30.2d	\n\t"\
		"fadd	v2.2d,v0.2d,v5.2d	\n\t	fmla	v16.2d,v23.2d,v31.2d	\n\t"\
		"fsub	v3.2d,v1.2d,v4.2d	\n\t	fmls	v17.2d,v22.2d,v31.2d	\n\t"\
		"ldp q28,q29,[x3,#0x60]		\n\t	ldp q20,q21,[x1,#0xa0]		\n\t"/* lcol: c4 */\
		"fmul	v8.2d,v2.2d,v28.2d	\n\t	trn1	v22.2d,v16.2d,v20.2d	\n\t"\
		"fmul	v9.2d,v3.2d,v28.2d	\n\t	trn1	v23.2d,v17.2d,v21.2d	\n\t"\
		"fmla	v8.2d,v3.2d,v29.2d	\n\t	trn2	v16.2d,v16.2d,v20.2d	\n\t"\
		"fmls	v9.2d,v2.2d,v29.2d	\n\t	trn2	v17.2d,v17.2d,v21.2d	\n\t"\
		"ldp q10,q11,[x1,#0x40]		\n\t	ldp q30,q31,[x3,#0xe0]		\n\t"/* rcol: c6 */\
		"trn1	v2.2d,v8.2d,v10.2d	\n\t	stp	q22,q23,[x0,#0xa0]		\n\t"\
		"trn1	v3.2d,v9.2d,v11.2d	\n\t	stp	q16,q17,[x1,#0xa0]		\n\t"\
		"trn2	v6.2d,v8.2d,v10.2d	\n\t	fsub	v16.2d,v12.2d,v19.2d	\n\t"/* v12,14 = xmm8,10 */\
		"trn2	v7.2d,v9.2d,v11.2d	\n\t	fsub	v17.2d,v14.2d,v18.2d	\n\t"\
		"stp	q2,q3,[x0,#0x40]	\n\t	fadd	v22.2d,v12.2d,v19.2d	\n\t"\
		"stp	q6,q7,[x1,#0x40]	\n\t	fadd	v23.2d,v14.2d,v18.2d	\n\t"\
		"ldp q20,q21,[x1,#0x60]		\n\t	fmul	v18.2d,v22.2d,v30.2d	\n\t"/* q20,21 are for rcol */\
		"ldp q28,q29,[x3,#0x80]		\n\t	fmul	v19.2d,v17.2d,v30.2d	\n\t"/* lcol: c12 */\
		"fsub	v8.2d,v0.2d,v5.2d	\n\t	fmla	v18.2d,v17.2d,v31.2d	\n\t"\
		"fadd	v9.2d,v1.2d,v4.2d	\n\t	fmls	v19.2d,v22.2d,v31.2d	\n\t"\
		"fmul	v0.2d,v8.2d,v28.2d	\n\t	trn1	v12.2d,v18.2d,v20.2d	\n\t"/* Not use v16,17 on lhs side here since v16,23 */\
		"fmul	v1.2d,v9.2d,v28.2d	\n\t	trn1	v13.2d,v19.2d,v21.2d	\n\t"/* still needed as counterparts of v17,22 below */\
		"fmla	v0.2d,v9.2d,v29.2d	\n\t	trn2	v18.2d,v18.2d,v20.2d	\n\t"\
		"fmls	v1.2d,v8.2d,v29.2d	\n\t	trn2	v19.2d,v19.2d,v21.2d	\n\t"\
		"ldp q10,q11,[x1,#0xc0]		\n\t	ldp q30,q31,[x3,#0x100]		\n\t"/* rcol: c14 */\
		"trn1	v2.2d,v0.2d,v10.2d	\n\t	stp	q12,q13,[x0,#0x60]		\n\t"\
		"trn1	v3.2d,v1.2d,v11.2d	\n\t	stp	q18,q19,[x1,#0x60]		\n\t"\
		"trn2	v6.2d,v0.2d,v10.2d	\n\t	fmul	v18.2d,v16.2d,v30.2d	\n\t"\
		"trn2	v7.2d,v1.2d,v11.2d	\n\t	fmul	v19.2d,v23.2d,v30.2d	\n\t"\
		"stp	q2,q3,[x0,#0xc0]	\n\t	fmla	v18.2d,v23.2d,v31.2d	\n\t"\
		"stp	q6,q7,[x1,#0xc0]	\n\t	fmls	v19.2d,v16.2d,v31.2d	\n\t"\
		"ldp q20,q21,[x1,#0xe0]		\n\t	trn1	v12.2d,v18.2d,v20.2d	\n\t"/* q20,21 are for rcol */\
		"									trn1	v13.2d,v19.2d,v21.2d	\n\t"\
		"									trn2	v18.2d,v18.2d,v20.2d	\n\t"\
		"									trn2	v19.2d,v19.2d,v21.2d	\n\t"\
		"									stp	q12,q13,[x0,#0xe0]		\n\t"\
		"									stp	q18,q19,[x1,#0xe0]		\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__add1] "m" (Xadd1)\
		 ,[__r1] "m" (Xr1)\
		 ,[__r9] "m" (Xr9)\
		 ,[__r17] "m" (Xr17)\
		 ,[__r25] "m" (Xr25)\
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__cc0] "m" (Xcc0)\
		 ,[__c1] "m" (Xc1)\
		 ,[__c2] "m" (Xc2)\
		 ,[__c3] "m" (Xc3)\
		 ,[__c4] "m" (Xc4)\
		 ,[__c5] "m" (Xc5)\
		 ,[__c6] "m" (Xc6)\
		 ,[__c7] "m" (Xc7)\
		 ,[__c8] "m" (Xc8)\
		 ,[__c9] "m" (Xc9)\
		 ,[__c10] "m" (Xc10)\
		 ,[__c11] "m" (Xc11)\
		 ,[__c12] "m" (Xc12)\
		 ,[__c13] "m" (Xc13)\
		 ,[__c14] "m" (Xc14)\
		 ,[__c15] "m" (Xc15)\
		: "cc","memory","x0","x1","x2","x3","v0","v1","v2","v3","v4","v5","v6","v7","v8","v9","v10","v11",\
				"v12","v13","v14","v15","v16","v17","v18","v19","v20","v21","v22","v23", "v28","v29","v30","v31"	/* Clobbered registers */\
	);\
	}

#elif defined(USE_AVX512)	// AVX512 implements a 512-bit-register version of the the AVX2 ALL_FMA-macro

	// Since the add0/1 block addresses advance in opposite directions, our prefetch scheme here is like so:
	// In fwd DIF prefetch  ahead 8 cache lines' worth of data w.r.to add0;
	// In inv DIT prefetch behind 8 cache lines' worth of data w.r.to add1;

  #ifdef USE_IMCI512	// 1st-gen Xeon Phi - Use modified 8x8 doubles-transpose algo [1a] from util.c:test_simd_transpose_8x8()

	#define SSE2_RADIX16_WRAPPER_DIF(Xadd0,Xadd1,Xadd2,Xadd3,Xadd4,Xadd5,Xadd6,Xadd7,Xr1,Xr9,Xr17,Xr25,Xisrt2,Xcc0,Xc1,Xc2,Xc3,Xc4,Xc5,Xc6,Xc7,Xc8,Xc9,Xc10,Xc11,Xc12,Xc13,Xc14,Xc15,Xpfetch_dist)\
	{\
	__asm__ volatile (\
	/* Cf. IMCI512 8x8 doubles-transpose algo [1a] in util.c ...
	do mask-register-setting first, so as to not clobber any GPRs used by compute section: */\
		"movl $0b10101010,%%ebx	\n\t movl $0b11001100,%%ecx	\n\t movl $0b11110000,%%edx	\n\t"\
		"kmov	%%ebx,%%k1		\n\t kmov	%%ecx,%%k3		\n\t kmov	%%edx,%%k5		\n\t"\
		"knot	%%k1 ,%%k2		\n\t knot	%%k3 ,%%k4		\n\t"\
		"movq	%[__add0],%%rax	\n\t"\
		"movq	%[__r1]  ,%%rsi	\n\t"\
		"leaq	0x10c0(%%rsi),%%r10		\n\t"/* &one */\
		"vmovaps	0x800(%%rsi),%%zmm29	\n\t"/* isrt2 */\
		"vmovaps	     (%%r10),%%zmm30	\n\t"/* 1.0 */\
		"vmovaps	0x040(%%r10),%%zmm31	\n\t"/* 2.0 */\
	/*************************************************************/\
	/*                  1st set of inputs:                       */\
	/*************************************************************/\
		"movq	%[__add1],%%rbx	\n\t"\
		"movq	%[__add2],%%rcx	\n\t"\
		"movq	%[__add3],%%rdx	\n\t"\
		"movq	%[__add4],%%r10	\n\t"\
		"movq	%[__add5],%%r11	\n\t"\
		"movq	%[__add6],%%r12	\n\t"\
		"movq	%[__add7],%%r13	\n\t"\
	/**** Start with 8-way interleaving - Cf. radix-32 wrapper-DFT macros for commented versions of in-register shuffle code: ****/\
	/* a[j+p0]: Inputs from add0+[0,0x100,0x200,0x300,0x400,0x500,0x600,0x700]. Outputs into local store at r1+[same byte offsets]: */\
		"vmovaps 		(%%rax),%%zmm0			\n\t	vmovaps		0x40(%%rax),%%zmm10	\n\t"\
		"vmovaps 		(%%rbx),%%zmm1			\n\t	vmovaps		0x40(%%rbx),%%zmm11	\n\t"\
		"vmovaps 		(%%rcx),%%zmm2			\n\t	vmovaps		0x40(%%rcx),%%zmm12	\n\t"\
		"vmovaps 		(%%rdx),%%zmm3			\n\t	vmovaps		0x40(%%rdx),%%zmm13	\n\t"\
		"vmovaps 		(%%r10),%%zmm4			\n\t	vmovaps		0x40(%%r10),%%zmm14	\n\t"\
		"vmovaps 		(%%r11),%%zmm5			\n\t	vmovaps		0x40(%%r11),%%zmm15	\n\t"\
		"vmovaps 		(%%r12),%%zmm6			\n\t	vmovaps		0x40(%%r12),%%zmm16	\n\t"\
		"vmovaps 		(%%r13),%%zmm7			\n\t	vmovaps		0x40(%%r13),%%zmm17	\n\t"\
		/* [1] First step is a quartet of [UNPCKLPD,UNPCKHPD] pairs to effect transposed 2x2 submatrices - */\
		/* under IMCI-512 use VBLENDMPD with 4-double swizzle of src1 to emulate AVX-512 [UNPCKLPD,UNPCKHPD]: */\
		"vblendmpd		%%zmm1%{cdab%},%%zmm0,%%zmm8%{%%k1%}	\n\t	vblendmpd		%%zmm11%{cdab%},%%zmm10,%%zmm9 %{%%k1%}	\n\t"\
		"vblendmpd		%%zmm0%{cdab%},%%zmm1,%%zmm1%{%%k2%}	\n\t	vblendmpd		%%zmm10%{cdab%},%%zmm11,%%zmm11%{%%k2%}	\n\t"\
		"vblendmpd		%%zmm3%{cdab%},%%zmm2,%%zmm0%{%%k1%}	\n\t	vblendmpd		%%zmm13%{cdab%},%%zmm12,%%zmm10%{%%k1%}	\n\t"\
		"vblendmpd		%%zmm2%{cdab%},%%zmm3,%%zmm3%{%%k2%}	\n\t	vblendmpd		%%zmm12%{cdab%},%%zmm13,%%zmm13%{%%k2%}	\n\t"\
		"vblendmpd		%%zmm5%{cdab%},%%zmm4,%%zmm2%{%%k1%}	\n\t	vblendmpd		%%zmm15%{cdab%},%%zmm14,%%zmm12%{%%k1%}	\n\t"\
		"vblendmpd		%%zmm4%{cdab%},%%zmm5,%%zmm5%{%%k2%}	\n\t	vblendmpd		%%zmm14%{cdab%},%%zmm15,%%zmm15%{%%k2%}	\n\t"\
		"vblendmpd		%%zmm7%{cdab%},%%zmm6,%%zmm4%{%%k1%}	\n\t	vblendmpd		%%zmm17%{cdab%},%%zmm16,%%zmm14%{%%k1%}	\n\t"\
		"vblendmpd		%%zmm6%{cdab%},%%zmm7,%%zmm7%{%%k2%}	\n\t	vblendmpd		%%zmm16%{cdab%},%%zmm17,%%zmm17%{%%k2%}	\n\t"\
		/* [2] Second layer of VBLENDMPD with pairwise-double swizzle of src1 gives us fully transposed double-quartets: */\
		"vblendmpd		%%zmm0%{badc%},%%zmm8,%%zmm6%{%%k3%}	\n\t	vblendmpd		%%zmm10%{badc%},%%zmm9 ,%%zmm16%{%%k3%}	\n\t"\
		"vblendmpd		%%zmm8%{badc%},%%zmm0,%%zmm0%{%%k4%}	\n\t	vblendmpd		%%zmm9 %{badc%},%%zmm10,%%zmm10%{%%k4%}	\n\t"\
		"vblendmpd		%%zmm3%{badc%},%%zmm1,%%zmm8%{%%k3%}	\n\t	vblendmpd		%%zmm13%{badc%},%%zmm11,%%zmm9 %{%%k3%}	\n\t"\
		"vblendmpd		%%zmm1%{badc%},%%zmm3,%%zmm3%{%%k4%}	\n\t	vblendmpd		%%zmm11%{badc%},%%zmm13,%%zmm13%{%%k4%}	\n\t"\
		"vblendmpd		%%zmm4%{badc%},%%zmm2,%%zmm1%{%%k3%}	\n\t	vblendmpd		%%zmm14%{badc%},%%zmm12,%%zmm11%{%%k3%}	\n\t"\
		"vblendmpd		%%zmm2%{badc%},%%zmm4,%%zmm4%{%%k4%}	\n\t	vblendmpd		%%zmm12%{badc%},%%zmm14,%%zmm14%{%%k4%}	\n\t"\
		"vblendmpd		%%zmm7%{badc%},%%zmm5,%%zmm2%{%%k3%}	\n\t	vblendmpd		%%zmm17%{badc%},%%zmm15,%%zmm12%{%%k3%}	\n\t"\
		"vblendmpd		%%zmm5%{badc%},%%zmm7,%%zmm7%{%%k4%}	\n\t	vblendmpd		%%zmm15%{badc%},%%zmm17,%%zmm17%{%%k4%}	\n\t"\
		/* [3] Swap/combine 256-bit register-halves across the register midline: */\
		"vpermf32x4	$78,%%zmm1,%%zmm1							\n\t	vpermf32x4	$78,%%zmm11,%%zmm11							\n\t"\
		"vpermf32x4	$78,%%zmm2,%%zmm2							\n\t	vpermf32x4	$78,%%zmm12,%%zmm12							\n\t"\
		"vpermf32x4	$78,%%zmm4,%%zmm4							\n\t	vpermf32x4	$78,%%zmm14,%%zmm14							\n\t"\
		"vpermf32x4	$78,%%zmm7,%%zmm7							\n\t	vpermf32x4	$78,%%zmm17,%%zmm17							\n\t"\
		"vblendmpd		%%zmm1,%%zmm6,%%zmm5%{%%k5%}			\n\t	vblendmpd		%%zmm11,%%zmm16,%%zmm15%{%%k5%}			\n\t"\
		"vblendmpd		%%zmm6,%%zmm1,%%zmm1%{%%k5%}			\n\t	vblendmpd		%%zmm16,%%zmm11,%%zmm11%{%%k5%}			\n\t"\
		"vblendmpd		%%zmm2,%%zmm8,%%zmm6%{%%k5%}			\n\t	vblendmpd		%%zmm12,%%zmm9 ,%%zmm16%{%%k5%}			\n\t"\
		"vblendmpd		%%zmm8,%%zmm2,%%zmm2%{%%k5%}			\n\t	vblendmpd		%%zmm9 ,%%zmm12,%%zmm12%{%%k5%}			\n\t"\
		"vblendmpd		%%zmm4,%%zmm0,%%zmm8%{%%k5%}			\n\t	vblendmpd		%%zmm14,%%zmm10,%%zmm9 %{%%k5%}			\n\t"\
		"vblendmpd		%%zmm0,%%zmm4,%%zmm4%{%%k5%}			\n\t	vblendmpd		%%zmm10,%%zmm14,%%zmm14%{%%k5%}			\n\t"\
		"vblendmpd		%%zmm7,%%zmm3,%%zmm0%{%%k5%}			\n\t	vblendmpd		%%zmm17,%%zmm13,%%zmm10%{%%k5%}			\n\t"\
		"vblendmpd		%%zmm3,%%zmm7,%%zmm7%{%%k5%}			\n\t	vblendmpd		%%zmm13,%%zmm17,%%zmm17%{%%k5%}			\n\t"\
		"vpermf32x4	$78,%%zmm1,%%zmm1							\n\t	vpermf32x4	$78,%%zmm11,%%zmm11							\n\t"\
		"vpermf32x4	$78,%%zmm2,%%zmm2							\n\t	vpermf32x4	$78,%%zmm12,%%zmm12							\n\t"\
		"vpermf32x4	$78,%%zmm4,%%zmm4							\n\t	vpermf32x4	$78,%%zmm14,%%zmm14							\n\t"\
		"vpermf32x4	$78,%%zmm7,%%zmm7							\n\t	vpermf32x4	$78,%%zmm17,%%zmm17							\n\t"\
	/************ Same register-index output pattern as in AVX-512 version: ************/\
		/**** rows 0,1,4,5,8,9,c,d contain original cols 0,4,2,6,1,5,3,7: ****/\
		"vmovaps		%%zmm5,     (%%rsi)		\n\t	vmovaps	%%zmm15,0x040(%%rsi)	\n\t"\
		"vmovaps		%%zmm6,0x400(%%rsi)		\n\t	vmovaps	%%zmm16,0x440(%%rsi)	\n\t"\
		"vmovaps		%%zmm8,0x200(%%rsi)		\n\t	vmovaps	%%zmm9 ,0x240(%%rsi)	\n\t"\
		"vmovaps		%%zmm0,0x600(%%rsi)		\n\t	vmovaps	%%zmm10,0x640(%%rsi)	\n\t"\
		"vmovaps		%%zmm1,0x080(%%rsi)		\n\t	vmovaps	%%zmm11,0x0c0(%%rsi)	\n\t"\
		"vmovaps		%%zmm2,0x480(%%rsi)		\n\t	vmovaps	%%zmm12,0x4c0(%%rsi)	\n\t"\
		"vmovaps		%%zmm4,0x280(%%rsi)		\n\t	vmovaps	%%zmm14,0x2c0(%%rsi)	\n\t"\
		"vmovaps		%%zmm7,0x680(%%rsi)		\n\t	vmovaps	%%zmm17,0x6c0(%%rsi)	\n\t"\
		"\n\t"\
	/* a[j+p4]: Inputs from add0+[0,0x100,0x200,0x300,0x400,0x500,0x600,0x700]+0x80. Outputs into r1+[same byte offsets]+0x40: */\
	"addq	$0x80,%%rax	\n\t"\
	"addq	$0x80,%%rbx	\n\t"\
	"addq	$0x80,%%rcx	\n\t"\
	"addq	$0x80,%%rdx	\n\t"\
	"addq	$0x80,%%r10	\n\t"\
	"addq	$0x80,%%r11	\n\t"\
	"addq	$0x80,%%r12	\n\t"\
	"addq	$0x80,%%r13	\n\t"\
	"addq	$0x100,%%rsi	\n\t"\
		"vmovaps 		(%%rax),%%zmm0			\n\t	vmovaps		0x40(%%rax),%%zmm10	\n\t"\
		"vmovaps 		(%%rbx),%%zmm1			\n\t	vmovaps		0x40(%%rbx),%%zmm11	\n\t"\
		"vmovaps 		(%%rcx),%%zmm2			\n\t	vmovaps		0x40(%%rcx),%%zmm12	\n\t"\
		"vmovaps 		(%%rdx),%%zmm3			\n\t	vmovaps		0x40(%%rdx),%%zmm13	\n\t"\
		"vmovaps 		(%%r10),%%zmm4			\n\t	vmovaps		0x40(%%r10),%%zmm14	\n\t"\
		"vmovaps 		(%%r11),%%zmm5			\n\t	vmovaps		0x40(%%r11),%%zmm15	\n\t"\
		"vmovaps 		(%%r12),%%zmm6			\n\t	vmovaps		0x40(%%r12),%%zmm16	\n\t"\
		"vmovaps 		(%%r13),%%zmm7			\n\t	vmovaps		0x40(%%r13),%%zmm17	\n\t"\
		"\n\t"\
		"vblendmpd		%%zmm1%{cdab%},%%zmm0,%%zmm8%{%%k1%}	\n\t	vblendmpd		%%zmm11%{cdab%},%%zmm10,%%zmm9 %{%%k1%}	\n\t"\
		"vblendmpd		%%zmm0%{cdab%},%%zmm1,%%zmm1%{%%k2%}	\n\t	vblendmpd		%%zmm10%{cdab%},%%zmm11,%%zmm11%{%%k2%}	\n\t"\
		"vblendmpd		%%zmm3%{cdab%},%%zmm2,%%zmm0%{%%k1%}	\n\t	vblendmpd		%%zmm13%{cdab%},%%zmm12,%%zmm10%{%%k1%}	\n\t"\
		"vblendmpd		%%zmm2%{cdab%},%%zmm3,%%zmm3%{%%k2%}	\n\t	vblendmpd		%%zmm12%{cdab%},%%zmm13,%%zmm13%{%%k2%}	\n\t"\
		"vblendmpd		%%zmm5%{cdab%},%%zmm4,%%zmm2%{%%k1%}	\n\t	vblendmpd		%%zmm15%{cdab%},%%zmm14,%%zmm12%{%%k1%}	\n\t"\
		"vblendmpd		%%zmm4%{cdab%},%%zmm5,%%zmm5%{%%k2%}	\n\t	vblendmpd		%%zmm14%{cdab%},%%zmm15,%%zmm15%{%%k2%}	\n\t"\
		"vblendmpd		%%zmm7%{cdab%},%%zmm6,%%zmm4%{%%k1%}	\n\t	vblendmpd		%%zmm17%{cdab%},%%zmm16,%%zmm14%{%%k1%}	\n\t"\
		"vblendmpd		%%zmm6%{cdab%},%%zmm7,%%zmm7%{%%k2%}	\n\t	vblendmpd		%%zmm16%{cdab%},%%zmm17,%%zmm17%{%%k2%}	\n\t"\
		"\n\t"\
		"vblendmpd		%%zmm0%{badc%},%%zmm8,%%zmm6%{%%k3%}	\n\t	vblendmpd		%%zmm10%{badc%},%%zmm9 ,%%zmm16%{%%k3%}	\n\t"\
		"vblendmpd		%%zmm8%{badc%},%%zmm0,%%zmm0%{%%k4%}	\n\t	vblendmpd		%%zmm9 %{badc%},%%zmm10,%%zmm10%{%%k4%}	\n\t"\
		"vblendmpd		%%zmm3%{badc%},%%zmm1,%%zmm8%{%%k3%}	\n\t	vblendmpd		%%zmm13%{badc%},%%zmm11,%%zmm9 %{%%k3%}	\n\t"\
		"vblendmpd		%%zmm1%{badc%},%%zmm3,%%zmm3%{%%k4%}	\n\t	vblendmpd		%%zmm11%{badc%},%%zmm13,%%zmm13%{%%k4%}	\n\t"\
		"vblendmpd		%%zmm4%{badc%},%%zmm2,%%zmm1%{%%k3%}	\n\t	vblendmpd		%%zmm14%{badc%},%%zmm12,%%zmm11%{%%k3%}	\n\t"\
		"vblendmpd		%%zmm2%{badc%},%%zmm4,%%zmm4%{%%k4%}	\n\t	vblendmpd		%%zmm12%{badc%},%%zmm14,%%zmm14%{%%k4%}	\n\t"\
		"vblendmpd		%%zmm7%{badc%},%%zmm5,%%zmm2%{%%k3%}	\n\t	vblendmpd		%%zmm17%{badc%},%%zmm15,%%zmm12%{%%k3%}	\n\t"\
		"vblendmpd		%%zmm5%{badc%},%%zmm7,%%zmm7%{%%k4%}	\n\t	vblendmpd		%%zmm15%{badc%},%%zmm17,%%zmm17%{%%k4%}	\n\t"\
		"\n\t"\
		"vpermf32x4	$78,%%zmm1,%%zmm1							\n\t	vpermf32x4	$78,%%zmm11,%%zmm11							\n\t"\
		"vpermf32x4	$78,%%zmm2,%%zmm2							\n\t	vpermf32x4	$78,%%zmm12,%%zmm12							\n\t"\
		"vpermf32x4	$78,%%zmm4,%%zmm4							\n\t	vpermf32x4	$78,%%zmm14,%%zmm14							\n\t"\
		"vpermf32x4	$78,%%zmm7,%%zmm7							\n\t	vpermf32x4	$78,%%zmm17,%%zmm17							\n\t"\
		"vblendmpd		%%zmm1,%%zmm6,%%zmm5%{%%k5%}			\n\t	vblendmpd		%%zmm11,%%zmm16,%%zmm15%{%%k5%}			\n\t"\
		"vblendmpd		%%zmm6,%%zmm1,%%zmm1%{%%k5%}			\n\t	vblendmpd		%%zmm16,%%zmm11,%%zmm11%{%%k5%}			\n\t"\
		"vblendmpd		%%zmm2,%%zmm8,%%zmm6%{%%k5%}			\n\t	vblendmpd		%%zmm12,%%zmm9 ,%%zmm16%{%%k5%}			\n\t"\
		"vblendmpd		%%zmm8,%%zmm2,%%zmm2%{%%k5%}			\n\t	vblendmpd		%%zmm9 ,%%zmm12,%%zmm12%{%%k5%}			\n\t"\
		"vblendmpd		%%zmm4,%%zmm0,%%zmm8%{%%k5%}			\n\t	vblendmpd		%%zmm14,%%zmm10,%%zmm9 %{%%k5%}			\n\t"\
		"vblendmpd		%%zmm0,%%zmm4,%%zmm4%{%%k5%}			\n\t	vblendmpd		%%zmm10,%%zmm14,%%zmm14%{%%k5%}			\n\t"\
		"vblendmpd		%%zmm7,%%zmm3,%%zmm0%{%%k5%}			\n\t	vblendmpd		%%zmm17,%%zmm13,%%zmm10%{%%k5%}			\n\t"\
		"vblendmpd		%%zmm3,%%zmm7,%%zmm7%{%%k5%}			\n\t	vblendmpd		%%zmm13,%%zmm17,%%zmm17%{%%k5%}			\n\t"\
		"vpermf32x4	$78,%%zmm1,%%zmm1							\n\t	vpermf32x4	$78,%%zmm11,%%zmm11							\n\t"\
		"vpermf32x4	$78,%%zmm2,%%zmm2							\n\t	vpermf32x4	$78,%%zmm12,%%zmm12							\n\t"\
		"vpermf32x4	$78,%%zmm4,%%zmm4							\n\t	vpermf32x4	$78,%%zmm14,%%zmm14							\n\t"\
		"vpermf32x4	$78,%%zmm7,%%zmm7							\n\t	vpermf32x4	$78,%%zmm17,%%zmm17							\n\t"\
		/**** rows 2,3,6,7,a,b,e,f contain original cols 8,c,a,e,9,d,b,f: ****/\
		"vmovaps		%%zmm5,     (%%rsi)		\n\t	vmovaps	%%zmm15,0x040(%%rsi)	\n\t"\
		"vmovaps		%%zmm6,0x400(%%rsi)		\n\t	vmovaps	%%zmm16,0x440(%%rsi)	\n\t"\
		"vmovaps		%%zmm8,0x200(%%rsi)		\n\t	vmovaps	%%zmm9 ,0x240(%%rsi)	\n\t"\
		"vmovaps		%%zmm0,0x600(%%rsi)		\n\t	vmovaps	%%zmm10,0x640(%%rsi)	\n\t"\
		"vmovaps		%%zmm1,0x080(%%rsi)		\n\t	vmovaps	%%zmm11,0x0c0(%%rsi)	\n\t"\
		"vmovaps		%%zmm2,0x480(%%rsi)		\n\t	vmovaps	%%zmm12,0x4c0(%%rsi)	\n\t"\
		"vmovaps		%%zmm4,0x280(%%rsi)		\n\t	vmovaps	%%zmm14,0x2c0(%%rsi)	\n\t"\
		"vmovaps		%%zmm7,0x680(%%rsi)		\n\t	vmovaps	%%zmm17,0x6c0(%%rsi)	\n\t"\
		"\n\t"\
	/* The data-patterning here is somewhat funky:
	Overall result rows 0-f contain original cols 0,4,8,c,2,6,a,e,1,5,9,d,3,7,b,f
	Bit-rev'd: row(br[0-f]) contain original cols 0,1,2,3,8,9,a,b,4,5,6,7,c,d,e,f, i.e. middle 2 quartets swapped.
	Ordering-by-cols we have that original main-array cols 0-f end up in local-mem rows 0,8,4,c,1,9,5,9,2,a,6,e,3,b,7,f.
	*/\
	/*****************/\
	/* Radix-16 DIF: */\
	/*****************/\
	"movq	%[__add0],%%rax	\n\t"/* Use for FMA-related spills */\
	"movq	%[__r1] ,%%rcx	\n\t"\
"leaq	0x1100(%%rcx),%%rsi	\n\t"/* two */\
	/*...Block 1: */													/*...Block 2: */\
	"leaq	0xac0(%%rcx),%%rdi	/* c2 */\n\t"\
	"leaq	0x9c0(%%rcx),%%rdx	/* c4 */\n\t"\
		"vmovaps	0x080(%%rcx),%%zmm0	/* zmm0 <-     a[jt+p4] */		\n\t		vmovaps		0x200(%%rcx),%%zmm8	/* zmm10 <-     a[jt+p2] */			\n\t"\
		"vmovaps	0x0c0(%%rcx),%%zmm1	/* zmm1 <-     a[jp+p4] */		\n\t		vmovaps		0x240(%%rcx),%%zmm9	/* zmm11 <-     a[jp+p2] */			\n\t"\
		"vmovaps	%%zmm0		,%%zmm2	/* zmm2 <- cpy a[jt+p4] */		\n\t		vmovaps		%%zmm8 	,%%zmm10	/* zmm10 <- cpy a[jt+p2] */			\n\t"\
		"vmovaps	%%zmm1		,%%zmm3	/* zmm3 <- cpy a[jp+p4] */		\n\t		vmovaps		%%zmm9 	,%%zmm11	/* zmm11 <- cpy a[jp+p2] */			\n\t"\
		/***************************************************************************/\
		/*** From hereon, things are identical to the code in radix16_dif_pass: ****/\
		/***************************************************************************/\
		"vmovaps		  (%%rdx),%%zmm16	\n\t	vmovaps		     (%%rdi),%%zmm20	\n\t"\
		"vmovaps	 0x080(%%rdx),%%zmm17	\n\t	vmovaps		0x080(%%rdi),%%zmm21	\n\t"\
		"vmovaps	 0x040(%%rdx),%%zmm18	\n\t	vmovaps		0x040(%%rdi),%%zmm22	\n\t"\
		"vmovaps	 0x0c0(%%rdx),%%zmm19	\n\t	vmovaps		0x0c0(%%rdi),%%zmm23	\n\t"\
		"vmovaps	0x180(%%rcx),%%zmm4			/* zmm4 <-     a[jt+p12] */	\n\t		vmovaps		0x300(%%rcx),%%zmm12			/* zmm12 <-     a[jt+p10] */\n\t"\
		"vmovaps	0x1c0(%%rcx),%%zmm5			/* zmm5 <-     a[jp+p12] */	\n\t		vmovaps		0x340(%%rcx),%%zmm13			/* zmm13 <-     a[jp+p10] */\n\t"\
		"vmovaps	%%zmm4		,%%zmm6			/* zmm6 <- cpy a[jt+p12] */	\n\t		vmovaps			%%zmm12	,%%zmm14			/* zmm14 <- cpy a[jt+p10] */\n\t"\
		"vmovaps	%%zmm5		,%%zmm7			/* zmm7 <- cpy a[jp+p12] */	\n\t		vmovaps			%%zmm13	,%%zmm15			/* zmm15 <- cpy a[jp+p10] */\n\t"\
		"vmulpd			%%zmm16,%%zmm0,%%zmm0		/* A: a[jt+ p4]* c4 */	\n\t		vmulpd			%%zmm20,%%zmm8 ,%%zmm8 	\n\t"/* E: a[jt+ p2]* c2 */\
		"vmulpd			%%zmm16,%%zmm1,%%zmm1		/* B: a[jp+ p4]* c4 */	\n\t		vmulpd			%%zmm20,%%zmm9 ,%%zmm9 	\n\t"/* F: a[jp+ p2]* c2 */\
		"vmulpd			%%zmm17,%%zmm4,%%zmm4		/* C: a[jt+p12]*c12 */	\n\t		vmulpd			%%zmm21,%%zmm12,%%zmm12	\n\t"/* G: a[jt+p10]*c10 */\
		"vmulpd			%%zmm17,%%zmm5,%%zmm5		/* D: a[jp+p12]*c12 */	\n\t		vmulpd			%%zmm21,%%zmm13,%%zmm13	\n\t"/* H: a[jp+p10]*c10 */\
	"vfnmadd231pd		%%zmm18,%%zmm3,%%zmm0		/* A -= a[jp+ p4]* s4 */\n\t	vfnmadd231pd		%%zmm22,%%zmm11,%%zmm8 	\n\t"/* E -= a[jp+ p2]* s2 */\
	" vfmadd231pd		%%zmm18,%%zmm2,%%zmm1		/* B += a[jt+ p4]* s4 */\n\t	 vfmadd231pd		%%zmm22,%%zmm10,%%zmm9 	\n\t"/* F += a[jt+ p2]* s2 */\
	"vfnmadd231pd		%%zmm19,%%zmm7,%%zmm4		/* C -= a[jp+p12]*s12 */\n\t	vfnmadd231pd		%%zmm23,%%zmm15,%%zmm12	\n\t"/* G -= a[jp+p10]*s10 */\
	" vfmadd231pd		%%zmm19,%%zmm6,%%zmm5		/* D += a[jt+p12]*s12 */\n\t	 vfmadd231pd		%%zmm23,%%zmm14,%%zmm13	\n\t"/* H += a[jt+p10]*s10 */\
		/* Since zmm11 is last-used oreg, place zmm9 vopy in memory and instead use zmm11 to store 1.0 needed by FMAs-in-place-of-ADD/SUB: */\
		"vmovaps	%%zmm1		,%%zmm3		/* zmm3 <- cpy t6 */			\n\t		vmovaps		%%zmm9 		,%%zmm11		/* zmm11 <- cpy t14 */\n\t"\
		"vmovaps	%%zmm0		,%%zmm2		/* zmm2 <- cpy t5 */			\n\t		vmovaps		%%zmm8 		,%%zmm10		/* zmm10 <- cpy t13 */\n\t"\
		"vfmadd132pd	%%zmm30,%%zmm4,%%zmm0		/* ~t5 <- t5 +rt */		\n\t		vfmadd132pd	%%zmm30,%%zmm12,%%zmm8  	/* ~t13<- t13+rt */	\n\t"\
		"vfmadd132pd	%%zmm30,%%zmm5,%%zmm1		/* ~t6 <- t6 +it */		\n\t		vfmadd132pd	%%zmm30,%%zmm13,%%zmm9  	/* ~t14<- t14+it */	\n\t"\
		"vfmsub132pd	%%zmm30,%%zmm4,%%zmm2		/* ~t7 <- t5 -rt */		\n\t		vfmsub132pd	%%zmm30,%%zmm12,%%zmm10 	/* ~t15<- t13-rt */	\n\t"\
		"vfmsub132pd	%%zmm30,%%zmm5,%%zmm3		/* ~t8 <- t6 -it */		\n\t		vfmsub132pd	%%zmm30,%%zmm13,%%zmm11 	/* ~t16<- t14-it	zmm12,13 free */\n\t"\
		"\n\t"\
		/* Now do the p0,8 combo: */												/* Do the p6,14 combo - do p14 first so registers come out in same order as for p2,10 */\
	"leaq	0x940(%%rcx),%%rdx	/* c8 */								\n\t		leaq	0xc40(%%rcx),%%rdi	/* c14 */\n\t"\
		"vmovaps		  (%%rdx),%%zmm16	\n\t	vmovaps		     (%%rdi),%%zmm20	\n\t"\
		"vmovaps	 0x040(%%rdx),%%zmm17	\n\t	vmovaps		0x040(%%rdi),%%zmm21	\n\t"\
		"vmovaps	0x100(%%rcx)	,%%zmm4		/* a[jt+p8 ] */				\n\t		vmovaps		0x380(%%rcx),%%zmm12		/* a[jt+p14] */				\n\t"\
		"vmovaps	0x140(%%rcx)	,%%zmm5		/* a[jp+p8 ] */				\n\t		vmovaps		0x3c0(%%rcx),%%zmm13		/* a[jp+p14] */				\n\t"\
		"vmovaps	%%zmm4		,%%zmm6	/* zmm6 <- cpy a[jt+p8] */			\n\t		vmovaps			%%zmm12	,%%zmm14		/* zmm14 <- cpy a[jt+p14] */\n\t"\
		"vmovaps	%%zmm5		,%%zmm7	/* zmm7 <- cpy a[jp+p8] */			\n\t		vmovaps			%%zmm13	,%%zmm15		/* zmm15 <- cpy a[jp+p14] */\n\t"\
		"vmulpd			%%zmm16,%%zmm4,%%zmm4		/* a[jt+p8]*c8 */		\n\t		vmulpd			%%zmm20,%%zmm12,%%zmm12		/* a[jt+p14]*c14 */			\n\t"\
		"vmulpd			%%zmm16,%%zmm5,%%zmm5		/* a[jp+p8]*c8 */		\n\t		vmulpd			%%zmm20,%%zmm13,%%zmm13		/* a[jp+p14]*c14 */			\n\t"\
	"vfnmadd231pd		%%zmm17,%%zmm7,%%zmm4		/* a[jp+p8]*s8 */		\n\t	vfnmadd231pd		%%zmm21,%%zmm15,%%zmm12		/* a[jp+p14]*s14 */			\n\t"\
	" vfmadd231pd		%%zmm17,%%zmm6,%%zmm5		/* a[jt+p8]*s8 */		\n\t	 vfmadd231pd		%%zmm21,%%zmm14,%%zmm13		/* a[jt+p14]*s14 */			\n\t"\
		"																				vmovaps		%%zmm13		,0x3c0(%%rcx)	/* Store it in t16*/		\n\t"\
		"																				vmovaps		%%zmm12		,0x380(%%rcx)	/* Store rt in t15*/		\n\t"\
		"																				subq	$0x080,%%rdi	/* c6  */	\n\t"\
		"																				vmovaps		     (%%rdi),%%zmm20	\n\t"\
		"																				vmovaps		0x040(%%rdi),%%zmm21	\n\t"\
		"vmovaps		 (%%rcx),%%zmm6		/* a[jt    ] */					\n\t		vmovaps		0x280(%%rcx),%%zmm12		/* a[jt+p6 ] */				\n\t"\
		"vmovaps	0x040(%%rcx),%%zmm7		/* a[jp    ] */					\n\t		vmovaps		0x2c0(%%rcx),%%zmm13		/* a[jp+p6 ] */				\n\t"\
		"																				vmovaps			%%zmm12	,%%zmm14		/* zmm14 <- cpy a[jt+p6] */		\n\t"\
		"																				vmovaps			%%zmm13	,%%zmm15		/* zmm15 <- cpy a[jp+p6] */		\n\t"\
		"vsubpd		%%zmm4,%%zmm6,%%zmm6	/* ~t3 <- t1 -rt */				\n\t		vmulpd			%%zmm20,%%zmm12,%%zmm12		/* a[jt+p6]*c6 */			\n\t"\
		"vsubpd		%%zmm5,%%zmm7,%%zmm7	/* ~t4 <- t2 -it */				\n\t		vmulpd			%%zmm20,%%zmm13,%%zmm13		/* a[jp+p6]*c6 */			\n\t"\
																				"	 vfmadd231pd		%%zmm21,%%zmm14,%%zmm13		/* a[jt+p6]*s6 */			\n\t"\
																				"	vfnmadd231pd		%%zmm21,%%zmm15,%%zmm12		/* a[jp+p6]*s6 */			\n\t"\
		" vfmadd132pd	%%zmm31,%%zmm6,%%zmm4	/* ~t1 <- t1 +rt */			\n\t		vmovaps		%%zmm13		,%%zmm15		/* zmm15 <- cpy t14*/			\n\t"\
		" vfmadd132pd	%%zmm31,%%zmm7,%%zmm5	/* ~t2 <- t2 +it */\n\t					vmovaps		%%zmm12		,%%zmm14		/* zmm14 <- cpy t13*/			\n\t"\
			"									/* zmm4,5 free */						vsubpd		0x380(%%rcx),%%zmm12,%%zmm12		/* ~t15<- t13-rt */			\n\t"\
			"																			vsubpd		0x3c0(%%rcx),%%zmm13,%%zmm13		/* ~t16<- t14-it */			\n\t"\
			"																			vaddpd		0x380(%%rcx),%%zmm14,%%zmm14		/* ~t13<- t13+rt */			\n\t"\
			"																			vaddpd		0x3c0(%%rcx),%%zmm15,%%zmm15		/* ~t14<- t14+it */			\n\t"\
		/* Finish radix-4 butterfly and store results into temporary-array slots: */\
		"vsubpd		%%zmm0,%%zmm4,%%zmm4	/*~t5 =t1 -t5 */			\n\t		vfmsub132pd	%%zmm30,%%zmm14,%%zmm8 	/*~t13*/						\n\t"\
		"vsubpd		%%zmm1,%%zmm5,%%zmm5	/*~t6 =t2 -t6 */			\n\t		vfmsub132pd	%%zmm30,%%zmm15,%%zmm9 	/*~t14*/						\n\t"\
		"vmovaps	%%zmm4		,0x100(%%rcx)	/* a[jt+p8 ] <- ~t5 */		\n\t		vmovaps		%%zmm8 		,0x300(%%rcx)	/* a[jt+p8 ] <- ~t13*/		\n\t"\
		"vmovaps	%%zmm5		,0x140(%%rcx)	/* a[jp+p8 ] <- ~t6 */		\n\t		vmovaps		%%zmm9 		,0x340(%%rcx)	/* a[jp+p8 ] <- ~t14*/		\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm4,%%zmm0		/* 2*t5 */				\n\t		vfmadd132pd	%%zmm31,%%zmm8 ,%%zmm14	/* 2*t13*/						\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm5,%%zmm1		/* 2*t6 */				\n\t		vfmadd132pd	%%zmm31,%%zmm9 ,%%zmm15	/* 2*t14*/						\n\t"\
		"vmovaps	%%zmm0		,     (%%rcx)	/* a[jt    ] <- ~t1 */		\n\t		vmovaps		%%zmm14		,0x200(%%rcx)	/* a[jt    ] <- ~t9 */		\n\t"\
		"vmovaps	%%zmm1		,0x040(%%rcx)	/* a[jp    ] <- ~t2 */		\n\t		vmovaps		%%zmm15		,0x240(%%rcx)	/* a[jp    ] <- ~t10*/		\n\t"\
		"\n\t"\
		"vsubpd		%%zmm3,%%zmm6,%%zmm6	/*~t3 =t3 -t8 */			\n\t		vfmsub132pd	%%zmm30,%%zmm13,%%zmm10	/*~t11*/				\n\t"\
		"vsubpd		%%zmm2,%%zmm7,%%zmm7	/*~t8 =t4 -t7 */			\n\t		vfmsub132pd	%%zmm30,%%zmm12,%%zmm11	/*~t16*/				\n\t"\
		"vmovaps	%%zmm6		,0x080(%%rcx)	/* a[jt+p4 ] <- ~t3 */		\n\t		vmovaps		%%zmm10		,0x280(%%rcx)	/* a[jt+p4 ] <- ~t11*/		\n\t"\
		"vmovaps	%%zmm7		,0x1c0(%%rcx)	/* a[jp+p12] <- ~t8 */		\n\t		vmovaps		%%zmm11		,0x3c0(%%rcx)	/* a[jp+p12] <- ~t16*/		\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm6,%%zmm3	/*~t7 =t3 +t8 */			\n\t		vfmadd132pd	%%zmm31,%%zmm10,%%zmm13			/*~t15*/				\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm7,%%zmm2	/*~t4 =t4 +t7 */			\n\t		vfmadd132pd	%%zmm31,%%zmm11,%%zmm12			/*~t12*/				\n\t"\
		"vmovaps	%%zmm3		,0x180(%%rcx)	/* a[jt+p12] <- ~t7 */		\n\t		vmovaps		%%zmm13		,0x380(%%rcx)	/* a[jt+p12] <- ~t15*/		\n\t"\
		"vmovaps	%%zmm2		,0x0c0(%%rcx)	/* a[jp+p4 ] <- ~t4 */		\n\t		vmovaps		%%zmm12		,0x2c0(%%rcx)	/* a[jp+p4 ] <- ~t12*/		\n\t"\
		"\n\t"\
	/*...Block 3: */															/*...Block 4: */\
	/*	SSE2_RADIX4_DIF_4TWIDDLE(r17,r21,r19,r23, r17,c1): */					/* SSE2_RADIX4_DIF_4TWIDDLE(r25,r29,r27,r31, r25,c3): */\
		/* Do the p0,p8 combo: */													/* Do the p0,p8 combo: */\
	"leaq	0xcc0(%%rcx),%%rbx	\n\t"/* c1 */	/* All __r and __c pointers incr by +0x200 in rcol w.r.to lcol: */\
	"addq	$0x400,%%rcx		\n\t"/* r17 */											/* c3, r25 */\
		"vmovaps		 (%%rcx),%%zmm0		/* a[jt   ] */					\n\t		vmovaps		0x200(%%rcx),%%zmm8 		/* a[jt    ] */\n\t"\
		"vmovaps	0x040(%%rcx),%%zmm1		/* a[jp   ] */					\n\t		vmovaps		0x240(%%rcx),%%zmm9 		/* a[jp    ] */\n\t"\
		"vmovaps	0x100(%%rcx),%%zmm4			/* zmm4 <-     a[jt+p12] */	\n\t		vmovaps		0x300(%%rcx),%%zmm12			/* zmm12 <-     a[jt+p10] */\n\t"\
		"vmovaps	0x140(%%rcx),%%zmm5			/* zmm5 <-     a[jp+p12] */	\n\t		vmovaps		0x340(%%rcx),%%zmm13			/* zmm13 <-     a[jp+p10] */\n\t"\
		"vmovaps	%%zmm0		,%%zmm2		/* zmm2 <- cpy a[jt   ] */		\n\t		vmovaps		%%zmm8 		,%%zmm10		/* zmm10 <- cpy a[jt   ] */\n\t"\
		"vmovaps	%%zmm1		,%%zmm3		/* zmm3 <- cpy a[jp   ] */		\n\t		vmovaps		%%zmm9 		,%%zmm11		/* zmm11 <- cpy a[jp   ] */\n\t"\
		"vmovaps	%%zmm4		,%%zmm6			/* zmm6 <- cpy a[jt+p12] */	\n\t		vmovaps			%%zmm12	,%%zmm14			/* zmm14 <- cpy a[jt+p10] */\n\t"\
		"vmovaps	%%zmm5		,%%zmm7			/* zmm7 <- cpy a[jp+p12] */	\n\t		vmovaps			%%zmm13	,%%zmm15			/* zmm15 <- cpy a[jp+p10] */\n\t"\
		"vmulpd		     (%%rbx),%%zmm0,%%zmm0		/* A: a[jt+ p4]* c4 */	\n\t		vmulpd		0x200(%%rbx),%%zmm8 ,%%zmm8 	\n\t"/* E: a[jt+ p2]* c2 */\
		"vmulpd		     (%%rbx),%%zmm1,%%zmm1		/* B: a[jp+ p4]* c4 */	\n\t		vmulpd		0x200(%%rbx),%%zmm9 ,%%zmm9 	\n\t"/* F: a[jp+ p2]* c2 */\
		"vmulpd		0x080(%%rbx),%%zmm4,%%zmm4		/* C: a[jt+p12]*c12 */	\n\t		vmulpd		0x280(%%rbx),%%zmm12,%%zmm12	\n\t"/* G: a[jt+p10]*c10 */\
		"vmulpd		0x080(%%rbx),%%zmm5,%%zmm5		/* D: a[jp+p12]*c12 */	\n\t		vmulpd		0x280(%%rbx),%%zmm13,%%zmm13	\n\t"/* H: a[jp+p10]*c10 */\
	"vfnmadd231pd	0x040(%%rbx),%%zmm3,%%zmm0		/* A -= a[jp+ p4]* s4 */\n\t	vfnmadd231pd	0x240(%%rbx),%%zmm11,%%zmm8 	\n\t"/* E -= a[jp+ p2]* s2 */\
	" vfmadd231pd	0x040(%%rbx),%%zmm2,%%zmm1		/* B += a[jt+ p4]* s4 */\n\t	 vfmadd231pd	0x240(%%rbx),%%zmm10,%%zmm9 	\n\t"/* F += a[jt+ p2]* s2 */\
	"vfnmadd231pd	0x0c0(%%rbx),%%zmm7,%%zmm4		/* C -= a[jp+p12]*s12 */\n\t	vfnmadd231pd	0x2c0(%%rbx),%%zmm15,%%zmm12	\n\t"/* G -= a[jp+p10]*s10 */\
	" vfmadd231pd	0x0c0(%%rbx),%%zmm6,%%zmm5		/* D += a[jt+p12]*s12 */\n\t	 vfmadd231pd	0x2c0(%%rbx),%%zmm14,%%zmm13	\n\t"/* H += a[jt+p10]*s10 */\
		"vmovaps	%%zmm0		,%%zmm2		/* zmm2 <- cpy t5 */			\n\t		vmovaps		%%zmm8 		,%%zmm10		/* zmm10 <- cpy t9 */		\n\t"\
		"vmovaps	%%zmm1		,%%zmm3		/* zmm3 <- cpy t6 */			\n\t		vmovaps		%%zmm9 		,%%zmm11		/* zmm11 <- cpy t10*/		\n\t"\
		"vfmadd132pd	%%zmm30,%%zmm4,%%zmm0	/* ~t1 <- t1 +rt */			\n\t		vfmadd132pd	%%zmm30,%%zmm12,%%zmm8 		/* ~t1 <- t1 +rt */\n\t"\
		"vfmadd132pd	%%zmm30,%%zmm5,%%zmm1	/* ~t2 <- t2 +it */			\n\t		vfmadd132pd	%%zmm30,%%zmm13,%%zmm9 		/* ~t2 <- t2 +it */\n\t"\
		"vfmsub132pd	%%zmm30,%%zmm4,%%zmm2	/* ~t3 <- t1 -rt */			\n\t		vfmsub132pd	%%zmm30,%%zmm12,%%zmm10		/* ~t3 <- t1 -rt */\n\t"\
		"vfmsub132pd	%%zmm30,%%zmm5,%%zmm3	/* ~t4 <- t2 -it zmm4,5 free*/\n\t		vfmsub132pd	%%zmm30,%%zmm13,%%zmm11		/* ~t4 <- t2 -it	zmm12,5 free */\n\t"\
		"\n\t"\
	/* Do the p4,12 combo: */														/* Do the p4,12 combo: */\
	"addq	$0x180 ,%%rbx	\n\t"/* c13 */											/* c15 */\
		"vmovaps		  (%%rbx),%%zmm16	\n\t	vmovaps		0x200(%%rbx),%%zmm20	\n\t"\
		"vmovaps	 0x040(%%rbx),%%zmm17	\n\t	vmovaps		0x240(%%rbx),%%zmm21	\n\t"\
		"vmovaps	0x180(%%rcx),%%zmm6		/* a[jt+p12] */					\n\t		vmovaps		0x380(%%rcx),%%zmm14		/* a[jt+p12] */\n\t"\
		"vmovaps	0x1c0(%%rcx),%%zmm7		/* a[jp+p12] */					\n\t		vmovaps		0x3c0(%%rcx),%%zmm15		/* a[jp+p12] */\n\t"\
		"vmovaps	%%zmm6		,%%zmm4		/* zmm4 <- cpy a[jt+p12] */		\n\t		vmovaps		%%zmm14		,%%zmm12		/* zmm12 <- cpy a[jt+p12] */\n\t"\
		"vmovaps	%%zmm7		,%%zmm5		/* zmm5 <- cpy a[jp+p12] */		\n\t		vmovaps		%%zmm15		,%%zmm13		/* zmm13 <- cpy a[jp+p12] */\n\t"\
		"vmulpd			%%zmm16,%%zmm4,%%zmm4		/* a[jt+p12]*c12 */		\n\t		vmulpd			%%zmm20,%%zmm12,%%zmm12		/* a[jt+p12]*c12 */\n\t"\
		"vmulpd			%%zmm16,%%zmm5,%%zmm5		/* a[jp+p12]*c12 */		\n\t		vmulpd			%%zmm20,%%zmm13,%%zmm13		/* a[jp+p12]*c12 */\n\t"\
	"vfnmadd231pd		%%zmm17,%%zmm7,%%zmm4		/* a[jp+p12]*s12 */		\n\t	vfnmadd231pd		%%zmm21,%%zmm15,%%zmm12		/* a[jp+p12]*s12 */\n\t"\
	" vfmadd231pd		%%zmm17,%%zmm6,%%zmm5		/* a[jt+p12]*s12 */		\n\t	 vfmadd231pd		%%zmm21,%%zmm14,%%zmm13		/* a[jt+p12]*s12 */\n\t"\
		"vmovaps	%%zmm5		,0x040(%%rcx)	/* store it */				\n\t		vmovaps		%%zmm13		,0x240(%%rcx)	/* store it */\n\t"\
		"vmovaps	%%zmm4		,     (%%rcx)	/* store rt */				\n\t		vmovaps		%%zmm12		,0x200(%%rcx)	/* store rt */\n\t"\
		"\n\t"\
	"subq	$0x080 ,%%rbx	\n\t"/* c5 */											/* c7  */\
		"vmovaps		  (%%rbx),%%zmm16	\n\t	vmovaps		0x200(%%rbx),%%zmm20	\n\t"\
		"vmovaps	 0x040(%%rbx),%%zmm17	\n\t	vmovaps		0x240(%%rbx),%%zmm21	\n\t"\
		"vmovaps	0x080(%%rcx),%%zmm4		/* a[jt+p4] */					\n\t		vmovaps		0x280(%%rcx),%%zmm12		/* a[jt+p4] */\n\t"\
		"vmovaps	0x0c0(%%rcx),%%zmm5		/* a[jp+p4] */					\n\t		vmovaps		0x2c0(%%rcx),%%zmm13		/* a[jp+p4] */\n\t"\
		"vmovaps		%%zmm4	,%%zmm6		/* zmm4 <- cpy a[jt+p4] */		\n\t		vmovaps			%%zmm12	,%%zmm14		/* zmm12 <- cpy a[jt+p4] */\n\t"\
		"vmovaps		%%zmm5	,%%zmm7		/* zmm5 <- cpy a[jp+p4] */		\n\t		vmovaps			%%zmm13	,%%zmm15		/* zmm13 <- cpy a[jp+p4] */\n\t"\
		"vmulpd			%%zmm16,%%zmm4,%%zmm4		/* a[jt+p4]*c4 */		\n\t		vmulpd			%%zmm20,%%zmm12,%%zmm12		/* a[jt+p4]*c4 */\n\t"\
		"vmulpd			%%zmm16,%%zmm5,%%zmm5		/* a[jp+p4]*c4 */		\n\t		vmulpd			%%zmm20,%%zmm13,%%zmm13		/* a[jp+p4]*c4 */\n\t"\
	"vfnmadd231pd		%%zmm17,%%zmm7,%%zmm4		/* a[jp+p4]*s4 */		\n\t	vfnmadd231pd		%%zmm21,%%zmm15,%%zmm12		/* a[jp+p4]*s4 */\n\t"\
	" vfmadd231pd		%%zmm17,%%zmm6,%%zmm5		/* a[jt+p4]*s4 */		\n\t	 vfmadd231pd		%%zmm21,%%zmm14,%%zmm13		/* a[jt+p4]*s4 */\n\t"\
		"vmovaps	%%zmm5		,%%zmm7		/* zmm7 <- cpy t6 */			\n\t		vmovaps		%%zmm13		,%%zmm15		/* zmm15 <- cpy t6 */\n\t"\
		"vmovaps	%%zmm4		,%%zmm6		/* zmm6 <- cpy t5 */			\n\t		vmovaps		%%zmm12		,%%zmm14		/* zmm14 <- cpy t5 */\n\t"\
		"vsubpd		     (%%rcx),%%zmm4,%%zmm4		/* ~t7 <- t5 -rt */		\n\t		vsubpd		0x200(%%rcx),%%zmm12,%%zmm12		/* ~t7 <- t5 -rt */\n\t"\
		"vsubpd		0x040(%%rcx),%%zmm5,%%zmm5		/* ~t8 <- t6 -it */		\n\t		vsubpd		0x240(%%rcx),%%zmm13,%%zmm13		/* ~t8 <- t6 -it */\n\t"\
		"vaddpd		     (%%rcx),%%zmm6,%%zmm6		/* ~t5 <- t5 +rt */		\n\t		vaddpd		0x200(%%rcx),%%zmm14,%%zmm14		/* ~t5 <- t5 +rt */\n\t"\
		"vaddpd		0x040(%%rcx),%%zmm7,%%zmm7		/* ~t6 <- t6 +it */		\n\t		vaddpd		0x240(%%rcx),%%zmm15,%%zmm15		/* ~t6 <- t6 +it */\n\t"\
		"\n\t"\
	/* Finish radix-4 butterfly and store results into temp-array slots: */\
		"vfmsub132pd	%%zmm30,%%zmm6,%%zmm0	/*~t5 */					\n\t		vfmsub132pd	%%zmm30,%%zmm14,%%zmm8 	/*~t5 */\n\t"\
		"vfmsub132pd	%%zmm30,%%zmm7,%%zmm1	/*~t6 */					\n\t		vfmsub132pd	%%zmm30,%%zmm15,%%zmm9 	/*~t6 */\n\t"\
		"vfmsub132pd	%%zmm30,%%zmm5,%%zmm2	/*~t3 */					\n\t		vfmsub132pd	%%zmm30,%%zmm13,%%zmm10	/*~t3 */\n\t"\
		"vfmsub132pd	%%zmm30,%%zmm4,%%zmm3	/*~t8 */					\n\t		vfmsub132pd	%%zmm30,%%zmm12,%%zmm11	/*~t8 */\n\t"\
		"vmovaps	%%zmm0,0x100(%%rcx)	/* a[jt+p8 ] <- ~t5 */				\n\t		vmovaps		%%zmm8 ,0x300(%%rcx)	/* a[jt+p8 ] <- ~t5 */\n\t"\
		"vmovaps	%%zmm1,0x140(%%rcx)	/* a[jp+p8 ] <- ~t6 */				\n\t		vmovaps		%%zmm9 ,0x340(%%rcx)	/* a[jp+p8 ] <- ~t6 */\n\t"\
		"vmovaps	%%zmm2,0x080(%%rcx)	/* a[jt+p4 ] <- ~t3 */				\n\t		vmovaps		%%zmm10,0x280(%%rcx)	/* a[jt+p4 ] <- ~t3 */\n\t"\
		"vmovaps	%%zmm3,0x1c0(%%rcx)	/* a[jp+p12] <- ~t8 */				\n\t		vmovaps		%%zmm11,0x3c0(%%rcx)	/* a[jp+p12] <- ~t8 */\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm0,%%zmm6	/*~t1 */					\n\t		vfmadd132pd	%%zmm31,%%zmm8 ,%%zmm14	/*~t1 */\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm1,%%zmm7	/*~t2 */					\n\t		vfmadd132pd	%%zmm31,%%zmm9 ,%%zmm15	/*~t2 */\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm2,%%zmm5	/*~t7 */					\n\t		vfmadd132pd	%%zmm31,%%zmm10,%%zmm13	/*~t7 */\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm3,%%zmm4	/*~t4 */					\n\t		vfmadd132pd	%%zmm31,%%zmm11,%%zmm12	/*~t4 */\n\t"\
		"vmovaps	%%zmm6,     (%%rcx)/* a[jt    ] <- ~t1 */				\n\t		vmovaps		%%zmm14,0x200(%%rcx)	/* a[jt    ] <- ~t1 */\n\t"\
		"vmovaps	%%zmm7,0x040(%%rcx)	/* a[jp    ] <- ~t2 */				\n\t		vmovaps		%%zmm15,0x240(%%rcx)	/* a[jp    ] <- ~t2 */\n\t"\
		"vmovaps	%%zmm5,0x180(%%rcx)	/* a[jt+p12] <- ~t7 */				\n\t		vmovaps		%%zmm13,0x380(%%rcx)	/* a[jt+p12] <- ~t7 */\n\t"\
		"vmovaps	%%zmm4,0x0c0(%%rcx)	/* a[jp+p4 ] <- ~t4 */				\n\t		vmovaps		%%zmm12,0x2c0(%%rcx)	/* a[jp+p4 ] <- ~t4 */\n\t"\
		"\n\t"\
	/**************************************************************************************/\
	/*...and now do four more radix-4 transforms, including the internal twiddle factors: */\
	/**************************************************************************************/\
		"\n\t"\
	"/*...Block 1: t1,9,17,25 */\n\t									/*...Block 3: t5,13,21,29: All rax-offsets incr +0x080 in rcol w.r.to lcol: */\n\t"\
		"movq	%[__r1],%%rax							\n\t				vmovaps		0x100(%%rax),%%zmm8 		/* t5  */\n\t"\
		"movq	%[__r9],%%rbx							\n\t				vmovaps		0x140(%%rax),%%zmm9 		/* t6  */\n\t"\
		"movq	%[__isrt2],%%rdi	\n\t"/* &two still in rsi */		"	vmovaps		0x340(%%rax),%%zmm11		/* t14 */\n\t"\
		"vmovaps		 (%%rax),%%zmm0		/* t1  */\n\t					vmovaps		0x300(%%rax),%%zmm10		/* t13 */\n\t"\
		"vmovaps	0x040(%%rax),%%zmm1		/* t2  */\n\t					vsubpd		%%zmm11,%%zmm8 ,%%zmm8 		/* t5 =t5 -t14 */\n\t"\
		"vmovaps	0x200(%%rax),%%zmm2		/* t9  */\n\t					vsubpd		%%zmm10,%%zmm9 ,%%zmm9 		/* t14=t6 -t13 */\n\t"\
		"vmovaps	0x240(%%rax),%%zmm3		/* t14 */\n\t					vfmadd132pd	%%zmm31,%%zmm8 ,%%zmm11	/* t13=t5 +t14 */\n\t"\
		"\n\t																vfmadd132pd	%%zmm31,%%zmm9 ,%%zmm10	/* t6 =t6 +t13 */\n\t"\
		"vsubpd		%%zmm2,%%zmm0,%%zmm0		/* t9 =t1 -t9  */\n\t		vmovaps		0x500(%%rax),%%zmm12		/* t21 */\n\t"\
		"vsubpd		%%zmm3,%%zmm1,%%zmm1		/* t14=t2 -t14 */\n\t		vmovaps		0x540(%%rax),%%zmm13		/* t22 */\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm0,%%zmm2		/* t1 =t1 +t9  */\n\t		vmovaps		0x700(%%rax),%%zmm14		/* t29 */\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm1,%%zmm3		/* t2 =t2 +t14 */\n\t		vmovaps		0x740(%%rax),%%zmm15		/* t30 */\n\t"\
		"\n\t																vsubpd		0x540(%%rax),%%zmm12,%%zmm12		/* t21-t22 */\n\t"\
		"vmovaps	0x400(%%rax),%%zmm4		/* t17 */\n\t					vaddpd		0x500(%%rax),%%zmm13,%%zmm13		/* t22+t21 */\n\t"\
		"vmovaps	0x440(%%rax),%%zmm5		/* t18 */\n\t					vmulpd		%%zmm29		,%%zmm12,%%zmm12	/* t21 = (t21-t22)*ISRT2 */\n\t"\
		"vmovaps	0x600(%%rax),%%zmm6		/* t25 */\n\t					vmulpd		%%zmm29		,%%zmm13,%%zmm13	/* t22 = (t22+t21)*ISRT2 */\n\t"\
		"vmovaps	0x640(%%rax),%%zmm7		/* t26 */\n\t					vaddpd		0x740(%%rax),%%zmm14,%%zmm14		/* t29+t30 */\n\t"\
		"\n\t																vsubpd		0x700(%%rax),%%zmm15,%%zmm15		/* t30-t29 */\n\t"\
																		"	vmulpd		%%zmm29,%%zmm14,%%zmm14	/*  rt = (t29+t30)*ISRT2 */\n\t"\
																		"	vmulpd		%%zmm29,%%zmm15,%%zmm15	/*  it = (t30-t29)*ISRT2 */\n\t"\
		"vsubpd		%%zmm6,%%zmm4,%%zmm4		/* t25=t17-t25 */\n\t		vsubpd		%%zmm14,%%zmm12,%%zmm12		/* t21=t21-rt */\n\t"\
		"vsubpd		%%zmm7,%%zmm5,%%zmm5		/* t26=t18-t26 */\n\t		vsubpd		%%zmm15,%%zmm13,%%zmm13		/* t22=t22-it */\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm4,%%zmm6		/* t17=t17+t25 */\n\t		vfmadd132pd	%%zmm31,%%zmm12,%%zmm14		/* t29=t21+rt */\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm5,%%zmm7		/* t18=t18+t26 */\n\t		vfmadd132pd	%%zmm31,%%zmm13,%%zmm15		/* t30=t22+it */\n\t"\
		"movq	%[__r17],%%rcx						\n\t		movq		%[__r25],%%rdx			\n\t"\
		"vsubpd		%%zmm6,%%zmm2,%%zmm2						\n\t		vsubpd		%%zmm12,%%zmm8 ,%%zmm8 			\n\t"\
		"vsubpd		%%zmm7,%%zmm3,%%zmm3						\n\t		vsubpd		%%zmm13,%%zmm10,%%zmm10			\n\t"\
		"vmovaps		%%zmm2,     (%%rcx)						\n\t		vmovaps		%%zmm8 ,0x100(%%rcx)		\n\t"\
		"vmovaps		%%zmm3,0x040(%%rcx)						\n\t		vmovaps		%%zmm10,0x140(%%rcx)		\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm2,%%zmm6						\n\t		vfmadd132pd	%%zmm31,%%zmm8 ,%%zmm12			\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm3,%%zmm7						\n\t		vfmadd132pd	%%zmm31,%%zmm10,%%zmm13			\n\t"\
		"vmovaps		%%zmm6,     (%%rax)						\n\t		vmovaps		%%zmm12,0x100(%%rax)		\n\t"\
		"vmovaps		%%zmm7,0x040(%%rax)						\n\t		vmovaps		%%zmm13,0x140(%%rax)		\n\t"\
		"vsubpd		%%zmm5,%%zmm0,%%zmm0						\n\t		vsubpd		%%zmm15,%%zmm11,%%zmm11			\n\t"\
		"vsubpd		%%zmm4,%%zmm1,%%zmm1						\n\t		vsubpd		%%zmm14,%%zmm9 ,%%zmm9 			\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm0,%%zmm5						\n\t		vfmadd132pd	%%zmm31,%%zmm11,%%zmm15			\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm1,%%zmm4						\n\t		vfmadd132pd	%%zmm31,%%zmm9 ,%%zmm14			\n\t"\
		"vmovaps		%%zmm0,     (%%rbx)						\n\t		vmovaps		%%zmm11,0x100(%%rbx)		\n\t"\
		"vmovaps		%%zmm1,0x040(%%rdx)						\n\t		vmovaps		%%zmm9 ,0x140(%%rdx)		\n\t"\
		"vmovaps		%%zmm5,     (%%rdx)						\n\t		vmovaps		%%zmm15,0x100(%%rdx)		\n\t"\
		"vmovaps		%%zmm4,0x040(%%rbx)						\n\t		vmovaps		%%zmm14,0x140(%%rbx)		\n\t"\
		"\n\t"\
	/*...Block 2: t3,11,19,27 */\
		"addq	$0x080,%%rax		\n\t"/* r3  */\
		"addq	$0x080,%%rbx							\n\t"\
		"addq	$0x080,%%rcx							\n\t"\
		"addq	$0x080,%%rdx							\n\t"\
																			/*...Block 4: t7,15,23,31 */\
		"vmovaps	0x400(%%rax),%%zmm4		/* t19 */		\n\t		vmovaps		0x500(%%rax),%%zmm12		/* t23 */\n\t"\
		"vmovaps	0x440(%%rax),%%zmm5		/* t20 */		\n\t		vmovaps		0x540(%%rax),%%zmm13		/* t24 */\n\t"\
		"vmovaps	0x600(%%rax),%%zmm6	/* t27 */			\n\t		vmovaps		0x700(%%rax),%%zmm14		/* t31 */\n\t"\
		"vmovaps	0x640(%%rax),%%zmm7	/* t28 */			\n\t		vmovaps		0x740(%%rax),%%zmm15		/* t32 */\n\t"\
		"vmovaps	%%zmm4,%%zmm0		/* copy t19 */		\n\t		vmovaps		%%zmm12,%%zmm8 		/* copy t23 */\n\t"\
		"vmovaps	%%zmm5,%%zmm1		/* copy t20 */		\n\t		vmovaps		%%zmm13,%%zmm9 		/* copy t24 */\n\t"\
	/*	"vmovaps	%%zmm6,%%zmm2	/x copy t27 x/			\n\t		vmovaps		%%zmm14,%%zmm10		/x copy t31 x/\n\t" */\
		"vmovaps	%%zmm7,%%zmm3	/* copy t28 */			\n\t		vmovaps		%%zmm15,%%zmm11		/* copy t32 */\n\t"\
		"\n\t"\
	"vmovaps	0x040(%%rdi),%%zmm2	\n\t	vmovaps	0x080(%%rdi),%%zmm10	\n\t"/* cc0,ss0 data shared by both columns */\
		"vmulpd			%%zmm2 ,%%zmm4,%%zmm4	/* t19*c */	\n\t		vmulpd			%%zmm10,%%zmm12,%%zmm12		/* t23*s */\n\t"\
		"vmulpd			%%zmm10,%%zmm6,%%zmm6	/* t27*s */	\n\t		vmulpd			%%zmm2 ,%%zmm14,%%zmm14		/* t31*c */\n\t"\
		"vmulpd			%%zmm2 ,%%zmm5,%%zmm5	/* t20*c */	\n\t		vmulpd			%%zmm10,%%zmm13,%%zmm13		/* t24*s */\n\t"\
		"vmulpd			%%zmm10,%%zmm7,%%zmm7	/* t28*s */	\n\t		vmulpd			%%zmm2 ,%%zmm15,%%zmm15		/* t32*c */\n\t"\
	"vfnmadd231pd		%%zmm10,%%zmm1,%%zmm4	/* ~t19 */	\n\t	vfnmadd231pd		%%zmm2 ,%%zmm9 ,%%zmm12		/* ~t23 */\n\t"\
	"vfnmadd231pd		%%zmm2 ,%%zmm3,%%zmm6	/* rt */	\n\t	vfnmadd231pd		%%zmm10,%%zmm11,%%zmm14		/* rt */\n\t"\
	" vfmadd231pd		%%zmm10,%%zmm0,%%zmm5	/* ~t20 */	\n\t	 vfmadd231pd		%%zmm2 ,%%zmm8 ,%%zmm13		/* ~t24 */\n\t"\
	" vfmadd231pd	0x600(%%rax),%%zmm2,%%zmm7	/* it */	\n\t	 vfmadd231pd	0x700(%%rax),%%zmm10,%%zmm15	/* it */\n\t"\
		"\n\t"\
		"vsubpd		%%zmm6,%%zmm4,%%zmm4		/*~t27=t19-rt */\n\t		vsubpd		%%zmm14,%%zmm12,%%zmm12		/*~t23=t23-rt */\n\t"\
		"vsubpd		%%zmm7,%%zmm5,%%zmm5		/*~t28=t20-it */\n\t		vsubpd		%%zmm15,%%zmm13,%%zmm13		/*~t24=t24-it */\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm4,%%zmm6		/*~t19=t19+rt */\n\t	vfmadd132pd	%%zmm31,%%zmm12,%%zmm14		/*~t31=t23+rt */\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm5,%%zmm7		/*~t20=t20+it */\n\t	vfmadd132pd	%%zmm31,%%zmm13,%%zmm15		/*~t32=t24+it */\n\t"\
		"\n\t"\
		"vmovaps	0x200(%%rax),%%zmm2		/* t11 */\n\t					vmovaps		0x300(%%rax),%%zmm10		/* t15 */\n\t"\
		"vmovaps	0x240(%%rax),%%zmm3		/* t12 */\n\t					vmovaps		0x340(%%rax),%%zmm11		/* t16 */\n\t"\
		"vsubpd		0x240(%%rax),%%zmm2,%%zmm2		/* t11-t12 */\n\t		vaddpd		0x340(%%rax),%%zmm10,%%zmm10	/* t15+t16 */\n\t"\
		"vaddpd		0x200(%%rax),%%zmm3,%%zmm3		/* t12+t11 */\n\t		vsubpd		0x300(%%rax),%%zmm11,%%zmm11	/* t16-t15 */\n\t"\
		"vmulpd		%%zmm29,%%zmm2,%%zmm2	/* rt = (t11-t12)*ISRT2 */\n\t	vmulpd		%%zmm29		,%%zmm10,%%zmm10	/* rt = (t15+t16)*ISRT2 */\n\t"\
		"vmulpd		%%zmm29,%%zmm3,%%zmm3	/* it = (t12+t11)*ISRT2 */\n\t	vmulpd		%%zmm29		,%%zmm11,%%zmm11	/* it = (t16-t15)*ISRT2 */\n\t"\
		"\n\t"\
		"vmovaps	     (%%rax),%%zmm0		/* t3  */\n\t					vmovaps		0x100(%%rax),%%zmm8 		/* t7  */\n\t"\
		"vmovaps	0x040(%%rax),%%zmm1		/* t4  */\n\t					vmovaps		0x140(%%rax),%%zmm9 		/* t8  */\n\t"\
		"\n\t"\
		"vsubpd		      %%zmm2,%%zmm0,%%zmm0		/*~t11=t3 -rt */\n\t	vsubpd		     %%zmm10,%%zmm8 ,%%zmm8 	/*~t7 =t7 -rt */\n\t"\
		"vsubpd		      %%zmm3,%%zmm1,%%zmm1		/*~t12=t4 -it */\n\t	vsubpd		     %%zmm11,%%zmm9 ,%%zmm9 	/*~t8 =t8 -it */\n\t"\
		"vaddpd		     (%%rax),%%zmm2,%%zmm2		/*~t3 =rt +t3 */\n\t	vaddpd		0x100(%%rax),%%zmm10,%%zmm10	/*~t15=rt +t7 */\n\t"\
		"vaddpd		0x040(%%rax),%%zmm3,%%zmm3		/*~t4 =it +t4 */\n\t	vaddpd		0x140(%%rax),%%zmm11,%%zmm11	/*~t16=it +t8 */\n\t"\
		"\n\t"\
		"vsubpd		%%zmm6,%%zmm2,%%zmm2		/* t3 -t19 */\n\t			vsubpd		%%zmm12		,%%zmm8 ,%%zmm8 	/* t7 -t23 */\n\t"\
		"vsubpd		%%zmm7,%%zmm3,%%zmm3		/* t4 -t20 */\n\t			vsubpd		%%zmm13		,%%zmm9 ,%%zmm9 	/* t8 -t24 */\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm2,%%zmm6		/* t3 +t19 */\n\t		vfmadd132pd	%%zmm31,%%zmm8,%%zmm12		/* t7 +t23 */\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm3,%%zmm7		/* t4 +t20 */\n\t		vfmadd132pd	%%zmm31,%%zmm9,%%zmm13		/* t8 +t24 */\n\t"\
		"vmovaps		%%zmm2,     (%%rcx)						\n\t		vmovaps		%%zmm8 ,0x100(%%rcx)		\n\t"\
		"vmovaps		%%zmm3,0x040(%%rcx)						\n\t		vmovaps		%%zmm9 ,0x140(%%rcx)		\n\t"\
		"vmovaps		%%zmm6,     (%%rax)						\n\t		vmovaps		%%zmm12,0x100(%%rax)		\n\t"\
		"vmovaps		%%zmm7,0x040(%%rax)						\n\t		vmovaps		%%zmm13,0x140(%%rax)		\n\t"\
		"vsubpd		%%zmm5,%%zmm0,%%zmm0						\n\t		vsubpd		%%zmm15,%%zmm10,%%zmm10			\n\t"\
		"vsubpd		%%zmm4,%%zmm1,%%zmm1						\n\t		vsubpd		%%zmm14,%%zmm11,%%zmm11			\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm0,%%zmm5	\n\t				vfmadd132pd	%%zmm31,%%zmm10,%%zmm15		\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm1,%%zmm4	\n\t				vfmadd132pd	%%zmm31,%%zmm11,%%zmm14		\n\t"\
		"vmovaps		%%zmm0,     (%%rbx)						\n\t		vmovaps		%%zmm10,0x100(%%rbx)		\n\t"\
		"vmovaps		%%zmm1,0x040(%%rdx)						\n\t		vmovaps		%%zmm11,0x140(%%rdx)		\n\t"\
		"vmovaps		%%zmm5,     (%%rdx)						\n\t		vmovaps		%%zmm15,0x100(%%rdx)		\n\t"\
		"vmovaps		%%zmm4,0x040(%%rbx)						\n\t		vmovaps		%%zmm14,0x140(%%rbx)		\n\t"\
		 :					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__add1] "m" (Xadd1)\
		 ,[__add2] "m" (Xadd2)\
		 ,[__add3] "m" (Xadd3)\
		 ,[__add4] "m" (Xadd4)\
		 ,[__add5] "m" (Xadd5)\
		 ,[__add6] "m" (Xadd6)\
		 ,[__add7] "m" (Xadd7)\
		 ,[__r1] "m" (Xr1)\
		 ,[__r9] "m" (Xr9)\
		 ,[__r17] "m" (Xr17)\
		 ,[__r25] "m" (Xr25)\
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__cc0] "m" (Xcc0)\
		 ,[__c1] "m" (Xc1)\
		 ,[__c2] "m" (Xc2)\
		 ,[__c3] "m" (Xc3)\
		 ,[__c4] "m" (Xc4)\
		 ,[__c5] "m" (Xc5)\
		 ,[__c6] "m" (Xc6)\
		 ,[__c7] "m" (Xc7)\
		 ,[__c8] "m" (Xc8)\
		 ,[__c9] "m" (Xc9)\
		 ,[__c10] "m" (Xc10)\
		 ,[__c11] "m" (Xc11)\
		 ,[__c12] "m" (Xc12)\
		 ,[__c13] "m" (Xc13)\
		 ,[__c14] "m" (Xc14)\
		 ,[__c15] "m" (Xc15)\
	 ,[__pfetch_dist] "m" (Xpfetch_dist)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r10","r11","r12","r13","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15","xmm16","xmm17","xmm18","xmm19","xmm20","xmm21","xmm22","xmm23", "xmm29","xmm30","xmm31"	/* Clobbered registers */\
	);\
	}

  #else	// AVX-512 version:

	#define SSE2_RADIX16_WRAPPER_DIF(Xadd0,Xadd1,Xadd2,Xadd3,Xadd4,Xadd5,Xadd6,Xadd7,Xr1,Xr9,Xr17,Xr25,Xisrt2,Xcc0,Xc1,Xc2,Xc3,Xc4,Xc5,Xc6,Xc7,Xc8,Xc9,Xc10,Xc11,Xc12,Xc13,Xc14,Xc15,Xpfetch_dist)\
	{\
	__asm__ volatile (\
		"movq	%[__add0],%%rax	\n\t"\
		"movq	%[__r1]  ,%%rsi	\n\t"\
		"leaq	0x10c0(%%rsi),%%r10		\n\t"/* &one */\
		"vmovaps	0x800(%%rsi),%%zmm29	\n\t"/* isrt2 */\
		"vmovaps	     (%%r10),%%zmm30	\n\t"/* 1.0 */\
		"vmovaps	0x040(%%r10),%%zmm31	\n\t"/* 2.0 */\
	"movslq	%[__pfetch_dist],%%r14	\n\t"\
	"leaq	(%%rax,%%r14,8),%%r14	\n\t"	/* Block 0 [base-address + data-fetch-ahead index] */\
	"prefetcht1	(%%r14)\n\t"\
	/*************************************************************/\
	/*                  1st set of inputs:                       */\
	/*************************************************************/\
		"movq	%[__add1],%%rbx	\n\t"\
		"movq	%[__add2],%%rcx	\n\t"\
		"movq	%[__add3],%%rdx	\n\t"\
		"movq	%[__add4],%%r10	\n\t"\
		"movq	%[__add5],%%r11	\n\t"\
		"movq	%[__add6],%%r12	\n\t"\
		"movq	%[__add7],%%r13	\n\t"\
	/**** Start with 8-way interleaving - Cf. radix-32 wrapper-DFT macros for commented versions of in-register shuffle code: ****/\
	/* a[j+p0]: Inputs from add0+[0,0x100,0x200,0x300,0x400,0x500,0x600,0x700]. Outputs into local store at r1+[same byte offsets]: */\
		"vmovaps 		(%%rax),%%zmm0			\n\t	vmovaps		0x40(%%rax),%%zmm10	\n\t"\
		"vmovaps 		(%%rbx),%%zmm1			\n\t	vmovaps		0x40(%%rbx),%%zmm11	\n\t"\
		"vmovaps 		(%%rcx),%%zmm2			\n\t	vmovaps		0x40(%%rcx),%%zmm12	\n\t"\
		"vmovaps 		(%%rdx),%%zmm3			\n\t	vmovaps		0x40(%%rdx),%%zmm13	\n\t"\
		"vmovaps 		(%%r10),%%zmm4			\n\t	vmovaps		0x40(%%r10),%%zmm14	\n\t"\
		"vmovaps 		(%%r11),%%zmm5			\n\t	vmovaps		0x40(%%r11),%%zmm15	\n\t"\
		"vmovaps 		(%%r12),%%zmm6			\n\t	vmovaps		0x40(%%r12),%%zmm16	\n\t"\
		"vmovaps 		(%%r13),%%zmm7			\n\t	vmovaps		0x40(%%r13),%%zmm17	\n\t"\
		"\n\t"\
		"vunpcklpd		%%zmm1,%%zmm0,%%zmm8	\n\t	vunpcklpd		%%zmm11,%%zmm10,%%zmm9 		\n\t"\
		"vunpckhpd		%%zmm1,%%zmm0,%%zmm1	\n\t	vunpckhpd		%%zmm11,%%zmm10,%%zmm11		\n\t"\
		"vunpcklpd		%%zmm3,%%zmm2,%%zmm0	\n\t	vunpcklpd		%%zmm13,%%zmm12,%%zmm10		\n\t"\
		"vunpckhpd		%%zmm3,%%zmm2,%%zmm3	\n\t	vunpckhpd		%%zmm13,%%zmm12,%%zmm13		\n\t"\
		"vunpcklpd		%%zmm5,%%zmm4,%%zmm2	\n\t	vunpcklpd		%%zmm15,%%zmm14,%%zmm12		\n\t"\
		"vunpckhpd		%%zmm5,%%zmm4,%%zmm5	\n\t	vunpckhpd		%%zmm15,%%zmm14,%%zmm15		\n\t"\
		"vunpcklpd		%%zmm7,%%zmm6,%%zmm4	\n\t	vunpcklpd		%%zmm17,%%zmm16,%%zmm14		\n\t"\
		"vunpckhpd		%%zmm7,%%zmm6,%%zmm7	\n\t	vunpckhpd		%%zmm17,%%zmm16,%%zmm17		\n\t"\
		"\n\t"\
		"vshuff64x2	$136,%%zmm0,%%zmm8,%%zmm6	\n\t	vshuff64x2	$136,%%zmm10,%%zmm9 ,%%zmm16	\n\t"\
		"vshuff64x2	$221,%%zmm0,%%zmm8,%%zmm0	\n\t	vshuff64x2	$221,%%zmm10,%%zmm9 ,%%zmm10	\n\t"\
		"vshuff64x2	$136,%%zmm3,%%zmm1,%%zmm8	\n\t	vshuff64x2	$136,%%zmm13,%%zmm11,%%zmm9 	\n\t"\
		"vshuff64x2	$221,%%zmm3,%%zmm1,%%zmm3	\n\t	vshuff64x2	$221,%%zmm13,%%zmm11,%%zmm13	\n\t"\
		"vshuff64x2	$136,%%zmm4,%%zmm2,%%zmm1	\n\t	vshuff64x2	$136,%%zmm14,%%zmm12,%%zmm11	\n\t"\
		"vshuff64x2	$221,%%zmm4,%%zmm2,%%zmm4	\n\t	vshuff64x2	$221,%%zmm14,%%zmm12,%%zmm14	\n\t"\
		"vshuff64x2	$136,%%zmm7,%%zmm5,%%zmm2	\n\t	vshuff64x2	$136,%%zmm17,%%zmm15,%%zmm12	\n\t"\
		"vshuff64x2	$221,%%zmm7,%%zmm5,%%zmm7	\n\t	vshuff64x2	$221,%%zmm17,%%zmm15,%%zmm17	\n\t"\
		"\n\t"\
		"vshuff64x2	$136,%%zmm1,%%zmm6,%%zmm5	\n\t	vshuff64x2	$136,%%zmm11,%%zmm16,%%zmm15	\n\t"\
		"vshuff64x2	$221,%%zmm1,%%zmm6,%%zmm1	\n\t	vshuff64x2	$221,%%zmm11,%%zmm16,%%zmm11	\n\t"\
		"vshuff64x2	$136,%%zmm2,%%zmm8,%%zmm6	\n\t	vshuff64x2	$136,%%zmm12,%%zmm9 ,%%zmm16	\n\t"\
		"vshuff64x2	$221,%%zmm2,%%zmm8,%%zmm2	\n\t	vshuff64x2	$221,%%zmm12,%%zmm9 ,%%zmm12	\n\t"\
		"vshuff64x2	$136,%%zmm4,%%zmm0,%%zmm8	\n\t	vshuff64x2	$136,%%zmm14,%%zmm10,%%zmm9 	\n\t"\
		"vshuff64x2	$221,%%zmm4,%%zmm0,%%zmm4	\n\t	vshuff64x2	$221,%%zmm14,%%zmm10,%%zmm14	\n\t"\
		"vshuff64x2	$136,%%zmm7,%%zmm3,%%zmm0	\n\t	vshuff64x2	$136,%%zmm17,%%zmm13,%%zmm10	\n\t"\
		"vshuff64x2	$221,%%zmm7,%%zmm3,%%zmm7	\n\t	vshuff64x2	$221,%%zmm17,%%zmm13,%%zmm17	\n\t"\
		/**** rows 0,1,4,5,8,9,c,d contain original cols 0,4,2,6,1,5,3,7: ****/\
		"vmovaps		%%zmm5,     (%%rsi)		\n\t	vmovaps	%%zmm15,0x040(%%rsi)	\n\t"\
		"vmovaps		%%zmm6,0x400(%%rsi)		\n\t	vmovaps	%%zmm16,0x440(%%rsi)	\n\t"\
		"vmovaps		%%zmm8,0x200(%%rsi)		\n\t	vmovaps	%%zmm9 ,0x240(%%rsi)	\n\t"\
		"vmovaps		%%zmm0,0x600(%%rsi)		\n\t	vmovaps	%%zmm10,0x640(%%rsi)	\n\t"\
		"vmovaps		%%zmm1,0x080(%%rsi)		\n\t	vmovaps	%%zmm11,0x0c0(%%rsi)	\n\t"\
		"vmovaps		%%zmm2,0x480(%%rsi)		\n\t	vmovaps	%%zmm12,0x4c0(%%rsi)	\n\t"\
		"vmovaps		%%zmm4,0x280(%%rsi)		\n\t	vmovaps	%%zmm14,0x2c0(%%rsi)	\n\t"\
		"vmovaps		%%zmm7,0x680(%%rsi)		\n\t	vmovaps	%%zmm17,0x6c0(%%rsi)	\n\t"\
		"\n\t"\
	/* a[j+p4]: Inputs from add0+[0,0x100,0x200,0x300,0x400,0x500,0x600,0x700]+0x80. Outputs into r1+[same byte offsets]+0x40: */\
	"addq	$0x80,%%rax	\n\t"\
	"addq	$0x80,%%rbx	\n\t"\
	"addq	$0x80,%%rcx	\n\t"\
	"addq	$0x80,%%rdx	\n\t"\
	"addq	$0x80,%%r10	\n\t"\
	"addq	$0x80,%%r11	\n\t"\
	"addq	$0x80,%%r12	\n\t"\
	"addq	$0x80,%%r13	\n\t"\
	"addq	$0x100,%%rsi	\n\t"\
		"vmovaps 		(%%rax),%%zmm0			\n\t	vmovaps		0x40(%%rax),%%zmm10	\n\t"\
		"vmovaps 		(%%rbx),%%zmm1			\n\t	vmovaps		0x40(%%rbx),%%zmm11	\n\t"\
		"vmovaps 		(%%rcx),%%zmm2			\n\t	vmovaps		0x40(%%rcx),%%zmm12	\n\t"\
		"vmovaps 		(%%rdx),%%zmm3			\n\t	vmovaps		0x40(%%rdx),%%zmm13	\n\t"\
		"vmovaps 		(%%r10),%%zmm4			\n\t	vmovaps		0x40(%%r10),%%zmm14	\n\t"\
		"vmovaps 		(%%r11),%%zmm5			\n\t	vmovaps		0x40(%%r11),%%zmm15	\n\t"\
		"vmovaps 		(%%r12),%%zmm6			\n\t	vmovaps		0x40(%%r12),%%zmm16	\n\t"\
		"vmovaps 		(%%r13),%%zmm7			\n\t	vmovaps		0x40(%%r13),%%zmm17	\n\t"\
		"\n\t"\
		"vunpcklpd		%%zmm1,%%zmm0,%%zmm8	\n\t	vunpcklpd		%%zmm11,%%zmm10,%%zmm9 		\n\t"\
		"vunpckhpd		%%zmm1,%%zmm0,%%zmm1	\n\t	vunpckhpd		%%zmm11,%%zmm10,%%zmm11		\n\t"\
		"vunpcklpd		%%zmm3,%%zmm2,%%zmm0	\n\t	vunpcklpd		%%zmm13,%%zmm12,%%zmm10		\n\t"\
		"vunpckhpd		%%zmm3,%%zmm2,%%zmm3	\n\t	vunpckhpd		%%zmm13,%%zmm12,%%zmm13		\n\t"\
		"vunpcklpd		%%zmm5,%%zmm4,%%zmm2	\n\t	vunpcklpd		%%zmm15,%%zmm14,%%zmm12		\n\t"\
		"vunpckhpd		%%zmm5,%%zmm4,%%zmm5	\n\t	vunpckhpd		%%zmm15,%%zmm14,%%zmm15		\n\t"\
		"vunpcklpd		%%zmm7,%%zmm6,%%zmm4	\n\t	vunpcklpd		%%zmm17,%%zmm16,%%zmm14		\n\t"\
		"vunpckhpd		%%zmm7,%%zmm6,%%zmm7	\n\t	vunpckhpd		%%zmm17,%%zmm16,%%zmm17		\n\t"\
		"\n\t"\
		"vshuff64x2	$136,%%zmm0,%%zmm8,%%zmm6	\n\t	vshuff64x2	$136,%%zmm10,%%zmm9 ,%%zmm16	\n\t"\
		"vshuff64x2	$221,%%zmm0,%%zmm8,%%zmm0	\n\t	vshuff64x2	$221,%%zmm10,%%zmm9 ,%%zmm10	\n\t"\
		"vshuff64x2	$136,%%zmm3,%%zmm1,%%zmm8	\n\t	vshuff64x2	$136,%%zmm13,%%zmm11,%%zmm9 	\n\t"\
		"vshuff64x2	$221,%%zmm3,%%zmm1,%%zmm3	\n\t	vshuff64x2	$221,%%zmm13,%%zmm11,%%zmm13	\n\t"\
		"vshuff64x2	$136,%%zmm4,%%zmm2,%%zmm1	\n\t	vshuff64x2	$136,%%zmm14,%%zmm12,%%zmm11	\n\t"\
		"vshuff64x2	$221,%%zmm4,%%zmm2,%%zmm4	\n\t	vshuff64x2	$221,%%zmm14,%%zmm12,%%zmm14	\n\t"\
		"vshuff64x2	$136,%%zmm7,%%zmm5,%%zmm2	\n\t	vshuff64x2	$136,%%zmm17,%%zmm15,%%zmm12	\n\t"\
		"vshuff64x2	$221,%%zmm7,%%zmm5,%%zmm7	\n\t	vshuff64x2	$221,%%zmm17,%%zmm15,%%zmm17	\n\t"\
		"\n\t"\
		"vshuff64x2	$136,%%zmm1,%%zmm6,%%zmm5	\n\t	vshuff64x2	$136,%%zmm11,%%zmm16,%%zmm15	\n\t"\
		"vshuff64x2	$221,%%zmm1,%%zmm6,%%zmm1	\n\t	vshuff64x2	$221,%%zmm11,%%zmm16,%%zmm11	\n\t"\
		"vshuff64x2	$136,%%zmm2,%%zmm8,%%zmm6	\n\t	vshuff64x2	$136,%%zmm12,%%zmm9 ,%%zmm16	\n\t"\
		"vshuff64x2	$221,%%zmm2,%%zmm8,%%zmm2	\n\t	vshuff64x2	$221,%%zmm12,%%zmm9 ,%%zmm12	\n\t"\
		"vshuff64x2	$136,%%zmm4,%%zmm0,%%zmm8	\n\t	vshuff64x2	$136,%%zmm14,%%zmm10,%%zmm9 	\n\t"\
		"vshuff64x2	$221,%%zmm4,%%zmm0,%%zmm4	\n\t	vshuff64x2	$221,%%zmm14,%%zmm10,%%zmm14	\n\t"\
		"vshuff64x2	$136,%%zmm7,%%zmm3,%%zmm0	\n\t	vshuff64x2	$136,%%zmm17,%%zmm13,%%zmm10	\n\t"\
		"vshuff64x2	$221,%%zmm7,%%zmm3,%%zmm7	\n\t	vshuff64x2	$221,%%zmm17,%%zmm13,%%zmm17	\n\t"\
		/**** rows 2,3,6,7,a,b,e,f contain original cols 8,c,a,e,9,d,b,f: ****/\
		"vmovaps		%%zmm5,     (%%rsi)		\n\t	vmovaps	%%zmm15,0x040(%%rsi)	\n\t"\
		"vmovaps		%%zmm6,0x400(%%rsi)		\n\t	vmovaps	%%zmm16,0x440(%%rsi)	\n\t"\
		"vmovaps		%%zmm8,0x200(%%rsi)		\n\t	vmovaps	%%zmm9 ,0x240(%%rsi)	\n\t"\
		"vmovaps		%%zmm0,0x600(%%rsi)		\n\t	vmovaps	%%zmm10,0x640(%%rsi)	\n\t"\
		"vmovaps		%%zmm1,0x080(%%rsi)		\n\t	vmovaps	%%zmm11,0x0c0(%%rsi)	\n\t"\
		"vmovaps		%%zmm2,0x480(%%rsi)		\n\t	vmovaps	%%zmm12,0x4c0(%%rsi)	\n\t"\
		"vmovaps		%%zmm4,0x280(%%rsi)		\n\t	vmovaps	%%zmm14,0x2c0(%%rsi)	\n\t"\
		"vmovaps		%%zmm7,0x680(%%rsi)		\n\t	vmovaps	%%zmm17,0x6c0(%%rsi)	\n\t"\
		"\n\t"\
	/* The data-patterning here is somewhat funky:
	Overall result rows 0-f contain original cols 0,4,8,c,2,6,a,e,1,5,9,d,3,7,b,f
	Bit-rev'd: row(br[0-f]) contain original cols 0,1,2,3,8,9,a,b,4,5,6,7,c,d,e,f, i.e. middle 2 quartets swapped.
	Ordering-by-cols we have that original main-array cols 0-f end up in local-mem rows 0,8,4,c,1,9,5,9,2,a,6,e,3,b,7,f.
	*/\
	/*****************/\
	/* Radix-16 DIF: */\
	/*****************/\
	"prefetcht1	0x100(%%r14)\n\t"\
	"movq	%[__add0],%%rax	\n\t"/* Use for FMA-related spills */\
	"movq	%[__r1] ,%%rcx	\n\t"\
"leaq	0x1100(%%rcx),%%rsi	\n\t"/* two */\
	/*...Block 1: */													/*...Block 2: */\
	"leaq	0xac0(%%rcx),%%rdi	/* c2 */\n\t"\
	"leaq	0x9c0(%%rcx),%%rdx	/* c4 */\n\t"\
		"vmovaps	0x080(%%rcx),%%zmm0	/* zmm0 <-     a[jt+p4] */		\n\t		vmovaps		0x200(%%rcx),%%zmm8	/* zmm10 <-     a[jt+p2] */			\n\t"\
		"vmovaps	0x0c0(%%rcx),%%zmm1	/* zmm1 <-     a[jp+p4] */		\n\t		vmovaps		0x240(%%rcx),%%zmm9	/* zmm11 <-     a[jp+p2] */			\n\t"\
		"vmovaps	%%zmm0		,%%zmm2	/* zmm2 <- cpy a[jt+p4] */		\n\t		vmovaps		%%zmm8 	,%%zmm10	/* zmm10 <- cpy a[jt+p2] */			\n\t"\
		"vmovaps	%%zmm1		,%%zmm3	/* zmm3 <- cpy a[jp+p4] */		\n\t		vmovaps		%%zmm9 	,%%zmm11	/* zmm11 <- cpy a[jp+p2] */			\n\t"\
		/***************************************************************************/\
		/*** From hereon, things are identical to the code in radix16_dif_pass: ****/\
		/***************************************************************************/\
		"vmovaps		  (%%rdx),%%zmm16	\n\t	vmovaps		     (%%rdi),%%zmm20	\n\t"\
		"vmovaps	 0x080(%%rdx),%%zmm17	\n\t	vmovaps		0x080(%%rdi),%%zmm21	\n\t"\
		"vmovaps	 0x040(%%rdx),%%zmm18	\n\t	vmovaps		0x040(%%rdi),%%zmm22	\n\t"\
		"vmovaps	 0x0c0(%%rdx),%%zmm19	\n\t	vmovaps		0x0c0(%%rdi),%%zmm23	\n\t"\
		"vmovaps	0x180(%%rcx),%%zmm4			/* zmm4 <-     a[jt+p12] */	\n\t		vmovaps		0x300(%%rcx),%%zmm12			/* zmm12 <-     a[jt+p10] */\n\t"\
		"vmovaps	0x1c0(%%rcx),%%zmm5			/* zmm5 <-     a[jp+p12] */	\n\t		vmovaps		0x340(%%rcx),%%zmm13			/* zmm13 <-     a[jp+p10] */\n\t"\
		"vmovaps	%%zmm4		,%%zmm6			/* zmm6 <- cpy a[jt+p12] */	\n\t		vmovaps			%%zmm12	,%%zmm14			/* zmm14 <- cpy a[jt+p10] */\n\t"\
		"vmovaps	%%zmm5		,%%zmm7			/* zmm7 <- cpy a[jp+p12] */	\n\t		vmovaps			%%zmm13	,%%zmm15			/* zmm15 <- cpy a[jp+p10] */\n\t"\
		"vmulpd			%%zmm16,%%zmm0,%%zmm0		/* A: a[jt+ p4]* c4 */	\n\t		vmulpd			%%zmm20,%%zmm8 ,%%zmm8 	\n\t"/* E: a[jt+ p2]* c2 */\
		"vmulpd			%%zmm16,%%zmm1,%%zmm1		/* B: a[jp+ p4]* c4 */	\n\t		vmulpd			%%zmm20,%%zmm9 ,%%zmm9 	\n\t"/* F: a[jp+ p2]* c2 */\
		"vmulpd			%%zmm17,%%zmm4,%%zmm4		/* C: a[jt+p12]*c12 */	\n\t		vmulpd			%%zmm21,%%zmm12,%%zmm12	\n\t"/* G: a[jt+p10]*c10 */\
		"vmulpd			%%zmm17,%%zmm5,%%zmm5		/* D: a[jp+p12]*c12 */	\n\t		vmulpd			%%zmm21,%%zmm13,%%zmm13	\n\t"/* H: a[jp+p10]*c10 */\
	"vfnmadd231pd		%%zmm18,%%zmm3,%%zmm0		/* A -= a[jp+ p4]* s4 */\n\t	vfnmadd231pd		%%zmm22,%%zmm11,%%zmm8 	\n\t"/* E -= a[jp+ p2]* s2 */\
	" vfmadd231pd		%%zmm18,%%zmm2,%%zmm1		/* B += a[jt+ p4]* s4 */\n\t	 vfmadd231pd		%%zmm22,%%zmm10,%%zmm9 	\n\t"/* F += a[jt+ p2]* s2 */\
	"vfnmadd231pd		%%zmm19,%%zmm7,%%zmm4		/* C -= a[jp+p12]*s12 */\n\t	vfnmadd231pd		%%zmm23,%%zmm15,%%zmm12	\n\t"/* G -= a[jp+p10]*s10 */\
	" vfmadd231pd		%%zmm19,%%zmm6,%%zmm5		/* D += a[jt+p12]*s12 */\n\t	 vfmadd231pd		%%zmm23,%%zmm14,%%zmm13	\n\t"/* H += a[jt+p10]*s10 */\
		/* Since zmm11 is last-used oreg, place zmm9 vopy in memory and instead use zmm11 to store 1.0 needed by FMAs-in-place-of-ADD/SUB: */\
		"vmovaps	%%zmm1		,%%zmm3		/* zmm3 <- cpy t6 */			\n\t		vmovaps		%%zmm9 		,%%zmm11		/* zmm11 <- cpy t14 */\n\t"\
		"vmovaps	%%zmm0		,%%zmm2		/* zmm2 <- cpy t5 */			\n\t		vmovaps		%%zmm8 		,%%zmm10		/* zmm10 <- cpy t13 */\n\t"\
		"vfmadd132pd	%%zmm30,%%zmm4,%%zmm0		/* ~t5 <- t5 +rt */		\n\t		vfmadd132pd	%%zmm30,%%zmm12,%%zmm8  	/* ~t13<- t13+rt */	\n\t"\
		"vfmadd132pd	%%zmm30,%%zmm5,%%zmm1		/* ~t6 <- t6 +it */		\n\t		vfmadd132pd	%%zmm30,%%zmm13,%%zmm9  	/* ~t14<- t14+it */	\n\t"\
		"vfmsub132pd	%%zmm30,%%zmm4,%%zmm2		/* ~t7 <- t5 -rt */		\n\t		vfmsub132pd	%%zmm30,%%zmm12,%%zmm10 	/* ~t15<- t13-rt */	\n\t"\
		"vfmsub132pd	%%zmm30,%%zmm5,%%zmm3		/* ~t8 <- t6 -it */		\n\t		vfmsub132pd	%%zmm30,%%zmm13,%%zmm11 	/* ~t16<- t14-it	zmm12,13 free */\n\t"\
	"prefetcht1	0x100(%%r14)\n\t"\
		"\n\t"\
		/* Now do the p0,8 combo: */												/* Do the p6,14 combo - do p14 first so registers come out in same order as for p2,10 */\
	"leaq	0x940(%%rcx),%%rdx	/* c8 */								\n\t		leaq	0xc40(%%rcx),%%rdi	/* c14 */\n\t"\
		"vmovaps		  (%%rdx),%%zmm16	\n\t	vmovaps		     (%%rdi),%%zmm20	\n\t"\
		"vmovaps	 0x040(%%rdx),%%zmm17	\n\t	vmovaps		0x040(%%rdi),%%zmm21	\n\t"\
		"vmovaps	0x100(%%rcx)	,%%zmm4		/* a[jt+p8 ] */				\n\t		vmovaps		0x380(%%rcx),%%zmm12		/* a[jt+p14] */				\n\t"\
		"vmovaps	0x140(%%rcx)	,%%zmm5		/* a[jp+p8 ] */				\n\t		vmovaps		0x3c0(%%rcx),%%zmm13		/* a[jp+p14] */				\n\t"\
		"vmovaps	%%zmm4		,%%zmm6	/* zmm6 <- cpy a[jt+p8] */			\n\t		vmovaps			%%zmm12	,%%zmm14		/* zmm14 <- cpy a[jt+p14] */\n\t"\
		"vmovaps	%%zmm5		,%%zmm7	/* zmm7 <- cpy a[jp+p8] */			\n\t		vmovaps			%%zmm13	,%%zmm15		/* zmm15 <- cpy a[jp+p14] */\n\t"\
		"vmulpd			%%zmm16,%%zmm4,%%zmm4		/* a[jt+p8]*c8 */		\n\t		vmulpd			%%zmm20,%%zmm12,%%zmm12		/* a[jt+p14]*c14 */			\n\t"\
		"vmulpd			%%zmm16,%%zmm5,%%zmm5		/* a[jp+p8]*c8 */		\n\t		vmulpd			%%zmm20,%%zmm13,%%zmm13		/* a[jp+p14]*c14 */			\n\t"\
	"vfnmadd231pd		%%zmm17,%%zmm7,%%zmm4		/* a[jp+p8]*s8 */		\n\t	vfnmadd231pd		%%zmm21,%%zmm15,%%zmm12		/* a[jp+p14]*s14 */			\n\t"\
	" vfmadd231pd		%%zmm17,%%zmm6,%%zmm5		/* a[jt+p8]*s8 */		\n\t	 vfmadd231pd		%%zmm21,%%zmm14,%%zmm13		/* a[jt+p14]*s14 */			\n\t"\
		"																				vmovaps		%%zmm13		,0x3c0(%%rcx)	/* Store it in t16*/		\n\t"\
		"																				vmovaps		%%zmm12		,0x380(%%rcx)	/* Store rt in t15*/		\n\t"\
		"																				subq	$0x080,%%rdi	/* c6  */	\n\t"\
		"																				vmovaps		     (%%rdi),%%zmm20	\n\t"\
		"																				vmovaps		0x040(%%rdi),%%zmm21	\n\t"\
		"vmovaps		 (%%rcx),%%zmm6		/* a[jt    ] */					\n\t		vmovaps		0x280(%%rcx),%%zmm12		/* a[jt+p6 ] */				\n\t"\
		"vmovaps	0x040(%%rcx),%%zmm7		/* a[jp    ] */					\n\t		vmovaps		0x2c0(%%rcx),%%zmm13		/* a[jp+p6 ] */				\n\t"\
		"																				vmovaps			%%zmm12	,%%zmm14		/* zmm14 <- cpy a[jt+p6] */		\n\t"\
		"																				vmovaps			%%zmm13	,%%zmm15		/* zmm15 <- cpy a[jp+p6] */		\n\t"\
		"vsubpd		%%zmm4,%%zmm6,%%zmm6	/* ~t3 <- t1 -rt */				\n\t		vmulpd			%%zmm20,%%zmm12,%%zmm12		/* a[jt+p6]*c6 */			\n\t"\
		"vsubpd		%%zmm5,%%zmm7,%%zmm7	/* ~t4 <- t2 -it */				\n\t		vmulpd			%%zmm20,%%zmm13,%%zmm13		/* a[jp+p6]*c6 */			\n\t"\
																				"	 vfmadd231pd		%%zmm21,%%zmm14,%%zmm13		/* a[jt+p6]*s6 */			\n\t"\
																				"	vfnmadd231pd		%%zmm21,%%zmm15,%%zmm12		/* a[jp+p6]*s6 */			\n\t"\
		" vfmadd132pd	%%zmm31,%%zmm6,%%zmm4	/* ~t1 <- t1 +rt */			\n\t		vmovaps		%%zmm13		,%%zmm15		/* zmm15 <- cpy t14*/			\n\t"\
		" vfmadd132pd	%%zmm31,%%zmm7,%%zmm5	/* ~t2 <- t2 +it */\n\t					vmovaps		%%zmm12		,%%zmm14		/* zmm14 <- cpy t13*/			\n\t"\
			"									/* zmm4,5 free */						vsubpd		0x380(%%rcx),%%zmm12,%%zmm12		/* ~t15<- t13-rt */			\n\t"\
			"																			vsubpd		0x3c0(%%rcx),%%zmm13,%%zmm13		/* ~t16<- t14-it */			\n\t"\
			"																			vaddpd		0x380(%%rcx),%%zmm14,%%zmm14		/* ~t13<- t13+rt */			\n\t"\
			"																			vaddpd		0x3c0(%%rcx),%%zmm15,%%zmm15		/* ~t14<- t14+it */			\n\t"\
		/* Finish radix-4 butterfly and store results into temporary-array slots: */\
		"vsubpd		%%zmm0,%%zmm4,%%zmm4	/*~t5 =t1 -t5 */			\n\t		vfmsub132pd	%%zmm30,%%zmm14,%%zmm8 	/*~t13*/						\n\t"\
		"vsubpd		%%zmm1,%%zmm5,%%zmm5	/*~t6 =t2 -t6 */			\n\t		vfmsub132pd	%%zmm30,%%zmm15,%%zmm9 	/*~t14*/						\n\t"\
		"vmovaps	%%zmm4		,0x100(%%rcx)	/* a[jt+p8 ] <- ~t5 */		\n\t		vmovaps		%%zmm8 		,0x300(%%rcx)	/* a[jt+p8 ] <- ~t13*/		\n\t"\
		"vmovaps	%%zmm5		,0x140(%%rcx)	/* a[jp+p8 ] <- ~t6 */		\n\t		vmovaps		%%zmm9 		,0x340(%%rcx)	/* a[jp+p8 ] <- ~t14*/		\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm4,%%zmm0		/* 2*t5 */				\n\t		vfmadd132pd	%%zmm31,%%zmm8 ,%%zmm14	/* 2*t13*/						\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm5,%%zmm1		/* 2*t6 */				\n\t		vfmadd132pd	%%zmm31,%%zmm9 ,%%zmm15	/* 2*t14*/						\n\t"\
		"vmovaps	%%zmm0		,     (%%rcx)	/* a[jt    ] <- ~t1 */		\n\t		vmovaps		%%zmm14		,0x200(%%rcx)	/* a[jt    ] <- ~t9 */		\n\t"\
		"vmovaps	%%zmm1		,0x040(%%rcx)	/* a[jp    ] <- ~t2 */		\n\t		vmovaps		%%zmm15		,0x240(%%rcx)	/* a[jp    ] <- ~t10*/		\n\t"\
		"\n\t"\
		"vsubpd		%%zmm3,%%zmm6,%%zmm6	/*~t3 =t3 -t8 */			\n\t		vfmsub132pd	%%zmm30,%%zmm13,%%zmm10	/*~t11*/				\n\t"\
		"vsubpd		%%zmm2,%%zmm7,%%zmm7	/*~t8 =t4 -t7 */			\n\t		vfmsub132pd	%%zmm30,%%zmm12,%%zmm11	/*~t16*/				\n\t"\
		"vmovaps	%%zmm6		,0x080(%%rcx)	/* a[jt+p4 ] <- ~t3 */		\n\t		vmovaps		%%zmm10		,0x280(%%rcx)	/* a[jt+p4 ] <- ~t11*/		\n\t"\
		"vmovaps	%%zmm7		,0x1c0(%%rcx)	/* a[jp+p12] <- ~t8 */		\n\t		vmovaps		%%zmm11		,0x3c0(%%rcx)	/* a[jp+p12] <- ~t16*/		\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm6,%%zmm3	/*~t7 =t3 +t8 */			\n\t		vfmadd132pd	%%zmm31,%%zmm10,%%zmm13			/*~t15*/				\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm7,%%zmm2	/*~t4 =t4 +t7 */			\n\t		vfmadd132pd	%%zmm31,%%zmm11,%%zmm12			/*~t12*/				\n\t"\
		"vmovaps	%%zmm3		,0x180(%%rcx)	/* a[jt+p12] <- ~t7 */		\n\t		vmovaps		%%zmm13		,0x380(%%rcx)	/* a[jt+p12] <- ~t15*/		\n\t"\
		"vmovaps	%%zmm2		,0x0c0(%%rcx)	/* a[jp+p4 ] <- ~t4 */		\n\t		vmovaps		%%zmm12		,0x2c0(%%rcx)	/* a[jp+p4 ] <- ~t12*/		\n\t"\
	"prefetcht1	0x100(%%r13)\n\t"\
		"\n\t"\
	/*...Block 3: */															/*...Block 4: */\
	/*	SSE2_RADIX4_DIF_4TWIDDLE(r17,r21,r19,r23, r17,c1): */					/* SSE2_RADIX4_DIF_4TWIDDLE(r25,r29,r27,r31, r25,c3): */\
		/* Do the p0,p8 combo: */													/* Do the p0,p8 combo: */\
	"leaq	0xcc0(%%rcx),%%rbx	\n\t"/* c1 */	/* All __r and __c pointers incr by +0x200 in rcol w.r.to lcol: */\
	"addq	$0x400,%%rcx		\n\t"/* r17 */											/* c3, r25 */\
		"vmovaps		 (%%rcx),%%zmm0		/* a[jt   ] */					\n\t		vmovaps		0x200(%%rcx),%%zmm8 		/* a[jt    ] */\n\t"\
		"vmovaps	0x040(%%rcx),%%zmm1		/* a[jp   ] */					\n\t		vmovaps		0x240(%%rcx),%%zmm9 		/* a[jp    ] */\n\t"\
		"vmovaps	0x100(%%rcx),%%zmm4			/* zmm4 <-     a[jt+p12] */	\n\t		vmovaps		0x300(%%rcx),%%zmm12			/* zmm12 <-     a[jt+p10] */\n\t"\
		"vmovaps	0x140(%%rcx),%%zmm5			/* zmm5 <-     a[jp+p12] */	\n\t		vmovaps		0x340(%%rcx),%%zmm13			/* zmm13 <-     a[jp+p10] */\n\t"\
		"vmovaps	%%zmm0		,%%zmm2		/* zmm2 <- cpy a[jt   ] */		\n\t		vmovaps		%%zmm8 		,%%zmm10		/* zmm10 <- cpy a[jt   ] */\n\t"\
		"vmovaps	%%zmm1		,%%zmm3		/* zmm3 <- cpy a[jp   ] */		\n\t		vmovaps		%%zmm9 		,%%zmm11		/* zmm11 <- cpy a[jp   ] */\n\t"\
		"vmovaps	%%zmm4		,%%zmm6			/* zmm6 <- cpy a[jt+p12] */	\n\t		vmovaps			%%zmm12	,%%zmm14			/* zmm14 <- cpy a[jt+p10] */\n\t"\
		"vmovaps	%%zmm5		,%%zmm7			/* zmm7 <- cpy a[jp+p12] */	\n\t		vmovaps			%%zmm13	,%%zmm15			/* zmm15 <- cpy a[jp+p10] */\n\t"\
		"vmulpd		     (%%rbx),%%zmm0,%%zmm0		/* A: a[jt+ p4]* c4 */	\n\t		vmulpd		0x200(%%rbx),%%zmm8 ,%%zmm8 	\n\t"/* E: a[jt+ p2]* c2 */\
		"vmulpd		     (%%rbx),%%zmm1,%%zmm1		/* B: a[jp+ p4]* c4 */	\n\t		vmulpd		0x200(%%rbx),%%zmm9 ,%%zmm9 	\n\t"/* F: a[jp+ p2]* c2 */\
		"vmulpd		0x080(%%rbx),%%zmm4,%%zmm4		/* C: a[jt+p12]*c12 */	\n\t		vmulpd		0x280(%%rbx),%%zmm12,%%zmm12	\n\t"/* G: a[jt+p10]*c10 */\
		"vmulpd		0x080(%%rbx),%%zmm5,%%zmm5		/* D: a[jp+p12]*c12 */	\n\t		vmulpd		0x280(%%rbx),%%zmm13,%%zmm13	\n\t"/* H: a[jp+p10]*c10 */\
	"vfnmadd231pd	0x040(%%rbx),%%zmm3,%%zmm0		/* A -= a[jp+ p4]* s4 */\n\t	vfnmadd231pd	0x240(%%rbx),%%zmm11,%%zmm8 	\n\t"/* E -= a[jp+ p2]* s2 */\
	" vfmadd231pd	0x040(%%rbx),%%zmm2,%%zmm1		/* B += a[jt+ p4]* s4 */\n\t	 vfmadd231pd	0x240(%%rbx),%%zmm10,%%zmm9 	\n\t"/* F += a[jt+ p2]* s2 */\
	"vfnmadd231pd	0x0c0(%%rbx),%%zmm7,%%zmm4		/* C -= a[jp+p12]*s12 */\n\t	vfnmadd231pd	0x2c0(%%rbx),%%zmm15,%%zmm12	\n\t"/* G -= a[jp+p10]*s10 */\
	" vfmadd231pd	0x0c0(%%rbx),%%zmm6,%%zmm5		/* D += a[jt+p12]*s12 */\n\t	 vfmadd231pd	0x2c0(%%rbx),%%zmm14,%%zmm13	\n\t"/* H += a[jt+p10]*s10 */\
		"vmovaps	%%zmm0		,%%zmm2		/* zmm2 <- cpy t5 */			\n\t		vmovaps		%%zmm8 		,%%zmm10		/* zmm10 <- cpy t9 */		\n\t"\
		"vmovaps	%%zmm1		,%%zmm3		/* zmm3 <- cpy t6 */			\n\t		vmovaps		%%zmm9 		,%%zmm11		/* zmm11 <- cpy t10*/		\n\t"\
		"vfmadd132pd	%%zmm30,%%zmm4,%%zmm0	/* ~t1 <- t1 +rt */			\n\t		vfmadd132pd	%%zmm30,%%zmm12,%%zmm8 		/* ~t1 <- t1 +rt */\n\t"\
		"vfmadd132pd	%%zmm30,%%zmm5,%%zmm1	/* ~t2 <- t2 +it */			\n\t		vfmadd132pd	%%zmm30,%%zmm13,%%zmm9 		/* ~t2 <- t2 +it */\n\t"\
		"vfmsub132pd	%%zmm30,%%zmm4,%%zmm2	/* ~t3 <- t1 -rt */			\n\t		vfmsub132pd	%%zmm30,%%zmm12,%%zmm10		/* ~t3 <- t1 -rt */\n\t"\
		"vfmsub132pd	%%zmm30,%%zmm5,%%zmm3	/* ~t4 <- t2 -it zmm4,5 free*/\n\t		vfmsub132pd	%%zmm30,%%zmm13,%%zmm11		/* ~t4 <- t2 -it	zmm12,5 free */\n\t"\
		"\n\t"\
	/* Do the p4,12 combo: */														/* Do the p4,12 combo: */\
	"addq	$0x180 ,%%rbx	\n\t"/* c13 */											/* c15 */\
		"vmovaps		  (%%rbx),%%zmm16	\n\t	vmovaps		0x200(%%rbx),%%zmm20	\n\t"\
		"vmovaps	 0x040(%%rbx),%%zmm17	\n\t	vmovaps		0x240(%%rbx),%%zmm21	\n\t"\
		"vmovaps	0x180(%%rcx),%%zmm6		/* a[jt+p12] */					\n\t		vmovaps		0x380(%%rcx),%%zmm14		/* a[jt+p12] */\n\t"\
		"vmovaps	0x1c0(%%rcx),%%zmm7		/* a[jp+p12] */					\n\t		vmovaps		0x3c0(%%rcx),%%zmm15		/* a[jp+p12] */\n\t"\
		"vmovaps	%%zmm6		,%%zmm4		/* zmm4 <- cpy a[jt+p12] */		\n\t		vmovaps		%%zmm14		,%%zmm12		/* zmm12 <- cpy a[jt+p12] */\n\t"\
		"vmovaps	%%zmm7		,%%zmm5		/* zmm5 <- cpy a[jp+p12] */		\n\t		vmovaps		%%zmm15		,%%zmm13		/* zmm13 <- cpy a[jp+p12] */\n\t"\
		"vmulpd			%%zmm16,%%zmm4,%%zmm4		/* a[jt+p12]*c12 */		\n\t		vmulpd			%%zmm20,%%zmm12,%%zmm12		/* a[jt+p12]*c12 */\n\t"\
		"vmulpd			%%zmm16,%%zmm5,%%zmm5		/* a[jp+p12]*c12 */		\n\t		vmulpd			%%zmm20,%%zmm13,%%zmm13		/* a[jp+p12]*c12 */\n\t"\
	"vfnmadd231pd		%%zmm17,%%zmm7,%%zmm4		/* a[jp+p12]*s12 */		\n\t	vfnmadd231pd		%%zmm21,%%zmm15,%%zmm12		/* a[jp+p12]*s12 */\n\t"\
	" vfmadd231pd		%%zmm17,%%zmm6,%%zmm5		/* a[jt+p12]*s12 */		\n\t	 vfmadd231pd		%%zmm21,%%zmm14,%%zmm13		/* a[jt+p12]*s12 */\n\t"\
		"vmovaps	%%zmm5		,0x040(%%rcx)	/* store it */				\n\t		vmovaps		%%zmm13		,0x240(%%rcx)	/* store it */\n\t"\
		"vmovaps	%%zmm4		,     (%%rcx)	/* store rt */				\n\t		vmovaps		%%zmm12		,0x200(%%rcx)	/* store rt */\n\t"\
		"\n\t"\
	"subq	$0x080 ,%%rbx	\n\t"/* c5 */											/* c7  */\
		"vmovaps		  (%%rbx),%%zmm16	\n\t	vmovaps		0x200(%%rbx),%%zmm20	\n\t"\
		"vmovaps	 0x040(%%rbx),%%zmm17	\n\t	vmovaps		0x240(%%rbx),%%zmm21	\n\t"\
		"vmovaps	0x080(%%rcx),%%zmm4		/* a[jt+p4] */					\n\t		vmovaps		0x280(%%rcx),%%zmm12		/* a[jt+p4] */\n\t"\
		"vmovaps	0x0c0(%%rcx),%%zmm5		/* a[jp+p4] */					\n\t		vmovaps		0x2c0(%%rcx),%%zmm13		/* a[jp+p4] */\n\t"\
		"vmovaps		%%zmm4	,%%zmm6		/* zmm4 <- cpy a[jt+p4] */		\n\t		vmovaps			%%zmm12	,%%zmm14		/* zmm12 <- cpy a[jt+p4] */\n\t"\
		"vmovaps		%%zmm5	,%%zmm7		/* zmm5 <- cpy a[jp+p4] */		\n\t		vmovaps			%%zmm13	,%%zmm15		/* zmm13 <- cpy a[jp+p4] */\n\t"\
		"vmulpd			%%zmm16,%%zmm4,%%zmm4		/* a[jt+p4]*c4 */		\n\t		vmulpd			%%zmm20,%%zmm12,%%zmm12		/* a[jt+p4]*c4 */\n\t"\
		"vmulpd			%%zmm16,%%zmm5,%%zmm5		/* a[jp+p4]*c4 */		\n\t		vmulpd			%%zmm20,%%zmm13,%%zmm13		/* a[jp+p4]*c4 */\n\t"\
	"vfnmadd231pd		%%zmm17,%%zmm7,%%zmm4		/* a[jp+p4]*s4 */		\n\t	vfnmadd231pd		%%zmm21,%%zmm15,%%zmm12		/* a[jp+p4]*s4 */\n\t"\
	" vfmadd231pd		%%zmm17,%%zmm6,%%zmm5		/* a[jt+p4]*s4 */		\n\t	 vfmadd231pd		%%zmm21,%%zmm14,%%zmm13		/* a[jt+p4]*s4 */\n\t"\
	"prefetcht1	0x200(%%r14)\n\t"\
		"vmovaps	%%zmm5		,%%zmm7		/* zmm7 <- cpy t6 */			\n\t		vmovaps		%%zmm13		,%%zmm15		/* zmm15 <- cpy t6 */\n\t"\
		"vmovaps	%%zmm4		,%%zmm6		/* zmm6 <- cpy t5 */			\n\t		vmovaps		%%zmm12		,%%zmm14		/* zmm14 <- cpy t5 */\n\t"\
		"vsubpd		     (%%rcx),%%zmm4,%%zmm4		/* ~t7 <- t5 -rt */		\n\t		vsubpd		0x200(%%rcx),%%zmm12,%%zmm12		/* ~t7 <- t5 -rt */\n\t"\
		"vsubpd		0x040(%%rcx),%%zmm5,%%zmm5		/* ~t8 <- t6 -it */		\n\t		vsubpd		0x240(%%rcx),%%zmm13,%%zmm13		/* ~t8 <- t6 -it */\n\t"\
		"vaddpd		     (%%rcx),%%zmm6,%%zmm6		/* ~t5 <- t5 +rt */		\n\t		vaddpd		0x200(%%rcx),%%zmm14,%%zmm14		/* ~t5 <- t5 +rt */\n\t"\
		"vaddpd		0x040(%%rcx),%%zmm7,%%zmm7		/* ~t6 <- t6 +it */		\n\t		vaddpd		0x240(%%rcx),%%zmm15,%%zmm15		/* ~t6 <- t6 +it */\n\t"\
		"\n\t"\
	/* Finish radix-4 butterfly and store results into temp-array slots: */\
		"vfmsub132pd	%%zmm30,%%zmm6,%%zmm0	/*~t5 */					\n\t		vfmsub132pd	%%zmm30,%%zmm14,%%zmm8 	/*~t5 */\n\t"\
		"vfmsub132pd	%%zmm30,%%zmm7,%%zmm1	/*~t6 */					\n\t		vfmsub132pd	%%zmm30,%%zmm15,%%zmm9 	/*~t6 */\n\t"\
		"vfmsub132pd	%%zmm30,%%zmm5,%%zmm2	/*~t3 */					\n\t		vfmsub132pd	%%zmm30,%%zmm13,%%zmm10	/*~t3 */\n\t"\
		"vfmsub132pd	%%zmm30,%%zmm4,%%zmm3	/*~t8 */					\n\t		vfmsub132pd	%%zmm30,%%zmm12,%%zmm11	/*~t8 */\n\t"\
		"vmovaps	%%zmm0,0x100(%%rcx)	/* a[jt+p8 ] <- ~t5 */				\n\t		vmovaps		%%zmm8 ,0x300(%%rcx)	/* a[jt+p8 ] <- ~t5 */\n\t"\
		"vmovaps	%%zmm1,0x140(%%rcx)	/* a[jp+p8 ] <- ~t6 */				\n\t		vmovaps		%%zmm9 ,0x340(%%rcx)	/* a[jp+p8 ] <- ~t6 */\n\t"\
		"vmovaps	%%zmm2,0x080(%%rcx)	/* a[jt+p4 ] <- ~t3 */				\n\t		vmovaps		%%zmm10,0x280(%%rcx)	/* a[jt+p4 ] <- ~t3 */\n\t"\
		"vmovaps	%%zmm3,0x1c0(%%rcx)	/* a[jp+p12] <- ~t8 */				\n\t		vmovaps		%%zmm11,0x3c0(%%rcx)	/* a[jp+p12] <- ~t8 */\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm0,%%zmm6	/*~t1 */					\n\t		vfmadd132pd	%%zmm31,%%zmm8 ,%%zmm14	/*~t1 */\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm1,%%zmm7	/*~t2 */					\n\t		vfmadd132pd	%%zmm31,%%zmm9 ,%%zmm15	/*~t2 */\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm2,%%zmm5	/*~t7 */					\n\t		vfmadd132pd	%%zmm31,%%zmm10,%%zmm13	/*~t7 */\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm3,%%zmm4	/*~t4 */					\n\t		vfmadd132pd	%%zmm31,%%zmm11,%%zmm12	/*~t4 */\n\t"\
		"vmovaps	%%zmm6,     (%%rcx)/* a[jt    ] <- ~t1 */				\n\t		vmovaps		%%zmm14,0x200(%%rcx)	/* a[jt    ] <- ~t1 */\n\t"\
		"vmovaps	%%zmm7,0x040(%%rcx)	/* a[jp    ] <- ~t2 */				\n\t		vmovaps		%%zmm15,0x240(%%rcx)	/* a[jp    ] <- ~t2 */\n\t"\
		"vmovaps	%%zmm5,0x180(%%rcx)	/* a[jt+p12] <- ~t7 */				\n\t		vmovaps		%%zmm13,0x380(%%rcx)	/* a[jt+p12] <- ~t7 */\n\t"\
		"vmovaps	%%zmm4,0x0c0(%%rcx)	/* a[jp+p4 ] <- ~t4 */				\n\t		vmovaps		%%zmm12,0x2c0(%%rcx)	/* a[jp+p4 ] <- ~t4 */\n\t"\
		"\n\t"\
	/**************************************************************************************/\
	/*...and now do four more radix-4 transforms, including the internal twiddle factors: */\
	/**************************************************************************************/\
		"\n\t"\
	"/*...Block 1: t1,9,17,25 */\n\t									/*...Block 3: t5,13,21,29: All rax-offsets incr +0x080 in rcol w.r.to lcol: */\n\t"\
		"movq	%[__r1],%%rax							\n\t				vmovaps		0x100(%%rax),%%zmm8 		/* t5  */\n\t"\
		"movq	%[__r9],%%rbx							\n\t				vmovaps		0x140(%%rax),%%zmm9 		/* t6  */\n\t"\
		"movq	%[__isrt2],%%rdi	\n\t"/* &two still in rsi */		"	vmovaps		0x340(%%rax),%%zmm11		/* t14 */\n\t"\
		"vmovaps		 (%%rax),%%zmm0		/* t1  */\n\t					vmovaps		0x300(%%rax),%%zmm10		/* t13 */\n\t"\
		"vmovaps	0x040(%%rax),%%zmm1		/* t2  */\n\t					vsubpd		%%zmm11,%%zmm8 ,%%zmm8 		/* t5 =t5 -t14 */\n\t"\
		"vmovaps	0x200(%%rax),%%zmm2		/* t9  */\n\t					vsubpd		%%zmm10,%%zmm9 ,%%zmm9 		/* t14=t6 -t13 */\n\t"\
		"vmovaps	0x240(%%rax),%%zmm3		/* t14 */\n\t					vfmadd132pd	%%zmm31,%%zmm8 ,%%zmm11	/* t13=t5 +t14 */\n\t"\
		"\n\t																vfmadd132pd	%%zmm31,%%zmm9 ,%%zmm10	/* t6 =t6 +t13 */\n\t"\
		"vsubpd		%%zmm2,%%zmm0,%%zmm0		/* t9 =t1 -t9  */\n\t		vmovaps		0x500(%%rax),%%zmm12		/* t21 */\n\t"\
		"vsubpd		%%zmm3,%%zmm1,%%zmm1		/* t14=t2 -t14 */\n\t		vmovaps		0x540(%%rax),%%zmm13		/* t22 */\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm0,%%zmm2		/* t1 =t1 +t9  */\n\t		vmovaps		0x700(%%rax),%%zmm14		/* t29 */\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm1,%%zmm3		/* t2 =t2 +t14 */\n\t		vmovaps		0x740(%%rax),%%zmm15		/* t30 */\n\t"\
		"\n\t																vsubpd		0x540(%%rax),%%zmm12,%%zmm12		/* t21-t22 */\n\t"\
		"vmovaps	0x400(%%rax),%%zmm4		/* t17 */\n\t					vaddpd		0x500(%%rax),%%zmm13,%%zmm13		/* t22+t21 */\n\t"\
		"vmovaps	0x440(%%rax),%%zmm5		/* t18 */\n\t					vmulpd		%%zmm29		,%%zmm12,%%zmm12	/* t21 = (t21-t22)*ISRT2 */\n\t"\
		"vmovaps	0x600(%%rax),%%zmm6		/* t25 */\n\t					vmulpd		%%zmm29		,%%zmm13,%%zmm13	/* t22 = (t22+t21)*ISRT2 */\n\t"\
		"vmovaps	0x640(%%rax),%%zmm7		/* t26 */\n\t					vaddpd		0x740(%%rax),%%zmm14,%%zmm14		/* t29+t30 */\n\t"\
		"\n\t																vsubpd		0x700(%%rax),%%zmm15,%%zmm15		/* t30-t29 */\n\t"\
	"prefetcht1	0x300(%%r14)\n\t"\
																		"	vmulpd		%%zmm29,%%zmm14,%%zmm14	/*  rt = (t29+t30)*ISRT2 */\n\t"\
																		"	vmulpd		%%zmm29,%%zmm15,%%zmm15	/*  it = (t30-t29)*ISRT2 */\n\t"\
		"vsubpd		%%zmm6,%%zmm4,%%zmm4		/* t25=t17-t25 */\n\t		vsubpd		%%zmm14,%%zmm12,%%zmm12		/* t21=t21-rt */\n\t"\
		"vsubpd		%%zmm7,%%zmm5,%%zmm5		/* t26=t18-t26 */\n\t		vsubpd		%%zmm15,%%zmm13,%%zmm13		/* t22=t22-it */\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm4,%%zmm6		/* t17=t17+t25 */\n\t		vfmadd132pd	%%zmm31,%%zmm12,%%zmm14		/* t29=t21+rt */\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm5,%%zmm7		/* t18=t18+t26 */\n\t		vfmadd132pd	%%zmm31,%%zmm13,%%zmm15		/* t30=t22+it */\n\t"\
		"movq	%[__r17],%%rcx						\n\t		movq		%[__r25],%%rdx			\n\t"\
		"vsubpd		%%zmm6,%%zmm2,%%zmm2						\n\t		vsubpd		%%zmm12,%%zmm8 ,%%zmm8 			\n\t"\
		"vsubpd		%%zmm7,%%zmm3,%%zmm3						\n\t		vsubpd		%%zmm13,%%zmm10,%%zmm10			\n\t"\
		"vmovaps		%%zmm2,     (%%rcx)						\n\t		vmovaps		%%zmm8 ,0x100(%%rcx)		\n\t"\
		"vmovaps		%%zmm3,0x040(%%rcx)						\n\t		vmovaps		%%zmm10,0x140(%%rcx)		\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm2,%%zmm6						\n\t		vfmadd132pd	%%zmm31,%%zmm8 ,%%zmm12			\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm3,%%zmm7						\n\t		vfmadd132pd	%%zmm31,%%zmm10,%%zmm13			\n\t"\
		"vmovaps		%%zmm6,     (%%rax)						\n\t		vmovaps		%%zmm12,0x100(%%rax)		\n\t"\
		"vmovaps		%%zmm7,0x040(%%rax)						\n\t		vmovaps		%%zmm13,0x140(%%rax)		\n\t"\
		"vsubpd		%%zmm5,%%zmm0,%%zmm0						\n\t		vsubpd		%%zmm15,%%zmm11,%%zmm11			\n\t"\
		"vsubpd		%%zmm4,%%zmm1,%%zmm1						\n\t		vsubpd		%%zmm14,%%zmm9 ,%%zmm9 			\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm0,%%zmm5						\n\t		vfmadd132pd	%%zmm31,%%zmm11,%%zmm15			\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm1,%%zmm4						\n\t		vfmadd132pd	%%zmm31,%%zmm9 ,%%zmm14			\n\t"\
		"vmovaps		%%zmm0,     (%%rbx)						\n\t		vmovaps		%%zmm11,0x100(%%rbx)		\n\t"\
		"vmovaps		%%zmm1,0x040(%%rdx)						\n\t		vmovaps		%%zmm9 ,0x140(%%rdx)		\n\t"\
		"vmovaps		%%zmm5,     (%%rdx)						\n\t		vmovaps		%%zmm15,0x100(%%rdx)		\n\t"\
		"vmovaps		%%zmm4,0x040(%%rbx)						\n\t		vmovaps		%%zmm14,0x140(%%rbx)		\n\t"\
		"\n\t"\
	/*...Block 2: t3,11,19,27 */\
		"addq	$0x080,%%rax		\n\t"/* r3  */\
		"addq	$0x080,%%rbx							\n\t"\
		"addq	$0x080,%%rcx							\n\t"\
		"addq	$0x080,%%rdx							\n\t"\
																			/*...Block 4: t7,15,23,31 */\
		"vmovaps	0x400(%%rax),%%zmm4		/* t19 */		\n\t		vmovaps		0x500(%%rax),%%zmm12		/* t23 */\n\t"\
		"vmovaps	0x440(%%rax),%%zmm5		/* t20 */		\n\t		vmovaps		0x540(%%rax),%%zmm13		/* t24 */\n\t"\
		"vmovaps	0x600(%%rax),%%zmm6	/* t27 */			\n\t		vmovaps		0x700(%%rax),%%zmm14		/* t31 */\n\t"\
		"vmovaps	0x640(%%rax),%%zmm7	/* t28 */			\n\t		vmovaps		0x740(%%rax),%%zmm15		/* t32 */\n\t"\
		"vmovaps	%%zmm4,%%zmm0		/* copy t19 */		\n\t		vmovaps		%%zmm12,%%zmm8 		/* copy t23 */\n\t"\
		"vmovaps	%%zmm5,%%zmm1		/* copy t20 */		\n\t		vmovaps		%%zmm13,%%zmm9 		/* copy t24 */\n\t"\
	/*	"vmovaps	%%zmm6,%%zmm2	/x copy t27 x/			\n\t		vmovaps		%%zmm14,%%zmm10		/x copy t31 x/\n\t" */\
		"vmovaps	%%zmm7,%%zmm3	/* copy t28 */			\n\t		vmovaps		%%zmm15,%%zmm11		/* copy t32 */\n\t"\
		"\n\t"\
	"vmovaps	0x040(%%rdi),%%zmm2	\n\t	vmovaps	0x080(%%rdi),%%zmm10	\n\t"/* cc0,ss0 data shared by both columns */\
		"vmulpd			%%zmm2 ,%%zmm4,%%zmm4	/* t19*c */	\n\t		vmulpd			%%zmm10,%%zmm12,%%zmm12		/* t23*s */\n\t"\
		"vmulpd			%%zmm10,%%zmm6,%%zmm6	/* t27*s */	\n\t		vmulpd			%%zmm2 ,%%zmm14,%%zmm14		/* t31*c */\n\t"\
		"vmulpd			%%zmm2 ,%%zmm5,%%zmm5	/* t20*c */	\n\t		vmulpd			%%zmm10,%%zmm13,%%zmm13		/* t24*s */\n\t"\
		"vmulpd			%%zmm10,%%zmm7,%%zmm7	/* t28*s */	\n\t		vmulpd			%%zmm2 ,%%zmm15,%%zmm15		/* t32*c */\n\t"\
	"vfnmadd231pd		%%zmm10,%%zmm1,%%zmm4	/* ~t19 */	\n\t	vfnmadd231pd		%%zmm2 ,%%zmm9 ,%%zmm12		/* ~t23 */\n\t"\
	"vfnmadd231pd		%%zmm2 ,%%zmm3,%%zmm6	/* rt */	\n\t	vfnmadd231pd		%%zmm10,%%zmm11,%%zmm14		/* rt */\n\t"\
	" vfmadd231pd		%%zmm10,%%zmm0,%%zmm5	/* ~t20 */	\n\t	 vfmadd231pd		%%zmm2 ,%%zmm8 ,%%zmm13		/* ~t24 */\n\t"\
	" vfmadd231pd	0x600(%%rax),%%zmm2,%%zmm7	/* it */	\n\t	 vfmadd231pd	0x700(%%rax),%%zmm10,%%zmm15	/* it */\n\t"\
	"prefetcht1	0x380(%%r14)\n\t"\
		"\n\t"\
		"vsubpd		%%zmm6,%%zmm4,%%zmm4		/*~t27=t19-rt */\n\t		vsubpd		%%zmm14,%%zmm12,%%zmm12		/*~t23=t23-rt */\n\t"\
		"vsubpd		%%zmm7,%%zmm5,%%zmm5		/*~t28=t20-it */\n\t		vsubpd		%%zmm15,%%zmm13,%%zmm13		/*~t24=t24-it */\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm4,%%zmm6		/*~t19=t19+rt */\n\t	vfmadd132pd	%%zmm31,%%zmm12,%%zmm14		/*~t31=t23+rt */\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm5,%%zmm7		/*~t20=t20+it */\n\t	vfmadd132pd	%%zmm31,%%zmm13,%%zmm15		/*~t32=t24+it */\n\t"\
		"\n\t"\
		"vmovaps	0x200(%%rax),%%zmm2		/* t11 */\n\t					vmovaps		0x300(%%rax),%%zmm10		/* t15 */\n\t"\
		"vmovaps	0x240(%%rax),%%zmm3		/* t12 */\n\t					vmovaps		0x340(%%rax),%%zmm11		/* t16 */\n\t"\
		"vsubpd		0x240(%%rax),%%zmm2,%%zmm2		/* t11-t12 */\n\t		vaddpd		0x340(%%rax),%%zmm10,%%zmm10	/* t15+t16 */\n\t"\
		"vaddpd		0x200(%%rax),%%zmm3,%%zmm3		/* t12+t11 */\n\t		vsubpd		0x300(%%rax),%%zmm11,%%zmm11	/* t16-t15 */\n\t"\
		"vmulpd		%%zmm29,%%zmm2,%%zmm2	/* rt = (t11-t12)*ISRT2 */\n\t	vmulpd		%%zmm29		,%%zmm10,%%zmm10	/* rt = (t15+t16)*ISRT2 */\n\t"\
		"vmulpd		%%zmm29,%%zmm3,%%zmm3	/* it = (t12+t11)*ISRT2 */\n\t	vmulpd		%%zmm29		,%%zmm11,%%zmm11	/* it = (t16-t15)*ISRT2 */\n\t"\
		"\n\t"\
		"vmovaps	     (%%rax),%%zmm0		/* t3  */\n\t					vmovaps		0x100(%%rax),%%zmm8 		/* t7  */\n\t"\
		"vmovaps	0x040(%%rax),%%zmm1		/* t4  */\n\t					vmovaps		0x140(%%rax),%%zmm9 		/* t8  */\n\t"\
		"\n\t"\
		"vsubpd		      %%zmm2,%%zmm0,%%zmm0		/*~t11=t3 -rt */\n\t	vsubpd		     %%zmm10,%%zmm8 ,%%zmm8 	/*~t7 =t7 -rt */\n\t"\
		"vsubpd		      %%zmm3,%%zmm1,%%zmm1		/*~t12=t4 -it */\n\t	vsubpd		     %%zmm11,%%zmm9 ,%%zmm9 	/*~t8 =t8 -it */\n\t"\
		"vaddpd		     (%%rax),%%zmm2,%%zmm2		/*~t3 =rt +t3 */\n\t	vaddpd		0x100(%%rax),%%zmm10,%%zmm10	/*~t15=rt +t7 */\n\t"\
		"vaddpd		0x040(%%rax),%%zmm3,%%zmm3		/*~t4 =it +t4 */\n\t	vaddpd		0x140(%%rax),%%zmm11,%%zmm11	/*~t16=it +t8 */\n\t"\
		"\n\t"\
		"vsubpd		%%zmm6,%%zmm2,%%zmm2		/* t3 -t19 */\n\t			vsubpd		%%zmm12		,%%zmm8 ,%%zmm8 	/* t7 -t23 */\n\t"\
		"vsubpd		%%zmm7,%%zmm3,%%zmm3		/* t4 -t20 */\n\t			vsubpd		%%zmm13		,%%zmm9 ,%%zmm9 	/* t8 -t24 */\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm2,%%zmm6		/* t3 +t19 */\n\t		vfmadd132pd	%%zmm31,%%zmm8,%%zmm12		/* t7 +t23 */\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm3,%%zmm7		/* t4 +t20 */\n\t		vfmadd132pd	%%zmm31,%%zmm9,%%zmm13		/* t8 +t24 */\n\t"\
		"vmovaps		%%zmm2,     (%%rcx)						\n\t		vmovaps		%%zmm8 ,0x100(%%rcx)		\n\t"\
		"vmovaps		%%zmm3,0x040(%%rcx)						\n\t		vmovaps		%%zmm9 ,0x140(%%rcx)		\n\t"\
		"vmovaps		%%zmm6,     (%%rax)						\n\t		vmovaps		%%zmm12,0x100(%%rax)		\n\t"\
		"vmovaps		%%zmm7,0x040(%%rax)						\n\t		vmovaps		%%zmm13,0x140(%%rax)		\n\t"\
		"vsubpd		%%zmm5,%%zmm0,%%zmm0						\n\t		vsubpd		%%zmm15,%%zmm10,%%zmm10			\n\t"\
		"vsubpd		%%zmm4,%%zmm1,%%zmm1						\n\t		vsubpd		%%zmm14,%%zmm11,%%zmm11			\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm0,%%zmm5	\n\t				vfmadd132pd	%%zmm31,%%zmm10,%%zmm15		\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm1,%%zmm4	\n\t				vfmadd132pd	%%zmm31,%%zmm11,%%zmm14		\n\t"\
		"vmovaps		%%zmm0,     (%%rbx)						\n\t		vmovaps		%%zmm10,0x100(%%rbx)		\n\t"\
		"vmovaps		%%zmm1,0x040(%%rdx)						\n\t		vmovaps		%%zmm11,0x140(%%rdx)		\n\t"\
		"vmovaps		%%zmm5,     (%%rdx)						\n\t		vmovaps		%%zmm15,0x100(%%rdx)		\n\t"\
		"vmovaps		%%zmm4,0x040(%%rbx)						\n\t		vmovaps		%%zmm14,0x140(%%rbx)		\n\t"\
		 :					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__add1] "m" (Xadd1)\
		 ,[__add2] "m" (Xadd2)\
		 ,[__add3] "m" (Xadd3)\
		 ,[__add4] "m" (Xadd4)\
		 ,[__add5] "m" (Xadd5)\
		 ,[__add6] "m" (Xadd6)\
		 ,[__add7] "m" (Xadd7)\
		 ,[__r1] "m" (Xr1)\
		 ,[__r9] "m" (Xr9)\
		 ,[__r17] "m" (Xr17)\
		 ,[__r25] "m" (Xr25)\
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__cc0] "m" (Xcc0)\
		 ,[__c1] "m" (Xc1)\
		 ,[__c2] "m" (Xc2)\
		 ,[__c3] "m" (Xc3)\
		 ,[__c4] "m" (Xc4)\
		 ,[__c5] "m" (Xc5)\
		 ,[__c6] "m" (Xc6)\
		 ,[__c7] "m" (Xc7)\
		 ,[__c8] "m" (Xc8)\
		 ,[__c9] "m" (Xc9)\
		 ,[__c10] "m" (Xc10)\
		 ,[__c11] "m" (Xc11)\
		 ,[__c12] "m" (Xc12)\
		 ,[__c13] "m" (Xc13)\
		 ,[__c14] "m" (Xc14)\
		 ,[__c15] "m" (Xc15)\
	 ,[__pfetch_dist] "m" (Xpfetch_dist)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r10","r11","r12","r13","r14","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15","xmm16","xmm17","xmm18","xmm19","xmm20","xmm21","xmm22","xmm23", "xmm29","xmm30","xmm31"	/* Clobbered registers */\
	);\
	}

  #endif	// (IMCI512 or AVX512?) toggle

  #ifdef USE_IMCI512	// 1st-gen Xeon Phi - Use modified 8x8 doubles-transpose algo [1a] from util.c:test_simd_transpose_8x8()

	#define SSE2_RADIX16_WRAPPER_DIT(Xadd0,Xadd1,Xadd2,Xadd3,Xadd4,Xadd5,Xadd6,Xadd7,Xr1,Xr9,Xr17,Xr25,Xisrt2,Xcc0,Xc1,Xc2,Xc3,Xc4,Xc5,Xc6,Xc7,Xc8,Xc9,Xc10,Xc11,Xc12,Xc13,Xc14,Xc15,Xpfetch_dist)\
	{\
	__asm__ volatile (\
	/* Cf. IMCI512 8x8 doubles-transpose algo [1a] in util.c ...
	do mask-register-setting first, so as to not clobber any GPRs used by compute section: */\
		"movl $0b10101010,%%ebx	\n\t movl $0b11001100,%%ecx	\n\t movl $0b11110000,%%edx	\n\t"\
		"kmov	%%ebx,%%k1		\n\t kmov	%%ecx,%%k3		\n\t kmov	%%edx,%%k5		\n\t"\
		"knot	%%k1 ,%%k2		\n\t knot	%%k3 ,%%k4		\n\t"\
	/*...Block 1: r1,9,17,25 */\
		"movq	%[__r1],%%rax					\n\t"\
		"leaq	0x10c0(%%rax),%%r10		\n\t"/* &one */\
		"vmovaps	0x800(%%rax),%%zmm29	\n\t"/* isrt2 */\
		"vmovaps	     (%%r10),%%zmm30	\n\t"/* 1.0 */\
		"vmovaps	0x040(%%r10),%%zmm31	\n\t"/* 2.0 */\
		"movq	%%rax,%%rbx						\n\t		leaq	0x1100(%%rax),%%rsi	\n\t"/* two */\
		"movq	%%rax,%%rcx						\n\t"	/*...Block 3: r3,11,19,27 */\
		"movq	%%rax,%%rdx						\n\t		leaq	0x080(%%rax),%%r10				\n\t"\
		"addq	$0x400,%%rbx					\n\t		leaq	0x080(%%rbx),%%r11				\n\t"\
		"addq	$0x200,%%rcx					\n\t		leaq	0x080(%%rcx),%%r12				\n\t"\
		"addq	$0x600,%%rdx					\n\t		leaq	0x080(%%rdx),%%r13				\n\t"\
		"/* SSE2_RADIX4_DIT_IN_PLACE(): */		\n\t		/* SSE2_RADIX4_DIT_IN_PLACE(): */		\n\t"\
		"vmovaps	     (%%rax),%%zmm2			\n\t		vmovaps	     (%%r10),%%zmm10			\n\t"\
		"vmovaps	0x040(%%rax),%%zmm3			\n\t		vmovaps	0x040(%%r10),%%zmm11			\n\t"\
		"vmovaps	     (%%rbx),%%zmm0			\n\t		vmovaps	     (%%r11),%%zmm8 			\n\t"\
		"vmovaps	0x040(%%rbx),%%zmm1			\n\t		vmovaps	0x040(%%r11),%%zmm9 			\n\t"\
		"vsubpd		%%zmm0,%%zmm2,%%zmm2		\n\t		vsubpd		%%zmm8 ,%%zmm10,%%zmm10		\n\t"\
		"vsubpd		%%zmm1,%%zmm3,%%zmm3		\n\t		vsubpd		%%zmm9 ,%%zmm11,%%zmm11		\n\t"\
		"vmovaps	     (%%rcx),%%zmm6			\n\t		vmovaps	     (%%r12),%%zmm14			\n\t"\
		"vmovaps	0x040(%%rcx),%%zmm7			\n\t		vmovaps	0x040(%%r12),%%zmm15			\n\t"\
		"vmovaps	     (%%rdx),%%zmm4			\n\t		vmovaps	     (%%r13),%%zmm12			\n\t"\
		"vmovaps	0x040(%%rdx),%%zmm5			\n\t		vmovaps	0x040(%%r13),%%zmm13			\n\t"\
		"vsubpd		%%zmm4,%%zmm6,%%zmm6		\n\t		vsubpd		 %%zmm12,%%zmm14,%%zmm14	\n\t"\
		"vsubpd		%%zmm5,%%zmm7,%%zmm7		\n\t		vsubpd	0x040(%%r13),%%zmm15,%%zmm15	\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm2,%%zmm0		\n\t	vfmadd132pd		%%zmm31,%%zmm10,%%zmm8 		\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm3,%%zmm1		\n\t	vfmadd132pd		%%zmm31,%%zmm11,%%zmm9 		\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm6,%%zmm4		\n\t	vfmadd132pd		%%zmm31,%%zmm14,%%zmm12		\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm7,%%zmm5		\n\t	vfmadd132pd		%%zmm31,%%zmm15,%%zmm13		\n\t"\
		"vsubpd	      %%zmm4,%%zmm0,%%zmm0		\n\t		vsubpd	      %%zmm12,%%zmm8 ,%%zmm8 	\n\t"\
		"vsubpd	      %%zmm5,%%zmm1,%%zmm1		\n\t		vsubpd	      %%zmm13,%%zmm9 ,%%zmm9 	\n\t"\
		"vsubpd	      %%zmm7,%%zmm2,%%zmm2		\n\t		vsubpd	      %%zmm15,%%zmm10,%%zmm10	\n\t"\
		"vsubpd	      %%zmm6,%%zmm3,%%zmm3		\n\t		vsubpd	      %%zmm14,%%zmm11,%%zmm11	\n\t"\
		"vmovaps	%%zmm0,     (%%rbx)			\n\t		vmovaps	%%zmm8 ,     (%%r11)			\n\t"\
		"vmovaps	%%zmm1,0x040(%%rbx)			\n\t		vmovaps	%%zmm9 ,0x040(%%r11)			\n\t"\
		"vmovaps	%%zmm2,     (%%rdx)			\n\t		vmovaps	%%zmm10,     (%%r13)			\n\t"\
		"vmovaps	%%zmm3,0x040(%%rcx)			\n\t		vmovaps	%%zmm11,0x040(%%r12)			\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm0,%%zmm4		\n\t	vfmadd132pd	%%zmm31,%%zmm8 ,%%zmm12		\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm1,%%zmm5		\n\t	vfmadd132pd	%%zmm31,%%zmm9 ,%%zmm13		\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm2,%%zmm7		\n\t	vfmadd132pd	%%zmm31,%%zmm10,%%zmm15		\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm3,%%zmm6		\n\t	vfmadd132pd	%%zmm31,%%zmm11,%%zmm14		\n\t"\
		"vmovaps	%%zmm4,     (%%rax)			\n\t		vmovaps	%%zmm12,     (%%r10)			\n\t"\
		"vmovaps	%%zmm5,0x040(%%rax)			\n\t		vmovaps	%%zmm13,0x040(%%r10)			\n\t"\
		"vmovaps	%%zmm7,     (%%rcx)			\n\t		vmovaps	%%zmm15,     (%%r12)			\n\t"\
		"vmovaps	%%zmm6,0x040(%%rdx)			\n\t		vmovaps	%%zmm14,0x040(%%r13)			\n\t"\
	/*...Block 2: r5,13,21,29 */						/*...Block 4: r7,15,23,31 */\
		"addq	$0x100,%%rax						\n\t		leaq	0x080(%%rax),%%r10				\n\t"\
		"addq	$0x100,%%rbx						\n\t		leaq	0x080(%%rbx),%%r11				\n\t"\
		"addq	$0x100,%%rcx						\n\t		leaq	0x080(%%rcx),%%r12				\n\t"\
		"addq	$0x100,%%rdx						\n\t		leaq	0x080(%%rdx),%%r13				\n\t"\
		"/* SSE2_RADIX4_DIT_IN_PLACE(): */		\n\t		/* SSE2_RADIX4_DIT_IN_PLACE(): */		\n\t"\
		"vmovaps	     (%%rax),%%zmm2			\n\t		vmovaps	     (%%r10),%%zmm10			\n\t"\
		"vmovaps	0x040(%%rax),%%zmm3			\n\t		vmovaps	0x040(%%r10),%%zmm11			\n\t"\
		"vmovaps	     (%%rbx),%%zmm0			\n\t		vmovaps	     (%%r11),%%zmm8 			\n\t"\
		"vmovaps	0x040(%%rbx),%%zmm1			\n\t		vmovaps	0x040(%%r11),%%zmm9 			\n\t"\
		"vsubpd		%%zmm0,%%zmm2,%%zmm2		\n\t		vsubpd		%%zmm8 ,%%zmm10,%%zmm10		\n\t"\
		"vsubpd		%%zmm1,%%zmm3,%%zmm3		\n\t		vsubpd		%%zmm9 ,%%zmm11,%%zmm11		\n\t"\
		"vmovaps	     (%%rcx),%%zmm6			\n\t		vmovaps	     (%%r12),%%zmm14			\n\t"\
		"vmovaps	0x040(%%rcx),%%zmm7			\n\t		vmovaps	0x040(%%r12),%%zmm15			\n\t"\
		"vmovaps	     (%%rdx),%%zmm4			\n\t		vmovaps	     (%%r13),%%zmm12			\n\t"\
		"vmovaps	0x040(%%rdx),%%zmm5			\n\t		vmovaps	0x040(%%r13),%%zmm13			\n\t"\
		"vsubpd		%%zmm4,%%zmm6,%%zmm6		\n\t		vsubpd		 %%zmm12,%%zmm14,%%zmm14	\n\t"\
		"vsubpd		%%zmm5,%%zmm7,%%zmm7		\n\t		vsubpd	0x040(%%r13),%%zmm15,%%zmm15	\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm2,%%zmm0		\n\t	vfmadd132pd		%%zmm31,%%zmm10,%%zmm8 		\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm3,%%zmm1		\n\t	vfmadd132pd		%%zmm31,%%zmm11,%%zmm9 		\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm6,%%zmm4		\n\t	vfmadd132pd		%%zmm31,%%zmm14,%%zmm12		\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm7,%%zmm5		\n\t	vfmadd132pd		%%zmm31,%%zmm15,%%zmm13		\n\t"\
		"vsubpd	      %%zmm4,%%zmm0,%%zmm0		\n\t		vsubpd	      %%zmm12,%%zmm8 ,%%zmm8 	\n\t"\
		"vsubpd	      %%zmm5,%%zmm1,%%zmm1		\n\t		vsubpd	      %%zmm13,%%zmm9 ,%%zmm9 	\n\t"\
		"vsubpd	      %%zmm7,%%zmm2,%%zmm2		\n\t		vsubpd	      %%zmm15,%%zmm10,%%zmm10	\n\t"\
		"vsubpd	      %%zmm6,%%zmm3,%%zmm3		\n\t		vsubpd	      %%zmm14,%%zmm11,%%zmm11	\n\t"\
		"vmovaps	%%zmm0,     (%%rbx)			\n\t		vmovaps	%%zmm8 ,     (%%r11)			\n\t"\
		"vmovaps	%%zmm1,0x040(%%rbx)			\n\t		vmovaps	%%zmm9 ,0x040(%%r11)			\n\t"\
		"vmovaps	%%zmm2,     (%%rdx)			\n\t		vmovaps	%%zmm10,     (%%r13)			\n\t"\
		"vmovaps	%%zmm3,0x040(%%rcx)			\n\t		vmovaps	%%zmm11,0x040(%%r12)			\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm0,%%zmm4		\n\t	vfmadd132pd	%%zmm31,%%zmm8 ,%%zmm12		\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm1,%%zmm5		\n\t	vfmadd132pd	%%zmm31,%%zmm9 ,%%zmm13		\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm2,%%zmm7		\n\t	vfmadd132pd	%%zmm31,%%zmm10,%%zmm15		\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm3,%%zmm6		\n\t	vfmadd132pd	%%zmm31,%%zmm11,%%zmm14		\n\t"\
		"vmovaps	%%zmm4,     (%%rax)			\n\t		vmovaps	%%zmm12,     (%%r10)			\n\t"\
		"vmovaps	%%zmm5,0x040(%%rax)			\n\t		vmovaps	%%zmm13,0x040(%%r10)			\n\t"\
		"vmovaps	%%zmm7,     (%%rcx)			\n\t		vmovaps	%%zmm15,     (%%r12)			\n\t"\
		"vmovaps	%%zmm6,0x040(%%rdx)			\n\t		vmovaps	%%zmm14,0x040(%%r13)			\n\t"\
		"\n\t"\
	/***************************************************************************************************/\
	/*..and now do four more radix-4 transforms, including the internal and external twiddle factors.  */\
	/*  Write even-index 16-byte output pairs to a(j1), odd-index to a(j2), unpack same as on inputs.  */\
	/*  We do the last 2 radix-4 blocks first, to make the unpack-interleaving of outputs smoother.    */\
	/***************************************************************************************************/\
	/* Main-array addresses still in add0,1; all rax-offsets incr +0x400 in rcol w.r.to lcol: */\
	/*...Block 3: t3,11,19,27 -> r9,13,11,15: */		/*...Block 4: t7,15,23,31 -> r25,29,27,31: */\
		"movq		%[__r9],%%rax				\n\t"\
		"movq		%[__cc0],%%rcx				\n\t"\
	"vmovaps	(%%rcx),%%zmm3	\n\t	vmovaps	0x040(%%rcx),%%zmm11	\n\t"/* cc0,ss0 data shared by both columns */\
		"vmovaps	0x080(%%rax),%%zmm4			\n\t		vmovaps		0x480(%%rax),%%zmm12			\n\t"\
		"vmovaps	0x180(%%rax),%%zmm0			\n\t		vmovaps		0x580(%%rax),%%zmm8 			\n\t"\
		"vmovaps	0x0c0(%%rax),%%zmm5			\n\t		vmovaps		0x4c0(%%rax),%%zmm13			\n\t"\
		"vmovaps	0x1c0(%%rax),%%zmm1			\n\t		vmovaps		0x5c0(%%rax),%%zmm9 			\n\t"\
		"vmulpd		     %%zmm3 ,%%zmm4,%%zmm4	\n\t		vmulpd		     %%zmm11,%%zmm12,%%zmm12	\n\t"\
		"vmulpd		     %%zmm11,%%zmm0,%%zmm0	\n\t		vmulpd		     %%zmm3 ,%%zmm8 ,%%zmm8 	\n\t"\
		"vmulpd		     %%zmm3 ,%%zmm5,%%zmm5	\n\t		vmulpd		     %%zmm11,%%zmm13,%%zmm13	\n\t"\
		"vmulpd		     %%zmm11,%%zmm1,%%zmm1	\n\t		vmulpd		     %%zmm3 ,%%zmm9 ,%%zmm9 	\n\t"\
		"vmovaps	0x080(%%rax),%%zmm6			\n\t		vmovaps		0x480(%%rax),%%zmm14			\n\t"\
		"vmovaps	0x180(%%rax),%%zmm2			\n\t		vmovaps		0x580(%%rax),%%zmm10			\n\t"\
		"vmovaps	0x0c0(%%rax),%%zmm7			\n\t		vmovaps		0x4c0(%%rax),%%zmm15			\n\t"\
	"vfnmadd231pd	     %%zmm11,%%zmm6,%%zmm5	\n\t	vfnmadd231pd	     %%zmm3 ,%%zmm14,%%zmm13	\n\t"\
	"vfnmadd231pd	     %%zmm3 ,%%zmm2,%%zmm1	\n\t	vfnmadd231pd	     %%zmm11,%%zmm10,%%zmm9 	\n\t"\
	" vfmadd231pd	     %%zmm11,%%zmm7,%%zmm4	\n\t	 vfmadd231pd	     %%zmm3 ,%%zmm15,%%zmm12	\n\t"\
	" vfmadd231pd	0x1c0(%%rax),%%zmm3,%%zmm0	\n\t	 vfmadd231pd	0x5c0(%%rax),%%zmm11,%%zmm8 	\n\t"\
		"vmovaps	%%zmm5,%%zmm7				\n\t		vmovaps		%%zmm13,%%zmm15				\n\t"\
		"vmovaps	%%zmm4,%%zmm6				\n\t		vmovaps		%%zmm12,%%zmm14				\n\t"\
		"vaddpd		%%zmm0,%%zmm4,%%zmm4		\n\t		vaddpd		%%zmm8 ,%%zmm12,%%zmm12		\n\t"\
		"vaddpd		%%zmm1,%%zmm5,%%zmm5		\n\t		vaddpd		%%zmm9 ,%%zmm13,%%zmm13		\n\t"\
		"vsubpd		%%zmm0,%%zmm6,%%zmm6		\n\t		vsubpd		%%zmm8 ,%%zmm14,%%zmm14		\n\t"\
		"vsubpd		%%zmm1,%%zmm7,%%zmm7		\n\t		vsubpd		%%zmm9 ,%%zmm15,%%zmm15		\n\t"\
		"vmovaps	0x100(%%rax),%%zmm2			\n\t		vmovaps		0x500(%%rax),%%zmm10		\n\t"\
		"vmovaps	0x140(%%rax),%%zmm3			\n\t		vmovaps		0x540(%%rax),%%zmm11		\n\t"\
		"vmovaps	     (%%rax),%%zmm0			\n\t		vmovaps		0x400(%%rax),%%zmm8 		\n\t"\
		"vmovaps	0x040(%%rax),%%zmm1			\n\t		vmovaps		0x440(%%rax),%%zmm9 		\n\t"\
		"vaddpd		0x140(%%rax),%%zmm2,%%zmm2	\n\t		vsubpd		0x540(%%rax),%%zmm10,%%zmm10\n\t"\
		"vsubpd		0x100(%%rax),%%zmm3,%%zmm3	\n\t		vaddpd		0x500(%%rax),%%zmm11,%%zmm11\n\t"\
		"vmulpd		     %%zmm29,%%zmm2,%%zmm2	\n\t		vmulpd		%%zmm29,%%zmm10,%%zmm10		\n\t"\
		"vmulpd		     %%zmm29,%%zmm3,%%zmm3	\n\t		vmulpd		%%zmm29,%%zmm11,%%zmm11		\n\t"\
		"vsubpd		%%zmm2,%%zmm0,%%zmm0		\n\t		vsubpd		%%zmm10,%%zmm8 ,%%zmm8 		\n\t"\
		"vsubpd		%%zmm3,%%zmm1,%%zmm1		\n\t		vsubpd		%%zmm11,%%zmm9 ,%%zmm9 		\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm0,%%zmm2		\n\t	vfmadd132pd	%%zmm31,%%zmm8 ,%%zmm10			\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm1,%%zmm3		\n\t	vfmadd132pd	%%zmm31,%%zmm9 ,%%zmm11			\n\t"\
	"addq	$0x480,%%rcx		\n\t"/* c1 from cc0 */\
	"leaq	0x80(%%rcx),%%rdx	\n\t"/* c9 from c1  */\
	/* Using the upper block(s) of the main array for temp-storage in the section below led to a nasty AVX bug to track down: */\
	/* [***The following details are the AVX-512-ification of the the analogous comment in the AVX version of this macro ***]
	In fermat-mod mode the 8 block addresses in ascending order are add0-7 with no 'gaps' between blocks, whereas for mersenne-mod
	the addresses in asc. order are add0,2,4,6,7,5,3,1 with a gap between contiguous-data-block quartets 0,2,4,6 and 7,5,3,1. Thus
	for fermat-mod we need [add4] as the base-address of the 'high-half' block for temp-storage; for mersenne-mod we need [add7]
	In both cases we have that add4 < add7 so instead use (add3 - add1): > 0 for fermat-mod, < 0 for mersenne - to differentiate: */\
		"movq	%[__add4],%%rdi		\n\t"/* destroyable copy of add4 */\
		"movq	%[__add4],%%rbx		\n\t"\
		"subq	%[__add3],%%rdi		\n\t"/* rdi = (add4 - add3); if this yields a borrow (i.e. sets CF) it's mersenne, else fermat. */\
	/*	"cmovcq %[__add7],%%rbx	\n\t" // CMOV not supported on k1om, so emulate using jump-if-CF-not-set: */\
	"jnc skip	\n\t"\
		"movq	%[__add7],%%rbx	\n\t" /* if CF set (i.e. h > l), copy source [add7] into dest (rbx), else leave dest = [add4]. */\
	"skip:	\n\t"\
		/* rbx shared between rcol/lcol; rcx/rdx-offsets incr +0x100 in rcol for rest of block: */\
		"vsubpd		%%zmm4,%%zmm2,%%zmm2			\n\t		vsubpd		%%zmm14,%%zmm8 ,%%zmm8 		\n\t"\
		"vsubpd		%%zmm5,%%zmm3,%%zmm3			\n\t		vsubpd		%%zmm15,%%zmm9 ,%%zmm9 		\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm2,%%zmm4			\n\t	vfmadd132pd	%%zmm31,%%zmm8 ,%%zmm14			\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm3,%%zmm5			\n\t	vfmadd132pd	%%zmm31,%%zmm9 ,%%zmm15			\n\t"\
		"vmovaps	%%zmm2,     (%%rax)				\n\t		vmovaps		%%zmm8 ,0x400(%%rax)		\n\t"\
		"vmovaps	%%zmm3,0x040(%%rax)				\n\t		vmovaps		%%zmm9 ,0x440(%%rax)		\n\t"\
		"vmovaps	%%zmm4,%%zmm2					\n\t		vmovaps		%%zmm14,%%zmm8 		\n\t"\
		"vmovaps	%%zmm5,%%zmm3					\n\t		vmovaps		%%zmm15,%%zmm9 		\n\t"\
		"vmulpd		     (%%rcx),%%zmm4,%%zmm4		\n\t		vmulpd		0x200(%%rcx),%%zmm14,%%zmm14\n\t"\
		"vmulpd		     (%%rcx),%%zmm5,%%zmm5		\n\t		vmulpd		0x200(%%rcx),%%zmm15,%%zmm15\n\t"\
		"vsubpd		%%zmm7,%%zmm0,%%zmm0			\n\t		vsubpd		%%zmm13,%%zmm10,%%zmm10		\n\t"\
		"vsubpd		%%zmm6,%%zmm1,%%zmm1			\n\t		vsubpd		%%zmm12,%%zmm11,%%zmm11		\n\t"\
	"vfnmadd231pd	0x040(%%rcx),%%zmm2,%%zmm5		\n\t	vfnmadd231pd	0x240(%%rcx),%%zmm8 ,%%zmm15	\n\t"\
	" vfmadd231pd	0x040(%%rcx),%%zmm3,%%zmm4		\n\t	 vfmadd231pd	0x240(%%rcx),%%zmm9 ,%%zmm14	\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm0,%%zmm7			\n\t	vfmadd132pd	%%zmm31,%%zmm10,%%zmm13			\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm1,%%zmm6			\n\t	vfmadd132pd	%%zmm31,%%zmm11,%%zmm12			\n\t"\
		"vmovaps	%%zmm5,0x040(%%rbx)				\n\t		vmovaps		%%zmm15,0x0c0(%%rbx)		\n\t"\
		"vmovaps	%%zmm4,     (%%rbx)				\n\t		vmovaps		%%zmm14,0x080(%%rbx)		\n\t"\
		"vmovaps		 (%%rax),%%zmm4				\n\t		vmovaps		0x400(%%rax),%%zmm14		\n\t"\
		"vmovaps	0x040(%%rax),%%zmm5				\n\t		vmovaps		0x440(%%rax),%%zmm15		\n\t"\
		"vmovaps	%%zmm4,%%zmm2					\n\t		vmovaps		%%zmm14,%%zmm8 		\n\t"\
		"vmovaps	%%zmm5,%%zmm3					\n\t		vmovaps		%%zmm15,%%zmm9 		\n\t"\
		"vmulpd		     (%%rdx),%%zmm4,%%zmm4		\n\t		vmulpd		0x200(%%rdx),%%zmm14,%%zmm14\n\t"\
		"vmulpd		     (%%rdx),%%zmm5,%%zmm5		\n\t		vmulpd		0x200(%%rdx),%%zmm15,%%zmm15\n\t"\
	"vfnmadd231pd	0x040(%%rdx),%%zmm2,%%zmm5		\n\t	vfnmadd231pd	0x240(%%rdx),%%zmm8 ,%%zmm15	\n\t"\
	" vfmadd231pd	0x040(%%rdx),%%zmm3,%%zmm4		\n\t	 vfmadd231pd	0x240(%%rdx),%%zmm9 ,%%zmm14	\n\t"\
		"vmovaps	%%zmm5,0x240(%%rbx)				\n\t		vmovaps		%%zmm15,0x2c0(%%rbx)		\n\t"\
		"vmovaps	%%zmm4,0x200(%%rbx)				\n\t		vmovaps		%%zmm14,0x280(%%rbx)		\n\t"\
		"addq	$0x100,%%rcx	\n\t"/* c1  from c1; rcol has c7 */\
		"addq	$0x100,%%rdx	\n\t"/* c9  from c9; rcol has c15*/\
		"vmovaps	%%zmm7,%%zmm4					\n\t		vmovaps		%%zmm13,%%zmm8 		\n\t"\
		"vmovaps	%%zmm1,%%zmm5					\n\t		vmovaps		%%zmm11,%%zmm9 		\n\t"\
		"vmulpd		     (%%rcx),%%zmm7,%%zmm7		\n\t		vmulpd		0x200(%%rcx),%%zmm13,%%zmm13\n\t"\
		"vmulpd		     (%%rcx),%%zmm1,%%zmm1		\n\t		vmulpd		0x200(%%rcx),%%zmm11,%%zmm11\n\t"\
	"vfnmadd231pd	0x040(%%rcx),%%zmm4,%%zmm1		\n\t	vfnmadd231pd	0x240(%%rcx),%%zmm8 ,%%zmm11	\n\t"\
	" vfmadd231pd	0x040(%%rcx),%%zmm5,%%zmm7		\n\t	 vfmadd231pd	0x240(%%rcx),%%zmm9 ,%%zmm13	\n\t"\
		"vmovaps	%%zmm0,%%zmm4					\n\t		vmovaps		%%zmm10,%%zmm8 		\n\t"\
		"vmovaps	%%zmm6,%%zmm5					\n\t		vmovaps		%%zmm12,%%zmm9 		\n\t"\
		"vmulpd		     (%%rdx),%%zmm0,%%zmm0		\n\t		vmulpd		0x200(%%rdx),%%zmm10,%%zmm10\n\t"\
		"vmulpd		     (%%rdx),%%zmm6,%%zmm6		\n\t		vmulpd		0x200(%%rdx),%%zmm12,%%zmm12\n\t"\
		"vmovaps	%%zmm1,0x140(%%rbx)				\n\t		vmovaps		%%zmm11,0x1c0(%%rbx)		\n\t"\
		"vmovaps	%%zmm7,0x100(%%rbx)				\n\t		vmovaps		%%zmm13,0x180(%%rbx)		\n\t"\
	"vfnmadd231pd	0x040(%%rdx),%%zmm4,%%zmm6		\n\t	vfnmadd231pd	0x240(%%rdx),%%zmm8 ,%%zmm12	\n\t"\
	" vfmadd231pd	0x040(%%rdx),%%zmm5,%%zmm0		\n\t	 vfmadd231pd	0x240(%%rdx),%%zmm9 ,%%zmm10	\n\t"\
		"vmovaps	%%zmm6,0x340(%%rbx)				\n\t		vmovaps		%%zmm12,0x3c0(%%rbx)		\n\t"\
		"vmovaps	%%zmm0,0x300(%%rbx)				\n\t		vmovaps		%%zmm10,0x380(%%rbx)		\n\t"\
		/*...Block 1: t1,9,17,25 -> r1,5,3,7: */			/*...Block 3: t5,13,21,29 -> r17,21,19,23: */\
		"movq	%[__r1],%%rax	\n\t"/* r17 in rcol */\
		"vmovaps	     (%%rax),%%zmm0				\n\t		vmovaps		0x480(%%rax),%%zmm12	\n\t"\
		"vmovaps	0x040(%%rax),%%zmm1				\n\t		vmovaps		0x4c0(%%rax),%%zmm13		\n\t"\
		"vmovaps	0x100(%%rax),%%zmm2				\n\t		vmovaps		0x580(%%rax),%%zmm8 		\n\t"\
		"vmovaps	0x140(%%rax),%%zmm3				\n\t		vmovaps		0x5c0(%%rax),%%zmm9 		\n\t"\
		"vsubpd		%%zmm2,%%zmm0,%%zmm0			\n\t		vaddpd		0x4c0(%%rax),%%zmm12,%%zmm12\n\t"\
		"vsubpd		%%zmm3,%%zmm1,%%zmm1			\n\t		vsubpd		0x480(%%rax),%%zmm13,%%zmm13\n\t"\
		"vmovaps	0x080(%%rax),%%zmm4				\n\t		vsubpd		0x5c0(%%rax),%%zmm8 ,%%zmm8 \n\t"\
		"vmovaps	0x0c0(%%rax),%%zmm5				\n\t		vaddpd		0x580(%%rax),%%zmm9 ,%%zmm9 \n\t"\
		"vmovaps	0x180(%%rax),%%zmm6				\n\t		vmulpd		%%zmm29,%%zmm12,%%zmm12		\n\t"\
		"vmovaps	0x1c0(%%rax),%%zmm7				\n\t		vmulpd		%%zmm29,%%zmm13,%%zmm13		\n\t"\
		"vsubpd		%%zmm6,%%zmm4,%%zmm4			\n\t		vmovaps		%%zmm12,%%zmm14		\n\t"\
		"vsubpd		%%zmm7,%%zmm5,%%zmm5			\n\t		vmovaps		%%zmm13,%%zmm15		\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm0,%%zmm2			\n\t	vfnmadd231pd	%%zmm29,%%zmm8 ,%%zmm12		\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm1,%%zmm3			\n\t	vfnmadd231pd	%%zmm29,%%zmm9 ,%%zmm13		\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm4,%%zmm6			\n\t	 vfmadd231pd	%%zmm29,%%zmm8 ,%%zmm14		\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm5,%%zmm7			\n\t	 vfmadd231pd	%%zmm29,%%zmm9 ,%%zmm15		\n\t"\
		"movq	%%rbx,%%rcx	\n\t"/* rbx saves 'high-half' block address; move that to rcx and put 'low-half' block address into rbx */\
		"movq	%[__add0],%%rbx						\n\t		vmovaps		0x400(%%rax),%%zmm8 		\n\t"\
		"subq	$0x500,%%rdx	/* c8 from c13 */	\n\t		vmovaps		0x440(%%rax),%%zmm9 		\n\t"\
		"vaddpd		%%zmm6,%%zmm2,%%zmm2			\n\t		vmovaps		0x500(%%rax),%%zmm10		\n\t"\
		"vaddpd		%%zmm7,%%zmm3,%%zmm3			\n\t		vmovaps		0x540(%%rax),%%zmm11		\n\t"\
		"vmovaps	%%zmm2,     (%%rbx)				\n\t		vsubpd		0x540(%%rax),%%zmm8 ,%%zmm8 \n\t"\
		"vmovaps	%%zmm3,0x040(%%rbx)				\n\t		vsubpd		0x500(%%rax),%%zmm9 ,%%zmm9 \n\t"\
	"vfnmadd231pd	%%zmm31,%%zmm6,%%zmm2			\n\t		vaddpd		0x400(%%rax),%%zmm11,%%zmm11\n\t"\
	"vfnmadd231pd	%%zmm31,%%zmm7,%%zmm3			\n\t		vaddpd		0x440(%%rax),%%zmm10,%%zmm10\n\t"\
		"vmovaps	%%zmm2,%%zmm6					\n\t		vsubpd		%%zmm12,%%zmm11,%%zmm11		\n\t"\
		"vmovaps	%%zmm3,%%zmm7					\n\t		vsubpd		%%zmm13,%%zmm9 ,%%zmm9 		\n\t"\
		"vmulpd		     (%%rdx),%%zmm2,%%zmm2		\n\t	vfmadd132pd	%%zmm31,%%zmm11,%%zmm12			\n\t"\
		"vmulpd		     (%%rdx),%%zmm3,%%zmm3		\n\t	vfmadd132pd	%%zmm31,%%zmm9 ,%%zmm13			\n\t"\
	"vfnmadd231pd	0x040(%%rdx),%%zmm6,%%zmm3		\n\t		vmovaps		%%zmm11,     (%%rax)		\n\t"\
	" vfmadd231pd	0x040(%%rdx),%%zmm7,%%zmm2		\n\t		vmovaps		%%zmm9 ,0x040(%%rax)		\n\t"\
		"vmovaps	0x240(%%rcx),%%zmm7				\n\t		vmovaps		%%zmm12,%%zmm11		\n\t"\
		"vmovaps	0x200(%%rcx),%%zmm6				\n\t		vmovaps		%%zmm13,%%zmm9 		\n\t"\
		"vmovaps	%%zmm3,0x240(%%rax)	/* r9 */	\n\t		vmulpd		0x180(%%rdx),%%zmm12,%%zmm12	\n\t"/* c2 */\
		"vmovaps	%%zmm2,0x200(%%rax)	/* r8 */	\n\t		vmulpd		0x180(%%rdx),%%zmm13,%%zmm13\n\t"\
		"vmovaps	%%zmm7,0x640(%%rax)	/* r25 */	\n\t	vfnmadd231pd	0x1c0(%%rdx),%%zmm11,%%zmm13\n\t"\
		"vmovaps	%%zmm6,0x600(%%rax)	/* r24 */	\n\t	 vfmadd231pd	0x1c0(%%rdx),%%zmm9 ,%%zmm12\n\t"\
		"vmovaps	0x040(%%rbx),%%zmm3				\n\t		vmovaps	0x0c0(%%rcx),%%zmm11	\n\t"\
		"vmovaps	     (%%rbx),%%zmm2				\n\t		vmovaps	0x080(%%rcx),%%zmm9		\n\t"\
		"vmovaps	0x040(%%rcx),%%zmm7				\n\t		vmovaps	%%zmm13,0x0c0(%%rax)	/* r3 */\n\t"\
		"vmovaps	     (%%rcx),%%zmm6				\n\t		vmovaps	%%zmm12,0x080(%%rax)	/* r2 */\n\t"\
	/* swap rcol/lcol in next 4 lines: */\
		"vmovaps		     (%%rax),%%zmm12		\n\t		vmovaps	%%zmm2 ,     (%%rax)	/* r0 */\n\t"\
		"vmovaps		0x040(%%rax),%%zmm13		\n\t		vmovaps	%%zmm3 ,0x040(%%rax)	/* r1 */\n\t"\
		"vmovaps	%%zmm11,0x4c0(%%rax)	/* r19*/\n\t		vmovaps		%%zmm12,%%zmm11		\n\t"\
		"vmovaps	%%zmm9 ,0x480(%%rax)	/* r18*/\n\t		vmovaps		%%zmm13,%%zmm9 		\n\t"\
		"														addq	$0x080,%%rdx		/* c4 in lcol; c10 in rcol*/		\n\t"\
		"vmovaps	%%zmm7,0x440(%%rax)	/* r17 */	\n\t		vmulpd		0x180(%%rdx),%%zmm12,%%zmm12		\n\t"\
		"vmovaps	%%zmm6,0x400(%%rax)	/* r16 */	\n\t		vmulpd		0x180(%%rdx),%%zmm13,%%zmm13		\n\t"\
		"vaddpd		%%zmm5,%%zmm0,%%zmm0			\n\t	vfnmadd231pd	0x1c0(%%rdx),%%zmm11,%%zmm13		\n\t"\
		"vsubpd		%%zmm4,%%zmm1,%%zmm1			\n\t	 vfmadd231pd	0x1c0(%%rdx),%%zmm9 ,%%zmm12		\n\t"\
		"vmovaps	%%zmm0,%%zmm2					\n\t		vmovaps	0x2c0(%%rcx),%%zmm11	\n\t"\
		"vmovaps	%%zmm1,%%zmm3					\n\t		vmovaps	0x280(%%rcx),%%zmm9		\n\t"\
		"			\n\t		vmovaps	%%zmm13,0x2c0(%%rax)	/* r11 */\n\t"\
		"			\n\t		vmovaps	%%zmm12,0x280(%%rax)	/* r10 */\n\t"\
		"vmovaps	%%zmm0,%%zmm6					\n\t		vmovaps	%%zmm11,0x6c0(%%rax)	/* r27 */\n\t"\
		"vmovaps	%%zmm1,%%zmm7					\n\t		vmovaps	%%zmm9 ,0x680(%%rax)	/* r26 */\n\t"\
	"vfnmadd231pd	%%zmm31,%%zmm5,%%zmm0			\n\t		vsubpd		%%zmm15,%%zmm8 ,%%zmm8 		\n\t"\
	" vfmadd231pd	%%zmm31,%%zmm4,%%zmm1			\n\t		vsubpd		%%zmm14,%%zmm10,%%zmm10		\n\t"\
		"vmulpd		     (%%rdx),%%zmm2,%%zmm2		\n\t	vfmadd132pd	%%zmm31,%%zmm8 ,%%zmm15		\n\t"\
		"vmulpd		     (%%rdx),%%zmm3,%%zmm3		\n\t	vfmadd132pd	%%zmm31,%%zmm10,%%zmm14		\n\t"\
	"vfnmadd231pd	0x040(%%rdx),%%zmm6,%%zmm3		\n\t		vmovaps		%%zmm15,%%zmm12		\n\t"\
	" vfmadd231pd	0x040(%%rdx),%%zmm7,%%zmm2		\n\t		vmovaps		%%zmm10,%%zmm13		\n\t"\
														"	addq	$0x080,%%rdx	\n\t"/* c12 in lcol; c6 in rcol*/\
		"vmovaps	0x140(%%rcx),%%zmm7				\n\t		vmulpd		0x180(%%rdx),%%zmm15,%%zmm15	\n\t"\
		"vmovaps	0x100(%%rcx),%%zmm6				\n\t		vmulpd		0x180(%%rdx),%%zmm10,%%zmm10	\n\t"\
		"vmovaps	%%zmm3,0x140(%%rax)	/* r5 */	\n\t	vfnmadd231pd	0x1c0(%%rdx),%%zmm12,%%zmm10	\n\t"\
		"vmovaps	%%zmm2,0x100(%%rax)	/* r4 */	\n\t	 vfmadd231pd	0x1c0(%%rdx),%%zmm13,%%zmm15	\n\t"\
		"vmovaps	%%zmm7,0x540(%%rax)	/* r21 */	\n\t		vmovaps	0x1c0(%%rcx),%%zmm13	\n\t"\
		"vmovaps	%%zmm6,0x500(%%rax)	/* r20 */	\n\t		vmovaps	0x180(%%rcx),%%zmm12	\n\t"\
		"vmovaps	%%zmm0,%%zmm6					\n\t		vmovaps	%%zmm10,0x1c0(%%rax)	/* r7 */\n\t"\
		"vmovaps	%%zmm1,%%zmm7					\n\t		vmovaps	%%zmm15,0x180(%%rax)	/* r6 */\n\t"\
		"vmulpd		     (%%rdx),%%zmm0,%%zmm0		\n\t		vmovaps	%%zmm13,0x5c0(%%rax)	/* r23 */\n\t"\
		"vmulpd		     (%%rdx),%%zmm1,%%zmm1		\n\t		vmovaps	%%zmm12,0x580(%%rax)	/* r22 */\n\t"\
	"vfnmadd231pd	0x040(%%rdx),%%zmm6,%%zmm1		\n\t		vmovaps		%%zmm8 ,%%zmm12		\n\t"\
	" vfmadd231pd	0x040(%%rdx),%%zmm7,%%zmm0		\n\t		vmovaps		%%zmm14,%%zmm13		\n\t"\
		"														vmulpd		0x200(%%rdx),%%zmm8 ,%%zmm8 \n\t"/* c14 */\
		"														vmulpd		0x200(%%rdx),%%zmm14,%%zmm14		\n\t"\
		"													vfnmadd231pd	0x240(%%rdx),%%zmm12,%%zmm14		\n\t"\
		"													 vfmadd231pd	0x240(%%rdx),%%zmm13,%%zmm8 		\n\t"\
		"vmovaps	0x340(%%rcx),%%zmm7				\n\t		vmovaps	0x3c0(%%rcx),%%zmm13	\n\t"\
		"vmovaps	0x300(%%rcx),%%zmm6				\n\t		vmovaps	0x380(%%rcx),%%zmm12	\n\t"\
		"vmovaps	%%zmm1,0x340(%%rax)	/* r13 */	\n\t		vmovaps	%%zmm14,0x3c0(%%rax)	/* r15 */\n\t"\
		"vmovaps	%%zmm0,0x300(%%rax)	/* r12 */	\n\t		vmovaps	%%zmm8 ,0x380(%%rax)	/* r14 */\n\t"\
		"vmovaps	%%zmm7,0x740(%%rax)	/* r29 */	\n\t		vmovaps	%%zmm13,0x7c0(%%rax)	/* r31 */\n\t"\
		"vmovaps	%%zmm6,0x700(%%rax)	/* r28 */	\n\t		vmovaps	%%zmm12,0x780(%%rax)	/* r30 */\n\t"\
	/*******************************************
	***** Finish with 8-way 'un'terleaving: ****
	Using the AVX-512 data layout, the rcol pattern is:
		a[ 0- 7] = re[ 0, 8, 1, 9, 2,10, 3,11].d0	a[ 8-15] = im[ 0, 8, 1, 9, 2,10, 3,11].d0
		a[16-23] = re[ 4,12, 5,13, 6,14, 7,15].d0	a[24-31] = im[ 4,12, 5,13, 6,14, 7,15].d0 ,
	and remaining seven 32-double blocks repeat same pattern with elts d1-d7 of the vector-doubles.
	*******************************************/\
		"movq	%[__add0],%%rax	\n\t"\
		"movq	%[__add1],%%rbx	\n\t"\
		"movq	%[__add2],%%rcx	\n\t"\
		"movq	%[__add3],%%rdx	\n\t"\
		"movq	%[__add4],%%r10	\n\t"\
		"movq	%[__add5],%%r11	\n\t"\
		"movq	%[__add6],%%r12	\n\t"\
		"movq	%[__add7],%%r13	\n\t"\
	"movq	%[__r1] ,%%rsi	\n\t"\
	/**** a[ 0- 7] = re[ 0, 8, 1, 9, 2,10, 3,11].d0	a[ 8-15] = im[ 0, 8, 1, 9, 2,10, 3,11].d0 : ****/\
		"vmovaps 	     (%%rsi),%%zmm0					\n\t		vmovaps 0x040(%%rsi),%%zmm10	\n\t"\
		"vmovaps 	0x400(%%rsi),%%zmm1					\n\t		vmovaps 0x440(%%rsi),%%zmm11	\n\t"\
		"vmovaps 	0x080(%%rsi),%%zmm2					\n\t		vmovaps 0x0c0(%%rsi),%%zmm12	\n\t"\
		"vmovaps 	0x480(%%rsi),%%zmm3					\n\t		vmovaps 0x4c0(%%rsi),%%zmm13	\n\t"\
		"vmovaps 	0x100(%%rsi),%%zmm4					\n\t		vmovaps 0x140(%%rsi),%%zmm14	\n\t"\
		"vmovaps 	0x500(%%rsi),%%zmm5					\n\t		vmovaps 0x540(%%rsi),%%zmm15	\n\t"\
		"vmovaps 	0x180(%%rsi),%%zmm6					\n\t		vmovaps 0x1c0(%%rsi),%%zmm16	\n\t"\
		"vmovaps 	0x580(%%rsi),%%zmm7					\n\t		vmovaps 0x5c0(%%rsi),%%zmm17	\n\t"\
		"\n\t"\
		"vblendmpd		%%zmm1%{cdab%},%%zmm0,%%zmm8%{%%k1%}	\n\t	vblendmpd		%%zmm11%{cdab%},%%zmm10,%%zmm9 %{%%k1%}	\n\t"\
		"vblendmpd		%%zmm0%{cdab%},%%zmm1,%%zmm1%{%%k2%}	\n\t	vblendmpd		%%zmm10%{cdab%},%%zmm11,%%zmm11%{%%k2%}	\n\t"\
		"vblendmpd		%%zmm3%{cdab%},%%zmm2,%%zmm0%{%%k1%}	\n\t	vblendmpd		%%zmm13%{cdab%},%%zmm12,%%zmm10%{%%k1%}	\n\t"\
		"vblendmpd		%%zmm2%{cdab%},%%zmm3,%%zmm3%{%%k2%}	\n\t	vblendmpd		%%zmm12%{cdab%},%%zmm13,%%zmm13%{%%k2%}	\n\t"\
		"vblendmpd		%%zmm5%{cdab%},%%zmm4,%%zmm2%{%%k1%}	\n\t	vblendmpd		%%zmm15%{cdab%},%%zmm14,%%zmm12%{%%k1%}	\n\t"\
		"vblendmpd		%%zmm4%{cdab%},%%zmm5,%%zmm5%{%%k2%}	\n\t	vblendmpd		%%zmm14%{cdab%},%%zmm15,%%zmm15%{%%k2%}	\n\t"\
		"vblendmpd		%%zmm7%{cdab%},%%zmm6,%%zmm4%{%%k1%}	\n\t	vblendmpd		%%zmm17%{cdab%},%%zmm16,%%zmm14%{%%k1%}	\n\t"\
		"vblendmpd		%%zmm6%{cdab%},%%zmm7,%%zmm7%{%%k2%}	\n\t	vblendmpd		%%zmm16%{cdab%},%%zmm17,%%zmm17%{%%k2%}	\n\t"\
		"\n\t"\
		"vblendmpd		%%zmm0%{badc%},%%zmm8,%%zmm6%{%%k3%}	\n\t	vblendmpd		%%zmm10%{badc%},%%zmm9 ,%%zmm16%{%%k3%}	\n\t"\
		"vblendmpd		%%zmm8%{badc%},%%zmm0,%%zmm0%{%%k4%}	\n\t	vblendmpd		%%zmm9 %{badc%},%%zmm10,%%zmm10%{%%k4%}	\n\t"\
		"vblendmpd		%%zmm3%{badc%},%%zmm1,%%zmm8%{%%k3%}	\n\t	vblendmpd		%%zmm13%{badc%},%%zmm11,%%zmm9 %{%%k3%}	\n\t"\
		"vblendmpd		%%zmm1%{badc%},%%zmm3,%%zmm3%{%%k4%}	\n\t	vblendmpd		%%zmm11%{badc%},%%zmm13,%%zmm13%{%%k4%}	\n\t"\
		"vblendmpd		%%zmm4%{badc%},%%zmm2,%%zmm1%{%%k3%}	\n\t	vblendmpd		%%zmm14%{badc%},%%zmm12,%%zmm11%{%%k3%}	\n\t"\
		"vblendmpd		%%zmm2%{badc%},%%zmm4,%%zmm4%{%%k4%}	\n\t	vblendmpd		%%zmm12%{badc%},%%zmm14,%%zmm14%{%%k4%}	\n\t"\
		"vblendmpd		%%zmm7%{badc%},%%zmm5,%%zmm2%{%%k3%}	\n\t	vblendmpd		%%zmm17%{badc%},%%zmm15,%%zmm12%{%%k3%}	\n\t"\
		"vblendmpd		%%zmm5%{badc%},%%zmm7,%%zmm7%{%%k4%}	\n\t	vblendmpd		%%zmm15%{badc%},%%zmm17,%%zmm17%{%%k4%}	\n\t"\
		"\n\t"\
		"vpermf32x4	$78,%%zmm1,%%zmm1							\n\t	vpermf32x4	$78,%%zmm11,%%zmm11							\n\t"\
		"vpermf32x4	$78,%%zmm2,%%zmm2							\n\t	vpermf32x4	$78,%%zmm12,%%zmm12							\n\t"\
		"vpermf32x4	$78,%%zmm4,%%zmm4							\n\t	vpermf32x4	$78,%%zmm14,%%zmm14							\n\t"\
		"vpermf32x4	$78,%%zmm7,%%zmm7							\n\t	vpermf32x4	$78,%%zmm17,%%zmm17							\n\t"\
		"vblendmpd		%%zmm1,%%zmm6,%%zmm5%{%%k5%}			\n\t	vblendmpd		%%zmm11,%%zmm16,%%zmm15%{%%k5%}			\n\t"\
		"vblendmpd		%%zmm6,%%zmm1,%%zmm1%{%%k5%}			\n\t	vblendmpd		%%zmm16,%%zmm11,%%zmm11%{%%k5%}			\n\t"\
		"vblendmpd		%%zmm2,%%zmm8,%%zmm6%{%%k5%}			\n\t	vblendmpd		%%zmm12,%%zmm9 ,%%zmm16%{%%k5%}			\n\t"\
		"vblendmpd		%%zmm8,%%zmm2,%%zmm2%{%%k5%}			\n\t	vblendmpd		%%zmm9 ,%%zmm12,%%zmm12%{%%k5%}			\n\t"\
		"vblendmpd		%%zmm4,%%zmm0,%%zmm8%{%%k5%}			\n\t	vblendmpd		%%zmm14,%%zmm10,%%zmm9 %{%%k5%}			\n\t"\
		"vblendmpd		%%zmm0,%%zmm4,%%zmm4%{%%k5%}			\n\t	vblendmpd		%%zmm10,%%zmm14,%%zmm14%{%%k5%}			\n\t"\
		"vblendmpd		%%zmm7,%%zmm3,%%zmm0%{%%k5%}			\n\t	vblendmpd		%%zmm17,%%zmm13,%%zmm10%{%%k5%}			\n\t"\
		"vblendmpd		%%zmm3,%%zmm7,%%zmm7%{%%k5%}			\n\t	vblendmpd		%%zmm13,%%zmm17,%%zmm17%{%%k5%}			\n\t"\
		"vpermf32x4	$78,%%zmm1,%%zmm1							\n\t	vpermf32x4	$78,%%zmm11,%%zmm11							\n\t"\
		"vpermf32x4	$78,%%zmm2,%%zmm2							\n\t	vpermf32x4	$78,%%zmm12,%%zmm12							\n\t"\
		"vpermf32x4	$78,%%zmm4,%%zmm4							\n\t	vpermf32x4	$78,%%zmm14,%%zmm14							\n\t"\
		"vpermf32x4	$78,%%zmm7,%%zmm7							\n\t	vpermf32x4	$78,%%zmm17,%%zmm17							\n\t"\
		"\n\t"\
		"vmovaps		%%zmm5,(%%rax)			\n\t	vmovaps	%%zmm15,0x40(%%rax)	\n\t"\
		"vmovaps		%%zmm6,(%%rbx)			\n\t	vmovaps	%%zmm16,0x40(%%rbx)	\n\t"\
		"vmovaps		%%zmm8,(%%rcx)			\n\t	vmovaps	%%zmm9 ,0x40(%%rcx)	\n\t"\
		"vmovaps		%%zmm0,(%%rdx)			\n\t	vmovaps	%%zmm10,0x40(%%rdx)	\n\t"\
		"vmovaps		%%zmm1,(%%r10)			\n\t	vmovaps	%%zmm11,0x40(%%r10)	\n\t"\
		"vmovaps		%%zmm2,(%%r11)			\n\t	vmovaps	%%zmm12,0x40(%%r11)	\n\t"\
		"vmovaps		%%zmm4,(%%r12)			\n\t	vmovaps	%%zmm14,0x40(%%r12)	\n\t"\
		"vmovaps		%%zmm7,(%%r13)			\n\t	vmovaps	%%zmm17,0x40(%%r13)	\n\t"\
		"\n\t"\
	/**** a[16-23] = re[ 4,12, 5,13, 6,14, 7,15].d0	a[24-31] = im[ 4,12, 5,13, 6,14, 7,15].d0 : ****/\
	"addq	$0x80,%%rax	\n\t"\
	"addq	$0x80,%%rbx	\n\t"\
	"addq	$0x80,%%rcx	\n\t"\
	"addq	$0x80,%%rdx	\n\t"\
	"addq	$0x80,%%r10	\n\t"\
	"addq	$0x80,%%r11	\n\t"\
	"addq	$0x80,%%r12	\n\t"\
	"addq	$0x80,%%r13	\n\t"\
	"addq	$0x200,%%rsi	\n\t"\
		"vmovaps 	     (%%rsi),%%zmm0					\n\t		vmovaps 0x040(%%rsi),%%zmm10	\n\t"\
		"vmovaps 	0x400(%%rsi),%%zmm1					\n\t		vmovaps 0x440(%%rsi),%%zmm11	\n\t"\
		"vmovaps 	0x080(%%rsi),%%zmm2					\n\t		vmovaps 0x0c0(%%rsi),%%zmm12	\n\t"\
		"vmovaps 	0x480(%%rsi),%%zmm3					\n\t		vmovaps 0x4c0(%%rsi),%%zmm13	\n\t"\
		"vmovaps 	0x100(%%rsi),%%zmm4					\n\t		vmovaps 0x140(%%rsi),%%zmm14	\n\t"\
		"vmovaps 	0x500(%%rsi),%%zmm5					\n\t		vmovaps 0x540(%%rsi),%%zmm15	\n\t"\
		"vmovaps 	0x180(%%rsi),%%zmm6					\n\t		vmovaps 0x1c0(%%rsi),%%zmm16	\n\t"\
		"vmovaps 	0x580(%%rsi),%%zmm7					\n\t		vmovaps 0x5c0(%%rsi),%%zmm17	\n\t"\
		"\n\t"\
		"vblendmpd		%%zmm1%{cdab%},%%zmm0,%%zmm8%{%%k1%}	\n\t	vblendmpd		%%zmm11%{cdab%},%%zmm10,%%zmm9 %{%%k1%}	\n\t"\
		"vblendmpd		%%zmm0%{cdab%},%%zmm1,%%zmm1%{%%k2%}	\n\t	vblendmpd		%%zmm10%{cdab%},%%zmm11,%%zmm11%{%%k2%}	\n\t"\
		"vblendmpd		%%zmm3%{cdab%},%%zmm2,%%zmm0%{%%k1%}	\n\t	vblendmpd		%%zmm13%{cdab%},%%zmm12,%%zmm10%{%%k1%}	\n\t"\
		"vblendmpd		%%zmm2%{cdab%},%%zmm3,%%zmm3%{%%k2%}	\n\t	vblendmpd		%%zmm12%{cdab%},%%zmm13,%%zmm13%{%%k2%}	\n\t"\
		"vblendmpd		%%zmm5%{cdab%},%%zmm4,%%zmm2%{%%k1%}	\n\t	vblendmpd		%%zmm15%{cdab%},%%zmm14,%%zmm12%{%%k1%}	\n\t"\
		"vblendmpd		%%zmm4%{cdab%},%%zmm5,%%zmm5%{%%k2%}	\n\t	vblendmpd		%%zmm14%{cdab%},%%zmm15,%%zmm15%{%%k2%}	\n\t"\
		"vblendmpd		%%zmm7%{cdab%},%%zmm6,%%zmm4%{%%k1%}	\n\t	vblendmpd		%%zmm17%{cdab%},%%zmm16,%%zmm14%{%%k1%}	\n\t"\
		"vblendmpd		%%zmm6%{cdab%},%%zmm7,%%zmm7%{%%k2%}	\n\t	vblendmpd		%%zmm16%{cdab%},%%zmm17,%%zmm17%{%%k2%}	\n\t"\
		"\n\t"\
		"vblendmpd		%%zmm0%{badc%},%%zmm8,%%zmm6%{%%k3%}	\n\t	vblendmpd		%%zmm10%{badc%},%%zmm9 ,%%zmm16%{%%k3%}	\n\t"\
		"vblendmpd		%%zmm8%{badc%},%%zmm0,%%zmm0%{%%k4%}	\n\t	vblendmpd		%%zmm9 %{badc%},%%zmm10,%%zmm10%{%%k4%}	\n\t"\
		"vblendmpd		%%zmm3%{badc%},%%zmm1,%%zmm8%{%%k3%}	\n\t	vblendmpd		%%zmm13%{badc%},%%zmm11,%%zmm9 %{%%k3%}	\n\t"\
		"vblendmpd		%%zmm1%{badc%},%%zmm3,%%zmm3%{%%k4%}	\n\t	vblendmpd		%%zmm11%{badc%},%%zmm13,%%zmm13%{%%k4%}	\n\t"\
		"vblendmpd		%%zmm4%{badc%},%%zmm2,%%zmm1%{%%k3%}	\n\t	vblendmpd		%%zmm14%{badc%},%%zmm12,%%zmm11%{%%k3%}	\n\t"\
		"vblendmpd		%%zmm2%{badc%},%%zmm4,%%zmm4%{%%k4%}	\n\t	vblendmpd		%%zmm12%{badc%},%%zmm14,%%zmm14%{%%k4%}	\n\t"\
		"vblendmpd		%%zmm7%{badc%},%%zmm5,%%zmm2%{%%k3%}	\n\t	vblendmpd		%%zmm17%{badc%},%%zmm15,%%zmm12%{%%k3%}	\n\t"\
		"vblendmpd		%%zmm5%{badc%},%%zmm7,%%zmm7%{%%k4%}	\n\t	vblendmpd		%%zmm15%{badc%},%%zmm17,%%zmm17%{%%k4%}	\n\t"\
		"\n\t"\
		"vpermf32x4	$78,%%zmm1,%%zmm1							\n\t	vpermf32x4	$78,%%zmm11,%%zmm11							\n\t"\
		"vpermf32x4	$78,%%zmm2,%%zmm2							\n\t	vpermf32x4	$78,%%zmm12,%%zmm12							\n\t"\
		"vpermf32x4	$78,%%zmm4,%%zmm4							\n\t	vpermf32x4	$78,%%zmm14,%%zmm14							\n\t"\
		"vpermf32x4	$78,%%zmm7,%%zmm7							\n\t	vpermf32x4	$78,%%zmm17,%%zmm17							\n\t"\
		"vblendmpd		%%zmm1,%%zmm6,%%zmm5%{%%k5%}			\n\t	vblendmpd		%%zmm11,%%zmm16,%%zmm15%{%%k5%}			\n\t"\
		"vblendmpd		%%zmm6,%%zmm1,%%zmm1%{%%k5%}			\n\t	vblendmpd		%%zmm16,%%zmm11,%%zmm11%{%%k5%}			\n\t"\
		"vblendmpd		%%zmm2,%%zmm8,%%zmm6%{%%k5%}			\n\t	vblendmpd		%%zmm12,%%zmm9 ,%%zmm16%{%%k5%}			\n\t"\
		"vblendmpd		%%zmm8,%%zmm2,%%zmm2%{%%k5%}			\n\t	vblendmpd		%%zmm9 ,%%zmm12,%%zmm12%{%%k5%}			\n\t"\
		"vblendmpd		%%zmm4,%%zmm0,%%zmm8%{%%k5%}			\n\t	vblendmpd		%%zmm14,%%zmm10,%%zmm9 %{%%k5%}			\n\t"\
		"vblendmpd		%%zmm0,%%zmm4,%%zmm4%{%%k5%}			\n\t	vblendmpd		%%zmm10,%%zmm14,%%zmm14%{%%k5%}			\n\t"\
		"vblendmpd		%%zmm7,%%zmm3,%%zmm0%{%%k5%}			\n\t	vblendmpd		%%zmm17,%%zmm13,%%zmm10%{%%k5%}			\n\t"\
		"vblendmpd		%%zmm3,%%zmm7,%%zmm7%{%%k5%}			\n\t	vblendmpd		%%zmm13,%%zmm17,%%zmm17%{%%k5%}			\n\t"\
		"vpermf32x4	$78,%%zmm1,%%zmm1							\n\t	vpermf32x4	$78,%%zmm11,%%zmm11							\n\t"\
		"vpermf32x4	$78,%%zmm2,%%zmm2							\n\t	vpermf32x4	$78,%%zmm12,%%zmm12							\n\t"\
		"vpermf32x4	$78,%%zmm4,%%zmm4							\n\t	vpermf32x4	$78,%%zmm14,%%zmm14							\n\t"\
		"vpermf32x4	$78,%%zmm7,%%zmm7							\n\t	vpermf32x4	$78,%%zmm17,%%zmm17							\n\t"\
		"\n\t"\
		"vmovaps		%%zmm5,(%%rax)			\n\t	vmovaps	%%zmm15,0x40(%%rax)	\n\t"\
		"vmovaps		%%zmm6,(%%rbx)			\n\t	vmovaps	%%zmm16,0x40(%%rbx)	\n\t"\
		"vmovaps		%%zmm8,(%%rcx)			\n\t	vmovaps	%%zmm9 ,0x40(%%rcx)	\n\t"\
		"vmovaps		%%zmm0,(%%rdx)			\n\t	vmovaps	%%zmm10,0x40(%%rdx)	\n\t"\
		"vmovaps		%%zmm1,(%%r10)			\n\t	vmovaps	%%zmm11,0x40(%%r10)	\n\t"\
		"vmovaps		%%zmm2,(%%r11)			\n\t	vmovaps	%%zmm12,0x40(%%r11)	\n\t"\
		"vmovaps		%%zmm4,(%%r12)			\n\t	vmovaps	%%zmm14,0x40(%%r12)	\n\t"\
		"vmovaps		%%zmm7,(%%r13)			\n\t	vmovaps	%%zmm17,0x40(%%r13)	\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__add1] "m" (Xadd1)\
		 ,[__add2] "m" (Xadd2)\
		 ,[__add3] "m" (Xadd3)\
		 ,[__add4] "m" (Xadd4)\
		 ,[__add5] "m" (Xadd5)\
		 ,[__add6] "m" (Xadd6)\
		 ,[__add7] "m" (Xadd7)\
		 ,[__r1] "m" (Xr1)\
		 ,[__r9] "m" (Xr9)\
		 ,[__r17] "m" (Xr17)\
		 ,[__r25] "m" (Xr25)\
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__cc0] "m" (Xcc0)\
		 ,[__c1] "m" (Xc1)\
		 ,[__c2] "m" (Xc2)\
		 ,[__c3] "m" (Xc3)\
		 ,[__c4] "m" (Xc4)\
		 ,[__c5] "m" (Xc5)\
		 ,[__c6] "m" (Xc6)\
		 ,[__c7] "m" (Xc7)\
		 ,[__c8] "m" (Xc8)\
		 ,[__c9] "m" (Xc9)\
		 ,[__c10] "m" (Xc10)\
		 ,[__c11] "m" (Xc11)\
		 ,[__c12] "m" (Xc12)\
		 ,[__c13] "m" (Xc13)\
		 ,[__c14] "m" (Xc14)\
		 ,[__c15] "m" (Xc15)\
	 ,[__pfetch_dist] "m" (Xpfetch_dist)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r10","r11","r12","r13","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15","xmm16","xmm17", "xmm29","xmm30","xmm31"	/* Clobbered registers */\
	);\
	}

  #else	// AVX-512 version:

	// Since the add0/1 block addresses advance in opposite directions, our prefetch scheme here is like so:
	// In fwd DIF prefetch  ahead 8 cache lines' worth of data w.r.to add0;
	// In inv DIT prefetch behind 8 cache lines' worth of data w.r.to add1;
	/****** To-Do: Need a 16-register version! ****************/
	#define SSE2_RADIX16_WRAPPER_DIT(Xadd0,Xadd1,Xadd2,Xadd3,Xadd4,Xadd5,Xadd6,Xadd7,Xr1,Xr9,Xr17,Xr25,Xisrt2,Xcc0,Xc1,Xc2,Xc3,Xc4,Xc5,Xc6,Xc7,Xc8,Xc9,Xc10,Xc11,Xc12,Xc13,Xc14,Xc15,Xpfetch_dist)\
	{\
	__asm__ volatile (\
	/*...Block 1: r1,9,17,25 */\
		"movq	%[__r1],%%rax					\n\t"\
		"leaq	0x10c0(%%rax),%%r10		\n\t"/* &one */\
		"vmovaps	0x800(%%rax),%%zmm29	\n\t"/* isrt2 */\
		"vmovaps	     (%%r10),%%zmm30	\n\t"/* 1.0 */\
		"vmovaps	0x040(%%r10),%%zmm31	\n\t"/* 2.0 */\
	"movq	%[__add1],%%r14	\n\t"\
	"movslq	%[__pfetch_dist],%%r14	\n\t"\
	"leaq	(%%rax,%%r14,8),%%r14	\n\t"	/* Block 1 [base-address + data-fetch-ahead index] */\
	"prefetcht1	-0x080(%%r14)\n\t"\
		"movq	%%rax,%%rbx						\n\t		leaq	0x1100(%%rax),%%rsi	\n\t"/* two */\
		"movq	%%rax,%%rcx						\n\t"	/*...Block 3: r3,11,19,27 */\
		"movq	%%rax,%%rdx						\n\t		leaq	0x080(%%rax),%%r10				\n\t"\
		"addq	$0x400,%%rbx					\n\t		leaq	0x080(%%rbx),%%r11				\n\t"\
		"addq	$0x200,%%rcx					\n\t		leaq	0x080(%%rcx),%%r12				\n\t"\
		"addq	$0x600,%%rdx					\n\t		leaq	0x080(%%rdx),%%r13				\n\t"\
		"/* SSE2_RADIX4_DIT_IN_PLACE(): */		\n\t		/* SSE2_RADIX4_DIT_IN_PLACE(): */		\n\t"\
		"vmovaps	     (%%rax),%%zmm2			\n\t		vmovaps	     (%%r10),%%zmm10			\n\t"\
		"vmovaps	0x040(%%rax),%%zmm3			\n\t		vmovaps	0x040(%%r10),%%zmm11			\n\t"\
		"vmovaps	     (%%rbx),%%zmm0			\n\t		vmovaps	     (%%r11),%%zmm8 			\n\t"\
		"vmovaps	0x040(%%rbx),%%zmm1			\n\t		vmovaps	0x040(%%r11),%%zmm9 			\n\t"\
		"vsubpd		%%zmm0,%%zmm2,%%zmm2		\n\t		vsubpd		%%zmm8 ,%%zmm10,%%zmm10		\n\t"\
		"vsubpd		%%zmm1,%%zmm3,%%zmm3		\n\t		vsubpd		%%zmm9 ,%%zmm11,%%zmm11		\n\t"\
		"vmovaps	     (%%rcx),%%zmm6			\n\t		vmovaps	     (%%r12),%%zmm14			\n\t"\
		"vmovaps	0x040(%%rcx),%%zmm7			\n\t		vmovaps	0x040(%%r12),%%zmm15			\n\t"\
		"vmovaps	     (%%rdx),%%zmm4			\n\t		vmovaps	     (%%r13),%%zmm12			\n\t"\
		"vmovaps	0x040(%%rdx),%%zmm5			\n\t		vmovaps	0x040(%%r13),%%zmm13			\n\t"\
		"vsubpd		%%zmm4,%%zmm6,%%zmm6		\n\t		vsubpd		 %%zmm12,%%zmm14,%%zmm14	\n\t"\
		"vsubpd		%%zmm5,%%zmm7,%%zmm7		\n\t		vsubpd	0x040(%%r13),%%zmm15,%%zmm15	\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm2,%%zmm0		\n\t	vfmadd132pd		%%zmm31,%%zmm10,%%zmm8 		\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm3,%%zmm1		\n\t	vfmadd132pd		%%zmm31,%%zmm11,%%zmm9 		\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm6,%%zmm4		\n\t	vfmadd132pd		%%zmm31,%%zmm14,%%zmm12		\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm7,%%zmm5		\n\t	vfmadd132pd		%%zmm31,%%zmm15,%%zmm13		\n\t"\
	"prefetcht1	-0x180(%%r14)\n\t"\
		"vsubpd	      %%zmm4,%%zmm0,%%zmm0		\n\t		vsubpd	      %%zmm12,%%zmm8 ,%%zmm8 	\n\t"\
		"vsubpd	      %%zmm5,%%zmm1,%%zmm1		\n\t		vsubpd	      %%zmm13,%%zmm9 ,%%zmm9 	\n\t"\
		"vsubpd	      %%zmm7,%%zmm2,%%zmm2		\n\t		vsubpd	      %%zmm15,%%zmm10,%%zmm10	\n\t"\
		"vsubpd	      %%zmm6,%%zmm3,%%zmm3		\n\t		vsubpd	      %%zmm14,%%zmm11,%%zmm11	\n\t"\
		"vmovaps	%%zmm0,     (%%rbx)			\n\t		vmovaps	%%zmm8 ,     (%%r11)			\n\t"\
		"vmovaps	%%zmm1,0x040(%%rbx)			\n\t		vmovaps	%%zmm9 ,0x040(%%r11)			\n\t"\
		"vmovaps	%%zmm2,     (%%rdx)			\n\t		vmovaps	%%zmm10,     (%%r13)			\n\t"\
		"vmovaps	%%zmm3,0x040(%%rcx)			\n\t		vmovaps	%%zmm11,0x040(%%r12)			\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm0,%%zmm4		\n\t	vfmadd132pd	%%zmm31,%%zmm8 ,%%zmm12		\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm1,%%zmm5		\n\t	vfmadd132pd	%%zmm31,%%zmm9 ,%%zmm13		\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm2,%%zmm7		\n\t	vfmadd132pd	%%zmm31,%%zmm10,%%zmm15		\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm3,%%zmm6		\n\t	vfmadd132pd	%%zmm31,%%zmm11,%%zmm14		\n\t"\
		"vmovaps	%%zmm4,     (%%rax)			\n\t		vmovaps	%%zmm12,     (%%r10)			\n\t"\
		"vmovaps	%%zmm5,0x040(%%rax)			\n\t		vmovaps	%%zmm13,0x040(%%r10)			\n\t"\
		"vmovaps	%%zmm7,     (%%rcx)			\n\t		vmovaps	%%zmm15,     (%%r12)			\n\t"\
		"vmovaps	%%zmm6,0x040(%%rdx)			\n\t		vmovaps	%%zmm14,0x040(%%r13)			\n\t"\
	/*...Block 2: r5,13,21,29 */						/*...Block 4: r7,15,23,31 */\
		"addq	$0x100,%%rax						\n\t		leaq	0x080(%%rax),%%r10				\n\t"\
		"addq	$0x100,%%rbx						\n\t		leaq	0x080(%%rbx),%%r11				\n\t"\
		"addq	$0x100,%%rcx						\n\t		leaq	0x080(%%rcx),%%r12				\n\t"\
		"addq	$0x100,%%rdx						\n\t		leaq	0x080(%%rdx),%%r13				\n\t"\
		"/* SSE2_RADIX4_DIT_IN_PLACE(): */		\n\t		/* SSE2_RADIX4_DIT_IN_PLACE(): */		\n\t"\
		"vmovaps	     (%%rax),%%zmm2			\n\t		vmovaps	     (%%r10),%%zmm10			\n\t"\
		"vmovaps	0x040(%%rax),%%zmm3			\n\t		vmovaps	0x040(%%r10),%%zmm11			\n\t"\
		"vmovaps	     (%%rbx),%%zmm0			\n\t		vmovaps	     (%%r11),%%zmm8 			\n\t"\
		"vmovaps	0x040(%%rbx),%%zmm1			\n\t		vmovaps	0x040(%%r11),%%zmm9 			\n\t"\
		"vsubpd		%%zmm0,%%zmm2,%%zmm2		\n\t		vsubpd		%%zmm8 ,%%zmm10,%%zmm10		\n\t"\
		"vsubpd		%%zmm1,%%zmm3,%%zmm3		\n\t		vsubpd		%%zmm9 ,%%zmm11,%%zmm11		\n\t"\
		"vmovaps	     (%%rcx),%%zmm6			\n\t		vmovaps	     (%%r12),%%zmm14			\n\t"\
		"vmovaps	0x040(%%rcx),%%zmm7			\n\t		vmovaps	0x040(%%r12),%%zmm15			\n\t"\
		"vmovaps	     (%%rdx),%%zmm4			\n\t		vmovaps	     (%%r13),%%zmm12			\n\t"\
		"vmovaps	0x040(%%rdx),%%zmm5			\n\t		vmovaps	0x040(%%r13),%%zmm13			\n\t"\
		"vsubpd		%%zmm4,%%zmm6,%%zmm6		\n\t		vsubpd		 %%zmm12,%%zmm14,%%zmm14	\n\t"\
		"vsubpd		%%zmm5,%%zmm7,%%zmm7		\n\t		vsubpd	0x040(%%r13),%%zmm15,%%zmm15	\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm2,%%zmm0		\n\t	vfmadd132pd		%%zmm31,%%zmm10,%%zmm8 		\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm3,%%zmm1		\n\t	vfmadd132pd		%%zmm31,%%zmm11,%%zmm9 		\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm6,%%zmm4		\n\t	vfmadd132pd		%%zmm31,%%zmm14,%%zmm12		\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm7,%%zmm5		\n\t	vfmadd132pd		%%zmm31,%%zmm15,%%zmm13		\n\t"\
	"prefetcht1	-0x100(%%r14)\n\t"\
		"vsubpd	      %%zmm4,%%zmm0,%%zmm0		\n\t		vsubpd	      %%zmm12,%%zmm8 ,%%zmm8 	\n\t"\
		"vsubpd	      %%zmm5,%%zmm1,%%zmm1		\n\t		vsubpd	      %%zmm13,%%zmm9 ,%%zmm9 	\n\t"\
		"vsubpd	      %%zmm7,%%zmm2,%%zmm2		\n\t		vsubpd	      %%zmm15,%%zmm10,%%zmm10	\n\t"\
		"vsubpd	      %%zmm6,%%zmm3,%%zmm3		\n\t		vsubpd	      %%zmm14,%%zmm11,%%zmm11	\n\t"\
		"vmovaps	%%zmm0,     (%%rbx)			\n\t		vmovaps	%%zmm8 ,     (%%r11)			\n\t"\
		"vmovaps	%%zmm1,0x040(%%rbx)			\n\t		vmovaps	%%zmm9 ,0x040(%%r11)			\n\t"\
		"vmovaps	%%zmm2,     (%%rdx)			\n\t		vmovaps	%%zmm10,     (%%r13)			\n\t"\
		"vmovaps	%%zmm3,0x040(%%rcx)			\n\t		vmovaps	%%zmm11,0x040(%%r12)			\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm0,%%zmm4		\n\t	vfmadd132pd	%%zmm31,%%zmm8 ,%%zmm12		\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm1,%%zmm5		\n\t	vfmadd132pd	%%zmm31,%%zmm9 ,%%zmm13		\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm2,%%zmm7		\n\t	vfmadd132pd	%%zmm31,%%zmm10,%%zmm15		\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm3,%%zmm6		\n\t	vfmadd132pd	%%zmm31,%%zmm11,%%zmm14		\n\t"\
		"vmovaps	%%zmm4,     (%%rax)			\n\t		vmovaps	%%zmm12,     (%%r10)			\n\t"\
		"vmovaps	%%zmm5,0x040(%%rax)			\n\t		vmovaps	%%zmm13,0x040(%%r10)			\n\t"\
		"vmovaps	%%zmm7,     (%%rcx)			\n\t		vmovaps	%%zmm15,     (%%r12)			\n\t"\
		"vmovaps	%%zmm6,0x040(%%rdx)			\n\t		vmovaps	%%zmm14,0x040(%%r13)			\n\t"\
		"\n\t"\
	/***************************************************************************************************/\
	/*..and now do four more radix-4 transforms, including the internal and external twiddle factors.  */\
	/*  Write even-index 16-byte output pairs to a(j1), odd-index to a(j2), unpack same as on inputs.  */\
	/*  We do the last 2 radix-4 blocks first, to make the unpack-interleaving of outputs smoother.    */\
	/***************************************************************************************************/\
	/* Main-array addresses still in add0,1; all rax-offsets incr +0x400 in rcol w.r.to lcol: */\
	"prefetcht1	-0x200(%%r14)\n\t"\
	/*...Block 3: t3,11,19,27 -> r9,13,11,15: */		/*...Block 4: t7,15,23,31 -> r25,29,27,31: */\
		"movq		%[__r9],%%rax				\n\t"\
		"movq		%[__cc0],%%rcx				\n\t"\
	"vmovaps	(%%rcx),%%zmm3	\n\t	vmovaps	0x040(%%rcx),%%zmm11	\n\t"/* cc0,ss0 data shared by both columns */\
		"vmovaps	0x080(%%rax),%%zmm4			\n\t		vmovaps		0x480(%%rax),%%zmm12			\n\t"\
		"vmovaps	0x180(%%rax),%%zmm0			\n\t		vmovaps		0x580(%%rax),%%zmm8 			\n\t"\
		"vmovaps	0x0c0(%%rax),%%zmm5			\n\t		vmovaps		0x4c0(%%rax),%%zmm13			\n\t"\
		"vmovaps	0x1c0(%%rax),%%zmm1			\n\t		vmovaps		0x5c0(%%rax),%%zmm9 			\n\t"\
		"vmulpd		     %%zmm3 ,%%zmm4,%%zmm4	\n\t		vmulpd		     %%zmm11,%%zmm12,%%zmm12	\n\t"\
		"vmulpd		     %%zmm11,%%zmm0,%%zmm0	\n\t		vmulpd		     %%zmm3 ,%%zmm8 ,%%zmm8 	\n\t"\
		"vmulpd		     %%zmm3 ,%%zmm5,%%zmm5	\n\t		vmulpd		     %%zmm11,%%zmm13,%%zmm13	\n\t"\
		"vmulpd		     %%zmm11,%%zmm1,%%zmm1	\n\t		vmulpd		     %%zmm3 ,%%zmm9 ,%%zmm9 	\n\t"\
		"vmovaps	0x080(%%rax),%%zmm6			\n\t		vmovaps		0x480(%%rax),%%zmm14			\n\t"\
		"vmovaps	0x180(%%rax),%%zmm2			\n\t		vmovaps		0x580(%%rax),%%zmm10			\n\t"\
		"vmovaps	0x0c0(%%rax),%%zmm7			\n\t		vmovaps		0x4c0(%%rax),%%zmm15			\n\t"\
	"vfnmadd231pd	     %%zmm11,%%zmm6,%%zmm5	\n\t	vfnmadd231pd	     %%zmm3 ,%%zmm14,%%zmm13	\n\t"\
	"vfnmadd231pd	     %%zmm3 ,%%zmm2,%%zmm1	\n\t	vfnmadd231pd	     %%zmm11,%%zmm10,%%zmm9 	\n\t"\
	" vfmadd231pd	     %%zmm11,%%zmm7,%%zmm4	\n\t	 vfmadd231pd	     %%zmm3 ,%%zmm15,%%zmm12	\n\t"\
	" vfmadd231pd	0x1c0(%%rax),%%zmm3,%%zmm0	\n\t	 vfmadd231pd	0x5c0(%%rax),%%zmm11,%%zmm8 	\n\t"\
		"vmovaps	%%zmm5,%%zmm7				\n\t		vmovaps		%%zmm13,%%zmm15				\n\t"\
		"vmovaps	%%zmm4,%%zmm6				\n\t		vmovaps		%%zmm12,%%zmm14				\n\t"\
		"vaddpd		%%zmm0,%%zmm4,%%zmm4		\n\t		vaddpd		%%zmm8 ,%%zmm12,%%zmm12		\n\t"\
		"vaddpd		%%zmm1,%%zmm5,%%zmm5		\n\t		vaddpd		%%zmm9 ,%%zmm13,%%zmm13		\n\t"\
		"vsubpd		%%zmm0,%%zmm6,%%zmm6		\n\t		vsubpd		%%zmm8 ,%%zmm14,%%zmm14		\n\t"\
		"vsubpd		%%zmm1,%%zmm7,%%zmm7		\n\t		vsubpd		%%zmm9 ,%%zmm15,%%zmm15		\n\t"\
		"vmovaps	0x100(%%rax),%%zmm2			\n\t		vmovaps		0x500(%%rax),%%zmm10		\n\t"\
		"vmovaps	0x140(%%rax),%%zmm3			\n\t		vmovaps		0x540(%%rax),%%zmm11		\n\t"\
		"vmovaps	     (%%rax),%%zmm0			\n\t		vmovaps		0x400(%%rax),%%zmm8 		\n\t"\
		"vmovaps	0x040(%%rax),%%zmm1			\n\t		vmovaps		0x440(%%rax),%%zmm9 		\n\t"\
		"vaddpd		0x140(%%rax),%%zmm2,%%zmm2	\n\t		vsubpd		0x540(%%rax),%%zmm10,%%zmm10\n\t"\
		"vsubpd		0x100(%%rax),%%zmm3,%%zmm3	\n\t		vaddpd		0x500(%%rax),%%zmm11,%%zmm11\n\t"\
		"vmulpd		     %%zmm29,%%zmm2,%%zmm2	\n\t		vmulpd		%%zmm29,%%zmm10,%%zmm10		\n\t"\
		"vmulpd		     %%zmm29,%%zmm3,%%zmm3	\n\t		vmulpd		%%zmm29,%%zmm11,%%zmm11		\n\t"\
		"vsubpd		%%zmm2,%%zmm0,%%zmm0		\n\t		vsubpd		%%zmm10,%%zmm8 ,%%zmm8 		\n\t"\
		"vsubpd		%%zmm3,%%zmm1,%%zmm1		\n\t		vsubpd		%%zmm11,%%zmm9 ,%%zmm9 		\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm0,%%zmm2		\n\t	vfmadd132pd	%%zmm31,%%zmm8 ,%%zmm10			\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm1,%%zmm3		\n\t	vfmadd132pd	%%zmm31,%%zmm9 ,%%zmm11			\n\t"\
	"addq	$0x480,%%rcx		\n\t"/* c1 from cc0 */\
	"leaq	0x80(%%rcx),%%rdx	\n\t"/* c9 from c1  */\
	/* Using the upper block(s) of the main array for temp-storage in the section below led to a nasty AVX bug to track down: */\
	/* [***The following details are the AVX-512-ification of the the analogous comment in the AVX version of this macro ***]
	In fermat-mod mode the 8 block addresses in ascending order are add0-7 with no 'gaps' between blocks, whereas for mersenne-mod
	the addresses in asc. order are add0,2,4,6,7,5,3,1 with a gap between contiguous-data-block quartets 0,2,4,6 and 7,5,3,1. Thus
	for fermat-mod we need [add4] as the base-address of the 'high-half' block for temp-storage; for mersenne-mod we need [add7]
	In both cases we have that add4 < add7 so instead use (add3 - add1): > 0 for fermat-mod, < 0 for mersenne - to differentiate: */\
		"movq	%[__add4],%%rdi		\n\t"/* destroyable copy of add4 */\
		"movq	%[__add4],%%rbx		\n\t"\
		"subq	%[__add3],%%rdi		\n\t"/* rdi = (add4 - add3); if this yields a borrow (i.e. sets CF) it's mersenne, else fermat. */\
		"cmovcq %[__add7],%%rbx	\n\t" /* if CF set (i.e. h > l), copy source [add7] into dest (rbx), else leave dest = [add4]. */\
		/* rbx shared between rcol/lcol; rcx/rdx-offsets incr +0x100 in rcol for rest of block: */\
	"prefetcht1	-0x280(%%r14)\n\t"\
		"vsubpd		%%zmm4,%%zmm2,%%zmm2			\n\t		vsubpd		%%zmm14,%%zmm8 ,%%zmm8 		\n\t"\
		"vsubpd		%%zmm5,%%zmm3,%%zmm3			\n\t		vsubpd		%%zmm15,%%zmm9 ,%%zmm9 		\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm2,%%zmm4			\n\t	vfmadd132pd	%%zmm31,%%zmm8 ,%%zmm14			\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm3,%%zmm5			\n\t	vfmadd132pd	%%zmm31,%%zmm9 ,%%zmm15			\n\t"\
		"vmovaps	%%zmm2,     (%%rax)				\n\t		vmovaps		%%zmm8 ,0x400(%%rax)		\n\t"\
		"vmovaps	%%zmm3,0x040(%%rax)				\n\t		vmovaps		%%zmm9 ,0x440(%%rax)		\n\t"\
		"vmovaps	%%zmm4,%%zmm2					\n\t		vmovaps		%%zmm14,%%zmm8 		\n\t"\
		"vmovaps	%%zmm5,%%zmm3					\n\t		vmovaps		%%zmm15,%%zmm9 		\n\t"\
		"vmulpd		     (%%rcx),%%zmm4,%%zmm4		\n\t		vmulpd		0x200(%%rcx),%%zmm14,%%zmm14\n\t"\
		"vmulpd		     (%%rcx),%%zmm5,%%zmm5		\n\t		vmulpd		0x200(%%rcx),%%zmm15,%%zmm15\n\t"\
		"vsubpd		%%zmm7,%%zmm0,%%zmm0			\n\t		vsubpd		%%zmm13,%%zmm10,%%zmm10		\n\t"\
		"vsubpd		%%zmm6,%%zmm1,%%zmm1			\n\t		vsubpd		%%zmm12,%%zmm11,%%zmm11		\n\t"\
	"vfnmadd231pd	0x040(%%rcx),%%zmm2,%%zmm5		\n\t	vfnmadd231pd	0x240(%%rcx),%%zmm8 ,%%zmm15	\n\t"\
	" vfmadd231pd	0x040(%%rcx),%%zmm3,%%zmm4		\n\t	 vfmadd231pd	0x240(%%rcx),%%zmm9 ,%%zmm14	\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm0,%%zmm7			\n\t	vfmadd132pd	%%zmm31,%%zmm10,%%zmm13			\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm1,%%zmm6			\n\t	vfmadd132pd	%%zmm31,%%zmm11,%%zmm12			\n\t"\
		"vmovaps	%%zmm5,0x040(%%rbx)				\n\t		vmovaps		%%zmm15,0x0c0(%%rbx)		\n\t"\
		"vmovaps	%%zmm4,     (%%rbx)				\n\t		vmovaps		%%zmm14,0x080(%%rbx)		\n\t"\
		"vmovaps		 (%%rax),%%zmm4				\n\t		vmovaps		0x400(%%rax),%%zmm14		\n\t"\
		"vmovaps	0x040(%%rax),%%zmm5				\n\t		vmovaps		0x440(%%rax),%%zmm15		\n\t"\
		"vmovaps	%%zmm4,%%zmm2					\n\t		vmovaps		%%zmm14,%%zmm8 		\n\t"\
		"vmovaps	%%zmm5,%%zmm3					\n\t		vmovaps		%%zmm15,%%zmm9 		\n\t"\
		"vmulpd		     (%%rdx),%%zmm4,%%zmm4		\n\t		vmulpd		0x200(%%rdx),%%zmm14,%%zmm14\n\t"\
		"vmulpd		     (%%rdx),%%zmm5,%%zmm5		\n\t		vmulpd		0x200(%%rdx),%%zmm15,%%zmm15\n\t"\
	"vfnmadd231pd	0x040(%%rdx),%%zmm2,%%zmm5		\n\t	vfnmadd231pd	0x240(%%rdx),%%zmm8 ,%%zmm15	\n\t"\
	" vfmadd231pd	0x040(%%rdx),%%zmm3,%%zmm4		\n\t	 vfmadd231pd	0x240(%%rdx),%%zmm9 ,%%zmm14	\n\t"\
		"vmovaps	%%zmm5,0x240(%%rbx)				\n\t		vmovaps		%%zmm15,0x2c0(%%rbx)		\n\t"\
		"vmovaps	%%zmm4,0x200(%%rbx)				\n\t		vmovaps		%%zmm14,0x280(%%rbx)		\n\t"\
		"addq	$0x100,%%rcx	\n\t"/* c1  from c1; rcol has c7 */\
		"addq	$0x100,%%rdx	\n\t"/* c9  from c9; rcol has c15*/\
		"vmovaps	%%zmm7,%%zmm4					\n\t		vmovaps		%%zmm13,%%zmm8 		\n\t"\
		"vmovaps	%%zmm1,%%zmm5					\n\t		vmovaps		%%zmm11,%%zmm9 		\n\t"\
		"vmulpd		     (%%rcx),%%zmm7,%%zmm7		\n\t		vmulpd		0x200(%%rcx),%%zmm13,%%zmm13\n\t"\
		"vmulpd		     (%%rcx),%%zmm1,%%zmm1		\n\t		vmulpd		0x200(%%rcx),%%zmm11,%%zmm11\n\t"\
	"vfnmadd231pd	0x040(%%rcx),%%zmm4,%%zmm1		\n\t	vfnmadd231pd	0x240(%%rcx),%%zmm8 ,%%zmm11	\n\t"\
	" vfmadd231pd	0x040(%%rcx),%%zmm5,%%zmm7		\n\t	 vfmadd231pd	0x240(%%rcx),%%zmm9 ,%%zmm13	\n\t"\
		"vmovaps	%%zmm0,%%zmm4					\n\t		vmovaps		%%zmm10,%%zmm8 		\n\t"\
		"vmovaps	%%zmm6,%%zmm5					\n\t		vmovaps		%%zmm12,%%zmm9 		\n\t"\
		"vmulpd		     (%%rdx),%%zmm0,%%zmm0		\n\t		vmulpd		0x200(%%rdx),%%zmm10,%%zmm10\n\t"\
		"vmulpd		     (%%rdx),%%zmm6,%%zmm6		\n\t		vmulpd		0x200(%%rdx),%%zmm12,%%zmm12\n\t"\
		"vmovaps	%%zmm1,0x140(%%rbx)				\n\t		vmovaps		%%zmm11,0x1c0(%%rbx)		\n\t"\
		"vmovaps	%%zmm7,0x100(%%rbx)				\n\t		vmovaps		%%zmm13,0x180(%%rbx)		\n\t"\
	"vfnmadd231pd	0x040(%%rdx),%%zmm4,%%zmm6		\n\t	vfnmadd231pd	0x240(%%rdx),%%zmm8 ,%%zmm12	\n\t"\
	" vfmadd231pd	0x040(%%rdx),%%zmm5,%%zmm0		\n\t	 vfmadd231pd	0x240(%%rdx),%%zmm9 ,%%zmm10	\n\t"\
		"vmovaps	%%zmm6,0x340(%%rbx)				\n\t		vmovaps		%%zmm12,0x3c0(%%rbx)		\n\t"\
		"vmovaps	%%zmm0,0x300(%%rbx)				\n\t		vmovaps		%%zmm10,0x380(%%rbx)		\n\t"\
	"prefetcht1	-0x300(%%r14)						\n\t"\
		/*...Block 1: t1,9,17,25 -> r1,5,3,7: */			/*...Block 3: t5,13,21,29 -> r17,21,19,23: */\
		"movq	%[__r1],%%rax	\n\t"/* r17 in rcol */\
		"vmovaps	     (%%rax),%%zmm0				\n\t		vmovaps		0x480(%%rax),%%zmm12	\n\t"\
		"vmovaps	0x040(%%rax),%%zmm1				\n\t		vmovaps		0x4c0(%%rax),%%zmm13		\n\t"\
		"vmovaps	0x100(%%rax),%%zmm2				\n\t		vmovaps		0x580(%%rax),%%zmm8 		\n\t"\
		"vmovaps	0x140(%%rax),%%zmm3				\n\t		vmovaps		0x5c0(%%rax),%%zmm9 		\n\t"\
		"vsubpd		%%zmm2,%%zmm0,%%zmm0			\n\t		vaddpd		0x4c0(%%rax),%%zmm12,%%zmm12\n\t"\
		"vsubpd		%%zmm3,%%zmm1,%%zmm1			\n\t		vsubpd		0x480(%%rax),%%zmm13,%%zmm13\n\t"\
		"vmovaps	0x080(%%rax),%%zmm4				\n\t		vsubpd		0x5c0(%%rax),%%zmm8 ,%%zmm8 \n\t"\
		"vmovaps	0x0c0(%%rax),%%zmm5				\n\t		vaddpd		0x580(%%rax),%%zmm9 ,%%zmm9 \n\t"\
		"vmovaps	0x180(%%rax),%%zmm6				\n\t		vmulpd		%%zmm29,%%zmm12,%%zmm12		\n\t"\
		"vmovaps	0x1c0(%%rax),%%zmm7				\n\t		vmulpd		%%zmm29,%%zmm13,%%zmm13		\n\t"\
		"vsubpd		%%zmm6,%%zmm4,%%zmm4			\n\t		vmovaps		%%zmm12,%%zmm14		\n\t"\
		"vsubpd		%%zmm7,%%zmm5,%%zmm5			\n\t		vmovaps		%%zmm13,%%zmm15		\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm0,%%zmm2			\n\t	vfnmadd231pd	%%zmm29,%%zmm8 ,%%zmm12		\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm1,%%zmm3			\n\t	vfnmadd231pd	%%zmm29,%%zmm9 ,%%zmm13		\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm4,%%zmm6			\n\t	 vfmadd231pd	%%zmm29,%%zmm8 ,%%zmm14		\n\t"\
	"vfmadd132pd	%%zmm31,%%zmm5,%%zmm7			\n\t	 vfmadd231pd	%%zmm29,%%zmm9 ,%%zmm15		\n\t"\
		"movq	%%rbx,%%rcx	\n\t"/* rbx saves 'high-half' block address; move that to rcx and put 'low-half' block address into rbx */\
		"movq	%[__add0],%%rbx						\n\t		vmovaps		0x400(%%rax),%%zmm8 		\n\t"\
		"subq	$0x500,%%rdx	/* c8 from c13 */	\n\t		vmovaps		0x440(%%rax),%%zmm9 		\n\t"\
		"vaddpd		%%zmm6,%%zmm2,%%zmm2			\n\t		vmovaps		0x500(%%rax),%%zmm10		\n\t"\
		"vaddpd		%%zmm7,%%zmm3,%%zmm3			\n\t		vmovaps		0x540(%%rax),%%zmm11		\n\t"\
		"vmovaps	%%zmm2,     (%%rbx)				\n\t		vsubpd		0x540(%%rax),%%zmm8 ,%%zmm8 \n\t"\
		"vmovaps	%%zmm3,0x040(%%rbx)				\n\t		vsubpd		0x500(%%rax),%%zmm9 ,%%zmm9 \n\t"\
	"vfnmadd231pd	%%zmm31,%%zmm6,%%zmm2			\n\t		vaddpd		0x400(%%rax),%%zmm11,%%zmm11\n\t"\
	"vfnmadd231pd	%%zmm31,%%zmm7,%%zmm3			\n\t		vaddpd		0x440(%%rax),%%zmm10,%%zmm10\n\t"\
		"vmovaps	%%zmm2,%%zmm6					\n\t		vsubpd		%%zmm12,%%zmm11,%%zmm11		\n\t"\
		"vmovaps	%%zmm3,%%zmm7					\n\t		vsubpd		%%zmm13,%%zmm9 ,%%zmm9 		\n\t"\
		"vmulpd		     (%%rdx),%%zmm2,%%zmm2		\n\t	vfmadd132pd	%%zmm31,%%zmm11,%%zmm12			\n\t"\
		"vmulpd		     (%%rdx),%%zmm3,%%zmm3		\n\t	vfmadd132pd	%%zmm31,%%zmm9 ,%%zmm13			\n\t"\
	"vfnmadd231pd	0x040(%%rdx),%%zmm6,%%zmm3		\n\t		vmovaps		%%zmm11,     (%%rax)		\n\t"\
	" vfmadd231pd	0x040(%%rdx),%%zmm7,%%zmm2		\n\t		vmovaps		%%zmm9 ,0x040(%%rax)		\n\t"\
		"vmovaps	0x240(%%rcx),%%zmm7				\n\t		vmovaps		%%zmm12,%%zmm11		\n\t"\
		"vmovaps	0x200(%%rcx),%%zmm6				\n\t		vmovaps		%%zmm13,%%zmm9 		\n\t"\
		"vmovaps	%%zmm3,0x240(%%rax)	/* r9 */	\n\t		vmulpd		0x180(%%rdx),%%zmm12,%%zmm12	\n\t"/* c2 */\
		"vmovaps	%%zmm2,0x200(%%rax)	/* r8 */	\n\t		vmulpd		0x180(%%rdx),%%zmm13,%%zmm13\n\t"\
		"vmovaps	%%zmm7,0x640(%%rax)	/* r25 */	\n\t	vfnmadd231pd	0x1c0(%%rdx),%%zmm11,%%zmm13\n\t"\
		"vmovaps	%%zmm6,0x600(%%rax)	/* r24 */	\n\t	 vfmadd231pd	0x1c0(%%rdx),%%zmm9 ,%%zmm12\n\t"\
		"vmovaps	0x040(%%rbx),%%zmm3				\n\t		vmovaps	0x0c0(%%rcx),%%zmm11	\n\t"\
		"vmovaps	     (%%rbx),%%zmm2				\n\t		vmovaps	0x080(%%rcx),%%zmm9		\n\t"\
		"vmovaps	0x040(%%rcx),%%zmm7				\n\t		vmovaps	%%zmm13,0x0c0(%%rax)	/* r3 */\n\t"\
		"vmovaps	     (%%rcx),%%zmm6				\n\t		vmovaps	%%zmm12,0x080(%%rax)	/* r2 */\n\t"\
	/* swap rcol/lcol in next 4 lines: */\
		"vmovaps		     (%%rax),%%zmm12		\n\t		vmovaps	%%zmm2 ,     (%%rax)	/* r0 */\n\t"\
		"vmovaps		0x040(%%rax),%%zmm13		\n\t		vmovaps	%%zmm3 ,0x040(%%rax)	/* r1 */\n\t"\
		"vmovaps	%%zmm11,0x4c0(%%rax)	/* r19*/\n\t		vmovaps		%%zmm12,%%zmm11		\n\t"\
		"vmovaps	%%zmm9 ,0x480(%%rax)	/* r18*/\n\t		vmovaps		%%zmm13,%%zmm9 		\n\t"\
		"														addq	$0x080,%%rdx		/* c4 in lcol; c10 in rcol*/		\n\t"\
		"vmovaps	%%zmm7,0x440(%%rax)	/* r17 */	\n\t		vmulpd		0x180(%%rdx),%%zmm12,%%zmm12		\n\t"\
		"vmovaps	%%zmm6,0x400(%%rax)	/* r16 */	\n\t		vmulpd		0x180(%%rdx),%%zmm13,%%zmm13		\n\t"\
		"vaddpd		%%zmm5,%%zmm0,%%zmm0			\n\t	vfnmadd231pd	0x1c0(%%rdx),%%zmm11,%%zmm13		\n\t"\
		"vsubpd		%%zmm4,%%zmm1,%%zmm1			\n\t	 vfmadd231pd	0x1c0(%%rdx),%%zmm9 ,%%zmm12		\n\t"\
		"vmovaps	%%zmm0,%%zmm2					\n\t		vmovaps	0x2c0(%%rcx),%%zmm11	\n\t"\
		"vmovaps	%%zmm1,%%zmm3					\n\t		vmovaps	0x280(%%rcx),%%zmm9		\n\t"\
		"			\n\t		vmovaps	%%zmm13,0x2c0(%%rax)	/* r11 */\n\t"\
		"			\n\t		vmovaps	%%zmm12,0x280(%%rax)	/* r10 */\n\t"\
		"vmovaps	%%zmm0,%%zmm6					\n\t		vmovaps	%%zmm11,0x6c0(%%rax)	/* r27 */\n\t"\
		"vmovaps	%%zmm1,%%zmm7					\n\t		vmovaps	%%zmm9 ,0x680(%%rax)	/* r26 */\n\t"\
	"vfnmadd231pd	%%zmm31,%%zmm5,%%zmm0			\n\t		vsubpd		%%zmm15,%%zmm8 ,%%zmm8 		\n\t"\
	" vfmadd231pd	%%zmm31,%%zmm4,%%zmm1			\n\t		vsubpd		%%zmm14,%%zmm10,%%zmm10		\n\t"\
		"vmulpd		     (%%rdx),%%zmm2,%%zmm2		\n\t	vfmadd132pd	%%zmm31,%%zmm8 ,%%zmm15		\n\t"\
		"vmulpd		     (%%rdx),%%zmm3,%%zmm3		\n\t	vfmadd132pd	%%zmm31,%%zmm10,%%zmm14		\n\t"\
	"vfnmadd231pd	0x040(%%rdx),%%zmm6,%%zmm3		\n\t		vmovaps		%%zmm15,%%zmm12		\n\t"\
	" vfmadd231pd	0x040(%%rdx),%%zmm7,%%zmm2		\n\t		vmovaps		%%zmm10,%%zmm13		\n\t"\
	"prefetcht1	-0x380(%%r14)						\n\t	addq	$0x080,%%rdx	\n\t"/* c12 in lcol; c6 in rcol*/\
		"vmovaps	0x140(%%rcx),%%zmm7				\n\t		vmulpd		0x180(%%rdx),%%zmm15,%%zmm15	\n\t"\
		"vmovaps	0x100(%%rcx),%%zmm6				\n\t		vmulpd		0x180(%%rdx),%%zmm10,%%zmm10	\n\t"\
		"vmovaps	%%zmm3,0x140(%%rax)	/* r5 */	\n\t	vfnmadd231pd	0x1c0(%%rdx),%%zmm12,%%zmm10	\n\t"\
		"vmovaps	%%zmm2,0x100(%%rax)	/* r4 */	\n\t	 vfmadd231pd	0x1c0(%%rdx),%%zmm13,%%zmm15	\n\t"\
		"vmovaps	%%zmm7,0x540(%%rax)	/* r21 */	\n\t		vmovaps	0x1c0(%%rcx),%%zmm13	\n\t"\
		"vmovaps	%%zmm6,0x500(%%rax)	/* r20 */	\n\t		vmovaps	0x180(%%rcx),%%zmm12	\n\t"\
		"vmovaps	%%zmm0,%%zmm6					\n\t		vmovaps	%%zmm10,0x1c0(%%rax)	/* r7 */\n\t"\
		"vmovaps	%%zmm1,%%zmm7					\n\t		vmovaps	%%zmm15,0x180(%%rax)	/* r6 */\n\t"\
		"vmulpd		     (%%rdx),%%zmm0,%%zmm0		\n\t		vmovaps	%%zmm13,0x5c0(%%rax)	/* r23 */\n\t"\
		"vmulpd		     (%%rdx),%%zmm1,%%zmm1		\n\t		vmovaps	%%zmm12,0x580(%%rax)	/* r22 */\n\t"\
	"vfnmadd231pd	0x040(%%rdx),%%zmm6,%%zmm1		\n\t		vmovaps		%%zmm8 ,%%zmm12		\n\t"\
	" vfmadd231pd	0x040(%%rdx),%%zmm7,%%zmm0		\n\t		vmovaps		%%zmm14,%%zmm13		\n\t"\
		"														vmulpd		0x200(%%rdx),%%zmm8 ,%%zmm8 \n\t"/* c14 */\
		"														vmulpd		0x200(%%rdx),%%zmm14,%%zmm14		\n\t"\
		"													vfnmadd231pd	0x240(%%rdx),%%zmm12,%%zmm14		\n\t"\
		"													 vfmadd231pd	0x240(%%rdx),%%zmm13,%%zmm8 		\n\t"\
		"vmovaps	0x340(%%rcx),%%zmm7				\n\t		vmovaps	0x3c0(%%rcx),%%zmm13	\n\t"\
		"vmovaps	0x300(%%rcx),%%zmm6				\n\t		vmovaps	0x380(%%rcx),%%zmm12	\n\t"\
		"vmovaps	%%zmm1,0x340(%%rax)	/* r13 */	\n\t		vmovaps	%%zmm14,0x3c0(%%rax)	/* r15 */\n\t"\
		"vmovaps	%%zmm0,0x300(%%rax)	/* r12 */	\n\t		vmovaps	%%zmm8 ,0x380(%%rax)	/* r14 */\n\t"\
		"vmovaps	%%zmm7,0x740(%%rax)	/* r29 */	\n\t		vmovaps	%%zmm13,0x7c0(%%rax)	/* r31 */\n\t"\
		"vmovaps	%%zmm6,0x700(%%rax)	/* r28 */	\n\t		vmovaps	%%zmm12,0x780(%%rax)	/* r30 */\n\t"\
	"prefetcht1	-0x400(%%r14)\n\t"\
	/*******************************************
	***** Finish with 8-way 'un'terleaving: ****
	Using the AVX-512 data layout, the rcol pattern is:
		a[ 0- 7] = re[ 0, 8, 1, 9, 2,10, 3,11].d0	a[ 8-15] = im[ 0, 8, 1, 9, 2,10, 3,11].d0
		a[16-23] = re[ 4,12, 5,13, 6,14, 7,15].d0	a[24-31] = im[ 4,12, 5,13, 6,14, 7,15].d0 ,
	and remaining seven 32-double blocks repeat same pattern with elts d1-d7 of the vector-doubles.
	*******************************************/\
		"movq	%[__add0],%%rax	\n\t"\
		"movq	%[__add1],%%rbx	\n\t"\
		"movq	%[__add2],%%rcx	\n\t"\
		"movq	%[__add3],%%rdx	\n\t"\
		"movq	%[__add4],%%r10	\n\t"\
		"movq	%[__add5],%%r11	\n\t"\
		"movq	%[__add6],%%r12	\n\t"\
		"movq	%[__add7],%%r13	\n\t"\
	"movq	%[__r1] ,%%rsi	\n\t"\
	/**** a[ 0- 7] = re[ 0, 8, 1, 9, 2,10, 3,11].d0	a[ 8-15] = im[ 0, 8, 1, 9, 2,10, 3,11].d0 : ****/\
		"vmovaps 	     (%%rsi),%%zmm0					\n\t		vmovaps 0x040(%%rsi),%%zmm10	\n\t"\
		"vmovaps 	0x400(%%rsi),%%zmm1					\n\t		vmovaps 0x440(%%rsi),%%zmm11	\n\t"\
		"vmovaps 	0x080(%%rsi),%%zmm2					\n\t		vmovaps 0x0c0(%%rsi),%%zmm12	\n\t"\
		"vmovaps 	0x480(%%rsi),%%zmm3					\n\t		vmovaps 0x4c0(%%rsi),%%zmm13	\n\t"\
		"vmovaps 	0x100(%%rsi),%%zmm4					\n\t		vmovaps 0x140(%%rsi),%%zmm14	\n\t"\
		"vmovaps 	0x500(%%rsi),%%zmm5					\n\t		vmovaps 0x540(%%rsi),%%zmm15	\n\t"\
		"vmovaps 	0x180(%%rsi),%%zmm6					\n\t		vmovaps 0x1c0(%%rsi),%%zmm16	\n\t"\
		"vmovaps 	0x580(%%rsi),%%zmm7					\n\t		vmovaps 0x5c0(%%rsi),%%zmm17	\n\t"\
		"\n\t"\
		"vunpcklpd		%%zmm1,%%zmm0,%%zmm8	\n\t	vunpcklpd		%%zmm11,%%zmm10,%%zmm9 		\n\t"\
		"vunpckhpd		%%zmm1,%%zmm0,%%zmm1	\n\t	vunpckhpd		%%zmm11,%%zmm10,%%zmm11		\n\t"\
		"vunpcklpd		%%zmm3,%%zmm2,%%zmm0	\n\t	vunpcklpd		%%zmm13,%%zmm12,%%zmm10		\n\t"\
		"vunpckhpd		%%zmm3,%%zmm2,%%zmm3	\n\t	vunpckhpd		%%zmm13,%%zmm12,%%zmm13		\n\t"\
		"vunpcklpd		%%zmm5,%%zmm4,%%zmm2	\n\t	vunpcklpd		%%zmm15,%%zmm14,%%zmm12		\n\t"\
		"vunpckhpd		%%zmm5,%%zmm4,%%zmm5	\n\t	vunpckhpd		%%zmm15,%%zmm14,%%zmm15		\n\t"\
		"vunpcklpd		%%zmm7,%%zmm6,%%zmm4	\n\t	vunpcklpd		%%zmm17,%%zmm16,%%zmm14		\n\t"\
		"vunpckhpd		%%zmm7,%%zmm6,%%zmm7	\n\t	vunpckhpd		%%zmm17,%%zmm16,%%zmm17		\n\t"\
		"\n\t"\
		"vshuff64x2	$136,%%zmm0,%%zmm8,%%zmm6	\n\t	vshuff64x2	$136,%%zmm10,%%zmm9 ,%%zmm16	\n\t"\
		"vshuff64x2	$221,%%zmm0,%%zmm8,%%zmm0	\n\t	vshuff64x2	$221,%%zmm10,%%zmm9 ,%%zmm10	\n\t"\
		"vshuff64x2	$136,%%zmm3,%%zmm1,%%zmm8	\n\t	vshuff64x2	$136,%%zmm13,%%zmm11,%%zmm9 	\n\t"\
		"vshuff64x2	$221,%%zmm3,%%zmm1,%%zmm3	\n\t	vshuff64x2	$221,%%zmm13,%%zmm11,%%zmm13	\n\t"\
		"vshuff64x2	$136,%%zmm4,%%zmm2,%%zmm1	\n\t	vshuff64x2	$136,%%zmm14,%%zmm12,%%zmm11	\n\t"\
		"vshuff64x2	$221,%%zmm4,%%zmm2,%%zmm4	\n\t	vshuff64x2	$221,%%zmm14,%%zmm12,%%zmm14	\n\t"\
		"vshuff64x2	$136,%%zmm7,%%zmm5,%%zmm2	\n\t	vshuff64x2	$136,%%zmm17,%%zmm15,%%zmm12	\n\t"\
		"vshuff64x2	$221,%%zmm7,%%zmm5,%%zmm7	\n\t	vshuff64x2	$221,%%zmm17,%%zmm15,%%zmm17	\n\t"\
		"\n\t"\
		"vshuff64x2	$136,%%zmm1,%%zmm6,%%zmm5	\n\t	vshuff64x2	$136,%%zmm11,%%zmm16,%%zmm15	\n\t"\
		"vshuff64x2	$221,%%zmm1,%%zmm6,%%zmm1	\n\t	vshuff64x2	$221,%%zmm11,%%zmm16,%%zmm11	\n\t"\
		"vshuff64x2	$136,%%zmm2,%%zmm8,%%zmm6	\n\t	vshuff64x2	$136,%%zmm12,%%zmm9 ,%%zmm16	\n\t"\
		"vshuff64x2	$221,%%zmm2,%%zmm8,%%zmm2	\n\t	vshuff64x2	$221,%%zmm12,%%zmm9 ,%%zmm12	\n\t"\
		"vshuff64x2	$136,%%zmm4,%%zmm0,%%zmm8	\n\t	vshuff64x2	$136,%%zmm14,%%zmm10,%%zmm9 	\n\t"\
		"vshuff64x2	$221,%%zmm4,%%zmm0,%%zmm4	\n\t	vshuff64x2	$221,%%zmm14,%%zmm10,%%zmm14	\n\t"\
		"vshuff64x2	$136,%%zmm7,%%zmm3,%%zmm0	\n\t	vshuff64x2	$136,%%zmm17,%%zmm13,%%zmm10	\n\t"\
		"vshuff64x2	$221,%%zmm7,%%zmm3,%%zmm7	\n\t	vshuff64x2	$221,%%zmm17,%%zmm13,%%zmm17	\n\t"\
		"\n\t"\
		"vmovaps		%%zmm5,(%%rax)			\n\t	vmovaps	%%zmm15,0x40(%%rax)	\n\t"\
		"vmovaps		%%zmm6,(%%rbx)			\n\t	vmovaps	%%zmm16,0x40(%%rbx)	\n\t"\
		"vmovaps		%%zmm8,(%%rcx)			\n\t	vmovaps	%%zmm9 ,0x40(%%rcx)	\n\t"\
		"vmovaps		%%zmm0,(%%rdx)			\n\t	vmovaps	%%zmm10,0x40(%%rdx)	\n\t"\
		"vmovaps		%%zmm1,(%%r10)			\n\t	vmovaps	%%zmm11,0x40(%%r10)	\n\t"\
		"vmovaps		%%zmm2,(%%r11)			\n\t	vmovaps	%%zmm12,0x40(%%r11)	\n\t"\
		"vmovaps		%%zmm4,(%%r12)			\n\t	vmovaps	%%zmm14,0x40(%%r12)	\n\t"\
		"vmovaps		%%zmm7,(%%r13)			\n\t	vmovaps	%%zmm17,0x40(%%r13)	\n\t"\
		"\n\t"\
	/**** a[16-23] = re[ 4,12, 5,13, 6,14, 7,15].d0	a[24-31] = im[ 4,12, 5,13, 6,14, 7,15].d0 : ****/\
	"addq	$0x80,%%rax	\n\t"\
	"addq	$0x80,%%rbx	\n\t"\
	"addq	$0x80,%%rcx	\n\t"\
	"addq	$0x80,%%rdx	\n\t"\
	"addq	$0x80,%%r10	\n\t"\
	"addq	$0x80,%%r11	\n\t"\
	"addq	$0x80,%%r12	\n\t"\
	"addq	$0x80,%%r13	\n\t"\
	"addq	$0x200,%%rsi	\n\t"\
		"vmovaps 	     (%%rsi),%%zmm0					\n\t		vmovaps 0x040(%%rsi),%%zmm10	\n\t"\
		"vmovaps 	0x400(%%rsi),%%zmm1					\n\t		vmovaps 0x440(%%rsi),%%zmm11	\n\t"\
		"vmovaps 	0x080(%%rsi),%%zmm2					\n\t		vmovaps 0x0c0(%%rsi),%%zmm12	\n\t"\
		"vmovaps 	0x480(%%rsi),%%zmm3					\n\t		vmovaps 0x4c0(%%rsi),%%zmm13	\n\t"\
		"vmovaps 	0x100(%%rsi),%%zmm4					\n\t		vmovaps 0x140(%%rsi),%%zmm14	\n\t"\
		"vmovaps 	0x500(%%rsi),%%zmm5					\n\t		vmovaps 0x540(%%rsi),%%zmm15	\n\t"\
		"vmovaps 	0x180(%%rsi),%%zmm6					\n\t		vmovaps 0x1c0(%%rsi),%%zmm16	\n\t"\
		"vmovaps 	0x580(%%rsi),%%zmm7					\n\t		vmovaps 0x5c0(%%rsi),%%zmm17	\n\t"\
		"\n\t"\
		"vunpcklpd		%%zmm1,%%zmm0,%%zmm8	\n\t	vunpcklpd		%%zmm11,%%zmm10,%%zmm9 		\n\t"\
		"vunpckhpd		%%zmm1,%%zmm0,%%zmm1	\n\t	vunpckhpd		%%zmm11,%%zmm10,%%zmm11		\n\t"\
		"vunpcklpd		%%zmm3,%%zmm2,%%zmm0	\n\t	vunpcklpd		%%zmm13,%%zmm12,%%zmm10		\n\t"\
		"vunpckhpd		%%zmm3,%%zmm2,%%zmm3	\n\t	vunpckhpd		%%zmm13,%%zmm12,%%zmm13		\n\t"\
		"vunpcklpd		%%zmm5,%%zmm4,%%zmm2	\n\t	vunpcklpd		%%zmm15,%%zmm14,%%zmm12		\n\t"\
		"vunpckhpd		%%zmm5,%%zmm4,%%zmm5	\n\t	vunpckhpd		%%zmm15,%%zmm14,%%zmm15		\n\t"\
		"vunpcklpd		%%zmm7,%%zmm6,%%zmm4	\n\t	vunpcklpd		%%zmm17,%%zmm16,%%zmm14		\n\t"\
		"vunpckhpd		%%zmm7,%%zmm6,%%zmm7	\n\t	vunpckhpd		%%zmm17,%%zmm16,%%zmm17		\n\t"\
		"\n\t"\
		"vshuff64x2	$136,%%zmm0,%%zmm8,%%zmm6	\n\t	vshuff64x2	$136,%%zmm10,%%zmm9 ,%%zmm16	\n\t"\
		"vshuff64x2	$221,%%zmm0,%%zmm8,%%zmm0	\n\t	vshuff64x2	$221,%%zmm10,%%zmm9 ,%%zmm10	\n\t"\
		"vshuff64x2	$136,%%zmm3,%%zmm1,%%zmm8	\n\t	vshuff64x2	$136,%%zmm13,%%zmm11,%%zmm9 	\n\t"\
		"vshuff64x2	$221,%%zmm3,%%zmm1,%%zmm3	\n\t	vshuff64x2	$221,%%zmm13,%%zmm11,%%zmm13	\n\t"\
		"vshuff64x2	$136,%%zmm4,%%zmm2,%%zmm1	\n\t	vshuff64x2	$136,%%zmm14,%%zmm12,%%zmm11	\n\t"\
		"vshuff64x2	$221,%%zmm4,%%zmm2,%%zmm4	\n\t	vshuff64x2	$221,%%zmm14,%%zmm12,%%zmm14	\n\t"\
		"vshuff64x2	$136,%%zmm7,%%zmm5,%%zmm2	\n\t	vshuff64x2	$136,%%zmm17,%%zmm15,%%zmm12	\n\t"\
		"vshuff64x2	$221,%%zmm7,%%zmm5,%%zmm7	\n\t	vshuff64x2	$221,%%zmm17,%%zmm15,%%zmm17	\n\t"\
		"\n\t"\
		"vshuff64x2	$136,%%zmm1,%%zmm6,%%zmm5	\n\t	vshuff64x2	$136,%%zmm11,%%zmm16,%%zmm15	\n\t"\
		"vshuff64x2	$221,%%zmm1,%%zmm6,%%zmm1	\n\t	vshuff64x2	$221,%%zmm11,%%zmm16,%%zmm11	\n\t"\
		"vshuff64x2	$136,%%zmm2,%%zmm8,%%zmm6	\n\t	vshuff64x2	$136,%%zmm12,%%zmm9 ,%%zmm16	\n\t"\
		"vshuff64x2	$221,%%zmm2,%%zmm8,%%zmm2	\n\t	vshuff64x2	$221,%%zmm12,%%zmm9 ,%%zmm12	\n\t"\
		"vshuff64x2	$136,%%zmm4,%%zmm0,%%zmm8	\n\t	vshuff64x2	$136,%%zmm14,%%zmm10,%%zmm9 	\n\t"\
		"vshuff64x2	$221,%%zmm4,%%zmm0,%%zmm4	\n\t	vshuff64x2	$221,%%zmm14,%%zmm10,%%zmm14	\n\t"\
		"vshuff64x2	$136,%%zmm7,%%zmm3,%%zmm0	\n\t	vshuff64x2	$136,%%zmm17,%%zmm13,%%zmm10	\n\t"\
		"vshuff64x2	$221,%%zmm7,%%zmm3,%%zmm7	\n\t	vshuff64x2	$221,%%zmm17,%%zmm13,%%zmm17	\n\t"\
		"\n\t"\
		"vmovaps		%%zmm5,(%%rax)			\n\t	vmovaps	%%zmm15,0x40(%%rax)	\n\t"\
		"vmovaps		%%zmm6,(%%rbx)			\n\t	vmovaps	%%zmm16,0x40(%%rbx)	\n\t"\
		"vmovaps		%%zmm8,(%%rcx)			\n\t	vmovaps	%%zmm9 ,0x40(%%rcx)	\n\t"\
		"vmovaps		%%zmm0,(%%rdx)			\n\t	vmovaps	%%zmm10,0x40(%%rdx)	\n\t"\
		"vmovaps		%%zmm1,(%%r10)			\n\t	vmovaps	%%zmm11,0x40(%%r10)	\n\t"\
		"vmovaps		%%zmm2,(%%r11)			\n\t	vmovaps	%%zmm12,0x40(%%r11)	\n\t"\
		"vmovaps		%%zmm4,(%%r12)			\n\t	vmovaps	%%zmm14,0x40(%%r12)	\n\t"\
		"vmovaps		%%zmm7,(%%r13)			\n\t	vmovaps	%%zmm17,0x40(%%r13)	\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__add1] "m" (Xadd1)\
		 ,[__add2] "m" (Xadd2)\
		 ,[__add3] "m" (Xadd3)\
		 ,[__add4] "m" (Xadd4)\
		 ,[__add5] "m" (Xadd5)\
		 ,[__add6] "m" (Xadd6)\
		 ,[__add7] "m" (Xadd7)\
		 ,[__r1] "m" (Xr1)\
		 ,[__r9] "m" (Xr9)\
		 ,[__r17] "m" (Xr17)\
		 ,[__r25] "m" (Xr25)\
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__cc0] "m" (Xcc0)\
		 ,[__c1] "m" (Xc1)\
		 ,[__c2] "m" (Xc2)\
		 ,[__c3] "m" (Xc3)\
		 ,[__c4] "m" (Xc4)\
		 ,[__c5] "m" (Xc5)\
		 ,[__c6] "m" (Xc6)\
		 ,[__c7] "m" (Xc7)\
		 ,[__c8] "m" (Xc8)\
		 ,[__c9] "m" (Xc9)\
		 ,[__c10] "m" (Xc10)\
		 ,[__c11] "m" (Xc11)\
		 ,[__c12] "m" (Xc12)\
		 ,[__c13] "m" (Xc13)\
		 ,[__c14] "m" (Xc14)\
		 ,[__c15] "m" (Xc15)\
	 ,[__pfetch_dist] "m" (Xpfetch_dist)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r10","r11","r12","r13","r14","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15","xmm16","xmm17", "xmm29","xmm30","xmm31"	/* Clobbered registers */\
	);\
	}

  #endif	// (IMCI512 or AVX512?) toggle

#elif defined(USE_AVX2)	// FMA-based versions of selected macros in this file for Intel AVX2/FMA3

	// Since the add0/1 block addresses advance in opposite directions, our prefetch scheme here is like so:
	// In fwd DIF prefetch  ahead 8 cache lines' worth of data w.r.to add0;
	// In inv DIT prefetch behind 8 cache lines' worth of data w.r.to add1;
	#define SSE2_RADIX16_WRAPPER_DIF(Xadd0,Xadd1,Xadd2,Xadd3,Xr1,Xr9,Xr17,Xr25,Xisrt2,Xcc0,Xc1,Xc2,Xc3,Xc4,Xc5,Xc6,Xc7,Xc8,Xc9,Xc10,Xc11,Xc12,Xc13,Xc14,Xc15,Xpfetch_dist)\
	{\
	__asm__ volatile (\
	"movq	%[__add0],%%rax						\n\t"\
	"movslq	%[__pfetch_dist],%%r14	\n\t"\
	"leaq	(%%rax,%%r14,8),%%r14	\n\t"	/* Block 0 [base-address + data-fetch-ahead index] */\
	"prefetcht1	(%%r14)\n\t"\
	/*************************************************************/\
	/*                  1st set of inputs:                       */\
	/*************************************************************/\
		"movq	%[__r1] ,%%rsi	\n\t"\
		"movq	%[__add1],%%rbx\n\t"\
		"movq	%[__add2],%%rcx\n\t"\
		"movq	%[__add3],%%rdx\n\t"\
	/**** Start with 4-way interleaving: ****/\
	"/* a[j+p0]: Inputs from [add0,add0+0x100,add1,add1+0x100]+0x0. Outputs into r1 +0/1, 8/9, 16/17, 24/25: */	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps	     (%%rax),%%ymm1						\n\t		vmovaps	0x020(%%rax),%%ymm3								\n\t"\
		"vmovaps	     (%%rbx),%%ymm5						\n\t		vmovaps	0x020(%%rbx),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps	     (%%rcx),%%ymm6						\n\t		vmovaps	0x020(%%rcx),%%ymm14							\n\t"\
		"vmovaps	     (%%rdx),%%ymm5						\n\t		vmovaps	0x020(%%rdx),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm6 ,%%ymm0				\n\t		vshufpd	$15 ,%%ymm13,%%ymm14 ,%%ymm2					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm6 ,%%ymm6				\n\t		vshufpd	$0  ,%%ymm13,%%ymm14 ,%%ymm14					\n\t"\
		"vperm2f128 $32 ,%%ymm0 ,%%ymm7 ,%%ymm5				\n\t		vperm2f128 $32 ,%%ymm2 ,%%ymm15,%%ymm13		/* outB	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm0 ,%%ymm7 ,%%ymm7				\n\t		vperm2f128 $49 ,%%ymm2 ,%%ymm15,%%ymm15		/* outD	*/	\n\t"\
		"vperm2f128 $32 ,%%ymm6 ,%%ymm1 ,%%ymm0				\n\t		vperm2f128 $32 ,%%ymm14,%%ymm3 ,%%ymm2 		/* outA	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm6 ,%%ymm1 ,%%ymm1				\n\t		vperm2f128 $49 ,%%ymm14,%%ymm3 ,%%ymm3 		/* outC	*/	\n\t"\
		"vmovaps 	%%ymm5 ,0x200(%%rsi)					\n\t		vmovaps %%ymm13,0x220(%%rsi)				/* outB	*/	\n\t"\
		"vmovaps 	%%ymm7 ,0x300(%%rsi)					\n\t		vmovaps %%ymm15,0x320(%%rsi)				/* outD	*/	\n\t"\
		"vmovaps 	%%ymm0 ,     (%%rsi)					\n\t		vmovaps %%ymm2 ,0x020(%%rsi)				/* outA	*/	\n\t"\
		"vmovaps 	%%ymm1 ,0x100(%%rsi)					\n\t		vmovaps %%ymm3 ,0x120(%%rsi)				/* outC	*/	\n\t"\
	"/* a[j+p2]: Inputs from [add0,add0+0x100,add1,add1+0x100]+0x40. Outputs into r3 +0/1, 8/9, 16/17, 24/25: */		\n\t"\
		"addq	$0x40,%%rax	\n\t"\
		"addq	$0x40,%%rbx	\n\t"\
		"addq	$0x40,%%rcx	\n\t"\
		"addq	$0x40,%%rdx	\n\t"\
		"addq	$0x40,%%rsi	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps	     (%%rax),%%ymm1						\n\t		vmovaps	0x020(%%rax),%%ymm3								\n\t"\
		"vmovaps	     (%%rbx),%%ymm5						\n\t		vmovaps	0x020(%%rbx),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps	     (%%rcx),%%ymm6						\n\t		vmovaps	0x020(%%rcx),%%ymm14							\n\t"\
		"vmovaps	     (%%rdx),%%ymm5						\n\t		vmovaps	0x020(%%rdx),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm6 ,%%ymm0				\n\t		vshufpd	$15 ,%%ymm13,%%ymm14 ,%%ymm2					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm6 ,%%ymm6				\n\t		vshufpd	$0  ,%%ymm13,%%ymm14 ,%%ymm14					\n\t"\
		"vperm2f128 $32 ,%%ymm0 ,%%ymm7 ,%%ymm5				\n\t		vperm2f128 $32 ,%%ymm2 ,%%ymm15,%%ymm13		/* outB	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm0 ,%%ymm7 ,%%ymm7				\n\t		vperm2f128 $49 ,%%ymm2 ,%%ymm15,%%ymm15		/* outD	*/	\n\t"\
		"vperm2f128 $32 ,%%ymm6 ,%%ymm1 ,%%ymm0				\n\t		vperm2f128 $32 ,%%ymm14,%%ymm3 ,%%ymm2 		/* outA	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm6 ,%%ymm1 ,%%ymm1				\n\t		vperm2f128 $49 ,%%ymm14,%%ymm3 ,%%ymm3 		/* outC	*/	\n\t"\
		"vmovaps 	%%ymm5 ,0x200(%%rsi)					\n\t		vmovaps %%ymm13,0x220(%%rsi)				/* outB	*/	\n\t"\
		"vmovaps 	%%ymm7 ,0x300(%%rsi)					\n\t		vmovaps %%ymm15,0x320(%%rsi)				/* outD	*/	\n\t"\
		"vmovaps 	%%ymm0 ,     (%%rsi)					\n\t		vmovaps %%ymm2 ,0x020(%%rsi)				/* outA	*/	\n\t"\
		"vmovaps 	%%ymm1 ,0x100(%%rsi)					\n\t		vmovaps %%ymm3 ,0x120(%%rsi)				/* outC	*/	\n\t"\
	"prefetcht1	0x40(%%r14)\n\t"\
	"/* a[j+p4]: Inputs from [add0,add0+0x100,add1,add1+0x100]+0x80. Outputs into r5 +0/1, 8/9, 16/17, 24/25: */		\n\t"\
		"addq	$0x40,%%rax	\n\t"\
		"addq	$0x40,%%rbx	\n\t"\
		"addq	$0x40,%%rcx	\n\t"\
		"addq	$0x40,%%rdx	\n\t"\
		"addq	$0x40,%%rsi	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps	     (%%rax),%%ymm1						\n\t		vmovaps	0x020(%%rax),%%ymm3								\n\t"\
		"vmovaps	     (%%rbx),%%ymm5						\n\t		vmovaps	0x020(%%rbx),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps	     (%%rcx),%%ymm6						\n\t		vmovaps	0x020(%%rcx),%%ymm14							\n\t"\
		"vmovaps	     (%%rdx),%%ymm5						\n\t		vmovaps	0x020(%%rdx),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm6 ,%%ymm0				\n\t		vshufpd	$15 ,%%ymm13,%%ymm14 ,%%ymm2					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm6 ,%%ymm6				\n\t		vshufpd	$0  ,%%ymm13,%%ymm14 ,%%ymm14					\n\t"\
		"vperm2f128 $32 ,%%ymm0 ,%%ymm7 ,%%ymm5				\n\t		vperm2f128 $32 ,%%ymm2 ,%%ymm15,%%ymm13		/* outB	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm0 ,%%ymm7 ,%%ymm7				\n\t		vperm2f128 $49 ,%%ymm2 ,%%ymm15,%%ymm15		/* outD	*/	\n\t"\
		"vperm2f128 $32 ,%%ymm6 ,%%ymm1 ,%%ymm0				\n\t		vperm2f128 $32 ,%%ymm14,%%ymm3 ,%%ymm2 		/* outA	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm6 ,%%ymm1 ,%%ymm1				\n\t		vperm2f128 $49 ,%%ymm14,%%ymm3 ,%%ymm3 		/* outC	*/	\n\t"\
		"vmovaps 	%%ymm5 ,0x200(%%rsi)					\n\t		vmovaps %%ymm13,0x220(%%rsi)				/* outB	*/	\n\t"\
		"vmovaps 	%%ymm7 ,0x300(%%rsi)					\n\t		vmovaps %%ymm15,0x320(%%rsi)				/* outD	*/	\n\t"\
		"vmovaps 	%%ymm0 ,     (%%rsi)					\n\t		vmovaps %%ymm2 ,0x020(%%rsi)				/* outA	*/	\n\t"\
		"vmovaps 	%%ymm1 ,0x100(%%rsi)					\n\t		vmovaps %%ymm3 ,0x120(%%rsi)				/* outC	*/	\n\t"\
	"/* a[j+p6]: Inputs from [add0,add0+0x100,add1,add1+0x100]+0xc0. Outputs into r7 +0/1, 8/9, 16/17, 24/25: */		\n\t"\
		"addq	$0x40,%%rax	\n\t"\
		"addq	$0x40,%%rbx	\n\t"\
		"addq	$0x40,%%rcx	\n\t"\
		"addq	$0x40,%%rdx	\n\t"\
		"addq	$0x40,%%rsi	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps	     (%%rax),%%ymm1						\n\t		vmovaps	0x020(%%rax),%%ymm3								\n\t"\
		"vmovaps	     (%%rbx),%%ymm5						\n\t		vmovaps	0x020(%%rbx),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps	     (%%rcx),%%ymm6						\n\t		vmovaps	0x020(%%rcx),%%ymm14							\n\t"\
		"vmovaps	     (%%rdx),%%ymm5						\n\t		vmovaps	0x020(%%rdx),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm6 ,%%ymm0				\n\t		vshufpd	$15 ,%%ymm13,%%ymm14 ,%%ymm2					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm6 ,%%ymm6				\n\t		vshufpd	$0  ,%%ymm13,%%ymm14 ,%%ymm14					\n\t"\
		"vperm2f128 $32 ,%%ymm0 ,%%ymm7 ,%%ymm5				\n\t		vperm2f128 $32 ,%%ymm2 ,%%ymm15,%%ymm13		/* outB	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm0 ,%%ymm7 ,%%ymm7				\n\t		vperm2f128 $49 ,%%ymm2 ,%%ymm15,%%ymm15		/* outD	*/	\n\t"\
		"vperm2f128 $32 ,%%ymm6 ,%%ymm1 ,%%ymm0				\n\t		vperm2f128 $32 ,%%ymm14,%%ymm3 ,%%ymm2 		/* outA	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm6 ,%%ymm1 ,%%ymm1				\n\t		vperm2f128 $49 ,%%ymm14,%%ymm3 ,%%ymm3 		/* outC	*/	\n\t"\
		"vmovaps 	%%ymm5 ,0x200(%%rsi)					\n\t		vmovaps %%ymm13,0x220(%%rsi)				/* outB	*/	\n\t"\
		"vmovaps 	%%ymm7 ,0x300(%%rsi)					\n\t		vmovaps %%ymm15,0x320(%%rsi)				/* outD	*/	\n\t"\
		"vmovaps 	%%ymm0 ,     (%%rsi)					\n\t		vmovaps %%ymm2 ,0x020(%%rsi)				/* outA	*/	\n\t"\
		"vmovaps 	%%ymm1 ,0x100(%%rsi)					\n\t		vmovaps %%ymm3 ,0x120(%%rsi)				/* outC	*/	\n\t"\
	/*****************/\
	/* Radix-16 DIF: */\
	/*****************/\
	"prefetcht1	0x80(%%r14)\n\t"\
	"movq	%[__add0],%%rax	\n\t"/* Use for FMA-related spills */\
	"movq	%[__r1] ,%%rcx	\n\t"\
	"leaq	0x880(%%rcx),%%rsi	\n\t"/* two */\
	/*...Block 1: */													/*...Block 2: */\
		"leaq	0x560(%%rcx),%%rdi	/* c2 */\n\t"\
		"leaq	0x4e0(%%rcx),%%rdx	/* c4 */\n\t"\
		"vmovaps	0x040(%%rcx),%%ymm0	/* ymm0 <-     a[jt+p4] */			\n\t		vmovaps		0x100(%%rcx),%%ymm8	/* ymm10 <-     a[jt+p2] */			\n\t"\
		"vmovaps	0x060(%%rcx),%%ymm1	/* ymm1 <-     a[jp+p4] */			\n\t		vmovaps		0x120(%%rcx),%%ymm9	/* ymm11 <-     a[jp+p2] */			\n\t"\
		"vmovaps	%%ymm0		,%%ymm2	/* ymm2 <- cpy a[jt+p4] */			\n\t		vmovaps		%%ymm8 	,%%ymm10	/* ymm10 <- cpy a[jt+p2] */			\n\t"\
		"vmovaps	%%ymm1		,%%ymm3	/* ymm3 <- cpy a[jp+p4] */			\n\t		vmovaps		%%ymm9 	,%%ymm11	/* ymm11 <- cpy a[jp+p2] */			\n\t"\
		/***************************************************************************/\
		/*** From hereon, things are identical to the code in radix16_dif_pass: ****/\
		/***************************************************************************/\
		"vmovaps	0x0c0(%%rcx),%%ymm4			/* ymm4 <-     a[jt+p12] */	\n\t		vmovaps		0x180(%%rcx),%%ymm12			/* ymm12 <-     a[jt+p10] */\n\t"\
		"vmovaps	0x0e0(%%rcx),%%ymm5			/* ymm5 <-     a[jp+p12] */	\n\t		vmovaps		0x1a0(%%rcx),%%ymm13			/* ymm13 <-     a[jp+p10] */\n\t"\
		"vmovaps	%%ymm4		,%%ymm6			/* ymm6 <- cpy a[jt+p12] */	\n\t		vmovaps			%%ymm12	,%%ymm14			/* ymm14 <- cpy a[jt+p10] */\n\t"\
		"vmovaps	%%ymm5		,%%ymm7			/* ymm7 <- cpy a[jp+p12] */	\n\t		vmovaps			%%ymm13	,%%ymm15			/* ymm15 <- cpy a[jp+p10] */\n\t"\
		"vmulpd		     (%%rdx),%%ymm0,%%ymm0		/* A: a[jt+ p4]* c4 */	\n\t		vmulpd		     (%%rdi),%%ymm8 ,%%ymm8 	\n\t"/* E: a[jt+ p2]* c2 */\
		"vmulpd		     (%%rdx),%%ymm1,%%ymm1		/* B: a[jp+ p4]* c4 */	\n\t		vmulpd		     (%%rdi),%%ymm9 ,%%ymm9 	\n\t"/* F: a[jp+ p2]* c2 */\
		"vmulpd		0x040(%%rdx),%%ymm4,%%ymm4		/* C: a[jt+p12]*c12 */	\n\t		vmulpd		0x040(%%rdi),%%ymm12,%%ymm12	\n\t"/* G: a[jt+p10]*c10 */\
		"vmulpd		0x040(%%rdx),%%ymm5,%%ymm5		/* D: a[jp+p12]*c12 */	\n\t		vmulpd		0x040(%%rdi),%%ymm13,%%ymm13	\n\t"/* H: a[jp+p10]*c10 */\
	"vfnmadd231pd	0x020(%%rdx),%%ymm3,%%ymm0		/* A -= a[jp+ p4]* s4 */\n\t	vfnmadd231pd	0x020(%%rdi),%%ymm11,%%ymm8 	\n\t"/* E -= a[jp+ p2]* s2 */\
	" vfmadd231pd	0x020(%%rdx),%%ymm2,%%ymm1		/* B += a[jt+ p4]* s4 */\n\t	 vfmadd231pd	0x020(%%rdi),%%ymm10,%%ymm9 	\n\t"/* F += a[jt+ p2]* s2 */\
	"vfnmadd231pd	0x060(%%rdx),%%ymm7,%%ymm4		/* C -= a[jp+p12]*s12 */\n\t	vfnmadd231pd	0x060(%%rdi),%%ymm15,%%ymm12	\n\t"/* G -= a[jp+p10]*s10 */\
	" vfmadd231pd	0x060(%%rdx),%%ymm6,%%ymm5		/* D += a[jt+p12]*s12 */\n\t	 vfmadd231pd	0x060(%%rdi),%%ymm14,%%ymm13	\n\t"/* H += a[jt+p10]*s10 */\
		"vmovaps	%%ymm1		,%%ymm3		/* ymm3 <- cpy t6 */			\n\t		vmovaps		%%ymm9 		,%%ymm11		/* ymm11 <- cpy t10*/		\n\t"\
		"vmovaps	%%ymm0		,%%ymm2		/* ymm2 <- cpy t5 */			\n\t		vmovaps		%%ymm8 		,%%ymm10		/* ymm10 <- cpy t9 */		\n\t"\
		"vaddpd		%%ymm4		,%%ymm0,%%ymm0		/* ~t5 <- t5 +rt */		\n\t		vaddpd		%%ymm12		,%%ymm8 ,%%ymm8  	/* ~t13<- t13+rt */						\n\t"\
		"vaddpd		%%ymm5		,%%ymm1,%%ymm1		/* ~t6 <- t6 +it */		\n\t		vaddpd		%%ymm13		,%%ymm9 ,%%ymm9  	/* ~t14<- t14+it */						\n\t"\
		"vsubpd		%%ymm4		,%%ymm2,%%ymm2		/* ~t7 <- t5 -rt */		\n\t		vsubpd		%%ymm12		,%%ymm10,%%ymm10 	/* ~t15<- t13-rt */						\n\t"\
		"vsubpd		%%ymm5		,%%ymm3,%%ymm3		/* ~t8 <- t6 -it */		\n\t		vsubpd		%%ymm13		,%%ymm11,%%ymm11 	/* ~t16<- t14-it	ymm12,13 free */	\n\t"\
	"prefetcht1	0xc0(%%r14)\n\t"\
		"\n\t"\
		"/* Now do the p0,8 combo: */\n\t									/* Do the p6,14 combo - do p14 first so registers come out in same order as for p2,10 */\n\t"\
	"leaq	0x4a0(%%rcx),%%rdx	/* c8 */								\n\t		leaq	0x620(%%rcx),%%rdi	/* c14 */\n\t"\
		"vmovaps	0x080(%%rcx)	,%%ymm4		/* a[jt+p8 ] */				\n\t		vmovaps		0x1c0(%%rcx),%%ymm12		/* a[jt+p14] */				\n\t"\
		"vmovaps	0x0a0(%%rcx)	,%%ymm5		/* a[jp+p8 ] */				\n\t		vmovaps		0x1e0(%%rcx),%%ymm13		/* a[jp+p14] */				\n\t"\
		"vmovaps	%%ymm4		,%%ymm6	/* ymm6 <- cpy a[jt+p8] */			\n\t		vmovaps			%%ymm12	,%%ymm14		/* ymm14 <- cpy a[jt+p14] */\n\t"\
		"vmovaps	%%ymm5		,%%ymm7	/* ymm7 <- cpy a[jp+p8] */			\n\t		vmovaps			%%ymm13	,%%ymm15		/* ymm15 <- cpy a[jp+p14] */\n\t"\
		"vmulpd		     (%%rdx),%%ymm4,%%ymm4		/* a[jt+p8]*c8 */		\n\t		vmulpd		    (%%rdi)	,%%ymm12,%%ymm12		/* a[jt+p14]*c14 */			\n\t"\
		"vmulpd		     (%%rdx),%%ymm5,%%ymm5		/* a[jp+p8]*c8 */		\n\t		vmulpd		    (%%rdi)	,%%ymm13,%%ymm13		/* a[jp+p14]*c14 */			\n\t"\
	"vfnmadd231pd	0x020(%%rdx),%%ymm7,%%ymm4		/* a[jp+p8]*s8 */		\n\t	vfnmadd231pd	0x020(%%rdi),%%ymm15,%%ymm12		/* a[jp+p14]*s14 */			\n\t"\
	" vfmadd231pd	0x020(%%rdx),%%ymm6,%%ymm5		/* a[jt+p8]*s8 */		\n\t	 vfmadd231pd	0x020(%%rdi),%%ymm14,%%ymm13		/* a[jt+p14]*s14 */			\n\t"\
		"																	\n\t		vmovaps		%%ymm13		,0x1e0(%%rcx)	/* Store it in t16*/		\n\t"\
	"/* Real parts:*/														\n\t		vmovaps		%%ymm12		,0x1c0(%%rcx)	/* Store rt in t15*/		\n\t"\
		"vmovaps		 (%%rcx),%%ymm6		/* a[jt    ] */					\n\t		subq	$0x040,%%rdi	/* c6  */	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm7		/* a[jp    ] */					\n\t	/* Real parts: */\n\t"\
		"vsubpd		%%ymm4 ,%%ymm6,%%ymm6	/* ~t3 <- t1 -rt */				\n\t		vmovaps		0x140(%%rcx),%%ymm12		/* a[jt+p6 ] */				\n\t"\
		"vsubpd		%%ymm5 ,%%ymm7,%%ymm7	/* ~t4 <- t2 -it */				\n\t		vmovaps		0x160(%%rcx),%%ymm13		/* a[jp+p6 ] */				\n\t"\
	"vmovaps	(%%rsi),%%ymm15	\n\t"/*two */\
	" vfmadd132pd	%%ymm15,%%ymm6,%%ymm4	/* ~t1 <- t1 +rt */				\n\t		vmovaps			%%ymm12	,%%ymm14		/* ymm14 <- cpy a[jt+p6] */		\n\t"\
	" vfmadd132pd	%%ymm15,%%ymm7,%%ymm5	/* ~t2 <- t2 +it */				\n\t		vmovaps			%%ymm13	,%%ymm15		/* ymm15 <- cpy a[jp+p6] */		\n\t"\
		/* Finish radix-4 butterfly and store results into temporary-array slots: */\
		"vsubpd		%%ymm0		,%%ymm4,%%ymm4	/*~t5 =t1 -t5 */			\n\t		vmulpd		    (%%rdi)	,%%ymm12,%%ymm12		/* a[jt+p6]*c6 */			\n\t"\
		"vsubpd		%%ymm1		,%%ymm5,%%ymm5	/*~t6 =t2 -t6 */			\n\t		vmulpd		    (%%rdi)	,%%ymm13,%%ymm13		/* a[jp+p6]*c6 */			\n\t"\
		"vsubpd		%%ymm3		,%%ymm6,%%ymm6	/*~t3 =t3 -t8 */			\n\t	 vfmadd231pd	0x020(%%rdi),%%ymm14,%%ymm13		/* a[jt+p6]*s6 */			\n\t"\
		"vsubpd		%%ymm2		,%%ymm7,%%ymm7	/*~t8 =t4 -t7 */			\n\t	vfnmadd231pd	0x020(%%rdi),%%ymm15,%%ymm12		/* a[jp+p6]*s6 */			\n\t"\
		"vmovaps	%%ymm4		,0x080(%%rcx)	/* a[jt+p8 ] <- ~t5 */		\n\t		vmovaps		%%ymm13		,%%ymm15		/* ymm15 <- cpy t14*/			\n\t"\
		"vmovaps	%%ymm5		,0x0a0(%%rcx)	/* a[jp+p8 ] <- ~t6 */		\n\t		vmovaps		%%ymm12		,%%ymm14		/* ymm14 <- cpy t13*/			\n\t"\
	"vfmadd132pd	(%%rsi),%%ymm4,%%ymm0		/* 2*t5 */				\n\t			vsubpd		0x1c0(%%rcx),%%ymm12,%%ymm12		/* ~t15<- t13-rt */			\n\t"\
	"vfmadd132pd	(%%rsi),%%ymm5,%%ymm1		/* 2*t6 */				\n\t			vsubpd		0x1e0(%%rcx),%%ymm13,%%ymm13		/* ~t16<- t14-it */			\n\t"\
	"vfmadd132pd	(%%rsi),%%ymm6,%%ymm3	/*~t7 =t3 +t8 */			\n\t			vaddpd		0x1c0(%%rcx),%%ymm14,%%ymm14		/* ~t13<- t13+rt */			\n\t"\
	"vfmadd132pd	(%%rsi),%%ymm7,%%ymm2	/*~t4 =t4 +t7 */			\n\t			vaddpd		0x1e0(%%rcx),%%ymm15,%%ymm15		/* ~t14<- t14+it */			\n\t"\
		"vmovaps	%%ymm0		,     (%%rcx)	/* a[jt    ] <- ~t1 */		\n\t		vsubpd		%%ymm14		,%%ymm8,%%ymm8 	/*~t13*/						\n\t"\
		"vmovaps	%%ymm1		,0x020(%%rcx)	/* a[jp    ] <- ~t2 */		\n\t		vsubpd		%%ymm15		,%%ymm9,%%ymm9 	/*~t14*/						\n\t"\
		"vmovaps	%%ymm6		,0x040(%%rcx)	/* a[jt+p4 ] <- ~t3 */		\n\t		vsubpd		%%ymm13		,%%ymm10,%%ymm10	/*~t11*/				\n\t"\
		"vmovaps	%%ymm7		,0x0e0(%%rcx)	/* a[jp+p12] <- ~t8 */		\n\t		vsubpd		%%ymm12		,%%ymm11,%%ymm11	/*~t16*/				\n\t"\
		"vmovaps	%%ymm3		,0x0c0(%%rcx)	/* a[jt+p12] <- ~t7 */		\n\t		vmovaps		%%ymm8 		,0x180(%%rcx)	/* a[jt+p8 ] <- ~t13*/		\n\t"\
		"vmovaps	%%ymm2		,0x060(%%rcx)	/* a[jp+p4 ] <- ~t4 */		\n\t		vmovaps		%%ymm9 		,0x1a0(%%rcx)	/* a[jp+p8 ] <- ~t14*/		\n\t"\
	"vmovaps	(%%rsi),%%ymm1	\n\t"/*two */\
		"																			vfmadd132pd	%%ymm1,%%ymm8 ,%%ymm14	/* 2*t13*/						\n\t"\
		"																			vfmadd132pd	%%ymm1,%%ymm9 ,%%ymm15	/* 2*t14*/						\n\t"\
		"																			vfmadd132pd	%%ymm1,%%ymm10,%%ymm13			/*~t15*/				\n\t"\
		"																			vfmadd132pd	%%ymm1,%%ymm11,%%ymm12			/*~t12*/				\n\t"\
		"																				vmovaps		%%ymm14		,0x100(%%rcx)	/* a[jt    ] <- ~t9 */		\n\t"\
		"																				vmovaps		%%ymm15		,0x120(%%rcx)	/* a[jp    ] <- ~t10*/		\n\t"\
		"																				vmovaps		%%ymm10		,0x140(%%rcx)	/* a[jt+p4 ] <- ~t11*/		\n\t"\
		"																				vmovaps		%%ymm11		,0x1e0(%%rcx)	/* a[jp+p12] <- ~t16*/		\n\t"\
		"																				vmovaps		%%ymm13		,0x1c0(%%rcx)	/* a[jt+p12] <- ~t15*/		\n\t"\
		"																				vmovaps		%%ymm12		,0x160(%%rcx)	/* a[jp+p4 ] <- ~t12*/		\n\t"\
	"prefetcht1	0x100(%%r14)\n\t"\
		"\n\t"\
	/*...Block 3: */															/*...Block 4: */\
	"/*	SSE2_RADIX4_DIF_4TWIDDLE(r17,r21,r19,r23, r17,c1): */					/* SSE2_RADIX4_DIF_4TWIDDLE(r25,r29,r27,r31, r25,c3): */\n\t"\
		"/* Do the p0,p8 combo: */											\n\t	/* Do the p0,p8 combo: */		\n\t"\
	"leaq	0x660(%%rcx),%%rbx	\n\t"/* c1 */	/* All __r and __c pointers incr by +0x100 in rcol w.r.to lcol: */\
	"addq	$0x200,%%rcx		\n\t"/* r17 */											/* c3, r25 */\
		"vmovaps		 (%%rcx),%%ymm0		/* a[jt   ] */					\n\t		vmovaps		0x100(%%rcx),%%ymm8 		/* a[jt    ] */\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm1		/* a[jp   ] */					\n\t		vmovaps		0x120(%%rcx),%%ymm9 		/* a[jp    ] */\n\t"\
		"vmovaps	0x080(%%rcx),%%ymm4			/* ymm4 <-     a[jt+p12] */	\n\t		vmovaps		0x180(%%rcx),%%ymm12			/* ymm12 <-     a[jt+p10] */\n\t"\
		"vmovaps	0x0a0(%%rcx),%%ymm5			/* ymm5 <-     a[jp+p12] */	\n\t		vmovaps		0x1a0(%%rcx),%%ymm13			/* ymm13 <-     a[jp+p10] */\n\t"\
		"vmovaps	%%ymm0		,%%ymm2		/* ymm2 <- cpy a[jt   ] */		\n\t		vmovaps		%%ymm8 		,%%ymm10		/* ymm10 <- cpy a[jt   ] */\n\t"\
		"vmovaps	%%ymm1		,%%ymm3		/* ymm3 <- cpy a[jp   ] */		\n\t		vmovaps		%%ymm9 		,%%ymm11		/* ymm11 <- cpy a[jp   ] */\n\t"\
		"vmovaps	%%ymm4		,%%ymm6			/* ymm6 <- cpy a[jt+p12] */	\n\t		vmovaps			%%ymm12	,%%ymm14			/* ymm14 <- cpy a[jt+p10] */\n\t"\
		"vmovaps	%%ymm5		,%%ymm7			/* ymm7 <- cpy a[jp+p12] */	\n\t		vmovaps			%%ymm13	,%%ymm15			/* ymm15 <- cpy a[jp+p10] */\n\t"\
		"vmulpd		     (%%rbx),%%ymm0,%%ymm0		/* A: a[jt+ p4]* c4 */	\n\t		vmulpd		0x100(%%rbx),%%ymm8 ,%%ymm8 	\n\t"/* E: a[jt+ p2]* c2 */\
		"vmulpd		     (%%rbx),%%ymm1,%%ymm1		/* B: a[jp+ p4]* c4 */	\n\t		vmulpd		0x100(%%rbx),%%ymm9 ,%%ymm9 	\n\t"/* F: a[jp+ p2]* c2 */\
		"vmulpd		0x040(%%rbx),%%ymm4,%%ymm4		/* C: a[jt+p12]*c12 */	\n\t		vmulpd		0x140(%%rbx),%%ymm12,%%ymm12	\n\t"/* G: a[jt+p10]*c10 */\
		"vmulpd		0x040(%%rbx),%%ymm5,%%ymm5		/* D: a[jp+p12]*c12 */	\n\t		vmulpd		0x140(%%rbx),%%ymm13,%%ymm13	\n\t"/* H: a[jp+p10]*c10 */\
	"vfnmadd231pd	0x020(%%rbx),%%ymm3,%%ymm0		/* A -= a[jp+ p4]* s4 */\n\t	vfnmadd231pd	0x120(%%rbx),%%ymm11,%%ymm8 	\n\t"/* E -= a[jp+ p2]* s2 */\
	" vfmadd231pd	0x020(%%rbx),%%ymm2,%%ymm1		/* B += a[jt+ p4]* s4 */\n\t	 vfmadd231pd	0x120(%%rbx),%%ymm10,%%ymm9 	\n\t"/* F += a[jt+ p2]* s2 */\
	"vfnmadd231pd	0x060(%%rbx),%%ymm7,%%ymm4		/* C -= a[jp+p12]*s12 */\n\t	vfnmadd231pd	0x160(%%rbx),%%ymm15,%%ymm12	\n\t"/* G -= a[jp+p10]*s10 */\
	" vfmadd231pd	0x060(%%rbx),%%ymm6,%%ymm5		/* D += a[jt+p12]*s12 */\n\t	 vfmadd231pd	0x160(%%rbx),%%ymm14,%%ymm13	\n\t"/* H += a[jt+p10]*s10 */\
		"vmovaps	%%ymm0		,%%ymm2		/* ymm2 <- cpy t5 */			\n\t		vmovaps		%%ymm8 		,%%ymm10		/* ymm10 <- cpy t9 */		\n\t"\
		"vmovaps	%%ymm1		,%%ymm3		/* ymm3 <- cpy t6 */			\n\t		vmovaps		%%ymm9 		,%%ymm11		/* ymm11 <- cpy t10*/		\n\t"\
		"vaddpd		%%ymm4	,%%ymm0,%%ymm0		/* ~t1 <- t1 +rt */			\n\t		vaddpd		%%ymm12	,%%ymm8 ,%%ymm8 		/* ~t1 <- t1 +rt */\n\t"\
		"vaddpd		%%ymm5	,%%ymm1,%%ymm1		/* ~t2 <- t2 +it */			\n\t		vaddpd		%%ymm13	,%%ymm9 ,%%ymm9 		/* ~t2 <- t2 +it */\n\t"\
		"vsubpd		%%ymm4	,%%ymm2,%%ymm2		/* ~t3 <- t1 -rt */			\n\t		vsubpd		%%ymm12	,%%ymm10,%%ymm10		/* ~t3 <- t1 -rt */\n\t"\
		"vsubpd		%%ymm5	,%%ymm3,%%ymm3	/* ~t4 <- t2 -it ymm4,5 free*/	\n\t		vsubpd		%%ymm13	,%%ymm11,%%ymm11		/* ~t4 <- t2 -it	ymm12,5 free */\n\t"\
		"\n\t"\
	"/* Do the p4,12 combo: */												\n\t	/* Do the p4,12 combo: */\n\t"\
	"addq	$0x0c0 ,%%rbx	\n\t"/* c13 */											/* c15 */\
		"vmovaps	0x0c0(%%rcx),%%ymm6		/* a[jt+p12] */					\n\t		vmovaps		0x1c0(%%rcx),%%ymm14		/* a[jt+p12] */\n\t"\
		"vmovaps	0x0e0(%%rcx),%%ymm7		/* a[jp+p12] */					\n\t		vmovaps		0x1e0(%%rcx),%%ymm15		/* a[jp+p12] */\n\t"\
		"vmovaps	%%ymm6		,%%ymm4		/* ymm4 <- cpy a[jt+p12] */		\n\t		vmovaps		%%ymm14		,%%ymm12		/* ymm12 <- cpy a[jt+p12] */\n\t"\
		"vmovaps	%%ymm7		,%%ymm5		/* ymm5 <- cpy a[jp+p12] */		\n\t		vmovaps		%%ymm15		,%%ymm13		/* ymm13 <- cpy a[jp+p12] */\n\t"\
		"vmulpd		(%%rbx)		,%%ymm4,%%ymm4		/* a[jt+p12]*c12 */		\n\t		vmulpd		0x100(%%rbx),%%ymm12,%%ymm12		/* a[jt+p12]*c12 */\n\t"\
		"vmulpd		(%%rbx)		,%%ymm5,%%ymm5		/* a[jp+p12]*c12 */		\n\t		vmulpd		0x100(%%rbx),%%ymm13,%%ymm13		/* a[jp+p12]*c12 */\n\t"\
	"vfnmadd231pd	0x020(%%rbx),%%ymm7,%%ymm4		/* a[jp+p12]*s12 */		\n\t	vfnmadd231pd	0x120(%%rbx),%%ymm15,%%ymm12		/* a[jp+p12]*s12 */\n\t"\
	" vfmadd231pd	0x020(%%rbx),%%ymm6,%%ymm5		/* a[jt+p12]*s12 */		\n\t	 vfmadd231pd	0x120(%%rbx),%%ymm14,%%ymm13		/* a[jt+p12]*s12 */\n\t"\
		"vmovaps	%%ymm5		,0x020(%%rcx)	/* store it */				\n\t		vmovaps		%%ymm13		,0x120(%%rcx)	/* store it */\n\t"\
		"vmovaps	%%ymm4		,     (%%rcx)	/* store rt */				\n\t		vmovaps		%%ymm12		,0x100(%%rcx)	/* store rt */\n\t"\
		"\n\t"\
	"subq	$0x040 ,%%rbx	\n\t"/* c5 */											/* c7  */\
		"vmovaps	0x040(%%rcx),%%ymm4		/* a[jt+p4] */					\n\t		vmovaps		0x140(%%rcx),%%ymm12		/* a[jt+p4] */\n\t"\
		"vmovaps	0x060(%%rcx),%%ymm5		/* a[jp+p4] */					\n\t		vmovaps		0x160(%%rcx),%%ymm13		/* a[jp+p4] */\n\t"\
		"vmovaps		%%ymm4	,%%ymm6		/* ymm4 <- cpy a[jt+p4] */		\n\t		vmovaps			%%ymm12	,%%ymm14		/* ymm12 <- cpy a[jt+p4] */\n\t"\
		"vmovaps		%%ymm5	,%%ymm7		/* ymm5 <- cpy a[jp+p4] */		\n\t		vmovaps			%%ymm13	,%%ymm15		/* ymm13 <- cpy a[jp+p4] */\n\t"\
		"vmulpd		     (%%rbx),%%ymm4,%%ymm4		/* a[jt+p4]*c4 */		\n\t		vmulpd		0x100(%%rbx),%%ymm12,%%ymm12		/* a[jt+p4]*c4 */\n\t"\
		"vmulpd		     (%%rbx),%%ymm5,%%ymm5		/* a[jp+p4]*c4 */		\n\t		vmulpd		0x100(%%rbx),%%ymm13,%%ymm13		/* a[jp+p4]*c4 */\n\t"\
	"vfnmadd231pd	0x020(%%rbx),%%ymm7,%%ymm4		/* a[jp+p4]*s4 */		\n\t	vfnmadd231pd	0x120(%%rbx),%%ymm15,%%ymm12		/* a[jp+p4]*s4 */\n\t"\
	" vfmadd231pd	0x020(%%rbx),%%ymm6,%%ymm5		/* a[jt+p4]*s4 */		\n\t	 vfmadd231pd	0x120(%%rbx),%%ymm14,%%ymm13		/* a[jt+p4]*s4 */\n\t"\
	"prefetcht1	0x140(%%r14)\n\t"\
		"vmovaps	%%ymm5		,%%ymm7		/* ymm7 <- cpy t6 */			\n\t		vmovaps		%%ymm13		,%%ymm15		/* ymm15 <- cpy t6 */\n\t"\
		"vmovaps	%%ymm4		,%%ymm6		/* ymm6 <- cpy t5 */			\n\t		vmovaps		%%ymm12		,%%ymm14		/* ymm14 <- cpy t5 */\n\t"\
		"vsubpd		     (%%rcx),%%ymm4,%%ymm4		/* ~t7 <- t5 -rt */		\n\t		vsubpd		0x100(%%rcx),%%ymm12,%%ymm12		/* ~t7 <- t5 -rt */\n\t"\
		"vsubpd		0x020(%%rcx),%%ymm5,%%ymm5		/* ~t8 <- t6 -it */		\n\t		vsubpd		0x120(%%rcx),%%ymm13,%%ymm13		/* ~t8 <- t6 -it */\n\t"\
		"vaddpd		     (%%rcx),%%ymm6,%%ymm6		/* ~t5 <- t5 +rt */		\n\t		vaddpd		0x100(%%rcx),%%ymm14,%%ymm14		/* ~t5 <- t5 +rt */\n\t"\
		"vaddpd		0x020(%%rcx),%%ymm7,%%ymm7		/* ~t6 <- t6 +it */		\n\t		vaddpd		0x120(%%rcx),%%ymm15,%%ymm15		/* ~t6 <- t6 +it */\n\t"\
		"\n\t"\
	/* Finish radix-4 butterfly and store results into temp-array slots: */\
		"vsubpd		%%ymm6,%%ymm0,%%ymm0	/*~t5 */						\n\t		vsubpd		%%ymm14,%%ymm8 ,%%ymm8 	/*~t5 */\n\t"\
		"vsubpd		%%ymm7,%%ymm1,%%ymm1	/*~t6 */						\n\t		vsubpd		%%ymm15,%%ymm9 ,%%ymm9 	/*~t6 */\n\t"\
		"vsubpd		%%ymm5,%%ymm2,%%ymm2	/*~t3 */						\n\t		vsubpd		%%ymm13,%%ymm10,%%ymm10	/*~t3 */\n\t"\
		"vsubpd		%%ymm4,%%ymm3,%%ymm3	/*~t8 */						\n\t		vsubpd		%%ymm12,%%ymm11,%%ymm11	/*~t8 */\n\t"\
		"vmovaps	%%ymm0,0x080(%%rcx)	/* a[jt+p8 ] <- ~t5 */				\n\t		vmovaps		%%ymm8 ,0x180(%%rcx)	/* a[jt+p8 ] <- ~t5 */\n\t"\
		"vmovaps	%%ymm1,0x0a0(%%rcx)	/* a[jp+p8 ] <- ~t6 */				\n\t		vmovaps		%%ymm9 ,0x1a0(%%rcx)	/* a[jp+p8 ] <- ~t6 */\n\t"\
		"vmovaps	%%ymm2,0x040(%%rcx)	/* a[jt+p4 ] <- ~t3 */				\n\t		vmovaps		%%ymm10,0x140(%%rcx)	/* a[jt+p4 ] <- ~t3 */\n\t"\
		"vmovaps	%%ymm3,0x0e0(%%rcx)	/* a[jp+p12] <- ~t8 */				\n\t		vmovaps		%%ymm11,0x1e0(%%rcx)	/* a[jp+p12] <- ~t8 */\n\t"\
	"vmovaps	%%ymm12,(%%rax)	\n\t"/* spill ymm12 to allow 2.0 to use a reg */\
	"vmovaps	(%%rsi),%%ymm12	\n\t"/* two */\
	"vfmadd132pd	%%ymm12,%%ymm0,%%ymm6	/*~t1 */					\n\t		vfmadd132pd	%%ymm12,%%ymm8 ,%%ymm14	/*~t1 */\n\t"\
	"vfmadd132pd	%%ymm12,%%ymm1,%%ymm7	/*~t2 */					\n\t		vfmadd132pd	%%ymm12,%%ymm9 ,%%ymm15	/*~t2 */\n\t"\
	"vfmadd132pd	%%ymm12,%%ymm2,%%ymm5	/*~t7 */					\n\t		vfmadd132pd	%%ymm12,%%ymm10,%%ymm13	/*~t7 */\n\t"\
	"vfmadd132pd	%%ymm12,%%ymm3,%%ymm4	/*~t4 */					\n\t		vfmadd132pd	(%%rax),%%ymm11,%%ymm12	/*~t4 */\n\t"\
		"vmovaps	%%ymm6,     (%%rcx)/* a[jt    ] <- ~t1 */				\n\t		vmovaps		%%ymm14,0x100(%%rcx)	/* a[jt    ] <- ~t1 */\n\t"\
		"vmovaps	%%ymm7,0x020(%%rcx)	/* a[jp    ] <- ~t2 */				\n\t		vmovaps		%%ymm15,0x120(%%rcx)	/* a[jp    ] <- ~t2 */\n\t"\
		"vmovaps	%%ymm5,0x0c0(%%rcx)	/* a[jt+p12] <- ~t7 */				\n\t		vmovaps		%%ymm13,0x1c0(%%rcx)	/* a[jt+p12] <- ~t7 */\n\t"\
		"vmovaps	%%ymm4,0x060(%%rcx)	/* a[jp+p4 ] <- ~t4 */				\n\t		vmovaps		%%ymm12,0x160(%%rcx)	/* a[jp+p4 ] <- ~t4 */\n\t"\
		"\n\t"\
	/**************************************************************************************/\
	/*...and now do four more radix-4 transforms, including the internal twiddle factors: */\
	/**************************************************************************************/\
		"\n\t"\
	"/*...Block 1: t1,9,17,25 */\n\t											/*...Block 3: t5,13,21,29: All rax-offsets incr +0x40 in rcol w.r.to lcol: */\n\t"\
		"movq	%[__r1],%%rax							\n\t				vmovaps		0x080(%%rax),%%ymm8 		/* t5  */\n\t"\
		"movq	%[__r9],%%rbx							\n\t				vmovaps		0x0a0(%%rax),%%ymm9 		/* t6  */\n\t"\
		"movq	%[__isrt2],%%rdi	\n\t"/* &two still in rsi */		"	vmovaps		0x1a0(%%rax),%%ymm11		/* t14 */\n\t"\
	"vmovaps	(%%rsi),%%ymm15	\n\t"/* two */\
		"vmovaps		 (%%rax),%%ymm0		/* t1  */\n\t					vmovaps		0x180(%%rax),%%ymm10		/* t13 */\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1		/* t2  */\n\t					vsubpd		%%ymm11,%%ymm8 ,%%ymm8 		/* t5 =t5 -t14 */\n\t"\
		"vmovaps	0x100(%%rax),%%ymm2		/* t9  */\n\t					vsubpd		%%ymm10,%%ymm9 ,%%ymm9 		/* t14=t6 -t13 */\n\t"\
		"vmovaps	0x120(%%rax),%%ymm3		/* t14 */\n\t					vfmadd132pd	%%ymm15,%%ymm8 ,%%ymm11	/* t13=t5 +t14 */\n\t"\
		"\n\t																vfmadd132pd	%%ymm15,%%ymm9 ,%%ymm10	/* t6 =t6 +t13 */\n\t"\
		"vsubpd		%%ymm2,%%ymm0,%%ymm0		/* t9 =t1 -t9  */\n\t		vmovaps		0x280(%%rax),%%ymm12		/* t21 */\n\t"\
		"vsubpd		%%ymm3,%%ymm1,%%ymm1		/* t14=t2 -t14 */\n\t		vmovaps		0x2a0(%%rax),%%ymm13		/* t22 */\n\t"\
	"vfmadd132pd	%%ymm15,%%ymm0,%%ymm2		/* t1 =t1 +t9  */\n\t		vmovaps		0x380(%%rax),%%ymm14		/* t29 */\n\t"\
	"vfmadd132pd	%%ymm15,%%ymm1,%%ymm3		/* t2 =t2 +t14 */\n\t		vmovaps		0x3a0(%%rax),%%ymm15		/* t30 */\n\t"\
		"\n\t																vsubpd		0x2a0(%%rax),%%ymm12,%%ymm12		/* t21-t22 */\n\t"\
		"vmovaps	0x200(%%rax),%%ymm4		/* t17 */\n\t					vaddpd		0x280(%%rax),%%ymm13,%%ymm13		/* t22+t21 */\n\t"\
		"vmovaps	0x220(%%rax),%%ymm5		/* t18 */\n\t					vmulpd		(%%rdi)		,%%ymm12,%%ymm12	/* t21 = (t21-t22)*ISRT2 */\n\t"\
		"vmovaps	0x300(%%rax),%%ymm6		/* t25 */\n\t					vmulpd		(%%rdi)		,%%ymm13,%%ymm13	/* t22 = (t22+t21)*ISRT2 */\n\t"\
		"vmovaps	0x320(%%rax),%%ymm7		/* t26 */\n\t					vaddpd		0x3a0(%%rax),%%ymm14,%%ymm14		/* t29+t30 */\n\t"\
		"\n\t																vsubpd		0x380(%%rax),%%ymm15,%%ymm15		/* t30-t29 */\n\t"\
	"prefetcht1	0x180(%%r14)\n\t"\
	"vmovaps	%%ymm0,(%%rax)	\n\t"/* spill ymm0 to allow 2.0-in-reg */"	vmulpd		(%%rdi),%%ymm14,%%ymm14	/*  rt = (t29+t30)*ISRT2 */\n\t"\
	"vmovaps	(%%rsi),%%ymm0	\n\t"/* two */							"	vmulpd		(%%rdi),%%ymm15,%%ymm15	/*  it = (t30-t29)*ISRT2 */\n\t"\
		"vsubpd		%%ymm6,%%ymm4,%%ymm4		/* t25=t17-t25 */\n\t		vsubpd		%%ymm14,%%ymm12,%%ymm12		/* t21=t21-rt */\n\t"\
		"vsubpd		%%ymm7,%%ymm5,%%ymm5		/* t26=t18-t26 */\n\t		vsubpd		%%ymm15,%%ymm13,%%ymm13		/* t22=t22-it */\n\t"\
	"vfmadd132pd	%%ymm0,%%ymm4,%%ymm6		/* t17=t17+t25 */\n\t		vfmadd132pd	%%ymm0,%%ymm12,%%ymm14		/* t29=t21+rt */\n\t"\
	"vfmadd132pd	%%ymm0,%%ymm5,%%ymm7		/* t18=t18+t26 */\n\t		vfmadd132pd	%%ymm0,%%ymm13,%%ymm15		/* t30=t22+it */\n\t"\
		"movq	%[__r17],%%rcx						\n\t		movq		%[__r25],%%rdx			\n\t"\
		"vsubpd		%%ymm6,%%ymm2,%%ymm2						\n\t		vsubpd		%%ymm12,%%ymm8 ,%%ymm8 			\n\t"\
		"vsubpd		%%ymm7,%%ymm3,%%ymm3						\n\t		vsubpd		%%ymm13,%%ymm10,%%ymm10			\n\t"\
		"vmovaps		%%ymm2,     (%%rcx)						\n\t		vmovaps		%%ymm8 ,0x080(%%rcx)		\n\t"\
		"vmovaps		%%ymm3,0x020(%%rcx)						\n\t		vmovaps		%%ymm10,0x0a0(%%rcx)		\n\t"\
	"vfmadd132pd	%%ymm0,%%ymm2,%%ymm6						\n\t		vfmadd132pd	%%ymm0,%%ymm8 ,%%ymm12			\n\t"\
	"vfmadd132pd	%%ymm0,%%ymm3,%%ymm7						\n\t		vfmadd132pd	%%ymm0,%%ymm10,%%ymm13			\n\t"\
	"vmovaps	%%ymm0,%%ymm10	\n\t"/* move two into ymm10 */\
	"vmovaps	(%%rax),%%ymm0	\n\t"/* restore spill */\
		"vmovaps		%%ymm6,     (%%rax)						\n\t		vmovaps		%%ymm12,0x080(%%rax)		\n\t"\
		"vmovaps		%%ymm7,0x020(%%rax)						\n\t		vmovaps		%%ymm13,0x0a0(%%rax)		\n\t"\
		"vsubpd		%%ymm5,%%ymm0,%%ymm0						\n\t		vsubpd		%%ymm15,%%ymm11,%%ymm11			\n\t"\
		"vsubpd		%%ymm4,%%ymm1,%%ymm1						\n\t		vsubpd		%%ymm14,%%ymm9 ,%%ymm9 			\n\t"\
	"vfmadd132pd	%%ymm10,%%ymm0,%%ymm5						\n\t		vfmadd132pd	%%ymm10,%%ymm11,%%ymm15			\n\t"\
	"vfmadd132pd	%%ymm10,%%ymm1,%%ymm4						\n\t		vfmadd132pd	%%ymm10,%%ymm9 ,%%ymm14			\n\t"\
		"vmovaps		%%ymm0,     (%%rbx)						\n\t		vmovaps		%%ymm11,0x080(%%rbx)		\n\t"\
		"vmovaps		%%ymm1,0x020(%%rdx)						\n\t		vmovaps		%%ymm9 ,0x0a0(%%rdx)		\n\t"\
		"vmovaps		%%ymm5,     (%%rdx)						\n\t		vmovaps		%%ymm15,0x080(%%rdx)		\n\t"\
		"vmovaps		%%ymm4,0x020(%%rbx)						\n\t		vmovaps		%%ymm14,0x0a0(%%rbx)		\n\t"\
		"\n\t"\
	/*...Block 2: t3,11,19,27 */\
		"addq	$0x040,%%rax		\n\t"/* r3  */\
		"addq	$0x040,%%rbx							\n\t"\
		"addq	$0x040,%%rcx							\n\t"\
		"addq	$0x040,%%rdx							\n\t"\
		"leaq	0x020(%%rdi),%%rdi	\n\t"/* cc0, from isrt2 */\
																			/*...Block 4: t7,15,23,31 */\
		"vmovaps	0x200(%%rax),%%ymm4		/* t19 */		\n\t		vmovaps		0x280(%%rax),%%ymm12		/* t23 */\n\t"\
		"vmovaps	0x220(%%rax),%%ymm5		/* t20 */		\n\t		vmovaps		0x2a0(%%rax),%%ymm13		/* t24 */\n\t"\
		"vmovaps	0x300(%%rax),%%ymm6	/* t27 */			\n\t		vmovaps		0x380(%%rax),%%ymm14		/* t31 */\n\t"\
		"vmovaps	0x320(%%rax),%%ymm7	/* t28 */			\n\t		vmovaps		0x3a0(%%rax),%%ymm15		/* t32 */\n\t"\
		"vmovaps	%%ymm4,%%ymm0		/* copy t19 */		\n\t		vmovaps		%%ymm12,%%ymm8 		/* copy t23 */\n\t"\
		"vmovaps	%%ymm5,%%ymm1		/* copy t20 */		\n\t		vmovaps		%%ymm13,%%ymm9 		/* copy t24 */\n\t"\
	/*	"vmovaps	%%ymm6,%%ymm2	/x copy t27 x/			\n\t		vmovaps		%%ymm14,%%ymm10		/x copy t31 x/\n\t" */\
		"vmovaps	%%ymm7,%%ymm3	/* copy t28 */			\n\t		vmovaps		%%ymm15,%%ymm11		/* copy t32 */\n\t"\
		"\n\t"\
	"vmovaps	(%%rdi),%%ymm2	\n\t	vmovaps	0x020(%%rdi),%%ymm10	\n\t"/* cc0,ss0 data shared by both columns */\
		"vmulpd			%%ymm2 ,%%ymm4,%%ymm4	/* t19*c */	\n\t		vmulpd			%%ymm10,%%ymm12,%%ymm12		/* t23*s */\n\t"\
		"vmulpd			%%ymm10,%%ymm6,%%ymm6	/* t27*s */	\n\t		vmulpd			%%ymm2 ,%%ymm14,%%ymm14		/* t31*c */\n\t"\
		"vmulpd			%%ymm2 ,%%ymm5,%%ymm5	/* t20*c */	\n\t		vmulpd			%%ymm10,%%ymm13,%%ymm13		/* t24*s */\n\t"\
		"vmulpd			%%ymm10,%%ymm7,%%ymm7	/* t28*s */	\n\t		vmulpd			%%ymm2 ,%%ymm15,%%ymm15		/* t32*c */\n\t"\
	"vfnmadd231pd		%%ymm10,%%ymm1,%%ymm4	/* ~t19 */	\n\t	vfnmadd231pd		%%ymm2 ,%%ymm9 ,%%ymm12		/* ~t23 */\n\t"\
	"vfnmadd231pd		%%ymm2 ,%%ymm3,%%ymm6	/* rt */	\n\t	vfnmadd231pd		%%ymm10,%%ymm11,%%ymm14		/* rt */\n\t"\
	" vfmadd231pd		%%ymm10,%%ymm0,%%ymm5	/* ~t20 */	\n\t	 vfmadd231pd		%%ymm2 ,%%ymm8 ,%%ymm13		/* ~t24 */\n\t"\
	" vfmadd231pd	0x300(%%rax),%%ymm2,%%ymm7	/* it */	\n\t	 vfmadd231pd	0x380(%%rax),%%ymm10,%%ymm15	/* it */\n\t"\
	"prefetcht1	0x1c0(%%r14)\n\t"\
		"\n\t"\
	"subq	$0x020,%%rdi	\n\t"/* isrt2, from cc0 */\
		"vsubpd		%%ymm6,%%ymm4,%%ymm4		/*~t27=t19-rt */\n\t		vsubpd		%%ymm14,%%ymm12,%%ymm12		/*~t23=t23-rt */\n\t"\
		"vsubpd		%%ymm7,%%ymm5,%%ymm5		/*~t28=t20-it */\n\t		vsubpd		%%ymm15,%%ymm13,%%ymm13		/*~t24=t24-it */\n\t"\
	"vmovaps	(%%rsi),%%ymm2	\n\t"/* two */\
	"vfmadd132pd	%%ymm2,%%ymm4,%%ymm6		/*~t19=t19+rt */\n\t	vfmadd132pd	%%ymm2,%%ymm12,%%ymm14		/*~t31=t23+rt */\n\t"\
	"vfmadd132pd	%%ymm2,%%ymm5,%%ymm7		/*~t20=t20+it */\n\t	vfmadd132pd	%%ymm2,%%ymm13,%%ymm15		/*~t32=t24+it */\n\t"\
		"\n\t"\
		"vmovaps	0x100(%%rax),%%ymm2		/* t11 */\n\t					vmovaps		0x180(%%rax),%%ymm10		/* t15 */\n\t"\
		"vmovaps	0x120(%%rax),%%ymm3		/* t12 */\n\t					vmovaps		0x1a0(%%rax),%%ymm11		/* t16 */\n\t"\
		"vsubpd		0x120(%%rax),%%ymm2,%%ymm2		/* t11-t12 */\n\t		vaddpd		0x1a0(%%rax),%%ymm10,%%ymm10	/* t15+t16 */\n\t"\
		"vaddpd		0x100(%%rax),%%ymm3,%%ymm3		/* t12+t11 */\n\t		vsubpd		0x180(%%rax),%%ymm11,%%ymm11	/* t16-t15 */\n\t"\
		"vmulpd		(%%rdi),%%ymm2,%%ymm2	/* rt = (t11-t12)*ISRT2 */\n\t	vmulpd		(%%rdi)		,%%ymm10,%%ymm10	/* rt = (t15+t16)*ISRT2 */\n\t"\
		"vmulpd		(%%rdi),%%ymm3,%%ymm3	/* it = (t12+t11)*ISRT2 */\n\t	vmulpd		(%%rdi)		,%%ymm11,%%ymm11	/* it = (t16-t15)*ISRT2 */\n\t"\
		"\n\t"\
		"vmovaps	     (%%rax),%%ymm0		/* t3  */\n\t					vmovaps		0x080(%%rax),%%ymm8 		/* t7  */\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1		/* t4  */\n\t					vmovaps		0x0a0(%%rax),%%ymm9 		/* t8  */\n\t"\
		"\n\t"\
		"vsubpd		      %%ymm2,%%ymm0,%%ymm0		/*~t11=t3 -rt */\n\t	vsubpd		     %%ymm10,%%ymm8 ,%%ymm8 	/*~t7 =t7 -rt */\n\t"\
		"vsubpd		      %%ymm3,%%ymm1,%%ymm1		/*~t12=t4 -it */\n\t	vsubpd		     %%ymm11,%%ymm9 ,%%ymm9 	/*~t8 =t8 -it */\n\t"\
		"vaddpd		     (%%rax),%%ymm2,%%ymm2		/*~t3 =rt +t3 */\n\t	vaddpd		0x080(%%rax),%%ymm10,%%ymm10	/*~t15=rt +t7 */\n\t"\
		"vaddpd		0x020(%%rax),%%ymm3,%%ymm3		/*~t4 =it +t4 */\n\t	vaddpd		0x0a0(%%rax),%%ymm11,%%ymm11	/*~t16=it +t8 */\n\t"\
		"\n\t"\
		"vsubpd		%%ymm6,%%ymm2,%%ymm2		/* t3 -t19 */\n\t			vsubpd		%%ymm12		,%%ymm8 ,%%ymm8 	/* t7 -t23 */\n\t"\
		"vsubpd		%%ymm7,%%ymm3,%%ymm3		/* t4 -t20 */\n\t			vsubpd		%%ymm13		,%%ymm9 ,%%ymm9 	/* t8 -t24 */\n\t"\
	"vfmadd132pd	(%%rsi),%%ymm2,%%ymm6		/* t3 +t19 */\n\t		vfmadd132pd	(%%rsi),%%ymm8,%%ymm12		/* t7 +t23 */\n\t"\
	"vfmadd132pd	(%%rsi),%%ymm3,%%ymm7		/* t4 +t20 */\n\t		vfmadd132pd	(%%rsi),%%ymm9,%%ymm13		/* t8 +t24 */\n\t"\
		"vmovaps		%%ymm2,     (%%rcx)						\n\t		vmovaps		%%ymm8 ,0x080(%%rcx)		\n\t"\
		"vmovaps		%%ymm3,0x020(%%rcx)						\n\t		vmovaps		%%ymm9 ,0x0a0(%%rcx)		\n\t"\
		"vmovaps		%%ymm6,     (%%rax)						\n\t		vmovaps		%%ymm12,0x080(%%rax)		\n\t"\
		"vmovaps		%%ymm7,0x020(%%rax)						\n\t		vmovaps		%%ymm13,0x0a0(%%rax)		\n\t"\
		"vsubpd		%%ymm5,%%ymm0,%%ymm0						\n\t		vsubpd		%%ymm15,%%ymm10,%%ymm10			\n\t"\
		"vsubpd		%%ymm4,%%ymm1,%%ymm1						\n\t		vsubpd		%%ymm14,%%ymm11,%%ymm11			\n\t"\
	"vfmadd132pd	(%%rsi),%%ymm0,%%ymm5	\n\t				vfmadd132pd	(%%rsi),%%ymm10,%%ymm15		\n\t"\
	"vfmadd132pd	(%%rsi),%%ymm1,%%ymm4	\n\t				vfmadd132pd	(%%rsi),%%ymm11,%%ymm14		\n\t"\
		"vmovaps		%%ymm0,     (%%rbx)						\n\t		vmovaps		%%ymm10,0x080(%%rbx)		\n\t"\
		"vmovaps		%%ymm1,0x020(%%rdx)						\n\t		vmovaps		%%ymm11,0x0a0(%%rdx)		\n\t"\
		"vmovaps		%%ymm5,     (%%rdx)						\n\t		vmovaps		%%ymm15,0x080(%%rdx)		\n\t"\
		"vmovaps		%%ymm4,0x020(%%rbx)						\n\t		vmovaps		%%ymm14,0x0a0(%%rbx)		\n\t"\
		 :					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__add1] "m" (Xadd1)\
		 ,[__add2] "m" (Xadd2)\
		 ,[__add3] "m" (Xadd3)\
		 ,[__r1] "m" (Xr1)\
		 ,[__r9] "m" (Xr9)\
		 ,[__r17] "m" (Xr17)\
		 ,[__r25] "m" (Xr25)\
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__cc0] "m" (Xcc0)\
		 ,[__c1] "m" (Xc1)\
		 ,[__c2] "m" (Xc2)\
		 ,[__c3] "m" (Xc3)\
		 ,[__c4] "m" (Xc4)\
		 ,[__c5] "m" (Xc5)\
		 ,[__c6] "m" (Xc6)\
		 ,[__c7] "m" (Xc7)\
		 ,[__c8] "m" (Xc8)\
		 ,[__c9] "m" (Xc9)\
		 ,[__c10] "m" (Xc10)\
		 ,[__c11] "m" (Xc11)\
		 ,[__c12] "m" (Xc12)\
		 ,[__c13] "m" (Xc13)\
		 ,[__c14] "m" (Xc14)\
		 ,[__c15] "m" (Xc15)\
	 ,[__pfetch_dist] "m" (Xpfetch_dist)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r14","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"		/* Clobbered registers */\
	);\
	}

	// Since the add0/1 block addresses advance in opposite directions, our prefetch scheme here is like so:
	// In fwd DIF prefetch  ahead 8 cache lines' worth of data w.r.to add0;
	// In inv DIT prefetch behind 8 cache lines' worth of data w.r.to add1;
	/****** To-Do: Need a 16-register version! ****************/
	#define SSE2_RADIX16_WRAPPER_DIT(Xadd0,Xadd1,Xadd2,Xadd3,Xr1,Xr9,Xr17,Xr25,Xisrt2,Xcc0,Xc1,Xc2,Xc3,Xc4,Xc5,Xc6,Xc7,Xc8,Xc9,Xc10,Xc11,Xc12,Xc13,Xc14,Xc15,Xpfetch_dist)\
	{\
	__asm__ volatile (\
	/*...Block 1: r1,9,17,25 */\
		"movq	%[__r1],%%rax					\n\t"\
	"movq	%[__add1],%%r14	\n\t"\
	"movslq	%[__pfetch_dist],%%r14	\n\t"\
	"leaq	(%%rax,%%r14,8),%%r14	\n\t"	/* Block 1 [base-address + data-fetch-ahead index] */\
	"prefetcht1	-0x40(%%r14)\n\t"\
		"movq	%%rax,%%rbx						\n\t		leaq	0x880(%%rax),%%rsi	\n\t"/* two */\
		"movq	%%rax,%%rcx						\n\t"	/*...Block 3: r3,11,19,27 */\
		"movq	%%rax,%%rdx						\n\t		leaq	0x40(%%rax),%%r10				\n\t"\
		"addq	$0x200,%%rbx					\n\t		leaq	0x40(%%rbx),%%r11				\n\t"\
		"addq	$0x100,%%rcx					\n\t		leaq	0x40(%%rcx),%%r12				\n\t"\
		"addq	$0x300,%%rdx					\n\t		leaq	0x40(%%rdx),%%r13				\n\t"\
		"/* SSE2_RADIX4_DIT_IN_PLACE(): */		\n\t		/* SSE2_RADIX4_DIT_IN_PLACE(): */		\n\t"\
	"vmovaps	(%%rsi),%%ymm13	\n\t"/* two */\
		"vmovaps	     (%%rax),%%ymm2			\n\t		vmovaps	     (%%r10),%%ymm10			\n\t"\
		"vmovaps	0x020(%%rax),%%ymm3			\n\t		vmovaps	0x020(%%r10),%%ymm11			\n\t"\
		"vmovaps	     (%%rbx),%%ymm0			\n\t		vmovaps	     (%%r11),%%ymm8 			\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm1			\n\t		vmovaps	0x020(%%r11),%%ymm9 			\n\t"\
		"vsubpd		%%ymm0,%%ymm2,%%ymm2		\n\t		vsubpd		%%ymm8 ,%%ymm10,%%ymm10		\n\t"\
		"vsubpd		%%ymm1,%%ymm3,%%ymm3		\n\t		vsubpd		%%ymm9 ,%%ymm11,%%ymm11		\n\t"\
		"vmovaps	     (%%rcx),%%ymm6			\n\t		vmovaps	     (%%r12),%%ymm14			\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm7			\n\t		vmovaps	0x020(%%r12),%%ymm15			\n\t"\
		"vmovaps	     (%%rdx),%%ymm4			\n\t		vmovaps	     (%%r13),%%ymm12			\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm5			\n\t"	/*	vmovaps	0x020(%%r13),%%ymm13	*/		\
		"vsubpd		%%ymm4,%%ymm6,%%ymm6		\n\t		vsubpd		 %%ymm12,%%ymm14,%%ymm14	\n\t"\
		"vsubpd		%%ymm5,%%ymm7,%%ymm7		\n\t		vsubpd	0x020(%%r13),%%ymm15,%%ymm15	\n\t"\
	"vfmadd132pd	%%ymm13,%%ymm2,%%ymm0		\n\t	vfmadd132pd		%%ymm13,%%ymm10,%%ymm8 		\n\t"\
	"vfmadd132pd	%%ymm13,%%ymm3,%%ymm1		\n\t	vfmadd132pd		%%ymm13,%%ymm11,%%ymm9 		\n\t"\
	"vfmadd132pd	%%ymm13,%%ymm6,%%ymm4		\n\t	vfmadd132pd		 %%ymm13,%%ymm14,%%ymm12		\n\t"\
	"vfmadd132pd	%%ymm13,%%ymm7,%%ymm5		\n\t	vfmadd132pd	0x020(%%r13),%%ymm15,%%ymm13		\n\t"\
	"prefetcht1	-0xc0(%%r14)\n\t"\
		"vsubpd	      %%ymm4,%%ymm0,%%ymm0		\n\t		vsubpd	      %%ymm12,%%ymm8 ,%%ymm8 	\n\t"\
		"vsubpd	      %%ymm5,%%ymm1,%%ymm1		\n\t		vsubpd	      %%ymm13,%%ymm9 ,%%ymm9 	\n\t"\
		"vsubpd	      %%ymm7,%%ymm2,%%ymm2		\n\t		vsubpd	      %%ymm15,%%ymm10,%%ymm10	\n\t"\
		"vsubpd	      %%ymm6,%%ymm3,%%ymm3		\n\t		vsubpd	      %%ymm14,%%ymm11,%%ymm11	\n\t"\
		"vmovaps	%%ymm0,     (%%rbx)			\n\t		vmovaps	%%ymm8 ,     (%%r11)			\n\t"\
		"vmovaps	%%ymm1,0x020(%%rbx)			\n\t		vmovaps	%%ymm9 ,0x020(%%r11)			\n\t"\
		"vmovaps	%%ymm2,     (%%rdx)			\n\t		vmovaps	%%ymm10,     (%%r13)			\n\t"\
		"vmovaps	%%ymm3,0x020(%%rcx)			\n\t		vmovaps	%%ymm11,0x020(%%r12)			\n\t"\
	"vmovaps	(%%rsi),%%ymm0	\n\t"/* two */\
	"vfmadd213pd	(%%rbx),%%ymm0,%%ymm4		\n\t	vfmadd132pd	%%ymm0,%%ymm8 ,%%ymm12		\n\t"\
	"vfmadd132pd	 %%ymm0,%%ymm1,%%ymm5		\n\t	vfmadd132pd	%%ymm0,%%ymm9 ,%%ymm13		\n\t"\
	"vfmadd132pd	 %%ymm0,%%ymm2,%%ymm7		\n\t	vfmadd132pd	%%ymm0,%%ymm10,%%ymm15		\n\t"\
	"vfmadd132pd	 %%ymm0,%%ymm3,%%ymm6		\n\t	vfmadd132pd	%%ymm0,%%ymm11,%%ymm14		\n\t"\
		"vmovaps	%%ymm4,     (%%rax)			\n\t		vmovaps	%%ymm12,     (%%r10)			\n\t"\
		"vmovaps	%%ymm5,0x020(%%rax)			\n\t		vmovaps	%%ymm13,0x020(%%r10)			\n\t"\
		"vmovaps	%%ymm7,     (%%rcx)			\n\t		vmovaps	%%ymm15,     (%%r12)			\n\t"\
		"vmovaps	%%ymm6,0x020(%%rdx)			\n\t		vmovaps	%%ymm14,0x020(%%r13)			\n\t"\
	/*...Block 2: r5,13,21,29 */						/*...Block 4: r7,15,23,31 */\
		"addq	$0x80,%%rax						\n\t		leaq	0x40(%%rax),%%r10				\n\t"\
		"addq	$0x80,%%rbx						\n\t		leaq	0x40(%%rbx),%%r11				\n\t"\
		"addq	$0x80,%%rcx						\n\t		leaq	0x40(%%rcx),%%r12				\n\t"\
		"addq	$0x80,%%rdx						\n\t		leaq	0x40(%%rdx),%%r13				\n\t"\
		"/* SSE2_RADIX4_DIT_IN_PLACE(): */		\n\t		/* SSE2_RADIX4_DIT_IN_PLACE(): */		\n\t"\
	"vmovaps	(%%rsi),%%ymm13	\n\t"/* two */\
		"vmovaps	     (%%rax),%%ymm2			\n\t		vmovaps	     (%%r10),%%ymm10			\n\t"\
		"vmovaps	0x020(%%rax),%%ymm3			\n\t		vmovaps	0x020(%%r10),%%ymm11			\n\t"\
		"vmovaps	     (%%rbx),%%ymm0			\n\t		vmovaps	     (%%r11),%%ymm8 			\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm1			\n\t		vmovaps	0x020(%%r11),%%ymm9 			\n\t"\
		"vsubpd		%%ymm0,%%ymm2,%%ymm2		\n\t		vsubpd		%%ymm8 ,%%ymm10,%%ymm10		\n\t"\
		"vsubpd		%%ymm1,%%ymm3,%%ymm3		\n\t		vsubpd		%%ymm9 ,%%ymm11,%%ymm11		\n\t"\
		"vmovaps	     (%%rcx),%%ymm6			\n\t		vmovaps	     (%%r12),%%ymm14			\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm7			\n\t		vmovaps	0x020(%%r12),%%ymm15			\n\t"\
		"vmovaps	     (%%rdx),%%ymm4			\n\t		vmovaps	     (%%r13),%%ymm12			\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm5			\n\t"	/*	vmovaps	0x020(%%r13),%%ymm13	*/		\
		"vsubpd		%%ymm4,%%ymm6,%%ymm6		\n\t		vsubpd		 %%ymm12,%%ymm14,%%ymm14	\n\t"\
		"vsubpd		%%ymm5,%%ymm7,%%ymm7		\n\t		vsubpd	0x020(%%r13),%%ymm15,%%ymm15	\n\t"\
	"vfmadd132pd	%%ymm13,%%ymm2,%%ymm0		\n\t	vfmadd132pd		%%ymm13,%%ymm10,%%ymm8 		\n\t"\
	"vfmadd132pd	%%ymm13,%%ymm3,%%ymm1		\n\t	vfmadd132pd		%%ymm13,%%ymm11,%%ymm9 		\n\t"\
	"vfmadd132pd	%%ymm13,%%ymm6,%%ymm4		\n\t	vfmadd132pd		 %%ymm13,%%ymm14,%%ymm12		\n\t"\
	"vfmadd132pd	%%ymm13,%%ymm7,%%ymm5		\n\t	vfmadd132pd	0x020(%%r13),%%ymm15,%%ymm13		\n\t"\
	"prefetcht1	-0x80(%%r14)\n\t"\
		"vsubpd	      %%ymm4,%%ymm0,%%ymm0		\n\t		vsubpd	      %%ymm12,%%ymm8 ,%%ymm8 	\n\t"\
		"vsubpd	      %%ymm5,%%ymm1,%%ymm1		\n\t		vsubpd	      %%ymm13,%%ymm9 ,%%ymm9 	\n\t"\
		"vsubpd	      %%ymm7,%%ymm2,%%ymm2		\n\t		vsubpd	      %%ymm15,%%ymm10,%%ymm10	\n\t"\
		"vsubpd	      %%ymm6,%%ymm3,%%ymm3		\n\t		vsubpd	      %%ymm14,%%ymm11,%%ymm11	\n\t"\
		"vmovaps	%%ymm0,     (%%rbx)			\n\t		vmovaps	%%ymm8 ,     (%%r11)			\n\t"\
		"vmovaps	%%ymm1,0x020(%%rbx)			\n\t		vmovaps	%%ymm9 ,0x020(%%r11)			\n\t"\
		"vmovaps	%%ymm2,     (%%rdx)			\n\t		vmovaps	%%ymm10,     (%%r13)			\n\t"\
		"vmovaps	%%ymm3,0x020(%%rcx)			\n\t		vmovaps	%%ymm11,0x020(%%r12)			\n\t"\
	"vmovaps	(%%rsi),%%ymm0	\n\t"/* two */\
	"vfmadd213pd	(%%rbx),%%ymm0,%%ymm4		\n\t	vfmadd132pd	%%ymm0,%%ymm8 ,%%ymm12		\n\t"\
	"vfmadd132pd	 %%ymm0,%%ymm1,%%ymm5		\n\t	vfmadd132pd	%%ymm0,%%ymm9 ,%%ymm13		\n\t"\
	"vfmadd132pd	 %%ymm0,%%ymm2,%%ymm7		\n\t	vfmadd132pd	%%ymm0,%%ymm10,%%ymm15		\n\t"\
	"vfmadd132pd	 %%ymm0,%%ymm3,%%ymm6		\n\t	vfmadd132pd	%%ymm0,%%ymm11,%%ymm14		\n\t"\
		"vmovaps	%%ymm4,     (%%rax)			\n\t		vmovaps	%%ymm12,     (%%r10)			\n\t"\
		"vmovaps	%%ymm5,0x020(%%rax)			\n\t		vmovaps	%%ymm13,0x020(%%r10)			\n\t"\
		"vmovaps	%%ymm7,     (%%rcx)			\n\t		vmovaps	%%ymm15,     (%%r12)			\n\t"\
		"vmovaps	%%ymm6,0x020(%%rdx)			\n\t		vmovaps	%%ymm14,0x020(%%r13)			\n\t"\
		"\n\t"\
	/***************************************************************************************************/\
	/*..and now do four more radix-4 transforms, including the internal and external twiddle factors.  */\
	/*  Write even-index 16-byte output pairs to a(j1), odd-index to a(j2), unpack same as on inputs.  */\
	/*  We do the last 2 radix-4 blocks first, to make the unpack-interleaving of outputs smoother.    */\
	/***************************************************************************************************/\
	"/* Main-array addresses still in add0,1; all rax-offsets incr +0x200 in rcol w.r.to lcol: */\n\t"\
	"prefetcht1	-0x100(%%r14)\n\t"\
	/*...Block 3: t3,11,19,27 -> r9,13,11,15: */		/*...Block 4: t7,15,23,31 -> r25,29,27,31: */\
		"movq		%[__r9],%%rax				\n\t"\
		"movq		%[__isrt2],%%rbx			\n\t"\
		"movq		%[__cc0],%%rcx				\n\t"\
	"vmovaps	(%%rcx),%%ymm3	\n\t	vmovaps	0x020(%%rcx),%%ymm11	\n\t"/* cc0,ss0 data shared by both columns */\
		"vmovaps	0x040(%%rax),%%ymm4			\n\t		vmovaps		0x240(%%rax),%%ymm12			\n\t"\
		"vmovaps	0x0c0(%%rax),%%ymm0			\n\t		vmovaps		0x2c0(%%rax),%%ymm8 			\n\t"\
		"vmovaps	0x060(%%rax),%%ymm5			\n\t		vmovaps		0x260(%%rax),%%ymm13			\n\t"\
		"vmovaps	0x0e0(%%rax),%%ymm1			\n\t		vmovaps		0x2e0(%%rax),%%ymm9 			\n\t"\
		"vmulpd		     %%ymm3 ,%%ymm4,%%ymm4	\n\t		vmulpd		     %%ymm11,%%ymm12,%%ymm12	\n\t"\
		"vmulpd		     %%ymm11,%%ymm0,%%ymm0	\n\t		vmulpd		     %%ymm3 ,%%ymm8 ,%%ymm8 	\n\t"\
		"vmulpd		     %%ymm3 ,%%ymm5,%%ymm5	\n\t		vmulpd		     %%ymm11,%%ymm13,%%ymm13	\n\t"\
		"vmulpd		     %%ymm11,%%ymm1,%%ymm1	\n\t		vmulpd		     %%ymm3 ,%%ymm9 ,%%ymm9 	\n\t"\
		"vmovaps	0x040(%%rax),%%ymm6			\n\t		vmovaps		0x240(%%rax),%%ymm14			\n\t"\
		"vmovaps	0x0c0(%%rax),%%ymm2			\n\t		vmovaps		0x2c0(%%rax),%%ymm10			\n\t"\
		"vmovaps	0x060(%%rax),%%ymm7			\n\t		vmovaps		0x260(%%rax),%%ymm15			\n\t"\
	"vfnmadd231pd	     %%ymm11,%%ymm6,%%ymm5	\n\t	vfnmadd231pd	     %%ymm3 ,%%ymm14,%%ymm13	\n\t"\
	"vfnmadd231pd	     %%ymm3 ,%%ymm2,%%ymm1	\n\t	vfnmadd231pd	     %%ymm11,%%ymm10,%%ymm9 	\n\t"\
	" vfmadd231pd	     %%ymm11,%%ymm7,%%ymm4	\n\t	 vfmadd231pd	     %%ymm3 ,%%ymm15,%%ymm12	\n\t"\
	" vfmadd231pd	0x0e0(%%rax),%%ymm3,%%ymm0	\n\t	 vfmadd231pd	0x2e0(%%rax),%%ymm11,%%ymm8 	\n\t"\
		"vmovaps	%%ymm5,%%ymm7				\n\t		vmovaps		%%ymm13,%%ymm15				\n\t"\
		"vmovaps	%%ymm4,%%ymm6				\n\t		vmovaps		%%ymm12,%%ymm14				\n\t"\
		"vaddpd		%%ymm0,%%ymm4,%%ymm4		\n\t		vaddpd		%%ymm8 ,%%ymm12,%%ymm12		\n\t"\
		"vaddpd		%%ymm1,%%ymm5,%%ymm5		\n\t		vaddpd		%%ymm9 ,%%ymm13,%%ymm13		\n\t"\
		"vsubpd		%%ymm0,%%ymm6,%%ymm6		\n\t		vsubpd		%%ymm8 ,%%ymm14,%%ymm14		\n\t"\
		"vsubpd		%%ymm1,%%ymm7,%%ymm7		\n\t		vsubpd		%%ymm9 ,%%ymm15,%%ymm15		\n\t"\
		"vmovaps	0x080(%%rax),%%ymm2			\n\t		vmovaps		0x280(%%rax),%%ymm10		\n\t"\
		"vmovaps	0x0a0(%%rax),%%ymm3			\n\t		vmovaps		0x2a0(%%rax),%%ymm11		\n\t"\
		"vmovaps	     (%%rax),%%ymm0			\n\t		vmovaps		0x200(%%rax),%%ymm8 		\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1			\n\t		vmovaps		0x220(%%rax),%%ymm9 		\n\t"\
		"vaddpd		0x0a0(%%rax),%%ymm2,%%ymm2	\n\t		vsubpd		0x2a0(%%rax),%%ymm10,%%ymm10\n\t"\
		"vsubpd		0x080(%%rax),%%ymm3,%%ymm3	\n\t		vaddpd		0x280(%%rax),%%ymm11,%%ymm11\n\t"\
		"vmulpd		     (%%rbx),%%ymm2,%%ymm2	\n\t		vmulpd		(%%rbx),%%ymm10,%%ymm10		\n\t"\
		"vmulpd		     (%%rbx),%%ymm3,%%ymm3	\n\t		vmulpd		(%%rbx),%%ymm11,%%ymm11		\n\t"\
		"vsubpd		%%ymm2,%%ymm0,%%ymm0		\n\t		vsubpd		%%ymm10,%%ymm8 ,%%ymm8 		\n\t"\
		"vsubpd		%%ymm3,%%ymm1,%%ymm1		\n\t		vsubpd		%%ymm11,%%ymm9 ,%%ymm9 		\n\t"\
	"vfmadd132pd	(%%rsi),%%ymm0,%%ymm2		\n\t	vfmadd132pd	(%%rsi),%%ymm8 ,%%ymm10			\n\t"\
	"vfmadd132pd	(%%rsi),%%ymm1,%%ymm3		\n\t	vfmadd132pd	(%%rsi),%%ymm9 ,%%ymm11			\n\t"\
	"addq	$0x240,%%rcx		\n\t"/* c1 from cc0  ; rcol has c1 */\
	"leaq	0x2a0(%%rbx),%%rdx	\n\t"/* c9 from isrt2; rcol has c9 */\
	/* Using the upper block(s) of the main array for temp-storage in the section below led to a nasty AVX bug to track down: */\
	/* In fermat-mod mode the 4 block addresses in ascending order are add0,1,2,3 with no 'gaps' between blocks, whereas for */\
	/* mersenne-mod the addresses in asc. order are add0,2,3,1 with a gap between contiguous-data-block pairs 0,2 and 3,1. Thus */\
	/* for fermat-mod we need [add2] as the base-address of the 'high-half' block for temp-storage; for mersenne-mod we need [add3]. */\
	/* In both cases we have that add2 < add3 so instead use (add2 - add1): > 0 for fermat-mod, < 0 for mersenne - to differentiate: */\
		"movq	%[__add2],%%rdi		\n\t"/* destroyable copy of add2 */\
		"movq	%[__add2],%%rbx		\n\t"\
		"subq	%[__add1],%%rdi		\n\t"/* rdi = (add2 - add1); if this yields a borrow (i.e. sets CF) it's mersenne, else fermat. */\
		"cmovcq %[__add3],%%rbx	\n\t" /* if CF set (i.e. h > l), copy source [add3] into dest (rbx), else leave dest = [add2]. */\
		/* rbx shared between rcol/lcol; rcx/rdx-offsets incr +0x80 in rcol for rest of block: */\
	"prefetcht1	-0x140(%%r14)\n\t"\
		"vsubpd		%%ymm4,%%ymm2,%%ymm2			\n\t		vsubpd		%%ymm14,%%ymm8 ,%%ymm8 		\n\t"\
		"vsubpd		%%ymm5,%%ymm3,%%ymm3			\n\t		vsubpd		%%ymm15,%%ymm9 ,%%ymm9 		\n\t"\
	"vfmadd132pd	(%%rsi),%%ymm2,%%ymm4			\n\t	vfmadd132pd	(%%rsi),%%ymm8 ,%%ymm14			\n\t"\
	"vfmadd132pd	(%%rsi),%%ymm3,%%ymm5			\n\t	vfmadd132pd	(%%rsi),%%ymm9 ,%%ymm15			\n\t"\
		"vmovaps	%%ymm2,     (%%rax)				\n\t		vmovaps		%%ymm8 ,0x200(%%rax)		\n\t"\
		"vmovaps	%%ymm3,0x020(%%rax)				\n\t		vmovaps		%%ymm9 ,0x220(%%rax)		\n\t"\
		"vmovaps	%%ymm4,%%ymm2					\n\t		vmovaps		%%ymm14,%%ymm8 		\n\t"\
		"vmovaps	%%ymm5,%%ymm3					\n\t		vmovaps		%%ymm15,%%ymm9 		\n\t"\
		"vmulpd		     (%%rcx),%%ymm4,%%ymm4		\n\t		vmulpd		0x100(%%rcx),%%ymm14,%%ymm14\n\t"\
		"vmulpd		     (%%rcx),%%ymm5,%%ymm5		\n\t		vmulpd		0x100(%%rcx),%%ymm15,%%ymm15\n\t"\
		"vsubpd		%%ymm7,%%ymm0,%%ymm0			\n\t		vsubpd		%%ymm13,%%ymm10,%%ymm10		\n\t"\
		"vsubpd		%%ymm6,%%ymm1,%%ymm1			\n\t		vsubpd		%%ymm12,%%ymm11,%%ymm11		\n\t"\
	"vfnmadd231pd	0x020(%%rcx),%%ymm2,%%ymm5		\n\t	vfnmadd231pd	0x120(%%rcx),%%ymm8 ,%%ymm15	\n\t"\
	" vfmadd231pd	0x020(%%rcx),%%ymm3,%%ymm4		\n\t	 vfmadd231pd	0x120(%%rcx),%%ymm9 ,%%ymm14	\n\t"\
	"vfmadd132pd	(%%rsi),%%ymm0,%%ymm7			\n\t	vfmadd132pd	(%%rsi),%%ymm10,%%ymm13			\n\t"\
	"vfmadd132pd	(%%rsi),%%ymm1,%%ymm6			\n\t	vfmadd132pd	(%%rsi),%%ymm11,%%ymm12			\n\t"\
		"vmovaps	%%ymm5,0x020(%%rbx)				\n\t		vmovaps		%%ymm15,0x060(%%rbx)		\n\t"\
		"vmovaps	%%ymm4,     (%%rbx)				\n\t		vmovaps		%%ymm14,0x040(%%rbx)		\n\t"\
		"vmovaps		 (%%rax),%%ymm4				\n\t		vmovaps		0x200(%%rax),%%ymm14		\n\t"\
		"vmovaps	0x020(%%rax),%%ymm5				\n\t		vmovaps		0x220(%%rax),%%ymm15		\n\t"\
		"vmovaps	%%ymm4,%%ymm2					\n\t		vmovaps		%%ymm14,%%ymm8 		\n\t"\
		"vmovaps	%%ymm5,%%ymm3					\n\t		vmovaps		%%ymm15,%%ymm9 		\n\t"\
		"vmulpd		     (%%rdx),%%ymm4,%%ymm4		\n\t		vmulpd		0x100(%%rdx),%%ymm14,%%ymm14\n\t"\
		"vmulpd		     (%%rdx),%%ymm5,%%ymm5		\n\t		vmulpd		0x100(%%rdx),%%ymm15,%%ymm15\n\t"\
	"vfnmadd231pd	0x020(%%rdx),%%ymm2,%%ymm5		\n\t	vfnmadd231pd	0x120(%%rdx),%%ymm8 ,%%ymm15	\n\t"\
	" vfmadd231pd	0x020(%%rdx),%%ymm3,%%ymm4		\n\t	 vfmadd231pd	0x120(%%rdx),%%ymm9 ,%%ymm14	\n\t"\
		"vmovaps	%%ymm5,0x120(%%rbx)				\n\t		vmovaps		%%ymm15,0x160(%%rbx)		\n\t"\
		"vmovaps	%%ymm4,0x100(%%rbx)				\n\t		vmovaps		%%ymm14,0x140(%%rbx)		\n\t"\
		"addq	$0x080,%%rcx	\n\t"/* c1  from c1; rcol has c7 */\
		"addq	$0x080,%%rdx	\n\t"/* c9  from c9; rcol has c15*/\
		"vmovaps	%%ymm7,%%ymm4					\n\t		vmovaps		%%ymm13,%%ymm8 		\n\t"\
		"vmovaps	%%ymm1,%%ymm5					\n\t		vmovaps		%%ymm11,%%ymm9 		\n\t"\
		"vmulpd		     (%%rcx),%%ymm7,%%ymm7		\n\t		vmulpd		0x100(%%rcx),%%ymm13,%%ymm13\n\t"\
		"vmulpd		     (%%rcx),%%ymm1,%%ymm1		\n\t		vmulpd		0x100(%%rcx),%%ymm11,%%ymm11\n\t"\
	"vfnmadd231pd	0x020(%%rcx),%%ymm4,%%ymm1		\n\t	vfnmadd231pd	0x120(%%rcx),%%ymm8 ,%%ymm11	\n\t"\
	" vfmadd231pd	0x020(%%rcx),%%ymm5,%%ymm7		\n\t	 vfmadd231pd	0x120(%%rcx),%%ymm9 ,%%ymm13	\n\t"\
		"vmovaps	%%ymm0,%%ymm4					\n\t		vmovaps		%%ymm10,%%ymm8 		\n\t"\
		"vmovaps	%%ymm6,%%ymm5					\n\t		vmovaps		%%ymm12,%%ymm9 		\n\t"\
		"vmulpd		     (%%rdx),%%ymm0,%%ymm0		\n\t		vmulpd		0x100(%%rdx),%%ymm10,%%ymm10\n\t"\
		"vmulpd		     (%%rdx),%%ymm6,%%ymm6		\n\t		vmulpd		0x100(%%rdx),%%ymm12,%%ymm12\n\t"\
		"vmovaps	%%ymm1,0x0a0(%%rbx)				\n\t		vmovaps		%%ymm11,0x0e0(%%rbx)		\n\t"\
		"vmovaps	%%ymm7,0x080(%%rbx)				\n\t		vmovaps		%%ymm13,0x0c0(%%rbx)		\n\t"\
	"vfnmadd231pd	0x020(%%rdx),%%ymm4,%%ymm6		\n\t	vfnmadd231pd	0x120(%%rdx),%%ymm8 ,%%ymm12	\n\t"\
	" vfmadd231pd	0x020(%%rdx),%%ymm5,%%ymm0		\n\t	 vfmadd231pd	0x120(%%rdx),%%ymm9 ,%%ymm10	\n\t"\
		"vmovaps	%%ymm6,0x1a0(%%rbx)				\n\t		vmovaps		%%ymm12,0x1e0(%%rbx)		\n\t"\
		"vmovaps	%%ymm0,0x180(%%rbx)				\n\t		vmovaps		%%ymm10,0x1c0(%%rbx)		\n\t"\
	"prefetcht1	-0x180(%%r14)						\n\t		movq		%[__isrt2],%%rdi		\n\t"\
		/*...Block 1: t1,9,17,25 -> r1,5,3,7: */			/*...Block 3: t5,13,21,29 -> r17,21,19,23: */\
		"movq	%[__r1],%%rax	\n\t"/* r17 in rcol */	"	vmovaps		(%%rdi),%%ymm10			\n\t"\
		"vmovaps	     (%%rax),%%ymm0				\n\t		vmovaps		0x240(%%rax),%%ymm12	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1				\n\t		vmovaps		0x260(%%rax),%%ymm13		\n\t"\
		"vmovaps	0x080(%%rax),%%ymm2				\n\t		vmovaps		0x2c0(%%rax),%%ymm8 		\n\t"\
		"vmovaps	0x0a0(%%rax),%%ymm3				\n\t		vmovaps		0x2e0(%%rax),%%ymm9 		\n\t"\
		"vsubpd		%%ymm2,%%ymm0,%%ymm0			\n\t		vaddpd		0x260(%%rax),%%ymm12,%%ymm12\n\t"\
		"vsubpd		%%ymm3,%%ymm1,%%ymm1			\n\t		vsubpd		0x240(%%rax),%%ymm13,%%ymm13\n\t"\
		"vmovaps	0x040(%%rax),%%ymm4				\n\t		vsubpd		0x2e0(%%rax),%%ymm8 ,%%ymm8 \n\t"\
		"vmovaps	0x060(%%rax),%%ymm5				\n\t		vaddpd		0x2c0(%%rax),%%ymm9 ,%%ymm9 \n\t"\
		"vmovaps	0x0c0(%%rax),%%ymm6				\n\t		vmulpd		%%ymm10,%%ymm12,%%ymm12		\n\t"\
	/*	vmovaps	0x0e0(%%rax),%%ymm7		*/		"	\n\t		vmulpd		%%ymm10,%%ymm13,%%ymm13		\n\t"\
	"vmovaps	(%%rsi),%%ymm7	\n\t"/* two */\
		"vsubpd		%%ymm6,%%ymm4,%%ymm4			\n\t		vmovaps		%%ymm12,%%ymm14		\n\t"\
		"vsubpd	0x0e0(%%rax),%%ymm5,%%ymm5			\n\t		vmovaps		%%ymm13,%%ymm15		\n\t"\
	"vfmadd132pd	%%ymm7,%%ymm0,%%ymm2			\n\t	vfnmadd231pd	%%ymm10,%%ymm8 ,%%ymm12		\n\t"\
	"vfmadd132pd	%%ymm7,%%ymm1,%%ymm3			\n\t	vfnmadd231pd	%%ymm10,%%ymm9 ,%%ymm13		\n\t"\
	"vfmadd132pd	%%ymm7,%%ymm4,%%ymm6			\n\t	 vfmadd231pd	%%ymm10,%%ymm8 ,%%ymm14		\n\t"\
	"vfmadd132pd 0x0e0(%%rax),%%ymm5,%%ymm7			\n\t	 vfmadd231pd	%%ymm10,%%ymm9 ,%%ymm15		\n\t"\
		"movq	%%rbx,%%rcx	\n\t"/* rbx saves 'high-half' block address; move that to rcx and put 'low-half' block address into rbx */\
		"movq	%[__add0],%%rbx						\n\t		vmovaps		0x200(%%rax),%%ymm8 		\n\t"\
		"subq	$0x280,%%rdx	/* c8 from c13 */	\n\t		vmovaps		0x220(%%rax),%%ymm9 		\n\t"\
		"vaddpd		%%ymm6,%%ymm2,%%ymm2			\n\t		vmovaps		0x280(%%rax),%%ymm10		\n\t"\
		"vaddpd		%%ymm7,%%ymm3,%%ymm3			\n\t		vmovaps		0x2a0(%%rax),%%ymm11		\n\t"\
		"vmovaps	%%ymm2,     (%%rbx)				\n\t		vsubpd		0x2a0(%%rax),%%ymm8 ,%%ymm8 \n\t"\
		"vmovaps	%%ymm3,0x020(%%rbx)				\n\t		vsubpd		0x280(%%rax),%%ymm9 ,%%ymm9 \n\t"\
	"vfnmadd231pd	(%%rsi),%%ymm6,%%ymm2			\n\t		vaddpd		0x200(%%rax),%%ymm11,%%ymm11\n\t"\
	"vfnmadd231pd	(%%rsi),%%ymm7,%%ymm3			\n\t		vaddpd		0x220(%%rax),%%ymm10,%%ymm10\n\t"\
		"vmovaps	%%ymm2,%%ymm6					\n\t		vsubpd		%%ymm12,%%ymm11,%%ymm11		\n\t"\
		"vmovaps	%%ymm3,%%ymm7					\n\t		vsubpd		%%ymm13,%%ymm9 ,%%ymm9 		\n\t"\
		"vmulpd		     (%%rdx),%%ymm2,%%ymm2		\n\t	vfmadd132pd	(%%rsi),%%ymm11,%%ymm12			\n\t"\
		"vmulpd		     (%%rdx),%%ymm3,%%ymm3		\n\t	vfmadd132pd	(%%rsi),%%ymm9 ,%%ymm13			\n\t"\
	"vfnmadd231pd	0x020(%%rdx),%%ymm6,%%ymm3		\n\t		vmovaps		%%ymm11,     (%%rax)		\n\t"\
	" vfmadd231pd	0x020(%%rdx),%%ymm7,%%ymm2		\n\t		vmovaps		%%ymm9 ,0x020(%%rax)		\n\t"\
		"vmovaps	0x120(%%rcx),%%ymm7				\n\t		vmovaps		%%ymm12,%%ymm11		\n\t"\
		"vmovaps	0x100(%%rcx),%%ymm6				\n\t		vmovaps		%%ymm13,%%ymm9 		\n\t"\
		"vmovaps	%%ymm3,0x120(%%rax)	/* r9 */	\n\t		vmulpd		0x0c0(%%rdx),%%ymm12,%%ymm12	\n\t"/* c2 */\
		"vmovaps	%%ymm2,0x100(%%rax)	/* r8 */	\n\t		vmulpd		0x0c0(%%rdx),%%ymm13,%%ymm13\n\t"\
		"vmovaps	%%ymm7,0x320(%%rax)	/* r25 */	\n\t	vfnmadd231pd	0x0e0(%%rdx),%%ymm11,%%ymm13\n\t"\
		"vmovaps	%%ymm6,0x300(%%rax)	/* r24 */	\n\t	 vfmadd231pd	0x0e0(%%rdx),%%ymm9 ,%%ymm12\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm3				\n\t		vmovaps	0x060(%%rcx),%%ymm11	\n\t"\
		"vmovaps	     (%%rbx),%%ymm2				\n\t		vmovaps	0x040(%%rcx),%%ymm9		\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm7				\n\t		vmovaps	%%ymm13,0x060(%%rax)	/* r3 */\n\t"\
		"vmovaps	     (%%rcx),%%ymm6				\n\t		vmovaps	%%ymm12,0x040(%%rax)	/* r2 */\n\t"\
	/* swap rcol/lcol in next 4 lines: */\
		"vmovaps		     (%%rax),%%ymm12		\n\t		vmovaps	%%ymm2 ,     (%%rax)	/* r0 */\n\t"\
		"vmovaps		0x020(%%rax),%%ymm13		\n\t		vmovaps	%%ymm3 ,0x020(%%rax)	/* r1 */\n\t"\
		"vmovaps	%%ymm11,0x260(%%rax)	/* r19*/\n\t		vmovaps		%%ymm12,%%ymm11		\n\t"\
		"vmovaps	%%ymm9 ,0x240(%%rax)	/* r18*/\n\t		vmovaps		%%ymm13,%%ymm9 		\n\t"\
		"														addq	$0x040,%%rdx		/* c4 in lcol; c10 in rcol*/		\n\t"\
		"vmovaps	%%ymm7,0x220(%%rax)	/* r17 */	\n\t		vmulpd		0x0c0(%%rdx),%%ymm12,%%ymm12		\n\t"\
		"vmovaps	%%ymm6,0x200(%%rax)	/* r16 */	\n\t		vmulpd		0x0c0(%%rdx),%%ymm13,%%ymm13		\n\t"\
		"vaddpd		%%ymm5,%%ymm0,%%ymm0			\n\t	vfnmadd231pd	0x0e0(%%rdx),%%ymm11,%%ymm13		\n\t"\
		"vsubpd		%%ymm4,%%ymm1,%%ymm1			\n\t	 vfmadd231pd	0x0e0(%%rdx),%%ymm9 ,%%ymm12		\n\t"\
		"vmovaps	%%ymm0,%%ymm2					\n\t		vmovaps	0x160(%%rcx),%%ymm11	\n\t"\
		"vmovaps	%%ymm1,%%ymm3					\n\t		vmovaps	0x140(%%rcx),%%ymm9		\n\t"\
		"			\n\t		vmovaps	%%ymm13,0x160(%%rax)	/* r11 */\n\t"\
		"			\n\t		vmovaps	%%ymm12,0x140(%%rax)	/* r10 */\n\t"\
		"vmovaps	%%ymm0,%%ymm6					\n\t		vmovaps	%%ymm11,0x360(%%rax)	/* r27 */\n\t"\
		"vmovaps	%%ymm1,%%ymm7					\n\t		vmovaps	%%ymm9 ,0x340(%%rax)	/* r26 */\n\t"\
	"vfnmadd231pd	(%%rsi),%%ymm5,%%ymm0			\n\t		vsubpd		%%ymm15,%%ymm8 ,%%ymm8 		\n\t"\
	" vfmadd231pd	(%%rsi),%%ymm4,%%ymm1			\n\t		vsubpd		%%ymm14,%%ymm10,%%ymm10		\n\t"\
		"vmulpd		     (%%rdx),%%ymm2,%%ymm2		\n\t	vfmadd132pd	(%%rsi),%%ymm8 ,%%ymm15		\n\t"\
		"vmulpd		     (%%rdx),%%ymm3,%%ymm3		\n\t	vfmadd132pd	(%%rsi),%%ymm10,%%ymm14		\n\t"\
	"vfnmadd231pd	0x020(%%rdx),%%ymm6,%%ymm3		\n\t		vmovaps		%%ymm15,%%ymm12		\n\t"\
	" vfmadd231pd	0x020(%%rdx),%%ymm7,%%ymm2		\n\t		vmovaps		%%ymm10,%%ymm13		\n\t"\
	"prefetcht1	-0x1c0(%%r14)						\n\t	addq	$0x040,%%rdx	\n\t"/* c12 in lcol; c6 in rcol*/\
		"vmovaps	0x0a0(%%rcx),%%ymm7				\n\t		vmulpd		0x0c0(%%rdx),%%ymm15,%%ymm15	\n\t"\
		"vmovaps	0x080(%%rcx),%%ymm6				\n\t		vmulpd		0x0c0(%%rdx),%%ymm10,%%ymm10	\n\t"\
		"vmovaps	%%ymm3,0x0a0(%%rax)	/* r5 */	\n\t	vfnmadd231pd	0x0e0(%%rdx),%%ymm12,%%ymm10	\n\t"\
		"vmovaps	%%ymm2,0x080(%%rax)	/* r4 */	\n\t	 vfmadd231pd	0x0e0(%%rdx),%%ymm13,%%ymm15	\n\t"\
		"vmovaps	%%ymm7,0x2a0(%%rax)	/* r21 */	\n\t		vmovaps	0x0e0(%%rcx),%%ymm13	\n\t"\
		"vmovaps	%%ymm6,0x280(%%rax)	/* r20 */	\n\t		vmovaps	0x0c0(%%rcx),%%ymm12	\n\t"\
		"vmovaps	%%ymm0,%%ymm6					\n\t		vmovaps	%%ymm10,0x0e0(%%rax)	/* r7 */\n\t"\
		"vmovaps	%%ymm1,%%ymm7					\n\t		vmovaps	%%ymm15,0x0c0(%%rax)	/* r6 */\n\t"\
		"vmulpd		     (%%rdx),%%ymm0,%%ymm0		\n\t		vmovaps	%%ymm13,0x2e0(%%rax)	/* r23 */\n\t"\
		"vmulpd		     (%%rdx),%%ymm1,%%ymm1		\n\t		vmovaps	%%ymm12,0x2c0(%%rax)	/* r22 */\n\t"\
	"vfnmadd231pd	0x020(%%rdx),%%ymm6,%%ymm1		\n\t		vmovaps		%%ymm8 ,%%ymm12		\n\t"\
	" vfmadd231pd	0x020(%%rdx),%%ymm7,%%ymm0		\n\t		vmovaps		%%ymm14,%%ymm13		\n\t"\
		"														vmulpd		0x100(%%rdx),%%ymm8 ,%%ymm8 \n\t"/* c14 */\
		"														vmulpd		0x100(%%rdx),%%ymm14,%%ymm14		\n\t"\
		"													vfnmadd231pd	0x120(%%rdx),%%ymm12,%%ymm14		\n\t"\
		"													 vfmadd231pd	0x120(%%rdx),%%ymm13,%%ymm8 		\n\t"\
		"vmovaps	0x1a0(%%rcx),%%ymm7				\n\t		vmovaps	0x1e0(%%rcx),%%ymm13	\n\t"\
		"vmovaps	0x180(%%rcx),%%ymm6				\n\t		vmovaps	0x1c0(%%rcx),%%ymm12	\n\t"\
		"vmovaps	%%ymm1,0x1a0(%%rax)	/* r13 */	\n\t		vmovaps	%%ymm14,0x1e0(%%rax)	/* r15 */\n\t"\
		"vmovaps	%%ymm0,0x180(%%rax)	/* r12 */	\n\t		vmovaps	%%ymm8 ,0x1c0(%%rax)	/* r14 */\n\t"\
		"vmovaps	%%ymm7,0x3a0(%%rax)	/* r29 */	\n\t		vmovaps	%%ymm13,0x3e0(%%rax)	/* r31 */\n\t"\
		"vmovaps	%%ymm6,0x380(%%rax)	/* r28 */	\n\t		vmovaps	%%ymm12,0x3c0(%%rax)	/* r30 */\n\t"\
	"prefetcht1	-0x200(%%r14)\n\t"\
	/**** Finish with 4-way 'un'terleaving: ****/\
		"movq	%[__r1] ,%%rsi	\n\t"\
		"movq	%[__add0],%%rax\n\t"\
		"movq	%[__add1],%%rbx\n\t"\
		"movq	%[__add2],%%rcx\n\t"\
		"movq	%[__add3],%%rdx\n\t"\
	"/* a[j+p0]: Inputs from r1 +0/1, 16/17, 2/3, 18/19. Outputs into [add0,add0+0x100,add1,add1+0x100]+0x0: */	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps	     (%%rsi),%%ymm1						\n\t		vmovaps	0x020(%%rsi),%%ymm3								\n\t"\
		"vmovaps	0x200(%%rsi),%%ymm5						\n\t		vmovaps	0x220(%%rsi),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps	0x040(%%rsi),%%ymm6						\n\t		vmovaps	0x060(%%rsi),%%ymm14							\n\t"\
		"vmovaps	0x240(%%rsi),%%ymm5						\n\t		vmovaps	0x260(%%rsi),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm6 ,%%ymm0				\n\t		vshufpd	$15 ,%%ymm13,%%ymm14 ,%%ymm2					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm6 ,%%ymm6				\n\t		vshufpd	$0  ,%%ymm13,%%ymm14 ,%%ymm14					\n\t"\
		"vperm2f128 $32 ,%%ymm0 ,%%ymm7 ,%%ymm5				\n\t		vperm2f128 $32 ,%%ymm2 ,%%ymm15,%%ymm13		/* outB	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm0 ,%%ymm7 ,%%ymm7				\n\t		vperm2f128 $49 ,%%ymm2 ,%%ymm15,%%ymm15		/* outD	*/	\n\t"\
		"vperm2f128 $32 ,%%ymm6 ,%%ymm1 ,%%ymm0				\n\t		vperm2f128 $32 ,%%ymm14,%%ymm3 ,%%ymm2 		/* outA	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm6 ,%%ymm1 ,%%ymm1				\n\t		vperm2f128 $49 ,%%ymm14,%%ymm3 ,%%ymm3 		/* outC	*/	\n\t"\
		"vmovaps 	%%ymm5 ,(%%rbx)							\n\t		vmovaps %%ymm13,0x020(%%rbx)				/* outB	*/	\n\t"\
		"vmovaps 	%%ymm7 ,(%%rdx)							\n\t		vmovaps %%ymm15,0x020(%%rdx)				/* outD	*/	\n\t"\
		"vmovaps 	%%ymm0 ,(%%rax)							\n\t		vmovaps %%ymm2 ,0x020(%%rax)				/* outA	*/	\n\t"\
		"vmovaps 	%%ymm1 ,(%%rcx)							\n\t		vmovaps %%ymm3 ,0x020(%%rcx)				/* outC	*/	\n\t"\
	"/* a[j+p2]: Inputs from r3 +0/1, 16/17, 2/3, 18/19. Outputs into [add0,add0+0x100,add1,add1+0x100]+0x40: */		\n\t"\
		"addq	$0x80,%%rsi	\n\t"\
		"addq	$0x40,%%rax	\n\t"\
		"addq	$0x40,%%rbx	\n\t"\
		"addq	$0x40,%%rcx	\n\t"\
		"addq	$0x40,%%rdx	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps	     (%%rsi),%%ymm1						\n\t		vmovaps	0x020(%%rsi),%%ymm3								\n\t"\
		"vmovaps	0x200(%%rsi),%%ymm5						\n\t		vmovaps	0x220(%%rsi),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps	0x040(%%rsi),%%ymm6						\n\t		vmovaps	0x060(%%rsi),%%ymm14							\n\t"\
		"vmovaps	0x240(%%rsi),%%ymm5						\n\t		vmovaps	0x260(%%rsi),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm6 ,%%ymm0				\n\t		vshufpd	$15 ,%%ymm13,%%ymm14 ,%%ymm2					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm6 ,%%ymm6				\n\t		vshufpd	$0  ,%%ymm13,%%ymm14 ,%%ymm14					\n\t"\
		"vperm2f128 $32 ,%%ymm0 ,%%ymm7 ,%%ymm5				\n\t		vperm2f128 $32 ,%%ymm2 ,%%ymm15,%%ymm13		/* outB	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm0 ,%%ymm7 ,%%ymm7				\n\t		vperm2f128 $49 ,%%ymm2 ,%%ymm15,%%ymm15		/* outD	*/	\n\t"\
		"vperm2f128 $32 ,%%ymm6 ,%%ymm1 ,%%ymm0				\n\t		vperm2f128 $32 ,%%ymm14,%%ymm3 ,%%ymm2 		/* outA	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm6 ,%%ymm1 ,%%ymm1				\n\t		vperm2f128 $49 ,%%ymm14,%%ymm3 ,%%ymm3 		/* outC	*/	\n\t"\
		"vmovaps 	%%ymm5 ,(%%rbx)							\n\t		vmovaps %%ymm13,0x020(%%rbx)				/* outB	*/	\n\t"\
		"vmovaps 	%%ymm7 ,(%%rdx)							\n\t		vmovaps %%ymm15,0x020(%%rdx)				/* outD	*/	\n\t"\
		"vmovaps 	%%ymm0 ,(%%rax)							\n\t		vmovaps %%ymm2 ,0x020(%%rax)				/* outA	*/	\n\t"\
		"vmovaps 	%%ymm1 ,(%%rcx)							\n\t		vmovaps %%ymm3 ,0x020(%%rcx)				/* outC	*/	\n\t"\
	"/* a[j+p4]: Inputs from r5 +0/1, 16/17, 2/3, 18/19. Outputs into [add0,add0+0x100,add1,add1+0x100]+0x80: */		\n\t"\
		"addq	$0x80,%%rsi	\n\t"\
		"addq	$0x40,%%rax	\n\t"\
		"addq	$0x40,%%rbx	\n\t"\
		"addq	$0x40,%%rcx	\n\t"\
		"addq	$0x40,%%rdx	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps	     (%%rsi),%%ymm1						\n\t		vmovaps	0x020(%%rsi),%%ymm3								\n\t"\
		"vmovaps	0x200(%%rsi),%%ymm5						\n\t		vmovaps	0x220(%%rsi),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps	0x040(%%rsi),%%ymm6						\n\t		vmovaps	0x060(%%rsi),%%ymm14							\n\t"\
		"vmovaps	0x240(%%rsi),%%ymm5						\n\t		vmovaps	0x260(%%rsi),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm6 ,%%ymm0				\n\t		vshufpd	$15 ,%%ymm13,%%ymm14 ,%%ymm2					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm6 ,%%ymm6				\n\t		vshufpd	$0  ,%%ymm13,%%ymm14 ,%%ymm14					\n\t"\
		"vperm2f128 $32 ,%%ymm0 ,%%ymm7 ,%%ymm5				\n\t		vperm2f128 $32 ,%%ymm2 ,%%ymm15,%%ymm13		/* outB	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm0 ,%%ymm7 ,%%ymm7				\n\t		vperm2f128 $49 ,%%ymm2 ,%%ymm15,%%ymm15		/* outD	*/	\n\t"\
		"vperm2f128 $32 ,%%ymm6 ,%%ymm1 ,%%ymm0				\n\t		vperm2f128 $32 ,%%ymm14,%%ymm3 ,%%ymm2 		/* outA	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm6 ,%%ymm1 ,%%ymm1				\n\t		vperm2f128 $49 ,%%ymm14,%%ymm3 ,%%ymm3 		/* outC	*/	\n\t"\
		"vmovaps 	%%ymm5 ,(%%rbx)							\n\t		vmovaps %%ymm13,0x020(%%rbx)				/* outB	*/	\n\t"\
		"vmovaps 	%%ymm7 ,(%%rdx)							\n\t		vmovaps %%ymm15,0x020(%%rdx)				/* outD	*/	\n\t"\
		"vmovaps 	%%ymm0 ,(%%rax)							\n\t		vmovaps %%ymm2 ,0x020(%%rax)				/* outA	*/	\n\t"\
		"vmovaps 	%%ymm1 ,(%%rcx)							\n\t		vmovaps %%ymm3 ,0x020(%%rcx)				/* outC	*/	\n\t"\
	"/* a[j+p6]: Inputs from r7 +0/1, 16/17, 2/3, 18/19. Outputs into [add0,add0+0x100,add1,add1+0x100]+0xc0: */		\n\t"\
		"addq	$0x80,%%rsi	\n\t"\
		"addq	$0x40,%%rax	\n\t"\
		"addq	$0x40,%%rbx	\n\t"\
		"addq	$0x40,%%rcx	\n\t"\
		"addq	$0x40,%%rdx	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps	     (%%rsi),%%ymm1						\n\t		vmovaps	0x020(%%rsi),%%ymm3								\n\t"\
		"vmovaps	0x200(%%rsi),%%ymm5						\n\t		vmovaps	0x220(%%rsi),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps	0x040(%%rsi),%%ymm6						\n\t		vmovaps	0x060(%%rsi),%%ymm14							\n\t"\
		"vmovaps	0x240(%%rsi),%%ymm5						\n\t		vmovaps	0x260(%%rsi),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm6 ,%%ymm0				\n\t		vshufpd	$15 ,%%ymm13,%%ymm14 ,%%ymm2					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm6 ,%%ymm6				\n\t		vshufpd	$0  ,%%ymm13,%%ymm14 ,%%ymm14					\n\t"\
		"vperm2f128 $32 ,%%ymm0 ,%%ymm7 ,%%ymm5				\n\t		vperm2f128 $32 ,%%ymm2 ,%%ymm15,%%ymm13		/* outB	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm0 ,%%ymm7 ,%%ymm7				\n\t		vperm2f128 $49 ,%%ymm2 ,%%ymm15,%%ymm15		/* outD	*/	\n\t"\
		"vperm2f128 $32 ,%%ymm6 ,%%ymm1 ,%%ymm0				\n\t		vperm2f128 $32 ,%%ymm14,%%ymm3 ,%%ymm2 		/* outA	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm6 ,%%ymm1 ,%%ymm1				\n\t		vperm2f128 $49 ,%%ymm14,%%ymm3 ,%%ymm3 		/* outC	*/	\n\t"\
		"vmovaps 	%%ymm5 ,(%%rbx)							\n\t		vmovaps %%ymm13,0x020(%%rbx)				/* outB	*/	\n\t"\
		"vmovaps 	%%ymm7 ,(%%rdx)							\n\t		vmovaps %%ymm15,0x020(%%rdx)				/* outD	*/	\n\t"\
		"vmovaps 	%%ymm0 ,(%%rax)							\n\t		vmovaps %%ymm2 ,0x020(%%rax)				/* outA	*/	\n\t"\
		"vmovaps 	%%ymm1 ,(%%rcx)							\n\t		vmovaps %%ymm3 ,0x020(%%rcx)				/* outC	*/	\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__add1] "m" (Xadd1)\
		 ,[__add2] "m" (Xadd2)\
		 ,[__add3] "m" (Xadd3)\
		 ,[__r1] "m" (Xr1)\
		 ,[__r9] "m" (Xr9)\
		 ,[__r17] "m" (Xr17)\
		 ,[__r25] "m" (Xr25)\
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__cc0] "m" (Xcc0)\
		 ,[__c1] "m" (Xc1)\
		 ,[__c2] "m" (Xc2)\
		 ,[__c3] "m" (Xc3)\
		 ,[__c4] "m" (Xc4)\
		 ,[__c5] "m" (Xc5)\
		 ,[__c6] "m" (Xc6)\
		 ,[__c7] "m" (Xc7)\
		 ,[__c8] "m" (Xc8)\
		 ,[__c9] "m" (Xc9)\
		 ,[__c10] "m" (Xc10)\
		 ,[__c11] "m" (Xc11)\
		 ,[__c12] "m" (Xc12)\
		 ,[__c13] "m" (Xc13)\
		 ,[__c14] "m" (Xc14)\
		 ,[__c15] "m" (Xc15)\
	 ,[__pfetch_dist] "m" (Xpfetch_dist)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r10","r11","r12","r13","r14","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
	);\
	}

#elif defined(USE_AVX)

	// Since the add0/1 block addresses advance in opposite directions, our prefetch scheme here is like so:
	// In fwd DIF prefetch  ahead 8 cache lines' worth of data w.r.to add0;
	// In inv DIT prefetch behind 8 cache lines' worth of data w.r.to add1;
	#define SSE2_RADIX16_WRAPPER_DIF(Xadd0,Xadd1,Xadd2,Xadd3,Xr1,Xr9,Xr17,Xr25,Xisrt2,Xcc0,Xc1,Xc2,Xc3,Xc4,Xc5,Xc6,Xc7,Xc8,Xc9,Xc10,Xc11,Xc12,Xc13,Xc14,Xc15,Xpfetch_dist)\
	{\
	__asm__ volatile (\
	"movq	%[__add0],%%rax						\n\t"\
	"movslq	%[__pfetch_dist],%%r14	\n\t"\
	"leaq	(%%rax,%%r14,8),%%r14	\n\t"	/* Block 0 [base-address + data-fetch-ahead index] */\
	"prefetcht1	(%%r14)\n\t"\
	/*************************************************************/\
	/*                  1st set of inputs:                       */\
	/*************************************************************/\
		"movq	%[__r1] ,%%rsi	\n\t"\
		"movq	%[__add1],%%rbx\n\t"\
		"movq	%[__add2],%%rcx\n\t"\
		"movq	%[__add3],%%rdx\n\t"\
	/**** Start with 4-way interleaving: ****/\
	"/* a[j+p0]: Inputs from [add0,add0+0x100,add1,add1+0x100]+0x0. Outputs into r1 +0/1, 8/9, 16/17, 24/25: */	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps	     (%%rax),%%ymm1						\n\t		vmovaps	0x020(%%rax),%%ymm3								\n\t"\
		"vmovaps	     (%%rbx),%%ymm5						\n\t		vmovaps	0x020(%%rbx),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps	     (%%rcx),%%ymm6						\n\t		vmovaps	0x020(%%rcx),%%ymm14							\n\t"\
		"vmovaps	     (%%rdx),%%ymm5						\n\t		vmovaps	0x020(%%rdx),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm6 ,%%ymm0				\n\t		vshufpd	$15 ,%%ymm13,%%ymm14 ,%%ymm2					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm6 ,%%ymm6				\n\t		vshufpd	$0  ,%%ymm13,%%ymm14 ,%%ymm14					\n\t"\
		"vperm2f128 $32 ,%%ymm0 ,%%ymm7 ,%%ymm5				\n\t		vperm2f128 $32 ,%%ymm2 ,%%ymm15,%%ymm13		/* outB	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm0 ,%%ymm7 ,%%ymm7				\n\t		vperm2f128 $49 ,%%ymm2 ,%%ymm15,%%ymm15		/* outD	*/	\n\t"\
		"vperm2f128 $32 ,%%ymm6 ,%%ymm1 ,%%ymm0				\n\t		vperm2f128 $32 ,%%ymm14,%%ymm3 ,%%ymm2 		/* outA	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm6 ,%%ymm1 ,%%ymm1				\n\t		vperm2f128 $49 ,%%ymm14,%%ymm3 ,%%ymm3 		/* outC	*/	\n\t"\
		"vmovaps 	%%ymm5 ,0x200(%%rsi)					\n\t		vmovaps %%ymm13,0x220(%%rsi)				/* outB	*/	\n\t"\
		"vmovaps 	%%ymm7 ,0x300(%%rsi)					\n\t		vmovaps %%ymm15,0x320(%%rsi)				/* outD	*/	\n\t"\
		"vmovaps 	%%ymm0 ,     (%%rsi)					\n\t		vmovaps %%ymm2 ,0x020(%%rsi)				/* outA	*/	\n\t"\
		"vmovaps 	%%ymm1 ,0x100(%%rsi)					\n\t		vmovaps %%ymm3 ,0x120(%%rsi)				/* outC	*/	\n\t"\
	"/* a[j+p2]: Inputs from [add0,add0+0x100,add1,add1+0x100]+0x40. Outputs into r3 +0/1, 8/9, 16/17, 24/25: */		\n\t"\
		"addq	$0x40,%%rax	\n\t"\
		"addq	$0x40,%%rbx	\n\t"\
		"addq	$0x40,%%rcx	\n\t"\
		"addq	$0x40,%%rdx	\n\t"\
		"addq	$0x40,%%rsi	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps	     (%%rax),%%ymm1						\n\t		vmovaps	0x020(%%rax),%%ymm3								\n\t"\
		"vmovaps	     (%%rbx),%%ymm5						\n\t		vmovaps	0x020(%%rbx),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps	     (%%rcx),%%ymm6						\n\t		vmovaps	0x020(%%rcx),%%ymm14							\n\t"\
		"vmovaps	     (%%rdx),%%ymm5						\n\t		vmovaps	0x020(%%rdx),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm6 ,%%ymm0				\n\t		vshufpd	$15 ,%%ymm13,%%ymm14 ,%%ymm2					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm6 ,%%ymm6				\n\t		vshufpd	$0  ,%%ymm13,%%ymm14 ,%%ymm14					\n\t"\
		"vperm2f128 $32 ,%%ymm0 ,%%ymm7 ,%%ymm5				\n\t		vperm2f128 $32 ,%%ymm2 ,%%ymm15,%%ymm13		/* outB	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm0 ,%%ymm7 ,%%ymm7				\n\t		vperm2f128 $49 ,%%ymm2 ,%%ymm15,%%ymm15		/* outD	*/	\n\t"\
		"vperm2f128 $32 ,%%ymm6 ,%%ymm1 ,%%ymm0				\n\t		vperm2f128 $32 ,%%ymm14,%%ymm3 ,%%ymm2 		/* outA	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm6 ,%%ymm1 ,%%ymm1				\n\t		vperm2f128 $49 ,%%ymm14,%%ymm3 ,%%ymm3 		/* outC	*/	\n\t"\
		"vmovaps 	%%ymm5 ,0x200(%%rsi)					\n\t		vmovaps %%ymm13,0x220(%%rsi)				/* outB	*/	\n\t"\
		"vmovaps 	%%ymm7 ,0x300(%%rsi)					\n\t		vmovaps %%ymm15,0x320(%%rsi)				/* outD	*/	\n\t"\
		"vmovaps 	%%ymm0 ,     (%%rsi)					\n\t		vmovaps %%ymm2 ,0x020(%%rsi)				/* outA	*/	\n\t"\
		"vmovaps 	%%ymm1 ,0x100(%%rsi)					\n\t		vmovaps %%ymm3 ,0x120(%%rsi)				/* outC	*/	\n\t"\
	"prefetcht1	0x40(%%r14)\n\t"\
	"/* a[j+p4]: Inputs from [add0,add0+0x100,add1,add1+0x100]+0x80. Outputs into r5 +0/1, 8/9, 16/17, 24/25: */		\n\t"\
		"addq	$0x40,%%rax	\n\t"\
		"addq	$0x40,%%rbx	\n\t"\
		"addq	$0x40,%%rcx	\n\t"\
		"addq	$0x40,%%rdx	\n\t"\
		"addq	$0x40,%%rsi	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps	     (%%rax),%%ymm1						\n\t		vmovaps	0x020(%%rax),%%ymm3								\n\t"\
		"vmovaps	     (%%rbx),%%ymm5						\n\t		vmovaps	0x020(%%rbx),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps	     (%%rcx),%%ymm6						\n\t		vmovaps	0x020(%%rcx),%%ymm14							\n\t"\
		"vmovaps	     (%%rdx),%%ymm5						\n\t		vmovaps	0x020(%%rdx),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm6 ,%%ymm0				\n\t		vshufpd	$15 ,%%ymm13,%%ymm14 ,%%ymm2					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm6 ,%%ymm6				\n\t		vshufpd	$0  ,%%ymm13,%%ymm14 ,%%ymm14					\n\t"\
		"vperm2f128 $32 ,%%ymm0 ,%%ymm7 ,%%ymm5				\n\t		vperm2f128 $32 ,%%ymm2 ,%%ymm15,%%ymm13		/* outB	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm0 ,%%ymm7 ,%%ymm7				\n\t		vperm2f128 $49 ,%%ymm2 ,%%ymm15,%%ymm15		/* outD	*/	\n\t"\
		"vperm2f128 $32 ,%%ymm6 ,%%ymm1 ,%%ymm0				\n\t		vperm2f128 $32 ,%%ymm14,%%ymm3 ,%%ymm2 		/* outA	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm6 ,%%ymm1 ,%%ymm1				\n\t		vperm2f128 $49 ,%%ymm14,%%ymm3 ,%%ymm3 		/* outC	*/	\n\t"\
		"vmovaps 	%%ymm5 ,0x200(%%rsi)					\n\t		vmovaps %%ymm13,0x220(%%rsi)				/* outB	*/	\n\t"\
		"vmovaps 	%%ymm7 ,0x300(%%rsi)					\n\t		vmovaps %%ymm15,0x320(%%rsi)				/* outD	*/	\n\t"\
		"vmovaps 	%%ymm0 ,     (%%rsi)					\n\t		vmovaps %%ymm2 ,0x020(%%rsi)				/* outA	*/	\n\t"\
		"vmovaps 	%%ymm1 ,0x100(%%rsi)					\n\t		vmovaps %%ymm3 ,0x120(%%rsi)				/* outC	*/	\n\t"\
	"/* a[j+p6]: Inputs from [add0,add0+0x100,add1,add1+0x100]+0xc0. Outputs into r7 +0/1, 8/9, 16/17, 24/25: */		\n\t"\
		"addq	$0x40,%%rax	\n\t"\
		"addq	$0x40,%%rbx	\n\t"\
		"addq	$0x40,%%rcx	\n\t"\
		"addq	$0x40,%%rdx	\n\t"\
		"addq	$0x40,%%rsi	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps	     (%%rax),%%ymm1						\n\t		vmovaps	0x020(%%rax),%%ymm3								\n\t"\
		"vmovaps	     (%%rbx),%%ymm5						\n\t		vmovaps	0x020(%%rbx),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps	     (%%rcx),%%ymm6						\n\t		vmovaps	0x020(%%rcx),%%ymm14							\n\t"\
		"vmovaps	     (%%rdx),%%ymm5						\n\t		vmovaps	0x020(%%rdx),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm6 ,%%ymm0				\n\t		vshufpd	$15 ,%%ymm13,%%ymm14 ,%%ymm2					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm6 ,%%ymm6				\n\t		vshufpd	$0  ,%%ymm13,%%ymm14 ,%%ymm14					\n\t"\
		"vperm2f128 $32 ,%%ymm0 ,%%ymm7 ,%%ymm5				\n\t		vperm2f128 $32 ,%%ymm2 ,%%ymm15,%%ymm13		/* outB	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm0 ,%%ymm7 ,%%ymm7				\n\t		vperm2f128 $49 ,%%ymm2 ,%%ymm15,%%ymm15		/* outD	*/	\n\t"\
		"vperm2f128 $32 ,%%ymm6 ,%%ymm1 ,%%ymm0				\n\t		vperm2f128 $32 ,%%ymm14,%%ymm3 ,%%ymm2 		/* outA	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm6 ,%%ymm1 ,%%ymm1				\n\t		vperm2f128 $49 ,%%ymm14,%%ymm3 ,%%ymm3 		/* outC	*/	\n\t"\
		"vmovaps 	%%ymm5 ,0x200(%%rsi)					\n\t		vmovaps %%ymm13,0x220(%%rsi)				/* outB	*/	\n\t"\
		"vmovaps 	%%ymm7 ,0x300(%%rsi)					\n\t		vmovaps %%ymm15,0x320(%%rsi)				/* outD	*/	\n\t"\
		"vmovaps 	%%ymm0 ,     (%%rsi)					\n\t		vmovaps %%ymm2 ,0x020(%%rsi)				/* outA	*/	\n\t"\
		"vmovaps 	%%ymm1 ,0x100(%%rsi)					\n\t		vmovaps %%ymm3 ,0x120(%%rsi)				/* outC	*/	\n\t"\
	/*****************/\
	/* Radix-16 DIF: */\
	/*****************/\
	"prefetcht1	0x80(%%r14)\n\t"\
	"movq	%[__add0],%%rax	\n\t"/* Use for FMA-related spills */\
		"movq	%[__r1] ,%%rcx	\n\t"\
	/*...Block 1: */													/*...Block 2: */\
		"leaq	0x560(%%rcx),%%rdi	/* c2 */\n\t"\
		"leaq	0x4e0(%%rcx),%%rdx	/* c4 */\n\t"\
		"vmovaps	0x040(%%rcx),%%ymm0	/* ymm0 <-     a[jt+p4] */			\n\t		vmovaps		0x100(%%rcx),%%ymm8	/* ymm10 <-     a[jt+p2] */			\n\t"\
		"vmovaps	0x060(%%rcx),%%ymm1	/* ymm1 <-     a[jp+p4] */			\n\t		vmovaps		0x120(%%rcx),%%ymm9	/* ymm11 <-     a[jp+p2] */			\n\t"\
		"vmovaps	%%ymm0		,%%ymm2	/* ymm2 <- cpy a[jt+p4] */			\n\t		vmovaps		%%ymm8 	,%%ymm10	/* ymm10 <- cpy a[jt+p2] */			\n\t"\
		"vmovaps	%%ymm1		,%%ymm3	/* ymm3 <- cpy a[jp+p4] */			\n\t		vmovaps		%%ymm9 	,%%ymm11	/* ymm11 <- cpy a[jp+p2] */			\n\t"\
		/***************************************************************************/\
		/*** From hereon, things are identical to the code in radix16_dif_pass: ****/\
		/***************************************************************************/\
		"vmulpd		(%%rdx)		,%%ymm0,%%ymm0		/* a[jt+p4]*c4 */		\n\t		vmulpd		     (%%rdi)	,%%ymm8 ,%%ymm8 		/* a[jt+p2]*c2 */	\n\t"\
		"vmulpd		(%%rdx)		,%%ymm1,%%ymm1		/* a[jp+p4]*c4 */		\n\t		vmulpd		     (%%rdi)	,%%ymm9 ,%%ymm9 		/* a[jp+p2]*c2 */	\n\t"\
		"vmulpd		0x020(%%rdx),%%ymm2,%%ymm2		/* a[jt+p4]*s4 */		\n\t		vmulpd		0x020(%%rdi),%%ymm10,%%ymm10		/* a[jt+p2]*s2 */	\n\t"\
		"vmulpd		0x020(%%rdx),%%ymm3,%%ymm3		/* a[jp+p4]*s4 */		\n\t		vmulpd		0x020(%%rdi),%%ymm11,%%ymm11		/* a[jp+p2]*s2 */	\n\t"\
		"vaddpd		%%ymm2		,%%ymm1,%%ymm1	/* ymm1 <- t6 */			\n\t		vaddpd		%%ymm10		,%%ymm9 ,%%ymm9 	/* ymm9  <- t10*/		\n\t"\
		"vsubpd		%%ymm3		,%%ymm0,%%ymm0	/* ymm0 <- t5 */			\n\t		vsubpd		%%ymm11		,%%ymm8 ,%%ymm8 	/* ymm8  <- t9 */		\n\t"\
		"vmovaps	%%ymm1		,%%ymm3		/* ymm3 <- cpy t6 */			\n\t		vmovaps		%%ymm9 		,%%ymm11		/* ymm11 <- cpy t10*/		\n\t"\
		"vmovaps	%%ymm0		,%%ymm2		/* ymm2 <- cpy t5 */			\n\t		vmovaps		%%ymm8 		,%%ymm10		/* ymm10 <- cpy t9 */		\n\t"\
		"addq		$0x040,%%rdx	/* c12 */								\n\t		addq	$0x040,%%rdi	/* c10 */		\n\t"\
		"vmovaps	0x0c0(%%rcx),%%ymm4			/* ymm4 <-     a[jt+p12] */	\n\t		vmovaps		0x180(%%rcx),%%ymm12			/* ymm12 <-     a[jt+p10] */\n\t"\
		"vmovaps	0x0e0(%%rcx),%%ymm5			/* ymm5 <-     a[jp+p12] */	\n\t		vmovaps		0x1a0(%%rcx),%%ymm13			/* ymm13 <-     a[jp+p10] */\n\t"\
		"vmovaps	%%ymm4		,%%ymm6			/* ymm6 <- cpy a[jt+p12] */	\n\t		vmovaps			%%ymm12	,%%ymm14			/* ymm14 <- cpy a[jt+p10] */\n\t"\
		"vmovaps	%%ymm5		,%%ymm7			/* ymm7 <- cpy a[jp+p12] */	\n\t		vmovaps			%%ymm13	,%%ymm15			/* ymm15 <- cpy a[jp+p10] */\n\t"\
		"vmulpd		(%%rdx)		,%%ymm4,%%ymm4		/* a[jt+p12]*c12 */		\n\t		vmulpd		     (%%rdi)	,%%ymm12,%%ymm12		/* a[jt+p10]*c10 */		\n\t"\
		"vmulpd		(%%rdx)		,%%ymm5,%%ymm5		/* a[jp+p12]*c12 */		\n\t		vmulpd		     (%%rdi)	,%%ymm13,%%ymm13		/* a[jp+p10]*c10 */		\n\t"\
		"vmulpd		0x020(%%rdx),%%ymm6,%%ymm6		/* a[jt+p12]*s12 */		\n\t		vmulpd		0x020(%%rdi),%%ymm14,%%ymm14		/* a[jt+p10]*s10 */		\n\t"\
		"vmulpd		0x020(%%rdx),%%ymm7,%%ymm7		/* a[jp+p12]*s12 */		\n\t		vmulpd		0x020(%%rdi),%%ymm15,%%ymm15		/* a[jp+p10]*s10 */		\n\t"\
		"vaddpd		%%ymm6		,%%ymm5,%%ymm5	/* ymm5 <- it */			\n\t		vaddpd		%%ymm14		,%%ymm13,%%ymm13	/* ymm13 <- it */						\n\t"\
		"vsubpd		%%ymm7		,%%ymm4,%%ymm4	/* ymm4 <- rt; ymm6,7 free*/\n\t		vsubpd		%%ymm15		,%%ymm12,%%ymm12	/* ymm12 <- rt		ymm14,7 free */		\n\t"\
		"vaddpd		%%ymm4		,%%ymm0,%%ymm0		/* ~t5 <- t5 +rt */		\n\t		vaddpd		%%ymm12		,%%ymm8 ,%%ymm8  	/* ~t13<- t13+rt */						\n\t"\
		"vaddpd		%%ymm5		,%%ymm1,%%ymm1		/* ~t6 <- t6 +it */		\n\t		vaddpd		%%ymm13		,%%ymm9 ,%%ymm9  	/* ~t14<- t14+it */						\n\t"\
		"vsubpd		%%ymm4		,%%ymm2,%%ymm2		/* ~t7 <- t5 -rt */		\n\t		vsubpd		%%ymm12		,%%ymm10,%%ymm10 	/* ~t15<- t13-rt */						\n\t"\
		"vsubpd		%%ymm5		,%%ymm3,%%ymm3		/* ~t8 <- t6 -it */		\n\t		vsubpd		%%ymm13		,%%ymm11,%%ymm11 	/* ~t16<- t14-it	ymm12,13 free */	\n\t"\
		"\n\t"\
		"/* Now do the p0,8 combo: */\n\t									/* Do the p6,14 combo - do p14 first so registers come out in same order as for p2,10 */\n\t"\
		"leaq	0x4a0(%%rcx),%%rdx	/* c8 */								\n\t		leaq	0x620(%%rcx),%%rdi	/* c14 */\n\t"\
		"vmovaps	0x080(%%rcx)	,%%ymm4		/* a[jt+p8 ] */				\n\t		vmovaps		0x1c0(%%rcx),%%ymm12		/* a[jt+p14] */				\n\t"\
		"vmovaps	0x0a0(%%rcx)	,%%ymm5		/* a[jp+p8 ] */				\n\t		vmovaps		0x1e0(%%rcx),%%ymm13		/* a[jp+p14] */				\n\t"\
		"																				vmovaps			%%ymm12	,%%ymm14		/* ymm14 <- cpy a[jt+p14] */\n\t"\
		"vmovaps	%%ymm4		,%%ymm6	/* ymm6 <- cpy a[jt+p8] */			\n\t		vmovaps			%%ymm13	,%%ymm15		/* ymm15 <- cpy a[jp+p14] */\n\t"\
		"vmovaps	%%ymm5		,%%ymm7	/* ymm7 <- cpy a[jp+p8] */			\n\t		vmulpd		     (%%rdi)	,%%ymm12,%%ymm12		/* a[jt+p14]*c14 */			\n\t"\
		"																				vmulpd		     (%%rdi)	,%%ymm13,%%ymm13		/* a[jp+p14]*c14 */			\n\t"\
		"vmulpd		(%%rdx)		,%%ymm4,%%ymm4		/* a[jt+p8]*c8 */		\n\t		vmulpd		0x020(%%rdi),%%ymm14,%%ymm14		/* a[jt+p14]*s14 */			\n\t"\
		"vmulpd		(%%rdx)		,%%ymm5,%%ymm5		/* a[jp+p8]*c8 */		\n\t		vmulpd		0x020(%%rdi),%%ymm15,%%ymm15		/* a[jp+p14]*s14 */			\n\t"\
		"vmulpd		0x020(%%rdx),%%ymm6,%%ymm6		/* a[jt+p8]*s8 */		\n\t		vaddpd		%%ymm14		,%%ymm13,%%ymm13		/* ymm13 <- it */			\n\t"\
		"vmulpd		0x020(%%rdx),%%ymm7,%%ymm7		/* a[jp+p8]*s8 */		\n\t		vsubpd		%%ymm15		,%%ymm12,%%ymm12		/* ymm12 <- rt ymm14,15 free */	\n\t"\
		"vaddpd		%%ymm6		,%%ymm5,%%ymm5	/* ymm5 <- it */			\n\t		vmovaps		%%ymm13		,0x1e0(%%rcx)	/* Store it in t16*/		\n\t"\
		"vsubpd		%%ymm7		,%%ymm4,%%ymm4	/* ymm4 <- rt; ymm6,7 free*/\n\t		vmovaps		%%ymm12		,0x1c0(%%rcx)	/* Store rt in t15*/		\n\t"\
		"																				subq	$0x040,%%rdi	/* c6  */	\n\t"\
	"prefetcht1	0xc0(%%r14)\n\t"\
	"/* Real parts:*/\n\t															/* Real parts: */\n\t"\
		"vmovaps		 (%%rcx),%%ymm6		/* a[jt    ] */					\n\t		vmovaps		0x140(%%rcx),%%ymm12		/* a[jt+p6 ] */				\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm7		/* a[jp    ] */					\n\t		vmovaps		0x160(%%rcx),%%ymm13		/* a[jp+p6 ] */				\n\t"\
		"																				vmovaps			%%ymm12	,%%ymm14		/* ymm14 <- cpy a[jt+p6] */		\n\t"\
		"																				vmovaps			%%ymm13	,%%ymm15		/* ymm15 <- cpy a[jp+p6] */		\n\t"\
		"vsubpd		%%ymm4	,%%ymm6,%%ymm6	/* ~t3 <- t1 -rt */				\n\t		vmulpd		     (%%rdi)	,%%ymm12,%%ymm12		/* a[jt+p6]*c6 */			\n\t"\
		"vsubpd		%%ymm5	,%%ymm7,%%ymm7	/* ~t4 <- t2 -it */				\n\t		vmulpd		     (%%rdi)	,%%ymm13,%%ymm13		/* a[jp+p6]*c6 */			\n\t"\
		"vaddpd		%%ymm4	,%%ymm4,%%ymm4	/*          2*rt */				\n\t		vmulpd		0x020(%%rdi),%%ymm14,%%ymm14		/* a[jt+p6]*s6 */			\n\t"\
		"vaddpd		%%ymm5	,%%ymm5,%%ymm5	/*          2*it */				\n\t		vmulpd		0x020(%%rdi),%%ymm15,%%ymm15		/* a[jp+p6]*s6 */			\n\t"\
		"vaddpd		%%ymm6	,%%ymm4,%%ymm4	/* ~t1 <- t1 +rt */				\n\t		vaddpd		%%ymm14		,%%ymm13,%%ymm13	/* ymm13 <- t14*/				\n\t"\
		"vaddpd		%%ymm7	,%%ymm5,%%ymm5	/* ~t2 <- t2 +it; ymm4,5 free */\n\t		vsubpd		%%ymm15		,%%ymm12,%%ymm12	/* ymm12 <- t13*/				\n\t"\
		"																				vmovaps		%%ymm13		,%%ymm15		/* ymm15 <- cpy t14*/			\n\t"\
		"																				vmovaps		%%ymm12		,%%ymm14		/* ymm14 <- cpy t13*/			\n\t"\
		"																				vsubpd		0x1c0(%%rcx),%%ymm12,%%ymm12		/* ~t15<- t13-rt */			\n\t"\
		"																				vsubpd		0x1e0(%%rcx),%%ymm13,%%ymm13		/* ~t16<- t14-it */			\n\t"\
		"																				vaddpd		0x1c0(%%rcx),%%ymm14,%%ymm14		/* ~t13<- t13+rt */			\n\t"\
		"																				vaddpd		0x1e0(%%rcx),%%ymm15,%%ymm15		/* ~t14<- t14+it */			\n\t"\
		"\n\t"\
	"/* Finish radix-4 butterfly and store results into temporary-array slots: */\n\t	/* Finish radix-4 butterfly and store results into temporary-array slots: */\n\t"\
		"vsubpd		%%ymm0		,%%ymm4,%%ymm4	/*~t5 =t1 -t5 */			\n\t		vsubpd		%%ymm14		,%%ymm8,%%ymm8 	/*~t13*/						\n\t"\
		"vsubpd		%%ymm1		,%%ymm5,%%ymm5	/*~t6 =t2 -t6 */			\n\t		vsubpd		%%ymm15		,%%ymm9,%%ymm9 	/*~t14*/						\n\t"\
		"vmovaps	%%ymm4		,0x080(%%rcx)	/* a[jt+p8 ] <- ~t5 */		\n\t		vmovaps		%%ymm8 		,0x180(%%rcx)	/* a[jt+p8 ] <- ~t13*/		\n\t"\
		"vmovaps	%%ymm5		,0x0a0(%%rcx)	/* a[jp+p8 ] <- ~t6 */		\n\t		vmovaps		%%ymm9 		,0x1a0(%%rcx)	/* a[jp+p8 ] <- ~t14*/		\n\t"\
		"vaddpd		%%ymm0		,%%ymm0,%%ymm0		/* 2*t5 */				\n\t		vaddpd		%%ymm14		,%%ymm14,%%ymm14	/* 2*t13*/						\n\t"\
		"vaddpd		%%ymm1		,%%ymm1,%%ymm1		/* 2*t6 */				\n\t		vaddpd		%%ymm15		,%%ymm15,%%ymm15	/* 2*t14*/						\n\t"\
		"vaddpd		%%ymm4		,%%ymm0,%%ymm0		/*~t1 =t1 +t5 */		\n\t		vaddpd		%%ymm8 		,%%ymm14,%%ymm14	/*~t9 */						\n\t"\
		"vaddpd		%%ymm5		,%%ymm1,%%ymm1		/*~t2 =t2 +t6 */		\n\t		vaddpd		%%ymm9 		,%%ymm15,%%ymm15	/*~t10*/						\n\t"\
		"vmovaps	%%ymm0		,     (%%rcx)	/* a[jt    ] <- ~t1 */		\n\t		vmovaps		%%ymm14		,0x100(%%rcx)	/* a[jt    ] <- ~t9 */		\n\t"\
		"vmovaps	%%ymm1		,0x020(%%rcx)	/* a[jp    ] <- ~t2 */		\n\t		vmovaps		%%ymm15		,0x120(%%rcx)	/* a[jp    ] <- ~t10*/		\n\t"\
		"\n\t"\
		"vsubpd		%%ymm3		,%%ymm6,%%ymm6	/*~t3 =t3 -t8 */			\n\t		vsubpd		%%ymm13		,%%ymm10,%%ymm10	/*~t11*/						\n\t"\
		"vsubpd		%%ymm2		,%%ymm7,%%ymm7	/*~t8 =t4 -t7 */			\n\t		vsubpd		%%ymm12		,%%ymm11,%%ymm11	/*~t16*/						\n\t"\
		"vmovaps	%%ymm6		,0x040(%%rcx)	/* a[jt+p4 ] <- ~t3 */		\n\t		vmovaps		%%ymm10		,0x140(%%rcx)	/* a[jt+p4 ] <- ~t11*/		\n\t"\
		"vmovaps	%%ymm7		,0x0e0(%%rcx)	/* a[jp+p12] <- ~t8 */		\n\t		vmovaps		%%ymm11		,0x1e0(%%rcx)	/* a[jp+p12] <- ~t16*/		\n\t"\
		"vaddpd		%%ymm3		,%%ymm3,%%ymm3		/* 2*t8 */				\n\t		vaddpd		%%ymm13		,%%ymm13,%%ymm13	/* 2*t16*/						\n\t"\
		"vaddpd		%%ymm2		,%%ymm2,%%ymm2		/* 2*t7 */				\n\t		vaddpd		%%ymm12		,%%ymm12,%%ymm12	/* 2*t15*/						\n\t"\
		"vaddpd		%%ymm6		,%%ymm3,%%ymm3		/*~t7 =t3 +t8 */		\n\t		vaddpd		%%ymm10		,%%ymm13,%%ymm13	/*~t15*/						\n\t"\
		"vaddpd		%%ymm7		,%%ymm2,%%ymm2		/*~t4 =t4 +t7 */		\n\t		vaddpd		%%ymm11		,%%ymm12,%%ymm12	/*~t12*/						\n\t"\
		"vmovaps	%%ymm3		,0x0c0(%%rcx)	/* a[jt+p12] <- ~t7 */		\n\t		vmovaps		%%ymm13		,0x1c0(%%rcx)	/* a[jt+p12] <- ~t15*/		\n\t"\
		"vmovaps	%%ymm2		,0x060(%%rcx)	/* a[jp+p4 ] <- ~t4 */		\n\t		vmovaps		%%ymm12		,0x160(%%rcx)	/* a[jp+p4 ] <- ~t12*/		\n\t"\
		"\n\t"\
	"/*...Block 3: */														\n\t	/*...Block 4: */\n\t"\
	"/*	SSE2_RADIX4_DIF_4TWIDDLE(r17,r21,r19,r23, r17,c1): */				\n\t	/* SSE2_RADIX4_DIF_4TWIDDLE(r25,r29,r27,r31, r25,c3): */\n\t"\
		"/* Do the p0,p8 combo: */											\n\t	/* Do the p0,p8 combo: */		\n\t"\
		"leaq	0x660(%%rcx),%%rbx	/* c1 */								\n\t/* All __r and __c pointers incr by +0x100 in rcol w.r.to lcol: */\n\t"\
		"addq	$0x200,%%rcx		/* r17 */								\n\t		/* r25 */\n\t"\
		"\n\t																\n\t		/* c3  */\n\t"\
		"vmovaps		 (%%rcx),%%ymm0		/* a[jt   ] */					\n\t		vmovaps		0x100(%%rcx),%%ymm8 		/* a[jt    ] */\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm1		/* a[jp   ] */					\n\t		vmovaps		0x120(%%rcx),%%ymm9 		/* a[jp    ] */\n\t"\
		"vmovaps	0x080(%%rcx),%%ymm4		/* a[jt+p8 ] */					\n\t		vmovaps		0x180(%%rcx),%%ymm12		/* a[jt+p8 ] */\n\t"\
		"vmovaps	0x0a0(%%rcx),%%ymm5		/* a[jp+p8 ] */					\n\t		vmovaps		0x1a0(%%rcx),%%ymm13		/* a[jp+p8 ] */\n\t"\
		"vmovaps		 (%%rbx),%%ymm6		/* c0 */						\n\t		vmovaps		0x100(%%rbx),%%ymm14		/* c0 */\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm7		/* s0 */						\n\t		vmovaps		0x120(%%rbx),%%ymm15		/* s0 */\n\t"\
		"vmovaps	%%ymm0		,%%ymm2		/* ymm2 <- cpy a[jt   ] */		\n\t		vmovaps		%%ymm8 		,%%ymm10		/* ymm10 <- cpy a[jt   ] */\n\t"\
		"vmovaps	%%ymm1		,%%ymm3		/* ymm3 <- cpy a[jp   ] */		\n\t		vmovaps		%%ymm9 		,%%ymm11		/* ymm11 <- cpy a[jp   ] */\n\t"\
	"prefetcht1	0x100(%%r14)\n\t"\
		"\n\t																\n\t"\
		"vmulpd   	%%ymm6		,%%ymm0,%%ymm0		/* a[jt   ]*c0 */		\n\t		vmulpd   	%%ymm14		,%%ymm8 ,%%ymm8 		/* a[jt   ]*c0 */\n\t"\
		"vmulpd   	%%ymm6		,%%ymm1,%%ymm1		/* a[jp   ]*c0 */		\n\t		vmulpd   	%%ymm14		,%%ymm9 ,%%ymm9 		/* a[jp   ]*c0 */\n\t"\
		"vmulpd   	%%ymm7		,%%ymm2,%%ymm2		/* a[jt   ]*s0 */		\n\t		vmulpd   	%%ymm15		,%%ymm10,%%ymm10		/* a[jt   ]*s0 */\n\t"\
		"vmulpd   	%%ymm7		,%%ymm3,%%ymm3	/* a[jp   ]*s0 ymm6,7 free*/\n\t		vmulpd   	%%ymm15		,%%ymm11,%%ymm11		/* a[jp   ]*s0	ymm14,7 free */\n\t"\
		"vmovaps	%%ymm4		,%%ymm6		/* ymm6 <- cpy a[jt+p8 ] */		\n\t		vmovaps		%%ymm12		,%%ymm14		/* ymm14 <- cpy a[jt+p8 ] */\n\t"\
		"vmovaps	%%ymm5		,%%ymm7		/* ymm7 <- cpy a[jp+p8 ] */		\n\t		vmovaps		%%ymm13		,%%ymm15		/* ymm15 <- cpy a[jp+p8 ] */\n\t"\
		"vaddpd   	%%ymm2		,%%ymm1,%%ymm1		/* ymm1 <- t2 */		\n\t		vaddpd   	%%ymm10		,%%ymm9 ,%%ymm9 		/* ymm9  <- t2 */			\n\t"\
		"vsubpd   	%%ymm3		,%%ymm0,%%ymm0		/* ymm0 <- t1 */		\n\t		vsubpd   	%%ymm11		,%%ymm8 ,%%ymm8 		/* ymm8  <- t1 */			\n\t"\
		"vmulpd		0x040(%%rbx),%%ymm4,%%ymm4		/* a[jt+p8 ]*c8 */		\n\t		vmulpd		0x140(%%rbx),%%ymm12,%%ymm12		/* a[jt+p8 ]*c8 */\n\t"\
		"vmulpd		0x040(%%rbx),%%ymm5,%%ymm5		/* a[jp+p8 ]*c8 */		\n\t		vmulpd		0x140(%%rbx),%%ymm13,%%ymm13		/* a[jp+p8 ]*c8 */\n\t"\
		"vmulpd		0x060(%%rbx),%%ymm6,%%ymm6		/* a[jt+p8 ]*s8 */		\n\t		vmulpd		0x160(%%rbx),%%ymm14,%%ymm14		/* a[jt+p8 ]*s8 */\n\t"\
		"vmulpd		0x060(%%rbx),%%ymm7,%%ymm7		/* a[jp+p8 ]*s8 */		\n\t		vmulpd		0x160(%%rbx),%%ymm15,%%ymm15		/* a[jp+p8 ]*s8 */\n\t"\
		"vmovaps	%%ymm0		,%%ymm2		/* ymm2 <- cpy t1 */			\n\t		vmovaps		%%ymm8 		,%%ymm10		/* ymm10 <- cpy t1 */		\n\t"\
		"vmovaps	%%ymm1		,%%ymm3		/* ymm3 <- cpy t2 */			\n\t		vmovaps		%%ymm9 		,%%ymm11		/* ymm11 <- cpy t2 */		\n\t"\
		"vaddpd		%%ymm6	    ,%%ymm5,%%ymm5		/* ymm5 <- it */		\n\t		vaddpd		%%ymm14	    ,%%ymm13,%%ymm13		/* ymm13 <- it */\n\t"\
		"vsubpd		%%ymm7		,%%ymm4,%%ymm4	/* ymm4 <- rt ymm6,7 free*/	\n\t		vsubpd		%%ymm15		,%%ymm12,%%ymm12		/* ymm12 <- rt    ymm14,7 free */\n\t"\
		"addq	$0x0c0 ,%%rbx	/* c13 */									\n\t		/* c15 */\n\t"\
		"vmovaps	0x0c0(%%rcx),%%ymm6		/* a[jt+p12] */					\n\t		vmovaps		0x1c0(%%rcx),%%ymm14		/* a[jt+p12] */\n\t"\
		"vmovaps	0x0e0(%%rcx),%%ymm7		/* a[jp+p12] */					\n\t		vmovaps		0x1e0(%%rcx),%%ymm15		/* a[jp+p12] */\n\t"\
		"vaddpd		%%ymm4	,%%ymm0,%%ymm0		/* ~t1 <- t1 +rt */			\n\t		vaddpd		%%ymm12	,%%ymm8 ,%%ymm8 		/* ~t1 <- t1 +rt */\n\t"\
		"vaddpd		%%ymm5	,%%ymm1,%%ymm1		/* ~t2 <- t2 +it */			\n\t		vaddpd		%%ymm13	,%%ymm9 ,%%ymm9 		/* ~t2 <- t2 +it */\n\t"\
		"vsubpd		%%ymm4	,%%ymm2,%%ymm2		/* ~t3 <- t1 -rt */			\n\t		vsubpd		%%ymm12	,%%ymm10,%%ymm10		/* ~t3 <- t1 -rt */\n\t"\
		"vsubpd		%%ymm5	,%%ymm3,%%ymm3	/* ~t4 <- t2 -it ymm4,5 free*/	\n\t		vsubpd		%%ymm13	,%%ymm11,%%ymm11		/* ~t4 <- t2 -it	ymm12,5 free */\n\t"\
		"\n\t"\
	"/* Do the p4,12 combo: */												\n\t	/* Do the p4,12 combo: */\n\t"\
		"vmovaps	%%ymm6		,%%ymm4		/* ymm4 <- cpy a[jt+p12] */		\n\t		vmovaps		%%ymm14		,%%ymm12		/* ymm12 <- cpy a[jt+p12] */\n\t"\
		"vmovaps	%%ymm7		,%%ymm5		/* ymm5 <- cpy a[jp+p12] */		\n\t		vmovaps		%%ymm15		,%%ymm13		/* ymm13 <- cpy a[jp+p12] */\n\t"\
		"vmulpd		(%%rbx)		,%%ymm4,%%ymm4		/* a[jt+p12]*c12 */		\n\t		vmulpd		0x100(%%rbx),%%ymm12,%%ymm12		/* a[jt+p12]*c12 */\n\t"\
		"vmulpd		(%%rbx)		,%%ymm5,%%ymm5		/* a[jp+p12]*c12 */		\n\t		vmulpd		0x100(%%rbx),%%ymm13,%%ymm13		/* a[jp+p12]*c12 */\n\t"\
		"vmulpd		0x020(%%rbx),%%ymm6,%%ymm6		/* a[jt+p12]*s12 */		\n\t		vmulpd		0x120(%%rbx),%%ymm14,%%ymm14		/* a[jt+p12]*s12 */\n\t"\
		"vmulpd		0x020(%%rbx),%%ymm7,%%ymm7		/* a[jp+p12]*s12 */		\n\t		vmulpd		0x120(%%rbx),%%ymm15,%%ymm15		/* a[jp+p12]*s12 */\n\t"\
		"vaddpd		%%ymm6		,%%ymm5,%%ymm5	/* ymm5 <- it */			\n\t		vaddpd		%%ymm14		,%%ymm13,%%ymm13	/* ymm13 <- it */\n\t"\
		"vsubpd		%%ymm7		,%%ymm4,%%ymm4	/* ymm4 <- rt */			\n\t		vsubpd		%%ymm15		,%%ymm12,%%ymm12	/* ymm12 <- rt */\n\t"\
		"vmovaps	%%ymm5		,0x020(%%rcx)	/* store it */				\n\t		vmovaps		%%ymm13		,0x120(%%rcx)	/* store it */\n\t"\
		"vmovaps	%%ymm4		,     (%%rcx)	/* store rt */				\n\t		vmovaps		%%ymm12		,0x100(%%rcx)	/* store rt */\n\t"\
		"\n\t"\
		"subq	$0x040 ,%%rbx	/* c5 */									\n\t		/* c7  */\n\t"\
		"vmovaps	0x040(%%rcx),%%ymm4		/* a[jt+p4] */					\n\t		vmovaps		0x140(%%rcx),%%ymm12		/* a[jt+p4] */\n\t"\
		"vmovaps	0x060(%%rcx),%%ymm5		/* a[jp+p4] */					\n\t		vmovaps		0x160(%%rcx),%%ymm13		/* a[jp+p4] */\n\t"\
		"vmovaps		%%ymm4	,%%ymm6		/* ymm4 <- cpy a[jt+p4] */		\n\t		vmovaps			%%ymm12	,%%ymm14		/* ymm12 <- cpy a[jt+p4] */\n\t"\
		"vmovaps		%%ymm5	,%%ymm7		/* ymm5 <- cpy a[jp+p4] */		\n\t		vmovaps			%%ymm13	,%%ymm15		/* ymm13 <- cpy a[jp+p4] */\n\t"\
		"vmulpd		     (%%rbx),%%ymm4,%%ymm4		/* a[jt+p4]*c4 */		\n\t		vmulpd		0x100(%%rbx),%%ymm12,%%ymm12		/* a[jt+p4]*c4 */\n\t"\
		"vmulpd		     (%%rbx),%%ymm5,%%ymm5		/* a[jp+p4]*c4 */		\n\t		vmulpd		0x100(%%rbx),%%ymm13,%%ymm13		/* a[jp+p4]*c4 */\n\t"\
		"vmulpd		0x020(%%rbx),%%ymm6,%%ymm6		/* a[jt+p4]*s4 */		\n\t		vmulpd		0x120(%%rbx),%%ymm14,%%ymm14		/* a[jt+p4]*s4 */\n\t"\
		"vmulpd		0x020(%%rbx),%%ymm7,%%ymm7		/* a[jp+p4]*s4 */		\n\t		vmulpd		0x120(%%rbx),%%ymm15,%%ymm15		/* a[jp+p4]*s4 */\n\t"\
		"vaddpd		%%ymm6		,%%ymm5,%%ymm5	/* ymm5 <- t6 */			\n\t		vaddpd		%%ymm14		,%%ymm13,%%ymm13	/* ymm13 <- t6 */\n\t"\
		"vsubpd		%%ymm7		,%%ymm4,%%ymm4	/* ymm4 <- t5 ymm6,7 free */\n\t		vsubpd		%%ymm15		,%%ymm12,%%ymm12	/* ymm12 <- t5 	ymm14,7 free */\n\t"\
		"vmovaps	%%ymm5		,%%ymm7		/* ymm7 <- cpy t6 */			\n\t		vmovaps		%%ymm13		,%%ymm15		/* ymm15 <- cpy t6 */\n\t"\
		"vmovaps	%%ymm4		,%%ymm6		/* ymm6 <- cpy t5 */			\n\t		vmovaps		%%ymm12		,%%ymm14		/* ymm14 <- cpy t5 */\n\t"\
		"vsubpd		     (%%rcx),%%ymm4,%%ymm4		/* ~t7 <- t5 -rt */		\n\t		vsubpd		0x100(%%rcx),%%ymm12,%%ymm12		/* ~t7 <- t5 -rt */\n\t"\
		"vsubpd		0x020(%%rcx),%%ymm5,%%ymm5		/* ~t8 <- t6 -it */		\n\t		vsubpd		0x120(%%rcx),%%ymm13,%%ymm13		/* ~t8 <- t6 -it */\n\t"\
		"vaddpd		     (%%rcx),%%ymm6,%%ymm6		/* ~t5 <- t5 +rt */		\n\t		vaddpd		0x100(%%rcx),%%ymm14,%%ymm14		/* ~t5 <- t5 +rt */\n\t"\
		"vaddpd		0x020(%%rcx),%%ymm7,%%ymm7		/* ~t6 <- t6 +it */		\n\t		vaddpd		0x120(%%rcx),%%ymm15,%%ymm15		/* ~t6 <- t6 +it */\n\t"\
	"prefetcht1	0x140(%%r14)\n\t"\
		"\n\t"\
	"/* Finish radix-4 butterfly and store results into temp-array slots: */\n\t"\
		"vsubpd		%%ymm6,%%ymm0,%%ymm0	/*~t5 */						\n\t		vsubpd		%%ymm14,%%ymm8 ,%%ymm8 	/*~t5 */						\n\t"\
		"vsubpd		%%ymm7,%%ymm1,%%ymm1	/*~t6 */						\n\t		vsubpd		%%ymm15,%%ymm9 ,%%ymm9 	/*~t6 */						\n\t"\
		"vsubpd		%%ymm5,%%ymm2,%%ymm2	/*~t3 */						\n\t		vsubpd		%%ymm13,%%ymm10,%%ymm10	/*~t3 */\n\t"\
		"vsubpd		%%ymm4,%%ymm3,%%ymm3	/*~t8 */						\n\t		vsubpd		%%ymm12,%%ymm11,%%ymm11	/*~t8 */\n\t"\
		"vmovaps	%%ymm0,0x080(%%rcx)	/* a[jt+p8 ] <- ~t5 */				\n\t		vmovaps		%%ymm8 ,0x180(%%rcx)	/* a[jt+p8 ] <- ~t5 */	\n\t"\
		"vmovaps	%%ymm1,0x0a0(%%rcx)	/* a[jp+p8 ] <- ~t6 */				\n\t		vmovaps		%%ymm9 ,0x1a0(%%rcx)	/* a[jp+p8 ] <- ~t6 */	\n\t"\
		"vmovaps	%%ymm2,0x040(%%rcx)	/* a[jt+p4 ] <- ~t3 */				\n\t		vmovaps		%%ymm10,0x140(%%rcx)	/* a[jt+p4 ] <- ~t3 */\n\t"\
		"vmovaps	%%ymm3,0x0e0(%%rcx)	/* a[jp+p12] <- ~t8 */				\n\t		vmovaps		%%ymm11,0x1e0(%%rcx)	/* a[jp+p12] <- ~t8 */\n\t"\
		"vaddpd		%%ymm6,%%ymm6,%%ymm6	/* 2*t5 */						\n\t		vaddpd		%%ymm14,%%ymm14,%%ymm14	/* 2*t5 */						\n\t"\
		"vaddpd		%%ymm7,%%ymm7,%%ymm7	/* 2*t6 */						\n\t		vaddpd		%%ymm15,%%ymm15,%%ymm15	/* 2*t6 */						\n\t"\
		"vaddpd		%%ymm5,%%ymm5,%%ymm5	/* 2*t8 */						\n\t		vaddpd		%%ymm13,%%ymm13,%%ymm13	/* 2*t8 */\n\t"\
		"vaddpd		%%ymm4,%%ymm4,%%ymm4	/* 2*t7 */						\n\t		vaddpd		%%ymm12,%%ymm12,%%ymm12	/* 2*t7 */\n\t"\
		"vaddpd		%%ymm0,%%ymm6,%%ymm6	/*~t1 */						\n\t		vaddpd		%%ymm8 ,%%ymm14,%%ymm14	/*~t1 */						\n\t"\
		"vaddpd		%%ymm1,%%ymm7,%%ymm7	/*~t2 */						\n\t		vaddpd		%%ymm9 ,%%ymm15,%%ymm15	/*~t2 */						\n\t"\
		"vaddpd		%%ymm2,%%ymm5,%%ymm5	/*~t7 */						\n\t		vaddpd		%%ymm10,%%ymm13,%%ymm13	/*~t7 */\n\t"\
		"vaddpd		%%ymm3,%%ymm4,%%ymm4	/*~t4 */						\n\t		vaddpd		%%ymm11,%%ymm12,%%ymm12	/*~t4 */\n\t"\
		"vmovaps	%%ymm6,     (%%rcx)/* a[jt    ] <- ~t1 */				\n\t		vmovaps		%%ymm14,0x100(%%rcx)	/* a[jt    ] <- ~t1 */	\n\t"\
		"vmovaps	%%ymm7,0x020(%%rcx)	/* a[jp    ] <- ~t2 */				\n\t		vmovaps		%%ymm15,0x120(%%rcx)	/* a[jp    ] <- ~t2 */	\n\t"\
		"vmovaps	%%ymm5,0x0c0(%%rcx)	/* a[jt+p12] <- ~t7 */				\n\t		vmovaps		%%ymm13,0x1c0(%%rcx)	/* a[jt+p12] <- ~t7 */\n\t"\
		"vmovaps	%%ymm4,0x060(%%rcx)	/* a[jp+p4 ] <- ~t4 */				\n\t		vmovaps		%%ymm12,0x160(%%rcx)	/* a[jp+p4 ] <- ~t4 */\n\t"\
		"\n\t"\
"/**************************************************************************************/\n\t"\
"/*...and now do four more radix-4 transforms, including the internal twiddle factors: */\n\t"\
"/**************************************************************************************/\n\t"\
		"\n\t"\
"/*...Block 1: t1,9,17,25 */\n\t											/*...Block 3: t5,13,21,29: All rax-offsets incr +0x40 in rcol w.r.to lcol: */\n\t"\
		"movq	%[__r1],%%rax							\n\t				vmovaps		0x080(%%rax),%%ymm8 		/* t5  */\n\t"\
		"movq	%[__r9],%%rbx							\n\t				vmovaps		0x0a0(%%rax),%%ymm9 		/* t6  */\n\t"\
		"movq	%[__isrt2],%%rsi						\n\t				vmovaps		0x1a0(%%rax),%%ymm11		/* t14 */\n\t"\
		"vmovaps		 (%%rax),%%ymm0		/* t1  */\n\t					vmovaps		0x180(%%rax),%%ymm10		/* t13 */\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1		/* t2  */\n\t					vsubpd		%%ymm11,%%ymm8 ,%%ymm8 		/* t5 =t5 -t14 */\n\t"\
		"vmovaps	0x100(%%rax),%%ymm2		/* t9  */\n\t					vsubpd		%%ymm10,%%ymm9 ,%%ymm9 		/* t14=t6 -t13 */\n\t"\
		"vmovaps	0x120(%%rax),%%ymm3		/* t14 */\n\t					vaddpd		%%ymm11,%%ymm11,%%ymm11		/* 2.t14 */\n\t"\
		"\n\t																vaddpd		%%ymm10,%%ymm10,%%ymm10		/* 2.t13 */\n\t"\
		"vsubpd		%%ymm2,%%ymm0,%%ymm0		/* t9 =t1 -t9  */\n\t		vaddpd		%%ymm8 ,%%ymm11,%%ymm11		/* t13=t5 +t14 */\n\t"\
		"vsubpd		%%ymm3,%%ymm1,%%ymm1		/* t14=t2 -t14 */\n\t		vaddpd		%%ymm9 ,%%ymm10,%%ymm10		/* t6 =t6 +t13 */\n\t"\
		"vaddpd		%%ymm2,%%ymm2,%%ymm2		/* 2.t9  */\n\t				vmovaps		0x280(%%rax),%%ymm12		/* t21 */\n\t"\
		"vaddpd		%%ymm3,%%ymm3,%%ymm3		/* 2.t14 */\n\t				vmovaps		0x2a0(%%rax),%%ymm13		/* t22 */\n\t"\
		"vaddpd		%%ymm0,%%ymm2,%%ymm2		/* t1 =t1 +t9  */\n\t		vmovaps		0x380(%%rax),%%ymm14		/* t29 */\n\t"\
		"vaddpd		%%ymm1,%%ymm3,%%ymm3		/* t2 =t2 +t14 */\n\t		vmovaps		0x3a0(%%rax),%%ymm15		/* t30 */\n\t"\
		"\n\t																vsubpd		0x2a0(%%rax),%%ymm12,%%ymm12		/* t21-t22 */\n\t"\
		"vmovaps	0x200(%%rax),%%ymm4		/* t17 */\n\t					vaddpd		0x280(%%rax),%%ymm13,%%ymm13		/* t22+t21 */\n\t"\
		"vmovaps	0x220(%%rax),%%ymm5		/* t18 */\n\t					vmulpd		(%%rsi)		,%%ymm12,%%ymm12	/* t21 = (t21-t22)*ISRT2 */\n\t"\
		"vmovaps	0x300(%%rax),%%ymm6		/* t25 */\n\t					vmulpd		(%%rsi)		,%%ymm13,%%ymm13	/* t22 = (t22+t21)*ISRT2 */\n\t"\
		"vmovaps	0x320(%%rax),%%ymm7		/* t26 */\n\t					vaddpd		0x3a0(%%rax),%%ymm14,%%ymm14		/* t29+t30 */\n\t"\
		"\n\t																vsubpd		0x380(%%rax),%%ymm15,%%ymm15		/* t30-t29 */\n\t"\
	"prefetcht1	0x180(%%r14)\n\t"\
		"vsubpd		%%ymm6,%%ymm4,%%ymm4		/* t25=t17-t25 */\n\t		vmulpd		(%%rsi),%%ymm14,%%ymm14	/*  rt = (t29+t30)*ISRT2 */\n\t"\
		"vsubpd		%%ymm7,%%ymm5,%%ymm5		/* t26=t18-t26 */\n\t		vmulpd		(%%rsi),%%ymm15,%%ymm15	/*  it = (t30-t29)*ISRT2 */\n\t"\
		"vaddpd		%%ymm6,%%ymm6,%%ymm6		/* 2.t25 */\n\t				vsubpd		%%ymm14,%%ymm12,%%ymm12		/* t21=t21-rt */\n\t"\
		"vaddpd		%%ymm7,%%ymm7,%%ymm7		/* 2.t26 */\n\t				vsubpd		%%ymm15,%%ymm13,%%ymm13		/* t22=t22-it */\n\t"\
		"vaddpd		%%ymm4,%%ymm6,%%ymm6		/* t17=t17+t25 */\n\t		vaddpd		%%ymm14,%%ymm14,%%ymm14		/*      2* rt */\n\t"\
		"vaddpd		%%ymm5,%%ymm7,%%ymm7		/* t18=t18+t26 */\n\t		vaddpd		%%ymm15,%%ymm15,%%ymm15		/*      2* it */\n\t"\
		"\n\t																vaddpd		%%ymm12,%%ymm14,%%ymm14		/* t29=t21+rt */\n\t"\
		"\n\t																vaddpd		%%ymm13,%%ymm15,%%ymm15		/* t30=t22+it */\n\t"\
		"movq	%[__r17],%%rcx						\n\t		movq		%[__r25],%%rdx			\n\t"\
		"vsubpd		%%ymm6,%%ymm2,%%ymm2						\n\t		vsubpd		%%ymm12,%%ymm8 ,%%ymm8 			\n\t"\
		"vsubpd		%%ymm7,%%ymm3,%%ymm3						\n\t		vsubpd		%%ymm13,%%ymm10,%%ymm10			\n\t"\
		"vaddpd		%%ymm6,%%ymm6,%%ymm6						\n\t		vaddpd		%%ymm12,%%ymm12,%%ymm12			\n\t"\
		"vaddpd		%%ymm7,%%ymm7,%%ymm7						\n\t		vaddpd		%%ymm13,%%ymm13,%%ymm13			\n\t"\
		"vmovaps		%%ymm2,     (%%rcx)						\n\t		vmovaps		%%ymm8 ,0x080(%%rcx)		\n\t"\
		"vmovaps		%%ymm3,0x020(%%rcx)						\n\t		vmovaps		%%ymm10,0x0a0(%%rcx)		\n\t"\
		"vaddpd		%%ymm2,%%ymm6,%%ymm6						\n\t		vaddpd		%%ymm8 ,%%ymm12,%%ymm12			\n\t"\
		"vaddpd		%%ymm3,%%ymm7,%%ymm7						\n\t		vaddpd		%%ymm10,%%ymm13,%%ymm13			\n\t"\
		"vmovaps		%%ymm6,     (%%rax)						\n\t		vmovaps		%%ymm12,0x080(%%rax)		\n\t"\
		"vmovaps		%%ymm7,0x020(%%rax)						\n\t		vmovaps		%%ymm13,0x0a0(%%rax)		\n\t"\
		"vsubpd		%%ymm5,%%ymm0,%%ymm0						\n\t		vsubpd		%%ymm15,%%ymm11,%%ymm11			\n\t"\
		"vsubpd		%%ymm4,%%ymm1,%%ymm1						\n\t		vsubpd		%%ymm14,%%ymm9 ,%%ymm9 			\n\t"\
		"vaddpd		%%ymm5,%%ymm5,%%ymm5						\n\t		vaddpd		%%ymm15,%%ymm15,%%ymm15			\n\t"\
		"vaddpd		%%ymm4,%%ymm4,%%ymm4						\n\t		vaddpd		%%ymm14,%%ymm14,%%ymm14			\n\t"\
		"vmovaps		%%ymm0,     (%%rbx)						\n\t		vmovaps		%%ymm11,0x080(%%rbx)		\n\t"\
		"vmovaps		%%ymm1,0x020(%%rdx)						\n\t		vmovaps		%%ymm9 ,0x0a0(%%rdx)		\n\t"\
		"vaddpd		%%ymm0,%%ymm5,%%ymm5						\n\t		vaddpd		%%ymm11,%%ymm15,%%ymm15			\n\t"\
		"vaddpd		%%ymm1,%%ymm4,%%ymm4						\n\t		vaddpd		%%ymm9 ,%%ymm14,%%ymm14			\n\t"\
		"vmovaps		%%ymm5,     (%%rdx)						\n\t		vmovaps		%%ymm15,0x080(%%rdx)		\n\t"\
		"vmovaps		%%ymm4,0x020(%%rbx)						\n\t		vmovaps		%%ymm14,0x0a0(%%rbx)		\n\t"\
		"\n\t"\
"/*...Block 2: t3,11,19,27 */\n\t"\
		"addq	$0x040,%%rax	/* r3  */\n\t"\
		"addq	$0x040,%%rbx							\n\t"\
		"addq	$0x040,%%rcx							\n\t"\
		"addq	$0x040,%%rdx							\n\t"\
		"leaq	0x020(%%rsi),%%rdi	/* cc0, from isrt2 */\n\t"\
		"\n\t																/*...Block 4: t7,15,23,31 */\n\t"\
		"vmovaps	0x200(%%rax),%%ymm4		/* t19 */			\n\t		vmovaps		0x280(%%rax),%%ymm12		/* t23 */\n\t"\
		"vmovaps	0x220(%%rax),%%ymm5		/* t20 */			\n\t		vmovaps		0x2a0(%%rax),%%ymm13		/* t24 */\n\t"\
		"vmovaps	0x300(%%rax),%%ymm6	/* t27 */\n\t						vmovaps		0x380(%%rax),%%ymm14		/* t31 */\n\t"\
		"vmovaps	0x320(%%rax),%%ymm7	/* t28 */\n\t						vmovaps		0x3a0(%%rax),%%ymm15		/* t32 */\n\t"\
		"vmovaps	0x200(%%rax),%%ymm0		/* copy t19 */		\n\t		vmovaps		0x280(%%rax),%%ymm8 		/* copy t23 */\n\t"\
		"vmovaps	0x220(%%rax),%%ymm1		/* copy t20 */		\n\t		vmovaps		0x2a0(%%rax),%%ymm9 		/* copy t24 */\n\t"\
		"vmovaps	0x300(%%rax),%%ymm2	/* copy t27 */\n\t					vmovaps		0x380(%%rax),%%ymm10		/* copy t31 */\n\t"\
		"vmovaps	0x320(%%rax),%%ymm3	/* copy t28 */\n\t					vmovaps		0x3a0(%%rax),%%ymm11		/* copy t32 */\n\t"\
		"\n\t"\
		"vmulpd		     (%%rdi),%%ymm4,%%ymm4		/* t19*c */	\n\t		vmulpd		0x020(%%rdi),%%ymm12,%%ymm12		/* t23*s */\n\t"\
		"vmulpd		0x020(%%rdi),%%ymm1,%%ymm1		/* t20*s */	\n\t		vmulpd		     (%%rdi),%%ymm9 ,%%ymm9 		/* t24*c */\n\t"\
		"vmulpd		0x020(%%rdi),%%ymm6,%%ymm6		/* t27*s */	\n\t		vmulpd		     (%%rdi),%%ymm14,%%ymm14		/* t31*c */\n\t"\
		"vmulpd		     (%%rdi),%%ymm3,%%ymm3		/* t28*c */	\n\t		vmulpd		0x020(%%rdi),%%ymm11,%%ymm11		/* t32*s */\n\t"\
		"vmulpd		     (%%rdi),%%ymm5,%%ymm5		/* t20*c */	\n\t		vmulpd		0x020(%%rdi),%%ymm13,%%ymm13		/* t24*s */\n\t"\
		"vmulpd		0x020(%%rdi),%%ymm7,%%ymm7		/* t28*s */	\n\t		vmulpd		     (%%rdi),%%ymm8 ,%%ymm8 		/* t23*c */\n\t"\
		"vmulpd		0x020(%%rdi),%%ymm0,%%ymm0		/* t19*s */	\n\t		vmulpd		     (%%rdi),%%ymm15,%%ymm15		/* t32*c */\n\t"\
		"vmulpd		     (%%rdi),%%ymm2,%%ymm2		/* t27*c */	\n\t		vmulpd		0x020(%%rdi),%%ymm10,%%ymm10		/* t31*s */\n\t"\
		"vsubpd		%%ymm1,%%ymm4,%%ymm4	/* ~t19 */			\n\t		vsubpd		%%ymm9 ,%%ymm12,%%ymm12		/* ~t23 */\n\t"\
		"vsubpd		%%ymm3,%%ymm6,%%ymm6	/* rt */\n\t					vaddpd		%%ymm8 ,%%ymm13,%%ymm13		/* ~t24 */\n\t"\
		"vaddpd		%%ymm0,%%ymm5,%%ymm5	/* ~t20 */			\n\t		vsubpd		%%ymm11,%%ymm14,%%ymm14		/* rt */\n\t"\
		"vaddpd		%%ymm2,%%ymm7,%%ymm7	/* it */\n\t					vaddpd		%%ymm10,%%ymm15,%%ymm15		/* it */\n\t"\
		"\n\t"\
		"vsubpd		%%ymm6,%%ymm4,%%ymm4		/*~t27=t19-rt */\n\t		vsubpd		%%ymm14,%%ymm12,%%ymm12		/*~t23=t23-rt */\n\t"\
		"vsubpd		%%ymm7,%%ymm5,%%ymm5		/*~t28=t20-it */\n\t		vsubpd		%%ymm15,%%ymm13,%%ymm13		/*~t24=t24-it */\n\t"\
		"vaddpd		%%ymm6,%%ymm6,%%ymm6		/*      2* rt */\n\t		vaddpd		%%ymm14,%%ymm14,%%ymm14		/*      2* rt */\n\t"\
		"vaddpd		%%ymm7,%%ymm7,%%ymm7		/*      2* it */\n\t		vaddpd		%%ymm15,%%ymm15,%%ymm15		/*      2* it */\n\t"\
		"vaddpd		%%ymm4,%%ymm6,%%ymm6		/*~t19=t19+rt */\n\t		vaddpd		%%ymm12,%%ymm14,%%ymm14		/*~t31=t23+rt */\n\t"\
		"vaddpd		%%ymm5,%%ymm7,%%ymm7		/*~t20=t20+it */\n\t		vaddpd		%%ymm13,%%ymm15,%%ymm15		/*~t32=t24+it */\n\t"\
	"prefetcht1	0x1c0(%%r14)\n\t"\
		"\n\t"\
		"vmovaps	0x100(%%rax),%%ymm2		/* t11 */\n\t					vmovaps		0x180(%%rax),%%ymm10		/* t15 */\n\t"\
		"vmovaps	0x120(%%rax),%%ymm3		/* t12 */\n\t					vmovaps		0x1a0(%%rax),%%ymm11		/* t16 */\n\t"\
		"vsubpd		0x120(%%rax),%%ymm2,%%ymm2		/* t11-t12 */\n\t		vaddpd		0x1a0(%%rax),%%ymm10,%%ymm10	/* t15+t16 */\n\t"\
		"vaddpd		0x100(%%rax),%%ymm3,%%ymm3		/* t12+t11 */\n\t		vsubpd		0x180(%%rax),%%ymm11,%%ymm11	/* t16-t15 */\n\t"\
		"vmulpd		(%%rsi),%%ymm2,%%ymm2	/* rt = (t11-t12)*ISRT2 */\n\t	vmulpd		(%%rsi)		,%%ymm10,%%ymm10	/* rt = (t15+t16)*ISRT2 */\n\t"\
		"vmulpd		(%%rsi),%%ymm3,%%ymm3	/* it = (t12+t11)*ISRT2 */\n\t	vmulpd		(%%rsi)		,%%ymm11,%%ymm11	/* it = (t16-t15)*ISRT2 */\n\t"\
		"\n\t"\
		"vmovaps	     (%%rax),%%ymm0		/* t3  */\n\t					vmovaps		0x080(%%rax),%%ymm8 		/* t7  */\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1		/* t4  */\n\t					vmovaps		0x0a0(%%rax),%%ymm9 		/* t8  */\n\t"\
		"\n\t"\
		"vsubpd		      %%ymm2,%%ymm0,%%ymm0		/*~t11=t3 -rt */\n\t	vsubpd		     %%ymm10,%%ymm8 ,%%ymm8 	/*~t7 =t7 -rt */\n\t"\
		"vsubpd		      %%ymm3,%%ymm1,%%ymm1		/*~t12=t4 -it */\n\t	vsubpd		     %%ymm11,%%ymm9 ,%%ymm9 	/*~t8 =t8 -it */\n\t"\
		"vaddpd		     (%%rax),%%ymm2,%%ymm2		/*~t3 =rt +t3 */\n\t	vaddpd		0x080(%%rax),%%ymm10,%%ymm10	/*~t15=rt +t7 */\n\t"\
		"vaddpd		0x020(%%rax),%%ymm3,%%ymm3		/*~t4 =it +t4 */\n\t	vaddpd		0x0a0(%%rax),%%ymm11,%%ymm11	/*~t16=it +t8 */\n\t"\
		"\n\t"\
		"vsubpd		%%ymm6,%%ymm2,%%ymm2		/* t3 -t19 */\n\t			vsubpd		%%ymm12		,%%ymm8 ,%%ymm8 	/* t7 -t23 */\n\t"\
		"vsubpd		%%ymm7,%%ymm3,%%ymm3		/* t4 -t20 */\n\t			vsubpd		%%ymm13		,%%ymm9 ,%%ymm9 	/* t8 -t24 */\n\t"\
		"vaddpd		%%ymm6,%%ymm6,%%ymm6		/*   2*t19 */\n\t			vaddpd		%%ymm12		,%%ymm12,%%ymm12	/*   2*t23 */\n\t"\
		"vaddpd		%%ymm7,%%ymm7,%%ymm7		/*   2*t20 */\n\t			vaddpd		%%ymm13		,%%ymm13,%%ymm13	/*   2*t24 */\n\t"\
		"vmovaps		%%ymm2,     (%%rcx)						\n\t		vmovaps		%%ymm8 ,0x080(%%rcx)		\n\t"\
		"vmovaps		%%ymm3,0x020(%%rcx)						\n\t		vmovaps		%%ymm9 ,0x0a0(%%rcx)		\n\t"\
		"vaddpd		%%ymm2,%%ymm6,%%ymm6						\n\t		vaddpd		%%ymm8 ,%%ymm12,%%ymm12			\n\t"\
		"vaddpd		%%ymm3,%%ymm7,%%ymm7						\n\t		vaddpd		%%ymm9 ,%%ymm13,%%ymm13			\n\t"\
		"vmovaps		%%ymm6,     (%%rax)						\n\t		vmovaps		%%ymm12,0x080(%%rax)		\n\t"\
		"vmovaps		%%ymm7,0x020(%%rax)						\n\t		vmovaps		%%ymm13,0x0a0(%%rax)		\n\t"\
		"vsubpd		%%ymm5,%%ymm0,%%ymm0						\n\t		vsubpd		%%ymm15,%%ymm10,%%ymm10			\n\t"\
		"vsubpd		%%ymm4,%%ymm1,%%ymm1						\n\t		vsubpd		%%ymm14,%%ymm11,%%ymm11			\n\t"\
		"vaddpd		%%ymm5,%%ymm5,%%ymm5						\n\t		vaddpd		%%ymm15,%%ymm15,%%ymm15			\n\t"\
		"vaddpd		%%ymm4,%%ymm4,%%ymm4						\n\t		vaddpd		%%ymm14,%%ymm14,%%ymm14			\n\t"\
		"vmovaps		%%ymm0,     (%%rbx)						\n\t		vmovaps		%%ymm10,0x080(%%rbx)		\n\t"\
		"vmovaps		%%ymm1,0x020(%%rdx)						\n\t		vmovaps		%%ymm11,0x0a0(%%rdx)		\n\t"\
		"vaddpd		%%ymm0,%%ymm5,%%ymm5						\n\t		vaddpd		%%ymm10,%%ymm15,%%ymm15			\n\t"\
		"vaddpd		%%ymm1,%%ymm4,%%ymm4						\n\t		vaddpd		%%ymm11,%%ymm14,%%ymm14			\n\t"\
		"vmovaps		%%ymm5,     (%%rdx)						\n\t		vmovaps		%%ymm15,0x080(%%rdx)		\n\t"\
		"vmovaps		%%ymm4,0x020(%%rbx)						\n\t		vmovaps		%%ymm14,0x0a0(%%rbx)		\n\t"\
		 :					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__add1] "m" (Xadd1)\
		 ,[__add2] "m" (Xadd2)\
		 ,[__add3] "m" (Xadd3)\
		 ,[__r1] "m" (Xr1)\
		 ,[__r9] "m" (Xr9)\
		 ,[__r17] "m" (Xr17)\
		 ,[__r25] "m" (Xr25)\
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__cc0] "m" (Xcc0)\
		 ,[__c1] "m" (Xc1)\
		 ,[__c2] "m" (Xc2)\
		 ,[__c3] "m" (Xc3)\
		 ,[__c4] "m" (Xc4)\
		 ,[__c5] "m" (Xc5)\
		 ,[__c6] "m" (Xc6)\
		 ,[__c7] "m" (Xc7)\
		 ,[__c8] "m" (Xc8)\
		 ,[__c9] "m" (Xc9)\
		 ,[__c10] "m" (Xc10)\
		 ,[__c11] "m" (Xc11)\
		 ,[__c12] "m" (Xc12)\
		 ,[__c13] "m" (Xc13)\
		 ,[__c14] "m" (Xc14)\
		 ,[__c15] "m" (Xc15)\
	 ,[__pfetch_dist] "m" (Xpfetch_dist)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r14","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"		/* Clobbered registers */\
	);\
	}

	// Since the add0/1 block addresses advance in opposite directions, our prefetch scheme here is like so:
	// In fwd DIF prefetch  ahead 8 cache lines' worth of data w.r.to add0;
	// In inv DIT prefetch behind 8 cache lines' worth of data w.r.to add1;
	/****** To-Do: Need a 16-register version! ****************/
	#define SSE2_RADIX16_WRAPPER_DIT(Xadd0,Xadd1,Xadd2,Xadd3,Xr1,Xr9,Xr17,Xr25,Xisrt2,Xcc0,Xc1,Xc2,Xc3,Xc4,Xc5,Xc6,Xc7,Xc8,Xc9,Xc10,Xc11,Xc12,Xc13,Xc14,Xc15,Xpfetch_dist)\
	{\
	__asm__ volatile (\
	/*...Block 1: r1,9,17,25 */\
		"movq	%[__r1],%%rax					\n\t"\
	"movq	%[__add1],%%r14	\n\t"\
	"movslq	%[__pfetch_dist],%%r14	\n\t"\
	"leaq	(%%rax,%%r14,8),%%r14	\n\t"	/* Block 1 [base-address + data-fetch-ahead index] */\
	"prefetcht1	-0x40(%%r14)\n\t"\
		"movq	%%rax,%%rbx						\n\t"\
		"movq	%%rax,%%rcx						\n\t"	/*...Block 3: r3,11,19,27 */\
		"movq	%%rax,%%rdx						\n\t		leaq	0x40(%%rax),%%r10				\n\t"\
		"addq	$0x200,%%rbx					\n\t		leaq	0x40(%%rbx),%%r11				\n\t"\
		"addq	$0x100,%%rcx					\n\t		leaq	0x40(%%rcx),%%r12				\n\t"\
		"addq	$0x300,%%rdx					\n\t		leaq	0x40(%%rdx),%%r13				\n\t"\
		"/* SSE2_RADIX4_DIT_IN_PLACE(): */		\n\t		/* SSE2_RADIX4_DIT_IN_PLACE(): */		\n\t"\
		"vmovaps	     (%%rax),%%ymm0			\n\t		vmovaps	     (%%r10),%%ymm8 			\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1			\n\t		vmovaps	0x020(%%r10),%%ymm9 			\n\t"\
		"vmovaps	     (%%rax),%%ymm2			\n\t		vmovaps	     (%%r10),%%ymm10			\n\t"\
		"vmovaps	0x020(%%rax),%%ymm3			\n\t		vmovaps	0x020(%%r10),%%ymm11			\n\t"\
		"vaddpd	     (%%rbx),%%ymm0,%%ymm0		\n\t		vaddpd	     (%%r11),%%ymm8 ,%%ymm8 	\n\t"\
		"vaddpd	0x020(%%rbx),%%ymm1,%%ymm1		\n\t		vaddpd	0x020(%%r11),%%ymm9 ,%%ymm9 	\n\t"\
		"vsubpd	     (%%rbx),%%ymm2,%%ymm2		\n\t		vsubpd	     (%%r11),%%ymm10,%%ymm10	\n\t"\
		"vsubpd	0x020(%%rbx),%%ymm3,%%ymm3		\n\t		vsubpd	0x020(%%r11),%%ymm11,%%ymm11	\n\t"\
		"vmovaps	     (%%rcx),%%ymm4			\n\t		vmovaps	     (%%r12),%%ymm12			\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5			\n\t		vmovaps	0x020(%%r12),%%ymm13			\n\t"\
		"vmovaps	     (%%rcx),%%ymm6			\n\t		vmovaps	     (%%r12),%%ymm14			\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm7			\n\t		vmovaps	0x020(%%r12),%%ymm15			\n\t"\
	"prefetcht1	-0xc0(%%r14)\n\t"\
		"vaddpd	     (%%rdx),%%ymm4,%%ymm4		\n\t		vaddpd	     (%%r13),%%ymm12,%%ymm12	\n\t"\
		"vaddpd	0x020(%%rdx),%%ymm5,%%ymm5		\n\t		vaddpd	0x020(%%r13),%%ymm13,%%ymm13	\n\t"\
		"vsubpd	     (%%rdx),%%ymm6,%%ymm6		\n\t		vsubpd	     (%%r13),%%ymm14,%%ymm14	\n\t"\
		"vsubpd	0x020(%%rdx),%%ymm7,%%ymm7		\n\t		vsubpd	0x020(%%r13),%%ymm15,%%ymm15	\n\t"\
		"vsubpd	      %%ymm4,%%ymm0,%%ymm0		\n\t		vsubpd	      %%ymm12,%%ymm8 ,%%ymm8 	\n\t"\
		"vsubpd	      %%ymm5,%%ymm1,%%ymm1		\n\t		vsubpd	      %%ymm13,%%ymm9 ,%%ymm9 	\n\t"\
		"vmovaps	%%ymm0,     (%%rbx)			\n\t		vmovaps	%%ymm8 ,     (%%r11)			\n\t"\
		"vmovaps	%%ymm1,0x020(%%rbx)			\n\t		vmovaps	%%ymm9 ,0x020(%%r11)			\n\t"\
		"vaddpd	      %%ymm4,%%ymm4,%%ymm4		\n\t		vaddpd	      %%ymm12,%%ymm12,%%ymm12	\n\t"\
		"vaddpd	      %%ymm5,%%ymm5,%%ymm5		\n\t		vaddpd	      %%ymm13,%%ymm13,%%ymm13	\n\t"\
		"vaddpd	      %%ymm0,%%ymm4,%%ymm4		\n\t		vaddpd	      %%ymm8 ,%%ymm12,%%ymm12	\n\t"\
		"vaddpd	      %%ymm1,%%ymm5,%%ymm5		\n\t		vaddpd	      %%ymm9 ,%%ymm13,%%ymm13	\n\t"\
		"vmovaps	%%ymm4,     (%%rax)			\n\t		vmovaps	%%ymm12,     (%%r10)			\n\t"\
		"vmovaps	%%ymm5,0x020(%%rax)			\n\t		vmovaps	%%ymm13,0x020(%%r10)			\n\t"\
		"vsubpd	      %%ymm7,%%ymm2,%%ymm2		\n\t		vsubpd	      %%ymm15,%%ymm10,%%ymm10	\n\t"\
		"vsubpd	      %%ymm6,%%ymm3,%%ymm3		\n\t		vsubpd	      %%ymm14,%%ymm11,%%ymm11	\n\t"\
		"vmovaps	%%ymm2,     (%%rdx)			\n\t		vmovaps	%%ymm10,     (%%r13)			\n\t"\
		"vmovaps	%%ymm3,0x020(%%rcx)			\n\t		vmovaps	%%ymm11,0x020(%%r12)			\n\t"\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7		\n\t		vaddpd	      %%ymm15,%%ymm15,%%ymm15	\n\t"\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6		\n\t		vaddpd	      %%ymm14,%%ymm14,%%ymm14	\n\t"\
		"vaddpd	      %%ymm2,%%ymm7,%%ymm7		\n\t		vaddpd	      %%ymm10,%%ymm15,%%ymm15	\n\t"\
		"vaddpd	      %%ymm3,%%ymm6,%%ymm6		\n\t		vaddpd	      %%ymm11,%%ymm14,%%ymm14	\n\t"\
		"vmovaps	%%ymm7,     (%%rcx)			\n\t		vmovaps	%%ymm15,     (%%r12)			\n\t"\
		"vmovaps	%%ymm6,0x020(%%rdx)			\n\t		vmovaps	%%ymm14,0x020(%%r13)			\n\t"\
	/*...Block 2: r5,13,21,29 */						/*...Block 4: r7,15,23,31 */\
		"addq	$0x80,%%rax						\n\t		leaq	0x40(%%rax),%%r10				\n\t"\
		"addq	$0x80,%%rbx						\n\t		leaq	0x40(%%rbx),%%r11				\n\t"\
		"addq	$0x80,%%rcx						\n\t		leaq	0x40(%%rcx),%%r12				\n\t"\
		"addq	$0x80,%%rdx						\n\t		leaq	0x40(%%rdx),%%r13				\n\t"\
		"/* SSE2_RADIX4_DIT_IN_PLACE(): */		\n\t		/* SSE2_RADIX4_DIT_IN_PLACE(): */		\n\t"\
		"vmovaps	     (%%rax),%%ymm0			\n\t		vmovaps	     (%%r10),%%ymm8 			\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1			\n\t		vmovaps	0x020(%%r10),%%ymm9 			\n\t"\
		"vmovaps	     (%%rax),%%ymm2			\n\t		vmovaps	     (%%r10),%%ymm10			\n\t"\
		"vmovaps	0x020(%%rax),%%ymm3			\n\t		vmovaps	0x020(%%r10),%%ymm11			\n\t"\
		"vaddpd	     (%%rbx),%%ymm0,%%ymm0		\n\t		vaddpd	     (%%r11),%%ymm8 ,%%ymm8 	\n\t"\
		"vaddpd	0x020(%%rbx),%%ymm1,%%ymm1		\n\t		vaddpd	0x020(%%r11),%%ymm9 ,%%ymm9 	\n\t"\
		"vsubpd	     (%%rbx),%%ymm2,%%ymm2		\n\t		vsubpd	     (%%r11),%%ymm10,%%ymm10	\n\t"\
		"vsubpd	0x020(%%rbx),%%ymm3,%%ymm3		\n\t		vsubpd	0x020(%%r11),%%ymm11,%%ymm11	\n\t"\
		"vmovaps	     (%%rcx),%%ymm4			\n\t		vmovaps	     (%%r12),%%ymm12			\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5			\n\t		vmovaps	0x020(%%r12),%%ymm13			\n\t"\
		"vmovaps	     (%%rcx),%%ymm6			\n\t		vmovaps	     (%%r12),%%ymm14			\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm7			\n\t		vmovaps	0x020(%%r12),%%ymm15			\n\t"\
		"vaddpd	     (%%rdx),%%ymm4,%%ymm4		\n\t		vaddpd	     (%%r13),%%ymm12,%%ymm12	\n\t"\
		"vaddpd	0x020(%%rdx),%%ymm5,%%ymm5		\n\t		vaddpd	0x020(%%r13),%%ymm13,%%ymm13	\n\t"\
		"vsubpd	     (%%rdx),%%ymm6,%%ymm6		\n\t		vsubpd	     (%%r13),%%ymm14,%%ymm14	\n\t"\
		"vsubpd	0x020(%%rdx),%%ymm7,%%ymm7		\n\t		vsubpd	0x020(%%r13),%%ymm15,%%ymm15	\n\t"\
		"vsubpd	      %%ymm4,%%ymm0,%%ymm0		\n\t		vsubpd	      %%ymm12,%%ymm8 ,%%ymm8 	\n\t"\
		"vsubpd	      %%ymm5,%%ymm1,%%ymm1		\n\t		vsubpd	      %%ymm13,%%ymm9 ,%%ymm9 	\n\t"\
	"prefetcht1	-0x80(%%r14)\n\t"\
		"vmovaps	%%ymm0,     (%%rbx)			\n\t		vmovaps	%%ymm8 ,     (%%r11)			\n\t"\
		"vmovaps	%%ymm1,0x020(%%rbx)			\n\t		vmovaps	%%ymm9 ,0x020(%%r11)			\n\t"\
		"vaddpd	      %%ymm4,%%ymm4,%%ymm4		\n\t		vaddpd	      %%ymm12,%%ymm12,%%ymm12	\n\t"\
		"vaddpd	      %%ymm5,%%ymm5,%%ymm5		\n\t		vaddpd	      %%ymm13,%%ymm13,%%ymm13	\n\t"\
		"vaddpd	      %%ymm0,%%ymm4,%%ymm4		\n\t		vaddpd	      %%ymm8 ,%%ymm12,%%ymm12	\n\t"\
		"vaddpd	      %%ymm1,%%ymm5,%%ymm5		\n\t		vaddpd	      %%ymm9 ,%%ymm13,%%ymm13	\n\t"\
		"vmovaps	%%ymm4,     (%%rax)			\n\t		vmovaps	%%ymm12,     (%%r10)			\n\t"\
		"vmovaps	%%ymm5,0x020(%%rax)			\n\t		vmovaps	%%ymm13,0x020(%%r10)			\n\t"\
		"vsubpd	      %%ymm7,%%ymm2,%%ymm2		\n\t		vsubpd	      %%ymm15,%%ymm10,%%ymm10	\n\t"\
		"vsubpd	      %%ymm6,%%ymm3,%%ymm3		\n\t		vsubpd	      %%ymm14,%%ymm11,%%ymm11	\n\t"\
		"vmovaps	%%ymm2,     (%%rdx)			\n\t		vmovaps	%%ymm10,     (%%r13)			\n\t"\
		"vmovaps	%%ymm3,0x020(%%rcx)			\n\t		vmovaps	%%ymm11,0x020(%%r12)			\n\t"\
		"vaddpd	      %%ymm7,%%ymm7,%%ymm7		\n\t		vaddpd	      %%ymm15,%%ymm15,%%ymm15	\n\t"\
		"vaddpd	      %%ymm6,%%ymm6,%%ymm6		\n\t		vaddpd	      %%ymm14,%%ymm14,%%ymm14	\n\t"\
		"vaddpd	      %%ymm2,%%ymm7,%%ymm7		\n\t		vaddpd	      %%ymm10,%%ymm15,%%ymm15	\n\t"\
		"vaddpd	      %%ymm3,%%ymm6,%%ymm6		\n\t		vaddpd	      %%ymm11,%%ymm14,%%ymm14	\n\t"\
		"vmovaps	%%ymm7,     (%%rcx)			\n\t		vmovaps	%%ymm15,     (%%r12)			\n\t"\
		"vmovaps	%%ymm6,0x020(%%rdx)			\n\t		vmovaps	%%ymm14,0x020(%%r13)			\n\t"\
		"\n\t"\
	/***************************************************************************************************/\
	/*..and now do four more radix-4 transforms, including the internal and external twiddle factors.  */\
	/*  Write even-index 16-byte output pairs to a(j1), odd-index to a(j2), unpack same as on inputs.  */\
	/*  We do the last 2 radix-4 blocks first, to make the unpack-interleaving of outputs smoother.    */\
	/***************************************************************************************************/\
	"/* Main-array addresses still in add0,1; all rax-offsets incr +0x200 in rcol w.r.to lcol: */\n\t"\
	"prefetcht1	-0x100(%%r14)\n\t"\
	/*...Block 3: t3,11,19,27 -> r9,13,11,15: */		/*...Block 4: t7,15,23,31 -> r25,29,27,31: */\
		"movq		%[__r9],%%rax				\n\t"\
		"movq		%[__isrt2],%%rbx			\n\t"\
		"movq		%[__cc0],%%rcx				\n\t"\
		"vmovaps	0x040(%%rax),%%ymm4			\n\t		vmovaps		0x240(%%rax),%%ymm12		\n\t"\
		"vmovaps	0x0c0(%%rax),%%ymm0			\n\t		vmovaps		0x2c0(%%rax),%%ymm8 		\n\t"\
		"vmovaps	0x060(%%rax),%%ymm5			\n\t		vmovaps		0x260(%%rax),%%ymm13		\n\t"\
		"vmovaps	0x0e0(%%rax),%%ymm1			\n\t		vmovaps		0x2e0(%%rax),%%ymm9 		\n\t"\
		"vmulpd		     (%%rcx),%%ymm4,%%ymm4	\n\t		vmulpd		0x020(%%rcx),%%ymm12,%%ymm12\n\t"\
		"vmulpd		0x020(%%rcx),%%ymm0,%%ymm0	\n\t		vmulpd		     (%%rcx),%%ymm8 ,%%ymm8 \n\t"\
		"vmulpd		     (%%rcx),%%ymm5,%%ymm5	\n\t		vmulpd		0x020(%%rcx),%%ymm13,%%ymm13\n\t"\
		"vmulpd		0x020(%%rcx),%%ymm1,%%ymm1	\n\t		vmulpd		     (%%rcx),%%ymm9 ,%%ymm9 \n\t"\
		"vmovaps	0x040(%%rax),%%ymm6			\n\t		vmovaps		0x240(%%rax),%%ymm14		\n\t"\
		"vmovaps	0x0c0(%%rax),%%ymm2			\n\t		vmovaps		0x2c0(%%rax),%%ymm10		\n\t"\
		"vmovaps	0x060(%%rax),%%ymm7			\n\t		vmovaps		0x260(%%rax),%%ymm15		\n\t"\
		"vmovaps	0x0e0(%%rax),%%ymm3			\n\t		vmovaps		0x2e0(%%rax),%%ymm11		\n\t"\
		"vmulpd		0x020(%%rcx),%%ymm6,%%ymm6	\n\t		vmulpd		     (%%rcx),%%ymm14,%%ymm14\n\t"\
		"vmulpd		     (%%rcx),%%ymm2,%%ymm2	\n\t		vmulpd		0x020(%%rcx),%%ymm10,%%ymm10\n\t"\
		"vmulpd		0x020(%%rcx),%%ymm7,%%ymm7	\n\t		vmulpd		     (%%rcx),%%ymm15,%%ymm15\n\t"\
		"vmulpd		     (%%rcx),%%ymm3,%%ymm3	\n\t		vmulpd		0x020(%%rcx),%%ymm11,%%ymm11\n\t"\
		"vsubpd		%%ymm6,%%ymm5,%%ymm5		\n\t		vsubpd		%%ymm14,%%ymm13,%%ymm13		\n\t"\
		"vsubpd		%%ymm2,%%ymm1,%%ymm1		\n\t		vsubpd		%%ymm10,%%ymm9 ,%%ymm9 		\n\t"\
		"vaddpd		%%ymm7,%%ymm4,%%ymm4		\n\t		vaddpd		%%ymm15,%%ymm12,%%ymm12		\n\t"\
		"vaddpd		%%ymm3,%%ymm0,%%ymm0		\n\t		vaddpd		%%ymm11,%%ymm8 ,%%ymm8 		\n\t"\
		"vmovaps	%%ymm5,%%ymm7				\n\t		vmovaps		%%ymm13,%%ymm15				\n\t"\
		"vmovaps	%%ymm4,%%ymm6				\n\t		vmovaps		%%ymm12,%%ymm14				\n\t"\
		"vaddpd		%%ymm0,%%ymm4,%%ymm4		\n\t		vaddpd		%%ymm8 ,%%ymm12,%%ymm12		\n\t"\
		"vaddpd		%%ymm1,%%ymm5,%%ymm5		\n\t		vaddpd		%%ymm9 ,%%ymm13,%%ymm13		\n\t"\
		"vsubpd		%%ymm0,%%ymm6,%%ymm6		\n\t		vsubpd		%%ymm8 ,%%ymm14,%%ymm14		\n\t"\
		"vsubpd		%%ymm1,%%ymm7,%%ymm7		\n\t		vsubpd		%%ymm9 ,%%ymm15,%%ymm15		\n\t"\
		"vmovaps	0x080(%%rax),%%ymm2			\n\t		vmovaps		0x280(%%rax),%%ymm10		\n\t"\
		"vmovaps	0x0a0(%%rax),%%ymm3			\n\t		vmovaps		0x2a0(%%rax),%%ymm11		\n\t"\
		"vmovaps	     (%%rax),%%ymm0			\n\t		vmovaps		0x200(%%rax),%%ymm8 		\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1			\n\t		vmovaps		0x220(%%rax),%%ymm9 		\n\t"\
		"vaddpd		0x0a0(%%rax),%%ymm2,%%ymm2	\n\t		vsubpd		0x2a0(%%rax),%%ymm10,%%ymm10\n\t"\
		"vsubpd		0x080(%%rax),%%ymm3,%%ymm3	\n\t		vaddpd		0x280(%%rax),%%ymm11,%%ymm11\n\t"\
		"vmulpd		     (%%rbx),%%ymm2,%%ymm2	\n\t		vmulpd		(%%rbx),%%ymm10,%%ymm10		\n\t"\
		"vmulpd		     (%%rbx),%%ymm3,%%ymm3	\n\t		vmulpd		(%%rbx),%%ymm11,%%ymm11		\n\t"\
		"vsubpd		%%ymm2,%%ymm0,%%ymm0		\n\t		vsubpd		%%ymm10,%%ymm8 ,%%ymm8 		\n\t"\
		"vsubpd		%%ymm3,%%ymm1,%%ymm1		\n\t		vsubpd		%%ymm11,%%ymm9 ,%%ymm9 		\n\t"\
		"vaddpd		%%ymm2,%%ymm2,%%ymm2		\n\t		vaddpd		%%ymm10,%%ymm10,%%ymm10		\n\t"\
		"vaddpd		%%ymm3,%%ymm3,%%ymm3		\n\t		vaddpd		%%ymm11,%%ymm11,%%ymm11		\n\t"\
		"vaddpd		%%ymm0,%%ymm2,%%ymm2		\n\t		vaddpd		%%ymm8 ,%%ymm10,%%ymm10		\n\t"\
		"vaddpd		%%ymm1,%%ymm3,%%ymm3		\n\t		vaddpd		%%ymm9 ,%%ymm11,%%ymm11		\n\t"\
	"addq	$0x240,%%rcx		\n\t"/* c1 from cc0   */"	vsubpd		%%ymm14,%%ymm8 ,%%ymm8 		\n\t"/* c1 */\
	"leaq	0x2a0(%%rbx),%%rdx	\n\t"/* c9 from isrt2 */"	vsubpd		%%ymm15,%%ymm9 ,%%ymm9 		\n\t"/* c9 */\
	/* Using the upper block(s) of the main array for temp-storage in the section below led to a nasty AVX bug to track down: */\
	/* In fermat-mod mode the 4 block addresses in ascending order are add0,1,2,3 with no 'gaps' between blocks, whereas for */\
	/* mersenne-mod the addresses in asc. order are add0,2,3,1 with a gap between contiguous-data-block pairs 0,2 and 3,1. Thus */\
	/* for fermat-mod we need [add2] as the base-address of the 'high-half' block for temp-storage; for mersenne-mod we need [add3]. */\
	/* In both cases we have that add2 < add3 so instead use (add2 - add1): > 0 for fermat-mod, < 0 for mersenne - to differentiate: */\
		"movq	%[__add2],%%rsi		\n\t"/* destroyable copy of add2 */\
		"movq	%[__add2],%%rbx		\n\t"\
		"subq	%[__add1],%%rsi		\n\t"/* rsi = (add2 - add1); if this yields a borrow (i.e. sets CF) it's mersenne, else fermat. */\
		"cmovcq %[__add3],%%rbx	\n\t" /* if CF set (i.e. h > l), copy source [add3] into dest (rbx), else leave dest = [add2]. */\
		/* rbx shared between rcol/lcol; rcx/rdx-offsets incr +0x80 in rcol for rest of block: */\
	"prefetcht1	-0x140(%%r14)\n\t"\
		"vsubpd		%%ymm4,%%ymm2,%%ymm2		\n\t		vaddpd		%%ymm14,%%ymm14,%%ymm14		\n\t"\
		"vsubpd		%%ymm5,%%ymm3,%%ymm3		\n\t		vaddpd		%%ymm15,%%ymm15,%%ymm15		\n\t"\
		"vaddpd		%%ymm4,%%ymm4,%%ymm4		\n\t		vaddpd		%%ymm8 ,%%ymm14,%%ymm14		\n\t"\
		"vaddpd		%%ymm5,%%ymm5,%%ymm5		\n\t		vaddpd		%%ymm9 ,%%ymm15,%%ymm15		\n\t"\
		"vaddpd		%%ymm2,%%ymm4,%%ymm4		\n\t		vmovaps		%%ymm8 ,0x200(%%rax)		\n\t"\
		"vaddpd		%%ymm3,%%ymm5,%%ymm5		\n\t		vmovaps		%%ymm9 ,0x220(%%rax)		\n\t"\
		"vmovaps	%%ymm2,     (%%rax)			\n\t		vmovaps		%%ymm14,%%ymm8 		\n\t"\
		"vmovaps	%%ymm3,0x020(%%rax)			\n\t		vmovaps		%%ymm15,%%ymm9 		\n\t"\
		"vmovaps	%%ymm4,%%ymm2				\n\t		vmulpd		0x100(%%rcx),%%ymm14,%%ymm14\n\t"\
		"vmovaps	%%ymm5,%%ymm3				\n\t		vmulpd		0x100(%%rcx),%%ymm15,%%ymm15\n\t"\
		"vmulpd		     (%%rcx),%%ymm4,%%ymm4	\n\t		vmulpd		0x120(%%rcx),%%ymm8 ,%%ymm8 \n\t"\
		"vmulpd		     (%%rcx),%%ymm5,%%ymm5	\n\t		vmulpd		0x120(%%rcx),%%ymm9 ,%%ymm9 \n\t"\
		"vmulpd		0x020(%%rcx),%%ymm2,%%ymm2	\n\t		vsubpd		%%ymm8      ,%%ymm15,%%ymm15\n\t"\
		"vmulpd		0x020(%%rcx),%%ymm3,%%ymm3	\n\t		vaddpd		%%ymm9      ,%%ymm14,%%ymm14\n\t"\
		"vsubpd		%%ymm2,%%ymm5,%%ymm5		\n\t		vmovaps		%%ymm15,0x060(%%rbx)		\n\t"\
		"vaddpd		%%ymm3,%%ymm4,%%ymm4		\n\t		vmovaps		%%ymm14,0x040(%%rbx)		\n\t"\
		"vmovaps	%%ymm5,0x020(%%rbx)			\n\t		vmovaps		0x200(%%rax),%%ymm14		\n\t"\
		"vmovaps	%%ymm4,     (%%rbx)			\n\t		vmovaps		0x220(%%rax),%%ymm15		\n\t"\
		"vmovaps		 (%%rax),%%ymm4			\n\t		vmovaps		%%ymm14,%%ymm8 		\n\t"\
		"vmovaps	0x020(%%rax),%%ymm5			\n\t		vmovaps		%%ymm15,%%ymm9 		\n\t"\
		"vmovaps	%%ymm4,%%ymm2				\n\t		vmulpd		0x100(%%rdx),%%ymm14,%%ymm14\n\t"\
		"vmovaps	%%ymm5,%%ymm3				\n\t		vmulpd		0x100(%%rdx),%%ymm15,%%ymm15\n\t"\
		"vmulpd		     (%%rdx),%%ymm4,%%ymm4	\n\t		vmulpd		0x120(%%rdx),%%ymm8 ,%%ymm8 \n\t"\
		"vmulpd		     (%%rdx),%%ymm5,%%ymm5	\n\t		vmulpd		0x120(%%rdx),%%ymm9 ,%%ymm9 \n\t"\
		"vmulpd		0x020(%%rdx),%%ymm2,%%ymm2	\n\t		vsubpd		%%ymm8      ,%%ymm15,%%ymm15\n\t"\
		"vmulpd		0x020(%%rdx),%%ymm3,%%ymm3	\n\t		vaddpd		%%ymm9      ,%%ymm14,%%ymm14\n\t"\
		"vsubpd		%%ymm2,%%ymm5,%%ymm5		\n\t		vmovaps		%%ymm15,0x160(%%rbx)		\n\t"\
		"vaddpd		%%ymm3,%%ymm4,%%ymm4		\n\t		vmovaps		%%ymm14,0x140(%%rbx)		\n\t"\
		"vmovaps	%%ymm5,0x120(%%rbx)			\n\t		vsubpd		%%ymm13,%%ymm10,%%ymm10		\n\t"\
		"vmovaps	%%ymm4,0x100(%%rbx)			\n\t		vsubpd		%%ymm12,%%ymm11,%%ymm11		\n\t"\
		"addq	$0x080,%%rcx	\n\t"/* c1  from c1 */	"	vaddpd		%%ymm13,%%ymm13,%%ymm13 	\n\t"/* c7  */\
		"addq	$0x080,%%rdx	\n\t"/* c9  from c9 */	"	vaddpd		%%ymm12,%%ymm12,%%ymm12 	\n\t"/* c15 */\
		"vsubpd		%%ymm7,%%ymm0,%%ymm0		\n\t		vaddpd		%%ymm10,%%ymm13,%%ymm13		\n\t"\
		"vsubpd		%%ymm6,%%ymm1,%%ymm1		\n\t		vaddpd		%%ymm11,%%ymm12,%%ymm12		\n\t"\
		"vaddpd		%%ymm7,%%ymm7,%%ymm7		\n\t		vmovaps		%%ymm13,%%ymm8 		\n\t"\
		"vaddpd		%%ymm6,%%ymm6,%%ymm6		\n\t		vmovaps		%%ymm11,%%ymm9 		\n\t"\
		"vaddpd		%%ymm0,%%ymm7,%%ymm7		\n\t		vmulpd		0x100(%%rcx),%%ymm13,%%ymm13\n\t"\
		"vaddpd		%%ymm1,%%ymm6,%%ymm6		\n\t		vmulpd		0x100(%%rcx),%%ymm11,%%ymm11\n\t"\
		"vmovaps	%%ymm7,%%ymm4				\n\t		vmulpd		0x120(%%rcx),%%ymm8 ,%%ymm8 \n\t"\
		"vmovaps	%%ymm1,%%ymm5				\n\t		vmulpd		0x120(%%rcx),%%ymm9 ,%%ymm9 \n\t"\
		"vmulpd		     (%%rcx),%%ymm7,%%ymm7	\n\t		vsubpd		%%ymm8      ,%%ymm11,%%ymm11\n\t"\
		"vmulpd		     (%%rcx),%%ymm1,%%ymm1	\n\t		vaddpd		%%ymm9      ,%%ymm13,%%ymm13\n\t"\
		"vmulpd		0x020(%%rcx),%%ymm4,%%ymm4	\n\t		vmovaps		%%ymm11,0x0e0(%%rbx)		\n\t"\
		"vmulpd		0x020(%%rcx),%%ymm5,%%ymm5	\n\t		vmovaps		%%ymm13,0x0c0(%%rbx)		\n\t"\
		"vsubpd		%%ymm4,%%ymm1,%%ymm1		\n\t		vmovaps		%%ymm10,%%ymm8 		\n\t"\
		"vaddpd		%%ymm5,%%ymm7,%%ymm7		\n\t		vmovaps		%%ymm12,%%ymm9 		\n\t"\
		"vmovaps	%%ymm1,0x0a0(%%rbx)			\n\t		vmulpd		0x100(%%rdx),%%ymm10,%%ymm10\n\t"\
		"vmovaps	%%ymm7,0x080(%%rbx)			\n\t		vmulpd		0x100(%%rdx),%%ymm12,%%ymm12\n\t"\
		"vmovaps	%%ymm0,%%ymm4				\n\t		vmulpd		0x120(%%rdx),%%ymm8 ,%%ymm8 \n\t"\
		"vmovaps	%%ymm6,%%ymm5				\n\t		vmulpd		0x120(%%rdx),%%ymm9 ,%%ymm9 \n\t"\
		"vmulpd		     (%%rdx),%%ymm0,%%ymm0	\n\t		vsubpd		%%ymm8      ,%%ymm12,%%ymm12\n\t"\
		"vmulpd		     (%%rdx),%%ymm6,%%ymm6	\n\t		vaddpd		%%ymm9      ,%%ymm10,%%ymm10\n\t"\
		"vmulpd		0x020(%%rdx),%%ymm4,%%ymm4	\n\t		vmovaps		%%ymm12,0x1e0(%%rbx)		\n\t"\
		"vmulpd		0x020(%%rdx),%%ymm5,%%ymm5	\n\t		vmovaps		%%ymm10,0x1c0(%%rbx)		\n\t"\
		"vsubpd		%%ymm4,%%ymm6,%%ymm6		\n\t"		/*...Block 3: t5,13,21,29 -> r17,21,19,23: */\
		"vaddpd		%%ymm5,%%ymm0,%%ymm0		\n\t		movq	%[__r1],%%rax		\n\t"/* r17 in rcol */\
		"vmovaps	%%ymm6,0x1a0(%%rbx)			\n\t		movq		%[__isrt2],%%rdi		\n\t"\
		"vmovaps	%%ymm0,0x180(%%rbx)			\n\t		vmovaps		(%%rdi),%%ymm10			\n\t"\
	"prefetcht1	-0x180(%%r14)					\n\t		vmovaps		0x240(%%rax),%%ymm12	\n\t"\
		/*...Block 1: t1,9,17,25 -> r1,5,3,7: */		"	vmovaps		0x260(%%rax),%%ymm13		\n\t"\
		"vmovaps	     (%%rax),%%ymm0			\n\t		vmovaps		0x2c0(%%rax),%%ymm8 		\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1			\n\t		vmovaps		0x2e0(%%rax),%%ymm9 		\n\t"\
		"vmovaps	0x080(%%rax),%%ymm2			\n\t		vaddpd		0x260(%%rax),%%ymm12,%%ymm12\n\t"\
		"vmovaps	0x0a0(%%rax),%%ymm3			\n\t		vsubpd		0x240(%%rax),%%ymm13,%%ymm13\n\t"\
		"vsubpd		0x080(%%rax),%%ymm0,%%ymm0	\n\t		vsubpd		0x2e0(%%rax),%%ymm8 ,%%ymm8 \n\t"\
		"vsubpd		0x0a0(%%rax),%%ymm1,%%ymm1	\n\t		vaddpd		0x2c0(%%rax),%%ymm9 ,%%ymm9 \n\t"\
		"vaddpd		     (%%rax),%%ymm2,%%ymm2	\n\t		vmulpd		%%ymm10,%%ymm12,%%ymm12		\n\t"\
		"vaddpd		0x020(%%rax),%%ymm3,%%ymm3	\n\t		vmulpd		%%ymm10,%%ymm13,%%ymm13		\n\t"\
		"vmovaps	0x040(%%rax),%%ymm4			\n\t		vmulpd		%%ymm10,%%ymm8 ,%%ymm8 		\n\t"\
		"vmovaps	0x060(%%rax),%%ymm5			\n\t		vmulpd		%%ymm10,%%ymm9 ,%%ymm9 		\n\t"\
		"vmovaps	0x0c0(%%rax),%%ymm6			\n\t		vmovaps		%%ymm12,%%ymm14		\n\t"\
		"vmovaps	0x0e0(%%rax),%%ymm7			\n\t		vmovaps		%%ymm13,%%ymm15		\n\t"\
		"vsubpd		0x0c0(%%rax),%%ymm4,%%ymm4	\n\t		vsubpd		%%ymm8 ,%%ymm12,%%ymm12		\n\t"\
		"vsubpd		0x0e0(%%rax),%%ymm5,%%ymm5	\n\t		vsubpd		%%ymm9 ,%%ymm13,%%ymm13		\n\t"\
		"vaddpd		0x040(%%rax),%%ymm6,%%ymm6	\n\t		vaddpd		%%ymm8 ,%%ymm14,%%ymm14		\n\t"\
		"vaddpd		0x060(%%rax),%%ymm7,%%ymm7	\n\t		vaddpd		%%ymm9 ,%%ymm15,%%ymm15		\n\t"\
		"movq	%%rbx,%%rcx	\n\t"/* rbx saves 'high-half' block address; move that to rcx and put 'low-half' block address into rbx */\
		"movq	%[__add0],%%rbx					\n\t		vmovaps		0x200(%%rax),%%ymm8 		\n\t"\
		"subq	$0x280,%%rdx	/* c8 from c13 */\n\t		vmovaps		0x220(%%rax),%%ymm9 		\n\t"\
		"vaddpd		%%ymm6,%%ymm2,%%ymm2		\n\t		vmovaps		0x280(%%rax),%%ymm10		\n\t"\
		"vaddpd		%%ymm7,%%ymm3,%%ymm3		\n\t		vmovaps		0x2a0(%%rax),%%ymm11		\n\t"\
		"vmovaps	%%ymm2,     (%%rbx)			\n\t		vsubpd		0x2a0(%%rax),%%ymm8 ,%%ymm8 \n\t"\
		"vmovaps	%%ymm3,0x020(%%rbx)			\n\t		vsubpd		0x280(%%rax),%%ymm9 ,%%ymm9 \n\t"\
		"vaddpd		%%ymm6,%%ymm6,%%ymm6		\n\t		vaddpd		0x200(%%rax),%%ymm11,%%ymm11\n\t"\
		"vaddpd		%%ymm7,%%ymm7,%%ymm7		\n\t		vaddpd		0x220(%%rax),%%ymm10,%%ymm10\n\t"\
		"vsubpd		%%ymm6,%%ymm2,%%ymm2		\n\t		vsubpd		%%ymm12,%%ymm11,%%ymm11		\n\t"\
		"vsubpd		%%ymm7,%%ymm3,%%ymm3		\n\t		vsubpd		%%ymm13,%%ymm9 ,%%ymm9 		\n\t"\
		"vmovaps	%%ymm2,%%ymm6				\n\t		vaddpd		%%ymm12,%%ymm12,%%ymm12		\n\t"\
		"vmovaps	%%ymm3,%%ymm7				\n\t		vaddpd		%%ymm13,%%ymm13,%%ymm13		\n\t"\
		"vmulpd		     (%%rdx),%%ymm2,%%ymm2	\n\t		vaddpd		%%ymm11,%%ymm12,%%ymm12		\n\t"\
		"vmulpd		     (%%rdx),%%ymm3,%%ymm3	\n\t		vaddpd		%%ymm9 ,%%ymm13,%%ymm13		\n\t"\
		"vmulpd		0x020(%%rdx),%%ymm6,%%ymm6	\n\t		vmovaps		%%ymm11,     (%%rax)		\n\t"\
		"vmulpd		0x020(%%rdx),%%ymm7,%%ymm7	\n\t		vmovaps		%%ymm9 ,0x020(%%rax)		\n\t"\
		"vsubpd		%%ymm6,%%ymm3,%%ymm3		\n\t		vmovaps		%%ymm12,%%ymm11		\n\t"\
		"vaddpd		%%ymm7,%%ymm2,%%ymm2		\n\t		vmovaps		%%ymm13,%%ymm9 		\n\t"\
/*==========================*/"\n\t"\
		"													vmulpd		0x0c0(%%rdx),%%ymm12,%%ymm12/* c2 */		\n\t"\
		"													vmulpd		0x0c0(%%rdx),%%ymm13,%%ymm13\n\t"\
		"													vmulpd		0x0e0(%%rdx),%%ymm11,%%ymm11\n\t"\
		"													vmulpd		0x0e0(%%rdx),%%ymm9 ,%%ymm9 \n\t"\
		"													vsubpd		%%ymm11     ,%%ymm13,%%ymm13\n\t"\
		"													vaddpd		%%ymm9      ,%%ymm12,%%ymm12\n\t"\
"vmovaps	0x120(%%rcx),%%ymm7		\n\t"\
"vmovaps	0x100(%%rcx),%%ymm6		\n\t"\
"vmovaps	%%ymm3,0x120(%%rax)	/* r9 */\n\t"\
"vmovaps	%%ymm2,0x100(%%rax)	/* r8 */\n\t"\
"vmovaps	%%ymm7,0x320(%%rax)	/* r25 */\n\t"\
"vmovaps	%%ymm6,0x300(%%rax)	/* r24 */\n\t"\
															"vmovaps	0x060(%%rcx),%%ymm11	\n\t"\
															"vmovaps	0x040(%%rcx),%%ymm9		\n\t"\
															"vmovaps	%%ymm13,0x060(%%rax)	/* r3 */\n\t"\
															"vmovaps	%%ymm12,0x040(%%rax)	/* r2 */\n\t"\
															"vmovaps	%%ymm11,0x260(%%rax)	/* r19 */\n\t"\
															"vmovaps	%%ymm9 ,0x240(%%rax)	/* r18 */\n\t"\
		"													addq	$0x040,%%rdx		/* c4 in lcol; c10 in rcol*/		\n\t"\
		"													vmovaps		     (%%rax),%%ymm12		\n\t"\
		"													vmovaps		0x020(%%rax),%%ymm13		\n\t"\
		"													vmovaps		%%ymm12,%%ymm11		\n\t"\
		"													vmovaps		%%ymm13,%%ymm9 		\n\t"\
		"													vmulpd		0x0c0(%%rdx),%%ymm12,%%ymm12		\n\t"\
		"													vmulpd		0x0c0(%%rdx),%%ymm13,%%ymm13		\n\t"\
		"													vmulpd		0x0e0(%%rdx),%%ymm11,%%ymm11		\n\t"\
		"													vmulpd		0x0e0(%%rdx),%%ymm9 ,%%ymm9 		\n\t"\
"vmovaps	0x020(%%rbx),%%ymm3		\n\t"\
"vmovaps	     (%%rbx),%%ymm2		\n\t"\
"vmovaps	0x020(%%rcx),%%ymm7		\n\t"\
"vmovaps	     (%%rcx),%%ymm6		\n\t"\
"vmovaps	%%ymm3,0x020(%%rax)	/* r1 */\n\t"\
"vmovaps	%%ymm2,     (%%rax)	/* r0 */\n\t"\
"vmovaps	%%ymm7,0x220(%%rax)	/* r17 */\n\t"\
"vmovaps	%%ymm6,0x200(%%rax)	/* r16 */\n\t"\
		"vaddpd		%%ymm5,%%ymm0,%%ymm0		\n\t		vsubpd		%%ymm11     ,%%ymm13,%%ymm13		\n\t"\
		"vsubpd		%%ymm4,%%ymm1,%%ymm1		\n\t		vaddpd		%%ymm9      ,%%ymm12,%%ymm12		\n\t"\
		"vmovaps	%%ymm0,%%ymm2				\n\t"\
		"vmovaps	%%ymm1,%%ymm3				\n\t"\
		"vaddpd		%%ymm5,%%ymm5,%%ymm5		\n\t"\
		"vaddpd		%%ymm4,%%ymm4,%%ymm4		\n\t"\
		"vmovaps	%%ymm0,%%ymm6				\n\t"\
		"vmovaps	%%ymm1,%%ymm7				\n\t"\
		"vmulpd		     (%%rdx),%%ymm2,%%ymm2	\n\t"\
		"vmulpd		     (%%rdx),%%ymm3,%%ymm3	\n\t"\
		"vmulpd		0x020(%%rdx),%%ymm6,%%ymm6	\n\t"\
		"vmulpd		0x020(%%rdx),%%ymm7,%%ymm7	\n\t"\
	"prefetcht1	-0x1c0(%%r14)\n\t"\
															"vmovaps	0x160(%%rcx),%%ymm11	\n\t"\
															"vmovaps	0x140(%%rcx),%%ymm9		\n\t"\
															"vmovaps	%%ymm13,0x160(%%rax)	/* r11 */\n\t"\
															"vmovaps	%%ymm12,0x140(%%rax)	/* r10 */\n\t"\
															"vmovaps	%%ymm11,0x360(%%rax)	/* r27 */\n\t"\
															"vmovaps	%%ymm9 ,0x340(%%rax)	/* r26 */\n\t"\
		"addq	$0x040,%%rdx		/* c12 in lcol; c6 in rcol*/\n\t		vsubpd		%%ymm15,%%ymm8 ,%%ymm8 		\n\t"\
		"vsubpd		%%ymm6,%%ymm3,%%ymm3		\n\t		vsubpd		%%ymm14,%%ymm10,%%ymm10		\n\t"\
		"vaddpd		%%ymm7,%%ymm2,%%ymm2		\n\t		vaddpd		%%ymm15,%%ymm15,%%ymm15		\n\t"\
		"													vaddpd		%%ymm14,%%ymm14,%%ymm14		\n\t"\
		"													vaddpd		%%ymm8 ,%%ymm15,%%ymm15		\n\t"\
		"													vaddpd		%%ymm10,%%ymm14,%%ymm14		\n\t"\
		"													vmovaps		%%ymm15,%%ymm12		\n\t"\
		"													vmovaps		%%ymm10,%%ymm13		\n\t"\
		"													vmulpd		0x0c0(%%rdx),%%ymm15,%%ymm15		\n\t"\
		"													vmulpd		0x0c0(%%rdx),%%ymm10,%%ymm10		\n\t"\
		"													vmulpd		0x0e0(%%rdx),%%ymm12,%%ymm12		\n\t"\
		"													vmulpd		0x0e0(%%rdx),%%ymm13,%%ymm13		\n\t"\
		"													vsubpd		%%ymm12     ,%%ymm10,%%ymm10		\n\t"\
"vmovaps	0x0a0(%%rcx),%%ymm7		\n\t"\
"vmovaps	0x080(%%rcx),%%ymm6		\n\t"\
"vmovaps	%%ymm3,0x0a0(%%rax)	/* r5 */\n\t"\
"vmovaps	%%ymm2,0x080(%%rax)	/* r4 */\n\t"\
"vmovaps	%%ymm7,0x2a0(%%rax)	/* r21 */\n\t"\
"vmovaps	%%ymm6,0x280(%%rax)	/* r20 */\n\t"\
		"vsubpd		%%ymm5,%%ymm0,%%ymm0		\n\t		vaddpd		%%ymm13     ,%%ymm15,%%ymm15		\n\t"\
		"vaddpd		%%ymm4,%%ymm1,%%ymm1		\n\t"\
		"vmovaps	%%ymm0,%%ymm6				\n\t"\
		"vmovaps	%%ymm1,%%ymm7				\n\t"\
		"vmulpd		     (%%rdx),%%ymm0,%%ymm0	\n\t"\
		"vmulpd		     (%%rdx),%%ymm1,%%ymm1	\n\t"\
		"vmulpd		0x020(%%rdx),%%ymm6,%%ymm6	\n\t"\
		"vmulpd		0x020(%%rdx),%%ymm7,%%ymm7	\n\t"\
		"vsubpd		%%ymm6,%%ymm1,%%ymm1		\n\t"\
		"vaddpd		%%ymm7,%%ymm0,%%ymm0		\n\t"\
															"vmovaps	0x0e0(%%rcx),%%ymm13	\n\t"\
															"vmovaps	0x0c0(%%rcx),%%ymm12	\n\t"\
															"vmovaps	%%ymm10,0x0e0(%%rax)	/* r7 */\n\t"\
															"vmovaps	%%ymm15,0x0c0(%%rax)	/* r6 */\n\t"\
															"vmovaps	%%ymm13,0x2e0(%%rax)	/* r23 */\n\t"\
															"vmovaps	%%ymm12,0x2c0(%%rax)	/* r22 */\n\t"\
		"													vmovaps		%%ymm8 ,%%ymm12		\n\t"\
		"													vmovaps		%%ymm14,%%ymm13		\n\t"\
		"													vmulpd		0x100(%%rdx),%%ymm8 ,%%ymm8 	/* c14 */		\n\t"\
		"													vmulpd		0x100(%%rdx),%%ymm14,%%ymm14		\n\t"\
		"													vmulpd		0x120(%%rdx),%%ymm12,%%ymm12		\n\t"\
		"													vmulpd		0x120(%%rdx),%%ymm13,%%ymm13		\n\t"\
		"													vsubpd		%%ymm12     ,%%ymm14,%%ymm14		\n\t"\
		"													vaddpd		%%ymm13     ,%%ymm8 ,%%ymm8 		\n\t"\
"vmovaps	0x1a0(%%rcx),%%ymm7				\n\t"\
"vmovaps	0x180(%%rcx),%%ymm6				\n\t"\
"vmovaps	%%ymm1,0x1a0(%%rax)	/* r13 */	\n\t"\
"vmovaps	%%ymm0,0x180(%%rax)	/* r12 */	\n\t"\
"vmovaps	%%ymm7,0x3a0(%%rax)	/* r29 */	\n\t"\
"vmovaps	%%ymm6,0x380(%%rax)	/* r28 */	\n\t"\
															"vmovaps	0x1e0(%%rcx),%%ymm13	\n\t"\
															"vmovaps	0x1c0(%%rcx),%%ymm12	\n\t"\
															"vmovaps	%%ymm14,0x1e0(%%rax)	/* r15 */\n\t"\
															"vmovaps	%%ymm8 ,0x1c0(%%rax)	/* r14 */\n\t"\
															"vmovaps	%%ymm13,0x3e0(%%rax)	/* r31 */\n\t"\
															"vmovaps	%%ymm12,0x3c0(%%rax)	/* r30 */\n\t"\
	"prefetcht1	-0x200(%%r14)\n\t"\
/*==========================*/"\n\t"\
	/**** Finish with 4-way 'un'terleaving: ****/\
		"movq	%[__r1] ,%%rsi	\n\t"\
		"movq	%[__add0],%%rax\n\t"\
		"movq	%[__add1],%%rbx\n\t"\
		"movq	%[__add2],%%rcx\n\t"\
		"movq	%[__add3],%%rdx\n\t"\
	"/* a[j+p0]: Inputs from r1 +0/1, 16/17, 2/3, 18/19. Outputs into [add0,add0+0x100,add1,add1+0x100]+0x0: */	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps	     (%%rsi),%%ymm1						\n\t		vmovaps	0x020(%%rsi),%%ymm3								\n\t"\
		"vmovaps	0x200(%%rsi),%%ymm5						\n\t		vmovaps	0x220(%%rsi),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps	0x040(%%rsi),%%ymm6						\n\t		vmovaps	0x060(%%rsi),%%ymm14							\n\t"\
		"vmovaps	0x240(%%rsi),%%ymm5						\n\t		vmovaps	0x260(%%rsi),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm6 ,%%ymm0				\n\t		vshufpd	$15 ,%%ymm13,%%ymm14 ,%%ymm2					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm6 ,%%ymm6				\n\t		vshufpd	$0  ,%%ymm13,%%ymm14 ,%%ymm14					\n\t"\
		"vperm2f128 $32 ,%%ymm0 ,%%ymm7 ,%%ymm5				\n\t		vperm2f128 $32 ,%%ymm2 ,%%ymm15,%%ymm13		/* outB	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm0 ,%%ymm7 ,%%ymm7				\n\t		vperm2f128 $49 ,%%ymm2 ,%%ymm15,%%ymm15		/* outD	*/	\n\t"\
		"vperm2f128 $32 ,%%ymm6 ,%%ymm1 ,%%ymm0				\n\t		vperm2f128 $32 ,%%ymm14,%%ymm3 ,%%ymm2 		/* outA	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm6 ,%%ymm1 ,%%ymm1				\n\t		vperm2f128 $49 ,%%ymm14,%%ymm3 ,%%ymm3 		/* outC	*/	\n\t"\
		"vmovaps 	%%ymm5 ,(%%rbx)							\n\t		vmovaps %%ymm13,0x020(%%rbx)				/* outB	*/	\n\t"\
		"vmovaps 	%%ymm7 ,(%%rdx)							\n\t		vmovaps %%ymm15,0x020(%%rdx)				/* outD	*/	\n\t"\
		"vmovaps 	%%ymm0 ,(%%rax)							\n\t		vmovaps %%ymm2 ,0x020(%%rax)				/* outA	*/	\n\t"\
		"vmovaps 	%%ymm1 ,(%%rcx)							\n\t		vmovaps %%ymm3 ,0x020(%%rcx)				/* outC	*/	\n\t"\
	"/* a[j+p2]: Inputs from r3 +0/1, 16/17, 2/3, 18/19. Outputs into [add0,add0+0x100,add1,add1+0x100]+0x40: */		\n\t"\
		"addq	$0x80,%%rsi	\n\t"\
		"addq	$0x40,%%rax	\n\t"\
		"addq	$0x40,%%rbx	\n\t"\
		"addq	$0x40,%%rcx	\n\t"\
		"addq	$0x40,%%rdx	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps	     (%%rsi),%%ymm1						\n\t		vmovaps	0x020(%%rsi),%%ymm3								\n\t"\
		"vmovaps	0x200(%%rsi),%%ymm5						\n\t		vmovaps	0x220(%%rsi),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps	0x040(%%rsi),%%ymm6						\n\t		vmovaps	0x060(%%rsi),%%ymm14							\n\t"\
		"vmovaps	0x240(%%rsi),%%ymm5						\n\t		vmovaps	0x260(%%rsi),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm6 ,%%ymm0				\n\t		vshufpd	$15 ,%%ymm13,%%ymm14 ,%%ymm2					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm6 ,%%ymm6				\n\t		vshufpd	$0  ,%%ymm13,%%ymm14 ,%%ymm14					\n\t"\
		"vperm2f128 $32 ,%%ymm0 ,%%ymm7 ,%%ymm5				\n\t		vperm2f128 $32 ,%%ymm2 ,%%ymm15,%%ymm13		/* outB	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm0 ,%%ymm7 ,%%ymm7				\n\t		vperm2f128 $49 ,%%ymm2 ,%%ymm15,%%ymm15		/* outD	*/	\n\t"\
		"vperm2f128 $32 ,%%ymm6 ,%%ymm1 ,%%ymm0				\n\t		vperm2f128 $32 ,%%ymm14,%%ymm3 ,%%ymm2 		/* outA	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm6 ,%%ymm1 ,%%ymm1				\n\t		vperm2f128 $49 ,%%ymm14,%%ymm3 ,%%ymm3 		/* outC	*/	\n\t"\
		"vmovaps 	%%ymm5 ,(%%rbx)							\n\t		vmovaps %%ymm13,0x020(%%rbx)				/* outB	*/	\n\t"\
		"vmovaps 	%%ymm7 ,(%%rdx)							\n\t		vmovaps %%ymm15,0x020(%%rdx)				/* outD	*/	\n\t"\
		"vmovaps 	%%ymm0 ,(%%rax)							\n\t		vmovaps %%ymm2 ,0x020(%%rax)				/* outA	*/	\n\t"\
		"vmovaps 	%%ymm1 ,(%%rcx)							\n\t		vmovaps %%ymm3 ,0x020(%%rcx)				/* outC	*/	\n\t"\
	"/* a[j+p4]: Inputs from r5 +0/1, 16/17, 2/3, 18/19. Outputs into [add0,add0+0x100,add1,add1+0x100]+0x80: */		\n\t"\
		"addq	$0x80,%%rsi	\n\t"\
		"addq	$0x40,%%rax	\n\t"\
		"addq	$0x40,%%rbx	\n\t"\
		"addq	$0x40,%%rcx	\n\t"\
		"addq	$0x40,%%rdx	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps	     (%%rsi),%%ymm1						\n\t		vmovaps	0x020(%%rsi),%%ymm3								\n\t"\
		"vmovaps	0x200(%%rsi),%%ymm5						\n\t		vmovaps	0x220(%%rsi),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps	0x040(%%rsi),%%ymm6						\n\t		vmovaps	0x060(%%rsi),%%ymm14							\n\t"\
		"vmovaps	0x240(%%rsi),%%ymm5						\n\t		vmovaps	0x260(%%rsi),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm6 ,%%ymm0				\n\t		vshufpd	$15 ,%%ymm13,%%ymm14 ,%%ymm2					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm6 ,%%ymm6				\n\t		vshufpd	$0  ,%%ymm13,%%ymm14 ,%%ymm14					\n\t"\
		"vperm2f128 $32 ,%%ymm0 ,%%ymm7 ,%%ymm5				\n\t		vperm2f128 $32 ,%%ymm2 ,%%ymm15,%%ymm13		/* outB	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm0 ,%%ymm7 ,%%ymm7				\n\t		vperm2f128 $49 ,%%ymm2 ,%%ymm15,%%ymm15		/* outD	*/	\n\t"\
		"vperm2f128 $32 ,%%ymm6 ,%%ymm1 ,%%ymm0				\n\t		vperm2f128 $32 ,%%ymm14,%%ymm3 ,%%ymm2 		/* outA	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm6 ,%%ymm1 ,%%ymm1				\n\t		vperm2f128 $49 ,%%ymm14,%%ymm3 ,%%ymm3 		/* outC	*/	\n\t"\
		"vmovaps 	%%ymm5 ,(%%rbx)							\n\t		vmovaps %%ymm13,0x020(%%rbx)				/* outB	*/	\n\t"\
		"vmovaps 	%%ymm7 ,(%%rdx)							\n\t		vmovaps %%ymm15,0x020(%%rdx)				/* outD	*/	\n\t"\
		"vmovaps 	%%ymm0 ,(%%rax)							\n\t		vmovaps %%ymm2 ,0x020(%%rax)				/* outA	*/	\n\t"\
		"vmovaps 	%%ymm1 ,(%%rcx)							\n\t		vmovaps %%ymm3 ,0x020(%%rcx)				/* outC	*/	\n\t"\
	"/* a[j+p6]: Inputs from r7 +0/1, 16/17, 2/3, 18/19. Outputs into [add0,add0+0x100,add1,add1+0x100]+0xc0: */		\n\t"\
		"addq	$0x80,%%rsi	\n\t"\
		"addq	$0x40,%%rax	\n\t"\
		"addq	$0x40,%%rbx	\n\t"\
		"addq	$0x40,%%rcx	\n\t"\
		"addq	$0x40,%%rdx	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps	     (%%rsi),%%ymm1						\n\t		vmovaps	0x020(%%rsi),%%ymm3								\n\t"\
		"vmovaps	0x200(%%rsi),%%ymm5						\n\t		vmovaps	0x220(%%rsi),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps	0x040(%%rsi),%%ymm6						\n\t		vmovaps	0x060(%%rsi),%%ymm14							\n\t"\
		"vmovaps	0x240(%%rsi),%%ymm5						\n\t		vmovaps	0x260(%%rsi),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm6 ,%%ymm0				\n\t		vshufpd	$15 ,%%ymm13,%%ymm14 ,%%ymm2					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm6 ,%%ymm6				\n\t		vshufpd	$0  ,%%ymm13,%%ymm14 ,%%ymm14					\n\t"\
		"vperm2f128 $32 ,%%ymm0 ,%%ymm7 ,%%ymm5				\n\t		vperm2f128 $32 ,%%ymm2 ,%%ymm15,%%ymm13		/* outB	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm0 ,%%ymm7 ,%%ymm7				\n\t		vperm2f128 $49 ,%%ymm2 ,%%ymm15,%%ymm15		/* outD	*/	\n\t"\
		"vperm2f128 $32 ,%%ymm6 ,%%ymm1 ,%%ymm0				\n\t		vperm2f128 $32 ,%%ymm14,%%ymm3 ,%%ymm2 		/* outA	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm6 ,%%ymm1 ,%%ymm1				\n\t		vperm2f128 $49 ,%%ymm14,%%ymm3 ,%%ymm3 		/* outC	*/	\n\t"\
		"vmovaps 	%%ymm5 ,(%%rbx)							\n\t		vmovaps %%ymm13,0x020(%%rbx)				/* outB	*/	\n\t"\
		"vmovaps 	%%ymm7 ,(%%rdx)							\n\t		vmovaps %%ymm15,0x020(%%rdx)				/* outD	*/	\n\t"\
		"vmovaps 	%%ymm0 ,(%%rax)							\n\t		vmovaps %%ymm2 ,0x020(%%rax)				/* outA	*/	\n\t"\
		"vmovaps 	%%ymm1 ,(%%rcx)							\n\t		vmovaps %%ymm3 ,0x020(%%rcx)				/* outC	*/	\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__add1] "m" (Xadd1)\
		 ,[__add2] "m" (Xadd2)\
		 ,[__add3] "m" (Xadd3)\
		 ,[__r1] "m" (Xr1)\
		 ,[__r9] "m" (Xr9)\
		 ,[__r17] "m" (Xr17)\
		 ,[__r25] "m" (Xr25)\
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__cc0] "m" (Xcc0)\
		 ,[__c1] "m" (Xc1)\
		 ,[__c2] "m" (Xc2)\
		 ,[__c3] "m" (Xc3)\
		 ,[__c4] "m" (Xc4)\
		 ,[__c5] "m" (Xc5)\
		 ,[__c6] "m" (Xc6)\
		 ,[__c7] "m" (Xc7)\
		 ,[__c8] "m" (Xc8)\
		 ,[__c9] "m" (Xc9)\
		 ,[__c10] "m" (Xc10)\
		 ,[__c11] "m" (Xc11)\
		 ,[__c12] "m" (Xc12)\
		 ,[__c13] "m" (Xc13)\
		 ,[__c14] "m" (Xc14)\
		 ,[__c15] "m" (Xc15)\
	 ,[__pfetch_dist] "m" (Xpfetch_dist)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r10","r11","r12","r13","r14","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
	);\
	}

#elif defined(USE_SSE2)

	// Since the add0/1 block addresses advance in opposite directions, our prefetch scheme here is like so:
	// In fwd DIF prefetch  ahead 8 cache lines' worth of data w.r.to add0;
	// In inv DIT prefetch behind 8 cache lines' worth of data w.r.to add1;
	#define SSE2_RADIX16_WRAPPER_DIF(Xadd0,Xadd1,Xr1,Xr9,Xr17,Xr25,Xisrt2,Xcc0,Xc1,Xc2,Xc3,Xc4,Xc5,Xc6,Xc7,Xc8,Xc9,Xc10,Xc11,Xc12,Xc13,Xc14,Xc15)\
	{\
	__asm__ volatile (\
	"movq	%[__add0],%%r14						\n\t"\
	"prefetcht1	0x40(%%r14)\n\t"\
	/*************************************************************/\
	/*                  1st set of inputs:                       */\
	/*************************************************************/\
		"/*...Block 1: */										\n\t	/*...Block 2: */					\n\t"\
		"movq		%[__add0],%%rax								\n\t	movq		%[__r9] ,%%r10			\n\t"\
		"movq		%[__add1],%%rbx								\n\t	movq		%[__c2] ,%%r11			\n\t"\
		"movq		%[__r1] ,%%rcx								\n\t	movaps		0x20(%%rax),%%xmm14		\n\t"\
		"movq		%[__c4] ,%%rdx								\n\t	movaps		0x20(%%rax),%%xmm8		\n\t"\
		"movaps		0x40(%%rax),%%xmm6							\n\t	unpckhpd	0x20(%%rbx),%%xmm14		\n\t"\
		"movaps		0x40(%%rbx),%%xmm2							\n\t	unpcklpd	0x20(%%rbx),%%xmm8		\n\t"\
		"movaps		%%xmm6,%%xmm0								\n\t	movaps			%%xmm14,0x100(%%r10)\n\t"\
		"movaps		%%xmm2,%%xmm3								\n\t	movaps		0x30(%%rax),%%xmm15		\n\t"\
		"unpckhpd	%%xmm2,%%xmm6								\n\t	movaps		0x30(%%rax),%%xmm9		\n\t"\
		"unpcklpd	%%xmm3,%%xmm0								\n\t	unpckhpd	0x30(%%rbx),%%xmm15		\n\t"\
		"movaps		%%xmm6,0x140(%%rcx)							\n\t	unpcklpd	0x30(%%rbx),%%xmm9		\n\t"\
		"movaps		0x50(%%rax),%%xmm7							\n\t	movaps		%%xmm15,0x110(%%r10)	\n\t"\
		"movaps		0x50(%%rbx),%%xmm4							\n\t	movaps		%%xmm8,%%xmm10			\n\t"\
		"movaps		%%xmm7,%%xmm1								\n\t	movaps		%%xmm9,%%xmm11			\n\t"\
		"movaps		%%xmm4,%%xmm5								\n\t	mulpd		    (%%r11),%%xmm8		\n\t"\
		"unpckhpd	%%xmm4,%%xmm7								\n\t	mulpd		    (%%r11),%%xmm9		\n\t"\
		"unpcklpd	%%xmm5,%%xmm1								\n\t	mulpd		0x10(%%r11),%%xmm10		\n\t"\
		"movaps		%%xmm7,0x150(%%rcx)							\n\t	mulpd		0x10(%%r11),%%xmm11		\n\t"\
		"movaps		%%xmm0,%%xmm2								\n\t	addpd		%%xmm10,%%xmm9			\n\t"\
		"movaps		%%xmm1,%%xmm3								\n\t	subpd		%%xmm11,%%xmm8			\n\t"\
		"/* Rest identical to code in radix16_dif_pass: */		\n\t	movaps		%%xmm9,%%xmm11			\n\t"\
		"mulpd		    (%%rdx),%%xmm0							\n\t	movaps		%%xmm8,%%xmm10			\n\t"\
		"mulpd		    (%%rdx),%%xmm1							\n\t	movq		%[__c10],%%r11			\n\t"\
		"mulpd		0x10(%%rdx),%%xmm2							\n\t	movaps		0xa0(%%rax),%%xmm14		\n\t"\
		"mulpd		0x10(%%rdx),%%xmm3							\n\t	movaps		0xa0(%%rax),%%xmm12		\n\t"\
		"addpd		%%xmm2,%%xmm1								\n\t	unpckhpd	0xa0(%%rbx),%%xmm14		\n\t"\
		"subpd		%%xmm3,%%xmm0								\n\t	unpcklpd	0xa0(%%rbx),%%xmm12		\n\t"\
		"movaps		%%xmm1,%%xmm3								\n\t	movaps			%%xmm14,0x120(%%r10)\n\t"\
		"movaps		%%xmm0,%%xmm2								\n\t	movaps		0xb0(%%rax),%%xmm15		\n\t"\
		"movq		%[__c12],%%rdx								\n\t	movaps		0xb0(%%rax),%%xmm13		\n\t"\
		"movaps		0xc0(%%rax),%%xmm6							\n\t	unpckhpd	0xb0(%%rbx),%%xmm15		\n\t"\
		"movaps		0xc0(%%rax),%%xmm4							\n\t	unpcklpd	0xb0(%%rbx),%%xmm13		\n\t"\
		"unpckhpd	0xc0(%%rbx),%%xmm6							\n\t	movaps			%%xmm15,0x130(%%r10)\n\t"\
		"unpcklpd	0xc0(%%rbx),%%xmm4							\n\t	movaps			%%xmm12,%%xmm14		\n\t"\
		"movaps		%%xmm6,0x160(%%rcx)							\n\t	movaps			%%xmm13,%%xmm15		\n\t"\
		"movaps		0xd0(%%rax),%%xmm7							\n\t	mulpd		    (%%r11),%%xmm12		\n\t"\
		"movaps		0xd0(%%rax),%%xmm5							\n\t	mulpd		    (%%r11),%%xmm13		\n\t"\
		"unpckhpd	0xd0(%%rbx),%%xmm7							\n\t	mulpd		0x10(%%r11),%%xmm14		\n\t"\
		"unpcklpd	0xd0(%%rbx),%%xmm5							\n\t	mulpd		0x10(%%r11),%%xmm15		\n\t"\
		"movaps		%%xmm7,0x170(%%rcx)							\n\t	addpd		%%xmm14,%%xmm13			\n\t"\
		"movaps		%%xmm4,%%xmm6								\n\t	subpd		%%xmm15,%%xmm12			\n\t"\
		"movaps		%%xmm5,%%xmm7								\n\t	addpd		%%xmm12,%%xmm8			\n\t"\
	"prefetcht1	0x80(%%r14)\n\t"\
		"mulpd		    (%%rdx),%%xmm4							\n\t	addpd		%%xmm13,%%xmm9			\n\t"\
		"mulpd		    (%%rdx),%%xmm5							\n\t	subpd		%%xmm12,%%xmm10			\n\t"\
		"mulpd		0x10(%%rdx),%%xmm6							\n\t	subpd		%%xmm13,%%xmm11			\n\t"\
		"mulpd		0x10(%%rdx),%%xmm7							\n\t	movq		%[__c14],%%r11			\n\t"\
		"addpd		%%xmm6,%%xmm5								\n\t	movaps		0xe0(%%rax),%%xmm14		\n\t"\
		"subpd		%%xmm7,%%xmm4								\n\t	movaps		0xe0(%%rax),%%xmm12		\n\t"\
		"addpd		%%xmm4,%%xmm0								\n\t	unpckhpd	0xe0(%%rbx),%%xmm14		\n\t"\
		"addpd		%%xmm5,%%xmm1								\n\t	unpcklpd	0xe0(%%rbx),%%xmm12		\n\t"\
		"subpd		%%xmm4,%%xmm2								\n\t	movaps			%%xmm14,0x160(%%r10)\n\t"\
		"subpd		%%xmm5,%%xmm3								\n\t	movaps		0xf0(%%rax),%%xmm15		\n\t"\
		"movq		%[__c8] ,%%rdx								\n\t	movaps		0xf0(%%rax),%%xmm13		\n\t"\
		"movaps		0x80(%%rax),%%xmm6							\n\t	unpckhpd	0xf0(%%rbx),%%xmm15		\n\t"\
		"movaps		0x80(%%rax),%%xmm4							\n\t	unpcklpd	0xf0(%%rbx),%%xmm13		\n\t"\
		"unpckhpd	0x80(%%rbx),%%xmm6							\n\t	movaps			%%xmm15,0x170(%%r10)\n\t"\
		"unpcklpd	0x80(%%rbx),%%xmm4							\n\t	movaps			%%xmm12,%%xmm14		\n\t"\
		"movaps			%%xmm6,0x120(%%rcx)						\n\t	movaps			%%xmm13,%%xmm15		\n\t"\
		"movaps		0x90(%%rax),%%xmm7							\n\t	mulpd		    (%%r11),%%xmm12		\n\t"\
		"movaps		0x90(%%rax),%%xmm5							\n\t	mulpd		    (%%r11),%%xmm13		\n\t"\
		"unpckhpd	0x90(%%rbx),%%xmm7							\n\t	mulpd		0x10(%%r11),%%xmm14		\n\t"\
		"unpcklpd	0x90(%%rbx),%%xmm5							\n\t	mulpd		0x10(%%r11),%%xmm15		\n\t"\
		"movaps		%%xmm7,0x130(%%rcx)							\n\t	addpd		%%xmm14,%%xmm13			\n\t"\
		"movaps		%%xmm4,%%xmm6								\n\t	subpd		%%xmm15,%%xmm12			\n\t"\
		"movaps		%%xmm5,%%xmm7								\n\t	movaps		%%xmm13,0x070(%%r10)	\n\t"\
		"mulpd		    (%%rdx),%%xmm4							\n\t	movaps		%%xmm12,0x060(%%r10)	\n\t"\
		"mulpd		    (%%rdx),%%xmm5							\n\t	movq		%[__c6] ,%%r11			\n\t"\
		"mulpd		0x10(%%rdx),%%xmm6							\n\t	movaps		0x60(%%rax),%%xmm14		\n\t"\
		"mulpd		0x10(%%rdx),%%xmm7							\n\t	movaps		0x60(%%rax),%%xmm12		\n\t"\
		"addpd		%%xmm6,%%xmm5								\n\t	unpckhpd	0x60(%%rbx),%%xmm14		\n\t"\
		"subpd		%%xmm7,%%xmm4								\n\t	unpcklpd	0x60(%%rbx),%%xmm12		\n\t"\
		"movaps		    (%%rax),%%xmm6							\n\t	movaps			%%xmm14,0x140(%%r10)\n\t"\
		"movaps		    (%%rax),%%xmm7							\n\t	movaps		0x70(%%rax),%%xmm15		\n\t"\
		"unpckhpd	    (%%rbx),%%xmm6							\n\t	movaps		0x70(%%rax),%%xmm13		\n\t"\
		"unpcklpd	    (%%rbx),%%xmm7							\n\t	unpckhpd	0x70(%%rbx),%%xmm15		\n\t"\
		"movaps		%%xmm6,0x100(%%rcx)							\n\t	unpcklpd	0x70(%%rbx),%%xmm13		\n\t"\
		"movaps		%%xmm7,     (%%rcx)							\n\t	movaps			%%xmm15,0x150(%%r10)\n\t"\
		"movaps		0x10(%%rax),%%xmm6							\n\t	movaps			%%xmm12,%%xmm14		\n\t"\
		"movaps		0x10(%%rax),%%xmm7							\n\t	movaps			%%xmm13,%%xmm15		\n\t"\
		"unpckhpd	0x10(%%rbx),%%xmm6							\n\t	mulpd		    (%%r11),%%xmm12		\n\t"\
		"unpcklpd	0x10(%%rbx),%%xmm7							\n\t	mulpd		    (%%r11),%%xmm13		\n\t"\
		"movaps		%%xmm6,0x110(%%rcx)							\n\t	mulpd		0x10(%%r11),%%xmm14		\n\t"\
		"movaps		    (%%rcx),%%xmm6							\n\t	mulpd		0x10(%%r11),%%xmm15		\n\t"\
		"subpd		%%xmm4,%%xmm6								\n\t	addpd		%%xmm14,%%xmm13			\n\t"\
		"subpd		%%xmm5,%%xmm7								\n\t	subpd		%%xmm15,%%xmm12			\n\t"\
		"addpd		%%xmm4,%%xmm4								\n\t	movaps		%%xmm13,%%xmm15			\n\t"\
		"addpd		%%xmm5,%%xmm5								\n\t	movaps		%%xmm12,%%xmm14			\n\t"\
		"addpd		%%xmm6,%%xmm4								\n\t	subpd		0x060(%%r10),%%xmm12	\n\t"\
		"addpd		%%xmm7,%%xmm5								\n\t	subpd		0x070(%%r10),%%xmm13	\n\t"\
		"/* Finish radix-4 butterfly: */						\n\t	addpd		0x060(%%r10),%%xmm14	\n\t"\
		"subpd		%%xmm0,%%xmm4								\n\t	addpd		0x070(%%r10),%%xmm15	\n\t"\
		"subpd		%%xmm1,%%xmm5								\n\t	/* Finish radix-4 butterfly: */		\n\t"\
	"prefetcht1	0xc0(%%r14)\n\t"\
		"movaps		%%xmm4,0x040(%%rcx)							\n\t	subpd		%%xmm14,%%xmm8			\n\t"\
		"movaps		%%xmm5,0x050(%%rcx)							\n\t	subpd		%%xmm15,%%xmm9			\n\t"\
		"addpd		%%xmm0,%%xmm0								\n\t	movaps		%%xmm8,0x040(%%r10)		\n\t"\
		"addpd		%%xmm1,%%xmm1								\n\t	movaps		%%xmm9,0x050(%%r10)		\n\t"\
		"addpd		%%xmm4,%%xmm0								\n\t	addpd		%%xmm14,%%xmm14			\n\t"\
		"addpd		%%xmm5,%%xmm1								\n\t	addpd		%%xmm15,%%xmm15			\n\t"\
		"movaps		%%xmm0,     (%%rcx)							\n\t	addpd		%%xmm8,%%xmm14			\n\t"\
		"movaps		%%xmm1,0x010(%%rcx)							\n\t	addpd		%%xmm9,%%xmm15			\n\t"\
		"subpd		%%xmm3,%%xmm6								\n\t	movaps		%%xmm14,     (%%r10)	\n\t"\
		"subpd		%%xmm2,%%xmm7								\n\t	movaps		%%xmm15,0x010(%%r10)	\n\t"\
		"movaps		%%xmm6,0x020(%%rcx)							\n\t	subpd		%%xmm13,%%xmm10			\n\t"\
		"movaps		%%xmm7,0x070(%%rcx)							\n\t	subpd		%%xmm12,%%xmm11			\n\t"\
		"addpd		%%xmm3,%%xmm3								\n\t	movaps		%%xmm10,0x020(%%r10)	\n\t"\
		"addpd		%%xmm2,%%xmm2								\n\t	movaps		%%xmm11,0x070(%%r10)	\n\t"\
		"addpd		%%xmm6,%%xmm3								\n\t	addpd		%%xmm13,%%xmm13			\n\t"\
		"addpd		%%xmm7,%%xmm2								\n\t	addpd		%%xmm12,%%xmm12			\n\t"\
		"movaps		%%xmm3,0x060(%%rcx)							\n\t	addpd		%%xmm10,%%xmm13			\n\t"\
		"movaps		%%xmm2,0x030(%%rcx)							\n\t	addpd		%%xmm11,%%xmm12			\n\t"\
		"														\n\t	movaps		%%xmm13,0x060(%%r10)	\n\t"\
		"														\n\t	movaps		%%xmm12,0x030(%%r10)	\n\t"\
	/****************************************************************************************************/\
	/* Each of the above 2 radix-4 blocks is unique, but the ensuing 2 blocks [operating on odd-indexed */\
	/* elements from the unpck*pd commands which were stored to temporaries] can use a common macro:    */\
	/****************************************************************************************************/\
		"/*...Block 3: */										\n\t	/*...Block 4: */					\n\t"\
		"/*	SSE2_RADIX4_DIF_4TWIDDLE(r17,r21,r19,r23, c1): */	\n\t	/* SSE2_RADIX4_DIF_4TWIDDLE(r25,r29,r27,r31, c3): */\n\t"\
		"/* Do the p0,p8 combo: */								\n\t	/* Do the p0,p8 combo: */			\n\t"\
		"movq		%[__r17],%%rax								\n\t	movq		%[__r25],%%r10			\n\t"\
		"movq		%[__c1] ,%%rbx								\n\t	movq		%[__c3] ,%%r11			\n\t"\
		"movq		%%rax   ,%%rcx								\n\t	movq		%%r10   ,%%r12			\n\t"\
		"addq		$0x20   ,%%rcx								\n\t	addq		$0x20   ,%%r12			\n\t"\
		"movaps		    (%%rax),%%xmm0							\n\t	movaps		    (%%r10),%%xmm8 		\n\t"\
		"movaps	        (%%rcx),%%xmm4							\n\t	movaps	        (%%r12),%%xmm12		\n\t"\
		"movaps		0x10(%%rax),%%xmm1							\n\t	movaps		0x10(%%r10),%%xmm9 		\n\t"\
		"movaps		0x10(%%rcx),%%xmm5							\n\t	movaps		0x10(%%r12),%%xmm13		\n\t"\
		"movaps		    (%%rbx),%%xmm6							\n\t	movaps		    (%%r11),%%xmm14		\n\t"\
		"movaps		0x10(%%rbx),%%xmm7							\n\t	movaps		0x10(%%r11),%%xmm15		\n\t"\
		"movaps		%%xmm0,%%xmm2								\n\t	movaps		%%xmm8 ,%%xmm10			\n\t"\
		"movaps		%%xmm1,%%xmm3								\n\t	movaps		%%xmm9 ,%%xmm11			\n\t"\
		"mulpd		%%xmm6,%%xmm0								\n\t	mulpd   	%%xmm14,%%xmm8 			\n\t"\
		"mulpd		%%xmm6,%%xmm1								\n\t	mulpd   	%%xmm14,%%xmm9 			\n\t"\
		"mulpd		%%xmm7,%%xmm2								\n\t	mulpd   	%%xmm15,%%xmm10			\n\t"\
		"mulpd		%%xmm7,%%xmm3								\n\t	mulpd   	%%xmm15,%%xmm11			\n\t"\
		"movaps		%%xmm4,%%xmm6								\n\t	movaps		%%xmm12,%%xmm14			\n\t"\
		"addpd		%%xmm2,%%xmm1								\n\t	addpd   	%%xmm10,%%xmm9 			\n\t"\
		"movaps		%%xmm5,%%xmm7								\n\t	movaps		%%xmm13,%%xmm15			\n\t"\
		"mulpd		0x20(%%rbx),%%xmm4							\n\t	mulpd		0x20(%%r11),%%xmm12		\n\t"\
		"subpd		%%xmm3,%%xmm0								\n\t	subpd		%%xmm11,%%xmm8 			\n\t"\
		"mulpd		0x20(%%rbx),%%xmm5							\n\t	mulpd		0x20(%%r11),%%xmm13		\n\t"\
		"mulpd		0x30(%%rbx),%%xmm6							\n\t	mulpd		0x30(%%r11),%%xmm14		\n\t"\
	"prefetcht1	0x100(%%r14)\n\t"\
		"movaps		%%xmm0,%%xmm2								\n\t	movaps		%%xmm8 ,%%xmm10			\n\t"\
		"mulpd		0x30(%%rbx),%%xmm7							\n\t	mulpd		0x30(%%r11),%%xmm15		\n\t"\
		"addpd		%%xmm6,%%xmm5								\n\t	addpd		%%xmm14,%%xmm13			\n\t"\
		"movaps		%%xmm1,%%xmm3								\n\t	movaps		%%xmm9 ,%%xmm11			\n\t"\
		"subpd		%%xmm7,%%xmm4								\n\t	subpd		%%xmm15,%%xmm12			\n\t"\
		"addq		$0x40,%%rcx									\n\t	addq		$0x40	,%%r12			\n\t"\
		"addq		$0x60,%%rbx									\n\t	addq		$0x60	,%%r11			\n\t"\
		"movaps			(%%rcx),%%xmm6							\n\t	movaps		    (%%r12),%%xmm14		\n\t"\
		"movaps		0x10(%%rcx),%%xmm7							\n\t	movaps		0x10(%%r12),%%xmm15		\n\t"\
		"addpd		%%xmm4,%%xmm0								\n\t	addpd		%%xmm12,%%xmm8 			\n\t"\
		"addpd		%%xmm5,%%xmm1								\n\t	addpd		%%xmm13,%%xmm9 			\n\t"\
		"subpd		%%xmm4,%%xmm2								\n\t	subpd		%%xmm12,%%xmm10			\n\t"\
		"subpd		%%xmm5,%%xmm3								\n\t	subpd		%%xmm13,%%xmm11			\n\t"\
		"/* Do the p4,12 combo: */								\n\t	/* Do the p4,12 combo: */			\n\t"\
		"movaps		%%xmm6,%%xmm4								\n\t	movaps		%%xmm14,%%xmm12			\n\t"\
		"movaps		%%xmm7,%%xmm5								\n\t	movaps		%%xmm15,%%xmm13			\n\t"\
		"mulpd		    (%%rbx),%%xmm4							\n\t	mulpd		    (%%r11),%%xmm12		\n\t"\
		"mulpd		    (%%rbx),%%xmm5							\n\t	mulpd		    (%%r11),%%xmm13		\n\t"\
		"mulpd		0x10(%%rbx),%%xmm6							\n\t	mulpd		0x10(%%r11),%%xmm14		\n\t"\
		"mulpd		0x10(%%rbx),%%xmm7							\n\t	mulpd		0x10(%%r11),%%xmm15		\n\t"\
		"movq		%%rax,%%rdx									\n\t	movq		%%r10,%%r13				\n\t"\
		"addpd		%%xmm6,%%xmm5								\n\t	addpd		%%xmm14,%%xmm13			\n\t"\
		"subpd		%%xmm7,%%xmm4								\n\t	subpd		%%xmm15,%%xmm12			\n\t"\
		"movaps		%%xmm5,0x010(%%rdx)							\n\t	movaps		%%xmm13,0x010(%%r13)	\n\t"\
		"movaps		%%xmm4,     (%%rdx)							\n\t	movaps		%%xmm12,     (%%r13)	\n\t"\
		"addq		$0x40 ,%%rax								\n\t	addq		$0x40 ,%%r10			\n\t"\
		"subq		$0x20 ,%%rbx								\n\t	subq		$0x20 ,%%r11			\n\t"\
		"movaps		    (%%rax),%%xmm4							\n\t	movaps		    (%%r10),%%xmm12		\n\t"\
		"movaps		0x10(%%rax),%%xmm5							\n\t	movaps		0x10(%%r10),%%xmm13		\n\t"\
		"movaps			%%xmm4,%%xmm6							\n\t	movaps			%%xmm12,%%xmm14		\n\t"\
		"movaps			%%xmm5,%%xmm7							\n\t	movaps			%%xmm13,%%xmm15		\n\t"\
		"mulpd		    (%%rbx),%%xmm4							\n\t	mulpd		    (%%r11),%%xmm12		\n\t"\
		"mulpd		    (%%rbx),%%xmm5							\n\t	mulpd		    (%%r11),%%xmm13		\n\t"\
		"mulpd		0x10(%%rbx),%%xmm6							\n\t	mulpd		0x10(%%r11),%%xmm14		\n\t"\
		"mulpd		0x10(%%rbx),%%xmm7							\n\t	mulpd		0x10(%%r11),%%xmm15		\n\t"\
		"addpd		%%xmm6,%%xmm5								\n\t	addpd		%%xmm14,%%xmm13			\n\t"\
		"subpd		%%xmm7,%%xmm4								\n\t	subpd		%%xmm15,%%xmm12			\n\t"\
		"movaps		%%xmm5,%%xmm7								\n\t	movaps		%%xmm13,%%xmm15			\n\t"\
		"movaps		%%xmm4,%%xmm6								\n\t	movaps		%%xmm12,%%xmm14			\n\t"\
		"subpd		     (%%rdx),%%xmm4							\n\t	subpd		     (%%r13),%%xmm12	\n\t"\
		"subpd		0x010(%%rdx),%%xmm5							\n\t	subpd		0x010(%%r13),%%xmm13	\n\t"\
		"addpd		     (%%rdx),%%xmm6							\n\t	addpd		     (%%r13),%%xmm14	\n\t"\
		"addpd		0x010(%%rdx),%%xmm7							\n\t	addpd		0x010(%%r13),%%xmm15	\n\t"\
		"/* Finish radix-4 butterfly: */						\n\t	/* Finish radix-4 butterfly: */		\n\t"\
		"subpd		%%xmm6,%%xmm0								\n\t	subpd		%%xmm14,%%xmm8 			\n\t"\
		"subpd		%%xmm5,%%xmm2								\n\t	subpd		%%xmm13,%%xmm10			\n\t"\
		"subpd		%%xmm7,%%xmm1								\n\t	subpd		%%xmm15,%%xmm9 			\n\t"\
		"subpd		%%xmm4,%%xmm3								\n\t	subpd		%%xmm12,%%xmm11			\n\t"\
	"prefetcht1	0x140(%%r14)\n\t"\
		"movaps		%%xmm0,0x040(%%rdx)							\n\t	movaps		%%xmm8 ,0x040(%%r13)	\n\t"\
		"movaps		%%xmm2,0x020(%%rdx)							\n\t	movaps		%%xmm10,0x020(%%r13)	\n\t"\
		"movaps		%%xmm1,0x050(%%rdx)							\n\t	movaps		%%xmm9 ,0x050(%%r13)	\n\t"\
		"movaps		%%xmm3,0x070(%%rdx)							\n\t	movaps		%%xmm11,0x070(%%r13)	\n\t"\
		"addpd		%%xmm6,%%xmm6								\n\t	addpd		%%xmm14,%%xmm14			\n\t"\
		"addpd		%%xmm5,%%xmm5								\n\t	addpd		%%xmm13,%%xmm13			\n\t"\
		"addpd		%%xmm7,%%xmm7								\n\t	addpd		%%xmm15,%%xmm15			\n\t"\
		"addpd		%%xmm4,%%xmm4								\n\t	addpd		%%xmm12,%%xmm12			\n\t"\
		"addpd		%%xmm0,%%xmm6								\n\t	addpd		%%xmm8 ,%%xmm14			\n\t"\
		"addpd		%%xmm2,%%xmm5								\n\t	addpd		%%xmm10,%%xmm13			\n\t"\
		"addpd		%%xmm1,%%xmm7								\n\t	addpd		%%xmm9 ,%%xmm15			\n\t"\
		"addpd		%%xmm3,%%xmm4								\n\t	addpd		%%xmm11,%%xmm12			\n\t"\
		"movaps		%%xmm6,     (%%rdx)							\n\t	movaps		%%xmm14,     (%%r13)	\n\t"\
		"movaps		%%xmm5,0x060(%%rdx)							\n\t	movaps		%%xmm13,0x060(%%r13)	\n\t"\
		"movaps		%%xmm7,0x010(%%rdx)							\n\t	movaps		%%xmm15,0x010(%%r13)	\n\t"\
		"movaps		%%xmm4,0x030(%%rdx)							\n\t	movaps		%%xmm12,0x030(%%r13)	\n\t"\
	/**************************************************************************************/\
	/*...and now do four more radix-4 transforms, including the internal twiddle factors: */\
	/**************************************************************************************/\
		"/*...Block 1: t1,9,17,25 */		\n\t		/*...Block 3: t5,13,21,29 */		\n\t"\
		"movq		%[__r1] ,%%rax						\n\t"\
		"movq		%[__r9] ,%%rbx						\n\t"\
		"movq		%[__r17],%%rcx						\n\t"\
		"movq		%[__r25],%%rdx			\n\t		movq		%[__isrt2],%%rsi		\n\t"\
		"movaps		     (%%rax),%%xmm0		\n\t		movaps		0x040(%%rax),%%xmm8		\n\t"\
		"movaps		0x010(%%rax),%%xmm1		\n\t		movaps		0x050(%%rax),%%xmm9		\n\t"\
		"movaps		     (%%rbx),%%xmm2		\n\t		movaps		0x0c0(%%rax),%%xmm10	\n\t"\
		"movaps		0x010(%%rbx),%%xmm3		\n\t		movaps		0x0d0(%%rax),%%xmm11	\n\t"\
		"subpd		     (%%rbx),%%xmm0		\n\t		subpd		0x0d0(%%rax),%%xmm8		\n\t"\
		"subpd		0x010(%%rbx),%%xmm1		\n\t		subpd		0x0c0(%%rax),%%xmm9		\n\t"\
		"addpd		     (%%rax),%%xmm2		\n\t		addpd		0x050(%%rax),%%xmm10	\n\t"\
		"addpd		0x010(%%rax),%%xmm3		\n\t		addpd		0x040(%%rax),%%xmm11	\n\t"\
		"movaps		     (%%rcx),%%xmm4		\n\t		movaps		0x140(%%rax),%%xmm12	\n\t"\
		"movaps		0x010(%%rcx),%%xmm5		\n\t		movaps		0x150(%%rax),%%xmm13	\n\t"\
		"movaps		     (%%rdx),%%xmm6		\n\t		movaps		0x1c0(%%rax),%%xmm14	\n\t"\
		"movaps		0x010(%%rdx),%%xmm7		\n\t		movaps		0x1d0(%%rax),%%xmm15	\n\t"\
		"subpd		     (%%rdx),%%xmm4		\n\t		subpd		0x150(%%rax),%%xmm12	\n\t"\
		"subpd		0x010(%%rdx),%%xmm5		\n\t		addpd		0x140(%%rax),%%xmm13	\n\t"\
		"addpd		     (%%rcx),%%xmm6		\n\t		mulpd		(%%rsi),%%xmm12			\n\t"\
		"addpd		0x010(%%rcx),%%xmm7		\n\t		mulpd		(%%rsi),%%xmm13			\n\t"\
		"subpd		%%xmm6,%%xmm2			\n\t		addpd		0x1d0(%%rax),%%xmm14	\n\t"\
		"subpd		%%xmm7,%%xmm3			\n\t		subpd		0x1c0(%%rax),%%xmm15	\n\t"\
		"addpd		%%xmm6,%%xmm6			\n\t		mulpd		(%%rsi),%%xmm14			\n\t"\
		"addpd		%%xmm7,%%xmm7			\n\t		mulpd		(%%rsi),%%xmm15			\n\t"\
		"movaps		%%xmm2,    (%%rcx)		\n\t		subpd		%%xmm14,%%xmm12			\n\t"\
		"movaps		%%xmm3,0x10(%%rcx)		\n\t		subpd		%%xmm15,%%xmm13			\n\t"\
		"addpd		%%xmm2,%%xmm6			\n\t		addpd		%%xmm14,%%xmm14			\n\t"\
		"addpd		%%xmm3,%%xmm7			\n\t		addpd		%%xmm15,%%xmm15			\n\t"\
		"movaps		%%xmm6,    (%%rax)		\n\t		addpd		%%xmm12,%%xmm14			\n\t"\
		"movaps		%%xmm7,0x10(%%rax)		\n\t		addpd		%%xmm13,%%xmm15			\n\t"\
	"prefetcht1	0x180(%%r14)\n\t"\
		"subpd		%%xmm5,%%xmm0			\n\t		subpd		%%xmm12,%%xmm8			\n\t"\
		"subpd		%%xmm4,%%xmm1			\n\t		subpd		%%xmm13,%%xmm10			\n\t"\
		"addpd		%%xmm5,%%xmm5			\n\t		addpd		%%xmm12,%%xmm12			\n\t"\
		"addpd		%%xmm4,%%xmm4			\n\t		addpd		%%xmm13,%%xmm13			\n\t"\
		"movaps		%%xmm0,    (%%rbx)		\n\t		movaps		%%xmm8,0x40(%%rcx)		\n\t"\
		"movaps		%%xmm1,0x10(%%rdx)		\n\t		movaps		%%xmm10,0x50(%%rcx)		\n\t"\
		"addpd		%%xmm0,%%xmm5			\n\t		addpd		%%xmm8,%%xmm12			\n\t"\
		"addpd		%%xmm1,%%xmm4			\n\t		addpd		%%xmm10,%%xmm13			\n\t"\
		"movaps		%%xmm5,    (%%rdx)		\n\t		movaps		%%xmm12,0x40(%%rax)		\n\t"\
		"movaps		%%xmm4,0x10(%%rbx)		\n\t		movaps		%%xmm13,0x50(%%rax)		\n\t"\
		"												subpd		%%xmm15,%%xmm11			\n\t"\
		"												subpd		%%xmm14,%%xmm9			\n\t"\
		"												addpd		%%xmm15,%%xmm15			\n\t"\
		"												addpd		%%xmm14,%%xmm14			\n\t"\
		"												movaps		%%xmm11,0x40(%%rbx)		\n\t"\
		"												movaps		%%xmm9,0x50(%%rdx)		\n\t"\
		"												addpd		%%xmm11,%%xmm15			\n\t"\
		"												addpd		%%xmm9,%%xmm14			\n\t"\
		"												movaps		%%xmm15,0x40(%%rdx)		\n\t"\
		"												movaps		%%xmm14,0x50(%%rbx)		\n\t"\
		"/*...Block 2: t3,11,19,27 */		\n\t		/*...Block 4: t7,15,23,31 */		\n\t"\
		"addq		$0x20,%%rax							\n\t"\
		"addq		$0x20,%%rbx							\n\t"\
		"addq		$0x20,%%rcx							\n\t"\
		"addq		$0x20,%%rdx							\n\t"\
		"movq		%[__cc0],%%rdi						\n\t"\
		"movaps		0x100(%%rax),%%xmm4		\n\t		movaps		0x140(%%rax),%%xmm12		\n\t"\
		"movaps		0x180(%%rax),%%xmm6		\n\t		movaps		0x1c0(%%rax),%%xmm14		\n\t"\
		"movaps		0x110(%%rax),%%xmm5		\n\t		movaps		0x150(%%rax),%%xmm13		\n\t"\
		"movaps		0x190(%%rax),%%xmm7		\n\t		movaps		0x1d0(%%rax),%%xmm15		\n\t"\
		"movaps		0x100(%%rax),%%xmm0		\n\t		movaps		0x140(%%rax),%%xmm8		\n\t"\
	/*	"movaps		0x180(%%rax),%%xmm2		\n\t		movaps		0x1c0(%%rax),%%xmm10		\n\t"*/\
		"movaps		0x110(%%rax),%%xmm1		\n\t		movaps		0x150(%%rax),%%xmm9		\n\t"\
		"movaps		0x190(%%rax),%%xmm3		\n\t		movaps		0x1d0(%%rax),%%xmm11		\n\t"\
	"movaps	    (%%rdi),%%xmm2			\n\t	movaps	0x10(%%rdi),%%xmm10			\n\t"\
		"mulpd		     %%xmm2 ,%%xmm4		\n\t		mulpd		     %%xmm10,%%xmm12		\n\t"\
		"mulpd		     %%xmm10,%%xmm6		\n\t		mulpd		     %%xmm2 ,%%xmm14		\n\t"\
		"mulpd		     %%xmm10,%%xmm1		\n\t		mulpd		     %%xmm2 ,%%xmm9			\n\t"\
		"mulpd		     %%xmm2 ,%%xmm3		\n\t		mulpd		     %%xmm10,%%xmm11		\n\t"\
		"mulpd		     %%xmm2 ,%%xmm5		\n\t		mulpd		     %%xmm10,%%xmm13		\n\t"\
		"mulpd		     %%xmm10,%%xmm7		\n\t		mulpd		     %%xmm2 ,%%xmm15		\n\t"\
		"mulpd		     %%xmm10,%%xmm0		\n\t		mulpd		     %%xmm2 ,%%xmm8			\n\t"\
		"mulpd		0x180(%%rax),%%xmm2 	\n\t		mulpd		0x1c0(%%rax),%%xmm10		\n\t"\
		"subpd		%%xmm1,%%xmm4			\n\t		subpd		%%xmm9 ,%%xmm12			\n\t"\
		"subpd		%%xmm3,%%xmm6			\n\t		subpd		%%xmm11,%%xmm14			\n\t"\
		"addpd		%%xmm0,%%xmm5			\n\t		addpd		%%xmm8 ,%%xmm13			\n\t"\
		"addpd		%%xmm2,%%xmm7			\n\t		addpd		%%xmm10,%%xmm15			\n\t"\
		"subpd		%%xmm6,%%xmm4			\n\t		subpd		%%xmm14,%%xmm12			\n\t"\
		"subpd		%%xmm7,%%xmm5			\n\t		subpd		%%xmm15,%%xmm13			\n\t"\
		"addpd		%%xmm6,%%xmm6			\n\t		addpd		%%xmm14,%%xmm14			\n\t"\
		"addpd		%%xmm7,%%xmm7			\n\t		addpd		%%xmm15,%%xmm15			\n\t"\
		"addpd		%%xmm4,%%xmm6			\n\t		addpd		%%xmm12,%%xmm14			\n\t"\
		"addpd		%%xmm5,%%xmm7			\n\t		addpd		%%xmm13,%%xmm15			\n\t"\
	"prefetcht1	0x1c0(%%r14)\n\t"\
		"movaps		0x080(%%rax),%%xmm2		\n\t		movaps		0x0c0(%%rax),%%xmm10		\n\t"\
		"movaps		0x090(%%rax),%%xmm3		\n\t		movaps		0x0d0(%%rax),%%xmm11		\n\t"\
		"subpd		0x090(%%rax),%%xmm2		\n\t		addpd		0x0d0(%%rax),%%xmm10		\n\t"\
		"addpd		0x080(%%rax),%%xmm3		\n\t		subpd		0x0c0(%%rax),%%xmm11		\n\t"\
		"mulpd		(%%rsi),%%xmm2			\n\t		mulpd		(%%rsi),%%xmm10			\n\t"\
		"mulpd		(%%rsi),%%xmm3			\n\t		mulpd		(%%rsi),%%xmm11			\n\t"\
		"movaps		     (%%rax),%%xmm0		\n\t		movaps		0x040(%%rax),%%xmm8		\n\t"\
		"movaps		0x010(%%rax),%%xmm1		\n\t		movaps		0x050(%%rax),%%xmm9		\n\t"\
		"subpd		%%xmm2,%%xmm0			\n\t		subpd		%%xmm10,%%xmm8			\n\t"\
		"subpd		%%xmm3,%%xmm1			\n\t		subpd		%%xmm11,%%xmm9			\n\t"\
		"addpd		     (%%rax),%%xmm2		\n\t		addpd		0x040(%%rax),%%xmm10		\n\t"\
		"addpd		0x010(%%rax),%%xmm3		\n\t		addpd		0x050(%%rax),%%xmm11		\n\t"\
		"subpd		%%xmm6,%%xmm2			\n\t		subpd		%%xmm12,%%xmm8			\n\t"\
		"subpd		%%xmm7,%%xmm3			\n\t		subpd		%%xmm13,%%xmm9			\n\t"\
		"addpd		%%xmm6,%%xmm6			\n\t		addpd		%%xmm12,%%xmm12			\n\t"\
		"addpd		%%xmm7,%%xmm7			\n\t		addpd		%%xmm13,%%xmm13			\n\t"\
		"movaps		%%xmm2,    (%%rcx)		\n\t		movaps		%%xmm8,0x40(%%rcx)		\n\t"\
		"movaps		%%xmm3,0x10(%%rcx)		\n\t		movaps		%%xmm9,0x50(%%rcx)		\n\t"\
		"addpd		%%xmm2,%%xmm6			\n\t		addpd		%%xmm8,%%xmm12			\n\t"\
		"addpd		%%xmm3,%%xmm7			\n\t		addpd		%%xmm9,%%xmm13			\n\t"\
		"movaps		%%xmm6,    (%%rax)		\n\t		movaps		%%xmm12,0x40(%%rax)		\n\t"\
		"movaps		%%xmm7,0x10(%%rax)		\n\t		movaps		%%xmm13,0x50(%%rax)		\n\t"\
		"subpd		%%xmm5,%%xmm0			\n\t		subpd		%%xmm15,%%xmm10			\n\t"\
		"subpd		%%xmm4,%%xmm1			\n\t		subpd		%%xmm14,%%xmm11			\n\t"\
		"addpd		%%xmm5,%%xmm5			\n\t		addpd		%%xmm15,%%xmm15			\n\t"\
		"addpd		%%xmm4,%%xmm4			\n\t		addpd		%%xmm14,%%xmm14			\n\t"\
		"movaps		%%xmm0,    (%%rbx)		\n\t		movaps		%%xmm10,0x40(%%rbx)		\n\t"\
		"movaps		%%xmm1,0x10(%%rdx)		\n\t		movaps		%%xmm11,0x50(%%rdx)		\n\t"\
		"addpd		%%xmm0,		%%xmm5		\n\t		addpd		%%xmm10,%%xmm15			\n\t"\
		"addpd		%%xmm1,		%%xmm4		\n\t		addpd		%%xmm11,%%xmm14			\n\t"\
		"movaps		%%xmm5,    (%%rdx)		\n\t		movaps		%%xmm15,0x40(%%rdx)		\n\t"\
		"movaps		%%xmm4,0x10(%%rbx)		\n\t		movaps		%%xmm14,0x50(%%rbx)		\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__add1] "m" (Xadd1)\
		 ,[__r1] "m" (Xr1)\
		 ,[__r9] "m" (Xr9)\
		 ,[__r17] "m" (Xr17)\
		 ,[__r25] "m" (Xr25)\
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__cc0] "m" (Xcc0)\
		 ,[__c1] "m" (Xc1)\
		 ,[__c2] "m" (Xc2)\
		 ,[__c3] "m" (Xc3)\
		 ,[__c4] "m" (Xc4)\
		 ,[__c5] "m" (Xc5)\
		 ,[__c6] "m" (Xc6)\
		 ,[__c7] "m" (Xc7)\
		 ,[__c8] "m" (Xc8)\
		 ,[__c9] "m" (Xc9)\
		 ,[__c10] "m" (Xc10)\
		 ,[__c11] "m" (Xc11)\
		 ,[__c12] "m" (Xc12)\
		 ,[__c13] "m" (Xc13)\
		 ,[__c14] "m" (Xc14)\
		 ,[__c15] "m" (Xc15)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r10","r11","r12","r13","r14","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"		/* Clobbered registers */\
	);\
	}

	// Since the add0/1 block addresses advance in opposite directions, our prefetch scheme here is like so:
	// In fwd DIF prefetch  ahead 8 cache lines' worth of data w.r.to add0;
	// In inv DIT prefetch behind 8 cache lines' worth of data w.r.to add1;
	#define SSE2_RADIX16_WRAPPER_DIT(Xadd0,Xadd1,Xr1,Xr9,Xr17,Xr25,Xisrt2,Xcc0,Xc1,Xc2,Xc3,Xc4,Xc5,Xc6,Xc7,Xc8,Xc9,Xc10,Xc11,Xc12,Xc13,Xc14,Xc15)\
	{\
	__asm__ volatile (\
	"movq	%[__add0],%%r14						\n\t"\
	"prefetcht1	-0x40(%%r14)\n\t"\
		/*...Block 1: */\
		"movq	%[__r1],%%rax					\n\t"\
		"movq	%%rax,%%rbx						\n\t"\
		"movq	%%rax,%%rcx						\n\t"\
		"movq	%%rax,%%rdx						\n\t"\
		"addq	$0x100,%%rbx					\n\t"\
		"addq	$0x080,%%rcx					\n\t"\
		"addq	$0x180,%%rdx					\n\t"		/*...Block 3: */\
		/* SSE2_RADIX4_DIT_IN_PLACE(): */					/* SSE2_RADIX4_DIT_IN_PLACE(): */\
		"movaps	    (%%rax),%%xmm0				\n\t		movaps	0x20(%%rax),%%xmm8 				\n\t"\
		"movaps	0x10(%%rax),%%xmm1				\n\t		movaps	0x30(%%rax),%%xmm9 				\n\t"\
		"movaps	    (%%rbx),%%xmm4				\n\t		movaps	0x20(%%rbx),%%xmm12				\n\t"\
		"movaps	0x10(%%rbx),%%xmm5				\n\t		movaps	0x30(%%rbx),%%xmm13				\n\t"\
		"movaps	     %%xmm0,%%xmm2				\n\t		movaps	    %%xmm8 ,%%xmm10				\n\t"\
		"movaps	     %%xmm1,%%xmm3				\n\t		movaps	    %%xmm9 ,%%xmm11				\n\t"\
		"addpd	     %%xmm4,%%xmm0				\n\t		addpd	    %%xmm12,%%xmm8 				\n\t"\
		"addpd	     %%xmm5,%%xmm1				\n\t		addpd	    %%xmm13,%%xmm9 				\n\t"\
		"subpd	     %%xmm4,%%xmm2				\n\t		subpd	    %%xmm12,%%xmm10				\n\t"\
		"subpd	     %%xmm5,%%xmm3				\n\t		subpd	    %%xmm13,%%xmm11				\n\t"\
		"movaps	    (%%rcx),%%xmm4				\n\t		movaps	0x20(%%rcx),%%xmm12				\n\t"\
		"movaps	0x10(%%rcx),%%xmm5				\n\t		movaps	0x30(%%rcx),%%xmm13				\n\t"\
		"movaps	     %%xmm4,%%xmm6				\n\t		movaps	    %%xmm12,%%xmm14				\n\t"\
		"movaps	     %%xmm5,%%xmm7				\n\t		movaps	    %%xmm13,%%xmm15				\n\t"\
		"addpd	    (%%rdx),%%xmm4				\n\t		addpd	0x20(%%rdx),%%xmm12				\n\t"\
		"addpd	0x10(%%rdx),%%xmm5				\n\t		addpd	0x30(%%rdx),%%xmm13				\n\t"\
		"subpd	    (%%rdx),%%xmm6				\n\t		subpd	0x20(%%rdx),%%xmm14				\n\t"\
		"subpd	0x10(%%rdx),%%xmm7				\n\t		subpd	0x30(%%rdx),%%xmm15				\n\t"\
		"subpd	%%xmm4,%%xmm0					\n\t		subpd	%%xmm12,%%xmm8 					\n\t"\
		"subpd	%%xmm5,%%xmm1					\n\t		subpd	%%xmm13,%%xmm9 					\n\t"\
		"movaps	%%xmm0,    (%%rbx)				\n\t		movaps	%%xmm8 ,0x20(%%rbx)				\n\t"\
		"movaps	%%xmm1,0x10(%%rbx)				\n\t		movaps	%%xmm9 ,0x30(%%rbx)				\n\t"\
		"addpd	%%xmm4,%%xmm4					\n\t		addpd	%%xmm12,%%xmm12					\n\t"\
		"addpd	%%xmm5,%%xmm5					\n\t		addpd	%%xmm13,%%xmm13					\n\t"\
		"addpd	%%xmm0,%%xmm4					\n\t		addpd	%%xmm8 ,%%xmm12					\n\t"\
		"addpd	%%xmm1,%%xmm5					\n\t		addpd	%%xmm9 ,%%xmm13					\n\t"\
		"movaps	%%xmm4,    (%%rax)				\n\t		movaps	%%xmm12,0x20(%%rax)				\n\t"\
		"movaps	%%xmm5,0x10(%%rax)				\n\t		movaps	%%xmm13,0x30(%%rax)				\n\t"\
		"subpd	%%xmm7,%%xmm2					\n\t		subpd	%%xmm15,%%xmm10					\n\t"\
		"subpd	%%xmm6,%%xmm3					\n\t		subpd	%%xmm14,%%xmm11					\n\t"\
		"movaps	%%xmm2,    (%%rdx)				\n\t		movaps	%%xmm10,0x20(%%rdx)				\n\t"\
		"movaps	%%xmm3,0x10(%%rcx)				\n\t		movaps	%%xmm11,0x30(%%rcx)				\n\t"\
		"addpd	%%xmm7,%%xmm7					\n\t		addpd	%%xmm15,%%xmm15					\n\t"\
		"addpd	%%xmm6,%%xmm6					\n\t		addpd	%%xmm14,%%xmm14					\n\t"\
		"addpd	%%xmm2,%%xmm7					\n\t		addpd	%%xmm10,%%xmm15					\n\t"\
		"addpd	%%xmm3,%%xmm6					\n\t		addpd	%%xmm11,%%xmm14					\n\t"\
		"movaps	%%xmm7,    (%%rcx)				\n\t		movaps	%%xmm15,0x20(%%rcx)				\n\t"\
		"movaps	%%xmm6,0x10(%%rdx)				\n\t		movaps	%%xmm14,0x30(%%rdx)				\n\t"\
		"/*...Block 2: */						\n\t		/*...Block 4: */						\n\t"\
	"prefetcht1	-0x80(%%r14)\n\t"\
		"addq	$0x040,%%rax					\n\t"\
		"addq	$0x040,%%rbx					\n\t"\
		"addq	$0x040,%%rcx					\n\t"\
		"addq	$0x040,%%rdx					\n\t"\
		"/* SSE2_RADIX4_DIT_IN_PLACE(): */		\n\t		/* SSE2_RADIX4_DIT_IN_PLACE(): */		\n\t"\
		"movaps	    (%%rax),%%xmm0				\n\t		movaps	0x20(%%rax),%%xmm8 				\n\t"\
		"movaps	0x10(%%rax),%%xmm1				\n\t		movaps	0x30(%%rax),%%xmm9 				\n\t"\
		"movaps	    (%%rbx),%%xmm4				\n\t		movaps	0x20(%%rbx),%%xmm12				\n\t"\
		"movaps	0x10(%%rbx),%%xmm5				\n\t		movaps	0x30(%%rbx),%%xmm13				\n\t"\
		"movaps	     %%xmm0,%%xmm2				\n\t		movaps	    %%xmm8 ,%%xmm10				\n\t"\
		"movaps	     %%xmm1,%%xmm3				\n\t		movaps	    %%xmm9 ,%%xmm11				\n\t"\
		"addpd	     %%xmm4,%%xmm0				\n\t		addpd	    %%xmm12,%%xmm8 				\n\t"\
		"addpd	     %%xmm5,%%xmm1				\n\t		addpd	    %%xmm13,%%xmm9 				\n\t"\
		"subpd	     %%xmm4,%%xmm2				\n\t		subpd	    %%xmm12,%%xmm10				\n\t"\
		"subpd	     %%xmm5,%%xmm3				\n\t		subpd	    %%xmm13,%%xmm11				\n\t"\
		"movaps	    (%%rcx),%%xmm4				\n\t		movaps	0x20(%%rcx),%%xmm12				\n\t"\
		"movaps	0x10(%%rcx),%%xmm5				\n\t		movaps	0x30(%%rcx),%%xmm13				\n\t"\
		"movaps	     %%xmm4,%%xmm6				\n\t		movaps	    %%xmm12,%%xmm14				\n\t"\
		"movaps	     %%xmm5,%%xmm7				\n\t		movaps	    %%xmm13,%%xmm15				\n\t"\
		"addpd	    (%%rdx),%%xmm4				\n\t		addpd	0x20(%%rdx),%%xmm12				\n\t"\
		"addpd	0x10(%%rdx),%%xmm5				\n\t		addpd	0x30(%%rdx),%%xmm13				\n\t"\
		"subpd	    (%%rdx),%%xmm6				\n\t		subpd	0x20(%%rdx),%%xmm14				\n\t"\
		"subpd	0x10(%%rdx),%%xmm7				\n\t		subpd	0x30(%%rdx),%%xmm15				\n\t"\
		"subpd	%%xmm4,%%xmm0					\n\t		subpd	%%xmm12,%%xmm8 					\n\t"\
		"subpd	%%xmm5,%%xmm1					\n\t		subpd	%%xmm13,%%xmm9 					\n\t"\
		"movaps	%%xmm0,    (%%rbx)				\n\t		movaps	%%xmm8 ,0x20(%%rbx)				\n\t"\
		"movaps	%%xmm1,0x10(%%rbx)				\n\t		movaps	%%xmm9 ,0x30(%%rbx)				\n\t"\
		"addpd	%%xmm4,%%xmm4					\n\t		addpd	%%xmm12,%%xmm12					\n\t"\
		"addpd	%%xmm5,%%xmm5					\n\t		addpd	%%xmm13,%%xmm13					\n\t"\
		"addpd	%%xmm0,%%xmm4					\n\t		addpd	%%xmm8 ,%%xmm12					\n\t"\
		"addpd	%%xmm1,%%xmm5					\n\t		addpd	%%xmm9 ,%%xmm13					\n\t"\
		"movaps	%%xmm4,    (%%rax)				\n\t		movaps	%%xmm12,0x20(%%rax)				\n\t"\
		"movaps	%%xmm5,0x10(%%rax)				\n\t		movaps	%%xmm13,0x30(%%rax)				\n\t"\
		"subpd	%%xmm7,%%xmm2					\n\t		subpd	%%xmm15,%%xmm10					\n\t"\
		"subpd	%%xmm6,%%xmm3					\n\t		subpd	%%xmm14,%%xmm11					\n\t"\
		"movaps	%%xmm2,    (%%rdx)				\n\t		movaps	%%xmm10,0x20(%%rdx)				\n\t"\
		"movaps	%%xmm3,0x10(%%rcx)				\n\t		movaps	%%xmm11,0x30(%%rcx)				\n\t"\
		"addpd	%%xmm7,%%xmm7					\n\t		addpd	%%xmm15,%%xmm15					\n\t"\
		"addpd	%%xmm6,%%xmm6					\n\t		addpd	%%xmm14,%%xmm14					\n\t"\
		"addpd	%%xmm2,%%xmm7					\n\t		addpd	%%xmm10,%%xmm15					\n\t"\
		"addpd	%%xmm3,%%xmm6					\n\t		addpd	%%xmm11,%%xmm14					\n\t"\
		"movaps	%%xmm7,    (%%rcx)				\n\t		movaps	%%xmm15,0x20(%%rcx)				\n\t"\
		"movaps	%%xmm6,0x10(%%rdx)				\n\t		movaps	%%xmm14,0x30(%%rdx)				\n\t"\
	/***************************************************************************************************/\
	/*..and now do four more radix-4 transforms, including the internal and external twiddle factors.  */\
	/*  Write even-index 16-byte output pairs to a(j1), odd-index to a(j2), unpack same as on inputs.  */\
	/*  We do the last 2 radix-4 blocks first, to make the unpack-interleaving of outputs smoother.    */\
	/***************************************************************************************************/\
		/* Main-array addresses still in add0,1, no need to re-init: */\
		/*...Block 3: t3,11,19,27 -> r9,13,11,15: */	/*...Block 4: t7,15,23,31 -> r25,29,27,31: */\
		"movq		%[__r9],%%rax				\n\t	"/*movq		%[__r25],%%rax	Instead incr rcol rax offsets by +0x100 */\
		"movq		%[__isrt2],%%rbx			\n\t"\
		"movq		%[__cc0],%%rcx				\n\t"\
	"prefetcht1	-0xc0(%%r14)\n\t"\
		"movaps		0x020(%%rax),%%xmm4			\n\t		movaps		0x120(%%rax),%%xmm12		\n\t"\
		"movaps		0x060(%%rax),%%xmm0			\n\t		movaps		0x160(%%rax),%%xmm8 		\n\t"\
		"movaps		0x030(%%rax),%%xmm5			\n\t		movaps		0x130(%%rax),%%xmm13		\n\t"\
		"movaps		0x070(%%rax),%%xmm1			\n\t		movaps		0x170(%%rax),%%xmm9 		\n\t"\
		"movaps		0x020(%%rax),%%xmm6			\n\t		movaps		0x120(%%rax),%%xmm14		\n\t"\
		"movaps		0x060(%%rax),%%xmm2			\n\t		movaps		0x160(%%rax),%%xmm10		\n\t"\
		"movaps		0x030(%%rax),%%xmm7			\n\t		movaps		0x130(%%rax),%%xmm15		\n\t"\
	/*	"movaps		0x070(%%rax),%%xmm3			\n\t		movaps		0x170(%%rax),%%xmm11		\n\t"*/\
	"movaps	    (%%rcx),%%xmm3			\n\t	movaps	0x10(%%rcx),%%xmm11			\n\t"\
		"mulpd		    %%xmm3 ,%%xmm4			\n\t		mulpd		    %%xmm11,%%xmm12			\n\t"\
		"mulpd		    %%xmm11,%%xmm0			\n\t		mulpd		    %%xmm3 ,%%xmm8 			\n\t"\
		"mulpd		    %%xmm3 ,%%xmm5			\n\t		mulpd		    %%xmm11,%%xmm13			\n\t"\
		"mulpd		    %%xmm11,%%xmm1			\n\t		mulpd		    %%xmm3 ,%%xmm9 			\n\t"\
		"mulpd		    %%xmm11,%%xmm6			\n\t		mulpd		    %%xmm3 ,%%xmm14			\n\t"\
		"mulpd		    %%xmm3 ,%%xmm2			\n\t		mulpd		    %%xmm11,%%xmm10			\n\t"\
		"mulpd		    %%xmm11,%%xmm7			\n\t		mulpd		    %%xmm3 ,%%xmm15			\n\t"\
		"mulpd	   0x070(%%rax),%%xmm3			\n\t		mulpd	   0x170(%%rax),%%xmm11			\n\t"\
		"subpd		%%xmm6,%%xmm5				\n\t		subpd		%%xmm14,%%xmm13				\n\t"\
		"subpd		%%xmm2,%%xmm1				\n\t		subpd		%%xmm10,%%xmm9 				\n\t"\
		"addpd		%%xmm7,%%xmm4				\n\t		addpd		%%xmm15,%%xmm12				\n\t"\
		"addpd		%%xmm3,%%xmm0				\n\t		addpd		%%xmm11,%%xmm8 				\n\t"\
		"movaps		%%xmm5,%%xmm7				\n\t		movaps		%%xmm13,%%xmm15				\n\t"\
		"movaps		%%xmm4,%%xmm6				\n\t		movaps		%%xmm12,%%xmm14				\n\t"\
		"addpd		%%xmm0,%%xmm4				\n\t		addpd		%%xmm8 ,%%xmm12				\n\t"\
		"addpd		%%xmm1,%%xmm5				\n\t		addpd		%%xmm9 ,%%xmm13				\n\t"\
		"subpd		%%xmm0,%%xmm6				\n\t		subpd		%%xmm8 ,%%xmm14				\n\t"\
		"subpd		%%xmm1,%%xmm7				\n\t		subpd		%%xmm9 ,%%xmm15				\n\t"\
		"movaps		0x040(%%rax),%%xmm2			\n\t		movaps		0x140(%%rax),%%xmm10		\n\t"\
		"movaps		0x050(%%rax),%%xmm3			\n\t		movaps		0x150(%%rax),%%xmm11		\n\t"\
		"movaps		     (%%rax),%%xmm0			\n\t		movaps		0x100(%%rax),%%xmm8 		\n\t"\
		"movaps		0x010(%%rax),%%xmm1			\n\t		movaps		0x110(%%rax),%%xmm9 		\n\t"\
		"addpd		0x050(%%rax),%%xmm2			\n\t		subpd		0x150(%%rax),%%xmm10		\n\t"\
		"subpd		0x040(%%rax),%%xmm3			\n\t		addpd		0x140(%%rax),%%xmm11		\n\t"\
		"mulpd		(%%rbx),%%xmm2				\n\t		mulpd		(%%rbx),%%xmm10				\n\t"\
		"mulpd		(%%rbx),%%xmm3				\n\t		mulpd		(%%rbx),%%xmm11				\n\t"\
		"subpd		%%xmm2,%%xmm0				\n\t		subpd		%%xmm10,%%xmm8 				\n\t"\
		"subpd		%%xmm3,%%xmm1				\n\t		subpd		%%xmm11,%%xmm9 				\n\t"\
		"addpd		%%xmm2,%%xmm2				\n\t		addpd		%%xmm10,%%xmm10				\n\t"\
		"addpd		%%xmm3,%%xmm3				\n\t		addpd		%%xmm11,%%xmm11				\n\t"\
		"addpd		%%xmm0,%%xmm2				\n\t		addpd		%%xmm8 ,%%xmm10				\n\t"\
		"addpd		%%xmm1,%%xmm3				\n\t		addpd		%%xmm9 ,%%xmm11				\n\t"\
		"movq		%[__add1],%%rbx				\n\t"\
		"movq		%[__c1],%%rcx				\n\t		"/*movq		%[__c3],%%rcx	c3  = c1 + 8 */\
		"movq		%[__c9],%%rdx				\n\t		"/*movq		%[__c11],%%rdx	c11 = c9 + 8 */\
		"subpd		%%xmm4,%%xmm2				\n\t		subpd		%%xmm14,%%xmm8 				\n\t"\
		"subpd		%%xmm5,%%xmm3				\n\t		subpd		%%xmm15,%%xmm9 				\n\t"\
		"addpd		%%xmm4,%%xmm4				\n\t		addpd		%%xmm14,%%xmm14				\n\t"\
		"addpd		%%xmm5,%%xmm5				\n\t		addpd		%%xmm15,%%xmm15				\n\t"\
		"addpd		%%xmm2,%%xmm4				\n\t		addpd		%%xmm8 ,%%xmm14				\n\t"\
		"addpd		%%xmm3,%%xmm5				\n\t		addpd		%%xmm9 ,%%xmm15				\n\t"\
	"prefetcht1	-0x100(%%r14)\n\t"\
		"movaps		%%xmm2,     (%%rax)			\n\t		movaps		%%xmm8 ,0x100(%%rax)		\n\t"\
		"movaps		%%xmm3,0x010(%%rax)			\n\t		movaps		%%xmm9 ,0x110(%%rax)		\n\t"\
		"movaps		%%xmm4,%%xmm2				\n\t		movaps		%%xmm14,%%xmm8 				\n\t"\
		"movaps		%%xmm5,%%xmm3				\n\t		movaps		%%xmm15,%%xmm9 				\n\t"\
		"mulpd		    (%%rcx),%%xmm4			\n\t		mulpd		0x80(%%rcx),%%xmm14			\n\t"\
		"mulpd		    (%%rcx),%%xmm5			\n\t		mulpd		0x80(%%rcx),%%xmm15			\n\t"\
		"mulpd		0x10(%%rcx),%%xmm2			\n\t		mulpd		0x90(%%rcx),%%xmm8 			\n\t"\
		"mulpd		0x10(%%rcx),%%xmm3			\n\t		mulpd		0x90(%%rcx),%%xmm9 			\n\t"\
		"subpd		%%xmm2,%%xmm5				\n\t		subpd		%%xmm8 ,%%xmm15				\n\t"\
		"addpd		%%xmm3,%%xmm4				\n\t		addpd		%%xmm9 ,%%xmm14				\n\t"\
/* Tmps for later unpack into add1+0,1,4,5,8,9,c,d: */	/* Tmps for later unpack into add1+2,3,6,7,a,b,e,f: */\
		"movaps		%%xmm5,0x10(%%rbx)			\n\t		movaps		%%xmm15,0x30(%%rbx)			\n\t"\
		"movaps		%%xmm4,    (%%rbx)			\n\t		movaps		%%xmm14,0x20(%%rbx)			\n\t"\
		"movaps		     (%%rax),%%xmm4			\n\t		movaps		0x100(%%rax),%%xmm14		\n\t"\
		"movaps		0x010(%%rax),%%xmm5			\n\t		movaps		0x110(%%rax),%%xmm15		\n\t"\
		"movaps		%%xmm4,%%xmm2				\n\t		movaps		%%xmm14,%%xmm8 				\n\t"\
		"movaps		%%xmm5,%%xmm3				\n\t		movaps		%%xmm15,%%xmm9 				\n\t"\
		"mulpd		    (%%rdx),%%xmm4			\n\t		mulpd		0x80(%%rdx),%%xmm14			\n\t"\
		"mulpd		    (%%rdx),%%xmm5			\n\t		mulpd		0x80(%%rdx),%%xmm15			\n\t"\
		"mulpd		0x10(%%rdx),%%xmm2			\n\t		mulpd		0x90(%%rdx),%%xmm8 			\n\t"\
		"mulpd		0x10(%%rdx),%%xmm3			\n\t		mulpd		0x90(%%rdx),%%xmm9 			\n\t"\
		"subpd		%%xmm2,%%xmm5				\n\t		subpd		%%xmm8 ,%%xmm15				\n\t"\
		"addpd		%%xmm3,%%xmm4				\n\t		addpd		%%xmm9 ,%%xmm14				\n\t"\
		"movaps		%%xmm5,0x90(%%rbx)			\n\t		movaps		%%xmm15,0xb0(%%rbx)			\n\t"\
		"movaps		%%xmm4,0x80(%%rbx)			\n\t		movaps		%%xmm14,0xa0(%%rbx)			\n\t"\
		"movq		%[__c5],%%rcx				\n\t		"/*movq		%[__c7],%%rcx	c7  = c5  + 8 */\
		"movq		%[__c13],%%rdx				\n\t		"/*movq		%[__c15],%%rdx	c15 = c13 + 8 */\
		"subpd		%%xmm7,%%xmm0				\n\t		subpd		%%xmm13,%%xmm10				\n\t"\
		"subpd		%%xmm6,%%xmm1				\n\t		subpd		%%xmm12,%%xmm11				\n\t"\
		"addpd		%%xmm7,%%xmm7				\n\t		addpd		%%xmm13,%%xmm13				\n\t"\
		"addpd		%%xmm6,%%xmm6				\n\t		addpd		%%xmm12,%%xmm12				\n\t"\
		"addpd		%%xmm0,%%xmm7				\n\t		addpd		%%xmm10,%%xmm13				\n\t"\
		"addpd		%%xmm1,%%xmm6				\n\t		addpd		%%xmm11,%%xmm12				\n\t"\
		"movaps		%%xmm7,%%xmm4				\n\t		movaps		%%xmm13,%%xmm8 				\n\t"\
		"movaps		%%xmm1,%%xmm5				\n\t		movaps		%%xmm11,%%xmm9 				\n\t"\
		"mulpd		    (%%rcx),%%xmm7			\n\t		mulpd		0x80(%%rcx),%%xmm13			\n\t"\
		"mulpd		    (%%rcx),%%xmm1			\n\t		mulpd		0x80(%%rcx),%%xmm11			\n\t"\
		"mulpd		0x10(%%rcx),%%xmm4			\n\t		mulpd		0x90(%%rcx),%%xmm8 			\n\t"\
		"mulpd		0x10(%%rcx),%%xmm5			\n\t		mulpd		0x90(%%rcx),%%xmm9 			\n\t"\
		"subpd		%%xmm4,%%xmm1				\n\t		subpd		%%xmm8 ,%%xmm11				\n\t"\
		"addpd		%%xmm5,%%xmm7				\n\t		addpd		%%xmm9 ,%%xmm13				\n\t"\
		"movaps		%%xmm1,0x50(%%rbx)			\n\t		movaps		%%xmm11,0x70(%%rbx)			\n\t"\
		"movaps		%%xmm7,0x40(%%rbx)			\n\t		movaps		%%xmm13,0x60(%%rbx)			\n\t"\
		"movaps		%%xmm0,%%xmm4				\n\t		movaps		%%xmm10,%%xmm8 				\n\t"\
		"movaps		%%xmm6,%%xmm5				\n\t		movaps		%%xmm12,%%xmm9 				\n\t"\
		"mulpd		    (%%rdx),%%xmm0			\n\t		mulpd		0x80(%%rdx),%%xmm10			\n\t"\
		"mulpd		    (%%rdx),%%xmm6			\n\t		mulpd		0x80(%%rdx),%%xmm12			\n\t"\
		"mulpd		0x10(%%rdx),%%xmm4			\n\t		mulpd		0x90(%%rdx),%%xmm8 			\n\t"\
		"mulpd		0x10(%%rdx),%%xmm5			\n\t		mulpd		0x90(%%rdx),%%xmm9 			\n\t"\
	"prefetcht1	-0x140(%%r14)\n\t"\
		"subpd		%%xmm4,%%xmm6				\n\t		subpd		%%xmm8 ,%%xmm12				\n\t"\
		"addpd		%%xmm5,%%xmm0				\n\t		addpd		%%xmm9 ,%%xmm10				\n\t"\
		"movaps		%%xmm6,0xd0(%%rbx)			\n\t		movaps		%%xmm12,0xf0(%%rbx)			\n\t"\
		"movaps		%%xmm0,0xc0(%%rbx)			\n\t		movaps		%%xmm10,0xe0(%%rbx)			\n\t"\
		/*...Block 1: t1,9,17,25 -> r1,5,3,7 */				/*...Block 3: t5,13,21,29 -> r17,21,19,23: */\
		"movq		%[__isrt2],%%rbx			\n\t"\
		"movq		%[__r1],%%rax				\n\t		movaps		0x120(%%rax),%%xmm12		\n\t"\
		"movaps		     (%%rax),%%xmm0			\n\t		movaps		0x130(%%rax),%%xmm13		\n\t"\
		"movaps		0x010(%%rax),%%xmm1			\n\t		movaps		0x160(%%rax),%%xmm8 		\n\t"\
		"movaps		0x040(%%rax),%%xmm2			\n\t		movaps		0x170(%%rax),%%xmm9 		\n\t"\
		"movaps		0x050(%%rax),%%xmm3			\n\t		movaps		(%%rbx),%%xmm10	/* isrt2 */	\n\t"\
		"subpd		0x040(%%rax),%%xmm0			\n\t		addpd		0x130(%%rax),%%xmm12		\n\t"\
		"subpd		0x050(%%rax),%%xmm1			\n\t		subpd		0x120(%%rax),%%xmm13		\n\t"\
		"addpd		     (%%rax),%%xmm2			\n\t		subpd		0x170(%%rax),%%xmm8 		\n\t"\
		"addpd		0x010(%%rax),%%xmm3			\n\t		addpd		0x160(%%rax),%%xmm9 		\n\t"\
		"movaps		0x020(%%rax),%%xmm4			\n\t		mulpd		%%xmm10,%%xmm12				\n\t"\
		"movaps		0x030(%%rax),%%xmm5			\n\t		mulpd		%%xmm10,%%xmm13				\n\t"\
		"movaps		0x060(%%rax),%%xmm6			\n\t		mulpd		%%xmm10,%%xmm8 				\n\t"\
		"movaps		0x070(%%rax),%%xmm7			\n\t		mulpd		%%xmm10,%%xmm9 				\n\t"\
		"subpd		0x060(%%rax),%%xmm4			\n\t		movaps		%%xmm12,%%xmm14				\n\t"\
		"subpd		0x070(%%rax),%%xmm5			\n\t		movaps		%%xmm13,%%xmm15				\n\t"\
		"addpd		0x020(%%rax),%%xmm6			\n\t		subpd		%%xmm8 ,%%xmm12				\n\t"\
		"addpd		0x030(%%rax),%%xmm7			\n\t		subpd		%%xmm9 ,%%xmm13				\n\t"\
		"movq		%[__add0],%%rbx				\n\t		addpd		%%xmm8 ,%%xmm14				\n\t"\
		"movq		%[__c8],%%rdx				\n\t		addpd		%%xmm9 ,%%xmm15				\n\t"\
		"addpd		%%xmm6,%%xmm2				\n\t		movaps		0x100(%%rax),%%xmm8 		\n\t"\
		"addpd		%%xmm7,%%xmm3				\n\t		movaps		0x110(%%rax),%%xmm9 		\n\t"\
		"movaps		%%xmm2,    (%%rbx)			\n\t		movaps		0x140(%%rax),%%xmm10		\n\t"\
		"movaps		%%xmm3,0x10(%%rbx)			\n\t		movaps		0x150(%%rax),%%xmm11		\n\t"\
		"addpd		%%xmm6,%%xmm6				\n\t		subpd		0x150(%%rax),%%xmm8 		\n\t"\
		"addpd		%%xmm7,%%xmm7				\n\t		subpd		0x140(%%rax),%%xmm9 		\n\t"\
		"subpd		%%xmm6,%%xmm2				\n\t		addpd		0x100(%%rax),%%xmm11		\n\t"\
		"subpd		%%xmm7,%%xmm3				\n\t		addpd		0x110(%%rax),%%xmm10		\n\t"\
		"movaps		%%xmm2,%%xmm6				\n\t		subpd		%%xmm12,%%xmm11				\n\t"\
		"movaps		%%xmm3,%%xmm7				\n\t		subpd		%%xmm13,%%xmm9 				\n\t"\
		"mulpd		    (%%rdx),%%xmm2			\n\t		addpd		%%xmm12,%%xmm12				\n\t"\
		"mulpd		    (%%rdx),%%xmm3			\n\t		addpd		%%xmm13,%%xmm13				\n\t"\
		"mulpd		0x10(%%rdx),%%xmm6			\n\t		addpd		%%xmm11,%%xmm12				\n\t"\
		"mulpd		0x10(%%rdx),%%xmm7			\n\t		addpd		%%xmm9 ,%%xmm13				\n\t"\
		"subpd		%%xmm6,%%xmm3				\n\t		movq %[__c10],%%rdx	\n\t"/*movq %[__c2],%%rcx	Use c2 = c10-2 since lcol needs rcx */\
		"addpd		%%xmm7,%%xmm2				\n\t		movaps		%%xmm11,0x100(%%rax)		\n\t"\
		"movq		%[__add1],%%rcx				\n\t		movaps		%%xmm9 ,0x110(%%rax)		\n\t"\
	/*********************************************************************************************/\
	/************ ALL STORES FROM HERE ON ARE FINAL **********************************************/\
	/*********************************************************************************************/\
		"movaps		%%xmm3,%%xmm7				\n\t		movaps		%%xmm12,%%xmm11				\n\t"\
		"movaps		%%xmm2,%%xmm6				\n\t		movaps		%%xmm13,%%xmm9 				\n\t"\
		"unpckhpd	0x90(%%rcx),%%xmm7			\n\t		mulpd		-0x20(%%rdx),%%xmm12		\n\t"\
		"unpcklpd	0x90(%%rcx),%%xmm3			\n\t		mulpd		-0x20(%%rdx),%%xmm13		\n\t"\
	"prefetcht1	-0x180(%%r14)\n\t"\
		"movaps		%%xmm7,0x90(%%rcx)			\n\t		mulpd		-0x10(%%rdx),%%xmm11		\n\t"\
		"unpckhpd	0x80(%%rcx),%%xmm6			\n\t		mulpd		-0x10(%%rdx),%%xmm9 		\n\t"\
		"unpcklpd	0x80(%%rcx),%%xmm2			\n\t		subpd		%%xmm11,%%xmm13				\n\t"\
		"movaps		%%xmm6,0x80(%%rcx)			\n\t		addpd		%%xmm9 ,%%xmm12				\n\t"\
		"movaps		%%xmm3,0x90(%%rbx)			\n\t		movaps		%%xmm13,%%xmm11				\n\t"\
		"movaps		%%xmm2,0x80(%%rbx)			\n\t		movaps		%%xmm12,%%xmm9 				\n\t"\
		"movaps		0x10(%%rbx),%%xmm3			\n\t		unpckhpd	0x30(%%rcx),%%xmm11			\n\t"\
		"movaps		    (%%rbx),%%xmm2			\n\t		unpcklpd	0x30(%%rcx),%%xmm13			\n\t"\
		"movaps		%%xmm3,%%xmm7				\n\t		movaps		%%xmm11,0x30(%%rcx)			\n\t"\
		"movaps		%%xmm2,%%xmm6				\n\t		unpckhpd	0x20(%%rcx),%%xmm9 			\n\t"\
		"unpckhpd	0x10(%%rcx),%%xmm7			\n\t		unpcklpd	0x20(%%rcx),%%xmm12			\n\t"\
		"unpcklpd	0x10(%%rcx),%%xmm3			\n\t		movaps		%%xmm9 ,0x20(%%rcx)			\n\t"\
		"movaps		%%xmm7,0x10(%%rcx)			\n\t		movaps		%%xmm13,0x30(%%rbx)			\n\t"\
		"unpckhpd	    (%%rcx),%%xmm6			\n\t		movaps		%%xmm12,0x20(%%rbx)			\n\t"\
		"unpcklpd	    (%%rcx),%%xmm2			\n\t		movaps		0x100(%%rax),%%xmm12		\n\t"\
		"movaps		%%xmm6,    (%%rcx)			\n\t		movaps		0x110(%%rax),%%xmm13		\n\t"\
		"movaps		%%xmm3,0x10(%%rbx)			\n\t		movaps		%%xmm12,%%xmm11				\n\t"\
		"movaps		%%xmm2,    (%%rbx)			\n\t		movaps		%%xmm13,%%xmm9 				\n\t"\
		"addpd		%%xmm5,%%xmm0				\n\t		mulpd		    (%%rdx),%%xmm12			\n\t"\
		"subpd		%%xmm4,%%xmm1				\n\t		mulpd		    (%%rdx),%%xmm13			\n\t"\
		"movaps		%%xmm0,%%xmm2				\n\t		mulpd		0x10(%%rdx),%%xmm11			\n\t"\
		"movaps		%%xmm1,%%xmm3				\n\t		mulpd		0x10(%%rdx),%%xmm9 			\n\t"\
		"movq %[__c12],%%rdx	/* c4 = c12-2 */\n\t		subpd		%%xmm11,%%xmm13				\n\t"\
		"addpd		%%xmm5,%%xmm5				\n\t		addpd		%%xmm9 ,%%xmm12				\n\t"\
		"addpd		%%xmm4,%%xmm4				\n\t		movaps		%%xmm13,%%xmm11				\n\t"\
		"movaps		%%xmm0,%%xmm6				\n\t		movaps		%%xmm12,%%xmm9 				\n\t"\
		"movaps		%%xmm1,%%xmm7				\n\t		unpckhpd	0xb0(%%rcx),%%xmm11			\n\t"\
		"mulpd		-0x20(%%rdx),%%xmm2			\n\t		unpcklpd	0xb0(%%rcx),%%xmm13			\n\t"\
		"mulpd		-0x20(%%rdx),%%xmm3			\n\t		movaps		%%xmm11,0xb0(%%rcx)			\n\t"\
		"mulpd		-0x10(%%rdx),%%xmm6			\n\t		unpckhpd	0xa0(%%rcx),%%xmm9 			\n\t"\
		"mulpd		-0x10(%%rdx),%%xmm7			\n\t		unpcklpd	0xa0(%%rcx),%%xmm12			\n\t"\
		"subpd		%%xmm6,%%xmm3				\n\t		movaps		%%xmm9 ,0xa0(%%rcx)			\n\t"\
		"addpd		%%xmm7,%%xmm2				\n\t		movaps		%%xmm13,0xb0(%%rbx)			\n\t"\
		"movaps		%%xmm3,%%xmm7				\n\t		movaps		%%xmm12,0xa0(%%rbx)			\n\t"\
		"movaps		%%xmm2,%%xmm6				\n\t		subpd		%%xmm15,%%xmm8 				\n\t"\
		"unpckhpd	0x50(%%rcx),%%xmm7			\n\t		subpd		%%xmm14,%%xmm10				\n\t"\
		"unpcklpd	0x50(%%rcx),%%xmm3			\n\t		addpd		%%xmm15,%%xmm15				\n\t"\
		"movaps		%%xmm7,0x50(%%rcx)			\n\t		addpd		%%xmm14,%%xmm14				\n\t"\
		"unpckhpd	0x40(%%rcx),%%xmm6			\n\t		addpd		%%xmm8 ,%%xmm15				\n\t"\
		"unpcklpd	0x40(%%rcx),%%xmm2			\n\t		addpd		%%xmm10,%%xmm14				\n\t"\
		"movaps		%%xmm6,0x40(%%rcx)			\n\t		movaps		%%xmm15,%%xmm12				\n\t"\
		"movaps		%%xmm3,0x50(%%rbx)			\n\t		movaps		%%xmm10,%%xmm13				\n\t"\
		"movaps		%%xmm2,0x40(%%rbx)			\n\t		mulpd		0x60(%%rdx),%%xmm15	/* c6  = c12 + 6 */\n\t"\
		"subpd		%%xmm5,%%xmm0				\n\t		mulpd		0x60(%%rdx),%%xmm10			\n\t"\
		"addpd		%%xmm4,%%xmm1				\n\t		mulpd		0x70(%%rdx),%%xmm12			\n\t"\
		"movaps		%%xmm0,%%xmm6				\n\t		mulpd		0x70(%%rdx),%%xmm13			\n\t"\
	"prefetcht1	-0x1c0(%%r14)\n\t"\
		"movaps		%%xmm1,%%xmm7				\n\t		subpd		%%xmm12,%%xmm10				\n\t"\
		"mulpd		    (%%rdx),%%xmm0			\n\t		addpd		%%xmm13,%%xmm15				\n\t"\
		"mulpd		    (%%rdx),%%xmm1			\n\t		movaps		%%xmm10,%%xmm13				\n\t"\
		"mulpd		0x10(%%rdx),%%xmm6			\n\t		movaps		%%xmm15,%%xmm12				\n\t"\
		"mulpd		0x10(%%rdx),%%xmm7			\n\t		unpckhpd	0x70(%%rcx),%%xmm13			\n\t"\
		"subpd		%%xmm6,%%xmm1				\n\t		unpcklpd	0x70(%%rcx),%%xmm10			\n\t"\
		"addpd		%%xmm7,%%xmm0				\n\t		movaps		%%xmm13,0x70(%%rcx)			\n\t"\
		"movaps		%%xmm1,%%xmm7				\n\t		unpckhpd	0x60(%%rcx),%%xmm12			\n\t"\
		"movaps		%%xmm0,%%xmm6				\n\t		unpcklpd	0x60(%%rcx),%%xmm15			\n\t"\
		"unpckhpd	0xd0(%%rcx),%%xmm7			\n\t		movaps		%%xmm12,0x60(%%rcx)			\n\t"\
		"unpcklpd	0xd0(%%rcx),%%xmm1			\n\t		movaps		%%xmm10,0x70(%%rbx)			\n\t"\
		"movaps		%%xmm7,0xd0(%%rcx)			\n\t		movaps		%%xmm15,0x60(%%rbx)			\n\t"\
		"unpckhpd	0xc0(%%rcx),%%xmm6			\n\t		movaps		%%xmm8 ,%%xmm12				\n\t"\
		"unpcklpd	0xc0(%%rcx),%%xmm0			\n\t		movaps		%%xmm14,%%xmm13				\n\t"\
		"movaps		%%xmm6,0xc0(%%rcx)			\n\t		mulpd		0x80(%%rdx),%%xmm8 	/* c14 = c12 + 8 */\n\t"\
		"movaps		%%xmm1,0xd0(%%rbx)			\n\t		mulpd		0x80(%%rdx),%%xmm14			\n\t"\
		"movaps		%%xmm0,0xc0(%%rbx)			\n\t		mulpd		0x90(%%rdx),%%xmm12			\n\t"\
		"													mulpd		0x90(%%rdx),%%xmm13			\n\t"\
		"													subpd		%%xmm12,%%xmm14				\n\t"\
		"													addpd		%%xmm13,%%xmm8 				\n\t"\
		"													movaps		%%xmm14,%%xmm13				\n\t"\
		"													movaps		%%xmm8 ,%%xmm12				\n\t"\
		"													unpckhpd	0xf0(%%rcx),%%xmm13			\n\t"\
		"													unpcklpd	0xf0(%%rcx),%%xmm14			\n\t"\
		"													movaps		%%xmm13,0xf0(%%rcx)			\n\t"\
		"													unpckhpd	0xe0(%%rcx),%%xmm12			\n\t"\
		"													unpcklpd	0xe0(%%rcx),%%xmm8 			\n\t"\
		"													movaps		%%xmm12,0xe0(%%rcx)			\n\t"\
		"													movaps		%%xmm14,0xf0(%%rbx)			\n\t"\
		"													movaps		%%xmm8 ,0xe0(%%rbx)			\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__add1] "m" (Xadd1)\
		 ,[__r1] "m" (Xr1)\
		 ,[__r9] "m" (Xr9)\
		 ,[__r17] "m" (Xr17)\
		 ,[__r25] "m" (Xr25)\
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__cc0] "m" (Xcc0)\
		 ,[__c1] "m" (Xc1)\
		 ,[__c2] "m" (Xc2)\
		 ,[__c3] "m" (Xc3)\
		 ,[__c4] "m" (Xc4)\
		 ,[__c5] "m" (Xc5)\
		 ,[__c6] "m" (Xc6)\
		 ,[__c7] "m" (Xc7)\
		 ,[__c8] "m" (Xc8)\
		 ,[__c9] "m" (Xc9)\
		 ,[__c10] "m" (Xc10)\
		 ,[__c11] "m" (Xc11)\
		 ,[__c12] "m" (Xc12)\
		 ,[__c13] "m" (Xc13)\
		 ,[__c14] "m" (Xc14)\
		 ,[__c15] "m" (Xc15)\
		: "cc","memory","rax","rbx","rcx","rdx","r14","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"		/* Clobbered registers */\
	);\
	}

#endif	// SSE2 or AVX?

#endif	/* radix16_wrapper_square_gcc_h_included */

