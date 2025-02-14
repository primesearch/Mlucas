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
#ifndef radix32_wrapper_square_gcc_h_included
#define radix32_wrapper_square_gcc_h_included

#ifdef USE_ARM_V8_SIMD

	#define SSE2_RADIX32_WRAPPER_DIF(Xadd0,Xadd1,Xr00,Xr10,Xr20,Xr30,Xisrt2,Xcc0,Xc00,Xc01,Xc02,Xc03,Xc05,Xc07)\
	{\
	__asm__ volatile (\
		/************************************************************************/\
		/* Forward DIF radix-32 pass on the interleaved block1 and block2 data: */\
		/************************************************************************/\
		"ldr	x2,%[__r00]			\n\t	ldr w11,=0x10				\n\t"/* NB: ~0, not -: */\
		"add	x10, x2,#0x480		\n\t	ld1r	{v29.2d},[x10],x11	\n\t"/* x10 has isrt2-ptr for LD1R, gets post-incr'ed by 0x10 to cc0 */\
		/* Block 0:					rcol vs lcol: Input-addresses += 0x40, twiddle-addrs += 0x80, tmp-store-addrs += 0x80: */\
		"ldr	x0,%[__add0]		\n\t	ldp	q8 ,q9 ,[x0,#0x040]	\n\t"\
		"ldr	x1,%[__add1]		\n\t	ldp	q24,q25,[x1,#0x040]	\n\t"\
		"ldr	x3,%[__c00]			\n\t	ldp	q22,q23,[x3,#0x80]	\n\t"\
		/* vector-data-to-be-interleaved not contiguous, but can load each of the 2 vector [Re,Im]-pairs via LDP and
		TRN1 (ARM analog of x86 UNPCKLPD) and TRN2 (ARM analog of x86 UNPCKHPD) to do the needed 2x2 transposes: */\
		"									trn1	v26.2d,v8.2d,v24.2d	\n\t"\
		"ldp	q0 ,q1 ,[x0       ]	\n\t	trn1	v27.2d,v9.2d,v25.2d	\n\t"\
		"ldp	q18,q19,[x1       ]	\n\t	trn2	v24.2d,v8.2d,v24.2d	\n\t"\
/* c00 */"ldp	q16,q17,[x3]		\n\t	trn2	v25.2d,v9.2d,v25.2d	\n\t"\
		"trn1	v20.2d,v0.2d,v18.2d	\n\t	fmul	v8.2d,v26.2d,v22.2d	\n\t"\
		"trn1	v21.2d,v1.2d,v19.2d	\n\t	fmul	v9.2d,v27.2d,v22.2d	\n\t"\
		"trn2	v18.2d,v0.2d,v18.2d	\n\t	fmls	v8.2d,v27.2d,v23.2d	\n\t"\
		"trn2	v19.2d,v1.2d,v19.2d	\n\t	fmla	v9.2d,v26.2d,v23.2d	\n\t"\
		"fmul	v0.2d,v20.2d,v16.2d	\n\t	stp	q24,q25,[x2,#0x280]		\n\t"\
		"fmul	v1.2d,v21.2d,v16.2d	\n\t	ldp	q10,q11,[x0,#0x140]		\n\t"\
		"fmls	v0.2d,v21.2d,v17.2d	\n\t	ldp	q24,q25,[x1,#0x140]		\n\t"\
		"fmla	v1.2d,v20.2d,v17.2d	\n\t	ldp	q22,q23,[x3,#0xa0]		\n\t"\
		"stp	q18,q19,[x2,#0x200]	\n\t	trn1	v26.2d,v10.2d,v24.2d\n\t"\
		"ldp	q2 ,q3 ,[x0,#0x100]	\n\t	trn1	v27.2d,v11.2d,v25.2d\n\t"\
		"ldp	q18,q19,[x1,#0x100]	\n\t	trn2	v24.2d,v10.2d,v24.2d\n\t"\
		"ldp	q16,q17,[x3,#0x20]	\n\t	trn2	v25.2d,v11.2d,v25.2d\n\t"\
		"trn1	v20.2d,v2.2d,v18.2d	\n\t	fmul	v10.2d,v26.2d,v22.2d\n\t"\
		"trn1	v21.2d,v3.2d,v19.2d	\n\t	fmul	v11.2d,v27.2d,v22.2d\n\t"\
		"trn2	v18.2d,v2.2d,v18.2d	\n\t	fmls	v10.2d,v27.2d,v23.2d\n\t"\
		"trn2	v19.2d,v3.2d,v19.2d	\n\t	fmla	v11.2d,v26.2d,v23.2d\n\t"\
		"fmul	v2.2d,v20.2d,v16.2d	\n\t	stp	q24,q25,[x2,#0x2a0]		\n\t"\
		"fmul	v3.2d,v21.2d,v16.2d	\n\t	ldp	q12,q13,[x0,#0x0c0]		\n\t"\
		"fmls	v2.2d,v21.2d,v17.2d	\n\t	ldp	q24,q25,[x1,#0x0c0]		\n\t"\
		"fmla	v3.2d,v20.2d,v17.2d	\n\t	ldp	q22,q23,[x3,#0xc0]		\n\t"\
		"stp	q18,q19,[x2,#0x220]	\n\t	trn1	v26.2d,v12.2d,v24.2d\n\t"\
		"ldp	q4 ,q5 ,[x0,#0x080]	\n\t	trn1	v27.2d,v13.2d,v25.2d\n\t"\
		"ldp	q18,q19,[x1,#0x080]	\n\t	trn2	v24.2d,v12.2d,v24.2d\n\t"\
		"ldp	q16,q17,[x3,#0x40]	\n\t	trn2	v25.2d,v13.2d,v25.2d\n\t"\
		"trn1	v20.2d,v4.2d,v18.2d	\n\t	fmul	v12.2d,v26.2d,v22.2d\n\t"\
		"trn1	v21.2d,v5.2d,v19.2d	\n\t	fmul	v13.2d,v27.2d,v22.2d\n\t"\
		"trn2	v18.2d,v4.2d,v18.2d	\n\t	fmls	v12.2d,v27.2d,v23.2d\n\t"\
		"trn2	v19.2d,v5.2d,v19.2d	\n\t	fmla	v13.2d,v26.2d,v23.2d\n\t"\
		"fmul	v4.2d,v20.2d,v16.2d	\n\t	stp	q24,q25,[x2,#0x2c0]		\n\t"\
		"fmul	v5.2d,v21.2d,v16.2d	\n\t	ldp	q14,q15,[x0,#0x1c0]		\n\t"\
		"fmls	v4.2d,v21.2d,v17.2d	\n\t	ldp	q24,q25,[x1,#0x1c0]		\n\t"\
		"fmla	v5.2d,v20.2d,v17.2d	\n\t	ldp	q22,q23,[x3,#0xe0]		\n\t"\
		"stp	q18,q19,[x2,#0x240]	\n\t	trn1	v26.2d,v14.2d,v24.2d\n\t"\
		"ldp	q6 ,q7 ,[x0,#0x180]	\n\t	trn1	v27.2d,v15.2d,v25.2d\n\t"\
		"ldp	q18,q19,[x1,#0x180]	\n\t	trn2	v24.2d,v14.2d,v24.2d\n\t"\
		"ldp	q16,q17,[x3,#0x60]	\n\t	trn2	v25.2d,v15.2d,v25.2d\n\t"\
		"trn1	v20.2d,v6.2d,v18.2d	\n\t	fmul	v14.2d,v26.2d,v22.2d\n\t"\
		"trn1	v21.2d,v7.2d,v19.2d	\n\t	fmul	v15.2d,v27.2d,v22.2d\n\t"\
		"trn2	v18.2d,v6.2d,v18.2d	\n\t	fmls	v14.2d,v27.2d,v23.2d\n\t"\
		"trn2	v19.2d,v7.2d,v19.2d	\n\t	fmla	v15.2d,v26.2d,v23.2d\n\t"\
		"fmul	v6.2d,v20.2d,v16.2d	\n\t	stp	q24,q25,[x2,#0x2e0]		\n\t"\
		"fmul	v7.2d,v21.2d,v16.2d	\n\t	fadd	v22.2d,v8.2d,v10.2d	\n\t"/* 2x2 butterflies: */\
		"fmls	v6.2d,v21.2d,v17.2d	\n\t	fadd	v23.2d,v9.2d,v11.2d	\n\t"\
		"fmla	v7.2d,v20.2d,v17.2d	\n\t	fsub	v26.2d,v8.2d,v10.2d	\n\t"\
		"stp	q18,q19,[x2,#0x260]	\n\t	fsub	v27.2d,v9.2d,v11.2d	\n\t"\
/*0*/	"fadd	v16.2d,v0.2d,v2.2d	\n\t	fsub	v24.2d,v12.2d,v14.2d\n\t"\
/*1*/	"fadd	v17.2d,v1.2d,v3.2d	\n\t	fsub	v25.2d,v13.2d,v15.2d\n\t"\
/*2*/	"fsub	v20.2d,v0.2d,v2.2d	\n\t	fadd	v14.2d,v12.2d,v14.2d\n\t"\
/*3*/	"fsub	v21.2d,v1.2d,v3.2d	\n\t	fadd	v15.2d,v13.2d,v15.2d\n\t"\
/*4*/	"fsub	v18.2d,v4.2d,v6.2d	\n\t	fsub	v8.2d ,v22.2d,v14.2d\n\t"/* combine to get 2 length-4 output subtransforms... */\
/*5*/	"fsub	v19.2d,v5.2d,v7.2d	\n\t	fsub	v9.2d ,v23.2d,v15.2d\n\t"\
		"fadd	v6.2d,v4.2d,v6.2d	\n\t	fadd	v14.2d,v22.2d,v14.2d\n\t"\
		"fadd	v7.2d,v5.2d,v7.2d	\n\t	fadd	v15.2d,v23.2d,v15.2d\n\t"\
		"fsub	v0.2d,v16.2d,v6.2d	\n\t	fsub	v10.2d,v26.2d,v25.2d\n\t"\
		"fsub	v1.2d,v17.2d,v7.2d	\n\t	fsub	v11.2d,v27.2d,v24.2d\n\t"\
		"fadd	v6.2d,v16.2d,v6.2d	\n\t	fadd	v25.2d,v26.2d,v25.2d\n\t"/* 13 */\
		"fadd	v7.2d,v17.2d,v7.2d	\n\t	fadd	v24.2d,v27.2d,v24.2d\n\t"/* 12 */\
		"fsub	v2.2d,v20.2d,v19.2d	\n\t	fadd	v12.2d,v10.2d,v24.2d\n\t"/* 12=10+12 */\
		"fsub	v3.2d,v21.2d,v18.2d	\n\t	fsub	v13.2d,v25.2d,v11.2d\n\t"/* 13=13-11 */\
		"fadd	v5.2d,v20.2d,v19.2d	\n\t	fsub	v10.2d,v10.2d,v24.2d\n\t"/* 10=10-12 */\
		"fadd	v4.2d,v21.2d,v18.2d	\n\t	fadd	v11.2d,v25.2d,v11.2d\n\t"/* 11=13+11 */\
/*7*/	"fsub	v16.2d,v7.2d,v15.2d	\n\t	fmul	v12.2d,v12.2d,v29.2d\n\t"/* 12 *= isrt2 */\
/*1*/	"fsub	v17.2d,v1.2d,v8.2d	\n\t	fmul	v13.2d,v13.2d,v29.2d\n\t"/* 13 *= isrt2 */\
/*6*/	"fsub	v18.2d,v6.2d,v14.2d	\n\t	fmul	v10.2d,v10.2d,v29.2d\n\t"/* 10 *= isrt2 */\
/*0*/	"fsub	v19.2d,v0.2d,v9.2d	\n\t	fmul	v11.2d,v11.2d,v29.2d\n\t"/* 11 *= isrt2 */\
		/* rcol: We are tempted to replace the [4 *=isrt2,8 addsub] by 8 FMA, but ARM fma's restricted
		addend-takes-output FMA syntax in conjunction with the scrambled output ordering wrecks that: */\
		"fadd	v15.2d,v7.2d,v15.2d	\n\t	fsub	v22.2d,v4.2d,v12.2d	\n\t"/*4*/\
		"fadd	v8.2d ,v1.2d,v8.2d	\n\t	fsub	v23.2d,v3.2d,v13.2d	\n\t"/*3*/\
		"fadd	v14.2d,v6.2d,v14.2d	\n\t	fsub	v24.2d,v2.2d,v10.2d	\n\t"/*2*/\
		"fadd	v9.2d ,v0.2d,v9.2d	\n\t	fsub	v25.2d,v5.2d,v11.2d	\n\t"/*5*/\
/*6,7*/	"stp	q18,q16,[x2,#0x80]	\n\t	fadd	v12.2d,v4.2d,v12.2d	\n\t"\
/*0,8*/	"stp	q19,q8 ,[x2,#0x40]	\n\t	fadd	v13.2d,v3.2d,v13.2d	\n\t"\
/*e,f*/	"stp	q14,q15,[x2      ]	\n\t	fadd	v10.2d,v2.2d,v10.2d	\n\t"\
/*9,1*/	"stp	q9 ,q17,[x2,#0xc0]	\n\t	fadd	v11.2d,v5.2d,v11.2d	\n\t"\
		"									stp	q24,q22,[x2,#0xa0]	\n\t"/*2,4*/\
		"									stp	q25,q13,[x2,#0x60]	\n\t"/*5,d*/\
		"									stp	q10,q12,[x2,#0x20]	\n\t"/*a,c*/\
		"									stp	q11,q23,[x2,#0xe0]	\n\t"/*b,3*/\
		/* Block 1: */\
		"add	x2, x2,#0x100		\n\t"\
		"add	x0, x0,#0x20		\n\t	ldp	q8 ,q9 ,[x0,#0x040]	\n\t"\
		"add	x1, x1,#0x20		\n\t	ldp	q24,q25,[x1,#0x040]	\n\t"\
		"add	x3, x3,#0x100		\n\t	ldp	q22,q23,[x3,#0x80]	\n\t"\
		"									trn1	v26.2d,v8.2d,v24.2d	\n\t"\
		"ldp	q0 ,q1 ,[x0       ]	\n\t	trn1	v27.2d,v9.2d,v25.2d	\n\t"\
		"ldp	q18,q19,[x1       ]	\n\t	trn2	v24.2d,v8.2d,v24.2d	\n\t"\
/* c02 */"ldp	q16,q17,[x3]		\n\t	trn2	v25.2d,v9.2d,v25.2d	\n\t"\
		"trn1	v20.2d,v0.2d,v18.2d	\n\t	fmul	v8.2d,v26.2d,v22.2d	\n\t"\
		"trn1	v21.2d,v1.2d,v19.2d	\n\t	fmul	v9.2d,v27.2d,v22.2d	\n\t"\
		"trn2	v18.2d,v0.2d,v18.2d	\n\t	fmls	v8.2d,v27.2d,v23.2d	\n\t"\
		"trn2	v19.2d,v1.2d,v19.2d	\n\t	fmla	v9.2d,v26.2d,v23.2d	\n\t"\
		"fmul	v0.2d,v20.2d,v16.2d	\n\t	stp	q24,q25,[x2,#0x280]		\n\t"\
		"fmul	v1.2d,v21.2d,v16.2d	\n\t	ldp	q10,q11,[x0,#0x140]		\n\t"\
		"fmls	v0.2d,v21.2d,v17.2d	\n\t	ldp	q24,q25,[x1,#0x140]		\n\t"\
		"fmla	v1.2d,v20.2d,v17.2d	\n\t	ldp	q22,q23,[x3,#0xa0]		\n\t"\
		"stp	q18,q19,[x2,#0x200]	\n\t	trn1	v26.2d,v10.2d,v24.2d\n\t"\
		"ldp	q2 ,q3 ,[x0,#0x100]	\n\t	trn1	v27.2d,v11.2d,v25.2d\n\t"\
		"ldp	q18,q19,[x1,#0x100]	\n\t	trn2	v24.2d,v10.2d,v24.2d\n\t"\
		"ldp	q16,q17,[x3,#0x20]	\n\t	trn2	v25.2d,v11.2d,v25.2d\n\t"\
		"trn1	v20.2d,v2.2d,v18.2d	\n\t	fmul	v10.2d,v26.2d,v22.2d\n\t"\
		"trn1	v21.2d,v3.2d,v19.2d	\n\t	fmul	v11.2d,v27.2d,v22.2d\n\t"\
		"trn2	v18.2d,v2.2d,v18.2d	\n\t	fmls	v10.2d,v27.2d,v23.2d\n\t"\
		"trn2	v19.2d,v3.2d,v19.2d	\n\t	fmla	v11.2d,v26.2d,v23.2d\n\t"\
		"fmul	v2.2d,v20.2d,v16.2d	\n\t	stp	q24,q25,[x2,#0x2a0]		\n\t"\
		"fmul	v3.2d,v21.2d,v16.2d	\n\t	ldp	q12,q13,[x0,#0x0c0]		\n\t"\
		"fmls	v2.2d,v21.2d,v17.2d	\n\t	ldp	q24,q25,[x1,#0x0c0]		\n\t"\
		"fmla	v3.2d,v20.2d,v17.2d	\n\t	ldp	q22,q23,[x3,#0xc0]		\n\t"\
		"stp	q18,q19,[x2,#0x220]	\n\t	trn1	v26.2d,v12.2d,v24.2d\n\t"\
		"ldp	q4 ,q5 ,[x0,#0x080]	\n\t	trn1	v27.2d,v13.2d,v25.2d\n\t"\
		"ldp	q18,q19,[x1,#0x080]	\n\t	trn2	v24.2d,v12.2d,v24.2d\n\t"\
		"ldp	q16,q17,[x3,#0x40]	\n\t	trn2	v25.2d,v13.2d,v25.2d\n\t"\
		"trn1	v20.2d,v4.2d,v18.2d	\n\t	fmul	v12.2d,v26.2d,v22.2d\n\t"\
		"trn1	v21.2d,v5.2d,v19.2d	\n\t	fmul	v13.2d,v27.2d,v22.2d\n\t"\
		"trn2	v18.2d,v4.2d,v18.2d	\n\t	fmls	v12.2d,v27.2d,v23.2d\n\t"\
		"trn2	v19.2d,v5.2d,v19.2d	\n\t	fmla	v13.2d,v26.2d,v23.2d\n\t"\
		"fmul	v4.2d,v20.2d,v16.2d	\n\t	stp	q24,q25,[x2,#0x2c0]		\n\t"\
		"fmul	v5.2d,v21.2d,v16.2d	\n\t	ldp	q14,q15,[x0,#0x1c0]		\n\t"\
		"fmls	v4.2d,v21.2d,v17.2d	\n\t	ldp	q24,q25,[x1,#0x1c0]		\n\t"\
		"fmla	v5.2d,v20.2d,v17.2d	\n\t	ldp	q22,q23,[x3,#0xe0]		\n\t"\
		"stp	q18,q19,[x2,#0x240]	\n\t	trn1	v26.2d,v14.2d,v24.2d\n\t"\
		"ldp	q6 ,q7 ,[x0,#0x180]	\n\t	trn1	v27.2d,v15.2d,v25.2d\n\t"\
		"ldp	q18,q19,[x1,#0x180]	\n\t	trn2	v24.2d,v14.2d,v24.2d\n\t"\
		"ldp	q16,q17,[x3,#0x60]	\n\t	trn2	v25.2d,v15.2d,v25.2d\n\t"\
		"trn1	v20.2d,v6.2d,v18.2d	\n\t	fmul	v14.2d,v26.2d,v22.2d\n\t"\
		"trn1	v21.2d,v7.2d,v19.2d	\n\t	fmul	v15.2d,v27.2d,v22.2d\n\t"\
		"trn2	v18.2d,v6.2d,v18.2d	\n\t	fmls	v14.2d,v27.2d,v23.2d\n\t"\
		"trn2	v19.2d,v7.2d,v19.2d	\n\t	fmla	v15.2d,v26.2d,v23.2d\n\t"\
		"fmul	v6.2d,v20.2d,v16.2d	\n\t	stp	q24,q25,[x2,#0x2e0]		\n\t"\
		"fmul	v7.2d,v21.2d,v16.2d	\n\t	fadd	v22.2d,v8.2d,v10.2d	\n\t"\
		"fmls	v6.2d,v21.2d,v17.2d	\n\t	fadd	v23.2d,v9.2d,v11.2d	\n\t"\
		"fmla	v7.2d,v20.2d,v17.2d	\n\t	fsub	v26.2d,v8.2d,v10.2d	\n\t"\
		"stp	q18,q19,[x2,#0x260]	\n\t	fsub	v27.2d,v9.2d,v11.2d	\n\t"\
		"fadd	v16.2d,v0.2d,v2.2d	\n\t	fsub	v24.2d,v12.2d,v14.2d\n\t"\
		"fadd	v17.2d,v1.2d,v3.2d	\n\t	fsub	v25.2d,v13.2d,v15.2d\n\t"\
		"fsub	v20.2d,v0.2d,v2.2d	\n\t	fadd	v14.2d,v12.2d,v14.2d\n\t"\
		"fsub	v21.2d,v1.2d,v3.2d	\n\t	fadd	v15.2d,v13.2d,v15.2d\n\t"\
		"fsub	v18.2d,v4.2d,v6.2d	\n\t	fsub	v8.2d ,v22.2d,v14.2d\n\t"\
		"fsub	v19.2d,v5.2d,v7.2d	\n\t	fsub	v9.2d ,v23.2d,v15.2d\n\t"\
		"fadd	v6.2d,v4.2d,v6.2d	\n\t	fadd	v14.2d,v22.2d,v14.2d\n\t"\
		"fadd	v7.2d,v5.2d,v7.2d	\n\t	fadd	v15.2d,v23.2d,v15.2d\n\t"\
		"fsub	v0.2d,v16.2d,v6.2d	\n\t	fsub	v10.2d,v26.2d,v25.2d\n\t"\
		"fsub	v1.2d,v17.2d,v7.2d	\n\t	fsub	v11.2d,v27.2d,v24.2d\n\t"\
		"fadd	v6.2d,v16.2d,v6.2d	\n\t	fadd	v25.2d,v26.2d,v25.2d\n\t"\
		"fadd	v7.2d,v17.2d,v7.2d	\n\t	fadd	v24.2d,v27.2d,v24.2d\n\t"\
		"fsub	v2.2d,v20.2d,v19.2d	\n\t	fadd	v12.2d,v10.2d,v24.2d\n\t"\
		"fsub	v3.2d,v21.2d,v18.2d	\n\t	fsub	v13.2d,v25.2d,v11.2d\n\t"\
		"fadd	v5.2d,v20.2d,v19.2d	\n\t	fsub	v10.2d,v10.2d,v24.2d\n\t"\
		"fadd	v4.2d,v21.2d,v18.2d	\n\t	fadd	v11.2d,v25.2d,v11.2d\n\t"\
		"fsub	v16.2d,v7.2d,v15.2d	\n\t	fmul	v12.2d,v12.2d,v29.2d\n\t"\
		"fsub	v17.2d,v1.2d,v8.2d	\n\t	fmul	v13.2d,v13.2d,v29.2d\n\t"\
		"fsub	v18.2d,v6.2d,v14.2d	\n\t	fmul	v10.2d,v10.2d,v29.2d\n\t"\
		"fsub	v19.2d,v0.2d,v9.2d	\n\t	fmul	v11.2d,v11.2d,v29.2d\n\t"\
		"fadd	v15.2d,v7.2d,v15.2d	\n\t	fsub	v22.2d,v4.2d,v12.2d	\n\t"\
		"fadd	v8.2d ,v1.2d,v8.2d	\n\t	fsub	v23.2d,v3.2d,v13.2d	\n\t"\
		"fadd	v14.2d,v6.2d,v14.2d	\n\t	fsub	v24.2d,v2.2d,v10.2d	\n\t"\
		"fadd	v9.2d ,v0.2d,v9.2d	\n\t	fsub	v25.2d,v5.2d,v11.2d	\n\t"\
		"stp	q18,q16,[x2,#0x80]	\n\t	fadd	v12.2d,v4.2d,v12.2d	\n\t"\
		"stp	q19,q8 ,[x2,#0x40]	\n\t	fadd	v13.2d,v3.2d,v13.2d	\n\t"\
		"stp	q14,q15,[x2      ]	\n\t	fadd	v10.2d,v2.2d,v10.2d	\n\t"\
		"stp	q9 ,q17,[x2,#0xc0]	\n\t	fadd	v11.2d,v5.2d,v11.2d	\n\t"\
		"									stp	q24,q22,[x2,#0xa0]	\n\t"\
		"									stp	q25,q13,[x2,#0x60]	\n\t"\
		"									stp	q10,q12,[x2,#0x20]	\n\t"\
		"									stp	q11,q23,[x2,#0xe0]	\n\t"\
		/************************************************************************************************************/\
		/* Next 2 blocks operate on odd-indexed elements from the unpck*pd commands which we stored to temporaries: */\
		/************************************************************************************************************/\
		/* Block 2: */\
		"add	x2, x2,#0x100		\n\t	ldp	q26,q27,[x2,#0x080]		\n\t"\
		"add	x3, x3,#0x100		\n\t	ldp	q22,q23,[x3,#0x80]		\n\t"\
		"ldp	q20,q21,[x2       ]	\n\t	fmul	v8.2d,v26.2d,v22.2d	\n\t"\
/* c01 */"ldp	q16,q17,[x3]		\n\t	fmul	v9.2d,v27.2d,v22.2d	\n\t"\
		"fmul	v0.2d,v20.2d,v16.2d	\n\t	fmls	v8.2d,v27.2d,v23.2d	\n\t"\
		"fmul	v1.2d,v21.2d,v16.2d	\n\t	fmla	v9.2d,v26.2d,v23.2d	\n\t"\
		"fmls	v0.2d,v21.2d,v17.2d	\n\t	ldp	q26,q27,[x2,#0x0a0]		\n\t"\
		"fmla	v1.2d,v20.2d,v17.2d	\n\t	ldp	q22,q23,[x3,#0xa0]		\n\t"\
		"ldp	q20,q21,[x2,#0x020]	\n\t	fmul	v10.2d,v26.2d,v22.2d\n\t"\
		"ldp	q16,q17,[x3,#0x20]	\n\t	fmul	v11.2d,v27.2d,v22.2d\n\t"\
		"fmul	v2.2d,v20.2d,v16.2d	\n\t	fmls	v10.2d,v27.2d,v23.2d\n\t"\
		"fmul	v3.2d,v21.2d,v16.2d	\n\t	fmla	v11.2d,v26.2d,v23.2d\n\t"\
		"fmls	v2.2d,v21.2d,v17.2d	\n\t	ldp	q26,q27,[x2,#0x0c0]		\n\t"\
		"fmla	v3.2d,v20.2d,v17.2d	\n\t	ldp	q22,q23,[x3,#0xc0]		\n\t"\
		"ldp	q20,q21,[x2,#0x040]	\n\t	fmul	v12.2d,v26.2d,v22.2d\n\t"\
		"ldp	q16,q17,[x3,#0x40]	\n\t	fmul	v13.2d,v27.2d,v22.2d\n\t"\
		"fmul	v4.2d,v20.2d,v16.2d	\n\t	fmls	v12.2d,v27.2d,v23.2d\n\t"\
		"fmul	v5.2d,v21.2d,v16.2d	\n\t	fmla	v13.2d,v26.2d,v23.2d\n\t"\
		"fmls	v4.2d,v21.2d,v17.2d	\n\t	ldp	q26,q27,[x2,#0x0e0]		\n\t"\
		"fmla	v5.2d,v20.2d,v17.2d	\n\t	ldp	q22,q23,[x3,#0xe0]		\n\t"\
		"ldp	q20,q21,[x2,#0x060]	\n\t	fmul	v14.2d,v26.2d,v22.2d\n\t"\
		"ldp	q16,q17,[x3,#0x60]	\n\t	fmul	v15.2d,v27.2d,v22.2d\n\t"\
		"fmul	v6.2d,v20.2d,v16.2d	\n\t	fmls	v14.2d,v27.2d,v23.2d\n\t"\
		"fmul	v7.2d,v21.2d,v16.2d	\n\t	fmla	v15.2d,v26.2d,v23.2d\n\t"\
		"fmls	v6.2d,v21.2d,v17.2d	\n\t	fadd	v22.2d,v8.2d,v10.2d	\n\t"/* 2x2 butterflies: */\
		"fmla	v7.2d,v20.2d,v17.2d	\n\t	fadd	v23.2d,v9.2d,v11.2d	\n\t"\
		"									fsub	v26.2d,v8.2d,v10.2d	\n\t"\
		"									fsub	v27.2d,v9.2d,v11.2d	\n\t"\
		"fadd	v16.2d,v0.2d,v2.2d	\n\t	fsub	v24.2d,v12.2d,v14.2d\n\t"\
		"fadd	v17.2d,v1.2d,v3.2d	\n\t	fsub	v25.2d,v13.2d,v15.2d\n\t"\
		"fsub	v20.2d,v0.2d,v2.2d	\n\t	fadd	v14.2d,v12.2d,v14.2d\n\t"\
		"fsub	v21.2d,v1.2d,v3.2d	\n\t	fadd	v15.2d,v13.2d,v15.2d\n\t"\
		"fsub	v18.2d,v4.2d,v6.2d	\n\t	fsub	v8.2d ,v22.2d,v14.2d\n\t"/* combine to get 2 length-4 output subtransforms... */\
		"fsub	v19.2d,v5.2d,v7.2d	\n\t	fsub	v9.2d ,v23.2d,v15.2d\n\t"\
		"fadd	v6.2d,v4.2d,v6.2d	\n\t	fadd	v14.2d,v22.2d,v14.2d\n\t"\
		"fadd	v7.2d,v5.2d,v7.2d	\n\t	fadd	v15.2d,v23.2d,v15.2d\n\t"\
		"fsub	v0.2d,v16.2d,v6.2d	\n\t	fsub	v10.2d,v26.2d,v25.2d\n\t"\
		"fsub	v1.2d,v17.2d,v7.2d	\n\t	fsub	v11.2d,v27.2d,v24.2d\n\t"\
		"fadd	v6.2d,v16.2d,v6.2d	\n\t	fadd	v25.2d,v26.2d,v25.2d\n\t"\
		"fadd	v7.2d,v17.2d,v7.2d	\n\t	fadd	v24.2d,v27.2d,v24.2d\n\t"\
		"fsub	v2.2d,v20.2d,v19.2d	\n\t	fadd	v12.2d,v10.2d,v24.2d\n\t"\
		"fsub	v3.2d,v21.2d,v18.2d	\n\t	fsub	v13.2d,v25.2d,v11.2d\n\t"\
		"fadd	v5.2d,v20.2d,v19.2d	\n\t	fsub	v10.2d,v10.2d,v24.2d\n\t"\
		"fadd	v4.2d,v21.2d,v18.2d	\n\t	fadd	v11.2d,v25.2d,v11.2d\n\t"\
		"fsub	v16.2d,v7.2d,v15.2d	\n\t	fmul	v12.2d,v12.2d,v29.2d\n\t"\
		"fsub	v17.2d,v1.2d,v8.2d	\n\t	fmul	v13.2d,v13.2d,v29.2d\n\t"\
		"fsub	v18.2d,v6.2d,v14.2d	\n\t	fmul	v10.2d,v10.2d,v29.2d\n\t"\
		"fsub	v19.2d,v0.2d,v9.2d	\n\t	fmul	v11.2d,v11.2d,v29.2d\n\t"\
		"fadd	v15.2d,v7.2d,v15.2d	\n\t	fsub	v22.2d,v4.2d,v12.2d	\n\t"\
		"fadd	v8.2d ,v1.2d,v8.2d	\n\t	fsub	v23.2d,v3.2d,v13.2d	\n\t"\
		"fadd	v14.2d,v6.2d,v14.2d	\n\t	fsub	v24.2d,v2.2d,v10.2d	\n\t"\
		"fadd	v9.2d ,v0.2d,v9.2d	\n\t	fsub	v25.2d,v5.2d,v11.2d	\n\t"\
		"stp	q18,q16,[x2,#0x80]	\n\t	fadd	v12.2d,v4.2d,v12.2d	\n\t"\
		"stp	q19,q8 ,[x2,#0x40]	\n\t	fadd	v13.2d,v3.2d,v13.2d	\n\t"\
		"stp	q14,q15,[x2      ]	\n\t	fadd	v10.2d,v2.2d,v10.2d	\n\t"\
		"stp	q9 ,q17,[x2,#0xc0]	\n\t	fadd	v11.2d,v5.2d,v11.2d	\n\t"\
		"									stp	q24,q22,[x2,#0xa0]	\n\t"\
		"									stp	q25,q13,[x2,#0x60]	\n\t"\
		"									stp	q10,q12,[x2,#0x20]	\n\t"\
		"									stp	q11,q23,[x2,#0xe0]	\n\t"\
		/* Block 3: */\
		"add	x2, x2,#0x100		\n\t	ldp	q26,q27,[x2,#0x080]		\n\t"\
		"add	x3, x3,#0x100		\n\t	ldp	q22,q23,[x3,#0x80]		\n\t"\
		"ldp	q20,q21,[x2       ]	\n\t	fmul	v8.2d,v26.2d,v22.2d	\n\t"\
/* c03 */"ldp	q16,q17,[x3]		\n\t	fmul	v9.2d,v27.2d,v22.2d	\n\t"\
		"fmul	v0.2d,v20.2d,v16.2d	\n\t	fmls	v8.2d,v27.2d,v23.2d	\n\t"\
		"fmul	v1.2d,v21.2d,v16.2d	\n\t	fmla	v9.2d,v26.2d,v23.2d	\n\t"\
		"fmls	v0.2d,v21.2d,v17.2d	\n\t	ldp	q26,q27,[x2,#0x0a0]		\n\t"\
		"fmla	v1.2d,v20.2d,v17.2d	\n\t	ldp	q22,q23,[x3,#0xa0]		\n\t"\
		"ldp	q20,q21,[x2,#0x020]	\n\t	fmul	v10.2d,v26.2d,v22.2d\n\t"\
		"ldp	q16,q17,[x3,#0x20]	\n\t	fmul	v11.2d,v27.2d,v22.2d\n\t"\
		"fmul	v2.2d,v20.2d,v16.2d	\n\t	fmls	v10.2d,v27.2d,v23.2d\n\t"\
		"fmul	v3.2d,v21.2d,v16.2d	\n\t	fmla	v11.2d,v26.2d,v23.2d\n\t"\
		"fmls	v2.2d,v21.2d,v17.2d	\n\t	ldp	q26,q27,[x2,#0x0c0]		\n\t"\
		"fmla	v3.2d,v20.2d,v17.2d	\n\t	ldp	q22,q23,[x3,#0xc0]		\n\t"\
		"ldp	q20,q21,[x2,#0x040]	\n\t	fmul	v12.2d,v26.2d,v22.2d\n\t"\
		"ldp	q16,q17,[x3,#0x40]	\n\t	fmul	v13.2d,v27.2d,v22.2d\n\t"\
		"fmul	v4.2d,v20.2d,v16.2d	\n\t	fmls	v12.2d,v27.2d,v23.2d\n\t"\
		"fmul	v5.2d,v21.2d,v16.2d	\n\t	fmla	v13.2d,v26.2d,v23.2d\n\t"\
		"fmls	v4.2d,v21.2d,v17.2d	\n\t	ldp	q26,q27,[x2,#0x0e0]		\n\t"\
		"fmla	v5.2d,v20.2d,v17.2d	\n\t	ldp	q22,q23,[x3,#0xe0]		\n\t"\
		"ldp	q20,q21,[x2,#0x060]	\n\t	fmul	v14.2d,v26.2d,v22.2d\n\t"\
		"ldp	q16,q17,[x3,#0x60]	\n\t	fmul	v15.2d,v27.2d,v22.2d\n\t"\
		"fmul	v6.2d,v20.2d,v16.2d	\n\t	fmls	v14.2d,v27.2d,v23.2d\n\t"\
		"fmul	v7.2d,v21.2d,v16.2d	\n\t	fmla	v15.2d,v26.2d,v23.2d\n\t"\
		"fmls	v6.2d,v21.2d,v17.2d	\n\t	fadd	v22.2d,v8.2d,v10.2d	\n\t"/* 2x2 butterflies: */\
		"fmla	v7.2d,v20.2d,v17.2d	\n\t	fadd	v23.2d,v9.2d,v11.2d	\n\t"\
		"									fsub	v26.2d,v8.2d,v10.2d	\n\t"\
		"									fsub	v27.2d,v9.2d,v11.2d	\n\t"\
		"fadd	v16.2d,v0.2d,v2.2d	\n\t	fsub	v24.2d,v12.2d,v14.2d\n\t"\
		"fadd	v17.2d,v1.2d,v3.2d	\n\t	fsub	v25.2d,v13.2d,v15.2d\n\t"\
		"fsub	v20.2d,v0.2d,v2.2d	\n\t	fadd	v14.2d,v12.2d,v14.2d\n\t"\
		"fsub	v21.2d,v1.2d,v3.2d	\n\t	fadd	v15.2d,v13.2d,v15.2d\n\t"\
		"fsub	v18.2d,v4.2d,v6.2d	\n\t	fsub	v8.2d ,v22.2d,v14.2d\n\t"/* combine to get 2 length-4 output subtransforms... */\
		"fsub	v19.2d,v5.2d,v7.2d	\n\t	fsub	v9.2d ,v23.2d,v15.2d\n\t"\
		"fadd	v6.2d,v4.2d,v6.2d	\n\t	fadd	v14.2d,v22.2d,v14.2d\n\t"\
		"fadd	v7.2d,v5.2d,v7.2d	\n\t	fadd	v15.2d,v23.2d,v15.2d\n\t"\
		"fsub	v0.2d,v16.2d,v6.2d	\n\t	fsub	v10.2d,v26.2d,v25.2d\n\t"\
		"fsub	v1.2d,v17.2d,v7.2d	\n\t	fsub	v11.2d,v27.2d,v24.2d\n\t"\
		"fadd	v6.2d,v16.2d,v6.2d	\n\t	fadd	v25.2d,v26.2d,v25.2d\n\t"\
		"fadd	v7.2d,v17.2d,v7.2d	\n\t	fadd	v24.2d,v27.2d,v24.2d\n\t"\
		"fsub	v2.2d,v20.2d,v19.2d	\n\t	fadd	v12.2d,v10.2d,v24.2d\n\t"\
		"fsub	v3.2d,v21.2d,v18.2d	\n\t	fsub	v13.2d,v25.2d,v11.2d\n\t"\
		"fadd	v5.2d,v20.2d,v19.2d	\n\t	fsub	v10.2d,v10.2d,v24.2d\n\t"\
		"fadd	v4.2d,v21.2d,v18.2d	\n\t	fadd	v11.2d,v25.2d,v11.2d\n\t"\
		"fsub	v16.2d,v7.2d,v15.2d	\n\t	fmul	v12.2d,v12.2d,v29.2d\n\t"\
		"fsub	v17.2d,v1.2d,v8.2d	\n\t	fmul	v13.2d,v13.2d,v29.2d\n\t"\
		"fsub	v18.2d,v6.2d,v14.2d	\n\t	fmul	v10.2d,v10.2d,v29.2d\n\t"\
		"fsub	v19.2d,v0.2d,v9.2d	\n\t	fmul	v11.2d,v11.2d,v29.2d\n\t"\
		"fadd	v15.2d,v7.2d,v15.2d	\n\t	fsub	v22.2d,v4.2d,v12.2d	\n\t"\
		"fadd	v8.2d ,v1.2d,v8.2d	\n\t	fsub	v23.2d,v3.2d,v13.2d	\n\t"\
		"fadd	v14.2d,v6.2d,v14.2d	\n\t	fsub	v24.2d,v2.2d,v10.2d	\n\t"\
		"fadd	v9.2d ,v0.2d,v9.2d	\n\t	fsub	v25.2d,v5.2d,v11.2d	\n\t"\
		"stp	q18,q16,[x2,#0x80]	\n\t	fadd	v12.2d,v4.2d,v12.2d	\n\t"\
		"stp	q19,q8 ,[x2,#0x40]	\n\t	fadd	v13.2d,v3.2d,v13.2d	\n\t"\
		"stp	q14,q15,[x2      ]	\n\t	fadd	v10.2d,v2.2d,v10.2d	\n\t"\
		"stp	q9 ,q17,[x2,#0xc0]	\n\t	fadd	v11.2d,v5.2d,v11.2d	\n\t"\
		"									stp	q24,q22,[x2,#0xa0]	\n\t"\
		"									stp	q25,q13,[x2,#0x60]	\n\t"\
		"									stp	q10,q12,[x2,#0x20]	\n\t"\
		"									stp	q11,q23,[x2,#0xe0]	\n\t"\
	/**********************************************************************************/\
	/*...and now do eight radix-4 transforms, including the internal twiddle factors: */\
	/**********************************************************************************/\
	/*...Block 1: t00,t10,t20,t30	*/		/*...Block 5: t08,t18,t28,t38	*/\
		"ldr	x2,%[__r00]				\n\t	ldr	x3,%[__cc0]			\n\t"\
		"ldp	q0,q1,[x2       ]		\n\t	ldp	q8 ,q9 ,[x2,#0x080]	\n\t"\
		"ldp	q2,q3,[x2,#0x100]		\n\t	ldp	q10,q11,[x2,#0x180]	\n\t"\
		"ldp	q4,q5,[x2,#0x200]		\n\t	ldp	q12,q13,[x2,#0x280]	\n\t"\
		"ldp	q6,q7,[x2,#0x300]		\n\t	ldp	q14,q15,[x2,#0x380]	\n\t"\
		"fadd	v16.2d,v0.2d,v2.2d		\n\t	fmul	v12.2d,v12.2d,v29.2d	\n\t"\
		"fadd	v17.2d,v1.2d,v3.2d		\n\t	fmul	v13.2d,v13.2d,v29.2d	\n\t"\
		"fadd	v18.2d,v4.2d,v6.2d		\n\t	fmul	v14.2d,v14.2d,v29.2d	\n\t"\
		"fadd	v19.2d,v5.2d,v7.2d		\n\t	fmul	v15.2d,v15.2d,v29.2d	\n\t"\
		"fsub	v0.2d ,v0.2d,v2.2d		\n\t	fadd	v20.2d,v8.2d ,v11.2d	\n\t"\
		"fsub	v1.2d ,v1.2d,v3.2d		\n\t	fadd	v21.2d,v12.2d,v13.2d	\n\t"\
		"fsub	v4.2d ,v4.2d,v6.2d		\n\t	fadd	v22.2d,v9.2d ,v10.2d	\n\t"\
		"fsub	v5.2d ,v5.2d,v7.2d		\n\t	fadd	v23.2d,v15.2d,v14.2d	\n\t"\
		"fsub	v2.2d,v16.2d,v18.2d		\n\t	fsub	v8.2d ,v8.2d ,v11.2d	\n\t"\
		"fsub	v3.2d,v17.2d,v19.2d		\n\t	fsub	v12.2d,v12.2d,v13.2d	\n\t"\
		"fadd	v6.2d,v16.2d,v18.2d		\n\t	fsub	v9.2d ,v9.2d ,v10.2d	\n\t"\
		"fadd	v7.2d,v17.2d,v19.2d		\n\t	fsub	v15.2d,v15.2d,v14.2d	\n\t"\
		"fsub	v16.2d,v0.2d,v5.2d		\n\t	fadd	v14.2d,v12.2d,v23.2d	\n\t"\
		"fsub	v17.2d,v1.2d,v4.2d		\n\t	fsub	v12.2d,v12.2d,v23.2d	\n\t"\
		"fadd	v5.2d ,v0.2d,v5.2d		\n\t	fsub	v13.2d,v21.2d,v15.2d	\n\t"\
		"fadd	v4.2d ,v1.2d,v4.2d		\n\t	fadd	v15.2d,v21.2d,v15.2d	\n\t"\
		/* Wrapper_square 32-DIF outputs to tmp-array rather than main-array, so needs middle 2 of each output quartet swapped: */\
		/* Lcol Outputs-to-main-array offset order is 0,1,2,3, so middle-2-swap (1<->2) changes that to 0,2,1,3: */\
		"stp	q6 ,q7 ,[x2       ]		\n\t	fsub	v21.2d,v8.2d ,v12.2d	\n\t"\
		"stp	q2, q3 ,[x2,#0x200]		\n\t	fsub	v10.2d,v22.2d,v13.2d	\n\t"\
		"stp	q16,q4 ,[x2,#0x100]		\n\t	fadd	v12.2d,v8.2d ,v12.2d	\n\t"\
		"stp	q5 ,q17,[x2,#0x300]		\n\t	fadd	v13.2d,v22.2d,v13.2d	\n\t"\
		/* Rcol Outputs-to-main-array offset order is 4,5,7,6, so middle-2-swap (5<->6) changes that to 4,6,7,5, and addresses now interleaved with Lcol ones rather than consecutive: */\
		"	fsub	v11.2d,v20.2d,v15.2d	\n\t	stp	q12,q13,[x2,#0x080]	\n\t"\
		"	fsub	v22.2d,v9.2d ,v14.2d	\n\t	stp	q21,q10,[x2,#0x280]	\n\t"\
		"	fadd	v15.2d,v20.2d,v15.2d	\n\t	stp	q15,q22,[x2,#0x380]	\n\t"\
		"	fadd	v14.2d,v9.2d ,v14.2d	\n\t	stp	q11,q14,[x2,#0x180]	\n\t"\
	/*...Block 3: t04,t14,t24,t34; outputs into add0 + p8-11: */\
		"add	x2,x2,#0x40	\n\t"/* r04 */"	ldp	q16,q17,[x10]			\n\t"/* cc0 */\
											/*...Block 7: t0C,t1C,t2C,t3C; outputs into add0 + p12-15: */\
		"ldp	q0 ,q1 ,[x2       ]		\n\t"\
		"ldp	q2 ,q3 ,[x2,#0x100]		\n\t"\
		"ldp	q4 ,q5 ,[x2,#0x200]		\n\t"\
		"ldp	q20,q21,[x2,#0x300]		\n\t"\
		"fmul	v18.2d,v4.2d,v16.2d		\n\t	ldp	q8 ,q9 ,[x2,#0x080]	\n\t"\
		"fmul	v19.2d,v5.2d,v16.2d		\n\t	ldp	q10,q11,[x2,#0x180]	\n\t"\
		"fmls	v18.2d,v5.2d,v17.2d		\n\t	ldp	q12,q13,[x2,#0x280]	\n\t"\
		"fmla	v19.2d,v4.2d,v17.2d		\n\t	ldp	q24,q25,[x2,#0x380]	\n\t"\
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
		/* Lcol Outputs-to-main-array offset order is 0,1,2,3, so middle-2-swap (1<->2) changes that to 0,2,1,3: */\
		"stp	q18,q19,[x2       ]		\n\t	fsub	v22.2d,v8.2d ,v12.2d	\n\t"/*  8= 8-12 */\
		"stp	q2 ,q3 ,[x2,#0x200]		\n\t	fsub	v23.2d,v9.2d ,v13.2d	\n\t"/*  9= 9-13 */\
		"stp	q20,q21,[x2,#0x100]		\n\t	fadd	v12.2d,v8.2d ,v12.2d	\n\t"/* 12= 8+12 */\
		"stp	q5 ,q1 ,[x2,#0x300]		\n\t	fadd	v13.2d,v9.2d ,v13.2d	\n\t"/* 13= 9+13 */\
		/* Rcol Outputs-to-main-array offset order is 6,7,5,4, so middle-2-swap (5<->6) changes that to 5,7,6,4: */\
		"	stp	q10,q24,[x2,#0x180]	\n\t"\
		"	stp	q25,q11,[x2,#0x380]	\n\t"\
		"	stp	q22,q23,[x2,#0x280]	\n\t"\
		"	stp	q12,q13,[x2,#0x080]	\n\t"\
	/*...Block 2: t02,t12,t22,t32: */		/*...Block 6: t0A,t1A,t2A,t3A; outputs into add0 + p20-23: */\
		"sub	x2,x2,#0x20				\n\t"/* r02 */\
		"ldp	q16,q17,[x10,#0x20]		\n\t	ldp	q14,q15,[x10,#0x40]		\n\t"/* cc1,cc3 */\
		"ldp	q0 ,q1 ,[x2       ]		\n\t"\
		"ldp	q2 ,q3 ,[x2,#0x100]		\n\t"\
		"ldp	q4 ,q5 ,[x2,#0x200]		\n\t"\
		"ldp	q20,q21,[x2,#0x300]		\n\t"\
		"fmul	v18.2d,v4.2d,v16.2d		\n\t	ldp	q8 ,q9 ,[x2,#0x080]	\n\t"\
		"fmul	v19.2d,v5.2d,v16.2d		\n\t	ldp	q10,q11,[x2,#0x180]	\n\t"\
		"fmls	v18.2d,v5.2d,v17.2d		\n\t	ldp	q12,q13,[x2,#0x280]	\n\t"\
		"fmla	v19.2d,v4.2d,v17.2d		\n\t	ldp	q24,q25,[x2,#0x380]	\n\t"\
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
		/* Lcol Outputs-to-main-array offset order is 0,1,2,3, so middle-2-swap (1<->2) changes that to 0,2,1,3: */\
		"stp	q18,q19,[x2       ]		\n\t	fsub	v22.2d,v8.2d ,v12.2d	\n\t"\
		"stp	q2 ,q3 ,[x2,#0x200]		\n\t	fsub	v23.2d,v9.2d ,v13.2d	\n\t"\
		"stp	q20,q21,[x2,#0x100]		\n\t	fadd	v12.2d,v8.2d ,v12.2d	\n\t"\
		"stp	q5 ,q1 ,[x2,#0x300]		\n\t	fadd	v13.2d,v9.2d ,v13.2d	\n\t"\
		/* Rcol Outputs-to-main-array offset order is 6,7,5,4, so middle-2-swap (5<->6) changes that to 5,7,6,4: */\
		"	stp	q10,q24,[x2,#0x180]	\n\t"\
		"	stp	q25,q11,[x2,#0x380]	\n\t"\
		"	stp	q22,q23,[x2,#0x280]	\n\t"\
		"	stp	q12,q13,[x2,#0x080]	\n\t"\
	/*...Block 4: t06,t16,t26,t36: */		/*...Block 8: t0E,t1E,t2E,t3E: */\
		"add	x2,x2,#0x40				\n\t"/* r06 */\
		"ldp	q16,q17,[x10,#0x40]		\n\t	ldp	q14,q15,[x10,#0x20]		\n\t"/* cc3,cc1 */\
		"ldp	q0 ,q1 ,[x2       ]		\n\t"\
		"ldp	q2 ,q3 ,[x2,#0x100]		\n\t"\
		"ldp	q4 ,q5 ,[x2,#0x200]		\n\t"\
		"ldp	q20,q21,[x2,#0x300]		\n\t"\
		"fmul	v18.2d,v4.2d,v16.2d		\n\t	ldp	q8 ,q9 ,[x2,#0x080]	\n\t"\
		"fmul	v19.2d,v5.2d,v16.2d		\n\t	ldp	q10,q11,[x2,#0x180]	\n\t"\
		"fmls	v18.2d,v5.2d,v17.2d		\n\t	ldp	q12,q13,[x2,#0x280]	\n\t"\
		"fmla	v19.2d,v4.2d,v17.2d		\n\t	ldp	q24,q25,[x2,#0x380]	\n\t"\
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
		/* Lcol Outputs-to-main-array offset order is 0,1,2,3, so middle-2-swap (1<->2) changes that to 0,2,1,3: */\
		"stp	q18,q19,[x2       ]		\n\t	fsub	v22.2d,v8.2d ,v12.2d	\n\t"\
		"stp	q2 ,q3 ,[x2,#0x200]		\n\t	fsub	v23.2d,v9.2d ,v13.2d	\n\t"\
		"stp	q20,q21,[x2,#0x100]		\n\t	fadd	v12.2d,v8.2d ,v12.2d	\n\t"\
		"stp	q5 ,q1 ,[x2,#0x300]		\n\t	fadd	v13.2d,v9.2d ,v13.2d	\n\t"\
		/* Rcol Outputs-to-main-array offset order is 6,7,5,4, so middle-2-swap (5<->6) changes that to 5,7,6,4: */\
		"	stp	q10,q24,[x2,#0x180]	\n\t"\
		"	stp	q25,q11,[x2,#0x380]	\n\t"\
		"	stp	q22,q23,[x2,#0x280]	\n\t"\
		"	stp	q12,q13,[x2,#0x080]	\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__add1] "m" (Xadd1)\
		 ,[__r00] "m" (Xr00)\
		 ,[__r10] "m" (Xr10)\
		 ,[__r20] "m" (Xr20)\
		 ,[__r30] "m" (Xr30)\
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__cc0] "m" (Xcc0)\
		 ,[__c00] "m" (Xc00)\
		 ,[__c01] "m" (Xc01)\
		 ,[__c02] "m" (Xc02)\
		 ,[__c03] "m" (Xc03)\
		 ,[__c05] "m" (Xc05)\
		 ,[__c07] "m" (Xc07)\
		: "cc","memory","x0","x1","x2","x3","x10","x11",\
		"v0","v1","v2","v3","v4","v5","v6","v7","v8","v9","v10","v11","v12","v13","v14","v15",\
		"v16","v17","v18","v19","v20","v21","v22","v23","v24","v25","v26","v27", "v29"	/* Clobbered registers */\
	);\
	}

	#define SSE2_RADIX32_WRAPPER_DIT(Xadd0,Xadd1,Xisrt2,Xr00,Xr08,Xr10,Xc01,Xc02,Xc04,Xc06,Xc08,Xc0A,Xc0C,Xc0E,Xc10,Xc12,Xc14,Xc16,Xc18,Xc1A,Xc1C,Xc1E)\
	{\
	__asm__ volatile (\
		"ldr	x0,%[__add0]		\n\t"\
		"ldr	x1,%[__add1]		\n\t"\
	/************************************************************************/\
	/*...And do an inverse DIT radix-32 pass on the squared-data blocks:	*/\
	/************************************************************************/\
		"ldr	x2,%[__r00]			\n\t	ldr w11,=0x10				\n\t"/* NB: ~0, not -: */\
		"add	x10, x2,#0x480		\n\t	ld1r	{v29.2d},[x10],x11	\n\t"/* x10 has isrt2-ptr for LD1R, gets post-incr'ed by 0x10 to cc0 */\
	/* Block 0:							rcol vs lcol: Input-addresses += 0x80: */\
		"ldp	q0,q1,[x2       ]	\n\t	ldp	q12,q13,[x2,#0x080]		\n\t"\
		"ldp	q8,q9,[x2,#0x200]	\n\t	ldp	q20,q21,[x2,#0x280]		\n\t"\
		"fsub	v2.2d,v0.2d,v8.2d	\n\t	fsub	v14.2d,v12.2d,v20.2d	\n\t"\
		"fsub	v3.2d,v1.2d,v9.2d	\n\t	fsub	v15.2d,v13.2d,v21.2d	\n\t"\
		"fadd	v0.2d,v0.2d,v8.2d	\n\t	fadd	v12.2d,v12.2d,v20.2d	\n\t"\
		"fadd	v1.2d,v1.2d,v9.2d	\n\t	fadd	v13.2d,v13.2d,v21.2d	\n\t"\
		"ldp	q4,q5,[x2,#0x100]	\n\t	ldp	q16,q17,[x2,#0x180]		\n\t"\
		"ldp	q8,q9,[x2,#0x300]	\n\t	ldp	q20,q21,[x2,#0x380]		\n\t"\
		"fsub	v6.2d,v4.2d,v8.2d	\n\t	fsub	v18.2d,v16.2d,v20.2d	\n\t"\
		"fsub	v7.2d,v5.2d,v9.2d	\n\t	fsub	v19.2d,v17.2d,v21.2d	\n\t"\
		"fadd	v4.2d,v4.2d,v8.2d	\n\t	fadd	v16.2d,v16.2d,v20.2d	\n\t"\
		"fadd	v5.2d,v5.2d,v9.2d	\n\t	fadd	v17.2d,v17.2d,v21.2d	\n\t"\
		/* combine to get the 2 length-4 transforms: */\
		"fsub	v8.2d,v0.2d,v4.2d	\n\t	fsub	v20.2d,v12.2d,v16.2d	\n\t"/* t5,13 */\
		"fsub	v9.2d,v1.2d,v5.2d	\n\t	fsub	v21.2d,v13.2d,v17.2d	\n\t"/* t6,14 */\
		"fadd	v0.2d,v0.2d,v4.2d	\n\t	fadd	v12.2d,v12.2d,v16.2d	\n\t"/* t1, 9 */\
		"fadd	v1.2d,v1.2d,v5.2d	\n\t	fadd	v13.2d,v13.2d,v17.2d	\n\t"/* t2,10 */\
		"fadd	v4.2d,v2.2d,v7.2d	\n\t	fadd	v16.2d,v14.2d,v19.2d	\n\t"/* t3,11 */\
		"fsub	v5.2d,v3.2d,v6.2d	\n\t	fsub	v17.2d,v15.2d,v18.2d	\n\t"/* t4,12 */\
		"fsub	v2.2d,v2.2d,v7.2d	\n\t	fsub	v14.2d,v14.2d,v19.2d	\n\t"/* t7,15 */\
		"fadd	v3.2d,v3.2d,v6.2d	\n\t	fadd	v15.2d,v15.2d,v18.2d	\n\t"/* t8,16 */\
		/* now combine the two half-transforms: */\
		"fadd	v6.2d,v0.2d,v12.2d	\n\t	fadd	v18.2d,v1.2d,v13.2d	\n\t"\
		"fsub	v7.2d,v0.2d,v12.2d	\n\t	fsub	v19.2d,v1.2d,v13.2d	\n\t"\
		"stp	q6,q18,[x2       ]	\n\t	stp	q7,q19,[x2,#0x080]		\n\t"\
		"fadd	v6.2d,v16.2d,v17.2d	\n\t	fsub	v7.2d,v16.2d,v17.2d	\n\t"/* t3,4 = v4,5; t11,12 = v16,17 */\
		"fmul	v16.2d,v29.2d,v6.2d	\n\t	fmul	v17.2d,v29.2d,v7.2d	\n\t"/* rt,it; */\
		"fadd	v6.2d,v4.2d,v16.2d	\n\t	fsub	v18.2d,v5.2d,v17.2d	\n\t"/* t11,12 = v6,18 */\
		"fsub	v7.2d,v4.2d,v16.2d	\n\t	fadd	v19.2d,v5.2d,v17.2d	\n\t"/* t3 ,4  = v7,19 */\
		"stp	q6,q18,[x2,#0x100]	\n\t	stp	q7,q19,[x2,#0x180]		\n\t"\
		"fadd	v6.2d,v8.2d,v21.2d	\n\t	fsub	v18.2d,v9.2d,v20.2d	\n\t"/* rt,it = v6,18 */\
		"fsub	v7.2d,v8.2d,v21.2d	\n\t	fadd	v19.2d,v9.2d,v20.2d	\n\t"/* t3,4  = v7,19 */\
		"stp	q6,q18,[x2,#0x200]	\n\t	stp	q7,q19,[x2,#0x280]		\n\t"\
		"fsub	v6.2d,v14.2d,v15.2d	\n\t	fadd	v7.2d,v14.2d,v15.2d	\n\t"/* t7,8 = v2,3; t15,16 = v14,15 */\
		"fmul	v14.2d,v29.2d,v6.2d	\n\t	fmul	v15.2d,v29.2d,v7.2d	\n\t"/* rt,it; */\
		"fsub	v6.2d,v2.2d,v14.2d	\n\t	fsub	v18.2d,v3.2d,v15.2d	\n\t"/* t15,16 = v6,18 */\
		"fadd	v7.2d,v2.2d,v14.2d	\n\t	fadd	v19.2d,v3.2d,v15.2d	\n\t"/* t7 ,8  = v7,19 */\
		"stp	q6,q18,[x2,#0x300]	\n\t	stp	q7,q19,[x2,#0x380]		\n\t"\
	/* Block 1: */\
		"add	x2, x2,#0x20	\n\t"\
		"ldp	q0,q1,[x2       ]	\n\t	ldp	q12,q13,[x2,#0x080]		\n\t"\
		"ldp	q8,q9,[x2,#0x200]	\n\t	ldp	q20,q21,[x2,#0x280]		\n\t"\
		"fsub	v2.2d,v0.2d,v8.2d	\n\t	fsub	v14.2d,v12.2d,v20.2d	\n\t"\
		"fsub	v3.2d,v1.2d,v9.2d	\n\t	fsub	v15.2d,v13.2d,v21.2d	\n\t"\
		"fadd	v0.2d,v0.2d,v8.2d	\n\t	fadd	v12.2d,v12.2d,v20.2d	\n\t"\
		"fadd	v1.2d,v1.2d,v9.2d	\n\t	fadd	v13.2d,v13.2d,v21.2d	\n\t"\
		"ldp	q4,q5,[x2,#0x100]	\n\t	ldp	q16,q17,[x2,#0x180]		\n\t"\
		"ldp	q8,q9,[x2,#0x300]	\n\t	ldp	q20,q21,[x2,#0x380]		\n\t"\
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
		"stp	q6,q18,[x2       ]	\n\t	stp	q7,q19,[x2,#0x080]		\n\t"\
		"fadd	v6.2d,v16.2d,v17.2d	\n\t	fsub	v7.2d,v16.2d,v17.2d	\n\t"\
		"fmul	v16.2d,v29.2d,v6.2d	\n\t	fmul	v17.2d,v29.2d,v7.2d	\n\t"\
		"fadd	v6.2d,v4.2d,v16.2d	\n\t	fsub	v18.2d,v5.2d,v17.2d	\n\t"\
		"fsub	v7.2d,v4.2d,v16.2d	\n\t	fadd	v19.2d,v5.2d,v17.2d	\n\t"\
		"stp	q6,q18,[x2,#0x100]	\n\t	stp	q7,q19,[x2,#0x180]		\n\t"\
		"fadd	v6.2d,v8.2d,v21.2d	\n\t	fsub	v18.2d,v9.2d,v20.2d	\n\t"\
		"fsub	v7.2d,v8.2d,v21.2d	\n\t	fadd	v19.2d,v9.2d,v20.2d	\n\t"\
		"stp	q6,q18,[x2,#0x200]	\n\t	stp	q7,q19,[x2,#0x280]		\n\t"\
		"fsub	v6.2d,v14.2d,v15.2d	\n\t	fadd	v7.2d,v14.2d,v15.2d	\n\t"\
		"fmul	v14.2d,v29.2d,v6.2d	\n\t	fmul	v15.2d,v29.2d,v7.2d	\n\t"\
		"fsub	v6.2d,v2.2d,v14.2d	\n\t	fsub	v18.2d,v3.2d,v15.2d	\n\t"\
		"fadd	v7.2d,v2.2d,v14.2d	\n\t	fadd	v19.2d,v3.2d,v15.2d	\n\t"\
		"stp	q6,q18,[x2,#0x300]	\n\t	stp	q7,q19,[x2,#0x380]		\n\t"\
	/* Block 2: */\
		"add	x2, x2,#0x20	\n\t"\
		"ldp	q0,q1,[x2       ]	\n\t	ldp	q12,q13,[x2,#0x080]		\n\t"\
		"ldp	q8,q9,[x2,#0x200]	\n\t	ldp	q20,q21,[x2,#0x280]		\n\t"\
		"fsub	v2.2d,v0.2d,v8.2d	\n\t	fsub	v14.2d,v12.2d,v20.2d	\n\t"\
		"fsub	v3.2d,v1.2d,v9.2d	\n\t	fsub	v15.2d,v13.2d,v21.2d	\n\t"\
		"fadd	v0.2d,v0.2d,v8.2d	\n\t	fadd	v12.2d,v12.2d,v20.2d	\n\t"\
		"fadd	v1.2d,v1.2d,v9.2d	\n\t	fadd	v13.2d,v13.2d,v21.2d	\n\t"\
		"ldp	q4,q5,[x2,#0x100]	\n\t	ldp	q16,q17,[x2,#0x180]		\n\t"\
		"ldp	q8,q9,[x2,#0x300]	\n\t	ldp	q20,q21,[x2,#0x380]		\n\t"\
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
		"stp	q6,q18,[x2       ]	\n\t	stp	q7,q19,[x2,#0x080]		\n\t"\
		"fadd	v6.2d,v16.2d,v17.2d	\n\t	fsub	v7.2d,v16.2d,v17.2d	\n\t"\
		"fmul	v16.2d,v29.2d,v6.2d	\n\t	fmul	v17.2d,v29.2d,v7.2d	\n\t"\
		"fadd	v6.2d,v4.2d,v16.2d	\n\t	fsub	v18.2d,v5.2d,v17.2d	\n\t"\
		"fsub	v7.2d,v4.2d,v16.2d	\n\t	fadd	v19.2d,v5.2d,v17.2d	\n\t"\
		"stp	q6,q18,[x2,#0x100]	\n\t	stp	q7,q19,[x2,#0x180]		\n\t"\
		"fadd	v6.2d,v8.2d,v21.2d	\n\t	fsub	v18.2d,v9.2d,v20.2d	\n\t"\
		"fsub	v7.2d,v8.2d,v21.2d	\n\t	fadd	v19.2d,v9.2d,v20.2d	\n\t"\
		"stp	q6,q18,[x2,#0x200]	\n\t	stp	q7,q19,[x2,#0x280]		\n\t"\
		"fsub	v6.2d,v14.2d,v15.2d	\n\t	fadd	v7.2d,v14.2d,v15.2d	\n\t"\
		"fmul	v14.2d,v29.2d,v6.2d	\n\t	fmul	v15.2d,v29.2d,v7.2d	\n\t"\
		"fsub	v6.2d,v2.2d,v14.2d	\n\t	fsub	v18.2d,v3.2d,v15.2d	\n\t"\
		"fadd	v7.2d,v2.2d,v14.2d	\n\t	fadd	v19.2d,v3.2d,v15.2d	\n\t"\
		"stp	q6,q18,[x2,#0x300]	\n\t	stp	q7,q19,[x2,#0x380]		\n\t"\
	/* Block 3: */\
		"add	x2, x2,#0x20	\n\t"\
		"ldp	q0,q1,[x2       ]	\n\t	ldp	q12,q13,[x2,#0x080]		\n\t"\
		"ldp	q8,q9,[x2,#0x200]	\n\t	ldp	q20,q21,[x2,#0x280]		\n\t"\
		"fsub	v2.2d,v0.2d,v8.2d	\n\t	fsub	v14.2d,v12.2d,v20.2d	\n\t"\
		"fsub	v3.2d,v1.2d,v9.2d	\n\t	fsub	v15.2d,v13.2d,v21.2d	\n\t"\
		"fadd	v0.2d,v0.2d,v8.2d	\n\t	fadd	v12.2d,v12.2d,v20.2d	\n\t"\
		"fadd	v1.2d,v1.2d,v9.2d	\n\t	fadd	v13.2d,v13.2d,v21.2d	\n\t"\
		"ldp	q4,q5,[x2,#0x100]	\n\t	ldp	q16,q17,[x2,#0x180]		\n\t"\
		"ldp	q8,q9,[x2,#0x300]	\n\t	ldp	q20,q21,[x2,#0x380]		\n\t"\
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
		"stp	q6,q18,[x2       ]	\n\t	stp	q7,q19,[x2,#0x080]		\n\t"\
		"fadd	v6.2d,v16.2d,v17.2d	\n\t	fsub	v7.2d,v16.2d,v17.2d	\n\t"\
		"fmul	v16.2d,v29.2d,v6.2d	\n\t	fmul	v17.2d,v29.2d,v7.2d	\n\t"\
		"fadd	v6.2d,v4.2d,v16.2d	\n\t	fsub	v18.2d,v5.2d,v17.2d	\n\t"\
		"fsub	v7.2d,v4.2d,v16.2d	\n\t	fadd	v19.2d,v5.2d,v17.2d	\n\t"\
		"stp	q6,q18,[x2,#0x100]	\n\t	stp	q7,q19,[x2,#0x180]		\n\t"\
		"fadd	v6.2d,v8.2d,v21.2d	\n\t	fsub	v18.2d,v9.2d,v20.2d	\n\t"\
		"fsub	v7.2d,v8.2d,v21.2d	\n\t	fadd	v19.2d,v9.2d,v20.2d	\n\t"\
		"stp	q6,q18,[x2,#0x200]	\n\t	stp	q7,q19,[x2,#0x280]		\n\t"\
		"fsub	v6.2d,v14.2d,v15.2d	\n\t	fadd	v7.2d,v14.2d,v15.2d	\n\t"\
		"fmul	v14.2d,v29.2d,v6.2d	\n\t	fmul	v15.2d,v29.2d,v7.2d	\n\t"\
		"fsub	v6.2d,v2.2d,v14.2d	\n\t	fsub	v18.2d,v3.2d,v15.2d	\n\t"\
		"fadd	v7.2d,v2.2d,v14.2d	\n\t	fadd	v19.2d,v3.2d,v15.2d	\n\t"\
		"stp	q6,q18,[x2,#0x300]	\n\t	stp	q7,q19,[x2,#0x380]		\n\t"\
	/***************************************************************************************/\
	/* Now do eight more radix-4 transforms, including the internal and external twiddles: */\
	/***************************************************************************************/\
/***
NOTE: Twiddles here shared between DIF/DIT portions, laid out like DIF, i.e. slightly differently than in DIT routine
	of radix32_dif_dit_pass ... address offsets within each columnwise quartet go as 0x[0,40,20,60], not 0x[0,20,40,60].
***/\
		/*...Block 2: t02,t12,t22,t32	*/		/*...Block 6: t0A,t1A,t2A,t3A*/\
		/* Lcol stores into add1+0x[000,080,100,180], Rcol into add1+0x[040,0c0,140,1c0]: */\
		"ldr	x2,%[__r10]				\n\t"/* In comments, let cc1=[c,s], cc3=[a,b]: */\
		"ldp	q0,q1,[x2,#0x20]		\n\t"\
		"ldp	q2,q3,[x2,#0x60]		\n\t"\
/*cc1:*/"ldp	q18,q19,[x10,#0x20]		\n\t"\
/*cc3:*/"ldp	q20,q21,[x10,#0x40]		\n\t"\
		"fmul	v4.2d,v0.2d,v18.2d		\n\t	ldp	q8 ,q9 ,[x2,#0xa0]		\n\t"\
		"fmul	v5.2d,v1.2d,v18.2d		\n\t	ldp	q10,q11,[x2,#0xe0]		\n\t"\
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
		"ldp	q0,q1,[x2      ]		\n\t	fsub	v12.2d,v12.2d,v24.2d	\n\t"\
		"ldp	q2,q3,[x2,#0x40]		\n\t	fsub	v13.2d,v13.2d,v25.2d	\n\t"\
		/* cc0 := [c,s]: */\
		"ldp	q20,q21,[x10]			\n\t"\
		"fmul	v16.2d,v2.2d,v20.2d		\n\t	ldp	q8 ,q9 ,[x2,#0x80]		\n\t"\
		"fmul	v17.2d,v3.2d,v20.2d		\n\t	ldp	q10,q11,[x2,#0xc0]		\n\t"\
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
		"ldp	q2,q3,[x10,#0x2a0]		\n\t"	/* Let c05 = [c,s], c0D = [a,b]: */\
		"fmul	v24.2d,v4.2d,v0.2d		\n\t	ldp	q8 ,q9 ,[x10,#0x2e0]		\n\t"\
		"fmul	v25.2d,v5.2d,v0.2d		\n\t	ldp	q10,q11,[x10,#0x320]		\n\t"\
/*4c+5s*/"fmla	v24.2d,v5.2d,v1.2d		\n\t	fmul	v26.2d,v12.2d,v8.2d		\n\t"\
/*5c-4s*/"fmls	v25.2d,v4.2d,v1.2d		\n\t	fmul	v27.2d,v13.2d,v8.2d		\n\t"\
		"fmul	v4.2d,v17.2d,v2.2d		\n\t	fmla	v26.2d,v13.2d,v9.2d		\n\t"/*12c+13s*/\
		"fmul	v5.2d,v7.2d ,v2.2d		\n\t	fmls	v27.2d,v12.2d,v9.2d		\n\t"/*13c-12s*/\
/*1a-7b*/"fmls	v4.2d,v7.2d ,v3.2d		\n\t	fmul	v12.2d,v21.2d,v10.2d	\n\t"\
/*7a+1b*/"fmla	v5.2d,v17.2d,v3.2d		\n\t	fmul	v13.2d,v15.2d,v10.2d	\n\t"\
		"stp	q24,q25,[x1       ]		\n\t	fmls	v12.2d,v15.2d,v11.2d	\n\t"/* 9a-15b*/\
		"stp	q5 ,q4 ,[x1,#0x080]		\n\t	fmla	v13.2d,v21.2d,v11.2d	\n\t"/*15a+ 9b*/\
		/* Let c11 = [c,s], c19 = [a,b]: */		/* Let c15 = [c,s], c1D = [a,b]: */\
		"ldp	q0,q1,[x10,#0x280]		\n\t"\
		"ldp	q2,q3,[x10,#0x2c0]		\n\t"\
		"fmul	v24.2d,v18.2d,v0.2d		\n\t	stp	q26,q27,[x1,#0x040]		\n\t"\
		"fmul	v25.2d,v19.2d,v0.2d		\n\t	stp	q13,q12,[x1,#0x0c0]		\n\t"\
/*2c+3s*/"fmla	v24.2d,v19.2d,v1.2d		\n\t	ldp	q8 ,q9 ,[x10,#0x300]		\n\t"\
/*3c-2s*/"fmls	v25.2d,v18.2d,v1.2d		\n\t	ldp	q10,q11,[x10,#0x340]		\n\t"\
		"fmul	v4.2d,v6.2d ,v2.2d		\n\t	fmul	v26.2d,v22.2d,v8.2d		\n\t"\
		"fmul	v5.2d,v16.2d,v2.2d		\n\t	fmul	v27.2d,v23.2d,v8.2d		\n\t"\
/*6a-0b*/"fmls	v4.2d,v16.2d,v3.2d		\n\t	fmla	v26.2d,v23.2d,v9.2d		\n\t"/*10c+11s*/\
/*0a+6b*/"fmla	v5.2d,v6.2d ,v3.2d		\n\t	fmls	v27.2d,v22.2d,v9.2d		\n\t"/*11c-10s*/\
		"stp	q24,q25,[x1,#0x100]		\n\t	fmul	v12.2d,v14.2d,v10.2d	\n\t"\
		"stp	q5 ,q4 ,[x1,#0x180]		\n\t	fmul	v13.2d,v20.2d,v10.2d	\n\t"\
		"								\n\t	fmls	v12.2d,v20.2d,v11.2d	\n\t"/*14a-8b*/\
		"								\n\t	fmla	v13.2d,v14.2d,v11.2d	\n\t"/*8a+14b*/\
		"								\n\t	stp	q26,q27,[x1,#0x140]		\n\t"\
		"								\n\t	stp	q13,q12,[x1,#0x1c0]		\n\t"\
		/*...Block 4: t06,t16,t26,t36*/			/*...Block 8: t0E,t1E,t2E,t3E*/\
		"add	x2,x2,#0x200			\n\t"/* r30 */\
		"ldp	q0,q1,[x2,#0x20]		\n\t"\
		"ldp	q2,q3,[x2,#0x60]		\n\t"\
/*cc1:=[c,s]*/"ldp	q18,q19,[x10,#0x20]	\n\t"\
/*cc3:=[a,b]*/"ldp	q20,q21,[x10,#0x40]	\n\t"\
		"fmul	v4.2d,v0.2d,v20.2d		\n\t	ldp	q8 ,q9 ,[x2,#0xa0]		\n\t"\
		"fmul	v5.2d,v1.2d,v20.2d		\n\t	ldp	q10,q11,[x2,#0xe0]		\n\t"\
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
		"ldp	q0,q1,[x2      ]		\n\t	fsub	v12.2d,v12.2d,v24.2d	\n\t"\
		"ldp	q2,q3,[x2,#0x40]		\n\t	fsub	v13.2d,v13.2d,v25.2d	\n\t"\
		/* cc0 := [c,s]: */\
		"ldp	q20,q21,[x10]			\n\t"\
		"fmul	v16.2d,v2.2d,v21.2d		\n\t	ldp	q8 ,q9 ,[x2,#0x80]		\n\t"\
		"fmul	v17.2d,v3.2d,v21.2d		\n\t	ldp	q10,q11,[x2,#0xc0]		\n\t"\
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
		"ldp	q2,q3,[x10,#0x40]		\n\t"	/* Let c05 = [c,s], c0F = [a,b]: */\
		"fmul	v24.2d,v4.2d,v0.2d		\n\t	ldp	q8 ,q9 ,[x10,#0x80]		\n\t"\
		"fmul	v25.2d,v5.2d,v0.2d		\n\t	ldp	q10,q11,[x10,#0xc0]		\n\t"\
/*4c+5s*/"fmla	v24.2d,v5.2d,v1.2d		\n\t	fmul	v26.2d,v12.2d,v8.2d		\n\t"\
/*5c-4s*/"fmls	v25.2d,v4.2d,v1.2d		\n\t	fmul	v27.2d,v13.2d,v8.2d		\n\t"\
		"fmul	v4.2d,v17.2d,v2.2d		\n\t	fmla	v26.2d,v13.2d,v9.2d		\n\t"/*12c+13s*/\
		"fmul	v5.2d,v7.2d ,v2.2d		\n\t	fmls	v27.2d,v12.2d,v9.2d		\n\t"/*13c-12s*/\
/*1a-7b*/"fmls	v4.2d,v7.2d ,v3.2d		\n\t	fmul	v12.2d,v21.2d,v10.2d	\n\t"\
/*7a+1b*/"fmla	v5.2d,v17.2d,v3.2d		\n\t	fmul	v13.2d,v15.2d,v10.2d	\n\t"\
		"stp	q24,q25,[x1,#0x020]		\n\t	fmls	v12.2d,v15.2d,v11.2d	\n\t"/* 9a-15b*/\
		"stp	q5 ,q4 ,[x1,#0x0a0]		\n\t	fmla	v13.2d,v21.2d,v11.2d	\n\t"/*15a+ 9b*/\
		/* Let c13 = [c,s], c1B = [a,b]: */		/* Let c17 = [c,s], c1F = [a,b]: */\
		"ldp	q0,q1,[x10,#0x20]		\n\t"\
		"ldp	q2,q3,[x10,#0x60]		\n\t"\
		"fmul	v24.2d,v18.2d,v0.2d		\n\t	stp	q26,q27,[x1,#0x060]		\n\t"\
		"fmul	v25.2d,v19.2d,v0.2d		\n\t	stp	q13,q12,[x1,#0x0e0]		\n\t"\
/*2c+3s*/"fmla	v24.2d,v19.2d,v1.2d		\n\t	ldp	q8 ,q9 ,[x10,#0xa0]		\n\t"\
/*3c-2s*/"fmls	v25.2d,v18.2d,v1.2d		\n\t	ldp	q10,q11,[x10,#0xe0]		\n\t"\
		"fmul	v4.2d,v6.2d ,v2.2d		\n\t	fmul	v26.2d,v22.2d,v8.2d		\n\t"\
		"fmul	v5.2d,v16.2d,v2.2d		\n\t	fmul	v27.2d,v23.2d,v8.2d		\n\t"\
/*6a-0b*/"fmls	v4.2d,v16.2d,v3.2d		\n\t	fmla	v26.2d,v23.2d,v9.2d		\n\t"/*10c+11s*/\
/*0a+6b*/"fmla	v5.2d,v6.2d ,v3.2d		\n\t	fmls	v27.2d,v22.2d,v9.2d		\n\t"/*11c-10s*/\
		"stp	q24,q25,[x1,#0x120]		\n\t	fmul	v12.2d,v14.2d,v10.2d	\n\t"\
		"stp	q5 ,q4 ,[x1,#0x1a0]		\n\t	fmul	v13.2d,v20.2d,v10.2d	\n\t"\
		"										fmls	v12.2d,v20.2d,v11.2d	\n\t"/*14a-8b*/\
		"										fmla	v13.2d,v14.2d,v11.2d	\n\t"/*8a+14b*/\
		"										stp	q26,q27,[x1,#0x160]		\n\t"\
		"										stp	q13,q12,[x1,#0x1e0]		\n\t"\
		"sub	x10,x10,#0x360	\n\t"\
		/*...Block 1: t00,t10,t20,t30	*/		/*...Block 5: t08,t18,t28,t38*/\
		"sub	x2,x2,#0x300			\n\t"/* r00 */\
		"ldp	q0,q1,[x2      ]		\n\t	ldp	q8 ,q9 ,[x2,#0x80]	\n\t"\
		"ldp	q2,q3,[x2,#0x40]		\n\t	ldp	q10,q11,[x2,#0xc0]	\n\t"\
		"ldp	q4,q5,[x2,#0x20]		\n\t	ldp	q12,q20,[x2,#0xa0]	\n\t"\
		"ldp	q6,q7,[x2,#0x60]		\n\t	ldp	q21,q15,[x2,#0xe0]	\n\t"\
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
/* c08 */"ldp	q18,q19,[x10,#0xa0]		\n\t	fadd	v15.2d,v13.2d,v15.2d	\n\t"/* f=d+f */\
/*0=0-5*/"fsub	v16.2d,v0.2d,v5.2d		\n\t	fsub	v8.2d ,v20.2d,v15.2d	\n\t"/* 8=8-f */\
/*1=1-4*/"fsub	v17.2d,v1.2d,v4.2d		\n\t	fsub	v9.2d ,v21.2d,v23.2d	\n\t"/* 9=9-d */\
/*5=0+5*/"fadd	v5.2d ,v0.2d,v5.2d		\n\t	fsub	v12.2d,v10.2d,v14.2d	\n\t"/* a=a-e */\
/*4=1+4*/"fadd	v4.2d ,v1.2d,v4.2d		\n\t	fsub	v13.2d,v11.2d,v22.2d	\n\t"/* b=b-c */\
		"										ldp	q24,q25,[x10,#0x0e0]	\n\t"/* c04 = [a,b] */\
		"										ldp	q26,q27,[x10,#0x120]	\n\t"/* c0C = [c,s] */\
		 "fmul	v0.2d,v17.2d,v18.2d		\n\t	fadd	v15.2d,v20.2d,v15.2d	\n\t"/* f=8+f */\
		 "fmul	v1.2d,v5.2d ,v18.2d		\n\t	fadd	v23.2d,v21.2d,v23.2d	\n\t"/* d=9+d */\
/*1c-5s*/"fmls	v0.2d,v5.2d ,v19.2d		\n\t	fadd	v14.2d,v10.2d,v14.2d	\n\t"/* e=a+e */\
/*1s+5c*/"fmla	v1.2d,v17.2d,v19.2d		\n\t	fadd	v22.2d,v11.2d,v22.2d	\n\t"/* c=b+c ... switch a-f to decimal 10-15 in ensuing twiddle-muls to avoid confusion with sincos data */\
		"ldp	q28,q29,[x1       ]		\n\t	fmul	v10.2d,v23.2d,v24.2d	\n\t"\
		"trn1	v30.2d,v6.2d,v28.2d		\n\t	fmul	v11.2d,v22.2d,v24.2d	\n\t"\
		"trn1	v31.2d,v7.2d,v29.2d		\n\t	fmls	v10.2d,v22.2d,v25.2d	\n\t"/* 13a-12b */\
		"trn2	v28.2d,v6.2d,v28.2d		\n\t	fmla	v11.2d,v23.2d,v25.2d	\n\t"/* 13b+12a */\
		"trn2	v29.2d,v7.2d,v29.2d		\n\t	fmul	v20.2d,v12.2d,v26.2d	\n\t"\
		"stp	q30,q31,[x0       ]		\n\t	fmul	v21.2d,v15.2d,v26.2d	\n\t"\
		"stp	q28,q29,[x1       ]		\n\t	fmls	v20.2d,v15.2d,v27.2d	\n\t"/* 10c-15s */\
	/* v6,7,18,19 free - use those as temps for next lcol TRN-block to free 28-31 up for rcol-TRNs: */\
		"ldp	q18,q19,[x1,#0x080]		\n\t	fmla	v21.2d,v12.2d,v27.2d	\n\t"/* 10s+15c */\
		"trn1	v6.2d ,v1.2d,v18.2d		\n\t	ldp	q24,q25,[x10,#0x100]	\n\t"/* c14 = [a,b] */\
		"trn1	v7.2d ,v0.2d,v19.2d		\n\t	ldp	q26,q27,[x10,#0x140]	\n\t"/* c1C = [c,s] */\
		"trn2	v18.2d,v1.2d,v18.2d		\n\t	ldp	q28,q29,[x1,#0x040]		\n\t"\
		"trn2	v19.2d,v0.2d,v19.2d		\n\t	trn1	v30.2d,v11.2d,v28.2d\n\t"\
		"stp	q6 ,q7 ,[x0,#0x080]		\n\t	trn1	v31.2d,v10.2d,v29.2d\n\t"\
		"stp	q18,q19,[x1,#0x080]		\n\t	trn2	v28.2d,v11.2d,v28.2d\n\t"\
/* c10=[a,b] */"ldp q5 ,q17,[x10,#0x80]	\n\t	trn2	v29.2d,v10.2d,v29.2d\n\t"\
/* c18=[c,s] */"ldp q18,q19,[x10,#0xc0]	\n\t	stp	q30,q31,[x0,#0x040]		\n\t"\
		"fmul	v0.2d,v3.2d ,v5.2d		\n\t	stp	q28,q29,[x1,#0x040]		\n\t"\
		"fmul	v1.2d,v2.2d ,v5.2d		\n\t	ldp	q28,q29,[x1,#0x0c0]		\n\t"\
/*3a-2b*/"fmls	v0.2d,v2.2d ,v17.2d		\n\t	trn1	v30.2d,v21.2d,v28.2d\n\t"\
/*3b+2a*/"fmla	v1.2d,v3.2d ,v17.2d		\n\t	trn1	v31.2d,v20.2d,v29.2d\n\t"\
		"fmul	v6.2d,v4.2d ,v18.2d		\n\t	trn2	v28.2d,v21.2d,v28.2d\n\t"\
		"fmul	v7.2d,v16.2d,v18.2d		\n\t	trn2	v29.2d,v20.2d,v29.2d\n\t"\
/*4c-0s*/"fmls	v6.2d,v16.2d,v19.2d		\n\t	stp	q30,q31,[x0,#0x0c0]		\n\t"\
/*4s+0c*/"fmla	v7.2d,v4.2d ,v19.2d		\n\t	stp	q28,q29,[x1,#0x0c0]		\n\t"\
		"										fmul	v10.2d,v9.2d ,v24.2d	\n\t"\
		"										fmul	v11.2d,v13.2d,v24.2d	\n\t"\
		"ldp	q28,q29,[x1,#0x100]		\n\t	fmls	v10.2d,v13.2d,v25.2d	\n\t"/*  9a-11b */\
		"trn1	v30.2d,v1.2d,v28.2d		\n\t	fmla	v11.2d,v9.2d ,v25.2d	\n\t"/*  9b+11a */\
		"trn1	v31.2d,v0.2d,v29.2d		\n\t	fmul	v20.2d,v14.2d,v26.2d	\n\t"/* lcol: add0 + p1  */\
		"trn2	v28.2d,v1.2d,v28.2d		\n\t	fmul	v21.2d,v8.2d ,v26.2d	\n\t"/* lcol: add0 + p9  */\
		"trn2	v29.2d,v0.2d,v29.2d		\n\t	fmls	v20.2d,v8.2d ,v27.2d	\n\t"/* 14c- 8s */\
		"stp	q30,q31,[x0,#0x100]		\n\t	fmla	v21.2d,v14.2d,v27.2d	\n\t"/* 14s+ 8c */\
		"stp	q28,q29,[x1,#0x100]		\n\t	ldp	q28,q29,[x1,#0x140]		\n\t"\
	/* v0,1,18,19 free - use those as temps for next lcol TRN-block to free 28-31 up for rcol-TRNs: */\
		"ldp	q18,q19,[x1,#0x180]		\n\t	trn1	v30.2d,v11.2d,v28.2d\n\t"\
		"trn1	v0.2d ,v7.2d,v18.2d		\n\t	trn1	v31.2d,v10.2d,v29.2d\n\t"\
		"trn1	v1.2d ,v6.2d,v19.2d		\n\t	trn2	v28.2d,v11.2d,v28.2d\n\t"\
		"trn2	v18.2d,v7.2d,v18.2d		\n\t	trn2	v29.2d,v10.2d,v29.2d\n\t"\
		"trn2	v19.2d,v6.2d,v19.2d		\n\t	stp	q30,q31,[x0,#0x140]		\n\t"\
		"stp	q0 ,q1 ,[x0,#0x180]		\n\t	stp	q28,q29,[x1,#0x140]		\n\t"\
		"stp	q18,q19,[x1,#0x180]		\n\t	ldp	q28,q29,[x1,#0x1c0]		\n\t"\
											"	trn1	v30.2d,v21.2d,v28.2d\n\t"\
											"	trn1	v31.2d,v20.2d,v29.2d\n\t"\
											"	trn2	v28.2d,v21.2d,v28.2d\n\t"\
											"	trn2	v29.2d,v20.2d,v29.2d\n\t"\
											"	stp	q30,q31,[x0,#0x1c0]		\n\t"\
											"	stp	q28,q29,[x1,#0x1c0]		\n\t"\
		"ldp q29,q30,[x10,#-0x10]		\n\t"/* Reload isrt2 into v29, will discard upper-register load data */\
		/*...Block 3: t04,t14,t24,t34*/			/*...Block 7: t0C,t1C,t2C,t3C*/\
		"add	x2,x2,#0x200			\n\t"/* r20 */\
		"ldp	q0,q1,[x2,#0x20]		\n\t"\
		"ldp	q2,q3,[x2,#0x60]		\n\t"\
/*cc0:*/"ldp	q18,q19,[x10]			\n\t"\
		"fmul	v4.2d,v0.2d,v18.2d		\n\t	ldp	q8 ,q9 ,[x2,#0xa0]		\n\t"\
		"fmul	v5.2d,v1.2d,v18.2d		\n\t	ldp	q10,q11,[x2,#0xe0]		\n\t"\
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
		"ldp	q0,q1,[x2      ]		\n\t	fsub	v12.2d,v12.2d,v24.2d	\n\t"\
		"ldp	q2,q3,[x2,#0x40]		\n\t	fsub	v13.2d,v13.2d,v25.2d	\n\t"\
		"fadd	v16.2d,v3.2d,v2.2d		\n\t	ldp	q8 ,q9 ,[x2,#0x80]		\n\t"\
		"fsub	v17.2d,v3.2d,v2.2d		\n\t	ldp	q10,q11,[x2,#0xc0]		\n\t"\
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
		"ldp	q2,q3,[x10,#0x1a0]		\n\t"	/* Let c06 = [c,s], c0E = [a,b]: */\
		"fmul	v24.2d,v4.2d,v0.2d		\n\t	ldp	q8 ,q9 ,[x10,#0x1e0]		\n\t"\
		"fmul	v25.2d,v5.2d,v0.2d		\n\t	ldp	q10,q11,[x10,#0x220]		\n\t"\
/*4c+5s*/"fmla	v24.2d,v5.2d,v1.2d		\n\t	fmul	v26.2d,v12.2d,v8.2d		\n\t"\
/*5c-4s*/"fmls	v25.2d,v4.2d,v1.2d		\n\t	fmul	v27.2d,v13.2d,v8.2d		\n\t"\
		"fmul	v4.2d,v17.2d,v2.2d		\n\t	fmla	v26.2d,v13.2d,v9.2d		\n\t"/*12c+13s*/\
		"fmul	v5.2d,v7.2d ,v2.2d		\n\t	fmls	v27.2d,v12.2d,v9.2d		\n\t"/*13c-12s*/\
/*1a-7b*/"fmls	v4.2d,v7.2d ,v3.2d		\n\t	fmul	v12.2d,v21.2d,v10.2d	\n\t"\
/*7a+1b*/"fmla	v5.2d,v17.2d,v3.2d		\n\t	fmul	v13.2d,v15.2d,v10.2d	\n\t"\
		"ldp	q28,q29,[x1,#0x020]		\n\t	fmls	v12.2d,v15.2d,v11.2d	\n\t"/* 9a-15b*/\
		"trn1	v30.2d,v24.2d,v28.2d	\n\t	fmla	v13.2d,v21.2d,v11.2d	\n\t"/*15a+ 9b*/\
		"trn1	v31.2d,v25.2d,v29.2d	\n\t"/* v8,9,10,11 free - use those as temps for TRN-block: */\
		"trn2	v28.2d,v24.2d,v28.2d	\n\t	ldp	q8 ,q9 ,[x1,#0x060]		\n\t"\
		"trn2	v29.2d,v25.2d,v29.2d	\n\t	trn1	v10.2d,v26.2d,v8.2d	\n\t"\
		"stp	q30,q31,[x0,#0x020]		\n\t	trn1	v11.2d,v27.2d,v9.2d	\n\t"\
		"stp	q28,q29,[x1,#0x020]		\n\t	trn2	v26.2d,v26.2d,v8.2d	\n\t"\
		"ldp	q28,q29,[x1,#0x0a0]		\n\t	trn2	v27.2d,v27.2d,v9.2d	\n\t"\
		"trn1	v30.2d,v5.2d ,v28.2d	\n\t	stp	q10,q11,[x0,#0x060]		\n\t"\
		"trn1	v31.2d,v4.2d ,v29.2d	\n\t	stp	q26,q27,[x1,#0x060]		\n\t"\
		"trn2	v28.2d,v5.2d ,v28.2d	\n\t"/* v8,9,10,11 free - use those as temps for TRN-block: */\
		"trn2	v29.2d,v4.2d ,v29.2d	\n\t	ldp	q8 ,q9 ,[x1,#0x0e0]		\n\t"\
		"stp	q30,q31,[x0,#0x0a0]		\n\t	trn1	v10.2d,v13.2d,v8.2d	\n\t"\
		"stp	q28,q29,[x1,#0x0a0]		\n\t	trn1	v11.2d,v12.2d,v9.2d	\n\t"\
		/* Let c12 = [c,s], c1A = [a,b]: */"	trn2	v13.2d,v13.2d,v8.2d	\n\t"\
		"ldp	q0,q1,[x10,#0x180]		\n\t	trn2	v12.2d,v12.2d,v9.2d	\n\t"\
		"ldp	q2,q3,[x10,#0x1c0]		\n\t	stp	q10,q11,[x0,#0x0e0]		\n\t"\
		"fmul	v24.2d,v18.2d,v0.2d		\n\t	stp	q13,q12,[x1,#0x0e0]		\n\t"\
		"fmul	v25.2d,v19.2d,v0.2d		\n\t"/* Let c16 = [c,s], c1E = [a,b]: */\
/*2c+3s*/"fmla	v24.2d,v19.2d,v1.2d		\n\t	ldp	q8 ,q9 ,[x10,#0x200]		\n\t"\
/*3c-2s*/"fmls	v25.2d,v18.2d,v1.2d		\n\t	ldp	q10,q11,[x10,#0x240]		\n\t"\
		"fmul	v4.2d,v6.2d ,v2.2d		\n\t	fmul	v26.2d,v22.2d,v8.2d		\n\t"\
		"fmul	v5.2d,v16.2d,v2.2d		\n\t	fmul	v27.2d,v23.2d,v8.2d		\n\t"\
/*6a-0b*/"fmls	v4.2d,v16.2d,v3.2d		\n\t	fmla	v26.2d,v23.2d,v9.2d		\n\t"/*10c+11s*/\
/*0a+6b*/"fmla	v5.2d,v6.2d ,v3.2d		\n\t	fmls	v27.2d,v22.2d,v9.2d		\n\t"/*11c-10s*/\
		"ldp	q28,q29,[x1,#0x120]		\n\t	fmul	v12.2d,v14.2d,v10.2d	\n\t"\
		"trn1	v30.2d,v24.2d,v28.2d	\n\t	fmul	v13.2d,v20.2d,v10.2d	\n\t"\
		"trn1	v31.2d,v25.2d,v29.2d	\n\t	fmls	v12.2d,v20.2d,v11.2d	\n\t"/*14a-8b*/\
		"trn2	v28.2d,v24.2d,v28.2d	\n\t	fmla	v13.2d,v14.2d,v11.2d	\n\t"/*8a+14b*/\
		"trn2	v29.2d,v25.2d,v29.2d	\n\t"/* v8,9,10,11 free - use those as temps for TRN-block: */\
		"stp	q30,q31,[x0,#0x120]		\n\t	ldp	q8 ,q9 ,[x1,#0x160]		\n\t"\
		"stp	q28,q29,[x1,#0x120]		\n\t	trn1	v10.2d,v26.2d,v8.2d	\n\t"\
		"ldp	q28,q29,[x1,#0x1a0]		\n\t	trn1	v11.2d,v27.2d,v9.2d	\n\t"\
		"trn1	v30.2d,v5.2d ,v28.2d	\n\t	trn2	v26.2d,v26.2d,v8.2d	\n\t"\
		"trn1	v31.2d,v4.2d ,v29.2d	\n\t	trn2	v27.2d,v27.2d,v9.2d	\n\t"\
		"trn2	v28.2d,v5.2d ,v28.2d	\n\t	stp	q10,q11,[x0,#0x160]		\n\t"\
		"trn2	v29.2d,v4.2d ,v29.2d	\n\t	stp	q26,q27,[x1,#0x160]		\n\t"\
		"stp	q30,q31,[x0,#0x1a0]		\n\t	ldp	q8 ,q9 ,[x1,#0x1e0]		\n\t"\
		"stp	q28,q29,[x1,#0x1a0]		\n\t	trn1	v10.2d,v13.2d,v8.2d	\n\t"\
		"										trn1	v11.2d,v12.2d,v9.2d	\n\t"\
		"										trn2	v13.2d,v13.2d,v8.2d	\n\t"\
		"										trn2	v12.2d,v12.2d,v9.2d	\n\t"\
		"										stp	q10,q11,[x0,#0x1e0]		\n\t"\
		"										stp	q13,q12,[x1,#0x1e0]		\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__add1] "m" (Xadd1)\
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__r00] "m" (Xr00)\
		 ,[__r08] "m" (Xr08)\
		 ,[__r10] "m" (Xr10)\
		 ,[__c01] "m" (Xc01)\
		 ,[__c02] "m" (Xc02)\
		 ,[__c04] "m" (Xc04)\
		 ,[__c06] "m" (Xc06)\
		 ,[__c08] "m" (Xc08)\
		 ,[__c0A] "m" (Xc0A)\
		 ,[__c0C] "m" (Xc0C)\
		 ,[__c0E] "m" (Xc0E)\
		 ,[__c10] "m" (Xc10)\
		 ,[__c12] "m" (Xc12)\
		 ,[__c14] "m" (Xc14)\
		 ,[__c16] "m" (Xc16)\
		 ,[__c18] "m" (Xc18)\
		 ,[__c1A] "m" (Xc1A)\
		 ,[__c1C] "m" (Xc1C)\
		 ,[__c1E] "m" (Xc1E)\
		: "cc","memory","x0","x1","x2","x10","x11",\
		"v0","v1","v2","v3","v4","v5","v6","v7","v8","v9","v10","v11","v12","v13","v14","v15",\
		"v16","v17","v18","v19","v20","v21","v22","v23","v24","v25","v26","v27","v28","v29","v30","v31"	/* Clobbered registers */\
	);\
	}

#elif defined(USE_AVX512)	// AVX512 implements a 512-bit-register version of the the AVX2/FMA3-macro (non-ALL_FMA version)

  #ifdef USE_16_REG 	// Default is 32-SIMD-register version further down; in my Jan 2017 tests on KNL found 1-2%
  						// overall gain [i.e. 8-10x that for this macro since that was only one that got toggled and
  						// accounts for ~10% of the overall runtime] from using all 32
	/***
		Dec 2018: NOTE that we keep the experimental-but-slow gather-based transposed in this version
		for reference purposes; alos, this version appears broken, get fatal ROE on iter 1 when use it.
	***/
	// Cost [vector-ops only]: 740 MEM (128 in transpose section, 610 in DFT proper), 540 ARITHMETIC
	#define SSE2_RADIX32_WRAPPER_DIF(Xadd0,Xadd1,Xadd2,Xadd3,Xadd4,Xadd5,Xadd6,Xadd7,Xr00,Xr10,Xr20,Xr30,Xisrt2,Xcc0,Xc00,Xc01,Xc02,Xc03,Xc05,Xc07)\
	{\
	__asm__ volatile (\
		"movq	%[__r00] ,%%rbx	\n\t	leaq	0x40(%%rbx),%%rcx	\n\t"\
		"movq	%[__add0],%%rax	\n\t"\
	/**** Start with 8-way interleaving: ****/\
		/* Auxiliary register data needed for columnwise loads: */\
		"movq	$0x0706050403020100,%%rsi	\n\t"/* 64-bit register w/byte offsets 0-7, bytes ordered left-to-right in decreasing significance */\
			"vmovq		%%rsi,%%xmm8 		\n\t"/* Copy byte pattern to low qword (64 bits) of ymm8 [NB: avx-512 only supports MOVQ to/from 128-bit vector regs] */\
			"vpmovzxbd	%%xmm8,%%ymm8		\n\t"/* vector-index offsets: ymm8 = [0,1,2,3,4,5,6,7] in 32-bit form in low 8 dwords */\
			"vpslld	$9,%%ymm8,%%ymm8		\n\t"/* The above bytewise offsets need scale *512 to get the needed ones - would include but
											e.g. 1<<9 overflows 1 byte - but x86 ISA only permits scale factors 1,2,4,8, so <<= 8 here. */\
		/* Mask-reg zmm9 = 11...11 - this is stupidly zeroed each time we do gather-load, so need to reinit: */\
		"movl	$-1,%%esi	\n\t"/* Init opmask k1 (Only need the low byte) */\
		/* Gather instruction sets mask-reg = 0, so must re-init opmask prior to each invocation */\
	/* a[j+p0]: Inputs from add0+[0,0x200,0x400,0x600,0x800,0xa00,0xc00,0xe00]. Outputs into local store at r0+[same byte offsets]: */\
		/* Real parts: */\
		"kmovw	%%esi,%%k1	\n\t	vgatherdpd 0x00(%%rax,%%ymm8),%%zmm0%{%%k1%}	\n\t"/* Col 0.re */\
		"kmovw	%%esi,%%k1	\n\t	vgatherdpd 0x10(%%rax,%%ymm8),%%zmm1%{%%k1%}	\n\t"/* Col 2.re */\
		"kmovw	%%esi,%%k1	\n\t	vgatherdpd 0x08(%%rax,%%ymm8),%%zmm2%{%%k1%}	\n\t"/* Col 1.re */\
		"kmovw	%%esi,%%k1	\n\t	vgatherdpd 0x18(%%rax,%%ymm8),%%zmm3%{%%k1%}	\n\t"/* Col 3.re */\
		"kmovw	%%esi,%%k1	\n\t	vgatherdpd 0x20(%%rax,%%ymm8),%%zmm4%{%%k1%}	\n\t"/* Col 4.re */\
		"kmovw	%%esi,%%k1	\n\t	vgatherdpd 0x30(%%rax,%%ymm8),%%zmm5%{%%k1%}	\n\t"/* Col 6.re */\
		"kmovw	%%esi,%%k1	\n\t	vgatherdpd 0x28(%%rax,%%ymm8),%%zmm6%{%%k1%}	\n\t"/* Col 5.re */\
		"kmovw	%%esi,%%k1	\n\t	vgatherdpd 0x38(%%rax,%%ymm8),%%zmm7%{%%k1%}	\n\t"/* Col 7.re */\
												/* Imag parts: */\
												"kmovw	%%esi,%%k1	\n\t	vgatherdpd 0x40(%%rax,%%ymm8),%%zmm10%{%%k1%}	\n\t"/* Col 0.im */\
												"kmovw	%%esi,%%k1	\n\t	vgatherdpd 0x50(%%rax,%%ymm8),%%zmm11%{%%k1%}	\n\t"/* Col 1.im */\
												"kmovw	%%esi,%%k1	\n\t	vgatherdpd 0x48(%%rax,%%ymm8),%%zmm12%{%%k1%}	\n\t"/* Col 2.im */\
												"kmovw	%%esi,%%k1	\n\t	vgatherdpd 0x58(%%rax,%%ymm8),%%zmm13%{%%k1%}	\n\t"/* Col 3.im */\
												"kmovw	%%esi,%%k1	\n\t	vgatherdpd 0x60(%%rax,%%ymm8),%%zmm14%{%%k1%}	\n\t"/* Col 4.im */\
												"kmovw	%%esi,%%k1	\n\t	vgatherdpd 0x70(%%rax,%%ymm8),%%zmm15%{%%k1%}	\n\t"/* Col 5.im */\
												"kmovw	%%esi,%%k1	\n\t	vgatherdpd 0x68(%%rax,%%ymm8),%%zmm16%{%%k1%}	\n\t"/* Col 6.im */\
												"kmovw	%%esi,%%k1	\n\t	vgatherdpd 0x78(%%rax,%%ymm8),%%zmm17%{%%k1%}	\n\t"/* Col 7.im */\
		/* Write original columns back as rows to the local store - sort code lines here by mem-offset: */\
	/**** Outs 0-7 [in terms of the rowwise ordering here, not the reg-number] go into slots 0,4,8,12,16,20,24,28: ****/\
		"vmovaps 	%%zmm0,     (%%rbx)			\n\t		vmovaps %%zmm10,     (%%rcx)	\n\t"\
		"vmovaps 	%%zmm4,0x200(%%rbx)			\n\t		vmovaps %%zmm14,0x200(%%rcx)	\n\t"\
		"vmovaps 	%%zmm1,0x400(%%rbx)			\n\t		vmovaps %%zmm11,0x400(%%rcx)	\n\t"\
		"vmovaps 	%%zmm5,0x600(%%rbx)			\n\t		vmovaps %%zmm15,0x600(%%rcx)	\n\t"\
		"vmovaps 	%%zmm2,0x800(%%rbx)			\n\t		vmovaps %%zmm12,0x800(%%rcx)	\n\t"\
		"vmovaps 	%%zmm6,0xa00(%%rbx)			\n\t		vmovaps %%zmm16,0xa00(%%rcx)	\n\t"\
		"vmovaps 	%%zmm3,0xc00(%%rbx)			\n\t		vmovaps %%zmm13,0xc00(%%rcx)	\n\t"\
		"vmovaps 	%%zmm7,0xe00(%%rbx)			\n\t		vmovaps %%zmm17,0xe00(%%rcx)	\n\t"\
	/* a[j+p2]: Inputs from add0+[0,0x100,0x200,0x300,0x400,0x500,0x600,0x700]+0x80. Outputs into r1+[same byte offsets]: */\
	"addq	$0x80,%%rax	\n\t"\
		/* Real parts: */\
		"kmovw	%%esi,%%k1	\n\t	vgatherdpd 0x00(%%rax,%%ymm8),%%zmm0%{%%k1%}	\n\t"/* Col 8.re */\
		"kmovw	%%esi,%%k1	\n\t	vgatherdpd 0x10(%%rax,%%ymm8),%%zmm1%{%%k1%}	\n\t"/* Col a.re */\
		"kmovw	%%esi,%%k1	\n\t	vgatherdpd 0x08(%%rax,%%ymm8),%%zmm2%{%%k1%}	\n\t"/* Col 9.re */\
		"kmovw	%%esi,%%k1	\n\t	vgatherdpd 0x18(%%rax,%%ymm8),%%zmm3%{%%k1%}	\n\t"/* Col b.re */\
		"kmovw	%%esi,%%k1	\n\t	vgatherdpd 0x20(%%rax,%%ymm8),%%zmm4%{%%k1%}	\n\t"/* Col c.re */\
		"kmovw	%%esi,%%k1	\n\t	vgatherdpd 0x30(%%rax,%%ymm8),%%zmm5%{%%k1%}	\n\t"/* Col e.re */\
		"kmovw	%%esi,%%k1	\n\t	vgatherdpd 0x28(%%rax,%%ymm8),%%zmm6%{%%k1%}	\n\t"/* Col d.re */\
		"kmovw	%%esi,%%k1	\n\t	vgatherdpd 0x38(%%rax,%%ymm8),%%zmm7%{%%k1%}	\n\t"/* Col f.re */\
												/* Imag parts: */\
												"kmovw	%%esi,%%k1	\n\t	vgatherdpd 0x40(%%rax,%%ymm8),%%zmm10%{%%k1%}	\n\t"/* Col 0.im */\
												"kmovw	%%esi,%%k1	\n\t	vgatherdpd 0x50(%%rax,%%ymm8),%%zmm11%{%%k1%}	\n\t"/* Col 1.im */\
												"kmovw	%%esi,%%k1	\n\t	vgatherdpd 0x48(%%rax,%%ymm8),%%zmm12%{%%k1%}	\n\t"/* Col 2.im */\
												"kmovw	%%esi,%%k1	\n\t	vgatherdpd 0x58(%%rax,%%ymm8),%%zmm13%{%%k1%}	\n\t"/* Col 3.im */\
												"kmovw	%%esi,%%k1	\n\t	vgatherdpd 0x60(%%rax,%%ymm8),%%zmm14%{%%k1%}	\n\t"/* Col 4.im */\
												"kmovw	%%esi,%%k1	\n\t	vgatherdpd 0x70(%%rax,%%ymm8),%%zmm15%{%%k1%}	\n\t"/* Col 5.im */\
												"kmovw	%%esi,%%k1	\n\t	vgatherdpd 0x68(%%rax,%%ymm8),%%zmm16%{%%k1%}	\n\t"/* Col 6.im */\
												"kmovw	%%esi,%%k1	\n\t	vgatherdpd 0x78(%%rax,%%ymm8),%%zmm17%{%%k1%}	\n\t"/* Col 7.im */\
		/* Write original columns back as rows to the local store: */\
	/**** Outs 0-7 [in terms of the rowwise ordering here, not the reg-number] go into slots 2,6,10,14,18,22,26,30: ****/\
	"addq	$0x100,%%rbx	\n\t	addq	$0x100,%%rcx	\n\t"\
		"vmovaps 	%%zmm0,     (%%rbx)			\n\t		vmovaps %%zmm10,     (%%rcx)	\n\t"\
		"vmovaps 	%%zmm4,0x200(%%rbx)			\n\t		vmovaps %%zmm14,0x200(%%rcx)	\n\t"\
		"vmovaps 	%%zmm1,0x400(%%rbx)			\n\t		vmovaps %%zmm11,0x400(%%rcx)	\n\t"\
		"vmovaps 	%%zmm5,0x600(%%rbx)			\n\t		vmovaps %%zmm15,0x600(%%rcx)	\n\t"\
		"vmovaps 	%%zmm2,0x800(%%rbx)			\n\t		vmovaps %%zmm12,0x800(%%rcx)	\n\t"\
		"vmovaps 	%%zmm6,0xa00(%%rbx)			\n\t		vmovaps %%zmm16,0xa00(%%rcx)	\n\t"\
		"vmovaps 	%%zmm3,0xc00(%%rbx)			\n\t		vmovaps %%zmm13,0xc00(%%rcx)	\n\t"\
		"vmovaps 	%%zmm7,0xe00(%%rbx)			\n\t		vmovaps %%zmm17,0xe00(%%rcx)	\n\t"\
	/* a[j+p4]: Inputs from add0+[0,0x100,0x200,0x300,0x400,0x500,0x600,0x700]+0x80. Outputs into r2+[same byte offsets]: */\
	"addq	$0x80,%%rax		\n\t"\
		/* Real parts: */\
		"kmovw	%%esi,%%k1	\n\t	vgatherdpd 0x00(%%rax,%%ymm8),%%zmm0%{%%k1%}	\n\t"/* Col 8.re */\
		"kmovw	%%esi,%%k1	\n\t	vgatherdpd 0x10(%%rax,%%ymm8),%%zmm1%{%%k1%}	\n\t"/* Col a.re */\
		"kmovw	%%esi,%%k1	\n\t	vgatherdpd 0x08(%%rax,%%ymm8),%%zmm2%{%%k1%}	\n\t"/* Col 9.re */\
		"kmovw	%%esi,%%k1	\n\t	vgatherdpd 0x18(%%rax,%%ymm8),%%zmm3%{%%k1%}	\n\t"/* Col b.re */\
		"kmovw	%%esi,%%k1	\n\t	vgatherdpd 0x20(%%rax,%%ymm8),%%zmm4%{%%k1%}	\n\t"/* Col c.re */\
		"kmovw	%%esi,%%k1	\n\t	vgatherdpd 0x30(%%rax,%%ymm8),%%zmm5%{%%k1%}	\n\t"/* Col e.re */\
		"kmovw	%%esi,%%k1	\n\t	vgatherdpd 0x28(%%rax,%%ymm8),%%zmm6%{%%k1%}	\n\t"/* Col d.re */\
		"kmovw	%%esi,%%k1	\n\t	vgatherdpd 0x38(%%rax,%%ymm8),%%zmm7%{%%k1%}	\n\t"/* Col f.re */\
												/* Imag parts: */\
												"kmovw	%%esi,%%k1	\n\t	vgatherdpd 0x40(%%rax,%%ymm8),%%zmm10%{%%k1%}	\n\t"/* Col 0.im */\
												"kmovw	%%esi,%%k1	\n\t	vgatherdpd 0x50(%%rax,%%ymm8),%%zmm11%{%%k1%}	\n\t"/* Col 1.im */\
												"kmovw	%%esi,%%k1	\n\t	vgatherdpd 0x48(%%rax,%%ymm8),%%zmm12%{%%k1%}	\n\t"/* Col 2.im */\
												"kmovw	%%esi,%%k1	\n\t	vgatherdpd 0x58(%%rax,%%ymm8),%%zmm13%{%%k1%}	\n\t"/* Col 3.im */\
												"kmovw	%%esi,%%k1	\n\t	vgatherdpd 0x60(%%rax,%%ymm8),%%zmm14%{%%k1%}	\n\t"/* Col 4.im */\
												"kmovw	%%esi,%%k1	\n\t	vgatherdpd 0x70(%%rax,%%ymm8),%%zmm15%{%%k1%}	\n\t"/* Col 5.im */\
												"kmovw	%%esi,%%k1	\n\t	vgatherdpd 0x68(%%rax,%%ymm8),%%zmm16%{%%k1%}	\n\t"/* Col 6.im */\
												"kmovw	%%esi,%%k1	\n\t	vgatherdpd 0x78(%%rax,%%ymm8),%%zmm17%{%%k1%}	\n\t"/* Col 7.im */\
		/* Write original columns back as rows to the local store: */\
	/**** Outs 0-7 [in terms of the rowwise ordering here, not the reg-number] go into slots 1,5,9,13,17,21,25,29: ****/\
	"subq	$0x80,%%rbx	\n\t	subq	$0x80,%%rcx	\n\t"\
		"vmovaps 	%%zmm0,     (%%rbx)			\n\t		vmovaps %%zmm10,     (%%rcx)	\n\t"\
		"vmovaps 	%%zmm4,0x200(%%rbx)			\n\t		vmovaps %%zmm14,0x200(%%rcx)	\n\t"\
		"vmovaps 	%%zmm1,0x400(%%rbx)			\n\t		vmovaps %%zmm11,0x400(%%rcx)	\n\t"\
		"vmovaps 	%%zmm5,0x600(%%rbx)			\n\t		vmovaps %%zmm15,0x600(%%rcx)	\n\t"\
		"vmovaps 	%%zmm2,0x800(%%rbx)			\n\t		vmovaps %%zmm12,0x800(%%rcx)	\n\t"\
		"vmovaps 	%%zmm6,0xa00(%%rbx)			\n\t		vmovaps %%zmm16,0xa00(%%rcx)	\n\t"\
		"vmovaps 	%%zmm3,0xc00(%%rbx)			\n\t		vmovaps %%zmm13,0xc00(%%rcx)	\n\t"\
		"vmovaps 	%%zmm7,0xe00(%%rbx)			\n\t		vmovaps %%zmm17,0xe00(%%rcx)	\n\t"\
	/* a[j+p6]: Inputs from add0+[0,0x100,0x200,0x300,0x400,0x500,0x600,0x700]+0x180. Outputs into r3+[same byte offsets]: */\
	"addq	$0x80,%%rax		\n\t"\
		/* Real parts: */\
		"kmovw	%%esi,%%k1	\n\t	vgatherdpd 0x00(%%rax,%%ymm8),%%zmm0%{%%k1%}	\n\t"/* Col 8.re */\
		"kmovw	%%esi,%%k1	\n\t	vgatherdpd 0x10(%%rax,%%ymm8),%%zmm1%{%%k1%}	\n\t"/* Col a.re */\
		"kmovw	%%esi,%%k1	\n\t	vgatherdpd 0x08(%%rax,%%ymm8),%%zmm2%{%%k1%}	\n\t"/* Col 9.re */\
		"kmovw	%%esi,%%k1	\n\t	vgatherdpd 0x18(%%rax,%%ymm8),%%zmm3%{%%k1%}	\n\t"/* Col b.re */\
		"kmovw	%%esi,%%k1	\n\t	vgatherdpd 0x20(%%rax,%%ymm8),%%zmm4%{%%k1%}	\n\t"/* Col c.re */\
		"kmovw	%%esi,%%k1	\n\t	vgatherdpd 0x30(%%rax,%%ymm8),%%zmm5%{%%k1%}	\n\t"/* Col e.re */\
		"kmovw	%%esi,%%k1	\n\t	vgatherdpd 0x28(%%rax,%%ymm8),%%zmm6%{%%k1%}	\n\t"/* Col d.re */\
		"kmovw	%%esi,%%k1	\n\t	vgatherdpd 0x38(%%rax,%%ymm8),%%zmm7%{%%k1%}	\n\t"/* Col f.re */\
												/* Imag parts: */\
												"kmovw	%%esi,%%k1	\n\t	vgatherdpd 0x40(%%rax,%%ymm8),%%zmm10%{%%k1%}	\n\t"/* Col 0.im */\
												"kmovw	%%esi,%%k1	\n\t	vgatherdpd 0x50(%%rax,%%ymm8),%%zmm11%{%%k1%}	\n\t"/* Col 1.im */\
												"kmovw	%%esi,%%k1	\n\t	vgatherdpd 0x48(%%rax,%%ymm8),%%zmm12%{%%k1%}	\n\t"/* Col 2.im */\
												"kmovw	%%esi,%%k1	\n\t	vgatherdpd 0x58(%%rax,%%ymm8),%%zmm13%{%%k1%}	\n\t"/* Col 3.im */\
												"kmovw	%%esi,%%k1	\n\t	vgatherdpd 0x60(%%rax,%%ymm8),%%zmm14%{%%k1%}	\n\t"/* Col 4.im */\
												"kmovw	%%esi,%%k1	\n\t	vgatherdpd 0x70(%%rax,%%ymm8),%%zmm15%{%%k1%}	\n\t"/* Col 5.im */\
												"kmovw	%%esi,%%k1	\n\t	vgatherdpd 0x68(%%rax,%%ymm8),%%zmm16%{%%k1%}	\n\t"/* Col 6.im */\
												"kmovw	%%esi,%%k1	\n\t	vgatherdpd 0x78(%%rax,%%ymm8),%%zmm17%{%%k1%}	\n\t"/* Col 7.im */\
		/* Write original columns back as rows to the local store: */\
	/**** Outs 0-7 [in terms of the rowwise ordering here, not the reg-number] go into slots 3,7,11,15,19,23,27,31: ****/\
	"addq	$0x100,%%rbx	\n\t	addq	$0x100,%%rcx	\n\t"\
		"vmovaps 	%%zmm0,     (%%rbx)			\n\t		vmovaps %%zmm10,     (%%rcx)	\n\t"\
		"vmovaps 	%%zmm4,0x200(%%rbx)			\n\t		vmovaps %%zmm14,0x200(%%rcx)	\n\t"\
		"vmovaps 	%%zmm1,0x400(%%rbx)			\n\t		vmovaps %%zmm11,0x400(%%rcx)	\n\t"\
		"vmovaps 	%%zmm5,0x600(%%rbx)			\n\t		vmovaps %%zmm15,0x600(%%rcx)	\n\t"\
		"vmovaps 	%%zmm2,0x800(%%rbx)			\n\t		vmovaps %%zmm12,0x800(%%rcx)	\n\t"\
		"vmovaps 	%%zmm6,0xa00(%%rbx)			\n\t		vmovaps %%zmm16,0xa00(%%rcx)	\n\t"\
		"vmovaps 	%%zmm3,0xc00(%%rbx)			\n\t		vmovaps %%zmm13,0xc00(%%rcx)	\n\t"\
		"vmovaps 	%%zmm7,0xe00(%%rbx)			\n\t		vmovaps %%zmm17,0xe00(%%rcx)	\n\t"\
		/************************************************************************/\
		/* Forward DIF radix-32 pass on the interleaved block1 and block2 data: */\
		/************************************************************************/\
	/*...Block 0: */\
		"movq	%[__isrt2],%%rsi				\n\t		leaq	0x11c0(%%rsi),%%rdi	/* two */	\n\t"\
	/* SSE2_RADIX4_DIF_4WRAPPER(c00,c08,c10,c18,r00) SSE2_RADIX4_DIF_4WRAPPER_2NDOFTWO(c04,c0C,c14,c1C,r08) */\
		"movq		%[__r00],%%rcx			\n\t		/*addq		$0x200,%%rcx // __r08 */	\n\t"\
		"movq		%[__c00],%%rdx			\n\t		/*addq		$0x200,%%rdx // __c04 */	\n\t"\
		"vmovaps			 (%%rcx),%%zmm0		\n\t		vmovaps		0x200(%%rcx),%%zmm8			\n\t"\
		"vmovaps		0x040(%%rcx),%%zmm1		\n\t		vmovaps		0x240(%%rcx),%%zmm9			\n\t"\
		"vmovaps		%%zmm0,%%zmm2			\n\t		vmovaps		%%zmm8,%%zmm10				\n\t"\
		"vmovaps		%%zmm1,%%zmm3			\n\t		vmovaps		%%zmm9,%%zmm11				\n\t"\
		"vmulpd			 (%%rdx),%%zmm0,%%zmm0	\n\t		vmulpd		0x200(%%rdx),%%zmm8,%%zmm8	\n\t"\
		"vmulpd			 (%%rdx),%%zmm1,%%zmm1	\n\t		vmulpd		0x200(%%rdx),%%zmm9,%%zmm9	\n\t"\
	"vfnmadd231pd	0x040(%%rdx),%%zmm7,%%zmm0	\n\t	vfnmadd231pd	0x240(%%rdx),%%zmm11,%%zmm8 \n\t"\
	" vfmadd231pd	0x040(%%rdx),%%zmm6,%%zmm1	\n\t	 vfmadd231pd	0x240(%%rdx),%%zmm10,%%zmm9 \n\t"\
		"vmovaps		%%zmm0,%%zmm2			\n\t		vmovaps		%%zmm8,%%zmm10				\n\t"\
		"vmovaps		%%zmm1,%%zmm3			\n\t		vmovaps		%%zmm9,%%zmm11				\n\t"\
		"addq		$0x080,%%rdx																\n\t"\
		"vmovaps		0x080(%%rcx),%%zmm4		\n\t		vmovaps		0x280(%%rcx),%%zmm12		\n\t"\
		"vmovaps		0x0c0(%%rcx),%%zmm5		\n\t		vmovaps		0x2c0(%%rcx),%%zmm13		\n\t"\
		"vmovaps		%%zmm4,%%zmm6			\n\t		vmovaps		%%zmm12,%%zmm14				\n\t"\
		"vmovaps		%%zmm5,%%zmm7			\n\t		vmovaps		%%zmm13,%%zmm15				\n\t"\
		"vmulpd			 (%%rdx),%%zmm4,%%zmm4	\n\t		vmulpd		0x200(%%rdx),%%zmm12,%%zmm12\n\t"\
		"vmulpd			 (%%rdx),%%zmm5,%%zmm5	\n\t		vmulpd		0x200(%%rdx),%%zmm13,%%zmm13\n\t"\
	"vfnmadd231pd	0x040(%%rdx),%%zmm7,%%zmm4	\n\t	vfnmadd231pd	0x240(%%rdx),%%zmm15,%%zmm12\n\t"\
	" vfmadd231pd	0x040(%%rdx),%%zmm6,%%zmm5	\n\t	 vfmadd231pd	0x240(%%rdx),%%zmm14,%%zmm13\n\t"\
		"vaddpd		%%zmm4,%%zmm0,%%zmm0	\n\t		vaddpd		%%zmm12,%%zmm8,%%zmm8	\n\t"\
		"vaddpd		%%zmm5,%%zmm1,%%zmm1	\n\t		vaddpd		%%zmm13,%%zmm9,%%zmm9	\n\t"\
		"vsubpd		%%zmm4,%%zmm2,%%zmm2	\n\t		vsubpd		%%zmm12,%%zmm10,%%zmm10\n\t"\
		"vsubpd		%%zmm5,%%zmm3,%%zmm3	\n\t		vsubpd		%%zmm13,%%zmm11,%%zmm11\n\t"\
		"addq		$0x100,%%rdx																\n\t"\
		"vmovaps		0x180(%%rcx),%%zmm4		\n\t		vmovaps		0x380(%%rcx),%%zmm12		\n\t"\
		"vmovaps		0x1c0(%%rcx),%%zmm5		\n\t		vmovaps		0x3c0(%%rcx),%%zmm13		\n\t"\
		"vmovaps		%%zmm4,%%zmm6			\n\t		vmovaps		%%zmm12,%%zmm14				\n\t"\
		"vmovaps		%%zmm5,%%zmm7			\n\t		vmovaps		%%zmm13,%%zmm15				\n\t"\
		"vmulpd			 (%%rdx),%%zmm4,%%zmm4	\n\t		vmulpd		0x200(%%rdx),%%zmm12,%%zmm12\n\t"\
		"vmulpd			 (%%rdx),%%zmm5,%%zmm5	\n\t		vmulpd		0x200(%%rdx),%%zmm13,%%zmm13\n\t"\
	"vfnmadd231pd	0x040(%%rdx),%%zmm7,%%zmm4	\n\t	vfnmadd231pd	0x240(%%rdx),%%zmm15,%%zmm12\n\t"\
	" vfmadd231pd	0x040(%%rdx),%%zmm6,%%zmm5	\n\t	 vfmadd231pd	0x240(%%rdx),%%zmm14,%%zmm13\n\t"\
		"vmovaps		%%zmm5,0x040(%%rcx)	\n\t		vmovaps		%%zmm13,0x240(%%rcx)		\n\t"\
		"vmovaps		%%zmm4,     (%%rcx)	\n\t		vmovaps		%%zmm12,0x200(%%rcx)		\n\t"\
		"subq		$0x080,%%rdx																\n\t"\
		"vmovaps		0x100(%%rcx),%%zmm4		\n\t		vmovaps		0x300(%%rcx),%%zmm12		\n\t"\
		"vmovaps		0x140(%%rcx),%%zmm5		\n\t		vmovaps		0x340(%%rcx),%%zmm13		\n\t"\
		"vmovaps		%%zmm4,%%zmm6			\n\t		vmovaps		%%zmm12,%%zmm14				\n\t"\
		"vmovaps		%%zmm5,%%zmm7			\n\t		vmovaps		%%zmm13,%%zmm15				\n\t"\
		"vmulpd			 (%%rdx),%%zmm4,%%zmm4	\n\t		vmulpd		0x200(%%rdx),%%zmm12,%%zmm12\n\t"\
		"vmulpd			 (%%rdx),%%zmm5,%%zmm5	\n\t		vmulpd		0x200(%%rdx),%%zmm13,%%zmm13\n\t"\
	"vfnmadd231pd	0x040(%%rdx),%%zmm7,%%zmm4	\n\t	vfnmadd231pd	0x240(%%rdx),%%zmm15,%%zmm12\n\t"\
	" vfmadd231pd	0x040(%%rdx),%%zmm6,%%zmm5	\n\t	 vfmadd231pd	0x240(%%rdx),%%zmm14,%%zmm13\n\t"\
		"vmovaps		%%zmm5,%%zmm7			\n\t		vmovaps		%%zmm13,%%zmm15				\n\t"\
		"vmovaps		%%zmm4,%%zmm6			\n\t		vmovaps		%%zmm12,%%zmm14				\n\t"\
		"vsubpd			 (%%rcx),%%zmm4,%%zmm4	\n\t		vsubpd		0x200(%%rcx),%%zmm12,%%zmm12\n\t"\
		"vsubpd		0x040(%%rcx),%%zmm5,%%zmm5	\n\t		vsubpd		0x240(%%rcx),%%zmm13,%%zmm13\n\t"\
		"vaddpd			 (%%rcx),%%zmm6,%%zmm6	\n\t		vaddpd		0x200(%%rcx),%%zmm14,%%zmm14\n\t"\
		"vaddpd		0x040(%%rcx),%%zmm7,%%zmm7	\n\t		vaddpd		0x240(%%rcx),%%zmm15,%%zmm15\n\t"\
		"vsubpd		%%zmm6,%%zmm0,%%zmm0	\n\t		vsubpd		%%zmm14,%%zmm8,%%zmm8	\n\t"\
		"vsubpd		%%zmm7,%%zmm1,%%zmm1	\n\t		vsubpd		%%zmm15,%%zmm9,%%zmm9	\n\t"\
		"vsubpd		%%zmm5,%%zmm2,%%zmm2	\n\t		vsubpd		%%zmm13,%%zmm10,%%zmm10\n\t"\
		"vsubpd		%%zmm4,%%zmm3,%%zmm3	\n\t		vsubpd		%%zmm12,%%zmm11,%%zmm11\n\t"\
	"vmovaps	%%zmm12,(%%rcx) \n\t"/* spill zmm12 to make room for two */"vmovaps	(%%rdi),%%zmm12 \n\t"/* two */\
	"vfmadd132pd	%%zmm12,%%zmm0,%%zmm6		\n\t	vfmadd132pd	%%zmm12,%%zmm8 ,%%zmm14			\n\t"\
	"vfmadd132pd	%%zmm12,%%zmm1,%%zmm7		\n\t	vfmadd132pd	%%zmm12,%%zmm9 ,%%zmm15			\n\t"\
	"vfmadd132pd	%%zmm12,%%zmm2,%%zmm5		\n\t	vfmadd132pd	%%zmm12,%%zmm10,%%zmm13			\n\t"\
	"vfmadd132pd	%%zmm12,%%zmm3,%%zmm4		\n\t	vfmadd132pd	(%%rcx),%%zmm11,%%zmm12			\n\t"\
		"												vmovaps		%%zmm14,0x200(%%rcx)			\n\t"\
		"												vmovaps		%%zmm15,0x240(%%rcx)			\n\t"\
		"												vmovaps		%%zmm10,%%zmm14					\n\t"\
		"												vmovaps		%%zmm13,%%zmm15					\n\t"\
		"												vsubpd		%%zmm12,%%zmm10,%%zmm10	\n\t"\
		"												vsubpd		%%zmm11,%%zmm13,%%zmm13	\n\t"\
		"												vaddpd		%%zmm12,%%zmm14,%%zmm14	\n\t"\
		"												vaddpd		%%zmm11,%%zmm15,%%zmm15	\n\t"\
		/*****	SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r00) *****/\
		"												vmovaps		0x200(%%rcx),%%zmm11			\n\t"\
		"												vmovaps		0x240(%%rcx),%%zmm12			\n\t"\
	"vmovaps %%zmm13,(%%rcx) \n\t"/* spill zmm13 to make room for isrt2 */"vmovaps (%%rsi),%%zmm13 	\n\t"/*isrt2*/\
		"vsubpd		%%zmm11,%%zmm6,%%zmm6		\n\t	vfnmadd231pd	%%zmm13,%%zmm10,%%zmm2		\n\t"\
		"vsubpd		%%zmm9,%%zmm0,%%zmm0		\n\t	vfnmadd231pd	%%zmm13,%%zmm15,%%zmm5		\n\t"\
		"vsubpd		%%zmm12,%%zmm7,%%zmm7		\n\t	vfnmadd231pd	%%zmm13,%%zmm14,%%zmm4		\n\t"\
		"vsubpd		%%zmm8 ,%%zmm1,%%zmm1		\n\t	vfnmadd231pd	(%%rcx),%%zmm13,%%zmm3		\n\t"\
		"vmovaps		%%zmm6,0x200(%%rcx)	\n\t		vmovaps		%%zmm2,0x280(%%rcx)	\n\t"\
		"vmovaps		%%zmm0,0x100(%%rcx)	\n\t		vmovaps		%%zmm5,0x180(%%rcx)	\n\t"\
		"vmovaps		%%zmm7,0x240(%%rcx)	\n\t		vmovaps		%%zmm4,0x2c0(%%rcx)	\n\t"\
		"vmovaps		%%zmm1,0x340(%%rcx)	\n\t		vmovaps		%%zmm3,0x3c0(%%rcx)	\n\t"\
	"vmovaps	0x80(%%rdi),%%zmm13 \n\t"/* sqrt2 *//* rcol: *2 => *sqrt2 due to above FMA leaving zmm10,13-15 unmultiplied-by-isrt2: */\
	"vfmadd132pd	(%%rdi),%%zmm6,%%zmm11	\n\t	 vfmadd132pd	%%zmm13,%%zmm2,%%zmm10	\n\t"\
	"vfmadd132pd	(%%rdi),%%zmm0,%%zmm9		\n\t	 vfmadd132pd	%%zmm13,%%zmm5,%%zmm15	\n\t"\
	"vfmadd132pd	(%%rdi),%%zmm7,%%zmm12	\n\t	 vfmadd132pd	%%zmm13,%%zmm4,%%zmm14	\n\t"\
	"vfmadd132pd	(%%rdi),%%zmm1,%%zmm8		\n\t	 vfmadd132pd	(%%rcx),%%zmm3,%%zmm13	\n\t"\
		"vmovaps		%%zmm11,     (%%rcx)	\n\t		vmovaps		%%zmm10,0x080(%%rcx)	\n\t"\
		"vmovaps		%%zmm9,0x300(%%rcx)	\n\t		vmovaps		%%zmm15,0x380(%%rcx)	\n\t"\
		"vmovaps		%%zmm12,0x040(%%rcx)	\n\t		vmovaps		%%zmm14,0x0c0(%%rcx)	\n\t"\
		"vmovaps		%%zmm8,0x140(%%rcx)	\n\t		vmovaps		%%zmm13,0x1c0(%%rcx)	\n\t"\
	/*...Block 2: */\
	/* SSE2_RADIX4_DIF_4WRAPPER(c02,c0A,c12,c1A,r10) SSE2_RADIX4_DIF_4WRAPPER_2NDOFTWO(c06,c0E,c16,c1E,r18) */\
		"movq		%[__r10],%%rcx			\n\t		/*addq		$0x200,%%rcx // __r18 */	\n\t"\
		"movq		%[__c02],%%rdx			\n\t		/*addq		$0x200,%%rdx // __c06 */	\n\t"\
		"vmovaps			 (%%rcx),%%zmm0		\n\t		vmovaps		0x200(%%rcx),%%zmm8			\n\t"\
		"vmovaps		0x040(%%rcx),%%zmm1		\n\t		vmovaps		0x240(%%rcx),%%zmm9			\n\t"\
		"vmovaps		%%zmm0,%%zmm2			\n\t		vmovaps		%%zmm8,%%zmm10				\n\t"\
		"vmovaps		%%zmm1,%%zmm3			\n\t		vmovaps		%%zmm9,%%zmm11				\n\t"\
		"vmulpd			 (%%rdx),%%zmm0,%%zmm0	\n\t		vmulpd		0x200(%%rdx),%%zmm8,%%zmm8	\n\t"\
		"vmulpd			 (%%rdx),%%zmm1,%%zmm1	\n\t		vmulpd		0x200(%%rdx),%%zmm9,%%zmm9	\n\t"\
	"vfnmadd231pd	0x040(%%rdx),%%zmm3,%%zmm0	\n\t	vfnmadd231pd	0x240(%%rdx),%%zmm11,%%zmm8 \n\t"\
	" vfmadd231pd	0x040(%%rdx),%%zmm2,%%zmm1	\n\t	 vfmadd231pd	0x240(%%rdx),%%zmm10,%%zmm9 \n\t"\
		"vmovaps		%%zmm0,%%zmm2			\n\t		vmovaps		%%zmm8,%%zmm10				\n\t"\
		"vmovaps		%%zmm1,%%zmm3			\n\t		vmovaps		%%zmm9,%%zmm11				\n\t"\
		"addq		$0x080,%%rdx																\n\t"\
		"vmovaps		0x080(%%rcx),%%zmm4		\n\t		vmovaps		0x280(%%rcx),%%zmm12		\n\t"\
		"vmovaps		0x0c0(%%rcx),%%zmm5		\n\t		vmovaps		0x2c0(%%rcx),%%zmm13		\n\t"\
		"vmovaps		%%zmm4,%%zmm6			\n\t		vmovaps		%%zmm12,%%zmm14				\n\t"\
		"vmovaps		%%zmm5,%%zmm7			\n\t		vmovaps		%%zmm13,%%zmm15				\n\t"\
		"vmulpd			 (%%rdx),%%zmm4,%%zmm4	\n\t		vmulpd		0x200(%%rdx),%%zmm12,%%zmm12\n\t"\
		"vmulpd			 (%%rdx),%%zmm5,%%zmm5	\n\t		vmulpd		0x200(%%rdx),%%zmm13,%%zmm13\n\t"\
	"vfnmadd231pd	0x040(%%rdx),%%zmm7,%%zmm4	\n\t	vfnmadd231pd	0x240(%%rdx),%%zmm15,%%zmm12\n\t"\
	" vfmadd231pd	0x040(%%rdx),%%zmm6,%%zmm5	\n\t	 vfmadd231pd	0x240(%%rdx),%%zmm14,%%zmm13\n\t"\
		"vaddpd		%%zmm4,%%zmm0,%%zmm0	\n\t		vaddpd		%%zmm12,%%zmm8,%%zmm8	\n\t"\
		"vaddpd		%%zmm5,%%zmm1,%%zmm1	\n\t		vaddpd		%%zmm13,%%zmm9,%%zmm9	\n\t"\
		"vsubpd		%%zmm4,%%zmm2,%%zmm2	\n\t		vsubpd		%%zmm12,%%zmm10,%%zmm10\n\t"\
		"vsubpd		%%zmm5,%%zmm3,%%zmm3	\n\t		vsubpd		%%zmm13,%%zmm11,%%zmm11\n\t"\
		"addq		$0x100,%%rdx																\n\t"\
		"vmovaps		0x180(%%rcx),%%zmm4		\n\t		vmovaps		0x380(%%rcx),%%zmm12		\n\t"\
		"vmovaps		0x1c0(%%rcx),%%zmm5		\n\t		vmovaps		0x3c0(%%rcx),%%zmm13		\n\t"\
		"vmovaps		%%zmm4,%%zmm6			\n\t		vmovaps		%%zmm12,%%zmm14				\n\t"\
		"vmovaps		%%zmm5,%%zmm7			\n\t		vmovaps		%%zmm13,%%zmm15				\n\t"\
		"vmulpd			 (%%rdx),%%zmm4,%%zmm4	\n\t		vmulpd		0x200(%%rdx),%%zmm12,%%zmm12\n\t"\
		"vmulpd			 (%%rdx),%%zmm5,%%zmm5	\n\t		vmulpd		0x200(%%rdx),%%zmm13,%%zmm13\n\t"\
	"vfnmadd231pd	0x040(%%rdx),%%zmm7,%%zmm4	\n\t	vfnmadd231pd	0x240(%%rdx),%%zmm15,%%zmm12\n\t"\
	" vfmadd231pd	0x040(%%rdx),%%zmm6,%%zmm5	\n\t	 vfmadd231pd	0x240(%%rdx),%%zmm14,%%zmm13\n\t"\
		"vmovaps		%%zmm5,0x040(%%rcx)	\n\t		vmovaps		%%zmm13,0x240(%%rcx)		\n\t"\
		"vmovaps		%%zmm4,     (%%rcx)	\n\t		vmovaps		%%zmm12,0x200(%%rcx)		\n\t"\
		"subq		$0x080,%%rdx																\n\t"\
		"vmovaps		0x100(%%rcx),%%zmm4		\n\t		vmovaps		0x300(%%rcx),%%zmm12		\n\t"\
		"vmovaps		0x140(%%rcx),%%zmm5		\n\t		vmovaps		0x340(%%rcx),%%zmm13		\n\t"\
		"vmovaps		%%zmm4,%%zmm6			\n\t		vmovaps		%%zmm12,%%zmm14				\n\t"\
		"vmovaps		%%zmm5,%%zmm7			\n\t		vmovaps		%%zmm13,%%zmm15				\n\t"\
		"vmulpd			 (%%rdx),%%zmm4,%%zmm4	\n\t		vmulpd		0x200(%%rdx),%%zmm12,%%zmm12\n\t"\
		"vmulpd			 (%%rdx),%%zmm5,%%zmm5	\n\t		vmulpd		0x200(%%rdx),%%zmm13,%%zmm13\n\t"\
	"vfnmadd231pd	0x040(%%rdx),%%zmm7,%%zmm4	\n\t	vfnmadd231pd	0x240(%%rdx),%%zmm15,%%zmm12\n\t"\
	" vfmadd231pd	0x040(%%rdx),%%zmm6,%%zmm5	\n\t	 vfmadd231pd	0x240(%%rdx),%%zmm14,%%zmm13\n\t"\
		"vmovaps		%%zmm5,%%zmm7			\n\t		vmovaps		%%zmm13,%%zmm15				\n\t"\
		"vmovaps		%%zmm4,%%zmm6			\n\t		vmovaps		%%zmm12,%%zmm14				\n\t"\
		"vsubpd			 (%%rcx),%%zmm4,%%zmm4	\n\t		vsubpd		0x200(%%rcx),%%zmm12,%%zmm12\n\t"\
		"vsubpd		0x040(%%rcx),%%zmm5,%%zmm5	\n\t		vsubpd		0x240(%%rcx),%%zmm13,%%zmm13\n\t"\
		"vaddpd			 (%%rcx),%%zmm6,%%zmm6	\n\t		vaddpd		0x200(%%rcx),%%zmm14,%%zmm14\n\t"\
		"vaddpd		0x040(%%rcx),%%zmm7,%%zmm7	\n\t		vaddpd		0x240(%%rcx),%%zmm15,%%zmm15\n\t"\
		"vsubpd		%%zmm6,%%zmm0,%%zmm0	\n\t		vsubpd		%%zmm14,%%zmm8,%%zmm8	\n\t"\
		"vsubpd		%%zmm7,%%zmm1,%%zmm1	\n\t		vsubpd		%%zmm15,%%zmm9,%%zmm9	\n\t"\
		"vsubpd		%%zmm5,%%zmm2,%%zmm2	\n\t		vsubpd		%%zmm13,%%zmm10,%%zmm10\n\t"\
		"vsubpd		%%zmm4,%%zmm3,%%zmm3	\n\t		vsubpd		%%zmm12,%%zmm11,%%zmm11\n\t"\
	"vmovaps	%%zmm12,(%%rcx) \n\t"/* spill zmm12 to make room for two */"vmovaps (%%rdi),%%zmm12	\n\t"/* two */\
	"vfmadd132pd	%%zmm12,%%zmm0,%%zmm6		\n\t	vfmadd132pd	%%zmm12,%%zmm8 ,%%zmm14			\n\t"\
	"vfmadd132pd	%%zmm12,%%zmm1,%%zmm7		\n\t	vfmadd132pd	%%zmm12,%%zmm9 ,%%zmm15			\n\t"\
	"vfmadd132pd	%%zmm12,%%zmm2,%%zmm5		\n\t	vfmadd132pd	%%zmm12,%%zmm10,%%zmm13			\n\t"\
	"vfmadd132pd	%%zmm12,%%zmm3,%%zmm4		\n\t	vfmadd132pd	(%%rcx),%%zmm11,%%zmm12			\n\t"\
		"												vmovaps		%%zmm14,0x200(%%rcx)			\n\t"\
		"												vmovaps		%%zmm15,0x240(%%rcx)			\n\t"\
		"												vmovaps		%%zmm10,%%zmm14					\n\t"\
		"												vmovaps		%%zmm13,%%zmm15					\n\t"\
		"												vsubpd		%%zmm12,%%zmm10,%%zmm10	\n\t"\
		"												vsubpd		%%zmm11,%%zmm13,%%zmm13	\n\t"\
		"												vaddpd		%%zmm12,%%zmm14,%%zmm14	\n\t"\
		"												vaddpd		%%zmm11,%%zmm15,%%zmm15	\n\t"\
		/*****	SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r10) *****/\
		"												vmovaps		0x200(%%rcx),%%zmm11			\n\t"\
		"												vmovaps		0x240(%%rcx),%%zmm12			\n\t"\
	"vmovaps %%zmm13,(%%rcx) \n\t"/* spill zmm13 to make room for isrt2 */"vmovaps (%%rsi),%%zmm13	\n\t"/*isrt2*/\
		"vsubpd		%%zmm11,%%zmm6,%%zmm6		\n\t	vfnmadd231pd	%%zmm13,%%zmm10,%%zmm2		\n\t"\
		"vsubpd		%%zmm9,%%zmm0,%%zmm0		\n\t	vfnmadd231pd	%%zmm13,%%zmm15,%%zmm5		\n\t"\
		"vsubpd		%%zmm12,%%zmm7,%%zmm7		\n\t	vfnmadd231pd	%%zmm13,%%zmm14,%%zmm4		\n\t"\
		"vsubpd		%%zmm8 ,%%zmm1,%%zmm1		\n\t	vfnmadd231pd	(%%rcx),%%zmm13,%%zmm3		\n\t"\
		"vmovaps		%%zmm6,0x200(%%rcx)	\n\t		vmovaps		%%zmm2,0x280(%%rcx)	\n\t"\
		"vmovaps		%%zmm0,0x100(%%rcx)	\n\t		vmovaps		%%zmm5,0x180(%%rcx)	\n\t"\
		"vmovaps		%%zmm7,0x240(%%rcx)	\n\t		vmovaps		%%zmm4,0x2c0(%%rcx)	\n\t"\
		"vmovaps		%%zmm1,0x340(%%rcx)	\n\t		vmovaps		%%zmm3,0x3c0(%%rcx)	\n\t"\
	"vmovaps	0x80(%%rdi),%%zmm13 \n\t"/* sqrt2 *//* rcol: *2 => *sqrt2 due to above FMA leaving zmm10,13-15 unmultiplied-by-isrt2: */\
	"vfmadd132pd	(%%rdi),%%zmm6,%%zmm11	\n\t	 vfmadd132pd	%%zmm13,%%zmm2,%%zmm10	\n\t"\
	"vfmadd132pd	(%%rdi),%%zmm0,%%zmm9		\n\t	 vfmadd132pd	%%zmm13,%%zmm5,%%zmm15	\n\t"\
	"vfmadd132pd	(%%rdi),%%zmm7,%%zmm12	\n\t	 vfmadd132pd	%%zmm13,%%zmm4,%%zmm14	\n\t"\
	"vfmadd132pd	(%%rdi),%%zmm1,%%zmm8		\n\t	 vfmadd132pd	(%%rcx),%%zmm3,%%zmm13	\n\t"\
		"vmovaps		%%zmm11,     (%%rcx)	\n\t		vmovaps		%%zmm10,0x080(%%rcx)	\n\t"\
		"vmovaps		%%zmm9,0x300(%%rcx)	\n\t		vmovaps		%%zmm15,0x380(%%rcx)	\n\t"\
		"vmovaps		%%zmm12,0x040(%%rcx)	\n\t		vmovaps		%%zmm14,0x0c0(%%rcx)	\n\t"\
		"vmovaps		%%zmm8,0x140(%%rcx)	\n\t		vmovaps		%%zmm13,0x1c0(%%rcx)	\n\t"\
		"\n\t"\
	/************************************************************************************************************/\
	/* Next 2 blocks operate on odd-indexed elements from the unpck*pd commands which we stored to temporaries: */\
	/************************************************************************************************************/\
	/*...Block 3: */\
	/*	SSE2_RADIX4_DIF_4TWIDDLE	 (r20,r24,r22,r26,r20,c01) */\
		"movq		%[__r20],%%rcx		\n\t"	/***	SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO(r28,r2C,r2A,r2E,r28,c05)	*****/\
		"movq		%[__c01],%%rbx			\n\t	/*	movq		%[__c05],%%rbx	*/		\n\t"\
		"movq		%%rcx,%%rax			\n\t	/*	addq		$0x100,%%rax	*/		\n\t"\
		"addq		$0x080,%%rcx			\n\t	/*	addq		$0x100,%%rcx	*/		\n\t"\
		"vmovaps		 (%%rax),%%zmm0	\n\t movq %%rax,%%rdx \n\t	vmovaps		0x200(%%rax),%%zmm8	\n\t"\
		"vmovaps		 (%%rcx),%%zmm4			\n\t		vmovaps		0x200(%%rcx),%%zmm12		\n\t"\
		"vmovaps	0x040(%%rax),%%zmm1			\n\t		vmovaps		0x240(%%rax),%%zmm9			\n\t"\
		"vmovaps	0x040(%%rcx),%%zmm5			\n\t		vmovaps		0x240(%%rcx),%%zmm13		\n\t"\
		"vmovaps		 (%%rbx),%%zmm6			\n\t		vmovaps		0x200(%%rbx),%%zmm14		\n\t"\
		"vmovaps	0x040(%%rbx),%%zmm7			\n\t		vmovaps		0x240(%%rbx),%%zmm15		\n\t"\
		"vmovaps	%%zmm0,%%zmm2			\n\t		vmovaps		%%zmm8,%%zmm10		\n\t"\
		"vmovaps	%%zmm1,%%zmm3			\n\t		vmovaps		%%zmm9,%%zmm11		\n\t"\
		"vmulpd		%%zmm6,%%zmm0,%%zmm0	\n\t		vmulpd		%%zmm14,%%zmm8,%%zmm8	\n\t"\
		"vmulpd		%%zmm6,%%zmm1,%%zmm1	\n\t		vmulpd		%%zmm14,%%zmm9,%%zmm9	\n\t"\
	"vfnmadd231pd	%%zmm7,%%zmm3,%%zmm0		\n\t	vfnmadd231pd	%%zmm15,%%zmm11,%%zmm8 		\n\t"\
	" vfmadd231pd	%%zmm7,%%zmm2,%%zmm1		\n\t	 vfmadd231pd	%%zmm15,%%zmm10,%%zmm9 		\n\t"\
		"vmovaps		%%zmm4,%%zmm6			\n\t		vmovaps		%%zmm12,%%zmm14		\n\t"\
		"vmovaps		%%zmm5,%%zmm7			\n\t		vmovaps		%%zmm13,%%zmm15		\n\t"\
		"vmulpd		0x080(%%rbx),%%zmm4,%%zmm4	\n\t		vmulpd		0x280(%%rbx),%%zmm12,%%zmm12\n\t"\
		"vmulpd		0x080(%%rbx),%%zmm5,%%zmm5	\n\t		vmulpd		0x280(%%rbx),%%zmm13,%%zmm13\n\t"\
		"vmovaps		%%zmm0,%%zmm2			\n\t		vmovaps		%%zmm8,%%zmm10		\n\t"\
		"vmovaps		%%zmm1,%%zmm3			\n\t		vmovaps		%%zmm9,%%zmm11		\n\t"\
	"vfnmadd231pd	0x0c0(%%rbx),%%zmm7,%%zmm4	\n\t	vfnmadd231pd	0x2c0(%%rbx),%%zmm15,%%zmm12\n\t"\
	" vfmadd231pd	0x0c0(%%rbx),%%zmm6,%%zmm5	\n\t	 vfmadd231pd	0x2c0(%%rbx),%%zmm14,%%zmm13\n\t"\
		"addq		$0x100,%%rcx			\n\t		addq		$0x180,%%rbx			\n\t"\
		"vmovaps			 (%%rcx),%%zmm6		\n\t		vaddpd		%%zmm12,%%zmm8,%%zmm8	\n\t"\
		"vmovaps		0x040(%%rcx),%%zmm7		\n\t		vaddpd		%%zmm13,%%zmm9,%%zmm9	\n\t"\
		"vaddpd		%%zmm4,%%zmm0,%%zmm0	\n\t		vsubpd		%%zmm12,%%zmm10,%%zmm10\n\t"\
		"vaddpd		%%zmm5,%%zmm1,%%zmm1	\n\t		vsubpd		%%zmm13,%%zmm11,%%zmm11\n\t"\
		"vsubpd		%%zmm4,%%zmm2,%%zmm2	\n\t		vmovaps		0x200(%%rcx),%%zmm12		\n\t"\
		"vsubpd		%%zmm5,%%zmm3,%%zmm3	\n\t		vmovaps		0x240(%%rcx),%%zmm13		\n\t"\
		"vmovaps		%%zmm6,%%zmm4			\n\t		vmovaps		0x200(%%rcx),%%zmm14		\n\t"\
		"vmovaps		%%zmm7,%%zmm5			\n\t		vmovaps		0x240(%%rcx),%%zmm15		\n\t"\
		"vmulpd			 (%%rbx),%%zmm4,%%zmm4	\n\t		vmulpd		0x200(%%rbx),%%zmm12,%%zmm12\n\t"\
		"vmulpd			 (%%rbx),%%zmm5,%%zmm5	\n\t		vmulpd		0x200(%%rbx),%%zmm13,%%zmm13\n\t"\
	"vfnmadd231pd	0x040(%%rbx),%%zmm7,%%zmm4	\n\t	vfnmadd231pd	0x240(%%rbx),%%zmm15,%%zmm12\n\t"\
	" vfmadd231pd	0x040(%%rbx),%%zmm6,%%zmm5	\n\t	 vfmadd231pd	0x240(%%rbx),%%zmm14,%%zmm13\n\t"\
		"vmovaps		%%zmm5,0x040(%%rdx)	\n\t		vmovaps		%%zmm13,0x240(%%rdx)		\n\t"\
		"vmovaps		%%zmm4,     (%%rdx)	\n\t		vmovaps		%%zmm12,0x200(%%rdx)		\n\t"\
		"										\n\t		addq	$0x100,%%rax					\n\t"\
		"										\n\t		subq	$0x080,%%rbx					\n\t"\
		"vmovaps			 (%%rax),%%zmm4		\n\t		vmovaps		0x200(%%rax),%%zmm12		\n\t"\
		"vmovaps		0x040(%%rax),%%zmm5		\n\t		vmovaps		0x240(%%rax),%%zmm13		\n\t"\
		"vmovaps		%%zmm4,%%zmm6			\n\t		vmovaps		0x200(%%rax),%%zmm14		\n\t"\
		"vmovaps		%%zmm5,%%zmm7			\n\t		vmovaps		0x240(%%rax),%%zmm15		\n\t"\
		"vmulpd			 (%%rbx),%%zmm4,%%zmm4	\n\t		vmulpd		0x200(%%rbx),%%zmm12,%%zmm12\n\t"\
		"vmulpd			 (%%rbx),%%zmm5,%%zmm5	\n\t		vmulpd		0x200(%%rbx),%%zmm13,%%zmm13\n\t"\
	"vfnmadd231pd	0x040(%%rbx),%%zmm7,%%zmm4	\n\t	vfnmadd231pd	0x240(%%rbx),%%zmm15,%%zmm12\n\t"\
	" vfmadd231pd	0x040(%%rbx),%%zmm6,%%zmm5	\n\t	 vfmadd231pd	0x240(%%rbx),%%zmm14,%%zmm13\n\t"\
		"vmovaps		%%zmm5,%%zmm7			\n\t		vmovaps		%%zmm13,%%zmm15		\n\t"\
		"vmovaps		%%zmm4,%%zmm6			\n\t		vmovaps		%%zmm12,%%zmm14		\n\t"\
		"vsubpd			 (%%rdx),%%zmm4,%%zmm4	\n\t		vsubpd		0x200(%%rdx),%%zmm12,%%zmm12\n\t"\
		"vsubpd		0x040(%%rdx),%%zmm5,%%zmm5	\n\t		vsubpd		0x240(%%rdx),%%zmm13,%%zmm13\n\t"\
		"vaddpd			 (%%rdx),%%zmm6,%%zmm6	\n\t		vaddpd		0x200(%%rdx),%%zmm14,%%zmm14\n\t"\
		"vaddpd		0x040(%%rdx),%%zmm7,%%zmm7	\n\t		vaddpd		0x240(%%rdx),%%zmm15,%%zmm15\n\t"\
		"vsubpd		%%zmm6,%%zmm0,%%zmm0	\n\t		vsubpd		%%zmm14,%%zmm8,%%zmm8	\n\t"\
		"vsubpd		%%zmm7,%%zmm1,%%zmm1	\n\t		vsubpd		%%zmm15,%%zmm9,%%zmm9	\n\t"\
		"vsubpd		%%zmm5,%%zmm2,%%zmm2	\n\t		vsubpd		%%zmm13,%%zmm10,%%zmm10\n\t"\
		"vsubpd		%%zmm4,%%zmm3,%%zmm3	\n\t		vsubpd		%%zmm12,%%zmm11,%%zmm11\n\t"\
	"vmovaps	%%zmm12,(%%rdx) \n\t"/* spill zmm12 to make room for two */"	vmovaps	(%%rdi),%%zmm12 	\n\t"/* two */\
	"vfmadd132pd	%%zmm12,%%zmm0,%%zmm6		\n\t	vfmadd132pd	%%zmm12,%%zmm8 ,%%zmm14			\n\t"\
	"vfmadd132pd	%%zmm12,%%zmm1,%%zmm7		\n\t	vfmadd132pd	%%zmm12,%%zmm9 ,%%zmm15			\n\t"\
	"vfmadd132pd	%%zmm12,%%zmm2,%%zmm5		\n\t	vfmadd132pd	%%zmm12,%%zmm10,%%zmm13			\n\t"\
	"vfmadd132pd	%%zmm12,%%zmm3,%%zmm4		\n\t	vfmadd132pd	(%%rdx),%%zmm11,%%zmm12			\n\t"\
		"												vmovaps		%%zmm14,0x200(%%rdx)			\n\t"\
		"												vmovaps		%%zmm15,0x240(%%rdx)			\n\t"\
		"												vmovaps		%%zmm10,%%zmm14					\n\t"\
		"												vmovaps		%%zmm13,%%zmm15					\n\t"\
		"												vsubpd		%%zmm12,%%zmm10,%%zmm10	\n\t"\
		"												vsubpd		%%zmm11,%%zmm13,%%zmm13	\n\t"\
		"												vaddpd		%%zmm12,%%zmm14,%%zmm14	\n\t"\
		"												vaddpd		%%zmm11,%%zmm15,%%zmm15	\n\t"\
		/*****	SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r20) *****/\
		"												vmovaps		0x200(%%rdx),%%zmm11			\n\t"\
		"												vmovaps		0x240(%%rdx),%%zmm12			\n\t"\
	"vmovaps	%%zmm13,(%%rdx) \n\t"/* spill zmm13 to make room for isrt2 */"	vmovaps	(%%rsi),%%zmm13 	\n\t"/*isrt2*/\
		"vsubpd		%%zmm11,%%zmm6,%%zmm6		\n\t	vfnmadd231pd	%%zmm13,%%zmm10,%%zmm2		\n\t"\
		"vsubpd		%%zmm9,%%zmm0,%%zmm0		\n\t	vfnmadd231pd	%%zmm13,%%zmm15,%%zmm5		\n\t"\
		"vsubpd		%%zmm12,%%zmm7,%%zmm7		\n\t	vfnmadd231pd	%%zmm13,%%zmm14,%%zmm4		\n\t"\
		"vsubpd		%%zmm8 ,%%zmm1,%%zmm1		\n\t	vfnmadd231pd	(%%rdx),%%zmm13,%%zmm3		\n\t"\
		"vmovaps		%%zmm6,0x200(%%rdx)	\n\t		vmovaps		%%zmm2,0x280(%%rdx)	\n\t"\
		"vmovaps		%%zmm0,0x100(%%rdx)	\n\t		vmovaps		%%zmm5,0x180(%%rdx)	\n\t"\
		"vmovaps		%%zmm7,0x240(%%rdx)	\n\t		vmovaps		%%zmm4,0x2c0(%%rdx)	\n\t"\
		"vmovaps		%%zmm1,0x340(%%rdx)	\n\t		vmovaps		%%zmm3,0x3c0(%%rdx)	\n\t"\
	"vmovaps	0x80(%%rdi),%%zmm13 \n\t"/* sqrt2 *//* rcol: *2 => *sqrt2 due to above FMA leaving zmm10,13-15 unmultiplied-by-isrt2: */\
	"vfmadd132pd	(%%rdi),%%zmm6,%%zmm11	\n\t	 vfmadd132pd	%%zmm13,%%zmm2,%%zmm10	\n\t"\
	"vfmadd132pd	(%%rdi),%%zmm0,%%zmm9		\n\t	 vfmadd132pd	%%zmm13,%%zmm5,%%zmm15	\n\t"\
	"vfmadd132pd	(%%rdi),%%zmm7,%%zmm12	\n\t	 vfmadd132pd	%%zmm13,%%zmm4,%%zmm14	\n\t"\
	"vfmadd132pd	(%%rdi),%%zmm1,%%zmm8		\n\t	 vfmadd132pd	(%%rdx),%%zmm3,%%zmm13	\n\t"\
		"vmovaps		%%zmm11,     (%%rdx)	\n\t		vmovaps		%%zmm10,0x080(%%rdx)	\n\t"\
		"vmovaps		%%zmm9,0x300(%%rdx)	\n\t		vmovaps		%%zmm15,0x380(%%rdx)	\n\t"\
		"vmovaps		%%zmm12,0x040(%%rdx)	\n\t		vmovaps		%%zmm14,0x0c0(%%rdx)	\n\t"\
		"vmovaps		%%zmm8,0x140(%%rdx)	\n\t		vmovaps		%%zmm13,0x1c0(%%rdx)	\n\t"\
		"\n\t"\
	/*...Block 4: */\
		"/*****	SSE2_RADIX4_DIF_4TWIDDLE        (r30,r34,r32,r36,r30,c03)	*****/\n\t"\
		"movq		%[__c03],%%rbx			\n\t"\
		"movq		%[__r30],%%rax			\n\t"\
		"movq		%%rax,%%rcx		/*****	SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO(r38,r3C,r3A,r3E,r38,c07)	*****/\n\t"\
		"addq		$0x080,%%rcx			\n\t"\
		"vmovaps		 (%%rax),%%zmm0	\n\t movq %%rax,%%rdx \n\t	vmovaps		0x200(%%rax),%%zmm8	\n\t"\
		"vmovaps		 (%%rcx),%%zmm4			\n\t		vmovaps		0x200(%%rcx),%%zmm12		\n\t"\
		"vmovaps	0x040(%%rax),%%zmm1			\n\t		vmovaps		0x240(%%rax),%%zmm9			\n\t"\
		"vmovaps	0x040(%%rcx),%%zmm5			\n\t		vmovaps		0x240(%%rcx),%%zmm13		\n\t"\
		"vmovaps		 (%%rbx),%%zmm6			\n\t		vmovaps		0x200(%%rbx),%%zmm14		\n\t"\
		"vmovaps	0x040(%%rbx),%%zmm7			\n\t		vmovaps		0x240(%%rbx),%%zmm15		\n\t"\
		"vmovaps	%%zmm0,%%zmm2			\n\t		vmovaps		%%zmm8,%%zmm10		\n\t"\
		"vmovaps	%%zmm1,%%zmm3			\n\t		vmovaps		%%zmm9,%%zmm11		\n\t"\
		"vmulpd		%%zmm6,%%zmm0,%%zmm0	\n\t		vmulpd		%%zmm14,%%zmm8,%%zmm8	\n\t"\
		"vmulpd		%%zmm6,%%zmm1,%%zmm1	\n\t		vmulpd		%%zmm14,%%zmm9,%%zmm9	\n\t"\
	"vfnmadd231pd	%%zmm7,%%zmm3,%%zmm0		\n\t	vfnmadd231pd	%%zmm15,%%zmm11,%%zmm8 		\n\t"\
	" vfmadd231pd	%%zmm7,%%zmm2,%%zmm1		\n\t	 vfmadd231pd	%%zmm15,%%zmm10,%%zmm9 		\n\t"\
		"vmovaps		%%zmm4,%%zmm6			\n\t		vmovaps		%%zmm12,%%zmm14		\n\t"\
		"vmovaps		%%zmm5,%%zmm7			\n\t		vmovaps		%%zmm13,%%zmm15		\n\t"\
		"vmulpd		0x080(%%rbx),%%zmm4,%%zmm4	\n\t		vmulpd		0x280(%%rbx),%%zmm12,%%zmm12\n\t"\
		"vmulpd		0x080(%%rbx),%%zmm5,%%zmm5	\n\t		vmulpd		0x280(%%rbx),%%zmm13,%%zmm13\n\t"\
		"vmovaps		%%zmm0,%%zmm2			\n\t		vmovaps		%%zmm8,%%zmm10		\n\t"\
		"vmovaps		%%zmm1,%%zmm3			\n\t		vmovaps		%%zmm9,%%zmm11		\n\t"\
	"vfnmadd231pd	0x0c0(%%rbx),%%zmm7,%%zmm4	\n\t	vfnmadd231pd	0x2c0(%%rbx),%%zmm15,%%zmm12\n\t"\
	" vfmadd231pd	0x0c0(%%rbx),%%zmm6,%%zmm5	\n\t	 vfmadd231pd	0x2c0(%%rbx),%%zmm14,%%zmm13\n\t"\
		"addq		$0x100,%%rcx			\n\t		addq		$0x180,%%rbx			\n\t"\
		"vmovaps			 (%%rcx),%%zmm6		\n\t		vaddpd		%%zmm12,%%zmm8,%%zmm8	\n\t"\
		"vmovaps		0x040(%%rcx),%%zmm7		\n\t		vaddpd		%%zmm13,%%zmm9,%%zmm9	\n\t"\
		"vaddpd		%%zmm4,%%zmm0,%%zmm0	\n\t		vsubpd		%%zmm12,%%zmm10,%%zmm10\n\t"\
		"vaddpd		%%zmm5,%%zmm1,%%zmm1	\n\t		vsubpd		%%zmm13,%%zmm11,%%zmm11\n\t"\
		"vsubpd		%%zmm4,%%zmm2,%%zmm2	\n\t		vmovaps		0x200(%%rcx),%%zmm12		\n\t"\
		"vsubpd		%%zmm5,%%zmm3,%%zmm3	\n\t		vmovaps		0x240(%%rcx),%%zmm13		\n\t"\
		"vmovaps		%%zmm6,%%zmm4			\n\t		vmovaps		0x200(%%rcx),%%zmm14		\n\t"\
		"vmovaps		%%zmm7,%%zmm5			\n\t		vmovaps		0x240(%%rcx),%%zmm15		\n\t"\
		"vmulpd			 (%%rbx),%%zmm4,%%zmm4	\n\t		vmulpd		0x200(%%rbx),%%zmm12,%%zmm12\n\t"\
		"vmulpd			 (%%rbx),%%zmm5,%%zmm5	\n\t		vmulpd		0x200(%%rbx),%%zmm13,%%zmm13\n\t"\
	"vfnmadd231pd	0x040(%%rbx),%%zmm7,%%zmm4	\n\t	vfnmadd231pd	0x240(%%rbx),%%zmm15,%%zmm12\n\t"\
	" vfmadd231pd	0x040(%%rbx),%%zmm6,%%zmm5	\n\t	 vfmadd231pd	0x240(%%rbx),%%zmm14,%%zmm13\n\t"\
		"vmovaps		%%zmm5,0x040(%%rdx)	\n\t		vmovaps		%%zmm13,0x240(%%rdx)		\n\t"\
		"vmovaps		%%zmm4,     (%%rdx)	\n\t		vmovaps		%%zmm12,0x200(%%rdx)		\n\t"\
		"										\n\t		addq	$0x100,%%rax					\n\t"\
		"										\n\t		subq	$0x080,%%rbx					\n\t"\
		"vmovaps			 (%%rax),%%zmm4		\n\t		vmovaps		0x200(%%rax),%%zmm12		\n\t"\
		"vmovaps		0x040(%%rax),%%zmm5		\n\t		vmovaps		0x240(%%rax),%%zmm13		\n\t"\
		"vmovaps		%%zmm4,%%zmm6			\n\t		vmovaps		0x200(%%rax),%%zmm14		\n\t"\
		"vmovaps		%%zmm5,%%zmm7			\n\t		vmovaps		0x240(%%rax),%%zmm15		\n\t"\
		"vmulpd			 (%%rbx),%%zmm4,%%zmm4	\n\t		vmulpd		0x200(%%rbx),%%zmm12,%%zmm12\n\t"\
		"vmulpd			 (%%rbx),%%zmm5,%%zmm5	\n\t		vmulpd		0x200(%%rbx),%%zmm13,%%zmm13\n\t"\
	"vfnmadd231pd	0x040(%%rbx),%%zmm7,%%zmm4	\n\t	vfnmadd231pd	0x240(%%rbx),%%zmm15,%%zmm12\n\t"\
	" vfmadd231pd	0x040(%%rbx),%%zmm6,%%zmm5	\n\t	 vfmadd231pd	0x240(%%rbx),%%zmm14,%%zmm13\n\t"\
		"vmovaps		%%zmm5,%%zmm7			\n\t		vmovaps		%%zmm13,%%zmm15		\n\t"\
		"vmovaps		%%zmm4,%%zmm6			\n\t		vmovaps		%%zmm12,%%zmm14		\n\t"\
		"vsubpd			 (%%rdx),%%zmm4,%%zmm4	\n\t		vsubpd		0x200(%%rdx),%%zmm12,%%zmm12\n\t"\
		"vsubpd		0x040(%%rdx),%%zmm5,%%zmm5	\n\t		vsubpd		0x240(%%rdx),%%zmm13,%%zmm13\n\t"\
		"vaddpd			 (%%rdx),%%zmm6,%%zmm6	\n\t		vaddpd		0x200(%%rdx),%%zmm14,%%zmm14\n\t"\
		"vaddpd		0x040(%%rdx),%%zmm7,%%zmm7	\n\t		vaddpd		0x240(%%rdx),%%zmm15,%%zmm15\n\t"\
		"vsubpd		%%zmm6,%%zmm0,%%zmm0	\n\t		vsubpd		%%zmm14,%%zmm8,%%zmm8	\n\t"\
		"vsubpd		%%zmm7,%%zmm1,%%zmm1	\n\t		vsubpd		%%zmm15,%%zmm9,%%zmm9	\n\t"\
		"vsubpd		%%zmm5,%%zmm2,%%zmm2	\n\t		vsubpd		%%zmm13,%%zmm10,%%zmm10\n\t"\
		"vsubpd		%%zmm4,%%zmm3,%%zmm3	\n\t		vsubpd		%%zmm12,%%zmm11,%%zmm11\n\t"\
	"vmovaps	%%zmm12,(%%rdx) \n\t"/* spill zmm12 to make room for two */"	vmovaps	(%%rdi),%%zmm12 	\n\t"/* two */\
	"vfmadd132pd	%%zmm12,%%zmm0,%%zmm6		\n\t	vfmadd132pd	%%zmm12,%%zmm8 ,%%zmm14			\n\t"\
	"vfmadd132pd	%%zmm12,%%zmm1,%%zmm7		\n\t	vfmadd132pd	%%zmm12,%%zmm9 ,%%zmm15			\n\t"\
	"vfmadd132pd	%%zmm12,%%zmm2,%%zmm5		\n\t	vfmadd132pd	%%zmm12,%%zmm10,%%zmm13			\n\t"\
	"vfmadd132pd	%%zmm12,%%zmm3,%%zmm4		\n\t	vfmadd132pd	(%%rdx),%%zmm11,%%zmm12			\n\t"\
		"												vmovaps		%%zmm14,0x200(%%rdx)			\n\t"\
		"												vmovaps		%%zmm15,0x240(%%rdx)			\n\t"\
		"												vmovaps		%%zmm10,%%zmm14					\n\t"\
		"												vmovaps		%%zmm13,%%zmm15					\n\t"\
		"												vsubpd		%%zmm12,%%zmm10,%%zmm10	\n\t"\
		"												vsubpd		%%zmm11,%%zmm13,%%zmm13	\n\t"\
		"												vaddpd		%%zmm12,%%zmm14,%%zmm14	\n\t"\
		"												vaddpd		%%zmm11,%%zmm15,%%zmm15	\n\t"\
		/*****	SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r30) *****/\
		"												vmovaps		0x200(%%rdx),%%zmm11			\n\t"\
		"												vmovaps		0x240(%%rdx),%%zmm12			\n\t"\
	"vmovaps	%%zmm13,(%%rdx) \n\t"/* spill zmm13 to make room for isrt2 */"	vmovaps	(%%rsi),%%zmm13 	\n\t"/*isrt2*/\
		"vsubpd		%%zmm11,%%zmm6,%%zmm6		\n\t	vfnmadd231pd	%%zmm13,%%zmm10,%%zmm2		\n\t"\
		"vsubpd		%%zmm9,%%zmm0,%%zmm0		\n\t	vfnmadd231pd	%%zmm13,%%zmm15,%%zmm5		\n\t"\
		"vsubpd		%%zmm12,%%zmm7,%%zmm7		\n\t	vfnmadd231pd	%%zmm13,%%zmm14,%%zmm4		\n\t"\
		"vsubpd		%%zmm8 ,%%zmm1,%%zmm1		\n\t	vfnmadd231pd	(%%rdx),%%zmm13,%%zmm3		\n\t"\
		"vmovaps		%%zmm6,0x200(%%rdx)	\n\t		vmovaps		%%zmm2,0x280(%%rdx)	\n\t"\
		"vmovaps		%%zmm0,0x100(%%rdx)	\n\t		vmovaps		%%zmm5,0x180(%%rdx)	\n\t"\
		"vmovaps		%%zmm7,0x240(%%rdx)	\n\t		vmovaps		%%zmm4,0x2c0(%%rdx)	\n\t"\
		"vmovaps		%%zmm1,0x340(%%rdx)	\n\t		vmovaps		%%zmm3,0x3c0(%%rdx)	\n\t"\
	"vmovaps	0x80(%%rdi),%%zmm13 \n\t"/* sqrt2 *//* rcol: *2 => *sqrt2 due to above FMA leaving zmm10,13-15 unmultiplied-by-isrt2: */\
	"vfmadd132pd	(%%rdi),%%zmm6,%%zmm11	\n\t	 vfmadd132pd	%%zmm13,%%zmm2,%%zmm10	\n\t"\
	"vfmadd132pd	(%%rdi),%%zmm0,%%zmm9		\n\t	 vfmadd132pd	%%zmm13,%%zmm5,%%zmm15	\n\t"\
	"vfmadd132pd	(%%rdi),%%zmm7,%%zmm12	\n\t	 vfmadd132pd	%%zmm13,%%zmm4,%%zmm14	\n\t"\
	"vfmadd132pd	(%%rdi),%%zmm1,%%zmm8		\n\t	 vfmadd132pd	(%%rdx),%%zmm3,%%zmm13	\n\t"\
		"vmovaps		%%zmm11,     (%%rdx)	\n\t		vmovaps		%%zmm10,0x080(%%rdx)	\n\t"\
		"vmovaps		%%zmm9,0x300(%%rdx)	\n\t		vmovaps		%%zmm15,0x380(%%rdx)	\n\t"\
		"vmovaps		%%zmm12,0x040(%%rdx)	\n\t		vmovaps		%%zmm14,0x0c0(%%rdx)	\n\t"\
		"vmovaps		%%zmm8,0x140(%%rdx)	\n\t		vmovaps		%%zmm13,0x1c0(%%rdx)	\n\t"\
		"\n\t"\
	/**********************************************************************************/\
	/*...and now do eight radix-4 transforms, including the internal twiddle factors: */\
	/**********************************************************************************/\
		"movq		%[__isrt2],%%rsi		\n\t"\
		/*...Block 1: t00,t10,t20,t30	*/				/*...Block 5: t08,t18,t28,t38	*/\
		"movq		%[__r00],%%rax			\n\t"\
		"movq		%[__r10],%%rbx			\n\t"\
		"movq		%[__r20],%%rcx			\n\t"\
		"movq		%[__r30],%%rdx			\n\t"\
		"vmovaps		 (%%rax),%%zmm0			\n\t		vmovaps		0x200(%%rax),%%zmm8			\n\t"\
		"vmovaps	0x040(%%rax),%%zmm1			\n\t		vmovaps		0x240(%%rax),%%zmm9			\n\t"\
		"vmovaps		 (%%rbx),%%zmm2			\n\t		vmovaps		0x200(%%rbx),%%zmm10		\n\t"\
		"vmovaps	0x040(%%rbx),%%zmm3			\n\t		vmovaps		0x240(%%rbx),%%zmm11		\n\t"\
		"vaddpd			%%zmm0,%%zmm2,%%zmm2	\n\t		vaddpd		%%zmm9 ,%%zmm10,%%zmm10		\n\t"\
		"vaddpd			%%zmm1,%%zmm3,%%zmm3	\n\t		vaddpd		%%zmm8 ,%%zmm11,%%zmm11		\n\t"\
		"vsubpd			 (%%rbx),%%zmm0,%%zmm0	\n\t		vsubpd		0x240(%%rbx),%%zmm8,%%zmm8	\n\t"\
		"vsubpd		0x040(%%rbx),%%zmm1,%%zmm1	\n\t		vsubpd		0x200(%%rbx),%%zmm9,%%zmm9	\n\t"\
		"vmovaps		 (%%rcx),%%zmm4			\n\t		vmovaps		0x200(%%rcx),%%zmm12		\n\t"\
		"vmovaps	0x040(%%rcx),%%zmm5			\n\t		vmovaps		0x240(%%rcx),%%zmm13		\n\t"\
		"vmovaps		 (%%rdx),%%zmm6			\n\t		vmovaps		0x200(%%rdx),%%zmm14		\n\t"\
		"vmovaps	0x040(%%rdx),%%zmm7			\n\t		vmovaps		0x240(%%rdx),%%zmm15		\n\t"\
		"vaddpd		%%zmm4,%%zmm6,%%zmm6		\n\t		vsubpd		%%zmm13,%%zmm12,%%zmm12		\n\t"\
		"vaddpd		%%zmm5,%%zmm7,%%zmm7		\n\t		vaddpd		0x200(%%rcx),%%zmm13,%%zmm13\n\t"\
		"vsubpd			 (%%rdx),%%zmm4,%%zmm4	\n\t		vaddpd			 %%zmm15,%%zmm14,%%zmm14\n\t"\
		"vsubpd		0x040(%%rdx),%%zmm5,%%zmm5	\n\t		vsubpd		0x200(%%rdx),%%zmm15,%%zmm15\n\t"\
		"vsubpd		%%zmm6,%%zmm2,%%zmm2	\n\t		vmulpd		(%%rsi),%%zmm12,%%zmm12\n\t"/* isrt2 */\
		"vsubpd		%%zmm7,%%zmm3,%%zmm3	\n\t		vmulpd		(%%rsi),%%zmm13,%%zmm13\n\t"\
	"vfmadd132pd	(%%rdi),%%zmm2,%%zmm6	\n\t	vfnmadd231pd	(%%rsi),%%zmm14,%%zmm12\n\t"\
	"vfmadd132pd	(%%rdi),%%zmm3,%%zmm7	\n\t	vfnmadd231pd	(%%rsi),%%zmm15,%%zmm13\n\t"\
		"vmovaps		%%zmm2,      (%%rcx)	\n\t	 vfmadd132pd	0x80(%%rdi),%%zmm12,%%zmm14\n\t"/* sqrt2 */\
		"vmovaps		%%zmm3, 0x040(%%rcx)	\n\t	 vfmadd132pd	0x80(%%rdi),%%zmm13,%%zmm15\n\t"\
		"										\n\t		vsubpd		%%zmm12,%%zmm8,%%zmm8	\n\t"\
		"										\n\t		vsubpd		%%zmm13,%%zmm10,%%zmm10\n\t"\
		"vmovaps		%%zmm6,      (%%rax)	\n\t	 vfmadd132pd	(%%rdi),%%zmm8,%%zmm12\n\t"\
		"vmovaps		%%zmm7, 0x040(%%rax)	\n\t	 vfmadd132pd	(%%rdi),%%zmm10,%%zmm13\n\t"\
		"vsubpd		%%zmm5,%%zmm0,%%zmm0	\n\t		vmovaps		%%zmm8 ,0x200(%%rcx)		\n\t"\
		"vsubpd		%%zmm4,%%zmm1,%%zmm1	\n\t		vmovaps		%%zmm10,0x240(%%rcx)		\n\t"\
	"vfmadd132pd	(%%rdi),%%zmm0,%%zmm5	\n\t		vmovaps		%%zmm12,0x200(%%rax)		\n\t"\
	"vfmadd132pd	(%%rdi),%%zmm1,%%zmm4	\n\t		vmovaps		%%zmm13,0x240(%%rax)		\n\t"\
		"vmovaps		%%zmm0,      (%%rbx)	\n\t		vsubpd		%%zmm15,%%zmm11,%%zmm11\n\t"\
		"vmovaps		%%zmm1, 0x040(%%rdx)	\n\t		vsubpd		%%zmm14,%%zmm9,%%zmm9	\n\t"\
		"										\n\t	 vfmadd132pd	(%%rdi),%%zmm11,%%zmm15\n\t"\
		"										\n\t	 vfmadd132pd	(%%rdi),%%zmm9 ,%%zmm14\n\t"\
		"vmovaps		%%zmm5,      (%%rdx)	\n\t		vmovaps		%%zmm11,0x200(%%rbx)		\n\t"\
		"vmovaps		%%zmm4, 0x040(%%rbx)	\n\t		vmovaps		%%zmm9 ,0x240(%%rdx)		\n\t"\
		"													vmovaps		%%zmm15,0x200(%%rdx)		\n\t"\
		"													vmovaps		%%zmm14,0x240(%%rbx)		\n\t"\
		/*...Block 3: t04,t14,t24,t34	*/			/*...Block 7: t0C,t1C,t2C,t3C	*/\
		"addq		$0x100,%%rax			\n\t"\
		"addq		$0x100,%%rbx			\n\t"\
		"addq		$0x100,%%rcx			\n\t"\
		"addq		$0x100,%%rdx			\n\t"\
	/* Instead of loading regs 0-2, 8-10 and copying to 4-6, 12-14, use the later MULs to effect copying: */\
		"vmovaps		 (%%rcx),%%zmm4			\n\t		vmovaps		0x200(%%rcx),%%zmm12		\n\t"\
		"vmovaps	0x040(%%rcx),%%zmm5			\n\t		vmovaps		0x240(%%rcx),%%zmm13		\n\t"\
		"vmovaps		 (%%rdx),%%zmm6			\n\t		vmovaps		0x200(%%rdx),%%zmm14		\n\t"\
		"vmovaps	0x040(%%rdx),%%zmm3			\n\t		vmovaps		0x240(%%rdx),%%zmm11		\n\t"\
	/* Strategic register-usage for sincos consts here cuts loads from 16 to 4: */\
		"vmovaps		0x080(%%rsi),%%zmm7 	\n\t		vmovaps		0x040(%%rsi),%%zmm15		\n\t"/* s,c */\
		"vmulpd		%%zmm7 ,%%zmm5,%%zmm1		\n\t		vmulpd		 %%zmm15,%%zmm13,%%zmm9		\n\t"\
		"vmulpd		%%zmm15,%%zmm3,%%zmm3		\n\t		vmulpd		 %%zmm7 ,%%zmm11,%%zmm11	\n\t"\
		"vmulpd		%%zmm7 ,%%zmm4,%%zmm0		\n\t		vmulpd		 %%zmm15,%%zmm12,%%zmm8		\n\t"\
		"vmulpd		%%zmm15,%%zmm6,%%zmm2		\n\t		vmulpd		 %%zmm7 ,%%zmm14,%%zmm10	\n\t"\
	"vfmsub132pd	%%zmm15,%%zmm1,%%zmm4		\n\t	vfmsub132pd		 %%zmm7 ,%%zmm9 ,%%zmm12	\n\t"\
	"vfmsub132pd	%%zmm7 ,%%zmm3,%%zmm6		\n\t	vfmsub132pd		 %%zmm15,%%zmm11,%%zmm14	\n\t"\
	"vfmadd132pd	%%zmm15,%%zmm0,%%zmm5		\n\t	vfmadd132pd		 %%zmm7 ,%%zmm8 ,%%zmm13	\n\t"\
	"vfmadd132pd 0x40(%%rdx),%%zmm2,%%zmm7		\n\t	vfmadd132pd 0x240(%%rdx),%%zmm10,%%zmm15	\n\t"\
		"vsubpd		%%zmm6,%%zmm4,%%zmm4	\n\t		vsubpd		%%zmm14,%%zmm12,%%zmm12\n\t"\
		"vsubpd		%%zmm7,%%zmm5,%%zmm5	\n\t		vsubpd		%%zmm15,%%zmm13,%%zmm13\n\t"\
		"vmovaps		 (%%rbx),%%zmm2			\n\t		vmovaps		0x200(%%rbx),%%zmm10		\n\t"\
		"vmovaps	0x040(%%rbx),%%zmm3			\n\t		vmovaps		0x240(%%rbx),%%zmm11		\n\t"\
	"vfmadd132pd	(%%rdi),%%zmm4,%%zmm6	\n\t	vfmadd132pd		(%%rdi),%%zmm12,%%zmm14\n\t"\
	"vfmadd132pd	(%%rdi),%%zmm5,%%zmm7	\n\t	vfmadd132pd		(%%rdi),%%zmm13,%%zmm15\n\t"\
		"vsubpd		0x040(%%rbx),%%zmm2,%%zmm2	\n\t		vaddpd		0x240(%%rbx),%%zmm10,%%zmm10\n\t"\
		"vaddpd			 (%%rbx),%%zmm3,%%zmm3	\n\t		vsubpd		0x200(%%rbx),%%zmm11,%%zmm11\n\t"\
	"vmovaps %%zmm6 ,(%%rcx) \n\t"/* spill zmm6 to make room for isrt2 */"vmovaps (%%rsi),%%zmm6 	\n\t"/* isrt2 */\
		"vmovaps		 (%%rax),%%zmm0			\n\t		vmovaps		 0x200(%%rax),%%zmm8		\n\t"\
		"vmovaps	0x040(%%rax),%%zmm1			\n\t		vmovaps		 0x240(%%rax),%%zmm9		\n\t"\
	"vfnmadd231pd	%%zmm6 ,%%zmm2,%%zmm0		\n\t	vfnmadd231pd	%%zmm6 ,%%zmm10,%%zmm8		\n\t"\
	"vfnmadd231pd	%%zmm6 ,%%zmm3,%%zmm1		\n\t	vfnmadd231pd	%%zmm6 ,%%zmm11,%%zmm9		\n\t"\
	" vfmadd213pd		 (%%rax),%%zmm6 ,%%zmm2	\n\t	 vfmadd213pd	0x200(%%rax),%%zmm6 ,%%zmm10\n\t"\
	" vfmadd213pd	0x040(%%rax),%%zmm6 ,%%zmm3	\n\t	 vfmadd213pd	0x240(%%rax),%%zmm6 ,%%zmm11\n\t"\
	"vmovaps	(%%rcx),%%zmm6 	\n\t"/* restore spill */\
		"vsubpd		%%zmm6,%%zmm2,%%zmm2	\n\t		vsubpd		%%zmm12,%%zmm8,%%zmm8	\n\t"\
		"vsubpd		%%zmm7,%%zmm3,%%zmm3	\n\t		vsubpd		%%zmm13,%%zmm9,%%zmm9	\n\t"\
	"vfmadd132pd	(%%rdi),%%zmm2,%%zmm6	\n\t	vfmadd132pd		(%%rdi),%%zmm8 ,%%zmm12\n\t"\
	"vfmadd132pd	(%%rdi),%%zmm3,%%zmm7	\n\t	vfmadd132pd		(%%rdi),%%zmm9 ,%%zmm13\n\t"\
		"vmovaps		%%zmm2,      (%%rcx)	\n\t		vmovaps		%%zmm8 ,0x200(%%rcx)		\n\t"\
		"vmovaps		%%zmm3, 0x040(%%rcx)	\n\t		vmovaps		%%zmm9 ,0x240(%%rcx)		\n\t"\
		"vmovaps		%%zmm6,      (%%rax)	\n\t		vmovaps		%%zmm12,0x200(%%rax)		\n\t"\
		"vmovaps		%%zmm7, 0x040(%%rax)	\n\t		vmovaps		%%zmm13,0x240(%%rax)		\n\t"\
		"vsubpd		%%zmm5,%%zmm0,%%zmm0	\n\t		vsubpd		%%zmm15,%%zmm10,%%zmm10\n\t"\
		"vsubpd		%%zmm4,%%zmm1,%%zmm1	\n\t		vsubpd		%%zmm14,%%zmm11,%%zmm11\n\t"\
	"vfmadd132pd	(%%rdi),%%zmm0,%%zmm5	\n\t	vfmadd132pd		(%%rdi),%%zmm10,%%zmm15\n\t"\
	"vfmadd132pd	(%%rdi),%%zmm1,%%zmm4	\n\t	vfmadd132pd		(%%rdi),%%zmm11,%%zmm14\n\t"\
		"vmovaps		%%zmm0,      (%%rbx)	\n\t		vmovaps		%%zmm10,0x200(%%rbx)		\n\t"\
		"vmovaps		%%zmm1, 0x040(%%rdx)	\n\t		vmovaps		%%zmm11,0x240(%%rdx)		\n\t"\
		"vmovaps		%%zmm5,      (%%rdx)	\n\t		vmovaps		%%zmm15,0x200(%%rdx)		\n\t"\
		"vmovaps		%%zmm4, 0x040(%%rbx)	\n\t		vmovaps		%%zmm14,0x240(%%rbx)		\n\t"\
		/*...Block 2: t02,t12,t22,t32	*/				/*...Block 6: t0A,t1A,t2A,t3A */\
		"subq		$0x080,%%rax			\n\t"\
		"subq		$0x080,%%rbx			\n\t"\
		"subq		$0x080,%%rcx			\n\t"\
		"subq		$0x080,%%rdx			\n\t"\
		"addq		$0x0c0,%%rsi	\n\t"/* cc1 */\
	/* Instead of loading regs 0-2, 8 and copying to 4-6, 12, use the later MULs to effect copying: */\
		"vmovaps		 (%%rcx),%%zmm4			\n\t		vmovaps		0x200(%%rcx),%%zmm12		\n\t"\
		"vmovaps		 (%%rdx),%%zmm6			\n\t		vmovaps		0x200(%%rdx),%%zmm14		\n\t"\
		"vmovaps	0x040(%%rcx),%%zmm5			\n\t		vmovaps		0x240(%%rcx),%%zmm9			\n\t"\
		"vmovaps	0x040(%%rdx),%%zmm3			\n\t		vmovaps		0x240(%%rdx),%%zmm11		\n\t"\
	/* Strategic register-usage for sincos consts here cuts loads from 16 to 4: */\
		"vmovaps		 (%%rsi),%%zmm15		\n\t		vmovaps		0x040(%%rsi),%%zmm10		\n\t"\
		"vmovaps	0x080(%%rsi),%%zmm7 		\n\t		vmovaps		0x0c0(%%rsi),%%zmm13		\n\t"\
		"vmulpd		%%zmm10,%%zmm5,%%zmm1		\n\t		vmulpd			%%zmm7 ,%%zmm9 ,%%zmm9	\n\t"\
		"vmulpd		%%zmm13,%%zmm3,%%zmm3		\n\t		vmulpd			%%zmm10,%%zmm11,%%zmm11	\n\t"\
		"vmulpd		%%zmm10,%%zmm4,%%zmm0		\n\t		vmulpd			%%zmm7 ,%%zmm12,%%zmm8	\n\t"\
		"vmulpd		%%zmm13,%%zmm6,%%zmm2		\n\t		vmulpd		0x200(%%rdx),%%zmm10,%%zmm10\n\t"\
	"vfmsub132pd	%%zmm15,%%zmm1,%%zmm4		\n\t	vfmsub132pd			%%zmm13,%%zmm9 ,%%zmm12	\n\t"\
	"vfmsub132pd	%%zmm7 ,%%zmm3,%%zmm6		\n\t	vfmadd132pd			%%zmm15,%%zmm11,%%zmm14	\n\t"\
	"vfmadd132pd	%%zmm15,%%zmm0,%%zmm5		\n\t	vfmadd132pd		0x240(%%rcx),%%zmm8 ,%%zmm13\n\t"\
	"vfmadd132pd	0x040(%%rdx),%%zmm2,%%zmm7	\n\t	vfmsub132pd		0x240(%%rdx),%%zmm10,%%zmm15\n\t"\
		"vsubpd		%%zmm6,%%zmm4,%%zmm4	\n\t		vsubpd		%%zmm14,%%zmm12,%%zmm12\n\t"\
		"vsubpd		%%zmm7,%%zmm5,%%zmm5	\n\t		vsubpd		%%zmm15,%%zmm13,%%zmm13\n\t"\
		"vmovaps		 (%%rbx),%%zmm2			\n\t		vmovaps		0x200(%%rbx),%%zmm10		\n\t"\
		"vmovaps	0x040(%%rbx),%%zmm0			\n\t		vmovaps		0x240(%%rbx),%%zmm8			\n\t"\
	"vfmadd132pd	(%%rdi),%%zmm4,%%zmm6	\n\t	vfmadd132pd		(%%rdi),%%zmm12,%%zmm14\n\t"\
	"vfmadd132pd	(%%rdi),%%zmm5,%%zmm7	\n\t	vfmadd132pd		(%%rdi),%%zmm13,%%zmm15\n\t"\
		"vmovaps		 (%%rbx),%%zmm1			\n\t		vmovaps		0x200(%%rbx),%%zmm9			\n\t"\
		"vmovaps	-0x80(%%rsi),%%zmm3			\n\t		vmovaps		-0x40(%%rsi),%%zmm11		\n\t"\
		"vmulpd		%%zmm11,%%zmm0,%%zmm0		\n\t		vmulpd			%%zmm3 ,%%zmm8 ,%%zmm8	\n\t"\
		"vmulpd		%%zmm11,%%zmm1,%%zmm1		\n\t		vmulpd			%%zmm3 ,%%zmm9 ,%%zmm9	\n\t"\
	"vfmsub132pd	%%zmm3 ,%%zmm0,%%zmm2		\n\t	vfmadd132pd			%%zmm11,%%zmm8 ,%%zmm10	\n\t"\
	"vfmadd132pd	0x40(%%rbx),%%zmm1,%%zmm3	\n\t	vfmsub132pd		0x240(%%rbx),%%zmm9,%%zmm11	\n\t"\
		"vmovaps		 (%%rax),%%zmm0			\n\t		vmovaps		0x200(%%rax),%%zmm8			\n\t"\
		"vmovaps	0x040(%%rax),%%zmm1			\n\t		vmovaps		0x240(%%rax),%%zmm9			\n\t"\
		"vsubpd		%%zmm2,%%zmm0,%%zmm0	\n\t		vsubpd		%%zmm10,%%zmm8,%%zmm8	\n\t"\
		"vsubpd		%%zmm3,%%zmm1,%%zmm1	\n\t		vsubpd		%%zmm11,%%zmm9,%%zmm9	\n\t"\
		"vaddpd			 (%%rax),%%zmm2,%%zmm2	\n\t		vaddpd		0x200(%%rax),%%zmm10,%%zmm10\n\t"\
		"vaddpd		0x040(%%rax),%%zmm3,%%zmm3	\n\t		vaddpd		0x240(%%rax),%%zmm11,%%zmm11\n\t"\
		"vsubpd		%%zmm6,%%zmm2,%%zmm2	\n\t		vsubpd		%%zmm12,%%zmm8,%%zmm8	\n\t"\
		"vsubpd		%%zmm7,%%zmm3,%%zmm3	\n\t		vsubpd		%%zmm13,%%zmm9,%%zmm9	\n\t"\
	"vfmadd132pd	(%%rdi),%%zmm2,%%zmm6	\n\t	vfmadd132pd		(%%rdi),%%zmm8,%%zmm12\n\t"\
	"vfmadd132pd	(%%rdi),%%zmm3,%%zmm7	\n\t	vfmadd132pd		(%%rdi),%%zmm9,%%zmm13\n\t"\
		"vmovaps		%%zmm2,      (%%rcx)	\n\t		vmovaps		%%zmm8 ,0x200(%%rcx)		\n\t"\
		"vmovaps		%%zmm3, 0x040(%%rcx)	\n\t		vmovaps		%%zmm9 ,0x240(%%rcx)		\n\t"\
		"vmovaps		%%zmm6,      (%%rax)	\n\t		vmovaps		%%zmm12,0x200(%%rax)		\n\t"\
		"vmovaps		%%zmm7, 0x040(%%rax)	\n\t		vmovaps		%%zmm13,0x240(%%rax)		\n\t"\
		"vsubpd		%%zmm5,%%zmm0,%%zmm0	\n\t		vsubpd		%%zmm15,%%zmm10,%%zmm10\n\t"\
		"vsubpd		%%zmm4,%%zmm1,%%zmm1	\n\t		vsubpd		%%zmm14,%%zmm11,%%zmm11\n\t"\
	"vfmadd132pd	(%%rdi),%%zmm0,%%zmm5	\n\t	vfmadd132pd		(%%rdi),%%zmm10,%%zmm15\n\t"\
	"vfmadd132pd	(%%rdi),%%zmm1,%%zmm4	\n\t	vfmadd132pd		(%%rdi),%%zmm11,%%zmm14\n\t"\
		"vmovaps		%%zmm0,      (%%rbx)	\n\t		vmovaps		%%zmm10,0x200(%%rbx)		\n\t"\
		"vmovaps		%%zmm1, 0x040(%%rdx)	\n\t		vmovaps		%%zmm11,0x240(%%rdx)		\n\t"\
		"vmovaps		%%zmm5,      (%%rdx)	\n\t		vmovaps		%%zmm15,0x200(%%rdx)		\n\t"\
		"vmovaps		%%zmm4, 0x040(%%rbx)	\n\t		vmovaps		%%zmm14,0x240(%%rbx)		\n\t"\
		/*...Block 4: t06,t16,t26,t36	*/			/*...Block 8: t0E,t1E,t2E,t3E	*/\
		"addq		$0x100,%%rax			\n\t"\
		"addq		$0x100,%%rbx			\n\t"\
		"addq		$0x100,%%rcx			\n\t"\
		"addq		$0x100,%%rdx			\n\t"\
	/* Instead of loading regs 0, 8-10 and copying to 4, 12-14, use the later MULs to effect copying: */\
		"vmovaps		 (%%rcx),%%zmm4			\n\t		vmovaps		0x200(%%rcx),%%zmm12		\n\t"\
		"vmovaps		 (%%rdx),%%zmm6			\n\t		vmovaps		0x200(%%rdx),%%zmm14		\n\t"\
		"vmovaps	0x040(%%rcx),%%zmm1			\n\t		vmovaps		0x240(%%rcx),%%zmm13		\n\t"\
		"vmovaps	0x040(%%rdx),%%zmm3			\n\t		vmovaps		0x240(%%rdx),%%zmm11		\n\t"\
	/* Strategic register-usage for sincos consts here cuts loads from 16 to 4: */\
		"vmovaps		 (%%rsi),%%zmm2 		\n\t		vmovaps		0x040(%%rsi),%%zmm7 		\n\t"\
		"vmovaps	0x080(%%rsi),%%zmm5 		\n\t		vmovaps		0x0c0(%%rsi),%%zmm15		\n\t"\
		"vmulpd		%%zmm15,%%zmm1,%%zmm1		\n\t		vmulpd			%%zmm2 ,%%zmm13,%%zmm9	\n\t"\
		"vmulpd		%%zmm2 ,%%zmm3,%%zmm3		\n\t		vmulpd			%%zmm5 ,%%zmm11,%%zmm11	\n\t"\
		"vmulpd		%%zmm15,%%zmm4,%%zmm0		\n\t		vmulpd			%%zmm2 ,%%zmm12,%%zmm8	\n\t"\
		"vmulpd		 (%%rdx),%%zmm2,%%zmm2		\n\t		vmulpd			%%zmm5 ,%%zmm14,%%zmm10	\n\t"\
	"vfmsub132pd	%%zmm5 ,%%zmm1,%%zmm4		\n\t	vfmsub132pd			%%zmm7 ,%%zmm9 ,%%zmm12	\n\t"\
	"vfmadd132pd	%%zmm7 ,%%zmm3,%%zmm6		\n\t	vfmsub132pd			%%zmm15,%%zmm11,%%zmm14	\n\t"\
	"vfmadd132pd	0x040(%%rcx),%%zmm0,%%zmm5	\n\t	vfmadd132pd			%%zmm7 ,%%zmm8 ,%%zmm13	\n\t"\
	"vfmsub132pd	0x040(%%rdx),%%zmm2,%%zmm7	\n\t	vfmadd132pd		0x240(%%rdx),%%zmm10,%%zmm15\n\t"\
		"vsubpd		%%zmm6,%%zmm4,%%zmm4	\n\t		vsubpd		%%zmm14,%%zmm12,%%zmm12\n\t"\
		"vsubpd		%%zmm7,%%zmm5,%%zmm5	\n\t		vsubpd		%%zmm15,%%zmm13,%%zmm13\n\t"\
		"vmovaps		 (%%rbx),%%zmm2			\n\t		vmovaps		0x200(%%rbx),%%zmm10		\n\t"\
		"vmovaps	0x040(%%rbx),%%zmm0			\n\t		vmovaps		0x240(%%rbx),%%zmm8			\n\t"\
	"vfmadd132pd	(%%rdi),%%zmm4,%%zmm6	\n\t	vfmadd132pd		(%%rdi),%%zmm12,%%zmm14\n\t"\
	"vfmadd132pd	(%%rdi),%%zmm5,%%zmm7	\n\t	vfmadd132pd		(%%rdi),%%zmm13,%%zmm15\n\t"\
		"vmovaps		 (%%rbx),%%zmm2			\n\t		vmovaps		0x200(%%rbx),%%zmm10		\n\t"\
		"vmovaps	0x040(%%rbx),%%zmm0			\n\t		vmovaps		0x240(%%rbx),%%zmm8			\n\t"\
		"vmovaps		 (%%rbx),%%zmm1			\n\t		vmovaps		0x200(%%rbx),%%zmm9			\n\t"\
		"vmovaps	0x040(%%rbx),%%zmm3			\n\t		vmovaps		0x240(%%rbx),%%zmm11		\n\t"\
		"vmovaps		 (%%rbx),%%zmm1			\n\t		vmovaps		0x200(%%rbx),%%zmm9			\n\t"\
		"vmovaps	-0x40(%%rsi),%%zmm3			\n\t		vmovaps		-0x80(%%rsi),%%zmm11		\n\t"\
		"vmulpd		%%zmm11,%%zmm0,%%zmm0		\n\t		vmulpd			%%zmm3 ,%%zmm8 ,%%zmm8	\n\t"\
		"vmulpd		%%zmm11,%%zmm1,%%zmm1		\n\t		vmulpd			%%zmm3 ,%%zmm9 ,%%zmm9	\n\t"\
	"vfmsub132pd	%%zmm3 ,%%zmm0,%%zmm2		\n\t	vfmadd132pd			%%zmm11,%%zmm8 ,%%zmm10	\n\t"\
	"vfmadd132pd	0x40(%%rbx),%%zmm1,%%zmm3	\n\t	vfmsub132pd		0x240(%%rbx),%%zmm9,%%zmm11	\n\t"\
		"vmovaps		 (%%rax),%%zmm0			\n\t		vmovaps		0x200(%%rax),%%zmm8			\n\t"\
		"vmovaps	0x040(%%rax),%%zmm1			\n\t		vmovaps		0x240(%%rax),%%zmm9			\n\t"\
		"vsubpd		%%zmm2,%%zmm0,%%zmm0	\n\t		vsubpd		%%zmm10,%%zmm8,%%zmm8	\n\t"\
		"vsubpd		%%zmm3,%%zmm1,%%zmm1	\n\t		vsubpd		%%zmm11,%%zmm9,%%zmm9	\n\t"\
		"vaddpd			 (%%rax),%%zmm2,%%zmm2	\n\t		vaddpd		0x200(%%rax),%%zmm10,%%zmm10\n\t"\
		"vaddpd		0x040(%%rax),%%zmm3,%%zmm3	\n\t		vaddpd		0x240(%%rax),%%zmm11,%%zmm11\n\t"\
		"vsubpd		%%zmm4,%%zmm2,%%zmm2	\n\t		vsubpd		%%zmm12,%%zmm8,%%zmm8	\n\t"\
		"vsubpd		%%zmm5,%%zmm3,%%zmm3	\n\t		vsubpd		%%zmm13,%%zmm9,%%zmm9	\n\t"\
	"vfmadd132pd	(%%rdi),%%zmm2,%%zmm4	\n\t	vfmadd132pd		(%%rdi),%%zmm8,%%zmm12\n\t"\
	"vfmadd132pd	(%%rdi),%%zmm3,%%zmm5	\n\t	vfmadd132pd		(%%rdi),%%zmm9,%%zmm13\n\t"\
		"vmovaps		%%zmm2,      (%%rcx)	\n\t		vmovaps		%%zmm8 ,0x200(%%rcx)		\n\t"\
		"vmovaps		%%zmm3, 0x040(%%rcx)	\n\t		vmovaps		%%zmm9 ,0x240(%%rcx)		\n\t"\
		"vmovaps		%%zmm4,      (%%rax)	\n\t		vmovaps		%%zmm12,0x200(%%rax)		\n\t"\
		"vmovaps		%%zmm5, 0x040(%%rax)	\n\t		vmovaps		%%zmm13,0x240(%%rax)		\n\t"\
		"vsubpd		%%zmm7,%%zmm0,%%zmm0	\n\t		vsubpd		%%zmm15,%%zmm10,%%zmm10\n\t"\
		"vsubpd		%%zmm6,%%zmm1,%%zmm1	\n\t		vsubpd		%%zmm14,%%zmm11,%%zmm11\n\t"\
	"vfmadd132pd	(%%rdi),%%zmm0,%%zmm7	\n\t	vfmadd132pd		(%%rdi),%%zmm10,%%zmm15\n\t"\
	"vfmadd132pd	(%%rdi),%%zmm1,%%zmm6	\n\t	vfmadd132pd		(%%rdi),%%zmm11,%%zmm14\n\t"\
		"vmovaps		%%zmm0,      (%%rbx)	\n\t		vmovaps		%%zmm10,0x200(%%rbx)		\n\t"\
		"vmovaps		%%zmm1, 0x040(%%rdx)	\n\t		vmovaps		%%zmm11,0x240(%%rdx)		\n\t"\
		"vmovaps		%%zmm7,      (%%rdx)	\n\t		vmovaps		%%zmm15,0x200(%%rdx)		\n\t"\
		"vmovaps		%%zmm6, 0x040(%%rbx)	\n\t		vmovaps		%%zmm14,0x240(%%rbx)		\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__add1] "m" (Xadd1)\
		 ,[__add2] "m" (Xadd2)\
		 ,[__add3] "m" (Xadd3)\
		 ,[__add4] "m" (Xadd4)\
		 ,[__add5] "m" (Xadd5)\
		 ,[__add6] "m" (Xadd6)\
		 ,[__add7] "m" (Xadd7)\
		 ,[__r00] "m" (Xr00)\
		 ,[__r10] "m" (Xr10)\
		 ,[__r20] "m" (Xr20)\
		 ,[__r30] "m" (Xr30)\
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__cc0] "m" (Xcc0)\
		 ,[__c00] "m" (Xc00)\
		 ,[__c01] "m" (Xc01)\
		 ,[__c02] "m" (Xc02)\
		 ,[__c03] "m" (Xc03)\
		 ,[__c05] "m" (Xc05)\
		 ,[__c07] "m" (Xc07)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
	);\
	}

  #else	// USE_16_REG = False, i.e. use default 32-SIMD-reg version of macro:

   #ifdef USE_IMCI512	// 1st-gen Xeon Phi - Use modified 8x8 doubles-transpose algo [1a] from util.c:test_simd_transpose_8x8()

	#define SSE2_RADIX32_WRAPPER_DIF(Xadd0,Xadd1,Xadd2,Xadd3,Xadd4,Xadd5,Xadd6,Xadd7,Xr00,Xr10,Xr20,Xr30,Xisrt2,Xcc0,Xc00,Xc01,Xc02,Xc03,Xc05,Xc07)\
	{\
	__asm__ volatile (\
	/* Cf. IMCI512 8x8 doubles-transpose algo [1a] in util.c ...
	do mask-register-setting first, so as to not clobber any GPRs used by compute section: */\
		"movl $0b10101010,%%ebx	\n\t movl $0b11001100,%%ecx	\n\t movl $0b11110000,%%edx	\n\t"\
		"kmov	%%ebx,%%k1		\n\t kmov	%%ecx,%%k3		\n\t kmov	%%edx,%%k5		\n\t"\
		"knot	%%k1 ,%%k2		\n\t knot	%%k3 ,%%k4		\n\t"\
	/**** Start with 8-way interleaving: ****/\
		"movq	%[__add0],%%rax	\n\t"\
		"movq	%[__add1],%%rbx	\n\t"\
		"movq	%[__add2],%%rcx	\n\t"\
		"movq	%[__add3],%%rdx	\n\t"\
		"movq	%[__add4],%%r10	\n\t"\
		"movq	%[__add5],%%r11	\n\t"\
		"movq	%[__add6],%%r12	\n\t"\
		"movq	%[__add7],%%r13	\n\t"\
		"movq	%[__r00] ,%%rsi	\n\t"\
	/* See above 16-register version for transpose using columnwise gather-loads. */\
	/* a[j+p0]: Inputs from add0+[0,0x200,0x400,0x600,0x800,0xa00,0xc00,0xe00]. Outputs into local store at r0+[same byte offsets]: */\
		/* Read contiguous (untransposed) row-data from the data array: */\
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
		/* Write transposed data to the local store: */\
		"vmovaps		%%zmm5,0x000(%%rsi)		\n\t	vmovaps	%%zmm15,0x040(%%rsi)	\n\t"/* Out 0 */\
		"vmovaps		%%zmm6,0x800(%%rsi)		\n\t	vmovaps	%%zmm16,0x840(%%rsi)	\n\t"/* Out 4 */\
		"vmovaps		%%zmm8,0x400(%%rsi)		\n\t	vmovaps	%%zmm9 ,0x440(%%rsi)	\n\t"/* Out 2 */\
		"vmovaps		%%zmm0,0xc00(%%rsi)		\n\t	vmovaps	%%zmm10,0xc40(%%rsi)	\n\t"/* Out 6 */\
		"vmovaps		%%zmm1,0x200(%%rsi)		\n\t	vmovaps	%%zmm11,0x240(%%rsi)	\n\t"/* Out 1 */\
		"vmovaps		%%zmm2,0xa00(%%rsi)		\n\t	vmovaps	%%zmm12,0xa40(%%rsi)	\n\t"/* Out 5 */\
		"vmovaps		%%zmm4,0x600(%%rsi)		\n\t	vmovaps	%%zmm14,0x640(%%rsi)	\n\t"/* Out 3 */\
		"vmovaps		%%zmm7,0xe00(%%rsi)		\n\t	vmovaps	%%zmm17,0xe40(%%rsi)	\n\t"/* Out 7 */\
		"\n\t"\
	/* a[j+p2]: Inputs from add0+[0,0x100,0x200,0x300,0x400,0x500,0x600,0x700]+0x80. Outputs into r1+[same byte offsets]: */\
	"addq	$0x80 ,%%rax	\n\t"\
	"addq	$0x80 ,%%rbx	\n\t"\
	"addq	$0x80 ,%%rcx	\n\t"\
	"addq	$0x80 ,%%rdx	\n\t"\
	"addq	$0x80 ,%%r10	\n\t"\
	"addq	$0x80 ,%%r11	\n\t"\
	"addq	$0x80 ,%%r12	\n\t"\
	"addq	$0x80 ,%%r13	\n\t"\
	"addq	$0x100,%%rsi	\n\t"\
		/* Read contiguous (untransposed) row-data from the data array: */\
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
		/* Write transposed data to the local store: */\
		"vmovaps		%%zmm5,0x000(%%rsi)		\n\t	vmovaps	%%zmm15,0x040(%%rsi)	\n\t"/* Out 0 */\
		"vmovaps		%%zmm6,0x800(%%rsi)		\n\t	vmovaps	%%zmm16,0x840(%%rsi)	\n\t"/* Out 4 */\
		"vmovaps		%%zmm8,0x400(%%rsi)		\n\t	vmovaps	%%zmm9 ,0x440(%%rsi)	\n\t"/* Out 2 */\
		"vmovaps		%%zmm0,0xc00(%%rsi)		\n\t	vmovaps	%%zmm10,0xc40(%%rsi)	\n\t"/* Out 6 */\
		"vmovaps		%%zmm1,0x200(%%rsi)		\n\t	vmovaps	%%zmm11,0x240(%%rsi)	\n\t"/* Out 1 */\
		"vmovaps		%%zmm2,0xa00(%%rsi)		\n\t	vmovaps	%%zmm12,0xa40(%%rsi)	\n\t"/* Out 5 */\
		"vmovaps		%%zmm4,0x600(%%rsi)		\n\t	vmovaps	%%zmm14,0x640(%%rsi)	\n\t"/* Out 3 */\
		"vmovaps		%%zmm7,0xe00(%%rsi)		\n\t	vmovaps	%%zmm17,0xe40(%%rsi)	\n\t"/* Out 7 */\
		"\n\t"\
	/* a[j+p4]: Inputs from add0+[0,0x100,0x200,0x300,0x400,0x500,0x600,0x700]+0x80. Outputs into r2+[same byte offsets]: */\
	"addq	$0x80 ,%%rax	\n\t"\
	"addq	$0x80 ,%%rbx	\n\t"\
	"addq	$0x80 ,%%rcx	\n\t"\
	"addq	$0x80 ,%%rdx	\n\t"\
	"addq	$0x80 ,%%r10	\n\t"\
	"addq	$0x80 ,%%r11	\n\t"\
	"addq	$0x80 ,%%r12	\n\t"\
	"addq	$0x80 ,%%r13	\n\t"\
	"subq	$0x80,%%rsi	\n\t"\
		/* Read contiguous (untransposed) row-data from the data array: */\
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
		/* Write transposed data to the local store: */\
		"vmovaps		%%zmm5,0x000(%%rsi)		\n\t	vmovaps	%%zmm15,0x040(%%rsi)	\n\t"/* Out 0 */\
		"vmovaps		%%zmm6,0x800(%%rsi)		\n\t	vmovaps	%%zmm16,0x840(%%rsi)	\n\t"/* Out 4 */\
		"vmovaps		%%zmm8,0x400(%%rsi)		\n\t	vmovaps	%%zmm9 ,0x440(%%rsi)	\n\t"/* Out 2 */\
		"vmovaps		%%zmm0,0xc00(%%rsi)		\n\t	vmovaps	%%zmm10,0xc40(%%rsi)	\n\t"/* Out 6 */\
		"vmovaps		%%zmm1,0x200(%%rsi)		\n\t	vmovaps	%%zmm11,0x240(%%rsi)	\n\t"/* Out 1 */\
		"vmovaps		%%zmm2,0xa00(%%rsi)		\n\t	vmovaps	%%zmm12,0xa40(%%rsi)	\n\t"/* Out 5 */\
		"vmovaps		%%zmm4,0x600(%%rsi)		\n\t	vmovaps	%%zmm14,0x640(%%rsi)	\n\t"/* Out 3 */\
		"vmovaps		%%zmm7,0xe00(%%rsi)		\n\t	vmovaps	%%zmm17,0xe40(%%rsi)	\n\t"/* Out 7 */\
		"\n\t"\
	/* a[j+p6]: Inputs from add0+[0,0x100,0x200,0x300,0x400,0x500,0x600,0x700]+0x180. Outputs into r3+[same byte offsets]: */\
	"addq	$0x80 ,%%rax	\n\t"\
	"addq	$0x80 ,%%rbx	\n\t"\
	"addq	$0x80 ,%%rcx	\n\t"\
	"addq	$0x80 ,%%rdx	\n\t"\
	"addq	$0x80 ,%%r10	\n\t"\
	"addq	$0x80 ,%%r11	\n\t"\
	"addq	$0x80 ,%%r12	\n\t"\
	"addq	$0x80 ,%%r13	\n\t"\
	"addq	$0x100,%%rsi	\n\t"\
		/* Read contiguous (untransposed) row-data from the data array: */\
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
		/* Write transposed data to the local store: */\
		"vmovaps		%%zmm5,0x000(%%rsi)		\n\t	vmovaps	%%zmm15,0x040(%%rsi)	\n\t"/* Out 0 */\
		"vmovaps		%%zmm6,0x800(%%rsi)		\n\t	vmovaps	%%zmm16,0x840(%%rsi)	\n\t"/* Out 4 */\
		"vmovaps		%%zmm8,0x400(%%rsi)		\n\t	vmovaps	%%zmm9 ,0x440(%%rsi)	\n\t"/* Out 2 */\
		"vmovaps		%%zmm0,0xc00(%%rsi)		\n\t	vmovaps	%%zmm10,0xc40(%%rsi)	\n\t"/* Out 6 */\
		"vmovaps		%%zmm1,0x200(%%rsi)		\n\t	vmovaps	%%zmm11,0x240(%%rsi)	\n\t"/* Out 1 */\
		"vmovaps		%%zmm2,0xa00(%%rsi)		\n\t	vmovaps	%%zmm12,0xa40(%%rsi)	\n\t"/* Out 5 */\
		"vmovaps		%%zmm4,0x600(%%rsi)		\n\t	vmovaps	%%zmm14,0x640(%%rsi)	\n\t"/* Out 3 */\
		"vmovaps		%%zmm7,0xe00(%%rsi)		\n\t	vmovaps	%%zmm17,0xe40(%%rsi)	\n\t"/* Out 7 */\
		"\n\t"\
		/************************************************************************/							/************************************************************************************************************/\
		/* Forward DIF radix-32 pass on the interleaved block1 and block2 data: */							/* Next 2 blocks operate on odd-indexed elements from the unpck*pd commands which we stored to temporaries: */\
		/************************************************************************/							/************************************************************************************************************/\
	/*...Block 0: */																						/*...Block 3: */\
		"movq	%[__isrt2],%%rsi				\n\t	leaq	0x11c0(%%rsi),%%rdi	/* two */	\n\t"\
	/* SSE2_RADIX4_DIF_4WRAPPER(c00,c08,c10,c18,r00)...DIF_4WRAPPER_2NDOFTWO(c04,c0C,c14,c1C,r08)	SSE2_RADIX4_DIF_4TWIDDLE(r20,r24,r22,r26,r20,c01)...DIF_4TWIDDLE_2NDOFTWO(r28,r2C,r2A,r2E,r28,c05) */\
		"movq		%[__r00],%%rcx				\n\t	/*addq		$0x200,%%rcx // __r08 */	\n\t	movq %[__r20],%%r10 \n\t movq %[__c01],%%r11 \n\t leaq 0x080(%%r10),%%r12 \n\t	movq %%r10,%%r13 \n\t"\
		"movq		%[__c00],%%rdx				\n\t	/*addq		$0x200,%%rdx // __c04 */	\n\t	vmovaps		 (%%r10),%%zmm18 			\n\t	vmovaps		0x200(%%r10),%%zmm26		\n\t"\
		"vmovaps			 (%%rcx),%%zmm0		\n\t	vmovaps		0x200(%%rcx),%%zmm8			\n\t	vmovaps		 (%%r12),%%zmm20			\n\t	vmovaps		0x200(%%r12),%%zmm28		\n\t"\
		"vmovaps		0x040(%%rcx),%%zmm1		\n\t	vmovaps		0x240(%%rcx),%%zmm9			\n\t	vmovaps	0x040(%%r10),%%zmm19			\n\t	vmovaps		0x240(%%r10),%%zmm27		\n\t"\
		"vmovaps		%%zmm0,%%zmm2			\n\t	vmovaps		%%zmm8,%%zmm10				\n\t	vmovaps	0x040(%%r12),%%zmm21			\n\t	vmovaps		0x240(%%r12),%%zmm29		\n\t"\
		"vmovaps		%%zmm1,%%zmm3			\n\t	vmovaps		%%zmm9,%%zmm11				\n\t	vmovaps		 (%%r11),%%zmm22			\n\t	vmovaps		0x200(%%r11),%%zmm30		\n\t"\
		"vmulpd			 (%%rdx),%%zmm0,%%zmm0	\n\t	vmulpd		0x200(%%rdx),%%zmm8 ,%%zmm8	\n\t	vmovaps	0x040(%%r11),%%zmm23			\n\t	vmovaps		0x240(%%r11),%%zmm31		\n\t"\
		"vmulpd			 (%%rdx),%%zmm1,%%zmm1	\n\t	vmulpd		0x200(%%rdx),%%zmm9 ,%%zmm9	\n\t"\
/********** Where do the zmm6,7 come from here? ***********/\
	"vfnmadd231pd	0x040(%%rdx),%%zmm7,%%zmm0	\n\t vfnmadd231pd	0x240(%%rdx),%%zmm11,%%zmm8 \n\t"\
	" vfmadd231pd	0x040(%%rdx),%%zmm6,%%zmm1	\n\t vfmadd231pd	0x240(%%rdx),%%zmm10,%%zmm9 \n\t	vmulpd		%%zmm22,%%zmm18,%%zmm16		\n\t	vmulpd		%%zmm30,%%zmm26,%%zmm24	\n\t"\
		"vmovaps		%%zmm0,%%zmm2			\n\t	vmovaps		%%zmm8,%%zmm10				\n\t	vmulpd		%%zmm22,%%zmm19,%%zmm17		\n\t	vmulpd		%%zmm30,%%zmm27,%%zmm25	\n\t"\
		"vmovaps		%%zmm1,%%zmm3			\n\t	vmovaps		%%zmm9,%%zmm11				\n\t vfnmadd231pd	%%zmm23,%%zmm19,%%zmm16		\n\t vfnmadd231pd	%%zmm31,%%zmm27,%%zmm24		\n\t"\
		"addq		$0x080,%%rdx																\n\t vfmadd231pd	%%zmm23,%%zmm18,%%zmm17		\n\t vfmadd231pd	%%zmm31,%%zmm26,%%zmm25		\n\t"\
		"vmovaps		0x080(%%rcx),%%zmm6		\n\t	vmovaps		0x280(%%rcx),%%zmm14		\n\t	vmovaps		%%zmm20,%%zmm22				\n\t	vmovaps		%%zmm28,%%zmm30		\n\t"\
		"vmovaps		0x0c0(%%rcx),%%zmm7		\n\t	vmovaps		0x2c0(%%rcx),%%zmm15		\n\t	vmovaps		%%zmm21,%%zmm23				\n\t	vmovaps		%%zmm29,%%zmm31		\n\t"\
		"																								vmulpd	  0x080(%%r11),%%zmm20,%%zmm20	\n\t	vmulpd		0x280(%%r11),%%zmm28,%%zmm28\n\t"\
		"																								vmulpd	  0x080(%%r11),%%zmm21,%%zmm21	\n\t	vmulpd		0x280(%%r11),%%zmm29,%%zmm29\n\t"\
		"vmulpd			 (%%rdx),%%zmm6,%%zmm4	\n\t	vmulpd		0x200(%%rdx),%%zmm14,%%zmm12\n\t vfnmadd231pd 0x0c0(%%r11),%%zmm23,%%zmm20	\n\t vfnmadd231pd	0x2c0(%%r11),%%zmm31,%%zmm28\n\t"\
		"vmulpd			 (%%rdx),%%zmm7,%%zmm5	\n\t	vmulpd		0x200(%%rdx),%%zmm15,%%zmm13\n\t vfmadd231pd  0x0c0(%%r11),%%zmm22,%%zmm21	\n\t vfmadd231pd	0x2c0(%%r11),%%zmm30,%%zmm29\n\t"\
	"vfnmadd231pd	0x040(%%rdx),%%zmm7,%%zmm4	\n\t vfnmadd231pd	0x240(%%rdx),%%zmm15,%%zmm12\n\t	vmovaps		%%zmm16,%%zmm18				\n\t	vmovaps		%%zmm24,%%zmm26		\n\t"\
	" vfmadd231pd	0x040(%%rdx),%%zmm6,%%zmm5	\n\t vfmadd231pd	0x240(%%rdx),%%zmm14,%%zmm13\n\t	vmovaps		%%zmm17,%%zmm19				\n\t	vmovaps		%%zmm25,%%zmm27		\n\t"\
		"vaddpd		%%zmm4,%%zmm0,%%zmm0		\n\t	vaddpd		%%zmm12,%%zmm8 ,%%zmm8		\n\t	addq		$0x100,%%r12				\n\t	addq		$0x180,%%r11			\n\t"\
		"vaddpd		%%zmm5,%%zmm1,%%zmm1		\n\t	vaddpd		%%zmm13,%%zmm9 ,%%zmm9		\n\t	vmovaps			 (%%r12),%%zmm22		\n\t	vmovaps		0x200(%%r12),%%zmm30		\n\t"\
		"vsubpd		%%zmm4,%%zmm2,%%zmm2		\n\t	vsubpd		%%zmm12,%%zmm10,%%zmm10		\n\t	vmovaps		0x040(%%r12),%%zmm23		\n\t	vmovaps		0x240(%%r12),%%zmm31		\n\t"\
		"vsubpd		%%zmm5,%%zmm3,%%zmm3		\n\t	vsubpd		%%zmm13,%%zmm11,%%zmm11		\n\t	vaddpd		%%zmm20,%%zmm16,%%zmm16		\n\t	vaddpd		%%zmm28,%%zmm24,%%zmm24		\n\t"\
		"addq		$0x100,%%rdx																\n\t	vaddpd		%%zmm21,%%zmm17,%%zmm17		\n\t	vaddpd		%%zmm29,%%zmm25,%%zmm25		\n\t"\
		"vmovaps		0x180(%%rcx),%%zmm6		\n\t	vmovaps		0x380(%%rcx),%%zmm14		\n\t	vsubpd		%%zmm20,%%zmm18,%%zmm18		\n\t	vsubpd		%%zmm28,%%zmm26,%%zmm26		\n\t"\
		"vmovaps		0x1c0(%%rcx),%%zmm7		\n\t	vmovaps		0x3c0(%%rcx),%%zmm15		\n\t	vsubpd		%%zmm21,%%zmm19,%%zmm19		\n\t	vsubpd		%%zmm29,%%zmm27,%%zmm27		\n\t"\
		"vmulpd			 (%%rdx),%%zmm6,%%zmm4	\n\t	vmulpd		0x200(%%rdx),%%zmm14,%%zmm12\n\t	vmulpd		   (%%r11),%%zmm22,%%zmm20	\n\t	vmulpd		0x200(%%r11),%%zmm30,%%zmm28\n\t"\
		"vmulpd			 (%%rdx),%%zmm7,%%zmm5	\n\t	vmulpd		0x200(%%rdx),%%zmm15,%%zmm13\n\t	vmulpd		   (%%r11),%%zmm23,%%zmm21	\n\t	vmulpd		0x200(%%r11),%%zmm31,%%zmm29\n\t"\
	"vfnmadd231pd	0x040(%%rdx),%%zmm7,%%zmm4	\n\t vfnmadd231pd	0x240(%%rdx),%%zmm15,%%zmm12\n\t vfnmadd231pd 0x040(%%r11),%%zmm23,%%zmm20	\n\t vfnmadd231pd	0x240(%%r11),%%zmm31,%%zmm28\n\t"\
	" vfmadd231pd	0x040(%%rdx),%%zmm6,%%zmm5	\n\t vfmadd231pd	0x240(%%rdx),%%zmm14,%%zmm13\n\t vfmadd231pd  0x040(%%r11),%%zmm22,%%zmm21	\n\t vfmadd231pd	0x240(%%r11),%%zmm30,%%zmm29\n\t"\
		"vmovaps		%%zmm5,0x040(%%rcx)		\n\t	vmovaps		%%zmm13,0x240(%%rcx)		\n\t	vmovaps		%%zmm21,0x040(%%r13)		\n\t	vmovaps		%%zmm29,0x240(%%r13)		\n\t"\
		"vmovaps		%%zmm4,     (%%rcx)		\n\t	vmovaps		%%zmm12,0x200(%%rcx)		\n\t	vmovaps		%%zmm20,     (%%r13)		\n\t	vmovaps		%%zmm28,0x200(%%r13)		\n\t"\
		"subq		$0x080,%%rdx																\n\t	addq	$0x100,%%r10					\n\t	subq	$0x080,%%r11					\n\t"\
		"vmovaps		0x100(%%rcx),%%zmm6		\n\t	vmovaps		0x300(%%rcx),%%zmm14		\n\t	vmovaps			 (%%r10),%%zmm22		\n\t	vmovaps		0x200(%%r10),%%zmm30		\n\t"\
		"vmovaps		0x140(%%rcx),%%zmm7		\n\t	vmovaps		0x340(%%rcx),%%zmm15		\n\t	vmovaps		0x040(%%r10),%%zmm23		\n\t	vmovaps		0x240(%%r10),%%zmm31		\n\t"\
		"vmulpd			 (%%rdx),%%zmm6,%%zmm4	\n\t	vmulpd		0x200(%%rdx),%%zmm14,%%zmm12\n\t	vmulpd		   (%%r11),%%zmm22,%%zmm20	\n\t	vmulpd		0x200(%%r11),%%zmm30,%%zmm28\n\t"\
		"vmulpd			 (%%rdx),%%zmm7,%%zmm5	\n\t	vmulpd		0x200(%%rdx),%%zmm15,%%zmm13\n\t	vmulpd		   (%%r11),%%zmm23,%%zmm21	\n\t	vmulpd		0x200(%%r11),%%zmm31,%%zmm29\n\t"\
	"vfnmadd231pd	0x040(%%rdx),%%zmm7,%%zmm4	\n\t vfnmadd231pd	0x240(%%rdx),%%zmm15,%%zmm12\n\t vfnmadd231pd 0x040(%%r11),%%zmm23,%%zmm20	\n\t vfnmadd231pd	0x240(%%r11),%%zmm31,%%zmm28\n\t"\
	" vfmadd231pd	0x040(%%rdx),%%zmm6,%%zmm5	\n\t vfmadd231pd	0x240(%%rdx),%%zmm14,%%zmm13\n\t vfmadd231pd  0x040(%%r11),%%zmm22,%%zmm21	\n\t vfmadd231pd	0x240(%%r11),%%zmm30,%%zmm29\n\t"\
		"vmovaps		%%zmm5,%%zmm7			\n\t	vmovaps		%%zmm13,%%zmm15				\n\t	vmovaps		%%zmm21,%%zmm23				\n\t	vmovaps		%%zmm29,%%zmm31		\n\t"\
		"vmovaps		%%zmm4,%%zmm6			\n\t	vmovaps		%%zmm12,%%zmm14				\n\t	vmovaps		%%zmm20,%%zmm22				\n\t	vmovaps		%%zmm28,%%zmm30		\n\t"\
		"vsubpd			 (%%rcx),%%zmm4,%%zmm4	\n\t	vsubpd		0x200(%%rcx),%%zmm12,%%zmm12\n\t	vsubpd			 (%%r13),%%zmm20,%%zmm20\n\t	vsubpd		0x200(%%r13),%%zmm28,%%zmm28\n\t"\
		"vsubpd		0x040(%%rcx),%%zmm5,%%zmm5	\n\t	vsubpd		0x240(%%rcx),%%zmm13,%%zmm13\n\t	vsubpd		0x040(%%r13),%%zmm21,%%zmm21\n\t	vsubpd		0x240(%%r13),%%zmm29,%%zmm29\n\t"\
		"vaddpd			 (%%rcx),%%zmm6,%%zmm6	\n\t	vaddpd		0x200(%%rcx),%%zmm14,%%zmm14\n\t	vaddpd			 (%%r13),%%zmm22,%%zmm22\n\t	vaddpd		0x200(%%r13),%%zmm30,%%zmm30\n\t"\
		"vaddpd		0x040(%%rcx),%%zmm7,%%zmm7	\n\t	vaddpd		0x240(%%rcx),%%zmm15,%%zmm15\n\t	vaddpd		0x040(%%r13),%%zmm23,%%zmm23\n\t	vaddpd		0x240(%%r13),%%zmm31,%%zmm31\n\t"\
		"vsubpd		%%zmm6,%%zmm0,%%zmm0		\n\t	vsubpd		%%zmm14,%%zmm8,%%zmm8		\n\t	vsubpd		%%zmm22,%%zmm16,%%zmm16		\n\t	vsubpd		%%zmm30,%%zmm24,%%zmm24\n\t"\
		"vsubpd		%%zmm7,%%zmm1,%%zmm1		\n\t	vsubpd		%%zmm15,%%zmm9,%%zmm9		\n\t	vsubpd		%%zmm23,%%zmm17,%%zmm17		\n\t	vsubpd		%%zmm31,%%zmm25,%%zmm25\n\t"\
		"vsubpd		%%zmm5,%%zmm2,%%zmm2		\n\t	vsubpd		%%zmm13,%%zmm10,%%zmm10		\n\t	vsubpd		%%zmm21,%%zmm18,%%zmm18		\n\t	vsubpd		%%zmm29,%%zmm26,%%zmm26\n\t"\
		"vsubpd		%%zmm4,%%zmm3,%%zmm3		\n\t	vsubpd		%%zmm12,%%zmm11,%%zmm11		\n\t	vsubpd		%%zmm20,%%zmm19,%%zmm19		\n\t	vsubpd		%%zmm28,%%zmm27,%%zmm27\n\t"\
	"vmovaps %%zmm28,(%%r13) \n\t"/* spill zmm28 to make room for two */"vmovaps (%%rdi),%%zmm28\n\t"/* two */\
	"vfmadd132pd	%%zmm28,%%zmm0,%%zmm6		\n\t	vfmadd132pd	%%zmm28,%%zmm8 ,%%zmm14		\n\t	vfmadd132pd	%%zmm28,%%zmm16,%%zmm22		\n\t	vfmadd132pd	%%zmm28,%%zmm24,%%zmm30			\n\t"\
	"vfmadd132pd	%%zmm28,%%zmm1,%%zmm7		\n\t	vfmadd132pd	%%zmm28,%%zmm9 ,%%zmm15		\n\t	vfmadd132pd	%%zmm28,%%zmm17,%%zmm23		\n\t	vfmadd132pd	%%zmm28,%%zmm25,%%zmm31			\n\t"\
	"vfmadd132pd	%%zmm28,%%zmm2,%%zmm5		\n\t	vfmadd132pd	%%zmm28,%%zmm10,%%zmm13		\n\t	vfmadd132pd	%%zmm28,%%zmm18,%%zmm21		\n\t	vfmadd132pd	%%zmm28,%%zmm26,%%zmm29			\n\t"\
	"vfmadd132pd	%%zmm28,%%zmm3,%%zmm4		\n\t	vfmadd132pd	%%zmm28,%%zmm11,%%zmm12		\n\t	vfmadd132pd	%%zmm28,%%zmm19,%%zmm20		\n\t	vfmadd132pd	(%%r13),%%zmm27,%%zmm28			\n\t"\
		"												vmovaps		%%zmm14,0x200(%%rcx)		\n\t														vmovaps		%%zmm30,0x200(%%r13)			\n\t"\
		"												vmovaps		%%zmm15,0x240(%%rcx)		\n\t														vmovaps		%%zmm31,0x240(%%r13)			\n\t"\
		"												vmovaps		%%zmm10,%%zmm14				\n\t														vmovaps		%%zmm26,%%zmm30					\n\t"\
		"												vmovaps		%%zmm13,%%zmm15				\n\t														vmovaps		%%zmm29,%%zmm31					\n\t"\
		"												vsubpd		%%zmm12,%%zmm10,%%zmm10		\n\t														vsubpd		%%zmm28,%%zmm26,%%zmm26	\n\t"\
		"												vsubpd		%%zmm11,%%zmm13,%%zmm13		\n\t														vsubpd		%%zmm27,%%zmm29,%%zmm29	\n\t"\
		"												vaddpd		%%zmm12,%%zmm14,%%zmm14		\n\t														vaddpd		%%zmm28,%%zmm30,%%zmm30	\n\t"\
		"												vaddpd		%%zmm11,%%zmm15,%%zmm15		\n\t														vaddpd		%%zmm27,%%zmm31,%%zmm31	\n\t"\
		/*****	SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r00) *****/												/*****	SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r20) *****/\
		"												vmovaps		0x200(%%rcx),%%zmm11			\n\t														vmovaps		0x200(%%r13),%%zmm27			\n\t"\
		"												vmovaps		0x240(%%rcx),%%zmm12			\n\t														vmovaps		0x240(%%r13),%%zmm28			\n\t"\
	"vmovaps %%zmm13,(%%rcx) \n\t"/* spill zmm13 to make room for isrt2 */"vmovaps (%%rsi),%%zmm13 	\n\t	vmovaps %%zmm29,(%%r13) \n\t"/* spill zmm29 to make room for isrt2 */"vmovaps (%%rsi),%%zmm29	\n\t"/*isrt2*/\
		"vsubpd		%%zmm11,%%zmm6,%%zmm6		\n\t	vfnmadd231pd	%%zmm13,%%zmm10,%%zmm2		\n\t		vsubpd		%%zmm27,%%zmm22,%%zmm22	\n\t	vfnmadd231pd	%%zmm29,%%zmm26,%%zmm18		\n\t"\
		"vsubpd		%%zmm9 ,%%zmm0,%%zmm0		\n\t	vfnmadd231pd	%%zmm13,%%zmm15,%%zmm5		\n\t		vsubpd		%%zmm25,%%zmm16,%%zmm16	\n\t	vfnmadd231pd	%%zmm29,%%zmm31,%%zmm21		\n\t"\
		"vsubpd		%%zmm12,%%zmm7,%%zmm7		\n\t	vfnmadd231pd	%%zmm13,%%zmm14,%%zmm4		\n\t		vsubpd		%%zmm28,%%zmm23,%%zmm23	\n\t	vfnmadd231pd	%%zmm29,%%zmm30,%%zmm20		\n\t"\
		"vsubpd		%%zmm8 ,%%zmm1,%%zmm1		\n\t	vfnmadd231pd	(%%rcx),%%zmm13,%%zmm3		\n\t		vsubpd		%%zmm24,%%zmm17,%%zmm17	\n\t	vfnmadd231pd	(%%r13),%%zmm29,%%zmm19		\n\t"\
		"vmovaps		%%zmm6,0x200(%%rcx)		\n\t		vmovaps		%%zmm2,0x280(%%rcx)	\n\t		vmovaps		%%zmm22,0x200(%%r13)	\n\t		vmovaps		%%zmm18,0x280(%%r13)	\n\t"\
		"vmovaps		%%zmm0,0x100(%%rcx)		\n\t		vmovaps		%%zmm5,0x180(%%rcx)	\n\t		vmovaps		%%zmm16,0x100(%%r13)	\n\t		vmovaps		%%zmm21,0x180(%%r13)	\n\t"\
		"vmovaps		%%zmm7,0x240(%%rcx)		\n\t		vmovaps		%%zmm4,0x2c0(%%rcx)	\n\t		vmovaps		%%zmm23,0x240(%%r13)	\n\t		vmovaps		%%zmm20,0x2c0(%%r13)	\n\t"\
		"vmovaps		%%zmm1,0x340(%%rcx)		\n\t		vmovaps		%%zmm3,0x3c0(%%rcx)	\n\t		vmovaps		%%zmm17,0x340(%%r13)	\n\t		vmovaps		%%zmm19,0x3c0(%%r13)	\n\t"\
	"vmovaps	0x80(%%rdi),%%zmm13																	\n\t	vmovaps	0x80(%%rdi),%%zmm29 \n\t"/* sqrt2 *//* rcol: *2 => *sqrt2 due to above FMA leaving zmm26,13-15 unmultiplied-by-isrt2: */\
	"vfmadd132pd	(%%rdi),%%zmm6,%%zmm11		\n\t	 vfmadd132pd	%%zmm13,%%zmm2,%%zmm10	\n\t	vfmadd132pd	(%%rdi),%%zmm22,%%zmm27	\n\t	 vfmadd132pd	%%zmm29,%%zmm18,%%zmm26	\n\t"\
	"vfmadd132pd	(%%rdi),%%zmm0,%%zmm9		\n\t	 vfmadd132pd	%%zmm13,%%zmm5,%%zmm15	\n\t	vfmadd132pd	(%%rdi),%%zmm16,%%zmm25	\n\t	 vfmadd132pd	%%zmm29,%%zmm21,%%zmm31	\n\t"\
	"vfmadd132pd	(%%rdi),%%zmm7,%%zmm12		\n\t	 vfmadd132pd	%%zmm13,%%zmm4,%%zmm14	\n\t	vfmadd132pd	(%%rdi),%%zmm23,%%zmm28	\n\t	 vfmadd132pd	%%zmm29,%%zmm20,%%zmm30	\n\t"\
	"vfmadd132pd	(%%rdi),%%zmm1,%%zmm8		\n\t	 vfmadd132pd	(%%rcx),%%zmm3,%%zmm13	\n\t	vfmadd132pd	(%%rdi),%%zmm17,%%zmm24	\n\t	 vfmadd132pd	(%%r13),%%zmm19,%%zmm29	\n\t"\
		"vmovaps		%%zmm11,     (%%rcx)	\n\t		vmovaps		%%zmm10,0x080(%%rcx)	\n\t		vmovaps		%%zmm27,     (%%r13)	\n\t		vmovaps		%%zmm26,0x080(%%r13)	\n\t"\
		"vmovaps		%%zmm9 ,0x300(%%rcx)	\n\t		vmovaps		%%zmm15,0x380(%%rcx)	\n\t		vmovaps		%%zmm25,0x300(%%r13)	\n\t		vmovaps		%%zmm31,0x380(%%r13)	\n\t"\
		"vmovaps		%%zmm12,0x040(%%rcx)	\n\t		vmovaps		%%zmm14,0x0c0(%%rcx)	\n\t		vmovaps		%%zmm28,0x040(%%r13)	\n\t		vmovaps		%%zmm30,0x0c0(%%r13)	\n\t"\
		"vmovaps		%%zmm8 ,0x140(%%rcx)	\n\t		vmovaps		%%zmm13,0x1c0(%%rcx)	\n\t		vmovaps		%%zmm24,0x140(%%r13)	\n\t		vmovaps		%%zmm29,0x1c0(%%r13)	\n\t"\
	/*...Block 2: */	/*...Block 4: */\
	/* SSE2_RADIX4_DIF_4WRAPPER(c02,c0A,c12,c1A,r10)	...DIF_4WRAPPER_2NDOFTWO(c06,c0E,c16,c1E,r18) */	/*	SSE2_RADIX4_DIF_4TWIDDLE(r30,r34,r32,r36,r30,c03)	...DIF_4TWIDDLE_2NDOFTWO(r38,r3C,r3A,r3E,r38,c07) */\
		"movq		%[__r10],%%rcx				\n\t		/*addq		$0x200,%%rcx // __r18 */	\n\t		movq %[__r30],%%r10 \n\t movq %[__c03],%%r11 \n\t movq %%r10,%%r12 \n\t addq $0x080,%%r12	\n\t"\
		"movq		%[__c02],%%rdx				\n\t		/*addq		$0x200,%%rdx // __c06 */	\n\t		vmovaps		 (%%r10),%%zmm16 \n\t movq %%r10,%%r13 \n\t vmovaps 0x200(%%r10),%%zmm24	\n\t"\
		"vmovaps			 (%%rcx),%%zmm0		\n\t		vmovaps		0x200(%%rcx),%%zmm8			\n\t		vmovaps		 (%%r12),%%zmm20		\n\t		vmovaps		0x200(%%r12),%%zmm28		\n\t"\
		"vmovaps		0x040(%%rcx),%%zmm1		\n\t		vmovaps		0x240(%%rcx),%%zmm9			\n\t		vmovaps	0x040(%%r10),%%zmm17		\n\t		vmovaps		0x240(%%r10),%%zmm25		\n\t"\
		"vmovaps		%%zmm0,%%zmm2			\n\t		vmovaps		%%zmm8,%%zmm10				\n\t		vmovaps	0x040(%%r12),%%zmm21		\n\t		vmovaps		0x240(%%r12),%%zmm29		\n\t"\
		"vmovaps		%%zmm1,%%zmm3			\n\t		vmovaps		%%zmm9,%%zmm11				\n\t		vmovaps		 (%%r11),%%zmm22		\n\t		vmovaps		0x200(%%r11),%%zmm30		\n\t"\
		"vmulpd			 (%%rdx),%%zmm0,%%zmm0	\n\t		vmulpd		0x200(%%rdx),%%zmm8,%%zmm8	\n\t		vmovaps	0x040(%%r11),%%zmm23		\n\t		vmovaps		0x240(%%r11),%%zmm31		\n\t"\
		"vmulpd			 (%%rdx),%%zmm1,%%zmm1	\n\t		vmulpd		0x200(%%rdx),%%zmm9,%%zmm9	\n\t		vmovaps	%%zmm16,%%zmm18		\n\t		vmovaps		%%zmm24,%%zmm26		\n\t"\
	"vfnmadd231pd	0x040(%%rdx),%%zmm3,%%zmm0	\n\t	vfnmadd231pd	0x240(%%rdx),%%zmm11,%%zmm8 \n\t		vmovaps	%%zmm17,%%zmm19		\n\t		vmovaps		%%zmm25,%%zmm27		\n\t"\
	" vfmadd231pd	0x040(%%rdx),%%zmm2,%%zmm1	\n\t	 vfmadd231pd	0x240(%%rdx),%%zmm10,%%zmm9 \n\t		vmulpd		%%zmm22,%%zmm16,%%zmm16\n\t		vmulpd		%%zmm30,%%zmm24,%%zmm24\n\t"\
		"vmovaps		%%zmm0,%%zmm2			\n\t		vmovaps		%%zmm8,%%zmm10				\n\t		vmulpd		%%zmm22,%%zmm17,%%zmm17\n\t		vmulpd		%%zmm30,%%zmm25,%%zmm25\n\t"\
		"vmovaps		%%zmm1,%%zmm3			\n\t		vmovaps		%%zmm9,%%zmm11				\n\t	vfnmadd231pd	%%zmm23,%%zmm19,%%zmm16		\n\t	vfnmadd231pd	%%zmm31,%%zmm27,%%zmm24		\n\t"\
		"addq		$0x080,%%rdx																\n\t	 vfmadd231pd	%%zmm23,%%zmm18,%%zmm17		\n\t	 vfmadd231pd	%%zmm31,%%zmm26,%%zmm25		\n\t"\
		"vmovaps		0x080(%%rcx),%%zmm4		\n\t		vmovaps		0x280(%%rcx),%%zmm12		\n\t		vmovaps		%%zmm20,%%zmm22		\n\t		vmovaps		%%zmm28,%%zmm30		\n\t"\
		"vmovaps		0x0c0(%%rcx),%%zmm5		\n\t		vmovaps		0x2c0(%%rcx),%%zmm13		\n\t		vmovaps		%%zmm21,%%zmm23		\n\t		vmovaps		%%zmm29,%%zmm31		\n\t"\
		"vmovaps		%%zmm4,%%zmm6			\n\t		vmovaps		%%zmm12,%%zmm14				\n\t		vmulpd		0x080(%%r11),%%zmm20,%%zmm20\n\t		vmulpd		0x280(%%r11),%%zmm28,%%zmm28\n\t"\
		"vmovaps		%%zmm5,%%zmm7			\n\t		vmovaps		%%zmm13,%%zmm15				\n\t		vmulpd		0x080(%%r11),%%zmm21,%%zmm21\n\t		vmulpd		0x280(%%r11),%%zmm29,%%zmm29\n\t"\
		"vmulpd			 (%%rdx),%%zmm4,%%zmm4	\n\t		vmulpd		0x200(%%rdx),%%zmm12,%%zmm12\n\t		vmovaps		%%zmm16,%%zmm18		\n\t		vmovaps		%%zmm24,%%zmm26		\n\t"\
		"vmulpd			 (%%rdx),%%zmm5,%%zmm5	\n\t		vmulpd		0x200(%%rdx),%%zmm13,%%zmm13\n\t		vmovaps		%%zmm17,%%zmm19		\n\t		vmovaps		%%zmm25,%%zmm27		\n\t"\
	"vfnmadd231pd	0x040(%%rdx),%%zmm7,%%zmm4	\n\t	vfnmadd231pd	0x240(%%rdx),%%zmm15,%%zmm12\n\t	vfnmadd231pd	0x0c0(%%r11),%%zmm23,%%zmm20\n\t	vfnmadd231pd	0x2c0(%%r11),%%zmm31,%%zmm28\n\t"\
	" vfmadd231pd	0x040(%%rdx),%%zmm6,%%zmm5	\n\t	 vfmadd231pd	0x240(%%rdx),%%zmm14,%%zmm13\n\t	 vfmadd231pd	0x0c0(%%r11),%%zmm22,%%zmm21\n\t	 vfmadd231pd	0x2c0(%%r11),%%zmm30,%%zmm29\n\t"\
		"vaddpd		%%zmm4,%%zmm0,%%zmm0		\n\t		vaddpd		%%zmm12,%%zmm8,%%zmm8	\n\t		addq		$0x100,%%r12			\n\t		addq		$0x180,%%r11			\n\t"\
		"vaddpd		%%zmm5,%%zmm1,%%zmm1		\n\t		vaddpd		%%zmm13,%%zmm9,%%zmm9	\n\t		vmovaps			 (%%r12),%%zmm22	\n\t		vaddpd		%%zmm28,%%zmm24,%%zmm24\n\t"\
		"vsubpd		%%zmm4,%%zmm2,%%zmm2		\n\t		vsubpd		%%zmm12,%%zmm10,%%zmm10\n\t		vmovaps		0x040(%%r12),%%zmm23	\n\t		vaddpd		%%zmm29,%%zmm25,%%zmm25\n\t"\
		"vsubpd		%%zmm5,%%zmm3,%%zmm3		\n\t		vsubpd		%%zmm13,%%zmm11,%%zmm11\n\t		vaddpd		%%zmm20,%%zmm16,%%zmm16\n\t		vsubpd		%%zmm28,%%zmm26,%%zmm26\n\t"\
		"addq		$0x100,%%rdx																\n\t		vaddpd		%%zmm21,%%zmm17,%%zmm17\n\t		vsubpd		%%zmm29,%%zmm27,%%zmm27\n\t"\
		"vmovaps		0x180(%%rcx),%%zmm4		\n\t		vmovaps		0x380(%%rcx),%%zmm12		\n\t		vsubpd		%%zmm20,%%zmm18,%%zmm18\n\t		vmovaps		0x200(%%r12),%%zmm28		\n\t"\
		"vmovaps		0x1c0(%%rcx),%%zmm5		\n\t		vmovaps		0x3c0(%%rcx),%%zmm13		\n\t		vsubpd		%%zmm21,%%zmm19,%%zmm19\n\t		vmovaps		0x240(%%r12),%%zmm29		\n\t"\
		"vmovaps		%%zmm4,%%zmm6			\n\t		vmovaps		%%zmm12,%%zmm14				\n\t		vmovaps		%%zmm22,%%zmm20		\n\t		vmovaps		0x200(%%r12),%%zmm30		\n\t"\
		"vmovaps		%%zmm5,%%zmm7			\n\t		vmovaps		%%zmm13,%%zmm15				\n\t		vmovaps		%%zmm23,%%zmm21		\n\t		vmovaps		0x240(%%r12),%%zmm31		\n\t"\
		"vmulpd			 (%%rdx),%%zmm4,%%zmm4	\n\t		vmulpd		0x200(%%rdx),%%zmm12,%%zmm12\n\t		vmulpd			 (%%r11),%%zmm20,%%zmm20\n\t		vmulpd		0x200(%%r11),%%zmm28,%%zmm28\n\t"\
		"vmulpd			 (%%rdx),%%zmm5,%%zmm5	\n\t		vmulpd		0x200(%%rdx),%%zmm13,%%zmm13\n\t		vmulpd			 (%%r11),%%zmm21,%%zmm21\n\t		vmulpd		0x200(%%r11),%%zmm29,%%zmm29\n\t"\
	"vfnmadd231pd	0x040(%%rdx),%%zmm7,%%zmm4	\n\t	vfnmadd231pd	0x240(%%rdx),%%zmm15,%%zmm12\n\t	vfnmadd231pd	0x040(%%r11),%%zmm23,%%zmm20\n\t	vfnmadd231pd	0x240(%%r11),%%zmm31,%%zmm28\n\t"\
	" vfmadd231pd	0x040(%%rdx),%%zmm6,%%zmm5	\n\t	 vfmadd231pd	0x240(%%rdx),%%zmm14,%%zmm13\n\t	 vfmadd231pd	0x040(%%r11),%%zmm22,%%zmm21\n\t	 vfmadd231pd	0x240(%%r11),%%zmm30,%%zmm29\n\t"\
		"vmovaps		%%zmm5,0x040(%%rcx)		\n\t		vmovaps		%%zmm13,0x240(%%rcx)		\n\t		vmovaps		%%zmm21,0x040(%%r13)	\n\t		vmovaps		%%zmm29,0x240(%%r13)		\n\t"\
		"vmovaps		%%zmm4,     (%%rcx)		\n\t		vmovaps		%%zmm12,0x200(%%rcx)		\n\t		vmovaps		%%zmm20,     (%%r13)	\n\t		vmovaps		%%zmm28,0x200(%%r13)		\n\t"\
		"subq		$0x080,%%rdx																\n\t		addq	$0x100,%%r10					\n\t		subq	$0x080,%%r11					\n\t"\
		"vmovaps		0x100(%%rcx),%%zmm4		\n\t		vmovaps		0x300(%%rcx),%%zmm12		\n\t		vmovaps			 (%%r10),%%zmm20	\n\t		vmovaps		0x200(%%r10),%%zmm28		\n\t"\
		"vmovaps		0x140(%%rcx),%%zmm5		\n\t		vmovaps		0x340(%%rcx),%%zmm13		\n\t		vmovaps		0x040(%%r10),%%zmm21	\n\t		vmovaps		0x240(%%r10),%%zmm29		\n\t"\
		"vmovaps		%%zmm4,%%zmm6			\n\t		vmovaps		%%zmm12,%%zmm14				\n\t		vmovaps		%%zmm20,%%zmm22		\n\t		vmovaps		0x200(%%r10),%%zmm30		\n\t"\
		"vmovaps		%%zmm5,%%zmm7			\n\t		vmovaps		%%zmm13,%%zmm15				\n\t		vmovaps		%%zmm21,%%zmm23		\n\t		vmovaps		0x240(%%r10),%%zmm31		\n\t"\
		"vmulpd			 (%%rdx),%%zmm4,%%zmm4	\n\t		vmulpd		0x200(%%rdx),%%zmm12,%%zmm12\n\t		vmulpd			 (%%r11),%%zmm20,%%zmm20\n\t		vmulpd		0x200(%%r11),%%zmm28,%%zmm28\n\t"\
		"vmulpd			 (%%rdx),%%zmm5,%%zmm5	\n\t		vmulpd		0x200(%%rdx),%%zmm13,%%zmm13\n\t		vmulpd			 (%%r11),%%zmm21,%%zmm21\n\t		vmulpd		0x200(%%r11),%%zmm29,%%zmm29\n\t"\
	"vfnmadd231pd	0x040(%%rdx),%%zmm7,%%zmm4	\n\t	vfnmadd231pd	0x240(%%rdx),%%zmm15,%%zmm12\n\t	vfnmadd231pd	0x040(%%r11),%%zmm23,%%zmm20\n\t	vfnmadd231pd	0x240(%%r11),%%zmm31,%%zmm28\n\t"\
	" vfmadd231pd	0x040(%%rdx),%%zmm6,%%zmm5	\n\t	 vfmadd231pd	0x240(%%rdx),%%zmm14,%%zmm13\n\t	 vfmadd231pd	0x040(%%r11),%%zmm22,%%zmm21\n\t	 vfmadd231pd	0x240(%%r11),%%zmm30,%%zmm29\n\t"\
		"vmovaps		%%zmm5,%%zmm7			\n\t		vmovaps		%%zmm13,%%zmm15				\n\t		vmovaps		%%zmm21,%%zmm23		\n\t		vmovaps		%%zmm29,%%zmm31		\n\t"\
		"vmovaps		%%zmm4,%%zmm6			\n\t		vmovaps		%%zmm12,%%zmm14				\n\t		vmovaps		%%zmm20,%%zmm22		\n\t		vmovaps		%%zmm28,%%zmm30		\n\t"\
		"vsubpd			 (%%rcx),%%zmm4,%%zmm4	\n\t		vsubpd		0x200(%%rcx),%%zmm12,%%zmm12\n\t		vsubpd			 (%%r13),%%zmm20,%%zmm20\n\t		vsubpd		0x200(%%r13),%%zmm28,%%zmm28\n\t"\
		"vsubpd		0x040(%%rcx),%%zmm5,%%zmm5	\n\t		vsubpd		0x240(%%rcx),%%zmm13,%%zmm13\n\t		vsubpd		0x040(%%r13),%%zmm21,%%zmm21\n\t		vsubpd		0x240(%%r13),%%zmm29,%%zmm29\n\t"\
		"vaddpd			 (%%rcx),%%zmm6,%%zmm6	\n\t		vaddpd		0x200(%%rcx),%%zmm14,%%zmm14\n\t		vaddpd			 (%%r13),%%zmm22,%%zmm22\n\t		vaddpd		0x200(%%r13),%%zmm30,%%zmm30\n\t"\
		"vaddpd		0x040(%%rcx),%%zmm7,%%zmm7	\n\t		vaddpd		0x240(%%rcx),%%zmm15,%%zmm15\n\t		vaddpd		0x040(%%r13),%%zmm23,%%zmm23\n\t		vaddpd		0x240(%%r13),%%zmm31,%%zmm31\n\t"\
		"vsubpd		%%zmm6,%%zmm0,%%zmm0		\n\t		vsubpd		%%zmm14,%%zmm8,%%zmm8	\n\t		vsubpd		%%zmm22,%%zmm16,%%zmm16\n\t		vsubpd		%%zmm30,%%zmm24,%%zmm24\n\t"\
		"vsubpd		%%zmm7,%%zmm1,%%zmm1		\n\t		vsubpd		%%zmm15,%%zmm9,%%zmm9	\n\t		vsubpd		%%zmm23,%%zmm17,%%zmm17\n\t		vsubpd		%%zmm31,%%zmm25,%%zmm25\n\t"\
		"vsubpd		%%zmm5,%%zmm2,%%zmm2		\n\t		vsubpd		%%zmm13,%%zmm10,%%zmm10\n\t		vsubpd		%%zmm21,%%zmm18,%%zmm18\n\t		vsubpd		%%zmm29,%%zmm26,%%zmm26\n\t"\
		"vsubpd		%%zmm4,%%zmm3,%%zmm3		\n\t		vsubpd		%%zmm12,%%zmm11,%%zmm11\n\t		vsubpd		%%zmm20,%%zmm19,%%zmm19\n\t		vsubpd		%%zmm28,%%zmm27,%%zmm27\n\t"\
	"vmovaps	%%zmm12,(%%rcx) \n\t"/* spill zmm12 to make room for two */"vmovaps (%%rdi),%%zmm12	\n\t	vmovaps	%%zmm28,(%%r13) \n\t"/* spill zmm28 to make room for two */"vmovaps (%%rdi),%%zmm28	\n\t"/* two */\
	"vfmadd132pd	%%zmm12,%%zmm0,%%zmm6		\n\t	vfmadd132pd	%%zmm12,%%zmm8 ,%%zmm14			\n\t	vfmadd132pd	%%zmm28,%%zmm16,%%zmm22		\n\t	vfmadd132pd	%%zmm28,%%zmm24,%%zmm30			\n\t"\
	"vfmadd132pd	%%zmm12,%%zmm1,%%zmm7		\n\t	vfmadd132pd	%%zmm12,%%zmm9 ,%%zmm15			\n\t	vfmadd132pd	%%zmm28,%%zmm17,%%zmm23		\n\t	vfmadd132pd	%%zmm28,%%zmm25,%%zmm31			\n\t"\
	"vfmadd132pd	%%zmm12,%%zmm2,%%zmm5		\n\t	vfmadd132pd	%%zmm12,%%zmm10,%%zmm13			\n\t	vfmadd132pd	%%zmm28,%%zmm18,%%zmm21		\n\t	vfmadd132pd	%%zmm28,%%zmm26,%%zmm29			\n\t"\
	"vfmadd132pd	%%zmm12,%%zmm3,%%zmm4		\n\t	vfmadd132pd	(%%rcx),%%zmm11,%%zmm12			\n\t	vfmadd132pd	%%zmm28,%%zmm19,%%zmm20		\n\t	vfmadd132pd	(%%r13),%%zmm27,%%zmm28			\n\t"\
		"												vmovaps		%%zmm14,0x200(%%rcx)			\n\t														vmovaps		%%zmm30,0x200(%%r13)			\n\t"\
		"												vmovaps		%%zmm15,0x240(%%rcx)			\n\t														vmovaps		%%zmm31,0x240(%%r13)			\n\t"\
		"												vmovaps		%%zmm10,%%zmm14					\n\t														vmovaps		%%zmm26,%%zmm30					\n\t"\
		"												vmovaps		%%zmm13,%%zmm15					\n\t														vmovaps		%%zmm29,%%zmm31					\n\t"\
		"												vsubpd		%%zmm12,%%zmm10,%%zmm10	\n\t														vsubpd		%%zmm28,%%zmm26,%%zmm26	\n\t"\
		"												vsubpd		%%zmm11,%%zmm13,%%zmm13	\n\t														vsubpd		%%zmm27,%%zmm29,%%zmm29	\n\t"\
		"												vaddpd		%%zmm12,%%zmm14,%%zmm14	\n\t														vaddpd		%%zmm28,%%zmm30,%%zmm30	\n\t"\
		"												vaddpd		%%zmm11,%%zmm15,%%zmm15	\n\t														vaddpd		%%zmm27,%%zmm31,%%zmm31	\n\t"\
		/*****	SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r10) *****/												/*****	SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r30) *****/\
		"												vmovaps		0x200(%%rcx),%%zmm11			\n\t														vmovaps		0x200(%%r13),%%zmm27			\n\t"\
		"												vmovaps		0x240(%%rcx),%%zmm12			\n\t														vmovaps		0x240(%%r13),%%zmm28			\n\t"\
	"vmovaps %%zmm13,(%%rcx) \n\t"/* spill zmm13 to make room for isrt2 */"vmovaps (%%rsi),%%zmm13	\n\t	vmovaps %%zmm29,(%%r13) \n\t"/* spill zmm29 to make room for isrt2 */"vmovaps (%%rsi),%%zmm29	\n\t"/*isrt2*/\
		"vsubpd		%%zmm11,%%zmm6,%%zmm6		\n\t	vfnmadd231pd	%%zmm13,%%zmm10,%%zmm2		\n\t		vsubpd		%%zmm27,%%zmm22,%%zmm22	\n\t	vfnmadd231pd	%%zmm29,%%zmm26,%%zmm18		\n\t"\
		"vsubpd		%%zmm9,%%zmm0,%%zmm0		\n\t	vfnmadd231pd	%%zmm13,%%zmm15,%%zmm5		\n\t		vsubpd		%%zmm25,%%zmm16,%%zmm16	\n\t	vfnmadd231pd	%%zmm29,%%zmm31,%%zmm21		\n\t"\
		"vsubpd		%%zmm12,%%zmm7,%%zmm7		\n\t	vfnmadd231pd	%%zmm13,%%zmm14,%%zmm4		\n\t		vsubpd		%%zmm28,%%zmm23,%%zmm23	\n\t	vfnmadd231pd	%%zmm29,%%zmm30,%%zmm20		\n\t"\
		"vsubpd		%%zmm8 ,%%zmm1,%%zmm1		\n\t	vfnmadd231pd	(%%rcx),%%zmm13,%%zmm3		\n\t		vsubpd		%%zmm24,%%zmm17,%%zmm17	\n\t	vfnmadd231pd	(%%r13),%%zmm29,%%zmm19		\n\t"\
		"vmovaps		%%zmm6,0x200(%%rcx)		\n\t		vmovaps		%%zmm2,0x280(%%rcx)	\n\t		vmovaps		%%zmm22,0x200(%%r13)	\n\t		vmovaps		%%zmm18,0x280(%%r13)	\n\t"\
		"vmovaps		%%zmm0,0x100(%%rcx)		\n\t		vmovaps		%%zmm5,0x180(%%rcx)	\n\t		vmovaps		%%zmm16,0x100(%%r13)	\n\t		vmovaps		%%zmm21,0x180(%%r13)	\n\t"\
		"vmovaps		%%zmm7,0x240(%%rcx)		\n\t		vmovaps		%%zmm4,0x2c0(%%rcx)	\n\t		vmovaps		%%zmm23,0x240(%%r13)	\n\t		vmovaps		%%zmm20,0x2c0(%%r13)	\n\t"\
		"vmovaps		%%zmm1,0x340(%%rcx)		\n\t		vmovaps		%%zmm3,0x3c0(%%rcx)	\n\t		vmovaps		%%zmm17,0x340(%%r13)	\n\t		vmovaps		%%zmm19,0x3c0(%%r13)	\n\t"\
	"vmovaps	0x80(%%rdi),%%zmm13																	\n\t	vmovaps	0x80(%%rdi),%%zmm29 \n\t"/* sqrt2 *//* rcol: *2 => *sqrt2 due to above FMA leaving zmm26,13-15 unmultiplied-by-isrt2: */\
	"vfmadd132pd	(%%rdi),%%zmm6,%%zmm11		\n\t	 vfmadd132pd	%%zmm13,%%zmm2,%%zmm10	\n\t	vfmadd132pd	(%%rdi),%%zmm22,%%zmm27	\n\t	 vfmadd132pd	%%zmm29,%%zmm18,%%zmm26	\n\t"\
	"vfmadd132pd	(%%rdi),%%zmm0,%%zmm9		\n\t	 vfmadd132pd	%%zmm13,%%zmm5,%%zmm15	\n\t	vfmadd132pd	(%%rdi),%%zmm16,%%zmm25	\n\t	 vfmadd132pd	%%zmm29,%%zmm21,%%zmm31	\n\t"\
	"vfmadd132pd	(%%rdi),%%zmm7,%%zmm12		\n\t	 vfmadd132pd	%%zmm13,%%zmm4,%%zmm14	\n\t	vfmadd132pd	(%%rdi),%%zmm23,%%zmm28	\n\t	 vfmadd132pd	%%zmm29,%%zmm20,%%zmm30	\n\t"\
	"vfmadd132pd	(%%rdi),%%zmm1,%%zmm8		\n\t	 vfmadd132pd	(%%rcx),%%zmm3,%%zmm13	\n\t	vfmadd132pd	(%%rdi),%%zmm17,%%zmm24	\n\t	 vfmadd132pd	(%%r13),%%zmm19,%%zmm29	\n\t"\
		"vmovaps		%%zmm11,     (%%rcx)	\n\t		vmovaps		%%zmm10,0x080(%%rcx)	\n\t		vmovaps		%%zmm27,     (%%r13)	\n\t		vmovaps		%%zmm26,0x080(%%r13)	\n\t"\
		"vmovaps		%%zmm9 ,0x300(%%rcx)	\n\t		vmovaps		%%zmm15,0x380(%%rcx)	\n\t		vmovaps		%%zmm25,0x300(%%r13)	\n\t		vmovaps		%%zmm31,0x380(%%r13)	\n\t"\
		"vmovaps		%%zmm12,0x040(%%rcx)	\n\t		vmovaps		%%zmm14,0x0c0(%%rcx)	\n\t		vmovaps		%%zmm28,0x040(%%r13)	\n\t		vmovaps		%%zmm30,0x0c0(%%r13)	\n\t"\
		"vmovaps		%%zmm8 ,0x140(%%rcx)	\n\t		vmovaps		%%zmm13,0x1c0(%%rcx)	\n\t		vmovaps		%%zmm24,0x140(%%r13)	\n\t		vmovaps		%%zmm29,0x1c0(%%r13)	\n\t"\
		"\n\t"\
	/**********************************************************************************/\
	/*...and now do eight radix-4 transforms, including the internal twiddle factors: */\
	/**********************************************************************************/\
		"movq		%[__isrt2],%%rsi			\n\t"\
	/*...Block 1: t00,t10,t20,t30	*/					/*...Block 5: t08,t18,t28,t38	*/					/*...Block 3: t04,t14,t24,t34	*/			/*...Block 7: t0C,t1C,t2C,t3C	*/\
		"movq		%[__r20],%%rcx				\n\t					leaq	0x100(%%rcx),%%r12	\n\t		vmovaps		 (%%r12),%%zmm20			\n\t		vmovaps		0x200(%%r12),%%zmm28		\n\t"\
		"movq		%[__r30],%%rdx				\n\t					leaq	0x100(%%rdx),%%r13	\n\t		vmovaps	0x040(%%r12),%%zmm21			\n\t		vmovaps		0x240(%%r12),%%zmm29		\n\t"\
		"movq		%[__r00],%%rax				\n\t					leaq	0x100(%%rax),%%r10	\n\t		vmovaps		 (%%r13),%%zmm22			\n\t		vmovaps		0x200(%%r13),%%zmm30		\n\t"\
		"movq		%[__r10],%%rbx				\n\t					leaq	0x100(%%rbx),%%r11	\n\t		vmovaps	0x040(%%r13),%%zmm19			\n\t		vmovaps		0x240(%%r13),%%zmm27		\n\t"\
		"vmovaps		 (%%rax),%%zmm0			\n\t		vmovaps		0x200(%%rax),%%zmm8			\n\t	vmovaps 0x080(%%rsi),%%zmm23\n\t vmovaps 0x040(%%rsi),%%zmm31 \n\t"/* Use 7,15 for s,c 2 cut loads */\
		"vmovaps	0x040(%%rax),%%zmm1			\n\t		vmovaps		0x240(%%rax),%%zmm9			\n\t		vmulpd		%%zmm23,%%zmm21,%%zmm17	\n\t		vmulpd		 %%zmm31,%%zmm29,%%zmm25	\n\t"\
		"vmovaps		 (%%rbx),%%zmm2			\n\t		vmovaps		0x200(%%rbx),%%zmm10		\n\t		vmulpd		%%zmm31,%%zmm19,%%zmm19	\n\t		vmulpd		 %%zmm23,%%zmm27,%%zmm27	\n\t"\
		"vmovaps	0x040(%%rbx),%%zmm3			\n\t		vmovaps		0x240(%%rbx),%%zmm11		\n\t		vmulpd		%%zmm23,%%zmm20,%%zmm16	\n\t		vmulpd		 %%zmm31,%%zmm28,%%zmm24	\n\t"\
		"vaddpd			  %%zmm0,%%zmm2,%%zmm2	\n\t		vaddpd		%%zmm9 ,%%zmm10,%%zmm10		\n\t		vmulpd		%%zmm31,%%zmm22,%%zmm18	\n\t		vmulpd		 %%zmm23,%%zmm30,%%zmm26	\n\t"\
		"vaddpd			  %%zmm1,%%zmm3,%%zmm3	\n\t		vaddpd		%%zmm8 ,%%zmm11,%%zmm11		\n\t	vfmsub132pd	%%zmm31,%%zmm17,%%zmm20		\n\t	vfmsub132pd		 %%zmm23,%%zmm25,%%zmm28	\n\t"\
		"vsubpd			 (%%rbx),%%zmm0,%%zmm0	\n\t		vsubpd		0x240(%%rbx),%%zmm8,%%zmm8	\n\t	vfmsub132pd	%%zmm23,%%zmm19,%%zmm22		\n\t	vfmsub132pd		 %%zmm31,%%zmm27,%%zmm30	\n\t"\
		"vsubpd		0x040(%%rbx),%%zmm1,%%zmm1	\n\t		vsubpd		0x200(%%rbx),%%zmm9,%%zmm9	\n\t	vfmadd132pd	%%zmm31,%%zmm16,%%zmm21		\n\t	vfmadd132pd		 %%zmm23,%%zmm24,%%zmm29	\n\t"\
		"vmovaps		 (%%rcx),%%zmm4			\n\t		vmovaps		0x200(%%rcx),%%zmm12		\n\t	vfmadd132pd 0x40(%%r13),%%zmm18,%%zmm23		\n\t	vfmadd132pd 0x240(%%r13),%%zmm26,%%zmm31	\n\t"\
		"vmovaps	0x040(%%rcx),%%zmm5			\n\t		vmovaps		0x240(%%rcx),%%zmm13		\n\t		vsubpd		%%zmm22,%%zmm20,%%zmm20	\n\t		vsubpd		%%zmm30,%%zmm28,%%zmm28\n\t"\
		"vmovaps		 (%%rdx),%%zmm6			\n\t		vmovaps		0x200(%%rdx),%%zmm14		\n\t		vsubpd		%%zmm23,%%zmm21,%%zmm21	\n\t		vsubpd		%%zmm31,%%zmm29,%%zmm29\n\t"\
		"vmovaps	0x040(%%rdx),%%zmm7			\n\t		vmovaps		0x240(%%rdx),%%zmm15		\n\t		vmovaps		 (%%r11),%%zmm18			\n\t		vmovaps		0x200(%%r11),%%zmm26		\n\t"\
		"vaddpd		%%zmm4,%%zmm6,%%zmm6		\n\t		vsubpd		%%zmm13,%%zmm12,%%zmm12		\n\t		vmovaps	0x040(%%r11),%%zmm19			\n\t		vmovaps		0x240(%%r11),%%zmm27		\n\t"\
		"vaddpd		%%zmm5,%%zmm7,%%zmm7		\n\t		vaddpd		0x200(%%rcx),%%zmm13,%%zmm13\n\t	vfmadd132pd	(%%rdi),%%zmm20,%%zmm22	\n\t	vfmadd132pd		(%%rdi),%%zmm28,%%zmm30\n\t"\
		"vsubpd			 (%%rdx),%%zmm4,%%zmm4	\n\t		vaddpd			 %%zmm15,%%zmm14,%%zmm14\n\t	vfmadd132pd	(%%rdi),%%zmm21,%%zmm23	\n\t	vfmadd132pd		(%%rdi),%%zmm29,%%zmm31\n\t"\
		"vsubpd		0x040(%%rdx),%%zmm5,%%zmm5	\n\t		vsubpd		0x200(%%rdx),%%zmm15,%%zmm15\n\t		vsubpd		0x040(%%r11),%%zmm18,%%zmm18\n\t		vaddpd		0x240(%%r11),%%zmm26,%%zmm26\n\t"\
		"vsubpd		%%zmm6,%%zmm2,%%zmm2		\n\t		vmulpd		(%%rsi),%%zmm12,%%zmm12\n\t		vaddpd			 (%%r11),%%zmm19,%%zmm19\n\t		vsubpd		0x200(%%r11),%%zmm27,%%zmm27\n\t"\
		"vsubpd		%%zmm7,%%zmm3,%%zmm3		\n\t		vmulpd		(%%rsi),%%zmm13,%%zmm13\n\t	vmovaps %%zmm22,(%%r12) \n\t"/* spill zmm22 to make room for isrt2 */"vmovaps (%%rsi),%%zmm22	\n\t"\
	"vfmadd132pd	(%%rdi),%%zmm2,%%zmm6		\n\t	vfnmadd231pd	(%%rsi),%%zmm14,%%zmm12\n\t		vmovaps		 (%%r10),%%zmm16			\n\t		vmovaps		 0x200(%%r10),%%zmm24		\n\t"\
	"vfmadd132pd	(%%rdi),%%zmm3,%%zmm7		\n\t	vfnmadd231pd	(%%rsi),%%zmm15,%%zmm13\n\t		vmovaps	0x040(%%r10),%%zmm17			\n\t		vmovaps		 0x240(%%r10),%%zmm25		\n\t"\
		"vmovaps		%%zmm2,      (%%rcx)	\n\t	 vfmadd132pd	0x80(%%rdi),%%zmm12,%%zmm14\n\t	vfnmadd231pd	%%zmm22,%%zmm18,%%zmm16		\n\t	vfnmadd231pd	%%zmm22,%%zmm26,%%zmm24		\n\t"\
		"vmovaps		%%zmm3, 0x040(%%rcx)	\n\t	 vfmadd132pd	0x80(%%rdi),%%zmm13,%%zmm15\n\t	vfnmadd231pd	%%zmm22,%%zmm19,%%zmm17		\n\t	vfnmadd231pd	%%zmm22,%%zmm27,%%zmm25		\n\t"\
		"										\n\t		vsubpd		%%zmm12,%%zmm8,%%zmm8	\n\t	 vfmadd213pd		 (%%r10),%%zmm22,%%zmm18\n\t	 vfmadd213pd	0x200(%%r10),%%zmm22,%%zmm26\n\t"\
		"										\n\t		vsubpd		%%zmm13,%%zmm10,%%zmm10\n\t	 vfmadd213pd	0x040(%%r10),%%zmm22,%%zmm19\n\t	 vfmadd213pd	0x240(%%r10),%%zmm22,%%zmm27\n\t"\
		"vmovaps		%%zmm6,      (%%rax)	\n\t	 vfmadd132pd	(%%rdi),%%zmm8,%%zmm12\n\t	vmovaps	(%%r12),%%zmm22	\n\t"/* restore spill */\
		"vmovaps		%%zmm7, 0x040(%%rax)	\n\t	 vfmadd132pd	(%%rdi),%%zmm10,%%zmm13\n\t		vsubpd		%%zmm22,%%zmm18,%%zmm18	\n\t		vsubpd		%%zmm28,%%zmm24,%%zmm24\n\t"\
		"vsubpd		%%zmm5,%%zmm0,%%zmm0		\n\t		vmovaps		%%zmm8 ,0x200(%%rcx)		\n\t		vsubpd		%%zmm23,%%zmm19,%%zmm19	\n\t		vsubpd		%%zmm29,%%zmm25,%%zmm25\n\t"\
		"vsubpd		%%zmm4,%%zmm1,%%zmm1		\n\t		vmovaps		%%zmm10,0x240(%%rcx)		\n\t	vfmadd132pd	(%%rdi),%%zmm18,%%zmm22	\n\t	vfmadd132pd		(%%rdi),%%zmm24,%%zmm28\n\t"\
	"vfmadd132pd	(%%rdi),%%zmm0,%%zmm5		\n\t		vmovaps		%%zmm12,0x200(%%rax)		\n\t	vfmadd132pd	(%%rdi),%%zmm19,%%zmm23	\n\t	vfmadd132pd		(%%rdi),%%zmm25,%%zmm29\n\t"\
	"vfmadd132pd	(%%rdi),%%zmm1,%%zmm4		\n\t		vmovaps		%%zmm13,0x240(%%rax)		\n\t		vmovaps		%%zmm18,      (%%r12)		\n\t		vmovaps		%%zmm24,0x200(%%r12)		\n\t"\
		"vmovaps		%%zmm0,      (%%rbx)	\n\t		vsubpd		%%zmm15,%%zmm11,%%zmm11\n\t		vmovaps		%%zmm19, 0x040(%%r12)		\n\t		vmovaps		%%zmm25,0x240(%%r12)		\n\t"\
		"vmovaps		%%zmm1, 0x040(%%rdx)	\n\t		vsubpd		%%zmm14,%%zmm9,%%zmm9	\n\t		vmovaps		%%zmm22,      (%%r10)		\n\t		vmovaps		%%zmm28,0x200(%%r10)		\n\t"\
		"										\n\t	 vfmadd132pd	(%%rdi),%%zmm11,%%zmm15\n\t		vmovaps		%%zmm23, 0x040(%%r10)		\n\t		vmovaps		%%zmm29,0x240(%%r10)		\n\t"\
		"										\n\t	 vfmadd132pd	(%%rdi),%%zmm9 ,%%zmm14\n\t		vsubpd		%%zmm21,%%zmm16,%%zmm16	\n\t		vsubpd		%%zmm31,%%zmm26,%%zmm26\n\t"\
		"vmovaps		%%zmm5,      (%%rdx)	\n\t		vmovaps		%%zmm11,0x200(%%rbx)		\n\t		vsubpd		%%zmm20,%%zmm17,%%zmm17	\n\t		vsubpd		%%zmm30,%%zmm27,%%zmm27\n\t"\
		"vmovaps		%%zmm4, 0x040(%%rbx)	\n\t		vmovaps		%%zmm9 ,0x240(%%rdx)		\n\t	vfmadd132pd	(%%rdi),%%zmm16,%%zmm21	\n\t	vfmadd132pd		(%%rdi),%%zmm26,%%zmm31\n\t"\
		"													vmovaps		%%zmm15,0x200(%%rdx)		\n\t	vfmadd132pd	(%%rdi),%%zmm17,%%zmm20	\n\t	vfmadd132pd		(%%rdi),%%zmm27,%%zmm30\n\t"\
		"													vmovaps		%%zmm14,0x240(%%rbx)		\n\t		vmovaps		%%zmm16,      (%%r11)		\n\t		vmovaps		%%zmm26,0x200(%%r11)		\n\t"\
		"																										vmovaps		%%zmm17, 0x040(%%r13)		\n\t		vmovaps		%%zmm27,0x240(%%r13)		\n\t"\
		"																										vmovaps		%%zmm21,      (%%r13)		\n\t		vmovaps		%%zmm31,0x200(%%r13)		\n\t"\
		"																										vmovaps		%%zmm20, 0x040(%%r11)		\n\t		vmovaps		%%zmm30,0x240(%%r11)		\n\t"\
		"addq		$0x0c0,%%rsi	\n\t"/* cc1 */\
	/*...Block 2: t02,t12,t22,t32	*/					/*...Block 6: t0A,t1A,t2A,t3A */					/*...Block 4: t06,t16,t26,t36	*/			/*...Block 8: t0E,t1E,t2E,t3E	*/\
		"addq		$0x080,%%rcx				\n\t					addq	$0x080,%%r12	\n\t		vmovaps		 (%%r12),%%zmm20			\n\t		vmovaps		0x200(%%r12),%%zmm28		\n\t"\
		"addq		$0x080,%%rdx				\n\t					addq	$0x080,%%r13	\n\t		vmovaps		 (%%r13),%%zmm22			\n\t		vmovaps		0x200(%%r13),%%zmm30		\n\t"\
		"addq		$0x080,%%rax				\n\t					addq	$0x080,%%r10	\n\t		vmovaps	0x040(%%r12),%%zmm17			\n\t		vmovaps		0x240(%%r12),%%zmm29		\n\t"\
		"addq		$0x080,%%rbx				\n\t					addq	$0x080,%%r11	\n\t		vmovaps	0x040(%%r13),%%zmm19			\n\t		vmovaps		0x240(%%r13),%%zmm27		\n\t"\
		"vmovaps		 (%%rcx),%%zmm4			\n\t		vmovaps		0x200(%%rcx),%%zmm12		\n\t		vmovaps		 (%%rsi),%%zmm18			\n\t		vmovaps		0x040(%%rsi),%%zmm23		\n\t"\
		"vmovaps		 (%%rdx),%%zmm6			\n\t		vmovaps		0x200(%%rdx),%%zmm14		\n\t		vmovaps	0x080(%%rsi),%%zmm21			\n\t		vmovaps		0x0c0(%%rsi),%%zmm31		\n\t"\
		"vmovaps	0x040(%%rcx),%%zmm5			\n\t		vmovaps		0x240(%%rcx),%%zmm9			\n\t		vmulpd		%%zmm31,%%zmm17,%%zmm17	\n\t		vmulpd			%%zmm18,%%zmm29,%%zmm25	\n\t"\
		"vmovaps	0x040(%%rdx),%%zmm3			\n\t		vmovaps		0x240(%%rdx),%%zmm11		\n\t		vmulpd		%%zmm18,%%zmm19,%%zmm19	\n\t		vmulpd			%%zmm21,%%zmm27,%%zmm27	\n\t"\
		"vmovaps		 (%%rsi),%%zmm15		\n\t		vmovaps		0x040(%%rsi),%%zmm10		\n\t		vmulpd		%%zmm31,%%zmm20,%%zmm16	\n\t		vmulpd			%%zmm18,%%zmm28,%%zmm24	\n\t"\
		"vmovaps	0x080(%%rsi),%%zmm7 		\n\t		vmovaps		0x0c0(%%rsi),%%zmm13		\n\t		vmulpd		 (%%r13),%%zmm18,%%zmm18	\n\t		vmulpd			%%zmm21,%%zmm30,%%zmm26	\n\t"\
		"vmulpd		%%zmm10,%%zmm5,%%zmm1		\n\t		vmulpd			%%zmm7 ,%%zmm9 ,%%zmm9	\n\t	vfmsub132pd	%%zmm21,%%zmm17,%%zmm20		\n\t	vfmsub132pd			%%zmm23,%%zmm25,%%zmm28	\n\t"\
		"vmulpd		%%zmm13,%%zmm3,%%zmm3		\n\t		vmulpd			%%zmm10,%%zmm11,%%zmm11	\n\t	vfmadd132pd	%%zmm23,%%zmm19,%%zmm22		\n\t	vfmsub132pd			%%zmm31,%%zmm27,%%zmm30	\n\t"\
		"vmulpd		%%zmm10,%%zmm4,%%zmm0		\n\t		vmulpd			%%zmm7 ,%%zmm12,%%zmm8	\n\t	vfmadd132pd	0x040(%%r12),%%zmm16,%%zmm21	\n\t	vfmadd132pd			%%zmm23,%%zmm24,%%zmm29	\n\t"\
		"vmulpd		%%zmm13,%%zmm6,%%zmm2		\n\t		vmulpd		0x200(%%rdx),%%zmm10,%%zmm10\n\t	vfmsub132pd	0x040(%%r13),%%zmm18,%%zmm23	\n\t	vfmadd132pd		0x240(%%r13),%%zmm26,%%zmm31\n\t"\
	"vfmsub132pd	%%zmm15,%%zmm1,%%zmm4		\n\t	vfmsub132pd			%%zmm13,%%zmm9 ,%%zmm12	\n\t		vsubpd		%%zmm22,%%zmm20,%%zmm20	\n\t		vsubpd		%%zmm30,%%zmm28,%%zmm28\n\t"\
	"vfmsub132pd	%%zmm7 ,%%zmm3,%%zmm6		\n\t	vfmadd132pd			%%zmm15,%%zmm11,%%zmm14	\n\t		vsubpd		%%zmm23,%%zmm21,%%zmm21	\n\t		vsubpd		%%zmm31,%%zmm29,%%zmm29\n\t"\
	"vfmadd132pd	%%zmm15,%%zmm0,%%zmm5		\n\t	vfmadd132pd		0x240(%%rcx),%%zmm8 ,%%zmm13\n\t		vmovaps		 (%%r11),%%zmm18			\n\t		vmovaps		0x200(%%r11),%%zmm26		\n\t"\
	"vfmadd132pd	0x040(%%rdx),%%zmm2,%%zmm7	\n\t	vfmsub132pd		0x240(%%rdx),%%zmm10,%%zmm15\n\t		vmovaps	0x040(%%r11),%%zmm16			\n\t		vmovaps		0x240(%%r11),%%zmm24		\n\t"\
		"vsubpd		%%zmm6,%%zmm4,%%zmm4		\n\t		vsubpd		%%zmm14,%%zmm12,%%zmm12\n\t	vfmadd132pd	(%%rdi),%%zmm20,%%zmm22	\n\t	vfmadd132pd		(%%rdi),%%zmm28,%%zmm30\n\t"\
		"vsubpd		%%zmm7,%%zmm5,%%zmm5		\n\t		vsubpd		%%zmm15,%%zmm13,%%zmm13\n\t	vfmadd132pd	(%%rdi),%%zmm21,%%zmm23	\n\t	vfmadd132pd		(%%rdi),%%zmm29,%%zmm31\n\t"\
		"vmovaps		 (%%rbx),%%zmm2			\n\t		vmovaps		0x200(%%rbx),%%zmm10		\n\t		vmovaps		 (%%r11),%%zmm18			\n\t		vmovaps		0x200(%%r11),%%zmm26		\n\t"\
		"vmovaps	0x040(%%rbx),%%zmm0			\n\t		vmovaps		0x240(%%rbx),%%zmm8			\n\t		vmovaps	0x040(%%r11),%%zmm16			\n\t		vmovaps		0x240(%%r11),%%zmm24		\n\t"\
	"vfmadd132pd	(%%rdi),%%zmm4,%%zmm6		\n\t	vfmadd132pd		(%%rdi),%%zmm12,%%zmm14\n\t		vmovaps		 (%%r11),%%zmm17			\n\t		vmovaps		0x200(%%r11),%%zmm25		\n\t"\
	"vfmadd132pd	(%%rdi),%%zmm5,%%zmm7		\n\t	vfmadd132pd		(%%rdi),%%zmm13,%%zmm15\n\t		vmovaps	0x040(%%r11),%%zmm19			\n\t		vmovaps		0x240(%%r11),%%zmm27		\n\t"\
		"vmovaps		 (%%rbx),%%zmm1			\n\t		vmovaps		0x200(%%rbx),%%zmm9			\n\t		vmovaps		 (%%r11),%%zmm17			\n\t		vmovaps		0x200(%%r11),%%zmm25		\n\t"\
		"vmovaps	-0x80(%%rsi),%%zmm3			\n\t		vmovaps		-0x40(%%rsi),%%zmm11		\n\t		vmovaps	-0x40(%%rsi),%%zmm19			\n\t		vmovaps		-0x80(%%rsi),%%zmm27		\n\t"\
		"vmulpd		%%zmm11,%%zmm0,%%zmm0		\n\t		vmulpd			%%zmm3 ,%%zmm8 ,%%zmm8	\n\t		vmulpd		%%zmm27,%%zmm16,%%zmm16	\n\t		vmulpd			%%zmm19,%%zmm24,%%zmm24	\n\t"\
		"vmulpd		%%zmm11,%%zmm1,%%zmm1		\n\t		vmulpd			%%zmm3 ,%%zmm9 ,%%zmm9	\n\t		vmulpd		%%zmm27,%%zmm17,%%zmm17	\n\t		vmulpd			%%zmm19,%%zmm25,%%zmm25	\n\t"\
	"vfmsub132pd	%%zmm3 ,%%zmm0,%%zmm2		\n\t	vfmadd132pd			%%zmm11,%%zmm8 ,%%zmm10	\n\t	vfmsub132pd		%%zmm19,%%zmm16,%%zmm18		\n\t	vfmadd132pd			%%zmm27,%%zmm24,%%zmm26	\n\t"\
	"vfmadd132pd	0x40(%%rbx),%%zmm1,%%zmm3	\n\t	vfmsub132pd		0x240(%%rbx),%%zmm9,%%zmm11	\n\t	vfmadd132pd	0x40(%%r11),%%zmm17,%%zmm19		\n\t	vfmsub132pd		0x240(%%r11),%%zmm25,%%zmm27\n\t"\
		"vmovaps		 (%%rax),%%zmm0			\n\t		vmovaps		0x200(%%rax),%%zmm8			\n\t		vmovaps		 (%%r10),%%zmm16			\n\t		vmovaps		0x200(%%r10),%%zmm24		\n\t"\
		"vmovaps	0x040(%%rax),%%zmm1			\n\t		vmovaps		0x240(%%rax),%%zmm9			\n\t		vmovaps	0x040(%%r10),%%zmm17			\n\t		vmovaps		0x240(%%r10),%%zmm25		\n\t"\
		"vsubpd		%%zmm2,%%zmm0,%%zmm0		\n\t		vsubpd		%%zmm10,%%zmm8,%%zmm8	\n\t		vsubpd		%%zmm18,%%zmm16,%%zmm16	\n\t		vsubpd		%%zmm26,%%zmm24,%%zmm24\n\t"\
		"vsubpd		%%zmm3,%%zmm1,%%zmm1		\n\t		vsubpd		%%zmm11,%%zmm9,%%zmm9	\n\t		vsubpd		%%zmm19,%%zmm17,%%zmm17	\n\t		vsubpd		%%zmm27,%%zmm25,%%zmm25\n\t"\
		"vaddpd			 (%%rax),%%zmm2,%%zmm2	\n\t		vaddpd		0x200(%%rax),%%zmm10,%%zmm10\n\t		vaddpd			 (%%r10),%%zmm18,%%zmm18\n\t		vaddpd		0x200(%%r10),%%zmm26,%%zmm26\n\t"\
		"vaddpd		0x040(%%rax),%%zmm3,%%zmm3	\n\t		vaddpd		0x240(%%rax),%%zmm11,%%zmm11\n\t		vaddpd		0x040(%%r10),%%zmm19,%%zmm19\n\t		vaddpd		0x240(%%r10),%%zmm27,%%zmm27\n\t"\
		"vsubpd		%%zmm6,%%zmm2,%%zmm2		\n\t		vsubpd		%%zmm12,%%zmm8,%%zmm8	\n\t		vsubpd		%%zmm20,%%zmm18,%%zmm18	\n\t		vsubpd		%%zmm28,%%zmm24,%%zmm24\n\t"\
		"vsubpd		%%zmm7,%%zmm3,%%zmm3		\n\t		vsubpd		%%zmm13,%%zmm9,%%zmm9	\n\t		vsubpd		%%zmm21,%%zmm19,%%zmm19	\n\t		vsubpd		%%zmm29,%%zmm25,%%zmm25\n\t"\
	"vfmadd132pd	(%%rdi),%%zmm2,%%zmm6		\n\t	vfmadd132pd		(%%rdi),%%zmm8,%%zmm12\n\t	vfmadd132pd	(%%rdi),%%zmm18,%%zmm20	\n\t	vfmadd132pd		(%%rdi),%%zmm24,%%zmm28\n\t"\
	"vfmadd132pd	(%%rdi),%%zmm3,%%zmm7		\n\t	vfmadd132pd		(%%rdi),%%zmm9,%%zmm13\n\t	vfmadd132pd	(%%rdi),%%zmm19,%%zmm21	\n\t	vfmadd132pd		(%%rdi),%%zmm25,%%zmm29\n\t"\
		"vmovaps		%%zmm2,      (%%rcx)	\n\t		vmovaps		%%zmm8 ,0x200(%%rcx)		\n\t		vmovaps		%%zmm18,      (%%r12)		\n\t		vmovaps		%%zmm24,0x200(%%r12)		\n\t"\
		"vmovaps		%%zmm3, 0x040(%%rcx)	\n\t		vmovaps		%%zmm9 ,0x240(%%rcx)		\n\t		vmovaps		%%zmm19, 0x040(%%r12)		\n\t		vmovaps		%%zmm25,0x240(%%r12)		\n\t"\
		"vmovaps		%%zmm6,      (%%rax)	\n\t		vmovaps		%%zmm12,0x200(%%rax)		\n\t		vmovaps		%%zmm20,      (%%r10)		\n\t		vmovaps		%%zmm28,0x200(%%r10)		\n\t"\
		"vmovaps		%%zmm7, 0x040(%%rax)	\n\t		vmovaps		%%zmm13,0x240(%%rax)		\n\t		vmovaps		%%zmm21, 0x040(%%r10)		\n\t		vmovaps		%%zmm29,0x240(%%r10)		\n\t"\
		"vsubpd		%%zmm5,%%zmm0,%%zmm0		\n\t		vsubpd		%%zmm15,%%zmm10,%%zmm10\n\t		vsubpd		%%zmm23,%%zmm16,%%zmm16	\n\t		vsubpd		%%zmm31,%%zmm26,%%zmm26\n\t"\
		"vsubpd		%%zmm4,%%zmm1,%%zmm1		\n\t		vsubpd		%%zmm14,%%zmm11,%%zmm11\n\t		vsubpd		%%zmm22,%%zmm17,%%zmm17	\n\t		vsubpd		%%zmm30,%%zmm27,%%zmm27\n\t"\
	"vfmadd132pd	(%%rdi),%%zmm0,%%zmm5		\n\t	vfmadd132pd		(%%rdi),%%zmm10,%%zmm15\n\t	vfmadd132pd	(%%rdi),%%zmm16,%%zmm23	\n\t	vfmadd132pd		(%%rdi),%%zmm26,%%zmm31\n\t"\
	"vfmadd132pd	(%%rdi),%%zmm1,%%zmm4		\n\t	vfmadd132pd		(%%rdi),%%zmm11,%%zmm14\n\t	vfmadd132pd	(%%rdi),%%zmm17,%%zmm22	\n\t	vfmadd132pd		(%%rdi),%%zmm27,%%zmm30\n\t"\
		"vmovaps		%%zmm0,      (%%rbx)	\n\t		vmovaps		%%zmm10,0x200(%%rbx)		\n\t		vmovaps		%%zmm16,      (%%r11)		\n\t		vmovaps		%%zmm26,0x200(%%r11)		\n\t"\
		"vmovaps		%%zmm1, 0x040(%%rdx)	\n\t		vmovaps		%%zmm11,0x240(%%rdx)		\n\t		vmovaps		%%zmm17, 0x040(%%r13)		\n\t		vmovaps		%%zmm27,0x240(%%r13)		\n\t"\
		"vmovaps		%%zmm5,      (%%rdx)	\n\t		vmovaps		%%zmm15,0x200(%%rdx)		\n\t		vmovaps		%%zmm23,      (%%r13)		\n\t		vmovaps		%%zmm31,0x200(%%r13)		\n\t"\
		"vmovaps		%%zmm4, 0x040(%%rbx)	\n\t		vmovaps		%%zmm14,0x240(%%rbx)		\n\t		vmovaps		%%zmm22, 0x040(%%r11)		\n\t		vmovaps		%%zmm30,0x240(%%r11)		\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__add1] "m" (Xadd1)\
		 ,[__add2] "m" (Xadd2)\
		 ,[__add3] "m" (Xadd3)\
		 ,[__add4] "m" (Xadd4)\
		 ,[__add5] "m" (Xadd5)\
		 ,[__add6] "m" (Xadd6)\
		 ,[__add7] "m" (Xadd7)\
		 ,[__r00] "m" (Xr00)\
		 ,[__r10] "m" (Xr10)\
		 ,[__r20] "m" (Xr20)\
		 ,[__r30] "m" (Xr30)\
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__cc0] "m" (Xcc0)\
		 ,[__c00] "m" (Xc00)\
		 ,[__c01] "m" (Xc01)\
		 ,[__c02] "m" (Xc02)\
		 ,[__c03] "m" (Xc03)\
		 ,[__c05] "m" (Xc05)\
		 ,[__c07] "m" (Xc07)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r10","r11","r12","r13","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15","xmm16","xmm17","xmm18","xmm19","xmm20","xmm21","xmm22","xmm23","xmm24","xmm25","xmm26","xmm27","xmm28","xmm29","xmm30","xmm31"	/* Clobbered registers */\
	);\
	}

   #else	// AVX-512 version:

	// Cost [vector-ops only]: 740 MEM (128 in transpose section, 610 in DFT proper), 540 ARITHMETIC
	#define SSE2_RADIX32_WRAPPER_DIF(Xadd0,Xadd1,Xadd2,Xadd3,Xadd4,Xadd5,Xadd6,Xadd7,Xr00,Xr10,Xr20,Xr30,Xisrt2,Xcc0,Xc00,Xc01,Xc02,Xc03,Xc05,Xc07)\
	{\
	__asm__ volatile (\
		"movl	$0b00110011,%%eax	\n\t"/* Constant for vblendmpd instructions goes into mask-reg k1 */\
		"kmovw	%%eax,%%k1			\n\t"\
		/* Init vector index-consts needed by vpermt2pd instructions - if regs were at a premium,
		could also init just prior to [3] and use zmm6,7 to hold index-consts: */\
		"movq	$0x0c040e0608000a02,%%rax	\n\t"/* zmm30 = 8+4 0+4 8+6 0+6 8+0 0+0 8+2 0+2 [msw at left] */\
		"movq	$0x0d050f0709010b03,%%rbx	\n\t"/* zmm31 = 8+5 0+5 8+7 0+7 8+1 0+1 8+3 0+3 */\
			"vmovq		%%rax,%%xmm0 		\n\t"\
			"vmovq		%%rbx,%%xmm1 		\n\t"\
			"vpmovzxbq	%%xmm0,%%zmm30		\n\t"\
			"vpmovzxbq	%%xmm1,%%zmm31		\n\t"\
	/**** Start with 8-way interleaving: ****/\
		"movq	%[__add0],%%rax	\n\t"\
		"movq	%[__add1],%%rbx	\n\t"\
		"movq	%[__add2],%%rcx	\n\t"\
		"movq	%[__add3],%%rdx	\n\t"\
		"movq	%[__add4],%%r10	\n\t"\
		"movq	%[__add5],%%r11	\n\t"\
		"movq	%[__add6],%%r12	\n\t"\
		"movq	%[__add7],%%r13	\n\t"\
		"movq	%[__r00] ,%%rsi	\n\t"\
	/* See above 16-register version for transpose using columnwise gather-loads. */\
	/* a[j+p0]: Inputs from add0+[0,0x200,0x400,0x600,0x800,0xa00,0xc00,0xe00]. Outputs into local store at r0+[same byte offsets]: */\
		/* Read contiguous (untransposed) row-data from the data array: */\
		"vmovaps 	(%%rax),%%zmm0			\n\t	vmovaps 0x40(%%rax),%%zmm12	\n\t"\
		"vmovaps 	(%%rbx),%%zmm1			\n\t	vmovaps 0x40(%%rbx),%%zmm13	\n\t"\
		"vmovaps 	(%%rcx),%%zmm2			\n\t	vmovaps 0x40(%%rcx),%%zmm14	\n\t"\
		"vmovaps 	(%%rdx),%%zmm3			\n\t	vmovaps 0x40(%%rdx),%%zmm15	\n\t"\
		"vmovaps 	(%%r10),%%zmm4			\n\t	vmovaps 0x40(%%r10),%%zmm16	\n\t"\
		"vmovaps 	(%%r11),%%zmm5			\n\t	vmovaps 0x40(%%r11),%%zmm17	\n\t"\
		"vmovaps 	(%%r12),%%zmm6			\n\t	vmovaps 0x40(%%r12),%%zmm18	\n\t"\
		"vmovaps 	(%%r13),%%zmm7			\n\t	vmovaps 0x40(%%r13),%%zmm19	\n\t"\
		/* [1] Shuffle the 4-aparts - note the different patterning of the first and second output quartet: */\
		"vshuff64x2	$0b01000100, %%zmm4,	%%zmm0,	%%zmm8 	\n\t	vshuff64x2	$0b01000100, %%zmm16,%%zmm12,	%%zmm20	\n\t"\
		"vshuff64x2	$0b11101110, %%zmm4,	%%zmm0,	%%zmm4 	\n\t	vshuff64x2	$0b11101110, %%zmm16,%%zmm12,	%%zmm16	\n\t"\
		"vshuff64x2	$0b01000100, %%zmm5,	%%zmm1,	%%zmm9	\n\t	vshuff64x2	$0b01000100, %%zmm17,%%zmm13,	%%zmm21	\n\t"\
		"vshuff64x2	$0b11101110, %%zmm5,	%%zmm1,	%%zmm5	\n\t	vshuff64x2	$0b11101110, %%zmm17,%%zmm13,	%%zmm17	\n\t"\
		"vshuff64x2	$0b00010001, %%zmm6,	%%zmm2,	%%zmm10	\n\t	vshuff64x2	$0b00010001, %%zmm18,%%zmm14,	%%zmm22	\n\t"\
		"vshuff64x2	$0b10111011, %%zmm6,	%%zmm2,	%%zmm6	\n\t	vshuff64x2	$0b10111011, %%zmm18,%%zmm14,	%%zmm18	\n\t"\
		"vshuff64x2	$0b00010001, %%zmm7,	%%zmm3,	%%zmm11	\n\t	vshuff64x2	$0b00010001, %%zmm19,%%zmm15,	%%zmm23	\n\t"\
		"vshuff64x2	$0b10111011, %%zmm7,	%%zmm3,	%%zmm7	\n\t	vshuff64x2	$0b10111011, %%zmm19,%%zmm15,	%%zmm19	\n\t"\
		/* [2] Blend in the 2-aparts */\
		"vblendmpd	%%zmm8 ,	%%zmm10,	%%zmm0%{%%k1%}	\n\t	vblendmpd	%%zmm20,	%%zmm22,	%%zmm12%{%%k1%}	\n\t"\
		"vblendmpd	%%zmm10,	%%zmm8 ,	%%zmm8%{%%k1%}	\n\t	vblendmpd	%%zmm22,	%%zmm20,	%%zmm20%{%%k1%}	\n\t"\
		"vblendmpd	%%zmm4 ,	%%zmm6 ,	%%zmm1%{%%k1%}	\n\t	vblendmpd	%%zmm16,	%%zmm18,	%%zmm13%{%%k1%}	\n\t"\
		"vblendmpd	%%zmm6 ,	%%zmm4 ,	%%zmm4%{%%k1%}	\n\t	vblendmpd	%%zmm18,	%%zmm16,	%%zmm16%{%%k1%}	\n\t"\
		"vblendmpd	%%zmm9 ,	%%zmm11,	%%zmm2%{%%k1%}	\n\t	vblendmpd	%%zmm21,	%%zmm23,	%%zmm14%{%%k1%}	\n\t"\
		"vblendmpd	%%zmm11,	%%zmm9 ,	%%zmm9%{%%k1%}	\n\t	vblendmpd	%%zmm23,	%%zmm21,	%%zmm21%{%%k1%}	\n\t"\
		"vblendmpd	%%zmm5 ,	%%zmm7 ,	%%zmm3%{%%k1%}	\n\t	vblendmpd	%%zmm17,	%%zmm19,	%%zmm15%{%%k1%}	\n\t"\
		"vblendmpd	%%zmm7 ,	%%zmm5 ,	%%zmm5%{%%k1%}	\n\t	vblendmpd	%%zmm19,	%%zmm17,	%%zmm17%{%%k1%}	\n\t"\
		/* [3] Shuffle or permute in the 1-aparts */\
		"vshufpd	$0b00000000,%%zmm2,		%%zmm0,%%zmm10 	\n\t	vshufpd	$0b00000000,%%zmm14,	%%zmm12,%%zmm22	\n\t"\
		"vshufpd	$0b11111111,%%zmm2,		%%zmm0,%%zmm11 	\n\t	vshufpd	$0b11111111,%%zmm14,	%%zmm12,%%zmm23	\n\t"\
		"vmovapd	%%zmm8,%%zmm2							\n\t	vmovapd	%%zmm20,%%zmm14	\n\t"\
		"vpermt2pd				%%zmm9,		%%zmm30,%%zmm2 	\n\t	vpermt2pd			%%zmm21,	%%zmm30,%%zmm14	\n\t"\
		"vpermt2pd				%%zmm9,		%%zmm31,%%zmm8	\n\t	vpermt2pd			%%zmm21,	%%zmm31,%%zmm20	\n\t"\
		"vshufpd	$0b00000000,%%zmm3,		%%zmm1,%%zmm0 	\n\t	vshufpd	$0b00000000,%%zmm15,	%%zmm13,%%zmm12	\n\t"\
		"vshufpd	$0b11111111,%%zmm3,		%%zmm1,%%zmm1 	\n\t	vshufpd	$0b11111111,%%zmm15,	%%zmm13,%%zmm13	\n\t"\
		"vmovapd	%%zmm4,%%zmm3							\n\t	vmovapd	%%zmm16,%%zmm15	\n\t"\
		"vpermt2pd				%%zmm5,		%%zmm30,%%zmm3 	\n\t	vpermt2pd			%%zmm17,	%%zmm30,%%zmm15	\n\t"\
		"vpermt2pd				%%zmm5,		%%zmm31,%%zmm4	\n\t	vpermt2pd			%%zmm17,	%%zmm31,%%zmm16	\n\t"\
		/* Write transposed data to the local store: */\
		"vmovaps		%%zmm10,0x000(%%rsi)		\n\t	vmovaps	%%zmm22,0x040(%%rsi)	\n\t"/* Out 0 */\
		"vmovaps		%%zmm11,0x800(%%rsi)		\n\t	vmovaps	%%zmm23,0x840(%%rsi)	\n\t"/* Out 4 */\
		"vmovaps		%%zmm2 ,0x400(%%rsi)		\n\t	vmovaps	%%zmm14,0x440(%%rsi)	\n\t"/* Out 2 */\
		"vmovaps		%%zmm8 ,0xc00(%%rsi)		\n\t	vmovaps	%%zmm20,0xc40(%%rsi)	\n\t"/* Out 6 */\
		"vmovaps		%%zmm0 ,0x200(%%rsi)		\n\t	vmovaps	%%zmm12,0x240(%%rsi)	\n\t"/* Out 1 */\
		"vmovaps		%%zmm1 ,0xa00(%%rsi)		\n\t	vmovaps	%%zmm13,0xa40(%%rsi)	\n\t"/* Out 5 */\
		"vmovaps		%%zmm3 ,0x600(%%rsi)		\n\t	vmovaps	%%zmm15,0x640(%%rsi)	\n\t"/* Out 3 */\
		"vmovaps		%%zmm4, 0xe00(%%rsi)		\n\t	vmovaps	%%zmm16,0xe40(%%rsi)	\n\t"/* Out 7 */\
		"\n\t"\
	/* a[j+p2]: Inputs from add0+[0,0x100,0x200,0x300,0x400,0x500,0x600,0x700]+0x80. Outputs into r1+[same byte offsets]: */\
	"addq	$0x80 ,%%rax	\n\t"\
	"addq	$0x80 ,%%rbx	\n\t"\
	"addq	$0x80 ,%%rcx	\n\t"\
	"addq	$0x80 ,%%rdx	\n\t"\
	"addq	$0x80 ,%%r10	\n\t"\
	"addq	$0x80 ,%%r11	\n\t"\
	"addq	$0x80 ,%%r12	\n\t"\
	"addq	$0x80 ,%%r13	\n\t"\
	"addq	$0x100,%%rsi	\n\t"\
		/* Read contiguous (untransposed) row-data from the data array: */\
		"vmovaps 	(%%rax),%%zmm0			\n\t	vmovaps 0x40(%%rax),%%zmm12	\n\t"\
		"vmovaps 	(%%rbx),%%zmm1			\n\t	vmovaps 0x40(%%rbx),%%zmm13	\n\t"\
		"vmovaps 	(%%rcx),%%zmm2			\n\t	vmovaps 0x40(%%rcx),%%zmm14	\n\t"\
		"vmovaps 	(%%rdx),%%zmm3			\n\t	vmovaps 0x40(%%rdx),%%zmm15	\n\t"\
		"vmovaps 	(%%r10),%%zmm4			\n\t	vmovaps 0x40(%%r10),%%zmm16	\n\t"\
		"vmovaps 	(%%r11),%%zmm5			\n\t	vmovaps 0x40(%%r11),%%zmm17	\n\t"\
		"vmovaps 	(%%r12),%%zmm6			\n\t	vmovaps 0x40(%%r12),%%zmm18	\n\t"\
		"vmovaps 	(%%r13),%%zmm7			\n\t	vmovaps 0x40(%%r13),%%zmm19	\n\t"\
		/* [1] Shuffle the 4-aparts - note the different patterning of the first and second output quartet: */\
		"vshuff64x2	$0b01000100, %%zmm4,	%%zmm0,	%%zmm8 	\n\t	vshuff64x2	$0b01000100, %%zmm16,%%zmm12,	%%zmm20	\n\t"\
		"vshuff64x2	$0b11101110, %%zmm4,	%%zmm0,	%%zmm4 	\n\t	vshuff64x2	$0b11101110, %%zmm16,%%zmm12,	%%zmm16	\n\t"\
		"vshuff64x2	$0b01000100, %%zmm5,	%%zmm1,	%%zmm9	\n\t	vshuff64x2	$0b01000100, %%zmm17,%%zmm13,	%%zmm21	\n\t"\
		"vshuff64x2	$0b11101110, %%zmm5,	%%zmm1,	%%zmm5	\n\t	vshuff64x2	$0b11101110, %%zmm17,%%zmm13,	%%zmm17	\n\t"\
		"vshuff64x2	$0b00010001, %%zmm6,	%%zmm2,	%%zmm10	\n\t	vshuff64x2	$0b00010001, %%zmm18,%%zmm14,	%%zmm22	\n\t"\
		"vshuff64x2	$0b10111011, %%zmm6,	%%zmm2,	%%zmm6	\n\t	vshuff64x2	$0b10111011, %%zmm18,%%zmm14,	%%zmm18	\n\t"\
		"vshuff64x2	$0b00010001, %%zmm7,	%%zmm3,	%%zmm11	\n\t	vshuff64x2	$0b00010001, %%zmm19,%%zmm15,	%%zmm23	\n\t"\
		"vshuff64x2	$0b10111011, %%zmm7,	%%zmm3,	%%zmm7	\n\t	vshuff64x2	$0b10111011, %%zmm19,%%zmm15,	%%zmm19	\n\t"\
		/* [2] Blend in the 2-aparts */\
		"vblendmpd	%%zmm8 ,	%%zmm10,	%%zmm0%{%%k1%}	\n\t	vblendmpd	%%zmm20,	%%zmm22,	%%zmm12%{%%k1%}	\n\t"\
		"vblendmpd	%%zmm10,	%%zmm8 ,	%%zmm8%{%%k1%}	\n\t	vblendmpd	%%zmm22,	%%zmm20,	%%zmm20%{%%k1%}	\n\t"\
		"vblendmpd	%%zmm4 ,	%%zmm6 ,	%%zmm1%{%%k1%}	\n\t	vblendmpd	%%zmm16,	%%zmm18,	%%zmm13%{%%k1%}	\n\t"\
		"vblendmpd	%%zmm6 ,	%%zmm4 ,	%%zmm4%{%%k1%}	\n\t	vblendmpd	%%zmm18,	%%zmm16,	%%zmm16%{%%k1%}	\n\t"\
		"vblendmpd	%%zmm9 ,	%%zmm11,	%%zmm2%{%%k1%}	\n\t	vblendmpd	%%zmm21,	%%zmm23,	%%zmm14%{%%k1%}	\n\t"\
		"vblendmpd	%%zmm11,	%%zmm9 ,	%%zmm9%{%%k1%}	\n\t	vblendmpd	%%zmm23,	%%zmm21,	%%zmm21%{%%k1%}	\n\t"\
		"vblendmpd	%%zmm5 ,	%%zmm7 ,	%%zmm3%{%%k1%}	\n\t	vblendmpd	%%zmm17,	%%zmm19,	%%zmm15%{%%k1%}	\n\t"\
		"vblendmpd	%%zmm7 ,	%%zmm5 ,	%%zmm5%{%%k1%}	\n\t	vblendmpd	%%zmm19,	%%zmm17,	%%zmm17%{%%k1%}	\n\t"\
		/* [3] Shuffle or permute in the 1-aparts */\
		"vshufpd	$0b00000000,%%zmm2,		%%zmm0,%%zmm10 	\n\t	vshufpd	$0b00000000,%%zmm14,	%%zmm12,%%zmm22	\n\t"\
		"vshufpd	$0b11111111,%%zmm2,		%%zmm0,%%zmm11 	\n\t	vshufpd	$0b11111111,%%zmm14,	%%zmm12,%%zmm23	\n\t"\
		"vmovapd	%%zmm8,%%zmm2							\n\t	vmovapd	%%zmm20,%%zmm14	\n\t"\
		"vpermt2pd				%%zmm9,		%%zmm30,%%zmm2 	\n\t	vpermt2pd			%%zmm21,	%%zmm30,%%zmm14	\n\t"\
		"vpermt2pd				%%zmm9,		%%zmm31,%%zmm8	\n\t	vpermt2pd			%%zmm21,	%%zmm31,%%zmm20	\n\t"\
		"vshufpd	$0b00000000,%%zmm3,		%%zmm1,%%zmm0 	\n\t	vshufpd	$0b00000000,%%zmm15,	%%zmm13,%%zmm12	\n\t"\
		"vshufpd	$0b11111111,%%zmm3,		%%zmm1,%%zmm1 	\n\t	vshufpd	$0b11111111,%%zmm15,	%%zmm13,%%zmm13	\n\t"\
		"vmovapd	%%zmm4,%%zmm3							\n\t	vmovapd	%%zmm16,%%zmm15	\n\t"\
		"vpermt2pd				%%zmm5,		%%zmm30,%%zmm3 	\n\t	vpermt2pd			%%zmm17,	%%zmm30,%%zmm15	\n\t"\
		"vpermt2pd				%%zmm5,		%%zmm31,%%zmm4	\n\t	vpermt2pd			%%zmm17,	%%zmm31,%%zmm16	\n\t"\
		/* Write transposed data to the local store: */\
		"vmovaps		%%zmm10,0x000(%%rsi)		\n\t	vmovaps	%%zmm22,0x040(%%rsi)	\n\t"/* Out 0 */\
		"vmovaps		%%zmm11,0x800(%%rsi)		\n\t	vmovaps	%%zmm23,0x840(%%rsi)	\n\t"/* Out 4 */\
		"vmovaps		%%zmm2 ,0x400(%%rsi)		\n\t	vmovaps	%%zmm14,0x440(%%rsi)	\n\t"/* Out 2 */\
		"vmovaps		%%zmm8 ,0xc00(%%rsi)		\n\t	vmovaps	%%zmm20,0xc40(%%rsi)	\n\t"/* Out 6 */\
		"vmovaps		%%zmm0 ,0x200(%%rsi)		\n\t	vmovaps	%%zmm12,0x240(%%rsi)	\n\t"/* Out 1 */\
		"vmovaps		%%zmm1 ,0xa00(%%rsi)		\n\t	vmovaps	%%zmm13,0xa40(%%rsi)	\n\t"/* Out 5 */\
		"vmovaps		%%zmm3 ,0x600(%%rsi)		\n\t	vmovaps	%%zmm15,0x640(%%rsi)	\n\t"/* Out 3 */\
		"vmovaps		%%zmm4, 0xe00(%%rsi)		\n\t	vmovaps	%%zmm16,0xe40(%%rsi)	\n\t"/* Out 7 */\
		"\n\t"\
	/* a[j+p4]: Inputs from add0+[0,0x100,0x200,0x300,0x400,0x500,0x600,0x700]+0x80. Outputs into r2+[same byte offsets]: */\
	"addq	$0x80 ,%%rax	\n\t"\
	"addq	$0x80 ,%%rbx	\n\t"\
	"addq	$0x80 ,%%rcx	\n\t"\
	"addq	$0x80 ,%%rdx	\n\t"\
	"addq	$0x80 ,%%r10	\n\t"\
	"addq	$0x80 ,%%r11	\n\t"\
	"addq	$0x80 ,%%r12	\n\t"\
	"addq	$0x80 ,%%r13	\n\t"\
	"subq	$0x80,%%rsi	\n\t"\
		/* Read contiguous (untransposed) row-data from the data array: */\
		"vmovaps 	(%%rax),%%zmm0			\n\t	vmovaps 0x40(%%rax),%%zmm12	\n\t"\
		"vmovaps 	(%%rbx),%%zmm1			\n\t	vmovaps 0x40(%%rbx),%%zmm13	\n\t"\
		"vmovaps 	(%%rcx),%%zmm2			\n\t	vmovaps 0x40(%%rcx),%%zmm14	\n\t"\
		"vmovaps 	(%%rdx),%%zmm3			\n\t	vmovaps 0x40(%%rdx),%%zmm15	\n\t"\
		"vmovaps 	(%%r10),%%zmm4			\n\t	vmovaps 0x40(%%r10),%%zmm16	\n\t"\
		"vmovaps 	(%%r11),%%zmm5			\n\t	vmovaps 0x40(%%r11),%%zmm17	\n\t"\
		"vmovaps 	(%%r12),%%zmm6			\n\t	vmovaps 0x40(%%r12),%%zmm18	\n\t"\
		"vmovaps 	(%%r13),%%zmm7			\n\t	vmovaps 0x40(%%r13),%%zmm19	\n\t"\
		/* [1] Shuffle the 4-aparts - note the different patterning of the first and second output quartet: */\
		"vshuff64x2	$0b01000100, %%zmm4,	%%zmm0,	%%zmm8 	\n\t	vshuff64x2	$0b01000100, %%zmm16,%%zmm12,	%%zmm20	\n\t"\
		"vshuff64x2	$0b11101110, %%zmm4,	%%zmm0,	%%zmm4 	\n\t	vshuff64x2	$0b11101110, %%zmm16,%%zmm12,	%%zmm16	\n\t"\
		"vshuff64x2	$0b01000100, %%zmm5,	%%zmm1,	%%zmm9	\n\t	vshuff64x2	$0b01000100, %%zmm17,%%zmm13,	%%zmm21	\n\t"\
		"vshuff64x2	$0b11101110, %%zmm5,	%%zmm1,	%%zmm5	\n\t	vshuff64x2	$0b11101110, %%zmm17,%%zmm13,	%%zmm17	\n\t"\
		"vshuff64x2	$0b00010001, %%zmm6,	%%zmm2,	%%zmm10	\n\t	vshuff64x2	$0b00010001, %%zmm18,%%zmm14,	%%zmm22	\n\t"\
		"vshuff64x2	$0b10111011, %%zmm6,	%%zmm2,	%%zmm6	\n\t	vshuff64x2	$0b10111011, %%zmm18,%%zmm14,	%%zmm18	\n\t"\
		"vshuff64x2	$0b00010001, %%zmm7,	%%zmm3,	%%zmm11	\n\t	vshuff64x2	$0b00010001, %%zmm19,%%zmm15,	%%zmm23	\n\t"\
		"vshuff64x2	$0b10111011, %%zmm7,	%%zmm3,	%%zmm7	\n\t	vshuff64x2	$0b10111011, %%zmm19,%%zmm15,	%%zmm19	\n\t"\
		/* [2] Blend in the 2-aparts */\
		"vblendmpd	%%zmm8 ,	%%zmm10,	%%zmm0%{%%k1%}	\n\t	vblendmpd	%%zmm20,	%%zmm22,	%%zmm12%{%%k1%}	\n\t"\
		"vblendmpd	%%zmm10,	%%zmm8 ,	%%zmm8%{%%k1%}	\n\t	vblendmpd	%%zmm22,	%%zmm20,	%%zmm20%{%%k1%}	\n\t"\
		"vblendmpd	%%zmm4 ,	%%zmm6 ,	%%zmm1%{%%k1%}	\n\t	vblendmpd	%%zmm16,	%%zmm18,	%%zmm13%{%%k1%}	\n\t"\
		"vblendmpd	%%zmm6 ,	%%zmm4 ,	%%zmm4%{%%k1%}	\n\t	vblendmpd	%%zmm18,	%%zmm16,	%%zmm16%{%%k1%}	\n\t"\
		"vblendmpd	%%zmm9 ,	%%zmm11,	%%zmm2%{%%k1%}	\n\t	vblendmpd	%%zmm21,	%%zmm23,	%%zmm14%{%%k1%}	\n\t"\
		"vblendmpd	%%zmm11,	%%zmm9 ,	%%zmm9%{%%k1%}	\n\t	vblendmpd	%%zmm23,	%%zmm21,	%%zmm21%{%%k1%}	\n\t"\
		"vblendmpd	%%zmm5 ,	%%zmm7 ,	%%zmm3%{%%k1%}	\n\t	vblendmpd	%%zmm17,	%%zmm19,	%%zmm15%{%%k1%}	\n\t"\
		"vblendmpd	%%zmm7 ,	%%zmm5 ,	%%zmm5%{%%k1%}	\n\t	vblendmpd	%%zmm19,	%%zmm17,	%%zmm17%{%%k1%}	\n\t"\
		/* [3] Shuffle or permute in the 1-aparts */\
		"vshufpd	$0b00000000,%%zmm2,		%%zmm0,%%zmm10 	\n\t	vshufpd	$0b00000000,%%zmm14,	%%zmm12,%%zmm22	\n\t"\
		"vshufpd	$0b11111111,%%zmm2,		%%zmm0,%%zmm11 	\n\t	vshufpd	$0b11111111,%%zmm14,	%%zmm12,%%zmm23	\n\t"\
		"vmovapd	%%zmm8,%%zmm2							\n\t	vmovapd	%%zmm20,%%zmm14	\n\t"\
		"vpermt2pd				%%zmm9,		%%zmm30,%%zmm2 	\n\t	vpermt2pd			%%zmm21,	%%zmm30,%%zmm14	\n\t"\
		"vpermt2pd				%%zmm9,		%%zmm31,%%zmm8	\n\t	vpermt2pd			%%zmm21,	%%zmm31,%%zmm20	\n\t"\
		"vshufpd	$0b00000000,%%zmm3,		%%zmm1,%%zmm0 	\n\t	vshufpd	$0b00000000,%%zmm15,	%%zmm13,%%zmm12	\n\t"\
		"vshufpd	$0b11111111,%%zmm3,		%%zmm1,%%zmm1 	\n\t	vshufpd	$0b11111111,%%zmm15,	%%zmm13,%%zmm13	\n\t"\
		"vmovapd	%%zmm4,%%zmm3							\n\t	vmovapd	%%zmm16,%%zmm15	\n\t"\
		"vpermt2pd				%%zmm5,		%%zmm30,%%zmm3 	\n\t	vpermt2pd			%%zmm17,	%%zmm30,%%zmm15	\n\t"\
		"vpermt2pd				%%zmm5,		%%zmm31,%%zmm4	\n\t	vpermt2pd			%%zmm17,	%%zmm31,%%zmm16	\n\t"\
		/* Write transposed data to the local store: */\
		"vmovaps		%%zmm10,0x000(%%rsi)		\n\t	vmovaps	%%zmm22,0x040(%%rsi)	\n\t"/* Out 0 */\
		"vmovaps		%%zmm11,0x800(%%rsi)		\n\t	vmovaps	%%zmm23,0x840(%%rsi)	\n\t"/* Out 4 */\
		"vmovaps		%%zmm2 ,0x400(%%rsi)		\n\t	vmovaps	%%zmm14,0x440(%%rsi)	\n\t"/* Out 2 */\
		"vmovaps		%%zmm8 ,0xc00(%%rsi)		\n\t	vmovaps	%%zmm20,0xc40(%%rsi)	\n\t"/* Out 6 */\
		"vmovaps		%%zmm0 ,0x200(%%rsi)		\n\t	vmovaps	%%zmm12,0x240(%%rsi)	\n\t"/* Out 1 */\
		"vmovaps		%%zmm1 ,0xa00(%%rsi)		\n\t	vmovaps	%%zmm13,0xa40(%%rsi)	\n\t"/* Out 5 */\
		"vmovaps		%%zmm3 ,0x600(%%rsi)		\n\t	vmovaps	%%zmm15,0x640(%%rsi)	\n\t"/* Out 3 */\
		"vmovaps		%%zmm4, 0xe00(%%rsi)		\n\t	vmovaps	%%zmm16,0xe40(%%rsi)	\n\t"/* Out 7 */\
		"\n\t"\
	/* a[j+p6]: Inputs from add0+[0,0x100,0x200,0x300,0x400,0x500,0x600,0x700]+0x180. Outputs into r3+[same byte offsets]: */\
	"addq	$0x80 ,%%rax	\n\t"\
	"addq	$0x80 ,%%rbx	\n\t"\
	"addq	$0x80 ,%%rcx	\n\t"\
	"addq	$0x80 ,%%rdx	\n\t"\
	"addq	$0x80 ,%%r10	\n\t"\
	"addq	$0x80 ,%%r11	\n\t"\
	"addq	$0x80 ,%%r12	\n\t"\
	"addq	$0x80 ,%%r13	\n\t"\
	"addq	$0x100,%%rsi	\n\t"\
		/* Read contiguous (untransposed) row-data from the data array: */\
		"vmovaps 	(%%rax),%%zmm0			\n\t	vmovaps 0x40(%%rax),%%zmm12	\n\t"\
		"vmovaps 	(%%rbx),%%zmm1			\n\t	vmovaps 0x40(%%rbx),%%zmm13	\n\t"\
		"vmovaps 	(%%rcx),%%zmm2			\n\t	vmovaps 0x40(%%rcx),%%zmm14	\n\t"\
		"vmovaps 	(%%rdx),%%zmm3			\n\t	vmovaps 0x40(%%rdx),%%zmm15	\n\t"\
		"vmovaps 	(%%r10),%%zmm4			\n\t	vmovaps 0x40(%%r10),%%zmm16	\n\t"\
		"vmovaps 	(%%r11),%%zmm5			\n\t	vmovaps 0x40(%%r11),%%zmm17	\n\t"\
		"vmovaps 	(%%r12),%%zmm6			\n\t	vmovaps 0x40(%%r12),%%zmm18	\n\t"\
		"vmovaps 	(%%r13),%%zmm7			\n\t	vmovaps 0x40(%%r13),%%zmm19	\n\t"\
		/* [1] Shuffle the 4-aparts - note the different patterning of the first and second output quartet: */\
		"vshuff64x2	$0b01000100, %%zmm4,	%%zmm0,	%%zmm8 	\n\t	vshuff64x2	$0b01000100, %%zmm16,%%zmm12,	%%zmm20	\n\t"\
		"vshuff64x2	$0b11101110, %%zmm4,	%%zmm0,	%%zmm4 	\n\t	vshuff64x2	$0b11101110, %%zmm16,%%zmm12,	%%zmm16	\n\t"\
		"vshuff64x2	$0b01000100, %%zmm5,	%%zmm1,	%%zmm9	\n\t	vshuff64x2	$0b01000100, %%zmm17,%%zmm13,	%%zmm21	\n\t"\
		"vshuff64x2	$0b11101110, %%zmm5,	%%zmm1,	%%zmm5	\n\t	vshuff64x2	$0b11101110, %%zmm17,%%zmm13,	%%zmm17	\n\t"\
		"vshuff64x2	$0b00010001, %%zmm6,	%%zmm2,	%%zmm10	\n\t	vshuff64x2	$0b00010001, %%zmm18,%%zmm14,	%%zmm22	\n\t"\
		"vshuff64x2	$0b10111011, %%zmm6,	%%zmm2,	%%zmm6	\n\t	vshuff64x2	$0b10111011, %%zmm18,%%zmm14,	%%zmm18	\n\t"\
		"vshuff64x2	$0b00010001, %%zmm7,	%%zmm3,	%%zmm11	\n\t	vshuff64x2	$0b00010001, %%zmm19,%%zmm15,	%%zmm23	\n\t"\
		"vshuff64x2	$0b10111011, %%zmm7,	%%zmm3,	%%zmm7	\n\t	vshuff64x2	$0b10111011, %%zmm19,%%zmm15,	%%zmm19	\n\t"\
		/* [2] Blend in the 2-aparts */\
		"vblendmpd	%%zmm8 ,	%%zmm10,	%%zmm0%{%%k1%}	\n\t	vblendmpd	%%zmm20,	%%zmm22,	%%zmm12%{%%k1%}	\n\t"\
		"vblendmpd	%%zmm10,	%%zmm8 ,	%%zmm8%{%%k1%}	\n\t	vblendmpd	%%zmm22,	%%zmm20,	%%zmm20%{%%k1%}	\n\t"\
		"vblendmpd	%%zmm4 ,	%%zmm6 ,	%%zmm1%{%%k1%}	\n\t	vblendmpd	%%zmm16,	%%zmm18,	%%zmm13%{%%k1%}	\n\t"\
		"vblendmpd	%%zmm6 ,	%%zmm4 ,	%%zmm4%{%%k1%}	\n\t	vblendmpd	%%zmm18,	%%zmm16,	%%zmm16%{%%k1%}	\n\t"\
		"vblendmpd	%%zmm9 ,	%%zmm11,	%%zmm2%{%%k1%}	\n\t	vblendmpd	%%zmm21,	%%zmm23,	%%zmm14%{%%k1%}	\n\t"\
		"vblendmpd	%%zmm11,	%%zmm9 ,	%%zmm9%{%%k1%}	\n\t	vblendmpd	%%zmm23,	%%zmm21,	%%zmm21%{%%k1%}	\n\t"\
		"vblendmpd	%%zmm5 ,	%%zmm7 ,	%%zmm3%{%%k1%}	\n\t	vblendmpd	%%zmm17,	%%zmm19,	%%zmm15%{%%k1%}	\n\t"\
		"vblendmpd	%%zmm7 ,	%%zmm5 ,	%%zmm5%{%%k1%}	\n\t	vblendmpd	%%zmm19,	%%zmm17,	%%zmm17%{%%k1%}	\n\t"\
		/* [3] Shuffle or permute in the 1-aparts */\
		"vshufpd	$0b00000000,%%zmm2,		%%zmm0,%%zmm10 	\n\t	vshufpd	$0b00000000,%%zmm14,	%%zmm12,%%zmm22	\n\t"\
		"vshufpd	$0b11111111,%%zmm2,		%%zmm0,%%zmm11 	\n\t	vshufpd	$0b11111111,%%zmm14,	%%zmm12,%%zmm23	\n\t"\
		"vmovapd	%%zmm8,%%zmm2							\n\t	vmovapd	%%zmm20,%%zmm14	\n\t"\
		"vpermt2pd				%%zmm9,		%%zmm30,%%zmm2 	\n\t	vpermt2pd			%%zmm21,	%%zmm30,%%zmm14	\n\t"\
		"vpermt2pd				%%zmm9,		%%zmm31,%%zmm8	\n\t	vpermt2pd			%%zmm21,	%%zmm31,%%zmm20	\n\t"\
		"vshufpd	$0b00000000,%%zmm3,		%%zmm1,%%zmm0 	\n\t	vshufpd	$0b00000000,%%zmm15,	%%zmm13,%%zmm12	\n\t"\
		"vshufpd	$0b11111111,%%zmm3,		%%zmm1,%%zmm1 	\n\t	vshufpd	$0b11111111,%%zmm15,	%%zmm13,%%zmm13	\n\t"\
		"vmovapd	%%zmm4,%%zmm3							\n\t	vmovapd	%%zmm16,%%zmm15	\n\t"\
		"vpermt2pd				%%zmm5,		%%zmm30,%%zmm3 	\n\t	vpermt2pd			%%zmm17,	%%zmm30,%%zmm15	\n\t"\
		"vpermt2pd				%%zmm5,		%%zmm31,%%zmm4	\n\t	vpermt2pd			%%zmm17,	%%zmm31,%%zmm16	\n\t"\
		/* Write transposed data to the local store: */\
		"vmovaps		%%zmm10,0x000(%%rsi)		\n\t	vmovaps	%%zmm22,0x040(%%rsi)	\n\t"/* Out 0 */\
		"vmovaps		%%zmm11,0x800(%%rsi)		\n\t	vmovaps	%%zmm23,0x840(%%rsi)	\n\t"/* Out 4 */\
		"vmovaps		%%zmm2 ,0x400(%%rsi)		\n\t	vmovaps	%%zmm14,0x440(%%rsi)	\n\t"/* Out 2 */\
		"vmovaps		%%zmm8 ,0xc00(%%rsi)		\n\t	vmovaps	%%zmm20,0xc40(%%rsi)	\n\t"/* Out 6 */\
		"vmovaps		%%zmm0 ,0x200(%%rsi)		\n\t	vmovaps	%%zmm12,0x240(%%rsi)	\n\t"/* Out 1 */\
		"vmovaps		%%zmm1 ,0xa00(%%rsi)		\n\t	vmovaps	%%zmm13,0xa40(%%rsi)	\n\t"/* Out 5 */\
		"vmovaps		%%zmm3 ,0x600(%%rsi)		\n\t	vmovaps	%%zmm15,0x640(%%rsi)	\n\t"/* Out 3 */\
		"vmovaps		%%zmm4, 0xe00(%%rsi)		\n\t	vmovaps	%%zmm16,0xe40(%%rsi)	\n\t"/* Out 7 */\
		"\n\t"\
		/************************************************************************/							/************************************************************************************************************/\
		/* Forward DIF radix-32 pass on the interleaved block1 and block2 data: */							/* Next 2 blocks operate on odd-indexed elements from the unpck*pd commands which we stored to temporaries: */\
		/************************************************************************/							/************************************************************************************************************/\
	/*...Block 0: */																						/*...Block 3: */\
		"movq	%[__isrt2],%%rsi				\n\t	leaq	0x11c0(%%rsi),%%rdi	/* two */	\n\t"\
	/* SSE2_RADIX4_DIF_4WRAPPER(c00,c08,c10,c18,r00)...DIF_4WRAPPER_2NDOFTWO(c04,c0C,c14,c1C,r08)	SSE2_RADIX4_DIF_4TWIDDLE(r20,r24,r22,r26,r20,c01)...DIF_4TWIDDLE_2NDOFTWO(r28,r2C,r2A,r2E,r28,c05) */\
		"movq		%[__r00],%%rcx				\n\t	/*addq		$0x200,%%rcx // __r08 */	\n\t	movq %[__r20],%%r10 \n\t movq %[__c01],%%r11 \n\t leaq 0x080(%%r10),%%r12 \n\t	movq %%r10,%%r13 \n\t"\
		"movq		%[__c00],%%rdx				\n\t	/*addq		$0x200,%%rdx // __c04 */	\n\t	vmovaps		 (%%r10),%%zmm18 			\n\t	vmovaps		0x200(%%r10),%%zmm26		\n\t"\
		"vmovaps			 (%%rcx),%%zmm0		\n\t	vmovaps		0x200(%%rcx),%%zmm8			\n\t	vmovaps		 (%%r12),%%zmm20			\n\t	vmovaps		0x200(%%r12),%%zmm28		\n\t"\
		"vmovaps		0x040(%%rcx),%%zmm1		\n\t	vmovaps		0x240(%%rcx),%%zmm9			\n\t	vmovaps	0x040(%%r10),%%zmm19			\n\t	vmovaps		0x240(%%r10),%%zmm27		\n\t"\
		"vmovaps		%%zmm0,%%zmm2			\n\t	vmovaps		%%zmm8,%%zmm10				\n\t	vmovaps	0x040(%%r12),%%zmm21			\n\t	vmovaps		0x240(%%r12),%%zmm29		\n\t"\
		"vmovaps		%%zmm1,%%zmm3			\n\t	vmovaps		%%zmm9,%%zmm11				\n\t	vmovaps		 (%%r11),%%zmm22			\n\t	vmovaps		0x200(%%r11),%%zmm30		\n\t"\
		"vmulpd			 (%%rdx),%%zmm0,%%zmm0	\n\t	vmulpd		0x200(%%rdx),%%zmm8 ,%%zmm8	\n\t	vmovaps	0x040(%%r11),%%zmm23			\n\t	vmovaps		0x240(%%r11),%%zmm31		\n\t"\
		"vmulpd			 (%%rdx),%%zmm1,%%zmm1	\n\t	vmulpd		0x200(%%rdx),%%zmm9 ,%%zmm9	\n\t"\
/********** Where do the zmm6,7 come from here? ***********/\
	"vfnmadd231pd	0x040(%%rdx),%%zmm7,%%zmm0	\n\t vfnmadd231pd	0x240(%%rdx),%%zmm11,%%zmm8 \n\t"\
	" vfmadd231pd	0x040(%%rdx),%%zmm6,%%zmm1	\n\t vfmadd231pd	0x240(%%rdx),%%zmm10,%%zmm9 \n\t	vmulpd		%%zmm22,%%zmm18,%%zmm16		\n\t	vmulpd		%%zmm30,%%zmm26,%%zmm24	\n\t"\
		"vmovaps		%%zmm0,%%zmm2			\n\t	vmovaps		%%zmm8,%%zmm10				\n\t	vmulpd		%%zmm22,%%zmm19,%%zmm17		\n\t	vmulpd		%%zmm30,%%zmm27,%%zmm25	\n\t"\
		"vmovaps		%%zmm1,%%zmm3			\n\t	vmovaps		%%zmm9,%%zmm11				\n\t vfnmadd231pd	%%zmm23,%%zmm19,%%zmm16		\n\t vfnmadd231pd	%%zmm31,%%zmm27,%%zmm24		\n\t"\
		"addq		$0x080,%%rdx																\n\t vfmadd231pd	%%zmm23,%%zmm18,%%zmm17		\n\t vfmadd231pd	%%zmm31,%%zmm26,%%zmm25		\n\t"\
		"vmovaps		0x080(%%rcx),%%zmm6		\n\t	vmovaps		0x280(%%rcx),%%zmm14		\n\t	vmovaps		%%zmm20,%%zmm22				\n\t	vmovaps		%%zmm28,%%zmm30		\n\t"\
		"vmovaps		0x0c0(%%rcx),%%zmm7		\n\t	vmovaps		0x2c0(%%rcx),%%zmm15		\n\t	vmovaps		%%zmm21,%%zmm23				\n\t	vmovaps		%%zmm29,%%zmm31		\n\t"\
		"																								vmulpd	  0x080(%%r11),%%zmm20,%%zmm20	\n\t	vmulpd		0x280(%%r11),%%zmm28,%%zmm28\n\t"\
		"																								vmulpd	  0x080(%%r11),%%zmm21,%%zmm21	\n\t	vmulpd		0x280(%%r11),%%zmm29,%%zmm29\n\t"\
		"vmulpd			 (%%rdx),%%zmm6,%%zmm4	\n\t	vmulpd		0x200(%%rdx),%%zmm14,%%zmm12\n\t vfnmadd231pd 0x0c0(%%r11),%%zmm23,%%zmm20	\n\t vfnmadd231pd	0x2c0(%%r11),%%zmm31,%%zmm28\n\t"\
		"vmulpd			 (%%rdx),%%zmm7,%%zmm5	\n\t	vmulpd		0x200(%%rdx),%%zmm15,%%zmm13\n\t vfmadd231pd  0x0c0(%%r11),%%zmm22,%%zmm21	\n\t vfmadd231pd	0x2c0(%%r11),%%zmm30,%%zmm29\n\t"\
	"vfnmadd231pd	0x040(%%rdx),%%zmm7,%%zmm4	\n\t vfnmadd231pd	0x240(%%rdx),%%zmm15,%%zmm12\n\t	vmovaps		%%zmm16,%%zmm18				\n\t	vmovaps		%%zmm24,%%zmm26		\n\t"\
	" vfmadd231pd	0x040(%%rdx),%%zmm6,%%zmm5	\n\t vfmadd231pd	0x240(%%rdx),%%zmm14,%%zmm13\n\t	vmovaps		%%zmm17,%%zmm19				\n\t	vmovaps		%%zmm25,%%zmm27		\n\t"\
		"vaddpd		%%zmm4,%%zmm0,%%zmm0		\n\t	vaddpd		%%zmm12,%%zmm8 ,%%zmm8		\n\t	addq		$0x100,%%r12				\n\t	addq		$0x180,%%r11			\n\t"\
		"vaddpd		%%zmm5,%%zmm1,%%zmm1		\n\t	vaddpd		%%zmm13,%%zmm9 ,%%zmm9		\n\t	vmovaps			 (%%r12),%%zmm22		\n\t	vmovaps		0x200(%%r12),%%zmm30		\n\t"\
		"vsubpd		%%zmm4,%%zmm2,%%zmm2		\n\t	vsubpd		%%zmm12,%%zmm10,%%zmm10		\n\t	vmovaps		0x040(%%r12),%%zmm23		\n\t	vmovaps		0x240(%%r12),%%zmm31		\n\t"\
		"vsubpd		%%zmm5,%%zmm3,%%zmm3		\n\t	vsubpd		%%zmm13,%%zmm11,%%zmm11		\n\t	vaddpd		%%zmm20,%%zmm16,%%zmm16		\n\t	vaddpd		%%zmm28,%%zmm24,%%zmm24		\n\t"\
		"addq		$0x100,%%rdx																\n\t	vaddpd		%%zmm21,%%zmm17,%%zmm17		\n\t	vaddpd		%%zmm29,%%zmm25,%%zmm25		\n\t"\
		"vmovaps		0x180(%%rcx),%%zmm6		\n\t	vmovaps		0x380(%%rcx),%%zmm14		\n\t	vsubpd		%%zmm20,%%zmm18,%%zmm18		\n\t	vsubpd		%%zmm28,%%zmm26,%%zmm26		\n\t"\
		"vmovaps		0x1c0(%%rcx),%%zmm7		\n\t	vmovaps		0x3c0(%%rcx),%%zmm15		\n\t	vsubpd		%%zmm21,%%zmm19,%%zmm19		\n\t	vsubpd		%%zmm29,%%zmm27,%%zmm27		\n\t"\
		"vmulpd			 (%%rdx),%%zmm6,%%zmm4	\n\t	vmulpd		0x200(%%rdx),%%zmm14,%%zmm12\n\t	vmulpd		   (%%r11),%%zmm22,%%zmm20	\n\t	vmulpd		0x200(%%r11),%%zmm30,%%zmm28\n\t"\
		"vmulpd			 (%%rdx),%%zmm7,%%zmm5	\n\t	vmulpd		0x200(%%rdx),%%zmm15,%%zmm13\n\t	vmulpd		   (%%r11),%%zmm23,%%zmm21	\n\t	vmulpd		0x200(%%r11),%%zmm31,%%zmm29\n\t"\
	"vfnmadd231pd	0x040(%%rdx),%%zmm7,%%zmm4	\n\t vfnmadd231pd	0x240(%%rdx),%%zmm15,%%zmm12\n\t vfnmadd231pd 0x040(%%r11),%%zmm23,%%zmm20	\n\t vfnmadd231pd	0x240(%%r11),%%zmm31,%%zmm28\n\t"\
	" vfmadd231pd	0x040(%%rdx),%%zmm6,%%zmm5	\n\t vfmadd231pd	0x240(%%rdx),%%zmm14,%%zmm13\n\t vfmadd231pd  0x040(%%r11),%%zmm22,%%zmm21	\n\t vfmadd231pd	0x240(%%r11),%%zmm30,%%zmm29\n\t"\
		"vmovaps		%%zmm5,0x040(%%rcx)		\n\t	vmovaps		%%zmm13,0x240(%%rcx)		\n\t	vmovaps		%%zmm21,0x040(%%r13)		\n\t	vmovaps		%%zmm29,0x240(%%r13)		\n\t"\
		"vmovaps		%%zmm4,     (%%rcx)		\n\t	vmovaps		%%zmm12,0x200(%%rcx)		\n\t	vmovaps		%%zmm20,     (%%r13)		\n\t	vmovaps		%%zmm28,0x200(%%r13)		\n\t"\
		"subq		$0x080,%%rdx																\n\t	addq	$0x100,%%r10					\n\t	subq	$0x080,%%r11					\n\t"\
		"vmovaps		0x100(%%rcx),%%zmm6		\n\t	vmovaps		0x300(%%rcx),%%zmm14		\n\t	vmovaps			 (%%r10),%%zmm22		\n\t	vmovaps		0x200(%%r10),%%zmm30		\n\t"\
		"vmovaps		0x140(%%rcx),%%zmm7		\n\t	vmovaps		0x340(%%rcx),%%zmm15		\n\t	vmovaps		0x040(%%r10),%%zmm23		\n\t	vmovaps		0x240(%%r10),%%zmm31		\n\t"\
		"vmulpd			 (%%rdx),%%zmm6,%%zmm4	\n\t	vmulpd		0x200(%%rdx),%%zmm14,%%zmm12\n\t	vmulpd		   (%%r11),%%zmm22,%%zmm20	\n\t	vmulpd		0x200(%%r11),%%zmm30,%%zmm28\n\t"\
		"vmulpd			 (%%rdx),%%zmm7,%%zmm5	\n\t	vmulpd		0x200(%%rdx),%%zmm15,%%zmm13\n\t	vmulpd		   (%%r11),%%zmm23,%%zmm21	\n\t	vmulpd		0x200(%%r11),%%zmm31,%%zmm29\n\t"\
	"vfnmadd231pd	0x040(%%rdx),%%zmm7,%%zmm4	\n\t vfnmadd231pd	0x240(%%rdx),%%zmm15,%%zmm12\n\t vfnmadd231pd 0x040(%%r11),%%zmm23,%%zmm20	\n\t vfnmadd231pd	0x240(%%r11),%%zmm31,%%zmm28\n\t"\
	" vfmadd231pd	0x040(%%rdx),%%zmm6,%%zmm5	\n\t vfmadd231pd	0x240(%%rdx),%%zmm14,%%zmm13\n\t vfmadd231pd  0x040(%%r11),%%zmm22,%%zmm21	\n\t vfmadd231pd	0x240(%%r11),%%zmm30,%%zmm29\n\t"\
		"vmovaps		%%zmm5,%%zmm7			\n\t	vmovaps		%%zmm13,%%zmm15				\n\t	vmovaps		%%zmm21,%%zmm23				\n\t	vmovaps		%%zmm29,%%zmm31		\n\t"\
		"vmovaps		%%zmm4,%%zmm6			\n\t	vmovaps		%%zmm12,%%zmm14				\n\t	vmovaps		%%zmm20,%%zmm22				\n\t	vmovaps		%%zmm28,%%zmm30		\n\t"\
		"vsubpd			 (%%rcx),%%zmm4,%%zmm4	\n\t	vsubpd		0x200(%%rcx),%%zmm12,%%zmm12\n\t	vsubpd			 (%%r13),%%zmm20,%%zmm20\n\t	vsubpd		0x200(%%r13),%%zmm28,%%zmm28\n\t"\
		"vsubpd		0x040(%%rcx),%%zmm5,%%zmm5	\n\t	vsubpd		0x240(%%rcx),%%zmm13,%%zmm13\n\t	vsubpd		0x040(%%r13),%%zmm21,%%zmm21\n\t	vsubpd		0x240(%%r13),%%zmm29,%%zmm29\n\t"\
		"vaddpd			 (%%rcx),%%zmm6,%%zmm6	\n\t	vaddpd		0x200(%%rcx),%%zmm14,%%zmm14\n\t	vaddpd			 (%%r13),%%zmm22,%%zmm22\n\t	vaddpd		0x200(%%r13),%%zmm30,%%zmm30\n\t"\
		"vaddpd		0x040(%%rcx),%%zmm7,%%zmm7	\n\t	vaddpd		0x240(%%rcx),%%zmm15,%%zmm15\n\t	vaddpd		0x040(%%r13),%%zmm23,%%zmm23\n\t	vaddpd		0x240(%%r13),%%zmm31,%%zmm31\n\t"\
		"vsubpd		%%zmm6,%%zmm0,%%zmm0		\n\t	vsubpd		%%zmm14,%%zmm8,%%zmm8		\n\t	vsubpd		%%zmm22,%%zmm16,%%zmm16		\n\t	vsubpd		%%zmm30,%%zmm24,%%zmm24\n\t"\
		"vsubpd		%%zmm7,%%zmm1,%%zmm1		\n\t	vsubpd		%%zmm15,%%zmm9,%%zmm9		\n\t	vsubpd		%%zmm23,%%zmm17,%%zmm17		\n\t	vsubpd		%%zmm31,%%zmm25,%%zmm25\n\t"\
		"vsubpd		%%zmm5,%%zmm2,%%zmm2		\n\t	vsubpd		%%zmm13,%%zmm10,%%zmm10		\n\t	vsubpd		%%zmm21,%%zmm18,%%zmm18		\n\t	vsubpd		%%zmm29,%%zmm26,%%zmm26\n\t"\
		"vsubpd		%%zmm4,%%zmm3,%%zmm3		\n\t	vsubpd		%%zmm12,%%zmm11,%%zmm11		\n\t	vsubpd		%%zmm20,%%zmm19,%%zmm19		\n\t	vsubpd		%%zmm28,%%zmm27,%%zmm27\n\t"\
	"vmovaps %%zmm28,(%%r13) \n\t"/* spill zmm28 to make room for two */"vmovaps (%%rdi),%%zmm28\n\t"/* two */\
	"vfmadd132pd	%%zmm28,%%zmm0,%%zmm6		\n\t	vfmadd132pd	%%zmm28,%%zmm8 ,%%zmm14		\n\t	vfmadd132pd	%%zmm28,%%zmm16,%%zmm22		\n\t	vfmadd132pd	%%zmm28,%%zmm24,%%zmm30			\n\t"\
	"vfmadd132pd	%%zmm28,%%zmm1,%%zmm7		\n\t	vfmadd132pd	%%zmm28,%%zmm9 ,%%zmm15		\n\t	vfmadd132pd	%%zmm28,%%zmm17,%%zmm23		\n\t	vfmadd132pd	%%zmm28,%%zmm25,%%zmm31			\n\t"\
	"vfmadd132pd	%%zmm28,%%zmm2,%%zmm5		\n\t	vfmadd132pd	%%zmm28,%%zmm10,%%zmm13		\n\t	vfmadd132pd	%%zmm28,%%zmm18,%%zmm21		\n\t	vfmadd132pd	%%zmm28,%%zmm26,%%zmm29			\n\t"\
	"vfmadd132pd	%%zmm28,%%zmm3,%%zmm4		\n\t	vfmadd132pd	%%zmm28,%%zmm11,%%zmm12		\n\t	vfmadd132pd	%%zmm28,%%zmm19,%%zmm20		\n\t	vfmadd132pd	(%%r13),%%zmm27,%%zmm28			\n\t"\
		"												vmovaps		%%zmm14,0x200(%%rcx)		\n\t														vmovaps		%%zmm30,0x200(%%r13)			\n\t"\
		"												vmovaps		%%zmm15,0x240(%%rcx)		\n\t														vmovaps		%%zmm31,0x240(%%r13)			\n\t"\
		"												vmovaps		%%zmm10,%%zmm14				\n\t														vmovaps		%%zmm26,%%zmm30					\n\t"\
		"												vmovaps		%%zmm13,%%zmm15				\n\t														vmovaps		%%zmm29,%%zmm31					\n\t"\
		"												vsubpd		%%zmm12,%%zmm10,%%zmm10		\n\t														vsubpd		%%zmm28,%%zmm26,%%zmm26	\n\t"\
		"												vsubpd		%%zmm11,%%zmm13,%%zmm13		\n\t														vsubpd		%%zmm27,%%zmm29,%%zmm29	\n\t"\
		"												vaddpd		%%zmm12,%%zmm14,%%zmm14		\n\t														vaddpd		%%zmm28,%%zmm30,%%zmm30	\n\t"\
		"												vaddpd		%%zmm11,%%zmm15,%%zmm15		\n\t														vaddpd		%%zmm27,%%zmm31,%%zmm31	\n\t"\
		/*****	SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r00) *****/												/*****	SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r20) *****/\
		"												vmovaps		0x200(%%rcx),%%zmm11			\n\t														vmovaps		0x200(%%r13),%%zmm27			\n\t"\
		"												vmovaps		0x240(%%rcx),%%zmm12			\n\t														vmovaps		0x240(%%r13),%%zmm28			\n\t"\
	"vmovaps %%zmm13,(%%rcx) \n\t"/* spill zmm13 to make room for isrt2 */"vmovaps (%%rsi),%%zmm13 	\n\t	vmovaps %%zmm29,(%%r13) \n\t"/* spill zmm29 to make room for isrt2 */"vmovaps (%%rsi),%%zmm29	\n\t"/*isrt2*/\
		"vsubpd		%%zmm11,%%zmm6,%%zmm6		\n\t	vfnmadd231pd	%%zmm13,%%zmm10,%%zmm2		\n\t		vsubpd		%%zmm27,%%zmm22,%%zmm22	\n\t	vfnmadd231pd	%%zmm29,%%zmm26,%%zmm18		\n\t"\
		"vsubpd		%%zmm9 ,%%zmm0,%%zmm0		\n\t	vfnmadd231pd	%%zmm13,%%zmm15,%%zmm5		\n\t		vsubpd		%%zmm25,%%zmm16,%%zmm16	\n\t	vfnmadd231pd	%%zmm29,%%zmm31,%%zmm21		\n\t"\
		"vsubpd		%%zmm12,%%zmm7,%%zmm7		\n\t	vfnmadd231pd	%%zmm13,%%zmm14,%%zmm4		\n\t		vsubpd		%%zmm28,%%zmm23,%%zmm23	\n\t	vfnmadd231pd	%%zmm29,%%zmm30,%%zmm20		\n\t"\
		"vsubpd		%%zmm8 ,%%zmm1,%%zmm1		\n\t	vfnmadd231pd	(%%rcx),%%zmm13,%%zmm3		\n\t		vsubpd		%%zmm24,%%zmm17,%%zmm17	\n\t	vfnmadd231pd	(%%r13),%%zmm29,%%zmm19		\n\t"\
		"vmovaps		%%zmm6,0x200(%%rcx)		\n\t		vmovaps		%%zmm2,0x280(%%rcx)	\n\t		vmovaps		%%zmm22,0x200(%%r13)	\n\t		vmovaps		%%zmm18,0x280(%%r13)	\n\t"\
		"vmovaps		%%zmm0,0x100(%%rcx)		\n\t		vmovaps		%%zmm5,0x180(%%rcx)	\n\t		vmovaps		%%zmm16,0x100(%%r13)	\n\t		vmovaps		%%zmm21,0x180(%%r13)	\n\t"\
		"vmovaps		%%zmm7,0x240(%%rcx)		\n\t		vmovaps		%%zmm4,0x2c0(%%rcx)	\n\t		vmovaps		%%zmm23,0x240(%%r13)	\n\t		vmovaps		%%zmm20,0x2c0(%%r13)	\n\t"\
		"vmovaps		%%zmm1,0x340(%%rcx)		\n\t		vmovaps		%%zmm3,0x3c0(%%rcx)	\n\t		vmovaps		%%zmm17,0x340(%%r13)	\n\t		vmovaps		%%zmm19,0x3c0(%%r13)	\n\t"\
	"vmovaps	0x80(%%rdi),%%zmm13																	\n\t	vmovaps	0x80(%%rdi),%%zmm29 \n\t"/* sqrt2 *//* rcol: *2 => *sqrt2 due to above FMA leaving zmm26,13-15 unmultiplied-by-isrt2: */\
	"vfmadd132pd	(%%rdi),%%zmm6,%%zmm11		\n\t	 vfmadd132pd	%%zmm13,%%zmm2,%%zmm10	\n\t	vfmadd132pd	(%%rdi),%%zmm22,%%zmm27	\n\t	 vfmadd132pd	%%zmm29,%%zmm18,%%zmm26	\n\t"\
	"vfmadd132pd	(%%rdi),%%zmm0,%%zmm9		\n\t	 vfmadd132pd	%%zmm13,%%zmm5,%%zmm15	\n\t	vfmadd132pd	(%%rdi),%%zmm16,%%zmm25	\n\t	 vfmadd132pd	%%zmm29,%%zmm21,%%zmm31	\n\t"\
	"vfmadd132pd	(%%rdi),%%zmm7,%%zmm12		\n\t	 vfmadd132pd	%%zmm13,%%zmm4,%%zmm14	\n\t	vfmadd132pd	(%%rdi),%%zmm23,%%zmm28	\n\t	 vfmadd132pd	%%zmm29,%%zmm20,%%zmm30	\n\t"\
	"vfmadd132pd	(%%rdi),%%zmm1,%%zmm8		\n\t	 vfmadd132pd	(%%rcx),%%zmm3,%%zmm13	\n\t	vfmadd132pd	(%%rdi),%%zmm17,%%zmm24	\n\t	 vfmadd132pd	(%%r13),%%zmm19,%%zmm29	\n\t"\
		"vmovaps		%%zmm11,     (%%rcx)	\n\t		vmovaps		%%zmm10,0x080(%%rcx)	\n\t		vmovaps		%%zmm27,     (%%r13)	\n\t		vmovaps		%%zmm26,0x080(%%r13)	\n\t"\
		"vmovaps		%%zmm9 ,0x300(%%rcx)	\n\t		vmovaps		%%zmm15,0x380(%%rcx)	\n\t		vmovaps		%%zmm25,0x300(%%r13)	\n\t		vmovaps		%%zmm31,0x380(%%r13)	\n\t"\
		"vmovaps		%%zmm12,0x040(%%rcx)	\n\t		vmovaps		%%zmm14,0x0c0(%%rcx)	\n\t		vmovaps		%%zmm28,0x040(%%r13)	\n\t		vmovaps		%%zmm30,0x0c0(%%r13)	\n\t"\
		"vmovaps		%%zmm8 ,0x140(%%rcx)	\n\t		vmovaps		%%zmm13,0x1c0(%%rcx)	\n\t		vmovaps		%%zmm24,0x140(%%r13)	\n\t		vmovaps		%%zmm29,0x1c0(%%r13)	\n\t"\
	/*...Block 2: */	/*...Block 4: */\
	/* SSE2_RADIX4_DIF_4WRAPPER(c02,c0A,c12,c1A,r10)	...DIF_4WRAPPER_2NDOFTWO(c06,c0E,c16,c1E,r18) */	/*	SSE2_RADIX4_DIF_4TWIDDLE(r30,r34,r32,r36,r30,c03)	...DIF_4TWIDDLE_2NDOFTWO(r38,r3C,r3A,r3E,r38,c07) */\
		"movq		%[__r10],%%rcx				\n\t		/*addq		$0x200,%%rcx // __r18 */	\n\t		movq %[__r30],%%r10 \n\t movq %[__c03],%%r11 \n\t movq %%r10,%%r12 \n\t addq $0x080,%%r12	\n\t"\
		"movq		%[__c02],%%rdx				\n\t		/*addq		$0x200,%%rdx // __c06 */	\n\t		vmovaps		 (%%r10),%%zmm16 \n\t movq %%r10,%%r13 \n\t vmovaps 0x200(%%r10),%%zmm24	\n\t"\
		"vmovaps			 (%%rcx),%%zmm0		\n\t		vmovaps		0x200(%%rcx),%%zmm8			\n\t		vmovaps		 (%%r12),%%zmm20		\n\t		vmovaps		0x200(%%r12),%%zmm28		\n\t"\
		"vmovaps		0x040(%%rcx),%%zmm1		\n\t		vmovaps		0x240(%%rcx),%%zmm9			\n\t		vmovaps	0x040(%%r10),%%zmm17		\n\t		vmovaps		0x240(%%r10),%%zmm25		\n\t"\
		"vmovaps		%%zmm0,%%zmm2			\n\t		vmovaps		%%zmm8,%%zmm10				\n\t		vmovaps	0x040(%%r12),%%zmm21		\n\t		vmovaps		0x240(%%r12),%%zmm29		\n\t"\
		"vmovaps		%%zmm1,%%zmm3			\n\t		vmovaps		%%zmm9,%%zmm11				\n\t		vmovaps		 (%%r11),%%zmm22		\n\t		vmovaps		0x200(%%r11),%%zmm30		\n\t"\
		"vmulpd			 (%%rdx),%%zmm0,%%zmm0	\n\t		vmulpd		0x200(%%rdx),%%zmm8,%%zmm8	\n\t		vmovaps	0x040(%%r11),%%zmm23		\n\t		vmovaps		0x240(%%r11),%%zmm31		\n\t"\
		"vmulpd			 (%%rdx),%%zmm1,%%zmm1	\n\t		vmulpd		0x200(%%rdx),%%zmm9,%%zmm9	\n\t		vmovaps	%%zmm16,%%zmm18		\n\t		vmovaps		%%zmm24,%%zmm26		\n\t"\
	"vfnmadd231pd	0x040(%%rdx),%%zmm3,%%zmm0	\n\t	vfnmadd231pd	0x240(%%rdx),%%zmm11,%%zmm8 \n\t		vmovaps	%%zmm17,%%zmm19		\n\t		vmovaps		%%zmm25,%%zmm27		\n\t"\
	" vfmadd231pd	0x040(%%rdx),%%zmm2,%%zmm1	\n\t	 vfmadd231pd	0x240(%%rdx),%%zmm10,%%zmm9 \n\t		vmulpd		%%zmm22,%%zmm16,%%zmm16\n\t		vmulpd		%%zmm30,%%zmm24,%%zmm24\n\t"\
		"vmovaps		%%zmm0,%%zmm2			\n\t		vmovaps		%%zmm8,%%zmm10				\n\t		vmulpd		%%zmm22,%%zmm17,%%zmm17\n\t		vmulpd		%%zmm30,%%zmm25,%%zmm25\n\t"\
		"vmovaps		%%zmm1,%%zmm3			\n\t		vmovaps		%%zmm9,%%zmm11				\n\t	vfnmadd231pd	%%zmm23,%%zmm19,%%zmm16		\n\t	vfnmadd231pd	%%zmm31,%%zmm27,%%zmm24		\n\t"\
		"addq		$0x080,%%rdx																\n\t	 vfmadd231pd	%%zmm23,%%zmm18,%%zmm17		\n\t	 vfmadd231pd	%%zmm31,%%zmm26,%%zmm25		\n\t"\
		"vmovaps		0x080(%%rcx),%%zmm4		\n\t		vmovaps		0x280(%%rcx),%%zmm12		\n\t		vmovaps		%%zmm20,%%zmm22		\n\t		vmovaps		%%zmm28,%%zmm30		\n\t"\
		"vmovaps		0x0c0(%%rcx),%%zmm5		\n\t		vmovaps		0x2c0(%%rcx),%%zmm13		\n\t		vmovaps		%%zmm21,%%zmm23		\n\t		vmovaps		%%zmm29,%%zmm31		\n\t"\
		"vmovaps		%%zmm4,%%zmm6			\n\t		vmovaps		%%zmm12,%%zmm14				\n\t		vmulpd		0x080(%%r11),%%zmm20,%%zmm20\n\t		vmulpd		0x280(%%r11),%%zmm28,%%zmm28\n\t"\
		"vmovaps		%%zmm5,%%zmm7			\n\t		vmovaps		%%zmm13,%%zmm15				\n\t		vmulpd		0x080(%%r11),%%zmm21,%%zmm21\n\t		vmulpd		0x280(%%r11),%%zmm29,%%zmm29\n\t"\
		"vmulpd			 (%%rdx),%%zmm4,%%zmm4	\n\t		vmulpd		0x200(%%rdx),%%zmm12,%%zmm12\n\t		vmovaps		%%zmm16,%%zmm18		\n\t		vmovaps		%%zmm24,%%zmm26		\n\t"\
		"vmulpd			 (%%rdx),%%zmm5,%%zmm5	\n\t		vmulpd		0x200(%%rdx),%%zmm13,%%zmm13\n\t		vmovaps		%%zmm17,%%zmm19		\n\t		vmovaps		%%zmm25,%%zmm27		\n\t"\
	"vfnmadd231pd	0x040(%%rdx),%%zmm7,%%zmm4	\n\t	vfnmadd231pd	0x240(%%rdx),%%zmm15,%%zmm12\n\t	vfnmadd231pd	0x0c0(%%r11),%%zmm23,%%zmm20\n\t	vfnmadd231pd	0x2c0(%%r11),%%zmm31,%%zmm28\n\t"\
	" vfmadd231pd	0x040(%%rdx),%%zmm6,%%zmm5	\n\t	 vfmadd231pd	0x240(%%rdx),%%zmm14,%%zmm13\n\t	 vfmadd231pd	0x0c0(%%r11),%%zmm22,%%zmm21\n\t	 vfmadd231pd	0x2c0(%%r11),%%zmm30,%%zmm29\n\t"\
		"vaddpd		%%zmm4,%%zmm0,%%zmm0		\n\t		vaddpd		%%zmm12,%%zmm8,%%zmm8	\n\t		addq		$0x100,%%r12			\n\t		addq		$0x180,%%r11			\n\t"\
		"vaddpd		%%zmm5,%%zmm1,%%zmm1		\n\t		vaddpd		%%zmm13,%%zmm9,%%zmm9	\n\t		vmovaps			 (%%r12),%%zmm22	\n\t		vaddpd		%%zmm28,%%zmm24,%%zmm24\n\t"\
		"vsubpd		%%zmm4,%%zmm2,%%zmm2		\n\t		vsubpd		%%zmm12,%%zmm10,%%zmm10\n\t		vmovaps		0x040(%%r12),%%zmm23	\n\t		vaddpd		%%zmm29,%%zmm25,%%zmm25\n\t"\
		"vsubpd		%%zmm5,%%zmm3,%%zmm3		\n\t		vsubpd		%%zmm13,%%zmm11,%%zmm11\n\t		vaddpd		%%zmm20,%%zmm16,%%zmm16\n\t		vsubpd		%%zmm28,%%zmm26,%%zmm26\n\t"\
		"addq		$0x100,%%rdx																\n\t		vaddpd		%%zmm21,%%zmm17,%%zmm17\n\t		vsubpd		%%zmm29,%%zmm27,%%zmm27\n\t"\
		"vmovaps		0x180(%%rcx),%%zmm4		\n\t		vmovaps		0x380(%%rcx),%%zmm12		\n\t		vsubpd		%%zmm20,%%zmm18,%%zmm18\n\t		vmovaps		0x200(%%r12),%%zmm28		\n\t"\
		"vmovaps		0x1c0(%%rcx),%%zmm5		\n\t		vmovaps		0x3c0(%%rcx),%%zmm13		\n\t		vsubpd		%%zmm21,%%zmm19,%%zmm19\n\t		vmovaps		0x240(%%r12),%%zmm29		\n\t"\
		"vmovaps		%%zmm4,%%zmm6			\n\t		vmovaps		%%zmm12,%%zmm14				\n\t		vmovaps		%%zmm22,%%zmm20		\n\t		vmovaps		0x200(%%r12),%%zmm30		\n\t"\
		"vmovaps		%%zmm5,%%zmm7			\n\t		vmovaps		%%zmm13,%%zmm15				\n\t		vmovaps		%%zmm23,%%zmm21		\n\t		vmovaps		0x240(%%r12),%%zmm31		\n\t"\
		"vmulpd			 (%%rdx),%%zmm4,%%zmm4	\n\t		vmulpd		0x200(%%rdx),%%zmm12,%%zmm12\n\t		vmulpd			 (%%r11),%%zmm20,%%zmm20\n\t		vmulpd		0x200(%%r11),%%zmm28,%%zmm28\n\t"\
		"vmulpd			 (%%rdx),%%zmm5,%%zmm5	\n\t		vmulpd		0x200(%%rdx),%%zmm13,%%zmm13\n\t		vmulpd			 (%%r11),%%zmm21,%%zmm21\n\t		vmulpd		0x200(%%r11),%%zmm29,%%zmm29\n\t"\
	"vfnmadd231pd	0x040(%%rdx),%%zmm7,%%zmm4	\n\t	vfnmadd231pd	0x240(%%rdx),%%zmm15,%%zmm12\n\t	vfnmadd231pd	0x040(%%r11),%%zmm23,%%zmm20\n\t	vfnmadd231pd	0x240(%%r11),%%zmm31,%%zmm28\n\t"\
	" vfmadd231pd	0x040(%%rdx),%%zmm6,%%zmm5	\n\t	 vfmadd231pd	0x240(%%rdx),%%zmm14,%%zmm13\n\t	 vfmadd231pd	0x040(%%r11),%%zmm22,%%zmm21\n\t	 vfmadd231pd	0x240(%%r11),%%zmm30,%%zmm29\n\t"\
		"vmovaps		%%zmm5,0x040(%%rcx)		\n\t		vmovaps		%%zmm13,0x240(%%rcx)		\n\t		vmovaps		%%zmm21,0x040(%%r13)	\n\t		vmovaps		%%zmm29,0x240(%%r13)		\n\t"\
		"vmovaps		%%zmm4,     (%%rcx)		\n\t		vmovaps		%%zmm12,0x200(%%rcx)		\n\t		vmovaps		%%zmm20,     (%%r13)	\n\t		vmovaps		%%zmm28,0x200(%%r13)		\n\t"\
		"subq		$0x080,%%rdx																\n\t		addq	$0x100,%%r10					\n\t		subq	$0x080,%%r11					\n\t"\
		"vmovaps		0x100(%%rcx),%%zmm4		\n\t		vmovaps		0x300(%%rcx),%%zmm12		\n\t		vmovaps			 (%%r10),%%zmm20	\n\t		vmovaps		0x200(%%r10),%%zmm28		\n\t"\
		"vmovaps		0x140(%%rcx),%%zmm5		\n\t		vmovaps		0x340(%%rcx),%%zmm13		\n\t		vmovaps		0x040(%%r10),%%zmm21	\n\t		vmovaps		0x240(%%r10),%%zmm29		\n\t"\
		"vmovaps		%%zmm4,%%zmm6			\n\t		vmovaps		%%zmm12,%%zmm14				\n\t		vmovaps		%%zmm20,%%zmm22		\n\t		vmovaps		0x200(%%r10),%%zmm30		\n\t"\
		"vmovaps		%%zmm5,%%zmm7			\n\t		vmovaps		%%zmm13,%%zmm15				\n\t		vmovaps		%%zmm21,%%zmm23		\n\t		vmovaps		0x240(%%r10),%%zmm31		\n\t"\
		"vmulpd			 (%%rdx),%%zmm4,%%zmm4	\n\t		vmulpd		0x200(%%rdx),%%zmm12,%%zmm12\n\t		vmulpd			 (%%r11),%%zmm20,%%zmm20\n\t		vmulpd		0x200(%%r11),%%zmm28,%%zmm28\n\t"\
		"vmulpd			 (%%rdx),%%zmm5,%%zmm5	\n\t		vmulpd		0x200(%%rdx),%%zmm13,%%zmm13\n\t		vmulpd			 (%%r11),%%zmm21,%%zmm21\n\t		vmulpd		0x200(%%r11),%%zmm29,%%zmm29\n\t"\
	"vfnmadd231pd	0x040(%%rdx),%%zmm7,%%zmm4	\n\t	vfnmadd231pd	0x240(%%rdx),%%zmm15,%%zmm12\n\t	vfnmadd231pd	0x040(%%r11),%%zmm23,%%zmm20\n\t	vfnmadd231pd	0x240(%%r11),%%zmm31,%%zmm28\n\t"\
	" vfmadd231pd	0x040(%%rdx),%%zmm6,%%zmm5	\n\t	 vfmadd231pd	0x240(%%rdx),%%zmm14,%%zmm13\n\t	 vfmadd231pd	0x040(%%r11),%%zmm22,%%zmm21\n\t	 vfmadd231pd	0x240(%%r11),%%zmm30,%%zmm29\n\t"\
		"vmovaps		%%zmm5,%%zmm7			\n\t		vmovaps		%%zmm13,%%zmm15				\n\t		vmovaps		%%zmm21,%%zmm23		\n\t		vmovaps		%%zmm29,%%zmm31		\n\t"\
		"vmovaps		%%zmm4,%%zmm6			\n\t		vmovaps		%%zmm12,%%zmm14				\n\t		vmovaps		%%zmm20,%%zmm22		\n\t		vmovaps		%%zmm28,%%zmm30		\n\t"\
		"vsubpd			 (%%rcx),%%zmm4,%%zmm4	\n\t		vsubpd		0x200(%%rcx),%%zmm12,%%zmm12\n\t		vsubpd			 (%%r13),%%zmm20,%%zmm20\n\t		vsubpd		0x200(%%r13),%%zmm28,%%zmm28\n\t"\
		"vsubpd		0x040(%%rcx),%%zmm5,%%zmm5	\n\t		vsubpd		0x240(%%rcx),%%zmm13,%%zmm13\n\t		vsubpd		0x040(%%r13),%%zmm21,%%zmm21\n\t		vsubpd		0x240(%%r13),%%zmm29,%%zmm29\n\t"\
		"vaddpd			 (%%rcx),%%zmm6,%%zmm6	\n\t		vaddpd		0x200(%%rcx),%%zmm14,%%zmm14\n\t		vaddpd			 (%%r13),%%zmm22,%%zmm22\n\t		vaddpd		0x200(%%r13),%%zmm30,%%zmm30\n\t"\
		"vaddpd		0x040(%%rcx),%%zmm7,%%zmm7	\n\t		vaddpd		0x240(%%rcx),%%zmm15,%%zmm15\n\t		vaddpd		0x040(%%r13),%%zmm23,%%zmm23\n\t		vaddpd		0x240(%%r13),%%zmm31,%%zmm31\n\t"\
		"vsubpd		%%zmm6,%%zmm0,%%zmm0		\n\t		vsubpd		%%zmm14,%%zmm8,%%zmm8	\n\t		vsubpd		%%zmm22,%%zmm16,%%zmm16\n\t		vsubpd		%%zmm30,%%zmm24,%%zmm24\n\t"\
		"vsubpd		%%zmm7,%%zmm1,%%zmm1		\n\t		vsubpd		%%zmm15,%%zmm9,%%zmm9	\n\t		vsubpd		%%zmm23,%%zmm17,%%zmm17\n\t		vsubpd		%%zmm31,%%zmm25,%%zmm25\n\t"\
		"vsubpd		%%zmm5,%%zmm2,%%zmm2		\n\t		vsubpd		%%zmm13,%%zmm10,%%zmm10\n\t		vsubpd		%%zmm21,%%zmm18,%%zmm18\n\t		vsubpd		%%zmm29,%%zmm26,%%zmm26\n\t"\
		"vsubpd		%%zmm4,%%zmm3,%%zmm3		\n\t		vsubpd		%%zmm12,%%zmm11,%%zmm11\n\t		vsubpd		%%zmm20,%%zmm19,%%zmm19\n\t		vsubpd		%%zmm28,%%zmm27,%%zmm27\n\t"\
	"vmovaps	%%zmm12,(%%rcx) \n\t"/* spill zmm12 to make room for two */"vmovaps (%%rdi),%%zmm12	\n\t	vmovaps	%%zmm28,(%%r13) \n\t"/* spill zmm28 to make room for two */"vmovaps (%%rdi),%%zmm28	\n\t"/* two */\
	"vfmadd132pd	%%zmm12,%%zmm0,%%zmm6		\n\t	vfmadd132pd	%%zmm12,%%zmm8 ,%%zmm14			\n\t	vfmadd132pd	%%zmm28,%%zmm16,%%zmm22		\n\t	vfmadd132pd	%%zmm28,%%zmm24,%%zmm30			\n\t"\
	"vfmadd132pd	%%zmm12,%%zmm1,%%zmm7		\n\t	vfmadd132pd	%%zmm12,%%zmm9 ,%%zmm15			\n\t	vfmadd132pd	%%zmm28,%%zmm17,%%zmm23		\n\t	vfmadd132pd	%%zmm28,%%zmm25,%%zmm31			\n\t"\
	"vfmadd132pd	%%zmm12,%%zmm2,%%zmm5		\n\t	vfmadd132pd	%%zmm12,%%zmm10,%%zmm13			\n\t	vfmadd132pd	%%zmm28,%%zmm18,%%zmm21		\n\t	vfmadd132pd	%%zmm28,%%zmm26,%%zmm29			\n\t"\
	"vfmadd132pd	%%zmm12,%%zmm3,%%zmm4		\n\t	vfmadd132pd	(%%rcx),%%zmm11,%%zmm12			\n\t	vfmadd132pd	%%zmm28,%%zmm19,%%zmm20		\n\t	vfmadd132pd	(%%r13),%%zmm27,%%zmm28			\n\t"\
		"												vmovaps		%%zmm14,0x200(%%rcx)			\n\t														vmovaps		%%zmm30,0x200(%%r13)			\n\t"\
		"												vmovaps		%%zmm15,0x240(%%rcx)			\n\t														vmovaps		%%zmm31,0x240(%%r13)			\n\t"\
		"												vmovaps		%%zmm10,%%zmm14					\n\t														vmovaps		%%zmm26,%%zmm30					\n\t"\
		"												vmovaps		%%zmm13,%%zmm15					\n\t														vmovaps		%%zmm29,%%zmm31					\n\t"\
		"												vsubpd		%%zmm12,%%zmm10,%%zmm10	\n\t														vsubpd		%%zmm28,%%zmm26,%%zmm26	\n\t"\
		"												vsubpd		%%zmm11,%%zmm13,%%zmm13	\n\t														vsubpd		%%zmm27,%%zmm29,%%zmm29	\n\t"\
		"												vaddpd		%%zmm12,%%zmm14,%%zmm14	\n\t														vaddpd		%%zmm28,%%zmm30,%%zmm30	\n\t"\
		"												vaddpd		%%zmm11,%%zmm15,%%zmm15	\n\t														vaddpd		%%zmm27,%%zmm31,%%zmm31	\n\t"\
		/*****	SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r10) *****/												/*****	SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r30) *****/\
		"												vmovaps		0x200(%%rcx),%%zmm11			\n\t														vmovaps		0x200(%%r13),%%zmm27			\n\t"\
		"												vmovaps		0x240(%%rcx),%%zmm12			\n\t														vmovaps		0x240(%%r13),%%zmm28			\n\t"\
	"vmovaps %%zmm13,(%%rcx) \n\t"/* spill zmm13 to make room for isrt2 */"vmovaps (%%rsi),%%zmm13	\n\t	vmovaps %%zmm29,(%%r13) \n\t"/* spill zmm29 to make room for isrt2 */"vmovaps (%%rsi),%%zmm29	\n\t"/*isrt2*/\
		"vsubpd		%%zmm11,%%zmm6,%%zmm6		\n\t	vfnmadd231pd	%%zmm13,%%zmm10,%%zmm2		\n\t		vsubpd		%%zmm27,%%zmm22,%%zmm22	\n\t	vfnmadd231pd	%%zmm29,%%zmm26,%%zmm18		\n\t"\
		"vsubpd		%%zmm9,%%zmm0,%%zmm0		\n\t	vfnmadd231pd	%%zmm13,%%zmm15,%%zmm5		\n\t		vsubpd		%%zmm25,%%zmm16,%%zmm16	\n\t	vfnmadd231pd	%%zmm29,%%zmm31,%%zmm21		\n\t"\
		"vsubpd		%%zmm12,%%zmm7,%%zmm7		\n\t	vfnmadd231pd	%%zmm13,%%zmm14,%%zmm4		\n\t		vsubpd		%%zmm28,%%zmm23,%%zmm23	\n\t	vfnmadd231pd	%%zmm29,%%zmm30,%%zmm20		\n\t"\
		"vsubpd		%%zmm8 ,%%zmm1,%%zmm1		\n\t	vfnmadd231pd	(%%rcx),%%zmm13,%%zmm3		\n\t		vsubpd		%%zmm24,%%zmm17,%%zmm17	\n\t	vfnmadd231pd	(%%r13),%%zmm29,%%zmm19		\n\t"\
		"vmovaps		%%zmm6,0x200(%%rcx)		\n\t		vmovaps		%%zmm2,0x280(%%rcx)	\n\t		vmovaps		%%zmm22,0x200(%%r13)	\n\t		vmovaps		%%zmm18,0x280(%%r13)	\n\t"\
		"vmovaps		%%zmm0,0x100(%%rcx)		\n\t		vmovaps		%%zmm5,0x180(%%rcx)	\n\t		vmovaps		%%zmm16,0x100(%%r13)	\n\t		vmovaps		%%zmm21,0x180(%%r13)	\n\t"\
		"vmovaps		%%zmm7,0x240(%%rcx)		\n\t		vmovaps		%%zmm4,0x2c0(%%rcx)	\n\t		vmovaps		%%zmm23,0x240(%%r13)	\n\t		vmovaps		%%zmm20,0x2c0(%%r13)	\n\t"\
		"vmovaps		%%zmm1,0x340(%%rcx)		\n\t		vmovaps		%%zmm3,0x3c0(%%rcx)	\n\t		vmovaps		%%zmm17,0x340(%%r13)	\n\t		vmovaps		%%zmm19,0x3c0(%%r13)	\n\t"\
	"vmovaps	0x80(%%rdi),%%zmm13																	\n\t	vmovaps	0x80(%%rdi),%%zmm29 \n\t"/* sqrt2 *//* rcol: *2 => *sqrt2 due to above FMA leaving zmm26,13-15 unmultiplied-by-isrt2: */\
	"vfmadd132pd	(%%rdi),%%zmm6,%%zmm11		\n\t	 vfmadd132pd	%%zmm13,%%zmm2,%%zmm10	\n\t	vfmadd132pd	(%%rdi),%%zmm22,%%zmm27	\n\t	 vfmadd132pd	%%zmm29,%%zmm18,%%zmm26	\n\t"\
	"vfmadd132pd	(%%rdi),%%zmm0,%%zmm9		\n\t	 vfmadd132pd	%%zmm13,%%zmm5,%%zmm15	\n\t	vfmadd132pd	(%%rdi),%%zmm16,%%zmm25	\n\t	 vfmadd132pd	%%zmm29,%%zmm21,%%zmm31	\n\t"\
	"vfmadd132pd	(%%rdi),%%zmm7,%%zmm12		\n\t	 vfmadd132pd	%%zmm13,%%zmm4,%%zmm14	\n\t	vfmadd132pd	(%%rdi),%%zmm23,%%zmm28	\n\t	 vfmadd132pd	%%zmm29,%%zmm20,%%zmm30	\n\t"\
	"vfmadd132pd	(%%rdi),%%zmm1,%%zmm8		\n\t	 vfmadd132pd	(%%rcx),%%zmm3,%%zmm13	\n\t	vfmadd132pd	(%%rdi),%%zmm17,%%zmm24	\n\t	 vfmadd132pd	(%%r13),%%zmm19,%%zmm29	\n\t"\
		"vmovaps		%%zmm11,     (%%rcx)	\n\t		vmovaps		%%zmm10,0x080(%%rcx)	\n\t		vmovaps		%%zmm27,     (%%r13)	\n\t		vmovaps		%%zmm26,0x080(%%r13)	\n\t"\
		"vmovaps		%%zmm9 ,0x300(%%rcx)	\n\t		vmovaps		%%zmm15,0x380(%%rcx)	\n\t		vmovaps		%%zmm25,0x300(%%r13)	\n\t		vmovaps		%%zmm31,0x380(%%r13)	\n\t"\
		"vmovaps		%%zmm12,0x040(%%rcx)	\n\t		vmovaps		%%zmm14,0x0c0(%%rcx)	\n\t		vmovaps		%%zmm28,0x040(%%r13)	\n\t		vmovaps		%%zmm30,0x0c0(%%r13)	\n\t"\
		"vmovaps		%%zmm8 ,0x140(%%rcx)	\n\t		vmovaps		%%zmm13,0x1c0(%%rcx)	\n\t		vmovaps		%%zmm24,0x140(%%r13)	\n\t		vmovaps		%%zmm29,0x1c0(%%r13)	\n\t"\
		"\n\t"\
	/**********************************************************************************/\
	/*...and now do eight radix-4 transforms, including the internal twiddle factors: */\
	/**********************************************************************************/\
		"movq		%[__isrt2],%%rsi			\n\t"\
	/*...Block 1: t00,t10,t20,t30	*/					/*...Block 5: t08,t18,t28,t38	*/					/*...Block 3: t04,t14,t24,t34	*/			/*...Block 7: t0C,t1C,t2C,t3C	*/\
		"movq		%[__r20],%%rcx				\n\t					leaq	0x100(%%rcx),%%r12	\n\t		vmovaps		 (%%r12),%%zmm20			\n\t		vmovaps		0x200(%%r12),%%zmm28		\n\t"\
		"movq		%[__r30],%%rdx				\n\t					leaq	0x100(%%rdx),%%r13	\n\t		vmovaps	0x040(%%r12),%%zmm21			\n\t		vmovaps		0x240(%%r12),%%zmm29		\n\t"\
		"movq		%[__r00],%%rax				\n\t					leaq	0x100(%%rax),%%r10	\n\t		vmovaps		 (%%r13),%%zmm22			\n\t		vmovaps		0x200(%%r13),%%zmm30		\n\t"\
		"movq		%[__r10],%%rbx				\n\t					leaq	0x100(%%rbx),%%r11	\n\t		vmovaps	0x040(%%r13),%%zmm19			\n\t		vmovaps		0x240(%%r13),%%zmm27		\n\t"\
		"vmovaps		 (%%rax),%%zmm0			\n\t		vmovaps		0x200(%%rax),%%zmm8			\n\t	vmovaps 0x080(%%rsi),%%zmm23\n\t vmovaps 0x040(%%rsi),%%zmm31 \n\t"/* Use 7,15 for s,c 2 cut loads */\
		"vmovaps	0x040(%%rax),%%zmm1			\n\t		vmovaps		0x240(%%rax),%%zmm9			\n\t		vmulpd		%%zmm23,%%zmm21,%%zmm17	\n\t		vmulpd		 %%zmm31,%%zmm29,%%zmm25	\n\t"\
		"vmovaps		 (%%rbx),%%zmm2			\n\t		vmovaps		0x200(%%rbx),%%zmm10		\n\t		vmulpd		%%zmm31,%%zmm19,%%zmm19	\n\t		vmulpd		 %%zmm23,%%zmm27,%%zmm27	\n\t"\
		"vmovaps	0x040(%%rbx),%%zmm3			\n\t		vmovaps		0x240(%%rbx),%%zmm11		\n\t		vmulpd		%%zmm23,%%zmm20,%%zmm16	\n\t		vmulpd		 %%zmm31,%%zmm28,%%zmm24	\n\t"\
		"vaddpd			  %%zmm0,%%zmm2,%%zmm2	\n\t		vaddpd		%%zmm9 ,%%zmm10,%%zmm10		\n\t		vmulpd		%%zmm31,%%zmm22,%%zmm18	\n\t		vmulpd		 %%zmm23,%%zmm30,%%zmm26	\n\t"\
		"vaddpd			  %%zmm1,%%zmm3,%%zmm3	\n\t		vaddpd		%%zmm8 ,%%zmm11,%%zmm11		\n\t	vfmsub132pd	%%zmm31,%%zmm17,%%zmm20		\n\t	vfmsub132pd		 %%zmm23,%%zmm25,%%zmm28	\n\t"\
		"vsubpd			 (%%rbx),%%zmm0,%%zmm0	\n\t		vsubpd		0x240(%%rbx),%%zmm8,%%zmm8	\n\t	vfmsub132pd	%%zmm23,%%zmm19,%%zmm22		\n\t	vfmsub132pd		 %%zmm31,%%zmm27,%%zmm30	\n\t"\
		"vsubpd		0x040(%%rbx),%%zmm1,%%zmm1	\n\t		vsubpd		0x200(%%rbx),%%zmm9,%%zmm9	\n\t	vfmadd132pd	%%zmm31,%%zmm16,%%zmm21		\n\t	vfmadd132pd		 %%zmm23,%%zmm24,%%zmm29	\n\t"\
		"vmovaps		 (%%rcx),%%zmm4			\n\t		vmovaps		0x200(%%rcx),%%zmm12		\n\t	vfmadd132pd 0x40(%%r13),%%zmm18,%%zmm23		\n\t	vfmadd132pd 0x240(%%r13),%%zmm26,%%zmm31	\n\t"\
		"vmovaps	0x040(%%rcx),%%zmm5			\n\t		vmovaps		0x240(%%rcx),%%zmm13		\n\t		vsubpd		%%zmm22,%%zmm20,%%zmm20	\n\t		vsubpd		%%zmm30,%%zmm28,%%zmm28\n\t"\
		"vmovaps		 (%%rdx),%%zmm6			\n\t		vmovaps		0x200(%%rdx),%%zmm14		\n\t		vsubpd		%%zmm23,%%zmm21,%%zmm21	\n\t		vsubpd		%%zmm31,%%zmm29,%%zmm29\n\t"\
		"vmovaps	0x040(%%rdx),%%zmm7			\n\t		vmovaps		0x240(%%rdx),%%zmm15		\n\t		vmovaps		 (%%r11),%%zmm18			\n\t		vmovaps		0x200(%%r11),%%zmm26		\n\t"\
		"vaddpd		%%zmm4,%%zmm6,%%zmm6		\n\t		vsubpd		%%zmm13,%%zmm12,%%zmm12		\n\t		vmovaps	0x040(%%r11),%%zmm19			\n\t		vmovaps		0x240(%%r11),%%zmm27		\n\t"\
		"vaddpd		%%zmm5,%%zmm7,%%zmm7		\n\t		vaddpd		0x200(%%rcx),%%zmm13,%%zmm13\n\t	vfmadd132pd	(%%rdi),%%zmm20,%%zmm22	\n\t	vfmadd132pd		(%%rdi),%%zmm28,%%zmm30\n\t"\
		"vsubpd			 (%%rdx),%%zmm4,%%zmm4	\n\t		vaddpd			 %%zmm15,%%zmm14,%%zmm14\n\t	vfmadd132pd	(%%rdi),%%zmm21,%%zmm23	\n\t	vfmadd132pd		(%%rdi),%%zmm29,%%zmm31\n\t"\
		"vsubpd		0x040(%%rdx),%%zmm5,%%zmm5	\n\t		vsubpd		0x200(%%rdx),%%zmm15,%%zmm15\n\t		vsubpd		0x040(%%r11),%%zmm18,%%zmm18\n\t		vaddpd		0x240(%%r11),%%zmm26,%%zmm26\n\t"\
		"vsubpd		%%zmm6,%%zmm2,%%zmm2		\n\t		vmulpd		(%%rsi),%%zmm12,%%zmm12\n\t		vaddpd			 (%%r11),%%zmm19,%%zmm19\n\t		vsubpd		0x200(%%r11),%%zmm27,%%zmm27\n\t"\
		"vsubpd		%%zmm7,%%zmm3,%%zmm3		\n\t		vmulpd		(%%rsi),%%zmm13,%%zmm13\n\t	vmovaps %%zmm22,(%%r12) \n\t"/* spill zmm22 to make room for isrt2 */"vmovaps (%%rsi),%%zmm22	\n\t"\
	"vfmadd132pd	(%%rdi),%%zmm2,%%zmm6		\n\t	vfnmadd231pd	(%%rsi),%%zmm14,%%zmm12\n\t		vmovaps		 (%%r10),%%zmm16			\n\t		vmovaps		 0x200(%%r10),%%zmm24		\n\t"\
	"vfmadd132pd	(%%rdi),%%zmm3,%%zmm7		\n\t	vfnmadd231pd	(%%rsi),%%zmm15,%%zmm13\n\t		vmovaps	0x040(%%r10),%%zmm17			\n\t		vmovaps		 0x240(%%r10),%%zmm25		\n\t"\
		"vmovaps		%%zmm2,      (%%rcx)	\n\t	 vfmadd132pd	0x80(%%rdi),%%zmm12,%%zmm14\n\t	vfnmadd231pd	%%zmm22,%%zmm18,%%zmm16		\n\t	vfnmadd231pd	%%zmm22,%%zmm26,%%zmm24		\n\t"\
		"vmovaps		%%zmm3, 0x040(%%rcx)	\n\t	 vfmadd132pd	0x80(%%rdi),%%zmm13,%%zmm15\n\t	vfnmadd231pd	%%zmm22,%%zmm19,%%zmm17		\n\t	vfnmadd231pd	%%zmm22,%%zmm27,%%zmm25		\n\t"\
		"										\n\t		vsubpd		%%zmm12,%%zmm8,%%zmm8	\n\t	 vfmadd213pd		 (%%r10),%%zmm22,%%zmm18\n\t	 vfmadd213pd	0x200(%%r10),%%zmm22,%%zmm26\n\t"\
		"										\n\t		vsubpd		%%zmm13,%%zmm10,%%zmm10\n\t	 vfmadd213pd	0x040(%%r10),%%zmm22,%%zmm19\n\t	 vfmadd213pd	0x240(%%r10),%%zmm22,%%zmm27\n\t"\
		"vmovaps		%%zmm6,      (%%rax)	\n\t	 vfmadd132pd	(%%rdi),%%zmm8,%%zmm12\n\t	vmovaps	(%%r12),%%zmm22	\n\t"/* restore spill */\
		"vmovaps		%%zmm7, 0x040(%%rax)	\n\t	 vfmadd132pd	(%%rdi),%%zmm10,%%zmm13\n\t		vsubpd		%%zmm22,%%zmm18,%%zmm18	\n\t		vsubpd		%%zmm28,%%zmm24,%%zmm24\n\t"\
		"vsubpd		%%zmm5,%%zmm0,%%zmm0		\n\t		vmovaps		%%zmm8 ,0x200(%%rcx)		\n\t		vsubpd		%%zmm23,%%zmm19,%%zmm19	\n\t		vsubpd		%%zmm29,%%zmm25,%%zmm25\n\t"\
		"vsubpd		%%zmm4,%%zmm1,%%zmm1		\n\t		vmovaps		%%zmm10,0x240(%%rcx)		\n\t	vfmadd132pd	(%%rdi),%%zmm18,%%zmm22	\n\t	vfmadd132pd		(%%rdi),%%zmm24,%%zmm28\n\t"\
	"vfmadd132pd	(%%rdi),%%zmm0,%%zmm5		\n\t		vmovaps		%%zmm12,0x200(%%rax)		\n\t	vfmadd132pd	(%%rdi),%%zmm19,%%zmm23	\n\t	vfmadd132pd		(%%rdi),%%zmm25,%%zmm29\n\t"\
	"vfmadd132pd	(%%rdi),%%zmm1,%%zmm4		\n\t		vmovaps		%%zmm13,0x240(%%rax)		\n\t		vmovaps		%%zmm18,      (%%r12)		\n\t		vmovaps		%%zmm24,0x200(%%r12)		\n\t"\
		"vmovaps		%%zmm0,      (%%rbx)	\n\t		vsubpd		%%zmm15,%%zmm11,%%zmm11\n\t		vmovaps		%%zmm19, 0x040(%%r12)		\n\t		vmovaps		%%zmm25,0x240(%%r12)		\n\t"\
		"vmovaps		%%zmm1, 0x040(%%rdx)	\n\t		vsubpd		%%zmm14,%%zmm9,%%zmm9	\n\t		vmovaps		%%zmm22,      (%%r10)		\n\t		vmovaps		%%zmm28,0x200(%%r10)		\n\t"\
		"										\n\t	 vfmadd132pd	(%%rdi),%%zmm11,%%zmm15\n\t		vmovaps		%%zmm23, 0x040(%%r10)		\n\t		vmovaps		%%zmm29,0x240(%%r10)		\n\t"\
		"										\n\t	 vfmadd132pd	(%%rdi),%%zmm9 ,%%zmm14\n\t		vsubpd		%%zmm21,%%zmm16,%%zmm16	\n\t		vsubpd		%%zmm31,%%zmm26,%%zmm26\n\t"\
		"vmovaps		%%zmm5,      (%%rdx)	\n\t		vmovaps		%%zmm11,0x200(%%rbx)		\n\t		vsubpd		%%zmm20,%%zmm17,%%zmm17	\n\t		vsubpd		%%zmm30,%%zmm27,%%zmm27\n\t"\
		"vmovaps		%%zmm4, 0x040(%%rbx)	\n\t		vmovaps		%%zmm9 ,0x240(%%rdx)		\n\t	vfmadd132pd	(%%rdi),%%zmm16,%%zmm21	\n\t	vfmadd132pd		(%%rdi),%%zmm26,%%zmm31\n\t"\
		"													vmovaps		%%zmm15,0x200(%%rdx)		\n\t	vfmadd132pd	(%%rdi),%%zmm17,%%zmm20	\n\t	vfmadd132pd		(%%rdi),%%zmm27,%%zmm30\n\t"\
		"													vmovaps		%%zmm14,0x240(%%rbx)		\n\t		vmovaps		%%zmm16,      (%%r11)		\n\t		vmovaps		%%zmm26,0x200(%%r11)		\n\t"\
		"																										vmovaps		%%zmm17, 0x040(%%r13)		\n\t		vmovaps		%%zmm27,0x240(%%r13)		\n\t"\
		"																										vmovaps		%%zmm21,      (%%r13)		\n\t		vmovaps		%%zmm31,0x200(%%r13)		\n\t"\
		"																										vmovaps		%%zmm20, 0x040(%%r11)		\n\t		vmovaps		%%zmm30,0x240(%%r11)		\n\t"\
		"addq		$0x0c0,%%rsi	\n\t"/* cc1 */\
	/*...Block 2: t02,t12,t22,t32	*/					/*...Block 6: t0A,t1A,t2A,t3A */					/*...Block 4: t06,t16,t26,t36	*/			/*...Block 8: t0E,t1E,t2E,t3E	*/\
		"addq		$0x080,%%rcx				\n\t					addq	$0x080,%%r12	\n\t		vmovaps		 (%%r12),%%zmm20			\n\t		vmovaps		0x200(%%r12),%%zmm28		\n\t"\
		"addq		$0x080,%%rdx				\n\t					addq	$0x080,%%r13	\n\t		vmovaps		 (%%r13),%%zmm22			\n\t		vmovaps		0x200(%%r13),%%zmm30		\n\t"\
		"addq		$0x080,%%rax				\n\t					addq	$0x080,%%r10	\n\t		vmovaps	0x040(%%r12),%%zmm17			\n\t		vmovaps		0x240(%%r12),%%zmm29		\n\t"\
		"addq		$0x080,%%rbx				\n\t					addq	$0x080,%%r11	\n\t		vmovaps	0x040(%%r13),%%zmm19			\n\t		vmovaps		0x240(%%r13),%%zmm27		\n\t"\
		"vmovaps		 (%%rcx),%%zmm4			\n\t		vmovaps		0x200(%%rcx),%%zmm12		\n\t		vmovaps		 (%%rsi),%%zmm18			\n\t		vmovaps		0x040(%%rsi),%%zmm23		\n\t"\
		"vmovaps		 (%%rdx),%%zmm6			\n\t		vmovaps		0x200(%%rdx),%%zmm14		\n\t		vmovaps	0x080(%%rsi),%%zmm21			\n\t		vmovaps		0x0c0(%%rsi),%%zmm31		\n\t"\
		"vmovaps	0x040(%%rcx),%%zmm5			\n\t		vmovaps		0x240(%%rcx),%%zmm9			\n\t		vmulpd		%%zmm31,%%zmm17,%%zmm17	\n\t		vmulpd			%%zmm18,%%zmm29,%%zmm25	\n\t"\
		"vmovaps	0x040(%%rdx),%%zmm3			\n\t		vmovaps		0x240(%%rdx),%%zmm11		\n\t		vmulpd		%%zmm18,%%zmm19,%%zmm19	\n\t		vmulpd			%%zmm21,%%zmm27,%%zmm27	\n\t"\
		"vmovaps		 (%%rsi),%%zmm15		\n\t		vmovaps		0x040(%%rsi),%%zmm10		\n\t		vmulpd		%%zmm31,%%zmm20,%%zmm16	\n\t		vmulpd			%%zmm18,%%zmm28,%%zmm24	\n\t"\
		"vmovaps	0x080(%%rsi),%%zmm7 		\n\t		vmovaps		0x0c0(%%rsi),%%zmm13		\n\t		vmulpd		 (%%r13),%%zmm18,%%zmm18	\n\t		vmulpd			%%zmm21,%%zmm30,%%zmm26	\n\t"\
		"vmulpd		%%zmm10,%%zmm5,%%zmm1		\n\t		vmulpd			%%zmm7 ,%%zmm9 ,%%zmm9	\n\t	vfmsub132pd	%%zmm21,%%zmm17,%%zmm20		\n\t	vfmsub132pd			%%zmm23,%%zmm25,%%zmm28	\n\t"\
		"vmulpd		%%zmm13,%%zmm3,%%zmm3		\n\t		vmulpd			%%zmm10,%%zmm11,%%zmm11	\n\t	vfmadd132pd	%%zmm23,%%zmm19,%%zmm22		\n\t	vfmsub132pd			%%zmm31,%%zmm27,%%zmm30	\n\t"\
		"vmulpd		%%zmm10,%%zmm4,%%zmm0		\n\t		vmulpd			%%zmm7 ,%%zmm12,%%zmm8	\n\t	vfmadd132pd	0x040(%%r12),%%zmm16,%%zmm21	\n\t	vfmadd132pd			%%zmm23,%%zmm24,%%zmm29	\n\t"\
		"vmulpd		%%zmm13,%%zmm6,%%zmm2		\n\t		vmulpd		0x200(%%rdx),%%zmm10,%%zmm10\n\t	vfmsub132pd	0x040(%%r13),%%zmm18,%%zmm23	\n\t	vfmadd132pd		0x240(%%r13),%%zmm26,%%zmm31\n\t"\
	"vfmsub132pd	%%zmm15,%%zmm1,%%zmm4		\n\t	vfmsub132pd			%%zmm13,%%zmm9 ,%%zmm12	\n\t		vsubpd		%%zmm22,%%zmm20,%%zmm20	\n\t		vsubpd		%%zmm30,%%zmm28,%%zmm28\n\t"\
	"vfmsub132pd	%%zmm7 ,%%zmm3,%%zmm6		\n\t	vfmadd132pd			%%zmm15,%%zmm11,%%zmm14	\n\t		vsubpd		%%zmm23,%%zmm21,%%zmm21	\n\t		vsubpd		%%zmm31,%%zmm29,%%zmm29\n\t"\
	"vfmadd132pd	%%zmm15,%%zmm0,%%zmm5		\n\t	vfmadd132pd		0x240(%%rcx),%%zmm8 ,%%zmm13\n\t		vmovaps		 (%%r11),%%zmm18			\n\t		vmovaps		0x200(%%r11),%%zmm26		\n\t"\
	"vfmadd132pd	0x040(%%rdx),%%zmm2,%%zmm7	\n\t	vfmsub132pd		0x240(%%rdx),%%zmm10,%%zmm15\n\t		vmovaps	0x040(%%r11),%%zmm16			\n\t		vmovaps		0x240(%%r11),%%zmm24		\n\t"\
		"vsubpd		%%zmm6,%%zmm4,%%zmm4		\n\t		vsubpd		%%zmm14,%%zmm12,%%zmm12\n\t	vfmadd132pd	(%%rdi),%%zmm20,%%zmm22	\n\t	vfmadd132pd		(%%rdi),%%zmm28,%%zmm30\n\t"\
		"vsubpd		%%zmm7,%%zmm5,%%zmm5		\n\t		vsubpd		%%zmm15,%%zmm13,%%zmm13\n\t	vfmadd132pd	(%%rdi),%%zmm21,%%zmm23	\n\t	vfmadd132pd		(%%rdi),%%zmm29,%%zmm31\n\t"\
		"vmovaps		 (%%rbx),%%zmm2			\n\t		vmovaps		0x200(%%rbx),%%zmm10		\n\t		vmovaps		 (%%r11),%%zmm18			\n\t		vmovaps		0x200(%%r11),%%zmm26		\n\t"\
		"vmovaps	0x040(%%rbx),%%zmm0			\n\t		vmovaps		0x240(%%rbx),%%zmm8			\n\t		vmovaps	0x040(%%r11),%%zmm16			\n\t		vmovaps		0x240(%%r11),%%zmm24		\n\t"\
	"vfmadd132pd	(%%rdi),%%zmm4,%%zmm6		\n\t	vfmadd132pd		(%%rdi),%%zmm12,%%zmm14\n\t		vmovaps		 (%%r11),%%zmm17			\n\t		vmovaps		0x200(%%r11),%%zmm25		\n\t"\
	"vfmadd132pd	(%%rdi),%%zmm5,%%zmm7		\n\t	vfmadd132pd		(%%rdi),%%zmm13,%%zmm15\n\t		vmovaps	0x040(%%r11),%%zmm19			\n\t		vmovaps		0x240(%%r11),%%zmm27		\n\t"\
		"vmovaps		 (%%rbx),%%zmm1			\n\t		vmovaps		0x200(%%rbx),%%zmm9			\n\t		vmovaps		 (%%r11),%%zmm17			\n\t		vmovaps		0x200(%%r11),%%zmm25		\n\t"\
		"vmovaps	-0x80(%%rsi),%%zmm3			\n\t		vmovaps		-0x40(%%rsi),%%zmm11		\n\t		vmovaps	-0x40(%%rsi),%%zmm19			\n\t		vmovaps		-0x80(%%rsi),%%zmm27		\n\t"\
		"vmulpd		%%zmm11,%%zmm0,%%zmm0		\n\t		vmulpd			%%zmm3 ,%%zmm8 ,%%zmm8	\n\t		vmulpd		%%zmm27,%%zmm16,%%zmm16	\n\t		vmulpd			%%zmm19,%%zmm24,%%zmm24	\n\t"\
		"vmulpd		%%zmm11,%%zmm1,%%zmm1		\n\t		vmulpd			%%zmm3 ,%%zmm9 ,%%zmm9	\n\t		vmulpd		%%zmm27,%%zmm17,%%zmm17	\n\t		vmulpd			%%zmm19,%%zmm25,%%zmm25	\n\t"\
	"vfmsub132pd	%%zmm3 ,%%zmm0,%%zmm2		\n\t	vfmadd132pd			%%zmm11,%%zmm8 ,%%zmm10	\n\t	vfmsub132pd		%%zmm19,%%zmm16,%%zmm18		\n\t	vfmadd132pd			%%zmm27,%%zmm24,%%zmm26	\n\t"\
	"vfmadd132pd	0x40(%%rbx),%%zmm1,%%zmm3	\n\t	vfmsub132pd		0x240(%%rbx),%%zmm9,%%zmm11	\n\t	vfmadd132pd	0x40(%%r11),%%zmm17,%%zmm19		\n\t	vfmsub132pd		0x240(%%r11),%%zmm25,%%zmm27\n\t"\
		"vmovaps		 (%%rax),%%zmm0			\n\t		vmovaps		0x200(%%rax),%%zmm8			\n\t		vmovaps		 (%%r10),%%zmm16			\n\t		vmovaps		0x200(%%r10),%%zmm24		\n\t"\
		"vmovaps	0x040(%%rax),%%zmm1			\n\t		vmovaps		0x240(%%rax),%%zmm9			\n\t		vmovaps	0x040(%%r10),%%zmm17			\n\t		vmovaps		0x240(%%r10),%%zmm25		\n\t"\
		"vsubpd		%%zmm2,%%zmm0,%%zmm0		\n\t		vsubpd		%%zmm10,%%zmm8,%%zmm8	\n\t		vsubpd		%%zmm18,%%zmm16,%%zmm16	\n\t		vsubpd		%%zmm26,%%zmm24,%%zmm24\n\t"\
		"vsubpd		%%zmm3,%%zmm1,%%zmm1		\n\t		vsubpd		%%zmm11,%%zmm9,%%zmm9	\n\t		vsubpd		%%zmm19,%%zmm17,%%zmm17	\n\t		vsubpd		%%zmm27,%%zmm25,%%zmm25\n\t"\
		"vaddpd			 (%%rax),%%zmm2,%%zmm2	\n\t		vaddpd		0x200(%%rax),%%zmm10,%%zmm10\n\t		vaddpd			 (%%r10),%%zmm18,%%zmm18\n\t		vaddpd		0x200(%%r10),%%zmm26,%%zmm26\n\t"\
		"vaddpd		0x040(%%rax),%%zmm3,%%zmm3	\n\t		vaddpd		0x240(%%rax),%%zmm11,%%zmm11\n\t		vaddpd		0x040(%%r10),%%zmm19,%%zmm19\n\t		vaddpd		0x240(%%r10),%%zmm27,%%zmm27\n\t"\
		"vsubpd		%%zmm6,%%zmm2,%%zmm2		\n\t		vsubpd		%%zmm12,%%zmm8,%%zmm8	\n\t		vsubpd		%%zmm20,%%zmm18,%%zmm18	\n\t		vsubpd		%%zmm28,%%zmm24,%%zmm24\n\t"\
		"vsubpd		%%zmm7,%%zmm3,%%zmm3		\n\t		vsubpd		%%zmm13,%%zmm9,%%zmm9	\n\t		vsubpd		%%zmm21,%%zmm19,%%zmm19	\n\t		vsubpd		%%zmm29,%%zmm25,%%zmm25\n\t"\
	"vfmadd132pd	(%%rdi),%%zmm2,%%zmm6		\n\t	vfmadd132pd		(%%rdi),%%zmm8,%%zmm12\n\t	vfmadd132pd	(%%rdi),%%zmm18,%%zmm20	\n\t	vfmadd132pd		(%%rdi),%%zmm24,%%zmm28\n\t"\
	"vfmadd132pd	(%%rdi),%%zmm3,%%zmm7		\n\t	vfmadd132pd		(%%rdi),%%zmm9,%%zmm13\n\t	vfmadd132pd	(%%rdi),%%zmm19,%%zmm21	\n\t	vfmadd132pd		(%%rdi),%%zmm25,%%zmm29\n\t"\
		"vmovaps		%%zmm2,      (%%rcx)	\n\t		vmovaps		%%zmm8 ,0x200(%%rcx)		\n\t		vmovaps		%%zmm18,      (%%r12)		\n\t		vmovaps		%%zmm24,0x200(%%r12)		\n\t"\
		"vmovaps		%%zmm3, 0x040(%%rcx)	\n\t		vmovaps		%%zmm9 ,0x240(%%rcx)		\n\t		vmovaps		%%zmm19, 0x040(%%r12)		\n\t		vmovaps		%%zmm25,0x240(%%r12)		\n\t"\
		"vmovaps		%%zmm6,      (%%rax)	\n\t		vmovaps		%%zmm12,0x200(%%rax)		\n\t		vmovaps		%%zmm20,      (%%r10)		\n\t		vmovaps		%%zmm28,0x200(%%r10)		\n\t"\
		"vmovaps		%%zmm7, 0x040(%%rax)	\n\t		vmovaps		%%zmm13,0x240(%%rax)		\n\t		vmovaps		%%zmm21, 0x040(%%r10)		\n\t		vmovaps		%%zmm29,0x240(%%r10)		\n\t"\
		"vsubpd		%%zmm5,%%zmm0,%%zmm0		\n\t		vsubpd		%%zmm15,%%zmm10,%%zmm10\n\t		vsubpd		%%zmm23,%%zmm16,%%zmm16	\n\t		vsubpd		%%zmm31,%%zmm26,%%zmm26\n\t"\
		"vsubpd		%%zmm4,%%zmm1,%%zmm1		\n\t		vsubpd		%%zmm14,%%zmm11,%%zmm11\n\t		vsubpd		%%zmm22,%%zmm17,%%zmm17	\n\t		vsubpd		%%zmm30,%%zmm27,%%zmm27\n\t"\
	"vfmadd132pd	(%%rdi),%%zmm0,%%zmm5		\n\t	vfmadd132pd		(%%rdi),%%zmm10,%%zmm15\n\t	vfmadd132pd	(%%rdi),%%zmm16,%%zmm23	\n\t	vfmadd132pd		(%%rdi),%%zmm26,%%zmm31\n\t"\
	"vfmadd132pd	(%%rdi),%%zmm1,%%zmm4		\n\t	vfmadd132pd		(%%rdi),%%zmm11,%%zmm14\n\t	vfmadd132pd	(%%rdi),%%zmm17,%%zmm22	\n\t	vfmadd132pd		(%%rdi),%%zmm27,%%zmm30\n\t"\
		"vmovaps		%%zmm0,      (%%rbx)	\n\t		vmovaps		%%zmm10,0x200(%%rbx)		\n\t		vmovaps		%%zmm16,      (%%r11)		\n\t		vmovaps		%%zmm26,0x200(%%r11)		\n\t"\
		"vmovaps		%%zmm1, 0x040(%%rdx)	\n\t		vmovaps		%%zmm11,0x240(%%rdx)		\n\t		vmovaps		%%zmm17, 0x040(%%r13)		\n\t		vmovaps		%%zmm27,0x240(%%r13)		\n\t"\
		"vmovaps		%%zmm5,      (%%rdx)	\n\t		vmovaps		%%zmm15,0x200(%%rdx)		\n\t		vmovaps		%%zmm23,      (%%r13)		\n\t		vmovaps		%%zmm31,0x200(%%r13)		\n\t"\
		"vmovaps		%%zmm4, 0x040(%%rbx)	\n\t		vmovaps		%%zmm14,0x240(%%rbx)		\n\t		vmovaps		%%zmm22, 0x040(%%r11)		\n\t		vmovaps		%%zmm30,0x240(%%r11)		\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__add1] "m" (Xadd1)\
		 ,[__add2] "m" (Xadd2)\
		 ,[__add3] "m" (Xadd3)\
		 ,[__add4] "m" (Xadd4)\
		 ,[__add5] "m" (Xadd5)\
		 ,[__add6] "m" (Xadd6)\
		 ,[__add7] "m" (Xadd7)\
		 ,[__r00] "m" (Xr00)\
		 ,[__r10] "m" (Xr10)\
		 ,[__r20] "m" (Xr20)\
		 ,[__r30] "m" (Xr30)\
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__cc0] "m" (Xcc0)\
		 ,[__c00] "m" (Xc00)\
		 ,[__c01] "m" (Xc01)\
		 ,[__c02] "m" (Xc02)\
		 ,[__c03] "m" (Xc03)\
		 ,[__c05] "m" (Xc05)\
		 ,[__c07] "m" (Xc07)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r10","r11","r12","r13","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15","xmm16","xmm17","xmm18","xmm19","xmm20","xmm21","xmm22","xmm23","xmm24","xmm25","xmm26","xmm27","xmm28","xmm29","xmm30","xmm31"	/* Clobbered registers */\
	);\
	}

   #endif	// (IMCI512 or AVX512?) toggle

  #endif	// USE_16_REG ?

  #ifdef USE_16_REG

	// Cost [vector-ops only]: 880 MEM, 352 MUL (236 FMA), 204 ADD/SUB
	#define SSE2_RADIX32_WRAPPER_DIT(Xadd0,Xadd1,Xadd2,Xadd3,Xadd4,Xadd5,Xadd6,Xadd7,Xisrt2,Xr00,Xr08,Xr10,Xr20,Xr28,Xr30,Xc01,Xc02,Xc04,Xc06,Xc0A,Xc10,Xc12,Xc1A)\
	{\
	__asm__ volatile (\
	/************************************************************************/\
	/*...And do an inverse DIT radix-32 pass on the squared-data blocks.	*/\
	/************************************************************************/\
	/* Block 1: */\
		"movq	%[__isrt2],%%rsi				\n\t	leaq	0x11c0(%%rsi),%%rdi		\n\t"/* two */\
		"movq		%[__r00],%%rax		\n\t"\
		"leaq		0x800(%%rax),%%rbx			\n\t"\
		"leaq		0x400(%%rax),%%rcx			\n\t"\
		"leaq		0xc00(%%rax),%%rdx			\n\t		vmovaps		 (%%rdi),%%zmm13	\n\t"/* 2.0*/\
	/*****	SSE2_RADIX4_DIT_IN_PLACE()	*****/	/* Block 2 has tmp-addr offset +0x100 w.r.to Block 1:  */\
		"vmovaps		(%%rbx),%%zmm0			\n\t		vmovaps		0x100(%%rbx),%%zmm8 		\n\t"\
		"vmovaps	0x40(%%rbx),%%zmm1			\n\t		vmovaps		0x140(%%rbx),%%zmm9 		\n\t"\
		"vmovaps		(%%rax),%%zmm2			\n\t		vmovaps		0x100(%%rax),%%zmm10		\n\t"\
		"vmovaps	0x40(%%rax),%%zmm3			\n\t		vmovaps		0x140(%%rax),%%zmm11		\n\t"\
		"vmovaps		(%%rdx),%%zmm4			\n\t		vmovaps		0x100(%%rdx),%%zmm12		\n\t"\
		"vmovaps	0x40(%%rdx),%%zmm5			\n\t	/*	vmovaps		0x140(%%rdx),%%zmm13	*/	\n\t"\
		"vmovaps		(%%rcx),%%zmm6			\n\t		vmovaps		0x100(%%rcx),%%zmm14		\n\t"\
		"vmovaps	0x40(%%rcx),%%zmm7			\n\t		vmovaps		0x140(%%rcx),%%zmm15		\n\t"\
		"vsubpd			%%zmm0,%%zmm2,%%zmm2	\n\t		vsubpd		%%zmm8 ,%%zmm10,%%zmm10		\n\t"\
		"vsubpd			%%zmm1,%%zmm3,%%zmm3	\n\t		vsubpd		%%zmm9 ,%%zmm11,%%zmm11		\n\t"\
		"vsubpd			%%zmm4,%%zmm6,%%zmm6	\n\t		vsubpd		%%zmm12,%%zmm14,%%zmm14		\n\t"\
		"vsubpd			%%zmm5,%%zmm7,%%zmm7	\n\t		vsubpd 0x140(%%rdx),%%zmm15,%%zmm15		\n\t"\
		"vfmadd132pd	%%zmm13,%%zmm2,%%zmm0	\n\t	 vfmadd132pd	%%zmm13,%%zmm10,%%zmm8 		\n\t"\
		"vfmadd132pd	%%zmm13,%%zmm3,%%zmm1	\n\t	 vfmadd132pd	%%zmm13,%%zmm11,%%zmm9 		\n\t"\
		"vfmadd132pd	%%zmm13,%%zmm6,%%zmm4	\n\t	 vfmadd132pd	%%zmm13,%%zmm14,%%zmm12		\n\t"\
		"vfmadd132pd	%%zmm13,%%zmm7,%%zmm5	\n\t	 vfmadd132pd 0x140(%%rdx),%%zmm15,%%zmm13	\n\t"\
		"vsubpd			%%zmm4 ,%%zmm0,%%zmm0	\n\t		vsubpd		%%zmm12,%%zmm8 ,%%zmm8 		\n\t"\
		"vsubpd			%%zmm5 ,%%zmm1,%%zmm1	\n\t		vsubpd		%%zmm13,%%zmm9 ,%%zmm9 		\n\t"\
		"vsubpd			%%zmm7 ,%%zmm2,%%zmm2	\n\t		vsubpd		%%zmm15,%%zmm10,%%zmm10		\n\t"\
		"vsubpd			%%zmm6 ,%%zmm3,%%zmm3	\n\t		vsubpd		%%zmm14,%%zmm11,%%zmm11		\n\t"\
		"vmovaps		%%zmm0 ,      (%%rbx)	\n\t		vmovaps		%%zmm8 , 0x100(%%rbx)		\n\t"\
		"vmovaps		%%zmm1 , 0x040(%%rbx)	\n\t		vmovaps		%%zmm9 , 0x140(%%rbx)		\n\t"\
		"vmovaps		%%zmm2 ,      (%%rdx)	\n\t		vmovaps		%%zmm10, 0x100(%%rdx)		\n\t"\
		"vmovaps		%%zmm3 , 0x040(%%rcx)	\n\t		vmovaps		%%zmm11, 0x140(%%rcx)		\n\t"\
	"vmovaps %%zmm14,0x140(%%rdx)\n\t"/* spill zmm14 to make room for two */"vmovaps (%%rdi),%%zmm14\n\t"\
	"vfmadd132pd	%%zmm14,%%zmm0,%%zmm4		\n\t	vfmadd132pd		%%zmm14,%%zmm8 ,%%zmm12		\n\t"\
	"vfmadd132pd	%%zmm14,%%zmm1,%%zmm5		\n\t	vfmadd132pd		%%zmm14,%%zmm9 ,%%zmm13		\n\t"\
	"vfmadd132pd	%%zmm14,%%zmm2,%%zmm7		\n\t	vfmadd132pd		%%zmm14,%%zmm10,%%zmm15		\n\t"\
	"vfmadd132pd	%%zmm14,%%zmm3,%%zmm6		\n\t	vfmadd132pd	0x140(%%rdx),%%zmm11,%%zmm14	\n\t"\
		"vmovaps	%%zmm4 ,    (%%rax)			\n\t		vmovaps		%%zmm12,0x100(%%rax)		\n\t"\
		"vmovaps	%%zmm5 ,0x40(%%rax)			\n\t		vmovaps		%%zmm13,0x140(%%rax)		\n\t"\
		"vmovaps	%%zmm7 ,    (%%rcx)			\n\t		vmovaps		%%zmm15,0x100(%%rcx)		\n\t"\
		"vmovaps	%%zmm6 ,0x40(%%rdx)			\n\t		vmovaps		%%zmm14,0x140(%%rdx)		\n\t"\
		"addq		$0x200,%%rax			\n\t"\
		"addq		$0x200,%%rbx			\n\t"\
		"addq		$0x200,%%rcx			\n\t"\
		"addq		$0x200,%%rdx			\n\t		vmovaps		 (%%rdi),%%zmm13	\n\t"/* 2.0*/\
		"vmovaps		(%%rbx),%%zmm0			\n\t		vmovaps		0x100(%%rbx),%%zmm8 		\n\t"\
		"vmovaps	0x40(%%rbx),%%zmm1			\n\t		vmovaps		0x140(%%rbx),%%zmm9 		\n\t"\
		"vmovaps		(%%rax),%%zmm2			\n\t		vmovaps		0x100(%%rax),%%zmm10		\n\t"\
		"vmovaps	0x40(%%rax),%%zmm3			\n\t		vmovaps		0x140(%%rax),%%zmm11		\n\t"\
		"vmovaps		(%%rdx),%%zmm4			\n\t		vmovaps		0x100(%%rdx),%%zmm12		\n\t"\
		"vmovaps	0x40(%%rdx),%%zmm5			\n\t	/*	vmovaps		0x140(%%rdx),%%zmm13	*/	\n\t"\
		"vmovaps		(%%rcx),%%zmm6			\n\t		vmovaps		0x100(%%rcx),%%zmm14		\n\t"\
		"vmovaps	0x40(%%rcx),%%zmm7			\n\t		vmovaps		0x140(%%rcx),%%zmm15		\n\t"\
		"vsubpd			%%zmm0,%%zmm2,%%zmm2	\n\t		vsubpd		%%zmm8 ,%%zmm10,%%zmm10		\n\t"\
		"vsubpd			%%zmm1,%%zmm3,%%zmm3	\n\t		vsubpd		%%zmm9 ,%%zmm11,%%zmm11		\n\t"\
		"vsubpd			%%zmm4,%%zmm6,%%zmm6	\n\t		vsubpd		%%zmm12,%%zmm14,%%zmm14		\n\t"\
		"vsubpd			%%zmm5,%%zmm7,%%zmm7	\n\t		vsubpd 0x140(%%rdx),%%zmm15,%%zmm15		\n\t"\
		"vfmadd132pd	%%zmm13,%%zmm2,%%zmm0	\n\t	 vfmadd132pd	%%zmm13,%%zmm10,%%zmm8 		\n\t"\
		"vfmadd132pd	%%zmm13,%%zmm3,%%zmm1	\n\t	 vfmadd132pd	%%zmm13,%%zmm11,%%zmm9 		\n\t"\
		"vfmadd132pd	%%zmm13,%%zmm6,%%zmm4	\n\t	 vfmadd132pd	%%zmm13,%%zmm14,%%zmm12		\n\t"\
		"vfmadd132pd	%%zmm13,%%zmm7,%%zmm5	\n\t	 vfmadd132pd 0x140(%%rdx),%%zmm15,%%zmm13	\n\t"\
		"vsubpd			%%zmm4 ,%%zmm0,%%zmm0	\n\t		vsubpd		%%zmm12,%%zmm8 ,%%zmm8 		\n\t"\
		"vsubpd			%%zmm5 ,%%zmm1,%%zmm1	\n\t		vsubpd		%%zmm13,%%zmm9 ,%%zmm9 		\n\t"\
		"vsubpd			%%zmm7 ,%%zmm2,%%zmm2	\n\t		vsubpd		%%zmm15,%%zmm10,%%zmm10		\n\t"\
		"vsubpd			%%zmm6 ,%%zmm3,%%zmm3	\n\t		vsubpd		%%zmm14,%%zmm11,%%zmm11		\n\t"\
	"vmovaps %%zmm14,0x140(%%rdx)\n\t"/* spill zmm14 to make room for two */"vmovaps (%%rdi),%%zmm14\n\t"\
	"vfmadd132pd	%%zmm14,%%zmm0,%%zmm4		\n\t	vfmadd132pd		%%zmm14,%%zmm8 ,%%zmm12		\n\t"\
	"vfmadd132pd	%%zmm14,%%zmm1,%%zmm5		\n\t	vfmadd132pd		%%zmm14,%%zmm9 ,%%zmm13		\n\t"\
	"vfmadd132pd	%%zmm14,%%zmm2,%%zmm7		\n\t	vfmadd132pd		%%zmm14,%%zmm10,%%zmm15		\n\t"\
	"vfmadd132pd	%%zmm14,%%zmm3,%%zmm6		\n\t	vfmadd132pd	0x140(%%rdx),%%zmm11,%%zmm14	\n\t"\
		"vmovaps	%%zmm0 ,      (%%rbx)		\n\t		vmovaps			%%zmm8 , 0x100(%%rbx)	\n\t"\
		"vmovaps	%%zmm1 , 0x040(%%rbx)		\n\t		vmovaps			%%zmm9 , 0x140(%%rbx)	\n\t"\
		"vmovaps	%%zmm4 ,      (%%rax)		\n\t		vmovaps			%%zmm12, 0x100(%%rax)	\n\t"\
		"vmovaps	%%zmm5 , 0x040(%%rax)		\n\t		vmovaps			%%zmm13, 0x140(%%rax)	\n\t"\
		"vaddpd		%%zmm7 ,%%zmm3,%%zmm0		\n\t		vaddpd			%%zmm15,%%zmm11,%%zmm8 	\n\t"\
		"vaddpd		%%zmm2 ,%%zmm6,%%zmm1		\n\t		vaddpd			%%zmm10,%%zmm14,%%zmm9 	\n\t"\
		"vsubpd		%%zmm7 ,%%zmm3,%%zmm3		\n\t		vsubpd			%%zmm15,%%zmm11,%%zmm11	\n\t"\
		"vsubpd		%%zmm2 ,%%zmm6,%%zmm6		\n\t		vsubpd			%%zmm10,%%zmm14,%%zmm14	\n\t"\
	"vmovaps	(%%rsi),%%zmm15	\n\t"/* isrt2 */\
		"vmulpd		%%zmm15,%%zmm3,%%zmm3		\n\t		vmulpd			%%zmm15,%%zmm11,%%zmm11	\n\t"\
		"vmulpd		%%zmm15,%%zmm6,%%zmm6		\n\t		vmulpd			%%zmm15,%%zmm14,%%zmm14	\n\t"\
		"vmulpd		%%zmm15,%%zmm0,%%zmm0		\n\t		vmulpd			%%zmm15,%%zmm8 ,%%zmm8 	\n\t"\
		"vmulpd		%%zmm15,%%zmm1,%%zmm1		\n\t		vmulpd			%%zmm15,%%zmm9 ,%%zmm9 	\n\t"\
		"vmovaps	%%zmm3 , 0x040(%%rcx)		\n\t		vmovaps			%%zmm11, 0x140(%%rcx)	\n\t"\
		"vmovaps	%%zmm6 , 0x040(%%rdx)		\n\t		vmovaps			%%zmm14, 0x140(%%rdx)	\n\t"\
		"vmovaps	%%zmm0 ,      (%%rcx)		\n\t		vmovaps			%%zmm8 , 0x100(%%rcx)	\n\t"\
		"vmovaps	%%zmm1 ,      (%%rdx)		\n\t		vmovaps			%%zmm9 , 0x100(%%rdx)	\n\t"\
	/***** SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS *****/		/***** SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS *****/\
	/*        (r00,r10,r20,r30,r08,r18,r28,r38): */		/*        (r04,r14,r24,r34,r0C,r1C,r2C,r3C): */\
		"vmovaps	-0x200(%%rax),%%zmm0		\n\t		vmovaps		-0x100(%%rax),%%zmm8 		\n\t"\
		"vmovaps	-0x200(%%rbx),%%zmm4		\n\t		vmovaps		-0x100(%%rbx),%%zmm12		\n\t"\
		"vmovaps	-0x1c0(%%rax),%%zmm1		\n\t		vmovaps		-0x0c0(%%rax),%%zmm9 		\n\t"\
		"vmovaps	-0x1c0(%%rbx),%%zmm5		\n\t		vmovaps		-0x0c0(%%rbx),%%zmm13		\n\t"\
		"vmovaps		  (%%rax),%%zmm2		\n\t		vmovaps		 0x100(%%rax),%%zmm10		\n\t"\
		"vmovaps	 0x040(%%rbx),%%zmm7		\n\t		vmovaps		 0x140(%%rbx),%%zmm15		\n\t"\
		"vmovaps	 0x040(%%rax),%%zmm3		\n\t		vmovaps		 0x140(%%rax),%%zmm11		\n\t"\
		"vmovaps		  (%%rbx),%%zmm6		\n\t		vmovaps		 0x100(%%rbx),%%zmm14		\n\t"\
		"vsubpd		%%zmm2 ,%%zmm0,%%zmm0		\n\t		vsubpd		%%zmm10,%%zmm8 ,%%zmm8 		\n\t"\
		"vsubpd		%%zmm7 ,%%zmm4,%%zmm4		\n\t		vsubpd		%%zmm15,%%zmm12,%%zmm12		\n\t"\
		"vsubpd		%%zmm3 ,%%zmm1,%%zmm1		\n\t		vsubpd		%%zmm11,%%zmm9 ,%%zmm9 		\n\t"\
		"vsubpd		%%zmm6 ,%%zmm5,%%zmm5		\n\t		vsubpd		%%zmm14,%%zmm13,%%zmm13		\n\t"\
	"vmovaps %%zmm14,0x140(%%rbx)\n\t"/* spill zmm14 to make room for two */"vmovaps (%%rdi),%%zmm14\n\t"\
	"vfmadd132pd	%%zmm14,%%zmm0,%%zmm2		\n\t	vfmadd132pd		%%zmm14,%%zmm8 ,%%zmm10		\n\t"\
	"vfmadd132pd	%%zmm14,%%zmm4,%%zmm7		\n\t	vfmadd132pd		%%zmm14,%%zmm12,%%zmm15		\n\t"\
	"vfmadd132pd	%%zmm14,%%zmm1,%%zmm3		\n\t	vfmadd132pd		%%zmm14,%%zmm9 ,%%zmm11		\n\t"\
	"vfmadd132pd	%%zmm14,%%zmm5,%%zmm6		\n\t	vfmadd132pd	0x140(%%rbx),%%zmm13,%%zmm14	\n\t"\
		"vmovaps	%%zmm0,      (%%rax)		\n\t		vmovaps		%%zmm8 , 0x100(%%rax)		\n\t"\
		"vmovaps	%%zmm4,      (%%rbx)		\n\t		vmovaps		%%zmm12, 0x100(%%rbx)		\n\t"\
		"vmovaps	%%zmm1, 0x040(%%rax)		\n\t		vmovaps		%%zmm9 , 0x140(%%rax)		\n\t"\
		"vmovaps	%%zmm5,-0x1c0(%%rbx)		\n\t		vmovaps		%%zmm13,-0x0c0(%%rbx)		\n\t"\
		"vmovaps	%%zmm2,-0x200(%%rax)		\n\t		vmovaps		%%zmm10,-0x100(%%rax)		\n\t"\
		"vmovaps	%%zmm7,-0x200(%%rbx)		\n\t		vmovaps		%%zmm15,-0x100(%%rbx)		\n\t"\
		"vmovaps	%%zmm3,-0x1c0(%%rax)		\n\t		vmovaps		%%zmm11,-0x0c0(%%rax)		\n\t"\
		"vmovaps	%%zmm6, 0x040(%%rbx)		\n\t		vmovaps		%%zmm14, 0x140(%%rbx)		\n\t"\
		"vmovaps	-0x200(%%rcx),%%zmm0		\n\t		vmovaps		-0x100(%%rcx),%%zmm8 		\n\t"\
		"vmovaps	-0x200(%%rdx),%%zmm4		\n\t		vmovaps		-0x100(%%rdx),%%zmm12		\n\t"\
		"vmovaps	-0x1c0(%%rcx),%%zmm1		\n\t		vmovaps		-0x0c0(%%rcx),%%zmm9 		\n\t"\
		"vmovaps	-0x1c0(%%rdx),%%zmm5		\n\t		vmovaps		-0x0c0(%%rdx),%%zmm13		\n\t"\
		"vmovaps		  (%%rcx),%%zmm2		\n\t		vmovaps		 0x100(%%rcx),%%zmm10		\n\t"\
		"vmovaps	 0x040(%%rdx),%%zmm7		\n\t		vmovaps		 0x140(%%rdx),%%zmm15		\n\t"\
		"vmovaps	 0x040(%%rcx),%%zmm3		\n\t		vmovaps		 0x140(%%rcx),%%zmm11		\n\t"\
		"vmovaps		  (%%rdx),%%zmm6		\n\t		vmovaps		 0x100(%%rdx),%%zmm14		\n\t"\
		"vsubpd		%%zmm2 ,%%zmm0,%%zmm0		\n\t		vsubpd		%%zmm10,%%zmm8 ,%%zmm8 		\n\t"\
		"vsubpd		%%zmm7 ,%%zmm4,%%zmm4		\n\t		vsubpd		%%zmm15,%%zmm12,%%zmm12		\n\t"\
		"vsubpd		%%zmm3 ,%%zmm1,%%zmm1		\n\t		vsubpd		%%zmm11,%%zmm9 ,%%zmm9 		\n\t"\
		"vsubpd		%%zmm6 ,%%zmm5,%%zmm5		\n\t		vsubpd		%%zmm14,%%zmm13,%%zmm13		\n\t"\
	"vmovaps %%zmm14,0x140(%%rdx)\n\t"/* spill zmm14 to make room for two */"vmovaps (%%rdi),%%zmm14\n\t"\
	"vfmadd132pd	%%zmm14,%%zmm0,%%zmm2		\n\t	vfmadd132pd		%%zmm14,%%zmm8 ,%%zmm10		\n\t"\
	"vfmadd132pd	%%zmm14,%%zmm4,%%zmm7		\n\t	vfmadd132pd		%%zmm14,%%zmm12,%%zmm15		\n\t"\
	"vfmadd132pd	%%zmm14,%%zmm1,%%zmm3		\n\t	vfmadd132pd		%%zmm14,%%zmm9 ,%%zmm11		\n\t"\
	"vfmadd132pd	%%zmm14,%%zmm5,%%zmm6		\n\t	vfmadd132pd	0x140(%%rdx),%%zmm13,%%zmm14	\n\t"\
		"vmovaps	%%zmm0,      (%%rcx)		\n\t		vmovaps		%%zmm8 , 0x100(%%rcx)		\n\t"\
		"vmovaps	%%zmm4,      (%%rdx)		\n\t		vmovaps		%%zmm12, 0x100(%%rdx)		\n\t"\
		"vmovaps	%%zmm1, 0x040(%%rcx)		\n\t		vmovaps		%%zmm9 , 0x140(%%rcx)		\n\t"\
		"vmovaps	%%zmm5,-0x1c0(%%rdx)		\n\t		vmovaps		%%zmm13,-0x0c0(%%rdx)		\n\t"\
		"vmovaps	%%zmm2,-0x200(%%rcx)		\n\t		vmovaps		%%zmm10,-0x100(%%rcx)		\n\t"\
		"vmovaps	%%zmm7,-0x200(%%rdx)		\n\t		vmovaps		%%zmm15,-0x100(%%rdx)		\n\t"\
		"vmovaps	%%zmm3,-0x1c0(%%rcx)		\n\t		vmovaps		%%zmm11,-0x0c0(%%rcx)		\n\t"\
		"vmovaps	%%zmm6, 0x040(%%rdx)		\n\t		vmovaps		%%zmm14, 0x140(%%rdx)		\n\t"\
	/* Blocks 3,4 have tmp-addresses offset +0x80 w.r.to Blocks 1,2, respectively (thus +0x200-0x180 = +0x080: */\
		"subq		$0x180,%%rax			\n\t"\
		"subq		$0x180,%%rbx			\n\t"\
		"subq		$0x180,%%rcx			\n\t"\
		"subq		$0x180,%%rdx			\n\t		vmovaps		 (%%rdi),%%zmm13	\n\t"/* 2.0*/\
		"vmovaps		(%%rbx),%%zmm0			\n\t		vmovaps		0x100(%%rbx),%%zmm8 		\n\t"\
		"vmovaps	0x40(%%rbx),%%zmm1			\n\t		vmovaps		0x140(%%rbx),%%zmm9 		\n\t"\
		"vmovaps		(%%rax),%%zmm2			\n\t		vmovaps		0x100(%%rax),%%zmm10		\n\t"\
		"vmovaps	0x40(%%rax),%%zmm3			\n\t		vmovaps		0x140(%%rax),%%zmm11		\n\t"\
		"vmovaps		(%%rdx),%%zmm4			\n\t		vmovaps		0x100(%%rdx),%%zmm12		\n\t"\
		"vmovaps	0x40(%%rdx),%%zmm5			\n\t	/*	vmovaps		0x140(%%rdx),%%zmm13	*/	\n\t"\
		"vmovaps		(%%rcx),%%zmm6			\n\t		vmovaps		0x100(%%rcx),%%zmm14		\n\t"\
		"vmovaps	0x40(%%rcx),%%zmm7			\n\t		vmovaps		0x140(%%rcx),%%zmm15		\n\t"\
		"vsubpd			%%zmm0,%%zmm2,%%zmm2	\n\t		vsubpd		%%zmm8 ,%%zmm10,%%zmm10		\n\t"\
		"vsubpd			%%zmm1,%%zmm3,%%zmm3	\n\t		vsubpd		%%zmm9 ,%%zmm11,%%zmm11		\n\t"\
		"vsubpd			%%zmm4,%%zmm6,%%zmm6	\n\t		vsubpd		%%zmm12,%%zmm14,%%zmm14		\n\t"\
		"vsubpd			%%zmm5,%%zmm7,%%zmm7	\n\t		vsubpd 0x140(%%rdx),%%zmm15,%%zmm15		\n\t"\
		"vfmadd132pd	%%zmm13,%%zmm2,%%zmm0	\n\t	 vfmadd132pd	%%zmm13,%%zmm10,%%zmm8 		\n\t"\
		"vfmadd132pd	%%zmm13,%%zmm3,%%zmm1	\n\t	 vfmadd132pd	%%zmm13,%%zmm11,%%zmm9 		\n\t"\
		"vfmadd132pd	%%zmm13,%%zmm6,%%zmm4	\n\t	 vfmadd132pd	%%zmm13,%%zmm14,%%zmm12		\n\t"\
		"vfmadd132pd	%%zmm13,%%zmm7,%%zmm5	\n\t	 vfmadd132pd 0x140(%%rdx),%%zmm15,%%zmm13	\n\t"\
		"vsubpd			%%zmm4 ,%%zmm0,%%zmm0	\n\t		vsubpd			%%zmm12,%%zmm8 ,%%zmm8 	\n\t"\
		"vsubpd			%%zmm5 ,%%zmm1,%%zmm1	\n\t		vsubpd			%%zmm13,%%zmm9 ,%%zmm9 	\n\t"\
		"vsubpd			%%zmm7 ,%%zmm2,%%zmm2	\n\t		vsubpd			%%zmm15,%%zmm10,%%zmm10	\n\t"\
		"vsubpd			%%zmm6 ,%%zmm3,%%zmm3	\n\t		vsubpd			%%zmm14,%%zmm11,%%zmm11	\n\t"\
		"vmovaps		%%zmm0 ,      (%%rbx)	\n\t		vmovaps			%%zmm8 , 0x100(%%rbx)	\n\t"\
		"vmovaps		%%zmm1 , 0x040(%%rbx)	\n\t		vmovaps			%%zmm9 , 0x140(%%rbx)	\n\t"\
		"vmovaps		%%zmm2 ,      (%%rdx)	\n\t		vmovaps			%%zmm10, 0x100(%%rdx)	\n\t"\
		"vmovaps		%%zmm3 , 0x040(%%rcx)	\n\t		vmovaps			%%zmm11, 0x140(%%rcx)	\n\t"\
	"vmovaps %%zmm14,0x140(%%rdx)\n\t"/* spill zmm14 to make room for two */"vmovaps (%%rdi),%%zmm14\n\t"\
	"vfmadd132pd	%%zmm14,%%zmm0,%%zmm4		\n\t	vfmadd132pd		%%zmm14,%%zmm8 ,%%zmm12		\n\t"\
	"vfmadd132pd	%%zmm14,%%zmm1,%%zmm5		\n\t	vfmadd132pd		%%zmm14,%%zmm9 ,%%zmm13		\n\t"\
	"vfmadd132pd	%%zmm14,%%zmm2,%%zmm7		\n\t	vfmadd132pd		%%zmm14,%%zmm10,%%zmm15		\n\t"\
	"vfmadd132pd	%%zmm14,%%zmm3,%%zmm6		\n\t	vfmadd132pd	0x140(%%rdx),%%zmm11,%%zmm14		\n\t"\
		"vmovaps	%%zmm4 ,    (%%rax)			\n\t		vmovaps		%%zmm12,0x100(%%rax)			\n\t"\
		"vmovaps	%%zmm5 ,0x40(%%rax)			\n\t		vmovaps		%%zmm13,0x140(%%rax)			\n\t"\
		"vmovaps	%%zmm7 ,    (%%rcx)			\n\t		vmovaps		%%zmm15,0x100(%%rcx)			\n\t"\
		"vmovaps	%%zmm6 ,0x40(%%rdx)			\n\t		vmovaps		%%zmm14,0x140(%%rdx)			\n\t"\
		"addq		$0x200,%%rax			\n\t"\
		"addq		$0x200,%%rbx			\n\t"\
		"addq		$0x200,%%rcx			\n\t"\
		"addq		$0x200,%%rdx			\n\t		vmovaps		 (%%rdi),%%zmm13	\n\t"/* 2.0*/\
		"vmovaps		(%%rbx),%%zmm0			\n\t		vmovaps		0x100(%%rbx),%%zmm8 		\n\t"\
		"vmovaps	0x40(%%rbx),%%zmm1			\n\t		vmovaps		0x140(%%rbx),%%zmm9 		\n\t"\
		"vmovaps		(%%rax),%%zmm2			\n\t		vmovaps		0x100(%%rax),%%zmm10		\n\t"\
		"vmovaps	0x40(%%rax),%%zmm3			\n\t		vmovaps		0x140(%%rax),%%zmm11		\n\t"\
		"vmovaps		(%%rdx),%%zmm4			\n\t		vmovaps		0x100(%%rdx),%%zmm12		\n\t"\
		"vmovaps	0x40(%%rdx),%%zmm5			\n\t	/*	vmovaps		0x140(%%rdx),%%zmm13	*/	\n\t"\
		"vmovaps		(%%rcx),%%zmm6			\n\t		vmovaps		0x100(%%rcx),%%zmm14		\n\t"\
		"vmovaps	0x40(%%rcx),%%zmm7			\n\t		vmovaps		0x140(%%rcx),%%zmm15		\n\t"\
		"vsubpd			%%zmm0,%%zmm2,%%zmm2	\n\t		vsubpd		%%zmm8 ,%%zmm10,%%zmm10		\n\t"\
		"vsubpd			%%zmm1,%%zmm3,%%zmm3	\n\t		vsubpd		%%zmm9 ,%%zmm11,%%zmm11		\n\t"\
		"vsubpd			%%zmm4,%%zmm6,%%zmm6	\n\t		vsubpd		%%zmm12,%%zmm14,%%zmm14		\n\t"\
		"vsubpd			%%zmm5,%%zmm7,%%zmm7	\n\t		vsubpd 0x140(%%rdx),%%zmm15,%%zmm15		\n\t"\
		"vfmadd132pd	%%zmm13,%%zmm2,%%zmm0	\n\t	 vfmadd132pd	%%zmm13,%%zmm10,%%zmm8 		\n\t"\
		"vfmadd132pd	%%zmm13,%%zmm3,%%zmm1	\n\t	 vfmadd132pd	%%zmm13,%%zmm11,%%zmm9 		\n\t"\
		"vfmadd132pd	%%zmm13,%%zmm6,%%zmm4	\n\t	 vfmadd132pd	%%zmm13,%%zmm14,%%zmm12		\n\t"\
		"vfmadd132pd	%%zmm13,%%zmm7,%%zmm5	\n\t	 vfmadd132pd 0x140(%%rdx),%%zmm15,%%zmm13	\n\t"\
		"vsubpd			%%zmm4 ,%%zmm0,%%zmm0	\n\t		vsubpd		%%zmm12,%%zmm8 ,%%zmm8 		\n\t"\
		"vsubpd			%%zmm5 ,%%zmm1,%%zmm1	\n\t		vsubpd		%%zmm13,%%zmm9 ,%%zmm9 		\n\t"\
		"vsubpd			%%zmm7 ,%%zmm2,%%zmm2	\n\t		vsubpd		%%zmm15,%%zmm10,%%zmm10		\n\t"\
		"vsubpd			%%zmm6 ,%%zmm3,%%zmm3	\n\t		vsubpd		%%zmm14,%%zmm11,%%zmm11		\n\t"\
	"vmovaps %%zmm14,0x140(%%rdx)\n\t"/* spill zmm14 to make room for two */"vmovaps (%%rdi),%%zmm14\n\t"\
	"vfmadd132pd	%%zmm14,%%zmm0,%%zmm4		\n\t	vfmadd132pd		%%zmm14,%%zmm8 ,%%zmm12		\n\t"\
	"vfmadd132pd	%%zmm14,%%zmm1,%%zmm5		\n\t	vfmadd132pd		%%zmm14,%%zmm9 ,%%zmm13		\n\t"\
	"vfmadd132pd	%%zmm14,%%zmm2,%%zmm7		\n\t	vfmadd132pd		%%zmm14,%%zmm10,%%zmm15		\n\t"\
	"vfmadd132pd	%%zmm14,%%zmm3,%%zmm6		\n\t	vfmadd132pd	0x140(%%rdx),%%zmm11,%%zmm14	\n\t"\
		"vmovaps	%%zmm0 ,      (%%rbx)		\n\t		vmovaps			%%zmm8 , 0x100(%%rbx)	\n\t"\
		"vmovaps	%%zmm1 , 0x040(%%rbx)		\n\t		vmovaps			%%zmm9 , 0x140(%%rbx)	\n\t"\
		"vmovaps	%%zmm4 ,      (%%rax)		\n\t		vmovaps			%%zmm12, 0x100(%%rax)	\n\t"\
		"vmovaps	%%zmm5 , 0x040(%%rax)		\n\t		vmovaps			%%zmm13, 0x140(%%rax)	\n\t"\
		"vaddpd		%%zmm7 ,%%zmm3,%%zmm0		\n\t		vaddpd			%%zmm15,%%zmm11,%%zmm8 	\n\t"\
		"vaddpd		%%zmm2 ,%%zmm6,%%zmm1		\n\t		vaddpd			%%zmm10,%%zmm14,%%zmm9 	\n\t"\
		"vsubpd		%%zmm7 ,%%zmm3,%%zmm3		\n\t		vsubpd			%%zmm15,%%zmm11,%%zmm11	\n\t"\
		"vsubpd		%%zmm2 ,%%zmm6,%%zmm6		\n\t		vsubpd			%%zmm10,%%zmm14,%%zmm14	\n\t"\
	"vmovaps	(%%rsi),%%zmm15	\n\t"/* isrt2 */\
		"vmulpd		%%zmm15,%%zmm3,%%zmm3		\n\t		vmulpd			%%zmm15,%%zmm11,%%zmm11	\n\t"\
		"vmulpd		%%zmm15,%%zmm6,%%zmm6		\n\t		vmulpd			%%zmm15,%%zmm14,%%zmm14	\n\t"\
		"vmulpd		%%zmm15,%%zmm0,%%zmm0		\n\t		vmulpd			%%zmm15,%%zmm8 ,%%zmm8 	\n\t"\
		"vmulpd		%%zmm15,%%zmm1,%%zmm1		\n\t		vmulpd			%%zmm15,%%zmm9 ,%%zmm9 	\n\t"\
		"vmovaps	%%zmm3 , 0x040(%%rcx)		\n\t		vmovaps			%%zmm11, 0x140(%%rcx)	\n\t"\
		"vmovaps	%%zmm6 , 0x040(%%rdx)		\n\t		vmovaps			%%zmm14, 0x140(%%rdx)	\n\t"\
		"vmovaps	%%zmm0 ,      (%%rcx)		\n\t		vmovaps			%%zmm8 , 0x100(%%rcx)	\n\t"\
		"vmovaps	%%zmm1 ,      (%%rdx)		\n\t		vmovaps			%%zmm9 , 0x100(%%rdx)	\n\t"\
	/***** SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS *****/		/***** SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS *****/\
	/*        (r02,r12,r22,r32,r0A,r1A,r2A,r3A): */		/*        (r06,r16,r26,r36,r0E,r1E,r2E,r3E): */\
		"vmovaps	-0x200(%%rax),%%zmm0		\n\t		vmovaps		-0x100(%%rax),%%zmm8 		\n\t"\
		"vmovaps	-0x200(%%rbx),%%zmm4		\n\t		vmovaps		-0x100(%%rbx),%%zmm12		\n\t"\
		"vmovaps	-0x1c0(%%rax),%%zmm1		\n\t		vmovaps		-0x0c0(%%rax),%%zmm9 		\n\t"\
		"vmovaps	-0x1c0(%%rbx),%%zmm5		\n\t		vmovaps		-0x0c0(%%rbx),%%zmm13		\n\t"\
		"vmovaps		  (%%rax),%%zmm2		\n\t		vmovaps		 0x100(%%rax),%%zmm10		\n\t"\
		"vmovaps	 0x040(%%rbx),%%zmm7		\n\t		vmovaps		 0x140(%%rbx),%%zmm15		\n\t"\
		"vmovaps	 0x040(%%rax),%%zmm3		\n\t		vmovaps		 0x140(%%rax),%%zmm11		\n\t"\
		"vmovaps		  (%%rbx),%%zmm6		\n\t		vmovaps		 0x100(%%rbx),%%zmm14		\n\t"\
		"vsubpd		%%zmm2 ,%%zmm0,%%zmm0		\n\t		vsubpd		%%zmm10,%%zmm8 ,%%zmm8 		\n\t"\
		"vsubpd		%%zmm7 ,%%zmm4,%%zmm4		\n\t		vsubpd		%%zmm15,%%zmm12,%%zmm12		\n\t"\
		"vsubpd		%%zmm3 ,%%zmm1,%%zmm1		\n\t		vsubpd		%%zmm11,%%zmm9 ,%%zmm9 		\n\t"\
		"vsubpd		%%zmm6 ,%%zmm5,%%zmm5		\n\t		vsubpd		%%zmm14,%%zmm13,%%zmm13		\n\t"\
	"vmovaps %%zmm14,0x140(%%rbx)\n\t"/* spill zmm14 to make room for two */"vmovaps (%%rdi),%%zmm14\n\t"\
	"vfmadd132pd	%%zmm14,%%zmm0,%%zmm2		\n\t	vfmadd132pd		%%zmm14,%%zmm8 ,%%zmm10		\n\t"\
	"vfmadd132pd	%%zmm14,%%zmm4,%%zmm7		\n\t	vfmadd132pd		%%zmm14,%%zmm12,%%zmm15		\n\t"\
	"vfmadd132pd	%%zmm14,%%zmm1,%%zmm3		\n\t	vfmadd132pd		%%zmm14,%%zmm9 ,%%zmm11		\n\t"\
	"vfmadd132pd	%%zmm14,%%zmm5,%%zmm6		\n\t	vfmadd132pd	0x140(%%rbx),%%zmm13,%%zmm14	\n\t"\
		"vmovaps	%%zmm0,      (%%rax)		\n\t		vmovaps		%%zmm8 , 0x100(%%rax)		\n\t"\
		"vmovaps	%%zmm4,      (%%rbx)		\n\t		vmovaps		%%zmm12, 0x100(%%rbx)		\n\t"\
		"vmovaps	%%zmm1, 0x040(%%rax)		\n\t		vmovaps		%%zmm9 , 0x140(%%rax)		\n\t"\
		"vmovaps	%%zmm5,-0x1c0(%%rbx)		\n\t		vmovaps		%%zmm13,-0x0c0(%%rbx)		\n\t"\
		"vmovaps	%%zmm2,-0x200(%%rax)		\n\t		vmovaps		%%zmm10,-0x100(%%rax)		\n\t"\
		"vmovaps	%%zmm7,-0x200(%%rbx)		\n\t		vmovaps		%%zmm15,-0x100(%%rbx)		\n\t"\
		"vmovaps	%%zmm3,-0x1c0(%%rax)		\n\t		vmovaps		%%zmm11,-0x0c0(%%rax)		\n\t"\
		"vmovaps	%%zmm6, 0x040(%%rbx)		\n\t		vmovaps		%%zmm14, 0x140(%%rbx)		\n\t"\
		"vmovaps	-0x200(%%rcx),%%zmm0		\n\t		vmovaps		-0x100(%%rcx),%%zmm8 		\n\t"\
		"vmovaps	-0x200(%%rdx),%%zmm4		\n\t		vmovaps		-0x100(%%rdx),%%zmm12		\n\t"\
		"vmovaps	-0x1c0(%%rcx),%%zmm1		\n\t		vmovaps		-0x0c0(%%rcx),%%zmm9 		\n\t"\
		"vmovaps	-0x1c0(%%rdx),%%zmm5		\n\t		vmovaps		-0x0c0(%%rdx),%%zmm13		\n\t"\
		"vmovaps		  (%%rcx),%%zmm2		\n\t		vmovaps		 0x100(%%rcx),%%zmm10		\n\t"\
		"vmovaps	 0x040(%%rdx),%%zmm7		\n\t		vmovaps		 0x140(%%rdx),%%zmm15		\n\t"\
		"vmovaps	 0x040(%%rcx),%%zmm3		\n\t		vmovaps		 0x140(%%rcx),%%zmm11		\n\t"\
		"vmovaps		  (%%rdx),%%zmm6		\n\t		vmovaps		 0x100(%%rdx),%%zmm14		\n\t"\
		"vsubpd		%%zmm2 ,%%zmm0,%%zmm0		\n\t		vsubpd		%%zmm10,%%zmm8 ,%%zmm8 		\n\t"\
		"vsubpd		%%zmm7 ,%%zmm4,%%zmm4		\n\t		vsubpd		%%zmm15,%%zmm12,%%zmm12		\n\t"\
		"vsubpd		%%zmm3 ,%%zmm1,%%zmm1		\n\t		vsubpd		%%zmm11,%%zmm9 ,%%zmm9 		\n\t"\
		"vsubpd		%%zmm6 ,%%zmm5,%%zmm5		\n\t		vsubpd		%%zmm14,%%zmm13,%%zmm13		\n\t"\
	"vmovaps %%zmm14,0x140(%%rdx)\n\t"/* spill zmm14 to make room for two */"vmovaps (%%rdi),%%zmm14\n\t"\
	"vfmadd132pd	%%zmm14,%%zmm0,%%zmm2		\n\t	vfmadd132pd		%%zmm14,%%zmm8 ,%%zmm10		\n\t"\
	"vfmadd132pd	%%zmm14,%%zmm4,%%zmm7		\n\t	vfmadd132pd		%%zmm14,%%zmm12,%%zmm15		\n\t"\
	"vfmadd132pd	%%zmm14,%%zmm1,%%zmm3		\n\t	vfmadd132pd		%%zmm14,%%zmm9 ,%%zmm11		\n\t"\
	"vfmadd132pd	%%zmm14,%%zmm5,%%zmm6		\n\t	vfmadd132pd	0x140(%%rdx),%%zmm13,%%zmm14	\n\t"\
		"vmovaps	%%zmm0,      (%%rcx)		\n\t		vmovaps		%%zmm8 , 0x100(%%rcx)		\n\t"\
		"vmovaps	%%zmm4,      (%%rdx)		\n\t		vmovaps		%%zmm12, 0x100(%%rdx)		\n\t"\
		"vmovaps	%%zmm1, 0x040(%%rcx)		\n\t		vmovaps		%%zmm9 , 0x140(%%rcx)		\n\t"\
		"vmovaps	%%zmm5,-0x1c0(%%rdx)		\n\t		vmovaps		%%zmm13,-0x0c0(%%rdx)		\n\t"\
		"vmovaps	%%zmm2,-0x200(%%rcx)		\n\t		vmovaps		%%zmm10,-0x100(%%rcx)		\n\t"\
		"vmovaps	%%zmm7,-0x200(%%rdx)		\n\t		vmovaps		%%zmm15,-0x100(%%rdx)		\n\t"\
		"vmovaps	%%zmm3,-0x1c0(%%rcx)		\n\t		vmovaps		%%zmm11,-0x0c0(%%rcx)		\n\t"\
		"vmovaps	%%zmm6, 0x040(%%rdx)		\n\t		vmovaps		%%zmm14, 0x140(%%rdx)		\n\t"\
		"\n\t"\
	/***************************************************************************************/\
	/* Now do eight more radix-4 transforms, including the internal and external twiddles: */\
	/***************************************************************************************/\
		"\n\t"\
	/*
	Using upper block(s) of main array for temp-storage in section below led to a nasty AVX bug to track down:
	In fermat-mod mode, 4 block addresses in ascending order are add0,1,2,3 with no 'gaps' between blocks; for
	mersenne-mod addresses in asc. order are add0,2,3,1 with gaps between contiguous-data-block pairs 0,2 & 3,1.
	Thus for ferm-mod need [add2] as base-address of 'high-half' block for temp-storage; for mers-mod need [add3].
	In both cases have (add2 < add3) so instead use (add2 - add1) to differentiate: > 0 for fermat, < 0 for mersenne.
	*/\
		"movq	%[__add2],%%rsi		\n\t"/* destroyable copy of add2 */\
		"movq	%[__add2],%%rbx		\n\t"\
		"subq	%[__add1],%%rsi		\n\t"/* rsi = (add2 - add1); if this yields a borrow (i.e. sets CF) it's mersenne, else fermat. */\
		"cmovcq %[__add3],%%rbx	\n\t" /* if CF set (i.e. h > l), copy source [add3] into dest (rbx), else leave dest = [add2]. */\
		"\n\t"\
	/********************************************/	/********************************************/\
	/* Block 2: t02,t12,t22,t32 -> r10,14,12,16 */	/* Block 6: t0A,t1A,t2A,t3A -> r18,1C,1A,1E */\
	/********************************************/	/********************************************/\
		"movq	%[__isrt2],%%rsi				\n\t"\
		"movq	%[__r10],%%rax	/* base-addr in rcol = c05/r18, so rax/r10 offset +0x200 vs lcol */\n\t"\
		"movq	%%rsi,%%rcx					\n\t"\
		"movq	%%rsi,%%rdx					\n\t"\
		"movq	%[__c01],%%r10					\n\t"\
		"addq	$0x040,%%rcx	/* cc0 */		\n\t"\
		"addq	$0x0c0,%%rdx	/* cc1 */		\n\t"\
		"vmovaps		0x080(%%rax),%%zmm4		\n\t		vmovaps		 0x280(%%rax),%%zmm12		\n\t"\
		"vmovaps		0x180(%%rax),%%zmm0		\n\t		vmovaps		 0x380(%%rax),%%zmm8 		\n\t"\
		"vmovaps		0x0c0(%%rax),%%zmm5		\n\t		vmovaps		 0x2c0(%%rax),%%zmm15		\n\t"\
		"vmovaps		0x1c0(%%rax),%%zmm3		\n\t"	/*	vmovaps		 0x3c0(%%rax),%%zmm9 		\n\t"*/\
	/* Using ensuing MULs to effect copying is a no-go here due to I*[...] 'twisting' - same for all 8 blocks: */\
		"vmovaps			%%zmm4 ,%%zmm6		\n\t		vmovaps		 	%%zmm12,%%zmm14			\n\t"\
		"vmovaps			%%zmm0 ,%%zmm2		\n\t"	/*	vmovaps		 	%%zmm8 ,%%zmm10			\n\t"*/\
		"vmovaps			%%zmm5 ,%%zmm7		\n\t"	/*	vmovaps		 	%%zmm13,%%zmm15			\n\t"*/\
	/*	"vmovaps			%%zmm1 ,%%zmm3	*/		"		vmovaps		 	%%zmm9 ,%%zmm11			\n\t"\
	/* Strategic register-usage for sincos consts here cuts loads from 16 to 4: */\
		"vmovaps			 (%%rdx),%%zmm9 	\n\t		vmovaps		0x040(%%rdx),%%zmm10		\n\t"\
		"vmovaps		0x080(%%rdx),%%zmm1 	\n\t		vmovaps		0x0c0(%%rdx),%%zmm13		\n\t"\
		"vmulpd			%%zmm10,%%zmm7,%%zmm7	\n\t		vmulpd			 %%zmm1 ,%%zmm15,%%zmm15\n\t"\
		"vmulpd			%%zmm13,%%zmm3,%%zmm3	\n\t		vmulpd			 %%zmm10,%%zmm11,%%zmm11\n\t"\
		"vmulpd			%%zmm10,%%zmm6,%%zmm6	\n\t		vmulpd			 %%zmm1 ,%%zmm14,%%zmm14\n\t"\
		"vmulpd			%%zmm13,%%zmm2,%%zmm2	\n\t		vmulpd		0x380(%%rax),%%zmm10,%%zmm10\n\t"\
	"vfmadd132pd		%%zmm9 ,%%zmm7,%%zmm4	\n\t	vfmadd132pd	 		 %%zmm13,%%zmm15,%%zmm12\n\t"\
	"vfmadd132pd		%%zmm1 ,%%zmm3,%%zmm0	\n\t	vfmsub132pd	 		 %%zmm9 ,%%zmm11,%%zmm8 \n\t"\
	"vfmsub132pd		%%zmm9 ,%%zmm6,%%zmm5	\n\t	vfmsub132pd	 	0x2c0(%%rax),%%zmm14,%%zmm13\n\t"\
	"vfmsub132pd   0x1c0(%%rax),%%zmm2,%%zmm1	\n\t	vfmadd132pd	 	0x3c0(%%rax),%%zmm10,%%zmm9 \n\t"\
		"vmovaps		%%zmm5,%%zmm7			\n\t		vmovaps		%%zmm13,%%zmm15				\n\t"\
		"vmovaps		%%zmm4,%%zmm6			\n\t		vmovaps		%%zmm12,%%zmm14				\n\t"\
		"vsubpd			%%zmm0,%%zmm4,%%zmm4	\n\t		vaddpd		%%zmm8 ,%%zmm12,%%zmm12		\n\t"\
		"vsubpd			%%zmm1,%%zmm5,%%zmm5	\n\t		vaddpd		%%zmm9 ,%%zmm13,%%zmm13		\n\t"\
		"vaddpd			%%zmm0,%%zmm6,%%zmm6	\n\t		vsubpd		%%zmm8 ,%%zmm14,%%zmm14		\n\t"\
		"vaddpd			%%zmm1,%%zmm7,%%zmm7	\n\t		vsubpd		%%zmm9 ,%%zmm15,%%zmm15		\n\t"\
		"vmovaps		 0x100(%%rax),%%zmm0	\n\t		vmovaps		0x300(%%rax),%%zmm8 		\n\t"\
		"vmovaps		 0x140(%%rax),%%zmm1	\n\t		vmovaps		0x340(%%rax),%%zmm9 		\n\t"\
		"vmovaps		 0x100(%%rax),%%zmm2	\n\t		vmovaps		0x300(%%rax),%%zmm10		\n\t"\
		"vmovaps		 0x140(%%rax),%%zmm3	\n\t		vmovaps		0x340(%%rax),%%zmm11		\n\t"\
		"vmulpd		 0x040(%%rcx),%%zmm0,%%zmm0	\n\t		vmulpd		(%%rcx),%%zmm8 ,%%zmm8 	\n\t"\
		"vmulpd		 0x040(%%rcx),%%zmm1,%%zmm1	\n\t		vmulpd		(%%rcx),%%zmm9 ,%%zmm9 	\n\t"\
	"vfmsub132pd		  (%%rcx),%%zmm0,%%zmm3	\n\t	vfmadd132pd		0x040(%%rcx),%%zmm8 ,%%zmm11\n\t"\
	"vfmadd132pd		  (%%rcx),%%zmm1,%%zmm2	\n\t	vfmsub132pd		0x040(%%rcx),%%zmm9 ,%%zmm10\n\t"\
		"vmovaps			  (%%rax),%%zmm0	\n\t		vmovaps		0x200(%%rax),%%zmm8 		\n\t"\
		"vmovaps		 0x040(%%rax),%%zmm1	\n\t		vmovaps		0x240(%%rax),%%zmm9 		\n\t"\
		"vsubpd			%%zmm2,%%zmm0,%%zmm0	\n\t		vsubpd		%%zmm10,%%zmm8 ,%%zmm8 		\n\t"\
		"vsubpd			%%zmm3,%%zmm1,%%zmm1	\n\t		vsubpd		%%zmm11,%%zmm9 ,%%zmm9 		\n\t"\
	"vfmadd132pd		(%%rdi),%%zmm0,%%zmm2	\n\t	vfmadd132pd		(%%rdi),%%zmm8 ,%%zmm10		\n\t"\
	"vfmadd132pd		(%%rdi),%%zmm1,%%zmm3	\n\t	vfmadd132pd		(%%rdi),%%zmm9 ,%%zmm11		\n\t"\
		"vsubpd			%%zmm6,%%zmm2,%%zmm2	\n\t		vsubpd		%%zmm14,%%zmm8 ,%%zmm8 		\n\t"\
		"vsubpd			%%zmm7,%%zmm3,%%zmm3	\n\t		vsubpd		%%zmm15,%%zmm9 ,%%zmm9 		\n\t"\
	"vfmadd132pd		(%%rdi),%%zmm2,%%zmm6	\n\t	vfmadd132pd		(%%rdi),%%zmm8 ,%%zmm14		\n\t"\
	"vfmadd132pd		(%%rdi),%%zmm3,%%zmm7	\n\t	vfmadd132pd		(%%rdi),%%zmm9 ,%%zmm15		\n\t"\
		"vmovaps		%%zmm2, 0x080(%%rax)	\n\t		vmovaps		%%zmm8 , 0x280(%%rax)		\n\t"\
		"vmovaps		%%zmm3, 0x0c0(%%rax)	\n\t		vmovaps		%%zmm9 , 0x2c0(%%rax)		\n\t"\
		"vmovaps		%%zmm6,%%zmm2			\n\t		vmovaps		%%zmm14,%%zmm8 				\n\t"\
		"vmovaps		%%zmm7,%%zmm3			\n\t		vmovaps		%%zmm15,%%zmm9 				\n\t"\
		"vmulpd			  (%%r10),%%zmm6,%%zmm6	\n\t		vmulpd		0x200(%%r10),%%zmm14,%%zmm14\n\t"\
		"vmulpd			  (%%r10),%%zmm7,%%zmm7	\n\t		vmulpd		0x200(%%r10),%%zmm15,%%zmm15\n\t"\
	" vfmadd231pd	 0x040(%%r10),%%zmm3,%%zmm6	\n\t	 vfmadd231pd	0x240(%%r10),%%zmm9 ,%%zmm14\n\t"\
	"vfnmadd231pd	 0x040(%%r10),%%zmm2,%%zmm7	\n\t	vfnmadd231pd	0x240(%%r10),%%zmm8 ,%%zmm15\n\t"\
		"vmovaps		%%zmm7, 0x040(%%rbx)	\n\t		vmovaps		%%zmm15, 0x140(%%rbx)		\n\t"\
		"vmovaps		%%zmm6,      (%%rbx)	\n\t		vmovaps		%%zmm14, 0x100(%%rbx)		\n\t"\
		"vmovaps		 0x080(%%rax),%%zmm6	\n\t		vmovaps		0x280(%%rax),%%zmm14		\n\t"\
		"vmovaps		 0x0c0(%%rax),%%zmm7	\n\t		vmovaps		0x2c0(%%rax),%%zmm15		\n\t"\
		"vmovaps		%%zmm6,%%zmm2			\n\t		vmovaps		%%zmm14,%%zmm8 				\n\t"\
		"vmovaps		%%zmm7,%%zmm3			\n\t		vmovaps		%%zmm15,%%zmm9 				\n\t"\
		"vmulpd		 0x080(%%r10),%%zmm6,%%zmm6	\n\t		vmulpd		0x280(%%r10),%%zmm14,%%zmm14\n\t"\
		"vmulpd		 0x080(%%r10),%%zmm7,%%zmm7	\n\t		vmulpd		0x280(%%r10),%%zmm15,%%zmm15\n\t"\
		"vsubpd			%%zmm5,%%zmm0,%%zmm0	\n\t		vsubpd		%%zmm13,%%zmm10,%%zmm10		\n\t"\
		"vsubpd			%%zmm4,%%zmm1,%%zmm1	\n\t		vsubpd		%%zmm12,%%zmm11,%%zmm11		\n\t"\
	" vfmadd231pd	 0x0c0(%%r10),%%zmm3,%%zmm6	\n\t	 vfmadd231pd	0x2c0(%%r10),%%zmm9 ,%%zmm14\n\t"\
	"vfnmadd231pd	 0x0c0(%%r10),%%zmm2,%%zmm7	\n\t	vfnmadd231pd	0x2c0(%%r10),%%zmm8 ,%%zmm15\n\t"\
	"vfmadd132pd		(%%rdi),%%zmm0,%%zmm5	\n\t	vfmadd132pd		(%%rdi),%%zmm10,%%zmm13		\n\t"\
	"vfmadd132pd		(%%rdi),%%zmm1,%%zmm4	\n\t	vfmadd132pd		(%%rdi),%%zmm11,%%zmm12		\n\t"\
		"vmovaps		%%zmm7, 0x440(%%rbx)	\n\t		vmovaps		%%zmm15, 0x540(%%rbx)		\n\t"\
		"vmovaps		%%zmm6, 0x400(%%rbx)	\n\t		vmovaps		%%zmm14, 0x500(%%rbx)		\n\t"\
		"addq		$0x100,%%r10				\n\t"\
		"vmovaps		%%zmm5,%%zmm2			\n\t		vmovaps		%%zmm13,%%zmm8 				\n\t"\
		"vmovaps		%%zmm1,%%zmm3			\n\t		vmovaps		%%zmm11,%%zmm9 				\n\t"\
		"vmovaps		%%zmm0,%%zmm6			\n\t		vmovaps		%%zmm10,%%zmm14				\n\t"\
		"vmovaps		%%zmm4,%%zmm7			\n\t		vmovaps		%%zmm12,%%zmm15				\n\t"\
		"vmulpd			  (%%r10),%%zmm5,%%zmm5	\n\t		vmulpd		0x200(%%r10),%%zmm13,%%zmm13\n\t"\
		"vmulpd			  (%%r10),%%zmm1,%%zmm1	\n\t		vmulpd		0x200(%%r10),%%zmm11,%%zmm11\n\t"\
		"vmulpd		 0x080(%%r10),%%zmm0,%%zmm0	\n\t		vmulpd		0x280(%%r10),%%zmm10,%%zmm10\n\t"\
		"vmulpd		 0x080(%%r10),%%zmm4,%%zmm4	\n\t		vmulpd		0x280(%%r10),%%zmm12,%%zmm12\n\t"\
	" vfmadd231pd	 0x040(%%r10),%%zmm3,%%zmm5	\n\t	 vfmadd231pd	0x240(%%r10),%%zmm9 ,%%zmm13\n\t"\
	"vfnmadd231pd	 0x040(%%r10),%%zmm2,%%zmm1	\n\t	vfnmadd231pd	0x240(%%r10),%%zmm8 ,%%zmm11\n\t"\
	" vfmadd231pd	 0x0c0(%%r10),%%zmm7,%%zmm0	\n\t	 vfmadd231pd	0x2c0(%%r10),%%zmm15,%%zmm10\n\t"\
	"vfnmadd231pd	 0x0c0(%%r10),%%zmm6,%%zmm4	\n\t	vfnmadd231pd	0x2c0(%%r10),%%zmm14,%%zmm12\n\t"\
		"vmovaps		%%zmm1, 0x240(%%rbx)	\n\t		vmovaps		%%zmm11, 0x340(%%rbx)		\n\t"\
		"vmovaps		%%zmm5, 0x200(%%rbx)	\n\t		vmovaps		%%zmm13, 0x300(%%rbx)		\n\t"\
		"vmovaps		%%zmm4, 0x640(%%rbx)	\n\t		vmovaps		%%zmm12, 0x740(%%rbx)		\n\t"\
		"vmovaps		%%zmm0, 0x600(%%rbx)	\n\t		vmovaps		%%zmm10, 0x700(%%rbx)		\n\t"\
		"\n\t"\
	/********************************************/	/********************************************/\
	/* Block 4: t06,t16,t26,t36 -> r30,34,32,36 */	/* Block 8: t0E,t1E,t2E,t3E -> r38,3C,3A,3E */\
	/********************************************/	/********************************************/\
		"addq		$0x800,%%rax				\n\t		addq		$0x200,%%r10				\n\t"\
		"vmovaps		0x080(%%rax),%%zmm4		\n\t		vmovaps		0x280(%%rax),%%zmm12		\n\t"\
		"vmovaps		0x180(%%rax),%%zmm0		\n\t		vmovaps		0x380(%%rax),%%zmm10		\n\t"\
		"vmovaps		0x0c0(%%rax),%%zmm7		\n\t		vmovaps		0x2c0(%%rax),%%zmm15		\n\t"\
		"vmovaps		0x1c0(%%rax),%%zmm1		\n\t		vmovaps		0x3c0(%%rax),%%zmm9 		\n\t"\
		"vmovaps			%%zmm4 ,%%zmm6		\n\t"	/*	vmovaps		 	%%zmm12,%%zmm14			\n\t"*/\
		"vmovaps			%%zmm0 ,%%zmm2		\n\t"	/*	vmovaps		 	%%zmm8 ,%%zmm10			\n\t"*/\
	/*	"vmovaps			%%zmm5 ,%%zmm7	*/			/*	vmovaps		 	%%zmm13,%%zmm15			\n\t"*/\
		"vmovaps			%%zmm1 ,%%zmm3		\n\t		vmovaps		 	%%zmm9 ,%%zmm11			\n\t"\
	/* Strategic register-usage for sincos consts here cuts loads from 16 to 4: */\
		"vmovaps			 (%%rdx),%%zmm14	\n\t		vmovaps		0x040(%%rdx),%%zmm13		\n\t"\
		"vmovaps		0x080(%%rdx),%%zmm5 	\n\t		vmovaps		0x0c0(%%rdx),%%zmm8 		\n\t"\
	/* Cyclic-permute 4 paired-MUL lines so as to put a '...(%%rax)' in the SRC3 slot of last line'srcol: */\
		"vmulpd		 	 %%zmm14,%%zmm2,%%zmm2	\n\t		vmulpd			 %%zmm5 ,%%zmm10,%%zmm10\n\t"\
		"vmulpd		 	 %%zmm8 ,%%zmm7,%%zmm7	\n\t		vmulpd			 %%zmm14,%%zmm15,%%zmm15\n\t"\
		"vmulpd		 	 %%zmm14,%%zmm3,%%zmm3	\n\t		vmulpd			 %%zmm5 ,%%zmm11,%%zmm11\n\t"\
		"vmulpd		 	 %%zmm8 ,%%zmm6,%%zmm6	\n\t		vmulpd		0x280(%%rax),%%zmm14,%%zmm14\n\t"\
	/* Also Cyclic-permute these 4 lines so middle (SRC2) op's reg-index order matches that of abovefour: */\
	"vfmadd132pd	 	 %%zmm13,%%zmm2,%%zmm1	\n\t	vfmsub132pd	 		 %%zmm8 ,%%zmm10,%%zmm9 \n\t"\
	"vfmadd132pd	 	 %%zmm5 ,%%zmm7,%%zmm4	\n\t	vfmadd132pd	 		 %%zmm13,%%zmm15,%%zmm12\n\t"\
	"vfmsub132pd	 	 %%zmm13,%%zmm3,%%zmm0	\n\t	vfmadd132pd	 	0x380(%%rax),%%zmm11,%%zmm8 \n\t"\
	"vfmsub132pd	 0xc0(%%rax),%%zmm6,%%zmm5	\n\t	vfmsub132pd	 	0x2c0(%%rax),%%zmm14,%%zmm13\n\t"\
		"vmovaps		%%zmm5,%%zmm7			\n\t		vmovaps		%%zmm13,%%zmm15				\n\t"\
		"vmovaps		%%zmm4,%%zmm6			\n\t		vmovaps		%%zmm12,%%zmm14				\n\t"\
		"vaddpd			%%zmm0,%%zmm4,%%zmm4	\n\t		vaddpd		%%zmm8 ,%%zmm12,%%zmm12		\n\t"\
		"vaddpd			%%zmm1,%%zmm5,%%zmm5	\n\t		vaddpd		%%zmm9 ,%%zmm13,%%zmm13		\n\t"\
		"vsubpd			%%zmm0,%%zmm6,%%zmm6	\n\t		vsubpd		%%zmm8 ,%%zmm14,%%zmm14		\n\t"\
		"vsubpd			%%zmm1,%%zmm7,%%zmm7	\n\t		vsubpd		%%zmm9 ,%%zmm15,%%zmm15		\n\t"\
		"vmovaps		 0x100(%%rax),%%zmm0	\n\t		vmovaps		0x300(%%rax),%%zmm8 		\n\t"\
		"vmovaps		 0x140(%%rax),%%zmm1	\n\t		vmovaps		0x340(%%rax),%%zmm9 		\n\t"\
		"vmovaps		 0x100(%%rax),%%zmm2	\n\t		vmovaps		0x300(%%rax),%%zmm10		\n\t"\
		"vmovaps		 0x140(%%rax),%%zmm3	\n\t		vmovaps		0x340(%%rax),%%zmm11		\n\t"\
		"vmulpd			  (%%rcx),%%zmm0,%%zmm0	\n\t		vmulpd		0x040(%%rcx),%%zmm8 ,%%zmm8 \n\t"\
		"vmulpd			  (%%rcx),%%zmm1,%%zmm1	\n\t		vmulpd		0x040(%%rcx),%%zmm9 ,%%zmm9 \n\t"\
	"vfmsub132pd	 0x040(%%rcx),%%zmm0,%%zmm3	\n\t	vfmadd132pd		(%%rcx),%%zmm8 ,%%zmm11		\n\t"\
	"vfmadd132pd	 0x040(%%rcx),%%zmm1,%%zmm2	\n\t	vfmsub132pd		(%%rcx),%%zmm9 ,%%zmm10		\n\t"\
		"vmovaps			  (%%rax),%%zmm0	\n\t		vmovaps		0x200(%%rax),%%zmm8 		\n\t"\
		"vmovaps		 0x040(%%rax),%%zmm1	\n\t		vmovaps		0x240(%%rax),%%zmm9 		\n\t"\
		"vsubpd			%%zmm2,%%zmm0,%%zmm0	\n\t		vsubpd		%%zmm10,%%zmm8 ,%%zmm8 		\n\t"\
		"vsubpd			%%zmm3,%%zmm1,%%zmm1	\n\t		vsubpd		%%zmm11,%%zmm9 ,%%zmm9 		\n\t"\
	"vfmadd132pd		(%%rdi),%%zmm0,%%zmm2	\n\t	vfmadd132pd		(%%rdi),%%zmm8 ,%%zmm10		\n\t"\
	"vfmadd132pd		(%%rdi),%%zmm1,%%zmm3	\n\t	vfmadd132pd		(%%rdi),%%zmm9 ,%%zmm11		\n\t"\
		"vsubpd			%%zmm6,%%zmm2,%%zmm2	\n\t		vsubpd		%%zmm14,%%zmm8 ,%%zmm8 		\n\t"\
		"vsubpd			%%zmm7,%%zmm3,%%zmm3	\n\t		vsubpd		%%zmm15,%%zmm9 ,%%zmm9 		\n\t"\
	"vfmadd132pd		(%%rdi),%%zmm2,%%zmm6	\n\t	vfmadd132pd		(%%rdi),%%zmm8 ,%%zmm14		\n\t"\
	"vfmadd132pd		(%%rdi),%%zmm3,%%zmm7	\n\t	vfmadd132pd		(%%rdi),%%zmm9 ,%%zmm15		\n\t"\
		"addq		$0x100,%%r10				\n\t"\
		"vmovaps		%%zmm2, 0x080(%%rax)	\n\t		vmovaps		%%zmm8 , 0x280(%%rax)		\n\t"\
		"vmovaps		%%zmm3, 0x0c0(%%rax)	\n\t		vmovaps		%%zmm9 , 0x2c0(%%rax)		\n\t"\
		"vmovaps		%%zmm6,%%zmm2			\n\t		vmovaps		%%zmm14,%%zmm8 				\n\t"\
		"vmovaps		%%zmm7,%%zmm3			\n\t		vmovaps		%%zmm15,%%zmm9 				\n\t"\
		"vmulpd			  (%%r10),%%zmm6,%%zmm6	\n\t		vmulpd		0x200(%%r10),%%zmm14,%%zmm14\n\t"\
		"vmulpd			  (%%r10),%%zmm7,%%zmm7	\n\t		vmulpd		0x200(%%r10),%%zmm15,%%zmm15\n\t"\
	" vfmadd231pd	 0x040(%%r10),%%zmm3,%%zmm6	\n\t	 vfmadd231pd	0x240(%%r10),%%zmm9 ,%%zmm14\n\t"\
	"vfnmadd231pd	 0x040(%%r10),%%zmm2,%%zmm7	\n\t	vfnmadd231pd	0x240(%%r10),%%zmm8 ,%%zmm15\n\t"\
		"vmovaps		%%zmm7, 0x0c0(%%rbx)	\n\t		vmovaps		%%zmm15, 0x1c0(%%rbx)		\n\t"\
		"vmovaps		%%zmm6, 0x080(%%rbx)	\n\t		vmovaps		%%zmm14, 0x180(%%rbx)		\n\t"\
		"vmovaps		 0x080(%%rax),%%zmm6	\n\t		vmovaps		0x280(%%rax),%%zmm14		\n\t"\
		"vmovaps		 0x0c0(%%rax),%%zmm7	\n\t		vmovaps		0x2c0(%%rax),%%zmm15		\n\t"\
		"vmovaps		%%zmm6,%%zmm2			\n\t		vmovaps		%%zmm14,%%zmm8 				\n\t"\
		"vmovaps		%%zmm7,%%zmm3			\n\t		vmovaps		%%zmm15,%%zmm9 				\n\t"\
		"vmulpd		 0x080(%%r10),%%zmm6,%%zmm6	\n\t		vmulpd		0x280(%%r10),%%zmm14,%%zmm14\n\t"\
		"vmulpd		 0x080(%%r10),%%zmm7,%%zmm7	\n\t		vmulpd		0x280(%%r10),%%zmm15,%%zmm15\n\t"\
		"vsubpd			%%zmm5,%%zmm0,%%zmm0	\n\t		vsubpd		%%zmm13,%%zmm10,%%zmm10		\n\t"\
		"vsubpd			%%zmm4,%%zmm1,%%zmm1	\n\t		vsubpd		%%zmm12,%%zmm11,%%zmm11		\n\t"\
	" vfmadd231pd	 0x0c0(%%r10),%%zmm3,%%zmm6	\n\t	 vfmadd231pd	0x2c0(%%r10),%%zmm9 ,%%zmm14\n\t"\
	"vfnmadd231pd	 0x0c0(%%r10),%%zmm2,%%zmm7	\n\t	vfnmadd231pd	0x2c0(%%r10),%%zmm8 ,%%zmm15\n\t"\
	"vfmadd132pd		(%%rdi),%%zmm0,%%zmm5	\n\t	vfmadd132pd		(%%rdi),%%zmm10,%%zmm13		\n\t"\
	"vfmadd132pd		(%%rdi),%%zmm1,%%zmm4	\n\t	vfmadd132pd		(%%rdi),%%zmm11,%%zmm12		\n\t"\
		"vmovaps		%%zmm7, 0x4c0(%%rbx)	\n\t		vmovaps		%%zmm15, 0x5c0(%%rbx)		\n\t"\
		"vmovaps		%%zmm6, 0x480(%%rbx)	\n\t		vmovaps		%%zmm14, 0x580(%%rbx)		\n\t"\
		"addq		$0x100,%%r10				\n\t"\
		"vmovaps		%%zmm5,%%zmm2			\n\t		vmovaps		%%zmm13,%%zmm8 				\n\t"\
		"vmovaps		%%zmm1,%%zmm3			\n\t		vmovaps		%%zmm11,%%zmm9 				\n\t"\
		"vmovaps		%%zmm0,%%zmm6			\n\t		vmovaps		%%zmm10,%%zmm14				\n\t"\
		"vmovaps		%%zmm4,%%zmm7			\n\t		vmovaps		%%zmm12,%%zmm15				\n\t"\
		"vmulpd			  (%%r10),%%zmm5,%%zmm5	\n\t		vmulpd		0x200(%%r10),%%zmm13,%%zmm13\n\t"\
		"vmulpd			  (%%r10),%%zmm1,%%zmm1	\n\t		vmulpd		0x200(%%r10),%%zmm11,%%zmm11\n\t"\
		"vmulpd		 0x080(%%r10),%%zmm0,%%zmm0	\n\t		vmulpd		0x280(%%r10),%%zmm10,%%zmm10\n\t"\
		"vmulpd		 0x080(%%r10),%%zmm4,%%zmm4	\n\t		vmulpd		0x280(%%r10),%%zmm12,%%zmm12\n\t"\
	" vfmadd231pd	 0x040(%%r10),%%zmm3,%%zmm5	\n\t	 vfmadd231pd	0x240(%%r10),%%zmm9 ,%%zmm13\n\t"\
	"vfnmadd231pd	 0x040(%%r10),%%zmm2,%%zmm1	\n\t	vfnmadd231pd	0x240(%%r10),%%zmm8 ,%%zmm11\n\t"\
	" vfmadd231pd	 0x0c0(%%r10),%%zmm7,%%zmm0	\n\t	 vfmadd231pd	0x2c0(%%r10),%%zmm15,%%zmm10\n\t"\
	"vfnmadd231pd	 0x0c0(%%r10),%%zmm6,%%zmm4	\n\t	vfnmadd231pd	0x2c0(%%r10),%%zmm14,%%zmm12\n\t"\
		"vmovaps		%%zmm1, 0x2c0(%%rbx)	\n\t		vmovaps		%%zmm11, 0x3c0(%%rbx)		\n\t"\
		"vmovaps		%%zmm5, 0x280(%%rbx)	\n\t		vmovaps		%%zmm13, 0x380(%%rbx)		\n\t"\
		"vmovaps		%%zmm4, 0x6c0(%%rbx)	\n\t		vmovaps		%%zmm12, 0x7c0(%%rbx)		\n\t"\
		"vmovaps		%%zmm0, 0x680(%%rbx)	\n\t		vmovaps		%%zmm10, 0x780(%%rbx)		\n\t"\
		"\n\t"\
	/********************************************/	/********************************************/\
	/* Block 1: t00,t10,t20,t30 -> r00,04,02,06 */	/* Block 5: t08,t18,t28,t38 -> r08,0C,0A,0E */\
	/********************************************/	/********************************************/\
		"movq	%[__r00],%%rdx					\n\t		vmovaps	(%%rsi),%%zmm10	\n\t"/* isrt2 */\
													/* base-addr in rcol = r08, thus rdx + 0x200 there */\
		"vmovaps			  (%%rdx),%%zmm0	\n\t		vmovaps		 0x280(%%rdx),%%zmm12		\n\t"\
		"vmovaps		 0x040(%%rdx),%%zmm1	\n\t		vmovaps		 0x2c0(%%rdx),%%zmm13		\n\t"\
		"vmovaps		 0x100(%%rdx),%%zmm2	\n\t		vmovaps		 0x380(%%rdx),%%zmm8 		\n\t"\
		"vmovaps		 0x140(%%rdx),%%zmm3	\n\t		vmovaps		 0x3c0(%%rdx),%%zmm9 		\n\t"\
		"vsubpd		 0x100(%%rdx),%%zmm0,%%zmm0	\n\t		vaddpd	 0x2c0(%%rdx),%%zmm12,%%zmm12	\n\t"\
		"vsubpd		 0x140(%%rdx),%%zmm1,%%zmm1	\n\t		vsubpd	 0x280(%%rdx),%%zmm13,%%zmm13	\n\t"\
		"vaddpd			  (%%rdx),%%zmm2,%%zmm2	\n\t		vsubpd	 0x3c0(%%rdx),%%zmm8 ,%%zmm8 	\n\t"\
		"vaddpd		 0x040(%%rdx),%%zmm3,%%zmm3	\n\t		vaddpd	 0x380(%%rdx),%%zmm9 ,%%zmm9 	\n\t"\
		"vmovaps		 0x080(%%rdx),%%zmm4	\n\t		vmulpd		%%zmm10,%%zmm8 ,%%zmm8 		\n\t"\
		"vmovaps		 0x0c0(%%rdx),%%zmm5	\n\t		vmulpd		%%zmm10,%%zmm9 ,%%zmm9 		\n\t"\
		"vmovaps		 0x180(%%rdx),%%zmm6	\n\t		vmovaps		%%zmm12,%%zmm14				\n\t"\
		"vmovaps		 0x1c0(%%rdx),%%zmm7	\n\t		vmovaps		%%zmm13,%%zmm15				\n\t"\
		"vsubpd		 0x180(%%rdx),%%zmm4,%%zmm4	\n\t	vfmsub132pd		%%zmm10,%%zmm8 ,%%zmm12		\n\t"\
		"vsubpd		 0x1c0(%%rdx),%%zmm5,%%zmm5	\n\t	vfmsub132pd		%%zmm10,%%zmm9 ,%%zmm13		\n\t"\
		"vaddpd		 0x080(%%rdx),%%zmm6,%%zmm6	\n\t	vfmadd132pd		%%zmm10,%%zmm8 ,%%zmm14		\n\t"\
		"vaddpd		 0x0c0(%%rdx),%%zmm7,%%zmm7	\n\t	vfmadd132pd		%%zmm10,%%zmm9 ,%%zmm15		\n\t"\
		"movq		%[__c10],%%rcx			\n\t"/* base-twiddle in l/rcol = c00/c04 ==> rcx+0x200 in rcol */\
		"vaddpd			%%zmm6,%%zmm2,%%zmm2	\n\t		vmovaps		 0x200(%%rdx),%%zmm8 		\n\t"\
		"vaddpd			%%zmm7,%%zmm3,%%zmm3	\n\t		vmovaps		 0x240(%%rdx),%%zmm9 		\n\t"\
		"vmovaps		%%zmm2,      (%%rdx)	\n\t		vmovaps		 0x300(%%rdx),%%zmm10		\n\t"\
		"vmovaps		%%zmm3, 0x040(%%rdx)	\n\t		vmovaps		 0x340(%%rdx),%%zmm11		\n\t"\
		"vaddpd			%%zmm6,%%zmm6,%%zmm6	\n\t		vsubpd	 0x340(%%rdx),%%zmm8 ,%%zmm8 	\n\t"\
		"vaddpd			%%zmm7,%%zmm7,%%zmm7	\n\t		vsubpd	 0x300(%%rdx),%%zmm9 ,%%zmm9 	\n\t"\
		"vsubpd			%%zmm6,%%zmm2,%%zmm2	\n\t		vaddpd	 0x200(%%rdx),%%zmm11,%%zmm11	\n\t"\
		"vsubpd			%%zmm7,%%zmm3,%%zmm3	\n\t		vaddpd	 0x240(%%rdx),%%zmm10,%%zmm10	\n\t"\
		"vmovaps		%%zmm2,%%zmm6			\n\t		vsubpd		%%zmm12,%%zmm11,%%zmm11		\n\t"\
		"vmovaps		%%zmm3,%%zmm7			\n\t		vsubpd		%%zmm13,%%zmm9 ,%%zmm9 		\n\t"\
		"vmulpd	(%%rcx),%%zmm2,%%zmm2	/* c10 */\n\t	vfmadd132pd		(%%rdi),%%zmm11,%%zmm12		\n\t"\
		"vmulpd	(%%rcx),%%zmm3,%%zmm3			\n\t	vfmadd132pd		(%%rdi),%%zmm9 ,%%zmm13		\n\t"\
	" vfmadd231pd	 0x040(%%rcx),%%zmm7,%%zmm2	\n\t		vmovaps		%%zmm11, 0x280(%%rdx)		\n\t"\
	"vfnmadd231pd	 0x040(%%rcx),%%zmm6,%%zmm3	\n\t		vmovaps		%%zmm9 , 0x2c0(%%rdx)		\n\t"\
		"subq $0x80,%%rcx		\n\t"/* put c00 in rcx to ease bookkeeping */\
	/* add0,1 in rax,rbx; __r00 in rdx.  For each */	"	vmovaps		%%zmm12,%%zmm11				\n\t"\
	/* complex output octet, complex pairs read from */"	vmovaps		%%zmm13,%%zmm9 	\n\t"/* c04:*/\
	/* offsets 0x0..,0x1..,0x2..,0x3.. go into local-*/"	vmulpd 0x200(%%rcx),%%zmm12,%%zmm12 	\n\t"\
	/* mem pairs rXY + 00/10,04/14,02/12,06/16. For  */"	vmulpd 0x200(%%rcx),%%zmm13,%%zmm13		\n\t"\
	/* 1st octet we read from offsets [0x2..,0x0..], */" vfmadd231pd	0x240(%%rcx),%%zmm9 ,%%zmm12\n\t"\
	/* [0x1..,0x3], other 3 octets order [0,2],[1,3] */"vfnmadd231pd	0x240(%%rcx),%%zmm11,%%zmm13\n\t"\
		"vmovaps	0x440(%%rbx),%%zmm7			\n\t		vmovaps	0x140(%%rbx),%%zmm11			\n\t"\
		"vmovaps	0x400(%%rbx),%%zmm6			\n\t		vmovaps	0x100(%%rbx),%%zmm9 			\n\t"\
		"vmovaps	%%zmm3,0x0c0(%%rdx)			\n\t		vmovaps	%%zmm13,0x240(%%rdx)			\n\t"\
		"vmovaps	%%zmm2,0x080(%%rdx)	/* r02,03 */\n\t	vmovaps	%%zmm12,0x200(%%rdx)			\n\t"/* r08,09 */\
		"vmovaps	%%zmm7,0x4c0(%%rdx)			\n\t		vmovaps	%%zmm11,0x640(%%rdx)			\n\t"\
		"vmovaps	%%zmm6,0x480(%%rdx)	/* r12,13 */\n\t	vmovaps	%%zmm9 ,0x600(%%rdx)			\n\t"/* r18,19 */\
		"vmovaps	0x040(%%rdx),%%zmm3			\n\t		vmovaps		0x280(%%rdx),%%zmm11		\n\t"\
		"vmovaps		 (%%rdx),%%zmm2			\n\t		vmovaps		0x2c0(%%rdx),%%zmm9 		\n\t"\
		"vmovaps	0x040(%%rbx),%%zmm7			\n\t		vmulpd	 0x280(%%rcx),%%zmm11,%%zmm12	\n\t"/* c14 */\
		"vmovaps	0x000(%%rbx),%%zmm6			\n\t		vmulpd	 0x280(%%rcx),%%zmm9 ,%%zmm13	\n\t"\
		"vmovaps	%%zmm3,0x040(%%rdx)			\n\t	 vfmadd231pd 0x2c0(%%rcx),%%zmm9 ,%%zmm12	\n\t"\
		"vmovaps	%%zmm2,     (%%rdx)	/* r00,01*/\n\t	vfnmadd231pd 0x2c0(%%rcx),%%zmm11,%%zmm13	\n\t"\
		"vmovaps	%%zmm7,0x440(%%rdx)			\n\t		vmovaps	0x540(%%rbx),%%zmm11			\n\t"\
		"vmovaps	%%zmm6,0x400(%%rdx)	/* r10,11*/\n\t		vmovaps	0x500(%%rbx),%%zmm9 			\n\t"\
		"vaddpd			%%zmm5,%%zmm0,%%zmm0	\n\t		vmovaps	%%zmm13,0x2c0(%%rdx)			\n\t"\
		"vsubpd			%%zmm4,%%zmm1,%%zmm1	\n\t		vmovaps	%%zmm12,0x280(%%rdx)			\n\t"/* r0a,0b */\
		"													vmovaps	%%zmm11,0x6c0(%%rdx)			\n\t"\
		"													vmovaps	%%zmm9 ,0x680(%%rdx)			\n\t"/* r1a,1b */\
	"vfnmadd132pd	(%%rdi),%%zmm0,%%zmm5		\n\t		vsubpd			%%zmm15,%%zmm8 ,%%zmm8 	\n\t"\
	" vfmadd132pd	(%%rdi),%%zmm1,%%zmm4		\n\t		vsubpd			%%zmm14,%%zmm10,%%zmm10	\n\t"\
		"vmulpd	0x100(%%rcx),%%zmm0,%%zmm2 /*c08*/\n\t	vfmadd132pd		(%%rdi),%%zmm8 ,%%zmm15		\n\t"\
		"vmulpd	0x100(%%rcx),%%zmm1,%%zmm3		\n\t	vfmadd132pd		(%%rdi),%%zmm10,%%zmm14		\n\t"\
	" vfmadd231pd	 0x140(%%rcx),%%zmm1,%%zmm2	\n\t		vmulpd	 0x300(%%rcx),%%zmm15,%%zmm12	\n\t"/* c0C */\
	"vfnmadd231pd	 0x140(%%rcx),%%zmm0,%%zmm3	\n\t		vmulpd	 0x300(%%rcx),%%zmm10,%%zmm13	\n\t"\
		"vmovaps	0x240(%%rbx),%%zmm7			\n\t	 vfmadd231pd 0x340(%%rcx),%%zmm10,%%zmm12	\n\t"\
		"vmovaps	0x200(%%rbx),%%zmm6			\n\t	vfnmadd231pd 0x340(%%rcx),%%zmm15,%%zmm13	\n\t"\
		"vmovaps	%%zmm3,0x140(%%rdx)			\n\t		vmovaps	0x340(%%rbx),%%zmm10			\n\t"\
		"vmovaps	%%zmm2,0x100(%%rdx)	/* r04,05 */\n\t	vmovaps	0x300(%%rbx),%%zmm15			\n\t"\
		"vmovaps	%%zmm7,0x540(%%rdx)			\n\t		vmovaps	%%zmm13,0x340(%%rdx)			\n\t"\
		"vmovaps	%%zmm6,0x500(%%rdx)	/* r14,15 */\n\t	vmovaps	%%zmm12,0x300(%%rdx)			\n\t"/* r0c,0d */\
		"													vmovaps	%%zmm10,0x740(%%rdx)			\n\t"\
		"													vmovaps	%%zmm15,0x700(%%rdx)			\n\t"/* r1c,1d */\
		"vmulpd	0x180(%%rcx),%%zmm5,%%zmm0	/* c18*/\n\t	vmulpd	 0x380(%%rcx),%%zmm8 ,%%zmm15	\n\t"/* c1C */\
		"vmulpd	0x180(%%rcx),%%zmm4,%%zmm1		\n\t		vmulpd	 0x380(%%rcx),%%zmm14,%%zmm10	\n\t"\
	" vfmadd231pd	 0x1c0(%%rcx),%%zmm4,%%zmm0	\n\t	 vfmadd231pd 0x3c0(%%rcx),%%zmm14,%%zmm15	\n\t"\
	"vfnmadd231pd	 0x1c0(%%rcx),%%zmm5,%%zmm1	\n\t	vfnmadd231pd 0x3c0(%%rcx),%%zmm8 ,%%zmm10	\n\t"\
	/*** Need to get rid of these read-from-main-array-spill-slots-and-write-to-local-array load/store pairs: ***/\
		"vmovaps	0x640(%%rbx),%%zmm7			\n\t		vmovaps	0x740(%%rbx),%%zmm14			\n\t"\
		"vmovaps	0x600(%%rbx),%%zmm6			\n\t		vmovaps	0x700(%%rbx),%%zmm8 			\n\t"\
		"vmovaps	%%zmm7,0x5c0(%%rdx)			\n\t		vmovaps	%%zmm14,0x7c0(%%rdx)			\n\t"\
		"vmovaps	%%zmm6,0x580(%%rdx)	/* r16,17 */\n\t	vmovaps	%%zmm8 ,0x780(%%rdx)			\n\t"/* r1e,1f */\
		"vmovaps	%%zmm1,0x1c0(%%rdx)			\n\t		vmovaps	%%zmm10,0x3c0(%%rdx)			\n\t"\
		"vmovaps	%%zmm0,0x180(%%rdx)	/* r06,07 */\n\t	vmovaps	%%zmm15,0x380(%%rdx)			\n\t"/* r0e,0f */\
		"\n\t"\
	/********************************************/	/********************************************/\
	/* Block 3: t04,t14,t24,t34 -> r20,24,22,26 */	/* Block 7: t0C,t1C,t2C,t3C -> r28,2C,2A,2E */\
	/********************************************/	/********************************************/\
		"movq		%[__r20],%%rdx	\n\t"			/* base-addr in rcol = r28, so rdx offset +0x200 vs lcol */\
		"leaq	0x040(%%rsi),%%rcx	\n\t"/* cc0; Note cc0/ss0 are shared between lcol/rcol, so no rcx-offset until get to twiddles*/\
		"vmovaps		0x080(%%rdx),%%zmm4		\n\t		vmovaps		 0x280(%%rdx),%%zmm12		\n\t"\
		"vmovaps		0x180(%%rdx),%%zmm0		\n\t		vmovaps		 0x380(%%rdx),%%zmm8 		\n\t"\
		"vmovaps		0x0c0(%%rdx),%%zmm5		\n\t		vmovaps		 0x2c0(%%rdx),%%zmm13		\n\t"\
		"vmovaps		0x1c0(%%rdx),%%zmm3		\n\t		vmovaps		 0x3c0(%%rdx),%%zmm11		\n\t"\
		"vmovaps			%%zmm4 ,%%zmm6		\n\t		vmovaps		 	%%zmm12,%%zmm14			\n\t"\
		"vmovaps			%%zmm0 ,%%zmm2		\n\t		vmovaps		 	%%zmm8 ,%%zmm10			\n\t"\
		"vmovaps			%%zmm5 ,%%zmm7		\n\t		vmovaps		 	%%zmm13,%%zmm15			\n\t"\
	/*	"vmovaps			%%zmm1 ,%%zmm3		\n\t		vmovaps		 	%%zmm9 ,%%zmm11			\n\t"*/\
	/* Strategic register-usage for sincos consts here cuts loads from 16 to 4: */\
		"vmovaps			 (%%rcx),%%zmm9 	\n\t		vmovaps		0x040(%%rcx),%%zmm1 		\n\t"\
		"vmulpd		 	 %%zmm1 ,%%zmm7,%%zmm7	\n\t		vmulpd			 %%zmm9 ,%%zmm15,%%zmm15\n\t"\
		"vmulpd		 	 %%zmm9 ,%%zmm3,%%zmm3	\n\t		vmulpd			 %%zmm1 ,%%zmm11,%%zmm11\n\t"\
		"vmulpd		 	 %%zmm1 ,%%zmm6,%%zmm6	\n\t		vmulpd			 %%zmm9 ,%%zmm14,%%zmm14\n\t"\
		"vmulpd		 	 %%zmm9 ,%%zmm2,%%zmm2	\n\t		vmulpd			 %%zmm1 ,%%zmm10,%%zmm10\n\t"\
	"vfmadd132pd	 	 %%zmm9 ,%%zmm7,%%zmm4	\n\t	vfmadd132pd			 %%zmm1 ,%%zmm15,%%zmm12\n\t"\
	"vfmadd132pd	 	 %%zmm1 ,%%zmm3,%%zmm0	\n\t	vfmadd132pd			 %%zmm9 ,%%zmm11,%%zmm8 \n\t"\
	"vfmsub132pd	 	 %%zmm9 ,%%zmm6,%%zmm5	\n\t	vfmsub132pd			 %%zmm1 ,%%zmm14,%%zmm13\n\t"\
	"vfmsub132pd	0x1c0(%%rdx),%%zmm2,%%zmm1	\n\t	vfmsub132pd		0x3c0(%%rdx),%%zmm10,%%zmm9 \n\t"\
		"vsubpd			%%zmm0,%%zmm4,%%zmm6	\n\t		vsubpd			%%zmm8 ,%%zmm12,%%zmm14	\n\t"\
		"vsubpd			%%zmm1,%%zmm5,%%zmm7	\n\t		vsubpd			%%zmm9 ,%%zmm13,%%zmm15	\n\t"\
		"vaddpd			%%zmm0,%%zmm4,%%zmm4	\n\t		vaddpd			%%zmm8 ,%%zmm12,%%zmm12	\n\t"\
		"vaddpd			%%zmm1,%%zmm5,%%zmm5	\n\t		vaddpd			%%zmm9 ,%%zmm13,%%zmm13	\n\t"\
		"vmovaps		 0x100(%%rdx),%%zmm2	\n\t		vmovaps		0x300(%%rdx),%%zmm10		\n\t"\
		"vmovaps		 0x140(%%rdx),%%zmm3	\n\t		vmovaps		0x340(%%rdx),%%zmm11		\n\t"\
		"vmovaps			  (%%rdx),%%zmm0	\n\t		vmovaps		0x200(%%rdx),%%zmm8 		\n\t"\
		"vmovaps		 0x040(%%rdx),%%zmm1	\n\t		vmovaps		0x240(%%rdx),%%zmm9 		\n\t"\
		"vaddpd		 0x140(%%rdx),%%zmm2,%%zmm2	\n\t		vsubpd		0x340(%%rdx),%%zmm10,%%zmm10\n\t"\
		"vsubpd		 0x100(%%rdx),%%zmm3,%%zmm3	\n\t		vaddpd		0x300(%%rdx),%%zmm11,%%zmm11\n\t"\
		"vmulpd			  (%%rsi),%%zmm2,%%zmm2	\n\t		vmulpd			 (%%rsi),%%zmm10,%%zmm10\n\t"\
		"vmulpd			  (%%rsi),%%zmm3,%%zmm3	\n\t		vmulpd			 (%%rsi),%%zmm11,%%zmm11\n\t"\
		"vsubpd			%%zmm2,%%zmm0,%%zmm0	\n\t		vsubpd			%%zmm10,%%zmm8 ,%%zmm8 	\n\t"\
		"vsubpd			%%zmm3,%%zmm1,%%zmm1	\n\t		vsubpd			%%zmm11,%%zmm9 ,%%zmm9 	\n\t"\
	"vfmadd132pd		(%%rdi),%%zmm0,%%zmm2	\n\t	vfmadd132pd		 (%%rdi),%%zmm8 ,%%zmm10	\n\t"\
	"vfmadd132pd		(%%rdi),%%zmm1,%%zmm3	\n\t	vfmadd132pd		 (%%rdi),%%zmm9 ,%%zmm11	\n\t"\
		"movq		%[__c02],%%rcx				\n\t"/* base-twiddle addr in rcol = c06, so rcx offset +0x200 vs lcol */\
		"vsubpd			%%zmm4,%%zmm2,%%zmm2	\n\t		vsubpd			%%zmm14,%%zmm8 ,%%zmm8 	\n\t"\
		"vsubpd			%%zmm5,%%zmm3,%%zmm3	\n\t		vsubpd			%%zmm15,%%zmm9 ,%%zmm9 	\n\t"\
	"vfmadd132pd		(%%rdi),%%zmm2,%%zmm4	\n\t	vfmadd132pd		(%%rdi),%%zmm8 ,%%zmm14	\n\t"\
	"vfmadd132pd		(%%rdi),%%zmm3,%%zmm5	\n\t	vfmadd132pd		(%%rdi),%%zmm9 ,%%zmm15	\n\t"\
		"vmovaps		%%zmm2, 0x080(%%rdx)	\n\t		vmovaps		%%zmm8 , 0x280(%%rdx)		\n\t"\
		"vmovaps		%%zmm3, 0x0c0(%%rdx)	\n\t		vmovaps		%%zmm9 , 0x2c0(%%rdx)		\n\t"\
		"vmulpd			  (%%rcx),%%zmm4,%%zmm2	\n\t		vmulpd		 0x200(%%rcx),%%zmm14,%%zmm8\n\t"\
		"vmulpd			  (%%rcx),%%zmm5,%%zmm3	\n\t		vmulpd		 0x200(%%rcx),%%zmm15,%%zmm9\n\t"\
	" vfmadd231pd	 0x040(%%rcx),%%zmm5,%%zmm2	\n\t	 vfmadd231pd	 0x240(%%rcx),%%zmm15,%%zmm8\n\t"\
	"vfnmadd231pd	 0x040(%%rcx),%%zmm4,%%zmm3	\n\t	vfnmadd231pd	 0x240(%%rcx),%%zmm14,%%zmm9\n\t"\
		"vmovaps	0x0c0(%%rbx),%%zmm5			\n\t		vmovaps	0x1c0(%%rbx),%%zmm15			\n\t"\
		"vmovaps	0x080(%%rbx),%%zmm4			\n\t		vmovaps	0x180(%%rbx),%%zmm14			\n\t"\
		"vmovaps	%%zmm3,0x040(%%rdx)			\n\t		vmovaps	%%zmm9 ,0x240(%%rdx)			\n\t"\
		"vmovaps	%%zmm2,     (%%rdx)	/* r20,21 */\n\t	vmovaps	%%zmm8 ,0x200(%%rdx)			\n\t"/* r28,29 */\
		"vmovaps	%%zmm5,0x440(%%rdx)			\n\t		vmovaps	%%zmm15,0x640(%%rdx)			\n\t"\
		"vmovaps	%%zmm4,0x400(%%rdx)	/* r30,31 */\n\t	vmovaps	%%zmm14,0x600(%%rdx)			\n\t"/* r38,39 */\
		"movq		%[__c12],%%rcx				\n\t"		/* rcol uses c16 */\
		"vmovaps		 0x080(%%rdx),%%zmm4	\n\t		vmovaps		 0x280(%%rdx),%%zmm14		\n\t"\
		"vmovaps		 0x0c0(%%rdx),%%zmm5	\n\t		vmovaps		 0x2c0(%%rdx),%%zmm15		\n\t"\
		"vmulpd			  (%%rcx),%%zmm4,%%zmm2	\n\t		vmulpd		 0x200(%%rcx),%%zmm14,%%zmm8\n\t"\
		"vmulpd			  (%%rcx),%%zmm5,%%zmm3	\n\t		vmulpd		 0x200(%%rcx),%%zmm15,%%zmm9\n\t"\
	" vfmadd231pd	 0x040(%%rcx),%%zmm5,%%zmm2	\n\t	 vfmadd231pd	 0x240(%%rcx),%%zmm15,%%zmm8\n\t"\
	"vfnmadd231pd	 0x040(%%rcx),%%zmm4,%%zmm3	\n\t	vfnmadd231pd	 0x240(%%rcx),%%zmm14,%%zmm9\n\t"\
		"vmovaps	0x4c0(%%rbx),%%zmm5			\n\t		vmovaps	0x5c0(%%rbx),%%zmm15			\n\t"\
		"vmovaps	0x480(%%rbx),%%zmm4			\n\t		vmovaps	0x580(%%rbx),%%zmm14			\n\t"\
		"vmovaps	%%zmm3,0x0c0(%%rdx)			\n\t		vmovaps	%%zmm9 ,0x2c0(%%rdx)			\n\t"\
		"vmovaps	%%zmm2,0x080(%%rdx)	/* r22,23 */\n\t	vmovaps	%%zmm8 ,0x280(%%rdx)			\n\t"/* r2a,2b */\
		"vmovaps	%%zmm5,0x4c0(%%rdx)			\n\t		vmovaps	%%zmm15,0x6c0(%%rdx)			\n\t"\
		"vmovaps	%%zmm4,0x480(%%rdx)	/* r32,33 */\n\t	vmovaps	%%zmm14,0x680(%%rdx)			\n\t"/* r3a,3b */\
		"movq		%[__c0A],%%rcx				\n\t"		/* rcol uses c0E */\
		"vsubpd			%%zmm7,%%zmm0,%%zmm0	\n\t		vsubpd			%%zmm13,%%zmm10,%%zmm10	\n\t"\
		"vsubpd			%%zmm6,%%zmm1,%%zmm1	\n\t		vsubpd			%%zmm12,%%zmm11,%%zmm11	\n\t"\
	"vfmadd132pd		  (%%rdi),%%zmm0,%%zmm7	\n\t	vfmadd132pd		(%%rdi),%%zmm10,%%zmm13	\n\t"\
	"vfmadd132pd		  (%%rdi),%%zmm1,%%zmm6	\n\t	vfmadd132pd		(%%rdi),%%zmm11,%%zmm12	\n\t"\
		"vmulpd			  (%%rcx),%%zmm7,%%zmm4	\n\t		vmulpd		 0x200(%%rcx),%%zmm13,%%zmm8\n\t"\
		"vmulpd			  (%%rcx),%%zmm1,%%zmm5	\n\t		vmulpd		 0x200(%%rcx),%%zmm11,%%zmm9\n\t"\
	" vfmadd231pd	 0x040(%%rcx),%%zmm1,%%zmm4	\n\t	 vfmadd231pd	 0x240(%%rcx),%%zmm11,%%zmm8\n\t"\
	"vfnmadd231pd	 0x040(%%rcx),%%zmm7,%%zmm5	\n\t	vfnmadd231pd	 0x240(%%rcx),%%zmm13,%%zmm9\n\t"\
		"vmovaps	0x2c0(%%rbx),%%zmm1			\n\t		vmovaps	0x3c0(%%rbx),%%zmm11			\n\t"\
		"vmovaps	0x280(%%rbx),%%zmm7			\n\t		vmovaps	0x380(%%rbx),%%zmm13			\n\t"\
		"vmovaps	%%zmm5,0x140(%%rdx)			\n\t		vmovaps	%%zmm9 ,0x340(%%rdx)			\n\t"\
		"vmovaps	%%zmm4,0x100(%%rdx)	/* r24,25 */\n\t	vmovaps	%%zmm8 ,0x300(%%rdx)			\n\t"/* r2c,2d */\
		"vmovaps	%%zmm1,0x540(%%rdx)			\n\t		vmovaps	%%zmm11,0x740(%%rdx)			\n\t"\
		"vmovaps	%%zmm7,0x500(%%rdx)	/* r34,35 */\n\t	vmovaps	%%zmm13,0x700(%%rdx)			\n\t"/* r3c,3d */\
		"movq		%[__c1A],%%rcx				\n\t"		/* rcol uses c1E */\
		"vmulpd			  (%%rcx),%%zmm0,%%zmm4	\n\t		vmulpd		 0x200(%%rcx),%%zmm10,%%zmm8\n\t"\
		"vmulpd			  (%%rcx),%%zmm6,%%zmm5	\n\t		vmulpd		 0x200(%%rcx),%%zmm12,%%zmm9\n\t"\
	" vfmadd231pd	 0x040(%%rcx),%%zmm6,%%zmm4	\n\t	 vfmadd231pd	 0x240(%%rcx),%%zmm12,%%zmm8\n\t"\
	"vfnmadd231pd	 0x040(%%rcx),%%zmm0,%%zmm5	\n\t	vfnmadd231pd	 0x240(%%rcx),%%zmm10,%%zmm9\n\t"\
		"vmovaps	0x6c0(%%rbx),%%zmm6			\n\t		vmovaps	0x7c0(%%rbx),%%zmm12			\n\t"\
		"vmovaps	0x680(%%rbx),%%zmm0			\n\t		vmovaps	0x780(%%rbx),%%zmm10			\n\t"\
		"vmovaps	%%zmm5,0x1c0(%%rdx)			\n\t		vmovaps	%%zmm9 ,0x3c0(%%rdx)			\n\t"\
		"vmovaps	%%zmm4,0x180(%%rdx)	/* r26,27 */\n\t	vmovaps	%%zmm8 ,0x380(%%rdx)			\n\t"/* r2e,2f */\
		"vmovaps	%%zmm6,0x5c0(%%rdx)			\n\t		vmovaps	%%zmm12,0x7c0(%%rdx)			\n\t"\
		"vmovaps	%%zmm0,0x580(%%rdx)	/* r36,37 */\n\t	vmovaps	%%zmm10,0x780(%%rdx)			\n\t"/* r3e,3f */\
	/*******************************************
	***** Finish with 8-way 'un'terleaving: ****
	Using the AVX-512 data layout, the rcol pattern is:
		a[ 0- 7] = re[ 0, 8,16,24, 4,12,20,28].d0	a[ 8-15] = im[ 0, 8,16,24, 4,12,20,28].d0
		a[16-23] = re[ 2,10,18,26, 6,14,22,30].d0	a[24-31] = im[ 2,10,18,26, 6,14,22,30].d0
		a[32-39] = re[ 1, 9,17,25, 5,13,21,29].d0	a[40-47] = im[ 1, 9,17,25, 5,13,21,29].d0
		a[48-55] = re[ 3,11,19,27, 7,15,23,31].d0	a[56-63] = im[ 3,11,19,27, 7,15,23,31].d0 ,
	and remaining seven 64-double blocks repeat same pattern with elts d1-d7 of the vector-doubles.
	*******************************************/\
	"movq	%[__add0],%%rax	\n\t"\
	"movq	%[__r00] ,%%rcx	\n\t"\
		/* Auxiliary register data needed for columnwise stores: */\
		"movq	$0x0706050403020100,%%rsi	\n\t"/* 64-bit register w/byte offsets 0-7, bytes ordered left-to-right in decreasing significance */\
			"vmovq		%%rsi,%%xmm8 		\n\t"/* Copy byte pattern to low qword (64 bits) of ymm8 [NB: avx-512 only supports MOVQ to/from 128-bit vector regs] */\
			"vpmovzxbd	%%xmm8,%%ymm8		\n\t"/* vector-index offsets: ymm8 = [0,1,2,3,4,5,6,7] in 32-bit form in low 8 dwords */\
			"vpslld	$9,%%ymm8,%%ymm8		\n\t"/* The above bytewise offsets need scale *512 to get the needed ones - would include but
											e.g. 1<<8 overflows 1 byte - but x86 ISA only permits scale factors 1,2,4,8, so <<= 8 here. */\
		/* Mask-reg zmm9 = 11...11 - opmask is stupidly zeroed each time we do scatter-store, so need to reinit prior to each invocation */\
		"movl	$-1,%%esi	\n\t"/* Init opmask k1 (Only need the low byte) */\
	/**** a[ 0- 7] = re[ 0, 8,16,24, 4,12,20,28].d0	a[ 8-15] = im[ 0, 8,16,24, 4,12,20,28].d0 : ****/\
		/* Read transposed row-data from the local store: */\
		"vmovaps 		 (%%rcx),%%zmm0					\n\t		vmovaps 0x040(%%rcx),%%zmm10	\n\t"\
		"vmovaps 	0x400(%%rcx),%%zmm1					\n\t		vmovaps 0x440(%%rcx),%%zmm11	\n\t"\
		"vmovaps 	0x800(%%rcx),%%zmm2					\n\t		vmovaps 0x840(%%rcx),%%zmm12	\n\t"\
		"vmovaps 	0xc00(%%rcx),%%zmm3					\n\t		vmovaps 0xc40(%%rcx),%%zmm13	\n\t"\
		"vmovaps 	0x200(%%rcx),%%zmm4					\n\t		vmovaps 0x240(%%rcx),%%zmm14	\n\t"\
		"vmovaps 	0x600(%%rcx),%%zmm5					\n\t		vmovaps 0x640(%%rcx),%%zmm15	\n\t"\
		"vmovaps 	0xa00(%%rcx),%%zmm6					\n\t		vmovaps 0xa40(%%rcx),%%zmm16	\n\t"\
		"vmovaps 	0xe00(%%rcx),%%zmm7					\n\t		vmovaps 0xe40(%%rcx),%%zmm17	\n\t"\
		/* Real parts: */\
		"kmovw	%%esi,%%k1	\n\t	vscatterdpd %%zmm0,0x00(%%rax,%%ymm8)%{%%k1%}	\n\t"\
		"kmovw	%%esi,%%k1	\n\t	vscatterdpd %%zmm1,0x08(%%rax,%%ymm8)%{%%k1%}	\n\t"\
		"kmovw	%%esi,%%k1	\n\t	vscatterdpd %%zmm2,0x10(%%rax,%%ymm8)%{%%k1%}	\n\t"\
		"kmovw	%%esi,%%k1	\n\t	vscatterdpd %%zmm3,0x18(%%rax,%%ymm8)%{%%k1%}	\n\t"\
		"kmovw	%%esi,%%k1	\n\t	vscatterdpd %%zmm4,0x20(%%rax,%%ymm8)%{%%k1%}	\n\t"\
		"kmovw	%%esi,%%k1	\n\t	vscatterdpd %%zmm5,0x28(%%rax,%%ymm8)%{%%k1%}	\n\t"\
		"kmovw	%%esi,%%k1	\n\t	vscatterdpd %%zmm6,0x30(%%rax,%%ymm8)%{%%k1%}	\n\t"\
		"kmovw	%%esi,%%k1	\n\t	vscatterdpd %%zmm7,0x38(%%rax,%%ymm8)%{%%k1%}	\n\t"\
													/* Imag parts: */\
													"kmovw	%%esi,%%k1	\n\t	vscatterdpd %%zmm10,0x40(%%rax,%%ymm8)%{%%k1%}	\n\t"\
													"kmovw	%%esi,%%k1	\n\t	vscatterdpd %%zmm11,0x48(%%rax,%%ymm8)%{%%k1%}	\n\t"\
													"kmovw	%%esi,%%k1	\n\t	vscatterdpd %%zmm12,0x50(%%rax,%%ymm8)%{%%k1%}	\n\t"\
													"kmovw	%%esi,%%k1	\n\t	vscatterdpd %%zmm13,0x58(%%rax,%%ymm8)%{%%k1%}	\n\t"\
													"kmovw	%%esi,%%k1	\n\t	vscatterdpd %%zmm14,0x60(%%rax,%%ymm8)%{%%k1%}	\n\t"\
													"kmovw	%%esi,%%k1	\n\t	vscatterdpd %%zmm15,0x68(%%rax,%%ymm8)%{%%k1%}	\n\t"\
													"kmovw	%%esi,%%k1	\n\t	vscatterdpd %%zmm16,0x70(%%rax,%%ymm8)%{%%k1%}	\n\t"\
													"kmovw	%%esi,%%k1	\n\t	vscatterdpd %%zmm17,0x78(%%rax,%%ymm8)%{%%k1%}	\n\t"\
	/**** a[16-23] = re[ 2,10,18,26, 6,14,22,30].d0	a[24-31] = im[ 2,10,18,26, 6,14,22,30].d0 : ****/\
	"addq	$0x80,%%rax	\n\t"\
	"addq	$0x100,%%rcx	\n\t"\
		/* Read transposed row-data from the local store: */\
		"vmovaps 		 (%%rcx),%%zmm0					\n\t		vmovaps 0x040(%%rcx),%%zmm10	\n\t"\
		"vmovaps 	0x400(%%rcx),%%zmm1					\n\t		vmovaps 0x440(%%rcx),%%zmm11	\n\t"\
		"vmovaps 	0x800(%%rcx),%%zmm2					\n\t		vmovaps 0x840(%%rcx),%%zmm12	\n\t"\
		"vmovaps 	0xc00(%%rcx),%%zmm3					\n\t		vmovaps 0xc40(%%rcx),%%zmm13	\n\t"\
		"vmovaps 	0x200(%%rcx),%%zmm4					\n\t		vmovaps 0x240(%%rcx),%%zmm14	\n\t"\
		"vmovaps 	0x600(%%rcx),%%zmm5					\n\t		vmovaps 0x640(%%rcx),%%zmm15	\n\t"\
		"vmovaps 	0xa00(%%rcx),%%zmm6					\n\t		vmovaps 0xa40(%%rcx),%%zmm16	\n\t"\
		"vmovaps 	0xe00(%%rcx),%%zmm7					\n\t		vmovaps 0xe40(%%rcx),%%zmm17	\n\t"\
		/* Real parts: */\
		"kmovw	%%esi,%%k1	\n\t	vscatterdpd %%zmm0,0x00(%%rax,%%ymm8)%{%%k1%}	\n\t"\
		"kmovw	%%esi,%%k1	\n\t	vscatterdpd %%zmm1,0x08(%%rax,%%ymm8)%{%%k1%}	\n\t"\
		"kmovw	%%esi,%%k1	\n\t	vscatterdpd %%zmm2,0x10(%%rax,%%ymm8)%{%%k1%}	\n\t"\
		"kmovw	%%esi,%%k1	\n\t	vscatterdpd %%zmm3,0x18(%%rax,%%ymm8)%{%%k1%}	\n\t"\
		"kmovw	%%esi,%%k1	\n\t	vscatterdpd %%zmm4,0x20(%%rax,%%ymm8)%{%%k1%}	\n\t"\
		"kmovw	%%esi,%%k1	\n\t	vscatterdpd %%zmm5,0x28(%%rax,%%ymm8)%{%%k1%}	\n\t"\
		"kmovw	%%esi,%%k1	\n\t	vscatterdpd %%zmm6,0x30(%%rax,%%ymm8)%{%%k1%}	\n\t"\
		"kmovw	%%esi,%%k1	\n\t	vscatterdpd %%zmm7,0x38(%%rax,%%ymm8)%{%%k1%}	\n\t"\
													/* Imag parts: */\
													"kmovw	%%esi,%%k1	\n\t	vscatterdpd %%zmm10,0x40(%%rax,%%ymm8)%{%%k1%}	\n\t"\
													"kmovw	%%esi,%%k1	\n\t	vscatterdpd %%zmm11,0x48(%%rax,%%ymm8)%{%%k1%}	\n\t"\
													"kmovw	%%esi,%%k1	\n\t	vscatterdpd %%zmm12,0x50(%%rax,%%ymm8)%{%%k1%}	\n\t"\
													"kmovw	%%esi,%%k1	\n\t	vscatterdpd %%zmm13,0x58(%%rax,%%ymm8)%{%%k1%}	\n\t"\
													"kmovw	%%esi,%%k1	\n\t	vscatterdpd %%zmm14,0x60(%%rax,%%ymm8)%{%%k1%}	\n\t"\
													"kmovw	%%esi,%%k1	\n\t	vscatterdpd %%zmm15,0x68(%%rax,%%ymm8)%{%%k1%}	\n\t"\
													"kmovw	%%esi,%%k1	\n\t	vscatterdpd %%zmm16,0x70(%%rax,%%ymm8)%{%%k1%}	\n\t"\
													"kmovw	%%esi,%%k1	\n\t	vscatterdpd %%zmm17,0x78(%%rax,%%ymm8)%{%%k1%}	\n\t"\
	/**** a[32-39] = re[ 1, 9,17,25, 5,13,21,29].d0	a[40-47] = im[ 1, 9,17,25, 5,13,21,29].d0 : ****/\
	"addq	$0x80,%%rax	\n\t"\
	"subq	$0x80,%%rcx	\n\t"\
		/* Read transposed row-data from the local store: */\
		"vmovaps 		 (%%rcx),%%zmm0					\n\t		vmovaps 0x040(%%rcx),%%zmm10	\n\t"\
		"vmovaps 	0x400(%%rcx),%%zmm1					\n\t		vmovaps 0x440(%%rcx),%%zmm11	\n\t"\
		"vmovaps 	0x800(%%rcx),%%zmm2					\n\t		vmovaps 0x840(%%rcx),%%zmm12	\n\t"\
		"vmovaps 	0xc00(%%rcx),%%zmm3					\n\t		vmovaps 0xc40(%%rcx),%%zmm13	\n\t"\
		"vmovaps 	0x200(%%rcx),%%zmm4					\n\t		vmovaps 0x240(%%rcx),%%zmm14	\n\t"\
		"vmovaps 	0x600(%%rcx),%%zmm5					\n\t		vmovaps 0x640(%%rcx),%%zmm15	\n\t"\
		"vmovaps 	0xa00(%%rcx),%%zmm6					\n\t		vmovaps 0xa40(%%rcx),%%zmm16	\n\t"\
		"vmovaps 	0xe00(%%rcx),%%zmm7					\n\t		vmovaps 0xe40(%%rcx),%%zmm17	\n\t"\
		/* Real parts: */\
		"kmovw	%%esi,%%k1	\n\t	vscatterdpd %%zmm0,0x00(%%rax,%%ymm8)%{%%k1%}	\n\t"\
		"kmovw	%%esi,%%k1	\n\t	vscatterdpd %%zmm1,0x08(%%rax,%%ymm8)%{%%k1%}	\n\t"\
		"kmovw	%%esi,%%k1	\n\t	vscatterdpd %%zmm2,0x10(%%rax,%%ymm8)%{%%k1%}	\n\t"\
		"kmovw	%%esi,%%k1	\n\t	vscatterdpd %%zmm3,0x18(%%rax,%%ymm8)%{%%k1%}	\n\t"\
		"kmovw	%%esi,%%k1	\n\t	vscatterdpd %%zmm4,0x20(%%rax,%%ymm8)%{%%k1%}	\n\t"\
		"kmovw	%%esi,%%k1	\n\t	vscatterdpd %%zmm5,0x28(%%rax,%%ymm8)%{%%k1%}	\n\t"\
		"kmovw	%%esi,%%k1	\n\t	vscatterdpd %%zmm6,0x30(%%rax,%%ymm8)%{%%k1%}	\n\t"\
		"kmovw	%%esi,%%k1	\n\t	vscatterdpd %%zmm7,0x38(%%rax,%%ymm8)%{%%k1%}	\n\t"\
													/* Imag parts: */\
													"kmovw	%%esi,%%k1	\n\t	vscatterdpd %%zmm10,0x40(%%rax,%%ymm8)%{%%k1%}	\n\t"\
													"kmovw	%%esi,%%k1	\n\t	vscatterdpd %%zmm11,0x48(%%rax,%%ymm8)%{%%k1%}	\n\t"\
													"kmovw	%%esi,%%k1	\n\t	vscatterdpd %%zmm12,0x50(%%rax,%%ymm8)%{%%k1%}	\n\t"\
													"kmovw	%%esi,%%k1	\n\t	vscatterdpd %%zmm13,0x58(%%rax,%%ymm8)%{%%k1%}	\n\t"\
													"kmovw	%%esi,%%k1	\n\t	vscatterdpd %%zmm14,0x60(%%rax,%%ymm8)%{%%k1%}	\n\t"\
													"kmovw	%%esi,%%k1	\n\t	vscatterdpd %%zmm15,0x68(%%rax,%%ymm8)%{%%k1%}	\n\t"\
													"kmovw	%%esi,%%k1	\n\t	vscatterdpd %%zmm16,0x70(%%rax,%%ymm8)%{%%k1%}	\n\t"\
													"kmovw	%%esi,%%k1	\n\t	vscatterdpd %%zmm17,0x78(%%rax,%%ymm8)%{%%k1%}	\n\t"\
	/**** a[48-55] = re[ 3,11,19,27, 7,15,23,31].d0	a[56-63] = im[ 3,11,19,27, 7,15,23,31].d0 : ****/\
	"addq	$0x80,%%rax	\n\t"\
	"addq	$0x100,%%rcx	\n\t"\
		/* Read transposed row-data from the local store: */\
		"vmovaps 		 (%%rcx),%%zmm0					\n\t		vmovaps 0x040(%%rcx),%%zmm10	\n\t"\
		"vmovaps 	0x400(%%rcx),%%zmm1					\n\t		vmovaps 0x440(%%rcx),%%zmm11	\n\t"\
		"vmovaps 	0x800(%%rcx),%%zmm2					\n\t		vmovaps 0x840(%%rcx),%%zmm12	\n\t"\
		"vmovaps 	0xc00(%%rcx),%%zmm3					\n\t		vmovaps 0xc40(%%rcx),%%zmm13	\n\t"\
		"vmovaps 	0x200(%%rcx),%%zmm4					\n\t		vmovaps 0x240(%%rcx),%%zmm14	\n\t"\
		"vmovaps 	0x600(%%rcx),%%zmm5					\n\t		vmovaps 0x640(%%rcx),%%zmm15	\n\t"\
		"vmovaps 	0xa00(%%rcx),%%zmm6					\n\t		vmovaps 0xa40(%%rcx),%%zmm16	\n\t"\
		"vmovaps 	0xe00(%%rcx),%%zmm7					\n\t		vmovaps 0xe40(%%rcx),%%zmm17	\n\t"\
		/* Real parts: */\
		"kmovw	%%esi,%%k1	\n\t	vscatterdpd %%zmm0,0x00(%%rax,%%ymm8)%{%%k1%}	\n\t"\
		"kmovw	%%esi,%%k1	\n\t	vscatterdpd %%zmm1,0x08(%%rax,%%ymm8)%{%%k1%}	\n\t"\
		"kmovw	%%esi,%%k1	\n\t	vscatterdpd %%zmm2,0x10(%%rax,%%ymm8)%{%%k1%}	\n\t"\
		"kmovw	%%esi,%%k1	\n\t	vscatterdpd %%zmm3,0x18(%%rax,%%ymm8)%{%%k1%}	\n\t"\
		"kmovw	%%esi,%%k1	\n\t	vscatterdpd %%zmm4,0x20(%%rax,%%ymm8)%{%%k1%}	\n\t"\
		"kmovw	%%esi,%%k1	\n\t	vscatterdpd %%zmm5,0x28(%%rax,%%ymm8)%{%%k1%}	\n\t"\
		"kmovw	%%esi,%%k1	\n\t	vscatterdpd %%zmm6,0x30(%%rax,%%ymm8)%{%%k1%}	\n\t"\
		"kmovw	%%esi,%%k1	\n\t	vscatterdpd %%zmm7,0x38(%%rax,%%ymm8)%{%%k1%}	\n\t"\
													/* Imag parts: */\
													"kmovw	%%esi,%%k1	\n\t	vscatterdpd %%zmm10,0x40(%%rax,%%ymm8)%{%%k1%}	\n\t"\
													"kmovw	%%esi,%%k1	\n\t	vscatterdpd %%zmm11,0x48(%%rax,%%ymm8)%{%%k1%}	\n\t"\
													"kmovw	%%esi,%%k1	\n\t	vscatterdpd %%zmm12,0x50(%%rax,%%ymm8)%{%%k1%}	\n\t"\
													"kmovw	%%esi,%%k1	\n\t	vscatterdpd %%zmm13,0x58(%%rax,%%ymm8)%{%%k1%}	\n\t"\
													"kmovw	%%esi,%%k1	\n\t	vscatterdpd %%zmm14,0x60(%%rax,%%ymm8)%{%%k1%}	\n\t"\
													"kmovw	%%esi,%%k1	\n\t	vscatterdpd %%zmm15,0x68(%%rax,%%ymm8)%{%%k1%}	\n\t"\
													"kmovw	%%esi,%%k1	\n\t	vscatterdpd %%zmm16,0x70(%%rax,%%ymm8)%{%%k1%}	\n\t"\
													"kmovw	%%esi,%%k1	\n\t	vscatterdpd %%zmm17,0x78(%%rax,%%ymm8)%{%%k1%}	\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__add1] "m" (Xadd1)\
		 ,[__add2] "m" (Xadd2)\
		 ,[__add3] "m" (Xadd3)\
		 ,[__add4] "m" (Xadd4)\
		 ,[__add5] "m" (Xadd5)\
		 ,[__add6] "m" (Xadd6)\
		 ,[__add7] "m" (Xadd7)\
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__r00] "m" (Xr00)\
		 ,[__r08] "m" (Xr08)\
		 ,[__r10] "m" (Xr10)\
		 ,[__r20] "m" (Xr20)\
		 ,[__r28] "m" (Xr28)\
		 ,[__r30] "m" (Xr30)\
		 ,[__c01] "m" (Xc01)\
		 ,[__c02] "m" (Xc02)\
		 ,[__c04] "m" (Xc04)\
		 ,[__c06] "m" (Xc06)\
		 ,[__c0A] "m" (Xc0A)\
		 ,[__c10] "m" (Xc10)\
		 ,[__c12] "m" (Xc12)\
		 ,[__c1A] "m" (Xc1A)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r10","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"		/* Clobbered registers */\
	);\
	}

  #else	// USE_16_REG = False, i.e. use default 32-SIMD-reg version of macro:

   #ifdef USE_IMCI512	// 1st-gen Xeon Phi - Use modified 8x8 doubles-transpose algo [1a] from util.c:test_simd_transpose_8x8()

	#define SSE2_RADIX32_WRAPPER_DIT(Xadd0,Xadd1,Xadd2,Xadd3,Xadd4,Xadd5,Xadd6,Xadd7,Xisrt2,Xr00,Xr08,Xr10,Xr20,Xr28,Xr30,Xc01,Xc02,Xc04,Xc06,Xc0A,Xc10,Xc12,Xc1A)\
	{\
	__asm__ volatile (\
	/* Cf. IMCI512 8x8 doubles-transpose algo [1a] in util.c ...
	do mask-register-setting first, so as to not clobber any GPRs used by compute section: */\
		"movl $0b10101010,%%ebx	\n\t movl $0b11001100,%%ecx	\n\t movl $0b11110000,%%edx	\n\t"\
		"kmov	%%ebx,%%k1		\n\t kmov	%%ecx,%%k3		\n\t kmov	%%edx,%%k5		\n\t"\
		"knot	%%k1 ,%%k2		\n\t knot	%%k3 ,%%k4		\n\t"\
	/************************************************************************/\
	/*...And do an inverse DIT radix-32 pass on the squared-data blocks.	*/\
	/************************************************************************/\
	/* Block 1: */\
		"movq	%[__isrt2],%%rsi			\n\t	leaq	0x11c0(%%rsi),%%rdi		\n\t"/* two */\
		"movq		%[__r00],%%rax			\n\t	vmovaps		 (%%rdi),%%zmm29	\n\t"/* 2.0*/	/* Next pair of 4-DFTs - Place these results in regs whose index is += 16 w.r.to the above 4-DFTs: */\
	/**** SSE2_RADIX4_DIT_IN_PLACE() ****/	/* Block 2 has addr-offsets +0x100 w.r.to Block 1: */	/* Byte-offsets all +0x200 larger in this 4-DFT pair */\
		"vmovaps	0x800(%%rax),%%zmm0		\n\t	vmovaps		0x900(%%rax),%%zmm8 		\n\t	vmovaps	0xa00(%%rax),%%zmm16		\n\t	vmovaps		0xb00(%%rax),%%zmm24		\n\t"\
		"vmovaps	0x840(%%rax),%%zmm1		\n\t	vmovaps		0x940(%%rax),%%zmm9 		\n\t	vmovaps	0xa40(%%rax),%%zmm17		\n\t	vmovaps		0xb40(%%rax),%%zmm25		\n\t"\
		"vmovaps		 (%%rax),%%zmm2		\n\t	vmovaps		0x100(%%rax),%%zmm10		\n\t	vmovaps	0x200(%%rax),%%zmm18		\n\t	vmovaps		0x300(%%rax),%%zmm26		\n\t"\
		"vmovaps	0x040(%%rax),%%zmm3		\n\t	vmovaps		0x140(%%rax),%%zmm11		\n\t	vmovaps	0x240(%%rax),%%zmm19		\n\t	vmovaps		0x340(%%rax),%%zmm27		\n\t"\
		"vsubpd		%%zmm0,%%zmm2,%%zmm2	\n\t	vsubpd		%%zmm8 ,%%zmm10,%%zmm10		\n\t	vsubpd		%%zmm16,%%zmm18,%%zmm18	\n\t	vsubpd		%%zmm24,%%zmm26,%%zmm26		\n\t"\
		"vsubpd		%%zmm1,%%zmm3,%%zmm3	\n\t	vsubpd		%%zmm9 ,%%zmm11,%%zmm11		\n\t	vsubpd		%%zmm17,%%zmm19,%%zmm19	\n\t	vsubpd		%%zmm25,%%zmm27,%%zmm27		\n\t"\
		"vmovaps	0xc00(%%rax),%%zmm4		\n\t	vmovaps		0xd00(%%rax),%%zmm12		\n\t	vmovaps	0xe00(%%rax),%%zmm20		\n\t	vmovaps		0xf00(%%rax),%%zmm28		\n\t"\
		"vmovaps	0xc40(%%rax),%%zmm5		\n\t	vmovaps		0xd40(%%rax),%%zmm13		\n\t	vmovaps	0xe40(%%rax),%%zmm21		\n\t /* vmovaps		0xf40(%%rax),%%zmm29	*/	\n\t"\
		"vmovaps	0x400(%%rax),%%zmm6		\n\t	vmovaps		0x500(%%rax),%%zmm14		\n\t	vmovaps	0x600(%%rax),%%zmm22		\n\t	vmovaps		0x700(%%rax),%%zmm30		\n\t"\
		"vmovaps	0x440(%%rax),%%zmm7		\n\t	vmovaps		0x540(%%rax),%%zmm15		\n\t	vmovaps	0x640(%%rax),%%zmm23		\n\t	vmovaps		0x740(%%rax),%%zmm31		\n\t"\
		"vsubpd		%%zmm4,%%zmm6,%%zmm6	\n\t	vsubpd		%%zmm12,%%zmm14,%%zmm14		\n\t	vsubpd		%%zmm20,%%zmm22,%%zmm22	\n\t	vsubpd		%%zmm28,%%zmm30,%%zmm30		\n\t"\
		"vsubpd		%%zmm5,%%zmm7,%%zmm7	\n\t	vsubpd		%%zmm13,%%zmm15,%%zmm15		\n\t	vsubpd		%%zmm21,%%zmm23,%%zmm23	\n\t	vsubpd 0xf40(%%rax),%%zmm31,%%zmm31		\n\t"\
	"vfmadd132pd	%%zmm29,%%zmm2,%%zmm0	\n\t vfmadd132pd	%%zmm29,%%zmm10,%%zmm8 		\n\t vfmadd132pd	%%zmm29,%%zmm18,%%zmm16	\n\t vfmadd132pd	%%zmm29,%%zmm26,%%zmm24		\n\t"\
	"vfmadd132pd	%%zmm29,%%zmm3,%%zmm1	\n\t vfmadd132pd	%%zmm29,%%zmm11,%%zmm9 		\n\t vfmadd132pd	%%zmm29,%%zmm19,%%zmm17	\n\t vfmadd132pd	%%zmm29,%%zmm27,%%zmm25		\n\t"\
	"vfmadd132pd	%%zmm29,%%zmm6,%%zmm4	\n\t vfmadd132pd	%%zmm29,%%zmm14,%%zmm12		\n\t vfmadd132pd	%%zmm29,%%zmm22,%%zmm20	\n\t vfmadd132pd	%%zmm29,%%zmm30,%%zmm28		\n\t"\
	"vfmadd132pd	%%zmm29,%%zmm7,%%zmm5	\n\t vfmadd132pd	%%zmm29,%%zmm15,%%zmm13		\n\t vfmadd132pd	%%zmm29,%%zmm23,%%zmm21	\n\t vfmadd132pd 0xf40(%%rax),%%zmm31,%%zmm29	\n\t"\
		"vsubpd		%%zmm4 ,%%zmm0,%%zmm0	\n\t	vsubpd		%%zmm12,%%zmm8 ,%%zmm8 		\n\t	vsubpd		%%zmm20,%%zmm16,%%zmm16	\n\t	vsubpd		%%zmm28,%%zmm24,%%zmm24		\n\t"\
		"vsubpd		%%zmm5 ,%%zmm1,%%zmm1	\n\t	vsubpd		%%zmm13,%%zmm9 ,%%zmm9 		\n\t	vsubpd		%%zmm21,%%zmm17,%%zmm17	\n\t	vsubpd		%%zmm29,%%zmm25,%%zmm25		\n\t"\
		"vsubpd		%%zmm7 ,%%zmm2,%%zmm2	\n\t	vsubpd		%%zmm15,%%zmm10,%%zmm10		\n\t	vsubpd		%%zmm23,%%zmm18,%%zmm18	\n\t	vsubpd		%%zmm31,%%zmm26,%%zmm26		\n\t"\
		"vsubpd		%%zmm6 ,%%zmm3,%%zmm3	\n\t	vsubpd		%%zmm14,%%zmm11,%%zmm11		\n\t	vsubpd		%%zmm22,%%zmm19,%%zmm19	\n\t	vsubpd		%%zmm30,%%zmm27,%%zmm27		\n\t"\
	"vmovaps %%zmm30,0xf40(%%rax)\n\t"/* spill zmm30 to make room for two */"vmovaps (%%rdi),%%zmm30\n\t"\
	"vfmadd132pd	%%zmm30,%%zmm0,%%zmm4	\n\t vfmadd132pd	%%zmm30,%%zmm8 ,%%zmm12		\n\t vfmadd132pd	%%zmm30,%%zmm16,%%zmm20	\n\t vfmadd132pd	%%zmm30,%%zmm24,%%zmm28		\n\t"\
	"vfmadd132pd	%%zmm30,%%zmm1,%%zmm5	\n\t vfmadd132pd	%%zmm30,%%zmm9 ,%%zmm13		\n\t vfmadd132pd	%%zmm30,%%zmm17,%%zmm21	\n\t vfmadd132pd	%%zmm30,%%zmm25,%%zmm29		\n\t"\
	"vfmadd132pd	%%zmm30,%%zmm2,%%zmm7	\n\t vfmadd132pd	%%zmm30,%%zmm10,%%zmm15		\n\t vfmadd132pd	%%zmm30,%%zmm18,%%zmm23	\n\t vfmadd132pd	%%zmm30,%%zmm26,%%zmm31		\n\t"\
	"vfmadd132pd	%%zmm30,%%zmm3,%%zmm6	\n\t vfmadd132pd	%%zmm30,%%zmm11,%%zmm14		\n\t vfmadd132pd	%%zmm30,%%zmm19,%%zmm22	\n\t vfmadd132pd 0xf40(%%rax),%%zmm27,%%zmm30	\n\t"\
	/**** SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS ****//**** SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS ****/	/* These are actually the final steps of the col3/4 4-DFTs, can do alongside */\
	/* 0x200(r00,r10,r20,r30,r08,r18,r28,r38): *//* 0x200(r04,r14,r24,r34,r0C,r1C,r2C,r3C): */	/* first steps of the combine-subs in cols 1/2 due to nonoverlapping registers: */\
		"vsubpd		%%zmm20,%%zmm4,%%zmm4	\n\t	vsubpd		%%zmm28,%%zmm12,%%zmm12		\n\t	vaddpd		%%zmm23,%%zmm19,%%zmm23	\n\t	vaddpd		%%zmm31,%%zmm27,%%zmm31		\n\t"\
		"vsubpd		%%zmm17,%%zmm0,%%zmm0	\n\t	vsubpd		%%zmm25,%%zmm8 ,%%zmm8 		\n\t	vaddpd		%%zmm18,%%zmm22,%%zmm18	\n\t	vaddpd		%%zmm26,%%zmm30,%%zmm26		\n\t"\
		"vsubpd		%%zmm21,%%zmm5,%%zmm5	\n\t	vsubpd		%%zmm29,%%zmm13,%%zmm13		\n\t vfmsub132pd	(%%rdi),%%zmm23,%%zmm19	\n\t vfmsub132pd	(%%rdi),%%zmm31,%%zmm27		\n\t"\
		"vsubpd		%%zmm16,%%zmm1,%%zmm1	\n\t	vsubpd		%%zmm24,%%zmm9 ,%%zmm9 		\n\t vfmsub132pd	(%%rdi),%%zmm18,%%zmm22	\n\t vfmsub132pd	(%%rdi),%%zmm26,%%zmm30		\n\t"\
		"vmovaps	%%zmm4 , 0x200(%%rax)	\n\t	vmovaps		%%zmm12, 0x300(%%rax)		\n\t"\
		"vmovaps	%%zmm0 , 0xa00(%%rax)	\n\t	vmovaps		%%zmm8 , 0xb00(%%rax)		\n\t"\
		"vmovaps	%%zmm5 , 0x240(%%rax)	\n\t	vmovaps		%%zmm13, 0x340(%%rax)		\n\t"\
		"vmovaps	%%zmm1 , 0x840(%%rax)	\n\t	vmovaps		%%zmm9 , 0x940(%%rax)		\n\t"\
	"vmovaps %%zmm24,0xb40(%%rax)\n\t"/* spill zmm24to make room for two */"vmovaps (%%rdi),%%zmm24\n\t"\
	"vfmadd132pd	%%zmm24,%%zmm4,%%zmm20	\n\t vfmadd132pd	 %%zmm24,%%zmm12,%%zmm28	\n\t"\
	"vfmadd132pd	%%zmm24,%%zmm0,%%zmm17	\n\t vfmadd132pd	 %%zmm24,%%zmm8 ,%%zmm25	\n\t"\
	"vfmadd132pd	%%zmm24,%%zmm5,%%zmm21	\n\t vfmadd132pd	 %%zmm24,%%zmm13,%%zmm29	\n\t"\
	"vfmadd132pd	%%zmm24,%%zmm1,%%zmm16	\n\t vfmadd132pd 0xb40(%%rax),%%zmm9 ,%%zmm24	\n\t"\
		"vmovaps	%%zmm20,      (%%rax)	\n\t	vmovaps		%%zmm28, 0x100(%%rax)		\n\t"\
		"vmovaps	%%zmm17, 0x800(%%rax)	\n\t	vmovaps		%%zmm25, 0x900(%%rax)		\n\t"\
		"vmovaps	%%zmm21, 0x040(%%rax)	\n\t	vmovaps		%%zmm29, 0x140(%%rax)		\n\t"\
		"vmovaps	%%zmm16, 0xa40(%%rax)	\n\t	vmovaps		%%zmm24, 0xb40(%%rax)		\n\t"\
	"vmovaps %%zmm14,0xf40(%%rax)\n\t"/* spill zmm14to make room for isrt2 */"vmovaps (%%rsi),%%zmm14\n\t"\
	"vfnmadd231pd	%%zmm23,%%zmm14,%%zmm7	\n\t vfnmadd231pd	%%zmm31,%%zmm14,%%zmm15		\n\t"\
	"vfnmadd231pd	%%zmm22,%%zmm14,%%zmm2	\n\t vfnmadd231pd	%%zmm30,%%zmm14,%%zmm10		\n\t"\
	"vfnmadd231pd	%%zmm19,%%zmm14,%%zmm3	\n\t vfnmadd231pd	%%zmm27,%%zmm14,%%zmm11		\n\t"\
	"vfnmadd231pd	%%zmm18,%%zmm14,%%zmm6	\n\t vfnmadd213pd 0xf40(%%rax),%%zmm26,%%zmm14	\n\t"\
		"vmovaps	%%zmm7 , 0x600(%%rax)	\n\t	vmovaps		%%zmm15, 0x700(%%rax)		\n\t"\
		"vmovaps	%%zmm2 , 0xe00(%%rax)	\n\t	vmovaps		%%zmm10, 0xf00(%%rax)		\n\t"\
		"vmovaps	%%zmm3 , 0x640(%%rax)	\n\t	vmovaps		%%zmm11, 0x740(%%rax)		\n\t"\
		"vmovaps	%%zmm6 , 0xc40(%%rax)	\n\t	vmovaps		%%zmm14, 0xd40(%%rax)		\n\t"\
	"vmovaps %%zmm26,0xf40(%%rax)\n\t"/* spill zmm26to make room for sqrt2 */"vmovaps 0x80(%%rdi),%%zmm26\n\t"\
	"vfmadd132pd	%%zmm26,%%zmm7,%%zmm23	\n\t vfmadd132pd	 %%zmm26,%%zmm15,%%zmm31	\n\t"\
	"vfmadd132pd	%%zmm26,%%zmm2,%%zmm22	\n\t vfmadd132pd	 %%zmm26,%%zmm10,%%zmm30	\n\t"\
	"vfmadd132pd	%%zmm26,%%zmm3,%%zmm19	\n\t vfmadd132pd	 %%zmm26,%%zmm11,%%zmm27	\n\t"\
	"vfmadd132pd	%%zmm26,%%zmm6,%%zmm18	\n\t vfmadd132pd 0xf40(%%rax),%%zmm14,%%zmm26	\n\t"\
		"vmovaps	%%zmm23, 0x400(%%rax)	\n\t	vmovaps		%%zmm31, 0x500(%%rax)		\n\t"\
		"vmovaps	%%zmm22, 0xc00(%%rax)	\n\t	vmovaps		%%zmm30, 0xd00(%%rax)		\n\t"\
		"vmovaps	%%zmm19, 0x440(%%rax)	\n\t	vmovaps		%%zmm27, 0x540(%%rax)		\n\t"\
		"vmovaps	%%zmm18, 0xe40(%%rax)	\n\t	vmovaps		%%zmm26, 0xf40(%%rax)		\n\t"\
	/*****	SSE2_RADIX4_DIT_IN_PLACE()	*****/	/* Rcol has tmp-addr offset +0x100 w.r.to lcol:  */\
	/* Blocks 3,4 have tmp-addresses offset +0x80 w.r.to Blocks 1,2, respectively. */\
		"											vmovaps		 (%%rdi),%%zmm29	\n\t"/* 2.0*/\
		"vmovaps	0x880(%%rax),%%zmm0		\n\t	vmovaps		0x980(%%rax),%%zmm8 		\n\t	vmovaps	0xa80(%%rax),%%zmm16		\n\t	vmovaps		0xb80(%%rax),%%zmm24		\n\t"\
		"vmovaps	0x8c0(%%rax),%%zmm1		\n\t	vmovaps		0x9c0(%%rax),%%zmm9 		\n\t	vmovaps	0xac0(%%rax),%%zmm17		\n\t	vmovaps		0xbc0(%%rax),%%zmm25		\n\t"\
		"vmovaps	0x080(%%rax),%%zmm2		\n\t	vmovaps		0x180(%%rax),%%zmm10		\n\t	vmovaps	0x280(%%rax),%%zmm18		\n\t	vmovaps		0x380(%%rax),%%zmm26		\n\t"\
		"vmovaps	0x0c0(%%rax),%%zmm3		\n\t	vmovaps		0x1c0(%%rax),%%zmm11		\n\t	vmovaps	0x2c0(%%rax),%%zmm19		\n\t	vmovaps		0x3c0(%%rax),%%zmm27		\n\t"\
		"vsubpd		%%zmm0,%%zmm2,%%zmm2	\n\t	vsubpd		%%zmm8 ,%%zmm10,%%zmm10		\n\t	vsubpd		%%zmm16,%%zmm18,%%zmm18	\n\t	vsubpd		%%zmm24,%%zmm26,%%zmm26		\n\t"\
		"vsubpd		%%zmm1,%%zmm3,%%zmm3	\n\t	vsubpd		%%zmm9 ,%%zmm11,%%zmm11		\n\t	vsubpd		%%zmm17,%%zmm19,%%zmm19	\n\t	vsubpd		%%zmm25,%%zmm27,%%zmm27		\n\t"\
		"vmovaps	0xc80(%%rax),%%zmm4		\n\t	vmovaps		0xd80(%%rax),%%zmm12		\n\t	vmovaps	0xe80(%%rax),%%zmm20		\n\t	vmovaps		0xf80(%%rax),%%zmm28		\n\t"\
		"vmovaps	0xcc0(%%rax),%%zmm5		\n\t	vmovaps		0xdc0(%%rax),%%zmm13		\n\t	vmovaps	0xec0(%%rax),%%zmm21		\n\t /* vmovaps		0xfc0(%%rax),%%zmm29	*/	\n\t"\
		"vmovaps	0x480(%%rax),%%zmm6		\n\t	vmovaps		0x580(%%rax),%%zmm14		\n\t	vmovaps	0x680(%%rax),%%zmm22		\n\t	vmovaps		0x780(%%rax),%%zmm30		\n\t"\
		"vmovaps	0x4c0(%%rax),%%zmm7		\n\t	vmovaps		0x5c0(%%rax),%%zmm15		\n\t	vmovaps	0x6c0(%%rax),%%zmm23		\n\t	vmovaps		0x7c0(%%rax),%%zmm31		\n\t"\
		"vsubpd		%%zmm4,%%zmm6,%%zmm6	\n\t	vsubpd		%%zmm12,%%zmm14,%%zmm14		\n\t	vsubpd		%%zmm20,%%zmm22,%%zmm22	\n\t	vsubpd		%%zmm28,%%zmm30,%%zmm30		\n\t"\
		"vsubpd		%%zmm5,%%zmm7,%%zmm7	\n\t	vsubpd		%%zmm13,%%zmm15,%%zmm15		\n\t	vsubpd		%%zmm21,%%zmm23,%%zmm23	\n\t	vsubpd 0xfc0(%%rax),%%zmm31,%%zmm31		\n\t"\
	"vfmadd132pd	%%zmm29,%%zmm2,%%zmm0	\n\t vfmadd132pd	%%zmm29,%%zmm10,%%zmm8 		\n\t vfmadd132pd	%%zmm29,%%zmm18,%%zmm16	\n\t vfmadd132pd	%%zmm29,%%zmm26,%%zmm24		\n\t"\
	"vfmadd132pd	%%zmm29,%%zmm3,%%zmm1	\n\t vfmadd132pd	%%zmm29,%%zmm11,%%zmm9 		\n\t vfmadd132pd	%%zmm29,%%zmm19,%%zmm17	\n\t vfmadd132pd	%%zmm29,%%zmm27,%%zmm25		\n\t"\
	"vfmadd132pd	%%zmm29,%%zmm6,%%zmm4	\n\t vfmadd132pd	%%zmm29,%%zmm14,%%zmm12		\n\t vfmadd132pd	%%zmm29,%%zmm22,%%zmm20	\n\t vfmadd132pd	%%zmm29,%%zmm30,%%zmm28		\n\t"\
	"vfmadd132pd	%%zmm29,%%zmm7,%%zmm5	\n\t vfmadd132pd	%%zmm29,%%zmm15,%%zmm13		\n\t vfmadd132pd	%%zmm29,%%zmm23,%%zmm21	\n\t vfmadd132pd 0xfc0(%%rax),%%zmm31,%%zmm29	\n\t"\
		"vsubpd		%%zmm4 ,%%zmm0,%%zmm0	\n\t	vsubpd		%%zmm12,%%zmm8 ,%%zmm8 		\n\t	vsubpd		%%zmm20,%%zmm16,%%zmm16	\n\t	vsubpd		%%zmm28,%%zmm24,%%zmm24		\n\t"\
		"vsubpd		%%zmm5 ,%%zmm1,%%zmm1	\n\t	vsubpd		%%zmm13,%%zmm9 ,%%zmm9 		\n\t	vsubpd		%%zmm21,%%zmm17,%%zmm17	\n\t	vsubpd		%%zmm29,%%zmm25,%%zmm25		\n\t"\
		"vsubpd		%%zmm7 ,%%zmm2,%%zmm2	\n\t	vsubpd		%%zmm15,%%zmm10,%%zmm10		\n\t	vsubpd		%%zmm23,%%zmm18,%%zmm18	\n\t	vsubpd		%%zmm31,%%zmm26,%%zmm26		\n\t"\
		"vsubpd		%%zmm6 ,%%zmm3,%%zmm3	\n\t	vsubpd		%%zmm14,%%zmm11,%%zmm11		\n\t	vsubpd		%%zmm22,%%zmm19,%%zmm19	\n\t	vsubpd		%%zmm30,%%zmm27,%%zmm27		\n\t"\
	"vmovaps %%zmm30,0xfc0(%%rax)\n\t"/* spill zmm30 to make room for two */"vmovaps (%%rdi),%%zmm30\n\t"\
	"vfmadd132pd	%%zmm30,%%zmm0,%%zmm4	\n\t vfmadd132pd	%%zmm30,%%zmm8 ,%%zmm12		\n\t vfmadd132pd	%%zmm30,%%zmm16,%%zmm20	\n\t vfmadd132pd	%%zmm30,%%zmm24,%%zmm28		\n\t"\
	"vfmadd132pd	%%zmm30,%%zmm1,%%zmm5	\n\t vfmadd132pd	%%zmm30,%%zmm9 ,%%zmm13		\n\t vfmadd132pd	%%zmm30,%%zmm17,%%zmm21	\n\t vfmadd132pd	%%zmm30,%%zmm25,%%zmm29		\n\t"\
	"vfmadd132pd	%%zmm30,%%zmm2,%%zmm7	\n\t vfmadd132pd	%%zmm30,%%zmm10,%%zmm15		\n\t vfmadd132pd	%%zmm30,%%zmm18,%%zmm23	\n\t vfmadd132pd	%%zmm30,%%zmm26,%%zmm31		\n\t"\
	"vfmadd132pd	%%zmm30,%%zmm3,%%zmm6	\n\t vfmadd132pd	%%zmm30,%%zmm11,%%zmm14		\n\t vfmadd132pd	%%zmm30,%%zmm19,%%zmm22	\n\t vfmadd132pd 0xfc0(%%rax),%%zmm27,%%zmm30	\n\t"\
	/**** SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS ****//**** SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS ****/	/* These are actually the final steps of the col3/4 4-DFTs, can do alongside */\
	/* 0x200(r00,r10,r20,r30,r08,r18,r28,r38): *//* 0x200(r04,r14,r24,r34,r0C,r1C,r2C,r3C): */	/* first steps of the combine-subs in cols 1/2 due to nonoverlapping registers: */\
		"vsubpd		%%zmm20,%%zmm4,%%zmm4	\n\t	vsubpd		%%zmm28,%%zmm12,%%zmm12		\n\t	vaddpd		%%zmm23,%%zmm19,%%zmm23	\n\t	vaddpd		%%zmm31,%%zmm27,%%zmm31		\n\t"\
		"vsubpd		%%zmm17,%%zmm0,%%zmm0	\n\t	vsubpd		%%zmm25,%%zmm8 ,%%zmm8 		\n\t	vaddpd		%%zmm18,%%zmm22,%%zmm18	\n\t	vaddpd		%%zmm26,%%zmm30,%%zmm26		\n\t"\
		"vsubpd		%%zmm21,%%zmm5,%%zmm5	\n\t	vsubpd		%%zmm29,%%zmm13,%%zmm13		\n\t vfmsub132pd	(%%rdi),%%zmm23,%%zmm19	\n\t vfmsub132pd	(%%rdi),%%zmm31,%%zmm27		\n\t"\
		"vsubpd		%%zmm16,%%zmm1,%%zmm1	\n\t	vsubpd		%%zmm24,%%zmm9 ,%%zmm9 		\n\t vfmsub132pd	(%%rdi),%%zmm18,%%zmm22	\n\t vfmsub132pd	(%%rdi),%%zmm26,%%zmm30		\n\t"\
		"vmovaps	%%zmm4 , 0x280(%%rax)	\n\t	vmovaps		%%zmm12, 0x380(%%rax)		\n\t"\
		"vmovaps	%%zmm0 , 0xa80(%%rax)	\n\t	vmovaps		%%zmm8 , 0xb80(%%rax)		\n\t"\
		"vmovaps	%%zmm5 , 0x2c0(%%rax)	\n\t	vmovaps		%%zmm13, 0x3c0(%%rax)		\n\t"\
		"vmovaps	%%zmm1 , 0x8c0(%%rax)	\n\t	vmovaps		%%zmm9 , 0x9c0(%%rax)		\n\t"\
	"vmovaps %%zmm24,0xbc0(%%rax)\n\t"/* spill zmm24to make room for two */"vmovaps (%%rdi),%%zmm24\n\t"\
	"vfmadd132pd	%%zmm24,%%zmm4,%%zmm20	\n\t vfmadd132pd	 %%zmm24,%%zmm12,%%zmm28	\n\t"\
	"vfmadd132pd	%%zmm24,%%zmm0,%%zmm17	\n\t vfmadd132pd	 %%zmm24,%%zmm8 ,%%zmm25	\n\t"\
	"vfmadd132pd	%%zmm24,%%zmm5,%%zmm21	\n\t vfmadd132pd	 %%zmm24,%%zmm13,%%zmm29	\n\t"\
	"vfmadd132pd	%%zmm24,%%zmm1,%%zmm16	\n\t vfmadd132pd 0xbc0(%%rax),%%zmm9 ,%%zmm24	\n\t"\
		"vmovaps	%%zmm20, 0x080(%%rax)	\n\t	vmovaps		%%zmm28, 0x180(%%rax)		\n\t"\
		"vmovaps	%%zmm17, 0x880(%%rax)	\n\t	vmovaps		%%zmm25, 0x980(%%rax)		\n\t"\
		"vmovaps	%%zmm21, 0x0c0(%%rax)	\n\t	vmovaps		%%zmm29, 0x1c0(%%rax)		\n\t"\
		"vmovaps	%%zmm16, 0xac0(%%rax)	\n\t	vmovaps		%%zmm24, 0xbc0(%%rax)		\n\t"\
	"vmovaps %%zmm14,0xfc0(%%rax)\n\t"/* spill zmm14to make room for isrt2 */"vmovaps (%%rsi),%%zmm14\n\t"\
	"vfnmadd231pd	%%zmm23,%%zmm14,%%zmm7	\n\t vfnmadd231pd	%%zmm31,%%zmm14,%%zmm15		\n\t"\
	"vfnmadd231pd	%%zmm22,%%zmm14,%%zmm2	\n\t vfnmadd231pd	%%zmm30,%%zmm14,%%zmm10		\n\t"\
	"vfnmadd231pd	%%zmm19,%%zmm14,%%zmm3	\n\t vfnmadd231pd	%%zmm27,%%zmm14,%%zmm11		\n\t"\
	"vfnmadd231pd	%%zmm18,%%zmm14,%%zmm6	\n\t vfnmadd213pd 0xfc0(%%rax),%%zmm26,%%zmm14	\n\t"\
		"vmovaps	%%zmm7 , 0x680(%%rax)	\n\t	vmovaps		%%zmm15, 0x780(%%rax)		\n\t"\
		"vmovaps	%%zmm2 , 0xe80(%%rax)	\n\t	vmovaps		%%zmm10, 0xf80(%%rax)		\n\t"\
		"vmovaps	%%zmm3 , 0x6c0(%%rax)	\n\t	vmovaps		%%zmm11, 0x7c0(%%rax)		\n\t"\
		"vmovaps	%%zmm6 , 0xcc0(%%rax)	\n\t	vmovaps		%%zmm14, 0xdc0(%%rax)		\n\t"\
	"vmovaps %%zmm26,0xfc0(%%rax)\n\t"/* spill zmm26to make room for sqrt2 */"vmovaps 0x80(%%rdi),%%zmm26\n\t"\
	"vfmadd132pd	%%zmm26,%%zmm7,%%zmm23	\n\t vfmadd132pd	 %%zmm26,%%zmm15,%%zmm31	\n\t"\
	"vfmadd132pd	%%zmm26,%%zmm2,%%zmm22	\n\t vfmadd132pd	 %%zmm26,%%zmm10,%%zmm30	\n\t"\
	"vfmadd132pd	%%zmm26,%%zmm3,%%zmm19	\n\t vfmadd132pd	 %%zmm26,%%zmm11,%%zmm27	\n\t"\
	"vfmadd132pd	%%zmm26,%%zmm6,%%zmm18	\n\t vfmadd132pd 0xfc0(%%rax),%%zmm14,%%zmm26	\n\t"\
		"vmovaps	%%zmm23, 0x480(%%rax)	\n\t	vmovaps		%%zmm31, 0x580(%%rax)		\n\t"\
		"vmovaps	%%zmm22, 0xc80(%%rax)	\n\t	vmovaps		%%zmm30, 0xd80(%%rax)		\n\t"\
		"vmovaps	%%zmm19, 0x4c0(%%rax)	\n\t	vmovaps		%%zmm27, 0x5c0(%%rax)		\n\t"\
		"vmovaps	%%zmm18, 0xec0(%%rax)	\n\t	vmovaps		%%zmm26, 0xfc0(%%rax)		\n\t"\
		"\n\t"\
	/***************************************************************************************/\
	/* Now do eight more radix-4 transforms, including the internal and external twiddles: */\
	/***************************************************************************************/\
		"\n\t"\
	/*
	Using upper block(s) of main array for temp-storage in section below led to a nasty AVX bug to track down:
	In fermat-mod mode, 4 block addresses in ascending order are add0,1,2,3 with no 'gaps' between blocks; for
	mersenne-mod addresses in asc. order are add0,2,3,1 with gaps between contiguous-data-block pairs 0,2 & 3,1.
	Thus for ferm-mod need [add2] as base-address of 'high-half' block for temp-storage; for mers-mod need [add3].
	In both cases have (add2 < add3) so instead use (add2 - add1) to differentiate: > 0 for fermat, < 0 for mersenne.
	*/\
		"movq	%[__add2],%%rsi		\n\t"/* destroyable copy of add2 */\
		"movq	%[__add2],%%rbx		\n\t"\
		"subq	%[__add1],%%rsi		\n\t"/* rsi = (add2 - add1); if this yields a borrow (i.e. sets CF) it's mersenne, else fermat. */\
	/*	"cmovcq %[__add3],%%rbx	\n\t" // CMOV not supported on k1om, so emulate using jump-if-CF-not-set: */\
	"jnc skip	\n\t"\
		"movq	%[__add3],%%rbx	\n\t" /* if CF set (i.e. h > l), copy source [add3] into dest (rbx), else leave dest = [add2]. */\
	"skip:	\n\t"\
		"movq	%[__isrt2],%%rsi				\n\t	vmovaps		 (%%rdi),%%zmm31	\n\t"/* 2.0*/\
													/* base-addr [rax] and sincos-ptr [r10] in rcol both offset +0x200 vs lcol */\
	/********************************************/	/********************************************/\
	/* Block 2: t02,t12,t22,t32 -> r10,14,12,16 */	/* Block 6: t0A,t1A,t2A,t3A -> r18,1C,1A,1E */\
	/********************************************/	/********************************************/\
		"movq	%[__r00],%%rax					\n\t"\
		"movq	%%rsi,%%rcx					\n\t"\
		"movq	%%rsi,%%rdx					\n\t"\
		"movq	%[__c01],%%r10					\n\t"\
		"addq	$0x040,%%rcx	/* cc0 */		\n\t"\
		"addq	$0x0c0,%%rdx	/* cc1 */		\n\t"\
		"vmovaps		0x480(%%rax),%%zmm4		\n\t		vmovaps		 0x680(%%rax),%%zmm12		\n\t"\
		"vmovaps		0x580(%%rax),%%zmm0		\n\t		vmovaps		 0x780(%%rax),%%zmm8 		\n\t"\
		"vmovaps		0x4c0(%%rax),%%zmm5		\n\t		vmovaps		 0x6c0(%%rax),%%zmm15		\n\t"\
		"vmovaps		0x5c0(%%rax),%%zmm3		\n\t		vmovaps		 0x7c0(%%rax),%%zmm11		\n\t"\
	/* Strategic register-usage for sincos consts here cuts loads from 16 to 4: */\
		"vmovaps			 (%%rdx),%%zmm9 	\n\t		vmovaps		0x040(%%rdx),%%zmm10		\n\t"\
		"vmovaps		0x080(%%rdx),%%zmm1 	\n\t		vmovaps		0x0c0(%%rdx),%%zmm13		\n\t"\
		"vmulpd			%%zmm10,%%zmm5,%%zmm7	\n\t		vmulpd			 %%zmm1 ,%%zmm15,%%zmm15\n\t"\
		"vmulpd			%%zmm13,%%zmm3,%%zmm3	\n\t		vmulpd			 %%zmm10,%%zmm11,%%zmm11\n\t"\
		"vmulpd			%%zmm10,%%zmm4,%%zmm6	\n\t		vmulpd			 %%zmm1 ,%%zmm12,%%zmm14\n\t"\
		"vmulpd			%%zmm13,%%zmm0,%%zmm2	\n\t		vmulpd		0x780(%%rax),%%zmm10,%%zmm10\n\t"\
	"vfmadd132pd		%%zmm9 ,%%zmm7,%%zmm4	\n\t	vfmadd132pd	 		 %%zmm13,%%zmm15,%%zmm12\n\t"\
	"vfmadd132pd		%%zmm1 ,%%zmm3,%%zmm0	\n\t	vfmsub132pd	 		 %%zmm9 ,%%zmm11,%%zmm8 \n\t"\
	"vfmsub132pd		%%zmm9 ,%%zmm6,%%zmm5	\n\t	vfmsub132pd	 	0x6c0(%%rax),%%zmm14,%%zmm13\n\t"\
	"vfmsub132pd   0x5c0(%%rax),%%zmm2,%%zmm1	\n\t	vfmadd132pd	 	0x7c0(%%rax),%%zmm10,%%zmm9 \n\t"\
		"vaddpd			%%zmm0,%%zmm4,%%zmm6	\n\t		vsubpd		%%zmm8 ,%%zmm12,%%zmm14		\n\t"\
		"vaddpd			%%zmm1,%%zmm5,%%zmm7	\n\t		vsubpd		%%zmm9 ,%%zmm13,%%zmm15		\n\t"\
		"vsubpd			%%zmm0,%%zmm4,%%zmm4	\n\t		vaddpd		%%zmm8 ,%%zmm12,%%zmm12		\n\t"\
		"vsubpd			%%zmm1,%%zmm5,%%zmm5	\n\t		vaddpd		%%zmm9 ,%%zmm13,%%zmm13		\n\t"\
		"vmovaps		 0x500(%%rax),%%zmm2	\n\t		vmovaps		0x700(%%rax),%%zmm10		\n\t"\
		"vmovaps		 0x540(%%rax),%%zmm3	\n\t		vmovaps		0x740(%%rax),%%zmm11		\n\t"\
		"vmovaps		  (%%rcx),%%zmm16		\n\t		vmovaps	 0x040(%%rcx),%%zmm17	\n\t"/* twiddle */\
		"vmulpd			%%zmm17,%%zmm2,%%zmm0	\n\t		vmulpd		%%zmm16,%%zmm10,%%zmm8	 	\n\t"\
		"vmulpd			%%zmm17,%%zmm3,%%zmm1	\n\t		vmulpd		%%zmm16,%%zmm11,%%zmm9 		\n\t"\
	"vfmsub132pd		%%zmm16,%%zmm0,%%zmm3	\n\t	vfmadd132pd		%%zmm17,%%zmm8 ,%%zmm11		\n\t"\
	"vfmadd132pd		%%zmm16,%%zmm1,%%zmm2	\n\t	vfmsub132pd		%%zmm17,%%zmm9 ,%%zmm10		\n\t"\
		"vmovaps		 0x400(%%rax),%%zmm0	\n\t		vmovaps		0x600(%%rax),%%zmm8 		\n\t"\
		"vmovaps		 0x440(%%rax),%%zmm1	\n\t		vmovaps		0x640(%%rax),%%zmm9 		\n\t"\
		"vsubpd			%%zmm2,%%zmm0,%%zmm0	\n\t		vsubpd		%%zmm10,%%zmm8 ,%%zmm8 		\n\t"\
		"vsubpd			%%zmm3,%%zmm1,%%zmm1	\n\t		vsubpd		%%zmm11,%%zmm9 ,%%zmm9 		\n\t"\
	"vfmadd132pd		%%zmm31,%%zmm0,%%zmm2	\n\t	vfmadd132pd		%%zmm31,%%zmm8 ,%%zmm10		\n\t"\
	"vfmadd132pd		%%zmm31,%%zmm1,%%zmm3	\n\t	vfmadd132pd		%%zmm31,%%zmm9 ,%%zmm11		\n\t"\
		"vsubpd			%%zmm6,%%zmm2,%%zmm2	\n\t		vsubpd		%%zmm14,%%zmm8 ,%%zmm8 		\n\t"\
		"vsubpd			%%zmm7,%%zmm3,%%zmm3	\n\t		vsubpd		%%zmm15,%%zmm9 ,%%zmm9 		\n\t"\
	"vfmadd132pd		%%zmm31,%%zmm2,%%zmm6	\n\t	vfmadd132pd		%%zmm31,%%zmm8 ,%%zmm14		\n\t"\
	"vfmadd132pd		%%zmm31,%%zmm3,%%zmm7	\n\t	vfmadd132pd		%%zmm31,%%zmm9 ,%%zmm15		\n\t"\
/*** Use regs m16-23 here to avoid 16-reg spill/fill crud: ***/\
		"vmovaps		  (%%r10),%%zmm20		\n\t		vmovaps		0x200(%%r10),%%zmm22		\n\t"\
		"vmovaps	 0x040(%%r10),%%zmm21		\n\t		vmovaps		0x240(%%r10),%%zmm23		\n\t"\
		"vmulpd			%%zmm20,%%zmm6,%%zmm16	\n\t		vmulpd		%%zmm22,%%zmm14,%%zmm18		\n\t"\
		"vmulpd			%%zmm20,%%zmm7,%%zmm17	\n\t		vmulpd		%%zmm22,%%zmm15,%%zmm19		\n\t"\
	" vfmadd231pd		%%zmm21,%%zmm7,%%zmm16	\n\t	 vfmadd231pd	%%zmm23,%%zmm15,%%zmm18		\n\t"\
	"vfnmadd231pd		%%zmm21,%%zmm6,%%zmm17	\n\t	vfnmadd231pd	%%zmm23,%%zmm14,%%zmm19		\n\t"\
		"vmovaps	 0x080(%%r10),%%zmm20		\n\t		vmovaps		0x280(%%r10),%%zmm22		\n\t"\
		"vmovaps	 0x0c0(%%r10),%%zmm21		\n\t		vmovaps		0x2c0(%%r10),%%zmm23		\n\t"\
		"vmovaps		%%zmm16, 0x400(%%rax)	\n\t		vmovaps		%%zmm18, 0x600(%%rax)		\n\t"\
		"vmovaps		%%zmm17, 0x440(%%rax)	\n\t		vmovaps		%%zmm19, 0x640(%%rax)		\n\t"\
		"vmulpd			%%zmm20,%%zmm2,%%zmm16	\n\t		vmulpd		%%zmm22,%%zmm8 ,%%zmm18		\n\t"\
		"vmulpd			%%zmm20,%%zmm3,%%zmm17	\n\t		vmulpd		%%zmm22,%%zmm9 ,%%zmm19		\n\t"\
		"vsubpd			%%zmm5,%%zmm0,%%zmm0	\n\t		vsubpd		%%zmm13,%%zmm10,%%zmm10		\n\t"\
		"vsubpd			%%zmm4,%%zmm1,%%zmm1	\n\t		vsubpd		%%zmm12,%%zmm11,%%zmm11		\n\t"\
	" vfmadd231pd		%%zmm21,%%zmm3,%%zmm16	\n\t	 vfmadd231pd	%%zmm23,%%zmm9 ,%%zmm18		\n\t"\
	"vfnmadd231pd		%%zmm21,%%zmm2,%%zmm17	\n\t	vfnmadd231pd	%%zmm23,%%zmm8 ,%%zmm19		\n\t"\
	"vfmadd132pd		%%zmm31,%%zmm0,%%zmm5	\n\t	vfmadd132pd		%%zmm31,%%zmm10,%%zmm13		\n\t"\
	"vfmadd132pd		%%zmm31,%%zmm1,%%zmm4	\n\t	vfmadd132pd		%%zmm31,%%zmm11,%%zmm12		\n\t"\
		"vmovaps		%%zmm16,0x480(%%rax)	\n\t		vmovaps		%%zmm18, 0x680(%%rax)		\n\t"\
		"vmovaps		%%zmm17,0x4c0(%%rax)	\n\t		vmovaps		%%zmm19, 0x6c0(%%rax)		\n\t"\
		"addq		$0x100,%%r10				\n\t"\
		"vmovaps		  (%%r10),%%zmm16		\n\t		vmovaps		0x200(%%r10),%%zmm20		\n\t"\
		"vmovaps	 0x080(%%r10),%%zmm17		\n\t		vmovaps		0x280(%%r10),%%zmm21		\n\t"\
		"vmovaps	 0x040(%%r10),%%zmm18		\n\t		vmovaps		0x240(%%r10),%%zmm22		\n\t"\
		"vmovaps	 0x0c0(%%r10),%%zmm19		\n\t		vmovaps		0x2c0(%%r10),%%zmm23		\n\t"\
		"vmulpd			%%zmm16,%%zmm5,%%zmm2	\n\t		vmulpd		%%zmm20,%%zmm13,%%zmm8 		\n\t"\
		"vmulpd			%%zmm17,%%zmm0,%%zmm6	\n\t		vmulpd		%%zmm21,%%zmm10,%%zmm14		\n\t"\
		"vmulpd			%%zmm16,%%zmm1,%%zmm3	\n\t		vmulpd		%%zmm20,%%zmm11,%%zmm9 		\n\t"\
		"vmulpd			%%zmm17,%%zmm4,%%zmm7	\n\t		vmulpd		%%zmm21,%%zmm12,%%zmm15		\n\t"\
	" vfmadd231pd		%%zmm18,%%zmm1,%%zmm2	\n\t	 vfmadd231pd	%%zmm22,%%zmm11,%%zmm8 		\n\t"\
	" vfmadd231pd		%%zmm19,%%zmm4,%%zmm6	\n\t	 vfmadd231pd	%%zmm23,%%zmm12,%%zmm14		\n\t"\
	"vfnmadd231pd		%%zmm18,%%zmm5,%%zmm3	\n\t	vfnmadd231pd	%%zmm22,%%zmm13,%%zmm9 		\n\t"\
	"vfnmadd231pd		%%zmm19,%%zmm0,%%zmm7	\n\t	vfnmadd231pd	%%zmm23,%%zmm10,%%zmm15		\n\t"\
		"vmovaps		%%zmm2, 0x500(%%rax)	\n\t		vmovaps		%%zmm8 , 0x700(%%rax)		\n\t"\
		"vmovaps		%%zmm3, 0x540(%%rax)	\n\t		vmovaps		%%zmm9 , 0x740(%%rax)		\n\t"\
		"vmovaps		%%zmm6, 0x580(%%rax)	\n\t		vmovaps		%%zmm14, 0x780(%%rax)		\n\t"\
		"vmovaps		%%zmm7, 0x5c0(%%rax)	\n\t		vmovaps		%%zmm15, 0x7c0(%%rax)		\n\t"\
		"\n\t"\
	/********************************************/	/********************************************/\
	/* Block 4: t06,t16,t26,t36 -> r30,34,32,36 */	/* Block 8: t0E,t1E,t2E,t3E -> r38,3C,3A,3E */\
	/********************************************/	/********************************************/\
		"													addq		$0x200,%%r10				\n\t"\
		"vmovaps		0xc80(%%rax),%%zmm4		\n\t		vmovaps		0xe80(%%rax),%%zmm12		\n\t"\
		"vmovaps		0xd80(%%rax),%%zmm0		\n\t		vmovaps		0xf80(%%rax),%%zmm10		\n\t"\
		"vmovaps		0xcc0(%%rax),%%zmm7		\n\t		vmovaps		0xec0(%%rax),%%zmm15		\n\t"\
		"vmovaps		0xdc0(%%rax),%%zmm1		\n\t		vmovaps		0xfc0(%%rax),%%zmm9 		\n\t"\
	/* Strategic register-usage for sincos consts here cuts loads from 16 to 4: */\
		"vmovaps			 (%%rdx),%%zmm14	\n\t		vmovaps		0x040(%%rdx),%%zmm13		\n\t"\
		"vmovaps		0x080(%%rdx),%%zmm5 	\n\t		vmovaps		0x0c0(%%rdx),%%zmm8 		\n\t"\
	/* Cyclic-permute 4 paired-MUL lines so as to put a '...(%%rax)' in the SRC3 slot of last line'srcol: */\
		"vmulpd		 	 %%zmm14,%%zmm0,%%zmm2	\n\t		vmulpd			 %%zmm5 ,%%zmm10,%%zmm10\n\t"\
		"vmulpd		 	 %%zmm8 ,%%zmm7,%%zmm7	\n\t		vmulpd			 %%zmm14,%%zmm15,%%zmm15\n\t"\
		"vmulpd		 	 %%zmm14,%%zmm1,%%zmm3	\n\t		vmulpd			 %%zmm5 ,%%zmm9 ,%%zmm11\n\t"\
		"vmulpd		 	 %%zmm8 ,%%zmm4,%%zmm6	\n\t		vmulpd		0xe80(%%rax),%%zmm14,%%zmm14\n\t"\
	/* Also Cyclic-permute these 4 lines so middle (SRC2) op's reg-index order matches that of abovefour: */\
	"vfmadd132pd	 	 %%zmm13,%%zmm2,%%zmm1	\n\t	vfmsub132pd	 		 %%zmm8 ,%%zmm10,%%zmm9 \n\t"\
	"vfmadd132pd	 	 %%zmm5 ,%%zmm7,%%zmm4	\n\t	vfmadd132pd	 		 %%zmm13,%%zmm15,%%zmm12\n\t"\
	"vfmsub132pd	 	 %%zmm13,%%zmm3,%%zmm0	\n\t	vfmadd132pd	 	0xf80(%%rax),%%zmm11,%%zmm8 \n\t"\
	"vfmsub132pd	0xcc0(%%rax),%%zmm6,%%zmm5	\n\t	vfmsub132pd	 	0xec0(%%rax),%%zmm14,%%zmm13\n\t"\
		"vsubpd			%%zmm0,%%zmm4,%%zmm6	\n\t		vsubpd		%%zmm8 ,%%zmm12,%%zmm14		\n\t"\
		"vsubpd			%%zmm1,%%zmm5,%%zmm7	\n\t		vsubpd		%%zmm9 ,%%zmm13,%%zmm15		\n\t"\
		"vaddpd			%%zmm0,%%zmm4,%%zmm4	\n\t		vaddpd		%%zmm8 ,%%zmm12,%%zmm12		\n\t"\
		"vaddpd			%%zmm1,%%zmm5,%%zmm5	\n\t		vaddpd		%%zmm9 ,%%zmm13,%%zmm13		\n\t"\
		"vmovaps		 0xd00(%%rax),%%zmm2	\n\t		vmovaps		0xf00(%%rax),%%zmm10		\n\t"\
		"vmovaps		 0xd40(%%rax),%%zmm3	\n\t		vmovaps		0xf40(%%rax),%%zmm11		\n\t"\
		"vmovaps		  (%%rcx),%%zmm16		\n\t		vmovaps	 0x040(%%rcx),%%zmm17	\n\t"/* twiddle */\
		"vmulpd			%%zmm16,%%zmm2,%%zmm0	\n\t		vmulpd		%%zmm17,%%zmm10,%%zmm8	 	\n\t"\
		"vmulpd			%%zmm16,%%zmm3,%%zmm1	\n\t		vmulpd		%%zmm17,%%zmm11,%%zmm9 		\n\t"\
	"vfmsub132pd		%%zmm17,%%zmm0,%%zmm3	\n\t	vfmadd132pd		%%zmm16,%%zmm8 ,%%zmm11		\n\t"\
	"vfmadd132pd		%%zmm17,%%zmm1,%%zmm2	\n\t	vfmsub132pd		%%zmm16,%%zmm9 ,%%zmm10		\n\t"\
		"vmovaps		 0xc00(%%rax),%%zmm0	\n\t		vmovaps		0xe00(%%rax),%%zmm8 		\n\t"\
		"vmovaps		 0xc40(%%rax),%%zmm1	\n\t		vmovaps		0xe40(%%rax),%%zmm9 		\n\t"\
		"vsubpd			%%zmm2,%%zmm0,%%zmm0	\n\t		vsubpd		%%zmm10,%%zmm8 ,%%zmm8 		\n\t"\
		"vsubpd			%%zmm3,%%zmm1,%%zmm1	\n\t		vsubpd		%%zmm11,%%zmm9 ,%%zmm9 		\n\t"\
	"vfmadd132pd		%%zmm31,%%zmm0,%%zmm2	\n\t	vfmadd132pd		%%zmm31,%%zmm8 ,%%zmm10		\n\t"\
	"vfmadd132pd		%%zmm31,%%zmm1,%%zmm3	\n\t	vfmadd132pd		%%zmm31,%%zmm9 ,%%zmm11		\n\t"\
		"vsubpd			%%zmm6,%%zmm2,%%zmm2	\n\t		vsubpd		%%zmm14,%%zmm8 ,%%zmm8 		\n\t"\
		"vsubpd			%%zmm7,%%zmm3,%%zmm3	\n\t		vsubpd		%%zmm15,%%zmm9 ,%%zmm9 		\n\t"\
	"vfmadd132pd		%%zmm31,%%zmm2,%%zmm6	\n\t	vfmadd132pd		%%zmm31,%%zmm8 ,%%zmm14		\n\t"\
	"vfmadd132pd		%%zmm31,%%zmm3,%%zmm7	\n\t	vfmadd132pd		%%zmm31,%%zmm9 ,%%zmm15		\n\t"\
		"addq		$0x100,%%r10				\n\t"\
		"vmovaps		  (%%r10),%%zmm20		\n\t		vmovaps		0x200(%%r10),%%zmm22		\n\t"\
		"vmovaps	 0x040(%%r10),%%zmm21		\n\t		vmovaps		0x240(%%r10),%%zmm23		\n\t"\
		"vmulpd			%%zmm20,%%zmm6,%%zmm16	\n\t		vmulpd		%%zmm22,%%zmm14,%%zmm18		\n\t"\
		"vmulpd			%%zmm20,%%zmm7,%%zmm17	\n\t		vmulpd		%%zmm22,%%zmm15,%%zmm19		\n\t"\
	" vfmadd231pd		%%zmm21,%%zmm7,%%zmm16	\n\t	 vfmadd231pd	%%zmm23,%%zmm15,%%zmm18		\n\t"\
	"vfnmadd231pd		%%zmm21,%%zmm6,%%zmm17	\n\t	vfnmadd231pd	%%zmm23,%%zmm14,%%zmm19		\n\t"\
		"vmovaps	 0x080(%%r10),%%zmm20		\n\t		vmovaps		0x280(%%r10),%%zmm22		\n\t"\
		"vmovaps	 0x0c0(%%r10),%%zmm21		\n\t		vmovaps		0x2c0(%%r10),%%zmm23		\n\t"\
		"vmulpd			%%zmm20,%%zmm2,%%zmm6	\n\t		vmulpd		%%zmm22,%%zmm8 ,%%zmm14		\n\t"\
		"vmulpd			%%zmm20,%%zmm3,%%zmm7	\n\t		vmulpd		%%zmm22,%%zmm9 ,%%zmm15		\n\t"\
		"vsubpd			%%zmm5,%%zmm0,%%zmm0	\n\t		vsubpd		%%zmm13,%%zmm10,%%zmm10		\n\t"\
		"vsubpd			%%zmm4,%%zmm1,%%zmm1	\n\t		vsubpd		%%zmm12,%%zmm11,%%zmm11		\n\t"\
	" vfmadd231pd		%%zmm21,%%zmm3,%%zmm6	\n\t	 vfmadd231pd	%%zmm23,%%zmm9 ,%%zmm14		\n\t"\
	"vfnmadd231pd		%%zmm21,%%zmm2,%%zmm7	\n\t	vfnmadd231pd	%%zmm23,%%zmm8 ,%%zmm15		\n\t"\
	"vfmadd132pd		%%zmm31,%%zmm0,%%zmm5	\n\t	vfmadd132pd		%%zmm31,%%zmm10,%%zmm13		\n\t"\
	"vfmadd132pd		%%zmm31,%%zmm1,%%zmm4	\n\t	vfmadd132pd		%%zmm31,%%zmm11,%%zmm12		\n\t"\
		"vmovaps		%%zmm16, 0xc00(%%rax)	\n\t		vmovaps		%%zmm18, 0xe00(%%rax)		\n\t"\
		"vmovaps		%%zmm17, 0xc40(%%rax)	\n\t		vmovaps		%%zmm19, 0xe40(%%rax)		\n\t"\
		"vmovaps		%%zmm6 , 0xc80(%%rax)	\n\t		vmovaps		%%zmm14, 0xe80(%%rax)		\n\t"\
		"vmovaps		%%zmm7 , 0xcc0(%%rax)	\n\t		vmovaps		%%zmm15, 0xec0(%%rax)		\n\t"\
		"addq		$0x100,%%r10				\n\t"\
		"vmovaps		  (%%r10),%%zmm16		\n\t		vmovaps		0x200(%%r10),%%zmm20		\n\t"\
		"vmovaps	 0x080(%%r10),%%zmm17		\n\t		vmovaps		0x280(%%r10),%%zmm21		\n\t"\
		"vmovaps	 0x040(%%r10),%%zmm18		\n\t		vmovaps		0x240(%%r10),%%zmm22		\n\t"\
		"vmovaps	 0x0c0(%%r10),%%zmm19		\n\t		vmovaps		0x2c0(%%r10),%%zmm23		\n\t"\
		"vmulpd			%%zmm16,%%zmm5,%%zmm2	\n\t		vmulpd		%%zmm20,%%zmm13,%%zmm8 		\n\t"\
		"vmulpd			%%zmm16,%%zmm1,%%zmm3	\n\t		vmulpd		%%zmm20,%%zmm11,%%zmm9 		\n\t"\
		"vmulpd			%%zmm17,%%zmm0,%%zmm6	\n\t		vmulpd		%%zmm21,%%zmm10,%%zmm14		\n\t"\
		"vmulpd			%%zmm17,%%zmm4,%%zmm7	\n\t		vmulpd		%%zmm21,%%zmm12,%%zmm15		\n\t"\
	" vfmadd231pd		%%zmm18,%%zmm1,%%zmm2	\n\t	 vfmadd231pd	%%zmm22,%%zmm11,%%zmm8 		\n\t"\
	"vfnmadd231pd		%%zmm18,%%zmm5,%%zmm3	\n\t	vfnmadd231pd	%%zmm22,%%zmm13,%%zmm9 		\n\t"\
	" vfmadd231pd		%%zmm19,%%zmm4,%%zmm6	\n\t	 vfmadd231pd	%%zmm23,%%zmm12,%%zmm14		\n\t"\
	"vfnmadd231pd		%%zmm19,%%zmm0,%%zmm7	\n\t	vfnmadd231pd	%%zmm23,%%zmm10,%%zmm15		\n\t"\
		"vmovaps		%%zmm2, 0xd00(%%rax)	\n\t		vmovaps		%%zmm8 , 0xf00(%%rax)		\n\t"\
		"vmovaps		%%zmm3, 0xd40(%%rax)	\n\t		vmovaps		%%zmm9 , 0xf40(%%rax)		\n\t"\
		"vmovaps		%%zmm6, 0xd80(%%rax)	\n\t		vmovaps		%%zmm14, 0xf80(%%rax)		\n\t"\
		"vmovaps		%%zmm7, 0xdc0(%%rax)	\n\t		vmovaps		%%zmm15, 0xfc0(%%rax)		\n\t"\
		"\n\t"\
	/********************************************/	/********************************************/\
	/* Block 1: t00,t10,t20,t30 -> r00,04,02,06 */	/* Block 5: t08,t18,t28,t38 -> r08,0C,0A,0E */\
	/********************************************/	/********************************************/\
		"													vmovaps	(%%rsi),%%zmm10	\n\t"/* isrt2 */\
		"vmovaps			  (%%rax),%%zmm0	\n\t		vmovaps		 0x280(%%rax),%%zmm12		\n\t"\
		"vmovaps		 0x040(%%rax),%%zmm1	\n\t		vmovaps		 0x2c0(%%rax),%%zmm13		\n\t"\
		"vmovaps		 0x100(%%rax),%%zmm2	\n\t		vmovaps		 0x380(%%rax),%%zmm8 		\n\t"\
		"vmovaps		 0x140(%%rax),%%zmm3	\n\t		vmovaps		 0x3c0(%%rax),%%zmm9 		\n\t"\
		"vsubpd			%%zmm2 ,%%zmm0,%%zmm0	\n\t		vsubpd		%%zmm12,%%zmm13,%%zmm13		\n\t"\
		"vsubpd			%%zmm3 ,%%zmm1,%%zmm1	\n\t	vfmadd132pd		%%zmm31,%%zmm13,%%zmm12		\n\t"\
	"vfmadd132pd		%%zmm31,%%zmm0,%%zmm2	\n\t		vsubpd		%%zmm9 ,%%zmm8 ,%%zmm8 		\n\t"\
	"vfmadd132pd		%%zmm31,%%zmm1,%%zmm3	\n\t	vfmadd132pd		%%zmm31,%%zmm8 ,%%zmm9 		\n\t"\
		"vmovaps		 0x080(%%rax),%%zmm4	\n\t		vmulpd		%%zmm10,%%zmm8 ,%%zmm8 		\n\t"\
		"vmovaps		 0x0c0(%%rax),%%zmm5	\n\t		vmulpd		%%zmm10,%%zmm9 ,%%zmm9 		\n\t"\
		"vmovaps		 0x180(%%rax),%%zmm6	\n\t		vmovaps		%%zmm12,%%zmm14				\n\t"\
		"vmovaps		 0x1c0(%%rax),%%zmm7	\n\t		vmovaps		%%zmm13,%%zmm15				\n\t"\
		"vsubpd			%%zmm6 ,%%zmm4,%%zmm4	\n\t	vfmsub132pd		%%zmm10,%%zmm8 ,%%zmm12		\n\t"\
		"vsubpd			%%zmm7 ,%%zmm5,%%zmm5	\n\t	vfmsub132pd		%%zmm10,%%zmm9 ,%%zmm13		\n\t"\
	"vfmadd132pd		%%zmm31,%%zmm4,%%zmm6	\n\t	vfmadd132pd		%%zmm10,%%zmm8 ,%%zmm14		\n\t"\
	"vfmadd132pd		%%zmm31,%%zmm5,%%zmm7	\n\t	vfmadd132pd		%%zmm10,%%zmm9 ,%%zmm15		\n\t"\
		"movq		%[__c10],%%rcx			\n\t"/* base-twiddle in l/rcol = c00/c04 ==> rcx+0x200 in rcol */\
		"vaddpd			%%zmm6,%%zmm2,%%zmm2	\n\t		vmovaps		 0x200(%%rax),%%zmm8 		\n\t"\
		"vaddpd			%%zmm7,%%zmm3,%%zmm3	\n\t		vmovaps		 0x240(%%rax),%%zmm9 		\n\t"\
		"vmovaps		%%zmm2,      (%%rax)	\n\t		vmovaps		 0x300(%%rax),%%zmm10		\n\t"\
		"vmovaps		%%zmm3, 0x040(%%rax)	\n\t		vmovaps		 0x340(%%rax),%%zmm11		\n\t"\
		"													vsubpd		%%zmm11,%%zmm8 ,%%zmm8 		\n\t"\
		"													vsubpd		%%zmm10,%%zmm9 ,%%zmm9 		\n\t"\
	"vfnmadd231pd		%%zmm31,%%zmm6,%%zmm2	\n\t	vfmadd132pd		%%zmm31,%%zmm8 ,%%zmm11		\n\t"\
	"vfnmadd231pd		%%zmm31,%%zmm7,%%zmm3	\n\t	vfmadd132pd		%%zmm31,%%zmm9 ,%%zmm10		\n\t"\
		"vmovaps		%%zmm2,%%zmm6			\n\t		vsubpd		%%zmm12,%%zmm11,%%zmm11		\n\t"\
		"vmovaps		%%zmm3,%%zmm7			\n\t		vsubpd		%%zmm13,%%zmm9 ,%%zmm9 		\n\t"\
		"vmovaps		  (%%rcx),%%zmm16		\n\t		vmovaps		0x180(%%rcx),%%zmm20		\n\t"\
		"vmovaps	 0x040(%%rcx),%%zmm17		\n\t		vmovaps		0x1c0(%%rcx),%%zmm21		\n\t"\
		"													vmovaps		0x200(%%rcx),%%zmm22		\n\t"\
		"													vmovaps		0x240(%%rcx),%%zmm23		\n\t"\
		"subq $0x80,%%rcx		\n\t"/* put c00 in rcx to ease bookkeeping */\
/*c10:*/"vmulpd			%%zmm16,%%zmm2,%%zmm2	\n\t	vfmadd132pd		%%zmm31,%%zmm11,%%zmm12		\n\t"\
		"vmulpd			%%zmm16,%%zmm3,%%zmm3	\n\t	vfmadd132pd		%%zmm31,%%zmm9 ,%%zmm13		\n\t"\
	" vfmadd231pd		%%zmm17,%%zmm7,%%zmm2	\n\t		vmulpd		%%zmm20,%%zmm12,%%zmm16		\n\t"/* c04:*/\
	"vfnmadd231pd		%%zmm17,%%zmm6,%%zmm3	\n\t		vmulpd		%%zmm20,%%zmm13,%%zmm17		\n\t"\
		"												 vfmadd231pd	%%zmm21,%%zmm13,%%zmm16		\n\t"\
		"												vfnmadd231pd	%%zmm21,%%zmm12,%%zmm17		\n\t"\
		"													vmulpd	 %%zmm22,%%zmm11,%%zmm12		\n\t"/* c14 */\
		"													vmulpd	 %%zmm22,%%zmm9 ,%%zmm13		\n\t"\
		"												 vfmadd231pd %%zmm23,%%zmm9 ,%%zmm12		\n\t"\
		"												vfnmadd231pd %%zmm23,%%zmm11,%%zmm13		\n\t"\
		"vaddpd			%%zmm5,%%zmm0,%%zmm0	\n\t		vmovaps	%%zmm12,0x280(%%rax)			\n\t"/* r0a,0b */\
		"vsubpd			%%zmm4,%%zmm1,%%zmm1	\n\t		vmovaps	%%zmm13,0x2c0(%%rax)			\n\t"\
		"vmovaps	%%zmm2,0x080(%%rax)	/* r02,03 */\n\t	vmovaps	%%zmm16,0x200(%%rax)			\n\t"/* r08,09 */\
		"vmovaps	%%zmm3,0x0c0(%%rax)				\n\t	vmovaps	%%zmm17,0x240(%%rax)			\n\t"\
		"vmovaps	 0x100(%%rcx),%%zmm16		\n\t		vmovaps		0x300(%%rcx),%%zmm20		\n\t"\
		"vmovaps	 0x140(%%rcx),%%zmm17		\n\t		vmovaps		0x340(%%rcx),%%zmm21		\n\t"\
		"vmovaps	 0x180(%%rcx),%%zmm18		\n\t		vmovaps		0x380(%%rcx),%%zmm22		\n\t"\
		"vmovaps	 0x1c0(%%rcx),%%zmm19		\n\t		vmovaps		0x3c0(%%rcx),%%zmm23		\n\t"\
	"vfnmadd132pd	%%zmm31,%%zmm0,%%zmm5		\n\t		vsubpd			%%zmm15,%%zmm8 ,%%zmm8 	\n\t"\
	" vfmadd132pd	%%zmm31,%%zmm1,%%zmm4		\n\t		vsubpd			%%zmm14,%%zmm10,%%zmm10	\n\t"\
		"vmulpd		%%zmm16,%%zmm0,%%zmm2 /*c08*/\n\t	vfmadd132pd		%%zmm31,%%zmm8 ,%%zmm15		\n\t"\
		"vmulpd		%%zmm16,%%zmm1,%%zmm3		\n\t	vfmadd132pd		%%zmm31,%%zmm10,%%zmm14		\n\t"\
	" vfmadd231pd	%%zmm17,%%zmm1,%%zmm2		\n\t		vmulpd		%%zmm20,%%zmm15,%%zmm12		\n\t"/* c0C */\
	"vfnmadd231pd	%%zmm17,%%zmm0,%%zmm3		\n\t		vmulpd		%%zmm20,%%zmm10,%%zmm13		\n\t"\
		"												 vfmadd231pd	%%zmm21,%%zmm10,%%zmm12		\n\t"\
		"												vfnmadd231pd	%%zmm21,%%zmm15,%%zmm13		\n\t"\
		"vmulpd		%%zmm18,%%zmm5,%%zmm0	/* c18*/\n\t	vmulpd		%%zmm22,%%zmm8 ,%%zmm15		\n\t"/* c1C */\
		"vmulpd		%%zmm18,%%zmm4,%%zmm1		\n\t		vmulpd		%%zmm22,%%zmm14,%%zmm10		\n\t"\
	" vfmadd231pd	%%zmm19,%%zmm4,%%zmm0		\n\t	 vfmadd231pd	%%zmm23,%%zmm14,%%zmm15		\n\t"\
	"vfnmadd231pd	%%zmm19,%%zmm5,%%zmm1		\n\t	vfnmadd231pd	%%zmm23,%%zmm8 ,%%zmm10		\n\t"\
		"vmovaps	%%zmm2,0x100(%%rax)	/* r04,05 */\n\t	vmovaps	%%zmm12,0x300(%%rax)			\n\t"/* r0c,0d */\
		"vmovaps	%%zmm3,0x140(%%rax)				\n\t	vmovaps	%%zmm13,0x340(%%rax)			\n\t"\
		"vmovaps	%%zmm0,0x180(%%rax)	/* r06,07 */\n\t	vmovaps	%%zmm15,0x380(%%rax)			\n\t"/* r0e,0f */\
		"vmovaps	%%zmm1,0x1c0(%%rax)			\n\t		vmovaps	%%zmm10,0x3c0(%%rax)			\n\t"\
		"\n\t"\
	/********************************************/	/********************************************/\
	/* Block 3: t04,t14,t24,t34 -> r20,24,22,26 */	/* Block 7: t0C,t1C,t2C,t3C -> r28,2C,2A,2E */\
	/********************************************/	/********************************************/\
		"leaq	0x040(%%rsi),%%rcx	\n\t"/* cc0; Note cc0/ss0 are shared between lcol/rcol, so no rcx-offset until get to twiddles*/\
		"vmovaps		0x880(%%rax),%%zmm4		\n\t		vmovaps		 0xa80(%%rax),%%zmm12		\n\t"\
		"vmovaps		0x980(%%rax),%%zmm0		\n\t		vmovaps		 0xb80(%%rax),%%zmm8 		\n\t"\
		"vmovaps		0x8c0(%%rax),%%zmm5		\n\t		vmovaps		 0xac0(%%rax),%%zmm13		\n\t"\
		"vmovaps		0x9c0(%%rax),%%zmm3		\n\t		vmovaps		 0xbc0(%%rax),%%zmm11		\n\t"\
		"vmovaps			%%zmm4 ,%%zmm6		\n\t		vmovaps		 	%%zmm12,%%zmm14			\n\t"\
		"vmovaps			%%zmm0 ,%%zmm2		\n\t		vmovaps		 	%%zmm8 ,%%zmm10			\n\t"\
		"vmovaps			%%zmm5 ,%%zmm7		\n\t		vmovaps		 	%%zmm13,%%zmm15			\n\t"\
	/*	"vmovaps			%%zmm1 ,%%zmm3		\n\t		vmovaps		 	%%zmm9 ,%%zmm11			\n\t"*/\
	/* Strategic register-usage for sincos consts here cuts loads from 16 to 4: */\
		"vmovaps			 (%%rcx),%%zmm9 	\n\t		vmovaps		0x040(%%rcx),%%zmm1 		\n\t"\
		"vmulpd		 	 %%zmm1 ,%%zmm7,%%zmm7	\n\t		vmulpd			 %%zmm9 ,%%zmm15,%%zmm15\n\t"\
		"vmulpd		 	 %%zmm9 ,%%zmm3,%%zmm3	\n\t		vmulpd			 %%zmm1 ,%%zmm11,%%zmm11\n\t"\
		"vmulpd		 	 %%zmm1 ,%%zmm6,%%zmm6	\n\t		vmulpd			 %%zmm9 ,%%zmm14,%%zmm14\n\t"\
		"vmulpd		 	 %%zmm9 ,%%zmm2,%%zmm2	\n\t		vmulpd			 %%zmm1 ,%%zmm10,%%zmm10\n\t"\
	"vfmadd132pd	 	 %%zmm9 ,%%zmm7,%%zmm4	\n\t	vfmadd132pd			 %%zmm1 ,%%zmm15,%%zmm12\n\t"\
	"vfmadd132pd	 	 %%zmm1 ,%%zmm3,%%zmm0	\n\t	vfmadd132pd			 %%zmm9 ,%%zmm11,%%zmm8 \n\t"\
	"vfmsub132pd	 	 %%zmm9 ,%%zmm6,%%zmm5	\n\t	vfmsub132pd			 %%zmm1 ,%%zmm14,%%zmm13\n\t"\
	"vfmsub132pd	0x9c0(%%rax),%%zmm2,%%zmm1	\n\t	vfmsub132pd		0xbc0(%%rax),%%zmm10,%%zmm9 \n\t"\
		"vsubpd			%%zmm0,%%zmm4,%%zmm6	\n\t		vsubpd			%%zmm8 ,%%zmm12,%%zmm14	\n\t"\
		"vsubpd			%%zmm1,%%zmm5,%%zmm7	\n\t		vsubpd			%%zmm9 ,%%zmm13,%%zmm15	\n\t"\
		"vaddpd			%%zmm0,%%zmm4,%%zmm4	\n\t		vaddpd			%%zmm8 ,%%zmm12,%%zmm12	\n\t"\
		"vaddpd			%%zmm1,%%zmm5,%%zmm5	\n\t		vaddpd			%%zmm9 ,%%zmm13,%%zmm13	\n\t"\
		"vmovaps		 0x900(%%rax),%%zmm2	\n\t		vmovaps		0xb00(%%rax),%%zmm10		\n\t"\
		"vmovaps		 0x940(%%rax),%%zmm3	\n\t		vmovaps		0xb40(%%rax),%%zmm11		\n\t"\
		"vmovaps		 0x800(%%rax),%%zmm0	\n\t		vmovaps		0xa00(%%rax),%%zmm8 		\n\t"\
		"vmovaps		 0x840(%%rax),%%zmm1	\n\t		vmovaps		0xa40(%%rax),%%zmm9 		\n\t"\
		"vaddpd		 0x940(%%rax),%%zmm2,%%zmm2	\n\t		vsubpd		0xb40(%%rax),%%zmm10,%%zmm10\n\t"\
		"vsubpd		 0x900(%%rax),%%zmm3,%%zmm3	\n\t		vaddpd		0xb00(%%rax),%%zmm11,%%zmm11\n\t"\
		"vmulpd			  (%%rsi),%%zmm2,%%zmm2	\n\t		vmulpd			 (%%rsi),%%zmm10,%%zmm10\n\t"\
		"vmulpd			  (%%rsi),%%zmm3,%%zmm3	\n\t		vmulpd			 (%%rsi),%%zmm11,%%zmm11\n\t"\
		"vsubpd			%%zmm2 ,%%zmm0,%%zmm0	\n\t		vsubpd			%%zmm10,%%zmm8 ,%%zmm8 	\n\t"\
		"vsubpd			%%zmm3 ,%%zmm1,%%zmm1	\n\t		vsubpd			%%zmm11,%%zmm9 ,%%zmm9 	\n\t"\
	"vfmadd132pd		%%zmm31,%%zmm0,%%zmm2	\n\t	vfmadd132pd		 %%zmm31,%%zmm8 ,%%zmm10	\n\t"\
	"vfmadd132pd		%%zmm31,%%zmm1,%%zmm3	\n\t	vfmadd132pd		 %%zmm31,%%zmm9 ,%%zmm11	\n\t"\
		"movq		%[__c02],%%rcx				\n\t"/* base-twiddle addr in rcol = c06, so rcx offset +0x200 vs lcol */\
		"vmovaps		  (%%rcx),%%zmm16		\n\t		vmovaps		0x200(%%rcx),%%zmm20		\n\t"\
		"vmovaps	 0x040(%%rcx),%%zmm17		\n\t		vmovaps		0x240(%%rcx),%%zmm21		\n\t"\
		"vsubpd			%%zmm4,%%zmm2,%%zmm2	\n\t		vsubpd			%%zmm14,%%zmm8 ,%%zmm8 	\n\t"\
		"vsubpd			%%zmm5,%%zmm3,%%zmm3	\n\t		vsubpd			%%zmm15,%%zmm9 ,%%zmm9 	\n\t"\
	"vfmadd132pd		%%zmm31,%%zmm2,%%zmm4	\n\t	vfmadd132pd		%%zmm31,%%zmm8 ,%%zmm14	\n\t"\
	"vfmadd132pd		%%zmm31,%%zmm3,%%zmm5	\n\t	vfmadd132pd		%%zmm31,%%zmm9 ,%%zmm15	\n\t"\
		"vmovaps		%%zmm2, 0x880(%%rax)	\n\t		vmovaps		%%zmm8 , 0xa80(%%rax)		\n\t"\
		"vmovaps		%%zmm3, 0x8c0(%%rax)	\n\t		vmovaps		%%zmm9 , 0xac0(%%rax)		\n\t"\
		"vmulpd			%%zmm16,%%zmm4,%%zmm2	\n\t		vmulpd			%%zmm20,%%zmm14,%%zmm8	\n\t"\
		"vmulpd			%%zmm16,%%zmm5,%%zmm3	\n\t		vmulpd			%%zmm20,%%zmm15,%%zmm9	\n\t"\
	" vfmadd231pd		%%zmm17,%%zmm5,%%zmm2	\n\t	 vfmadd231pd		%%zmm21,%%zmm15,%%zmm8	\n\t"\
	"vfnmadd231pd		%%zmm17,%%zmm4,%%zmm3	\n\t	vfnmadd231pd		%%zmm21,%%zmm14,%%zmm9	\n\t"\
		"movq		%[__c12],%%rcx				\n\t"		/* rcol uses c16 */\
		"vmovaps		  (%%rcx),%%zmm16		\n\t		vmovaps		0x200(%%rcx),%%zmm20		\n\t"\
		"vmovaps	 0x040(%%rcx),%%zmm17		\n\t		vmovaps		0x240(%%rcx),%%zmm21		\n\t"\
		"vmovaps		 0x880(%%rax),%%zmm4	\n\t		vmovaps		 0xa80(%%rax),%%zmm14		\n\t"\
		"vmovaps		 0x8c0(%%rax),%%zmm5	\n\t		vmovaps		 0xac0(%%rax),%%zmm15		\n\t"\
		"vmovaps	%%zmm2,0x800(%%rax)	/* r20,21 */\n\t	vmovaps	%%zmm8 ,0xa00(%%rax)			\n\t"/* r28,29 */\
		"vmovaps	%%zmm3,0x840(%%rax)			\n\t		vmovaps	%%zmm9 ,0xa40(%%rax)			\n\t"\
		"vmulpd			%%zmm16,%%zmm4,%%zmm2	\n\t		vmulpd			%%zmm20,%%zmm14,%%zmm8	\n\t"\
		"vmulpd			%%zmm16,%%zmm5,%%zmm3	\n\t		vmulpd			%%zmm20,%%zmm15,%%zmm9	\n\t"\
	" vfmadd231pd		%%zmm17,%%zmm5,%%zmm2	\n\t	 vfmadd231pd		%%zmm21,%%zmm15,%%zmm8	\n\t"\
	"vfnmadd231pd		%%zmm17,%%zmm4,%%zmm3	\n\t	vfnmadd231pd		%%zmm21,%%zmm14,%%zmm9	\n\t"\
		"movq		%[__c0A],%%rcx				\n\t"		/* rcol uses c0E */\
		"vmovaps		  (%%rcx),%%zmm16		\n\t		vmovaps		0x200(%%rcx),%%zmm20		\n\t"\
		"vmovaps	 0x040(%%rcx),%%zmm17		\n\t		vmovaps		0x240(%%rcx),%%zmm21		\n\t"\
		"vsubpd			%%zmm7,%%zmm0,%%zmm0	\n\t		vsubpd			%%zmm13,%%zmm10,%%zmm10	\n\t"\
		"vsubpd			%%zmm6,%%zmm1,%%zmm1	\n\t		vsubpd			%%zmm12,%%zmm11,%%zmm11	\n\t"\
	"vfmadd132pd		%%zmm31,%%zmm0,%%zmm7	\n\t	vfmadd132pd		%%zmm31,%%zmm10,%%zmm13		\n\t"\
	"vfmadd132pd		%%zmm31,%%zmm1,%%zmm6	\n\t	vfmadd132pd		%%zmm31,%%zmm11,%%zmm12		\n\t"\
		"vmovaps	%%zmm2,0x880(%%rax)	/* r22,23 */\n\t	vmovaps	%%zmm8 ,0xa80(%%rax)			\n\t"/* r2a,2b */\
		"vmovaps	%%zmm3,0x8c0(%%rax)			\n\t		vmovaps	%%zmm9 ,0xac0(%%rax)			\n\t"\
		"vmulpd			%%zmm16,%%zmm7,%%zmm4	\n\t		vmulpd			%%zmm20,%%zmm13,%%zmm8	\n\t"\
		"vmulpd			%%zmm16,%%zmm1,%%zmm5	\n\t		vmulpd			%%zmm20,%%zmm11,%%zmm9	\n\t"\
	" vfmadd231pd		%%zmm17,%%zmm1,%%zmm4	\n\t	 vfmadd231pd		%%zmm21,%%zmm11,%%zmm8	\n\t"\
	"vfnmadd231pd		%%zmm17,%%zmm7,%%zmm5	\n\t	vfnmadd231pd		%%zmm21,%%zmm13,%%zmm9	\n\t"\
		"movq		%[__c1A],%%rcx				\n\t"		/* rcol uses c1E */\
		"vmovaps		  (%%rcx),%%zmm16		\n\t		vmovaps		0x200(%%rcx),%%zmm20		\n\t"\
		"vmovaps	 0x040(%%rcx),%%zmm17		\n\t		vmovaps		0x240(%%rcx),%%zmm21		\n\t"\
		"vmovaps	%%zmm4,0x900(%%rax)	/* r24,25 */\n\t	vmovaps	%%zmm8 ,0xb00(%%rax)			\n\t"/* r2c,2d */\
		"vmovaps	%%zmm5,0x940(%%rax)			\n\t		vmovaps	%%zmm9 ,0xb40(%%rax)			\n\t"\
		"vmulpd			%%zmm16,%%zmm0,%%zmm4	\n\t		vmulpd			%%zmm20,%%zmm10,%%zmm8	\n\t"\
		"vmulpd			%%zmm16,%%zmm6,%%zmm5	\n\t		vmulpd			%%zmm20,%%zmm12,%%zmm9	\n\t"\
	" vfmadd231pd		%%zmm17,%%zmm6,%%zmm4	\n\t	 vfmadd231pd		%%zmm21,%%zmm12,%%zmm8	\n\t"\
	"vfnmadd231pd		%%zmm17,%%zmm0,%%zmm5	\n\t	vfnmadd231pd		%%zmm21,%%zmm10,%%zmm9	\n\t"\
		"vmovaps	%%zmm4,0x980(%%rax)	/* r26,27 */\n\t	vmovaps	%%zmm8 ,0xb80(%%rax)			\n\t"/* r2e,2f */\
		"vmovaps	%%zmm5,0x9c0(%%rax)			\n\t		vmovaps	%%zmm9 ,0xbc0(%%rax)			\n\t"\
	/*******************************************
	***** Finish with 8-way 'un'terleaving: ****
	Using the AVX-512 data layout, the rcol pattern is:
		a[ 0- 7] = re[ 0, 8,16,24, 4,12,20,28].d0	a[ 8-15] = im[ 0, 8,16,24, 4,12,20,28].d0
		a[16-23] = re[ 2,10,18,26, 6,14,22,30].d0	a[24-31] = im[ 2,10,18,26, 6,14,22,30].d0
		a[32-39] = re[ 1, 9,17,25, 5,13,21,29].d0	a[40-47] = im[ 1, 9,17,25, 5,13,21,29].d0
		a[48-55] = re[ 3,11,19,27, 7,15,23,31].d0	a[56-63] = im[ 3,11,19,27, 7,15,23,31].d0 ,
	and remaining seven 64-double blocks repeat same pattern with elts d1-d7 of the vector-doubles.
	*******************************************/\
		"movq	%[__add0],%%rax	\n\t"\
		"movq	%[__add1],%%rbx	\n\t"\
		"movq	%[__add2],%%rcx	\n\t"\
		"movq	%[__add3],%%rdx	\n\t"\
		"movq	%[__add4],%%r10	\n\t"\
		"movq	%[__add5],%%r11	\n\t"\
		"movq	%[__add6],%%r12	\n\t"\
		"movq	%[__add7],%%r13	\n\t"\
	"movq	%[__r00] ,%%rsi	\n\t"\
	/**** See above USE_16_REG version of this macro for closing transposition using columnwise scatter-stores. ****/\
	/**** a[ 0- 7] = re[ 0, 8,16,24, 4,12,20,28].d0	a[ 8-15] = im[ 0, 8,16,24, 4,12,20,28].d0 : ****/\
		/* Read transposed row-data from the local store: */\
		"vmovaps 		 (%%rsi),%%zmm0			\n\t	vmovaps 0x040(%%rsi),%%zmm10	\n\t"\
		"vmovaps 	0x400(%%rsi),%%zmm1			\n\t	vmovaps 0x440(%%rsi),%%zmm11	\n\t"\
		"vmovaps 	0x800(%%rsi),%%zmm2			\n\t	vmovaps 0x840(%%rsi),%%zmm12	\n\t"\
		"vmovaps 	0xc00(%%rsi),%%zmm3			\n\t	vmovaps 0xc40(%%rsi),%%zmm13	\n\t"\
		"vmovaps 	0x200(%%rsi),%%zmm4			\n\t	vmovaps 0x240(%%rsi),%%zmm14	\n\t"\
		"vmovaps 	0x600(%%rsi),%%zmm5			\n\t	vmovaps 0x640(%%rsi),%%zmm15	\n\t"\
		"vmovaps 	0xa00(%%rsi),%%zmm6			\n\t	vmovaps 0xa40(%%rsi),%%zmm16	\n\t"\
		"vmovaps 	0xe00(%%rsi),%%zmm7			\n\t	vmovaps 0xe40(%%rsi),%%zmm17	\n\t"\
		/* [1] Shuffle the 4-aparts - note the different patterning of the first and second output quartet: */\
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
		/* Write resulting 8 vec-cmplx data to eight 512[0x200]-byte-separated slots of main data array: */\
		"vmovaps		%%zmm5,(%%rax)		\n\t	vmovaps	%%zmm15,0x40(%%rax)	\n\t"\
		"vmovaps		%%zmm6,(%%rbx)		\n\t	vmovaps	%%zmm16,0x40(%%rbx)	\n\t"\
		"vmovaps		%%zmm8,(%%rcx)		\n\t	vmovaps	%%zmm9 ,0x40(%%rcx)	\n\t"\
		"vmovaps		%%zmm0,(%%rdx)		\n\t	vmovaps	%%zmm10,0x40(%%rdx)	\n\t"\
		"vmovaps		%%zmm1,(%%r10)		\n\t	vmovaps	%%zmm11,0x40(%%r10)	\n\t"\
		"vmovaps		%%zmm2,(%%r11)		\n\t	vmovaps	%%zmm12,0x40(%%r11)	\n\t"\
		"vmovaps		%%zmm4,(%%r12)		\n\t	vmovaps	%%zmm14,0x40(%%r12)	\n\t"\
		"vmovaps		%%zmm7,(%%r13)		\n\t	vmovaps	%%zmm17,0x40(%%r13)	\n\t"\
	/**** a[16-23] = re[ 2,10,18,26, 6,14,22,30].d0	a[24-31] = im[ 2,10,18,26, 6,14,22,30].d0 : ****/\
	"addq	$0x80 ,%%rax	\n\t"\
	"addq	$0x80 ,%%rbx	\n\t"\
	"addq	$0x80 ,%%rcx	\n\t"\
	"addq	$0x80 ,%%rdx	\n\t"\
	"addq	$0x80 ,%%r10	\n\t"\
	"addq	$0x80 ,%%r11	\n\t"\
	"addq	$0x80 ,%%r12	\n\t"\
	"addq	$0x80 ,%%r13	\n\t"\
	"addq	$0x100,%%rsi	\n\t"\
		/* Read transposed row-data from the local store: */\
		"vmovaps 		 (%%rsi),%%zmm0			\n\t	vmovaps 0x040(%%rsi),%%zmm10	\n\t"\
		"vmovaps 	0x400(%%rsi),%%zmm1			\n\t	vmovaps 0x440(%%rsi),%%zmm11	\n\t"\
		"vmovaps 	0x800(%%rsi),%%zmm2			\n\t	vmovaps 0x840(%%rsi),%%zmm12	\n\t"\
		"vmovaps 	0xc00(%%rsi),%%zmm3			\n\t	vmovaps 0xc40(%%rsi),%%zmm13	\n\t"\
		"vmovaps 	0x200(%%rsi),%%zmm4			\n\t	vmovaps 0x240(%%rsi),%%zmm14	\n\t"\
		"vmovaps 	0x600(%%rsi),%%zmm5			\n\t	vmovaps 0x640(%%rsi),%%zmm15	\n\t"\
		"vmovaps 	0xa00(%%rsi),%%zmm6			\n\t	vmovaps 0xa40(%%rsi),%%zmm16	\n\t"\
		"vmovaps 	0xe00(%%rsi),%%zmm7			\n\t	vmovaps 0xe40(%%rsi),%%zmm17	\n\t"\
		/* [1] Shuffle the 4-aparts - note the different patterning of the first and second output quartet: */\
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
		/* Write resulting 8 vec-cmplx data to eight 512[0x200]-byte-separated slots of main data array: */\
		"vmovaps		%%zmm5,(%%rax)		\n\t	vmovaps	%%zmm15,0x40(%%rax)	\n\t"\
		"vmovaps		%%zmm6,(%%rbx)		\n\t	vmovaps	%%zmm16,0x40(%%rbx)	\n\t"\
		"vmovaps		%%zmm8,(%%rcx)		\n\t	vmovaps	%%zmm9 ,0x40(%%rcx)	\n\t"\
		"vmovaps		%%zmm0,(%%rdx)		\n\t	vmovaps	%%zmm10,0x40(%%rdx)	\n\t"\
		"vmovaps		%%zmm1,(%%r10)		\n\t	vmovaps	%%zmm11,0x40(%%r10)	\n\t"\
		"vmovaps		%%zmm2,(%%r11)		\n\t	vmovaps	%%zmm12,0x40(%%r11)	\n\t"\
		"vmovaps		%%zmm4,(%%r12)		\n\t	vmovaps	%%zmm14,0x40(%%r12)	\n\t"\
		"vmovaps		%%zmm7,(%%r13)		\n\t	vmovaps	%%zmm17,0x40(%%r13)	\n\t"\
	/**** a[32-39] = re[ 1, 9,17,25, 5,13,21,29].d0	a[40-47] = im[ 1, 9,17,25, 5,13,21,29].d0 : ****/\
	"addq	$0x80 ,%%rax	\n\t"\
	"addq	$0x80 ,%%rbx	\n\t"\
	"addq	$0x80 ,%%rcx	\n\t"\
	"addq	$0x80 ,%%rdx	\n\t"\
	"addq	$0x80 ,%%r10	\n\t"\
	"addq	$0x80 ,%%r11	\n\t"\
	"addq	$0x80 ,%%r12	\n\t"\
	"addq	$0x80 ,%%r13	\n\t"\
	"subq	$0x80,%%rsi	\n\t"\
		/* Read transposed row-data from the local store: */\
		"vmovaps 		 (%%rsi),%%zmm0			\n\t	vmovaps 0x040(%%rsi),%%zmm10	\n\t"\
		"vmovaps 	0x400(%%rsi),%%zmm1			\n\t	vmovaps 0x440(%%rsi),%%zmm11	\n\t"\
		"vmovaps 	0x800(%%rsi),%%zmm2			\n\t	vmovaps 0x840(%%rsi),%%zmm12	\n\t"\
		"vmovaps 	0xc00(%%rsi),%%zmm3			\n\t	vmovaps 0xc40(%%rsi),%%zmm13	\n\t"\
		"vmovaps 	0x200(%%rsi),%%zmm4			\n\t	vmovaps 0x240(%%rsi),%%zmm14	\n\t"\
		"vmovaps 	0x600(%%rsi),%%zmm5			\n\t	vmovaps 0x640(%%rsi),%%zmm15	\n\t"\
		"vmovaps 	0xa00(%%rsi),%%zmm6			\n\t	vmovaps 0xa40(%%rsi),%%zmm16	\n\t"\
		"vmovaps 	0xe00(%%rsi),%%zmm7			\n\t	vmovaps 0xe40(%%rsi),%%zmm17	\n\t"\
		/* [1] Shuffle the 4-aparts - note the different patterning of the first and second output quartet: */\
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
		/* Write resulting 8 vec-cmplx data to eight 512[0x200]-byte-separated slots of main data array: */\
		"vmovaps		%%zmm5,(%%rax)		\n\t	vmovaps	%%zmm15,0x40(%%rax)	\n\t"\
		"vmovaps		%%zmm6,(%%rbx)		\n\t	vmovaps	%%zmm16,0x40(%%rbx)	\n\t"\
		"vmovaps		%%zmm8,(%%rcx)		\n\t	vmovaps	%%zmm9 ,0x40(%%rcx)	\n\t"\
		"vmovaps		%%zmm0,(%%rdx)		\n\t	vmovaps	%%zmm10,0x40(%%rdx)	\n\t"\
		"vmovaps		%%zmm1,(%%r10)		\n\t	vmovaps	%%zmm11,0x40(%%r10)	\n\t"\
		"vmovaps		%%zmm2,(%%r11)		\n\t	vmovaps	%%zmm12,0x40(%%r11)	\n\t"\
		"vmovaps		%%zmm4,(%%r12)		\n\t	vmovaps	%%zmm14,0x40(%%r12)	\n\t"\
		"vmovaps		%%zmm7,(%%r13)		\n\t	vmovaps	%%zmm17,0x40(%%r13)	\n\t"\
	/**** a[48-55] = re[ 3,11,19,27, 7,15,23,31].d0	a[56-63] = im[ 3,11,19,27, 7,15,23,31].d0 : ****/\
	"addq	$0x80 ,%%rax	\n\t"\
	"addq	$0x80 ,%%rbx	\n\t"\
	"addq	$0x80 ,%%rcx	\n\t"\
	"addq	$0x80 ,%%rdx	\n\t"\
	"addq	$0x80 ,%%r10	\n\t"\
	"addq	$0x80 ,%%r11	\n\t"\
	"addq	$0x80 ,%%r12	\n\t"\
	"addq	$0x80 ,%%r13	\n\t"\
	"addq	$0x100,%%rsi	\n\t"\
		/* Read transposed row-data from the local store: */\
		"vmovaps 		 (%%rsi),%%zmm0			\n\t	vmovaps 0x040(%%rsi),%%zmm10	\n\t"\
		"vmovaps 	0x400(%%rsi),%%zmm1			\n\t	vmovaps 0x440(%%rsi),%%zmm11	\n\t"\
		"vmovaps 	0x800(%%rsi),%%zmm2			\n\t	vmovaps 0x840(%%rsi),%%zmm12	\n\t"\
		"vmovaps 	0xc00(%%rsi),%%zmm3			\n\t	vmovaps 0xc40(%%rsi),%%zmm13	\n\t"\
		"vmovaps 	0x200(%%rsi),%%zmm4			\n\t	vmovaps 0x240(%%rsi),%%zmm14	\n\t"\
		"vmovaps 	0x600(%%rsi),%%zmm5			\n\t	vmovaps 0x640(%%rsi),%%zmm15	\n\t"\
		"vmovaps 	0xa00(%%rsi),%%zmm6			\n\t	vmovaps 0xa40(%%rsi),%%zmm16	\n\t"\
		"vmovaps 	0xe00(%%rsi),%%zmm7			\n\t	vmovaps 0xe40(%%rsi),%%zmm17	\n\t"\
		/* [1] Shuffle the 4-aparts - note the different patterning of the first and second output quartet: */\
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
		/* Write resulting 8 vec-cmplx data to eight 512[0x200]-byte-separated slots of main data array: */\
		"vmovaps		%%zmm5,(%%rax)		\n\t	vmovaps	%%zmm15,0x40(%%rax)	\n\t"\
		"vmovaps		%%zmm6,(%%rbx)		\n\t	vmovaps	%%zmm16,0x40(%%rbx)	\n\t"\
		"vmovaps		%%zmm8,(%%rcx)		\n\t	vmovaps	%%zmm9 ,0x40(%%rcx)	\n\t"\
		"vmovaps		%%zmm0,(%%rdx)		\n\t	vmovaps	%%zmm10,0x40(%%rdx)	\n\t"\
		"vmovaps		%%zmm1,(%%r10)		\n\t	vmovaps	%%zmm11,0x40(%%r10)	\n\t"\
		"vmovaps		%%zmm2,(%%r11)		\n\t	vmovaps	%%zmm12,0x40(%%r11)	\n\t"\
		"vmovaps		%%zmm4,(%%r12)		\n\t	vmovaps	%%zmm14,0x40(%%r12)	\n\t"\
		"vmovaps		%%zmm7,(%%r13)		\n\t	vmovaps	%%zmm17,0x40(%%r13)	\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__add1] "m" (Xadd1)\
		 ,[__add2] "m" (Xadd2)\
		 ,[__add3] "m" (Xadd3)\
		 ,[__add4] "m" (Xadd4)\
		 ,[__add5] "m" (Xadd5)\
		 ,[__add6] "m" (Xadd6)\
		 ,[__add7] "m" (Xadd7)\
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__r00] "m" (Xr00)\
		 ,[__r08] "m" (Xr08)\
		 ,[__r10] "m" (Xr10)\
		 ,[__r20] "m" (Xr20)\
		 ,[__r28] "m" (Xr28)\
		 ,[__r30] "m" (Xr30)\
		 ,[__c01] "m" (Xc01)\
		 ,[__c02] "m" (Xc02)\
		 ,[__c04] "m" (Xc04)\
		 ,[__c06] "m" (Xc06)\
		 ,[__c0A] "m" (Xc0A)\
		 ,[__c10] "m" (Xc10)\
		 ,[__c12] "m" (Xc12)\
		 ,[__c1A] "m" (Xc1A)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r10","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15","xmm16","xmm17","xmm18","xmm19","xmm20","xmm21","xmm22","xmm23","xmm24","xmm25","xmm26","xmm27","xmm28","xmm29","xmm30","xmm31"	/* Clobbered registers */\
	);\
	}

   #else	// AVX-512 version:

	// Cost [vector-ops only] : 528 MEM (400 in DFT proper), 370 MUL (270 FMA), 168 ADD/SUB	<*** 32-DIT MEM-ops reduced by more than half vs 16-reg version ***
	// Compare to AVX2 version: 992 MEM (864 in DFT proper), 352 MUL (236 FMA), 204 ADD/SUB
	#define SSE2_RADIX32_WRAPPER_DIT(Xadd0,Xadd1,Xadd2,Xadd3,Xadd4,Xadd5,Xadd6,Xadd7,Xisrt2,Xr00,Xr08,Xr10,Xr20,Xr28,Xr30,Xc01,Xc02,Xc04,Xc06,Xc0A,Xc10,Xc12,Xc1A)\
	{\
	__asm__ volatile (\
	/************************************************************************/\
	/*...And do an inverse DIT radix-32 pass on the squared-data blocks.	*/\
	/************************************************************************/\
	/* Block 1: */\
		"movq	%[__isrt2],%%rsi			\n\t	leaq	0x11c0(%%rsi),%%rdi		\n\t"/* two */\
		"movq		%[__r00],%%rax			\n\t	vmovaps		 (%%rdi),%%zmm29	\n\t"/* 2.0*/	/* Next pair of 4-DFTs - Place these results in regs whose index is += 16 w.r.to the above 4-DFTs: */\
	/**** SSE2_RADIX4_DIT_IN_PLACE() ****/	/* Block 2 has addr-offsets +0x100 w.r.to Block 1: */	/* Byte-offsets all +0x200 larger in this 4-DFT pair */\
		"vmovaps	0x800(%%rax),%%zmm0		\n\t	vmovaps		0x900(%%rax),%%zmm8 		\n\t	vmovaps	0xa00(%%rax),%%zmm16		\n\t	vmovaps		0xb00(%%rax),%%zmm24		\n\t"\
		"vmovaps	0x840(%%rax),%%zmm1		\n\t	vmovaps		0x940(%%rax),%%zmm9 		\n\t	vmovaps	0xa40(%%rax),%%zmm17		\n\t	vmovaps		0xb40(%%rax),%%zmm25		\n\t"\
		"vmovaps		 (%%rax),%%zmm2		\n\t	vmovaps		0x100(%%rax),%%zmm10		\n\t	vmovaps	0x200(%%rax),%%zmm18		\n\t	vmovaps		0x300(%%rax),%%zmm26		\n\t"\
		"vmovaps	0x040(%%rax),%%zmm3		\n\t	vmovaps		0x140(%%rax),%%zmm11		\n\t	vmovaps	0x240(%%rax),%%zmm19		\n\t	vmovaps		0x340(%%rax),%%zmm27		\n\t"\
		"vsubpd		%%zmm0,%%zmm2,%%zmm2	\n\t	vsubpd		%%zmm8 ,%%zmm10,%%zmm10		\n\t	vsubpd		%%zmm16,%%zmm18,%%zmm18	\n\t	vsubpd		%%zmm24,%%zmm26,%%zmm26		\n\t"\
		"vsubpd		%%zmm1,%%zmm3,%%zmm3	\n\t	vsubpd		%%zmm9 ,%%zmm11,%%zmm11		\n\t	vsubpd		%%zmm17,%%zmm19,%%zmm19	\n\t	vsubpd		%%zmm25,%%zmm27,%%zmm27		\n\t"\
		"vmovaps	0xc00(%%rax),%%zmm4		\n\t	vmovaps		0xd00(%%rax),%%zmm12		\n\t	vmovaps	0xe00(%%rax),%%zmm20		\n\t	vmovaps		0xf00(%%rax),%%zmm28		\n\t"\
		"vmovaps	0xc40(%%rax),%%zmm5		\n\t	vmovaps		0xd40(%%rax),%%zmm13		\n\t	vmovaps	0xe40(%%rax),%%zmm21		\n\t /* vmovaps		0xf40(%%rax),%%zmm29	*/	\n\t"\
		"vmovaps	0x400(%%rax),%%zmm6		\n\t	vmovaps		0x500(%%rax),%%zmm14		\n\t	vmovaps	0x600(%%rax),%%zmm22		\n\t	vmovaps		0x700(%%rax),%%zmm30		\n\t"\
		"vmovaps	0x440(%%rax),%%zmm7		\n\t	vmovaps		0x540(%%rax),%%zmm15		\n\t	vmovaps	0x640(%%rax),%%zmm23		\n\t	vmovaps		0x740(%%rax),%%zmm31		\n\t"\
		"vsubpd		%%zmm4,%%zmm6,%%zmm6	\n\t	vsubpd		%%zmm12,%%zmm14,%%zmm14		\n\t	vsubpd		%%zmm20,%%zmm22,%%zmm22	\n\t	vsubpd		%%zmm28,%%zmm30,%%zmm30		\n\t"\
		"vsubpd		%%zmm5,%%zmm7,%%zmm7	\n\t	vsubpd		%%zmm13,%%zmm15,%%zmm15		\n\t	vsubpd		%%zmm21,%%zmm23,%%zmm23	\n\t	vsubpd 0xf40(%%rax),%%zmm31,%%zmm31		\n\t"\
	"vfmadd132pd	%%zmm29,%%zmm2,%%zmm0	\n\t vfmadd132pd	%%zmm29,%%zmm10,%%zmm8 		\n\t vfmadd132pd	%%zmm29,%%zmm18,%%zmm16	\n\t vfmadd132pd	%%zmm29,%%zmm26,%%zmm24		\n\t"\
	"vfmadd132pd	%%zmm29,%%zmm3,%%zmm1	\n\t vfmadd132pd	%%zmm29,%%zmm11,%%zmm9 		\n\t vfmadd132pd	%%zmm29,%%zmm19,%%zmm17	\n\t vfmadd132pd	%%zmm29,%%zmm27,%%zmm25		\n\t"\
	"vfmadd132pd	%%zmm29,%%zmm6,%%zmm4	\n\t vfmadd132pd	%%zmm29,%%zmm14,%%zmm12		\n\t vfmadd132pd	%%zmm29,%%zmm22,%%zmm20	\n\t vfmadd132pd	%%zmm29,%%zmm30,%%zmm28		\n\t"\
	"vfmadd132pd	%%zmm29,%%zmm7,%%zmm5	\n\t vfmadd132pd	%%zmm29,%%zmm15,%%zmm13		\n\t vfmadd132pd	%%zmm29,%%zmm23,%%zmm21	\n\t vfmadd132pd 0xf40(%%rax),%%zmm31,%%zmm29	\n\t"\
		"vsubpd		%%zmm4 ,%%zmm0,%%zmm0	\n\t	vsubpd		%%zmm12,%%zmm8 ,%%zmm8 		\n\t	vsubpd		%%zmm20,%%zmm16,%%zmm16	\n\t	vsubpd		%%zmm28,%%zmm24,%%zmm24		\n\t"\
		"vsubpd		%%zmm5 ,%%zmm1,%%zmm1	\n\t	vsubpd		%%zmm13,%%zmm9 ,%%zmm9 		\n\t	vsubpd		%%zmm21,%%zmm17,%%zmm17	\n\t	vsubpd		%%zmm29,%%zmm25,%%zmm25		\n\t"\
		"vsubpd		%%zmm7 ,%%zmm2,%%zmm2	\n\t	vsubpd		%%zmm15,%%zmm10,%%zmm10		\n\t	vsubpd		%%zmm23,%%zmm18,%%zmm18	\n\t	vsubpd		%%zmm31,%%zmm26,%%zmm26		\n\t"\
		"vsubpd		%%zmm6 ,%%zmm3,%%zmm3	\n\t	vsubpd		%%zmm14,%%zmm11,%%zmm11		\n\t	vsubpd		%%zmm22,%%zmm19,%%zmm19	\n\t	vsubpd		%%zmm30,%%zmm27,%%zmm27		\n\t"\
	"vmovaps %%zmm30,0xf40(%%rax)\n\t"/* spill zmm30 to make room for two */"vmovaps (%%rdi),%%zmm30\n\t"\
	"vfmadd132pd	%%zmm30,%%zmm0,%%zmm4	\n\t vfmadd132pd	%%zmm30,%%zmm8 ,%%zmm12		\n\t vfmadd132pd	%%zmm30,%%zmm16,%%zmm20	\n\t vfmadd132pd	%%zmm30,%%zmm24,%%zmm28		\n\t"\
	"vfmadd132pd	%%zmm30,%%zmm1,%%zmm5	\n\t vfmadd132pd	%%zmm30,%%zmm9 ,%%zmm13		\n\t vfmadd132pd	%%zmm30,%%zmm17,%%zmm21	\n\t vfmadd132pd	%%zmm30,%%zmm25,%%zmm29		\n\t"\
	"vfmadd132pd	%%zmm30,%%zmm2,%%zmm7	\n\t vfmadd132pd	%%zmm30,%%zmm10,%%zmm15		\n\t vfmadd132pd	%%zmm30,%%zmm18,%%zmm23	\n\t vfmadd132pd	%%zmm30,%%zmm26,%%zmm31		\n\t"\
	"vfmadd132pd	%%zmm30,%%zmm3,%%zmm6	\n\t vfmadd132pd	%%zmm30,%%zmm11,%%zmm14		\n\t vfmadd132pd	%%zmm30,%%zmm19,%%zmm22	\n\t vfmadd132pd 0xf40(%%rax),%%zmm27,%%zmm30	\n\t"\
	/**** SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS ****//**** SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS ****/	/* These are actually the final steps of the col3/4 4-DFTs, can do alongside */\
	/* 0x200(r00,r10,r20,r30,r08,r18,r28,r38): *//* 0x200(r04,r14,r24,r34,r0C,r1C,r2C,r3C): */	/* first steps of the combine-subs in cols 1/2 due to nonoverlapping registers: */\
		"vsubpd		%%zmm20,%%zmm4,%%zmm4	\n\t	vsubpd		%%zmm28,%%zmm12,%%zmm12		\n\t	vaddpd		%%zmm23,%%zmm19,%%zmm23	\n\t	vaddpd		%%zmm31,%%zmm27,%%zmm31		\n\t"\
		"vsubpd		%%zmm17,%%zmm0,%%zmm0	\n\t	vsubpd		%%zmm25,%%zmm8 ,%%zmm8 		\n\t	vaddpd		%%zmm18,%%zmm22,%%zmm18	\n\t	vaddpd		%%zmm26,%%zmm30,%%zmm26		\n\t"\
		"vsubpd		%%zmm21,%%zmm5,%%zmm5	\n\t	vsubpd		%%zmm29,%%zmm13,%%zmm13		\n\t vfmsub132pd	(%%rdi),%%zmm23,%%zmm19	\n\t vfmsub132pd	(%%rdi),%%zmm31,%%zmm27		\n\t"\
		"vsubpd		%%zmm16,%%zmm1,%%zmm1	\n\t	vsubpd		%%zmm24,%%zmm9 ,%%zmm9 		\n\t vfmsub132pd	(%%rdi),%%zmm18,%%zmm22	\n\t vfmsub132pd	(%%rdi),%%zmm26,%%zmm30		\n\t"\
		"vmovaps	%%zmm4 , 0x200(%%rax)	\n\t	vmovaps		%%zmm12, 0x300(%%rax)		\n\t"\
		"vmovaps	%%zmm0 , 0xa00(%%rax)	\n\t	vmovaps		%%zmm8 , 0xb00(%%rax)		\n\t"\
		"vmovaps	%%zmm5 , 0x240(%%rax)	\n\t	vmovaps		%%zmm13, 0x340(%%rax)		\n\t"\
		"vmovaps	%%zmm1 , 0x840(%%rax)	\n\t	vmovaps		%%zmm9 , 0x940(%%rax)		\n\t"\
	"vmovaps %%zmm24,0xb40(%%rax)\n\t"/* spill zmm24to make room for two */"vmovaps (%%rdi),%%zmm24\n\t"\
	"vfmadd132pd	%%zmm24,%%zmm4,%%zmm20	\n\t vfmadd132pd	 %%zmm24,%%zmm12,%%zmm28	\n\t"\
	"vfmadd132pd	%%zmm24,%%zmm0,%%zmm17	\n\t vfmadd132pd	 %%zmm24,%%zmm8 ,%%zmm25	\n\t"\
	"vfmadd132pd	%%zmm24,%%zmm5,%%zmm21	\n\t vfmadd132pd	 %%zmm24,%%zmm13,%%zmm29	\n\t"\
	"vfmadd132pd	%%zmm24,%%zmm1,%%zmm16	\n\t vfmadd132pd 0xb40(%%rax),%%zmm9 ,%%zmm24	\n\t"\
		"vmovaps	%%zmm20,      (%%rax)	\n\t	vmovaps		%%zmm28, 0x100(%%rax)		\n\t"\
		"vmovaps	%%zmm17, 0x800(%%rax)	\n\t	vmovaps		%%zmm25, 0x900(%%rax)		\n\t"\
		"vmovaps	%%zmm21, 0x040(%%rax)	\n\t	vmovaps		%%zmm29, 0x140(%%rax)		\n\t"\
		"vmovaps	%%zmm16, 0xa40(%%rax)	\n\t	vmovaps		%%zmm24, 0xb40(%%rax)		\n\t"\
	"vmovaps %%zmm14,0xf40(%%rax)\n\t"/* spill zmm14to make room for isrt2 */"vmovaps (%%rsi),%%zmm14\n\t"\
	"vfnmadd231pd	%%zmm23,%%zmm14,%%zmm7	\n\t vfnmadd231pd	%%zmm31,%%zmm14,%%zmm15		\n\t"\
	"vfnmadd231pd	%%zmm22,%%zmm14,%%zmm2	\n\t vfnmadd231pd	%%zmm30,%%zmm14,%%zmm10		\n\t"\
	"vfnmadd231pd	%%zmm19,%%zmm14,%%zmm3	\n\t vfnmadd231pd	%%zmm27,%%zmm14,%%zmm11		\n\t"\
	"vfnmadd231pd	%%zmm18,%%zmm14,%%zmm6	\n\t vfnmadd213pd 0xf40(%%rax),%%zmm26,%%zmm14	\n\t"\
		"vmovaps	%%zmm7 , 0x600(%%rax)	\n\t	vmovaps		%%zmm15, 0x700(%%rax)		\n\t"\
		"vmovaps	%%zmm2 , 0xe00(%%rax)	\n\t	vmovaps		%%zmm10, 0xf00(%%rax)		\n\t"\
		"vmovaps	%%zmm3 , 0x640(%%rax)	\n\t	vmovaps		%%zmm11, 0x740(%%rax)		\n\t"\
		"vmovaps	%%zmm6 , 0xc40(%%rax)	\n\t	vmovaps		%%zmm14, 0xd40(%%rax)		\n\t"\
	"vmovaps %%zmm26,0xf40(%%rax)\n\t"/* spill zmm26to make room for sqrt2 */"vmovaps 0x80(%%rdi),%%zmm26\n\t"\
	"vfmadd132pd	%%zmm26,%%zmm7,%%zmm23	\n\t vfmadd132pd	 %%zmm26,%%zmm15,%%zmm31	\n\t"\
	"vfmadd132pd	%%zmm26,%%zmm2,%%zmm22	\n\t vfmadd132pd	 %%zmm26,%%zmm10,%%zmm30	\n\t"\
	"vfmadd132pd	%%zmm26,%%zmm3,%%zmm19	\n\t vfmadd132pd	 %%zmm26,%%zmm11,%%zmm27	\n\t"\
	"vfmadd132pd	%%zmm26,%%zmm6,%%zmm18	\n\t vfmadd132pd 0xf40(%%rax),%%zmm14,%%zmm26	\n\t"\
		"vmovaps	%%zmm23, 0x400(%%rax)	\n\t	vmovaps		%%zmm31, 0x500(%%rax)		\n\t"\
		"vmovaps	%%zmm22, 0xc00(%%rax)	\n\t	vmovaps		%%zmm30, 0xd00(%%rax)		\n\t"\
		"vmovaps	%%zmm19, 0x440(%%rax)	\n\t	vmovaps		%%zmm27, 0x540(%%rax)		\n\t"\
		"vmovaps	%%zmm18, 0xe40(%%rax)	\n\t	vmovaps		%%zmm26, 0xf40(%%rax)		\n\t"\
	/*****	SSE2_RADIX4_DIT_IN_PLACE()	*****/	/* Rcol has tmp-addr offset +0x100 w.r.to lcol:  */\
	/* Blocks 3,4 have tmp-addresses offset +0x80 w.r.to Blocks 1,2, respectively. */\
		"											vmovaps		 (%%rdi),%%zmm29	\n\t"/* 2.0*/\
		"vmovaps	0x880(%%rax),%%zmm0		\n\t	vmovaps		0x980(%%rax),%%zmm8 		\n\t	vmovaps	0xa80(%%rax),%%zmm16		\n\t	vmovaps		0xb80(%%rax),%%zmm24		\n\t"\
		"vmovaps	0x8c0(%%rax),%%zmm1		\n\t	vmovaps		0x9c0(%%rax),%%zmm9 		\n\t	vmovaps	0xac0(%%rax),%%zmm17		\n\t	vmovaps		0xbc0(%%rax),%%zmm25		\n\t"\
		"vmovaps	0x080(%%rax),%%zmm2		\n\t	vmovaps		0x180(%%rax),%%zmm10		\n\t	vmovaps	0x280(%%rax),%%zmm18		\n\t	vmovaps		0x380(%%rax),%%zmm26		\n\t"\
		"vmovaps	0x0c0(%%rax),%%zmm3		\n\t	vmovaps		0x1c0(%%rax),%%zmm11		\n\t	vmovaps	0x2c0(%%rax),%%zmm19		\n\t	vmovaps		0x3c0(%%rax),%%zmm27		\n\t"\
		"vsubpd		%%zmm0,%%zmm2,%%zmm2	\n\t	vsubpd		%%zmm8 ,%%zmm10,%%zmm10		\n\t	vsubpd		%%zmm16,%%zmm18,%%zmm18	\n\t	vsubpd		%%zmm24,%%zmm26,%%zmm26		\n\t"\
		"vsubpd		%%zmm1,%%zmm3,%%zmm3	\n\t	vsubpd		%%zmm9 ,%%zmm11,%%zmm11		\n\t	vsubpd		%%zmm17,%%zmm19,%%zmm19	\n\t	vsubpd		%%zmm25,%%zmm27,%%zmm27		\n\t"\
		"vmovaps	0xc80(%%rax),%%zmm4		\n\t	vmovaps		0xd80(%%rax),%%zmm12		\n\t	vmovaps	0xe80(%%rax),%%zmm20		\n\t	vmovaps		0xf80(%%rax),%%zmm28		\n\t"\
		"vmovaps	0xcc0(%%rax),%%zmm5		\n\t	vmovaps		0xdc0(%%rax),%%zmm13		\n\t	vmovaps	0xec0(%%rax),%%zmm21		\n\t /* vmovaps		0xfc0(%%rax),%%zmm29	*/	\n\t"\
		"vmovaps	0x480(%%rax),%%zmm6		\n\t	vmovaps		0x580(%%rax),%%zmm14		\n\t	vmovaps	0x680(%%rax),%%zmm22		\n\t	vmovaps		0x780(%%rax),%%zmm30		\n\t"\
		"vmovaps	0x4c0(%%rax),%%zmm7		\n\t	vmovaps		0x5c0(%%rax),%%zmm15		\n\t	vmovaps	0x6c0(%%rax),%%zmm23		\n\t	vmovaps		0x7c0(%%rax),%%zmm31		\n\t"\
		"vsubpd		%%zmm4,%%zmm6,%%zmm6	\n\t	vsubpd		%%zmm12,%%zmm14,%%zmm14		\n\t	vsubpd		%%zmm20,%%zmm22,%%zmm22	\n\t	vsubpd		%%zmm28,%%zmm30,%%zmm30		\n\t"\
		"vsubpd		%%zmm5,%%zmm7,%%zmm7	\n\t	vsubpd		%%zmm13,%%zmm15,%%zmm15		\n\t	vsubpd		%%zmm21,%%zmm23,%%zmm23	\n\t	vsubpd 0xfc0(%%rax),%%zmm31,%%zmm31		\n\t"\
	"vfmadd132pd	%%zmm29,%%zmm2,%%zmm0	\n\t vfmadd132pd	%%zmm29,%%zmm10,%%zmm8 		\n\t vfmadd132pd	%%zmm29,%%zmm18,%%zmm16	\n\t vfmadd132pd	%%zmm29,%%zmm26,%%zmm24		\n\t"\
	"vfmadd132pd	%%zmm29,%%zmm3,%%zmm1	\n\t vfmadd132pd	%%zmm29,%%zmm11,%%zmm9 		\n\t vfmadd132pd	%%zmm29,%%zmm19,%%zmm17	\n\t vfmadd132pd	%%zmm29,%%zmm27,%%zmm25		\n\t"\
	"vfmadd132pd	%%zmm29,%%zmm6,%%zmm4	\n\t vfmadd132pd	%%zmm29,%%zmm14,%%zmm12		\n\t vfmadd132pd	%%zmm29,%%zmm22,%%zmm20	\n\t vfmadd132pd	%%zmm29,%%zmm30,%%zmm28		\n\t"\
	"vfmadd132pd	%%zmm29,%%zmm7,%%zmm5	\n\t vfmadd132pd	%%zmm29,%%zmm15,%%zmm13		\n\t vfmadd132pd	%%zmm29,%%zmm23,%%zmm21	\n\t vfmadd132pd 0xfc0(%%rax),%%zmm31,%%zmm29	\n\t"\
		"vsubpd		%%zmm4 ,%%zmm0,%%zmm0	\n\t	vsubpd		%%zmm12,%%zmm8 ,%%zmm8 		\n\t	vsubpd		%%zmm20,%%zmm16,%%zmm16	\n\t	vsubpd		%%zmm28,%%zmm24,%%zmm24		\n\t"\
		"vsubpd		%%zmm5 ,%%zmm1,%%zmm1	\n\t	vsubpd		%%zmm13,%%zmm9 ,%%zmm9 		\n\t	vsubpd		%%zmm21,%%zmm17,%%zmm17	\n\t	vsubpd		%%zmm29,%%zmm25,%%zmm25		\n\t"\
		"vsubpd		%%zmm7 ,%%zmm2,%%zmm2	\n\t	vsubpd		%%zmm15,%%zmm10,%%zmm10		\n\t	vsubpd		%%zmm23,%%zmm18,%%zmm18	\n\t	vsubpd		%%zmm31,%%zmm26,%%zmm26		\n\t"\
		"vsubpd		%%zmm6 ,%%zmm3,%%zmm3	\n\t	vsubpd		%%zmm14,%%zmm11,%%zmm11		\n\t	vsubpd		%%zmm22,%%zmm19,%%zmm19	\n\t	vsubpd		%%zmm30,%%zmm27,%%zmm27		\n\t"\
	"vmovaps %%zmm30,0xfc0(%%rax)\n\t"/* spill zmm30 to make room for two */"vmovaps (%%rdi),%%zmm30\n\t"\
	"vfmadd132pd	%%zmm30,%%zmm0,%%zmm4	\n\t vfmadd132pd	%%zmm30,%%zmm8 ,%%zmm12		\n\t vfmadd132pd	%%zmm30,%%zmm16,%%zmm20	\n\t vfmadd132pd	%%zmm30,%%zmm24,%%zmm28		\n\t"\
	"vfmadd132pd	%%zmm30,%%zmm1,%%zmm5	\n\t vfmadd132pd	%%zmm30,%%zmm9 ,%%zmm13		\n\t vfmadd132pd	%%zmm30,%%zmm17,%%zmm21	\n\t vfmadd132pd	%%zmm30,%%zmm25,%%zmm29		\n\t"\
	"vfmadd132pd	%%zmm30,%%zmm2,%%zmm7	\n\t vfmadd132pd	%%zmm30,%%zmm10,%%zmm15		\n\t vfmadd132pd	%%zmm30,%%zmm18,%%zmm23	\n\t vfmadd132pd	%%zmm30,%%zmm26,%%zmm31		\n\t"\
	"vfmadd132pd	%%zmm30,%%zmm3,%%zmm6	\n\t vfmadd132pd	%%zmm30,%%zmm11,%%zmm14		\n\t vfmadd132pd	%%zmm30,%%zmm19,%%zmm22	\n\t vfmadd132pd 0xfc0(%%rax),%%zmm27,%%zmm30	\n\t"\
	/**** SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS ****//**** SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS ****/	/* These are actually the final steps of the col3/4 4-DFTs, can do alongside */\
	/* 0x200(r00,r10,r20,r30,r08,r18,r28,r38): *//* 0x200(r04,r14,r24,r34,r0C,r1C,r2C,r3C): */	/* first steps of the combine-subs in cols 1/2 due to nonoverlapping registers: */\
		"vsubpd		%%zmm20,%%zmm4,%%zmm4	\n\t	vsubpd		%%zmm28,%%zmm12,%%zmm12		\n\t	vaddpd		%%zmm23,%%zmm19,%%zmm23	\n\t	vaddpd		%%zmm31,%%zmm27,%%zmm31		\n\t"\
		"vsubpd		%%zmm17,%%zmm0,%%zmm0	\n\t	vsubpd		%%zmm25,%%zmm8 ,%%zmm8 		\n\t	vaddpd		%%zmm18,%%zmm22,%%zmm18	\n\t	vaddpd		%%zmm26,%%zmm30,%%zmm26		\n\t"\
		"vsubpd		%%zmm21,%%zmm5,%%zmm5	\n\t	vsubpd		%%zmm29,%%zmm13,%%zmm13		\n\t vfmsub132pd	(%%rdi),%%zmm23,%%zmm19	\n\t vfmsub132pd	(%%rdi),%%zmm31,%%zmm27		\n\t"\
		"vsubpd		%%zmm16,%%zmm1,%%zmm1	\n\t	vsubpd		%%zmm24,%%zmm9 ,%%zmm9 		\n\t vfmsub132pd	(%%rdi),%%zmm18,%%zmm22	\n\t vfmsub132pd	(%%rdi),%%zmm26,%%zmm30		\n\t"\
		"vmovaps	%%zmm4 , 0x280(%%rax)	\n\t	vmovaps		%%zmm12, 0x380(%%rax)		\n\t"\
		"vmovaps	%%zmm0 , 0xa80(%%rax)	\n\t	vmovaps		%%zmm8 , 0xb80(%%rax)		\n\t"\
		"vmovaps	%%zmm5 , 0x2c0(%%rax)	\n\t	vmovaps		%%zmm13, 0x3c0(%%rax)		\n\t"\
		"vmovaps	%%zmm1 , 0x8c0(%%rax)	\n\t	vmovaps		%%zmm9 , 0x9c0(%%rax)		\n\t"\
	"vmovaps %%zmm24,0xbc0(%%rax)\n\t"/* spill zmm24to make room for two */"vmovaps (%%rdi),%%zmm24\n\t"\
	"vfmadd132pd	%%zmm24,%%zmm4,%%zmm20	\n\t vfmadd132pd	 %%zmm24,%%zmm12,%%zmm28	\n\t"\
	"vfmadd132pd	%%zmm24,%%zmm0,%%zmm17	\n\t vfmadd132pd	 %%zmm24,%%zmm8 ,%%zmm25	\n\t"\
	"vfmadd132pd	%%zmm24,%%zmm5,%%zmm21	\n\t vfmadd132pd	 %%zmm24,%%zmm13,%%zmm29	\n\t"\
	"vfmadd132pd	%%zmm24,%%zmm1,%%zmm16	\n\t vfmadd132pd 0xbc0(%%rax),%%zmm9 ,%%zmm24	\n\t"\
		"vmovaps	%%zmm20, 0x080(%%rax)	\n\t	vmovaps		%%zmm28, 0x180(%%rax)		\n\t"\
		"vmovaps	%%zmm17, 0x880(%%rax)	\n\t	vmovaps		%%zmm25, 0x980(%%rax)		\n\t"\
		"vmovaps	%%zmm21, 0x0c0(%%rax)	\n\t	vmovaps		%%zmm29, 0x1c0(%%rax)		\n\t"\
		"vmovaps	%%zmm16, 0xac0(%%rax)	\n\t	vmovaps		%%zmm24, 0xbc0(%%rax)		\n\t"\
	"vmovaps %%zmm14,0xfc0(%%rax)\n\t"/* spill zmm14to make room for isrt2 */"vmovaps (%%rsi),%%zmm14\n\t"\
	"vfnmadd231pd	%%zmm23,%%zmm14,%%zmm7	\n\t vfnmadd231pd	%%zmm31,%%zmm14,%%zmm15		\n\t"\
	"vfnmadd231pd	%%zmm22,%%zmm14,%%zmm2	\n\t vfnmadd231pd	%%zmm30,%%zmm14,%%zmm10		\n\t"\
	"vfnmadd231pd	%%zmm19,%%zmm14,%%zmm3	\n\t vfnmadd231pd	%%zmm27,%%zmm14,%%zmm11		\n\t"\
	"vfnmadd231pd	%%zmm18,%%zmm14,%%zmm6	\n\t vfnmadd213pd 0xfc0(%%rax),%%zmm26,%%zmm14	\n\t"\
		"vmovaps	%%zmm7 , 0x680(%%rax)	\n\t	vmovaps		%%zmm15, 0x780(%%rax)		\n\t"\
		"vmovaps	%%zmm2 , 0xe80(%%rax)	\n\t	vmovaps		%%zmm10, 0xf80(%%rax)		\n\t"\
		"vmovaps	%%zmm3 , 0x6c0(%%rax)	\n\t	vmovaps		%%zmm11, 0x7c0(%%rax)		\n\t"\
		"vmovaps	%%zmm6 , 0xcc0(%%rax)	\n\t	vmovaps		%%zmm14, 0xdc0(%%rax)		\n\t"\
	"vmovaps %%zmm26,0xfc0(%%rax)\n\t"/* spill zmm26to make room for sqrt2 */"vmovaps 0x80(%%rdi),%%zmm26\n\t"\
	"vfmadd132pd	%%zmm26,%%zmm7,%%zmm23	\n\t vfmadd132pd	 %%zmm26,%%zmm15,%%zmm31	\n\t"\
	"vfmadd132pd	%%zmm26,%%zmm2,%%zmm22	\n\t vfmadd132pd	 %%zmm26,%%zmm10,%%zmm30	\n\t"\
	"vfmadd132pd	%%zmm26,%%zmm3,%%zmm19	\n\t vfmadd132pd	 %%zmm26,%%zmm11,%%zmm27	\n\t"\
	"vfmadd132pd	%%zmm26,%%zmm6,%%zmm18	\n\t vfmadd132pd 0xfc0(%%rax),%%zmm14,%%zmm26	\n\t"\
		"vmovaps	%%zmm23, 0x480(%%rax)	\n\t	vmovaps		%%zmm31, 0x580(%%rax)		\n\t"\
		"vmovaps	%%zmm22, 0xc80(%%rax)	\n\t	vmovaps		%%zmm30, 0xd80(%%rax)		\n\t"\
		"vmovaps	%%zmm19, 0x4c0(%%rax)	\n\t	vmovaps		%%zmm27, 0x5c0(%%rax)		\n\t"\
		"vmovaps	%%zmm18, 0xec0(%%rax)	\n\t	vmovaps		%%zmm26, 0xfc0(%%rax)		\n\t"\
		"\n\t"\
	/***************************************************************************************/\
	/* Now do eight more radix-4 transforms, including the internal and external twiddles: */\
	/***************************************************************************************/\
		"\n\t"\
	/*
	Using upper block(s) of main array for temp-storage in section below led to a nasty AVX bug to track down:
	In fermat-mod mode, 4 block addresses in ascending order are add0,1,2,3 with no 'gaps' between blocks; for
	mersenne-mod addresses in asc. order are add0,2,3,1 with gaps between contiguous-data-block pairs 0,2 & 3,1.
	Thus for ferm-mod need [add2] as base-address of 'high-half' block for temp-storage; for mers-mod need [add3].
	In both cases have (add2 < add3) so instead use (add2 - add1) to differentiate: > 0 for fermat, < 0 for mersenne.
	*/\
		"movq	%[__add2],%%rsi		\n\t"/* destroyable copy of add2 */\
		"movq	%[__add2],%%rbx		\n\t"\
		"subq	%[__add1],%%rsi		\n\t"/* rsi = (add2 - add1); if this yields a borrow (i.e. sets CF) it's mersenne, else fermat. */\
		"cmovcq %[__add3],%%rbx	\n\t" /* if CF set (i.e. h > l), copy source [add3] into dest (rbx), else leave dest = [add2]. */\
		"\n\t"\
		"movq	%[__isrt2],%%rsi				\n\t	vmovaps		 (%%rdi),%%zmm31	\n\t"/* 2.0*/\
													/* base-addr [rax] and sincos-ptr [r10] in rcol both offset +0x200 vs lcol */\
	/********************************************/	/********************************************/\
	/* Block 2: t02,t12,t22,t32 -> r10,14,12,16 */	/* Block 6: t0A,t1A,t2A,t3A -> r18,1C,1A,1E */\
	/********************************************/	/********************************************/\
		"movq	%[__r00],%%rax					\n\t"\
		"movq	%%rsi,%%rcx					\n\t"\
		"movq	%%rsi,%%rdx					\n\t"\
		"movq	%[__c01],%%r10					\n\t"\
		"addq	$0x040,%%rcx	/* cc0 */		\n\t"\
		"addq	$0x0c0,%%rdx	/* cc1 */		\n\t"\
		"vmovaps		0x480(%%rax),%%zmm4		\n\t		vmovaps		 0x680(%%rax),%%zmm12		\n\t"\
		"vmovaps		0x580(%%rax),%%zmm0		\n\t		vmovaps		 0x780(%%rax),%%zmm8 		\n\t"\
		"vmovaps		0x4c0(%%rax),%%zmm5		\n\t		vmovaps		 0x6c0(%%rax),%%zmm15		\n\t"\
		"vmovaps		0x5c0(%%rax),%%zmm3		\n\t		vmovaps		 0x7c0(%%rax),%%zmm11		\n\t"\
	/* Strategic register-usage for sincos consts here cuts loads from 16 to 4: */\
		"vmovaps			 (%%rdx),%%zmm9 	\n\t		vmovaps		0x040(%%rdx),%%zmm10		\n\t"\
		"vmovaps		0x080(%%rdx),%%zmm1 	\n\t		vmovaps		0x0c0(%%rdx),%%zmm13		\n\t"\
		"vmulpd			%%zmm10,%%zmm5,%%zmm7	\n\t		vmulpd			 %%zmm1 ,%%zmm15,%%zmm15\n\t"\
		"vmulpd			%%zmm13,%%zmm3,%%zmm3	\n\t		vmulpd			 %%zmm10,%%zmm11,%%zmm11\n\t"\
		"vmulpd			%%zmm10,%%zmm4,%%zmm6	\n\t		vmulpd			 %%zmm1 ,%%zmm12,%%zmm14\n\t"\
		"vmulpd			%%zmm13,%%zmm0,%%zmm2	\n\t		vmulpd		0x780(%%rax),%%zmm10,%%zmm10\n\t"\
	"vfmadd132pd		%%zmm9 ,%%zmm7,%%zmm4	\n\t	vfmadd132pd	 		 %%zmm13,%%zmm15,%%zmm12\n\t"\
	"vfmadd132pd		%%zmm1 ,%%zmm3,%%zmm0	\n\t	vfmsub132pd	 		 %%zmm9 ,%%zmm11,%%zmm8 \n\t"\
	"vfmsub132pd		%%zmm9 ,%%zmm6,%%zmm5	\n\t	vfmsub132pd	 	0x6c0(%%rax),%%zmm14,%%zmm13\n\t"\
	"vfmsub132pd   0x5c0(%%rax),%%zmm2,%%zmm1	\n\t	vfmadd132pd	 	0x7c0(%%rax),%%zmm10,%%zmm9 \n\t"\
		"vaddpd			%%zmm0,%%zmm4,%%zmm6	\n\t		vsubpd		%%zmm8 ,%%zmm12,%%zmm14		\n\t"\
		"vaddpd			%%zmm1,%%zmm5,%%zmm7	\n\t		vsubpd		%%zmm9 ,%%zmm13,%%zmm15		\n\t"\
		"vsubpd			%%zmm0,%%zmm4,%%zmm4	\n\t		vaddpd		%%zmm8 ,%%zmm12,%%zmm12		\n\t"\
		"vsubpd			%%zmm1,%%zmm5,%%zmm5	\n\t		vaddpd		%%zmm9 ,%%zmm13,%%zmm13		\n\t"\
		"vmovaps		 0x500(%%rax),%%zmm2	\n\t		vmovaps		0x700(%%rax),%%zmm10		\n\t"\
		"vmovaps		 0x540(%%rax),%%zmm3	\n\t		vmovaps		0x740(%%rax),%%zmm11		\n\t"\
		"vmovaps		  (%%rcx),%%zmm16		\n\t		vmovaps	 0x040(%%rcx),%%zmm17	\n\t"/* twiddle */\
		"vmulpd			%%zmm17,%%zmm2,%%zmm0	\n\t		vmulpd		%%zmm16,%%zmm10,%%zmm8	 	\n\t"\
		"vmulpd			%%zmm17,%%zmm3,%%zmm1	\n\t		vmulpd		%%zmm16,%%zmm11,%%zmm9 		\n\t"\
	"vfmsub132pd		%%zmm16,%%zmm0,%%zmm3	\n\t	vfmadd132pd		%%zmm17,%%zmm8 ,%%zmm11		\n\t"\
	"vfmadd132pd		%%zmm16,%%zmm1,%%zmm2	\n\t	vfmsub132pd		%%zmm17,%%zmm9 ,%%zmm10		\n\t"\
		"vmovaps		 0x400(%%rax),%%zmm0	\n\t		vmovaps		0x600(%%rax),%%zmm8 		\n\t"\
		"vmovaps		 0x440(%%rax),%%zmm1	\n\t		vmovaps		0x640(%%rax),%%zmm9 		\n\t"\
		"vsubpd			%%zmm2,%%zmm0,%%zmm0	\n\t		vsubpd		%%zmm10,%%zmm8 ,%%zmm8 		\n\t"\
		"vsubpd			%%zmm3,%%zmm1,%%zmm1	\n\t		vsubpd		%%zmm11,%%zmm9 ,%%zmm9 		\n\t"\
	"vfmadd132pd		%%zmm31,%%zmm0,%%zmm2	\n\t	vfmadd132pd		%%zmm31,%%zmm8 ,%%zmm10		\n\t"\
	"vfmadd132pd		%%zmm31,%%zmm1,%%zmm3	\n\t	vfmadd132pd		%%zmm31,%%zmm9 ,%%zmm11		\n\t"\
		"vsubpd			%%zmm6,%%zmm2,%%zmm2	\n\t		vsubpd		%%zmm14,%%zmm8 ,%%zmm8 		\n\t"\
		"vsubpd			%%zmm7,%%zmm3,%%zmm3	\n\t		vsubpd		%%zmm15,%%zmm9 ,%%zmm9 		\n\t"\
	"vfmadd132pd		%%zmm31,%%zmm2,%%zmm6	\n\t	vfmadd132pd		%%zmm31,%%zmm8 ,%%zmm14		\n\t"\
	"vfmadd132pd		%%zmm31,%%zmm3,%%zmm7	\n\t	vfmadd132pd		%%zmm31,%%zmm9 ,%%zmm15		\n\t"\
/*** Use regs m16-23 here to avoid 16-reg spill/fill crud: ***/\
		"vmovaps		  (%%r10),%%zmm20		\n\t		vmovaps		0x200(%%r10),%%zmm22		\n\t"\
		"vmovaps	 0x040(%%r10),%%zmm21		\n\t		vmovaps		0x240(%%r10),%%zmm23		\n\t"\
		"vmulpd			%%zmm20,%%zmm6,%%zmm16	\n\t		vmulpd		%%zmm22,%%zmm14,%%zmm18		\n\t"\
		"vmulpd			%%zmm20,%%zmm7,%%zmm17	\n\t		vmulpd		%%zmm22,%%zmm15,%%zmm19		\n\t"\
	" vfmadd231pd		%%zmm21,%%zmm7,%%zmm16	\n\t	 vfmadd231pd	%%zmm23,%%zmm15,%%zmm18		\n\t"\
	"vfnmadd231pd		%%zmm21,%%zmm6,%%zmm17	\n\t	vfnmadd231pd	%%zmm23,%%zmm14,%%zmm19		\n\t"\
		"vmovaps	 0x080(%%r10),%%zmm20		\n\t		vmovaps		0x280(%%r10),%%zmm22		\n\t"\
		"vmovaps	 0x0c0(%%r10),%%zmm21		\n\t		vmovaps		0x2c0(%%r10),%%zmm23		\n\t"\
		"vmovaps		%%zmm16, 0x400(%%rax)	\n\t		vmovaps		%%zmm18, 0x600(%%rax)		\n\t"\
		"vmovaps		%%zmm17, 0x440(%%rax)	\n\t		vmovaps		%%zmm19, 0x640(%%rax)		\n\t"\
		"vmulpd			%%zmm20,%%zmm2,%%zmm16	\n\t		vmulpd		%%zmm22,%%zmm8 ,%%zmm18		\n\t"\
		"vmulpd			%%zmm20,%%zmm3,%%zmm17	\n\t		vmulpd		%%zmm22,%%zmm9 ,%%zmm19		\n\t"\
		"vsubpd			%%zmm5,%%zmm0,%%zmm0	\n\t		vsubpd		%%zmm13,%%zmm10,%%zmm10		\n\t"\
		"vsubpd			%%zmm4,%%zmm1,%%zmm1	\n\t		vsubpd		%%zmm12,%%zmm11,%%zmm11		\n\t"\
	" vfmadd231pd		%%zmm21,%%zmm3,%%zmm16	\n\t	 vfmadd231pd	%%zmm23,%%zmm9 ,%%zmm18		\n\t"\
	"vfnmadd231pd		%%zmm21,%%zmm2,%%zmm17	\n\t	vfnmadd231pd	%%zmm23,%%zmm8 ,%%zmm19		\n\t"\
	"vfmadd132pd		%%zmm31,%%zmm0,%%zmm5	\n\t	vfmadd132pd		%%zmm31,%%zmm10,%%zmm13		\n\t"\
	"vfmadd132pd		%%zmm31,%%zmm1,%%zmm4	\n\t	vfmadd132pd		%%zmm31,%%zmm11,%%zmm12		\n\t"\
		"vmovaps		%%zmm16,0x480(%%rax)	\n\t		vmovaps		%%zmm18, 0x680(%%rax)		\n\t"\
		"vmovaps		%%zmm17,0x4c0(%%rax)	\n\t		vmovaps		%%zmm19, 0x6c0(%%rax)		\n\t"\
		"addq		$0x100,%%r10				\n\t"\
		"vmovaps		  (%%r10),%%zmm16		\n\t		vmovaps		0x200(%%r10),%%zmm20		\n\t"\
		"vmovaps	 0x080(%%r10),%%zmm17		\n\t		vmovaps		0x280(%%r10),%%zmm21		\n\t"\
		"vmovaps	 0x040(%%r10),%%zmm18		\n\t		vmovaps		0x240(%%r10),%%zmm22		\n\t"\
		"vmovaps	 0x0c0(%%r10),%%zmm19		\n\t		vmovaps		0x2c0(%%r10),%%zmm23		\n\t"\
		"vmulpd			%%zmm16,%%zmm5,%%zmm2	\n\t		vmulpd		%%zmm20,%%zmm13,%%zmm8 		\n\t"\
		"vmulpd			%%zmm17,%%zmm0,%%zmm6	\n\t		vmulpd		%%zmm21,%%zmm10,%%zmm14		\n\t"\
		"vmulpd			%%zmm16,%%zmm1,%%zmm3	\n\t		vmulpd		%%zmm20,%%zmm11,%%zmm9 		\n\t"\
		"vmulpd			%%zmm17,%%zmm4,%%zmm7	\n\t		vmulpd		%%zmm21,%%zmm12,%%zmm15		\n\t"\
	" vfmadd231pd		%%zmm18,%%zmm1,%%zmm2	\n\t	 vfmadd231pd	%%zmm22,%%zmm11,%%zmm8 		\n\t"\
	" vfmadd231pd		%%zmm19,%%zmm4,%%zmm6	\n\t	 vfmadd231pd	%%zmm23,%%zmm12,%%zmm14		\n\t"\
	"vfnmadd231pd		%%zmm18,%%zmm5,%%zmm3	\n\t	vfnmadd231pd	%%zmm22,%%zmm13,%%zmm9 		\n\t"\
	"vfnmadd231pd		%%zmm19,%%zmm0,%%zmm7	\n\t	vfnmadd231pd	%%zmm23,%%zmm10,%%zmm15		\n\t"\
		"vmovaps		%%zmm2, 0x500(%%rax)	\n\t		vmovaps		%%zmm8 , 0x700(%%rax)		\n\t"\
		"vmovaps		%%zmm3, 0x540(%%rax)	\n\t		vmovaps		%%zmm9 , 0x740(%%rax)		\n\t"\
		"vmovaps		%%zmm6, 0x580(%%rax)	\n\t		vmovaps		%%zmm14, 0x780(%%rax)		\n\t"\
		"vmovaps		%%zmm7, 0x5c0(%%rax)	\n\t		vmovaps		%%zmm15, 0x7c0(%%rax)		\n\t"\
		"\n\t"\
	/********************************************/	/********************************************/\
	/* Block 4: t06,t16,t26,t36 -> r30,34,32,36 */	/* Block 8: t0E,t1E,t2E,t3E -> r38,3C,3A,3E */\
	/********************************************/	/********************************************/\
		"													addq		$0x200,%%r10				\n\t"\
		"vmovaps		0xc80(%%rax),%%zmm4		\n\t		vmovaps		0xe80(%%rax),%%zmm12		\n\t"\
		"vmovaps		0xd80(%%rax),%%zmm0		\n\t		vmovaps		0xf80(%%rax),%%zmm10		\n\t"\
		"vmovaps		0xcc0(%%rax),%%zmm7		\n\t		vmovaps		0xec0(%%rax),%%zmm15		\n\t"\
		"vmovaps		0xdc0(%%rax),%%zmm1		\n\t		vmovaps		0xfc0(%%rax),%%zmm9 		\n\t"\
	/* Strategic register-usage for sincos consts here cuts loads from 16 to 4: */\
		"vmovaps			 (%%rdx),%%zmm14	\n\t		vmovaps		0x040(%%rdx),%%zmm13		\n\t"\
		"vmovaps		0x080(%%rdx),%%zmm5 	\n\t		vmovaps		0x0c0(%%rdx),%%zmm8 		\n\t"\
	/* Cyclic-permute 4 paired-MUL lines so as to put a '...(%%rax)' in the SRC3 slot of last line'srcol: */\
		"vmulpd		 	 %%zmm14,%%zmm0,%%zmm2	\n\t		vmulpd			 %%zmm5 ,%%zmm10,%%zmm10\n\t"\
		"vmulpd		 	 %%zmm8 ,%%zmm7,%%zmm7	\n\t		vmulpd			 %%zmm14,%%zmm15,%%zmm15\n\t"\
		"vmulpd		 	 %%zmm14,%%zmm1,%%zmm3	\n\t		vmulpd			 %%zmm5 ,%%zmm9 ,%%zmm11\n\t"\
		"vmulpd		 	 %%zmm8 ,%%zmm4,%%zmm6	\n\t		vmulpd		0xe80(%%rax),%%zmm14,%%zmm14\n\t"\
	/* Also Cyclic-permute these 4 lines so middle (SRC2) op's reg-index order matches that of abovefour: */\
	"vfmadd132pd	 	 %%zmm13,%%zmm2,%%zmm1	\n\t	vfmsub132pd	 		 %%zmm8 ,%%zmm10,%%zmm9 \n\t"\
	"vfmadd132pd	 	 %%zmm5 ,%%zmm7,%%zmm4	\n\t	vfmadd132pd	 		 %%zmm13,%%zmm15,%%zmm12\n\t"\
	"vfmsub132pd	 	 %%zmm13,%%zmm3,%%zmm0	\n\t	vfmadd132pd	 	0xf80(%%rax),%%zmm11,%%zmm8 \n\t"\
	"vfmsub132pd	0xcc0(%%rax),%%zmm6,%%zmm5	\n\t	vfmsub132pd	 	0xec0(%%rax),%%zmm14,%%zmm13\n\t"\
		"vsubpd			%%zmm0,%%zmm4,%%zmm6	\n\t		vsubpd		%%zmm8 ,%%zmm12,%%zmm14		\n\t"\
		"vsubpd			%%zmm1,%%zmm5,%%zmm7	\n\t		vsubpd		%%zmm9 ,%%zmm13,%%zmm15		\n\t"\
		"vaddpd			%%zmm0,%%zmm4,%%zmm4	\n\t		vaddpd		%%zmm8 ,%%zmm12,%%zmm12		\n\t"\
		"vaddpd			%%zmm1,%%zmm5,%%zmm5	\n\t		vaddpd		%%zmm9 ,%%zmm13,%%zmm13		\n\t"\
		"vmovaps		 0xd00(%%rax),%%zmm2	\n\t		vmovaps		0xf00(%%rax),%%zmm10		\n\t"\
		"vmovaps		 0xd40(%%rax),%%zmm3	\n\t		vmovaps		0xf40(%%rax),%%zmm11		\n\t"\
		"vmovaps		  (%%rcx),%%zmm16		\n\t		vmovaps	 0x040(%%rcx),%%zmm17	\n\t"/* twiddle */\
		"vmulpd			%%zmm16,%%zmm2,%%zmm0	\n\t		vmulpd		%%zmm17,%%zmm10,%%zmm8	 	\n\t"\
		"vmulpd			%%zmm16,%%zmm3,%%zmm1	\n\t		vmulpd		%%zmm17,%%zmm11,%%zmm9 		\n\t"\
	"vfmsub132pd		%%zmm17,%%zmm0,%%zmm3	\n\t	vfmadd132pd		%%zmm16,%%zmm8 ,%%zmm11		\n\t"\
	"vfmadd132pd		%%zmm17,%%zmm1,%%zmm2	\n\t	vfmsub132pd		%%zmm16,%%zmm9 ,%%zmm10		\n\t"\
		"vmovaps		 0xc00(%%rax),%%zmm0	\n\t		vmovaps		0xe00(%%rax),%%zmm8 		\n\t"\
		"vmovaps		 0xc40(%%rax),%%zmm1	\n\t		vmovaps		0xe40(%%rax),%%zmm9 		\n\t"\
		"vsubpd			%%zmm2,%%zmm0,%%zmm0	\n\t		vsubpd		%%zmm10,%%zmm8 ,%%zmm8 		\n\t"\
		"vsubpd			%%zmm3,%%zmm1,%%zmm1	\n\t		vsubpd		%%zmm11,%%zmm9 ,%%zmm9 		\n\t"\
	"vfmadd132pd		%%zmm31,%%zmm0,%%zmm2	\n\t	vfmadd132pd		%%zmm31,%%zmm8 ,%%zmm10		\n\t"\
	"vfmadd132pd		%%zmm31,%%zmm1,%%zmm3	\n\t	vfmadd132pd		%%zmm31,%%zmm9 ,%%zmm11		\n\t"\
		"vsubpd			%%zmm6,%%zmm2,%%zmm2	\n\t		vsubpd		%%zmm14,%%zmm8 ,%%zmm8 		\n\t"\
		"vsubpd			%%zmm7,%%zmm3,%%zmm3	\n\t		vsubpd		%%zmm15,%%zmm9 ,%%zmm9 		\n\t"\
	"vfmadd132pd		%%zmm31,%%zmm2,%%zmm6	\n\t	vfmadd132pd		%%zmm31,%%zmm8 ,%%zmm14		\n\t"\
	"vfmadd132pd		%%zmm31,%%zmm3,%%zmm7	\n\t	vfmadd132pd		%%zmm31,%%zmm9 ,%%zmm15		\n\t"\
		"addq		$0x100,%%r10				\n\t"\
		"vmovaps		  (%%r10),%%zmm20		\n\t		vmovaps		0x200(%%r10),%%zmm22		\n\t"\
		"vmovaps	 0x040(%%r10),%%zmm21		\n\t		vmovaps		0x240(%%r10),%%zmm23		\n\t"\
		"vmulpd			%%zmm20,%%zmm6,%%zmm16	\n\t		vmulpd		%%zmm22,%%zmm14,%%zmm18		\n\t"\
		"vmulpd			%%zmm20,%%zmm7,%%zmm17	\n\t		vmulpd		%%zmm22,%%zmm15,%%zmm19		\n\t"\
	" vfmadd231pd		%%zmm21,%%zmm7,%%zmm16	\n\t	 vfmadd231pd	%%zmm23,%%zmm15,%%zmm18		\n\t"\
	"vfnmadd231pd		%%zmm21,%%zmm6,%%zmm17	\n\t	vfnmadd231pd	%%zmm23,%%zmm14,%%zmm19		\n\t"\
		"vmovaps	 0x080(%%r10),%%zmm20		\n\t		vmovaps		0x280(%%r10),%%zmm22		\n\t"\
		"vmovaps	 0x0c0(%%r10),%%zmm21		\n\t		vmovaps		0x2c0(%%r10),%%zmm23		\n\t"\
		"vmulpd			%%zmm20,%%zmm2,%%zmm6	\n\t		vmulpd		%%zmm22,%%zmm8 ,%%zmm14		\n\t"\
		"vmulpd			%%zmm20,%%zmm3,%%zmm7	\n\t		vmulpd		%%zmm22,%%zmm9 ,%%zmm15		\n\t"\
		"vsubpd			%%zmm5,%%zmm0,%%zmm0	\n\t		vsubpd		%%zmm13,%%zmm10,%%zmm10		\n\t"\
		"vsubpd			%%zmm4,%%zmm1,%%zmm1	\n\t		vsubpd		%%zmm12,%%zmm11,%%zmm11		\n\t"\
	" vfmadd231pd		%%zmm21,%%zmm3,%%zmm6	\n\t	 vfmadd231pd	%%zmm23,%%zmm9 ,%%zmm14		\n\t"\
	"vfnmadd231pd		%%zmm21,%%zmm2,%%zmm7	\n\t	vfnmadd231pd	%%zmm23,%%zmm8 ,%%zmm15		\n\t"\
	"vfmadd132pd		%%zmm31,%%zmm0,%%zmm5	\n\t	vfmadd132pd		%%zmm31,%%zmm10,%%zmm13		\n\t"\
	"vfmadd132pd		%%zmm31,%%zmm1,%%zmm4	\n\t	vfmadd132pd		%%zmm31,%%zmm11,%%zmm12		\n\t"\
		"vmovaps		%%zmm16, 0xc00(%%rax)	\n\t		vmovaps		%%zmm18, 0xe00(%%rax)		\n\t"\
		"vmovaps		%%zmm17, 0xc40(%%rax)	\n\t		vmovaps		%%zmm19, 0xe40(%%rax)		\n\t"\
		"vmovaps		%%zmm6 , 0xc80(%%rax)	\n\t		vmovaps		%%zmm14, 0xe80(%%rax)		\n\t"\
		"vmovaps		%%zmm7 , 0xcc0(%%rax)	\n\t		vmovaps		%%zmm15, 0xec0(%%rax)		\n\t"\
		"addq		$0x100,%%r10				\n\t"\
		"vmovaps		  (%%r10),%%zmm16		\n\t		vmovaps		0x200(%%r10),%%zmm20		\n\t"\
		"vmovaps	 0x080(%%r10),%%zmm17		\n\t		vmovaps		0x280(%%r10),%%zmm21		\n\t"\
		"vmovaps	 0x040(%%r10),%%zmm18		\n\t		vmovaps		0x240(%%r10),%%zmm22		\n\t"\
		"vmovaps	 0x0c0(%%r10),%%zmm19		\n\t		vmovaps		0x2c0(%%r10),%%zmm23		\n\t"\
		"vmulpd			%%zmm16,%%zmm5,%%zmm2	\n\t		vmulpd		%%zmm20,%%zmm13,%%zmm8 		\n\t"\
		"vmulpd			%%zmm16,%%zmm1,%%zmm3	\n\t		vmulpd		%%zmm20,%%zmm11,%%zmm9 		\n\t"\
		"vmulpd			%%zmm17,%%zmm0,%%zmm6	\n\t		vmulpd		%%zmm21,%%zmm10,%%zmm14		\n\t"\
		"vmulpd			%%zmm17,%%zmm4,%%zmm7	\n\t		vmulpd		%%zmm21,%%zmm12,%%zmm15		\n\t"\
	" vfmadd231pd		%%zmm18,%%zmm1,%%zmm2	\n\t	 vfmadd231pd	%%zmm22,%%zmm11,%%zmm8 		\n\t"\
	"vfnmadd231pd		%%zmm18,%%zmm5,%%zmm3	\n\t	vfnmadd231pd	%%zmm22,%%zmm13,%%zmm9 		\n\t"\
	" vfmadd231pd		%%zmm19,%%zmm4,%%zmm6	\n\t	 vfmadd231pd	%%zmm23,%%zmm12,%%zmm14		\n\t"\
	"vfnmadd231pd		%%zmm19,%%zmm0,%%zmm7	\n\t	vfnmadd231pd	%%zmm23,%%zmm10,%%zmm15		\n\t"\
		"vmovaps		%%zmm2, 0xd00(%%rax)	\n\t		vmovaps		%%zmm8 , 0xf00(%%rax)		\n\t"\
		"vmovaps		%%zmm3, 0xd40(%%rax)	\n\t		vmovaps		%%zmm9 , 0xf40(%%rax)		\n\t"\
		"vmovaps		%%zmm6, 0xd80(%%rax)	\n\t		vmovaps		%%zmm14, 0xf80(%%rax)		\n\t"\
		"vmovaps		%%zmm7, 0xdc0(%%rax)	\n\t		vmovaps		%%zmm15, 0xfc0(%%rax)		\n\t"\
		"\n\t"\
	/********************************************/	/********************************************/\
	/* Block 1: t00,t10,t20,t30 -> r00,04,02,06 */	/* Block 5: t08,t18,t28,t38 -> r08,0C,0A,0E */\
	/********************************************/	/********************************************/\
		"													vmovaps	(%%rsi),%%zmm10	\n\t"/* isrt2 */\
		"vmovaps			  (%%rax),%%zmm0	\n\t		vmovaps		 0x280(%%rax),%%zmm12		\n\t"\
		"vmovaps		 0x040(%%rax),%%zmm1	\n\t		vmovaps		 0x2c0(%%rax),%%zmm13		\n\t"\
		"vmovaps		 0x100(%%rax),%%zmm2	\n\t		vmovaps		 0x380(%%rax),%%zmm8 		\n\t"\
		"vmovaps		 0x140(%%rax),%%zmm3	\n\t		vmovaps		 0x3c0(%%rax),%%zmm9 		\n\t"\
		"vsubpd			%%zmm2 ,%%zmm0,%%zmm0	\n\t		vsubpd		%%zmm12,%%zmm13,%%zmm13		\n\t"\
		"vsubpd			%%zmm3 ,%%zmm1,%%zmm1	\n\t	vfmadd132pd		%%zmm31,%%zmm13,%%zmm12		\n\t"\
	"vfmadd132pd		%%zmm31,%%zmm0,%%zmm2	\n\t		vsubpd		%%zmm9 ,%%zmm8 ,%%zmm8 		\n\t"\
	"vfmadd132pd		%%zmm31,%%zmm1,%%zmm3	\n\t	vfmadd132pd		%%zmm31,%%zmm8 ,%%zmm9 		\n\t"\
		"vmovaps		 0x080(%%rax),%%zmm4	\n\t		vmulpd		%%zmm10,%%zmm8 ,%%zmm8 		\n\t"\
		"vmovaps		 0x0c0(%%rax),%%zmm5	\n\t		vmulpd		%%zmm10,%%zmm9 ,%%zmm9 		\n\t"\
		"vmovaps		 0x180(%%rax),%%zmm6	\n\t		vmovaps		%%zmm12,%%zmm14				\n\t"\
		"vmovaps		 0x1c0(%%rax),%%zmm7	\n\t		vmovaps		%%zmm13,%%zmm15				\n\t"\
		"vsubpd			%%zmm6 ,%%zmm4,%%zmm4	\n\t	vfmsub132pd		%%zmm10,%%zmm8 ,%%zmm12		\n\t"\
		"vsubpd			%%zmm7 ,%%zmm5,%%zmm5	\n\t	vfmsub132pd		%%zmm10,%%zmm9 ,%%zmm13		\n\t"\
	"vfmadd132pd		%%zmm31,%%zmm4,%%zmm6	\n\t	vfmadd132pd		%%zmm10,%%zmm8 ,%%zmm14		\n\t"\
	"vfmadd132pd		%%zmm31,%%zmm5,%%zmm7	\n\t	vfmadd132pd		%%zmm10,%%zmm9 ,%%zmm15		\n\t"\
		"movq		%[__c10],%%rcx			\n\t"/* base-twiddle in l/rcol = c00/c04 ==> rcx+0x200 in rcol */\
		"vaddpd			%%zmm6,%%zmm2,%%zmm2	\n\t		vmovaps		 0x200(%%rax),%%zmm8 		\n\t"\
		"vaddpd			%%zmm7,%%zmm3,%%zmm3	\n\t		vmovaps		 0x240(%%rax),%%zmm9 		\n\t"\
		"vmovaps		%%zmm2,      (%%rax)	\n\t		vmovaps		 0x300(%%rax),%%zmm10		\n\t"\
		"vmovaps		%%zmm3, 0x040(%%rax)	\n\t		vmovaps		 0x340(%%rax),%%zmm11		\n\t"\
		"													vsubpd		%%zmm11,%%zmm8 ,%%zmm8 		\n\t"\
		"													vsubpd		%%zmm10,%%zmm9 ,%%zmm9 		\n\t"\
	"vfnmadd231pd		%%zmm31,%%zmm6,%%zmm2	\n\t	vfmadd132pd		%%zmm31,%%zmm8 ,%%zmm11		\n\t"\
	"vfnmadd231pd		%%zmm31,%%zmm7,%%zmm3	\n\t	vfmadd132pd		%%zmm31,%%zmm9 ,%%zmm10		\n\t"\
		"vmovaps		%%zmm2,%%zmm6			\n\t		vsubpd		%%zmm12,%%zmm11,%%zmm11		\n\t"\
		"vmovaps		%%zmm3,%%zmm7			\n\t		vsubpd		%%zmm13,%%zmm9 ,%%zmm9 		\n\t"\
		"vmovaps		  (%%rcx),%%zmm16		\n\t		vmovaps		0x180(%%rcx),%%zmm20		\n\t"\
		"vmovaps	 0x040(%%rcx),%%zmm17		\n\t		vmovaps		0x1c0(%%rcx),%%zmm21		\n\t"\
		"													vmovaps		0x200(%%rcx),%%zmm22		\n\t"\
		"													vmovaps		0x240(%%rcx),%%zmm23		\n\t"\
		"subq $0x80,%%rcx		\n\t"/* put c00 in rcx to ease bookkeeping */\
/*c10:*/"vmulpd			%%zmm16,%%zmm2,%%zmm2	\n\t	vfmadd132pd		%%zmm31,%%zmm11,%%zmm12		\n\t"\
		"vmulpd			%%zmm16,%%zmm3,%%zmm3	\n\t	vfmadd132pd		%%zmm31,%%zmm9 ,%%zmm13		\n\t"\
	" vfmadd231pd		%%zmm17,%%zmm7,%%zmm2	\n\t		vmulpd		%%zmm20,%%zmm12,%%zmm16		\n\t"/* c04:*/\
	"vfnmadd231pd		%%zmm17,%%zmm6,%%zmm3	\n\t		vmulpd		%%zmm20,%%zmm13,%%zmm17		\n\t"\
		"												 vfmadd231pd	%%zmm21,%%zmm13,%%zmm16		\n\t"\
		"												vfnmadd231pd	%%zmm21,%%zmm12,%%zmm17		\n\t"\
		"													vmulpd	 %%zmm22,%%zmm11,%%zmm12		\n\t"/* c14 */\
		"													vmulpd	 %%zmm22,%%zmm9 ,%%zmm13		\n\t"\
		"												 vfmadd231pd %%zmm23,%%zmm9 ,%%zmm12		\n\t"\
		"												vfnmadd231pd %%zmm23,%%zmm11,%%zmm13		\n\t"\
		"vaddpd			%%zmm5,%%zmm0,%%zmm0	\n\t		vmovaps	%%zmm12,0x280(%%rax)			\n\t"/* r0a,0b */\
		"vsubpd			%%zmm4,%%zmm1,%%zmm1	\n\t		vmovaps	%%zmm13,0x2c0(%%rax)			\n\t"\
		"vmovaps	%%zmm2,0x080(%%rax)	/* r02,03 */\n\t	vmovaps	%%zmm16,0x200(%%rax)			\n\t"/* r08,09 */\
		"vmovaps	%%zmm3,0x0c0(%%rax)				\n\t	vmovaps	%%zmm17,0x240(%%rax)			\n\t"\
		"vmovaps	 0x100(%%rcx),%%zmm16		\n\t		vmovaps		0x300(%%rcx),%%zmm20		\n\t"\
		"vmovaps	 0x140(%%rcx),%%zmm17		\n\t		vmovaps		0x340(%%rcx),%%zmm21		\n\t"\
		"vmovaps	 0x180(%%rcx),%%zmm18		\n\t		vmovaps		0x380(%%rcx),%%zmm22		\n\t"\
		"vmovaps	 0x1c0(%%rcx),%%zmm19		\n\t		vmovaps		0x3c0(%%rcx),%%zmm23		\n\t"\
	"vfnmadd132pd	%%zmm31,%%zmm0,%%zmm5		\n\t		vsubpd			%%zmm15,%%zmm8 ,%%zmm8 	\n\t"\
	" vfmadd132pd	%%zmm31,%%zmm1,%%zmm4		\n\t		vsubpd			%%zmm14,%%zmm10,%%zmm10	\n\t"\
		"vmulpd		%%zmm16,%%zmm0,%%zmm2 /*c08*/\n\t	vfmadd132pd		%%zmm31,%%zmm8 ,%%zmm15		\n\t"\
		"vmulpd		%%zmm16,%%zmm1,%%zmm3		\n\t	vfmadd132pd		%%zmm31,%%zmm10,%%zmm14		\n\t"\
	" vfmadd231pd	%%zmm17,%%zmm1,%%zmm2		\n\t		vmulpd		%%zmm20,%%zmm15,%%zmm12		\n\t"/* c0C */\
	"vfnmadd231pd	%%zmm17,%%zmm0,%%zmm3		\n\t		vmulpd		%%zmm20,%%zmm10,%%zmm13		\n\t"\
		"												 vfmadd231pd	%%zmm21,%%zmm10,%%zmm12		\n\t"\
		"												vfnmadd231pd	%%zmm21,%%zmm15,%%zmm13		\n\t"\
		"vmulpd		%%zmm18,%%zmm5,%%zmm0	/* c18*/\n\t	vmulpd		%%zmm22,%%zmm8 ,%%zmm15		\n\t"/* c1C */\
		"vmulpd		%%zmm18,%%zmm4,%%zmm1		\n\t		vmulpd		%%zmm22,%%zmm14,%%zmm10		\n\t"\
	" vfmadd231pd	%%zmm19,%%zmm4,%%zmm0		\n\t	 vfmadd231pd	%%zmm23,%%zmm14,%%zmm15		\n\t"\
	"vfnmadd231pd	%%zmm19,%%zmm5,%%zmm1		\n\t	vfnmadd231pd	%%zmm23,%%zmm8 ,%%zmm10		\n\t"\
		"vmovaps	%%zmm2,0x100(%%rax)	/* r04,05 */\n\t	vmovaps	%%zmm12,0x300(%%rax)			\n\t"/* r0c,0d */\
		"vmovaps	%%zmm3,0x140(%%rax)				\n\t	vmovaps	%%zmm13,0x340(%%rax)			\n\t"\
		"vmovaps	%%zmm0,0x180(%%rax)	/* r06,07 */\n\t	vmovaps	%%zmm15,0x380(%%rax)			\n\t"/* r0e,0f */\
		"vmovaps	%%zmm1,0x1c0(%%rax)			\n\t		vmovaps	%%zmm10,0x3c0(%%rax)			\n\t"\
		"\n\t"\
	/********************************************/	/********************************************/\
	/* Block 3: t04,t14,t24,t34 -> r20,24,22,26 */	/* Block 7: t0C,t1C,t2C,t3C -> r28,2C,2A,2E */\
	/********************************************/	/********************************************/\
		"leaq	0x040(%%rsi),%%rcx	\n\t"/* cc0; Note cc0/ss0 are shared between lcol/rcol, so no rcx-offset until get to twiddles*/\
		"vmovaps		0x880(%%rax),%%zmm4		\n\t		vmovaps		 0xa80(%%rax),%%zmm12		\n\t"\
		"vmovaps		0x980(%%rax),%%zmm0		\n\t		vmovaps		 0xb80(%%rax),%%zmm8 		\n\t"\
		"vmovaps		0x8c0(%%rax),%%zmm5		\n\t		vmovaps		 0xac0(%%rax),%%zmm13		\n\t"\
		"vmovaps		0x9c0(%%rax),%%zmm3		\n\t		vmovaps		 0xbc0(%%rax),%%zmm11		\n\t"\
		"vmovaps			%%zmm4 ,%%zmm6		\n\t		vmovaps		 	%%zmm12,%%zmm14			\n\t"\
		"vmovaps			%%zmm0 ,%%zmm2		\n\t		vmovaps		 	%%zmm8 ,%%zmm10			\n\t"\
		"vmovaps			%%zmm5 ,%%zmm7		\n\t		vmovaps		 	%%zmm13,%%zmm15			\n\t"\
	/*	"vmovaps			%%zmm1 ,%%zmm3		\n\t		vmovaps		 	%%zmm9 ,%%zmm11			\n\t"*/\
	/* Strategic register-usage for sincos consts here cuts loads from 16 to 4: */\
		"vmovaps			 (%%rcx),%%zmm9 	\n\t		vmovaps		0x040(%%rcx),%%zmm1 		\n\t"\
		"vmulpd		 	 %%zmm1 ,%%zmm7,%%zmm7	\n\t		vmulpd			 %%zmm9 ,%%zmm15,%%zmm15\n\t"\
		"vmulpd		 	 %%zmm9 ,%%zmm3,%%zmm3	\n\t		vmulpd			 %%zmm1 ,%%zmm11,%%zmm11\n\t"\
		"vmulpd		 	 %%zmm1 ,%%zmm6,%%zmm6	\n\t		vmulpd			 %%zmm9 ,%%zmm14,%%zmm14\n\t"\
		"vmulpd		 	 %%zmm9 ,%%zmm2,%%zmm2	\n\t		vmulpd			 %%zmm1 ,%%zmm10,%%zmm10\n\t"\
	"vfmadd132pd	 	 %%zmm9 ,%%zmm7,%%zmm4	\n\t	vfmadd132pd			 %%zmm1 ,%%zmm15,%%zmm12\n\t"\
	"vfmadd132pd	 	 %%zmm1 ,%%zmm3,%%zmm0	\n\t	vfmadd132pd			 %%zmm9 ,%%zmm11,%%zmm8 \n\t"\
	"vfmsub132pd	 	 %%zmm9 ,%%zmm6,%%zmm5	\n\t	vfmsub132pd			 %%zmm1 ,%%zmm14,%%zmm13\n\t"\
	"vfmsub132pd	0x9c0(%%rax),%%zmm2,%%zmm1	\n\t	vfmsub132pd		0xbc0(%%rax),%%zmm10,%%zmm9 \n\t"\
		"vsubpd			%%zmm0,%%zmm4,%%zmm6	\n\t		vsubpd			%%zmm8 ,%%zmm12,%%zmm14	\n\t"\
		"vsubpd			%%zmm1,%%zmm5,%%zmm7	\n\t		vsubpd			%%zmm9 ,%%zmm13,%%zmm15	\n\t"\
		"vaddpd			%%zmm0,%%zmm4,%%zmm4	\n\t		vaddpd			%%zmm8 ,%%zmm12,%%zmm12	\n\t"\
		"vaddpd			%%zmm1,%%zmm5,%%zmm5	\n\t		vaddpd			%%zmm9 ,%%zmm13,%%zmm13	\n\t"\
		"vmovaps		 0x900(%%rax),%%zmm2	\n\t		vmovaps		0xb00(%%rax),%%zmm10		\n\t"\
		"vmovaps		 0x940(%%rax),%%zmm3	\n\t		vmovaps		0xb40(%%rax),%%zmm11		\n\t"\
		"vmovaps		 0x800(%%rax),%%zmm0	\n\t		vmovaps		0xa00(%%rax),%%zmm8 		\n\t"\
		"vmovaps		 0x840(%%rax),%%zmm1	\n\t		vmovaps		0xa40(%%rax),%%zmm9 		\n\t"\
		"vaddpd		 0x940(%%rax),%%zmm2,%%zmm2	\n\t		vsubpd		0xb40(%%rax),%%zmm10,%%zmm10\n\t"\
		"vsubpd		 0x900(%%rax),%%zmm3,%%zmm3	\n\t		vaddpd		0xb00(%%rax),%%zmm11,%%zmm11\n\t"\
		"vmulpd			  (%%rsi),%%zmm2,%%zmm2	\n\t		vmulpd			 (%%rsi),%%zmm10,%%zmm10\n\t"\
		"vmulpd			  (%%rsi),%%zmm3,%%zmm3	\n\t		vmulpd			 (%%rsi),%%zmm11,%%zmm11\n\t"\
		"vsubpd			%%zmm2 ,%%zmm0,%%zmm0	\n\t		vsubpd			%%zmm10,%%zmm8 ,%%zmm8 	\n\t"\
		"vsubpd			%%zmm3 ,%%zmm1,%%zmm1	\n\t		vsubpd			%%zmm11,%%zmm9 ,%%zmm9 	\n\t"\
	"vfmadd132pd		%%zmm31,%%zmm0,%%zmm2	\n\t	vfmadd132pd		 %%zmm31,%%zmm8 ,%%zmm10	\n\t"\
	"vfmadd132pd		%%zmm31,%%zmm1,%%zmm3	\n\t	vfmadd132pd		 %%zmm31,%%zmm9 ,%%zmm11	\n\t"\
		"movq		%[__c02],%%rcx				\n\t"/* base-twiddle addr in rcol = c06, so rcx offset +0x200 vs lcol */\
		"vmovaps		  (%%rcx),%%zmm16		\n\t		vmovaps		0x200(%%rcx),%%zmm20		\n\t"\
		"vmovaps	 0x040(%%rcx),%%zmm17		\n\t		vmovaps		0x240(%%rcx),%%zmm21		\n\t"\
		"vsubpd			%%zmm4,%%zmm2,%%zmm2	\n\t		vsubpd			%%zmm14,%%zmm8 ,%%zmm8 	\n\t"\
		"vsubpd			%%zmm5,%%zmm3,%%zmm3	\n\t		vsubpd			%%zmm15,%%zmm9 ,%%zmm9 	\n\t"\
	"vfmadd132pd		%%zmm31,%%zmm2,%%zmm4	\n\t	vfmadd132pd		%%zmm31,%%zmm8 ,%%zmm14	\n\t"\
	"vfmadd132pd		%%zmm31,%%zmm3,%%zmm5	\n\t	vfmadd132pd		%%zmm31,%%zmm9 ,%%zmm15	\n\t"\
		"vmovaps		%%zmm2, 0x880(%%rax)	\n\t		vmovaps		%%zmm8 , 0xa80(%%rax)		\n\t"\
		"vmovaps		%%zmm3, 0x8c0(%%rax)	\n\t		vmovaps		%%zmm9 , 0xac0(%%rax)		\n\t"\
		"vmulpd			%%zmm16,%%zmm4,%%zmm2	\n\t		vmulpd			%%zmm20,%%zmm14,%%zmm8	\n\t"\
		"vmulpd			%%zmm16,%%zmm5,%%zmm3	\n\t		vmulpd			%%zmm20,%%zmm15,%%zmm9	\n\t"\
	" vfmadd231pd		%%zmm17,%%zmm5,%%zmm2	\n\t	 vfmadd231pd		%%zmm21,%%zmm15,%%zmm8	\n\t"\
	"vfnmadd231pd		%%zmm17,%%zmm4,%%zmm3	\n\t	vfnmadd231pd		%%zmm21,%%zmm14,%%zmm9	\n\t"\
		"movq		%[__c12],%%rcx				\n\t"		/* rcol uses c16 */\
		"vmovaps		  (%%rcx),%%zmm16		\n\t		vmovaps		0x200(%%rcx),%%zmm20		\n\t"\
		"vmovaps	 0x040(%%rcx),%%zmm17		\n\t		vmovaps		0x240(%%rcx),%%zmm21		\n\t"\
		"vmovaps		 0x880(%%rax),%%zmm4	\n\t		vmovaps		 0xa80(%%rax),%%zmm14		\n\t"\
		"vmovaps		 0x8c0(%%rax),%%zmm5	\n\t		vmovaps		 0xac0(%%rax),%%zmm15		\n\t"\
		"vmovaps	%%zmm2,0x800(%%rax)	/* r20,21 */\n\t	vmovaps	%%zmm8 ,0xa00(%%rax)			\n\t"/* r28,29 */\
		"vmovaps	%%zmm3,0x840(%%rax)			\n\t		vmovaps	%%zmm9 ,0xa40(%%rax)			\n\t"\
		"vmulpd			%%zmm16,%%zmm4,%%zmm2	\n\t		vmulpd			%%zmm20,%%zmm14,%%zmm8	\n\t"\
		"vmulpd			%%zmm16,%%zmm5,%%zmm3	\n\t		vmulpd			%%zmm20,%%zmm15,%%zmm9	\n\t"\
	" vfmadd231pd		%%zmm17,%%zmm5,%%zmm2	\n\t	 vfmadd231pd		%%zmm21,%%zmm15,%%zmm8	\n\t"\
	"vfnmadd231pd		%%zmm17,%%zmm4,%%zmm3	\n\t	vfnmadd231pd		%%zmm21,%%zmm14,%%zmm9	\n\t"\
		"movq		%[__c0A],%%rcx				\n\t"		/* rcol uses c0E */\
		"vmovaps		  (%%rcx),%%zmm16		\n\t		vmovaps		0x200(%%rcx),%%zmm20		\n\t"\
		"vmovaps	 0x040(%%rcx),%%zmm17		\n\t		vmovaps		0x240(%%rcx),%%zmm21		\n\t"\
		"vsubpd			%%zmm7,%%zmm0,%%zmm0	\n\t		vsubpd			%%zmm13,%%zmm10,%%zmm10	\n\t"\
		"vsubpd			%%zmm6,%%zmm1,%%zmm1	\n\t		vsubpd			%%zmm12,%%zmm11,%%zmm11	\n\t"\
	"vfmadd132pd		%%zmm31,%%zmm0,%%zmm7	\n\t	vfmadd132pd		%%zmm31,%%zmm10,%%zmm13		\n\t"\
	"vfmadd132pd		%%zmm31,%%zmm1,%%zmm6	\n\t	vfmadd132pd		%%zmm31,%%zmm11,%%zmm12		\n\t"\
		"vmovaps	%%zmm2,0x880(%%rax)	/* r22,23 */\n\t	vmovaps	%%zmm8 ,0xa80(%%rax)			\n\t"/* r2a,2b */\
		"vmovaps	%%zmm3,0x8c0(%%rax)			\n\t		vmovaps	%%zmm9 ,0xac0(%%rax)			\n\t"\
		"vmulpd			%%zmm16,%%zmm7,%%zmm4	\n\t		vmulpd			%%zmm20,%%zmm13,%%zmm8	\n\t"\
		"vmulpd			%%zmm16,%%zmm1,%%zmm5	\n\t		vmulpd			%%zmm20,%%zmm11,%%zmm9	\n\t"\
	" vfmadd231pd		%%zmm17,%%zmm1,%%zmm4	\n\t	 vfmadd231pd		%%zmm21,%%zmm11,%%zmm8	\n\t"\
	"vfnmadd231pd		%%zmm17,%%zmm7,%%zmm5	\n\t	vfnmadd231pd		%%zmm21,%%zmm13,%%zmm9	\n\t"\
		"movq		%[__c1A],%%rcx				\n\t"		/* rcol uses c1E */\
		"vmovaps		  (%%rcx),%%zmm16		\n\t		vmovaps		0x200(%%rcx),%%zmm20		\n\t"\
		"vmovaps	 0x040(%%rcx),%%zmm17		\n\t		vmovaps		0x240(%%rcx),%%zmm21		\n\t"\
		"vmovaps	%%zmm4,0x900(%%rax)	/* r24,25 */\n\t	vmovaps	%%zmm8 ,0xb00(%%rax)			\n\t"/* r2c,2d */\
		"vmovaps	%%zmm5,0x940(%%rax)			\n\t		vmovaps	%%zmm9 ,0xb40(%%rax)			\n\t"\
		"vmulpd			%%zmm16,%%zmm0,%%zmm4	\n\t		vmulpd			%%zmm20,%%zmm10,%%zmm8	\n\t"\
		"vmulpd			%%zmm16,%%zmm6,%%zmm5	\n\t		vmulpd			%%zmm20,%%zmm12,%%zmm9	\n\t"\
	" vfmadd231pd		%%zmm17,%%zmm6,%%zmm4	\n\t	 vfmadd231pd		%%zmm21,%%zmm12,%%zmm8	\n\t"\
	"vfnmadd231pd		%%zmm17,%%zmm0,%%zmm5	\n\t	vfnmadd231pd		%%zmm21,%%zmm10,%%zmm9	\n\t"\
		"vmovaps	%%zmm4,0x980(%%rax)	/* r26,27 */\n\t	vmovaps	%%zmm8 ,0xb80(%%rax)			\n\t"/* r2e,2f */\
		"vmovaps	%%zmm5,0x9c0(%%rax)			\n\t		vmovaps	%%zmm9 ,0xbc0(%%rax)			\n\t"\
	/*******************************************
	***** Finish with 8-way 'un'terleaving: ****
	Using the AVX-512 data layout, the rcol pattern is:
		a[ 0- 7] = re[ 0, 8,16,24, 4,12,20,28].d0	a[ 8-15] = im[ 0, 8,16,24, 4,12,20,28].d0
		a[16-23] = re[ 2,10,18,26, 6,14,22,30].d0	a[24-31] = im[ 2,10,18,26, 6,14,22,30].d0
		a[32-39] = re[ 1, 9,17,25, 5,13,21,29].d0	a[40-47] = im[ 1, 9,17,25, 5,13,21,29].d0
		a[48-55] = re[ 3,11,19,27, 7,15,23,31].d0	a[56-63] = im[ 3,11,19,27, 7,15,23,31].d0 ,
	and remaining seven 64-double blocks repeat same pattern with elts d1-d7 of the vector-doubles.
	*******************************************/\
		"movl	$0b00110011,%%eax	\n\t"/* Constant for vblendmpd instructions goes into mask-reg k1 */\
		"kmovw	%%eax,%%k1			\n\t"\
		/* Init vector index-consts needed by vpermt2pd instructions - if regs were at a premium,
		could also init just prior to [3] and use zmm6,7 to hold index-consts: */\
		"movq	$0x0c040e0608000a02,%%rax	\n\t"/* zmm30 = 8+4 0+4 8+6 0+6 8+0 0+0 8+2 0+2 [msw at left] */\
		"movq	$0x0d050f0709010b03,%%rbx	\n\t"/* zmm31 = 8+5 0+5 8+7 0+7 8+1 0+1 8+3 0+3 */\
			"vmovq		%%rax,%%xmm0 		\n\t"\
			"vmovq		%%rbx,%%xmm1 		\n\t"\
			"vpmovzxbq	%%xmm0,%%zmm30		\n\t"\
			"vpmovzxbq	%%xmm1,%%zmm31		\n\t"\
		"movq	%[__add0],%%rax	\n\t"\
		"movq	%[__add1],%%rbx	\n\t"\
		"movq	%[__add2],%%rcx	\n\t"\
		"movq	%[__add3],%%rdx	\n\t"\
		"movq	%[__add4],%%r10	\n\t"\
		"movq	%[__add5],%%r11	\n\t"\
		"movq	%[__add6],%%r12	\n\t"\
		"movq	%[__add7],%%r13	\n\t"\
	"movq	%[__r00] ,%%rsi	\n\t"\
	/**** See above USE_16_REG version of this macro for closing transposition using columnwise scatter-stores. ****/\
	/**** a[ 0- 7] = re[ 0, 8,16,24, 4,12,20,28].d0	a[ 8-15] = im[ 0, 8,16,24, 4,12,20,28].d0 : ****/\
		/* Read transposed row-data from the local store: */\
		"vmovaps 		 (%%rsi),%%zmm0			\n\t	vmovaps 0x040(%%rsi),%%zmm12	\n\t"\
		"vmovaps 	0x400(%%rsi),%%zmm1			\n\t	vmovaps 0x440(%%rsi),%%zmm13	\n\t"\
		"vmovaps 	0x800(%%rsi),%%zmm2			\n\t	vmovaps 0x840(%%rsi),%%zmm14	\n\t"\
		"vmovaps 	0xc00(%%rsi),%%zmm3			\n\t	vmovaps 0xc40(%%rsi),%%zmm15	\n\t"\
		"vmovaps 	0x200(%%rsi),%%zmm4			\n\t	vmovaps 0x240(%%rsi),%%zmm16	\n\t"\
		"vmovaps 	0x600(%%rsi),%%zmm5			\n\t	vmovaps 0x640(%%rsi),%%zmm17	\n\t"\
		"vmovaps 	0xa00(%%rsi),%%zmm6			\n\t	vmovaps 0xa40(%%rsi),%%zmm18	\n\t"\
		"vmovaps 	0xe00(%%rsi),%%zmm7			\n\t	vmovaps 0xe40(%%rsi),%%zmm19	\n\t"\
		/* [1] Shuffle the 4-aparts - note the different patterning of the first and second output quartet: */\
		"vshuff64x2	$0b01000100, %%zmm4,	%%zmm0,	%%zmm8 	\n\t	vshuff64x2	$0b01000100, %%zmm16,%%zmm12,	%%zmm20	\n\t"\
		"vshuff64x2	$0b11101110, %%zmm4,	%%zmm0,	%%zmm4 	\n\t	vshuff64x2	$0b11101110, %%zmm16,%%zmm12,	%%zmm16	\n\t"\
		"vshuff64x2	$0b01000100, %%zmm5,	%%zmm1,	%%zmm9	\n\t	vshuff64x2	$0b01000100, %%zmm17,%%zmm13,	%%zmm21	\n\t"\
		"vshuff64x2	$0b11101110, %%zmm5,	%%zmm1,	%%zmm5	\n\t	vshuff64x2	$0b11101110, %%zmm17,%%zmm13,	%%zmm17	\n\t"\
		"vshuff64x2	$0b00010001, %%zmm6,	%%zmm2,	%%zmm10	\n\t	vshuff64x2	$0b00010001, %%zmm18,%%zmm14,	%%zmm22	\n\t"\
		"vshuff64x2	$0b10111011, %%zmm6,	%%zmm2,	%%zmm6	\n\t	vshuff64x2	$0b10111011, %%zmm18,%%zmm14,	%%zmm18	\n\t"\
		"vshuff64x2	$0b00010001, %%zmm7,	%%zmm3,	%%zmm11	\n\t	vshuff64x2	$0b00010001, %%zmm19,%%zmm15,	%%zmm23	\n\t"\
		"vshuff64x2	$0b10111011, %%zmm7,	%%zmm3,	%%zmm7	\n\t	vshuff64x2	$0b10111011, %%zmm19,%%zmm15,	%%zmm19	\n\t"\
		/* [2] Blend in the 2-aparts */\
		"vblendmpd	%%zmm8 ,	%%zmm10,	%%zmm0%{%%k1%}	\n\t	vblendmpd	%%zmm20,	%%zmm22,	%%zmm12%{%%k1%}	\n\t"\
		"vblendmpd	%%zmm10,	%%zmm8 ,	%%zmm8%{%%k1%}	\n\t	vblendmpd	%%zmm22,	%%zmm20,	%%zmm20%{%%k1%}	\n\t"\
		"vblendmpd	%%zmm4 ,	%%zmm6 ,	%%zmm1%{%%k1%}	\n\t	vblendmpd	%%zmm16,	%%zmm18,	%%zmm13%{%%k1%}	\n\t"\
		"vblendmpd	%%zmm6 ,	%%zmm4 ,	%%zmm4%{%%k1%}	\n\t	vblendmpd	%%zmm18,	%%zmm16,	%%zmm16%{%%k1%}	\n\t"\
		"vblendmpd	%%zmm9 ,	%%zmm11,	%%zmm2%{%%k1%}	\n\t	vblendmpd	%%zmm21,	%%zmm23,	%%zmm14%{%%k1%}	\n\t"\
		"vblendmpd	%%zmm11,	%%zmm9 ,	%%zmm9%{%%k1%}	\n\t	vblendmpd	%%zmm23,	%%zmm21,	%%zmm21%{%%k1%}	\n\t"\
		"vblendmpd	%%zmm5 ,	%%zmm7 ,	%%zmm3%{%%k1%}	\n\t	vblendmpd	%%zmm17,	%%zmm19,	%%zmm15%{%%k1%}	\n\t"\
		"vblendmpd	%%zmm7 ,	%%zmm5 ,	%%zmm5%{%%k1%}	\n\t	vblendmpd	%%zmm19,	%%zmm17,	%%zmm17%{%%k1%}	\n\t"\
		/* [3] Shuffle or permute in the 1-aparts */\
		"vshufpd	$0b00000000,%%zmm2,		%%zmm0,%%zmm10 	\n\t	vshufpd	$0b00000000,%%zmm14,	%%zmm12,%%zmm22	\n\t"\
		"vshufpd	$0b11111111,%%zmm2,		%%zmm0,%%zmm11 	\n\t	vshufpd	$0b11111111,%%zmm14,	%%zmm12,%%zmm23	\n\t"\
		"vmovapd	%%zmm8,%%zmm2							\n\t	vmovapd	%%zmm20,%%zmm14	\n\t"\
		"vpermt2pd				%%zmm9,		%%zmm30,%%zmm2 	\n\t	vpermt2pd			%%zmm21,	%%zmm30,%%zmm14	\n\t"\
		"vpermt2pd				%%zmm9,		%%zmm31,%%zmm8	\n\t	vpermt2pd			%%zmm21,	%%zmm31,%%zmm20	\n\t"\
		"vshufpd	$0b00000000,%%zmm3,		%%zmm1,%%zmm0 	\n\t	vshufpd	$0b00000000,%%zmm15,	%%zmm13,%%zmm12	\n\t"\
		"vshufpd	$0b11111111,%%zmm3,		%%zmm1,%%zmm1 	\n\t	vshufpd	$0b11111111,%%zmm15,	%%zmm13,%%zmm13	\n\t"\
		"vmovapd	%%zmm4,%%zmm3							\n\t	vmovapd	%%zmm16,%%zmm15	\n\t"\
		"vpermt2pd				%%zmm5,		%%zmm30,%%zmm3 	\n\t	vpermt2pd			%%zmm17,	%%zmm30,%%zmm15	\n\t"\
		"vpermt2pd				%%zmm5,		%%zmm31,%%zmm4	\n\t	vpermt2pd			%%zmm17,	%%zmm31,%%zmm16	\n\t"\
		/* Write resulting 8 vec-cmplx data to eight 512[0x200]-byte-separated slots of main data array: */\
		"vmovaps		%%zmm10,(%%rax)		\n\t	vmovaps	%%zmm22,0x40(%%rax)	\n\t"\
		"vmovaps		%%zmm11,(%%rbx)		\n\t	vmovaps	%%zmm23,0x40(%%rbx)	\n\t"\
		"vmovaps		%%zmm2 ,(%%rcx)		\n\t	vmovaps	%%zmm14,0x40(%%rcx)	\n\t"\
		"vmovaps		%%zmm8 ,(%%rdx)		\n\t	vmovaps	%%zmm20,0x40(%%rdx)	\n\t"\
		"vmovaps		%%zmm0 ,(%%r10)		\n\t	vmovaps	%%zmm12,0x40(%%r10)	\n\t"\
		"vmovaps		%%zmm1 ,(%%r11)		\n\t	vmovaps	%%zmm13,0x40(%%r11)	\n\t"\
		"vmovaps		%%zmm3 ,(%%r12)		\n\t	vmovaps	%%zmm15,0x40(%%r12)	\n\t"\
		"vmovaps		%%zmm4 ,(%%r13)		\n\t	vmovaps	%%zmm16,0x40(%%r13)	\n\t"\
	/**** a[16-23] = re[ 2,10,18,26, 6,14,22,30].d0	a[24-31] = im[ 2,10,18,26, 6,14,22,30].d0 : ****/\
	"addq	$0x80 ,%%rax	\n\t"\
	"addq	$0x80 ,%%rbx	\n\t"\
	"addq	$0x80 ,%%rcx	\n\t"\
	"addq	$0x80 ,%%rdx	\n\t"\
	"addq	$0x80 ,%%r10	\n\t"\
	"addq	$0x80 ,%%r11	\n\t"\
	"addq	$0x80 ,%%r12	\n\t"\
	"addq	$0x80 ,%%r13	\n\t"\
	"addq	$0x100,%%rsi	\n\t"\
		/* Read transposed row-data from the local store: */\
		"vmovaps 		 (%%rsi),%%zmm0			\n\t	vmovaps 0x040(%%rsi),%%zmm12	\n\t"\
		"vmovaps 	0x400(%%rsi),%%zmm1			\n\t	vmovaps 0x440(%%rsi),%%zmm13	\n\t"\
		"vmovaps 	0x800(%%rsi),%%zmm2			\n\t	vmovaps 0x840(%%rsi),%%zmm14	\n\t"\
		"vmovaps 	0xc00(%%rsi),%%zmm3			\n\t	vmovaps 0xc40(%%rsi),%%zmm15	\n\t"\
		"vmovaps 	0x200(%%rsi),%%zmm4			\n\t	vmovaps 0x240(%%rsi),%%zmm16	\n\t"\
		"vmovaps 	0x600(%%rsi),%%zmm5			\n\t	vmovaps 0x640(%%rsi),%%zmm17	\n\t"\
		"vmovaps 	0xa00(%%rsi),%%zmm6			\n\t	vmovaps 0xa40(%%rsi),%%zmm18	\n\t"\
		"vmovaps 	0xe00(%%rsi),%%zmm7			\n\t	vmovaps 0xe40(%%rsi),%%zmm19	\n\t"\
		/* [1] Shuffle the 4-aparts - note the different patterning of the first and second output quartet: */\
		"vshuff64x2	$0b01000100, %%zmm4,	%%zmm0,	%%zmm8 	\n\t	vshuff64x2	$0b01000100, %%zmm16,%%zmm12,	%%zmm20	\n\t"\
		"vshuff64x2	$0b11101110, %%zmm4,	%%zmm0,	%%zmm4 	\n\t	vshuff64x2	$0b11101110, %%zmm16,%%zmm12,	%%zmm16	\n\t"\
		"vshuff64x2	$0b01000100, %%zmm5,	%%zmm1,	%%zmm9	\n\t	vshuff64x2	$0b01000100, %%zmm17,%%zmm13,	%%zmm21	\n\t"\
		"vshuff64x2	$0b11101110, %%zmm5,	%%zmm1,	%%zmm5	\n\t	vshuff64x2	$0b11101110, %%zmm17,%%zmm13,	%%zmm17	\n\t"\
		"vshuff64x2	$0b00010001, %%zmm6,	%%zmm2,	%%zmm10	\n\t	vshuff64x2	$0b00010001, %%zmm18,%%zmm14,	%%zmm22	\n\t"\
		"vshuff64x2	$0b10111011, %%zmm6,	%%zmm2,	%%zmm6	\n\t	vshuff64x2	$0b10111011, %%zmm18,%%zmm14,	%%zmm18	\n\t"\
		"vshuff64x2	$0b00010001, %%zmm7,	%%zmm3,	%%zmm11	\n\t	vshuff64x2	$0b00010001, %%zmm19,%%zmm15,	%%zmm23	\n\t"\
		"vshuff64x2	$0b10111011, %%zmm7,	%%zmm3,	%%zmm7	\n\t	vshuff64x2	$0b10111011, %%zmm19,%%zmm15,	%%zmm19	\n\t"\
		/* [2] Blend in the 2-aparts */\
		"vblendmpd	%%zmm8 ,	%%zmm10,	%%zmm0%{%%k1%}	\n\t	vblendmpd	%%zmm20,	%%zmm22,	%%zmm12%{%%k1%}	\n\t"\
		"vblendmpd	%%zmm10,	%%zmm8 ,	%%zmm8%{%%k1%}	\n\t	vblendmpd	%%zmm22,	%%zmm20,	%%zmm20%{%%k1%}	\n\t"\
		"vblendmpd	%%zmm4 ,	%%zmm6 ,	%%zmm1%{%%k1%}	\n\t	vblendmpd	%%zmm16,	%%zmm18,	%%zmm13%{%%k1%}	\n\t"\
		"vblendmpd	%%zmm6 ,	%%zmm4 ,	%%zmm4%{%%k1%}	\n\t	vblendmpd	%%zmm18,	%%zmm16,	%%zmm16%{%%k1%}	\n\t"\
		"vblendmpd	%%zmm9 ,	%%zmm11,	%%zmm2%{%%k1%}	\n\t	vblendmpd	%%zmm21,	%%zmm23,	%%zmm14%{%%k1%}	\n\t"\
		"vblendmpd	%%zmm11,	%%zmm9 ,	%%zmm9%{%%k1%}	\n\t	vblendmpd	%%zmm23,	%%zmm21,	%%zmm21%{%%k1%}	\n\t"\
		"vblendmpd	%%zmm5 ,	%%zmm7 ,	%%zmm3%{%%k1%}	\n\t	vblendmpd	%%zmm17,	%%zmm19,	%%zmm15%{%%k1%}	\n\t"\
		"vblendmpd	%%zmm7 ,	%%zmm5 ,	%%zmm5%{%%k1%}	\n\t	vblendmpd	%%zmm19,	%%zmm17,	%%zmm17%{%%k1%}	\n\t"\
		/* [3] Shuffle or permute in the 1-aparts */\
		"vshufpd	$0b00000000,%%zmm2,		%%zmm0,%%zmm10 	\n\t	vshufpd	$0b00000000,%%zmm14,	%%zmm12,%%zmm22	\n\t"\
		"vshufpd	$0b11111111,%%zmm2,		%%zmm0,%%zmm11 	\n\t	vshufpd	$0b11111111,%%zmm14,	%%zmm12,%%zmm23	\n\t"\
		"vmovapd	%%zmm8,%%zmm2							\n\t	vmovapd	%%zmm20,%%zmm14	\n\t"\
		"vpermt2pd				%%zmm9,		%%zmm30,%%zmm2 	\n\t	vpermt2pd			%%zmm21,	%%zmm30,%%zmm14	\n\t"\
		"vpermt2pd				%%zmm9,		%%zmm31,%%zmm8	\n\t	vpermt2pd			%%zmm21,	%%zmm31,%%zmm20	\n\t"\
		"vshufpd	$0b00000000,%%zmm3,		%%zmm1,%%zmm0 	\n\t	vshufpd	$0b00000000,%%zmm15,	%%zmm13,%%zmm12	\n\t"\
		"vshufpd	$0b11111111,%%zmm3,		%%zmm1,%%zmm1 	\n\t	vshufpd	$0b11111111,%%zmm15,	%%zmm13,%%zmm13	\n\t"\
		"vmovapd	%%zmm4,%%zmm3							\n\t	vmovapd	%%zmm16,%%zmm15	\n\t"\
		"vpermt2pd				%%zmm5,		%%zmm30,%%zmm3 	\n\t	vpermt2pd			%%zmm17,	%%zmm30,%%zmm15	\n\t"\
		"vpermt2pd				%%zmm5,		%%zmm31,%%zmm4	\n\t	vpermt2pd			%%zmm17,	%%zmm31,%%zmm16	\n\t"\
		/* Write resulting 8 vec-cmplx data to eight 512[0x200]-byte-separated slots of main data array: */\
		"vmovaps		%%zmm10,(%%rax)		\n\t	vmovaps	%%zmm22,0x40(%%rax)	\n\t"\
		"vmovaps		%%zmm11,(%%rbx)		\n\t	vmovaps	%%zmm23,0x40(%%rbx)	\n\t"\
		"vmovaps		%%zmm2 ,(%%rcx)		\n\t	vmovaps	%%zmm14,0x40(%%rcx)	\n\t"\
		"vmovaps		%%zmm8 ,(%%rdx)		\n\t	vmovaps	%%zmm20,0x40(%%rdx)	\n\t"\
		"vmovaps		%%zmm0 ,(%%r10)		\n\t	vmovaps	%%zmm12,0x40(%%r10)	\n\t"\
		"vmovaps		%%zmm1 ,(%%r11)		\n\t	vmovaps	%%zmm13,0x40(%%r11)	\n\t"\
		"vmovaps		%%zmm3 ,(%%r12)		\n\t	vmovaps	%%zmm15,0x40(%%r12)	\n\t"\
		"vmovaps		%%zmm4 ,(%%r13)		\n\t	vmovaps	%%zmm16,0x40(%%r13)	\n\t"\
		"\n\t"\
	/**** a[32-39] = re[ 1, 9,17,25, 5,13,21,29].d0	a[40-47] = im[ 1, 9,17,25, 5,13,21,29].d0 : ****/\
	"addq	$0x80 ,%%rax	\n\t"\
	"addq	$0x80 ,%%rbx	\n\t"\
	"addq	$0x80 ,%%rcx	\n\t"\
	"addq	$0x80 ,%%rdx	\n\t"\
	"addq	$0x80 ,%%r10	\n\t"\
	"addq	$0x80 ,%%r11	\n\t"\
	"addq	$0x80 ,%%r12	\n\t"\
	"addq	$0x80 ,%%r13	\n\t"\
	"subq	$0x80,%%rsi	\n\t"\
		/* Read transposed row-data from the local store: */\
		"vmovaps 		 (%%rsi),%%zmm0			\n\t	vmovaps 0x040(%%rsi),%%zmm12	\n\t"\
		"vmovaps 	0x400(%%rsi),%%zmm1			\n\t	vmovaps 0x440(%%rsi),%%zmm13	\n\t"\
		"vmovaps 	0x800(%%rsi),%%zmm2			\n\t	vmovaps 0x840(%%rsi),%%zmm14	\n\t"\
		"vmovaps 	0xc00(%%rsi),%%zmm3			\n\t	vmovaps 0xc40(%%rsi),%%zmm15	\n\t"\
		"vmovaps 	0x200(%%rsi),%%zmm4			\n\t	vmovaps 0x240(%%rsi),%%zmm16	\n\t"\
		"vmovaps 	0x600(%%rsi),%%zmm5			\n\t	vmovaps 0x640(%%rsi),%%zmm17	\n\t"\
		"vmovaps 	0xa00(%%rsi),%%zmm6			\n\t	vmovaps 0xa40(%%rsi),%%zmm18	\n\t"\
		"vmovaps 	0xe00(%%rsi),%%zmm7			\n\t	vmovaps 0xe40(%%rsi),%%zmm19	\n\t"\
		/* [1] Shuffle the 4-aparts - note the different patterning of the first and second output quartet: */\
		"vshuff64x2	$0b01000100, %%zmm4,	%%zmm0,	%%zmm8 	\n\t	vshuff64x2	$0b01000100, %%zmm16,%%zmm12,	%%zmm20	\n\t"\
		"vshuff64x2	$0b11101110, %%zmm4,	%%zmm0,	%%zmm4 	\n\t	vshuff64x2	$0b11101110, %%zmm16,%%zmm12,	%%zmm16	\n\t"\
		"vshuff64x2	$0b01000100, %%zmm5,	%%zmm1,	%%zmm9	\n\t	vshuff64x2	$0b01000100, %%zmm17,%%zmm13,	%%zmm21	\n\t"\
		"vshuff64x2	$0b11101110, %%zmm5,	%%zmm1,	%%zmm5	\n\t	vshuff64x2	$0b11101110, %%zmm17,%%zmm13,	%%zmm17	\n\t"\
		"vshuff64x2	$0b00010001, %%zmm6,	%%zmm2,	%%zmm10	\n\t	vshuff64x2	$0b00010001, %%zmm18,%%zmm14,	%%zmm22	\n\t"\
		"vshuff64x2	$0b10111011, %%zmm6,	%%zmm2,	%%zmm6	\n\t	vshuff64x2	$0b10111011, %%zmm18,%%zmm14,	%%zmm18	\n\t"\
		"vshuff64x2	$0b00010001, %%zmm7,	%%zmm3,	%%zmm11	\n\t	vshuff64x2	$0b00010001, %%zmm19,%%zmm15,	%%zmm23	\n\t"\
		"vshuff64x2	$0b10111011, %%zmm7,	%%zmm3,	%%zmm7	\n\t	vshuff64x2	$0b10111011, %%zmm19,%%zmm15,	%%zmm19	\n\t"\
		/* [2] Blend in the 2-aparts */\
		"vblendmpd	%%zmm8 ,	%%zmm10,	%%zmm0%{%%k1%}	\n\t	vblendmpd	%%zmm20,	%%zmm22,	%%zmm12%{%%k1%}	\n\t"\
		"vblendmpd	%%zmm10,	%%zmm8 ,	%%zmm8%{%%k1%}	\n\t	vblendmpd	%%zmm22,	%%zmm20,	%%zmm20%{%%k1%}	\n\t"\
		"vblendmpd	%%zmm4 ,	%%zmm6 ,	%%zmm1%{%%k1%}	\n\t	vblendmpd	%%zmm16,	%%zmm18,	%%zmm13%{%%k1%}	\n\t"\
		"vblendmpd	%%zmm6 ,	%%zmm4 ,	%%zmm4%{%%k1%}	\n\t	vblendmpd	%%zmm18,	%%zmm16,	%%zmm16%{%%k1%}	\n\t"\
		"vblendmpd	%%zmm9 ,	%%zmm11,	%%zmm2%{%%k1%}	\n\t	vblendmpd	%%zmm21,	%%zmm23,	%%zmm14%{%%k1%}	\n\t"\
		"vblendmpd	%%zmm11,	%%zmm9 ,	%%zmm9%{%%k1%}	\n\t	vblendmpd	%%zmm23,	%%zmm21,	%%zmm21%{%%k1%}	\n\t"\
		"vblendmpd	%%zmm5 ,	%%zmm7 ,	%%zmm3%{%%k1%}	\n\t	vblendmpd	%%zmm17,	%%zmm19,	%%zmm15%{%%k1%}	\n\t"\
		"vblendmpd	%%zmm7 ,	%%zmm5 ,	%%zmm5%{%%k1%}	\n\t	vblendmpd	%%zmm19,	%%zmm17,	%%zmm17%{%%k1%}	\n\t"\
		/* [3] Shuffle or permute in the 1-aparts */\
		"vshufpd	$0b00000000,%%zmm2,		%%zmm0,%%zmm10 	\n\t	vshufpd	$0b00000000,%%zmm14,	%%zmm12,%%zmm22	\n\t"\
		"vshufpd	$0b11111111,%%zmm2,		%%zmm0,%%zmm11 	\n\t	vshufpd	$0b11111111,%%zmm14,	%%zmm12,%%zmm23	\n\t"\
		"vmovapd	%%zmm8,%%zmm2							\n\t	vmovapd	%%zmm20,%%zmm14	\n\t"\
		"vpermt2pd				%%zmm9,		%%zmm30,%%zmm2 	\n\t	vpermt2pd			%%zmm21,	%%zmm30,%%zmm14	\n\t"\
		"vpermt2pd				%%zmm9,		%%zmm31,%%zmm8	\n\t	vpermt2pd			%%zmm21,	%%zmm31,%%zmm20	\n\t"\
		"vshufpd	$0b00000000,%%zmm3,		%%zmm1,%%zmm0 	\n\t	vshufpd	$0b00000000,%%zmm15,	%%zmm13,%%zmm12	\n\t"\
		"vshufpd	$0b11111111,%%zmm3,		%%zmm1,%%zmm1 	\n\t	vshufpd	$0b11111111,%%zmm15,	%%zmm13,%%zmm13	\n\t"\
		"vmovapd	%%zmm4,%%zmm3							\n\t	vmovapd	%%zmm16,%%zmm15	\n\t"\
		"vpermt2pd				%%zmm5,		%%zmm30,%%zmm3 	\n\t	vpermt2pd			%%zmm17,	%%zmm30,%%zmm15	\n\t"\
		"vpermt2pd				%%zmm5,		%%zmm31,%%zmm4	\n\t	vpermt2pd			%%zmm17,	%%zmm31,%%zmm16	\n\t"\
		/* Write resulting 8 vec-cmplx data to eight 512[0x200]-byte-separated slots of main data array: */\
		"vmovaps		%%zmm10,(%%rax)		\n\t	vmovaps	%%zmm22,0x40(%%rax)	\n\t"\
		"vmovaps		%%zmm11,(%%rbx)		\n\t	vmovaps	%%zmm23,0x40(%%rbx)	\n\t"\
		"vmovaps		%%zmm2 ,(%%rcx)		\n\t	vmovaps	%%zmm14,0x40(%%rcx)	\n\t"\
		"vmovaps		%%zmm8 ,(%%rdx)		\n\t	vmovaps	%%zmm20,0x40(%%rdx)	\n\t"\
		"vmovaps		%%zmm0 ,(%%r10)		\n\t	vmovaps	%%zmm12,0x40(%%r10)	\n\t"\
		"vmovaps		%%zmm1 ,(%%r11)		\n\t	vmovaps	%%zmm13,0x40(%%r11)	\n\t"\
		"vmovaps		%%zmm3 ,(%%r12)		\n\t	vmovaps	%%zmm15,0x40(%%r12)	\n\t"\
		"vmovaps		%%zmm4 ,(%%r13)		\n\t	vmovaps	%%zmm16,0x40(%%r13)	\n\t"\
		"\n\t"\
	/**** a[48-55] = re[ 3,11,19,27, 7,15,23,31].d0	a[56-63] = im[ 3,11,19,27, 7,15,23,31].d0 : ****/\
	"addq	$0x80 ,%%rax	\n\t"\
	"addq	$0x80 ,%%rbx	\n\t"\
	"addq	$0x80 ,%%rcx	\n\t"\
	"addq	$0x80 ,%%rdx	\n\t"\
	"addq	$0x80 ,%%r10	\n\t"\
	"addq	$0x80 ,%%r11	\n\t"\
	"addq	$0x80 ,%%r12	\n\t"\
	"addq	$0x80 ,%%r13	\n\t"\
	"addq	$0x100,%%rsi	\n\t"\
		/* Read transposed row-data from the local store: */\
		"vmovaps 		 (%%rsi),%%zmm0			\n\t	vmovaps 0x040(%%rsi),%%zmm12	\n\t"\
		"vmovaps 	0x400(%%rsi),%%zmm1			\n\t	vmovaps 0x440(%%rsi),%%zmm13	\n\t"\
		"vmovaps 	0x800(%%rsi),%%zmm2			\n\t	vmovaps 0x840(%%rsi),%%zmm14	\n\t"\
		"vmovaps 	0xc00(%%rsi),%%zmm3			\n\t	vmovaps 0xc40(%%rsi),%%zmm15	\n\t"\
		"vmovaps 	0x200(%%rsi),%%zmm4			\n\t	vmovaps 0x240(%%rsi),%%zmm16	\n\t"\
		"vmovaps 	0x600(%%rsi),%%zmm5			\n\t	vmovaps 0x640(%%rsi),%%zmm17	\n\t"\
		"vmovaps 	0xa00(%%rsi),%%zmm6			\n\t	vmovaps 0xa40(%%rsi),%%zmm18	\n\t"\
		"vmovaps 	0xe00(%%rsi),%%zmm7			\n\t	vmovaps 0xe40(%%rsi),%%zmm19	\n\t"\
		/* [1] Shuffle the 4-aparts - note the different patterning of the first and second output quartet: */\
		"vshuff64x2	$0b01000100, %%zmm4,	%%zmm0,	%%zmm8 	\n\t	vshuff64x2	$0b01000100, %%zmm16,%%zmm12,	%%zmm20	\n\t"\
		"vshuff64x2	$0b11101110, %%zmm4,	%%zmm0,	%%zmm4 	\n\t	vshuff64x2	$0b11101110, %%zmm16,%%zmm12,	%%zmm16	\n\t"\
		"vshuff64x2	$0b01000100, %%zmm5,	%%zmm1,	%%zmm9	\n\t	vshuff64x2	$0b01000100, %%zmm17,%%zmm13,	%%zmm21	\n\t"\
		"vshuff64x2	$0b11101110, %%zmm5,	%%zmm1,	%%zmm5	\n\t	vshuff64x2	$0b11101110, %%zmm17,%%zmm13,	%%zmm17	\n\t"\
		"vshuff64x2	$0b00010001, %%zmm6,	%%zmm2,	%%zmm10	\n\t	vshuff64x2	$0b00010001, %%zmm18,%%zmm14,	%%zmm22	\n\t"\
		"vshuff64x2	$0b10111011, %%zmm6,	%%zmm2,	%%zmm6	\n\t	vshuff64x2	$0b10111011, %%zmm18,%%zmm14,	%%zmm18	\n\t"\
		"vshuff64x2	$0b00010001, %%zmm7,	%%zmm3,	%%zmm11	\n\t	vshuff64x2	$0b00010001, %%zmm19,%%zmm15,	%%zmm23	\n\t"\
		"vshuff64x2	$0b10111011, %%zmm7,	%%zmm3,	%%zmm7	\n\t	vshuff64x2	$0b10111011, %%zmm19,%%zmm15,	%%zmm19	\n\t"\
		/* [2] Blend in the 2-aparts */\
		"vblendmpd	%%zmm8 ,	%%zmm10,	%%zmm0%{%%k1%}	\n\t	vblendmpd	%%zmm20,	%%zmm22,	%%zmm12%{%%k1%}	\n\t"\
		"vblendmpd	%%zmm10,	%%zmm8 ,	%%zmm8%{%%k1%}	\n\t	vblendmpd	%%zmm22,	%%zmm20,	%%zmm20%{%%k1%}	\n\t"\
		"vblendmpd	%%zmm4 ,	%%zmm6 ,	%%zmm1%{%%k1%}	\n\t	vblendmpd	%%zmm16,	%%zmm18,	%%zmm13%{%%k1%}	\n\t"\
		"vblendmpd	%%zmm6 ,	%%zmm4 ,	%%zmm4%{%%k1%}	\n\t	vblendmpd	%%zmm18,	%%zmm16,	%%zmm16%{%%k1%}	\n\t"\
		"vblendmpd	%%zmm9 ,	%%zmm11,	%%zmm2%{%%k1%}	\n\t	vblendmpd	%%zmm21,	%%zmm23,	%%zmm14%{%%k1%}	\n\t"\
		"vblendmpd	%%zmm11,	%%zmm9 ,	%%zmm9%{%%k1%}	\n\t	vblendmpd	%%zmm23,	%%zmm21,	%%zmm21%{%%k1%}	\n\t"\
		"vblendmpd	%%zmm5 ,	%%zmm7 ,	%%zmm3%{%%k1%}	\n\t	vblendmpd	%%zmm17,	%%zmm19,	%%zmm15%{%%k1%}	\n\t"\
		"vblendmpd	%%zmm7 ,	%%zmm5 ,	%%zmm5%{%%k1%}	\n\t	vblendmpd	%%zmm19,	%%zmm17,	%%zmm17%{%%k1%}	\n\t"\
		/* [3] Shuffle or permute in the 1-aparts */\
		"vshufpd	$0b00000000,%%zmm2,		%%zmm0,%%zmm10 	\n\t	vshufpd	$0b00000000,%%zmm14,	%%zmm12,%%zmm22	\n\t"\
		"vshufpd	$0b11111111,%%zmm2,		%%zmm0,%%zmm11 	\n\t	vshufpd	$0b11111111,%%zmm14,	%%zmm12,%%zmm23	\n\t"\
		"vmovapd	%%zmm8,%%zmm2							\n\t	vmovapd	%%zmm20,%%zmm14	\n\t"\
		"vpermt2pd				%%zmm9,		%%zmm30,%%zmm2 	\n\t	vpermt2pd			%%zmm21,	%%zmm30,%%zmm14	\n\t"\
		"vpermt2pd				%%zmm9,		%%zmm31,%%zmm8	\n\t	vpermt2pd			%%zmm21,	%%zmm31,%%zmm20	\n\t"\
		"vshufpd	$0b00000000,%%zmm3,		%%zmm1,%%zmm0 	\n\t	vshufpd	$0b00000000,%%zmm15,	%%zmm13,%%zmm12	\n\t"\
		"vshufpd	$0b11111111,%%zmm3,		%%zmm1,%%zmm1 	\n\t	vshufpd	$0b11111111,%%zmm15,	%%zmm13,%%zmm13	\n\t"\
		"vmovapd	%%zmm4,%%zmm3							\n\t	vmovapd	%%zmm16,%%zmm15	\n\t"\
		"vpermt2pd				%%zmm5,		%%zmm30,%%zmm3 	\n\t	vpermt2pd			%%zmm17,	%%zmm30,%%zmm15	\n\t"\
		"vpermt2pd				%%zmm5,		%%zmm31,%%zmm4	\n\t	vpermt2pd			%%zmm17,	%%zmm31,%%zmm16	\n\t"\
		/* Write resulting 8 vec-cmplx data to eight 512[0x200]-byte-separated slots of main data array: */\
		"vmovaps		%%zmm10,(%%rax)		\n\t	vmovaps	%%zmm22,0x40(%%rax)	\n\t"\
		"vmovaps		%%zmm11,(%%rbx)		\n\t	vmovaps	%%zmm23,0x40(%%rbx)	\n\t"\
		"vmovaps		%%zmm2 ,(%%rcx)		\n\t	vmovaps	%%zmm14,0x40(%%rcx)	\n\t"\
		"vmovaps		%%zmm8 ,(%%rdx)		\n\t	vmovaps	%%zmm20,0x40(%%rdx)	\n\t"\
		"vmovaps		%%zmm0 ,(%%r10)		\n\t	vmovaps	%%zmm12,0x40(%%r10)	\n\t"\
		"vmovaps		%%zmm1 ,(%%r11)		\n\t	vmovaps	%%zmm13,0x40(%%r11)	\n\t"\
		"vmovaps		%%zmm3 ,(%%r12)		\n\t	vmovaps	%%zmm15,0x40(%%r12)	\n\t"\
		"vmovaps		%%zmm4 ,(%%r13)		\n\t	vmovaps	%%zmm16,0x40(%%r13)	\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__add1] "m" (Xadd1)\
		 ,[__add2] "m" (Xadd2)\
		 ,[__add3] "m" (Xadd3)\
		 ,[__add4] "m" (Xadd4)\
		 ,[__add5] "m" (Xadd5)\
		 ,[__add6] "m" (Xadd6)\
		 ,[__add7] "m" (Xadd7)\
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__r00] "m" (Xr00)\
		 ,[__r08] "m" (Xr08)\
		 ,[__r10] "m" (Xr10)\
		 ,[__r20] "m" (Xr20)\
		 ,[__r28] "m" (Xr28)\
		 ,[__r30] "m" (Xr30)\
		 ,[__c01] "m" (Xc01)\
		 ,[__c02] "m" (Xc02)\
		 ,[__c04] "m" (Xc04)\
		 ,[__c06] "m" (Xc06)\
		 ,[__c0A] "m" (Xc0A)\
		 ,[__c10] "m" (Xc10)\
		 ,[__c12] "m" (Xc12)\
		 ,[__c1A] "m" (Xc1A)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r10","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15","xmm16","xmm17","xmm18","xmm19","xmm20","xmm21","xmm22","xmm23","xmm24","xmm25","xmm26","xmm27","xmm28","xmm29","xmm30","xmm31"	/* Clobbered registers */\
	);\
	}

   #endif	// (IMCI512 or AVX512?) toggle

  #endif	// USE_16_REG ?

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
  #ifdef ALL_FMA	// Aggressive-FMA version: Replace most ADD/SUB by FMA-with-one-unity-multiplicand

	// Aggressive-FMA: replace [58 ADD, 158 SUB, 98 MUL, 232 FMA, 844 memref] ==> [10 ADD, 2 SUB, 428 FMA (232 nontrivial), 98 MUL, 844 memref].
	//
	#define SSE2_RADIX32_WRAPPER_DIF(Xadd0,Xadd1,Xadd2,Xadd3,Xr00,Xr10,Xr20,Xr30,Xisrt2,Xcc0,Xc00,Xc01,Xc02,Xc03,Xc05,Xc07)\
	{\
	__asm__ volatile (\
		"movq	%[__r00] ,%%rsi	\n\t"\
		"movq	%[__add0],%%rax\n\t"\
		"movq	%[__add1],%%rbx\n\t"\
		"movq	%[__add2],%%rcx\n\t"\
		"movq	%[__add3],%%rdx\n\t"\
	"/**** Start with 4-way interleaving: ****/\n\t"\
	"/* a[j+p0]: Inputs from [add0,add0+0x200,add1,add1+0x200]+0x0. Outputs into r0 +0/1, 16/17, 32/33, 48/49: */	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps		 (%%rax),%%ymm1						\n\t		vmovaps	0x020(%%rax),%%ymm3								\n\t"\
		"vmovaps		 (%%rbx),%%ymm5						\n\t		vmovaps	0x020(%%rbx),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps		 (%%rcx),%%ymm6						\n\t		vmovaps	0x020(%%rcx),%%ymm14							\n\t"\
		"vmovaps		 (%%rdx),%%ymm5						\n\t		vmovaps	0x020(%%rdx),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm6 ,%%ymm0				\n\t		vshufpd	$15 ,%%ymm13,%%ymm14 ,%%ymm2					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm6 ,%%ymm6				\n\t		vshufpd	$0  ,%%ymm13,%%ymm14 ,%%ymm14					\n\t"\
		"vperm2f128 $32 ,%%ymm0 ,%%ymm7 ,%%ymm5				\n\t		vperm2f128 $32 ,%%ymm2 ,%%ymm15,%%ymm13		/* outB	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm0 ,%%ymm7 ,%%ymm7				\n\t		vperm2f128 $49 ,%%ymm2 ,%%ymm15,%%ymm15		/* outD	*/	\n\t"\
		"vperm2f128 $32 ,%%ymm6 ,%%ymm1 ,%%ymm0				\n\t		vperm2f128 $32 ,%%ymm14,%%ymm3 ,%%ymm2 		/* outA	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm6 ,%%ymm1 ,%%ymm1				\n\t		vperm2f128 $49 ,%%ymm14,%%ymm3 ,%%ymm3 		/* outC	*/	\n\t"\
		"vmovaps 	%%ymm5 ,0x400(%%rsi)					\n\t		vmovaps %%ymm13,0x420(%%rsi)				/* outB	*/	\n\t"\
		"vmovaps 	%%ymm7 ,0x600(%%rsi)					\n\t		vmovaps %%ymm15,0x620(%%rsi)				/* outD	*/	\n\t"\
		"vmovaps 	%%ymm0 ,     (%%rsi)					\n\t		vmovaps %%ymm2 ,0x020(%%rsi)				/* outA	*/	\n\t"\
		"vmovaps 	%%ymm1 ,0x200(%%rsi)					\n\t		vmovaps %%ymm3 ,0x220(%%rsi)				/* outC	*/	\n\t"\
	"/* a[j+p2]: Inputs from [add0,add0+0x200,add1,add1+0x200]+0x40. Outputs into **r8** +0/1, 16/17, 32/33, 48/49: */	\n\t"\
		"addq	$0x40,%%rax	\n\t"\
		"addq	$0x40,%%rbx	\n\t"\
		"addq	$0x40,%%rcx	\n\t"\
		"addq	$0x40,%%rdx	\n\t"\
		"addq	$0x100,%%rsi	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps		 (%%rax),%%ymm1						\n\t		vmovaps	0x020(%%rax),%%ymm3								\n\t"\
		"vmovaps		 (%%rbx),%%ymm5						\n\t		vmovaps	0x020(%%rbx),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps		 (%%rcx),%%ymm6						\n\t		vmovaps	0x020(%%rcx),%%ymm14							\n\t"\
		"vmovaps		 (%%rdx),%%ymm5						\n\t		vmovaps	0x020(%%rdx),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm6 ,%%ymm0				\n\t		vshufpd	$15 ,%%ymm13,%%ymm14 ,%%ymm2					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm6 ,%%ymm6				\n\t		vshufpd	$0  ,%%ymm13,%%ymm14 ,%%ymm14					\n\t"\
		"vperm2f128 $32 ,%%ymm0 ,%%ymm7 ,%%ymm5				\n\t		vperm2f128 $32 ,%%ymm2 ,%%ymm15,%%ymm13		/* outB	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm0 ,%%ymm7 ,%%ymm7				\n\t		vperm2f128 $49 ,%%ymm2 ,%%ymm15,%%ymm15		/* outD	*/	\n\t"\
		"vperm2f128 $32 ,%%ymm6 ,%%ymm1 ,%%ymm0				\n\t		vperm2f128 $32 ,%%ymm14,%%ymm3 ,%%ymm2 		/* outA	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm6 ,%%ymm1 ,%%ymm1				\n\t		vperm2f128 $49 ,%%ymm14,%%ymm3 ,%%ymm3 		/* outC	*/	\n\t"\
		"vmovaps 	%%ymm5 ,0x400(%%rsi)					\n\t		vmovaps %%ymm13,0x420(%%rsi)				/* outB	*/	\n\t"\
		"vmovaps 	%%ymm7 ,0x600(%%rsi)					\n\t		vmovaps %%ymm15,0x620(%%rsi)				/* outD	*/	\n\t"\
		"vmovaps 	%%ymm0 ,     (%%rsi)					\n\t		vmovaps %%ymm2 ,0x020(%%rsi)				/* outA	*/	\n\t"\
		"vmovaps 	%%ymm1 ,0x200(%%rsi)					\n\t		vmovaps %%ymm3 ,0x220(%%rsi)				/* outC	*/	\n\t"\
	"/* a[j+p4]: Inputs from [add0,add0+0x200,add1,add1+0x200]+0x80. Outputs into r4 +0/1, 16/17, 32/33, 48/49: */	\n\t"\
		"addq	$0x40,%%rax	\n\t"\
		"addq	$0x40,%%rbx	\n\t"\
		"addq	$0x40,%%rcx	\n\t"\
		"addq	$0x40,%%rdx	\n\t"\
		"subq	$0x80,%%rsi	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps		 (%%rax),%%ymm1						\n\t		vmovaps	0x020(%%rax),%%ymm3								\n\t"\
		"vmovaps		 (%%rbx),%%ymm5						\n\t		vmovaps	0x020(%%rbx),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps		 (%%rcx),%%ymm6						\n\t		vmovaps	0x020(%%rcx),%%ymm14							\n\t"\
		"vmovaps		 (%%rdx),%%ymm5						\n\t		vmovaps	0x020(%%rdx),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm6 ,%%ymm0				\n\t		vshufpd	$15 ,%%ymm13,%%ymm14 ,%%ymm2					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm6 ,%%ymm6				\n\t		vshufpd	$0  ,%%ymm13,%%ymm14 ,%%ymm14					\n\t"\
		"vperm2f128 $32 ,%%ymm0 ,%%ymm7 ,%%ymm5				\n\t		vperm2f128 $32 ,%%ymm2 ,%%ymm15,%%ymm13		/* outB	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm0 ,%%ymm7 ,%%ymm7				\n\t		vperm2f128 $49 ,%%ymm2 ,%%ymm15,%%ymm15		/* outD	*/	\n\t"\
		"vperm2f128 $32 ,%%ymm6 ,%%ymm1 ,%%ymm0				\n\t		vperm2f128 $32 ,%%ymm14,%%ymm3 ,%%ymm2 		/* outA	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm6 ,%%ymm1 ,%%ymm1				\n\t		vperm2f128 $49 ,%%ymm14,%%ymm3 ,%%ymm3 		/* outC	*/	\n\t"\
		"vmovaps 	%%ymm5 ,0x400(%%rsi)					\n\t		vmovaps %%ymm13,0x420(%%rsi)				/* outB	*/	\n\t"\
		"vmovaps 	%%ymm7 ,0x600(%%rsi)					\n\t		vmovaps %%ymm15,0x620(%%rsi)				/* outD	*/	\n\t"\
		"vmovaps 	%%ymm0 ,     (%%rsi)					\n\t		vmovaps %%ymm2 ,0x020(%%rsi)				/* outA	*/	\n\t"\
		"vmovaps 	%%ymm1 ,0x200(%%rsi)					\n\t		vmovaps %%ymm3 ,0x220(%%rsi)				/* outC	*/	\n\t"\
	"/* a[j+p6]: Inputs from [add0,add0+0x200,add1,add1+0x200]+0xc0. Outputs into **r12** +0/1, 16/17, 32/33, 48/49: */	\n\t"\
		"addq	$0x40,%%rax	\n\t"\
		"addq	$0x40,%%rbx	\n\t"\
		"addq	$0x40,%%rcx	\n\t"\
		"addq	$0x40,%%rdx	\n\t"\
		"addq	$0x100,%%rsi	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps		 (%%rax),%%ymm1						\n\t		vmovaps	0x020(%%rax),%%ymm3								\n\t"\
		"vmovaps		 (%%rbx),%%ymm5						\n\t		vmovaps	0x020(%%rbx),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps		 (%%rcx),%%ymm6						\n\t		vmovaps	0x020(%%rcx),%%ymm14							\n\t"\
		"vmovaps		 (%%rdx),%%ymm5						\n\t		vmovaps	0x020(%%rdx),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm6 ,%%ymm0				\n\t		vshufpd	$15 ,%%ymm13,%%ymm14 ,%%ymm2					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm6 ,%%ymm6				\n\t		vshufpd	$0  ,%%ymm13,%%ymm14 ,%%ymm14					\n\t"\
		"vperm2f128 $32 ,%%ymm0 ,%%ymm7 ,%%ymm5				\n\t		vperm2f128 $32 ,%%ymm2 ,%%ymm15,%%ymm13		/* outB	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm0 ,%%ymm7 ,%%ymm7				\n\t		vperm2f128 $49 ,%%ymm2 ,%%ymm15,%%ymm15		/* outD	*/	\n\t"\
		"vperm2f128 $32 ,%%ymm6 ,%%ymm1 ,%%ymm0				\n\t		vperm2f128 $32 ,%%ymm14,%%ymm3 ,%%ymm2 		/* outA	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm6 ,%%ymm1 ,%%ymm1				\n\t		vperm2f128 $49 ,%%ymm14,%%ymm3 ,%%ymm3 		/* outC	*/	\n\t"\
		"vmovaps 	%%ymm5 ,0x400(%%rsi)					\n\t		vmovaps %%ymm13,0x420(%%rsi)				/* outB	*/	\n\t"\
		"vmovaps 	%%ymm7 ,0x600(%%rsi)					\n\t		vmovaps %%ymm15,0x620(%%rsi)				/* outD	*/	\n\t"\
		"vmovaps 	%%ymm0 ,     (%%rsi)					\n\t		vmovaps %%ymm2 ,0x020(%%rsi)				/* outA	*/	\n\t"\
		"vmovaps 	%%ymm1 ,0x200(%%rsi)					\n\t		vmovaps %%ymm3 ,0x220(%%rsi)				/* outC	*/	\n\t"\
	"/* a[j+p8]: Inputs from [add0,add0+0x200,add1,add1+0x200]+0x100. Outputs into **r2** +0/1, 16/17, 32/33, 48/49: */	\n\t"\
		"addq	$0x40,%%rax	\n\t"\
		"addq	$0x40,%%rbx	\n\t"\
		"addq	$0x40,%%rcx	\n\t"\
		"addq	$0x40,%%rdx	\n\t"\
		"subq	$0x140,%%rsi	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps		 (%%rax),%%ymm1						\n\t		vmovaps	0x020(%%rax),%%ymm3								\n\t"\
		"vmovaps		 (%%rbx),%%ymm5						\n\t		vmovaps	0x020(%%rbx),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps		 (%%rcx),%%ymm6						\n\t		vmovaps	0x020(%%rcx),%%ymm14							\n\t"\
		"vmovaps		 (%%rdx),%%ymm5						\n\t		vmovaps	0x020(%%rdx),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm6 ,%%ymm0				\n\t		vshufpd	$15 ,%%ymm13,%%ymm14 ,%%ymm2					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm6 ,%%ymm6				\n\t		vshufpd	$0  ,%%ymm13,%%ymm14 ,%%ymm14					\n\t"\
		"vperm2f128 $32 ,%%ymm0 ,%%ymm7 ,%%ymm5				\n\t		vperm2f128 $32 ,%%ymm2 ,%%ymm15,%%ymm13		/* outB	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm0 ,%%ymm7 ,%%ymm7				\n\t		vperm2f128 $49 ,%%ymm2 ,%%ymm15,%%ymm15		/* outD	*/	\n\t"\
		"vperm2f128 $32 ,%%ymm6 ,%%ymm1 ,%%ymm0				\n\t		vperm2f128 $32 ,%%ymm14,%%ymm3 ,%%ymm2 		/* outA	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm6 ,%%ymm1 ,%%ymm1				\n\t		vperm2f128 $49 ,%%ymm14,%%ymm3 ,%%ymm3 		/* outC	*/	\n\t"\
		"vmovaps 	%%ymm5 ,0x400(%%rsi)					\n\t		vmovaps %%ymm13,0x420(%%rsi)				/* outB	*/	\n\t"\
		"vmovaps 	%%ymm7 ,0x600(%%rsi)					\n\t		vmovaps %%ymm15,0x620(%%rsi)				/* outD	*/	\n\t"\
		"vmovaps 	%%ymm0 ,     (%%rsi)					\n\t		vmovaps %%ymm2 ,0x020(%%rsi)				/* outA	*/	\n\t"\
		"vmovaps 	%%ymm1 ,0x200(%%rsi)					\n\t		vmovaps %%ymm3 ,0x220(%%rsi)				/* outC	*/	\n\t"\
	"/* a[j+p10]: Inputs from [add0,add0+0x200,add1,add1+0x200]+0x140. Outputs into r10 +0/1, 16/17, 32/33, 48/49: */	\n\t"\
		"addq	$0x40,%%rax	\n\t"\
		"addq	$0x40,%%rbx	\n\t"\
		"addq	$0x40,%%rcx	\n\t"\
		"addq	$0x40,%%rdx	\n\t"\
		"addq	$0x100,%%rsi	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps		 (%%rax),%%ymm1						\n\t		vmovaps	0x020(%%rax),%%ymm3								\n\t"\
		"vmovaps		 (%%rbx),%%ymm5						\n\t		vmovaps	0x020(%%rbx),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps		 (%%rcx),%%ymm6						\n\t		vmovaps	0x020(%%rcx),%%ymm14							\n\t"\
		"vmovaps		 (%%rdx),%%ymm5						\n\t		vmovaps	0x020(%%rdx),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm6 ,%%ymm0				\n\t		vshufpd	$15 ,%%ymm13,%%ymm14 ,%%ymm2					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm6 ,%%ymm6				\n\t		vshufpd	$0  ,%%ymm13,%%ymm14 ,%%ymm14					\n\t"\
		"vperm2f128 $32 ,%%ymm0 ,%%ymm7 ,%%ymm5				\n\t		vperm2f128 $32 ,%%ymm2 ,%%ymm15,%%ymm13		/* outB	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm0 ,%%ymm7 ,%%ymm7				\n\t		vperm2f128 $49 ,%%ymm2 ,%%ymm15,%%ymm15		/* outD	*/	\n\t"\
		"vperm2f128 $32 ,%%ymm6 ,%%ymm1 ,%%ymm0				\n\t		vperm2f128 $32 ,%%ymm14,%%ymm3 ,%%ymm2 		/* outA	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm6 ,%%ymm1 ,%%ymm1				\n\t		vperm2f128 $49 ,%%ymm14,%%ymm3 ,%%ymm3 		/* outC	*/	\n\t"\
		"vmovaps 	%%ymm5 ,0x400(%%rsi)					\n\t		vmovaps %%ymm13,0x420(%%rsi)				/* outB	*/	\n\t"\
		"vmovaps 	%%ymm7 ,0x600(%%rsi)					\n\t		vmovaps %%ymm15,0x620(%%rsi)				/* outD	*/	\n\t"\
		"vmovaps 	%%ymm0 ,     (%%rsi)					\n\t		vmovaps %%ymm2 ,0x020(%%rsi)				/* outA	*/	\n\t"\
		"vmovaps 	%%ymm1 ,0x200(%%rsi)					\n\t		vmovaps %%ymm3 ,0x220(%%rsi)				/* outC	*/	\n\t"\
	"/* a[j+p12]: Inputs from [add0,add0+0x200,add1,add1+0x200]+0x180. Outputs into **r6** +0/1, 16/17, 32/33, 48/49: */	\n\t"\
		"addq	$0x40,%%rax	\n\t"\
		"addq	$0x40,%%rbx	\n\t"\
		"addq	$0x40,%%rcx	\n\t"\
		"addq	$0x40,%%rdx	\n\t"\
		"subq	$0x80,%%rsi	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps		 (%%rax),%%ymm1						\n\t		vmovaps	0x020(%%rax),%%ymm3								\n\t"\
		"vmovaps		 (%%rbx),%%ymm5						\n\t		vmovaps	0x020(%%rbx),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps		 (%%rcx),%%ymm6						\n\t		vmovaps	0x020(%%rcx),%%ymm14							\n\t"\
		"vmovaps		 (%%rdx),%%ymm5						\n\t		vmovaps	0x020(%%rdx),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm6 ,%%ymm0				\n\t		vshufpd	$15 ,%%ymm13,%%ymm14 ,%%ymm2					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm6 ,%%ymm6				\n\t		vshufpd	$0  ,%%ymm13,%%ymm14 ,%%ymm14					\n\t"\
		"vperm2f128 $32 ,%%ymm0 ,%%ymm7 ,%%ymm5				\n\t		vperm2f128 $32 ,%%ymm2 ,%%ymm15,%%ymm13		/* outB	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm0 ,%%ymm7 ,%%ymm7				\n\t		vperm2f128 $49 ,%%ymm2 ,%%ymm15,%%ymm15		/* outD	*/	\n\t"\
		"vperm2f128 $32 ,%%ymm6 ,%%ymm1 ,%%ymm0				\n\t		vperm2f128 $32 ,%%ymm14,%%ymm3 ,%%ymm2 		/* outA	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm6 ,%%ymm1 ,%%ymm1				\n\t		vperm2f128 $49 ,%%ymm14,%%ymm3 ,%%ymm3 		/* outC	*/	\n\t"\
		"vmovaps 	%%ymm5 ,0x400(%%rsi)					\n\t		vmovaps %%ymm13,0x420(%%rsi)				/* outB	*/	\n\t"\
		"vmovaps 	%%ymm7 ,0x600(%%rsi)					\n\t		vmovaps %%ymm15,0x620(%%rsi)				/* outD	*/	\n\t"\
		"vmovaps 	%%ymm0 ,     (%%rsi)					\n\t		vmovaps %%ymm2 ,0x020(%%rsi)				/* outA	*/	\n\t"\
		"vmovaps 	%%ymm1 ,0x200(%%rsi)					\n\t		vmovaps %%ymm3 ,0x220(%%rsi)				/* outC	*/	\n\t"\
	"/* a[j+p14]: Inputs from [add0,add0+0x200,add1,add1+0x200]+0x1c0. Outputs into r14 +0/1, 16/17, 32/33, 48/49: */	\n\t"\
		"addq	$0x40,%%rax	\n\t"\
		"addq	$0x40,%%rbx	\n\t"\
		"addq	$0x40,%%rcx	\n\t"\
		"addq	$0x40,%%rdx	\n\t"\
		"addq	$0x100,%%rsi	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps		 (%%rax),%%ymm1						\n\t		vmovaps	0x020(%%rax),%%ymm3								\n\t"\
		"vmovaps		 (%%rbx),%%ymm5						\n\t		vmovaps	0x020(%%rbx),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps		 (%%rcx),%%ymm6						\n\t		vmovaps	0x020(%%rcx),%%ymm14							\n\t"\
		"vmovaps		 (%%rdx),%%ymm5						\n\t		vmovaps	0x020(%%rdx),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm6 ,%%ymm0				\n\t		vshufpd	$15 ,%%ymm13,%%ymm14 ,%%ymm2					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm6 ,%%ymm6				\n\t		vshufpd	$0  ,%%ymm13,%%ymm14 ,%%ymm14					\n\t"\
		"vperm2f128 $32 ,%%ymm0 ,%%ymm7 ,%%ymm5				\n\t		vperm2f128 $32 ,%%ymm2 ,%%ymm15,%%ymm13		/* outB	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm0 ,%%ymm7 ,%%ymm7				\n\t		vperm2f128 $49 ,%%ymm2 ,%%ymm15,%%ymm15		/* outD	*/	\n\t"\
		"vperm2f128 $32 ,%%ymm6 ,%%ymm1 ,%%ymm0				\n\t		vperm2f128 $32 ,%%ymm14,%%ymm3 ,%%ymm2 		/* outA	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm6 ,%%ymm1 ,%%ymm1				\n\t		vperm2f128 $49 ,%%ymm14,%%ymm3 ,%%ymm3 		/* outC	*/	\n\t"\
		"vmovaps 	%%ymm5 ,0x400(%%rsi)					\n\t		vmovaps %%ymm13,0x420(%%rsi)				/* outB	*/	\n\t"\
		"vmovaps 	%%ymm7 ,0x600(%%rsi)					\n\t		vmovaps %%ymm15,0x620(%%rsi)				/* outD	*/	\n\t"\
		"vmovaps 	%%ymm0 ,     (%%rsi)					\n\t		vmovaps %%ymm2 ,0x020(%%rsi)				/* outA	*/	\n\t"\
		"vmovaps 	%%ymm1 ,0x200(%%rsi)					\n\t		vmovaps %%ymm3 ,0x220(%%rsi)				/* outC	*/	\n\t"\
		"\n\t"\
		/************************************************************************/\
		/* Forward DIF radix-32 pass on the interleaved block1 and block2 data: */\
		/************************************************************************/\
	/*...Block 0: */\
		"movq	%[__isrt2],%%rsi	\n\t	leaq	0x940(%%rsi),%%r8	/* one */\n\t	leaq	0x8e0(%%rsi),%%rdi	/* two */\n\t"\
	/*****	SSE2_RADIX4_DIF_4WRAPPER(c00,c08,c10,c18,r00) ****	SSE2_RADIX4_DIF_4WRAPPER_2NDOFTWO(c04,c0C,c14,c1C,r08)	*****/\
		"movq		%[__r00],%%rcx						\n\t		/*addq		$0x100,%%rcx // __r08 */	\n\t"\
		"movq		%[__c00],%%rdx						\n\t		/*addq		$0x100,%%rdx // __c04 */	\n\t"\
		"vmovaps			 (%%rcx),%%ymm0					\n\t		vmovaps		0x100(%%rcx),%%ymm8			\n\t"\
		"vmovaps		0x020(%%rcx),%%ymm1					\n\t		vmovaps		0x120(%%rcx),%%ymm9			\n\t"\
		"vmovaps		%%ymm0,%%ymm2					\n\t		vmovaps		%%ymm8,%%ymm10				\n\t"\
		"vmovaps		%%ymm1,%%ymm3					\n\t		vmovaps		%%ymm9,%%ymm11				\n\t"\
		"vmulpd			 (%%rdx),%%ymm0,%%ymm0				\n\t		vmulpd		0x100(%%rdx),%%ymm8,%%ymm8			\n\t"\
		"vmulpd			 (%%rdx),%%ymm1,%%ymm1				\n\t		vmulpd		0x100(%%rdx),%%ymm9,%%ymm9			\n\t"\
	"vfnmadd231pd	0x020(%%rdx),%%ymm7,%%ymm0				\n\t	vfnmadd231pd	0x120(%%rdx),%%ymm11,%%ymm8 	\n\t"\
	" vfmadd231pd	0x020(%%rdx),%%ymm6,%%ymm1				\n\t	 vfmadd231pd	0x120(%%rdx),%%ymm10,%%ymm9 	\n\t"\
		"vmovaps		0x040(%%rcx),%%ymm4					\n\t		vmovaps		0x140(%%rcx),%%ymm12				\n\t"\
		"vmovaps		0x060(%%rcx),%%ymm5					\n\t		vmovaps		0x160(%%rcx),%%ymm13				\n\t"\
		"vmovaps		%%ymm4,%%ymm6					\n\t		vmovaps		%%ymm12,%%ymm14						\n\t"\
		"vmovaps		%%ymm5,%%ymm7					\n\t		vmovaps		%%ymm13,%%ymm15						\n\t"\
		"addq		$0x040,%%rdx						\n\t"\
		"vmovaps		%%ymm0,%%ymm2					\n\t		vmovaps		%%ymm8,%%ymm10						\n\t"\
		"vmovaps		%%ymm1,%%ymm3					\n\t		vmovaps		%%ymm9,%%ymm11						\n\t"\
	"vmovaps	%%ymm11,(%%rbx)	\n\t	vmovaps	(%%r8),%%ymm11	\n\t"/* spill to make room for 1.0 */\
		"vmulpd			 (%%rdx),%%ymm4,%%ymm4				\n\t		vmulpd		0x100(%%rdx),%%ymm12,%%ymm12		\n\t"\
		"vmulpd			 (%%rdx),%%ymm5,%%ymm5				\n\t		vmulpd		0x100(%%rdx),%%ymm13,%%ymm13		\n\t"\
	"vfnmadd231pd	0x020(%%rdx),%%ymm7,%%ymm4				\n\t	vfnmadd231pd	0x120(%%rdx),%%ymm15,%%ymm12	\n\t"\
	" vfmadd231pd	0x020(%%rdx),%%ymm6,%%ymm5				\n\t	 vfmadd231pd	0x120(%%rdx),%%ymm14,%%ymm13	\n\t"\
		"vfmadd132pd	%%ymm11,%%ymm4,%%ymm0				\n\t		vfmadd132pd	%%ymm11,%%ymm12,%%ymm8				\n\t"\
		"vfmadd132pd	%%ymm11,%%ymm5,%%ymm1				\n\t		vfmadd132pd	%%ymm11,%%ymm13,%%ymm9				\n\t"\
		"vfmsub132pd	%%ymm11,%%ymm4,%%ymm2				\n\t		vfmsub132pd	%%ymm11,%%ymm12,%%ymm10				\n\t"\
		"vfmsub132pd	%%ymm11,%%ymm5,%%ymm3				\n\t		vfmsub132pd	(%%rbx),%%ymm13,%%ymm11				\n\t"\
		"addq		$0x080,%%rdx						\n\t"\
		"vmovaps		0x0c0(%%rcx),%%ymm4					\n\t		vmovaps		0x1c0(%%rcx),%%ymm12				\n\t"\
		"vmovaps		0x0e0(%%rcx),%%ymm5					\n\t		vmovaps		0x1e0(%%rcx),%%ymm13				\n\t"\
		"vmovaps		%%ymm4,%%ymm6					\n\t		vmovaps		%%ymm12,%%ymm14						\n\t"\
		"vmovaps		%%ymm5,%%ymm7					\n\t		vmovaps		%%ymm13,%%ymm15						\n\t"\
		"vmulpd			 (%%rdx),%%ymm4,%%ymm4				\n\t		vmulpd		0x100(%%rdx),%%ymm12,%%ymm12		\n\t"\
		"vmulpd			 (%%rdx),%%ymm5,%%ymm5				\n\t		vmulpd		0x100(%%rdx),%%ymm13,%%ymm13		\n\t"\
	"vfnmadd231pd	0x020(%%rdx),%%ymm7,%%ymm4				\n\t	vfnmadd231pd	0x120(%%rdx),%%ymm15,%%ymm12	\n\t"\
	" vfmadd231pd	0x020(%%rdx),%%ymm6,%%ymm5				\n\t	 vfmadd231pd	0x120(%%rdx),%%ymm14,%%ymm13	\n\t"\
		"vmovaps		%%ymm5,0x020(%%rcx)			\n\t		vmovaps		%%ymm13,0x120(%%rcx)				\n\t"\
		"vmovaps		%%ymm4,     (%%rcx)			\n\t		vmovaps		%%ymm12,0x100(%%rcx)				\n\t"\
		"subq		$0x040,%%rdx						\n\t"\
		"vmovaps		0x080(%%rcx),%%ymm4					\n\t		vmovaps		0x180(%%rcx),%%ymm12				\n\t"\
		"vmovaps		0x0a0(%%rcx),%%ymm5					\n\t		vmovaps		0x1a0(%%rcx),%%ymm13				\n\t"\
		"vmovaps		%%ymm4,%%ymm6					\n\t		vmovaps		%%ymm12,%%ymm14						\n\t"\
		"vmovaps		%%ymm5,%%ymm7					\n\t		vmovaps		%%ymm13,%%ymm15						\n\t"\
		"vmulpd			 (%%rdx),%%ymm4,%%ymm4				\n\t		vmulpd		0x100(%%rdx),%%ymm12,%%ymm12		\n\t"\
		"vmulpd			 (%%rdx),%%ymm5,%%ymm5				\n\t		vmulpd		0x100(%%rdx),%%ymm13,%%ymm13		\n\t"\
	"vfnmadd231pd	0x020(%%rdx),%%ymm7,%%ymm4				\n\t	vfnmadd231pd	0x120(%%rdx),%%ymm15,%%ymm12	\n\t"\
	" vfmadd231pd	0x020(%%rdx),%%ymm6,%%ymm5				\n\t	 vfmadd231pd	0x120(%%rdx),%%ymm14,%%ymm13	\n\t"\
		"vmovaps		%%ymm5,%%ymm7					\n\t		vmovaps		%%ymm13,%%ymm15					\n\t"\
		"vmovaps		%%ymm4,%%ymm6					\n\t		vmovaps		%%ymm12,%%ymm14					\n\t"\
	"vmovaps	%%ymm11,(%%rbx)	\n\t	vmovaps	(%%r8),%%ymm11	\n\t"/* spill to make room for 1.0 */\
		"vfnmadd231pd		 (%%rcx),%%ymm11,%%ymm4			\n\t		vfnmadd231pd	0x100(%%rcx),%%ymm11,%%ymm12	\n\t"\
		"vfnmadd231pd	0x020(%%rcx),%%ymm11,%%ymm5			\n\t		vfnmadd231pd	0x120(%%rcx),%%ymm11,%%ymm13	\n\t"\
		" vfmadd231pd		 (%%rcx),%%ymm11,%%ymm6			\n\t		 vfmadd231pd	0x100(%%rcx),%%ymm11,%%ymm14	\n\t"\
		" vfmadd231pd	0x020(%%rcx),%%ymm11,%%ymm7			\n\t		 vfmadd231pd	0x120(%%rcx),%%ymm11,%%ymm15	\n\t"\
		"vfmsub132pd	%%ymm11,%%ymm6,%%ymm0				\n\t		vfmsub132pd	%%ymm11,%%ymm14,%%ymm8		\n\t"\
		"vfmsub132pd	%%ymm11,%%ymm7,%%ymm1				\n\t		vfmsub132pd	%%ymm11,%%ymm15,%%ymm9		\n\t"\
		"vfmsub132pd	%%ymm11,%%ymm5,%%ymm2				\n\t		vfmsub132pd	%%ymm11,%%ymm13,%%ymm10	\n\t"\
		"vfmsub132pd	%%ymm11,%%ymm4,%%ymm3				\n\t		vfmsub132pd	(%%rbx),%%ymm12,%%ymm11	\n\t"\
	"vmovaps	%%ymm12,(%%rcx) 	\n\t"/* spill ymm12 to make room for two */"	vmovaps	(%%rdi),%%ymm12 	\n\t"/* two */\
	"vfmadd132pd	%%ymm12,%%ymm0,%%ymm6					\n\t	vfmadd132pd	%%ymm12,%%ymm8 ,%%ymm14				\n\t"\
	"vfmadd132pd	%%ymm12,%%ymm1,%%ymm7					\n\t	vfmadd132pd	%%ymm12,%%ymm9 ,%%ymm15				\n\t"\
	"vfmadd132pd	%%ymm12,%%ymm2,%%ymm5					\n\t	vfmadd132pd	%%ymm12,%%ymm10,%%ymm13				\n\t"\
	"vfmadd132pd	%%ymm12,%%ymm3,%%ymm4					\n\t	vfmadd132pd	(%%rcx),%%ymm11,%%ymm12				\n\t"\
		"																vmovaps		%%ymm14,0x100(%%rcx)			\n\t"\
		"																vmovaps		%%ymm15,0x120(%%rcx)			\n\t"\
		"																vmovaps		%%ymm10,%%ymm14					\n\t"\
		"																vmovaps		%%ymm13,%%ymm15					\n\t"\
		"																vfmsub132pd	(%%r8),%%ymm12,%%ymm10	\n\t"\
		"																vfmsub132pd	(%%r8),%%ymm11,%%ymm13	\n\t"\
		"																vfmadd132pd	(%%r8),%%ymm12,%%ymm14	\n\t"\
		"																vfmadd132pd	(%%r8),%%ymm11,%%ymm15	\n\t"\
		/*****	SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r00) *****/\
		"																vmovaps		0x100(%%rcx),%%ymm11		\n\t"\
		"																vmovaps		0x120(%%rcx),%%ymm12		\n\t"\
	"vmovaps	%%ymm13,(%%rcx) 	\n\t"/* spill ymm13 to make room for isrt2 */"	vmovaps	(%%rsi),%%ymm13 	\n\t"/*isrt2*/\
		"vfmsub132pd	(%%r8),%%ymm11,%%ymm6				\n\t	vfnmadd231pd	%%ymm13,%%ymm10,%%ymm2			\n\t"\
		"vfmsub132pd	(%%r8),%%ymm9,%%ymm0				\n\t	vfnmadd231pd	%%ymm13,%%ymm15,%%ymm5			\n\t"\
		"vfmsub132pd	(%%r8),%%ymm12,%%ymm7				\n\t	vfnmadd231pd	%%ymm13,%%ymm14,%%ymm4			\n\t"\
		"vfmsub132pd	(%%r8),%%ymm8 ,%%ymm1				\n\t	vfnmadd231pd	(%%rcx),%%ymm13,%%ymm3			\n\t"\
		"vmovaps		%%ymm6,0x100(%%rcx)			\n\t		vmovaps		%%ymm2,0x140(%%rcx)	\n\t"\
		"vmovaps		%%ymm0,0x080(%%rcx)			\n\t		vmovaps		%%ymm5,0x0c0(%%rcx)	\n\t"\
		"vmovaps		%%ymm7,0x120(%%rcx)			\n\t		vmovaps		%%ymm4,0x160(%%rcx)	\n\t"\
		"vmovaps		%%ymm1,0x1a0(%%rcx)			\n\t		vmovaps		%%ymm3,0x1e0(%%rcx)	\n\t"\
	"vmovaps	0x40(%%rdi),%%ymm13 \n\t"/* sqrt2 *//* rcol: *2 => *sqrt2 due to above FMA leaving ymm10,13-15 unmultiplied-by-isrt2: */\
	"vfmadd132pd	(%%rdi),%%ymm6,%%ymm11				\n\t	 vfmadd132pd	%%ymm13,%%ymm2,%%ymm10	\n\t"\
	"vfmadd132pd	(%%rdi),%%ymm0,%%ymm9					\n\t	 vfmadd132pd	%%ymm13,%%ymm5,%%ymm15	\n\t"\
	"vfmadd132pd	(%%rdi),%%ymm7,%%ymm12				\n\t	 vfmadd132pd	%%ymm13,%%ymm4,%%ymm14	\n\t"\
	"vfmadd132pd	(%%rdi),%%ymm1,%%ymm8					\n\t	 vfmadd132pd	(%%rcx),%%ymm3,%%ymm13	\n\t"\
		"vmovaps		%%ymm11,     (%%rcx)			\n\t		vmovaps		%%ymm10,0x040(%%rcx)	\n\t"\
		"vmovaps		%%ymm9,0x180(%%rcx)			\n\t		vmovaps		%%ymm15,0x1c0(%%rcx)	\n\t"\
		"vmovaps		%%ymm12,0x020(%%rcx)			\n\t		vmovaps		%%ymm14,0x060(%%rcx)	\n\t"\
		"vmovaps		%%ymm8,0x0a0(%%rcx)			\n\t		vmovaps		%%ymm13,0x0e0(%%rcx)	\n\t"\
		"\n\t"\
	/*...Block 2: */\
	/*****	SSE2_RADIX4_DIF_4WRAPPER(c02,c0A,c12,c1A,r10) ****	SSE2_RADIX4_DIF_4WRAPPER_2NDOFTWO(c06,c0E,c16,c1E,r18)	*****/\
		"movq		%[__r10],%%rcx						\n\t		/*addq		$0x100,%%rcx // __r18 */	\n\t"\
		"movq		%[__c02],%%rdx						\n\t		/*addq		$0x100,%%rdx // __c06 */	\n\t"\
		"vmovaps			 (%%rcx),%%ymm0					\n\t		vmovaps		0x100(%%rcx),%%ymm8			\n\t"\
		"vmovaps		0x020(%%rcx),%%ymm1					\n\t		vmovaps		0x120(%%rcx),%%ymm9			\n\t"\
		"vmovaps		%%ymm0,%%ymm2					\n\t		vmovaps		%%ymm8,%%ymm10				\n\t"\
		"vmovaps		%%ymm1,%%ymm3					\n\t		vmovaps		%%ymm9,%%ymm11				\n\t"\
		"vmulpd			 (%%rdx),%%ymm0,%%ymm0				\n\t		vmulpd		0x100(%%rdx),%%ymm8,%%ymm8			\n\t"\
		"vmulpd			 (%%rdx),%%ymm1,%%ymm1				\n\t		vmulpd		0x100(%%rdx),%%ymm9,%%ymm9			\n\t"\
	"vfnmadd231pd	0x020(%%rdx),%%ymm3,%%ymm0				\n\t	vfnmadd231pd	0x120(%%rdx),%%ymm11,%%ymm8 	\n\t"\
	" vfmadd231pd	0x020(%%rdx),%%ymm2,%%ymm1				\n\t	 vfmadd231pd	0x120(%%rdx),%%ymm10,%%ymm9 	\n\t"\
		"vmovaps		0x040(%%rcx),%%ymm4					\n\t		vmovaps		0x140(%%rcx),%%ymm12				\n\t"\
		"vmovaps		0x060(%%rcx),%%ymm5					\n\t		vmovaps		0x160(%%rcx),%%ymm13				\n\t"\
		"vmovaps		%%ymm4,%%ymm6					\n\t		vmovaps		%%ymm12,%%ymm14						\n\t"\
		"vmovaps		%%ymm5,%%ymm7					\n\t		vmovaps		%%ymm13,%%ymm15						\n\t"\
		"addq		$0x040,%%rdx						\n\t"\
		"vmovaps		%%ymm0,%%ymm2					\n\t		vmovaps		%%ymm8,%%ymm10						\n\t"\
		"vmovaps		%%ymm1,%%ymm3					\n\t		vmovaps		%%ymm9,%%ymm11						\n\t"\
	"vmovaps	%%ymm11,(%%rbx)	\n\t	vmovaps	(%%r8),%%ymm11	\n\t"/* spill to make room for 1.0 */\
		"vmulpd			 (%%rdx),%%ymm4,%%ymm4				\n\t		vmulpd		0x100(%%rdx),%%ymm12,%%ymm12		\n\t"\
		"vmulpd			 (%%rdx),%%ymm5,%%ymm5				\n\t		vmulpd		0x100(%%rdx),%%ymm13,%%ymm13		\n\t"\
	"vfnmadd231pd	0x020(%%rdx),%%ymm7,%%ymm4				\n\t	vfnmadd231pd	0x120(%%rdx),%%ymm15,%%ymm12	\n\t"\
	" vfmadd231pd	0x020(%%rdx),%%ymm6,%%ymm5				\n\t	 vfmadd231pd	0x120(%%rdx),%%ymm14,%%ymm13	\n\t"\
		"vfmadd132pd	%%ymm11,%%ymm4,%%ymm0				\n\t		vfmadd132pd	%%ymm11,%%ymm12,%%ymm8				\n\t"\
		"vfmadd132pd	%%ymm11,%%ymm5,%%ymm1				\n\t		vfmadd132pd	%%ymm11,%%ymm13,%%ymm9				\n\t"\
		"vfmsub132pd	%%ymm11,%%ymm4,%%ymm2				\n\t		vfmsub132pd	%%ymm11,%%ymm12,%%ymm10				\n\t"\
		"vfmsub132pd	%%ymm11,%%ymm5,%%ymm3				\n\t		vfmsub132pd	(%%rbx),%%ymm13,%%ymm11				\n\t"\
		"addq		$0x080,%%rdx						\n\t"\
		"vmovaps		0x0c0(%%rcx),%%ymm4					\n\t		vmovaps		0x1c0(%%rcx),%%ymm12				\n\t"\
		"vmovaps		0x0e0(%%rcx),%%ymm5					\n\t		vmovaps		0x1e0(%%rcx),%%ymm13				\n\t"\
		"vmovaps		%%ymm4,%%ymm6					\n\t		vmovaps		%%ymm12,%%ymm14						\n\t"\
		"vmovaps		%%ymm5,%%ymm7					\n\t		vmovaps		%%ymm13,%%ymm15						\n\t"\
		"vmulpd			 (%%rdx),%%ymm4,%%ymm4				\n\t		vmulpd		0x100(%%rdx),%%ymm12,%%ymm12		\n\t"\
		"vmulpd			 (%%rdx),%%ymm5,%%ymm5				\n\t		vmulpd		0x100(%%rdx),%%ymm13,%%ymm13		\n\t"\
	"vfnmadd231pd	0x020(%%rdx),%%ymm7,%%ymm4				\n\t	vfnmadd231pd	0x120(%%rdx),%%ymm15,%%ymm12	\n\t"\
	" vfmadd231pd	0x020(%%rdx),%%ymm6,%%ymm5				\n\t	 vfmadd231pd	0x120(%%rdx),%%ymm14,%%ymm13	\n\t"\
		"vmovaps		%%ymm5,0x020(%%rcx)			\n\t		vmovaps		%%ymm13,0x120(%%rcx)				\n\t"\
		"vmovaps		%%ymm4,     (%%rcx)			\n\t		vmovaps		%%ymm12,0x100(%%rcx)				\n\t"\
		"subq		$0x040,%%rdx						\n\t"\
		"vmovaps		0x080(%%rcx),%%ymm4					\n\t		vmovaps		0x180(%%rcx),%%ymm12				\n\t"\
		"vmovaps		0x0a0(%%rcx),%%ymm5					\n\t		vmovaps		0x1a0(%%rcx),%%ymm13				\n\t"\
		"vmovaps		%%ymm4,%%ymm6					\n\t		vmovaps		%%ymm12,%%ymm14						\n\t"\
		"vmovaps		%%ymm5,%%ymm7					\n\t		vmovaps		%%ymm13,%%ymm15						\n\t"\
		"vmulpd			 (%%rdx),%%ymm4,%%ymm4				\n\t		vmulpd		0x100(%%rdx),%%ymm12,%%ymm12		\n\t"\
		"vmulpd			 (%%rdx),%%ymm5,%%ymm5				\n\t		vmulpd		0x100(%%rdx),%%ymm13,%%ymm13		\n\t"\
	"vfnmadd231pd	0x020(%%rdx),%%ymm7,%%ymm4				\n\t	vfnmadd231pd	0x120(%%rdx),%%ymm15,%%ymm12	\n\t"\
	" vfmadd231pd	0x020(%%rdx),%%ymm6,%%ymm5				\n\t	 vfmadd231pd	0x120(%%rdx),%%ymm14,%%ymm13	\n\t"\
		"vmovaps		%%ymm5,%%ymm7					\n\t		vmovaps		%%ymm13,%%ymm15						\n\t"\
		"vmovaps		%%ymm4,%%ymm6					\n\t		vmovaps		%%ymm12,%%ymm14						\n\t"\
	"vmovaps	%%ymm11,(%%rbx)	\n\t	vmovaps	(%%r8),%%ymm11	\n\t"/* spill to make room for 1.0 */\
		"vfnmadd231pd		 (%%rcx),%%ymm11,%%ymm4			\n\t		vfnmadd231pd	0x100(%%rcx),%%ymm11,%%ymm12	\n\t"\
		"vfnmadd231pd	0x020(%%rcx),%%ymm11,%%ymm5			\n\t		vfnmadd231pd	0x120(%%rcx),%%ymm11,%%ymm13	\n\t"\
		" vfmadd231pd		 (%%rcx),%%ymm11,%%ymm6			\n\t		 vfmadd231pd	0x100(%%rcx),%%ymm11,%%ymm14	\n\t"\
		" vfmadd231pd	0x020(%%rcx),%%ymm11,%%ymm7			\n\t		 vfmadd231pd	0x120(%%rcx),%%ymm11,%%ymm15	\n\t"\
		"vfmsub132pd	%%ymm11,%%ymm6,%%ymm0				\n\t		vfmsub132pd	%%ymm11,%%ymm14,%%ymm8		\n\t"\
		"vfmsub132pd	%%ymm11,%%ymm7,%%ymm1				\n\t		vfmsub132pd	%%ymm11,%%ymm15,%%ymm9		\n\t"\
		"vfmsub132pd	%%ymm11,%%ymm5,%%ymm2				\n\t		vfmsub132pd	%%ymm11,%%ymm13,%%ymm10	\n\t"\
		"vfmsub132pd	%%ymm11,%%ymm4,%%ymm3				\n\t		vfmsub132pd	(%%rbx),%%ymm12,%%ymm11	\n\t"\
	"vmovaps	%%ymm12,(%%rcx) 	\n\t"/* spill ymm12 to make room for two */"	vmovaps	(%%rdi),%%ymm12 	\n\t"/* two */\
	"vfmadd132pd	%%ymm12,%%ymm0,%%ymm6					\n\t	vfmadd132pd	%%ymm12,%%ymm8 ,%%ymm14				\n\t"\
	"vfmadd132pd	%%ymm12,%%ymm1,%%ymm7					\n\t	vfmadd132pd	%%ymm12,%%ymm9 ,%%ymm15				\n\t"\
	"vfmadd132pd	%%ymm12,%%ymm2,%%ymm5					\n\t	vfmadd132pd	%%ymm12,%%ymm10,%%ymm13				\n\t"\
	"vfmadd132pd	%%ymm12,%%ymm3,%%ymm4					\n\t	vfmadd132pd	(%%rcx),%%ymm11,%%ymm12				\n\t"\
		"																vmovaps		%%ymm14,0x100(%%rcx)			\n\t"\
		"																vmovaps		%%ymm15,0x120(%%rcx)			\n\t"\
		"																vmovaps		%%ymm10,%%ymm14					\n\t"\
		"																vmovaps		%%ymm13,%%ymm15					\n\t"\
		"																vfmsub132pd	(%%r8),%%ymm12,%%ymm10	\n\t"\
		"																vfmsub132pd	(%%r8),%%ymm11,%%ymm13	\n\t"\
		"																vfmadd132pd	(%%r8),%%ymm12,%%ymm14	\n\t"\
		"																vfmadd132pd	(%%r8),%%ymm11,%%ymm15	\n\t"\
		/*****	SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r10) *****/\
		"																vmovaps		0x100(%%rcx),%%ymm11		\n\t"\
		"																vmovaps		0x120(%%rcx),%%ymm12		\n\t"\
	"vmovaps	%%ymm13,(%%rcx) 	\n\t"/* spill ymm13 to make room for isrt2 */"	vmovaps	(%%rsi),%%ymm13 	\n\t"/*isrt2*/\
		"vfmsub132pd	(%%r8),%%ymm11,%%ymm6				\n\t	vfnmadd231pd	%%ymm13,%%ymm10,%%ymm2			\n\t"\
		"vfmsub132pd	(%%r8),%%ymm9 ,%%ymm0				\n\t	vfnmadd231pd	%%ymm13,%%ymm15,%%ymm5			\n\t"\
		"vfmsub132pd	(%%r8),%%ymm12,%%ymm7				\n\t	vfnmadd231pd	%%ymm13,%%ymm14,%%ymm4			\n\t"\
		"vfmsub132pd	(%%r8),%%ymm8 ,%%ymm1				\n\t	vfnmadd231pd	(%%rcx),%%ymm13,%%ymm3			\n\t"\
		"vmovaps		%%ymm6,0x100(%%rcx)			\n\t		vmovaps		%%ymm2,0x140(%%rcx)	\n\t"\
		"vmovaps		%%ymm0,0x080(%%rcx)			\n\t		vmovaps		%%ymm5,0x0c0(%%rcx)	\n\t"\
		"vmovaps		%%ymm7,0x120(%%rcx)			\n\t		vmovaps		%%ymm4,0x160(%%rcx)	\n\t"\
		"vmovaps		%%ymm1,0x1a0(%%rcx)			\n\t		vmovaps		%%ymm3,0x1e0(%%rcx)	\n\t"\
	"vmovaps	0x40(%%rdi),%%ymm13 \n\t"/* sqrt2 *//* rcol: *2 => *sqrt2 due to above FMA leaving ymm10,13-15 unmultiplied-by-isrt2: */\
	"vfmadd132pd	(%%rdi),%%ymm6,%%ymm11				\n\t	 vfmadd132pd	%%ymm13,%%ymm2,%%ymm10	\n\t"\
	"vfmadd132pd	(%%rdi),%%ymm0,%%ymm9					\n\t	 vfmadd132pd	%%ymm13,%%ymm5,%%ymm15	\n\t"\
	"vfmadd132pd	(%%rdi),%%ymm7,%%ymm12				\n\t	 vfmadd132pd	%%ymm13,%%ymm4,%%ymm14	\n\t"\
	"vfmadd132pd	(%%rdi),%%ymm1,%%ymm8					\n\t	 vfmadd132pd	(%%rcx),%%ymm3,%%ymm13	\n\t"\
		"vmovaps		%%ymm11,     (%%rcx)			\n\t		vmovaps		%%ymm10,0x040(%%rcx)	\n\t"\
		"vmovaps		%%ymm9,0x180(%%rcx)			\n\t		vmovaps		%%ymm15,0x1c0(%%rcx)	\n\t"\
		"vmovaps		%%ymm12,0x020(%%rcx)			\n\t		vmovaps		%%ymm14,0x060(%%rcx)	\n\t"\
		"vmovaps		%%ymm8,0x0a0(%%rcx)			\n\t		vmovaps		%%ymm13,0x0e0(%%rcx)	\n\t"\
		"\n\t"\
	/************************************************************************************************************/\
	/* Next 2 blocks operate on odd-indexed elements from the unpck*pd commands which we stored to temporaries: */\
	/************************************************************************************************************/\
	/*...Block 3: */\
	/*	SSE2_RADIX4_DIF_4TWIDDLE         (r20,r24,r22,r26,r20,c01) */\
		"addq		$0x200,%%rcx		\n\t"	/***	SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO(r28,r2C,r2A,r2E,r28,c05)	*****/\
		"movq		%[__c01],%%rbx						\n\t	/*	movq		%[__c05],%%rbx	*/		\n\t"\
		"movq		%%rcx,%%rax						\n\t	/*	addq		$0x080,%%rax	*/		\n\t"\
		"addq		$0x040,%%rcx						\n\t	/*	addq		$0x080,%%rcx	*/		\n\t"\
		"vmovaps		 (%%rax),%%ymm0	\n\t movq %%rax,%%rdx \n\t	vmovaps		0x100(%%rax),%%ymm8			\n\t"\
		"vmovaps		 (%%rcx),%%ymm4						\n\t		vmovaps		0x100(%%rcx),%%ymm12		\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1						\n\t		vmovaps		0x120(%%rax),%%ymm9			\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5						\n\t		vmovaps		0x120(%%rcx),%%ymm13		\n\t"\
		"vmovaps		 (%%rbx),%%ymm6						\n\t		vmovaps		0x100(%%rbx),%%ymm14		\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm7						\n\t		vmovaps		0x120(%%rbx),%%ymm15		\n\t"\
		"vmovaps	%%ymm0,%%ymm2						\n\t		vmovaps		%%ymm8,%%ymm10		\n\t"\
		"vmovaps	%%ymm1,%%ymm3						\n\t		vmovaps		%%ymm9,%%ymm11		\n\t"\
		"vmulpd		%%ymm6,%%ymm0,%%ymm0				\n\t		vmulpd		%%ymm14,%%ymm8,%%ymm8			\n\t"\
		"vmulpd		%%ymm6,%%ymm1,%%ymm1				\n\t		vmulpd		%%ymm14,%%ymm9,%%ymm9			\n\t"\
	"vfnmadd231pd	%%ymm7,%%ymm3,%%ymm0					\n\t	vfnmadd231pd	%%ymm15,%%ymm11,%%ymm8 	\n\t"\
	" vfmadd231pd	%%ymm7,%%ymm2,%%ymm1					\n\t	 vfmadd231pd	%%ymm15,%%ymm10,%%ymm9 	\n\t"\
		"vmovaps		%%ymm4,%%ymm6					\n\t		vmovaps		%%ymm12,%%ymm14				\n\t"\
		"vmovaps		%%ymm5,%%ymm7					\n\t		vmovaps		%%ymm13,%%ymm15				\n\t"\
		"vmulpd		0x040(%%rbx),%%ymm4,%%ymm4				\n\t		vmulpd		0x140(%%rbx),%%ymm12,%%ymm12		\n\t"\
		"vmulpd		0x040(%%rbx),%%ymm5,%%ymm5				\n\t		vmulpd		0x140(%%rbx),%%ymm13,%%ymm13		\n\t"\
		"vmovaps		%%ymm0,%%ymm2					\n\t		vmovaps		%%ymm8,%%ymm10				\n\t"\
		"vmovaps		%%ymm1,%%ymm3					\n\t		vmovaps		%%ymm9,%%ymm11				\n\t"\
	"vfnmadd231pd	0x060(%%rbx),%%ymm7,%%ymm4				\n\t	vfnmadd231pd	0x160(%%rbx),%%ymm15,%%ymm12	\n\t"\
	" vfmadd231pd	0x060(%%rbx),%%ymm6,%%ymm5				\n\t	 vfmadd231pd	0x160(%%rbx),%%ymm14,%%ymm13	\n\t"\
		"addq		$0x080,%%rcx						\n\t		addq		$0x0c0,%%rbx				\n\t"\
	"vmovaps	%%ymm11,(%%rdx)	\n\t	vmovaps	(%%r8),%%ymm11	\n\t"/* spill to make room for 1.0 */\
		"vmovaps			 (%%rcx),%%ymm6					\n\t"\
		"vmovaps		0x020(%%rcx),%%ymm7					\n\t"\
		"vfmadd132pd	%%ymm11,%%ymm4,%%ymm0				\n\t		vfmadd132pd	%%ymm11,%%ymm12,%%ymm8		\n\t"\
		"vfmadd132pd	%%ymm11,%%ymm5,%%ymm1				\n\t		vfmadd132pd	%%ymm11,%%ymm13,%%ymm9		\n\t"\
		"vfmsub132pd	%%ymm11,%%ymm4,%%ymm2				\n\t		vfmsub132pd	%%ymm11,%%ymm12,%%ymm10	\n\t"\
		"vfmsub132pd	%%ymm11,%%ymm5,%%ymm3				\n\t		vfmsub132pd	(%%rdx),%%ymm13,%%ymm11	\n\t"\
		"vmovaps		%%ymm6,%%ymm4					\n\t		vmovaps		0x100(%%rcx),%%ymm12				\n\t"\
		"vmovaps		%%ymm7,%%ymm5					\n\t		vmovaps		0x120(%%rcx),%%ymm13				\n\t"\
		"vmulpd			 (%%rbx),%%ymm4,%%ymm4				\n\t		vmovaps		0x100(%%rcx),%%ymm14				\n\t"\
		"vmulpd			 (%%rbx),%%ymm5,%%ymm5				\n\t		vmovaps		0x120(%%rcx),%%ymm15				\n\t"\
	"vfnmadd231pd	0x020(%%rbx),%%ymm7,%%ymm4				\n\t		vmulpd		0x100(%%rbx),%%ymm12,%%ymm12		\n\t"\
	" vfmadd231pd	0x020(%%rbx),%%ymm6,%%ymm5				\n\t		vmulpd		0x100(%%rbx),%%ymm13,%%ymm13		\n\t"\
		"vmovaps		%%ymm5,0x020(%%rdx)			\n\t	vfnmadd231pd	0x120(%%rbx),%%ymm15,%%ymm12	\n\t"\
		"vmovaps		%%ymm4,     (%%rdx)			\n\t	 vfmadd231pd	0x120(%%rbx),%%ymm14,%%ymm13	\n\t"\
		"addq	$0x080,%%rax								\n\t		vmovaps		%%ymm13,0x120(%%rdx)				\n\t"\
		"subq	$0x040,%%rbx								\n\t		vmovaps		%%ymm12,0x100(%%rdx)				\n\t"\
		"vmovaps			 (%%rax),%%ymm4					\n\t		vmovaps		0x100(%%rax),%%ymm12				\n\t"\
		"vmovaps		0x020(%%rax),%%ymm5					\n\t		vmovaps		0x120(%%rax),%%ymm13				\n\t"\
		"vmovaps		%%ymm4,%%ymm6					\n\t		vmovaps		0x100(%%rax),%%ymm14				\n\t"\
		"vmovaps		%%ymm5,%%ymm7					\n\t		vmovaps		0x120(%%rax),%%ymm15				\n\t"\
		"vmulpd			 (%%rbx),%%ymm4,%%ymm4				\n\t		vmulpd		0x100(%%rbx),%%ymm12,%%ymm12		\n\t"\
		"vmulpd			 (%%rbx),%%ymm5,%%ymm5				\n\t		vmulpd		0x100(%%rbx),%%ymm13,%%ymm13		\n\t"\
	"vfnmadd231pd	0x020(%%rbx),%%ymm7,%%ymm4				\n\t	vfnmadd231pd	0x120(%%rbx),%%ymm15,%%ymm12	\n\t"\
	" vfmadd231pd	0x020(%%rbx),%%ymm6,%%ymm5				\n\t	 vfmadd231pd	0x120(%%rbx),%%ymm14,%%ymm13	\n\t"\
		"vmovaps		%%ymm5,%%ymm7					\n\t		vmovaps		%%ymm13,%%ymm15				\n\t"\
		"vmovaps		%%ymm4,%%ymm6					\n\t		vmovaps		%%ymm12,%%ymm14				\n\t"\
	"vmovaps	%%ymm11,(%%rax)	\n\t	vmovaps	(%%r8),%%ymm11	\n\t"/* spill to make room for 1.0 */\
		"vfnmadd231pd		 (%%rdx),%%ymm11,%%ymm4			\n\t		vfnmadd231pd	0x100(%%rdx),%%ymm11,%%ymm12	\n\t"\
		"vfnmadd231pd	0x020(%%rdx),%%ymm11,%%ymm5			\n\t		vfnmadd231pd	0x120(%%rdx),%%ymm11,%%ymm13	\n\t"\
		" vfmadd231pd		 (%%rdx),%%ymm11,%%ymm6			\n\t		 vfmadd231pd	0x100(%%rdx),%%ymm11,%%ymm14	\n\t"\
		" vfmadd231pd	0x020(%%rdx),%%ymm11,%%ymm7			\n\t		 vfmadd231pd	0x120(%%rdx),%%ymm11,%%ymm15	\n\t"\
		"vfmsub132pd	%%ymm11,%%ymm6,%%ymm0				\n\t		vfmsub132pd	%%ymm11,%%ymm14,%%ymm8		\n\t"\
		"vfmsub132pd	%%ymm11,%%ymm7,%%ymm1				\n\t		vfmsub132pd	%%ymm11,%%ymm15,%%ymm9		\n\t"\
		"vfmsub132pd	%%ymm11,%%ymm5,%%ymm2				\n\t		vfmsub132pd	%%ymm11,%%ymm13,%%ymm10	\n\t"\
		"vfmsub132pd	%%ymm11,%%ymm4,%%ymm3				\n\t		vfmsub132pd	(%%rax),%%ymm12,%%ymm11	\n\t"\
	"vmovaps	%%ymm12,(%%rdx) 	\n\t"/* spill ymm12 to make room for two */"	vmovaps	(%%rdi),%%ymm12 	\n\t"/* two */\
	"vfmadd132pd	%%ymm12,%%ymm0,%%ymm6					\n\t	vfmadd132pd	%%ymm12,%%ymm8 ,%%ymm14				\n\t"\
	"vfmadd132pd	%%ymm12,%%ymm1,%%ymm7					\n\t	vfmadd132pd	%%ymm12,%%ymm9 ,%%ymm15				\n\t"\
	"vfmadd132pd	%%ymm12,%%ymm2,%%ymm5					\n\t	vfmadd132pd	%%ymm12,%%ymm10,%%ymm13				\n\t"\
	"vfmadd132pd	%%ymm12,%%ymm3,%%ymm4					\n\t	vfmadd132pd	(%%rdx),%%ymm11,%%ymm12				\n\t"\
		"																vmovaps		%%ymm14,0x100(%%rdx)			\n\t"\
		"																vmovaps		%%ymm15,0x120(%%rdx)			\n\t"\
		"																vmovaps		%%ymm10,%%ymm14					\n\t"\
		"																vmovaps		%%ymm13,%%ymm15					\n\t"\
		"																vfmsub132pd	(%%r8),%%ymm12,%%ymm10	\n\t"\
		"																vfmsub132pd	(%%r8),%%ymm11,%%ymm13	\n\t"\
		"																vfmadd132pd	(%%r8),%%ymm12,%%ymm14	\n\t"\
		"																vfmadd132pd	(%%r8),%%ymm11,%%ymm15	\n\t"\
		/*****	SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r20) *****/\
		"																vmovaps		0x100(%%rdx),%%ymm11		\n\t"\
		"																vmovaps		0x120(%%rdx),%%ymm12		\n\t"\
	"vmovaps	%%ymm13,(%%rdx) 	\n\t"/* spill ymm13 to make room for isrt2 */"	vmovaps	(%%rsi),%%ymm13 	\n\t"/*isrt2*/\
		"vfmsub132pd	(%%r8),%%ymm11,%%ymm6				\n\t	vfnmadd231pd	%%ymm13,%%ymm10,%%ymm2			\n\t"\
		"vfmsub132pd	(%%r8),%%ymm9,%%ymm0				\n\t	vfnmadd231pd	%%ymm13,%%ymm15,%%ymm5			\n\t"\
		"vfmsub132pd	(%%r8),%%ymm12,%%ymm7				\n\t	vfnmadd231pd	%%ymm13,%%ymm14,%%ymm4			\n\t"\
		"vfmsub132pd	(%%r8),%%ymm8 ,%%ymm1				\n\t	vfnmadd231pd	(%%rdx),%%ymm13,%%ymm3			\n\t"\
		"vmovaps		%%ymm6,0x100(%%rdx)			\n\t		vmovaps		%%ymm2,0x140(%%rdx)	\n\t"\
		"vmovaps		%%ymm0,0x080(%%rdx)			\n\t		vmovaps		%%ymm5,0x0c0(%%rdx)	\n\t"\
		"vmovaps		%%ymm7,0x120(%%rdx)			\n\t		vmovaps		%%ymm4,0x160(%%rdx)	\n\t"\
		"vmovaps		%%ymm1,0x1a0(%%rdx)			\n\t		vmovaps		%%ymm3,0x1e0(%%rdx)	\n\t"\
	"vmovaps	0x40(%%rdi),%%ymm13 \n\t"/* sqrt2 *//* rcol: *2 => *sqrt2 due to above FMA leaving ymm10,13-15 unmultiplied-by-isrt2: */\
	"vfmadd132pd	(%%rdi),%%ymm6,%%ymm11				\n\t	 vfmadd132pd	%%ymm13,%%ymm2,%%ymm10	\n\t"\
	"vfmadd132pd	(%%rdi),%%ymm0,%%ymm9					\n\t	 vfmadd132pd	%%ymm13,%%ymm5,%%ymm15	\n\t"\
	"vfmadd132pd	(%%rdi),%%ymm7,%%ymm12				\n\t	 vfmadd132pd	%%ymm13,%%ymm4,%%ymm14	\n\t"\
	"vfmadd132pd	(%%rdi),%%ymm1,%%ymm8					\n\t	 vfmadd132pd	(%%rdx),%%ymm3,%%ymm13	\n\t"\
		"vmovaps		%%ymm11,     (%%rdx)			\n\t		vmovaps		%%ymm10,0x040(%%rdx)	\n\t"\
		"vmovaps		%%ymm9,0x180(%%rdx)			\n\t		vmovaps		%%ymm15,0x1c0(%%rdx)	\n\t"\
		"vmovaps		%%ymm12,0x020(%%rdx)			\n\t		vmovaps		%%ymm14,0x060(%%rdx)	\n\t"\
		"vmovaps		%%ymm8,0x0a0(%%rdx)			\n\t		vmovaps		%%ymm13,0x0e0(%%rdx)	\n\t"\
		"\n\t"\
	/*...Block 4: */\
		"/*****	SSE2_RADIX4_DIF_4TWIDDLE         (r30,r34,r32,r36,r30,c03)	*****/\n\t"\
		"movq		%[__c03],%%rbx					\n\t"\
		"movq		%[__r30],%%rax					\n\t"\
		"movq		%%rax,%%rcx		/*****	SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO(r38,r3C,r3A,r3E,r38,c07)	*****/\n\t"\
		"addq		$0x040,%%rcx					\n\t"\
		"vmovaps		 (%%rax),%%ymm0	\n\t movq %%rax,%%rdx \n\t	vmovaps		0x100(%%rax),%%ymm8			\n\t"\
		"vmovaps		 (%%rcx),%%ymm4						\n\t		vmovaps		0x100(%%rcx),%%ymm12		\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1						\n\t		vmovaps		0x120(%%rax),%%ymm9			\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5						\n\t		vmovaps		0x120(%%rcx),%%ymm13		\n\t"\
		"vmovaps		 (%%rbx),%%ymm6						\n\t		vmovaps		0x100(%%rbx),%%ymm14		\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm7						\n\t		vmovaps		0x120(%%rbx),%%ymm15		\n\t"\
		"vmovaps	%%ymm0,%%ymm2						\n\t		vmovaps		%%ymm8,%%ymm10		\n\t"\
		"vmovaps	%%ymm1,%%ymm3						\n\t		vmovaps		%%ymm9,%%ymm11		\n\t"\
		"vmulpd		%%ymm6,%%ymm0,%%ymm0				\n\t		vmulpd		%%ymm14,%%ymm8,%%ymm8			\n\t"\
		"vmulpd		%%ymm6,%%ymm1,%%ymm1				\n\t		vmulpd		%%ymm14,%%ymm9,%%ymm9			\n\t"\
	"vfnmadd231pd	%%ymm7,%%ymm3,%%ymm0			\n\t	vfnmadd231pd	%%ymm15,%%ymm11,%%ymm8 	\n\t"\
	" vfmadd231pd	%%ymm7,%%ymm2,%%ymm1			\n\t	 vfmadd231pd	%%ymm15,%%ymm10,%%ymm9 	\n\t"\
		"vmovaps		%%ymm4,%%ymm6					\n\t		vmovaps		%%ymm12,%%ymm14				\n\t"\
		"vmovaps		%%ymm5,%%ymm7					\n\t		vmovaps		%%ymm13,%%ymm15				\n\t"\
		"vmulpd		0x040(%%rbx),%%ymm4,%%ymm4				\n\t		vmulpd		0x140(%%rbx),%%ymm12,%%ymm12		\n\t"\
		"vmulpd		0x040(%%rbx),%%ymm5,%%ymm5				\n\t		vmulpd		0x140(%%rbx),%%ymm13,%%ymm13		\n\t"\
		"vmovaps		%%ymm0,%%ymm2					\n\t		vmovaps		%%ymm8,%%ymm10				\n\t"\
		"vmovaps		%%ymm1,%%ymm3					\n\t		vmovaps		%%ymm9,%%ymm11				\n\t"\
	"vfnmadd231pd	0x060(%%rbx),%%ymm7,%%ymm4				\n\t	vfnmadd231pd	0x160(%%rbx),%%ymm15,%%ymm12	\n\t"\
	" vfmadd231pd	0x060(%%rbx),%%ymm6,%%ymm5				\n\t	 vfmadd231pd	0x160(%%rbx),%%ymm14,%%ymm13	\n\t"\
		"addq		$0x080,%%rcx						\n\t		addq		$0x0c0,%%rbx				\n\t"\
	"vmovaps	%%ymm11,(%%rdx)	\n\t	vmovaps	(%%r8),%%ymm11	\n\t"/* spill to make room for 1.0 */\
		"vmovaps			 (%%rcx),%%ymm6					\n\t"\
		"vmovaps		0x020(%%rcx),%%ymm7					\n\t"\
		"vfmadd132pd	%%ymm11,%%ymm4,%%ymm0				\n\t		vfmadd132pd	%%ymm11,%%ymm12,%%ymm8		\n\t"\
		"vfmadd132pd	%%ymm11,%%ymm5,%%ymm1				\n\t		vfmadd132pd	%%ymm11,%%ymm13,%%ymm9		\n\t"\
		"vfmsub132pd	%%ymm11,%%ymm4,%%ymm2				\n\t		vfmsub132pd	%%ymm11,%%ymm12,%%ymm10	\n\t"\
		"vfmsub132pd	%%ymm11,%%ymm5,%%ymm3				\n\t		vfmsub132pd	(%%rdx),%%ymm13,%%ymm11	\n\t"\
		"vmovaps		%%ymm6,%%ymm4					\n\t		vmovaps		0x100(%%rcx),%%ymm12				\n\t"\
		"vmovaps		%%ymm7,%%ymm5					\n\t		vmovaps		0x120(%%rcx),%%ymm13				\n\t"\
		"vmulpd			 (%%rbx),%%ymm4,%%ymm4				\n\t		vmovaps		0x100(%%rcx),%%ymm14				\n\t"\
		"vmulpd			 (%%rbx),%%ymm5,%%ymm5				\n\t		vmovaps		0x120(%%rcx),%%ymm15				\n\t"\
	"vfnmadd231pd	0x020(%%rbx),%%ymm7,%%ymm4				\n\t		vmulpd		0x100(%%rbx),%%ymm12,%%ymm12		\n\t"\
	" vfmadd231pd	0x020(%%rbx),%%ymm6,%%ymm5				\n\t		vmulpd		0x100(%%rbx),%%ymm13,%%ymm13		\n\t"\
		"vmovaps		%%ymm5,0x020(%%rdx)			\n\t	vfnmadd231pd	0x120(%%rbx),%%ymm15,%%ymm12	\n\t"\
		"vmovaps		%%ymm4,     (%%rdx)			\n\t	 vfmadd231pd	0x120(%%rbx),%%ymm14,%%ymm13	\n\t"\
		"addq	$0x080,%%rax								\n\t		vmovaps		%%ymm13,0x120(%%rdx)				\n\t"\
		"subq	$0x040,%%rbx								\n\t		vmovaps		%%ymm12,0x100(%%rdx)				\n\t"\
		"vmovaps			 (%%rax),%%ymm4					\n\t		vmovaps		0x100(%%rax),%%ymm12				\n\t"\
		"vmovaps		0x020(%%rax),%%ymm5					\n\t		vmovaps		0x120(%%rax),%%ymm13				\n\t"\
		"vmovaps		%%ymm4,%%ymm6					\n\t		vmovaps		0x100(%%rax),%%ymm14				\n\t"\
		"vmovaps		%%ymm5,%%ymm7					\n\t		vmovaps		0x120(%%rax),%%ymm15				\n\t"\
		"vmulpd			 (%%rbx),%%ymm4,%%ymm4				\n\t		vmulpd		0x100(%%rbx),%%ymm12,%%ymm12		\n\t"\
		"vmulpd			 (%%rbx),%%ymm5,%%ymm5				\n\t		vmulpd		0x100(%%rbx),%%ymm13,%%ymm13		\n\t"\
	"vfnmadd231pd	0x020(%%rbx),%%ymm7,%%ymm4				\n\t	vfnmadd231pd	0x120(%%rbx),%%ymm15,%%ymm12	\n\t"\
	" vfmadd231pd	0x020(%%rbx),%%ymm6,%%ymm5				\n\t	 vfmadd231pd	0x120(%%rbx),%%ymm14,%%ymm13	\n\t"\
		"vmovaps		%%ymm5,%%ymm7					\n\t		vmovaps		%%ymm13,%%ymm15				\n\t"\
		"vmovaps		%%ymm4,%%ymm6					\n\t		vmovaps		%%ymm12,%%ymm14				\n\t"\
	"vmovaps	%%ymm11,(%%rax)	\n\t	vmovaps	(%%r8),%%ymm11	\n\t"/* spill to make room for 1.0 */\
		"vfnmadd231pd		 (%%rdx),%%ymm11,%%ymm4			\n\t		vfnmadd231pd	0x100(%%rdx),%%ymm11,%%ymm12	\n\t"\
		"vfnmadd231pd	0x020(%%rdx),%%ymm11,%%ymm5			\n\t		vfnmadd231pd	0x120(%%rdx),%%ymm11,%%ymm13	\n\t"\
		" vfmadd231pd		 (%%rdx),%%ymm11,%%ymm6			\n\t		 vfmadd231pd	0x100(%%rdx),%%ymm11,%%ymm14	\n\t"\
		" vfmadd231pd	0x020(%%rdx),%%ymm11,%%ymm7			\n\t		 vfmadd231pd	0x120(%%rdx),%%ymm11,%%ymm15	\n\t"\
		"vfmsub132pd	%%ymm11,%%ymm6,%%ymm0				\n\t		vfmsub132pd	%%ymm11,%%ymm14,%%ymm8		\n\t"\
		"vfmsub132pd	%%ymm11,%%ymm7,%%ymm1				\n\t		vfmsub132pd	%%ymm11,%%ymm15,%%ymm9		\n\t"\
		"vfmsub132pd	%%ymm11,%%ymm5,%%ymm2				\n\t		vfmsub132pd	%%ymm11,%%ymm13,%%ymm10	\n\t"\
		"vfmsub132pd	%%ymm11,%%ymm4,%%ymm3				\n\t		vfmsub132pd	(%%rax),%%ymm12,%%ymm11	\n\t"\
	"vmovaps	%%ymm12,(%%rdx) 	\n\t"/* spill ymm12 to make room for two */"	vmovaps	(%%rdi),%%ymm12 	\n\t"/* two */\
	"vfmadd132pd	%%ymm12,%%ymm0,%%ymm6					\n\t	vfmadd132pd	%%ymm12,%%ymm8 ,%%ymm14				\n\t"\
	"vfmadd132pd	%%ymm12,%%ymm1,%%ymm7					\n\t	vfmadd132pd	%%ymm12,%%ymm9 ,%%ymm15				\n\t"\
	"vfmadd132pd	%%ymm12,%%ymm2,%%ymm5					\n\t	vfmadd132pd	%%ymm12,%%ymm10,%%ymm13				\n\t"\
	"vfmadd132pd	%%ymm12,%%ymm3,%%ymm4					\n\t	vfmadd132pd	(%%rdx),%%ymm11,%%ymm12				\n\t"\
		"																vmovaps		%%ymm14,0x100(%%rdx)			\n\t"\
		"																vmovaps		%%ymm15,0x120(%%rdx)			\n\t"\
		"																vmovaps		%%ymm10,%%ymm14					\n\t"\
		"																vmovaps		%%ymm13,%%ymm15					\n\t"\
		"																vfmsub132pd	(%%r8),%%ymm12,%%ymm10	\n\t"\
		"																vfmsub132pd	(%%r8),%%ymm11,%%ymm13	\n\t"\
		"																vfmadd132pd	(%%r8),%%ymm12,%%ymm14	\n\t"\
		"																vfmadd132pd	(%%r8),%%ymm11,%%ymm15	\n\t"\
		/*****	SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r30) *****/\
		"																vmovaps		0x100(%%rdx),%%ymm11		\n\t"\
		"																vmovaps		0x120(%%rdx),%%ymm12		\n\t"\
	"vmovaps	%%ymm13,(%%rdx) 	\n\t"/* spill ymm13 to make room for isrt2 */"	vmovaps	(%%rsi),%%ymm13 	\n\t"/*isrt2*/\
		"vfmsub132pd	(%%r8),%%ymm11,%%ymm6				\n\t	vfnmadd231pd	%%ymm13,%%ymm10,%%ymm2			\n\t"\
		"vfmsub132pd	(%%r8),%%ymm9,%%ymm0				\n\t	vfnmadd231pd	%%ymm13,%%ymm15,%%ymm5			\n\t"\
		"vfmsub132pd	(%%r8),%%ymm12,%%ymm7				\n\t	vfnmadd231pd	%%ymm13,%%ymm14,%%ymm4			\n\t"\
		"vfmsub132pd	(%%r8),%%ymm8 ,%%ymm1				\n\t	vfnmadd231pd	(%%rdx),%%ymm13,%%ymm3			\n\t"\
		"vmovaps		%%ymm6,0x100(%%rdx)			\n\t		vmovaps		%%ymm2,0x140(%%rdx)	\n\t"\
		"vmovaps		%%ymm0,0x080(%%rdx)			\n\t		vmovaps		%%ymm5,0x0c0(%%rdx)	\n\t"\
		"vmovaps		%%ymm7,0x120(%%rdx)			\n\t		vmovaps		%%ymm4,0x160(%%rdx)	\n\t"\
		"vmovaps		%%ymm1,0x1a0(%%rdx)			\n\t		vmovaps		%%ymm3,0x1e0(%%rdx)	\n\t"\
	"vmovaps	0x40(%%rdi),%%ymm13 \n\t"/* sqrt2 *//* rcol: *2 => *sqrt2 due to above FMA leaving ymm10,13-15 unmultiplied-by-isrt2: */\
	"vfmadd132pd	(%%rdi),%%ymm6,%%ymm11				\n\t	 vfmadd132pd	%%ymm13,%%ymm2,%%ymm10	\n\t"\
	"vfmadd132pd	(%%rdi),%%ymm0,%%ymm9					\n\t	 vfmadd132pd	%%ymm13,%%ymm5,%%ymm15	\n\t"\
	"vfmadd132pd	(%%rdi),%%ymm7,%%ymm12				\n\t	 vfmadd132pd	%%ymm13,%%ymm4,%%ymm14	\n\t"\
	"vfmadd132pd	(%%rdi),%%ymm1,%%ymm8					\n\t	 vfmadd132pd	(%%rdx),%%ymm3,%%ymm13	\n\t"\
		"vmovaps		%%ymm11,     (%%rdx)			\n\t		vmovaps		%%ymm10,0x040(%%rdx)	\n\t"\
		"vmovaps		%%ymm9,0x180(%%rdx)			\n\t		vmovaps		%%ymm15,0x1c0(%%rdx)	\n\t"\
		"vmovaps		%%ymm12,0x020(%%rdx)			\n\t		vmovaps		%%ymm14,0x060(%%rdx)	\n\t"\
		"vmovaps		%%ymm8,0x0a0(%%rdx)			\n\t		vmovaps		%%ymm13,0x0e0(%%rdx)	\n\t"\
		"\n\t"\
	/**********************************************************************************/\
	/*...and now do eight radix-4 transforms, including the internal twiddle factors: */\
	/**********************************************************************************/\
		"movq		%[__isrt2],%%rsi		\n\t"\
		/*...Block 1: t00,t10,t20,t30	*/							/*...Block 5: t08,t18,t28,t38	*/\
		"movq		%[__r00],%%rax					\n\t		leaq	0x100(%%rax),%%r10	\n\t"\
		"movq		%[__r10],%%rbx					\n\t		leaq	0x100(%%rbx),%%r11	\n\t"\
		"movq		%[__r20],%%rcx					\n\t		leaq	0x100(%%rcx),%%r12	\n\t"\
		"movq		%[__r30],%%rdx					\n\t		leaq	0x100(%%rdx),%%r13	\n\t"\
		"vmovaps			  (%%rax),%%ymm0				\n\t		vmovaps			(%%r10),%%ymm8			\n\t"\
		"vmovaps		 0x020(%%rax),%%ymm1				\n\t		vmovaps		0x20(%%r10),%%ymm9			\n\t"\
		"vmovaps			  (%%rbx),%%ymm2				\n\t		vmovaps			(%%r11),%%ymm10			\n\t"\
		"vmovaps		 0x020(%%rbx),%%ymm3				\n\t		vmovaps		0x20(%%r11),%%ymm11			\n\t"\
		"vmovaps			  (%%rcx),%%ymm4				\n\t		vmovaps			(%%r12),%%ymm12			\n\t"\
		"vmovaps		 0x020(%%rcx),%%ymm5				\n\t		vmovaps		0x20(%%r12),%%ymm13			\n\t"\
		"vmovaps			  (%%rdx),%%ymm6				\n\t		vmovaps			(%%r13),%%ymm14			\n\t"\
		"vmovaps		 0x020(%%rdx),%%ymm7				\n\t"	/*	vmovaps		0x20(%%r13),%%ymm15			\n\t"*/\
	"vmovaps	(%%r8),%%ymm15	\n\t"/* spill to make room for 1.0 */\
		"vfmsub132pd	%%ymm15,%%ymm2 ,%%ymm0 				\n\t		vfmsub132pd		%%ymm15,%%ymm11,%%ymm8 	\n\t"\
		"vfmsub132pd	%%ymm15,%%ymm3 ,%%ymm1 				\n\t		vfmsub132pd		%%ymm15,%%ymm10,%%ymm9 	\n\t"\
		"vfmsub132pd	%%ymm15,%%ymm6 ,%%ymm4 				\n\t		vfmsub132pd		%%ymm15,%%ymm13,%%ymm12	\n\t"\
		"vfmsub132pd	%%ymm15,%%ymm7 ,%%ymm5 				\n\t		vfmsub132pd	0x20(%%r13),%%ymm14,%%ymm15	\n\t"\
	"vmovaps	%%ymm14,(%%rcx) 	\n\t"/* spill ymm14 to make room for two */"	vmovaps	(%%rdi),%%ymm14 	\n\t"/* two */\
	"vfmadd132pd		%%ymm14,%%ymm0,%%ymm2 				\n\t	 vfmadd132pd		%%ymm14,%%ymm8,%%ymm11	\n\t"\
	"vfmadd132pd		%%ymm14,%%ymm1,%%ymm3 				\n\t	 vfmadd132pd		%%ymm14,%%ymm9,%%ymm10	\n\t"\
	"vfmadd132pd		%%ymm14,%%ymm4,%%ymm6 				\n\t	 vfmadd132pd		%%ymm14,%%ymm12,%%ymm13	\n\t"\
	"vfmadd132pd		%%ymm14,%%ymm5,%%ymm7 				\n\t	 vfmadd132pd		(%%rcx),%%ymm15,%%ymm14	\n\t"\
		"vfmsub132pd	(%%r8),%%ymm6,%%ymm2				\n\t		vmulpd		(%%rsi),%%ymm12,%%ymm12		\n\t"/* isrt2 */\
		"vfmsub132pd	(%%r8),%%ymm7,%%ymm3				\n\t		vmulpd		(%%rsi),%%ymm13,%%ymm13		\n\t"\
	"vfmadd132pd	(%%rdi),%%ymm2,%%ymm6			\n\t	vfnmadd231pd	(%%rsi),%%ymm14,%%ymm12		\n\t"\
	"vfmadd132pd	(%%rdi),%%ymm3,%%ymm7			\n\t	vfnmadd231pd	(%%rsi),%%ymm15,%%ymm13		\n\t"\
		"vmovaps		%%ymm2,      (%%rcx)			\n\t	 vfmadd132pd	0x40(%%rdi),%%ymm12,%%ymm14		\n\t"/* sqrt2 */\
		"vmovaps		%%ymm3, 0x020(%%rcx)			\n\t	 vfmadd132pd	0x40(%%rdi),%%ymm13,%%ymm15		\n\t"\
		"													\n\t		vfmsub132pd	(%%r8),%%ymm12,%%ymm8			\n\t"\
		"													\n\t		vfmsub132pd	(%%r8),%%ymm13,%%ymm10		\n\t"\
		"vmovaps		%%ymm6,      (%%rax)			\n\t	 vfmadd132pd	(%%rdi),%%ymm8,%%ymm12		\n\t"\
		"vmovaps		%%ymm7, 0x020(%%rax)			\n\t	 vfmadd132pd	(%%rdi),%%ymm10,%%ymm13		\n\t"\
		"vfmsub132pd	(%%r8),%%ymm5,%%ymm0				\n\t		vmovaps		%%ymm8 ,     (%%r12)				\n\t"\
		"vfmsub132pd	(%%r8),%%ymm4,%%ymm1				\n\t		vmovaps		%%ymm10,0x020(%%r12)				\n\t"\
	"vfmadd132pd	(%%rdi),%%ymm0,%%ymm5			\n\t		vmovaps		%%ymm12,     (%%r10)				\n\t"\
	"vfmadd132pd	(%%rdi),%%ymm1,%%ymm4			\n\t		vmovaps		%%ymm13,0x020(%%r10)				\n\t"\
		"vmovaps		%%ymm0,      (%%rbx)			\n\t		vfmsub132pd	(%%r8),%%ymm15,%%ymm11		\n\t"\
		"vmovaps		%%ymm1, 0x020(%%rdx)			\n\t		vfmsub132pd	(%%r8),%%ymm14,%%ymm9			\n\t"\
		"													\n\t	 vfmadd132pd	(%%rdi),%%ymm11,%%ymm15		\n\t"\
		"													\n\t	 vfmadd132pd	(%%rdi),%%ymm9 ,%%ymm14		\n\t"\
		"vmovaps		%%ymm5,      (%%rdx)			\n\t		vmovaps		%%ymm11,     (%%r11)				\n\t"\
		"vmovaps		%%ymm4, 0x020(%%rbx)			\n\t		vmovaps		%%ymm9 ,0x020(%%r13)				\n\t"\
		"																vmovaps		%%ymm15,     (%%r13)				\n\t"\
		"																vmovaps		%%ymm14,0x020(%%r11)				\n\t"\
		/*...Block 3: t04,t14,t24,t34	*/							/*...Block 7: t0C,t1C,t2C,t3C	*/\
		"addq		$0x080,%%rax						\n\t		addq		$0x080,%%r10					\n\t"\
		"addq		$0x080,%%rbx						\n\t		addq		$0x080,%%r11					\n\t"\
		"addq		$0x080,%%rcx						\n\t		addq		$0x080,%%r12					\n\t"\
		"addq		$0x080,%%rdx						\n\t		addq		$0x080,%%r13					\n\t"\
		"vmovaps			 (%%rcx),%%ymm0				\n\t		vmovaps			 (%%r12),%%ymm8 				\n\t"\
		"vmovaps			 (%%rdx),%%ymm2				\n\t		vmovaps			 (%%r13),%%ymm10				\n\t"\
		"vmovaps		0x020(%%rcx),%%ymm1				\n\t		vmovaps		0x020(%%r12),%%ymm9					\n\t"\
		"vmovaps		0x020(%%rdx),%%ymm3				\n\t		vmovaps		0x020(%%r13),%%ymm11				\n\t"\
		"vmovaps			  %%ymm0,%%ymm4				\n\t		vmovaps			 %%ymm8 ,%%ymm12				\n\t"\
		"vmovaps			  %%ymm2,%%ymm6				\n\t		vmovaps			 %%ymm10,%%ymm14				\n\t"\
		"vmovaps			  %%ymm1,%%ymm5				\n\t		vmovaps			 %%ymm9,%%ymm13				\n\t"\
	/*	"vmovaps			  %%ymm3,%%ymm7				\n\t		vmovaps			 %%ymm11,%%ymm15				\n\t"*/\
	/* Strategic register-usage for sincos consts here cuts loads from 16 to 4: */\
		"vmovaps		0x040(%%rsi),%%ymm7 	\n\t"/* s */"		vmovaps		0x020(%%rsi),%%ymm15	\n\t"/* c */\
		"vmulpd		%%ymm7 ,%%ymm1,%%ymm1					\n\t		vmulpd		 %%ymm15,%%ymm9,%%ymm9				\n\t"\
		"vmulpd		%%ymm15,%%ymm3,%%ymm3					\n\t		vmulpd		 %%ymm7 ,%%ymm11,%%ymm11			\n\t"\
		"vmulpd		%%ymm7 ,%%ymm0,%%ymm0					\n\t		vmulpd		 %%ymm15,%%ymm8,%%ymm8				\n\t"\
		"vmulpd		%%ymm15,%%ymm2,%%ymm2					\n\t		vmulpd		 %%ymm7 ,%%ymm10,%%ymm10			\n\t"\
	"vfmsub132pd	%%ymm15,%%ymm1,%%ymm4					\n\t	vfmsub132pd		 %%ymm7 ,%%ymm9 ,%%ymm12			\n\t"\
	"vfmsub132pd	%%ymm7 ,%%ymm3,%%ymm6					\n\t	vfmsub132pd		 %%ymm15,%%ymm11,%%ymm14			\n\t"\
	"vfmadd132pd	%%ymm15,%%ymm0,%%ymm5					\n\t	vfmadd132pd		 %%ymm7 ,%%ymm8 ,%%ymm13			\n\t"\
	"vfmadd132pd 0x20(%%rdx),%%ymm2,%%ymm7					\n\t	vfmadd132pd 0x20(%%r13) ,%%ymm10,%%ymm15			\n\t"\
		"vfmsub132pd	(%%r8),%%ymm6,%%ymm4				\n\t		vfmsub132pd	(%%r8),%%ymm14,%%ymm12		\n\t"\
		"vfmsub132pd	(%%r8),%%ymm7,%%ymm5				\n\t		vfmsub132pd	(%%r8),%%ymm15,%%ymm13		\n\t"\
		"vmovaps			  (%%rbx),%%ymm2				\n\t		vmovaps			 (%%r11),%%ymm10				\n\t"\
		"vmovaps		 0x020(%%rbx),%%ymm3				\n\t		vmovaps		0x020(%%r11),%%ymm11				\n\t"\
	"vfmadd132pd	(%%rdi),%%ymm4,%%ymm6				\n\t	vfmadd132pd		(%%rdi),%%ymm12,%%ymm14		\n\t"\
	"vfmadd132pd	(%%rdi),%%ymm5,%%ymm7				\n\t	vfmadd132pd		(%%rdi),%%ymm13,%%ymm15		\n\t"\
		"vsubpd		0x020(%%rbx),%%ymm2,%%ymm2				\n\t		vaddpd		0x020(%%r11),%%ymm10,%%ymm10	\n\t"\
		"vaddpd			 (%%rbx),%%ymm3,%%ymm3				\n\t		vsubpd			 (%%r11),%%ymm11,%%ymm11	\n\t"\
	"vmovaps	%%ymm6 ,(%%rcx) 	\n\t"/* spill ymm6  to make room for isrt2 */"	vmovaps	(%%rsi),%%ymm6  	\n\t"/* isrt2 */\
		"vmovaps			  (%%rax),%%ymm0				\n\t		vmovaps			  (%%r10),%%ymm8				\n\t"\
		"vmovaps		 0x020(%%rax),%%ymm1				\n\t		vmovaps		 0x020(%%r10),%%ymm9				\n\t"\
	"vfnmadd231pd	%%ymm6 ,%%ymm2,%%ymm0					\n\t	vfnmadd231pd	%%ymm6 ,%%ymm10,%%ymm8			\n\t"\
	"vfnmadd231pd	%%ymm6 ,%%ymm3,%%ymm1					\n\t	vfnmadd231pd	%%ymm6 ,%%ymm11,%%ymm9			\n\t"\
	" vfmadd213pd		 (%%rax),%%ymm6 ,%%ymm2				\n\t	 vfmadd213pd		 (%%r10),%%ymm6 ,%%ymm10	\n\t"\
	" vfmadd213pd	0x020(%%rax),%%ymm6 ,%%ymm3				\n\t	 vfmadd213pd	0x020(%%r10),%%ymm6 ,%%ymm11	\n\t"\
	"vmovaps	(%%rcx),%%ymm6 	\n\t"/* restore spill */\
		"vfmsub132pd	(%%r8),%%ymm6,%%ymm2				\n\t		vfmsub132pd	(%%r8),%%ymm12,%%ymm8			\n\t"\
		"vfmsub132pd	(%%r8),%%ymm7,%%ymm3				\n\t		vfmsub132pd	(%%r8),%%ymm13,%%ymm9			\n\t"\
	"vfmadd132pd	(%%rdi),%%ymm2,%%ymm6				\n\t	vfmadd132pd		(%%rdi),%%ymm8 ,%%ymm12		\n\t"\
	"vfmadd132pd	(%%rdi),%%ymm3,%%ymm7				\n\t	vfmadd132pd		(%%rdi),%%ymm9 ,%%ymm13		\n\t"\
		"vmovaps		%%ymm2,      (%%rcx)			\n\t		vmovaps		%%ymm8 ,     (%%r12)				\n\t"\
		"vmovaps		%%ymm3, 0x020(%%rcx)			\n\t		vmovaps		%%ymm9 ,0x020(%%r12)				\n\t"\
		"vmovaps		%%ymm6,      (%%rax)			\n\t		vmovaps		%%ymm12,     (%%r10)				\n\t"\
		"vmovaps		%%ymm7, 0x020(%%rax)			\n\t		vmovaps		%%ymm13,0x020(%%r10)				\n\t"\
		"vfmsub132pd	(%%r8),%%ymm5,%%ymm0				\n\t		vfmsub132pd	(%%r8),%%ymm15,%%ymm10		\n\t"\
		"vfmsub132pd	(%%r8),%%ymm4,%%ymm1				\n\t		vfmsub132pd	(%%r8),%%ymm14,%%ymm11		\n\t"\
	"vfmadd132pd	(%%rdi),%%ymm0,%%ymm5				\n\t	vfmadd132pd		(%%rdi),%%ymm10,%%ymm15		\n\t"\
	"vfmadd132pd	(%%rdi),%%ymm1,%%ymm4				\n\t	vfmadd132pd		(%%rdi),%%ymm11,%%ymm14		\n\t"\
		"vmovaps		%%ymm0,      (%%rbx)			\n\t		vmovaps		%%ymm10,     (%%r11)				\n\t"\
		"vmovaps		%%ymm1, 0x020(%%rdx)			\n\t		vmovaps		%%ymm11,0x020(%%r13)				\n\t"\
		"vmovaps	%%ymm5,      (%%rdx)				\n\t		vmovaps		%%ymm15,     (%%r13)		\n\t"\
		"vmovaps	%%ymm4, 0x020(%%rbx)				\n\t		vmovaps		%%ymm14,0x020(%%r11)		\n\t"\
		/*...Block 2: t02,t12,t22,t32	*/				"	\n\t"	/*...Block 6: t0A,t1A,t2A,t3A */\
		"subq		$0x040,%%rax						\n\t		subq		$0x040,%%r10			\n\t"\
		"subq		$0x040,%%rbx						\n\t		subq		$0x040,%%r11			\n\t"\
		"subq		$0x040,%%rcx						\n\t		subq		$0x040,%%r12			\n\t"\
		"subq		$0x040,%%rdx						\n\t		subq		$0x040,%%r13			\n\t"\
		"addq		$0x060,%%rsi	\n\t"/* cc1 */\
		"vmovaps			 (%%rcx),%%ymm0				\n\t		vmovaps			 (%%r12),%%ymm8 				\n\t"\
		"vmovaps			 (%%rdx),%%ymm2				\n\t		vmovaps			 (%%r13),%%ymm14				\n\t"\
		"vmovaps		0x020(%%rcx),%%ymm1				\n\t		vmovaps		0x020(%%r12),%%ymm9					\n\t"\
		"vmovaps		0x020(%%rdx),%%ymm3				\n\t		vmovaps		0x020(%%r13),%%ymm11				\n\t"\
		"vmovaps			  %%ymm0,%%ymm4				\n\t		vmovaps			 %%ymm8 ,%%ymm12				\n\t"\
		"vmovaps			  %%ymm2,%%ymm6				\n\t"	/*	vmovaps			 %%ymm10,%%ymm14				\n\t"*/\
		"vmovaps			  %%ymm1,%%ymm5				\n\t"	/*	vmovaps			 %%ymm9,%%ymm13				\n\t"*/\
	/*	"vmovaps			  %%ymm3,%%ymm7				\n\t		vmovaps			 %%ymm11,%%ymm15				\n\t"*/\
	/* Strategic register-usage for sincos consts here cuts loads from 16 to 4: */\
		"vmovaps			 (%%rsi),%%ymm15			\n\t		vmovaps		0x020(%%rsi),%%ymm10			\n\t"\
		"vmovaps		0x040(%%rsi),%%ymm7 			\n\t		vmovaps		0x060(%%rsi),%%ymm13			\n\t"\
		"vmulpd			%%ymm10,%%ymm1,%%ymm1				\n\t		vmulpd			%%ymm7 ,%%ymm9 ,%%ymm9		\n\t"\
		"vmulpd			%%ymm13,%%ymm3,%%ymm3				\n\t		vmulpd			%%ymm10,%%ymm11,%%ymm11		\n\t"\
		"vmulpd			%%ymm10,%%ymm0,%%ymm0				\n\t		vmulpd			%%ymm7 ,%%ymm8 ,%%ymm8		\n\t"\
		"vmulpd			%%ymm13,%%ymm2,%%ymm2				\n\t		vmulpd			(%%r13),%%ymm10,%%ymm10		\n\t"\
	"vfmsub132pd		%%ymm15,%%ymm1,%%ymm4				\n\t	vfmsub132pd			%%ymm13,%%ymm9 ,%%ymm12		\n\t"\
	"vfmsub132pd		%%ymm7 ,%%ymm3,%%ymm6				\n\t	vfmadd132pd			%%ymm15,%%ymm11,%%ymm14		\n\t"\
	"vfmadd132pd		%%ymm15,%%ymm0,%%ymm5				\n\t	vfmadd132pd		0x020(%%r12),%%ymm8 ,%%ymm13	\n\t"\
	"vfmadd132pd	0x020(%%rdx),%%ymm2,%%ymm7				\n\t	vfmsub132pd		0x020(%%r13),%%ymm10,%%ymm15	\n\t"\
		"vfmsub132pd	(%%r8),%%ymm6,%%ymm4				\n\t		vfmsub132pd	(%%r8),%%ymm14,%%ymm12		\n\t"\
		"vfmsub132pd	(%%r8),%%ymm7,%%ymm5				\n\t		vfmsub132pd	(%%r8),%%ymm15,%%ymm13		\n\t"\
		"vmovaps			  (%%rbx),%%ymm2				\n\t		vmovaps			 (%%r11),%%ymm10				\n\t"\
		"vmovaps		 0x020(%%rbx),%%ymm0				\n\t		vmovaps		0x020(%%r11),%%ymm8					\n\t"\
	"vfmadd132pd	(%%rdi),%%ymm4,%%ymm6				\n\t	vfmadd132pd		(%%rdi),%%ymm12,%%ymm14		\n\t"\
	"vfmadd132pd	(%%rdi),%%ymm5,%%ymm7				\n\t	vfmadd132pd		(%%rdi),%%ymm13,%%ymm15		\n\t"\
		"vmovaps			  (%%rbx),%%ymm1				\n\t		vmovaps			 (%%r11),%%ymm9					\n\t"\
		"vmovaps		-0x040(%%rsi),%%ymm3				\n\t		vmovaps		-0x20(%%rsi),%%ymm11				\n\t"\
		"vmulpd			%%ymm11,%%ymm0,%%ymm0				\n\t		vmulpd			%%ymm3 ,%%ymm8 ,%%ymm8		\n\t"\
		"vmulpd			%%ymm11,%%ymm1,%%ymm1				\n\t		vmulpd			%%ymm3 ,%%ymm9 ,%%ymm9		\n\t"\
	"vfmsub132pd		%%ymm3 ,%%ymm0,%%ymm2				\n\t	vfmadd132pd			%%ymm11,%%ymm8 ,%%ymm10		\n\t"\
	"vfmadd132pd	 0x20(%%rbx),%%ymm1,%%ymm3				\n\t	vfmsub132pd		0x20(%%r11),%%ymm9 ,%%ymm11		\n\t"\
		"vmovaps			  (%%rax),%%ymm0				\n\t		vmovaps			 (%%r10),%%ymm8					\n\t"\
		"vmovaps		 0x020(%%rax),%%ymm1				\n\t		vmovaps		0x020(%%r10),%%ymm9					\n\t"\
		"vfmsub132pd	(%%r8),%%ymm2,%%ymm0				\n\t		vfmsub132pd	(%%r8),%%ymm10,%%ymm8			\n\t"\
		"vfmsub132pd	(%%r8),%%ymm3,%%ymm1				\n\t		vfmsub132pd	(%%r8),%%ymm11,%%ymm9			\n\t"\
		"vaddpd			  (%%rax),%%ymm2,%%ymm2			\n\t		vaddpd			 (%%r10),%%ymm10,%%ymm10	\n\t"\
		"vaddpd		 0x020(%%rax),%%ymm3,%%ymm3			\n\t		vaddpd		0x020(%%r10),%%ymm11,%%ymm11	\n\t"\
		"vfmsub132pd	(%%r8),%%ymm6,%%ymm2				\n\t		vfmsub132pd	(%%r8),%%ymm12,%%ymm8			\n\t"\
		"vfmsub132pd	(%%r8),%%ymm7,%%ymm3				\n\t		vfmsub132pd	(%%r8),%%ymm13,%%ymm9			\n\t"\
	"vfmadd132pd	(%%rdi),%%ymm2,%%ymm6			\n\t	vfmadd132pd		(%%rdi),%%ymm8,%%ymm12		\n\t"\
	"vfmadd132pd	(%%rdi),%%ymm3,%%ymm7			\n\t	vfmadd132pd		(%%rdi),%%ymm9,%%ymm13		\n\t"\
		"vmovaps		%%ymm2,      (%%rcx)			\n\t		vmovaps		%%ymm8 ,     (%%r12)				\n\t"\
		"vmovaps		%%ymm3, 0x020(%%rcx)			\n\t		vmovaps		%%ymm9 ,0x020(%%r12)				\n\t"\
		"vmovaps		%%ymm6,      (%%rax)			\n\t		vmovaps		%%ymm12,     (%%r10)				\n\t"\
		"vmovaps		%%ymm7, 0x020(%%rax)			\n\t		vmovaps		%%ymm13,0x020(%%r10)				\n\t"\
		"vfmsub132pd	(%%r8),%%ymm5,%%ymm0				\n\t		vfmsub132pd	(%%r8),%%ymm15,%%ymm10		\n\t"\
		"vfmsub132pd	(%%r8),%%ymm4,%%ymm1				\n\t		vfmsub132pd	(%%r8),%%ymm14,%%ymm11		\n\t"\
	"vfmadd132pd	(%%rdi),%%ymm0,%%ymm5			\n\t	vfmadd132pd		(%%rdi),%%ymm10,%%ymm15		\n\t"\
	"vfmadd132pd	(%%rdi),%%ymm1,%%ymm4			\n\t	vfmadd132pd		(%%rdi),%%ymm11,%%ymm14		\n\t"\
		"vmovaps		%%ymm0,      (%%rbx)			\n\t		vmovaps		%%ymm10,     (%%r11)				\n\t"\
		"vmovaps		%%ymm1, 0x020(%%rdx)			\n\t		vmovaps		%%ymm11,0x020(%%r13)				\n\t"\
		"vmovaps		%%ymm5,      (%%rdx)			\n\t		vmovaps		%%ymm15,     (%%r13)		\n\t"\
		"vmovaps		%%ymm4, 0x020(%%rbx)			\n\t		vmovaps		%%ymm14,0x020(%%r11)		\n\t"\
		/*...Block 4: t06,t16,t26,t36	*/							/*...Block 8: t0E,t1E,t2E,t3E	*/\
		"addq		$0x080,%%rax						\n\t		addq		$0x080,%%r10			\n\t"\
		"addq		$0x080,%%rbx						\n\t		addq		$0x080,%%r11			\n\t"\
		"addq		$0x080,%%rcx						\n\t		addq		$0x080,%%r12			\n\t"\
		"addq		$0x080,%%rdx						\n\t		addq		$0x080,%%r13			\n\t"\
		"vmovaps			 (%%rcx),%%ymm0				\n\t		vmovaps			 (%%r12),%%ymm8 				\n\t"\
		"vmovaps			 (%%rdx),%%ymm6				\n\t		vmovaps			 (%%r13),%%ymm10				\n\t"\
		"vmovaps		0x020(%%rcx),%%ymm1				\n\t		vmovaps		0x020(%%r12),%%ymm9					\n\t"\
		"vmovaps		0x020(%%rdx),%%ymm3				\n\t		vmovaps		0x020(%%r13),%%ymm11				\n\t"\
		"vmovaps			  %%ymm0,%%ymm4				\n\t		vmovaps			 %%ymm8 ,%%ymm12				\n\t"\
	/*	"vmovaps			  %%ymm2,%%ymm6			*/"	\n\t		vmovaps			 %%ymm10,%%ymm14				\n\t"\
	/*	"vmovaps			  %%ymm1,%%ymm5			*/"	\n\t		vmovaps			 %%ymm9,%%ymm13				\n\t"\
	/*	"vmovaps			  %%ymm3,%%ymm7				\n\t		vmovaps			 %%ymm11,%%ymm15				\n\t"*/\
	/* Strategic register-usage for sincos consts here cuts loads from 16 to 4: */\
		"vmovaps			 (%%rsi),%%ymm2 			\n\t		vmovaps		0x020(%%rsi),%%ymm7 			\n\t"\
		"vmovaps		0x040(%%rsi),%%ymm5 			\n\t		vmovaps		0x060(%%rsi),%%ymm15			\n\t"\
		"vmulpd		 	%%ymm15,%%ymm1,%%ymm1				\n\t		vmulpd			%%ymm2 ,%%ymm9,%%ymm9		\n\t"\
		"vmulpd		 	%%ymm2 ,%%ymm3,%%ymm3				\n\t		vmulpd			%%ymm5 ,%%ymm11,%%ymm11		\n\t"\
		"vmulpd		 	%%ymm15,%%ymm0,%%ymm0				\n\t		vmulpd			%%ymm2 ,%%ymm8,%%ymm8		\n\t"\
		"vmulpd			 (%%rdx),%%ymm2,%%ymm2				\n\t		vmulpd			%%ymm5 ,%%ymm10,%%ymm10		\n\t"\
	"vfmsub132pd	 	%%ymm5 ,%%ymm1,%%ymm4				\n\t	vfmsub132pd			%%ymm7 ,%%ymm9 ,%%ymm12		\n\t"\
	"vfmadd132pd	 	%%ymm7 ,%%ymm3,%%ymm6				\n\t	vfmsub132pd			%%ymm15,%%ymm11,%%ymm14		\n\t"\
	"vfmadd132pd	0x020(%%rcx),%%ymm0,%%ymm5				\n\t	vfmadd132pd			%%ymm7 ,%%ymm8 ,%%ymm13		\n\t"\
	"vfmsub132pd	0x020(%%rdx),%%ymm2,%%ymm7				\n\t	vfmadd132pd		0x20(%%r13),%%ymm10,%%ymm15		\n\t"\
		"vfmsub132pd	(%%r8),%%ymm6,%%ymm4				\n\t		vfmsub132pd	(%%r8),%%ymm14,%%ymm12		\n\t"\
		"vfmsub132pd	(%%r8),%%ymm7,%%ymm5				\n\t		vfmsub132pd	(%%r8),%%ymm15,%%ymm13		\n\t"\
		"vmovaps			  (%%rbx),%%ymm2				\n\t		vmovaps			 (%%r11),%%ymm10				\n\t"\
		"vmovaps		 0x020(%%rbx),%%ymm0				\n\t		vmovaps		0x020(%%r11),%%ymm8					\n\t"\
	"vfmadd132pd	(%%rdi),%%ymm4,%%ymm6				\n\t	vfmadd132pd		(%%rdi),%%ymm12,%%ymm14		\n\t"\
	"vfmadd132pd	(%%rdi),%%ymm5,%%ymm7				\n\t	vfmadd132pd		(%%rdi),%%ymm13,%%ymm15		\n\t"\
		"vmovaps			  (%%rbx),%%ymm1				\n\t		vmovaps			 (%%r11),%%ymm9					\n\t"\
		"vmovaps		-0x020(%%rsi),%%ymm3				\n\t		vmovaps		-0x40(%%rsi),%%ymm11				\n\t"\
		"vmulpd			%%ymm11,%%ymm0,%%ymm0				\n\t		vmulpd			%%ymm3 ,%%ymm8 ,%%ymm8		\n\t"\
		"vmulpd			%%ymm11,%%ymm1,%%ymm1				\n\t		vmulpd			%%ymm3 ,%%ymm9 ,%%ymm9		\n\t"\
	"vfmsub132pd		%%ymm3 ,%%ymm0,%%ymm2				\n\t	vfmadd132pd			%%ymm11,%%ymm8 ,%%ymm10		\n\t"\
	"vfmadd132pd	 0x20(%%rbx),%%ymm1,%%ymm3				\n\t	vfmsub132pd		0x20(%%r11),%%ymm9 ,%%ymm11		\n\t"\
		"vmovaps			  (%%rax),%%ymm0				\n\t		vmovaps			 (%%r10),%%ymm8					\n\t"\
		"vmovaps		 0x020(%%rax),%%ymm1				\n\t		vmovaps		0x020(%%r10),%%ymm9					\n\t"\
		"vfmsub132pd	(%%r8),%%ymm2,%%ymm0				\n\t		vfmsub132pd	(%%r8),%%ymm10,%%ymm8			\n\t"\
		"vfmsub132pd	(%%r8),%%ymm3,%%ymm1				\n\t		vfmsub132pd	(%%r8),%%ymm11,%%ymm9			\n\t"\
		"vaddpd			  (%%rax),%%ymm2,%%ymm2			\n\t		vaddpd			 (%%r10),%%ymm10,%%ymm10	\n\t"\
		"vaddpd		 0x020(%%rax),%%ymm3,%%ymm3			\n\t		vaddpd		0x020(%%r10),%%ymm11,%%ymm11	\n\t"\
		"vfmsub132pd	(%%r8),%%ymm4,%%ymm2				\n\t		vfmsub132pd	(%%r8),%%ymm12,%%ymm8			\n\t"\
		"vfmsub132pd	(%%r8),%%ymm5,%%ymm3				\n\t		vfmsub132pd	(%%r8),%%ymm13,%%ymm9			\n\t"\
	"vfmadd132pd	(%%rdi),%%ymm2,%%ymm4			\n\t	vfmadd132pd		(%%rdi),%%ymm8,%%ymm12		\n\t"\
	"vfmadd132pd	(%%rdi),%%ymm3,%%ymm5			\n\t	vfmadd132pd		(%%rdi),%%ymm9,%%ymm13		\n\t"\
		"vmovaps		%%ymm2,      (%%rcx)			\n\t		vmovaps		%%ymm8 ,     (%%r12)				\n\t"\
		"vmovaps		%%ymm3, 0x020(%%rcx)			\n\t		vmovaps		%%ymm9 ,0x020(%%r12)				\n\t"\
		"vmovaps		%%ymm4,      (%%rax)			\n\t		vmovaps		%%ymm12,     (%%r10)				\n\t"\
		"vmovaps		%%ymm5, 0x020(%%rax)			\n\t		vmovaps		%%ymm13,0x020(%%r10)				\n\t"\
		"vfmsub132pd	(%%r8),%%ymm7,%%ymm0				\n\t		vfmsub132pd	(%%r8),%%ymm15,%%ymm10		\n\t"\
		"vfmsub132pd	(%%r8),%%ymm6,%%ymm1				\n\t		vfmsub132pd	(%%r8),%%ymm14,%%ymm11		\n\t"\
	"vfmadd132pd	(%%rdi),%%ymm0,%%ymm7			\n\t	vfmadd132pd		(%%rdi),%%ymm10,%%ymm15		\n\t"\
	"vfmadd132pd	(%%rdi),%%ymm1,%%ymm6			\n\t	vfmadd132pd		(%%rdi),%%ymm11,%%ymm14		\n\t"\
		"vmovaps		%%ymm0,      (%%rbx)			\n\t		vmovaps		%%ymm10,     (%%r11)				\n\t"\
		"vmovaps		%%ymm1, 0x020(%%rdx)			\n\t		vmovaps		%%ymm11,0x020(%%r13)				\n\t"\
		"vmovaps		%%ymm7,      (%%rdx)			\n\t		vmovaps		%%ymm15,     (%%r13)		\n\t"\
		"vmovaps		%%ymm6, 0x020(%%rbx)			\n\t		vmovaps		%%ymm14,0x020(%%r11)		\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__add1] "m" (Xadd1)\
		 ,[__add2] "m" (Xadd2)\
		 ,[__add3] "m" (Xadd3)\
		 ,[__r00] "m" (Xr00)\
		 ,[__r10] "m" (Xr10)\
		 ,[__r20] "m" (Xr20)\
		 ,[__r30] "m" (Xr30)\
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__cc0] "m" (Xcc0)\
		 ,[__c00] "m" (Xc00)\
		 ,[__c01] "m" (Xc01)\
		 ,[__c02] "m" (Xc02)\
		 ,[__c03] "m" (Xc03)\
		 ,[__c05] "m" (Xc05)\
		 ,[__c07] "m" (Xc07)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r8","r10","r11","r12","r13","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"		/* Clobbered registers */\
	);\
	}

	// Aggressive-FMA: replace [70 ADD, 170 SUB, 116 MUL, 202 FMA, 996 memref] ==> [6 ADD, 6 SUB, 426 FMA (202 nontrivial), 116 MUL, 996 memref].
	//
	#define SSE2_RADIX32_WRAPPER_DIT(Xadd0,Xadd1,Xadd2,Xadd3,Xisrt2,Xr00,Xr08,Xr10,Xr20,Xr28,Xr30,Xc01,Xc02,Xc04,Xc06,Xc08,Xc0A,Xc0C,Xc0E,Xc10,Xc12,Xc14,Xc16,Xc18,Xc1A,Xc1C,Xc1E)\
	{\
	__asm__ volatile (\
	/************************************************************************/\
	/*...And do an inverse DIT radix-32 pass on the squared-data blocks:	*/\
	/************************************************************************/\
	/*...Block 1: */\
		"movq	%[__isrt2],%%rsi	\n\t	leaq	0x940(%%rsi),%%r8	/* one */\n\t	leaq	0x8e0(%%rsi),%%rdi	/* two */\n\t"\
		"movq		%[__r00],%%rax			\n\t"\
		"leaq		0x400(%%rax),%%rbx				\n\t"\
		"leaq		0x200(%%rax),%%rcx				\n\t"\
		"leaq		0x600(%%rax),%%rdx				\n\t"\
		/*****	SSE2_RADIX4_DIT_IN_PLACE()	*****/				/*...Block 2 has tmp-addresses offset +0x80 w.r.to Block 1:	*/\
		"vmovaps		(%%rax),%%ymm0				\n\t		vmovaps		0x80(%%rax),%%ymm8 			\n\t"\
		"vmovaps	0x20(%%rax),%%ymm1				\n\t		vmovaps		0xa0(%%rax),%%ymm9 			\n\t"\
		"vmovaps		(%%rax),%%ymm2				\n\t		vmovaps		0x80(%%rax),%%ymm10			\n\t"\
		"vmovaps	0x20(%%rax),%%ymm3				\n\t		vmovaps		0xa0(%%rax),%%ymm11			\n\t"\
	"vmovaps	(%%r8),%%ymm15	\n\t"/* 1.0 */\
		" vfmadd231pd		(%%rbx),%%ymm15,%%ymm0	\n\t		 vfmadd231pd	0x80(%%rbx),%%ymm15,%%ymm8 	\n\t"\
		" vfmadd231pd	0x20(%%rbx),%%ymm15,%%ymm1	\n\t		 vfmadd231pd	0xa0(%%rbx),%%ymm15,%%ymm9 	\n\t"\
		"vfnmadd231pd		(%%rbx),%%ymm15,%%ymm2	\n\t		vfnmadd231pd	0x80(%%rbx),%%ymm15,%%ymm10	\n\t"\
		"vfnmadd231pd	0x20(%%rbx),%%ymm15,%%ymm3	\n\t		vfnmadd231pd	0xa0(%%rbx),%%ymm15,%%ymm11	\n\t"\
		"vmovaps		(%%rcx),%%ymm4				\n\t		vmovaps		0x80(%%rcx),%%ymm12			\n\t"\
		"vmovaps	0x20(%%rcx),%%ymm5				\n\t		vmovaps		0xa0(%%rcx),%%ymm13			\n\t"\
		"vmovaps		(%%rcx),%%ymm6				\n\t		vmovaps		0x80(%%rcx),%%ymm14			\n\t"\
		"vmovaps	0x20(%%rcx),%%ymm7				\n\t		vmovaps		0xa0(%%rcx),%%ymm15			\n\t"\
	"vmovaps	%%ymm11,(%%rbx)	\n\t	vmovaps	(%%r8),%%ymm11	\n\t"/* spill to make room for 1.0 */\
		" vfmadd231pd		(%%rdx),%%ymm11,%%ymm4	\n\t		 vfmadd231pd	0x80(%%rdx),%%ymm11,%%ymm12	\n\t"\
		" vfmadd231pd	0x20(%%rdx),%%ymm11,%%ymm5	\n\t		 vfmadd231pd	0xa0(%%rdx),%%ymm11,%%ymm13	\n\t"\
		"vfnmadd231pd		(%%rdx),%%ymm11,%%ymm6	\n\t		vfnmadd231pd	0x80(%%rdx),%%ymm11,%%ymm14	\n\t"\
		"vfnmadd231pd	0x20(%%rdx),%%ymm11,%%ymm7	\n\t		vfnmadd231pd	0xa0(%%rdx),%%ymm11,%%ymm15	\n\t"\
		"vfmsub132pd		%%ymm11,%%ymm4,%%ymm0	\n\t		vfmsub132pd		%%ymm11,%%ymm12,%%ymm8 	\n\t"\
		"vfmsub132pd		%%ymm11,%%ymm5,%%ymm1	\n\t		vfmsub132pd		%%ymm11,%%ymm13,%%ymm9 	\n\t"\
		"vfmsub132pd		%%ymm11,%%ymm7,%%ymm2	\n\t		vfmsub132pd		%%ymm11,%%ymm15,%%ymm10	\n\t"\
		"vfmsub132pd		%%ymm11,%%ymm6,%%ymm3	\n\t		vfmsub132pd		(%%rbx),%%ymm14,%%ymm11	\n\t"\
		"vmovaps		%%ymm0 ,      (%%rbx)		\n\t		vmovaps			%%ymm8 , 0x080(%%rbx)	\n\t"\
		"vmovaps		%%ymm1 , 0x020(%%rbx)		\n\t		vmovaps			%%ymm9 , 0x0a0(%%rbx)	\n\t"\
		"vmovaps		%%ymm2 ,      (%%rdx)		\n\t		vmovaps			%%ymm10, 0x080(%%rdx)	\n\t"\
		"vmovaps		%%ymm3 , 0x020(%%rcx)		\n\t		vmovaps			%%ymm11, 0x0a0(%%rcx)	\n\t"\
	"vmovaps	%%ymm14,0xa0(%%rdx)	\n\t"/* spill ymm14 to make room for two */"	vmovaps	(%%rdi),%%ymm14	\n\t"/* two */\
	"vfmadd132pd	%%ymm14,%%ymm0,%%ymm4			\n\t	vfmadd132pd		%%ymm14,%%ymm8 ,%%ymm12		\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm1,%%ymm5			\n\t	vfmadd132pd		%%ymm14,%%ymm9 ,%%ymm13		\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm2,%%ymm7			\n\t	vfmadd132pd		%%ymm14,%%ymm10,%%ymm15		\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm3,%%ymm6			\n\t	vfmadd132pd	0xa0(%%rdx),%%ymm11,%%ymm14		\n\t"\
		"vmovaps	%%ymm4 ,    (%%rax)				\n\t		vmovaps		%%ymm12,0x80(%%rax)			\n\t"\
		"vmovaps	%%ymm5 ,0x20(%%rax)				\n\t		vmovaps		%%ymm13,0xa0(%%rax)			\n\t"\
		"vmovaps	%%ymm7 ,    (%%rcx)				\n\t		vmovaps		%%ymm15,0x80(%%rcx)			\n\t"\
		"vmovaps	%%ymm6 ,0x20(%%rdx)				\n\t		vmovaps		%%ymm14,0xa0(%%rdx)			\n\t"\
		"addq		$0x100,%%rax				\n\t"\
		"addq		$0x100,%%rbx				\n\t"\
		"addq		$0x100,%%rcx				\n\t"\
		"addq		$0x100,%%rdx				\n\t"\
		"vmovaps		(%%rax),%%ymm0				\n\t		vmovaps		0x80(%%rax),%%ymm8 			\n\t"\
		"vmovaps	0x20(%%rax),%%ymm1				\n\t		vmovaps		0xa0(%%rax),%%ymm9 			\n\t"\
		"vmovaps		(%%rax),%%ymm2				\n\t		vmovaps		0x80(%%rax),%%ymm10			\n\t"\
		"vmovaps	0x20(%%rax),%%ymm3				\n\t		vmovaps		0xa0(%%rax),%%ymm11			\n\t"\
	"vmovaps	(%%r8),%%ymm15	\n\t"/* 1.0 */\
		" vfmadd231pd		(%%rbx),%%ymm15,%%ymm0	\n\t		 vfmadd231pd	0x80(%%rbx),%%ymm15,%%ymm8 	\n\t"\
		" vfmadd231pd	0x20(%%rbx),%%ymm15,%%ymm1	\n\t		 vfmadd231pd	0xa0(%%rbx),%%ymm15,%%ymm9 	\n\t"\
		"vfnmadd231pd		(%%rbx),%%ymm15,%%ymm2	\n\t		vfnmadd231pd	0x80(%%rbx),%%ymm15,%%ymm10	\n\t"\
		"vfnmadd231pd	0x20(%%rbx),%%ymm15,%%ymm3	\n\t		vfnmadd231pd	0xa0(%%rbx),%%ymm15,%%ymm11	\n\t"\
		"vmovaps		(%%rcx),%%ymm4				\n\t		vmovaps		0x80(%%rcx),%%ymm12			\n\t"\
		"vmovaps	0x20(%%rcx),%%ymm5				\n\t		vmovaps		0xa0(%%rcx),%%ymm13			\n\t"\
		"vmovaps		(%%rcx),%%ymm6				\n\t		vmovaps		0x80(%%rcx),%%ymm14			\n\t"\
		"vmovaps	0x20(%%rcx),%%ymm7				\n\t		vmovaps		0xa0(%%rcx),%%ymm15			\n\t"\
	"vmovaps	%%ymm11,(%%rbx)	\n\t	vmovaps	(%%r8),%%ymm11	\n\t"/* spill to make room for 1.0 */\
		" vfmadd231pd		(%%rdx),%%ymm11,%%ymm4	\n\t		 vfmadd231pd	0x80(%%rdx),%%ymm11,%%ymm12	\n\t"\
		" vfmadd231pd	0x20(%%rdx),%%ymm11,%%ymm5	\n\t		 vfmadd231pd	0xa0(%%rdx),%%ymm11,%%ymm13	\n\t"\
		"vfnmadd231pd		(%%rdx),%%ymm11,%%ymm6	\n\t		vfnmadd231pd	0x80(%%rdx),%%ymm11,%%ymm14	\n\t"\
		"vfnmadd231pd	0x20(%%rdx),%%ymm11,%%ymm7	\n\t		vfnmadd231pd	0xa0(%%rdx),%%ymm11,%%ymm15	\n\t"\
		"vfmsub132pd		%%ymm11,%%ymm4,%%ymm0		\n\t		vfmsub132pd		%%ymm11,%%ymm12,%%ymm8 	\n\t"\
		"vfmsub132pd		%%ymm11,%%ymm5,%%ymm1		\n\t		vfmsub132pd		%%ymm11,%%ymm13,%%ymm9 	\n\t"\
		"vfmsub132pd		%%ymm11,%%ymm7,%%ymm2		\n\t		vfmsub132pd		%%ymm11,%%ymm15,%%ymm10	\n\t"\
		"vfmsub132pd		%%ymm11,%%ymm6,%%ymm3		\n\t		vfmsub132pd		(%%rbx),%%ymm14,%%ymm11	\n\t"\
	"vmovaps	%%ymm14,0xa0(%%rdx)	\n\t"/* spill ymm14 to make room for two */"	vmovaps	(%%rdi),%%ymm14	\n\t"/* two */\
	"vfmadd132pd	%%ymm14,%%ymm0,%%ymm4			\n\t	vfmadd132pd		%%ymm14,%%ymm8 ,%%ymm12		\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm1,%%ymm5			\n\t	vfmadd132pd		%%ymm14,%%ymm9 ,%%ymm13		\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm2,%%ymm7			\n\t	vfmadd132pd		%%ymm14,%%ymm10,%%ymm15		\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm3,%%ymm6			\n\t	vfmadd132pd	0xa0(%%rdx),%%ymm11,%%ymm14		\n\t"\
		"vmovaps	%%ymm0 ,      (%%rbx)			\n\t		vmovaps			%%ymm8 , 0x080(%%rbx)	\n\t"\
		"vmovaps	%%ymm1 , 0x020(%%rbx)			\n\t		vmovaps			%%ymm9 , 0x0a0(%%rbx)	\n\t"\
		"vmovaps	%%ymm4 ,      (%%rax)			\n\t		vmovaps			%%ymm12, 0x080(%%rax)	\n\t"\
		"vmovaps	%%ymm5 , 0x020(%%rax)			\n\t		vmovaps			%%ymm13, 0x0a0(%%rax)	\n\t"\
		"vmovaps	%%ymm3 ,%%ymm0					\n\t		vmovaps			%%ymm11,%%ymm8 			\n\t"\
		"vmovaps	%%ymm6 ,%%ymm1					\n\t"	/*	vmovaps			%%ymm14,%%ymm9 			\n\t"*/\
	"vmovaps	%%ymm14,(%%rdx)	\n\t	vmovaps	(%%r8),%%ymm9	\n\t"/* spill to make room for 1.0 */\
		"vfmsub132pd	%%ymm9 ,%%ymm7,%%ymm3			\n\t		vfmsub132pd		%%ymm9 ,%%ymm15,%%ymm11	\n\t"\
		"vfmsub132pd	%%ymm9 ,%%ymm2,%%ymm6			\n\t		vfmsub132pd		%%ymm9 ,%%ymm10,%%ymm14	\n\t"\
		"vfmadd132pd	%%ymm9 ,%%ymm7,%%ymm0			\n\t		vfmadd132pd		%%ymm9 ,%%ymm15,%%ymm8 	\n\t"\
		"vfmadd132pd	%%ymm9 ,%%ymm2,%%ymm1			\n\t		vfmadd132pd		(%%rdx),%%ymm10,%%ymm9 	\n\t"\
	"vmovaps	(%%rsi),%%ymm15	\n\t"/* isrt2 */\
		"vmulpd		%%ymm15,%%ymm3,%%ymm3			\n\t		vmulpd			%%ymm15,%%ymm11,%%ymm11	\n\t"\
		"vmulpd		%%ymm15,%%ymm6,%%ymm6			\n\t		vmulpd			%%ymm15,%%ymm14,%%ymm14	\n\t"\
		"vmulpd		%%ymm15,%%ymm0,%%ymm0			\n\t		vmulpd			%%ymm15,%%ymm8 ,%%ymm8 	\n\t"\
		"vmulpd		%%ymm15,%%ymm1,%%ymm1			\n\t		vmulpd			%%ymm15,%%ymm9 ,%%ymm9 	\n\t"\
		"vmovaps	%%ymm3 , 0x020(%%rcx)			\n\t		vmovaps			%%ymm11, 0x0a0(%%rcx)	\n\t"\
		"vmovaps	%%ymm6 , 0x020(%%rdx)			\n\t		vmovaps			%%ymm14, 0x0a0(%%rdx)	\n\t"\
		"vmovaps	%%ymm0 ,      (%%rcx)			\n\t		vmovaps			%%ymm8 , 0x080(%%rcx)	\n\t"\
		"vmovaps	%%ymm1 ,      (%%rdx)			\n\t		vmovaps			%%ymm9 , 0x080(%%rdx)	\n\t"\
		/*****	SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS *****/			/*****	SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS *****/\
		/*        (r00,r10,r20,r30,r08,r18,r28,r38):  */			/*        (r04,r14,r24,r34,r0C,r1C,r2C,r3C):  */\
		"vmovaps	-0x100(%%rax),%%ymm0			\n\t		vmovaps		-0x080(%%rax),%%ymm8 		\n\t"\
		"vmovaps	-0x100(%%rbx),%%ymm4			\n\t		vmovaps		-0x080(%%rbx),%%ymm12		\n\t"\
		"vmovaps	-0x0e0(%%rax),%%ymm1			\n\t		vmovaps		-0x060(%%rax),%%ymm9 		\n\t"\
		"vmovaps	-0x0e0(%%rbx),%%ymm5			\n\t"	/*	vmovaps		-0x060(%%rbx),%%ymm13		\n\t"*/	"vmovaps	(%%r8),%%ymm13	\n\t"/* spill to make room for 1.0 */\
		"vmovaps		  (%%rax),%%ymm2			\n\t		vmovaps		 0x080(%%rax),%%ymm10		\n\t"\
		"vmovaps	 0x020(%%rbx),%%ymm7			\n\t		vmovaps		 0x0a0(%%rbx),%%ymm15		\n\t"\
		"vmovaps	 0x020(%%rax),%%ymm3			\n\t		vmovaps		 0x0a0(%%rax),%%ymm11		\n\t"\
		"vmovaps		  (%%rbx),%%ymm6			\n\t		vmovaps		 0x080(%%rbx),%%ymm14		\n\t"\
		"vfmsub132pd	%%ymm13,%%ymm2,%%ymm0			\n\t		vfmsub132pd		%%ymm13 ,%%ymm10,%%ymm8 		\n\t"\
		"vfmsub132pd	%%ymm13,%%ymm7,%%ymm4			\n\t		vfmsub132pd		%%ymm13 ,%%ymm15,%%ymm12		\n\t"\
		"vfmsub132pd	%%ymm13,%%ymm3,%%ymm1			\n\t		vfmsub132pd		%%ymm13 ,%%ymm11,%%ymm9 		\n\t"\
		"vfmsub132pd	%%ymm13,%%ymm6,%%ymm5			\n\t		vfmsub132pd	-0x60(%%rbx),%%ymm14,%%ymm13		\n\t"\
	"vmovaps	%%ymm14,0xa0(%%rbx)	\n\t"/* spill ymm14 to make room for two */"	vmovaps	(%%rdi),%%ymm14	\n\t"/* two */\
	"vfmadd132pd	%%ymm14,%%ymm0,%%ymm2			\n\t	vfmadd132pd		%%ymm14,%%ymm8 ,%%ymm10		\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm4,%%ymm7			\n\t	vfmadd132pd		%%ymm14,%%ymm12,%%ymm15		\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm1,%%ymm3			\n\t	vfmadd132pd		%%ymm14,%%ymm9 ,%%ymm11		\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm5,%%ymm6			\n\t	vfmadd132pd	0xa0(%%rbx),%%ymm13,%%ymm14		\n\t"\
		"vmovaps	%%ymm0,      (%%rax)			\n\t		vmovaps		%%ymm8 , 0x080(%%rax)		\n\t"\
		"vmovaps	%%ymm4,      (%%rbx)			\n\t		vmovaps		%%ymm12, 0x080(%%rbx)		\n\t"\
		"vmovaps	%%ymm1, 0x020(%%rax)			\n\t		vmovaps		%%ymm9 , 0x0a0(%%rax)		\n\t"\
		"vmovaps	%%ymm5,-0x0e0(%%rbx)			\n\t		vmovaps		%%ymm13,-0x060(%%rbx)		\n\t"\
		"vmovaps	%%ymm2,-0x100(%%rax)			\n\t		vmovaps		%%ymm10,-0x080(%%rax)		\n\t"\
		"vmovaps	%%ymm7,-0x100(%%rbx)			\n\t		vmovaps		%%ymm15,-0x080(%%rbx)		\n\t"\
		"vmovaps	%%ymm3,-0x0e0(%%rax)			\n\t		vmovaps		%%ymm11,-0x060(%%rax)		\n\t"\
		"vmovaps	%%ymm6, 0x020(%%rbx)			\n\t		vmovaps		%%ymm14, 0x0a0(%%rbx)		\n\t"\
		"vmovaps	-0x100(%%rcx),%%ymm0			\n\t		vmovaps		-0x080(%%rcx),%%ymm8 		\n\t"\
		"vmovaps	-0x100(%%rdx),%%ymm4			\n\t		vmovaps		-0x080(%%rdx),%%ymm12		\n\t"\
		"vmovaps	-0x0e0(%%rcx),%%ymm1			\n\t		vmovaps		-0x060(%%rcx),%%ymm9 		\n\t"\
		"vmovaps	-0x0e0(%%rdx),%%ymm5			\n\t"	/*	vmovaps		-0x060(%%rdx),%%ymm13		\n\t"*/	"vmovaps	(%%r8),%%ymm13	\n\t"/* spill to make room for 1.0 */\
		"vmovaps		  (%%rcx),%%ymm2			\n\t		vmovaps		 0x080(%%rcx),%%ymm10		\n\t"\
		"vmovaps	 0x020(%%rdx),%%ymm7			\n\t		vmovaps		 0x0a0(%%rdx),%%ymm15		\n\t"\
		"vmovaps	 0x020(%%rcx),%%ymm3			\n\t		vmovaps		 0x0a0(%%rcx),%%ymm11		\n\t"\
		"vmovaps		  (%%rdx),%%ymm6			\n\t		vmovaps		 0x080(%%rdx),%%ymm14		\n\t"\
		"vfmsub132pd	%%ymm13,%%ymm2,%%ymm0			\n\t		vfmsub132pd		%%ymm13 ,%%ymm10,%%ymm8 		\n\t"\
		"vfmsub132pd	%%ymm13,%%ymm7,%%ymm4			\n\t		vfmsub132pd		%%ymm13 ,%%ymm15,%%ymm12		\n\t"\
		"vfmsub132pd	%%ymm13,%%ymm3,%%ymm1			\n\t		vfmsub132pd		%%ymm13 ,%%ymm11,%%ymm9 		\n\t"\
		"vfmsub132pd	%%ymm13,%%ymm6,%%ymm5			\n\t		vfmsub132pd	-0x60(%%rdx),%%ymm14,%%ymm13		\n\t"\
	"vmovaps	%%ymm14,0xa0(%%rdx)	\n\t"/* spill ymm14 to make room for two */"	vmovaps	(%%rdi),%%ymm14	\n\t"/* two */\
	"vfmadd132pd	%%ymm14,%%ymm0,%%ymm2			\n\t	vfmadd132pd		%%ymm14,%%ymm8 ,%%ymm10		\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm4,%%ymm7			\n\t	vfmadd132pd		%%ymm14,%%ymm12,%%ymm15		\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm1,%%ymm3			\n\t	vfmadd132pd		%%ymm14,%%ymm9 ,%%ymm11		\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm5,%%ymm6			\n\t	vfmadd132pd	0xa0(%%rdx),%%ymm13,%%ymm14		\n\t"\
		"vmovaps	%%ymm0,      (%%rcx)			\n\t		vmovaps		%%ymm8 , 0x080(%%rcx)		\n\t"\
		"vmovaps	%%ymm4,      (%%rdx)			\n\t		vmovaps		%%ymm12, 0x080(%%rdx)		\n\t"\
		"vmovaps	%%ymm1, 0x020(%%rcx)			\n\t		vmovaps		%%ymm9 , 0x0a0(%%rcx)		\n\t"\
		"vmovaps	%%ymm5,-0x0e0(%%rdx)			\n\t		vmovaps		%%ymm13,-0x060(%%rdx)		\n\t"\
		"vmovaps	%%ymm2,-0x100(%%rcx)			\n\t		vmovaps		%%ymm10,-0x080(%%rcx)		\n\t"\
		"vmovaps	%%ymm7,-0x100(%%rdx)			\n\t		vmovaps		%%ymm15,-0x080(%%rdx)		\n\t"\
		"vmovaps	%%ymm3,-0x0e0(%%rcx)			\n\t		vmovaps		%%ymm11,-0x060(%%rcx)		\n\t"\
		"vmovaps	%%ymm6, 0x020(%%rdx)			\n\t		vmovaps		%%ymm14, 0x0a0(%%rdx)		\n\t"\
	/*...Blocks 3,4 have tmp-addresses offset +0x40 w.r.to Blocks 1,2, respectively (thus +0x100-0x0c0 = +0x040: */\
		"subq		$0xc0,%%rax				\n\t"\
		"subq		$0xc0,%%rbx				\n\t"\
		"subq		$0xc0,%%rcx				\n\t"\
		"subq		$0xc0,%%rdx				\n\t"\
		"vmovaps		(%%rax),%%ymm0				\n\t		vmovaps		0x80(%%rax),%%ymm8 			\n\t"\
		"vmovaps	0x20(%%rax),%%ymm1				\n\t		vmovaps		0xa0(%%rax),%%ymm9 			\n\t"\
		"vmovaps		(%%rax),%%ymm2				\n\t		vmovaps		0x80(%%rax),%%ymm10			\n\t"\
		"vmovaps	0x20(%%rax),%%ymm3				\n\t		vmovaps		0xa0(%%rax),%%ymm11			\n\t"\
	"vmovaps	(%%r8),%%ymm15	\n\t"/* 1.0 */\
		" vfmadd231pd		(%%rbx),%%ymm15,%%ymm0	\n\t		 vfmadd231pd	0x80(%%rbx),%%ymm15,%%ymm8 	\n\t"\
		" vfmadd231pd	0x20(%%rbx),%%ymm15,%%ymm1	\n\t		 vfmadd231pd	0xa0(%%rbx),%%ymm15,%%ymm9 	\n\t"\
		"vfnmadd231pd		(%%rbx),%%ymm15,%%ymm2	\n\t		vfnmadd231pd	0x80(%%rbx),%%ymm15,%%ymm10	\n\t"\
		"vfnmadd231pd	0x20(%%rbx),%%ymm15,%%ymm3	\n\t		vfnmadd231pd	0xa0(%%rbx),%%ymm15,%%ymm11	\n\t"\
		"vmovaps		(%%rcx),%%ymm4				\n\t		vmovaps		0x80(%%rcx),%%ymm12			\n\t"\
		"vmovaps	0x20(%%rcx),%%ymm5				\n\t		vmovaps		0xa0(%%rcx),%%ymm13			\n\t"\
		"vmovaps		(%%rcx),%%ymm6				\n\t		vmovaps		0x80(%%rcx),%%ymm14			\n\t"\
		"vmovaps	0x20(%%rcx),%%ymm7				\n\t		vmovaps		0xa0(%%rcx),%%ymm15			\n\t"\
	"vmovaps	%%ymm11,(%%rbx)	\n\t	vmovaps	(%%r8),%%ymm11	\n\t"/* spill to make room for 1.0 */\
		" vfmadd231pd		(%%rdx),%%ymm11,%%ymm4	\n\t		 vfmadd231pd	0x80(%%rdx),%%ymm11,%%ymm12	\n\t"\
		" vfmadd231pd	0x20(%%rdx),%%ymm11,%%ymm5	\n\t		 vfmadd231pd	0xa0(%%rdx),%%ymm11,%%ymm13	\n\t"\
		"vfnmadd231pd		(%%rdx),%%ymm11,%%ymm6	\n\t		vfnmadd231pd	0x80(%%rdx),%%ymm11,%%ymm14	\n\t"\
		"vfnmadd231pd	0x20(%%rdx),%%ymm11,%%ymm7	\n\t		vfnmadd231pd	0xa0(%%rdx),%%ymm11,%%ymm15	\n\t"\
		"vfmsub132pd		%%ymm11,%%ymm4,%%ymm0		\n\t		vfmsub132pd		%%ymm11,%%ymm12,%%ymm8 	\n\t"\
		"vfmsub132pd		%%ymm11,%%ymm5,%%ymm1		\n\t		vfmsub132pd		%%ymm11,%%ymm13,%%ymm9 	\n\t"\
		"vfmsub132pd		%%ymm11,%%ymm7,%%ymm2		\n\t		vfmsub132pd		%%ymm11,%%ymm15,%%ymm10	\n\t"\
		"vfmsub132pd		%%ymm11,%%ymm6,%%ymm3		\n\t		vfmsub132pd		(%%rbx),%%ymm14,%%ymm11	\n\t"\
		"vmovaps		%%ymm0 ,      (%%rbx)		\n\t		vmovaps			%%ymm8 , 0x080(%%rbx)	\n\t"\
		"vmovaps		%%ymm1 , 0x020(%%rbx)		\n\t		vmovaps			%%ymm9 , 0x0a0(%%rbx)	\n\t"\
		"vmovaps		%%ymm2 ,      (%%rdx)		\n\t		vmovaps			%%ymm10, 0x080(%%rdx)	\n\t"\
		"vmovaps		%%ymm3 , 0x020(%%rcx)		\n\t		vmovaps			%%ymm11, 0x0a0(%%rcx)	\n\t"\
	"vmovaps	%%ymm14,0xa0(%%rdx)	\n\t"/* spill ymm14 to make room for two */"	vmovaps	(%%rdi),%%ymm14	\n\t"/* two */\
	"vfmadd132pd	%%ymm14,%%ymm0,%%ymm4			\n\t	vfmadd132pd		%%ymm14,%%ymm8 ,%%ymm12		\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm1,%%ymm5			\n\t	vfmadd132pd		%%ymm14,%%ymm9 ,%%ymm13		\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm2,%%ymm7			\n\t	vfmadd132pd		%%ymm14,%%ymm10,%%ymm15		\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm3,%%ymm6			\n\t	vfmadd132pd	0xa0(%%rdx),%%ymm11,%%ymm14		\n\t"\
		"vmovaps	%%ymm4 ,    (%%rax)				\n\t		vmovaps		%%ymm12,0x80(%%rax)			\n\t"\
		"vmovaps	%%ymm5 ,0x20(%%rax)				\n\t		vmovaps		%%ymm13,0xa0(%%rax)			\n\t"\
		"vmovaps	%%ymm7 ,    (%%rcx)				\n\t		vmovaps		%%ymm15,0x80(%%rcx)			\n\t"\
		"vmovaps	%%ymm6 ,0x20(%%rdx)				\n\t		vmovaps		%%ymm14,0xa0(%%rdx)			\n\t"\
		"addq		$0x100,%%rax				\n\t"\
		"addq		$0x100,%%rbx				\n\t"\
		"addq		$0x100,%%rcx				\n\t"\
		"addq		$0x100,%%rdx				\n\t"\
		"vmovaps		(%%rax),%%ymm0				\n\t		vmovaps		0x80(%%rax),%%ymm8 			\n\t"\
		"vmovaps	0x20(%%rax),%%ymm1				\n\t		vmovaps		0xa0(%%rax),%%ymm9 			\n\t"\
		"vmovaps		(%%rax),%%ymm2				\n\t		vmovaps		0x80(%%rax),%%ymm10			\n\t"\
		"vmovaps	0x20(%%rax),%%ymm3				\n\t		vmovaps		0xa0(%%rax),%%ymm11			\n\t"\
	"vmovaps	(%%r8),%%ymm15	\n\t"/* 1.0 */\
		" vfmadd231pd		(%%rbx),%%ymm15,%%ymm0	\n\t		 vfmadd231pd	0x80(%%rbx),%%ymm15,%%ymm8 	\n\t"\
		" vfmadd231pd	0x20(%%rbx),%%ymm15,%%ymm1	\n\t		 vfmadd231pd	0xa0(%%rbx),%%ymm15,%%ymm9 	\n\t"\
		"vfnmadd231pd		(%%rbx),%%ymm15,%%ymm2	\n\t		vfnmadd231pd	0x80(%%rbx),%%ymm15,%%ymm10	\n\t"\
		"vfnmadd231pd	0x20(%%rbx),%%ymm15,%%ymm3	\n\t		vfnmadd231pd	0xa0(%%rbx),%%ymm15,%%ymm11	\n\t"\
		"vmovaps		(%%rcx),%%ymm4				\n\t		vmovaps		0x80(%%rcx),%%ymm12			\n\t"\
		"vmovaps	0x20(%%rcx),%%ymm5				\n\t		vmovaps		0xa0(%%rcx),%%ymm13			\n\t"\
		"vmovaps		(%%rcx),%%ymm6				\n\t		vmovaps		0x80(%%rcx),%%ymm14			\n\t"\
		"vmovaps	0x20(%%rcx),%%ymm7				\n\t		vmovaps		0xa0(%%rcx),%%ymm15			\n\t"\
	"vmovaps	%%ymm11,(%%rbx)	\n\t	vmovaps	(%%r8),%%ymm11	\n\t"/* spill to make room for 1.0 */\
		" vfmadd231pd		(%%rdx),%%ymm11,%%ymm4	\n\t		 vfmadd231pd	0x80(%%rdx),%%ymm11,%%ymm12	\n\t"\
		" vfmadd231pd	0x20(%%rdx),%%ymm11,%%ymm5	\n\t		 vfmadd231pd	0xa0(%%rdx),%%ymm11,%%ymm13	\n\t"\
		"vfnmadd231pd		(%%rdx),%%ymm11,%%ymm6	\n\t		vfnmadd231pd	0x80(%%rdx),%%ymm11,%%ymm14	\n\t"\
		"vfnmadd231pd	0x20(%%rdx),%%ymm11,%%ymm7	\n\t		vfnmadd231pd	0xa0(%%rdx),%%ymm11,%%ymm15	\n\t"\
		"vfmsub132pd		%%ymm11,%%ymm4,%%ymm0		\n\t		vfmsub132pd		%%ymm11,%%ymm12,%%ymm8 	\n\t"\
		"vfmsub132pd		%%ymm11,%%ymm5,%%ymm1		\n\t		vfmsub132pd		%%ymm11,%%ymm13,%%ymm9 	\n\t"\
		"vfmsub132pd		%%ymm11,%%ymm7,%%ymm2		\n\t		vfmsub132pd		%%ymm11,%%ymm15,%%ymm10	\n\t"\
		"vfmsub132pd		%%ymm11,%%ymm6,%%ymm3		\n\t		vfmsub132pd		(%%rbx),%%ymm14,%%ymm11	\n\t"\
	"vmovaps	%%ymm14,0xa0(%%rdx)	\n\t"/* spill ymm14 to make room for two */"	vmovaps	(%%rdi),%%ymm14	\n\t"/* two */\
	"vfmadd132pd	%%ymm14,%%ymm0,%%ymm4			\n\t	vfmadd132pd		%%ymm14,%%ymm8 ,%%ymm12		\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm1,%%ymm5			\n\t	vfmadd132pd		%%ymm14,%%ymm9 ,%%ymm13		\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm2,%%ymm7			\n\t	vfmadd132pd		%%ymm14,%%ymm10,%%ymm15		\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm3,%%ymm6			\n\t	vfmadd132pd	0xa0(%%rdx),%%ymm11,%%ymm14		\n\t"\
		"vmovaps	%%ymm0 ,      (%%rbx)			\n\t		vmovaps			%%ymm8 , 0x080(%%rbx)	\n\t"\
		"vmovaps	%%ymm1 , 0x020(%%rbx)			\n\t		vmovaps			%%ymm9 , 0x0a0(%%rbx)	\n\t"\
		"vmovaps	%%ymm4 ,      (%%rax)			\n\t		vmovaps			%%ymm12, 0x080(%%rax)	\n\t"\
		"vmovaps	%%ymm5 , 0x020(%%rax)			\n\t		vmovaps			%%ymm13, 0x0a0(%%rax)	\n\t"\
		"vmovaps	%%ymm3 ,%%ymm0					\n\t		vmovaps			%%ymm11,%%ymm8 			\n\t"\
		"vmovaps	%%ymm6 ,%%ymm1					\n\t"	/*	vmovaps			%%ymm14,%%ymm9 			\n\t"*/\
	"vmovaps	%%ymm14,(%%rdx)	\n\t	vmovaps	(%%r8),%%ymm9	\n\t"/* spill to make room for 1.0 */\
		"vfmsub132pd	%%ymm9 ,%%ymm7,%%ymm3			\n\t		vfmsub132pd		%%ymm9 ,%%ymm15,%%ymm11	\n\t"\
		"vfmsub132pd	%%ymm9 ,%%ymm2,%%ymm6			\n\t		vfmsub132pd		%%ymm9 ,%%ymm10,%%ymm14	\n\t"\
		"vfmadd132pd	%%ymm9 ,%%ymm7,%%ymm0			\n\t		vfmadd132pd		%%ymm9 ,%%ymm15,%%ymm8 	\n\t"\
		"vfmadd132pd	%%ymm9 ,%%ymm2,%%ymm1			\n\t		vfmadd132pd		(%%rdx),%%ymm10,%%ymm9 	\n\t"\
	"vmovaps	(%%rsi),%%ymm15	\n\t"/* isrt2 */\
		"vmulpd		%%ymm15,%%ymm3,%%ymm3			\n\t		vmulpd			%%ymm15,%%ymm11,%%ymm11	\n\t"\
		"vmulpd		%%ymm15,%%ymm6,%%ymm6			\n\t		vmulpd			%%ymm15,%%ymm14,%%ymm14	\n\t"\
		"vmulpd		%%ymm15,%%ymm0,%%ymm0			\n\t		vmulpd			%%ymm15,%%ymm8 ,%%ymm8 	\n\t"\
		"vmulpd		%%ymm15,%%ymm1,%%ymm1			\n\t		vmulpd			%%ymm15,%%ymm9 ,%%ymm9 	\n\t"\
		"vmovaps	%%ymm3 , 0x020(%%rcx)			\n\t		vmovaps			%%ymm11, 0x0a0(%%rcx)	\n\t"\
		"vmovaps	%%ymm6 , 0x020(%%rdx)			\n\t		vmovaps			%%ymm14, 0x0a0(%%rdx)	\n\t"\
		"vmovaps	%%ymm0 ,      (%%rcx)			\n\t		vmovaps			%%ymm8 , 0x080(%%rcx)	\n\t"\
		"vmovaps	%%ymm1 ,      (%%rdx)			\n\t		vmovaps			%%ymm9 , 0x080(%%rdx)	\n\t"\
		/*****	SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS *****/			/*****	SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS *****/\
		/*        (r02,r12,r22,r32,r0A,r1A,r2A,r3A):  */			/*        (r06,r16,r26,r36,r0E,r1E,r2E,r3E):  */\
		"vmovaps	-0x100(%%rax),%%ymm0			\n\t		vmovaps		-0x080(%%rax),%%ymm8 		\n\t"\
		"vmovaps	-0x100(%%rbx),%%ymm4			\n\t		vmovaps		-0x080(%%rbx),%%ymm12		\n\t"\
		"vmovaps	-0x0e0(%%rax),%%ymm1			\n\t		vmovaps		-0x060(%%rax),%%ymm9 		\n\t"\
		"vmovaps	-0x0e0(%%rbx),%%ymm5			\n\t		vmovaps		-0x060(%%rbx),%%ymm13		\n\t"\
		"vmovaps		  (%%rax),%%ymm2			\n\t		vmovaps		 0x080(%%rax),%%ymm10		\n\t"\
		"vmovaps	 0x020(%%rbx),%%ymm7			\n\t		vmovaps		 0x0a0(%%rbx),%%ymm15		\n\t"\
		"vmovaps	 0x020(%%rax),%%ymm3			\n\t		vmovaps		 0x0a0(%%rax),%%ymm11		\n\t"\
		"vmovaps		  (%%rbx),%%ymm6			\n\t		vmovaps		 0x080(%%rbx),%%ymm14		\n\t"\
	"vmovaps	%%ymm13,(%%rbx)	\n\t	vmovaps	(%%r8),%%ymm13	\n\t"/* spill to make room for 1.0 */\
		"vfmsub132pd	%%ymm13,%%ymm2,%%ymm0			\n\t		vfmsub132pd	%%ymm13,%%ymm10,%%ymm8 		\n\t"\
		"vfmsub132pd	%%ymm13,%%ymm7,%%ymm4			\n\t		vfmsub132pd	%%ymm13,%%ymm15,%%ymm12		\n\t"\
		"vfmsub132pd	%%ymm13,%%ymm3,%%ymm1			\n\t		vfmsub132pd	%%ymm13,%%ymm11,%%ymm9 		\n\t"\
		"vfmsub132pd	%%ymm13,%%ymm6,%%ymm5			\n\t		vfmsub132pd	(%%rbx),%%ymm14,%%ymm13		\n\t"\
	"vmovaps	%%ymm14,0xa0(%%rbx)	\n\t"/* spill ymm14 to make room for two */"	vmovaps	(%%rdi),%%ymm14	\n\t"/* two */\
	"vfmadd132pd	%%ymm14,%%ymm0,%%ymm2			\n\t	vfmadd132pd		%%ymm14,%%ymm8 ,%%ymm10		\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm4,%%ymm7			\n\t	vfmadd132pd		%%ymm14,%%ymm12,%%ymm15		\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm1,%%ymm3			\n\t	vfmadd132pd		%%ymm14,%%ymm9 ,%%ymm11		\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm5,%%ymm6			\n\t	vfmadd132pd	0xa0(%%rbx),%%ymm13,%%ymm14		\n\t"\
		"vmovaps	%%ymm0,      (%%rax)			\n\t		vmovaps		%%ymm8 , 0x080(%%rax)		\n\t"\
		"vmovaps	%%ymm4,      (%%rbx)			\n\t		vmovaps		%%ymm12, 0x080(%%rbx)		\n\t"\
		"vmovaps	%%ymm1, 0x020(%%rax)			\n\t		vmovaps		%%ymm9 , 0x0a0(%%rax)		\n\t"\
		"vmovaps	%%ymm5,-0x0e0(%%rbx)			\n\t		vmovaps		%%ymm13,-0x060(%%rbx)		\n\t"\
		"vmovaps	%%ymm2,-0x100(%%rax)			\n\t		vmovaps		%%ymm10,-0x080(%%rax)		\n\t"\
		"vmovaps	%%ymm7,-0x100(%%rbx)			\n\t		vmovaps		%%ymm15,-0x080(%%rbx)		\n\t"\
		"vmovaps	%%ymm3,-0x0e0(%%rax)			\n\t		vmovaps		%%ymm11,-0x060(%%rax)		\n\t"\
		"vmovaps	%%ymm6, 0x020(%%rbx)			\n\t		vmovaps		%%ymm14, 0x0a0(%%rbx)		\n\t"\
		"vmovaps	-0x100(%%rcx),%%ymm0			\n\t		vmovaps		-0x080(%%rcx),%%ymm8 		\n\t"\
		"vmovaps	-0x100(%%rdx),%%ymm4			\n\t		vmovaps		-0x080(%%rdx),%%ymm12		\n\t"\
		"vmovaps	-0x0e0(%%rcx),%%ymm1			\n\t		vmovaps		-0x060(%%rcx),%%ymm9 		\n\t"\
		"vmovaps	-0x0e0(%%rdx),%%ymm5			\n\t		vmovaps		-0x060(%%rdx),%%ymm13		\n\t"\
		"vmovaps		  (%%rcx),%%ymm2			\n\t		vmovaps		 0x080(%%rcx),%%ymm10		\n\t"\
		"vmovaps	 0x020(%%rdx),%%ymm7			\n\t		vmovaps		 0x0a0(%%rdx),%%ymm15		\n\t"\
		"vmovaps	 0x020(%%rcx),%%ymm3			\n\t		vmovaps		 0x0a0(%%rcx),%%ymm11		\n\t"\
		"vmovaps		  (%%rdx),%%ymm6			\n\t		vmovaps		 0x080(%%rdx),%%ymm14		\n\t"\
	"vmovaps	%%ymm13,(%%rdx)	\n\t	vmovaps	(%%r8),%%ymm13	\n\t"/* spill to make room for 1.0 */\
		"vfmsub132pd	%%ymm13,%%ymm2,%%ymm0			\n\t		vfmsub132pd	%%ymm13,%%ymm10,%%ymm8 		\n\t"\
		"vfmsub132pd	%%ymm13,%%ymm7,%%ymm4			\n\t		vfmsub132pd	%%ymm13,%%ymm15,%%ymm12		\n\t"\
		"vfmsub132pd	%%ymm13,%%ymm3,%%ymm1			\n\t		vfmsub132pd	%%ymm13,%%ymm11,%%ymm9 		\n\t"\
		"vfmsub132pd	%%ymm13,%%ymm6,%%ymm5			\n\t		vfmsub132pd	(%%rdx),%%ymm14,%%ymm13		\n\t"\
	"vmovaps	%%ymm14,0xa0(%%rdx)	\n\t"/* spill ymm14 to make room for two */"	vmovaps	(%%rdi),%%ymm14	\n\t"/* two */\
	"vfmadd132pd	%%ymm14,%%ymm0,%%ymm2			\n\t	vfmadd132pd		%%ymm14,%%ymm8 ,%%ymm10		\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm4,%%ymm7			\n\t	vfmadd132pd		%%ymm14,%%ymm12,%%ymm15		\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm1,%%ymm3			\n\t	vfmadd132pd		%%ymm14,%%ymm9 ,%%ymm11		\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm5,%%ymm6			\n\t	vfmadd132pd	0xa0(%%rdx),%%ymm13,%%ymm14		\n\t"\
		"vmovaps	%%ymm0,      (%%rcx)			\n\t		vmovaps		%%ymm8 , 0x080(%%rcx)		\n\t"\
		"vmovaps	%%ymm4,      (%%rdx)			\n\t		vmovaps		%%ymm12, 0x080(%%rdx)		\n\t"\
		"vmovaps	%%ymm1, 0x020(%%rcx)			\n\t		vmovaps		%%ymm9 , 0x0a0(%%rcx)		\n\t"\
		"vmovaps	%%ymm5,-0x0e0(%%rdx)			\n\t		vmovaps		%%ymm13,-0x060(%%rdx)		\n\t"\
		"vmovaps	%%ymm2,-0x100(%%rcx)			\n\t		vmovaps		%%ymm10,-0x080(%%rcx)		\n\t"\
		"vmovaps	%%ymm7,-0x100(%%rdx)			\n\t		vmovaps		%%ymm15,-0x080(%%rdx)		\n\t"\
		"vmovaps	%%ymm3,-0x0e0(%%rcx)			\n\t		vmovaps		%%ymm11,-0x060(%%rcx)		\n\t"\
		"vmovaps	%%ymm6, 0x020(%%rdx)			\n\t		vmovaps		%%ymm14, 0x0a0(%%rdx)		\n\t"\
		"\n\t"\
	/***************************************************************************************/\
	/* Now do eight more radix-4 transforms, including the internal and external twiddles: */\
	/***************************************************************************************/\
		"\n\t"\
		/* Using the upper block(s) of the main array for temp-storage in the section below led to a nasty AVX bug to track down: */\
		/* In fermat-mod mode the 4 block addresses in ascending order are add0,1,2,3 with no 'gaps' between blocks, whereas for */\
		/* mersenne-mod the addresses in asc. order are add0,2,3,1 with a gap between contiguous-data-block pairs 0,2 and 3,1. Thus */\
		/* for fermat-mod we need [add2] as the base-address of the 'high-half' block for temp-storage; for mersenne-mod we need [add3]. */\
		/* In both cases we have that add2 < add3 so instead use (add2 - add1): > 0 for fermat-mod, < 0 for mersenne - to differentiate: */\
		"movq	%[__add2],%%rsi		\n\t"/* destroyable copy of add2 */\
		"movq	%[__add2],%%rbx		\n\t"\
		"subq	%[__add1],%%rsi		\n\t"/* rsi = (add2 - add1); if this yields a borrow (i.e. sets CF) it's mersenne, else fermat. */\
		"cmovcq %[__add3],%%rbx	\n\t" /* if CF set (i.e. h > l), copy source [add3] into dest (rbx), else leave dest = [add2]. */\
		"\n\t"\
	/************************************************/				/************************************************/\
	/*...Block 2: t02,t12,t22,t32 -> r10,14,12,16:  */				/*...Block 6: t0A,t1A,t2A,t3A -> r18,1C,1A,1E:  */\
	/************************************************/				/************************************************/\
		"movq		%[__isrt2],%%rsi					\n\t"\
		"movq		%[__r10],%%rax	/* base-addr in rcol = c05/r18, so rax/r10 offset +0x100 vs lcol */\n\t"\
		"movq		%%rsi,%%rcx					\n\t"\
		"movq		%%rsi,%%rdx					\n\t"\
		"movq		%[__c01],%%r10					\n\t"\
		"addq		$0x020,%%rcx	/* cc0 */		\n\t"\
		"addq		$0x060,%%rdx	/* cc1 */		\n\t"\
		"vmovaps		0x040(%%rax),%%ymm4					\n\t		vmovaps		 0x140(%%rax),%%ymm12			\n\t"\
		"vmovaps		0x0c0(%%rax),%%ymm0					\n\t		vmovaps		 0x1c0(%%rax),%%ymm8 			\n\t"\
		"vmovaps		0x060(%%rax),%%ymm5					\n\t		vmovaps		 0x160(%%rax),%%ymm15			\n\t"\
		"vmovaps		0x0e0(%%rax),%%ymm3					\n\t"	/*	vmovaps		 0x1e0(%%rax),%%ymm9 			\n\t"*/\
		"vmovaps			%%ymm4 ,%%ymm6					\n\t		vmovaps		 	%%ymm12,%%ymm14				\n\t"\
		"vmovaps			%%ymm0 ,%%ymm2					\n\t"	/*	vmovaps		 	%%ymm8 ,%%ymm10				\n\t"*/\
		"vmovaps			%%ymm5 ,%%ymm7					\n\t"	/*	vmovaps		 	%%ymm13,%%ymm15				\n\t"*/\
	/*	"vmovaps			%%ymm1 ,%%ymm3			*/			"		vmovaps		 	%%ymm9 ,%%ymm11				\n\t"\
	/* Strategic register-usage for sincos consts here cuts loads from 16 to 4: */\
		"vmovaps			 (%%rdx),%%ymm9 			\n\t		vmovaps		0x020(%%rdx),%%ymm10			\n\t"\
		"vmovaps		0x040(%%rdx),%%ymm1 			\n\t		vmovaps		0x060(%%rdx),%%ymm13			\n\t"\
		"vmulpd			%%ymm10,%%ymm7,%%ymm7				\n\t		vmulpd			 %%ymm1 ,%%ymm15,%%ymm15	\n\t"\
		"vmulpd			%%ymm13,%%ymm3,%%ymm3				\n\t		vmulpd			 %%ymm10,%%ymm11,%%ymm11	\n\t"\
		"vmulpd			%%ymm10,%%ymm6,%%ymm6				\n\t		vmulpd			 %%ymm1 ,%%ymm14,%%ymm14	\n\t"\
		"vmulpd			%%ymm13,%%ymm2,%%ymm2				\n\t		vmulpd		0x1c0(%%rax),%%ymm10,%%ymm10	\n\t"\
	"vfmadd132pd		%%ymm9 ,%%ymm7,%%ymm4				\n\t	vfmadd132pd	 		 %%ymm13,%%ymm15,%%ymm12	\n\t"\
	"vfmadd132pd		%%ymm1 ,%%ymm3,%%ymm0				\n\t	vfmsub132pd	 		 %%ymm9 ,%%ymm11,%%ymm8 	\n\t"\
	"vfmsub132pd		%%ymm9 ,%%ymm6,%%ymm5				\n\t	vfmsub132pd	 	0x160(%%rax),%%ymm14,%%ymm13	\n\t"\
	"vfmsub132pd	0xe0(%%rax),%%ymm2,%%ymm1				\n\t	vfmadd132pd	 	0x1e0(%%rax),%%ymm10,%%ymm9 	\n\t"\
		"vmovaps		%%ymm5,%%ymm7					\n\t		vmovaps		%%ymm13,%%ymm15					\n\t"\
		"vmovaps		%%ymm4,%%ymm6					\n\t		vmovaps		%%ymm12,%%ymm14					\n\t"\
	"vmovaps	%%ymm15,(%%rbx)	\n\t	vmovaps	(%%r8),%%ymm15	\n\t"/* spill to make room for 1.0 */\
		"vfmsub132pd	%%ymm15,%%ymm0,%%ymm4			\n\t		vfmadd132pd	%%ymm15,%%ymm8,%%ymm12			\n\t"\
		"vfmsub132pd	%%ymm15,%%ymm1,%%ymm5			\n\t		vfmadd132pd	%%ymm15,%%ymm9,%%ymm13			\n\t"\
		"vfmadd132pd	%%ymm15,%%ymm0,%%ymm6			\n\t		vfmsub132pd	%%ymm15,%%ymm8,%%ymm14			\n\t"\
		"vfmadd132pd	%%ymm15,%%ymm1,%%ymm7			\n\t		vfmsub132pd	(%%rbx),%%ymm9,%%ymm15			\n\t"\
		"vmovaps		 0x080(%%rax),%%ymm0				\n\t		vmovaps		 0x180(%%rax),%%ymm8 				\n\t"\
		"vmovaps		 0x0a0(%%rax),%%ymm1				\n\t		vmovaps		 0x1a0(%%rax),%%ymm9 				\n\t"\
		"vmovaps		 0x080(%%rax),%%ymm2				\n\t		vmovaps		 0x180(%%rax),%%ymm10				\n\t"\
		"vmovaps		 0x0a0(%%rax),%%ymm3				\n\t		vmovaps		 0x1a0(%%rax),%%ymm11				\n\t"\
		"vmulpd		 0x020(%%rcx),%%ymm0,%%ymm0			\n\t		vmulpd			  (%%rcx),%%ymm8 ,%%ymm8 			\n\t"\
		"vmulpd		 0x020(%%rcx),%%ymm1,%%ymm1			\n\t		vmulpd			  (%%rcx),%%ymm9 ,%%ymm9 			\n\t"\
	"vfmsub132pd		  (%%rcx),%%ymm0,%%ymm3			\n\t	vfmadd132pd		 0x020(%%rcx),%%ymm8 ,%%ymm11			\n\t"\
	"vfmadd132pd		  (%%rcx),%%ymm1,%%ymm2			\n\t	vfmsub132pd		 0x020(%%rcx),%%ymm9 ,%%ymm10			\n\t"\
		"vmovaps			  (%%rax),%%ymm0				\n\t		vmovaps		 0x100(%%rax),%%ymm8 				\n\t"\
		"vmovaps		 0x020(%%rax),%%ymm1				\n\t		vmovaps		 0x120(%%rax),%%ymm9 				\n\t"\
	"vmovaps	%%ymm12,    (%%rbx)	\n\t	vmovaps	(%%r8) ,%%ymm12	\n\t"/* spill to make room for 1.0 */\
	"vmovaps	%%ymm13,0x20(%%rbx)	\n\t	vmovaps	(%%rdi),%%ymm13	\n\t"/* spill to make room for 2.0 */\
		"vfmsub132pd	%%ymm12,%%ymm2,%%ymm0				\n\t		vfmsub132pd	%%ymm12,%%ymm10,%%ymm8 			\n\t"\
		"vfmsub132pd	%%ymm12,%%ymm3,%%ymm1				\n\t		vfmsub132pd	%%ymm12,%%ymm11,%%ymm9 			\n\t"\
	"vfmadd132pd		%%ymm13,%%ymm0,%%ymm2			\n\t	vfmadd132pd		%%ymm13,%%ymm8 ,%%ymm10		\n\t"\
	"vfmadd132pd		%%ymm13,%%ymm1,%%ymm3			\n\t	vfmadd132pd		%%ymm13,%%ymm9 ,%%ymm11		\n\t"\
		"vfmsub132pd	%%ymm12,%%ymm6,%%ymm2				\n\t		vfmsub132pd	%%ymm12,%%ymm14,%%ymm8 			\n\t"\
		"vfmsub132pd	%%ymm12,%%ymm7,%%ymm3				\n\t		vfmsub132pd	%%ymm12,%%ymm15,%%ymm9 			\n\t"\
	"vfmadd132pd		%%ymm13,%%ymm2,%%ymm6			\n\t	vfmadd132pd		%%ymm13,%%ymm8 ,%%ymm14		\n\t"\
	"vfmadd132pd		%%ymm13,%%ymm3,%%ymm7			\n\t	vfmadd132pd		%%ymm13,%%ymm9 ,%%ymm15		\n\t"\
	"vmovaps		(%%rbx),%%ymm12	\n\t"/* restore spill */\
	"vmovaps	0x20(%%rbx),%%ymm13	\n\t"/* restore spill */\
		"vmovaps		%%ymm2, 0x040(%%rax)			\n\t		vmovaps		%%ymm8 , 0x140(%%rax)			\n\t"\
		"vmovaps		%%ymm3, 0x060(%%rax)			\n\t		vmovaps		%%ymm9 , 0x160(%%rax)			\n\t"\
		"vmovaps		%%ymm6,%%ymm2					\n\t		vmovaps		%%ymm14,%%ymm8 					\n\t"\
		"vmovaps		%%ymm7,%%ymm3					\n\t		vmovaps		%%ymm15,%%ymm9 					\n\t"\
		"vmulpd			  (%%r10),%%ymm6,%%ymm6			\n\t		vmulpd		 0x100(%%r10),%%ymm14,%%ymm14			\n\t"\
		"vmulpd			  (%%r10),%%ymm7,%%ymm7			\n\t		vmulpd		 0x100(%%r10),%%ymm15,%%ymm15			\n\t"\
	" vfmadd231pd	 0x020(%%r10),%%ymm3,%%ymm6			\n\t	 vfmadd231pd	 0x120(%%r10),%%ymm9 ,%%ymm14			\n\t"\
	"vfnmadd231pd	 0x020(%%r10),%%ymm2,%%ymm7			\n\t	vfnmadd231pd	 0x120(%%r10),%%ymm8 ,%%ymm15			\n\t"\
		"vmovaps		%%ymm7, 0x020(%%rbx)			\n\t		vmovaps		%%ymm15, 0x0a0(%%rbx)			\n\t"\
		"vmovaps		%%ymm6,      (%%rbx)			\n\t		vmovaps		%%ymm14, 0x080(%%rbx)			\n\t"\
		"vmovaps		 0x040(%%rax),%%ymm6				\n\t		vmovaps		 0x140(%%rax),%%ymm14				\n\t"\
		"vmovaps		 0x060(%%rax),%%ymm7				\n\t		vmovaps		 0x160(%%rax),%%ymm15				\n\t"\
		"vmovaps		%%ymm6,%%ymm2					\n\t		vmovaps		%%ymm14,%%ymm8 					\n\t"\
		"vmovaps		%%ymm7,%%ymm3					\n\t		vmovaps		%%ymm15,%%ymm9 					\n\t"\
		"vmulpd		 0x040(%%r10),%%ymm6,%%ymm6			\n\t		vmulpd		 0x140(%%r10),%%ymm14,%%ymm14			\n\t"\
		"vmulpd		 0x040(%%r10),%%ymm7,%%ymm7			\n\t		vmulpd		 0x140(%%r10),%%ymm15,%%ymm15			\n\t"\
		"vfmsub132pd	(%%r8),%%ymm5,%%ymm0			\n\t		vfmsub132pd	(%%r8),%%ymm13,%%ymm10			\n\t"\
		"vfmsub132pd	(%%r8),%%ymm4,%%ymm1			\n\t		vfmsub132pd	(%%r8),%%ymm12,%%ymm11			\n\t"\
	" vfmadd231pd	 0x060(%%r10),%%ymm3,%%ymm6			\n\t	 vfmadd231pd	 0x160(%%r10),%%ymm9 ,%%ymm14			\n\t"\
	"vfnmadd231pd	 0x060(%%r10),%%ymm2,%%ymm7			\n\t	vfnmadd231pd	 0x160(%%r10),%%ymm8 ,%%ymm15			\n\t"\
	"vfmadd132pd		  (%%rdi),%%ymm0,%%ymm5			\n\t	vfmadd132pd		(%%rdi),%%ymm10,%%ymm13			\n\t"\
	"vfmadd132pd		  (%%rdi),%%ymm1,%%ymm4			\n\t	vfmadd132pd		(%%rdi),%%ymm11,%%ymm12			\n\t"\
		"vmovaps		%%ymm7, 0x220(%%rbx)			\n\t		vmovaps		%%ymm15, 0x2a0(%%rbx)			\n\t"\
		"vmovaps		%%ymm6, 0x200(%%rbx)			\n\t		vmovaps		%%ymm14, 0x280(%%rbx)			\n\t"\
		"addq		$0x080,%%r10					\n\t"\
		"vmovaps		%%ymm5,%%ymm2					\n\t		vmovaps		%%ymm13,%%ymm8 					\n\t"\
		"vmovaps		%%ymm1,%%ymm3					\n\t		vmovaps		%%ymm11,%%ymm9 					\n\t"\
		"vmovaps		%%ymm0,%%ymm6					\n\t		vmovaps		%%ymm10,%%ymm14					\n\t"\
		"vmovaps		%%ymm4,%%ymm7					\n\t		vmovaps		%%ymm12,%%ymm15					\n\t"\
		"vmulpd			  (%%r10),%%ymm5,%%ymm5			\n\t		vmulpd		 0x100(%%r10),%%ymm13,%%ymm13			\n\t"\
		"vmulpd			  (%%r10),%%ymm1,%%ymm1			\n\t		vmulpd		 0x100(%%r10),%%ymm11,%%ymm11			\n\t"\
		"vmulpd		 0x040(%%r10),%%ymm0,%%ymm0			\n\t		vmulpd		 0x140(%%r10),%%ymm10,%%ymm10			\n\t"\
		"vmulpd		 0x040(%%r10),%%ymm4,%%ymm4			\n\t		vmulpd		 0x140(%%r10),%%ymm12,%%ymm12			\n\t"\
	" vfmadd231pd	 0x020(%%r10),%%ymm3,%%ymm5			\n\t	 vfmadd231pd	 0x120(%%r10),%%ymm9 ,%%ymm13			\n\t"\
	"vfnmadd231pd	 0x020(%%r10),%%ymm2,%%ymm1			\n\t	vfnmadd231pd	 0x120(%%r10),%%ymm8 ,%%ymm11			\n\t"\
	" vfmadd231pd	 0x060(%%r10),%%ymm7,%%ymm0			\n\t	 vfmadd231pd	 0x160(%%r10),%%ymm15,%%ymm10			\n\t"\
	"vfnmadd231pd	 0x060(%%r10),%%ymm6,%%ymm4			\n\t	vfnmadd231pd	 0x160(%%r10),%%ymm14,%%ymm12			\n\t"\
		"vmovaps		%%ymm1, 0x120(%%rbx)			\n\t		vmovaps		%%ymm11, 0x1a0(%%rbx)			\n\t"\
		"vmovaps		%%ymm5, 0x100(%%rbx)			\n\t		vmovaps		%%ymm13, 0x180(%%rbx)			\n\t"\
		"vmovaps		%%ymm4, 0x320(%%rbx)			\n\t		vmovaps		%%ymm12, 0x3a0(%%rbx)			\n\t"\
		"vmovaps		%%ymm0, 0x300(%%rbx)			\n\t		vmovaps		%%ymm10, 0x380(%%rbx)			\n\t"\
		"\n\t"\
	/************************************************/				/************************************************/\
	/*...Block 4: t06,t16,t26,t36 -> r30,34,32,36:  */				/*...Block 8: t0E,t1E,t2E,t3E -> r38,3C,3A,3E:  */\
	/************************************************/				/************************************************/\
		"addq		$0x400,%%rax						\n\t		addq		$0x100,%%r10						\n\t"\
		"vmovaps		0x040(%%rax),%%ymm4					\n\t		vmovaps		 0x140(%%rax),%%ymm12			\n\t"\
		"vmovaps		0x0c0(%%rax),%%ymm0					\n\t		vmovaps		 0x1c0(%%rax),%%ymm10			\n\t"\
		"vmovaps		0x060(%%rax),%%ymm7					\n\t		vmovaps		 0x160(%%rax),%%ymm15			\n\t"\
		"vmovaps		0x0e0(%%rax),%%ymm1					\n\t		vmovaps		 0x1e0(%%rax),%%ymm9 			\n\t"\
		"vmovaps			%%ymm4 ,%%ymm6					\n\t"	/*	vmovaps		 	%%ymm12,%%ymm14				\n\t"*/\
		"vmovaps			%%ymm0 ,%%ymm2					\n\t"	/*	vmovaps		 	%%ymm8 ,%%ymm10				\n\t"*/\
	/*	"vmovaps			%%ymm5 ,%%ymm7	*/						/*	vmovaps		 	%%ymm13,%%ymm15				\n\t"*/\
		"vmovaps			%%ymm1 ,%%ymm3					\n\t		vmovaps		 	%%ymm9 ,%%ymm11				\n\t"\
	/* Strategic register-usage for sincos consts here cuts loads from 16 to 4: */\
		"vmovaps			 (%%rdx),%%ymm14			\n\t		vmovaps		0x020(%%rdx),%%ymm13			\n\t"\
		"vmovaps		0x040(%%rdx),%%ymm5 			\n\t		vmovaps		0x060(%%rdx),%%ymm8 			\n\t"\
	/* We cyc-permuted the 4 paired-MUL lines so as to put a '     (%%rdx)' in the SRC3 slot of the last line's rcol: */\
		"vmulpd		 	 %%ymm14,%%ymm2,%%ymm2				\n\t		vmulpd			 %%ymm5 ,%%ymm10,%%ymm10	\n\t"\
		"vmulpd		 	 %%ymm8 ,%%ymm7,%%ymm7				\n\t		vmulpd			 %%ymm14,%%ymm15,%%ymm15	\n\t"\
		"vmulpd		 	 %%ymm14,%%ymm3,%%ymm3				\n\t		vmulpd			 %%ymm5 ,%%ymm11,%%ymm11	\n\t"\
		"vmulpd		 	 %%ymm8 ,%%ymm6,%%ymm6				\n\t		vmulpd		0x140(%%rax),%%ymm14,%%ymm14	\n\t"\
	/* ...and similarly cyc-permute these 4 lines so middle (SRC2) op's reg-index order matches that of the above 4: */\
	"vfmadd132pd	 	 %%ymm13,%%ymm2,%%ymm1				\n\t	vfmsub132pd	 		 %%ymm8 ,%%ymm10,%%ymm9 	\n\t"\
	"vfmadd132pd	 	 %%ymm5 ,%%ymm7,%%ymm4				\n\t	vfmadd132pd	 		 %%ymm13,%%ymm15,%%ymm12	\n\t"\
	"vfmsub132pd	 	 %%ymm13,%%ymm3,%%ymm0				\n\t	vfmadd132pd	 	0x1c0(%%rax),%%ymm11,%%ymm8 	\n\t"\
	"vfmsub132pd	 0x60(%%rax),%%ymm6,%%ymm5				\n\t	vfmsub132pd	 	0x160(%%rax),%%ymm14,%%ymm13	\n\t"\
		"vmovaps		%%ymm5,%%ymm7					\n\t		vmovaps		%%ymm13,%%ymm15					\n\t"\
		"vmovaps		%%ymm4,%%ymm6					\n\t		vmovaps		%%ymm12,%%ymm14					\n\t"\
	"vmovaps	%%ymm15,0x40(%%rax)	\n\t	vmovaps	(%%r8),%%ymm15	\n\t"/* spill to make room for 1.0 */\
		"vfmadd132pd	%%ymm15,%%ymm0,%%ymm4			\n\t		vfmadd132pd		%%ymm15,%%ymm8,%%ymm12			\n\t"\
		"vfmadd132pd	%%ymm15,%%ymm1,%%ymm5			\n\t		vfmadd132pd		%%ymm15,%%ymm9,%%ymm13			\n\t"\
		"vfmsub132pd	%%ymm15,%%ymm0,%%ymm6			\n\t		vfmsub132pd		%%ymm15,%%ymm8,%%ymm14			\n\t"\
		"vfmsub132pd	%%ymm15,%%ymm1,%%ymm7			\n\t		vfmsub132pd	0x40(%%rax),%%ymm9,%%ymm15			\n\t"\
		"vmovaps		 0x080(%%rax),%%ymm0				\n\t		vmovaps		 0x180(%%rax),%%ymm8 				\n\t"\
		"vmovaps		 0x0a0(%%rax),%%ymm1				\n\t		vmovaps		 0x1a0(%%rax),%%ymm9 				\n\t"\
		"vmovaps		 0x080(%%rax),%%ymm2				\n\t		vmovaps		 0x180(%%rax),%%ymm10				\n\t"\
		"vmovaps		 0x0a0(%%rax),%%ymm3				\n\t		vmovaps		 0x1a0(%%rax),%%ymm11				\n\t"\
		"vmulpd			  (%%rcx),%%ymm0,%%ymm0			\n\t		vmulpd		 0x020(%%rcx),%%ymm8 ,%%ymm8 			\n\t"\
		"vmulpd			  (%%rcx),%%ymm1,%%ymm1			\n\t		vmulpd		 0x020(%%rcx),%%ymm9 ,%%ymm9 			\n\t"\
	"vfmsub132pd	 0x020(%%rcx),%%ymm0,%%ymm3			\n\t	vfmadd132pd			  (%%rcx),%%ymm8 ,%%ymm11			\n\t"\
	"vfmadd132pd	 0x020(%%rcx),%%ymm1,%%ymm2			\n\t	vfmsub132pd			  (%%rcx),%%ymm9 ,%%ymm10			\n\t"\
		"vmovaps			  (%%rax),%%ymm0				\n\t		vmovaps		 0x100(%%rax),%%ymm8 				\n\t"\
		"vmovaps		 0x020(%%rax),%%ymm1				\n\t		vmovaps		 0x120(%%rax),%%ymm9 				\n\t"\
	"vmovaps	%%ymm12,0x40(%%rbx)	\n\t	vmovaps	(%%r8) ,%%ymm12	\n\t"/* spill to make room for 1.0 */\
	"vmovaps	%%ymm13,0x60(%%rbx)	\n\t	vmovaps	(%%rdi),%%ymm13	\n\t"/* spill to make room for 2.0 */\
		"vfmsub132pd	%%ymm12,%%ymm2,%%ymm0				\n\t		vfmsub132pd	%%ymm12,%%ymm10,%%ymm8 			\n\t"\
		"vfmsub132pd	%%ymm12,%%ymm3,%%ymm1				\n\t		vfmsub132pd	%%ymm12,%%ymm11,%%ymm9 			\n\t"\
	"vfmadd132pd		%%ymm13,%%ymm0,%%ymm2			\n\t	vfmadd132pd		%%ymm13,%%ymm8 ,%%ymm10		\n\t"\
	"vfmadd132pd		%%ymm13,%%ymm1,%%ymm3			\n\t	vfmadd132pd		%%ymm13,%%ymm9 ,%%ymm11		\n\t"\
		"vfmsub132pd	%%ymm12,%%ymm6,%%ymm2				\n\t		vfmsub132pd	%%ymm12,%%ymm14,%%ymm8 			\n\t"\
		"vfmsub132pd	%%ymm12,%%ymm7,%%ymm3				\n\t		vfmsub132pd	%%ymm12,%%ymm15,%%ymm9 			\n\t"\
	"vfmadd132pd		%%ymm13,%%ymm2,%%ymm6			\n\t	vfmadd132pd		%%ymm13,%%ymm8 ,%%ymm14		\n\t"\
	"vfmadd132pd		%%ymm13,%%ymm3,%%ymm7			\n\t	vfmadd132pd		%%ymm13,%%ymm9 ,%%ymm15		\n\t"\
	"vmovaps	0x40(%%rbx),%%ymm12	\n\t"/* restore spill */\
	"vmovaps	0x60(%%rbx),%%ymm13	\n\t"/* restore spill */\
		"addq		$0x080,%%r10						\n\t"\
		"vmovaps		%%ymm2, 0x040(%%rax)			\n\t		vmovaps		%%ymm8 , 0x140(%%rax)			\n\t"\
		"vmovaps		%%ymm3, 0x060(%%rax)			\n\t		vmovaps		%%ymm9 , 0x160(%%rax)			\n\t"\
		"vmovaps		%%ymm6,%%ymm2					\n\t		vmovaps		%%ymm14,%%ymm8 					\n\t"\
		"vmovaps		%%ymm7,%%ymm3					\n\t		vmovaps		%%ymm15,%%ymm9 					\n\t"\
		"vmulpd			  (%%r10),%%ymm6,%%ymm6			\n\t		vmulpd		 0x100(%%r10),%%ymm14,%%ymm14			\n\t"\
		"vmulpd			  (%%r10),%%ymm7,%%ymm7			\n\t		vmulpd		 0x100(%%r10),%%ymm15,%%ymm15			\n\t"\
	" vfmadd231pd	 0x020(%%r10),%%ymm3,%%ymm6			\n\t	 vfmadd231pd	 0x120(%%r10),%%ymm9 ,%%ymm14			\n\t"\
	"vfnmadd231pd	 0x020(%%r10),%%ymm2,%%ymm7			\n\t	vfnmadd231pd	 0x120(%%r10),%%ymm8 ,%%ymm15			\n\t"\
		"vmovaps		%%ymm7, 0x060(%%rbx)			\n\t		vmovaps		%%ymm15, 0x0e0(%%rbx)			\n\t"\
		"vmovaps		%%ymm6, 0x040(%%rbx)			\n\t		vmovaps		%%ymm14, 0x0c0(%%rbx)			\n\t"\
		"vmovaps		 0x040(%%rax),%%ymm6				\n\t		vmovaps		 0x140(%%rax),%%ymm14				\n\t"\
		"vmovaps		 0x060(%%rax),%%ymm7				\n\t		vmovaps		 0x160(%%rax),%%ymm15				\n\t"\
		"vmovaps		%%ymm6,%%ymm2					\n\t		vmovaps		%%ymm14,%%ymm8 					\n\t"\
		"vmovaps		%%ymm7,%%ymm3					\n\t		vmovaps		%%ymm15,%%ymm9 					\n\t"\
		"vmulpd		 0x040(%%r10),%%ymm6,%%ymm6			\n\t		vmulpd		 0x140(%%r10),%%ymm14,%%ymm14			\n\t"\
		"vmulpd		 0x040(%%r10),%%ymm7,%%ymm7			\n\t		vmulpd		 0x140(%%r10),%%ymm15,%%ymm15			\n\t"\
		"vfmsub132pd	(%%r8),%%ymm5,%%ymm0			\n\t		vfmsub132pd	(%%r8),%%ymm13,%%ymm10			\n\t"\
		"vfmsub132pd	(%%r8),%%ymm4,%%ymm1			\n\t		vfmsub132pd	(%%r8),%%ymm12,%%ymm11			\n\t"\
	" vfmadd231pd	 0x060(%%r10),%%ymm3,%%ymm6			\n\t	 vfmadd231pd	 0x160(%%r10),%%ymm9 ,%%ymm14			\n\t"\
	"vfnmadd231pd	 0x060(%%r10),%%ymm2,%%ymm7			\n\t	vfnmadd231pd	 0x160(%%r10),%%ymm8 ,%%ymm15			\n\t"\
	"vfmadd132pd		  (%%rdi),%%ymm0,%%ymm5			\n\t	vfmadd132pd		(%%rdi),%%ymm10,%%ymm13			\n\t"\
	"vfmadd132pd		  (%%rdi),%%ymm1,%%ymm4			\n\t	vfmadd132pd		(%%rdi),%%ymm11,%%ymm12			\n\t"\
		"vmovaps		%%ymm7, 0x260(%%rbx)			\n\t		vmovaps		%%ymm15, 0x2e0(%%rbx)			\n\t"\
		"vmovaps		%%ymm6, 0x240(%%rbx)			\n\t		vmovaps		%%ymm14, 0x2c0(%%rbx)			\n\t"\
		"addq		$0x080,%%r10					\n\t"\
		"vmovaps		%%ymm5,%%ymm2					\n\t		vmovaps		%%ymm13,%%ymm8 					\n\t"\
		"vmovaps		%%ymm1,%%ymm3					\n\t		vmovaps		%%ymm11,%%ymm9 					\n\t"\
		"vmovaps		%%ymm0,%%ymm6					\n\t		vmovaps		%%ymm10,%%ymm14					\n\t"\
		"vmovaps		%%ymm4,%%ymm7					\n\t		vmovaps		%%ymm12,%%ymm15					\n\t"\
		"vmulpd			  (%%r10),%%ymm5,%%ymm5			\n\t		vmulpd		 0x100(%%r10),%%ymm13,%%ymm13			\n\t"\
		"vmulpd			  (%%r10),%%ymm1,%%ymm1			\n\t		vmulpd		 0x100(%%r10),%%ymm11,%%ymm11			\n\t"\
		"vmulpd		 0x040(%%r10),%%ymm0,%%ymm0			\n\t		vmulpd		 0x140(%%r10),%%ymm10,%%ymm10			\n\t"\
		"vmulpd		 0x040(%%r10),%%ymm4,%%ymm4			\n\t		vmulpd		 0x140(%%r10),%%ymm12,%%ymm12			\n\t"\
	" vfmadd231pd	 0x020(%%r10),%%ymm3,%%ymm5			\n\t	 vfmadd231pd	 0x120(%%r10),%%ymm9 ,%%ymm13			\n\t"\
	"vfnmadd231pd	 0x020(%%r10),%%ymm2,%%ymm1			\n\t	vfnmadd231pd	 0x120(%%r10),%%ymm8 ,%%ymm11			\n\t"\
	" vfmadd231pd	 0x060(%%r10),%%ymm7,%%ymm0			\n\t	 vfmadd231pd	 0x160(%%r10),%%ymm15,%%ymm10			\n\t"\
	"vfnmadd231pd	 0x060(%%r10),%%ymm6,%%ymm4			\n\t	vfnmadd231pd	 0x160(%%r10),%%ymm14,%%ymm12			\n\t"\
		"vmovaps		%%ymm1, 0x160(%%rbx)			\n\t		vmovaps		%%ymm11, 0x1e0(%%rbx)			\n\t"\
		"vmovaps		%%ymm5, 0x140(%%rbx)			\n\t		vmovaps		%%ymm13, 0x1c0(%%rbx)			\n\t"\
		"vmovaps		%%ymm4, 0x360(%%rbx)			\n\t		vmovaps		%%ymm12, 0x3e0(%%rbx)			\n\t"\
		"vmovaps		%%ymm0, 0x340(%%rbx)			\n\t		vmovaps		%%ymm10, 0x3c0(%%rbx)			\n\t"\
		"\n\t"\
	/************************************************/				/************************************************/\
	/*...Block 1: t00,t10,t20,t30 -> r00,04,02,06:  */				/*...Block 5: t08,t18,t28,t38 -> r08,0C,0A,0E:  */\
	/************************************************/				/************************************************/\
		"movq	%[__r00],%%rdx	/* base-addr in rcol = r08, so rdx+0x100 in rcol */	\n\t	vmovaps	(%%rsi),%%ymm10	\n\t"/* isrt2 */\
		"vmovaps			  (%%rdx),%%ymm0				\n\t		vmovaps		 0x140(%%rdx),%%ymm12				\n\t"\
		"vmovaps		 0x020(%%rdx),%%ymm1				\n\t		vmovaps		 0x160(%%rdx),%%ymm13				\n\t"\
		"vmovaps		 0x080(%%rdx),%%ymm2				\n\t		vmovaps		 0x1c0(%%rdx),%%ymm8 				\n\t"\
		"vmovaps		 0x0a0(%%rdx),%%ymm3				\n\t		vmovaps		 0x1e0(%%rdx),%%ymm9 				\n\t"\
	"vmovaps	(%%r8),%%ymm15	\n\t"/* 1.0 */\
		"vfnmadd231pd	0x080(%%rdx),%%ymm15,%%ymm0			\n\t		 vfmadd231pd	0x160(%%rdx),%%ymm15,%%ymm12	\n\t"\
		"vfnmadd231pd	0x0a0(%%rdx),%%ymm15,%%ymm1			\n\t		vfnmadd231pd	0x140(%%rdx),%%ymm15,%%ymm13	\n\t"\
		" vfmadd231pd		 (%%rdx),%%ymm15,%%ymm2			\n\t		vfnmadd231pd	0x1e0(%%rdx),%%ymm15,%%ymm8 	\n\t"\
		" vfmadd231pd	0x020(%%rdx),%%ymm15,%%ymm3			\n\t		 vfmadd231pd	0x1c0(%%rdx),%%ymm15,%%ymm9 	\n\t"\
		"vmovaps		 0x040(%%rdx),%%ymm4				\n\t		vmulpd			%%ymm10,%%ymm8 ,%%ymm8 			\n\t"\
		"vmovaps		 0x060(%%rdx),%%ymm5				\n\t		vmulpd			%%ymm10,%%ymm9 ,%%ymm9 			\n\t"\
		"vmovaps		 0x0c0(%%rdx),%%ymm6				\n\t		vmovaps		%%ymm12,%%ymm14					\n\t"\
		"vmovaps		 0x0e0(%%rdx),%%ymm7				\n\t		vmovaps		%%ymm13,%%ymm15					\n\t"\
		"vsubpd		 0x0c0(%%rdx),%%ymm4,%%ymm4			\n\t	vfmsub132pd		%%ymm10,%%ymm8 ,%%ymm12			\n\t"\
		"vsubpd		 0x0e0(%%rdx),%%ymm5,%%ymm5			\n\t	vfmsub132pd		%%ymm10,%%ymm9 ,%%ymm13			\n\t"\
		"vaddpd		 0x040(%%rdx),%%ymm6,%%ymm6			\n\t	vfmadd132pd		%%ymm10,%%ymm8 ,%%ymm14			\n\t"\
		"vaddpd		 0x060(%%rdx),%%ymm7,%%ymm7			\n\t	vfmadd132pd		%%ymm10,%%ymm9 ,%%ymm15			\n\t"\
		/* base-twiddle in l/rcol = c00/c04, so rcx+0x100 in rcol*/"	vmovaps		 0x100(%%rdx),%%ymm8 				\n\t"\
		"movq		%[__c10],%%rcx						\n\t		vmovaps		 0x120(%%rdx),%%ymm9 				\n\t"\
		"vfmadd132pd		(%%r8),%%ymm6,%%ymm2			\n\t		vmovaps		 0x180(%%rdx),%%ymm10				\n\t"\
		"vfmadd132pd		(%%r8),%%ymm7,%%ymm3			\n\t		vmovaps		 0x1a0(%%rdx),%%ymm11				\n\t"\
		"vmovaps		%%ymm2,      (%%rdx)			\n\t		vsubpd		 0x1a0(%%rdx),%%ymm8 ,%%ymm8 			\n\t"\
		"vmovaps		%%ymm3, 0x020(%%rdx)			\n\t		vsubpd		 0x180(%%rdx),%%ymm9 ,%%ymm9 			\n\t"\
	"vfnmadd231pd		(%%rdi),%%ymm6,%%ymm2				\n\t		vaddpd		 0x100(%%rdx),%%ymm11,%%ymm11			\n\t"\
	"vfnmadd231pd		(%%rdi),%%ymm7,%%ymm3				\n\t		vaddpd		 0x120(%%rdx),%%ymm10,%%ymm10			\n\t"\
		"vmovaps		%%ymm2,%%ymm6					\n\t		vfmsub132pd		(%%r8),%%ymm12,%%ymm11			\n\t"\
		"vmovaps		%%ymm3,%%ymm7					\n\t		vfmsub132pd		(%%r8),%%ymm13,%%ymm9 			\n\t"\
		"vmulpd			  (%%rcx),%%ymm2,%%ymm2	\n\t"/* c10 */"	vfmadd132pd		(%%rdi),%%ymm11,%%ymm12			\n\t"\
		"vmulpd			  (%%rcx),%%ymm3,%%ymm3			\n\t	vfmadd132pd		(%%rdi),%%ymm9 ,%%ymm13			\n\t"\
	" vfmadd231pd	 0x020(%%rcx),%%ymm7,%%ymm2			\n\t		vmovaps		%%ymm11, 0x140(%%rdx)			\n\t"\
	"vfnmadd231pd	 0x020(%%rcx),%%ymm6,%%ymm3			\n\t		vmovaps		%%ymm9 , 0x160(%%rdx)			\n\t"\
		"subq $0x40,%%rcx	/* put c00 in rcx to ease bookkeeping*/\n\t	vmovaps		%%ymm12,%%ymm11					\n\t"\
		/* add0,1 in rax,rbx; __r00 in rdx: */						"	vmovaps		%%ymm13,%%ymm9 					\n\t"\
		/* For each complex output octet, complex pairs having */	"	vmulpd		 0x100(%%rcx),%%ymm12,%%ymm12	/* c04 */\n\t"\
		/* reads from offsets 0x0..,0x1..,0x2..,0x3.. go into  */	"	vmulpd		 0x100(%%rcx),%%ymm13,%%ymm13			\n\t"\
		/* local-mem pairs rXY + 00/10, 04/14, 02/12, 06/16.   */	" vfmadd231pd	 0x120(%%rcx),%%ymm9 ,%%ymm12			\n\t"\
		/* For 1st octet we read from offsets [0x2..,0x0..],   */	"vfnmadd231pd	 0x120(%%rcx),%%ymm11,%%ymm13			\n\t"\
		/* [0x1..,0x3], other 3 octets use order [0,2],[1,3].  */\
		"vmovaps	0x220(%%rbx),%%ymm7						\n\t		vmovaps	0x0a0(%%rbx),%%ymm11						\n\t"\
		"vmovaps	0x200(%%rbx),%%ymm6						\n\t		vmovaps	0x080(%%rbx),%%ymm9 						\n\t"\
		"vmovaps	%%ymm3,0x060(%%rdx)						\n\t		vmovaps	%%ymm13,0x120(%%rdx)						\n\t"\
		"vmovaps	%%ymm2,0x040(%%rdx)		\n\t/* r02,03 */			vmovaps	%%ymm12,0x100(%%rdx)			\n\t"/* r08,09 */\
		"vmovaps	%%ymm7,0x260(%%rdx)						\n\t		vmovaps	%%ymm11,0x320(%%rdx)						\n\t"\
		"vmovaps	%%ymm6,0x240(%%rdx)		\n\t/* r12,13 */			vmovaps	%%ymm9 ,0x300(%%rdx)			\n\t"/* r18,19 */\
		"vmovaps		 0x020(%%rdx),%%ymm3				\n\t		vmovaps		0x140(%%rdx),%%ymm12				\n\t"\
		"vmovaps			  (%%rdx),%%ymm2				\n\t		vmovaps		0x160(%%rdx),%%ymm13				\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm7	\n\t vmovaps %%ymm12,%%ymm11 \n\t	vmulpd	 0x140(%%rcx),%%ymm12,%%ymm12	/* c14 */\n\t"\
		"vmovaps	0x000(%%rbx),%%ymm6	\n\t vmovaps %%ymm13,%%ymm9  \n\t	vmulpd	 0x140(%%rcx),%%ymm13,%%ymm13			\n\t"\
		"vmovaps	%%ymm3,0x020(%%rdx)						\n\t	 vfmadd231pd	 0x160(%%rcx),%%ymm9 ,%%ymm12			\n\t"\
		"vmovaps	%%ymm2,     (%%rdx)			/* r00,01 */\n\t	vfnmadd231pd	 0x160(%%rcx),%%ymm11,%%ymm13			\n\t"\
		"vmovaps	%%ymm7,0x220(%%rdx)						\n\t		vmovaps	0x2a0(%%rbx),%%ymm11						\n\t"\
		"vmovaps	%%ymm6,0x200(%%rdx)			/* r10,11 */\n\t		vmovaps	0x280(%%rbx),%%ymm9 						\n\t"\
		"vfmadd132pd		(%%r8),%%ymm5,%%ymm0			\n\t		vmovaps	%%ymm13,0x160(%%rdx)						\n\t"\
		"vfmsub132pd		(%%r8),%%ymm4,%%ymm1			\n\t		vmovaps	%%ymm12,0x140(%%rdx)			\n\t"/* r0a,0b */\
		"vmovaps		%%ymm0,%%ymm2					\n\t		vfmsub132pd		(%%r8),%%ymm15,%%ymm8 			\n\t"\
		"vmovaps		%%ymm1,%%ymm3					\n\t		vfmsub132pd		(%%r8),%%ymm14,%%ymm10			\n\t"\
		"vmovaps		%%ymm0,%%ymm6					\n\t		vmovaps	%%ymm11,0x360(%%rdx)						\n\t"\
		"vmovaps		%%ymm1,%%ymm7					\n\t		vmovaps	%%ymm9 ,0x340(%%rdx)			\n\t"/* r1a,1b */\
		"vfnmadd231pd		(%%rdi),%%ymm5,%%ymm0			\n\t	vfmadd132pd		(%%rdi),%%ymm8 ,%%ymm15			\n\t"\
		" vfmadd231pd		(%%rdi),%%ymm4,%%ymm1			\n\t	vfmadd132pd		(%%rdi),%%ymm10,%%ymm14			\n\t"\
		"vmulpd		 0x080(%%rcx),%%ymm2,%%ymm2	/* c08*/\n\t		vmovaps		%%ymm15,%%ymm12					\n\t"\
		"vmulpd		 0x080(%%rcx),%%ymm3,%%ymm3			\n\t		vmovaps		%%ymm10,%%ymm13					\n\t"\
	" vfmadd231pd	 0x0a0(%%rcx),%%ymm7,%%ymm2			\n\t		vmulpd		 0x180(%%rcx),%%ymm15,%%ymm15	/* c0C */\n\t"\
	"vfnmadd231pd	 0x0a0(%%rcx),%%ymm6,%%ymm3			\n\t		vmulpd		 0x180(%%rcx),%%ymm10,%%ymm10			\n\t"\
		"vmovaps	0x120(%%rbx),%%ymm7						\n\t	 vfmadd231pd	 0x1a0(%%rcx),%%ymm13,%%ymm15			\n\t"\
		"vmovaps	0x100(%%rbx),%%ymm6						\n\t	vfnmadd231pd	 0x1a0(%%rcx),%%ymm12,%%ymm10			\n\t"\
		"vmovaps	%%ymm3,0x0a0(%%rdx)						\n\t		vmovaps	0x1a0(%%rbx),%%ymm13						\n\t"\
		"vmovaps	%%ymm2,0x080(%%rdx)		\n\t/* r04,05 */\n\t		vmovaps	0x180(%%rbx),%%ymm12						\n\t"\
		"vmovaps	%%ymm7,0x2a0(%%rdx)						\n\t		vmovaps	%%ymm10,0x1a0(%%rdx)						\n\t"\
		"vmovaps	%%ymm6,0x280(%%rdx)		\n\t/* r14,15 */\n\t		vmovaps	%%ymm15,0x180(%%rdx)			\n\t"/* r0c,0d */\
		"													\n\t		vmovaps	%%ymm13,0x3a0(%%rdx)						\n\t"\
		"													\n\t		vmovaps	%%ymm12,0x380(%%rdx)			\n\t"/* r1c,1d */\
		"vmovaps		%%ymm0,%%ymm6					\n\t		vmovaps		%%ymm8 ,%%ymm12					\n\t"\
		"vmovaps		%%ymm1,%%ymm7					\n\t		vmovaps		%%ymm14,%%ymm13					\n\t"\
		"vmulpd		 0x0c0(%%rcx),%%ymm0,%%ymm0	/* c18*/\n\t		vmulpd		 0x1c0(%%rcx),%%ymm8 ,%%ymm8 	/* c1C */\n\t"\
		"vmulpd		 0x0c0(%%rcx),%%ymm1,%%ymm1			\n\t		vmulpd		 0x1c0(%%rcx),%%ymm14,%%ymm14			\n\t"\
	" vfmadd231pd	 0x0e0(%%rcx),%%ymm7,%%ymm0			\n\t	 vfmadd231pd	 0x1e0(%%rcx),%%ymm13,%%ymm8 			\n\t"\
	"vfnmadd231pd	 0x0e0(%%rcx),%%ymm6,%%ymm1			\n\t	vfnmadd231pd	 0x1e0(%%rcx),%%ymm12,%%ymm14			\n\t"\
		"vmovaps	0x320(%%rbx),%%ymm7						\n\t		vmovaps	0x3a0(%%rbx),%%ymm13						\n\t"\
		"vmovaps	0x300(%%rbx),%%ymm6						\n\t		vmovaps	0x380(%%rbx),%%ymm12						\n\t"\
		"vmovaps	%%ymm7,0x2e0(%%rdx)						\n\t		vmovaps	%%ymm13,0x3e0(%%rdx)						\n\t"\
		"vmovaps	%%ymm6,0x2c0(%%rdx)		/* r16,17 */	\n\t		vmovaps	%%ymm12,0x3c0(%%rdx)			\n\t"/* r1e,1f */\
		"vmovaps	%%ymm1,0x0e0(%%rdx)						\n\t		vmovaps	%%ymm14,0x1e0(%%rdx)						\n\t"\
		"vmovaps	%%ymm0,0x0c0(%%rdx)		/* r06,07 */	\n\t		vmovaps	%%ymm8 ,0x1c0(%%rdx)			\n\t"/* r0e,0f */\
		"\n\t"\
	/************************************************/				/************************************************/\
	/*...Block 3: t04,t14,t24,t34 -> r20,24,22,26:  */				/*...Block 7: t0C,t1C,t2C,t3C -> r28,2C,2A,2E:  */\
	/************************************************/				/************************************************/\
		"movq		%[__r20],%%rdx	\n\t"							/* base-addr in rcol = r28, so rdx offset +0x100 vs lcol */\
		"leaq	0x020(%%rsi),%%rcx	\n\t"/* cc0; Note cc0/ss0 are shared between lcol/rcol, so no rcx-offset until get to twiddles*/\
		"vmovaps		0x040(%%rdx),%%ymm4					\n\t		vmovaps		 0x140(%%rdx),%%ymm12			\n\t"\
		"vmovaps		0x0c0(%%rdx),%%ymm0					\n\t		vmovaps		 0x1c0(%%rdx),%%ymm8 			\n\t"\
		"vmovaps		0x060(%%rdx),%%ymm5					\n\t		vmovaps		 0x160(%%rdx),%%ymm13			\n\t"\
		"vmovaps		0x0e0(%%rdx),%%ymm3					\n\t		vmovaps		 0x1e0(%%rdx),%%ymm11			\n\t"\
		"vmovaps			%%ymm4 ,%%ymm6					\n\t		vmovaps		 	%%ymm12,%%ymm14				\n\t"\
		"vmovaps			%%ymm0 ,%%ymm2					\n\t		vmovaps		 	%%ymm8 ,%%ymm10				\n\t"\
		"vmovaps			%%ymm5 ,%%ymm7					\n\t		vmovaps		 	%%ymm13,%%ymm15				\n\t"\
	/*	"vmovaps			%%ymm1 ,%%ymm3					\n\t		vmovaps		 	%%ymm9 ,%%ymm11				\n\t"*/\
	/* Strategic register-usage for sincos consts here cuts loads from 16 to 4: */\
		"vmovaps			 (%%rcx),%%ymm9 			\n\t		vmovaps		0x020(%%rcx),%%ymm1 			\n\t"\
		"vmulpd		 	 %%ymm1 ,%%ymm7,%%ymm7				\n\t		vmulpd			 %%ymm9 ,%%ymm15,%%ymm15	\n\t"\
		"vmulpd		 	 %%ymm9 ,%%ymm3,%%ymm3				\n\t		vmulpd			 %%ymm1 ,%%ymm11,%%ymm11	\n\t"\
		"vmulpd		 	 %%ymm1 ,%%ymm6,%%ymm6				\n\t		vmulpd			 %%ymm9 ,%%ymm14,%%ymm14	\n\t"\
		"vmulpd		 	 %%ymm9 ,%%ymm2,%%ymm2				\n\t		vmulpd			 %%ymm1 ,%%ymm10,%%ymm10	\n\t"\
	"vfmadd132pd	 	 %%ymm9 ,%%ymm7,%%ymm4				\n\t	vfmadd132pd			 %%ymm1 ,%%ymm15,%%ymm12	\n\t"\
	"vfmadd132pd	 	 %%ymm1 ,%%ymm3,%%ymm0				\n\t	vfmadd132pd			 %%ymm9 ,%%ymm11,%%ymm8 	\n\t"\
	"vfmsub132pd	 	 %%ymm9 ,%%ymm6,%%ymm5				\n\t	vfmsub132pd			 %%ymm1 ,%%ymm14,%%ymm13	\n\t"\
	"vfmsub132pd	 0xe0(%%rdx),%%ymm2,%%ymm1				\n\t	vfmsub132pd		0x1e0(%%rdx),%%ymm10,%%ymm9 	\n\t"\
		"vmovaps		%%ymm5,%%ymm7					\n\t		vmovaps		%%ymm13,%%ymm15					\n\t"\
		"vmovaps		%%ymm4,%%ymm6					\n\t		vmovaps		%%ymm12,%%ymm14					\n\t"\
	"vmovaps	%%ymm15,0x40(%%rdx)	\n\t	vmovaps	(%%r8),%%ymm15	\n\t"/* spill to make room for 1.0 */\
		"vfmadd132pd	%%ymm15,%%ymm0,%%ymm4				\n\t		vfmadd132pd		%%ymm15,%%ymm8,%%ymm12			\n\t"\
		"vfmadd132pd	%%ymm15,%%ymm1,%%ymm5				\n\t		vfmadd132pd		%%ymm15,%%ymm9,%%ymm13			\n\t"\
		"vfmsub132pd	%%ymm15,%%ymm0,%%ymm6				\n\t		vfmsub132pd		%%ymm15,%%ymm8,%%ymm14			\n\t"\
		"vfmsub132pd	%%ymm15,%%ymm1,%%ymm7				\n\t		vfmsub132pd	0x40(%%rdx),%%ymm9,%%ymm15			\n\t"\
		"vmovaps		 0x080(%%rdx),%%ymm2				\n\t		vmovaps		 0x180(%%rdx),%%ymm10				\n\t"\
		"vmovaps		 0x0a0(%%rdx),%%ymm3				\n\t		vmovaps		 0x1a0(%%rdx),%%ymm11				\n\t"\
		"vmovaps			  (%%rdx),%%ymm0				\n\t		vmovaps		 0x100(%%rdx),%%ymm8 				\n\t"\
		"vmovaps		 0x020(%%rdx),%%ymm1				\n\t		vmovaps		 0x120(%%rdx),%%ymm9 				\n\t"\
		"vaddpd		 0x0a0(%%rdx),%%ymm2,%%ymm2			\n\t		vsubpd		 0x1a0(%%rdx),%%ymm10,%%ymm10			\n\t"\
		"vsubpd		 0x080(%%rdx),%%ymm3,%%ymm3			\n\t		vaddpd		 0x180(%%rdx),%%ymm11,%%ymm11			\n\t"\
		"vmulpd			  (%%rsi),%%ymm2,%%ymm2			\n\t		vmulpd			  (%%rsi),%%ymm10,%%ymm10			\n\t"\
		"vmulpd			  (%%rsi),%%ymm3,%%ymm3			\n\t		vmulpd			  (%%rsi),%%ymm11,%%ymm11			\n\t"\
		"movq		%[__c02],%%rcx				/* base-twiddle addr in rcol = c06, so rcx offset +0x100 vs lcol */	\n\t"\
	"vmovaps	%%ymm12,0x40(%%rdx)	\n\t	vmovaps	(%%r8) ,%%ymm12	\n\t"/* spill to make room for 1.0 */\
	"vmovaps	%%ymm13,0x60(%%rdx)	\n\t	vmovaps	(%%rdi),%%ymm13	\n\t"/* spill to make room for 2.0 */\
		"vfmsub132pd	%%ymm12,%%ymm2,%%ymm0				\n\t		vfmsub132pd	%%ymm12,%%ymm10,%%ymm8 			\n\t"\
		"vfmsub132pd	%%ymm12,%%ymm3,%%ymm1				\n\t		vfmsub132pd	%%ymm12,%%ymm11,%%ymm9 			\n\t"\
	"vfmadd132pd		%%ymm13,%%ymm0,%%ymm2			\n\t	vfmadd132pd		%%ymm13,%%ymm8 ,%%ymm10		\n\t"\
	"vfmadd132pd		%%ymm13,%%ymm1,%%ymm3			\n\t	vfmadd132pd		%%ymm13,%%ymm9 ,%%ymm11		\n\t"\
		"vfmsub132pd	%%ymm12,%%ymm4,%%ymm2				\n\t		vfmsub132pd	%%ymm12,%%ymm14,%%ymm8 			\n\t"\
		"vfmsub132pd	%%ymm12,%%ymm5,%%ymm3				\n\t		vfmsub132pd	%%ymm12,%%ymm15,%%ymm9 			\n\t"\
	"vfmadd132pd		%%ymm13,%%ymm2,%%ymm4			\n\t	vfmadd132pd		%%ymm13,%%ymm8 ,%%ymm14		\n\t"\
	"vfmadd132pd		%%ymm13,%%ymm3,%%ymm5			\n\t	vfmadd132pd		%%ymm13,%%ymm9 ,%%ymm15		\n\t"\
	"vmovaps	0x40(%%rdx),%%ymm12	\n\t"/* restore spill */\
	"vmovaps	0x60(%%rdx),%%ymm13	\n\t"/* restore spill */\
		"vmovaps		%%ymm2, 0x040(%%rdx)			\n\t		vmovaps		%%ymm8 , 0x140(%%rdx)			\n\t"\
		"vmovaps		%%ymm3, 0x060(%%rdx)			\n\t		vmovaps		%%ymm9 , 0x160(%%rdx)			\n\t"\
		"vmovaps		%%ymm4,%%ymm2					\n\t		vmovaps		%%ymm14,%%ymm8 					\n\t"\
		"vmovaps		%%ymm5,%%ymm3					\n\t		vmovaps		%%ymm15,%%ymm9 					\n\t"\
		"vmulpd			  (%%rcx),%%ymm4,%%ymm4			\n\t		vmulpd		 0x100(%%rcx),%%ymm14,%%ymm14			\n\t"\
		"vmulpd			  (%%rcx),%%ymm5,%%ymm5			\n\t		vmulpd		 0x100(%%rcx),%%ymm15,%%ymm15			\n\t"\
	" vfmadd231pd	 0x020(%%rcx),%%ymm3,%%ymm4			\n\t	 vfmadd231pd	 0x120(%%rcx),%%ymm9 ,%%ymm14			\n\t"\
	"vfnmadd231pd	 0x020(%%rcx),%%ymm2,%%ymm5			\n\t	vfnmadd231pd	 0x120(%%rcx),%%ymm8 ,%%ymm15			\n\t"\
		"vmovaps	0x060(%%rbx),%%ymm3						\n\t		vmovaps	0x0e0(%%rbx),%%ymm9 						\n\t"\
		"vmovaps	0x040(%%rbx),%%ymm2						\n\t		vmovaps	0x0c0(%%rbx),%%ymm8 						\n\t"\
		"vmovaps	%%ymm5,0x020(%%rdx)						\n\t		vmovaps	%%ymm15,0x120(%%rdx)						\n\t"\
		"vmovaps	%%ymm4,     (%%rdx)			\n\t/* r20,21 */		vmovaps	%%ymm14,0x100(%%rdx)			\n\t"/* r28,29 */\
		"vmovaps	%%ymm3,0x220(%%rdx)						\n\t		vmovaps	%%ymm9 ,0x320(%%rdx)						\n\t"\
		"vmovaps	%%ymm2,0x200(%%rdx)			\n\t/* r30,31 */		vmovaps	%%ymm8 ,0x300(%%rdx)			\n\t"/* r38,39 */\
		"movq		%[__c12],%%rcx					\n\t"		/* rcol uses c16 */\
		"vmovaps		 0x040(%%rdx),%%ymm4				\n\t		vmovaps		 0x140(%%rdx),%%ymm14				\n\t"\
		"vmovaps		 0x060(%%rdx),%%ymm5				\n\t		vmovaps		 0x160(%%rdx),%%ymm15				\n\t"\
		"vmovaps		%%ymm4,%%ymm2					\n\t		vmovaps		%%ymm14,%%ymm8 					\n\t"\
		"vmovaps		%%ymm5,%%ymm3					\n\t		vmovaps		%%ymm15,%%ymm9 					\n\t"\
		"vmulpd			  (%%rcx),%%ymm4,%%ymm4			\n\t		vmulpd		 0x100(%%rcx),%%ymm14,%%ymm14			\n\t"\
		"vmulpd			  (%%rcx),%%ymm5,%%ymm5			\n\t		vmulpd		 0x100(%%rcx),%%ymm15,%%ymm15			\n\t"\
	" vfmadd231pd	 0x020(%%rcx),%%ymm3,%%ymm4			\n\t	 vfmadd231pd	 0x120(%%rcx),%%ymm9 ,%%ymm14			\n\t"\
	"vfnmadd231pd	 0x020(%%rcx),%%ymm2,%%ymm5			\n\t	vfnmadd231pd	 0x120(%%rcx),%%ymm8 ,%%ymm15			\n\t"\
		"vmovaps	0x260(%%rbx),%%ymm3						\n\t		vmovaps	0x2e0(%%rbx),%%ymm9 						\n\t"\
		"vmovaps	0x240(%%rbx),%%ymm2						\n\t		vmovaps	0x2c0(%%rbx),%%ymm8 						\n\t"\
		"vmovaps	%%ymm5,0x060(%%rdx)						\n\t		vmovaps	%%ymm15,0x160(%%rdx)						\n\t"\
		"vmovaps	%%ymm4,0x040(%%rdx)			\n\t/* r22,23 */		vmovaps	%%ymm14,0x140(%%rdx)			\n\t"/* r2a,2b */\
		"vmovaps	%%ymm3,0x260(%%rdx)						\n\t		vmovaps	%%ymm9 ,0x360(%%rdx)						\n\t"\
		"vmovaps	%%ymm2,0x240(%%rdx)			\n\t/* r32,33 */		vmovaps	%%ymm8 ,0x340(%%rdx)			\n\t"/* r3a,3b */\
		"movq		%[__c0A],%%rcx					\n\t"		/* rcol uses c0E */\
	"vfmsub132pd		(%%r8) ,%%ymm7,%%ymm0			\n\t	vfmsub132pd		(%%r8) ,%%ymm13,%%ymm10			\n\t"\
	"vfmsub132pd		(%%r8) ,%%ymm6,%%ymm1			\n\t	vfmsub132pd		(%%r8) ,%%ymm12,%%ymm11			\n\t"\
	"vfmadd132pd		(%%rdi),%%ymm0,%%ymm7			\n\t	vfmadd132pd		(%%rdi),%%ymm10,%%ymm13			\n\t"\
	"vfmadd132pd		(%%rdi),%%ymm1,%%ymm6			\n\t	vfmadd132pd		(%%rdi),%%ymm11,%%ymm12			\n\t"\
		"vmovaps		%%ymm7,%%ymm4					\n\t		vmovaps		%%ymm13,%%ymm8 					\n\t"\
		"vmovaps		%%ymm1,%%ymm5					\n\t		vmovaps		%%ymm11,%%ymm9 					\n\t"\
		"vmulpd			  (%%rcx),%%ymm7,%%ymm7			\n\t		vmulpd		 0x100(%%rcx),%%ymm13,%%ymm13			\n\t"\
		"vmulpd			  (%%rcx),%%ymm1,%%ymm1			\n\t		vmulpd		 0x100(%%rcx),%%ymm11,%%ymm11			\n\t"\
	" vfmadd231pd	 0x020(%%rcx),%%ymm5,%%ymm7			\n\t	 vfmadd231pd	 0x120(%%rcx),%%ymm9 ,%%ymm13			\n\t"\
	"vfnmadd231pd	 0x020(%%rcx),%%ymm4,%%ymm1			\n\t	vfnmadd231pd	 0x120(%%rcx),%%ymm8 ,%%ymm11			\n\t"\
		"vmovaps	0x160(%%rbx),%%ymm5						\n\t		vmovaps	0x1e0(%%rbx),%%ymm9 						\n\t"\
		"vmovaps	0x140(%%rbx),%%ymm4						\n\t		vmovaps	0x1c0(%%rbx),%%ymm8 						\n\t"\
		"vmovaps	%%ymm1,0x0a0(%%rdx)						\n\t		vmovaps	%%ymm11,0x1a0(%%rdx)						\n\t"\
		"vmovaps	%%ymm7,0x080(%%rdx)			\n\t/* r24,25 */		vmovaps	%%ymm13,0x180(%%rdx)			\n\t"/* r2c,2d */\
		"vmovaps	%%ymm5,0x2a0(%%rdx)						\n\t		vmovaps	%%ymm9 ,0x3a0(%%rdx)						\n\t"\
		"vmovaps	%%ymm4,0x280(%%rdx)			\n\t/* r34,35 */		vmovaps	%%ymm8 ,0x380(%%rdx)			\n\t"/* r3c,3d */\
		"movq		%[__c1A],%%rcx					\n\t"		/* rcol uses c1E */\
		"vmovaps		%%ymm0,%%ymm4					\n\t		vmovaps		%%ymm10,%%ymm8 					\n\t"\
		"vmovaps		%%ymm6,%%ymm5					\n\t		vmovaps		%%ymm12,%%ymm9 					\n\t"\
		"vmulpd			  (%%rcx),%%ymm0,%%ymm0			\n\t		vmulpd		 0x100(%%rcx),%%ymm10,%%ymm10			\n\t"\
		"vmulpd			  (%%rcx),%%ymm6,%%ymm6			\n\t		vmulpd		 0x100(%%rcx),%%ymm12,%%ymm12			\n\t"\
	" vfmadd231pd	 0x020(%%rcx),%%ymm5,%%ymm0			\n\t	 vfmadd231pd	 0x120(%%rcx),%%ymm9 ,%%ymm10			\n\t"\
	"vfnmadd231pd	 0x020(%%rcx),%%ymm4,%%ymm6			\n\t	vfnmadd231pd	 0x120(%%rcx),%%ymm8 ,%%ymm12			\n\t"\
		"vmovaps	0x360(%%rbx),%%ymm5						\n\t		vmovaps	0x3e0(%%rbx),%%ymm9 						\n\t"\
		"vmovaps	0x340(%%rbx),%%ymm4						\n\t		vmovaps	0x3c0(%%rbx),%%ymm8 						\n\t"\
		"vmovaps	%%ymm6,0x0e0(%%rdx)						\n\t		vmovaps	%%ymm12,0x1e0(%%rdx)						\n\t"\
		"vmovaps	%%ymm0,0x0c0(%%rdx)			\n\t/* r26,27 */		vmovaps	%%ymm10,0x1c0(%%rdx)			\n\t"/* r2e,2f */\
		"vmovaps	%%ymm5,0x2e0(%%rdx)						\n\t		vmovaps	%%ymm9 ,0x3e0(%%rdx)						\n\t"\
		"vmovaps	%%ymm4,0x2c0(%%rdx)			\n\t/* r36,37 */		vmovaps	%%ymm8 ,0x3c0(%%rdx)			\n\t"/* r3e,3f */\
	/*******************************************/\
	/**** Finish with 4-way 'un'terleaving: ****/\
	/*******************************************/\
		"movq	%[__r00] ,%%rsi\n\t"\
		"movq	%[__add0],%%rax\n\t"\
		"movq	%[__add1],%%rbx\n\t"\
		"movq	%[__add2],%%rcx\n\t"\
		"movq	%[__add3],%%rdx\n\t"\
	"/* a[j+p0]: Inputs from r00 +0/1, 8/9, 16/17, 24/25. Outputs into [add0,add0+0x100,add1,add1+0x100]+0x0: */	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps		 (%%rsi),%%ymm1						\n\t		vmovaps	0x020(%%rsi),%%ymm3								\n\t"\
		"vmovaps	0x200(%%rsi),%%ymm5						\n\t		vmovaps	0x220(%%rsi),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps	0x400(%%rsi),%%ymm6						\n\t		vmovaps	0x420(%%rsi),%%ymm14							\n\t"\
		"vmovaps	0x600(%%rsi),%%ymm5						\n\t		vmovaps	0x620(%%rsi),%%ymm13							\n\t"\
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
	"/* a[j+p2]: Inputs from r08 +0/1, 8/9, 16/17, 24/25. Outputs into [add0,add0+0x100,add1,add1+0x100]+0x0: */	\n\t"\
		"addq	$0x100,%%rsi	\n\t"\
		"addq	$0x40,%%rax	\n\t"\
		"addq	$0x40,%%rbx	\n\t"\
		"addq	$0x40,%%rcx	\n\t"\
		"addq	$0x40,%%rdx	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps		 (%%rsi),%%ymm1						\n\t		vmovaps	0x020(%%rsi),%%ymm3								\n\t"\
		"vmovaps	0x200(%%rsi),%%ymm5						\n\t		vmovaps	0x220(%%rsi),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps	0x400(%%rsi),%%ymm6						\n\t		vmovaps	0x420(%%rsi),%%ymm14							\n\t"\
		"vmovaps	0x600(%%rsi),%%ymm5						\n\t		vmovaps	0x620(%%rsi),%%ymm13							\n\t"\
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
	"/* a[j+p4]: Inputs from r04 +0/1, 8/9, 16/17, 24/25. Outputs into [add0,add0+0x100,add1,add1+0x100]+0x0: */	\n\t"\
		"subq	$0x80,%%rsi	\n\t"\
		"addq	$0x40,%%rax	\n\t"\
		"addq	$0x40,%%rbx	\n\t"\
		"addq	$0x40,%%rcx	\n\t"\
		"addq	$0x40,%%rdx	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps		 (%%rsi),%%ymm1						\n\t		vmovaps	0x020(%%rsi),%%ymm3								\n\t"\
		"vmovaps	0x200(%%rsi),%%ymm5						\n\t		vmovaps	0x220(%%rsi),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps	0x400(%%rsi),%%ymm6						\n\t		vmovaps	0x420(%%rsi),%%ymm14							\n\t"\
		"vmovaps	0x600(%%rsi),%%ymm5						\n\t		vmovaps	0x620(%%rsi),%%ymm13							\n\t"\
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
	"/* a[j+p6]: Inputs from r0c +0/1, 8/9, 16/17, 24/25. Outputs into [add0,add0+0x100,add1,add1+0x100]+0x0: */	\n\t"\
		"addq	$0x100,%%rsi	\n\t"\
		"addq	$0x40,%%rax	\n\t"\
		"addq	$0x40,%%rbx	\n\t"\
		"addq	$0x40,%%rcx	\n\t"\
		"addq	$0x40,%%rdx	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps		 (%%rsi),%%ymm1						\n\t		vmovaps	0x020(%%rsi),%%ymm3								\n\t"\
		"vmovaps	0x200(%%rsi),%%ymm5						\n\t		vmovaps	0x220(%%rsi),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps	0x400(%%rsi),%%ymm6						\n\t		vmovaps	0x420(%%rsi),%%ymm14							\n\t"\
		"vmovaps	0x600(%%rsi),%%ymm5						\n\t		vmovaps	0x620(%%rsi),%%ymm13							\n\t"\
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
	"/* a[j+p8]: Inputs from r02 +0/1, 8/9, 16/17, 24/25. Outputs into [add0,add0+0x100,add1,add1+0x100]+0x0: */	\n\t"\
		"subq	$0x140,%%rsi	\n\t"\
		"addq	$0x40,%%rax	\n\t"\
		"addq	$0x40,%%rbx	\n\t"\
		"addq	$0x40,%%rcx	\n\t"\
		"addq	$0x40,%%rdx	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps		 (%%rsi),%%ymm1						\n\t		vmovaps	0x020(%%rsi),%%ymm3								\n\t"\
		"vmovaps	0x200(%%rsi),%%ymm5						\n\t		vmovaps	0x220(%%rsi),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps	0x400(%%rsi),%%ymm6						\n\t		vmovaps	0x420(%%rsi),%%ymm14							\n\t"\
		"vmovaps	0x600(%%rsi),%%ymm5						\n\t		vmovaps	0x620(%%rsi),%%ymm13							\n\t"\
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
	"/* a[j+p10]: Inputs from r0a +0/1, 8/9, 16/17, 24/25. Outputs into [add0,add0+0x100,add1,add1+0x100]+0x0: */	\n\t"\
		"addq	$0x100,%%rsi	\n\t"\
		"addq	$0x40,%%rax	\n\t"\
		"addq	$0x40,%%rbx	\n\t"\
		"addq	$0x40,%%rcx	\n\t"\
		"addq	$0x40,%%rdx	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps		 (%%rsi),%%ymm1						\n\t		vmovaps	0x020(%%rsi),%%ymm3								\n\t"\
		"vmovaps	0x200(%%rsi),%%ymm5						\n\t		vmovaps	0x220(%%rsi),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps	0x400(%%rsi),%%ymm6						\n\t		vmovaps	0x420(%%rsi),%%ymm14							\n\t"\
		"vmovaps	0x600(%%rsi),%%ymm5						\n\t		vmovaps	0x620(%%rsi),%%ymm13							\n\t"\
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
	"/* a[j+p12]: Inputs from r06 +0/1, 8/9, 16/17, 24/25. Outputs into [add0,add0+0x100,add1,add1+0x100]+0x0: */	\n\t"\
		"subq	$0x80,%%rsi	\n\t"\
		"addq	$0x40,%%rax	\n\t"\
		"addq	$0x40,%%rbx	\n\t"\
		"addq	$0x40,%%rcx	\n\t"\
		"addq	$0x40,%%rdx	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps		 (%%rsi),%%ymm1						\n\t		vmovaps	0x020(%%rsi),%%ymm3								\n\t"\
		"vmovaps	0x200(%%rsi),%%ymm5						\n\t		vmovaps	0x220(%%rsi),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps	0x400(%%rsi),%%ymm6						\n\t		vmovaps	0x420(%%rsi),%%ymm14							\n\t"\
		"vmovaps	0x600(%%rsi),%%ymm5						\n\t		vmovaps	0x620(%%rsi),%%ymm13							\n\t"\
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
	"/* a[j+p14]: Inputs from r0e +0/1, 8/9, 16/17, 24/25. Outputs into [add0,add0+0x100,add1,add1+0x100]+0x0: */	\n\t"\
		"addq	$0x100,%%rsi	\n\t"\
		"addq	$0x40,%%rax	\n\t"\
		"addq	$0x40,%%rbx	\n\t"\
		"addq	$0x40,%%rcx	\n\t"\
		"addq	$0x40,%%rdx	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps		 (%%rsi),%%ymm1						\n\t		vmovaps	0x020(%%rsi),%%ymm3								\n\t"\
		"vmovaps	0x200(%%rsi),%%ymm5						\n\t		vmovaps	0x220(%%rsi),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps	0x400(%%rsi),%%ymm6						\n\t		vmovaps	0x420(%%rsi),%%ymm14							\n\t"\
		"vmovaps	0x600(%%rsi),%%ymm5						\n\t		vmovaps	0x620(%%rsi),%%ymm13							\n\t"\
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
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__r00] "m" (Xr00)\
		 ,[__r08] "m" (Xr08)\
		 ,[__r10] "m" (Xr10)\
		 ,[__r20] "m" (Xr20)\
		 ,[__r28] "m" (Xr28)\
		 ,[__r30] "m" (Xr30)\
		 ,[__c01] "m" (Xc01)\
		 ,[__c02] "m" (Xc02)\
		 ,[__c04] "m" (Xc04)\
		 ,[__c06] "m" (Xc06)\
		 ,[__c08] "m" (Xc08)\
		 ,[__c0A] "m" (Xc0A)\
		 ,[__c0C] "m" (Xc0C)\
		 ,[__c0E] "m" (Xc0E)\
		 ,[__c10] "m" (Xc10)\
		 ,[__c12] "m" (Xc12)\
		 ,[__c14] "m" (Xc14)\
		 ,[__c16] "m" (Xc16)\
		 ,[__c18] "m" (Xc18)\
		 ,[__c1A] "m" (Xc1A)\
		 ,[__c1C] "m" (Xc1C)\
		 ,[__c1E] "m" (Xc1E)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r8","r10","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"		/* Clobbered registers */\
	);\
	}

  #else	// ALL_FMA = False: All FMAs have non-unity multiplicands:

	// FMA version: replace [222 ADD, 220 SUB, 322 MUL, 840 memref] ==> [58 ADD, 158 SUB, 98 MUL, 232 FMA, 844 memref].
	//
	#define SSE2_RADIX32_WRAPPER_DIF(Xadd0,Xadd1,Xadd2,Xadd3,Xr00,Xr10,Xr20,Xr30,Xisrt2,Xcc0,Xc00,Xc01,Xc02,Xc03,Xc05,Xc07)\
	{\
	__asm__ volatile (\
		"movq	%[__r00] ,%%rsi	\n\t"\
		"movq	%[__add0],%%rax\n\t"\
		"movq	%[__add1],%%rbx\n\t"\
		"movq	%[__add2],%%rcx\n\t"\
		"movq	%[__add3],%%rdx\n\t"\
	"/**** Start with 4-way interleaving: ****/\n\t"\
	"/* a[j+p0]: Inputs from [add0,add0+0x200,add1,add1+0x200]+0x0. Outputs into r0 +0/1, 16/17, 32/33, 48/49: */	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps		 (%%rax),%%ymm1						\n\t		vmovaps	0x020(%%rax),%%ymm3								\n\t"\
		"vmovaps		 (%%rbx),%%ymm5						\n\t		vmovaps	0x020(%%rbx),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps		 (%%rcx),%%ymm6						\n\t		vmovaps	0x020(%%rcx),%%ymm14							\n\t"\
		"vmovaps		 (%%rdx),%%ymm5						\n\t		vmovaps	0x020(%%rdx),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm6 ,%%ymm0				\n\t		vshufpd	$15 ,%%ymm13,%%ymm14 ,%%ymm2					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm6 ,%%ymm6				\n\t		vshufpd	$0  ,%%ymm13,%%ymm14 ,%%ymm14					\n\t"\
		"vperm2f128 $32 ,%%ymm0 ,%%ymm7 ,%%ymm5				\n\t		vperm2f128 $32 ,%%ymm2 ,%%ymm15,%%ymm13		/* outB	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm0 ,%%ymm7 ,%%ymm7				\n\t		vperm2f128 $49 ,%%ymm2 ,%%ymm15,%%ymm15		/* outD	*/	\n\t"\
		"vperm2f128 $32 ,%%ymm6 ,%%ymm1 ,%%ymm0				\n\t		vperm2f128 $32 ,%%ymm14,%%ymm3 ,%%ymm2 		/* outA	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm6 ,%%ymm1 ,%%ymm1				\n\t		vperm2f128 $49 ,%%ymm14,%%ymm3 ,%%ymm3 		/* outC	*/	\n\t"\
		"vmovaps 	%%ymm5 ,0x400(%%rsi)					\n\t		vmovaps %%ymm13,0x420(%%rsi)				/* outB	*/	\n\t"\
		"vmovaps 	%%ymm7 ,0x600(%%rsi)					\n\t		vmovaps %%ymm15,0x620(%%rsi)				/* outD	*/	\n\t"\
		"vmovaps 	%%ymm0 ,     (%%rsi)					\n\t		vmovaps %%ymm2 ,0x020(%%rsi)				/* outA	*/	\n\t"\
		"vmovaps 	%%ymm1 ,0x200(%%rsi)					\n\t		vmovaps %%ymm3 ,0x220(%%rsi)				/* outC	*/	\n\t"\
	"/* a[j+p2]: Inputs from [add0,add0+0x200,add1,add1+0x200]+0x40. Outputs into **r8** +0/1, 16/17, 32/33, 48/49: */	\n\t"\
		"addq	$0x40,%%rax	\n\t"\
		"addq	$0x40,%%rbx	\n\t"\
		"addq	$0x40,%%rcx	\n\t"\
		"addq	$0x40,%%rdx	\n\t"\
		"addq	$0x100,%%rsi	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps		 (%%rax),%%ymm1						\n\t		vmovaps	0x020(%%rax),%%ymm3								\n\t"\
		"vmovaps		 (%%rbx),%%ymm5						\n\t		vmovaps	0x020(%%rbx),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps		 (%%rcx),%%ymm6						\n\t		vmovaps	0x020(%%rcx),%%ymm14							\n\t"\
		"vmovaps		 (%%rdx),%%ymm5						\n\t		vmovaps	0x020(%%rdx),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm6 ,%%ymm0				\n\t		vshufpd	$15 ,%%ymm13,%%ymm14 ,%%ymm2					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm6 ,%%ymm6				\n\t		vshufpd	$0  ,%%ymm13,%%ymm14 ,%%ymm14					\n\t"\
		"vperm2f128 $32 ,%%ymm0 ,%%ymm7 ,%%ymm5				\n\t		vperm2f128 $32 ,%%ymm2 ,%%ymm15,%%ymm13		/* outB	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm0 ,%%ymm7 ,%%ymm7				\n\t		vperm2f128 $49 ,%%ymm2 ,%%ymm15,%%ymm15		/* outD	*/	\n\t"\
		"vperm2f128 $32 ,%%ymm6 ,%%ymm1 ,%%ymm0				\n\t		vperm2f128 $32 ,%%ymm14,%%ymm3 ,%%ymm2 		/* outA	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm6 ,%%ymm1 ,%%ymm1				\n\t		vperm2f128 $49 ,%%ymm14,%%ymm3 ,%%ymm3 		/* outC	*/	\n\t"\
		"vmovaps 	%%ymm5 ,0x400(%%rsi)					\n\t		vmovaps %%ymm13,0x420(%%rsi)				/* outB	*/	\n\t"\
		"vmovaps 	%%ymm7 ,0x600(%%rsi)					\n\t		vmovaps %%ymm15,0x620(%%rsi)				/* outD	*/	\n\t"\
		"vmovaps 	%%ymm0 ,     (%%rsi)					\n\t		vmovaps %%ymm2 ,0x020(%%rsi)				/* outA	*/	\n\t"\
		"vmovaps 	%%ymm1 ,0x200(%%rsi)					\n\t		vmovaps %%ymm3 ,0x220(%%rsi)				/* outC	*/	\n\t"\
	"/* a[j+p4]: Inputs from [add0,add0+0x200,add1,add1+0x200]+0x80. Outputs into r4 +0/1, 16/17, 32/33, 48/49: */	\n\t"\
		"addq	$0x40,%%rax	\n\t"\
		"addq	$0x40,%%rbx	\n\t"\
		"addq	$0x40,%%rcx	\n\t"\
		"addq	$0x40,%%rdx	\n\t"\
		"subq	$0x80,%%rsi	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps		 (%%rax),%%ymm1						\n\t		vmovaps	0x020(%%rax),%%ymm3								\n\t"\
		"vmovaps		 (%%rbx),%%ymm5						\n\t		vmovaps	0x020(%%rbx),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps		 (%%rcx),%%ymm6						\n\t		vmovaps	0x020(%%rcx),%%ymm14							\n\t"\
		"vmovaps		 (%%rdx),%%ymm5						\n\t		vmovaps	0x020(%%rdx),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm6 ,%%ymm0				\n\t		vshufpd	$15 ,%%ymm13,%%ymm14 ,%%ymm2					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm6 ,%%ymm6				\n\t		vshufpd	$0  ,%%ymm13,%%ymm14 ,%%ymm14					\n\t"\
		"vperm2f128 $32 ,%%ymm0 ,%%ymm7 ,%%ymm5				\n\t		vperm2f128 $32 ,%%ymm2 ,%%ymm15,%%ymm13		/* outB	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm0 ,%%ymm7 ,%%ymm7				\n\t		vperm2f128 $49 ,%%ymm2 ,%%ymm15,%%ymm15		/* outD	*/	\n\t"\
		"vperm2f128 $32 ,%%ymm6 ,%%ymm1 ,%%ymm0				\n\t		vperm2f128 $32 ,%%ymm14,%%ymm3 ,%%ymm2 		/* outA	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm6 ,%%ymm1 ,%%ymm1				\n\t		vperm2f128 $49 ,%%ymm14,%%ymm3 ,%%ymm3 		/* outC	*/	\n\t"\
		"vmovaps 	%%ymm5 ,0x400(%%rsi)					\n\t		vmovaps %%ymm13,0x420(%%rsi)				/* outB	*/	\n\t"\
		"vmovaps 	%%ymm7 ,0x600(%%rsi)					\n\t		vmovaps %%ymm15,0x620(%%rsi)				/* outD	*/	\n\t"\
		"vmovaps 	%%ymm0 ,     (%%rsi)					\n\t		vmovaps %%ymm2 ,0x020(%%rsi)				/* outA	*/	\n\t"\
		"vmovaps 	%%ymm1 ,0x200(%%rsi)					\n\t		vmovaps %%ymm3 ,0x220(%%rsi)				/* outC	*/	\n\t"\
	"/* a[j+p6]: Inputs from [add0,add0+0x200,add1,add1+0x200]+0xc0. Outputs into **r12** +0/1, 16/17, 32/33, 48/49: */	\n\t"\
		"addq	$0x40,%%rax	\n\t"\
		"addq	$0x40,%%rbx	\n\t"\
		"addq	$0x40,%%rcx	\n\t"\
		"addq	$0x40,%%rdx	\n\t"\
		"addq	$0x100,%%rsi	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps		 (%%rax),%%ymm1						\n\t		vmovaps	0x020(%%rax),%%ymm3								\n\t"\
		"vmovaps		 (%%rbx),%%ymm5						\n\t		vmovaps	0x020(%%rbx),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps		 (%%rcx),%%ymm6						\n\t		vmovaps	0x020(%%rcx),%%ymm14							\n\t"\
		"vmovaps		 (%%rdx),%%ymm5						\n\t		vmovaps	0x020(%%rdx),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm6 ,%%ymm0				\n\t		vshufpd	$15 ,%%ymm13,%%ymm14 ,%%ymm2					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm6 ,%%ymm6				\n\t		vshufpd	$0  ,%%ymm13,%%ymm14 ,%%ymm14					\n\t"\
		"vperm2f128 $32 ,%%ymm0 ,%%ymm7 ,%%ymm5				\n\t		vperm2f128 $32 ,%%ymm2 ,%%ymm15,%%ymm13		/* outB	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm0 ,%%ymm7 ,%%ymm7				\n\t		vperm2f128 $49 ,%%ymm2 ,%%ymm15,%%ymm15		/* outD	*/	\n\t"\
		"vperm2f128 $32 ,%%ymm6 ,%%ymm1 ,%%ymm0				\n\t		vperm2f128 $32 ,%%ymm14,%%ymm3 ,%%ymm2 		/* outA	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm6 ,%%ymm1 ,%%ymm1				\n\t		vperm2f128 $49 ,%%ymm14,%%ymm3 ,%%ymm3 		/* outC	*/	\n\t"\
		"vmovaps 	%%ymm5 ,0x400(%%rsi)					\n\t		vmovaps %%ymm13,0x420(%%rsi)				/* outB	*/	\n\t"\
		"vmovaps 	%%ymm7 ,0x600(%%rsi)					\n\t		vmovaps %%ymm15,0x620(%%rsi)				/* outD	*/	\n\t"\
		"vmovaps 	%%ymm0 ,     (%%rsi)					\n\t		vmovaps %%ymm2 ,0x020(%%rsi)				/* outA	*/	\n\t"\
		"vmovaps 	%%ymm1 ,0x200(%%rsi)					\n\t		vmovaps %%ymm3 ,0x220(%%rsi)				/* outC	*/	\n\t"\
	"/* a[j+p8]: Inputs from [add0,add0+0x200,add1,add1+0x200]+0x100. Outputs into **r2** +0/1, 16/17, 32/33, 48/49: */	\n\t"\
		"addq	$0x40,%%rax	\n\t"\
		"addq	$0x40,%%rbx	\n\t"\
		"addq	$0x40,%%rcx	\n\t"\
		"addq	$0x40,%%rdx	\n\t"\
		"subq	$0x140,%%rsi	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps		 (%%rax),%%ymm1						\n\t		vmovaps	0x020(%%rax),%%ymm3								\n\t"\
		"vmovaps		 (%%rbx),%%ymm5						\n\t		vmovaps	0x020(%%rbx),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps		 (%%rcx),%%ymm6						\n\t		vmovaps	0x020(%%rcx),%%ymm14							\n\t"\
		"vmovaps		 (%%rdx),%%ymm5						\n\t		vmovaps	0x020(%%rdx),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm6 ,%%ymm0				\n\t		vshufpd	$15 ,%%ymm13,%%ymm14 ,%%ymm2					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm6 ,%%ymm6				\n\t		vshufpd	$0  ,%%ymm13,%%ymm14 ,%%ymm14					\n\t"\
		"vperm2f128 $32 ,%%ymm0 ,%%ymm7 ,%%ymm5				\n\t		vperm2f128 $32 ,%%ymm2 ,%%ymm15,%%ymm13		/* outB	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm0 ,%%ymm7 ,%%ymm7				\n\t		vperm2f128 $49 ,%%ymm2 ,%%ymm15,%%ymm15		/* outD	*/	\n\t"\
		"vperm2f128 $32 ,%%ymm6 ,%%ymm1 ,%%ymm0				\n\t		vperm2f128 $32 ,%%ymm14,%%ymm3 ,%%ymm2 		/* outA	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm6 ,%%ymm1 ,%%ymm1				\n\t		vperm2f128 $49 ,%%ymm14,%%ymm3 ,%%ymm3 		/* outC	*/	\n\t"\
		"vmovaps 	%%ymm5 ,0x400(%%rsi)					\n\t		vmovaps %%ymm13,0x420(%%rsi)				/* outB	*/	\n\t"\
		"vmovaps 	%%ymm7 ,0x600(%%rsi)					\n\t		vmovaps %%ymm15,0x620(%%rsi)				/* outD	*/	\n\t"\
		"vmovaps 	%%ymm0 ,     (%%rsi)					\n\t		vmovaps %%ymm2 ,0x020(%%rsi)				/* outA	*/	\n\t"\
		"vmovaps 	%%ymm1 ,0x200(%%rsi)					\n\t		vmovaps %%ymm3 ,0x220(%%rsi)				/* outC	*/	\n\t"\
	"/* a[j+p10]: Inputs from [add0,add0+0x200,add1,add1+0x200]+0x140. Outputs into r10 +0/1, 16/17, 32/33, 48/49: */	\n\t"\
		"addq	$0x40,%%rax	\n\t"\
		"addq	$0x40,%%rbx	\n\t"\
		"addq	$0x40,%%rcx	\n\t"\
		"addq	$0x40,%%rdx	\n\t"\
		"addq	$0x100,%%rsi	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps		 (%%rax),%%ymm1						\n\t		vmovaps	0x020(%%rax),%%ymm3								\n\t"\
		"vmovaps		 (%%rbx),%%ymm5						\n\t		vmovaps	0x020(%%rbx),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps		 (%%rcx),%%ymm6						\n\t		vmovaps	0x020(%%rcx),%%ymm14							\n\t"\
		"vmovaps		 (%%rdx),%%ymm5						\n\t		vmovaps	0x020(%%rdx),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm6 ,%%ymm0				\n\t		vshufpd	$15 ,%%ymm13,%%ymm14 ,%%ymm2					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm6 ,%%ymm6				\n\t		vshufpd	$0  ,%%ymm13,%%ymm14 ,%%ymm14					\n\t"\
		"vperm2f128 $32 ,%%ymm0 ,%%ymm7 ,%%ymm5				\n\t		vperm2f128 $32 ,%%ymm2 ,%%ymm15,%%ymm13		/* outB	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm0 ,%%ymm7 ,%%ymm7				\n\t		vperm2f128 $49 ,%%ymm2 ,%%ymm15,%%ymm15		/* outD	*/	\n\t"\
		"vperm2f128 $32 ,%%ymm6 ,%%ymm1 ,%%ymm0				\n\t		vperm2f128 $32 ,%%ymm14,%%ymm3 ,%%ymm2 		/* outA	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm6 ,%%ymm1 ,%%ymm1				\n\t		vperm2f128 $49 ,%%ymm14,%%ymm3 ,%%ymm3 		/* outC	*/	\n\t"\
		"vmovaps 	%%ymm5 ,0x400(%%rsi)					\n\t		vmovaps %%ymm13,0x420(%%rsi)				/* outB	*/	\n\t"\
		"vmovaps 	%%ymm7 ,0x600(%%rsi)					\n\t		vmovaps %%ymm15,0x620(%%rsi)				/* outD	*/	\n\t"\
		"vmovaps 	%%ymm0 ,     (%%rsi)					\n\t		vmovaps %%ymm2 ,0x020(%%rsi)				/* outA	*/	\n\t"\
		"vmovaps 	%%ymm1 ,0x200(%%rsi)					\n\t		vmovaps %%ymm3 ,0x220(%%rsi)				/* outC	*/	\n\t"\
	"/* a[j+p12]: Inputs from [add0,add0+0x200,add1,add1+0x200]+0x180. Outputs into **r6** +0/1, 16/17, 32/33, 48/49: */	\n\t"\
		"addq	$0x40,%%rax	\n\t"\
		"addq	$0x40,%%rbx	\n\t"\
		"addq	$0x40,%%rcx	\n\t"\
		"addq	$0x40,%%rdx	\n\t"\
		"subq	$0x80,%%rsi	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps		 (%%rax),%%ymm1						\n\t		vmovaps	0x020(%%rax),%%ymm3								\n\t"\
		"vmovaps		 (%%rbx),%%ymm5						\n\t		vmovaps	0x020(%%rbx),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps		 (%%rcx),%%ymm6						\n\t		vmovaps	0x020(%%rcx),%%ymm14							\n\t"\
		"vmovaps		 (%%rdx),%%ymm5						\n\t		vmovaps	0x020(%%rdx),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm6 ,%%ymm0				\n\t		vshufpd	$15 ,%%ymm13,%%ymm14 ,%%ymm2					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm6 ,%%ymm6				\n\t		vshufpd	$0  ,%%ymm13,%%ymm14 ,%%ymm14					\n\t"\
		"vperm2f128 $32 ,%%ymm0 ,%%ymm7 ,%%ymm5				\n\t		vperm2f128 $32 ,%%ymm2 ,%%ymm15,%%ymm13		/* outB	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm0 ,%%ymm7 ,%%ymm7				\n\t		vperm2f128 $49 ,%%ymm2 ,%%ymm15,%%ymm15		/* outD	*/	\n\t"\
		"vperm2f128 $32 ,%%ymm6 ,%%ymm1 ,%%ymm0				\n\t		vperm2f128 $32 ,%%ymm14,%%ymm3 ,%%ymm2 		/* outA	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm6 ,%%ymm1 ,%%ymm1				\n\t		vperm2f128 $49 ,%%ymm14,%%ymm3 ,%%ymm3 		/* outC	*/	\n\t"\
		"vmovaps 	%%ymm5 ,0x400(%%rsi)					\n\t		vmovaps %%ymm13,0x420(%%rsi)				/* outB	*/	\n\t"\
		"vmovaps 	%%ymm7 ,0x600(%%rsi)					\n\t		vmovaps %%ymm15,0x620(%%rsi)				/* outD	*/	\n\t"\
		"vmovaps 	%%ymm0 ,     (%%rsi)					\n\t		vmovaps %%ymm2 ,0x020(%%rsi)				/* outA	*/	\n\t"\
		"vmovaps 	%%ymm1 ,0x200(%%rsi)					\n\t		vmovaps %%ymm3 ,0x220(%%rsi)				/* outC	*/	\n\t"\
	"/* a[j+p14]: Inputs from [add0,add0+0x200,add1,add1+0x200]+0x1c0. Outputs into r14 +0/1, 16/17, 32/33, 48/49: */	\n\t"\
		"addq	$0x40,%%rax	\n\t"\
		"addq	$0x40,%%rbx	\n\t"\
		"addq	$0x40,%%rcx	\n\t"\
		"addq	$0x40,%%rdx	\n\t"\
		"addq	$0x100,%%rsi	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps		 (%%rax),%%ymm1						\n\t		vmovaps	0x020(%%rax),%%ymm3								\n\t"\
		"vmovaps		 (%%rbx),%%ymm5						\n\t		vmovaps	0x020(%%rbx),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps		 (%%rcx),%%ymm6						\n\t		vmovaps	0x020(%%rcx),%%ymm14							\n\t"\
		"vmovaps		 (%%rdx),%%ymm5						\n\t		vmovaps	0x020(%%rdx),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm6 ,%%ymm0				\n\t		vshufpd	$15 ,%%ymm13,%%ymm14 ,%%ymm2					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm6 ,%%ymm6				\n\t		vshufpd	$0  ,%%ymm13,%%ymm14 ,%%ymm14					\n\t"\
		"vperm2f128 $32 ,%%ymm0 ,%%ymm7 ,%%ymm5				\n\t		vperm2f128 $32 ,%%ymm2 ,%%ymm15,%%ymm13		/* outB	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm0 ,%%ymm7 ,%%ymm7				\n\t		vperm2f128 $49 ,%%ymm2 ,%%ymm15,%%ymm15		/* outD	*/	\n\t"\
		"vperm2f128 $32 ,%%ymm6 ,%%ymm1 ,%%ymm0				\n\t		vperm2f128 $32 ,%%ymm14,%%ymm3 ,%%ymm2 		/* outA	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm6 ,%%ymm1 ,%%ymm1				\n\t		vperm2f128 $49 ,%%ymm14,%%ymm3 ,%%ymm3 		/* outC	*/	\n\t"\
		"vmovaps 	%%ymm5 ,0x400(%%rsi)					\n\t		vmovaps %%ymm13,0x420(%%rsi)				/* outB	*/	\n\t"\
		"vmovaps 	%%ymm7 ,0x600(%%rsi)					\n\t		vmovaps %%ymm15,0x620(%%rsi)				/* outD	*/	\n\t"\
		"vmovaps 	%%ymm0 ,     (%%rsi)					\n\t		vmovaps %%ymm2 ,0x020(%%rsi)				/* outA	*/	\n\t"\
		"vmovaps 	%%ymm1 ,0x200(%%rsi)					\n\t		vmovaps %%ymm3 ,0x220(%%rsi)				/* outC	*/	\n\t"\
		"\n\t"\
		/************************************************************************/\
		/* Forward DIF radix-32 pass on the interleaved block1 and block2 data: */\
		/************************************************************************/\
	/*...Block 0: */\
		"movq	%[__isrt2],%%rsi							\n\t		leaq	0x8e0(%%rsi),%%rdi	/* two */	\n\t"\
	/*****	SSE2_RADIX4_DIF_4WRAPPER(c00,c08,c10,c18,r00) ****	SSE2_RADIX4_DIF_4WRAPPER_2NDOFTWO(c04,c0C,c14,c1C,r08)	*****/\
		"movq		%[__r00],%%rcx						\n\t		/*addq		$0x100,%%rcx // __r08 */	\n\t"\
		"movq		%[__c00],%%rdx						\n\t		/*addq		$0x100,%%rdx // __c04 */	\n\t"\
		"vmovaps			 (%%rcx),%%ymm0					\n\t		vmovaps		0x100(%%rcx),%%ymm8			\n\t"\
		"vmovaps		0x020(%%rcx),%%ymm1					\n\t		vmovaps		0x120(%%rcx),%%ymm9			\n\t"\
		"vmovaps		%%ymm0,%%ymm2					\n\t		vmovaps		%%ymm8,%%ymm10				\n\t"\
		"vmovaps		%%ymm1,%%ymm3					\n\t		vmovaps		%%ymm9,%%ymm11				\n\t"\
		"vmulpd			 (%%rdx),%%ymm0,%%ymm0				\n\t		vmulpd		0x100(%%rdx),%%ymm8,%%ymm8			\n\t"\
		"vmulpd			 (%%rdx),%%ymm1,%%ymm1				\n\t		vmulpd		0x100(%%rdx),%%ymm9,%%ymm9			\n\t"\
	"vfnmadd231pd	0x020(%%rdx),%%ymm7,%%ymm0				\n\t	vfnmadd231pd	0x120(%%rdx),%%ymm11,%%ymm8 	\n\t"\
	" vfmadd231pd	0x020(%%rdx),%%ymm6,%%ymm1				\n\t	 vfmadd231pd	0x120(%%rdx),%%ymm10,%%ymm9 	\n\t"\
		"vmovaps		%%ymm0,%%ymm2					\n\t		vmovaps		%%ymm8,%%ymm10						\n\t"\
		"vmovaps		%%ymm1,%%ymm3					\n\t		vmovaps		%%ymm9,%%ymm11						\n\t"\
		"addq		$0x040,%%rdx						\n\t"\
		"vmovaps		0x040(%%rcx),%%ymm4					\n\t		vmovaps		0x140(%%rcx),%%ymm12				\n\t"\
		"vmovaps		0x060(%%rcx),%%ymm5					\n\t		vmovaps		0x160(%%rcx),%%ymm13				\n\t"\
		"vmovaps		%%ymm4,%%ymm6					\n\t		vmovaps		%%ymm12,%%ymm14						\n\t"\
		"vmovaps		%%ymm5,%%ymm7					\n\t		vmovaps		%%ymm13,%%ymm15						\n\t"\
		"vmulpd			 (%%rdx),%%ymm4,%%ymm4				\n\t		vmulpd		0x100(%%rdx),%%ymm12,%%ymm12		\n\t"\
		"vmulpd			 (%%rdx),%%ymm5,%%ymm5				\n\t		vmulpd		0x100(%%rdx),%%ymm13,%%ymm13		\n\t"\
	"vfnmadd231pd	0x020(%%rdx),%%ymm7,%%ymm4				\n\t	vfnmadd231pd	0x120(%%rdx),%%ymm15,%%ymm12	\n\t"\
	" vfmadd231pd	0x020(%%rdx),%%ymm6,%%ymm5				\n\t	 vfmadd231pd	0x120(%%rdx),%%ymm14,%%ymm13	\n\t"\
		"vaddpd		%%ymm4,%%ymm0,%%ymm0				\n\t		vaddpd		%%ymm12,%%ymm8,%%ymm8				\n\t"\
		"vaddpd		%%ymm5,%%ymm1,%%ymm1				\n\t		vaddpd		%%ymm13,%%ymm9,%%ymm9				\n\t"\
		"vsubpd		%%ymm4,%%ymm2,%%ymm2				\n\t		vsubpd		%%ymm12,%%ymm10,%%ymm10				\n\t"\
		"vsubpd		%%ymm5,%%ymm3,%%ymm3				\n\t		vsubpd		%%ymm13,%%ymm11,%%ymm11				\n\t"\
		"addq		$0x080,%%rdx						\n\t"\
		"vmovaps		0x0c0(%%rcx),%%ymm4					\n\t		vmovaps		0x1c0(%%rcx),%%ymm12				\n\t"\
		"vmovaps		0x0e0(%%rcx),%%ymm5					\n\t		vmovaps		0x1e0(%%rcx),%%ymm13				\n\t"\
		"vmovaps		%%ymm4,%%ymm6					\n\t		vmovaps		%%ymm12,%%ymm14						\n\t"\
		"vmovaps		%%ymm5,%%ymm7					\n\t		vmovaps		%%ymm13,%%ymm15						\n\t"\
		"vmulpd			 (%%rdx),%%ymm4,%%ymm4				\n\t		vmulpd		0x100(%%rdx),%%ymm12,%%ymm12		\n\t"\
		"vmulpd			 (%%rdx),%%ymm5,%%ymm5				\n\t		vmulpd		0x100(%%rdx),%%ymm13,%%ymm13		\n\t"\
	"vfnmadd231pd	0x020(%%rdx),%%ymm7,%%ymm4				\n\t	vfnmadd231pd	0x120(%%rdx),%%ymm15,%%ymm12	\n\t"\
	" vfmadd231pd	0x020(%%rdx),%%ymm6,%%ymm5				\n\t	 vfmadd231pd	0x120(%%rdx),%%ymm14,%%ymm13	\n\t"\
		"vmovaps		%%ymm5,0x020(%%rcx)			\n\t		vmovaps		%%ymm13,0x120(%%rcx)				\n\t"\
		"vmovaps		%%ymm4,     (%%rcx)			\n\t		vmovaps		%%ymm12,0x100(%%rcx)				\n\t"\
		"subq		$0x040,%%rdx						\n\t"\
		"vmovaps		0x080(%%rcx),%%ymm4					\n\t		vmovaps		0x180(%%rcx),%%ymm12				\n\t"\
		"vmovaps		0x0a0(%%rcx),%%ymm5					\n\t		vmovaps		0x1a0(%%rcx),%%ymm13				\n\t"\
		"vmovaps		%%ymm4,%%ymm6					\n\t		vmovaps		%%ymm12,%%ymm14						\n\t"\
		"vmovaps		%%ymm5,%%ymm7					\n\t		vmovaps		%%ymm13,%%ymm15						\n\t"\
		"vmulpd			 (%%rdx),%%ymm4,%%ymm4				\n\t		vmulpd		0x100(%%rdx),%%ymm12,%%ymm12		\n\t"\
		"vmulpd			 (%%rdx),%%ymm5,%%ymm5				\n\t		vmulpd		0x100(%%rdx),%%ymm13,%%ymm13		\n\t"\
	"vfnmadd231pd	0x020(%%rdx),%%ymm7,%%ymm4				\n\t	vfnmadd231pd	0x120(%%rdx),%%ymm15,%%ymm12	\n\t"\
	" vfmadd231pd	0x020(%%rdx),%%ymm6,%%ymm5				\n\t	 vfmadd231pd	0x120(%%rdx),%%ymm14,%%ymm13	\n\t"\
		"vmovaps		%%ymm5,%%ymm7					\n\t		vmovaps		%%ymm13,%%ymm15					\n\t"\
		"vmovaps		%%ymm4,%%ymm6					\n\t		vmovaps		%%ymm12,%%ymm14					\n\t"\
		"vsubpd			 (%%rcx),%%ymm4,%%ymm4				\n\t		vsubpd		0x100(%%rcx),%%ymm12,%%ymm12	\n\t"\
		"vsubpd		0x020(%%rcx),%%ymm5,%%ymm5				\n\t		vsubpd		0x120(%%rcx),%%ymm13,%%ymm13	\n\t"\
		"vaddpd			 (%%rcx),%%ymm6,%%ymm6				\n\t		vaddpd		0x100(%%rcx),%%ymm14,%%ymm14	\n\t"\
		"vaddpd		0x020(%%rcx),%%ymm7,%%ymm7				\n\t		vaddpd		0x120(%%rcx),%%ymm15,%%ymm15	\n\t"\
		"vsubpd		%%ymm6,%%ymm0,%%ymm0				\n\t		vsubpd		%%ymm14,%%ymm8,%%ymm8		\n\t"\
		"vsubpd		%%ymm7,%%ymm1,%%ymm1				\n\t		vsubpd		%%ymm15,%%ymm9,%%ymm9		\n\t"\
		"vsubpd		%%ymm5,%%ymm2,%%ymm2				\n\t		vsubpd		%%ymm13,%%ymm10,%%ymm10	\n\t"\
		"vsubpd		%%ymm4,%%ymm3,%%ymm3				\n\t		vsubpd		%%ymm12,%%ymm11,%%ymm11	\n\t"\
	"vmovaps	%%ymm12,(%%rcx) 	\n\t"/* spill ymm12 to make room for two */"	vmovaps	(%%rdi),%%ymm12 	\n\t"/* two */\
	"vfmadd132pd	%%ymm12,%%ymm0,%%ymm6					\n\t	vfmadd132pd	%%ymm12,%%ymm8 ,%%ymm14				\n\t"\
	"vfmadd132pd	%%ymm12,%%ymm1,%%ymm7					\n\t	vfmadd132pd	%%ymm12,%%ymm9 ,%%ymm15				\n\t"\
	"vfmadd132pd	%%ymm12,%%ymm2,%%ymm5					\n\t	vfmadd132pd	%%ymm12,%%ymm10,%%ymm13				\n\t"\
	"vfmadd132pd	%%ymm12,%%ymm3,%%ymm4					\n\t	vfmadd132pd	(%%rcx),%%ymm11,%%ymm12				\n\t"\
		"																vmovaps		%%ymm14,0x100(%%rcx)			\n\t"\
		"																vmovaps		%%ymm15,0x120(%%rcx)			\n\t"\
		"																vmovaps		%%ymm10,%%ymm14					\n\t"\
		"																vmovaps		%%ymm13,%%ymm15					\n\t"\
		"																vsubpd		%%ymm12,%%ymm10,%%ymm10	\n\t"\
		"																vsubpd		%%ymm11,%%ymm13,%%ymm13	\n\t"\
		"																vaddpd		%%ymm12,%%ymm14,%%ymm14	\n\t"\
		"																vaddpd		%%ymm11,%%ymm15,%%ymm15	\n\t"\
		/*****	SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r00) *****/\
		"																vmovaps		0x100(%%rcx),%%ymm11		\n\t"\
		"																vmovaps		0x120(%%rcx),%%ymm12		\n\t"\
	"vmovaps	%%ymm13,(%%rcx) 	\n\t"/* spill ymm13 to make room for isrt2 */"	vmovaps	(%%rsi),%%ymm13 	\n\t"/*isrt2*/\
		"vsubpd		%%ymm11,%%ymm6,%%ymm6				\n\t	vfnmadd231pd	%%ymm13,%%ymm10,%%ymm2			\n\t"\
		"vsubpd		%%ymm9,%%ymm0,%%ymm0				\n\t	vfnmadd231pd	%%ymm13,%%ymm15,%%ymm5			\n\t"\
		"vsubpd		%%ymm12,%%ymm7,%%ymm7				\n\t	vfnmadd231pd	%%ymm13,%%ymm14,%%ymm4			\n\t"\
		"vsubpd		%%ymm8 ,%%ymm1,%%ymm1				\n\t	vfnmadd231pd	(%%rcx),%%ymm13,%%ymm3			\n\t"\
		"vmovaps		%%ymm6,0x100(%%rcx)			\n\t		vmovaps		%%ymm2,0x140(%%rcx)	\n\t"\
		"vmovaps		%%ymm0,0x080(%%rcx)			\n\t		vmovaps		%%ymm5,0x0c0(%%rcx)	\n\t"\
		"vmovaps		%%ymm7,0x120(%%rcx)			\n\t		vmovaps		%%ymm4,0x160(%%rcx)	\n\t"\
		"vmovaps		%%ymm1,0x1a0(%%rcx)			\n\t		vmovaps		%%ymm3,0x1e0(%%rcx)	\n\t"\
	"vmovaps	0x40(%%rdi),%%ymm13 \n\t"/* sqrt2 *//* rcol: *2 => *sqrt2 due to above FMA leaving ymm10,13-15 unmultiplied-by-isrt2: */\
	"vfmadd132pd	(%%rdi),%%ymm6,%%ymm11				\n\t	 vfmadd132pd	%%ymm13,%%ymm2,%%ymm10	\n\t"\
	"vfmadd132pd	(%%rdi),%%ymm0,%%ymm9					\n\t	 vfmadd132pd	%%ymm13,%%ymm5,%%ymm15	\n\t"\
	"vfmadd132pd	(%%rdi),%%ymm7,%%ymm12				\n\t	 vfmadd132pd	%%ymm13,%%ymm4,%%ymm14	\n\t"\
	"vfmadd132pd	(%%rdi),%%ymm1,%%ymm8					\n\t	 vfmadd132pd	(%%rcx),%%ymm3,%%ymm13	\n\t"\
		"vmovaps		%%ymm11,     (%%rcx)			\n\t		vmovaps		%%ymm10,0x040(%%rcx)	\n\t"\
		"vmovaps		%%ymm9,0x180(%%rcx)			\n\t		vmovaps		%%ymm15,0x1c0(%%rcx)	\n\t"\
		"vmovaps		%%ymm12,0x020(%%rcx)			\n\t		vmovaps		%%ymm14,0x060(%%rcx)	\n\t"\
		"vmovaps		%%ymm8,0x0a0(%%rcx)			\n\t		vmovaps		%%ymm13,0x0e0(%%rcx)	\n\t"\
		"\n\t"\
	/*...Block 2: */\
	/*****	SSE2_RADIX4_DIF_4WRAPPER(c02,c0A,c12,c1A,r10) ****	SSE2_RADIX4_DIF_4WRAPPER_2NDOFTWO(c06,c0E,c16,c1E,r18)	*****/\
		"movq		%[__r10],%%rcx						\n\t		/*addq		$0x100,%%rcx // __r18 */	\n\t"\
		"movq		%[__c02],%%rdx						\n\t		/*addq		$0x100,%%rdx // __c06 */	\n\t"\
		"vmovaps			 (%%rcx),%%ymm0					\n\t		vmovaps		0x100(%%rcx),%%ymm8			\n\t"\
		"vmovaps		0x020(%%rcx),%%ymm1					\n\t		vmovaps		0x120(%%rcx),%%ymm9			\n\t"\
		"vmovaps		%%ymm0,%%ymm2					\n\t		vmovaps		%%ymm8,%%ymm10				\n\t"\
		"vmovaps		%%ymm1,%%ymm3					\n\t		vmovaps		%%ymm9,%%ymm11				\n\t"\
		"vmulpd			 (%%rdx),%%ymm0,%%ymm0				\n\t		vmulpd		0x100(%%rdx),%%ymm8,%%ymm8			\n\t"\
		"vmulpd			 (%%rdx),%%ymm1,%%ymm1				\n\t		vmulpd		0x100(%%rdx),%%ymm9,%%ymm9			\n\t"\
	"vfnmadd231pd	0x020(%%rdx),%%ymm3,%%ymm0				\n\t	vfnmadd231pd	0x120(%%rdx),%%ymm11,%%ymm8 	\n\t"\
	" vfmadd231pd	0x020(%%rdx),%%ymm2,%%ymm1				\n\t	 vfmadd231pd	0x120(%%rdx),%%ymm10,%%ymm9 	\n\t"\
		"vmovaps		%%ymm0,%%ymm2					\n\t		vmovaps		%%ymm8,%%ymm10						\n\t"\
		"vmovaps		%%ymm1,%%ymm3					\n\t		vmovaps		%%ymm9,%%ymm11						\n\t"\
		"addq		$0x040,%%rdx						\n\t"\
		"vmovaps		0x040(%%rcx),%%ymm4					\n\t		vmovaps		0x140(%%rcx),%%ymm12				\n\t"\
		"vmovaps		0x060(%%rcx),%%ymm5					\n\t		vmovaps		0x160(%%rcx),%%ymm13				\n\t"\
		"vmovaps		%%ymm4,%%ymm6					\n\t		vmovaps		%%ymm12,%%ymm14						\n\t"\
		"vmovaps		%%ymm5,%%ymm7					\n\t		vmovaps		%%ymm13,%%ymm15						\n\t"\
		"vmulpd			 (%%rdx),%%ymm4,%%ymm4				\n\t		vmulpd		0x100(%%rdx),%%ymm12,%%ymm12		\n\t"\
		"vmulpd			 (%%rdx),%%ymm5,%%ymm5				\n\t		vmulpd		0x100(%%rdx),%%ymm13,%%ymm13		\n\t"\
	"vfnmadd231pd	0x020(%%rdx),%%ymm7,%%ymm4				\n\t	vfnmadd231pd	0x120(%%rdx),%%ymm15,%%ymm12	\n\t"\
	" vfmadd231pd	0x020(%%rdx),%%ymm6,%%ymm5				\n\t	 vfmadd231pd	0x120(%%rdx),%%ymm14,%%ymm13	\n\t"\
		"vaddpd		%%ymm4,%%ymm0,%%ymm0				\n\t		vaddpd		%%ymm12,%%ymm8,%%ymm8				\n\t"\
		"vaddpd		%%ymm5,%%ymm1,%%ymm1				\n\t		vaddpd		%%ymm13,%%ymm9,%%ymm9				\n\t"\
		"vsubpd		%%ymm4,%%ymm2,%%ymm2				\n\t		vsubpd		%%ymm12,%%ymm10,%%ymm10				\n\t"\
		"vsubpd		%%ymm5,%%ymm3,%%ymm3				\n\t		vsubpd		%%ymm13,%%ymm11,%%ymm11				\n\t"\
		"addq		$0x080,%%rdx						\n\t"\
		"vmovaps		0x0c0(%%rcx),%%ymm4					\n\t		vmovaps		0x1c0(%%rcx),%%ymm12				\n\t"\
		"vmovaps		0x0e0(%%rcx),%%ymm5					\n\t		vmovaps		0x1e0(%%rcx),%%ymm13				\n\t"\
		"vmovaps		%%ymm4,%%ymm6					\n\t		vmovaps		%%ymm12,%%ymm14						\n\t"\
		"vmovaps		%%ymm5,%%ymm7					\n\t		vmovaps		%%ymm13,%%ymm15						\n\t"\
		"vmulpd			 (%%rdx),%%ymm4,%%ymm4				\n\t		vmulpd		0x100(%%rdx),%%ymm12,%%ymm12		\n\t"\
		"vmulpd			 (%%rdx),%%ymm5,%%ymm5				\n\t		vmulpd		0x100(%%rdx),%%ymm13,%%ymm13		\n\t"\
	"vfnmadd231pd	0x020(%%rdx),%%ymm7,%%ymm4				\n\t	vfnmadd231pd	0x120(%%rdx),%%ymm15,%%ymm12	\n\t"\
	" vfmadd231pd	0x020(%%rdx),%%ymm6,%%ymm5				\n\t	 vfmadd231pd	0x120(%%rdx),%%ymm14,%%ymm13	\n\t"\
		"vmovaps		%%ymm5,0x020(%%rcx)			\n\t		vmovaps		%%ymm13,0x120(%%rcx)				\n\t"\
		"vmovaps		%%ymm4,     (%%rcx)			\n\t		vmovaps		%%ymm12,0x100(%%rcx)				\n\t"\
		"subq		$0x040,%%rdx						\n\t"\
		"vmovaps		0x080(%%rcx),%%ymm4					\n\t		vmovaps		0x180(%%rcx),%%ymm12				\n\t"\
		"vmovaps		0x0a0(%%rcx),%%ymm5					\n\t		vmovaps		0x1a0(%%rcx),%%ymm13				\n\t"\
		"vmovaps		%%ymm4,%%ymm6					\n\t		vmovaps		%%ymm12,%%ymm14						\n\t"\
		"vmovaps		%%ymm5,%%ymm7					\n\t		vmovaps		%%ymm13,%%ymm15						\n\t"\
		"vmulpd			 (%%rdx),%%ymm4,%%ymm4				\n\t		vmulpd		0x100(%%rdx),%%ymm12,%%ymm12		\n\t"\
		"vmulpd			 (%%rdx),%%ymm5,%%ymm5				\n\t		vmulpd		0x100(%%rdx),%%ymm13,%%ymm13		\n\t"\
	"vfnmadd231pd	0x020(%%rdx),%%ymm7,%%ymm4				\n\t	vfnmadd231pd	0x120(%%rdx),%%ymm15,%%ymm12	\n\t"\
	" vfmadd231pd	0x020(%%rdx),%%ymm6,%%ymm5				\n\t	 vfmadd231pd	0x120(%%rdx),%%ymm14,%%ymm13	\n\t"\
		"vmovaps		%%ymm5,%%ymm7					\n\t		vmovaps		%%ymm13,%%ymm15						\n\t"\
		"vmovaps		%%ymm4,%%ymm6					\n\t		vmovaps		%%ymm12,%%ymm14						\n\t"\
		"vsubpd			 (%%rcx),%%ymm4,%%ymm4				\n\t		vsubpd		0x100(%%rcx),%%ymm12,%%ymm12		\n\t"\
		"vsubpd		0x020(%%rcx),%%ymm5,%%ymm5				\n\t		vsubpd		0x120(%%rcx),%%ymm13,%%ymm13		\n\t"\
		"vaddpd			 (%%rcx),%%ymm6,%%ymm6				\n\t		vaddpd		0x100(%%rcx),%%ymm14,%%ymm14		\n\t"\
		"vaddpd		0x020(%%rcx),%%ymm7,%%ymm7				\n\t		vaddpd		0x120(%%rcx),%%ymm15,%%ymm15		\n\t"\
		"vsubpd		%%ymm6,%%ymm0,%%ymm0				\n\t		vsubpd		%%ymm14,%%ymm8,%%ymm8		\n\t"\
		"vsubpd		%%ymm7,%%ymm1,%%ymm1				\n\t		vsubpd		%%ymm15,%%ymm9,%%ymm9		\n\t"\
		"vsubpd		%%ymm5,%%ymm2,%%ymm2				\n\t		vsubpd		%%ymm13,%%ymm10,%%ymm10	\n\t"\
		"vsubpd		%%ymm4,%%ymm3,%%ymm3				\n\t		vsubpd		%%ymm12,%%ymm11,%%ymm11	\n\t"\
	"vmovaps	%%ymm12,(%%rcx) 	\n\t"/* spill ymm12 to make room for two */"	vmovaps	(%%rdi),%%ymm12 	\n\t"/* two */\
	"vfmadd132pd	%%ymm12,%%ymm0,%%ymm6					\n\t	vfmadd132pd	%%ymm12,%%ymm8 ,%%ymm14				\n\t"\
	"vfmadd132pd	%%ymm12,%%ymm1,%%ymm7					\n\t	vfmadd132pd	%%ymm12,%%ymm9 ,%%ymm15				\n\t"\
	"vfmadd132pd	%%ymm12,%%ymm2,%%ymm5					\n\t	vfmadd132pd	%%ymm12,%%ymm10,%%ymm13				\n\t"\
	"vfmadd132pd	%%ymm12,%%ymm3,%%ymm4					\n\t	vfmadd132pd	(%%rcx),%%ymm11,%%ymm12				\n\t"\
		"																vmovaps		%%ymm14,0x100(%%rcx)			\n\t"\
		"																vmovaps		%%ymm15,0x120(%%rcx)			\n\t"\
		"																vmovaps		%%ymm10,%%ymm14					\n\t"\
		"																vmovaps		%%ymm13,%%ymm15					\n\t"\
		"																vsubpd		%%ymm12,%%ymm10,%%ymm10	\n\t"\
		"																vsubpd		%%ymm11,%%ymm13,%%ymm13	\n\t"\
		"																vaddpd		%%ymm12,%%ymm14,%%ymm14	\n\t"\
		"																vaddpd		%%ymm11,%%ymm15,%%ymm15	\n\t"\
		/*****	SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r10) *****/\
		"																vmovaps		0x100(%%rcx),%%ymm11		\n\t"\
		"																vmovaps		0x120(%%rcx),%%ymm12		\n\t"\
	"vmovaps	%%ymm13,(%%rcx) 	\n\t"/* spill ymm13 to make room for isrt2 */"	vmovaps	(%%rsi),%%ymm13 	\n\t"/*isrt2*/\
		"vsubpd		%%ymm11,%%ymm6,%%ymm6				\n\t	vfnmadd231pd	%%ymm13,%%ymm10,%%ymm2			\n\t"\
		"vsubpd		%%ymm9,%%ymm0,%%ymm0				\n\t	vfnmadd231pd	%%ymm13,%%ymm15,%%ymm5			\n\t"\
		"vsubpd		%%ymm12,%%ymm7,%%ymm7				\n\t	vfnmadd231pd	%%ymm13,%%ymm14,%%ymm4			\n\t"\
		"vsubpd		%%ymm8 ,%%ymm1,%%ymm1				\n\t	vfnmadd231pd	(%%rcx),%%ymm13,%%ymm3			\n\t"\
		"vmovaps		%%ymm6,0x100(%%rcx)			\n\t		vmovaps		%%ymm2,0x140(%%rcx)	\n\t"\
		"vmovaps		%%ymm0,0x080(%%rcx)			\n\t		vmovaps		%%ymm5,0x0c0(%%rcx)	\n\t"\
		"vmovaps		%%ymm7,0x120(%%rcx)			\n\t		vmovaps		%%ymm4,0x160(%%rcx)	\n\t"\
		"vmovaps		%%ymm1,0x1a0(%%rcx)			\n\t		vmovaps		%%ymm3,0x1e0(%%rcx)	\n\t"\
	"vmovaps	0x40(%%rdi),%%ymm13 \n\t"/* sqrt2 *//* rcol: *2 => *sqrt2 due to above FMA leaving ymm10,13-15 unmultiplied-by-isrt2: */\
	"vfmadd132pd	(%%rdi),%%ymm6,%%ymm11				\n\t	 vfmadd132pd	%%ymm13,%%ymm2,%%ymm10	\n\t"\
	"vfmadd132pd	(%%rdi),%%ymm0,%%ymm9					\n\t	 vfmadd132pd	%%ymm13,%%ymm5,%%ymm15	\n\t"\
	"vfmadd132pd	(%%rdi),%%ymm7,%%ymm12				\n\t	 vfmadd132pd	%%ymm13,%%ymm4,%%ymm14	\n\t"\
	"vfmadd132pd	(%%rdi),%%ymm1,%%ymm8					\n\t	 vfmadd132pd	(%%rcx),%%ymm3,%%ymm13	\n\t"\
		"vmovaps		%%ymm11,     (%%rcx)			\n\t		vmovaps		%%ymm10,0x040(%%rcx)	\n\t"\
		"vmovaps		%%ymm9,0x180(%%rcx)			\n\t		vmovaps		%%ymm15,0x1c0(%%rcx)	\n\t"\
		"vmovaps		%%ymm12,0x020(%%rcx)			\n\t		vmovaps		%%ymm14,0x060(%%rcx)	\n\t"\
		"vmovaps		%%ymm8,0x0a0(%%rcx)			\n\t		vmovaps		%%ymm13,0x0e0(%%rcx)	\n\t"\
		"\n\t"\
	/************************************************************************************************************/\
	/* Next 2 blocks operate on odd-indexed elements from the unpck*pd commands which we stored to temporaries: */\
	/************************************************************************************************************/\
	/*...Block 3: */\
	/*	SSE2_RADIX4_DIF_4TWIDDLE         (r20,r24,r22,r26,r20,c01) */\
		"addq		$0x200,%%rcx		\n\t"	/***	SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO(r28,r2C,r2A,r2E,r28,c05)	*****/\
		"movq		%[__c01],%%rbx						\n\t	/*	movq		%[__c05],%%rbx	*/		\n\t"\
		"movq		%%rcx,%%rax						\n\t	/*	addq		$0x080,%%rax	*/		\n\t"\
		"addq		$0x040,%%rcx						\n\t	/*	addq		$0x080,%%rcx	*/		\n\t"\
		"vmovaps		 (%%rax),%%ymm0	\n\t movq %%rax,%%rdx \n\t	vmovaps		0x100(%%rax),%%ymm8			\n\t"\
		"vmovaps		 (%%rcx),%%ymm4						\n\t		vmovaps		0x100(%%rcx),%%ymm12		\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1						\n\t		vmovaps		0x120(%%rax),%%ymm9			\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5						\n\t		vmovaps		0x120(%%rcx),%%ymm13		\n\t"\
		"vmovaps		 (%%rbx),%%ymm6						\n\t		vmovaps		0x100(%%rbx),%%ymm14		\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm7						\n\t		vmovaps		0x120(%%rbx),%%ymm15		\n\t"\
		"vmovaps	%%ymm0,%%ymm2						\n\t		vmovaps		%%ymm8,%%ymm10		\n\t"\
		"vmovaps	%%ymm1,%%ymm3						\n\t		vmovaps		%%ymm9,%%ymm11		\n\t"\
		"vmulpd		%%ymm6,%%ymm0,%%ymm0				\n\t		vmulpd		%%ymm14,%%ymm8,%%ymm8			\n\t"\
		"vmulpd		%%ymm6,%%ymm1,%%ymm1				\n\t		vmulpd		%%ymm14,%%ymm9,%%ymm9			\n\t"\
	"vfnmadd231pd	%%ymm7,%%ymm3,%%ymm0			\n\t	vfnmadd231pd	%%ymm15,%%ymm11,%%ymm8 	\n\t"\
	" vfmadd231pd	%%ymm7,%%ymm2,%%ymm1			\n\t	 vfmadd231pd	%%ymm15,%%ymm10,%%ymm9 	\n\t"\
		"vmovaps		%%ymm4,%%ymm6					\n\t		vmovaps		%%ymm12,%%ymm14				\n\t"\
		"vmovaps		%%ymm5,%%ymm7					\n\t		vmovaps		%%ymm13,%%ymm15				\n\t"\
		"vmulpd		0x040(%%rbx),%%ymm4,%%ymm4				\n\t		vmulpd		0x140(%%rbx),%%ymm12,%%ymm12		\n\t"\
		"vmulpd		0x040(%%rbx),%%ymm5,%%ymm5				\n\t		vmulpd		0x140(%%rbx),%%ymm13,%%ymm13		\n\t"\
		"vmovaps		%%ymm0,%%ymm2					\n\t		vmovaps		%%ymm8,%%ymm10				\n\t"\
		"vmovaps		%%ymm1,%%ymm3					\n\t		vmovaps		%%ymm9,%%ymm11				\n\t"\
	"vfnmadd231pd	0x060(%%rbx),%%ymm7,%%ymm4				\n\t	vfnmadd231pd	0x160(%%rbx),%%ymm15,%%ymm12	\n\t"\
	" vfmadd231pd	0x060(%%rbx),%%ymm6,%%ymm5				\n\t	 vfmadd231pd	0x160(%%rbx),%%ymm14,%%ymm13	\n\t"\
		"addq		$0x080,%%rcx						\n\t		addq		$0x0c0,%%rbx				\n\t"\
		"vmovaps			 (%%rcx),%%ymm6					\n\t		vaddpd		%%ymm12,%%ymm8,%%ymm8			\n\t"\
		"vmovaps		0x020(%%rcx),%%ymm7					\n\t		vaddpd		%%ymm13,%%ymm9,%%ymm9			\n\t"\
		"vaddpd		%%ymm4,%%ymm0,%%ymm0				\n\t		vsubpd		%%ymm12,%%ymm10,%%ymm10		\n\t"\
		"vaddpd		%%ymm5,%%ymm1,%%ymm1				\n\t		vsubpd		%%ymm13,%%ymm11,%%ymm11		\n\t"\
		"vsubpd		%%ymm4,%%ymm2,%%ymm2				\n\t		vmovaps		0x100(%%rcx),%%ymm12				\n\t"\
		"vsubpd		%%ymm5,%%ymm3,%%ymm3				\n\t		vmovaps		0x120(%%rcx),%%ymm13				\n\t"\
		"vmovaps		%%ymm6,%%ymm4					\n\t		vmovaps		0x100(%%rcx),%%ymm14				\n\t"\
		"vmovaps		%%ymm7,%%ymm5					\n\t		vmovaps		0x120(%%rcx),%%ymm15				\n\t"\
		"vmulpd			 (%%rbx),%%ymm4,%%ymm4				\n\t		vmulpd		0x100(%%rbx),%%ymm12,%%ymm12		\n\t"\
		"vmulpd			 (%%rbx),%%ymm5,%%ymm5				\n\t		vmulpd		0x100(%%rbx),%%ymm13,%%ymm13		\n\t"\
	"vfnmadd231pd	0x020(%%rbx),%%ymm7,%%ymm4				\n\t	vfnmadd231pd	0x120(%%rbx),%%ymm15,%%ymm12	\n\t"\
	" vfmadd231pd	0x020(%%rbx),%%ymm6,%%ymm5				\n\t	 vfmadd231pd	0x120(%%rbx),%%ymm14,%%ymm13	\n\t"\
		"vmovaps		%%ymm5,0x020(%%rdx)			\n\t		vmovaps		%%ymm13,0x120(%%rdx)				\n\t"\
		"vmovaps		%%ymm4,     (%%rdx)			\n\t		vmovaps		%%ymm12,0x100(%%rdx)				\n\t"\
		"													\n\t		addq	$0x080,%%rax							\n\t"\
		"													\n\t		subq	$0x040,%%rbx							\n\t"\
		"vmovaps			 (%%rax),%%ymm4					\n\t		vmovaps		0x100(%%rax),%%ymm12				\n\t"\
		"vmovaps		0x020(%%rax),%%ymm5					\n\t		vmovaps		0x120(%%rax),%%ymm13				\n\t"\
		"vmovaps		%%ymm4,%%ymm6					\n\t		vmovaps		0x100(%%rax),%%ymm14				\n\t"\
		"vmovaps		%%ymm5,%%ymm7					\n\t		vmovaps		0x120(%%rax),%%ymm15				\n\t"\
		"vmulpd			 (%%rbx),%%ymm4,%%ymm4				\n\t		vmulpd		0x100(%%rbx),%%ymm12,%%ymm12		\n\t"\
		"vmulpd			 (%%rbx),%%ymm5,%%ymm5				\n\t		vmulpd		0x100(%%rbx),%%ymm13,%%ymm13		\n\t"\
	"vfnmadd231pd	0x020(%%rbx),%%ymm7,%%ymm4				\n\t	vfnmadd231pd	0x120(%%rbx),%%ymm15,%%ymm12	\n\t"\
	" vfmadd231pd	0x020(%%rbx),%%ymm6,%%ymm5				\n\t	 vfmadd231pd	0x120(%%rbx),%%ymm14,%%ymm13	\n\t"\
		"vmovaps		%%ymm5,%%ymm7					\n\t		vmovaps		%%ymm13,%%ymm15				\n\t"\
		"vmovaps		%%ymm4,%%ymm6					\n\t		vmovaps		%%ymm12,%%ymm14				\n\t"\
		"vsubpd			 (%%rdx),%%ymm4,%%ymm4				\n\t		vsubpd		0x100(%%rdx),%%ymm12,%%ymm12		\n\t"\
		"vsubpd		0x020(%%rdx),%%ymm5,%%ymm5				\n\t		vsubpd		0x120(%%rdx),%%ymm13,%%ymm13		\n\t"\
		"vaddpd			 (%%rdx),%%ymm6,%%ymm6				\n\t		vaddpd		0x100(%%rdx),%%ymm14,%%ymm14		\n\t"\
		"vaddpd		0x020(%%rdx),%%ymm7,%%ymm7				\n\t		vaddpd		0x120(%%rdx),%%ymm15,%%ymm15		\n\t"\
		"vsubpd		%%ymm6,%%ymm0,%%ymm0				\n\t		vsubpd		%%ymm14,%%ymm8,%%ymm8		\n\t"\
		"vsubpd		%%ymm7,%%ymm1,%%ymm1				\n\t		vsubpd		%%ymm15,%%ymm9,%%ymm9		\n\t"\
		"vsubpd		%%ymm5,%%ymm2,%%ymm2				\n\t		vsubpd		%%ymm13,%%ymm10,%%ymm10	\n\t"\
		"vsubpd		%%ymm4,%%ymm3,%%ymm3				\n\t		vsubpd		%%ymm12,%%ymm11,%%ymm11	\n\t"\
	"vmovaps	%%ymm12,(%%rdx) 	\n\t"/* spill ymm12 to make room for two */"	vmovaps	(%%rdi),%%ymm12 	\n\t"/* two */\
	"vfmadd132pd	%%ymm12,%%ymm0,%%ymm6					\n\t	vfmadd132pd	%%ymm12,%%ymm8 ,%%ymm14				\n\t"\
	"vfmadd132pd	%%ymm12,%%ymm1,%%ymm7					\n\t	vfmadd132pd	%%ymm12,%%ymm9 ,%%ymm15				\n\t"\
	"vfmadd132pd	%%ymm12,%%ymm2,%%ymm5					\n\t	vfmadd132pd	%%ymm12,%%ymm10,%%ymm13				\n\t"\
	"vfmadd132pd	%%ymm12,%%ymm3,%%ymm4					\n\t	vfmadd132pd	(%%rdx),%%ymm11,%%ymm12				\n\t"\
		"																vmovaps		%%ymm14,0x100(%%rdx)			\n\t"\
		"																vmovaps		%%ymm15,0x120(%%rdx)			\n\t"\
		"																vmovaps		%%ymm10,%%ymm14					\n\t"\
		"																vmovaps		%%ymm13,%%ymm15					\n\t"\
		"																vsubpd		%%ymm12,%%ymm10,%%ymm10	\n\t"\
		"																vsubpd		%%ymm11,%%ymm13,%%ymm13	\n\t"\
		"																vaddpd		%%ymm12,%%ymm14,%%ymm14	\n\t"\
		"																vaddpd		%%ymm11,%%ymm15,%%ymm15	\n\t"\
		/*****	SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r20) *****/\
		"																vmovaps		0x100(%%rdx),%%ymm11		\n\t"\
		"																vmovaps		0x120(%%rdx),%%ymm12		\n\t"\
	"vmovaps	%%ymm13,(%%rdx) 	\n\t"/* spill ymm13 to make room for isrt2 */"	vmovaps	(%%rsi),%%ymm13 	\n\t"/*isrt2*/\
		"vsubpd		%%ymm11,%%ymm6,%%ymm6				\n\t	vfnmadd231pd	%%ymm13,%%ymm10,%%ymm2			\n\t"\
		"vsubpd		%%ymm9,%%ymm0,%%ymm0				\n\t	vfnmadd231pd	%%ymm13,%%ymm15,%%ymm5			\n\t"\
		"vsubpd		%%ymm12,%%ymm7,%%ymm7				\n\t	vfnmadd231pd	%%ymm13,%%ymm14,%%ymm4			\n\t"\
		"vsubpd		%%ymm8 ,%%ymm1,%%ymm1				\n\t	vfnmadd231pd	(%%rdx),%%ymm13,%%ymm3			\n\t"\
		"vmovaps		%%ymm6,0x100(%%rdx)			\n\t		vmovaps		%%ymm2,0x140(%%rdx)	\n\t"\
		"vmovaps		%%ymm0,0x080(%%rdx)			\n\t		vmovaps		%%ymm5,0x0c0(%%rdx)	\n\t"\
		"vmovaps		%%ymm7,0x120(%%rdx)			\n\t		vmovaps		%%ymm4,0x160(%%rdx)	\n\t"\
		"vmovaps		%%ymm1,0x1a0(%%rdx)			\n\t		vmovaps		%%ymm3,0x1e0(%%rdx)	\n\t"\
	"vmovaps	0x40(%%rdi),%%ymm13 \n\t"/* sqrt2 *//* rcol: *2 => *sqrt2 due to above FMA leaving ymm10,13-15 unmultiplied-by-isrt2: */\
	"vfmadd132pd	(%%rdi),%%ymm6,%%ymm11				\n\t	 vfmadd132pd	%%ymm13,%%ymm2,%%ymm10	\n\t"\
	"vfmadd132pd	(%%rdi),%%ymm0,%%ymm9					\n\t	 vfmadd132pd	%%ymm13,%%ymm5,%%ymm15	\n\t"\
	"vfmadd132pd	(%%rdi),%%ymm7,%%ymm12				\n\t	 vfmadd132pd	%%ymm13,%%ymm4,%%ymm14	\n\t"\
	"vfmadd132pd	(%%rdi),%%ymm1,%%ymm8					\n\t	 vfmadd132pd	(%%rdx),%%ymm3,%%ymm13	\n\t"\
		"vmovaps		%%ymm11,     (%%rdx)			\n\t		vmovaps		%%ymm10,0x040(%%rdx)	\n\t"\
		"vmovaps		%%ymm9,0x180(%%rdx)			\n\t		vmovaps		%%ymm15,0x1c0(%%rdx)	\n\t"\
		"vmovaps		%%ymm12,0x020(%%rdx)			\n\t		vmovaps		%%ymm14,0x060(%%rdx)	\n\t"\
		"vmovaps		%%ymm8,0x0a0(%%rdx)			\n\t		vmovaps		%%ymm13,0x0e0(%%rdx)	\n\t"\
		"\n\t"\
	/*...Block 4: */\
		"/*****	SSE2_RADIX4_DIF_4TWIDDLE         (r30,r34,r32,r36,r30,c03)	*****/\n\t"\
		"movq		%[__c03],%%rbx					\n\t"\
		"movq		%[__r30],%%rax					\n\t"\
		"movq		%%rax,%%rcx		/*****	SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO(r38,r3C,r3A,r3E,r38,c07)	*****/\n\t"\
		"addq		$0x040,%%rcx					\n\t"\
		"vmovaps		 (%%rax),%%ymm0	\n\t movq %%rax,%%rdx \n\t	vmovaps		0x100(%%rax),%%ymm8			\n\t"\
		"vmovaps		 (%%rcx),%%ymm4						\n\t		vmovaps		0x100(%%rcx),%%ymm12		\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1						\n\t		vmovaps		0x120(%%rax),%%ymm9			\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5						\n\t		vmovaps		0x120(%%rcx),%%ymm13		\n\t"\
		"vmovaps		 (%%rbx),%%ymm6						\n\t		vmovaps		0x100(%%rbx),%%ymm14		\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm7						\n\t		vmovaps		0x120(%%rbx),%%ymm15		\n\t"\
		"vmovaps	%%ymm0,%%ymm2						\n\t		vmovaps		%%ymm8,%%ymm10		\n\t"\
		"vmovaps	%%ymm1,%%ymm3						\n\t		vmovaps		%%ymm9,%%ymm11		\n\t"\
		"vmulpd		%%ymm6,%%ymm0,%%ymm0				\n\t		vmulpd		%%ymm14,%%ymm8,%%ymm8			\n\t"\
		"vmulpd		%%ymm6,%%ymm1,%%ymm1				\n\t		vmulpd		%%ymm14,%%ymm9,%%ymm9			\n\t"\
	"vfnmadd231pd	%%ymm7,%%ymm3,%%ymm0			\n\t	vfnmadd231pd	%%ymm15,%%ymm11,%%ymm8 	\n\t"\
	" vfmadd231pd	%%ymm7,%%ymm2,%%ymm1			\n\t	 vfmadd231pd	%%ymm15,%%ymm10,%%ymm9 	\n\t"\
		"vmovaps		%%ymm4,%%ymm6					\n\t		vmovaps		%%ymm12,%%ymm14				\n\t"\
		"vmovaps		%%ymm5,%%ymm7					\n\t		vmovaps		%%ymm13,%%ymm15				\n\t"\
		"vmulpd		0x040(%%rbx),%%ymm4,%%ymm4				\n\t		vmulpd		0x140(%%rbx),%%ymm12,%%ymm12		\n\t"\
		"vmulpd		0x040(%%rbx),%%ymm5,%%ymm5				\n\t		vmulpd		0x140(%%rbx),%%ymm13,%%ymm13		\n\t"\
		"vmovaps		%%ymm0,%%ymm2					\n\t		vmovaps		%%ymm8,%%ymm10				\n\t"\
		"vmovaps		%%ymm1,%%ymm3					\n\t		vmovaps		%%ymm9,%%ymm11				\n\t"\
	"vfnmadd231pd	0x060(%%rbx),%%ymm7,%%ymm4				\n\t	vfnmadd231pd	0x160(%%rbx),%%ymm15,%%ymm12	\n\t"\
	" vfmadd231pd	0x060(%%rbx),%%ymm6,%%ymm5				\n\t	 vfmadd231pd	0x160(%%rbx),%%ymm14,%%ymm13	\n\t"\
		"addq		$0x080,%%rcx						\n\t		addq		$0x0c0,%%rbx				\n\t"\
		"vmovaps			 (%%rcx),%%ymm6					\n\t		vaddpd		%%ymm12,%%ymm8,%%ymm8			\n\t"\
		"vmovaps		0x020(%%rcx),%%ymm7					\n\t		vaddpd		%%ymm13,%%ymm9,%%ymm9			\n\t"\
		"vaddpd		%%ymm4,%%ymm0,%%ymm0				\n\t		vsubpd		%%ymm12,%%ymm10,%%ymm10		\n\t"\
		"vaddpd		%%ymm5,%%ymm1,%%ymm1				\n\t		vsubpd		%%ymm13,%%ymm11,%%ymm11		\n\t"\
		"vsubpd		%%ymm4,%%ymm2,%%ymm2				\n\t		vmovaps		0x100(%%rcx),%%ymm12				\n\t"\
		"vsubpd		%%ymm5,%%ymm3,%%ymm3				\n\t		vmovaps		0x120(%%rcx),%%ymm13				\n\t"\
		"vmovaps		%%ymm6,%%ymm4					\n\t		vmovaps		0x100(%%rcx),%%ymm14				\n\t"\
		"vmovaps		%%ymm7,%%ymm5					\n\t		vmovaps		0x120(%%rcx),%%ymm15				\n\t"\
		"vmulpd			 (%%rbx),%%ymm4,%%ymm4				\n\t		vmulpd		0x100(%%rbx),%%ymm12,%%ymm12		\n\t"\
		"vmulpd			 (%%rbx),%%ymm5,%%ymm5				\n\t		vmulpd		0x100(%%rbx),%%ymm13,%%ymm13		\n\t"\
	"vfnmadd231pd	0x020(%%rbx),%%ymm7,%%ymm4				\n\t	vfnmadd231pd	0x120(%%rbx),%%ymm15,%%ymm12	\n\t"\
	" vfmadd231pd	0x020(%%rbx),%%ymm6,%%ymm5				\n\t	 vfmadd231pd	0x120(%%rbx),%%ymm14,%%ymm13	\n\t"\
		"vmovaps		%%ymm5,0x020(%%rdx)			\n\t		vmovaps		%%ymm13,0x120(%%rdx)				\n\t"\
		"vmovaps		%%ymm4,     (%%rdx)			\n\t		vmovaps		%%ymm12,0x100(%%rdx)				\n\t"\
		"													\n\t		addq	$0x080,%%rax							\n\t"\
		"													\n\t		subq	$0x040,%%rbx							\n\t"\
		"vmovaps			 (%%rax),%%ymm4					\n\t		vmovaps		0x100(%%rax),%%ymm12				\n\t"\
		"vmovaps		0x020(%%rax),%%ymm5					\n\t		vmovaps		0x120(%%rax),%%ymm13				\n\t"\
		"vmovaps		%%ymm4,%%ymm6					\n\t		vmovaps		0x100(%%rax),%%ymm14				\n\t"\
		"vmovaps		%%ymm5,%%ymm7					\n\t		vmovaps		0x120(%%rax),%%ymm15				\n\t"\
		"vmulpd			 (%%rbx),%%ymm4,%%ymm4				\n\t		vmulpd		0x100(%%rbx),%%ymm12,%%ymm12		\n\t"\
		"vmulpd			 (%%rbx),%%ymm5,%%ymm5				\n\t		vmulpd		0x100(%%rbx),%%ymm13,%%ymm13		\n\t"\
	"vfnmadd231pd	0x020(%%rbx),%%ymm7,%%ymm4				\n\t	vfnmadd231pd	0x120(%%rbx),%%ymm15,%%ymm12	\n\t"\
	" vfmadd231pd	0x020(%%rbx),%%ymm6,%%ymm5				\n\t	 vfmadd231pd	0x120(%%rbx),%%ymm14,%%ymm13	\n\t"\
		"vmovaps		%%ymm5,%%ymm7					\n\t		vmovaps		%%ymm13,%%ymm15				\n\t"\
		"vmovaps		%%ymm4,%%ymm6					\n\t		vmovaps		%%ymm12,%%ymm14				\n\t"\
		"vsubpd			 (%%rdx),%%ymm4,%%ymm4				\n\t		vsubpd		0x100(%%rdx),%%ymm12,%%ymm12		\n\t"\
		"vsubpd		0x020(%%rdx),%%ymm5,%%ymm5				\n\t		vsubpd		0x120(%%rdx),%%ymm13,%%ymm13		\n\t"\
		"vaddpd			 (%%rdx),%%ymm6,%%ymm6				\n\t		vaddpd		0x100(%%rdx),%%ymm14,%%ymm14		\n\t"\
		"vaddpd		0x020(%%rdx),%%ymm7,%%ymm7				\n\t		vaddpd		0x120(%%rdx),%%ymm15,%%ymm15		\n\t"\
		"vsubpd		%%ymm6,%%ymm0,%%ymm0				\n\t		vsubpd		%%ymm14,%%ymm8,%%ymm8		\n\t"\
		"vsubpd		%%ymm7,%%ymm1,%%ymm1				\n\t		vsubpd		%%ymm15,%%ymm9,%%ymm9		\n\t"\
		"vsubpd		%%ymm5,%%ymm2,%%ymm2				\n\t		vsubpd		%%ymm13,%%ymm10,%%ymm10	\n\t"\
		"vsubpd		%%ymm4,%%ymm3,%%ymm3				\n\t		vsubpd		%%ymm12,%%ymm11,%%ymm11	\n\t"\
	"vmovaps	%%ymm12,(%%rdx) 	\n\t"/* spill ymm12 to make room for two */"	vmovaps	(%%rdi),%%ymm12 	\n\t"/* two */\
	"vfmadd132pd	%%ymm12,%%ymm0,%%ymm6					\n\t	vfmadd132pd	%%ymm12,%%ymm8 ,%%ymm14				\n\t"\
	"vfmadd132pd	%%ymm12,%%ymm1,%%ymm7					\n\t	vfmadd132pd	%%ymm12,%%ymm9 ,%%ymm15				\n\t"\
	"vfmadd132pd	%%ymm12,%%ymm2,%%ymm5					\n\t	vfmadd132pd	%%ymm12,%%ymm10,%%ymm13				\n\t"\
	"vfmadd132pd	%%ymm12,%%ymm3,%%ymm4					\n\t	vfmadd132pd	(%%rdx),%%ymm11,%%ymm12				\n\t"\
		"																vmovaps		%%ymm14,0x100(%%rdx)			\n\t"\
		"																vmovaps		%%ymm15,0x120(%%rdx)			\n\t"\
		"																vmovaps		%%ymm10,%%ymm14					\n\t"\
		"																vmovaps		%%ymm13,%%ymm15					\n\t"\
		"																vsubpd		%%ymm12,%%ymm10,%%ymm10	\n\t"\
		"																vsubpd		%%ymm11,%%ymm13,%%ymm13	\n\t"\
		"																vaddpd		%%ymm12,%%ymm14,%%ymm14	\n\t"\
		"																vaddpd		%%ymm11,%%ymm15,%%ymm15	\n\t"\
		/*****	SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r30) *****/\
		"																vmovaps		0x100(%%rdx),%%ymm11		\n\t"\
		"																vmovaps		0x120(%%rdx),%%ymm12		\n\t"\
	"vmovaps	%%ymm13,(%%rdx) 	\n\t"/* spill ymm13 to make room for isrt2 */"	vmovaps	(%%rsi),%%ymm13 	\n\t"/*isrt2*/\
		"vsubpd		%%ymm11,%%ymm6,%%ymm6				\n\t	vfnmadd231pd	%%ymm13,%%ymm10,%%ymm2			\n\t"\
		"vsubpd		%%ymm9,%%ymm0,%%ymm0				\n\t	vfnmadd231pd	%%ymm13,%%ymm15,%%ymm5			\n\t"\
		"vsubpd		%%ymm12,%%ymm7,%%ymm7				\n\t	vfnmadd231pd	%%ymm13,%%ymm14,%%ymm4			\n\t"\
		"vsubpd		%%ymm8 ,%%ymm1,%%ymm1				\n\t	vfnmadd231pd	(%%rdx),%%ymm13,%%ymm3			\n\t"\
		"vmovaps		%%ymm6,0x100(%%rdx)			\n\t		vmovaps		%%ymm2,0x140(%%rdx)	\n\t"\
		"vmovaps		%%ymm0,0x080(%%rdx)			\n\t		vmovaps		%%ymm5,0x0c0(%%rdx)	\n\t"\
		"vmovaps		%%ymm7,0x120(%%rdx)			\n\t		vmovaps		%%ymm4,0x160(%%rdx)	\n\t"\
		"vmovaps		%%ymm1,0x1a0(%%rdx)			\n\t		vmovaps		%%ymm3,0x1e0(%%rdx)	\n\t"\
	"vmovaps	0x40(%%rdi),%%ymm13 \n\t"/* sqrt2 *//* rcol: *2 => *sqrt2 due to above FMA leaving ymm10,13-15 unmultiplied-by-isrt2: */\
	"vfmadd132pd	(%%rdi),%%ymm6,%%ymm11				\n\t	 vfmadd132pd	%%ymm13,%%ymm2,%%ymm10	\n\t"\
	"vfmadd132pd	(%%rdi),%%ymm0,%%ymm9					\n\t	 vfmadd132pd	%%ymm13,%%ymm5,%%ymm15	\n\t"\
	"vfmadd132pd	(%%rdi),%%ymm7,%%ymm12				\n\t	 vfmadd132pd	%%ymm13,%%ymm4,%%ymm14	\n\t"\
	"vfmadd132pd	(%%rdi),%%ymm1,%%ymm8					\n\t	 vfmadd132pd	(%%rdx),%%ymm3,%%ymm13	\n\t"\
		"vmovaps		%%ymm11,     (%%rdx)			\n\t		vmovaps		%%ymm10,0x040(%%rdx)	\n\t"\
		"vmovaps		%%ymm9,0x180(%%rdx)			\n\t		vmovaps		%%ymm15,0x1c0(%%rdx)	\n\t"\
		"vmovaps		%%ymm12,0x020(%%rdx)			\n\t		vmovaps		%%ymm14,0x060(%%rdx)	\n\t"\
		"vmovaps		%%ymm8,0x0a0(%%rdx)			\n\t		vmovaps		%%ymm13,0x0e0(%%rdx)	\n\t"\
		"\n\t"\
	/**********************************************************************************/\
	/*...and now do eight radix-4 transforms, including the internal twiddle factors: */\
	/**********************************************************************************/\
		"movq		%[__isrt2],%%rsi		\n\t"\
		/*...Block 1: t00,t10,t20,t30	*/							/*...Block 5: t08,t18,t28,t38	*/\
		"movq		%[__r00],%%rax					\n\t		leaq	0x100(%%rax),%%r10	\n\t"\
		"movq		%[__r10],%%rbx					\n\t		leaq	0x100(%%rbx),%%r11	\n\t"\
		"movq		%[__r20],%%rcx					\n\t		leaq	0x100(%%rcx),%%r12	\n\t"\
		"movq		%[__r30],%%rdx					\n\t		leaq	0x100(%%rdx),%%r13	\n\t"\
		"vmovaps			  (%%rax),%%ymm0				\n\t		vmovaps			(%%r10),%%ymm8			\n\t"\
		"vmovaps		 0x020(%%rax),%%ymm1				\n\t		vmovaps		0x20(%%r10),%%ymm9			\n\t"\
		"vmovaps			  (%%rbx),%%ymm2				\n\t		vmovaps			(%%r11),%%ymm10			\n\t"\
		"vmovaps		 0x020(%%rbx),%%ymm3				\n\t		vmovaps		0x20(%%r11),%%ymm11			\n\t"\
		"vaddpd			   %%ymm0,%%ymm2,%%ymm2			\n\t		vaddpd			 %%ymm9 ,%%ymm10,%%ymm10	\n\t"\
		"vaddpd			   %%ymm1,%%ymm3,%%ymm3			\n\t		vaddpd			 %%ymm8 ,%%ymm11,%%ymm11	\n\t"\
		"vsubpd			  (%%rbx),%%ymm0,%%ymm0			\n\t		vsubpd		0x020(%%r11),%%ymm8,%%ymm8		\n\t"\
		"vsubpd		 0x020(%%rbx),%%ymm1,%%ymm1			\n\t		vsubpd			 (%%r11),%%ymm9,%%ymm9		\n\t"\
		"vmovaps			  (%%rcx),%%ymm4				\n\t		vmovaps			(%%r12),%%ymm12					\n\t"\
		"vmovaps		 0x020(%%rcx),%%ymm5				\n\t		vmovaps		0x20(%%r12),%%ymm13					\n\t"\
		"vmovaps			  (%%rdx),%%ymm6				\n\t		vmovaps			(%%r13),%%ymm14					\n\t"\
		"vmovaps		 0x020(%%rdx),%%ymm7				\n\t		vmovaps		0x20(%%r13),%%ymm15					\n\t"\
		"vaddpd			   %%ymm4,%%ymm6,%%ymm6			\n\t		vsubpd			 %%ymm13,%%ymm12,%%ymm12	\n\t"\
		"vaddpd			   %%ymm5,%%ymm7,%%ymm7			\n\t		vaddpd			(%%r12),%%ymm13,%%ymm13	\n\t"\
		"vsubpd			  (%%rdx),%%ymm4,%%ymm4			\n\t		vaddpd			 %%ymm15,%%ymm14,%%ymm14	\n\t"\
		"vsubpd		 0x020(%%rdx),%%ymm5,%%ymm5			\n\t		vsubpd			(%%r13),%%ymm15,%%ymm15	\n\t"\
		"vsubpd		%%ymm6,%%ymm2,%%ymm2			\n\t		vmulpd		(%%rsi),%%ymm12,%%ymm12		\n\t"/* isrt2 */\
		"vsubpd		%%ymm7,%%ymm3,%%ymm3			\n\t		vmulpd		(%%rsi),%%ymm13,%%ymm13		\n\t"\
	"vfmadd132pd	(%%rdi),%%ymm2,%%ymm6			\n\t	vfnmadd231pd	(%%rsi),%%ymm14,%%ymm12		\n\t"\
	"vfmadd132pd	(%%rdi),%%ymm3,%%ymm7			\n\t	vfnmadd231pd	(%%rsi),%%ymm15,%%ymm13		\n\t"\
		"vmovaps		%%ymm2,      (%%rcx)			\n\t	 vfmadd132pd	0x40(%%rdi),%%ymm12,%%ymm14		\n\t"/* sqrt2 */\
		"vmovaps		%%ymm3, 0x020(%%rcx)			\n\t	 vfmadd132pd	0x40(%%rdi),%%ymm13,%%ymm15		\n\t"\
		"													\n\t		vsubpd		%%ymm12,%%ymm8,%%ymm8			\n\t"\
		"													\n\t		vsubpd		%%ymm13,%%ymm10,%%ymm10		\n\t"\
		"vmovaps		%%ymm6,      (%%rax)			\n\t	 vfmadd132pd	(%%rdi),%%ymm8,%%ymm12		\n\t"\
		"vmovaps		%%ymm7, 0x020(%%rax)			\n\t	 vfmadd132pd	(%%rdi),%%ymm10,%%ymm13		\n\t"\
		"vsubpd		%%ymm5,%%ymm0,%%ymm0			\n\t		vmovaps		%%ymm8 ,     (%%r12)				\n\t"\
		"vsubpd		%%ymm4,%%ymm1,%%ymm1			\n\t		vmovaps		%%ymm10,0x020(%%r12)				\n\t"\
	"vfmadd132pd	(%%rdi),%%ymm0,%%ymm5			\n\t		vmovaps		%%ymm12,     (%%r10)				\n\t"\
	"vfmadd132pd	(%%rdi),%%ymm1,%%ymm4			\n\t		vmovaps		%%ymm13,0x020(%%r10)				\n\t"\
		"vmovaps		%%ymm0,      (%%rbx)			\n\t		vsubpd		%%ymm15,%%ymm11,%%ymm11		\n\t"\
		"vmovaps		%%ymm1, 0x020(%%rdx)			\n\t		vsubpd		%%ymm14,%%ymm9,%%ymm9			\n\t"\
		"													\n\t	 vfmadd132pd	(%%rdi),%%ymm11,%%ymm15		\n\t"\
		"													\n\t	 vfmadd132pd	(%%rdi),%%ymm9 ,%%ymm14		\n\t"\
		"vmovaps		%%ymm5,      (%%rdx)			\n\t		vmovaps		%%ymm11,     (%%r11)				\n\t"\
		"vmovaps		%%ymm4, 0x020(%%rbx)			\n\t		vmovaps		%%ymm9 ,0x020(%%r13)				\n\t"\
		"																vmovaps		%%ymm15,     (%%r13)				\n\t"\
		"																vmovaps		%%ymm14,0x020(%%r11)				\n\t"\
		/*...Block 3: t04,t14,t24,t34	*/							/*...Block 7: t0C,t1C,t2C,t3C	*/\
		"addq		$0x080,%%rax						\n\t		addq		$0x080,%%r10					\n\t"\
		"addq		$0x080,%%rbx						\n\t		addq		$0x080,%%r11					\n\t"\
		"addq		$0x080,%%rcx						\n\t		addq		$0x080,%%r12					\n\t"\
		"addq		$0x080,%%rdx						\n\t		addq		$0x080,%%r13					\n\t"\
		"vmovaps			 (%%rcx),%%ymm0				\n\t		vmovaps			 (%%r12),%%ymm8 				\n\t"\
		"vmovaps			 (%%rdx),%%ymm2				\n\t		vmovaps			 (%%r13),%%ymm10				\n\t"\
		"vmovaps		0x020(%%rcx),%%ymm1				\n\t		vmovaps		0x020(%%r12),%%ymm9					\n\t"\
		"vmovaps		0x020(%%rdx),%%ymm3				\n\t		vmovaps		0x020(%%r13),%%ymm11				\n\t"\
		"vmovaps			  %%ymm0,%%ymm4				\n\t		vmovaps			 %%ymm8 ,%%ymm12				\n\t"\
		"vmovaps			  %%ymm2,%%ymm6				\n\t		vmovaps			 %%ymm10,%%ymm14				\n\t"\
		"vmovaps			  %%ymm1,%%ymm5				\n\t		vmovaps			 %%ymm9,%%ymm13				\n\t"\
	/*	"vmovaps			  %%ymm3,%%ymm7				\n\t		vmovaps			 %%ymm11,%%ymm15				\n\t"*/\
	/* Strategic register-usage for sincos consts here cuts loads from 16 to 4: */\
		"vmovaps		0x040(%%rsi),%%ymm7 	\n\t"/* s */"		vmovaps		0x020(%%rsi),%%ymm15	\n\t"/* c */\
		"vmulpd		%%ymm7 ,%%ymm1,%%ymm1					\n\t		vmulpd		 %%ymm15,%%ymm9,%%ymm9				\n\t"\
		"vmulpd		%%ymm15,%%ymm3,%%ymm3					\n\t		vmulpd		 %%ymm7 ,%%ymm11,%%ymm11			\n\t"\
		"vmulpd		%%ymm7 ,%%ymm0,%%ymm0					\n\t		vmulpd		 %%ymm15,%%ymm8,%%ymm8				\n\t"\
		"vmulpd		%%ymm15,%%ymm2,%%ymm2					\n\t		vmulpd		 %%ymm7 ,%%ymm10,%%ymm10			\n\t"\
	"vfmsub132pd	%%ymm15,%%ymm1,%%ymm4					\n\t	vfmsub132pd		 %%ymm7 ,%%ymm9 ,%%ymm12			\n\t"\
	"vfmsub132pd	%%ymm7 ,%%ymm3,%%ymm6					\n\t	vfmsub132pd		 %%ymm15,%%ymm11,%%ymm14			\n\t"\
	"vfmadd132pd	%%ymm15,%%ymm0,%%ymm5					\n\t	vfmadd132pd		 %%ymm7 ,%%ymm8 ,%%ymm13			\n\t"\
	"vfmadd132pd 0x20(%%rdx),%%ymm2,%%ymm7					\n\t	vfmadd132pd 0x20(%%r13) ,%%ymm10,%%ymm15			\n\t"\
		"vsubpd		%%ymm6,%%ymm4,%%ymm4				\n\t		vsubpd		%%ymm14,%%ymm12,%%ymm12		\n\t"\
		"vsubpd		%%ymm7,%%ymm5,%%ymm5				\n\t		vsubpd		%%ymm15,%%ymm13,%%ymm13		\n\t"\
		"vmovaps			  (%%rbx),%%ymm2				\n\t		vmovaps			 (%%r11),%%ymm10				\n\t"\
		"vmovaps		 0x020(%%rbx),%%ymm3				\n\t		vmovaps		0x020(%%r11),%%ymm11				\n\t"\
	"vfmadd132pd	(%%rdi),%%ymm4,%%ymm6				\n\t	vfmadd132pd		(%%rdi),%%ymm12,%%ymm14		\n\t"\
	"vfmadd132pd	(%%rdi),%%ymm5,%%ymm7				\n\t	vfmadd132pd		(%%rdi),%%ymm13,%%ymm15		\n\t"\
		"vsubpd		0x020(%%rbx),%%ymm2,%%ymm2				\n\t		vaddpd		0x020(%%r11),%%ymm10,%%ymm10	\n\t"\
		"vaddpd			 (%%rbx),%%ymm3,%%ymm3				\n\t		vsubpd			 (%%r11),%%ymm11,%%ymm11	\n\t"\
	"vmovaps	%%ymm6 ,(%%rcx) 	\n\t"/* spill ymm6  to make room for isrt2 */"	vmovaps	(%%rsi),%%ymm6  	\n\t"/* isrt2 */\
		"vmovaps			  (%%rax),%%ymm0				\n\t		vmovaps			  (%%r10),%%ymm8				\n\t"\
		"vmovaps		 0x020(%%rax),%%ymm1				\n\t		vmovaps		 0x020(%%r10),%%ymm9				\n\t"\
	"vfnmadd231pd	%%ymm6 ,%%ymm2,%%ymm0					\n\t	vfnmadd231pd	%%ymm6 ,%%ymm10,%%ymm8			\n\t"\
	"vfnmadd231pd	%%ymm6 ,%%ymm3,%%ymm1					\n\t	vfnmadd231pd	%%ymm6 ,%%ymm11,%%ymm9			\n\t"\
	" vfmadd213pd		 (%%rax),%%ymm6 ,%%ymm2				\n\t	 vfmadd213pd		 (%%r10),%%ymm6 ,%%ymm10	\n\t"\
	" vfmadd213pd	0x020(%%rax),%%ymm6 ,%%ymm3				\n\t	 vfmadd213pd	0x020(%%r10),%%ymm6 ,%%ymm11	\n\t"\
	"vmovaps	(%%rcx),%%ymm6 	\n\t"/* restore spill */\
		"vsubpd		%%ymm6,%%ymm2,%%ymm2				\n\t		vsubpd		%%ymm12,%%ymm8,%%ymm8			\n\t"\
		"vsubpd		%%ymm7,%%ymm3,%%ymm3				\n\t		vsubpd		%%ymm13,%%ymm9,%%ymm9			\n\t"\
	"vfmadd132pd	(%%rdi),%%ymm2,%%ymm6				\n\t	vfmadd132pd		(%%rdi),%%ymm8 ,%%ymm12		\n\t"\
	"vfmadd132pd	(%%rdi),%%ymm3,%%ymm7				\n\t	vfmadd132pd		(%%rdi),%%ymm9 ,%%ymm13		\n\t"\
		"vmovaps		%%ymm2,      (%%rcx)			\n\t		vmovaps		%%ymm8 ,     (%%r12)				\n\t"\
		"vmovaps		%%ymm3, 0x020(%%rcx)			\n\t		vmovaps		%%ymm9 ,0x020(%%r12)				\n\t"\
		"vmovaps		%%ymm6,      (%%rax)			\n\t		vmovaps		%%ymm12,     (%%r10)				\n\t"\
		"vmovaps		%%ymm7, 0x020(%%rax)			\n\t		vmovaps		%%ymm13,0x020(%%r10)				\n\t"\
		"vsubpd		%%ymm5,%%ymm0,%%ymm0				\n\t		vsubpd		%%ymm15,%%ymm10,%%ymm10		\n\t"\
		"vsubpd		%%ymm4,%%ymm1,%%ymm1				\n\t		vsubpd		%%ymm14,%%ymm11,%%ymm11		\n\t"\
	"vfmadd132pd	(%%rdi),%%ymm0,%%ymm5				\n\t	vfmadd132pd		(%%rdi),%%ymm10,%%ymm15		\n\t"\
	"vfmadd132pd	(%%rdi),%%ymm1,%%ymm4				\n\t	vfmadd132pd		(%%rdi),%%ymm11,%%ymm14		\n\t"\
		"vmovaps		%%ymm0,      (%%rbx)			\n\t		vmovaps		%%ymm10,     (%%r11)				\n\t"\
		"vmovaps		%%ymm1, 0x020(%%rdx)			\n\t		vmovaps		%%ymm11,0x020(%%r13)				\n\t"\
		"vmovaps	%%ymm5,      (%%rdx)				\n\t		vmovaps		%%ymm15,     (%%r13)		\n\t"\
		"vmovaps	%%ymm4, 0x020(%%rbx)				\n\t		vmovaps		%%ymm14,0x020(%%r11)		\n\t"\
		/*...Block 2: t02,t12,t22,t32	*/				"	\n\t"	/*...Block 6: t0A,t1A,t2A,t3A */\
		"subq		$0x040,%%rax						\n\t		subq		$0x040,%%r10			\n\t"\
		"subq		$0x040,%%rbx						\n\t		subq		$0x040,%%r11			\n\t"\
		"subq		$0x040,%%rcx						\n\t		subq		$0x040,%%r12			\n\t"\
		"subq		$0x040,%%rdx						\n\t		subq		$0x040,%%r13			\n\t"\
		"addq		$0x060,%%rsi	\n\t"/* cc1 */\
		"vmovaps			 (%%rcx),%%ymm0				\n\t		vmovaps			 (%%r12),%%ymm8 				\n\t"\
		"vmovaps			 (%%rdx),%%ymm2				\n\t		vmovaps			 (%%r13),%%ymm14				\n\t"\
		"vmovaps		0x020(%%rcx),%%ymm1				\n\t		vmovaps		0x020(%%r12),%%ymm9					\n\t"\
		"vmovaps		0x020(%%rdx),%%ymm3				\n\t		vmovaps		0x020(%%r13),%%ymm11				\n\t"\
		"vmovaps			  %%ymm0,%%ymm4				\n\t		vmovaps			 %%ymm8 ,%%ymm12				\n\t"\
		"vmovaps			  %%ymm2,%%ymm6				\n\t"	/*	vmovaps			 %%ymm10,%%ymm14				\n\t"*/\
		"vmovaps			  %%ymm1,%%ymm5				\n\t"	/*	vmovaps			 %%ymm9,%%ymm13				\n\t"*/\
	/*	"vmovaps			  %%ymm3,%%ymm7				\n\t		vmovaps			 %%ymm11,%%ymm15				\n\t"*/\
	/* Strategic register-usage for sincos consts here cuts loads from 16 to 4: */\
		"vmovaps			 (%%rsi),%%ymm15			\n\t		vmovaps		0x020(%%rsi),%%ymm10			\n\t"\
		"vmovaps		0x040(%%rsi),%%ymm7 			\n\t		vmovaps		0x060(%%rsi),%%ymm13			\n\t"\
		"vmulpd			%%ymm10,%%ymm1,%%ymm1				\n\t		vmulpd			%%ymm7 ,%%ymm9 ,%%ymm9		\n\t"\
		"vmulpd			%%ymm13,%%ymm3,%%ymm3				\n\t		vmulpd			%%ymm10,%%ymm11,%%ymm11		\n\t"\
		"vmulpd			%%ymm10,%%ymm0,%%ymm0				\n\t		vmulpd			%%ymm7 ,%%ymm8 ,%%ymm8		\n\t"\
		"vmulpd			%%ymm13,%%ymm2,%%ymm2				\n\t		vmulpd			(%%r13),%%ymm10,%%ymm10		\n\t"\
	"vfmsub132pd		%%ymm15,%%ymm1,%%ymm4				\n\t	vfmsub132pd			%%ymm13,%%ymm9 ,%%ymm12		\n\t"\
	"vfmsub132pd		%%ymm7 ,%%ymm3,%%ymm6				\n\t	vfmadd132pd			%%ymm15,%%ymm11,%%ymm14		\n\t"\
	"vfmadd132pd		%%ymm15,%%ymm0,%%ymm5				\n\t	vfmadd132pd		0x020(%%r12),%%ymm8 ,%%ymm13	\n\t"\
	"vfmadd132pd	0x020(%%rdx),%%ymm2,%%ymm7				\n\t	vfmsub132pd		0x020(%%r13),%%ymm10,%%ymm15	\n\t"\
		"vsubpd		%%ymm6,%%ymm4,%%ymm4			\n\t		vsubpd		%%ymm14,%%ymm12,%%ymm12		\n\t"\
		"vsubpd		%%ymm7,%%ymm5,%%ymm5			\n\t		vsubpd		%%ymm15,%%ymm13,%%ymm13		\n\t"\
		"vmovaps			  (%%rbx),%%ymm2				\n\t		vmovaps			 (%%r11),%%ymm10				\n\t"\
		"vmovaps		 0x020(%%rbx),%%ymm0				\n\t		vmovaps		0x020(%%r11),%%ymm8					\n\t"\
	"vfmadd132pd	(%%rdi),%%ymm4,%%ymm6				\n\t	vfmadd132pd		(%%rdi),%%ymm12,%%ymm14		\n\t"\
	"vfmadd132pd	(%%rdi),%%ymm5,%%ymm7				\n\t	vfmadd132pd		(%%rdi),%%ymm13,%%ymm15		\n\t"\
		"vmovaps			  (%%rbx),%%ymm1				\n\t		vmovaps			 (%%r11),%%ymm9					\n\t"\
		"vmovaps		-0x040(%%rsi),%%ymm3				\n\t		vmovaps		-0x20(%%rsi),%%ymm11				\n\t"\
		"vmulpd			%%ymm11,%%ymm0,%%ymm0			\n\t		vmulpd			%%ymm3 ,%%ymm8 ,%%ymm8		\n\t"\
		"vmulpd			%%ymm11,%%ymm1,%%ymm1			\n\t		vmulpd			%%ymm3 ,%%ymm9 ,%%ymm9		\n\t"\
	"vfmsub132pd		%%ymm3 ,%%ymm0,%%ymm2			\n\t	vfmadd132pd			%%ymm11,%%ymm8 ,%%ymm10		\n\t"\
	"vfmadd132pd	 0x20(%%rbx),%%ymm1,%%ymm3			\n\t	vfmsub132pd		0x20(%%r11),%%ymm9 ,%%ymm11		\n\t"\
		"vmovaps			  (%%rax),%%ymm0				\n\t		vmovaps			 (%%r10),%%ymm8					\n\t"\
		"vmovaps		 0x020(%%rax),%%ymm1				\n\t		vmovaps		0x020(%%r10),%%ymm9					\n\t"\
		"vsubpd		%%ymm2,%%ymm0,%%ymm0			\n\t		vsubpd		%%ymm10,%%ymm8,%%ymm8			\n\t"\
		"vsubpd		%%ymm3,%%ymm1,%%ymm1			\n\t		vsubpd		%%ymm11,%%ymm9,%%ymm9			\n\t"\
		"vaddpd			  (%%rax),%%ymm2,%%ymm2			\n\t		vaddpd			 (%%r10),%%ymm10,%%ymm10	\n\t"\
		"vaddpd		 0x020(%%rax),%%ymm3,%%ymm3			\n\t		vaddpd		0x020(%%r10),%%ymm11,%%ymm11	\n\t"\
		"vsubpd		%%ymm6,%%ymm2,%%ymm2			\n\t		vsubpd		%%ymm12,%%ymm8,%%ymm8			\n\t"\
		"vsubpd		%%ymm7,%%ymm3,%%ymm3			\n\t		vsubpd		%%ymm13,%%ymm9,%%ymm9			\n\t"\
	"vfmadd132pd	(%%rdi),%%ymm2,%%ymm6			\n\t	vfmadd132pd		(%%rdi),%%ymm8,%%ymm12		\n\t"\
	"vfmadd132pd	(%%rdi),%%ymm3,%%ymm7			\n\t	vfmadd132pd		(%%rdi),%%ymm9,%%ymm13		\n\t"\
		"vmovaps		%%ymm2,      (%%rcx)			\n\t		vmovaps		%%ymm8 ,     (%%r12)				\n\t"\
		"vmovaps		%%ymm3, 0x020(%%rcx)			\n\t		vmovaps		%%ymm9 ,0x020(%%r12)				\n\t"\
		"vmovaps		%%ymm6,      (%%rax)			\n\t		vmovaps		%%ymm12,     (%%r10)				\n\t"\
		"vmovaps		%%ymm7, 0x020(%%rax)			\n\t		vmovaps		%%ymm13,0x020(%%r10)				\n\t"\
		"vsubpd		%%ymm5,%%ymm0,%%ymm0			\n\t		vsubpd		%%ymm15,%%ymm10,%%ymm10		\n\t"\
		"vsubpd		%%ymm4,%%ymm1,%%ymm1			\n\t		vsubpd		%%ymm14,%%ymm11,%%ymm11		\n\t"\
	"vfmadd132pd	(%%rdi),%%ymm0,%%ymm5			\n\t	vfmadd132pd		(%%rdi),%%ymm10,%%ymm15		\n\t"\
	"vfmadd132pd	(%%rdi),%%ymm1,%%ymm4			\n\t	vfmadd132pd		(%%rdi),%%ymm11,%%ymm14		\n\t"\
		"vmovaps		%%ymm0,      (%%rbx)			\n\t		vmovaps		%%ymm10,     (%%r11)				\n\t"\
		"vmovaps		%%ymm1, 0x020(%%rdx)			\n\t		vmovaps		%%ymm11,0x020(%%r13)				\n\t"\
		"vmovaps		%%ymm5,      (%%rdx)			\n\t		vmovaps		%%ymm15,     (%%r13)		\n\t"\
		"vmovaps		%%ymm4, 0x020(%%rbx)			\n\t		vmovaps		%%ymm14,0x020(%%r11)		\n\t"\
		/*...Block 4: t06,t16,t26,t36	*/							/*...Block 8: t0E,t1E,t2E,t3E	*/\
		"addq		$0x080,%%rax						\n\t		addq		$0x080,%%r10			\n\t"\
		"addq		$0x080,%%rbx						\n\t		addq		$0x080,%%r11			\n\t"\
		"addq		$0x080,%%rcx						\n\t		addq		$0x080,%%r12			\n\t"\
		"addq		$0x080,%%rdx						\n\t		addq		$0x080,%%r13			\n\t"\
		"vmovaps			 (%%rcx),%%ymm0				\n\t		vmovaps			 (%%r12),%%ymm8 				\n\t"\
		"vmovaps			 (%%rdx),%%ymm6				\n\t		vmovaps			 (%%r13),%%ymm10				\n\t"\
		"vmovaps		0x020(%%rcx),%%ymm1				\n\t		vmovaps		0x020(%%r12),%%ymm9					\n\t"\
		"vmovaps		0x020(%%rdx),%%ymm3				\n\t		vmovaps		0x020(%%r13),%%ymm11				\n\t"\
		"vmovaps			  %%ymm0,%%ymm4				\n\t		vmovaps			 %%ymm8 ,%%ymm12				\n\t"\
	/*	"vmovaps			  %%ymm2,%%ymm6			*/"	\n\t		vmovaps			 %%ymm10,%%ymm14				\n\t"\
	/*	"vmovaps			  %%ymm1,%%ymm5			*/"	\n\t		vmovaps			 %%ymm9,%%ymm13				\n\t"\
	/*	"vmovaps			  %%ymm3,%%ymm7				\n\t		vmovaps			 %%ymm11,%%ymm15				\n\t"*/\
	/* Strategic register-usage for sincos consts here cuts loads from 16 to 4: */\
		"vmovaps			 (%%rsi),%%ymm2 			\n\t		vmovaps		0x020(%%rsi),%%ymm7 			\n\t"\
		"vmovaps		0x040(%%rsi),%%ymm5 			\n\t		vmovaps		0x060(%%rsi),%%ymm15			\n\t"\
		"vmulpd		 	%%ymm15,%%ymm1,%%ymm1				\n\t		vmulpd			%%ymm2 ,%%ymm9,%%ymm9		\n\t"\
		"vmulpd		 	%%ymm2 ,%%ymm3,%%ymm3				\n\t		vmulpd			%%ymm5 ,%%ymm11,%%ymm11		\n\t"\
		"vmulpd		 	%%ymm15,%%ymm0,%%ymm0				\n\t		vmulpd			%%ymm2 ,%%ymm8,%%ymm8		\n\t"\
		"vmulpd			 (%%rdx),%%ymm2,%%ymm2				\n\t		vmulpd			%%ymm5 ,%%ymm10,%%ymm10		\n\t"\
	"vfmsub132pd	 	%%ymm5 ,%%ymm1,%%ymm4				\n\t	vfmsub132pd			%%ymm7 ,%%ymm9 ,%%ymm12		\n\t"\
	"vfmadd132pd	 	%%ymm7 ,%%ymm3,%%ymm6				\n\t	vfmsub132pd			%%ymm15,%%ymm11,%%ymm14		\n\t"\
	"vfmadd132pd	0x020(%%rcx),%%ymm0,%%ymm5				\n\t	vfmadd132pd			%%ymm7 ,%%ymm8 ,%%ymm13		\n\t"\
	"vfmsub132pd	0x020(%%rdx),%%ymm2,%%ymm7				\n\t	vfmadd132pd		0x20(%%r13),%%ymm10,%%ymm15		\n\t"\
		"vsubpd		%%ymm6,%%ymm4,%%ymm4			\n\t		vsubpd		%%ymm14,%%ymm12,%%ymm12		\n\t"\
		"vsubpd		%%ymm7,%%ymm5,%%ymm5			\n\t		vsubpd		%%ymm15,%%ymm13,%%ymm13		\n\t"\
		"vmovaps			  (%%rbx),%%ymm2				\n\t		vmovaps			 (%%r11),%%ymm10				\n\t"\
		"vmovaps		 0x020(%%rbx),%%ymm0				\n\t		vmovaps		0x020(%%r11),%%ymm8					\n\t"\
	"vfmadd132pd	(%%rdi),%%ymm4,%%ymm6				\n\t	vfmadd132pd		(%%rdi),%%ymm12,%%ymm14		\n\t"\
	"vfmadd132pd	(%%rdi),%%ymm5,%%ymm7				\n\t	vfmadd132pd		(%%rdi),%%ymm13,%%ymm15		\n\t"\
		"vmovaps			  (%%rbx),%%ymm2				\n\t		vmovaps			 (%%r11),%%ymm10				\n\t"\
		"vmovaps		 0x020(%%rbx),%%ymm0				\n\t		vmovaps		0x020(%%r11),%%ymm8					\n\t"\
		"vmovaps			  (%%rbx),%%ymm1				\n\t		vmovaps			 (%%r11),%%ymm9					\n\t"\
		"vmovaps		 0x020(%%rbx),%%ymm3				\n\t		vmovaps		0x020(%%r11),%%ymm11				\n\t"\
		"vmovaps			  (%%rbx),%%ymm1				\n\t		vmovaps			 (%%r11),%%ymm9					\n\t"\
		"vmovaps		-0x020(%%rsi),%%ymm3				\n\t		vmovaps		-0x40(%%rsi),%%ymm11				\n\t"\
		"vmulpd			%%ymm11,%%ymm0,%%ymm0			\n\t		vmulpd			%%ymm3 ,%%ymm8 ,%%ymm8		\n\t"\
		"vmulpd			%%ymm11,%%ymm1,%%ymm1			\n\t		vmulpd			%%ymm3 ,%%ymm9 ,%%ymm9		\n\t"\
	"vfmsub132pd		%%ymm3 ,%%ymm0,%%ymm2			\n\t	vfmadd132pd			%%ymm11,%%ymm8 ,%%ymm10		\n\t"\
	"vfmadd132pd	 0x20(%%rbx),%%ymm1,%%ymm3			\n\t	vfmsub132pd		0x20(%%r11),%%ymm9 ,%%ymm11		\n\t"\
		"vmovaps			  (%%rax),%%ymm0				\n\t		vmovaps			 (%%r10),%%ymm8					\n\t"\
		"vmovaps		 0x020(%%rax),%%ymm1				\n\t		vmovaps		0x020(%%r10),%%ymm9					\n\t"\
		"vsubpd		%%ymm2,%%ymm0,%%ymm0			\n\t		vsubpd		%%ymm10,%%ymm8,%%ymm8			\n\t"\
		"vsubpd		%%ymm3,%%ymm1,%%ymm1			\n\t		vsubpd		%%ymm11,%%ymm9,%%ymm9			\n\t"\
		"vaddpd			  (%%rax),%%ymm2,%%ymm2			\n\t		vaddpd			 (%%r10),%%ymm10,%%ymm10	\n\t"\
		"vaddpd		 0x020(%%rax),%%ymm3,%%ymm3			\n\t		vaddpd		0x020(%%r10),%%ymm11,%%ymm11	\n\t"\
		"vsubpd		%%ymm4,%%ymm2,%%ymm2			\n\t		vsubpd		%%ymm12,%%ymm8,%%ymm8			\n\t"\
		"vsubpd		%%ymm5,%%ymm3,%%ymm3			\n\t		vsubpd		%%ymm13,%%ymm9,%%ymm9			\n\t"\
	"vfmadd132pd	(%%rdi),%%ymm2,%%ymm4			\n\t	vfmadd132pd		(%%rdi),%%ymm8,%%ymm12		\n\t"\
	"vfmadd132pd	(%%rdi),%%ymm3,%%ymm5			\n\t	vfmadd132pd		(%%rdi),%%ymm9,%%ymm13		\n\t"\
		"vmovaps		%%ymm2,      (%%rcx)			\n\t		vmovaps		%%ymm8 ,     (%%r12)				\n\t"\
		"vmovaps		%%ymm3, 0x020(%%rcx)			\n\t		vmovaps		%%ymm9 ,0x020(%%r12)				\n\t"\
		"vmovaps		%%ymm4,      (%%rax)			\n\t		vmovaps		%%ymm12,     (%%r10)				\n\t"\
		"vmovaps		%%ymm5, 0x020(%%rax)			\n\t		vmovaps		%%ymm13,0x020(%%r10)				\n\t"\
		"vsubpd		%%ymm7,%%ymm0,%%ymm0			\n\t		vsubpd		%%ymm15,%%ymm10,%%ymm10		\n\t"\
		"vsubpd		%%ymm6,%%ymm1,%%ymm1			\n\t		vsubpd		%%ymm14,%%ymm11,%%ymm11		\n\t"\
	"vfmadd132pd	(%%rdi),%%ymm0,%%ymm7			\n\t	vfmadd132pd		(%%rdi),%%ymm10,%%ymm15		\n\t"\
	"vfmadd132pd	(%%rdi),%%ymm1,%%ymm6			\n\t	vfmadd132pd		(%%rdi),%%ymm11,%%ymm14		\n\t"\
		"vmovaps		%%ymm0,      (%%rbx)			\n\t		vmovaps		%%ymm10,     (%%r11)				\n\t"\
		"vmovaps		%%ymm1, 0x020(%%rdx)			\n\t		vmovaps		%%ymm11,0x020(%%r13)				\n\t"\
		"vmovaps		%%ymm7,      (%%rdx)			\n\t		vmovaps		%%ymm15,     (%%r13)		\n\t"\
		"vmovaps		%%ymm6, 0x020(%%rbx)			\n\t		vmovaps		%%ymm14,0x020(%%r11)		\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__add1] "m" (Xadd1)\
		 ,[__add2] "m" (Xadd2)\
		 ,[__add3] "m" (Xadd3)\
		 ,[__r00] "m" (Xr00)\
		 ,[__r10] "m" (Xr10)\
		 ,[__r20] "m" (Xr20)\
		 ,[__r30] "m" (Xr30)\
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__cc0] "m" (Xcc0)\
		 ,[__c00] "m" (Xc00)\
		 ,[__c01] "m" (Xc01)\
		 ,[__c02] "m" (Xc02)\
		 ,[__c03] "m" (Xc03)\
		 ,[__c05] "m" (Xc05)\
		 ,[__c07] "m" (Xc07)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r10","r11","r12","r13","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"		/* Clobbered registers */\
	);\
	}

	// FMA version: replace [327 ADD, 219 SUB, 212 MUL, 934 memref] ==> [70 ADD, 170 SUB, 116 MUL, 202 FMA, 996 memref].
	//
	#define SSE2_RADIX32_WRAPPER_DIT(Xadd0,Xadd1,Xadd2,Xadd3,Xisrt2,Xr00,Xr08,Xr10,Xr20,Xr28,Xr30,Xc01,Xc02,Xc04,Xc06,Xc08,Xc0A,Xc0C,Xc0E,Xc10,Xc12,Xc14,Xc16,Xc18,Xc1A,Xc1C,Xc1E)\
	{\
	__asm__ volatile (\
	/************************************************************************/\
	/*...And do an inverse DIT radix-32 pass on the squared-data blocks.	*/\
	/************************************************************************/\
	/*...Block 1: */\
		"movq	%[__isrt2],%%rsi					\n\t		leaq	0x8e0(%%rsi),%%rdi	/* two */	\n\t"\
		"movq		%[__r00],%%rax			\n\t"\
		"leaq		0x400(%%rax),%%rbx				\n\t"\
		"leaq		0x200(%%rax),%%rcx				\n\t"\
		"leaq		0x600(%%rax),%%rdx				\n\t"\
		/*****	SSE2_RADIX4_DIT_IN_PLACE()	*****/				/*...Block 2 has tmp-addresses offset +0x80 w.r.to Block 1:	*/\
		"vmovaps		(%%rax),%%ymm0				\n\t		vmovaps		0x80(%%rax),%%ymm8 			\n\t"\
		"vmovaps	0x20(%%rax),%%ymm1				\n\t		vmovaps		0xa0(%%rax),%%ymm9 			\n\t"\
		"vmovaps		(%%rax),%%ymm2				\n\t		vmovaps		0x80(%%rax),%%ymm10			\n\t"\
		"vmovaps	0x20(%%rax),%%ymm3				\n\t		vmovaps		0xa0(%%rax),%%ymm11			\n\t"\
		"vaddpd			(%%rbx),%%ymm0,%%ymm0		\n\t		vaddpd		0x80(%%rbx),%%ymm8 ,%%ymm8 	\n\t"\
		"vaddpd		0x20(%%rbx),%%ymm1,%%ymm1		\n\t		vaddpd		0xa0(%%rbx),%%ymm9 ,%%ymm9 	\n\t"\
		"vsubpd			(%%rbx),%%ymm2,%%ymm2		\n\t		vsubpd		0x80(%%rbx),%%ymm10,%%ymm10	\n\t"\
		"vsubpd		0x20(%%rbx),%%ymm3,%%ymm3		\n\t		vsubpd		0xa0(%%rbx),%%ymm11,%%ymm11	\n\t"\
		"vmovaps		(%%rcx),%%ymm4				\n\t		vmovaps		0x80(%%rcx),%%ymm12			\n\t"\
		"vmovaps	0x20(%%rcx),%%ymm5				\n\t		vmovaps		0xa0(%%rcx),%%ymm13			\n\t"\
		"vmovaps		(%%rcx),%%ymm6				\n\t		vmovaps		0x80(%%rcx),%%ymm14			\n\t"\
		"vmovaps	0x20(%%rcx),%%ymm7				\n\t		vmovaps		0xa0(%%rcx),%%ymm15			\n\t"\
		"vaddpd			(%%rdx),%%ymm4,%%ymm4		\n\t		vaddpd		0x80(%%rdx),%%ymm12,%%ymm12	\n\t"\
		"vaddpd		0x20(%%rdx),%%ymm5,%%ymm5		\n\t		vaddpd		0xa0(%%rdx),%%ymm13,%%ymm13	\n\t"\
		"vsubpd			(%%rdx),%%ymm6,%%ymm6		\n\t		vsubpd		0x80(%%rdx),%%ymm14,%%ymm14	\n\t"\
		"vsubpd		0x20(%%rdx),%%ymm7,%%ymm7		\n\t		vsubpd		0xa0(%%rdx),%%ymm15,%%ymm15	\n\t"\
		"vsubpd			%%ymm4 ,%%ymm0,%%ymm0		\n\t		vsubpd			%%ymm12,%%ymm8 ,%%ymm8 	\n\t"\
		"vsubpd			%%ymm5 ,%%ymm1,%%ymm1		\n\t		vsubpd			%%ymm13,%%ymm9 ,%%ymm9 	\n\t"\
		"vsubpd			%%ymm7 ,%%ymm2,%%ymm2		\n\t		vsubpd			%%ymm15,%%ymm10,%%ymm10	\n\t"\
		"vsubpd			%%ymm6 ,%%ymm3,%%ymm3		\n\t		vsubpd			%%ymm14,%%ymm11,%%ymm11	\n\t"\
		"vmovaps		%%ymm0 ,      (%%rbx)		\n\t		vmovaps			%%ymm8 , 0x080(%%rbx)	\n\t"\
		"vmovaps		%%ymm1 , 0x020(%%rbx)		\n\t		vmovaps			%%ymm9 , 0x0a0(%%rbx)	\n\t"\
		"vmovaps		%%ymm2 ,      (%%rdx)		\n\t		vmovaps			%%ymm10, 0x080(%%rdx)	\n\t"\
		"vmovaps		%%ymm3 , 0x020(%%rcx)		\n\t		vmovaps			%%ymm11, 0x0a0(%%rcx)	\n\t"\
	"vmovaps	%%ymm14,0xa0(%%rdx)	\n\t"/* spill ymm14 to make room for two */"	vmovaps	(%%rdi),%%ymm14	\n\t"/* two */\
	"vfmadd132pd	%%ymm14,%%ymm0,%%ymm4			\n\t	vfmadd132pd		%%ymm14,%%ymm8 ,%%ymm12		\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm1,%%ymm5			\n\t	vfmadd132pd		%%ymm14,%%ymm9 ,%%ymm13		\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm2,%%ymm7			\n\t	vfmadd132pd		%%ymm14,%%ymm10,%%ymm15		\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm3,%%ymm6			\n\t	vfmadd132pd	0xa0(%%rdx),%%ymm11,%%ymm14		\n\t"\
		"vmovaps	%%ymm4 ,    (%%rax)				\n\t		vmovaps		%%ymm12,0x80(%%rax)			\n\t"\
		"vmovaps	%%ymm5 ,0x20(%%rax)				\n\t		vmovaps		%%ymm13,0xa0(%%rax)			\n\t"\
		"vmovaps	%%ymm7 ,    (%%rcx)				\n\t		vmovaps		%%ymm15,0x80(%%rcx)			\n\t"\
		"vmovaps	%%ymm6 ,0x20(%%rdx)				\n\t		vmovaps		%%ymm14,0xa0(%%rdx)			\n\t"\
		"addq		$0x100,%%rax				\n\t"\
		"addq		$0x100,%%rbx				\n\t"\
		"addq		$0x100,%%rcx				\n\t"\
		"addq		$0x100,%%rdx				\n\t"\
		"vmovaps		(%%rax),%%ymm0				\n\t		vmovaps		0x80(%%rax),%%ymm8 			\n\t"\
		"vmovaps	0x20(%%rax),%%ymm1				\n\t		vmovaps		0xa0(%%rax),%%ymm9 			\n\t"\
		"vmovaps		(%%rax),%%ymm2				\n\t		vmovaps		0x80(%%rax),%%ymm10			\n\t"\
		"vmovaps	0x20(%%rax),%%ymm3				\n\t		vmovaps		0xa0(%%rax),%%ymm11			\n\t"\
		"vaddpd			(%%rbx),%%ymm0,%%ymm0		\n\t		vaddpd		0x80(%%rbx),%%ymm8 ,%%ymm8 	\n\t"\
		"vaddpd		0x20(%%rbx),%%ymm1,%%ymm1		\n\t		vaddpd		0xa0(%%rbx),%%ymm9 ,%%ymm9 	\n\t"\
		"vsubpd			(%%rbx),%%ymm2,%%ymm2		\n\t		vsubpd		0x80(%%rbx),%%ymm10,%%ymm10	\n\t"\
		"vsubpd		0x20(%%rbx),%%ymm3,%%ymm3		\n\t		vsubpd		0xa0(%%rbx),%%ymm11,%%ymm11	\n\t"\
		"vmovaps		(%%rcx),%%ymm4				\n\t		vmovaps		0x80(%%rcx),%%ymm12			\n\t"\
		"vmovaps	0x20(%%rcx),%%ymm5				\n\t		vmovaps		0xa0(%%rcx),%%ymm13			\n\t"\
		"vmovaps		(%%rcx),%%ymm6				\n\t		vmovaps		0x80(%%rcx),%%ymm14			\n\t"\
		"vmovaps	0x20(%%rcx),%%ymm7				\n\t		vmovaps		0xa0(%%rcx),%%ymm15			\n\t"\
		"vaddpd			(%%rdx),%%ymm4,%%ymm4		\n\t		vaddpd		0x80(%%rdx),%%ymm12,%%ymm12	\n\t"\
		"vaddpd		0x20(%%rdx),%%ymm5,%%ymm5		\n\t		vaddpd		0xa0(%%rdx),%%ymm13,%%ymm13	\n\t"\
		"vsubpd			(%%rdx),%%ymm6,%%ymm6		\n\t		vsubpd		0x80(%%rdx),%%ymm14,%%ymm14	\n\t"\
		"vsubpd		0x20(%%rdx),%%ymm7,%%ymm7		\n\t		vsubpd		0xa0(%%rdx),%%ymm15,%%ymm15	\n\t"\
		"vsubpd			%%ymm4 ,%%ymm0,%%ymm0		\n\t		vsubpd			%%ymm12,%%ymm8 ,%%ymm8 	\n\t"\
		"vsubpd			%%ymm5 ,%%ymm1,%%ymm1		\n\t		vsubpd			%%ymm13,%%ymm9 ,%%ymm9 	\n\t"\
		"vsubpd			%%ymm7 ,%%ymm2,%%ymm2		\n\t		vsubpd			%%ymm15,%%ymm10,%%ymm10	\n\t"\
		"vsubpd			%%ymm6 ,%%ymm3,%%ymm3		\n\t		vsubpd			%%ymm14,%%ymm11,%%ymm11	\n\t"\
	"vmovaps	%%ymm14,0xa0(%%rdx)	\n\t"/* spill ymm14 to make room for two */"	vmovaps	(%%rdi),%%ymm14	\n\t"/* two */\
	"vfmadd132pd	%%ymm14,%%ymm0,%%ymm4			\n\t	vfmadd132pd		%%ymm14,%%ymm8 ,%%ymm12		\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm1,%%ymm5			\n\t	vfmadd132pd		%%ymm14,%%ymm9 ,%%ymm13		\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm2,%%ymm7			\n\t	vfmadd132pd		%%ymm14,%%ymm10,%%ymm15		\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm3,%%ymm6			\n\t	vfmadd132pd	0xa0(%%rdx),%%ymm11,%%ymm14		\n\t"\
		"vmovaps	%%ymm0 ,      (%%rbx)			\n\t		vmovaps			%%ymm8 , 0x080(%%rbx)	\n\t"\
		"vmovaps	%%ymm1 , 0x020(%%rbx)			\n\t		vmovaps			%%ymm9 , 0x0a0(%%rbx)	\n\t"\
		"vmovaps	%%ymm4 ,      (%%rax)			\n\t		vmovaps			%%ymm12, 0x080(%%rax)	\n\t"\
		"vmovaps	%%ymm5 , 0x020(%%rax)			\n\t		vmovaps			%%ymm13, 0x0a0(%%rax)	\n\t"\
		"vmovaps	%%ymm3 ,%%ymm0					\n\t		vmovaps			%%ymm11,%%ymm8 			\n\t"\
		"vmovaps	%%ymm6 ,%%ymm1					\n\t		vmovaps			%%ymm14,%%ymm9 			\n\t"\
		"vsubpd		%%ymm7 ,%%ymm3,%%ymm3			\n\t		vsubpd			%%ymm15,%%ymm11,%%ymm11	\n\t"\
		"vsubpd		%%ymm2 ,%%ymm6,%%ymm6			\n\t		vsubpd			%%ymm10,%%ymm14,%%ymm14	\n\t"\
		"vaddpd		%%ymm7 ,%%ymm0,%%ymm0			\n\t		vaddpd			%%ymm15,%%ymm8 ,%%ymm8 	\n\t"\
		"vaddpd		%%ymm2 ,%%ymm1,%%ymm1			\n\t		vaddpd			%%ymm10,%%ymm9 ,%%ymm9 	\n\t"\
	"vmovaps	(%%rsi),%%ymm15	\n\t"/* isrt2 */\
		"vmulpd		%%ymm15,%%ymm3,%%ymm3			\n\t		vmulpd			%%ymm15,%%ymm11,%%ymm11	\n\t"\
		"vmulpd		%%ymm15,%%ymm6,%%ymm6			\n\t		vmulpd			%%ymm15,%%ymm14,%%ymm14	\n\t"\
		"vmulpd		%%ymm15,%%ymm0,%%ymm0			\n\t		vmulpd			%%ymm15,%%ymm8 ,%%ymm8 	\n\t"\
		"vmulpd		%%ymm15,%%ymm1,%%ymm1			\n\t		vmulpd			%%ymm15,%%ymm9 ,%%ymm9 	\n\t"\
		"vmovaps	%%ymm3 , 0x020(%%rcx)			\n\t		vmovaps			%%ymm11, 0x0a0(%%rcx)	\n\t"\
		"vmovaps	%%ymm6 , 0x020(%%rdx)			\n\t		vmovaps			%%ymm14, 0x0a0(%%rdx)	\n\t"\
		"vmovaps	%%ymm0 ,      (%%rcx)			\n\t		vmovaps			%%ymm8 , 0x080(%%rcx)	\n\t"\
		"vmovaps	%%ymm1 ,      (%%rdx)			\n\t		vmovaps			%%ymm9 , 0x080(%%rdx)	\n\t"\
		/*****	SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS *****/			/*****	SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS *****/\
		/*        (r00,r10,r20,r30,r08,r18,r28,r38):  */			/*        (r04,r14,r24,r34,r0C,r1C,r2C,r3C):  */\
		"vmovaps	-0x100(%%rax),%%ymm0			\n\t		vmovaps		-0x080(%%rax),%%ymm8 		\n\t"\
		"vmovaps	-0x100(%%rbx),%%ymm4			\n\t		vmovaps		-0x080(%%rbx),%%ymm12		\n\t"\
		"vmovaps	-0x0e0(%%rax),%%ymm1			\n\t		vmovaps		-0x060(%%rax),%%ymm9 		\n\t"\
		"vmovaps	-0x0e0(%%rbx),%%ymm5			\n\t		vmovaps		-0x060(%%rbx),%%ymm13		\n\t"\
		"vmovaps		  (%%rax),%%ymm2			\n\t		vmovaps		 0x080(%%rax),%%ymm10		\n\t"\
		"vmovaps	 0x020(%%rbx),%%ymm7			\n\t		vmovaps		 0x0a0(%%rbx),%%ymm15		\n\t"\
		"vmovaps	 0x020(%%rax),%%ymm3			\n\t		vmovaps		 0x0a0(%%rax),%%ymm11		\n\t"\
		"vmovaps		  (%%rbx),%%ymm6			\n\t		vmovaps		 0x080(%%rbx),%%ymm14		\n\t"\
		"vsubpd		%%ymm2 ,%%ymm0,%%ymm0			\n\t		vsubpd		%%ymm10,%%ymm8 ,%%ymm8 		\n\t"\
		"vsubpd		%%ymm7 ,%%ymm4,%%ymm4			\n\t		vsubpd		%%ymm15,%%ymm12,%%ymm12		\n\t"\
		"vsubpd		%%ymm3 ,%%ymm1,%%ymm1			\n\t		vsubpd		%%ymm11,%%ymm9 ,%%ymm9 		\n\t"\
		"vsubpd		%%ymm6 ,%%ymm5,%%ymm5			\n\t		vsubpd		%%ymm14,%%ymm13,%%ymm13		\n\t"\
	"vmovaps	%%ymm14,0xa0(%%rbx)	\n\t"/* spill ymm14 to make room for two */"	vmovaps	(%%rdi),%%ymm14	\n\t"/* two */\
	"vfmadd132pd	%%ymm14,%%ymm0,%%ymm2			\n\t	vfmadd132pd		%%ymm14,%%ymm8 ,%%ymm10		\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm4,%%ymm7			\n\t	vfmadd132pd		%%ymm14,%%ymm12,%%ymm15		\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm1,%%ymm3			\n\t	vfmadd132pd		%%ymm14,%%ymm9 ,%%ymm11		\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm5,%%ymm6			\n\t	vfmadd132pd	0xa0(%%rbx),%%ymm13,%%ymm14		\n\t"\
		"vmovaps	%%ymm0,      (%%rax)			\n\t		vmovaps		%%ymm8 , 0x080(%%rax)		\n\t"\
		"vmovaps	%%ymm4,      (%%rbx)			\n\t		vmovaps		%%ymm12, 0x080(%%rbx)		\n\t"\
		"vmovaps	%%ymm1, 0x020(%%rax)			\n\t		vmovaps		%%ymm9 , 0x0a0(%%rax)		\n\t"\
		"vmovaps	%%ymm5,-0x0e0(%%rbx)			\n\t		vmovaps		%%ymm13,-0x060(%%rbx)		\n\t"\
		"vmovaps	%%ymm2,-0x100(%%rax)			\n\t		vmovaps		%%ymm10,-0x080(%%rax)		\n\t"\
		"vmovaps	%%ymm7,-0x100(%%rbx)			\n\t		vmovaps		%%ymm15,-0x080(%%rbx)		\n\t"\
		"vmovaps	%%ymm3,-0x0e0(%%rax)			\n\t		vmovaps		%%ymm11,-0x060(%%rax)		\n\t"\
		"vmovaps	%%ymm6, 0x020(%%rbx)			\n\t		vmovaps		%%ymm14, 0x0a0(%%rbx)		\n\t"\
		"vmovaps	-0x100(%%rcx),%%ymm0			\n\t		vmovaps		-0x080(%%rcx),%%ymm8 		\n\t"\
		"vmovaps	-0x100(%%rdx),%%ymm4			\n\t		vmovaps		-0x080(%%rdx),%%ymm12		\n\t"\
		"vmovaps	-0x0e0(%%rcx),%%ymm1			\n\t		vmovaps		-0x060(%%rcx),%%ymm9 		\n\t"\
		"vmovaps	-0x0e0(%%rdx),%%ymm5			\n\t		vmovaps		-0x060(%%rdx),%%ymm13		\n\t"\
		"vmovaps		  (%%rcx),%%ymm2			\n\t		vmovaps		 0x080(%%rcx),%%ymm10		\n\t"\
		"vmovaps	 0x020(%%rdx),%%ymm7			\n\t		vmovaps		 0x0a0(%%rdx),%%ymm15		\n\t"\
		"vmovaps	 0x020(%%rcx),%%ymm3			\n\t		vmovaps		 0x0a0(%%rcx),%%ymm11		\n\t"\
		"vmovaps		  (%%rdx),%%ymm6			\n\t		vmovaps		 0x080(%%rdx),%%ymm14		\n\t"\
		"vsubpd		%%ymm2 ,%%ymm0,%%ymm0			\n\t		vsubpd		%%ymm10,%%ymm8 ,%%ymm8 		\n\t"\
		"vsubpd		%%ymm7 ,%%ymm4,%%ymm4			\n\t		vsubpd		%%ymm15,%%ymm12,%%ymm12		\n\t"\
		"vsubpd		%%ymm3 ,%%ymm1,%%ymm1			\n\t		vsubpd		%%ymm11,%%ymm9 ,%%ymm9 		\n\t"\
		"vsubpd		%%ymm6 ,%%ymm5,%%ymm5			\n\t		vsubpd		%%ymm14,%%ymm13,%%ymm13		\n\t"\
	"vmovaps	%%ymm14,0xa0(%%rdx)	\n\t"/* spill ymm14 to make room for two */"	vmovaps	(%%rdi),%%ymm14	\n\t"/* two */\
	"vfmadd132pd	%%ymm14,%%ymm0,%%ymm2			\n\t	vfmadd132pd		%%ymm14,%%ymm8 ,%%ymm10		\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm4,%%ymm7			\n\t	vfmadd132pd		%%ymm14,%%ymm12,%%ymm15		\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm1,%%ymm3			\n\t	vfmadd132pd		%%ymm14,%%ymm9 ,%%ymm11		\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm5,%%ymm6			\n\t	vfmadd132pd	0xa0(%%rdx),%%ymm13,%%ymm14		\n\t"\
		"vmovaps	%%ymm0,      (%%rcx)			\n\t		vmovaps		%%ymm8 , 0x080(%%rcx)		\n\t"\
		"vmovaps	%%ymm4,      (%%rdx)			\n\t		vmovaps		%%ymm12, 0x080(%%rdx)		\n\t"\
		"vmovaps	%%ymm1, 0x020(%%rcx)			\n\t		vmovaps		%%ymm9 , 0x0a0(%%rcx)		\n\t"\
		"vmovaps	%%ymm5,-0x0e0(%%rdx)			\n\t		vmovaps		%%ymm13,-0x060(%%rdx)		\n\t"\
		"vmovaps	%%ymm2,-0x100(%%rcx)			\n\t		vmovaps		%%ymm10,-0x080(%%rcx)		\n\t"\
		"vmovaps	%%ymm7,-0x100(%%rdx)			\n\t		vmovaps		%%ymm15,-0x080(%%rdx)		\n\t"\
		"vmovaps	%%ymm3,-0x0e0(%%rcx)			\n\t		vmovaps		%%ymm11,-0x060(%%rcx)		\n\t"\
		"vmovaps	%%ymm6, 0x020(%%rdx)			\n\t		vmovaps		%%ymm14, 0x0a0(%%rdx)		\n\t"\
	/*...Blocks 3,4 have tmp-addresses offset +0x40 w.r.to Blocks 1,2, respectively (thus +0x100-0x0c0 = +0x040: */\
		"subq		$0xc0,%%rax				\n\t"\
		"subq		$0xc0,%%rbx				\n\t"\
		"subq		$0xc0,%%rcx				\n\t"\
		"subq		$0xc0,%%rdx				\n\t"\
		"vmovaps		(%%rax),%%ymm0				\n\t		vmovaps		0x80(%%rax),%%ymm8 			\n\t"\
		"vmovaps	0x20(%%rax),%%ymm1				\n\t		vmovaps		0xa0(%%rax),%%ymm9 			\n\t"\
		"vmovaps		(%%rax),%%ymm2				\n\t		vmovaps		0x80(%%rax),%%ymm10			\n\t"\
		"vmovaps	0x20(%%rax),%%ymm3				\n\t		vmovaps		0xa0(%%rax),%%ymm11			\n\t"\
		"vaddpd			(%%rbx),%%ymm0,%%ymm0		\n\t		vaddpd		0x80(%%rbx),%%ymm8 ,%%ymm8 	\n\t"\
		"vaddpd		0x20(%%rbx),%%ymm1,%%ymm1		\n\t		vaddpd		0xa0(%%rbx),%%ymm9 ,%%ymm9 	\n\t"\
		"vsubpd			(%%rbx),%%ymm2,%%ymm2		\n\t		vsubpd		0x80(%%rbx),%%ymm10,%%ymm10	\n\t"\
		"vsubpd		0x20(%%rbx),%%ymm3,%%ymm3		\n\t		vsubpd		0xa0(%%rbx),%%ymm11,%%ymm11	\n\t"\
		"vmovaps		(%%rcx),%%ymm4				\n\t		vmovaps		0x80(%%rcx),%%ymm12			\n\t"\
		"vmovaps	0x20(%%rcx),%%ymm5				\n\t		vmovaps		0xa0(%%rcx),%%ymm13			\n\t"\
		"vmovaps		(%%rcx),%%ymm6				\n\t		vmovaps		0x80(%%rcx),%%ymm14			\n\t"\
		"vmovaps	0x20(%%rcx),%%ymm7				\n\t		vmovaps		0xa0(%%rcx),%%ymm15			\n\t"\
		"vaddpd			(%%rdx),%%ymm4,%%ymm4		\n\t		vaddpd		0x80(%%rdx),%%ymm12,%%ymm12	\n\t"\
		"vaddpd		0x20(%%rdx),%%ymm5,%%ymm5		\n\t		vaddpd		0xa0(%%rdx),%%ymm13,%%ymm13	\n\t"\
		"vsubpd			(%%rdx),%%ymm6,%%ymm6		\n\t		vsubpd		0x80(%%rdx),%%ymm14,%%ymm14	\n\t"\
		"vsubpd		0x20(%%rdx),%%ymm7,%%ymm7		\n\t		vsubpd		0xa0(%%rdx),%%ymm15,%%ymm15	\n\t"\
		"vsubpd			%%ymm4 ,%%ymm0,%%ymm0		\n\t		vsubpd			%%ymm12,%%ymm8 ,%%ymm8 	\n\t"\
		"vsubpd			%%ymm5 ,%%ymm1,%%ymm1		\n\t		vsubpd			%%ymm13,%%ymm9 ,%%ymm9 	\n\t"\
		"vsubpd			%%ymm7 ,%%ymm2,%%ymm2		\n\t		vsubpd			%%ymm15,%%ymm10,%%ymm10	\n\t"\
		"vsubpd			%%ymm6 ,%%ymm3,%%ymm3		\n\t		vsubpd			%%ymm14,%%ymm11,%%ymm11	\n\t"\
		"vmovaps		%%ymm0 ,      (%%rbx)		\n\t		vmovaps			%%ymm8 , 0x080(%%rbx)	\n\t"\
		"vmovaps		%%ymm1 , 0x020(%%rbx)		\n\t		vmovaps			%%ymm9 , 0x0a0(%%rbx)	\n\t"\
		"vmovaps		%%ymm2 ,      (%%rdx)		\n\t		vmovaps			%%ymm10, 0x080(%%rdx)	\n\t"\
		"vmovaps		%%ymm3 , 0x020(%%rcx)		\n\t		vmovaps			%%ymm11, 0x0a0(%%rcx)	\n\t"\
	"vmovaps	%%ymm14,0xa0(%%rdx)	\n\t"/* spill ymm14 to make room for two */"	vmovaps	(%%rdi),%%ymm14	\n\t"/* two */\
	"vfmadd132pd	%%ymm14,%%ymm0,%%ymm4			\n\t	vfmadd132pd		%%ymm14,%%ymm8 ,%%ymm12		\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm1,%%ymm5			\n\t	vfmadd132pd		%%ymm14,%%ymm9 ,%%ymm13		\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm2,%%ymm7			\n\t	vfmadd132pd		%%ymm14,%%ymm10,%%ymm15		\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm3,%%ymm6			\n\t	vfmadd132pd	0xa0(%%rdx),%%ymm11,%%ymm14		\n\t"\
		"vmovaps	%%ymm4 ,    (%%rax)				\n\t		vmovaps		%%ymm12,0x80(%%rax)			\n\t"\
		"vmovaps	%%ymm5 ,0x20(%%rax)				\n\t		vmovaps		%%ymm13,0xa0(%%rax)			\n\t"\
		"vmovaps	%%ymm7 ,    (%%rcx)				\n\t		vmovaps		%%ymm15,0x80(%%rcx)			\n\t"\
		"vmovaps	%%ymm6 ,0x20(%%rdx)				\n\t		vmovaps		%%ymm14,0xa0(%%rdx)			\n\t"\
		"addq		$0x100,%%rax				\n\t"\
		"addq		$0x100,%%rbx				\n\t"\
		"addq		$0x100,%%rcx				\n\t"\
		"addq		$0x100,%%rdx				\n\t"\
		"vmovaps		(%%rax),%%ymm0				\n\t		vmovaps		0x80(%%rax),%%ymm8 			\n\t"\
		"vmovaps	0x20(%%rax),%%ymm1				\n\t		vmovaps		0xa0(%%rax),%%ymm9 			\n\t"\
		"vmovaps		(%%rax),%%ymm2				\n\t		vmovaps		0x80(%%rax),%%ymm10			\n\t"\
		"vmovaps	0x20(%%rax),%%ymm3				\n\t		vmovaps		0xa0(%%rax),%%ymm11			\n\t"\
		"vaddpd			(%%rbx),%%ymm0,%%ymm0		\n\t		vaddpd		0x80(%%rbx),%%ymm8 ,%%ymm8 	\n\t"\
		"vaddpd		0x20(%%rbx),%%ymm1,%%ymm1		\n\t		vaddpd		0xa0(%%rbx),%%ymm9 ,%%ymm9 	\n\t"\
		"vsubpd			(%%rbx),%%ymm2,%%ymm2		\n\t		vsubpd		0x80(%%rbx),%%ymm10,%%ymm10	\n\t"\
		"vsubpd		0x20(%%rbx),%%ymm3,%%ymm3		\n\t		vsubpd		0xa0(%%rbx),%%ymm11,%%ymm11	\n\t"\
		"vmovaps		(%%rcx),%%ymm4				\n\t		vmovaps		0x80(%%rcx),%%ymm12			\n\t"\
		"vmovaps	0x20(%%rcx),%%ymm5				\n\t		vmovaps		0xa0(%%rcx),%%ymm13			\n\t"\
		"vmovaps		(%%rcx),%%ymm6				\n\t		vmovaps		0x80(%%rcx),%%ymm14			\n\t"\
		"vmovaps	0x20(%%rcx),%%ymm7				\n\t		vmovaps		0xa0(%%rcx),%%ymm15			\n\t"\
		"vaddpd			(%%rdx),%%ymm4,%%ymm4		\n\t		vaddpd		0x80(%%rdx),%%ymm12,%%ymm12	\n\t"\
		"vaddpd		0x20(%%rdx),%%ymm5,%%ymm5		\n\t		vaddpd		0xa0(%%rdx),%%ymm13,%%ymm13	\n\t"\
		"vsubpd			(%%rdx),%%ymm6,%%ymm6		\n\t		vsubpd		0x80(%%rdx),%%ymm14,%%ymm14	\n\t"\
		"vsubpd		0x20(%%rdx),%%ymm7,%%ymm7		\n\t		vsubpd		0xa0(%%rdx),%%ymm15,%%ymm15	\n\t"\
		"vsubpd			%%ymm4 ,%%ymm0,%%ymm0		\n\t		vsubpd			%%ymm12,%%ymm8 ,%%ymm8 	\n\t"\
		"vsubpd			%%ymm5 ,%%ymm1,%%ymm1		\n\t		vsubpd			%%ymm13,%%ymm9 ,%%ymm9 	\n\t"\
		"vsubpd			%%ymm7 ,%%ymm2,%%ymm2		\n\t		vsubpd			%%ymm15,%%ymm10,%%ymm10	\n\t"\
		"vsubpd			%%ymm6 ,%%ymm3,%%ymm3		\n\t		vsubpd			%%ymm14,%%ymm11,%%ymm11	\n\t"\
	"vmovaps	%%ymm14,0xa0(%%rdx)	\n\t"/* spill ymm14 to make room for two */"	vmovaps	(%%rdi),%%ymm14	\n\t"/* two */\
	"vfmadd132pd	%%ymm14,%%ymm0,%%ymm4			\n\t	vfmadd132pd		%%ymm14,%%ymm8 ,%%ymm12		\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm1,%%ymm5			\n\t	vfmadd132pd		%%ymm14,%%ymm9 ,%%ymm13		\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm2,%%ymm7			\n\t	vfmadd132pd		%%ymm14,%%ymm10,%%ymm15		\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm3,%%ymm6			\n\t	vfmadd132pd	0xa0(%%rdx),%%ymm11,%%ymm14		\n\t"\
		"vmovaps	%%ymm0 ,      (%%rbx)			\n\t		vmovaps			%%ymm8 , 0x080(%%rbx)	\n\t"\
		"vmovaps	%%ymm1 , 0x020(%%rbx)			\n\t		vmovaps			%%ymm9 , 0x0a0(%%rbx)	\n\t"\
		"vmovaps	%%ymm4 ,      (%%rax)			\n\t		vmovaps			%%ymm12, 0x080(%%rax)	\n\t"\
		"vmovaps	%%ymm5 , 0x020(%%rax)			\n\t		vmovaps			%%ymm13, 0x0a0(%%rax)	\n\t"\
		"vmovaps	%%ymm3 ,%%ymm0					\n\t		vmovaps			%%ymm11,%%ymm8 			\n\t"\
		"vmovaps	%%ymm6 ,%%ymm1					\n\t		vmovaps			%%ymm14,%%ymm9 			\n\t"\
		"vsubpd		%%ymm7 ,%%ymm3,%%ymm3			\n\t		vsubpd			%%ymm15,%%ymm11,%%ymm11	\n\t"\
		"vsubpd		%%ymm2 ,%%ymm6,%%ymm6			\n\t		vsubpd			%%ymm10,%%ymm14,%%ymm14	\n\t"\
		"vaddpd		%%ymm7 ,%%ymm0,%%ymm0			\n\t		vaddpd			%%ymm15,%%ymm8 ,%%ymm8 	\n\t"\
		"vaddpd		%%ymm2 ,%%ymm1,%%ymm1			\n\t		vaddpd			%%ymm10,%%ymm9 ,%%ymm9 	\n\t"\
	"vmovaps	(%%rsi),%%ymm15	\n\t"/* isrt2 */\
		"vmulpd		%%ymm15,%%ymm3,%%ymm3			\n\t		vmulpd			%%ymm15,%%ymm11,%%ymm11	\n\t"\
		"vmulpd		%%ymm15,%%ymm6,%%ymm6			\n\t		vmulpd			%%ymm15,%%ymm14,%%ymm14	\n\t"\
		"vmulpd		%%ymm15,%%ymm0,%%ymm0			\n\t		vmulpd			%%ymm15,%%ymm8 ,%%ymm8 	\n\t"\
		"vmulpd		%%ymm15,%%ymm1,%%ymm1			\n\t		vmulpd			%%ymm15,%%ymm9 ,%%ymm9 	\n\t"\
		"vmovaps	%%ymm3 , 0x020(%%rcx)			\n\t		vmovaps			%%ymm11, 0x0a0(%%rcx)	\n\t"\
		"vmovaps	%%ymm6 , 0x020(%%rdx)			\n\t		vmovaps			%%ymm14, 0x0a0(%%rdx)	\n\t"\
		"vmovaps	%%ymm0 ,      (%%rcx)			\n\t		vmovaps			%%ymm8 , 0x080(%%rcx)	\n\t"\
		"vmovaps	%%ymm1 ,      (%%rdx)			\n\t		vmovaps			%%ymm9 , 0x080(%%rdx)	\n\t"\
		/*****	SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS *****/			/*****	SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS *****/\
		/*        (r02,r12,r22,r32,r0A,r1A,r2A,r3A):  */			/*        (r06,r16,r26,r36,r0E,r1E,r2E,r3E):  */\
		"vmovaps	-0x100(%%rax),%%ymm0			\n\t		vmovaps		-0x080(%%rax),%%ymm8 		\n\t"\
		"vmovaps	-0x100(%%rbx),%%ymm4			\n\t		vmovaps		-0x080(%%rbx),%%ymm12		\n\t"\
		"vmovaps	-0x0e0(%%rax),%%ymm1			\n\t		vmovaps		-0x060(%%rax),%%ymm9 		\n\t"\
		"vmovaps	-0x0e0(%%rbx),%%ymm5			\n\t		vmovaps		-0x060(%%rbx),%%ymm13		\n\t"\
		"vmovaps		  (%%rax),%%ymm2			\n\t		vmovaps		 0x080(%%rax),%%ymm10		\n\t"\
		"vmovaps	 0x020(%%rbx),%%ymm7			\n\t		vmovaps		 0x0a0(%%rbx),%%ymm15		\n\t"\
		"vmovaps	 0x020(%%rax),%%ymm3			\n\t		vmovaps		 0x0a0(%%rax),%%ymm11		\n\t"\
		"vmovaps		  (%%rbx),%%ymm6			\n\t		vmovaps		 0x080(%%rbx),%%ymm14		\n\t"\
		"vsubpd		%%ymm2 ,%%ymm0,%%ymm0			\n\t		vsubpd		%%ymm10,%%ymm8 ,%%ymm8 		\n\t"\
		"vsubpd		%%ymm7 ,%%ymm4,%%ymm4			\n\t		vsubpd		%%ymm15,%%ymm12,%%ymm12		\n\t"\
		"vsubpd		%%ymm3 ,%%ymm1,%%ymm1			\n\t		vsubpd		%%ymm11,%%ymm9 ,%%ymm9 		\n\t"\
		"vsubpd		%%ymm6 ,%%ymm5,%%ymm5			\n\t		vsubpd		%%ymm14,%%ymm13,%%ymm13		\n\t"\
	"vmovaps	%%ymm14,0xa0(%%rbx)	\n\t"/* spill ymm14 to make room for two */"	vmovaps	(%%rdi),%%ymm14	\n\t"/* two */\
	"vfmadd132pd	%%ymm14,%%ymm0,%%ymm2			\n\t	vfmadd132pd		%%ymm14,%%ymm8 ,%%ymm10		\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm4,%%ymm7			\n\t	vfmadd132pd		%%ymm14,%%ymm12,%%ymm15		\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm1,%%ymm3			\n\t	vfmadd132pd		%%ymm14,%%ymm9 ,%%ymm11		\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm5,%%ymm6			\n\t	vfmadd132pd	0xa0(%%rbx),%%ymm13,%%ymm14		\n\t"\
		"vmovaps	%%ymm0,      (%%rax)			\n\t		vmovaps		%%ymm8 , 0x080(%%rax)		\n\t"\
		"vmovaps	%%ymm4,      (%%rbx)			\n\t		vmovaps		%%ymm12, 0x080(%%rbx)		\n\t"\
		"vmovaps	%%ymm1, 0x020(%%rax)			\n\t		vmovaps		%%ymm9 , 0x0a0(%%rax)		\n\t"\
		"vmovaps	%%ymm5,-0x0e0(%%rbx)			\n\t		vmovaps		%%ymm13,-0x060(%%rbx)		\n\t"\
		"vmovaps	%%ymm2,-0x100(%%rax)			\n\t		vmovaps		%%ymm10,-0x080(%%rax)		\n\t"\
		"vmovaps	%%ymm7,-0x100(%%rbx)			\n\t		vmovaps		%%ymm15,-0x080(%%rbx)		\n\t"\
		"vmovaps	%%ymm3,-0x0e0(%%rax)			\n\t		vmovaps		%%ymm11,-0x060(%%rax)		\n\t"\
		"vmovaps	%%ymm6, 0x020(%%rbx)			\n\t		vmovaps		%%ymm14, 0x0a0(%%rbx)		\n\t"\
		"vmovaps	-0x100(%%rcx),%%ymm0			\n\t		vmovaps		-0x080(%%rcx),%%ymm8 		\n\t"\
		"vmovaps	-0x100(%%rdx),%%ymm4			\n\t		vmovaps		-0x080(%%rdx),%%ymm12		\n\t"\
		"vmovaps	-0x0e0(%%rcx),%%ymm1			\n\t		vmovaps		-0x060(%%rcx),%%ymm9 		\n\t"\
		"vmovaps	-0x0e0(%%rdx),%%ymm5			\n\t		vmovaps		-0x060(%%rdx),%%ymm13		\n\t"\
		"vmovaps		  (%%rcx),%%ymm2			\n\t		vmovaps		 0x080(%%rcx),%%ymm10		\n\t"\
		"vmovaps	 0x020(%%rdx),%%ymm7			\n\t		vmovaps		 0x0a0(%%rdx),%%ymm15		\n\t"\
		"vmovaps	 0x020(%%rcx),%%ymm3			\n\t		vmovaps		 0x0a0(%%rcx),%%ymm11		\n\t"\
		"vmovaps		  (%%rdx),%%ymm6			\n\t		vmovaps		 0x080(%%rdx),%%ymm14		\n\t"\
		"vsubpd		%%ymm2 ,%%ymm0,%%ymm0			\n\t		vsubpd		%%ymm10,%%ymm8 ,%%ymm8 		\n\t"\
		"vsubpd		%%ymm7 ,%%ymm4,%%ymm4			\n\t		vsubpd		%%ymm15,%%ymm12,%%ymm12		\n\t"\
		"vsubpd		%%ymm3 ,%%ymm1,%%ymm1			\n\t		vsubpd		%%ymm11,%%ymm9 ,%%ymm9 		\n\t"\
		"vsubpd		%%ymm6 ,%%ymm5,%%ymm5			\n\t		vsubpd		%%ymm14,%%ymm13,%%ymm13		\n\t"\
	"vmovaps	%%ymm14,0xa0(%%rdx)	\n\t"/* spill ymm14 to make room for two */"	vmovaps	(%%rdi),%%ymm14	\n\t"/* two */\
	"vfmadd132pd	%%ymm14,%%ymm0,%%ymm2			\n\t	vfmadd132pd		%%ymm14,%%ymm8 ,%%ymm10		\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm4,%%ymm7			\n\t	vfmadd132pd		%%ymm14,%%ymm12,%%ymm15		\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm1,%%ymm3			\n\t	vfmadd132pd		%%ymm14,%%ymm9 ,%%ymm11		\n\t"\
	"vfmadd132pd	%%ymm14,%%ymm5,%%ymm6			\n\t	vfmadd132pd	0xa0(%%rdx),%%ymm13,%%ymm14		\n\t"\
		"vmovaps	%%ymm0,      (%%rcx)			\n\t		vmovaps		%%ymm8 , 0x080(%%rcx)		\n\t"\
		"vmovaps	%%ymm4,      (%%rdx)			\n\t		vmovaps		%%ymm12, 0x080(%%rdx)		\n\t"\
		"vmovaps	%%ymm1, 0x020(%%rcx)			\n\t		vmovaps		%%ymm9 , 0x0a0(%%rcx)		\n\t"\
		"vmovaps	%%ymm5,-0x0e0(%%rdx)			\n\t		vmovaps		%%ymm13,-0x060(%%rdx)		\n\t"\
		"vmovaps	%%ymm2,-0x100(%%rcx)			\n\t		vmovaps		%%ymm10,-0x080(%%rcx)		\n\t"\
		"vmovaps	%%ymm7,-0x100(%%rdx)			\n\t		vmovaps		%%ymm15,-0x080(%%rdx)		\n\t"\
		"vmovaps	%%ymm3,-0x0e0(%%rcx)			\n\t		vmovaps		%%ymm11,-0x060(%%rcx)		\n\t"\
		"vmovaps	%%ymm6, 0x020(%%rdx)			\n\t		vmovaps		%%ymm14, 0x0a0(%%rdx)		\n\t"\
		"\n\t"\
	/***************************************************************************************/\
	/* Now do eight more radix-4 transforms, including the internal and external twiddles: */\
	/***************************************************************************************/\
		"\n\t"\
		/* Using the upper block(s) of the main array for temp-storage in the section below led to a nasty AVX bug to track down: */\
		/* In fermat-mod mode the 4 block addresses in ascending order are add0,1,2,3 with no 'gaps' between blocks, whereas for */\
		/* mersenne-mod the addresses in asc. order are add0,2,3,1 with a gap between contiguous-data-block pairs 0,2 and 3,1. Thus */\
		/* for fermat-mod we need [add2] as the base-address of the 'high-half' block for temp-storage; for mersenne-mod we need [add3]. */\
		/* In both cases we have that add2 < add3 so instead use (add2 - add1): > 0 for fermat-mod, < 0 for mersenne - to differentiate: */\
		"movq	%[__add2],%%rsi		\n\t"/* destroyable copy of add2 */\
		"movq	%[__add2],%%rbx		\n\t"\
		"subq	%[__add1],%%rsi		\n\t"/* rsi = (add2 - add1); if this yields a borrow (i.e. sets CF) it's mersenne, else fermat. */\
		"cmovcq %[__add3],%%rbx	\n\t" /* if CF set (i.e. h > l), copy source [add3] into dest (rbx), else leave dest = [add2]. */\
		"\n\t"\
	/************************************************/				/************************************************/\
	/*...Block 2: t02,t12,t22,t32 -> r10,14,12,16:  */				/*...Block 6: t0A,t1A,t2A,t3A -> r18,1C,1A,1E:  */\
	/************************************************/				/************************************************/\
		"movq		%[__isrt2],%%rsi					\n\t"\
		"movq		%[__r10],%%rax	/* base-addr in rcol = c05/r18, so rax/r10 offset +0x100 vs lcol */\n\t"\
		"movq		%%rsi,%%rcx					\n\t"\
		"movq		%%rsi,%%rdx					\n\t"\
		"movq		%[__c01],%%r10					\n\t"\
		"addq		$0x020,%%rcx	/* cc0 */		\n\t"\
		"addq		$0x060,%%rdx	/* cc1 */		\n\t"\
		"vmovaps		0x040(%%rax),%%ymm4					\n\t		vmovaps		 0x140(%%rax),%%ymm12			\n\t"\
		"vmovaps		0x0c0(%%rax),%%ymm0					\n\t		vmovaps		 0x1c0(%%rax),%%ymm8 			\n\t"\
		"vmovaps		0x060(%%rax),%%ymm5					\n\t		vmovaps		 0x160(%%rax),%%ymm15			\n\t"\
		"vmovaps		0x0e0(%%rax),%%ymm3					\n\t"	/*	vmovaps		 0x1e0(%%rax),%%ymm9 			\n\t"*/\
		"vmovaps			%%ymm4 ,%%ymm6					\n\t		vmovaps		 	%%ymm12,%%ymm14				\n\t"\
		"vmovaps			%%ymm0 ,%%ymm2					\n\t"	/*	vmovaps		 	%%ymm8 ,%%ymm10				\n\t"*/\
		"vmovaps			%%ymm5 ,%%ymm7					\n\t"	/*	vmovaps		 	%%ymm13,%%ymm15				\n\t"*/\
	/*	"vmovaps			%%ymm1 ,%%ymm3			*/			"		vmovaps		 	%%ymm9 ,%%ymm11				\n\t"\
	/* Strategic register-usage for sincos consts here cuts loads from 16 to 4: */\
		"vmovaps			 (%%rdx),%%ymm9 			\n\t		vmovaps		0x020(%%rdx),%%ymm10			\n\t"\
		"vmovaps		0x040(%%rdx),%%ymm1 			\n\t		vmovaps		0x060(%%rdx),%%ymm13			\n\t"\
		"vmulpd			%%ymm10,%%ymm7,%%ymm7				\n\t		vmulpd			 %%ymm1 ,%%ymm15,%%ymm15	\n\t"\
		"vmulpd			%%ymm13,%%ymm3,%%ymm3				\n\t		vmulpd			 %%ymm10,%%ymm11,%%ymm11	\n\t"\
		"vmulpd			%%ymm10,%%ymm6,%%ymm6				\n\t		vmulpd			 %%ymm1 ,%%ymm14,%%ymm14	\n\t"\
		"vmulpd			%%ymm13,%%ymm2,%%ymm2				\n\t		vmulpd		0x1c0(%%rax),%%ymm10,%%ymm10	\n\t"\
	"vfmadd132pd		%%ymm9 ,%%ymm7,%%ymm4				\n\t	vfmadd132pd	 		 %%ymm13,%%ymm15,%%ymm12	\n\t"\
	"vfmadd132pd		%%ymm1 ,%%ymm3,%%ymm0				\n\t	vfmsub132pd	 		 %%ymm9 ,%%ymm11,%%ymm8 	\n\t"\
	"vfmsub132pd		%%ymm9 ,%%ymm6,%%ymm5				\n\t	vfmsub132pd	 	0x160(%%rax),%%ymm14,%%ymm13	\n\t"\
	"vfmsub132pd	0xe0(%%rax),%%ymm2,%%ymm1				\n\t	vfmadd132pd	 	0x1e0(%%rax),%%ymm10,%%ymm9 	\n\t"\
		"vmovaps		%%ymm5,%%ymm7					\n\t		vmovaps		%%ymm13,%%ymm15					\n\t"\
		"vmovaps		%%ymm4,%%ymm6					\n\t		vmovaps		%%ymm12,%%ymm14					\n\t"\
		"vsubpd		%%ymm0,%%ymm4,%%ymm4			\n\t		vaddpd		%%ymm8 ,%%ymm12,%%ymm12			\n\t"\
		"vsubpd		%%ymm1,%%ymm5,%%ymm5			\n\t		vaddpd		%%ymm9 ,%%ymm13,%%ymm13			\n\t"\
		"vaddpd		%%ymm0,%%ymm6,%%ymm6			\n\t		vsubpd		%%ymm8 ,%%ymm14,%%ymm14			\n\t"\
		"vaddpd		%%ymm1,%%ymm7,%%ymm7			\n\t		vsubpd		%%ymm9 ,%%ymm15,%%ymm15			\n\t"\
		"vmovaps		 0x080(%%rax),%%ymm0				\n\t		vmovaps		 0x180(%%rax),%%ymm8 				\n\t"\
		"vmovaps		 0x0a0(%%rax),%%ymm1				\n\t		vmovaps		 0x1a0(%%rax),%%ymm9 				\n\t"\
		"vmovaps		 0x080(%%rax),%%ymm2				\n\t		vmovaps		 0x180(%%rax),%%ymm10				\n\t"\
		"vmovaps		 0x0a0(%%rax),%%ymm3				\n\t		vmovaps		 0x1a0(%%rax),%%ymm11				\n\t"\
		"vmulpd		 0x020(%%rcx),%%ymm0,%%ymm0			\n\t		vmulpd			  (%%rcx),%%ymm8 ,%%ymm8 			\n\t"\
		"vmulpd		 0x020(%%rcx),%%ymm1,%%ymm1			\n\t		vmulpd			  (%%rcx),%%ymm9 ,%%ymm9 			\n\t"\
	"vfmsub132pd		  (%%rcx),%%ymm0,%%ymm3			\n\t	vfmadd132pd		 0x020(%%rcx),%%ymm8 ,%%ymm11			\n\t"\
	"vfmadd132pd		  (%%rcx),%%ymm1,%%ymm2			\n\t	vfmsub132pd		 0x020(%%rcx),%%ymm9 ,%%ymm10			\n\t"\
		"vmovaps			  (%%rax),%%ymm0				\n\t		vmovaps		 0x100(%%rax),%%ymm8 				\n\t"\
		"vmovaps		 0x020(%%rax),%%ymm1				\n\t		vmovaps		 0x120(%%rax),%%ymm9 				\n\t"\
		"vsubpd		%%ymm2,%%ymm0,%%ymm0			\n\t		vsubpd		%%ymm10,%%ymm8 ,%%ymm8 			\n\t"\
		"vsubpd		%%ymm3,%%ymm1,%%ymm1			\n\t		vsubpd		%%ymm11,%%ymm9 ,%%ymm9 			\n\t"\
	"vfmadd132pd		(%%rdi),%%ymm0,%%ymm2			\n\t	vfmadd132pd		(%%rdi),%%ymm8 ,%%ymm10			\n\t"\
	"vfmadd132pd		(%%rdi),%%ymm1,%%ymm3			\n\t	vfmadd132pd		(%%rdi),%%ymm9 ,%%ymm11			\n\t"\
		"vsubpd		%%ymm6,%%ymm2,%%ymm2			\n\t		vsubpd		%%ymm14,%%ymm8 ,%%ymm8 			\n\t"\
		"vsubpd		%%ymm7,%%ymm3,%%ymm3			\n\t		vsubpd		%%ymm15,%%ymm9 ,%%ymm9 			\n\t"\
	"vfmadd132pd		(%%rdi),%%ymm2,%%ymm6			\n\t	vfmadd132pd		(%%rdi),%%ymm8 ,%%ymm14			\n\t"\
	"vfmadd132pd		(%%rdi),%%ymm3,%%ymm7			\n\t	vfmadd132pd		(%%rdi),%%ymm9 ,%%ymm15			\n\t"\
		"vmovaps		%%ymm2, 0x040(%%rax)			\n\t		vmovaps		%%ymm8 , 0x140(%%rax)			\n\t"\
		"vmovaps		%%ymm3, 0x060(%%rax)			\n\t		vmovaps		%%ymm9 , 0x160(%%rax)			\n\t"\
		"vmovaps		%%ymm6,%%ymm2					\n\t		vmovaps		%%ymm14,%%ymm8 					\n\t"\
		"vmovaps		%%ymm7,%%ymm3					\n\t		vmovaps		%%ymm15,%%ymm9 					\n\t"\
		"vmulpd			  (%%r10),%%ymm6,%%ymm6			\n\t		vmulpd		 0x100(%%r10),%%ymm14,%%ymm14			\n\t"\
		"vmulpd			  (%%r10),%%ymm7,%%ymm7			\n\t		vmulpd		 0x100(%%r10),%%ymm15,%%ymm15			\n\t"\
	" vfmadd231pd	 0x020(%%r10),%%ymm3,%%ymm6			\n\t	 vfmadd231pd	 0x120(%%r10),%%ymm9 ,%%ymm14			\n\t"\
	"vfnmadd231pd	 0x020(%%r10),%%ymm2,%%ymm7			\n\t	vfnmadd231pd	 0x120(%%r10),%%ymm8 ,%%ymm15			\n\t"\
		"vmovaps		%%ymm7, 0x020(%%rbx)			\n\t		vmovaps		%%ymm15, 0x0a0(%%rbx)			\n\t"\
		"vmovaps		%%ymm6,      (%%rbx)			\n\t		vmovaps		%%ymm14, 0x080(%%rbx)			\n\t"\
		"vmovaps		 0x040(%%rax),%%ymm6				\n\t		vmovaps		 0x140(%%rax),%%ymm14				\n\t"\
		"vmovaps		 0x060(%%rax),%%ymm7				\n\t		vmovaps		 0x160(%%rax),%%ymm15				\n\t"\
		"vmovaps		%%ymm6,%%ymm2					\n\t		vmovaps		%%ymm14,%%ymm8 					\n\t"\
		"vmovaps		%%ymm7,%%ymm3					\n\t		vmovaps		%%ymm15,%%ymm9 					\n\t"\
		"vmulpd		 0x040(%%r10),%%ymm6,%%ymm6			\n\t		vmulpd		 0x140(%%r10),%%ymm14,%%ymm14			\n\t"\
		"vmulpd		 0x040(%%r10),%%ymm7,%%ymm7			\n\t		vmulpd		 0x140(%%r10),%%ymm15,%%ymm15			\n\t"\
		"vsubpd		%%ymm5,%%ymm0,%%ymm0			\n\t		vsubpd		%%ymm13,%%ymm10,%%ymm10			\n\t"\
		"vsubpd		%%ymm4,%%ymm1,%%ymm1			\n\t		vsubpd		%%ymm12,%%ymm11,%%ymm11			\n\t"\
	" vfmadd231pd	 0x060(%%r10),%%ymm3,%%ymm6			\n\t	 vfmadd231pd	 0x160(%%r10),%%ymm9 ,%%ymm14			\n\t"\
	"vfnmadd231pd	 0x060(%%r10),%%ymm2,%%ymm7			\n\t	vfnmadd231pd	 0x160(%%r10),%%ymm8 ,%%ymm15			\n\t"\
	"vfmadd132pd		  (%%rdi),%%ymm0,%%ymm5			\n\t	vfmadd132pd		(%%rdi),%%ymm10,%%ymm13			\n\t"\
	"vfmadd132pd		  (%%rdi),%%ymm1,%%ymm4			\n\t	vfmadd132pd		(%%rdi),%%ymm11,%%ymm12			\n\t"\
		"vmovaps		%%ymm7, 0x220(%%rbx)			\n\t		vmovaps		%%ymm15, 0x2a0(%%rbx)			\n\t"\
		"vmovaps		%%ymm6, 0x200(%%rbx)			\n\t		vmovaps		%%ymm14, 0x280(%%rbx)			\n\t"\
		"addq		$0x080,%%r10					\n\t"\
		"vmovaps		%%ymm5,%%ymm2					\n\t		vmovaps		%%ymm13,%%ymm8 					\n\t"\
		"vmovaps		%%ymm1,%%ymm3					\n\t		vmovaps		%%ymm11,%%ymm9 					\n\t"\
		"vmovaps		%%ymm0,%%ymm6					\n\t		vmovaps		%%ymm10,%%ymm14					\n\t"\
		"vmovaps		%%ymm4,%%ymm7					\n\t		vmovaps		%%ymm12,%%ymm15					\n\t"\
		"vmulpd			  (%%r10),%%ymm5,%%ymm5			\n\t		vmulpd		 0x100(%%r10),%%ymm13,%%ymm13			\n\t"\
		"vmulpd			  (%%r10),%%ymm1,%%ymm1			\n\t		vmulpd		 0x100(%%r10),%%ymm11,%%ymm11			\n\t"\
		"vmulpd		 0x040(%%r10),%%ymm0,%%ymm0			\n\t		vmulpd		 0x140(%%r10),%%ymm10,%%ymm10			\n\t"\
		"vmulpd		 0x040(%%r10),%%ymm4,%%ymm4			\n\t		vmulpd		 0x140(%%r10),%%ymm12,%%ymm12			\n\t"\
	" vfmadd231pd	 0x020(%%r10),%%ymm3,%%ymm5			\n\t	 vfmadd231pd	 0x120(%%r10),%%ymm9 ,%%ymm13			\n\t"\
	"vfnmadd231pd	 0x020(%%r10),%%ymm2,%%ymm1			\n\t	vfnmadd231pd	 0x120(%%r10),%%ymm8 ,%%ymm11			\n\t"\
	" vfmadd231pd	 0x060(%%r10),%%ymm7,%%ymm0			\n\t	 vfmadd231pd	 0x160(%%r10),%%ymm15,%%ymm10			\n\t"\
	"vfnmadd231pd	 0x060(%%r10),%%ymm6,%%ymm4			\n\t	vfnmadd231pd	 0x160(%%r10),%%ymm14,%%ymm12			\n\t"\
		"vmovaps		%%ymm1, 0x120(%%rbx)			\n\t		vmovaps		%%ymm11, 0x1a0(%%rbx)			\n\t"\
		"vmovaps		%%ymm5, 0x100(%%rbx)			\n\t		vmovaps		%%ymm13, 0x180(%%rbx)			\n\t"\
		"vmovaps		%%ymm4, 0x320(%%rbx)			\n\t		vmovaps		%%ymm12, 0x3a0(%%rbx)			\n\t"\
		"vmovaps		%%ymm0, 0x300(%%rbx)			\n\t		vmovaps		%%ymm10, 0x380(%%rbx)			\n\t"\
		"\n\t"\
	/************************************************/				/************************************************/\
	/*...Block 4: t06,t16,t26,t36 -> r30,34,32,36:  */				/*...Block 8: t0E,t1E,t2E,t3E -> r38,3C,3A,3E:  */\
	/************************************************/				/************************************************/\
		"addq		$0x400,%%rax						\n\t		addq		$0x100,%%r10						\n\t"\
		"vmovaps		0x040(%%rax),%%ymm4					\n\t		vmovaps		 0x140(%%rax),%%ymm12			\n\t"\
		"vmovaps		0x0c0(%%rax),%%ymm0					\n\t		vmovaps		 0x1c0(%%rax),%%ymm10			\n\t"\
		"vmovaps		0x060(%%rax),%%ymm7					\n\t		vmovaps		 0x160(%%rax),%%ymm15			\n\t"\
		"vmovaps		0x0e0(%%rax),%%ymm1					\n\t		vmovaps		 0x1e0(%%rax),%%ymm9 			\n\t"\
		"vmovaps			%%ymm4 ,%%ymm6					\n\t"	/*	vmovaps		 	%%ymm12,%%ymm14				\n\t"*/\
		"vmovaps			%%ymm0 ,%%ymm2					\n\t"	/*	vmovaps		 	%%ymm8 ,%%ymm10				\n\t"*/\
	/*	"vmovaps			%%ymm5 ,%%ymm7	*/						/*	vmovaps		 	%%ymm13,%%ymm15				\n\t"*/\
		"vmovaps			%%ymm1 ,%%ymm3					\n\t		vmovaps		 	%%ymm9 ,%%ymm11				\n\t"\
	/* Strategic register-usage for sincos consts here cuts loads from 16 to 4: */\
		"vmovaps			 (%%rdx),%%ymm14			\n\t		vmovaps		0x020(%%rdx),%%ymm13			\n\t"\
		"vmovaps		0x040(%%rdx),%%ymm5 			\n\t		vmovaps		0x060(%%rdx),%%ymm8 			\n\t"\
	/* We cyc-permuted the 4 paired-MUL lines so as to put a '     (%%rdx)' in the SRC3 slot of the last line's rcol: */\
		"vmulpd		 	 %%ymm14,%%ymm2,%%ymm2				\n\t		vmulpd			 %%ymm5 ,%%ymm10,%%ymm10	\n\t"\
		"vmulpd		 	 %%ymm8 ,%%ymm7,%%ymm7				\n\t		vmulpd			 %%ymm14,%%ymm15,%%ymm15	\n\t"\
		"vmulpd		 	 %%ymm14,%%ymm3,%%ymm3				\n\t		vmulpd			 %%ymm5 ,%%ymm11,%%ymm11	\n\t"\
		"vmulpd		 	 %%ymm8 ,%%ymm6,%%ymm6				\n\t		vmulpd		0x140(%%rax),%%ymm14,%%ymm14	\n\t"\
	/* ...and similarly cyc-permute these 4 lines so middle (SRC2) op's reg-index order matches that of the above 4: */\
	"vfmadd132pd	 	 %%ymm13,%%ymm2,%%ymm1				\n\t	vfmsub132pd	 		 %%ymm8 ,%%ymm10,%%ymm9 	\n\t"\
	"vfmadd132pd	 	 %%ymm5 ,%%ymm7,%%ymm4				\n\t	vfmadd132pd	 		 %%ymm13,%%ymm15,%%ymm12	\n\t"\
	"vfmsub132pd	 	 %%ymm13,%%ymm3,%%ymm0				\n\t	vfmadd132pd	 	0x1c0(%%rax),%%ymm11,%%ymm8 	\n\t"\
	"vfmsub132pd	 0x60(%%rax),%%ymm6,%%ymm5				\n\t	vfmsub132pd	 	0x160(%%rax),%%ymm14,%%ymm13	\n\t"\
		"vmovaps		%%ymm5,%%ymm7					\n\t		vmovaps		%%ymm13,%%ymm15					\n\t"\
		"vmovaps		%%ymm4,%%ymm6					\n\t		vmovaps		%%ymm12,%%ymm14					\n\t"\
		"vaddpd		%%ymm0,%%ymm4,%%ymm4			\n\t		vaddpd		%%ymm8 ,%%ymm12,%%ymm12			\n\t"\
		"vaddpd		%%ymm1,%%ymm5,%%ymm5			\n\t		vaddpd		%%ymm9 ,%%ymm13,%%ymm13			\n\t"\
		"vsubpd		%%ymm0,%%ymm6,%%ymm6			\n\t		vsubpd		%%ymm8 ,%%ymm14,%%ymm14			\n\t"\
		"vsubpd		%%ymm1,%%ymm7,%%ymm7			\n\t		vsubpd		%%ymm9 ,%%ymm15,%%ymm15			\n\t"\
		"vmovaps		 0x080(%%rax),%%ymm0				\n\t		vmovaps		 0x180(%%rax),%%ymm8 				\n\t"\
		"vmovaps		 0x0a0(%%rax),%%ymm1				\n\t		vmovaps		 0x1a0(%%rax),%%ymm9 				\n\t"\
		"vmovaps		 0x080(%%rax),%%ymm2				\n\t		vmovaps		 0x180(%%rax),%%ymm10				\n\t"\
		"vmovaps		 0x0a0(%%rax),%%ymm3				\n\t		vmovaps		 0x1a0(%%rax),%%ymm11				\n\t"\
		"vmulpd			  (%%rcx),%%ymm0,%%ymm0			\n\t		vmulpd		 0x020(%%rcx),%%ymm8 ,%%ymm8 			\n\t"\
		"vmulpd			  (%%rcx),%%ymm1,%%ymm1			\n\t		vmulpd		 0x020(%%rcx),%%ymm9 ,%%ymm9 			\n\t"\
	"vfmsub132pd	 0x020(%%rcx),%%ymm0,%%ymm3			\n\t	vfmadd132pd			  (%%rcx),%%ymm8 ,%%ymm11			\n\t"\
	"vfmadd132pd	 0x020(%%rcx),%%ymm1,%%ymm2			\n\t	vfmsub132pd			  (%%rcx),%%ymm9 ,%%ymm10			\n\t"\
		"vmovaps			  (%%rax),%%ymm0				\n\t		vmovaps		 0x100(%%rax),%%ymm8 				\n\t"\
		"vmovaps		 0x020(%%rax),%%ymm1				\n\t		vmovaps		 0x120(%%rax),%%ymm9 				\n\t"\
		"vsubpd		%%ymm2,%%ymm0,%%ymm0			\n\t		vsubpd		%%ymm10,%%ymm8 ,%%ymm8 			\n\t"\
		"vsubpd		%%ymm3,%%ymm1,%%ymm1			\n\t		vsubpd		%%ymm11,%%ymm9 ,%%ymm9 			\n\t"\
	"vfmadd132pd		(%%rdi),%%ymm0,%%ymm2			\n\t	vfmadd132pd		(%%rdi),%%ymm8 ,%%ymm10			\n\t"\
	"vfmadd132pd		(%%rdi),%%ymm1,%%ymm3			\n\t	vfmadd132pd		(%%rdi),%%ymm9 ,%%ymm11			\n\t"\
		"vsubpd		%%ymm6,%%ymm2,%%ymm2			\n\t		vsubpd		%%ymm14,%%ymm8 ,%%ymm8 			\n\t"\
		"vsubpd		%%ymm7,%%ymm3,%%ymm3			\n\t		vsubpd		%%ymm15,%%ymm9 ,%%ymm9 			\n\t"\
	"vfmadd132pd		(%%rdi),%%ymm2,%%ymm6			\n\t	vfmadd132pd		(%%rdi),%%ymm8 ,%%ymm14			\n\t"\
	"vfmadd132pd		(%%rdi),%%ymm3,%%ymm7			\n\t	vfmadd132pd		(%%rdi),%%ymm9 ,%%ymm15			\n\t"\
		"addq		$0x080,%%r10						\n\t"\
		"vmovaps		%%ymm2, 0x040(%%rax)			\n\t		vmovaps		%%ymm8 , 0x140(%%rax)			\n\t"\
		"vmovaps		%%ymm3, 0x060(%%rax)			\n\t		vmovaps		%%ymm9 , 0x160(%%rax)			\n\t"\
		"vmovaps		%%ymm6,%%ymm2					\n\t		vmovaps		%%ymm14,%%ymm8 					\n\t"\
		"vmovaps		%%ymm7,%%ymm3					\n\t		vmovaps		%%ymm15,%%ymm9 					\n\t"\
		"vmulpd			  (%%r10),%%ymm6,%%ymm6			\n\t		vmulpd		 0x100(%%r10),%%ymm14,%%ymm14			\n\t"\
		"vmulpd			  (%%r10),%%ymm7,%%ymm7			\n\t		vmulpd		 0x100(%%r10),%%ymm15,%%ymm15			\n\t"\
	" vfmadd231pd	 0x020(%%r10),%%ymm3,%%ymm6			\n\t	 vfmadd231pd	 0x120(%%r10),%%ymm9 ,%%ymm14			\n\t"\
	"vfnmadd231pd	 0x020(%%r10),%%ymm2,%%ymm7			\n\t	vfnmadd231pd	 0x120(%%r10),%%ymm8 ,%%ymm15			\n\t"\
		"vmovaps		%%ymm7, 0x060(%%rbx)			\n\t		vmovaps		%%ymm15, 0x0e0(%%rbx)			\n\t"\
		"vmovaps		%%ymm6, 0x040(%%rbx)			\n\t		vmovaps		%%ymm14, 0x0c0(%%rbx)			\n\t"\
		"vmovaps		 0x040(%%rax),%%ymm6				\n\t		vmovaps		 0x140(%%rax),%%ymm14				\n\t"\
		"vmovaps		 0x060(%%rax),%%ymm7				\n\t		vmovaps		 0x160(%%rax),%%ymm15				\n\t"\
		"vmovaps		%%ymm6,%%ymm2					\n\t		vmovaps		%%ymm14,%%ymm8 					\n\t"\
		"vmovaps		%%ymm7,%%ymm3					\n\t		vmovaps		%%ymm15,%%ymm9 					\n\t"\
		"vmulpd		 0x040(%%r10),%%ymm6,%%ymm6			\n\t		vmulpd		 0x140(%%r10),%%ymm14,%%ymm14			\n\t"\
		"vmulpd		 0x040(%%r10),%%ymm7,%%ymm7			\n\t		vmulpd		 0x140(%%r10),%%ymm15,%%ymm15			\n\t"\
		"vsubpd		%%ymm5,%%ymm0,%%ymm0			\n\t		vsubpd		%%ymm13,%%ymm10,%%ymm10			\n\t"\
		"vsubpd		%%ymm4,%%ymm1,%%ymm1			\n\t		vsubpd		%%ymm12,%%ymm11,%%ymm11			\n\t"\
	" vfmadd231pd	 0x060(%%r10),%%ymm3,%%ymm6			\n\t	 vfmadd231pd	 0x160(%%r10),%%ymm9 ,%%ymm14			\n\t"\
	"vfnmadd231pd	 0x060(%%r10),%%ymm2,%%ymm7			\n\t	vfnmadd231pd	 0x160(%%r10),%%ymm8 ,%%ymm15			\n\t"\
	"vfmadd132pd		  (%%rdi),%%ymm0,%%ymm5			\n\t	vfmadd132pd		(%%rdi),%%ymm10,%%ymm13			\n\t"\
	"vfmadd132pd		  (%%rdi),%%ymm1,%%ymm4			\n\t	vfmadd132pd		(%%rdi),%%ymm11,%%ymm12			\n\t"\
		"vmovaps		%%ymm7, 0x260(%%rbx)			\n\t		vmovaps		%%ymm15, 0x2e0(%%rbx)			\n\t"\
		"vmovaps		%%ymm6, 0x240(%%rbx)			\n\t		vmovaps		%%ymm14, 0x2c0(%%rbx)			\n\t"\
		"addq		$0x080,%%r10					\n\t"\
		"vmovaps		%%ymm5,%%ymm2					\n\t		vmovaps		%%ymm13,%%ymm8 					\n\t"\
		"vmovaps		%%ymm1,%%ymm3					\n\t		vmovaps		%%ymm11,%%ymm9 					\n\t"\
		"vmovaps		%%ymm0,%%ymm6					\n\t		vmovaps		%%ymm10,%%ymm14					\n\t"\
		"vmovaps		%%ymm4,%%ymm7					\n\t		vmovaps		%%ymm12,%%ymm15					\n\t"\
		"vmulpd			  (%%r10),%%ymm5,%%ymm5			\n\t		vmulpd		 0x100(%%r10),%%ymm13,%%ymm13			\n\t"\
		"vmulpd			  (%%r10),%%ymm1,%%ymm1			\n\t		vmulpd		 0x100(%%r10),%%ymm11,%%ymm11			\n\t"\
		"vmulpd		 0x040(%%r10),%%ymm0,%%ymm0			\n\t		vmulpd		 0x140(%%r10),%%ymm10,%%ymm10			\n\t"\
		"vmulpd		 0x040(%%r10),%%ymm4,%%ymm4			\n\t		vmulpd		 0x140(%%r10),%%ymm12,%%ymm12			\n\t"\
	" vfmadd231pd	 0x020(%%r10),%%ymm3,%%ymm5			\n\t	 vfmadd231pd	 0x120(%%r10),%%ymm9 ,%%ymm13			\n\t"\
	"vfnmadd231pd	 0x020(%%r10),%%ymm2,%%ymm1			\n\t	vfnmadd231pd	 0x120(%%r10),%%ymm8 ,%%ymm11			\n\t"\
	" vfmadd231pd	 0x060(%%r10),%%ymm7,%%ymm0			\n\t	 vfmadd231pd	 0x160(%%r10),%%ymm15,%%ymm10			\n\t"\
	"vfnmadd231pd	 0x060(%%r10),%%ymm6,%%ymm4			\n\t	vfnmadd231pd	 0x160(%%r10),%%ymm14,%%ymm12			\n\t"\
		"vmovaps		%%ymm1, 0x160(%%rbx)			\n\t		vmovaps		%%ymm11, 0x1e0(%%rbx)			\n\t"\
		"vmovaps		%%ymm5, 0x140(%%rbx)			\n\t		vmovaps		%%ymm13, 0x1c0(%%rbx)			\n\t"\
		"vmovaps		%%ymm4, 0x360(%%rbx)			\n\t		vmovaps		%%ymm12, 0x3e0(%%rbx)			\n\t"\
		"vmovaps		%%ymm0, 0x340(%%rbx)			\n\t		vmovaps		%%ymm10, 0x3c0(%%rbx)			\n\t"\
		"\n\t"\
	/************************************************/				/************************************************/\
	/*...Block 1: t00,t10,t20,t30 -> r00,04,02,06:  */				/*...Block 5: t08,t18,t28,t38 -> r08,0C,0A,0E:  */\
	/************************************************/				/************************************************/\
		"movq	%[__r00],%%rdx	/* base-addr in rcol = r08, so rdx+0x100 in rcol */	\n\t	vmovaps	(%%rsi),%%ymm10	\n\t"/* isrt2 */\
		"vmovaps			  (%%rdx),%%ymm0				\n\t		vmovaps		 0x140(%%rdx),%%ymm12				\n\t"\
		"vmovaps		 0x020(%%rdx),%%ymm1				\n\t		vmovaps		 0x160(%%rdx),%%ymm13				\n\t"\
		"vmovaps		 0x080(%%rdx),%%ymm2				\n\t		vmovaps		 0x1c0(%%rdx),%%ymm8 				\n\t"\
		"vmovaps		 0x0a0(%%rdx),%%ymm3				\n\t		vmovaps		 0x1e0(%%rdx),%%ymm9 				\n\t"\
		"vsubpd		 0x080(%%rdx),%%ymm0,%%ymm0			\n\t		vaddpd		 0x160(%%rdx),%%ymm12,%%ymm12			\n\t"\
		"vsubpd		 0x0a0(%%rdx),%%ymm1,%%ymm1			\n\t		vsubpd		 0x140(%%rdx),%%ymm13,%%ymm13			\n\t"\
		"vaddpd			  (%%rdx),%%ymm2,%%ymm2			\n\t		vsubpd		 0x1e0(%%rdx),%%ymm8 ,%%ymm8 			\n\t"\
		"vaddpd		 0x020(%%rdx),%%ymm3,%%ymm3			\n\t		vaddpd		 0x1c0(%%rdx),%%ymm9 ,%%ymm9 			\n\t"\
		"vmovaps		 0x040(%%rdx),%%ymm4				\n\t		vmulpd			%%ymm10,%%ymm8 ,%%ymm8 			\n\t"\
		"vmovaps		 0x060(%%rdx),%%ymm5				\n\t		vmulpd			%%ymm10,%%ymm9 ,%%ymm9 			\n\t"\
		"vmovaps		 0x0c0(%%rdx),%%ymm6				\n\t"\
		"vmovaps		 0x0e0(%%rdx),%%ymm7				\n\t"\
		"vsubpd		 0x0c0(%%rdx),%%ymm4,%%ymm4			\n\t		vmovaps		%%ymm12,%%ymm14					\n\t"\
		"vsubpd		 0x0e0(%%rdx),%%ymm5,%%ymm5			\n\t		vmovaps		%%ymm13,%%ymm15					\n\t"\
		"vaddpd		 0x040(%%rdx),%%ymm6,%%ymm6			\n\t	vfmsub132pd		%%ymm10,%%ymm8 ,%%ymm12			\n\t"\
		"vaddpd		 0x060(%%rdx),%%ymm7,%%ymm7			\n\t	vfmsub132pd		%%ymm10,%%ymm9 ,%%ymm13			\n\t"\
		/* base-twiddle in l/rcol = c00/c04, so rcx+0x100 in rcol*/"vfmadd132pd		%%ymm10,%%ymm8 ,%%ymm14			\n\t"\
		"movq		%[__c10],%%rcx						\n\t	vfmadd132pd		%%ymm10,%%ymm9 ,%%ymm15			\n\t"\
		"vaddpd			%%ymm6,%%ymm2,%%ymm2			\n\t		vmovaps		 0x100(%%rdx),%%ymm8 				\n\t"\
		"vaddpd			%%ymm7,%%ymm3,%%ymm3			\n\t		vmovaps		 0x120(%%rdx),%%ymm9 				\n\t"\
		"vmovaps		%%ymm2,      (%%rdx)			\n\t		vmovaps		 0x180(%%rdx),%%ymm10				\n\t"\
		"vmovaps		%%ymm3, 0x020(%%rdx)			\n\t		vmovaps		 0x1a0(%%rdx),%%ymm11				\n\t"\
		"vaddpd			%%ymm6,%%ymm6,%%ymm6			\n\t		vsubpd		 0x1a0(%%rdx),%%ymm8 ,%%ymm8 			\n\t"\
		"vaddpd			%%ymm7,%%ymm7,%%ymm7			\n\t		vsubpd		 0x180(%%rdx),%%ymm9 ,%%ymm9 			\n\t"\
		"vsubpd			%%ymm6,%%ymm2,%%ymm2			\n\t		vaddpd		 0x100(%%rdx),%%ymm11,%%ymm11			\n\t"\
		"vsubpd			%%ymm7,%%ymm3,%%ymm3			\n\t		vaddpd		 0x120(%%rdx),%%ymm10,%%ymm10			\n\t"\
		"vmovaps		%%ymm2,%%ymm6					\n\t		vsubpd			%%ymm12,%%ymm11,%%ymm11			\n\t"\
		"vmovaps		%%ymm3,%%ymm7					\n\t		vsubpd			%%ymm13,%%ymm9 ,%%ymm9 			\n\t"\
		"vmulpd			  (%%rcx),%%ymm2,%%ymm2	\n\t"/* c10 */"	vfmadd132pd		(%%rdi),%%ymm11,%%ymm12			\n\t"\
		"vmulpd			  (%%rcx),%%ymm3,%%ymm3			\n\t	vfmadd132pd		(%%rdi),%%ymm9 ,%%ymm13			\n\t"\
	" vfmadd231pd	 0x020(%%rcx),%%ymm7,%%ymm2			\n\t		vmovaps		%%ymm11, 0x140(%%rdx)			\n\t"\
	"vfnmadd231pd	 0x020(%%rcx),%%ymm6,%%ymm3			\n\t		vmovaps		%%ymm9 , 0x160(%%rdx)			\n\t"\
		"subq $0x40,%%rcx	/* put c00 in rcx to ease bookkeeping*/\n\t	vmovaps		%%ymm12,%%ymm11					\n\t"\
		/* add0,1 in rax,rbx; __r00 in rdx: */						"	vmovaps		%%ymm13,%%ymm9 					\n\t"\
		/* For each complex output octet, complex pairs having */	"	vmulpd		 0x100(%%rcx),%%ymm12,%%ymm12	/* c04 */\n\t"\
		/* reads from offsets 0x0..,0x1..,0x2..,0x3.. go into  */	"	vmulpd		 0x100(%%rcx),%%ymm13,%%ymm13			\n\t"\
		/* local-mem pairs rXY + 00/10, 04/14, 02/12, 06/16.   */	" vfmadd231pd	 0x120(%%rcx),%%ymm9 ,%%ymm12			\n\t"\
		/* For 1st octet we read from offsets [0x2..,0x0..],   */	"vfnmadd231pd	 0x120(%%rcx),%%ymm11,%%ymm13			\n\t"\
		/* [0x1..,0x3], other 3 octets use order [0,2],[1,3].  */\
		"vmovaps	0x220(%%rbx),%%ymm7						\n\t		vmovaps	0x0a0(%%rbx),%%ymm11						\n\t"\
		"vmovaps	0x200(%%rbx),%%ymm6						\n\t		vmovaps	0x080(%%rbx),%%ymm9 						\n\t"\
		"vmovaps	%%ymm3,0x060(%%rdx)						\n\t		vmovaps	%%ymm13,0x120(%%rdx)						\n\t"\
		"vmovaps	%%ymm2,0x040(%%rdx)		\n\t/* r02,03 */			vmovaps	%%ymm12,0x100(%%rdx)			\n\t"/* r08,09 */\
		"vmovaps	%%ymm7,0x260(%%rdx)						\n\t		vmovaps	%%ymm11,0x320(%%rdx)						\n\t"\
		"vmovaps	%%ymm6,0x240(%%rdx)		\n\t/* r12,13 */			vmovaps	%%ymm9 ,0x300(%%rdx)			\n\t"/* r18,19 */\
		"vmovaps		 0x020(%%rdx),%%ymm3				\n\t		vmovaps		0x140(%%rdx),%%ymm12				\n\t"\
		"vmovaps			  (%%rdx),%%ymm2				\n\t		vmovaps		0x160(%%rdx),%%ymm13				\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm7	\n\t vmovaps %%ymm12,%%ymm11 \n\t	vmulpd	 0x140(%%rcx),%%ymm12,%%ymm12	/* c14 */\n\t"\
		"vmovaps	0x000(%%rbx),%%ymm6	\n\t vmovaps %%ymm13,%%ymm9  \n\t	vmulpd	 0x140(%%rcx),%%ymm13,%%ymm13			\n\t"\
		"vmovaps	%%ymm3,0x020(%%rdx)						\n\t	 vfmadd231pd	 0x160(%%rcx),%%ymm9 ,%%ymm12			\n\t"\
		"vmovaps	%%ymm2,     (%%rdx)			/* r00,01 */\n\t	vfnmadd231pd	 0x160(%%rcx),%%ymm11,%%ymm13			\n\t"\
		"vmovaps	%%ymm7,0x220(%%rdx)						\n\t		vmovaps	0x2a0(%%rbx),%%ymm11						\n\t"\
		"vmovaps	%%ymm6,0x200(%%rdx)			/* r10,11 */\n\t		vmovaps	0x280(%%rbx),%%ymm9 						\n\t"\
		"vaddpd			%%ymm5,%%ymm0,%%ymm0			\n\t		vmovaps	%%ymm13,0x160(%%rdx)						\n\t"\
		"vsubpd			%%ymm4,%%ymm1,%%ymm1			\n\t		vmovaps	%%ymm12,0x140(%%rdx)			\n\t"/* r0a,0b */\
		"vmovaps		%%ymm0,%%ymm2					\n\t		vmovaps	%%ymm11,0x360(%%rdx)						\n\t"\
		"vmovaps		%%ymm1,%%ymm3					\n\t		vmovaps	%%ymm9 ,0x340(%%rdx)			\n\t"/* r1a,1b */\
		"vaddpd			%%ymm5,%%ymm5,%%ymm5			\n\t		vsubpd			%%ymm15,%%ymm8 ,%%ymm8 			\n\t"\
		"vaddpd			%%ymm4,%%ymm4,%%ymm4			\n\t		vsubpd			%%ymm14,%%ymm10,%%ymm10			\n\t"\
		"vmovaps		%%ymm0,%%ymm6					\n\t	vfmadd132pd		(%%rdi),%%ymm8 ,%%ymm15			\n\t"\
		"vmovaps		%%ymm1,%%ymm7					\n\t	vfmadd132pd		(%%rdi),%%ymm10,%%ymm14			\n\t"\
		"vmulpd		 0x080(%%rcx),%%ymm2,%%ymm2	/* c08*/\n\t		vmovaps		%%ymm15,%%ymm12					\n\t"\
		"vmulpd		 0x080(%%rcx),%%ymm3,%%ymm3			\n\t		vmovaps		%%ymm10,%%ymm13					\n\t"\
	" vfmadd231pd	 0x0a0(%%rcx),%%ymm7,%%ymm2			\n\t		vmulpd		 0x180(%%rcx),%%ymm15,%%ymm15	/* c0C */\n\t"\
	"vfnmadd231pd	 0x0a0(%%rcx),%%ymm6,%%ymm3			\n\t		vmulpd		 0x180(%%rcx),%%ymm10,%%ymm10			\n\t"\
		"vmovaps	0x120(%%rbx),%%ymm7						\n\t	 vfmadd231pd	 0x1a0(%%rcx),%%ymm13,%%ymm15			\n\t"\
		"vmovaps	0x100(%%rbx),%%ymm6						\n\t	vfnmadd231pd	 0x1a0(%%rcx),%%ymm12,%%ymm10			\n\t"\
		"vmovaps	%%ymm3,0x0a0(%%rdx)						\n\t		vmovaps	0x1a0(%%rbx),%%ymm13						\n\t"\
		"vmovaps	%%ymm2,0x080(%%rdx)		\n\t/* r04,05 */\n\t		vmovaps	0x180(%%rbx),%%ymm12						\n\t"\
		"vmovaps	%%ymm7,0x2a0(%%rdx)						\n\t		vmovaps	%%ymm10,0x1a0(%%rdx)						\n\t"\
		"vmovaps	%%ymm6,0x280(%%rdx)		\n\t/* r14,15 */\n\t		vmovaps	%%ymm15,0x180(%%rdx)			\n\t"/* r0c,0d */\
		"vsubpd			%%ymm5,%%ymm0,%%ymm0			\n\t		vmovaps	%%ymm13,0x3a0(%%rdx)						\n\t"\
		"vaddpd			%%ymm4,%%ymm1,%%ymm1			\n\t		vmovaps	%%ymm12,0x380(%%rdx)			\n\t"/* r1c,1d */\
		"vmovaps		%%ymm0,%%ymm6					\n\t		vmovaps		%%ymm8 ,%%ymm12					\n\t"\
		"vmovaps		%%ymm1,%%ymm7					\n\t		vmovaps		%%ymm14,%%ymm13					\n\t"\
		"vmulpd		 0x0c0(%%rcx),%%ymm0,%%ymm0	/* c18*/\n\t		vmulpd		 0x1c0(%%rcx),%%ymm8 ,%%ymm8 	/* c1C */\n\t"\
		"vmulpd		 0x0c0(%%rcx),%%ymm1,%%ymm1			\n\t		vmulpd		 0x1c0(%%rcx),%%ymm14,%%ymm14			\n\t"\
	" vfmadd231pd	 0x0e0(%%rcx),%%ymm7,%%ymm0			\n\t	 vfmadd231pd	 0x1e0(%%rcx),%%ymm13,%%ymm8 			\n\t"\
	"vfnmadd231pd	 0x0e0(%%rcx),%%ymm6,%%ymm1			\n\t	vfnmadd231pd	 0x1e0(%%rcx),%%ymm12,%%ymm14			\n\t"\
		"vmovaps	0x320(%%rbx),%%ymm7						\n\t		vmovaps	0x3a0(%%rbx),%%ymm13						\n\t"\
		"vmovaps	0x300(%%rbx),%%ymm6						\n\t		vmovaps	0x380(%%rbx),%%ymm12						\n\t"\
		"vmovaps	%%ymm7,0x2e0(%%rdx)						\n\t		vmovaps	%%ymm13,0x3e0(%%rdx)						\n\t"\
		"vmovaps	%%ymm6,0x2c0(%%rdx)		/* r16,17 */	\n\t		vmovaps	%%ymm12,0x3c0(%%rdx)			\n\t"/* r1e,1f */\
		"vmovaps	%%ymm1,0x0e0(%%rdx)						\n\t		vmovaps	%%ymm14,0x1e0(%%rdx)						\n\t"\
		"vmovaps	%%ymm0,0x0c0(%%rdx)		/* r06,07 */	\n\t		vmovaps	%%ymm8 ,0x1c0(%%rdx)			\n\t"/* r0e,0f */\
		"\n\t"\
	/************************************************/				/************************************************/\
	/*...Block 3: t04,t14,t24,t34 -> r20,24,22,26:  */				/*...Block 7: t0C,t1C,t2C,t3C -> r28,2C,2A,2E:  */\
	/************************************************/				/************************************************/\
		"movq		%[__r20],%%rdx	\n\t"							/* base-addr in rcol = r28, so rdx offset +0x100 vs lcol */\
		"leaq	0x020(%%rsi),%%rcx	\n\t"/* cc0; Note cc0/ss0 are shared between lcol/rcol, so no rcx-offset until get to twiddles*/\
		"vmovaps		0x040(%%rdx),%%ymm4					\n\t		vmovaps		 0x140(%%rdx),%%ymm12			\n\t"\
		"vmovaps		0x0c0(%%rdx),%%ymm0					\n\t		vmovaps		 0x1c0(%%rdx),%%ymm8 			\n\t"\
		"vmovaps		0x060(%%rdx),%%ymm5					\n\t		vmovaps		 0x160(%%rdx),%%ymm13			\n\t"\
		"vmovaps		0x0e0(%%rdx),%%ymm3					\n\t		vmovaps		 0x1e0(%%rdx),%%ymm11			\n\t"\
		"vmovaps			%%ymm4 ,%%ymm6					\n\t		vmovaps		 	%%ymm12,%%ymm14				\n\t"\
		"vmovaps			%%ymm0 ,%%ymm2					\n\t		vmovaps		 	%%ymm8 ,%%ymm10				\n\t"\
		"vmovaps			%%ymm5 ,%%ymm7					\n\t		vmovaps		 	%%ymm13,%%ymm15				\n\t"\
	/*	"vmovaps			%%ymm1 ,%%ymm3					\n\t		vmovaps		 	%%ymm9 ,%%ymm11				\n\t"*/\
	/* Strategic register-usage for sincos consts here cuts loads from 16 to 4: */\
		"vmovaps			 (%%rcx),%%ymm9 			\n\t		vmovaps		0x020(%%rcx),%%ymm1 			\n\t"\
		"vmulpd		 	 %%ymm1 ,%%ymm7,%%ymm7				\n\t		vmulpd			 %%ymm9 ,%%ymm15,%%ymm15	\n\t"\
		"vmulpd		 	 %%ymm9 ,%%ymm3,%%ymm3				\n\t		vmulpd			 %%ymm1 ,%%ymm11,%%ymm11	\n\t"\
		"vmulpd		 	 %%ymm1 ,%%ymm6,%%ymm6				\n\t		vmulpd			 %%ymm9 ,%%ymm14,%%ymm14	\n\t"\
		"vmulpd		 	 %%ymm9 ,%%ymm2,%%ymm2				\n\t		vmulpd			 %%ymm1 ,%%ymm10,%%ymm10	\n\t"\
	"vfmadd132pd	 	 %%ymm9 ,%%ymm7,%%ymm4				\n\t	vfmadd132pd			 %%ymm1 ,%%ymm15,%%ymm12	\n\t"\
	"vfmadd132pd	 	 %%ymm1 ,%%ymm3,%%ymm0				\n\t	vfmadd132pd			 %%ymm9 ,%%ymm11,%%ymm8 	\n\t"\
	"vfmsub132pd	 	 %%ymm9 ,%%ymm6,%%ymm5				\n\t	vfmsub132pd			 %%ymm1 ,%%ymm14,%%ymm13	\n\t"\
	"vfmsub132pd	 0xe0(%%rdx),%%ymm2,%%ymm1				\n\t	vfmsub132pd		0x1e0(%%rdx),%%ymm10,%%ymm9 	\n\t"\
		"vmovaps		%%ymm5,%%ymm7					\n\t		vmovaps		%%ymm13,%%ymm15					\n\t"\
		"vmovaps		%%ymm4,%%ymm6					\n\t		vmovaps		%%ymm12,%%ymm14					\n\t"\
		"vaddpd			%%ymm0,%%ymm4,%%ymm4			\n\t		vaddpd			%%ymm8 ,%%ymm12,%%ymm12			\n\t"\
		"vaddpd			%%ymm1,%%ymm5,%%ymm5			\n\t		vaddpd			%%ymm9 ,%%ymm13,%%ymm13			\n\t"\
		"vsubpd			%%ymm0,%%ymm6,%%ymm6			\n\t		vsubpd			%%ymm8 ,%%ymm14,%%ymm14			\n\t"\
		"vsubpd			%%ymm1,%%ymm7,%%ymm7			\n\t		vsubpd			%%ymm9 ,%%ymm15,%%ymm15			\n\t"\
		"vmovaps		 0x080(%%rdx),%%ymm2				\n\t		vmovaps		 0x180(%%rdx),%%ymm10				\n\t"\
		"vmovaps		 0x0a0(%%rdx),%%ymm3				\n\t		vmovaps		 0x1a0(%%rdx),%%ymm11				\n\t"\
		"vmovaps			  (%%rdx),%%ymm0				\n\t		vmovaps		 0x100(%%rdx),%%ymm8 				\n\t"\
		"vmovaps		 0x020(%%rdx),%%ymm1				\n\t		vmovaps		 0x120(%%rdx),%%ymm9 				\n\t"\
		"vaddpd		 0x0a0(%%rdx),%%ymm2,%%ymm2			\n\t		vsubpd		 0x1a0(%%rdx),%%ymm10,%%ymm10			\n\t"\
		"vsubpd		 0x080(%%rdx),%%ymm3,%%ymm3			\n\t		vaddpd		 0x180(%%rdx),%%ymm11,%%ymm11			\n\t"\
		"vmulpd			  (%%rsi),%%ymm2,%%ymm2			\n\t		vmulpd			  (%%rsi),%%ymm10,%%ymm10			\n\t"\
		"vmulpd			  (%%rsi),%%ymm3,%%ymm3			\n\t		vmulpd			  (%%rsi),%%ymm11,%%ymm11			\n\t"\
		"vsubpd			%%ymm2,%%ymm0,%%ymm0			\n\t		vsubpd			%%ymm10,%%ymm8 ,%%ymm8 			\n\t"\
		"vsubpd			%%ymm3,%%ymm1,%%ymm1			\n\t		vsubpd			%%ymm11,%%ymm9 ,%%ymm9 			\n\t"\
	"vfmadd132pd		(%%rdi),%%ymm0,%%ymm2			\n\t	vfmadd132pd		(%%rdi),%%ymm8 ,%%ymm10			\n\t"\
	"vfmadd132pd		(%%rdi),%%ymm1,%%ymm3			\n\t	vfmadd132pd		(%%rdi),%%ymm9 ,%%ymm11			\n\t"\
		"movq		%[__c02],%%rcx				/* base-twiddle addr in rcol = c06, so rcx offset +0x100 vs lcol */	\n\t"\
		"vsubpd			%%ymm4,%%ymm2,%%ymm2			\n\t		vsubpd			%%ymm14,%%ymm8 ,%%ymm8 			\n\t"\
		"vsubpd			%%ymm5,%%ymm3,%%ymm3			\n\t		vsubpd			%%ymm15,%%ymm9 ,%%ymm9 			\n\t"\
	"vfmadd132pd		(%%rdi),%%ymm2,%%ymm4			\n\t	vfmadd132pd		(%%rdi),%%ymm8 ,%%ymm14			\n\t"\
	"vfmadd132pd		(%%rdi),%%ymm3,%%ymm5			\n\t	vfmadd132pd		(%%rdi),%%ymm9 ,%%ymm15			\n\t"\
		"vmovaps		%%ymm2, 0x040(%%rdx)			\n\t		vmovaps		%%ymm8 , 0x140(%%rdx)			\n\t"\
		"vmovaps		%%ymm3, 0x060(%%rdx)			\n\t		vmovaps		%%ymm9 , 0x160(%%rdx)			\n\t"\
		"vmovaps		%%ymm4,%%ymm2					\n\t		vmovaps		%%ymm14,%%ymm8 					\n\t"\
		"vmovaps		%%ymm5,%%ymm3					\n\t		vmovaps		%%ymm15,%%ymm9 					\n\t"\
		"vmulpd			  (%%rcx),%%ymm4,%%ymm4			\n\t		vmulpd		 0x100(%%rcx),%%ymm14,%%ymm14			\n\t"\
		"vmulpd			  (%%rcx),%%ymm5,%%ymm5			\n\t		vmulpd		 0x100(%%rcx),%%ymm15,%%ymm15			\n\t"\
	" vfmadd231pd	 0x020(%%rcx),%%ymm3,%%ymm4			\n\t	 vfmadd231pd	 0x120(%%rcx),%%ymm9 ,%%ymm14			\n\t"\
	"vfnmadd231pd	 0x020(%%rcx),%%ymm2,%%ymm5			\n\t	vfnmadd231pd	 0x120(%%rcx),%%ymm8 ,%%ymm15			\n\t"\
		"vmovaps	0x060(%%rbx),%%ymm3						\n\t		vmovaps	0x0e0(%%rbx),%%ymm9 						\n\t"\
		"vmovaps	0x040(%%rbx),%%ymm2						\n\t		vmovaps	0x0c0(%%rbx),%%ymm8 						\n\t"\
		"vmovaps	%%ymm5,0x020(%%rdx)						\n\t		vmovaps	%%ymm15,0x120(%%rdx)						\n\t"\
		"vmovaps	%%ymm4,     (%%rdx)			\n\t/* r20,21 */		vmovaps	%%ymm14,0x100(%%rdx)			\n\t"/* r28,29 */\
		"vmovaps	%%ymm3,0x220(%%rdx)						\n\t		vmovaps	%%ymm9 ,0x320(%%rdx)						\n\t"\
		"vmovaps	%%ymm2,0x200(%%rdx)			\n\t/* r30,31 */		vmovaps	%%ymm8 ,0x300(%%rdx)			\n\t"/* r38,39 */\
		"movq		%[__c12],%%rcx					\n\t"		/* rcol uses c16 */\
		"vmovaps		 0x040(%%rdx),%%ymm4				\n\t		vmovaps		 0x140(%%rdx),%%ymm14				\n\t"\
		"vmovaps		 0x060(%%rdx),%%ymm5				\n\t		vmovaps		 0x160(%%rdx),%%ymm15				\n\t"\
		"vmovaps		%%ymm4,%%ymm2					\n\t		vmovaps		%%ymm14,%%ymm8 					\n\t"\
		"vmovaps		%%ymm5,%%ymm3					\n\t		vmovaps		%%ymm15,%%ymm9 					\n\t"\
		"vmulpd			  (%%rcx),%%ymm4,%%ymm4			\n\t		vmulpd		 0x100(%%rcx),%%ymm14,%%ymm14			\n\t"\
		"vmulpd			  (%%rcx),%%ymm5,%%ymm5			\n\t		vmulpd		 0x100(%%rcx),%%ymm15,%%ymm15			\n\t"\
	" vfmadd231pd	 0x020(%%rcx),%%ymm3,%%ymm4			\n\t	 vfmadd231pd	 0x120(%%rcx),%%ymm9 ,%%ymm14			\n\t"\
	"vfnmadd231pd	 0x020(%%rcx),%%ymm2,%%ymm5			\n\t	vfnmadd231pd	 0x120(%%rcx),%%ymm8 ,%%ymm15			\n\t"\
		"vmovaps	0x260(%%rbx),%%ymm3						\n\t		vmovaps	0x2e0(%%rbx),%%ymm9 						\n\t"\
		"vmovaps	0x240(%%rbx),%%ymm2						\n\t		vmovaps	0x2c0(%%rbx),%%ymm8 						\n\t"\
		"vmovaps	%%ymm5,0x060(%%rdx)						\n\t		vmovaps	%%ymm15,0x160(%%rdx)						\n\t"\
		"vmovaps	%%ymm4,0x040(%%rdx)			\n\t/* r22,23 */		vmovaps	%%ymm14,0x140(%%rdx)			\n\t"/* r2a,2b */\
		"vmovaps	%%ymm3,0x260(%%rdx)						\n\t		vmovaps	%%ymm9 ,0x360(%%rdx)						\n\t"\
		"vmovaps	%%ymm2,0x240(%%rdx)			\n\t/* r32,33 */		vmovaps	%%ymm8 ,0x340(%%rdx)			\n\t"/* r3a,3b */\
		"movq		%[__c0A],%%rcx					\n\t"		/* rcol uses c0E */\
		"vsubpd			%%ymm7,%%ymm0,%%ymm0			\n\t		vsubpd			%%ymm13,%%ymm10,%%ymm10			\n\t"\
		"vsubpd			%%ymm6,%%ymm1,%%ymm1			\n\t		vsubpd			%%ymm12,%%ymm11,%%ymm11			\n\t"\
	"vfmadd132pd		(%%rdi),%%ymm0,%%ymm7			\n\t	vfmadd132pd		(%%rdi),%%ymm10,%%ymm13			\n\t"\
	"vfmadd132pd		(%%rdi),%%ymm1,%%ymm6			\n\t	vfmadd132pd		(%%rdi),%%ymm11,%%ymm12			\n\t"\
		"vmovaps		%%ymm7,%%ymm4					\n\t		vmovaps		%%ymm13,%%ymm8 					\n\t"\
		"vmovaps		%%ymm1,%%ymm5					\n\t		vmovaps		%%ymm11,%%ymm9 					\n\t"\
		"vmulpd			  (%%rcx),%%ymm7,%%ymm7			\n\t		vmulpd		 0x100(%%rcx),%%ymm13,%%ymm13			\n\t"\
		"vmulpd			  (%%rcx),%%ymm1,%%ymm1			\n\t		vmulpd		 0x100(%%rcx),%%ymm11,%%ymm11			\n\t"\
	" vfmadd231pd	 0x020(%%rcx),%%ymm5,%%ymm7			\n\t	 vfmadd231pd	 0x120(%%rcx),%%ymm9 ,%%ymm13			\n\t"\
	"vfnmadd231pd	 0x020(%%rcx),%%ymm4,%%ymm1			\n\t	vfnmadd231pd	 0x120(%%rcx),%%ymm8 ,%%ymm11			\n\t"\
		"vmovaps	0x160(%%rbx),%%ymm5						\n\t		vmovaps	0x1e0(%%rbx),%%ymm9 						\n\t"\
		"vmovaps	0x140(%%rbx),%%ymm4						\n\t		vmovaps	0x1c0(%%rbx),%%ymm8 						\n\t"\
		"vmovaps	%%ymm1,0x0a0(%%rdx)						\n\t		vmovaps	%%ymm11,0x1a0(%%rdx)						\n\t"\
		"vmovaps	%%ymm7,0x080(%%rdx)			\n\t/* r24,25 */		vmovaps	%%ymm13,0x180(%%rdx)			\n\t"/* r2c,2d */\
		"vmovaps	%%ymm5,0x2a0(%%rdx)						\n\t		vmovaps	%%ymm9 ,0x3a0(%%rdx)						\n\t"\
		"vmovaps	%%ymm4,0x280(%%rdx)			\n\t/* r34,35 */		vmovaps	%%ymm8 ,0x380(%%rdx)			\n\t"/* r3c,3d */\
		"movq		%[__c1A],%%rcx					\n\t"		/* rcol uses c1E */\
		"vmovaps		%%ymm0,%%ymm4					\n\t		vmovaps		%%ymm10,%%ymm8 					\n\t"\
		"vmovaps		%%ymm6,%%ymm5					\n\t		vmovaps		%%ymm12,%%ymm9 					\n\t"\
		"vmulpd			  (%%rcx),%%ymm0,%%ymm0			\n\t		vmulpd		 0x100(%%rcx),%%ymm10,%%ymm10			\n\t"\
		"vmulpd			  (%%rcx),%%ymm6,%%ymm6			\n\t		vmulpd		 0x100(%%rcx),%%ymm12,%%ymm12			\n\t"\
	" vfmadd231pd	 0x020(%%rcx),%%ymm5,%%ymm0			\n\t	 vfmadd231pd	 0x120(%%rcx),%%ymm9 ,%%ymm10			\n\t"\
	"vfnmadd231pd	 0x020(%%rcx),%%ymm4,%%ymm6			\n\t	vfnmadd231pd	 0x120(%%rcx),%%ymm8 ,%%ymm12			\n\t"\
		"vmovaps	0x360(%%rbx),%%ymm5						\n\t		vmovaps	0x3e0(%%rbx),%%ymm9 						\n\t"\
		"vmovaps	0x340(%%rbx),%%ymm4						\n\t		vmovaps	0x3c0(%%rbx),%%ymm8 						\n\t"\
		"vmovaps	%%ymm6,0x0e0(%%rdx)						\n\t		vmovaps	%%ymm12,0x1e0(%%rdx)						\n\t"\
		"vmovaps	%%ymm0,0x0c0(%%rdx)			\n\t/* r26,27 */		vmovaps	%%ymm10,0x1c0(%%rdx)			\n\t"/* r2e,2f */\
		"vmovaps	%%ymm5,0x2e0(%%rdx)						\n\t		vmovaps	%%ymm9 ,0x3e0(%%rdx)						\n\t"\
		"vmovaps	%%ymm4,0x2c0(%%rdx)			\n\t/* r36,37 */		vmovaps	%%ymm8 ,0x3c0(%%rdx)			\n\t"/* r3e,3f */\
	/*******************************************/\
	/**** Finish with 4-way 'un'terleaving: ****/\
	/*******************************************/\
		"movq	%[__r00] ,%%rsi\n\t"\
		"movq	%[__add0],%%rax\n\t"\
		"movq	%[__add1],%%rbx\n\t"\
		"movq	%[__add2],%%rcx\n\t"\
		"movq	%[__add3],%%rdx\n\t"\
	"/* a[j+p0]: Inputs from r00 +0/1, 8/9, 16/17, 24/25. Outputs into [add0,add0+0x100,add1,add1+0x100]+0x0: */	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps		 (%%rsi),%%ymm1						\n\t		vmovaps	0x020(%%rsi),%%ymm3								\n\t"\
		"vmovaps	0x200(%%rsi),%%ymm5						\n\t		vmovaps	0x220(%%rsi),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps	0x400(%%rsi),%%ymm6						\n\t		vmovaps	0x420(%%rsi),%%ymm14							\n\t"\
		"vmovaps	0x600(%%rsi),%%ymm5						\n\t		vmovaps	0x620(%%rsi),%%ymm13							\n\t"\
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
	"/* a[j+p2]: Inputs from r08 +0/1, 8/9, 16/17, 24/25. Outputs into [add0,add0+0x100,add1,add1+0x100]+0x0: */	\n\t"\
		"addq	$0x100,%%rsi	\n\t"\
		"addq	$0x40,%%rax	\n\t"\
		"addq	$0x40,%%rbx	\n\t"\
		"addq	$0x40,%%rcx	\n\t"\
		"addq	$0x40,%%rdx	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps		 (%%rsi),%%ymm1						\n\t		vmovaps	0x020(%%rsi),%%ymm3								\n\t"\
		"vmovaps	0x200(%%rsi),%%ymm5						\n\t		vmovaps	0x220(%%rsi),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps	0x400(%%rsi),%%ymm6						\n\t		vmovaps	0x420(%%rsi),%%ymm14							\n\t"\
		"vmovaps	0x600(%%rsi),%%ymm5						\n\t		vmovaps	0x620(%%rsi),%%ymm13							\n\t"\
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
	"/* a[j+p4]: Inputs from r04 +0/1, 8/9, 16/17, 24/25. Outputs into [add0,add0+0x100,add1,add1+0x100]+0x0: */	\n\t"\
		"subq	$0x80,%%rsi	\n\t"\
		"addq	$0x40,%%rax	\n\t"\
		"addq	$0x40,%%rbx	\n\t"\
		"addq	$0x40,%%rcx	\n\t"\
		"addq	$0x40,%%rdx	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps		 (%%rsi),%%ymm1						\n\t		vmovaps	0x020(%%rsi),%%ymm3								\n\t"\
		"vmovaps	0x200(%%rsi),%%ymm5						\n\t		vmovaps	0x220(%%rsi),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps	0x400(%%rsi),%%ymm6						\n\t		vmovaps	0x420(%%rsi),%%ymm14							\n\t"\
		"vmovaps	0x600(%%rsi),%%ymm5						\n\t		vmovaps	0x620(%%rsi),%%ymm13							\n\t"\
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
	"/* a[j+p6]: Inputs from r0c +0/1, 8/9, 16/17, 24/25. Outputs into [add0,add0+0x100,add1,add1+0x100]+0x0: */	\n\t"\
		"addq	$0x100,%%rsi	\n\t"\
		"addq	$0x40,%%rax	\n\t"\
		"addq	$0x40,%%rbx	\n\t"\
		"addq	$0x40,%%rcx	\n\t"\
		"addq	$0x40,%%rdx	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps		 (%%rsi),%%ymm1						\n\t		vmovaps	0x020(%%rsi),%%ymm3								\n\t"\
		"vmovaps	0x200(%%rsi),%%ymm5						\n\t		vmovaps	0x220(%%rsi),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps	0x400(%%rsi),%%ymm6						\n\t		vmovaps	0x420(%%rsi),%%ymm14							\n\t"\
		"vmovaps	0x600(%%rsi),%%ymm5						\n\t		vmovaps	0x620(%%rsi),%%ymm13							\n\t"\
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
	"/* a[j+p8]: Inputs from r02 +0/1, 8/9, 16/17, 24/25. Outputs into [add0,add0+0x100,add1,add1+0x100]+0x0: */	\n\t"\
		"subq	$0x140,%%rsi	\n\t"\
		"addq	$0x40,%%rax	\n\t"\
		"addq	$0x40,%%rbx	\n\t"\
		"addq	$0x40,%%rcx	\n\t"\
		"addq	$0x40,%%rdx	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps		 (%%rsi),%%ymm1						\n\t		vmovaps	0x020(%%rsi),%%ymm3								\n\t"\
		"vmovaps	0x200(%%rsi),%%ymm5						\n\t		vmovaps	0x220(%%rsi),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps	0x400(%%rsi),%%ymm6						\n\t		vmovaps	0x420(%%rsi),%%ymm14							\n\t"\
		"vmovaps	0x600(%%rsi),%%ymm5						\n\t		vmovaps	0x620(%%rsi),%%ymm13							\n\t"\
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
	"/* a[j+p10]: Inputs from r0a +0/1, 8/9, 16/17, 24/25. Outputs into [add0,add0+0x100,add1,add1+0x100]+0x0: */	\n\t"\
		"addq	$0x100,%%rsi	\n\t"\
		"addq	$0x40,%%rax	\n\t"\
		"addq	$0x40,%%rbx	\n\t"\
		"addq	$0x40,%%rcx	\n\t"\
		"addq	$0x40,%%rdx	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps		 (%%rsi),%%ymm1						\n\t		vmovaps	0x020(%%rsi),%%ymm3								\n\t"\
		"vmovaps	0x200(%%rsi),%%ymm5						\n\t		vmovaps	0x220(%%rsi),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps	0x400(%%rsi),%%ymm6						\n\t		vmovaps	0x420(%%rsi),%%ymm14							\n\t"\
		"vmovaps	0x600(%%rsi),%%ymm5						\n\t		vmovaps	0x620(%%rsi),%%ymm13							\n\t"\
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
	"/* a[j+p12]: Inputs from r06 +0/1, 8/9, 16/17, 24/25. Outputs into [add0,add0+0x100,add1,add1+0x100]+0x0: */	\n\t"\
		"subq	$0x80,%%rsi	\n\t"\
		"addq	$0x40,%%rax	\n\t"\
		"addq	$0x40,%%rbx	\n\t"\
		"addq	$0x40,%%rcx	\n\t"\
		"addq	$0x40,%%rdx	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps		 (%%rsi),%%ymm1						\n\t		vmovaps	0x020(%%rsi),%%ymm3								\n\t"\
		"vmovaps	0x200(%%rsi),%%ymm5						\n\t		vmovaps	0x220(%%rsi),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps	0x400(%%rsi),%%ymm6						\n\t		vmovaps	0x420(%%rsi),%%ymm14							\n\t"\
		"vmovaps	0x600(%%rsi),%%ymm5						\n\t		vmovaps	0x620(%%rsi),%%ymm13							\n\t"\
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
	"/* a[j+p14]: Inputs from r0e +0/1, 8/9, 16/17, 24/25. Outputs into [add0,add0+0x100,add1,add1+0x100]+0x0: */	\n\t"\
		"addq	$0x100,%%rsi	\n\t"\
		"addq	$0x40,%%rax	\n\t"\
		"addq	$0x40,%%rbx	\n\t"\
		"addq	$0x40,%%rcx	\n\t"\
		"addq	$0x40,%%rdx	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps		 (%%rsi),%%ymm1						\n\t		vmovaps	0x020(%%rsi),%%ymm3								\n\t"\
		"vmovaps	0x200(%%rsi),%%ymm5						\n\t		vmovaps	0x220(%%rsi),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps	0x400(%%rsi),%%ymm6						\n\t		vmovaps	0x420(%%rsi),%%ymm14							\n\t"\
		"vmovaps	0x600(%%rsi),%%ymm5						\n\t		vmovaps	0x620(%%rsi),%%ymm13							\n\t"\
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
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__r00] "m" (Xr00)\
		 ,[__r08] "m" (Xr08)\
		 ,[__r10] "m" (Xr10)\
		 ,[__r20] "m" (Xr20)\
		 ,[__r28] "m" (Xr28)\
		 ,[__r30] "m" (Xr30)\
		 ,[__c01] "m" (Xc01)\
		 ,[__c02] "m" (Xc02)\
		 ,[__c04] "m" (Xc04)\
		 ,[__c06] "m" (Xc06)\
		 ,[__c08] "m" (Xc08)\
		 ,[__c0A] "m" (Xc0A)\
		 ,[__c0C] "m" (Xc0C)\
		 ,[__c0E] "m" (Xc0E)\
		 ,[__c10] "m" (Xc10)\
		 ,[__c12] "m" (Xc12)\
		 ,[__c14] "m" (Xc14)\
		 ,[__c16] "m" (Xc16)\
		 ,[__c18] "m" (Xc18)\
		 ,[__c1A] "m" (Xc1A)\
		 ,[__c1C] "m" (Xc1C)\
		 ,[__c1E] "m" (Xc1E)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r10","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"		/* Clobbered registers */\
	);\
	}

  #endif	// ALL_FMA ?

#elif defined(USE_AVX)	// AVX and AVX2 both use 256-bit registers

	// Our non-FMA DIF and DIT implementations below have the following arithmetic opcounts - The disparate ADD/MUL balances between
	// DIF and DIT reflect use of 2 different variants for the middle step of the radix-2 in-place butterfly, shown in the right column:
	//		DIF: [222 ADD, 220 SUB, 322 MUL]		x -= y; y *= 2; y += x
	//		DIT: [327 ADD, 219 SUB, 212 MUL]		x -= y; y += y; y += x
	//
	#define SSE2_RADIX32_WRAPPER_DIF(Xadd0,Xadd1,Xadd2,Xadd3,Xr00,Xr10,Xr20,Xr30,Xisrt2,Xcc0,Xc00,Xc01,Xc02,Xc03,Xc05,Xc07)\
	{\
	__asm__ volatile (\
		"movq	%[__r00] ,%%rsi	\n\t"\
		"movq	%[__add0],%%rax\n\t"\
		"movq	%[__add1],%%rbx\n\t"\
		"movq	%[__add2],%%rcx\n\t"\
		"movq	%[__add3],%%rdx\n\t"\
	"/**** Start with 4-way interleaving: ****/\n\t"\
	"/* a[j+p0]: Inputs from [add0,add0+0x200,add1,add1+0x200]+0x0. Outputs into r0 +0/1, 16/17, 32/33, 48/49: */	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps		 (%%rax),%%ymm1						\n\t		vmovaps	0x020(%%rax),%%ymm3								\n\t"\
		"vmovaps		 (%%rbx),%%ymm5						\n\t		vmovaps	0x020(%%rbx),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps		 (%%rcx),%%ymm6						\n\t		vmovaps	0x020(%%rcx),%%ymm14							\n\t"\
		"vmovaps		 (%%rdx),%%ymm5						\n\t		vmovaps	0x020(%%rdx),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm6 ,%%ymm0				\n\t		vshufpd	$15 ,%%ymm13,%%ymm14 ,%%ymm2					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm6 ,%%ymm6				\n\t		vshufpd	$0  ,%%ymm13,%%ymm14 ,%%ymm14					\n\t"\
		"vperm2f128 $32 ,%%ymm0 ,%%ymm7 ,%%ymm5				\n\t		vperm2f128 $32 ,%%ymm2 ,%%ymm15,%%ymm13		/* outB	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm0 ,%%ymm7 ,%%ymm7				\n\t		vperm2f128 $49 ,%%ymm2 ,%%ymm15,%%ymm15		/* outD	*/	\n\t"\
		"vperm2f128 $32 ,%%ymm6 ,%%ymm1 ,%%ymm0				\n\t		vperm2f128 $32 ,%%ymm14,%%ymm3 ,%%ymm2 		/* outA	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm6 ,%%ymm1 ,%%ymm1				\n\t		vperm2f128 $49 ,%%ymm14,%%ymm3 ,%%ymm3 		/* outC	*/	\n\t"\
		"vmovaps 	%%ymm5 ,0x400(%%rsi)					\n\t		vmovaps %%ymm13,0x420(%%rsi)				/* outB	*/	\n\t"\
		"vmovaps 	%%ymm7 ,0x600(%%rsi)					\n\t		vmovaps %%ymm15,0x620(%%rsi)				/* outD	*/	\n\t"\
		"vmovaps 	%%ymm0 ,     (%%rsi)					\n\t		vmovaps %%ymm2 ,0x020(%%rsi)				/* outA	*/	\n\t"\
		"vmovaps 	%%ymm1 ,0x200(%%rsi)					\n\t		vmovaps %%ymm3 ,0x220(%%rsi)				/* outC	*/	\n\t"\
	"/* a[j+p2]: Inputs from [add0,add0+0x200,add1,add1+0x200]+0x40. Outputs into **r8** +0/1, 16/17, 32/33, 48/49: */	\n\t"\
		"addq	$0x40,%%rax	\n\t"\
		"addq	$0x40,%%rbx	\n\t"\
		"addq	$0x40,%%rcx	\n\t"\
		"addq	$0x40,%%rdx	\n\t"\
		"addq	$0x100,%%rsi	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps		 (%%rax),%%ymm1						\n\t		vmovaps	0x020(%%rax),%%ymm3								\n\t"\
		"vmovaps		 (%%rbx),%%ymm5						\n\t		vmovaps	0x020(%%rbx),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps		 (%%rcx),%%ymm6						\n\t		vmovaps	0x020(%%rcx),%%ymm14							\n\t"\
		"vmovaps		 (%%rdx),%%ymm5						\n\t		vmovaps	0x020(%%rdx),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm6 ,%%ymm0				\n\t		vshufpd	$15 ,%%ymm13,%%ymm14 ,%%ymm2					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm6 ,%%ymm6				\n\t		vshufpd	$0  ,%%ymm13,%%ymm14 ,%%ymm14					\n\t"\
		"vperm2f128 $32 ,%%ymm0 ,%%ymm7 ,%%ymm5				\n\t		vperm2f128 $32 ,%%ymm2 ,%%ymm15,%%ymm13		/* outB	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm0 ,%%ymm7 ,%%ymm7				\n\t		vperm2f128 $49 ,%%ymm2 ,%%ymm15,%%ymm15		/* outD	*/	\n\t"\
		"vperm2f128 $32 ,%%ymm6 ,%%ymm1 ,%%ymm0				\n\t		vperm2f128 $32 ,%%ymm14,%%ymm3 ,%%ymm2 		/* outA	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm6 ,%%ymm1 ,%%ymm1				\n\t		vperm2f128 $49 ,%%ymm14,%%ymm3 ,%%ymm3 		/* outC	*/	\n\t"\
		"vmovaps 	%%ymm5 ,0x400(%%rsi)					\n\t		vmovaps %%ymm13,0x420(%%rsi)				/* outB	*/	\n\t"\
		"vmovaps 	%%ymm7 ,0x600(%%rsi)					\n\t		vmovaps %%ymm15,0x620(%%rsi)				/* outD	*/	\n\t"\
		"vmovaps 	%%ymm0 ,     (%%rsi)					\n\t		vmovaps %%ymm2 ,0x020(%%rsi)				/* outA	*/	\n\t"\
		"vmovaps 	%%ymm1 ,0x200(%%rsi)					\n\t		vmovaps %%ymm3 ,0x220(%%rsi)				/* outC	*/	\n\t"\
	"/* a[j+p4]: Inputs from [add0,add0+0x200,add1,add1+0x200]+0x80. Outputs into r4 +0/1, 16/17, 32/33, 48/49: */	\n\t"\
		"addq	$0x40,%%rax	\n\t"\
		"addq	$0x40,%%rbx	\n\t"\
		"addq	$0x40,%%rcx	\n\t"\
		"addq	$0x40,%%rdx	\n\t"\
		"subq	$0x80,%%rsi	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps		 (%%rax),%%ymm1						\n\t		vmovaps	0x020(%%rax),%%ymm3								\n\t"\
		"vmovaps		 (%%rbx),%%ymm5						\n\t		vmovaps	0x020(%%rbx),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps		 (%%rcx),%%ymm6						\n\t		vmovaps	0x020(%%rcx),%%ymm14							\n\t"\
		"vmovaps		 (%%rdx),%%ymm5						\n\t		vmovaps	0x020(%%rdx),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm6 ,%%ymm0				\n\t		vshufpd	$15 ,%%ymm13,%%ymm14 ,%%ymm2					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm6 ,%%ymm6				\n\t		vshufpd	$0  ,%%ymm13,%%ymm14 ,%%ymm14					\n\t"\
		"vperm2f128 $32 ,%%ymm0 ,%%ymm7 ,%%ymm5				\n\t		vperm2f128 $32 ,%%ymm2 ,%%ymm15,%%ymm13		/* outB	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm0 ,%%ymm7 ,%%ymm7				\n\t		vperm2f128 $49 ,%%ymm2 ,%%ymm15,%%ymm15		/* outD	*/	\n\t"\
		"vperm2f128 $32 ,%%ymm6 ,%%ymm1 ,%%ymm0				\n\t		vperm2f128 $32 ,%%ymm14,%%ymm3 ,%%ymm2 		/* outA	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm6 ,%%ymm1 ,%%ymm1				\n\t		vperm2f128 $49 ,%%ymm14,%%ymm3 ,%%ymm3 		/* outC	*/	\n\t"\
		"vmovaps 	%%ymm5 ,0x400(%%rsi)					\n\t		vmovaps %%ymm13,0x420(%%rsi)				/* outB	*/	\n\t"\
		"vmovaps 	%%ymm7 ,0x600(%%rsi)					\n\t		vmovaps %%ymm15,0x620(%%rsi)				/* outD	*/	\n\t"\
		"vmovaps 	%%ymm0 ,     (%%rsi)					\n\t		vmovaps %%ymm2 ,0x020(%%rsi)				/* outA	*/	\n\t"\
		"vmovaps 	%%ymm1 ,0x200(%%rsi)					\n\t		vmovaps %%ymm3 ,0x220(%%rsi)				/* outC	*/	\n\t"\
	"/* a[j+p6]: Inputs from [add0,add0+0x200,add1,add1+0x200]+0xc0. Outputs into **r12** +0/1, 16/17, 32/33, 48/49: */	\n\t"\
		"addq	$0x40,%%rax	\n\t"\
		"addq	$0x40,%%rbx	\n\t"\
		"addq	$0x40,%%rcx	\n\t"\
		"addq	$0x40,%%rdx	\n\t"\
		"addq	$0x100,%%rsi	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps		 (%%rax),%%ymm1						\n\t		vmovaps	0x020(%%rax),%%ymm3								\n\t"\
		"vmovaps		 (%%rbx),%%ymm5						\n\t		vmovaps	0x020(%%rbx),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps		 (%%rcx),%%ymm6						\n\t		vmovaps	0x020(%%rcx),%%ymm14							\n\t"\
		"vmovaps		 (%%rdx),%%ymm5						\n\t		vmovaps	0x020(%%rdx),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm6 ,%%ymm0				\n\t		vshufpd	$15 ,%%ymm13,%%ymm14 ,%%ymm2					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm6 ,%%ymm6				\n\t		vshufpd	$0  ,%%ymm13,%%ymm14 ,%%ymm14					\n\t"\
		"vperm2f128 $32 ,%%ymm0 ,%%ymm7 ,%%ymm5				\n\t		vperm2f128 $32 ,%%ymm2 ,%%ymm15,%%ymm13		/* outB	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm0 ,%%ymm7 ,%%ymm7				\n\t		vperm2f128 $49 ,%%ymm2 ,%%ymm15,%%ymm15		/* outD	*/	\n\t"\
		"vperm2f128 $32 ,%%ymm6 ,%%ymm1 ,%%ymm0				\n\t		vperm2f128 $32 ,%%ymm14,%%ymm3 ,%%ymm2 		/* outA	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm6 ,%%ymm1 ,%%ymm1				\n\t		vperm2f128 $49 ,%%ymm14,%%ymm3 ,%%ymm3 		/* outC	*/	\n\t"\
		"vmovaps 	%%ymm5 ,0x400(%%rsi)					\n\t		vmovaps %%ymm13,0x420(%%rsi)				/* outB	*/	\n\t"\
		"vmovaps 	%%ymm7 ,0x600(%%rsi)					\n\t		vmovaps %%ymm15,0x620(%%rsi)				/* outD	*/	\n\t"\
		"vmovaps 	%%ymm0 ,     (%%rsi)					\n\t		vmovaps %%ymm2 ,0x020(%%rsi)				/* outA	*/	\n\t"\
		"vmovaps 	%%ymm1 ,0x200(%%rsi)					\n\t		vmovaps %%ymm3 ,0x220(%%rsi)				/* outC	*/	\n\t"\
	"/* a[j+p8]: Inputs from [add0,add0+0x200,add1,add1+0x200]+0x100. Outputs into **r2** +0/1, 16/17, 32/33, 48/49: */	\n\t"\
		"addq	$0x40,%%rax	\n\t"\
		"addq	$0x40,%%rbx	\n\t"\
		"addq	$0x40,%%rcx	\n\t"\
		"addq	$0x40,%%rdx	\n\t"\
		"subq	$0x140,%%rsi	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps		 (%%rax),%%ymm1						\n\t		vmovaps	0x020(%%rax),%%ymm3								\n\t"\
		"vmovaps		 (%%rbx),%%ymm5						\n\t		vmovaps	0x020(%%rbx),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps		 (%%rcx),%%ymm6						\n\t		vmovaps	0x020(%%rcx),%%ymm14							\n\t"\
		"vmovaps		 (%%rdx),%%ymm5						\n\t		vmovaps	0x020(%%rdx),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm6 ,%%ymm0				\n\t		vshufpd	$15 ,%%ymm13,%%ymm14 ,%%ymm2					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm6 ,%%ymm6				\n\t		vshufpd	$0  ,%%ymm13,%%ymm14 ,%%ymm14					\n\t"\
		"vperm2f128 $32 ,%%ymm0 ,%%ymm7 ,%%ymm5				\n\t		vperm2f128 $32 ,%%ymm2 ,%%ymm15,%%ymm13		/* outB	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm0 ,%%ymm7 ,%%ymm7				\n\t		vperm2f128 $49 ,%%ymm2 ,%%ymm15,%%ymm15		/* outD	*/	\n\t"\
		"vperm2f128 $32 ,%%ymm6 ,%%ymm1 ,%%ymm0				\n\t		vperm2f128 $32 ,%%ymm14,%%ymm3 ,%%ymm2 		/* outA	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm6 ,%%ymm1 ,%%ymm1				\n\t		vperm2f128 $49 ,%%ymm14,%%ymm3 ,%%ymm3 		/* outC	*/	\n\t"\
		"vmovaps 	%%ymm5 ,0x400(%%rsi)					\n\t		vmovaps %%ymm13,0x420(%%rsi)				/* outB	*/	\n\t"\
		"vmovaps 	%%ymm7 ,0x600(%%rsi)					\n\t		vmovaps %%ymm15,0x620(%%rsi)				/* outD	*/	\n\t"\
		"vmovaps 	%%ymm0 ,     (%%rsi)					\n\t		vmovaps %%ymm2 ,0x020(%%rsi)				/* outA	*/	\n\t"\
		"vmovaps 	%%ymm1 ,0x200(%%rsi)					\n\t		vmovaps %%ymm3 ,0x220(%%rsi)				/* outC	*/	\n\t"\
	"/* a[j+p10]: Inputs from [add0,add0+0x200,add1,add1+0x200]+0x140. Outputs into r10 +0/1, 16/17, 32/33, 48/49: */	\n\t"\
		"addq	$0x40,%%rax	\n\t"\
		"addq	$0x40,%%rbx	\n\t"\
		"addq	$0x40,%%rcx	\n\t"\
		"addq	$0x40,%%rdx	\n\t"\
		"addq	$0x100,%%rsi	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps		 (%%rax),%%ymm1						\n\t		vmovaps	0x020(%%rax),%%ymm3								\n\t"\
		"vmovaps		 (%%rbx),%%ymm5						\n\t		vmovaps	0x020(%%rbx),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps		 (%%rcx),%%ymm6						\n\t		vmovaps	0x020(%%rcx),%%ymm14							\n\t"\
		"vmovaps		 (%%rdx),%%ymm5						\n\t		vmovaps	0x020(%%rdx),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm6 ,%%ymm0				\n\t		vshufpd	$15 ,%%ymm13,%%ymm14 ,%%ymm2					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm6 ,%%ymm6				\n\t		vshufpd	$0  ,%%ymm13,%%ymm14 ,%%ymm14					\n\t"\
		"vperm2f128 $32 ,%%ymm0 ,%%ymm7 ,%%ymm5				\n\t		vperm2f128 $32 ,%%ymm2 ,%%ymm15,%%ymm13		/* outB	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm0 ,%%ymm7 ,%%ymm7				\n\t		vperm2f128 $49 ,%%ymm2 ,%%ymm15,%%ymm15		/* outD	*/	\n\t"\
		"vperm2f128 $32 ,%%ymm6 ,%%ymm1 ,%%ymm0				\n\t		vperm2f128 $32 ,%%ymm14,%%ymm3 ,%%ymm2 		/* outA	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm6 ,%%ymm1 ,%%ymm1				\n\t		vperm2f128 $49 ,%%ymm14,%%ymm3 ,%%ymm3 		/* outC	*/	\n\t"\
		"vmovaps 	%%ymm5 ,0x400(%%rsi)					\n\t		vmovaps %%ymm13,0x420(%%rsi)				/* outB	*/	\n\t"\
		"vmovaps 	%%ymm7 ,0x600(%%rsi)					\n\t		vmovaps %%ymm15,0x620(%%rsi)				/* outD	*/	\n\t"\
		"vmovaps 	%%ymm0 ,     (%%rsi)					\n\t		vmovaps %%ymm2 ,0x020(%%rsi)				/* outA	*/	\n\t"\
		"vmovaps 	%%ymm1 ,0x200(%%rsi)					\n\t		vmovaps %%ymm3 ,0x220(%%rsi)				/* outC	*/	\n\t"\
	"/* a[j+p12]: Inputs from [add0,add0+0x200,add1,add1+0x200]+0x180. Outputs into **r6** +0/1, 16/17, 32/33, 48/49: */	\n\t"\
		"addq	$0x40,%%rax	\n\t"\
		"addq	$0x40,%%rbx	\n\t"\
		"addq	$0x40,%%rcx	\n\t"\
		"addq	$0x40,%%rdx	\n\t"\
		"subq	$0x80,%%rsi	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps		 (%%rax),%%ymm1						\n\t		vmovaps	0x020(%%rax),%%ymm3								\n\t"\
		"vmovaps		 (%%rbx),%%ymm5						\n\t		vmovaps	0x020(%%rbx),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps		 (%%rcx),%%ymm6						\n\t		vmovaps	0x020(%%rcx),%%ymm14							\n\t"\
		"vmovaps		 (%%rdx),%%ymm5						\n\t		vmovaps	0x020(%%rdx),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm6 ,%%ymm0				\n\t		vshufpd	$15 ,%%ymm13,%%ymm14 ,%%ymm2					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm6 ,%%ymm6				\n\t		vshufpd	$0  ,%%ymm13,%%ymm14 ,%%ymm14					\n\t"\
		"vperm2f128 $32 ,%%ymm0 ,%%ymm7 ,%%ymm5				\n\t		vperm2f128 $32 ,%%ymm2 ,%%ymm15,%%ymm13		/* outB	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm0 ,%%ymm7 ,%%ymm7				\n\t		vperm2f128 $49 ,%%ymm2 ,%%ymm15,%%ymm15		/* outD	*/	\n\t"\
		"vperm2f128 $32 ,%%ymm6 ,%%ymm1 ,%%ymm0				\n\t		vperm2f128 $32 ,%%ymm14,%%ymm3 ,%%ymm2 		/* outA	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm6 ,%%ymm1 ,%%ymm1				\n\t		vperm2f128 $49 ,%%ymm14,%%ymm3 ,%%ymm3 		/* outC	*/	\n\t"\
		"vmovaps 	%%ymm5 ,0x400(%%rsi)					\n\t		vmovaps %%ymm13,0x420(%%rsi)				/* outB	*/	\n\t"\
		"vmovaps 	%%ymm7 ,0x600(%%rsi)					\n\t		vmovaps %%ymm15,0x620(%%rsi)				/* outD	*/	\n\t"\
		"vmovaps 	%%ymm0 ,     (%%rsi)					\n\t		vmovaps %%ymm2 ,0x020(%%rsi)				/* outA	*/	\n\t"\
		"vmovaps 	%%ymm1 ,0x200(%%rsi)					\n\t		vmovaps %%ymm3 ,0x220(%%rsi)				/* outC	*/	\n\t"\
	"/* a[j+p14]: Inputs from [add0,add0+0x200,add1,add1+0x200]+0x1c0. Outputs into r14 +0/1, 16/17, 32/33, 48/49: */	\n\t"\
		"addq	$0x40,%%rax	\n\t"\
		"addq	$0x40,%%rbx	\n\t"\
		"addq	$0x40,%%rcx	\n\t"\
		"addq	$0x40,%%rdx	\n\t"\
		"addq	$0x100,%%rsi	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps		 (%%rax),%%ymm1						\n\t		vmovaps	0x020(%%rax),%%ymm3								\n\t"\
		"vmovaps		 (%%rbx),%%ymm5						\n\t		vmovaps	0x020(%%rbx),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps		 (%%rcx),%%ymm6						\n\t		vmovaps	0x020(%%rcx),%%ymm14							\n\t"\
		"vmovaps		 (%%rdx),%%ymm5						\n\t		vmovaps	0x020(%%rdx),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm6 ,%%ymm0				\n\t		vshufpd	$15 ,%%ymm13,%%ymm14 ,%%ymm2					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm6 ,%%ymm6				\n\t		vshufpd	$0  ,%%ymm13,%%ymm14 ,%%ymm14					\n\t"\
		"vperm2f128 $32 ,%%ymm0 ,%%ymm7 ,%%ymm5				\n\t		vperm2f128 $32 ,%%ymm2 ,%%ymm15,%%ymm13		/* outB	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm0 ,%%ymm7 ,%%ymm7				\n\t		vperm2f128 $49 ,%%ymm2 ,%%ymm15,%%ymm15		/* outD	*/	\n\t"\
		"vperm2f128 $32 ,%%ymm6 ,%%ymm1 ,%%ymm0				\n\t		vperm2f128 $32 ,%%ymm14,%%ymm3 ,%%ymm2 		/* outA	*/	\n\t"\
		"vperm2f128 $49 ,%%ymm6 ,%%ymm1 ,%%ymm1				\n\t		vperm2f128 $49 ,%%ymm14,%%ymm3 ,%%ymm3 		/* outC	*/	\n\t"\
		"vmovaps 	%%ymm5 ,0x400(%%rsi)					\n\t		vmovaps %%ymm13,0x420(%%rsi)				/* outB	*/	\n\t"\
		"vmovaps 	%%ymm7 ,0x600(%%rsi)					\n\t		vmovaps %%ymm15,0x620(%%rsi)				/* outD	*/	\n\t"\
		"vmovaps 	%%ymm0 ,     (%%rsi)					\n\t		vmovaps %%ymm2 ,0x020(%%rsi)				/* outA	*/	\n\t"\
		"vmovaps 	%%ymm1 ,0x200(%%rsi)					\n\t		vmovaps %%ymm3 ,0x220(%%rsi)				/* outC	*/	\n\t"\
		"\n\t"\
		/************************************************************************/\
		/* Forward DIF radix-32 pass on the interleaved block1 and block2 data: */\
		/************************************************************************/\
		"/*...Block 0:	*/									\n\t"\
		"movq	%[__isrt2],%%rsi							\n\t		leaq	0x8e0(%%rsi),%%rdi	/* two */	\n\t"\
		"/*****	SSE2_RADIX4_DIF_4WRAPPER(c00,c08,c10,c18,r00) ****	SSE2_RADIX4_DIF_4WRAPPER_2NDOFTWO(c04,c0C,c14,c1C,r08)	*****/\n\t"\
		"movq		%[__r00],%%rcx						\n\t		/*addq		$0x100,%%rcx // __r08 */	\n\t"\
		"movq		%[__c00],%%rdx						\n\t		/*addq		$0x100,%%rdx // __c04 */	\n\t"\
		"vmovaps			 (%%rcx),%%ymm0					\n\t		vmovaps		0x100(%%rcx),%%ymm8			\n\t"\
		"vmovaps		0x020(%%rcx),%%ymm1					\n\t		vmovaps		0x120(%%rcx),%%ymm9			\n\t"\
		"vmovaps		%%ymm0,%%ymm2					\n\t		vmovaps		%%ymm8,%%ymm10				\n\t"\
		"vmovaps		%%ymm1,%%ymm3					\n\t		vmovaps		%%ymm9,%%ymm11				\n\t"\
		"vmulpd			 (%%rdx),%%ymm0,%%ymm0				\n\t		vmulpd		0x100(%%rdx),%%ymm8,%%ymm8			\n\t"\
		"vmulpd			 (%%rdx),%%ymm1,%%ymm1				\n\t		vmulpd		0x100(%%rdx),%%ymm9,%%ymm9			\n\t"\
		"vmulpd		0x020(%%rdx),%%ymm2,%%ymm2				\n\t		vmulpd		0x120(%%rdx),%%ymm10,%%ymm10		\n\t"\
		"vmulpd		0x020(%%rdx),%%ymm3,%%ymm3				\n\t		vmulpd		0x120(%%rdx),%%ymm11,%%ymm11		\n\t"\
		"vaddpd		%%ymm2,%%ymm1,%%ymm1				\n\t		vaddpd		%%ymm10,%%ymm9,%%ymm9				\n\t"\
		"vsubpd		%%ymm3,%%ymm0,%%ymm0				\n\t		vsubpd		%%ymm11,%%ymm8,%%ymm8				\n\t"\
		"vmovaps		0x040(%%rcx),%%ymm4					\n\t		vmovaps		0x140(%%rcx),%%ymm12				\n\t"\
		"vmovaps		0x060(%%rcx),%%ymm5					\n\t		vmovaps		0x160(%%rcx),%%ymm13				\n\t"\
		"vmovaps		%%ymm4,%%ymm6					\n\t		vmovaps		%%ymm12,%%ymm14						\n\t"\
		"vmovaps		%%ymm5,%%ymm7					\n\t		vmovaps		%%ymm13,%%ymm15						\n\t"\
		"addq		$0x040,%%rdx						\n\t"\
		"vmovaps		%%ymm0,%%ymm2					\n\t		vmovaps		%%ymm8,%%ymm10						\n\t"\
		"vmovaps		%%ymm1,%%ymm3					\n\t		vmovaps		%%ymm9,%%ymm11						\n\t"\
		"vmulpd			 (%%rdx),%%ymm4,%%ymm4				\n\t		vmulpd		0x100(%%rdx),%%ymm12,%%ymm12		\n\t"\
		"vmulpd			 (%%rdx),%%ymm5,%%ymm5				\n\t		vmulpd		0x100(%%rdx),%%ymm13,%%ymm13		\n\t"\
		"vmulpd		0x020(%%rdx),%%ymm6,%%ymm6				\n\t		vmulpd		0x120(%%rdx),%%ymm14,%%ymm14		\n\t"\
		"vmulpd		0x020(%%rdx),%%ymm7,%%ymm7				\n\t		vmulpd		0x120(%%rdx),%%ymm15,%%ymm15		\n\t"\
		"vaddpd		%%ymm6,%%ymm5,%%ymm5				\n\t		vaddpd		%%ymm14,%%ymm13,%%ymm13				\n\t"\
		"vsubpd		%%ymm7,%%ymm4,%%ymm4				\n\t		vsubpd		%%ymm15,%%ymm12,%%ymm12				\n\t"\
		"vaddpd		%%ymm4,%%ymm0,%%ymm0				\n\t		vaddpd		%%ymm12,%%ymm8,%%ymm8				\n\t"\
		"vaddpd		%%ymm5,%%ymm1,%%ymm1				\n\t		vaddpd		%%ymm13,%%ymm9,%%ymm9				\n\t"\
		"vsubpd		%%ymm4,%%ymm2,%%ymm2				\n\t		vsubpd		%%ymm12,%%ymm10,%%ymm10				\n\t"\
		"vsubpd		%%ymm5,%%ymm3,%%ymm3				\n\t		vsubpd		%%ymm13,%%ymm11,%%ymm11				\n\t"\
		"addq		$0x080,%%rdx						\n\t"\
		"vmovaps		0x0c0(%%rcx),%%ymm4					\n\t		vmovaps		0x1c0(%%rcx),%%ymm12				\n\t"\
		"vmovaps		0x0e0(%%rcx),%%ymm5					\n\t		vmovaps		0x1e0(%%rcx),%%ymm13				\n\t"\
		"vmovaps		%%ymm4,%%ymm6					\n\t		vmovaps		%%ymm12,%%ymm14						\n\t"\
		"vmovaps		%%ymm5,%%ymm7					\n\t		vmovaps		%%ymm13,%%ymm15						\n\t"\
		"vmulpd			 (%%rdx),%%ymm4,%%ymm4				\n\t		vmulpd		0x100(%%rdx),%%ymm12,%%ymm12		\n\t"\
		"vmulpd			 (%%rdx),%%ymm5,%%ymm5				\n\t		vmulpd		0x100(%%rdx),%%ymm13,%%ymm13		\n\t"\
		"vmulpd		0x020(%%rdx),%%ymm6,%%ymm6				\n\t		vmulpd		0x120(%%rdx),%%ymm14,%%ymm14		\n\t"\
		"vmulpd		0x020(%%rdx),%%ymm7,%%ymm7				\n\t		vmulpd		0x120(%%rdx),%%ymm15,%%ymm15		\n\t"\
		"vaddpd		%%ymm6,%%ymm5,%%ymm5				\n\t		vaddpd		%%ymm14,%%ymm13,%%ymm13				\n\t"\
		"vsubpd		%%ymm7,%%ymm4,%%ymm4				\n\t		vsubpd		%%ymm15,%%ymm12,%%ymm12				\n\t"\
		"vmovaps		%%ymm5,0x020(%%rcx)			\n\t		vmovaps		%%ymm13,0x120(%%rcx)				\n\t"\
		"vmovaps		%%ymm4,     (%%rcx)			\n\t		vmovaps		%%ymm12,0x100(%%rcx)				\n\t"\
		"subq		$0x040,%%rdx						\n\t"\
		"vmovaps		0x080(%%rcx),%%ymm4					\n\t		vmovaps		0x180(%%rcx),%%ymm12				\n\t"\
		"vmovaps		0x0a0(%%rcx),%%ymm5					\n\t		vmovaps		0x1a0(%%rcx),%%ymm13				\n\t"\
		"vmovaps		%%ymm4,%%ymm6					\n\t		vmovaps		%%ymm12,%%ymm14						\n\t"\
		"vmovaps		%%ymm5,%%ymm7					\n\t		vmovaps		%%ymm13,%%ymm15						\n\t"\
		"vmulpd			 (%%rdx),%%ymm4,%%ymm4				\n\t		vmulpd		0x100(%%rdx),%%ymm12,%%ymm12		\n\t"\
		"vmulpd			 (%%rdx),%%ymm5,%%ymm5				\n\t		vmulpd		0x100(%%rdx),%%ymm13,%%ymm13		\n\t"\
		"vmulpd		0x020(%%rdx),%%ymm6,%%ymm6				\n\t		vmulpd		0x120(%%rdx),%%ymm14,%%ymm14		\n\t"\
		"vmulpd		0x020(%%rdx),%%ymm7,%%ymm7				\n\t		vmulpd		0x120(%%rdx),%%ymm15,%%ymm15		\n\t"\
		"vaddpd		%%ymm6,%%ymm5,%%ymm5				\n\t		vaddpd		%%ymm14,%%ymm13,%%ymm13				\n\t"\
		"vsubpd		%%ymm7,%%ymm4,%%ymm4				\n\t		vsubpd		%%ymm15,%%ymm12,%%ymm12				\n\t"\
		"vmovaps		%%ymm5,%%ymm7					\n\t		vmovaps		%%ymm13,%%ymm15						\n\t"\
		"vmovaps		%%ymm4,%%ymm6					\n\t		vmovaps		%%ymm12,%%ymm14						\n\t"\
		"vsubpd			 (%%rcx),%%ymm4,%%ymm4				\n\t		vsubpd		0x100(%%rcx),%%ymm12,%%ymm12		\n\t"\
		"vsubpd		0x020(%%rcx),%%ymm5,%%ymm5				\n\t		vsubpd		0x120(%%rcx),%%ymm13,%%ymm13		\n\t"\
		"vaddpd			 (%%rcx),%%ymm6,%%ymm6				\n\t		vaddpd		0x100(%%rcx),%%ymm14,%%ymm14		\n\t"\
		"vaddpd		0x020(%%rcx),%%ymm7,%%ymm7				\n\t		vaddpd		0x120(%%rcx),%%ymm15,%%ymm15		\n\t"\
		"vsubpd		%%ymm6,%%ymm0,%%ymm0				\n\t		vsubpd		%%ymm14,%%ymm8,%%ymm8				\n\t"\
		"vsubpd		%%ymm7,%%ymm1,%%ymm1				\n\t		vsubpd		%%ymm15,%%ymm9,%%ymm9				\n\t"\
		"vmulpd		(%%rdi),%%ymm6,%%ymm6				\n\t		vmulpd		(%%rdi),%%ymm14,%%ymm14				\n\t"\
		"vmulpd		(%%rdi),%%ymm7,%%ymm7				\n\t		vmulpd		(%%rdi),%%ymm15,%%ymm15				\n\t"\
		"vaddpd		%%ymm0,%%ymm6,%%ymm6				\n\t		vaddpd		%%ymm8,%%ymm14,%%ymm14				\n\t"\
		"vaddpd		%%ymm1,%%ymm7,%%ymm7				\n\t		vaddpd		%%ymm9,%%ymm15,%%ymm15				\n\t"\
		"/*vmovaps	%%ymm6,     (%%rcx)	*/			\n\t		vmovaps		%%ymm14,0x100(%%rcx)				\n\t"\
		"/*vmovaps	%%ymm7,0x020(%%rcx)	*/			\n\t		vmovaps		%%ymm15,0x120(%%rcx)				\n\t"\
		"vsubpd		%%ymm5,%%ymm2,%%ymm2				\n\t		vsubpd		%%ymm13,%%ymm10,%%ymm10				\n\t"\
		"vsubpd		%%ymm4,%%ymm3,%%ymm3				\n\t		vsubpd		%%ymm12,%%ymm11,%%ymm11				\n\t"\
		"/*vmovaps	%%ymm2,0x040(%%rcx)	*/			\n\t		vmulpd		(%%rdi),%%ymm13,%%ymm13				\n\t"\
		"/*vmovaps	%%ymm3,0x0e0(%%rcx)	*/			\n\t		vmulpd		(%%rdi),%%ymm12,%%ymm12				\n\t"\
		"vmulpd		(%%rdi),%%ymm5,%%ymm5				\n\t		vaddpd		%%ymm10,%%ymm13,%%ymm13				\n\t"\
		"vmulpd		(%%rdi),%%ymm4,%%ymm4				\n\t		vaddpd		%%ymm11,%%ymm12,%%ymm12				\n\t"\
		"vaddpd		%%ymm2,%%ymm5,%%ymm5				\n\t		vmovaps		%%ymm10,%%ymm14						\n\t"\
		"vaddpd		%%ymm3,%%ymm4,%%ymm4				\n\t		vmovaps		%%ymm13,%%ymm15						\n\t"\
		"/*vmovaps	%%ymm5,0x0c0(%%rcx)	*/			\n\t		vsubpd		%%ymm12,%%ymm10,%%ymm10				\n\t"\
		"/*vmovaps	%%ymm4,0x060(%%rcx)	*/						vsubpd		%%ymm11,%%ymm13,%%ymm13				\n\t"\
		"																vaddpd		%%ymm12,%%ymm14,%%ymm14				\n\t"\
		"																vaddpd		%%ymm11,%%ymm15,%%ymm15				\n\t"\
		"																vmulpd		(%%rsi),%%ymm10,%%ymm10				\n\t"\
		"																vmulpd		(%%rsi),%%ymm13,%%ymm13				\n\t"\
		"																vmulpd		(%%rsi),%%ymm14,%%ymm14				\n\t"\
		"																vmulpd		(%%rsi),%%ymm15,%%ymm15				\n\t"\
		"/*****	SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r00) *****/\n\t"\
		"																vmovaps		0x100(%%rcx),%%ymm11		\n\t"\
		"																vmovaps		0x120(%%rcx),%%ymm12		\n\t"\
		"vsubpd		%%ymm11,%%ymm6,%%ymm6				\n\t		vsubpd		%%ymm10,%%ymm2,%%ymm2			\n\t"\
		"vsubpd		%%ymm9,%%ymm0,%%ymm0				\n\t		vsubpd		%%ymm15,%%ymm5,%%ymm5			\n\t"\
		"vsubpd		%%ymm12,%%ymm7,%%ymm7				\n\t		vsubpd		%%ymm14,%%ymm4,%%ymm4			\n\t"\
		"vsubpd		%%ymm8,%%ymm1,%%ymm1				\n\t		vsubpd		%%ymm13,%%ymm3,%%ymm3			\n\t"\
		"vmulpd		(%%rdi),%%ymm11,%%ymm11			\n\t		vmulpd		(%%rdi),%%ymm10,%%ymm10		\n\t"\
		"vmulpd		(%%rdi),%%ymm9,%%ymm9				\n\t		vmulpd		(%%rdi),%%ymm15,%%ymm15		\n\t"\
		"vmulpd		(%%rdi),%%ymm12,%%ymm12			\n\t		vmulpd		(%%rdi),%%ymm14,%%ymm14		\n\t"\
		"vmulpd		(%%rdi),%%ymm8,%%ymm8				\n\t		vmulpd		(%%rdi),%%ymm13,%%ymm13		\n\t"\
		"vaddpd		%%ymm6,%%ymm11,%%ymm11			\n\t		vaddpd		%%ymm2,%%ymm10,%%ymm10		\n\t"\
		"vaddpd		%%ymm0,%%ymm9,%%ymm9				\n\t		vaddpd		%%ymm5,%%ymm15,%%ymm15		\n\t"\
		"vaddpd		%%ymm7,%%ymm12,%%ymm12			\n\t		vaddpd		%%ymm4,%%ymm14,%%ymm14		\n\t"\
		"vaddpd		%%ymm1,%%ymm8,%%ymm8				\n\t		vaddpd		%%ymm3,%%ymm13,%%ymm13		\n\t"\
		"vmovaps		%%ymm6,0x100(%%rcx)			\n\t		vmovaps		%%ymm2,0x140(%%rcx)	\n\t"\
		"vmovaps		%%ymm0,0x080(%%rcx)			\n\t		vmovaps		%%ymm5,0x0c0(%%rcx)	\n\t"\
		"vmovaps		%%ymm7,0x120(%%rcx)			\n\t		vmovaps		%%ymm4,0x160(%%rcx)	\n\t"\
		"vmovaps		%%ymm1,0x1a0(%%rcx)			\n\t		vmovaps		%%ymm3,0x1e0(%%rcx)	\n\t"\
		"vmovaps		%%ymm11,     (%%rcx)			\n\t		vmovaps		%%ymm10,0x040(%%rcx)	\n\t"\
		"vmovaps		%%ymm9,0x180(%%rcx)			\n\t		vmovaps		%%ymm15,0x1c0(%%rcx)	\n\t"\
		"vmovaps		%%ymm12,0x020(%%rcx)			\n\t		vmovaps		%%ymm14,0x060(%%rcx)	\n\t"\
		"vmovaps		%%ymm8,0x0a0(%%rcx)			\n\t		vmovaps		%%ymm13,0x0e0(%%rcx)	\n\t"\
		"\n\t"\
		"/*...Block 2:	*/\n\t"\
		"/*****	SSE2_RADIX4_DIF_4WRAPPER(c02,c0A,c12,c1A,r10) ****	SSE2_RADIX4_DIF_4WRAPPER_2NDOFTWO(c06,c0E,c16,c1E,r18)	*****/\n\t"\
		"movq		%[__r10],%%rcx						\n\t		/*addq		$0x100,%%rcx // __r18 */	\n\t"\
		"movq		%[__c02],%%rdx						\n\t		/*addq		$0x100,%%rdx // __c06 */	\n\t"\
		"vmovaps			 (%%rcx),%%ymm0					\n\t		vmovaps		0x100(%%rcx),%%ymm8			\n\t"\
		"vmovaps		0x020(%%rcx),%%ymm1					\n\t		vmovaps		0x120(%%rcx),%%ymm9			\n\t"\
		"vmovaps		%%ymm0,%%ymm2					\n\t		vmovaps		%%ymm8,%%ymm10				\n\t"\
		"vmovaps		%%ymm1,%%ymm3					\n\t		vmovaps		%%ymm9,%%ymm11				\n\t"\
		"vmulpd			 (%%rdx),%%ymm0,%%ymm0				\n\t		vmulpd		0x100(%%rdx),%%ymm8,%%ymm8			\n\t"\
		"vmulpd			 (%%rdx),%%ymm1,%%ymm1				\n\t		vmulpd		0x100(%%rdx),%%ymm9,%%ymm9			\n\t"\
		"vmulpd		0x020(%%rdx),%%ymm2,%%ymm2				\n\t		vmulpd		0x120(%%rdx),%%ymm10,%%ymm10		\n\t"\
		"vmulpd		0x020(%%rdx),%%ymm3,%%ymm3				\n\t		vmulpd		0x120(%%rdx),%%ymm11,%%ymm11		\n\t"\
		"vaddpd		%%ymm2,%%ymm1,%%ymm1				\n\t		vaddpd		%%ymm10,%%ymm9,%%ymm9				\n\t"\
		"vsubpd		%%ymm3,%%ymm0,%%ymm0				\n\t		vsubpd		%%ymm11,%%ymm8,%%ymm8				\n\t"\
		"vmovaps		%%ymm0,%%ymm2					\n\t		vmovaps		%%ymm8,%%ymm10						\n\t"\
		"vmovaps		%%ymm1,%%ymm3					\n\t		vmovaps		%%ymm9,%%ymm11						\n\t"\
		"addq		$0x040,%%rdx						\n\t"\
		"vmovaps		0x040(%%rcx),%%ymm4					\n\t		vmovaps		0x140(%%rcx),%%ymm12				\n\t"\
		"vmovaps		0x060(%%rcx),%%ymm5					\n\t		vmovaps		0x160(%%rcx),%%ymm13				\n\t"\
		"vmovaps		%%ymm4,%%ymm6					\n\t		vmovaps		%%ymm12,%%ymm14						\n\t"\
		"vmovaps		%%ymm5,%%ymm7					\n\t		vmovaps		%%ymm13,%%ymm15						\n\t"\
		"vmulpd			 (%%rdx),%%ymm4,%%ymm4				\n\t		vmulpd		0x100(%%rdx),%%ymm12,%%ymm12		\n\t"\
		"vmulpd			 (%%rdx),%%ymm5,%%ymm5				\n\t		vmulpd		0x100(%%rdx),%%ymm13,%%ymm13		\n\t"\
		"vmulpd		0x020(%%rdx),%%ymm6,%%ymm6				\n\t		vmulpd		0x120(%%rdx),%%ymm14,%%ymm14		\n\t"\
		"vmulpd		0x020(%%rdx),%%ymm7,%%ymm7				\n\t		vmulpd		0x120(%%rdx),%%ymm15,%%ymm15		\n\t"\
		"vaddpd		%%ymm6,%%ymm5,%%ymm5				\n\t		vaddpd		%%ymm14,%%ymm13,%%ymm13				\n\t"\
		"vsubpd		%%ymm7,%%ymm4,%%ymm4				\n\t		vsubpd		%%ymm15,%%ymm12,%%ymm12				\n\t"\
		"vaddpd		%%ymm4,%%ymm0,%%ymm0				\n\t		vaddpd		%%ymm12,%%ymm8,%%ymm8				\n\t"\
		"vaddpd		%%ymm5,%%ymm1,%%ymm1				\n\t		vaddpd		%%ymm13,%%ymm9,%%ymm9				\n\t"\
		"vsubpd		%%ymm4,%%ymm2,%%ymm2				\n\t		vsubpd		%%ymm12,%%ymm10,%%ymm10				\n\t"\
		"vsubpd		%%ymm5,%%ymm3,%%ymm3				\n\t		vsubpd		%%ymm13,%%ymm11,%%ymm11				\n\t"\
		"addq		$0x080,%%rdx						\n\t"\
		"vmovaps		0x0c0(%%rcx),%%ymm4					\n\t		vmovaps		0x1c0(%%rcx),%%ymm12				\n\t"\
		"vmovaps		0x0e0(%%rcx),%%ymm5					\n\t		vmovaps		0x1e0(%%rcx),%%ymm13				\n\t"\
		"vmovaps		%%ymm4,%%ymm6					\n\t		vmovaps		%%ymm12,%%ymm14						\n\t"\
		"vmovaps		%%ymm5,%%ymm7					\n\t		vmovaps		%%ymm13,%%ymm15						\n\t"\
		"vmulpd			 (%%rdx),%%ymm4,%%ymm4				\n\t		vmulpd		0x100(%%rdx),%%ymm12,%%ymm12		\n\t"\
		"vmulpd			 (%%rdx),%%ymm5,%%ymm5				\n\t		vmulpd		0x100(%%rdx),%%ymm13,%%ymm13		\n\t"\
		"vmulpd		0x020(%%rdx),%%ymm6,%%ymm6				\n\t		vmulpd		0x120(%%rdx),%%ymm14,%%ymm14		\n\t"\
		"vmulpd		0x020(%%rdx),%%ymm7,%%ymm7				\n\t		vmulpd		0x120(%%rdx),%%ymm15,%%ymm15		\n\t"\
		"vaddpd		%%ymm6,%%ymm5,%%ymm5				\n\t		vaddpd		%%ymm14,%%ymm13,%%ymm13				\n\t"\
		"vsubpd		%%ymm7,%%ymm4,%%ymm4				\n\t		vsubpd		%%ymm15,%%ymm12,%%ymm12				\n\t"\
		"vmovaps		%%ymm5,0x020(%%rcx)			\n\t		vmovaps		%%ymm13,0x120(%%rcx)				\n\t"\
		"vmovaps		%%ymm4,     (%%rcx)			\n\t		vmovaps		%%ymm12,0x100(%%rcx)				\n\t"\
		"subq		$0x040,%%rdx						\n\t"\
		"vmovaps		0x080(%%rcx),%%ymm4					\n\t		vmovaps		0x180(%%rcx),%%ymm12				\n\t"\
		"vmovaps		0x0a0(%%rcx),%%ymm5					\n\t		vmovaps		0x1a0(%%rcx),%%ymm13				\n\t"\
		"vmovaps		%%ymm4,%%ymm6					\n\t		vmovaps		%%ymm12,%%ymm14						\n\t"\
		"vmovaps		%%ymm5,%%ymm7					\n\t		vmovaps		%%ymm13,%%ymm15						\n\t"\
		"vmulpd			 (%%rdx),%%ymm4,%%ymm4				\n\t		vmulpd		0x100(%%rdx),%%ymm12,%%ymm12		\n\t"\
		"vmulpd			 (%%rdx),%%ymm5,%%ymm5				\n\t		vmulpd		0x100(%%rdx),%%ymm13,%%ymm13		\n\t"\
		"vmulpd		0x020(%%rdx),%%ymm6,%%ymm6				\n\t		vmulpd		0x120(%%rdx),%%ymm14,%%ymm14		\n\t"\
		"vmulpd		0x020(%%rdx),%%ymm7,%%ymm7				\n\t		vmulpd		0x120(%%rdx),%%ymm15,%%ymm15		\n\t"\
		"vaddpd		%%ymm6,%%ymm5,%%ymm5				\n\t		vaddpd		%%ymm14,%%ymm13,%%ymm13				\n\t"\
		"vsubpd		%%ymm7,%%ymm4,%%ymm4				\n\t		vsubpd		%%ymm15,%%ymm12,%%ymm12				\n\t"\
		"vmovaps		%%ymm5,%%ymm7					\n\t		vmovaps		%%ymm13,%%ymm15						\n\t"\
		"vmovaps		%%ymm4,%%ymm6					\n\t		vmovaps		%%ymm12,%%ymm14						\n\t"\
		"vsubpd			 (%%rcx),%%ymm4,%%ymm4				\n\t		vsubpd		0x100(%%rcx),%%ymm12,%%ymm12		\n\t"\
		"vsubpd		0x020(%%rcx),%%ymm5,%%ymm5				\n\t		vsubpd		0x120(%%rcx),%%ymm13,%%ymm13		\n\t"\
		"vaddpd			 (%%rcx),%%ymm6,%%ymm6				\n\t		vaddpd		0x100(%%rcx),%%ymm14,%%ymm14		\n\t"\
		"vaddpd		0x020(%%rcx),%%ymm7,%%ymm7				\n\t		vaddpd		0x120(%%rcx),%%ymm15,%%ymm15		\n\t"\
		"vsubpd		%%ymm6,%%ymm0,%%ymm0				\n\t		vsubpd		%%ymm14,%%ymm8,%%ymm8				\n\t"\
		"vsubpd		%%ymm7,%%ymm1,%%ymm1				\n\t		vsubpd		%%ymm15,%%ymm9,%%ymm9				\n\t"\
		"vmulpd		(%%rdi),%%ymm6,%%ymm6				\n\t		vmulpd		(%%rdi),%%ymm14,%%ymm14				\n\t"\
		"vmulpd		(%%rdi),%%ymm7,%%ymm7				\n\t		vmulpd		(%%rdi),%%ymm15,%%ymm15				\n\t"\
		"vaddpd		%%ymm0,%%ymm6,%%ymm6				\n\t		vaddpd		%%ymm8 ,%%ymm14,%%ymm14				\n\t"\
		"vaddpd		%%ymm1,%%ymm7,%%ymm7				\n\t		vaddpd		%%ymm9 ,%%ymm15,%%ymm15				\n\t"\
		"/*vmovaps	%%ymm6,     (%%rcx)	*/			\n\t		vmovaps		%%ymm14,0x100(%%rcx)				\n\t"\
		"/*vmovaps	%%ymm7,0x020(%%rcx)	*/			\n\t		vmovaps		%%ymm15,0x120(%%rcx)				\n\t"\
		"vsubpd		%%ymm5,%%ymm2,%%ymm2				\n\t		vsubpd		%%ymm13,%%ymm10,%%ymm10				\n\t"\
		"vsubpd		%%ymm4,%%ymm3,%%ymm3				\n\t		vsubpd		%%ymm12,%%ymm11,%%ymm11				\n\t"\
		"/*vmovaps	%%ymm2,0x040(%%rcx)	*/			\n\t		vmulpd		(%%rdi),%%ymm13,%%ymm13				\n\t"\
		"/*vmovaps	%%ymm3,0x0e0(%%rcx)	*/			\n\t		vmulpd		(%%rdi),%%ymm12,%%ymm12				\n\t"\
		"vmulpd		(%%rdi),%%ymm5,%%ymm5				\n\t		vaddpd		%%ymm10,%%ymm13,%%ymm13				\n\t"\
		"vmulpd		(%%rdi),%%ymm4,%%ymm4				\n\t		vaddpd		%%ymm11,%%ymm12,%%ymm12				\n\t"\
		"vaddpd		%%ymm2,%%ymm5,%%ymm5				\n\t		vmovaps		%%ymm10,%%ymm14						\n\t"\
		"vaddpd		%%ymm3,%%ymm4,%%ymm4				\n\t		vmovaps		%%ymm13,%%ymm15						\n\t"\
		"/*vmovaps	%%ymm5,0x0c0(%%rcx)	*/			\n\t		vsubpd		%%ymm12,%%ymm10,%%ymm10				\n\t"\
		"/*vmovaps	%%ymm4,0x060(%%rcx)	*/			\n\t		vsubpd		%%ymm11,%%ymm13,%%ymm13				\n\t"\
		"																vaddpd		%%ymm12,%%ymm14,%%ymm14				\n\t"\
		"																vaddpd		%%ymm11,%%ymm15,%%ymm15				\n\t"\
		"																vmulpd		(%%rsi),%%ymm10,%%ymm10	/* isrt2 */	\n\t"\
		"																vmulpd		(%%rsi),%%ymm13,%%ymm13				\n\t"\
		"																vmulpd		(%%rsi),%%ymm14,%%ymm14				\n\t"\
		"																vmulpd		(%%rsi),%%ymm15,%%ymm15				\n\t"\
		"/*****	SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r10) *****/\n\t"\
		"																vmovaps		0x100(%%rcx),%%ymm11		\n\t"\
		"																vmovaps		0x120(%%rcx),%%ymm12		\n\t"\
		"vsubpd		%%ymm11,%%ymm6,%%ymm6				\n\t		vsubpd		%%ymm10,%%ymm2,%%ymm2			\n\t"\
		"vsubpd		%%ymm9,%%ymm0,%%ymm0				\n\t		vsubpd		%%ymm15,%%ymm5,%%ymm5			\n\t"\
		"vsubpd		%%ymm12,%%ymm7,%%ymm7				\n\t		vsubpd		%%ymm14,%%ymm4,%%ymm4			\n\t"\
		"vsubpd		%%ymm8,%%ymm1,%%ymm1				\n\t		vsubpd		%%ymm13,%%ymm3,%%ymm3			\n\t"\
		"vmulpd		(%%rdi),%%ymm11,%%ymm11			\n\t		vmulpd		(%%rdi),%%ymm10,%%ymm10		\n\t"\
		"vmulpd		(%%rdi),%%ymm9,%%ymm9				\n\t		vmulpd		(%%rdi),%%ymm15,%%ymm15		\n\t"\
		"vmulpd		(%%rdi),%%ymm12,%%ymm12			\n\t		vmulpd		(%%rdi),%%ymm14,%%ymm14		\n\t"\
		"vmulpd		(%%rdi),%%ymm8,%%ymm8				\n\t		vmulpd		(%%rdi),%%ymm13,%%ymm13		\n\t"\
		"vaddpd		%%ymm6,%%ymm11,%%ymm11			\n\t		vaddpd		%%ymm2,%%ymm10,%%ymm10		\n\t"\
		"vaddpd		%%ymm0,%%ymm9,%%ymm9				\n\t		vaddpd		%%ymm5,%%ymm15,%%ymm15		\n\t"\
		"vaddpd		%%ymm7,%%ymm12,%%ymm12			\n\t		vaddpd		%%ymm4,%%ymm14,%%ymm14		\n\t"\
		"vaddpd		%%ymm1,%%ymm8,%%ymm8				\n\t		vaddpd		%%ymm3,%%ymm13,%%ymm13		\n\t"\
		"vmovaps		%%ymm6,0x100(%%rcx)			\n\t		vmovaps		%%ymm2,0x140(%%rcx)	\n\t"\
		"vmovaps		%%ymm0,0x080(%%rcx)			\n\t		vmovaps		%%ymm5,0x0c0(%%rcx)	\n\t"\
		"vmovaps		%%ymm7,0x120(%%rcx)			\n\t		vmovaps		%%ymm4,0x160(%%rcx)	\n\t"\
		"vmovaps		%%ymm1,0x1a0(%%rcx)			\n\t		vmovaps		%%ymm3,0x1e0(%%rcx)	\n\t"\
		"vmovaps		%%ymm11,     (%%rcx)			\n\t		vmovaps		%%ymm10,0x040(%%rcx)	\n\t"\
		"vmovaps		%%ymm9,0x180(%%rcx)			\n\t		vmovaps		%%ymm15,0x1c0(%%rcx)	\n\t"\
		"vmovaps		%%ymm12,0x020(%%rcx)			\n\t		vmovaps		%%ymm14,0x060(%%rcx)	\n\t"\
		"vmovaps		%%ymm8,0x0a0(%%rcx)			\n\t		vmovaps		%%ymm13,0x0e0(%%rcx)	\n\t"\
		/********************************************************************************************************\
		 Next 2 blocks operate on odd-indexed elements from the unpck*pd commands which we stored to temporaries:\
		********************************************************************************************************/\
		"/*...Block 3:	*/\n\t"\
		"/*	SSE2_RADIX4_DIF_4TWIDDLE         (r20,r24,r22,r26,r20,c01) */\n\t"\
		"addq		$0x200,%%rcx		/*****	SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO(r28,r2C,r2A,r2E,r28,c05)	*****/\n\t"\
		"movq		%[__c01],%%rbx						\n\t	/*	movq		%[__c05],%%rbx	*/		\n\t"\
		"movq		%%rcx,%%rax						\n\t	/*	addq		$0x080,%%rax	*/		\n\t"\
		"addq		$0x040,%%rcx						\n\t	/*	addq		$0x080,%%rcx	*/		\n\t"\
		"vmovaps		 (%%rax),%%ymm0	\n\t movq %%rax,%%rdx \n\t	vmovaps		0x100(%%rax),%%ymm8			\n\t"\
		"vmovaps		 (%%rcx),%%ymm4						\n\t		vmovaps		0x100(%%rcx),%%ymm12		\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1						\n\t		vmovaps		0x120(%%rax),%%ymm9			\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5						\n\t		vmovaps		0x120(%%rcx),%%ymm13		\n\t"\
		"vmovaps		 (%%rbx),%%ymm6						\n\t		vmovaps		0x100(%%rbx),%%ymm14		\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm7						\n\t		vmovaps		0x120(%%rbx),%%ymm15		\n\t"\
		"vmovaps	%%ymm0,%%ymm2						\n\t		vmovaps		%%ymm8,%%ymm10		\n\t"\
		"vmovaps	%%ymm1,%%ymm3						\n\t		vmovaps		%%ymm9,%%ymm11		\n\t"\
		"vmulpd		%%ymm6,%%ymm0,%%ymm0				\n\t		vmulpd		%%ymm14,%%ymm8,%%ymm8			\n\t"\
		"vmulpd		%%ymm6,%%ymm1,%%ymm1				\n\t		vmulpd		%%ymm14,%%ymm9,%%ymm9			\n\t"\
		"vmulpd		%%ymm7,%%ymm2,%%ymm2				\n\t		vmulpd		%%ymm15,%%ymm10,%%ymm10		\n\t"\
		"vmulpd		%%ymm7,%%ymm3,%%ymm3				\n\t		vmulpd		%%ymm15,%%ymm11,%%ymm11		\n\t"\
		"vmovaps		%%ymm4,%%ymm6					\n\t		vmovaps		%%ymm12,%%ymm14				\n\t"\
		"vaddpd		%%ymm2,%%ymm1,%%ymm1				\n\t		vaddpd		%%ymm10,%%ymm9,%%ymm9			\n\t"\
		"vmovaps		%%ymm5,%%ymm7					\n\t		vmovaps		%%ymm13,%%ymm15				\n\t"\
		"vmulpd		0x040(%%rbx),%%ymm4,%%ymm4				\n\t		vmulpd		0x140(%%rbx),%%ymm12,%%ymm12		\n\t"\
		"vsubpd		%%ymm3,%%ymm0,%%ymm0				\n\t		vsubpd		%%ymm11,%%ymm8,%%ymm8			\n\t"\
		"vmulpd		0x040(%%rbx),%%ymm5,%%ymm5				\n\t		vmulpd		0x140(%%rbx),%%ymm13,%%ymm13		\n\t"\
		"vmulpd		0x060(%%rbx),%%ymm6,%%ymm6				\n\t		vmulpd		0x160(%%rbx),%%ymm14,%%ymm14		\n\t"\
		"vmovaps		%%ymm0,%%ymm2					\n\t		vmovaps		%%ymm8,%%ymm10				\n\t"\
		"vmulpd		0x060(%%rbx),%%ymm7,%%ymm7				\n\t		vmulpd		0x160(%%rbx),%%ymm15,%%ymm15		\n\t"\
		"vaddpd		%%ymm6,%%ymm5,%%ymm5				\n\t		vaddpd		%%ymm14,%%ymm13,%%ymm13		\n\t"\
		"vmovaps		%%ymm1,%%ymm3					\n\t		vmovaps		%%ymm9,%%ymm11				\n\t"\
		"vsubpd		%%ymm7,%%ymm4,%%ymm4				\n\t		vsubpd		%%ymm15,%%ymm12,%%ymm12		\n\t"\
		"addq		$0x080,%%rcx						\n\t		vaddpd		%%ymm12,%%ymm8,%%ymm8			\n\t"\
		"addq		$0x0c0,%%rbx						\n\t		vaddpd		%%ymm13,%%ymm9,%%ymm9			\n\t"\
		"vmovaps			 (%%rcx),%%ymm6					\n\t		vsubpd		%%ymm12,%%ymm10,%%ymm10		\n\t"\
		"vmovaps		0x020(%%rcx),%%ymm7					\n\t		vsubpd		%%ymm13,%%ymm11,%%ymm11		\n\t"\
		"vaddpd		%%ymm4,%%ymm0,%%ymm0				\n\t		vmovaps		0x100(%%rcx),%%ymm12				\n\t"\
		"vaddpd		%%ymm5,%%ymm1,%%ymm1				\n\t		vmovaps		0x120(%%rcx),%%ymm13				\n\t"\
		"vsubpd		%%ymm4,%%ymm2,%%ymm2				\n\t		vmovaps		0x100(%%rcx),%%ymm14				\n\t"\
		"vsubpd		%%ymm5,%%ymm3,%%ymm3				\n\t		vmovaps		0x120(%%rcx),%%ymm15				\n\t"\
		"vmovaps		%%ymm6,%%ymm4					\n\t		vmulpd		0x100(%%rbx),%%ymm12,%%ymm12		\n\t"\
		"vmovaps		%%ymm7,%%ymm5					\n\t		vmulpd		0x100(%%rbx),%%ymm13,%%ymm13		\n\t"\
		"vmulpd			 (%%rbx),%%ymm4,%%ymm4				\n\t		vmulpd		0x120(%%rbx),%%ymm14,%%ymm14		\n\t"\
		"vmulpd			 (%%rbx),%%ymm5,%%ymm5				\n\t		vmulpd		0x120(%%rbx),%%ymm15,%%ymm15		\n\t"\
		"vmulpd		0x020(%%rbx),%%ymm6,%%ymm6				\n\t		vaddpd		%%ymm14,%%ymm13,%%ymm13		\n\t"\
		"vmulpd		0x020(%%rbx),%%ymm7,%%ymm7				\n\t		vsubpd		%%ymm15,%%ymm12,%%ymm12		\n\t"\
		"/*movq		%%rax,%%rdx		*/				\n\t		vmovaps		%%ymm13,0x120(%%rdx)				\n\t"\
		"vaddpd		%%ymm6,%%ymm5,%%ymm5				\n\t		vmovaps		%%ymm12,0x100(%%rdx)				\n\t"\
		"vsubpd		%%ymm7,%%ymm4,%%ymm4				\n\t"\
		"vmovaps		%%ymm5,0x020(%%rdx)			\n\t		addq	$0x080,%%rax							\n\t"\
		"vmovaps		%%ymm4,     (%%rdx)			\n\t		subq	$0x040,%%rbx							\n\t"\
		"vmovaps			 (%%rax),%%ymm4					\n\t		vmovaps		0x100(%%rax),%%ymm12				\n\t"\
		"vmovaps		0x020(%%rax),%%ymm5					\n\t		vmovaps		0x120(%%rax),%%ymm13				\n\t"\
		"vmovaps		%%ymm4,%%ymm6					\n\t		vmovaps		0x100(%%rax),%%ymm14				\n\t"\
		"vmovaps		%%ymm5,%%ymm7					\n\t		vmovaps		0x120(%%rax),%%ymm15				\n\t"\
		"vmulpd			 (%%rbx),%%ymm4,%%ymm4				\n\t		vmulpd		0x100(%%rbx),%%ymm12,%%ymm12		\n\t"\
		"vmulpd			 (%%rbx),%%ymm5,%%ymm5				\n\t		vmulpd		0x100(%%rbx),%%ymm13,%%ymm13		\n\t"\
		"vmulpd		0x020(%%rbx),%%ymm6,%%ymm6				\n\t		vmulpd		0x120(%%rbx),%%ymm14,%%ymm14		\n\t"\
		"vmulpd		0x020(%%rbx),%%ymm7,%%ymm7				\n\t		vmulpd		0x120(%%rbx),%%ymm15,%%ymm15		\n\t"\
		"vaddpd		%%ymm6,%%ymm5,%%ymm5				\n\t		vaddpd		%%ymm14,%%ymm13,%%ymm13		\n\t"\
		"vsubpd		%%ymm7,%%ymm4,%%ymm4				\n\t		vsubpd		%%ymm15,%%ymm12,%%ymm12		\n\t"\
		"vmovaps		%%ymm5,%%ymm7					\n\t		vmovaps		%%ymm13,%%ymm15				\n\t"\
		"vmovaps		%%ymm4,%%ymm6					\n\t		vmovaps		%%ymm12,%%ymm14				\n\t"\
		"vsubpd			 (%%rdx),%%ymm4,%%ymm4				\n\t		vsubpd		0x100(%%rdx),%%ymm12,%%ymm12		\n\t"\
		"vsubpd		0x020(%%rdx),%%ymm5,%%ymm5				\n\t		vsubpd		0x120(%%rdx),%%ymm13,%%ymm13		\n\t"\
		"vaddpd			 (%%rdx),%%ymm6,%%ymm6				\n\t		vaddpd		0x100(%%rdx),%%ymm14,%%ymm14		\n\t"\
		"vaddpd		0x020(%%rdx),%%ymm7,%%ymm7				\n\t		vaddpd		0x120(%%rdx),%%ymm15,%%ymm15		\n\t"\
		"vsubpd		%%ymm6,%%ymm0,%%ymm0				\n\t		vsubpd		%%ymm14,%%ymm8,%%ymm8			\n\t"\
		"vsubpd		%%ymm5,%%ymm2,%%ymm2				\n\t		vsubpd		%%ymm15,%%ymm9,%%ymm9			\n\t"\
		"vsubpd		%%ymm7,%%ymm1,%%ymm1				\n\t"\
		"vsubpd		%%ymm4,%%ymm3,%%ymm3				\n\t		vsubpd		%%ymm13,%%ymm10,%%ymm10		\n\t"\
		"/*vmovaps		%%ymm0,0x080(%%rdx)	*/		\n\t"\
		"/*vmovaps		%%ymm2,0x040(%%rdx)	*/		\n\t		vsubpd		%%ymm12,%%ymm11,%%ymm11		\n\t"\
		"/*vmovaps		%%ymm1,0x0a0(%%rdx)	*/		\n\t		vmulpd		(%%rdi),%%ymm14,%%ymm14		\n\t"\
		"/*vmovaps		%%ymm3,0x0e0(%%rdx)	*/		\n\t		vmulpd		(%%rdi),%%ymm13,%%ymm13		\n\t"\
		"vmulpd		(%%rdi),%%ymm6,%%ymm6				\n\t		vmulpd		(%%rdi),%%ymm15,%%ymm15		\n\t"\
		"vmulpd		(%%rdi),%%ymm5,%%ymm5				\n\t		vmulpd		(%%rdi),%%ymm12,%%ymm12		\n\t"\
		"vmulpd		(%%rdi),%%ymm7,%%ymm7				\n\t		vaddpd		%%ymm8,%%ymm14,%%ymm14		\n\t"\
		"vmulpd		(%%rdi),%%ymm4,%%ymm4				\n\t		vaddpd		%%ymm10,%%ymm13,%%ymm13		\n\t"\
		"vaddpd		%%ymm0,%%ymm6,%%ymm6				\n\t		vaddpd		%%ymm9,%%ymm15,%%ymm15		\n\t"\
		"vaddpd		%%ymm2,%%ymm5,%%ymm5				\n\t		vaddpd		%%ymm11,%%ymm12,%%ymm12		\n\t"\
		"vaddpd		%%ymm1,%%ymm7,%%ymm7				\n\t		vmovaps		%%ymm14,0x100(%%rdx)				\n\t"\
		"vaddpd		%%ymm3,%%ymm4,%%ymm4				\n\t		vmovaps		%%ymm15,0x120(%%rdx)				\n\t"\
		"/*vmovaps		%%ymm6,     (%%rdx)	*/		\n\t		vmovaps		%%ymm10,%%ymm14						\n\t"\
		"/*vmovaps		%%ymm5,0x0c0(%%rdx)	*/		\n\t		vmovaps		%%ymm13,%%ymm15						\n\t"\
		"/*vmovaps		%%ymm7,0x020(%%rdx)	*/		\n\t		vsubpd		%%ymm12,%%ymm10,%%ymm10		\n\t"\
		"/*vmovaps		%%ymm4,0x060(%%rdx)	*/		\n\t		vsubpd		%%ymm11,%%ymm13,%%ymm13		\n\t"\
		"													\n\t		vaddpd		%%ymm12,%%ymm14,%%ymm14		\n\t"\
		"													\n\t		vaddpd		%%ymm11,%%ymm15,%%ymm15		\n\t"\
		"													\n\t		vmulpd		(%%rsi),%%ymm10,%%ymm10		\n\t"\
		"													\n\t		vmulpd		(%%rsi),%%ymm13,%%ymm13		\n\t"\
		"													\n\t		vmulpd		(%%rsi),%%ymm14,%%ymm14		\n\t"\
		"													\n\t		vmulpd		(%%rsi),%%ymm15,%%ymm15		\n\t"\
		"/*****	SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r20) *****/\n\t"\
		"																vmovaps		0x100(%%rdx),%%ymm11		\n\t"\
		"																vmovaps		0x120(%%rdx),%%ymm12		\n\t"\
		"vsubpd		%%ymm11,%%ymm6,%%ymm6				\n\t		vsubpd		%%ymm10,%%ymm2,%%ymm2			\n\t"\
		"vsubpd		%%ymm9,%%ymm0,%%ymm0				\n\t		vsubpd		%%ymm15,%%ymm5,%%ymm5			\n\t"\
		"vsubpd		%%ymm12,%%ymm7,%%ymm7				\n\t		vsubpd		%%ymm14,%%ymm4,%%ymm4			\n\t"\
		"vsubpd		%%ymm8,%%ymm1,%%ymm1				\n\t		vsubpd		%%ymm13,%%ymm3,%%ymm3			\n\t"\
		"vmulpd		(%%rdi),%%ymm11,%%ymm11			\n\t		vmulpd		(%%rdi),%%ymm10,%%ymm10		\n\t"\
		"vmulpd		(%%rdi),%%ymm9,%%ymm9				\n\t		vmulpd		(%%rdi),%%ymm15,%%ymm15		\n\t"\
		"vmulpd		(%%rdi),%%ymm12,%%ymm12			\n\t		vmulpd		(%%rdi),%%ymm14,%%ymm14		\n\t"\
		"vmulpd		(%%rdi),%%ymm8,%%ymm8				\n\t		vmulpd		(%%rdi),%%ymm13,%%ymm13		\n\t"\
		"vaddpd		%%ymm6,%%ymm11,%%ymm11			\n\t		vaddpd		%%ymm2,%%ymm10,%%ymm10		\n\t"\
		"vaddpd		%%ymm0,%%ymm9,%%ymm9				\n\t		vaddpd		%%ymm5,%%ymm15,%%ymm15		\n\t"\
		"vaddpd		%%ymm7,%%ymm12,%%ymm12			\n\t		vaddpd		%%ymm4,%%ymm14,%%ymm14		\n\t"\
		"vaddpd		%%ymm1,%%ymm8,%%ymm8				\n\t		vaddpd		%%ymm3,%%ymm13,%%ymm13		\n\t"\
		"vmovaps		%%ymm6,0x100(%%rdx)			\n\t		vmovaps		%%ymm2,0x140(%%rdx)	\n\t"\
		"vmovaps		%%ymm0,0x080(%%rdx)			\n\t		vmovaps		%%ymm5,0x0c0(%%rdx)	\n\t"\
		"vmovaps		%%ymm7,0x120(%%rdx)			\n\t		vmovaps		%%ymm4,0x160(%%rdx)	\n\t"\
		"vmovaps		%%ymm1,0x1a0(%%rdx)			\n\t		vmovaps		%%ymm3,0x1e0(%%rdx)	\n\t"\
		"vmovaps		%%ymm11,     (%%rdx)			\n\t		vmovaps		%%ymm10,0x040(%%rdx)	\n\t"\
		"vmovaps		%%ymm9,0x180(%%rdx)			\n\t		vmovaps		%%ymm15,0x1c0(%%rdx)	\n\t"\
		"vmovaps		%%ymm12,0x020(%%rdx)			\n\t		vmovaps		%%ymm14,0x060(%%rdx)	\n\t"\
		"vmovaps		%%ymm8,0x0a0(%%rdx)			\n\t		vmovaps		%%ymm13,0x0e0(%%rdx)	\n\t"\
		"/*...Block 4:	*/\n\t"\
		"/*****	SSE2_RADIX4_DIF_4TWIDDLE         (r30,r34,r32,r36,r30,c03)	*****/\n\t"\
		"movq		%[__c03],%%rbx					\n\t"\
		"movq		%[__r30],%%rax					\n\t"\
		"movq		%%rax,%%rcx		/*****	SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO(r38,r3C,r3A,r3E,r38,c07)	*****/\n\t"\
		"addq		$0x040,%%rcx					\n\t"\
		"vmovaps		 (%%rax),%%ymm0	\n\t	movq %%rax,%%rdx \n\t	vmovaps		0x100(%%rax),%%ymm8			\n\t"\
		"vmovaps		 (%%rcx),%%ymm4						\n\t		vmovaps		0x100(%%rcx),%%ymm12		\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1						\n\t		vmovaps		0x120(%%rax),%%ymm9			\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5						\n\t		vmovaps		0x120(%%rcx),%%ymm13		\n\t"\
		"vmovaps		 (%%rbx),%%ymm6						\n\t		vmovaps		0x100(%%rbx),%%ymm14		\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm7						\n\t		vmovaps		0x120(%%rbx),%%ymm15		\n\t"\
		"vmovaps	%%ymm0,%%ymm2						\n\t		vmovaps		%%ymm8,%%ymm10		\n\t"\
		"vmovaps	%%ymm1,%%ymm3						\n\t		vmovaps		%%ymm9,%%ymm11		\n\t"\
		"vmulpd		%%ymm6,%%ymm0,%%ymm0				\n\t		vmulpd		%%ymm14,%%ymm8,%%ymm8			\n\t"\
		"vmulpd		%%ymm6,%%ymm1,%%ymm1				\n\t		vmulpd		%%ymm14,%%ymm9,%%ymm9			\n\t"\
		"vmulpd		%%ymm7,%%ymm2,%%ymm2				\n\t		vmulpd		%%ymm15,%%ymm10,%%ymm10		\n\t"\
		"vmulpd		%%ymm7,%%ymm3,%%ymm3				\n\t		vmulpd		%%ymm15,%%ymm11,%%ymm11		\n\t"\
		"vmovaps		%%ymm4,%%ymm6					\n\t		vmovaps		%%ymm12,%%ymm14				\n\t"\
		"vaddpd		%%ymm2,%%ymm1,%%ymm1				\n\t		vaddpd		%%ymm10,%%ymm9,%%ymm9			\n\t"\
		"vmovaps		%%ymm5,%%ymm7					\n\t		vmovaps		%%ymm13,%%ymm15				\n\t"\
		"vmulpd		0x040(%%rbx),%%ymm4,%%ymm4				\n\t		vmulpd		0x140(%%rbx),%%ymm12,%%ymm12		\n\t"\
		"vsubpd		%%ymm3,%%ymm0,%%ymm0				\n\t		vsubpd		%%ymm11,%%ymm8,%%ymm8			\n\t"\
		"vmulpd		0x040(%%rbx),%%ymm5,%%ymm5				\n\t		vmulpd		0x140(%%rbx),%%ymm13,%%ymm13		\n\t"\
		"vmulpd		0x060(%%rbx),%%ymm6,%%ymm6				\n\t		vmulpd		0x160(%%rbx),%%ymm14,%%ymm14		\n\t"\
		"vmovaps		%%ymm0,%%ymm2					\n\t		vmovaps		%%ymm8,%%ymm10				\n\t"\
		"vmulpd		0x060(%%rbx),%%ymm7,%%ymm7				\n\t		vmulpd		0x160(%%rbx),%%ymm15,%%ymm15		\n\t"\
		"vaddpd		%%ymm6,%%ymm5,%%ymm5				\n\t		vaddpd		%%ymm14,%%ymm13,%%ymm13		\n\t"\
		"vmovaps		%%ymm1,%%ymm3					\n\t		vmovaps		%%ymm9,%%ymm11				\n\t"\
		"vsubpd		%%ymm7,%%ymm4,%%ymm4				\n\t		vsubpd		%%ymm15,%%ymm12,%%ymm12		\n\t"\
		"addq		$0x080,%%rcx						\n\t		vaddpd		%%ymm12,%%ymm8,%%ymm8			\n\t"\
		"addq		$0x0c0,%%rbx						\n\t		vaddpd		%%ymm13,%%ymm9,%%ymm9			\n\t"\
		"vmovaps			 (%%rcx),%%ymm6					\n\t		vsubpd		%%ymm12,%%ymm10,%%ymm10		\n\t"\
		"vmovaps		0x020(%%rcx),%%ymm7					\n\t		vsubpd		%%ymm13,%%ymm11,%%ymm11		\n\t"\
		"vaddpd		%%ymm4,%%ymm0,%%ymm0				\n\t		vmovaps		0x100(%%rcx),%%ymm12				\n\t"\
		"vaddpd		%%ymm5,%%ymm1,%%ymm1				\n\t		vmovaps		0x120(%%rcx),%%ymm13				\n\t"\
		"vsubpd		%%ymm4,%%ymm2,%%ymm2				\n\t		vmovaps		0x100(%%rcx),%%ymm14				\n\t"\
		"vsubpd		%%ymm5,%%ymm3,%%ymm3				\n\t		vmovaps		0x120(%%rcx),%%ymm15				\n\t"\
		"vmovaps		%%ymm6,%%ymm4					\n\t		vmulpd		0x100(%%rbx),%%ymm12,%%ymm12		\n\t"\
		"vmovaps		%%ymm7,%%ymm5					\n\t		vmulpd		0x100(%%rbx),%%ymm13,%%ymm13		\n\t"\
		"vmulpd			 (%%rbx),%%ymm4,%%ymm4				\n\t		vmulpd		0x120(%%rbx),%%ymm14,%%ymm14		\n\t"\
		"vmulpd			 (%%rbx),%%ymm5,%%ymm5				\n\t		vmulpd		0x120(%%rbx),%%ymm15,%%ymm15		\n\t"\
		"vmulpd		0x020(%%rbx),%%ymm6,%%ymm6				\n\t		vaddpd		%%ymm14,%%ymm13,%%ymm13		\n\t"\
		"vmulpd		0x020(%%rbx),%%ymm7,%%ymm7				\n\t		vsubpd		%%ymm15,%%ymm12,%%ymm12		\n\t"\
		"/*movq		%%rax,%%rdx,%%rdx		*/		\n\t		vmovaps		%%ymm13,0x120(%%rdx)				\n\t"\
		"vaddpd		%%ymm6,%%ymm5,%%ymm5				\n\t		vmovaps		%%ymm12,0x100(%%rdx)				\n\t"\
		"vsubpd		%%ymm7,%%ymm4,%%ymm4				\n\t"\
		"vmovaps		%%ymm5,0x020(%%rdx)			\n\t		addq	$0x080,%%rax							\n\t"\
		"vmovaps		%%ymm4,     (%%rdx)			\n\t		subq	$0x040,%%rbx							\n\t"\
		"vmovaps			 (%%rax),%%ymm4					\n\t		vmovaps		0x100(%%rax),%%ymm12				\n\t"\
		"vmovaps		0x020(%%rax),%%ymm5					\n\t		vmovaps		0x120(%%rax),%%ymm13				\n\t"\
		"vmovaps		%%ymm4,%%ymm6					\n\t		vmovaps		0x100(%%rax),%%ymm14				\n\t"\
		"vmovaps		%%ymm5,%%ymm7					\n\t		vmovaps		0x120(%%rax),%%ymm15				\n\t"\
		"vmulpd			 (%%rbx),%%ymm4,%%ymm4				\n\t		vmulpd		0x100(%%rbx),%%ymm12,%%ymm12		\n\t"\
		"vmulpd			 (%%rbx),%%ymm5,%%ymm5				\n\t		vmulpd		0x100(%%rbx),%%ymm13,%%ymm13		\n\t"\
		"vmulpd		0x020(%%rbx),%%ymm6,%%ymm6				\n\t		vmulpd		0x120(%%rbx),%%ymm14,%%ymm14		\n\t"\
		"vmulpd		0x020(%%rbx),%%ymm7,%%ymm7				\n\t		vmulpd		0x120(%%rbx),%%ymm15,%%ymm15		\n\t"\
		"vaddpd		%%ymm6,%%ymm5,%%ymm5				\n\t		vaddpd		%%ymm14,%%ymm13,%%ymm13		\n\t"\
		"vsubpd		%%ymm7,%%ymm4,%%ymm4				\n\t		vsubpd		%%ymm15,%%ymm12,%%ymm12		\n\t"\
		"vmovaps		%%ymm5,%%ymm7					\n\t		vmovaps		%%ymm13,%%ymm15				\n\t"\
		"vmovaps		%%ymm4,%%ymm6					\n\t		vmovaps		%%ymm12,%%ymm14				\n\t"\
		"vsubpd			 (%%rdx),%%ymm4,%%ymm4				\n\t		vsubpd		0x100(%%rdx),%%ymm12,%%ymm12		\n\t"\
		"vsubpd		0x020(%%rdx),%%ymm5,%%ymm5				\n\t		vsubpd		0x120(%%rdx),%%ymm13,%%ymm13		\n\t"\
		"vaddpd			 (%%rdx),%%ymm6,%%ymm6				\n\t		vaddpd		0x100(%%rdx),%%ymm14,%%ymm14		\n\t"\
		"vaddpd		0x020(%%rdx),%%ymm7,%%ymm7				\n\t		vaddpd		0x120(%%rdx),%%ymm15,%%ymm15		\n\t"\
		"vsubpd		%%ymm6,%%ymm0,%%ymm0				\n\t		vsubpd		%%ymm14,%%ymm8,%%ymm8			\n\t"\
		"vsubpd		%%ymm5,%%ymm2,%%ymm2				\n\t		vsubpd		%%ymm15,%%ymm9,%%ymm9			\n\t"\
		"vsubpd		%%ymm7,%%ymm1,%%ymm1				\n\t"\
		"vsubpd		%%ymm4,%%ymm3,%%ymm3				\n\t		vsubpd		%%ymm13,%%ymm10,%%ymm10		\n\t"\
		"/*vmovaps		%%ymm0,0x080(%%rdx)	*/		\n\t"\
		"/*vmovaps		%%ymm2,0x040(%%rdx)	*/		\n\t		vsubpd		%%ymm12,%%ymm11,%%ymm11		\n\t"\
		"/*vmovaps		%%ymm1,0x0a0(%%rdx)	*/		\n\t		vmulpd		(%%rdi),%%ymm14,%%ymm14		\n\t"\
		"/*vmovaps		%%ymm3,0x0e0(%%rdx)	*/		\n\t		vmulpd		(%%rdi),%%ymm13,%%ymm13		\n\t"\
		"vmulpd		(%%rdi),%%ymm6,%%ymm6				\n\t		vmulpd		(%%rdi),%%ymm15,%%ymm15		\n\t"\
		"vmulpd		(%%rdi),%%ymm5,%%ymm5				\n\t		vmulpd		(%%rdi),%%ymm12,%%ymm12		\n\t"\
		"vmulpd		(%%rdi),%%ymm7,%%ymm7				\n\t		vaddpd		%%ymm8,%%ymm14,%%ymm14		\n\t"\
		"vmulpd		(%%rdi),%%ymm4,%%ymm4				\n\t		vaddpd		%%ymm10,%%ymm13,%%ymm13		\n\t"\
		"vaddpd		%%ymm0,%%ymm6,%%ymm6				\n\t		vaddpd		%%ymm9,%%ymm15,%%ymm15		\n\t"\
		"vaddpd		%%ymm2,%%ymm5,%%ymm5				\n\t		vaddpd		%%ymm11,%%ymm12,%%ymm12		\n\t"\
		"vaddpd		%%ymm1,%%ymm7,%%ymm7				\n\t		vmovaps		%%ymm14,0x100(%%rdx)				\n\t"\
		"vaddpd		%%ymm3,%%ymm4,%%ymm4				\n\t		vmovaps		%%ymm15,0x120(%%rdx)				\n\t"\
		"/*vmovaps		%%ymm6,     (%%rdx)	*/		\n\t		vmovaps		%%ymm10,%%ymm14						\n\t"\
		"/*vmovaps		%%ymm5,0x0c0(%%rdx)	*/		\n\t		vmovaps		%%ymm13,%%ymm15						\n\t"\
		"/*vmovaps		%%ymm7,0x020(%%rdx)	*/		\n\t		vsubpd		%%ymm12,%%ymm10,%%ymm10		\n\t"\
		"/*vmovaps		%%ymm4,0x060(%%rdx)	*/		\n\t		vsubpd		%%ymm11,%%ymm13,%%ymm13		\n\t"\
		"													\n\t		vaddpd		%%ymm12,%%ymm14,%%ymm14		\n\t"\
		"													\n\t		vaddpd		%%ymm11,%%ymm15,%%ymm15		\n\t"\
		"													\n\t		vmulpd		(%%rsi),%%ymm10,%%ymm10		\n\t"\
		"													\n\t		vmulpd		(%%rsi),%%ymm13,%%ymm13		\n\t"\
		"													\n\t		vmulpd		(%%rsi),%%ymm14,%%ymm14		\n\t"\
		"													\n\t		vmulpd		(%%rsi),%%ymm15,%%ymm15		\n\t"\
		"/*****	SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r20) *****/\n\t"\
		"																vmovaps		0x100(%%rdx),%%ymm11		\n\t"\
		"																vmovaps		0x120(%%rdx),%%ymm12		\n\t"\
		"vsubpd		%%ymm11,%%ymm6,%%ymm6				\n\t		vsubpd		%%ymm10,%%ymm2,%%ymm2			\n\t"\
		"vsubpd		%%ymm9,%%ymm0,%%ymm0				\n\t		vsubpd		%%ymm15,%%ymm5,%%ymm5			\n\t"\
		"vsubpd		%%ymm12,%%ymm7,%%ymm7				\n\t		vsubpd		%%ymm14,%%ymm4,%%ymm4			\n\t"\
		"vsubpd		%%ymm8,%%ymm1,%%ymm1				\n\t		vsubpd		%%ymm13,%%ymm3,%%ymm3			\n\t"\
		"vmulpd		(%%rdi),%%ymm11,%%ymm11			\n\t		vmulpd		(%%rdi),%%ymm10,%%ymm10		\n\t"\
		"vmulpd		(%%rdi),%%ymm9,%%ymm9				\n\t		vmulpd		(%%rdi),%%ymm15,%%ymm15		\n\t"\
		"vmulpd		(%%rdi),%%ymm12,%%ymm12			\n\t		vmulpd		(%%rdi),%%ymm14,%%ymm14		\n\t"\
		"vmulpd		(%%rdi),%%ymm8,%%ymm8				\n\t		vmulpd		(%%rdi),%%ymm13,%%ymm13		\n\t"\
		"vaddpd		%%ymm6,%%ymm11,%%ymm11			\n\t		vaddpd		%%ymm2,%%ymm10,%%ymm10		\n\t"\
		"vaddpd		%%ymm0,%%ymm9,%%ymm9				\n\t		vaddpd		%%ymm5,%%ymm15,%%ymm15		\n\t"\
		"vaddpd		%%ymm7,%%ymm12,%%ymm12			\n\t		vaddpd		%%ymm4,%%ymm14,%%ymm14		\n\t"\
		"vaddpd		%%ymm1,%%ymm8,%%ymm8				\n\t		vaddpd		%%ymm3,%%ymm13,%%ymm13		\n\t"\
		"vmovaps		%%ymm6,0x100(%%rdx)			\n\t		vmovaps		%%ymm2,0x140(%%rdx)	\n\t"\
		"vmovaps		%%ymm0,0x080(%%rdx)			\n\t		vmovaps		%%ymm5,0x0c0(%%rdx)	\n\t"\
		"vmovaps		%%ymm7,0x120(%%rdx)			\n\t		vmovaps		%%ymm4,0x160(%%rdx)	\n\t"\
		"vmovaps		%%ymm1,0x1a0(%%rdx)			\n\t		vmovaps		%%ymm3,0x1e0(%%rdx)	\n\t"\
		"vmovaps		%%ymm11,     (%%rdx)			\n\t		vmovaps		%%ymm10,0x040(%%rdx)	\n\t"\
		"vmovaps		%%ymm9,0x180(%%rdx)			\n\t		vmovaps		%%ymm15,0x1c0(%%rdx)	\n\t"\
		"vmovaps		%%ymm12,0x020(%%rdx)			\n\t		vmovaps		%%ymm14,0x060(%%rdx)	\n\t"\
		"vmovaps		%%ymm8,0x0a0(%%rdx)			\n\t		vmovaps		%%ymm13,0x0e0(%%rdx)	\n\t"\
	/**********************************************************************************/\
	/*...and now do eight radix-4 transforms, including the internal twiddle factors: */\
	/**********************************************************************************/\
		"movq		%[__isrt2],%%rsi		\n\t"\
		"/*...Block 1: t00,t10,t20,t30	*/					\n\t		/*...Block 5: t08,t18,t28,t38	*/		\n\t"\
		"movq		%[__r00],%%rax					\n\t		movq		$0x100,%%r10			\n\t"\
		"movq		%[__r10],%%rbx					\n\t		movq		$0x100,%%r11			\n\t"\
		"movq		%[__r20],%%rcx					\n\t		movq		$0x100,%%r12			\n\t"\
		"movq		%[__r30],%%rdx					\n\t		movq		$0x100,%%r13			\n\t"\
		"vmovaps			  (%%rax),%%ymm0				\n\t		addq		%%rax ,%%r10			\n\t"\
		"vmovaps		 0x020(%%rax),%%ymm1				\n\t		addq		%%rbx ,%%r11			\n\t"\
		"vmovaps			  (%%rbx),%%ymm2				\n\t		addq		%%rcx ,%%r12			\n\t"\
		"vmovaps		 0x020(%%rbx),%%ymm3				\n\t		addq		%%rdx ,%%r13			\n\t"\
		"vaddpd			   %%ymm0,%%ymm2,%%ymm2			\n\t		vmovaps			(%%r10),%%ymm8			\n\t"\
		"vaddpd			   %%ymm1,%%ymm3,%%ymm3			\n\t		vmovaps		0x20(%%r10),%%ymm9			\n\t"\
		"vsubpd			  (%%rbx),%%ymm0,%%ymm0			\n\t		vmovaps			(%%r11),%%ymm10			\n\t"\
		"vsubpd		 0x020(%%rbx),%%ymm1,%%ymm1			\n\t		vmovaps		0x20(%%r11),%%ymm11			\n\t"\
		"vmovaps			  (%%rcx),%%ymm4				\n\t		vaddpd			 %%ymm9 ,%%ymm10,%%ymm10	\n\t"\
		"vmovaps		 0x020(%%rcx),%%ymm5				\n\t		vaddpd			 %%ymm8 ,%%ymm11,%%ymm11	\n\t"\
		"vmovaps			  (%%rdx),%%ymm6				\n\t		vsubpd		0x020(%%r11),%%ymm8,%%ymm8		\n\t"\
		"vmovaps		 0x020(%%rdx),%%ymm7				\n\t		vsubpd			 (%%r11),%%ymm9,%%ymm9		\n\t"\
		"vaddpd			   %%ymm4,%%ymm6,%%ymm6			\n\t		vmovaps			(%%r12),%%ymm12					\n\t"\
		"vaddpd			   %%ymm5,%%ymm7,%%ymm7			\n\t		vmovaps		0x20(%%r12),%%ymm13					\n\t"\
		"vsubpd			  (%%rdx),%%ymm4,%%ymm4			\n\t		vmovaps			(%%r13),%%ymm14					\n\t"\
		"vsubpd		 0x020(%%rdx),%%ymm5,%%ymm5			\n\t		vmovaps		0x20(%%r13),%%ymm15					\n\t"\
		"vsubpd		%%ymm6,%%ymm2,%%ymm2			\n\t		vsubpd			 %%ymm13,%%ymm12,%%ymm12	\n\t"\
		"vsubpd		%%ymm7,%%ymm3,%%ymm3			\n\t		vaddpd			(%%r12),%%ymm13,%%ymm13	\n\t"\
		"vmulpd		(%%rdi),%%ymm6,%%ymm6			\n\t		vaddpd			 %%ymm15,%%ymm14,%%ymm14	\n\t"\
		"vmulpd		(%%rdi),%%ymm7,%%ymm7			\n\t		vsubpd			(%%r13),%%ymm15,%%ymm15	\n\t"\
		"vmovaps		%%ymm2,      (%%rcx)			\n\t		vmulpd		(%%rsi),%%ymm12,%%ymm12		\n\t"\
		"vmovaps		%%ymm3, 0x020(%%rcx)			\n\t		vmulpd		(%%rsi),%%ymm13,%%ymm13		\n\t"\
		"vaddpd		%%ymm2,%%ymm6,%%ymm6			\n\t		vmulpd		(%%rsi),%%ymm14,%%ymm14		\n\t"\
		"vaddpd		%%ymm3,%%ymm7,%%ymm7			\n\t		vmulpd		(%%rsi),%%ymm15,%%ymm15		\n\t"\
		"vmovaps		%%ymm6,      (%%rax)			\n\t		vsubpd		%%ymm14,%%ymm12,%%ymm12		\n\t"\
		"vmovaps		%%ymm7, 0x020(%%rax)			\n\t		vsubpd		%%ymm15,%%ymm13,%%ymm13		\n\t"\
		"vsubpd		%%ymm5,%%ymm0,%%ymm0			\n\t		vmulpd		(%%rdi),%%ymm14,%%ymm14		\n\t"\
		"vsubpd		%%ymm4,%%ymm1,%%ymm1			\n\t		vmulpd		(%%rdi),%%ymm15,%%ymm15		\n\t"\
		"vmulpd		(%%rdi),%%ymm5,%%ymm5			\n\t		vaddpd		%%ymm12,%%ymm14,%%ymm14		\n\t"\
		"vmulpd		(%%rdi),%%ymm4,%%ymm4			\n\t		vaddpd		%%ymm13,%%ymm15,%%ymm15		\n\t"\
		"vmovaps		%%ymm0,      (%%rbx)			\n\t		vsubpd		%%ymm12,%%ymm8,%%ymm8			\n\t"\
		"vmovaps		%%ymm1, 0x020(%%rdx)			\n\t		vsubpd		%%ymm13,%%ymm10,%%ymm10		\n\t"\
		"vaddpd		%%ymm0,%%ymm5,%%ymm5			\n\t		vmulpd		(%%rdi),%%ymm12,%%ymm12		\n\t"\
		"vaddpd		%%ymm1,%%ymm4,%%ymm4			\n\t		vmulpd		(%%rdi),%%ymm13,%%ymm13		\n\t"\
		"vmovaps		%%ymm5,      (%%rdx)			\n\t		vmovaps		%%ymm8 ,     (%%r12)				\n\t"\
		"vmovaps		%%ymm4, 0x020(%%rbx)			\n\t		vmovaps		%%ymm10,0x020(%%r12)				\n\t"\
		"/*...Block 3: t04,t14,t24,t34	*/					\n\t"\
		"addq		$0x080,%%rax						\n\t		vaddpd		%%ymm8,%%ymm12,%%ymm12		\n\t"\
		"addq		$0x080,%%rbx						\n\t		vaddpd		%%ymm10,%%ymm13,%%ymm13		\n\t"\
		"addq		$0x080,%%rcx						\n\t		vmovaps		%%ymm12,     (%%r10)				\n\t"\
		"addq		$0x080,%%rdx						\n\t		vmovaps		%%ymm13,0x020(%%r10)				\n\t"\
		"vmovaps		0x020(%%rsi),%%ymm8		/* c */	\n\t		vsubpd		%%ymm15,%%ymm11,%%ymm11		\n\t"\
		"vmovaps		0x040(%%rsi),%%ymm10	/* s */	\n\t		vsubpd		%%ymm14,%%ymm9,%%ymm9			\n\t"\
		"vmovaps			 (%%rcx),%%ymm4				\n\t		vmulpd		(%%rdi),%%ymm15,%%ymm15		\n\t"\
		"vmovaps			 (%%rdx),%%ymm6				\n\t		vmulpd		(%%rdi),%%ymm14,%%ymm14		\n\t"\
		"vmovaps		0x020(%%rcx),%%ymm5				\n\t		vmovaps		%%ymm11,     (%%r11)				\n\t"\
		"vmovaps		0x020(%%rdx),%%ymm7				\n\t		vmovaps		%%ymm9 ,0x020(%%r13)				\n\t"\
		"vmovaps			  %%ymm4,%%ymm0				\n\t		vaddpd		%%ymm11,%%ymm15,%%ymm15		\n\t"\
		"vmovaps			  %%ymm6,%%ymm2				\n\t		vaddpd		%%ymm9,%%ymm14,%%ymm14		\n\t"\
		"vmovaps			  %%ymm5,%%ymm1				\n\t		vmovaps		%%ymm15,     (%%r13)				\n\t"\
		"vmovaps			  %%ymm7,%%ymm3				\n\t		vmovaps		%%ymm14,0x020(%%r11)				\n\t"\
		"													/*...Block 7: t0C,t1C,t2C,t3C	*/					\n\t"\
		"vmulpd		%%ymm8,%%ymm4,%%ymm4				\n\t		addq		$0x080,%%r10					\n\t"\
		"vmulpd		%%ymm10,%%ymm6,%%ymm6				\n\t		addq		$0x080,%%r11					\n\t"\
		"vmulpd		%%ymm10,%%ymm1,%%ymm1				\n\t		addq		$0x080,%%r12					\n\t"\
		"vmulpd		%%ymm8,%%ymm3,%%ymm3				\n\t		addq		$0x080,%%r13					\n\t"\
		"vmulpd		%%ymm8,%%ymm5,%%ymm5				\n\t		vmovaps			 (%%r12),%%ymm12				\n\t"\
		"vmulpd		%%ymm10,%%ymm7,%%ymm7				\n\t		vmovaps			 (%%r13),%%ymm14				\n\t"\
		"vmulpd		%%ymm10,%%ymm0,%%ymm0				\n\t		vmovaps		0x020(%%r12),%%ymm13				\n\t"\
		"vmulpd		%%ymm8,%%ymm2,%%ymm2				\n\t		vmovaps		0x020(%%r13),%%ymm15				\n\t"\
		"vsubpd		%%ymm1,%%ymm4,%%ymm4				\n\t		vmovaps			 %%ymm13,%%ymm9					\n\t"\
		"vsubpd		%%ymm3,%%ymm6,%%ymm6				\n\t		vmovaps			 %%ymm15,%%ymm11				\n\t"\
		"vaddpd		%%ymm0,%%ymm5,%%ymm5				\n\t		vmulpd		 %%ymm10,%%ymm12,%%ymm12			\n\t"\
		"vaddpd		%%ymm2,%%ymm7,%%ymm7				\n\t		vmulpd		 %%ymm8,%%ymm14,%%ymm14			\n\t"\
		"vsubpd		%%ymm6,%%ymm4,%%ymm4				\n\t		vmulpd		 %%ymm8,%%ymm9,%%ymm9				\n\t"\
		"vsubpd		%%ymm7,%%ymm5,%%ymm5				\n\t		vmulpd		 %%ymm10,%%ymm11,%%ymm11			\n\t"\
		"vmulpd		(%%rdi),%%ymm6,%%ymm6				\n\t		vmulpd		 %%ymm10,%%ymm13,%%ymm13			\n\t"\
		"vmulpd		(%%rdi),%%ymm7,%%ymm7				\n\t		vmulpd		 %%ymm8,%%ymm15,%%ymm15			\n\t"\
		"vaddpd		%%ymm4,%%ymm6,%%ymm6				\n\t		vmulpd		(%%r12),%%ymm8,%%ymm8				\n\t"\
		"vaddpd		%%ymm5,%%ymm7,%%ymm7				\n\t		vmulpd		(%%r13),%%ymm10,%%ymm10			\n\t"\
		"vmovaps			  (%%rbx),%%ymm2				\n\t		vsubpd		%%ymm9,%%ymm12,%%ymm12		\n\t"\
		"vmovaps		 0x020(%%rbx),%%ymm3				\n\t		vsubpd		%%ymm11,%%ymm14,%%ymm14		\n\t"\
		"vsubpd		0x020(%%rbx),%%ymm2,%%ymm2				\n\t		vaddpd		%%ymm8,%%ymm13,%%ymm13		\n\t"\
		"vaddpd			 (%%rbx),%%ymm3,%%ymm3				\n\t		vaddpd		%%ymm10,%%ymm15,%%ymm15		\n\t"\
		"vmulpd			 (%%rsi),%%ymm2,%%ymm2				\n\t		vsubpd		%%ymm14,%%ymm12,%%ymm12		\n\t"\
		"vmulpd			 (%%rsi),%%ymm3,%%ymm3				\n\t		vsubpd		%%ymm15,%%ymm13,%%ymm13		\n\t"\
		"vmovaps			  (%%rax),%%ymm0				\n\t		vmulpd		(%%rdi),%%ymm14,%%ymm14		\n\t"\
		"vmovaps		 0x020(%%rax),%%ymm1				\n\t		vmulpd		(%%rdi),%%ymm15,%%ymm15		\n\t"\
		"vsubpd		%%ymm2,%%ymm0,%%ymm0				\n\t		vaddpd		%%ymm12,%%ymm14,%%ymm14		\n\t"\
		"vsubpd		%%ymm3,%%ymm1,%%ymm1				\n\t		vaddpd		%%ymm13,%%ymm15,%%ymm15		\n\t"\
		"vaddpd			 (%%rax),%%ymm2,%%ymm2				\n\t		vmovaps			 (%%r11),%%ymm10				\n\t"\
		"vaddpd		0x020(%%rax),%%ymm3,%%ymm3				\n\t		vmovaps		0x020(%%r11),%%ymm11				\n\t"\
		"vsubpd		%%ymm6,%%ymm2,%%ymm2				\n\t		vaddpd		0x020(%%r11),%%ymm10,%%ymm10	\n\t"\
		"vsubpd		%%ymm7,%%ymm3,%%ymm3				\n\t		vsubpd			 (%%r11),%%ymm11,%%ymm11	\n\t"\
		"vmulpd		(%%rdi),%%ymm6,%%ymm6				\n\t		vmulpd			 (%%rsi),%%ymm10,%%ymm10	\n\t"\
		"vmulpd		(%%rdi),%%ymm7,%%ymm7				\n\t		vmulpd			 (%%rsi),%%ymm11,%%ymm11	\n\t"\
		"vmovaps		%%ymm2,      (%%rcx)			\n\t		vmovaps			  (%%r10),%%ymm8				\n\t"\
		"vmovaps		%%ymm3, 0x020(%%rcx)			\n\t		vmovaps		 0x020(%%r10),%%ymm9				\n\t"\
		"vaddpd		%%ymm2,%%ymm6,%%ymm6				\n\t		vsubpd		%%ymm10,%%ymm8,%%ymm8			\n\t"\
		"vaddpd		%%ymm3,%%ymm7,%%ymm7				\n\t		vsubpd		%%ymm11,%%ymm9,%%ymm9			\n\t"\
		"vmovaps		%%ymm6,      (%%rax)			\n\t		vaddpd			 (%%r10),%%ymm10,%%ymm10	\n\t"\
		"vmovaps		%%ymm7, 0x020(%%rax)			\n\t		vaddpd		0x020(%%r10),%%ymm11,%%ymm11	\n\t"\
		"vsubpd		%%ymm5,%%ymm0,%%ymm0				\n\t		vsubpd		%%ymm12,%%ymm8,%%ymm8			\n\t"\
		"vsubpd		%%ymm4,%%ymm1,%%ymm1				\n\t		vsubpd		%%ymm13,%%ymm9,%%ymm9			\n\t"\
		"vmulpd		(%%rdi),%%ymm5,%%ymm5				\n\t		vmulpd		(%%rdi),%%ymm12,%%ymm12		\n\t"\
		"vmulpd		(%%rdi),%%ymm4,%%ymm4				\n\t		vmulpd		(%%rdi),%%ymm13,%%ymm13		\n\t"\
		"vmovaps		%%ymm0,      (%%rbx)			\n\t		vmovaps		%%ymm8 ,     (%%r12)				\n\t"\
		"vmovaps		%%ymm1, 0x020(%%rdx)			\n\t		vmovaps		%%ymm9 ,0x020(%%r12)				\n\t"\
		"vaddpd		%%ymm0,%%ymm5,%%ymm5				\n\t		vaddpd		%%ymm8,%%ymm12,%%ymm12		\n\t"\
		"vaddpd		%%ymm1,%%ymm4,%%ymm4				\n\t		vaddpd		%%ymm9,%%ymm13,%%ymm13		\n\t"\
		"vmovaps	%%ymm5,      (%%rdx)				\n\t		vmovaps		%%ymm12,     (%%r10)				\n\t"\
		"vmovaps	%%ymm4, 0x020(%%rbx)				\n\t		vmovaps		%%ymm13,0x020(%%r10)				\n\t"\
		"/*...Block 2: t02,t12,t22,t32	*/\n\t"\
		"subq		$0x040,%%rax						\n\t		vsubpd		%%ymm15,%%ymm10,%%ymm10		\n\t"\
		"subq		$0x040,%%rbx						\n\t		vsubpd		%%ymm14,%%ymm11,%%ymm11		\n\t"\
		"subq		$0x040,%%rcx						\n\t		vmulpd		(%%rdi),%%ymm15,%%ymm15		\n\t"\
		"subq		$0x040,%%rdx						\n\t		vmulpd		(%%rdi),%%ymm14,%%ymm14		\n\t"\
		"addq		$0x060,%%rsi	/* cc1 */			\n\t		vmovaps		%%ymm10,     (%%r11)				\n\t"\
		"vmovaps		 (%%rcx),%%ymm4					\n\t		vmovaps		%%ymm11,0x020(%%r13)				\n\t"\
		"vmovaps		 (%%rdx),%%ymm6					\n\t		vaddpd		%%ymm10,%%ymm15,%%ymm15		\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5					\n\t		vaddpd		%%ymm11,%%ymm14,%%ymm14		\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm7					\n\t		vmovaps		%%ymm15,     (%%r13)		\n\t"\
		"vmovaps		 (%%rcx),%%ymm0					\n\t		vmovaps		%%ymm14,0x020(%%r11)		\n\t"\
		"																/*...Block 6: t0A,t1A,t2A,t3A	*/\n\t"\
		"vmovaps		 (%%rdx),%%ymm2					\n\t		subq		$0x040,%%r10			\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm1					\n\t		subq		$0x040,%%r11			\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm3					\n\t		subq		$0x040,%%r12			\n\t"\
		"vmulpd			 (%%rsi),%%ymm4,%%ymm4			\n\t		subq		$0x040,%%r13			\n\t"\
		"vmulpd		0x040(%%rsi),%%ymm6,%%ymm6			\n\t		vmovaps			  (%%r12),%%ymm12	\n\t"\
		"vmulpd		0x020(%%rsi),%%ymm1,%%ymm1			\n\t		vmovaps			  (%%r13),%%ymm14	\n\t"\
		"vmulpd		0x060(%%rsi),%%ymm3,%%ymm3			\n\t		vmovaps		 0x020(%%r12),%%ymm13	\n\t"\
		"vmulpd			 (%%rsi),%%ymm5,%%ymm5			\n\t		vmovaps		 0x020(%%r13),%%ymm15	\n\t"\
		"vmulpd		0x040(%%rsi),%%ymm7,%%ymm7			\n\t		vmovaps			  (%%r12),%%ymm8		\n\t"\
		"vmulpd		0x020(%%rsi),%%ymm0,%%ymm0			\n\t		vmovaps			  (%%r13),%%ymm10	\n\t"\
		"vmulpd		0x060(%%rsi),%%ymm2,%%ymm2			\n\t		vmovaps		 0x020(%%r12),%%ymm9		\n\t"\
		"vsubpd		%%ymm1,%%ymm4,%%ymm4			\n\t		vmovaps		 0x020(%%r13),%%ymm11	\n\t"\
		"vsubpd		%%ymm3,%%ymm6,%%ymm6			\n\t		vmulpd		0x060(%%rsi),%%ymm12,%%ymm12	\n\t"\
		"vaddpd		%%ymm0,%%ymm5,%%ymm5			\n\t		vmulpd			 (%%rsi),%%ymm14,%%ymm14	\n\t"\
		"vaddpd		%%ymm2,%%ymm7,%%ymm7			\n\t		vmulpd		0x040(%%rsi),%%ymm9,%%ymm9		\n\t"\
		"vsubpd		%%ymm6,%%ymm4,%%ymm4			\n\t		vmulpd		0x020(%%rsi),%%ymm11,%%ymm11	\n\t"\
		"vsubpd		%%ymm7,%%ymm5,%%ymm5			\n\t		vmulpd		0x060(%%rsi),%%ymm13,%%ymm13	\n\t"\
		"vmulpd		(%%rdi),%%ymm6,%%ymm6			\n\t		vmulpd			 (%%rsi),%%ymm15,%%ymm15	\n\t"\
		"vmulpd		(%%rdi),%%ymm7,%%ymm7			\n\t		vmulpd		0x040(%%rsi),%%ymm8,%%ymm8		\n\t"\
		"vaddpd		%%ymm4,%%ymm6,%%ymm6			\n\t		vmulpd		0x020(%%rsi),%%ymm10,%%ymm10	\n\t"\
		"vaddpd		%%ymm5,%%ymm7,%%ymm7			\n\t		vsubpd		%%ymm9,%%ymm12,%%ymm12		\n\t"\
		"vmovaps			  (%%rbx),%%ymm2				\n\t		vaddpd		%%ymm11,%%ymm14,%%ymm14		\n\t"\
		"vmovaps		 0x020(%%rbx),%%ymm0				\n\t		vaddpd		%%ymm8,%%ymm13,%%ymm13		\n\t"\
		"vmovaps			  (%%rbx),%%ymm1				\n\t		vsubpd		%%ymm10,%%ymm15,%%ymm15		\n\t"\
		"vmovaps		 0x020(%%rbx),%%ymm3				\n\t		vsubpd		%%ymm14,%%ymm12,%%ymm12		\n\t"\
		"vmulpd		-0x040(%%rsi),%%ymm2,%%ymm2			\n\t		vsubpd		%%ymm15,%%ymm13,%%ymm13		\n\t"\
		"vmulpd		-0x020(%%rsi),%%ymm0,%%ymm0			\n\t		vmulpd		(%%rdi),%%ymm14,%%ymm14		\n\t"\
		"vmulpd		-0x040(%%rsi),%%ymm3,%%ymm3			\n\t		vmulpd		(%%rdi),%%ymm15,%%ymm15		\n\t"\
		"vmulpd		-0x020(%%rsi),%%ymm1,%%ymm1			\n\t		vaddpd		%%ymm12,%%ymm14,%%ymm14		\n\t"\
		"vsubpd		%%ymm0,%%ymm2,%%ymm2			\n\t		vaddpd		%%ymm13,%%ymm15,%%ymm15		\n\t"\
		"vaddpd		%%ymm1,%%ymm3,%%ymm3			\n\t		vmovaps			 (%%r11),%%ymm10				\n\t"\
		"vmovaps			  (%%rax),%%ymm0				\n\t		vmovaps		0x020(%%r11),%%ymm8					\n\t"\
		"vmovaps		 0x020(%%rax),%%ymm1				\n\t		vmovaps			 (%%r11),%%ymm9					\n\t"\
		"vsubpd		%%ymm2,%%ymm0,%%ymm0			\n\t		vmovaps		0x020(%%r11),%%ymm11				\n\t"\
		"vsubpd		%%ymm3,%%ymm1,%%ymm1			\n\t		vmulpd		-0x20(%%rsi),%%ymm10,%%ymm10	\n\t"\
		"vaddpd			  (%%rax),%%ymm2,%%ymm2			\n\t		vmulpd		-0x40(%%rsi),%%ymm8,%%ymm8		\n\t"\
		"vaddpd		 0x020(%%rax),%%ymm3,%%ymm3			\n\t		vmulpd		-0x20(%%rsi),%%ymm11,%%ymm11	\n\t"\
		"vsubpd		%%ymm6,%%ymm2,%%ymm2			\n\t		vmulpd		-0x40(%%rsi),%%ymm9,%%ymm9		\n\t"\
		"vsubpd		%%ymm7,%%ymm3,%%ymm3			\n\t		vaddpd		%%ymm8,%%ymm10,%%ymm10		\n\t"\
		"vmulpd		(%%rdi),%%ymm6,%%ymm6			\n\t		vsubpd		%%ymm9,%%ymm11,%%ymm11		\n\t"\
		"vmulpd		(%%rdi),%%ymm7,%%ymm7			\n\t		vmovaps			 (%%r10),%%ymm8					\n\t"\
		"vmovaps		%%ymm2,      (%%rcx)			\n\t		vmovaps		0x020(%%r10),%%ymm9					\n\t"\
		"vmovaps		%%ymm3, 0x020(%%rcx)			\n\t		vsubpd		%%ymm10,%%ymm8,%%ymm8			\n\t"\
		"vaddpd		%%ymm2,%%ymm6,%%ymm6			\n\t		vsubpd		%%ymm11,%%ymm9,%%ymm9			\n\t"\
		"vaddpd		%%ymm3,%%ymm7,%%ymm7			\n\t		vaddpd			 (%%r10),%%ymm10,%%ymm10	\n\t"\
		"vmovaps		%%ymm6,      (%%rax)			\n\t		vaddpd		0x020(%%r10),%%ymm11,%%ymm11	\n\t"\
		"vmovaps		%%ymm7, 0x020(%%rax)			\n\t		vsubpd		%%ymm12,%%ymm8,%%ymm8			\n\t"\
		"vsubpd		%%ymm5,%%ymm0,%%ymm0			\n\t		vsubpd		%%ymm13,%%ymm9,%%ymm9			\n\t"\
		"vsubpd		%%ymm4,%%ymm1,%%ymm1			\n\t		vmulpd		(%%rdi),%%ymm12,%%ymm12		\n\t"\
		"vmulpd		(%%rdi),%%ymm5,%%ymm5			\n\t		vmulpd		(%%rdi),%%ymm13,%%ymm13		\n\t"\
		"vmulpd		(%%rdi),%%ymm4,%%ymm4			\n\t		vmovaps		%%ymm8 ,     (%%r12)				\n\t"\
		"vmovaps		%%ymm0,      (%%rbx)			\n\t		vmovaps		%%ymm9 ,0x020(%%r12)				\n\t"\
		"vmovaps		%%ymm1, 0x020(%%rdx)			\n\t		vaddpd		%%ymm8,%%ymm12,%%ymm12		\n\t"\
		"vaddpd		%%ymm0,%%ymm5,%%ymm5			\n\t		vaddpd		%%ymm9,%%ymm13,%%ymm13		\n\t"\
		"vaddpd		%%ymm1,%%ymm4,%%ymm4			\n\t		vmovaps		%%ymm12,     (%%r10)				\n\t"\
		"vmovaps		%%ymm5,      (%%rdx)			\n\t		vmovaps		%%ymm13,0x020(%%r10)				\n\t"\
		"vmovaps		%%ymm4, 0x020(%%rbx)			\n\t		vsubpd		%%ymm15,%%ymm10,%%ymm10		\n\t"\
		"/*...Block 4: t06,t16,t26,t36	*/\n\t"\
		"addq		$0x080,%%rax						\n\t		vsubpd		%%ymm14,%%ymm11,%%ymm11		\n\t"\
		"addq		$0x080,%%rbx						\n\t		vaddpd		%%ymm15,%%ymm15,%%ymm15		\n\t"\
		"addq		$0x080,%%rcx						\n\t		vaddpd		%%ymm14,%%ymm14,%%ymm14		\n\t"\
		"addq		$0x080,%%rdx						\n\t		vmovaps		%%ymm10,     (%%r11)				\n\t"\
		"vmovaps			  (%%rcx),%%ymm4				\n\t		vmovaps		%%ymm11,0x020(%%r13)				\n\t"\
		"vmovaps			  (%%rdx),%%ymm6				\n\t		vaddpd		%%ymm10,%%ymm15,%%ymm15		\n\t"\
		"vmovaps		 0x020(%%rcx),%%ymm5				\n\t		vaddpd		%%ymm11,%%ymm14,%%ymm14		\n\t"\
		"vmovaps		 0x020(%%rdx),%%ymm7				\n\t		vmovaps		%%ymm15,     (%%r13)		\n\t"\
		"vmovaps			  (%%rcx),%%ymm0				\n\t		vmovaps		%%ymm14,0x020(%%r11)		\n\t"\
		"																/*...Block 8: t0E,t1E,t2E,t3E	*/		\n\t"\
		"vmovaps			  (%%rdx),%%ymm2				\n\t		addq		$0x080,%%r10			\n\t"\
		"vmovaps		 0x020(%%rcx),%%ymm1				\n\t		addq		$0x080,%%r11			\n\t"\
		"vmovaps		 0x020(%%rdx),%%ymm3				\n\t		addq		$0x080,%%r12			\n\t"\
		"vmulpd		 0x040(%%rsi),%%ymm4,%%ymm4			\n\t		addq		$0x080,%%r13			\n\t"\
		"vmulpd		 0x020(%%rsi),%%ymm6,%%ymm6			\n\t		vmovaps			  (%%r12),%%ymm12	\n\t"\
		"vmulpd		 0x060(%%rsi),%%ymm1,%%ymm1			\n\t		vmovaps			  (%%r13),%%ymm14	\n\t"\
		"vmulpd			  (%%rsi),%%ymm3,%%ymm3			\n\t		vmovaps		 0x020(%%r12),%%ymm13	\n\t"\
		"vmulpd		 0x040(%%rsi),%%ymm5,%%ymm5			\n\t		vmovaps		 0x020(%%r13),%%ymm15	\n\t"\
		"vmulpd		 0x020(%%rsi),%%ymm7,%%ymm7			\n\t		vmovaps			  (%%r12),%%ymm8		\n\t"\
		"vmulpd		 0x060(%%rsi),%%ymm0,%%ymm0			\n\t		vmovaps			  (%%r13),%%ymm10	\n\t"\
		"vmulpd			  (%%rsi),%%ymm2,%%ymm2			\n\t		vmovaps		 0x020(%%r12),%%ymm9		\n\t"\
		"vsubpd		%%ymm1,%%ymm4,%%ymm4			\n\t		vmovaps		 0x020(%%r13),%%ymm11	\n\t"\
		"vaddpd		%%ymm3,%%ymm6,%%ymm6			\n\t		vmulpd		0x020(%%rsi),%%ymm12,%%ymm12	\n\t"\
		"vaddpd		%%ymm0,%%ymm5,%%ymm5			\n\t		vmulpd		0x060(%%rsi),%%ymm14,%%ymm14	\n\t"\
		"vsubpd		%%ymm2,%%ymm7,%%ymm7			\n\t		vmulpd			 (%%rsi),%%ymm9,%%ymm9		\n\t"\
		"vsubpd		%%ymm6,%%ymm4,%%ymm4			\n\t		vmulpd		0x040(%%rsi),%%ymm11,%%ymm11	\n\t"\
		"vsubpd		%%ymm7,%%ymm5,%%ymm5			\n\t		vmulpd		0x020(%%rsi),%%ymm13,%%ymm13	\n\t"\
		"vmulpd		(%%rdi),%%ymm6,%%ymm6			\n\t		vmulpd		0x060(%%rsi),%%ymm15,%%ymm15	\n\t"\
		"vmulpd		(%%rdi),%%ymm7,%%ymm7			\n\t		vmulpd			 (%%rsi),%%ymm8,%%ymm8		\n\t"\
		"vaddpd		%%ymm4,%%ymm6,%%ymm6			\n\t		vmulpd		0x040(%%rsi),%%ymm10,%%ymm10	\n\t"\
		"vaddpd		%%ymm5,%%ymm7,%%ymm7			\n\t		vsubpd		%%ymm9,%%ymm12,%%ymm12		\n\t"\
		"vmovaps			  (%%rbx),%%ymm2				\n\t		vsubpd		%%ymm11,%%ymm14,%%ymm14		\n\t"\
		"vmovaps		 0x020(%%rbx),%%ymm0				\n\t		vaddpd		%%ymm8,%%ymm13,%%ymm13		\n\t"\
		"vmovaps			  (%%rbx),%%ymm1				\n\t		vaddpd		%%ymm10,%%ymm15,%%ymm15		\n\t"\
		"vmovaps		 0x020(%%rbx),%%ymm3				\n\t		vsubpd		%%ymm14,%%ymm12,%%ymm12		\n\t"\
		"vmulpd		-0x020(%%rsi),%%ymm2,%%ymm2			\n\t		vsubpd		%%ymm15,%%ymm13,%%ymm13		\n\t"\
		"vmulpd		-0x040(%%rsi),%%ymm0,%%ymm0			\n\t		vmulpd		(%%rdi),%%ymm14,%%ymm14		\n\t"\
		"vmulpd		-0x020(%%rsi),%%ymm3,%%ymm3			\n\t		vmulpd		(%%rdi),%%ymm15,%%ymm15		\n\t"\
		"vmulpd		-0x040(%%rsi),%%ymm1,%%ymm1			\n\t		vaddpd		%%ymm12,%%ymm14,%%ymm14		\n\t"\
		"vsubpd		%%ymm0,%%ymm2,%%ymm2			\n\t		vaddpd		%%ymm13,%%ymm15,%%ymm15		\n\t"\
		"vaddpd		%%ymm1,%%ymm3,%%ymm3			\n\t		vmovaps			 (%%r11),%%ymm10				\n\t"\
		"vmovaps			  (%%rax),%%ymm0				\n\t		vmovaps		0x020(%%r11),%%ymm8					\n\t"\
		"vmovaps		 0x020(%%rax),%%ymm1				\n\t		vmovaps			 (%%r11),%%ymm9					\n\t"\
		"vsubpd		%%ymm2,%%ymm0,%%ymm0			\n\t		vmovaps		0x020(%%r11),%%ymm11				\n\t"\
		"vsubpd		%%ymm3,%%ymm1,%%ymm1			\n\t		vmulpd		-0x40(%%rsi),%%ymm10,%%ymm10	\n\t"\
		"vaddpd			  (%%rax),%%ymm2,%%ymm2			\n\t		vmulpd		-0x20(%%rsi),%%ymm8,%%ymm8		\n\t"\
		"vaddpd		 0x020(%%rax),%%ymm3,%%ymm3			\n\t		vmulpd		-0x40(%%rsi),%%ymm11,%%ymm11	\n\t"\
		"vsubpd		%%ymm4,%%ymm2,%%ymm2			\n\t		vmulpd		-0x20(%%rsi),%%ymm9,%%ymm9		\n\t"\
		"vsubpd		%%ymm5,%%ymm3,%%ymm3			\n\t		vaddpd		%%ymm8,%%ymm10,%%ymm10		\n\t"\
		"vmulpd		(%%rdi),%%ymm4,%%ymm4			\n\t		vsubpd		%%ymm9,%%ymm11,%%ymm11		\n\t"\
		"vmulpd		(%%rdi),%%ymm5,%%ymm5			\n\t		vmovaps			 (%%r10),%%ymm8					\n\t"\
		"vmovaps		%%ymm2,      (%%rcx)			\n\t		vmovaps		0x020(%%r10),%%ymm9					\n\t"\
		"vmovaps		%%ymm3, 0x020(%%rcx)			\n\t		vsubpd		%%ymm10,%%ymm8,%%ymm8			\n\t"\
		"vaddpd		%%ymm2,%%ymm4,%%ymm4			\n\t		vsubpd		%%ymm11,%%ymm9,%%ymm9			\n\t"\
		"vaddpd		%%ymm3,%%ymm5,%%ymm5			\n\t		vaddpd			 (%%r10),%%ymm10,%%ymm10	\n\t"\
		"vmovaps		%%ymm4,      (%%rax)			\n\t		vaddpd		0x020(%%r10),%%ymm11,%%ymm11	\n\t"\
		"vmovaps		%%ymm5, 0x020(%%rax)			\n\t		vsubpd		%%ymm12,%%ymm8,%%ymm8			\n\t"\
		"vsubpd		%%ymm7,%%ymm0,%%ymm0			\n\t		vsubpd		%%ymm13,%%ymm9,%%ymm9			\n\t"\
		"vsubpd		%%ymm6,%%ymm1,%%ymm1			\n\t		vmulpd		(%%rdi),%%ymm12,%%ymm12		\n\t"\
		"vmulpd		(%%rdi),%%ymm7,%%ymm7			\n\t		vmulpd		(%%rdi),%%ymm13,%%ymm13		\n\t"\
		"vmulpd		(%%rdi),%%ymm6,%%ymm6			\n\t		vmovaps		%%ymm8 ,     (%%r12)				\n\t"\
		"vmovaps		%%ymm0,      (%%rbx)			\n\t		vmovaps		%%ymm9 ,0x020(%%r12)				\n\t"\
		"vmovaps		%%ymm1, 0x020(%%rdx)			\n\t		vaddpd		%%ymm8,%%ymm12,%%ymm12		\n\t"\
		"vaddpd		%%ymm0,%%ymm7,%%ymm7			\n\t		vaddpd		%%ymm9,%%ymm13,%%ymm13		\n\t"\
		"vaddpd		%%ymm1,%%ymm6,%%ymm6			\n\t		vmovaps		%%ymm12,     (%%r10)				\n\t"\
		"vmovaps		%%ymm7,      (%%rdx)			\n\t		vmovaps		%%ymm13,0x020(%%r10)				\n\t"\
		"vmovaps		%%ymm6, 0x020(%%rbx)			\n\t		vsubpd		%%ymm15,%%ymm10,%%ymm10		\n\t"\
		"																vsubpd		%%ymm14,%%ymm11,%%ymm11		\n\t"\
		"																vmulpd		(%%rdi),%%ymm15,%%ymm15		\n\t"\
		"																vmulpd		(%%rdi),%%ymm14,%%ymm14		\n\t"\
		"																vmovaps		%%ymm10,     (%%r11)				\n\t"\
		"																vmovaps		%%ymm11,0x020(%%r13)				\n\t"\
		"																vaddpd		%%ymm10,%%ymm15,%%ymm15		\n\t"\
		"																vaddpd		%%ymm11,%%ymm14,%%ymm14		\n\t"\
		"																vmovaps		%%ymm15,     (%%r13)		\n\t"\
		"																vmovaps		%%ymm14,0x020(%%r11)		\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__add1] "m" (Xadd1)\
		 ,[__add2] "m" (Xadd2)\
		 ,[__add3] "m" (Xadd3)\
		 ,[__r00] "m" (Xr00)\
		 ,[__r10] "m" (Xr10)\
		 ,[__r20] "m" (Xr20)\
		 ,[__r30] "m" (Xr30)\
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__cc0] "m" (Xcc0)\
		 ,[__c00] "m" (Xc00)\
		 ,[__c01] "m" (Xc01)\
		 ,[__c02] "m" (Xc02)\
		 ,[__c03] "m" (Xc03)\
		 ,[__c05] "m" (Xc05)\
		 ,[__c07] "m" (Xc07)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r10","r11","r12","r13","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"		/* Clobbered registers */\
	);\
	}

	#define SSE2_RADIX32_WRAPPER_DIT(Xadd0,Xadd1,Xadd2,Xadd3,Xisrt2,Xr00,Xr08,Xr10,Xr20,Xr28,Xr30,Xc01,Xc02,Xc04,Xc06,Xc08,Xc0A,Xc0C,Xc0E,Xc10,Xc12,Xc14,Xc16,Xc18,Xc1A,Xc1C,Xc1E)\
	{\
	__asm__ volatile (\
	/************************************************************************/\
	/*...And do an inverse DIT radix-32 pass on the squared-data blocks:	*/\
	/************************************************************************/\
		"/*...Block 1: */								\n\t"\
		"movq		%[__isrt2],%%rsi				\n\t"\
		"movq		%[__r00],%%rax				\n\t"\
		"movq		%%rax,%%rbx					\n\t"\
		"movq		%%rax,%%rcx					\n\t"\
		"movq		%%rax,%%rdx					\n\t"\
		"addq		$0x400,%%rbx					\n\t"\
		"addq		$0x200,%%rcx					\n\t"\
		"addq		$0x600,%%rdx					\n\t"\
		"/*****	SSE2_RADIX4_DIT_IN_PLACE()	*****/		\n\t		/*...Block 2 has tmp-addresses offset +0x80 w.r.to Block 1:	*/\n\t"\
		"vmovaps		  (%%rax),%%ymm0				\n\t		vmovaps		 0x080(%%rax),%%ymm8 			\n\t"\
		"vmovaps	 0x020(%%rax),%%ymm1				\n\t		vmovaps		 0x0a0(%%rax),%%ymm9 			\n\t"\
		"vmovaps		  (%%rax),%%ymm2				\n\t		vmovaps		 0x080(%%rax),%%ymm10			\n\t"\
		"vmovaps	 0x020(%%rax),%%ymm3				\n\t		vmovaps		 0x0a0(%%rax),%%ymm11			\n\t"\
		"vaddpd			  (%%rbx),%%ymm0,%%ymm0		\n\t		vaddpd		 0x080(%%rbx),%%ymm8 ,%%ymm8 	\n\t"\
		"vaddpd		 0x020(%%rbx),%%ymm1,%%ymm1		\n\t		vaddpd		 0x0a0(%%rbx),%%ymm9 ,%%ymm9 	\n\t"\
		"vsubpd			  (%%rbx),%%ymm2,%%ymm2		\n\t		vsubpd		 0x080(%%rbx),%%ymm10,%%ymm10	\n\t"\
		"vsubpd		 0x020(%%rbx),%%ymm3,%%ymm3		\n\t		vsubpd		 0x0a0(%%rbx),%%ymm11,%%ymm11	\n\t"\
		"vmovaps		  (%%rcx),%%ymm4				\n\t		vmovaps		 0x080(%%rcx),%%ymm12			\n\t"\
		"vmovaps	 0x020(%%rcx),%%ymm5				\n\t		vmovaps		 0x0a0(%%rcx),%%ymm13			\n\t"\
		"vmovaps		  (%%rcx),%%ymm6				\n\t		vmovaps		 0x080(%%rcx),%%ymm14			\n\t"\
		"vmovaps	 0x020(%%rcx),%%ymm7				\n\t		vmovaps		 0x0a0(%%rcx),%%ymm15			\n\t"\
		"vaddpd			  (%%rdx),%%ymm4,%%ymm4		\n\t		vaddpd		 0x080(%%rdx),%%ymm12,%%ymm12	\n\t"\
		"vaddpd		 0x020(%%rdx),%%ymm5,%%ymm5		\n\t		vaddpd		 0x0a0(%%rdx),%%ymm13,%%ymm13	\n\t"\
		"vsubpd			  (%%rdx),%%ymm6,%%ymm6		\n\t		vsubpd		 0x080(%%rdx),%%ymm14,%%ymm14	\n\t"\
		"vsubpd		 0x020(%%rdx),%%ymm7,%%ymm7		\n\t		vsubpd		 0x0a0(%%rdx),%%ymm15,%%ymm15	\n\t"\
		"vsubpd			%%ymm4,%%ymm0,%%ymm0		\n\t		vsubpd			%%ymm12,%%ymm8 ,%%ymm8 	\n\t"\
		"vsubpd			%%ymm5,%%ymm1,%%ymm1		\n\t		vsubpd			%%ymm13,%%ymm9 ,%%ymm9 	\n\t"\
		"vmovaps		%%ymm0,      (%%rbx)		\n\t		vmovaps			%%ymm8 , 0x080(%%rbx)			\n\t"\
		"vmovaps		%%ymm1, 0x020(%%rbx)		\n\t		vmovaps			%%ymm9 , 0x0a0(%%rbx)			\n\t"\
		"vaddpd			%%ymm4,%%ymm4,%%ymm4		\n\t		vaddpd			%%ymm12,%%ymm12,%%ymm12	\n\t"\
		"vaddpd			%%ymm5,%%ymm5,%%ymm5		\n\t		vaddpd			%%ymm13,%%ymm13,%%ymm13	\n\t"\
		"vaddpd			%%ymm0,%%ymm4,%%ymm4		\n\t		vaddpd			%%ymm8 ,%%ymm12,%%ymm12	\n\t"\
		"vaddpd			%%ymm1,%%ymm5,%%ymm5		\n\t		vaddpd			%%ymm9 ,%%ymm13,%%ymm13	\n\t"\
		"vmovaps		%%ymm4,      (%%rax)		\n\t		vmovaps			%%ymm12, 0x080(%%rax)			\n\t"\
		"vmovaps		%%ymm5, 0x020(%%rax)		\n\t		vmovaps			%%ymm13, 0x0a0(%%rax)			\n\t"\
		"vsubpd			%%ymm7,%%ymm2,%%ymm2		\n\t		vsubpd			%%ymm15,%%ymm10,%%ymm10	\n\t"\
		"vsubpd			%%ymm6,%%ymm3,%%ymm3		\n\t		vsubpd			%%ymm14,%%ymm11,%%ymm11	\n\t"\
		"vmovaps		%%ymm2,      (%%rdx)		\n\t		vmovaps			%%ymm10, 0x080(%%rdx)			\n\t"\
		"vmovaps		%%ymm3, 0x020(%%rcx)		\n\t		vmovaps			%%ymm11, 0x0a0(%%rcx)			\n\t"\
		"vaddpd			%%ymm7,%%ymm7,%%ymm7		\n\t		vaddpd			%%ymm15,%%ymm15,%%ymm15	\n\t"\
		"vaddpd			%%ymm6,%%ymm6,%%ymm6		\n\t		vaddpd			%%ymm14,%%ymm14,%%ymm14	\n\t"\
		"vaddpd			%%ymm2,%%ymm7,%%ymm7		\n\t		vaddpd			%%ymm10,%%ymm15,%%ymm15	\n\t"\
		"vaddpd			%%ymm3,%%ymm6,%%ymm6		\n\t		vaddpd			%%ymm11,%%ymm14,%%ymm14	\n\t"\
		"vmovaps		%%ymm7,      (%%rcx)		\n\t		vmovaps			%%ymm15, 0x080(%%rcx)			\n\t"\
		"vmovaps		%%ymm6, 0x020(%%rdx)		\n\t		vmovaps			%%ymm14, 0x0a0(%%rdx)			\n\t"\
		"addq		$0x100,%%rax					\n\t"\
		"addq		$0x100,%%rbx					\n\t"\
		"addq		$0x100,%%rcx					\n\t"\
		"addq		$0x100,%%rdx					\n\t"\
		"vmovaps		  (%%rax),%%ymm0				\n\t		vmovaps		 0x080(%%rax),%%ymm8 			\n\t"\
		"vmovaps	 0x020(%%rax),%%ymm1				\n\t		vmovaps		 0x0a0(%%rax),%%ymm9 			\n\t"\
		"vmovaps		  (%%rax),%%ymm2				\n\t		vmovaps		 0x080(%%rax),%%ymm10			\n\t"\
		"vmovaps	 0x020(%%rax),%%ymm3				\n\t		vmovaps		 0x0a0(%%rax),%%ymm11			\n\t"\
		"vaddpd			  (%%rbx),%%ymm0,%%ymm0		\n\t		vaddpd		 0x080(%%rbx),%%ymm8 ,%%ymm8 	\n\t"\
		"vaddpd		 0x020(%%rbx),%%ymm1,%%ymm1		\n\t		vaddpd		 0x0a0(%%rbx),%%ymm9 ,%%ymm9 	\n\t"\
		"vsubpd			  (%%rbx),%%ymm2,%%ymm2		\n\t		vsubpd		 0x080(%%rbx),%%ymm10,%%ymm10	\n\t"\
		"vsubpd		 0x020(%%rbx),%%ymm3,%%ymm3		\n\t		vsubpd		 0x0a0(%%rbx),%%ymm11,%%ymm11	\n\t"\
		"vmovaps		  (%%rcx),%%ymm4				\n\t		vmovaps		 0x080(%%rcx),%%ymm12			\n\t"\
		"vmovaps	 0x020(%%rcx),%%ymm5				\n\t		vmovaps		 0x0a0(%%rcx),%%ymm13			\n\t"\
		"vmovaps		  (%%rcx),%%ymm6				\n\t		vmovaps		 0x080(%%rcx),%%ymm14			\n\t"\
		"vmovaps	 0x020(%%rcx),%%ymm7				\n\t		vmovaps		 0x0a0(%%rcx),%%ymm15			\n\t"\
		"vaddpd			  (%%rdx),%%ymm4,%%ymm4		\n\t		vaddpd		 0x080(%%rdx),%%ymm12,%%ymm12	\n\t"\
		"vaddpd		 0x020(%%rdx),%%ymm5,%%ymm5		\n\t		vaddpd		 0x0a0(%%rdx),%%ymm13,%%ymm13	\n\t"\
		"vsubpd			  (%%rdx),%%ymm6,%%ymm6		\n\t		vsubpd		 0x080(%%rdx),%%ymm14,%%ymm14	\n\t"\
		"vsubpd		 0x020(%%rdx),%%ymm7,%%ymm7		\n\t		vsubpd		 0x0a0(%%rdx),%%ymm15,%%ymm15	\n\t"\
		"vsubpd			%%ymm4,%%ymm0,%%ymm0		\n\t		vsubpd			%%ymm12,%%ymm8 ,%%ymm8 	\n\t"\
		"vsubpd			%%ymm5,%%ymm1,%%ymm1		\n\t		vsubpd			%%ymm13,%%ymm9 ,%%ymm9 	\n\t"\
		"vmovaps		%%ymm0,      (%%rbx)		\n\t		vmovaps			%%ymm8 , 0x080(%%rbx)			\n\t"\
		"vmovaps		%%ymm1, 0x020(%%rbx)		\n\t		vmovaps			%%ymm9 , 0x0a0(%%rbx)			\n\t"\
		"vaddpd			%%ymm4,%%ymm4,%%ymm4		\n\t		vaddpd			%%ymm12,%%ymm12,%%ymm12	\n\t"\
		"vaddpd			%%ymm5,%%ymm5,%%ymm5		\n\t		vaddpd			%%ymm13,%%ymm13,%%ymm13	\n\t"\
		"vaddpd			%%ymm0,%%ymm4,%%ymm4		\n\t		vaddpd			%%ymm8 ,%%ymm12,%%ymm12	\n\t"\
		"vaddpd			%%ymm1,%%ymm5,%%ymm5		\n\t		vaddpd			%%ymm9 ,%%ymm13,%%ymm13	\n\t"\
		"vmovaps		%%ymm4,      (%%rax)		\n\t		vmovaps			%%ymm12, 0x080(%%rax)			\n\t"\
		"vmovaps		%%ymm5, 0x020(%%rax)		\n\t		vmovaps			%%ymm13, 0x0a0(%%rax)			\n\t"\
		"vsubpd			%%ymm7,%%ymm2,%%ymm2		\n\t		vsubpd			%%ymm15,%%ymm10,%%ymm10	\n\t"\
		"vsubpd			%%ymm6,%%ymm3,%%ymm3		\n\t		vsubpd			%%ymm14,%%ymm11,%%ymm11	\n\t"\
		"vaddpd			%%ymm7,%%ymm7,%%ymm7		\n\t		vaddpd			%%ymm15,%%ymm15,%%ymm15	\n\t"\
		"vaddpd			%%ymm6,%%ymm6,%%ymm6		\n\t		vaddpd			%%ymm14,%%ymm14,%%ymm14	\n\t"\
		"vaddpd			%%ymm2,%%ymm7,%%ymm7		\n\t		vaddpd			%%ymm10,%%ymm15,%%ymm15	\n\t"\
		"vaddpd			%%ymm3,%%ymm6,%%ymm6		\n\t		vaddpd			%%ymm11,%%ymm14,%%ymm14	\n\t"\
		"vmovaps		%%ymm3,%%ymm0				\n\t		vmovaps			%%ymm11,%%ymm8 				\n\t"\
		"vmovaps		%%ymm6,%%ymm1				\n\t		vmovaps			%%ymm14,%%ymm9 				\n\t"\
		"vsubpd			%%ymm7,%%ymm3,%%ymm3		\n\t		vsubpd			%%ymm15,%%ymm11,%%ymm11	\n\t"\
		"vsubpd			%%ymm2,%%ymm6,%%ymm6		\n\t		vsubpd			%%ymm10,%%ymm14,%%ymm14	\n\t"\
		"vaddpd			%%ymm7,%%ymm0,%%ymm0		\n\t		vaddpd			%%ymm15,%%ymm8 ,%%ymm8 	\n\t"\
		"vaddpd			%%ymm2,%%ymm1,%%ymm1		\n\t		vaddpd			%%ymm10,%%ymm9 ,%%ymm9 	\n\t"\
		"vmulpd			  (%%rsi),%%ymm3,%%ymm3		\n\t		vmulpd			  (%%rsi),%%ymm11,%%ymm11	\n\t"\
		"vmulpd			  (%%rsi),%%ymm6,%%ymm6		\n\t		vmulpd			  (%%rsi),%%ymm14,%%ymm14	\n\t"\
		"vmulpd			  (%%rsi),%%ymm0,%%ymm0		\n\t		vmulpd			  (%%rsi),%%ymm8 ,%%ymm8 	\n\t"\
		"vmulpd			  (%%rsi),%%ymm1,%%ymm1		\n\t		vmulpd			  (%%rsi),%%ymm9 ,%%ymm9 	\n\t"\
		"vmovaps		%%ymm3, 0x020(%%rcx)		\n\t		vmovaps			%%ymm11, 0x0a0(%%rcx)			\n\t"\
		"vmovaps		%%ymm6, 0x020(%%rdx)		\n\t		vmovaps			%%ymm14, 0x0a0(%%rdx)			\n\t"\
		"vmovaps		%%ymm0,      (%%rcx)		\n\t		vmovaps			%%ymm8 , 0x080(%%rcx)			\n\t"\
		"vmovaps		%%ymm1,      (%%rdx)		\n\t		vmovaps			%%ymm9 , 0x080(%%rdx)			\n\t"\
		"/*****	SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS *****/\n\t		/*****	SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS *****/\n\t"\
		"/*        (r00,r10,r20,r30,r08,r18,r28,r38): */\n\t		/*        (r04,r14,r24,r34,r0C,r1C,r2C,r3C):  */\n\t"\
		"vmovaps		-0x100(%%rax),%%ymm0			\n\t		vmovaps		-0x080(%%rax),%%ymm8 			\n\t"\
		"vmovaps		-0x100(%%rbx),%%ymm4			\n\t		vmovaps		-0x080(%%rbx),%%ymm12			\n\t"\
		"vmovaps		-0x0e0(%%rax),%%ymm1			\n\t		vmovaps		-0x060(%%rax),%%ymm9 			\n\t"\
		"vmovaps		-0x0e0(%%rbx),%%ymm5			\n\t		vmovaps		-0x060(%%rbx),%%ymm13			\n\t"\
		"vmovaps			  (%%rax),%%ymm2			\n\t		vmovaps		 0x080(%%rax),%%ymm10			\n\t"\
		"vmovaps		 0x020(%%rbx),%%ymm7			\n\t		vmovaps		 0x0a0(%%rbx),%%ymm15			\n\t"\
		"vmovaps		 0x020(%%rax),%%ymm3			\n\t		vmovaps		 0x0a0(%%rax),%%ymm11			\n\t"\
		"vmovaps			  (%%rbx),%%ymm6			\n\t		vmovaps		 0x080(%%rbx),%%ymm14			\n\t"\
		"vsubpd			%%ymm2,%%ymm0,%%ymm0		\n\t		vsubpd		%%ymm10,%%ymm8 ,%%ymm8 		\n\t"\
		"vsubpd			%%ymm7,%%ymm4,%%ymm4		\n\t		vsubpd		%%ymm15,%%ymm12,%%ymm12		\n\t"\
		"vsubpd			%%ymm3,%%ymm1,%%ymm1		\n\t		vsubpd		%%ymm11,%%ymm9 ,%%ymm9 		\n\t"\
		"vsubpd			%%ymm6,%%ymm5,%%ymm5		\n\t		vsubpd		%%ymm14,%%ymm13,%%ymm13		\n\t"\
		"vaddpd			%%ymm2,%%ymm2,%%ymm2		\n\t		vaddpd		%%ymm10,%%ymm10,%%ymm10		\n\t"\
		"vaddpd			%%ymm7,%%ymm7,%%ymm7		\n\t		vaddpd		%%ymm15,%%ymm15,%%ymm15		\n\t"\
		"vaddpd			%%ymm3,%%ymm3,%%ymm3		\n\t		vaddpd		%%ymm11,%%ymm11,%%ymm11		\n\t"\
		"vaddpd			%%ymm6,%%ymm6,%%ymm6		\n\t		vaddpd		%%ymm14,%%ymm14,%%ymm14		\n\t"\
		"vaddpd			%%ymm0,%%ymm2,%%ymm2		\n\t		vaddpd		%%ymm8 ,%%ymm10,%%ymm10		\n\t"\
		"vaddpd			%%ymm4,%%ymm7,%%ymm7		\n\t		vaddpd		%%ymm12,%%ymm15,%%ymm15		\n\t"\
		"vaddpd			%%ymm1,%%ymm3,%%ymm3		\n\t		vaddpd		%%ymm9 ,%%ymm11,%%ymm11		\n\t"\
		"vaddpd			%%ymm5,%%ymm6,%%ymm6		\n\t		vaddpd		%%ymm13,%%ymm14,%%ymm14		\n\t"\
		"vmovaps		%%ymm0,      (%%rax)		\n\t		vmovaps		%%ymm8 , 0x080(%%rax)			\n\t"\
		"vmovaps		%%ymm4,      (%%rbx)		\n\t		vmovaps		%%ymm12, 0x080(%%rbx)			\n\t"\
		"vmovaps		%%ymm1, 0x020(%%rax)		\n\t		vmovaps		%%ymm9 , 0x0a0(%%rax)			\n\t"\
		"vmovaps		%%ymm5,-0x0e0(%%rbx)		\n\t		vmovaps		%%ymm13,-0x060(%%rbx)			\n\t"\
		"vmovaps		%%ymm2,-0x100(%%rax)		\n\t		vmovaps		%%ymm10,-0x080(%%rax)			\n\t"\
		"vmovaps		%%ymm7,-0x100(%%rbx)		\n\t		vmovaps		%%ymm15,-0x080(%%rbx)			\n\t"\
		"vmovaps		%%ymm3,-0x0e0(%%rax)		\n\t		vmovaps		%%ymm11,-0x060(%%rax)			\n\t"\
		"vmovaps		%%ymm6, 0x020(%%rbx)		\n\t		vmovaps		%%ymm14, 0x0a0(%%rbx)			\n\t"\
		"vmovaps		-0x100(%%rcx),%%ymm0			\n\t		vmovaps		-0x080(%%rcx),%%ymm8 			\n\t"\
		"vmovaps		-0x100(%%rdx),%%ymm4			\n\t		vmovaps		-0x080(%%rdx),%%ymm12			\n\t"\
		"vmovaps		-0x0e0(%%rcx),%%ymm1			\n\t		vmovaps		-0x060(%%rcx),%%ymm9 			\n\t"\
		"vmovaps		-0x0e0(%%rdx),%%ymm5			\n\t		vmovaps		-0x060(%%rdx),%%ymm13			\n\t"\
		"vmovaps			  (%%rcx),%%ymm2			\n\t		vmovaps		 0x080(%%rcx),%%ymm10			\n\t"\
		"vmovaps		 0x020(%%rdx),%%ymm7			\n\t		vmovaps		 0x0a0(%%rdx),%%ymm15			\n\t"\
		"vmovaps		 0x020(%%rcx),%%ymm3			\n\t		vmovaps		 0x0a0(%%rcx),%%ymm11			\n\t"\
		"vmovaps			  (%%rdx),%%ymm6			\n\t		vmovaps		 0x080(%%rdx),%%ymm14			\n\t"\
		"vsubpd			%%ymm2,%%ymm0,%%ymm0		\n\t		vsubpd		%%ymm10,%%ymm8 ,%%ymm8 		\n\t"\
		"vsubpd			%%ymm7,%%ymm4,%%ymm4		\n\t		vsubpd		%%ymm15,%%ymm12,%%ymm12		\n\t"\
		"vsubpd			%%ymm3,%%ymm1,%%ymm1		\n\t		vsubpd		%%ymm11,%%ymm9 ,%%ymm9 		\n\t"\
		"vsubpd			%%ymm6,%%ymm5,%%ymm5		\n\t		vsubpd		%%ymm14,%%ymm13,%%ymm13		\n\t"\
		"vaddpd			%%ymm2,%%ymm2,%%ymm2		\n\t		vaddpd		%%ymm10,%%ymm10,%%ymm10		\n\t"\
		"vaddpd			%%ymm7,%%ymm7,%%ymm7		\n\t		vaddpd		%%ymm15,%%ymm15,%%ymm15		\n\t"\
		"vaddpd			%%ymm3,%%ymm3,%%ymm3		\n\t		vaddpd		%%ymm11,%%ymm11,%%ymm11		\n\t"\
		"vaddpd			%%ymm6,%%ymm6,%%ymm6		\n\t		vaddpd		%%ymm14,%%ymm14,%%ymm14		\n\t"\
		"vaddpd			%%ymm0,%%ymm2,%%ymm2		\n\t		vaddpd		%%ymm8 ,%%ymm10,%%ymm10		\n\t"\
		"vaddpd			%%ymm4,%%ymm7,%%ymm7		\n\t		vaddpd		%%ymm12,%%ymm15,%%ymm15		\n\t"\
		"vaddpd			%%ymm1,%%ymm3,%%ymm3		\n\t		vaddpd		%%ymm9 ,%%ymm11,%%ymm11		\n\t"\
		"vaddpd			%%ymm5,%%ymm6,%%ymm6		\n\t		vaddpd		%%ymm13,%%ymm14,%%ymm14		\n\t"\
		"vmovaps		%%ymm0,      (%%rcx)		\n\t		vmovaps		%%ymm8 , 0x080(%%rcx)			\n\t"\
		"vmovaps		%%ymm4,      (%%rdx)		\n\t		vmovaps		%%ymm12, 0x080(%%rdx)			\n\t"\
		"vmovaps		%%ymm1, 0x020(%%rcx)		\n\t		vmovaps		%%ymm9 , 0x0a0(%%rcx)			\n\t"\
		"vmovaps		%%ymm5,-0x0e0(%%rdx)		\n\t		vmovaps		%%ymm13,-0x060(%%rdx)			\n\t"\
		"vmovaps		%%ymm2,-0x100(%%rcx)		\n\t		vmovaps		%%ymm10,-0x080(%%rcx)			\n\t"\
		"vmovaps		%%ymm7,-0x100(%%rdx)		\n\t		vmovaps		%%ymm15,-0x080(%%rdx)			\n\t"\
		"vmovaps		%%ymm3,-0x0e0(%%rcx)		\n\t		vmovaps		%%ymm11,-0x060(%%rcx)			\n\t"\
		"vmovaps		%%ymm6, 0x020(%%rdx)		\n\t		vmovaps		%%ymm14, 0x0a0(%%rdx)			\n\t"\
		"/*...Blocks 3,4 have tmp-addresses offset +0x40 w.r.to Blocks 1,2, respectively (thus +0x100-0x0c0 = +0x040: */\n\t"\
		"subq		$0xc0,%%rax					\n\t"\
		"subq		$0xc0,%%rbx					\n\t"\
		"subq		$0xc0,%%rcx					\n\t"\
		"subq		$0xc0,%%rdx					\n\t"\
		"vmovaps		  (%%rax),%%ymm0				\n\t		vmovaps		 0x080(%%rax),%%ymm8 			\n\t"\
		"vmovaps	 0x020(%%rax),%%ymm1				\n\t		vmovaps		 0x0a0(%%rax),%%ymm9 			\n\t"\
		"vmovaps		  (%%rax),%%ymm2				\n\t		vmovaps		 0x080(%%rax),%%ymm10			\n\t"\
		"vmovaps	 0x020(%%rax),%%ymm3				\n\t		vmovaps		 0x0a0(%%rax),%%ymm11			\n\t"\
		"vaddpd			  (%%rbx),%%ymm0,%%ymm0		\n\t		vaddpd		 0x080(%%rbx),%%ymm8 ,%%ymm8 	\n\t"\
		"vaddpd		 0x020(%%rbx),%%ymm1,%%ymm1		\n\t		vaddpd		 0x0a0(%%rbx),%%ymm9 ,%%ymm9 	\n\t"\
		"vsubpd			  (%%rbx),%%ymm2,%%ymm2		\n\t		vsubpd		 0x080(%%rbx),%%ymm10,%%ymm10	\n\t"\
		"vsubpd		 0x020(%%rbx),%%ymm3,%%ymm3		\n\t		vsubpd		 0x0a0(%%rbx),%%ymm11,%%ymm11	\n\t"\
		"vmovaps		  (%%rcx),%%ymm4				\n\t		vmovaps		 0x080(%%rcx),%%ymm12			\n\t"\
		"vmovaps	 0x020(%%rcx),%%ymm5				\n\t		vmovaps		 0x0a0(%%rcx),%%ymm13			\n\t"\
		"vmovaps		  (%%rcx),%%ymm6				\n\t		vmovaps		 0x080(%%rcx),%%ymm14			\n\t"\
		"vmovaps	 0x020(%%rcx),%%ymm7				\n\t		vmovaps		 0x0a0(%%rcx),%%ymm15			\n\t"\
		"vaddpd			  (%%rdx),%%ymm4,%%ymm4		\n\t		vaddpd		 0x080(%%rdx),%%ymm12,%%ymm12	\n\t"\
		"vaddpd		 0x020(%%rdx),%%ymm5,%%ymm5		\n\t		vaddpd		 0x0a0(%%rdx),%%ymm13,%%ymm13	\n\t"\
		"vsubpd			  (%%rdx),%%ymm6,%%ymm6		\n\t		vsubpd		 0x080(%%rdx),%%ymm14,%%ymm14	\n\t"\
		"vsubpd		 0x020(%%rdx),%%ymm7,%%ymm7		\n\t		vsubpd		 0x0a0(%%rdx),%%ymm15,%%ymm15	\n\t"\
		"vsubpd			%%ymm4,%%ymm0,%%ymm0		\n\t		vsubpd			%%ymm12,%%ymm8 ,%%ymm8 	\n\t"\
		"vsubpd			%%ymm5,%%ymm1,%%ymm1		\n\t		vsubpd			%%ymm13,%%ymm9 ,%%ymm9 	\n\t"\
		"vmovaps		%%ymm0,      (%%rbx)		\n\t		vmovaps			%%ymm8 , 0x080(%%rbx)			\n\t"\
		"vmovaps		%%ymm1, 0x020(%%rbx)		\n\t		vmovaps			%%ymm9 , 0x0a0(%%rbx)			\n\t"\
		"vaddpd			%%ymm4,%%ymm4,%%ymm4		\n\t		vaddpd			%%ymm12,%%ymm12,%%ymm12	\n\t"\
		"vaddpd			%%ymm5,%%ymm5,%%ymm5		\n\t		vaddpd			%%ymm13,%%ymm13,%%ymm13	\n\t"\
		"vaddpd			%%ymm0,%%ymm4,%%ymm4		\n\t		vaddpd			%%ymm8 ,%%ymm12,%%ymm12	\n\t"\
		"vaddpd			%%ymm1,%%ymm5,%%ymm5		\n\t		vaddpd			%%ymm9 ,%%ymm13,%%ymm13	\n\t"\
		"vmovaps		%%ymm4,      (%%rax)		\n\t		vmovaps			%%ymm12, 0x080(%%rax)			\n\t"\
		"vmovaps		%%ymm5, 0x020(%%rax)		\n\t		vmovaps			%%ymm13, 0x0a0(%%rax)			\n\t"\
		"vsubpd			%%ymm7,%%ymm2,%%ymm2		\n\t		vsubpd			%%ymm15,%%ymm10,%%ymm10	\n\t"\
		"vsubpd			%%ymm6,%%ymm3,%%ymm3		\n\t		vsubpd			%%ymm14,%%ymm11,%%ymm11	\n\t"\
		"vmovaps		%%ymm2,      (%%rdx)		\n\t		vmovaps			%%ymm10, 0x080(%%rdx)			\n\t"\
		"vmovaps		%%ymm3, 0x020(%%rcx)		\n\t		vmovaps			%%ymm11, 0x0a0(%%rcx)			\n\t"\
		"vaddpd			%%ymm7,%%ymm7,%%ymm7		\n\t		vaddpd			%%ymm15,%%ymm15,%%ymm15	\n\t"\
		"vaddpd			%%ymm6,%%ymm6,%%ymm6		\n\t		vaddpd			%%ymm14,%%ymm14,%%ymm14	\n\t"\
		"vaddpd			%%ymm2,%%ymm7,%%ymm7		\n\t		vaddpd			%%ymm10,%%ymm15,%%ymm15	\n\t"\
		"vaddpd			%%ymm3,%%ymm6,%%ymm6		\n\t		vaddpd			%%ymm11,%%ymm14,%%ymm14	\n\t"\
		"vmovaps		%%ymm7,      (%%rcx)		\n\t		vmovaps			%%ymm15, 0x080(%%rcx)			\n\t"\
		"vmovaps		%%ymm6, 0x020(%%rdx)		\n\t		vmovaps			%%ymm14, 0x0a0(%%rdx)			\n\t"\
		"addq		$0x100,%%rax					\n\t"\
		"addq		$0x100,%%rbx					\n\t"\
		"addq		$0x100,%%rcx					\n\t"\
		"addq		$0x100,%%rdx					\n\t"\
		"vmovaps		  (%%rax),%%ymm0				\n\t		vmovaps		 0x080(%%rax),%%ymm8 			\n\t"\
		"vmovaps	 0x020(%%rax),%%ymm1				\n\t		vmovaps		 0x0a0(%%rax),%%ymm9 			\n\t"\
		"vmovaps		  (%%rax),%%ymm2				\n\t		vmovaps		 0x080(%%rax),%%ymm10			\n\t"\
		"vmovaps	 0x020(%%rax),%%ymm3				\n\t		vmovaps		 0x0a0(%%rax),%%ymm11			\n\t"\
		"vaddpd			  (%%rbx),%%ymm0,%%ymm0		\n\t		vaddpd		 0x080(%%rbx),%%ymm8 ,%%ymm8 	\n\t"\
		"vaddpd		 0x020(%%rbx),%%ymm1,%%ymm1		\n\t		vaddpd		 0x0a0(%%rbx),%%ymm9 ,%%ymm9 	\n\t"\
		"vsubpd			  (%%rbx),%%ymm2,%%ymm2		\n\t		vsubpd		 0x080(%%rbx),%%ymm10,%%ymm10	\n\t"\
		"vsubpd		 0x020(%%rbx),%%ymm3,%%ymm3		\n\t		vsubpd		 0x0a0(%%rbx),%%ymm11,%%ymm11	\n\t"\
		"vmovaps		  (%%rcx),%%ymm4				\n\t		vmovaps		 0x080(%%rcx),%%ymm12			\n\t"\
		"vmovaps	 0x020(%%rcx),%%ymm5				\n\t		vmovaps		 0x0a0(%%rcx),%%ymm13			\n\t"\
		"vmovaps		  (%%rcx),%%ymm6				\n\t		vmovaps		 0x080(%%rcx),%%ymm14			\n\t"\
		"vmovaps	 0x020(%%rcx),%%ymm7				\n\t		vmovaps		 0x0a0(%%rcx),%%ymm15			\n\t"\
		"vaddpd			  (%%rdx),%%ymm4,%%ymm4		\n\t		vaddpd		 0x080(%%rdx),%%ymm12,%%ymm12	\n\t"\
		"vaddpd		 0x020(%%rdx),%%ymm5,%%ymm5		\n\t		vaddpd		 0x0a0(%%rdx),%%ymm13,%%ymm13	\n\t"\
		"vsubpd			  (%%rdx),%%ymm6,%%ymm6		\n\t		vsubpd		 0x080(%%rdx),%%ymm14,%%ymm14	\n\t"\
		"vsubpd		 0x020(%%rdx),%%ymm7,%%ymm7		\n\t		vsubpd		 0x0a0(%%rdx),%%ymm15,%%ymm15	\n\t"\
		"vsubpd			%%ymm4,%%ymm0,%%ymm0		\n\t		vsubpd			%%ymm12,%%ymm8 ,%%ymm8 	\n\t"\
		"vsubpd			%%ymm5,%%ymm1,%%ymm1		\n\t		vsubpd			%%ymm13,%%ymm9 ,%%ymm9 	\n\t"\
		"vmovaps		%%ymm0,      (%%rbx)		\n\t		vmovaps			%%ymm8 , 0x080(%%rbx)			\n\t"\
		"vmovaps		%%ymm1, 0x020(%%rbx)		\n\t		vmovaps			%%ymm9 , 0x0a0(%%rbx)			\n\t"\
		"vaddpd			%%ymm4,%%ymm4,%%ymm4		\n\t		vaddpd			%%ymm12,%%ymm12,%%ymm12	\n\t"\
		"vaddpd			%%ymm5,%%ymm5,%%ymm5		\n\t		vaddpd			%%ymm13,%%ymm13,%%ymm13	\n\t"\
		"vaddpd			%%ymm0,%%ymm4,%%ymm4		\n\t		vaddpd			%%ymm8 ,%%ymm12,%%ymm12	\n\t"\
		"vaddpd			%%ymm1,%%ymm5,%%ymm5		\n\t		vaddpd			%%ymm9 ,%%ymm13,%%ymm13	\n\t"\
		"vmovaps		%%ymm4,      (%%rax)		\n\t		vmovaps			%%ymm12, 0x080(%%rax)			\n\t"\
		"vmovaps		%%ymm5, 0x020(%%rax)		\n\t		vmovaps			%%ymm13, 0x0a0(%%rax)			\n\t"\
		"vsubpd			%%ymm7,%%ymm2,%%ymm2		\n\t		vsubpd			%%ymm15,%%ymm10,%%ymm10	\n\t"\
		"vsubpd			%%ymm6,%%ymm3,%%ymm3		\n\t		vsubpd			%%ymm14,%%ymm11,%%ymm11	\n\t"\
		"vaddpd			%%ymm7,%%ymm7,%%ymm7		\n\t		vaddpd			%%ymm15,%%ymm15,%%ymm15	\n\t"\
		"vaddpd			%%ymm6,%%ymm6,%%ymm6		\n\t		vaddpd			%%ymm14,%%ymm14,%%ymm14	\n\t"\
		"vaddpd			%%ymm2,%%ymm7,%%ymm7		\n\t		vaddpd			%%ymm10,%%ymm15,%%ymm15	\n\t"\
		"vaddpd			%%ymm3,%%ymm6,%%ymm6		\n\t		vaddpd			%%ymm11,%%ymm14,%%ymm14	\n\t"\
		"vmovaps		%%ymm3,%%ymm0				\n\t		vmovaps			%%ymm11,%%ymm8 				\n\t"\
		"vmovaps		%%ymm6,%%ymm1				\n\t		vmovaps			%%ymm14,%%ymm9 				\n\t"\
		"vsubpd			%%ymm7,%%ymm3,%%ymm3		\n\t		vsubpd			%%ymm15,%%ymm11,%%ymm11	\n\t"\
		"vsubpd			%%ymm2,%%ymm6,%%ymm6		\n\t		vsubpd			%%ymm10,%%ymm14,%%ymm14	\n\t"\
		"vaddpd			%%ymm7,%%ymm0,%%ymm0		\n\t		vaddpd			%%ymm15,%%ymm8 ,%%ymm8 	\n\t"\
		"vaddpd			%%ymm2,%%ymm1,%%ymm1		\n\t		vaddpd			%%ymm10,%%ymm9 ,%%ymm9 	\n\t"\
		"vmulpd			  (%%rsi),%%ymm3,%%ymm3		\n\t		vmulpd			  (%%rsi),%%ymm11,%%ymm11	\n\t"\
		"vmulpd			  (%%rsi),%%ymm6,%%ymm6		\n\t		vmulpd			  (%%rsi),%%ymm14,%%ymm14	\n\t"\
		"vmulpd			  (%%rsi),%%ymm0,%%ymm0		\n\t		vmulpd			  (%%rsi),%%ymm8 ,%%ymm8 	\n\t"\
		"vmulpd			  (%%rsi),%%ymm1,%%ymm1		\n\t		vmulpd			  (%%rsi),%%ymm9 ,%%ymm9 	\n\t"\
		"vmovaps		%%ymm3, 0x020(%%rcx)		\n\t		vmovaps			%%ymm11, 0x0a0(%%rcx)			\n\t"\
		"vmovaps		%%ymm6, 0x020(%%rdx)		\n\t		vmovaps			%%ymm14, 0x0a0(%%rdx)			\n\t"\
		"vmovaps		%%ymm0,      (%%rcx)		\n\t		vmovaps			%%ymm8 , 0x080(%%rcx)			\n\t"\
		"vmovaps		%%ymm1,      (%%rdx)		\n\t		vmovaps			%%ymm9 , 0x080(%%rdx)			\n\t"\
		"/*****	SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS *****/\n\t		/*****	SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS *****/\n\t"\
		"/*        (r02,r12,r22,r32,r0A,r1A,r2A,r3A): */\n\t		/*        (r06,r16,r26,r36,r0E,r1E,r2E,r3E):  */\n\t"\
		"vmovaps		-0x100(%%rax),%%ymm0			\n\t		vmovaps		-0x080(%%rax),%%ymm8 			\n\t"\
		"vmovaps		-0x100(%%rbx),%%ymm4			\n\t		vmovaps		-0x080(%%rbx),%%ymm12			\n\t"\
		"vmovaps		-0x0e0(%%rax),%%ymm1			\n\t		vmovaps		-0x060(%%rax),%%ymm9 			\n\t"\
		"vmovaps		-0x0e0(%%rbx),%%ymm5			\n\t		vmovaps		-0x060(%%rbx),%%ymm13			\n\t"\
		"vmovaps			  (%%rax),%%ymm2			\n\t		vmovaps		 0x080(%%rax),%%ymm10			\n\t"\
		"vmovaps		 0x020(%%rbx),%%ymm7			\n\t		vmovaps		 0x0a0(%%rbx),%%ymm15			\n\t"\
		"vmovaps		 0x020(%%rax),%%ymm3			\n\t		vmovaps		 0x0a0(%%rax),%%ymm11			\n\t"\
		"vmovaps			  (%%rbx),%%ymm6			\n\t		vmovaps		 0x080(%%rbx),%%ymm14			\n\t"\
		"vsubpd			%%ymm2,%%ymm0,%%ymm0		\n\t		vsubpd		%%ymm10,%%ymm8 ,%%ymm8 		\n\t"\
		"vsubpd			%%ymm7,%%ymm4,%%ymm4		\n\t		vsubpd		%%ymm15,%%ymm12,%%ymm12		\n\t"\
		"vsubpd			%%ymm3,%%ymm1,%%ymm1		\n\t		vsubpd		%%ymm11,%%ymm9 ,%%ymm9 		\n\t"\
		"vsubpd			%%ymm6,%%ymm5,%%ymm5		\n\t		vsubpd		%%ymm14,%%ymm13,%%ymm13		\n\t"\
		"vaddpd			%%ymm2,%%ymm2,%%ymm2		\n\t		vaddpd		%%ymm10,%%ymm10,%%ymm10		\n\t"\
		"vaddpd			%%ymm7,%%ymm7,%%ymm7		\n\t		vaddpd		%%ymm15,%%ymm15,%%ymm15		\n\t"\
		"vaddpd			%%ymm3,%%ymm3,%%ymm3		\n\t		vaddpd		%%ymm11,%%ymm11,%%ymm11		\n\t"\
		"vaddpd			%%ymm6,%%ymm6,%%ymm6		\n\t		vaddpd		%%ymm14,%%ymm14,%%ymm14		\n\t"\
		"vaddpd			%%ymm0,%%ymm2,%%ymm2		\n\t		vaddpd		%%ymm8 ,%%ymm10,%%ymm10		\n\t"\
		"vaddpd			%%ymm4,%%ymm7,%%ymm7		\n\t		vaddpd		%%ymm12,%%ymm15,%%ymm15		\n\t"\
		"vaddpd			%%ymm1,%%ymm3,%%ymm3		\n\t		vaddpd		%%ymm9 ,%%ymm11,%%ymm11		\n\t"\
		"vaddpd			%%ymm5,%%ymm6,%%ymm6		\n\t		vaddpd		%%ymm13,%%ymm14,%%ymm14		\n\t"\
		"vmovaps		%%ymm0,      (%%rax)		\n\t		vmovaps		%%ymm8 , 0x080(%%rax)			\n\t"\
		"vmovaps		%%ymm4,      (%%rbx)		\n\t		vmovaps		%%ymm12, 0x080(%%rbx)			\n\t"\
		"vmovaps		%%ymm1, 0x020(%%rax)		\n\t		vmovaps		%%ymm9 , 0x0a0(%%rax)			\n\t"\
		"vmovaps		%%ymm5,-0x0e0(%%rbx)		\n\t		vmovaps		%%ymm13,-0x060(%%rbx)			\n\t"\
		"vmovaps		%%ymm2,-0x100(%%rax)		\n\t		vmovaps		%%ymm10,-0x080(%%rax)			\n\t"\
		"vmovaps		%%ymm7,-0x100(%%rbx)		\n\t		vmovaps		%%ymm15,-0x080(%%rbx)			\n\t"\
		"vmovaps		%%ymm3,-0x0e0(%%rax)		\n\t		vmovaps		%%ymm11,-0x060(%%rax)			\n\t"\
		"vmovaps		%%ymm6, 0x020(%%rbx)		\n\t		vmovaps		%%ymm14, 0x0a0(%%rbx)			\n\t"\
		"vmovaps		-0x100(%%rcx),%%ymm0			\n\t		vmovaps		-0x080(%%rcx),%%ymm8 			\n\t"\
		"vmovaps		-0x100(%%rdx),%%ymm4			\n\t		vmovaps		-0x080(%%rdx),%%ymm12			\n\t"\
		"vmovaps		-0x0e0(%%rcx),%%ymm1			\n\t		vmovaps		-0x060(%%rcx),%%ymm9 			\n\t"\
		"vmovaps		-0x0e0(%%rdx),%%ymm5			\n\t		vmovaps		-0x060(%%rdx),%%ymm13			\n\t"\
		"vmovaps			  (%%rcx),%%ymm2			\n\t		vmovaps		 0x080(%%rcx),%%ymm10			\n\t"\
		"vmovaps		 0x020(%%rdx),%%ymm7			\n\t		vmovaps		 0x0a0(%%rdx),%%ymm15			\n\t"\
		"vmovaps		 0x020(%%rcx),%%ymm3			\n\t		vmovaps		 0x0a0(%%rcx),%%ymm11			\n\t"\
		"vmovaps			  (%%rdx),%%ymm6			\n\t		vmovaps		 0x080(%%rdx),%%ymm14			\n\t"\
		"vsubpd			%%ymm2,%%ymm0,%%ymm0		\n\t		vsubpd		%%ymm10,%%ymm8 ,%%ymm8 		\n\t"\
		"vsubpd			%%ymm7,%%ymm4,%%ymm4		\n\t		vsubpd		%%ymm15,%%ymm12,%%ymm12		\n\t"\
		"vsubpd			%%ymm3,%%ymm1,%%ymm1		\n\t		vsubpd		%%ymm11,%%ymm9 ,%%ymm9 		\n\t"\
		"vsubpd			%%ymm6,%%ymm5,%%ymm5		\n\t		vsubpd		%%ymm14,%%ymm13,%%ymm13		\n\t"\
		"vaddpd			%%ymm2,%%ymm2,%%ymm2		\n\t		vaddpd		%%ymm10,%%ymm10,%%ymm10		\n\t"\
		"vaddpd			%%ymm7,%%ymm7,%%ymm7		\n\t		vaddpd		%%ymm15,%%ymm15,%%ymm15		\n\t"\
		"vaddpd			%%ymm3,%%ymm3,%%ymm3		\n\t		vaddpd		%%ymm11,%%ymm11,%%ymm11		\n\t"\
		"vaddpd			%%ymm6,%%ymm6,%%ymm6		\n\t		vaddpd		%%ymm14,%%ymm14,%%ymm14		\n\t"\
		"vaddpd			%%ymm0,%%ymm2,%%ymm2		\n\t		vaddpd		%%ymm8 ,%%ymm10,%%ymm10		\n\t"\
		"vaddpd			%%ymm4,%%ymm7,%%ymm7		\n\t		vaddpd		%%ymm12,%%ymm15,%%ymm15		\n\t"\
		"vaddpd			%%ymm1,%%ymm3,%%ymm3		\n\t		vaddpd		%%ymm9 ,%%ymm11,%%ymm11		\n\t"\
		"vaddpd			%%ymm5,%%ymm6,%%ymm6		\n\t		vaddpd		%%ymm13,%%ymm14,%%ymm14		\n\t"\
		"vmovaps		%%ymm0,      (%%rcx)		\n\t		vmovaps		%%ymm8 , 0x080(%%rcx)			\n\t"\
		"vmovaps		%%ymm4,      (%%rdx)		\n\t		vmovaps		%%ymm12, 0x080(%%rdx)			\n\t"\
		"vmovaps		%%ymm1, 0x020(%%rcx)		\n\t		vmovaps		%%ymm9 , 0x0a0(%%rcx)			\n\t"\
		"vmovaps		%%ymm5,-0x0e0(%%rdx)		\n\t		vmovaps		%%ymm13,-0x060(%%rdx)			\n\t"\
		"vmovaps		%%ymm2,-0x100(%%rcx)		\n\t		vmovaps		%%ymm10,-0x080(%%rcx)			\n\t"\
		"vmovaps		%%ymm7,-0x100(%%rdx)		\n\t		vmovaps		%%ymm15,-0x080(%%rdx)			\n\t"\
		"vmovaps		%%ymm3,-0x0e0(%%rcx)		\n\t		vmovaps		%%ymm11,-0x060(%%rcx)			\n\t"\
		"vmovaps		%%ymm6, 0x020(%%rdx)		\n\t		vmovaps		%%ymm14, 0x0a0(%%rdx)			\n\t"\
		"/***************************************************************************************/\n\t"\
		"/* Now do eight more radix-4 transforms, including the internal and external twiddles: */\n\t"\
		"/***************************************************************************************/\n\t"\
		/* Using the upper block(s) of the main array for temp-storage in the section below led to a nasty AVX bug to track down: */\
		/* In fermat-mod mode the 4 block addresses in ascending order are add0,1,2,3 with no 'gaps' between blocks, whereas for */\
		/* mersenne-mod the addresses in asc. order are add0,2,3,1 with a gap between contiguous-data-block pairs 0,2 and 3,1. Thus */\
		/* for fermat-mod we need [add2] as the base-address of the 'high-half' block for temp-storage; for mersenne-mod we need [add3]. */\
		/* In both cases we have that add2 < add3 so instead use (add2 - add1): > 0 for fermat-mod, < 0 for mersenne - to differentiate: */\
		"movq	%[__add2],%%rsi		\n\t"/* destroyable copy of add2 */\
		"movq	%[__add2],%%rbx		\n\t"\
		"subq	%[__add1],%%rsi		\n\t"/* rsi = (add2 - add1); if this yields a borrow (i.e. sets CF) it's mersenne, else fermat. */\
		"cmovcq %[__add3],%%rbx	\n\t" /* if CF set (i.e. h > l), copy source [add3] into dest (rbx), else leave dest = [add2]. */\
		"\n\t"\
		"/************************************************/	\n\t		/************************************************/	\n\t"\
		"/*...Block 2: t02,t12,t22,t32 -> r10,14,12,16:  */	\n\t		/*...Block 6: t0A,t1A,t2A,t3A -> r18,1C,1A,1E:  */	\n\t"\
		"/************************************************/	\n\t		/************************************************/	\n\t"\
		"movq		%[__isrt2],%%rsi					\n\t"\
		"movq		%[__r10],%%rax	/* base-addr in rcol = c05/r18, so rax/rdi offset +0x100 vs lcol */\n\t"\
		"movq		%%rsi,%%rcx					\n\t"\
		"movq		%%rsi,%%rdx					\n\t"\
		"movq		%[__c01],%%rdi					\n\t"\
		"addq		$0x020,%%rcx	/* cc0 */		\n\t"\
		"addq		$0x060,%%rdx	/* cc1 */		\n\t"\
		"vmovaps		 0x040(%%rax),%%ymm4				\n\t		vmovaps		 0x140(%%rax),%%ymm12				\n\t"\
		"vmovaps		 0x0c0(%%rax),%%ymm0				\n\t		vmovaps		 0x1c0(%%rax),%%ymm8 				\n\t"\
		"vmovaps		 0x060(%%rax),%%ymm5				\n\t		vmovaps		 0x160(%%rax),%%ymm13				\n\t"\
		"vmovaps		 0x0e0(%%rax),%%ymm1				\n\t		vmovaps		 0x1e0(%%rax),%%ymm9 				\n\t"\
		"vmovaps		 0x040(%%rax),%%ymm6				\n\t		vmovaps		 0x140(%%rax),%%ymm14				\n\t"\
		"vmovaps		 0x0c0(%%rax),%%ymm2				\n\t		vmovaps		 0x1c0(%%rax),%%ymm10				\n\t"\
		"vmovaps		 0x060(%%rax),%%ymm7				\n\t		vmovaps		 0x160(%%rax),%%ymm15				\n\t"\
		"vmovaps		 0x0e0(%%rax),%%ymm3				\n\t		vmovaps		 0x1e0(%%rax),%%ymm11				\n\t"\
		"vmulpd			  (%%rdx),%%ymm4,%%ymm4			\n\t		vmulpd		 0x060(%%rdx),%%ymm12,%%ymm12			\n\t"\
		"vmulpd		 0x040(%%rdx),%%ymm0,%%ymm0			\n\t		vmulpd			  (%%rdx),%%ymm8 ,%%ymm8 			\n\t"\
		"vmulpd			  (%%rdx),%%ymm5,%%ymm5			\n\t		vmulpd		 0x060(%%rdx),%%ymm13,%%ymm13			\n\t"\
		"vmulpd		 0x040(%%rdx),%%ymm1,%%ymm1			\n\t		vmulpd			  (%%rdx),%%ymm9 ,%%ymm9 			\n\t"\
		"vmulpd		 0x020(%%rdx),%%ymm6,%%ymm6			\n\t		vmulpd		 0x040(%%rdx),%%ymm14,%%ymm14			\n\t"\
		"vmulpd		 0x060(%%rdx),%%ymm2,%%ymm2			\n\t		vmulpd		 0x020(%%rdx),%%ymm10,%%ymm10			\n\t"\
		"vmulpd		 0x020(%%rdx),%%ymm7,%%ymm7			\n\t		vmulpd		 0x040(%%rdx),%%ymm15,%%ymm15			\n\t"\
		"vmulpd		 0x060(%%rdx),%%ymm3,%%ymm3			\n\t		vmulpd		 0x020(%%rdx),%%ymm11,%%ymm11			\n\t"\
		"vsubpd		%%ymm6,%%ymm5,%%ymm5			\n\t		vsubpd		%%ymm14,%%ymm13,%%ymm13			\n\t"\
		"vsubpd		%%ymm2,%%ymm1,%%ymm1			\n\t		vaddpd		%%ymm10,%%ymm9 ,%%ymm9 			\n\t"\
		"vaddpd		%%ymm7,%%ymm4,%%ymm4			\n\t		vaddpd		%%ymm15,%%ymm12,%%ymm12			\n\t"\
		"vaddpd		%%ymm3,%%ymm0,%%ymm0			\n\t		vsubpd		%%ymm11,%%ymm8 ,%%ymm8 			\n\t"\
		"vmovaps		%%ymm5,%%ymm7					\n\t		vmovaps		%%ymm13,%%ymm15					\n\t"\
		"vmovaps		%%ymm4,%%ymm6					\n\t		vmovaps		%%ymm12,%%ymm14					\n\t"\
		"vsubpd		%%ymm0,%%ymm4,%%ymm4			\n\t		vaddpd		%%ymm8 ,%%ymm12,%%ymm12			\n\t"\
		"vsubpd		%%ymm1,%%ymm5,%%ymm5			\n\t		vaddpd		%%ymm9 ,%%ymm13,%%ymm13			\n\t"\
		"vaddpd		%%ymm0,%%ymm6,%%ymm6			\n\t		vsubpd		%%ymm8 ,%%ymm14,%%ymm14			\n\t"\
		"vaddpd		%%ymm1,%%ymm7,%%ymm7			\n\t		vsubpd		%%ymm9 ,%%ymm15,%%ymm15			\n\t"\
		"vmovaps		 0x080(%%rax),%%ymm2				\n\t		vmovaps		 0x180(%%rax),%%ymm10				\n\t"\
		"vmovaps		 0x0a0(%%rax),%%ymm3				\n\t		vmovaps		 0x1a0(%%rax),%%ymm11				\n\t"\
		"vmovaps		 0x080(%%rax),%%ymm0				\n\t		vmovaps		 0x180(%%rax),%%ymm8 				\n\t"\
		"vmovaps		 0x0a0(%%rax),%%ymm1				\n\t		vmovaps		 0x1a0(%%rax),%%ymm9 				\n\t"\
		"vmulpd			  (%%rcx),%%ymm2,%%ymm2			\n\t		vmulpd		 0x020(%%rcx),%%ymm10,%%ymm10			\n\t"\
		"vmulpd		 0x020(%%rcx),%%ymm1,%%ymm1			\n\t		vmulpd			  (%%rcx),%%ymm9 ,%%ymm9 			\n\t"\
		"vmulpd			  (%%rcx),%%ymm3,%%ymm3			\n\t		vmulpd		 0x020(%%rcx),%%ymm11,%%ymm11			\n\t"\
		"vmulpd		 0x020(%%rcx),%%ymm0,%%ymm0			\n\t		vmulpd			  (%%rcx),%%ymm8 ,%%ymm8 			\n\t"\
		"vaddpd		%%ymm1,%%ymm2,%%ymm2			\n\t		vsubpd		%%ymm9 ,%%ymm10,%%ymm10			\n\t"\
		"vsubpd		%%ymm0,%%ymm3,%%ymm3			\n\t		vaddpd		%%ymm8 ,%%ymm11,%%ymm11			\n\t"\
		"vmovaps			  (%%rax),%%ymm0				\n\t		vmovaps		 0x100(%%rax),%%ymm8 				\n\t"\
		"vmovaps		 0x020(%%rax),%%ymm1				\n\t		vmovaps		 0x120(%%rax),%%ymm9 				\n\t"\
		"vsubpd		%%ymm2,%%ymm0,%%ymm0			\n\t		vsubpd		%%ymm10,%%ymm8 ,%%ymm8 			\n\t"\
		"vsubpd		%%ymm3,%%ymm1,%%ymm1			\n\t		vsubpd		%%ymm11,%%ymm9 ,%%ymm9 			\n\t"\
		"vaddpd		%%ymm2,%%ymm2,%%ymm2			\n\t		vaddpd		%%ymm10,%%ymm10,%%ymm10			\n\t"\
		"vaddpd		%%ymm3,%%ymm3,%%ymm3			\n\t		vaddpd		%%ymm11,%%ymm11,%%ymm11			\n\t"\
		"vaddpd		%%ymm0,%%ymm2,%%ymm2			\n\t		vaddpd		%%ymm8 ,%%ymm10,%%ymm10			\n\t"\
		"vaddpd		%%ymm1,%%ymm3,%%ymm3			\n\t		vaddpd		%%ymm9 ,%%ymm11,%%ymm11			\n\t"\
		"vsubpd		%%ymm6,%%ymm2,%%ymm2			\n\t		vsubpd		%%ymm14,%%ymm8 ,%%ymm8 			\n\t"\
		"vsubpd		%%ymm7,%%ymm3,%%ymm3			\n\t		vsubpd		%%ymm15,%%ymm9 ,%%ymm9 			\n\t"\
		"vaddpd		%%ymm6,%%ymm6,%%ymm6			\n\t		vaddpd		%%ymm14,%%ymm14,%%ymm14			\n\t"\
		"vaddpd		%%ymm7,%%ymm7,%%ymm7			\n\t		vaddpd		%%ymm15,%%ymm15,%%ymm15			\n\t"\
		"vaddpd		%%ymm2,%%ymm6,%%ymm6			\n\t		vaddpd		%%ymm8 ,%%ymm14,%%ymm14			\n\t"\
		"vaddpd		%%ymm3,%%ymm7,%%ymm7			\n\t		vaddpd		%%ymm9 ,%%ymm15,%%ymm15			\n\t"\
		"vmovaps		%%ymm2, 0x040(%%rax)			\n\t		vmovaps		%%ymm8 , 0x140(%%rax)			\n\t"\
		"vmovaps		%%ymm3, 0x060(%%rax)			\n\t		vmovaps		%%ymm9 , 0x160(%%rax)			\n\t"\
		"vmovaps		%%ymm6,%%ymm2					\n\t		vmovaps		%%ymm14,%%ymm8 					\n\t"\
		"vmovaps		%%ymm7,%%ymm3					\n\t		vmovaps		%%ymm15,%%ymm9 					\n\t"\
		"vmulpd			  (%%rdi),%%ymm6,%%ymm6			\n\t		vmulpd		 0x100(%%rdi),%%ymm14,%%ymm14			\n\t"\
		"vmulpd			  (%%rdi),%%ymm7,%%ymm7			\n\t		vmulpd		 0x100(%%rdi),%%ymm15,%%ymm15			\n\t"\
		"vmulpd		 0x020(%%rdi),%%ymm2,%%ymm2			\n\t		vmulpd		 0x120(%%rdi),%%ymm8 ,%%ymm8 			\n\t"\
		"vmulpd		 0x020(%%rdi),%%ymm3,%%ymm3			\n\t		vmulpd		 0x120(%%rdi),%%ymm9 ,%%ymm9 			\n\t"\
		"vsubpd		%%ymm2,%%ymm7,%%ymm7			\n\t		vsubpd		%%ymm8 ,%%ymm15,%%ymm15			\n\t"\
		"vaddpd		%%ymm3,%%ymm6,%%ymm6			\n\t		vaddpd		%%ymm9 ,%%ymm14,%%ymm14			\n\t"\
		"vmovaps		%%ymm7, 0x020(%%rbx)			\n\t		vmovaps		%%ymm15, 0x0a0(%%rbx)			\n\t"\
		"vmovaps		%%ymm6,      (%%rbx)			\n\t		vmovaps		%%ymm14, 0x080(%%rbx)			\n\t"\
		"vmovaps		 0x040(%%rax),%%ymm6				\n\t		vmovaps		 0x140(%%rax),%%ymm14				\n\t"\
		"vmovaps		 0x060(%%rax),%%ymm7				\n\t		vmovaps		 0x160(%%rax),%%ymm15				\n\t"\
		"vmovaps		%%ymm6,%%ymm2					\n\t		vmovaps		%%ymm14,%%ymm8 					\n\t"\
		"vmovaps		%%ymm7,%%ymm3					\n\t		vmovaps		%%ymm15,%%ymm9 					\n\t"\
		"vmulpd		 0x040(%%rdi),%%ymm6,%%ymm6			\n\t		vmulpd		 0x140(%%rdi),%%ymm14,%%ymm14			\n\t"\
		"vmulpd		 0x040(%%rdi),%%ymm7,%%ymm7			\n\t		vmulpd		 0x140(%%rdi),%%ymm15,%%ymm15			\n\t"\
		"vmulpd		 0x060(%%rdi),%%ymm2,%%ymm2			\n\t		vmulpd		 0x160(%%rdi),%%ymm8 ,%%ymm8 			\n\t"\
		"vmulpd		 0x060(%%rdi),%%ymm3,%%ymm3			\n\t		vmulpd		 0x160(%%rdi),%%ymm9 ,%%ymm9 			\n\t"\
		"vsubpd		%%ymm2,%%ymm7,%%ymm7			\n\t		vsubpd		%%ymm8 ,%%ymm15,%%ymm15			\n\t"\
		"vaddpd		%%ymm3,%%ymm6,%%ymm6			\n\t		vaddpd		%%ymm9 ,%%ymm14,%%ymm14			\n\t"\
		"vmovaps		%%ymm7, 0x220(%%rbx)			\n\t		vmovaps		%%ymm15, 0x2a0(%%rbx)			\n\t"\
		"vmovaps		%%ymm6, 0x200(%%rbx)			\n\t		vmovaps		%%ymm14, 0x280(%%rbx)			\n\t"\
		"addq		$0x080,%%rdi					\n\t"\
		"vsubpd		%%ymm5,%%ymm0,%%ymm0			\n\t		vsubpd		%%ymm13,%%ymm10,%%ymm10			\n\t"\
		"vsubpd		%%ymm4,%%ymm1,%%ymm1			\n\t		vsubpd		%%ymm12,%%ymm11,%%ymm11			\n\t"\
		"vaddpd		%%ymm5,%%ymm5,%%ymm5			\n\t		vaddpd		%%ymm13,%%ymm13,%%ymm13			\n\t"\
		"vaddpd		%%ymm4,%%ymm4,%%ymm4			\n\t		vaddpd		%%ymm12,%%ymm12,%%ymm12			\n\t"\
		"vaddpd		%%ymm0,%%ymm5,%%ymm5			\n\t		vaddpd		%%ymm10,%%ymm13,%%ymm13			\n\t"\
		"vaddpd		%%ymm1,%%ymm4,%%ymm4			\n\t		vaddpd		%%ymm11,%%ymm12,%%ymm12			\n\t"\
		"vmovaps		%%ymm5,%%ymm2					\n\t		vmovaps		%%ymm13,%%ymm8 					\n\t"\
		"vmovaps		%%ymm1,%%ymm3					\n\t		vmovaps		%%ymm11,%%ymm9 					\n\t"\
		"vmulpd			  (%%rdi),%%ymm5,%%ymm5			\n\t		vmulpd		 0x100(%%rdi),%%ymm13,%%ymm13			\n\t"\
		"vmulpd			  (%%rdi),%%ymm1,%%ymm1			\n\t		vmulpd		 0x100(%%rdi),%%ymm11,%%ymm11			\n\t"\
		"vmulpd		 0x020(%%rdi),%%ymm2,%%ymm2			\n\t		vmulpd		 0x120(%%rdi),%%ymm8 ,%%ymm8 			\n\t"\
		"vmulpd		 0x020(%%rdi),%%ymm3,%%ymm3			\n\t		vmulpd		 0x120(%%rdi),%%ymm9 ,%%ymm9 			\n\t"\
		"vsubpd		%%ymm2,%%ymm1,%%ymm1			\n\t		vsubpd		%%ymm8 ,%%ymm11,%%ymm11			\n\t"\
		"vaddpd		%%ymm3,%%ymm5,%%ymm5			\n\t		vaddpd		%%ymm9 ,%%ymm13,%%ymm13			\n\t"\
		"vmovaps		%%ymm1, 0x120(%%rbx)			\n\t		vmovaps		%%ymm11, 0x1a0(%%rbx)			\n\t"\
		"vmovaps		%%ymm5, 0x100(%%rbx)			\n\t		vmovaps		%%ymm13, 0x180(%%rbx)			\n\t"\
		"vmovaps		%%ymm0,%%ymm2					\n\t		vmovaps		%%ymm10,%%ymm8 					\n\t"\
		"vmovaps		%%ymm4,%%ymm3					\n\t		vmovaps		%%ymm12,%%ymm9 					\n\t"\
		"vmulpd		 0x040(%%rdi),%%ymm0,%%ymm0			\n\t		vmulpd		 0x140(%%rdi),%%ymm10,%%ymm10			\n\t"\
		"vmulpd		 0x040(%%rdi),%%ymm4,%%ymm4			\n\t		vmulpd		 0x140(%%rdi),%%ymm12,%%ymm12			\n\t"\
		"vmulpd		 0x060(%%rdi),%%ymm2,%%ymm2			\n\t		vmulpd		 0x160(%%rdi),%%ymm8 ,%%ymm8 			\n\t"\
		"vmulpd		 0x060(%%rdi),%%ymm3,%%ymm3			\n\t		vmulpd		 0x160(%%rdi),%%ymm9 ,%%ymm9 			\n\t"\
		"vsubpd		%%ymm2,%%ymm4,%%ymm4			\n\t		vsubpd		%%ymm8 ,%%ymm12,%%ymm12			\n\t"\
		"vaddpd		%%ymm3,%%ymm0,%%ymm0			\n\t		vaddpd		%%ymm9 ,%%ymm10,%%ymm10			\n\t"\
		"vmovaps		%%ymm4, 0x320(%%rbx)			\n\t		vmovaps		%%ymm12, 0x3a0(%%rbx)			\n\t"\
		"vmovaps		%%ymm0, 0x300(%%rbx)			\n\t		vmovaps		%%ymm10, 0x380(%%rbx)			\n\t"\
		"\n\t"\
		"/************************************************/	\n\t		/************************************************/	\n\t"\
		"/*...Block 4: t06,t16,t26,t36 -> r30,34,32,36:  */	\n\t		/*...Block 8: t0E,t1E,t2E,t3E -> r38,3C,3A,3E:  */	\n\t"\
		"/************************************************/	\n\t		/************************************************/	\n\t"\
		"addq		$0x400,%%rax						\n\t		addq		$0x100,%%rdi						\n\t"\
		"vmovaps		 0x040(%%rax),%%ymm4				\n\t		vmovaps		 0x140(%%rax),%%ymm12				\n\t"\
		"vmovaps		 0x0c0(%%rax),%%ymm0				\n\t		vmovaps		 0x1c0(%%rax),%%ymm8 				\n\t"\
		"vmovaps		 0x060(%%rax),%%ymm5				\n\t		vmovaps		 0x160(%%rax),%%ymm13				\n\t"\
		"vmovaps		 0x0e0(%%rax),%%ymm1				\n\t		vmovaps		 0x1e0(%%rax),%%ymm9 				\n\t"\
		"vmovaps		 0x040(%%rax),%%ymm6				\n\t		vmovaps		 0x140(%%rax),%%ymm14				\n\t"\
		"vmovaps		 0x0c0(%%rax),%%ymm2				\n\t		vmovaps		 0x1c0(%%rax),%%ymm10				\n\t"\
		"vmovaps		 0x060(%%rax),%%ymm7				\n\t		vmovaps		 0x160(%%rax),%%ymm15				\n\t"\
		"vmovaps		 0x0e0(%%rax),%%ymm3				\n\t		vmovaps		 0x1e0(%%rax),%%ymm11				\n\t"\
		"vmulpd		 0x040(%%rdx),%%ymm4,%%ymm4			\n\t		vmulpd		 0x020(%%rdx),%%ymm12,%%ymm12			\n\t"\
		"vmulpd		 0x020(%%rdx),%%ymm0,%%ymm0			\n\t		vmulpd		 0x060(%%rdx),%%ymm8 ,%%ymm8 			\n\t"\
		"vmulpd		 0x040(%%rdx),%%ymm5,%%ymm5			\n\t		vmulpd		 0x020(%%rdx),%%ymm13,%%ymm13			\n\t"\
		"vmulpd		 0x020(%%rdx),%%ymm1,%%ymm1			\n\t		vmulpd		 0x060(%%rdx),%%ymm9 ,%%ymm9 			\n\t"\
		"vmulpd		 0x060(%%rdx),%%ymm6,%%ymm6			\n\t		vmulpd			  (%%rdx),%%ymm14,%%ymm14			\n\t"\
		"vmulpd			  (%%rdx),%%ymm2,%%ymm2			\n\t		vmulpd		 0x040(%%rdx),%%ymm10,%%ymm10			\n\t"\
		"vmulpd		 0x060(%%rdx),%%ymm7,%%ymm7			\n\t		vmulpd			  (%%rdx),%%ymm15,%%ymm15			\n\t"\
		"vmulpd			  (%%rdx),%%ymm3,%%ymm3			\n\t		vmulpd		 0x040(%%rdx),%%ymm11,%%ymm11			\n\t"\
		"vsubpd		%%ymm6,%%ymm5,%%ymm5			\n\t		vsubpd		%%ymm14,%%ymm13,%%ymm13			\n\t"\
		"vaddpd		%%ymm2,%%ymm1,%%ymm1			\n\t		vsubpd		%%ymm10,%%ymm9 ,%%ymm9 			\n\t"\
		"vaddpd		%%ymm7,%%ymm4,%%ymm4			\n\t		vaddpd		%%ymm15,%%ymm12,%%ymm12			\n\t"\
		"vsubpd		%%ymm3,%%ymm0,%%ymm0			\n\t		vaddpd		%%ymm11,%%ymm8 ,%%ymm8 			\n\t"\
		"vmovaps		%%ymm5,%%ymm7					\n\t		vmovaps		%%ymm13,%%ymm15					\n\t"\
		"vmovaps		%%ymm4,%%ymm6					\n\t		vmovaps		%%ymm12,%%ymm14					\n\t"\
		"vaddpd		%%ymm0,%%ymm4,%%ymm4			\n\t		vaddpd		%%ymm8 ,%%ymm12,%%ymm12			\n\t"\
		"vaddpd		%%ymm1,%%ymm5,%%ymm5			\n\t		vaddpd		%%ymm9 ,%%ymm13,%%ymm13			\n\t"\
		"vsubpd		%%ymm0,%%ymm6,%%ymm6			\n\t		vsubpd		%%ymm8 ,%%ymm14,%%ymm14			\n\t"\
		"vsubpd		%%ymm1,%%ymm7,%%ymm7			\n\t		vsubpd		%%ymm9 ,%%ymm15,%%ymm15			\n\t"\
		"vmovaps		 0x080(%%rax),%%ymm2				\n\t		vmovaps		 0x180(%%rax),%%ymm10				\n\t"\
		"vmovaps		 0x0a0(%%rax),%%ymm3				\n\t		vmovaps		 0x1a0(%%rax),%%ymm11				\n\t"\
		"vmovaps		 0x080(%%rax),%%ymm0				\n\t		vmovaps		 0x180(%%rax),%%ymm8 				\n\t"\
		"vmovaps		 0x0a0(%%rax),%%ymm1				\n\t		vmovaps		 0x1a0(%%rax),%%ymm9 				\n\t"\
		"vmulpd		 0x020(%%rcx),%%ymm2,%%ymm2			\n\t		vmulpd			  (%%rcx),%%ymm10,%%ymm10			\n\t"\
		"vmulpd			  (%%rcx),%%ymm1,%%ymm1			\n\t		vmulpd		 0x020(%%rcx),%%ymm9 ,%%ymm9 			\n\t"\
		"vmulpd		 0x020(%%rcx),%%ymm3,%%ymm3			\n\t		vmulpd			  (%%rcx),%%ymm11,%%ymm11			\n\t"\
		"vmulpd			  (%%rcx),%%ymm0,%%ymm0			\n\t		vmulpd		 0x020(%%rcx),%%ymm8 ,%%ymm8 			\n\t"\
		"vaddpd		%%ymm1,%%ymm2,%%ymm2			\n\t		vsubpd		%%ymm9 ,%%ymm10,%%ymm10			\n\t"\
		"vsubpd		%%ymm0,%%ymm3,%%ymm3			\n\t		vaddpd		%%ymm8 ,%%ymm11,%%ymm11			\n\t"\
		"vmovaps			  (%%rax),%%ymm0				\n\t		vmovaps		 0x100(%%rax),%%ymm8 				\n\t"\
		"vmovaps		 0x020(%%rax),%%ymm1				\n\t		vmovaps		 0x120(%%rax),%%ymm9 				\n\t"\
		"vsubpd		%%ymm2,%%ymm0,%%ymm0			\n\t		vsubpd		%%ymm10,%%ymm8 ,%%ymm8 			\n\t"\
		"vsubpd		%%ymm3,%%ymm1,%%ymm1			\n\t		vsubpd		%%ymm11,%%ymm9 ,%%ymm9 			\n\t"\
		"vaddpd		%%ymm2,%%ymm2,%%ymm2			\n\t		vaddpd		%%ymm10,%%ymm10,%%ymm10			\n\t"\
		"vaddpd		%%ymm3,%%ymm3,%%ymm3			\n\t		vaddpd		%%ymm11,%%ymm11,%%ymm11			\n\t"\
		"vaddpd		%%ymm0,%%ymm2,%%ymm2			\n\t		vaddpd		%%ymm8 ,%%ymm10,%%ymm10			\n\t"\
		"vaddpd		%%ymm1,%%ymm3,%%ymm3			\n\t		vaddpd		%%ymm9 ,%%ymm11,%%ymm11			\n\t"\
		"addq		$0x080,%%rdi						\n\t"\
		"vsubpd		%%ymm6,%%ymm2,%%ymm2			\n\t		vsubpd		%%ymm14,%%ymm8 ,%%ymm8 			\n\t"\
		"vsubpd		%%ymm7,%%ymm3,%%ymm3			\n\t		vsubpd		%%ymm15,%%ymm9 ,%%ymm9 			\n\t"\
		"vaddpd		%%ymm6,%%ymm6,%%ymm6			\n\t		vaddpd		%%ymm14,%%ymm14,%%ymm14			\n\t"\
		"vaddpd		%%ymm7,%%ymm7,%%ymm7			\n\t		vaddpd		%%ymm15,%%ymm15,%%ymm15			\n\t"\
		"vaddpd		%%ymm2,%%ymm6,%%ymm6			\n\t		vaddpd		%%ymm8 ,%%ymm14,%%ymm14			\n\t"\
		"vaddpd		%%ymm3,%%ymm7,%%ymm7			\n\t		vaddpd		%%ymm9 ,%%ymm15,%%ymm15			\n\t"\
		"vmovaps		%%ymm2, 0x040(%%rax)			\n\t		vmovaps		%%ymm8 , 0x140(%%rax)			\n\t"\
		"vmovaps		%%ymm3, 0x060(%%rax)			\n\t		vmovaps		%%ymm9 , 0x160(%%rax)			\n\t"\
		"vmovaps		%%ymm6,%%ymm2					\n\t		vmovaps		%%ymm14,%%ymm8 					\n\t"\
		"vmovaps		%%ymm7,%%ymm3					\n\t		vmovaps		%%ymm15,%%ymm9 					\n\t"\
		"vmulpd			  (%%rdi),%%ymm6,%%ymm6			\n\t		vmulpd		 0x100(%%rdi),%%ymm14,%%ymm14			\n\t"\
		"vmulpd			  (%%rdi),%%ymm7,%%ymm7			\n\t		vmulpd		 0x100(%%rdi),%%ymm15,%%ymm15			\n\t"\
		"vmulpd		 0x020(%%rdi),%%ymm2,%%ymm2			\n\t		vmulpd		 0x120(%%rdi),%%ymm8 ,%%ymm8 			\n\t"\
		"vmulpd		 0x020(%%rdi),%%ymm3,%%ymm3			\n\t		vmulpd		 0x120(%%rdi),%%ymm9 ,%%ymm9 			\n\t"\
		"vsubpd		%%ymm2,%%ymm7,%%ymm7			\n\t		vsubpd		%%ymm8 ,%%ymm15,%%ymm15			\n\t"\
		"vaddpd		%%ymm3,%%ymm6,%%ymm6			\n\t		vaddpd		%%ymm9 ,%%ymm14,%%ymm14			\n\t"\
		"vmovaps		%%ymm7, 0x060(%%rbx)			\n\t		vmovaps		%%ymm15, 0x0e0(%%rbx)			\n\t"\
		"vmovaps		%%ymm6, 0x040(%%rbx)			\n\t		vmovaps		%%ymm14, 0x0c0(%%rbx)			\n\t"\
		"vmovaps		 0x040(%%rax),%%ymm6				\n\t		vmovaps		 0x140(%%rax),%%ymm14				\n\t"\
		"vmovaps		 0x060(%%rax),%%ymm7				\n\t		vmovaps		 0x160(%%rax),%%ymm15				\n\t"\
		"vmovaps		%%ymm6,%%ymm2					\n\t		vmovaps		%%ymm14,%%ymm8 					\n\t"\
		"vmovaps		%%ymm7,%%ymm3					\n\t		vmovaps		%%ymm15,%%ymm9 					\n\t"\
		"vmulpd		 0x040(%%rdi),%%ymm6,%%ymm6			\n\t		vmulpd		 0x140(%%rdi),%%ymm14,%%ymm14			\n\t"\
		"vmulpd		 0x040(%%rdi),%%ymm7,%%ymm7			\n\t		vmulpd		 0x140(%%rdi),%%ymm15,%%ymm15			\n\t"\
		"vmulpd		 0x060(%%rdi),%%ymm2,%%ymm2			\n\t		vmulpd		 0x160(%%rdi),%%ymm8 ,%%ymm8 			\n\t"\
		"vmulpd		 0x060(%%rdi),%%ymm3,%%ymm3			\n\t		vmulpd		 0x160(%%rdi),%%ymm9 ,%%ymm9 			\n\t"\
		"vsubpd		%%ymm2,%%ymm7,%%ymm7			\n\t		vsubpd		%%ymm8 ,%%ymm15,%%ymm15			\n\t"\
		"vaddpd		%%ymm3,%%ymm6,%%ymm6			\n\t		vaddpd		%%ymm9 ,%%ymm14,%%ymm14			\n\t"\
		"vmovaps		%%ymm7, 0x260(%%rbx)			\n\t		vmovaps		%%ymm15, 0x2e0(%%rbx)			\n\t"\
		"vmovaps		%%ymm6, 0x240(%%rbx)			\n\t		vmovaps		%%ymm14, 0x2c0(%%rbx)			\n\t"\
		"addq		$0x080,%%rdi						\n\t"\
		"vsubpd		%%ymm5,%%ymm0,%%ymm0			\n\t		vsubpd		%%ymm13,%%ymm10,%%ymm10			\n\t"\
		"vsubpd		%%ymm4,%%ymm1,%%ymm1			\n\t		vsubpd		%%ymm12,%%ymm11,%%ymm11			\n\t"\
		"vaddpd		%%ymm5,%%ymm5,%%ymm5			\n\t		vaddpd		%%ymm13,%%ymm13,%%ymm13			\n\t"\
		"vaddpd		%%ymm4,%%ymm4,%%ymm4			\n\t		vaddpd		%%ymm12,%%ymm12,%%ymm12			\n\t"\
		"vaddpd		%%ymm0,%%ymm5,%%ymm5			\n\t		vaddpd		%%ymm10,%%ymm13,%%ymm13			\n\t"\
		"vaddpd		%%ymm1,%%ymm4,%%ymm4			\n\t		vaddpd		%%ymm11,%%ymm12,%%ymm12			\n\t"\
		"vmovaps		%%ymm5,%%ymm2					\n\t		vmovaps		%%ymm13,%%ymm8 					\n\t"\
		"vmovaps		%%ymm1,%%ymm3					\n\t		vmovaps		%%ymm11,%%ymm9 					\n\t"\
		"vmulpd			  (%%rdi),%%ymm5,%%ymm5			\n\t		vmulpd		 0x100(%%rdi),%%ymm13,%%ymm13			\n\t"\
		"vmulpd			  (%%rdi),%%ymm1,%%ymm1			\n\t		vmulpd		 0x100(%%rdi),%%ymm11,%%ymm11			\n\t"\
		"vmulpd		 0x020(%%rdi),%%ymm2,%%ymm2			\n\t		vmulpd		 0x120(%%rdi),%%ymm8 ,%%ymm8 			\n\t"\
		"vmulpd		 0x020(%%rdi),%%ymm3,%%ymm3			\n\t		vmulpd		 0x120(%%rdi),%%ymm9 ,%%ymm9 			\n\t"\
		"vsubpd		%%ymm2,%%ymm1,%%ymm1			\n\t		vsubpd		%%ymm8 ,%%ymm11,%%ymm11			\n\t"\
		"vaddpd		%%ymm3,%%ymm5,%%ymm5			\n\t		vaddpd		%%ymm9 ,%%ymm13,%%ymm13			\n\t"\
		"vmovaps		%%ymm1, 0x160(%%rbx)			\n\t		vmovaps		%%ymm11, 0x1e0(%%rbx)			\n\t"\
		"vmovaps		%%ymm5, 0x140(%%rbx)			\n\t		vmovaps		%%ymm13, 0x1c0(%%rbx)			\n\t"\
		"vmovaps		%%ymm0,%%ymm2					\n\t		vmovaps		%%ymm10,%%ymm8 					\n\t"\
		"vmovaps		%%ymm4,%%ymm3					\n\t		vmovaps		%%ymm12,%%ymm9 					\n\t"\
		"vmulpd		 0x040(%%rdi),%%ymm0,%%ymm0			\n\t		vmulpd		 0x140(%%rdi),%%ymm10,%%ymm10			\n\t"\
		"vmulpd		 0x040(%%rdi),%%ymm4,%%ymm4			\n\t		vmulpd		 0x140(%%rdi),%%ymm12,%%ymm12			\n\t"\
		"vmulpd		 0x060(%%rdi),%%ymm2,%%ymm2			\n\t		vmulpd		 0x160(%%rdi),%%ymm8 ,%%ymm8 			\n\t"\
		"vmulpd		 0x060(%%rdi),%%ymm3,%%ymm3			\n\t		vmulpd		 0x160(%%rdi),%%ymm9 ,%%ymm9 			\n\t"\
		"vsubpd		%%ymm2,%%ymm4,%%ymm4			\n\t		vsubpd			%%ymm8 ,%%ymm12,%%ymm12			\n\t"\
		"vaddpd		%%ymm3,%%ymm0,%%ymm0			\n\t		vaddpd			%%ymm9 ,%%ymm10,%%ymm10			\n\t"\
		"vmovaps		%%ymm4, 0x360(%%rbx)			\n\t		vmovaps		%%ymm12, 0x3e0(%%rbx)			\n\t"\
		"vmovaps		%%ymm0, 0x340(%%rbx)			\n\t		vmovaps		%%ymm10, 0x3c0(%%rbx)			\n\t"\
		"\n\t"\
		"/************************************************/	\n\t		/************************************************/	\n\t"\
		"/*...Block 1: t00,t10,t20,t30 -> r00,04,02,06:  */	\n\t		/*...Block 5: t08,t18,t28,t38 -> r08,0C,0A,0E:  */	\n\t"\
		"/************************************************/	\n\t		/************************************************/	\n\t"\
		"movq	%[__r00],%%rdx	/* base-addr in rcol = r08, so rdx+0x100 in rcol */	\n\t	vmovaps	(%%rsi),%%ymm10	/* isrt2 */	\n\t"\
		"vmovaps			  (%%rdx),%%ymm0				\n\t		vmovaps		 0x140(%%rdx),%%ymm12				\n\t"\
		"vmovaps		 0x020(%%rdx),%%ymm1				\n\t		vmovaps		 0x160(%%rdx),%%ymm13				\n\t"\
		"vmovaps		 0x080(%%rdx),%%ymm2				\n\t		vmovaps		 0x1c0(%%rdx),%%ymm8 				\n\t"\
		"vmovaps		 0x0a0(%%rdx),%%ymm3				\n\t		vmovaps		 0x1e0(%%rdx),%%ymm9 				\n\t"\
		"vsubpd		 0x080(%%rdx),%%ymm0,%%ymm0			\n\t		vaddpd		 0x160(%%rdx),%%ymm12,%%ymm12			\n\t"\
		"vsubpd		 0x0a0(%%rdx),%%ymm1,%%ymm1			\n\t		vsubpd		 0x140(%%rdx),%%ymm13,%%ymm13			\n\t"\
		"vaddpd			  (%%rdx),%%ymm2,%%ymm2			\n\t		vsubpd		 0x1e0(%%rdx),%%ymm8 ,%%ymm8 			\n\t"\
		"vaddpd		 0x020(%%rdx),%%ymm3,%%ymm3			\n\t		vaddpd		 0x1c0(%%rdx),%%ymm9 ,%%ymm9 			\n\t"\
		"vmovaps		 0x040(%%rdx),%%ymm4				\n\t		vmulpd			%%ymm10,%%ymm12,%%ymm12			\n\t"\
		"vmovaps		 0x060(%%rdx),%%ymm5				\n\t		vmulpd			%%ymm10,%%ymm13,%%ymm13			\n\t"\
		"vmovaps		 0x0c0(%%rdx),%%ymm6				\n\t		vmulpd			%%ymm10,%%ymm8 ,%%ymm8 			\n\t"\
		"vmovaps		 0x0e0(%%rdx),%%ymm7				\n\t		vmulpd			%%ymm10,%%ymm9 ,%%ymm9 			\n\t"\
		"vsubpd		 0x0c0(%%rdx),%%ymm4,%%ymm4			\n\t		vmovaps		%%ymm12,%%ymm14					\n\t"\
		"vsubpd		 0x0e0(%%rdx),%%ymm5,%%ymm5			\n\t		vmovaps		%%ymm13,%%ymm15					\n\t"\
		"vaddpd		 0x040(%%rdx),%%ymm6,%%ymm6			\n\t		vsubpd			%%ymm8 ,%%ymm12,%%ymm12			\n\t"\
		"vaddpd		 0x060(%%rdx),%%ymm7,%%ymm7			\n\t		vsubpd			%%ymm9 ,%%ymm13,%%ymm13			\n\t"\
		"/* base-twiddle in l/rcol = c00/c04, so rcx+0x100 in rcol */	vaddpd			%%ymm8 ,%%ymm14,%%ymm14			\n\t"\
		"movq		%[__c10],%%rcx						\n\t		vaddpd			%%ymm9 ,%%ymm15,%%ymm15			\n\t"\
		"vaddpd			%%ymm6,%%ymm2,%%ymm2			\n\t		vmovaps		 0x100(%%rdx),%%ymm8 				\n\t"\
		"vaddpd			%%ymm7,%%ymm3,%%ymm3			\n\t		vmovaps		 0x120(%%rdx),%%ymm9 				\n\t"\
		"vmovaps		%%ymm2,      (%%rdx)			\n\t		vmovaps		 0x180(%%rdx),%%ymm10				\n\t"\
		"vmovaps		%%ymm3, 0x020(%%rdx)			\n\t		vmovaps		 0x1a0(%%rdx),%%ymm11				\n\t"\
		"vaddpd			%%ymm6,%%ymm6,%%ymm6			\n\t		vsubpd		 0x1a0(%%rdx),%%ymm8 ,%%ymm8 			\n\t"\
		"vaddpd			%%ymm7,%%ymm7,%%ymm7			\n\t		vsubpd		 0x180(%%rdx),%%ymm9 ,%%ymm9 			\n\t"\
		"vsubpd			%%ymm6,%%ymm2,%%ymm2			\n\t		vaddpd		 0x100(%%rdx),%%ymm11,%%ymm11			\n\t"\
		"vsubpd			%%ymm7,%%ymm3,%%ymm3			\n\t		vaddpd		 0x120(%%rdx),%%ymm10,%%ymm10			\n\t"\
		"vmovaps		%%ymm2,%%ymm6					\n\t		vsubpd			%%ymm12,%%ymm11,%%ymm11			\n\t"\
		"vmovaps		%%ymm3,%%ymm7					\n\t		vsubpd			%%ymm13,%%ymm9 ,%%ymm9 			\n\t"\
		"vmulpd			  (%%rcx),%%ymm2,%%ymm2	/* c10 */	\n\t	vaddpd			%%ymm12,%%ymm12,%%ymm12			\n\t"\
		"vmulpd			  (%%rcx),%%ymm3,%%ymm3			\n\t		vaddpd			%%ymm13,%%ymm13,%%ymm13			\n\t"\
		"vmulpd		 0x020(%%rcx),%%ymm6,%%ymm6			\n\t		vaddpd			%%ymm11,%%ymm12,%%ymm12			\n\t"\
		"vmulpd		 0x020(%%rcx),%%ymm7,%%ymm7			\n\t		vaddpd			%%ymm9 ,%%ymm13,%%ymm13			\n\t"\
		"vsubpd			%%ymm6,%%ymm3,%%ymm3			\n\t		vmovaps		%%ymm11, 0x140(%%rdx)			\n\t"\
		"vaddpd			%%ymm7,%%ymm2,%%ymm2			\n\t		vmovaps		%%ymm9 , 0x160(%%rdx)			\n\t"\
		"subq $0x40,%%rcx	/* put c00 in rcx to ease bookkeeping*/\n\t	vmovaps		%%ymm12,%%ymm11					\n\t"\
		"/* add0,1 in rax,rbx; __r00 in rdx: */							vmovaps		%%ymm13,%%ymm9 					\n\t"\
		"/* For each complex output octet, complex pairs having */		vmulpd		 0x100(%%rcx),%%ymm12,%%ymm12	/* c04 */\n\t"\
		"/* reads from offsets 0x0..,0x1..,0x2..,0x3.. go into  */		vmulpd		 0x100(%%rcx),%%ymm13,%%ymm13			\n\t"\
		"/* local-mem pairs rXY + 00/10, 04/14, 02/12, 06/16.   */		vmulpd		 0x120(%%rcx),%%ymm11,%%ymm11			\n\t"\
		"/* For 1st octet we read from offsets [0x2..,0x0..],   */		vmulpd		 0x120(%%rcx),%%ymm9 ,%%ymm9 			\n\t"\
		"/* [0x1..,0x3], other 3 octets use order [0,2],[1,3].  */		vsubpd			%%ymm11,%%ymm13,%%ymm13			\n\t"\
		"vmovaps	0x220(%%rbx),%%ymm7						\n\t		vaddpd			%%ymm9 ,%%ymm12,%%ymm12			\n\t"\
		"vmovaps	0x200(%%rbx),%%ymm6						\n\t		vmovaps	0x0a0(%%rbx),%%ymm11						\n\t"\
		"vmovaps	%%ymm3,0x060(%%rdx)						\n\t		vmovaps	0x080(%%rbx),%%ymm9 						\n\t"\
		"vmovaps	%%ymm2,0x040(%%rdx)		\n\t/* r02,03 */			vmovaps	%%ymm13,0x120(%%rdx)						\n\t"\
		"vmovaps	%%ymm7,0x260(%%rdx)						\n\t		vmovaps	%%ymm12,0x100(%%rdx)			\n\t"/* r08,09 */\
		"vmovaps	%%ymm6,0x240(%%rdx)		\n\t/* r12,13 */			vmovaps	%%ymm11,0x320(%%rdx)						\n\t"\
		"vmovaps		 0x020(%%rdx),%%ymm3				\n\t		vmovaps	%%ymm9 ,0x300(%%rdx)			\n\t"/* r18,19 */\
		"vmovaps			  (%%rdx),%%ymm2				\n\t		vmovaps		0x140(%%rdx),%%ymm12				\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm7						\n\t		vmovaps		0x160(%%rdx),%%ymm13				\n\t"\
		"vmovaps	0x000(%%rbx),%%ymm6						\n\t		vmovaps		%%ymm12,%%ymm11					\n\t"\
		"vmovaps	%%ymm3,0x020(%%rdx)						\n\t		vmovaps		%%ymm13,%%ymm9 					\n\t"\
		"vmovaps	%%ymm2,     (%%rdx)		\n\t/* r00,01 */			vmulpd		 0x140(%%rcx),%%ymm12,%%ymm12	/* c14 */\n\t"\
		"vmovaps	%%ymm7,0x220(%%rdx)						\n\t		vmulpd		 0x140(%%rcx),%%ymm13,%%ymm13			\n\t"\
		"vmovaps	%%ymm6,0x200(%%rdx)		\n\t/* r10,11 */			vmulpd		 0x160(%%rcx),%%ymm11,%%ymm11			\n\t"\
		"vaddpd			%%ymm5,%%ymm0,%%ymm0			\n\t		vmulpd		 0x160(%%rcx),%%ymm9 ,%%ymm9 			\n\t"\
		"vsubpd			%%ymm4,%%ymm1,%%ymm1			\n\t		vsubpd			%%ymm11,%%ymm13,%%ymm13			\n\t"\
		"vmovaps		%%ymm0,%%ymm2					\n\t		vaddpd			%%ymm9 ,%%ymm12,%%ymm12			\n\t"\
		"vmovaps		%%ymm1,%%ymm3					\n\t		vmovaps	0x2a0(%%rbx),%%ymm11						\n\t"\
		"vaddpd			%%ymm5,%%ymm5,%%ymm5			\n\t		vmovaps	0x280(%%rbx),%%ymm9 						\n\t"\
		"vaddpd			%%ymm4,%%ymm4,%%ymm4			\n\t		vmovaps	%%ymm13,0x160(%%rdx)						\n\t"\
		"vmovaps		%%ymm0,%%ymm6					\n\t		vmovaps	%%ymm12,0x140(%%rdx)			\n\t"/* r0a,0b */\
		"vmovaps		%%ymm1,%%ymm7					\n\t		vmovaps	%%ymm11,0x360(%%rdx)						\n\t"\
		"vmulpd		 0x080(%%rcx),%%ymm2,%%ymm2	/* c08 */	\n\t	vmovaps	%%ymm9 ,0x340(%%rdx)			\n\t"/* r1a,1b */\
		"vmulpd		 0x080(%%rcx),%%ymm3,%%ymm3			\n\t		vsubpd			%%ymm15,%%ymm8 ,%%ymm8 			\n\t"\
		"vmulpd		 0x0a0(%%rcx),%%ymm6,%%ymm6			\n\t		vsubpd			%%ymm14,%%ymm10,%%ymm10			\n\t"\
		"vmulpd		 0x0a0(%%rcx),%%ymm7,%%ymm7			\n\t		vaddpd			%%ymm15,%%ymm15,%%ymm15			\n\t"\
		"vsubpd			%%ymm6,%%ymm3,%%ymm3			\n\t		vaddpd			%%ymm14,%%ymm14,%%ymm14			\n\t"\
		"vaddpd			%%ymm7,%%ymm2,%%ymm2			\n\t		vaddpd			%%ymm8 ,%%ymm15,%%ymm15			\n\t"\
		"vmovaps	0x120(%%rbx),%%ymm7						\n\t		vaddpd			%%ymm10,%%ymm14,%%ymm14			\n\t"\
		"vmovaps	0x100(%%rbx),%%ymm6						\n\t		vmovaps		%%ymm15,%%ymm12					\n\t"\
		"vmovaps	%%ymm3,0x0a0(%%rdx)						\n\t		vmovaps		%%ymm10,%%ymm13					\n\t"\
		"vmovaps	%%ymm2,0x080(%%rdx)		\n\t/* r04,05 */			vmulpd		 0x180(%%rcx),%%ymm15,%%ymm15	/* c0C */\n\t"\
		"vmovaps	%%ymm7,0x2a0(%%rdx)						\n\t		vmulpd		 0x180(%%rcx),%%ymm10,%%ymm10			\n\t"\
		"vmovaps	%%ymm6,0x280(%%rdx)		\n\t/* r14,15 */			vmulpd		 0x1a0(%%rcx),%%ymm12,%%ymm12			\n\t"\
		"vsubpd			%%ymm5,%%ymm0,%%ymm0			\n\t		vmulpd		 0x1a0(%%rcx),%%ymm13,%%ymm13			\n\t"\
		"vaddpd			%%ymm4,%%ymm1,%%ymm1			\n\t		vsubpd			%%ymm12,%%ymm10,%%ymm10			\n\t"\
		"vmovaps		%%ymm0,%%ymm6					\n\t		vaddpd			%%ymm13,%%ymm15,%%ymm15			\n\t"\
		"vmovaps		%%ymm1,%%ymm7					\n\t		vmovaps	0x1a0(%%rbx),%%ymm13						\n\t"\
		"vmulpd		 0x0c0(%%rcx),%%ymm0,%%ymm0	/* c18 */	\n\t	vmovaps	0x180(%%rbx),%%ymm12						\n\t"\
		"vmulpd		 0x0c0(%%rcx),%%ymm1,%%ymm1			\n\t		vmovaps	%%ymm10,0x1a0(%%rdx)						\n\t"\
		"vmulpd		 0x0e0(%%rcx),%%ymm6,%%ymm6			\n\t		vmovaps	%%ymm15,0x180(%%rdx)			\n\t"/* r0c,0d */\
		"vmulpd		 0x0e0(%%rcx),%%ymm7,%%ymm7			\n\t		vmovaps	%%ymm13,0x3a0(%%rdx)						\n\t"\
		"vsubpd			%%ymm6,%%ymm1,%%ymm1			\n\t		vmovaps	%%ymm12,0x380(%%rdx)			\n\t"/* r1c,1d */\
		"vaddpd			%%ymm7,%%ymm0,%%ymm0			\n\t		vmovaps		%%ymm8 ,%%ymm12					\n\t"\
		"vmovaps	0x320(%%rbx),%%ymm7						\n\t		vmovaps		%%ymm14,%%ymm13					\n\t"\
		"vmovaps	0x300(%%rbx),%%ymm6						\n\t		vmulpd		 0x1c0(%%rcx),%%ymm8 ,%%ymm8 	/* c1C */\n\t"\
		"vmovaps	%%ymm1,0x0e0(%%rdx)						\n\t		vmulpd		 0x1c0(%%rcx),%%ymm14,%%ymm14			\n\t"\
		"vmovaps	%%ymm0,0x0c0(%%rdx)		\n\t/* r06,07 */			vmulpd		 0x1e0(%%rcx),%%ymm12,%%ymm12			\n\t"\
		"vmovaps	%%ymm7,0x2e0(%%rdx)						\n\t		vmulpd		 0x1e0(%%rcx),%%ymm13,%%ymm13			\n\t"\
		"vmovaps	%%ymm6,0x2c0(%%rdx)		\n\t/* r16,17 */			vsubpd			%%ymm12,%%ymm14,%%ymm14			\n\t"\
		"																vaddpd			%%ymm13,%%ymm8 ,%%ymm8 			\n\t"\
		"																vmovaps	0x3a0(%%rbx),%%ymm13						\n\t"\
		"																vmovaps	0x380(%%rbx),%%ymm12						\n\t"\
		"																vmovaps	%%ymm14,0x1e0(%%rdx)						\n\t"\
		"																vmovaps	%%ymm8 ,0x1c0(%%rdx)			\n\t"/* r0e,0f */\
		"																vmovaps	%%ymm13,0x3e0(%%rdx)						\n\t"\
		"																vmovaps	%%ymm12,0x3c0(%%rdx)			\n\t"/* r1e,1f */\
		"\n\t"\
		"/************************************************/	\n\t		/************************************************/	\n\t"\
		"/*...Block 3: t04,t14,t24,t34 -> r20,24,22,26:  */	\n\t		/*...Block 7: t0C,t1C,t2C,t3C -> r28,2C,2A,2E:  */	\n\t"\
		"/************************************************/	\n\t		/************************************************/	\n\t"\
		"movq		%[__r20],%%rdx								/* base-addr in rcol = r28, so rdx offset +0x100 vs lcol */	\n\t"\
		"leaq	0x020(%%rsi),%%rcx	/* cc0; Note cc0/ss0 are shared between lcol/rcol, so no rcx-offset until get to twiddles*/\n\t"\
		"vmovaps		 0x040(%%rdx),%%ymm4				\n\t		vmovaps		 0x140(%%rdx),%%ymm12				\n\t"\
		"vmovaps		 0x0c0(%%rdx),%%ymm0				\n\t		vmovaps		 0x1c0(%%rdx),%%ymm8 				\n\t"\
		"vmovaps		 0x060(%%rdx),%%ymm5				\n\t		vmovaps		 0x160(%%rdx),%%ymm13				\n\t"\
		"vmovaps		 0x0e0(%%rdx),%%ymm1				\n\t		vmovaps		 0x1e0(%%rdx),%%ymm9 				\n\t"\
		"vmovaps		 0x040(%%rdx),%%ymm6				\n\t		vmovaps		 0x140(%%rdx),%%ymm14				\n\t"\
		"vmovaps		 0x0c0(%%rdx),%%ymm2				\n\t		vmovaps		 0x1c0(%%rdx),%%ymm10				\n\t"\
		"vmovaps		 0x060(%%rdx),%%ymm7				\n\t		vmovaps		 0x160(%%rdx),%%ymm15				\n\t"\
		"vmovaps		 0x0e0(%%rdx),%%ymm3				\n\t		vmovaps		 0x1e0(%%rdx),%%ymm11				\n\t"\
		"vmulpd			  (%%rcx),%%ymm4,%%ymm4			\n\t		vmulpd		 0x020(%%rcx),%%ymm12,%%ymm12			\n\t"\
		"vmulpd		 0x020(%%rcx),%%ymm0,%%ymm0			\n\t		vmulpd			  (%%rcx),%%ymm8 ,%%ymm8 			\n\t"\
		"vmulpd			  (%%rcx),%%ymm5,%%ymm5			\n\t		vmulpd		 0x020(%%rcx),%%ymm13,%%ymm13			\n\t"\
		"vmulpd		 0x020(%%rcx),%%ymm1,%%ymm1			\n\t		vmulpd			  (%%rcx),%%ymm9 ,%%ymm9 			\n\t"\
		"vmulpd		 0x020(%%rcx),%%ymm6,%%ymm6			\n\t		vmulpd			  (%%rcx),%%ymm14,%%ymm14			\n\t"\
		"vmulpd			  (%%rcx),%%ymm2,%%ymm2			\n\t		vmulpd		 0x020(%%rcx),%%ymm10,%%ymm10			\n\t"\
		"vmulpd		 0x020(%%rcx),%%ymm7,%%ymm7			\n\t		vmulpd			  (%%rcx),%%ymm15,%%ymm15			\n\t"\
		"vmulpd			  (%%rcx),%%ymm3,%%ymm3			\n\t		vmulpd		 0x020(%%rcx),%%ymm11,%%ymm11			\n\t"\
		"vsubpd			%%ymm6,%%ymm5,%%ymm5			\n\t		vsubpd			%%ymm14,%%ymm13,%%ymm13			\n\t"\
		"vsubpd			%%ymm2,%%ymm1,%%ymm1			\n\t		vsubpd			%%ymm10,%%ymm9 ,%%ymm9 			\n\t"\
		"vaddpd			%%ymm7,%%ymm4,%%ymm4			\n\t		vaddpd			%%ymm15,%%ymm12,%%ymm12			\n\t"\
		"vaddpd			%%ymm3,%%ymm0,%%ymm0			\n\t		vaddpd			%%ymm11,%%ymm8 ,%%ymm8 			\n\t"\
		"vmovaps		%%ymm5,%%ymm7					\n\t		vmovaps		%%ymm13,%%ymm15					\n\t"\
		"vmovaps		%%ymm4,%%ymm6					\n\t		vmovaps		%%ymm12,%%ymm14					\n\t"\
		"vaddpd			%%ymm0,%%ymm4,%%ymm4			\n\t		vaddpd			%%ymm8 ,%%ymm12,%%ymm12			\n\t"\
		"vaddpd			%%ymm1,%%ymm5,%%ymm5			\n\t		vaddpd			%%ymm9 ,%%ymm13,%%ymm13			\n\t"\
		"vsubpd			%%ymm0,%%ymm6,%%ymm6			\n\t		vsubpd			%%ymm8 ,%%ymm14,%%ymm14			\n\t"\
		"vsubpd			%%ymm1,%%ymm7,%%ymm7			\n\t		vsubpd			%%ymm9 ,%%ymm15,%%ymm15			\n\t"\
		"vmovaps		 0x080(%%rdx),%%ymm2				\n\t		vmovaps		 0x180(%%rdx),%%ymm10				\n\t"\
		"vmovaps		 0x0a0(%%rdx),%%ymm3				\n\t		vmovaps		 0x1a0(%%rdx),%%ymm11				\n\t"\
		"vmovaps			  (%%rdx),%%ymm0				\n\t		vmovaps		 0x100(%%rdx),%%ymm8 				\n\t"\
		"vmovaps		 0x020(%%rdx),%%ymm1				\n\t		vmovaps		 0x120(%%rdx),%%ymm9 				\n\t"\
		"vaddpd		 0x0a0(%%rdx),%%ymm2,%%ymm2			\n\t		vsubpd		 0x1a0(%%rdx),%%ymm10,%%ymm10			\n\t"\
		"vsubpd		 0x080(%%rdx),%%ymm3,%%ymm3			\n\t		vaddpd		 0x180(%%rdx),%%ymm11,%%ymm11			\n\t"\
		"vmulpd			  (%%rsi),%%ymm2,%%ymm2			\n\t		vmulpd			  (%%rsi),%%ymm10,%%ymm10			\n\t"\
		"vmulpd			  (%%rsi),%%ymm3,%%ymm3			\n\t		vmulpd			  (%%rsi),%%ymm11,%%ymm11			\n\t"\
		"vsubpd			%%ymm2,%%ymm0,%%ymm0			\n\t		vsubpd			%%ymm10,%%ymm8 ,%%ymm8 			\n\t"\
		"vsubpd			%%ymm3,%%ymm1,%%ymm1			\n\t		vsubpd			%%ymm11,%%ymm9 ,%%ymm9 			\n\t"\
		"vaddpd			%%ymm2,%%ymm2,%%ymm2			\n\t		vaddpd			%%ymm10,%%ymm10,%%ymm10			\n\t"\
		"vaddpd			%%ymm3,%%ymm3,%%ymm3			\n\t		vaddpd			%%ymm11,%%ymm11,%%ymm11			\n\t"\
		"vaddpd			%%ymm0,%%ymm2,%%ymm2			\n\t		vaddpd			%%ymm8 ,%%ymm10,%%ymm10			\n\t"\
		"vaddpd			%%ymm1,%%ymm3,%%ymm3			\n\t		vaddpd			%%ymm9 ,%%ymm11,%%ymm11			\n\t"\
		"movq		%[__c02],%%rcx				/* base-twiddle addr in rcol = c06, so rcx offset +0x100 vs lcol */	\n\t"\
		"vsubpd			%%ymm4,%%ymm2,%%ymm2			\n\t		vsubpd			%%ymm14,%%ymm8 ,%%ymm8 			\n\t"\
		"vsubpd			%%ymm5,%%ymm3,%%ymm3			\n\t		vsubpd			%%ymm15,%%ymm9 ,%%ymm9 			\n\t"\
		"vaddpd			%%ymm4,%%ymm4,%%ymm4			\n\t		vaddpd			%%ymm14,%%ymm14,%%ymm14			\n\t"\
		"vaddpd			%%ymm5,%%ymm5,%%ymm5			\n\t		vaddpd			%%ymm15,%%ymm15,%%ymm15			\n\t"\
		"vaddpd			%%ymm2,%%ymm4,%%ymm4			\n\t		vaddpd			%%ymm8 ,%%ymm14,%%ymm14			\n\t"\
		"vaddpd			%%ymm3,%%ymm5,%%ymm5			\n\t		vaddpd			%%ymm9 ,%%ymm15,%%ymm15			\n\t"\
		"vmovaps		%%ymm2, 0x040(%%rdx)			\n\t		vmovaps		%%ymm8 , 0x140(%%rdx)			\n\t"\
		"vmovaps		%%ymm3, 0x060(%%rdx)			\n\t		vmovaps		%%ymm9 , 0x160(%%rdx)			\n\t"\
		"vmovaps		%%ymm4,%%ymm2					\n\t		vmovaps		%%ymm14,%%ymm8 					\n\t"\
		"vmovaps		%%ymm5,%%ymm3					\n\t		vmovaps		%%ymm15,%%ymm9 					\n\t"\
		"vmulpd			  (%%rcx),%%ymm4,%%ymm4			\n\t		vmulpd		 0x100(%%rcx),%%ymm14,%%ymm14			\n\t"\
		"vmulpd			  (%%rcx),%%ymm5,%%ymm5			\n\t		vmulpd		 0x100(%%rcx),%%ymm15,%%ymm15			\n\t"\
		"vmulpd		 0x020(%%rcx),%%ymm2,%%ymm2			\n\t		vmulpd		 0x120(%%rcx),%%ymm8 ,%%ymm8 			\n\t"\
		"vmulpd		 0x020(%%rcx),%%ymm3,%%ymm3			\n\t		vmulpd		 0x120(%%rcx),%%ymm9 ,%%ymm9 			\n\t"\
		"vsubpd			%%ymm2,%%ymm5,%%ymm5			\n\t		vsubpd			%%ymm8 ,%%ymm15,%%ymm15			\n\t"\
		"vaddpd			%%ymm3,%%ymm4,%%ymm4			\n\t		vaddpd			%%ymm9 ,%%ymm14,%%ymm14			\n\t"\
		"vmovaps	0x060(%%rbx),%%ymm3						\n\t		vmovaps	0x0e0(%%rbx),%%ymm9 						\n\t"\
		"vmovaps	0x040(%%rbx),%%ymm2						\n\t		vmovaps	0x0c0(%%rbx),%%ymm8 						\n\t"\
		"vmovaps	%%ymm5,0x020(%%rdx)						\n\t		vmovaps	%%ymm15,0x120(%%rdx)						\n\t"\
		"vmovaps	%%ymm4,     (%%rdx)			\n\t/* r20,21 */		vmovaps	%%ymm14,0x100(%%rdx)			\n\t"/* r28,29 */\
		"vmovaps	%%ymm3,0x220(%%rdx)						\n\t		vmovaps	%%ymm9 ,0x320(%%rdx)						\n\t"\
		"vmovaps	%%ymm2,0x200(%%rdx)			\n\t/* r30,31 */		vmovaps	%%ymm8 ,0x300(%%rdx)			\n\t"/* r38,39 */\
		"movq		%[__c12],%%rcx					\n\t"		/* rcol uses c16 */\
		"vmovaps		 0x040(%%rdx),%%ymm4				\n\t		vmovaps		 0x140(%%rdx),%%ymm14				\n\t"\
		"vmovaps		 0x060(%%rdx),%%ymm5				\n\t		vmovaps		 0x160(%%rdx),%%ymm15				\n\t"\
		"vmovaps		%%ymm4,%%ymm2					\n\t		vmovaps		%%ymm14,%%ymm8 					\n\t"\
		"vmovaps		%%ymm5,%%ymm3					\n\t		vmovaps		%%ymm15,%%ymm9 					\n\t"\
		"vmulpd			  (%%rcx),%%ymm4,%%ymm4			\n\t		vmulpd		 0x100(%%rcx),%%ymm14,%%ymm14			\n\t"\
		"vmulpd			  (%%rcx),%%ymm5,%%ymm5			\n\t		vmulpd		 0x100(%%rcx),%%ymm15,%%ymm15			\n\t"\
		"vmulpd		 0x020(%%rcx),%%ymm2,%%ymm2			\n\t		vmulpd		 0x120(%%rcx),%%ymm8 ,%%ymm8 			\n\t"\
		"vmulpd		 0x020(%%rcx),%%ymm3,%%ymm3			\n\t		vmulpd		 0x120(%%rcx),%%ymm9 ,%%ymm9 			\n\t"\
		"vsubpd			%%ymm2,%%ymm5,%%ymm5			\n\t		vsubpd			%%ymm8 ,%%ymm15,%%ymm15			\n\t"\
		"vaddpd			%%ymm3,%%ymm4,%%ymm4			\n\t		vaddpd			%%ymm9 ,%%ymm14,%%ymm14			\n\t"\
		"vmovaps	0x260(%%rbx),%%ymm3						\n\t		vmovaps	0x2e0(%%rbx),%%ymm9 						\n\t"\
		"vmovaps	0x240(%%rbx),%%ymm2						\n\t		vmovaps	0x2c0(%%rbx),%%ymm8 						\n\t"\
		"vmovaps	%%ymm5,0x060(%%rdx)						\n\t		vmovaps	%%ymm15,0x160(%%rdx)						\n\t"\
		"vmovaps	%%ymm4,0x040(%%rdx)			\n\t/* r22,23 */		vmovaps	%%ymm14,0x140(%%rdx)			\n\t"/* r2a,2b */\
		"vmovaps	%%ymm3,0x260(%%rdx)						\n\t		vmovaps	%%ymm9 ,0x360(%%rdx)						\n\t"\
		"vmovaps	%%ymm2,0x240(%%rdx)			\n\t/* r32,33 */		vmovaps	%%ymm8 ,0x340(%%rdx)			\n\t"/* r3a,3b */\
		"movq		%[__c0A],%%rcx					\n\t"		/* rcol uses c0E */\
		"vsubpd			%%ymm7,%%ymm0,%%ymm0			\n\t		vsubpd			%%ymm13,%%ymm10,%%ymm10			\n\t"\
		"vsubpd			%%ymm6,%%ymm1,%%ymm1			\n\t		vsubpd			%%ymm12,%%ymm11,%%ymm11			\n\t"\
		"vaddpd			%%ymm7,%%ymm7,%%ymm7			\n\t		vaddpd			%%ymm13,%%ymm13,%%ymm13			\n\t"\
		"vaddpd			%%ymm6,%%ymm6,%%ymm6			\n\t		vaddpd			%%ymm12,%%ymm12,%%ymm12			\n\t"\
		"vaddpd			%%ymm0,%%ymm7,%%ymm7			\n\t		vaddpd			%%ymm10,%%ymm13,%%ymm13			\n\t"\
		"vaddpd			%%ymm1,%%ymm6,%%ymm6			\n\t		vaddpd			%%ymm11,%%ymm12,%%ymm12			\n\t"\
		"vmovaps		%%ymm7,%%ymm4					\n\t		vmovaps		%%ymm13,%%ymm8 					\n\t"\
		"vmovaps		%%ymm1,%%ymm5					\n\t		vmovaps		%%ymm11,%%ymm9 					\n\t"\
		"vmulpd			  (%%rcx),%%ymm7,%%ymm7			\n\t		vmulpd		 0x100(%%rcx),%%ymm13,%%ymm13			\n\t"\
		"vmulpd			  (%%rcx),%%ymm1,%%ymm1			\n\t		vmulpd		 0x100(%%rcx),%%ymm11,%%ymm11			\n\t"\
		"vmulpd		 0x020(%%rcx),%%ymm4,%%ymm4			\n\t		vmulpd		 0x120(%%rcx),%%ymm8 ,%%ymm8 			\n\t"\
		"vmulpd		 0x020(%%rcx),%%ymm5,%%ymm5			\n\t		vmulpd		 0x120(%%rcx),%%ymm9 ,%%ymm9 			\n\t"\
		"vsubpd			%%ymm4,%%ymm1,%%ymm1			\n\t		vsubpd			%%ymm8 ,%%ymm11,%%ymm11			\n\t"\
		"vaddpd			%%ymm5,%%ymm7,%%ymm7			\n\t		vaddpd			%%ymm9 ,%%ymm13,%%ymm13			\n\t"\
		"vmovaps	0x160(%%rbx),%%ymm5						\n\t		vmovaps	0x1e0(%%rbx),%%ymm9 						\n\t"\
		"vmovaps	0x140(%%rbx),%%ymm4						\n\t		vmovaps	0x1c0(%%rbx),%%ymm8 						\n\t"\
		"vmovaps	%%ymm1,0x0a0(%%rdx)						\n\t		vmovaps	%%ymm11,0x1a0(%%rdx)						\n\t"\
		"vmovaps	%%ymm7,0x080(%%rdx)			\n\t/* r24,25 */		vmovaps	%%ymm13,0x180(%%rdx)			\n\t"/* r2c,2d */\
		"vmovaps	%%ymm5,0x2a0(%%rdx)						\n\t		vmovaps	%%ymm9 ,0x3a0(%%rdx)						\n\t"\
		"vmovaps	%%ymm4,0x280(%%rdx)			\n\t/* r34,35 */		vmovaps	%%ymm8 ,0x380(%%rdx)			\n\t"/* r3c,3d */\
		"movq		%[__c1A],%%rcx					\n\t"		/* rcol uses c1E */\
		"vmovaps		%%ymm0,%%ymm4					\n\t		vmovaps		%%ymm10,%%ymm8 					\n\t"\
		"vmovaps		%%ymm6,%%ymm5					\n\t		vmovaps		%%ymm12,%%ymm9 					\n\t"\
		"vmulpd			  (%%rcx),%%ymm0,%%ymm0			\n\t		vmulpd		 0x100(%%rcx),%%ymm10,%%ymm10			\n\t"\
		"vmulpd			  (%%rcx),%%ymm6,%%ymm6			\n\t		vmulpd		 0x100(%%rcx),%%ymm12,%%ymm12			\n\t"\
		"vmulpd		 0x020(%%rcx),%%ymm4,%%ymm4			\n\t		vmulpd		 0x120(%%rcx),%%ymm8 ,%%ymm8 			\n\t"\
		"vmulpd		 0x020(%%rcx),%%ymm5,%%ymm5			\n\t		vmulpd		 0x120(%%rcx),%%ymm9 ,%%ymm9 			\n\t"\
		"vsubpd			%%ymm4,%%ymm6,%%ymm6			\n\t		vsubpd			%%ymm8 ,%%ymm12,%%ymm12			\n\t"\
		"vaddpd			%%ymm5,%%ymm0,%%ymm0			\n\t		vaddpd			%%ymm9 ,%%ymm10,%%ymm10			\n\t"\
		"vmovaps	0x360(%%rbx),%%ymm5						\n\t		vmovaps	0x3e0(%%rbx),%%ymm9 						\n\t"\
		"vmovaps	0x340(%%rbx),%%ymm4						\n\t		vmovaps	0x3c0(%%rbx),%%ymm8 						\n\t"\
		"vmovaps	%%ymm6,0x0e0(%%rdx)						\n\t		vmovaps	%%ymm12,0x1e0(%%rdx)						\n\t"\
		"vmovaps	%%ymm0,0x0c0(%%rdx)			\n\t/* r26,27 */		vmovaps	%%ymm10,0x1c0(%%rdx)			\n\t"/* r2e,2f */\
		"vmovaps	%%ymm5,0x2e0(%%rdx)						\n\t		vmovaps	%%ymm9 ,0x3e0(%%rdx)						\n\t"\
		"vmovaps	%%ymm4,0x2c0(%%rdx)			\n\t/* r36,37 */		vmovaps	%%ymm8 ,0x3c0(%%rdx)			\n\t"/* r3e,3f */\
	/*******************************************/\
	/**** Finish with 4-way 'un'terleaving: ****/\
	/*******************************************/\
		"movq	%[__r00] ,%%rsi\n\t"\
		"movq	%[__add0],%%rax\n\t"\
		"movq	%[__add1],%%rbx\n\t"\
		"movq	%[__add2],%%rcx\n\t"\
		"movq	%[__add3],%%rdx\n\t"\
	"/* a[j+p0]: Inputs from r00 +0/1, 8/9, 16/17, 24/25. Outputs into [add0,add0+0x100,add1,add1+0x100]+0x0: */	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps		 (%%rsi),%%ymm1						\n\t		vmovaps	0x020(%%rsi),%%ymm3								\n\t"\
		"vmovaps	0x200(%%rsi),%%ymm5						\n\t		vmovaps	0x220(%%rsi),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps	0x400(%%rsi),%%ymm6						\n\t		vmovaps	0x420(%%rsi),%%ymm14							\n\t"\
		"vmovaps	0x600(%%rsi),%%ymm5						\n\t		vmovaps	0x620(%%rsi),%%ymm13							\n\t"\
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
	"/* a[j+p2]: Inputs from r08 +0/1, 8/9, 16/17, 24/25. Outputs into [add0,add0+0x100,add1,add1+0x100]+0x0: */	\n\t"\
		"addq	$0x100,%%rsi	\n\t"\
		"addq	$0x40,%%rax	\n\t"\
		"addq	$0x40,%%rbx	\n\t"\
		"addq	$0x40,%%rcx	\n\t"\
		"addq	$0x40,%%rdx	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps		 (%%rsi),%%ymm1						\n\t		vmovaps	0x020(%%rsi),%%ymm3								\n\t"\
		"vmovaps	0x200(%%rsi),%%ymm5						\n\t		vmovaps	0x220(%%rsi),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps	0x400(%%rsi),%%ymm6						\n\t		vmovaps	0x420(%%rsi),%%ymm14							\n\t"\
		"vmovaps	0x600(%%rsi),%%ymm5						\n\t		vmovaps	0x620(%%rsi),%%ymm13							\n\t"\
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
	"/* a[j+p4]: Inputs from r04 +0/1, 8/9, 16/17, 24/25. Outputs into [add0,add0+0x100,add1,add1+0x100]+0x0: */	\n\t"\
		"subq	$0x80,%%rsi	\n\t"\
		"addq	$0x40,%%rax	\n\t"\
		"addq	$0x40,%%rbx	\n\t"\
		"addq	$0x40,%%rcx	\n\t"\
		"addq	$0x40,%%rdx	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps		 (%%rsi),%%ymm1						\n\t		vmovaps	0x020(%%rsi),%%ymm3								\n\t"\
		"vmovaps	0x200(%%rsi),%%ymm5						\n\t		vmovaps	0x220(%%rsi),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps	0x400(%%rsi),%%ymm6						\n\t		vmovaps	0x420(%%rsi),%%ymm14							\n\t"\
		"vmovaps	0x600(%%rsi),%%ymm5						\n\t		vmovaps	0x620(%%rsi),%%ymm13							\n\t"\
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
	"/* a[j+p6]: Inputs from r0c +0/1, 8/9, 16/17, 24/25. Outputs into [add0,add0+0x100,add1,add1+0x100]+0x0: */	\n\t"\
		"addq	$0x100,%%rsi	\n\t"\
		"addq	$0x40,%%rax	\n\t"\
		"addq	$0x40,%%rbx	\n\t"\
		"addq	$0x40,%%rcx	\n\t"\
		"addq	$0x40,%%rdx	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps		 (%%rsi),%%ymm1						\n\t		vmovaps	0x020(%%rsi),%%ymm3								\n\t"\
		"vmovaps	0x200(%%rsi),%%ymm5						\n\t		vmovaps	0x220(%%rsi),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps	0x400(%%rsi),%%ymm6						\n\t		vmovaps	0x420(%%rsi),%%ymm14							\n\t"\
		"vmovaps	0x600(%%rsi),%%ymm5						\n\t		vmovaps	0x620(%%rsi),%%ymm13							\n\t"\
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
	"/* a[j+p8]: Inputs from r02 +0/1, 8/9, 16/17, 24/25. Outputs into [add0,add0+0x100,add1,add1+0x100]+0x0: */	\n\t"\
		"subq	$0x140,%%rsi	\n\t"\
		"addq	$0x40,%%rax	\n\t"\
		"addq	$0x40,%%rbx	\n\t"\
		"addq	$0x40,%%rcx	\n\t"\
		"addq	$0x40,%%rdx	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps		 (%%rsi),%%ymm1						\n\t		vmovaps	0x020(%%rsi),%%ymm3								\n\t"\
		"vmovaps	0x200(%%rsi),%%ymm5						\n\t		vmovaps	0x220(%%rsi),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps	0x400(%%rsi),%%ymm6						\n\t		vmovaps	0x420(%%rsi),%%ymm14							\n\t"\
		"vmovaps	0x600(%%rsi),%%ymm5						\n\t		vmovaps	0x620(%%rsi),%%ymm13							\n\t"\
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
	"/* a[j+p10]: Inputs from r0a +0/1, 8/9, 16/17, 24/25. Outputs into [add0,add0+0x100,add1,add1+0x100]+0x0: */	\n\t"\
		"addq	$0x100,%%rsi	\n\t"\
		"addq	$0x40,%%rax	\n\t"\
		"addq	$0x40,%%rbx	\n\t"\
		"addq	$0x40,%%rcx	\n\t"\
		"addq	$0x40,%%rdx	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps		 (%%rsi),%%ymm1						\n\t		vmovaps	0x020(%%rsi),%%ymm3								\n\t"\
		"vmovaps	0x200(%%rsi),%%ymm5						\n\t		vmovaps	0x220(%%rsi),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps	0x400(%%rsi),%%ymm6						\n\t		vmovaps	0x420(%%rsi),%%ymm14							\n\t"\
		"vmovaps	0x600(%%rsi),%%ymm5						\n\t		vmovaps	0x620(%%rsi),%%ymm13							\n\t"\
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
	"/* a[j+p12]: Inputs from r06 +0/1, 8/9, 16/17, 24/25. Outputs into [add0,add0+0x100,add1,add1+0x100]+0x0: */	\n\t"\
		"subq	$0x80,%%rsi	\n\t"\
		"addq	$0x40,%%rax	\n\t"\
		"addq	$0x40,%%rbx	\n\t"\
		"addq	$0x40,%%rcx	\n\t"\
		"addq	$0x40,%%rdx	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps		 (%%rsi),%%ymm1						\n\t		vmovaps	0x020(%%rsi),%%ymm3								\n\t"\
		"vmovaps	0x200(%%rsi),%%ymm5						\n\t		vmovaps	0x220(%%rsi),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps	0x400(%%rsi),%%ymm6						\n\t		vmovaps	0x420(%%rsi),%%ymm14							\n\t"\
		"vmovaps	0x600(%%rsi),%%ymm5						\n\t		vmovaps	0x620(%%rsi),%%ymm13							\n\t"\
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
	"/* a[j+p14]: Inputs from r0e +0/1, 8/9, 16/17, 24/25. Outputs into [add0,add0+0x100,add1,add1+0x100]+0x0: */	\n\t"\
		"addq	$0x100,%%rsi	\n\t"\
		"addq	$0x40,%%rax	\n\t"\
		"addq	$0x40,%%rbx	\n\t"\
		"addq	$0x40,%%rcx	\n\t"\
		"addq	$0x40,%%rdx	\n\t"\
	"/* Real parts: */										\n\t	/* Imag parts: */\n\t"\
		"vmovaps		 (%%rsi),%%ymm1						\n\t		vmovaps	0x020(%%rsi),%%ymm3								\n\t"\
		"vmovaps	0x200(%%rsi),%%ymm5						\n\t		vmovaps	0x220(%%rsi),%%ymm13							\n\t"\
		"vshufpd	$15 ,%%ymm5 ,%%ymm1 ,%%ymm7				\n\t		vshufpd	$15 ,%%ymm13,%%ymm3 ,%%ymm15					\n\t"\
		"vshufpd	$0  ,%%ymm5 ,%%ymm1 ,%%ymm1				\n\t		vshufpd	$0  ,%%ymm13,%%ymm3 ,%%ymm3						\n\t"\
		"vmovaps	0x400(%%rsi),%%ymm6						\n\t		vmovaps	0x420(%%rsi),%%ymm14							\n\t"\
		"vmovaps	0x600(%%rsi),%%ymm5						\n\t		vmovaps	0x620(%%rsi),%%ymm13							\n\t"\
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
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__r00] "m" (Xr00)\
		 ,[__r08] "m" (Xr08)\
		 ,[__r10] "m" (Xr10)\
		 ,[__r20] "m" (Xr20)\
		 ,[__r28] "m" (Xr28)\
		 ,[__r30] "m" (Xr30)\
		 ,[__c01] "m" (Xc01)\
		 ,[__c02] "m" (Xc02)\
		 ,[__c04] "m" (Xc04)\
		 ,[__c06] "m" (Xc06)\
		 ,[__c08] "m" (Xc08)\
		 ,[__c0A] "m" (Xc0A)\
		 ,[__c0C] "m" (Xc0C)\
		 ,[__c0E] "m" (Xc0E)\
		 ,[__c10] "m" (Xc10)\
		 ,[__c12] "m" (Xc12)\
		 ,[__c14] "m" (Xc14)\
		 ,[__c16] "m" (Xc16)\
		 ,[__c18] "m" (Xc18)\
		 ,[__c1A] "m" (Xc1A)\
		 ,[__c1C] "m" (Xc1C)\
		 ,[__c1E] "m" (Xc1E)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"		/* Clobbered registers */\
	);\
	}

#elif defined(USE_SSE2)

	#define SSE2_RADIX32_WRAPPER_DIF(Xadd0,Xadd1,Xr00,Xr10,Xr20,Xr30,Xisrt2,Xcc0,Xc00,Xc01,Xc02,Xc03,Xc05,Xc07)\
	{\
	__asm__ volatile (\
		/************************************************************************/\
		/* Forward DIF radix-32 pass on the interleaved block1 and block2 data: */\
		/************************************************************************/\
		/*...Block 0:	*/\
		"movq		%[__add0],%%rax						\n\t"\
		"movq	%[__isrt2],%%rsi							\n\t		movq		%%rsi,%%rdi			\n\t"\
		"movq		%[__add1],%%rbx						\n\t		addq	$0x470,%%rdi	/* two */	\n\t"\
		/*****	SSE2_RADIX4_DIF_4WRAPPER(c00,c08,c10,c18,r00) ****	SSE2_RADIX4_DIF_4WRAPPER_2NDOFTWO(c04,c0C,c14,c1C,r08)	*****/\
		"movq		%[__r00],%%rcx						\n\t		/*addq		$0x080,%%rcx // __r08 */	\n\t"\
		"movq		%[__c00],%%rdx						\n\t		/*addq		$0x080,%%rdx // __c04 */	\n\t"\
		"movaps			 (%%rax),%%xmm6						\n\t		movaps		0x40 (%%rax),%%xmm14		\n\t"\
		"movaps			 (%%rax),%%xmm0						\n\t		movaps		0x40 (%%rax),%%xmm8			\n\t"\
		"unpckhpd		 (%%rbx),%%xmm6						\n\t		unpckhpd	0x40 (%%rbx),%%xmm14		\n\t"\
		"unpcklpd		 (%%rbx),%%xmm0						\n\t		unpcklpd	0x40 (%%rbx),%%xmm8			\n\t"\
		"movaps		%%xmm6,0x200(%%rcx)				\n\t		movaps		%%xmm14,0x280(%%rcx)		\n\t"\
		"movaps		0x010(%%rax),%%xmm7						\n\t		movaps		0x50 (%%rax),%%xmm15		\n\t"\
		"movaps		0x010(%%rax),%%xmm1						\n\t		movaps		0x50 (%%rax),%%xmm9			\n\t"\
		"unpckhpd	0x010(%%rbx),%%xmm7						\n\t		unpckhpd	0x50 (%%rbx),%%xmm15		\n\t"\
		"unpcklpd	0x010(%%rbx),%%xmm1						\n\t		unpcklpd	0x50 (%%rbx),%%xmm9			\n\t"\
		"movaps		%%xmm7,0x210(%%rcx)				\n\t		movaps		%%xmm15,0x290(%%rcx)		\n\t"\
		"movaps		%%xmm0,%%xmm2						\n\t		movaps		%%xmm8,%%xmm10				\n\t"\
		"movaps		%%xmm1,%%xmm3						\n\t		movaps		%%xmm9,%%xmm11				\n\t"\
		"mulpd			 (%%rdx),%%xmm0						\n\t		mulpd		0x080(%%rdx),%%xmm8			\n\t"\
		"mulpd			 (%%rdx),%%xmm1						\n\t		mulpd		0x080(%%rdx),%%xmm9			\n\t"\
		"mulpd		0x010(%%rdx),%%xmm2						\n\t		mulpd		0x090(%%rdx),%%xmm10		\n\t"\
		"mulpd		0x010(%%rdx),%%xmm3						\n\t		mulpd		0x090(%%rdx),%%xmm11		\n\t"\
		"addpd		%%xmm2,%%xmm1						\n\t		addpd		%%xmm10,%%xmm9				\n\t"\
		"subpd		%%xmm3,%%xmm0						\n\t		subpd		%%xmm11,%%xmm8				\n\t"\
		"movaps		%%xmm0,%%xmm2						\n\t		movaps		%%xmm8,%%xmm10				\n\t"\
		"movaps		%%xmm1,%%xmm3						\n\t		movaps		%%xmm9,%%xmm11				\n\t"\
		"addq		$0x020,%%rdx																			\n\t"\
		"movaps		0x100(%%rax),%%xmm6						\n\t		movaps		0x140(%%rax),%%xmm14		\n\t"\
		"movaps		0x100(%%rax),%%xmm4						\n\t		movaps		0x140(%%rax),%%xmm12		\n\t"\
		"unpckhpd	0x100(%%rbx),%%xmm6						\n\t		unpckhpd	0x140(%%rbx),%%xmm14		\n\t"\
		"unpcklpd	0x100(%%rbx),%%xmm4						\n\t		unpcklpd	0x140(%%rbx),%%xmm12		\n\t"\
		"movaps		%%xmm6,0x220(%%rcx)				\n\t		movaps		%%xmm14,0x2a0(%%rcx)		\n\t"\
		"movaps		0x110(%%rax),%%xmm7						\n\t		movaps		0x150(%%rax),%%xmm15		\n\t"\
		"movaps		0x110(%%rax),%%xmm5						\n\t		movaps		0x150(%%rax),%%xmm13		\n\t"\
		"unpckhpd	0x110(%%rbx),%%xmm7						\n\t		unpckhpd	0x150(%%rbx),%%xmm15		\n\t"\
		"unpcklpd	0x110(%%rbx),%%xmm5						\n\t		unpcklpd	0x150(%%rbx),%%xmm13		\n\t"\
		"movaps		%%xmm7,0x230(%%rcx)				\n\t		movaps		%%xmm15,0x2b0(%%rcx)		\n\t"\
		"movaps		%%xmm4,%%xmm6						\n\t		movaps		%%xmm12,%%xmm14				\n\t"\
		"movaps		%%xmm5,%%xmm7						\n\t		movaps		%%xmm13,%%xmm15				\n\t"\
		"mulpd			 (%%rdx),%%xmm4						\n\t		mulpd		0x080(%%rdx),%%xmm12		\n\t"\
		"mulpd			 (%%rdx),%%xmm5						\n\t		mulpd		0x080(%%rdx),%%xmm13		\n\t"\
		"mulpd		0x010(%%rdx),%%xmm6						\n\t		mulpd		0x090(%%rdx),%%xmm14		\n\t"\
		"mulpd		0x010(%%rdx),%%xmm7						\n\t		mulpd		0x090(%%rdx),%%xmm15		\n\t"\
		"addpd		%%xmm6,%%xmm5						\n\t		addpd		%%xmm14,%%xmm13				\n\t"\
		"subpd		%%xmm7,%%xmm4						\n\t		subpd		%%xmm15,%%xmm12				\n\t"\
		"addpd		%%xmm4,%%xmm0						\n\t		addpd		%%xmm12,%%xmm8				\n\t"\
		"addpd		%%xmm5,%%xmm1						\n\t		addpd		%%xmm13,%%xmm9				\n\t"\
		"subpd		%%xmm4,%%xmm2						\n\t		subpd		%%xmm12,%%xmm10				\n\t"\
		"subpd		%%xmm5,%%xmm3						\n\t		subpd		%%xmm13,%%xmm11				\n\t"\
		"addq		$0x040,%%rdx																			\n\t"\
		"movaps		0x180(%%rax),%%xmm6						\n\t		movaps		0x1c0(%%rax),%%xmm14		\n\t"\
		"movaps		0x180(%%rax),%%xmm4						\n\t		movaps		0x1c0(%%rax),%%xmm12		\n\t"\
		"unpckhpd	0x180(%%rbx),%%xmm6						\n\t		unpckhpd	0x1c0(%%rbx),%%xmm14		\n\t"\
		"unpcklpd	0x180(%%rbx),%%xmm4						\n\t		unpcklpd	0x1c0(%%rbx),%%xmm12		\n\t"\
		"movaps		%%xmm6,0x260(%%rcx)				\n\t		movaps		%%xmm14,0x2e0(%%rcx)		\n\t"\
		"movaps		0x190(%%rax),%%xmm7						\n\t		movaps		0x1d0(%%rax),%%xmm15		\n\t"\
		"movaps		0x190(%%rax),%%xmm5						\n\t		movaps		0x1d0(%%rax),%%xmm13		\n\t"\
		"unpckhpd	0x190(%%rbx),%%xmm7						\n\t		unpckhpd	0x1d0(%%rbx),%%xmm15		\n\t"\
		"unpcklpd	0x190(%%rbx),%%xmm5						\n\t		unpcklpd	0x1d0(%%rbx),%%xmm13		\n\t"\
		"movaps		%%xmm7,0x270(%%rcx)				\n\t		movaps		%%xmm15,0x2f0(%%rcx)		\n\t"\
		"movaps		%%xmm4,%%xmm6						\n\t		movaps		%%xmm12,%%xmm14				\n\t"\
		"movaps		%%xmm5,%%xmm7						\n\t		movaps		%%xmm13,%%xmm15				\n\t"\
		"mulpd			 (%%rdx),%%xmm4						\n\t		mulpd		0x080(%%rdx),%%xmm12		\n\t"\
		"mulpd			 (%%rdx),%%xmm5						\n\t		mulpd		0x080(%%rdx),%%xmm13		\n\t"\
		"mulpd		0x010(%%rdx),%%xmm6						\n\t		mulpd		0x090(%%rdx),%%xmm14		\n\t"\
		"mulpd		0x010(%%rdx),%%xmm7						\n\t		mulpd		0x090(%%rdx),%%xmm15		\n\t"\
		"addpd		%%xmm6,%%xmm5						\n\t		addpd		%%xmm14,%%xmm13				\n\t"\
		"subpd		%%xmm7,%%xmm4						\n\t		subpd		%%xmm15,%%xmm12				\n\t"\
		"movaps		%%xmm5,0x010(%%rcx)				\n\t		movaps		%%xmm13,0x090(%%rcx)		\n\t"\
		"movaps		%%xmm4,     (%%rcx)				\n\t		movaps		%%xmm12,0x080(%%rcx)		\n\t"\
		"subq		$0x020,%%rdx																			\n\t"\
		"movaps		0x080(%%rax),%%xmm6						\n\t		movaps		0x0c0(%%rax),%%xmm14		\n\t"\
		"movaps		0x080(%%rax),%%xmm4						\n\t		movaps		0x0c0(%%rax),%%xmm12		\n\t"\
		"unpckhpd	0x080(%%rbx),%%xmm6						\n\t		unpckhpd	0x0c0(%%rbx),%%xmm14		\n\t"\
		"unpcklpd	0x080(%%rbx),%%xmm4						\n\t		unpcklpd	0x0c0(%%rbx),%%xmm12		\n\t"\
		"movaps		%%xmm6,0x240(%%rcx)				\n\t		movaps		%%xmm14,0x2c0(%%rcx)		\n\t"\
		"movaps		0x090(%%rax),%%xmm7						\n\t		movaps		0x0d0(%%rax),%%xmm15		\n\t"\
		"movaps		0x090(%%rax),%%xmm5						\n\t		movaps		0x0d0(%%rax),%%xmm13		\n\t"\
		"unpckhpd	0x090(%%rbx),%%xmm7						\n\t		unpckhpd	0x0d0(%%rbx),%%xmm15		\n\t"\
		"unpcklpd	0x090(%%rbx),%%xmm5						\n\t		unpcklpd	0x0d0(%%rbx),%%xmm13		\n\t"\
		"movaps		%%xmm7,0x250(%%rcx)				\n\t		movaps		%%xmm15,0x2d0(%%rcx)		\n\t"\
		"movaps		%%xmm4,%%xmm6						\n\t		movaps		%%xmm12,%%xmm14				\n\t"\
		"movaps		%%xmm5,%%xmm7						\n\t		movaps		%%xmm13,%%xmm15				\n\t"\
		"mulpd			 (%%rdx),%%xmm4						\n\t		mulpd		0x080(%%rdx),%%xmm12		\n\t"\
		"mulpd			 (%%rdx),%%xmm5						\n\t		mulpd		0x080(%%rdx),%%xmm13		\n\t"\
		"mulpd		0x010(%%rdx),%%xmm6						\n\t		mulpd		0x090(%%rdx),%%xmm14		\n\t"\
		"mulpd		0x010(%%rdx),%%xmm7						\n\t		mulpd		0x090(%%rdx),%%xmm15		\n\t"\
		"addpd		%%xmm6,%%xmm5						\n\t		addpd		%%xmm14,%%xmm13				\n\t"\
		"subpd		%%xmm7,%%xmm4						\n\t		subpd		%%xmm15,%%xmm12				\n\t"\
		"movaps		%%xmm5,%%xmm7						\n\t		movaps		%%xmm13,%%xmm15				\n\t"\
		"movaps		%%xmm4,%%xmm6						\n\t		movaps		%%xmm12,%%xmm14				\n\t"\
		"subpd			 (%%rcx),%%xmm4						\n\t		subpd		0x080(%%rcx),%%xmm12		\n\t"\
		"subpd		0x010(%%rcx),%%xmm5						\n\t		subpd		0x090(%%rcx),%%xmm13		\n\t"\
		"addpd			 (%%rcx),%%xmm6						\n\t		addpd		0x080(%%rcx),%%xmm14		\n\t"\
		"addpd		0x010(%%rcx),%%xmm7						\n\t		addpd		0x090(%%rcx),%%xmm15		\n\t"\
		"subpd		%%xmm6,%%xmm0						\n\t		subpd		%%xmm14,%%xmm8				\n\t"\
		"subpd		%%xmm7,%%xmm1						\n\t		subpd		%%xmm15,%%xmm9				\n\t"\
		"/*movaps	%%xmm0,0x040(%%rcx)	*/			\n\t	/*	movaps		%%xmm8,0x0c0(%%rcx)		*/	\n\t"\
		"/*movaps	%%xmm1,0x050(%%rcx)	*/			\n\t	/*	movaps		%%xmm9,0x0d0(%%rcx)		*/	\n\t"\
		"mulpd		(%%rdi),%%xmm6						\n\t		mulpd		(%%rdi),%%xmm14				\n\t"\
		"mulpd		(%%rdi),%%xmm7						\n\t		mulpd		(%%rdi),%%xmm15				\n\t"\
		"addpd		%%xmm0,%%xmm6						\n\t		addpd		%%xmm8,%%xmm14				\n\t"\
		"addpd		%%xmm1,%%xmm7						\n\t		addpd		%%xmm9,%%xmm15				\n\t"\
		"/*movaps	%%xmm6,     (%%rcx)	*/			\n\t		movaps		%%xmm14,0x080(%%rcx)		\n\t"\
		"/*movaps	%%xmm7,0x010(%%rcx)	*/			\n\t		movaps		%%xmm15,0x090(%%rcx)		\n\t"\
		"subpd		%%xmm5,%%xmm2						\n\t		subpd		%%xmm13,%%xmm10				\n\t"\
		"subpd		%%xmm4,%%xmm3						\n\t		subpd		%%xmm12,%%xmm11				\n\t"\
		"/*movaps	%%xmm2,0x020(%%rcx)	*/			\n\t		mulpd		(%%rdi),%%xmm13				\n\t"\
		"/*movaps	%%xmm3,0x070(%%rcx)	*/			\n\t		mulpd		(%%rdi),%%xmm12				\n\t"\
		"mulpd		(%%rdi),%%xmm5						\n\t		addpd		%%xmm10,%%xmm13				\n\t"\
		"mulpd		(%%rdi),%%xmm4						\n\t		addpd		%%xmm11,%%xmm12				\n\t"\
		"addpd		%%xmm2,%%xmm5						\n\t		movaps		%%xmm10,%%xmm14				\n\t"\
		"addpd		%%xmm3,%%xmm4						\n\t		movaps		%%xmm13,%%xmm15				\n\t"\
		"/*movaps	%%xmm5,0x060(%%rcx)	*/			\n\t		subpd		%%xmm12,%%xmm10				\n\t"\
		"/*movaps	%%xmm4,0x030(%%rcx)	*/						subpd		%%xmm11,%%xmm13				\n\t"\
		"																addpd		%%xmm12,%%xmm14				\n\t"\
		"																addpd		%%xmm11,%%xmm15				\n\t"\
		"																mulpd		(%%rsi),%%xmm10				\n\t"\
		"																mulpd		(%%rsi),%%xmm13				\n\t"\
		"																mulpd		(%%rsi),%%xmm14				\n\t"\
		"																mulpd		(%%rsi),%%xmm15				\n\t"\
		"															/*	movaps		%%xmm10,0x0a0(%%rcx)	*/	\n\t"\
		"															/*	movaps		%%xmm13,0x0e0(%%rcx)	*/	\n\t"\
		"															/*	movaps		%%xmm14,0x0b0(%%rcx)	*/	\n\t"\
		"															/*	movaps		%%xmm15,0x0f0(%%rcx)	*/	\n\t"\
		/*****	SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r00)	*****/\
		"																movaps		0x080(%%rcx),%%xmm11		\n\t"\
		"																movaps		0x090(%%rcx),%%xmm12		\n\t"\
		"subpd		%%xmm11,%%xmm6						\n\t		subpd		%%xmm10,%%xmm2			\n\t"\
		"subpd		%%xmm9,%%xmm0						\n\t		subpd		%%xmm15,%%xmm5			\n\t"\
		"subpd		%%xmm12,%%xmm7						\n\t		subpd		%%xmm14,%%xmm4			\n\t"\
		"subpd		%%xmm8,%%xmm1						\n\t		subpd		%%xmm13,%%xmm3			\n\t"\
		"mulpd		(%%rdi),%%xmm11					\n\t		mulpd		(%%rdi),%%xmm10		\n\t"\
		"mulpd		(%%rdi),%%xmm9						\n\t		mulpd		(%%rdi),%%xmm15		\n\t"\
		"mulpd		(%%rdi),%%xmm12					\n\t		mulpd		(%%rdi),%%xmm14		\n\t"\
		"mulpd		(%%rdi),%%xmm8						\n\t		mulpd		(%%rdi),%%xmm13		\n\t"\
		"addpd		%%xmm6,%%xmm11					\n\t		addpd		%%xmm2,%%xmm10		\n\t"\
		"addpd		%%xmm0,%%xmm9						\n\t		addpd		%%xmm5,%%xmm15		\n\t"\
		"addpd		%%xmm7,%%xmm12					\n\t		addpd		%%xmm4,%%xmm14		\n\t"\
		"addpd		%%xmm1,%%xmm8						\n\t		addpd		%%xmm3,%%xmm13		\n\t"\
		"movaps		%%xmm6,0x080(%%rcx)				\n\t		movaps		%%xmm2,0x0a0(%%rcx)	\n\t"\
		"movaps		%%xmm0,0x040(%%rcx)				\n\t		movaps		%%xmm5,0x060(%%rcx)	\n\t"\
		"movaps		%%xmm7,0x090(%%rcx)				\n\t		movaps		%%xmm4,0x0b0(%%rcx)	\n\t"\
		"movaps		%%xmm1,0x0d0(%%rcx)				\n\t		movaps		%%xmm3,0x0f0(%%rcx)	\n\t"\
		"movaps		%%xmm11,     (%%rcx)				\n\t		movaps		%%xmm10,0x020(%%rcx)	\n\t"\
		"movaps		%%xmm9,0x0c0(%%rcx)				\n\t		movaps		%%xmm15,0x0e0(%%rcx)	\n\t"\
		"movaps		%%xmm12,0x010(%%rcx)				\n\t		movaps		%%xmm14,0x030(%%rcx)	\n\t"\
		"movaps		%%xmm8,0x050(%%rcx)				\n\t		movaps		%%xmm13,0x070(%%rcx)	\n\t"\
		"\n\t"\
		/*...Block 2:	*/\
		"addq		$0x20,%%rax\n\t"\
		"addq		$0x20,%%rbx\n\t"\
		/*****	SSE2_RADIX4_DIF_4WRAPPER(c02,c0A,c12,c1A,r10) ****	SSE2_RADIX4_DIF_4WRAPPER_2NDOFTWO(c06,c0E,c16,c1E,r18)	*****/\
		"movq		%[__r10],%%rcx						\n\t		/*addq		$0x080,%%rcx // __r18 */	\n\t"\
		"movq		%[__c02],%%rdx						\n\t		/*addq		$0x080,%%rdx // __c06 */	\n\t"\
		"movaps			 (%%rax),%%xmm6						\n\t		movaps		0x40 (%%rax),%%xmm14		\n\t"\
		"movaps			 (%%rax),%%xmm0						\n\t		movaps		0x40 (%%rax),%%xmm8			\n\t"\
		"unpckhpd		 (%%rbx),%%xmm6						\n\t		unpckhpd	0x40 (%%rbx),%%xmm14		\n\t"\
		"unpcklpd		 (%%rbx),%%xmm0						\n\t		unpcklpd	0x40 (%%rbx),%%xmm8			\n\t"\
		"movaps		%%xmm6,0x200(%%rcx)				\n\t		movaps		%%xmm14,0x280(%%rcx)		\n\t"\
		"movaps		0x010(%%rax),%%xmm7						\n\t		movaps		0x50 (%%rax),%%xmm15		\n\t"\
		"movaps		0x010(%%rax),%%xmm1						\n\t		movaps		0x50 (%%rax),%%xmm9			\n\t"\
		"unpckhpd	0x010(%%rbx),%%xmm7						\n\t		unpckhpd	0x50 (%%rbx),%%xmm15		\n\t"\
		"unpcklpd	0x010(%%rbx),%%xmm1						\n\t		unpcklpd	0x50 (%%rbx),%%xmm9			\n\t"\
		"movaps		%%xmm7,0x210(%%rcx)				\n\t		movaps		%%xmm15,0x290(%%rcx)		\n\t"\
		"movaps		%%xmm0,%%xmm2						\n\t		movaps		%%xmm8,%%xmm10				\n\t"\
		"movaps		%%xmm1,%%xmm3						\n\t		movaps		%%xmm9,%%xmm11				\n\t"\
		"mulpd			 (%%rdx),%%xmm0						\n\t		mulpd		0x080(%%rdx),%%xmm8			\n\t"\
		"mulpd			 (%%rdx),%%xmm1						\n\t		mulpd		0x080(%%rdx),%%xmm9			\n\t"\
		"mulpd		0x010(%%rdx),%%xmm2						\n\t		mulpd		0x090(%%rdx),%%xmm10		\n\t"\
		"mulpd		0x010(%%rdx),%%xmm3						\n\t		mulpd		0x090(%%rdx),%%xmm11		\n\t"\
		"addpd		%%xmm2,%%xmm1						\n\t		addpd		%%xmm10,%%xmm9				\n\t"\
		"subpd		%%xmm3,%%xmm0						\n\t		subpd		%%xmm11,%%xmm8				\n\t"\
		"movaps		%%xmm0,%%xmm2						\n\t		movaps		%%xmm8,%%xmm10				\n\t"\
		"movaps		%%xmm1,%%xmm3						\n\t		movaps		%%xmm9,%%xmm11				\n\t"\
		"addq		$0x020,%%rdx																			\n\t"\
		"movaps		0x100(%%rax),%%xmm6						\n\t		movaps		0x140(%%rax),%%xmm14		\n\t"\
		"movaps		0x100(%%rax),%%xmm4						\n\t		movaps		0x140(%%rax),%%xmm12		\n\t"\
		"unpckhpd	0x100(%%rbx),%%xmm6						\n\t		unpckhpd	0x140(%%rbx),%%xmm14		\n\t"\
		"unpcklpd	0x100(%%rbx),%%xmm4						\n\t		unpcklpd	0x140(%%rbx),%%xmm12		\n\t"\
		"movaps		%%xmm6,0x220(%%rcx)				\n\t		movaps		%%xmm14,0x2a0(%%rcx)		\n\t"\
		"movaps		0x110(%%rax),%%xmm7						\n\t		movaps		0x150(%%rax),%%xmm15		\n\t"\
		"movaps		0x110(%%rax),%%xmm5						\n\t		movaps		0x150(%%rax),%%xmm13		\n\t"\
		"unpckhpd	0x110(%%rbx),%%xmm7						\n\t		unpckhpd	0x150(%%rbx),%%xmm15		\n\t"\
		"unpcklpd	0x110(%%rbx),%%xmm5						\n\t		unpcklpd	0x150(%%rbx),%%xmm13		\n\t"\
		"movaps		%%xmm7,0x230(%%rcx)				\n\t		movaps		%%xmm15,0x2b0(%%rcx)		\n\t"\
		"movaps		%%xmm4,%%xmm6						\n\t		movaps		%%xmm12,%%xmm14				\n\t"\
		"movaps		%%xmm5,%%xmm7						\n\t		movaps		%%xmm13,%%xmm15				\n\t"\
		"mulpd			 (%%rdx),%%xmm4						\n\t		mulpd		0x080(%%rdx),%%xmm12		\n\t"\
		"mulpd			 (%%rdx),%%xmm5						\n\t		mulpd		0x080(%%rdx),%%xmm13		\n\t"\
		"mulpd		0x010(%%rdx),%%xmm6						\n\t		mulpd		0x090(%%rdx),%%xmm14		\n\t"\
		"mulpd		0x010(%%rdx),%%xmm7						\n\t		mulpd		0x090(%%rdx),%%xmm15		\n\t"\
		"addpd		%%xmm6,%%xmm5						\n\t		addpd		%%xmm14,%%xmm13				\n\t"\
		"subpd		%%xmm7,%%xmm4						\n\t		subpd		%%xmm15,%%xmm12				\n\t"\
		"addpd		%%xmm4,%%xmm0						\n\t		addpd		%%xmm12,%%xmm8				\n\t"\
		"addpd		%%xmm5,%%xmm1						\n\t		addpd		%%xmm13,%%xmm9				\n\t"\
		"subpd		%%xmm4,%%xmm2						\n\t		subpd		%%xmm12,%%xmm10				\n\t"\
		"subpd		%%xmm5,%%xmm3						\n\t		subpd		%%xmm13,%%xmm11				\n\t"\
		"addq		$0x040,%%rdx																			\n\t"\
		"movaps		0x180(%%rax),%%xmm6						\n\t		movaps		0x1c0(%%rax),%%xmm14		\n\t"\
		"movaps		0x180(%%rax),%%xmm4						\n\t		movaps		0x1c0(%%rax),%%xmm12		\n\t"\
		"unpckhpd	0x180(%%rbx),%%xmm6						\n\t		unpckhpd	0x1c0(%%rbx),%%xmm14		\n\t"\
		"unpcklpd	0x180(%%rbx),%%xmm4						\n\t		unpcklpd	0x1c0(%%rbx),%%xmm12		\n\t"\
		"movaps		%%xmm6,0x260(%%rcx)				\n\t		movaps		%%xmm14,0x2e0(%%rcx)		\n\t"\
		"movaps		0x190(%%rax),%%xmm7						\n\t		movaps		0x1d0(%%rax),%%xmm15		\n\t"\
		"movaps		0x190(%%rax),%%xmm5						\n\t		movaps		0x1d0(%%rax),%%xmm13		\n\t"\
		"unpckhpd	0x190(%%rbx),%%xmm7						\n\t		unpckhpd	0x1d0(%%rbx),%%xmm15		\n\t"\
		"unpcklpd	0x190(%%rbx),%%xmm5						\n\t		unpcklpd	0x1d0(%%rbx),%%xmm13		\n\t"\
		"movaps		%%xmm7,0x270(%%rcx)				\n\t		movaps		%%xmm15,0x2f0(%%rcx)		\n\t"\
		"movaps		%%xmm4,%%xmm6						\n\t		movaps		%%xmm12,%%xmm14				\n\t"\
		"movaps		%%xmm5,%%xmm7						\n\t		movaps		%%xmm13,%%xmm15				\n\t"\
		"mulpd			 (%%rdx),%%xmm4						\n\t		mulpd		0x080(%%rdx),%%xmm12		\n\t"\
		"mulpd			 (%%rdx),%%xmm5						\n\t		mulpd		0x080(%%rdx),%%xmm13		\n\t"\
		"mulpd		0x010(%%rdx),%%xmm6						\n\t		mulpd		0x090(%%rdx),%%xmm14		\n\t"\
		"mulpd		0x010(%%rdx),%%xmm7						\n\t		mulpd		0x090(%%rdx),%%xmm15		\n\t"\
		"addpd		%%xmm6,%%xmm5						\n\t		addpd		%%xmm14,%%xmm13				\n\t"\
		"subpd		%%xmm7,%%xmm4						\n\t		subpd		%%xmm15,%%xmm12				\n\t"\
		"movaps		%%xmm5,0x010(%%rcx)				\n\t		movaps		%%xmm13,0x090(%%rcx)		\n\t"\
		"movaps		%%xmm4,     (%%rcx)				\n\t		movaps		%%xmm12,0x080(%%rcx)		\n\t"\
		"subq		$0x020,%%rdx																			\n\t"\
		"movaps		0x080(%%rax),%%xmm6						\n\t		movaps		0x0c0(%%rax),%%xmm14		\n\t"\
		"movaps		0x080(%%rax),%%xmm4						\n\t		movaps		0x0c0(%%rax),%%xmm12		\n\t"\
		"unpckhpd	0x080(%%rbx),%%xmm6						\n\t		unpckhpd	0x0c0(%%rbx),%%xmm14		\n\t"\
		"unpcklpd	0x080(%%rbx),%%xmm4						\n\t		unpcklpd	0x0c0(%%rbx),%%xmm12		\n\t"\
		"movaps		%%xmm6,0x240(%%rcx)				\n\t		movaps		%%xmm14,0x2c0(%%rcx)		\n\t"\
		"movaps		0x090(%%rax),%%xmm7						\n\t		movaps		0x0d0(%%rax),%%xmm15		\n\t"\
		"movaps		0x090(%%rax),%%xmm5						\n\t		movaps		0x0d0(%%rax),%%xmm13		\n\t"\
		"unpckhpd	0x090(%%rbx),%%xmm7						\n\t		unpckhpd	0x0d0(%%rbx),%%xmm15		\n\t"\
		"unpcklpd	0x090(%%rbx),%%xmm5						\n\t		unpcklpd	0x0d0(%%rbx),%%xmm13		\n\t"\
		"movaps		%%xmm7,0x250(%%rcx)				\n\t		movaps		%%xmm15,0x2d0(%%rcx)		\n\t"\
		"movaps		%%xmm4,%%xmm6						\n\t		movaps		%%xmm12,%%xmm14				\n\t"\
		"movaps		%%xmm5,%%xmm7						\n\t		movaps		%%xmm13,%%xmm15				\n\t"\
		"mulpd			 (%%rdx),%%xmm4						\n\t		mulpd		0x080(%%rdx),%%xmm12		\n\t"\
		"mulpd			 (%%rdx),%%xmm5						\n\t		mulpd		0x080(%%rdx),%%xmm13		\n\t"\
		"mulpd		0x010(%%rdx),%%xmm6						\n\t		mulpd		0x090(%%rdx),%%xmm14		\n\t"\
		"mulpd		0x010(%%rdx),%%xmm7						\n\t		mulpd		0x090(%%rdx),%%xmm15		\n\t"\
		"addpd		%%xmm6,%%xmm5						\n\t		addpd		%%xmm14,%%xmm13				\n\t"\
		"subpd		%%xmm7,%%xmm4						\n\t		subpd		%%xmm15,%%xmm12				\n\t"\
		"movaps		%%xmm5,%%xmm7						\n\t		movaps		%%xmm13,%%xmm15				\n\t"\
		"movaps		%%xmm4,%%xmm6						\n\t		movaps		%%xmm12,%%xmm14				\n\t"\
		"subpd			 (%%rcx),%%xmm4						\n\t		subpd		0x080(%%rcx),%%xmm12		\n\t"\
		"subpd		0x010(%%rcx),%%xmm5						\n\t		subpd		0x090(%%rcx),%%xmm13		\n\t"\
		"addpd			 (%%rcx),%%xmm6						\n\t		addpd		0x080(%%rcx),%%xmm14		\n\t"\
		"addpd		0x010(%%rcx),%%xmm7						\n\t		addpd		0x090(%%rcx),%%xmm15		\n\t"\
		"subpd		%%xmm6,%%xmm0						\n\t		subpd		%%xmm14,%%xmm8				\n\t"\
		"subpd		%%xmm7,%%xmm1						\n\t		subpd		%%xmm15,%%xmm9				\n\t"\
		"/*movaps	%%xmm0,0x040(%%rcx)	*/			\n\t	/*	movaps		%%xmm8,0x0c0(%%rcx)		*/	\n\t"\
		"/*movaps	%%xmm1,0x050(%%rcx)	*/			\n\t	/*	movaps		%%xmm9,0x0d0(%%rcx)		*/	\n\t"\
		"mulpd		(%%rdi),%%xmm6						\n\t		mulpd		(%%rdi),%%xmm14				\n\t"\
		"mulpd		(%%rdi),%%xmm7						\n\t		mulpd		(%%rdi),%%xmm15				\n\t"\
		"addpd		%%xmm0,%%xmm6						\n\t		addpd		%%xmm8,%%xmm14				\n\t"\
		"addpd		%%xmm1,%%xmm7						\n\t		addpd		%%xmm9,%%xmm15				\n\t"\
		"/*movaps	%%xmm6,     (%%rcx)	*/			\n\t		movaps		%%xmm14,0x080(%%rcx)		\n\t"\
		"/*movaps	%%xmm7,0x010(%%rcx)	*/			\n\t		movaps		%%xmm15,0x090(%%rcx)		\n\t"\
		"subpd		%%xmm5,%%xmm2						\n\t		subpd		%%xmm13,%%xmm10				\n\t"\
		"subpd		%%xmm4,%%xmm3						\n\t		subpd		%%xmm12,%%xmm11				\n\t"\
		"/*movaps	%%xmm2,0x020(%%rcx)	*/			\n\t		mulpd		(%%rdi),%%xmm13				\n\t"\
		"/*movaps	%%xmm3,0x070(%%rcx)	*/			\n\t		mulpd		(%%rdi),%%xmm12				\n\t"\
		"mulpd		(%%rdi),%%xmm5						\n\t		addpd		%%xmm10,%%xmm13				\n\t"\
		"mulpd		(%%rdi),%%xmm4						\n\t		addpd		%%xmm11,%%xmm12				\n\t"\
		"addpd		%%xmm2,%%xmm5						\n\t		movaps		%%xmm10,%%xmm14				\n\t"\
		"addpd		%%xmm3,%%xmm4						\n\t		movaps		%%xmm13,%%xmm15				\n\t"\
		"/*movaps	%%xmm5,0x060(%%rcx)	*/			\n\t		subpd		%%xmm12,%%xmm10				\n\t"\
		"/*movaps	%%xmm4,0x030(%%rcx)	*/			\n\t		subpd		%%xmm11,%%xmm13				\n\t"\
		"																addpd		%%xmm12,%%xmm14				\n\t"\
		"																addpd		%%xmm11,%%xmm15				\n\t"\
		"																mulpd		(%%rsi),%%xmm10	/* isrt2 */	\n\t"\
		"																mulpd		(%%rsi),%%xmm13				\n\t"\
		"																mulpd		(%%rsi),%%xmm14				\n\t"\
		"																mulpd		(%%rsi),%%xmm15				\n\t"\
		"															/*	movaps		%%xmm10,0x0a0(%%rcx)	*/	\n\t"\
		"															/*	movaps		%%xmm13,0x0e0(%%rcx)	*/	\n\t"\
		"															/*	movaps		%%xmm14,0x0b0(%%rcx)	*/	\n\t"\
		"															/*	movaps		%%xmm15,0x0f0(%%rcx)	*/	\n\t"\
		/*****	SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r10)	*****/\
		"																movaps		0x080(%%rcx),%%xmm11		\n\t"\
		"																movaps		0x090(%%rcx),%%xmm12		\n\t"\
		"subpd		%%xmm11,%%xmm6						\n\t		subpd		%%xmm10,%%xmm2			\n\t"\
		"subpd		%%xmm9,%%xmm0						\n\t		subpd		%%xmm15,%%xmm5			\n\t"\
		"subpd		%%xmm12,%%xmm7						\n\t		subpd		%%xmm14,%%xmm4			\n\t"\
		"subpd		%%xmm8,%%xmm1						\n\t		subpd		%%xmm13,%%xmm3			\n\t"\
		"mulpd		(%%rdi),%%xmm11						\n\t		mulpd		(%%rdi),%%xmm10		\n\t"\
		"mulpd		(%%rdi),%%xmm9							\n\t		mulpd		(%%rdi),%%xmm15		\n\t"\
		"mulpd		(%%rdi),%%xmm12						\n\t		mulpd		(%%rdi),%%xmm14		\n\t"\
		"mulpd		(%%rdi),%%xmm8							\n\t		mulpd		(%%rdi),%%xmm13		\n\t"\
		"addpd		%%xmm6,%%xmm11					\n\t		addpd		%%xmm2,%%xmm10		\n\t"\
		"addpd		%%xmm0,%%xmm9						\n\t		addpd		%%xmm5,%%xmm15		\n\t"\
		"addpd		%%xmm7,%%xmm12					\n\t		addpd		%%xmm4,%%xmm14		\n\t"\
		"addpd		%%xmm1,%%xmm8						\n\t		addpd		%%xmm3,%%xmm13		\n\t"\
		"movaps		%%xmm6,0x080(%%rcx)				\n\t		movaps		%%xmm2,0x0a0(%%rcx)	\n\t"\
		"movaps		%%xmm0,0x040(%%rcx)				\n\t		movaps		%%xmm5,0x060(%%rcx)	\n\t"\
		"movaps		%%xmm7,0x090(%%rcx)				\n\t		movaps		%%xmm4,0x0b0(%%rcx)	\n\t"\
		"movaps		%%xmm1,0x0d0(%%rcx)				\n\t		movaps		%%xmm3,0x0f0(%%rcx)	\n\t"\
		"movaps		%%xmm11,     (%%rcx)				\n\t		movaps		%%xmm10,0x020(%%rcx)	\n\t"\
		"movaps		%%xmm9,0x0c0(%%rcx)				\n\t		movaps		%%xmm15,0x0e0(%%rcx)	\n\t"\
		"movaps		%%xmm12,0x010(%%rcx)				\n\t		movaps		%%xmm14,0x030(%%rcx)	\n\t"\
		"movaps		%%xmm8,0x050(%%rcx)				\n\t		movaps		%%xmm13,0x070(%%rcx)	\n\t"\
		/************************************************************************************************************/\
		/* Next 2 blocks operate on odd-indexed elements from the unpck*pd commands which we stored to temporaries: */\
		/************************************************************************************************************/\
		/*...Block 3:	*/\
		/*	SSE2_RADIX4_DIF_4TWIDDLE         (r20,r24,r22,r26,r20,c01) */\
		"addq		$0x100,%%rcx		\n\t"	/*****	SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO(r28,r2C,r2A,r2E,r28,c05)	*****/\
		"movq		%[__c01],%%rbx						\n\t	/*	movq		%[__c05],%%rbx	*/		\n\t"\
		"movq		%%rcx,%%rax						\n\t	/*	addq		$0x040,%%rax	*/		\n\t"\
		"addq		$0x020,%%rcx						\n\t	/*	addq		$0x040,%%rcx	*/		\n\t"\
		"movaps			 (%%rax),%%xmm0	\n\t	movq %%rax,%%rdx \n\t	movaps		0x080(%%rax),%%xmm8			\n\t"\
		"movaps			 (%%rcx),%%xmm4						\n\t		movaps		0x080(%%rcx),%%xmm12		\n\t"\
		"movaps		0x010(%%rax),%%xmm1						\n\t		movaps		0x090(%%rax),%%xmm9			\n\t"\
		"movaps		0x010(%%rcx),%%xmm5						\n\t		movaps		0x090(%%rcx),%%xmm13		\n\t"\
		"movaps			 (%%rbx),%%xmm6						\n\t		movaps		0x080(%%rbx),%%xmm14		\n\t"\
		"movaps		0x010(%%rbx),%%xmm7						\n\t		movaps		0x090(%%rbx),%%xmm15		\n\t"\
		"movaps		%%xmm0,%%xmm2						\n\t		movaps		%%xmm8,%%xmm10		\n\t"\
		"movaps		%%xmm1,%%xmm3						\n\t		movaps		%%xmm9,%%xmm11		\n\t"\
		"mulpd		%%xmm6,%%xmm0						\n\t		mulpd		%%xmm14,%%xmm8			\n\t"\
		"mulpd		%%xmm6,%%xmm1						\n\t		mulpd		%%xmm14,%%xmm9			\n\t"\
		"mulpd		%%xmm7,%%xmm2						\n\t		mulpd		%%xmm15,%%xmm10		\n\t"\
		"mulpd		%%xmm7,%%xmm3						\n\t		mulpd		%%xmm15,%%xmm11		\n\t"\
		"movaps		%%xmm4,%%xmm6						\n\t		movaps		%%xmm12,%%xmm14		\n\t"\
		"addpd		%%xmm2,%%xmm1						\n\t		addpd		%%xmm10,%%xmm9			\n\t"\
		"movaps		%%xmm5,%%xmm7						\n\t		movaps		%%xmm13,%%xmm15		\n\t"\
		"mulpd		0x20 (%%rbx),%%xmm4						\n\t		mulpd		0xa0 (%%rbx),%%xmm12		\n\t"\
		"subpd		%%xmm3,%%xmm0						\n\t		subpd		%%xmm11,%%xmm8			\n\t"\
		"mulpd		0x20 (%%rbx),%%xmm5						\n\t		mulpd		0xa0 (%%rbx),%%xmm13		\n\t"\
		"mulpd		0x30 (%%rbx),%%xmm6						\n\t		mulpd		0xb0 (%%rbx),%%xmm14		\n\t"\
		"movaps		%%xmm0,%%xmm2						\n\t		movaps		%%xmm8,%%xmm10		\n\t"\
		"mulpd		0x30 (%%rbx),%%xmm7						\n\t		mulpd		0xb0 (%%rbx),%%xmm15		\n\t"\
		"addpd		%%xmm6,%%xmm5						\n\t		addpd		%%xmm14,%%xmm13		\n\t"\
		"movaps		%%xmm1,%%xmm3						\n\t		movaps		%%xmm9,%%xmm11		\n\t"\
		"subpd		%%xmm7,%%xmm4						\n\t		subpd		%%xmm15,%%xmm12		\n\t"\
		"addq		$0x040,%%rcx						\n\t		addpd		%%xmm12,%%xmm8			\n\t"\
		"addq		$0x60,%%rbx						\n\t		addpd		%%xmm13,%%xmm9			\n\t"\
		"movaps			 (%%rcx),%%xmm6						\n\t		subpd		%%xmm12,%%xmm10		\n\t"\
		"movaps		0x010(%%rcx),%%xmm7						\n\t		subpd		%%xmm13,%%xmm11		\n\t"\
		"addpd		%%xmm4,%%xmm0						\n\t		movaps		0x080(%%rcx),%%xmm12		\n\t"\
		"addpd		%%xmm5,%%xmm1						\n\t		movaps		0x090(%%rcx),%%xmm13		\n\t"\
		"subpd		%%xmm4,%%xmm2						\n\t		movaps		0x080(%%rcx),%%xmm14		\n\t"\
		"subpd		%%xmm5,%%xmm3						\n\t		movaps		0x090(%%rcx),%%xmm15		\n\t"\
		"movaps		%%xmm6,%%xmm4						\n\t		mulpd		0x080(%%rbx),%%xmm12		\n\t"\
		"movaps		%%xmm7,%%xmm5						\n\t		mulpd		0x080(%%rbx),%%xmm13		\n\t"\
		"mulpd			 (%%rbx),%%xmm4						\n\t		mulpd		0x090(%%rbx),%%xmm14		\n\t"\
		"mulpd			 (%%rbx),%%xmm5						\n\t		mulpd		0x090(%%rbx),%%xmm15		\n\t"\
		"mulpd		0x010(%%rbx),%%xmm6						\n\t		addpd		%%xmm14,%%xmm13		\n\t"\
		"mulpd		0x010(%%rbx),%%xmm7						\n\t		subpd		%%xmm15,%%xmm12		\n\t"\
		"/*movq		%%rax,%%rdx		*/				\n\t		movaps		%%xmm13,0x090(%%rdx)	\n\t"\
		"addpd		%%xmm6,%%xmm5						\n\t		movaps		%%xmm12,0x080(%%rdx)	\n\t"\
		"subpd		%%xmm7,%%xmm4						\n\t	/*	subq	$0x20,%%rbx */				\n\t"\
		"movaps		%%xmm5,0x010(%%rdx)	\n\t	addq	$0x040,%%rax								\n\t"\
		"movaps		%%xmm4,     (%%rdx)	\n\t	subq	$0x20,%%rbx								\n\t"\
		"movaps			 (%%rax),%%xmm4						\n\t		movaps		0x080(%%rax),%%xmm12		\n\t"\
		"movaps		0x010(%%rax),%%xmm5						\n\t		movaps		0x090(%%rax),%%xmm13		\n\t"\
		"movaps		%%xmm4,%%xmm6						\n\t		movaps		0x080(%%rax),%%xmm14		\n\t"\
		"movaps		%%xmm5,%%xmm7						\n\t		movaps		0x090(%%rax),%%xmm15		\n\t"\
		"mulpd			 (%%rbx),%%xmm4						\n\t		mulpd		0x080(%%rbx),%%xmm12		\n\t"\
		"mulpd			 (%%rbx),%%xmm5						\n\t		mulpd		0x080(%%rbx),%%xmm13		\n\t"\
		"mulpd		0x010(%%rbx),%%xmm6						\n\t		mulpd		0x090(%%rbx),%%xmm14		\n\t"\
		"mulpd		0x010(%%rbx),%%xmm7						\n\t		mulpd		0x090(%%rbx),%%xmm15		\n\t"\
		"addpd		%%xmm6,%%xmm5						\n\t		addpd		%%xmm14,%%xmm13		\n\t"\
		"subpd		%%xmm7,%%xmm4						\n\t		subpd		%%xmm15,%%xmm12		\n\t"\
		"movaps		%%xmm5,%%xmm7						\n\t		movaps		%%xmm13,%%xmm15		\n\t"\
		"movaps		%%xmm4,%%xmm6						\n\t		movaps		%%xmm12,%%xmm14		\n\t"\
		"subpd			 (%%rdx),%%xmm4						\n\t		subpd		0x080(%%rdx),%%xmm12		\n\t"\
		"subpd		0x010(%%rdx),%%xmm5						\n\t		subpd		0x090(%%rdx),%%xmm13		\n\t"\
		"addpd			 (%%rdx),%%xmm6						\n\t		addpd		0x080(%%rdx),%%xmm14		\n\t"\
		"addpd		0x010(%%rdx),%%xmm7						\n\t		addpd		0x090(%%rdx),%%xmm15		\n\t"\
		"subpd		%%xmm6,%%xmm0						\n\t		subpd		%%xmm14,%%xmm8			\n\t"\
		"subpd		%%xmm5,%%xmm2						\n\t		subpd		%%xmm15,%%xmm9			\n\t"\
		"subpd		%%xmm7,%%xmm1						\n\t	/*	movaps		%%xmm8,0x0c0(%%rdx)*/	\n\t"\
		"subpd		%%xmm4,%%xmm3						\n\t		subpd		%%xmm13,%%xmm10		\n\t"\
		"/*movaps		%%xmm0,0x040(%%rdx)	*/		\n\t	/*	movaps		%%xmm9,0x0d0(%%rdx)*/	\n\t"\
		"/*movaps		%%xmm2,0x020(%%rdx)	*/		\n\t		subpd		%%xmm12,%%xmm11		\n\t"\
		"/*movaps		%%xmm1,0x050(%%rdx)	*/		\n\t		mulpd		(%%rdi),%%xmm14		\n\t"\
		"/*movaps		%%xmm3,0x070(%%rdx)	*/		\n\t		mulpd		(%%rdi),%%xmm13		\n\t"\
		"mulpd		(%%rdi),%%xmm6						\n\t		mulpd		(%%rdi),%%xmm15		\n\t"\
		"mulpd		(%%rdi),%%xmm5						\n\t		mulpd		(%%rdi),%%xmm12		\n\t"\
		"mulpd		(%%rdi),%%xmm7						\n\t		addpd		%%xmm8,%%xmm14		\n\t"\
		"mulpd		(%%rdi),%%xmm4						\n\t		addpd		%%xmm10,%%xmm13		\n\t"\
		"addpd		%%xmm0,%%xmm6						\n\t		addpd		%%xmm9,%%xmm15		\n\t"\
		"addpd		%%xmm2,%%xmm5						\n\t		addpd		%%xmm11,%%xmm12		\n\t"\
		"addpd		%%xmm1,%%xmm7						\n\t		movaps		%%xmm14,0x080(%%rdx)	\n\t"\
		"addpd		%%xmm3,%%xmm4						\n\t		movaps		%%xmm15,0x090(%%rdx)	\n\t"\
		"/*movaps		%%xmm6,     (%%rdx)	*/		\n\t		movaps		%%xmm10,%%xmm14		\n\t"\
		"/*movaps		%%xmm5,0x060(%%rdx)	*/		\n\t		movaps		%%xmm13,%%xmm15		\n\t"\
		"/*movaps		%%xmm7,0x010(%%rdx)	*/		\n\t		subpd		%%xmm12,%%xmm10		\n\t"\
		"/*movaps		%%xmm4,0x030(%%rdx)	*/		\n\t		subpd		%%xmm11,%%xmm13		\n\t"\
		"													\n\t		addpd		%%xmm12,%%xmm14		\n\t"\
		"													\n\t		addpd		%%xmm11,%%xmm15		\n\t"\
		"													\n\t		mulpd		(%%rsi),%%xmm10		\n\t"\
		"													\n\t		mulpd		(%%rsi),%%xmm13		\n\t"\
		"													\n\t		mulpd		(%%rsi),%%xmm14		\n\t"\
		"													\n\t		mulpd		(%%rsi),%%xmm15		\n\t"\
		"													\n\t	/*	movaps		%%xmm10,0x0a0(%%rdx)*/	\n\t"\
		"													\n\t	/*	movaps		%%xmm13,0x0e0(%%rdx)*/	\n\t"\
		"													\n\t	/*	movaps		%%xmm14,0x0b0(%%rdx)*/	\n\t"\
		"													\n\t	/*	movaps		%%xmm15,0x0f0(%%rdx)*/	\n\t"\
		/*****	SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r20)	*****/\
		"																movaps		0x080(%%rdx),%%xmm11		\n\t"\
		"																movaps		0x090(%%rdx),%%xmm12		\n\t"\
		"subpd		%%xmm11,%%xmm6						\n\t		subpd		%%xmm10,%%xmm2			\n\t"\
		"subpd		%%xmm9,%%xmm0						\n\t		subpd		%%xmm15,%%xmm5			\n\t"\
		"subpd		%%xmm12,%%xmm7						\n\t		subpd		%%xmm14,%%xmm4			\n\t"\
		"subpd		%%xmm8,%%xmm1						\n\t		subpd		%%xmm13,%%xmm3			\n\t"\
		"mulpd		(%%rdi),%%xmm11					\n\t		mulpd		(%%rdi),%%xmm10		\n\t"\
		"mulpd		(%%rdi),%%xmm9						\n\t		mulpd		(%%rdi),%%xmm15		\n\t"\
		"mulpd		(%%rdi),%%xmm12					\n\t		mulpd		(%%rdi),%%xmm14		\n\t"\
		"mulpd		(%%rdi),%%xmm8						\n\t		mulpd		(%%rdi),%%xmm13		\n\t"\
		"addpd		%%xmm6,%%xmm11					\n\t		addpd		%%xmm2,%%xmm10		\n\t"\
		"addpd		%%xmm0,%%xmm9						\n\t		addpd		%%xmm5,%%xmm15		\n\t"\
		"addpd		%%xmm7,%%xmm12					\n\t		addpd		%%xmm4,%%xmm14		\n\t"\
		"addpd		%%xmm1,%%xmm8						\n\t		addpd		%%xmm3,%%xmm13		\n\t"\
		"movaps		%%xmm6,0x080(%%rdx)				\n\t		movaps		%%xmm2,0x0a0(%%rdx)	\n\t"\
		"movaps		%%xmm0,0x040(%%rdx)				\n\t		movaps		%%xmm5,0x060(%%rdx)	\n\t"\
		"movaps		%%xmm7,0x090(%%rdx)				\n\t		movaps		%%xmm4,0x0b0(%%rdx)	\n\t"\
		"movaps		%%xmm1,0x0d0(%%rdx)				\n\t		movaps		%%xmm3,0x0f0(%%rdx)	\n\t"\
		"movaps		%%xmm11,     (%%rdx)				\n\t		movaps		%%xmm10,0x020(%%rdx)	\n\t"\
		"movaps		%%xmm9,0x0c0(%%rdx)				\n\t		movaps		%%xmm15,0x0e0(%%rdx)	\n\t"\
		"movaps		%%xmm12,0x010(%%rdx)				\n\t		movaps		%%xmm14,0x030(%%rdx)	\n\t"\
		"movaps		%%xmm8,0x050(%%rdx)				\n\t		movaps		%%xmm13,0x070(%%rdx)	\n\t"\
		/*...Block 4:	*/\
		/*****	SSE2_RADIX4_DIF_4TWIDDLE         (r30,r34,r32,r36,r30,c03)	*****/\
		"movq		%[__c03],%%rbx					\n\t"\
		"movq		%[__r30],%%rax					\n\t"\
		"movq		%%rax,%%rcx		\n\t"		/*****	SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO(r38,r3C,r3A,r3E,r38,c07)	*****/\
		"addq		$0x020,%%rcx					\n\t"\
		"movaps			 (%%rax),%%xmm0	\n\t	movq %%rax,%%rdx \n\t	movaps		0x080(%%rax),%%xmm8			\n\t"\
		"movaps			 (%%rcx),%%xmm4						\n\t		movaps		0x080(%%rcx),%%xmm12		\n\t"\
		"movaps		0x010(%%rax),%%xmm1						\n\t		movaps		0x090(%%rax),%%xmm9			\n\t"\
		"movaps		0x010(%%rcx),%%xmm5						\n\t		movaps		0x090(%%rcx),%%xmm13		\n\t"\
		"movaps			 (%%rbx),%%xmm6						\n\t		movaps		0x080(%%rbx),%%xmm14		\n\t"\
		"movaps		0x010(%%rbx),%%xmm7						\n\t		movaps		0x090(%%rbx),%%xmm15		\n\t"\
		"movaps		%%xmm0,%%xmm2						\n\t		movaps		%%xmm8,%%xmm10		\n\t"\
		"movaps		%%xmm1,%%xmm3						\n\t		movaps		%%xmm9,%%xmm11		\n\t"\
		"mulpd		%%xmm6,%%xmm0						\n\t		mulpd		%%xmm14,%%xmm8			\n\t"\
		"mulpd		%%xmm6,%%xmm1						\n\t		mulpd		%%xmm14,%%xmm9			\n\t"\
		"mulpd		%%xmm7,%%xmm2						\n\t		mulpd		%%xmm15,%%xmm10		\n\t"\
		"mulpd		%%xmm7,%%xmm3						\n\t		mulpd		%%xmm15,%%xmm11		\n\t"\
		"movaps		%%xmm4,%%xmm6						\n\t		movaps		%%xmm12,%%xmm14		\n\t"\
		"addpd		%%xmm2,%%xmm1						\n\t		addpd		%%xmm10,%%xmm9			\n\t"\
		"movaps		%%xmm5,%%xmm7						\n\t		movaps		%%xmm13,%%xmm15		\n\t"\
		"mulpd		0x20 (%%rbx),%%xmm4						\n\t		mulpd		0xa0 (%%rbx),%%xmm12		\n\t"\
		"subpd		%%xmm3,%%xmm0						\n\t		subpd		%%xmm11,%%xmm8			\n\t"\
		"mulpd		0x20 (%%rbx),%%xmm5						\n\t		mulpd		0xa0 (%%rbx),%%xmm13		\n\t"\
		"mulpd		0x30 (%%rbx),%%xmm6						\n\t		mulpd		0xb0 (%%rbx),%%xmm14		\n\t"\
		"movaps		%%xmm0,%%xmm2						\n\t		movaps		%%xmm8,%%xmm10		\n\t"\
		"mulpd		0x30 (%%rbx),%%xmm7						\n\t		mulpd		0xb0 (%%rbx),%%xmm15		\n\t"\
		"addpd		%%xmm6,%%xmm5						\n\t		addpd		%%xmm14,%%xmm13		\n\t"\
		"movaps		%%xmm1,%%xmm3						\n\t		movaps		%%xmm9,%%xmm11		\n\t"\
		"subpd		%%xmm7,%%xmm4						\n\t		subpd		%%xmm15,%%xmm12		\n\t"\
		"addq		$0x040,%%rcx						\n\t		addpd		%%xmm12,%%xmm8			\n\t"\
		"addq		$0x60,%%rbx						\n\t		addpd		%%xmm13,%%xmm9			\n\t"\
		"movaps			 (%%rcx),%%xmm6						\n\t		subpd		%%xmm12,%%xmm10		\n\t"\
		"movaps		0x010(%%rcx),%%xmm7						\n\t		subpd		%%xmm13,%%xmm11		\n\t"\
		"addpd		%%xmm4,%%xmm0						\n\t		movaps		0x080(%%rcx),%%xmm12		\n\t"\
		"addpd		%%xmm5,%%xmm1						\n\t		movaps		0x090(%%rcx),%%xmm13		\n\t"\
		"subpd		%%xmm4,%%xmm2						\n\t		movaps		0x080(%%rcx),%%xmm14		\n\t"\
		"subpd		%%xmm5,%%xmm3						\n\t		movaps		0x090(%%rcx),%%xmm15		\n\t"\
		"movaps		%%xmm6,%%xmm4						\n\t		mulpd		0x080(%%rbx),%%xmm12		\n\t"\
		"movaps		%%xmm7,%%xmm5						\n\t		mulpd		0x080(%%rbx),%%xmm13		\n\t"\
		"mulpd			 (%%rbx),%%xmm4						\n\t		mulpd		0x090(%%rbx),%%xmm14		\n\t"\
		"mulpd			 (%%rbx),%%xmm5						\n\t		mulpd		0x090(%%rbx),%%xmm15		\n\t"\
		"mulpd		0x010(%%rbx),%%xmm6						\n\t		addpd		%%xmm14,%%xmm13		\n\t"\
		"mulpd		0x010(%%rbx),%%xmm7						\n\t		subpd		%%xmm15,%%xmm12		\n\t"\
		"/*movq		%%rax,%%rdx		*/				\n\t		movaps		%%xmm13,0x090(%%rdx)	\n\t"\
		"addpd		%%xmm6,%%xmm5						\n\t		movaps		%%xmm12,0x080(%%rdx)	\n\t"\
		"subpd		%%xmm7,%%xmm4						\n\t	/*	subq	$0x20,%%rbx */				\n\t"\
		"movaps		%%xmm5,0x010(%%rdx)	\n\t	addq	$0x040,%%rax								\n\t"\
		"movaps		%%xmm4,     (%%rdx)	\n\t	subq	$0x20,%%rbx								\n\t"\
		"movaps			 (%%rax),%%xmm4						\n\t		movaps		0x080(%%rax),%%xmm12		\n\t"\
		"movaps		0x010(%%rax),%%xmm5						\n\t		movaps		0x090(%%rax),%%xmm13		\n\t"\
		"movaps		%%xmm4,%%xmm6						\n\t		movaps		0x080(%%rax),%%xmm14		\n\t"\
		"movaps		%%xmm5,%%xmm7						\n\t		movaps		0x090(%%rax),%%xmm15		\n\t"\
		"mulpd			 (%%rbx),%%xmm4						\n\t		mulpd		0x080(%%rbx),%%xmm12		\n\t"\
		"mulpd			 (%%rbx),%%xmm5						\n\t		mulpd		0x080(%%rbx),%%xmm13		\n\t"\
		"mulpd		0x010(%%rbx),%%xmm6						\n\t		mulpd		0x090(%%rbx),%%xmm14		\n\t"\
		"mulpd		0x010(%%rbx),%%xmm7						\n\t		mulpd		0x090(%%rbx),%%xmm15		\n\t"\
		"addpd		%%xmm6,%%xmm5						\n\t		addpd		%%xmm14,%%xmm13		\n\t"\
		"subpd		%%xmm7,%%xmm4						\n\t		subpd		%%xmm15,%%xmm12		\n\t"\
		"movaps		%%xmm5,%%xmm7						\n\t		movaps		%%xmm13,%%xmm15		\n\t"\
		"movaps		%%xmm4,%%xmm6						\n\t		movaps		%%xmm12,%%xmm14		\n\t"\
		"subpd			 (%%rdx),%%xmm4						\n\t		subpd		0x080(%%rdx),%%xmm12		\n\t"\
		"subpd		0x010(%%rdx),%%xmm5						\n\t		subpd		0x090(%%rdx),%%xmm13		\n\t"\
		"addpd			 (%%rdx),%%xmm6						\n\t		addpd		0x080(%%rdx),%%xmm14		\n\t"\
		"addpd		0x010(%%rdx),%%xmm7						\n\t		addpd		0x090(%%rdx),%%xmm15		\n\t"\
		"subpd		%%xmm6,%%xmm0						\n\t		subpd		%%xmm14,%%xmm8			\n\t"\
		"subpd		%%xmm5,%%xmm2						\n\t		subpd		%%xmm15,%%xmm9			\n\t"\
		"subpd		%%xmm7,%%xmm1						\n\t	/*	movaps		%%xmm8,0x0c0(%%rdx)*/	\n\t"\
		"subpd		%%xmm4,%%xmm3						\n\t		subpd		%%xmm13,%%xmm10		\n\t"\
		"/*movaps		%%xmm0,0x040(%%rdx)	*/		\n\t	/*	movaps		%%xmm9,0x0d0(%%rdx)*/	\n\t"\
		"/*movaps		%%xmm2,0x020(%%rdx)	*/		\n\t		subpd		%%xmm12,%%xmm11		\n\t"\
		"/*movaps		%%xmm1,0x050(%%rdx)	*/		\n\t		mulpd		(%%rdi),%%xmm14		\n\t"\
		"/*movaps		%%xmm3,0x070(%%rdx)	*/		\n\t		mulpd		(%%rdi),%%xmm13		\n\t"\
		"mulpd		(%%rdi),%%xmm6						\n\t		mulpd		(%%rdi),%%xmm15		\n\t"\
		"mulpd		(%%rdi),%%xmm5						\n\t		mulpd		(%%rdi),%%xmm12		\n\t"\
		"mulpd		(%%rdi),%%xmm7						\n\t		addpd		%%xmm8,%%xmm14		\n\t"\
		"mulpd		(%%rdi),%%xmm4						\n\t		addpd		%%xmm10,%%xmm13		\n\t"\
		"addpd		%%xmm0,%%xmm6						\n\t		addpd		%%xmm9,%%xmm15		\n\t"\
		"addpd		%%xmm2,%%xmm5						\n\t		addpd		%%xmm11,%%xmm12		\n\t"\
		"addpd		%%xmm1,%%xmm7						\n\t		movaps		%%xmm14,0x080(%%rdx)	\n\t"\
		"addpd		%%xmm3,%%xmm4						\n\t		movaps		%%xmm15,0x090(%%rdx)	\n\t"\
		"/*movaps		%%xmm6,     (%%rdx)	*/		\n\t		movaps		%%xmm10,%%xmm14		\n\t"\
		"/*movaps		%%xmm5,0x060(%%rdx)	*/		\n\t		movaps		%%xmm13,%%xmm15		\n\t"\
		"/*movaps		%%xmm7,0x010(%%rdx)	*/		\n\t		subpd		%%xmm12,%%xmm10		\n\t"\
		"/*movaps		%%xmm4,0x030(%%rdx)	*/		\n\t		subpd		%%xmm11,%%xmm13		\n\t"\
		"													\n\t		addpd		%%xmm12,%%xmm14		\n\t"\
		"													\n\t		addpd		%%xmm11,%%xmm15		\n\t"\
		"													\n\t		mulpd		(%%rsi),%%xmm10		\n\t"\
		"													\n\t		mulpd		(%%rsi),%%xmm13		\n\t"\
		"													\n\t		mulpd		(%%rsi),%%xmm14		\n\t"\
		"													\n\t		mulpd		(%%rsi),%%xmm15		\n\t"\
		"													\n\t	/*	movaps		%%xmm10,0x0a0(%%rdx)*/	\n\t"\
		"													\n\t	/*	movaps		%%xmm13,0x0e0(%%rdx)*/	\n\t"\
		"													\n\t	/*	movaps		%%xmm14,0x0b0(%%rdx)*/	\n\t"\
		"													\n\t	/*	movaps		%%xmm15,0x0f0(%%rdx)*/	\n\t"\
		/*****	SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r20)	*****/\
		"																movaps		0x080(%%rdx),%%xmm11		\n\t"\
		"																movaps		0x090(%%rdx),%%xmm12		\n\t"\
		"subpd		%%xmm11,%%xmm6						\n\t		subpd		%%xmm10,%%xmm2			\n\t"\
		"subpd		%%xmm9,%%xmm0						\n\t		subpd		%%xmm15,%%xmm5			\n\t"\
		"subpd		%%xmm12,%%xmm7						\n\t		subpd		%%xmm14,%%xmm4			\n\t"\
		"subpd		%%xmm8,%%xmm1						\n\t		subpd		%%xmm13,%%xmm3			\n\t"\
		"mulpd		(%%rdi),%%xmm11					\n\t		mulpd		(%%rdi),%%xmm10		\n\t"\
		"mulpd		(%%rdi),%%xmm9						\n\t		mulpd		(%%rdi),%%xmm15		\n\t"\
		"mulpd		(%%rdi),%%xmm12					\n\t		mulpd		(%%rdi),%%xmm14		\n\t"\
		"mulpd		(%%rdi),%%xmm8						\n\t		mulpd		(%%rdi),%%xmm13		\n\t"\
		"addpd		%%xmm6,%%xmm11					\n\t		addpd		%%xmm2,%%xmm10		\n\t"\
		"addpd		%%xmm0,%%xmm9						\n\t		addpd		%%xmm5,%%xmm15		\n\t"\
		"addpd		%%xmm7,%%xmm12					\n\t		addpd		%%xmm4,%%xmm14		\n\t"\
		"addpd		%%xmm1,%%xmm8						\n\t		addpd		%%xmm3,%%xmm13		\n\t"\
		"movaps		%%xmm6,0x080(%%rdx)				\n\t		movaps		%%xmm2,0x0a0(%%rdx)	\n\t"\
		"movaps		%%xmm0,0x040(%%rdx)				\n\t		movaps		%%xmm5,0x060(%%rdx)	\n\t"\
		"movaps		%%xmm7,0x090(%%rdx)				\n\t		movaps		%%xmm4,0x0b0(%%rdx)	\n\t"\
		"movaps		%%xmm1,0x0d0(%%rdx)				\n\t		movaps		%%xmm3,0x0f0(%%rdx)	\n\t"\
		"movaps		%%xmm11,     (%%rdx)				\n\t		movaps		%%xmm10,0x020(%%rdx)	\n\t"\
		"movaps		%%xmm9,0x0c0(%%rdx)				\n\t		movaps		%%xmm15,0x0e0(%%rdx)	\n\t"\
		"movaps		%%xmm12,0x010(%%rdx)				\n\t		movaps		%%xmm14,0x030(%%rdx)	\n\t"\
		"movaps		%%xmm8,0x050(%%rdx)				\n\t		movaps		%%xmm13,0x070(%%rdx)	\n\t"\
		/**********************************************************************************/\
		/*...and now do eight radix-4 transforms, including the internal twiddle factots: */\
		/**********************************************************************************/\
		"movq		%[__isrt2],%%rsi		\n\t"\
		/*...Block 1: t00,t10,t20,t30	*/					/*...Block 5: t08,t18,t28,t38	*/\
		"movq		%[__r00],%%rax			\n\t	movq		$0x080,%%r10			\n\t"\
		"movq		%[__r10],%%rbx			\n\t	movq		$0x080,%%r11			\n\t"\
		"movq		%[__r20],%%rcx			\n\t	movq		$0x080,%%r12			\n\t"\
		"movq		%[__r30],%%rdx			\n\t	movq		$0x080,%%r13			\n\t"\
		"movaps			  (%%rax),%%xmm0			\n\t	addq		%%rax,%%r10			\n\t"\
		"movaps		 0x010(%%rax),%%xmm1			\n\t	addq		%%rbx,%%r11			\n\t"\
		"movaps			  (%%rbx),%%xmm2			\n\t	addq		%%rcx,%%r12			\n\t"\
		"movaps		 0x010(%%rbx),%%xmm3			\n\t	addq		%%rdx,%%r13			\n\t"\
		"addpd			   %%xmm0,%%xmm2			\n\t	movaps			  (%%r10),%%xmm8		\n\t"\
		"addpd			   %%xmm1,%%xmm3			\n\t	movaps		 0x010(%%r10),%%xmm9		\n\t"\
		"subpd			  (%%rbx),%%xmm0			\n\t	movaps			  (%%r11),%%xmm10	\n\t"\
		"subpd		 0x010(%%rbx),%%xmm1			\n\t	movaps		 0x010(%%r11),%%xmm11	\n\t"\
		"movaps			  (%%rcx),%%xmm4			\n\t	addpd			   %%xmm9,%%xmm10	\n\t"\
		"movaps		 0x010(%%rcx),%%xmm5			\n\t	addpd			   %%xmm8,%%xmm11	\n\t"\
		"movaps			  (%%rdx),%%xmm6			\n\t	subpd		 0x010(%%r11),%%xmm8		\n\t"\
		"movaps		 0x010(%%rdx),%%xmm7			\n\t	subpd			  (%%r11),%%xmm9		\n\t"\
		"addpd			   %%xmm4,%%xmm6			\n\t	movaps			  (%%r12),%%xmm12	\n\t"\
		"addpd			   %%xmm5,%%xmm7			\n\t	movaps		 0x010(%%r12),%%xmm13	\n\t"\
		"subpd			  (%%rdx),%%xmm4			\n\t	movaps			  (%%r13),%%xmm14	\n\t"\
		"subpd		 0x010(%%rdx),%%xmm5			\n\t	movaps		 0x010(%%r13),%%xmm15	\n\t"\
		"subpd		%%xmm6,%%xmm2				\n\t	subpd			   %%xmm13,%%xmm12	\n\t"\
		"subpd		%%xmm7,%%xmm3				\n\t	addpd			  (%%r12),%%xmm13	\n\t"\
		"mulpd		(%%rdi),%%xmm6				\n\t	addpd			   %%xmm15,%%xmm14	\n\t"\
		"mulpd		(%%rdi),%%xmm7				\n\t	subpd			  (%%r13),%%xmm15	\n\t"\
		"movaps		%%xmm2,      (%%rcx)		\n\t	mulpd		(%%rsi),%%xmm12		\n\t"\
		"movaps		%%xmm3, 0x010(%%rcx)		\n\t	mulpd		(%%rsi),%%xmm13		\n\t"\
		"addpd		%%xmm2,%%xmm6				\n\t	mulpd		(%%rsi),%%xmm14		\n\t"\
		"addpd		%%xmm3,%%xmm7				\n\t	mulpd		(%%rsi),%%xmm15		\n\t"\
		"movaps		%%xmm6,      (%%rax)		\n\t	subpd		%%xmm14,%%xmm12		\n\t"\
		"movaps		%%xmm7, 0x010(%%rax)		\n\t	subpd		%%xmm15,%%xmm13		\n\t"\
		"subpd		%%xmm5,%%xmm0				\n\t	mulpd		(%%rdi),%%xmm14		\n\t"\
		"subpd		%%xmm4,%%xmm1				\n\t	mulpd		(%%rdi),%%xmm15		\n\t"\
		"mulpd		(%%rdi),%%xmm5				\n\t	addpd		%%xmm12,%%xmm14		\n\t"\
		"mulpd		(%%rdi),%%xmm4				\n\t	addpd		%%xmm13,%%xmm15		\n\t"\
		"movaps		%%xmm0,      (%%rbx)		\n\t	subpd		%%xmm12,%%xmm8			\n\t"\
		"movaps		%%xmm1, 0x010(%%rdx)		\n\t	subpd		%%xmm13,%%xmm10		\n\t"\
		"addpd		%%xmm0,%%xmm5				\n\t	mulpd		(%%rdi),%%xmm12		\n\t"\
		"addpd		%%xmm1,%%xmm4				\n\t	mulpd		(%%rdi),%%xmm13		\n\t"\
		"movaps		%%xmm5,      (%%rdx)		\n\t	movaps		%%xmm8,      (%%r12)	\n\t"\
		"movaps		%%xmm4, 0x010(%%rbx)		\n\t	movaps		%%xmm10, 0x010(%%r12)	\n\t"\
		/*...Block 3: t04,t14,t24,t34	*/\
		"addq		$0x040,%%rax				\n\t	addpd		%%xmm8,%%xmm12		\n\t"\
		"addq		$0x040,%%rbx				\n\t	addpd		%%xmm10,%%xmm13		\n\t"\
		"addq		$0x040,%%rcx				\n\t	movaps		%%xmm12,      (%%r10)	\n\t"\
		"addq		$0x040,%%rdx				\n\t	movaps		%%xmm13, 0x010(%%r10)	\n\t"\
		"movaps		0x10(%%rsi),%%xmm8		/* c */	\n\t	subpd		%%xmm15,%%xmm11		\n\t"\
		"movaps		0x20(%%rsi),%%xmm10	/* s */	\n\t	subpd		%%xmm14,%%xmm9			\n\t"\
		"movaps			  (%%rcx),%%xmm4			\n\t	mulpd		(%%rdi),%%xmm15		\n\t"\
		"movaps			  (%%rdx),%%xmm6			\n\t	mulpd		(%%rdi),%%xmm14		\n\t"\
		"movaps		 0x010(%%rcx),%%xmm5			\n\t	movaps		%%xmm11,      (%%r11)	\n\t"\
		"movaps		 0x010(%%rdx),%%xmm7			\n\t	movaps		%%xmm9, 0x010(%%r13)	\n\t"\
		"movaps			   %%xmm4,%%xmm0			\n\t	addpd		%%xmm11,%%xmm15		\n\t"\
		"movaps			   %%xmm6,%%xmm2			\n\t	addpd		%%xmm9,%%xmm14		\n\t"\
		"movaps			   %%xmm5,%%xmm1			\n\t	movaps		%%xmm15,      (%%r13)	\n\t"\
		"movaps			   %%xmm7,%%xmm3			\n\t	movaps		%%xmm14, 0x010(%%r11)	\n\t"\
															/*...Block 7: t0C,t1C,t2C,t3C	*/\
		"mulpd		 %%xmm8,%%xmm4					\n\t	addq		$0x040,%%r10			\n\t"\
		"mulpd		 %%xmm10,%%xmm6					\n\t	addq		$0x040,%%r11			\n\t"\
		"mulpd		 %%xmm10,%%xmm1					\n\t	addq		$0x040,%%r12			\n\t"\
		"mulpd		 %%xmm8,%%xmm3					\n\t	addq		$0x040,%%r13			\n\t"\
		"mulpd		 %%xmm8,%%xmm5					\n\t	movaps			  (%%r12),%%xmm12	\n\t"\
		"mulpd		 %%xmm10,%%xmm7					\n\t	movaps			  (%%r13),%%xmm14	\n\t"\
		"mulpd		 %%xmm10,%%xmm0					\n\t	movaps		 0x010(%%r12),%%xmm13	\n\t"\
		"mulpd		 %%xmm8,%%xmm2					\n\t	movaps		 0x010(%%r13),%%xmm15	\n\t"\
		"subpd		%%xmm1,%%xmm4				\n\t	movaps			   %%xmm13,%%xmm9		\n\t"\
		"subpd		%%xmm3,%%xmm6				\n\t	movaps			   %%xmm15,%%xmm11	\n\t"\
		"addpd		%%xmm0,%%xmm5				\n\t	mulpd		 %%xmm10,%%xmm12			\n\t"\
		"addpd		%%xmm2,%%xmm7				\n\t	mulpd		 %%xmm8,%%xmm14			\n\t"\
		"subpd		%%xmm6,%%xmm4				\n\t	mulpd		 %%xmm8,%%xmm9				\n\t"\
		"subpd		%%xmm7,%%xmm5				\n\t	mulpd		 %%xmm10,%%xmm11			\n\t"\
		"mulpd		(%%rdi),%%xmm6				\n\t	mulpd		 %%xmm10,%%xmm13			\n\t"\
		"mulpd		(%%rdi),%%xmm7				\n\t	mulpd		 %%xmm8,%%xmm15			\n\t"\
		"addpd		%%xmm4,%%xmm6				\n\t	mulpd		(%%r12),%%xmm8				\n\t"\
		"addpd		%%xmm5,%%xmm7				\n\t	mulpd		(%%r13),%%xmm10			\n\t"\
		"movaps			  (%%rbx),%%xmm2			\n\t	subpd		%%xmm9,%%xmm12		\n\t"\
		"movaps		 0x010(%%rbx),%%xmm3			\n\t	subpd		%%xmm11,%%xmm14		\n\t"\
		"subpd		 0x010(%%rbx),%%xmm2			\n\t	addpd		%%xmm8,%%xmm13		\n\t"\
		"addpd			  (%%rbx),%%xmm3			\n\t	addpd		%%xmm10,%%xmm15		\n\t"\
		"mulpd			  (%%rsi),%%xmm2			\n\t	subpd		%%xmm14,%%xmm12		\n\t"\
		"mulpd			  (%%rsi),%%xmm3			\n\t	subpd		%%xmm15,%%xmm13		\n\t"\
		"movaps			  (%%rax),%%xmm0			\n\t	mulpd		(%%rdi),%%xmm14		\n\t"\
		"movaps		 0x010(%%rax),%%xmm1			\n\t	mulpd		(%%rdi),%%xmm15		\n\t"\
		"subpd		%%xmm2,%%xmm0				\n\t	addpd		%%xmm12,%%xmm14		\n\t"\
		"subpd		%%xmm3,%%xmm1				\n\t	addpd		%%xmm13,%%xmm15		\n\t"\
		"addpd			  (%%rax),%%xmm2			\n\t	movaps			  (%%r11),%%xmm10	\n\t"\
		"addpd		 0x010(%%rax),%%xmm3			\n\t	movaps		 0x010(%%r11),%%xmm11	\n\t"\
		"subpd		%%xmm6,%%xmm2				\n\t	addpd		 0x010(%%r11),%%xmm10	\n\t"\
		"subpd		%%xmm7,%%xmm3				\n\t	subpd			  (%%r11),%%xmm11	\n\t"\
		"mulpd		(%%rdi),%%xmm6				\n\t	mulpd			  (%%rsi),%%xmm10	\n\t"\
		"mulpd		(%%rdi),%%xmm7				\n\t	mulpd			  (%%rsi),%%xmm11	\n\t"\
		"movaps		%%xmm2,      (%%rcx)		\n\t	movaps			  (%%r10),%%xmm8		\n\t"\
		"movaps		%%xmm3, 0x010(%%rcx)		\n\t	movaps		 0x010(%%r10),%%xmm9		\n\t"\
		"addpd		%%xmm2,%%xmm6				\n\t	subpd		%%xmm10,%%xmm8			\n\t"\
		"addpd		%%xmm3,%%xmm7				\n\t	subpd		%%xmm11,%%xmm9			\n\t"\
		"movaps		%%xmm6,      (%%rax)		\n\t	addpd			  (%%r10),%%xmm10	\n\t"\
		"movaps		%%xmm7, 0x010(%%rax)		\n\t	addpd		 0x010(%%r10),%%xmm11	\n\t"\
		"subpd		%%xmm5,%%xmm0				\n\t	subpd		%%xmm12,%%xmm8			\n\t"\
		"subpd		%%xmm4,%%xmm1				\n\t	subpd		%%xmm13,%%xmm9			\n\t"\
		"mulpd		(%%rdi),%%xmm5				\n\t	mulpd		(%%rdi),%%xmm12		\n\t"\
		"mulpd		(%%rdi),%%xmm4				\n\t	mulpd		(%%rdi),%%xmm13		\n\t"\
		"movaps		%%xmm0,      (%%rbx)		\n\t	movaps		%%xmm8,      (%%r12)	\n\t"\
		"movaps		%%xmm1, 0x010(%%rdx)		\n\t	movaps		%%xmm9, 0x010(%%r12)	\n\t"\
		"addpd		%%xmm0,%%xmm5				\n\t	addpd		%%xmm8,%%xmm12		\n\t"\
		"addpd		%%xmm1,%%xmm4				\n\t	addpd		%%xmm9,%%xmm13		\n\t"\
		"movaps		%%xmm5,      (%%rdx)		\n\t	movaps		%%xmm12,      (%%r10)	\n\t"\
		"movaps		%%xmm4, 0x010(%%rbx)		\n\t	movaps		%%xmm13, 0x010(%%r10)	\n\t"\
		/*...Block 2: t02,t12,t22,t32	*/\
		"subq		$0x020,%%rax				\n\t	subpd		%%xmm15,%%xmm10		\n\t"\
		"subq		$0x020,%%rbx				\n\t	subpd		%%xmm14,%%xmm11		\n\t"\
		"subq		$0x020,%%rcx				\n\t	mulpd		(%%rdi),%%xmm15		\n\t"\
		"subq		$0x020,%%rdx				\n\t	mulpd		(%%rdi),%%xmm14		\n\t"\
		"addq		$0x30,%%rsi	/* cc1 */	\n\t	movaps		%%xmm10,      (%%r11)	\n\t"\
		"movaps			  (%%rcx),%%xmm4			\n\t	movaps		%%xmm11, 0x010(%%r13)	\n\t"\
		"movaps			  (%%rdx),%%xmm6			\n\t	addpd		%%xmm10,%%xmm15		\n\t"\
		"movaps		 0x010(%%rcx),%%xmm5			\n\t	addpd		%%xmm11,%%xmm14		\n\t"\
		"movaps		 0x010(%%rdx),%%xmm7			\n\t	movaps		%%xmm15,      (%%r13)	\n\t"\
		"movaps			  (%%rcx),%%xmm0			\n\t	movaps		%%xmm14, 0x010(%%r11)	\n\t"\
															/*...Block 6: t0A,t1A,t2A,t3A	*/\
		"movaps			  (%%rdx),%%xmm2			\n\t	subq		$0x020,%%r10			\n\t"\
		"movaps		 0x010(%%rcx),%%xmm1			\n\t	subq		$0x020,%%r11			\n\t"\
		"movaps		 0x010(%%rdx),%%xmm3			\n\t	subq		$0x020,%%r12			\n\t"\
		"mulpd			  (%%rsi),%%xmm4			\n\t	subq		$0x020,%%r13			\n\t"\
		"mulpd		 0x20 (%%rsi),%%xmm6			\n\t	movaps			  (%%r12),%%xmm12	\n\t"\
		"mulpd		 0x010(%%rsi),%%xmm1			\n\t	movaps			  (%%r13),%%xmm14	\n\t"\
		"mulpd		 0x30 (%%rsi),%%xmm3			\n\t	movaps		 0x010(%%r12),%%xmm13	\n\t"\
		"mulpd			  (%%rsi),%%xmm5			\n\t	movaps		 0x010(%%r13),%%xmm15	\n\t"\
		"mulpd		 0x20 (%%rsi),%%xmm7			\n\t	movaps			  (%%r12),%%xmm8		\n\t"\
		"mulpd		 0x010(%%rsi),%%xmm0			\n\t	movaps			  (%%r13),%%xmm10	\n\t"\
		"mulpd		 0x30 (%%rsi),%%xmm2			\n\t	movaps		 0x010(%%r12),%%xmm9		\n\t"\
		"subpd		%%xmm1,%%xmm4				\n\t	movaps		 0x010(%%r13),%%xmm11	\n\t"\
		"subpd		%%xmm3,%%xmm6				\n\t	mulpd		 0x30 (%%rsi),%%xmm12	\n\t"\
		"addpd		%%xmm0,%%xmm5				\n\t	mulpd			  (%%rsi),%%xmm14	\n\t"\
		"addpd		%%xmm2,%%xmm7				\n\t	mulpd		 0x20 (%%rsi),%%xmm9		\n\t"\
		"subpd		%%xmm6,%%xmm4				\n\t	mulpd		 0x010(%%rsi),%%xmm11	\n\t"\
		"subpd		%%xmm7,%%xmm5				\n\t	mulpd		 0x30 (%%rsi),%%xmm13	\n\t"\
		"mulpd		(%%rdi),%%xmm6				\n\t	mulpd			  (%%rsi),%%xmm15	\n\t"\
		"mulpd		(%%rdi),%%xmm7				\n\t	mulpd		 0x20 (%%rsi),%%xmm8		\n\t"\
		"addpd		%%xmm4,%%xmm6				\n\t	mulpd		 0x010(%%rsi),%%xmm10	\n\t"\
		"addpd		%%xmm5,%%xmm7				\n\t	subpd		%%xmm9,%%xmm12		\n\t"\
		"movaps			  (%%rbx),%%xmm2			\n\t	addpd		%%xmm11,%%xmm14		\n\t"\
		"movaps		 0x010(%%rbx),%%xmm0			\n\t	addpd		%%xmm8,%%xmm13		\n\t"\
		"movaps			  (%%rbx),%%xmm1			\n\t	subpd		%%xmm10,%%xmm15		\n\t"\
		"movaps		 0x010(%%rbx),%%xmm3			\n\t	subpd		%%xmm14,%%xmm12		\n\t"\
		"mulpd		-0x020(%%rsi),%%xmm2			\n\t	subpd		%%xmm15,%%xmm13		\n\t"\
		"mulpd		-0x010(%%rsi),%%xmm0			\n\t	mulpd		(%%rdi),%%xmm14		\n\t"\
		"mulpd		-0x020(%%rsi),%%xmm3			\n\t	mulpd		(%%rdi),%%xmm15		\n\t"\
		"mulpd		-0x010(%%rsi),%%xmm1			\n\t	addpd		%%xmm12,%%xmm14		\n\t"\
		"subpd		%%xmm0,%%xmm2				\n\t	addpd		%%xmm13,%%xmm15		\n\t"\
		"addpd		%%xmm1,%%xmm3				\n\t	movaps			  (%%r11),%%xmm10	\n\t"\
		"movaps			  (%%rax),%%xmm0			\n\t	movaps		 0x010(%%r11),%%xmm8		\n\t"\
		"movaps		 0x010(%%rax),%%xmm1			\n\t	movaps			  (%%r11),%%xmm9		\n\t"\
		"subpd		%%xmm2,%%xmm0				\n\t	movaps		 0x010(%%r11),%%xmm11	\n\t"\
		"subpd		%%xmm3,%%xmm1				\n\t	mulpd		-0x010(%%rsi),%%xmm10	\n\t"\
		"addpd			  (%%rax),%%xmm2			\n\t	mulpd		-0x020(%%rsi),%%xmm8		\n\t"\
		"addpd		 0x010(%%rax),%%xmm3			\n\t	mulpd		-0x010(%%rsi),%%xmm11	\n\t"\
		"subpd		%%xmm6,%%xmm2				\n\t	mulpd		-0x020(%%rsi),%%xmm9		\n\t"\
		"subpd		%%xmm7,%%xmm3				\n\t	addpd		%%xmm8,%%xmm10		\n\t"\
		"mulpd		(%%rdi),%%xmm6				\n\t	subpd		%%xmm9,%%xmm11		\n\t"\
		"mulpd		(%%rdi),%%xmm7				\n\t	movaps			  (%%r10),%%xmm8		\n\t"\
		"movaps		%%xmm2,      (%%rcx)		\n\t	movaps		 0x010(%%r10),%%xmm9		\n\t"\
		"movaps		%%xmm3, 0x010(%%rcx)		\n\t	subpd		%%xmm10,%%xmm8			\n\t"\
		"addpd		%%xmm2,%%xmm6				\n\t	subpd		%%xmm11,%%xmm9			\n\t"\
		"addpd		%%xmm3,%%xmm7				\n\t	addpd			  (%%r10),%%xmm10	\n\t"\
		"movaps		%%xmm6,      (%%rax)		\n\t	addpd		 0x010(%%r10),%%xmm11	\n\t"\
		"movaps		%%xmm7, 0x010(%%rax)		\n\t	subpd		%%xmm12,%%xmm8			\n\t"\
		"subpd		%%xmm5,%%xmm0				\n\t	subpd		%%xmm13,%%xmm9			\n\t"\
		"subpd		%%xmm4,%%xmm1				\n\t	mulpd		(%%rdi),%%xmm12		\n\t"\
		"mulpd		(%%rdi),%%xmm5				\n\t	mulpd		(%%rdi),%%xmm13		\n\t"\
		"mulpd		(%%rdi),%%xmm4				\n\t	movaps		%%xmm8,      (%%r12)	\n\t"\
		"movaps		%%xmm0,      (%%rbx)		\n\t	movaps		%%xmm9, 0x010(%%r12)	\n\t"\
		"movaps		%%xmm1, 0x010(%%rdx)		\n\t	addpd		%%xmm8,%%xmm12		\n\t"\
		"addpd		%%xmm0,%%xmm5				\n\t	addpd		%%xmm9,%%xmm13		\n\t"\
		"addpd		%%xmm1,%%xmm4				\n\t	movaps		%%xmm12,      (%%r10)	\n\t"\
		"movaps		%%xmm5,      (%%rdx)		\n\t	movaps		%%xmm13, 0x010(%%r10)	\n\t"\
		"movaps		%%xmm4, 0x010(%%rbx)		\n\t	subpd		%%xmm15,%%xmm10		\n\t"\
		/*...Block 4: t06,t16,t26,t36	*/\
		"addq		$0x040,%%rax				\n\t	subpd		%%xmm14,%%xmm11		\n\t"\
		"addq		$0x040,%%rbx				\n\t	addpd		%%xmm15,%%xmm15		\n\t"\
		"addq		$0x040,%%rcx				\n\t	addpd		%%xmm14,%%xmm14		\n\t"\
		"addq		$0x040,%%rdx				\n\t	movaps		%%xmm10,      (%%r11)	\n\t"\
		"movaps			  (%%rcx),%%xmm4			\n\t	movaps		%%xmm11, 0x010(%%r13)	\n\t"\
		"movaps			  (%%rdx),%%xmm6			\n\t	addpd		%%xmm10,%%xmm15		\n\t"\
		"movaps		 0x010(%%rcx),%%xmm5			\n\t	addpd		%%xmm11,%%xmm14		\n\t"\
		"movaps		 0x010(%%rdx),%%xmm7			\n\t	movaps		%%xmm15,      (%%r13)	\n\t"\
		"movaps			  (%%rcx),%%xmm0			\n\t	movaps		%%xmm14, 0x010(%%r11)	\n\t"\
															/*...Block 8: t0E,t1E,t2E,t3E	*/\
		"movaps			  (%%rdx),%%xmm2			\n\t	addq		$0x040,%%r10			\n\t"\
		"movaps		 0x010(%%rcx),%%xmm1			\n\t	addq		$0x040,%%r11			\n\t"\
		"movaps		 0x010(%%rdx),%%xmm3			\n\t	addq		$0x040,%%r12			\n\t"\
		"mulpd		 0x20 (%%rsi),%%xmm4			\n\t	addq		$0x040,%%r13			\n\t"\
		"mulpd		 0x010(%%rsi),%%xmm6			\n\t	movaps			  (%%r12),%%xmm12	\n\t"\
		"mulpd		 0x30 (%%rsi),%%xmm1			\n\t	movaps			  (%%r13),%%xmm14	\n\t"\
		"mulpd			  (%%rsi),%%xmm3			\n\t	movaps		 0x010(%%r12),%%xmm13	\n\t"\
		"mulpd		 0x20 (%%rsi),%%xmm5			\n\t	movaps		 0x010(%%r13),%%xmm15	\n\t"\
		"mulpd		 0x010(%%rsi),%%xmm7			\n\t	movaps			  (%%r12),%%xmm8		\n\t"\
		"mulpd		 0x30 (%%rsi),%%xmm0			\n\t	movaps			  (%%r13),%%xmm10	\n\t"\
		"mulpd			  (%%rsi),%%xmm2			\n\t	movaps		 0x010(%%r12),%%xmm9		\n\t"\
		"subpd		%%xmm1,%%xmm4				\n\t	movaps		 0x010(%%r13),%%xmm11	\n\t"\
		"addpd		%%xmm3,%%xmm6				\n\t	mulpd		 0x010(%%rsi),%%xmm12	\n\t"\
		"addpd		%%xmm0,%%xmm5				\n\t	mulpd		 0x30 (%%rsi),%%xmm14	\n\t"\
		"subpd		%%xmm2,%%xmm7				\n\t	mulpd			  (%%rsi),%%xmm9		\n\t"\
		"subpd		%%xmm6,%%xmm4				\n\t	mulpd		 0x20 (%%rsi),%%xmm11	\n\t"\
		"subpd		%%xmm7,%%xmm5				\n\t	mulpd		 0x010(%%rsi),%%xmm13	\n\t"\
		"mulpd		(%%rdi),%%xmm6				\n\t	mulpd		 0x30 (%%rsi),%%xmm15	\n\t"\
		"mulpd		(%%rdi),%%xmm7				\n\t	mulpd			  (%%rsi),%%xmm8		\n\t"\
		"addpd		%%xmm4,%%xmm6				\n\t	mulpd		 0x20 (%%rsi),%%xmm10	\n\t"\
		"addpd		%%xmm5,%%xmm7				\n\t	subpd		%%xmm9,%%xmm12		\n\t"\
		"movaps			  (%%rbx),%%xmm2			\n\t	subpd		%%xmm11,%%xmm14		\n\t"\
		"movaps		 0x010(%%rbx),%%xmm0			\n\t	addpd		%%xmm8,%%xmm13		\n\t"\
		"movaps			  (%%rbx),%%xmm1			\n\t	addpd		%%xmm10,%%xmm15		\n\t"\
		"movaps		 0x010(%%rbx),%%xmm3			\n\t	subpd		%%xmm14,%%xmm12		\n\t"\
		"mulpd		-0x010(%%rsi),%%xmm2			\n\t	subpd		%%xmm15,%%xmm13		\n\t"\
		"mulpd		-0x020(%%rsi),%%xmm0			\n\t	mulpd		(%%rdi),%%xmm14		\n\t"\
		"mulpd		-0x010(%%rsi),%%xmm3			\n\t	mulpd		(%%rdi),%%xmm15		\n\t"\
		"mulpd		-0x020(%%rsi),%%xmm1			\n\t	addpd		%%xmm12,%%xmm14		\n\t"\
		"subpd		%%xmm0,%%xmm2				\n\t	addpd		%%xmm13,%%xmm15		\n\t"\
		"addpd		%%xmm1,%%xmm3				\n\t	movaps			  (%%r11),%%xmm10	\n\t"\
		"movaps			  (%%rax),%%xmm0			\n\t	movaps		 0x010(%%r11),%%xmm8		\n\t"\
		"movaps		 0x010(%%rax),%%xmm1			\n\t	movaps			  (%%r11),%%xmm9		\n\t"\
		"subpd		%%xmm2,%%xmm0				\n\t	movaps		 0x010(%%r11),%%xmm11	\n\t"\
		"subpd		%%xmm3,%%xmm1				\n\t	mulpd		-0x020(%%rsi),%%xmm10	\n\t"\
		"addpd			  (%%rax),%%xmm2			\n\t	mulpd		-0x010(%%rsi),%%xmm8		\n\t"\
		"addpd		 0x010(%%rax),%%xmm3			\n\t	mulpd		-0x020(%%rsi),%%xmm11	\n\t"\
		"subpd		%%xmm4,%%xmm2				\n\t	mulpd		-0x010(%%rsi),%%xmm9		\n\t"\
		"subpd		%%xmm5,%%xmm3				\n\t	addpd		%%xmm8,%%xmm10		\n\t"\
		"mulpd		(%%rdi),%%xmm4				\n\t	subpd		%%xmm9,%%xmm11		\n\t"\
		"mulpd		(%%rdi),%%xmm5				\n\t	movaps			  (%%r10),%%xmm8		\n\t"\
		"movaps		%%xmm2,      (%%rcx)		\n\t	movaps		 0x010(%%r10),%%xmm9		\n\t"\
		"movaps		%%xmm3, 0x010(%%rcx)		\n\t	subpd		%%xmm10,%%xmm8			\n\t"\
		"addpd		%%xmm2,%%xmm4				\n\t	subpd		%%xmm11,%%xmm9			\n\t"\
		"addpd		%%xmm3,%%xmm5				\n\t	addpd			  (%%r10),%%xmm10	\n\t"\
		"movaps		%%xmm4,      (%%rax)		\n\t	addpd		 0x010(%%r10),%%xmm11	\n\t"\
		"movaps		%%xmm5, 0x010(%%rax)		\n\t	subpd		%%xmm12,%%xmm8			\n\t"\
		"subpd		%%xmm7,%%xmm0				\n\t	subpd		%%xmm13,%%xmm9			\n\t"\
		"subpd		%%xmm6,%%xmm1				\n\t	mulpd		(%%rdi),%%xmm12		\n\t"\
		"mulpd		(%%rdi),%%xmm7				\n\t	mulpd		(%%rdi),%%xmm13		\n\t"\
		"mulpd		(%%rdi),%%xmm6				\n\t	movaps		%%xmm8,      (%%r12)	\n\t"\
		"movaps		%%xmm0,      (%%rbx)		\n\t	movaps		%%xmm9, 0x010(%%r12)	\n\t"\
		"movaps		%%xmm1, 0x010(%%rdx)		\n\t	addpd		%%xmm8,%%xmm12		\n\t"\
		"addpd		%%xmm0,%%xmm7				\n\t	addpd		%%xmm9,%%xmm13		\n\t"\
		"addpd		%%xmm1,%%xmm6				\n\t	movaps		%%xmm12,      (%%r10)	\n\t"\
		"movaps		%%xmm7,      (%%rdx)		\n\t	movaps		%%xmm13, 0x010(%%r10)	\n\t"\
		"movaps		%%xmm6, 0x010(%%rbx)		\n\t	subpd		%%xmm15,%%xmm10		\n\t"\
		"													subpd		%%xmm14,%%xmm11		\n\t"\
		"													mulpd		(%%rdi),%%xmm15		\n\t"\
		"													mulpd		(%%rdi),%%xmm14		\n\t"\
		"													movaps		%%xmm10,      (%%r11)	\n\t"\
		"													movaps		%%xmm11, 0x010(%%r13)	\n\t"\
		"													addpd		%%xmm10,%%xmm15		\n\t"\
		"													addpd		%%xmm11,%%xmm14		\n\t"\
		"													movaps		%%xmm15,      (%%r13)	\n\t"\
		"													movaps		%%xmm14, 0x010(%%r11)	\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__add1] "m" (Xadd1)\
		 ,[__r00] "m" (Xr00)\
		 ,[__r10] "m" (Xr10)\
		 ,[__r20] "m" (Xr20)\
		 ,[__r30] "m" (Xr30)\
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__cc0] "m" (Xcc0)\
		 ,[__c00] "m" (Xc00)\
		 ,[__c01] "m" (Xc01)\
		 ,[__c02] "m" (Xc02)\
		 ,[__c03] "m" (Xc03)\
		 ,[__c05] "m" (Xc05)\
		 ,[__c07] "m" (Xc07)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r10","r11","r12","r13","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"		/* Clobbered registers */\
	);\
	}

	#define SSE2_RADIX32_WRAPPER_DIT(Xadd0,Xadd1,Xisrt2,Xr00,Xr08,Xr10,Xc01,Xc02,Xc04,Xc06,Xc08,Xc0A,Xc0C,Xc0E,Xc10,Xc12,Xc14,Xc16,Xc18,Xc1A,Xc1C,Xc1E)\
	{\
	__asm__ volatile (\
	/************************************************************************/\
	/*...And do an inverse DIT radix-32 pass on the squared-data blocks:	*/\
	/************************************************************************/\
		/*...Block 1: */\
		"movq		%[__isrt2],%%rsi	\n\t"\
		"movq		%[__r00],%%rax	\n\t"\
		"leaq		0x200(%%rax),%%rbx	\n\t"\
		"leaq		0x100(%%rax),%%rcx	\n\t"\
		"leaq		0x300(%%rax),%%rdx	\n\t"\
		"/*****	SSE2_RADIX4_DIT_IN_PLACE()	*****/	\n\t		/*...Block 2 has tmp-addresses offset +0x40 w.r.to Block 1:	*/\n\t"\
		"movaps			  (%%rax),%%xmm0			\n\t		movaps		 0x040(%%rax),%%xmm8 \n\t"\
		"movaps		 0x010(%%rax),%%xmm1			\n\t		movaps		 0x050(%%rax),%%xmm9 \n\t"\
		"movaps			  (%%rax),%%xmm2			\n\t		movaps		 0x040(%%rax),%%xmm10\n\t"\
		"movaps		 0x010(%%rax),%%xmm3			\n\t		movaps		 0x050(%%rax),%%xmm11\n\t"\
		"addpd			  (%%rbx),%%xmm0			\n\t		addpd		 0x040(%%rbx),%%xmm8 \n\t"\
		"addpd		 0x010(%%rbx),%%xmm1			\n\t		addpd		 0x050(%%rbx),%%xmm9 \n\t"\
		"subpd			  (%%rbx),%%xmm2			\n\t		subpd		 0x040(%%rbx),%%xmm10\n\t"\
		"subpd		 0x010(%%rbx),%%xmm3			\n\t		subpd		 0x050(%%rbx),%%xmm11\n\t"\
		"movaps			  (%%rcx),%%xmm4			\n\t		movaps		 0x040(%%rcx),%%xmm12\n\t"\
		"movaps		 0x010(%%rcx),%%xmm5			\n\t		movaps		 0x050(%%rcx),%%xmm13\n\t"\
		"movaps			  (%%rcx),%%xmm6			\n\t		movaps		 0x040(%%rcx),%%xmm14\n\t"\
		"movaps		 0x010(%%rcx),%%xmm7			\n\t		movaps		 0x050(%%rcx),%%xmm15\n\t"\
		"addpd			  (%%rdx),%%xmm4			\n\t		addpd		 0x040(%%rdx),%%xmm12\n\t"\
		"addpd		 0x010(%%rdx),%%xmm5			\n\t		addpd		 0x050(%%rdx),%%xmm13\n\t"\
		"subpd			  (%%rdx),%%xmm6			\n\t		subpd		 0x040(%%rdx),%%xmm14\n\t"\
		"subpd		 0x010(%%rdx),%%xmm7			\n\t		subpd		 0x050(%%rdx),%%xmm15\n\t"\
		"subpd		%%xmm4,%%xmm0				\n\t		subpd		%%xmm12,%%xmm8 \n\t"\
		"subpd		%%xmm5,%%xmm1				\n\t		subpd		%%xmm13,%%xmm9 \n\t"\
		"movaps		%%xmm0,      (%%rbx)		\n\t		movaps		%%xmm8 , 0x040(%%rbx)\n\t"\
		"movaps		%%xmm1, 0x010(%%rbx)		\n\t		movaps		%%xmm9 , 0x050(%%rbx)\n\t"\
		"addpd		%%xmm4,%%xmm4				\n\t		addpd		%%xmm12,%%xmm12\n\t"\
		"addpd		%%xmm5,%%xmm5				\n\t		addpd		%%xmm13,%%xmm13\n\t"\
		"addpd		%%xmm0,%%xmm4				\n\t		addpd		%%xmm8 ,%%xmm12\n\t"\
		"addpd		%%xmm1,%%xmm5				\n\t		addpd		%%xmm9 ,%%xmm13\n\t"\
		"movaps		%%xmm4,      (%%rax)		\n\t		movaps		%%xmm12, 0x040(%%rax)\n\t"\
		"movaps		%%xmm5, 0x010(%%rax)		\n\t		movaps		%%xmm13, 0x050(%%rax)\n\t"\
		"subpd		%%xmm7,%%xmm2				\n\t		subpd		%%xmm15,%%xmm10\n\t"\
		"subpd		%%xmm6,%%xmm3				\n\t		subpd		%%xmm14,%%xmm11\n\t"\
		"movaps		%%xmm2,      (%%rdx)		\n\t		movaps		%%xmm10, 0x040(%%rdx)\n\t"\
		"movaps		%%xmm3, 0x010(%%rcx)		\n\t		movaps		%%xmm11, 0x050(%%rcx)\n\t"\
		"addpd		%%xmm7,%%xmm7				\n\t		addpd		%%xmm15,%%xmm15\n\t"\
		"addpd		%%xmm6,%%xmm6				\n\t		addpd		%%xmm14,%%xmm14\n\t"\
		"addpd		%%xmm2,%%xmm7				\n\t		addpd		%%xmm10,%%xmm15\n\t"\
		"addpd		%%xmm3,%%xmm6				\n\t		addpd		%%xmm11,%%xmm14\n\t"\
		"movaps		%%xmm7,      (%%rcx)		\n\t		movaps		%%xmm15, 0x040(%%rcx)\n\t"\
		"movaps		%%xmm6, 0x010(%%rdx)		\n\t		movaps		%%xmm14, 0x050(%%rdx)\n\t"\
		"addq		$0x080,%%rax				\n\t"\
		"addq		$0x080,%%rbx				\n\t"\
		"addq		$0x080,%%rcx				\n\t"\
		"addq		$0x080,%%rdx				\n\t"\
		"movaps			  (%%rax),%%xmm0			\n\t		movaps		 0x040(%%rax),%%xmm8 \n\t"\
		"movaps		 0x010(%%rax),%%xmm1			\n\t		movaps		 0x050(%%rax),%%xmm9 \n\t"\
		"movaps			  (%%rax),%%xmm2			\n\t		movaps		 0x040(%%rax),%%xmm10\n\t"\
		"movaps		 0x010(%%rax),%%xmm3			\n\t		movaps		 0x050(%%rax),%%xmm11\n\t"\
		"addpd			  (%%rbx),%%xmm0			\n\t		addpd		 0x040(%%rbx),%%xmm8 \n\t"\
		"addpd		 0x010(%%rbx),%%xmm1			\n\t		addpd		 0x050(%%rbx),%%xmm9 \n\t"\
		"subpd			  (%%rbx),%%xmm2			\n\t		subpd		 0x040(%%rbx),%%xmm10\n\t"\
		"subpd		 0x010(%%rbx),%%xmm3			\n\t		subpd		 0x050(%%rbx),%%xmm11\n\t"\
		"movaps			  (%%rcx),%%xmm4			\n\t		movaps		 0x040(%%rcx),%%xmm12\n\t"\
		"movaps		 0x010(%%rcx),%%xmm5			\n\t		movaps		 0x050(%%rcx),%%xmm13\n\t"\
		"movaps			  (%%rcx),%%xmm6			\n\t		movaps		 0x040(%%rcx),%%xmm14\n\t"\
		"movaps		 0x010(%%rcx),%%xmm7			\n\t		movaps		 0x050(%%rcx),%%xmm15\n\t"\
		"addpd			  (%%rdx),%%xmm4			\n\t		addpd		 0x040(%%rdx),%%xmm12\n\t"\
		"addpd		 0x010(%%rdx),%%xmm5			\n\t		addpd		 0x050(%%rdx),%%xmm13\n\t"\
		"subpd			  (%%rdx),%%xmm6			\n\t		subpd		 0x040(%%rdx),%%xmm14\n\t"\
		"subpd		 0x010(%%rdx),%%xmm7			\n\t		subpd		 0x050(%%rdx),%%xmm15\n\t"\
		"subpd		%%xmm4,%%xmm0				\n\t		subpd		%%xmm12,%%xmm8 \n\t"\
		"subpd		%%xmm5,%%xmm1				\n\t		subpd		%%xmm13,%%xmm9 \n\t"\
		"movaps		%%xmm0,      (%%rbx)		\n\t		movaps		%%xmm8 , 0x040(%%rbx)\n\t"\
		"movaps		%%xmm1, 0x010(%%rbx)		\n\t		movaps		%%xmm9 , 0x050(%%rbx)\n\t"\
		"addpd		%%xmm4,%%xmm4				\n\t		addpd		%%xmm12,%%xmm12\n\t"\
		"addpd		%%xmm5,%%xmm5				\n\t		addpd		%%xmm13,%%xmm13\n\t"\
		"addpd		%%xmm0,%%xmm4				\n\t		addpd		%%xmm8 ,%%xmm12\n\t"\
		"addpd		%%xmm1,%%xmm5				\n\t		addpd		%%xmm9 ,%%xmm13\n\t"\
		"movaps		%%xmm4,      (%%rax)		\n\t		movaps		%%xmm12, 0x040(%%rax)\n\t"\
		"movaps		%%xmm5, 0x010(%%rax)		\n\t		movaps		%%xmm13, 0x050(%%rax)\n\t"\
		"subpd		%%xmm7,%%xmm2				\n\t		subpd		%%xmm15,%%xmm10\n\t"\
		"subpd		%%xmm6,%%xmm3				\n\t		subpd		%%xmm14,%%xmm11\n\t"\
		"addpd		%%xmm7,%%xmm7				\n\t		addpd		%%xmm15,%%xmm15\n\t"\
		"addpd		%%xmm6,%%xmm6				\n\t		addpd		%%xmm14,%%xmm14\n\t"\
		"addpd		%%xmm2,%%xmm7				\n\t		addpd		%%xmm10,%%xmm15\n\t"\
		"addpd		%%xmm3,%%xmm6				\n\t		addpd		%%xmm11,%%xmm14\n\t"\
		"movaps		%%xmm3,%%xmm0				\n\t		movaps		%%xmm11,%%xmm8 \n\t"\
		"movaps		%%xmm6,%%xmm1				\n\t		movaps		%%xmm14,%%xmm9 \n\t"\
		"subpd		%%xmm7,%%xmm3				\n\t		subpd		%%xmm15,%%xmm11\n\t"\
		"subpd		%%xmm2,%%xmm6				\n\t		subpd		%%xmm10,%%xmm14\n\t"\
		"addpd		%%xmm7,%%xmm0				\n\t		addpd		%%xmm15,%%xmm8 \n\t"\
		"addpd		%%xmm2,%%xmm1				\n\t		addpd		%%xmm10,%%xmm9 \n\t"\
		"mulpd			  (%%rsi),%%xmm3			\n\t		mulpd			  (%%rsi),%%xmm11\n\t"\
		"mulpd			  (%%rsi),%%xmm6			\n\t		mulpd			  (%%rsi),%%xmm14\n\t"\
		"mulpd			  (%%rsi),%%xmm0			\n\t		mulpd			  (%%rsi),%%xmm8 \n\t"\
		"mulpd			  (%%rsi),%%xmm1			\n\t		mulpd			  (%%rsi),%%xmm9 \n\t"\
		"movaps		%%xmm3, 0x010(%%rcx)		\n\t		movaps		%%xmm11, 0x050(%%rcx)\n\t"\
		"movaps		%%xmm6, 0x010(%%rdx)		\n\t		movaps		%%xmm14, 0x050(%%rdx)\n\t"\
		"movaps		%%xmm0,      (%%rcx)		\n\t		movaps		%%xmm8 , 0x040(%%rcx)\n\t"\
		"movaps		%%xmm1,      (%%rdx)		\n\t		movaps		%%xmm9 , 0x040(%%rdx)\n\t"\
		"/***************************	SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS ************************************/\n\t"\
		"movaps		-0x080(%%rax),%%xmm0			\n\t		movaps		-0x040(%%rax),%%xmm8 \n\t"\
		"movaps		-0x080(%%rbx),%%xmm4			\n\t		movaps		-0x040(%%rbx),%%xmm12\n\t"\
		"movaps		-0x070(%%rax),%%xmm1			\n\t		movaps		-0x030(%%rax),%%xmm9 \n\t"\
		"movaps		-0x070(%%rbx),%%xmm5			\n\t		movaps		-0x030(%%rbx),%%xmm13\n\t"\
		"movaps			  (%%rax),%%xmm2			\n\t		movaps		 0x040(%%rax),%%xmm10\n\t"\
		"movaps		 0x010(%%rbx),%%xmm7			\n\t		movaps		 0x050(%%rbx),%%xmm15\n\t"\
		"movaps		 0x010(%%rax),%%xmm3			\n\t		movaps		 0x050(%%rax),%%xmm11\n\t"\
		"movaps			  (%%rbx),%%xmm6			\n\t		movaps		 0x040(%%rbx),%%xmm14\n\t"\
		"subpd		%%xmm2,%%xmm0				\n\t		subpd		%%xmm10,%%xmm8 \n\t"\
		"subpd		%%xmm7,%%xmm4				\n\t		subpd		%%xmm15,%%xmm12\n\t"\
		"subpd		%%xmm3,%%xmm1				\n\t		subpd		%%xmm11,%%xmm9 \n\t"\
		"subpd		%%xmm6,%%xmm5				\n\t		subpd		%%xmm14,%%xmm13\n\t"\
		"addpd		%%xmm2,%%xmm2				\n\t		addpd		%%xmm10,%%xmm10\n\t"\
		"addpd		%%xmm7,%%xmm7				\n\t		addpd		%%xmm15,%%xmm15\n\t"\
		"addpd		%%xmm3,%%xmm3				\n\t		addpd		%%xmm11,%%xmm11\n\t"\
		"addpd		%%xmm6,%%xmm6				\n\t		addpd		%%xmm14,%%xmm14\n\t"\
		"addpd		%%xmm0,%%xmm2				\n\t		addpd		%%xmm8 ,%%xmm10\n\t"\
		"addpd		%%xmm4,%%xmm7				\n\t		addpd		%%xmm12,%%xmm15\n\t"\
		"addpd		%%xmm1,%%xmm3				\n\t		addpd		%%xmm9 ,%%xmm11\n\t"\
		"addpd		%%xmm5,%%xmm6				\n\t		addpd		%%xmm13,%%xmm14\n\t"\
		"movaps		%%xmm0,      (%%rax)		\n\t		movaps		%%xmm8 , 0x040(%%rax)\n\t"\
		"movaps		%%xmm4,      (%%rbx)		\n\t		movaps		%%xmm12, 0x040(%%rbx)\n\t"\
		"movaps		%%xmm1, 0x010(%%rax)		\n\t		movaps		%%xmm9 , 0x050(%%rax)\n\t"\
		"movaps		%%xmm5,-0x070(%%rbx)		\n\t		movaps		%%xmm13,-0x030(%%rbx)\n\t"\
		"movaps		%%xmm2,-0x080(%%rax)		\n\t		movaps		%%xmm10,-0x040(%%rax)\n\t"\
		"movaps		%%xmm7,-0x080(%%rbx)		\n\t		movaps		%%xmm15,-0x040(%%rbx)\n\t"\
		"movaps		%%xmm3,-0x070(%%rax)		\n\t		movaps		%%xmm11,-0x030(%%rax)\n\t"\
		"movaps		%%xmm6, 0x010(%%rbx)		\n\t		movaps		%%xmm14, 0x050(%%rbx)\n\t"\
		"movaps		-0x080(%%rcx),%%xmm0			\n\t		movaps		-0x040(%%rcx),%%xmm8 \n\t"\
		"movaps		-0x080(%%rdx),%%xmm4			\n\t		movaps		-0x040(%%rdx),%%xmm12\n\t"\
		"movaps		-0x070(%%rcx),%%xmm1			\n\t		movaps		-0x030(%%rcx),%%xmm9 \n\t"\
		"movaps		-0x070(%%rdx),%%xmm5			\n\t		movaps		-0x030(%%rdx),%%xmm13\n\t"\
		"movaps			  (%%rcx),%%xmm2			\n\t		movaps		 0x040(%%rcx),%%xmm10\n\t"\
		"movaps		 0x010(%%rdx),%%xmm7			\n\t		movaps		 0x050(%%rdx),%%xmm15\n\t"\
		"movaps		 0x010(%%rcx),%%xmm3			\n\t		movaps		 0x050(%%rcx),%%xmm11\n\t"\
		"movaps			  (%%rdx),%%xmm6			\n\t		movaps		 0x040(%%rdx),%%xmm14\n\t"\
		"subpd		%%xmm2,%%xmm0				\n\t		subpd		%%xmm10,%%xmm8 \n\t"\
		"subpd		%%xmm7,%%xmm4				\n\t		subpd		%%xmm15,%%xmm12\n\t"\
		"subpd		%%xmm3,%%xmm1				\n\t		subpd		%%xmm11,%%xmm9 \n\t"\
		"subpd		%%xmm6,%%xmm5				\n\t		subpd		%%xmm14,%%xmm13\n\t"\
		"addpd		%%xmm2,%%xmm2				\n\t		addpd		%%xmm10,%%xmm10\n\t"\
		"addpd		%%xmm7,%%xmm7				\n\t		addpd		%%xmm15,%%xmm15\n\t"\
		"addpd		%%xmm3,%%xmm3				\n\t		addpd		%%xmm11,%%xmm11\n\t"\
		"addpd		%%xmm6,%%xmm6				\n\t		addpd		%%xmm14,%%xmm14\n\t"\
		"addpd		%%xmm0,%%xmm2				\n\t		addpd		%%xmm8 ,%%xmm10\n\t"\
		"addpd		%%xmm4,%%xmm7				\n\t		addpd		%%xmm12,%%xmm15\n\t"\
		"addpd		%%xmm1,%%xmm3				\n\t		addpd		%%xmm9 ,%%xmm11\n\t"\
		"addpd		%%xmm5,%%xmm6				\n\t		addpd		%%xmm13,%%xmm14\n\t"\
		"movaps		%%xmm0,      (%%rcx)		\n\t		movaps		%%xmm8 , 0x040(%%rcx)\n\t"\
		"movaps		%%xmm4,      (%%rdx)		\n\t		movaps		%%xmm12, 0x040(%%rdx)\n\t"\
		"movaps		%%xmm1, 0x010(%%rcx)		\n\t		movaps		%%xmm9 , 0x050(%%rcx)\n\t"\
		"movaps		%%xmm5,-0x070(%%rdx)		\n\t		movaps		%%xmm13,-0x030(%%rdx)\n\t"\
		"movaps		%%xmm2,-0x080(%%rcx)		\n\t		movaps		%%xmm10,-0x040(%%rcx)\n\t"\
		"movaps		%%xmm7,-0x080(%%rdx)		\n\t		movaps		%%xmm15,-0x040(%%rdx)\n\t"\
		"movaps		%%xmm3,-0x070(%%rcx)		\n\t		movaps		%%xmm11,-0x030(%%rcx)\n\t"\
		"movaps		%%xmm6, 0x010(%%rdx)		\n\t		movaps		%%xmm14, 0x050(%%rdx)\n\t"\
		"\n\t"\
		"/*...Blocks 3,4 have tmp-addresses offset +0x20 w.r.to Blocks 1,2, respectively (thus +0x80-0x060 = +0x20: */\n\t"\
		"subq		$0x60,%%rax				\n\t"\
		"subq		$0x60,%%rbx				\n\t"\
		"subq		$0x60,%%rcx				\n\t"\
		"subq		$0x60,%%rdx				\n\t"\
		"movaps			  (%%rax),%%xmm0			\n\t		movaps		 0x040(%%rax),%%xmm8 \n\t"\
		"movaps		 0x010(%%rax),%%xmm1			\n\t		movaps		 0x050(%%rax),%%xmm9 \n\t"\
		"movaps			  (%%rax),%%xmm2			\n\t		movaps		 0x040(%%rax),%%xmm10\n\t"\
		"movaps		 0x010(%%rax),%%xmm3			\n\t		movaps		 0x050(%%rax),%%xmm11\n\t"\
		"addpd			  (%%rbx),%%xmm0			\n\t		addpd		 0x040(%%rbx),%%xmm8 \n\t"\
		"addpd		 0x010(%%rbx),%%xmm1			\n\t		addpd		 0x050(%%rbx),%%xmm9 \n\t"\
		"subpd			  (%%rbx),%%xmm2			\n\t		subpd		 0x040(%%rbx),%%xmm10\n\t"\
		"subpd		 0x010(%%rbx),%%xmm3			\n\t		subpd		 0x050(%%rbx),%%xmm11\n\t"\
		"movaps			  (%%rcx),%%xmm4			\n\t		movaps		 0x040(%%rcx),%%xmm12\n\t"\
		"movaps		 0x010(%%rcx),%%xmm5			\n\t		movaps		 0x050(%%rcx),%%xmm13\n\t"\
		"movaps			  (%%rcx),%%xmm6			\n\t		movaps		 0x040(%%rcx),%%xmm14\n\t"\
		"movaps		 0x010(%%rcx),%%xmm7			\n\t		movaps		 0x050(%%rcx),%%xmm15\n\t"\
		"addpd			  (%%rdx),%%xmm4			\n\t		addpd		 0x040(%%rdx),%%xmm12\n\t"\
		"addpd		 0x010(%%rdx),%%xmm5			\n\t		addpd		 0x050(%%rdx),%%xmm13\n\t"\
		"subpd			  (%%rdx),%%xmm6			\n\t		subpd		 0x040(%%rdx),%%xmm14\n\t"\
		"subpd		 0x010(%%rdx),%%xmm7			\n\t		subpd		 0x050(%%rdx),%%xmm15\n\t"\
		"subpd		%%xmm4,%%xmm0				\n\t		subpd		%%xmm12,%%xmm8 \n\t"\
		"subpd		%%xmm5,%%xmm1				\n\t		subpd		%%xmm13,%%xmm9 \n\t"\
		"movaps		%%xmm0,      (%%rbx)		\n\t		movaps		%%xmm8 , 0x040(%%rbx)\n\t"\
		"movaps		%%xmm1, 0x010(%%rbx)		\n\t		movaps		%%xmm9 , 0x050(%%rbx)\n\t"\
		"addpd		%%xmm4,%%xmm4				\n\t		addpd		%%xmm12,%%xmm12\n\t"\
		"addpd		%%xmm5,%%xmm5				\n\t		addpd		%%xmm13,%%xmm13\n\t"\
		"addpd		%%xmm0,%%xmm4				\n\t		addpd		%%xmm8 ,%%xmm12\n\t"\
		"addpd		%%xmm1,%%xmm5				\n\t		addpd		%%xmm9 ,%%xmm13\n\t"\
		"movaps		%%xmm4,      (%%rax)		\n\t		movaps		%%xmm12, 0x040(%%rax)\n\t"\
		"movaps		%%xmm5, 0x010(%%rax)		\n\t		movaps		%%xmm13, 0x050(%%rax)\n\t"\
		"subpd		%%xmm7,%%xmm2				\n\t		subpd		%%xmm15,%%xmm10\n\t"\
		"subpd		%%xmm6,%%xmm3				\n\t		subpd		%%xmm14,%%xmm11\n\t"\
		"movaps		%%xmm2,      (%%rdx)		\n\t		movaps		%%xmm10, 0x040(%%rdx)\n\t"\
		"movaps		%%xmm3, 0x010(%%rcx)		\n\t		movaps		%%xmm11, 0x050(%%rcx)\n\t"\
		"addpd		%%xmm7,%%xmm7				\n\t		addpd		%%xmm15,%%xmm15\n\t"\
		"addpd		%%xmm6,%%xmm6				\n\t		addpd		%%xmm14,%%xmm14\n\t"\
		"addpd		%%xmm2,%%xmm7				\n\t		addpd		%%xmm10,%%xmm15\n\t"\
		"addpd		%%xmm3,%%xmm6				\n\t		addpd		%%xmm11,%%xmm14\n\t"\
		"movaps		%%xmm7,      (%%rcx)		\n\t		movaps		%%xmm15, 0x040(%%rcx)\n\t"\
		"movaps		%%xmm6, 0x010(%%rdx)		\n\t		movaps		%%xmm14, 0x050(%%rdx)\n\t"\
		"addq		$0x080,%%rax				\n\t"\
		"addq		$0x080,%%rbx				\n\t"\
		"addq		$0x080,%%rcx				\n\t"\
		"addq		$0x080,%%rdx				\n\t"\
		"movaps			  (%%rax),%%xmm0			\n\t		movaps		 0x040(%%rax),%%xmm8 \n\t"\
		"movaps		 0x010(%%rax),%%xmm1			\n\t		movaps		 0x050(%%rax),%%xmm9 \n\t"\
		"movaps			  (%%rax),%%xmm2			\n\t		movaps		 0x040(%%rax),%%xmm10\n\t"\
		"movaps		 0x010(%%rax),%%xmm3			\n\t		movaps		 0x050(%%rax),%%xmm11\n\t"\
		"addpd			  (%%rbx),%%xmm0			\n\t		addpd		 0x040(%%rbx),%%xmm8 \n\t"\
		"addpd		 0x010(%%rbx),%%xmm1			\n\t		addpd		 0x050(%%rbx),%%xmm9 \n\t"\
		"subpd			  (%%rbx),%%xmm2			\n\t		subpd		 0x040(%%rbx),%%xmm10\n\t"\
		"subpd		 0x010(%%rbx),%%xmm3			\n\t		subpd		 0x050(%%rbx),%%xmm11\n\t"\
		"movaps			  (%%rcx),%%xmm4			\n\t		movaps		 0x040(%%rcx),%%xmm12\n\t"\
		"movaps		 0x010(%%rcx),%%xmm5			\n\t		movaps		 0x050(%%rcx),%%xmm13\n\t"\
		"movaps			  (%%rcx),%%xmm6			\n\t		movaps		 0x040(%%rcx),%%xmm14\n\t"\
		"movaps		 0x010(%%rcx),%%xmm7			\n\t		movaps		 0x050(%%rcx),%%xmm15\n\t"\
		"addpd			  (%%rdx),%%xmm4			\n\t		addpd		 0x040(%%rdx),%%xmm12\n\t"\
		"addpd		 0x010(%%rdx),%%xmm5			\n\t		addpd		 0x050(%%rdx),%%xmm13\n\t"\
		"subpd			  (%%rdx),%%xmm6			\n\t		subpd		 0x040(%%rdx),%%xmm14\n\t"\
		"subpd		 0x010(%%rdx),%%xmm7			\n\t		subpd		 0x050(%%rdx),%%xmm15\n\t"\
		"subpd		%%xmm4,%%xmm0				\n\t		subpd		%%xmm12,%%xmm8 \n\t"\
		"subpd		%%xmm5,%%xmm1				\n\t		subpd		%%xmm13,%%xmm9 \n\t"\
		"movaps		%%xmm0,      (%%rbx)		\n\t		movaps		%%xmm8 , 0x040(%%rbx)\n\t"\
		"movaps		%%xmm1, 0x010(%%rbx)		\n\t		movaps		%%xmm9 , 0x050(%%rbx)\n\t"\
		"addpd		%%xmm4,%%xmm4				\n\t		addpd		%%xmm12,%%xmm12\n\t"\
		"addpd		%%xmm5,%%xmm5				\n\t		addpd		%%xmm13,%%xmm13\n\t"\
		"addpd		%%xmm0,%%xmm4				\n\t		addpd		%%xmm8 ,%%xmm12\n\t"\
		"addpd		%%xmm1,%%xmm5				\n\t		addpd		%%xmm9 ,%%xmm13\n\t"\
		"movaps		%%xmm4,      (%%rax)		\n\t		movaps		%%xmm12, 0x040(%%rax)\n\t"\
		"movaps		%%xmm5, 0x010(%%rax)		\n\t		movaps		%%xmm13, 0x050(%%rax)\n\t"\
		"subpd		%%xmm7,%%xmm2				\n\t		subpd		%%xmm15,%%xmm10\n\t"\
		"subpd		%%xmm6,%%xmm3				\n\t		subpd		%%xmm14,%%xmm11\n\t"\
		"addpd		%%xmm7,%%xmm7				\n\t		addpd		%%xmm15,%%xmm15\n\t"\
		"addpd		%%xmm6,%%xmm6				\n\t		addpd		%%xmm14,%%xmm14\n\t"\
		"addpd		%%xmm2,%%xmm7				\n\t		addpd		%%xmm10,%%xmm15\n\t"\
		"addpd		%%xmm3,%%xmm6				\n\t		addpd		%%xmm11,%%xmm14\n\t"\
		"movaps		%%xmm3,%%xmm0				\n\t		movaps		%%xmm11,%%xmm8 \n\t"\
		"movaps		%%xmm6,%%xmm1				\n\t		movaps		%%xmm14,%%xmm9 \n\t"\
		"subpd		%%xmm7,%%xmm3				\n\t		subpd		%%xmm15,%%xmm11\n\t"\
		"subpd		%%xmm2,%%xmm6				\n\t		subpd		%%xmm10,%%xmm14\n\t"\
		"addpd		%%xmm7,%%xmm0				\n\t		addpd		%%xmm15,%%xmm8 \n\t"\
		"addpd		%%xmm2,%%xmm1				\n\t		addpd		%%xmm10,%%xmm9 \n\t"\
		"mulpd			  (%%rsi),%%xmm3			\n\t		mulpd			  (%%rsi),%%xmm11\n\t"\
		"mulpd			  (%%rsi),%%xmm6			\n\t		mulpd			  (%%rsi),%%xmm14\n\t"\
		"mulpd			  (%%rsi),%%xmm0			\n\t		mulpd			  (%%rsi),%%xmm8 \n\t"\
		"mulpd			  (%%rsi),%%xmm1			\n\t		mulpd			  (%%rsi),%%xmm9 \n\t"\
		"movaps		%%xmm3, 0x010(%%rcx)		\n\t		movaps		%%xmm11, 0x050(%%rcx)\n\t"\
		"movaps		%%xmm6, 0x010(%%rdx)		\n\t		movaps		%%xmm14, 0x050(%%rdx)\n\t"\
		"movaps		%%xmm0,      (%%rcx)		\n\t		movaps		%%xmm8 , 0x040(%%rcx)\n\t"\
		"movaps		%%xmm1,      (%%rdx)		\n\t		movaps		%%xmm9 , 0x040(%%rdx)\n\t"\
		"/***************************	SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS ************************************/\n\t"\
		"movaps		-0x080(%%rax),%%xmm0			\n\t		movaps		-0x040(%%rax),%%xmm8 \n\t"\
		"movaps		-0x080(%%rbx),%%xmm4			\n\t		movaps		-0x040(%%rbx),%%xmm12\n\t"\
		"movaps		-0x070(%%rax),%%xmm1			\n\t		movaps		-0x030(%%rax),%%xmm9 \n\t"\
		"movaps		-0x070(%%rbx),%%xmm5			\n\t		movaps		-0x030(%%rbx),%%xmm13\n\t"\
		"movaps			  (%%rax),%%xmm2			\n\t		movaps		 0x040(%%rax),%%xmm10\n\t"\
		"movaps		 0x010(%%rbx),%%xmm7			\n\t		movaps		 0x050(%%rbx),%%xmm15\n\t"\
		"movaps		 0x010(%%rax),%%xmm3			\n\t		movaps		 0x050(%%rax),%%xmm11\n\t"\
		"movaps			  (%%rbx),%%xmm6			\n\t		movaps		 0x040(%%rbx),%%xmm14\n\t"\
		"subpd		%%xmm2,%%xmm0				\n\t		subpd		%%xmm10,%%xmm8 \n\t"\
		"subpd		%%xmm7,%%xmm4				\n\t		subpd		%%xmm15,%%xmm12\n\t"\
		"subpd		%%xmm3,%%xmm1				\n\t		subpd		%%xmm11,%%xmm9 \n\t"\
		"subpd		%%xmm6,%%xmm5				\n\t		subpd		%%xmm14,%%xmm13\n\t"\
		"addpd		%%xmm2,%%xmm2				\n\t		addpd		%%xmm10,%%xmm10\n\t"\
		"addpd		%%xmm7,%%xmm7				\n\t		addpd		%%xmm15,%%xmm15\n\t"\
		"addpd		%%xmm3,%%xmm3				\n\t		addpd		%%xmm11,%%xmm11\n\t"\
		"addpd		%%xmm6,%%xmm6				\n\t		addpd		%%xmm14,%%xmm14\n\t"\
		"addpd		%%xmm0,%%xmm2				\n\t		addpd		%%xmm8 ,%%xmm10\n\t"\
		"addpd		%%xmm4,%%xmm7				\n\t		addpd		%%xmm12,%%xmm15\n\t"\
		"addpd		%%xmm1,%%xmm3				\n\t		addpd		%%xmm9 ,%%xmm11\n\t"\
		"addpd		%%xmm5,%%xmm6				\n\t		addpd		%%xmm13,%%xmm14\n\t"\
		"movaps		%%xmm0,      (%%rax)		\n\t		movaps		%%xmm8 , 0x040(%%rax)\n\t"\
		"movaps		%%xmm4,      (%%rbx)		\n\t		movaps		%%xmm12, 0x040(%%rbx)\n\t"\
		"movaps		%%xmm1, 0x010(%%rax)		\n\t		movaps		%%xmm9 , 0x050(%%rax)\n\t"\
		"movaps		%%xmm5,-0x070(%%rbx)		\n\t		movaps		%%xmm13,-0x030(%%rbx)\n\t"\
		"movaps		%%xmm2,-0x080(%%rax)		\n\t		movaps		%%xmm10,-0x040(%%rax)\n\t"\
		"movaps		%%xmm7,-0x080(%%rbx)		\n\t		movaps		%%xmm15,-0x040(%%rbx)\n\t"\
		"movaps		%%xmm3,-0x070(%%rax)		\n\t		movaps		%%xmm11,-0x030(%%rax)\n\t"\
		"movaps		%%xmm6, 0x010(%%rbx)		\n\t		movaps		%%xmm14, 0x050(%%rbx)\n\t"\
		"movaps		-0x080(%%rcx),%%xmm0			\n\t		movaps		-0x040(%%rcx),%%xmm8 \n\t"\
		"movaps		-0x080(%%rdx),%%xmm4			\n\t		movaps		-0x040(%%rdx),%%xmm12\n\t"\
		"movaps		-0x070(%%rcx),%%xmm1			\n\t		movaps		-0x030(%%rcx),%%xmm9 \n\t"\
		"movaps		-0x070(%%rdx),%%xmm5			\n\t		movaps		-0x030(%%rdx),%%xmm13\n\t"\
		"movaps			  (%%rcx),%%xmm2			\n\t		movaps		 0x040(%%rcx),%%xmm10\n\t"\
		"movaps		 0x010(%%rdx),%%xmm7			\n\t		movaps		 0x050(%%rdx),%%xmm15\n\t"\
		"movaps		 0x010(%%rcx),%%xmm3			\n\t		movaps		 0x050(%%rcx),%%xmm11\n\t"\
		"movaps			  (%%rdx),%%xmm6			\n\t		movaps		 0x040(%%rdx),%%xmm14\n\t"\
		"subpd		%%xmm2,%%xmm0				\n\t		subpd		%%xmm10,%%xmm8 \n\t"\
		"subpd		%%xmm7,%%xmm4				\n\t		subpd		%%xmm15,%%xmm12\n\t"\
		"subpd		%%xmm3,%%xmm1				\n\t		subpd		%%xmm11,%%xmm9 \n\t"\
		"subpd		%%xmm6,%%xmm5				\n\t		subpd		%%xmm14,%%xmm13\n\t"\
		"addpd		%%xmm2,%%xmm2				\n\t		addpd		%%xmm10,%%xmm10\n\t"\
		"addpd		%%xmm7,%%xmm7				\n\t		addpd		%%xmm15,%%xmm15\n\t"\
		"addpd		%%xmm3,%%xmm3				\n\t		addpd		%%xmm11,%%xmm11\n\t"\
		"addpd		%%xmm6,%%xmm6				\n\t		addpd		%%xmm14,%%xmm14\n\t"\
		"addpd		%%xmm0,%%xmm2				\n\t		addpd		%%xmm8 ,%%xmm10\n\t"\
		"addpd		%%xmm4,%%xmm7				\n\t		addpd		%%xmm12,%%xmm15\n\t"\
		"addpd		%%xmm1,%%xmm3				\n\t		addpd		%%xmm9 ,%%xmm11\n\t"\
		"addpd		%%xmm5,%%xmm6				\n\t		addpd		%%xmm13,%%xmm14\n\t"\
		"movaps		%%xmm0,      (%%rcx)		\n\t		movaps		%%xmm8 , 0x040(%%rcx)\n\t"\
		"movaps		%%xmm4,      (%%rdx)		\n\t		movaps		%%xmm12, 0x040(%%rdx)\n\t"\
		"movaps		%%xmm1, 0x010(%%rcx)		\n\t		movaps		%%xmm9 , 0x050(%%rcx)\n\t"\
		"movaps		%%xmm5,-0x070(%%rdx)		\n\t		movaps		%%xmm13,-0x030(%%rdx)\n\t"\
		"movaps		%%xmm2,-0x080(%%rcx)		\n\t		movaps		%%xmm10,-0x040(%%rcx)\n\t"\
		"movaps		%%xmm7,-0x080(%%rdx)		\n\t		movaps		%%xmm15,-0x040(%%rdx)\n\t"\
		"movaps		%%xmm3,-0x070(%%rcx)		\n\t		movaps		%%xmm11,-0x030(%%rcx)\n\t"\
		"movaps		%%xmm6, 0x010(%%rdx)		\n\t		movaps		%%xmm14, 0x050(%%rdx)\n\t"\
	/***************************************************************************************/\
	/* Now do eight more radix-4 transforms, including the internal and external twiddles: */\
	/***************************************************************************************/\
		/***********************************************/			/************************************************/\
		/*...Block 2: t02,t12,t22,t32 -> r10,14,12,16: */			/*...Block 6: t0A,t1A,t2A,t3A -> r18,1C,1A,1E:  */\
		/***********************************************/			/************************************************/\
		"movq		%[__isrt2],%%rsi	\n\t"\
		"movq		%[__r10],%%rax	/* base-addr in rcol = c05/r18, so rax/rdi offset +0x80 vs lcol */\n\t"\
		"movq		%[__add1],%%rbx	\n\t"\
		"movq		%%rsi,%%rcx	\n\t"\
		"movq		%%rsi,%%rdx	\n\t"\
		"movq		%[__c01],%%rdi	\n\t"\
		"addq		$0x010,%%rcx	/* cc0 */	\n\t"\
		"addq		$0x030,%%rdx	/* cc1 */	\n\t"\
		"movaps		 0x020(%%rax),%%xmm4				\n\t		movaps		 0x0a0(%%rax),%%xmm12	\n\t"\
		"movaps		 0x060(%%rax),%%xmm0				\n\t		movaps		 0x0e0(%%rax),%%xmm8 	\n\t"\
		"movaps		 0x030(%%rax),%%xmm5				\n\t		movaps		 0x0b0(%%rax),%%xmm13	\n\t"\
		"movaps		 0x070(%%rax),%%xmm1				\n\t		movaps		 0x0f0(%%rax),%%xmm9 	\n\t"\
		"movaps		 0x020(%%rax),%%xmm6				\n\t		movaps		 0x0a0(%%rax),%%xmm14	\n\t"\
		"movaps		 0x060(%%rax),%%xmm2				\n\t		movaps		 0x0e0(%%rax),%%xmm10	\n\t"\
		"movaps		 0x030(%%rax),%%xmm7				\n\t		movaps		 0x0b0(%%rax),%%xmm15	\n\t"\
		"movaps		 0x070(%%rax),%%xmm3				\n\t		movaps		 0x0f0(%%rax),%%xmm11	\n\t"\
		"mulpd			  (%%rdx),%%xmm4				\n\t		mulpd		 0x30 (%%rdx),%%xmm12	\n\t"\
		"mulpd		 0x20 (%%rdx),%%xmm0				\n\t		mulpd			  (%%rdx),%%xmm8 	\n\t"\
		"mulpd			  (%%rdx),%%xmm5				\n\t		mulpd		 0x30 (%%rdx),%%xmm13	\n\t"\
		"mulpd		 0x20 (%%rdx),%%xmm1				\n\t		mulpd			  (%%rdx),%%xmm9 	\n\t"\
		"mulpd		 0x010(%%rdx),%%xmm6				\n\t		mulpd		 0x20 (%%rdx),%%xmm14	\n\t"\
		"mulpd		 0x30 (%%rdx),%%xmm2				\n\t		mulpd		 0x010(%%rdx),%%xmm10	\n\t"\
		"mulpd		 0x010(%%rdx),%%xmm7				\n\t		mulpd		 0x20 (%%rdx),%%xmm15	\n\t"\
		"mulpd		 0x30 (%%rdx),%%xmm3				\n\t		mulpd		 0x010(%%rdx),%%xmm11	\n\t"\
		"subpd		%%xmm6,%%xmm5					\n\t		subpd		%%xmm14,%%xmm13		\n\t"\
		"subpd		%%xmm2,%%xmm1					\n\t		addpd		%%xmm10,%%xmm9 		\n\t"\
		"addpd		%%xmm7,%%xmm4					\n\t		addpd		%%xmm15,%%xmm12		\n\t"\
		"addpd		%%xmm3,%%xmm0					\n\t		subpd		%%xmm11,%%xmm8 		\n\t"\
		"movaps		%%xmm5,%%xmm7					\n\t		movaps		%%xmm13,%%xmm15		\n\t"\
		"movaps		%%xmm4,%%xmm6					\n\t		movaps		%%xmm12,%%xmm14		\n\t"\
		"subpd		%%xmm0,%%xmm4					\n\t		addpd		%%xmm8 ,%%xmm12		\n\t"\
		"subpd		%%xmm1,%%xmm5					\n\t		addpd		%%xmm9 ,%%xmm13		\n\t"\
		"addpd		%%xmm0,%%xmm6					\n\t		subpd		%%xmm8 ,%%xmm14		\n\t"\
		"addpd		%%xmm1,%%xmm7					\n\t		subpd		%%xmm9 ,%%xmm15		\n\t"\
		"movaps		 0x040(%%rax),%%xmm2				\n\t		movaps		 0x0c0(%%rax),%%xmm10	\n\t"\
		"movaps		 0x050(%%rax),%%xmm3				\n\t		movaps		 0x0d0(%%rax),%%xmm11	\n\t"\
		"movaps		 0x040(%%rax),%%xmm0				\n\t		movaps		 0x0c0(%%rax),%%xmm8 	\n\t"\
		"movaps		 0x050(%%rax),%%xmm1				\n\t		movaps		 0x0d0(%%rax),%%xmm9 	\n\t"\
		"mulpd			  (%%rcx),%%xmm2				\n\t		mulpd		 0x010(%%rcx),%%xmm10	\n\t"\
		"mulpd		 0x010(%%rcx),%%xmm1				\n\t		mulpd			  (%%rcx),%%xmm9 	\n\t"\
		"mulpd			  (%%rcx),%%xmm3				\n\t		mulpd		 0x010(%%rcx),%%xmm11	\n\t"\
		"mulpd		 0x010(%%rcx),%%xmm0				\n\t		mulpd			  (%%rcx),%%xmm8 	\n\t"\
		"addpd		%%xmm1,%%xmm2					\n\t		subpd		%%xmm9 ,%%xmm10		\n\t"\
		"subpd		%%xmm0,%%xmm3					\n\t		addpd		%%xmm8 ,%%xmm11		\n\t"\
		"movaps			  (%%rax),%%xmm0				\n\t		movaps		 0x080(%%rax),%%xmm8 	\n\t"\
		"movaps		 0x010(%%rax),%%xmm1				\n\t		movaps		 0x090(%%rax),%%xmm9 	\n\t"\
		"subpd		%%xmm2,%%xmm0					\n\t		subpd		%%xmm10,%%xmm8 		\n\t"\
		"subpd		%%xmm3,%%xmm1					\n\t		subpd		%%xmm11,%%xmm9 		\n\t"\
		"addpd		%%xmm2,%%xmm2					\n\t		addpd		%%xmm10,%%xmm10		\n\t"\
		"addpd		%%xmm3,%%xmm3					\n\t		addpd		%%xmm11,%%xmm11		\n\t"\
		"addpd		%%xmm0,%%xmm2					\n\t		addpd		%%xmm8 ,%%xmm10		\n\t"\
		"addpd		%%xmm1,%%xmm3					\n\t		addpd		%%xmm9 ,%%xmm11		\n\t"\
		"subpd		%%xmm6,%%xmm2					\n\t		subpd		%%xmm14,%%xmm8 		\n\t"\
		"subpd		%%xmm7,%%xmm3					\n\t		subpd		%%xmm15,%%xmm9 		\n\t"\
		"addpd		%%xmm6,%%xmm6					\n\t		addpd		%%xmm14,%%xmm14		\n\t"\
		"addpd		%%xmm7,%%xmm7					\n\t		addpd		%%xmm15,%%xmm15		\n\t"\
		"addpd		%%xmm2,%%xmm6					\n\t		addpd		%%xmm8 ,%%xmm14		\n\t"\
		"addpd		%%xmm3,%%xmm7					\n\t		addpd		%%xmm9 ,%%xmm15		\n\t"\
		"movaps		%%xmm2, 0x020(%%rax)			\n\t		movaps		%%xmm8 , 0x0a0(%%rax)	\n\t"\
		"movaps		%%xmm3, 0x030(%%rax)			\n\t		movaps		%%xmm9 , 0x0b0(%%rax)	\n\t"\
		"movaps		%%xmm6,%%xmm2					\n\t		movaps		%%xmm14,%%xmm8 		\n\t"\
		"movaps		%%xmm7,%%xmm3					\n\t		movaps		%%xmm15,%%xmm9 		\n\t"\
		"mulpd			  (%%rdi),%%xmm6				\n\t		mulpd		 0x080(%%rdi),%%xmm14	\n\t"\
		"mulpd			  (%%rdi),%%xmm7				\n\t		mulpd		 0x080(%%rdi),%%xmm15	\n\t"\
		"mulpd		 0x010(%%rdi),%%xmm2				\n\t		mulpd		 0x090(%%rdi),%%xmm8 	\n\t"\
		"mulpd		 0x010(%%rdi),%%xmm3				\n\t		mulpd		 0x090(%%rdi),%%xmm9 	\n\t"\
		"subpd		%%xmm2,%%xmm7					\n\t		subpd		%%xmm8 ,%%xmm15		\n\t"\
		"addpd		%%xmm3,%%xmm6					\n\t		addpd		%%xmm9 ,%%xmm14		\n\t"\
		"movaps		%%xmm7, 0x010(%%rbx)			\n\t		movaps		%%xmm15, 0x050(%%rbx)	\n\t"\
		"movaps		%%xmm6,      (%%rbx)			\n\t		movaps		%%xmm14, 0x040(%%rbx)	\n\t"\
		"movaps		 0x020(%%rax),%%xmm6				\n\t		movaps		 0x0a0(%%rax),%%xmm14	\n\t"\
		"movaps		 0x030(%%rax),%%xmm7				\n\t		movaps		 0x0b0(%%rax),%%xmm15	\n\t"\
		"movaps		%%xmm6,%%xmm2					\n\t		movaps		%%xmm14,%%xmm8 		\n\t"\
		"movaps		%%xmm7,%%xmm3					\n\t		movaps		%%xmm15,%%xmm9 		\n\t"\
		"mulpd		 0x20 (%%rdi),%%xmm6				\n\t		mulpd		 0xa0 (%%rdi),%%xmm14	\n\t"\
		"mulpd		 0x20 (%%rdi),%%xmm7				\n\t		mulpd		 0xa0 (%%rdi),%%xmm15	\n\t"\
		"mulpd		 0x30 (%%rdi),%%xmm2				\n\t		mulpd		 0xb0 (%%rdi),%%xmm8 	\n\t"\
		"mulpd		 0x30 (%%rdi),%%xmm3				\n\t		mulpd		 0xb0 (%%rdi),%%xmm9 	\n\t"\
		"subpd		%%xmm2,%%xmm7					\n\t		subpd		%%xmm8 ,%%xmm15		\n\t"\
		"addpd		%%xmm3,%%xmm6					\n\t		addpd		%%xmm9 ,%%xmm14		\n\t"\
		"movaps		%%xmm7, 0x110(%%rbx)			\n\t		movaps		%%xmm15, 0x150(%%rbx)	\n\t"\
		"movaps		%%xmm6, 0x100(%%rbx)			\n\t		movaps		%%xmm14, 0x140(%%rbx)	\n\t"\
		"addq		$0x040,%%rdi					\n\t"\
		"subpd		%%xmm5,%%xmm0					\n\t		subpd		%%xmm13,%%xmm10		\n\t"\
		"subpd		%%xmm4,%%xmm1					\n\t		subpd		%%xmm12,%%xmm11		\n\t"\
		"addpd		%%xmm5,%%xmm5					\n\t		addpd		%%xmm13,%%xmm13		\n\t"\
		"addpd		%%xmm4,%%xmm4					\n\t		addpd		%%xmm12,%%xmm12		\n\t"\
		"addpd		%%xmm0,%%xmm5					\n\t		addpd		%%xmm10,%%xmm13		\n\t"\
		"addpd		%%xmm1,%%xmm4					\n\t		addpd		%%xmm11,%%xmm12		\n\t"\
		"movaps		%%xmm5,%%xmm2					\n\t		movaps		%%xmm13,%%xmm8 		\n\t"\
		"movaps		%%xmm1,%%xmm3					\n\t		movaps		%%xmm11,%%xmm9 		\n\t"\
		"mulpd			  (%%rdi),%%xmm5				\n\t		mulpd		 0x080(%%rdi),%%xmm13	\n\t"\
		"mulpd			  (%%rdi),%%xmm1				\n\t		mulpd		 0x080(%%rdi),%%xmm11	\n\t"\
		"mulpd		 0x010(%%rdi),%%xmm2				\n\t		mulpd		 0x090(%%rdi),%%xmm8 	\n\t"\
		"mulpd		 0x010(%%rdi),%%xmm3				\n\t		mulpd		 0x090(%%rdi),%%xmm9 	\n\t"\
		"subpd		%%xmm2,%%xmm1					\n\t		subpd		%%xmm8 ,%%xmm11		\n\t"\
		"addpd		%%xmm3,%%xmm5					\n\t		addpd		%%xmm9 ,%%xmm13		\n\t"\
		"movaps		%%xmm1, 0x090(%%rbx)			\n\t		movaps		%%xmm11, 0x0d0(%%rbx)	\n\t"\
		"movaps		%%xmm5, 0x080(%%rbx)			\n\t		movaps		%%xmm13, 0x0c0(%%rbx)	\n\t"\
		"movaps		%%xmm0,%%xmm2					\n\t		movaps		%%xmm10,%%xmm8 		\n\t"\
		"movaps		%%xmm4,%%xmm3					\n\t		movaps		%%xmm12,%%xmm9 		\n\t"\
		"mulpd		 0x20 (%%rdi),%%xmm0				\n\t		mulpd		 0xa0 (%%rdi),%%xmm10	\n\t"\
		"mulpd		 0x20 (%%rdi),%%xmm4				\n\t		mulpd		 0xa0 (%%rdi),%%xmm12	\n\t"\
		"mulpd		 0x30 (%%rdi),%%xmm2				\n\t		mulpd		 0xb0 (%%rdi),%%xmm8 	\n\t"\
		"mulpd		 0x30 (%%rdi),%%xmm3				\n\t		mulpd		 0xb0 (%%rdi),%%xmm9 	\n\t"\
		"subpd		%%xmm2,%%xmm4					\n\t		subpd		%%xmm8 ,%%xmm12		\n\t"\
		"addpd		%%xmm3,%%xmm0					\n\t		addpd		%%xmm9 ,%%xmm10		\n\t"\
		"movaps		%%xmm4, 0x190(%%rbx)			\n\t		movaps		%%xmm12, 0x1d0(%%rbx)	\n\t"\
		"movaps		%%xmm0, 0x180(%%rbx)			\n\t		movaps		%%xmm10, 0x1c0(%%rbx)	\n\t"\
		"\n\t"\
		"/************************************************/	\n\t	/************************************************/	\n\t"\
		"/*...Block 4: t06,t16,t26,t36 -> r30,34,32,36:  */	\n\t	/*...Block 8: t0E,t1E,t2E,t3E -> r38,3C,3A,3E:  */	\n\t"\
		"/************************************************/	\n\t	/************************************************/	\n\t"\
		"addq		$0x200,%%rax					\n\t		addq		$0x080,%%rdi			\n\t"\
		"movaps		 0x020(%%rax),%%xmm4				\n\t		movaps		 0x0a0(%%rax),%%xmm12	\n\t"\
		"movaps		 0x060(%%rax),%%xmm0				\n\t		movaps		 0x0e0(%%rax),%%xmm8 	\n\t"\
		"movaps		 0x030(%%rax),%%xmm5				\n\t		movaps		 0x0b0(%%rax),%%xmm13	\n\t"\
		"movaps		 0x070(%%rax),%%xmm1				\n\t		movaps		 0x0f0(%%rax),%%xmm9 	\n\t"\
		"movaps		 0x020(%%rax),%%xmm6				\n\t		movaps		 0x0a0(%%rax),%%xmm14	\n\t"\
		"movaps		 0x060(%%rax),%%xmm2				\n\t		movaps		 0x0e0(%%rax),%%xmm10	\n\t"\
		"movaps		 0x030(%%rax),%%xmm7				\n\t		movaps		 0x0b0(%%rax),%%xmm15	\n\t"\
		"movaps		 0x070(%%rax),%%xmm3				\n\t		movaps		 0x0f0(%%rax),%%xmm11	\n\t"\
		"mulpd		 0x20 (%%rdx),%%xmm4				\n\t		mulpd		 0x010(%%rdx),%%xmm12	\n\t"\
		"mulpd		 0x010(%%rdx),%%xmm0				\n\t		mulpd		 0x30 (%%rdx),%%xmm8 	\n\t"\
		"mulpd		 0x20 (%%rdx),%%xmm5				\n\t		mulpd		 0x010(%%rdx),%%xmm13	\n\t"\
		"mulpd		 0x010(%%rdx),%%xmm1				\n\t		mulpd		 0x30 (%%rdx),%%xmm9 	\n\t"\
		"mulpd		 0x30 (%%rdx),%%xmm6				\n\t		mulpd			  (%%rdx),%%xmm14	\n\t"\
		"mulpd			  (%%rdx),%%xmm2				\n\t		mulpd		 0x20 (%%rdx),%%xmm10	\n\t"\
		"mulpd		 0x30 (%%rdx),%%xmm7				\n\t		mulpd			  (%%rdx),%%xmm15	\n\t"\
		"mulpd			  (%%rdx),%%xmm3				\n\t		mulpd		 0x20 (%%rdx),%%xmm11	\n\t"\
		"subpd		%%xmm6,%%xmm5					\n\t		subpd		%%xmm14,%%xmm13		\n\t"\
		"addpd		%%xmm2,%%xmm1					\n\t		subpd		%%xmm10,%%xmm9 		\n\t"\
		"addpd		%%xmm7,%%xmm4					\n\t		addpd		%%xmm15,%%xmm12		\n\t"\
		"subpd		%%xmm3,%%xmm0					\n\t		addpd		%%xmm11,%%xmm8 		\n\t"\
		"movaps		%%xmm5,%%xmm7					\n\t		movaps		%%xmm13,%%xmm15		\n\t"\
		"movaps		%%xmm4,%%xmm6					\n\t		movaps		%%xmm12,%%xmm14		\n\t"\
		"addpd		%%xmm0,%%xmm4					\n\t		addpd		%%xmm8 ,%%xmm12		\n\t"\
		"addpd		%%xmm1,%%xmm5					\n\t		addpd		%%xmm9 ,%%xmm13		\n\t"\
		"subpd		%%xmm0,%%xmm6					\n\t		subpd		%%xmm8 ,%%xmm14		\n\t"\
		"subpd		%%xmm1,%%xmm7					\n\t		subpd		%%xmm9 ,%%xmm15		\n\t"\
		"movaps		 0x040(%%rax),%%xmm2				\n\t		movaps		 0x0c0(%%rax),%%xmm10	\n\t"\
		"movaps		 0x050(%%rax),%%xmm3				\n\t		movaps		 0x0d0(%%rax),%%xmm11	\n\t"\
		"movaps		 0x040(%%rax),%%xmm0				\n\t		movaps		 0x0c0(%%rax),%%xmm8 	\n\t"\
		"movaps		 0x050(%%rax),%%xmm1				\n\t		movaps		 0x0d0(%%rax),%%xmm9 	\n\t"\
		"mulpd		 0x010(%%rcx),%%xmm2				\n\t		mulpd			  (%%rcx),%%xmm10	\n\t"\
		"mulpd			  (%%rcx),%%xmm1				\n\t		mulpd		 0x010(%%rcx),%%xmm9 	\n\t"\
		"mulpd		 0x010(%%rcx),%%xmm3				\n\t		mulpd			  (%%rcx),%%xmm11	\n\t"\
		"mulpd			  (%%rcx),%%xmm0				\n\t		mulpd		 0x010(%%rcx),%%xmm8 	\n\t"\
		"addpd		%%xmm1,%%xmm2					\n\t		subpd		%%xmm9 ,%%xmm10		\n\t"\
		"subpd		%%xmm0,%%xmm3					\n\t		addpd		%%xmm8 ,%%xmm11		\n\t"\
		"movaps			  (%%rax),%%xmm0				\n\t		movaps		 0x080(%%rax),%%xmm8 	\n\t"\
		"movaps		 0x010(%%rax),%%xmm1				\n\t		movaps		 0x090(%%rax),%%xmm9 	\n\t"\
		"subpd		%%xmm2,%%xmm0					\n\t		subpd		%%xmm10,%%xmm8 		\n\t"\
		"subpd		%%xmm3,%%xmm1					\n\t		subpd		%%xmm11,%%xmm9 		\n\t"\
		"addpd		%%xmm2,%%xmm2					\n\t		addpd		%%xmm10,%%xmm10		\n\t"\
		"addpd		%%xmm3,%%xmm3					\n\t		addpd		%%xmm11,%%xmm11		\n\t"\
		"addpd		%%xmm0,%%xmm2					\n\t		addpd		%%xmm8 ,%%xmm10		\n\t"\
		"addpd		%%xmm1,%%xmm3					\n\t		addpd		%%xmm9 ,%%xmm11		\n\t"\
		"addq		$0x040,%%rdi					\n\t"\
		"subpd		%%xmm6,%%xmm2					\n\t		subpd		%%xmm14,%%xmm8 		\n\t"\
		"subpd		%%xmm7,%%xmm3					\n\t		subpd		%%xmm15,%%xmm9 		\n\t"\
		"addpd		%%xmm6,%%xmm6					\n\t		addpd		%%xmm14,%%xmm14		\n\t"\
		"addpd		%%xmm7,%%xmm7					\n\t		addpd		%%xmm15,%%xmm15		\n\t"\
		"addpd		%%xmm2,%%xmm6					\n\t		addpd		%%xmm8 ,%%xmm14		\n\t"\
		"addpd		%%xmm3,%%xmm7					\n\t		addpd		%%xmm9 ,%%xmm15		\n\t"\
		"movaps		%%xmm2, 0x020(%%rax)			\n\t		movaps		%%xmm8 , 0x0a0(%%rax)	\n\t"\
		"movaps		%%xmm3, 0x030(%%rax)			\n\t		movaps		%%xmm9 , 0x0b0(%%rax)	\n\t"\
		"movaps		%%xmm6,%%xmm2					\n\t		movaps		%%xmm14,%%xmm8 		\n\t"\
		"movaps		%%xmm7,%%xmm3					\n\t		movaps		%%xmm15,%%xmm9 		\n\t"\
		"mulpd			  (%%rdi),%%xmm6				\n\t		mulpd		 0x080(%%rdi),%%xmm14	\n\t"\
		"mulpd			  (%%rdi),%%xmm7				\n\t		mulpd		 0x080(%%rdi),%%xmm15	\n\t"\
		"mulpd		 0x010(%%rdi),%%xmm2				\n\t		mulpd		 0x090(%%rdi),%%xmm8 	\n\t"\
		"mulpd		 0x010(%%rdi),%%xmm3				\n\t		mulpd		 0x090(%%rdi),%%xmm9 	\n\t"\
		"subpd		%%xmm2,%%xmm7					\n\t		subpd		%%xmm8 ,%%xmm15		\n\t"\
		"addpd		%%xmm3,%%xmm6					\n\t		addpd		%%xmm9 ,%%xmm14		\n\t"\
		"movaps		%%xmm7, 0x030(%%rbx)			\n\t		movaps		%%xmm15, 0x070(%%rbx)	\n\t"\
		"movaps		%%xmm6, 0x020(%%rbx)			\n\t		movaps		%%xmm14, 0x060(%%rbx)	\n\t"\
		"movaps		 0x020(%%rax),%%xmm6				\n\t		movaps		 0x0a0(%%rax),%%xmm14	\n\t"\
		"movaps		 0x030(%%rax),%%xmm7				\n\t		movaps		 0x0b0(%%rax),%%xmm15	\n\t"\
		"movaps		%%xmm6,%%xmm2					\n\t		movaps		%%xmm14,%%xmm8 		\n\t"\
		"movaps		%%xmm7,%%xmm3					\n\t		movaps		%%xmm15,%%xmm9 		\n\t"\
		"mulpd		 0x20 (%%rdi),%%xmm6				\n\t		mulpd		 0xa0 (%%rdi),%%xmm14	\n\t"\
		"mulpd		 0x20 (%%rdi),%%xmm7				\n\t		mulpd		 0xa0 (%%rdi),%%xmm15	\n\t"\
		"mulpd		 0x30 (%%rdi),%%xmm2				\n\t		mulpd		 0xb0 (%%rdi),%%xmm8 	\n\t"\
		"mulpd		 0x30 (%%rdi),%%xmm3				\n\t		mulpd		 0xb0 (%%rdi),%%xmm9 	\n\t"\
		"subpd		%%xmm2,%%xmm7					\n\t		subpd		%%xmm8 ,%%xmm15		\n\t"\
		"addpd		%%xmm3,%%xmm6					\n\t		addpd		%%xmm9 ,%%xmm14		\n\t"\
		"movaps		%%xmm7, 0x130(%%rbx)			\n\t		movaps		%%xmm15, 0x170(%%rbx)	\n\t"\
		"movaps		%%xmm6, 0x120(%%rbx)			\n\t		movaps		%%xmm14, 0x160(%%rbx)	\n\t"\
		"addq		$0x040,%%rdi					\n\t"\
		"subpd		%%xmm5,%%xmm0					\n\t		subpd		%%xmm13,%%xmm10		\n\t"\
		"subpd		%%xmm4,%%xmm1					\n\t		subpd		%%xmm12,%%xmm11		\n\t"\
		"addpd		%%xmm5,%%xmm5					\n\t		addpd		%%xmm13,%%xmm13		\n\t"\
		"addpd		%%xmm4,%%xmm4					\n\t		addpd		%%xmm12,%%xmm12		\n\t"\
		"addpd		%%xmm0,%%xmm5					\n\t		addpd		%%xmm10,%%xmm13		\n\t"\
		"addpd		%%xmm1,%%xmm4					\n\t		addpd		%%xmm11,%%xmm12		\n\t"\
		"movaps		%%xmm5,%%xmm2					\n\t		movaps		%%xmm13,%%xmm8 		\n\t"\
		"movaps		%%xmm1,%%xmm3					\n\t		movaps		%%xmm11,%%xmm9 		\n\t"\
		"mulpd			  (%%rdi),%%xmm5				\n\t		mulpd		 0x080(%%rdi),%%xmm13	\n\t"\
		"mulpd			  (%%rdi),%%xmm1				\n\t		mulpd		 0x080(%%rdi),%%xmm11	\n\t"\
		"mulpd		 0x010(%%rdi),%%xmm2				\n\t		mulpd		 0x090(%%rdi),%%xmm8 	\n\t"\
		"mulpd		 0x010(%%rdi),%%xmm3				\n\t		mulpd		 0x090(%%rdi),%%xmm9 	\n\t"\
		"subpd		%%xmm2,%%xmm1					\n\t		subpd		%%xmm8 ,%%xmm11		\n\t"\
		"addpd		%%xmm3,%%xmm5					\n\t		addpd		%%xmm9 ,%%xmm13		\n\t"\
		"movaps		%%xmm1, 0x0b0(%%rbx)			\n\t		movaps		%%xmm11, 0x0f0(%%rbx)	\n\t"\
		"movaps		%%xmm5, 0x0a0(%%rbx)			\n\t		movaps		%%xmm13, 0x0e0(%%rbx)	\n\t"\
		"movaps		%%xmm0,%%xmm2					\n\t		movaps		%%xmm10,%%xmm8 		\n\t"\
		"movaps		%%xmm4,%%xmm3					\n\t		movaps		%%xmm12,%%xmm9 		\n\t"\
		"mulpd		 0x20 (%%rdi),%%xmm0				\n\t		mulpd		 0xa0 (%%rdi),%%xmm10	\n\t"\
		"mulpd		 0x20 (%%rdi),%%xmm4				\n\t		mulpd		 0xa0 (%%rdi),%%xmm12	\n\t"\
		"mulpd		 0x30 (%%rdi),%%xmm2				\n\t		mulpd		 0xb0 (%%rdi),%%xmm8 	\n\t"\
		"mulpd		 0x30 (%%rdi),%%xmm3				\n\t		mulpd		 0xb0 (%%rdi),%%xmm9 	\n\t"\
		"subpd		%%xmm2,%%xmm4					\n\t		subpd		%%xmm8 ,%%xmm12		\n\t"\
		"addpd		%%xmm3,%%xmm0					\n\t		addpd		%%xmm9 ,%%xmm10		\n\t"\
		"movaps		%%xmm4, 0x1b0(%%rbx)			\n\t		movaps		%%xmm12, 0x1f0(%%rbx)	\n\t"\
		"movaps		%%xmm0, 0x1a0(%%rbx)			\n\t		movaps		%%xmm10, 0x1e0(%%rbx)	\n\t"\
		"/************************************************/	\n\t		/************************************************/	\n\t"\
		"/*...Block 1: t00,t10,t20,t30 -> r00,04,02,06:  */	\n\t		/*...Block 5: t08,t18,t28,t38 -> r08,0C,0A,0E:  */	\n\t"\
		"/************************************************/	\n\t		/************************************************/	\n\t"\
		"movq	%[__r00],%%rdx	/* base-addr in rcol = c04/r08, so rcx,rdx offset +0x80 in rcol */	\n\t	movaps	(%%rsi),%%xmm10	/* isrt2 */	\n\t"\
		"movaps			  (%%rdx),%%xmm0					\n\t		movaps		 0x0a0(%%rdx),%%xmm12				\n\t"\
		"movaps		 0x010(%%rdx),%%xmm1					\n\t		movaps		 0x0b0(%%rdx),%%xmm13				\n\t"\
		"movaps		 0x040(%%rdx),%%xmm2					\n\t		movaps		 0x0e0(%%rdx),%%xmm8 				\n\t"\
		"movaps		 0x050(%%rdx),%%xmm3					\n\t		movaps		 0x0f0(%%rdx),%%xmm9 				\n\t"\
		"subpd		 0x040(%%rdx),%%xmm0					\n\t		addpd		 0x0b0(%%rdx),%%xmm12				\n\t"\
		"subpd		 0x050(%%rdx),%%xmm1					\n\t		subpd		 0x0a0(%%rdx),%%xmm13				\n\t"\
		"addpd			  (%%rdx),%%xmm2					\n\t		subpd		 0x0f0(%%rdx),%%xmm8 				\n\t"\
		"addpd		 0x010(%%rdx),%%xmm3					\n\t		addpd		 0x0e0(%%rdx),%%xmm9 				\n\t"\
		"movaps		 0x020(%%rdx),%%xmm4					\n\t		mulpd		%%xmm10,%%xmm12					\n\t"\
		"movaps		 0x030(%%rdx),%%xmm5					\n\t		mulpd		%%xmm10,%%xmm13					\n\t"\
		"movaps		 0x060(%%rdx),%%xmm6					\n\t		mulpd		%%xmm10,%%xmm8 					\n\t"\
		"movaps		 0x070(%%rdx),%%xmm7					\n\t		mulpd		%%xmm10,%%xmm9 					\n\t"\
		"subpd		 0x060(%%rdx),%%xmm4					\n\t		movaps		%%xmm12,%%xmm14					\n\t"\
		"subpd		 0x070(%%rdx),%%xmm5					\n\t		movaps		%%xmm13,%%xmm15					\n\t"\
		"addpd		 0x020(%%rdx),%%xmm6					\n\t		subpd		%%xmm8 ,%%xmm12					\n\t"\
		"addpd		 0x030(%%rdx),%%xmm7					\n\t		subpd		%%xmm9 ,%%xmm13					\n\t"\
		"movq		%[__add0],%%rax					\n\t		addpd		%%xmm8 ,%%xmm14					\n\t"\
		"movq		%[__c10],%%rcx					\n\t		addpd		%%xmm9 ,%%xmm15					\n\t"\
		"addpd		%%xmm6,%%xmm2						\n\t		movaps		 0x080(%%rdx),%%xmm8 				\n\t"\
		"addpd		%%xmm7,%%xmm3						\n\t		movaps		 0x090(%%rdx),%%xmm9 				\n\t"\
		"movaps		%%xmm2,      (%%rdx)				\n\t		movaps		 0x0c0(%%rdx),%%xmm10				\n\t"\
		"movaps		%%xmm3, 0x010(%%rdx)				\n\t		movaps		 0x0d0(%%rdx),%%xmm11				\n\t"\
		"addpd		%%xmm6,%%xmm6						\n\t		subpd		 0x0d0(%%rdx),%%xmm8 				\n\t"\
		"addpd		%%xmm7,%%xmm7						\n\t		subpd		 0x0c0(%%rdx),%%xmm9 				\n\t"\
		"subpd		%%xmm6,%%xmm2						\n\t		addpd		 0x080(%%rdx),%%xmm11				\n\t"\
		"subpd		%%xmm7,%%xmm3						\n\t		addpd		 0x090(%%rdx),%%xmm10				\n\t"\
		"movaps		%%xmm2,%%xmm6						\n\t		subpd		%%xmm12,%%xmm11					\n\t"\
		"movaps		%%xmm3,%%xmm7						\n\t		subpd		%%xmm13,%%xmm9 					\n\t"\
		"mulpd			  (%%rcx),%%xmm2					\n\t		addpd		%%xmm12,%%xmm12					\n\t"\
		"mulpd			  (%%rcx),%%xmm3					\n\t		addpd		%%xmm13,%%xmm13					\n\t"\
		"mulpd		 0x010(%%rcx),%%xmm6					\n\t		addpd		%%xmm11,%%xmm12					\n\t"\
		"mulpd		 0x010(%%rcx),%%xmm7					\n\t		addpd		%%xmm9 ,%%xmm13					\n\t"\
		"subpd		%%xmm6,%%xmm3						\n\t		movaps		%%xmm11, 0x080(%%rdx)				\n\t"\
		"addpd		%%xmm7,%%xmm2						\n\t		movaps		%%xmm9 , 0x090(%%rdx)				\n\t"\
		"subq	$0x20,%%rcx	/* put c00 in rcx to ease bookkeeping*/\n\t	movaps		%%xmm12,%%xmm11					\n\t"\
		"movq		%[__add1],%%rbx					\n\t		movaps		%%xmm13,%%xmm9 					\n\t"\
		"movaps		%%xmm3,%%xmm7						\n\t		mulpd		 0x080(%%rcx),%%xmm12	/* c04 */	\n\t"\
		"movaps		%%xmm2,%%xmm6						\n\t		mulpd		 0x080(%%rcx),%%xmm13				\n\t"\
		"unpckhpd	 0x110(%%rbx),%%xmm7					\n\t		mulpd		 0x090(%%rcx),%%xmm11				\n\t"\
		"unpcklpd	 0x110(%%rbx),%%xmm3					\n\t		mulpd		 0x090(%%rcx),%%xmm9 				\n\t"\
		"movaps		%%xmm7, 0x110(%%rbx)				\n\t		subpd		%%xmm11,%%xmm13					\n\t"\
		"unpckhpd	 0x100(%%rbx),%%xmm6					\n\t		addpd		%%xmm9 ,%%xmm12					\n\t"\
		"unpcklpd	 0x100(%%rbx),%%xmm2					\n\t		movaps		%%xmm13,%%xmm11					\n\t"\
		"movaps		%%xmm6, 0x100(%%rbx)				\n\t		movaps		%%xmm12,%%xmm9 					\n\t"\
		"movaps		%%xmm3, 0x110(%%rax)				\n\t		unpckhpd	 0x050(%%rbx),%%xmm11				\n\t"\
		"movaps		%%xmm2, 0x100(%%rax)				\n\t		unpcklpd	 0x050(%%rbx),%%xmm13				\n\t"\
		"movaps		 0x010(%%rdx),%%xmm3					\n\t		movaps		%%xmm11, 0x050(%%rbx)				\n\t"\
		"movaps			  (%%rdx),%%xmm2					\n\t		unpckhpd	 0x040(%%rbx),%%xmm9 				\n\t"\
		"movaps		%%xmm3,%%xmm7						\n\t		unpcklpd	 0x040(%%rbx),%%xmm12				\n\t"\
		"movaps		%%xmm2,%%xmm6						\n\t		movaps		%%xmm9 , 0x040(%%rbx)				\n\t"\
		"unpckhpd	 0x010(%%rbx),%%xmm7					\n\t		movaps		%%xmm13, 0x050(%%rax)				\n\t"\
		"unpcklpd	 0x010(%%rbx),%%xmm3					\n\t		movaps		%%xmm12, 0x040(%%rax)				\n\t"\
		"movaps		%%xmm7, 0x010(%%rbx)				\n\t		movaps		 0x080(%%rdx),%%xmm12				\n\t"\
		"unpckhpd		(%%rbx),%%xmm6					\n\t		movaps		 0x090(%%rdx),%%xmm13				\n\t"\
		"unpcklpd		(%%rbx),%%xmm2					\n\t		movaps		%%xmm12,%%xmm11					\n\t"\
		"movaps		%%xmm6,      (%%rbx)				\n\t		movaps		%%xmm13,%%xmm9 					\n\t"\
		"movaps		%%xmm3, 0x010(%%rax)				\n\t		mulpd		 0x0a0(%%rcx),%%xmm12	/* c14 */	\n\t"\
		"movaps		%%xmm2,      (%%rax)				\n\t		mulpd		 0x0a0(%%rcx),%%xmm13				\n\t"\
		"addpd		%%xmm5,%%xmm0						\n\t		mulpd		 0x0b0(%%rcx),%%xmm11				\n\t"\
		"subpd		%%xmm4,%%xmm1						\n\t		mulpd		 0x0b0(%%rcx),%%xmm9 				\n\t"\
		"movaps		%%xmm0,%%xmm2						\n\t		subpd		%%xmm11,%%xmm13					\n\t"\
		"movaps		%%xmm1,%%xmm3						\n\t		addpd		%%xmm9 ,%%xmm12					\n\t"\
		"addpd		%%xmm5,%%xmm5						\n\t		movaps		%%xmm13,%%xmm11					\n\t"\
		"addpd		%%xmm4,%%xmm4						\n\t		movaps		%%xmm12,%%xmm9 					\n\t"\
		"movaps		%%xmm0,%%xmm6						\n\t		unpckhpd	 0x150(%%rbx),%%xmm11				\n\t"\
		"movaps		%%xmm1,%%xmm7						\n\t		unpcklpd	 0x150(%%rbx),%%xmm13				\n\t"\
		"mulpd		 0x040(%%rcx),%%xmm2		/* c08 */	\n\t		movaps		%%xmm11, 0x150(%%rbx)				\n\t"\
		"mulpd		 0x040(%%rcx),%%xmm3					\n\t		unpckhpd	 0x140(%%rbx),%%xmm9 				\n\t"\
		"mulpd		 0x050(%%rcx),%%xmm6					\n\t		unpcklpd	 0x140(%%rbx),%%xmm12				\n\t"\
		"mulpd		 0x050(%%rcx),%%xmm7					\n\t		movaps		%%xmm9 , 0x140(%%rbx)				\n\t"\
		"subpd		%%xmm6,%%xmm3						\n\t		movaps		%%xmm13, 0x150(%%rax)				\n\t"\
		"addpd		%%xmm7,%%xmm2						\n\t		movaps		%%xmm12, 0x140(%%rax)				\n\t"\
		"movaps		%%xmm3,%%xmm7						\n\t		subpd		%%xmm15,%%xmm8 					\n\t"\
		"movaps		%%xmm2,%%xmm6						\n\t		subpd		%%xmm14,%%xmm10					\n\t"\
		"unpckhpd	 0x090(%%rbx),%%xmm7					\n\t		addpd		%%xmm15,%%xmm15					\n\t"\
		"unpcklpd	 0x090(%%rbx),%%xmm3					\n\t		addpd		%%xmm14,%%xmm14					\n\t"\
		"movaps		%%xmm7, 0x090(%%rbx)				\n\t		addpd		%%xmm8 ,%%xmm15					\n\t"\
		"unpckhpd	 0x080(%%rbx),%%xmm6					\n\t		addpd		%%xmm10,%%xmm14					\n\t"\
		"unpcklpd	 0x080(%%rbx),%%xmm2					\n\t		movaps		%%xmm15,%%xmm12					\n\t"\
		"movaps		%%xmm6, 0x080(%%rbx)				\n\t		movaps		%%xmm10,%%xmm13					\n\t"\
		"movaps		%%xmm3, 0x090(%%rax)				\n\t		mulpd		 0x0c0(%%rcx),%%xmm15	/* c0C */	\n\t"\
		"movaps		%%xmm2, 0x080(%%rax)				\n\t		mulpd		 0x0c0(%%rcx),%%xmm10				\n\t"\
		"subpd		%%xmm5,%%xmm0						\n\t		mulpd		 0x0d0(%%rcx),%%xmm12				\n\t"\
		"addpd		%%xmm4,%%xmm1						\n\t		mulpd		 0x0d0(%%rcx),%%xmm13				\n\t"\
		"movaps		%%xmm0,%%xmm6						\n\t		subpd		%%xmm12,%%xmm10					\n\t"\
		"movaps		%%xmm1,%%xmm7						\n\t		addpd		%%xmm13,%%xmm15					\n\t"\
		"mulpd		 0x060(%%rcx),%%xmm0		/* c18 */	\n\t		movaps		%%xmm10,%%xmm13					\n\t"\
		"mulpd		 0x060(%%rcx),%%xmm1					\n\t		movaps		%%xmm15,%%xmm12					\n\t"\
		"mulpd		 0x070(%%rcx),%%xmm6					\n\t		unpckhpd	 0x0d0(%%rbx),%%xmm13				\n\t"\
		"mulpd		 0x070(%%rcx),%%xmm7					\n\t		unpcklpd	 0x0d0(%%rbx),%%xmm10				\n\t"\
		"subpd		%%xmm6,%%xmm1						\n\t		movaps		%%xmm13, 0x0d0(%%rbx)				\n\t"\
		"addpd		%%xmm7,%%xmm0						\n\t		unpckhpd	 0x0c0(%%rbx),%%xmm12				\n\t"\
		"movaps		%%xmm1,%%xmm7						\n\t		unpcklpd	 0x0c0(%%rbx),%%xmm15				\n\t"\
		"movaps		%%xmm0,%%xmm6						\n\t		movaps		%%xmm12, 0x0c0(%%rbx)				\n\t"\
		"unpckhpd	 0x190(%%rbx),%%xmm7					\n\t		movaps		%%xmm10, 0x0d0(%%rax)				\n\t"\
		"unpcklpd	 0x190(%%rbx),%%xmm1					\n\t		movaps		%%xmm15, 0x0c0(%%rax)				\n\t"\
		"movaps		%%xmm7, 0x190(%%rbx)				\n\t		movaps		%%xmm8 ,%%xmm12					\n\t"\
		"unpckhpd	 0x180(%%rbx),%%xmm6					\n\t		movaps		%%xmm14,%%xmm13					\n\t"\
		"unpcklpd	 0x180(%%rbx),%%xmm0					\n\t		mulpd		 0x0e0(%%rcx),%%xmm8 	/* c1C */	\n\t"\
		"movaps		%%xmm6, 0x180(%%rbx)				\n\t		mulpd		 0x0e0(%%rcx),%%xmm14				\n\t"\
		"movaps		%%xmm1, 0x190(%%rax)				\n\t		mulpd		 0x0f0(%%rcx),%%xmm12				\n\t"\
		"movaps		%%xmm0, 0x180(%%rax)				\n\t		mulpd		 0x0f0(%%rcx),%%xmm13				\n\t"\
		"																subpd		%%xmm12,%%xmm14					\n\t"\
		"																addpd		%%xmm13,%%xmm8 					\n\t"\
		"																movaps		%%xmm14,%%xmm13					\n\t"\
		"																movaps		%%xmm8 ,%%xmm12					\n\t"\
		"																unpckhpd	 0x1d0(%%rbx),%%xmm13				\n\t"\
		"																unpcklpd	 0x1d0(%%rbx),%%xmm14				\n\t"\
		"																movaps		%%xmm13, 0x1d0(%%rbx)				\n\t"\
		"																unpckhpd	 0x1c0(%%rbx),%%xmm12				\n\t"\
		"																unpcklpd	 0x1c0(%%rbx),%%xmm8 				\n\t"\
		"																movaps		%%xmm12, 0x1c0(%%rbx)				\n\t"\
		"																movaps		%%xmm14, 0x1d0(%%rax)				\n\t"\
		"																movaps		%%xmm8 , 0x1c0(%%rax)				\n\t"\
		"/************************************************/	\n\t		/************************************************/	\n\t"\
		"/*...Block 3: t04,t14,t24,t34 -> r20,24,22,26:  */	\n\t		/*...Block 7: t0C,t1C,t2C,t3C -> r28,2C,2A,2E:  */	\n\t"\
		"/************************************************/	\n\t		/************************************************/	\n\t"\
		"addq	$0x200,%%rdx	/* __r20 */					/* base-addr in rcol = r28, so rdx offset +0x80 vs lcol */	\n\t"\
		"leaq	0x010(%%rsi),%%rcx	/* cc0; Note cc0/ss0 are shared between lcol/rcol, so no rcx-offset until get to twiddles*/\n\t"\
		"movaps		 0x020(%%rdx),%%xmm4					\n\t		movaps		 0x0a0(%%rdx),%%xmm12				\n\t"\
		"movaps		 0x060(%%rdx),%%xmm0					\n\t		movaps		 0x0e0(%%rdx),%%xmm8 				\n\t"\
		"movaps		 0x030(%%rdx),%%xmm5					\n\t		movaps		 0x0b0(%%rdx),%%xmm13				\n\t"\
		"movaps		 0x070(%%rdx),%%xmm1					\n\t		movaps		 0x0f0(%%rdx),%%xmm9 				\n\t"\
		"movaps		 0x020(%%rdx),%%xmm6					\n\t		movaps		 0x0a0(%%rdx),%%xmm14				\n\t"\
		"movaps		 0x060(%%rdx),%%xmm2					\n\t		movaps		 0x0e0(%%rdx),%%xmm10				\n\t"\
		"movaps		 0x030(%%rdx),%%xmm7					\n\t		movaps		 0x0b0(%%rdx),%%xmm15				\n\t"\
		"movaps		 0x070(%%rdx),%%xmm3					\n\t		movaps		 0x0f0(%%rdx),%%xmm11				\n\t"\
		"mulpd			  (%%rcx),%%xmm4					\n\t		mulpd		 0x010(%%rcx),%%xmm12				\n\t"\
		"mulpd		 0x010(%%rcx),%%xmm0					\n\t		mulpd			  (%%rcx),%%xmm8 				\n\t"\
		"mulpd			  (%%rcx),%%xmm5					\n\t		mulpd		 0x010(%%rcx),%%xmm13				\n\t"\
		"mulpd		 0x010(%%rcx),%%xmm1					\n\t		mulpd			  (%%rcx),%%xmm9 				\n\t"\
		"mulpd		 0x010(%%rcx),%%xmm6					\n\t		mulpd			  (%%rcx),%%xmm14				\n\t"\
		"mulpd			  (%%rcx),%%xmm2					\n\t		mulpd		 0x010(%%rcx),%%xmm10				\n\t"\
		"mulpd		 0x010(%%rcx),%%xmm7					\n\t		mulpd			  (%%rcx),%%xmm15				\n\t"\
		"mulpd			  (%%rcx),%%xmm3					\n\t		mulpd		 0x010(%%rcx),%%xmm11				\n\t"\
		"subpd		%%xmm6,%%xmm5						\n\t		subpd		%%xmm14,%%xmm13					\n\t"\
		"subpd		%%xmm2,%%xmm1						\n\t		subpd		%%xmm10,%%xmm9 					\n\t"\
		"addpd		%%xmm7,%%xmm4						\n\t		addpd		%%xmm15,%%xmm12					\n\t"\
		"addpd		%%xmm3,%%xmm0						\n\t		addpd		%%xmm11,%%xmm8 					\n\t"\
		"movaps		%%xmm5,%%xmm7						\n\t		movaps		%%xmm13,%%xmm15					\n\t"\
		"movaps		%%xmm4,%%xmm6						\n\t		movaps		%%xmm12,%%xmm14					\n\t"\
		"addpd		%%xmm0,%%xmm4						\n\t		addpd		%%xmm8 ,%%xmm12					\n\t"\
		"addpd		%%xmm1,%%xmm5						\n\t		addpd		%%xmm9 ,%%xmm13					\n\t"\
		"subpd		%%xmm0,%%xmm6						\n\t		subpd		%%xmm8 ,%%xmm14					\n\t"\
		"subpd		%%xmm1,%%xmm7						\n\t		subpd		%%xmm9 ,%%xmm15					\n\t"\
		"movaps		 0x040(%%rdx),%%xmm2					\n\t		movaps		 0x0c0(%%rdx),%%xmm10				\n\t"\
		"movaps		 0x050(%%rdx),%%xmm3					\n\t		movaps		 0x0d0(%%rdx),%%xmm11				\n\t"\
		"movaps			  (%%rdx),%%xmm0					\n\t		movaps		 0x080(%%rdx),%%xmm8 				\n\t"\
		"movaps		 0x010(%%rdx),%%xmm1					\n\t		movaps		 0x090(%%rdx),%%xmm9 				\n\t"\
		"addpd		 0x050(%%rdx),%%xmm2					\n\t		subpd		 0x0d0(%%rdx),%%xmm10				\n\t"\
		"subpd		 0x040(%%rdx),%%xmm3					\n\t		addpd		 0x0c0(%%rdx),%%xmm11				\n\t"\
		"mulpd			  (%%rsi),%%xmm2					\n\t		mulpd			  (%%rsi),%%xmm10				\n\t"\
		"mulpd			  (%%rsi),%%xmm3					\n\t		mulpd			  (%%rsi),%%xmm11				\n\t"\
		"subpd		%%xmm2,%%xmm0						\n\t		subpd		%%xmm10,%%xmm8 					\n\t"\
		"subpd		%%xmm3,%%xmm1						\n\t		subpd		%%xmm11,%%xmm9 					\n\t"\
		"addpd		%%xmm2,%%xmm2						\n\t		addpd		%%xmm10,%%xmm10					\n\t"\
		"addpd		%%xmm3,%%xmm3						\n\t		addpd		%%xmm11,%%xmm11					\n\t"\
		"addpd		%%xmm0,%%xmm2						\n\t		addpd		%%xmm8 ,%%xmm10					\n\t"\
		"addpd		%%xmm1,%%xmm3						\n\t		addpd		%%xmm9 ,%%xmm11					\n\t"\
		"movq		%[__add0],%%rax					\n\t"\
		"movq		%[__c02],%%rcx					\n\t"\
		"subpd		%%xmm4,%%xmm2						\n\t		subpd		%%xmm14,%%xmm8 					\n\t"\
		"subpd		%%xmm5,%%xmm3						\n\t		subpd		%%xmm15,%%xmm9 					\n\t"\
		"addpd		%%xmm4,%%xmm4						\n\t		addpd		%%xmm14,%%xmm14					\n\t"\
		"addpd		%%xmm5,%%xmm5						\n\t		addpd		%%xmm15,%%xmm15					\n\t"\
		"addpd		%%xmm2,%%xmm4						\n\t		addpd		%%xmm8 ,%%xmm14					\n\t"\
		"addpd		%%xmm3,%%xmm5						\n\t		addpd		%%xmm9 ,%%xmm15					\n\t"\
		"movaps		%%xmm2,      (%%rdx)				\n\t		movaps		%%xmm8 , 0x080(%%rdx)				\n\t"\
		"movaps		%%xmm3, 0x010(%%rdx)				\n\t		movaps		%%xmm9 , 0x090(%%rdx)				\n\t"\
		"movaps		%%xmm4,%%xmm2						\n\t		movaps		%%xmm14,%%xmm8 					\n\t"\
		"movaps		%%xmm5,%%xmm3						\n\t		movaps		%%xmm15,%%xmm9 					\n\t"\
		"mulpd			  (%%rcx),%%xmm4					\n\t		mulpd		 0x080(%%rcx),%%xmm14	/* c06 */	\n\t"\
		"mulpd			  (%%rcx),%%xmm5					\n\t		mulpd		 0x080(%%rcx),%%xmm15				\n\t"\
		"mulpd		 0x010(%%rcx),%%xmm2					\n\t		mulpd		 0x090(%%rcx),%%xmm8 				\n\t"\
		"mulpd		 0x010(%%rcx),%%xmm3					\n\t		mulpd		 0x090(%%rcx),%%xmm9 				\n\t"\
		"subpd		%%xmm2,%%xmm5						\n\t		subpd		%%xmm8 ,%%xmm15					\n\t"\
		"addpd		%%xmm3,%%xmm4						\n\t		addpd		%%xmm9 ,%%xmm14					\n\t"\
		"movq		%[__add1],%%rbx					\n\t"\
		"movaps		%%xmm5,%%xmm3						\n\t		movaps		%%xmm15,%%xmm9 					\n\t"\
		"movaps		%%xmm4,%%xmm2						\n\t		movaps		%%xmm14,%%xmm8 					\n\t"\
		"unpckhpd	 0x030(%%rbx),%%xmm3					\n\t		unpckhpd	 0x070(%%rbx),%%xmm9 				\n\t"\
		"unpcklpd	 0x030(%%rbx),%%xmm5					\n\t		unpcklpd	 0x070(%%rbx),%%xmm15				\n\t"\
		"movaps		%%xmm3, 0x030(%%rbx)				\n\t		movaps		%%xmm9 , 0x070(%%rbx)				\n\t"\
		"unpckhpd	 0x020(%%rbx),%%xmm2					\n\t		unpckhpd	 0x060(%%rbx),%%xmm8 				\n\t"\
		"unpcklpd	 0x020(%%rbx),%%xmm4					\n\t		unpcklpd	 0x060(%%rbx),%%xmm14				\n\t"\
		"movaps		%%xmm2, 0x020(%%rbx)				\n\t		movaps		%%xmm8 , 0x060(%%rbx)				\n\t"\
		"movaps		%%xmm5, 0x030(%%rax)				\n\t		movaps		%%xmm15, 0x070(%%rax)				\n\t"\
		"movaps		%%xmm4, 0x020(%%rax)				\n\t		movaps		%%xmm14, 0x060(%%rax)				\n\t"\
		"movq		%[__c12],%%rcx					\n\t"\
		"movaps			  (%%rdx),%%xmm4					\n\t		movaps		 0x080(%%rdx),%%xmm14				\n\t"\
		"movaps		 0x010(%%rdx),%%xmm5					\n\t		movaps		 0x090(%%rdx),%%xmm15				\n\t"\
		"movaps		%%xmm4,%%xmm2						\n\t		movaps		%%xmm14,%%xmm8 					\n\t"\
		"movaps		%%xmm5,%%xmm3						\n\t		movaps		%%xmm15,%%xmm9 					\n\t"\
		"mulpd			  (%%rcx),%%xmm4					\n\t		mulpd		 0x080(%%rcx),%%xmm14	/* c16 */	\n\t"\
		"mulpd			  (%%rcx),%%xmm5					\n\t		mulpd		 0x080(%%rcx),%%xmm15				\n\t"\
		"mulpd		 0x010(%%rcx),%%xmm2					\n\t		mulpd		 0x090(%%rcx),%%xmm8 				\n\t"\
		"mulpd		 0x010(%%rcx),%%xmm3					\n\t		mulpd		 0x090(%%rcx),%%xmm9 				\n\t"\
		"subpd		%%xmm2,%%xmm5						\n\t		subpd		%%xmm8 ,%%xmm15					\n\t"\
		"addpd		%%xmm3,%%xmm4						\n\t		addpd		%%xmm9 ,%%xmm14					\n\t"\
		"movaps		%%xmm5,%%xmm3						\n\t		movaps		%%xmm15,%%xmm9 					\n\t"\
		"movaps		%%xmm4,%%xmm2						\n\t		movaps		%%xmm14,%%xmm8 					\n\t"\
		"unpckhpd	 0x130(%%rbx),%%xmm3					\n\t		unpckhpd	 0x170(%%rbx),%%xmm9 				\n\t"\
		"unpcklpd	 0x130(%%rbx),%%xmm5					\n\t		unpcklpd	 0x170(%%rbx),%%xmm15				\n\t"\
		"movaps		%%xmm3, 0x130(%%rbx)				\n\t		movaps		%%xmm9 , 0x170(%%rbx)				\n\t"\
		"unpckhpd	 0x120(%%rbx),%%xmm2					\n\t		unpckhpd	 0x160(%%rbx),%%xmm8 				\n\t"\
		"unpcklpd	 0x120(%%rbx),%%xmm4					\n\t		unpcklpd	 0x160(%%rbx),%%xmm14				\n\t"\
		"movaps		%%xmm2, 0x120(%%rbx)				\n\t		movaps		%%xmm8 , 0x160(%%rbx)				\n\t"\
		"movaps		%%xmm5, 0x130(%%rax)				\n\t		movaps		%%xmm15, 0x170(%%rax)				\n\t"\
		"movaps		%%xmm4, 0x120(%%rax)				\n\t		movaps		%%xmm14, 0x160(%%rax)				\n\t"\
		"movq		%[__c0A],%%rcx					\n\t"\
		"subpd		%%xmm7,%%xmm0						\n\t		subpd		%%xmm13,%%xmm10					\n\t"\
		"subpd		%%xmm6,%%xmm1						\n\t		subpd		%%xmm12,%%xmm11					\n\t"\
		"addpd		%%xmm7,%%xmm7						\n\t		addpd		%%xmm13,%%xmm13					\n\t"\
		"addpd		%%xmm6,%%xmm6						\n\t		addpd		%%xmm12,%%xmm12					\n\t"\
		"addpd		%%xmm0,%%xmm7						\n\t		addpd		%%xmm10,%%xmm13					\n\t"\
		"addpd		%%xmm1,%%xmm6						\n\t		addpd		%%xmm11,%%xmm12					\n\t"\
		"movaps		%%xmm7,%%xmm4						\n\t		movaps		%%xmm13,%%xmm8 					\n\t"\
		"movaps		%%xmm1,%%xmm5						\n\t		movaps		%%xmm11,%%xmm9 					\n\t"\
		"mulpd			  (%%rcx),%%xmm7					\n\t		mulpd		 0x080(%%rcx),%%xmm13	/* c0E */	\n\t"\
		"mulpd			  (%%rcx),%%xmm1					\n\t		mulpd		 0x080(%%rcx),%%xmm11				\n\t"\
		"mulpd		 0x010(%%rcx),%%xmm4					\n\t		mulpd		 0x090(%%rcx),%%xmm8 				\n\t"\
		"mulpd		 0x010(%%rcx),%%xmm5					\n\t		mulpd		 0x090(%%rcx),%%xmm9 				\n\t"\
		"subpd		%%xmm4,%%xmm1						\n\t		subpd		%%xmm8 ,%%xmm11					\n\t"\
		"addpd		%%xmm5,%%xmm7						\n\t		addpd		%%xmm9 ,%%xmm13					\n\t"\
		"movaps		%%xmm1,%%xmm5						\n\t		movaps		%%xmm11,%%xmm9 					\n\t"\
		"movaps		%%xmm7,%%xmm4						\n\t		movaps		%%xmm13,%%xmm8 					\n\t"\
		"unpckhpd	 0x0b0(%%rbx),%%xmm5					\n\t		unpckhpd	 0x0f0(%%rbx),%%xmm9 				\n\t"\
		"unpcklpd	 0x0b0(%%rbx),%%xmm1					\n\t		unpcklpd	 0x0f0(%%rbx),%%xmm11				\n\t"\
		"movaps		%%xmm5, 0x0b0(%%rbx)				\n\t		movaps		%%xmm9 , 0x0f0(%%rbx)				\n\t"\
		"unpckhpd	 0x0a0(%%rbx),%%xmm4					\n\t		unpckhpd	 0x0e0(%%rbx),%%xmm8 				\n\t"\
		"unpcklpd	 0x0a0(%%rbx),%%xmm7					\n\t		unpcklpd	 0x0e0(%%rbx),%%xmm13				\n\t"\
		"movaps		%%xmm4, 0x0a0(%%rbx)				\n\t		movaps		%%xmm8 , 0x0e0(%%rbx)				\n\t"\
		"movaps		%%xmm1, 0x0b0(%%rax)				\n\t		movaps		%%xmm11, 0x0f0(%%rax)				\n\t"\
		"movaps		%%xmm7, 0x0a0(%%rax)				\n\t		movaps		%%xmm13, 0x0e0(%%rax)				\n\t"\
		"movq		%[__c1A],%%rcx					\n\t"\
		"movaps		%%xmm0,%%xmm4						\n\t		movaps		%%xmm10,%%xmm8 					\n\t"\
		"movaps		%%xmm6,%%xmm5						\n\t		movaps		%%xmm12,%%xmm9 					\n\t"\
		"mulpd			  (%%rcx),%%xmm0					\n\t		mulpd		 0x080(%%rcx),%%xmm10	/* c1E */	\n\t"\
		"mulpd			  (%%rcx),%%xmm6					\n\t		mulpd		 0x080(%%rcx),%%xmm12				\n\t"\
		"mulpd		 0x010(%%rcx),%%xmm4					\n\t		mulpd		 0x090(%%rcx),%%xmm8 				\n\t"\
		"mulpd		 0x010(%%rcx),%%xmm5					\n\t		mulpd		 0x090(%%rcx),%%xmm9 				\n\t"\
		"subpd		%%xmm4,%%xmm6						\n\t		subpd		%%xmm8 ,%%xmm12					\n\t"\
		"addpd		%%xmm5,%%xmm0						\n\t		addpd		%%xmm9 ,%%xmm10					\n\t"\
		"movaps		%%xmm6,%%xmm5						\n\t		movaps		%%xmm12,%%xmm9 					\n\t"\
		"movaps		%%xmm0,%%xmm4						\n\t		movaps		%%xmm10,%%xmm8 					\n\t"\
		"unpckhpd	 0x1b0(%%rbx),%%xmm5					\n\t		unpckhpd	 0x1f0(%%rbx),%%xmm9 				\n\t"\
		"unpcklpd	 0x1b0(%%rbx),%%xmm6					\n\t		unpcklpd	 0x1f0(%%rbx),%%xmm12				\n\t"\
		"movaps		%%xmm5, 0x1b0(%%rbx)				\n\t		movaps		%%xmm9 , 0x1f0(%%rbx)				\n\t"\
		"unpckhpd	 0x1a0(%%rbx),%%xmm4					\n\t		unpckhpd	 0x1e0(%%rbx),%%xmm8 				\n\t"\
		"unpcklpd	 0x1a0(%%rbx),%%xmm0					\n\t		unpcklpd	 0x1e0(%%rbx),%%xmm10				\n\t"\
		"movaps		%%xmm4, 0x1a0(%%rbx)				\n\t		movaps		%%xmm8 , 0x1e0(%%rbx)				\n\t"\
		"movaps		%%xmm6, 0x1b0(%%rax)				\n\t		movaps		%%xmm12, 0x1f0(%%rax)				\n\t"\
		"movaps		%%xmm0, 0x1a0(%%rax)				\n\t		movaps		%%xmm10, 0x1e0(%%rax)				\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__add1] "m" (Xadd1)\
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__r00] "m" (Xr00)\
		 ,[__r08] "m" (Xr08)\
		 ,[__r10] "m" (Xr10)\
		 ,[__c01] "m" (Xc01)\
		 ,[__c02] "m" (Xc02)\
		 ,[__c04] "m" (Xc04)\
		 ,[__c06] "m" (Xc06)\
		 ,[__c08] "m" (Xc08)\
		 ,[__c0A] "m" (Xc0A)\
		 ,[__c0C] "m" (Xc0C)\
		 ,[__c0E] "m" (Xc0E)\
		 ,[__c10] "m" (Xc10)\
		 ,[__c12] "m" (Xc12)\
		 ,[__c14] "m" (Xc14)\
		 ,[__c16] "m" (Xc16)\
		 ,[__c18] "m" (Xc18)\
		 ,[__c1A] "m" (Xc1A)\
		 ,[__c1C] "m" (Xc1C)\
		 ,[__c1E] "m" (Xc1E)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"		/* Clobbered registers */\
	);\
	}

#endif	// SSE2 or AVX?

#endif	/* radix32_wrapper_square_gcc_h_included */

