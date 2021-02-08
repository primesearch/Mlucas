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
#ifndef radix8_dif_dit_pass_asm_h_included
#define radix8_dif_dit_pass_asm_h_included

#ifdef USE_ARM_V8_SIMD

	// Made the mistake here of using the SSE2 macro as a coding template - that was so munged based on that ISA
	// (and 16 vector registers) that it proved a royal PITA, ended up 'unscrambling' most of the trickier parts
	// of the SSE2-asm and having to track all the associated spill/fills and reg-indexing ugliness:
	#define SSE2_RADIX8_DIF_TWIDDLE(Xadd0,Xadd1,Xadd2,Xadd3,Xadd4,Xadd5,Xadd6,Xadd7,Xisrt2,Xc1,Xc2,Xc3,Xc4,Xc5,Xc6,Xc7)\
	{\
	__asm__ volatile (\
		"ldr	x10,%[isrt2]		\n\t	ld1r	{v29.2d},[x10]	\n\t"/* isrt2 */\
										"	ldr	x1,%[add1]			\n\t"\
										"	ldr	x5,%[add5]			\n\t"\
		"ldr	x0,%[add0]			\n\t	ldr	x12,%[c1]			\n\t"\
		"ldr	x4,%[add4]			\n\t	ldr	x13,%[c5]			\n\t"\
		"ldr	x10,%[c4]			\n\t	ldp	q12,q13,[x1]			\n\t"\
		"ldp	q0,q1,[x0]			\n\t	ldp	q14,q15,[x5]			\n\t"\
		"ldp	q2,q3,[x4]			\n\t	ldp	q16,q17,[x12]			\n\t"/* c1 */\
/* c4 */"ldp	q8,q9,[x10]			\n\t	fmul	v18.2d,v12.2d,v16.2d	\n\t"\
		"fmul	v4.2d,v2.2d,v8.2d	\n\t	fmul	v19.2d,v13.2d,v16.2d	\n\t"\
		"fmul	v5.2d,v3.2d,v8.2d	\n\t	fmls	v18.2d,v13.2d,v17.2d	\n\t"\
		"fmls	v4.2d,v3.2d,v9.2d	\n\t	fmla	v19.2d,v12.2d,v17.2d	\n\t"\
		"fmla	v5.2d,v2.2d,v9.2d	\n\t	ldp	q16,q17,[x13]			\n\t"/* c5 */\
		"fsub	v2.2d,v0.2d,v4.2d	\n\t	fmul	v12.2d,v14.2d,v16.2d	\n\t"\
		"fsub	v3.2d,v1.2d,v5.2d	\n\t	fmul	v13.2d,v15.2d,v16.2d	\n\t"\
		"fadd	v0.2d,v0.2d,v4.2d	\n\t	fmls	v12.2d,v15.2d,v17.2d	\n\t"\
		"fadd	v1.2d,v1.2d,v5.2d	\n\t	fmla	v13.2d,v14.2d,v17.2d	\n\t"\
		"ldr	x2,%[add2]			\n\t	fsub	v14.2d,v18.2d,v12.2d	\n\t"\
		"ldr	x6,%[add6]			\n\t	fsub	v15.2d,v19.2d,v13.2d	\n\t"\
		"ldr	x10,%[c2]			\n\t	fadd	v12.2d,v18.2d,v12.2d	\n\t"\
		"ldr	x11,%[c6]			\n\t	fadd	v13.2d,v19.2d,v13.2d	\n\t"\
		"ldp	q4,q5,[x2]			\n\t	ldr	x3,%[add3]			\n\t"\
		"ldp	q6,q7,[x6]			\n\t	ldr	x7,%[add7]			\n\t"\
/* c2 */"ldp	q8,q9,[x10]			\n\t	ldr	x12,%[c3]			\n\t"\
		"fmul	v10.2d,v4.2d,v8.2d	\n\t	ldr	x13,%[c7]			\n\t"\
		"fmul	v11.2d,v5.2d,v8.2d	\n\t	ldp	q16,q17,[x3]			\n\t"\
		"fmls	v10.2d,v5.2d,v9.2d	\n\t	ldp	q18,q19,[x7]			\n\t"\
		"fmla	v11.2d,v4.2d,v9.2d	\n\t	ldp	q20,q21,[x12]			\n\t"/* c3 */\
/* c6 */"ldp	q8,q9,[x11]			\n\t	fmul	v22.2d,v16.2d,v20.2d	\n\t"\
		"fmul	v4.2d,v6.2d,v8.2d	\n\t	fmul	v23.2d,v17.2d,v20.2d	\n\t"\
		"fmul	v5.2d,v7.2d,v8.2d	\n\t	fmls	v22.2d,v17.2d,v21.2d	\n\t"\
		"fmls	v4.2d,v7.2d,v9.2d	\n\t	fmla	v23.2d,v16.2d,v21.2d	\n\t"\
		"fmla	v5.2d,v6.2d,v9.2d	\n\t	ldp	q20,q21,[x13]			\n\t"/* c7 */\
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
		"fsub	v6.2d,v0.2d,v12.2d	\n\t"/* 1r: xmm0 = 0-12 */	"fsub	v18.2d,v14.2d,v15.2d	\n\t"/* 18 =14-15 */\
		"fsub	v7.2d,v1.2d,v13.2d	\n\t"/* 1i: xmm1 = 1-13 */	"fadd	v19.2d,v14.2d,v15.2d	\n\t"/* 19 =14+15 */\
		"fadd	v0.2d,v0.2d,v12.2d	\n\t"/* 0r: xmm2 = 0+12 */	"fmul	v18.2d,v18.2d,v29.2d	\n\t"/* 18 *=isrt2 */\
		"fadd	v1.2d,v1.2d,v13.2d	\n\t"/* 0i: xmm3 = 1+13 */	"fmul	v19.2d,v19.2d,v29.2d	\n\t"/* 19 *=isrt2 */\
		"fsub	v12.2d,v8.2d,v21.2d	\n\t"/* 2r: xmm4 = 8-21 */	"fadd	v14.2d,v16.2d,v17.2d	\n\t"/* 14=16+17 */\
		"fsub	v13.2d,v9.2d,v20.2d	\n\t"/* 3i: xmm5 = 9-20 */	"fsub	v15.2d,v16.2d,v17.2d	\n\t"/* 15=16-17 */\
		"fadd	v8.2d,v8.2d,v21.2d	\n\t"/* 3r: xmm7 = 8+21 */	"fmul	v14.2d,v14.2d,v29.2d	\n\t"/* 14 *=isrt2 */\
		"fadd	v9.2d,v9.2d,v20.2d	\n\t"/* 2i: xmm6 = 9+20 */	"fmul	v15.2d,v15.2d,v29.2d	\n\t"/* 15 *=isrt2 */\
																"fsub	v16.2d,v2.2d ,v18.2d	\n\t"/* 5r: xmm8 =2-18 */\
																"fsub	v17.2d,v3.2d ,v19.2d	\n\t"/* 5i: xmm9 =3-19 */\
																"fadd	v2.2d ,v2.2d ,v18.2d	\n\t"/* 4r: xmm10=2+18 */\
																"fadd	v3.2d ,v3.2d ,v19.2d	\n\t"/* 4i: xmm13=3+19 */\
																"fsub	v18.2d,v4.2d ,v14.2d	\n\t"/* 6r: xmm14=4-14 */\
																"fsub	v19.2d,v5.2d ,v15.2d	\n\t"/* 6i: xmm11=5-15 */\
																"fadd	v4.2d ,v4.2d ,v14.2d	\n\t"/* 7r: xmm12=4+14 */\
																"fadd	v5.2d ,v5.2d ,v15.2d	\n\t"/* 7i: xmm15=5+15 */\
		"stp	q0 ,q1 ,[x0]		\n\t	stp	q2 ,q3 ,[x4]			\n\t"\
		"stp	q6 ,q7 ,[x1]		\n\t	stp	q16,q17,[x5]			\n\t"\
		"stp	q12,q9 ,[x2]		\n\t	stp	q18,q19,[x6]			\n\t"\
		"stp	q8 ,q13,[x3]		\n\t	stp	q4 ,q5 ,[x7]			\n\t"\
		:					/* outputs: none */\
		: [add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[add1] "m" (Xadd1)\
		 ,[add2] "m" (Xadd2)\
		 ,[add3] "m" (Xadd3)\
		 ,[add4] "m" (Xadd4)\
		 ,[add5] "m" (Xadd5)\
		 ,[add6] "m" (Xadd6)\
		 ,[add7] "m" (Xadd7)\
		 ,[isrt2] "m" (Xisrt2)\
		 ,[c1] "m" (Xc1)\
		 ,[c2] "m" (Xc2)\
		 ,[c3] "m" (Xc3)\
		 ,[c4] "m" (Xc4)\
		 ,[c5] "m" (Xc5)\
		 ,[c6] "m" (Xc6)\
		 ,[c7] "m" (Xc7)\
		: "cc","memory","x0","x1","x2","x3","x4","x5","x6","x7","x10","x11","x12","x13","v0","v1","v2","v3","v4","v5","v6","v7","v8","v9",\
		"v10","v11","v12","v13","v14","v15","v16","v17","v18","v19","v20","v21","v22","v23","v24","v25", "v29"	/* Clobbered registers */\
	);\
	}

	// Based on the aforementioned experience of using the SSE2 macro as a coding template for the radix-8 DIF macro,
	// for DIT instead to directly using the C code as a template, i.e. doing it from scratch, as it were:
	#define SSE2_RADIX8_DIT_TWIDDLE(Xadd0,Xadd1,Xadd2,Xadd3,Xadd4,Xadd5,Xadd6,Xadd7,Xisrt2,Xc1,Xc2,Xc3,Xc4,Xc5,Xc6,Xc7)\
	{\
	__asm__ volatile (\
/*vvv*********** in terms of the C-code's t-indexes **********/\
		"ldr	x10,%[isrt2]		\n\t	ld1r	{v29.2d},[x10]	\n\t"/* isrt2 */\
		"ldr	x0,%[add0]			\n\t	ldr	x4,%[add4]			\n\t"\
		"ldr	x1,%[add1]			\n\t	ldr	x5,%[add5]			\n\t"\
		"ldr	x2,%[add2]			\n\t	ldr	x6,%[add6]			\n\t"\
		"ldr	x3,%[add3]			\n\t	ldr	x7,%[add7]			\n\t"\
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
	/* combine to get the 2 length-4 transform:
		rt =t5;	t5 =t1 -rt;	t1 =t1 +rt;		rt =t13; t13=t9 -rt; t9 =t9 +rt;
		it =t6;	t6 =t2 -it;	t2 =t2 +it;		it =t14; t14=t10-it; t10=t10+it; */\
		"fsub	v8.2d,v0.2d,v4.2d	\n\t	fsub	v20.2d,v12.2d,v16.2d	\n\t"/* t5,13 */ \
		"fsub	v9.2d,v1.2d,v5.2d	\n\t	fsub	v21.2d,v13.2d,v17.2d	\n\t"/* t6,14 */ \
		"fadd	v0.2d,v0.2d,v4.2d	\n\t	fadd	v12.2d,v12.2d,v16.2d	\n\t"/* t1,9  */ \
		"fadd	v1.2d,v1.2d,v5.2d	\n\t	fadd	v13.2d,v13.2d,v17.2d	\n\t"/* t2,10 */ \
	/*	rt =t7;	t7 =t3 -t8;	t3 =t3 +t8;		rt =t15; t15=t11-t16; t11=t11+t16;
			t8 =t4 +rt;	t4 =t4 -rt;					 t16=t12+rt ; t12=t12-rt ; */\
		"fadd	v4.2d,v2.2d,v7.2d	\n\t	fadd	v16.2d,v14.2d,v19.2d	\n\t"/* t3,11 */\
		"fsub	v5.2d,v3.2d,v6.2d	\n\t	fsub	v17.2d,v15.2d,v18.2d	\n\t"/* t4,12 */\
		"fsub	v2.2d,v2.2d,v7.2d	\n\t	fsub	v14.2d,v14.2d,v19.2d	\n\t"/* t7,15 */\
		"fadd	v3.2d,v3.2d,v6.2d	\n\t	fadd	v15.2d,v15.2d,v18.2d	\n\t"/* t8,16 */\
	/* now combine the two half-transforms: */\
		"ldr	x10,%[c4]			\n\t	ldp	q10,q11,[x10]			\n\t"/* c4 */\
	  /* [0]	a[j1   ]=t1+t9;				a[j2   ]=t2+t10;
				t1      =t1-t9;				t2      =t2-t10; */\
		"fadd	v6.2d,v0.2d,v12.2d	\n\t	fadd	v18.2d,v1.2d,v13.2d	\n\t"\
		"fsub	v7.2d,v0.2d,v12.2d	\n\t	fsub	v19.2d,v1.2d,v13.2d	\n\t"/* t1,2 = v7,19 */\
	  /*	a[j1+p4]=t1 *c4 +t2 *s4;	a[j2+p4]=t2 *c4 -t1 *s4; */\
		"fmul	v0.2d,v7.2d,v10.2d	\n\t	fmul	v1.2d,v19.2d,v10.2d	\n\t"\
		"fmla	v0.2d,v19.2d,v11.2d	\n\t	fmls	v1.2d,v7.2d,v11.2d	\n\t"\
		"stp	q6,q18,[x0]			\n\t	stp	q0,q1,[x4]			\n\t"\
	  /* [1]	rt=(t11+t12)*ISRT2;			it=(t11-t12)*ISRT2;
				t11     =t3+rt;				t12       =t4-it;
				t3      =t3-rt;				t4        =t4+it; */\
		"ldr	x10,%[c1]			\n\t	ldp	q10,q11,[x10]			\n\t"/* c1 */\
		"ldr	x11,%[c5]			\n\t	ldp	q12,q13,[x11]			\n\t"/* c5 */\
		"fadd	v6.2d,v16.2d,v17.2d	\n\t	fsub	v7.2d,v16.2d,v17.2d	\n\t"/* t3,4 = v4,5; t11,12 = v16,17 */\
		"fmul	v16.2d,v29.2d,v6.2d	\n\t	fmul	v17.2d,v29.2d,v7.2d	\n\t"/* rt,it; */\
		"fadd	v6.2d,v4.2d,v16.2d	\n\t	fsub	v18.2d,v5.2d,v17.2d	\n\t"/* t11,12 = v6,18 */\
		"fsub	v7.2d,v4.2d,v16.2d	\n\t	fadd	v19.2d,v5.2d,v17.2d	\n\t"/* t3 ,4  = v7,19 */\
	  /*	a[j1+p1]=t11*c1 +t12*s1;	a[j2+p1]=t12*c1 -t11*s1;
			a[j1+p5]=t3 *c5 +t4 *s5;	a[j2+p5]=t4 *c5 -t3 *s5; */\
		"fmul	v0.2d,v6.2d,v10.2d	\n\t	fmul	v1.2d,v18.2d,v10.2d	\n\t"\
		"fmul	v4.2d,v7.2d,v12.2d	\n\t	fmul	v5.2d,v19.2d,v12.2d	\n\t"\
		"fmla	v0.2d,v18.2d,v11.2d	\n\t	fmls	v1.2d,v6.2d,v11.2d	\n\t"\
		"fmla	v4.2d,v19.2d,v13.2d	\n\t	fmls	v5.2d,v7.2d,v13.2d	\n\t"\
		"stp	q0,q1,[x1]			\n\t	stp	q4,q5,[x5]			\n\t"\
	  /*rt      =t5+t14;			it        =t6-t13;
		t5      =t5-t14;			t6        =t6+t13;
		a[j1+p2]=rt *c2 +it *s2;	a[j2+p2]=it *c2 -rt *s2;
		a[j1+p6]=t5 *c6 +t6 *s6;	a[j2+p6]=t6 *c6 -t5 *s6;	ARM: On input, t5,6 = v8,9; t13,14 = v20,21: */\
		"ldr	x10,%[c2]			\n\t	ldp	q10,q11,[x10]			\n\t"/* c2 */\
		"ldr	x11,%[c6]			\n\t	ldp	q12,q13,[x11]			\n\t"/* c6 */\
		"fadd	v6.2d,v8.2d,v21.2d	\n\t	fsub	v18.2d,v9.2d,v20.2d	\n\t"/* rt,it = v6,18 */\
		"fsub	v7.2d,v8.2d,v21.2d	\n\t	fadd	v19.2d,v9.2d,v20.2d	\n\t"/* t3,4  = v7,19 */\
		"fmul	v0.2d,v6.2d,v10.2d	\n\t	fmul	v1.2d,v18.2d,v10.2d	\n\t"\
		"fmul	v4.2d,v7.2d,v12.2d	\n\t	fmul	v5.2d,v19.2d,v12.2d	\n\t"\
		"fmla	v0.2d,v18.2d,v11.2d	\n\t	fmls	v1.2d,v6.2d,v11.2d	\n\t"\
		"fmla	v4.2d,v19.2d,v13.2d	\n\t	fmls	v5.2d,v7.2d,v13.2d	\n\t"\
		"stp	q0,q1,[x2]			\n\t	stp	q4,q5,[x6]			\n\t"\
	  /*rt=(t15-t16)*ISRT2;			it=(t15+t16)*ISRT2;
		t15     =t7-rt;				t16       =t8-it;
		t7      =t7+rt;				t8        =t8+it; */\
		"ldr	x10,%[c3]			\n\t	ldp	q10,q11,[x10]			\n\t"/* c3 */\
		"ldr	x11,%[c7]			\n\t	ldp	q12,q13,[x11]			\n\t"/* c7 */\
		"fsub	v6.2d,v14.2d,v15.2d	\n\t	fadd	v7.2d,v14.2d,v15.2d	\n\t"/* t7,8 = v2,3; t15,16 = v14,15 */\
		"fmul	v14.2d,v29.2d,v6.2d	\n\t	fmul	v15.2d,v29.2d,v7.2d	\n\t"/* rt,it; */\
		"fsub	v6.2d,v2.2d,v14.2d	\n\t	fsub	v18.2d,v3.2d,v15.2d	\n\t"/* t15,16 = v6,18 */\
		"fadd	v7.2d,v2.2d,v14.2d	\n\t	fadd	v19.2d,v3.2d,v15.2d	\n\t"/* t7 ,8  = v7,19 */\
	  /*a[j1+p3]=t15*c3 +t16*s3;	a[j2+p3]=t16*c3 -t15*s3;
		a[j1+p7]=t7 *c7 +t8 *s7;	a[j2+p7]=t8 *c7 -t7 *s7; */\
		"fmul	v0.2d,v6.2d,v10.2d	\n\t	fmul	v1.2d,v18.2d,v10.2d	\n\t"\
		"fmul	v4.2d,v7.2d,v12.2d	\n\t	fmul	v5.2d,v19.2d,v12.2d	\n\t"\
		"fmla	v0.2d,v18.2d,v11.2d	\n\t	fmls	v1.2d,v6.2d,v11.2d	\n\t"\
		"fmla	v4.2d,v19.2d,v13.2d	\n\t	fmls	v5.2d,v7.2d,v13.2d	\n\t"\
		"stp	q0,q1,[x3]			\n\t	stp	q4,q5,[x7]			\n\t"\
		:					/* outputs: none */\
		: [add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[add1] "m" (Xadd1)\
		 ,[add2] "m" (Xadd2)\
		 ,[add3] "m" (Xadd3)\
		 ,[add4] "m" (Xadd4)\
		 ,[add5] "m" (Xadd5)\
		 ,[add6] "m" (Xadd6)\
		 ,[add7] "m" (Xadd7)\
		 ,[isrt2] "m" (Xisrt2)\
		 ,[c1] "m" (Xc1)\
		 ,[c2] "m" (Xc2)\
		 ,[c3] "m" (Xc3)\
		 ,[c4] "m" (Xc4)\
		 ,[c5] "m" (Xc5)\
		 ,[c6] "m" (Xc6)\
		 ,[c7] "m" (Xc7)\
		: "cc","memory","x0","x1","x2","x3","x4","x5","x6","x7","x10","x11","x12","x13","v0","v1","v2","v3","v4","v5","v6","v7","v8","v9",\
		"v10","v11","v12","v13","v14","v15","v16","v17","v18","v19","v20","v21","v22","v23","v24","v25", "v29"	/* Clobbered registers */\
	);\
	}

#elif defined(USE_AVX512)	// Cf. avx2 code for commentary

	#define SSE2_RADIX8_DIF_TWIDDLE(Xadd0,Xadd1,Xadd2,Xadd3,Xadd4,Xadd5,Xadd6,Xadd7,Xisrt2,Xc1,Xc2,Xc3,Xc4,Xc5,Xc6,Xc7)\
	{\
	__asm__ volatile (\
	"movq	%[isrt2],%%rsi	\n\t"\
	"leaq	0x40(%%rsi),%%rsi	\n\t"/* 1.0 */\
		"							movq	%[add1]	,%%r10	\n\t	movq	%[c1],%%r12				\n\t"\
		"							movq	%[add5]	,%%r11	\n\t	movq	%[c5],%%r13				\n\t"\
		"										\n\t		vmovaps	    (%%r10)	,%%zmm8 			\n\t"\
		"										\n\t		vmovaps	0x40(%%r10)	,%%zmm10			\n\t"\
		"										\n\t		vmovaps	     %%zmm8	,%%zmm9 			\n\t"\
		"movq		%[add0]	,%%rax				\n\t		vmovaps	     %%zmm10,%%zmm11			\n\t"\
		"movq		%[add4]	,%%rbx				\n\t		vmulpd	     (%%r12),%%zmm8 ,%%zmm8 	\n\t"\
		"movq		%[c4]	,%%rcx				\n\t		vmulpd	 0x40(%%r12),%%zmm9 ,%%zmm9 	\n\t"\
		"vmovaps	    (%%rax)	,%%zmm0			\n\t	vfnmadd231pd 0x40(%%r12),%%zmm10,%%zmm8 	\n\t"\
		"vmovaps	0x40(%%rax)	,%%zmm1			\n\t	 vfmadd231pd     (%%r12),%%zmm11,%%zmm9 	\n\t"\
		"vmovaps	     %%zmm0	,%%zmm6			\n\t		vmovaps	    (%%r11)	,%%zmm10			\n\t"\
		"vmovaps	     %%zmm1	,%%zmm7			\n\t		vmovaps	0x40(%%r11)	,%%zmm11			\n\t"\
		"vmovaps	    (%%rbx)	,%%zmm2			\n\t		vmovaps	     %%zmm10,%%zmm12			\n\t"\
		"vmovaps	0x40(%%rbx)	,%%zmm3			\n\t		vmovaps	     %%zmm11,%%zmm13			\n\t"\
		"vmovaps	     %%zmm2	,%%zmm4			\n\t		vmulpd	     (%%r13),%%zmm10,%%zmm10	\n\t"\
		"vmovaps	     %%zmm3	,%%zmm5			\n\t		vmulpd	 0x40(%%r13),%%zmm12,%%zmm12	\n\t"\
		"vmulpd	        (%%rcx)	,%%zmm2,%%zmm2	\n\t	vfnmadd231pd 0x40(%%r13),%%zmm11,%%zmm10	\n\t"\
		"vmulpd	        (%%rcx)	,%%zmm3,%%zmm3	\n\t	 vfmadd231pd     (%%r13),%%zmm13,%%zmm12	\n\t"\
	"vfnmadd231pd	0x40(%%rcx)	,%%zmm5,%%zmm2	\n\t		vmovaps	%%zmm10		,%%zmm11			\n\t"\
	" vfmadd231pd	0x40(%%rcx)	,%%zmm4,%%zmm3	\n\t		vmovaps	%%zmm12		,%%zmm13			\n\t"\
	"vmovaps	(%%rsi),%%zmm4	\n\t"/* 1.0 */\
		"vfmadd132pd	%%zmm4,%%zmm2,%%zmm0	\n\t		vfmadd132pd	%%zmm4,%%zmm8 ,%%zmm10	\n\t"\
		"vfmadd132pd	%%zmm4,%%zmm3,%%zmm1	\n\t		vfmsub132pd	%%zmm4,%%zmm11,%%zmm8 	\n\t"\
		"vfmsub132pd	%%zmm4,%%zmm2,%%zmm6	\n\t		vfmadd132pd	%%zmm4,%%zmm9 ,%%zmm12	\n\t"\
		"vfmsub132pd	%%zmm4,%%zmm3,%%zmm7	\n\t		vfmsub132pd	%%zmm4,%%zmm13,%%zmm9 	\n\t"\
		"vmovaps	%%zmm0		,    (%%rax)	\n\t		vmovaps	%%zmm10		,    (%%r10)		\n\t"\
		"vmovaps	%%zmm1		,0x40(%%rax)	\n\t		vmovaps	%%zmm12		,0x40(%%r10)		\n\t"\
		"vmovaps	%%zmm6		,    (%%rbx)	\n\t		vmovaps	%%zmm8 		,    (%%r11)		\n\t"\
		"vmovaps	%%zmm7		,0x40(%%rbx)	\n\t		vmovaps	%%zmm9 		,0x40(%%r11)		\n\t"\
		"movq		%[add2]	,%%rax				\n\t		movq		%[add3]	,%%r10				\n\t"\
		"movq		%[add6]	,%%rbx				\n\t		movq		%[add7]	,%%r11				\n\t"\
		"movq		%[c2]	,%%rcx				\n\t		movq		%[c3]	,%%r12				\n\t"\
		"movq		%[c6]	,%%rdx				\n\t		movq		%[c7]	,%%r13				\n\t"\
		"vmovaps	    (%%rax)	,%%zmm0			\n\t		vmovaps	    (%%r10)	,%%zmm8 			\n\t"\
		"vmovaps	0x40(%%rax)	,%%zmm2			\n\t		vmovaps	0x40(%%r10)	,%%zmm10			\n\t"\
		"vmovaps	     %%zmm0	,%%zmm1			\n\t		vmovaps	     %%zmm8	,%%zmm9 			\n\t"\
		"vmovaps	     %%zmm2	,%%zmm3			\n\t		vmovaps	     %%zmm10,%%zmm11			\n\t"\
		"vmulpd	      (%%rcx),%%zmm0,%%zmm0		\n\t		vmulpd	     (%%r12),%%zmm8 ,%%zmm8 	\n\t"\
		"vmulpd	  0x40(%%rcx),%%zmm1,%%zmm1		\n\t		vmulpd	 0x40(%%r12),%%zmm9 ,%%zmm9 	\n\t"\
	"vfnmadd231pd 0x40(%%rcx),%%zmm2,%%zmm0		\n\t	vfnmadd231pd 0x40(%%r12),%%zmm10,%%zmm8 	\n\t"\
	" vfmadd231pd     (%%rcx),%%zmm3,%%zmm1		\n\t	 vfmadd231pd     (%%r12),%%zmm11,%%zmm9 	\n\t"\
		"vmovaps	    (%%rbx)	,%%zmm2			\n\t		vmovaps	    (%%r11)	,%%zmm10			\n\t"\
		"vmovaps	0x40(%%rbx)	,%%zmm3			\n\t		vmovaps	0x40(%%r11)	,%%zmm11			\n\t"\
		"vmovaps	     %%zmm2	,%%zmm4			\n\t		vmovaps	     %%zmm10,%%zmm12			\n\t"\
		"vmovaps	     %%zmm3	,%%zmm5			\n\t		vmovaps	     %%zmm11,%%zmm13			\n\t"\
		"vmulpd	      (%%rdx),%%zmm2,%%zmm2		\n\t		vmulpd	     (%%r13),%%zmm10,%%zmm10	\n\t"\
		"vmulpd	  0x40(%%rdx),%%zmm4,%%zmm4		\n\t		vmulpd	 0x40(%%r13),%%zmm12,%%zmm12	\n\t"\
	"vfnmadd231pd 0x40(%%rdx),%%zmm3,%%zmm2		\n\t	vfnmadd231pd 0x40(%%r13),%%zmm11,%%zmm10	\n\t"\
	" vfmadd231pd     (%%rdx),%%zmm5,%%zmm4		\n\t	 vfmadd231pd     (%%r13),%%zmm13,%%zmm12	\n\t"\
		"vmovaps	%%zmm2		,%%zmm3			\n\t		vmovaps	%%zmm10		,%%zmm11			\n\t"\
		"vmovaps	%%zmm4		,%%zmm5			\n\t		vmovaps	%%zmm12		,%%zmm13			\n\t"\
	"vmovaps	(%%rsi),%%zmm6	\n\t"/* 1.0 */\
		"vfmadd132pd	%%zmm6,%%zmm0,%%zmm2	\n\t		vfmadd132pd	%%zmm6,%%zmm8 ,%%zmm10	\n\t"\
		"vfmsub132pd	%%zmm6,%%zmm3,%%zmm0	\n\t		vfmsub132pd	%%zmm6,%%zmm11,%%zmm8 	\n\t"\
		"vfmadd132pd	%%zmm6,%%zmm1,%%zmm4	\n\t		vfmadd132pd	%%zmm6,%%zmm9 ,%%zmm12	\n\t"\
		"vfmsub132pd	%%zmm6,%%zmm5,%%zmm1	\n\t		vfmsub132pd	%%zmm6,%%zmm13,%%zmm9 	\n\t"\
		"vmovaps	%%zmm2		,    (%%rax)	\n\t		vmovaps	%%zmm10		,    (%%r10)		\n\t"\
		"vmovaps	%%zmm4		,0x40(%%rax)	\n\t		vmovaps	%%zmm12		,0x40(%%r10)		\n\t"\
		"vmovaps	%%zmm0		,    (%%rbx)	\n\t		vmovaps	%%zmm8 		,    (%%r11)		\n\t"\
		"vmovaps	%%zmm1		,0x40(%%rbx)	\n\t		vmovaps	%%zmm9 		,0x40(%%r11)		\n\t"\
	/* combine to get 2 length-4 output subtransforms... */\
		"movq		%[add0]	,%%rax				\n\t		movq		%[add4]	,%%r10				\n\t"\
		"movq		%[add2]	,%%rbx				\n\t		movq		%[add6]	,%%r11				\n\t"\
		"vmovaps	    (%%rax)	,%%zmm0			\n\t		vmovaps	    (%%r10)	,%%zmm8 			\n\t"\
		"vmovaps	0x40(%%rax)	,%%zmm1			\n\t		vmovaps	0x40(%%r10)	,%%zmm9 			\n\t"\
		"vmovaps	%%zmm0		,%%zmm4			\n\t		vmovaps	%%zmm8 		,%%zmm12			\n\t"\
		"vmovaps	%%zmm1		,%%zmm5			\n\t		vmovaps	%%zmm9 		,%%zmm13			\n\t"\
		"vfnmadd231pd	    (%%rbx),%%zmm6,%%zmm4	\n\t	vfnmadd231pd	0x40(%%r11),%%zmm6,%%zmm8 	\n\t"\
		"vfnmadd231pd	0x40(%%rbx),%%zmm6,%%zmm5	\n\t	vfnmadd231pd	    (%%r11),%%zmm6,%%zmm13	\n\t"\
		"vfmadd231pd	    (%%rbx),%%zmm6,%%zmm0	\n\t	vfmadd231pd		    (%%r11),%%zmm6,%%zmm9 	\n\t"\
		"vfmadd231pd	0x40(%%rbx),%%zmm6,%%zmm1	\n\t	vfmadd231pd	    0x40(%%r11),%%zmm6,%%zmm12	\n\t"\
		"vmovaps	%%zmm4		,    (%%rbx)	\n\t		vmovaps	%%zmm8 		,    (%%r10)	\n\t"\
		"vmovaps	%%zmm5		,0x40(%%rbx)	\n\t		vmovaps	%%zmm13		,0x40(%%r11)	\n\t"\
		"vmovaps	%%zmm0		,    (%%rax)	\n\t		vmovaps	%%zmm9 		,0x40(%%r10)	\n\t"\
		"vmovaps	%%zmm1		,0x40(%%rax)	\n\t		vmovaps	%%zmm12		,    (%%r11)	\n\t"\
		"movq		%[add1]	,%%rcx				\n\t		movq		%[add5]	,%%r12				\n\t"\
		"movq		%[add3]	,%%rdx				\n\t		movq		%[add7]	,%%r13				\n\t"\
		"vmovaps	    (%%rcx)	,%%zmm2			\n\t		vmovaps	    (%%r12)	,%%zmm10			\n\t"\
		"vmovaps	0x40(%%rcx)	,%%zmm3			\n\t		vmovaps	0x40(%%r12)	,%%zmm11			\n\t"\
		"vmovaps	%%zmm2		,%%zmm6			\n\t		vmovaps	%%zmm10		,%%zmm14			\n\t"\
		"vmovaps	%%zmm3		,%%zmm7			\n\t		vmovaps	%%zmm11		,%%zmm15			\n\t"\
	"vmovaps	(%%rsi),%%zmm8	\n\t"/* 1.0 */\
		"vfmadd231pd	    (%%rdx),%%zmm8,%%zmm2	\n\t	vfnmadd231pd	0x40(%%r13),%%zmm8,%%zmm10	\n\t"\
		"vfnmadd231pd	    (%%rdx),%%zmm8,%%zmm6	\n\t	vfmadd231pd		0x40(%%r13),%%zmm8,%%zmm14	\n\t"\
		"vfnmadd231pd	0x40(%%rdx),%%zmm8,%%zmm7	\n\t	vfnmadd231pd	    (%%r13),%%zmm8,%%zmm15	\n\t"\
		"vfmadd231pd	0x40(%%rdx),%%zmm8,%%zmm3	\n\t	vfmadd231pd	    	(%%r13),%%zmm8,%%zmm11	\n\t"\
		"vfmsub132pd	%%zmm8,%%zmm2,%%zmm0	\n\t		vmovaps	%%zmm12		,    (%%r13)	\n\t"\
		"vfmsub132pd	%%zmm8,%%zmm3,%%zmm1	\n\t		vmovaps	%%zmm13		,0x40(%%r13)	\n\t"\
		"vfmsub132pd	%%zmm8,%%zmm7,%%zmm4	\n\t		vmovaps	%%zmm10		,%%zmm13			\n\t"\
		"vfmsub132pd	%%zmm8,%%zmm6,%%zmm5	\n\t		vfmsub132pd	%%zmm8,%%zmm11,%%zmm10	\n\t"\
		"vfmadd231pd	    (%%rax),%%zmm8,%%zmm2	\n\t	vfmadd132pd	%%zmm8,%%zmm11,%%zmm13	\n\t"\
		"vfmadd231pd	0x40(%%rax),%%zmm8,%%zmm3	\n\t	vmulpd	-0x40(%%rsi),%%zmm10,%%zmm10	\n\t"/* *= isrt2 */\
		"vfmadd231pd	    (%%rbx),%%zmm8,%%zmm7	\n\t	vmulpd	-0x40(%%rsi),%%zmm13,%%zmm13	\n\t"\
		"vfmadd231pd	0x40(%%rbx),%%zmm8,%%zmm6	\n\t	vmovaps	0x40(%%r13)	,%%zmm11			\n\t"\
		"vmovaps	%%zmm0		,    (%%rcx)	\n\t		vmovaps	%%zmm15		,%%zmm12			\n\t"\
		"vmovaps	%%zmm1		,0x40(%%rcx)	\n\t		vfmadd132pd	%%zmm8,%%zmm14,%%zmm12	\n\t"\
		"vmovaps	%%zmm4		,    (%%rbx)	\n\t		vfmsub132pd	%%zmm8,%%zmm14,%%zmm15	\n\t"\
		"vmovaps	%%zmm5		,0x40(%%rdx)	\n\t		vmulpd	-0x40(%%rsi),%%zmm12,%%zmm12	\n\t"\
		"vmovaps	%%zmm2		,    (%%rax)	\n\t		vmulpd	-0x40(%%rsi),%%zmm15,%%zmm15	\n\t"\
		"vmovaps	%%zmm3		,0x40(%%rax)	\n\t		vmovaps		(%%r13)	,%%zmm14			\n\t"\
												"vmovaps	%%zmm8,%%zmm0	\n\t	vmovaps	(%%r10),%%zmm8	\n\t"/* Move 1.0 into zmm0 and restore earlier spill of zmm8 */\
		"vmovaps	%%zmm7		,    (%%rdx)	\n\t		vfmsub132pd	%%zmm0,%%zmm13,%%zmm9 	\n\t"\
		"vmovaps	%%zmm6		,0x40(%%rbx)	\n\t		vfmsub132pd	%%zmm0,%%zmm12,%%zmm14	\n\t"\
		"													vfmsub132pd	%%zmm0,%%zmm15,%%zmm11	\n\t"\
		"													vfmsub132pd	%%zmm0,%%zmm10,%%zmm8 	\n\t"\
		"													vfmadd231pd	    (%%r10),%%zmm0,%%zmm10	\n\t"\
		"													vfmadd231pd	0x40(%%r10),%%zmm0,%%zmm13	\n\t"\
		"													vfmadd231pd	    (%%r11),%%zmm0,%%zmm12	\n\t"\
		"													vfmadd231pd	0x40(%%r11),%%zmm0,%%zmm15	\n\t"\
		"													vmovaps	%%zmm8 		,    (%%r12)	\n\t"\
		"													vmovaps	%%zmm9 		,0x40(%%r12)	\n\t"\
		"													vmovaps	%%zmm14		,    (%%r11)	\n\t"\
		"													vmovaps	%%zmm11		,0x40(%%r11)	\n\t"\
		"													vmovaps	%%zmm10		,    (%%r10)	\n\t"\
		"													vmovaps	%%zmm13		,0x40(%%r10)	\n\t"\
		"													vmovaps	%%zmm12		,    (%%r13)	\n\t"\
		"													vmovaps	%%zmm15		,0x40(%%r13)	\n\t"\
		:					/* outputs: none */\
		: [add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[add1] "m" (Xadd1)\
		 ,[add2] "m" (Xadd2)\
		 ,[add3] "m" (Xadd3)\
		 ,[add4] "m" (Xadd4)\
		 ,[add5] "m" (Xadd5)\
		 ,[add6] "m" (Xadd6)\
		 ,[add7] "m" (Xadd7)\
		 ,[isrt2] "m" (Xisrt2)\
		 ,[c1] "m" (Xc1)\
		 ,[c2] "m" (Xc2)\
		 ,[c3] "m" (Xc3)\
		 ,[c4] "m" (Xc4)\
		 ,[c5] "m" (Xc5)\
		 ,[c6] "m" (Xc6)\
		 ,[c7] "m" (Xc7)\
		: "cc","memory","rax","rbx","rcx","rdx","rsi","r10","r11","r12","r13","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
	);\
	}

	// In our FMA-ized version, we replace the original mix of [35 ADD, 33 SUB, 32 MUL] with [4 SUB, 62 FMA, 18 MUL].
	//
	// FMA version assumes 1.0,2.0 in slots just above isrt2, i.e. vec_dbl *one = isrt2+1, *two = isrt2+2, or +0x40,0x80 in terms of byte offsets:
	#define SSE2_RADIX8_DIT_TWIDDLE(Xadd0,Xadd1,Xadd2,Xadd3,Xadd4,Xadd5,Xadd6,Xadd7,Xisrt2,Xc1,Xc2,Xc3,Xc4,Xc5,Xc6,Xc7)\
	{\
	__asm__ volatile (\
	"movq	%[isrt2],%%rsi	\n\t"\
	"leaq	0x40(%%rsi),%%rsi	\n\t"/* 1.0 */\
		/*** 2nd of 2 length-4 subtransforms gets done first: ***/\
		/* SSE2_RADIX4_DIT_0TWIDDLE_2NDOFTWO(add4-7):	SSE2_RADIX4_DIT_0TWIDDLE(add0-3): */\
		"movq		%[add4]		,%%rax			\n\t		movq		%[add0]		,%%r10		\n\t"\
		"movq		%[add5]		,%%rbx			\n\t		movq		%[add1]		,%%r11		\n\t"\
		"vmovaps		(%%rax)	,%%zmm0			\n\t		vmovaps			(%%r10)	,%%zmm8 	\n\t"\
		"vmovaps	0x40(%%rax)	,%%zmm1			\n\t		vmovaps		0x40(%%r10)	,%%zmm9 	\n\t"\
		"vmovaps		(%%rbx)	,%%zmm2			\n\t		vmovaps			(%%r11)	,%%zmm10	\n\t"\
		"vmovaps	0x40(%%rbx)	,%%zmm3			\n\t		vmovaps		0x40(%%r11)	,%%zmm11	\n\t"\
	"vmovaps	(%%rsi),%%zmm4 	\n\t"/* 1.0 */\
		"vfmsub132pd	%%zmm4 ,%%zmm2,%%zmm0	\n\t		vfmsub132pd	%%zmm4 ,%%zmm10,%%zmm8	\n\t"\
		"vfmsub132pd	%%zmm4 ,%%zmm3,%%zmm1	\n\t		vfmsub132pd	%%zmm4 ,%%zmm11,%%zmm9	\n\t"\
		"movq		%[add6]		,%%rcx			\n\t		movq		%[add2]		,%%r12		\n\t"\
		"movq		%[add7]		,%%rdx			\n\t		movq		%[add3]		,%%r13		\n\t"\
		"vmovaps			(%%rcx)	,%%zmm4		\n\t		vmovaps			(%%r12)	,%%zmm12	\n\t"\
		"vmovaps		0x40(%%rcx)	,%%zmm5		\n\t		vmovaps		0x40(%%r12)	,%%zmm13	\n\t"\
		"vmovaps			(%%rdx)	,%%zmm6		\n\t		vmovaps			(%%r13)	,%%zmm14	\n\t"/* t9  in zmm6, needed below in RHS! */\
		"vmovaps		0x40(%%rdx)	,%%zmm7		\n\t		vmovaps		0x40(%%r13)	,%%zmm15	\n\t"/* t10 in zmm7, needed below in RHS! */\
		"vfmsub132pd	(%%rsi),%%zmm6,%%zmm4	\n\t		vfmsub132pd	(%%rsi),%%zmm14,%%zmm12	\n\t"\
		"vfmsub132pd	(%%rsi),%%zmm7,%%zmm5	\n\t		vfmsub132pd	(%%rsi),%%zmm15,%%zmm13	\n\t"\
	"vmovaps	%%zmm15,(%%r12)	\n\t"/* spill zmm15 to make room for 2.0 */"	vmovaps	0x40(%%rsi),%%zmm15	\n\t"/* 2.0 */\
	"vfmadd132pd	%%zmm15,%%zmm0,%%zmm2		\n\t	vfmadd132pd	%%zmm15,%%zmm8 ,%%zmm10		\n\t"\
	"vfmadd132pd	%%zmm15,%%zmm1,%%zmm3		\n\t	vfmadd132pd	%%zmm15,%%zmm9 ,%%zmm11		\n\t"\
	"vfmadd132pd	%%zmm15,%%zmm4,%%zmm6		\n\t	vfmadd132pd	%%zmm15,%%zmm12,%%zmm14		\n\t"\
	"vfmadd132pd	%%zmm15,%%zmm5,%%zmm7		\n\t	vfmadd132pd	(%%r12),%%zmm13,%%zmm15		\n\t"\
		"vmovaps		%%zmm6	,    (%%rcx)	\n\t		vmovaps		%%zmm14		,    (%%r12)\n\t"\
		"vmovaps		%%zmm7	,0x40(%%rcx)	\n\t		vmovaps		%%zmm15		,0x40(%%r12)\n\t"\
	"vmovaps	%%zmm0 ,(%%r10)	\n\t"/* spill zmm0  to make room for 1.0 */"	vmovaps	(%%rsi),%%zmm0 	\n\t"/* 1.0 */\
		"vmovaps		%%zmm4	,    (%%rdx)	\n\t		vfmadd132pd	%%zmm0 ,%%zmm10,%%zmm14	\n\t"/* zmm14 <- ~t1 */\
		"vmovaps		%%zmm5	,0x40(%%rdx)	\n\t		vfmadd132pd	%%zmm0 ,%%zmm11,%%zmm15	\n\t"/* zmm15 <- ~t2 */\
		"vfmadd132pd	%%zmm0 ,%%zmm2,%%zmm6	\n\t		vsubpd		(%%r12),%%zmm10,%%zmm10	\n\t"/* zmm10 <- ~t5 */\
		"vfmadd132pd	%%zmm0 ,%%zmm3,%%zmm7	\n\t		vsubpd	0x40(%%r12),%%zmm11,%%zmm11	\n\t"/* zmm11 <- ~t6 */\
		"vsubpd		(%%rcx)	,%%zmm2,%%zmm2		\n\t		vfmadd132pd	%%zmm0 ,%%zmm6 ,%%zmm14	\n\t"/* t1+t9 */\
		"vsubpd	0x40(%%rcx)	,%%zmm3,%%zmm3		\n\t		vfmadd132pd	%%zmm0 ,%%zmm7 ,%%zmm15	\n\t"/* t2+t10*/\
	"vmovaps	%%zmm10,(%%rax)	\n\t"/* spill zmm10 to make room for 1.0 */"	vmovaps	%%zmm0 ,%%zmm10	\n\t"/* 1.0 */\
	"vmovaps	(%%r10),%%zmm0	\n\t"/* restore zmm0 */\
		"vmovaps		%%zmm2	,    (%%rcx)	\n\t"\
		"vmovaps		%%zmm3	,0x40(%%rcx)	\n\t"\
		"vmovaps		%%zmm4	,%%zmm2			\n\t		vmovaps		%%zmm14		,    (%%r10)\n\t"/* a[j1   ], DONE. */\
		"vmovaps		%%zmm5	,%%zmm3			\n\t		vmovaps		%%zmm15		,0x40(%%r10)\n\t"/* a[j2   ], DONE. */\
		"vfmadd132pd	%%zmm10,%%zmm0,%%zmm5	\n\t	vfnmadd231pd	0x40(%%rsi),%%zmm6 ,%%zmm14		\n\t"/* t1-t9  = [t1+t9 ] - 2*t9  */\
		"vfmsub132pd	%%zmm10,%%zmm3,%%zmm0	\n\t	vfnmadd231pd	0x40(%%rsi),%%zmm7 ,%%zmm15		\n\t"/* t2-t10 = [t2+t10] - 2*t10 */\
		"vfmadd132pd	%%zmm10,%%zmm1,%%zmm4	\n\t		vmovaps		%%zmm12		,%%zmm6		\n\t"/* zmm6<- copy of t7 */\
		"vfmsub132pd	%%zmm10,%%zmm2,%%zmm1	\n\t		vmovaps		%%zmm13		,%%zmm7		\n\t"/* zmm7<- copy of t8 */\
		"vmovaps		%%zmm5	,%%zmm2			\n\t		vfmadd132pd	%%zmm10,%%zmm8,%%zmm13	\n\t"/* zmm13<- ~t3 */\
		"vmovaps		%%zmm1	,%%zmm3			\n\t		vfmsub132pd	%%zmm10,%%zmm7,%%zmm8	\n\t"/* zmm8 <- ~t7 */\
		"vfmadd132pd	%%zmm10,%%zmm1,%%zmm5	\n\t		vfmadd132pd	%%zmm10,%%zmm9,%%zmm12	\n\t"/* zmm12<- ~t8 */\
		"vfmsub132pd	%%zmm10,%%zmm3,%%zmm2	\n\t		vfmsub132pd	%%zmm10,%%zmm6,%%zmm9 	\n\t"/* zmm9 <- ~t4 */\
		"vmovaps	-0x40(%%rsi),%%zmm1	\n\t"/* isrt2 */	/* Combine Outputs 0,4 of two half-transforms, as these are ready: */\
		"vmulpd	%%zmm1		,%%zmm5,%%zmm5		\n\t		movq		%[c4]		,%%rdi		\n\t"\
		"vmulpd	%%zmm1		,%%zmm2,%%zmm2		\n\t		vmovaps		%%zmm14,%%zmm6			\n\t"\
		"vmovaps		%%zmm0	,%%zmm3			\n\t		vmovaps		%%zmm15,%%zmm7			\n\t"\
		"vmovaps		%%zmm5	,    (%%rbx)	\n\t		vmulpd	0x40(%%rdi)	,%%zmm6	,%%zmm6 \n\t"\
		"vmovaps		%%zmm2	,0x40(%%rbx)	\n\t		vmulpd	0x40(%%rdi)	,%%zmm7	,%%zmm7 \n\t"\
		"vmovaps		%%zmm4	,%%zmm5			\n\t"\
		"vfmadd132pd	%%zmm10,%%zmm4,%%zmm0	\n\t"\
		"vfmsub132pd	%%zmm10,%%zmm5,%%zmm3	\n\t	vfmadd132pd		(%%rdi)	,%%zmm7 ,%%zmm14\n\t"\
		"vmulpd	%%zmm1		,%%zmm0,%%zmm0		\n\t	vfmsub132pd		(%%rdi)	,%%zmm6 ,%%zmm15\n\t"\
	"vmovaps	(%%rax),%%zmm10	\n\t"/* restore zmm10 */\
		"vmulpd	%%zmm1		,%%zmm3,%%zmm3		\n\t		vmovaps		%%zmm14,    (%%rax)		\n\t"\
		"vmovaps		%%zmm0	,    (%%rdx)	\n\t		vmovaps		%%zmm15,0x40(%%rax)		\n\t"\
		"vmovaps		%%zmm3	,0x40(%%rdx)	\n\t"\
	/* Now combine the two half-transforms & store outputs back into original array slots. */\
	/* add0-7 in r10,11,12,13,ax,bx,cx,dx; SIMD Registers 0-7,14-15 FREE: */\
		/* Outputs 1,5: Use zmm 4,5,9,13,14,15			Outputs 2,6: Use zmm 0,1,2,3,6,7,10,11 : */\
		"movq		%[c2],%%rdi		\n\t"/* c6 = c2+2, c1 = c2+4, c5 = c2+6 */\
		"vmovaps		    (%%rbx)	,%%zmm4 	\n\t		vmovaps		    (%%rcx)	,%%zmm2		\n\t"\
		"vmovaps		0x40(%%rbx)	,%%zmm5 	\n\t		vmovaps		0x40(%%rcx)	,%%zmm3		\n\t"\
	"vmovaps	(%%rsi),%%zmm0 	\n\t"/* 1.0 */\
		"vmovaps		%%zmm13	,%%zmm14		\n\t		vmovaps		%%zmm10		,%%zmm6		\n\t"\
		"vmovaps		%%zmm9 	,%%zmm15		\n\t		vmovaps		%%zmm11		,%%zmm7		\n\t"\
		"vfmadd132pd	%%zmm0 ,%%zmm4,%%zmm13	\n\t		vfmadd132pd	%%zmm0 ,%%zmm3,%%zmm10	\n\t"\
		"vfmsub132pd	%%zmm0 ,%%zmm5,%%zmm9 	\n\t		vfmsub132pd	%%zmm0 ,%%zmm2,%%zmm11	\n\t"\
		"vfmsub132pd	%%zmm0 ,%%zmm4,%%zmm14	\n\t"\
		"vfmadd132pd	%%zmm0 ,%%zmm5,%%zmm15	\n\t"\
		"vmovaps		%%zmm14	,%%zmm4 		\n\t		vfmsub132pd	%%zmm0 ,%%zmm3,%%zmm6	\n\t"\
		"vmovaps		%%zmm15	,%%zmm5 		\n\t		vfmadd132pd	%%zmm0 ,%%zmm2,%%zmm7	\n\t"/* zmm2,3 free */\
		"vmovaps		%%zmm13	,%%zmm14		\n\t		vmovaps		%%zmm10		,%%zmm0 	\n\t"\
		"vmovaps		%%zmm9 	,%%zmm15		\n\t		vmovaps		%%zmm11		,%%zmm1		\n\t"\
		"vmulpd	0x100(%%rdi)	,%%zmm13,%%zmm13	\n\t		vmulpd		(%%rdi)	,%%zmm10,%%zmm10\n\t"\
		"vmulpd	0x100(%%rdi)	,%%zmm9 ,%%zmm9 	\n\t		vmulpd		(%%rdi)	,%%zmm11,%%zmm11\n\t"\
	"vfnmadd231pd	0x140(%%rdi),%%zmm14,%%zmm9 	\n\t	vfnmadd231pd	0x40(%%rdi)	,%%zmm0 ,%%zmm11\n\t"\
	" vfmadd231pd	0x140(%%rdi),%%zmm15,%%zmm13	\n\t	 vfmadd231pd	0x40(%%rdi)	,%%zmm1	,%%zmm10\n\t"\
		"vmovaps		%%zmm9 	,0x40(%%r11)	\n\t		vmovaps		0x80(%%rdi)	,%%zmm2		\n\t"\
		"vmovaps		%%zmm13	,    (%%r11)	\n\t		vmovaps		0xc0(%%rdi)	,%%zmm3		\n\t"\
		"vmovaps		%%zmm4 	,%%zmm13		\n\t		vmovaps		%%zmm11		,0x40(%%r12)\n\t"\
		"vmovaps		%%zmm5 	,%%zmm9 		\n\t		vmovaps		%%zmm10		,    (%%r12)\n\t"\
		"vmovaps		%%zmm13	,%%zmm14		\n\t		vmovaps		%%zmm6		,%%zmm0 	\n\t"\
		"vmovaps		%%zmm9 	,%%zmm15		\n\t		vmovaps		%%zmm7		,%%zmm1		\n\t"\
		"vmulpd	0x180(%%rdi)	,%%zmm13,%%zmm13	\n\t		vmulpd	%%zmm2	,%%zmm6	,%%zmm6 	\n\t"\
		"vmulpd	0x180(%%rdi)	,%%zmm9 ,%%zmm9 	\n\t		vmulpd	%%zmm2	,%%zmm7	,%%zmm7		\n\t"\
	"vfnmadd231pd	0x1c0(%%rdi),%%zmm14,%%zmm9 	\n\t	vfnmadd231pd	%%zmm3	,%%zmm0 ,%%zmm7	\n\t"\
	" vfmadd231pd	0x1c0(%%rdi),%%zmm15,%%zmm13	\n\t	 vfmadd231pd	%%zmm3	,%%zmm1	,%%zmm6	\n\t"\
		"vmovaps		%%zmm9 	,0x40(%%rbx)	\n\t		vmovaps		%%zmm7		,0x40(%%rcx)\n\t"\
		"vmovaps		%%zmm13	,    (%%rbx)	\n\t		vmovaps		%%zmm6		,    (%%rcx)\n\t"\
		/* Outputs 3,7: All SIMD regs except for zmm8,12 (containing one of the 2 complex data being butterflied) free: */\
	"vmovaps	(%%rsi),%%zmm10	\n\t"/* 1.0 */\
		"movq		%[c3]		,%%rdi			\n\t"/* c7 = c3+2 */\
		"vmovaps		    (%%rdx)	,%%zmm4		\n\t"\
		"vmovaps		0x40(%%rdx)	,%%zmm5		\n\t"\
		"vmovaps		%%zmm8 	,%%zmm14		\n\t"\
		"vmovaps		%%zmm12	,%%zmm15		\n\t"\
		"vfmsub132pd	%%zmm10,%%zmm5,%%zmm8	\n\t		vfmsub132pd	%%zmm10,%%zmm4,%%zmm12	\n\t"\
		"vmovaps		    (%%rdi)	,%%zmm0		\n\t"\
		"vmovaps		0x40(%%rdi)	,%%zmm1		\n\t"\
		"vfmadd132pd	%%zmm10,%%zmm5,%%zmm14	\n\t		vmovaps		0x80(%%rdi)	,%%zmm2		\n\t"\
		"vfmadd132pd	%%zmm10,%%zmm4,%%zmm15	\n\t		vmovaps		0xc0(%%rdi)	,%%zmm3		\n\t"\
		"vmovaps		%%zmm8 	,%%zmm9 		\n\t		vmovaps		%%zmm14		,%%zmm6 	\n\t"\
		"vmovaps		%%zmm12	,%%zmm13		\n\t		vmovaps		%%zmm15		,%%zmm7		\n\t"\
		"vmulpd	%%zmm0	,%%zmm8 ,%%zmm8 		\n\t		vmulpd	%%zmm2	,%%zmm14,%%zmm14	\n\t"\
		"vmulpd	%%zmm0	,%%zmm12,%%zmm12		\n\t		vmulpd	%%zmm2	,%%zmm15,%%zmm15	\n\t"\
	"vfnmadd231pd	%%zmm1,%%zmm9 ,%%zmm12		\n\t	vfnmadd231pd	%%zmm3,%%zmm6,%%zmm15	\n\t"\
	" vfmadd231pd	%%zmm1,%%zmm13,%%zmm8 		\n\t	 vfmadd231pd	%%zmm3,%%zmm7,%%zmm14	\n\t"\
		"vmovaps		%%zmm12	,0x40(%%r13)	\n\t		vmovaps		%%zmm15		,0x40(%%rdx)\n\t"\
		"vmovaps		%%zmm8 	,    (%%r13)	\n\t		vmovaps		%%zmm14		,    (%%rdx)\n\t"\
		:					/* outputs: none */\
		: [add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[add1] "m" (Xadd1)\
		 ,[add2] "m" (Xadd2)\
		 ,[add3] "m" (Xadd3)\
		 ,[add4] "m" (Xadd4)\
		 ,[add5] "m" (Xadd5)\
		 ,[add6] "m" (Xadd6)\
		 ,[add7] "m" (Xadd7)\
		 ,[isrt2] "m" (Xisrt2)\
		 ,[c1] "m" (Xc1)\
		 ,[c2] "m" (Xc2)\
		 ,[c3] "m" (Xc3)\
		 ,[c4] "m" (Xc4)\
		 ,[c5] "m" (Xc5)\
		 ,[c6] "m" (Xc6)\
		 ,[c7] "m" (Xc7)\
		: "cc","memory","rax","rbx","rcx","rdx","rsi","rdi","r10","r11","r12","r13","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
	);\
	}

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
	// In our FMA-ized version, we replace the original mix of [33 ADD, 33 SUB, 32 MUL] with [66 FMA, 18 MUL].
	//
	// FMA version assumes 1.0,2.0 in slots just above isrt2, i.e. vec_dbl *one = isrt2+1, *two = isrt2+2, or +0x20,0x40 in terms of byte offsets:
	#define SSE2_RADIX8_DIF_TWIDDLE(Xadd0,Xadd1,Xadd2,Xadd3,Xadd4,Xadd5,Xadd6,Xadd7,Xisrt2,Xc1,Xc2,Xc3,Xc4,Xc5,Xc6,Xc7)\
	{\
	__asm__ volatile (\
	"movq	%[isrt2],%%rsi	\n\t"\
	"leaq	0x20(%%rsi),%%rsi	\n\t"/* 1.0 */\
		"							movq	%[add1]	,%%r10	\n\t	movq	%[c1],%%r12				\n\t"\
		"							movq	%[add5]	,%%r11	\n\t	movq	%[c5],%%r13				\n\t"\
		"										\n\t		vmovaps	    (%%r10)	,%%ymm8 			\n\t"\
		"										\n\t		vmovaps	0x20(%%r10)	,%%ymm10			\n\t"\
		"										\n\t		vmovaps	     %%ymm8	,%%ymm9 			\n\t"\
		"movq		%[add0]	,%%rax				\n\t		vmovaps	     %%ymm10,%%ymm11			\n\t"\
		"movq		%[add4]	,%%rbx				\n\t		vmulpd	     (%%r12),%%ymm8 ,%%ymm8 	\n\t"\
		"movq		%[c4]	,%%rcx				\n\t		vmulpd	 0x20(%%r12),%%ymm9 ,%%ymm9 	\n\t"\
		"vmovaps	    (%%rax)	,%%ymm0			\n\t	vfnmadd231pd 0x20(%%r12),%%ymm10,%%ymm8 	\n\t"\
		"vmovaps	0x20(%%rax)	,%%ymm1			\n\t	 vfmadd231pd     (%%r12),%%ymm11,%%ymm9 	\n\t"\
		"vmovaps	     %%ymm0	,%%ymm6			\n\t		vmovaps	    (%%r11)	,%%ymm10			\n\t"\
		"vmovaps	     %%ymm1	,%%ymm7			\n\t		vmovaps	0x20(%%r11)	,%%ymm11			\n\t"\
		"vmovaps	    (%%rbx)	,%%ymm2			\n\t		vmovaps	     %%ymm10,%%ymm12			\n\t"\
		"vmovaps	0x20(%%rbx)	,%%ymm3			\n\t		vmovaps	     %%ymm11,%%ymm13			\n\t"\
		"vmovaps	     %%ymm2	,%%ymm4			\n\t		vmulpd	     (%%r13),%%ymm10,%%ymm10	\n\t"\
		"vmovaps	     %%ymm3	,%%ymm5			\n\t		vmulpd	 0x20(%%r13),%%ymm12,%%ymm12	\n\t"\
		"vmulpd	        (%%rcx)	,%%ymm2,%%ymm2	\n\t	vfnmadd231pd 0x20(%%r13),%%ymm11,%%ymm10	\n\t"\
		"vmulpd	        (%%rcx)	,%%ymm3,%%ymm3	\n\t	 vfmadd231pd     (%%r13),%%ymm13,%%ymm12	\n\t"\
	"vfnmadd231pd	0x20(%%rcx)	,%%ymm5,%%ymm2	\n\t		vmovaps	%%ymm10		,%%ymm11			\n\t"\
	" vfmadd231pd	0x20(%%rcx)	,%%ymm4,%%ymm3	\n\t		vmovaps	%%ymm12		,%%ymm13			\n\t"\
	"vmovaps	(%%rsi),%%ymm4	\n\t"/* 1.0 */\
		"vfmadd132pd	%%ymm4,%%ymm2,%%ymm0	\n\t		vfmadd132pd	%%ymm4,%%ymm8 ,%%ymm10	\n\t"\
		"vfmadd132pd	%%ymm4,%%ymm3,%%ymm1	\n\t		vfmsub132pd	%%ymm4,%%ymm11,%%ymm8 	\n\t"\
		"vfmsub132pd	%%ymm4,%%ymm2,%%ymm6	\n\t		vfmadd132pd	%%ymm4,%%ymm9 ,%%ymm12	\n\t"\
		"vfmsub132pd	%%ymm4,%%ymm3,%%ymm7	\n\t		vfmsub132pd	%%ymm4,%%ymm13,%%ymm9 	\n\t"\
		"vmovaps	%%ymm0		,    (%%rax)	\n\t		vmovaps	%%ymm10		,    (%%r10)		\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rax)	\n\t		vmovaps	%%ymm12		,0x20(%%r10)		\n\t"\
		"vmovaps	%%ymm6		,    (%%rbx)	\n\t		vmovaps	%%ymm8 		,    (%%r11)		\n\t"\
		"vmovaps	%%ymm7		,0x20(%%rbx)	\n\t		vmovaps	%%ymm9 		,0x20(%%r11)		\n\t"\
		"movq		%[add2]	,%%rax				\n\t		movq		%[add3]	,%%r10				\n\t"\
		"movq		%[add6]	,%%rbx				\n\t		movq		%[add7]	,%%r11				\n\t"\
		"movq		%[c2]	,%%rcx				\n\t		movq		%[c3]	,%%r12				\n\t"\
		"movq		%[c6]	,%%rdx				\n\t		movq		%[c7]	,%%r13				\n\t"\
		"vmovaps	    (%%rax)	,%%ymm0			\n\t		vmovaps	    (%%r10)	,%%ymm8 			\n\t"\
		"vmovaps	0x20(%%rax)	,%%ymm2			\n\t		vmovaps	0x20(%%r10)	,%%ymm10			\n\t"\
		"vmovaps	     %%ymm0	,%%ymm1			\n\t		vmovaps	     %%ymm8	,%%ymm9 			\n\t"\
		"vmovaps	     %%ymm2	,%%ymm3			\n\t		vmovaps	     %%ymm10,%%ymm11			\n\t"\
		"vmulpd	      (%%rcx),%%ymm0,%%ymm0		\n\t		vmulpd	     (%%r12),%%ymm8 ,%%ymm8 	\n\t"\
		"vmulpd	  0x20(%%rcx),%%ymm1,%%ymm1		\n\t		vmulpd	 0x20(%%r12),%%ymm9 ,%%ymm9 	\n\t"\
	"vfnmadd231pd 0x20(%%rcx),%%ymm2,%%ymm0		\n\t	vfnmadd231pd 0x20(%%r12),%%ymm10,%%ymm8 	\n\t"\
	" vfmadd231pd     (%%rcx),%%ymm3,%%ymm1		\n\t	 vfmadd231pd     (%%r12),%%ymm11,%%ymm9 	\n\t"\
		"vmovaps	    (%%rbx)	,%%ymm2			\n\t		vmovaps	    (%%r11)	,%%ymm10			\n\t"\
		"vmovaps	0x20(%%rbx)	,%%ymm3			\n\t		vmovaps	0x20(%%r11)	,%%ymm11			\n\t"\
		"vmovaps	     %%ymm2	,%%ymm4			\n\t		vmovaps	     %%ymm10,%%ymm12			\n\t"\
		"vmovaps	     %%ymm3	,%%ymm5			\n\t		vmovaps	     %%ymm11,%%ymm13			\n\t"\
		"vmulpd	      (%%rdx),%%ymm2,%%ymm2		\n\t		vmulpd	     (%%r13),%%ymm10,%%ymm10	\n\t"\
		"vmulpd	  0x20(%%rdx),%%ymm4,%%ymm4		\n\t		vmulpd	 0x20(%%r13),%%ymm12,%%ymm12	\n\t"\
	"vfnmadd231pd 0x20(%%rdx),%%ymm3,%%ymm2		\n\t	vfnmadd231pd 0x20(%%r13),%%ymm11,%%ymm10	\n\t"\
	" vfmadd231pd     (%%rdx),%%ymm5,%%ymm4		\n\t	 vfmadd231pd     (%%r13),%%ymm13,%%ymm12	\n\t"\
		"vmovaps	%%ymm2		,%%ymm3			\n\t		vmovaps	%%ymm10		,%%ymm11			\n\t"\
		"vmovaps	%%ymm4		,%%ymm5			\n\t		vmovaps	%%ymm12		,%%ymm13			\n\t"\
	"vmovaps	(%%rsi),%%ymm6	\n\t"/* 1.0 */\
		"vfmadd132pd	%%ymm6,%%ymm0,%%ymm2	\n\t		vfmadd132pd	%%ymm6,%%ymm8 ,%%ymm10	\n\t"\
		"vfmsub132pd	%%ymm6,%%ymm3,%%ymm0	\n\t		vfmsub132pd	%%ymm6,%%ymm11,%%ymm8 	\n\t"\
		"vfmadd132pd	%%ymm6,%%ymm1,%%ymm4	\n\t		vfmadd132pd	%%ymm6,%%ymm9 ,%%ymm12	\n\t"\
		"vfmsub132pd	%%ymm6,%%ymm5,%%ymm1	\n\t		vfmsub132pd	%%ymm6,%%ymm13,%%ymm9 	\n\t"\
		"vmovaps	%%ymm2		,    (%%rax)	\n\t		vmovaps	%%ymm10		,    (%%r10)		\n\t"\
		"vmovaps	%%ymm4		,0x20(%%rax)	\n\t		vmovaps	%%ymm12		,0x20(%%r10)		\n\t"\
		"vmovaps	%%ymm0		,    (%%rbx)	\n\t		vmovaps	%%ymm8 		,    (%%r11)		\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rbx)	\n\t		vmovaps	%%ymm9 		,0x20(%%r11)		\n\t"\
	/* combine to get 2 length-4 output subtransforms... */\
		"movq		%[add0]	,%%rax				\n\t		movq		%[add4]	,%%r10				\n\t"\
		"movq		%[add2]	,%%rbx				\n\t		movq		%[add6]	,%%r11				\n\t"\
		"vmovaps	    (%%rax)	,%%ymm0			\n\t		vmovaps	    (%%r10)	,%%ymm8 			\n\t"\
		"vmovaps	0x20(%%rax)	,%%ymm1			\n\t		vmovaps	0x20(%%r10)	,%%ymm9 			\n\t"\
		"vmovaps	%%ymm0		,%%ymm4			\n\t		vmovaps	%%ymm8 		,%%ymm12			\n\t"\
		"vmovaps	%%ymm1		,%%ymm5			\n\t		vmovaps	%%ymm9 		,%%ymm13			\n\t"\
		"vfnmadd231pd	    (%%rbx),%%ymm6,%%ymm4	\n\t	vfnmadd231pd	0x20(%%r11),%%ymm6,%%ymm8 	\n\t"\
		"vfnmadd231pd	0x20(%%rbx),%%ymm6,%%ymm5	\n\t	vfnmadd231pd	    (%%r11),%%ymm6,%%ymm13	\n\t"\
		"vfmadd231pd	    (%%rbx),%%ymm6,%%ymm0	\n\t	vfmadd231pd		    (%%r11),%%ymm6,%%ymm9 	\n\t"\
		"vfmadd231pd	0x20(%%rbx),%%ymm6,%%ymm1	\n\t	vfmadd231pd	    0x20(%%r11),%%ymm6,%%ymm12	\n\t"\
		"vmovaps	%%ymm4		,    (%%rbx)	\n\t		vmovaps	%%ymm8 		,    (%%r10)	\n\t"\
		"vmovaps	%%ymm5		,0x20(%%rbx)	\n\t		vmovaps	%%ymm13		,0x20(%%r11)	\n\t"\
		"vmovaps	%%ymm0		,    (%%rax)	\n\t		vmovaps	%%ymm9 		,0x20(%%r10)	\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rax)	\n\t		vmovaps	%%ymm12		,    (%%r11)	\n\t"\
		"movq		%[add1]	,%%rcx				\n\t		movq		%[add5]	,%%r12				\n\t"\
		"movq		%[add3]	,%%rdx				\n\t		movq		%[add7]	,%%r13				\n\t"\
		"vmovaps	    (%%rcx)	,%%ymm2			\n\t		vmovaps	    (%%r12)	,%%ymm10			\n\t"\
		"vmovaps	0x20(%%rcx)	,%%ymm3			\n\t		vmovaps	0x20(%%r12)	,%%ymm11			\n\t"\
		"vmovaps	%%ymm2		,%%ymm6			\n\t		vmovaps	%%ymm10		,%%ymm14			\n\t"\
		"vmovaps	%%ymm3		,%%ymm7			\n\t		vmovaps	%%ymm11		,%%ymm15			\n\t"\
	"vmovaps	(%%rsi),%%ymm8	\n\t"/* 1.0 */\
		"vfmadd231pd	    (%%rdx),%%ymm8,%%ymm2	\n\t	vfnmadd231pd	0x20(%%r13),%%ymm8,%%ymm10	\n\t"\
		"vfnmadd231pd	    (%%rdx),%%ymm8,%%ymm6	\n\t	vfmadd231pd		0x20(%%r13),%%ymm8,%%ymm14	\n\t"\
		"vfnmadd231pd	0x20(%%rdx),%%ymm8,%%ymm7	\n\t	vfnmadd231pd	    (%%r13),%%ymm8,%%ymm15	\n\t"\
		"vfmadd231pd	0x20(%%rdx),%%ymm8,%%ymm3	\n\t	vfmadd231pd	    	(%%r13),%%ymm8,%%ymm11	\n\t"\
		"vfmsub132pd	%%ymm8,%%ymm2,%%ymm0	\n\t		vmovaps	%%ymm12		,    (%%r13)	\n\t"\
		"vfmsub132pd	%%ymm8,%%ymm3,%%ymm1	\n\t		vmovaps	%%ymm13		,0x20(%%r13)	\n\t"\
		"vfmsub132pd	%%ymm8,%%ymm7,%%ymm4	\n\t		vmovaps	%%ymm10		,%%ymm13			\n\t"\
		"vfmsub132pd	%%ymm8,%%ymm6,%%ymm5	\n\t		vfmsub132pd	%%ymm8,%%ymm11,%%ymm10	\n\t"\
		"vfmadd231pd	    (%%rax),%%ymm8,%%ymm2	\n\t	vfmadd132pd	%%ymm8,%%ymm11,%%ymm13	\n\t"\
		"vfmadd231pd	0x20(%%rax),%%ymm8,%%ymm3	\n\t	vmulpd	-0x20(%%rsi),%%ymm10,%%ymm10	\n\t"/* *= isrt2 */\
		"vfmadd231pd	    (%%rbx),%%ymm8,%%ymm7	\n\t	vmulpd	-0x20(%%rsi),%%ymm13,%%ymm13	\n\t"\
		"vfmadd231pd	0x20(%%rbx),%%ymm8,%%ymm6	\n\t	vmovaps	0x20(%%r13)	,%%ymm11			\n\t"\
		"vmovaps	%%ymm0		,    (%%rcx)	\n\t		vmovaps	%%ymm15		,%%ymm12			\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rcx)	\n\t		vfmadd132pd	%%ymm8,%%ymm14,%%ymm12	\n\t"\
		"vmovaps	%%ymm4		,    (%%rbx)	\n\t		vfmsub132pd	%%ymm8,%%ymm14,%%ymm15	\n\t"\
		"vmovaps	%%ymm5		,0x20(%%rdx)	\n\t		vmulpd	-0x20(%%rsi),%%ymm12,%%ymm12	\n\t"\
		"vmovaps	%%ymm2		,    (%%rax)	\n\t		vmulpd	-0x20(%%rsi),%%ymm15,%%ymm15	\n\t"\
		"vmovaps	%%ymm3		,0x20(%%rax)	\n\t		vmovaps		(%%r13)	,%%ymm14			\n\t"\
												"vmovaps	%%ymm8,%%ymm0	\n\t	vmovaps	(%%r10),%%ymm8	\n\t"/* Move 1.0 into ymm0 and restore earlier spill of ymm8 */\
		"vmovaps	%%ymm7		,    (%%rdx)	\n\t		vfmsub132pd	%%ymm0,%%ymm13,%%ymm9 	\n\t"\
		"vmovaps	%%ymm6		,0x20(%%rbx)	\n\t		vfmsub132pd	%%ymm0,%%ymm12,%%ymm14	\n\t"\
		"													vfmsub132pd	%%ymm0,%%ymm15,%%ymm11	\n\t"\
		"													vfmsub132pd	%%ymm0,%%ymm10,%%ymm8 	\n\t"\
		"													vfmadd231pd	    (%%r10),%%ymm0,%%ymm10	\n\t"\
		"													vfmadd231pd	0x20(%%r10),%%ymm0,%%ymm13	\n\t"\
		"													vfmadd231pd	    (%%r11),%%ymm0,%%ymm12	\n\t"\
		"													vfmadd231pd	0x20(%%r11),%%ymm0,%%ymm15	\n\t"\
		"													vmovaps	%%ymm8 		,    (%%r12)	\n\t"\
		"													vmovaps	%%ymm9 		,0x20(%%r12)	\n\t"\
		"													vmovaps	%%ymm14		,    (%%r11)	\n\t"\
		"													vmovaps	%%ymm11		,0x20(%%r11)	\n\t"\
		"													vmovaps	%%ymm10		,    (%%r10)	\n\t"\
		"													vmovaps	%%ymm13		,0x20(%%r10)	\n\t"\
		"													vmovaps	%%ymm12		,    (%%r13)	\n\t"\
		"													vmovaps	%%ymm15		,0x20(%%r13)	\n\t"\
		:					/* outputs: none */\
		: [add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[add1] "m" (Xadd1)\
		 ,[add2] "m" (Xadd2)\
		 ,[add3] "m" (Xadd3)\
		 ,[add4] "m" (Xadd4)\
		 ,[add5] "m" (Xadd5)\
		 ,[add6] "m" (Xadd6)\
		 ,[add7] "m" (Xadd7)\
		 ,[isrt2] "m" (Xisrt2)\
		 ,[c1] "m" (Xc1)\
		 ,[c2] "m" (Xc2)\
		 ,[c3] "m" (Xc3)\
		 ,[c4] "m" (Xc4)\
		 ,[c5] "m" (Xc5)\
		 ,[c6] "m" (Xc6)\
		 ,[c7] "m" (Xc7)\
		: "cc","memory","rax","rbx","rcx","rdx","rsi","r10","r11","r12","r13","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
	);\
	}

	// In our FMA-ized version, we replace the original mix of [35 ADD, 33 SUB, 32 MUL] with [4 SUB, 62 FMA, 18 MUL].
	//
	// FMA version assumes 1.0,2.0 in slots just above isrt2, i.e. vec_dbl *one = isrt2+1, *two = isrt2+2, or +0x20,0x40 in terms of byte offsets:
	#define SSE2_RADIX8_DIT_TWIDDLE(Xadd0,Xadd1,Xadd2,Xadd3,Xadd4,Xadd5,Xadd6,Xadd7,Xisrt2,Xc1,Xc2,Xc3,Xc4,Xc5,Xc6,Xc7)\
	{\
	__asm__ volatile (\
	"movq	%[isrt2],%%rsi	\n\t"\
	"leaq	0x20(%%rsi),%%rsi	\n\t"/* 1.0 */\
		/*** 2nd of 2 length-4 subtransforms gets done first: ***/\
		/* SSE2_RADIX4_DIT_0TWIDDLE_2NDOFTWO(add4-7):	SSE2_RADIX4_DIT_0TWIDDLE(add0-3): */\
		"movq		%[add4]		,%%rax			\n\t		movq		%[add0]		,%%r10		\n\t"\
		"movq		%[add5]		,%%rbx			\n\t		movq		%[add1]		,%%r11		\n\t"\
		"vmovaps		(%%rax)	,%%ymm0			\n\t		vmovaps			(%%r10)	,%%ymm8 	\n\t"\
		"vmovaps	0x20(%%rax)	,%%ymm1			\n\t		vmovaps		0x20(%%r10)	,%%ymm9 	\n\t"\
		"vmovaps		(%%rbx)	,%%ymm2			\n\t		vmovaps			(%%r11)	,%%ymm10	\n\t"\
		"vmovaps	0x20(%%rbx)	,%%ymm3			\n\t		vmovaps		0x20(%%r11)	,%%ymm11	\n\t"\
	"vmovaps	(%%rsi),%%ymm4 	\n\t"/* 1.0 */\
		"vfmsub132pd	%%ymm4 ,%%ymm2,%%ymm0	\n\t		vfmsub132pd	%%ymm4 ,%%ymm10,%%ymm8	\n\t"\
		"vfmsub132pd	%%ymm4 ,%%ymm3,%%ymm1	\n\t		vfmsub132pd	%%ymm4 ,%%ymm11,%%ymm9	\n\t"\
		"movq		%[add6]		,%%rcx			\n\t		movq		%[add2]		,%%r12		\n\t"\
		"movq		%[add7]		,%%rdx			\n\t		movq		%[add3]		,%%r13		\n\t"\
		"vmovaps			(%%rcx)	,%%ymm4		\n\t		vmovaps			(%%r12)	,%%ymm12	\n\t"\
		"vmovaps		0x20(%%rcx)	,%%ymm5		\n\t		vmovaps		0x20(%%r12)	,%%ymm13	\n\t"\
		"vmovaps			(%%rdx)	,%%ymm6		\n\t		vmovaps			(%%r13)	,%%ymm14	\n\t"/* t9  in ymm6, needed below in RHS! */\
		"vmovaps		0x20(%%rdx)	,%%ymm7		\n\t		vmovaps		0x20(%%r13)	,%%ymm15	\n\t"/* t10 in ymm7, needed below in RHS! */\
		"vfmsub132pd	(%%rsi),%%ymm6,%%ymm4	\n\t		vfmsub132pd	(%%rsi),%%ymm14,%%ymm12	\n\t"\
		"vfmsub132pd	(%%rsi),%%ymm7,%%ymm5	\n\t		vfmsub132pd	(%%rsi),%%ymm15,%%ymm13	\n\t"\
	"vmovaps	%%ymm15,(%%r12)	\n\t"/* spill ymm15 to make room for 2.0 */"	vmovaps	0x20(%%rsi),%%ymm15	\n\t"/* 2.0 */\
	"vfmadd132pd	%%ymm15,%%ymm0,%%ymm2		\n\t	vfmadd132pd	%%ymm15,%%ymm8 ,%%ymm10		\n\t"\
	"vfmadd132pd	%%ymm15,%%ymm1,%%ymm3		\n\t	vfmadd132pd	%%ymm15,%%ymm9 ,%%ymm11		\n\t"\
	"vfmadd132pd	%%ymm15,%%ymm4,%%ymm6		\n\t	vfmadd132pd	%%ymm15,%%ymm12,%%ymm14		\n\t"\
	"vfmadd132pd	%%ymm15,%%ymm5,%%ymm7		\n\t	vfmadd132pd	(%%r12),%%ymm13,%%ymm15		\n\t"\
		"vmovaps		%%ymm6	,    (%%rcx)	\n\t		vmovaps		%%ymm14		,    (%%r12)\n\t"\
		"vmovaps		%%ymm7	,0x20(%%rcx)	\n\t		vmovaps		%%ymm15		,0x20(%%r12)\n\t"\
	"vmovaps	%%ymm0 ,(%%r10)	\n\t"/* spill ymm0  to make room for 1.0 */"	vmovaps	(%%rsi),%%ymm0 	\n\t"/* 1.0 */\
		"vmovaps		%%ymm4	,    (%%rdx)	\n\t		vfmadd132pd	%%ymm0 ,%%ymm10,%%ymm14	\n\t"/* ymm14 <- ~t1 */\
		"vmovaps		%%ymm5	,0x20(%%rdx)	\n\t		vfmadd132pd	%%ymm0 ,%%ymm11,%%ymm15	\n\t"/* ymm15 <- ~t2 */\
		"vfmadd132pd	%%ymm0 ,%%ymm2,%%ymm6	\n\t		vsubpd		(%%r12),%%ymm10,%%ymm10	\n\t"/* ymm10 <- ~t5 */\
		"vfmadd132pd	%%ymm0 ,%%ymm3,%%ymm7	\n\t		vsubpd	0x20(%%r12),%%ymm11,%%ymm11	\n\t"/* ymm11 <- ~t6 */\
		"vsubpd		(%%rcx)	,%%ymm2,%%ymm2		\n\t		vfmadd132pd	%%ymm0 ,%%ymm6 ,%%ymm14	\n\t"/* t1+t9 */\
		"vsubpd	0x20(%%rcx)	,%%ymm3,%%ymm3		\n\t		vfmadd132pd	%%ymm0 ,%%ymm7 ,%%ymm15	\n\t"/* t2+t10*/\
	"vmovaps	%%ymm10,(%%rax)	\n\t"/* spill ymm10 to make room for 1.0 */"	vmovaps	%%ymm0 ,%%ymm10	\n\t"/* 1.0 */\
	"vmovaps	(%%r10),%%ymm0	\n\t"/* restore ymm0 */\
		"vmovaps		%%ymm2	,    (%%rcx)	\n\t"\
		"vmovaps		%%ymm3	,0x20(%%rcx)	\n\t"\
		"vmovaps		%%ymm4	,%%ymm2			\n\t		vmovaps		%%ymm14		,    (%%r10)\n\t"/* a[j1   ], DONE. */\
		"vmovaps		%%ymm5	,%%ymm3			\n\t		vmovaps		%%ymm15		,0x20(%%r10)\n\t"/* a[j2   ], DONE. */\
		"vfmadd132pd	%%ymm10,%%ymm0,%%ymm5	\n\t	vfnmadd231pd	0x20(%%rsi),%%ymm6 ,%%ymm14		\n\t"/* t1-t9  = [t1+t9 ] - 2*t9  */\
		"vfmsub132pd	%%ymm10,%%ymm3,%%ymm0	\n\t	vfnmadd231pd	0x20(%%rsi),%%ymm7 ,%%ymm15		\n\t"/* t2-t10 = [t2+t10] - 2*t10 */\
		"vfmadd132pd	%%ymm10,%%ymm1,%%ymm4	\n\t		vmovaps		%%ymm12		,%%ymm6		\n\t"/* ymm6<- copy of t7 */\
		"vfmsub132pd	%%ymm10,%%ymm2,%%ymm1	\n\t		vmovaps		%%ymm13		,%%ymm7		\n\t"/* ymm7<- copy of t8 */\
		"vmovaps		%%ymm5	,%%ymm2			\n\t		vfmadd132pd	%%ymm10,%%ymm8,%%ymm13	\n\t"/* ymm13<- ~t3 */\
		"vmovaps		%%ymm1	,%%ymm3			\n\t		vfmsub132pd	%%ymm10,%%ymm7,%%ymm8	\n\t"/* ymm8 <- ~t7 */\
		"vfmadd132pd	%%ymm10,%%ymm1,%%ymm5	\n\t		vfmadd132pd	%%ymm10,%%ymm9,%%ymm12	\n\t"/* ymm12<- ~t8 */\
		"vfmsub132pd	%%ymm10,%%ymm3,%%ymm2	\n\t		vfmsub132pd	%%ymm10,%%ymm6,%%ymm9 	\n\t"/* ymm9 <- ~t4 */\
		"vmovaps	-0x20(%%rsi),%%ymm1	\n\t"/* isrt2 */	/* Combine Outputs 0,4 of two half-transforms, as these are ready: */\
		"vmulpd	%%ymm1		,%%ymm5,%%ymm5		\n\t		movq		%[c4]		,%%rdi		\n\t"\
		"vmulpd	%%ymm1		,%%ymm2,%%ymm2		\n\t		vmovaps		%%ymm14,%%ymm6			\n\t"\
		"vmovaps		%%ymm0	,%%ymm3			\n\t		vmovaps		%%ymm15,%%ymm7			\n\t"\
		"vmovaps		%%ymm5	,    (%%rbx)	\n\t		vmulpd	0x20(%%rdi)	,%%ymm6	,%%ymm6 \n\t"\
		"vmovaps		%%ymm2	,0x20(%%rbx)	\n\t		vmulpd	0x20(%%rdi)	,%%ymm7	,%%ymm7 \n\t"\
		"vmovaps		%%ymm4	,%%ymm5			\n\t"\
		"vfmadd132pd	%%ymm10,%%ymm4,%%ymm0	\n\t"\
		"vfmsub132pd	%%ymm10,%%ymm5,%%ymm3	\n\t	vfmadd132pd		(%%rdi)	,%%ymm7 ,%%ymm14\n\t"\
		"vmulpd	%%ymm1		,%%ymm0,%%ymm0		\n\t	vfmsub132pd		(%%rdi)	,%%ymm6 ,%%ymm15\n\t"\
	"vmovaps	(%%rax),%%ymm10	\n\t"/* restore ymm10 */\
		"vmulpd	%%ymm1		,%%ymm3,%%ymm3		\n\t		vmovaps		%%ymm14,    (%%rax)		\n\t"\
		"vmovaps		%%ymm0	,    (%%rdx)	\n\t		vmovaps		%%ymm15,0x20(%%rax)		\n\t"\
		"vmovaps		%%ymm3	,0x20(%%rdx)	\n\t"\
	/* Now combine the two half-transforms & store outputs back into original array slots. */\
	/* add0-7 in r10,11,12,13,ax,bx,cx,dx; SIMD Registers 0-7,14-15 FREE: */\
		/* Outputs 1,5: Use ymm 4,5,9,13,14,15			Outputs 2,6: Use ymm 0,1,2,3,6,7,10,11 : */\
		"movq		%[c2],%%rdi		\n\t"/* c6 = c2+2, c1 = c2+4, c5 = c2+6 */\
		"vmovaps		    (%%rbx)	,%%ymm4 	\n\t		vmovaps		    (%%rcx)	,%%ymm2		\n\t"\
		"vmovaps		0x20(%%rbx)	,%%ymm5 	\n\t		vmovaps		0x20(%%rcx)	,%%ymm3		\n\t"\
	"vmovaps	(%%rsi),%%ymm0 	\n\t"/* 1.0 */\
		"vmovaps		%%ymm13	,%%ymm14		\n\t		vmovaps		%%ymm10		,%%ymm6		\n\t"\
		"vmovaps		%%ymm9 	,%%ymm15		\n\t		vmovaps		%%ymm11		,%%ymm7		\n\t"\
		"vfmadd132pd	%%ymm0 ,%%ymm4,%%ymm13	\n\t		vfmadd132pd	%%ymm0 ,%%ymm3,%%ymm10	\n\t"\
		"vfmsub132pd	%%ymm0 ,%%ymm5,%%ymm9 	\n\t		vfmsub132pd	%%ymm0 ,%%ymm2,%%ymm11	\n\t"\
		"vfmsub132pd	%%ymm0 ,%%ymm4,%%ymm14	\n\t"\
		"vfmadd132pd	%%ymm0 ,%%ymm5,%%ymm15	\n\t"\
		"vmovaps		%%ymm14	,%%ymm4 		\n\t		vfmsub132pd	%%ymm0 ,%%ymm3,%%ymm6	\n\t"\
		"vmovaps		%%ymm15	,%%ymm5 		\n\t		vfmadd132pd	%%ymm0 ,%%ymm2,%%ymm7	\n\t"/* ymm2,3 free */\
		"vmovaps		%%ymm13	,%%ymm14		\n\t		vmovaps		%%ymm10		,%%ymm0 	\n\t"\
		"vmovaps		%%ymm9 	,%%ymm15		\n\t		vmovaps		%%ymm11		,%%ymm1		\n\t"\
		"vmulpd	0x80(%%rdi)	,%%ymm13,%%ymm13	\n\t		vmulpd		(%%rdi)	,%%ymm10,%%ymm10\n\t"\
		"vmulpd	0x80(%%rdi)	,%%ymm9 ,%%ymm9 	\n\t		vmulpd		(%%rdi)	,%%ymm11,%%ymm11\n\t"\
	"vfnmadd231pd	0xa0(%%rdi),%%ymm14,%%ymm9 	\n\t	vfnmadd231pd	0x20(%%rdi)	,%%ymm0 ,%%ymm11\n\t"\
	" vfmadd231pd	0xa0(%%rdi),%%ymm15,%%ymm13	\n\t	 vfmadd231pd	0x20(%%rdi)	,%%ymm1	,%%ymm10\n\t"\
		"vmovaps		%%ymm9 	,0x20(%%r11)	\n\t		vmovaps		0x40(%%rdi)	,%%ymm2		\n\t"\
		"vmovaps		%%ymm13	,    (%%r11)	\n\t		vmovaps		0x60(%%rdi)	,%%ymm3		\n\t"\
		"vmovaps		%%ymm4 	,%%ymm13		\n\t		vmovaps		%%ymm11		,0x20(%%r12)\n\t"\
		"vmovaps		%%ymm5 	,%%ymm9 		\n\t		vmovaps		%%ymm10		,    (%%r12)\n\t"\
		"vmovaps		%%ymm13	,%%ymm14		\n\t		vmovaps		%%ymm6		,%%ymm0 	\n\t"\
		"vmovaps		%%ymm9 	,%%ymm15		\n\t		vmovaps		%%ymm7		,%%ymm1		\n\t"\
		"vmulpd	0xc0(%%rdi)	,%%ymm13,%%ymm13	\n\t		vmulpd	%%ymm2	,%%ymm6	,%%ymm6 	\n\t"\
		"vmulpd	0xc0(%%rdi)	,%%ymm9 ,%%ymm9 	\n\t		vmulpd	%%ymm2	,%%ymm7	,%%ymm7		\n\t"\
	"vfnmadd231pd	0xe0(%%rdi),%%ymm14,%%ymm9 	\n\t	vfnmadd231pd	%%ymm3	,%%ymm0 ,%%ymm7	\n\t"\
	" vfmadd231pd	0xe0(%%rdi),%%ymm15,%%ymm13	\n\t	 vfmadd231pd	%%ymm3	,%%ymm1	,%%ymm6	\n\t"\
		"vmovaps		%%ymm9 	,0x20(%%rbx)	\n\t		vmovaps		%%ymm7		,0x20(%%rcx)\n\t"\
		"vmovaps		%%ymm13	,    (%%rbx)	\n\t		vmovaps		%%ymm6		,    (%%rcx)\n\t"\
		/* Outputs 3,7: All SIMD regs except for ymm8,12 (containing one of the 2 complex data being butterflied) free: */\
	"vmovaps	(%%rsi),%%ymm10	\n\t"/* 1.0 */\
		"movq		%[c3]		,%%rdi			\n\t"/* c7 = c3+2 */\
		"vmovaps		    (%%rdx)	,%%ymm4		\n\t"\
		"vmovaps		0x20(%%rdx)	,%%ymm5		\n\t"\
		"vmovaps		%%ymm8 	,%%ymm14		\n\t"\
		"vmovaps		%%ymm12	,%%ymm15		\n\t"\
		"vfmsub132pd	%%ymm10,%%ymm5,%%ymm8	\n\t		vfmsub132pd	%%ymm10,%%ymm4,%%ymm12	\n\t"\
		"vmovaps		    (%%rdi)	,%%ymm0		\n\t"\
		"vmovaps		0x20(%%rdi)	,%%ymm1		\n\t"\
		"vfmadd132pd	%%ymm10,%%ymm5,%%ymm14	\n\t		vmovaps		0x40(%%rdi)	,%%ymm2		\n\t"\
		"vfmadd132pd	%%ymm10,%%ymm4,%%ymm15	\n\t		vmovaps		0x60(%%rdi)	,%%ymm3		\n\t"\
		"vmovaps		%%ymm8 	,%%ymm9 		\n\t		vmovaps		%%ymm14		,%%ymm6 	\n\t"\
		"vmovaps		%%ymm12	,%%ymm13		\n\t		vmovaps		%%ymm15		,%%ymm7		\n\t"\
		"vmulpd	%%ymm0	,%%ymm8 ,%%ymm8 		\n\t		vmulpd	%%ymm2	,%%ymm14,%%ymm14	\n\t"\
		"vmulpd	%%ymm0	,%%ymm12,%%ymm12		\n\t		vmulpd	%%ymm2	,%%ymm15,%%ymm15	\n\t"\
	"vfnmadd231pd	%%ymm1,%%ymm9 ,%%ymm12		\n\t	vfnmadd231pd	%%ymm3,%%ymm6,%%ymm15	\n\t"\
	" vfmadd231pd	%%ymm1,%%ymm13,%%ymm8 		\n\t	 vfmadd231pd	%%ymm3,%%ymm7,%%ymm14	\n\t"\
		"vmovaps		%%ymm12	,0x20(%%r13)	\n\t		vmovaps		%%ymm15		,0x20(%%rdx)\n\t"\
		"vmovaps		%%ymm8 	,    (%%r13)	\n\t		vmovaps		%%ymm14		,    (%%rdx)\n\t"\
		:					/* outputs: none */\
		: [add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[add1] "m" (Xadd1)\
		 ,[add2] "m" (Xadd2)\
		 ,[add3] "m" (Xadd3)\
		 ,[add4] "m" (Xadd4)\
		 ,[add5] "m" (Xadd5)\
		 ,[add6] "m" (Xadd6)\
		 ,[add7] "m" (Xadd7)\
		 ,[isrt2] "m" (Xisrt2)\
		 ,[c1] "m" (Xc1)\
		 ,[c2] "m" (Xc2)\
		 ,[c3] "m" (Xc3)\
		 ,[c4] "m" (Xc4)\
		 ,[c5] "m" (Xc5)\
		 ,[c6] "m" (Xc6)\
		 ,[c7] "m" (Xc7)\
		: "cc","memory","rax","rbx","rcx","rdx","rsi","rdi","r10","r11","r12","r13","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
	);\
	}

#elif defined(USE_AVX)	// AVX and AVX2 both use 256-bit registers

	#define SSE2_RADIX8_DIF_TWIDDLE(Xadd0,Xadd1,Xadd2,Xadd3,Xadd4,Xadd5,Xadd6,Xadd7,Xisrt2,Xc1,Xc2,Xc3,Xc4,Xc5,Xc6,Xc7)\
	{\
	__asm__ volatile (\
		"													movq		%[add1]	,%%r10				\n\t"\
		"													movq		%[add5]	,%%r11				\n\t"\
		"													movq		%[c1]	,%%r12				\n\t"\
		"													movq		%[c5]	,%%r13				\n\t"\
		"													vmovaps	    (%%r10)	,%%ymm8 			\n\t"\
		"movq		%[add0]	,%%rax				\n\t		vmovaps	0x20(%%r10)	,%%ymm10			\n\t"\
		"movq		%[add4]	,%%rbx				\n\t		vmovaps	    (%%r10)	,%%ymm9 			\n\t"\
		"movq		%[c4]	,%%rcx				\n\t		vmovaps	0x20(%%r10)	,%%ymm11			\n\t"\
		"vmovaps	    (%%rax)	,%%ymm0			\n\t		vmulpd	    (%%r12)	,%%ymm8 ,%%ymm8 	\n\t"\
		"vmovaps	0x20(%%rax)	,%%ymm1			\n\t		vmulpd	0x20(%%r12)	,%%ymm10,%%ymm10	\n\t"\
		"vmovaps	    (%%rax)	,%%ymm6			\n\t		vmulpd	0x20(%%r12)	,%%ymm9 ,%%ymm9 	\n\t"\
		"vmovaps	0x20(%%rax)	,%%ymm7			\n\t		vmulpd	    (%%r12)	,%%ymm11,%%ymm11	\n\t"\
		"vmovaps	    (%%rbx)	,%%ymm2			\n\t		vsubpd	%%ymm10		,%%ymm8 ,%%ymm8 	\n\t"\
		"vmovaps	0x20(%%rbx)	,%%ymm3			\n\t		vaddpd	%%ymm11		,%%ymm9 ,%%ymm9 	\n\t"\
		"vmovaps	    (%%rbx)	,%%ymm4			\n\t		vmovaps	    (%%r11)	,%%ymm10			\n\t"\
		"vmovaps	0x20(%%rbx)	,%%ymm5			\n\t		vmovaps	0x20(%%r11)	,%%ymm11			\n\t"\
		"vmulpd	    (%%rcx)	,%%ymm2,%%ymm2		\n\t		vmovaps	    (%%r11)	,%%ymm12			\n\t"\
		"vmulpd	    (%%rcx)	,%%ymm3,%%ymm3		\n\t		vmovaps	0x20(%%r11)	,%%ymm13			\n\t"\
		"vmulpd	0x20(%%rcx)	,%%ymm4,%%ymm4		\n\t		vmulpd	    (%%r13)	,%%ymm10,%%ymm10	\n\t"\
		"vmulpd	0x20(%%rcx)	,%%ymm5,%%ymm5		\n\t		vmulpd	0x20(%%r13)	,%%ymm11,%%ymm11	\n\t"\
		"vsubpd	%%ymm5		,%%ymm2,%%ymm2		\n\t		vmulpd	0x20(%%r13)	,%%ymm12,%%ymm12	\n\t"\
		"vaddpd	%%ymm4		,%%ymm3,%%ymm3		\n\t		vmulpd	    (%%r13)	,%%ymm13,%%ymm13	\n\t"\
		"vaddpd	%%ymm2		,%%ymm0,%%ymm0		\n\t		vsubpd	%%ymm11		,%%ymm10,%%ymm10	\n\t"\
		"vaddpd	%%ymm3		,%%ymm1,%%ymm1		\n\t		vaddpd	%%ymm13		,%%ymm12,%%ymm12	\n\t"\
		"vsubpd	%%ymm2		,%%ymm6,%%ymm6		\n\t		vmovaps	%%ymm10		,%%ymm11			\n\t"\
		"vsubpd	%%ymm3		,%%ymm7,%%ymm7		\n\t		vmovaps	%%ymm12		,%%ymm13			\n\t"\
		"vmovaps	%%ymm0		,    (%%rax)	\n\t		vaddpd	%%ymm8 		,%%ymm10,%%ymm10	\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rax)	\n\t		vsubpd	%%ymm11		,%%ymm8 ,%%ymm8 	\n\t"\
		"vmovaps	%%ymm6		,    (%%rbx)	\n\t		vaddpd	%%ymm9 		,%%ymm12,%%ymm12	\n\t"\
		"vmovaps	%%ymm7		,0x20(%%rbx)	\n\t		vsubpd	%%ymm13		,%%ymm9 ,%%ymm9 	\n\t"\
		"movq		%[add2]	,%%rax				\n\t		vmovaps	%%ymm10		,    (%%r10)	\n\t"\
		"movq		%[add6]	,%%rbx				\n\t		vmovaps	%%ymm12		,0x20(%%r10)	\n\t"\
		"movq		%[c2]	,%%rcx				\n\t		vmovaps	%%ymm8 		,    (%%r11)	\n\t"\
		"movq		%[c6]	,%%rdx				\n\t		vmovaps	%%ymm9 		,0x20(%%r11)	\n\t"\
		"vmovaps	    (%%rax)	,%%ymm0			\n\t		movq		%[add3]	,%%r10				\n\t"\
		"vmovaps	0x20(%%rax)	,%%ymm2			\n\t		movq		%[add7]	,%%r11				\n\t"\
		"vmovaps	    (%%rax)	,%%ymm1			\n\t		movq		%[c3]	,%%r12				\n\t"\
		"vmovaps	0x20(%%rax)	,%%ymm3			\n\t		movq		%[c7]	,%%r13				\n\t"\
		"vmulpd	    (%%rcx)	,%%ymm0,%%ymm0		\n\t		vmovaps	    (%%r10)	,%%ymm8 			\n\t"\
		"vmulpd	0x20(%%rcx)	,%%ymm2,%%ymm2		\n\t		vmovaps	0x20(%%r10)	,%%ymm10			\n\t"\
		"vmulpd	0x20(%%rcx)	,%%ymm1,%%ymm1		\n\t		vmovaps	    (%%r10)	,%%ymm9 			\n\t"\
		"vmulpd	    (%%rcx)	,%%ymm3,%%ymm3		\n\t		vmovaps	0x20(%%r10)	,%%ymm11			\n\t"\
		"vsubpd	%%ymm2		,%%ymm0,%%ymm0		\n\t		vmulpd	    (%%r12)	,%%ymm8 ,%%ymm8 	\n\t"\
		"vaddpd	%%ymm3		,%%ymm1,%%ymm1		\n\t		vmulpd	0x20(%%r12)	,%%ymm10,%%ymm10	\n\t"\
		"vmovaps	    (%%rbx)	,%%ymm2			\n\t		vmulpd	0x20(%%r12)	,%%ymm9 ,%%ymm9 	\n\t"\
		"vmovaps	0x20(%%rbx)	,%%ymm3			\n\t		vmulpd	    (%%r12)	,%%ymm11,%%ymm11	\n\t"\
		"vmovaps	    (%%rbx)	,%%ymm4			\n\t		vsubpd	%%ymm10		,%%ymm8 ,%%ymm8 	\n\t"\
		"vmovaps	0x20(%%rbx)	,%%ymm5			\n\t		vaddpd	%%ymm11		,%%ymm9 ,%%ymm9 	\n\t"\
		"vmulpd	    (%%rdx)	,%%ymm2,%%ymm2		\n\t		vmovaps	    (%%r11)	,%%ymm10			\n\t"\
		"vmulpd	0x20(%%rdx)	,%%ymm3,%%ymm3		\n\t		vmovaps	0x20(%%r11)	,%%ymm11			\n\t"\
		"vmulpd	0x20(%%rdx)	,%%ymm4,%%ymm4		\n\t		vmovaps	    (%%r11)	,%%ymm12			\n\t"\
		"vmulpd	    (%%rdx)	,%%ymm5,%%ymm5		\n\t		vmovaps	0x20(%%r11)	,%%ymm13			\n\t"\
		"vsubpd	%%ymm3		,%%ymm2,%%ymm2		\n\t		vmulpd	    (%%r13)	,%%ymm10,%%ymm10	\n\t"\
		"vaddpd	%%ymm5		,%%ymm4,%%ymm4		\n\t		vmulpd	0x20(%%r13)	,%%ymm11,%%ymm11	\n\t"\
		"vmovaps	%%ymm2		,%%ymm3			\n\t		vmulpd	0x20(%%r13)	,%%ymm12,%%ymm12	\n\t"\
		"vmovaps	%%ymm4		,%%ymm5			\n\t		vmulpd	    (%%r13)	,%%ymm13,%%ymm13	\n\t"\
		"vaddpd	%%ymm0		,%%ymm2,%%ymm2		\n\t		vsubpd	%%ymm11		,%%ymm10,%%ymm10	\n\t"\
		"vsubpd	%%ymm3		,%%ymm0,%%ymm0		\n\t		vaddpd	%%ymm13		,%%ymm12,%%ymm12	\n\t"\
		"vaddpd	%%ymm1		,%%ymm4,%%ymm4		\n\t		vmovaps	%%ymm10		,%%ymm11			\n\t"\
		"vsubpd	%%ymm5		,%%ymm1,%%ymm1		\n\t		vmovaps	%%ymm12		,%%ymm13			\n\t"\
		"vmovaps	%%ymm2		,    (%%rax)	\n\t		vaddpd	%%ymm8 		,%%ymm10,%%ymm10	\n\t"\
		"vmovaps	%%ymm4		,0x20(%%rax)	\n\t		vsubpd	%%ymm11		,%%ymm8 ,%%ymm8 	\n\t"\
		"vmovaps	%%ymm0		,    (%%rbx)	\n\t		vaddpd	%%ymm9 		,%%ymm12,%%ymm12	\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rbx)	\n\t		vsubpd	%%ymm13		,%%ymm9 ,%%ymm9 	\n\t"\
		"													vmovaps	%%ymm10		,    (%%r10)	\n\t"\
		"													vmovaps	%%ymm12		,0x20(%%r10)	\n\t"\
		"													vmovaps	%%ymm8 		,    (%%r11)	\n\t"\
		"													vmovaps	%%ymm9 		,0x20(%%r11)	\n\t"\
/* combine to get 2 length-4 output subtransforms... */\
		"movq		%[add0]	,%%rax				\n\t		movq		%[add4]	,%%r10				\n\t"\
		"movq		%[add2]	,%%rbx				\n\t		movq		%[add6]	,%%r11				\n\t"\
		"vmovaps	    (%%rax)	,%%ymm0			\n\t		vmovaps	    (%%r10)	,%%ymm8 			\n\t"\
		"vmovaps	0x20(%%rax)	,%%ymm1			\n\t		vmovaps	0x20(%%r10)	,%%ymm9 			\n\t"\
		"vmovaps	%%ymm0		,%%ymm4			\n\t		vmovaps	%%ymm8 		,%%ymm12			\n\t"\
		"vmovaps	%%ymm1		,%%ymm5			\n\t		vmovaps	%%ymm9 		,%%ymm13			\n\t"\
		"vaddpd	    (%%rbx)	,%%ymm0,%%ymm0		\n\t		vsubpd	0x20(%%r11)	,%%ymm8 ,%%ymm8 	\n\t"\
		"vsubpd	    (%%rbx)	,%%ymm4,%%ymm4		\n\t		vaddpd	0x20(%%r11)	,%%ymm12,%%ymm12	\n\t"\
		"vaddpd	0x20(%%rbx)	,%%ymm1,%%ymm1		\n\t		vaddpd	    (%%r11)	,%%ymm9 ,%%ymm9 	\n\t"\
		"vsubpd	0x20(%%rbx)	,%%ymm5,%%ymm5		\n\t		vsubpd	    (%%r11)	,%%ymm13,%%ymm13	\n\t"\
		"vmovaps	%%ymm0		,    (%%rax)	\n\t		vmovaps	%%ymm8 		,    (%%r10)	\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rax)	\n\t		vmovaps	%%ymm9 		,0x20(%%r10)	\n\t"\
		"vmovaps	%%ymm4		,    (%%rbx)	\n\t		vmovaps	%%ymm12		,    (%%r11)	\n\t"\
		"vmovaps	%%ymm5		,0x20(%%rbx)	\n\t		vmovaps	%%ymm13		,0x20(%%r11)	\n\t"\
		"movq		%[add1]	,%%rcx				\n\t		movq		%[add5]	,%%r12				\n\t"\
		"movq		%[add3]	,%%rdx				\n\t		movq		%[add7]	,%%r13				\n\t"\
		"vmovaps	    (%%rcx)	,%%ymm2			\n\t		vmovaps	    (%%r12)	,%%ymm10			\n\t"\
		"vmovaps	0x20(%%rcx)	,%%ymm3			\n\t		vmovaps	0x20(%%r12)	,%%ymm11			\n\t"\
		"vmovaps	%%ymm2		,%%ymm6			\n\t		vmovaps	%%ymm10		,%%ymm14			\n\t"\
		"vmovaps	%%ymm3		,%%ymm7			\n\t		vmovaps	%%ymm11		,%%ymm15			\n\t"\
		"vaddpd	    (%%rdx)	,%%ymm2,%%ymm2		\n\t		vsubpd	0x20(%%r13)	,%%ymm10,%%ymm10	\n\t"\
		"vsubpd	    (%%rdx)	,%%ymm6,%%ymm6		\n\t		vaddpd	0x20(%%r13)	,%%ymm14,%%ymm14	\n\t"\
		"vaddpd	0x20(%%rdx)	,%%ymm3,%%ymm3		\n\t		vaddpd	    (%%r13)	,%%ymm11,%%ymm11	\n\t"\
		"vsubpd	0x20(%%rdx)	,%%ymm7,%%ymm7		\n\t		vsubpd	    (%%r13)	,%%ymm15,%%ymm15	\n\t"\
		"vsubpd	%%ymm2		,%%ymm0,%%ymm0		\n\t		vmovaps	%%ymm12		,    (%%r13)	\n\t"\
		"vsubpd	%%ymm3		,%%ymm1,%%ymm1		\n\t		vmovaps	%%ymm13		,0x20(%%r13)	\n\t"\
		"vsubpd	%%ymm7		,%%ymm4,%%ymm4		\n\t		movq		%[isrt2]	,%%r10			\n\t"\
		"vsubpd	%%ymm6		,%%ymm5,%%ymm5		\n\t		vmovaps	%%ymm10		,%%ymm13			\n\t"\
		"vaddpd	    (%%rax)	,%%ymm2,%%ymm2		\n\t		vsubpd	%%ymm11		,%%ymm10,%%ymm10	\n\t"\
		"vaddpd	0x20(%%rax)	,%%ymm3,%%ymm3		\n\t		vaddpd	%%ymm11		,%%ymm13,%%ymm13	\n\t"\
		"vaddpd	    (%%rbx)	,%%ymm7,%%ymm7		\n\t		vmulpd	    (%%r10)	,%%ymm10,%%ymm10	\n\t"\
		"vaddpd	0x20(%%rbx)	,%%ymm6,%%ymm6		\n\t		vmulpd	    (%%r10)	,%%ymm13,%%ymm13	\n\t"\
		"vmovaps	%%ymm2		,    (%%rax)	\n\t		vmovaps	0x20(%%r13)	,%%ymm11			\n\t"\
		"vmovaps	%%ymm3		,0x20(%%rax)	\n\t		vmovaps	%%ymm15		,%%ymm12			\n\t"\
		"vmovaps	%%ymm4		,    (%%rbx)	\n\t		vaddpd	%%ymm14		,%%ymm12,%%ymm12	\n\t"\
		"vmovaps	%%ymm6		,0x20(%%rbx)	\n\t		vsubpd	%%ymm14		,%%ymm15,%%ymm15	\n\t"\
		"vmovaps	%%ymm0		,    (%%rcx)	\n\t		vmulpd	(%%r10)		,%%ymm12,%%ymm12	\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rcx)	\n\t		vmulpd	(%%r10)		,%%ymm15,%%ymm15	\n\t"\
		"vmovaps	%%ymm7		,    (%%rdx)	\n\t		vmovaps		(%%r13)	,%%ymm14			\n\t"\
		"vmovaps	%%ymm5		,0x20(%%rdx)	\n\t		movq		%[add4]	,%%r10				\n\t"\
		"													vsubpd	%%ymm10		,%%ymm8 ,%%ymm8 	\n\t"\
		"													vsubpd	%%ymm13		,%%ymm9 ,%%ymm9 	\n\t"\
		"													vsubpd	%%ymm12		,%%ymm14,%%ymm14	\n\t"\
		"													vsubpd	%%ymm15		,%%ymm11,%%ymm11	\n\t"\
		"													vaddpd	    (%%r10)	,%%ymm10,%%ymm10	\n\t"\
		"													vaddpd	0x20(%%r10)	,%%ymm13,%%ymm13	\n\t"\
		"													vaddpd	    (%%r11)	,%%ymm12,%%ymm12	\n\t"\
		"													vaddpd	0x20(%%r11)	,%%ymm15,%%ymm15	\n\t"\
		"													vmovaps	%%ymm10		,    (%%r10)	\n\t"\
		"													vmovaps	%%ymm13		,0x20(%%r10)	\n\t"\
		"													vmovaps	%%ymm14		,    (%%r11)	\n\t"\
		"													vmovaps	%%ymm11		,0x20(%%r11)	\n\t"\
		"													vmovaps	%%ymm8 		,    (%%r12)	\n\t"\
		"													vmovaps	%%ymm9 		,0x20(%%r12)	\n\t"\
		"													vmovaps	%%ymm12		,    (%%r13)	\n\t"\
		"													vmovaps	%%ymm15		,0x20(%%r13)	\n\t"\
		:					/* outputs: none */\
		: [add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[add1] "m" (Xadd1)\
		 ,[add2] "m" (Xadd2)\
		 ,[add3] "m" (Xadd3)\
		 ,[add4] "m" (Xadd4)\
		 ,[add5] "m" (Xadd5)\
		 ,[add6] "m" (Xadd6)\
		 ,[add7] "m" (Xadd7)\
		 ,[isrt2] "m" (Xisrt2)\
		 ,[c1] "m" (Xc1)\
		 ,[c2] "m" (Xc2)\
		 ,[c3] "m" (Xc3)\
		 ,[c4] "m" (Xc4)\
		 ,[c5] "m" (Xc5)\
		 ,[c6] "m" (Xc6)\
		 ,[c7] "m" (Xc7)\
		: "cc","memory","rax","rbx","rcx","rdx","r10","r11","r12","r13","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
	);\
	}

	#define SSE2_RADIX8_DIT_TWIDDLE(Xadd0,Xadd1,Xadd2,Xadd3,Xadd4,Xadd5,Xadd6,Xadd7,Xisrt2,Xc1,Xc2,Xc3,Xc4,Xc5,Xc6,Xc7)\
	{\
	__asm__ volatile (\
		/*** 2nd of 2 length-4 subtransforms gets done first: ***/\
		/* SSE2_RADIX4_DIT_0TWIDDLE_2NDOFTWO(add4-7):	SSE2_RADIX4_DIT_0TWIDDLE(add0-3): */\
		"movq		%[add4]		,%%rax			\n\t		movq		%[add0]		,%%r10		\n\t"\
		"movq		%[add5]		,%%rbx			\n\t		movq		%[add1]		,%%r11		\n\t"\
		"vmovaps		(%%rax)	,%%ymm0			\n\t		vmovaps			(%%r10)	,%%ymm8 	\n\t"\
		"vmovaps	0x20(%%rax)	,%%ymm1			\n\t		vmovaps		0x20(%%r10)	,%%ymm9 	\n\t"\
		"vmovaps	%%ymm0		,%%ymm2			\n\t		vmovaps		%%ymm8 		,%%ymm10	\n\t"\
		"vmovaps	%%ymm1		,%%ymm3			\n\t		vmovaps		%%ymm9 		,%%ymm11	\n\t"\
		"vaddpd		(%%rbx)	,%%ymm2,%%ymm2		\n\t		vaddpd		(%%r11)	,%%ymm10,%%ymm10\n\t"\
		"vaddpd	0x20(%%rbx)	,%%ymm3,%%ymm3		\n\t		vaddpd	0x20(%%r11)	,%%ymm11,%%ymm11\n\t"\
		"vsubpd		(%%rbx)	,%%ymm0,%%ymm0		\n\t		vsubpd		(%%r11)	,%%ymm8 ,%%ymm8 \n\t"\
		"vsubpd	0x20(%%rbx)	,%%ymm1,%%ymm1		\n\t		vsubpd	0x20(%%r11)	,%%ymm9 ,%%ymm9 \n\t"\
		"movq		%[add6]		,%%rcx			\n\t		movq		%[add2]		,%%r12		\n\t"\
		"movq		%[add7]		,%%rdx			\n\t		movq		%[add3]		,%%r13		\n\t"\
		"vmovaps			(%%rcx)	,%%ymm4		\n\t		vmovaps			(%%r12)	,%%ymm12	\n\t"\
		"vmovaps		0x20(%%rcx)	,%%ymm5		\n\t		vmovaps		0x20(%%r12)	,%%ymm13	\n\t"/* t9  in ymm6, needed below in RHS! */\
		"vmovaps		%%ymm4		,%%ymm6		\n\t		vmovaps		%%ymm12		,%%ymm14	\n\t"/* t10 in ymm7, needed below in RHS! */\
		"vmovaps		%%ymm5		,%%ymm7		\n\t		vmovaps		%%ymm13		,%%ymm15	\n\t"\
		"vaddpd		(%%rdx)	,%%ymm6,%%ymm6		\n\t		vaddpd		(%%r13)	,%%ymm14,%%ymm14\n\t"\
		"vaddpd	0x20(%%rdx)	,%%ymm7,%%ymm7		\n\t		vaddpd	0x20(%%r13)	,%%ymm15,%%ymm15\n\t"\
		"vsubpd		(%%rdx)	,%%ymm4,%%ymm4		\n\t		vsubpd		(%%r13)	,%%ymm12,%%ymm12\n\t"\
		"vsubpd	0x20(%%rdx)	,%%ymm5,%%ymm5		\n\t		vsubpd	0x20(%%r13)	,%%ymm13,%%ymm13\n\t"\
		"vmovaps		%%ymm6	,    (%%rcx)	\n\t		vmovaps		%%ymm14		,    (%%r12)\n\t"\
		"vmovaps		%%ymm7	,0x20(%%rcx)	\n\t		vmovaps		%%ymm15		,0x20(%%r12)\n\t"\
		"vmovaps		%%ymm4	,    (%%rdx)	\n\t		vaddpd	%%ymm10		,%%ymm14,%%ymm14\n\t"/* ymm14 <- ~t1 */\
		"vmovaps		%%ymm5	,0x20(%%rdx)	\n\t		vaddpd	%%ymm11		,%%ymm15,%%ymm15\n\t"/* ymm15 <- ~t2 */\
		"vaddpd	%%ymm2		,%%ymm6,%%ymm6		\n\t		vsubpd		(%%r12)	,%%ymm10,%%ymm10\n\t"/* ymm10 <- ~t5 */\
		"vaddpd	%%ymm3		,%%ymm7,%%ymm7		\n\t		vsubpd	0x20(%%r12)	,%%ymm11,%%ymm11\n\t"/* ymm11 <- ~t6 */\
		"vsubpd		(%%rcx)	,%%ymm2,%%ymm2		\n\t		vaddpd	%%ymm6,%%ymm14,%%ymm14		\n\t"/* t1+t9 */\
		"vsubpd	0x20(%%rcx)	,%%ymm3,%%ymm3		\n\t		vaddpd	%%ymm7,%%ymm15,%%ymm15		\n\t"/* t2+t10*/\
		"vmovaps		%%ymm2	,    (%%rcx)	\n\t		vaddpd	%%ymm6,%%ymm6 ,%%ymm6		\n\t"/* 2*t9  */\
		"vmovaps		%%ymm3	,0x20(%%rcx)	\n\t		vaddpd	%%ymm7,%%ymm7 ,%%ymm7		\n\t"/* 2*t10 */\
		"vmovaps		%%ymm4	,%%ymm2			\n\t		vmovaps		%%ymm14		,    (%%r10)\n\t"/* a[j1   ], DONE. */\
		"vmovaps		%%ymm5	,%%ymm3			\n\t		vmovaps		%%ymm15		,0x20(%%r10)\n\t"/* a[j2   ], DONE. */\
		"vaddpd	%%ymm0		,%%ymm5,%%ymm5		\n\t		vsubpd	%%ymm6,%%ymm14,%%ymm14		\n\t"/* t1-t9  = [t1+t9 ] - 2*t9  */\
		"vsubpd	%%ymm3		,%%ymm0,%%ymm0		\n\t		vsubpd	%%ymm7,%%ymm15,%%ymm15		\n\t"/* t2-t10 = [t2+t10] - 2*t10 */\
		"vaddpd	%%ymm1		,%%ymm4,%%ymm4		\n\t		vmovaps		%%ymm12		,%%ymm6		\n\t"/* ymm6<- copy of t7 */\
		"vsubpd	%%ymm2		,%%ymm1,%%ymm1		\n\t		vmovaps		%%ymm13		,%%ymm7		\n\t"/* ymm7<- copy of t8 */\
		"vmovaps		%%ymm5	,%%ymm2			\n\t		vaddpd	%%ymm8 		,%%ymm13,%%ymm13\n\t"/* ymm13<- ~t3 */\
		"vmovaps		%%ymm1	,%%ymm3			\n\t		vsubpd	%%ymm7		,%%ymm8 ,%%ymm8	\n\t"/* ymm8 <- ~t7 */\
		"vaddpd	%%ymm1		,%%ymm5,%%ymm5		\n\t		vaddpd	%%ymm9 		,%%ymm12,%%ymm12\n\t"/* ymm12<- ~t8 */\
		"movq		%[isrt2]	,%%rsi			\n\t		vsubpd	%%ymm6		,%%ymm9 ,%%ymm9 \n\t"/* ymm9 <- ~t4 */\
		"vmovaps			(%%rsi)	,%%ymm1		\n\t"\
		"vsubpd	%%ymm3		,%%ymm2,%%ymm2		\n\t"		/* Combine Outputs 0,4 of two half-transforms, as these are ready: */\
		"vmulpd	%%ymm1		,%%ymm5,%%ymm5		\n\t		movq		%[c4]		,%%rsi		\n\t"\
		"vmulpd	%%ymm1		,%%ymm2,%%ymm2		\n\t		vmovaps		%%ymm14,%%ymm6			\n\t"\
		"vmovaps		%%ymm0	,%%ymm3			\n\t		vmovaps		%%ymm15,%%ymm7			\n\t"\
		"vmovaps		%%ymm5	,    (%%rbx)	\n\t		vmulpd	0x20(%%rsi)	,%%ymm6	,%%ymm6 \n\t"\
		"vmovaps		%%ymm2	,0x20(%%rbx)	\n\t		vmulpd	0x20(%%rsi)	,%%ymm7	,%%ymm7 \n\t"\
		"vmovaps		%%ymm4	,%%ymm5			\n\t		vmulpd		(%%rsi)	,%%ymm14,%%ymm14\n\t"\
		"vaddpd	%%ymm4		,%%ymm0,%%ymm0		\n\t		vmulpd		(%%rsi)	,%%ymm15,%%ymm15\n\t"\
		"vsubpd	%%ymm5		,%%ymm3,%%ymm3		\n\t		vaddpd	%%ymm7 ,%%ymm14,%%ymm14		\n\t"\
		"vmulpd	%%ymm1		,%%ymm0,%%ymm0		\n\t		vsubpd	%%ymm6 ,%%ymm15,%%ymm15		\n\t"\
		"vmulpd	%%ymm1		,%%ymm3,%%ymm3		\n\t		vmovaps		%%ymm14,    (%%rax)		\n\t"\
		"vmovaps		%%ymm0	,    (%%rdx)	\n\t		vmovaps		%%ymm15,0x20(%%rax)		\n\t"\
		"vmovaps		%%ymm3	,0x20(%%rdx)	\n\t"\
	/* Now combine the two half-transforms & store outputs back into original array slots. */\
	/* add0-7 in r10,11,12,13,ax,bx,cx,dx; SIMD Registers 0-7,14-15 FREE: */\
		/* Outputs 1,5: Use ymm 4,5,9,13,14,15			Outputs 2,6: Use ymm 0,1,2,3,6,7,10,11 : */\
		"movq		%[c2],%%rsi		\n\t"/* c6 = c2+2, c1 = c2+4, c5 = c2+6 */\
		"vmovaps		    (%%rbx)	,%%ymm4 	\n\t		vmovaps		    (%%rcx)	,%%ymm2		\n\t"\
		"vmovaps		0x20(%%rbx)	,%%ymm5 	\n\t		vmovaps		0x20(%%rcx)	,%%ymm3		\n\t"\
		"vmovaps		%%ymm13	,%%ymm14		\n\t		vmovaps		%%ymm10		,%%ymm6		\n\t"\
		"vmovaps		%%ymm9 	,%%ymm15		\n\t		vmovaps		%%ymm11		,%%ymm7		\n\t"\
		"vaddpd	%%ymm4 	,%%ymm13,%%ymm13		\n\t		vaddpd	%%ymm3	,%%ymm10,%%ymm10	\n\t"\
		"vsubpd	%%ymm5 	,%%ymm9 ,%%ymm9 		\n\t		vsubpd	%%ymm2	,%%ymm11,%%ymm11	\n\t"\
		"vsubpd	%%ymm4 	,%%ymm14,%%ymm14		\n\t"	/*	vmovaps		    (%%rsi)	,%%ymm0	Need ymm0 to replace ymm9  below */\
		"vaddpd	%%ymm5 	,%%ymm15,%%ymm15		\n\t"	/*	vmovaps		0x20(%%rsi)	,%%ymm1	Need ymm1 to replace ymm13 below */\
		"vmovaps		%%ymm14	,%%ymm4 		\n\t		vsubpd	%%ymm3	,%%ymm6	,%%ymm6		\n\t"\
		"vmovaps		%%ymm15	,%%ymm5 		\n\t		vaddpd	%%ymm2	,%%ymm7	,%%ymm7		\n\t"/* ymm2,3 free */\
		"vmovaps		%%ymm13	,%%ymm14		\n\t		vmovaps		%%ymm10		,%%ymm0 	\n\t"\
		"vmovaps		%%ymm9 	,%%ymm15		\n\t		vmovaps		%%ymm11		,%%ymm1		\n\t"\
		"vmulpd	0x80(%%rsi)	,%%ymm13,%%ymm13	\n\t		vmulpd		(%%rsi)	,%%ymm10,%%ymm10\n\t"\
		"vmulpd	0x80(%%rsi)	,%%ymm9 ,%%ymm9 	\n\t		vmulpd		(%%rsi)	,%%ymm11,%%ymm11\n\t"\
		"vmulpd	0xa0(%%rsi)	,%%ymm14,%%ymm14	\n\t		vmulpd	0x20(%%rsi)	,%%ymm0 ,%%ymm0 \n\t"\
		"vmulpd	0xa0(%%rsi)	,%%ymm15,%%ymm15	\n\t		vmulpd	0x20(%%rsi)	,%%ymm1	,%%ymm1	\n\t"\
		"vsubpd	%%ymm14		,%%ymm9 ,%%ymm9 	\n\t		vsubpd	%%ymm0 		,%%ymm11,%%ymm11\n\t"\
		"vaddpd	%%ymm15		,%%ymm13,%%ymm13	\n\t		vaddpd	%%ymm1		,%%ymm10,%%ymm10\n\t"\
		"vmovaps		%%ymm9 	,0x20(%%r11)	\n\t		vmovaps		0x40(%%rsi)	,%%ymm2		\n\t"\
		"vmovaps		%%ymm13	,    (%%r11)	\n\t		vmovaps		0x60(%%rsi)	,%%ymm3		\n\t"\
		"vmovaps		%%ymm4 	,%%ymm13		\n\t		vmovaps		%%ymm11		,0x20(%%r12)\n\t"\
		"vmovaps		%%ymm5 	,%%ymm9 		\n\t		vmovaps		%%ymm10		,    (%%r12)\n\t"\
		"vmovaps		%%ymm13	,%%ymm14		\n\t		vmovaps		%%ymm6		,%%ymm0 	\n\t"\
		"vmovaps		%%ymm9 	,%%ymm15		\n\t		vmovaps		%%ymm7		,%%ymm1		\n\t"\
		"vmulpd	0xc0(%%rsi)	,%%ymm13,%%ymm13	\n\t		vmulpd	%%ymm2	,%%ymm6	,%%ymm6 	\n\t"\
		"vmulpd	0xc0(%%rsi)	,%%ymm9 ,%%ymm9 	\n\t		vmulpd	%%ymm2	,%%ymm7	,%%ymm7		\n\t"\
		"vmulpd	0xe0(%%rsi)	,%%ymm14,%%ymm14	\n\t		vmulpd	%%ymm3	,%%ymm0 ,%%ymm0		\n\t"\
		"vmulpd	0xe0(%%rsi)	,%%ymm15,%%ymm15	\n\t		vmulpd	%%ymm3	,%%ymm1	,%%ymm1		\n\t"\
		"vsubpd	%%ymm14		,%%ymm9 ,%%ymm9 	\n\t		vsubpd	%%ymm0 	,%%ymm7	,%%ymm7		\n\t"\
		"vaddpd	%%ymm15		,%%ymm13,%%ymm13	\n\t		vaddpd	%%ymm1	,%%ymm6	,%%ymm6		\n\t"\
		"vmovaps		%%ymm9 	,0x20(%%rbx)	\n\t		vmovaps		%%ymm7		,0x20(%%rcx)\n\t"\
		"vmovaps		%%ymm13	,    (%%rbx)	\n\t		vmovaps		%%ymm6		,    (%%rcx)\n\t"\
		/* Outputs 3,7: All SIMD regs except for ymm8,12 (containing one of the 2 complex data being butterflied) free: */\
		"movq		%[c3]		,%%rsi			\n\t"/* c7 = c3+2 */\
		"vmovaps		    (%%rdx)	,%%ymm4		\n\t"\
		"vmovaps		0x20(%%rdx)	,%%ymm5		\n\t"\
		"vmovaps		%%ymm8 	,%%ymm14		\n\t"\
		"vmovaps		%%ymm12	,%%ymm15		\n\t"\
		"vsubpd	%%ymm5	,%%ymm8 ,%%ymm8			\n\t		vsubpd	%%ymm4	,%%ymm12,%%ymm12	\n\t"\
		"vmovaps		    (%%rsi)	,%%ymm0		\n\t"\
		"vmovaps		0x20(%%rsi)	,%%ymm1		\n\t"\
		"vaddpd	%%ymm5	,%%ymm14,%%ymm14		\n\t		vmovaps		0x40(%%rsi)	,%%ymm2		\n\t"\
		"vaddpd	%%ymm4	,%%ymm15,%%ymm15		\n\t		vmovaps		0x60(%%rsi)	,%%ymm3		\n\t"\
		"vmovaps		%%ymm8 	,%%ymm9 		\n\t		vmovaps		%%ymm14		,%%ymm6 	\n\t"\
		"vmovaps		%%ymm12	,%%ymm13		\n\t		vmovaps		%%ymm15		,%%ymm7		\n\t"\
		"vmulpd	%%ymm0	,%%ymm8 ,%%ymm8 		\n\t		vmulpd	%%ymm2	,%%ymm14,%%ymm14	\n\t"\
		"vmulpd	%%ymm0	,%%ymm12,%%ymm12		\n\t		vmulpd	%%ymm2	,%%ymm15,%%ymm15	\n\t"\
		"vmulpd	%%ymm1	,%%ymm9 ,%%ymm9 		\n\t		vmulpd	%%ymm3	,%%ymm6 ,%%ymm6 	\n\t"\
		"vmulpd	%%ymm1	,%%ymm13,%%ymm13		\n\t		vmulpd	%%ymm3	,%%ymm7	,%%ymm7		\n\t"\
		"vsubpd	%%ymm9 		,%%ymm12,%%ymm12	\n\t		vsubpd	%%ymm6 	,%%ymm15,%%ymm15	\n\t"\
		"vaddpd	%%ymm13		,%%ymm8 ,%%ymm8 	\n\t		vaddpd	%%ymm7	,%%ymm14,%%ymm14	\n\t"\
		"vmovaps		%%ymm12	,0x20(%%r13)	\n\t		vmovaps		%%ymm15		,0x20(%%rdx)\n\t"\
		"vmovaps		%%ymm8 	,    (%%r13)	\n\t		vmovaps		%%ymm14		,    (%%rdx)\n\t"\
		:					/* outputs: none */\
		: [add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[add1] "m" (Xadd1)\
		 ,[add2] "m" (Xadd2)\
		 ,[add3] "m" (Xadd3)\
		 ,[add4] "m" (Xadd4)\
		 ,[add5] "m" (Xadd5)\
		 ,[add6] "m" (Xadd6)\
		 ,[add7] "m" (Xadd7)\
		 ,[isrt2] "m" (Xisrt2)\
		 ,[c1] "m" (Xc1)\
		 ,[c2] "m" (Xc2)\
		 ,[c3] "m" (Xc3)\
		 ,[c4] "m" (Xc4)\
		 ,[c5] "m" (Xc5)\
		 ,[c6] "m" (Xc6)\
		 ,[c7] "m" (Xc7)\
		: "cc","memory","rax","rbx","rcx","rdx","rsi","r10","r11","r12","r13","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
	);\
	}

#elif defined(USE_SSE2) && (OS_BITS == 64)

	// To-do: propagate load-saving tweaks from 8-reg version to this one:
	#define SSE2_RADIX8_DIF_TWIDDLE(Xadd0,Xadd1,Xadd2,Xadd3,Xadd4,Xadd5,Xadd6,Xadd7,Xisrt2,Xc1,Xc2,Xc3,Xc4,Xc5,Xc6,Xc7)\
	{\
	__asm__ volatile (\
		"												movq		%[add1]	,%%r10			\n\t"\
		"												movq		%[add5]	,%%r11			\n\t"\
		"												movq		%[c1]	,%%r12			\n\t"\
		"												movq		%[c5]	,%%r13			\n\t"\
		"												movaps		    (%%r10)	,%%xmm8 	\n\t"\
		"movq		%[add0]	,%%rax			\n\t		movaps		0x10(%%r10)	,%%xmm10	\n\t"\
		"movq		%[add4]	,%%rbx			\n\t		movaps		    (%%r10)	,%%xmm9 	\n\t"\
		"movq		%[c4]	,%%rcx			\n\t		movaps		0x10(%%r10)	,%%xmm11	\n\t"\
		"movaps		    (%%rax)	,%%xmm0		\n\t		mulpd			(%%r12)	,%%xmm8 	\n\t"\
		"movaps		0x10(%%rax)	,%%xmm1		\n\t		mulpd		0x10(%%r12)	,%%xmm10	\n\t"\
		"movaps		    (%%rax)	,%%xmm6		\n\t		mulpd		0x10(%%r12)	,%%xmm9 	\n\t"\
		"movaps		0x10(%%rax)	,%%xmm7		\n\t		mulpd			(%%r12)	,%%xmm11	\n\t"\
		"movaps		    (%%rbx)	,%%xmm2		\n\t		subpd		%%xmm10		,%%xmm8 	\n\t"\
		"movaps		0x10(%%rbx)	,%%xmm3		\n\t		addpd		%%xmm11		,%%xmm9 	\n\t"\
		"movaps		    (%%rbx)	,%%xmm4		\n\t		movaps		    (%%r11)	,%%xmm10	\n\t"\
		"movaps		0x10(%%rbx)	,%%xmm5		\n\t		movaps		0x10(%%r11)	,%%xmm11	\n\t"\
		"mulpd			(%%rcx)	,%%xmm2		\n\t		movaps		    (%%r11)	,%%xmm12	\n\t"\
		"mulpd			(%%rcx)	,%%xmm3		\n\t		movaps		0x10(%%r11)	,%%xmm13	\n\t"\
		"mulpd		0x10(%%rcx)	,%%xmm4		\n\t		mulpd			(%%r13)	,%%xmm10	\n\t"\
		"mulpd		0x10(%%rcx)	,%%xmm5		\n\t		mulpd		0x10(%%r13)	,%%xmm11	\n\t"\
		"subpd		%%xmm5		,%%xmm2		\n\t		mulpd		0x10(%%r13)	,%%xmm12	\n\t"\
		"addpd		%%xmm4		,%%xmm3		\n\t		mulpd			(%%r13)	,%%xmm13	\n\t"\
		"addpd		%%xmm2		,%%xmm0		\n\t		subpd		%%xmm11		,%%xmm10	\n\t"\
		"addpd		%%xmm3		,%%xmm1		\n\t		addpd		%%xmm13		,%%xmm12	\n\t"\
		"subpd		%%xmm2		,%%xmm6		\n\t		movaps		%%xmm10		,%%xmm11	\n\t"\
		"subpd		%%xmm3		,%%xmm7		\n\t		movaps		%%xmm12		,%%xmm13	\n\t"\
		"movaps		%%xmm0	,    (%%rax)	\n\t		addpd		%%xmm8 		,%%xmm10	\n\t"\
		"movaps		%%xmm1	,0x10(%%rax)	\n\t		subpd		%%xmm11		,%%xmm8 	\n\t"\
		"movaps		%%xmm6	,    (%%rbx)	\n\t		addpd		%%xmm9 		,%%xmm12	\n\t"\
		"movaps		%%xmm7	,0x10(%%rbx)	\n\t		subpd		%%xmm13		,%%xmm9 	\n\t"\
		"movq		%[add2]	,%%rax			\n\t		movaps		%%xmm10	,    (%%r10)	\n\t"\
		"movq		%[add6]	,%%rbx			\n\t		movaps		%%xmm12	,0x10(%%r10)	\n\t"\
		"movq		%[c2]	,%%rcx			\n\t		movaps		%%xmm8 	,    (%%r11)	\n\t"\
		"movq		%[c6]	,%%rdx			\n\t		movaps		%%xmm9 	,0x10(%%r11)	\n\t"\
		"movaps		    (%%rax)	,%%xmm0		\n\t		movq		%[add3]	,%%r10			\n\t"\
		"movaps		0x10(%%rax)	,%%xmm2		\n\t		movq		%[add7]	,%%r11			\n\t"\
		"movaps		    (%%rax)	,%%xmm1		\n\t		movq		%[c3]	,%%r12			\n\t"\
		"movaps		0x10(%%rax)	,%%xmm3		\n\t		movq		%[c7]	,%%r13			\n\t"\
		"mulpd			(%%rcx)	,%%xmm0		\n\t		movaps		    (%%r10)	,%%xmm8 	\n\t"\
		"mulpd		0x10(%%rcx)	,%%xmm2		\n\t		movaps		0x10(%%r10)	,%%xmm10	\n\t"\
		"mulpd		0x10(%%rcx)	,%%xmm1		\n\t		movaps		    (%%r10)	,%%xmm9 	\n\t"\
		"mulpd			(%%rcx)	,%%xmm3		\n\t		movaps		0x10(%%r10)	,%%xmm11	\n\t"\
		"subpd		%%xmm2		,%%xmm0		\n\t		mulpd			(%%r12)	,%%xmm8 	\n\t"\
		"addpd		%%xmm3		,%%xmm1		\n\t		mulpd		0x10(%%r12)	,%%xmm10	\n\t"\
		"movaps		    (%%rbx)	,%%xmm2		\n\t		mulpd		0x10(%%r12)	,%%xmm9 	\n\t"\
		"movaps		0x10(%%rbx)	,%%xmm3		\n\t		mulpd			(%%r12)	,%%xmm11	\n\t"\
		"movaps		    (%%rbx)	,%%xmm4		\n\t		subpd		%%xmm10		,%%xmm8 	\n\t"\
		"movaps		0x10(%%rbx)	,%%xmm5		\n\t		addpd		%%xmm11		,%%xmm9 	\n\t"\
		"mulpd			(%%rdx)	,%%xmm2		\n\t		movaps		    (%%r11)	,%%xmm10	\n\t"\
		"mulpd		0x10(%%rdx)	,%%xmm3		\n\t		movaps		0x10(%%r11)	,%%xmm11	\n\t"\
		"mulpd		0x10(%%rdx)	,%%xmm4		\n\t		movaps		    (%%r11)	,%%xmm12	\n\t"\
		"mulpd			(%%rdx)	,%%xmm5		\n\t		movaps		0x10(%%r11)	,%%xmm13	\n\t"\
		"subpd		%%xmm3		,%%xmm2		\n\t		mulpd			(%%r13)	,%%xmm10	\n\t"\
		"addpd		%%xmm5		,%%xmm4		\n\t		mulpd		0x10(%%r13)	,%%xmm11	\n\t"\
		"movaps		%%xmm2		,%%xmm3		\n\t		mulpd		0x10(%%r13)	,%%xmm12	\n\t"\
		"movaps		%%xmm4		,%%xmm5		\n\t		mulpd			(%%r13)	,%%xmm13	\n\t"\
		"addpd		%%xmm0		,%%xmm2		\n\t		subpd		%%xmm11		,%%xmm10	\n\t"\
		"subpd		%%xmm3		,%%xmm0		\n\t		addpd		%%xmm13		,%%xmm12	\n\t"\
		"addpd		%%xmm1		,%%xmm4		\n\t		movaps		%%xmm10		,%%xmm11	\n\t"\
		"subpd		%%xmm5		,%%xmm1		\n\t		movaps		%%xmm12		,%%xmm13	\n\t"\
		"movaps		%%xmm2	,    (%%rax)	\n\t		addpd		%%xmm8 		,%%xmm10	\n\t"\
		"movaps		%%xmm4	,0x10(%%rax)	\n\t		subpd		%%xmm11		,%%xmm8 	\n\t"\
		"movaps		%%xmm0	,    (%%rbx)	\n\t		addpd		%%xmm9 		,%%xmm12	\n\t"\
		"movaps		%%xmm1	,0x10(%%rbx)	\n\t		subpd		%%xmm13		,%%xmm9 	\n\t"\
		"												movaps		%%xmm10	,    (%%r10)	\n\t"\
		"												movaps		%%xmm12	,0x10(%%r10)	\n\t"\
		"												movaps		%%xmm8 	,    (%%r11)	\n\t"\
		"												movaps		%%xmm9 	,0x10(%%r11)	\n\t"\
/* combine to get 2 length-4 output subtransforms... */\
		"movq		%[add0]	,%%rax			\n\t		movq		%[add4]	,%%r10			\n\t"\
		"movq		%[add2]	,%%rbx			\n\t		movq		%[add6]	,%%r11			\n\t"\
		"movaps		    (%%rax)	,%%xmm0		\n\t		movaps		    (%%r10)	,%%xmm8 	\n\t"\
		"movaps		0x10(%%rax)	,%%xmm1		\n\t		movaps		0x10(%%r10)	,%%xmm9 	\n\t"\
		"movaps		%%xmm0		,%%xmm4		\n\t		movaps		%%xmm8 		,%%xmm12	\n\t"\
		"movaps		%%xmm1		,%%xmm5		\n\t		movaps		%%xmm9 		,%%xmm13	\n\t"\
		"addpd			(%%rbx)	,%%xmm0		\n\t		subpd		0x10(%%r11)	,%%xmm8 	\n\t"\
		"subpd			(%%rbx)	,%%xmm4		\n\t		addpd		0x10(%%r11)	,%%xmm12	\n\t"\
		"addpd		0x10(%%rbx)	,%%xmm1		\n\t		addpd			(%%r11)	,%%xmm9 	\n\t"\
		"subpd		0x10(%%rbx)	,%%xmm5		\n\t		subpd			(%%r11)	,%%xmm13	\n\t"\
		"movaps		%%xmm0	,    (%%rax)	\n\t		movaps		%%xmm8 	,    (%%r10)	\n\t"\
		"movaps		%%xmm1	,0x10(%%rax)	\n\t		movaps		%%xmm9 	,0x10(%%r10)	\n\t"\
		"movaps		%%xmm4	,    (%%rbx)	\n\t		movaps		%%xmm12	,    (%%r11)	\n\t"\
		"movaps		%%xmm5	,0x10(%%rbx)	\n\t		movaps		%%xmm13	,0x10(%%r11)	\n\t"\
		"movq		%[add1]	,%%rcx			\n\t		movq		%[add5]	,%%r12			\n\t"\
		"movq		%[add3]	,%%rdx			\n\t		movq		%[add7]	,%%r13			\n\t"\
		"movaps		    (%%rcx)	,%%xmm2		\n\t		movaps		    (%%r12)	,%%xmm10	\n\t"\
		"movaps		0x10(%%rcx)	,%%xmm3		\n\t		movaps		0x10(%%r12)	,%%xmm11	\n\t"\
		"movaps		%%xmm2		,%%xmm6		\n\t		movaps		%%xmm10		,%%xmm14	\n\t"\
		"movaps		%%xmm3		,%%xmm7		\n\t		movaps		%%xmm11		,%%xmm15	\n\t"\
		"addpd			(%%rdx)	,%%xmm2		\n\t		subpd		0x10(%%r13)	,%%xmm10	\n\t"\
		"subpd			(%%rdx)	,%%xmm6		\n\t		addpd		0x10(%%r13)	,%%xmm14	\n\t"\
		"addpd		0x10(%%rdx)	,%%xmm3		\n\t		addpd			(%%r13)	,%%xmm11	\n\t"\
		"subpd		0x10(%%rdx)	,%%xmm7		\n\t		subpd			(%%r13)	,%%xmm15	\n\t"\
		"subpd		%%xmm2		,%%xmm0		\n\t		movaps		%%xmm12	,    (%%r13)	\n\t"\
		"subpd		%%xmm3		,%%xmm1		\n\t		movaps		%%xmm13	,0x10(%%r13)	\n\t"\
		"subpd		%%xmm7		,%%xmm4		\n\t		movq		%[isrt2]	,%%r10		\n\t"\
		"subpd		%%xmm6		,%%xmm5		\n\t		movaps		%%xmm10		,%%xmm13	\n\t"\
		"addpd			(%%rax)	,%%xmm2		\n\t		subpd		%%xmm11		,%%xmm10	\n\t"\
		"addpd		0x10(%%rax)	,%%xmm3		\n\t		addpd		%%xmm11		,%%xmm13	\n\t"\
		"addpd			(%%rbx)	,%%xmm7		\n\t		mulpd			(%%r10)	,%%xmm10	\n\t"\
		"addpd		0x10(%%rbx)	,%%xmm6		\n\t		mulpd			(%%r10)	,%%xmm13	\n\t"\
		"movaps		%%xmm2	,    (%%rax)	\n\t		movaps		0x10(%%r13)	,%%xmm11	\n\t"\
		"movaps		%%xmm3	,0x10(%%rax)	\n\t		movaps		%%xmm15		,%%xmm12	\n\t"\
		"movaps		%%xmm4	,    (%%rbx)	\n\t		addpd		%%xmm14		,%%xmm12	\n\t"\
		"movaps		%%xmm6	,0x10(%%rbx)	\n\t		subpd		%%xmm14		,%%xmm15	\n\t"\
		"movaps		%%xmm0	,    (%%rcx)	\n\t		mulpd			(%%r10)	,%%xmm12	\n\t"\
		"movaps		%%xmm1	,0x10(%%rcx)	\n\t		mulpd			(%%r10)	,%%xmm15	\n\t"\
		"movaps		%%xmm7	,    (%%rdx)	\n\t		movaps			(%%r13)	,%%xmm14	\n\t"\
		"movaps		%%xmm5	,0x10(%%rdx)	\n\t		movq		%[add4]	,%%r10			\n\t"\
		"												subpd		%%xmm10		,%%xmm8 	\n\t"\
		"												subpd		%%xmm13		,%%xmm9 	\n\t"\
		"												subpd		%%xmm12		,%%xmm14	\n\t"\
		"												subpd		%%xmm15		,%%xmm11	\n\t"\
		"												addpd			(%%r10)	,%%xmm10	\n\t"\
		"												addpd		0x10(%%r10)	,%%xmm13	\n\t"\
		"												addpd			(%%r11)	,%%xmm12	\n\t"\
		"												addpd		0x10(%%r11)	,%%xmm15	\n\t"\
		"												movaps		%%xmm10	,    (%%r10)	\n\t"\
		"												movaps		%%xmm13	,0x10(%%r10)	\n\t"\
		"												movaps		%%xmm14	,    (%%r11)	\n\t"\
		"												movaps		%%xmm11	,0x10(%%r11)	\n\t"\
		"												movaps		%%xmm8 	,    (%%r12)	\n\t"\
		"												movaps		%%xmm9 	,0x10(%%r12)	\n\t"\
		"												movaps		%%xmm12	,    (%%r13)	\n\t"\
		"												movaps		%%xmm15	,0x10(%%r13)	\n\t"\
		:					/* outputs: none */\
		: [add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[add1] "m" (Xadd1)\
		 ,[add2] "m" (Xadd2)\
		 ,[add3] "m" (Xadd3)\
		 ,[add4] "m" (Xadd4)\
		 ,[add5] "m" (Xadd5)\
		 ,[add6] "m" (Xadd6)\
		 ,[add7] "m" (Xadd7)\
		 ,[isrt2] "m" (Xisrt2)\
		 ,[c1] "m" (Xc1)\
		 ,[c2] "m" (Xc2)\
		 ,[c3] "m" (Xc3)\
		 ,[c4] "m" (Xc4)\
		 ,[c5] "m" (Xc5)\
		 ,[c6] "m" (Xc6)\
		 ,[c7] "m" (Xc7)\
		: "cc","memory","rax","rbx","rcx","rdx","r10","r11","r12","r13","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
	);\
	}

	// Note: The SIMD radix-8 DIF is a very basic one, no opts re. memory movement, etc. Thus DIF has 65% more load/stores than DIT
	// DIF SIMD opcount: 140 load/store [56 implicit], 66 add/sub, 32 mul
	// DIT SIMD opcount:  85 load/store [36 implicit], 68 add/sub, 32 mul	<*** To-Do: propagate same optimizations used in DIT to DIF!

	#define SSE2_RADIX8_DIT_TWIDDLE(Xadd0,Xadd1,Xadd2,Xadd3,Xadd4,Xadd5,Xadd6,Xadd7,Xisrt2,Xc1,Xc2,Xc3,Xc4,Xc5,Xc6,Xc7)\
	{\
	__asm__ volatile (\
		/*** 2nd of 2 length-4 subtransforms gets done first: ***/\
		/* SSE2_RADIX4_DIT_0TWIDDLE_2NDOFTWO(add4-7):	SSE2_RADIX4_DIT_0TWIDDLE(add0-3): */\
		"movq		%[add4]		,%%rax		\n\t		movq		%[add0]		,%%r10		\n\t"\
		"movq		%[add5]		,%%rbx		\n\t		movq		%[add1]		,%%r11		\n\t"\
		"movaps		(%%rax)	,%%xmm0			\n\t		movaps			(%%r10)	,%%xmm8 		\n\t"\
		"movaps	0x10(%%rax)	,%%xmm1			\n\t		movaps		0x10(%%r10)	,%%xmm9 		\n\t"\
		"movaps	%%xmm0		,%%xmm2			\n\t		movaps		%%xmm8 		,%%xmm10		\n\t"\
		"movaps	%%xmm1		,%%xmm3			\n\t		movaps		%%xmm9 		,%%xmm11		\n\t"\
		"addpd			(%%rbx)	,%%xmm2		\n\t		addpd			(%%r11)	,%%xmm10		\n\t"\
		"addpd		0x10(%%rbx)	,%%xmm3		\n\t		addpd		0x10(%%r11)	,%%xmm11		\n\t"\
		"subpd			(%%rbx)	,%%xmm0		\n\t		subpd			(%%r11)	,%%xmm8 		\n\t"\
		"subpd		0x10(%%rbx)	,%%xmm1		\n\t		subpd		0x10(%%r11)	,%%xmm9 		\n\t"\
		"movq		%[add6]		,%%rcx		\n\t		movq		%[add2]		,%%r12			\n\t"\
		"movq		%[add7]		,%%rdx		\n\t		movq		%[add3]		,%%r13			\n\t"\
		"movaps			(%%rcx)	,%%xmm4		\n\t		movaps			(%%r12)	,%%xmm12		\n\t"\
		"movaps		0x10(%%rcx)	,%%xmm5		\n\t		movaps		0x10(%%r12)	,%%xmm13		\n\t"/* t9  in xmm6, needed below in RHS! */\
		"movaps		%%xmm4		,%%xmm6		\n\t		movaps		%%xmm12		,%%xmm14		\n\t"/* t10 in xmm7, needed below in RHS! */\
		"movaps		%%xmm5		,%%xmm7		\n\t		movaps		%%xmm13		,%%xmm15		\n\t"\
		"addpd			(%%rdx)	,%%xmm6		\n\t		addpd			(%%r13)	,%%xmm14		\n\t"\
		"addpd		0x10(%%rdx)	,%%xmm7		\n\t		addpd		0x10(%%r13)	,%%xmm15		\n\t"\
		"subpd			(%%rdx)	,%%xmm4		\n\t		subpd			(%%r13)	,%%xmm12		\n\t"\
		"subpd		0x10(%%rdx)	,%%xmm5		\n\t		subpd		0x10(%%r13)	,%%xmm13		\n\t"\
		"movaps		%%xmm6		,    (%%rcx)\n\t		movaps		%%xmm14		,    (%%r12)	\n\t"\
		"movaps		%%xmm7		,0x10(%%rcx)\n\t		movaps		%%xmm15		,0x10(%%r12)	\n\t"\
		"movaps		%%xmm4		,    (%%rdx)\n\t		addpd		%%xmm10		,%%xmm14		\n\t"/* xmm14 <- ~t1 */\
		"movaps		%%xmm5		,0x10(%%rdx)\n\t		addpd		%%xmm11		,%%xmm15		\n\t"/* xmm15 <- ~t2 */\
		"addpd		%%xmm2		,%%xmm6		\n\t		subpd			(%%r12)	,%%xmm10		\n\t"/* xmm10 <- ~t5 */\
		"addpd		%%xmm3		,%%xmm7		\n\t		subpd		0x10(%%r12)	,%%xmm11		\n\t"/* xmm11 <- ~t6 */\
		"subpd			(%%rcx)	,%%xmm2		\n\t		addpd		%%xmm6,%%xmm14		\n\t"/* t1+t9 */\
		"subpd		0x10(%%rcx)	,%%xmm3		\n\t		addpd		%%xmm7,%%xmm15		\n\t"/* t2+t10*/\
		"movaps		%%xmm2		,    (%%rcx)\n\t		addpd		%%xmm6,%%xmm6		\n\t"/* 2*t9  */\
		"movaps		%%xmm3		,0x10(%%rcx)\n\t		addpd		%%xmm7,%%xmm7		\n\t"/* 2*t10 */\
		"movaps		%%xmm4		,%%xmm2		\n\t		movaps		%%xmm14		,    (%%r10)	\n\t"/* a[j1   ], DONE. */\
		"movaps		%%xmm5		,%%xmm3		\n\t		movaps		%%xmm15		,0x10(%%r10)	\n\t"/* a[j2   ], DONE. */\
		"addpd		%%xmm0		,%%xmm5		\n\t		subpd		%%xmm6,%%xmm14		\n\t"/* t1-t9  = [t1+t9 ] - 2*t9  */\
		"subpd		%%xmm3		,%%xmm0		\n\t		subpd		%%xmm7,%%xmm15		\n\t"/* t2-t10 = [t2+t10] - 2*t10 */\
		"addpd		%%xmm1		,%%xmm4		\n\t		movaps		%%xmm12		,%%xmm6		\n\t"/* xmm6<- copy of t7 */\
		"subpd		%%xmm2		,%%xmm1		\n\t		movaps		%%xmm13		,%%xmm7		\n\t"/* xmm7<- copy of t8 */\
		"movaps		%%xmm5		,%%xmm2		\n\t		addpd		%%xmm8 		,%%xmm13		\n\t"/* xmm13<- ~t3 */\
		"movaps		%%xmm1		,%%xmm3		\n\t		subpd		%%xmm7		,%%xmm8 		\n\t"/* xmm8 <- ~t7 */\
		"addpd		%%xmm1		,%%xmm5		\n\t		addpd		%%xmm9 		,%%xmm12		\n\t"/* xmm12<- ~t8 */\
		"movq		%[isrt2]	,%%rsi		\n\t		subpd		%%xmm6		,%%xmm9 		\n\t"/* xmm9 <- ~t4 */\
		"movaps			(%%rsi)	,%%xmm1		\n\t"\
		"subpd		%%xmm3		,%%xmm2		\n\t"		/* Combine Outputs 0,4 of two half-transforms, as these are ready: */\
		"mulpd		%%xmm1		,%%xmm5		\n\t		movq		%[c4]		,%%rsi		\n\t"\
		"mulpd		%%xmm1		,%%xmm2		\n\t		movaps		%%xmm14,%%xmm6			\n\t"\
		"movaps		%%xmm0		,%%xmm3		\n\t		movaps		%%xmm15,%%xmm7			\n\t"\
		"movaps		%%xmm5		,    (%%rbx)\n\t		mulpd		0x10(%%rsi)	,%%xmm6		\n\t"\
		"movaps		%%xmm2		,0x10(%%rbx)\n\t		mulpd		0x10(%%rsi)	,%%xmm7		\n\t"\
		"movaps		%%xmm4		,%%xmm5		\n\t		mulpd			(%%rsi)	,%%xmm14	\n\t"\
		"addpd		%%xmm4		,%%xmm0		\n\t		mulpd			(%%rsi)	,%%xmm15	\n\t"\
		"subpd		%%xmm5		,%%xmm3		\n\t		addpd		%%xmm7 ,%%xmm14			\n\t"\
		"mulpd		%%xmm1		,%%xmm0		\n\t		subpd		%%xmm6 ,%%xmm15			\n\t"\
		"mulpd		%%xmm1		,%%xmm3		\n\t		movaps		%%xmm14,    (%%rax)		\n\t"\
		"movaps		%%xmm0		,    (%%rdx)\n\t		movaps		%%xmm15,0x10(%%rax)		\n\t"\
		"movaps		%%xmm3		,0x10(%%rdx)\n\t"\
	/* Now combine the two half-transforms & store outputs back into original array slots. */\
	/* add0-7 in r10,11,12,13,ax,bx,cx,dx; SIMD Registers 0-7,14-15 FREE: */\
		/* Outputs 1,5: Use xmm 4,5,9,13,14,15			Outputs 2,6: Use xmm 0,1,2,3,6,7,10,11 : */\
		"movq		%[c2],%%rsi		\n\t"/* c6 = c2+2, c1 = c2+4, c5 = c2+6 */\
		"movaps		    (%%rbx)	,%%xmm4 	\n\t		movaps		    (%%rcx)	,%%xmm2		\n\t"\
		"movaps		0x10(%%rbx)	,%%xmm5 	\n\t		movaps		0x10(%%rcx)	,%%xmm3		\n\t"\
		"movaps		%%xmm13		,%%xmm14	\n\t		movaps		%%xmm10		,%%xmm6		\n\t"\
		"movaps		%%xmm9 		,%%xmm15	\n\t		movaps		%%xmm11		,%%xmm7		\n\t"\
		"addpd		%%xmm4 	,%%xmm13		\n\t		addpd		%%xmm3	,%%xmm10		\n\t"\
		"subpd		%%xmm5 	,%%xmm9 		\n\t		subpd		%%xmm2	,%%xmm11		\n\t"\
		"subpd		%%xmm4 	,%%xmm14		\n\t"	/*	movaps		    (%%rsi)	,%%xmm0	Need xmm0 to replace xmm9  below */\
		"addpd		%%xmm5 	,%%xmm15		\n\t"	/*	movaps		0x10(%%rsi)	,%%xmm1	Need xmm1 to replace xmm13 below */\
		"movaps		%%xmm14		,%%xmm4 	\n\t		subpd		%%xmm3	,%%xmm6		\n\t"\
		"movaps		%%xmm15		,%%xmm5 	\n\t		addpd		%%xmm2	,%%xmm7		\n\t"/* xmm2,3 free */\
		"movaps		%%xmm13		,%%xmm14	\n\t		movaps		%%xmm10		,%%xmm0 		\n\t"\
		"movaps		%%xmm9 		,%%xmm15	\n\t		movaps		%%xmm11		,%%xmm1		\n\t"\
		"mulpd		0x40(%%rsi)	,%%xmm13	\n\t		mulpd			(%%rsi)	,%%xmm10		\n\t"\
		"mulpd		0x40(%%rsi)	,%%xmm9 	\n\t		mulpd			(%%rsi)	,%%xmm11		\n\t"\
		"mulpd		0x50(%%rsi)	,%%xmm14	\n\t		mulpd		0x10(%%rsi)	,%%xmm0 		\n\t"\
		"mulpd		0x50(%%rsi)	,%%xmm15	\n\t		mulpd		0x10(%%rsi)	,%%xmm1		\n\t"\
		"subpd		%%xmm14		,%%xmm9 	\n\t		subpd		%%xmm0 		,%%xmm11		\n\t"\
		"addpd		%%xmm15		,%%xmm13	\n\t		addpd		%%xmm1		,%%xmm10		\n\t"\
		"movaps		%%xmm9 		,0x10(%%r11)\n\t		movaps		0x20(%%rsi)	,%%xmm2		\n\t"\
		"movaps		%%xmm13		,    (%%r11)\n\t		movaps		0x30(%%rsi)	,%%xmm3		\n\t"\
		"movaps		%%xmm4 	,%%xmm13		\n\t		movaps		%%xmm11		,0x10(%%r12)	\n\t"\
		"movaps		%%xmm5 	,%%xmm9 		\n\t		movaps		%%xmm10		,    (%%r12)	\n\t"\
		"movaps		%%xmm13		,%%xmm14	\n\t		movaps		%%xmm6		,%%xmm0 		\n\t"\
		"movaps		%%xmm9 		,%%xmm15	\n\t		movaps		%%xmm7		,%%xmm1		\n\t"\
		"mulpd		0x60(%%rsi)	,%%xmm13	\n\t		mulpd		%%xmm2	,%%xmm6		\n\t"\
		"mulpd		0x60(%%rsi)	,%%xmm9 	\n\t		mulpd		%%xmm2	,%%xmm7		\n\t"\
		"mulpd		0x70(%%rsi)	,%%xmm14	\n\t		mulpd		%%xmm3	,%%xmm0 		\n\t"\
		"mulpd		0x70(%%rsi)	,%%xmm15	\n\t		mulpd		%%xmm3	,%%xmm1		\n\t"\
		"subpd		%%xmm14		,%%xmm9 	\n\t		subpd		%%xmm0 		,%%xmm7		\n\t"\
		"addpd		%%xmm15		,%%xmm13	\n\t		addpd		%%xmm1		,%%xmm6		\n\t"\
		"movaps		%%xmm9 		,0x10(%%rbx)\n\t		movaps		%%xmm7		,0x10(%%rcx)	\n\t"\
		"movaps		%%xmm13		,    (%%rbx)\n\t		movaps		%%xmm6		,    (%%rcx)	\n\t"\
		/* Outputs 3,7: All SIMD regs except for xmm8,12 (containing one of the 2 complex data being butterflied) free: */\
		"movq		%[c3]		,%%rsi		\n\t"/* c7 = c3+2 */\
		"movaps		    (%%rdx)	,%%xmm4		\n\t"\
		"movaps		0x10(%%rdx)	,%%xmm5		\n\t"\
		"movaps		%%xmm8 		,%%xmm14	\n\t"\
		"movaps		%%xmm12		,%%xmm15	\n\t"\
		"subpd		%%xmm5	,%%xmm8 		\n\t		subpd		%%xmm4	,%%xmm12		\n\t"\
		"movaps		    (%%rsi)	,%%xmm0		\n\t"\
		"movaps		0x10(%%rsi)	,%%xmm1		\n\t"\
		"addpd		%%xmm5	,%%xmm14		\n\t		movaps		0x20(%%rsi)	,%%xmm2		\n\t"\
		"addpd		%%xmm4	,%%xmm15		\n\t		movaps		0x30(%%rsi)	,%%xmm3		\n\t"\
		"movaps		%%xmm8 		,%%xmm9 	\n\t		movaps		%%xmm14		,%%xmm6 	\n\t"\
		"movaps		%%xmm12		,%%xmm13	\n\t		movaps		%%xmm15		,%%xmm7		\n\t"\
		"mulpd		%%xmm0	,%%xmm8 		\n\t		mulpd		%%xmm2	,%%xmm14		\n\t"\
		"mulpd		%%xmm0	,%%xmm12		\n\t		mulpd		%%xmm2	,%%xmm15		\n\t"\
		"mulpd		%%xmm1	,%%xmm9 		\n\t		mulpd		%%xmm3	,%%xmm6 		\n\t"\
		"mulpd		%%xmm1	,%%xmm13		\n\t		mulpd		%%xmm3	,%%xmm7			\n\t"\
		"subpd		%%xmm9 		,%%xmm12	\n\t		subpd		%%xmm6 		,%%xmm15	\n\t"\
		"addpd		%%xmm13		,%%xmm8 	\n\t		addpd		%%xmm7		,%%xmm14	\n\t"\
		"movaps		%%xmm12		,0x10(%%r13)\n\t		movaps		%%xmm15		,0x10(%%rdx)\n\t"\
		"movaps		%%xmm8 		,    (%%r13)\n\t		movaps		%%xmm14		,    (%%rdx)\n\t"\
		:					/* outputs: none */\
		: [add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[add1] "m" (Xadd1)\
		 ,[add2] "m" (Xadd2)\
		 ,[add3] "m" (Xadd3)\
		 ,[add4] "m" (Xadd4)\
		 ,[add5] "m" (Xadd5)\
		 ,[add6] "m" (Xadd6)\
		 ,[add7] "m" (Xadd7)\
		 ,[isrt2] "m" (Xisrt2)\
		 ,[c1] "m" (Xc1)\
		 ,[c2] "m" (Xc2)\
		 ,[c3] "m" (Xc3)\
		 ,[c4] "m" (Xc4)\
		 ,[c5] "m" (Xc5)\
		 ,[c6] "m" (Xc6)\
		 ,[c7] "m" (Xc7)\
		: "cc","memory","rax","rbx","rcx","rdx","rsi","r10","r11","r12","r13","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
	);\
	}

#elif defined(USE_SSE2) && (OS_BITS == 32)

	#error 32-bit OSes no longer supported for SIMD builds!

#else

	#error Unhandled combination of preprocessr flags!

#endif	// x86 simd version ?

#endif	/* radix8_dif_dit_pass_asm_h_included */

