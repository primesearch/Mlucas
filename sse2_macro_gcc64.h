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
#ifndef sse2_macro_gcc_h_included
#define sse2_macro_gcc_h_included

/*
SSE2-ified version of PAIR_SQUARE_4. Data enter in [__tAr, ~tAr], [__tAi, ~tAi], pairs, where the imaginary part
of each input pair is assumed offset +0x10 in memory from the real part, i.e. needs no explicit pointer reference.

For the sincos twiddles: using the notation of the scalar PAIR_SQUARE_4() macro,"__c" means [c0,s1], "__s" means [s0,c1].
For these, due to the buterfly indexing pattern, we cannot assume that __s = __c + 0x10, so feed both pointers explicitly.

We use shufpd xmm, xmm, 1 to swap lo and hi doubles of an xmm register for the various operations with one swapped input.
*/
	/* Complex multiply of 2 roots of unity - use e.g. for "multiply up" of sincos twiddles. */
	#define SSE2_CMUL_EXPO(XcA,XcB,XcAmB,XcApB)\
	{\
	__asm__ volatile (\
		"movq	%[__cA]		,%%rax\n\t"\
		"movq	%[__cB]		,%%rbx\n\t"\
		"movq	%[__cAmB]	,%%rcx\n\t"\
		"movq	%[__cApB]	,%%rdx\n\t"\
		"\n\t"\
		"movaps	    (%%rax),%%xmm0\n\t"\
		"movaps	0x10(%%rax),%%xmm2\n\t"\
		"movaps	    (%%rbx),%%xmm4\n\t"\
		"movaps	0x10(%%rbx),%%xmm5\n\t"\
		"movaps	%%xmm0,%%xmm1\n\t"\
		"movaps	%%xmm2,%%xmm3\n\t"\
		"\n\t"\
		"mulpd	%%xmm4,%%xmm0\n\t"\
		"mulpd	%%xmm5,%%xmm1\n\t"\
		"mulpd	%%xmm4,%%xmm2\n\t"\
		"mulpd	%%xmm5,%%xmm3\n\t"\
		"movaps	%%xmm0,%%xmm4\n\t"\
		"movaps	%%xmm1,%%xmm5\n\t"\
		"addpd	%%xmm3,%%xmm0\n\t"\
		"subpd	%%xmm2,%%xmm1\n\t"\
		"subpd	%%xmm3,%%xmm4\n\t"\
		"addpd	%%xmm2,%%xmm5\n\t"\
		"movaps	%%xmm0,    (%%rcx)\n\t"\
		"movaps	%%xmm1,0x10(%%rcx)\n\t"\
		"movaps	%%xmm4,    (%%rdx)\n\t"\
		"movaps	%%xmm5,0x10(%%rdx)\n\t"\
		:					/* outputs: none */\
		: [__cA]  "m" (XcA)	/* All inputs from memory addresses here */\
		 ,[__cB]  "m" (XcB)\
		 ,[__cAmB] "m" (XcAmB)\
		 ,[__cApB] "m" (XcApB)\
		: "rax","rbx","rcx","rdx","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5"		/* Clobbered registers */\
	);\
	}

	#define PAIR_SQUARE_4_SSE2(XtAr, XtBr, XtCr, XtDr, Xc, Xs, Xforth)\
	{\
	__asm__ volatile (\
		"/*   calculate cross-product terms...\n\t"\
		"	__rt=__tAr* ~tDr+__tAi* ~tDi; __rt=__rt+__rt;\n\t"\
		"	__it=__tAi* ~tDr-__tAr* ~tDi; __it=__it+__it;\n\t"\
		"*/\n\t"\
		"movq	%[__tDr]	,%%rdx\n\t"\
		"movq	%[__tAr]	,%%rax\n\t"\
		"\n\t"\
		"movaps	    (%%rdx)	,%%xmm6		/* tDr */\n\t"\
		"movaps	0x10(%%rdx)	,%%xmm7		/* tDi */\n\t"\
		"movaps	    (%%rax)	,%%xmm0		/* tAr */\n\t"\
		"movaps	0x10(%%rax)	,%%xmm3		/* tAi */\n\t"\
		"shufpd	$1	,%%xmm6	,%%xmm6	/*~tDr */\n\t"\
		"shufpd	$1	,%%xmm7	,%%xmm7	/*~tDi */\n\t"\
		"movaps	    (%%rax)	,%%xmm2		/* cpy tAr */\n\t"\
		"movaps	0x10(%%rax)	,%%xmm1		/* cpy tAi */\n\t"\
		"\n\t"\
		"mulpd	%%xmm6		,%%xmm0	/* tAr*~tDr */\n\t"\
		"mulpd	%%xmm7		,%%xmm3	/* tAi*~tDi */\n\t"\
		"mulpd	%%xmm6		,%%xmm1	/* tAi*~tDr */\n\t"\
		"mulpd	%%xmm7		,%%xmm2	/* tAr*~tDi */\n\t"\
		"addpd	%%xmm3		,%%xmm0	/* rt */\n\t"\
		"subpd	%%xmm2		,%%xmm1	/* it */\n\t"\
		"addpd	%%xmm0		,%%xmm0	/* rt=rt+rt */\n\t"\
		"addpd	%%xmm1		,%%xmm1	/* it=it+it; xmm2-7 free */\n\t"\
		"/*\n\t"\
		"	__st=__tBr* ~tCr+__tBi* ~tCi; __st=__st+__st;\n\t"\
		"	__jt=__tBi* ~tCr-__tBr* ~tCi; __jt=__jt+__jt;\n\t"\
		"*/\n\t"\
		"movq	%[__tCr]	,%%rcx\n\t"\
		"movq	%[__tBr]	,%%rbx\n\t"\
		"\n\t"\
		"movaps	    (%%rcx)	,%%xmm6		/* tCr */\n\t"\
		"movaps	0x10(%%rcx)	,%%xmm7		/* tCi */\n\t"\
		"movaps	    (%%rbx)	,%%xmm2		/* tBr */\n\t"\
		"movaps	0x10(%%rbx)	,%%xmm5		/* tBi */\n\t"\
		"shufpd	$1	,%%xmm6	,%%xmm6	/*~tCr */\n\t"\
		"shufpd	$1	,%%xmm7	,%%xmm7	/*~tCi */\n\t"\
		"movaps	    (%%rbx)	,%%xmm4		/* cpy tBr */\n\t"\
		"movaps	0x10(%%rbx)	,%%xmm3		/* cpy tBi */\n\t"\
		"\n\t"\
		"mulpd	%%xmm6		,%%xmm2	/* tBr*~tCr */\n\t"\
		"mulpd	%%xmm7		,%%xmm5	/* tBi*~tCi */\n\t"\
		"mulpd	%%xmm6		,%%xmm3	/* tBi*~tCr */\n\t"\
		"mulpd	%%xmm7		,%%xmm4	/* tBr*~tCi */\n\t"\
		"addpd	%%xmm5		,%%xmm2	/* st */\n\t"\
		"subpd	%%xmm4		,%%xmm3	/* jt */\n\t"\
		"addpd	%%xmm2		,%%xmm2	/* st=st+st */\n\t"\
		"addpd	%%xmm3		,%%xmm3	/* jt=jt+jt; xmm4-7 free */\n\t"\
		"\n\t"\
		"/*   now calculate square terms and __store back in the same temporaries:	*/\n\t"\
		"/*	__tmp=(__tAr+__tAi)*(__tAr-__tAi); __tAi=__tAr*__tAi; __tAi=__tAi+__tAi; __tAr=__tmp;	*/\n\t"\
		"\n\t"\
		"movaps	    (%%rax)	,%%xmm4		/* __tAr */\n\t"\
		"movaps	0x10(%%rax)	,%%xmm5		/* __tAi */\n\t"\
		"subpd	%%xmm5		,%%xmm4		/* (__tAr-__tAi) */\n\t"\
		"addpd	%%xmm5		,%%xmm5		/*      2*__tAi  */\n\t"\
		"addpd	%%xmm4		,%%xmm5		/* (__tAr+__tAi) */\n\t"\
		"mulpd	%%xmm5		,%%xmm4		/*>__tAr */\n\t"\
		"\n\t"\
		"movaps	    (%%rax)	,%%xmm5		/* __tAr */\n\t"\
		"mulpd	0x10(%%rax)	,%%xmm5		/* __tAr*__tAi */\n\t"\
		"addpd	%%xmm5		,%%xmm5		/*>__tAi */\n\t"\
		"movaps	%%xmm4	,    (%%rax)	/* tmp store >__tAr */\n\t"\
		"movaps	%%xmm5	,0x10(%%rax)	/* tmp store >__tAi */\n\t"\
		"\n\t"\
		"subpd	%%xmm4		,%%xmm0	/* rt-__tAr */\n\t"\
		"subpd	%%xmm5		,%%xmm1	/* it-__tAi; xmm4-7 free */\n\t"\
		"\n\t"\
		"/*	__tmp=(__tBr+__tBi)*(__tBr-__tBi); __tBi=__tBr*__tBi; __tBi=__tBi+__tBi; __tBr=__tmp;	*/\n\t"\
		"/*** [Can be done in parallel with above segment] ***/\n\t"\
		"\n\t"\
		"movaps	    (%%rbx)	,%%xmm6		/* __tBr */\n\t"\
		"movaps	0x10(%%rbx)	,%%xmm7		/* __tBi */\n\t"\
		"subpd	%%xmm7		,%%xmm6		/* (__tBr-__tBi) */\n\t"\
		"addpd	%%xmm7		,%%xmm7		/*      2*__tBi  */\n\t"\
		"addpd	%%xmm6		,%%xmm7		/* (__tBr+__tBi) */\n\t"\
		"mulpd	%%xmm7		,%%xmm6		/*>__tBr */\n\t"\
		"\n\t"\
		"movaps	    (%%rbx)	,%%xmm7		/* __tBr */\n\t"\
		"mulpd	0x10(%%rbx)	,%%xmm7		/* __tBr*__tBi */\n\t"\
		"addpd	%%xmm7		,%%xmm7		/*>__tBi */\n\t"\
		"movaps	%%xmm6	,    (%%rbx)	/* tmp store >__tBr */\n\t"\
		"movaps	%%xmm7	,0x10(%%rbx)	/* tmp store >__tBi */\n\t"\
		"\n\t"\
		"subpd	%%xmm6		,%%xmm2	/* st-__tBr */\n\t"\
		"subpd	%%xmm7		,%%xmm3	/* jt-__tBi; xmm4-7 free */\n\t"\
		"\n\t"\
		"/*	__tmp=(__tDr+__tDi)*(__tDr-__tDi); __tDi=__tDr*__tDi; __tDi=__tDi+__tDi; __tDr=__tmp;	*/\n\t"\
		"\n\t"\
		"movaps	    (%%rdx)	,%%xmm4		/* __tDr */\n\t"\
		"movaps	0x10(%%rdx)	,%%xmm5		/* __tDi */\n\t"\
		"subpd	%%xmm5		,%%xmm4		/* (__tDr-__tDi) */\n\t"\
		"addpd	%%xmm5		,%%xmm5		/*      2*__tDi  */\n\t"\
		"addpd	%%xmm4		,%%xmm5		/* (__tDr+__tDi) */\n\t"\
		"mulpd	%%xmm5		,%%xmm4		/*>__tDr */\n\t"\
		"\n\t"\
		"movaps	    (%%rdx)	,%%xmm5		/* __tDr */\n\t"\
		"mulpd	0x10(%%rdx)	,%%xmm5		/* __tDr*__tDi */\n\t"\
		"addpd	%%xmm5		,%%xmm5		/*>__tDi */\n\t"\
		"movaps	%%xmm4	,    (%%rdx)	/* tmp store ~tDr */\n\t"\
		"movaps	%%xmm5	,0x10(%%rdx)	/* tmp store ~tDi */\n\t"\
		"shufpd	$1	,%%xmm4	,%%xmm4	/*~tDr */\n\t"\
		"shufpd	$1	,%%xmm5	,%%xmm5	/*~tDi */\n\t"\
		"\n\t"\
		"subpd	%%xmm4		,%%xmm0	/* rt-__tAr- ~tDr */\n\t"\
		"addpd	%%xmm5		,%%xmm1	/* it-__tAi+ ~tDi; xmm4-7 free */\n\t"\
		"\n\t"\
		"/*	__tmp=(__tCr+__tCi)*(__tCr-__tCi); __tCi=__tCr*__tCi; __tCi=__tCi+__tCi; __tCr=__tmp;	*/\n\t"\
		"/*** [Can be done in parallel with above segment] ***/\n\t"\
		"\n\t"\
		"movaps	    (%%rcx)	,%%xmm6		/* __tCr */\n\t"\
		"movaps	0x10(%%rcx)	,%%xmm7		/* __tCi */\n\t"\
		"subpd	%%xmm7		,%%xmm6		/* (__tCr-__tCi) */\n\t"\
		"addpd	%%xmm7		,%%xmm7		/*      2*__tCi  */\n\t"\
		"addpd	%%xmm6		,%%xmm7		/* (__tCr+__tCi) */\n\t"\
		"mulpd	%%xmm7		,%%xmm6		/*>__tCr */\n\t"\
		"\n\t"\
		"movaps	    (%%rcx)	,%%xmm7		/* __tCr */\n\t"\
		"mulpd	0x10(%%rcx)	,%%xmm7		/* __tCr*__tCi */\n\t"\
		"addpd	%%xmm7		,%%xmm7		/*>__tCi */\n\t"\
		"movaps	%%xmm6	,    (%%rcx)	/* tmp store ~tCr */\n\t"\
		"movaps	%%xmm7	,0x10(%%rcx)	/* tmp store ~tCi */\n\t"\
		"shufpd	$1	,%%xmm6	,%%xmm6	/*~tCr */\n\t"\
		"shufpd	$1	,%%xmm7	,%%xmm7	/*~tCi */\n\t"\
		"\n\t"\
		"subpd	%%xmm6		,%%xmm2	/* st-__tBr- ~tCr */\n\t"\
		"addpd	%%xmm7		,%%xmm3	/* jt-__tBi+ ~tCi; xmm4-7 free */\n\t"\
		"/*\n\t"\
		"	__tmp=((1.0+__c)*__rt-__s*__it)*0.25;\n\t"\
		"	__it =((1.0+__c)*__it+__s*__rt)*0.25;	__rt=__tmp;\n\t"\
		"*/\n\t"\
		"/*** [Can be done in parallel with above segment] ***/\n\t"\
		"movq	%[__c]		,%%rax\n\t"\
		"movq	%[__s]		,%%rbx\n\t"\
		"movq	%[__forth]	,%%rdx\n\t"\
		"movaps	%%xmm0		,%%xmm4		/* cpy rt */\n\t"\
		"movaps	%%xmm1		,%%xmm5		/* cpy it */\n\t"\
		"mulpd	(%%rax)		,%%xmm0		/* c*rt */\n\t"\
		"mulpd	(%%rax)		,%%xmm1		/* c*it */\n\t"\
		"addpd	%%xmm4		,%%xmm0		/* (c+1.0)*rt */\n\t"\
		"addpd	%%xmm5		,%%xmm1		/* (c+1.0)*it */\n\t"\
		"mulpd	(%%rbx)		,%%xmm4		/* s*rt */\n\t"\
		"mulpd	(%%rbx)		,%%xmm5		/* s*it */\n\t"\
		"subpd	%%xmm5		,%%xmm0		/* (c+1.0)*rt-s*it */\n\t"\
		"addpd	%%xmm4		,%%xmm1		/* (c+1.0)*it+s*rt; xmm4,5 free */\n\t"\
		"mulpd	(%%rdx)		,%%xmm0	/* -rt Both of these inherit the sign flip [w.r.to the non-SSE2 PAIR_SQUARE_4 macro] */\n\t"\
		"mulpd	(%%rdx)		,%%xmm1	/* -it that resulted from the in-place-friendlier (rt-__tAr- ~tDr) reordering above. */\n\t"\
		"/*\n\t"\
		"	__tmp=((1.0-__s)*__st-__c*__jt)*0.25;\n\t"\
		"	__jt =((1.0-__s)*__jt+__c*__st)*0.25	__st=__tmp;\n\t"\
		"*/\n\t"\
		"/*** [Can be done in parallel wjth above segment] ***/\n\t"\
		"movaps	%%xmm2		,%%xmm6		/* cpy st */\n\t"\
		"movaps	%%xmm3		,%%xmm7		/* cpy jt */\n\t"\
		"mulpd	(%%rbx)		,%%xmm2		/* s*st */\n\t"\
		"mulpd	(%%rbx)		,%%xmm3		/* s*jt */\n\t"\
		"subpd	%%xmm6		,%%xmm2		/* (s-1.0)*st, note sign flip! */\n\t"\
		"subpd	%%xmm7		,%%xmm3		/* (s-1.0)*jt, note sign flip! */\n\t"\
		"mulpd	(%%rax)		,%%xmm6		/* c*st */\n\t"\
		"mulpd	(%%rax)		,%%xmm7		/* c*jt */\n\t"\
		"addpd	%%xmm7		,%%xmm2		/* -[(1.0-s)*st-c*jt] */\n\t"\
		"subpd	%%xmm6		,%%xmm3		/* -[(1.0-s)*jt+c*st]; xmm6,7 free */\n\t"\
		"mulpd	(%%rdx)		,%%xmm2	/* +st Sign flip due to (s-1.0) reordering here */\n\t"\
		"mulpd	(%%rdx)		,%%xmm3	/* +jt cancels earlier one due to in-place-friendlier (st-__tBr- ~tCr) reordering above. */\n\t"\
		"/*...and now complete and store the results. We flip the signs on st and jt here to undo the above -st,-jt negations. */\n\t"\
		"/*	__tAr = (__tAr+__rt);\n\t"\
		"	__tAi = (__tAi+__it);\n\t"\
		"	__tBr = (__tBr-__st);\n\t"\
		"	__tBi = (__tBi-__jt);\n\t"\
		"*/\n\t"\
		"movq	%[__tAr]	,%%rax\n\t"\
		"movq	%[__tBr]	,%%rbx\n\t"\
		"\n\t"\
		"movaps	    (%%rax)	,%%xmm4		/* __tAr */\n\t"\
		"movaps	0x10(%%rax)	,%%xmm5		/* __tAi */\n\t"\
		"movaps	    (%%rbx)	,%%xmm6		/* __tBr */\n\t"\
		"movaps	0x10(%%rbx)	,%%xmm7		/* __tBi */\n\t"\
		"addpd	%%xmm0		,%%xmm4		/* (__tAr+__rt) */\n\t"\
		"addpd	%%xmm1		,%%xmm5		/* (__tAi+__it) */\n\t"\
		"subpd	%%xmm2		,%%xmm6		/* (__tBr-__st) */\n\t"\
		"subpd	%%xmm3		,%%xmm7		/* (__tBi-__jt) */\n\t"\
		"movaps	%%xmm4	,    (%%rax)	/* store >__tAr */\n\t"\
		"movaps	%%xmm5	,0x10(%%rax)	/* store >__tAi */\n\t"\
		"movaps	%%xmm6	,    (%%rbx)	/* store >__tBr */\n\t"\
		"movaps	%%xmm7	,0x10(%%rbx)	/* store >__tBi */\n\t"\
		"/*...N-j terms are as above, but with the replacements: __tAr<--> ~tDr, __tAi<--> ~tDi, __it|-->-__it. */\n\t"\
		"/*	__tDr = (__tDr+ ~rt);\n\t"\
		"	__tDi = (__tDi- ~it);\n\t"\
		"	__tCr = (__tCr- ~st);\n\t"\
		"	__tCi = (__tCi+ ~jt);\n\t"\
		"*/\n\t"\
		"movq	%[__tCr]	,%%rcx\n\t"\
		"movq	%[__tDr]	,%%rdx\n\t"\
		"\n\t"\
		"shufpd	$1	,%%xmm0	,%%xmm0		/* ~rt */\n\t"\
		"shufpd	$1	,%%xmm1	,%%xmm1		/* ~it */\n\t"\
		"shufpd	$1	,%%xmm2	,%%xmm2		/* ~st */\n\t"\
		"shufpd	$1	,%%xmm3	,%%xmm3		/* ~jt */\n\t"\
		"\n\t"\
		"movaps	    (%%rdx)	,%%xmm4		/* __tDr */\n\t"\
		"movaps	0x10(%%rdx)	,%%xmm5		/* __tDi */\n\t"\
		"movaps	    (%%rcx)	,%%xmm6		/* __tCr */\n\t"\
		"movaps	0x10(%%rcx)	,%%xmm7		/* __tCi */\n\t"\
		"addpd	%%xmm0		,%%xmm4		/* (__tDr+ ~rt) */\n\t"\
		"subpd	%%xmm1		,%%xmm5		/* (__tDi- ~it) */\n\t"\
		"subpd	%%xmm2		,%%xmm6		/* (__tCr- ~st) */\n\t"\
		"addpd	%%xmm3		,%%xmm7		/* (__tCi+ ~jt) */\n\t"\
		"movaps	%%xmm4	,    (%%rdx)	/* store >__tDr */\n\t"\
		"movaps	%%xmm5	,0x10(%%rdx)	/* store >__tDi */\n\t"\
		"movaps	%%xmm6	,    (%%rcx)	/* store >__tCr */\n\t"\
		"movaps	%%xmm7	,0x10(%%rcx)	/* store >__tCi */\n\t"\
		:					/* outputs: none */\
		: [__tAr] "m" (XtAr)	/* All inputs from memory addresses here */\
		 ,[__tBr] "m" (XtBr)\
		 ,[__tCr] "m" (XtCr)\
		 ,[__tDr] "m" (XtDr)\
		 ,[__c] "m" (Xc)\
		 ,[__s] "m" (Xs)\
		 ,[__forth] "m" (Xforth)\
		: "rax","rbx","rcx","rdx","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
	);\
	}

	#define	SSE2_RADIX8_DIT_0TWIDDLE(Xin0,Xin1,Xin2,Xin3,Xin4,Xin5,Xin6,Xin7, Xout, Xisrt2)\
	{\
	__asm__ volatile (\
		"movq	%[__in4],%%rax	\n\t"\
		"movq	%[__in5],%%rbx	\n\t"\
		"movq	%[__in6],%%rcx	\n\t"\
		"movq	%[__in7],%%rdx	\n\t"\
		"movq	%[__out],%%rsi	\n\t"\
		"movaps	    (%%rax),%%xmm0			\n\t"\
		"movaps	0x10(%%rax),%%xmm1			\n\t"\
		"movaps	%%xmm0,%%xmm2				\n\t"\
		"movaps	%%xmm1,%%xmm3				\n\t"\
		"addpd	    (%%rbx),%%xmm2			\n\t"\
		"addpd	0x10(%%rbx),%%xmm3			\n\t"\
		"subpd	    (%%rbx),%%xmm0			\n\t"\
		"subpd	0x10(%%rbx),%%xmm1			\n\t"\
		"\n\t"\
		"movaps	    (%%rcx),%%xmm4			\n\t"\
		"movaps	0x10(%%rcx),%%xmm5			\n\t"\
		"movaps	%%xmm4,%%xmm6				\n\t"\
		"movaps	%%xmm5,%%xmm7				\n\t"\
		"addpd	    (%%rdx),%%xmm6			\n\t"\
		"addpd	0x10(%%rdx),%%xmm7			\n\t"\
		"subpd	    (%%rdx),%%xmm4			\n\t"\
		"subpd	0x10(%%rdx),%%xmm5			\n\t"\
		"/* Copy t6r,i into main-array slot add6 */\n\t"\
		"movaps	%%xmm6,0xc0(%%rsi)			\n\t"\
		"movaps	%%xmm7,0xd0(%%rsi)			\n\t"\
		"/* Copy t7r,i into main-array slot add7 */\n\t"\
		"movaps	%%xmm4,0xe0(%%rsi)			\n\t"\
		"movaps	%%xmm5,0xf0(%%rsi)			\n\t"\
		"\n\t"\
		"/* Move outputs t5r,i into a(j1,j2+p5), first doing the addsub and mul by ISRT2: */\n\t"\
		"/* Move outputs t7r,i into a(j1,j2+p7), first doing the addsub and mul by ISRT2: */\n\t"\
		"\n\t"\
		"addpd	%%xmm2,%%xmm6				\n\t"\
		"addpd	%%xmm3,%%xmm7				\n\t"\
		"subpd	0xc0(%%rsi),%%xmm2			\n\t"\
		"subpd	0xd0(%%rsi),%%xmm3			\n\t"\
		"/* Move t6r,i into a(j1,j2+p6) */	\n\t"\
		"movaps	%%xmm2,0xc0(%%rsi)			\n\t"\
		"movaps	%%xmm3,0xd0(%%rsi)			\n\t"\
		"\n\t"\
		"movaps	%%xmm4,%%xmm2				\n\t"\
		"movaps	%%xmm5,%%xmm3				\n\t"\
		"addpd	%%xmm0,%%xmm5				\n\t"\
		"subpd	%%xmm3,%%xmm0				\n\t"\
		"addpd	%%xmm1,%%xmm4				\n\t"\
		"subpd	%%xmm2,%%xmm1				\n\t"\
		"\n\t"\
		"movaps	%%xmm5,%%xmm2				\n\t"\
		"movaps	%%xmm1,%%xmm3				\n\t"\
		"addpd	%%xmm1,%%xmm5				\n\t"\
		"movq	%[__isrt2],%%rdi			\n\t"\
		"movaps	(%%rdi),%%xmm1	/* isrt2 */	\n\t"\
		"subpd	%%xmm3,%%xmm2				\n\t"\
		"mulpd	%%xmm1,%%xmm5				\n\t"\
		"mulpd	%%xmm1,%%xmm2				\n\t"\
		"movaps	%%xmm0,%%xmm3				\n\t"\
		"\n\t"\
		"movaps	%%xmm5,0xa0(%%rsi)			\n\t"\
		"movaps	%%xmm4,%%xmm5				\n\t"\
		"addpd	%%xmm4,%%xmm0				\n\t"\
		"movaps	%%xmm2,0xb0(%%rsi)			\n\t"\
		"subpd	%%xmm5,%%xmm3				\n\t"\
		"mulpd	%%xmm1,%%xmm0				\n\t"\
		"mulpd	%%xmm1,%%xmm3				\n\t"\
		"movaps	%%xmm0,0xe0(%%rsi)			\n\t"\
		"movaps	%%xmm3,0xf0(%%rsi)			\n\t"\
		"\n\t"\
		"/**************** 1st of the 2 length-4 subtransforms... **************/\n\t"\
		"movq	%[__in0],%%rax	\n\t"\
		"movq	%[__in1],%%rbx	\n\t"\
		"movq	%[__in2],%%rcx	\n\t"\
		"movq	%[__in3],%%rdx	\n\t"\
		"\n\t"\
		"movaps	    (%%rax),%%xmm0			\n\t"\
		"movaps	0x10(%%rax),%%xmm1			\n\t"\
		"movaps	%%xmm0,%%xmm2				\n\t"\
		"movaps	%%xmm1,%%xmm3				\n\t"\
		"addpd	    (%%rbx),%%xmm2			\n\t"\
		"addpd	0x10(%%rbx),%%xmm3			\n\t"\
		"subpd	    (%%rbx),%%xmm0			\n\t"\
		"subpd	0x10(%%rbx),%%xmm1			\n\t"\
		"\n\t"\
		"/* Move   t4r,i into temp(j1+p0) in anticipation of final outputs (t0+t4)r,i which will go there: */\n\t"\
		"movaps	%%xmm6,    (%%rsi)			\n\t"\
		"movaps	%%xmm7,0x10(%%rsi)			\n\t"\
		"addpd	%%xmm6,%%xmm6				\n\t"\
		"addpd	%%xmm7,%%xmm7				\n\t"\
		"/* Move 2*t4r,i into temp(j1+p4) in anticipation of final outputs (t0-t4)r,i which will go there: */\n\t"\
		"movaps	%%xmm6,0x80(%%rsi)			\n\t"\
		"movaps	%%xmm7,0x90(%%rsi)			\n\t"\
		"\n\t"\
		"movaps	    (%%rcx),%%xmm4			\n\t"\
		"movaps	0x10(%%rcx),%%xmm5			\n\t"\
		"movaps	%%xmm4,%%xmm6				\n\t"\
		"movaps	%%xmm5,%%xmm7				\n\t"\
		"addpd	    (%%rdx),%%xmm6			\n\t"\
		"addpd	0x10(%%rdx),%%xmm7			\n\t"\
		"subpd	    (%%rdx),%%xmm4			\n\t"\
		"subpd	0x10(%%rdx),%%xmm5			\n\t"\
		"/* Copy t5,6 into temp(j1+p2) */	\n\t"\
		"movaps	%%xmm6,0x40(%%rsi)			\n\t"\
		"movaps	%%xmm7,0x50(%%rsi)			\n\t"\
		"addpd	%%xmm2,%%xmm6				\n\t"\
		"addpd	%%xmm3,%%xmm7				\n\t"\
		"subpd	0x40(%%rsi),%%xmm2			\n\t"\
		"subpd	0x50(%%rsi),%%xmm3			\n\t"\
		"\n\t"\
		"/* Compute and dump first 2 outputs now, in order to free up 2 registers: */\n\t"\
		"addpd	    (%%rsi),%%xmm6			\n\t"\
		"addpd	0x10(%%rsi),%%xmm7			\n\t"\
		"movaps	%%xmm6,    (%%rsi)			\n\t"\
		"movaps	%%xmm7,0x10(%%rsi)			\n\t"\
		"\n\t"\
		"subpd	0x80(%%rsi),%%xmm6			\n\t"\
		"subpd	0x90(%%rsi),%%xmm7			\n\t"\
		"movaps	%%xmm6,0x80(%%rsi)			\n\t"\
		"movaps	%%xmm7,0x90(%%rsi)			\n\t"\
		"\n\t"\
		"movaps	%%xmm4,%%xmm6				\n\t"\
		"movaps	%%xmm5,%%xmm7				\n\t"\
		"addpd	%%xmm0,%%xmm5				\n\t"\
		"subpd	%%xmm7,%%xmm0				\n\t"\
		"addpd	%%xmm1,%%xmm4				\n\t"\
		"subpd	%%xmm6,%%xmm1				\n\t"\
		"\n\t"\
		"movaps	%%xmm2,%%xmm6				\n\t"\
		"movaps	%%xmm3,%%xmm7				\n\t"\
		"addpd	0xd0(%%rsi),%%xmm2			\n\t"\
		"subpd	0xc0(%%rsi),%%xmm3			\n\t"\
		"subpd	0xd0(%%rsi),%%xmm6			\n\t"\
		"addpd	0xc0(%%rsi),%%xmm7			\n\t"\
		"movaps	%%xmm2,0xc0(%%rsi)			\n\t"\
		"movaps	%%xmm3,0xd0(%%rsi)			\n\t"\
		"movaps	%%xmm6,0x40(%%rsi)			\n\t"\
		"movaps	%%xmm7,0x50(%%rsi)			\n\t"\
		"\n\t"\
		"movaps	%%xmm5,%%xmm2				\n\t"\
		"movaps	%%xmm0,%%xmm6				\n\t"\
		"movaps	%%xmm1,%%xmm3				\n\t"\
		"movaps	%%xmm4,%%xmm7				\n\t"\
		"addpd	0xa0(%%rsi),%%xmm5			\n\t"\
		"subpd	0xf0(%%rsi),%%xmm0			\n\t"\
		"subpd	0xb0(%%rsi),%%xmm1			\n\t"\
		"subpd	0xe0(%%rsi),%%xmm4			\n\t"\
		"subpd	0xa0(%%rsi),%%xmm2			\n\t"\
		"addpd	0xf0(%%rsi),%%xmm6			\n\t"\
		"addpd	0xb0(%%rsi),%%xmm3			\n\t"\
		"addpd	0xe0(%%rsi),%%xmm7			\n\t"\
		"movaps	%%xmm5,0xe0(%%rsi)			\n\t"\
		"movaps	%%xmm0,0xa0(%%rsi)			\n\t"\
		"movaps	%%xmm1,0xf0(%%rsi)			\n\t"\
		"movaps	%%xmm4,0xb0(%%rsi)			\n\t"\
		"movaps	%%xmm2,0x60(%%rsi)			\n\t"\
		"movaps	%%xmm6,0x20(%%rsi)			\n\t"\
		"movaps	%%xmm3,0x70(%%rsi)			\n\t"\
		"movaps	%%xmm7,0x30(%%rsi)			\n\t"\
		"\n\t"\
		:					/* outputs: none */\
		: [__in0] "m" (Xin0)	/* All inputs from memory addresses here */\
		 ,[__in1] "m" (Xin1)\
		 ,[__in2] "m" (Xin2)\
		 ,[__in3] "m" (Xin3)\
		 ,[__in4] "m" (Xin4)\
		 ,[__in5] "m" (Xin5)\
		 ,[__in6] "m" (Xin6)\
		 ,[__in7] "m" (Xin7)\
		 ,[__out] "m" (Xout)\
		 ,[__isrt2] "m" (Xisrt2)\
		: "rax","rbx","rcx","rdx","rdi","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
	);\
	}

	#define SSE2_RADIX_05_DFT_0TWIDDLE(Xi0,Xi1,Xi2,Xi3,Xi4, Xcc1, Xo0,Xo1,Xo2,Xo3,Xo4)\
	{\
	__asm__ volatile (\
		"movq	%[__i0],%%rsi		\n\t"\
		"movq	%[__i1],%%rax		\n\t"\
		"movq	%[__i2],%%rbx		\n\t"\
		"movq	%[__i3],%%rcx		\n\t"\
		"movq	%[__i4],%%rdx		\n\t"\
		"movq	%[__o0],%%rdi		\n\t"\
		"movaps	    (%%rax),%%xmm0	\n\t"\
		"movaps	0x10(%%rax),%%xmm1	\n\t"\
		"movaps	    (%%rbx),%%xmm2	\n\t"\
		"movaps	0x10(%%rbx),%%xmm3	\n\t"\
		"movaps	    (%%rcx),%%xmm4	\n\t"\
		"movaps	0x10(%%rcx),%%xmm5	\n\t"\
		"movaps	    (%%rdx),%%xmm6	\n\t"\
		"movaps	0x10(%%rdx),%%xmm7	\n\t"\
		"subpd	%%xmm6,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm1		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm0,%%xmm6		\n\t"\
		"addpd	%%xmm1,%%xmm7		\n\t"\
		"subpd	%%xmm4,%%xmm2		\n\t"\
		"subpd	%%xmm5,%%xmm3		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm2,%%xmm4		\n\t"\
		"addpd	%%xmm3,%%xmm5		\n\t"\
	"movq	%[__cc1],%%rax		\n\t"\
		"subpd	%%xmm4,%%xmm6		\n\t"\
		"subpd	%%xmm5,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm5		\n\t"\
		"addpd	    (%%rsi),%%xmm4	\n\t"\
		"addpd	0x10(%%rsi),%%xmm5	\n\t"\
		"movaps	%%xmm4,    (%%rdi)	\n\t"\
		"movaps	%%xmm5,0x10(%%rdi)	\n\t"\
		"mulpd	0x10(%%rax),%%xmm6	\n\t"\
		"mulpd	0x10(%%rax),%%xmm7	\n\t"\
		"subpd	     (%%rsi),%%xmm4	\n\t"\
		"subpd	0x010(%%rsi),%%xmm5	\n\t"\
		"mulpd	    (%%rax),%%xmm4	\n\t"\
		"mulpd	    (%%rax),%%xmm5	\n\t"\
		"addpd	     (%%rdi),%%xmm4	\n\t"\
		"addpd	0x010(%%rdi),%%xmm5	\n\t"\
		"subpd	%%xmm6,%%xmm4		\n\t"\
		"subpd	%%xmm7,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm6		\n\t"\
		"addpd	%%xmm5,%%xmm7		\n\t"\
		"movaps	%%xmm4,    (%%rsi)	\n\t"\
		"movaps	%%xmm5,0x10(%%rsi)	\n\t"\
		"movaps	%%xmm0,%%xmm4		\n\t"\
		"movaps	%%xmm1,%%xmm5		\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"mulpd	0x20(%%rax),%%xmm0	\n\t"\
		"mulpd	0x20(%%rax),%%xmm1	\n\t"\
		"mulpd	0x30(%%rax),%%xmm2	\n\t"\
		"mulpd	0x30(%%rax),%%xmm3	\n\t"\
		"mulpd	0x40(%%rax),%%xmm4	\n\t"\
		"mulpd	0x40(%%rax),%%xmm5	\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"movaps	    (%%rsi),%%xmm4	\n\t"\
		"movaps	0x10(%%rsi),%%xmm5	\n\t"\
		"movq	%[__o1],%%rax		\n\t"\
		"movq	%[__o4],%%rdx		\n\t"\
		"subpd	%%xmm3,%%xmm6		\n\t"\
		"subpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"movaps	%%xmm6,    (%%rax)	\n\t"\
		"movaps	%%xmm7,0x10(%%rdx)	\n\t"\
		"addpd	%%xmm6,%%xmm3		\n\t"\
		"addpd	%%xmm7,%%xmm2		\n\t"\
		"movaps	%%xmm3,    (%%rdx)	\n\t"\
		"movaps	%%xmm2,0x10(%%rax)	\n\t"\
		"movq	%[__o2],%%rbx		\n\t"\
		"movq	%[__o3],%%rcx		\n\t"\
		"subpd	%%xmm1,%%xmm4		\n\t"\
		"subpd	%%xmm0,%%xmm5		\n\t"\
		"addpd	%%xmm1,%%xmm1		\n\t"\
		"addpd	%%xmm0,%%xmm0		\n\t"\
		"movaps	%%xmm4,    (%%rbx)	\n\t"\
		"movaps	%%xmm5,0x10(%%rcx)	\n\t"\
		"addpd	%%xmm4,%%xmm1		\n\t"\
		"addpd	%%xmm5,%%xmm0		\n\t"\
		"movaps	%%xmm1,    (%%rcx)	\n\t"\
		"movaps	%%xmm0,0x10(%%rbx)	\n\t"\
		"							\n\t"\
		:					/* outputs: none */\
		: [__i0] "m" (Xi0)	/* All inputs from memory addresses here */\
		 ,[__i1] "m" (Xi1)\
		 ,[__i2] "m" (Xi2)\
		 ,[__i3] "m" (Xi3)\
		 ,[__i4] "m" (Xi4)\
		 ,[__cc1] "m" (Xcc1)\
		 ,[__o0] "m" (Xo0)\
		 ,[__o1] "m" (Xo1)\
		 ,[__o2] "m" (Xo2)\
		 ,[__o3] "m" (Xo3)\
		 ,[__o4] "m" (Xo4)\
		: "rax","rbx","rcx","rdx","rdi","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
	);\
	}

	/* Twiddleless version of SSE2_RADIX8_DIF_TWIDDLE. Inputs enter in memory locations __r0 + [__i1,__i2,__i3,__i4,__i5,__i6,__i7],;
	where r0 is a memory address and the i's are LITERAL [BYTE] OFFSETS. Outputs go into memory locations __o0,__o1,__o2,__o3,__o4,__o5,__o6,__o7, assumed disjoint with inputs:\
	*/
	#define SSE2_RADIX8_DIF_0TWIDDLE(Xr0, Xi1,Xi2,Xi3,Xi4,Xi5,Xi6,Xi7, Xo0,Xo1,Xo2,Xo3,Xo4,Xo5,Xo6,Xo7, Xisrt2)\
	{\
	__asm__ volatile (\
		"movq	%[__r0],%%rax	/* i0 = r00 */\n\t"\
		"movq	%[__i2],%%rbx	/* i2 */	\n\t"\
		"movq	%[__i4],%%rcx	/* i4 */	\n\t"\
		"movq	%[__i6],%%rdx	/* i6 */	\n\t"\
		"addq	%%rax,%%rbx					\n\t"\
		"addq	%%rax,%%rcx					\n\t"\
		"addq	%%rax,%%rdx					\n\t"\
		"\n\t"\
		"/* Do the p0,p4 combo: */			\n\t"\
		"\n\t"\
		"movaps	    (%%rax),%%xmm0			\n\t"\
		"movaps	    (%%rcx),%%xmm4			\n\t"\
		"movaps	0x10(%%rax),%%xmm1			\n\t"\
		"movaps	0x10(%%rcx),%%xmm5			\n\t"\
		"movaps	%%xmm0,%%xmm2				\n\t"\
		"movaps	%%xmm1,%%xmm3				\n\t"\
		"\n\t"\
		"addpd	%%xmm4,%%xmm0				\n\t"\
		"addpd	%%xmm5,%%xmm1				\n\t"\
		"subpd	%%xmm4,%%xmm2				\n\t"\
		"subpd	%%xmm5,%%xmm3				\n\t"\
		"\n\t"\
		"/* Do the p2,6 combo: */			\n\t"\
		"\n\t"\
		"movaps	    (%%rbx),%%xmm4			\n\t"\
		"movaps	0x10(%%rbx),%%xmm5			\n\t"\
		"movaps	%%xmm4,%%xmm6				\n\t"\
		"movaps	%%xmm5,%%xmm7				\n\t"\
		"\n\t"\
		"addpd	     (%%rdx),%%xmm4			\n\t"\
		"addpd	0x010(%%rdx),%%xmm5			\n\t"\
		"subpd	     (%%rdx),%%xmm6			\n\t"\
		"subpd	0x010(%%rdx),%%xmm7			\n\t"\
		"\n\t"\
		"subpd	%%xmm4,%%xmm0				\n\t"\
		"subpd	%%xmm7,%%xmm2				\n\t"\
		"subpd	%%xmm5,%%xmm1				\n\t"\
		"subpd	%%xmm6,%%xmm3				\n\t"\
		"movaps	%%xmm0,     (%%rcx)			\n\t"\
		"movaps	%%xmm2,     (%%rbx)			\n\t"\
		"movaps	%%xmm1,0x010(%%rcx)			\n\t"\
		"movaps	%%xmm3,0x010(%%rdx)			\n\t"\
		"addpd	%%xmm4,%%xmm4				\n\t"\
		"addpd	%%xmm7,%%xmm7				\n\t"\
		"addpd	%%xmm5,%%xmm5				\n\t"\
		"addpd	%%xmm6,%%xmm6				\n\t"\
		"addpd	%%xmm0,%%xmm4				\n\t"\
		"addpd	%%xmm2,%%xmm7				\n\t"\
		"addpd	%%xmm1,%%xmm5				\n\t"\
		"addpd	%%xmm3,%%xmm6				\n\t"\
		"movaps	%%xmm4,     (%%rax)			\n\t"\
		"movaps	%%xmm7,     (%%rdx)			\n\t"\
		"movaps	%%xmm5,0x010(%%rax)			\n\t"\
		"movaps	%%xmm6,0x010(%%rbx)			\n\t"\
		"\n\t"\
		"/* Inline of SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO_B(r8,c1): Do 2nd of 2 radix-4 subtransforms, but now keep outputs in registers: */\n\t"\
		"\n\t"\
		"movq	%[__i3],%%rbx		/* i3 */	\n\t"\
		"movq	%[__i5],%%rcx		/* i5 */	\n\t"\
		"movq	%[__i7],%%rdx		/* i7 */	\n\t"\
		"addq	%%rax,%%rbx					\n\t"\
		"addq	%%rax,%%rcx					\n\t"\
		"addq	%%rax,%%rdx					\n\t"\
		"addq	%[__i1],%%rax		/* i1 */	\n\t"\
		"\n\t"\
		"/* Do the p1,p5 combo: */			\n\t"\
		"\n\t"\
		"movaps	    (%%rax),%%xmm0			\n\t"\
		"movaps	    (%%rcx),%%xmm4			\n\t"\
		"movaps	0x10(%%rax),%%xmm1			\n\t"\
		"movaps	0x10(%%rcx),%%xmm5			\n\t"\
		"movaps	%%xmm0,%%xmm2				\n\t"\
		"movaps	%%xmm1,%%xmm3				\n\t"\
		"\n\t"\
		"addpd	%%xmm4,%%xmm0				\n\t"\
		"addpd	%%xmm5,%%xmm1				\n\t"\
		"subpd	%%xmm4,%%xmm2				\n\t"\
		"subpd	%%xmm5,%%xmm3				\n\t"\
		"\n\t"\
		"/* Do the p3,7 combo: */			\n\t"\
		"\n\t"\
		"movaps	    (%%rbx),%%xmm4			\n\t"\
		"movaps	0x10(%%rbx),%%xmm5			\n\t"\
		"movaps	%%xmm4,%%xmm6				\n\t"\
		"movaps	%%xmm5,%%xmm7				\n\t"\
		"\n\t"\
		"subpd	    (%%rdx),%%xmm4			\n\t"\
		"subpd	0x10(%%rdx),%%xmm5			\n\t"\
		"addpd	    (%%rdx),%%xmm6			\n\t"\
		"addpd	0x10(%%rdx),%%xmm7			\n\t"\
		"\n\t"\
		"/* Finish radix-4 butterfly and store just the 1st of the 8 outputs into output-array slots: */\n\t"\
		"\n\t"\
		"subpd	%%xmm6,%%xmm0				\n\t"\
		"subpd	%%xmm7,%%xmm1				\n\t"\
		"subpd	%%xmm5,%%xmm2				\n\t"\
		"subpd	%%xmm4,%%xmm3				\n\t"\
		"addpd	%%xmm6,%%xmm6				\n\t"\
		"addpd	%%xmm5,%%xmm5				\n\t"\
		"addpd	%%xmm7,%%xmm7				\n\t"\
		"addpd	%%xmm4,%%xmm4				\n\t"\
		"addpd	%%xmm0,%%xmm6				\n\t"\
		"addpd	%%xmm2,%%xmm5				\n\t"\
		"addpd	%%xmm1,%%xmm7				\n\t"\
		"addpd	%%xmm3,%%xmm4				\n\t"\
		"movaps	%%xmm6,    (%%rax)			\n\t"\
		"movaps	%%xmm7,0x10(%%rax)			\n\t"\
		"\n\t"\
		"movaps	%%xmm2,%%xmm6				\n\t"\
		"movaps	%%xmm5,%%xmm7				\n\t"\
		"subpd	%%xmm4,%%xmm2				\n\t"\
		"subpd	%%xmm3,%%xmm5				\n\t"\
		"addpd	%%xmm6,%%xmm4				\n\t"\
		"addpd	%%xmm7,%%xmm3				\n\t"\
		"\n\t"\
		"movq	%[__isrt2],%%rsi			\n\t"\
		"movaps	(%%rsi),%%xmm6	/* isrt2 */	\n\t"\
		"mulpd	%%xmm6,%%xmm2				\n\t"\
		"mulpd	%%xmm6,%%xmm5				\n\t"\
		"mulpd	%%xmm6,%%xmm4				\n\t"\
		"mulpd	%%xmm6,%%xmm3				\n\t"\
		"\n\t"\
		"/* Inline the contents of SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_A(r0): Combine the 2 radix-4 subtransforms and write outputs to temp array: */\n\t"\
		"/***** t0,1,2,3,4,5,6,7 in __i0,2,4,6 (eax,ebx,ecx,edx) *****/\n\t"\
		"/***** t8,9,a,b,c,d,e,f in __i1,xmm2,4,0,1,5,3; xmm6,7 FREE */\n\t"\
		"movq	%[__r0],%%rdi			\n\t"\
		"movq	%[__r0],%%rsi			\n\t"\
		"addq	$0x10,%%rsi				\n\t"\
		"movq	%[__o4],%%rax	\n\t"\
		"movq	%[__o5],%%rbx	\n\t"\
		"movq	%[__o6],%%rcx	\n\t"\
		"movq	%[__o7],%%rdx	\n\t"\
		"\n\t"\
		"movaps	%c[__i2](%%rdi),%%xmm6		\n\t"\
		"movaps	%c[__i2](%%rsi),%%xmm7		\n\t"\
		"subpd   %%xmm2,%%xmm6				\n\t"\
		"subpd   %%xmm4,%%xmm7				\n\t"\
		"addpd   %%xmm2,%%xmm2				\n\t"\
		"addpd   %%xmm4,%%xmm4				\n\t"\
		"addpd   %%xmm6,%%xmm2				\n\t"\
		"addpd   %%xmm7,%%xmm4				\n\t"\
		"\n\t"\
		"movaps	%%xmm6,    (%%rbx)	/* o5r */\n\t"\
		"movaps	%%xmm7,0x10(%%rbx)	/* o5i */\n\t"\
		"movaps	%%xmm2,    (%%rax)	/* o4r */\n\t"\
		"movaps	%%xmm4,0x10(%%rax)	/* o4i */\n\t"\
		"\n\t"\
		"/* Can't do 2nd set side-by-side with 1st, since only 2 registers (xmm6,7) available for temp storage: */\n\t"\
		"movaps	%c[__i6](%%rdi),%%xmm6		\n\t"\
		"movaps	%c[__i6](%%rsi),%%xmm7		\n\t"\
		"subpd   %%xmm3,%%xmm6			\n\t"\
		"subpd   %%xmm5,%%xmm7			\n\t"\
		"addpd   %%xmm3,%%xmm3			\n\t"\
		"addpd   %%xmm5,%%xmm5			\n\t"\
		"addpd   %%xmm6,%%xmm3			\n\t"\
		"addpd   %%xmm7,%%xmm5			\n\t"\
		"\n\t"\
		"movaps	%%xmm6,    (%%rcx)	/* o6r */\n\t"\
		"movaps	%%xmm7,0x10(%%rdx)	/* o7i */\n\t"\
		"movaps	%%xmm3,    (%%rdx)	/* o7r */\n\t"\
		"movaps	%%xmm5,0x10(%%rcx)	/* o6i */\n\t"\
		"\n\t"\
		"movq	%[__o0],%%rax	\n\t"\
		"movq	%[__o1],%%rbx	\n\t"\
		"movq	%[__o2],%%rcx	\n\t"\
		"movq	%[__o3],%%rdx	\n\t"\
		"\n\t"\
		"movaps	    (%%rdi),%%xmm6			\n\t"\
		"movaps	%c[__i4](%%rdi),%%xmm4		\n\t"\
		"movaps	    (%%rsi),%%xmm7			\n\t"\
		"movaps	%c[__i4](%%rsi),%%xmm5		\n\t"\
		"movaps	%c[__i1](%%rdi),%%xmm2		\n\t"\
		"movaps	%c[__i1](%%rsi),%%xmm3		\n\t"\
		"subpd   %%xmm2,%%xmm6				\n\t"\
		"subpd   %%xmm1,%%xmm4				\n\t"\
		"subpd   %%xmm3,%%xmm7				\n\t"\
		"subpd   %%xmm0,%%xmm5				\n\t"\
		"addpd   %%xmm2,%%xmm2				\n\t"\
		"addpd   %%xmm1,%%xmm1				\n\t"\
		"addpd   %%xmm3,%%xmm3				\n\t"\
		"addpd   %%xmm0,%%xmm0				\n\t"\
		"addpd   %%xmm6,%%xmm2				\n\t"\
		"addpd   %%xmm4,%%xmm1				\n\t"\
		"addpd   %%xmm7,%%xmm3				\n\t"\
		"addpd   %%xmm5,%%xmm0				\n\t"\
		"\n\t"\
		"movaps	%%xmm6,    (%%rbx)	/* o1r */\n\t"\
		"movaps	%%xmm4,    (%%rcx)	/* o2r */\n\t"\
		"movaps	%%xmm7,0x10(%%rbx)	/* o1i */\n\t"\
		"movaps	%%xmm5,0x10(%%rdx)	/* o3i */\n\t"\
		"movaps	%%xmm2,    (%%rax)	/* o0r */\n\t"\
		"movaps	%%xmm1,    (%%rdx)	/* o3r */\n\t"\
		"movaps	%%xmm3,0x10(%%rax)	/* o0r */\n\t"\
		"movaps	%%xmm0,0x10(%%rcx)	/* o2i */\n\t"\
		"\n\t"\
		:					/* outputs: none */\
		: [__r0] "m" (Xr0)	/* All inputs from memory addresses here */\
		 ,[__i1] "e" (Xi1)\
		 ,[__i2] "e" (Xi2)\
		 ,[__i3] "e" (Xi3)\
		 ,[__i4] "e" (Xi4)\
		 ,[__i5] "e" (Xi5)\
		 ,[__i6] "e" (Xi6)\
		 ,[__i7] "e" (Xi7)\
		 ,[__o0] "m" (Xo0)\
		 ,[__o1] "m" (Xo1)\
		 ,[__o2] "m" (Xo2)\
		 ,[__o3] "m" (Xo3)\
		 ,[__o4] "m" (Xo4)\
		 ,[__o5] "m" (Xo5)\
		 ,[__o6] "m" (Xo6)\
		 ,[__o7] "m" (Xo7)\
		 ,[__isrt2] "m" (Xisrt2)\
		: "rax","rbx","rcx","rdx","rdi","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
	);\
	}

#endif	/* sse2_macro_gcc_h_included */

