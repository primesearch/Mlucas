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

	#define SSE2_RADIX_03_DFT(Xi0,Xi1,Xi2, Xcc1, Xo0,Xo1,Xo2)\
	{\
	__asm__ volatile (\
			"movq	%[__i0],%%rax		\n\t"\
			"movq	%[__i1],%%rbx		\n\t"\
			"movq	%[__i2],%%rcx		\n\t"\
			"movq	%[__cc1],%%rdx		\n\t"\
			"\n\t"\
			"movaps	    (%%rbx),%%xmm2	\n\t"\
			"movaps	0x10(%%rbx),%%xmm3	\n\t"\
			"movaps	    (%%rax),%%xmm0	\n\t"\
			"movaps	0x10(%%rax),%%xmm1	\n\t"\
			"movaps	    (%%rcx),%%xmm6	\n\t"\
			"movaps	0x10(%%rcx),%%xmm7	\n\t"\
			"movaps	%%xmm2,%%xmm4		\n\t"\
			"movaps	%%xmm3,%%xmm5		\n\t"\
			"\n\t"\
			"movq	%[__o0],%%rax		\n\t"\
			"movq	%[__o1],%%rbx		\n\t"\
			"movq	%[__o2],%%rcx		\n\t"\
			"addpd	%%xmm6,%%xmm2		\n\t"\
			"addpd	%%xmm7,%%xmm3		\n\t"\
			"subpd	%%xmm6,%%xmm4		\n\t"\
			"subpd	%%xmm7,%%xmm5		\n\t"\
			"addpd	%%xmm2,%%xmm0		\n\t"\
			"addpd	%%xmm3,%%xmm1		\n\t"\
			"movaps	    (%%rdx),%%xmm6	\n\t"\
			"movaps	0x10(%%rdx),%%xmm7	\n\t"\
			"movaps	%%xmm0,     (%%rax)	\n\t"\
			"movaps	%%xmm1,0x010(%%rax)	\n\t"\
			"\n\t"\
			"mulpd	%%xmm6,%%xmm2		\n\t"\
			"mulpd	%%xmm6,%%xmm3		\n\t"\
			"mulpd	%%xmm7,%%xmm4		\n\t"\
			"mulpd	%%xmm7,%%xmm5		\n\t"\
			"addpd	%%xmm0,%%xmm2		\n\t"\
			"addpd	%%xmm1,%%xmm3		\n\t"\
			"\n\t"\
			"movaps	%%xmm2,%%xmm0		\n\t"\
			"movaps	%%xmm3,%%xmm1		\n\t"\
			"\n\t"\
			"subpd	%%xmm5,%%xmm2		\n\t"\
			"addpd	%%xmm4,%%xmm3		\n\t"\
			"addpd	%%xmm5,%%xmm0		\n\t"\
			"subpd	%%xmm4,%%xmm1		\n\t"\
			"\n\t"\
			"movaps	%%xmm2,     (%%rbx)	\n\t"\
			"movaps	%%xmm3,0x010(%%rbx)	\n\t"\
			"movaps	%%xmm0,     (%%rcx)	\n\t"\
			"movaps	%%xmm1,0x010(%%rcx)	\n\t"\
		:					/* outputs: none */\
		: [__i0] "m" (Xi0)	/* All inputs from memory addresses here */\
		 ,[__i1] "m" (Xi1)\
		 ,[__i2] "m" (Xi2)\
		 ,[__cc1] "m" (Xcc1)\
		 ,[__o0] "m" (Xo0)\
		 ,[__o1] "m" (Xo1)\
		 ,[__o2] "m" (Xo2)\
		: "rax","rbx","rcx","rdx","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
	);\
	}

	#define SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(Xadd0, Xadd1, Xadd2, Xadd3, Xtmp, Xstride)\
	{\
	__asm__ volatile (\
		"movq	%[__tmp]   ,%%rax	\n\t"\
		"movq	%[__stride],%%rsi	\n\t"\
		"movq	%%rax,%%rbx			\n\t"\
		"addq	%%rsi,%%rbx			/* add_in1  */\n\t"\
		"shlq	$1,%%rsi			/* stride*2 */\n\t"\
		"movaps	    (%%rax),%%xmm0	\n\t"\
		"movaps	    (%%rbx),%%xmm2	\n\t"\
		"movaps	0x10(%%rax),%%xmm1	\n\t"\
		"movaps	0x10(%%rbx),%%xmm3	\n\t"\
		"movaps	    (%%rax),%%xmm4	\n\t"\
		"movaps	    (%%rbx),%%xmm6	\n\t"\
		"movaps	0x10(%%rax),%%xmm5	\n\t"\
		"movaps	0x10(%%rbx),%%xmm7	\n\t"\
		"addq	%%rsi,%%rax			/* add_in2  */\n\t"\
		"addq	%%rsi,%%rbx			/* add_in3  */\n\t"\
		"addpd	    (%%rax),%%xmm0	\n\t"\
		"addpd	    (%%rbx),%%xmm2	\n\t"\
		"addpd	0x10(%%rax),%%xmm1	\n\t"\
		"addpd	0x10(%%rbx),%%xmm3	\n\t"\
		"subpd	    (%%rax),%%xmm4	\n\t"\
		"subpd	    (%%rbx),%%xmm6	\n\t"\
		"subpd	0x10(%%rax),%%xmm5	\n\t"\
		"subpd	0x10(%%rbx),%%xmm7	\n\t"\
		"/* Finish radix-4 butterfly and store results into main-array slots: */\n\t"\
		"movq	%[__add0],%%rax		\n\t"\
		"movq	%[__add1],%%rbx		\n\t"\
		"movq	%[__add2],%%rcx		\n\t"\
		"movq	%[__add3],%%rdx		\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm4		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"subpd	%%xmm6,%%xmm5		\n\t"\
		"movaps	%%xmm0,     (%%rbx)	\n\t"\
		"movaps	%%xmm4,     (%%rcx)	\n\t"\
		"movaps	%%xmm1,0x010(%%rbx)	\n\t"\
		"movaps	%%xmm5,0x010(%%rdx)	\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm4,%%xmm7		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"addpd	%%xmm5,%%xmm6		\n\t"\
		"movaps	%%xmm2,     (%%rax)	\n\t"\
		"movaps	%%xmm7,     (%%rdx)	\n\t"\
		"movaps	%%xmm3,0x010(%%rax)	\n\t"\
		"movaps	%%xmm6,0x010(%%rcx)	\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__add1] "m" (Xadd1)\
		 ,[__add2] "m" (Xadd2)\
		 ,[__add3] "m" (Xadd3)\
		 ,[__tmp] "m" (Xtmp)\
		 ,[__stride] "e" (Xstride)\
		: "rax","rbx","rcx","rdx","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
	);\
	}

	#define SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(Xadd0, Xadd1, Xadd2, Xadd3, Xtmp, Xstride)\
	{\
	__asm__ volatile (\
		"movq	%[__add0],%%rax		\n\t"\
		"movq	%[__add1],%%rbx		\n\t"\
		"movq	%[__add2],%%rcx		\n\t"\
		"movq	%[__add3],%%rdx		\n\t"\
		"movaps	    (%%rax),%%xmm0	\n\t"\
		"movaps	    (%%rcx),%%xmm4	\n\t"\
		"movaps	0x10(%%rax),%%xmm1	\n\t"\
		"movaps	0x10(%%rcx),%%xmm5	\n\t"\
		"movaps	%%xmm0,%%xmm2			\n\t"\
		"movaps	%%xmm4,%%xmm6			\n\t"\
		"movaps	%%xmm1,%%xmm3			\n\t"\
		"movaps	%%xmm5,%%xmm7			\n\t"\
		"movq	%[__tmp]   ,%%rax	\n\t"\
		"movq	%[__stride],%%rcx	\n\t"\
		"addpd	    (%%rbx),%%xmm0	\n\t"\
		"addpd	    (%%rdx),%%xmm4	\n\t"\
		"addpd	0x10(%%rbx),%%xmm1	\n\t"\
		"addpd	0x10(%%rdx),%%xmm5	\n\t"\
		"subpd	    (%%rbx),%%xmm2	\n\t"\
		"subpd	    (%%rdx),%%xmm6	\n\t"\
		"subpd	0x10(%%rbx),%%xmm3	\n\t"\
		"subpd	0x10(%%rdx),%%xmm7	\n\t"\
		"movq	%%rax,%%rbx			\n\t"\
		"addq	%%rcx,%%rbx			\n\t"\
		"movq	%%rbx,%%rdx			\n\t"\
		"addq	%%rcx,%%rcx			\n\t"\
		"addq	%%rcx,%%rdx			\n\t"\
		"addq	%%rax,%%rcx			\n\t"\
		"/* Finish radix-4 butterfly and store results into temp-array slots: */\n\t"\
		"subpd	%%xmm4,%%xmm0			\n\t"\
		"subpd	%%xmm7,%%xmm2			\n\t"\
		"subpd	%%xmm5,%%xmm1			\n\t"\
		"subpd	%%xmm6,%%xmm3			\n\t"\
		"movaps	%%xmm0,     (%%rcx)	\n\t"\
		"movaps	%%xmm2,     (%%rdx)	\n\t"\
		"movaps	%%xmm1,0x010(%%rcx)	\n\t"\
		"movaps	%%xmm3,0x010(%%rbx)	\n\t"\
		"addpd	%%xmm4,%%xmm4			\n\t"\
		"addpd	%%xmm7,%%xmm7			\n\t"\
		"addpd	%%xmm5,%%xmm5			\n\t"\
		"addpd	%%xmm6,%%xmm6			\n\t"\
		"addpd	%%xmm0,%%xmm4			\n\t"\
		"addpd	%%xmm2,%%xmm7			\n\t"\
		"addpd	%%xmm1,%%xmm5			\n\t"\
		"addpd	%%xmm3,%%xmm6			\n\t"\
		"movaps	%%xmm4,     (%%rax)	\n\t"\
		"movaps	%%xmm7,     (%%rbx)	\n\t"\
		"movaps	%%xmm5,0x010(%%rax)	\n\t"\
		"movaps	%%xmm6,0x010(%%rdx)	\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__add1] "m" (Xadd1)\
		 ,[__add2] "m" (Xadd2)\
		 ,[__add3] "m" (Xadd3)\
		 ,[__tmp] "m" (Xtmp)\
		 ,[__stride] "e" (Xstride)\
		: "rax","rbx","rcx","rdx","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
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

	/* 16-xmm-register version does 2 of the above side-by-side: */
	#define SSE2_RADIX_05_DFT_0TWID_X2(Xcc1, Xi0,Xi1,Xi2,Xi3,Xi4, Xo0,Xo1,Xo2,Xo3,Xo4, Yi0,Yi1,Yi2,Yi3,Yi4, Yo0,Yo1,Yo2,Yo3,Yo4)\
	{\
	__asm__ volatile (\
		"movq	%[__i0],%%rsi		\n\t"\
		"movq	%[__i1],%%rax		\n\t"\
		"movq	%[__i2],%%rbx		\n\t"\
		"movq	%[__i3],%%rcx		\n\t"\
		"movq	%[__i4],%%rdx		\n\t"\
		"movq	%[__o0],%%rdi		\n\t"\
		"movaps	    (%%rax),%%xmm0	\n\t		movq	%[__i5],%%r10		\n\t"\
		"movaps	0x10(%%rax),%%xmm1	\n\t		movq	%[__i6],%%r11		\n\t"\
		"movaps	    (%%rbx),%%xmm2	\n\t		movq	%[__i7],%%r12		\n\t"\
		"movaps	0x10(%%rbx),%%xmm3	\n\t		movq	%[__i8],%%r13		\n\t"\
		"movaps	    (%%rcx),%%xmm4	\n\t		movq	%[__i9],%%r14		\n\t"\
		"movaps	0x10(%%rcx),%%xmm5	\n\t		movq	%[__o5],%%r15		\n\t"\
		"movaps	    (%%rdx),%%xmm6	\n\t		movaps	    (%%r11),%%xmm8 	\n\t"\
		"movaps	0x10(%%rdx),%%xmm7	\n\t		movaps	0x10(%%r11),%%xmm9 	\n\t"\
		"subpd	%%xmm6,%%xmm0		\n\t		movaps	    (%%r12),%%xmm10	\n\t"\
		"subpd	%%xmm7,%%xmm1		\n\t		movaps	0x10(%%r12),%%xmm11	\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t		movaps	    (%%r13),%%xmm12	\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t		movaps	0x10(%%r13),%%xmm13	\n\t"\
		"addpd	%%xmm0,%%xmm6		\n\t		movaps	    (%%r14),%%xmm14	\n\t"\
		"addpd	%%xmm1,%%xmm7		\n\t		movaps	0x10(%%r14),%%xmm15	\n\t"\
		"subpd	%%xmm4,%%xmm2		\n\t		subpd	%%xmm14,%%xmm8 		\n\t"\
		"subpd	%%xmm5,%%xmm3		\n\t		subpd	%%xmm15,%%xmm9 		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t		addpd	%%xmm14,%%xmm14		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t		addpd	%%xmm15,%%xmm15		\n\t"\
		"addpd	%%xmm2,%%xmm4		\n\t		addpd	%%xmm8 ,%%xmm14		\n\t"\
		"addpd	%%xmm3,%%xmm5		\n\t		addpd	%%xmm9 ,%%xmm15		\n\t"\
		"movq	%[__cc1],%%rax		/* Sincos data shared by both streams */\n\t"\
		"subpd	%%xmm4,%%xmm6		\n\t		subpd	%%xmm12,%%xmm10		\n\t"\
		"subpd	%%xmm5,%%xmm7		\n\t		subpd	%%xmm13,%%xmm11		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t		addpd	%%xmm12,%%xmm12		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t		addpd	%%xmm13,%%xmm13		\n\t"\
		"addpd	%%xmm6,%%xmm4		\n\t		addpd	%%xmm10,%%xmm12		\n\t"\
		"addpd	%%xmm7,%%xmm5		\n\t		addpd	%%xmm11,%%xmm13		\n\t"\
		"addpd	    (%%rsi),%%xmm4	\n\t		subpd	%%xmm12,%%xmm14		\n\t"\
		"addpd	0x10(%%rsi),%%xmm5	\n\t		subpd	%%xmm13,%%xmm15		\n\t"\
		"movaps	%%xmm4,    (%%rdi)	\n\t		addpd	%%xmm12,%%xmm12		\n\t"\
		"movaps	%%xmm5,0x10(%%rdi)	\n\t		addpd	%%xmm13,%%xmm13		\n\t"\
		"mulpd	0x10(%%rax),%%xmm6	\n\t		addpd	%%xmm14,%%xmm12		\n\t"\
		"mulpd	0x10(%%rax),%%xmm7	\n\t		addpd	%%xmm15,%%xmm13		\n\t"\
		"subpd	    (%%rsi),%%xmm4	\n\t		addpd	    (%%r10),%%xmm12	\n\t"\
		"subpd	0x10(%%rsi),%%xmm5	\n\t		addpd	0x10(%%r10),%%xmm13	\n\t"\
		"mulpd	    (%%rax),%%xmm4	\n\t		movaps	%%xmm12,    (%%r15)	\n\t"\
		"mulpd	    (%%rax),%%xmm5	\n\t		movaps	%%xmm13,0x10(%%r15)	\n\t"\
		"addpd	    (%%rdi),%%xmm4	\n\t		mulpd	0x10(%%rax),%%xmm14	\n\t"\
		"addpd	0x10(%%rdi),%%xmm5	\n\t		mulpd	0x10(%%rax),%%xmm15	\n\t"\
		"subpd	%%xmm6,%%xmm4		\n\t		subpd	    (%%r10),%%xmm12	\n\t"\
		"subpd	%%xmm7,%%xmm5		\n\t		subpd	0x10(%%r10),%%xmm13	\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t		mulpd	    (%%rax),%%xmm12	\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t		mulpd	    (%%rax),%%xmm13	\n\t"\
		"addpd	%%xmm4,%%xmm6		\n\t		addpd	    (%%r15),%%xmm12	\n\t"\
		"addpd	%%xmm5,%%xmm7		\n\t		addpd	0x10(%%r15),%%xmm13	\n\t"\
		"movaps	%%xmm4,    (%%rsi)	\n\t		subpd	%%xmm14,%%xmm12		\n\t"\
		"movaps	%%xmm5,0x10(%%rsi)	\n\t		subpd	%%xmm15,%%xmm13		\n\t"\
		"movaps	%%xmm0,%%xmm4		\n\t		addpd	%%xmm14,%%xmm14		\n\t"\
		"movaps	%%xmm1,%%xmm5		\n\t		addpd	%%xmm15,%%xmm15		\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t		addpd	%%xmm12,%%xmm14		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t		addpd	%%xmm13,%%xmm15		\n\t"\
		"mulpd	0x20(%%rax),%%xmm0	\n\t		movaps	%%xmm12,    (%%r10)	\n\t"\
		"mulpd	0x20(%%rax),%%xmm1	\n\t		movaps	%%xmm13,0x10(%%r10)	\n\t"\
		"mulpd	0x30(%%rax),%%xmm2	\n\t		movaps	%%xmm8 ,%%xmm12		\n\t"\
		"mulpd	0x30(%%rax),%%xmm3	\n\t		movaps	%%xmm9 ,%%xmm13		\n\t"\
		"mulpd	0x40(%%rax),%%xmm4	\n\t		subpd	%%xmm10,%%xmm8 		\n\t"\
		"mulpd	0x40(%%rax),%%xmm5	\n\t		subpd	%%xmm11,%%xmm9 		\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t		mulpd	0x20(%%rax),%%xmm8 	\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t		mulpd	0x20(%%rax),%%xmm9 	\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t		mulpd	0x30(%%rax),%%xmm10	\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t		mulpd	0x30(%%rax),%%xmm11	\n\t"\
		"movaps	    (%%rsi),%%xmm4	\n\t		mulpd	0x40(%%rax),%%xmm12	\n\t"\
		"movaps	0x10(%%rsi),%%xmm5	\n\t		mulpd	0x40(%%rax),%%xmm13	\n\t"\
		"/*------ End of shared-trig-data accessed via rax register ------*/\n\t"\
		"movq	%[__o1],%%rax		\n\t		addpd	%%xmm8 ,%%xmm10		\n\t"\
		"movq	%[__o4],%%rdx		\n\t		addpd	%%xmm9 ,%%xmm11		\n\t"\
		"subpd	%%xmm3,%%xmm6		\n\t		subpd	%%xmm12,%%xmm8 		\n\t"\
		"subpd	%%xmm2,%%xmm7		\n\t		subpd	%%xmm13,%%xmm9 		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t		movaps	    (%%r10),%%xmm12	\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t		movaps	0x10(%%r10),%%xmm13	\n\t"\
		"movaps	%%xmm6,    (%%rax)	\n\t		movq	%[__o6],%%r11		\n\t"\
		"movaps	%%xmm7,0x10(%%rdx)	\n\t		movq	%[__o9],%%r14		\n\t"\
		"addpd	%%xmm6,%%xmm3		\n\t		subpd	%%xmm11,%%xmm14		\n\t"\
		"addpd	%%xmm7,%%xmm2		\n\t		subpd	%%xmm10,%%xmm15		\n\t"\
		"movaps	%%xmm3,    (%%rdx)	\n\t		addpd	%%xmm11,%%xmm11		\n\t"\
		"movaps	%%xmm2,0x10(%%rax)	\n\t		addpd	%%xmm10,%%xmm10		\n\t"\
		"movq	%[__o2],%%rbx		\n\t		movaps	%%xmm14,    (%%r11)	\n\t"\
		"movq	%[__o3],%%rcx		\n\t		movaps	%%xmm15,0x10(%%r14)	\n\t"\
		"subpd	%%xmm1,%%xmm4		\n\t		addpd	%%xmm14,%%xmm11		\n\t"\
		"subpd	%%xmm0,%%xmm5		\n\t		addpd	%%xmm15,%%xmm10		\n\t"\
		"addpd	%%xmm1,%%xmm1		\n\t		movaps	%%xmm11,    (%%r14)	\n\t"\
		"addpd	%%xmm0,%%xmm0		\n\t		movaps	%%xmm10,0x10(%%r11)	\n\t"\
		"movaps	%%xmm4,    (%%rbx)	\n\t		movq	%[__o7],%%r12		\n\t"\
		"movaps	%%xmm5,0x10(%%rcx)	\n\t		movq	%[__o8],%%r13		\n\t"\
		"addpd	%%xmm4,%%xmm1		\n\t		subpd	%%xmm9 ,%%xmm12		\n\t"\
		"addpd	%%xmm5,%%xmm0		\n\t		subpd	%%xmm8 ,%%xmm13		\n\t"\
		"movaps	%%xmm1,    (%%rcx)	\n\t		addpd	%%xmm9 ,%%xmm9 		\n\t"\
		"movaps	%%xmm0,0x10(%%rbx)	\n\t		addpd	%%xmm8 ,%%xmm8 		\n\t"\
		"										movaps	%%xmm12,    (%%r12)	\n\t"\
		"										movaps	%%xmm13,0x10(%%r13)	\n\t"\
		"										addpd	%%xmm12,%%xmm9 		\n\t"\
		"										addpd	%%xmm13,%%xmm8 		\n\t"\
		"										movaps	%%xmm9 ,    (%%r13)	\n\t"\
		"										movaps	%%xmm8 ,0x10(%%r12)	\n\t"\
		"							\n\t"\
		:					/* outputs: none */\
		: [__cc1] "m" (Xcc1)	/* All inputs from memory addresses here */\
		 ,[__i0] "m" (Xi0)\
		 ,[__i1] "m" (Xi1)\
		 ,[__i2] "m" (Xi2)\
		 ,[__i3] "m" (Xi3)\
		 ,[__i4] "m" (Xi4)\
		 ,[__o0] "m" (Xo0)\
		 ,[__o1] "m" (Xo1)\
		 ,[__o2] "m" (Xo2)\
		 ,[__o3] "m" (Xo3)\
		 ,[__o4] "m" (Xo4)\
		 ,[__i5] "m" (Yi0)\
		 ,[__i6] "m" (Yi1)\
		 ,[__i7] "m" (Yi2)\
		 ,[__i8] "m" (Yi3)\
		 ,[__i9] "m" (Yi4)\
		 ,[__o5] "m" (Yo0)\
		 ,[__o6] "m" (Yo1)\
		 ,[__o7] "m" (Yo2)\
		 ,[__o8] "m" (Yo3)\
		 ,[__o9] "m" (Yo4)\
		: "cc","memory","rax","rbx","rcx","rdx","rsi","rdi","r10","r11","r12","r13","r14","r15","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"		/* Clobbered registers */\
	);\
	}

	/* Twiddleless version of SSE2_RADIX8_DIF_TWIDDLE. Inputs enter in memory locations __r0 + [__i1,__i2,__i3,__i4,__i5,__i6,__i7],;
	where r0 is a memory address and the i's are LITERAL [BYTE] OFFSETS. Outputs go into memory locations __o0,__o1,__o2,__o3,__o4,__o5,__o6,__o7, assumed disjoint with inputs:\
	*/
	#define SSE2_RADIX8_DIF_0TWIDDLE(Xr0, Xi1,Xi2,Xi3,Xi4,Xi5,Xi6,Xi7, Xo0,Xo1,Xo2,Xo3,Xo4,Xo5,Xo6,Xo7, Xisrt2)\
	{\
	__asm__ volatile (\
		"/* 1st of 2 radix-4 subtransforms, data in xmm0-7: */\n\t	/* 2nd of 2 radix-4 subtransforms, data in xmm8-15: */\n\t"\
		"movq	%[__r0],%%rax	/* i0 = r00 */	\n\t			movq	%[__i1],%%r10		/* i1 */	\n\t"\
		"movq	%[__i2],%%rbx	/* i2 */		\n\t			movq	%[__i3],%%r11		/* i3 */	\n\t"\
		"movq	%[__i4],%%rcx	/* i4 */		\n\t			movq	%[__i5],%%r12		/* i5 */	\n\t"\
		"movq	%[__i6],%%rdx	/* i6 */		\n\t			movq	%[__i7],%%r13		/* i7 */	\n\t"\
		"addq	%%rax,%%rbx						\n\t			addq	%%rax,%%r10						\n\t"\
		"addq	%%rax,%%rcx						\n\t			addq	%%rax,%%r11						\n\t"\
		"addq	%%rax,%%rdx						\n\t			addq	%%rax,%%r12						\n\t"\
		"movq	%[__isrt2],%%rsi				\n\t			addq	%%rax,%%r13						\n\t"\
		"										\n\t			/* p1,5 combo: x+y into xmm8 /1, x-y in xmm10/3: */	\n\t"\
		"/* p0,4 combo: x+y into xmm0/1, x-y in xmm2/3: */\n\t	movaps	    (%%r12),%%xmm8 				\n\t"\
		"										\n\t			movaps	0x10(%%r12),%%xmm9 				\n\t"\
		"movaps	    (%%rcx),%%xmm0				\n\t			movaps	    (%%r10),%%xmm10				\n\t"\
		"movaps	0x10(%%rcx),%%xmm1				\n\t			movaps	0x10(%%r10),%%xmm11				\n\t"\
		"movaps	    (%%rax),%%xmm2				\n\t			subpd	%%xmm8 ,%%xmm10					\n\t"\
		"movaps	0x10(%%rax),%%xmm3				\n\t			subpd	%%xmm9 ,%%xmm11					\n\t"\
		"subpd	%%xmm0,%%xmm2					\n\t			addpd	%%xmm8 ,%%xmm8 					\n\t"\
		"subpd	%%xmm1,%%xmm3					\n\t			addpd	%%xmm9 ,%%xmm9 					\n\t"\
		"addpd	%%xmm0,%%xmm0					\n\t			addpd	%%xmm10,%%xmm8 					\n\t"\
		"addpd	%%xmm1,%%xmm1					\n\t			addpd	%%xmm11,%%xmm9 					\n\t"\
		"addpd	%%xmm2,%%xmm0					\n\t			/* p3,7 combo: x+y into xmm14/7, x-y in xmm12/5: */	\n\t"\
		"addpd	%%xmm3,%%xmm1					\n\t			movaps	    (%%r11),%%xmm12				\n\t"\
		"										\n\t			movaps	0x10(%%r11),%%xmm13				\n\t"\
		"/* p2,6 combo: x+y into xmm4/5, x-y in xmm6/7: */\n\t	movaps	    (%%r13),%%xmm14				\n\t"\
		"										\n\t			movaps	0x10(%%r13),%%xmm15				\n\t"\
		"movaps	    (%%rdx),%%xmm4				\n\t			subpd	%%xmm14,%%xmm12					\n\t"\
		"movaps	0x10(%%rdx),%%xmm5				\n\t			subpd	%%xmm15,%%xmm13					\n\t"\
		"movaps	    (%%rbx),%%xmm6				\n\t			addpd	%%xmm14,%%xmm14					\n\t"\
		"movaps	0x10(%%rbx),%%xmm7				\n\t			addpd	%%xmm15,%%xmm15					\n\t"\
		"subpd	%%xmm4,%%xmm6					\n\t			addpd	%%xmm12,%%xmm14					\n\t"\
		"subpd	%%xmm5,%%xmm7					\n\t			addpd	%%xmm13,%%xmm15					\n\t"\
		"addpd	%%xmm4,%%xmm4					\n\t			/* Finish radix-4 butterfly, tmp-store 1st of 4 outputs to free up 2 registers: */\n\t"\
		"addpd	%%xmm5,%%xmm5					\n\t			subpd	%%xmm14,%%xmm8 					\n\t"\
		"addpd	%%xmm6,%%xmm4					\n\t			subpd	%%xmm15,%%xmm9 					\n\t"\
		"addpd	%%xmm7,%%xmm5					\n\t			subpd	%%xmm13,%%xmm10					\n\t"\
		"										\n\t			subpd	%%xmm12,%%xmm11					\n\t"\
		"subpd	%%xmm4,%%xmm0					\n\t			addpd	%%xmm14,%%xmm14					\n\t"\
		"subpd	%%xmm7,%%xmm2					\n\t			addpd	%%xmm13,%%xmm13					\n\t"\
		"subpd	%%xmm5,%%xmm1					\n\t			addpd	%%xmm15,%%xmm15					\n\t"\
		"subpd	%%xmm6,%%xmm3					\n\t			addpd	%%xmm12,%%xmm12					\n\t"\
		"														addpd	%%xmm8 ,%%xmm14					\n\t"\
		"														addpd	%%xmm10,%%xmm13					\n\t"\
		"														addpd	%%xmm9 ,%%xmm15					\n\t"\
		"														addpd	%%xmm11,%%xmm12					\n\t"\
		"addpd	%%xmm4,%%xmm4					\n\t			movaps	%%xmm14,    (%%r10)				\n\t"\
		"addpd	%%xmm7,%%xmm7					\n\t			movaps	%%xmm15,0x10(%%r10)				\n\t"\
		"addpd	%%xmm5,%%xmm5					\n\t			movaps	%%xmm10,%%xmm14					\n\t"\
		"addpd	%%xmm6,%%xmm6					\n\t			movaps	%%xmm13,%%xmm15					\n\t"\
		"addpd	%%xmm0,%%xmm4					\n\t			subpd	%%xmm12,%%xmm10					\n\t"\
		"addpd	%%xmm2,%%xmm7					\n\t			subpd	%%xmm11,%%xmm13					\n\t"\
		"addpd	%%xmm1,%%xmm5					\n\t			addpd	%%xmm14,%%xmm12					\n\t"\
		"addpd	%%xmm3,%%xmm6					\n\t			addpd	%%xmm15,%%xmm11					\n\t"\
		"														movaps	(%%rsi),%%xmm14	/* isrt2 */		\n\t"\
		"														mulpd	%%xmm14,%%xmm10					\n\t"\
		"														mulpd	%%xmm14,%%xmm13					\n\t"\
		"														mulpd	%%xmm14,%%xmm12					\n\t"\
		"										\n\t			mulpd	%%xmm14,%%xmm11					\n\t"\
		"movaps	    (%%r10),%%xmm14	/* restore spilled */\n\t"\
		"movaps	0x10(%%r10),%%xmm15	/* restore spilled */\n\t"\
		"										\n\t"\
		"/* Inline of SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_A(r0): Combine radix-4 subtransforms and write outputs: */\n\t"\
		"/***** t0,1,2,3,4,5,6,7 in xmm[ 4, 5| 2, 6| 0, 1| 7, 3] *****/\n\t"\
		"/***** t8,9,a,b,c,d,e,f in xmm[14,15|10,12| 8, 9|13,11] */\n\t"\
		"movq	%[__o4],%%rax					\n\t			subpd   %%xmm10,%%xmm2					\n\t"\
		"movq	%[__o5],%%rbx					\n\t			subpd   %%xmm12,%%xmm6					\n\t"\
		"movq	%[__o6],%%rcx					\n\t			addpd   %%xmm10,%%xmm10					\n\t"\
		"movq	%[__o7],%%rdx					\n\t			addpd   %%xmm12,%%xmm12					\n\t"\
		"										\n\t			addpd   %%xmm2,%%xmm10					\n\t"\
		"subpd   %%xmm11,%%xmm7					\n\t			addpd   %%xmm6,%%xmm12					\n\t"\
		"subpd   %%xmm13,%%xmm3					\n\t													\n\t"\
		"addpd   %%xmm11,%%xmm11				\n\t			movaps	%%xmm2 ,    (%%rbx)	/* o5r */	\n\t"\
		"addpd   %%xmm13,%%xmm13				\n\t			movaps	%%xmm6 ,0x10(%%rbx)	/* o5i */	\n\t"\
		"addpd   %%xmm7,%%xmm11					\n\t			movaps	%%xmm10,    (%%rax)	/* o4r */	\n\t"\
		"addpd   %%xmm3,%%xmm13					\n\t			movaps	%%xmm12,0x10(%%rax)	/* o4i */	\n\t"\
		"										\n\t"\
		"movaps	%%xmm7 ,    (%%rcx)	/* o6r */	\n\t"\
		"movaps	%%xmm3 ,0x10(%%rdx)	/* o7i */	\n\t"\
		"movaps	%%xmm11,    (%%rdx)	/* o7r */	\n\t"\
		"movaps	%%xmm13,0x10(%%rcx)	/* o6i */	\n\t"\
		"										\n\t"\
		"movq	%[__o0],%%rax					\n\t"\
		"movq	%[__o1],%%rbx					\n\t"\
		"movq	%[__o2],%%rcx					\n\t"\
		"movq	%[__o3],%%rdx					\n\t"\
		"										\n\t"\
		"subpd	%%xmm14,%%xmm4 					\n\t"\
		"subpd	%%xmm15,%%xmm5 					\n\t"\
		"subpd	%%xmm9 ,%%xmm0 					\n\t"\
		"subpd	%%xmm8 ,%%xmm1 					\n\t"\
		"addpd	%%xmm14,%%xmm14					\n\t			movaps	%%xmm4 ,    (%%rbx)	/* o1r */	\n\t"\
		"addpd	%%xmm15,%%xmm15					\n\t			movaps	%%xmm5 ,0x10(%%rbx)	/* o1i */	\n\t"\
		"addpd	%%xmm9 ,%%xmm9 					\n\t			movaps	%%xmm0 ,    (%%rcx)	/* o2r */	\n\t"\
		"addpd	%%xmm8 ,%%xmm8 					\n\t			movaps	%%xmm1 ,0x10(%%rdx)	/* o3i */	\n\t"\
		"addpd	%%xmm4 ,%%xmm14					\n\t"\
		"addpd	%%xmm5 ,%%xmm15					\n\t"\
		"addpd	%%xmm0 ,%%xmm9 					\n\t"\
		"addpd	%%xmm1 ,%%xmm8 					\n\t"\
		"										\n\t"\
		"movaps	%%xmm14,    (%%rax)	/* o0r */	\n\t"\
		"movaps	%%xmm15,0x10(%%rax)	/* o0r */	\n\t"\
		"movaps	%%xmm9 ,    (%%rdx)	/* o3r */	\n\t"\
		"movaps	%%xmm8 ,0x10(%%rcx)	/* o2i */	\n\t"\
		"										\n\t"\
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
		: "cc","memory","rax","rbx","rcx","rdx","rsi","r10","r11","r12","r13","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"		/* Clobbered registers */\
	);\
	}

	/* Twiddleless version of SSE2_RADIX8_DIT_TWIDDLE. Inputs enter in memory locations __i0,__i1,__i2,__i3,__i4,__i5,__i6,__i7.
	Outputs go into 16 contiguous 32-byte memory locations starting at __out and assumed disjoint with inputs.
	This macro built on the same code template as SSE2_RADIX8_DIF_TWIDDLE0, but with the I/O-location indices mutually bit reversed:
	01234567 <--> 04261537, which can be effected via the pairwise swaps 1 <--> 4 and 3 <--> 6.
	*/
	#define	SSE2_RADIX8_DIT_0TWIDDLE(Xi0,Xi1,Xi2,Xi3,Xi4,Xi5,Xi6,Xi7, Xout, Xisrt2)\
	{\
	__asm__ volatile (\
		"/* 1st of 2 radix-4 subtransforms, data in xmm0-7: */\n\t	/* 2nd of 2 radix-4 subtransforms, data in xmm8-15: */\n\t"\
		"movq	%[__i0],%%rax					\n\t			movq	%[__i4],%%r10					\n\t"\
		"movq	%[__i1],%%rbx					\n\t			movq	%[__i5],%%r11					\n\t"\
		"movq	%[__i2],%%rcx					\n\t			movq	%[__i6],%%r12					\n\t"\
		"movq	%[__i3],%%rdx					\n\t			movq	%[__i7],%%r13					\n\t"\
		"										\n\t			/* p1,5 combo: x+y into xmm8 /1, x-y in xmm10/3: */	\n\t"\
		"/* p0,4 combo: x+y into xmm0/1, x-y in xmm2/3: */\n\t	movaps	    (%%r11),%%xmm8 				\n\t"\
		"										\n\t			movaps	0x10(%%r11),%%xmm9 				\n\t"\
		"movaps	    (%%rbx),%%xmm0				\n\t			movaps	    (%%r10),%%xmm10				\n\t"\
		"movaps	0x10(%%rbx),%%xmm1				\n\t			movaps	0x10(%%r10),%%xmm11				\n\t"\
		"movaps	    (%%rax),%%xmm2				\n\t			subpd	%%xmm8 ,%%xmm10					\n\t"\
		"movaps	0x10(%%rax),%%xmm3				\n\t			subpd	%%xmm9 ,%%xmm11					\n\t"\
		"subpd	%%xmm0,%%xmm2					\n\t			addpd	%%xmm8 ,%%xmm8 					\n\t"\
		"subpd	%%xmm1,%%xmm3					\n\t			addpd	%%xmm9 ,%%xmm9 					\n\t"\
		"addpd	%%xmm0,%%xmm0					\n\t			addpd	%%xmm10,%%xmm8 					\n\t"\
		"addpd	%%xmm1,%%xmm1					\n\t			addpd	%%xmm11,%%xmm9 					\n\t"\
		"addpd	%%xmm2,%%xmm0					\n\t			/* p3,7 combo: x+y into xmm14/7, x-y in xmm12/5: */	\n\t"\
		"addpd	%%xmm3,%%xmm1					\n\t			movaps	    (%%r12),%%xmm12				\n\t"\
		"										\n\t			movaps	0x10(%%r12),%%xmm13				\n\t"\
		"/* p2,6 combo: x+y into xmm4/5, x-y in xmm6/7: */\n\t	movaps	    (%%r13),%%xmm14				\n\t"\
		"										\n\t			movaps	0x10(%%r13),%%xmm15				\n\t"\
		"movaps	    (%%rdx),%%xmm4				\n\t			subpd	%%xmm14,%%xmm12					\n\t"\
		"movaps	0x10(%%rdx),%%xmm5				\n\t			subpd	%%xmm15,%%xmm13					\n\t"\
		"movaps	    (%%rcx),%%xmm6				\n\t			addpd	%%xmm14,%%xmm14					\n\t"\
		"movaps	0x10(%%rcx),%%xmm7				\n\t			addpd	%%xmm15,%%xmm15					\n\t"\
		"subpd	%%xmm4,%%xmm6					\n\t			addpd	%%xmm12,%%xmm14					\n\t"\
		"subpd	%%xmm5,%%xmm7					\n\t			addpd	%%xmm13,%%xmm15					\n\t"\
		"addpd	%%xmm4,%%xmm4					\n\t			/* Finish radix-4 butterfly, tmp-store 1st of 4 outputs to free up 2 registers: */\n\t"\
		"addpd	%%xmm5,%%xmm5					\n\t			subpd	%%xmm14,%%xmm8 					\n\t"\
		"addpd	%%xmm6,%%xmm4					\n\t			subpd	%%xmm15,%%xmm9 					\n\t"\
		"addpd	%%xmm7,%%xmm5					\n\t			subpd	%%xmm13,%%xmm10					\n\t"\
		"										\n\t			subpd	%%xmm12,%%xmm11					\n\t"\
		"subpd	%%xmm4,%%xmm0					\n\t			addpd	%%xmm14,%%xmm14					\n\t"\
		"subpd	%%xmm7,%%xmm2					\n\t			addpd	%%xmm13,%%xmm13					\n\t"\
		"subpd	%%xmm5,%%xmm1					\n\t			addpd	%%xmm15,%%xmm15					\n\t"\
		"subpd	%%xmm6,%%xmm3					\n\t			addpd	%%xmm12,%%xmm12					\n\t"\
		"														addpd	%%xmm8 ,%%xmm14					\n\t"\
		"														addpd	%%xmm10,%%xmm13					\n\t"\
		"														addpd	%%xmm9 ,%%xmm15					\n\t"\
		"														addpd	%%xmm11,%%xmm12					\n\t"\
		"														movq	%[__isrt2],%%rsi	/* isrt2 */	\n\t"\
		"addpd	%%xmm4,%%xmm4					\n\t			movaps	%%xmm14,    (%%rax)	/* spill */	\n\t"\
		"addpd	%%xmm7,%%xmm7					\n\t			movaps	%%xmm15,0x10(%%rax)	/* spill */	\n\t"\
		"addpd	%%xmm5,%%xmm5					\n\t			movaps	%%xmm10,%%xmm14					\n\t"\
		"addpd	%%xmm6,%%xmm6					\n\t			movaps	%%xmm13,%%xmm15					\n\t"\
		"addpd	%%xmm0,%%xmm4					\n\t			subpd	%%xmm12,%%xmm10					\n\t"\
		"addpd	%%xmm2,%%xmm7					\n\t			subpd	%%xmm11,%%xmm13					\n\t"\
		"addpd	%%xmm1,%%xmm5					\n\t			addpd	%%xmm14,%%xmm12					\n\t"\
		"addpd	%%xmm3,%%xmm6					\n\t			addpd	%%xmm15,%%xmm11					\n\t"\
		"														movaps	(%%rsi),%%xmm14		/* isrt2 */	\n\t"\
		"														mulpd	%%xmm14,%%xmm10					\n\t"\
		"														mulpd	%%xmm14,%%xmm13					\n\t"\
		"														mulpd	%%xmm14,%%xmm12					\n\t"\
		"														mulpd	%%xmm14,%%xmm11					\n\t"\
		"/* Combine radix-4 subtransforms and write outputs: */\n\t"\
		"\n\t"\
		"movaps	    (%%rax),%%xmm14	/* restore spilled */\n\t	subpd   %%xmm10,%%xmm2					\n\t"\
		"movaps	0x10(%%rax),%%xmm15	/* restore spilled */\n\t	subpd   %%xmm12,%%xmm6					\n\t"\
		"														addpd   %%xmm10,%%xmm10					\n\t"\
		"movq	%[__out],%%rax					\n\t			addpd   %%xmm12,%%xmm12					\n\t"\
		"										\n\t			addpd   %%xmm2,%%xmm10					\n\t"\
		"subpd   %%xmm11,%%xmm7					\n\t			addpd   %%xmm6,%%xmm12					\n\t"\
		"subpd   %%xmm13,%%xmm3					\n\t													\n\t"\
		"addpd   %%xmm11,%%xmm11				\n\t			movaps	%%xmm2 ,0xa0(%%rax)	/* o5r */	\n\t"\
		"addpd   %%xmm13,%%xmm13				\n\t			movaps	%%xmm6 ,0xb0(%%rax)	/* o5i */	\n\t"\
		"addpd   %%xmm7,%%xmm11					\n\t			movaps	%%xmm10,0x20(%%rax)	/* o1r */	\n\t"\
		"addpd   %%xmm3,%%xmm13					\n\t			movaps	%%xmm12,0x30(%%rax)	/* o1i */	\n\t"\
		"										\n\t"\
		"movaps	%%xmm7 ,0x60(%%rax)	/* o3r */	\n\t"\
		"movaps	%%xmm3 ,0xf0(%%rax)	/* o7i */	\n\t"\
		"movaps	%%xmm11,0xe0(%%rax)	/* o7r */	\n\t"\
		"movaps	%%xmm13,0x70(%%rax)	/* o3i */	\n\t"\
		"										\n\t"\
		"subpd	%%xmm14,%%xmm4 					\n\t"\
		"subpd	%%xmm15,%%xmm5 					\n\t"\
		"subpd	%%xmm9 ,%%xmm0 					\n\t"\
		"subpd	%%xmm8 ,%%xmm1 					\n\t"\
		"addpd	%%xmm14,%%xmm14					\n\t			movaps	%%xmm4 ,0x80(%%rax)	/* o4r */	\n\t"\
		"addpd	%%xmm15,%%xmm15					\n\t			movaps	%%xmm5 ,0x90(%%rax)	/* o4i */	\n\t"\
		"addpd	%%xmm9 ,%%xmm9 					\n\t			movaps	%%xmm0 ,0x40(%%rax)	/* o2r */	\n\t"\
		"addpd	%%xmm8 ,%%xmm8 					\n\t			movaps	%%xmm1 ,0xd0(%%rax)	/* o6i */	\n\t"\
		"addpd	%%xmm4 ,%%xmm14					\n\t"\
		"addpd	%%xmm5 ,%%xmm15					\n\t"\
		"addpd	%%xmm0 ,%%xmm9 					\n\t"\
		"addpd	%%xmm1 ,%%xmm8 					\n\t"\
		"										\n\t"\
		"movaps	%%xmm14,    (%%rax)	/* o0r */	\n\t"\
		"movaps	%%xmm15,0x10(%%rax)	/* o0r */	\n\t"\
		"movaps	%%xmm9 ,0xc0(%%rax)	/* o6r */	\n\t"\
		"movaps	%%xmm8 ,0x50(%%rax)	/* o2i */	\n\t"\
		"										\n\t"\
		:					/* outputs: none */\
		: [__i0] "m" (Xi0)	/* All iputs from memory addresses here */\
		 ,[__i1] "m" (Xi1)\
		 ,[__i2] "m" (Xi2)\
		 ,[__i3] "m" (Xi3)\
		 ,[__i4] "m" (Xi4)\
		 ,[__i5] "m" (Xi5)\
		 ,[__i6] "m" (Xi6)\
		 ,[__i7] "m" (Xi7)\
		 ,[__out] "m" (Xout)\
		 ,[__isrt2] "m" (Xisrt2)\
		: "cc","memory","rax","rbx","rcx","rdx","rsi","r10","r11","r12","r13","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"		/* Clobbered registers */\
	);\
	}

#endif	/* sse2_macro_gcc_h_included */

