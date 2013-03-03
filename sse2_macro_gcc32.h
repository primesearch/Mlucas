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
		"movl	%[__cA]		,%%eax\n\t"\
		"movl	%[__cB]		,%%ebx\n\t"\
		"movl	%[__cAmB]	,%%ecx\n\t"\
		"movl	%[__cApB]	,%%edx\n\t"\
		"\n\t"\
		"movaps	    (%%eax),%%xmm0\n\t"\
		"movaps	0x10(%%eax),%%xmm2\n\t"\
		"movaps	    (%%ebx),%%xmm4\n\t"\
		"movaps	0x10(%%ebx),%%xmm5\n\t"\
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
		"movaps	%%xmm0,    (%%ecx)\n\t"\
		"movaps	%%xmm1,0x10(%%ecx)\n\t"\
		"movaps	%%xmm4,    (%%edx)\n\t"\
		"movaps	%%xmm5,0x10(%%edx)\n\t"\
		:					/* outputs: none */\
		: [__cA]  "m" (XcA)	/* All inputs from memory addresses here */\
		 ,[__cB]  "m" (XcB)\
		 ,[__cAmB] "m" (XcAmB)\
		 ,[__cApB] "m" (XcApB)\
		: "eax","ebx","ecx","edx"		/* Clobbered registers */\
	);\
	}

	#define PAIR_SQUARE_4_SSE2(XtAr, XtBr, XtCr, XtDr, Xc, Xs, Xforth)\
	{\
	__asm__ volatile (\
		"/*   calculate cross-product terms...\n\t"\
		"	__rt=__tAr* ~tDr+__tAi* ~tDi; __rt=__rt+__rt;\n\t"\
		"	__it=__tAi* ~tDr-__tAr* ~tDi; __it=__it+__it;\n\t"\
		"*/\n\t"\
		"movl	%[__tDr]	,%%edx\n\t"\
		"movl	%[__tAr]	,%%eax\n\t"\
		"\n\t"\
		"movaps	    (%%edx)	,%%xmm6		/* tDr */\n\t"\
		"movaps	0x10(%%edx)	,%%xmm7		/* tDi */\n\t"\
		"movaps	    (%%eax)	,%%xmm0		/* tAr */\n\t"\
		"movaps	0x10(%%eax)	,%%xmm3		/* tAi */\n\t"\
		"shufpd	$1	,%%xmm6	,%%xmm6	/*~tDr */\n\t"\
		"shufpd	$1	,%%xmm7	,%%xmm7	/*~tDi */\n\t"\
		"movaps	    (%%eax)	,%%xmm2		/* cpy tAr */\n\t"\
		"movaps	0x10(%%eax)	,%%xmm1		/* cpy tAi */\n\t"\
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
		"movl	%[__tCr]	,%%ecx\n\t"\
		"movl	%[__tBr]	,%%ebx\n\t"\
		"\n\t"\
		"movaps	    (%%ecx)	,%%xmm6		/* tCr */\n\t"\
		"movaps	0x10(%%ecx)	,%%xmm7		/* tCi */\n\t"\
		"movaps	    (%%ebx)	,%%xmm2		/* tBr */\n\t"\
		"movaps	0x10(%%ebx)	,%%xmm5		/* tBi */\n\t"\
		"shufpd	$1	,%%xmm6	,%%xmm6	/*~tCr */\n\t"\
		"shufpd	$1	,%%xmm7	,%%xmm7	/*~tCi */\n\t"\
		"movaps	    (%%ebx)	,%%xmm4		/* cpy tBr */\n\t"\
		"movaps	0x10(%%ebx)	,%%xmm3		/* cpy tBi */\n\t"\
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
		"movaps	    (%%eax)	,%%xmm4		/* __tAr */\n\t"\
		"movaps	0x10(%%eax)	,%%xmm5		/* __tAi */\n\t"\
		"subpd	%%xmm5		,%%xmm4		/* (__tAr-__tAi) */\n\t"\
		"addpd	%%xmm5		,%%xmm5		/*      2*__tAi  */\n\t"\
		"addpd	%%xmm4		,%%xmm5		/* (__tAr+__tAi) */\n\t"\
		"mulpd	%%xmm5		,%%xmm4		/*>__tAr */\n\t"\
		"\n\t"\
		"movaps	    (%%eax)	,%%xmm5		/* __tAr */\n\t"\
		"mulpd	0x10(%%eax)	,%%xmm5		/* __tAr*__tAi */\n\t"\
		"addpd	%%xmm5		,%%xmm5		/*>__tAi */\n\t"\
		"movaps	%%xmm4	,    (%%eax)	/* tmp store >__tAr */\n\t"\
		"movaps	%%xmm5	,0x10(%%eax)	/* tmp store >__tAi */\n\t"\
		"\n\t"\
		"subpd	%%xmm4		,%%xmm0	/* rt-__tAr */\n\t"\
		"subpd	%%xmm5		,%%xmm1	/* it-__tAi; xmm4-7 free */\n\t"\
		"\n\t"\
		"/*	__tmp=(__tBr+__tBi)*(__tBr-__tBi); __tBi=__tBr*__tBi; __tBi=__tBi+__tBi; __tBr=__tmp;	*/\n\t"\
		"/*** [Can be done in parallel with above segment] ***/\n\t"\
		"\n\t"\
		"movaps	    (%%ebx)	,%%xmm6		/* __tBr */\n\t"\
		"movaps	0x10(%%ebx)	,%%xmm7		/* __tBi */\n\t"\
		"subpd	%%xmm7		,%%xmm6		/* (__tBr-__tBi) */\n\t"\
		"addpd	%%xmm7		,%%xmm7		/*      2*__tBi  */\n\t"\
		"addpd	%%xmm6		,%%xmm7		/* (__tBr+__tBi) */\n\t"\
		"mulpd	%%xmm7		,%%xmm6		/*>__tBr */\n\t"\
		"\n\t"\
		"movaps	    (%%ebx)	,%%xmm7		/* __tBr */\n\t"\
		"mulpd	0x10(%%ebx)	,%%xmm7		/* __tBr*__tBi */\n\t"\
		"addpd	%%xmm7		,%%xmm7		/*>__tBi */\n\t"\
		"movaps	%%xmm6	,    (%%ebx)	/* tmp store >__tBr */\n\t"\
		"movaps	%%xmm7	,0x10(%%ebx)	/* tmp store >__tBi */\n\t"\
		"\n\t"\
		"subpd	%%xmm6		,%%xmm2	/* st-__tBr */\n\t"\
		"subpd	%%xmm7		,%%xmm3	/* jt-__tBi; xmm4-7 free */\n\t"\
		"\n\t"\
		"/*	__tmp=(__tDr+__tDi)*(__tDr-__tDi); __tDi=__tDr*__tDi; __tDi=__tDi+__tDi; __tDr=__tmp;	*/\n\t"\
		"\n\t"\
		"movaps	    (%%edx)	,%%xmm4		/* __tDr */\n\t"\
		"movaps	0x10(%%edx)	,%%xmm5		/* __tDi */\n\t"\
		"subpd	%%xmm5		,%%xmm4		/* (__tDr-__tDi) */\n\t"\
		"addpd	%%xmm5		,%%xmm5		/*      2*__tDi  */\n\t"\
		"addpd	%%xmm4		,%%xmm5		/* (__tDr+__tDi) */\n\t"\
		"mulpd	%%xmm5		,%%xmm4		/*>__tDr */\n\t"\
		"\n\t"\
		"movaps	    (%%edx)	,%%xmm5		/* __tDr */\n\t"\
		"mulpd	0x10(%%edx)	,%%xmm5		/* __tDr*__tDi */\n\t"\
		"addpd	%%xmm5		,%%xmm5		/*>__tDi */\n\t"\
		"movaps	%%xmm4	,    (%%edx)	/* tmp store ~tDr */\n\t"\
		"movaps	%%xmm5	,0x10(%%edx)	/* tmp store ~tDi */\n\t"\
		"shufpd	$1	,%%xmm4	,%%xmm4	/*~tDr */\n\t"\
		"shufpd	$1	,%%xmm5	,%%xmm5	/*~tDi */\n\t"\
		"\n\t"\
		"subpd	%%xmm4		,%%xmm0	/* rt-__tAr- ~tDr */\n\t"\
		"addpd	%%xmm5		,%%xmm1	/* it-__tAi+ ~tDi; xmm4-7 free */\n\t"\
		"\n\t"\
		"/*	__tmp=(__tCr+__tCi)*(__tCr-__tCi); __tCi=__tCr*__tCi; __tCi=__tCi+__tCi; __tCr=__tmp;	*/\n\t"\
		"/*** [Can be done in parallel with above segment] ***/\n\t"\
		"\n\t"\
		"movaps	    (%%ecx)	,%%xmm6		/* __tCr */\n\t"\
		"movaps	0x10(%%ecx)	,%%xmm7		/* __tCi */\n\t"\
		"subpd	%%xmm7		,%%xmm6		/* (__tCr-__tCi) */\n\t"\
		"addpd	%%xmm7		,%%xmm7		/*      2*__tCi  */\n\t"\
		"addpd	%%xmm6		,%%xmm7		/* (__tCr+__tCi) */\n\t"\
		"mulpd	%%xmm7		,%%xmm6		/*>__tCr */\n\t"\
		"\n\t"\
		"movaps	    (%%ecx)	,%%xmm7		/* __tCr */\n\t"\
		"mulpd	0x10(%%ecx)	,%%xmm7		/* __tCr*__tCi */\n\t"\
		"addpd	%%xmm7		,%%xmm7		/*>__tCi */\n\t"\
		"movaps	%%xmm6	,    (%%ecx)	/* tmp store ~tCr */\n\t"\
		"movaps	%%xmm7	,0x10(%%ecx)	/* tmp store ~tCi */\n\t"\
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
		"movl	%[__c]		,%%eax\n\t"\
		"movl	%[__s]		,%%ebx\n\t"\
		"movl	%[__forth]	,%%edx\n\t"\
		"movaps	%%xmm0		,%%xmm4		/* cpy rt */\n\t"\
		"movaps	%%xmm1		,%%xmm5		/* cpy it */\n\t"\
		"mulpd	(%%eax)		,%%xmm0		/* c*rt */\n\t"\
		"mulpd	(%%eax)		,%%xmm1		/* c*it */\n\t"\
		"addpd	%%xmm4		,%%xmm0		/* (c+1.0)*rt */\n\t"\
		"addpd	%%xmm5		,%%xmm1		/* (c+1.0)*it */\n\t"\
		"mulpd	(%%ebx)		,%%xmm4		/* s*rt */\n\t"\
		"mulpd	(%%ebx)		,%%xmm5		/* s*it */\n\t"\
		"subpd	%%xmm5		,%%xmm0		/* (c+1.0)*rt-s*it */\n\t"\
		"addpd	%%xmm4		,%%xmm1		/* (c+1.0)*it+s*rt; xmm4,5 free */\n\t"\
		"mulpd	(%%edx)		,%%xmm0	/* -rt Both of these inherit the sign flip [w.r.to the non-SSE2 PAIR_SQUARE_4 macro] */\n\t"\
		"mulpd	(%%edx)		,%%xmm1	/* -it that resulted from the in-place-friendlier (rt-__tAr- ~tDr) reordering above. */\n\t"\
		"/*\n\t"\
		"	__tmp=((1.0-__s)*__st-__c*__jt)*0.25;\n\t"\
		"	__jt =((1.0-__s)*__jt+__c*__st)*0.25	__st=__tmp;\n\t"\
		"*/\n\t"\
		"/*** [Can be done in parallel wjth above segment] ***/\n\t"\
		"movaps	%%xmm2		,%%xmm6		/* cpy st */\n\t"\
		"movaps	%%xmm3		,%%xmm7		/* cpy jt */\n\t"\
		"mulpd	(%%ebx)		,%%xmm2		/* s*st */\n\t"\
		"mulpd	(%%ebx)		,%%xmm3		/* s*jt */\n\t"\
		"subpd	%%xmm6		,%%xmm2		/* (s-1.0)*st, note sign flip! */\n\t"\
		"subpd	%%xmm7		,%%xmm3		/* (s-1.0)*jt, note sign flip! */\n\t"\
		"mulpd	(%%eax)		,%%xmm6		/* c*st */\n\t"\
		"mulpd	(%%eax)		,%%xmm7		/* c*jt */\n\t"\
		"addpd	%%xmm7		,%%xmm2		/* -[(1.0-s)*st-c*jt] */\n\t"\
		"subpd	%%xmm6		,%%xmm3		/* -[(1.0-s)*jt+c*st]; xmm6,7 free */\n\t"\
		"mulpd	(%%edx)		,%%xmm2	/* +st Sign flip due to (s-1.0) reordering here */\n\t"\
		"mulpd	(%%edx)		,%%xmm3	/* +jt cancels earlier one due to in-place-friendlier (st-__tBr- ~tCr) reordering above. */\n\t"\
		"/*...and now complete and store the results. We flip the signs on st and jt here to undo the above -st,-jt negations. */\n\t"\
		"/*	__tAr = (__tAr+__rt);\n\t"\
		"	__tAi = (__tAi+__it);\n\t"\
		"	__tBr = (__tBr-__st);\n\t"\
		"	__tBi = (__tBi-__jt);\n\t"\
		"*/\n\t"\
		"movl	%[__tAr]	,%%eax\n\t"\
		"movl	%[__tBr]	,%%ebx\n\t"\
		"\n\t"\
		"movaps	    (%%eax)	,%%xmm4		/* __tAr */\n\t"\
		"movaps	0x10(%%eax)	,%%xmm5		/* __tAi */\n\t"\
		"movaps	    (%%ebx)	,%%xmm6		/* __tBr */\n\t"\
		"movaps	0x10(%%ebx)	,%%xmm7		/* __tBi */\n\t"\
		"addpd	%%xmm0		,%%xmm4		/* (__tAr+__rt) */\n\t"\
		"addpd	%%xmm1		,%%xmm5		/* (__tAi+__it) */\n\t"\
		"subpd	%%xmm2		,%%xmm6		/* (__tBr-__st) */\n\t"\
		"subpd	%%xmm3		,%%xmm7		/* (__tBi-__jt) */\n\t"\
		"movaps	%%xmm4	,    (%%eax)	/* store >__tAr */\n\t"\
		"movaps	%%xmm5	,0x10(%%eax)	/* store >__tAi */\n\t"\
		"movaps	%%xmm6	,    (%%ebx)	/* store >__tBr */\n\t"\
		"movaps	%%xmm7	,0x10(%%ebx)	/* store >__tBi */\n\t"\
		"/*...N-j terms are as above, but with the replacements: __tAr<--> ~tDr, __tAi<--> ~tDi, __it|-->-__it. */\n\t"\
		"/*	__tDr = (__tDr+ ~rt);\n\t"\
		"	__tDi = (__tDi- ~it);\n\t"\
		"	__tCr = (__tCr- ~st);\n\t"\
		"	__tCi = (__tCi+ ~jt);\n\t"\
		"*/\n\t"\
		"movl	%[__tCr]	,%%ecx\n\t"\
		"movl	%[__tDr]	,%%edx\n\t"\
		"\n\t"\
		"shufpd	$1	,%%xmm0	,%%xmm0		/* ~rt */\n\t"\
		"shufpd	$1	,%%xmm1	,%%xmm1		/* ~it */\n\t"\
		"shufpd	$1	,%%xmm2	,%%xmm2		/* ~st */\n\t"\
		"shufpd	$1	,%%xmm3	,%%xmm3		/* ~jt */\n\t"\
		"\n\t"\
		"movaps	    (%%edx)	,%%xmm4		/* __tDr */\n\t"\
		"movaps	0x10(%%edx)	,%%xmm5		/* __tDi */\n\t"\
		"movaps	    (%%ecx)	,%%xmm6		/* __tCr */\n\t"\
		"movaps	0x10(%%ecx)	,%%xmm7		/* __tCi */\n\t"\
		"addpd	%%xmm0		,%%xmm4		/* (__tDr+ ~rt) */\n\t"\
		"subpd	%%xmm1		,%%xmm5		/* (__tDi- ~it) */\n\t"\
		"subpd	%%xmm2		,%%xmm6		/* (__tCr- ~st) */\n\t"\
		"addpd	%%xmm3		,%%xmm7		/* (__tCi+ ~jt) */\n\t"\
		"movaps	%%xmm4	,    (%%edx)	/* store >__tDr */\n\t"\
		"movaps	%%xmm5	,0x10(%%edx)	/* store >__tDi */\n\t"\
		"movaps	%%xmm6	,    (%%ecx)	/* store >__tCr */\n\t"\
		"movaps	%%xmm7	,0x10(%%ecx)	/* store >__tCi */\n\t"\
		:					/* outputs: none */\
		: [__tAr] "m" (XtAr)	/* All inputs from memory addresses here */\
		 ,[__tBr] "m" (XtBr)\
		 ,[__tCr] "m" (XtCr)\
		 ,[__tDr] "m" (XtDr)\
		 ,[__c] "m" (Xc)\
		 ,[__s] "m" (Xs)\
		 ,[__forth] "m" (Xforth)\
		: "eax","ebx","ecx","edx"		/* Clobbered registers */\
	);\
	}

	#define SSE2_RADIX_03_DFT(Xi0,Xi1,Xi2, Xcc1, Xo0,Xo1,Xo2)\
	{\
	__asm__ volatile (\
			"movl	%[__i0],%%eax		\n\t"\
			"movl	%[__i1],%%ebx		\n\t"\
			"movl	%[__i2],%%ecx		\n\t"\
			"movl	%[__cc1],%%edx		\n\t"\
			"\n\t"\
			"movaps	    (%%ebx),%%xmm2	\n\t"\
			"movaps	0x10(%%ebx),%%xmm3	\n\t"\
			"movaps	    (%%eax),%%xmm0	\n\t"\
			"movaps	0x10(%%eax),%%xmm1	\n\t"\
			"movaps	    (%%ecx),%%xmm6	\n\t"\
			"movaps	0x10(%%ecx),%%xmm7	\n\t"\
			"movaps	%%xmm2,%%xmm4		\n\t"\
			"movaps	%%xmm3,%%xmm5		\n\t"\
			"\n\t"\
			"movl	%[__o0],%%eax		\n\t"\
			"movl	%[__o1],%%ebx		\n\t"\
			"movl	%[__o2],%%ecx		\n\t"\
			"addpd	%%xmm6,%%xmm2		\n\t"\
			"addpd	%%xmm7,%%xmm3		\n\t"\
			"subpd	%%xmm6,%%xmm4		\n\t"\
			"subpd	%%xmm7,%%xmm5		\n\t"\
			"addpd	%%xmm2,%%xmm0		\n\t"\
			"addpd	%%xmm3,%%xmm1		\n\t"\
			"movaps	    (%%edx),%%xmm6	\n\t"\
			"movaps	0x10(%%edx),%%xmm7	\n\t"\
			"movaps	%%xmm0,     (%%eax)	\n\t"\
			"movaps	%%xmm1,0x010(%%eax)	\n\t"\
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
			"movaps	%%xmm2,     (%%ebx)	\n\t"\
			"movaps	%%xmm3,0x010(%%ebx)	\n\t"\
			"movaps	%%xmm0,     (%%ecx)	\n\t"\
			"movaps	%%xmm1,0x010(%%ecx)	\n\t"\
		:					/* outputs: none */\
		: [__i0] "m" (Xi0)	/* All inputs from memory addresses here */\
		 ,[__i1] "m" (Xi1)\
		 ,[__i2] "m" (Xi2)\
		 ,[__cc1] "m" (Xcc1)\
		 ,[__o0] "m" (Xo0)\
		 ,[__o1] "m" (Xo1)\
		 ,[__o2] "m" (Xo2)\
		: "rax","rbx","rcx","rdx"		/* Clobbered registers */\
	);\
	}

	#define SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(Xadd0, Xadd1, Xadd2, Xadd3, Xtmp, Xstride)\
	{\
	__asm__ volatile (\
		"movl	%[__tmp]   ,%%eax	\n\t"\
		"movl	%[__stride],%%esi	\n\t"\
		"movl	%%eax,%%ebx			\n\t"\
		"addl	%%esi,%%ebx			/* add_in1  */\n\t"\
		"shll	$1,%%esi			/* stride*2 */\n\t"\
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	    (%%ebx),%%xmm2	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"movaps	0x10(%%ebx),%%xmm3	\n\t"\
		"movaps	    (%%eax),%%xmm4	\n\t"\
		"movaps	    (%%ebx),%%xmm6	\n\t"\
		"movaps	0x10(%%eax),%%xmm5	\n\t"\
		"movaps	0x10(%%ebx),%%xmm7	\n\t"\
		"addl	%%esi,%%eax			/* add_in2  */\n\t"\
		"addl	%%esi,%%ebx			/* add_in3  */\n\t"\
		"addpd	    (%%eax),%%xmm0	\n\t"\
		"addpd	    (%%ebx),%%xmm2	\n\t"\
		"addpd	0x10(%%eax),%%xmm1	\n\t"\
		"addpd	0x10(%%ebx),%%xmm3	\n\t"\
		"subpd	    (%%eax),%%xmm4	\n\t"\
		"subpd	    (%%ebx),%%xmm6	\n\t"\
		"subpd	0x10(%%eax),%%xmm5	\n\t"\
		"subpd	0x10(%%ebx),%%xmm7	\n\t"\
		"/* Finish radix-4 butterfly and store results into main-array slots: */\n\t"\
		"movl	%[__add0],%%eax		\n\t"\
		"movl	%[__add1],%%ebx		\n\t"\
		"movl	%[__add2],%%ecx		\n\t"\
		"movl	%[__add3],%%edx		\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm4		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"subpd	%%xmm6,%%xmm5		\n\t"\
		"movaps	%%xmm0,     (%%ebx)	\n\t"\
		"movaps	%%xmm4,     (%%ecx)	\n\t"\
		"movaps	%%xmm1,0x010(%%ebx)	\n\t"\
		"movaps	%%xmm5,0x010(%%edx)	\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm4,%%xmm7		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"addpd	%%xmm5,%%xmm6		\n\t"\
		"movaps	%%xmm2,     (%%eax)	\n\t"\
		"movaps	%%xmm7,     (%%edx)	\n\t"\
		"movaps	%%xmm3,0x010(%%eax)	\n\t"\
		"movaps	%%xmm6,0x010(%%ecx)	\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__add1] "m" (Xadd1)\
		 ,[__add2] "m" (Xadd2)\
		 ,[__add3] "m" (Xadd3)\
		 ,[__tmp] "m" (Xtmp)\
		 ,[__stride] "e" (Xstride)\
		: "eax","ebx","ecx","edx","rsi"		/* Clobbered registers */\
	);\
	}

	#define SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(Xadd0, Xadd1, Xadd2, Xadd3, Xtmp, Xstride)\
	{\
	__asm__ volatile (\
		"movl	%[__add0],%%eax		\n\t"\
		"movl	%[__add1],%%ebx		\n\t"\
		"movl	%[__add2],%%ecx		\n\t"\
		"movl	%[__add3],%%edx		\n\t"\
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	    (%%ecx),%%xmm4	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"movaps	0x10(%%ecx),%%xmm5	\n\t"\
		"movaps	%%xmm0,%%xmm2			\n\t"\
		"movaps	%%xmm4,%%xmm6			\n\t"\
		"movaps	%%xmm1,%%xmm3			\n\t"\
		"movaps	%%xmm5,%%xmm7			\n\t"\
		"movl	%[__tmp]   ,%%eax	\n\t"\
		"movl	%[__stride],%%ecx	\n\t"\
		"addpd	    (%%ebx),%%xmm0	\n\t"\
		"addpd	    (%%edx),%%xmm4	\n\t"\
		"addpd	0x10(%%ebx),%%xmm1	\n\t"\
		"addpd	0x10(%%edx),%%xmm5	\n\t"\
		"subpd	    (%%ebx),%%xmm2	\n\t"\
		"subpd	    (%%edx),%%xmm6	\n\t"\
		"subpd	0x10(%%ebx),%%xmm3	\n\t"\
		"subpd	0x10(%%edx),%%xmm7	\n\t"\
		"movl	%%eax,%%ebx			\n\t"\
		"addl	%%ecx,%%ebx			\n\t"\
		"movl	%%ebx,%%edx			\n\t"\
		"addl	%%ecx,%%ecx			\n\t"\
		"addl	%%ecx,%%edx			\n\t"\
		"addl	%%eax,%%ecx			\n\t"\
		"/* Finish radix-4 butterfly and store results into temp-array slots: */\n\t"\
		"subpd	%%xmm4,%%xmm0			\n\t"\
		"subpd	%%xmm7,%%xmm2			\n\t"\
		"subpd	%%xmm5,%%xmm1			\n\t"\
		"subpd	%%xmm6,%%xmm3			\n\t"\
		"movaps	%%xmm0,     (%%ecx)	\n\t"\
		"movaps	%%xmm2,     (%%edx)	\n\t"\
		"movaps	%%xmm1,0x010(%%ecx)	\n\t"\
		"movaps	%%xmm3,0x010(%%ebx)	\n\t"\
		"addpd	%%xmm4,%%xmm4			\n\t"\
		"addpd	%%xmm7,%%xmm7			\n\t"\
		"addpd	%%xmm5,%%xmm5			\n\t"\
		"addpd	%%xmm6,%%xmm6			\n\t"\
		"addpd	%%xmm0,%%xmm4			\n\t"\
		"addpd	%%xmm2,%%xmm7			\n\t"\
		"addpd	%%xmm1,%%xmm5			\n\t"\
		"addpd	%%xmm3,%%xmm6			\n\t"\
		"movaps	%%xmm4,     (%%eax)	\n\t"\
		"movaps	%%xmm7,     (%%ebx)	\n\t"\
		"movaps	%%xmm5,0x010(%%eax)	\n\t"\
		"movaps	%%xmm6,0x010(%%edx)	\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__add1] "m" (Xadd1)\
		 ,[__add2] "m" (Xadd2)\
		 ,[__add3] "m" (Xadd3)\
		 ,[__tmp] "m" (Xtmp)\
		 ,[__stride] "e" (Xstride)\
		: "eax","ebx","ecx","edx","esi"		/* Clobbered registers */\
	);\
	}

	#define SSE2_RADIX_05_DFT_0TWIDDLE(Xi0,Xi1,Xi2,Xi3,Xi4, Xcc1, Xo0,Xo1,Xo2,Xo3,Xo4)\
	{\
	__asm__ volatile (\
		"movl	%[__i0],%%esi		\n\t"\
		"movl	%[__i1],%%eax		\n\t"\
		"movl	%[__i2],%%ebx		\n\t"\
		"movl	%[__i3],%%ecx		\n\t"\
		"movl	%[__i4],%%edx		\n\t"\
		"movl	%[__o0],%%edi		\n\t"\
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"movaps	    (%%ebx),%%xmm2	\n\t"\
		"movaps	0x10(%%ebx),%%xmm3	\n\t"\
		"movaps	    (%%ecx),%%xmm4	\n\t"\
		"movaps	0x10(%%ecx),%%xmm5	\n\t"\
		"movaps	    (%%edx),%%xmm6	\n\t"\
		"movaps	0x10(%%edx),%%xmm7	\n\t"\
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
		"movl	%[__cc1],%%eax		\n\t"\
		"subpd	%%xmm4,%%xmm6		\n\t"\
		"subpd	%%xmm5,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm5		\n\t"\
		"addpd	    (%%esi),%%xmm4	\n\t"\
		"addpd	0x10(%%esi),%%xmm5	\n\t"\
		"movaps	%%xmm4,    (%%edi)	\n\t"\
		"movaps	%%xmm5,0x10(%%edi)	\n\t"\
		"mulpd	0x10(%%eax),%%xmm6	\n\t"\
		"mulpd	0x10(%%eax),%%xmm7	\n\t"\
		"subpd	     (%%esi),%%xmm4	\n\t"\
		"subpd	0x010(%%esi),%%xmm5	\n\t"\
		"mulpd	    (%%eax),%%xmm4	\n\t"\
		"mulpd	    (%%eax),%%xmm5	\n\t"\
		"addpd	     (%%edi),%%xmm4	\n\t"\
		"addpd	0x010(%%edi),%%xmm5	\n\t"\
		"subpd	%%xmm6,%%xmm4		\n\t"\
		"subpd	%%xmm7,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm6		\n\t"\
		"addpd	%%xmm5,%%xmm7		\n\t"\
		"movaps	%%xmm4,    (%%esi)	\n\t"\
		"movaps	%%xmm5,0x10(%%esi)	\n\t"\
		"movaps	%%xmm0,%%xmm4		\n\t"\
		"movaps	%%xmm1,%%xmm5		\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"mulpd	0x20(%%eax),%%xmm0	\n\t"\
		"mulpd	0x20(%%eax),%%xmm1	\n\t"\
		"mulpd	0x30(%%eax),%%xmm2	\n\t"\
		"mulpd	0x30(%%eax),%%xmm3	\n\t"\
		"mulpd	0x40(%%eax),%%xmm4	\n\t"\
		"mulpd	0x40(%%eax),%%xmm5	\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"movaps	    (%%esi),%%xmm4	\n\t"\
		"movaps	0x10(%%esi),%%xmm5	\n\t"\
		"movl	%[__o1],%%eax		\n\t"\
		"movl	%[__o4],%%edx		\n\t"\
		"subpd	%%xmm3,%%xmm6		\n\t"\
		"subpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"movaps	%%xmm6,    (%%eax)	\n\t"\
		"movaps	%%xmm7,0x10(%%edx)	\n\t"\
		"addpd	%%xmm6,%%xmm3		\n\t"\
		"addpd	%%xmm7,%%xmm2		\n\t"\
		"movaps	%%xmm3,    (%%edx)	\n\t"\
		"movaps	%%xmm2,0x10(%%eax)	\n\t"\
		"movl	%[__o2],%%ebx		\n\t"\
		"movl	%[__o3],%%ecx		\n\t"\
		"subpd	%%xmm1,%%xmm4		\n\t"\
		"subpd	%%xmm0,%%xmm5		\n\t"\
		"addpd	%%xmm1,%%xmm1		\n\t"\
		"addpd	%%xmm0,%%xmm0		\n\t"\
		"movaps	%%xmm4,    (%%ebx)	\n\t"\
		"movaps	%%xmm5,0x10(%%ecx)	\n\t"\
		"addpd	%%xmm4,%%xmm1		\n\t"\
		"addpd	%%xmm5,%%xmm0		\n\t"\
		"movaps	%%xmm1,    (%%ecx)	\n\t"\
		"movaps	%%xmm0,0x10(%%ebx)	\n\t"\
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
		: "eax","ebx","ecx","edx","edi","esi"		/* Clobbered registers */\
	);\
	}

	/* Twiddleless version of SSE2_RADIX8_DIF_TWIDDLE. Inputs enter in memory locations __r0 + [__i1,__i2,__i3,__i4,__i5,__i6,__i7],;
	where r0 is a memory address and the i's are LITERAL [BYTE] OFFSETS. Outputs go into memory locations __o0,__o1,__o2,__o3,__o4,__o5,__o6,__o7, assumed disjoint with inputs:\
	NOTE: The 64-bit version of this macro uses a non-stringified GCC-inline synatx for the literal address offsets __i1-7, by way of reference.
	*/
	#define SSE2_RADIX8_DIF_0TWIDDLE(Xr0, __i1,__i2,__i3,__i4,__i5,__i6,__i7, Xo0,Xo1,Xo2,Xo3,Xo4,Xo5,Xo6,Xo7, Xisrt2)\
	{\
	__asm__ volatile (\
		"movl	%[__r0],%%edi			\n\t"\
		"movl	%[__r0],%%esi			\n\t"\
		"addl	$0x10,%%esi				\n\t"\
		"\n\t"\
		"/* Do the p0,p4 combo: */		\n\t"\
		"\n\t"\
		"movaps	       (%%edi),%%xmm0	\n\t"\
		"movaps	"#__i4"(%%edi),%%xmm4	\n\t"\
		"movaps	       (%%esi),%%xmm1	\n\t"\
		"movaps	"#__i4"(%%esi),%%xmm5	\n\t"\
		"movaps	%%xmm0,%%xmm2			\n\t"\
		"movaps	%%xmm1,%%xmm3			\n\t"\
		"\n\t"\
		"addpd	%%xmm4,%%xmm0			\n\t"\
		"addpd	%%xmm5,%%xmm1			\n\t"\
		"subpd	%%xmm4,%%xmm2			\n\t"\
		"subpd	%%xmm5,%%xmm3			\n\t"\
		"\n\t"\
		"/* Do the p2,6 combo: */		\n\t"\
		"\n\t"\
		"movaps	"#__i2"(%%edi),%%xmm4	\n\t"\
		"movaps	"#__i2"(%%esi),%%xmm5	\n\t"\
		"movaps	%%xmm4,%%xmm6			\n\t"\
		"movaps	%%xmm5,%%xmm7			\n\t"\
		"\n\t"\
		"addpd	"#__i6"(%%edi),%%xmm4	\n\t"\
		"addpd	"#__i6"(%%esi),%%xmm5	\n\t"\
		"subpd	"#__i6"(%%edi),%%xmm6	\n\t"\
		"subpd	"#__i6"(%%esi),%%xmm7	\n\t"\
		"\n\t"\
		"subpd	%%xmm4,%%xmm0			\n\t"\
		"subpd	%%xmm7,%%xmm2			\n\t"\
		"subpd	%%xmm5,%%xmm1			\n\t"\
		"subpd	%%xmm6,%%xmm3			\n\t"\
		"movaps	%%xmm0,"#__i4"(%%edi)	\n\t"\
		"movaps	%%xmm2,"#__i2"(%%edi)	\n\t"\
		"movaps	%%xmm1,"#__i4"(%%esi)	\n\t"\
		"movaps	%%xmm3,"#__i6"(%%esi)	\n\t"\
		"addpd	%%xmm4,%%xmm4			\n\t"\
		"addpd	%%xmm7,%%xmm7			\n\t"\
		"addpd	%%xmm5,%%xmm5			\n\t"\
		"addpd	%%xmm6,%%xmm6			\n\t"\
		"addpd	%%xmm0,%%xmm4			\n\t"\
		"addpd	%%xmm2,%%xmm7			\n\t"\
		"addpd	%%xmm1,%%xmm5			\n\t"\
		"addpd	%%xmm3,%%xmm6			\n\t"\
		"movaps	%%xmm4,       (%%edi)	\n\t"\
		"movaps	%%xmm7,"#__i6"(%%edi)	\n\t"\
		"movaps	%%xmm5,       (%%esi)	\n\t"\
		"movaps	%%xmm6,"#__i2"(%%esi)	\n\t"\
		"\n\t"\
		"/* Inline of SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO_B(r8,c1): Do 2nd of 2 radix-4 subtransforms, but now keep outputs in registers: */\n\t"\
		"\n\t"\
		"/* Do the p1,p5 combo: */		\n\t"\
		"\n\t"\
		"movaps	"#__i1"(%%edi),%%xmm0	\n\t"\
		"movaps	"#__i5"(%%edi),%%xmm4	\n\t"\
		"movaps	"#__i1"(%%esi),%%xmm1	\n\t"\
		"movaps	"#__i5"(%%esi),%%xmm5	\n\t"\
		"movaps	%%xmm0,%%xmm2			\n\t"\
		"movaps	%%xmm1,%%xmm3			\n\t"\
		"\n\t"\
		"addpd	%%xmm4,%%xmm0			\n\t"\
		"addpd	%%xmm5,%%xmm1			\n\t"\
		"subpd	%%xmm4,%%xmm2			\n\t"\
		"subpd	%%xmm5,%%xmm3			\n\t"\
		"\n\t"\
		"/* Do the p3,7 combo: */		\n\t"\
		"\n\t"\
		"movaps	"#__i3"(%%edi),%%xmm4	\n\t"\
		"movaps	"#__i3"(%%esi),%%xmm5	\n\t"\
		"movaps	%%xmm4,%%xmm6			\n\t"\
		"movaps	%%xmm5,%%xmm7			\n\t"\
		"\n\t"\
		"subpd	"#__i7"(%%edi),%%xmm4	\n\t"\
		"subpd	"#__i7"(%%esi),%%xmm5	\n\t"\
		"addpd	"#__i7"(%%edi),%%xmm6	\n\t"\
		"addpd	"#__i7"(%%esi),%%xmm7	\n\t"\
		"\n\t"\
		"/* Finish radix-4 butterfly and store just the 1st of the 8 outputs into output-array slots: */\n\t"\
		"\n\t"\
		"subpd	%%xmm6,%%xmm0			\n\t"\
		"subpd	%%xmm7,%%xmm1			\n\t"\
		"subpd	%%xmm5,%%xmm2			\n\t"\
		"subpd	%%xmm4,%%xmm3			\n\t"\
		"addpd	%%xmm6,%%xmm6			\n\t"\
		"addpd	%%xmm5,%%xmm5			\n\t"\
		"addpd	%%xmm7,%%xmm7			\n\t"\
		"addpd	%%xmm4,%%xmm4			\n\t"\
		"addpd	%%xmm0,%%xmm6			\n\t"\
		"addpd	%%xmm2,%%xmm5			\n\t"\
		"addpd	%%xmm1,%%xmm7			\n\t"\
		"addpd	%%xmm3,%%xmm4			\n\t"\
		"movaps	%%xmm6,"#__i1"(%%edi)	\n\t"\
		"movaps	%%xmm7,"#__i1"(%%esi)	\n\t"\
		"\n\t"\
		"movaps	%%xmm2,%%xmm6			\n\t"\
		"movaps	%%xmm5,%%xmm7			\n\t"\
		"subpd	%%xmm4,%%xmm2			\n\t"\
		"subpd	%%xmm3,%%xmm5			\n\t"\
		"addpd	%%xmm6,%%xmm4			\n\t"\
		"addpd	%%xmm7,%%xmm3			\n\t"\
		"\n\t"\
		"movl	%[__isrt2],%%eax		\n\t"\
		"movaps	(%%eax),%%xmm6	/* isrt2 */\n\t"\
		"mulpd	%%xmm6,%%xmm2			\n\t"\
		"mulpd	%%xmm6,%%xmm5			\n\t"\
		"mulpd	%%xmm6,%%xmm4			\n\t"\
		"mulpd	%%xmm6,%%xmm3			\n\t"\
		"\n\t"\
		"/* Inline the contents of SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_A(r0): Combine the 2 radix-4 subtransforms and write outputs to temp array: */\n\t"\
		"/***** t0,1,2,3,4,5,6,7 in __i0,2,4,6 (eax,ebx,ecx,edx) *****/\n\t"\
		"/***** t8,9,a,b,c,d,e,f in __i1,xmm2,4,0,1,5,3; xmm6,7 FREE */\n\t"\
		"movl	%[__o4],%%eax	\n\t"\
		"movl	%[__o5],%%ebx	\n\t"\
		"movl	%[__o6],%%ecx	\n\t"\
		"movl	%[__o7],%%edx	\n\t"\
		"\n\t"\
		"movaps	"#__i2"(%%edi),%%xmm6		\n\t"\
		"movaps	"#__i2"(%%esi),%%xmm7		\n\t"\
		"subpd   %%xmm2,%%xmm6				\n\t"\
		"subpd   %%xmm4,%%xmm7				\n\t"\
		"addpd   %%xmm2,%%xmm2				\n\t"\
		"addpd   %%xmm4,%%xmm4				\n\t"\
		"addpd   %%xmm6,%%xmm2				\n\t"\
		"addpd   %%xmm7,%%xmm4				\n\t"\
		"\n\t"\
		"movaps	%%xmm6,    (%%ebx)	/* o5r */\n\t"\
		"movaps	%%xmm7,0x10(%%ebx)	/* o5i */\n\t"\
		"movaps	%%xmm2,    (%%eax)	/* o4r */\n\t"\
		"movaps	%%xmm4,0x10(%%eax)	/* o4i */\n\t"\
		"\n\t"\
		"/* Can't do 2nd set side-by-side with 1st, since only 2 registers (xmm6,7) available for temp storage: */\n\t"\
		"movaps	"#__i6"(%%edi),%%xmm6		\n\t"\
		"movaps	"#__i6"(%%esi),%%xmm7		\n\t"\
		"subpd   %%xmm3,%%xmm6				\n\t"\
		"subpd   %%xmm5,%%xmm7				\n\t"\
		"addpd   %%xmm3,%%xmm3				\n\t"\
		"addpd   %%xmm5,%%xmm5				\n\t"\
		"addpd   %%xmm6,%%xmm3				\n\t"\
		"addpd   %%xmm7,%%xmm5				\n\t"\
		"\n\t"\
		"movaps	%%xmm6,    (%%ecx)	/* o6r */\n\t"\
		"movaps	%%xmm7,0x10(%%edx)	/* o7i */\n\t"\
		"movaps	%%xmm3,    (%%edx)	/* o7r */\n\t"\
		"movaps	%%xmm5,0x10(%%ecx)	/* o6i */\n\t"\
		"\n\t"\
		"movl	%[__o0],%%eax	\n\t"\
		"movl	%[__o1],%%ebx	\n\t"\
		"movl	%[__o2],%%ecx	\n\t"\
		"movl	%[__o3],%%edx	\n\t"\
		"\n\t"\
		"movaps	       (%%edi),%%xmm6		\n\t"\
		"movaps	"#__i4"(%%edi),%%xmm4		\n\t"\
		"movaps	       (%%esi),%%xmm7		\n\t"\
		"movaps	"#__i4"(%%esi),%%xmm5		\n\t"\
		"movaps	"#__i1"(%%edi),%%xmm2		\n\t"\
		"movaps	"#__i1"(%%esi),%%xmm3		\n\t"\
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
		"movaps	%%xmm6,    (%%ebx)	/* o1r */\n\t"\
		"movaps	%%xmm4,    (%%ecx)	/* o2r */\n\t"\
		"movaps	%%xmm7,0x10(%%ebx)	/* o1i */\n\t"\
		"movaps	%%xmm5,0x10(%%edx)	/* o3i */\n\t"\
		"movaps	%%xmm2,    (%%eax)	/* o0r */\n\t"\
		"movaps	%%xmm1,    (%%edx)	/* o3r */\n\t"\
		"movaps	%%xmm3,0x10(%%eax)	/* o0r */\n\t"\
		"movaps	%%xmm0,0x10(%%ecx)	/* o2i */\n\t"\
		"\n\t"\
		:					/* outputs: none */\
		: [__r0] "m" (Xr0)	/* All inputs from memory addresses here */\
		 ,[__o0] "m" (Xo0)\
		 ,[__o1] "m" (Xo1)\
		 ,[__o2] "m" (Xo2)\
		 ,[__o3] "m" (Xo3)\
		 ,[__o4] "m" (Xo4)\
		 ,[__o5] "m" (Xo5)\
		 ,[__o6] "m" (Xo6)\
		 ,[__o7] "m" (Xo7)\
		 ,[__isrt2] "m" (Xisrt2)\
		: "eax","ebx","ecx","edx","edi","esi"		/* Clobbered registers */\
	);\
	}

	#define	SSE2_RADIX8_DIT_0TWIDDLE(Xin0,Xin1,Xin2,Xin3,Xin4,Xin5,Xin6,Xin7, Xout, Xisrt2)\
	{\
	__asm__ volatile (\
		"movl	%[__in4],%%eax	\n\t"\
		"movl	%[__in5],%%ebx	\n\t"\
		"movl	%[__in6],%%ecx	\n\t"\
		"movl	%[__in7],%%edx	\n\t"\
		"movl	%[__out],%%esi	\n\t"\
		"movaps	    (%%eax),%%xmm0			\n\t"\
		"movaps	0x10(%%eax),%%xmm1			\n\t"\
		"movaps	%%xmm0,%%xmm2				\n\t"\
		"movaps	%%xmm1,%%xmm3				\n\t"\
		"addpd	    (%%ebx),%%xmm2			\n\t"\
		"addpd	0x10(%%ebx),%%xmm3			\n\t"\
		"subpd	    (%%ebx),%%xmm0			\n\t"\
		"subpd	0x10(%%ebx),%%xmm1			\n\t"\
		"\n\t"\
		"movaps	    (%%ecx),%%xmm4			\n\t"\
		"movaps	0x10(%%ecx),%%xmm5			\n\t"\
		"movaps	%%xmm4,%%xmm6				\n\t"\
		"movaps	%%xmm5,%%xmm7				\n\t"\
		"addpd	    (%%edx),%%xmm6			\n\t"\
		"addpd	0x10(%%edx),%%xmm7			\n\t"\
		"subpd	    (%%edx),%%xmm4			\n\t"\
		"subpd	0x10(%%edx),%%xmm5			\n\t"\
		"/* Copy t6r,i into main-array slot add6 */\n\t"\
		"movaps	%%xmm6,0xc0(%%esi)			\n\t"\
		"movaps	%%xmm7,0xd0(%%esi)			\n\t"\
		"/* Copy t7r,i into main-array slot add7 */\n\t"\
		"movaps	%%xmm4,0xe0(%%esi)			\n\t"\
		"movaps	%%xmm5,0xf0(%%esi)			\n\t"\
		"\n\t"\
		"/* Move outputs t5r,i into a(j1,j2+p5), first doing the addsub and mul by ISRT2: */\n\t"\
		"/* Move outputs t7r,i into a(j1,j2+p7), first doing the addsub and mul by ISRT2: */\n\t"\
		"\n\t"\
		"addpd	%%xmm2,%%xmm6				\n\t"\
		"addpd	%%xmm3,%%xmm7				\n\t"\
		"subpd	0xc0(%%esi),%%xmm2			\n\t"\
		"subpd	0xd0(%%esi),%%xmm3			\n\t"\
		"/* Move t6r,i into a(j1,j2+p6) */	\n\t"\
		"movaps	%%xmm2,0xc0(%%esi)			\n\t"\
		"movaps	%%xmm3,0xd0(%%esi)			\n\t"\
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
		"movl	%[__isrt2],%%edi			\n\t"\
		"movaps	(%%edi),%%xmm1	/* isrt2 */	\n\t"\
		"subpd	%%xmm3,%%xmm2				\n\t"\
		"mulpd	%%xmm1,%%xmm5				\n\t"\
		"mulpd	%%xmm1,%%xmm2				\n\t"\
		"movaps	%%xmm0,%%xmm3				\n\t"\
		"\n\t"\
		"movaps	%%xmm5,0xa0(%%esi)			\n\t"\
		"movaps	%%xmm4,%%xmm5				\n\t"\
		"addpd	%%xmm4,%%xmm0				\n\t"\
		"movaps	%%xmm2,0xb0(%%esi)			\n\t"\
		"subpd	%%xmm5,%%xmm3				\n\t"\
		"mulpd	%%xmm1,%%xmm0				\n\t"\
		"mulpd	%%xmm1,%%xmm3				\n\t"\
		"movaps	%%xmm0,0xe0(%%esi)			\n\t"\
		"movaps	%%xmm3,0xf0(%%esi)			\n\t"\
		"\n\t"\
		"/**************** 1st of the 2 length-4 subtransforms... **************/\n\t"\
		"movl	%[__in0],%%eax	\n\t"\
		"movl	%[__in1],%%ebx	\n\t"\
		"movl	%[__in2],%%ecx	\n\t"\
		"movl	%[__in3],%%edx	\n\t"\
		"\n\t"\
		"movaps	    (%%eax),%%xmm0			\n\t"\
		"movaps	0x10(%%eax),%%xmm1			\n\t"\
		"movaps	%%xmm0,%%xmm2				\n\t"\
		"movaps	%%xmm1,%%xmm3				\n\t"\
		"addpd	    (%%ebx),%%xmm2			\n\t"\
		"addpd	0x10(%%ebx),%%xmm3			\n\t"\
		"subpd	    (%%ebx),%%xmm0			\n\t"\
		"subpd	0x10(%%ebx),%%xmm1			\n\t"\
		"\n\t"\
		"/* Move   t4r,i into temp(j1+p0) in anticipation of final outputs (t0+t4)r,i which will go there: */\n\t"\
		"movaps	%%xmm6,    (%%esi)			\n\t"\
		"movaps	%%xmm7,0x10(%%esi)			\n\t"\
		"addpd	%%xmm6,%%xmm6				\n\t"\
		"addpd	%%xmm7,%%xmm7				\n\t"\
		"/* Move 2*t4r,i into temp(j1+p4) in anticipation of final outputs (t0-t4)r,i which will go there: */\n\t"\
		"movaps	%%xmm6,0x80(%%esi)			\n\t"\
		"movaps	%%xmm7,0x90(%%esi)			\n\t"\
		"\n\t"\
		"movaps	    (%%ecx),%%xmm4			\n\t"\
		"movaps	0x10(%%ecx),%%xmm5			\n\t"\
		"movaps	%%xmm4,%%xmm6				\n\t"\
		"movaps	%%xmm5,%%xmm7				\n\t"\
		"addpd	    (%%edx),%%xmm6			\n\t"\
		"addpd	0x10(%%edx),%%xmm7			\n\t"\
		"subpd	    (%%edx),%%xmm4			\n\t"\
		"subpd	0x10(%%edx),%%xmm5			\n\t"\
		"/* Copy t5,6 into temp(j1+p2) */	\n\t"\
		"movaps	%%xmm6,0x40(%%esi)			\n\t"\
		"movaps	%%xmm7,0x50(%%esi)			\n\t"\
		"addpd	%%xmm2,%%xmm6				\n\t"\
		"addpd	%%xmm3,%%xmm7				\n\t"\
		"subpd	0x40(%%esi),%%xmm2			\n\t"\
		"subpd	0x50(%%esi),%%xmm3			\n\t"\
		"\n\t"\
		"/* Compute and dump first 2 outputs now, in order to free up 2 registers: */\n\t"\
		"addpd	    (%%esi),%%xmm6			\n\t"\
		"addpd	0x10(%%esi),%%xmm7			\n\t"\
		"movaps	%%xmm6,    (%%esi)			\n\t"\
		"movaps	%%xmm7,0x10(%%esi)			\n\t"\
		"\n\t"\
		"subpd	0x80(%%esi),%%xmm6			\n\t"\
		"subpd	0x90(%%esi),%%xmm7			\n\t"\
		"movaps	%%xmm6,0x80(%%esi)			\n\t"\
		"movaps	%%xmm7,0x90(%%esi)			\n\t"\
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
		"addpd	0xd0(%%esi),%%xmm2			\n\t"\
		"subpd	0xc0(%%esi),%%xmm3			\n\t"\
		"subpd	0xd0(%%esi),%%xmm6			\n\t"\
		"addpd	0xc0(%%esi),%%xmm7			\n\t"\
		"movaps	%%xmm2,0xc0(%%esi)			\n\t"\
		"movaps	%%xmm3,0xd0(%%esi)			\n\t"\
		"movaps	%%xmm6,0x40(%%esi)			\n\t"\
		"movaps	%%xmm7,0x50(%%esi)			\n\t"\
		"\n\t"\
		"movaps	%%xmm5,%%xmm2				\n\t"\
		"movaps	%%xmm0,%%xmm6				\n\t"\
		"movaps	%%xmm1,%%xmm3				\n\t"\
		"movaps	%%xmm4,%%xmm7				\n\t"\
		"addpd	0xa0(%%esi),%%xmm5			\n\t"\
		"subpd	0xf0(%%esi),%%xmm0			\n\t"\
		"subpd	0xb0(%%esi),%%xmm1			\n\t"\
		"subpd	0xe0(%%esi),%%xmm4			\n\t"\
		"subpd	0xa0(%%esi),%%xmm2			\n\t"\
		"addpd	0xf0(%%esi),%%xmm6			\n\t"\
		"addpd	0xb0(%%esi),%%xmm3			\n\t"\
		"addpd	0xe0(%%esi),%%xmm7			\n\t"\
		"movaps	%%xmm5,0xe0(%%esi)			\n\t"\
		"movaps	%%xmm0,0xa0(%%esi)			\n\t"\
		"movaps	%%xmm1,0xf0(%%esi)			\n\t"\
		"movaps	%%xmm4,0xb0(%%esi)			\n\t"\
		"movaps	%%xmm2,0x60(%%esi)			\n\t"\
		"movaps	%%xmm6,0x20(%%esi)			\n\t"\
		"movaps	%%xmm3,0x70(%%esi)			\n\t"\
		"movaps	%%xmm7,0x30(%%esi)			\n\t"\
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
		: "eax","ebx","ecx","edx","edi","esi"		/* Clobbered registers */\
	);\
	}

#endif	/* sse2_macro_gcc_h_included */

