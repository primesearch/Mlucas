/*******************************************************************************
*                                                                              *
*   (C) 1997-2013 by Ernst W. Mayer.                                           *
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
		"movl	%[__cB]		,%%esi\n\t"\
		"movl	%[__cAmB]	,%%ecx\n\t"\
		"movl	%[__cApB]	,%%edx\n\t"\
		"movaps	    (%%eax),%%xmm0\n\t"\
		"movaps	0x10(%%eax),%%xmm2\n\t"\
		"movaps	    (%%esi),%%xmm4\n\t"\
		"movaps	0x10(%%esi),%%xmm5\n\t"\
		"movaps	%%xmm0,%%xmm1\n\t"\
		"movaps	%%xmm2,%%xmm3\n\t"\
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
		: "cc","memory","eax","esi","ecx","edx","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5"		/* Clobbered registers */\
	);\
	}

	#define PAIR_SQUARE_4_SSE2(XtAr, XtBr, XtCr, XtDr, Xc, Xs, Xforth)\
	{\
	__asm__ volatile (\
		/*   calculate cross-product terms...\
			__rt=__tAr* ~tDr+__tAi* ~tDi; __rt=__rt+__rt;\
			__it=__tAi* ~tDr-__tAr* ~tDi; __it=__it+__it;\
		*/\
		"movl	%[__tDr]	,%%edx\n\t"\
		"movl	%[__tAr]	,%%eax\n\t"\
		"movaps	    (%%edx)	,%%xmm6		/* tDr */\n\t"\
		"movaps	0x10(%%edx)	,%%xmm7		/* tDi */\n\t"\
		"movaps	    (%%eax)	,%%xmm0		/* tAr */\n\t"\
		"movaps	0x10(%%eax)	,%%xmm3		/* tAi */\n\t"\
		"shufpd	$1	,%%xmm6	,%%xmm6	/*~tDr */\n\t"\
		"shufpd	$1	,%%xmm7	,%%xmm7	/*~tDi */\n\t"\
		"movaps	    (%%eax)	,%%xmm2		/* cpy tAr */\n\t"\
		"movaps	0x10(%%eax)	,%%xmm1		/* cpy tAi */\n\t"\
		"mulpd	%%xmm6		,%%xmm0	/* tAr*~tDr */\n\t"\
		"mulpd	%%xmm7		,%%xmm3	/* tAi*~tDi */\n\t"\
		"mulpd	%%xmm6		,%%xmm1	/* tAi*~tDr */\n\t"\
		"mulpd	%%xmm7		,%%xmm2	/* tAr*~tDi */\n\t"\
		"addpd	%%xmm3		,%%xmm0	/* rt */\n\t"\
		"subpd	%%xmm2		,%%xmm1	/* it */\n\t"\
		"addpd	%%xmm0		,%%xmm0	/* rt=rt+rt */\n\t"\
		"addpd	%%xmm1		,%%xmm1	/* it=it+it; xmm2-7 free */\n\t"\
		/*\
			__st=__tBr* ~tCr+__tBi* ~tCi; __st=__st+__st;\
			__jt=__tBi* ~tCr-__tBr* ~tCi; __jt=__jt+__jt;\
		*/\
		"movl	%[__tCr]	,%%ecx\n\t"\
		"movl	%[__tBr]	,%%esi\n\t"\
		"movaps	    (%%ecx)	,%%xmm6		/* tCr */\n\t"\
		"movaps	0x10(%%ecx)	,%%xmm7		/* tCi */\n\t"\
		"movaps	    (%%esi)	,%%xmm2		/* tBr */\n\t"\
		"movaps	0x10(%%esi)	,%%xmm5		/* tBi */\n\t"\
		"shufpd	$1	,%%xmm6	,%%xmm6	/*~tCr */\n\t"\
		"shufpd	$1	,%%xmm7	,%%xmm7	/*~tCi */\n\t"\
		"movaps	    (%%esi)	,%%xmm4		/* cpy tBr */\n\t"\
		"movaps	0x10(%%esi)	,%%xmm3		/* cpy tBi */\n\t"\
		"mulpd	%%xmm6		,%%xmm2	/* tBr*~tCr */\n\t"\
		"mulpd	%%xmm7		,%%xmm5	/* tBi*~tCi */\n\t"\
		"mulpd	%%xmm6		,%%xmm3	/* tBi*~tCr */\n\t"\
		"mulpd	%%xmm7		,%%xmm4	/* tBr*~tCi */\n\t"\
		"addpd	%%xmm5		,%%xmm2	/* st */\n\t"\
		"subpd	%%xmm4		,%%xmm3	/* jt */\n\t"\
		"addpd	%%xmm2		,%%xmm2	/* st=st+st */\n\t"\
		"addpd	%%xmm3		,%%xmm3	/* jt=jt+jt; xmm4-7 free */\n\t"\
	/*   now calculate square terms and __store back in the same temporaries:	*/\
	/*	__tmp=(__tAr+__tAi)*(__tAr-__tAi); __tAi=__tAr*__tAi; __tAi=__tAi+__tAi; __tAr=__tmp;	*/\
		"movaps	    (%%eax)	,%%xmm4		/* __tAr */\n\t"\
		"movaps	0x10(%%eax)	,%%xmm5		/* __tAi */\n\t"\
		"subpd	%%xmm5		,%%xmm4		/* (__tAr-__tAi) */\n\t"\
		"addpd	%%xmm5		,%%xmm5		/*      2*__tAi  */\n\t"\
		"addpd	%%xmm4		,%%xmm5		/* (__tAr+__tAi) */\n\t"\
		"mulpd	%%xmm5		,%%xmm4		/*>__tAr */\n\t"\
		"movaps	    (%%eax)	,%%xmm5		/* __tAr */\n\t"\
		"mulpd	0x10(%%eax)	,%%xmm5		/* __tAr*__tAi */\n\t"\
		"addpd	%%xmm5		,%%xmm5		/*>__tAi */\n\t"\
		"movaps	%%xmm4	,    (%%eax)	/* tmp store >__tAr */\n\t"\
		"movaps	%%xmm5	,0x10(%%eax)	/* tmp store >__tAi */\n\t"\
		"subpd	%%xmm4		,%%xmm0	/* rt-__tAr */\n\t"\
		"subpd	%%xmm5		,%%xmm1	/* it-__tAi; xmm4-7 free */\n\t"\
	/*	__tmp=(__tBr+__tBi)*(__tBr-__tBi); __tBi=__tBr*__tBi; __tBi=__tBi+__tBi; __tBr=__tmp;	*/\
	/* [Can be done in parallel with above segment] */\
		"movaps	    (%%esi)	,%%xmm6		/* __tBr */\n\t"\
		"movaps	0x10(%%esi)	,%%xmm7		/* __tBi */\n\t"\
		"subpd	%%xmm7		,%%xmm6		/* (__tBr-__tBi) */\n\t"\
		"addpd	%%xmm7		,%%xmm7		/*      2*__tBi  */\n\t"\
		"addpd	%%xmm6		,%%xmm7		/* (__tBr+__tBi) */\n\t"\
		"mulpd	%%xmm7		,%%xmm6		/*>__tBr */\n\t"\
		"movaps	    (%%esi)	,%%xmm7		/* __tBr */\n\t"\
		"mulpd	0x10(%%esi)	,%%xmm7		/* __tBr*__tBi */\n\t"\
		"addpd	%%xmm7		,%%xmm7		/*>__tBi */\n\t"\
		"movaps	%%xmm6	,    (%%esi)	/* tmp store >__tBr */\n\t"\
		"movaps	%%xmm7	,0x10(%%esi)	/* tmp store >__tBi */\n\t"\
		"subpd	%%xmm6		,%%xmm2	/* st-__tBr */\n\t"\
		"subpd	%%xmm7		,%%xmm3	/* jt-__tBi; xmm4-7 free */\n\t"\
	/*	__tmp=(__tDr+__tDi)*(__tDr-__tDi); __tDi=__tDr*__tDi; __tDi=__tDi+__tDi; __tDr=__tmp;	*/\
		"movaps	    (%%edx)	,%%xmm4		/* __tDr */\n\t"\
		"movaps	0x10(%%edx)	,%%xmm5		/* __tDi */\n\t"\
		"subpd	%%xmm5		,%%xmm4		/* (__tDr-__tDi) */\n\t"\
		"addpd	%%xmm5		,%%xmm5		/*      2*__tDi  */\n\t"\
		"addpd	%%xmm4		,%%xmm5		/* (__tDr+__tDi) */\n\t"\
		"mulpd	%%xmm5		,%%xmm4		/*>__tDr */\n\t"\
		"movaps	    (%%edx)	,%%xmm5		/* __tDr */\n\t"\
		"mulpd	0x10(%%edx)	,%%xmm5		/* __tDr*__tDi */\n\t"\
		"addpd	%%xmm5		,%%xmm5		/*>__tDi */\n\t"\
		"movaps	%%xmm4	,    (%%edx)	/* tmp store ~tDr */\n\t"\
		"movaps	%%xmm5	,0x10(%%edx)	/* tmp store ~tDi */\n\t"\
		"shufpd	$1	,%%xmm4	,%%xmm4	/*~tDr */\n\t"\
		"shufpd	$1	,%%xmm5	,%%xmm5	/*~tDi */\n\t"\
		"subpd	%%xmm4		,%%xmm0	/* rt-__tAr- ~tDr */\n\t"\
		"addpd	%%xmm5		,%%xmm1	/* it-__tAi+ ~tDi; xmm4-7 free */\n\t"\
	/*	__tmp=(__tCr+__tCi)*(__tCr-__tCi); __tCi=__tCr*__tCi; __tCi=__tCi+__tCi; __tCr=__tmp;	*/\
	/* [Can be done in parallel with above segment] */\
		"movaps	    (%%ecx)	,%%xmm6		/* __tCr */\n\t"\
		"movaps	0x10(%%ecx)	,%%xmm7		/* __tCi */\n\t"\
		"subpd	%%xmm7		,%%xmm6		/* (__tCr-__tCi) */\n\t"\
		"addpd	%%xmm7		,%%xmm7		/*      2*__tCi  */\n\t"\
		"addpd	%%xmm6		,%%xmm7		/* (__tCr+__tCi) */\n\t"\
		"mulpd	%%xmm7		,%%xmm6		/*>__tCr */\n\t"\
		"movaps	    (%%ecx)	,%%xmm7		/* __tCr */\n\t"\
		"mulpd	0x10(%%ecx)	,%%xmm7		/* __tCr*__tCi */\n\t"\
		"addpd	%%xmm7		,%%xmm7		/*>__tCi */\n\t"\
		"movaps	%%xmm6	,    (%%ecx)	/* tmp store ~tCr */\n\t"\
		"movaps	%%xmm7	,0x10(%%ecx)	/* tmp store ~tCi */\n\t"\
		"shufpd	$1	,%%xmm6	,%%xmm6	/*~tCr */\n\t"\
		"shufpd	$1	,%%xmm7	,%%xmm7	/*~tCi */\n\t"\
		"subpd	%%xmm6		,%%xmm2	/* st-__tBr- ~tCr */\n\t"\
		"addpd	%%xmm7		,%%xmm3	/* jt-__tBi+ ~tCi; xmm4-7 free */\n\t"\
	/*\
		__tmp=((1.0+__c)*__rt-__s*__it)*0.25;\
		__it =((1.0+__c)*__it+__s*__rt)*0.25;	__rt=__tmp;\
	*/\
	/* [Can be done in parallel with above segment] */\
		"movl	%[__c]		,%%eax\n\t"\
		"movl	%[__s]		,%%esi\n\t"\
		"movl	%[__forth]	,%%edx\n\t"\
		"movaps	%%xmm0		,%%xmm4		/* cpy rt */\n\t"\
		"movaps	%%xmm1		,%%xmm5		/* cpy it */\n\t"\
		"mulpd	(%%eax)		,%%xmm0		/* c*rt */\n\t"\
		"mulpd	(%%eax)		,%%xmm1		/* c*it */\n\t"\
		"addpd	%%xmm4		,%%xmm0		/* (c+1.0)*rt */\n\t"\
		"addpd	%%xmm5		,%%xmm1		/* (c+1.0)*it */\n\t"\
		"mulpd	(%%esi)		,%%xmm4		/* s*rt */\n\t"\
		"mulpd	(%%esi)		,%%xmm5		/* s*it */\n\t"\
		"subpd	%%xmm5		,%%xmm0		/* (c+1.0)*rt-s*it */\n\t"\
		"addpd	%%xmm4		,%%xmm1		/* (c+1.0)*it+s*rt; xmm4,5 free */\n\t"\
		"mulpd	(%%edx)		,%%xmm0	/* -rt Both of these inherit the sign flip [w.r.to the non-SSE2 PAIR_SQUARE_4 macro] */\n\t"\
		"mulpd	(%%edx)		,%%xmm1	/* -it that resulted from the in-place-friendlier (rt-__tAr- ~tDr) reordering above. */\n\t"\
	/*\
		__tmp=((1.0-__s)*__st-__c*__jt)*0.25;\
		__jt =((1.0-__s)*__jt+__c*__st)*0.25	__st=__tmp;\
	*/\
	/* [Can be done in parallel wjth above segment] */\
		"movaps	%%xmm2		,%%xmm6		/* cpy st */\n\t"\
		"movaps	%%xmm3		,%%xmm7		/* cpy jt */\n\t"\
		"mulpd	(%%esi)		,%%xmm2		/* s*st */\n\t"\
		"mulpd	(%%esi)		,%%xmm3		/* s*jt */\n\t"\
		"subpd	%%xmm6		,%%xmm2		/* (s-1.0)*st, note sign flip! */\n\t"\
		"subpd	%%xmm7		,%%xmm3		/* (s-1.0)*jt, note sign flip! */\n\t"\
		"mulpd	(%%eax)		,%%xmm6		/* c*st */\n\t"\
		"mulpd	(%%eax)		,%%xmm7		/* c*jt */\n\t"\
		"addpd	%%xmm7		,%%xmm2		/* -[(1.0-s)*st-c*jt] */\n\t"\
		"subpd	%%xmm6		,%%xmm3		/* -[(1.0-s)*jt+c*st]; xmm6,7 free */\n\t"\
		"mulpd	(%%edx)		,%%xmm2	/* +st Sign flip due to (s-1.0) reordering here */\n\t"\
		"mulpd	(%%edx)		,%%xmm3	/* +jt cancels earlier one due to in-place-friendlier (st-__tBr- ~tCr) reordering above. */\n\t"\
	/*...and now complete and store the results. We flip the signs on st and jt here to undo the above -st,-jt negations. */\
	/*	__tAr = (__tAr+__rt);\
		__tAi = (__tAi+__it);\
		__tBr = (__tBr-__st);\
		__tBi = (__tBi-__jt);\
	*/\
		"movl	%[__tAr]	,%%eax\n\t"\
		"movl	%[__tBr]	,%%esi\n\t"\
		"movaps	    (%%eax)	,%%xmm4		/* __tAr */\n\t"\
		"movaps	0x10(%%eax)	,%%xmm5		/* __tAi */\n\t"\
		"movaps	    (%%esi)	,%%xmm6		/* __tBr */\n\t"\
		"movaps	0x10(%%esi)	,%%xmm7		/* __tBi */\n\t"\
		"addpd	%%xmm0		,%%xmm4		/* (__tAr+__rt) */\n\t"\
		"addpd	%%xmm1		,%%xmm5		/* (__tAi+__it) */\n\t"\
		"subpd	%%xmm2		,%%xmm6		/* (__tBr-__st) */\n\t"\
		"subpd	%%xmm3		,%%xmm7		/* (__tBi-__jt) */\n\t"\
		"movaps	%%xmm4	,    (%%eax)	/* store >__tAr */\n\t"\
		"movaps	%%xmm5	,0x10(%%eax)	/* store >__tAi */\n\t"\
		"movaps	%%xmm6	,    (%%esi)	/* store >__tBr */\n\t"\
		"movaps	%%xmm7	,0x10(%%esi)	/* store >__tBi */\n\t"\
	/*...N-j terms are as above, but with the replacements: __tAr<--> ~tDr, __tAi<--> ~tDi, __it|-->-__it. */\
	/*	__tDr = (__tDr+ ~rt);\
		__tDi = (__tDi- ~it);\
		__tCr = (__tCr- ~st);\
		__tCi = (__tCi+ ~jt);\
	*/\
		"movl	%[__tCr]	,%%ecx\n\t"\
		"movl	%[__tDr]	,%%edx\n\t"\
		"shufpd	$1	,%%xmm0	,%%xmm0		/* ~rt */\n\t"\
		"shufpd	$1	,%%xmm1	,%%xmm1		/* ~it */\n\t"\
		"shufpd	$1	,%%xmm2	,%%xmm2		/* ~st */\n\t"\
		"shufpd	$1	,%%xmm3	,%%xmm3		/* ~jt */\n\t"\
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
		: "cc","memory","eax","esi","ecx","edx","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
	);\
	}

	#define SSE2_RADIX_03_DFT(Xi0,Xi1,Xi2, Xcc1, Xo0,Xo1,Xo2)\
	{\
	__asm__ volatile (\
	"pushl %%ebx	\n\t"/* Explicit save/restore of PIC register */\
		"movl	%[__i0],%%eax		\n\t"\
		"movl	%[__i1],%%ebx		\n\t"\
		"movl	%[__i2],%%ecx		\n\t"\
		"movl	%[__cc1],%%edx		\n\t"\
		"movaps	    (%%ebx),%%xmm2	\n\t"\
		"movaps	0x10(%%ebx),%%xmm3	\n\t"\
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"movaps	    (%%ecx),%%xmm6	\n\t"\
		"movaps	0x10(%%ecx),%%xmm7	\n\t"\
		"movaps	%%xmm2,%%xmm4		\n\t"\
		"movaps	%%xmm3,%%xmm5		\n\t"\
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
		"mulpd	%%xmm6,%%xmm2		\n\t"\
		"mulpd	%%xmm6,%%xmm3		\n\t"\
		"mulpd	%%xmm7,%%xmm4		\n\t"\
		"mulpd	%%xmm7,%%xmm5		\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"movaps	%%xmm2,%%xmm0		\n\t"\
		"movaps	%%xmm3,%%xmm1		\n\t"\
		"subpd	%%xmm5,%%xmm2		\n\t"\
		"addpd	%%xmm4,%%xmm3		\n\t"\
		"addpd	%%xmm5,%%xmm0		\n\t"\
		"subpd	%%xmm4,%%xmm1		\n\t"\
		"movaps	%%xmm2,     (%%ebx)	\n\t"\
		"movaps	%%xmm3,0x010(%%ebx)	\n\t"\
		"movaps	%%xmm0,     (%%ecx)	\n\t"\
		"movaps	%%xmm1,0x010(%%ecx)	\n\t"\
	"popl %%ebx	\n\t"\
		:					/* outputs: none */\
		: [__i0] "m" (Xi0)	/* All inputs from memory addresses here */\
		 ,[__i1] "m" (Xi1)\
		 ,[__i2] "m" (Xi2)\
		 ,[__cc1] "m" (Xcc1)\
		 ,[__o0] "m" (Xo0)\
		 ,[__o1] "m" (Xo1)\
		 ,[__o2] "m" (Xo2)\
		: "cc","memory","eax",/*"ebx",*/"ecx","edx","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"	/* Clobbered registers */\
	);\
	}

	#define SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(Xadd0, Xadd1, Xadd2, Xadd3, Xtmp, Xstride)\
	{\
	__asm__ volatile (\
	"pushl %%ebx	\n\t"/* Explicit save/restore of PIC register */\
		"movl	%[__tmp]   ,%%eax	\n\t"\
		"movl	%[__stride],%%esi	\n\t"\
		"movl	%%eax,%%ebx			\n\t"\
		"addl	%%esi,%%ebx			\n\t"/* add_in1  */\
		"shll	$1,%%esi			\n\t"/* stride*2 */\
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	    (%%ebx),%%xmm2	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"movaps	0x10(%%ebx),%%xmm3	\n\t"\
		"movaps	    (%%eax),%%xmm4	\n\t"\
		"movaps	    (%%ebx),%%xmm6	\n\t"\
		"movaps	0x10(%%eax),%%xmm5	\n\t"\
		"movaps	0x10(%%ebx),%%xmm7	\n\t"\
		"addl	%%esi,%%eax			\n\t"/* add_in2  */\
		"addl	%%esi,%%ebx			\n\t"/* add_in3  */\
		"addpd	    (%%eax),%%xmm0	\n\t"\
		"addpd	    (%%ebx),%%xmm2	\n\t"\
		"addpd	0x10(%%eax),%%xmm1	\n\t"\
		"addpd	0x10(%%ebx),%%xmm3	\n\t"\
		"subpd	    (%%eax),%%xmm4	\n\t"\
		"subpd	    (%%ebx),%%xmm6	\n\t"\
		"subpd	0x10(%%eax),%%xmm5	\n\t"\
		"subpd	0x10(%%ebx),%%xmm7	\n\t"\
	/* Finish radix-4 butterfly and store results into main-array slots: */\
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
	"popl %%ebx	\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__add1] "m" (Xadd1)\
		 ,[__add2] "m" (Xadd2)\
		 ,[__add3] "m" (Xadd3)\
		 ,[__tmp] "m" (Xtmp)\
		 ,[__stride] "e" (Xstride)\
		: "cc","memory","eax",/*"ebx",*/"ecx","edx","esi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"	/* Clobbered registers */\
	);\
	}

	/* DIF radix-4 subconvolution, sans twiddles, inputs in __i0-3, outputs in __o0-3, possibly coincident with inputs: */
	#define SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_E(Xi0,Xi1,Xi2,Xi3, Xo0,Xo1,Xo2,Xo3)\
	{\
	__asm__ volatile (\
	"pushl %%ebx	\n\t"/* Explicit save/restore of PIC register */\
		"movl	%[__i0],%%eax		\n\t"\
		"movl	%[__i1],%%ebx		\n\t"\
		"movl	%[__i2],%%ecx		\n\t"\
		"movl	%[__i3],%%edx		\n\t"\
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	    (%%ebx),%%xmm4	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"movaps	0x10(%%ebx),%%xmm5	\n\t"\
		"movaps	%%xmm0,%%xmm2	\n\t"\
		"movaps	%%xmm4,%%xmm6	\n\t"\
		"movaps	%%xmm1,%%xmm3	\n\t"\
		"movaps	%%xmm5,%%xmm7	\n\t"\
		"addpd	    (%%ecx),%%xmm0	\n\t"\
		"addpd	    (%%edx),%%xmm4	\n\t"\
		"addpd	0x10(%%ecx),%%xmm1	\n\t"\
		"addpd	0x10(%%edx),%%xmm5	\n\t"\
		"subpd	    (%%ecx),%%xmm2	\n\t"\
		"subpd	    (%%edx),%%xmm6	\n\t"\
		"subpd	0x10(%%ecx),%%xmm3	\n\t"\
		"subpd	0x10(%%edx),%%xmm7	\n\t"\
	/* Finish radix-4 butterfly and store results into main-array slots: */\
		"movl	%[__o0],%%eax		\n\t"\
		"movl	%[__o1],%%ebx		\n\t"\
		"movl	%[__o2],%%ecx		\n\t"\
		"movl	%[__o3],%%edx		\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm2		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"subpd	%%xmm6,%%xmm3		\n\t"\
		"movaps	%%xmm0,    (%%ebx)	\n\t"\
		"movaps	%%xmm2,    (%%ecx)	\n\t"\
		"movaps	%%xmm1,0x10(%%ebx)	\n\t"\
		"movaps	%%xmm3,0x10(%%edx)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"addpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	%%xmm4,    (%%eax)	\n\t"\
		"movaps	%%xmm7,    (%%edx)	\n\t"\
		"movaps	%%xmm5,0x10(%%eax)	\n\t"\
		"movaps	%%xmm6,0x10(%%ecx)	\n\t"\
	"popl %%ebx	\n\t"\
		:					/* outputs: none */\
		: [__i0] "m" (Xi0)	/* All inputs from memory addresses here */\
		 ,[__i1] "m" (Xi1)\
		 ,[__i2] "m" (Xi2)\
		 ,[__i3] "m" (Xi3)\
		 ,[__o0] "m" (Xo0)\
		 ,[__o1] "m" (Xo1)\
		 ,[__o2] "m" (Xo2)\
		 ,[__o3] "m" (Xo3)\
		: "cc","memory","eax",/*"ebx",*/"ecx","edx","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"	/* Clobbered registers */\
	);\
	}

	#define SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(Xadd0, Xadd1, Xadd2, Xadd3, Xtmp, Xstride)\
	{\
	__asm__ volatile (\
	"pushl %%ebx	\n\t"/* Explicit save/restore of PIC register */\
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
	/* Finish radix-4 butterfly and store results into temp-array slots: */\
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
	"popl %%ebx	\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__add1] "m" (Xadd1)\
		 ,[__add2] "m" (Xadd2)\
		 ,[__add3] "m" (Xadd3)\
		 ,[__tmp] "m" (Xtmp)\
		 ,[__stride] "e" (Xstride)\
		: "cc","memory","eax",/*"ebx",*/"ecx","edx","esi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"	/* Clobbered registers */\
	);\
	}

	/* DIT radix-4 subconvolution, sans twiddles, inputs in __i0-3, outputs in __o0-3, possibly coincident with inputs: */
	#define SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_E(Xi0,Xi1,Xi2,Xi3, Xo0,Xo1,Xo2,Xo3)\
	{\
	__asm__ volatile (\
	"pushl %%ebx	\n\t"/* Explicit save/restore of PIC register */\
		"movl	%[__i0],%%eax		\n\t"\
		"movl	%[__i1],%%ebx		\n\t"\
		"movl	%[__i2],%%ecx		\n\t"\
		"movl	%[__i3],%%edx		\n\t"\
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	    (%%ecx),%%xmm4	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"movaps	0x10(%%ecx),%%xmm5	\n\t"\
		"movaps	%%xmm0,%%xmm2			\n\t"\
		"movaps	%%xmm4,%%xmm6			\n\t"\
		"movaps	%%xmm1,%%xmm3			\n\t"\
		"movaps	%%xmm5,%%xmm7			\n\t"\
		"addpd	    (%%ebx),%%xmm0	\n\t"\
		"addpd	    (%%edx),%%xmm4	\n\t"\
		"addpd	0x10(%%ebx),%%xmm1	\n\t"\
		"addpd	0x10(%%edx),%%xmm5	\n\t"\
		"subpd	    (%%ebx),%%xmm2	\n\t"\
		"subpd	    (%%edx),%%xmm6	\n\t"\
		"subpd	0x10(%%ebx),%%xmm3	\n\t"\
		"subpd	0x10(%%edx),%%xmm7	\n\t"\
	/* Finish radix-4 butterfly and store results into output-array slots: */\
		"movl	%[__o0],%%eax		\n\t"\
		"movl	%[__o1],%%ebx		\n\t"\
		"movl	%[__o2],%%ecx		\n\t"\
		"movl	%[__o3],%%edx		\n\t"\
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
	"popl %%ebx	\n\t"\
		:					/* outputs: none */\
		: [__i0] "m" (Xi0)	/* All inputs from memory addresses here */\
		 ,[__i1] "m" (Xi1)\
		 ,[__i2] "m" (Xi2)\
		 ,[__i3] "m" (Xi3)\
		 ,[__o0] "m" (Xo0)\
		 ,[__o1] "m" (Xo1)\
		 ,[__o2] "m" (Xo2)\
		 ,[__o3] "m" (Xo3)\
		: "cc","memory","eax",/*"ebx",*/"ecx","edx","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"	/* Clobbered registers */\
	);\
	}

	#define SSE2_RADIX_05_DFT_0TWIDDLE(Xi0,Xi1,Xi2,Xi3,Xi4, Xcc1, Xo0,Xo1,Xo2,Xo3,Xo4)\
	{\
	__asm__ volatile (\
	"pushl %%ebx	\n\t"/* Explicit save/restore of PIC register */\
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
	"popl %%ebx	\n\t"\
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
		: "cc","memory","eax",/*"ebx",*/"ecx","edx","edi","esi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"	/* Clobbered registers */\
	);\
	}

	/*...Radix-7 DFT: Inputs in memlocs __i0-6, outputs into __o0-6, possibly coincident with inputs:\ */\
	#define SSE2_RADIX_07_DFT(Xi0,Xi1,Xi2,Xi3,Xi4,Xi5,Xi6, Xcc, Xo0,Xo1,Xo2,Xo3,Xo4,Xo5,Xo6)\
	{\
	__asm__ volatile (\
	"pushl %%ebx	\n\t"/* Explicit save/restore of PIC register */\
		"movl	%[__i1],%%eax	\n\t"\
		"movl	%[__i2],%%ebx	\n\t"\
		"movl	%[__i3],%%ecx	\n\t"\
		"movl	%[__i4],%%edx	\n\t"\
		"movl	%[__i5],%%esi	\n\t"\
		"movl	%[__i6],%%edi	\n\t"\
		"movaps	(%%eax),%%xmm6	\n\t"/* A1r */\
		"movaps	(%%edi),%%xmm1	\n\t"/* A6r */\
		"movaps	(%%ebx),%%xmm5	\n\t"/* A2r */\
		"movaps	(%%esi),%%xmm2	\n\t"/* A5r */\
		"movaps	(%%ecx),%%xmm4	\n\t"/* A3r */\
		"movaps	(%%edx),%%xmm3	\n\t"/* A4r */\
		"movl	%[__i0],%%ebx	\n\t"\
		"subpd	%%xmm1,%%xmm6	\n\t"/* t6r = A1r-A6r */\
		"addpd	%%xmm1,%%xmm1	\n\t"/*         2*A6r */\
		"addpd	%%xmm6,%%xmm1	\n\t"/* t1r = A1r+A6r */\
		"subpd	%%xmm2,%%xmm5	\n\t"/* t5r = A2r-A5r */\
		"addpd	%%xmm2,%%xmm2	\n\t"/*         2*A5r */\
		"addpd	%%xmm5,%%xmm2	\n\t"/* t2r = A2r+A5r */\
		"movaps	(%%ebx),%%xmm0	\n\t"/* Ar0           */\
		"subpd	%%xmm3,%%xmm4	\n\t"/* t4r = A3r-A4r */\
		"addpd	%%xmm3,%%xmm3	\n\t"/*         2*A4r */\
		"addpd	%%xmm4,%%xmm3	\n\t"/* t3r = A3r+A4r */\
		"movl	%[__o0],%%ecx	/* This might be same address as any of i0-i6 */\n\t"\
		"movl	%[__cc],%%esi		\n\t"\
		"movaps	%%xmm0,0x80(%%esi)	/* cpy t0 into scratch sincos slot */\n\t"\
		"movaps	%%xmm6,0x90(%%esi)	/* cpy t6 into scratch sincos slot */\n\t"\
		"addpd	%%xmm1,%%xmm0		/*~A0 = A0+t1 */\n\t"\
		"movaps	%%xmm5,%%xmm7		/* cpy t5 */\n\t"\
		"addpd	%%xmm2,%%xmm3		/*~t3 = t3+t2 */\n\t"\
		"subpd	%%xmm4,%%xmm5		/*~t5 = t5-t4 */\n\t"\
		"subpd	%%xmm2,%%xmm1		/*~t1 = t1-t2 */\n\t"\
		"subpd	%%xmm7,%%xmm6		/*~t6 = t6-t5 */\n\t"\
		"addpd	%%xmm2,%%xmm2		/* 2*t2 */\n\t"\
		"addpd	%%xmm7,%%xmm4		/*~t5 = t4+t5 */\n\t"\
		"addpd	%%xmm3,%%xmm0		/* B0 */\n\t"\
		"addpd	0x90(%%esi),%%xmm5	/* t3 = (t5-t4)+t6 */\n\t"\
		"subpd	%%xmm2,%%xmm3		/*~t2 =  (t2+t3) - 2*t2 = t3-t2 */\n\t"\
		"movaps	%%xmm4,%%xmm7		/* cpy t5 */\n\t"\
		"movaps	%%xmm0,    (%%ecx)	/* <-B0; %%xmm0 FREE */\n\t"\
		"subpd	%%xmm6,%%xmm4		/* t4 = ~t5-~t6 */\n\t"\
		"movaps	%%xmm1,%%xmm2		/* cpy ~t1 */\n\t"\
		"subpd	0x80(%%esi),%%xmm0	/* r = B0 - t0 */\n\t"\
		"mulpd	0x10(%%esi),%%xmm5	/*~t3 = t3*sx0 */\n\t"\
		"addpd	%%xmm3,%%xmm2		/* ~t1+~t2 */\n\t"\
		"mulpd	0x40(%%esi),%%xmm3	/* t2 = t2*cx2 */\n\t"\
		"mulpd	0x70(%%esi),%%xmm4	/*~t4 = t4*sx3 */\n\t"\
		"mulpd	0x20(%%esi),%%xmm1	/* t1 = t1*cx1 */\n\t"\
		"mulpd	0x30(%%esi),%%xmm6	/*~t6 = t6*sx1 */\n\t"\
		"mulpd	    (%%esi),%%xmm0	/* ~r = r*(cx0-1) */\n\t"\
		"mulpd	0x50(%%esi),%%xmm7	/*~t5 = t5*sx2 */\n\t"\
		"mulpd	0x60(%%esi),%%xmm2	/* t3 =(t1+t2)*cx3 */\n\t"\
		"addpd	    (%%ecx),%%xmm0	/* t0 =~r + B0 */\n\t"\
		"addpd	%%xmm4,%%xmm6		/*~t6 = t4+t6 */\n\t"\
		"subpd	%%xmm2,%%xmm1		/* tt = t1-t3 */\n\t"\
		"subpd	%%xmm7,%%xmm4		/*~t5 = t4-t5; %%xmm7 FREE */\n\t"\
		"subpd	%%xmm2,%%xmm3		/* t2 = t2-t3; %%xmm2 FREE */\n\t"\
		"movl	%[__o1],%%eax		\n\t"\
		"movl	%[__o2],%%ebx		\n\t"\
		"movl	%[__o3],%%ecx		\n\t"\
		"movl	%[__o4],%%edx		\n\t"\
		"movl	%[__o5],%%esi		\n\t"\
		"movl	%[__o6],%%edi		\n\t"\
		"movaps	%%xmm0,%%xmm2		/* cpy t0 */\n\t"\
		"movaps	%%xmm5,%%xmm7		/* cpy t3 */\n\t"\
		"addpd	%%xmm1,%%xmm0		/*~t0 = t0+tt */\n\t"\
		"addpd	%%xmm6,%%xmm5		/*~t3 = t3+t6 */\n\t"\
		"addpd	%%xmm3,%%xmm1		/*~tt = tt+t2 */\n\t"\
		"addpd	%%xmm4,%%xmm6		/*      t6+t5 */\n\t"\
		"addpd	%%xmm2,%%xmm3		/*~t2 = t2+t0 */\n\t"\
		"addpd	%%xmm7,%%xmm4		/*~t5 = t5+t3 */\n\t"\
		"subpd	%%xmm1,%%xmm2		/*~t1 = t0-tt-t2 */\n\t"\
		"subpd	%%xmm6,%%xmm7		/*~t4 = t3-t6-t5 */\n\t"\
		"movaps	%%xmm0,(%%eax)		\n\t"/* B1 <- t0 */\
		"movaps	%%xmm5,(%%edi)		\n\t"/* B6 <- t3 */\
		"movaps	%%xmm2,(%%ebx)		\n\t"/* B2 <- t1 */\
		"movaps	%%xmm7,(%%esi)		\n\t"/* B5 <- t4 */\
		"movaps	%%xmm3,(%%ecx)		\n\t"/* B3 <- t2 */\
		"movaps	%%xmm4,(%%edx)		\n\t"/* B4 <- t5 */\
	/* Imaginary Parts: */\
		"movl	%[__i1],%%eax		\n\t"\
		"movl	%[__i2],%%ebx		\n\t"\
		"movl	%[__i3],%%ecx		\n\t"\
		"movl	%[__i4],%%edx		\n\t"\
		"movl	%[__i5],%%esi		\n\t"\
		"movl	%[__i6],%%edi		\n\t"\
		"movaps	0x10(%%eax),%%xmm6	\n\t"/* A1i */\
		"movaps	0x10(%%edi),%%xmm1	\n\t"/* A6i */\
		"movaps	0x10(%%ebx),%%xmm5	\n\t"/* A2i */\
		"movaps	0x10(%%esi),%%xmm2	\n\t"/* A5i */\
		"movaps	0x10(%%ecx),%%xmm4	\n\t"/* A3i */\
		"movaps	0x10(%%edx),%%xmm3	\n\t"/* A4i */\
		"movl	%[__i0],%%ebx		\n\t"\
		"subpd	%%xmm1,%%xmm6		\n\t"/* t6i = A1i-A6i */\
		"addpd	%%xmm1,%%xmm1		\n\t"/*         2*A6i */\
		"addpd	%%xmm6,%%xmm1		\n\t"/* t1i = A1i+A6i */\
		"subpd	%%xmm2,%%xmm5		\n\t"/* t5i = A2i-A5i */\
		"addpd	%%xmm2,%%xmm2		\n\t"/*         2*A5i */\
		"addpd	%%xmm5,%%xmm2		\n\t"/* t2i = A2i+A5i */\
		"movaps	0x10(%%ebx),%%xmm0	\n\t"/* Ai0           */\
		"subpd	%%xmm3,%%xmm4		\n\t"/* t4i = A3i-A4i */\
		"addpd	%%xmm3,%%xmm3		\n\t"/*         2*A4i */\
		"addpd	%%xmm4,%%xmm3		\n\t"/* t3i = A3i+A4i */\
		"movl	%[__o0],%%ecx		\n\t"\
		"movl	%[__cc],%%esi		\n\t"\
		"movaps	%%xmm0,0x80(%%esi)	/* cpy t0 into scratch sincos slot */\n\t"\
		"movaps	%%xmm6,0x90(%%esi)	/* cpy t6 into scratch sincos slot */\n\t"\
		"addpd	%%xmm1,%%xmm0		/*~A0 = A0+t1 */\n\t"\
		"movaps	%%xmm5,%%xmm7		/* cpy t5 */\n\t"\
		"addpd	%%xmm2,%%xmm3		/*~t3 = t3+t2 */\n\t"\
		"subpd	%%xmm4,%%xmm5		/*~t5 = t5-t4 */\n\t"\
		"subpd	%%xmm2,%%xmm1		/*~t1 = t1-t2 */\n\t"\
		"subpd	%%xmm7,%%xmm6		/*~t6 = t6-t5 */\n\t"\
		"addpd	%%xmm2,%%xmm2		/* 2*t2 */\n\t"\
		"addpd	%%xmm7,%%xmm4		/*~t5 = t4+t5 */\n\t"\
		"addpd	%%xmm3,%%xmm0		/* B0 */\n\t"\
		"addpd	0x90(%%esi),%%xmm5	/* t3 = (t5-t4)+t6 */\n\t"\
		"subpd	%%xmm2,%%xmm3		/*~t2 =  (t2+t3) - 2*t2 = t3-t2 */\n\t"\
		"movaps	%%xmm4,%%xmm7		/* cpy t5 */\n\t"\
		"movaps	%%xmm0,0x10(%%ecx)	/* <-B0; %%xmm0 FREE */\n\t"\
		"subpd	%%xmm6,%%xmm4		/* t4 = ~t5-~t6 */\n\t"\
		"movaps	%%xmm1,%%xmm2		/* cpy ~t1 */\n\t"\
		"subpd	0x80(%%esi),%%xmm0	/* r = B0 - t0 */\n\t"\
		"mulpd	0x10(%%esi),%%xmm5	/*~t3 = t3*sx0 */\n\t"\
		"addpd	%%xmm3,%%xmm2		/* ~t1+~t2 */\n\t"\
		"mulpd	0x40(%%esi),%%xmm3	/* t2 = t2*cx2 */\n\t"\
		"mulpd	0x70(%%esi),%%xmm4	/*~t4 = t4*sx3 */\n\t"\
		"mulpd	0x20(%%esi),%%xmm1	/* t1 = t1*cx1 */\n\t"\
		"mulpd	0x30(%%esi),%%xmm6	/*~t6 = t6*sx1 */\n\t"\
		"mulpd	    (%%esi),%%xmm0/* ~r = r*(cx0-1) */\n\t"\
		"mulpd	0x50(%%esi),%%xmm7	/*~t5 = t5*sx2 */\n\t"\
		"mulpd	0x60(%%esi),%%xmm2	/* t3 =(t1+t2)*cx3 */\n\t"\
		"addpd	0x10(%%ecx),%%xmm0	/* t0 =~r + B0 */\n\t"\
		"addpd	%%xmm4,%%xmm6		/*~t6 = t4+t6 */\n\t"\
		"subpd	%%xmm2,%%xmm1		/* tt = t1-t3 */\n\t"\
		"subpd	%%xmm7,%%xmm4		/*~t5 = t4-t5; %%xmm7 FREE */\n\t"\
		"subpd	%%xmm2,%%xmm3		/* t2 = t2-t3; %%xmm2 FREE */\n\t"\
		"movl	%[__o1],%%eax		\n\t"\
		"movl	%[__o2],%%ebx		\n\t"\
		"movl	%[__o3],%%ecx		\n\t"\
		"movaps	%%xmm0,%%xmm2		/* cpy t0 */\n\t"\
		"movaps	%%xmm5,%%xmm7		/* cpy t3 */\n\t"\
		"addpd	%%xmm1,%%xmm0		/*~t0 = t0+tt */\n\t"\
		"addpd	%%xmm6,%%xmm5		/*~t3 = t3+t6 */\n\t"\
		"addpd	%%xmm3,%%xmm1		/*~tt = tt+t2 */\n\t"\
		"addpd	%%xmm4,%%xmm6		/*      t6+t5 */\n\t"\
		"addpd	%%xmm2,%%xmm3		/*~t2 = t2+t0 */\n\t"\
		"addpd	%%xmm7,%%xmm4		/*~t5 = t5+t3 */\n\t"\
		"subpd	%%xmm1,%%xmm2		/*~t1 = t0-tt-t2; %%xmm1 FREE */\n\t"\
		"subpd	%%xmm6,%%xmm7		/*~t4 = t3-t6-t5; %%xmm6 FREE */\n\t"\
		"movl	%[__o4],%%edx		\n\t"\
		"movl	%[__o5],%%esi		\n\t"\
		"movl	%[__o6],%%edi		\n\t"/* %%xmm1;6 FREE */\
		"movaps	(%%eax),%%xmm1		/* t0r */\n\t"\
		"movaps	(%%edi),%%xmm6		/* t3r */\n\t"\
		"subpd	%%xmm5,%%xmm1		/* B1r =t0r-t3i */\n\t"\
		"subpd	%%xmm6,%%xmm0		/* B6i =t0i-t3r */\n\t"\
		"addpd	%%xmm5,%%xmm5		/*        2*t3i */\n\t"\
		"addpd	%%xmm6,%%xmm6		/*        2*t3r */\n\t"\
		"addpd	%%xmm1,%%xmm5		/* B6r =t0r+t3i */\n\t"\
		"addpd	%%xmm0,%%xmm6		/* B1i =t0i+t3r */\n\t"\
		"movaps	%%xmm1,    (%%eax)	\n\t"/* <-B1r */\
		"movaps	%%xmm0,0x10(%%edi)	\n\t"/* <-B6i */\
		"movaps	%%xmm5,    (%%edi)	\n\t"/* <-B6r */\
		"movaps	%%xmm6,0x10(%%eax)	\n\t"/* <-B1i */\
		"movaps	(%%ebx),%%xmm1		/* t1r */\n\t"\
		"movaps	(%%esi),%%xmm6		/* t4r */\n\t"\
		"subpd	%%xmm7,%%xmm1		/* B2r =t1r-t4i */\n\t"\
		"subpd	%%xmm6,%%xmm2		/* B5i =t1i-t4r */\n\t"\
		"addpd	%%xmm7,%%xmm7		/*        2*t4i */\n\t"\
		"addpd	%%xmm6,%%xmm6		/*        2*t4r */\n\t"\
		"addpd	%%xmm1,%%xmm7		/* B5r =t1r+t4i */\n\t"\
		"addpd	%%xmm2,%%xmm6		/* B2i =t1i+t4r */\n\t"\
		"movaps	%%xmm1,    (%%ebx)	\n\t"/* <-B2r */\
		"movaps	%%xmm2,0x10(%%esi)	\n\t"/* <-B5i */\
		"movaps	%%xmm7,    (%%esi)	\n\t"/* <-B5r */\
		"movaps	%%xmm6,0x10(%%ebx)	\n\t"/* <-B2i */\
	/* Note the order reversal onthis pair of outputs: */\
		"movaps	(%%ecx),%%xmm0		/* t2r */\n\t"\
		"movaps	(%%edx),%%xmm5		/* t5r */\n\t"\
		"subpd	%%xmm4,%%xmm0		/* B4r =t2r-t5i */\n\t"\
		"subpd	%%xmm5,%%xmm3		/* B3i =t2i-t5r */\n\t"\
		"addpd	%%xmm4,%%xmm4		/*        2*t5i */\n\t"\
		"addpd	%%xmm5,%%xmm5		/*        2*t5r */\n\t"\
		"addpd	%%xmm0,%%xmm4		/* B3r =t2r+t5i */\n\t"\
		"addpd	%%xmm3,%%xmm5		/* B4i =t2i+t5r */\n\t"\
		"movaps	%%xmm0,    (%%edx)	\n\t"/* <-B4r */\
		"movaps	%%xmm3,0x10(%%ecx)	\n\t"/* <-B3i */\
		"movaps	%%xmm4,    (%%ecx)	\n\t"/* <-B3r */\
		"movaps	%%xmm5,0x10(%%edx)	\n\t"/* <-B4i */\
	"popl %%ebx	\n\t"\
		:					/* outputs: none */\
		: [__i0] "m" (Xi0)	/* All inputs from memory addresses here */\
		 ,[__i1] "m" (Xi1)\
		 ,[__i2] "m" (Xi2)\
		 ,[__i3] "m" (Xi3)\
		 ,[__i4] "m" (Xi4)\
		 ,[__i5] "m" (Xi5)\
		 ,[__i6] "m" (Xi6)\
		 ,[__cc] "m" (Xcc)\
		 ,[__o0] "m" (Xo0)\
		 ,[__o1] "m" (Xo1)\
		 ,[__o2] "m" (Xo2)\
		 ,[__o3] "m" (Xo3)\
		 ,[__o4] "m" (Xo4)\
		 ,[__o5] "m" (Xo5)\
		 ,[__o6] "m" (Xo6)\
		: "cc","memory","eax","ecx","edx","edi","esi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"	/* Clobbered registers */\
	);\
	}

	/* Twiddleless version of SSE2_RADIX8_DIF_TWIDDLE. Inputs enter in memory locations __r0 + [__i1,__i2,__i3,__i4,__i5,__i6,__i7],;
	where r0 is a memory address and the i's are LITERAL [BYTE] OFFSETS. Outputs go into memory locations __o0,__o1,__o2,__o3,__o4,__o5,__o6,__o7, assumed disjoint with inputs:\
	*/
	#define SSE2_RADIX8_DIF_0TWIDDLE(Xr0, Xi1,Xi2,Xi3,Xi4,Xi5,Xi6,Xi7, Xo0,Xo1,Xo2,Xo3,Xo4,Xo5,Xo6,Xo7, Xisrt2)\
	{\
	__asm__ volatile (\
	"pushl %%ebx	\n\t"/* Explicit save/restore of PIC register */\
		"movl	%[__r0],%%eax	/* i0 = r00 */	\n\t"\
		"leal	%c[__i2](%%eax),%%ebx			\n\t"\
		"leal	%c[__i4](%%eax),%%ecx			\n\t"\
		"leal	%c[__i6](%%eax),%%edx			\n\t"\
		"/* Do the p0,p4 combo: */		\n\t"\
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	    (%%ecx),%%xmm4	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"movaps	0x10(%%ecx),%%xmm5	\n\t"\
		"movaps	%%xmm0,%%xmm2			\n\t"\
		"movaps	%%xmm1,%%xmm3			\n\t"\
		"addpd	%%xmm4,%%xmm0			\n\t"\
		"addpd	%%xmm5,%%xmm1			\n\t"\
		"subpd	%%xmm4,%%xmm2			\n\t"\
		"subpd	%%xmm5,%%xmm3			\n\t"\
		"/* Do the p2,6 combo: */		\n\t"\
		"movaps	    (%%ebx),%%xmm4	\n\t"\
		"movaps	0x10(%%ebx),%%xmm5	\n\t"\
		"movaps	%%xmm4,%%xmm6			\n\t"\
		"movaps	%%xmm5,%%xmm7			\n\t"\
		"addpd	    (%%edx),%%xmm4	\n\t"\
		"addpd	0x10(%%edx),%%xmm5	\n\t"\
		"subpd	    (%%edx),%%xmm6	\n\t"\
		"subpd	0x10(%%edx),%%xmm7	\n\t"\
		"subpd	%%xmm4,%%xmm0			\n\t"\
		"subpd	%%xmm7,%%xmm2			\n\t"\
		"subpd	%%xmm5,%%xmm1			\n\t"\
		"subpd	%%xmm6,%%xmm3			\n\t"\
		"movaps	%%xmm0,    (%%ecx)	\n\t"\
		"movaps	%%xmm2,    (%%ebx)	\n\t"\
		"movaps	%%xmm1,0x10(%%ecx)	\n\t"\
		"movaps	%%xmm3,0x10(%%edx)	\n\t"\
		"addpd	%%xmm4,%%xmm4			\n\t"\
		"addpd	%%xmm7,%%xmm7			\n\t"\
		"addpd	%%xmm5,%%xmm5			\n\t"\
		"addpd	%%xmm6,%%xmm6			\n\t"\
		"addpd	%%xmm0,%%xmm4			\n\t"\
		"addpd	%%xmm2,%%xmm7			\n\t"\
		"addpd	%%xmm1,%%xmm5			\n\t"\
		"addpd	%%xmm3,%%xmm6			\n\t"\
		"movaps	%%xmm4,    (%%eax)	\n\t"\
		"movaps	%%xmm7,    (%%edx)	\n\t"\
		"movaps	%%xmm5,0x10(%%eax)	\n\t"\
		"movaps	%%xmm6,0x10(%%ebx)	\n\t"\
	/* Inline of SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO_B(r8,c1):\
	Do 2nd of 2 radix-4 subtransforms, but now keep outputs in registers: */\
		"leal	%c[__i3](%%eax),%%ebx			\n\t"\
		"leal	%c[__i5](%%eax),%%ecx			\n\t"\
		"leal	%c[__i7](%%eax),%%edx			\n\t"\
		"leal	%c[__i1](%%eax),%%eax	\n\t"/* Must do this last here */\
		"/* Do the p1,p5 combo: */		\n\t"\
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	    (%%ecx),%%xmm4	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"movaps	0x10(%%ecx),%%xmm5	\n\t"\
		"movaps	%%xmm0,%%xmm2			\n\t"\
		"movaps	%%xmm1,%%xmm3			\n\t"\
		"addpd	%%xmm4,%%xmm0			\n\t"\
		"addpd	%%xmm5,%%xmm1			\n\t"\
		"subpd	%%xmm4,%%xmm2			\n\t"\
		"subpd	%%xmm5,%%xmm3			\n\t"\
		"/* Do the p3,7 combo: */		\n\t"\
		"movaps	    (%%ebx),%%xmm4	\n\t"\
		"movaps	0x10(%%ebx),%%xmm5	\n\t"\
		"movaps	%%xmm4,%%xmm6			\n\t"\
		"movaps	%%xmm5,%%xmm7			\n\t"\
		"subpd	    (%%edx),%%xmm4	\n\t"\
		"subpd	0x10(%%edx),%%xmm5	\n\t"\
		"addpd	    (%%edx),%%xmm6	\n\t"\
		"addpd	0x10(%%edx),%%xmm7	\n\t"\
	/* Finish radix-4 butterfly and store just the 1st of the 8 outputs into output-array slots: */\
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
		"movaps	%%xmm6,    (%%eax)	\n\t"\
		"movaps	%%xmm7,0x10(%%eax)	\n\t"\
		"movaps	%%xmm2,%%xmm6			\n\t"\
		"movaps	%%xmm5,%%xmm7			\n\t"\
		"subpd	%%xmm4,%%xmm2			\n\t"\
		"subpd	%%xmm3,%%xmm5			\n\t"\
		"addpd	%%xmm6,%%xmm4			\n\t"\
		"addpd	%%xmm7,%%xmm3			\n\t"\
		"movl	%[__isrt2],%%eax		\n\t"\
		"movaps	(%%eax),%%xmm6	/* isrt2 */\n\t"\
		"mulpd	%%xmm6,%%xmm2			\n\t"\
		"mulpd	%%xmm6,%%xmm5			\n\t"\
		"mulpd	%%xmm6,%%xmm4			\n\t"\
		"mulpd	%%xmm6,%%xmm3			\n\t"\
	/* Inline the contents of SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_A(r0): Combine the 2 radix-4 subtransforms and write outputs to temp array: */\
	/* t0,1,2,3,4,5,6,7 in __i0,2,4,6 (eax,ebx,ecx,edx) */\
	/* t8,9,a,b,c,d,e,f in __i1,xmm2,4,0,1,5,3; xmm6,7 FREE */\
		"movl	%[__r0],%%edi			\n\t"\
		"leal	%c[__i2](%%edi),%%esi	\n\t"\
		"movl	%[__o4],%%eax	\n\t"\
		"movl	%[__o5],%%ebx	\n\t"\
		"movl	%[__o6],%%ecx	\n\t"\
		"movl	%[__o7],%%edx	\n\t"\
		"movaps	    (%%esi),%%xmm6		\n\t"\
		"movaps	0x10(%%esi),%%xmm7		\n\t"\
		"subpd   %%xmm2,%%xmm6				\n\t"\
		"subpd   %%xmm4,%%xmm7				\n\t"\
		"addpd   %%xmm2,%%xmm2				\n\t"\
		"addpd   %%xmm4,%%xmm4				\n\t"\
		"addpd   %%xmm6,%%xmm2				\n\t"\
		"addpd   %%xmm7,%%xmm4				\n\t"\
		"movaps	%%xmm6,    (%%ebx)	/* o5r */\n\t"\
		"movaps	%%xmm7,0x10(%%ebx)	/* o5i */\n\t"\
		"movaps	%%xmm2,    (%%eax)	/* o4r */\n\t"\
		"movaps	%%xmm4,0x10(%%eax)	/* o4i */\n\t"\
	/* Can't do 2nd set side-by-side with 1st, since only 2 registers (xmm6,7) available for temp storage: */\
		"leal	%c[__i6](%%edi),%%esi			\n\t"\
		"movaps	    (%%esi),%%xmm6		\n\t"\
		"movaps	0x10(%%esi),%%xmm7		\n\t"\
		"subpd   %%xmm3,%%xmm6				\n\t"\
		"subpd   %%xmm5,%%xmm7				\n\t"\
		"addpd   %%xmm3,%%xmm3				\n\t"\
		"addpd   %%xmm5,%%xmm5				\n\t"\
		"addpd   %%xmm6,%%xmm3				\n\t"\
		"addpd   %%xmm7,%%xmm5				\n\t"\
		"movaps	%%xmm6,    (%%ecx)	/* o6r */\n\t"\
		"movaps	%%xmm7,0x10(%%edx)	/* o7i */\n\t"\
		"movaps	%%xmm3,    (%%edx)	/* o7r */\n\t"\
		"movaps	%%xmm5,0x10(%%ecx)	/* o6i */\n\t"\
		"movl	%[__o0],%%eax	\n\t"\
		"movl	%[__o1],%%ebx	\n\t"\
		"movl	%[__o2],%%ecx	\n\t"\
		"movl	%[__o3],%%edx	\n\t"\
		"leal	%c[__i4](%%edi),%%esi	\n\t"\
		"movaps	    (%%edi),%%xmm6		\n\t"\
		"movaps	    (%%esi),%%xmm4		\n\t"\
		"movaps	0x10(%%edi),%%xmm7		\n\t"\
		"movaps	0x10(%%esi),%%xmm5		\n\t"\
		"leal	%c[__i1](%%edi),%%esi	\n\t"\
		"movaps	    (%%esi),%%xmm2		\n\t"\
		"movaps	0x10(%%esi),%%xmm3		\n\t"\
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
		"movaps	%%xmm6,    (%%ebx)	/* o1r */\n\t"\
		"movaps	%%xmm4,    (%%ecx)	/* o2r */\n\t"\
		"movaps	%%xmm7,0x10(%%ebx)	/* o1i */\n\t"\
		"movaps	%%xmm5,0x10(%%edx)	/* o3i */\n\t"\
		"movaps	%%xmm2,    (%%eax)	/* o0r */\n\t"\
		"movaps	%%xmm1,    (%%edx)	/* o3r */\n\t"\
		"movaps	%%xmm3,0x10(%%eax)	/* o0r */\n\t"\
		"movaps	%%xmm0,0x10(%%ecx)	/* o2i */\n\t"\
	"popl %%ebx	\n\t"\
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
		: "cc","memory","eax",/*"ebx",*/"ecx","edx","edi","esi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"	/* Clobbered registers */\
	);\
	}

	// Need a 2nd version of above which takes the i-strides as intvars rather than literal bytes:
	#define SSE2_RADIX8_DIF_0TWIDDLE_B(Xr0, Xi1,Xi2,Xi3,Xi4,Xi5,Xi6,Xi7, Xo0,Xo1,Xo2,Xo3,Xo4,Xo5,Xo6,Xo7, Xisrt2)\
	{\
	__asm__ volatile (\
	"pushl %%ebx	\n\t"/* Explicit save/restore of PIC register */\
		"movl	%[__r0],%%eax	/* i0 = r00 */	\n\t"\
		"movl	%[__i2],%%ebx	/* i2 */		\n\t"\
		"movl	%[__i4],%%ecx	/* i4 */		\n\t"\
		"movl	%[__i6],%%edx	/* i6 */		\n\t"\
		"addl	%%eax,%%ebx						\n\t"\
		"addl	%%eax,%%ecx						\n\t"\
		"addl	%%eax,%%edx						\n\t"\
		/* Do the p0,p4 combo: */\
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	    (%%ecx),%%xmm4	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"movaps	0x10(%%ecx),%%xmm5	\n\t"\
		"movaps	%%xmm0,%%xmm2			\n\t"\
		"movaps	%%xmm1,%%xmm3			\n\t"\
		"addpd	%%xmm4,%%xmm0			\n\t"\
		"addpd	%%xmm5,%%xmm1			\n\t"\
		"subpd	%%xmm4,%%xmm2			\n\t"\
		"subpd	%%xmm5,%%xmm3			\n\t"\
		/* Do the p2,6 combo: */\
		"movaps	    (%%ebx),%%xmm4	\n\t"\
		"movaps	0x10(%%ebx),%%xmm5	\n\t"\
		"movaps	%%xmm4,%%xmm6			\n\t"\
		"movaps	%%xmm5,%%xmm7			\n\t"\
		"addpd	    (%%edx),%%xmm4	\n\t"\
		"addpd	0x10(%%edx),%%xmm5	\n\t"\
		"subpd	    (%%edx),%%xmm6	\n\t"\
		"subpd	0x10(%%edx),%%xmm7	\n\t"\
		"subpd	%%xmm4,%%xmm0			\n\t"\
		"subpd	%%xmm7,%%xmm2			\n\t"\
		"subpd	%%xmm5,%%xmm1			\n\t"\
		"subpd	%%xmm6,%%xmm3			\n\t"\
		"movaps	%%xmm0,    (%%ecx)	\n\t"\
		"movaps	%%xmm2,    (%%ebx)	\n\t"\
		"movaps	%%xmm1,0x10(%%ecx)	\n\t"\
		"movaps	%%xmm3,0x10(%%edx)	\n\t"\
		"addpd	%%xmm4,%%xmm4			\n\t"\
		"addpd	%%xmm7,%%xmm7			\n\t"\
		"addpd	%%xmm5,%%xmm5			\n\t"\
		"addpd	%%xmm6,%%xmm6			\n\t"\
		"addpd	%%xmm0,%%xmm4			\n\t"\
		"addpd	%%xmm2,%%xmm7			\n\t"\
		"addpd	%%xmm1,%%xmm5			\n\t"\
		"addpd	%%xmm3,%%xmm6			\n\t"\
		"movaps	%%xmm4,    (%%eax)	\n\t"\
		"movaps	%%xmm7,    (%%edx)	\n\t"\
		"movaps	%%xmm5,0x10(%%eax)	\n\t"\
		"movaps	%%xmm6,0x10(%%ebx)	\n\t"\
	/* Inline of SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO_B(r8,c1):\
	Do 2nd of 2 radix-4 subtransforms, but now keep outputs in registers: */\
		"movl	%[__i3],%%ebx	/* i2 */		\n\t"\
		"movl	%[__i5],%%ecx	/* i4 */		\n\t"\
		"movl	%[__i7],%%edx	/* i6 */		\n\t"\
		"addl	%%eax,%%ebx						\n\t"\
		"addl	%%eax,%%ecx						\n\t"\
		"addl	%%eax,%%edx						\n\t"\
		"addl	%[__i1],%%eax	\n\t"/* Must do this last here */\
		"/* Do the p1,p5 combo: */		\n\t"\
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	    (%%ecx),%%xmm4	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"movaps	0x10(%%ecx),%%xmm5	\n\t"\
		"movaps	%%xmm0,%%xmm2			\n\t"\
		"movaps	%%xmm1,%%xmm3			\n\t"\
		"addpd	%%xmm4,%%xmm0			\n\t"\
		"addpd	%%xmm5,%%xmm1			\n\t"\
		"subpd	%%xmm4,%%xmm2			\n\t"\
		"subpd	%%xmm5,%%xmm3			\n\t"\
		"/* Do the p3,7 combo: */		\n\t"\
		"movaps	    (%%ebx),%%xmm4	\n\t"\
		"movaps	0x10(%%ebx),%%xmm5	\n\t"\
		"movaps	%%xmm4,%%xmm6			\n\t"\
		"movaps	%%xmm5,%%xmm7			\n\t"\
		"subpd	    (%%edx),%%xmm4	\n\t"\
		"subpd	0x10(%%edx),%%xmm5	\n\t"\
		"addpd	    (%%edx),%%xmm6	\n\t"\
		"addpd	0x10(%%edx),%%xmm7	\n\t"\
	/* Finish radix-4 butterfly and store just the 1st of the 8 outputs into output-array slots: */\
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
		"movaps	%%xmm6,    (%%eax)	\n\t"\
		"movaps	%%xmm7,0x10(%%eax)	\n\t"\
		"movaps	%%xmm2,%%xmm6			\n\t"\
		"movaps	%%xmm5,%%xmm7			\n\t"\
		"subpd	%%xmm4,%%xmm2			\n\t"\
		"subpd	%%xmm3,%%xmm5			\n\t"\
		"addpd	%%xmm6,%%xmm4			\n\t"\
		"addpd	%%xmm7,%%xmm3			\n\t"\
		"movl	%[__isrt2],%%eax		\n\t"\
		"movaps	(%%eax),%%xmm6	/* isrt2 */\n\t"\
		"mulpd	%%xmm6,%%xmm2			\n\t"\
		"mulpd	%%xmm6,%%xmm5			\n\t"\
		"mulpd	%%xmm6,%%xmm4			\n\t"\
		"mulpd	%%xmm6,%%xmm3			\n\t"\
	/* Inline the contents of SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_A(r0): Combine the 2 radix-4 subtransforms and write outputs to temp array: */\
	/* t0,1,2,3,4,5,6,7 in __i0,2,4,6 (eax,ebx,ecx,edx) */\
	/* t8,9,a,b,c,d,e,f in __i1,xmm2,4,0,1,5,3; xmm6,7 FREE */\
		"movl	%[__r0],%%edi			\n\t"\
		"movl	%[__i2],%%esi	/* i2 */	\n\t"\
		"addl	%%edi,%%esi					\n\t"\
		"movl	%[__o4],%%eax	\n\t"\
		"movl	%[__o5],%%ebx	\n\t"\
		"movl	%[__o6],%%ecx	\n\t"\
		"movl	%[__o7],%%edx	\n\t"\
		"movaps	    (%%esi),%%xmm6		\n\t"\
		"movaps	0x10(%%esi),%%xmm7		\n\t"\
		"subpd   %%xmm2,%%xmm6				\n\t"\
		"subpd   %%xmm4,%%xmm7				\n\t"\
		"addpd   %%xmm2,%%xmm2				\n\t"\
		"addpd   %%xmm4,%%xmm4				\n\t"\
		"addpd   %%xmm6,%%xmm2				\n\t"\
		"addpd   %%xmm7,%%xmm4				\n\t"\
		"movaps	%%xmm6,    (%%ebx)	/* o5r */\n\t"\
		"movaps	%%xmm7,0x10(%%ebx)	/* o5i */\n\t"\
		"movaps	%%xmm2,    (%%eax)	/* o4r */\n\t"\
		"movaps	%%xmm4,0x10(%%eax)	/* o4i */\n\t"\
	/* Can't do 2nd set side-by-side with 1st, since only 2 registers (xmm6,7) available for temp storage: */\
		"movl	%[__i6],%%esi	/* i6 */	\n\t"\
		"addl	%%edi,%%esi					\n\t"\
		"movaps	    (%%esi),%%xmm6		\n\t"\
		"movaps	0x10(%%esi),%%xmm7		\n\t"\
		"subpd   %%xmm3,%%xmm6				\n\t"\
		"subpd   %%xmm5,%%xmm7				\n\t"\
		"addpd   %%xmm3,%%xmm3				\n\t"\
		"addpd   %%xmm5,%%xmm5				\n\t"\
		"addpd   %%xmm6,%%xmm3				\n\t"\
		"addpd   %%xmm7,%%xmm5				\n\t"\
		"movaps	%%xmm6,    (%%ecx)	/* o6r */\n\t"\
		"movaps	%%xmm7,0x10(%%edx)	/* o7i */\n\t"\
		"movaps	%%xmm3,    (%%edx)	/* o7r */\n\t"\
		"movaps	%%xmm5,0x10(%%ecx)	/* o6i */\n\t"\
		"movl	%[__o0],%%eax	\n\t"\
		"movl	%[__o1],%%ebx	\n\t"\
		"movl	%[__o2],%%ecx	\n\t"\
		"movl	%[__o3],%%edx	\n\t"\
		"movl	%[__i4],%%esi	/* i4 */	\n\t"\
		"addl	%%edi,%%esi					\n\t"\
		"movaps	    (%%edi),%%xmm6		\n\t"\
		"movaps	    (%%esi),%%xmm4		\n\t"\
		"movaps	0x10(%%edi),%%xmm7		\n\t"\
		"movaps	0x10(%%esi),%%xmm5		\n\t"\
		"movl	%[__i1],%%esi	/* i1 */	\n\t"\
		"addl	%%edi,%%esi					\n\t"\
		"movaps	    (%%esi),%%xmm2		\n\t"\
		"movaps	0x10(%%esi),%%xmm3		\n\t"\
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
		"movaps	%%xmm6,    (%%ebx)	/* o1r */\n\t"\
		"movaps	%%xmm4,    (%%ecx)	/* o2r */\n\t"\
		"movaps	%%xmm7,0x10(%%ebx)	/* o1i */\n\t"\
		"movaps	%%xmm5,0x10(%%edx)	/* o3i */\n\t"\
		"movaps	%%xmm2,    (%%eax)	/* o0r */\n\t"\
		"movaps	%%xmm1,    (%%edx)	/* o3r */\n\t"\
		"movaps	%%xmm3,0x10(%%eax)	/* o0r */\n\t"\
		"movaps	%%xmm0,0x10(%%ecx)	/* o2i */\n\t"\
	"popl %%ebx	\n\t"\
		:					/* outputs: none */\
		: [__r0] "m" (Xr0)	/* All inputs from memory addresses here */\
		 ,[__i1] "m" (Xi1)\
		 ,[__i2] "m" (Xi2)\
		 ,[__i3] "m" (Xi3)\
		 ,[__i4] "m" (Xi4)\
		 ,[__i5] "m" (Xi5)\
		 ,[__i6] "m" (Xi6)\
		 ,[__i7] "m" (Xi7)\
		 ,[__o0] "m" (Xo0)\
		 ,[__o1] "m" (Xo1)\
		 ,[__o2] "m" (Xo2)\
		 ,[__o3] "m" (Xo3)\
		 ,[__o4] "m" (Xo4)\
		 ,[__o5] "m" (Xo5)\
		 ,[__o6] "m" (Xo6)\
		 ,[__o7] "m" (Xo7)\
		 ,[__isrt2] "m" (Xisrt2)\
		: "cc","memory","eax",/*"ebx",*/"ecx","edx","edi","esi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"	/* Clobbered registers */\
	);\
	}

	/* Twiddleless version of SSE2_RADIX8_DIT_TWIDDLE. Inputs enter in memory locations __i0,__i1,__i2,__i3,__i4,__i5,__i6,__i7.
	Outputs go into 16 contiguous 32-byte memory locations starting at __out and assumed disjoint with inputs.
	This macro built on the same code template as SSE2_RADIX8_DIF_TWIDDLE0, but with the I/O-location indices mutually bit reversed:
	01234567 <--> 04261537, which can be effected via the pairwise swaps 1 <--> 4 and 3 <--> 6.
	*/
	#define	SSE2_RADIX8_DIT_0TWIDDLE(Xin0,Xin1,Xin2,Xin3,Xin4,Xin5,Xin6,Xin7, Xout, Xisrt2)\
	{\
	__asm__ volatile (\
	"pushl %%ebx	\n\t"/* Explicit save/restore of PIC register */\
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
		"/* Move outputs t5r,i into a(j1,j2+p5), first doing the addsub and mul by ISRT2: */\n\t"\
		"/* Move outputs t7r,i into a(j1,j2+p7), first doing the addsub and mul by ISRT2: */\n\t"\
		"addpd	%%xmm2,%%xmm6				\n\t"\
		"addpd	%%xmm3,%%xmm7				\n\t"\
		"subpd	0xc0(%%esi),%%xmm2			\n\t"\
		"subpd	0xd0(%%esi),%%xmm3			\n\t"\
		"/* Move t6r,i into a(j1,j2+p6) */	\n\t"\
		"movaps	%%xmm2,0xc0(%%esi)			\n\t"\
		"movaps	%%xmm3,0xd0(%%esi)			\n\t"\
		"movaps	%%xmm4,%%xmm2				\n\t"\
		"movaps	%%xmm5,%%xmm3				\n\t"\
		"addpd	%%xmm0,%%xmm5				\n\t"\
		"subpd	%%xmm3,%%xmm0				\n\t"\
		"addpd	%%xmm1,%%xmm4				\n\t"\
		"subpd	%%xmm2,%%xmm1				\n\t"\
		"movaps	%%xmm5,%%xmm2				\n\t"\
		"movaps	%%xmm1,%%xmm3				\n\t"\
		"addpd	%%xmm1,%%xmm5				\n\t"\
		"movl	%[__isrt2],%%edi			\n\t"\
		"movaps	(%%edi),%%xmm1	/* isrt2 */	\n\t"\
		"subpd	%%xmm3,%%xmm2				\n\t"\
		"mulpd	%%xmm1,%%xmm5				\n\t"\
		"mulpd	%%xmm1,%%xmm2				\n\t"\
		"movaps	%%xmm0,%%xmm3				\n\t"\
		"movaps	%%xmm5,0xa0(%%esi)			\n\t"\
		"movaps	%%xmm4,%%xmm5				\n\t"\
		"addpd	%%xmm4,%%xmm0				\n\t"\
		"movaps	%%xmm2,0xb0(%%esi)			\n\t"\
		"subpd	%%xmm5,%%xmm3				\n\t"\
		"mulpd	%%xmm1,%%xmm0				\n\t"\
		"mulpd	%%xmm1,%%xmm3				\n\t"\
		"movaps	%%xmm0,0xe0(%%esi)			\n\t"\
		"movaps	%%xmm3,0xf0(%%esi)			\n\t"\
	/**************** 1st of the 2 length-4 subtransforms... **************/\
		"movl	%[__in0],%%eax	\n\t"\
		"movl	%[__in1],%%ebx	\n\t"\
		"movl	%[__in2],%%ecx	\n\t"\
		"movl	%[__in3],%%edx	\n\t"\
		"movaps	    (%%eax),%%xmm0			\n\t"\
		"movaps	0x10(%%eax),%%xmm1			\n\t"\
		"movaps	%%xmm0,%%xmm2				\n\t"\
		"movaps	%%xmm1,%%xmm3				\n\t"\
		"addpd	    (%%ebx),%%xmm2			\n\t"\
		"addpd	0x10(%%ebx),%%xmm3			\n\t"\
		"subpd	    (%%ebx),%%xmm0			\n\t"\
		"subpd	0x10(%%ebx),%%xmm1			\n\t"\
		"/* Move   t4r,i into temp(j1+p0) in anticipation of final outputs (t0+t4)r,i which will go there: */\n\t"\
		"movaps	%%xmm6,    (%%esi)			\n\t"\
		"movaps	%%xmm7,0x10(%%esi)			\n\t"\
		"addpd	%%xmm6,%%xmm6				\n\t"\
		"addpd	%%xmm7,%%xmm7				\n\t"\
		"/* Move 2*t4r,i into temp(j1+p4) in anticipation of final outputs (t0-t4)r,i which will go there: */\n\t"\
		"movaps	%%xmm6,0x80(%%esi)			\n\t"\
		"movaps	%%xmm7,0x90(%%esi)			\n\t"\
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
	/* Compute and dump first 2 outputs now, in order to free up 2 registers: */\
		"addpd	    (%%esi),%%xmm6			\n\t"\
		"addpd	0x10(%%esi),%%xmm7			\n\t"\
		"movaps	%%xmm6,    (%%esi)			\n\t"/* o0r */\
		"movaps	%%xmm7,0x10(%%esi)			\n\t"/* o0i */\
		"subpd	0x80(%%esi),%%xmm6			\n\t"\
		"subpd	0x90(%%esi),%%xmm7			\n\t"\
		"movaps	%%xmm6,0x80(%%esi)			\n\t"/* o4r */\
		"movaps	%%xmm7,0x90(%%esi)			\n\t"/* o4i */\
		"movaps	%%xmm4,%%xmm6				\n\t"\
		"movaps	%%xmm5,%%xmm7				\n\t"\
		"addpd	%%xmm0,%%xmm5				\n\t"\
		"subpd	%%xmm7,%%xmm0				\n\t"\
		"addpd	%%xmm1,%%xmm4				\n\t"\
		"subpd	%%xmm6,%%xmm1				\n\t"\
		"movaps	%%xmm2,%%xmm6				\n\t"\
		"movaps	%%xmm3,%%xmm7				\n\t"\
		"addpd	0xd0(%%esi),%%xmm2			\n\t"\
		"subpd	0xc0(%%esi),%%xmm3			\n\t"\
		"subpd	0xd0(%%esi),%%xmm6			\n\t"\
		"addpd	0xc0(%%esi),%%xmm7			\n\t"\
		"movaps	%%xmm2,0xc0(%%esi)			\n\t"/* o6r */\
		"movaps	%%xmm3,0xd0(%%esi)			\n\t"/* o6r */\
		"movaps	%%xmm6,0x40(%%esi)			\n\t"/* o2r */\
		"movaps	%%xmm7,0x50(%%esi)			\n\t"/* o2i */\
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
		"movaps	%%xmm5,0xe0(%%esi)			\n\t"/* o7r */\
		"movaps	%%xmm0,0xa0(%%esi)			\n\t"/* o5r */\
		"movaps	%%xmm1,0xf0(%%esi)			\n\t"/* o7i */\
		"movaps	%%xmm4,0xb0(%%esi)			\n\t"/* o5i */\
		"movaps	%%xmm2,0x60(%%esi)			\n\t"/* o3r */\
		"movaps	%%xmm6,0x20(%%esi)			\n\t"/* o1r */\
		"movaps	%%xmm3,0x70(%%esi)			\n\t"/* o3i */\
		"movaps	%%xmm7,0x30(%%esi)			\n\t"/* o1i */\
	"popl %%ebx	\n\t"\
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
		: "cc","memory","eax",/*"ebx",*/"ecx","edx","edi","esi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"	/* Clobbered registers */\
	);\
	}

	// Same as SSE2_RADIX8_DIT_0TWIDDLE but with user-specifiable [i.e. not nec. contiguous] output addresses:
	#define	SSE2_RADIX8_DIT_0TWIDDLE_OOP(Xi0,Xi1,Xi2,Xi3,Xi4,Xi5,Xi6,Xi7, Xo0,Xo1,Xo2,Xo3,Xo4,Xo5,Xo6,Xo7, Xisrt2)\
	{\
	__asm__ volatile (\
	"pushl %%ebx	\n\t"/* Explicit save/restore of PIC register */\
		"movl	%[__i4],%%eax	\n\t"\
		"movl	%[__i5],%%ebx	\n\t"\
		"movl	%[__i6],%%ecx	\n\t"\
		"movl	%[__i7],%%edx	\n\t"\
		"movaps	    (%%eax),%%xmm0			\n\t"\
		"movaps	0x10(%%eax),%%xmm1			\n\t"\
		"movaps	%%xmm0,%%xmm2				\n\t"\
		"movaps	%%xmm1,%%xmm3				\n\t"\
		"addpd	    (%%ebx),%%xmm2			\n\t"\
		"addpd	0x10(%%ebx),%%xmm3			\n\t"\
		"subpd	    (%%ebx),%%xmm0			\n\t"\
		"subpd	0x10(%%ebx),%%xmm1			\n\t"\
		"movaps	    (%%ecx),%%xmm4			\n\t"\
		"movaps	0x10(%%ecx),%%xmm5			\n\t"\
		"movaps	%%xmm4,%%xmm6				\n\t"\
		"movaps	%%xmm5,%%xmm7				\n\t"\
		"addpd	    (%%edx),%%xmm6			\n\t"\
		"addpd	0x10(%%edx),%%xmm7			\n\t"\
		"subpd	    (%%edx),%%xmm4			\n\t"\
		"subpd	0x10(%%edx),%%xmm5			\n\t"\
		"/* Copy t6r,i into main-array slot add6 */\n\t"\
		"movl	%[__o6],%%esi	\n\t"\
		"/*movl	%[__o7],%%edi	*/\n\t"\
		"movaps	%%xmm6,    (%%esi)			\n\t"\
		"movaps	%%xmm7,0x10(%%esi)			\n\t"\
		"/* Copy t7r,i into main-array slot add7 */\n\t"\
		"/* movaps	%%xmm4,    (%%edi)		wtf? these 2 stores are useless, since we rewrite */	\n\t"\
		"/* movaps	%%xmm5,0x10(%%edi)		same memlocs below, sans an intervening read-from */	\n\t"\
		"/* Move outputs t5r,i into a(j1,j2+p5), first doing the addsub and mul by ISRT2: */\n\t"\
		"/* Move outputs t7r,i into a(j1,j2+p7), first doing the addsub and mul by ISRT2: */\n\t"\
		"addpd	%%xmm2,%%xmm6				\n\t"\
		"addpd	%%xmm3,%%xmm7				\n\t"\
		"subpd	    (%%esi),%%xmm2			\n\t"\
		"subpd	0x10(%%esi),%%xmm3			\n\t"\
		"/* Move t6r,i into a(j1,j2+p6) */	\n\t"\
		"movaps	%%xmm2,    (%%esi)			\n\t"\
		"movaps	%%xmm3,0x10(%%esi)			\n\t"\
		"movaps	%%xmm4,%%xmm2				\n\t"\
		"movaps	%%xmm5,%%xmm3				\n\t"\
		"addpd	%%xmm0,%%xmm5				\n\t"\
		"subpd	%%xmm3,%%xmm0				\n\t"\
		"addpd	%%xmm1,%%xmm4				\n\t"\
		"subpd	%%xmm2,%%xmm1				\n\t"\
		"movaps	%%xmm5,%%xmm2				\n\t"\
		"movaps	%%xmm1,%%xmm3				\n\t"\
		"addpd	%%xmm1,%%xmm5				\n\t"\
		"movl	%[__isrt2],%%edi			\n\t"\
		"movaps	(%%edi),%%xmm1	/* isrt2 */	\n\t"\
		"subpd	%%xmm3,%%xmm2				\n\t"\
		"mulpd	%%xmm1,%%xmm5				\n\t"\
		"mulpd	%%xmm1,%%xmm2				\n\t"\
		"movaps	%%xmm0,%%xmm3				\n\t"\
		"movl	%[__o5],%%esi	\n\t"\
		"movl	%[__o7],%%edi	\n\t"\
		"movaps	%%xmm5,    (%%esi)			\n\t"\
		"movaps	%%xmm4,%%xmm5				\n\t"\
		"addpd	%%xmm4,%%xmm0				\n\t"\
		"movaps	%%xmm2,0x10(%%esi)			\n\t"\
		"subpd	%%xmm5,%%xmm3				\n\t"\
		"mulpd	%%xmm1,%%xmm0				\n\t"\
		"mulpd	%%xmm1,%%xmm3				\n\t"\
		"movaps	%%xmm0,    (%%edi)			\n\t"\
		"movaps	%%xmm3,0x10(%%edi)			\n\t"\
	/**************** 1st of the 2 length-4 subtransforms... **************/\
		"movl	%[__i0],%%eax	\n\t"\
		"movl	%[__i1],%%ebx	\n\t"\
		"movl	%[__i2],%%ecx	\n\t"\
		"movl	%[__i3],%%edx	\n\t"\
		"movaps	    (%%eax),%%xmm0			\n\t"\
		"movaps	0x10(%%eax),%%xmm1			\n\t"\
		"movaps	%%xmm0,%%xmm2				\n\t"\
		"movaps	%%xmm1,%%xmm3				\n\t"\
		"addpd	    (%%ebx),%%xmm2			\n\t"\
		"addpd	0x10(%%ebx),%%xmm3			\n\t"\
		"subpd	    (%%ebx),%%xmm0			\n\t"\
		"subpd	0x10(%%ebx),%%xmm1			\n\t"\
		"/* Move   t4r,i into temp(j1+p0) in anticipation of final outputs (t0+t4)r,i which will go there: */\n\t"\
		"movl	%[__o0],%%esi	\n\t"\
		"movl	%[__o4],%%edi	\n\t"\
		"movaps	%%xmm6,    (%%esi)			\n\t"\
		"movaps	%%xmm7,0x10(%%esi)			\n\t"\
		"addpd	%%xmm6,%%xmm6				\n\t"\
		"addpd	%%xmm7,%%xmm7				\n\t"\
		"/* Move 2*t4r,i into temp(j1+p4) in anticipation of final outputs (t0-t4)r,i which will go there: */\n\t"\
		"movaps	%%xmm6,    (%%edi)			\n\t"\
		"movaps	%%xmm7,0x10(%%edi)			\n\t"\
		"movaps	    (%%ecx),%%xmm4			\n\t"\
		"movaps	0x10(%%ecx),%%xmm5			\n\t"\
		"movaps	%%xmm4,%%xmm6				\n\t"\
		"movaps	%%xmm5,%%xmm7				\n\t"\
		"addpd	    (%%edx),%%xmm6			\n\t"\
		"addpd	0x10(%%edx),%%xmm7			\n\t"\
		"subpd	    (%%edx),%%xmm4			\n\t"\
		"subpd	0x10(%%edx),%%xmm5			\n\t"\
		"/* Copy t5,6 into temp(j1+p2) */	\n\t"\
		"movl	%[__o2],%%esi	\n\t"\
		"movaps	%%xmm6,    (%%esi)			\n\t"\
		"movaps	%%xmm7,0x10(%%esi)			\n\t"\
		"addpd	%%xmm2,%%xmm6				\n\t"\
		"addpd	%%xmm3,%%xmm7				\n\t"\
		"subpd	    (%%esi),%%xmm2			\n\t"\
		"subpd	0x10(%%esi),%%xmm3			\n\t"\
	/* Compute and dump first 2 outputs now, in order to free up 2 registers: */\
		"movl	%[__o0],%%esi	\n\t"\
		"addpd	    (%%esi),%%xmm6			\n\t"\
		"addpd	0x10(%%esi),%%xmm7			\n\t"\
		"movaps	%%xmm6,    (%%esi)			\n\t"\
		"movaps	%%xmm7,0x10(%%esi)			\n\t"\
		"subpd	    (%%edi),%%xmm6			\n\t"\
		"subpd	0x10(%%edi),%%xmm7			\n\t"\
		"movaps	%%xmm6,    (%%edi)			\n\t"\
		"movaps	%%xmm7,0x10(%%edi)			\n\t"\
		"movaps	%%xmm4,%%xmm6				\n\t"\
		"movaps	%%xmm5,%%xmm7				\n\t"\
		"addpd	%%xmm0,%%xmm5				\n\t"\
		"subpd	%%xmm7,%%xmm0				\n\t"\
		"addpd	%%xmm1,%%xmm4				\n\t"\
		"subpd	%%xmm6,%%xmm1				\n\t"\
		"movaps	%%xmm2,%%xmm6				\n\t"\
		"movaps	%%xmm3,%%xmm7				\n\t"\
		"movl	%[__o2],%%esi	\n\t"\
		"movl	%[__o6],%%edi	\n\t"\
		"addpd	0x10(%%edi),%%xmm2			\n\t"\
		"subpd	    (%%edi),%%xmm3			\n\t"\
		"subpd	0x10(%%edi),%%xmm6			\n\t"\
		"addpd	    (%%edi),%%xmm7			\n\t"\
		"movaps	%%xmm2,    (%%edi)			\n\t"\
		"movaps	%%xmm3,0x10(%%edi)			\n\t"\
		"movaps	%%xmm6,    (%%esi)			\n\t"\
		"movaps	%%xmm7,0x10(%%esi)			\n\t"\
		"movl	%[__o1],%%eax	\n\t"\
		"movl	%[__o3],%%ebx	\n\t"\
		"movl	%[__o5],%%ecx	\n\t"\
		"movl	%[__o7],%%edx	\n\t"\
		"movaps	%%xmm5,%%xmm2				\n\t"\
		"movaps	%%xmm0,%%xmm6				\n\t"\
		"movaps	%%xmm1,%%xmm3				\n\t"\
		"movaps	%%xmm4,%%xmm7				\n\t"\
		"addpd	    (%%ecx),%%xmm5			\n\t"\
		"subpd	0x10(%%edx),%%xmm0			\n\t"\
		"subpd	0x10(%%ecx),%%xmm1			\n\t"\
		"subpd	    (%%edx),%%xmm4			\n\t"\
		"subpd	    (%%ecx),%%xmm2			\n\t"\
		"addpd	0x10(%%edx),%%xmm6			\n\t"\
		"addpd	0x10(%%ecx),%%xmm3			\n\t"\
		"addpd	    (%%edx),%%xmm7			\n\t"\
		"movaps	%%xmm5,    (%%edx)			\n\t"\
		"movaps	%%xmm0,    (%%ecx)			\n\t"\
		"movaps	%%xmm1,0x10(%%edx)			\n\t"\
		"movaps	%%xmm4,0x10(%%ecx)			\n\t"\
		"movaps	%%xmm2,    (%%ebx)			\n\t"\
		"movaps	%%xmm6,    (%%eax)			\n\t"\
		"movaps	%%xmm3,0x10(%%ebx)			\n\t"\
		"movaps	%%xmm7,0x10(%%eax)			\n\t"\
	"popl %%ebx	\n\t"\
		:					/* outputs: none */\
		: [__i0] "m" (Xi0)	/* All inputs from memory addresses here */\
		 ,[__i1] "m" (Xi1)\
		 ,[__i2] "m" (Xi2)\
		 ,[__i3] "m" (Xi3)\
		 ,[__i4] "m" (Xi4)\
		 ,[__i5] "m" (Xi5)\
		 ,[__i6] "m" (Xi6)\
		 ,[__i7] "m" (Xi7)\
		 ,[__o0] "m" (Xo0)\
		 ,[__o1] "m" (Xo1)\
		 ,[__o2] "m" (Xo2)\
		 ,[__o3] "m" (Xo3)\
		 ,[__o4] "m" (Xo4)\
		 ,[__o5] "m" (Xo5)\
		 ,[__o6] "m" (Xo6)\
		 ,[__o7] "m" (Xo7)\
		 ,[__isrt2] "m" (Xisrt2)\
		: "cc","memory","eax",/*"ebx",*/"ecx","edx","edi","esi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"	/* Clobbered registers */\
	);\
	}

	#define SSE2_RADIX8_DIF_TWIDDLE_OOP(Xi0,Xi1,Xi2,Xi3,Xi4,Xi5,Xi6,Xi7, Xo0,Xo1,Xo2,Xo3,Xo4,Xo5,Xo6,Xo7 ,Xc1,Xs1,Xc2,Xs2,Xc3,Xs3,Xc4,Xs4,Xc5,Xs5,Xc6,Xs6,Xc7,Xs7)\
	{\
	__asm__ volatile (\
	"pushl %%ebx	\n\t"/* Explicit save/restore of PIC register */\
		/* Block 0,4: */\
		"movl		%[i0]	,%%eax			\n\t"\
		"movl		%[i4]	,%%ebx			\n\t"\
		"movl		%[c4]	,%%ecx			\n\t"\
		"movl		%[s4]	,%%esi			\n\t"\
	/* [esi] (and if needed rdi) points to sine components of each sincos pair, which is not really a pair here in terms of relative addressing: */\
		"movaps		    (%%ebx)	,%%xmm2		\n\t"\
		"movaps		0x10(%%ebx)	,%%xmm3		\n\t"\
		"movaps		(%%ecx)	,%%xmm0		\n\t"\
		"movaps		(%%esi)	,%%xmm1		\n\t"\
		"movaps		%%xmm2	,%%xmm4		\n\t"\
		"movaps		%%xmm3	,%%xmm5		\n\t"\
		"mulpd		%%xmm0	,%%xmm2		\n\t"\
		"mulpd		%%xmm0	,%%xmm3		\n\t"\
		"mulpd		%%xmm1	,%%xmm4		\n\t"\
		"mulpd		%%xmm1	,%%xmm5		\n\t"\
		"movaps		    (%%eax)	,%%xmm0		\n\t"\
		"movaps		0x10(%%eax)	,%%xmm1		\n\t"\
		"subpd		%%xmm5		,%%xmm2		\n\t"\
		"addpd		%%xmm4		,%%xmm3		\n\t"\
		"movaps		%%xmm0	,%%xmm6		\n\t"\
		"movaps		%%xmm1	,%%xmm7		\n\t"\
		"addpd		%%xmm2		,%%xmm0		\n\t"\
		"addpd		%%xmm3		,%%xmm1		\n\t"\
		"subpd		%%xmm2		,%%xmm6		\n\t"\
		"subpd		%%xmm3		,%%xmm7		\n\t"\
		"movaps		%%xmm0		,    (%%eax)	\n\t"/* __tr0 */\
		"movaps		%%xmm1		,0x10(%%eax)	\n\t"/* __ti0 */\
		"movaps		%%xmm6		,    (%%ebx)	\n\t"/* __tr1 */\
		"movaps		%%xmm7		,0x10(%%ebx)	\n\t"/* __ti1 */\
		/* Block 2,6: */\
		"movl		%[i2]	,%%eax			\n\t"\
		"movl		%[i6]	,%%ebx			\n\t"\
		"movl		%[c2]	,%%ecx			\n\t"\
		"movl		%[c6]	,%%edx			\n\t"\
		"movl		%[s2]	,%%esi			\n\t"\
		"movl		%[s6]	,%%edi			\n\t"\
		"movaps		    (%%eax)	,%%xmm0		\n\t"\
		"movaps		0x10(%%eax)	,%%xmm2		\n\t"\
		"movaps		(%%ecx)	,%%xmm6		\n\t"\
		"movaps		(%%esi)	,%%xmm7		\n\t"\
		"movaps		%%xmm0	,%%xmm1		\n\t"\
		"movaps		%%xmm2	,%%xmm3		\n\t"\
		"mulpd		%%xmm6	,%%xmm0		\n\t"\
		"mulpd		%%xmm7	,%%xmm2		\n\t"\
		"mulpd		%%xmm7	,%%xmm1		\n\t"\
		"mulpd		%%xmm6	,%%xmm3		\n\t"\
		"subpd		%%xmm2		,%%xmm0		\n\t"\
		"addpd		%%xmm3		,%%xmm1		\n\t"\
		"movaps		    (%%ebx)	,%%xmm2		\n\t"\
		"movaps		0x10(%%ebx)	,%%xmm3		\n\t"\
		"movaps		(%%edx)	,%%xmm6		\n\t"\
		"movaps		(%%edi)	,%%xmm7		\n\t"\
		"movaps		%%xmm2	,%%xmm4		\n\t"\
		"movaps		%%xmm3	,%%xmm5		\n\t"\
		"mulpd		%%xmm6	,%%xmm2		\n\t"\
		"mulpd		%%xmm7	,%%xmm3		\n\t"\
		"mulpd		%%xmm7	,%%xmm4		\n\t"\
		"mulpd		%%xmm6	,%%xmm5		\n\t"\
		"subpd		%%xmm3		,%%xmm2		\n\t"\
		"addpd		%%xmm5		,%%xmm4		\n\t"\
		"movaps		%%xmm2		,%%xmm3		\n\t"\
		"movaps		%%xmm4		,%%xmm5		\n\t"\
		"addpd		%%xmm0		,%%xmm2		\n\t"\
		"subpd		%%xmm3		,%%xmm0		\n\t"\
		"addpd		%%xmm1		,%%xmm4		\n\t"\
		"subpd		%%xmm5		,%%xmm1		\n\t"\
		"movaps		%%xmm2		,    (%%eax)	\n\t"/* __tr2 */\
		"movaps		%%xmm4		,0x10(%%eax)	\n\t"/* __ti2 */\
		"movaps		%%xmm0		,    (%%ebx)	\n\t"/* __tr3 */\
		"movaps		%%xmm1		,0x10(%%ebx)	\n\t"/* __ti3 */\
		/* Block 1,5: */\
		"movl		%[i1]	,%%eax			\n\t"\
		"movl		%[i5]	,%%ebx			\n\t"\
		"movl		%[c1]	,%%ecx			\n\t"\
		"movl		%[c5]	,%%edx			\n\t"\
		"movl		%[s1]	,%%esi			\n\t"\
		"movl		%[s5]	,%%edi			\n\t"\
		"movaps		    (%%eax)	,%%xmm0		\n\t"\
		"movaps		0x10(%%eax)	,%%xmm2		\n\t"\
		"movaps		(%%ecx)	,%%xmm6		\n\t"\
		"movaps		(%%esi)	,%%xmm7		\n\t"\
		"movaps		%%xmm0	,%%xmm1		\n\t"\
		"movaps		%%xmm2	,%%xmm3		\n\t"\
		"mulpd		%%xmm6	,%%xmm0		\n\t"\
		"mulpd		%%xmm7	,%%xmm2		\n\t"\
		"mulpd		%%xmm7	,%%xmm1		\n\t"\
		"mulpd		%%xmm6	,%%xmm3		\n\t"\
		"subpd		%%xmm2		,%%xmm0		\n\t"\
		"addpd		%%xmm3		,%%xmm1		\n\t"\
		"movaps		    (%%ebx)	,%%xmm2		\n\t"\
		"movaps		0x10(%%ebx)	,%%xmm3		\n\t"\
		"movaps		(%%edx)	,%%xmm6		\n\t"\
		"movaps		(%%edi)	,%%xmm7		\n\t"\
		"movaps		%%xmm2	,%%xmm4		\n\t"\
		"movaps		%%xmm3	,%%xmm5		\n\t"\
		"mulpd		%%xmm6	,%%xmm2		\n\t"\
		"mulpd		%%xmm7	,%%xmm3		\n\t"\
		"mulpd		%%xmm7	,%%xmm4		\n\t"\
		"mulpd		%%xmm6	,%%xmm5		\n\t"\
		"subpd		%%xmm3		,%%xmm2		\n\t"\
		"addpd		%%xmm5		,%%xmm4		\n\t"\
		"movaps		%%xmm2		,%%xmm3		\n\t"\
		"movaps		%%xmm4		,%%xmm5		\n\t"\
		"addpd		%%xmm0		,%%xmm2		\n\t"\
		"subpd		%%xmm3		,%%xmm0		\n\t"\
		"addpd		%%xmm1		,%%xmm4		\n\t"\
		"subpd		%%xmm5		,%%xmm1		\n\t"\
		"movaps		%%xmm2		,    (%%eax)	\n\t"/* __tr4 */\
		"movaps		%%xmm4		,0x10(%%eax)	\n\t"/* __ti4 */\
		"movaps		%%xmm0		,    (%%ebx)	\n\t"/* __tr5 */\
		"movaps		%%xmm1		,0x10(%%ebx)	\n\t"/* __ti5 */\
		/* Block 3,7: */\
		"movl		%[i3]	,%%eax			\n\t"\
		"movl		%[i7]	,%%ebx			\n\t"\
		"movl		%[c3]	,%%ecx			\n\t"\
		"movl		%[c7]	,%%edx			\n\t"\
		"movl		%[s3]	,%%esi			\n\t"\
		"movl		%[s7]	,%%edi			\n\t"\
		"movaps		    (%%eax)	,%%xmm0		\n\t"\
		"movaps		0x10(%%eax)	,%%xmm2		\n\t"\
		"movaps		(%%ecx)	,%%xmm6		\n\t"\
		"movaps		(%%esi)	,%%xmm7		\n\t"\
		"movaps		%%xmm0	,%%xmm1		\n\t"\
		"movaps		%%xmm2	,%%xmm3		\n\t"\
		"mulpd		%%xmm6	,%%xmm0		\n\t"\
		"mulpd		%%xmm7	,%%xmm2		\n\t"\
		"mulpd		%%xmm7	,%%xmm1		\n\t"\
		"mulpd		%%xmm6	,%%xmm3		\n\t"\
		"subpd		%%xmm2		,%%xmm0		\n\t"\
		"addpd		%%xmm3		,%%xmm1		\n\t"\
		"movaps		    (%%ebx)	,%%xmm2		\n\t"\
		"movaps		0x10(%%ebx)	,%%xmm3		\n\t"\
		"movaps		(%%edx)	,%%xmm6		\n\t"\
		"movaps		(%%edi)	,%%xmm7		\n\t"\
		"movaps		%%xmm2	,%%xmm4		\n\t"\
		"movaps		%%xmm3	,%%xmm5		\n\t"\
		"mulpd		%%xmm6	,%%xmm2		\n\t"\
		"mulpd		%%xmm7	,%%xmm3		\n\t"\
		"mulpd		%%xmm7	,%%xmm4		\n\t"\
		"mulpd		%%xmm6	,%%xmm5		\n\t"\
		"subpd		%%xmm3		,%%xmm2		\n\t"\
		"addpd		%%xmm5		,%%xmm4		\n\t"\
		"movaps		%%xmm2		,%%xmm3		\n\t"\
		"movaps		%%xmm4		,%%xmm5		\n\t"\
		"addpd		%%xmm0		,%%xmm2		\n\t"\
		"subpd		%%xmm3		,%%xmm0		\n\t"\
		"addpd		%%xmm1		,%%xmm4		\n\t"\
		"subpd		%%xmm5		,%%xmm1		\n\t"\
		"movaps		%%xmm2		,    (%%eax)	\n\t"/* __tr6 */\
		"movaps		%%xmm4		,0x10(%%eax)	\n\t"/* __ti6 */\
		"movaps		%%xmm0		,    (%%ebx)	\n\t"/* __tr7 */\
		"movaps		%%xmm1		,0x10(%%ebx)	\n\t"/* __ti7 */\
	/***********************************************/\
	/* Combine to get the 2 length-4 transforms... */\
	/***********************************************/\
		/* Combo 1: */\
		"movl		%[i0]	,%%eax			\n\t"\
		"movl		%[i2]	,%%ebx			\n\t"\
		"movaps		    (%%eax)	,%%xmm0		\n\t"\
		"movaps		0x10(%%eax)	,%%xmm1		\n\t"\
		"movaps		    (%%ebx)	,%%xmm2		\n\t"/* __rt = __tr2 */\
		"movaps		0x10(%%ebx)	,%%xmm3		\n\t"/* __it = __ti2 */\
		"movaps		%%xmm0		,%%xmm4		\n\t"/* copy __tr0; */\
		"movaps		%%xmm1		,%%xmm5		\n\t"/* copy __ti0; */\
		"addpd		%%xmm2	,%%xmm0			\n\t"/* __tr0 = __tr0 + __rt; */\
		"subpd		%%xmm2	,%%xmm4			\n\t"/* __tr2 = __tr0 - __rt; */\
		"addpd		%%xmm3	,%%xmm1			\n\t"/* __ti0 = __ti0 + __it; */\
		"subpd		%%xmm3	,%%xmm5			\n\t"/* __ti2 = __ti0 - __it; */\
		"movaps		%%xmm0	,    (%%eax)	\n\t"\
		"movaps		%%xmm1	,0x10(%%eax)	\n\t"\
		"movaps		%%xmm4	,    (%%ebx)	\n\t"\
		"movaps		%%xmm5	,0x10(%%ebx)	\n\t"\
		"movl		%[i1]	,%%ecx			\n\t"\
		"movl		%[i3]	,%%edx			\n\t"\
		"movaps		    (%%ecx)	,%%xmm2		\n\t"\
		"movaps		0x10(%%ecx)	,%%xmm3		\n\t"\
		"movaps		%%xmm2		,%%xmm6		\n\t"/* copy __tr4; [__rt = __tr6] */\
		"movaps		%%xmm3		,%%xmm7		\n\t"/* copy __ti4; [__it = __ti6] */\
		"addpd		    (%%edx)	,%%xmm2		\n\t"/* __tr4 = __tr4 + __rt; */\
		"subpd		    (%%edx)	,%%xmm6		\n\t"/* __tr6 = __tr4 - __rt; */\
		"addpd		0x10(%%edx)	,%%xmm3		\n\t"\
		"subpd		0x10(%%edx)	,%%xmm7		\n\t"\
		"subpd		%%xmm2		,%%xmm0		\n\t"/* __Br1 = __tr0 - __tr4; */\
		"subpd		%%xmm3		,%%xmm1		\n\t"/* __Bi1 = __ti0 - __ti4; */\
		"subpd		%%xmm7		,%%xmm4		\n\t"/* __Br2 = __tr2 - __ti6; */\
		"subpd		%%xmm6		,%%xmm5		\n\t"/* __Bi3 = __ti2 - __tr6; */\
		"addpd			(%%eax)	,%%xmm2		\n\t"/* __Br0 = __tr0 + __tr4; */\
		"addpd		0x10(%%eax)	,%%xmm3		\n\t"/* __Bi0 = __ti0 + __ti4; */\
		"addpd			(%%ebx)	,%%xmm7		\n\t"/* __Br3 = __tr2 + __ti6; */\
		"addpd		0x10(%%ebx)	,%%xmm6		\n\t"/* __Bi2 = __ti2 + __tr6; */\
		"movl		%[o0]	,%%eax			\n\t"\
		"movl		%[o2]	,%%ebx			\n\t"\
		"movl		%[o1]	,%%ecx			\n\t"\
		"movl		%[o3]	,%%edx			\n\t"\
		"movaps		%%xmm2		,    (%%eax)	\n\t"/* [o0].re */\
		"movaps		%%xmm3		,0x10(%%eax)	\n\t"/* [o0].im */\
		"movaps		%%xmm4		,    (%%ebx)	\n\t"/* [o2].re */\
		"movaps		%%xmm6		,0x10(%%ebx)	\n\t"/* [o2].im */\
		"movaps		%%xmm0		,    (%%ecx)	\n\t"/* [o1].re */\
		"movaps		%%xmm1		,0x10(%%ecx)	\n\t"/* [o1].im */\
		"movaps		%%xmm7		,    (%%edx)	\n\t"/* [o3].re */\
		"movaps		%%xmm5		,0x10(%%edx)	\n\t"/* [o3].im */\
		/* Combo 2: */\
		"movl		%[i4]	,%%eax			\n\t"\
		"movl		%[i6]	,%%ebx			\n\t"\
		"movaps		    (%%eax)	,%%xmm0		\n\t"\
		"movaps		0x10(%%eax)	,%%xmm1		\n\t"\
		"movaps		    (%%ebx)	,%%xmm2		\n\t"/* __rt = __tr3 */\
		"movaps		0x10(%%ebx)	,%%xmm3		\n\t"/* __it = __ti3 */\
		"movaps		%%xmm0		,%%xmm4		\n\t"/* copy __tr1; */\
		"movaps		%%xmm1		,%%xmm5		\n\t"/* copy __ti1; */\
		"subpd		%%xmm3	,%%xmm0			\n\t"/* __tr1 = __tr1 - __it; */\
		"addpd		%%xmm3	,%%xmm4			\n\t"/* __tr3 = __tr1 + __it; */\
		"addpd		%%xmm2	,%%xmm1			\n\t"/* __ti1 = __ti1 + __rt; */\
		"subpd		%%xmm2	,%%xmm5			\n\t"/* __ti3 = __ti1 - __rt; */\
		"movaps		%%xmm0	,    (%%eax)	\n\t"/* __tr1 */\
		"movaps		%%xmm1	,0x10(%%eax)	\n\t"/* __ti1 */\
		"movaps		%%xmm4	,    (%%ebx)	\n\t"/* __tr3 */\
		"movaps		%%xmm5	,0x10(%%ebx)	\n\t"/* __ti3 */\
		"movl		%[i5]	,%%ecx			\n\t"\
		"movl		%[i7]	,%%edx			\n\t"\
		"movaps		    (%%ecx)	,%%xmm2		\n\t"\
		"movaps		0x10(%%ecx)	,%%xmm3		\n\t"\
		"movaps		%%xmm2		,%%xmm6		\n\t"/* copy __tr5; [__rt = __tr7] */\
		"movaps		%%xmm3		,%%xmm7		\n\t"/* copy __ti5; [__it = __ti7] */\
		"subpd		0x10(%%edx)	,%%xmm2		\n\t"/* __tr5 = __tr5 - __rt; */\
		"addpd		0x10(%%edx)	,%%xmm6		\n\t"/* __tr7 = __tr5 + __rt; */\
		"addpd		    (%%edx)	,%%xmm3		\n\t"/* __ti5 = __ti5 + __it; */\
		"subpd		    (%%edx)	,%%xmm7		\n\t"/* __ti7 = __ti5 - __it; */\
	/* t:01,23,45,67 in [i]:04,26,15,37 */\
	/***********************************************/\
	/* Now combine the two half-transforms:        */\
	/***********************************************/\
	/* Use the cosine term of the [c1,s1] pair, which is the *middle* [4th of 7] of our 7 input pairs, in terms \
	of the input-arg bit-reversal reordering defined in the __X[c,s] --> [c,s] mapping below and happens to \
	always in fact *be* a true cosine term, which is a requirement for our "decr 1 gives isrt2" data-copy scheme: */\
		"movl		%[c1],%%esi			\n\t"/* isrt2 in [c1]-1 */\
		"movaps		%%xmm2		,%%xmm5		\n\t"/* copy __tr5 */\
		"subpd		%%xmm3		,%%xmm2		\n\t"/* __it  = __tr5 + __ti5; */\
		"addpd		%%xmm3		,%%xmm5		\n\t"/* __rt  = __tr5 - __ti5; */\
		"mulpd	-0x10(%%esi)	,%%xmm2		\n\t"/* __rt *= ISRT2; */\
		"mulpd	-0x10(%%esi)	,%%xmm5		\n\t"/* __it *= ISRT2; */\
		"movaps		0x10(%%ebx)	,%%xmm3		\n\t"/* __tr3 */\
		"movaps		%%xmm7		,%%xmm4		\n\t"/* copy __ti7 */\
		"addpd		%%xmm6		,%%xmm4		\n\t"/* __rt  = __tr7 + __ti7; */\
		"subpd		%%xmm6		,%%xmm7		\n\t"/* __it  = __ti7 - __tr7; */\
		"mulpd	-0x10(%%esi)	,%%xmm4		\n\t"/* __rt *= ISRT2; */\
		"mulpd	-0x10(%%esi)	,%%xmm7		\n\t"/* __it *= ISRT2; */\
		"movaps			(%%ebx)	,%%xmm6		\n\t"/* __ti3; last ref to edx-for-in-address */\
		"movl		%[o5]	,%%ecx			\n\t"\
		"movl		%[o7]	,%%edx			\n\t"\
		"subpd		%%xmm2		,%%xmm0		\n\t"/* __Br5 = __tr1 - __rt; */\
		"subpd		%%xmm5		,%%xmm1		\n\t"/* __Bi5 = __ti1 - __it; */\
		"subpd		%%xmm4		,%%xmm6		\n\t"/* __Br6 = __tr3 - __rt; */\
		"subpd		%%xmm7		,%%xmm3		\n\t"/* __Bi6 = __ti3 - __it; */\
		"addpd			(%%eax)	,%%xmm2		\n\t"/* __Br4 = __tr1 + __rt; */\
		"addpd		0x10(%%eax)	,%%xmm5		\n\t"/* __Bi4 = __ti1 + __it; */\
		"addpd			(%%ebx)	,%%xmm4		\n\t"/* __Br7 = __tr3 + __rt; */\
		"addpd		0x10(%%ebx)	,%%xmm7		\n\t"/* __Bi7 = __ti3 + __it; */\
		"movl		%[o4]	,%%eax			\n\t"\
		"movl		%[o6]	,%%ebx			\n\t"\
		"movaps		%%xmm2		,    (%%eax)	\n\t"/* [o4].re */\
		"movaps		%%xmm5		,0x10(%%eax)	\n\t"/* [o4].im */\
		"movaps		%%xmm6		,    (%%ebx)	\n\t"/* [o6].re */\
		"movaps		%%xmm3		,0x10(%%ebx)	\n\t"/* [o6].im */\
		"movaps		%%xmm0		,    (%%ecx)	\n\t"/* [o5].re */\
		"movaps		%%xmm1		,0x10(%%ecx)	\n\t"/* [o5].im */\
		"movaps		%%xmm4		,    (%%edx)	\n\t"/* [o7].re */\
		"movaps		%%xmm7		,0x10(%%edx)	\n\t"/* [o7].im */\
	"popl %%ebx	\n\t"\
		:					/* outputs: none */\
		: [i0] "m" (Xi0)	/* All inputs from memory addresses here */\
		 ,[i1] "m" (Xi1)\
		 ,[i2] "m" (Xi2)\
		 ,[i3] "m" (Xi3)\
		 ,[i4] "m" (Xi4)\
		 ,[i5] "m" (Xi5)\
		 ,[i6] "m" (Xi6)\
		 ,[i7] "m" (Xi7)\
		 ,[o0] "m" (Xo0)\
		 ,[o1] "m" (Xo1)\
		 ,[o2] "m" (Xo2)\
		 ,[o3] "m" (Xo3)\
		 ,[o4] "m" (Xo4)\
		 ,[o5] "m" (Xo5)\
		 ,[o6] "m" (Xo6)\
		 ,[o7] "m" (Xo7)\
		 ,[c4] "m" (Xc1),[s4] "m" (Xs1)\
		 ,[c2] "m" (Xc2),[s2] "m" (Xs2)\
		 ,[c6] "m" (Xc3),[s6] "m" (Xs3)\
		 ,[c1] "m" (Xc4),[s1] "m" (Xs4)\
		 ,[c5] "m" (Xc5),[s5] "m" (Xs5)\
		 ,[c3] "m" (Xc6),[s3] "m" (Xs6)\
		 ,[c7] "m" (Xc7),[s7] "m" (Xs7)\
		: "cc","memory","eax",/*"ebx",*/"ecx","edx","esi","edi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"	/* Clobbered registers */\
	);\
	}

	#define SSE2_RADIX8_DIT_TWIDDLE_OOP(Xi0,Xi1,Xi2,Xi3,Xi4,Xi5,Xi6,Xi7, Xo0,Xo1,Xo2,Xo3,Xo4,Xo5,Xo6,Xo7 ,Xc1,Xs1,Xc2,Xs2,Xc3,Xs3,Xc4,Xs4,Xc5,Xs5,Xc6,Xs6,Xc7,Xs7)\
	{\
	__asm__ volatile (\
	"pushl %%ebx	\n\t"/* Explicit save/restore of PIC register */\
		/* Block 0,1: */\
		"movl		%[i0]	,%%eax			\n\t"\
		"movl		%[i1]	,%%ebx			\n\t"\
		"movl		%[c1]	,%%ecx			\n\t"\
		"movl		%[s1]	,%%esi			\n\t"\
	/* [esi] (and if needed rdi) points to sine components of each sincos pair, which is not really a pair here in terms of relative addressing: */\
		"movaps		    (%%ebx)	,%%xmm2		\n\t"\
		"movaps		0x10(%%ebx)	,%%xmm3		\n\t"\
		"movaps		(%%ecx)	,%%xmm0		\n\t"\
		"movaps		(%%esi)	,%%xmm1		\n\t"\
		"movaps		%%xmm2	,%%xmm4		\n\t"\
		"movaps		%%xmm3	,%%xmm5		\n\t"\
		"mulpd		%%xmm0	,%%xmm2		\n\t"/* tr*c */\
		"mulpd		%%xmm0	,%%xmm3		\n\t"/* ti*c */\
		"mulpd		%%xmm1	,%%xmm4		\n\t"/* tr*s */\
		"mulpd		%%xmm1	,%%xmm5		\n\t"/* ti*s */\
		"movaps		    (%%eax)	,%%xmm0		\n\t"\
		"movaps		0x10(%%eax)	,%%xmm1		\n\t"\
		"addpd		%%xmm5		,%%xmm2		\n\t"\
		"subpd		%%xmm4		,%%xmm3		\n\t"/* ti*c - tr*s ... sign correc, no need to negate xmm3 below */\
		"movaps		%%xmm0	,%%xmm6		\n\t"\
		"movaps		%%xmm1	,%%xmm7		\n\t"\
		"addpd		%%xmm2		,%%xmm0		\n\t"\
		"addpd		%%xmm3		,%%xmm1		\n\t"\
		"subpd		%%xmm2		,%%xmm6		\n\t"\
		"subpd		%%xmm3		,%%xmm7		\n\t"\
		"movaps		%%xmm0		,    (%%eax)	\n\t"/* __tr0 */\
		"movaps		%%xmm1		,0x10(%%eax)	\n\t"/* __ti0 */\
		"movaps		%%xmm6		,    (%%ebx)	\n\t"/* __tr1 */\
		"movaps		%%xmm7		,0x10(%%ebx)	\n\t"/* __ti1 */\
		/* Block 2,3: */\
		"movl		%[i2]	,%%eax			\n\t"\
		"movl		%[i3]	,%%ebx			\n\t"\
		"movl		%[c2]	,%%ecx			\n\t"\
		"movl		%[c3]	,%%edx			\n\t"\
		"movl		%[s2]	,%%esi			\n\t"\
		"movl		%[s3]	,%%edi			\n\t"\
		"movaps		    (%%eax)	,%%xmm0		\n\t"\
		"movaps		0x10(%%eax)	,%%xmm2		\n\t"\
		"movaps		(%%ecx)	,%%xmm6		\n\t"\
		"movaps		(%%esi)	,%%xmm7		\n\t"\
		"movaps		%%xmm0	,%%xmm1		\n\t"\
		"movaps		%%xmm2	,%%xmm3		\n\t"\
		"mulpd		%%xmm6	,%%xmm0		\n\t"/* tr*c */\
		"mulpd		%%xmm7	,%%xmm2		\n\t"/* ti*s */\
		"mulpd		%%xmm7	,%%xmm1		\n\t"/* tr*s */\
		"mulpd		%%xmm6	,%%xmm3		\n\t"/* ti*c */\
		"addpd		%%xmm2		,%%xmm0		\n\t"/* tr*c + ti*s */\
		"subpd		%%xmm3		,%%xmm1		\n\t"/* tr*s - ti*c ... need to flip sign of xmm1... */\
		"movaps		    (%%ebx)	,%%xmm2		\n\t"\
		"movaps		0x10(%%ebx)	,%%xmm3		\n\t"\
		"movaps		(%%edx)	,%%xmm6		\n\t"\
		"movaps		(%%edi)	,%%xmm7		\n\t"\
		"movaps		%%xmm2	,%%xmm4		\n\t"\
		"movaps		%%xmm3	,%%xmm5		\n\t"\
		"mulpd		%%xmm6	,%%xmm2		\n\t"/* tr*c */\
		"mulpd		%%xmm7	,%%xmm3		\n\t"/* ti*s */\
		"mulpd		%%xmm7	,%%xmm4		\n\t"/* tr*s */\
		"mulpd		%%xmm6	,%%xmm5		\n\t"/* ti*c */\
		"addpd		%%xmm3		,%%xmm2		\n\t"/* tr*c + ti*s */\
		"subpd		%%xmm5		,%%xmm4		\n\t"/* tr*s - ti*c ... need to flip sign of xmm4... */\
		"movaps		%%xmm2		,%%xmm3		\n\t"\
		"movaps		%%xmm4		,%%xmm5		\n\t"/* copy xmm4, thus need to flip sign of xmm4,5 below */\
		"addpd		%%xmm0		,%%xmm2		\n\t"\
		"subpd		%%xmm3		,%%xmm0		\n\t"\
		"addpd		%%xmm1		,%%xmm4		\n\t"/* need to flip sign of xmm1,4, thus of the sum (==> xmm4)... */\
		"subpd		%%xmm5		,%%xmm1		\n\t"/* ditto for (xmm1-xmm5) ==> xmm1... */\
		"xorpd		%%xmm6		,%%xmm6		\n\t"/* 0.0 */\
		"movaps		%%xmm6		,%%xmm7		\n\t"/* 0.0 */\
		"subpd		%%xmm4		,%%xmm6		\n\t"/* flip sign of xmm4... */\
		"subpd		%%xmm1		,%%xmm7		\n\t"/* flip sign of xmm1... */\
		"movaps		%%xmm2		,    (%%eax)	\n\t"/* __tr2 */\
		"movaps		%%xmm6		,0x10(%%eax)	\n\t"/* __ti2 */\
		"movaps		%%xmm0		,    (%%ebx)	\n\t"/* __tr3 */\
		"movaps		%%xmm7		,0x10(%%ebx)	\n\t"/* __ti3 */\
		/* Block 4,5: */\
		"movl		%[i4]	,%%eax			\n\t"\
		"movl		%[i5]	,%%ebx			\n\t"\
		"movl		%[c4]	,%%ecx			\n\t"\
		"movl		%[c5]	,%%edx			\n\t"\
		"movl		%[s4]	,%%esi			\n\t"\
		"movl		%[s5]	,%%edi			\n\t"\
		"movaps		    (%%eax)	,%%xmm0		\n\t"\
		"movaps		0x10(%%eax)	,%%xmm2		\n\t"\
		"movaps		(%%ecx)	,%%xmm6		\n\t"\
		"movaps		(%%esi)	,%%xmm7		\n\t"\
		"movaps		%%xmm0	,%%xmm1		\n\t"\
		"movaps		%%xmm2	,%%xmm3		\n\t"\
		"mulpd		%%xmm6	,%%xmm0		\n\t"\
		"mulpd		%%xmm7	,%%xmm2		\n\t"\
		"mulpd		%%xmm7	,%%xmm1		\n\t"\
		"mulpd		%%xmm6	,%%xmm3		\n\t"\
		"addpd		%%xmm2		,%%xmm0		\n\t"\
		"subpd		%%xmm3		,%%xmm1		\n\t"\
		"movaps		    (%%ebx)	,%%xmm2		\n\t"\
		"movaps		0x10(%%ebx)	,%%xmm3		\n\t"\
		"movaps		(%%edx)	,%%xmm6		\n\t"\
		"movaps		(%%edi)	,%%xmm7		\n\t"\
		"movaps		%%xmm2	,%%xmm4		\n\t"\
		"movaps		%%xmm3	,%%xmm5		\n\t"\
		"mulpd		%%xmm6	,%%xmm2		\n\t"\
		"mulpd		%%xmm7	,%%xmm3		\n\t"\
		"mulpd		%%xmm7	,%%xmm4		\n\t"\
		"mulpd		%%xmm6	,%%xmm5		\n\t"\
		"addpd		%%xmm3		,%%xmm2		\n\t"\
		"subpd		%%xmm5		,%%xmm4		\n\t"\
		"movaps		%%xmm2		,%%xmm3		\n\t"\
		"movaps		%%xmm4		,%%xmm5		\n\t"\
		"addpd		%%xmm0		,%%xmm2		\n\t"\
		"subpd		%%xmm3		,%%xmm0		\n\t"\
		"addpd		%%xmm1		,%%xmm4		\n\t"\
		"subpd		%%xmm5		,%%xmm1		\n\t"\
		"xorpd		%%xmm6		,%%xmm6		\n\t"/* 0.0 */\
		"movaps		%%xmm6		,%%xmm7		\n\t"/* 0.0 */\
		"subpd		%%xmm4		,%%xmm6		\n\t"/* flip sign of xmm4... */\
		"subpd		%%xmm1		,%%xmm7		\n\t"/* flip sign of xmm1... */\
		"movaps		%%xmm2		,    (%%eax)	\n\t"/* __tr4 */\
		"movaps		%%xmm6		,0x10(%%eax)	\n\t"/* __ti4 */\
		"movaps		%%xmm0		,    (%%ebx)	\n\t"/* __tr5 */\
		"movaps		%%xmm7		,0x10(%%ebx)	\n\t"/* __ti5 */\
		/* Block 6,7: */\
		"movl		%[i6]	,%%eax			\n\t"\
		"movl		%[i7]	,%%ebx			\n\t"\
		"movl		%[c6]	,%%ecx			\n\t"\
		"movl		%[c7]	,%%edx			\n\t"\
		"movl		%[s6]	,%%esi			\n\t"\
		"movl		%[s7]	,%%edi			\n\t"\
		"movaps		    (%%eax)	,%%xmm0		\n\t"\
		"movaps		0x10(%%eax)	,%%xmm2		\n\t"\
		"movaps		(%%ecx)	,%%xmm6		\n\t"\
		"movaps		(%%esi)	,%%xmm7		\n\t"\
		"movaps		%%xmm0	,%%xmm1		\n\t"\
		"movaps		%%xmm2	,%%xmm3		\n\t"\
		"mulpd		%%xmm6	,%%xmm0		\n\t"\
		"mulpd		%%xmm7	,%%xmm2		\n\t"\
		"mulpd		%%xmm7	,%%xmm1		\n\t"\
		"mulpd		%%xmm6	,%%xmm3		\n\t"\
		"addpd		%%xmm2		,%%xmm0		\n\t"\
		"subpd		%%xmm3		,%%xmm1		\n\t"\
		"movaps		    (%%ebx)	,%%xmm2		\n\t"\
		"movaps		0x10(%%ebx)	,%%xmm3		\n\t"\
		"movaps		(%%edx)	,%%xmm6		\n\t"\
		"movaps		(%%edi)	,%%xmm7		\n\t"\
		"movaps		%%xmm2	,%%xmm4		\n\t"\
		"movaps		%%xmm3	,%%xmm5		\n\t"\
		"mulpd		%%xmm6	,%%xmm2		\n\t"\
		"mulpd		%%xmm7	,%%xmm3		\n\t"\
		"mulpd		%%xmm7	,%%xmm4		\n\t"\
		"mulpd		%%xmm6	,%%xmm5		\n\t"\
		"addpd		%%xmm3		,%%xmm2		\n\t"\
		"subpd		%%xmm5		,%%xmm4		\n\t"\
		"movaps		%%xmm2		,%%xmm3		\n\t"\
		"movaps		%%xmm4		,%%xmm5		\n\t"\
		"addpd		%%xmm0		,%%xmm2		\n\t"\
		"subpd		%%xmm3		,%%xmm0		\n\t"\
		"addpd		%%xmm1		,%%xmm4		\n\t"\
		"subpd		%%xmm5		,%%xmm1		\n\t"\
		"xorpd		%%xmm6		,%%xmm6		\n\t"/* 0.0 */\
		"movaps		%%xmm6		,%%xmm7		\n\t"/* 0.0 */\
		"subpd		%%xmm4		,%%xmm6		\n\t"/* flip sign of xmm4... */\
		"subpd		%%xmm1		,%%xmm7		\n\t"/* flip sign of xmm1... */\
		"movaps		%%xmm2		,    (%%eax)	\n\t"/* __tr6 */\
		"movaps		%%xmm6		,0x10(%%eax)	\n\t"/* __ti6 */\
		"movaps		%%xmm0		,    (%%ebx)	\n\t"/* __tr7 */\
		"movaps		%%xmm7		,0x10(%%ebx)	\n\t"/* __ti7 */\
	/***********************************************/\
	/* Combine to get the 2 length-4 transforms... */\
	/***********************************************/\
		/* Combo 1: */\
		"movl		%[i0]	,%%eax			\n\t"\
		"movl		%[i2]	,%%ebx			\n\t"\
		"movaps		    (%%eax)	,%%xmm0		\n\t"\
		"movaps		0x10(%%eax)	,%%xmm1		\n\t"\
		"movaps		    (%%ebx)	,%%xmm2		\n\t"/* __rt = __tr2 */\
		"movaps		0x10(%%ebx)	,%%xmm3		\n\t"/* __it = __ti2 */\
		"movaps		%%xmm0		,%%xmm4		\n\t"/* copy __tr0; */\
		"movaps		%%xmm1		,%%xmm5		\n\t"/* copy __ti0; */\
		"addpd		%%xmm2	,%%xmm0			\n\t"/* __tr0 = __tr0 + __rt; */\
		"subpd		%%xmm2	,%%xmm4			\n\t"/* __tr2 = __tr0 - __rt; */\
		"addpd		%%xmm3	,%%xmm1			\n\t"/* __ti0 = __ti0 + __it; */\
		"subpd		%%xmm3	,%%xmm5			\n\t"/* __ti2 = __ti0 - __it; */\
		"movaps		%%xmm0	,    (%%eax)	\n\t"\
		"movaps		%%xmm1	,0x10(%%eax)	\n\t"\
		"movaps		%%xmm4	,    (%%ebx)	\n\t"\
		"movaps		%%xmm5	,0x10(%%ebx)	\n\t"\
		"movl		%[i4]	,%%ecx			\n\t"\
		"movl		%[i6]	,%%edx			\n\t"\
		"movaps		    (%%ecx)	,%%xmm2		\n\t"\
		"movaps		0x10(%%ecx)	,%%xmm3		\n\t"\
		"movaps		%%xmm2		,%%xmm6		\n\t"/* copy __tr4; [__rt = __tr6] */\
		"movaps		%%xmm3		,%%xmm7		\n\t"/* copy __ti4; [__it = __ti6] */\
		"addpd		    (%%edx)	,%%xmm2		\n\t"/* __tr4 = __tr4 + __rt; */\
		"subpd		    (%%edx)	,%%xmm6		\n\t"/* __tr6 = __tr4 - __rt; */\
		"addpd		0x10(%%edx)	,%%xmm3		\n\t"/* __ti4 = __ti4 + __it; */\
		"subpd		0x10(%%edx)	,%%xmm7		\n\t"/* __ti6 = __ti4 - __it; */\
		"subpd		%%xmm2		,%%xmm0		\n\t"/* __Br1 = __tr0 - __tr4; */\
		"subpd		%%xmm3		,%%xmm1		\n\t"/* __Bi1 = __ti0 - __ti4; */\
		"subpd		%%xmm7		,%%xmm4		\n\t"/* __Br2 = __tr2 - __ti6; */\
		"subpd		%%xmm6		,%%xmm5		\n\t"/* __Bi3 = __ti2 - __tr6; */\
		"addpd			(%%eax)	,%%xmm2		\n\t"/* __Br0 = __tr0 + __tr4; */\
		"addpd		0x10(%%eax)	,%%xmm3		\n\t"/* __Bi0 = __ti0 + __ti4; */\
		"addpd			(%%ebx)	,%%xmm7		\n\t"/* __Br3 = __tr2 + __ti6; */\
		"addpd		0x10(%%ebx)	,%%xmm6		\n\t"/* __Bi2 = __ti2 + __tr6; */\
		"movl		%[o0]	,%%eax			\n\t"\
		"movl		%[o2]	,%%ebx			\n\t"/* Swap o2/3 here vs DIF, by swapping regs b/d */\
		"movl		%[o1]	,%%ecx			\n\t"\
		"movl		%[o3]	,%%edx			\n\t"\
		"movaps		%%xmm2		,    (%%eax)	\n\t"/* [o0].re */\
		"movaps		%%xmm3		,0x10(%%eax)	\n\t"/* [o0].im */\
		"movaps		%%xmm4		,    (%%edx)	\n\t"/* [o3].re */\
		"movaps		%%xmm6		,0x10(%%edx)	\n\t"/* [o3].im */\
		"movaps		%%xmm0		,    (%%ecx)	\n\t"/* [o1].re */\
		"movaps		%%xmm1		,0x10(%%ecx)	\n\t"/* [o1].im */\
		"movaps		%%xmm7		,    (%%ebx)	\n\t"/* [o2].re */\
		"movaps		%%xmm5		,0x10(%%ebx)	\n\t"/* [o2].im */\
		/* Combo 2: */\
		"movl		%[i1]	,%%eax			\n\t"\
		"movl		%[i3]	,%%ebx			\n\t"\
		"movaps		    (%%eax)	,%%xmm0		\n\t"\
		"movaps		0x10(%%eax)	,%%xmm1		\n\t"\
		"movaps		    (%%ebx)	,%%xmm2		\n\t"/* __rt = __tr3 */\
		"movaps		0x10(%%ebx)	,%%xmm3		\n\t"/* __it = __ti3 */\
		"movaps		%%xmm0		,%%xmm4		\n\t"/* copy __tr1; */\
		"movaps		%%xmm1		,%%xmm5		\n\t"/* copy __ti1; */\
		"addpd		%%xmm3	,%%xmm0			\n\t"/* __tr1 = __tr1 + __it; */\
		"subpd		%%xmm3	,%%xmm4			\n\t"/* __tr3 = __tr1 - __it; */\
		"subpd		%%xmm2	,%%xmm1			\n\t"/* __ti1 = __ti1 - __rt; */\
		"addpd		%%xmm2	,%%xmm5			\n\t"/* __ti3 = __ti1 + __rt; */\
		"movaps		%%xmm0	,    (%%eax)	\n\t"/* __tr1 */\
		"movaps		%%xmm1	,0x10(%%eax)	\n\t"/* __ti1 */\
		"movaps		%%xmm4	,    (%%ebx)	\n\t"/* __tr3 */\
		"movaps		%%xmm5	,0x10(%%ebx)	\n\t"/* __ti3 */\
		"movl		%[i5]	,%%ecx			\n\t"\
		"movl		%[i7]	,%%edx			\n\t"\
		"movaps		    (%%ecx)	,%%xmm2		\n\t"\
		"movaps		0x10(%%ecx)	,%%xmm3		\n\t"\
		"movaps		%%xmm2		,%%xmm6		\n\t"/* copy __tr5; [__rt = __tr7] */\
		"movaps		%%xmm3		,%%xmm7		\n\t"/* copy __ti5; [__it = __ti7] */\
		"addpd		0x10(%%edx)	,%%xmm2		\n\t"/* __tr5 = __tr5 + __rt; */\
		"subpd		0x10(%%edx)	,%%xmm6		\n\t"/* __tr7 = __tr5 - __rt; */\
		"subpd		    (%%edx)	,%%xmm3		\n\t"/* __ti5 = __ti5 - __it; */\
		"addpd		    (%%edx)	,%%xmm7		\n\t"/* __ti7 = __ti5 + __it; */\
	/* t:01,23,45,67 in [i]:01,23,45,67 */\
	/***********************************************/\
	/* Now combine the two half-transforms:        */\
	/***********************************************/\
	/* Use the cosine term of the [c4,s4] pair, which is the *middle* [4th of 7] of our 7 input pairs, in terms \
	of the input-arg bit-reversal reordering defined in the __X[c,s] --> [c,s] mapping below and happens to \
	always in fact *be* a true cosine term, which is a requirement for our "decr 1 gives isrt2" data-copy scheme: */\
		"movl		%[c4],%%esi			\n\t"/* isrt2 in [c4]-1 */\
		"movaps		%%xmm2		,%%xmm5		\n\t"/* copy __tr5; flip +- below to effect it/rt swap vs DIF */\
		"addpd		%%xmm3		,%%xmm2		\n\t"/* __rt5  = __tr5 + __ti5; */\
		"subpd		%%xmm3		,%%xmm5		\n\t"/* __it5  = __tr5 - __ti5; */\
		"mulpd	-0x10(%%esi)	,%%xmm2		\n\t"/* __rt5 *= ISRT2; */\
		"mulpd	-0x10(%%esi)	,%%xmm5		\n\t"/* __it5 *= ISRT2; */\
		"movaps		0x10(%%ebx)	,%%xmm3		\n\t"/* __ti3 */\
		"movaps		%%xmm6		,%%xmm4		\n\t"/* copy __tr7; flip +- and xmm6/7 below to effect it/rt swap vs DIF */\
		"subpd		%%xmm7		,%%xmm4		\n\t"/* __rt7  = __tr7 - __ti7; */\
		"addpd		%%xmm7		,%%xmm6		\n\t"/* __it7  = __tr7 + __ti7; */\
		"mulpd	-0x10(%%esi)	,%%xmm4		\n\t"/* __rt7 *= ISRT2; */\
		"mulpd	-0x10(%%esi)	,%%xmm6		\n\t"/* __it7 *= ISRT2; */\
		"movaps			(%%ebx)	,%%xmm7		\n\t"/* __tr3; last ref to edx-for-in-address */\
		"movl		%[o5]	,%%ecx			\n\t"\
		"movl		%[o7]	,%%edx			\n\t"/* __Bi4/5 swapped vs DIF: */\
		"subpd		%%xmm2		,%%xmm0		\n\t"/* __Br5 = __tr1 - __rt5; */\
		"subpd		%%xmm5		,%%xmm1		\n\t"/* __Bi4 = __ti1 - __it5; */\
		"subpd		%%xmm4		,%%xmm7		\n\t"/* __Br6 = __tr3 - __rt7; */\
		"subpd		%%xmm6		,%%xmm3		\n\t"/* __Bi6 = __ti3 - __it7; */\
		"addpd			(%%eax)	,%%xmm2		\n\t"/* __Br4 = __tr1 + __rt5; */\
		"addpd		0x10(%%eax)	,%%xmm5		\n\t"/* __Bi5 = __ti1 + __it5; */\
		"addpd			(%%ebx)	,%%xmm4		\n\t"/* __Br7 = __tr3 + __rt7; */\
		"addpd		0x10(%%ebx)	,%%xmm6		\n\t"/* __Bi7 = __ti3 + __it7; */\
		"movl		%[o4]	,%%eax			\n\t"\
		"movl		%[o6]	,%%ebx			\n\t"/* __Bi4/5 swapped vs DIF: */\
		"movaps		%%xmm2		,    (%%eax)	\n\t"/* [o4].re */\
		"movaps		%%xmm5		,0x10(%%ecx)	\n\t"/* [o5].im */\
		"movaps		%%xmm7		,    (%%ebx)	\n\t"/* [o6].re */\
		"movaps		%%xmm3		,0x10(%%ebx)	\n\t"/* [o6].im */\
		"movaps		%%xmm0		,    (%%ecx)	\n\t"/* [o5].re */\
		"movaps		%%xmm1		,0x10(%%eax)	\n\t"/* [o4].im */\
		"movaps		%%xmm4		,    (%%edx)	\n\t"/* [o7].re */\
		"movaps		%%xmm6		,0x10(%%edx)	\n\t"/* [o7].im */\
	"popl %%ebx	\n\t"\
		:					/* outputs: none */\
		: [i0] "m" (Xi0)	/* All inputs from memory addresses here */\
		 ,[i1] "m" (Xi1)\
		 ,[i2] "m" (Xi2)\
		 ,[i3] "m" (Xi3)\
		 ,[i4] "m" (Xi4)\
		 ,[i5] "m" (Xi5)\
		 ,[i6] "m" (Xi6)\
		 ,[i7] "m" (Xi7)\
		 ,[o0] "m" (Xo0)\
		 ,[o1] "m" (Xo1)\
		 ,[o2] "m" (Xo2)\
		 ,[o3] "m" (Xo3)\
		 ,[o4] "m" (Xo4)\
		 ,[o5] "m" (Xo5)\
		 ,[o6] "m" (Xo6)\
		 ,[o7] "m" (Xo7)\
		 ,[c1] "m" (Xc1),[s1] "m" (Xs1)\
		 ,[c2] "m" (Xc2),[s2] "m" (Xs2)\
		 ,[c3] "m" (Xc3),[s3] "m" (Xs3)\
		 ,[c4] "m" (Xc4),[s4] "m" (Xs4)\
		 ,[c5] "m" (Xc5),[s5] "m" (Xs5)\
		 ,[c6] "m" (Xc6),[s6] "m" (Xs6)\
		 ,[c7] "m" (Xc7),[s7] "m" (Xs7)\
		: "cc","memory","eax",/*"ebx",*/"ecx","edx","esi","edi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"	/* Clobbered registers */\
	);\
	}

	// Based on the SSE2_RADIX16_DIT_NOTWIDDLE macro in radix16_ditN_cy_dif1_gcc64.h, but with completely
	// specifiable 16-input addressing required for usage as the power-of-2 component of a twiddleless
	// radix = [odd*2^n] DFT routine.
	// We use just a single output base-pointer plus literal ostrides which are [1,2,3,4]-multiples of
	// __01; this allows us to cut GP-register usage, which is absolutely a must for the 32-bit version
	// of the macro, and is a benefit to the 64-bit versions which code-fold to yield 2 side-by-side
	// streams of independently executable instructions, one for data in xmm0-7, the other using xmm8-15.
	#define SSE2_RADIX16_DIT_0TWIDDLE(Xin0,Xin1,Xin2,Xin3,Xin4,Xin5,Xin6,Xin7,Xin8,Xin9,Xina,Xinb,Xinc,Xind,Xine,Xinf, Xisrt2, Xout0,Xo1,Xo2,Xo3,Xo4)\
	{\
	__asm__ volatile (\
	"pushl %%ebx	\n\t"/* Explicit save/restore of PIC register */\
		"movl	%[__in0],%%eax		\n\t"\
		"movl	%[__in1],%%ebx		\n\t"\
		"movl	%[__in2],%%ecx		\n\t"\
		"movl	%[__in3],%%edx		\n\t"\
		/* SSE2_RADIX4_DIT_0TWIDDLE_B(r0 ): */\
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
		"movl	%[__out0],%%esi		\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm2		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"subpd	%%xmm6,%%xmm3		\n\t"\
		"leal	0x10(%%esi),%%edi	\n\t"/* Need separate address Im parts of outputs due to literal-offsets below */\
		"movaps	%%xmm0,%c[__o2](%%esi)	\n\t"\
		"movaps	%%xmm2,%c[__o3](%%esi)	\n\t"\
		"movaps	%%xmm1,%c[__o2](%%edi)	\n\t"\
		"movaps	%%xmm3,%c[__o1](%%edi)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"addpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	%%xmm4,        (%%esi)	\n\t"\
		"movaps	%%xmm7,%c[__o1](%%esi)	\n\t"\
		"movaps	%%xmm5,        (%%edi)	\n\t"\
		"movaps	%%xmm6,%c[__o3](%%edi)	\n\t"\
		/* SSE2_RADIX4_DIT_0TWIDDLE_B(r8 ): */\
		"movl	%[__in4],%%eax		\n\t	leal	%c[__o4](%%esi),%%esi	\n\t"/* __out0 + 4*ostride */\
		"movl	%[__in5],%%ebx		\n\t"\
		"movl	%[__in6],%%ecx		\n\t"\
		"movl	%[__in7],%%edx		\n\t"\
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
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm2		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"subpd	%%xmm6,%%xmm3		\n\t"\
		"leal	0x10(%%esi),%%edi	\n\t"\
		"movaps	%%xmm0,%c[__o2](%%esi)	\n\t"\
		"movaps	%%xmm2,%c[__o3](%%esi)	\n\t"\
		"movaps	%%xmm1,%c[__o2](%%edi)	\n\t"\
		"movaps	%%xmm3,%c[__o1](%%edi)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"addpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	%%xmm4,        (%%esi)	\n\t"\
		"movaps	%%xmm7,%c[__o1](%%esi)	\n\t"\
		"movaps	%%xmm5,        (%%edi)	\n\t"\
		"movaps	%%xmm6,%c[__o3](%%edi)	\n\t"\
		/* SSE2_RADIX4_DIT_0TWIDDLE_B(r16): */\
		"movl	%[__in8],%%eax			\n\t	leal	%c[__o4](%%esi),%%esi	\n\t"/* __out0 + 8*ostride */\
		"movl	%[__in9],%%ebx			\n\t"\
		"movl	%[__ina],%%ecx			\n\t"\
		"movl	%[__inb],%%edx			\n\t"\
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
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm2		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"subpd	%%xmm6,%%xmm3		\n\t"\
		"leal	0x10(%%esi),%%edi	\n\t"\
		"movaps	%%xmm0,%c[__o2](%%esi)	\n\t"\
		"movaps	%%xmm2,%c[__o3](%%esi)	\n\t"\
		"movaps	%%xmm1,%c[__o2](%%edi)	\n\t"\
		"movaps	%%xmm3,%c[__o1](%%edi)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"addpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	%%xmm4,        (%%esi)	\n\t"\
		"movaps	%%xmm7,%c[__o1](%%esi)	\n\t"\
		"movaps	%%xmm5,        (%%edi)	\n\t"\
		"movaps	%%xmm6,%c[__o3](%%edi)	\n\t"\
		/* SSE2_RADIX4_DIT_0TWIDDLE_B(r24): */\
		"movl	%[__inc],%%eax			\n\t	leal	%c[__o4](%%esi),%%esi	\n\t"/* __out0 + c*ostride */\
		"movl	%[__ind],%%ebx			\n\t"\
		"movl	%[__ine],%%ecx			\n\t"\
		"movl	%[__inf],%%edx			\n\t"\
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
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm2		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"subpd	%%xmm6,%%xmm3		\n\t"\
		"leal	0x10(%%esi),%%edi	\n\t"\
		"movaps	%%xmm0,%c[__o2](%%esi)	\n\t"\
		"movaps	%%xmm2,%c[__o3](%%esi)	\n\t"\
		"movaps	%%xmm1,%c[__o2](%%edi)	\n\t"\
		"movaps	%%xmm3,%c[__o1](%%edi)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm2,%%xmm7		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"addpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	%%xmm4,        (%%esi)	\n\t"\
		"movaps	%%xmm7,%c[__o1](%%esi)	\n\t"\
		"movaps	%%xmm5,        (%%edi)	\n\t"\
		"movaps	%%xmm6,%c[__o3](%%edi)	\n\t"\
	/*** Now do 4 DFTs with internal twiddles on the 4*stride - separated data: ***/\
		"movl	%[__out0],%%eax		\n\t"\
		"leal	%c[__o4](%%eax),%%ebx	\n\t"/* __out0 +   [4*ostride] */\
		"leal	%c[__o4](%%ebx),%%ecx	\n\t"/* __out0 + 2*[4*ostride] */\
		"leal	%c[__o4](%%ecx),%%edx	\n\t"/* __out0 + 3*[4*ostride] */\
		/* Block 0: */\
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"movaps	    (%%ebx),%%xmm2	\n\t"\
		"movaps	0x10(%%ebx),%%xmm3	\n\t"\
		"subpd	    (%%ebx),%%xmm0	\n\t"\
		"subpd	0x10(%%ebx),%%xmm1	\n\t"\
		"addpd	    (%%eax),%%xmm2	\n\t"\
		"addpd	0x10(%%eax),%%xmm3	\n\t"\
		"movaps	    (%%ecx),%%xmm4	\n\t"\
		"movaps	0x10(%%ecx),%%xmm5	\n\t"\
		"movaps	    (%%edx),%%xmm6	\n\t"\
		"movaps	0x10(%%edx),%%xmm7	\n\t"\
		"subpd	    (%%edx),%%xmm4	\n\t"\
		"subpd	0x10(%%edx),%%xmm5	\n\t"\
		"addpd	    (%%ecx),%%xmm6	\n\t"\
		"addpd	0x10(%%ecx),%%xmm7	\n\t"\
		"subpd	%%xmm6,%%xmm2		\n\t"\
		"subpd	%%xmm7,%%xmm3		\n\t"\
		"movaps	%%xmm2,    (%%ecx)	\n\t"\
		"movaps	%%xmm3,0x10(%%ecx)	\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm2,%%xmm6		\n\t"\
		"addpd	%%xmm3,%%xmm7		\n\t"\
		"movaps	%%xmm6,    (%%eax)	\n\t"\
		"movaps	%%xmm7,0x10(%%eax)	\n\t"\
		"subpd	%%xmm5,%%xmm0		\n\t"\
		"subpd	%%xmm4,%%xmm1		\n\t"\
		"movaps	%%xmm0,    (%%edx)	\n\t"\
		"movaps	%%xmm1,0x10(%%ebx)	\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm0,%%xmm5		\n\t"\
		"addpd	%%xmm1,%%xmm4		\n\t"\
		"movaps	%%xmm5,    (%%ebx)	\n\t"\
		"movaps	%%xmm4,0x10(%%edx)	\n\t"\
		/* Block 2: */\
		"leal	%c[__o2](%%eax),%%eax	\n\t"/* All addresses += 2*ostride */\
		"leal	%c[__o2](%%ebx),%%ebx	\n\t"\
		"leal	%c[__o2](%%ecx),%%ecx	\n\t"\
		"leal	%c[__o2](%%edx),%%edx	\n\t"\
		"movl	%[__isrt2],%%edi	\n\t"\
		"movaps	(%%edi),%%xmm2		\n\t"\
		"movaps	    (%%ecx),%%xmm4	\n\t"\
		"movaps	0x10(%%ecx),%%xmm5	\n\t"\
		"movaps	    (%%edx),%%xmm0	\n\t"\
		"movaps	0x10(%%edx),%%xmm1	\n\t"\
		"addpd	0x10(%%ecx),%%xmm4	\n\t"\
		"subpd	    (%%ecx),%%xmm5	\n\t"\
		"subpd	0x10(%%edx),%%xmm0	\n\t"\
		"addpd	    (%%edx),%%xmm1	\n\t"\
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
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"movaps	    (%%ebx),%%xmm2	\n\t"\
		"movaps	0x10(%%ebx),%%xmm3	\n\t"\
		"subpd	0x10(%%ebx),%%xmm0	\n\t"\
		"subpd	    (%%ebx),%%xmm1	\n\t"\
		"addpd	    (%%eax),%%xmm3	\n\t"\
		"addpd	0x10(%%eax),%%xmm2	\n\t"\
		"subpd	%%xmm4,%%xmm3		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"movaps	%%xmm3,    (%%ecx)	\n\t"\
		"movaps	%%xmm1,0x10(%%ecx)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm3,%%xmm4		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"movaps	%%xmm4,    (%%eax)	\n\t"\
		"movaps	%%xmm5,0x10(%%eax)	\n\t"\
		"subpd	%%xmm7,%%xmm0		\n\t"\
		"subpd	%%xmm6,%%xmm2		\n\t"\
		"movaps	%%xmm0,    (%%edx)	\n\t"\
		"movaps	%%xmm2,0x10(%%ebx)	\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm7		\n\t"\
		"addpd	%%xmm2,%%xmm6		\n\t"\
		"movaps	%%xmm7,    (%%ebx)	\n\t"\
		"movaps	%%xmm6,0x10(%%edx)	\n\t"\
		/* Block 1: */\
		"leal	-%c[__o1](%%eax),%%eax	\n\t"/* All addresses -= 1*ostride */\
		"leal	-%c[__o1](%%ebx),%%ebx	\n\t"\
		"leal	-%c[__o1](%%ecx),%%ecx	\n\t"\
		"leal	-%c[__o1](%%edx),%%edx	\n\t"\
		"leal	0x10(%%edi),%%esi	\n\t"/* cc0 */\
		"movaps	    (%%edx),%%xmm0	\n\t"\
		"movaps	0x10(%%edx),%%xmm1	\n\t"\
		"movaps	    (%%edx),%%xmm2	\n\t"\
		"movaps	0x10(%%edx),%%xmm3	\n\t"\
		"mulpd	0x10(%%esi),%%xmm0	\n\t"\
		"mulpd	0x10(%%esi),%%xmm1	\n\t"\
		"mulpd	    (%%esi),%%xmm2	\n\t"\
		"mulpd	    (%%esi),%%xmm3	\n\t"\
		"subpd	%%xmm2,%%xmm1		\n\t"\
		"addpd	%%xmm3,%%xmm0		\n\t"\
		"movaps	    (%%ecx),%%xmm4	\n\t"\
		"movaps	0x10(%%ecx),%%xmm5	\n\t"\
		"movaps	    (%%ecx),%%xmm6	\n\t"\
		"movaps	0x10(%%ecx),%%xmm7	\n\t"\
		"mulpd	    (%%esi),%%xmm4	\n\t"\
		"mulpd	    (%%esi),%%xmm5	\n\t"\
		"mulpd	0x10(%%esi),%%xmm6	\n\t"\
		"mulpd	0x10(%%esi),%%xmm7	\n\t"\
		"subpd	%%xmm6,%%xmm5		\n\t"\
		"addpd	%%xmm7,%%xmm4		\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t"\
		"movaps	%%xmm4,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"subpd	%%xmm0,%%xmm6		\n\t"\
		"subpd	%%xmm1,%%xmm7		\n\t"\
		"movaps	    (%%ebx),%%xmm2	\n\t"\
		"movaps	0x10(%%ebx),%%xmm3	\n\t"\
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"addpd	0x10(%%ebx),%%xmm2	\n\t"\
		"subpd	    (%%ebx),%%xmm3	\n\t"\
		"mulpd	    (%%edi),%%xmm2	\n\t"\
		"mulpd	    (%%edi),%%xmm3	\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"subpd	%%xmm4,%%xmm2		\n\t"\
		"subpd	%%xmm5,%%xmm3		\n\t"\
		"movaps	%%xmm2,    (%%ecx)	\n\t"\
		"movaps	%%xmm3,0x10(%%ecx)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm2,%%xmm4		\n\t"\
		"addpd	%%xmm3,%%xmm5		\n\t"\
		"movaps	%%xmm4,    (%%eax)	\n\t"\
		"movaps	%%xmm5,0x10(%%eax)	\n\t"\
		"subpd	%%xmm7,%%xmm0		\n\t"\
		"subpd	%%xmm6,%%xmm1		\n\t"\
		"movaps	%%xmm0,    (%%edx)	\n\t"\
		"movaps	%%xmm1,0x10(%%ebx)	\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm7		\n\t"\
		"addpd	%%xmm1,%%xmm6		\n\t"\
		"movaps	%%xmm7,    (%%ebx)	\n\t"\
		"movaps	%%xmm6,0x10(%%edx)	\n\t"\
		/* Block 3: */\
		"leal	%c[__o2](%%eax),%%eax	\n\t"/* All addresses += 2*ostride */\
		"leal	%c[__o2](%%ebx),%%ebx	\n\t"\
		"leal	%c[__o2](%%ecx),%%ecx	\n\t"\
		"leal	%c[__o2](%%edx),%%edx	\n\t"\
		"movaps	    (%%edx),%%xmm0	\n\t"\
		"movaps	0x10(%%edx),%%xmm1	\n\t"\
		"movaps	    (%%edx),%%xmm2	\n\t"\
		"movaps	0x10(%%edx),%%xmm3	\n\t"\
		"mulpd	    (%%esi),%%xmm0	\n\t"\
		"mulpd	    (%%esi),%%xmm1	\n\t"\
		"mulpd	0x10(%%esi),%%xmm2	\n\t"\
		"mulpd	0x10(%%esi),%%xmm3	\n\t"\
		"subpd	%%xmm2,%%xmm1		\n\t"\
		"addpd	%%xmm3,%%xmm0		\n\t"\
		"movaps	    (%%ecx),%%xmm4	\n\t"\
		"movaps	0x10(%%ecx),%%xmm5	\n\t"\
		"movaps	    (%%ecx),%%xmm6	\n\t"\
		"movaps	0x10(%%ecx),%%xmm7	\n\t"\
		"mulpd	0x10(%%esi),%%xmm4	\n\t"\
		"mulpd	0x10(%%esi),%%xmm5	\n\t"\
		"mulpd	    (%%esi),%%xmm6	\n\t"\
		"mulpd	    (%%esi),%%xmm7	\n\t"\
		"subpd	%%xmm6,%%xmm5		\n\t"\
		"addpd	%%xmm7,%%xmm4		\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t"\
		"movaps	%%xmm4,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"subpd	%%xmm0,%%xmm6		\n\t"\
		"subpd	%%xmm1,%%xmm7		\n\t"\
		"movaps	    (%%ebx),%%xmm2	\n\t"\
		"movaps	0x10(%%ebx),%%xmm3	\n\t"\
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"subpd	0x10(%%ebx),%%xmm2	\n\t"\
		"addpd	    (%%ebx),%%xmm3	\n\t"\
		"mulpd	    (%%edi),%%xmm2	\n\t"\
		"mulpd	    (%%edi),%%xmm3	\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"subpd	%%xmm6,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm1		\n\t"\
		"movaps	%%xmm0,    (%%ecx)	\n\t"\
		"movaps	%%xmm1,0x10(%%ecx)	\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm0,%%xmm6		\n\t"\
		"addpd	%%xmm1,%%xmm7		\n\t"\
		"movaps	%%xmm6,    (%%eax)	\n\t"\
		"movaps	%%xmm7,0x10(%%eax)	\n\t"\
		"subpd	%%xmm5,%%xmm2		\n\t"\
		"subpd	%%xmm4,%%xmm3		\n\t"\
		"movaps	%%xmm2,    (%%edx)	\n\t"\
		"movaps	%%xmm3,0x10(%%ebx)	\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm2,%%xmm5		\n\t"\
		"addpd	%%xmm3,%%xmm4		\n\t"\
		"movaps	%%xmm5,    (%%ebx)	\n\t"\
		"movaps	%%xmm4,0x10(%%edx)	\n\t"\
	"popl %%ebx	\n\t"\
		:					/* outputs: none */\
		:[__in0] "m" (Xin0)	/* All inputs from memory addresses here */\
		,[__in1] "m" (Xin1)\
		,[__in2] "m" (Xin2)\
		,[__in3] "m" (Xin3)\
		,[__in4] "m" (Xin4)\
		,[__in5] "m" (Xin5)\
		,[__in6] "m" (Xin6)\
		,[__in7] "m" (Xin7)\
		,[__in8] "m" (Xin8)\
		,[__in9] "m" (Xin9)\
		,[__ina] "m" (Xina)\
		,[__inb] "m" (Xinb)\
		,[__inc] "m" (Xinc)\
		,[__ind] "m" (Xind)\
		,[__ine] "m" (Xine)\
		,[__inf] "m" (Xinf)\
		,[__isrt2] "m" (Xisrt2)\
		,[__out0] "m" (Xout0)\
		,[__o1] "e" (Xo1)\
		,[__o2] "e" (Xo2)\
		,[__o3] "e" (Xo3)\
		,[__o4] "e" (Xo4)\
		: "cc","memory","eax",/*"ebx",*/"ecx","edx","edi","esi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"	/* Clobbered registers */\
	);\
	}

	// Based on the SSE2_RADIX16_DIF_NOTWIDDLE macro in radix16_ditN_cy_dif1_gcc64.h, but with completely
	// specifiable 16-output addressing required for usage as the power-of-2 component of a twiddleless
	// radix = [odd*2^n] DFT routine.
	#define SSE2_RADIX16_DIF_0TWIDDLE(Xin0,Xi1,Xi2,Xi3,Xi4, Xisrt2, Xout0,Xout1,Xout2,Xout3,Xout4,Xout5,Xout6,Xout7,Xout8,Xout9,Xouta,Xoutb,Xoutc,Xoutd,Xoute,Xoutf)\
	{\
	__asm__ volatile (\
	"pushl %%ebx	\n\t"/* Explicit save/restore of PIC register */\
		/* SSE2_RADIX4_DIF_IN_PLACE(r1 , r17, r9 , r25): */\
		"movl	%[__in0],%%eax		\n\t"\
		"leal	%c[__i4](%%eax),%%ecx	\n\t"/* __in0 +   [4*istride]; note BR of [a,b,c,d]-ptrs, i.e. b/c swap */\
		"leal	%c[__i4](%%ecx),%%ebx	\n\t"/* __in0 + 2*[4*istride] */\
		"leal	%c[__i4](%%ebx),%%edx	\n\t"/* __in0 + 3*[4*istride] */\
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
		/* SSE2_RADIX4_DIF_IN_PLACE(r5 , r21, rdx, r29): */\
		"leal	%c[__i2](%%eax),%%eax	\n\t"/* All addresses += 2*ostride */\
		"leal	%c[__i2](%%ebx),%%ebx	\n\t"\
		"leal	%c[__i2](%%ecx),%%ecx	\n\t"\
		"leal	%c[__i2](%%edx),%%edx	\n\t"\
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
		/* SSE2_RADIX4_DIF_IN_PLACE(r3 , r19, rbx, r27): */\
		"leal	-%c[__i1](%%eax),%%eax	\n\t"/* All addresses -= 1*ostride */\
		"leal	-%c[__i1](%%ebx),%%ebx	\n\t"\
		"leal	-%c[__i1](%%ecx),%%ecx	\n\t"\
		"leal	-%c[__i1](%%edx),%%edx	\n\t"\
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
		/* SSE2_RADIX4_DIF_IN_PLACE(r7 , r23, r15, r31): */\
		"leal	%c[__i2](%%eax),%%eax	\n\t"/* All addresses += 2*ostride */\
		"leal	%c[__i2](%%ebx),%%ebx	\n\t"\
		"leal	%c[__i2](%%ecx),%%ecx	\n\t"\
		"leal	%c[__i2](%%edx),%%edx	\n\t"\
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
	/*** Now do 4 DFTs with internal twiddles on the 1*stride - separated data. Do blocks in order 0,2,1,3 to allow increment-only of rsi-datum from 1 block to the next: ***/\
		/* Block 0: r0-3 */\
		"movl	%[__in0],%%esi	\n\t"\
		"movl	%[__out0],%%eax		\n\t"\
		"movl	%[__out1],%%ebx		\n\t"\
		"movl	%[__out2],%%ecx		\n\t"\
		"movl	%[__out3],%%edx		\n\t"\
		"leal	0x10(%%esi),%%edi	\n\t"/* Need separate address Im parts of outputs due to literal-offsets below */\
		"movaps	        (%%esi),%%xmm0	\n\t"\
		"movaps	        (%%edi),%%xmm1	\n\t"\
		"movaps	%c[__i2](%%esi),%%xmm2	\n\t"\
		"movaps	%c[__i2](%%edi),%%xmm3	\n\t"\
		"subpd	%c[__i2](%%esi),%%xmm0	\n\t"\
		"subpd	%c[__i2](%%edi),%%xmm1	\n\t"\
		"addpd	        (%%esi),%%xmm2	\n\t"\
		"addpd	        (%%edi),%%xmm3	\n\t"\
		"movaps	%c[__i1](%%esi),%%xmm4	\n\t"\
		"movaps	%c[__i1](%%edi),%%xmm5	\n\t"\
		"movaps	%c[__i3](%%esi),%%xmm6	\n\t"\
		"movaps	%c[__i3](%%edi),%%xmm7	\n\t"\
		"subpd	%c[__i3](%%esi),%%xmm4	\n\t"\
		"subpd	%c[__i3](%%edi),%%xmm5	\n\t"\
		"addpd	%c[__i1](%%esi),%%xmm6	\n\t"\
		"addpd	%c[__i1](%%edi),%%xmm7	\n\t"\
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
		/* Block 2: */\
		"movl	%[__out8],%%eax		\n\t	leal	%c[__i4](%%esi),%%esi	\n\t"/* __in0 + 4*ostride */\
		"movl	%[__out9],%%ebx		\n\t"\
		"movl	%[__outa],%%ecx		\n\t"\
		"movl	%[__outb],%%edx		\n\t"\
		"leal	0x10(%%esi),%%edi	\n\t"\
		"movaps	%c[__i1](%%esi),%%xmm4	\n\t"\
		"movaps	%c[__i3](%%esi),%%xmm6	\n\t"\
		"movaps	%c[__i1](%%edi),%%xmm5	\n\t"\
		"movaps	%c[__i3](%%edi),%%xmm7	\n\t"\
		"movaps	%%xmm4,%%xmm0		\n\t"\
		"movaps	%%xmm6,%%xmm2		\n\t"\
		"movaps	%%xmm5,%%xmm1		\n\t"\
		"movaps	%%xmm7,%%xmm3		\n\t"\
		"movl	%[__isrt2],%%edi	\n\t"\
		"addl	$0x10,%%edi	\n\t"/* cc0 */\
		"mulpd	    (%%edi),%%xmm4	\n\t"\
		"mulpd	0x10(%%edi),%%xmm6	\n\t"\
		"mulpd	0x10(%%edi),%%xmm1	\n\t"\
		"mulpd	    (%%edi),%%xmm3	\n\t"\
		"mulpd	    (%%edi),%%xmm5	\n\t"\
		"mulpd	0x10(%%edi),%%xmm7	\n\t"\
		"mulpd	0x10(%%edi),%%xmm0	\n\t"\
		"mulpd	    (%%edi),%%xmm2	\n\t"\
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
		"leal	0x10(%%esi),%%edi	\n\t"\
		"movaps	%c[__i2](%%esi),%%xmm2	\n\t"\
		"movaps	%c[__i2](%%edi),%%xmm3	\n\t"\
		"subpd	%c[__i2](%%edi),%%xmm2	\n\t"\
		"addpd	%c[__i2](%%esi),%%xmm3	\n\t"\
		"movl	%[__isrt2],%%edi	\n\t"\
		"mulpd	(%%edi),%%xmm2	\n\t"/* mul by isrt2 */\
		"mulpd	(%%edi),%%xmm3	\n\t"\
		"leal	0x10(%%esi),%%edi	\n\t"\
		"movaps	        (%%esi),%%xmm0	\n\t"\
		"movaps	        (%%edi),%%xmm1	\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"addpd	        (%%esi),%%xmm2	\n\t"\
		"addpd	        (%%edi),%%xmm3	\n\t"\
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
		/* Block 1: r8-b */\
		"movl	%[__out4],%%eax		\n\t	leal	%c[__i4](%%esi),%%esi	\n\t"/* __in0 + 8*ostride */\
		"movl	%[__out5],%%ebx		\n\t"\
		"movl	%[__out6],%%ecx		\n\t"\
		"movl	%[__out7],%%edx		\n\t"\
		"leal	0x10(%%esi),%%edi	\n\t"\
		"movaps	        (%%esi),%%xmm0	\n\t"\
		"movaps	        (%%edi),%%xmm1	\n\t"\
		"movaps	%c[__i2](%%esi),%%xmm2	\n\t"\
		"movaps	%c[__i2](%%edi),%%xmm3	\n\t"\
		"subpd	%c[__i2](%%edi),%%xmm0	\n\t"\
		"subpd	%c[__i2](%%esi),%%xmm1	\n\t"\
		"addpd	        (%%edi),%%xmm2	\n\t"\
		"addpd	        (%%esi),%%xmm3	\n\t"\
		"movaps	%c[__i1](%%esi),%%xmm4	\n\t"\
		"movaps	%c[__i1](%%edi),%%xmm5	\n\t"\
		"movaps	%c[__i3](%%esi),%%xmm6	\n\t"\
		"movaps	%c[__i3](%%edi),%%xmm7	\n\t"\
		"subpd	%c[__i1](%%edi),%%xmm4	\n\t"\
		"addpd	%c[__i1](%%esi),%%xmm5	\n\t"\
		"addpd	%c[__i3](%%edi),%%xmm6	\n\t"\
		"subpd	%c[__i3](%%esi),%%xmm7	\n\t"\
		"movl	%[__isrt2],%%edi	\n\t"\
		"mulpd	(%%edi),%%xmm4		\n\t"\
		"mulpd	(%%edi),%%xmm5		\n\t"\
		"mulpd	(%%edi),%%xmm6		\n\t"\
		"mulpd	(%%edi),%%xmm7		\n\t"\
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
		/* Block 3: */\
		"movl	%[__outc],%%eax		\n\t	leal	%c[__i4](%%esi),%%esi	\n\t"/* __in0 + c*ostride */\
		"movl	%[__outd],%%ebx		\n\t"\
		"movl	%[__oute],%%ecx		\n\t"\
		"movl	%[__outf],%%edx		\n\t"\
		"leal	0x10(%%esi),%%edi	\n\t"\
		"movaps	%c[__i1](%%esi),%%xmm4	\n\t"\
		"movaps	%c[__i3](%%esi),%%xmm6	\n\t"\
		"movaps	%c[__i1](%%edi),%%xmm5	\n\t"\
		"movaps	%c[__i3](%%edi),%%xmm7	\n\t"\
		"movaps	%%xmm4,%%xmm0		\n\t"\
		"movaps	%%xmm6,%%xmm2		\n\t"\
		"movaps	%%xmm5,%%xmm1		\n\t"\
		"movaps	%%xmm7,%%xmm3		\n\t"\
		"movl	%[__isrt2],%%edi	\n\t"\
		"addl	$0x10,%%edi	\n\t"/* cc0 */\
		"mulpd	0x10(%%edi),%%xmm4	\n\t"\
		"mulpd	    (%%edi),%%xmm6	\n\t"\
		"mulpd	    (%%edi),%%xmm1	\n\t"\
		"mulpd	0x10(%%edi),%%xmm3	\n\t"\
		"mulpd	0x10(%%edi),%%xmm5	\n\t"\
		"mulpd	    (%%edi),%%xmm7	\n\t"\
		"mulpd	    (%%edi),%%xmm0	\n\t"\
		"mulpd	0x10(%%edi),%%xmm2	\n\t"\
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
		"leal	0x10(%%esi),%%edi	\n\t"\
		"movaps	%c[__i2](%%esi),%%xmm2	\n\t"\
		"movaps	%c[__i2](%%edi),%%xmm3	\n\t"\
		"addpd	%c[__i2](%%edi),%%xmm2	\n\t"\
		"subpd	%c[__i2](%%esi),%%xmm3	\n\t"\
		"movl	%[__isrt2],%%edi	\n\t"\
		"mulpd	(%%edi),%%xmm2	\n\t"/* mul by isrt2 */\
		"mulpd	(%%edi),%%xmm3	\n\t"\
		"leal	0x10(%%esi),%%edi	\n\t"\
		"movaps	        (%%esi),%%xmm0	\n\t"\
		"movaps	        (%%edi),%%xmm1	\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"addpd	        (%%esi),%%xmm2	\n\t"\
		"addpd	        (%%edi),%%xmm3	\n\t"\
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
		:[__in0] "m" (Xin0)	/* All inputs from memory addresses here */\
		,[__i1] "e" (Xi1)\
		,[__i2] "e" (Xi2)\
		,[__i3] "e" (Xi3)\
		,[__i4] "e" (Xi4)\
		,[__isrt2] "m" (Xisrt2)\
		,[__out0] "m" (Xout0)\
		,[__out1] "m" (Xout1)\
		,[__out2] "m" (Xout2)\
		,[__out3] "m" (Xout3)\
		,[__out4] "m" (Xout4)\
		,[__out5] "m" (Xout5)\
		,[__out6] "m" (Xout6)\
		,[__out7] "m" (Xout7)\
		,[__out8] "m" (Xout8)\
		,[__out9] "m" (Xout9)\
		,[__outa] "m" (Xouta)\
		,[__outb] "m" (Xoutb)\
		,[__outc] "m" (Xoutc)\
		,[__outd] "m" (Xoutd)\
		,[__oute] "m" (Xoute)\
		,[__outf] "m" (Xoutf)\
		: "cc","memory","eax",/*"ebx",*/"ecx","edx","edi","esi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"	/* Clobbered registers */\
	);\
	}

	/* With-twiddles out-of-place analog of above twiddleless DIT macro: 15 nontrivial complex input twiddles E1-f [E0 assumed = 1],
	The DIT version of this macro processes the twiddles in-order.
	NOTE: SINCE THIS MACRO IS SPECIFICALLY DESIGNED AS THE 2ND-PASS OF LARGE-POWER-OF-2-TWIDDLELESS DFT SYNTHESIS, THE
	"TWIDDLES" HERE ARE PURE OF THE DFT-INTERNAL VARIETY, AND THUS APPLIED TO THE INPUTS, JUST AS FOR THE ABOVE DIF COUNTERPART.

	Sincos layout: Two portions:

	Radix-16 shared consts anchored at isrt2:

	  isrt2 + 0x000;	cc0 + 0x010;	ss0 + 0x020;

	Per-block-specific set of 15 complex twiddles anchored at c1:

		c1  + 0x000;	s1  + 0x010;
		c2  + 0x020;	s2  + 0x030;
		c3  + 0x040;	s3  + 0x050;
		c4  + 0x060;	s4  + 0x070;
		c5  + 0x080;	s5  + 0x090;
		c6  + 0x0a0;	s6  + 0x0b0;
		c7  + 0x0c0;	s7  + 0x0d0;
		c8  + 0x0e0;	s8  + 0x0f0;
		c9  + 0x100;	s9  + 0x110;
		c10 + 0x120;	s10 + 0x130;
		c11 + 0x140;	s11 + 0x150;
		c12 + 0x160;	s12 + 0x170;
		c13 + 0x180;	s13 + 0x190;
		c14 + 0x1a0;	s14 + 0x1b0;
		c15 + 0x1c0;	s15 + 0x1d0;

	Use radix-16 DIF as template for DIT/OOP here, since need a pre-twiddles algorithm:
	*/
	#define SSE2_RADIX16_DIT_TWIDDLE_OOP(Xin0,Xi1,Xi2,Xi3,Xi4, Xout0,Xo1,Xo2,Xo3,Xo4, Xisrt2,Xc1)\
	{\
	__asm__ volatile (\
	"pushl %%ebx	\n\t"/* Explicit save/restore of PIC register */\
	/*...Block 0: Do in-place, i.e. outputs into __in0 + [0,1,2,3]*istride: */\
		"movl	%[__in0],%%eax		\n\t"\
		"leal	%c[__i1](%%eax),%%ecx	\n\t"/* __in0 +   istride */\
		"leal	%c[__i2](%%eax),%%ebx	\n\t"/* __in0 + 2*istride */\
		"leal	%c[__i3](%%eax),%%edx	\n\t"/* __in0 + 3*istride */\
		"/* Do	the p0,1 combo: */	\n\t"\
		"movl	%[__c1],%%esi 	/* c1 */\n\t"\
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	    (%%ecx),%%xmm4	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"movaps	0x10(%%ecx),%%xmm5	\n\t"\
		"movaps	%%xmm4,%%xmm6		\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t"\
		"mulpd	    (%%esi),%%xmm4	\n\t"\
		"mulpd	    (%%esi),%%xmm5	\n\t"\
		"mulpd	0x10(%%esi),%%xmm6	\n\t"\
		"movaps	%%xmm0,%%xmm2		\n\t"\
		"mulpd	0x10(%%esi),%%xmm7	\n\t"\
		"subpd	%%xmm6,%%xmm5		\n\t"\
		"movaps	%%xmm1,%%xmm3		\n\t"\
		"addpd	%%xmm7,%%xmm4		\n\t"\
		"addpd	%%xmm4,%%xmm0		\n\t"\
		"addpd	%%xmm5,%%xmm1		\n\t"\
		"subpd	%%xmm4,%%xmm2		\n\t"\
		"subpd	%%xmm5,%%xmm3		\n\t"\
		"/* Do	the p2,3 combo: */	\n\t"\
		"addl	$0x20,%%esi 	/* c2,3 */\n\t"\
		"movaps	    (%%edx),%%xmm4	\n\t"\
		"movaps	0x10(%%edx),%%xmm5	\n\t"\
		"movaps	    (%%edx),%%xmm6	\n\t"\
		"movaps	0x10(%%edx),%%xmm7	\n\t"\
		"mulpd	0x20(%%esi),%%xmm4	/* c3 */\n\t"\
		"mulpd	0x20(%%esi),%%xmm5	\n\t"\
		"mulpd	0x30(%%esi),%%xmm6	\n\t"\
		"mulpd	0x30(%%esi),%%xmm7	\n\t"\
		"subpd	%%xmm6,%%xmm5		\n\t"\
		"addpd	%%xmm7,%%xmm4		\n\t"\
		"movaps	%%xmm5,0x10(%%eax)	\n\t"\
		"movaps	%%xmm4,    (%%eax)	\n\t"\
		"movaps	    (%%ebx),%%xmm4	\n\t"\
		"movaps	0x10(%%ebx),%%xmm5	\n\t"\
		"movaps	    (%%ebx),%%xmm6	\n\t"\
		"movaps	0x10(%%ebx),%%xmm7	\n\t"\
		"mulpd	    (%%esi),%%xmm4	/* c2 */\n\t"\
		"mulpd	    (%%esi),%%xmm5	\n\t"\
		"mulpd	0x10(%%esi),%%xmm6	\n\t"\
		"mulpd	0x10(%%esi),%%xmm7	\n\t"\
		"subpd	%%xmm6,%%xmm5		\n\t"\
		"addpd	%%xmm7,%%xmm4		\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t"\
		"movaps	%%xmm4,%%xmm6		\n\t"\
		"subpd	    (%%eax),%%xmm4	\n\t"\
		"subpd	0x10(%%eax),%%xmm5	\n\t"\
		"addpd	    (%%eax),%%xmm6	\n\t"\
		"addpd	0x10(%%eax),%%xmm7	\n\t"\
		"/* Finish radix-4 butterfly and store results into temporary-array slots: */\n\t"\
		"subpd	%%xmm6,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm2		\n\t"\
		"subpd	%%xmm7,%%xmm1		\n\t"\
		"subpd	%%xmm4,%%xmm3		\n\t"\
		/* DIT has outputs (indexed in real-temp form as 0-7) 2/6,3/7 swapped, i.e. swap oregs c/d vs DIF: */\
		"movaps	%%xmm0,    (%%ebx)	\n\t"\
		"movaps	%%xmm2,    (%%edx)	\n\t"\
		"movaps	%%xmm1,0x10(%%ebx)	\n\t"\
		"movaps	%%xmm3,0x10(%%ecx)	\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm0,%%xmm6		\n\t"\
		"addpd	%%xmm2,%%xmm5		\n\t"\
		"addpd	%%xmm1,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm4		\n\t"\
		"movaps	%%xmm6,    (%%eax)	\n\t"\
		"movaps	%%xmm5,    (%%ecx)	\n\t"\
		"movaps	%%xmm7,0x10(%%eax)	\n\t"\
		"movaps	%%xmm4,0x10(%%edx)	\n\t"\
		"\n\t"\
	/*...Block 1: outputs into __in0 + [4,5,6,7]*istride: */\
		"addl	$%c[__i4],%%eax	\n\t"/* __in0 + 4*istride */\
		"addl	$%c[__i4],%%ecx	\n\t"/* __in0 + 5*istride */\
		"addl	$%c[__i4],%%ebx	\n\t"/* __in0 + 6*istride */\
		"addl	$%c[__i4],%%edx	\n\t"/* __in0 + 7*istride */\
		"/* Do	the p0,1 combo: */	\n\t"\
		"addl	$0x40,%%esi 	/* c4,5 */\n\t"\
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	    (%%ecx),%%xmm4	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"movaps	0x10(%%ecx),%%xmm5	\n\t"\
		"movaps	    (%%esi),%%xmm6	\n\t"\
		"movaps	0x10(%%esi),%%xmm7	\n\t"\
		"movaps	%%xmm0,%%xmm2		\n\t"\
		"movaps	%%xmm1,%%xmm3		\n\t"\
		"mulpd	%%xmm6,%%xmm0		/* c4 */\n\t"\
		"mulpd	%%xmm6,%%xmm1		\n\t"\
		"mulpd	%%xmm7,%%xmm2		\n\t"\
		"mulpd	%%xmm7,%%xmm3		\n\t"\
		"movaps	%%xmm4,%%xmm6		\n\t"\
		"subpd	%%xmm2,%%xmm1		\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t"\
		"mulpd	0x20(%%esi),%%xmm4	/* c5 */\n\t"\
		"addpd	%%xmm3,%%xmm0		\n\t"\
		"mulpd	0x20(%%esi),%%xmm5	\n\t"\
		"mulpd	0x30(%%esi),%%xmm6	\n\t"\
		"movaps	%%xmm0,%%xmm2		\n\t"\
		"mulpd	0x30(%%esi),%%xmm7	\n\t"\
		"subpd	%%xmm6,%%xmm5		\n\t"\
		"movaps	%%xmm1,%%xmm3		\n\t"\
		"addpd	%%xmm7,%%xmm4		\n\t"\
		"addpd	%%xmm4,%%xmm0		\n\t"\
		"addpd	%%xmm5,%%xmm1		\n\t"\
		"subpd	%%xmm4,%%xmm2		\n\t"\
		"subpd	%%xmm5,%%xmm3		\n\t"\
		"/* Do	the p2,3 combo: */	\n\t"\
		"addl	$0x40,%%esi 	/* c6,7 */\n\t"\
		"movaps	    (%%edx),%%xmm4	\n\t"\
		"movaps	0x10(%%edx),%%xmm5	\n\t"\
		"movaps	    (%%edx),%%xmm6	\n\t"\
		"movaps	0x10(%%edx),%%xmm7	\n\t"\
		"mulpd	0x20(%%esi),%%xmm4	/* c7 */\n\t"\
		"mulpd	0x20(%%esi),%%xmm5	\n\t"\
		"mulpd	0x30(%%esi),%%xmm6	\n\t"\
		"mulpd	0x30(%%esi),%%xmm7	\n\t"\
		"subpd	%%xmm6,%%xmm5		\n\t"\
		"addpd	%%xmm7,%%xmm4		\n\t"\
		"movaps	%%xmm5,0x10(%%eax)	\n\t"\
		"movaps	%%xmm4,    (%%eax)	\n\t"\
		"movaps	    (%%ebx),%%xmm4	\n\t"\
		"movaps	0x10(%%ebx),%%xmm5	\n\t"\
		"movaps	    (%%ebx),%%xmm6	\n\t"\
		"movaps	0x10(%%ebx),%%xmm7	\n\t"\
		"mulpd	    (%%esi),%%xmm4	/* c6 */\n\t"\
		"mulpd	    (%%esi),%%xmm5	\n\t"\
		"mulpd	0x10(%%esi),%%xmm6	\n\t"\
		"mulpd	0x10(%%esi),%%xmm7	\n\t"\
		"subpd	%%xmm6,%%xmm5		\n\t"\
		"addpd	%%xmm7,%%xmm4		\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t"\
		"movaps	%%xmm4,%%xmm6		\n\t"\
		"subpd	    (%%eax),%%xmm4	\n\t"\
		"subpd	0x10(%%eax),%%xmm5	\n\t"\
		"addpd	    (%%eax),%%xmm6	\n\t"\
		"addpd	0x10(%%eax),%%xmm7	\n\t"\
		"/* Finish radix-4 butterfly and store results into temporary-array slots: */\n\t"\
		"subpd	%%xmm6,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm2		\n\t"\
		"subpd	%%xmm7,%%xmm1		\n\t"\
		"subpd	%%xmm4,%%xmm3		\n\t"\
		"movaps	%%xmm0,    (%%ebx)	\n\t"\
		"movaps	%%xmm2,    (%%edx)	\n\t"\
		"movaps	%%xmm1,0x10(%%ebx)	\n\t"\
		"movaps	%%xmm3,0x10(%%ecx)	\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm0,%%xmm6		\n\t"\
		"addpd	%%xmm2,%%xmm5		\n\t"\
		"addpd	%%xmm1,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm4		\n\t"\
		"movaps	%%xmm6,    (%%eax)	\n\t"\
		"movaps	%%xmm5,    (%%ecx)	\n\t"\
		"movaps	%%xmm7,0x10(%%eax)	\n\t"\
		"movaps	%%xmm4,0x10(%%edx)	\n\t"\
		"\n\t"\
	/*...Block 2: outputs into __in0 + [8,9,a,b]*istride: */\
		"addl	$%c[__i4],%%eax	\n\t"/* __in0 + 8*istride */\
		"addl	$%c[__i4],%%ecx	\n\t"/* __in0 + 9*istride */\
		"addl	$%c[__i4],%%ebx	\n\t"/* __in0 + a*istride */\
		"addl	$%c[__i4],%%edx	\n\t"/* __in0 + b*istride */\
		"/* Do	the p0,1 combo: */	\n\t"\
		"addl	$0x40,%%esi 	/* c8,9 */\n\t"\
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	    (%%ecx),%%xmm4	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"movaps	0x10(%%ecx),%%xmm5	\n\t"\
		"movaps	    (%%esi),%%xmm6	\n\t"\
		"movaps	0x10(%%esi),%%xmm7	\n\t"\
		"movaps	%%xmm0,%%xmm2		\n\t"\
		"movaps	%%xmm1,%%xmm3		\n\t"\
		"mulpd	%%xmm6,%%xmm0		/* c8 */\n\t"\
		"mulpd	%%xmm6,%%xmm1		\n\t"\
		"mulpd	%%xmm7,%%xmm2		\n\t"\
		"mulpd	%%xmm7,%%xmm3		\n\t"\
		"movaps	%%xmm4,%%xmm6		\n\t"\
		"subpd	%%xmm2,%%xmm1		\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t"\
		"mulpd	0x20(%%esi),%%xmm4	/* c9 */\n\t"\
		"addpd	%%xmm3,%%xmm0		\n\t"\
		"mulpd	0x20(%%esi),%%xmm5	\n\t"\
		"mulpd	0x30(%%esi),%%xmm6	\n\t"\
		"movaps	%%xmm0,%%xmm2		\n\t"\
		"mulpd	0x30(%%esi),%%xmm7	\n\t"\
		"subpd	%%xmm6,%%xmm5		\n\t"\
		"movaps	%%xmm1,%%xmm3		\n\t"\
		"addpd	%%xmm7,%%xmm4		\n\t"\
		"addpd	%%xmm4,%%xmm0		\n\t"\
		"addpd	%%xmm5,%%xmm1		\n\t"\
		"subpd	%%xmm4,%%xmm2		\n\t"\
		"subpd	%%xmm5,%%xmm3		\n\t"\
		"/* Do	the p2,3 combo: */	\n\t"\
		"addl	$0x40,%%esi 	/* ca,b */\n\t"\
		"movaps	    (%%edx),%%xmm4	\n\t"\
		"movaps	0x10(%%edx),%%xmm5	\n\t"\
		"movaps	    (%%edx),%%xmm6	\n\t"\
		"movaps	0x10(%%edx),%%xmm7	\n\t"\
		"mulpd	0x20(%%esi),%%xmm4	/* cb */\n\t"\
		"mulpd	0x20(%%esi),%%xmm5	\n\t"\
		"mulpd	0x30(%%esi),%%xmm6	\n\t"\
		"mulpd	0x30(%%esi),%%xmm7	\n\t"\
		"subpd	%%xmm6,%%xmm5		\n\t"\
		"addpd	%%xmm7,%%xmm4		\n\t"\
		"movaps	%%xmm5,0x10(%%eax)	\n\t"\
		"movaps	%%xmm4,    (%%eax)	\n\t"\
		"movaps	    (%%ebx),%%xmm4	\n\t"\
		"movaps	0x10(%%ebx),%%xmm5	\n\t"\
		"movaps	    (%%ebx),%%xmm6	\n\t"\
		"movaps	0x10(%%ebx),%%xmm7	\n\t"\
		"mulpd	    (%%esi),%%xmm4	/* ca */\n\t"\
		"mulpd	    (%%esi),%%xmm5	\n\t"\
		"mulpd	0x10(%%esi),%%xmm6	\n\t"\
		"mulpd	0x10(%%esi),%%xmm7	\n\t"\
		"subpd	%%xmm6,%%xmm5		\n\t"\
		"addpd	%%xmm7,%%xmm4		\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t"\
		"movaps	%%xmm4,%%xmm6		\n\t"\
		"subpd	    (%%eax),%%xmm4	\n\t"\
		"subpd	0x10(%%eax),%%xmm5	\n\t"\
		"addpd	    (%%eax),%%xmm6	\n\t"\
		"addpd	0x10(%%eax),%%xmm7	\n\t"\
		"/* Finish radix-4 butterfly and store results into temporary-array slots: */\n\t"\
		"subpd	%%xmm6,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm2		\n\t"\
		"subpd	%%xmm7,%%xmm1		\n\t"\
		"subpd	%%xmm4,%%xmm3		\n\t"\
		"movaps	%%xmm0,    (%%ebx)	\n\t"\
		"movaps	%%xmm2,    (%%edx)	\n\t"\
		"movaps	%%xmm1,0x10(%%ebx)	\n\t"\
		"movaps	%%xmm3,0x10(%%ecx)	\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm0,%%xmm6		\n\t"\
		"addpd	%%xmm2,%%xmm5		\n\t"\
		"addpd	%%xmm1,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm4		\n\t"\
		"movaps	%%xmm6,    (%%eax)	\n\t"\
		"movaps	%%xmm5,    (%%ecx)	\n\t"\
		"movaps	%%xmm7,0x10(%%eax)	\n\t"\
		"movaps	%%xmm4,0x10(%%edx)	\n\t"\
		"\n\t"\
	/*...Block 3: outputs into __in0 + [c,d,e,f]*istride: */\
		"addl	$%c[__i4],%%eax	\n\t"/* __in0 + c*istride */\
		"addl	$%c[__i4],%%ecx	\n\t"/* __in0 + d*istride */\
		"addl	$%c[__i4],%%ebx	\n\t"/* __in0 + e*istride */\
		"addl	$%c[__i4],%%edx	\n\t"/* __in0 + f*istride */\
		"/* Do	the p0,1 combo: */	\n\t"\
		"addl	$0x40,%%esi 	/* cc,d */\n\t"\
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	    (%%ecx),%%xmm4	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"movaps	0x10(%%ecx),%%xmm5	\n\t"\
		"movaps	    (%%esi),%%xmm6	\n\t"\
		"movaps	0x10(%%esi),%%xmm7	\n\t"\
		"movaps	%%xmm0,%%xmm2		\n\t"\
		"movaps	%%xmm1,%%xmm3		\n\t"\
		"mulpd	%%xmm6,%%xmm0		/* cc */\n\t"\
		"mulpd	%%xmm6,%%xmm1		\n\t"\
		"mulpd	%%xmm7,%%xmm2		\n\t"\
		"mulpd	%%xmm7,%%xmm3		\n\t"\
		"movaps	%%xmm4,%%xmm6		\n\t"\
		"subpd	%%xmm2,%%xmm1		\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t"\
		"mulpd	0x20(%%esi),%%xmm4	/* cd */\n\t"\
		"addpd	%%xmm3,%%xmm0		\n\t"\
		"mulpd	0x20(%%esi),%%xmm5	\n\t"\
		"mulpd	0x30(%%esi),%%xmm6	\n\t"\
		"movaps	%%xmm0,%%xmm2		\n\t"\
		"mulpd	0x30(%%esi),%%xmm7	\n\t"\
		"subpd	%%xmm6,%%xmm5		\n\t"\
		"movaps	%%xmm1,%%xmm3		\n\t"\
		"addpd	%%xmm7,%%xmm4		\n\t"\
		"addpd	%%xmm4,%%xmm0		\n\t"\
		"addpd	%%xmm5,%%xmm1		\n\t"\
		"subpd	%%xmm4,%%xmm2		\n\t"\
		"subpd	%%xmm5,%%xmm3		\n\t"\
		"/* Do	the p2,3 combo: */	\n\t"\
		"addl	$0x40,%%esi 	/* ce,f */\n\t"\
		"movaps	    (%%edx),%%xmm4	\n\t"\
		"movaps	0x10(%%edx),%%xmm5	\n\t"\
		"movaps	    (%%edx),%%xmm6	\n\t"\
		"movaps	0x10(%%edx),%%xmm7	\n\t"\
		"mulpd	0x20(%%esi),%%xmm4	/* cf */\n\t"\
		"mulpd	0x20(%%esi),%%xmm5	\n\t"\
		"mulpd	0x30(%%esi),%%xmm6	\n\t"\
		"mulpd	0x30(%%esi),%%xmm7	\n\t"\
		"subpd	%%xmm6,%%xmm5		\n\t"\
		"addpd	%%xmm7,%%xmm4		\n\t"\
		"movaps	%%xmm5,0x10(%%eax)	\n\t"\
		"movaps	%%xmm4,    (%%eax)	\n\t"\
		"movaps	    (%%ebx),%%xmm4	\n\t"\
		"movaps	0x10(%%ebx),%%xmm5	\n\t"\
		"movaps	    (%%ebx),%%xmm6	\n\t"\
		"movaps	0x10(%%ebx),%%xmm7	\n\t"\
		"mulpd	    (%%esi),%%xmm4	/* ce */\n\t"\
		"mulpd	    (%%esi),%%xmm5	\n\t"\
		"mulpd	0x10(%%esi),%%xmm6	\n\t"\
		"mulpd	0x10(%%esi),%%xmm7	\n\t"\
		"subpd	%%xmm6,%%xmm5		\n\t"\
		"addpd	%%xmm7,%%xmm4		\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t"\
		"movaps	%%xmm4,%%xmm6		\n\t"\
		"subpd	    (%%eax),%%xmm4	\n\t"\
		"subpd	0x10(%%eax),%%xmm5	\n\t"\
		"addpd	    (%%eax),%%xmm6	\n\t"\
		"addpd	0x10(%%eax),%%xmm7	\n\t"\
		"/* Finish radix-4 butterfly and store results into temporary-array slots: */\n\t"\
		"subpd	%%xmm6,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm2		\n\t"\
		"subpd	%%xmm7,%%xmm1		\n\t"\
		"subpd	%%xmm4,%%xmm3		\n\t"\
		"movaps	%%xmm0,    (%%ebx)	\n\t"\
		"movaps	%%xmm2,    (%%edx)	\n\t"\
		"movaps	%%xmm1,0x10(%%ebx)	\n\t"\
		"movaps	%%xmm3,0x10(%%ecx)	\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm0,%%xmm6		\n\t"\
		"addpd	%%xmm2,%%xmm5		\n\t"\
		"addpd	%%xmm1,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm4		\n\t"\
		"movaps	%%xmm6,    (%%eax)	\n\t"\
		"movaps	%%xmm5,    (%%ecx)	\n\t"\
		"movaps	%%xmm7,0x10(%%eax)	\n\t"\
		"movaps	%%xmm4,0x10(%%edx)	\n\t"\
	/*************************************************************************************/\
	/*  And now do four more radix-4 transforms, including the internal twiddle factors: */\
	/*************************************************************************************/\
		"movl	%[__isrt2],%%esi 	\n\t"\
	/* Block 0: Combine 0-output of each radix-4, i.e. inputs from __in0 + [0,4,8,c]*istride: */\
		"movl	%[__in0],%%eax		\n\t"\
		"leal	%c[__i4](%%eax),%%ebx	\n\t"/* __in0 +   [4*istride] */\
		"leal	%c[__i4](%%ebx),%%ecx	\n\t"/* __in0 + 2*[4*istride] */\
		"leal	%c[__i4](%%ecx),%%edx	\n\t"/* __in0 + 3*[4*istride] */\
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	    (%%ecx),%%xmm4	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"movaps	0x10(%%ecx),%%xmm5	\n\t"\
		"movaps	    (%%ebx),%%xmm2	\n\t"\
		"movaps	    (%%edx),%%xmm6	\n\t"\
		"movaps	0x10(%%ebx),%%xmm3	\n\t"\
		"movaps	0x10(%%edx),%%xmm7	\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm6,%%xmm4		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"subpd	%%xmm7,%%xmm5		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm4,%%xmm6		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"addpd	%%xmm5,%%xmm7		\n\t"\
		"subpd	%%xmm6,%%xmm2		\n\t"\
		"subpd	%%xmm7,%%xmm3		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"movl	%[__out0],%%eax		\n\t"\
		"leal	%c[__o4](%%eax),%%ebx	\n\t"/* __out0 + 4*ostride */\
		"leal	%c[__o4](%%ebx),%%ecx	\n\t"/* __out0 + 8*ostride */\
		"leal	%c[__o4](%%ecx),%%edx	\n\t"/* __out0 + c*ostride */\
		"movaps	%%xmm2,    (%%ecx)	\n\t"\
		"movaps	%%xmm3,0x10(%%ecx)	\n\t"\
		"addpd	%%xmm2,	%%xmm6	\n\t"\
		"addpd	%%xmm3,	%%xmm7	\n\t"\
		"movaps	%%xmm6,    (%%eax)	\n\t"\
		"movaps	%%xmm7,0x10(%%eax)	\n\t"\
		"subpd	%%xmm5,	%%xmm0	\n\t"\
		"subpd	%%xmm4,	%%xmm1	\n\t"\
		"addpd	%%xmm5,	%%xmm5	\n\t"\
		"addpd	%%xmm4,	%%xmm4	\n\t"\
		"movaps	%%xmm0,    (%%edx)	\n\t"/* These 2 outputs [4/c] swapped w.r.to dif [2/3] due to +-I sign diff */\
		"movaps	%%xmm1,0x10(%%ebx)	\n\t"\
		"addpd	%%xmm0,	%%xmm5	\n\t"\
		"addpd	%%xmm1,	%%xmm4	\n\t"\
		"movaps	%%xmm5,    (%%ebx)	\n\t"\
		"movaps	%%xmm4,0x10(%%edx)	\n\t"\
	/* Block 1: Combine 1-output of each radix-4, i.e. inputs from __in0 + [1,5,9,d]*istride: */\
		"movl	%[__in0],%%eax		\n\t"\
		"leal	%c[__i4](%%eax),%%ebx	\n\t"/* __in0 +   [4*istride] */\
		"leal	%c[__i4](%%ebx),%%ecx	\n\t"/* __in0 + 2*[4*istride] */\
		"leal	%c[__i4](%%ecx),%%edx	\n\t"/* __in0 + 3*[4*istride] */\
		"addl	$%c[__i1],%%eax	\n\t"/* __in0 + 1*istride */\
		"addl	$%c[__i1],%%ecx	\n\t"/* __in0 + 5*istride */\
		"addl	$%c[__i1],%%ebx	\n\t"/* __in0 + 9*istride */\
		"addl	$%c[__i1],%%edx	\n\t"/* __in0 + d*istride */\
		"movaps	    (%%edx),%%xmm0	\n\t"\
		"movaps	0x10(%%edx),%%xmm1	\n\t"\
		"movaps	    (%%edx),%%xmm2	\n\t"\
		"movaps	0x10(%%edx),%%xmm3	\n\t"\
		"mulpd	0x20(%%esi),%%xmm0	\n\t"/* ss0 */\
		"mulpd	0x20(%%esi),%%xmm1	\n\t"\
		"mulpd	0x10(%%esi),%%xmm2	\n\t"/* cc0 */\
		"mulpd	0x10(%%esi),%%xmm3	\n\t"\
		"subpd	%%xmm2,%%xmm1		\n\t"\
		"addpd	%%xmm3,%%xmm0		\n\t"\
		"movaps	    (%%ecx),%%xmm4	\n\t"\
		"movaps	0x10(%%ecx),%%xmm5	\n\t"\
		"movaps	    (%%ecx),%%xmm6	\n\t"\
		"movaps	0x10(%%ecx),%%xmm7	\n\t"\
		"mulpd	0x10(%%esi),%%xmm4	\n\t"/* cc0 */\
		"mulpd	0x10(%%esi),%%xmm5	\n\t"\
		"mulpd	0x20(%%esi),%%xmm6	\n\t"/* ss0 */\
		"mulpd	0x20(%%esi),%%xmm7	\n\t"\
		"subpd	%%xmm6,%%xmm5		\n\t"\
		"addpd	%%xmm7,%%xmm4		\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t"\
		"movaps	%%xmm4,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"subpd	%%xmm0,%%xmm6		\n\t"\
		"subpd	%%xmm1,%%xmm7		\n\t"\
		"movaps	    (%%ebx),%%xmm2	\n\t"\
		"movaps	0x10(%%ebx),%%xmm3	\n\t"\
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"addpd	0x10(%%ebx),%%xmm2	\n\t"\
		"subpd	    (%%ebx),%%xmm3	\n\t"\
		"mulpd	    (%%esi),%%xmm2	\n\t"/* isrt2 */\
		"mulpd	    (%%esi),%%xmm3	\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"movl	%[__out0],%%eax		\n\t"\
		"leal	%c[__o4](%%eax),%%ebx	\n\t"/* __out0 + 4*ostride */\
		"leal	%c[__o4](%%ebx),%%ecx	\n\t"/* __out0 + 8*ostride */\
		"leal	%c[__o4](%%ecx),%%edx	\n\t"/* __out0 + c*ostride */\
		"addl	$%c[__o1],%%eax	\n\t"/* __out0 + 1*ostride */\
		"addl	$%c[__o1],%%ebx	\n\t"/* __out0 + 5*ostride */\
		"addl	$%c[__o1],%%ecx	\n\t"/* __out0 + 9*ostride */\
		"addl	$%c[__o1],%%edx	\n\t"/* __out0 + d*ostride */\
		"subpd	%%xmm4,%%xmm2		\n\t"\
		"subpd	%%xmm5,%%xmm3		\n\t"\
		"movaps	%%xmm2,    (%%ecx)	\n\t"\
		"movaps	%%xmm3,0x10(%%ecx)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm2,%%xmm4		\n\t"\
		"addpd	%%xmm3,%%xmm5		\n\t"\
		"movaps	%%xmm4,    (%%eax)	\n\t"\
		"movaps	%%xmm5,0x10(%%eax)	\n\t"\
		"subpd	%%xmm7,%%xmm0		\n\t"\
		"subpd	%%xmm6,%%xmm1		\n\t"\
		"movaps	%%xmm0,    (%%edx)	\n\t"\
		"movaps	%%xmm1,0x10(%%ebx)	\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm7		\n\t"\
		"addpd	%%xmm1,%%xmm6		\n\t"\
		"movaps	%%xmm7,    (%%ebx)	\n\t"\
		"movaps	%%xmm6,0x10(%%edx)	\n\t"\
	/* Block 2: Combine 2-output of each radix-4, i.e. inputs from __in0 + [2,6,a,e]*istride: */\
		"movaps	(%%esi),%%xmm2	/* isrt2 */\n\t"\
		"movl	%[__in0],%%eax		\n\t"\
		"leal	%c[__i4](%%eax),%%ebx	\n\t"/* __in0 +   [4*istride] */\
		"leal	%c[__i4](%%ebx),%%ecx	\n\t"/* __in0 + 2*[4*istride] */\
		"leal	%c[__i4](%%ecx),%%edx	\n\t"/* __in0 + 3*[4*istride] */\
		"addl	$%c[__i2],%%eax	\n\t"/* __in0 + 2*istride */\
		"addl	$%c[__i2],%%ebx	\n\t"/* __in0 + 6*istride */\
		"addl	$%c[__i2],%%ecx	\n\t"/* __in0 + a*istride */\
		"addl	$%c[__i2],%%edx	\n\t"/* __in0 + e*istride */\
		"movaps	    (%%ecx),%%xmm4	\n\t"\
		"movaps	0x10(%%ecx),%%xmm5	\n\t"\
		"movaps	    (%%edx),%%xmm0	\n\t"\
		"movaps	0x10(%%edx),%%xmm1	\n\t"\
		"addpd	0x10(%%ecx),%%xmm4	\n\t"\
		"subpd	    (%%ecx),%%xmm5	\n\t"\
		"subpd	0x10(%%edx),%%xmm0	\n\t"\
		"addpd	    (%%edx),%%xmm1	\n\t"\
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
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"movaps	    (%%ebx),%%xmm2	\n\t"\
		"movaps	0x10(%%ebx),%%xmm3	\n\t"\
		"subpd	0x10(%%ebx),%%xmm0	\n\t"\
		"subpd	    (%%ebx),%%xmm1	\n\t"\
		"addpd	    (%%eax),%%xmm3	\n\t"\
		"addpd	0x10(%%eax),%%xmm2	\n\t"\
		"movl	%[__out0],%%eax		\n\t"\
		"leal	%c[__o4](%%eax),%%ebx	\n\t"/* __out0 + 4*ostride */\
		"leal	%c[__o4](%%ebx),%%ecx	\n\t"/* __out0 + 8*ostride */\
		"leal	%c[__o4](%%ecx),%%edx	\n\t"/* __out0 + c*ostride */\
		"addl	$%c[__o2],%%eax	\n\t"/* __out0 + 2*ostride */\
		"addl	$%c[__o2],%%ebx	\n\t"/* __out0 + 6*ostride */\
		"addl	$%c[__o2],%%ecx	\n\t"/* __out0 + a*ostride */\
		"addl	$%c[__o2],%%edx	\n\t"/* __out0 + e*ostride */\
		"subpd	%%xmm4,%%xmm3		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"movaps	%%xmm3,    (%%ecx)	\n\t"\
		"movaps	%%xmm1,0x10(%%ecx)	\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm3,%%xmm4		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"movaps	%%xmm4,    (%%eax)	\n\t"\
		"movaps	%%xmm5,0x10(%%eax)	\n\t"\
		"subpd	%%xmm7,%%xmm0		\n\t"\
		"subpd	%%xmm6,%%xmm2		\n\t"\
		"movaps	%%xmm0,    (%%edx)	\n\t"\
		"movaps	%%xmm2,0x10(%%ebx)	\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm7		\n\t"\
		"addpd	%%xmm2,%%xmm6		\n\t"\
		"movaps	%%xmm7,    (%%ebx)	\n\t"\
		"movaps	%%xmm6,0x10(%%edx)	\n\t"\
	/* Block 3: Combine 3-output of each radix-4, i.e. inputs from __in0 + [3,7,b,f]*istride: */\
		"movl	%[__in0],%%eax		\n\t"\
		"leal	%c[__i4](%%eax),%%ebx	\n\t"/* __in0 +   [4*istride] */\
		"leal	%c[__i4](%%ebx),%%ecx	\n\t"/* __in0 + 2*[4*istride] */\
		"leal	%c[__i4](%%ecx),%%edx	\n\t"/* __in0 + 3*[4*istride] */\
		"addl	$%c[__i3],%%eax	\n\t"/* __in0 + 3*istride */\
		"addl	$%c[__i3],%%ecx	\n\t"/* __in0 + 7*istride */\
		"addl	$%c[__i3],%%ebx	\n\t"/* __in0 + b*istride */\
		"addl	$%c[__i3],%%edx	\n\t"/* __in0 + f*istride */\
		"movaps	    (%%edx),%%xmm0	\n\t"\
		"movaps	0x10(%%edx),%%xmm1	\n\t"\
		"movaps	    (%%edx),%%xmm2	\n\t"\
		"movaps	0x10(%%edx),%%xmm3	\n\t"\
		"mulpd	0x10(%%esi),%%xmm0	\n\t"/* cc0 */\
		"mulpd	0x10(%%esi),%%xmm1	\n\t"\
		"mulpd	0x20(%%esi),%%xmm2	\n\t"/* ss0 */\
		"mulpd	0x20(%%esi),%%xmm3	\n\t"\
		"subpd	%%xmm2,%%xmm1		\n\t"\
		"addpd	%%xmm3,%%xmm0		\n\t"\
		"movaps	    (%%ecx),%%xmm4	\n\t"\
		"movaps	0x10(%%ecx),%%xmm5	\n\t"\
		"movaps	    (%%ecx),%%xmm6	\n\t"\
		"movaps	0x10(%%ecx),%%xmm7	\n\t"\
		"mulpd	0x20(%%esi),%%xmm4	\n\t"/* ss0 */\
		"mulpd	0x20(%%esi),%%xmm5	\n\t"\
		"mulpd	0x10(%%esi),%%xmm6	\n\t"/* cc0 */\
		"mulpd	0x10(%%esi),%%xmm7	\n\t"\
		"subpd	%%xmm6,%%xmm5		\n\t"\
		"addpd	%%xmm7,%%xmm4		\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t"\
		"movaps	%%xmm4,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm4		\n\t"\
		"addpd	%%xmm1,%%xmm5		\n\t"\
		"subpd	%%xmm0,%%xmm6		\n\t"\
		"subpd	%%xmm1,%%xmm7		\n\t"\
		"movaps	    (%%ebx),%%xmm2	\n\t"\
		"movaps	0x10(%%ebx),%%xmm3	\n\t"\
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"subpd	0x10(%%ebx),%%xmm2	\n\t"\
		"addpd	    (%%ebx),%%xmm3	\n\t"\
		"mulpd	    (%%esi),%%xmm2	\n\t"/* isrt2 */\
		"mulpd	    (%%esi),%%xmm3	\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"movl	%[__out0],%%eax		\n\t"\
		"leal	%c[__o4](%%eax),%%ebx	\n\t"/* __out0 + 4*ostride */\
		"leal	%c[__o4](%%ebx),%%ecx	\n\t"/* __out0 + 8*ostride */\
		"leal	%c[__o4](%%ecx),%%edx	\n\t"/* __out0 + c*ostride */\
		"addl	$%c[__o3],%%eax	\n\t"/* __out0 + 3*ostride */\
		"addl	$%c[__o3],%%ecx	\n\t"/* __out0 + 7*ostride */\
		"addl	$%c[__o3],%%ebx	\n\t"/* __out0 + b*ostride */\
		"addl	$%c[__o3],%%edx	\n\t"/* __out0 + f*ostride */\
		"subpd	%%xmm6,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm1		\n\t"\
		"movaps	%%xmm0,    (%%ecx)	\n\t"\
		"movaps	%%xmm1,0x10(%%ecx)	\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm0,%%xmm6		\n\t"\
		"addpd	%%xmm1,%%xmm7		\n\t"\
		"movaps	%%xmm6,    (%%eax)	\n\t"\
		"movaps	%%xmm7,0x10(%%eax)	\n\t"\
		"subpd	%%xmm5,%%xmm2		\n\t"\
		"subpd	%%xmm4,%%xmm3		\n\t"\
		"movaps	%%xmm2,    (%%edx)	\n\t"\
		"movaps	%%xmm3,0x10(%%ebx)	\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm2,%%xmm5		\n\t"\
		"addpd	%%xmm3,%%xmm4		\n\t"\
		"movaps	%%xmm5,    (%%ebx)	\n\t"\
		"movaps	%%xmm4,0x10(%%edx)	\n\t"\
	"popl %%ebx	\n\t"\
		:					/* outputs: none */\
		: [__in0] "m" (Xin0)	/* All inputs from memory addresses here */\
		 ,[__i1] "e" (Xi1)\
		 ,[__i2] "e" (Xi2)\
		 ,[__i3] "e" (Xi3)\
		 ,[__i4] "e" (Xi4)\
		 ,[__out0] "m" (Xout0)\
		 ,[__o1] "e" (Xo1)\
		 ,[__o2] "e" (Xo2)\
		 ,[__o3] "e" (Xo3)\
		 ,[__o4] "e" (Xo4)\
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__c1] "m" (Xc1)\
		: "cc","memory","eax",/*"ebx",*/"ecx","edx","esi","edi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
	);\
	}

	// DIF version of above shares same sincos layout & data:
	#define SSE2_RADIX16_DIF_TWIDDLE_OOP(Xin0,Xi1,Xi2,Xi3,Xi4, Xout0,Xout1,Xout2,Xout3,Xout4,Xout5,Xout6,Xout7,Xout8,Xout9,Xouta,Xoutb,Xoutc,Xoutd,Xoute,Xoutf, Xisrt2,Xc1)\
	{\
	__asm__ volatile (\
	"pushl %%ebx	\n\t"/* Explicit save/restore of PIC register */\
	/*...Block 0: Do in-place, i.e. outputs into __in0 + [0,1,2,3]*istride: */\
		"movl	%[__in0],%%eax		\n\t"\
		"leal	%c[__i1](%%eax),%%ecx	\n\t"/* __in0 +   istride */\
		"leal	%c[__i2](%%eax),%%ebx	\n\t"/* __in0 + 2*istride */\
		"leal	%c[__i3](%%eax),%%edx	\n\t"/* __in0 + 3*istride */\
		"/* Do	the p0,1 combo: */	\n\t"\
		"movl	%[__c1],%%esi 	/* Roots sets c1-15 same as for DIT, w/c1 as base-ptr */\n\t"\
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	    (%%ecx),%%xmm4	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"movaps	0x10(%%ecx),%%xmm5	\n\t"\
		"movaps	%%xmm4,%%xmm6		\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t"\
		"mulpd	    (%%esi),%%xmm4	\n\t"\
		"mulpd	    (%%esi),%%xmm5	\n\t"\
		"mulpd	0x10(%%esi),%%xmm6	\n\t"\
		"movaps	%%xmm0,%%xmm2		\n\t"\
		"mulpd	0x10(%%esi),%%xmm7	\n\t"\
		"addpd	%%xmm6,%%xmm5		\n\t"\
		"movaps	%%xmm1,%%xmm3		\n\t"\
		"subpd	%%xmm7,%%xmm4		\n\t"\
		"addpd	%%xmm4,%%xmm0		\n\t"\
		"addpd	%%xmm5,%%xmm1		\n\t"\
		"subpd	%%xmm4,%%xmm2		\n\t"\
		"subpd	%%xmm5,%%xmm3		\n\t"\
		"/* Do	the p2,3 combo: */	\n\t"\
		"addl	$0x20,%%esi 	/* c2,3 */\n\t"\
		"movaps	    (%%edx),%%xmm4	\n\t"\
		"movaps	0x10(%%edx),%%xmm5	\n\t"\
		"movaps	    (%%edx),%%xmm6	\n\t"\
		"movaps	0x10(%%edx),%%xmm7	\n\t"\
		"mulpd	0x20(%%esi),%%xmm4	/* c3 */\n\t"\
		"mulpd	0x20(%%esi),%%xmm5	\n\t"\
		"mulpd	0x30(%%esi),%%xmm6	\n\t"\
		"mulpd	0x30(%%esi),%%xmm7	\n\t"\
		"addpd	%%xmm6,%%xmm5		\n\t"\
		"subpd	%%xmm7,%%xmm4		\n\t"\
		"movaps	%%xmm5,0x10(%%eax)	\n\t"\
		"movaps	%%xmm4,    (%%eax)	\n\t"\
		"movaps	    (%%ebx),%%xmm4	\n\t"\
		"movaps	0x10(%%ebx),%%xmm5	\n\t"\
		"movaps	    (%%ebx),%%xmm6	\n\t"\
		"movaps	0x10(%%ebx),%%xmm7	\n\t"\
		"mulpd	    (%%esi),%%xmm4	/* c2 */\n\t"\
		"mulpd	    (%%esi),%%xmm5	\n\t"\
		"mulpd	0x10(%%esi),%%xmm6	\n\t"\
		"mulpd	0x10(%%esi),%%xmm7	\n\t"\
		"addpd	%%xmm6,%%xmm5		\n\t"\
		"subpd	%%xmm7,%%xmm4		\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t"\
		"movaps	%%xmm4,%%xmm6		\n\t"\
		"subpd	    (%%eax),%%xmm4	\n\t"\
		"subpd	0x10(%%eax),%%xmm5	\n\t"\
		"addpd	    (%%eax),%%xmm6	\n\t"\
		"addpd	0x10(%%eax),%%xmm7	\n\t"\
		"/* Finish radix-4 butterfly and store results into temporary-array slots: */\n\t"\
		"subpd	%%xmm6,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm2		\n\t"\
		"subpd	%%xmm7,%%xmm1		\n\t"\
		"subpd	%%xmm4,%%xmm3		\n\t"\
		"movaps	%%xmm0,    (%%ecx)	\n\t"\
		"movaps	%%xmm2,    (%%ebx)	\n\t"\
		"movaps	%%xmm1,0x10(%%ecx)	\n\t"\
		"movaps	%%xmm3,0x10(%%edx)	\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm0,%%xmm6		\n\t"\
		"addpd	%%xmm2,%%xmm5		\n\t"\
		"addpd	%%xmm1,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm4		\n\t"\
		"movaps	%%xmm6,    (%%eax)	\n\t"\
		"movaps	%%xmm5,    (%%edx)	\n\t"\
		"movaps	%%xmm7,0x10(%%eax)	\n\t"\
		"movaps	%%xmm4,0x10(%%ebx)	\n\t"\
		"\n\t"\
	/*...Block 1: outputs into __in0 + [4,5,6,7]*istride: */\
		"addl	$%c[__i4],%%eax	\n\t"/* __in0 + 4*istride */\
		"addl	$%c[__i4],%%ecx	\n\t"/* __in0 + 5*istride */\
		"addl	$%c[__i4],%%ebx	\n\t"/* __in0 + 6*istride */\
		"addl	$%c[__i4],%%edx	\n\t"/* __in0 + 7*istride */\
		"/* Do	the p0,1 combo: */	\n\t"\
		"addl	$0x40,%%esi 	/* c4,5 */\n\t"\
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	    (%%ecx),%%xmm4	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"movaps	0x10(%%ecx),%%xmm5	\n\t"\
		"movaps	    (%%esi),%%xmm6	\n\t"\
		"movaps	0x10(%%esi),%%xmm7	\n\t"\
		"movaps	%%xmm0,%%xmm2		\n\t"\
		"movaps	%%xmm1,%%xmm3		\n\t"\
		"mulpd	%%xmm6,%%xmm0		/* c4 */\n\t"\
		"mulpd	%%xmm6,%%xmm1		\n\t"\
		"mulpd	%%xmm7,%%xmm2		\n\t"\
		"mulpd	%%xmm7,%%xmm3		\n\t"\
		"movaps	%%xmm4,%%xmm6		\n\t"\
		"addpd	%%xmm2,%%xmm1		\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t"\
		"mulpd	0x20(%%esi),%%xmm4	/* c5 */\n\t"\
		"subpd	%%xmm3,%%xmm0		\n\t"\
		"mulpd	0x20(%%esi),%%xmm5	\n\t"\
		"mulpd	0x30(%%esi),%%xmm6	\n\t"\
		"movaps	%%xmm0,%%xmm2		\n\t"\
		"mulpd	0x30(%%esi),%%xmm7	\n\t"\
		"addpd	%%xmm6,%%xmm5		\n\t"\
		"movaps	%%xmm1,%%xmm3		\n\t"\
		"subpd	%%xmm7,%%xmm4		\n\t"\
		"addpd	%%xmm4,%%xmm0		\n\t"\
		"addpd	%%xmm5,%%xmm1		\n\t"\
		"subpd	%%xmm4,%%xmm2		\n\t"\
		"subpd	%%xmm5,%%xmm3		\n\t"\
		"/* Do	the p2,3 combo: */	\n\t"\
		"addl	$0x40,%%esi 	/* c6,7 */\n\t"\
		"movaps	    (%%edx),%%xmm4	\n\t"\
		"movaps	0x10(%%edx),%%xmm5	\n\t"\
		"movaps	    (%%edx),%%xmm6	\n\t"\
		"movaps	0x10(%%edx),%%xmm7	\n\t"\
		"mulpd	0x20(%%esi),%%xmm4	/* c7 */\n\t"\
		"mulpd	0x20(%%esi),%%xmm5	\n\t"\
		"mulpd	0x30(%%esi),%%xmm6	\n\t"\
		"mulpd	0x30(%%esi),%%xmm7	\n\t"\
		"addpd	%%xmm6,%%xmm5		\n\t"\
		"subpd	%%xmm7,%%xmm4		\n\t"\
		"movaps	%%xmm5,0x10(%%eax)	\n\t"\
		"movaps	%%xmm4,    (%%eax)	\n\t"\
		"movaps	    (%%ebx),%%xmm4	\n\t"\
		"movaps	0x10(%%ebx),%%xmm5	\n\t"\
		"movaps	    (%%ebx),%%xmm6	\n\t"\
		"movaps	0x10(%%ebx),%%xmm7	\n\t"\
		"mulpd	    (%%esi),%%xmm4	/* c6 */\n\t"\
		"mulpd	    (%%esi),%%xmm5	\n\t"\
		"mulpd	0x10(%%esi),%%xmm6	\n\t"\
		"mulpd	0x10(%%esi),%%xmm7	\n\t"\
		"addpd	%%xmm6,%%xmm5		\n\t"\
		"subpd	%%xmm7,%%xmm4		\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t"\
		"movaps	%%xmm4,%%xmm6		\n\t"\
		"subpd	    (%%eax),%%xmm4	\n\t"\
		"subpd	0x10(%%eax),%%xmm5	\n\t"\
		"addpd	    (%%eax),%%xmm6	\n\t"\
		"addpd	0x10(%%eax),%%xmm7	\n\t"\
		"/* Finish radix-4 butterfly and store results into temporary-array slots: */\n\t"\
		"subpd	%%xmm6,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm2		\n\t"\
		"subpd	%%xmm7,%%xmm1		\n\t"\
		"subpd	%%xmm4,%%xmm3		\n\t"\
		"movaps	%%xmm0,    (%%ecx)	\n\t"\
		"movaps	%%xmm2,    (%%ebx)	\n\t"\
		"movaps	%%xmm1,0x10(%%ecx)	\n\t"\
		"movaps	%%xmm3,0x10(%%edx)	\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm0,%%xmm6		\n\t"\
		"addpd	%%xmm2,%%xmm5		\n\t"\
		"addpd	%%xmm1,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm4		\n\t"\
		"movaps	%%xmm6,    (%%eax)	\n\t"\
		"movaps	%%xmm5,    (%%edx)	\n\t"\
		"movaps	%%xmm7,0x10(%%eax)	\n\t"\
		"movaps	%%xmm4,0x10(%%ebx)	\n\t"\
		"\n\t"\
	/*...Block 2: outputs into __in0 + [8,9,a,b]*istride: */\
		"addl	$%c[__i4],%%eax	\n\t"/* __in0 + 8*istride */\
		"addl	$%c[__i4],%%ecx	\n\t"/* __in0 + 9*istride */\
		"addl	$%c[__i4],%%ebx	\n\t"/* __in0 + a*istride */\
		"addl	$%c[__i4],%%edx	\n\t"/* __in0 + b*istride */\
		"/* Do	the p0,1 combo: */	\n\t"\
		"addl	$0x40,%%esi 	/* c8,9 */\n\t"\
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	    (%%ecx),%%xmm4	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"movaps	0x10(%%ecx),%%xmm5	\n\t"\
		"movaps	    (%%esi),%%xmm6	\n\t"\
		"movaps	0x10(%%esi),%%xmm7	\n\t"\
		"movaps	%%xmm0,%%xmm2		\n\t"\
		"movaps	%%xmm1,%%xmm3		\n\t"\
		"mulpd	%%xmm6,%%xmm0		/* c8 */\n\t"\
		"mulpd	%%xmm6,%%xmm1		\n\t"\
		"mulpd	%%xmm7,%%xmm2		\n\t"\
		"mulpd	%%xmm7,%%xmm3		\n\t"\
		"movaps	%%xmm4,%%xmm6		\n\t"\
		"addpd	%%xmm2,%%xmm1		\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t"\
		"mulpd	0x20(%%esi),%%xmm4	/* c9 */\n\t"\
		"subpd	%%xmm3,%%xmm0		\n\t"\
		"mulpd	0x20(%%esi),%%xmm5	\n\t"\
		"mulpd	0x30(%%esi),%%xmm6	\n\t"\
		"movaps	%%xmm0,%%xmm2		\n\t"\
		"mulpd	0x30(%%esi),%%xmm7	\n\t"\
		"addpd	%%xmm6,%%xmm5		\n\t"\
		"movaps	%%xmm1,%%xmm3		\n\t"\
		"subpd	%%xmm7,%%xmm4		\n\t"\
		"addpd	%%xmm4,%%xmm0		\n\t"\
		"addpd	%%xmm5,%%xmm1		\n\t"\
		"subpd	%%xmm4,%%xmm2		\n\t"\
		"subpd	%%xmm5,%%xmm3		\n\t"\
		"/* Do	the p2,3 combo: */	\n\t"\
		"addl	$0x40,%%esi 	/* ca,b */\n\t"\
		"movaps	    (%%edx),%%xmm4	\n\t"\
		"movaps	0x10(%%edx),%%xmm5	\n\t"\
		"movaps	    (%%edx),%%xmm6	\n\t"\
		"movaps	0x10(%%edx),%%xmm7	\n\t"\
		"mulpd	0x20(%%esi),%%xmm4	/* cb */\n\t"\
		"mulpd	0x20(%%esi),%%xmm5	\n\t"\
		"mulpd	0x30(%%esi),%%xmm6	\n\t"\
		"mulpd	0x30(%%esi),%%xmm7	\n\t"\
		"addpd	%%xmm6,%%xmm5		\n\t"\
		"subpd	%%xmm7,%%xmm4		\n\t"\
		"movaps	%%xmm5,0x10(%%eax)	\n\t"\
		"movaps	%%xmm4,    (%%eax)	\n\t"\
		"movaps	    (%%ebx),%%xmm4	\n\t"\
		"movaps	0x10(%%ebx),%%xmm5	\n\t"\
		"movaps	    (%%ebx),%%xmm6	\n\t"\
		"movaps	0x10(%%ebx),%%xmm7	\n\t"\
		"mulpd	    (%%esi),%%xmm4	/* ca */\n\t"\
		"mulpd	    (%%esi),%%xmm5	\n\t"\
		"mulpd	0x10(%%esi),%%xmm6	\n\t"\
		"mulpd	0x10(%%esi),%%xmm7	\n\t"\
		"addpd	%%xmm6,%%xmm5		\n\t"\
		"subpd	%%xmm7,%%xmm4		\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t"\
		"movaps	%%xmm4,%%xmm6		\n\t"\
		"subpd	    (%%eax),%%xmm4	\n\t"\
		"subpd	0x10(%%eax),%%xmm5	\n\t"\
		"addpd	    (%%eax),%%xmm6	\n\t"\
		"addpd	0x10(%%eax),%%xmm7	\n\t"\
		"/* Finish radix-4 butterfly and store results into temporary-array slots: */\n\t"\
		"subpd	%%xmm6,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm2		\n\t"\
		"subpd	%%xmm7,%%xmm1		\n\t"\
		"subpd	%%xmm4,%%xmm3		\n\t"\
		"movaps	%%xmm0,    (%%ecx)	\n\t"\
		"movaps	%%xmm2,    (%%ebx)	\n\t"\
		"movaps	%%xmm1,0x10(%%ecx)	\n\t"\
		"movaps	%%xmm3,0x10(%%edx)	\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm0,%%xmm6		\n\t"\
		"addpd	%%xmm2,%%xmm5		\n\t"\
		"addpd	%%xmm1,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm4		\n\t"\
		"movaps	%%xmm6,    (%%eax)	\n\t"\
		"movaps	%%xmm5,    (%%edx)	\n\t"\
		"movaps	%%xmm7,0x10(%%eax)	\n\t"\
		"movaps	%%xmm4,0x10(%%ebx)	\n\t"\
		"\n\t"\
	/*...Block 3: outputs into __in0 + [c,d,e,f]*istride: */\
		"addl	$%c[__i4],%%eax	\n\t"/* __in0 + c*istride */\
		"addl	$%c[__i4],%%ecx	\n\t"/* __in0 + d*istride */\
		"addl	$%c[__i4],%%ebx	\n\t"/* __in0 + e*istride */\
		"addl	$%c[__i4],%%edx	\n\t"/* __in0 + f*istride */\
		"/* Do	the p0,1 combo: */	\n\t"\
		"addl	$0x40,%%esi 	/* cc,d */\n\t"\
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	    (%%ecx),%%xmm4	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"movaps	0x10(%%ecx),%%xmm5	\n\t"\
		"movaps	    (%%esi),%%xmm6	\n\t"\
		"movaps	0x10(%%esi),%%xmm7	\n\t"\
		"movaps	%%xmm0,%%xmm2		\n\t"\
		"movaps	%%xmm1,%%xmm3		\n\t"\
		"mulpd	%%xmm6,%%xmm0		/* cc */\n\t"\
		"mulpd	%%xmm6,%%xmm1		\n\t"\
		"mulpd	%%xmm7,%%xmm2		\n\t"\
		"mulpd	%%xmm7,%%xmm3		\n\t"\
		"movaps	%%xmm4,%%xmm6		\n\t"\
		"addpd	%%xmm2,%%xmm1		\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t"\
		"mulpd	0x20(%%esi),%%xmm4	/* cd */\n\t"\
		"subpd	%%xmm3,%%xmm0		\n\t"\
		"mulpd	0x20(%%esi),%%xmm5	\n\t"\
		"mulpd	0x30(%%esi),%%xmm6	\n\t"\
		"movaps	%%xmm0,%%xmm2		\n\t"\
		"mulpd	0x30(%%esi),%%xmm7	\n\t"\
		"addpd	%%xmm6,%%xmm5		\n\t"\
		"movaps	%%xmm1,%%xmm3		\n\t"\
		"subpd	%%xmm7,%%xmm4		\n\t"\
		"addpd	%%xmm4,%%xmm0		\n\t"\
		"addpd	%%xmm5,%%xmm1		\n\t"\
		"subpd	%%xmm4,%%xmm2		\n\t"\
		"subpd	%%xmm5,%%xmm3		\n\t"\
		"/* Do	the p2,3 combo: */	\n\t"\
		"addl	$0x40,%%esi 	/* ce,f */\n\t"\
		"movaps	    (%%edx),%%xmm4	\n\t"\
		"movaps	0x10(%%edx),%%xmm5	\n\t"\
		"movaps	    (%%edx),%%xmm6	\n\t"\
		"movaps	0x10(%%edx),%%xmm7	\n\t"\
		"mulpd	0x20(%%esi),%%xmm4	/* cf */\n\t"\
		"mulpd	0x20(%%esi),%%xmm5	\n\t"\
		"mulpd	0x30(%%esi),%%xmm6	\n\t"\
		"mulpd	0x30(%%esi),%%xmm7	\n\t"\
		"addpd	%%xmm6,%%xmm5		\n\t"\
		"subpd	%%xmm7,%%xmm4		\n\t"\
		"movaps	%%xmm5,0x10(%%eax)	\n\t"\
		"movaps	%%xmm4,    (%%eax)	\n\t"\
		"movaps	    (%%ebx),%%xmm4	\n\t"\
		"movaps	0x10(%%ebx),%%xmm5	\n\t"\
		"movaps	    (%%ebx),%%xmm6	\n\t"\
		"movaps	0x10(%%ebx),%%xmm7	\n\t"\
		"mulpd	    (%%esi),%%xmm4	/* ce */\n\t"\
		"mulpd	    (%%esi),%%xmm5	\n\t"\
		"mulpd	0x10(%%esi),%%xmm6	\n\t"\
		"mulpd	0x10(%%esi),%%xmm7	\n\t"\
		"addpd	%%xmm6,%%xmm5		\n\t"\
		"subpd	%%xmm7,%%xmm4		\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t"\
		"movaps	%%xmm4,%%xmm6		\n\t"\
		"subpd	    (%%eax),%%xmm4	\n\t"\
		"subpd	0x10(%%eax),%%xmm5	\n\t"\
		"addpd	    (%%eax),%%xmm6	\n\t"\
		"addpd	0x10(%%eax),%%xmm7	\n\t"\
		"/* Finish radix-4 butterfly and store results into temporary-array slots: */\n\t"\
		"subpd	%%xmm6,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm2		\n\t"\
		"subpd	%%xmm7,%%xmm1		\n\t"\
		"subpd	%%xmm4,%%xmm3		\n\t"\
		"movaps	%%xmm0,    (%%ecx)	\n\t"\
		"movaps	%%xmm2,    (%%ebx)	\n\t"\
		"movaps	%%xmm1,0x10(%%ecx)	\n\t"\
		"movaps	%%xmm3,0x10(%%edx)	\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm0,%%xmm6		\n\t"\
		"addpd	%%xmm2,%%xmm5		\n\t"\
		"addpd	%%xmm1,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm4		\n\t"\
		"movaps	%%xmm6,    (%%eax)	\n\t"\
		"movaps	%%xmm5,    (%%edx)	\n\t"\
		"movaps	%%xmm7,0x10(%%eax)	\n\t"\
		"movaps	%%xmm4,0x10(%%ebx)	\n\t"\
	/*************************************************************************************/\
	/*  And now do four more radix-4 transforms, including the internal twiddle factors: */\
	/*************************************************************************************/\
	/* Block 0: Combine 0-output of each radix-4, i.e. inputs from __in0 + [0,4,8,c]*istride: */\
		"movl	%[__in0],%%eax		\n\t"\
		"leal	%c[__i4](%%eax),%%ebx	\n\t"/* __in0 +   [4*istride] */\
		"leal	%c[__i4](%%ebx),%%ecx	\n\t"/* __in0 + 2*[4*istride] */\
		"leal	%c[__i4](%%ecx),%%edx	\n\t"/* __in0 + 3*[4*istride] */\
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	    (%%ecx),%%xmm4	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"movaps	0x10(%%ecx),%%xmm5	\n\t"\
		"movaps	    (%%ebx),%%xmm2	\n\t"\
		"movaps	    (%%edx),%%xmm6	\n\t"\
		"movaps	0x10(%%ebx),%%xmm3	\n\t"\
		"movaps	0x10(%%edx),%%xmm7	\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm6,%%xmm4		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"subpd	%%xmm7,%%xmm5		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm4,%%xmm6		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"addpd	%%xmm5,%%xmm7		\n\t"\
		"subpd	%%xmm6,%%xmm2		\n\t"\
		"subpd	%%xmm7,%%xmm3		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"movl	%[__out0],%%eax		\n\t"\
		"movl	%[__out1],%%ebx		\n\t"\
		"movl	%[__out2],%%ecx		\n\t"\
		"movl	%[__out3],%%edx		\n\t"\
	"prefetcht1	0x100(%%ebx)\n\t"\
		"movaps	%%xmm2,    (%%ebx)	\n\t"\
		"movaps	%%xmm3,0x10(%%ebx)	\n\t"\
		"addpd	%%xmm2,	%%xmm6	\n\t"\
		"addpd	%%xmm3,	%%xmm7	\n\t"\
		"movaps	%%xmm6,    (%%eax)	\n\t"\
		"movaps	%%xmm7,0x10(%%eax)	\n\t"\
		"subpd	%%xmm5,	%%xmm0	\n\t"\
		"subpd	%%xmm4,	%%xmm1	\n\t"\
		"addpd	%%xmm5,	%%xmm5	\n\t"\
		"addpd	%%xmm4,	%%xmm4	\n\t"\
		"movaps	%%xmm0,    (%%ecx)	\n\t"\
		"movaps	%%xmm1,0x10(%%edx)	\n\t"\
		"addpd	%%xmm0,	%%xmm5	\n\t"\
		"addpd	%%xmm1,	%%xmm4	\n\t"\
		"movaps	%%xmm5,    (%%edx)	\n\t"\
		"movaps	%%xmm4,0x10(%%ecx)	\n\t"\
	/* Block 2: Combine 2-output of each radix-4, i.e. inputs from __in0 + [4,5,6,7]*istride: */\
		"movl	%[__isrt2],%%esi 	\n\t"\
		"movaps	(%%esi),%%xmm3	/* isrt2 */\n\t"\
		"movl	%[__in0],%%eax		\n\t"\
		"leal	%c[__i4](%%eax),%%ebx	\n\t"/* __in0 +   [4*istride] */\
		"leal	%c[__i4](%%ebx),%%ecx	\n\t"/* __in0 + 2*[4*istride] */\
		"leal	%c[__i4](%%ecx),%%edx	\n\t"/* __in0 + 3*[4*istride] */\
		"addl	$%c[__i1],%%eax	\n\t"/* __in0 + 1*istride */\
		"addl	$%c[__i1],%%ebx	\n\t"/* __in0 + 5*istride */\
		"addl	$%c[__i1],%%ecx	\n\t"/* __in0 + 9*istride */\
		"addl	$%c[__i1],%%edx	\n\t"/* __in0 + d*istride */\
	"movl	%[__out3],%%edi		\n\t"\
	"prefetcht1	0x100(%%edi)\n\t"\
		"movaps	    (%%ecx),%%xmm4	\n\t"\
		"movaps	0x10(%%ecx),%%xmm5	\n\t"\
		"movaps	    (%%edx),%%xmm6	\n\t"\
		"movaps	0x10(%%edx),%%xmm7	\n\t"\
		"mulpd	%%xmm3,%%xmm4		\n\t"\
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"mulpd	%%xmm3,%%xmm5		\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"mulpd	%%xmm3,%%xmm6		\n\t"\
		"movaps	    (%%ebx),%%xmm2	\n\t"\
		"mulpd	%%xmm3,%%xmm7		\n\t"\
		"movaps	0x10(%%ebx),%%xmm3	\n\t"\
		"subpd	%%xmm3,%%xmm0		\n\t"\
		"subpd	%%xmm5,%%xmm4		\n\t"\
		"subpd	%%xmm2,%%xmm1		\n\t"\
		"subpd	%%xmm6,%%xmm7		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm0,%%xmm3		\n\t"\
		"addpd	%%xmm4,%%xmm5		\n\t"\
		"addpd	%%xmm1,%%xmm2		\n\t"\
		"addpd	%%xmm7,%%xmm6		\n\t"\
		"subpd	%%xmm6,%%xmm4		\n\t"\
		"subpd	%%xmm7,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm6		\n\t"\
		"addpd	%%xmm5,%%xmm7		\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm3		\n\t"\
		"subpd	%%xmm5,%%xmm2		\n\t"\
		"subpd	%%xmm6,%%xmm1		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"movl	%[__out4],%%eax		\n\t"\
		"movl	%[__out5],%%ebx		\n\t"\
		"movl	%[__out6],%%ecx		\n\t"\
		"movl	%[__out7],%%edx		\n\t"\
	"prefetcht1	0x100(%%ebx)\n\t"\
		"movaps	%%xmm0,    (%%ebx)	\n\t"\
		"movaps	%%xmm3,    (%%ecx)	\n\t"\
		"movaps	%%xmm2,0x10(%%ebx)	\n\t"\
		"movaps	%%xmm1,0x10(%%edx)	\n\t"\
		"addpd	%%xmm0,%%xmm4	\n\t"\
		"addpd	%%xmm3,%%xmm7	\n\t"\
		"addpd	%%xmm2,%%xmm5	\n\t"\
		"addpd	%%xmm1,%%xmm6	\n\t"\
		"movaps	%%xmm4,    (%%eax)	\n\t"\
		"movaps	%%xmm7,    (%%edx)	\n\t"\
		"movaps	%%xmm5,0x10(%%eax)	\n\t"\
		"movaps	%%xmm6,0x10(%%ecx)	\n\t"\
	/* Block 1: Combine 1-output of each radix-4, i.e. inputs from __in0 + [8,9,a,b]*istride: */\
		"movl	%[__in0],%%eax		\n\t"\
		"leal	%c[__i4](%%eax),%%ebx	\n\t"/* __in0 +   [4*istride] */\
		"leal	%c[__i4](%%ebx),%%ecx	\n\t"/* __in0 + 2*[4*istride] */\
		"leal	%c[__i4](%%ecx),%%edx	\n\t"/* __in0 + 3*[4*istride] */\
		"addl	$%c[__i2],%%eax	\n\t"/* __in0 + 2*istride */\
		"addl	$%c[__i2],%%ebx	\n\t"/* __in0 + 6*istride */\
		"addl	$%c[__i2],%%ecx	\n\t"/* __in0 + a*istride */\
		"addl	$%c[__i2],%%edx	\n\t"/* __in0 + e*istride */\
	"movl	%[__out7],%%edi		\n\t"\
	"prefetcht1	0x100(%%edi)\n\t"\
		"movaps	    (%%ecx),%%xmm4	\n\t"\
		"movaps	0x10(%%ecx),%%xmm5	\n\t"\
		"movaps	0x10(%%esi),%%xmm3	/* cc0, using isrt2 as base-ptr */\n\t"\
		"movaps	0x20(%%esi),%%xmm2	/* ss0, using isrt2 as base-ptr */\n\t"\
		"movaps	%%xmm4,%%xmm6		\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t"\
		"mulpd	%%xmm3,%%xmm4		\n\t"\
		"mulpd	%%xmm3,%%xmm5		\n\t"\
		"mulpd	%%xmm2,%%xmm6		\n\t"\
		"movaps	    (%%edx),%%xmm0	\n\t"\
		"mulpd	%%xmm2,%%xmm7		\n\t"\
		"movaps	0x10(%%edx),%%xmm1	\n\t"\
		"addpd	%%xmm6,%%xmm5		\n\t"\
		"movaps	%%xmm0,%%xmm6		\n\t"\
		"subpd	%%xmm7,%%xmm4		\n\t"\
		"movaps	%%xmm1,%%xmm7		\n\t"\
		"mulpd	%%xmm2,%%xmm6		\n\t"\
		"mulpd	%%xmm2,%%xmm7		\n\t"\
		"mulpd	%%xmm3,%%xmm0		\n\t"\
		"mulpd	%%xmm3,%%xmm1		\n\t"\
		"addpd	%%xmm0,%%xmm7		\n\t"\
		"subpd	%%xmm1,%%xmm6		\n\t"\
		"movaps	%%xmm4,%%xmm2		\n\t"\
		"movaps	%%xmm5,%%xmm3		\n\t"\
		"subpd	%%xmm6,%%xmm4		\n\t"\
		"subpd	%%xmm7,%%xmm5		\n\t"\
		"addpd	%%xmm2,%%xmm6		\n\t"\
		"addpd	%%xmm3,%%xmm7		\n\t"\
		"movaps	    (%%ebx),%%xmm2	\n\t"\
		"movaps	0x10(%%ebx),%%xmm3	\n\t"\
		"movaps	    (%%esi),%%xmm1	/* isrt2 */\n\t"\
		"movaps	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm2		\n\t"\
		"addpd	%%xmm0,%%xmm3		\n\t"\
		"mulpd	%%xmm1,%%xmm2		\n\t"\
		"mulpd	%%xmm1,%%xmm3		\n\t"\
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"subpd	%%xmm6,%%xmm2		\n\t"\
		"subpd	%%xmm5,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm3		\n\t"\
		"subpd	%%xmm4,%%xmm1		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"movl	%[__out8],%%eax		\n\t"\
		"movl	%[__out9],%%ebx		\n\t"\
		"movl	%[__outa],%%ecx		\n\t"\
		"movl	%[__outb],%%edx		\n\t"\
	"prefetcht1	0x100(%%ebx)\n\t"\
		"movaps	%%xmm2,    (%%ebx)	\n\t"\
		"movaps	%%xmm0,    (%%ecx)	\n\t"\
		"movaps	%%xmm3,0x10(%%ebx)	\n\t"\
		"movaps	%%xmm1,0x10(%%edx)	\n\t"\
		"addpd	%%xmm2,%%xmm6	\n\t"\
		"addpd	%%xmm0,%%xmm5	\n\t"\
		"addpd	%%xmm3,%%xmm7	\n\t"\
		"addpd	%%xmm1,%%xmm4	\n\t"\
		"movaps	%%xmm6,    (%%eax)	\n\t"\
		"movaps	%%xmm5,    (%%edx)	\n\t"\
		"movaps	%%xmm7,0x10(%%eax)	\n\t"\
		"movaps	%%xmm4,0x10(%%ecx)	\n\t"\
	/* Block 3: Combine 3-output of each radix-4, i.e. inputs from __in0 + [c,d,e,f]*istride: */\
		"movl	%[__in0],%%eax		\n\t"\
		"leal	%c[__i4](%%eax),%%ebx	\n\t"/* __in0 +   [4*istride] */\
		"leal	%c[__i4](%%ebx),%%ecx	\n\t"/* __in0 + 2*[4*istride] */\
		"leal	%c[__i4](%%ecx),%%edx	\n\t"/* __in0 + 3*[4*istride] */\
		"addl	$%c[__i3],%%eax	\n\t"/* __in0 + 3*istride */\
		"addl	$%c[__i3],%%ebx	\n\t"/* __in0 + 7*istride */\
		"addl	$%c[__i3],%%ecx	\n\t"/* __in0 + b*istride */\
		"addl	$%c[__i3],%%edx	\n\t"/* __in0 + f*istride */\
	"movl	%[__outb],%%edi		\n\t"\
	"prefetcht1	0x100(%%edi)\n\t"\
		"movaps	    (%%ecx),%%xmm4	\n\t"\
		"movaps	0x10(%%ecx),%%xmm5	\n\t"\
		"movaps	0x10(%%esi),%%xmm2	/* cc0, using isrt2 as base-ptr */\n\t"\
		"movaps	0x20(%%esi),%%xmm3	/* ss0, using isrt2 as base-ptr */\n\t"\
		"movaps	%%xmm4,%%xmm6		\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t"\
		"mulpd	%%xmm3,%%xmm4		\n\t"\
		"mulpd	%%xmm3,%%xmm5		\n\t"\
		"mulpd	%%xmm2,%%xmm6		\n\t"\
		"movaps	    (%%edx),%%xmm0	\n\t"\
		"mulpd	%%xmm2,%%xmm7		\n\t"\
		"movaps	0x10(%%edx),%%xmm1	\n\t"\
		"addpd	%%xmm6,%%xmm5		\n\t"\
		"movaps	%%xmm0,%%xmm6		\n\t"\
		"subpd	%%xmm7,%%xmm4		\n\t"\
		"movaps	%%xmm1,%%xmm7		\n\t"\
		"mulpd	%%xmm2,%%xmm6		\n\t"\
		"mulpd	%%xmm2,%%xmm7		\n\t"\
		"mulpd	%%xmm3,%%xmm0		\n\t"\
		"mulpd	%%xmm3,%%xmm1		\n\t"\
		"addpd	%%xmm0,%%xmm7		\n\t"\
		"subpd	%%xmm1,%%xmm6		\n\t"\
		"movaps	%%xmm4,%%xmm2		\n\t"\
		"movaps	%%xmm5,%%xmm3		\n\t"\
		"subpd	%%xmm6,%%xmm4		\n\t"\
		"subpd	%%xmm7,%%xmm5		\n\t"\
		"addpd	%%xmm2,%%xmm6		\n\t"\
		"addpd	%%xmm3,%%xmm7		\n\t"\
		"movaps	    (%%ebx),%%xmm2	\n\t"\
		"movaps	0x10(%%ebx),%%xmm3	\n\t"\
		"movaps	    (%%esi),%%xmm1		/* isrt2 */\n\t"\
		"movaps	%%xmm2,%%xmm0		\n\t"\
		"addpd	%%xmm3,%%xmm2		\n\t"\
		"subpd	%%xmm0,%%xmm3		\n\t"\
		"mulpd	%%xmm1,%%xmm2		\n\t"\
		"mulpd	%%xmm1,%%xmm3		\n\t"\
		"movaps	    (%%eax),%%xmm0	\n\t"\
		"movaps	0x10(%%eax),%%xmm1	\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t"\
		"addpd	%%xmm2,%%xmm2		\n\t"\
		"addpd	%%xmm3,%%xmm3		\n\t"\
		"addpd	%%xmm0,%%xmm2		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t"\
		"subpd	%%xmm4,%%xmm0		\n\t"\
		"subpd	%%xmm7,%%xmm2		\n\t"\
		"subpd	%%xmm5,%%xmm1		\n\t"\
		"subpd	%%xmm6,%%xmm3		\n\t"\
		"addpd	%%xmm4,%%xmm4		\n\t"\
		"addpd	%%xmm7,%%xmm7		\n\t"\
		"addpd	%%xmm5,%%xmm5		\n\t"\
		"addpd	%%xmm6,%%xmm6		\n\t"\
		"movl	%[__outc],%%eax		\n\t"\
		"movl	%[__outd],%%ebx		\n\t"\
		"movl	%[__oute],%%ecx		\n\t"\
		"movl	%[__outf],%%edx		\n\t"\
	"prefetcht1	0x100(%%ebx)\n\t"\
		"movaps	%%xmm0,    (%%ebx)	\n\t"\
		"movaps	%%xmm2,    (%%ecx)	\n\t"\
		"movaps	%%xmm1,0x10(%%ebx)	\n\t"\
		"movaps	%%xmm3,0x10(%%edx)	\n\t"\
		"addpd	%%xmm0,%%xmm4	\n\t"\
		"addpd	%%xmm2,%%xmm7	\n\t"\
		"addpd	%%xmm1,%%xmm5	\n\t"\
		"addpd	%%xmm3,%%xmm6	\n\t"\
		"movaps	%%xmm4,    (%%eax)	\n\t"\
		"movaps	%%xmm7,    (%%edx)	\n\t"\
		"movaps	%%xmm5,0x10(%%eax)	\n\t"\
		"movaps	%%xmm6,0x10(%%ecx)	\n\t"\
	"prefetcht1	0x100(%%edx)\n\t"\
	"popl %%ebx	\n\t"\
		:					/* outputs: none */\
		: [__in0] "m" (Xin0)	/* All inputs from memory addresses here */\
		 ,[__i1] "e" (Xi1)\
		 ,[__i2] "e" (Xi2)\
		 ,[__i3] "e" (Xi3)\
		 ,[__i4] "e" (Xi4)\
		 ,[__out0] "m" (Xout0)\
		 ,[__out1] "m" (Xout1)\
		 ,[__out2] "m" (Xout2)\
		 ,[__out3] "m" (Xout3)\
		 ,[__out4] "m" (Xout4)\
		 ,[__out5] "m" (Xout5)\
		 ,[__out6] "m" (Xout6)\
		 ,[__out7] "m" (Xout7)\
		 ,[__out8] "m" (Xout8)\
		 ,[__out9] "m" (Xout9)\
		 ,[__outa] "m" (Xouta)\
		 ,[__outb] "m" (Xoutb)\
		 ,[__outc] "m" (Xoutc)\
		 ,[__outd] "m" (Xoutd)\
		 ,[__oute] "m" (Xoute)\
		 ,[__outf] "m" (Xoutf)\
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__c1] "m" (Xc1)\
		: "cc","memory","eax",/*"ebx",*/"ecx","edx","esi","edi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
	);\
	}

#endif	/* sse2_macro_gcc_h_included */

