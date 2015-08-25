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
#ifndef radix16_wrapper_square_gcc_h_included
#define radix16_wrapper_square_gcc_h_included

	// Since the add0/1 block addresses advance in opposite directions, our prefetch scheme here is like so:
	// In fwd DIF prefetch  ahead 8 cache lines' worth of data w.r.to add0;
	// In inv DIT prefetch behind 8 cache lines' worth of data w.r.to add1;
	#define SSE2_RADIX16_WRAPPER_DIF(Xadd0,Xadd1,Xr1,Xr9,Xr17,Xr25,Xisrt2,Xcc0,Xc1,Xc2,Xc3,Xc4,Xc5,Xc6,Xc7,Xc8,Xc9,Xc10,Xc11,Xc12,Xc13,Xc14,Xc15)\
	{\
	__asm__ volatile (\
	"movl	%[__add0],%%edi	\n\t"\
	"prefetcht0	0x00(%%edi)\n\t"\
	"/*************************************************************/\n\t"\
	"/*                  1st set of inputs:                       */\n\t"\
	"/*************************************************************/\n\t"\
	"/*...Block 1: */\n\t"\
		"movl		%[__add0],%%eax\n\t"\
		"movl		%[__add1],%%esi\n\t"\
		"movl		%[__r1] ,%%ecx\n\t"\
		"movl		%[__c4] ,%%edx\n\t"\
		"\n\t"\
	"/* For interleaved [j1,j2] version, replace e.g.\n\t"\
		"\n\t"\
		"movaps	xmm0,[eax+0x40]	// a[jt+p4]\n\t"\
		"movaps	xmm1,[eax+0x50]	// a[jp+p4]\n\t"\
		"movaps	xmm2,[eax+0x40]	// xmm2 <- cpy a[jt+p4]\n\t"\
		"movaps	xmm3,[eax+0x50]	// xmm3 <- cpy a[jp+p4]\n\t"\
		"\n\t"\
	"by the following:\n\t"\
	"*/\n\t"\
	"/* Real parts: */\n\t"\
		"movaps		0x40(%%eax)	,%%xmm6		/* a[j1+p4], this is the scratch xmm register  */\n\t"\
		"movaps		0x40(%%esi)	,%%xmm2		/* a[j2+p4] */\n\t"\
		"movaps		%%xmm6		,%%xmm0		/* a[j1+p4] copy, his is the active  xmm register */\n\t"\
		"movaps		%%xmm2		,%%xmm3		/* a[j2+p4] copy */\n\t"\
		"unpckhpd	%%xmm2,%%xmm6\n\t"\
		"unpcklpd	%%xmm3,%%xmm0\n\t"\
		"movaps		%%xmm6		,0x140(%%ecx)	/* Store hi real in t21 */\n\t"\
	"/* Imag parts: */\n\t"\
		"movaps		0x50(%%eax)	,%%xmm7\n\t"\
		"movaps		0x50(%%esi)	,%%xmm4\n\t"\
		"movaps		%%xmm7,%%xmm1\n\t"\
		"movaps		%%xmm4,%%xmm5\n\t"\
		"unpckhpd	%%xmm4,%%xmm7\n\t"\
		"unpcklpd	%%xmm5,%%xmm1\n\t"\
		"movaps		%%xmm7		,0x150(%%ecx)	/* Store hi imag in t22 */\n\t"\
		"\n\t"\
		"movaps		%%xmm0		,%%xmm2	/* xmm2 <- cpy a[jt+p4] */\n\t"\
		"movaps		%%xmm1		,%%xmm3	/* xmm3 <- cpy a[jp+p4] */\n\t"\
	"/*************************************************************************************/\n\t"\
	"/******** From here on, things are identical to the code in radix16_dif_pass: ********/\n\t"\
	"/*************************************************************************************/\n\t"\
		"mulpd		    (%%edx)	,%%xmm0		/* a[jt+p4]*c4 */\n\t"\
		"mulpd		    (%%edx)	,%%xmm1		/* a[jp+p4]*c4 */\n\t"\
		"mulpd		0x10(%%edx)	,%%xmm2		/* a[jt+p4]*s4 */\n\t"\
		"mulpd		0x10(%%edx)	,%%xmm3		/* a[jp+p4]*s4 */\n\t"\
		"addpd		%%xmm2		,%%xmm1	/* xmm1 <- t6 */\n\t"\
		"subpd		%%xmm3		,%%xmm0	/* xmm0 <- t5 */\n\t"\
		"movaps		%%xmm1		,%%xmm3	/* xmm3 <- cpy t6 */\n\t"\
		"movaps		%%xmm0		,%%xmm2	/* xmm2 <- cpy t5 */\n\t"\
		"\n\t"\
		"movl		%[__c12],%%edx\n\t"\
	"/* Real parts: */\n\t"\
		"movaps		0xc0(%%eax)	,%%xmm6		/* a[j1+p12], this is the scratch xmm register  */\n\t"\
		"movaps		0xc0(%%eax)	,%%xmm4		/* a[j1+p12], this is the active  xmm register */\n\t"\
		"unpckhpd	0xc0(%%esi)	,%%xmm6		/* a[j2+p12] gets read twice */\n\t"\
		"unpcklpd	0xc0(%%esi)	,%%xmm4		/* a[jt+p12] */\n\t"\
		"movaps		%%xmm6		,0x160(%%ecx)	/* Store hi real in t23 */\n\t"\
	"/* Imag parts: */\n\t"\
		"movaps		0xd0(%%eax)	,%%xmm7\n\t"\
		"movaps		0xd0(%%eax)	,%%xmm5\n\t"\
		"unpckhpd	0xd0(%%esi)	,%%xmm7\n\t"\
		"unpcklpd	0xd0(%%esi)	,%%xmm5		/* a[jp+p12] */\n\t"\
		"movaps		%%xmm7		,0x170(%%ecx)	/* Store hi imag in t24 */\n\t"\
		"\n\t"\
		"movaps		%%xmm4		,%%xmm6	/* xmm6 <- cpy a[jt+p12] */\n\t"\
		"movaps		%%xmm5		,%%xmm7	/* xmm7 <- cpy a[jp+p12] */\n\t"\
		"\n\t"\
		"mulpd		    (%%edx)	,%%xmm4		/* a[jt+p12]*c12 */\n\t"\
		"mulpd		    (%%edx)	,%%xmm5		/* a[jp+p12]*c12 */\n\t"\
		"mulpd		0x10(%%edx)	,%%xmm6		/* a[jt+p12]*s12 */\n\t"\
		"mulpd		0x10(%%edx)	,%%xmm7		/* a[jp+p12]*s12 */\n\t"\
		"addpd		%%xmm6		,%%xmm5	/* xmm5 <- it */\n\t"\
		"subpd		%%xmm7		,%%xmm4	/* xmm4 <- rt		xmm6,7 free */\n\t"\
		"\n\t"\
		"addpd		%%xmm4		,%%xmm0	/* ~t5 <- t5 +rt */\n\t"\
		"addpd		%%xmm5		,%%xmm1	/* ~t6 <- t6 +it */\n\t"\
		"subpd		%%xmm4		,%%xmm2	/* ~t7 <- t5 -rt */\n\t"\
		"subpd		%%xmm5		,%%xmm3	/* ~t8 <- t6 -it	xmm4,5 free */\n\t"\
		"\n\t"\
		"/* Now do the p0,8 combo: */\n\t"\
		"movl		%[__c8] ,%%edx\n\t"\
	"/* Real parts: */\n\t"\
	"prefetcht0	0x20(%%edi)\n\t"\
		"movaps		0x80(%%eax)	,%%xmm6		/* a[j1+p8 ], this is the scratch xmm register  */\n\t"\
		"movaps		0x80(%%eax)	,%%xmm4		/* a[j1+p8 ], this is the active  xmm register */\n\t"\
		"unpckhpd	0x80(%%esi)	,%%xmm6		/* a[j2+p8 ] gets read twice */\n\t"\
		"unpcklpd	0x80(%%esi)	,%%xmm4		/* a[jt+p8 ] */\n\t"\
		"movaps			%%xmm6	,0x120(%%ecx)	/* Store hi real in t19 */\n\t"\
	"/* Imag parts: */\n\t"\
		"movaps		0x90(%%eax)	,%%xmm7\n\t"\
		"movaps		0x90(%%eax)	,%%xmm5\n\t"\
		"unpckhpd	0x90(%%esi)	,%%xmm7\n\t"\
		"unpcklpd	0x90(%%esi)	,%%xmm5		/* a[jp+p8 ] */\n\t"\
		"movaps		%%xmm7		,0x130(%%ecx)	/* Store hi imag in t20 */\n\t"\
		"\n\t"\
		"movaps		%%xmm4		,%%xmm6	/* xmm6 <- cpy a[jt+p8] */\n\t"\
		"movaps		%%xmm5		,%%xmm7	/* xmm7 <- cpy a[jp+p8] */\n\t"\
		"\n\t"\
		"mulpd		    (%%edx)	,%%xmm4		/* a[jt+p8]*c8 */\n\t"\
		"mulpd		    (%%edx)	,%%xmm5		/* a[jp+p8]*c8 */\n\t"\
		"mulpd		0x10(%%edx)	,%%xmm6		/* a[jt+p8]*s8 */\n\t"\
		"mulpd		0x10(%%edx)	,%%xmm7		/* a[jp+p8]*s8 */\n\t"\
		"addpd		%%xmm6		,%%xmm5	/* xmm5 <- it */\n\t"\
		"subpd		%%xmm7		,%%xmm4	/* xmm4 <- rt 	xmm6,7 free - stick t1,t2 in those */\n\t"\
		"\n\t"\
	"/* Real parts: */\n\t"\
		"movaps		    (%%eax)	,%%xmm6		/* a[j1    ], this is the scratch xmm register  */\n\t"\
		"movaps		    (%%eax)	,%%xmm7		/* a[j1    ], this is the active  xmm register */\n\t"\
		"unpckhpd	    (%%esi)	,%%xmm6		/* a[j2    ] gets read twice */\n\t"\
		"unpcklpd	    (%%esi)	,%%xmm7		/* a[jt] = t1*/\n\t"\
		"movaps		%%xmm6		,0x100(%%ecx)	/* Store hi real in t17 */\n\t"\
		"movaps		%%xmm7		,     (%%ecx)	/* Store active  in t1  */\n\t"\
	"/* Imag parts: */\n\t"\
		"movaps		0x10(%%eax)	,%%xmm6\n\t"\
		"movaps		0x10(%%eax)	,%%xmm7\n\t"\
		"unpckhpd	0x10(%%esi)	,%%xmm6\n\t"\
		"unpcklpd	0x10(%%esi)	,%%xmm7		/* a[jp] = t2*/\n\t"\
		"movaps		%%xmm6		,0x110(%%ecx)	/* Store hi imag in t18... */\n\t"\
		"movaps		    (%%ecx)	,%%xmm6		/* ...and reload t1. */\n\t"\
		"\n\t"\
		"subpd		%%xmm4		,%%xmm6	/* ~t3 <- t1 -rt */\n\t"\
		"subpd		%%xmm5		,%%xmm7	/* ~t4 <- t2 -it */\n\t"\
		"addpd		%%xmm4		,%%xmm4	/*          2*rt */\n\t"\
		"addpd		%%xmm5		,%%xmm5	/*          2*it */\n\t"\
		"addpd		%%xmm6		,%%xmm4	/* ~t1 <- t1 +rt */\n\t"\
		"addpd		%%xmm7		,%%xmm5	/* ~t2 <- t2 +it	xmm4,5 free */\n\t"\
		"\n\t"\
	"/* Finish radix-4 butterfly and store results into temporary-array slots: */\n\t"\
		"movl		%[__r1] ,%%eax\n\t"\
	"/*\n\t"\
	"~t5 =t1 -t5;		~t1 =t1 +t5;\n\t"\
	"~t6 =t2 -t6;		~t2 =t2 +t6;\n\t"\
	"*/\n\t"\
		"subpd		%%xmm0		,%%xmm4	/*~t5 =t1 -t5 */\n\t"\
		"subpd		%%xmm1		,%%xmm5	/*~t6 =t2 -t6 */\n\t"\
		"movaps		%%xmm4		,0x040(%%eax)	/* a[jt+p8 ] <- ~t5 */\n\t"\
		"movaps		%%xmm5		,0x050(%%eax)	/* a[jp+p8 ] <- ~t6 */\n\t"\
		"addpd		%%xmm0		,%%xmm0			/* 2*t5 */\n\t"\
		"addpd		%%xmm1		,%%xmm1			/* 2*t6 */\n\t"\
		"addpd		%%xmm4		,%%xmm0			/*~t1 =t1 +t5 */\n\t"\
		"addpd		%%xmm5		,%%xmm1			/*~t2 =t2 +t6 */\n\t"\
		"movaps		%%xmm0		,     (%%eax)	/* a[jt    ] <- ~t1 */\n\t"\
		"movaps		%%xmm1		,0x010(%%eax)	/* a[jp    ] <- ~t2 */\n\t"\
		"\n\t"\
	"/*\n\t"\
	"~t7 =t3 +t8;		~t3 =t3 -t8;\n\t"\
	"~t8 =t4 -t7;		~t4 =t4 +t7;\n\t"\
	"*/\n\t"\
		"subpd		%%xmm3		,%%xmm6	/*~t3 =t3 -t8 */\n\t"\
		"subpd		%%xmm2		,%%xmm7	/*~t8 =t4 -t7 */\n\t"\
		"movaps		%%xmm6		,0x020(%%eax)	/* a[jt+p4 ] <- ~t3 */\n\t"\
		"movaps		%%xmm7		,0x070(%%eax)	/* a[jp+p12] <- ~t8 */\n\t"\
		"addpd		%%xmm3		,%%xmm3			/* 2*t8 */\n\t"\
		"addpd		%%xmm2		,%%xmm2			/* 2*t7 */\n\t"\
		"addpd		%%xmm6		,%%xmm3			/*~t7 =t3 +t8 */\n\t"\
		"addpd		%%xmm7		,%%xmm2			/*~t4 =t4 +t7 */\n\t"\
		"movaps		%%xmm3		,0x060(%%eax)	/* a[jt+p12] <- ~t7 */\n\t"\
		"movaps		%%xmm2		,0x030(%%eax)	/* a[jp+p4 ] <- ~t4 */\n\t"\
		"\n\t"\
"/*...Block 2:		Cost: 46 MOVapd, 16 UNPCKHPD, 28 ADD/SUBpd, 16 MULpd */\n\t"\
		"\n\t"\
	"prefetcht0	0x40(%%edi)\n\t"\
		"movl		%[__add0],%%eax\n\t"\
		"movl		%[__add1],%%esi\n\t"\
		"movl		%[__r9] ,%%ecx\n\t"\
		"\n\t"\
"/* Do the p2,10 combo: */\n\t"\
		"movl		%[__c2] ,%%edx\n\t"\
	"/* Real parts: */\n\t"\
		"movaps		0x20(%%eax)	,%%xmm6		/* a[j1+p2 ], this is the scratch xmm register */\n\t"\
		"movaps		0x20(%%eax)	,%%xmm0		/* a[j1+p2 ], this is the active  xmm register */\n\t"\
		"unpckhpd	0x20(%%esi)	,%%xmm6		/* a[j2+p2 ] gets read twice */\n\t"\
		"unpcklpd	0x20(%%esi)	,%%xmm0		/* a[jt+p2 ] */\n\t"\
		"movaps			%%xmm6	,0x100(%%ecx)	/* Store hi real in t9 +16 */\n\t"\
	"/* Imag parts: */\n\t"\
		"movaps		0x30(%%eax)	,%%xmm7\n\t"\
		"movaps		0x30(%%eax)	,%%xmm1\n\t"\
		"unpckhpd	0x30(%%esi)	,%%xmm7\n\t"\
		"unpcklpd	0x30(%%esi)	,%%xmm1		/* a[jp+p2 ] */\n\t"\
		"movaps		%%xmm7		,0x110(%%ecx)	/* Store hi imag in t10+16 */\n\t"\
		"\n\t"\
		"movaps		%%xmm0		,%%xmm2	/* xmm2 <- cpy a[jt+p2] */\n\t"\
		"movaps		%%xmm1		,%%xmm3	/* xmm3 <- cpy a[jp+p2] */\n\t"\
		"\n\t"\
		"mulpd		    (%%edx)	,%%xmm0		/* a[jt+p2]*c2 */\n\t"\
		"mulpd		    (%%edx)	,%%xmm1		/* a[jp+p2]*c2 */\n\t"\
		"mulpd		0x10(%%edx)	,%%xmm2		/* a[jt+p2]*s2 */\n\t"\
		"mulpd		0x10(%%edx)	,%%xmm3		/* a[jp+p2]*s2 */\n\t"\
		"addpd		%%xmm2		,%%xmm1	/* xmm1 <- t10*/\n\t"\
		"subpd		%%xmm3		,%%xmm0	/* xmm0 <- t9 */\n\t"\
		"movaps		%%xmm1		,%%xmm3	/* xmm3 <- cpy t10*/\n\t"\
		"movaps		%%xmm0		,%%xmm2	/* xmm2 <- cpy t9 */\n\t"\
		"\n\t"\
		"movl		%[__c10],%%edx\n\t"\
	"/* Real parts: */\n\t"\
		"movaps		0xa0(%%eax)	,%%xmm6		/* a[j1+p10], this is the scratch xmm register  */\n\t"\
		"movaps		0xa0(%%eax)	,%%xmm4		/* a[j1+p10], this is the active  xmm register */\n\t"\
		"unpckhpd	0xa0(%%esi)	,%%xmm6		/* a[j2+p10] gets read twice */\n\t"\
		"unpcklpd	0xa0(%%esi)	,%%xmm4		/* a[jt+p10] */\n\t"\
		"movaps			%%xmm6	,0x120(%%ecx)	/* Store hi real in t11+16 */\n\t"\
	"/* Imag parts: */\n\t"\
		"movaps		0xb0(%%eax)	,%%xmm7\n\t"\
		"movaps		0xb0(%%eax)	,%%xmm5\n\t"\
		"unpckhpd	0xb0(%%esi)	,%%xmm7\n\t"\
		"unpcklpd	0xb0(%%esi)	,%%xmm5		/* a[jp+p10] */\n\t"\
		"movaps			%%xmm7	,0x130(%%ecx)	/* Store hi imag in t12+16 */\n\t"\
		"\n\t"\
		"movaps			%%xmm4	,%%xmm6	/* xmm6 <- cpy a[jt+p10] */\n\t"\
		"movaps			%%xmm5	,%%xmm7	/* xmm7 <- cpy a[jp+p10] */\n\t"\
		"\n\t"\
		"mulpd		    (%%edx)	,%%xmm4		/* a[jt+p10]*c10 */\n\t"\
		"mulpd		    (%%edx)	,%%xmm5		/* a[jp+p10]*c10 */\n\t"\
		"mulpd		0x10(%%edx)	,%%xmm6		/* a[jt+p10]*s10 */\n\t"\
		"mulpd		0x10(%%edx)	,%%xmm7		/* a[jp+p10]*s10 */\n\t"\
		"addpd		%%xmm6		,%%xmm5	/* xmm5 <- it */\n\t"\
		"subpd		%%xmm7		,%%xmm4	/* xmm4 <- rt		xmm6,7 free */\n\t"\
		"\n\t"\
		"addpd		%%xmm4		,%%xmm0	/* ~t13<- t13+rt */\n\t"\
		"addpd		%%xmm5		,%%xmm1	/* ~t14<- t14+it */\n\t"\
		"subpd		%%xmm4		,%%xmm2	/* ~t15<- t13-rt */\n\t"\
		"subpd		%%xmm5		,%%xmm3	/* ~t16<- t14-it	xmm4,5 free */\n\t"\
		"\n\t"\
"/* Do the p6,14 combo - do p14 first so register assignments come out in same relative order as for p2,10 */\n\t"\
		"movl		%[__c14],%%edx\n\t"\
	"/* Real parts: */\n\t"\
		"movaps		0xe0(%%eax)	,%%xmm6		/* a[j1+p14], this is the scratch xmm register  */\n\t"\
		"movaps		0xe0(%%eax)	,%%xmm4		/* a[j1+p14], this is the active  xmm register */\n\t"\
		"unpckhpd	0xe0(%%esi)	,%%xmm6		/* a[j2+p14] gets read twice */\n\t"\
		"unpcklpd	0xe0(%%esi)	,%%xmm4		/* a[jt+p14] */\n\t"\
		"movaps			%%xmm6	,0x160(%%ecx)	/* Store hi real in t15+16 */\n\t"\
	"/* Imag parts: */\n\t"\
		"movaps		0xf0(%%eax)	,%%xmm7\n\t"\
		"movaps		0xf0(%%eax)	,%%xmm5\n\t"\
		"unpckhpd	0xf0(%%esi)	,%%xmm7\n\t"\
		"unpcklpd	0xf0(%%esi)	,%%xmm5		/* a[jp+p14] */\n\t"\
		"movaps			%%xmm7	,0x170(%%ecx)	/* Store hi imag in t16+16 */\n\t"\
		"\n\t"\
		"movaps			%%xmm4	,%%xmm6		/* xmm6 <- cpy a[jt+p14] */\n\t"\
		"movaps			%%xmm5	,%%xmm7		/* xmm7 <- cpy a[jp+p14] */\n\t"\
		"\n\t"\
		"mulpd		    (%%edx)	,%%xmm4		/* a[jt+p14]*c14 */\n\t"\
		"mulpd		    (%%edx)	,%%xmm5		/* a[jp+p14]*c14 */\n\t"\
		"mulpd		0x10(%%edx)	,%%xmm6		/* a[jt+p14]*s14 */\n\t"\
		"mulpd		0x10(%%edx)	,%%xmm7		/* a[jp+p14]*s14 */\n\t"\
		"addpd		%%xmm6		,%%xmm5		/* xmm5 <- it */\n\t"\
		"subpd		%%xmm7		,%%xmm4		/* xmm4 <- rt		xmm6,7 free */\n\t"\
		"movaps		%%xmm5		,0x070(%%ecx)	/* Store it in t16*/\n\t"\
		"movaps		%%xmm4		,0x060(%%ecx)	/* Store rt in t15*/\n\t"\
		"\n\t"\
	"prefetcht0	0x60(%%edi)\n\t"\
		"movl		%[__c6] ,%%edx\n\t"\
	"/* Real parts: */\n\t"\
		"movaps		0x60(%%eax)	,%%xmm6		/* a[j1+p6 ], this is the scratch xmm register  */\n\t"\
		"movaps		0x60(%%eax)	,%%xmm4		/* a[j1+p6 ], this is the active  xmm register */\n\t"\
		"unpckhpd	0x60(%%esi)	,%%xmm6		/* a[j2+p6 ] gets read twice */\n\t"\
		"unpcklpd	0x60(%%esi)	,%%xmm4		/* a[jt+p6 ] */\n\t"\
		"movaps			%%xmm6	,0x140(%%ecx)	/* Store hi real in t13+16 */\n\t"\
	"/* Imag parts: */\n\t"\
		"movaps		0x70(%%eax)	,%%xmm7\n\t"\
		"movaps		0x70(%%eax)	,%%xmm5\n\t"\
		"unpckhpd	0x70(%%esi)	,%%xmm7\n\t"\
		"unpcklpd	0x70(%%esi)	,%%xmm5		/* a[jp+p6 ] */\n\t"\
		"movaps			%%xmm7	,0x150(%%ecx)	/* Store hi imag in t14+16 */\n\t"\
		"\n\t"\
		"movaps			%%xmm4	,%%xmm6	/* xmm6 <- cpy a[jt+p6] */\n\t"\
		"movaps			%%xmm5	,%%xmm7	/* xmm7 <- cpy a[jp+p6] */\n\t"\
		"\n\t"\
		"mulpd		    (%%edx)	,%%xmm4		/* a[jt+p6]*c6 */\n\t"\
		"mulpd		    (%%edx)	,%%xmm5		/* a[jp+p6]*c6 */\n\t"\
		"mulpd		0x10(%%edx)	,%%xmm6		/* a[jt+p6]*s6 */\n\t"\
		"mulpd		0x10(%%edx)	,%%xmm7		/* a[jp+p6]*s6 */\n\t"\
		"addpd		%%xmm6		,%%xmm5	/* xmm5 <- t14*/\n\t"\
		"subpd		%%xmm7		,%%xmm4	/* xmm4 <- t13*/\n\t"\
		"movaps		%%xmm5		,%%xmm7	/* xmm7 <- cpy t14*/\n\t"\
		"movaps		%%xmm4		,%%xmm6	/* xmm6 <- cpy t13*/\n\t"\
		"\n\t"\
		"subpd		0x060(%%ecx)	,%%xmm4		/* ~t15<- t13-rt */\n\t"\
		"subpd		0x070(%%ecx)	,%%xmm5		/* ~t16<- t14-it */\n\t"\
		"addpd		0x060(%%ecx)	,%%xmm6		/* ~t13<- t13+rt */\n\t"\
		"addpd		0x070(%%ecx)	,%%xmm7		/* ~t14<- t14+it */\n\t"\
		"\n\t"\
	"/* Finish radix-4 butterfly and store results into temporary-array slots: */\n\t"\
	"/*\n\t"\
	"~t13=t9 -t5;		~t9 =t9 +t5;\n\t"\
	"~t14=t10-t6;		~t10=t10+t6;\n\t"\
	"*/\n\t"\
		"subpd		%%xmm6		,%%xmm0	/*~t13*/\n\t"\
		"subpd		%%xmm7		,%%xmm1	/*~t14*/\n\t"\
		"movaps		%%xmm0		,0x040(%%ecx)	/* a[jt+p8 ] <- ~t13*/\n\t"\
		"movaps		%%xmm1		,0x050(%%ecx)	/* a[jp+p8 ] <- ~t14*/\n\t"\
		"addpd		%%xmm6		,%%xmm6	/* 2*t13*/\n\t"\
		"addpd		%%xmm7		,%%xmm7	/* 2*t14*/\n\t"\
		"addpd		%%xmm0		,%%xmm6	/*~t9 */\n\t"\
		"addpd		%%xmm1		,%%xmm7	/*~t10*/\n\t"\
		"movaps		%%xmm6		,     (%%ecx)	/* a[jt    ] <- ~t9 */\n\t"\
		"movaps		%%xmm7		,0x010(%%ecx)	/* a[jp    ] <- ~t10*/\n\t"\
		"\n\t"\
	"/*\n\t"\
	"~t15=t11+t8;		~t11=t11-t8;\n\t"\
	"~t16=t12-t7;		~t12=t12+t7;\n\t"\
	"*/\n\t"\
		"subpd		%%xmm5		,%%xmm2	/*~t11*/\n\t"\
		"subpd		%%xmm4		,%%xmm3	/*~t16*/\n\t"\
		"movaps		%%xmm2		,0x020(%%ecx)	/* a[jt+p4 ] <- ~t11*/\n\t"\
		"movaps		%%xmm3		,0x070(%%ecx)	/* a[jp+p12] <- ~t16*/\n\t"\
		"addpd		%%xmm5		,%%xmm5	/* 2*t16*/\n\t"\
		"addpd		%%xmm4		,%%xmm4	/* 2*t15*/\n\t"\
		"addpd		%%xmm2		,%%xmm5	/*~t15*/\n\t"\
		"addpd		%%xmm3		,%%xmm4	/*~t12*/\n\t"\
		"movaps		%%xmm5		,0x060(%%ecx)	/* a[jt+p12] <- ~t15*/\n\t"\
		"movaps		%%xmm4		,0x030(%%ecx)	/* a[jp+p4 ] <- ~t12*/\n\t"\
		"\n\t"\
"/******************************************************************************************************************************/\n\t"\
"/* Each of the above 2 radix-4 blocks is unique, but the ensuing 2 blocks                                                     */\n\t"\
"/* [operating on the odd-indexed elements from the unpck*pd commands which were stored to temporaries can use a common macro: */\n\t"\
"/******************************************************************************************************************************/\n\t"\
		"\n\t"\
"/*...Block 3: */\n\t"\
		"\n\t"\
"/*	SSE2_RADIX4_DIF_4TWIDDLE(r17,r21,r19,r23, r17,c1): */\n\t"\
	"/* Do the p0,p8 combo: */\n\t"\
		"movl		%[__r17],%%eax\n\t"\
		"movl		%[__c1] ,%%esi\n\t"\
		"movl		%%eax   ,%%ecx\n\t"\
		"addl		$0x20   ,%%ecx	/* r19 */\n\t"\
		"\n\t"\
	"prefetcht0	0x80(%%edi)\n\t"\
		"movaps		    (%%eax)	,%%xmm0		/* a[jt   ] */				\n\t	movaps	    (%%ecx),%%xmm4		/* a[jt+p8 ] */\n\t"\
		"movaps		0x10(%%eax)	,%%xmm1		/* a[jp   ] */				\n\t	movaps	0x10(%%ecx),%%xmm5		/* a[jp+p8 ] */\n\t"\
		"movaps		    (%%esi)	,%%xmm6		/* c0 */\n\t"\
		"movaps		0x10(%%esi)	,%%xmm7		/* s0 */\n\t"\
		"movaps		%%xmm0		,%%xmm2		/* xmm2 <- cpy a[jt   ] */\n\t"\
		"movaps		%%xmm1		,%%xmm3		/* xmm3 <- cpy a[jp   ] */\n\t"\
		"\n\t"\
		"mulpd   	%%xmm6		,%%xmm0		/* a[jt   ]*c0 */\n\t"\
		"mulpd   	%%xmm6		,%%xmm1		/* a[jp   ]*c0 */\n\t"\
		"mulpd   	%%xmm7		,%%xmm2		/* a[jt   ]*s0 */\n\t"\
		"mulpd   	%%xmm7		,%%xmm3		/* a[jp   ]*s0	xmm6,7 free */\n\t"\
		"																\n\t	movaps	%%xmm4		,%%xmm6		/* xmm6 <- cpy a[jt+p8 ] */\n\t"\
		"addpd   	%%xmm2		,%%xmm1		/* xmm1 <- t2 */			\n\t	movaps	%%xmm5		,%%xmm7		/* xmm7 <- cpy a[jp+p8 ] */\n\t"\
		"																\n\t	mulpd   0x20(%%esi)	,%%xmm4		/* a[jt+p8 ]*c8 */\n\t"\
		"subpd   	%%xmm3		,%%xmm0		/* xmm0 <- t1 */			\n\t	mulpd   0x20(%%esi)	,%%xmm5		/* a[jp+p8 ]*c8 */\n\t"\
		"																\n\t	mulpd   0x30(%%esi)	,%%xmm6		/* a[jt+p8 ]*s8 */\n\t"\
		"movaps		%%xmm0		,%%xmm2		/* xmm2 <- cpy t1 */		\n\t	mulpd	0x30(%%esi)	,%%xmm7		/* a[jp+p8 ]*s8 */\n\t"\
		"																\n\t	addpd   %%xmm6	    ,%%xmm5		/* xmm5 <- it */\n\t"\
		"movaps		%%xmm1		,%%xmm3		/* xmm3 <- cpy t2 */		\n\t	subpd	%%xmm7		,%%xmm4		/* xmm4 <- rt    xmm6,7 free */\n\t"\
		"\n\t"\
			"															\n\t	addl	$0x40 ,%%ecx	/* r23 */\n\t"\
		"																\n\t	addl	$0x60 ,%%esi\n\t"\
		"																\n\t	movaps	    (%%ecx),%%xmm6		/* a[jt+p12] */\n\t"\
		"																\n\t	movaps	0x10(%%ecx),%%xmm7		/* a[jp+p12] */\n\t"\
		"\n\t"\
		"addpd		%%xmm4		,%%xmm0		/* ~t1 <- t1 +rt */\n\t"\
		"addpd		%%xmm5		,%%xmm1		/* ~t2 <- t2 +it */\n\t"\
		"subpd		%%xmm4		,%%xmm2		/* ~t3 <- t1 -rt */\n\t"\
		"subpd		%%xmm5		,%%xmm3		/* ~t4 <- t2 -it	xmm4,5 free */\n\t"\
		"\n\t"\
	"/* Do the p4,12 combo: */\n\t"\
		"movaps		%%xmm6		,%%xmm4		/* xmm4 <- cpy a[jt+p12] */\n\t"\
		"movaps		%%xmm7		,%%xmm5		/* xmm5 <- cpy a[jp+p12] */\n\t"\
		"\n\t"\
		"mulpd		    (%%esi)	,%%xmm4		/* a[jt+p12]*c12 */\n\t"\
		"mulpd		    (%%esi)	,%%xmm5		/* a[jp+p12]*c12 */\n\t"\
		"mulpd		0x10(%%esi)	,%%xmm6		/* a[jt+p12]*s12 */\n\t"\
		"mulpd		0x10(%%esi)	,%%xmm7		/* a[jp+p12]*s12 */\n\t"\
		"movl		%%eax,%%edx	/* r17 */\n\t"\
		"addpd		%%xmm6		,%%xmm5	/* xmm5 <- it */\n\t"\
		"subpd		%%xmm7		,%%xmm4	/* xmm4 <- rt */\n\t"\
		"movaps		%%xmm5		,0x010(%%edx)	/* store it */\n\t"\
		"movaps		%%xmm4		,     (%%edx)	/* store rt */\n\t"\
		"\n\t"\
		"addl		$0x40 ,%%eax	/* r21 */\n\t"\
		"subl		$0x20 ,%%esi\n\t"\
		"movaps		    (%%eax)	,%%xmm4		/* a[jt+p4] */\n\t"\
		"movaps		0x10(%%eax)	,%%xmm5		/* a[jp+p4] */\n\t"\
		"movaps			%%xmm4	,%%xmm6		/* xmm4 <- cpy a[jt+p4] */\n\t"\
		"movaps			%%xmm5	,%%xmm7		/* xmm5 <- cpy a[jp+p4] */\n\t"\
		"\n\t"\
		"mulpd		    (%%esi)	,%%xmm4		/* a[jt+p4]*c4 */\n\t"\
		"mulpd		    (%%esi)	,%%xmm5		/* a[jp+p4]*c4 */\n\t"\
		"mulpd		0x10(%%esi)	,%%xmm6		/* a[jt+p4]*s4 */\n\t"\
		"mulpd		0x10(%%esi)	,%%xmm7		/* a[jp+p4]*s4 */\n\t"\
		"addpd		%%xmm6		,%%xmm5	/* xmm5 <- t6 */\n\t"\
		"subpd		%%xmm7		,%%xmm4	/* xmm4 <- t5 	xmm6,7 free */\n\t"\
		"movaps		%%xmm5		,%%xmm7	/* xmm7 <- cpy t6 */\n\t"\
		"movaps		%%xmm4		,%%xmm6	/* xmm6 <- cpy t5 */\n\t"\
		"\n\t"\
		"subpd		     (%%edx),%%xmm4		/* ~t7 <- t5 -rt */\n\t"\
		"subpd		0x010(%%edx),%%xmm5		/* ~t8 <- t6 -it */\n\t"\
		"addpd		     (%%edx),%%xmm6		/* ~t5 <- t5 +rt */\n\t"\
		"addpd		0x010(%%edx),%%xmm7		/* ~t6 <- t6 +it */\n\t"\
		"\n\t"\
	"/* Finish radix-4 butterfly and store results into temporary-array slots: */\n\t"\
		"subpd		%%xmm6		,%%xmm0	/*~t5 */						\n\t	subpd	%%xmm5,%%xmm2			/*~t3 */\n\t"\
		"subpd		%%xmm7		,%%xmm1	/*~t6 */						\n\t	subpd	%%xmm4,%%xmm3			/*~t8 */\n\t"\
		"movaps		%%xmm0		,0x040(%%edx)	/* a[jt+p8 ] <- ~t5 */	\n\t	movaps	%%xmm2,0x020(%%edx)		/* a[jt+p4 ] <- ~t3 */\n\t"\
		"movaps		%%xmm1		,0x050(%%edx)	/* a[jp+p8 ] <- ~t6 */	\n\t	movaps	%%xmm3,0x070(%%edx)		/* a[jp+p12] <- ~t8 */\n\t"\
		"addpd		%%xmm6		,%%xmm6	/* 2*t5 */						\n\t	addpd	%%xmm5,%%xmm5			/* 2*t8 */\n\t"\
		"addpd		%%xmm7		,%%xmm7	/* 2*t6 */						\n\t	addpd	%%xmm4,%%xmm4			/* 2*t7 */\n\t"\
		"addpd		%%xmm0		,%%xmm6	/*~t1 */						\n\t	addpd	%%xmm2,%%xmm5			/*~t7 */\n\t"\
		"addpd		%%xmm1		,%%xmm7	/*~t2 */						\n\t	addpd	%%xmm3,%%xmm4			/*~t4 */\n\t"\
		"movaps		%%xmm6		,     (%%edx)	/* a[jt    ] <- ~t1 */	\n\t	movaps	%%xmm5,0x060(%%edx)		/* a[jt+p12] <- ~t7 */\n\t"\
		"movaps		%%xmm7		,0x010(%%edx)	/* a[jp    ] <- ~t2 */	\n\t	movaps	%%xmm4,0x030(%%edx)		/* a[jp+p4 ] <- ~t4 */\n\t"\
		"\n\t"\
"/*...Block 4: */\n\t"\
		"\n\t"\
	"prefetcht0	0xa0(%%edi)\n\t"\
"/* SSE2_RADIX4_DIF_4TWIDDLE(r25,r29,r27,r31, r25,c3): */\n\t"\
	"/* Do the p0,p8 combo: */\n\t"\
		"movl		%[__r25],%%eax\n\t"\
		"movl		%[__c3] ,%%esi\n\t"\
		"movl		%%eax   ,%%ecx\n\t"\
		"addl		$0x20   ,%%ecx	/* r27 */\n\t"\
		"\n\t"\
		"movaps		    (%%eax)	,%%xmm0		/* a[jt   ] */				\n\t	movaps	    (%%ecx),%%xmm4		/* a[jt+p8 ] */\n\t"\
		"movaps		0x10(%%eax)	,%%xmm1		/* a[jp   ] */				\n\t	movaps	0x10(%%ecx),%%xmm5		/* a[jp+p8 ] */\n\t"\
		"movaps		    (%%esi)	,%%xmm6		/* c0 */\n\t"\
		"movaps		0x10(%%esi)	,%%xmm7		/* s0 */\n\t"\
		"movaps		%%xmm0		,%%xmm2		/* xmm2 <- cpy a[jt   ] */\n\t"\
		"movaps		%%xmm1		,%%xmm3		/* xmm3 <- cpy a[jp   ] */\n\t"\
		"\n\t"\
		"mulpd   	%%xmm6		,%%xmm0		/* a[jt   ]*c0 */\n\t"\
		"mulpd   	%%xmm6		,%%xmm1		/* a[jp   ]*c0 */\n\t"\
		"mulpd   	%%xmm7		,%%xmm2		/* a[jt   ]*s0 */\n\t"\
		"mulpd   	%%xmm7		,%%xmm3		/* a[jp   ]*s0	xmm6,7 free */\n\t"\
		"																\n\t	movaps	%%xmm4		,%%xmm6		/* xmm6 <- cpy a[jt+p8 ] */\n\t"\
		"addpd   	%%xmm2		,%%xmm1		/* xmm1 <- t2 */			\n\t	movaps	%%xmm5		,%%xmm7		/* xmm7 <- cpy a[jp+p8 ] */\n\t"\
		"																\n\t	mulpd   0x20(%%esi)	,%%xmm4		/* a[jt+p8 ]*c8 */\n\t"\
		"subpd   	%%xmm3		,%%xmm0		/* xmm0 <- t1 */			\n\t	mulpd   0x20(%%esi)	,%%xmm5		/* a[jp+p8 ]*c8 */\n\t"\
		"																\n\t	mulpd   0x30(%%esi)	,%%xmm6		/* a[jt+p8 ]*s8 */\n\t"\
		"movaps		%%xmm0		,%%xmm2		/* xmm2 <- cpy t1 */		\n\t	mulpd	0x30(%%esi)	,%%xmm7		/* a[jp+p8 ]*s8 */\n\t"\
		"																\n\t	addpd   %%xmm6	    ,%%xmm5		/* xmm5 <- it */\n\t"\
		"movaps		%%xmm1		,%%xmm3		/* xmm3 <- cpy t2 */		\n\t	subpd	%%xmm7		,%%xmm4		/* xmm4 <- rt    xmm6,7 free */\n\t"\
		"\n\t"\
			"															\n\t	addl	$0x40 ,%%ecx	/* r31 */\n\t"\
		"																\n\t	addl	$0x60 ,%%esi\n\t"\
		"																\n\t	movaps	    (%%ecx),%%xmm6		/* a[jt+p12] */\n\t"\
		"																\n\t	movaps	0x10(%%ecx),%%xmm7		/* a[jp+p12] */\n\t"\
		"\n\t"\
		"addpd		%%xmm4		,%%xmm0		/* ~t1 <- t1 +rt */\n\t"\
		"addpd		%%xmm5		,%%xmm1		/* ~t2 <- t2 +it */\n\t"\
		"subpd		%%xmm4		,%%xmm2		/* ~t3 <- t1 -rt */\n\t"\
		"subpd		%%xmm5		,%%xmm3		/* ~t4 <- t2 -it	xmm4,5 free */\n\t"\
		"\n\t"\
	"/* Do the p4,12 combo: */\n\t"\
		"movaps		%%xmm6		,%%xmm4		/* xmm4 <- cpy a[jt+p12] */\n\t"\
		"movaps		%%xmm7		,%%xmm5		/* xmm5 <- cpy a[jp+p12] */\n\t"\
		"\n\t"\
		"mulpd		    (%%esi)	,%%xmm4		/* a[jt+p12]*c12 */\n\t"\
		"mulpd		    (%%esi)	,%%xmm5		/* a[jp+p12]*c12 */\n\t"\
		"mulpd		0x10(%%esi)	,%%xmm6		/* a[jt+p12]*s12 */\n\t"\
		"mulpd		0x10(%%esi)	,%%xmm7		/* a[jp+p12]*s12 */\n\t"\
		"movl		%%eax,%%edx	/* r25 */\n\t"\
		"addpd		%%xmm6		,%%xmm5	/* xmm5 <- it */\n\t"\
		"subpd		%%xmm7		,%%xmm4	/* xmm4 <- rt */\n\t"\
		"movaps		%%xmm5		,0x010(%%edx)	/* store it */\n\t"\
		"movaps		%%xmm4		,     (%%edx)	/* store rt */\n\t"\
		"\n\t"\
		"addl		$0x40 ,%%eax	/* r29 */\n\t"\
		"subl		$0x20 ,%%esi\n\t"\
		"movaps		    (%%eax)	,%%xmm4		/* a[jt+p4] */\n\t"\
		"movaps		0x10(%%eax)	,%%xmm5		/* a[jp+p4] */\n\t"\
		"movaps			%%xmm4	,%%xmm6		/* xmm4 <- cpy a[jt+p4] */\n\t"\
		"movaps			%%xmm5	,%%xmm7		/* xmm5 <- cpy a[jp+p4] */\n\t"\
		"\n\t"\
		"mulpd		    (%%esi)	,%%xmm4		/* a[jt+p4]*c4 */\n\t"\
		"mulpd		    (%%esi)	,%%xmm5		/* a[jp+p4]*c4 */\n\t"\
		"mulpd		0x10(%%esi)	,%%xmm6		/* a[jt+p4]*s4 */\n\t"\
		"mulpd		0x10(%%esi)	,%%xmm7		/* a[jp+p4]*s4 */\n\t"\
		"addpd		%%xmm6		,%%xmm5	/* xmm5 <- t6 */\n\t"\
		"subpd		%%xmm7		,%%xmm4	/* xmm4 <- t5 	xmm6,7 free */\n\t"\
		"movaps		%%xmm5		,%%xmm7	/* xmm7 <- cpy t6 */\n\t"\
		"movaps		%%xmm4		,%%xmm6	/* xmm6 <- cpy t5 */\n\t"\
		"\n\t"\
		"subpd		     (%%edx),%%xmm4		/* ~t7 <- t5 -rt */\n\t"\
		"subpd		0x010(%%edx),%%xmm5		/* ~t8 <- t6 -it */\n\t"\
		"addpd		     (%%edx),%%xmm6		/* ~t5 <- t5 +rt */\n\t"\
		"addpd		0x010(%%edx),%%xmm7		/* ~t6 <- t6 +it */\n\t"\
		"\n\t"\
	"/* Finish radix-4 butterfly and store results into temporary-array slots: */\n\t"\
	"prefetcht0	0xc0(%%edi)\n\t"\
		"subpd		%%xmm6		,%%xmm0	/*~t5 */						\n\t	subpd	%%xmm5,%%xmm2			/*~t3 */\n\t"\
		"subpd		%%xmm7		,%%xmm1	/*~t6 */						\n\t	subpd	%%xmm4,%%xmm3			/*~t8 */\n\t"\
		"movaps		%%xmm0		,0x040(%%edx)	/* a[jt+p8 ] <- ~t5 */	\n\t	movaps	%%xmm2,0x020(%%edx)		/* a[jt+p4 ] <- ~t3 */\n\t"\
		"movaps		%%xmm1		,0x050(%%edx)	/* a[jp+p8 ] <- ~t6 */	\n\t	movaps	%%xmm3,0x070(%%edx)		/* a[jp+p12] <- ~t8 */\n\t"\
		"addpd		%%xmm6		,%%xmm6	/* 2*t5 */						\n\t	addpd	%%xmm5,%%xmm5			/* 2*t8 */\n\t"\
		"addpd		%%xmm7		,%%xmm7	/* 2*t6 */						\n\t	addpd	%%xmm4,%%xmm4			/* 2*t7 */\n\t"\
		"addpd		%%xmm0		,%%xmm6	/*~t1 */						\n\t	addpd	%%xmm2,%%xmm5			/*~t7 */\n\t"\
		"addpd		%%xmm1		,%%xmm7	/*~t2 */						\n\t	addpd	%%xmm3,%%xmm4			/*~t4 */\n\t"\
		"movaps		%%xmm6		,     (%%edx)	/* a[jt    ] <- ~t1 */	\n\t	movaps	%%xmm5,0x060(%%edx)		/* a[jt+p12] <- ~t7 */\n\t"\
		"movaps		%%xmm7		,0x010(%%edx)	/* a[jp    ] <- ~t2 */	\n\t	movaps	%%xmm4,0x030(%%edx)		/* a[jp+p4 ] <- ~t4 */\n\t"\
		"\n\t"\
/**************************************************************************************/\
/*...and now do four more radix-4 transforms, including the internal twiddle factors: */\
/**************************************************************************************/\
		"\n\t"\
"/*...Block 1: t1,9,17,25 */\n\t"\
		"movl		%[__r1] ,%%eax\n\t"\
		"movl		%[__r9] ,%%esi\n\t"\
		"movl		%[__r17],%%ecx\n\t"\
		"movl		%[__r25],%%edx\n\t"\
		"\n\t"\
		"movaps		     (%%eax)	,%%xmm0		/* t1  */\n\t"\
		"movaps		0x010(%%eax)	,%%xmm1		/* t2  */\n\t"\
		"movaps		     (%%esi)	,%%xmm2		/* t9  */\n\t"\
		"movaps		0x010(%%esi)	,%%xmm3		/* t10 */\n\t"\
		"\n\t"\
		"subpd		     (%%esi)	,%%xmm0		/* t9 =t1 -rt */\n\t"\
		"subpd		0x010(%%esi)	,%%xmm1		/* t10=t2 -it */\n\t"\
		"addpd		     (%%eax)	,%%xmm2		/* t1 =t1 +rt */\n\t"\
		"addpd		0x010(%%eax)	,%%xmm3		/* t2 =t2 +it */\n\t"\
		"\n\t"\
		"movaps		     (%%ecx)	,%%xmm4		/* t17 */\n\t"\
		"movaps		0x010(%%ecx)	,%%xmm5		/* t18 */\n\t"\
		"movaps		     (%%edx)	,%%xmm6		/* t25 */\n\t"\
		"movaps		0x010(%%edx)	,%%xmm7		/* t26 */\n\t"\
		"\n\t"\
		"subpd		     (%%edx)	,%%xmm4		/* t25=t17-rt */\n\t"\
		"subpd		0x010(%%edx)	,%%xmm5		/* t26=t18-it */\n\t"\
		"addpd		     (%%ecx)	,%%xmm6		/* t17=t17+rt */\n\t"\
		"addpd		0x010(%%ecx)	,%%xmm7		/* t18=t18+it */\n\t"\
		"\n\t"\
		"subpd		%%xmm6		,%%xmm2		/* t1  <- t1 -t17 */\n\t"\
		"subpd		%%xmm7		,%%xmm3		/* t2  <- t2 -t18 */\n\t"\
		"addpd		%%xmm6		,%%xmm6		/*          2*t17 */\n\t"\
		"addpd		%%xmm7		,%%xmm7		/*          2*t18 */\n\t"\
		"movaps		%%xmm2		,    (%%ecx)/* a[jt+p1 ], store in t17 */\n\t"\
		"movaps		%%xmm3		,0x10(%%ecx)/* a[jp+p1 ], store in t18 */\n\t"\
		"addpd		%%xmm2		,%%xmm6		/* t17 <- t1 +t17 */\n\t"\
		"addpd		%%xmm3		,%%xmm7		/* t18 <- t2 +t18 */\n\t"\
		"movaps		%%xmm6		,    (%%eax)/* a[jt+p0 ], store in t0  */\n\t"\
		"movaps		%%xmm7		,0x10(%%eax)/* a[jp+p0 ], store in t1  */\n\t"\
		"\n\t"\
		"subpd		%%xmm5		,%%xmm0		/* t9  <- t9 -t26 */\n\t"\
		"subpd		%%xmm4		,%%xmm1		/* t10 <- t10-t25 */\n\t"\
		"addpd		%%xmm5		,%%xmm5		/*          2*t26 */\n\t"\
		"addpd		%%xmm4		,%%xmm4		/*          2*t25 */\n\t"\
		"movaps		%%xmm0		,    (%%esi)/* a[jt+p2 ], store in t9  */\n\t"\
		"movaps		%%xmm1		,0x10(%%edx)/* a[jp+p3 ], store in t26 */\n\t"\
		"addpd		%%xmm0		,%%xmm5		/* t26 <- t9 +t26 */\n\t"\
		"addpd		%%xmm1		,%%xmm4		/* t25 <- t10+t25 */\n\t"\
		"movaps		%%xmm5		,    (%%edx)/* a[jt+p3 ], store in t25 */\n\t"\
		"movaps		%%xmm4		,0x10(%%esi)/* a[jp+p2 ], store in t10 */\n\t"\
		"\n\t"\
"/*...Block 3: t5,13,21,29 */\n\t"\
		"addl		$0x40,%%eax	/* r5  */	\n\t"\
		"addl		$0x40,%%esi	/* r13 */	\n\t"\
		"addl		$0x40,%%ecx	/* r21 */	\n\t"\
		"addl		$0x40,%%edx	/* r29 */	\n\t"\
		"\n\t"\
		"movaps		     (%%eax)	,%%xmm0		/* t5  */\n\t"\
		"movaps		0x010(%%eax)	,%%xmm1		/* t6  */\n\t"\
		"movaps		0x080(%%eax)	,%%xmm2		/* t13 */\n\t"\
		"movaps		0x090(%%eax)	,%%xmm3		/* t14 */\n\t"\
		"\n\t"\
		"subpd		0x090(%%eax)	,%%xmm0		/* t5 =t5 -t14*/\n\t"\
		"subpd		0x080(%%eax)	,%%xmm1		/* t14=t6 -t13*/\n\t"\
		"addpd		0x010(%%eax)	,%%xmm2		/* t6 =t13+t6 */\n\t"\
		"addpd		     (%%eax)	,%%xmm3		/* t13=t14+t5 */\n\t"\
	"prefetcht0	0xe0(%%edi)\n\t"\
		"movl		%[__isrt2],%%edi\n\t"\
		"\n\t"\
		"movaps		0x100(%%eax)	,%%xmm4		/* t21 */\n\t"\
		"movaps		0x110(%%eax)	,%%xmm5		/* t22 */\n\t"\
		"movaps		0x180(%%eax)	,%%xmm6		/* t29 */\n\t"\
		"movaps		0x190(%%eax)	,%%xmm7		/* t30 */\n\t"\
		"\n\t"\
		"subpd		0x110(%%eax)	,%%xmm4		/* t21-t22 */\n\t"\
		"addpd		0x100(%%eax)	,%%xmm5		/* t22+t21 */\n\t"\
		"mulpd		(%%edi)		,%%xmm4	/* t21 = (t21-t22)*ISRT2 */\n\t"\
		"mulpd		(%%edi)		,%%xmm5	/* t22 = (t22+t21)*ISRT2 */\n\t"\
		"\n\t"\
		"addpd		0x190(%%eax)	,%%xmm6		/* t29+t30 */\n\t"\
		"subpd		0x180(%%eax)	,%%xmm7		/* t30-t29 */\n\t"\
		"mulpd		(%%edi)		,%%xmm6	/*  rt = (t29+t30)*ISRT2 */\n\t"\
		"mulpd		(%%edi)		,%%xmm7	/*  it = (t30-t29)*ISRT2 */\n\t"\
		"\n\t"\
		"subpd		%%xmm6		,%%xmm4		/* t21=t21-rt */\n\t"\
		"subpd		%%xmm7		,%%xmm5		/* t22=t22-it */\n\t"\
		"addpd		%%xmm6		,%%xmm6		/*      2* rt */\n\t"\
		"addpd		%%xmm7		,%%xmm7		/*      2* it */\n\t"\
		"addpd		%%xmm4		,%%xmm6		/* t29=t21+rt */\n\t"\
		"addpd		%%xmm5		,%%xmm7		/* t30=t22+it */\n\t"\
		"\n\t"\
		"subpd		%%xmm4		,%%xmm0		/* t5 -t21 */\n\t"\
		"subpd		%%xmm5		,%%xmm2		/* t6 -t22 */\n\t"\
		"addpd		%%xmm4		,%%xmm4		/*   2*t21 */\n\t"\
		"addpd		%%xmm5		,%%xmm5		/*   2*t22 */\n\t"\
		"\n\t"\
		"movaps		%%xmm0		,    (%%ecx)/* a[jt+p1 ] */\n\t"\
		"movaps		%%xmm2		,0x10(%%ecx)/* a[jp+p1 ] */\n\t"\
		"addpd		%%xmm0		,%%xmm4		/* t5 +t21 */\n\t"\
		"addpd		%%xmm2		,%%xmm5		/* t6 +t22 */\n\t"\
		"movaps		%%xmm4		,    (%%eax)/* a[jt+p0 ] */\n\t"\
		"movaps		%%xmm5		,0x10(%%eax)/* a[jp+p0 ] */\n\t"\
		"\n\t"\
		"subpd		%%xmm7		,%%xmm3		/* t13-t30 */\n\t"\
		"subpd		%%xmm6		,%%xmm1		/* t14-t29 */\n\t"\
		"addpd		%%xmm7		,%%xmm7		/*   2*t30 */\n\t"\
		"addpd		%%xmm6		,%%xmm6		/*   2*t29 */\n\t"\
		"movaps		%%xmm3		,    (%%esi)/* a[jt+p2 ] */\n\t"\
		"movaps		%%xmm1		,0x10(%%edx)/* a[jp+p3 ] */\n\t"\
		"addpd		%%xmm3		,%%xmm7		/* t13+t30 */\n\t"\
		"addpd		%%xmm1		,%%xmm6		/* t14+t29 */\n\t"\
		"movaps		%%xmm7		,    (%%edx)/* a[jt+p3 ] */\n\t"\
		"movaps		%%xmm6		,0x10(%%esi)/* a[jp+p2 ] */\n\t"\
		"\n\t"\
"/*...Block 2: t3,11,19,27 */\n\t"\
		"subl		$0x20,%%eax	/* r3  */	\n\t"\
		"subl		$0x20,%%esi	/* r11 */	\n\t"\
		"subl		$0x20,%%ecx	/* r19 */	\n\t"\
		"subl		$0x20,%%edx	/* r27 */	\n\t"\
		"movl		%[__cc0],%%edi\n\t"\
		"\n\t"\
		"movaps		0x100(%%eax),%%xmm4		/* t19 */			\n\t	movaps	0x180(%%eax),%%xmm6	/* t27 */\n\t"\
		"movaps		0x110(%%eax),%%xmm5		/* t20 */			\n\t	movaps	0x190(%%eax),%%xmm7	/* t28 */\n\t"\
		"movaps		0x100(%%eax),%%xmm0		/* copy t19 */		\n\t	movaps	0x180(%%eax),%%xmm2	/* copy t27 */\n\t"\
		"movaps		0x110(%%eax),%%xmm1		/* copy t20 */		\n\t	movaps	0x190(%%eax),%%xmm3	/* copy t28 */\n\t"\
		"\n\t"\
		"mulpd		    (%%edi)	,%%xmm4		/* t19*c */			\n\t	mulpd	0x10(%%edi)	,%%xmm6	/* t27*s */\n\t"\
		"mulpd		0x10(%%edi)	,%%xmm1		/* t20*s */			\n\t	mulpd	    (%%edi)	,%%xmm3	/* t28*c */\n\t"\
		"mulpd		    (%%edi)	,%%xmm5		/* t20*c */			\n\t	mulpd	0x10(%%edi)	,%%xmm7	/* t28*s */\n\t"\
		"mulpd		0x10(%%edi)	,%%xmm0		/* t19*s */			\n\t	mulpd	    (%%edi)	,%%xmm2	/* t27*c */\n\t"\
		"subpd		%%xmm1		,%%xmm4	/* ~t19 */				\n\t	subpd	%%xmm3,%%xmm6	/* rt */\n\t"\
		"addpd		%%xmm0		,%%xmm5	/* ~t20 */				\n\t	addpd	%%xmm2,%%xmm7	/* it */\n\t"\
		"\n\t"\
		"subpd		%%xmm6		,%%xmm4		/*~t27=t19-rt */\n\t"\
		"subpd		%%xmm7		,%%xmm5		/*~t28=t20-it */\n\t"\
		"addpd		%%xmm6		,%%xmm6		/*      2* rt */\n\t"\
		"addpd		%%xmm7		,%%xmm7		/*      2* it */\n\t"\
		"addpd		%%xmm4		,%%xmm6		/*~t19=t19+rt */\n\t"\
		"addpd		%%xmm5		,%%xmm7		/*~t20=t20+it */\n\t"\
		"\n\t"\
		"movaps		0x080(%%eax)	,%%xmm2		/* t11 */\n\t"\
		"movaps		0x090(%%eax)	,%%xmm3		/* t12 */\n\t"\
		"subpd		0x090(%%eax)	,%%xmm2		/* t11-t12 */\n\t"\
		"addpd		0x080(%%eax)	,%%xmm3		/* t12+t11 */\n\t"\
		"movl		%[__isrt2],%%edi\n\t"\
		"mulpd		(%%edi)		,%%xmm2	/* rt = (t11-t12)*ISRT2 */\n\t"\
		"mulpd		(%%edi)		,%%xmm3	/* it = (t12+t11)*ISRT2 */\n\t"\
		"\n\t"\
		"movaps		     (%%eax)	,%%xmm0		/* t3  */\n\t"\
		"movaps		0x010(%%eax)	,%%xmm1		/* t4  */\n\t"\
		"\n\t"\
		"subpd		%%xmm2,%%xmm0					/*~t11=t3 -rt */\n\t"\
		"subpd		%%xmm3,%%xmm1					/*~t12=t4 -it */\n\t"\
		"addpd		     (%%eax)	,%%xmm2		/*~t3 =rt +t3 */\n\t"\
		"addpd		0x010(%%eax)	,%%xmm3		/*~t4 =it +t4 */\n\t"\
		"\n\t"\
		"\n\t"\
		"subpd		%%xmm6		,%%xmm2		/* t3 -t19 */\n\t"\
		"subpd		%%xmm7		,%%xmm3		/* t4 -t20 */\n\t"\
		"addpd		%%xmm6		,%%xmm6		/*   2*t19 */\n\t"\
		"addpd		%%xmm7		,%%xmm7		/*   2*t20 */\n\t"\
		"movaps		%%xmm2		,    (%%ecx)/* a[jt+p1 ] */\n\t"\
		"movaps		%%xmm3		,0x10(%%ecx)/* a[jp+p1 ] */\n\t"\
		"addpd		%%xmm2		,%%xmm6		/* t3 +t19 */\n\t"\
		"addpd		%%xmm3		,%%xmm7		/* t4 +t20 */\n\t"\
		"movaps		%%xmm6		,    (%%eax)/* a[jt+p0 ] */\n\t"\
		"movaps		%%xmm7		,0x10(%%eax)/* a[jp+p0 ] */\n\t"\
		"\n\t"\
		"subpd		%%xmm5		,%%xmm0		/* t11-t28 */\n\t"\
		"subpd		%%xmm4		,%%xmm1		/* t12-t27 */\n\t"\
		"addpd		%%xmm5		,%%xmm5		/*          2*t28 */\n\t"\
		"addpd		%%xmm4		,%%xmm4		/*          2*t27 */\n\t"\
		"movaps		%%xmm0		,    (%%esi)/* a[jt+p2 ] */\n\t"\
		"movaps		%%xmm1		,0x10(%%edx)/* a[jp+p3 ] */\n\t"\
		"addpd		%%xmm0		,		%%xmm5/* t11+t28 */\n\t"\
		"addpd		%%xmm1		,		%%xmm4/* t12+t27 */\n\t"\
		"movaps		%%xmm5		,    (%%edx)/* a[jt+p3 ] */\n\t"\
		"movaps		%%xmm4		,0x10(%%esi)/* a[jp+p2 ] */\n\t"\
		"\n\t"\
"/*...Block 4: t7,15,23,31 */\n\t"\
		"addl		$0x40,%%eax	/* r7  */	\n\t"\
		"addl		$0x40,%%esi	/* r15 */	\n\t"\
		"addl		$0x40,%%ecx	/* r23 */	\n\t"\
		"addl		$0x40,%%edx	/* r31 */	\n\t"\
		"movl		%[__cc0],%%edi\n\t"\
		"\n\t"\
		"movaps		0x100(%%eax),%%xmm4		/* t23 */			\n\t	movaps	0x180(%%eax),%%xmm6		/* t31 */\n\t"\
		"movaps		0x110(%%eax),%%xmm5		/* t24 */			\n\t	movaps	0x190(%%eax),%%xmm7		/* t32 */\n\t"\
		"movaps		0x100(%%eax),%%xmm0		/* copy t23 */		\n\t	movaps	0x180(%%eax),%%xmm2		/* copy t31 */\n\t"\
		"movaps		0x110(%%eax),%%xmm1		/* copy t24 */		\n\t	movaps	0x190(%%eax),%%xmm3		/* copy t32 */\n\t"\
		"\n\t"\
		"mulpd		0x10(%%edi)	,%%xmm4		/* t23*s */			\n\t	mulpd	    (%%edi)	,%%xmm6		/* t31*c */\n\t"\
		"mulpd		    (%%edi)	,%%xmm1		/* t24*c */			\n\t	mulpd	0x10(%%edi)	,%%xmm3		/* t32*s */\n\t"\
		"mulpd		0x10(%%edi)	,%%xmm5		/* t24*s */			\n\t	mulpd	    (%%edi)	,%%xmm7		/* t32*c */\n\t"\
		"mulpd		    (%%edi)	,%%xmm0		/* t23*c */			\n\t	mulpd	0x10(%%edi)	,%%xmm2		/* t31*s */\n\t"\
		"subpd		%%xmm1		,%%xmm4	/* ~t23 */				\n\t	subpd	%%xmm3		,%%xmm6		/* rt */\n\t"\
		"addpd		%%xmm0		,%%xmm5	/* ~t24 */				\n\t	addpd	%%xmm2		,%%xmm7		/* it */\n\t"\
		"\n\t"\
		"subpd		%%xmm6		,%%xmm4		/*~t23=t23-rt */\n\t"\
		"subpd		%%xmm7		,%%xmm5		/*~t24=t24-it */\n\t"\
		"addpd		%%xmm6		,%%xmm6		/*      2* rt */\n\t"\
		"addpd		%%xmm7		,%%xmm7		/*      2* it */\n\t"\
		"addpd		%%xmm4		,%%xmm6		/*~t31=t23+rt */\n\t"\
		"addpd		%%xmm5		,%%xmm7		/*~t32=t24+it */\n\t"\
		"\n\t"\
		"movaps		0x080(%%eax)	,%%xmm2		/* t15 */\n\t"\
		"movaps		0x090(%%eax)	,%%xmm3		/* t16 */\n\t"\
		"addpd		0x090(%%eax)	,%%xmm2		/* t15+t16 */\n\t"\
		"subpd		0x080(%%eax)	,%%xmm3		/* t16-t15 */\n\t"\
		"movl		%[__isrt2],%%edi\n\t"\
		"mulpd		(%%edi)		,%%xmm2	/* rt = (t15+t16)*ISRT2 */\n\t"\
		"mulpd		(%%edi)		,%%xmm3	/* it = (t16-t15)*ISRT2 */\n\t"\
		"\n\t"\
		"movaps		     (%%eax)	,%%xmm0		/* t7  */\n\t"\
		"movaps		0x010(%%eax)	,%%xmm1		/* t8  */\n\t"\
		"\n\t"\
		"subpd		%%xmm2,%%xmm0					/*~t7 =t7 -rt */\n\t"\
		"subpd		%%xmm3,%%xmm1					/*~t8 =t8 -it */\n\t"\
		"addpd		     (%%eax)	,%%xmm2		/*~t15=rt +t7 */\n\t"\
		"addpd		0x010(%%eax)	,%%xmm3		/*~t16=it +t8 */\n\t"\
		"\n\t"\
		"subpd		%%xmm4		,%%xmm0		/* t7 -t23 */\n\t"\
		"subpd		%%xmm5		,%%xmm1		/* t8 -t24 */\n\t"\
		"addpd		%%xmm4		,%%xmm4		/*   2*t23 */\n\t"\
		"addpd		%%xmm5		,%%xmm5		/*   2*t24 */\n\t"\
		"movaps		%%xmm0		,    (%%ecx)/* a[jt+p1 ] */\n\t"\
		"movaps		%%xmm1		,0x10(%%ecx)/* a[jp+p1 ] */\n\t"\
		"addpd		%%xmm0		,%%xmm4		/* t7 +t23 */\n\t"\
		"addpd		%%xmm1		,%%xmm5		/* t8 +t24 */\n\t"\
		"movaps		%%xmm4		,    (%%eax)/* a[jt+p0 ] */\n\t"\
		"movaps		%%xmm5		,0x10(%%eax)/* a[jp+p0 ] */\n\t"\
		"\n\t"\
		"subpd		%%xmm7		,%%xmm2		/* t15-t32 */\n\t"\
		"subpd		%%xmm6		,%%xmm3		/* t16-t31 */\n\t"\
		"addpd		%%xmm7		,%%xmm7		/*   2*t32 */\n\t"\
		"addpd		%%xmm6		,%%xmm6		/*   2*t31 */\n\t"\
		"movaps		%%xmm2		,    (%%esi)/* a[jt+p2 ] */\n\t"\
		"movaps		%%xmm3		,0x10(%%edx)/* a[jp+p3 ] */\n\t"\
		"addpd		%%xmm2		,%%xmm7		/* t15+t32 */\n\t"\
		"addpd		%%xmm3		,%%xmm6		/* t16+t31 */\n\t"\
		"movaps		%%xmm7		,    (%%edx)/* a[jt+p3 ] */\n\t"\
		"movaps		%%xmm6		,0x10(%%esi)/* a[jp+p2 ] */\n\t"\
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
		: "cc","memory","eax","esi","ecx","edx","edi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"	/* Clobbered registers */\
	);\
	}

	// Since the add0/1 block addresses advance in opposite directions, our prefetch scheme here is like so:
	// In fwd DIF prefetch  ahead 8 cache lines' worth of data w.r.to add0;
	// In inv DIT prefetch behind 8 cache lines' worth of data w.r.to add1;
	#define SSE2_RADIX16_WRAPPER_DIT(Xadd0,Xadd1,Xr1,Xr9,Xr17,Xr25,Xisrt2,Xcc0,Xc1,Xc2,Xc3,Xc4,Xc5,Xc6,Xc7,Xc8,Xc9,Xc10,Xc11,Xc12,Xc13,Xc14,Xc15)\
	{\
	__asm__ volatile (\
		"/*...Block 1: */\n\t"\
		"movl	%[__r1],%%eax\n\t"\
		"movl	%%eax,%%esi\n\t"\
		"movl	%%eax,%%ecx\n\t"\
		"movl	%%eax,%%edx\n\t"\
		"addl	$0x100,%%esi\n\t"\
		"addl	$0x080,%%ecx\n\t"\
		"addl	$0x180,%%edx\n\t"\
		"/* SSE2_RADIX4_DIT_IN_PLACE(): */\n\t"\
		"movaps	    (%%eax),%%xmm0\n\t"\
		"movaps	0x10(%%eax),%%xmm1\n\t"\
		"movaps	    (%%eax),%%xmm2\n\t"\
		"movaps	0x10(%%eax),%%xmm3\n\t"\
		"\n\t"\
		"addpd	    (%%esi),%%xmm0\n\t"\
		"addpd	0x10(%%esi),%%xmm1\n\t"\
		"subpd	    (%%esi),%%xmm2\n\t"\
		"subpd	0x10(%%esi),%%xmm3\n\t"\
		"\n\t"\
		"movaps	    (%%ecx),%%xmm4\n\t"\
		"movaps	0x10(%%ecx),%%xmm5\n\t"\
		"movaps	    (%%ecx),%%xmm6\n\t"\
		"movaps	0x10(%%ecx),%%xmm7\n\t"\
		"\n\t"\
		"addpd	    (%%edx),%%xmm4\n\t"\
		"addpd	0x10(%%edx),%%xmm5\n\t"\
		"subpd	    (%%edx),%%xmm6\n\t"\
		"subpd	0x10(%%edx),%%xmm7\n\t"\
		"\n\t"\
		"subpd	%%xmm4,%%xmm0\n\t"\
		"subpd	%%xmm5,%%xmm1\n\t"\
		"movaps	%%xmm0,     (%%esi)\n\t"\
		"movaps	%%xmm1,0x010(%%esi)\n\t"\
		"addpd	%%xmm4,%%xmm4\n\t"\
		"addpd	%%xmm5,%%xmm5\n\t"\
		"addpd	%%xmm0,%%xmm4\n\t"\
		"addpd	%%xmm1,%%xmm5\n\t"\
		"movaps	%%xmm4,     (%%eax)\n\t"\
		"movaps	%%xmm5,0x010(%%eax)\n\t"\
		"\n\t"\
		"subpd	%%xmm7,%%xmm2\n\t"\
		"subpd	%%xmm6,%%xmm3\n\t"\
		"movaps	%%xmm2,     (%%edx)\n\t"\
		"movaps	%%xmm3,0x010(%%ecx)\n\t"\
		"addpd	%%xmm7,%%xmm7\n\t"\
		"addpd	%%xmm6,%%xmm6\n\t"\
		"addpd	%%xmm2,%%xmm7\n\t"\
		"addpd	%%xmm3,%%xmm6\n\t"\
		"movaps	%%xmm7,     (%%ecx)\n\t"\
		"movaps	%%xmm6,0x010(%%edx)\n\t"\
		"\n\t"\
		"/*...Block 2: */\n\t"\
		"addl	$0x040,%%eax\n\t"\
		"addl	$0x040,%%esi\n\t"\
		"addl	$0x040,%%ecx\n\t"\
		"addl	$0x040,%%edx\n\t"\
		"/* SSE2_RADIX4_DIT_IN_PLACE(): */\n\t"\
		"movaps	    (%%eax),%%xmm0\n\t"\
		"movaps	0x10(%%eax),%%xmm1\n\t"\
		"movaps	    (%%eax),%%xmm2\n\t"\
		"movaps	0x10(%%eax),%%xmm3\n\t"\
		"\n\t"\
		"addpd	    (%%esi),%%xmm0\n\t"\
		"addpd	0x10(%%esi),%%xmm1\n\t"\
		"subpd	    (%%esi),%%xmm2\n\t"\
		"subpd	0x10(%%esi),%%xmm3\n\t"\
		"\n\t"\
		"movaps	    (%%ecx),%%xmm4\n\t"\
		"movaps	0x10(%%ecx),%%xmm5\n\t"\
		"movaps	    (%%ecx),%%xmm6\n\t"\
		"movaps	0x10(%%ecx),%%xmm7\n\t"\
		"\n\t"\
		"addpd	    (%%edx),%%xmm4\n\t"\
		"addpd	0x10(%%edx),%%xmm5\n\t"\
		"subpd	    (%%edx),%%xmm6\n\t"\
		"subpd	0x10(%%edx),%%xmm7\n\t"\
	"movl	%[__add1],%%edi	\n\t"\
	"prefetcht0	-0x20(%%edi)\n\t"\
		"\n\t"\
		"subpd	%%xmm4,%%xmm0\n\t"\
		"subpd	%%xmm5,%%xmm1\n\t"\
		"movaps	%%xmm0,     (%%esi)\n\t"\
		"movaps	%%xmm1,0x010(%%esi)\n\t"\
		"addpd	%%xmm4,%%xmm4\n\t"\
		"addpd	%%xmm5,%%xmm5\n\t"\
		"addpd	%%xmm0,%%xmm4\n\t"\
		"addpd	%%xmm1,%%xmm5\n\t"\
		"movaps	%%xmm4,     (%%eax)\n\t"\
		"movaps	%%xmm5,0x010(%%eax)\n\t"\
		"\n\t"\
		"subpd	%%xmm7,%%xmm2\n\t"\
		"subpd	%%xmm6,%%xmm3\n\t"\
		"movaps	%%xmm2,     (%%edx)\n\t"\
		"movaps	%%xmm3,0x010(%%ecx)\n\t"\
		"addpd	%%xmm7,%%xmm7\n\t"\
		"addpd	%%xmm6,%%xmm6\n\t"\
		"addpd	%%xmm2,%%xmm7\n\t"\
		"addpd	%%xmm3,%%xmm6\n\t"\
		"movaps	%%xmm7,     (%%ecx)\n\t"\
		"movaps	%%xmm6,0x010(%%edx)\n\t"\
		"\n\t"\
		"/*...Block 3: */\n\t"\
		"subl	$0x020,%%eax\n\t"\
		"subl	$0x020,%%esi\n\t"\
		"subl	$0x020,%%ecx\n\t"\
		"subl	$0x020,%%edx\n\t"\
		"/* SSE2_RADIX4_DIT_IN_PLACE(): */\n\t"\
		"movaps	    (%%eax),%%xmm0\n\t"\
		"movaps	0x10(%%eax),%%xmm1\n\t"\
		"movaps	    (%%eax),%%xmm2\n\t"\
		"movaps	0x10(%%eax),%%xmm3\n\t"\
		"\n\t"\
		"addpd	    (%%esi),%%xmm0\n\t"\
		"addpd	0x10(%%esi),%%xmm1\n\t"\
		"subpd	    (%%esi),%%xmm2\n\t"\
		"subpd	0x10(%%esi),%%xmm3\n\t"\
		"\n\t"\
		"movaps	    (%%ecx),%%xmm4\n\t"\
		"movaps	0x10(%%ecx),%%xmm5\n\t"\
		"movaps	    (%%ecx),%%xmm6\n\t"\
		"movaps	0x10(%%ecx),%%xmm7\n\t"\
		"\n\t"\
		"addpd	    (%%edx),%%xmm4\n\t"\
		"addpd	0x10(%%edx),%%xmm5\n\t"\
		"subpd	    (%%edx),%%xmm6\n\t"\
		"subpd	0x10(%%edx),%%xmm7\n\t"\
		"\n\t"\
		"subpd	%%xmm4,%%xmm0\n\t"\
		"subpd	%%xmm5,%%xmm1\n\t"\
		"movaps	%%xmm0,     (%%esi)\n\t"\
		"movaps	%%xmm1,0x010(%%esi)\n\t"\
		"addpd	%%xmm4,%%xmm4\n\t"\
		"addpd	%%xmm5,%%xmm5\n\t"\
		"addpd	%%xmm0,%%xmm4\n\t"\
		"addpd	%%xmm1,%%xmm5\n\t"\
		"movaps	%%xmm4,     (%%eax)\n\t"\
		"movaps	%%xmm5,0x010(%%eax)\n\t"\
		"\n\t"\
		"subpd	%%xmm7,%%xmm2\n\t"\
		"subpd	%%xmm6,%%xmm3\n\t"\
		"movaps	%%xmm2,     (%%edx)\n\t"\
		"movaps	%%xmm3,0x010(%%ecx)\n\t"\
		"addpd	%%xmm7,%%xmm7\n\t"\
		"addpd	%%xmm6,%%xmm6\n\t"\
		"addpd	%%xmm2,%%xmm7\n\t"\
		"addpd	%%xmm3,%%xmm6\n\t"\
		"movaps	%%xmm7,     (%%ecx)\n\t"\
		"movaps	%%xmm6,0x010(%%edx)\n\t"\
		"\n\t"\
		"/*...Block 4: */\n\t"\
		"addl	$0x040,%%eax\n\t"\
		"addl	$0x040,%%esi\n\t"\
		"addl	$0x040,%%ecx\n\t"\
		"addl	$0x040,%%edx\n\t"\
		"/* SSE2_RADIX4_DIT_IN_PLACE(): */\n\t"\
		"movaps	    (%%eax),%%xmm0\n\t"\
		"movaps	0x10(%%eax),%%xmm1\n\t"\
		"movaps	    (%%eax),%%xmm2\n\t"\
		"movaps	0x10(%%eax),%%xmm3\n\t"\
		"\n\t"\
		"addpd	    (%%esi),%%xmm0\n\t"\
		"addpd	0x10(%%esi),%%xmm1\n\t"\
		"subpd	    (%%esi),%%xmm2\n\t"\
		"subpd	0x10(%%esi),%%xmm3\n\t"\
	"prefetcht0	-0x40(%%edi)\n\t"\
		"\n\t"\
		"movaps	    (%%ecx),%%xmm4\n\t"\
		"movaps	0x10(%%ecx),%%xmm5\n\t"\
		"movaps	    (%%ecx),%%xmm6\n\t"\
		"movaps	0x10(%%ecx),%%xmm7\n\t"\
		"\n\t"\
		"addpd	    (%%edx),%%xmm4\n\t"\
		"addpd	0x10(%%edx),%%xmm5\n\t"\
		"subpd	    (%%edx),%%xmm6\n\t"\
		"subpd	0x10(%%edx),%%xmm7\n\t"\
		"\n\t"\
		"subpd	%%xmm4,%%xmm0\n\t"\
		"subpd	%%xmm5,%%xmm1\n\t"\
		"movaps	%%xmm0,     (%%esi)\n\t"\
		"movaps	%%xmm1,0x010(%%esi)\n\t"\
		"addpd	%%xmm4,%%xmm4\n\t"\
		"addpd	%%xmm5,%%xmm5\n\t"\
		"addpd	%%xmm0,%%xmm4\n\t"\
		"addpd	%%xmm1,%%xmm5\n\t"\
		"movaps	%%xmm4,     (%%eax)\n\t"\
		"movaps	%%xmm5,0x010(%%eax)\n\t"\
		"\n\t"\
		"subpd	%%xmm7,%%xmm2\n\t"\
		"subpd	%%xmm6,%%xmm3\n\t"\
		"movaps	%%xmm2,     (%%edx)\n\t"\
		"movaps	%%xmm3,0x010(%%ecx)\n\t"\
		"addpd	%%xmm7,%%xmm7\n\t"\
		"addpd	%%xmm6,%%xmm6\n\t"\
		"addpd	%%xmm2,%%xmm7\n\t"\
		"addpd	%%xmm3,%%xmm6\n\t"\
		"movaps	%%xmm7,     (%%ecx)\n\t"\
		"movaps	%%xmm6,0x010(%%edx)\n\t"\
		"\n\t"\
	/***************************************************************************************************/\
	/*..and now do four more radix-4 transforms, including the internal and external twiddle factors.  */\
	/*  Write even-index 16-byte output pairs to a(j1), odd-index to a(j2), unpack same as on inputs.  */\
	/*  We do the last 2 radix-4 blocks first, to make the unpack-interleaving of outputs smoother.    */\
	/***************************************************************************************************/\
		"\n\t"\
		"/* Main-array addresses still in add0,1, no need to re-init: */\n\t"\
		"\n\t"\
		"/*...Block 3: t3,11,19,27 -> r9,13,11,15: */\n\t"\
		"movl		%[__r9],%%eax\n\t"\
		"movl		%[__isrt2],%%esi\n\t"\
		"movl		%[__cc0],%%ecx\n\t"\
		"\n\t"\
		"movaps		0x020(%%eax),%%xmm4\n\t"\
		"movaps		0x060(%%eax),%%xmm0\n\t"\
		"movaps		0x030(%%eax),%%xmm5\n\t"\
		"movaps		0x070(%%eax),%%xmm1\n\t"\
		"movaps		0x020(%%eax),%%xmm6\n\t"\
		"movaps		0x060(%%eax),%%xmm2\n\t"\
		"movaps		0x030(%%eax),%%xmm7\n\t"\
		"movaps		0x070(%%eax),%%xmm3\n\t"\
		"\n\t"\
		"mulpd		    (%%ecx),%%xmm4\n\t"\
		"mulpd		0x10(%%ecx),%%xmm0\n\t"\
		"mulpd		    (%%ecx),%%xmm5\n\t"\
		"mulpd		0x10(%%ecx),%%xmm1\n\t"\
		"mulpd		0x10(%%ecx),%%xmm6\n\t"\
		"mulpd		    (%%ecx),%%xmm2\n\t"\
		"mulpd		0x10(%%ecx),%%xmm7\n\t"\
		"mulpd		    (%%ecx),%%xmm3\n\t"\
		"subpd		%%xmm6,%%xmm5\n\t"\
		"subpd		%%xmm2,%%xmm1\n\t"\
		"addpd		%%xmm7,%%xmm4\n\t"\
		"addpd		%%xmm3,%%xmm0\n\t"\
		"movaps		%%xmm5,%%xmm7\n\t"\
		"movaps		%%xmm4,%%xmm6\n\t"\
		"\n\t"\
		"addpd		%%xmm0,%%xmm4\n\t"\
		"addpd		%%xmm1,%%xmm5\n\t"\
		"subpd		%%xmm0,%%xmm6\n\t"\
		"subpd		%%xmm1,%%xmm7\n\t"\
		"\n\t"\
		"movaps		0x040(%%eax),%%xmm2\n\t"\
		"movaps		0x050(%%eax),%%xmm3\n\t"\
		"movaps		     (%%eax),%%xmm0\n\t"\
		"movaps		0x010(%%eax),%%xmm1\n\t"\
		"addpd		0x050(%%eax),%%xmm2\n\t"\
		"subpd		0x040(%%eax),%%xmm3\n\t"\
		"mulpd		(%%esi),%%xmm2\n\t"\
		"mulpd		(%%esi),%%xmm3\n\t"\
	"prefetcht0	-0x60(%%edi)\n\t"\
		"\n\t"\
		"subpd		%%xmm2,%%xmm0\n\t"\
		"subpd		%%xmm3,%%xmm1\n\t"\
		"addpd		%%xmm2,%%xmm2\n\t"\
		"addpd		%%xmm3,%%xmm3\n\t"\
		"addpd		%%xmm0,%%xmm2\n\t"\
		"addpd		%%xmm1,%%xmm3\n\t"\
		"\n\t"\
		"movl		%[__add1],%%esi\n\t"\
		"movl		%[__c1],%%ecx\n\t"\
		"movl		%[__c9],%%edx\n\t"\
		"\n\t"\
		"subpd		%%xmm4,%%xmm2\n\t"\
		"subpd		%%xmm5,%%xmm3\n\t"\
		"addpd		%%xmm4,%%xmm4\n\t"\
		"addpd		%%xmm5,%%xmm5\n\t"\
		"addpd		%%xmm2,%%xmm4\n\t"\
		"addpd		%%xmm3,%%xmm5\n\t"\
		"movaps		%%xmm2,     (%%eax)\n\t"\
		"movaps		%%xmm3,0x010(%%eax)\n\t"\
		"movaps		%%xmm4,%%xmm2\n\t"\
		"movaps		%%xmm5,%%xmm3\n\t"\
		"mulpd		    (%%ecx),%%xmm4\n\t"\
		"mulpd		    (%%ecx),%%xmm5\n\t"\
		"mulpd		0x10(%%ecx),%%xmm2\n\t"\
		"mulpd		0x10(%%ecx),%%xmm3\n\t"\
		"subpd		%%xmm2,%%xmm5\n\t"\
		"addpd		%%xmm3,%%xmm4\n\t"\
		"\n\t"\
		"movaps		%%xmm5,0x10(%%esi)\n\t"\
		"movaps		%%xmm4,    (%%esi)\n\t"\
		"\n\t"\
		"movaps		     (%%eax),%%xmm4\n\t"\
		"movaps		0x010(%%eax),%%xmm5\n\t"\
		"movaps		%%xmm4,%%xmm2\n\t"\
		"movaps		%%xmm5,%%xmm3\n\t"\
		"mulpd		    (%%edx),%%xmm4\n\t"\
		"mulpd		    (%%edx),%%xmm5\n\t"\
		"mulpd		0x10(%%edx),%%xmm2\n\t"\
		"mulpd		0x10(%%edx),%%xmm3\n\t"\
		"subpd		%%xmm2,%%xmm5\n\t"\
		"addpd		%%xmm3,%%xmm4\n\t"\
		"\n\t"\
		"movaps		%%xmm5,0x90(%%esi)\n\t"\
		"movaps		%%xmm4,0x80(%%esi)\n\t"\
		"\n\t"\
		"movl		%[__c5],%%ecx\n\t"\
		"movl		%[__c13],%%edx\n\t"\
		"\n\t"\
		"subpd		%%xmm7,%%xmm0\n\t"\
		"subpd		%%xmm6,%%xmm1\n\t"\
		"addpd		%%xmm7,%%xmm7\n\t"\
		"addpd		%%xmm6,%%xmm6\n\t"\
		"addpd		%%xmm0,%%xmm7\n\t"\
		"addpd		%%xmm1,%%xmm6\n\t"\
		"movaps		%%xmm7,%%xmm4\n\t"\
		"movaps		%%xmm1,%%xmm5\n\t"\
		"mulpd		    (%%ecx),%%xmm7\n\t"\
		"mulpd		    (%%ecx),%%xmm1\n\t"\
		"mulpd		0x10(%%ecx),%%xmm4\n\t"\
		"mulpd		0x10(%%ecx),%%xmm5\n\t"\
		"subpd		%%xmm4,%%xmm1\n\t"\
		"addpd		%%xmm5,%%xmm7\n\t"\
		"\n\t"\
		"movaps		%%xmm1,0x50(%%esi)\n\t"\
		"movaps		%%xmm7,0x40(%%esi)\n\t"\
		"\n\t"\
		"movaps		%%xmm0,%%xmm4\n\t"\
		"movaps		%%xmm6,%%xmm5\n\t"\
		"mulpd		    (%%edx),%%xmm0\n\t"\
		"mulpd		    (%%edx),%%xmm6\n\t"\
		"mulpd		0x10(%%edx),%%xmm4\n\t"\
		"mulpd		0x10(%%edx),%%xmm5\n\t"\
		"subpd		%%xmm4,%%xmm6\n\t"\
		"addpd		%%xmm5,%%xmm0\n\t"\
		"\n\t"\
		"movaps		%%xmm6,0xd0(%%esi)\n\t"\
		"movaps		%%xmm0,0xc0(%%esi)\n\t"\
	"prefetcht0	-0x80(%%edi)\n\t"\
		"\n\t"\
		"/*...Block 4: t7,15,23,31 -> r25,29,27,31: */\n\t"\
		"movl		%[__r25],%%eax\n\t"\
		"movl		%[__isrt2],%%esi\n\t"\
		"movl		%[__cc0],%%ecx\n\t"\
		"\n\t"\
		"movaps		0x020(%%eax),%%xmm4\n\t"\
		"movaps		0x060(%%eax),%%xmm0\n\t"\
		"movaps		0x030(%%eax),%%xmm5\n\t"\
		"movaps		0x070(%%eax),%%xmm1\n\t"\
		"movaps		0x020(%%eax),%%xmm6\n\t"\
		"movaps		0x060(%%eax),%%xmm2\n\t"\
		"movaps		0x030(%%eax),%%xmm7\n\t"\
		"movaps		0x070(%%eax),%%xmm3\n\t"\
		"\n\t"\
		"mulpd		0x10(%%ecx),%%xmm4\n\t"\
		"mulpd		    (%%ecx),%%xmm0\n\t"\
		"mulpd		0x10(%%ecx),%%xmm5\n\t"\
		"mulpd		    (%%ecx),%%xmm1\n\t"\
		"mulpd		    (%%ecx),%%xmm6\n\t"\
		"mulpd		0x10(%%ecx),%%xmm2\n\t"\
		"mulpd		    (%%ecx),%%xmm7\n\t"\
		"mulpd		0x10(%%ecx),%%xmm3\n\t"\
		"subpd		%%xmm6,%%xmm5\n\t"\
		"subpd		%%xmm2,%%xmm1\n\t"\
		"addpd		%%xmm7,%%xmm4\n\t"\
		"addpd		%%xmm3,%%xmm0\n\t"\
		"movaps		%%xmm5,%%xmm7\n\t"\
		"movaps		%%xmm4,%%xmm6\n\t"\
		"\n\t"\
		"addpd		%%xmm0,%%xmm4\n\t"\
		"addpd		%%xmm1,%%xmm5\n\t"\
		"subpd		%%xmm0,%%xmm6\n\t"\
		"subpd		%%xmm1,%%xmm7\n\t"\
		"\n\t"\
		"movaps		0x040(%%eax),%%xmm2\n\t"\
		"movaps		0x050(%%eax),%%xmm3\n\t"\
		"movaps		     (%%eax),%%xmm0\n\t"\
		"movaps		0x010(%%eax),%%xmm1\n\t"\
		"subpd		0x050(%%eax),%%xmm2\n\t"\
		"addpd		0x040(%%eax),%%xmm3\n\t"\
		"mulpd		(%%esi),%%xmm2\n\t"\
		"mulpd		(%%esi),%%xmm3\n\t"\
		"\n\t"\
		"subpd		%%xmm2,%%xmm0\n\t"\
		"subpd		%%xmm3,%%xmm1\n\t"\
		"addpd		%%xmm2,%%xmm2\n\t"\
		"addpd		%%xmm3,%%xmm3\n\t"\
		"addpd		%%xmm0,%%xmm2\n\t"\
		"addpd		%%xmm1,%%xmm3\n\t"\
		"\n\t"\
		"movl		%[__add1],%%esi\n\t"\
		"movl		%[__c3],%%ecx\n\t"\
		"movl		%[__c11],%%edx\n\t"\
		"\n\t"\
		"subpd		%%xmm6,%%xmm0\n\t"\
		"subpd		%%xmm7,%%xmm1\n\t"\
		"addpd		%%xmm6,%%xmm6\n\t"\
		"addpd		%%xmm7,%%xmm7\n\t"\
		"addpd		%%xmm0,%%xmm6\n\t"\
		"addpd		%%xmm1,%%xmm7\n\t"\
		"movaps		%%xmm0,     (%%eax)\n\t"\
		"movaps		%%xmm1,0x010(%%eax)\n\t"\
		"movaps		%%xmm6,%%xmm0\n\t"\
		"movaps		%%xmm7,%%xmm1\n\t"\
		"mulpd		    (%%ecx),%%xmm6\n\t"\
		"mulpd		    (%%ecx),%%xmm7\n\t"\
		"mulpd		0x10(%%ecx),%%xmm0\n\t"\
		"mulpd		0x10(%%ecx),%%xmm1\n\t"\
		"subpd		%%xmm0,%%xmm7\n\t"\
		"addpd		%%xmm1,%%xmm6\n\t"\
		"\n\t"\
		"movaps		%%xmm7,0x30(%%esi)\n\t"\
		"movaps		%%xmm6,0x20(%%esi)\n\t"\
		"\n\t"\
		"movaps		     (%%eax),%%xmm6\n\t"\
		"movaps		0x010(%%eax),%%xmm7\n\t"\
		"movaps		%%xmm6,%%xmm0\n\t"\
		"movaps		%%xmm7,%%xmm1\n\t"\
		"mulpd		    (%%edx),%%xmm6\n\t"\
		"mulpd		    (%%edx),%%xmm7\n\t"\
		"mulpd		0x10(%%edx),%%xmm0\n\t"\
		"mulpd		0x10(%%edx),%%xmm1\n\t"\
		"subpd		%%xmm0,%%xmm7\n\t"\
		"addpd		%%xmm1,%%xmm6\n\t"\
	"prefetcht0	-0xa0(%%edi)\n\t"\
		"\n\t"\
		"movaps		%%xmm7,0xb0(%%esi)\n\t"\
		"movaps		%%xmm6,0xa0(%%esi)\n\t"\
		"\n\t"\
		"movl		%[__c7],%%ecx\n\t"\
		"movl		%[__c15],%%edx\n\t"\
		"\n\t"\
		"subpd		%%xmm5,%%xmm2\n\t"\
		"subpd		%%xmm4,%%xmm3\n\t"\
		"addpd		%%xmm5,%%xmm5\n\t"\
		"addpd		%%xmm4,%%xmm4\n\t"\
		"addpd		%%xmm2,%%xmm5\n\t"\
		"addpd		%%xmm3,%%xmm4\n\t"\
		"movaps		%%xmm5,%%xmm0\n\t"\
		"movaps		%%xmm3,%%xmm1\n\t"\
		"mulpd		    (%%ecx),%%xmm5\n\t"\
		"mulpd		    (%%ecx),%%xmm3\n\t"\
		"mulpd		0x10(%%ecx),%%xmm0\n\t"\
		"mulpd		0x10(%%ecx),%%xmm1\n\t"\
		"subpd		%%xmm0,%%xmm3\n\t"\
		"addpd		%%xmm1,%%xmm5\n\t"\
		"\n\t"\
		"movaps		%%xmm3,0x70(%%esi)\n\t"\
		"movaps		%%xmm5,0x60(%%esi)\n\t"\
		"\n\t"\
		"movaps		%%xmm2,%%xmm0\n\t"\
		"movaps		%%xmm4,%%xmm1\n\t"\
		"mulpd		    (%%edx),%%xmm2\n\t"\
		"mulpd		    (%%edx),%%xmm4\n\t"\
		"mulpd		0x10(%%edx),%%xmm0\n\t"\
		"mulpd		0x10(%%edx),%%xmm1\n\t"\
		"subpd		%%xmm0,%%xmm4\n\t"\
		"addpd		%%xmm1,%%xmm2\n\t"\
		"\n\t"\
		"movaps		%%xmm4,0xf0(%%esi)\n\t"\
		"movaps		%%xmm2,0xe0(%%esi)\n\t"\
		"\n\t"\
		"/*...Block 1: t1,9,17,25 -> r1,5,3,7: */\n\t"\
		"movl		%[__r1],%%eax\n\t"\
		"\n\t"\
		"movaps		     (%%eax),%%xmm0\n\t"\
		"movaps		0x010(%%eax),%%xmm1\n\t"\
		"movaps		0x040(%%eax),%%xmm2\n\t"\
		"movaps		0x050(%%eax),%%xmm3\n\t"\
		"\n\t"\
		"subpd		0x040(%%eax),%%xmm0\n\t"\
		"subpd		0x050(%%eax),%%xmm1\n\t"\
		"addpd		     (%%eax),%%xmm2\n\t"\
		"addpd		0x010(%%eax),%%xmm3\n\t"\
		"\n\t"\
		"movaps		0x020(%%eax),%%xmm4\n\t"\
		"movaps		0x030(%%eax),%%xmm5\n\t"\
		"movaps		0x060(%%eax),%%xmm6\n\t"\
		"movaps		0x070(%%eax),%%xmm7\n\t"\
		"\n\t"\
		"subpd		0x060(%%eax),%%xmm4\n\t"\
		"subpd		0x070(%%eax),%%xmm5\n\t"\
		"addpd		0x020(%%eax),%%xmm6\n\t"\
		"addpd		0x030(%%eax),%%xmm7\n\t"\
		"\n\t"\
		"movl		%[__add0],%%eax\n\t"\
		"movl		%[__c8],%%edx\n\t"\
		"\n\t"\
		"addpd		%%xmm6,%%xmm2\n\t"\
		"addpd		%%xmm7,%%xmm3\n\t"\
		"\n\t"\
		"movaps		%%xmm2,    (%%eax)\n\t"\
		"movaps		%%xmm3,0x10(%%eax)\n\t"\
		"addpd		%%xmm6,%%xmm6\n\t"\
		"addpd		%%xmm7,%%xmm7\n\t"\
		"subpd		%%xmm6,%%xmm2\n\t"\
		"subpd		%%xmm7,%%xmm3\n\t"\
		"movaps		%%xmm2,%%xmm6\n\t"\
		"movaps		%%xmm3,%%xmm7\n\t"\
		"mulpd		    (%%edx),%%xmm2\n\t"\
		"mulpd		    (%%edx),%%xmm3\n\t"\
		"mulpd		0x10(%%edx),%%xmm6\n\t"\
		"mulpd		0x10(%%edx),%%xmm7\n\t"\
		"subpd		%%xmm6,%%xmm3\n\t"\
		"addpd		%%xmm7,%%xmm2\n\t"\
	"prefetcht0	-0xc0(%%edi)\n\t"\
		"\n\t"\
		"movl		%[__add1],%%ecx\n\t"\
		"movaps		%%xmm3,%%xmm7\n\t"\
		"movaps		%%xmm2,%%xmm6\n\t"\
		"unpckhpd	0x90(%%ecx),%%xmm7\n\t"\
		"unpcklpd	0x90(%%ecx),%%xmm3\n\t"\
		"movaps		%%xmm7,0x90(%%ecx)\n\t"\
		"unpckhpd	0x80(%%ecx),%%xmm6\n\t"\
		"unpcklpd	0x80(%%ecx),%%xmm2\n\t"\
		"movaps		%%xmm6,0x80(%%ecx)\n\t"\
		"\n\t"\
		"movaps		%%xmm3,0x90(%%eax)\n\t"\
		"movaps		%%xmm2,0x80(%%eax)\n\t"\
		"\n\t"\
		"movaps		0x10(%%eax),%%xmm3\n\t"\
		"movaps		(%%eax),%%xmm2\n\t"\
		"movaps		%%xmm3,%%xmm7\n\t"\
		"movaps		%%xmm2,%%xmm6\n\t"\
		"unpckhpd	0x10(%%ecx),%%xmm7\n\t"\
		"unpcklpd	0x10(%%ecx),%%xmm3\n\t"\
		"movaps		%%xmm7,0x10(%%ecx)\n\t"\
		"unpckhpd	(%%ecx),%%xmm6\n\t"\
		"unpcklpd	(%%ecx),%%xmm2\n\t"\
		"movaps		%%xmm6,    (%%ecx)\n\t"\
		"\n\t"\
		"movaps		%%xmm3,0x10(%%eax)\n\t"\
		"movaps		%%xmm2,    (%%eax)\n\t"\
		"\n\t"\
		"movl		%[__c4],%%ecx\n\t"\
		"movl		%[__c12],%%edx\n\t"\
		"\n\t"\
		"addpd		%%xmm5,%%xmm0\n\t"\
		"subpd		%%xmm4,%%xmm1\n\t"\
		"movaps		%%xmm0,%%xmm2\n\t"\
		"movaps		%%xmm1,%%xmm3\n\t"\
		"addpd		%%xmm5,%%xmm5\n\t"\
		"addpd		%%xmm4,%%xmm4\n\t"\
		"movaps		%%xmm0,%%xmm6\n\t"\
		"movaps		%%xmm1,%%xmm7\n\t"\
		"mulpd		    (%%ecx),%%xmm2\n\t"\
		"mulpd		    (%%ecx),%%xmm3\n\t"\
		"mulpd		0x10(%%ecx),%%xmm6\n\t"\
		"mulpd		0x10(%%ecx),%%xmm7\n\t"\
		"subpd		%%xmm6,%%xmm3\n\t"\
		"addpd		%%xmm7,%%xmm2\n\t"\
		"\n\t"\
		"movl		%[__add1],%%ecx\n\t"\
		"movaps		%%xmm3,%%xmm7\n\t"\
		"movaps		%%xmm2,%%xmm6\n\t"\
		"unpckhpd	0x50(%%ecx),%%xmm7\n\t"\
		"unpcklpd	0x50(%%ecx),%%xmm3\n\t"\
		"movaps		%%xmm7,0x50(%%ecx)\n\t"\
		"unpckhpd	0x40(%%ecx),%%xmm6\n\t"\
		"unpcklpd	0x40(%%ecx),%%xmm2\n\t"\
		"movaps		%%xmm6,0x40(%%ecx)\n\t"\
		"\n\t"\
		"movaps		%%xmm3,0x50(%%eax)\n\t"\
		"movaps		%%xmm2,0x40(%%eax)\n\t"\
		"\n\t"\
		"subpd		%%xmm5,%%xmm0\n\t"\
		"addpd		%%xmm4,%%xmm1\n\t"\
		"movaps		%%xmm0,%%xmm6\n\t"\
		"movaps		%%xmm1,%%xmm7\n\t"\
		"mulpd		    (%%edx),%%xmm0\n\t"\
		"mulpd		    (%%edx),%%xmm1\n\t"\
		"mulpd		0x10(%%edx),%%xmm6\n\t"\
		"mulpd		0x10(%%edx),%%xmm7\n\t"\
		"subpd		%%xmm6,%%xmm1\n\t"\
		"addpd		%%xmm7,%%xmm0\n\t"\
		"\n\t"\
		"movl		%[__add1],%%ecx\n\t"\
		"movaps		%%xmm1,%%xmm7\n\t"\
		"movaps		%%xmm0,%%xmm6\n\t"\
		"unpckhpd	0xd0(%%ecx),%%xmm7\n\t"\
		"unpcklpd	0xd0(%%ecx),%%xmm1\n\t"\
		"movaps		%%xmm7,0xd0(%%ecx)\n\t"\
		"unpckhpd	0xc0(%%ecx),%%xmm6\n\t"\
		"unpcklpd	0xc0(%%ecx),%%xmm0\n\t"\
		"movaps		%%xmm6,0xc0(%%ecx)\n\t"\
		"\n\t"\
	"prefetcht0	-0xe0(%%edi)\n\t"\
		"movaps		%%xmm1,0xd0(%%eax)\n\t"\
		"movaps		%%xmm0,0xc0(%%eax)\n\t"\
		"\n\t"\
		"/*...Block 3: t5,13,21,29 -> r17,21,19,23: */\n\t"\
		"movl		%[__r17],%%eax\n\t"\
		"movl		%[__isrt2],%%esi\n\t"\
		"movaps		(%%esi),%%xmm2\n\t"\
		"\n\t"\
		"movaps		0x020(%%eax),%%xmm4\n\t"\
		"movaps		0x030(%%eax),%%xmm5\n\t"\
		"movaps		0x060(%%eax),%%xmm0\n\t"\
		"movaps		0x070(%%eax),%%xmm1\n\t"\
		"\n\t"\
		"addpd		0x030(%%eax),%%xmm4\n\t"\
		"subpd		0x020(%%eax),%%xmm5\n\t"\
		"subpd		0x070(%%eax),%%xmm0\n\t"\
		"addpd		0x060(%%eax),%%xmm1\n\t"\
		"mulpd		%%xmm2,%%xmm4\n\t"\
		"mulpd		%%xmm2,%%xmm5\n\t"\
		"mulpd		%%xmm2,%%xmm0\n\t"\
		"mulpd		%%xmm2,%%xmm1\n\t"\
		"movaps		%%xmm4,%%xmm6\n\t"\
		"movaps		%%xmm5,%%xmm7\n\t"\
		"\n\t"\
		"subpd		%%xmm0,%%xmm4\n\t"\
		"subpd		%%xmm1,%%xmm5\n\t"\
		"addpd		%%xmm0,%%xmm6\n\t"\
		"addpd		%%xmm1,%%xmm7\n\t"\
		"\n\t"\
		"movaps		     (%%eax),%%xmm0\n\t"\
		"movaps		0x010(%%eax),%%xmm1\n\t"\
		"movaps		0x040(%%eax),%%xmm2\n\t"\
		"movaps		0x050(%%eax),%%xmm3\n\t"\
		"\n\t"\
		"subpd		0x050(%%eax),%%xmm0\n\t"\
		"subpd		0x040(%%eax),%%xmm1\n\t"\
		"addpd		     (%%eax),%%xmm3\n\t"\
		"addpd		0x010(%%eax),%%xmm2\n\t"\
		"\n\t"\
		"movl		%[__add0],%%esi\n\t"\
		"movl		%[__c2],%%ecx\n\t"\
		"movl		%[__c10],%%edx\n\t"\
		"\n\t"\
		"subpd		%%xmm4,%%xmm3\n\t"\
		"subpd		%%xmm5,%%xmm1\n\t"\
		"addpd		%%xmm4,%%xmm4\n\t"\
		"addpd		%%xmm5,%%xmm5\n\t"\
		"addpd		%%xmm3,%%xmm4\n\t"\
		"addpd		%%xmm1,%%xmm5\n\t"\
		"movaps		%%xmm3,     (%%eax)\n\t"\
		"movaps		%%xmm1,0x010(%%eax)\n\t"\
		"movaps		%%xmm4,%%xmm3\n\t"\
		"movaps		%%xmm5,%%xmm1\n\t"\
		"mulpd		    (%%ecx),%%xmm4\n\t"\
		"mulpd		    (%%ecx),%%xmm5\n\t"\
		"mulpd		0x10(%%ecx),%%xmm3\n\t"\
		"mulpd		0x10(%%ecx),%%xmm1\n\t"\
		"subpd		%%xmm3,%%xmm5\n\t"\
		"addpd		%%xmm1,%%xmm4\n\t"\
		"\n\t"\
		"movl		%[__add1],%%ecx\n\t"\
		"movaps		%%xmm5,%%xmm3\n\t"\
		"movaps		%%xmm4,%%xmm1\n\t"\
		"unpckhpd	0x30(%%ecx),%%xmm3\n\t"\
		"unpcklpd	0x30(%%ecx),%%xmm5\n\t"\
		"movaps		%%xmm3,0x30(%%ecx)\n\t"\
		"unpckhpd	0x20(%%ecx),%%xmm1\n\t"\
		"unpcklpd	0x20(%%ecx),%%xmm4\n\t"\
		"movaps		%%xmm1,0x20(%%ecx)\n\t"\
		"\n\t"\
		"movaps		%%xmm5,0x30(%%esi)\n\t"\
		"movaps		%%xmm4,0x20(%%esi)\n\t"\
		"\n\t"\
		"movaps		     (%%eax),%%xmm4\n\t"\
		"movaps		0x010(%%eax),%%xmm5\n\t"\
		"movaps		%%xmm4,%%xmm3\n\t"\
		"movaps		%%xmm5,%%xmm1\n\t"\
	"prefetcht0	-0x100(%%edi)\n\t"\
		"mulpd		    (%%edx),%%xmm4\n\t"\
		"mulpd		    (%%edx),%%xmm5\n\t"\
		"mulpd		0x10(%%edx),%%xmm3\n\t"\
		"mulpd		0x10(%%edx),%%xmm1\n\t"\
		"subpd		%%xmm3,%%xmm5\n\t"\
		"addpd		%%xmm1,%%xmm4\n\t"\
		"\n\t"\
		"movaps		%%xmm5,%%xmm3\n\t"\
		"movaps		%%xmm4,%%xmm1\n\t"\
		"unpckhpd	0xb0(%%ecx),%%xmm3\n\t"\
		"unpcklpd	0xb0(%%ecx),%%xmm5\n\t"\
		"movaps		%%xmm3,0xb0(%%ecx)\n\t"\
		"unpckhpd	0xa0(%%ecx),%%xmm1\n\t"\
		"unpcklpd	0xa0(%%ecx),%%xmm4\n\t"\
		"movaps		%%xmm1,0xa0(%%ecx)\n\t"\
		"\n\t"\
		"movaps		%%xmm5,0xb0(%%esi)\n\t"\
		"movaps		%%xmm4,0xa0(%%esi)\n\t"\
		"\n\t"\
		"movl		%[__c6],%%ecx\n\t"\
		"movl		%[__c14],%%edx\n\t"\
		"\n\t"\
		"subpd		%%xmm7,%%xmm0\n\t"\
		"subpd		%%xmm6,%%xmm2\n\t"\
		"addpd		%%xmm7,%%xmm7\n\t"\
		"addpd		%%xmm6,%%xmm6\n\t"\
		"addpd		%%xmm0,%%xmm7\n\t"\
		"addpd		%%xmm2,%%xmm6\n\t"\
		"movaps		%%xmm7,%%xmm4\n\t"\
		"movaps		%%xmm2,%%xmm5\n\t"\
		"mulpd		    (%%ecx),%%xmm7\n\t"\
		"mulpd		    (%%ecx),%%xmm2\n\t"\
		"mulpd		0x10(%%ecx),%%xmm4\n\t"\
		"mulpd		0x10(%%ecx),%%xmm5\n\t"\
		"subpd		%%xmm4,%%xmm2\n\t"\
		"addpd		%%xmm5,%%xmm7\n\t"\
		"\n\t"\
		"movl		%[__add1],%%ecx\n\t"\
		"movaps		%%xmm2,%%xmm5\n\t"\
		"movaps		%%xmm7,%%xmm4\n\t"\
		"unpckhpd	0x70(%%ecx),%%xmm5\n\t"\
		"unpcklpd	0x70(%%ecx),%%xmm2\n\t"\
		"movaps		%%xmm5,0x70(%%ecx)\n\t"\
		"unpckhpd	0x60(%%ecx),%%xmm4\n\t"\
		"unpcklpd	0x60(%%ecx),%%xmm7\n\t"\
		"movaps		%%xmm4,0x60(%%ecx)\n\t"\
		"\n\t"\
		"movaps		%%xmm2,0x70(%%esi)\n\t"\
		"movaps		%%xmm7,0x60(%%esi)\n\t"\
		"\n\t"\
		"movaps		%%xmm0,%%xmm4\n\t"\
		"movaps		%%xmm6,%%xmm5\n\t"\
		"mulpd		    (%%edx),%%xmm0\n\t"\
		"mulpd		    (%%edx),%%xmm6\n\t"\
		"mulpd		0x10(%%edx),%%xmm4\n\t"\
		"mulpd		0x10(%%edx),%%xmm5\n\t"\
		"subpd		%%xmm4,%%xmm6\n\t"\
		"addpd		%%xmm5,%%xmm0\n\t"\
		"\n\t"\
		"movaps		%%xmm6,%%xmm5\n\t"\
		"movaps		%%xmm0,%%xmm4\n\t"\
		"unpckhpd	0xf0(%%ecx),%%xmm5\n\t"\
		"unpcklpd	0xf0(%%ecx),%%xmm6\n\t"\
		"movaps		%%xmm5,0xf0(%%ecx)\n\t"\
		"unpckhpd	0xe0(%%ecx),%%xmm4\n\t"\
		"unpcklpd	0xe0(%%ecx),%%xmm0\n\t"\
		"movaps		%%xmm4,0xe0(%%ecx)\n\t"\
		"\n\t"\
		"movaps		%%xmm6,0xf0(%%esi)\n\t"\
		"movaps		%%xmm0,0xe0(%%esi)\n\t"\
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
		: "cc","memory","eax","esi","ecx","edx","edi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"	/* Clobbered registers */\
	);\
	}

#endif	/* radix16_wrapper_square_gcc_h_included */

