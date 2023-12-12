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

/****************************************************************************
 * We now include this header file if it was not included before.
 ****************************************************************************/
#ifndef pair_square_included
#define pair_square_included

/***************/

/*
!   Given complex scalars H[j] = (x1,y1) and H[N-j] = (x2,y2) along with complex exponential E = (c,s),
!   calculates I[j] = H[j]^2 + {1 + exp(2*pi*I*j/N)}*{H[j]-H~[N-j]}^2/4 and its complex conjugate I~,
!   returns the former in H[j] and the latter in H[N-j].
*/
#define PAIR_SQUARE2A(__x1, __y1, __x2, __y2, __x3, __y3, __x4, __y4, __c, __s)\
{\
	double __rt0,__rt3,__it3, __st0,__st3,__jt3;\
\
/*   calculate cross-product terms... */\
	__rt3=__x1*__x2+__y1*__y2; __rt3=__rt3+__rt3;\
	__st3=__x3*__x4+__y3*__y4; __st3=__st3+__st3;\
\
	__it3=__y1*__x2-__x1*__y2; __it3=__it3+__it3;\
	__jt3=__y3*__x4-__x3*__y4; __jt3=__jt3+__jt3;\
\
/*   now calculate square terms and __store back in the same temporaries. */\
	__rt0=(__x1+__y1)*(__x1-__y1); __y1=__x1*__y1; __y1=__y1+__y1; __x1=__rt0;\
	__st0=(__x3+__y3)*(__x3-__y3); __y3=__x3*__y3; __y3=__y3+__y3; __x3=__st0;\
\
	__rt0=(__x2+__y2)*(__x2-__y2); __y2=__x2*__y2; __y2=__y2+__y2; __x2=__rt0;\
	__st0=(__x4+__y4)*(__x4-__y4); __y4=__x4*__y4; __y4=__y4+__y4; __x4=__st0;\
\
/*   use that (H[j] - H~[N-j])^2 = H(j)^2 - 2*H(j)*H~(N-j) + H~(N-j)^2... */\
	__rt3=__x1+__x2-__rt3;\
	__st3=__x3+__x4-__st3;\
\
	__it3=__y1-__y2-__it3;\
	__jt3=__y3-__y4-__jt3;\
\
	__rt0=((1.0+__c)*__rt3-__s*__it3)*0.25;\
	__st0=((1.0-__s)*__st3-__c*__jt3)*0.25;	/* [c,s] --> [-c,-s] here */\
\
	__it3=((1.0+__c)*__it3+__s*__rt3)*0.25;\
	__jt3=((1.0-__s)*__jt3+__c*__st3)*0.25;	/* [c,s] --> [-c,-s] here */\
\
/*...and now complete and store the results. */\
	__x1 = (__x1-__rt0);\
	__x3 = (__x3-__st0);\
\
	__y1 = (__y1-__it3);\
	__y3 = (__y3-__jt3);\
\
/*...N-j terms are as above, but with the replacements: __x1<-->__x2, __y1<-->__y2, __it3|-->-__it3. */\
	__x2 = (__x2-__rt0);\
	__x4 = (__x4-__st0);\
\
	__y2 = (__y2+__it3);\
	__y4 = (__y4+__jt3);\
}

/***************/

/*
!   Same as PAIR_SQUARE2A, except signs flipped on [c,s] in computation of st0 and jt3.
*/
#define PAIR_SQUARE2B(__x1, __y1, __x2, __y2, __x3, __y3, __x4, __y4, __c, __s)\
{\
	double __rt0,__rt3,__it3, __st0,__st3,__jt3;\
\
/*   calculate cross-product terms... */\
	__rt3=__x1*__x2+__y1*__y2; __rt3=__rt3+__rt3;\
	__st3=__x3*__x4+__y3*__y4; __st3=__st3+__st3;\
\
	__it3=__y1*__x2-__x1*__y2; __it3=__it3+__it3;\
	__jt3=__y3*__x4-__x3*__y4; __jt3=__jt3+__jt3;\
\
/*   now calculate square terms and __store back in the same temporaries. */\
	__rt0=(__x1+__y1)*(__x1-__y1); __y1=__x1*__y1; __y1=__y1+__y1; __x1=__rt0;\
	__st0=(__x3+__y3)*(__x3-__y3); __y3=__x3*__y3; __y3=__y3+__y3; __x3=__st0;\
\
	__rt0=(__x2+__y2)*(__x2-__y2); __y2=__x2*__y2; __y2=__y2+__y2; __x2=__rt0;\
	__st0=(__x4+__y4)*(__x4-__y4); __y4=__x4*__y4; __y4=__y4+__y4; __x4=__st0;\
\
/*   use that (H[j] - H~[N-j])^2 = H(j)^2 - 2*H(j)*H~(N-j) + H~(N-j)^2... */\
	__rt3=__x1+__x2-__rt3;\
	__st3=__x3+__x4-__st3;\
\
	__it3=__y1-__y2-__it3;\
	__jt3=__y3-__y4-__jt3;\
\
	__rt0=((1.0-__c)*__rt3-__s*__it3)*0.25;\
	__st0=((1.0+__s)*__st3-__c*__jt3)*0.25;	/* [1+c,1-s] --> [1-c,1+s] here */\
\
	__it3=((1.0-__c)*__it3+__s*__rt3)*0.25;\
	__jt3=((1.0+__s)*__jt3+__c*__st3)*0.25;	/* [1+c,1-s] --> [1-c,1+s] here */\
\
/*...and now complete and store the results. */\
	__x1 = (__x1-__rt0);\
	__x3 = (__x3-__st0);\
\
	__y1 = (__y1-__it3);\
	__y3 = (__y3-__jt3);\
\
/*...N-j terms are as above, but with the replacements: __x1<-->__x2, __y1<-->__y2, __it3|-->-__it3. */\
	__x2 = (__x2-__rt0);\
	__x4 = (__x4-__st0);\
\
	__y2 = (__y2+__it3);\
	__y4 = (__y4+__jt3);\
}

/***************/

#define PAIR_SQUARE_4(__aj1pAr, __aj1pAi, __aj2pDr, __aj2pDi, __aj1pBr, __aj1pBi, __aj2pCr, __aj2pCi, __c0, __s0 \
					, __aj2pBr, __aj2pBi, __aj1pCr, __aj1pCi, __aj2pAr, __aj2pAi, __aj1pDr, __aj1pDi, __c1, __s1)\
{\
	double __rt0,__rt1, __it0,__it1, __st0,__st1, __jt0,__jt1;\
	double __tmp, __rt0_tmp,__rt1_tmp, __st0_tmp,__st1_tmp;\
\
/*   calculate cross-product terms... */\
	__rt0=__aj1pAr*__aj2pDr+__aj1pAi*__aj2pDi; __rt0=__rt0+__rt0;	/* <*** ~D*/\
	__st1=__aj2pAr*__aj1pDr+__aj2pAi*__aj1pDi; __st1=__st1+__st1;\
\
	__it0=__aj1pAi*__aj2pDr-__aj1pAr*__aj2pDi; __it0=__it0+__it0;	/* <*** ~D*/\
	__jt1=__aj2pAi*__aj1pDr-__aj2pAr*__aj1pDi; __jt1=__jt1+__jt1;\
\
	__st0=__aj1pBr*__aj2pCr+__aj1pBi*__aj2pCi; __st0=__st0+__st0;	/* <*** ~C*/\
	__rt1=__aj2pBr*__aj1pCr+__aj2pBi*__aj1pCi; __rt1=__rt1+__rt1;\
\
	__jt0=__aj1pBi*__aj2pCr-__aj1pBr*__aj2pCi; __jt0=__jt0+__jt0;	/* <*** ~C*/\
	__it1=__aj2pBi*__aj1pCr-__aj2pBr*__aj1pCi; __it1=__it1+__it1;\
\
/*   now calculate square terms and __store back in the same temporaries. */\
	__tmp=(__aj1pAr+__aj1pAi)*(__aj1pAr-__aj1pAi); __aj1pAi=__aj1pAr*__aj1pAi; __aj1pAi=__aj1pAi+__aj1pAi; __aj1pAr=__tmp;\
	__tmp=(__aj2pAr+__aj2pAi)*(__aj2pAr-__aj2pAi); __aj2pAi=__aj2pAr*__aj2pAi; __aj2pAi=__aj2pAi+__aj2pAi; __aj2pAr=__tmp;\
\
	__tmp=(__aj1pBr+__aj1pBi)*(__aj1pBr-__aj1pBi); __aj1pBi=__aj1pBr*__aj1pBi; __aj1pBi=__aj1pBi+__aj1pBi; __aj1pBr=__tmp;\
	__tmp=(__aj2pBr+__aj2pBi)*(__aj2pBr-__aj2pBi); __aj2pBi=__aj2pBr*__aj2pBi; __aj2pBi=__aj2pBi+__aj2pBi; __aj2pBr=__tmp;\
\
	__tmp=(__aj2pDr+__aj2pDi)*(__aj2pDr-__aj2pDi); __aj2pDi=__aj2pDr*__aj2pDi; __aj2pDi=__aj2pDi+__aj2pDi; __aj2pDr=__tmp;\
	__tmp=(__aj1pDr+__aj1pDi)*(__aj1pDr-__aj1pDi); __aj1pDi=__aj1pDr*__aj1pDi; __aj1pDi=__aj1pDi+__aj1pDi; __aj1pDr=__tmp;\
\
	__tmp=(__aj2pCr+__aj2pCi)*(__aj2pCr-__aj2pCi); __aj2pCi=__aj2pCr*__aj2pCi; __aj2pCi=__aj2pCi+__aj2pCi; __aj2pCr=__tmp;\
	__tmp=(__aj1pCr+__aj1pCi)*(__aj1pCr-__aj1pCi); __aj1pCi=__aj1pCr*__aj1pCi; __aj1pCi=__aj1pCi+__aj1pCi; __aj1pCr=__tmp;\
\
/*   use that (H[j] - H~[N-j])^2 = H(j)^2 - 2*H(j)*H~(N-j) + H~(N-j)^2... */\
	__rt0=__aj1pAr+__aj2pDr-__rt0;	/* <*** ~D*/\
	__st1=__aj2pAr+__aj1pDr-__st1;\
\
	__it0=__aj1pAi-__aj2pDi-__it0;	/* <*** ~D*/\
	__jt1=__aj2pAi-__aj1pDi-__jt1;\
\
	__st0=__aj1pBr+__aj2pCr-__st0;	/* <*** ~C*/\
	__rt1=__aj2pBr+__aj1pCr-__rt1;\
\
	__jt0=__aj1pBi-__aj2pCi-__jt0;	/* <*** ~C*/\
	__it1=__aj2pBi-__aj1pCi-__it1;\
\
	__rt0_tmp=((1.0+__c0)*__rt0-__s0*__it0)*0.25;\
	__st1_tmp=((1.0+__s1)*__st1-__c1*__jt1)*0.25;\
\
	__it0=((1.0+__c0)*__it0+__s0*__rt0)*0.25;	__rt0=__rt0_tmp;\
	__jt1=((1.0+__s1)*__jt1+__c1*__st1)*0.25;	__st1=__st1_tmp;\
\
	__st0_tmp=((1.0-__s0)*__st0-__c0*__jt0)*0.25;\
	__rt1_tmp=((1.0-__c1)*__rt1-__s1*__it1)*0.25;\
\
	__jt0=((1.0-__s0)*__jt0+__c0*__st0)*0.25;	__st0=__st0_tmp;\
	__it1=((1.0-__c1)*__it1+__s1*__rt1)*0.25;	__rt1=__rt1_tmp;\
\
/*...and now complete and store the results. */\
	__aj1pAr = (__aj1pAr-__rt0);\
	__aj2pAr = (__aj2pAr-__st1);\
\
	__aj1pAi = (__aj1pAi-__it0);\
	__aj2pAi = (__aj2pAi-__jt1);\
\
	__aj1pBr = (__aj1pBr-__st0);\
	__aj2pBr = (__aj2pBr-__rt1);\
\
	__aj1pBi = (__aj1pBi-__jt0);\
	__aj2pBi = (__aj2pBi-__it1);\
\
/*...N-j terms are as above, but with the replacements: __aj1pAr<-->__aj2pDr, __aj1pAi<-->__aj2pDi, __it1|-->-__it1. */\
	__aj2pDr = (__aj2pDr-__rt0);	/* <*** ~RT0*/\
	__aj1pDr = (__aj1pDr-__st1);\
\
	__aj2pDi = (__aj2pDi+__it0);	/* <*** ~IT1*/\
	__aj1pDi = (__aj1pDi+__jt1);\
\
	__aj2pCr = (__aj2pCr-__st0);	/* <*** ~ST0*/\
	__aj1pCr = (__aj1pCr-__rt1);\
\
	__aj2pCi = (__aj2pCi+__jt0);	/* <*** ~JT1*/\
	__aj1pCi = (__aj1pCi+__it1);\
}

/* 4-way version of pair_mul() macro in radix16_wrapper_square.c, equivalent to interleaving the 4 calls
	pair_mul(*ax1, *ay1, *ax2, *ay2, ax3, ay3, ax4, ay4, c0, s0);
	pair_mul(*bx1, *by1, *bx2, *by2, bx3, by3, bx4, by4,-s0, c0);
	pair_mul(*cx1, *cy1, *cx2, *cy2, cx3, cy3, cx4, cy4,-c1, s1);
	pair_mul(*dx1, *dy1, *dx2, *dy2, dx3, dy3, dx4, dy4, s1, c1);
...and we rearrange the input args in our implementation below to group pairs of inputs in the same manner which said
inputs will experience in the SSE2 vector version of the code, whose inline-asm pair-mul macro is derived from this one.
Specifically, in the SSE2 version the inputs would be paired as [a|b][0,1,2,3][d0,d1] and we swap the d0,d1 components
of [a|b][2,3] after loading those data into xmm registers:
*/
#define PAIR_MUL_4(\
	__a0rd0, __a0id0, __a3rd1, __a3id1, __b0rd0, __b0id0, __b3rd1, __b3id1,\
	__a0rd1, __a0id1, __a3rd0, __a3id0, __b0rd1, __b0id1, __b3rd0, __b3id0,\
	__a1rd0, __a1id0, __a2rd1, __a2id1, __b1rd0, __b1id0, __b2rd1, __b2id1,\
	__a1rd1, __a1id1, __a2rd0, __a2id0, __b1rd1, __b1id1, __b2rd0, __b2id0, __c0, __s0, __c1, __s1\
)\
{\
	double __cc0 = __c0*.25,__ss0 = __s0*.25, __cc1 = __c1*.25,__ss1 = __s1*.25;\
	double __t0d0,__t0d1, __t1d0,__t1d1, __t2d0,__t2d1, __t3d0,__t3d1, __t4d0,__t4d1, __t5d0,__t5d1;\
/* calculate 2nd square-like term and store in temp...these don't care if the [a,b][2|3] data are d0,1-swapped or not, but do those
   swaps because the following 'calculate difference terms' *does* need such a [lo,hi] double-data swap in the vector version: */\
	__t0d0 = __a3rd1*__b3rd1 - __a3id1*__b3id1;	__t2d0 = __a3rd1*__b3id1 + __a3id1*__b3rd1;\
	__t0d1 = __a3rd0*__b3rd0 - __a3id0*__b3id0;	__t2d1 = __a3rd0*__b3id0 + __a3id0*__b3rd0;\
	__t1d0 = __a2rd1*__b2rd1 - __a2id1*__b2id1;	__t3d0 = __a2rd1*__b2id1 + __a2id1*__b2rd1;\
	__t1d1 = __a2rd0*__b2rd0 - __a2id0*__b2id0;	__t3d1 = __a2rd0*__b2id0 + __a2id0*__b2rd0;\
/* calculate difference terms...these need the [a,b][2|3] vector-data to be d0,1-swapped: */\
	__a3rd1 = __a0rd0 - __a3rd1;	__a3id1 = __a0id0 + __a3id1;	__b3rd1 = __b0rd0 - __b3rd1;	__b3id1 = __b0id0 + __b3id1;\
	__a3rd0 = __a0rd1 - __a3rd0;	__a3id0 = __a0id1 + __a3id0;	__b3rd0 = __b0rd1 - __b3rd0;	__b3id0 = __b0id1 + __b3id0;\
	__a2rd1 = __a1rd0 - __a2rd1;	__a2id1 = __a1id0 + __a2id1;	__b2rd1 = __b1rd0 - __b2rd1;	__b2id1 = __b1id0 + __b2id1;\
	__a2rd0 = __a1rd1 - __a2rd0;	__a2id0 = __a1id1 + __a2id0;	__b2rd0 = __b1rd1 - __b2rd0;	__b2id0 = __b1id1 + __b2id0;\
/* now calculate 1st square-like term and store back in H(j) slot... */\
	__t4d0 = __a0rd0*__b0rd0 - __a0id0*__b0id0;	__a0id0 = __a0rd0*__b0id0 + __a0id0*__b0rd0; __a0rd0 = __t4d0;\
	__t4d1 = __a0rd1*__b0rd1 - __a0id1*__b0id1;	__a0id1 = __a0rd1*__b0id1 + __a0id1*__b0rd1; __a0rd1 = __t4d1;\
	__t5d0 = __a1rd0*__b1rd0 - __a1id0*__b1id0;	__a1id0 = __a1rd0*__b1id0 + __a1id0*__b1rd0; __a1rd0 = __t5d0;\
	__t5d1 = __a1rd1*__b1rd1 - __a1id1*__b1id1;	__a1id1 = __a1rd1*__b1id1 + __a1id1*__b1rd1; __a1rd1 = __t5d1;\
/* calculate the complex products to build the second term... */\
	__t4d0 = __a3rd1*__b3rd1 - __a3id1*__b3id1;	__a3id1 = __a3rd1*__b3id1 + __a3id1*__b3rd1;\
	__t4d1 = __a3rd0*__b3rd0 - __a3id0*__b3id0;	__a3id0 = __a3rd0*__b3id0 + __a3id0*__b3rd0;\
	__t5d0 = __a2rd1*__b2rd1 - __a2id1*__b2id1;	__a2id1 = __a2rd1*__b2id1 + __a2id1*__b2rd1;\
	__t5d1 = __a2rd0*__b2rd0 - __a2id0*__b2id0;	__a2id0 = __a2rd0*__b2id0 + __a2id0*__b2rd0;\
/* Note: precomputing the +-0.25 combos here means trading 2 adds here for (1 add, 1 mul) of precomputation, i.e. no point.
   Note also the mixing of the [cc0,ss1],[ss0,cc1] sincos data in the vector version of this sequence: */\
	__a3rd1 = ((__cc0+0.25)*__t4d0 - __ss0*__a3id1);	__a3id1 = (__ss0*__t4d0 + (__cc0+0.25)*__a3id1);	/*  __cc0,__ss0 */\
	__a3rd0 = ((__ss1+0.25)*__t4d1 - __cc1*__a3id0);	__a3id0 = (__cc1*__t4d1 + (__ss1+0.25)*__a3id0);	/*  __ss1,__cc1 */\
	__a2rd1 = ((0.25-__ss0)*__t5d0 - __cc0*__a2id1);	__a2id1 = (__cc0*__t5d0 + (0.25-__ss0)*__a2id1);	/* -__ss0,__cc0 */\
	__a2rd0 = ((0.25-__cc1)*__t5d1 - __ss1*__a2id0);	__a2id0 = (__ss1*__t5d1 + (0.25-__cc1)*__a2id0);	/* -__cc1,__ss1 */\
/* and now complete and store the results. */\
	__a0rd0 = (__a0rd0-__a3rd1);	__a0id0 = (__a0id0-__a3id1);\
	__a0rd1 = (__a0rd1-__a3rd0);	__a0id1 = (__a0id1-__a3id0);\
	__a1rd0 = (__a1rd0-__a2rd1);	__a1id0 = (__a1id0-__a2id1);\
	__a1rd1 = (__a1rd1-__a2rd0);	__a1id1 = (__a1id1-__a2id0);\
/* N-j terms: */\
	__a3rd1 = (__t0d0-__a3rd1);	__a3id1 = (__t2d0+__a3id1);\
	__a3rd0 = (__t0d1-__a3rd0);	__a3id0 = (__t2d1+__a3id0);\
	__a2rd1 = (__t1d0-__a2rd1);	__a2id1 = (__t3d0+__a2id1);\
	__a2rd0 = (__t1d1-__a2rd0);	__a2id0 = (__t3d1+__a2id0);\
}

#endif	/* #ifndef pair_square_included */
