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

/****************************************************************************
 * We now include this header file if it was not included before.
 ****************************************************************************/
#ifndef imul_macro256_h_included
#define imul_macro256_h_included

#include "types.h"
#include "util.h"

/*
Wide integer multiply macros, with ASM fragments to access efficient non-standard-C
machine instructions wherever feasible. Also includes macros for wide-integer add,
subtract, < > <= == comparison, and bitwise logical shift.
*/

#ifdef __cplusplus
extern "C" {
#endif

#undef	MAX_IMUL_WORDS
#define MAX_IMUL_WORDS	4	/* Max. # of 64-bit words allowed for multiplicands;
							   64 x this should match the bit count of the largest MUL macro. */

/*****************************************************/
/*                    256 x 256                      */
/*****************************************************/

/* Full 512 bits of the product of uint256 x and y.

Let x = x0 + x1*2^64 + x2*2^128 + x3*2^192,  with x0, x1, x2 < 2^64,
and y = y0 + y1*2^64 + y2*2^128 + y3*2^192,  with y0, y1, y2 < 2^64.

Then,

 x*y = 2^  0 * (x0*y0                        )
     + 2^ 64 * (x0*y1 + x1*y0                )
     + 2^128 * (x0*y2 + x1*y1 + x2*y0        )
     + 2^192 * (x0*y3 + x1*y2 + x2*y1 + x3*y0)
     + 2^256 * (        x1*y3 + x2*y2 + x3*y1)
     + 2^320 * (                x2*y3 + x3*y2)
     + 2^384 * (                        x3*y3) .

In terms of the 8 output coefficients, here is what goes into each of those:

	w0 = (x0*y0).lo
	w1 = (x0*y0).hi + (x0*y1 + x1*y0).lo
	w2 =              (x0*y1 + x1*y0).hi + (x0*y2 + x1*y1 + x2*y0).lo
	w3 =                                   (x0*y2 + x1*y1 + x2*y0).hi + (x0*y3 + x1*y2 + x2*y1 + x3*y0).lo
	w4 =                                                                (x0*y3 + x1*y2 + x2*y1 + x3*y0).hi + (x1*y3 + x2*y2 + x3*y1).lo
	w5 =                                                                                                     (x1*y3 + x2*y2 + x3*y1).hi + (x2*y3 + x3*y2).lo
	w6 =                                                                                                                                  (x2*y3 + x3*y2).hi + (x3*y3).lo
	w7 =                                                                                                                                                       (x3*y3).hi

On Alpha, this needs a total of 32 MUL instructions and 82 ALU ops.
*/
#ifdef MUL_LOHI64_SUBROUTINE

	#define MUL_LOHI256(__x, __y, __lo, __hi)\
	{\
		uint64 __w0, __w1, __w2, __w3, __w4, __w5, __w6, __w7\
						,__cy2,__cy3,__cy4,__cy5,__cy6,__cy7\
		,__a,__b,__c,__d,__e,__f,__g,__h,__i,__j,__k,__l\
		,__m,__n,__o,__p,__q,__r,__s,__t,__u,__v,__w,__z;\
		\
		MUL_LOHI64(__x.d0,__y.d0,&__w0,&__w1);	/*   x0*y0 */\
		MUL_LOHI64(__x.d1,__y.d1,&__w2,&__w3);	/*   x1*y1 */\
		MUL_LOHI64(__x.d2,__y.d2,&__w4,&__w5);	/*   x2*y2 */\
		MUL_LOHI64(__x.d3,__y.d3,&__w6,&__w7);	/*   x3*y3 */\
		\
		MUL_LOHI64(__x.d0,__y.d1,&__a ,&__b );	/*   x0*y1 */\
		MUL_LOHI64(__x.d1,__y.d0,&__c ,&__d );	/*   x1*y0 */\
		\
		MUL_LOHI64(__x.d0,__y.d2,&__e ,&__f );	/*   x0*y2 */\
		MUL_LOHI64(__x.d2,__y.d0,&__g ,&__h );	/*   x2*y0 */\
		\
		MUL_LOHI64(__x.d0,__y.d3,&__i ,&__j );	/*   x0*y3 */\
		MUL_LOHI64(__x.d1,__y.d2,&__k ,&__l );	/*   x1*y2 */\
		MUL_LOHI64(__x.d2,__y.d1,&__m ,&__n );	/*   x2*y1 */\
		MUL_LOHI64(__x.d3,__y.d0,&__o ,&__p );	/*   x3*y0 */\
		\
		MUL_LOHI64(__x.d1,__y.d3,&__q ,&__r );	/*   x1*y3 */\
		MUL_LOHI64(__x.d3,__y.d1,&__s ,&__t );	/*   x3*y1 */\
		\
		MUL_LOHI64(__x.d2,__y.d3,&__u ,&__v );	/*   x2*y3 */\
		MUL_LOHI64(__x.d3,__y.d2,&__w ,&__z );	/*   x3*y2 */\
		\
		/* Now add cross terms: */\
		/* Add x0*y1 to w1-2: */\
		__w1 += __a;	__cy2  = (__w1 < __a);\
		__w2 += __b;	__cy3  = (__w2 < __b);\
		\
		/* Add x1*y0 to w1-2: */\
		__w1 += __c;	__cy2 += (__w1 < __c);\
		__w2 += __d;	__cy3 += (__w2 < __d);\
		\
		/* Add x0*y2 to w2-3: */\
		__w2 += __e;	__cy3 += (__w2 < __e);\
		__w3 += __f;	__cy4  = (__w3 < __f);\
		\
		/* Add x2*y0 to w2-3: */\
		__w2 += __g;	__cy3 += (__w2 < __g);\
		__w3 += __h;	__cy4 += (__w3 < __h);\
		\
		/* Add x0*y3 to w3-4: */\
		__w3 += __i;	__cy4 += (__w3 < __i);\
		__w4 += __j;	__cy5  = (__w4 < __j);\
		\
		/* Add x1*y2 to w3-4: */\
		__w3 += __k;	__cy4 += (__w3 < __k);\
		__w4 += __l;	__cy5 += (__w4 < __l);\
		\
		/* Add x2*y1 to w3-4: */\
		__w3 += __m;	__cy4 += (__w3 < __m);\
		__w4 += __n;	__cy5 += (__w4 < __n);\
		\
		/* Add x3*y0 to w3-4: */\
		__w3 += __o;	__cy4 += (__w3 < __o);\
		__w4 += __p;	__cy5 += (__w4 < __p);\
		\
		/* Add x1*y3 to w4-5: */\
		__w4 += __q;	__cy5 += (__w4 < __q);\
		__w5 += __r;	__cy6  = (__w5 < __r);\
		\
		/* Add x3*y1 to w4-5: */\
		__w4 += __s;	__cy5 += (__w4 < __s);\
		__w5 += __t;	__cy6 += (__w5 < __t);\
		\
		/* Add x2*y3 to w5-6: */\
		__w5 += __u;	__cy6 += (__w5 < __u);\
		__w6 += __v;	__cy7  = (__w6 < __v);\
		\
		/* Add x3*y2 to w5-6: */\
		__w5 += __w;	__cy6 += (__w5 < __w);\
		__w6 += __z;	__cy7 += (__w6 < __z);\
		\
		/* Now process carries: */\
		__w2 += __cy2;	__cy3 += (__w2 < __cy2);\
		__w3 += __cy3;	__cy4 += (__w3 < __cy3);\
		__w4 += __cy4;	__cy5 += (__w4 < __cy4);\
		__w5 += __cy5;	__cy6 += (__w5 < __cy5);\
		__w6 += __cy6;	__cy7 += (__w6 < __cy6);\
		__w7 += __cy7;\
		\
		/* Now split the result between __lo and __hi: */\
		__lo.d0 = __w0; __lo.d1 = __w1; __lo.d2 = __w2; __lo.d3 = __w3;\
		__hi.d0 = __w4;	__hi.d1 = __w5;	__hi.d2 = __w6;	__hi.d3 = __w7;\
	}

#else

	#define MUL_LOHI256(__x, __y, __lo, __hi)\
	{\
		uint64 __w0, __w1, __w2, __w3, __w4, __w5, __w6, __w7\
						,__cy2,__cy3,__cy4,__cy5,__cy6,__cy7\
		,__a,__b,__c,__d,__e,__f,__g,__h,__i,__j,__k,__l\
		,__m,__n,__o,__p,__q,__r,__s,__t,__u,__v,__w,__z;\
		\
		MUL_LOHI64(__x.d0,__y.d0, __w0, __w1);	/*   x0*y0 */\
		MUL_LOHI64(__x.d1,__y.d1, __w2, __w3);	/*   x1*y1 */\
		MUL_LOHI64(__x.d2,__y.d2, __w4, __w5);	/*   x2*y2 */\
		MUL_LOHI64(__x.d3,__y.d3, __w6, __w7);	/*   x3*y3 */\
		\
		MUL_LOHI64(__x.d0,__y.d1, __a , __b );	/*   x0*y1 */\
		MUL_LOHI64(__x.d1,__y.d0, __c , __d );	/*   x1*y0 */\
		\
		MUL_LOHI64(__x.d0,__y.d2, __e , __f );	/*   x0*y2 */\
		MUL_LOHI64(__x.d2,__y.d0, __g , __h );	/*   x2*y0 */\
		\
		MUL_LOHI64(__x.d0,__y.d3, __i , __j );	/*   x0*y3 */\
		MUL_LOHI64(__x.d1,__y.d2, __k , __l );	/*   x1*y2 */\
		MUL_LOHI64(__x.d2,__y.d1, __m , __n );	/*   x2*y1 */\
		MUL_LOHI64(__x.d3,__y.d0, __o , __p );	/*   x3*y0 */\
		\
		MUL_LOHI64(__x.d1,__y.d3, __q , __r );	/*   x1*y3 */\
		MUL_LOHI64(__x.d3,__y.d1, __s , __t );	/*   x3*y1 */\
		\
		MUL_LOHI64(__x.d2,__y.d3, __u , __v );	/*   x2*y3 */\
		MUL_LOHI64(__x.d3,__y.d2, __w , __z );	/*   x3*y2 */\
		\
		/* Now add cross terms: */\
		/* Add x0*y1 to w1-2: */\
		__w1 += __a;	__cy2  = (__w1 < __a);\
		__w2 += __b;	__cy3  = (__w2 < __b);\
		\
		/* Add x1*y0 to w1-2: */\
		__w1 += __c;	__cy2 += (__w1 < __c);\
		__w2 += __d;	__cy3 += (__w2 < __d);\
		\
		/* Add x0*y2 to w2-3: */\
		__w2 += __e;	__cy3 += (__w2 < __e);\
		__w3 += __f;	__cy4  = (__w3 < __f);\
		\
		/* Add x2*y0 to w2-3: */\
		__w2 += __g;	__cy3 += (__w2 < __g);\
		__w3 += __h;	__cy4 += (__w3 < __h);\
		\
		/* Add x0*y3 to w3-4: */\
		__w3 += __i;	__cy4 += (__w3 < __i);\
		__w4 += __j;	__cy5  = (__w4 < __j);\
		\
		/* Add x1*y2 to w3-4: */\
		__w3 += __k;	__cy4 += (__w3 < __k);\
		__w4 += __l;	__cy5 += (__w4 < __l);\
		\
		/* Add x2*y1 to w3-4: */\
		__w3 += __m;	__cy4 += (__w3 < __m);\
		__w4 += __n;	__cy5 += (__w4 < __n);\
		\
		/* Add x3*y0 to w3-4: */\
		__w3 += __o;	__cy4 += (__w3 < __o);\
		__w4 += __p;	__cy5 += (__w4 < __p);\
		\
		/* Add x1*y3 to w4-5: */\
		__w4 += __q;	__cy5 += (__w4 < __q);\
		__w5 += __r;	__cy6  = (__w5 < __r);\
		\
		/* Add x3*y1 to w4-5: */\
		__w4 += __s;	__cy5 += (__w4 < __s);\
		__w5 += __t;	__cy6 += (__w5 < __t);\
		\
		/* Add x2*y3 to w5-6: */\
		__w5 += __u;	__cy6 += (__w5 < __u);\
		__w6 += __v;	__cy7  = (__w6 < __v);\
		\
		/* Add x3*y2 to w5-6: */\
		__w5 += __w;	__cy6 += (__w5 < __w);\
		__w6 += __z;	__cy7 += (__w6 < __z);\
		\
		/* Now process carries: */\
		__w2 += __cy2;	__cy3 += (__w2 < __cy2);\
		__w3 += __cy3;	__cy4 += (__w3 < __cy3);\
		__w4 += __cy4;	__cy5 += (__w4 < __cy4);\
		__w5 += __cy5;	__cy6 += (__w5 < __cy5);\
		__w6 += __cy6;	__cy7 += (__w6 < __cy6);\
		__w7 += __cy7;\
		\
		/* Now split the result between __lo and __hi: */\
		__lo.d0 = __w0; __lo.d1 = __w1; __lo.d2 = __w2; __lo.d3 = __w3;\
		__hi.d0 = __w4;	__hi.d1 = __w5;	__hi.d2 = __w6;	__hi.d3 = __w7;\
	}

#endif


/* 512-bit square of uint256 x. Result is returned in a pair of uint256.

On Alpha, this needs a total of 20 MUL instructions and 82 ALU ops.
*/
#ifdef MUL_LOHI64_SUBROUTINE

	#define SQR_LOHI256(__x, __lo, __hi)\
	{\
		uint64 __w0, __w1, __w2, __w3, __w4, __w5, __w6, __w7\
						,__cy2,__cy3,__cy4,__cy5,__cy6,__cy7\
		,__a,__b,__e,__f,__i,__j,__k,__l\
		,__q,__r,__u,__v;\
		\
		SQR_LOHI64(__x.d0       ,&__w0,&__w1);	/*   x0^2  */\
		SQR_LOHI64(__x.d1       ,&__w2,&__w3);	/*   x1^2  */\
		SQR_LOHI64(__x.d2       ,&__w4,&__w5);	/*   x2^2  */\
		SQR_LOHI64(__x.d3       ,&__w6,&__w7);	/*   x3^2  */\
		\
		MUL_LOHI64(__x.d0,__x.d1,&__a ,&__b );	/*   x0*x1 */\
		\
		MUL_LOHI64(__x.d0,__x.d2,&__e ,&__f );	/*   x0*x2 */\
		\
		MUL_LOHI64(__x.d0,__x.d3,&__i ,&__j );	/*   x0*x3 */\
		MUL_LOHI64(__x.d1,__x.d2,&__k ,&__l );	/*   x1*x2 */\
		\
		MUL_LOHI64(__x.d1,__x.d3,&__q ,&__r );	/*   x1*x3 */\
		\
		MUL_LOHI64(__x.d2,__x.d3,&__u ,&__v );	/*   x2*x3 */\
		\
		/* Now add cross terms: */\
		/* Add x0*x1 twice to w1-2: */\
		__w1 += __a;	__cy2  = (__w1 < __a);\
		__w2 += __b;	__cy3  = (__w2 < __b);\
		\
		__w1 += __a;	__cy2 += (__w1 < __a);\
		__w2 += __b;	__cy3 += (__w2 < __b);\
		\
		/* Add x0*x2 twice to w2-3: */\
		__w2 += __e;	__cy3 += (__w2 < __e);\
		__w3 += __f;	__cy4  = (__w3 < __f);\
		\
		__w2 += __e;	__cy3 += (__w2 < __e);\
		__w3 += __f;	__cy4 += (__w3 < __f);\
		\
		/* Add x0*x3 twice to w3-4: */\
		__w3 += __i;	__cy4 += (__w3 < __i);\
		__w4 += __j;	__cy5  = (__w4 < __j);\
		\
		__w3 += __i;	__cy4 += (__w3 < __i);\
		__w4 += __j;	__cy5 += (__w4 < __j);\
		\
		/* Add x1*x2 twice to w3-4: */\
		__w3 += __k;	__cy4 += (__w3 < __k);\
		__w4 += __l;	__cy5 += (__w4 < __l);\
		\
		__w3 += __k;	__cy4 += (__w3 < __k);\
		__w4 += __l;	__cy5 += (__w4 < __l);\
		\
		/* Add x1*x3 twice to w4-5: */\
		__w4 += __q;	__cy5 += (__w4 < __q);\
		__w5 += __r;	__cy6  = (__w5 < __r);\
		\
		__w4 += __q;	__cy5 += (__w4 < __q);\
		__w5 += __r;	__cy6 += (__w5 < __r);\
		\
		/* Add x2*x3 twice to w5-6: */\
		__w5 += __u;	__cy6 += (__w5 < __u);\
		__w6 += __v;	__cy7  = (__w6 < __v);\
		\
		__w5 += __u;	__cy6 += (__w5 < __u);\
		__w6 += __v;	__cy7 += (__w6 < __v);\
		\
		/* Now process carries: */\
		__w2 += __cy2;	__cy3 += (__w2 < __cy2);\
		__w3 += __cy3;	__cy4 += (__w3 < __cy3);\
		__w4 += __cy4;	__cy5 += (__w4 < __cy4);\
		__w5 += __cy5;	__cy6 += (__w5 < __cy5);\
		__w6 += __cy6;	__cy7 += (__w6 < __cy6);\
		__w7 += __cy7;\
		\
		/* Now split the result between __lo and __hi: */\
		__lo.d0 = __w0; __lo.d1 = __w1; __lo.d2 = __w2; __lo.d3 = __w3;\
		__hi.d0 = __w4;	__hi.d1 = __w5;	__hi.d2 = __w6;	__hi.d3 = __w7;\
	}

#else

	#define SQR_LOHI256(__x, __lo, __hi)\
	{\
		uint64 __w0, __w1, __w2, __w3, __w4, __w5, __w6, __w7\
						,__cy2,__cy3,__cy4,__cy5,__cy6,__cy7\
		,__a,__b,__e,__f,__i,__j,__k,__l\
		,__q,__r,__u,__v;\
		\
		SQR_LOHI64(__x.d0       , __w0, __w1);	/*   x0^2  */\
		SQR_LOHI64(__x.d1       , __w2, __w3);	/*   x1^2  */\
		SQR_LOHI64(__x.d2       , __w4, __w5);	/*   x2^2  */\
		SQR_LOHI64(__x.d3       , __w6, __w7);	/*   x3^2  */\
		\
		MUL_LOHI64(__x.d0,__x.d1, __a , __b );	/*   x0*x1 */\
		\
		MUL_LOHI64(__x.d0,__x.d2, __e , __f );	/*   x0*x2 */\
		\
		MUL_LOHI64(__x.d0,__x.d3, __i , __j );	/*   x0*x3 */\
		MUL_LOHI64(__x.d1,__x.d2, __k , __l );	/*   x1*x2 */\
		\
		MUL_LOHI64(__x.d1,__x.d3, __q , __r );	/*   x1*x3 */\
		\
		MUL_LOHI64(__x.d2,__x.d3, __u , __v );	/*   x2*x3 */\
		\
		/* Now add cross terms: */\
		/* Add x0*x1 twice to w1-2: */\
		__w1 += __a;	__cy2  = (__w1 < __a);\
		__w2 += __b;	__cy3  = (__w2 < __b);\
		\
		__w1 += __a;	__cy2 += (__w1 < __a);\
		__w2 += __b;	__cy3 += (__w2 < __b);\
		\
		/* Add x0*x2 twice to w2-3: */\
		__w2 += __e;	__cy3 += (__w2 < __e);\
		__w3 += __f;	__cy4  = (__w3 < __f);\
		\
		__w2 += __e;	__cy3 += (__w2 < __e);\
		__w3 += __f;	__cy4 += (__w3 < __f);\
		\
		/* Add x0*x3 twice to w3-4: */\
		__w3 += __i;	__cy4 += (__w3 < __i);\
		__w4 += __j;	__cy5  = (__w4 < __j);\
		\
		__w3 += __i;	__cy4 += (__w3 < __i);\
		__w4 += __j;	__cy5 += (__w4 < __j);\
		\
		/* Add x1*x2 twice to w3-4: */\
		__w3 += __k;	__cy4 += (__w3 < __k);\
		__w4 += __l;	__cy5 += (__w4 < __l);\
		\
		__w3 += __k;	__cy4 += (__w3 < __k);\
		__w4 += __l;	__cy5 += (__w4 < __l);\
		\
		/* Add x1*x3 twice to w4-5: */\
		__w4 += __q;	__cy5 += (__w4 < __q);\
		__w5 += __r;	__cy6  = (__w5 < __r);\
		\
		__w4 += __q;	__cy5 += (__w4 < __q);\
		__w5 += __r;	__cy6 += (__w5 < __r);\
		\
		/* Add x2*x3 twice to w5-6: */\
		__w5 += __u;	__cy6 += (__w5 < __u);\
		__w6 += __v;	__cy7  = (__w6 < __v);\
		\
		__w5 += __u;	__cy6 += (__w5 < __u);\
		__w6 += __v;	__cy7 += (__w6 < __v);\
		\
		/* Now process carries: */\
		__w2 += __cy2;	__cy3 += (__w2 < __cy2);\
		__w3 += __cy3;	__cy4 += (__w3 < __cy3);\
		__w4 += __cy4;	__cy5 += (__w4 < __cy4);\
		__w5 += __cy5;	__cy6 += (__w5 < __cy5);\
		__w6 += __cy6;	__cy7 += (__w6 < __cy6);\
		__w7 += __cy7;\
		\
		/* Now split the result between __lo and __hi: */\
		__lo.d0 = __w0; __lo.d1 = __w1; __lo.d2 = __w2; __lo.d3 = __w3;\
		__hi.d0 = __w4;	__hi.d1 = __w5;	__hi.d2 = __w6;	__hi.d3 = __w7;\
	}

#endif

	/* 4-operand-pipelined version: */
	#define SQR_LOHI256_q4(\
	  __x0, __lo0, __hi0\
	, __x1, __lo1, __hi1\
	, __x2, __lo2, __hi2\
	, __x3, __lo3, __hi3)\
	{\
		SQR_LOHI256(__x0, __lo0, __hi0);\
		SQR_LOHI256(__x1, __lo1, __hi1);\
		SQR_LOHI256(__x2, __lo2, __hi2);\
		SQR_LOHI256(__x3, __lo3, __hi3);\
	}

	/* 8-operand-pipelined version: */
	#define SQR_LOHI256_q8(\
	  __x0, __lo0, __hi0\
	, __x1, __lo1, __hi1\
	, __x2, __lo2, __hi2\
	, __x3, __lo3, __hi3\
	, __x4, __lo4, __hi4\
	, __x5, __lo5, __hi5\
	, __x6, __lo6, __hi6\
	, __x7, __lo7, __hi7)\
	{\
		SQR_LOHI256(__x0, __lo0, __hi0);\
		SQR_LOHI256(__x1, __lo1, __hi1);\
		SQR_LOHI256(__x2, __lo2, __hi2);\
		SQR_LOHI256(__x3, __lo3, __hi3);\
		SQR_LOHI256(__x4, __lo4, __hi4);\
		SQR_LOHI256(__x5, __lo5, __hi5);\
		SQR_LOHI256(__x6, __lo6, __hi6);\
		SQR_LOHI256(__x7, __lo7, __hi7);\
	}

/* Lower 256 bits of the product of uint256 x and y.

On Alpha, this needs a total of 10 MUL instructions and 26 ALU ops.
*/
#ifdef MUL_LOHI64_SUBROUTINE

	#define MULL256(__x, __y, __lo)\
	{\
		uint64 __w0, __w1, __w2, __w3\
						,__cy2,__cy3\
		,__a,__b,__c,__d,__e,__f,__g,__h,__i,__j,__k,__l;\
		\
		MUL_LOHI64(__x.d0,__y.d0,&__w0,&__w1);	/*   x0*y0 */\
		MUL_LOHI64(__x.d1,__y.d1,&__w2,&__w3);	/*   x1*y1 */\
		\
		MUL_LOHI64(__x.d0,__y.d1,&__a ,&__b );	/*   x0*y1 */\
		MUL_LOHI64(__x.d1,__y.d0,&__c ,&__d );	/*   x1*y0 */\
		\
		MUL_LOHI64(__x.d0,__y.d2,&__e ,&__f );	/*   x0*y2 */\
		MUL_LOHI64(__x.d2,__y.d0,&__g ,&__h );	/*   x2*y0 */\
		\
		__i = __x.d0*__y.d3;				/* (x0*y3).lo */\
		__j = __x.d1*__y.d2;				/* (x1*y2).lo */\
		__k = __x.d2*__y.d1;				/* (x2*y1).lo */\
		__l = __x.d3*__y.d0;				/* (x3*y0).lo */\
		\
		/* Now add cross terms: */\
		/* Add x0*y1 to w1-2: */\
		__w1 += __a;	__cy2  = (__w1 < __a);\
		__w2 += __b;	__cy3  = (__w2 < __b);\
		\
		/* Add x1*y0 to w1-2: */\
		__w1 += __c;	__cy2 += (__w1 < __c);\
		__w2 += __d;	__cy3 += (__w2 < __d);\
		\
		/* Add x0*y2 to w2-3: */\
		__w2 += __e;	__cy3 += (__w2 < __e);\
		__w3 += __f;\
		\
		/* Add x2*y0 to w2-3: */\
		__w2 += __g;	__cy3 += (__w2 < __g);\
		__w3 += __h;\
		\
		/* Add (x0*y3 + x1*y2 + x2*y1 + x3*y0).lo to w3: */\
		__w3 += __i+__j+__k+__l;\
		\
		/* Now process carries: */\
		__w2 += __cy2;	__cy3 += (__w2 < __cy2);\
		__w3 += __cy3;\
		\
		/* Now return the result in __lo: */\
		__lo.d0 = __w0; __lo.d1 = __w1; __lo.d2 = __w2; __lo.d3 = __w3;\
	}

#else

	#define MULL256(__x, __y, __lo)\
	{\
		uint64 __w0, __w1, __w2, __w3\
						,__cy2,__cy3\
		,__a,__b,__c,__d,__e,__f,__g,__h,__i,__j,__k,__l;\
		\
		MUL_LOHI64(__x.d0,__y.d0, __w0, __w1);	/*   x0*y0 */\
		MUL_LOHI64(__x.d1,__y.d1, __w2, __w3);	/*   x1*y1 */\
		\
		MUL_LOHI64(__x.d0,__y.d1, __a , __b );	/*   x0*y1 */\
		MUL_LOHI64(__x.d1,__y.d0, __c , __d );	/*   x1*y0 */\
		\
		MUL_LOHI64(__x.d0,__y.d2, __e , __f );	/*   x0*y2 */\
		MUL_LOHI64(__x.d2,__y.d0, __g , __h );	/*   x2*y0 */\
		\
		__i = __x.d0*__y.d3;				/* (x0*y3).lo */\
		__j = __x.d1*__y.d2;				/* (x1*y2).lo */\
		__k = __x.d2*__y.d1;				/* (x2*y1).lo */\
		__l = __x.d3*__y.d0;				/* (x3*y0).lo */\
		\
		/* Now add cross terms: */\
		/* Add x0*y1 to w1-2: */\
		__w1 += __a;	__cy2  = (__w1 < __a);\
		__w2 += __b;	__cy3  = (__w2 < __b);\
		\
		/* Add x1*y0 to w1-2: */\
		__w1 += __c;	__cy2 += (__w1 < __c);\
		__w2 += __d;	__cy3 += (__w2 < __d);\
		\
		/* Add x0*y2 to w2-3: */\
		__w2 += __e;	__cy3 += (__w2 < __e);\
		__w3 += __f;\
		\
		/* Add x2*y0 to w2-3: */\
		__w2 += __g;	__cy3 += (__w2 < __g);\
		__w3 += __h;\
		\
		/* Add (x0*y3 + x1*y2 + x2*y1 + x3*y0).lo to w3: */\
		__w3 += __i+__j+__k+__l;\
		\
		/* Now process carries: */\
		__w2 += __cy2;	__cy3 += (__w2 < __cy2);\
		__w3 += __cy3;\
		\
		/* Now return the result in __lo: */\
		__lo.d0 = __w0; __lo.d1 = __w1; __lo.d2 = __w2; __lo.d3 = __w3;\
	}

#endif

	/* 4-operand-pipelined version: */
	#define MULL256_q4(\
	  __x0, __y0, __lo0\
	, __x1, __y1, __lo1\
	, __x2, __y2, __lo2\
	, __x3, __y3, __lo3)\
	{\
		MULL256(__x0, __y0, __lo0);\
		MULL256(__x1, __y1, __lo1);\
		MULL256(__x2, __y2, __lo2);\
		MULL256(__x3, __y3, __lo3);\
	}

	/* 8-operand-pipelined version: */
	#define MULL256_q8(\
	  __x0, __y0, __lo0\
	, __x1, __y1, __lo1\
	, __x2, __y2, __lo2\
	, __x3, __y3, __lo3\
	, __x4, __y4, __lo4\
	, __x5, __y5, __lo5\
	, __x6, __y6, __lo6\
	, __x7, __y7, __lo7)\
	{\
		MULL256(__x0, __y0, __lo0);\
		MULL256(__x1, __y1, __lo1);\
		MULL256(__x2, __y2, __lo2);\
		MULL256(__x3, __y3, __lo3);\
		MULL256(__x4, __y4, __lo4);\
		MULL256(__x5, __y5, __lo5);\
		MULL256(__x6, __y6, __lo6);\
		MULL256(__x7, __y7, __lo7);\
	}

/* Upper 256 bits of the product of uint256 x and y.

On Alpha, this needs a total of 32 MUL, 82 ALU ops.
*/
#ifdef MUL_LOHI64_SUBROUTINE

	#define MULH256(__x, __y, __hi)\
	{\
		uint256 __lo;\
		MUL_LOHI256(__x, __y, __lo, __hi);\
	}

#else

	#define MULH256(__x, __y, __hi)\
	{\
		uint64 __w0, __w1, __w2, __w3, __w4, __w5, __w6, __w7\
						,__cy2,__cy3,__cy4,__cy5,__cy6,__cy7\
		,__a,__b,__c,__d,__e,__f,__g,__h,__i,__j,__k,__l\
		,__m,__n,__o,__p,__q,__r,__s,__t,__u,__v,__w,__z;\
		\
		MUL_LOHI64(__x.d0,__y.d0, __w0, __w1);	/*   x0*y0 */\
		MUL_LOHI64(__x.d1,__y.d1, __w2, __w3);	/*   x1*y1 */\
		MUL_LOHI64(__x.d2,__y.d2, __w4, __w5);	/*   x2*y2 */\
		MUL_LOHI64(__x.d3,__y.d3, __w6, __w7);	/*   x3*y3 */\
		\
		MUL_LOHI64(__x.d0,__y.d1, __a , __b );	/*   x0*y1 */\
		MUL_LOHI64(__x.d1,__y.d0, __c , __d );	/*   x1*y0 */\
		\
		MUL_LOHI64(__x.d0,__y.d2, __e , __f );	/*   x0*y2 */\
		MUL_LOHI64(__x.d2,__y.d0, __g , __h );	/*   x2*y0 */\
		\
		MUL_LOHI64(__x.d0,__y.d3, __i , __j );	/*   x0*y3 */\
		MUL_LOHI64(__x.d1,__y.d2, __k , __l );	/*   x1*y2 */\
		MUL_LOHI64(__x.d2,__y.d1, __m , __n );	/*   x2*y1 */\
		MUL_LOHI64(__x.d3,__y.d0, __o , __p );	/*   x3*y0 */\
		\
		MUL_LOHI64(__x.d1,__y.d3, __q , __r );	/*   x1*y3 */\
		MUL_LOHI64(__x.d3,__y.d1, __s , __t );	/*   x3*y1 */\
		\
		MUL_LOHI64(__x.d2,__y.d3, __u , __v );	/*   x2*y3 */\
		MUL_LOHI64(__x.d3,__y.d2, __w , __z );	/*   x3*y2 */\
		\
		/* Now add cross terms: */\
		/* Add x0*y1 to w1-2: */\
		__w1 += __a;	__cy2  = (__w1 < __a);\
		__w2 += __b;	__cy3  = (__w2 < __b);\
		\
		/* Add x1*y0 to w1-2: */\
		__w1 += __c;	__cy2 += (__w1 < __c);\
		__w2 += __d;	__cy3 += (__w2 < __d);\
		\
		/* Add x0*y2 to w2-3: */\
		__w2 += __e;	__cy3 += (__w2 < __e);\
		__w3 += __f;	__cy4  = (__w3 < __f);\
		\
		/* Add x2*y0 to w2-3: */\
		__w2 += __g;	__cy3 += (__w2 < __g);\
		__w3 += __h;	__cy4 += (__w3 < __h);\
		\
		/* Add x0*y3 to w3-4: */\
		__w3 += __i;	__cy4 += (__w3 < __i);\
		__w4 += __j;	__cy5  = (__w4 < __j);\
		\
		/* Add x1*y2 to w3-4: */\
		__w3 += __k;	__cy4 += (__w3 < __k);\
		__w4 += __l;	__cy5 += (__w4 < __l);\
		\
		/* Add x2*y1 to w3-4: */\
		__w3 += __m;	__cy4 += (__w3 < __m);\
		__w4 += __n;	__cy5 += (__w4 < __n);\
		\
		/* Add x3*y0 to w3-4: */\
		__w3 += __o;	__cy4 += (__w3 < __o);\
		__w4 += __p;	__cy5 += (__w4 < __p);\
		\
		/* Add x1*y3 to w4-5: */\
		__w4 += __q;	__cy5 += (__w4 < __q);\
		__w5 += __r;	__cy6  = (__w5 < __r);\
		\
		/* Add x3*y1 to w4-5: */\
		__w4 += __s;	__cy5 += (__w4 < __s);\
		__w5 += __t;	__cy6 += (__w5 < __t);\
		\
		/* Add x2*y3 to w5-6: */\
		__w5 += __u;	__cy6 += (__w5 < __u);\
		__w6 += __v;	__cy7  = (__w6 < __v);\
		\
		/* Add x3*y2 to w5-6: */\
		__w5 += __w;	__cy6 += (__w5 < __w);\
		__w6 += __z;	__cy7 += (__w6 < __z);\
		\
		/* Now process carries: */\
		__w2 += __cy2;	__cy3 += (__w2 < __cy2);\
		__w3 += __cy3;	__cy4 += (__w3 < __cy3);\
		__w4 += __cy4;	__cy5 += (__w4 < __cy4);\
		__w5 += __cy5;	__cy6 += (__w5 < __cy5);\
		__w6 += __cy6;	__cy7 += (__w6 < __cy6);\
		__w7 += __cy7;\
		\
		/* Now return the upper half of the product in __hi: */\
		__hi.d0 = __w4;	__hi.d1 = __w5;	__hi.d2 = __w6;	__hi.d3 = __w7;\
	}

#endif

	/* 4-operand-pipelined version: */
	#define MULH256_q4(\
	  __x0, __y0, __hi0\
	, __x1, __y1, __hi1\
	, __x2, __y2, __hi2\
	, __x3, __y3, __hi3)\
	{\
		MULH256(__x0, __y0, __hi0);\
		MULH256(__x1, __y1, __hi1);\
		MULH256(__x2, __y2, __hi2);\
		MULH256(__x3, __y3, __hi3);\
	}

	/* 8-operand-pipelined version: */
	#define MULH256_q8(\
	  __x0, __y0, __hi0\
	, __x1, __y1, __hi1\
	, __x2, __y2, __hi2\
	, __x3, __y3, __hi3\
	, __x4, __y4, __hi4\
	, __x5, __y5, __hi5\
	, __x6, __y6, __hi6\
	, __x7, __y7, __hi7)\
	{\
		MULH256(__x0, __y0, __hi0);\
		MULH256(__x1, __y1, __hi1);\
		MULH256(__x2, __y2, __hi2);\
		MULH256(__x3, __y3, __hi3);\
		MULH256(__x4, __y4, __hi4);\
		MULH256(__x5, __y5, __hi5);\
		MULH256(__x6, __y6, __hi6);\
		MULH256(__x7, __y7, __hi7);\
	}

#ifdef __cplusplus
}
#endif

#endif	/* imul_macro256_h_included */

