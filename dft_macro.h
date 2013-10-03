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

/****************************************************************************
 * We now include this header file if it was not included before.
 ****************************************************************************/
#ifndef dft_macro_included
#define dft_macro_included

/* All these radix-N macros are structured so that any 2 of the 3 sets of arguments (A, t, B) may point to the same set of
   addresses, i.e. so the DFT may be done in place (A == B) or re-using a single set of temporaries (A == t or B == t).
   By "in place" we mean that A == B up to an arbitrary permutation of array indices (or element subscripts). */

/****** RADIX = 3: ALLOWS IN-PLACE ******/

/* Totals: 12 ADD, 4 MUL	*/
#define RADIX_03_DFT(\
	__s,__c3m1,\
	__Ar0,__Ai0,\
	__Ar1,__Ai1,\
	__Ar2,__Ai2,\
	__tr0,__ti0,\
	__tr1,__ti1,\
	__tr2,__ti2,\
	__Br0,__Bi0,\
	__Br1,__Bi1,\
	__Br2,__Bi2)\
{\
	__tr0 =__Ar2;					__ti0 =__Ai2;\
	__tr2 =__Ar1 - __tr0;			__ti2 =__Ai1 - __ti0;\
	__tr1 =__Ar1 + __tr0;			__ti1 =__Ai1 + __ti0;\
	__Br0 =__Ar0 + __tr1;			__Bi0 =__Ai0 + __ti1;\
	__Br1 =__Br0 + __tr1*__c3m1;	__Bi1 =__Bi0 + __ti1*__c3m1;\
	__tr0 =__s*__tr2;				__ti0 =__s*__ti2;\
	__Br2 =__Br1 + __ti0;			__Bi2 =__Bi1 - __tr0;\
	__Br1 =__Br1 - __ti0;			__Bi1 =__Bi1 + __tr0;\
}

/* Totals: 12 ADD, 4 MUL	*/
#define RADIX_03_DFT_PFETCH(\
	__s,__c3m1,\
	__Ar0,__Ai0,\
	__Ar1,__Ai1,\
	__Ar2,__Ai2,\
	__tr0,__ti0,\
	__tr1,__ti1,\
	__tr2,__ti2,\
	__Br0,__Bi0,\
	__Br1,__Bi1,\
	__Br2,__Bi2,\
	__add_offset)\
{\
	__tr0 =__Ar2;					__ti0 =__Ai2;\
	__tr2 =__Ar1 - __tr0;			__ti2 =__Ai1 - __ti0;\
	__tr1 =__Ar1 + __tr0;			__ti1 =__Ai1 + __ti0;\
	__Br0 =__Ar0 + __tr1;			__Bi0 =__Ai0 + __ti1;\
	__Br1 =__Br0 + __tr1*__c3m1;	__Bi1 =__Bi0 + __ti1*__c3m1;\
\
	addr = add0 + __add_offset;	prefetch_p_doubles(addr);\
\
	__tr0 =__s*__tr2;				__ti0 =__s*__ti2;\
	__Br2 =__Br1 + __ti0;			__Bi2 =__Bi1 - __tr0;\
	__Br1 =__Br1 - __ti0;			__Bi1 =__Bi1 + __tr0;\
}

/****** RADIX = 4: ******/

/* Totals: 16 ADD, 0 MUL	*/
#define RADIX_04_DIF(\
	__Ar0,__Ai0,\
	__Ar1,__Ai1,\
	__Ar2,__Ai2,\
	__Ar3,__Ai3,\
	__Br0,__Bi0,\
	__Br1,__Bi1,\
	__Br2,__Bi2,\
	__Br3,__Bi3,\
	__rt,__it)\
{\
	__rt = __Ar2;	__Ar2 =__Ar0 - __rt;	__Ar0 = __Ar0 + __rt;\
	__it = __Ai2;	__Ai2 =__Ai0 - __it;	__Ai0 = __Ai0 + __it;\
	__rt = __Ar3;	__Ar3 =__Ar1 - __rt;	__Ar1 = __Ar1 + __rt;\
	__it = __Ai3;	__Ai3 =__Ai1 - __it;	__Ai1 = __Ai1 + __it;\
\
	__rt  = __Ar1;	__Br1 = __Ar0 - __rt;	__Br0 = __Ar0 + __rt;\
	__it  = __Ai1;	__Bi1 = __Ai0 - __it;	__Bi0 = __Ai0 + __it;\
	/* mpy by I is inlined here...	*/\
	__rt  = __Ar3;	__Br3 = __Ar2 + __Ai3;	__Br2 = __Ar2 - __Ai3;\
			__Bi3 = __Ai2 - __rt;			__Bi2 = __Ai2 + __rt;\
}

/* Totals: 16 ADD, 0 MUL	*/
#define RADIX_04_DIF_PFETCH(\
	__Ar0,__Ai0,\
	__Ar1,__Ai1,\
	__Ar2,__Ai2,\
	__Ar3,__Ai3,\
	__Br0,__Bi0,\
	__Br1,__Bi1,\
	__Br2,__Bi2,\
	__Br3,__Bi3,\
	__rt,__it,\
	__addr,__addp,__add_offset0,__add_offset1)\
{\
	__rt = __Ar2;	__Ar2 =__Ar0 - __rt;	__Ar0 = __Ar0 + __rt;\
	__it = __Ai2;	__Ai2 =__Ai0 - __it;	__Ai0 = __Ai0 + __it;\
	__rt = __Ar3;	__Ar3 =__Ar1 - __rt;	__Ar1 = __Ar1 + __rt;\
	__it = __Ai3;	__Ai3 =__Ai1 - __it;	__Ai1 = __Ai1 + __it;\
\
	__addp = __addr + __add_offset0;	prefetch_p_doubles(__addp);\
\
	__rt  = __Ar1;	__Br1 = __Ar0 - __rt;	__Br0 = __Ar0 + __rt;\
	__it  = __Ai1;	__Bi1 = __Ai0 - __it;	__Bi0 = __Ai0 + __it;\
	/* mpy by I is inlined here...	*/\
	__rt  = __Ar3;	__Br3 = __Ar2 + __Ai3;	__Br2 = __Ar2 - __Ai3;\
			__Bi3 = __Ai2 - __rt;			__Bi2 = __Ai2 + __rt;\
\
	__addp = __addr + __add_offset1;	prefetch_p_doubles(__addp);\
}

/* Totals: 16 ADD, 0 MUL	*/
#define RADIX_04_DIT(\
	__Ar0,__Ai0,\
	__Ar1,__Ai1,\
	__Ar2,__Ai2,\
	__Ar3,__Ai3,\
	__Br0,__Bi0,\
	__Br1,__Bi1,\
	__Br2,__Bi2,\
	__Br3,__Bi3,\
	__rt,__it)\
{\
	__rt = __Ar1;	__Ar1 =__Ar0 - __rt;	__Ar0 = __Ar0 + __rt;\
	__it = __Ai1;	__Ai1 =__Ai0 - __it;	__Ai0 = __Ai0 + __it;\
	__rt = __Ar3;	__Ar3 =__Ar2 - __rt;	__Ar2 = __Ar2 + __rt;\
	__it = __Ai3;	__Ai3 =__Ai2 - __it;	__Ai2 = __Ai2 + __it;\
\
	__rt  = __Ar2;	__Br2 = __Ar0 - __rt;	__Br0 = __Ar0 + __rt;\
	__it  = __Ai2;	__Bi2 = __Ai0 - __it;	__Bi0 = __Ai0 + __it;\
	/* mpy by I is inlined here...	*/\
	__rt  = __Ar3;	__Br3 = __Ar1 - __Ai3;	__Br1 = __Ar1 + __Ai3;\
					__Bi3 = __Ai1 + __rt;	__Bi1 = __Ai1 - __rt;\
}

/****** RADIX = 5: ******/

/* Totals: 34 ADD, 10 MUL	*/
#define RADIX_05_DFT(\
	__cc1, __cc2, __s2, __ss1, __ss2,\
	__Ar0,__Ai0,\
	__Ar1,__Ai1,\
	__Ar2,__Ai2,\
	__Ar3,__Ai3,\
	__Ar4,__Ai4,\
	__Br0,__Bi0,\
	__Br1,__Bi1,\
	__Br2,__Bi2,\
	__Br3,__Bi3,\
	__Br4,__Bi4,\
	__rt,__it)\
{\
	__Br0 = __Ar0;					__Bi0 = __Ai0;\
	__rt  = __Ar4;					__it  = __Ai4;\
	__Br4 = __Ar1 - __rt;			__Bi4 = __Ai1 - __it;\
	__Br1 = __Ar1 + __rt;			__Bi1 = __Ai1 + __it;\
	__rt  = __Ar3;					__it  = __Ai3;\
	__Br3 = __Ar2 - __rt;			__Bi3 = __Ai2 - __it;\
	__Br2 = __Ar2 + __rt;			__Bi2 = __Ai2 + __it;\
	__rt  = __Br1 + __Br2;			__it  = __Bi1 + __Bi2;\
	__Br0 = __Br0 + __rt;			__Bi0 = __Bi0 + __it;		/* y0 */\
	__rt  = __Br0 + __rt  *__cc1;	__it  = __Bi0 + __it  *__cc1;\
	__Br2 =(__Br1 - __Br2)*__cc2;	__Bi2 =(__Bi1 - __Bi2)*__cc2;\
	__Br1 = __rt  + __Br2;			__Bi1 = __it  + __Bi2;\
	__Br2 = __rt  - __Br2;			__Bi2 = __it  - __Bi2;\
	__rt  =(__Br4 - __Br3)* __s2;	__it  =(__Bi4 - __Bi3)* __s2;\
	__Br3 = __rt  + __Br3 *__ss1;	__Bi3 = __it  + __Bi3 *__ss1;\
	__Br4 = __rt  - __Br4 *__ss2;	__Bi4 = __it  - __Bi4 *__ss2;\
	__rt  = __Br4;					__it  = __Bi4;\
	__Br4 = __Br1 + __Bi3;			__Bi4 = __Bi1 - __Br3;		/* y4 - note swap with y3 below! */\
	__Br1 = __Br1 - __Bi3;			__Bi1 = __Bi1 + __Br3;		/* y1 */\
	__Br3 = __Br2 + __it;			__Bi3 = __Bi2 - __rt;		/* y3 */\
	__Br2 = __Br2 - __it;			__Bi2 = __Bi2 + __rt;		/* y2 */\
}

/* Totals: 34 ADD, 10 MUL	*/
#define RADIX_05_DFT_PFETCH(\
	__cc1, __cc2, __s2, __ss1, __ss2,\
	__Ar0,__Ai0,\
	__Ar1,__Ai1,\
	__Ar2,__Ai2,\
	__Ar3,__Ai3,\
	__Ar4,__Ai4,\
	__Br0,__Bi0,\
	__Br1,__Bi1,\
	__Br2,__Bi2,\
	__Br3,__Bi3,\
	__Br4,__Bi4,\
	__rt,__it,\
	__addr,__addp,__add_offset0,__add_offset1)\
{\
	__Br0 = __Ar0;					__Bi0 = __Ai0;\
	__rt  = __Ar4;					__it  = __Ai4;\
	__Br4 = __Ar1 - __rt;			__Bi4 = __Ai1 - __it;\
	__Br1 = __Ar1 + __rt;			__Bi1 = __Ai1 + __it;\
	__rt  = __Ar3;					__it  = __Ai3;\
	__Br3 = __Ar2 - __rt;			__Bi3 = __Ai2 - __it;\
	__Br2 = __Ar2 + __rt;			__Bi2 = __Ai2 + __it;\
	__rt  = __Br1 + __Br2;			__it  = __Bi1 + __Bi2;\
	__Br0 = __Br0 + __rt;			__Bi0 = __Bi0 + __it;		/* y0 */\
	__rt  = __Br0 + __rt  *__cc1;	__it  = __Bi0 + __it  *__cc1;\
	__Br2 =(__Br1 - __Br2)*__cc2;	__Bi2 =(__Bi1 - __Bi2)*__cc2;\
\
	__addp = __addr + __add_offset0;	prefetch_p_doubles(__addp);\
\
	__Br1 = __rt  + __Br2;			__Bi1 = __it  + __Bi2;\
	__Br2 = __rt  - __Br2;			__Bi2 = __it  - __Bi2;\
	__rt  =(__Br4 - __Br3)* __s2;	__it  =(__Bi4 - __Bi3)* __s2;\
	__Br3 = __rt  + __Br3 *__ss1;	__Bi3 = __it  + __Bi3 *__ss1;\
	__Br4 = __rt  - __Br4 *__ss2;	__Bi4 = __it  - __Bi4 *__ss2;\
	__rt  = __Br4;					__it  = __Bi4;\
	__Br4 = __Br1 + __Bi3;			__Bi4 = __Bi1 - __Br3;		/* y4 - note swap with y3 below! */\
	__Br1 = __Br1 - __Bi3;			__Bi1 = __Bi1 + __Br3;		/* y1 */\
	__Br3 = __Br2 + __it;			__Bi3 = __Bi2 - __rt;		/* y3 */\
	__Br2 = __Br2 - __it;			__Bi2 = __Bi2 + __rt;		/* y2 */\
\
	__addp = __addr + __add_offset1;	prefetch_p_doubles(__addp);\
}

/****** RADIX = 7: ******/

/* Simple version may be best if ADD throughput is limiting factor.
Totals: 60 ADD, 36 MUL	*/
#define RADIX_07_DFT(\
	__Ar0,__Ai0,\
	__Ar1,__Ai1,\
	__Ar2,__Ai2,\
	__Ar3,__Ai3,\
	__Ar4,__Ai4,\
	__Ar5,__Ai5,\
	__Ar6,__Ai6,\
	__tr0,__ti0,\
	__tr1,__ti1,\
	__tr2,__ti2,\
	__tr3,__ti3,\
	__tr4,__ti4,\
	__tr5,__ti5,\
	__tr6,__ti6,\
	__Br0,__Bi0,\
	__Br1,__Bi1,\
	__Br2,__Bi2,\
	__Br3,__Bi3,\
	__Br4,__Bi4,\
	__Br5,__Bi5,\
	__Br6,__Bi6,\
	__cc1,__ss1,\
	__cc2,__ss2,\
	__cc3,__ss3,\
	__rt,__it,\
	__re,__im)\
{\
	__tr0 = __Ar0;												__ti0 = __Ai0;\
	__tr6 = __Ar1 - __Ar6;										__ti6 = __Ai1 - __Ai6;	/* x1 - x6 */\
	__tr1 = __Ar1 + __Ar6;										__ti1 = __Ai1 + __Ai6;	/* x1 + x6 */\
\
	__tr5 = __Ar2 - __Ar5;										__ti5 = __Ai2 - __Ai5;	/* x2 - x5 */\
	__tr2 = __Ar2 + __Ar5;										__ti2 = __Ai2 + __Ai5;	/* x2 + x5 */\
\
	__tr4 = __Ar3 - __Ar4;										__ti4 = __Ai3 - __Ai4;	/* x3 - x4 */\
	__tr3 = __Ar3 + __Ar4;										__ti3 = __Ai3 + __Ai4;	/* x3 + x4 */\
\
	__Br0 = __tr0 + __tr1 + __tr2 + __tr3;						__Bi0 = __ti0 + __ti1 + __ti2 + __ti3;\
	__rt  = __tr0 + __cc1*__tr1 + __cc2*__tr2 + __cc3*__tr3;	__it  = __ti0 + __cc1*__ti1 + __cc2*__ti2 + __cc3*__ti3;\
	__re  = __tr0 + __cc2*__tr1 + __cc3*__tr2 + __cc1*__tr3;	__im  = __ti0 + __cc2*__ti1 + __cc3*__ti2 + __cc1*__ti3;\
	__tr0 = __tr0 + __cc3*__tr1 + __cc1*__tr2 + __cc2*__tr3;	__ti0 = __ti0 + __cc3*__ti1 + __cc1*__ti2 + __cc2*__ti3;\
\
	__tr1 = __ss1*__tr6 + __ss2*__tr5 + __ss3*__tr4;			__ti1 = __ss1*__ti6 + __ss2*__ti5 + __ss3*__ti4;\
	__tr2 = __ss2*__tr6 - __ss3*__tr5 - __ss1*__tr4;			__ti2 = __ss2*__ti6 - __ss3*__ti5 - __ss1*__ti4;\
	__tr3 = __ss3*__tr6 - __ss1*__tr5 + __ss2*__tr4;			__ti3 = __ss3*__ti6 - __ss1*__ti5 + __ss2*__ti4;\
	/* Output permutation causes signs to get flipped here: */\
	__Br1 = __rt  - __ti1;										__Bi1 = __it  + __tr1;\
	__Br2 = __re  - __ti2;										__Bi2 = __im  + __tr2;\
	__Br3 = __tr0 - __ti3;										__Bi3 = __ti0 + __tr3;\
	__Br4 = __tr0 + __ti3;										__Bi4 = __ti0 - __tr3;\
	__Br5 = __re  + __ti2;										__Bi5 = __im  - __tr2;\
	__Br6 = __rt  + __ti1;										__Bi6 = __it  - __tr1;\
}

/* Low-mul Nussbaumer-style DFT implementation has 20 fewer MUL but at cost of 12 more ADD.
Totals: 72 ADD, 16 MUL	*/
#define RADIX_07_DFT_NUSS(\
	__A0r,__A0i,\
	__A1r,__A1i,\
	__A2r,__A2i,\
	__A3r,__A3i,\
	__A4r,__A4i,\
	__A5r,__A5i,\
	__A6r,__A6i,\
	__t0r,__t0i,\
	__t1r,__t1i,\
	__t2r,__t2i,\
	__t3r,__t3i,\
	__t4r,__t4i,\
	__t5r,__t5i,\
	__t6r,__t6i,\
	__B0r,__B0i,\
	__B1r,__B1i,\
	__B2r,__B2i,\
	__B3r,__B3i,\
	__B4r,__B4i,\
	__B5r,__B5i,\
	__B6r,__B6i,\
	__cx0,__sx0,\
	__cx1,__sx1,\
	__cx2,__sx2,\
	__cx3,__sx3,\
	__rt,__it)\
{\
	__t0r = __A0r;						__t0i = __A0i;\
	__t6r = __A1r - __A6r;				__t6i = __A1i - __A6i;	/* x1 - x6 */\
	__t1r = __A1r + __A6r;				__t1i = __A1i + __A6i;	/* x1 + x6 */\
\
	__t5r = __A2r - __A5r;				__t5i = __A2i - __A5i;	/* x2 - x5 */\
	__t2r = __A2r + __A5r;				__t2i = __A2i + __A5i;	/* x2 + x5 */\
\
	__t4r = __A3r - __A4r;				__t4i = __A3i - __A4i;	/* x3 - x4 */\
	__t3r = __A3r + __A4r;				__t3i = __A3i + __A4i;	/* x3 + x4 */\
\
	__rt = __t1r+__t2r+__t3r;			__it = __t1i+__t2i+__t3i;\
	__B0r= __rt+__t0r;					__B0i= __it+__t0i;\
	__t0r= __rt*__cx0+__t0r;			__t0i= __it*__cx0+__t0i;\
	__t1r= __t1r-__t2r;					__t1i= __t1i-__t2i;\
	__t2r= __t3r-__t2r;					__t2i= __t3i-__t2i;\
	__t3r=(__t1r+__t2r)*__cx3;			__t3i=(__t1i+__t2i)*__cx3;\
	__t1r= __t1r*__cx1;					__t1i= __t1i*__cx1;\
	__t2r= __t2r*__cx2;					__t2i= __t2i*__cx2;\
	__rt = __t1r-__t3r;					__it = __t1i-__t3i;\
	__t2r= __t2r-__t3r;					__t2i= __t2i-__t3i;\
\
	__t1r= __t0r-__rt-__t2r;			__t1i= __t0i-__it-__t2i;\
	__t2r= __t0r+__t2r;					__t2i= __t0i+__t2i;\
	__t0r= __t0r+__rt;					__t0i= __t0i+__it;\
\
	__t3r=(__t6r-__t4r+__t5r)*__sx0;	__t3i=(__t6i-__t4i+__t5i)*__sx0;\
	__t6r= __t6r-__t5r;					__t6i= __t6i-__t5i;\
	__t5r= __t4r+__t5r;					__t5i= __t4i+__t5i;\
	__t4r=(__t5r-__t6r)*__sx3;			__t4i=(__t5i-__t6i)*__sx3;\
	__t6r= __t6r*__sx1;					__t6i= __t6i*__sx1;\
	__t5r= __t5r*__sx2;					__t5i= __t5i*__sx2;\
	__t6r= __t4r+__t6r;					__t6i= __t4i+__t6i;\
	__t5r= __t4r-__t5r;					__t5i= __t4i-__t5i;\
\
	__t4r= __t3r-__t6r-__t5r;			__t4i= __t3i-__t6i-__t5i;\
	__t5r= __t3r+__t5r;					__t5i= __t3i+__t5i;\
	__t3r= __t3r+__t6r;					__t3i= __t3i+__t6i;\
\
	__B1r =__t0r-__t3i;					__B1i =__t0i+__t3r;\
	__B2r =__t1r-__t4i;					__B2i =__t1i+__t4r;\
	__B3r =__t2r+__t5i;					__B3i =__t2i-__t5r;\
	__B4r =__t2r-__t5i;					__B4i =__t2i+__t5r;\
	__B5r =__t1r+__t4i;					__B5i =__t1i-__t4r;\
	__B6r =__t0r+__t3i;					__B6i =__t0i-__t3r;\
}

/* Totals: 60 ADD, 36 MUL	*/
#define RADIX_07_DFT_PFETCH(\
	__Ar0,__Ai0,\
	__Ar1,__Ai1,\
	__Ar2,__Ai2,\
	__Ar3,__Ai3,\
	__Ar4,__Ai4,\
	__Ar5,__Ai5,\
	__Ar6,__Ai6,\
	__tr0,__ti0,\
	__tr1,__ti1,\
	__tr2,__ti2,\
	__tr3,__ti3,\
	__tr4,__ti4,\
	__tr5,__ti5,\
	__tr6,__ti6,\
	__Br0,__Bi0,\
	__Br1,__Bi1,\
	__Br2,__Bi2,\
	__Br3,__Bi3,\
	__Br4,__Bi4,\
	__Br5,__Bi5,\
	__Br6,__Bi6,\
	__cc1,__ss1,\
	__cc2,__ss2,\
	__cc3,__ss3,\
	__rt,__it,\
	__re,__im,\
	__addr,__addp,__add_offset0,__add_offset1,__add_offset2)\
{\
	__tr0 = __Ar0;												__ti0 = __Ai0;\
	__tr6 = __Ar1 - __Ar6;										__ti6 = __Ai1 - __Ai6;	/* x1 - x6 */\
	__tr1 = __Ar1 + __Ar6;										__ti1 = __Ai1 + __Ai6;	/* x1 + x6 */\
\
	__tr5 = __Ar2 - __Ar5;										__ti5 = __Ai2 - __Ai5;	/* x2 - x5 */\
	__tr2 = __Ar2 + __Ar5;										__ti2 = __Ai2 + __Ai5;	/* x2 + x5 */\
\
	__tr4 = __Ar3 - __Ar4;										__ti4 = __Ai3 - __Ai4;	/* x3 - x4 */\
	__tr3 = __Ar3 + __Ar4;										__ti3 = __Ai3 + __Ai4;	/* x3 + x4 */\
\
	__addp = __addr + __add_offset0;	prefetch_p_doubles(__addp);\
\
	__Br0 = __tr0 + __tr1 + __tr2 + __tr3;						__Bi0 = __ti0 + __ti1 + __ti2 + __ti3;\
	__rt  = __tr0 + __cc1*__tr1 + __cc2*__tr2 + __cc3*__tr3;	__it  = __ti0 + __cc1*__ti1 + __cc2*__ti2 + __cc3*__ti3;\
	__re  = __tr0 + __cc2*__tr1 + __cc3*__tr2 + __cc1*__tr3;	__im  = __ti0 + __cc2*__ti1 + __cc3*__ti2 + __cc1*__ti3;\
	__tr0 = __tr0 + __cc3*__tr1 + __cc1*__tr2 + __cc2*__tr3;	__ti0 = __ti0 + __cc3*__ti1 + __cc1*__ti2 + __cc2*__ti3;\
\
	__addp = __addr + __add_offset1;	prefetch_p_doubles(__addp);\
\
	__tr1 = __ss1*__tr6 + __ss2*__tr5 + __ss3*__tr4;			__ti1 = __ss1*__ti6 + __ss2*__ti5 + __ss3*__ti4;\
	__tr2 = __ss2*__tr6 - __ss3*__tr5 - __ss1*__tr4;			__ti2 = __ss2*__ti6 - __ss3*__ti5 - __ss1*__ti4;\
	__tr3 = __ss3*__tr6 - __ss1*__tr5 + __ss2*__tr4;			__ti3 = __ss3*__ti6 - __ss1*__ti5 + __ss2*__ti4;\
\
	__addp = __addr + __add_offset2;	prefetch_p_doubles(__addp);\
\
	/* Output permutation causes signs to get flipped here: */\
	__Br1 = __rt  - __ti1;										__Bi1 = __it  + __tr1;\
	__Br2 = __re  - __ti2;										__Bi2 = __im  + __tr2;\
	__Br3 = __tr0 - __ti3;										__Bi3 = __ti0 + __tr3;\
	__Br4 = __tr0 + __ti3;										__Bi4 = __ti0 - __tr3;\
	__Br5 = __re  + __ti2;										__Bi5 = __im  - __tr2;\
	__Br6 = __rt  + __ti1;										__Bi6 = __it  - __tr1;\
}

/****** RADIX = 8: ******/

/* Totals: 52 ADD, 4 MUL	*/
#define RADIX_08_DIF(\
	__Ar0,__Ai0,\
	__Ar1,__Ai1,\
	__Ar2,__Ai2,\
	__Ar3,__Ai3,\
	__Ar4,__Ai4,\
	__Ar5,__Ai5,\
	__Ar6,__Ai6,\
	__Ar7,__Ai7,\
	__tr0,__ti0,\
	__tr1,__ti1,\
	__tr2,__ti2,\
	__tr3,__ti3,\
	__tr4,__ti4,\
	__tr5,__ti5,\
	__tr6,__ti6,\
	__tr7,__ti7,\
	__Br0,__Bi0,\
	__Br1,__Bi1,\
	__Br2,__Bi2,\
	__Br3,__Bi3,\
	__Br4,__Bi4,\
	__Br5,__Bi5,\
	__Br6,__Bi6,\
	__Br7,__Bi7,\
	__rt,__it)\
{\
	__tr0 = __Ar0;			__ti0 = __Ai0;\
	__rt  = __Ar4;			__it  = __Ai4;\
	__tr1 = __tr0 - __rt;		__ti1 = __ti0 - __it;\
	__tr0 = __tr0 + __rt;		__ti0 = __ti0 + __it;\
\
	__tr2 = __Ar2;			__ti2 = __Ai2;\
	__rt  = __Ar6;			__it  = __Ai6;\
	__tr3 = __tr2 - __rt;		__ti3 = __ti2 - __it;\
	__tr2 = __tr2 + __rt;		__ti2 = __ti2 + __it;\
\
	__tr4 = __Ar1;			__ti4 = __Ai1;\
	__rt  = __Ar5;			__it  = __Ai5;\
	__tr5 = __tr4 - __rt;		__ti5 = __ti4 - __it;\
	__tr4 = __tr4 + __rt;		__ti4 = __ti4 + __it;\
\
	__tr6 = __Ar3;			__ti6 = __Ai3;\
	__rt  = __Ar7;			__it  = __Ai7;\
	__tr7 = __tr6 - __rt;		__ti7 = __ti6 - __it;\
	__tr6 = __tr6 + __rt;		__ti6 = __ti6 + __it;\
/*       combine to get the 2 length-4 transforms...	*/\
	__rt  = __tr2;			__it  = __ti2;\
	__tr2 = __tr0 - __rt;		__ti2 = __ti0 - __it;\
	__tr0 = __tr0 + __rt;		__ti0 = __ti0 + __it;\
\
	__rt  = __tr3;			__it  = __ti3;\
	__tr3 = __tr1 + __it;		__ti3 = __ti1 - __rt;\
	__tr1 = __tr1 - __it;		__ti1 = __ti1 + __rt;\
\
	__rt  = __tr6;			__it  = __ti6;\
	__tr6 = __tr4 - __rt;		__ti6 = __ti4 - __it;\
	__tr4 = __tr4 + __rt;		__ti4 = __ti4 + __it;\
\
	__rt  = __tr7;			__it  = __ti7;\
	__tr7 = __tr5 + __it;		__ti7 = __ti5 - __rt;\
	__tr5 = __tr5 - __it;		__ti5 = __ti5 + __rt;\
/*       now combine the two half-transforms	*/\
	__Br0 = __tr0 + __tr4;		__Bi0 = __ti0 + __ti4;\
	__Br1 = __tr0 - __tr4;		__Bi1 = __ti0 - __ti4;\
\
	__Br2 = __tr2 - __ti6;		__Bi2 = __ti2 + __tr6;\
	__Br3 = __tr2 + __ti6;		__Bi3 = __ti2 - __tr6;\
\
	__rt = (__tr5 - __ti5)*ISRT2;	__it = (__tr5 + __ti5)*ISRT2;\
	__Br4 = __tr1 + __rt;		__Bi4 = __ti1 + __it;\
	__Br5 = __tr1 - __rt;		__Bi5 = __ti1 - __it;\
\
	__rt = (__tr7 + __ti7)*ISRT2;	__it = (__ti7 - __tr7)*ISRT2;\
	__Br6 = __tr3 - __rt;		__Bi6 = __ti3 - __it;\
	__Br7 = __tr3 + __rt;		__Bi7 = __ti3 + __it;\
}

// Radix-8 DIF subtransform macro for use by larger radix-8*k macros.
// OOP = out of place, i.e. Assumes output locs != input locs.
// For the scalar-data versino of this macro it doesn't matter whether the in/outputs are local
// scalars or array locations, but the __A-naming indicates that this macro is intended to serve
// as a prototype for a SIMD/ASM macro in which the __A-inputs are read from an array with arbitrary
// index stride and the __t-outputs go into a block of contiguous local-allocated storage:
#define RADIX_08_DIF_OOP(\
	__Ar0,__Ai0,\
	__Ar1,__Ai1,\
	__Ar2,__Ai2,\
	__Ar3,__Ai3,\
	__Ar4,__Ai4,\
	__Ar5,__Ai5,\
	__Ar6,__Ai6,\
	__Ar7,__Ai7,\
	__tr0,__ti0,\
	__tr1,__ti1,\
	__tr2,__ti2,\
	__tr3,__ti3,\
	__tr4,__ti4,\
	__tr5,__ti5,\
	__tr6,__ti6,\
	__tr7,__ti7)\
{\
	double __rt,__it;\
	__tr0 = __Ar0;			__ti0 = __Ai0;\
	__rt  = __Ar4;			__it  = __Ai4;\
	__tr1 = __tr0 - __rt;	__ti1 = __ti0 - __it;\
	__tr0 = __tr0 + __rt;	__ti0 = __ti0 + __it;\
\
	__tr2 = __Ar2;			__ti2 = __Ai2;\
	__rt  = __Ar6;			__it  = __Ai6;\
	__tr3 = __tr2 - __rt;	__ti3 = __ti2 - __it;\
	__tr2 = __tr2 + __rt;	__ti2 = __ti2 + __it;\
\
	__tr4 = __Ar1;			__ti4 = __Ai1;\
	__rt  = __Ar5;			__it  = __Ai5;\
	__tr5 = __tr4 - __rt;	__ti5 = __ti4 - __it;\
	__tr4 = __tr4 + __rt;	__ti4 = __ti4 + __it;\
\
	__tr6 = __Ar3;			__ti6 = __Ai3;\
	__rt  = __Ar7;			__it  = __Ai7;\
	__tr7 = __tr6 - __rt;	__ti7 = __ti6 - __it;\
	__tr6 = __tr6 + __rt;	__ti6 = __ti6 + __it;\
/* Combine to get the 2 length-4 transforms... */\
	__rt  = __tr2;			__it  = __ti2;\
	__tr2 = __tr0 - __rt;	__ti2 = __ti0 - __it;\
	__tr0 = __tr0 + __rt;	__ti0 = __ti0 + __it;\
\
	__rt  = __tr3;			__it  = __ti3;\
	__tr3 = __tr1 + __it;	__ti3 = __ti1 - __rt;\
	__tr1 = __tr1 - __it;	__ti1 = __ti1 + __rt;\
\
	__rt  = __tr6;			__it  = __ti6;\
	__tr6 = __tr4 - __rt;	__ti6 = __ti4 - __it;\
	__tr4 = __tr4 + __rt;	__ti4 = __ti4 + __it;\
\
	__rt  = __tr7;			__it  = __ti7;\
	__tr7 = __tr5 + __it;	__ti7 = __ti5 - __rt;\
	__tr5 = __tr5 - __it;	__ti5 = __ti5 + __rt;\
/* Now combine the two half-transforms: */\
	__rt  = __tr4;			__it  = __ti4;\
	__tr4 = __tr0 - __rt;	__ti4 = __ti0 - __it;\
	__tr0 = __tr0 + __rt;	__ti0 = __ti0 + __it;\
\
	__rt  = __tr6;			__it  = __ti6;\
	__tr6 = __tr2 + __it;	__ti6 = __ti2 - __rt;\
	__tr2 = __tr2 - __it;	__ti2 = __ti2 + __rt;\
\
	__rt = __tr5 - __ti5;	__it = __tr5 + __ti5;\
	__rt *= ISRT2;			__it *= ISRT2;\
	__tr5 = __tr1 - __rt;	__ti5 = __ti1 - __it;\
	__tr1 = __tr1 + __rt;	__ti1 = __ti1 + __it;\
\
	__rt = __tr7 + __ti7;	__it = __ti7 - __tr7;\
	__rt *= ISRT2;			__it *= ISRT2;\
	__tr7 = __tr3 + __rt;	__ti7 = __ti3 + __it;\
	__tr3 = __tr3 - __rt;	__ti3 = __ti3 - __it;\
}

// With-twiddles analog of above DIF macro: 7 nontrivial complex input twiddles E1-7 [E0 assumed = 1],
// The DIF version of this macro processes the twiddles in bit-reversed order: 0,4,2,6,1,5,3,7.
// For the scalar-data versino of this macro it doesn't matter whether the in/outputs are local
// scalars or array locations, but the __B-naming indicates that this macro is intended to serve
// as a prototype for a SIMD/ASM macro in which the __t-inputs are read from a block of contiguous
// local-allocated storage and the __B-outputs go into an array with arbitrary index stride:
#define RADIX_08_DIF_TWIDDLE_OOP(\
	__tr0,__ti0,/* Inputs */\
	__tr1,__ti1,\
	__tr2,__ti2,\
	__tr3,__ti3,\
	__tr4,__ti4,\
	__tr5,__ti5,\
	__tr6,__ti6,\
	__tr7,__ti7,\
	__Br0,__Bi0,/* Outputs */\
	__Br1,__Bi1,\
	__Br2,__Bi2,\
	__Br3,__Bi3,\
	__Br4,__Bi4,\
	__Br5,__Bi5,\
	__Br6,__Bi6,\
	__Br7,__Bi7,\
	__Wr4,__Wi4,/* Roots of unity - 'W' is for 'Wurzel' */\
	__Wr2,__Wi2,\
	__Wr6,__Wi6,\
	__Wr1,__Wi1,\
	__Wr5,__Wi5,\
	__Wr3,__Wi3,\
	__Wr7,__Wi7)\
{\
	double __rt,__it;\
	__rt  = __tr1*__Wr4 - __ti1*__Wi4;	__it  = __ti1*__Wr4 + __tr1*__Wi4;\
	__tr1 = __tr0 - __rt;	__ti1 = __ti0 - __it;\
	__tr0 = __tr0 + __rt;	__ti0 = __ti0 + __it;\
\
	__rt  = __tr2*__Wr2 - __ti2*__Wi2;	__ti2 = __ti2*__Wr2 + __tr2*__Wi2;	__tr2 = __rt;\
	__rt  = __tr3*__Wr6 - __ti3*__Wi6;	__it  = __ti3*__Wr6 + __tr3*__Wi6;\
	__tr3 = __tr2 - __rt;	__ti3 = __ti2 - __it;\
	__tr2 = __tr2 + __rt;	__ti2 = __ti2 + __it;\
\
	__rt  = __tr4*__Wr1 - __ti4*__Wi1;	__ti4 = __ti4*__Wr1 + __tr4*__Wi1;	__tr4 = __rt;\
	__rt  = __tr5*__Wr5 - __ti5*__Wi5;	__it  = __ti5*__Wr5 + __tr5*__Wi5;\
	__tr5 = __tr4 - __rt;	__ti5 = __ti4 - __it;\
	__tr4 = __tr4 + __rt;	__ti4 = __ti4 + __it;\
\
	__rt  = __tr6*__Wr3 - __ti6*__Wi3;	__ti6 = __ti6*__Wr3 + __tr6*__Wi3;	__tr6 = __rt;\
	__rt  = __tr7*__Wr7 - __ti7*__Wi7;	__it  = __ti7*__Wr7 + __tr7*__Wi7;\
	__tr7 = __tr6 - __rt;	__ti7 = __ti6 - __it;\
	__tr6 = __tr6 + __rt;	__ti6 = __ti6 + __it;\
/* Combine to get the 2 length-4 transforms... */\
	__rt  = __tr2;			__it  = __ti2;\
	__tr2 = __tr0 - __rt;	__ti2 = __ti0 - __it;\
	__tr0 = __tr0 + __rt;	__ti0 = __ti0 + __it;\
\
	__rt  = __tr3;			__it  = __ti3;\
	__tr3 = __tr1 + __it;	__ti3 = __ti1 - __rt;\
	__tr1 = __tr1 - __it;	__ti1 = __ti1 + __rt;\
\
	__rt  = __tr6;			__it  = __ti6;\
	__tr6 = __tr4 - __rt;	__ti6 = __ti4 - __it;\
	__tr4 = __tr4 + __rt;	__ti4 = __ti4 + __it;\
\
	__rt  = __tr7;			__it  = __ti7;\
	__tr7 = __tr5 + __it;	__ti7 = __ti5 - __rt;\
	__tr5 = __tr5 - __it;	__ti5 = __ti5 + __rt;\
/* Now combine the two half-transforms: */\
	__rt  = __tr4;			__it  = __ti4;\
	__Br0 = __tr0 + __rt;	__Bi0 = __ti0 + __it;\
	__Br1 = __tr0 - __rt;	__Bi1 = __ti0 - __it;\
\
	__rt  = __tr6;			__it  = __ti6;\
	__Br2 = __tr2 - __it;	__Bi2 = __ti2 + __rt;\
	__Br3 = __tr2 + __it;	__Bi3 = __ti2 - __rt;\
\
	__rt  = __tr5 - __ti5;	__it  = __tr5 + __ti5;\
	__rt *= ISRT2;			__it *= ISRT2;\
	__Br4 = __tr1 + __rt;	__Bi4 = __ti1 + __it;\
	__Br5 = __tr1 - __rt;	__Bi5 = __ti1 - __it;\
\
	__rt  = __tr7 + __ti7;	__it  = __ti7 - __tr7;\
	__rt *= ISRT2;			__it *= ISRT2;\
	__Br6 = __tr3 - __rt;	__Bi6 = __ti3 - __it;\
	__Br7 = __tr3 + __rt;	__Bi7 = __ti3 + __it;\
}


// Radix-8 DIT subtransform macro for use by larger radix-8*k macros.
// OOP = out of place, i.e. Assumes output locs != input locs.
// For the scalar-data versino of this macro it doesn't matter whether the in/outputs are local
// scalars or array locations, but the __A-naming indicates that this macro is intended to serve
// as a prototype for a SIMD/ASM macro in which the __A-inputs are read from an array with arbitrary
// index stride and the __t-outputs go into a block of contiguous local-allocated storage:
#define RADIX_08_DIT_OOP(\
	__Ar0,__Ai0,\
	__Ar1,__Ai1,\
	__Ar2,__Ai2,\
	__Ar3,__Ai3,\
	__Ar4,__Ai4,\
	__Ar5,__Ai5,\
	__Ar6,__Ai6,\
	__Ar7,__Ai7,\
	__tr0,__ti0,\
	__tr1,__ti1,\
	__tr2,__ti2,\
	__tr3,__ti3,\
	__tr4,__ti4,\
	__tr5,__ti5,\
	__tr6,__ti6,\
	__tr7,__ti7)\
{\
	double __rt,__it;\
	__tr0 = __Ar0;			__ti0 = __Ai0;\
	__rt  = __Ar1;			__it  = __Ai1;\
	__tr1 = __tr0 - __rt;	__ti1 = __ti0 - __it;\
	__tr0 = __tr0 + __rt;	__ti0 = __ti0 + __it;\
\
	__tr2 = __Ar2;			__ti2 = __Ai2;\
	__rt  = __Ar3;			__it  = __Ai3;\
	__tr3 = __tr2 - __rt;	__ti3 = __ti2 - __it;\
	__tr2 = __tr2 + __rt;	__ti2 = __ti2 + __it;\
\
	__tr4 = __Ar4;			__ti4 = __Ai4;\
	__rt  = __Ar5;			__it  = __Ai5;\
	__tr5 = __tr4 - __rt;	__ti5 = __ti4 - __it;\
	__tr4 = __tr4 + __rt;	__ti4 = __ti4 + __it;\
\
	__tr6 = __Ar6;			__ti6 = __Ai6;\
	__rt  = __Ar7;			__it  = __Ai7;\
	__tr7 = __tr6 - __rt;	__ti7 = __ti6 - __it;\
	__tr6 = __tr6 + __rt;	__ti6 = __ti6 + __it;\
/* Combine to get the 2 length-4 transforms... */\
	__rt  = __tr2;			__it  = __ti2;\
	__tr2 = __tr0 - __rt;	__ti2 = __ti0 - __it;\
	__tr0 = __tr0 + __rt;	__ti0 = __ti0 + __it;\
\
	__rt  = __tr3;			__it  = __ti3;\
	__tr3 = __tr1 - __it;	__ti3 = __ti1 + __rt;\
	__tr1 = __tr1 + __it;	__ti1 = __ti1 - __rt;\
\
	__rt  = __tr6;			__it  = __ti6;\
	__tr6 = __tr4 - __rt;	__ti6 = __ti4 - __it;\
	__tr4 = __tr4 + __rt;	__ti4 = __ti4 + __it;\
\
	__rt  = __tr7;			__it  = __ti7;\
	__tr7 = __tr5 - __it;	__ti7 = __ti5 + __rt;\
	__tr5 = __tr5 + __it;	__ti5 = __ti5 - __rt;\
/* Now combine the two half-transforms: */\
	__rt  = __tr4;			__it  = __ti4;\
	__tr4 = __tr0 - __rt;	__ti4 = __ti0 - __it;\
	__tr0 = __tr0 + __rt;	__ti0 = __ti0 + __it;\
\
	__rt  = __tr6;			__it  = __ti6;\
	__tr6 = __tr2 - __it;	__ti6 = __ti2 + __rt;\
	__tr2 = __tr2 + __it;	__ti2 = __ti2 - __rt;\
\
	__rt  = __tr5 + __ti5;	__it  = __tr5 - __ti5;\
	__rt *= ISRT2;			__it *= ISRT2;\
	__tr5 = __tr1 - __rt;	__ti5 = __ti1 + __it;\
	__tr1 = __tr1 + __rt;	__ti1 = __ti1 - __it;\
\
	__rt  = __tr7 - __ti7;	__it  = __ti7 + __tr7;\
	__rt *= ISRT2;			__it *= ISRT2;\
	__tr7 = __tr3 + __rt;	__ti7 = __ti3 + __it;\
	__tr3 = __tr3 - __rt;	__ti3 = __ti3 - __it;\
}

// With-twiddles analog of above DIT macro: 7 nontrivial complex input twiddles E1-7 [E0 assumed = 1],
// The DIT version of this macro processes the twiddles in bit-reversed order: 0,4,2,6,1,5,3,7.
// For the scalar-data versino of this macro it doesn't matter whether the in/outputs are local
// scalars or array locations, but the __B-naming indicates that this macro is intended to serve
// as a prototype for a SIMD/ASM macro in which the __t-inputs are read from a block of contiguous
// local-allocated storage and the __B-outputs go into an array with arbitrary index stride:
#define RADIX_08_DIT_TWIDDLE_OOP(\
	__tr0,__ti0,/* Inputs */\
	__tr1,__ti1,\
	__tr2,__ti2,\
	__tr3,__ti3,\
	__tr4,__ti4,\
	__tr5,__ti5,\
	__tr6,__ti6,\
	__tr7,__ti7,\
	__Br0,__Bi0,/* Outputs */\
	__Br1,__Bi1,\
	__Br2,__Bi2,\
	__Br3,__Bi3,\
	__Br4,__Bi4,\
	__Br5,__Bi5,\
	__Br6,__Bi6,\
	__Br7,__Bi7,\
	__Wr1,__Wi1,/* Roots of unity - 'W' is for 'Wurzel' */\
	__Wr2,__Wi2,\
	__Wr3,__Wi3,\
	__Wr4,__Wi4,\
	__Wr5,__Wi5,\
	__Wr6,__Wi6,\
	__Wr7,__Wi7)\
{\
/* Original-version code:\
	double __rt,__it;\
	__rt  = __tr1*__Wr1 + __ti1*__Wi1;	__it  = __ti1*__Wr1 - __tr1*__Wi1;\
	__tr1 = __tr0 - __rt;	__ti1 = __ti0 - __it;\
	__tr0 = __tr0 + __rt;	__ti0 = __ti0 + __it;\
\
	__rt  = __tr2*__Wr2 + __ti2*__Wi2;	__ti2 = __ti2*__Wr2 - __tr2*__Wi2;	__tr2 = __rt;\
	__rt  = __tr3*__Wr3 + __ti3*__Wi3;	__it  = __ti3*__Wr3 - __tr3*__Wi3;\
	__tr3 = __tr2 - __rt;	__ti3 = __ti2 - __it;\
	__tr2 = __tr2 + __rt;	__ti2 = __ti2 + __it;\
\
	__rt  = __tr4*__Wr4 + __ti4*__Wi4;	__ti4 = __ti4*__Wr4 - __tr4*__Wi4;	__tr4 = __rt;\
	__rt  = __tr5*__Wr5 + __ti5*__Wi5;	__it  = __ti5*__Wr5 - __tr5*__Wi5;\
	__tr5 = __tr4 - __rt;	__ti5 = __ti4 - __it;\
	__tr4 = __tr4 + __rt;	__ti4 = __ti4 + __it;\
\
	__rt  = __tr6*__Wr6 + __ti6*__Wi6;	__ti6 = __ti6*__Wr6 - __tr6*__Wi6;	__tr6 = __rt;\
	__rt  = __tr7*__Wr7 + __ti7*__Wi7;	__it  = __ti7*__Wr7 - __tr7*__Wi7;\
	__tr7 = __tr6 - __rt;	__ti7 = __ti6 - __it;\
	__tr6 = __tr6 + __rt;	__ti6 = __ti6 + __it;\
** Combine to get the 2 length-4 transforms... **\
	__rt  = __tr2;			__it  = __ti2;\
	__tr2 = __tr0 - __rt;	__ti2 = __ti0 - __it;\
	__tr0 = __tr0 + __rt;	__ti0 = __ti0 + __it;\
\
	__rt  = __tr3;			__it  = __ti3;\
	__tr3 = __tr1 - __it;	__ti3 = __ti1 + __rt;\
	__tr1 = __tr1 + __it;	__ti1 = __ti1 - __rt;\
\
	__rt  = __tr6;			__it  = __ti6;\
	__tr6 = __tr4 - __rt;	__ti6 = __ti4 - __it;\
	__tr4 = __tr4 + __rt;	__ti4 = __ti4 + __it;\
\
	__rt  = __tr7;			__it  = __ti7;\
	__tr7 = __tr5 - __it;	__ti7 = __ti5 + __rt;\
	__tr5 = __tr5 + __it;	__ti5 = __ti5 - __rt;\
** Now combine the two half-transforms: **\
	__rt  = __tr4;			__it  = __ti4;\
	__Br0 = __tr0 + __rt;	__Bi0 = __ti0 + __it;\
	__Br1 = __tr0 - __rt;	__Bi1 = __ti0 - __it;\
\
	__rt  = __tr6;			__it  = __ti6;\
	__Br2 = __tr2 + __it;	__Bi2 = __ti2 - __rt;\
	__Br3 = __tr2 - __it;	__Bi3 = __ti2 + __rt;\
\
	__rt  = __tr5 + __ti5;	__it  = __tr5 - __ti5;\
	__rt *= ISRT2;			__it *= ISRT2;\
	__Br4 = __tr1 + __rt;	__Bi4 = __ti1 - __it;\
	__Br5 = __tr1 - __rt;	__Bi5 = __ti1 + __it;\
\
	__rt  = __tr7 - __ti7;	__it  = __ti7 + __tr7;\
	__rt *= ISRT2;			__it *= ISRT2;\
	__Br6 = __tr3 - __rt;	__Bi6 = __ti3 - __it;\
	__Br7 = __tr3 + __rt;	__Bi7 = __ti3 + __it;\
*/\
/* Rewrite with x86 ISA and SSE2/16-register target: Assume 16 registers _r0-f available, all arithmetic is on these: */\
/* In 16-register mode Block pairs 0/1,2/3 and 4/5,6/7 can be done at same time; 8-register mode just replace _r8-f with _r0-7. */\
	double _r0,_r1,_r2,_r3,_r4,_r5,_r6,_r7,_r8,_r9,_ra,_rb,_rc,_rd,_re,_rf;\
	const double _two = 2;\
	/* Block 0/1 has just one twiddle-CMUL: */\
	_r0  = __tr0;	_r1  = __ti0;\
	/* [r4,r5] = CMUL(__t1,__W1): */\
	_r4  = __tr1;	_r5  = __ti1;\
	_r6  = _r5;		_r7  = _r4;\
	_r4 *= __Wr1;	_r5 *= __Wr1;\
	_r6 *= __Wi1;	_r7 *= __Wi1;\
	_r4 += _r6;		_r5 -= _r7;\
	_r2  = _r0;		_r3  = _r1;\
	_r0 += _r4;		_r1 += _r5;\
	_r2 -= _r4;		_r3 -= _r5;\
	__tr0 = _r0;	__ti0 = _r1;\
	__tr1 = _r2;	__ti1 = _r3;\
\
	/* Blocks 2/3: */\
	/* [r8,r9] = CMUL(__t2,__W2): */\
	_r8  = __tr2;	_r9  = __ti2;\
	_ra  = _r9;		_rb  = _r8;\
	_r8 *= __Wr2;	_r9 *= __Wr2;\
	_ra *= __Wi2;	_rb *= __Wi2;\
	_r8 += _ra;		_r9 -= _rb;\
	/* [rc,rd] = CMUL(__t3,__W3): */\
	_rc  = __tr3;	_rd  = __ti3;\
	_re  = _rd;		_rf  = _rc;\
	_rc *= __Wr3;	_rd *= __Wr3;\
	_re *= __Wi3;	_rf *= __Wi3;\
	_rc += _re;		_rd -= _rf;\
	/* Now do radix-2 butterfly: */\
	_ra  = _r8;		_rb  = _r9;\
	_r8 += _rc;		_r9 += _rd;\
	_ra -= _rc;		_rb -= _rd;\
	__tr2 = _r8;	__ti2 = _r9;\
	__tr3 = _ra;	__ti3 = _rb;\
\
	/* Blocks 4/5: */\
	/* [r0,r1] = CMUL(__t4,__W4): */\
	_r0  = __tr4;	_r1  = __ti4;\
	_r2  = _r1;		_r3  = _r0;\
	_r0 *= __Wr4;	_r1 *= __Wr4;\
	_r2 *= __Wi4;	_r3 *= __Wi4;\
	_r0 += _r2;		_r1 -= _r3;\
	/* [r4,r5] = CMUL(__t5,__W5): */\
	_r4  = __tr5;	_r5  = __ti5;\
	_r6  = _r5;		_r7  = _r4;\
	_r4 *= __Wr5;	_r5 *= __Wr5;\
	_r6 *= __Wi5;	_r7 *= __Wi5;\
	_r4 += _r6;		_r5 -= _r7;\
	/* Now do radix-2 butterfly: */\
	_r2  = _r0;		_r3  = _r1;\
	_r0 += _r4;		_r1 += _r5;\
	_r2 -= _r4;		_r3 -= _r5;\
\
	/* Blocks 6/7: */\
	/* [r8,r9] = CMUL(__t6,__W6): */\
	_r8  = __tr6;	_r9  = __ti6;\
	_ra  = _r9;		_rb  = _r8;\
	_r8 *= __Wr6;	_r9 *= __Wr6;\
	_ra *= __Wi6;	_rb *= __Wi6;\
	_r8 += _ra;		_r9 -= _rb;\
	/* [rc,rd] = CMUL(__t7,__W7): */\
	_rc  = __tr7;	_rd  = __ti7;\
	_re  = _rd;		_rf  = _rc;\
	_rc *= __Wr7;	_rd *= __Wr7;\
	_re *= __Wi7;	_rf *= __Wi7;\
	_rc += _re;		_rd -= _rf;\
	/* Now do radix-2 butterfly: */\
	_ra  = _r8;		_rb  = _r9;\
	_r8 += _rc;		_r9 += _rd;\
	_ra -= _rc;		_rb -= _rd;\
\
/* Reload Block 0-3 outputs into r4-7,c-f, combine to get the 2 length-4 subtransform... */\
	_r4 = __tr0;	_r5 = __ti0;\
	_r6 = __tr1;	_r7 = __ti1;\
	_rc = __tr2;	_rd = __ti2;\
	_re = __tr3;	_rf = __ti3;\
\
	_r4 -= _rc;		_r5 -= _rd;\
	_r6 -= _rf;		_r7 -= _re;\
	_r0 -= _r8;		_r1 -= _r9;\
	_r2 -= _rb;		_r3 -= _ra;\
	_rc *= _two;	_rd *= _two;\
	_rf *= _two;	_re *= _two;\
	_r8 *= _two;	_r9 *= _two;\
	_rb *= _two;	_ra *= _two;\
	_rc += _r4;		_rd += _r5;\
	_rf += _r6;		_re += _r7;\
	_r8 += _r0;		_r9 += _r1;\
	_rb += _r2;		_ra += _r3;\
	/* In terms of our original scalar-code prototyping macro, here are where the data are:
	__tr0 = _rc;	__ti0 = _rd;\
	__tr1 = _rf;	__ti1 = _r7;\
	__tr2 = _r4;	__ti2 = _r5;\
	__tr3 = _r6;	__ti3 = _re;\
	__tr4 = _r8;	__ti4 = _r9;\
	__tr5 = _rb;	__ti5 = _r3;\
	__tr6 = _r0;	__ti6 = _r1;\
	__tr7 = _r2;	__ti7 = _ra;\
	*/\
\
/* Now combine the two half-transforms: */\
	/* Need r2/3 +- a/b combos for the *ISRT2 preceding the output 4-7 radix-2 butterflies, so start them first: */\
	_rb -= _r3;		_r2 -= _ra;\
	_rc -= _r8;		_rd -= _r9;\
	_r4 -= _r1;		_r5 -= _r0;	\
	_r3 *= _two;	_ra *= _two;\
	_r8 *= _two;	_r9 *= _two;\
	_r1 *= _two;	_r0 *= _two;\
	_r3 += _rb;		_ra += _r2;\
	_r8 += _rc;		_r9 += _rd;\
	_r1 += _r4;		_r0 += _r5;	\
\
	__Br0 = _r8;	__Bi0 = _r9;\
	__Br1 = _rc;	__Bi1 = _rd;\
	__Br2 = _r1;	__Bi2 = _r5;\
	__Br3 = _r4;	__Bi3 = _r0;\
\
	_r3 *= ISRT2;	_rb *= ISRT2;\
	_r2 *= ISRT2;	_ra *= ISRT2;\
\
	_rf -= _r3;		_r7 -= _rb;\
	_r6 -= _r2;		_re -= _ra;\
	_r3 *= _two;	_rb *= _two;\
	_r2 *= _two;	_ra *= _two;\
	_r3 += _rf;		_rb += _r7;\
	_r2 += _r6;		_ra += _re;\
\
	__Br4 = _r3;	__Bi4 = _r7;\
	__Br5 = _rf;	__Bi5 = _rb;\
	__Br6 = _r6;	__Bi6 = _re;\
	__Br7 = _r2;	__Bi7 = _ra;\
}

/* Totals: 52 ADD, 4 MUL	*/
#define RADIX_08_DIF_PFETCH(\
	__Ar0,__Ai0,\
	__Ar1,__Ai1,\
	__Ar2,__Ai2,\
	__Ar3,__Ai3,\
	__Ar4,__Ai4,\
	__Ar5,__Ai5,\
	__Ar6,__Ai6,\
	__Ar7,__Ai7,\
	__tr0,__ti0,\
	__tr1,__ti1,\
	__tr2,__ti2,\
	__tr3,__ti3,\
	__tr4,__ti4,\
	__tr5,__ti5,\
	__tr6,__ti6,\
	__tr7,__ti7,\
	__Br0,__Bi0,\
	__Br1,__Bi1,\
	__Br2,__Bi2,\
	__Br3,__Bi3,\
	__Br4,__Bi4,\
	__Br5,__Bi5,\
	__Br6,__Bi6,\
	__Br7,__Bi7,\
	__rt,__it,\
	__add_offset0,__add_offset1,__add_offset2,__add_offset3,__add_offset4)\
{\
	__tr0 = __Ar0;			__ti0 = __Ai0;\
	__rt  = __Ar4;			__it  = __Ai4;\
	__tr1 = __tr0 - __rt;		__ti1 = __ti0 - __it;\
	__tr0 = __tr0 + __rt;		__ti0 = __ti0 + __it;\
\
	__tr2 = __Ar2;			__ti2 = __Ai2;\
	__rt  = __Ar6;			__it  = __Ai6;\
	__tr3 = __tr2 - __rt;		__ti3 = __ti2 - __it;\
	__tr2 = __tr2 + __rt;		__ti2 = __ti2 + __it;\
\
addr = add0 + __add_offset0;	prefetch_p_doubles(addr);\
\
	__tr4 = __Ar1;			__ti4 = __Ai1;\
	__rt  = __Ar5;			__it  = __Ai5;\
	__tr5 = __tr4 - __rt;		__ti5 = __ti4 - __it;\
	__tr4 = __tr4 + __rt;		__ti4 = __ti4 + __it;\
\
	__tr6 = __Ar3;			__ti6 = __Ai3;\
	__rt  = __Ar7;			__it  = __Ai7;\
	__tr7 = __tr6 - __rt;		__ti7 = __ti6 - __it;\
	__tr6 = __tr6 + __rt;		__ti6 = __ti6 + __it;\
\
addr = add0 + __add_offset1;	prefetch_p_doubles(addr);\
\
/*       combine to get the 2 length-4 transforms...	*/\
	__rt  = __tr2;			__it  = __ti2;\
	__tr2 = __tr0 - __rt;		__ti2 = __ti0 - __it;\
	__tr0 = __tr0 + __rt;		__ti0 = __ti0 + __it;\
\
	__rt  = __tr3;			__it  = __ti3;\
	__tr3 = __tr1 + __it;		__ti3 = __ti1 - __rt;\
	__tr1 = __tr1 - __it;		__ti1 = __ti1 + __rt;\
\
addr = add0 + __add_offset2;	prefetch_p_doubles(addr);\
\
	__rt  = __tr6;			__it  = __ti6;\
	__tr6 = __tr4 - __rt;		__ti6 = __ti4 - __it;\
	__tr4 = __tr4 + __rt;		__ti4 = __ti4 + __it;\
\
	__rt  = __tr7;			__it  = __ti7;\
	__tr7 = __tr5 + __it;		__ti7 = __ti5 - __rt;\
	__tr5 = __tr5 - __it;		__ti5 = __ti5 + __rt;\
\
addr = add0 + __add_offset3;	prefetch_p_doubles(addr);\
\
/*       now combine the two half-transforms	*/\
	__Br0 = __tr0 + __tr4;		__Bi0 = __ti0 + __ti4;\
	__Br1 = __tr0 - __tr4;		__Bi1 = __ti0 - __ti4;\
\
	__Br2 = __tr2 - __ti6;		__Bi2 = __ti2 + __tr6;\
	__Br3 = __tr2 + __ti6;		__Bi3 = __ti2 - __tr6;\
\
addr = add0 + __add_offset4;	prefetch_p_doubles(addr);\
\
	__rt = (__tr5 - __ti5)*ISRT2;	__it = (__tr5 + __ti5)*ISRT2;\
	__Br4 = __tr1 + __rt;		__Bi4 = __ti1 + __it;\
	__Br5 = __tr1 - __rt;		__Bi5 = __ti1 - __it;\
\
	__rt = (__tr7 + __ti7)*ISRT2;	__it = (__ti7 - __tr7)*ISRT2;\
	__Br6 = __tr3 - __rt;		__Bi6 = __ti3 - __it;\
	__Br7 = __tr3 + __rt;		__Bi7 = __ti3 + __it;\
}

/* Totals: 52 ADD, 4 MUL	*/
#define RADIX_08_DIT(\
	__Ar0,__Ai0,\
	__Ar1,__Ai1,\
	__Ar2,__Ai2,\
	__Ar3,__Ai3,\
	__Ar4,__Ai4,\
	__Ar5,__Ai5,\
	__Ar6,__Ai6,\
	__Ar7,__Ai7,\
	__tr0,__ti0,\
	__tr1,__ti1,\
	__tr2,__ti2,\
	__tr3,__ti3,\
	__tr4,__ti4,\
	__tr5,__ti5,\
	__tr6,__ti6,\
	__tr7,__ti7,\
	__Br0,__Bi0,\
	__Br1,__Bi1,\
	__Br2,__Bi2,\
	__Br3,__Bi3,\
	__Br4,__Bi4,\
	__Br5,__Bi5,\
	__Br6,__Bi6,\
	__Br7,__Bi7,\
	__rt,__it)\
{\
	__tr0 = __Ar0;			__ti0 = __Ai0;\
	__rt  = __Ar1;			__it  = __Ai1;\
	__tr1 = __tr0 - __rt;		__ti1 = __ti0 - __it;\
	__tr0 = __tr0 + __rt;		__ti0 = __ti0 + __it;\
\
	__tr2 = __Ar2;			__ti2 = __Ai2;\
	__rt  = __Ar3;			__it  = __Ai3;\
	__tr3 = __tr2 - __rt;		__ti3 = __ti2 - __it;\
	__tr2 = __tr2 + __rt;		__ti2 = __ti2 + __it;\
\
	__tr4 = __Ar4;			__ti4 = __Ai4;\
	__rt  = __Ar5;			__it  = __Ai5;\
	__tr5 = __tr4 - __rt;		__ti5 = __ti4 - __it;\
	__tr4 = __tr4 + __rt;		__ti4 = __ti4 + __it;\
\
	__tr6 = __Ar6;			__ti6 = __Ai6;\
	__rt  = __Ar7;			__it  = __Ai7;\
	__tr7 = __tr6 - __rt;		__ti7 = __ti6 - __it;\
	__tr6 = __tr6 + __rt;		__ti6 = __ti6 + __it;\
/*       combine to get the 2 length-4 transforms...	*/\
	__rt  = __tr2;			__it  = __ti2;\
	__tr2 = __tr0 - __rt;		__ti2 = __ti0 - __it;\
	__tr0 = __tr0 + __rt;		__ti0 = __ti0 + __it;\
\
	__rt  = __tr3;			__it  = __ti3;\
	__tr3 = __tr1 + __it;		__ti3 = __ti1 - __rt;\
	__tr1 = __tr1 - __it;		__ti1 = __ti1 + __rt;\
\
	__rt  = __tr6;			__it  = __ti6;\
	__tr6 = __tr4 - __rt;		__ti6 = __ti4 - __it;\
	__tr4 = __tr4 + __rt;		__ti4 = __ti4 + __it;\
\
	__rt  = __tr7;			__it  = __ti7;\
	__tr7 = __tr5 + __it;		__ti7 = __ti5 - __rt;\
	__tr5 = __tr5 - __it;		__ti5 = __ti5 + __rt;\
/*       now combine the two half-transforms	*/\
	__Br0 = __tr0 + __tr4;		__Bi0 = __ti0 + __ti4;\
	__Br4 = __tr0 - __tr4;		__Bi4 = __ti0 - __ti4;\
\
	__Br2 = __tr2 - __ti6;		__Bi2 = __ti2 + __tr6;\
	__Br6 = __tr2 + __ti6;		__Bi6 = __ti2 - __tr6;\
\
	__rt = (__tr5 - __ti5)*ISRT2;	__it = (__tr5 + __ti5)*ISRT2;\
	__Br1 = __tr1 + __rt;		__Bi1 = __ti1 + __it;\
	__Br5 = __tr1 - __rt;		__Bi5 = __ti1 - __it;\
\
	__rt = (__tr7 + __ti7)*ISRT2;	__it = (__ti7 - __tr7)*ISRT2;\
	__Br3 = __tr3 - __rt;		__Bi3 = __ti3 - __it;\
	__Br7 = __tr3 + __rt;		__Bi7 = __ti3 + __it;\
}


/****** RADIX = 9: ******/

/*...Radix-9 DIF: A/Bs are in/outputs and t's are temporaries (all doubles): */
/* Totals: 80 ADD, 40 MUL	*/
#define RADIX_09_DIF(__A0r,__A0i,__A1r,__A1i,__A2r,__A2i,__A3r,__A3i,__A4r,__A4i,__A5r,__A5i,__A6r,__A6i,__A7r,__A7i,__A8r,__A8i\
                    ,__t00,__t01,__t02,__t03,__t04,__t05,__t06,__t07,__t08,__t09,__t0a,__t0b,__t0c,__t0d,__t0e,__t0f,__t0g,__t0h\
                    ,__rt,__it,__tt)\
{\
/*...gather the needed data (9 64-bit complex, i.e. 18 64-bit reals) and do three radix-3 transforms...	*/\
	__t00 =__A0r;				__t01 =__A0i;\
	__t02 =__A3r+__A6r;			__t03 =__A3i+__A6i;\
	__t04 =__A3r-__A6r;			__t05 =__A3i-__A6i;\
	__t00 =__t00+__t02;			__t01 =__t01+__t03;\
	__t02 =__t00+c3m1*__t02;	__t03 =__t01+c3m1*__t03;\
	__rt  =s3*__t04;			__it  =s3*__t05;\
	__t04 =__t02+__it;			__t05 =__t03-__rt;\
	__t02 =__t02-__it;			__t03 =__t03+__rt;\
\
	__t06 =__A1r;				__t07 =__A1i;\
	__t08 =__A4r+__A7r;			__t09 =__A4i+__A7i;\
	__t0a =__A4r-__A7r;			__t0b =__A4i-__A7i;\
	__t06 =__t06+__t08;			__t07 =__t07+__t09;\
	__t08 =__t06+c3m1*__t08;	__t09 =__t07+c3m1*__t09;\
	__rt  =s3*__t0a;			__it  =s3*__t0b;\
	__t0a =__t08+__it;			__t0b =__t09-__rt;\
	__t08 =__t08-__it;			__t09 =__t09+__rt;\
\
	__t0c =__A2r;				__t0d =__A2i;\
	__t0e =__A5r+__A8r;			__t0f =__A5i+__A8i;\
	__t0g =__A5r-__A8r;			__t0h =__A5i-__A8i;\
	__t0c =__t0c+__t0e;			__t0d =__t0d+__t0f;\
	__t0e =__t0c+c3m1*__t0e;	__t0f =__t0d+c3m1*__t0f;\
	__rt  =s3*__t0g;			__it  =s3*__t0h;\
	__t0g =__t0e+__it;			__t0h =__t0f-__rt;\
	__t0e =__t0e-__it;			__t0f =__t0f+__rt;\
/*...and now do three more radix-3 transforms, including the twiddle factors: */\
	__rt  =__t06;				__it  =__t07;\
	__t06 =__rt+__t0c;			__t07 =__it+__t0d;\
	__t0c =__rt-__t0c;			__t0d =__it-__t0d;\
	__t00 =__t00+__t06;			__t01 =__t01+__t07;\
	__t06 =__t00+c3m1*__t06;	__t07 =__t01+c3m1*__t07;\
	__rt  =s3*__t0c;			__it  =s3*__t0d;\
	__t0c =__t06+__it;			__t0d =__t07-__rt;\
	__t06 =__t06-__it;			__t07 =__t07+__rt;\
\
	__rt  =__t08*c -__t09*s;	__it  =__t08*s +__t09*c;\
	__tt  =__t0e*c2-__t0f*s2;	__t0f =__t0e*s2+__t0f*c2;	__t0e=__tt;\
	__t08 =__rt+__t0e;			__t09 =__it+__t0f;\
	__t0e =__rt-__t0e;			__t0f =__it-__t0f;\
	__t02 =__t02+__t08;			__t03 =__t03+__t09;\
	__t08 =__t02+c3m1*__t08;	__t09 =__t03+c3m1*__t09;\
	__rt  =s3*__t0e;			__it  =s3*__t0f;\
	__t0e =__t08+__it;			__t0f =__t09-__rt;\
	__t08 =__t08-__it;			__t09 =__t09+__rt;\
\
	__rt  =__t0a*c2-__t0b*s2;	__it  =__t0a*s2+__t0b*c2;\
	__tt  =__t0g*c4-__t0h*s4;	__t0h =__t0g*s4+__t0h*c4;	__t0g=__tt;\
	__t0a =__rt+__t0g;			__t0b =__it+__t0h;\
	__t0g =__rt-__t0g;			__t0h =__it-__t0h;\
	__t04 =__t04+__t0a;			__t05 =__t05+__t0b;\
	__t0a =__t04+c3m1*__t0a;	__t0b =__t05+c3m1*__t0b;\
	__rt  =s3*__t0g;			__it  =s3*__t0h;\
	__t0g =__t0a+__it;			__t0h =__t0b-__rt;\
	__t0a =__t0a-__it;			__t0b =__t0b+__rt;\
}


/*...Radix-9 DIT: t's are temporaries and As are outputs (all doubles) */
/* Totals: 68 ADD, 40 MUL	*/
#define RADIX_09_DIT(__t00,__t01,__t02,__t03,__t04,__t05,__t06,__t07,__t08,__t09,__t0a,__t0b,__t0c,__t0d,__t0e,__t0f,__t0g,__t0h\
					,__A0r,__A0i,__A1r,__A1i,__A2r,__A2i,__A3r,__A3i,__A4r,__A4i,__A5r,__A5i,__A6r,__A6i,__A7r,__A7i,__A8r,__A8i\
                    ,__rt,__it,__tt)\
{\
/*...gather the needed data (9 64-b__it complex, i.e. 18 64-b__it reals) and do three radix-3 transforms...	*/\
	__rt  =__t02;				__it  =__t03;\
	__t02 =__rt+__t04;			__t03 =__it+__t05;\
	__t04 =__rt-__t04;			__t05 =__it-__t05;\
	__t00 =__t00+__t02;			__t01 =__t01+__t03;\
	__t02 =__t00+c3m1*__t02;	__t03 =__t01+c3m1*__t03;\
	__rt  =s3*__t04;			__it  =s3*__t05;\
	__t04 =__t02-__it;			__t05 =__t03+__rt;\
	__t02 =__t02+__it;			__t03 =__t03-__rt;\
\
	__rt  =__t08;				__it  =__t09;\
	__t08 =__rt+__t0a;			__t09 =__it+__t0b;\
	__t0a =__rt-__t0a;			__t0b =__it-__t0b;\
	__t06 =__t06+__t08;			__t07 =__t07+__t09;\
	__t08 =__t06+c3m1*__t08;	__t09 =__t07+c3m1*__t09;\
	__rt  =s3*__t0a;			__it  =s3*__t0b;\
	__t0a =__t08-__it;			__t0b =__t09+__rt;\
	__t08 =__t08+__it;			__t09 =__t09-__rt;\
\
	__rt  =__t0e;				__it  =__t0f;\
	__t0e =__rt+__t0g;			__t0f =__it+__t0h;\
	__t0g =__rt-__t0g;			__t0h =__it-__t0h;\
	__t0c =__t0c+__t0e;			__t0d =__t0d+__t0f;\
	__t0e =__t0c+c3m1*__t0e;	__t0f =__t0d+c3m1*__t0f;\
	__rt  =s3*__t0g;			__it  =s3*__t0h;\
	__t0g =__t0e-__it;			__t0h =__t0f+__rt;\
	__t0e =__t0e+__it;			__t0f =__t0f-__rt;\
/*...and now do three more radix-3 transforms, including the twiddle factors: */\
	__rt  =__t06;				__it  =__t07;\
	__t06 =__rt+__t0c;			__t07 =__it+__t0d;\
	__t0c =__rt-__t0c;			__t0d =__it-__t0d;\
	__t00 =__t00+__t06;			__t01 =__t01+__t07;\
	__A0r =__t00;				__A0i =__t01;\
	__t06 =__t00+c3m1*__t06;	__t07 =__t01+c3m1*__t07;\
	__rt  =s3*__t0c;			__it  =s3*__t0d;\
	__A3r =__t06+__it;			__A3i =__t07-__rt;\
	__A6r =__t06-__it;			__A6i =__t07+__rt;\
\
	__rt  =__t08*c +__t09*s;	__it  =__t09*c -__t08*s;\
	__tt  =__t0e*c2+__t0f*s2;	__t0f =__t0f*c2-__t0e*s2;	__t0e=__tt;\
	__t08 =__rt+__t0e;			__t09 =__it+__t0f;\
	__t0e =__rt-__t0e;			__t0f =__it-__t0f;\
	__t02 =__t02+__t08;			__t03 =__t03+__t09;\
	__A7r =__t02;				__A7i =__t03;\
	__t08 =__t02+c3m1*__t08;	__t09 =__t03+c3m1*__t09;\
	__rt  =s3*__t0e;			__it  =s3*__t0f;\
	__A1r =__t08+__it;			__A1i =__t09-__rt;\
	__A4r =__t08-__it;			__A4i =__t09+__rt;\
\
	__rt  =__t0a*c2+__t0b*s2;	__it  =__t0b*c2-__t0a*s2;\
	__tt  =__t0g*c4+__t0h*s4;	__t0h =__t0h*c4-__t0g*s4;	__t0g=__tt;\
	__t0a =__rt+__t0g;			__t0b =__it+__t0h;\
	__t0g =__rt-__t0g;			__t0h =__it-__t0h;\
	__t04 =__t04+__t0a;			__t05 =__t05+__t0b;\
	__A5r =__t04;				__A5i =__t05;\
	__t0a =__t04+c3m1*__t0a;	__t0b =__t05+c3m1*__t0b;\
	__rt  =s3*__t0g;			__it  =s3*__t0h;\
	__A8r =__t0a+__it;			__A8i =__t0b-__rt;\
	__A2r =__t0a-__it;			__A2i =__t0b+__rt;\
}

/****** RADIX = 11: ******/

/*...Simple-algo Radix-11 DFT.  Totals: 140 ADD, 100 MUL	*/
#define RADIX_11_DFT_BASIC(__A0r,__A0i,__A1r,__A1i,__A2r,__A2i,__A3r,__A3i,__A4r,__A4i,__A5r,__A5i,__A6r,__A6i,__A7r,__A7i,__A8r,__A8i,__A9r,__A9i,__Aar,__Aai\
					,__B0r,__B0i,__B1r,__B1i,__B2r,__B2i,__B3r,__B3i,__B4r,__B4i,__B5r,__B5i,__B6r,__B6i,__B7r,__B7i,__B8r,__B8i,__B9r,__B9i,__Bar,__Bai\
					)\
{\
	double t0r,t0i,t1r,t1i,t2r,t2i,t3r,t3i,t4r,t4i,t5r,t5i,t6r,t6i,t7r,t7i,t8r,t8i,t9r,t9i,tar,tai;\
	double cr1,cr2,cr3,cr4,cr5,ci1,ci2,ci3,ci4,ci5,sr1,sr2,sr3,sr4,sr5,si1,si2,si3,si4,si5;\
\
	t0r = __A0r       ;			t0i = __A0i        ;	/* x0		*/\
	t1r = __A1r +__Aar;			t1i = __A1i + __Aai;	/* x1 + x10	*/\
	t2r = __A2r +__A9r;			t2i = __A2i + __A9i;	/* x2 + x9	*/\
	t3r = __A3r +__A8r;			t3i = __A3i + __A8i;	/* x3 + x8	*/\
	t4r = __A4r +__A7r;			t4i = __A4i + __A7i;	/* x4 + x7	*/\
	t5r = __A5r +__A6r;			t5i = __A5i + __A6i;	/* x5 + x6	*/\
	t6r = __A5r -__A6r;			t6i = __A5i - __A6i;	/* x5 - x6	*/\
	t7r = __A4r -__A7r;			t7i = __A4i - __A7i;	/* x4 - x7	*/\
	t8r = __A3r -__A8r;			t8i = __A3i - __A8i;	/* x3 - x8	*/\
	t9r = __A2r -__A9r;			t9i = __A2i - __A9i;	/* x2 - x9	*/\
	tar = __A1r -__Aar;			tai = __A1i - __Aai;	/* x1 - x10	*/\
\
	cr1= t0r+cc1*t1r+cc2*t2r+cc3*t3r+cc4*t4r+cc5*t5r;	ci1= t0i+cc1*t1i+cc2*t2i+cc3*t3i+cc4*t4i+cc5*t5i;	/* C1	*/\
	cr2= t0r+cc2*t1r+cc4*t2r+cc5*t3r+cc3*t4r+cc1*t5r;	ci2= t0i+cc2*t1i+cc4*t2i+cc5*t3i+cc3*t4i+cc1*t5i;	/* C2	*/\
	cr3= t0r+cc3*t1r+cc5*t2r+cc2*t3r+cc1*t4r+cc4*t5r;	ci3= t0i+cc3*t1i+cc5*t2i+cc2*t3i+cc1*t4i+cc4*t5i;	/* C3	*/\
	cr4= t0r+cc4*t1r+cc3*t2r+cc1*t3r+cc5*t4r+cc2*t5r;	ci4= t0i+cc4*t1i+cc3*t2i+cc1*t3i+cc5*t4i+cc2*t5i;	/* C4	*/\
	cr5= t0r+cc5*t1r+cc1*t2r+cc4*t3r+cc2*t4r+cc3*t5r;	ci5= t0i+cc5*t1i+cc1*t2i+cc4*t3i+cc2*t4i+cc3*t5i;	/* C5	*/\
\
	sr1=     ss1*tar+ss2*t9r+ss3*t8r+ss4*t7r+ss5*t6r;	si1=     ss1*tai+ss2*t9i+ss3*t8i+ss4*t7i+ss5*t6i;	/* S1	*/\
	sr2=     ss2*tar+ss4*t9r-ss5*t8r-ss3*t7r-ss1*t6r;	si2=     ss2*tai+ss4*t9i-ss5*t8i-ss3*t7i-ss1*t6i;	/* S2	*/\
	sr3=     ss3*tar-ss5*t9r-ss2*t8r+ss1*t7r+ss4*t6r;	si3=     ss3*tai-ss5*t9i-ss2*t8i+ss1*t7i+ss4*t6i;	/* S3	*/\
	sr4=     ss4*tar-ss3*t9r+ss1*t8r+ss5*t7r-ss2*t6r;	si4=     ss4*tai-ss3*t9i+ss1*t8i+ss5*t7i-ss2*t6i;	/* S4	*/\
	sr5=     ss5*tar-ss1*t9r+ss4*t8r-ss2*t7r+ss3*t6r;	si5=     ss5*tai-ss1*t9i+ss4*t8i-ss2*t7i+ss3*t6i;	/* S5	*/\
\
	__B0r = t0r+t1r+t2r+t3r+t4r+t5r;__B0i = t0i+t1i+t2i+t3i+t4i+t5i;	/* X0	*/\
	__B1r = cr1 - si1;				__B1i = ci1 + sr1;	/* X1 = C1 + I*S1	*/\
	__B2r = cr2 - si2;				__B2i = ci2 + sr2;	/* X2 = C2 + I*S2	*/\
	__B3r = cr3 - si3;				__B3i = ci3 + sr3;	/* X3 = C3 + I*S3	*/\
	__B4r = cr4 - si4;				__B4i = ci4 + sr4;	/* X4 = C4 + I*S4	*/\
	__B5r = cr5 - si5;				__B5i = ci5 + sr5;	/* X5 = C5 + I*S5	*/\
	__B6r = cr5 + si5;				__B6i = ci5 - sr5;	/* X6 =	C5 - I*S5	*/\
	__B7r = cr4 + si4;				__B7i = ci4 - sr4;	/* X7 =	C4 - I*S4	*/\
	__B8r = cr3 + si3;				__B8i = ci3 - sr3;	/* X8 =	C3 - I*S3	*/\
	__B9r = cr2 + si2;				__B9i = ci2 - sr2;	/* X9 =	C2 - I*S2	*/\
	__Bar = cr1 + si1;				__Bai = ci1 - sr1;	/* X10=	C1 - I*S1	*/\
}

/*...Radix-11 DFT using length-5 cyclic convolution scheme for the 5 x 5 matrix submultiplies.  Totals: 160 ADD,  44 MUL	*/
#define RADIX_11_DFT(__A0r,__A0i,__A1r,__A1i,__A2r,__A2i,__A3r,__A3i,__A4r,__A4i,__A5r,__A5i,__A6r,__A6i,__A7r,__A7i,__A8r,__A8i,__A9r,__A9i,__Aar,__Aai\
					,__B0r,__B0i,__B1r,__B1i,__B2r,__B2i,__B3r,__B3i,__B4r,__B4i,__B5r,__B5i,__B6r,__B6i,__B7r,__B7i,__B8r,__B8i,__B9r,__B9i,__Bar,__Bai\
					)\
{\
	double t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22;\
	double cr1,cr2,cr3,cr4,cr5,ci1,ci2,ci3,ci4,ci5,sr1,sr2,sr3,sr4,sr5,si1,si2,si3,si4,si5;\
	double c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10;\
\
	__B0r = __A0r;				__B0i = __A0i;			/* x0		; store in t2 */\
	t3  = __A1r +__Aar;			t4  = __A1i + __Aai;	/* x1 + x10	; load a1-a4,a7-aa into xmm0-7, compute t4-t10 (save in xmm0-3) */\
	t5  = __A2r +__A9r;			t6  = __A2i + __A9i;	/* x2 + x9	; and t16-22 (dump to memory), load a0,5,6 into xmm4-6, compute t12,14, */\
	t7  = __A3r +__A8r;			t8  = __A3i + __A8i;	/* x3 + x8	; dump t14 to memory. */\
	t9  = __A4r +__A7r;			t10 = __A4i + __A7i;	/* x4 + x7	*/\
	t11 = __A5r +__A6r;			t12 = __A5i + __A6i;	/* x5 + x6	*/\
	t13 = __A5r -__A6r;			t14 = __A5i - __A6i;	/* x5 - x6	*/\
	t15 = __A4r -__A7r;			t16 = __A4i - __A7i;	/* x4 - x7	*/\
	t17 = __A3r -__A8r;			t18 = __A3i - __A8i;	/* x3 - x8	*/\
	t19 = __A2r -__A9r;			t20 = __A2i - __A9i;	/* x2 - x9	*/\
	t21 = __A1r -__Aar;			t22 = __A1i - __Aai;	/* x1 - x10	*/\
\
/*  Here are the 5 cosine terms:\
	let y0 = (x1+x10) = t3,4, y4 = (x2+x9) = t5,6, y2 = (x3+x8) = t7,8, y3 = (x4+x7) = t9,10, y1 = (x5+x6) = t11,12, then form\
*/\
	/* b0/t2,t4,t6,t8,t10,t12 in xmm0-5, xmm6,7 free */\
	c1 = t3-t5;					s1 = t4-t6;		/* store in s1/t4,  make 2 copies (call them t0,t14) */\
	c2 = t11-t5;				s2 = t12-t6;	/* store in s2/t12  */\
	c4 = t7-t5;					s4 = t8-t6;		/* store in s4/t8 */\
	c5 = t9-t5;					s5 = t10-t6;	/* store in s5/t10, then mul t6*5 */\
	c3 = c1+c2;					s3 = s1+s2;	/* t4+t12-2*t6; store in s3/t0 */\
	c7 = c1-c4;					s7 = s1-s4;	/* t4-t8; store in t14, mul by a6 and store. */\
	c6 = c4+c5;					s6 = s4+s5;	/* t8+t10-2*t6; store in t8 (do after s1-s4, since this will overwrite s4) */\
	c8 = c2-c5;					s8 = s2-s5;	/* t12-t10; copy s2 to t14, mul s2-s5 by a7 and store. */\
	c9 = c3-c6;					s9 = s3-s6;	/* t4+t12-(t8+t10); copy s3 to t14, store result in t14. */\
	c10= c3+c6+5*t5;			s10= s3+s6+5*t6;	/* c10= t3+t5+t7+t9+t11;	 s10= t4+t6+t8+t10+t12; store in t6. */\
\
	__B0r += c10;				__B0i += s10;	/* X0/t2, store result to memory */\
\
	c3 = a2*c3;					s3 = a2*s3;	/* t0  */\
	c6 = a5*c6;					s6 = a5*s6;	/* t8  */\
	c9 = a8*c9;					s9 = a8*s9;	/* t14 */\
	c10= a9*c10+__B0r;			s10= a9*s10+__B0i;	/* t6  */\
\
	c1 = a0*c1+c3;				s1 = a0*s1+s3;\
	c2 = a1*c2+c3;				s2 = a1*s2+s3;\
	c4 = a3*c4+c6;				s4 = a3*s4+s6;\
	c5 = a4*c5+c6;				s5 = a4*s5+s6;\
	c7 = a6*c7+c9;				s7 = a6*s7+s9;\
	c8 = a7*c8+c9;				s8 = a7*s8+s9;\
\
	cr1 = c10+c1-c7;			ci1 = s10+s1-s7;\
	cr2 = c10-c1-c2-c4-c5;		ci2 = s10-s1-s2-s4-s5;\
	cr3 = c10+c4+c7;			ci3 = s10+s4+s7;\
	cr4 = c10+c5+c8;			ci4 = s10+s5+s8;\
	cr5 = c10+c2-c8;			ci5 = s10+s2-s8;\
\
/*  Here are the 5 sine terms:\
let y0 = (x1-x10) = t21,22, y4 = (x9-x2) = -t19,20, y2 = (x3-x8) = t17,18, y3 = (x4-x7) = t15,16, y1 = (x5-x6) = t13,t14, then form\
*/\
	c1 = t21+t19;				s1 = t22+t20;\
	c2 = t13+t19;				s2 = t14+t20;\
	c3 = c1+c2;					s3 = s1+s2;\
	c4 = t17+t19;				s4 = t18+t20;\
	c5 = t15+t19;				s5 = t16+t20;\
	c6 = c4+c5;					s6 = s4+s5;\
	c7 = c1-c4;					s7 = s1-s4;\
	c8 = c2-c5;					s8 = s2-s5;\
	c9 = c3-c6;					s9 = s3-s6;\
	c10= c3+c6-5*t19;			s10= s3+s6-5*t20;	/* c10= t21-t19+t17+t15+t13;	s10= t22-t20+t18+t16+t14; */\
\
	c3 = b2*c3;					s3 = b2*s3;\
	c6 = b5*c6;					s6 = b5*s6;\
	c9 = b8*c9;					s9 = b8*s9;\
	c10= b9*c10;				s10= b9*s10;\
\
	c1 = b0*c1+c3;				s1 = b0*s1+s3;\
	c2 = b1*c2+c3;				s2 = b1*s2+s3;\
	c3 = b3*c4+c6;				s3 = b3*s4+s6;\
	c4 = b4*c5+c6;				s4 = b4*s5+s6;\
	c5 = b6*c7+c9;				s5 = b6*s7+s9;\
	c6 = b7*c8+c9;				s6 = b7*s8+s9;\
\
	sr1 = c10+c1-c5;			si1 = s10+s1-s5;\
	sr2 = c1+c2+c3+c4-c10;		si2 = s1+s2+s3+s4-s10;\
	sr3 = c10+c3+c5;			si3 = s10+s3+s5;\
	sr4 = c10+c4+c6;			si4 = s10+s4+s6;\
	sr5 = c10+c2-c6;			si5 = s10+s2-s6;\
\
	__B1r = cr1 - si1;				__B1i = ci1 + sr1;	/* X1 = C1 + I*S1	*/\
	__B2r = cr2 - si2;				__B2i = ci2 + sr2;	/* X2 = C2 + I*S2	*/\
	__B3r = cr3 - si3;				__B3i = ci3 + sr3;	/* X3 = C3 + I*S3	*/\
	__B4r = cr4 - si4;				__B4i = ci4 + sr4;	/* X4 = C4 + I*S4	*/\
	__B5r = cr5 - si5;				__B5i = ci5 + sr5;	/* X5 = C5 + I*S5	*/\
	__B6r = cr5 + si5;				__B6i = ci5 - sr5;	/* X6 =	C5 - I*S5	*/\
	__B7r = cr4 + si4;				__B7i = ci4 - sr4;	/* X7 =	C4 - I*S4	*/\
	__B8r = cr3 + si3;				__B8i = ci3 - sr3;	/* X8 =	C3 - I*S3	*/\
	__B9r = cr2 + si2;				__B9i = ci2 - sr2;	/* X9 =	C2 - I*S2	*/\
	__Bar = cr1 + si1;				__Bai = ci1 - sr1;	/* X10=	C1 - I*S1	*/\
}

/****** RADIX = 13: ******/

/* This version implements a straight tangent-based radix-13 DFT and only tries to minimize the number
of temporaries used in each of the major phases, with view toward a 16-register-optimized version.
(That is why we do the x and y-pieces side-by-side).
For an 8-register-optimized version which moreover mimics x86-style destructive (2-input, one overwritten with output)
register arithmetic, see the RADIX_13_XYZ version of the same macro in radix13_ditN_cy_dif1.c .

Output ordering is as for DIF, so for DIT caller must reverse order of last 12 B-output args.

ASSUMES that the needed DC and DS-constants have been defined in the calling routine via inclusion of the radix13.h header file..
*/
#define RADIX_13_DFT(__A0r,__A0i,__A1r,__A1i,__A2r,__A2i,__A3r,__A3i,__A4r,__A4i,__A5r,__A5i,__A6r,__A6i,__A7r,__A7i,__A8r,__A8i,__A9r,__A9i,__Aar,__Aai,__Abr,__Abi,__Acr,__Aci\
					,__B0r,__B0i,__B1r,__B1i,__B2r,__B2i,__B3r,__B3i,__B4r,__B4i,__B5r,__B5i,__B6r,__B6i,__B7r,__B7i,__B8r,__B8i,__B9r,__B9i,__Bar,__Bai,__Bbr,__Bbi,__Bcr,__Bci\
					)\
{\
	double xr0, xi0, xr1, xi1, xr2, xi2, xr3, xi3, xr4, xi4, xr5, xi5, xr6, xi6, xr7, xi7;\
	double yr1, yi1, yr2, yi2, yr3, yi3, yr4, yi4, yr5, yi5, yr6, yi6;\
\
		xr7 = __A6r + __A7r;			xi7 = __A6i + __A7i;\
		xr5 = __A5r + __A8r;			xi5 = __A5i + __A8i;\
		xr4 = __A4r + __A9r;			xi4 = __A4i + __A9i;\
		xr6 = __A3r + __Aar;			xi6 = __A3i + __Aai;\
		xr3 = __A2r + __Abr;			xi3 = __A2i + __Abi;\
		xr1 = __A1r + __Acr;			xi1 = __A1i + __Aci;\
\
		yr1 = __A1r - __Acr;			yi1 = __A1i - __Aci;\
		yr2 = __A2r - __Abr;			yi2 = __A2i - __Abi;\
		yr5 = __A3r - __Aar;			yi5 = __A3i - __Aai;\
		yr3 = __A4r - __A9r;			yi3 = __A4i - __A9i;\
		yr4 = __A5r - __A8r;			yi4 = __A5i - __A8i;\
		yr6 = __A6r - __A7r;			yi6 = __A6i - __A7i;\
\
/* x-terms need 8 registers for each side: */\
		xr2 = xr1+xr5;					xi2 = xi1+xi5;\
		xr5 = xr1-xr5;					xi5 = xi1-xi5;\
		xr0 = xr3+xr6;					xi0 = xi3+xi6;\
		xr6 = xr3-xr6;					xi6 = xi3-xi6;\
		xr3 = xr4+xr7;					xi3 = xi4+xi7;\
		xr7 = xr4-xr7;					xi7 = xi4-xi7;\
		xr1 = xr2+xr0+xr3;				xi1 = xi2+xi0+xi3;\
		xr4 = __A0r+xr1;				xi4 = __A0i+xi1;\
		__B0r = xr4;					__B0i = xi4;\
		xr4 +=    DC1*xr1;				xi4 +=    DC1*xi1;\
		xr1 = xr4+DC2*xr2+DC3*xr0;		xi1 = xi4+DC2*xi2+DC3*xi0;\
		xr2 = xr4+DC2*xr3+DC3*xr2;		xi2 = xi4+DC2*xi3+DC3*xi2;\
		xr0 = xr4+DC2*xr0+DC3*xr3;		xi0 = xi4+DC2*xi0+DC3*xi3;\
		xr4 = DC4*xr5+DC5*xr6+DC6*xr7;	xi4 = DC4*xi5+DC5*xi6+DC6*xi7;\
		xr3 = DC4*xr6+DC5*xr7-DC6*xr5;	xi3 = DC4*xi6+DC5*xi7-DC6*xi5;\
		xr7 = DC4*xr7-DC5*xr5-DC6*xr6;	xi7 = DC4*xi7-DC5*xi5-DC6*xi6;\
		xr5 = xr0+xr4;					xi5 = xi0+xi4;\
		xr0 = xr0-xr4;					xi0 = xi0-xi4;\
		xr6 = xr2+xr3;					xi6 = xi2+xi3;\
		xr2 = xr2-xr3;					xi2 = xi2-xi3;\
		xr3 = xr1+xr7;					xi3 = xi1+xi7;\
		xr7 = xr1-xr7;					xi7 = xi1-xi7;\
/* x4,1 free...now do y-terms: */\
		xr4 = yr1-yr3+yr5;				xi4 = yi1-yi3+yi5;\
		yr1 = yr1+yr3;					yi1 = yi1+yi3;\
		yr5 = yr3+yr5;		 			yi5 = yi3+yi5;\
		yr3 = yr2+yr4+yr6;				yi3 = yi2+yi4+yi6;\
		yr2 = yr2-yr4;					yi2 = yi2-yi4;\
		yr4 = yr6-yr4;		 			yi4 = yi6-yi4;\
		xr1 = DS1*xr4+DS2*yr3;			xi1 = DS1*xi4+DS2*yi3;\
/* Need to use one output as a temporary here if have just 8 registers: */\
		__B1i = DS1*yr3-DS2*xr4;		__B1r = DS1*yi3-DS2*xi4;\
		yr3 = yr1+yr2;		 			yi3 = yi1+yi2;\
		yr6 = yr5+yr4;					yi6 = yi5+yi4;\
		xr4 = DS3*yr3-DS6*yr6;			xi4 = DS3*yi3-DS6*yi6;\
		yr3 = DS9*yr3-DS3*yr6;			yi3 = DS9*yi3-DS3*yi6;\
		yr6 = xr4+DS4*yr2-DS7*yr4;		yi6 = xi4+DS4*yi2-DS7*yi4;\
		yr2 = yr3+DSa*yr2-DS4*yr4;		yi2 = yi3+DSa*yi2-DS4*yi4;\
		yr4 = xr4+DS5*yr1-DS8*yr5;		yi4 = xi4+DS5*yi1-DS8*yi5;\
		yr3 = yr3+DSb*yr1-DS5*yr5;		yi3 = yi3+DSb*yi1-DS5*yi5;\
		xr4 = __B1i;					xi4 = __B1r;\
		yr5 = xr1+yr6;					yi5 = xi1+yi6;\
		yr6 = yr6-xr1-yr2;				yi6 = yi6-xi1-yi2;\
		yr2 = xr1-yr2;					yi2 = xi1-yi2;\
		yr1 = xr4+yr4;					yi1 = xi4+yi4;\
		yr4 = yr4-xr4-yr3;				yi4 = yi4-xi4-yi3;\
		xr4 = xr4-yr3;					xi4 = xi4-yi3;\
/* In ASM, do xr and yi-terms and combine to get real parts of outputs,\
   then do xi and yr-terms and combine to get imaginary parts: */\
		__B1r = xr7-xi4;				__B1i = xi7+xr4;\
		__Bcr = xr7+xi4;	 			__Bci = xi7-xr4;\
		__B2r = xr2-yi2;	 			__B2i = xi2+yr2;\
		__Bbr = xr2+yi2;				__Bbi = xi2-yr2;\
		__B3r = xr6-yi1;	 			__B3i = xi6+yr1;\
		__Bar = xr6+yi1;				__Bai = xi6-yr1;\
		__B4r = xr0-yi4;	 			__B4i = xi0+yr4;\
		__B9r = xr0+yi4;				__B9i = xi0-yr4;\
		__B5r = xr3+yi6;				__B5i = xi3-yr6;\
		__B8r = xr3-yi6;				__B8i = xi3+yr6;\
		__B6r = xr5-yi5;				__B6i = xi5+yr5;\
		__B7r = xr5+yi5;				__B7i = xi5-yr5;\
/* Totals: 164 ADD, 64 MUL. */\
}

/****** RADIX = 15: ******/

/* Totals: 162 ADD, 50 MUL	*/
#define RADIX_15_DIF(\
	__Ar00,__Ai00,\
	__Ar01,__Ai01,\
	__Ar02,__Ai02,\
	__Ar03,__Ai03,\
	__Ar04,__Ai04,\
	__Ar05,__Ai05,\
	__Ar06,__Ai06,\
	__Ar07,__Ai07,\
	__Ar08,__Ai08,\
	__Ar09,__Ai09,\
	__Ar10,__Ai10,\
	__Ar11,__Ai11,\
	__Ar12,__Ai12,\
	__Ar13,__Ai13,\
	__Ar14,__Ai14,\
	__tr00,__ti00,\
	__tr01,__ti01,\
	__tr02,__ti02,\
	__tr03,__ti03,\
	__tr04,__ti04,\
	__tr05,__ti05,\
	__tr06,__ti06,\
	__tr07,__ti07,\
	__tr08,__ti08,\
	__tr09,__ti09,\
	__tr10,__ti10,\
	__tr11,__ti11,\
	__tr12,__ti12,\
	__tr13,__ti13,\
	__tr14,__ti14,\
	__Br00,__Bi00,\
	__Br01,__Bi01,\
	__Br02,__Bi02,\
	__Br03,__Bi03,\
	__Br04,__Bi04,\
	__Br05,__Bi05,\
	__Br06,__Bi06,\
	__Br07,__Bi07,\
	__Br08,__Bi08,\
	__Br09,__Bi09,\
	__Br10,__Bi10,\
	__Br11,__Bi11,\
	__Br12,__Bi12,\
	__Br13,__Bi13,\
	__Br14,__Bi14,\
	__rt,__it)\
{\
/* Default input order for a twiddleless radix-15 DIF (cf.radix15_dif_pass1()) is [0,10,5,12,7,2,9,4,14,6,1,11,3,13,8]. */\
/*...gather the needed data (15 64-bit complex, i.e. 30 64-bit reals) and do 3 radix-5 transforms...*/\
/*...Block 1: Swap inputs as 03 <-> 12, 06 <-> 09: */\
	__tr00 = __Ar00;				__ti00 = __Ai00;\
	__tr01 = __Ar12;				__ti01 = __Ai12;\
	__rt   = __Ar03;				__it   = __Ai03;\
	__tr03 = __tr01-__rt;			__ti03 = __ti01-__it;\
	__tr01 = __tr01+__rt;			__ti01 = __ti01+__it;\
	__tr02 = __Ar09;				__ti02 = __Ai09;\
	__rt   = __Ar06;				__it   = __Ai06;\
	__tr04 = __tr02-__rt;			__ti04 = __ti02-__it;\
	__tr02 = __tr02+__rt;			__ti02 = __ti02+__it;\
\
	__rt   = __tr01+__tr02;			__it   = __ti01+__ti02;\
	__tr00 = __tr00+__rt;			__ti00 = __ti00+__it;\
	__rt   = __tr00+cn1*__rt;		__it   = __ti00+cn1*__it;\
	__tr02 = cn2*(__tr01-__tr02);	__ti02 = cn2*(__ti01-__ti02);\
	__tr01 = __rt+__tr02;			__ti01 = __it+__ti02;\
	__tr02 = __rt-__tr02;			__ti02 = __it-__ti02;\
	__rt   = ss3*(__tr03-__tr04);	__it   = ss3*(__ti03-__ti04);\
	__tr04 = __rt+sn1*__tr04;		__ti04 = __it+sn1*__ti04;\
	__tr03 = __rt-sn2*__tr03;		__ti03 = __it-sn2*__ti03;\
	__rt   = __tr04;				__it   = __ti04;\
	__tr04 = __tr01+__it;			__ti04 = __ti01-__rt;		/*<==prefer these to be stored in __tr,i04 */\
	__tr01 = __tr01-__it;			__ti01 = __ti01+__rt;\
	__rt   = __tr03;				__it   = __ti03;\
	__tr03 = __tr02+__it;			__ti03 = __ti02-__rt;		/*<==prefer these to be stored in __tr,i03 */\
	__tr02 = __tr02-__it;			__ti02 = __ti02+__rt;\
\
/*...Block 2: Swap inputs as 01 <-> 10, 04 <-> 07: */\
	__tr05 = __Ar10;				__ti05 = __Ai10;\
	__tr06 = __Ar07;				__ti06 = __Ai07;\
	__rt   = __Ar13;				__it   = __Ai13;\
	__tr08 = __tr06 -__rt;			__ti08 = __ti06 -__it;\
	__tr06 = __tr06 +__rt;			__ti06 = __ti06 +__it;\
	__tr07 = __Ar04;				__ti07 = __Ai04;\
	__rt   = __Ar01;				__it   = __Ai01;\
	__tr09 = __tr07 -__rt;			__ti09 = __ti07 -__it;\
	__tr07 = __tr07 +__rt;			__ti07 = __ti07 +__it;\
\
	__rt   = __tr06+__tr07;			__it   = __ti06+__ti07;\
	__tr05 = __tr05+__rt;			__ti05 = __ti05+__it;\
	__rt   = __tr05+cn1*__rt;		__it   = __ti05+cn1*__it;\
	__tr07 = cn2*(__tr06-__tr07);	__ti07 = cn2*(__ti06-__ti07);\
	__tr06 = __rt+__tr07;			__ti06 = __it+__ti07;\
	__tr07 = __rt-__tr07;			__ti07 = __it-__ti07;\
	__rt   = ss3*(__tr08-__tr09);	__it   = ss3*(__ti08-__ti09);\
	__tr09 = __rt+sn1*__tr09;		__ti09 = __it+sn1*__ti09;\
	__tr08 = __rt-sn2*__tr08;		__ti08 = __it-sn2*__ti08;\
	__rt   = __tr09;				__it   = __ti09;\
	__tr09 = __tr06+__it;			__ti09 = __ti06-__rt;\
	__tr06 = __tr06-__it;			__ti06 = __ti06+__rt;\
	__rt   = __tr08;				__it   = __ti08;\
	__tr08 = __tr07+__it;			__ti08 = __ti07-__rt;\
	__tr07 = __tr07-__it;			__ti07 = __ti07+__rt;\
\
/*...Block 3: Swap inputs as 02 <-> 05, 08 <-> 14: */\
	__tr10 = __Ar05;				__ti10 = __Ai05;\
	__tr11 = __Ar02;				__ti11 = __Ai02;\
	__rt   = __Ar08;				__it   = __Ai08;\
	__tr13 = __tr11 -__rt;			__ti13 = __ti11 -__it;\
	__tr11 = __tr11 +__rt;			__ti11 = __ti11 +__it;\
	__tr12 = __Ar14;				__ti12 = __Ai14;\
	__rt   = __Ar11;				__it   = __Ai11;\
	__tr14 = __tr12 -__rt;			__ti14 = __ti12 -__it;\
	__tr12 = __tr12 +__rt;			__ti12 = __ti12 +__it;\
\
	__rt   = __tr11+__tr12;			__it   = __ti11+__ti12;\
	__tr10 = __tr10+__rt;			__ti10 = __ti10+__it;\
	__rt   = __tr10+cn1*__rt;		__it   = __ti10+cn1*__it;\
	__tr12 = cn2*(__tr11-__tr12);	__ti12 = cn2*(__ti11-__ti12);\
	__tr11 = __rt+__tr12;			__ti11 = __it+__ti12;\
	__tr12 = __rt-__tr12;			__ti12 = __it-__ti12;\
	__rt   = ss3*(__tr13-__tr14);	__it   = ss3*(__ti13-__ti14);\
	__tr14 = __rt+sn1*__tr14;		__ti14 = __it+sn1*__ti14;\
	__tr13 = __rt-sn2*__tr13;		__ti13 = __it-sn2*__ti13;\
	__rt   = __tr14;				__it   = __ti14;\
	__tr14 = __tr11+__it;			__ti14 = __ti11-__rt;\
	__tr11 = __tr11-__it;			__ti11 = __ti11+__rt;\
	__rt   = __tr13;				__it   = __ti13;\
	__tr13 = __tr12+__it;			__ti13 = __ti12-__rt;\
	__tr12 = __tr12-__it;			__ti12 = __ti12+__rt;\
\
/*...and now do five radix-3 transforms.\
	The required output permutation is [0,1,2,13,14,12,9,10,11,8,6,7,4,5,3]. */\
/*...Block 1: */\
	__rt   = __tr10;				__it   = __ti10;\
	__tr10 = __tr05-__rt;			__ti10 = __ti05-__it;\
	__tr05 = __tr05+__rt;			__ti05 = __ti05+__it;\
	__tr00 = __tr00+__tr05;			__ti00 = __ti00+__ti05;\
	__Br00 = __tr00;				__Bi00 = __ti00;\
	__tr05 = __tr00+c3m1*__tr05;	__ti05 = __ti00+c3m1*__ti05;\
	__rt   = s*__tr10;				__it   = s*__ti10;\
	__Br01 = __tr05-__it;			__Bi01 = __ti05+__rt;\
	__Br02 = __tr05+__it;			__Bi02 = __ti05-__rt;\
\
/*...Block 2: */\
	__rt   = __tr11;				__it   = __ti11;\
	__tr11 = __tr06-__rt;			__ti11 = __ti06-__it;\
	__tr06 = __tr06+__rt;			__ti06 = __ti06+__it;\
	__tr01 = __tr01+__tr06;			__ti01 = __ti01+__ti06;\
	__Br13 = __tr01;				__Bi13 = __ti01;\
	__tr06 = __tr01+c3m1*__tr06;	__ti06 = __ti01+c3m1*__ti06;\
	__rt   = s*__tr11;				__it   = s*__ti11;\
	__Br14 = __tr06-__it;			__Bi14 = __ti06+__rt;\
	__Br12 = __tr06+__it;			__Bi12 = __ti06-__rt;\
\
/*...Block 3: */\
	__rt   = __tr12;				__it   = __ti12;\
	__tr12 = __tr07-__rt;			__ti12 = __ti07-__it;\
	__tr07 = __tr07+__rt;			__ti07 = __ti07+__it;\
	__tr02 = __tr02+__tr07;			__ti02 = __ti02+__ti07;\
	__Br09 = __tr02;				__Bi09 = __ti02;\
	__tr07 = __tr02+c3m1*__tr07;	__ti07 = __ti02+c3m1*__ti07;\
	__rt   = s*__tr12;				__it   = s*__ti12;\
	__Br10 = __tr07-__it;			__Bi10 = __ti07+__rt;\
	__Br11 = __tr07+__it;			__Bi11 = __ti07-__rt;\
\
/*...Block 4: */\
	__rt   = __tr13;				__it   = __ti13;\
	__tr13 = __tr08-__rt;			__ti13 = __ti08-__it;\
	__tr08 = __tr08+__rt;			__ti08 = __ti08+__it;\
	__tr03 = __tr03+__tr08;			__ti03 = __ti03+__ti08;\
	__Br08 = __tr03;				__Bi08 = __ti03;\
	__tr08 = __tr03+c3m1*__tr08;	__ti08 = __ti03+c3m1*__ti08;\
	__rt   = s*__tr13;				__it   = s*__ti13;\
	__Br06 = __tr08-__it;			__Bi06 = __ti08+__rt;\
	__Br07 = __tr08+__it;			__Bi07 = __ti08-__rt;\
\
/*...Block 5: */\
	__rt   = __tr14;				__it   = __ti14;\
	__tr14 = __tr09-__rt;			__ti14 = __ti09-__it;\
	__tr09 = __tr09+__rt;			__ti09 = __ti09+__it;\
	__tr04 = __tr04+__tr09;			__ti04 = __ti04+__ti09;\
	__Br04 = __tr04;				__Bi04 = __ti04;\
	__tr09 = __tr04+c3m1*__tr09;	__ti09 = __ti04+c3m1*__ti09;\
	__rt   = s*__tr14;				__it   = s*__ti14;\
	__Br05 = __tr09-__it;			__Bi05 = __ti09+__rt;\
	__Br03 = __tr09+__it;			__Bi03 = __ti09-__rt;\
}

/* Totals: 162 ADD, 50 MUL	*/
#define RADIX_15_DIT(\
	__Ar00,__Ai00,\
	__Ar01,__Ai01,\
	__Ar02,__Ai02,\
	__Ar03,__Ai03,\
	__Ar04,__Ai04,\
	__Ar05,__Ai05,\
	__Ar06,__Ai06,\
	__Ar07,__Ai07,\
	__Ar08,__Ai08,\
	__Ar09,__Ai09,\
	__Ar10,__Ai10,\
	__Ar11,__Ai11,\
	__Ar12,__Ai12,\
	__Ar13,__Ai13,\
	__Ar14,__Ai14,\
	__tr00,__ti00,\
	__tr01,__ti01,\
	__tr02,__ti02,\
	__tr03,__ti03,\
	__tr04,__ti04,\
	__tr05,__ti05,\
	__tr06,__ti06,\
	__tr07,__ti07,\
	__tr08,__ti08,\
	__tr09,__ti09,\
	__tr10,__ti10,\
	__tr11,__ti11,\
	__tr12,__ti12,\
	__tr13,__ti13,\
	__tr14,__ti14,\
	__Br00,__Bi00,\
	__Br01,__Bi01,\
	__Br02,__Bi02,\
	__Br03,__Bi03,\
	__Br04,__Bi04,\
	__Br05,__Bi05,\
	__Br06,__Bi06,\
	__Br07,__Bi07,\
	__Br08,__Bi08,\
	__Br09,__Bi09,\
	__Br10,__Bi10,\
	__Br11,__Bi11,\
	__Br12,__Bi12,\
	__Br13,__Bi13,\
	__Br14,__Bi14,\
	__rt,__it)\
{\
/* Default input order for a twiddleless radix-15 DIT (cf.radix15_dit_pass1()) is [0,2,1,8,7,6,13,12,14,4,3,5,9,11,10]. */\
/*...gather the needed data (15 64-bit complex, i.e. 30 64-bit reals) and do 5 radix-3 transforms...*/\
/*...Block 1: Swap inputs 01 <-> 02: */\
	    __tr00 = __Ar00;				__ti00 = __Ai00;\
	    __tr01 = __Ar02;				__ti01 = __Ai02;\
	    __rt   = __Ar01;				__it   = __Ai01;\
	    __tr02 = __tr01-__rt;			__ti02 = __ti01-__it;\
	    __tr01 = __tr01+__rt;			__ti01 = __ti01+__it;\
	    __tr00 = __tr00+__tr01;			__ti00 = __ti00+__ti01;\
	    __tr01 = __tr00+c3m1*__tr01;	__ti01 = __ti00+c3m1*__ti01;\
	    __rt   = s*__tr02;				__it   = s*__ti02;\
	    __tr02 = __tr01-__it;			__ti02 = __ti01+__rt;\
	    __tr01 = __tr01+__it;			__ti01 = __ti01-__rt;\
\
/*...Block 2: Swap inputs 3,4,5 -> 8,7,6: */\
	    __tr03 = __Ar08;				__ti03 = __Ai08;\
	    __tr04 = __Ar07;				__ti04 = __Ai07;\
	    __rt   = __Ar06;				__it   = __Ai06;\
	    __tr05 = __tr04-__rt;			__ti05 = __ti04-__it;\
	    __tr04 = __tr04+__rt;			__ti04 = __ti04+__it;\
	    __tr03 = __tr03+__tr04;			__ti03 = __ti03+__ti04;\
	    __tr04 = __tr03+c3m1*__tr04;	__ti04 = __ti03+c3m1*__ti04;\
	    __rt   = s*__tr05;				__it   = s*__ti05;\
	    __tr05 = __tr04-__it;			__ti05 = __ti04+__rt;\
	    __tr04 = __tr04+__it;			__ti04 = __ti04-__rt;\
\
/*...Block 3: Swap inputs 6,7,8 -> 13,12,14: */\
	    __tr06 = __Ar13;				__ti06 = __Ai13;\
	    __tr07 = __Ar12;				__ti07 = __Ai12;\
	    __rt   = __Ar14;				__it   = __Ai14;\
	    __tr08 = __tr07-__rt;			__ti08 = __ti07-__it;\
	    __tr07 = __tr07+__rt;			__ti07 = __ti07+__it;\
	    __tr06 = __tr06+__tr07;			__ti06 = __ti06+__ti07;\
	    __tr07 = __tr06+c3m1*__tr07;	__ti07 = __ti06+c3m1*__ti07;\
	    __rt   = s*__tr08;				__it   = s*__ti08;\
	    __tr08 = __tr07-__it;			__ti08 = __ti07+__rt;\
	    __tr07 = __tr07+__it;			__ti07 = __ti07-__rt;\
\
/*...Block 4: Swap inputs 9,10,11 -> 4,3,5: */\
	    __tr09 = __Ar04;				__ti09 = __Ai04;\
	    __tr10 = __Ar03;				__ti10 = __Ai03;\
	    __rt   = __Ar05;				__it   = __Ai05;\
	    __tr11 = __tr10-__rt;			__ti11 = __ti10-__it;\
	    __tr10 = __tr10+__rt;			__ti10 = __ti10+__it;\
	    __tr09 = __tr09+__tr10;			__ti09 = __ti09+__ti10;\
	    __tr10 = __tr09+c3m1*__tr10;	__ti10 = __ti09+c3m1*__ti10;\
	    __rt   = s*__tr11;				__it   = s*__ti11;\
	    __tr11 = __tr10-__it;			__ti11 = __ti10+__rt;\
	    __tr10 = __tr10+__it;			__ti10 = __ti10-__rt;\
\
/*...Block 2: Swap inputs 12,13,14 -> 9,11,10: */\
	    __tr12 = __Ar09;				__ti12 = __Ai09;\
	    __tr13 = __Ar11;				__ti13 = __Ai11;\
	    __rt   = __Ar10;				__it   = __Ai10;\
	    __tr14 = __tr13-__rt;			__ti14 = __ti13-__it;\
	    __tr13 = __tr13+__rt;			__ti13 = __ti13+__it;\
	    __tr12 = __tr12+__tr13;			__ti12 = __ti12+__ti13;\
	    __tr13 = __tr12+c3m1*__tr13;	__ti13 = __ti12+c3m1*__ti13;\
	    __rt   = s*__tr14;				__it   = s*__ti14;\
	    __tr14 = __tr13-__it;			__ti14 = __ti13+__rt;\
	    __tr13 = __tr13+__it;			__ti13 = __ti13-__rt;\
\
/*...and now do three radix-5 transforms.\
	The required output permutation is [0,5,10,9,14,4,3,8,13,12,2,7,6,11,1]. */\
/*...Block 1: output permutation is 0,9,3,12,6 */\
	    __rt   = __tr12;				__it   = __ti12;\
	    __tr12 = __tr03-__rt;			__ti12 = __ti03-__it;\
	    __tr03 = __tr03+__rt;			__ti03 = __ti03+__it;\
	    __rt   = __tr09;				__it   = __ti09;\
	    __tr09 = __tr06-__rt;			__ti09 = __ti06-__it;\
	    __tr06 = __tr06+__rt;			__ti06 = __ti06+__it;\
\
	    __rt   = __tr03+__tr06;			__it   = __ti03+__ti06;\
	    __tr00 = __tr00+__rt;			__ti00 = __ti00+__it;\
	    __rt   = __tr00+cn1*__rt;		__it   = __ti00+cn1*__it;\
	    __tr06 = cn2*(__tr03-__tr06);	__ti06 = cn2*(__ti03-__ti06);\
	    __tr03 = __rt+__tr06;			__ti03 = __it+__ti06;\
	    __tr06 = __rt-__tr06;			__ti06 = __it-__ti06;\
	    __rt   = ss3*(__tr09-__tr12);	__it   = ss3*(__ti09-__ti12);\
	    __tr09 = __rt-sn1*__tr09;		__ti09 = __it-sn1*__ti09;\
	    __tr12 = __rt+sn2*__tr12;		__ti12 = __it+sn2*__ti12;\
\
	    __Br00 = __tr00;				__Bi00 = __ti00;\
	    __Br09 = __tr03-__ti09;			__Bi09 = __ti03+__tr09;\
	    __Br03 = __tr06-__ti12;			__Bi03 = __ti06+__tr12;\
	    __Br12 = __tr06+__ti12;			__Bi12 = __ti06-__tr12;\
	    __Br06 = __tr03+__ti09;			__Bi06 = __ti03-__tr09;\
\
/*...Block 2: output permutation is 5,14,8,2,11 */\
	    __rt   = __tr13;				__it   = __ti13;\
	    __tr13 = __tr04-__rt;			__ti13 = __ti04-__it;\
	    __tr04 = __tr04+__rt;			__ti04 = __ti04+__it;\
	    __rt   = __tr10;				__it   = __ti10;\
	    __tr10 = __tr07-__rt;			__ti10 = __ti07-__it;\
	    __tr07 = __tr07+__rt;			__ti07 = __ti07+__it;\
\
	    __rt   = __tr04+__tr07;			__it   = __ti04+__ti07;\
	    __tr01 = __tr01+__rt;			__ti01 = __ti01+__it;\
	    __rt   = __tr01+cn1*__rt;		__it   = __ti01+cn1*__it;\
	    __tr07 = cn2*(__tr04-__tr07);	__ti07 = cn2*(__ti04-__ti07);\
	    __tr04 = __rt+__tr07;			__ti04 = __it+__ti07;\
	    __tr07 = __rt-__tr07;			__ti07 = __it-__ti07;\
	    __rt   = ss3*(__tr10-__tr13);	__it   = ss3*(__ti10-__ti13);\
	    __tr10 = __rt-sn1*__tr10;		__ti10 = __it-sn1*__ti10;\
	    __tr13 = __rt+sn2*__tr13;		__ti13 = __it+sn2*__ti13;\
\
	    __Br05 = __tr01;				__Bi05 = __ti01;\
	    __Br14 = __tr04-__ti10;			__Bi14 = __ti04+__tr10;\
	    __Br08 = __tr07-__ti13;			__Bi08 = __ti07+__tr13;\
	    __Br02 = __tr07+__ti13;			__Bi02 = __ti07-__tr13;\
	    __Br11 = __tr04+__ti10;			__Bi11 = __ti04-__tr10;\
\
/*...Block 3: output permutation is 10,4,13,7,1 */\
	    __rt   = __tr14;				__it   = __ti14;\
	    __tr14 = __tr05-__rt;			__ti14 = __ti05-__it;\
	    __tr05 = __tr05+__rt;			__ti05 = __ti05+__it;\
	    __rt   = __tr11;				__it   = __ti11;\
	    __tr11 = __tr08-__rt;			__ti11 = __ti08-__it;\
	    __tr08 = __tr08+__rt;			__ti08 = __ti08+__it;\
\
	    __rt   = __tr05+__tr08;			__it   = __ti05+__ti08;\
	    __tr02 = __tr02+__rt;			__ti02 = __ti02+__it;\
	    __rt   = __tr02+cn1*__rt;		__it   = __ti02+cn1*__it;\
	    __tr08 = cn2*(__tr05-__tr08);	__ti08 = cn2*(__ti05-__ti08);\
	    __tr05 = __rt+__tr08;			__ti05 = __it+__ti08;\
	    __tr08 = __rt-__tr08;			__ti08 = __it-__ti08;\
	    __rt   = ss3*(__tr11-__tr14);	__it   = ss3*(__ti11-__ti14);\
	    __tr11 = __rt-sn1*__tr11;		__ti11 = __it-sn1*__ti11;\
	    __tr14 = __rt+sn2*__tr14;		__ti14 = __it+sn2*__ti14;\
\
	    __Br10 = __tr02;				__Bi10 = __ti02;\
	    __Br04 = __tr05-__ti11;			__Bi04 = __ti05+__tr11;\
	    __Br13 = __tr08-__ti14;			__Bi13 = __ti08+__tr14;\
	    __Br07 = __tr08+__ti14;			__Bi07 = __ti08-__tr14;\
	    __Br01 = __tr05+__ti11;			__Bi01 = __ti05-__tr11;\
}

/* Twiddleless radix-16 DFTs: */
/* Totals: 144 ADD, 24 MUL. With twiddles adds 30 CMUL @[4 MUL, 2 ADD] each and thus needs 174 ADD, 84 MUL. */
#define RADIX_16_DIF(\
	__A0r,__A0i,__A1r,__A1i,__A2r,__A2i,__A3r,__A3i,__A4r,__A4i,__A5r,__A5i,__A6r,__A6i,__A7r,__A7i,__A8r,__A8i,__A9r,__A9i,__AAr,__AAi,__ABr,__ABi,__ACr,__ACi,__ADr,__ADi,__AEr,__AEi,__AFr,__AFi,\
	__B0r,__B0i,__B1r,__B1i,__B2r,__B2i,__B3r,__B3i,__B4r,__B4i,__B5r,__B5i,__B6r,__B6i,__B7r,__B7i,__B8r,__B8i,__B9r,__B9i,__BAr,__BAi,__BBr,__BBi,__BCr,__BCi,__BDr,__BDi,__BEr,__BEi,__BFr,__BFi,\
	__c,__s)\
{\
	double t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32;\
	double _r,_i;\
	/*       gather the needed data (16 64-bit complex, i.e. 32 64-bit reals) and do the first set of four length-4 transforms.	*/\
	/*...Block 1: */\
	t3 =__A0r -__A8r;	t1 =__A0r +__A8r;\
	t4 =__A0i -__A8i;	t2 =__A0i +__A8i;\
\
	t7 =__A4r -__ACr;	t5 =__A4r +__ACr;\
	t8 =__A4i -__ACi;	t6 =__A4i +__ACi;\
\
	_r =t5;		t5 =t1 -_r;	t1 =t1 +_r;\
	_i =t6;		t6 =t2 -_i;	t2 =t2 +_i;\
\
	_r =t7;		t7 =t3 +t8;	t3 =t3 -t8;\
				t8 =t4 -_r;	t4 =t4 +_r;\
\
	/*...Block 2: */\
	t11=__A2r -__AAr;	t9 =__A2r +__AAr;\
	t12=__A2i -__AAi;	t10=__A2i +__AAi;\
\
	t15=__A6r -__AEr;	t13=__A6r +__AEr;\
	t16=__A6i -__AEi;	t14=__A6i +__AEi;\
\
	_r =t13;	t13=t9 -_r;		t9 =t9 +_r;\
	_i =t14;	t14=t10-_i;		t10=t10+_i;\
\
	_r =t15;	t15=t11+t16;	t11=t11-t16;\
				t16=t12-_r;		t12=t12+_r;\
\
	/*...Block 3: */\
	t19=__A1r -__A9r;	t17=__A1r +__A9r;\
	t20=__A1i -__A9i;	t18=__A1i +__A9i;\
\
	t23=__A5r -__ADr;	t21=__A5r +__ADr;\
	t24=__A5i -__ADi;	t22=__A5i +__ADi;\
\
	_r =t21;	t21=t17-_r;		t17=t17+_r;\
	_i =t22;	t22=t18-_i;		t18=t18+_i;\
\
	_r =t23;	t23=t19+t24;	t19=t19-t24;\
				t24=t20-_r;		t20=t20+_r;\
\
	/*...Block 4: */\
	t27=__A3r -__ABr;	t25=__A3r +__ABr;\
	t28=__A3i -__ABi;	t26=__A3i +__ABi;\
\
	t31=__A7r -__AFr;	t29=__A7r +__AFr;\
	t32=__A7i -__AFi;	t30=__A7i +__AFi;\
\
	_r =t29;	t29=t25-_r;		t25=t25+_r;\
	_i =t30;	t30=t26-_i;		t26=t26+_i;\
\
	_r =t31;	t31=t27+t32;	t27=t27-t32;\
				t32=t28-_r;		t28=t28+_r;\
\
	/*...and now do four more radix-4 transforms, including the internal twiddle factors: */\
\
	/*...Block 1: t1,9,17,25	*/\
	_r =t9 ;	t9 =t1 -_r;	t1 =t1 +_r;\
	_i =t10;	t10=t2 -_i;	t2 =t2 +_i;\
\
	_r =t25;	t25=t17-_r;	t17=t17+_r;\
	_i =t26;	t26=t18-_i;	t18=t18+_i;\
\
	__B0r =t1+t17;			__B0i =t2+t18;\
	__B1r =t1-t17;			__B1i =t2-t18;\
\
	__B2r =t9 -t26;			__B2i =t10+t25;	/* mpy by E^4=i is inlined here...	*/\
	__B3r =t9 +t26;			__B3i =t10-t25;\
\
	/*...Block 3: t5,13,21,29	*/\
	_r =t13;t13=t5 +t14;	t5 =t5 -t14;	/* twiddle mpy by E^4 = I	*/\
			t14=t6 -_r;		t6 =t6 +_r;\
\
	_r =(t21-t22)*ISRT2;	t22=(t21+t22)*ISRT2;	t21=_r;	/* twiddle mpy by E^2; could also do mpy by ISRT2 directly on t21,22,29,30	*/\
	_r =(t30+t29)*ISRT2;	_i =(t30-t29)*ISRT2;	/* twiddle mpy by -E^6 is here...	*/\
	t29=t21+_r;				t21=t21-_r;				/* ...and get E^6=(i-1)/sq_r by flipping signs here.	*/\
	t30=t22+_i;				t22=t22-_i;\
\
	__B4r =t5+t21;			__B4i =t6+t22;\
	__B5r =t5-t21;			__B5i =t6-t22;\
\
	__B6r =t13-t30;			__B6i =t14+t29;			/* mpy by E^4=i is inlined here...	*/\
	__B7r =t13+t30;			__B7i =t14-t29;\
\
	/*...Block 2: t3,11,19,27	*/\
	_r =(t11-t12)*ISRT2;	_i =(t11+t12)*ISRT2;	/* twiddle mpy by E^2; could also do mpy by ISRT2 directly on t11,12	*/\
	t11=t3 -_r;				t3 =t3 +_r;\
	t12=t4 -_i;				t4 =t4 +_i;\
\
	_r =t19*__c - t20*__s;	t20=t20*__c + t19*__s;	t19=_r;	/* twiddle mpy by E^1	*/\
	_r =t27*__s - t28*__c;	_i =t28*__s + t27*__c;	/* twiddle mpy by E^3	*/\
	t27=t19-_r;				t19=t19+_r;\
	t28=t20-_i;				t20=t20+_i;\
\
	__B8r =t3+t19;			__B8i =t4+t20;\
	__B9r =t3-t19;			__B9i =t4-t20;\
\
	__BAr =t11-t28;			__BAi =t12+t27;			/* mpy by E^4=i is inlined here...	*/\
	__BBr =t11+t28;			__BBi =t12-t27;\
\
	/*...Block 4: t7,15,23,31	*/\
	_r =(t16+t15)*ISRT2;	_i =(t16-t15)*ISRT2;	/* twiddle mpy by -E^6 is here...	*/\
	t15=t7 +_r;				t7 =t7 -_r;			/* ...and get E^6=(i-1)/sqrt2 by flipping signs here.	*/\
	t16=t8 +_i;				t8 =t8 -_i;\
\
	_r =t23*__s - t24*__c;	t24=t24*__s + t23*__c;	t23=_r;	/* twiddle mpy by E^3	*/\
	_r =t31*__c - t32*__s;	_i =t32*__c + t31*__s;	/* twiddle mpy by E^1 = -E^9...	*/\
	t31=t23+_r;				t23=t23-_r;			/* ...and get E^9 by flipping signs here.	*/\
	t32=t24+_i;				t24=t24-_i;			/* Note: t23+_r = t23*(s+1)	*/\
\
	__BCr =t7+t23;			__BCi =t8+t24;\
	__BDr =t7-t23;			__BDi =t8-t24;\
\
	__BEr =t15-t32;			__BEi =t16+t31;			/* mpy by E^4=i is inlined here...	*/\
	__BFr =t15+t32;			__BFi =t16-t31;\
}

/* Totals: 144 ADD, 24 MUL	*/
#define RADIX_16_DIT(\
	__A0r,__A0i,__A1r,__A1i,__A2r,__A2i,__A3r,__A3i,__A4r,__A4i,__A5r,__A5i,__A6r,__A6i,__A7r,__A7i,__A8r,__A8i,__A9r,__A9i,__AAr,__AAi,__ABr,__ABi,__ACr,__ACi,__ADr,__ADi,__AEr,__AEi,__AFr,__AFi,\
	__B0r,__B0i,__B1r,__B1i,__B2r,__B2i,__B3r,__B3i,__B4r,__B4i,__B5r,__B5i,__B6r,__B6i,__B7r,__B7i,__B8r,__B8i,__B9r,__B9i,__BAr,__BAi,__BBr,__BBi,__BCr,__BCi,__BDr,__BDi,__BEr,__BEi,__BFr,__BFi,\
	__c,__s)\
{\
	double t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32;\
	double _r,_i;\
	/*...Block 1: */\
	t3 =__A0r -__A1r;	t1 =__A0r +__A1r;\
	t4 =__A0i -__A1i;	t2 =__A0i +__A1i;\
\
	t7 =__A2r -__A3r;	t5 =__A2r +__A3r;\
	t8 =__A2i -__A3i;	t6 =__A2i +__A3i;\
\
	_r =t5;		t5 =t1 -_r;	t1 =t1 +_r;\
	_i =t6;		t6 =t2 -_i;	t2 =t2 +_i;\
\
	_r =t7;		t7 =t3 -t8;	t3 =t3 +t8;\
				t8 =t4 +_r;	t4 =t4 -_r;\
\
	/*...Block 2: */\
	t11=__A4r -__A5r;	t9 =__A4r +__A5r;\
	t12=__A4i -__A5i;	t10=__A4i +__A5i;\
\
	t15=__A6r -__A7r;	t13=__A6r +__A7r;\
	t16=__A6i -__A7i;	t14=__A6i +__A7i;\
\
	_r =t13;	t13=t9 -_r;		t9 =t9 +_r;\
	_i =t14;	t14=t10-_i;		t10=t10+_i;\
\
	_r =t15;	t15=t11-t16;	t11=t11+t16;\
				t16=t12+_r;		t12=t12-_r;\
\
	/*...Block 3: */\
	t19=__A8r -__A9r;	t17=__A8r +__A9r;\
	t20=__A8i -__A9i;	t18=__A8i +__A9i;\
\
	t23=__AAr -__ABr;	t21=__AAr +__ABr;\
	t24=__AAi -__ABi;	t22=__AAi +__ABi;\
\
	_r =t21;	t21=t17-_r;		t17=t17+_r;\
	_i =t22;	t22=t18-_i;		t18=t18+_i;\
\
	_r =t23;	t23=t19-t24;	t19=t19+t24;\
				t24=t20+_r;		t20=t20-_r;\
\
	/*...Block 4: */\
	t27=__ACr -__ADr;	t25=__ACr +__ADr; \
	t28=__ACi -__ADi;	t26=__ACi +__ADi;\
\
	t31=__AEr -__AFr;	t29=__AEr +__AFr;\
	t32=__AEi -__AFi;	t30=__AEi +__AFi;\
\
	_r =t29;	t29=t25-_r;		t25=t25+_r;\
	_i =t30;	t30=t26-_i;		t26=t26+_i;\
\
	_r =t31;	t31=t27-t32;	t27=t27+t32;\
				t32=t28+_r;		t28=t28-_r;\
\
	/*...and now do four more radix-4 transforms, including the internal twiddle factors: */\
\
	/*...Block 1: t1,9,17,25	*/\
	_r =t9 ;	t9 =t1 -_r;	t1 =t1 +_r;\
	_i =t10;	t10=t2 -_i;	t2 =t2 +_i;\
\
	_r =t25;	t25=t17-_r;	t17=t17+_r;\
	_i =t26;	t26=t18-_i;	t18=t18+_i;\
\
	__B0r=t1+t17;			__B0i =t2+t18;\
	__B8r=t1-t17;			__B8i =t2-t18;\
\
	__B4r=t9 +t26;			__B4i =t10-t25;	/* mpy by E^-4 = -I is inlined here...	*/\
	__BCr=t9 -t26;			__BCi =t10+t25;\
\
	/*...Block 3: t5,13,21,29	*/\
	_r =t13;t13=t5 -t14;	t5 =t5 +t14;	/* twiddle mpy by E^4 =-I	*/\
			t14=t6 +_r;		t6 =t6 -_r;\
\
	_r =(t22+t21)*ISRT2;	t22=(t22-t21)*ISRT2;	t21=_r;	/* twiddle mpy by E^-2	*/\
	_r =(t29-t30)*ISRT2;	_i =(t29+t30)*ISRT2;	/* twiddle mpy by E^2 = -E^-6 is here...	*/\
	t29=t21+_r;				t21=t21-_r;				/* ...and get E^-6 by flipping signs here.	*/\
	t30=t22+_i;				t22=t22-_i;\
\
	__B2r=t5+t21;			__B2i =t6+t22;\
	__BAr=t5-t21;			__BAi =t6-t22;\
\
	__B6r=t13+t30;			__B6i =t14-t29;			/* mpy by E^-4 =-I is inlined here...	*/\
	__BEr=t13-t30;			__BEi =t14+t29;\
\
	/*...Block 2: t3,11,19,27	*/\
	_r =(t12+t11)*ISRT2;	_i =(t12-t11)*ISRT2;	/* twiddle mpy by E^-2	*/\
	t11=t3 -_r;				t3 =t3 +_r;\
	t12=t4 -_i;				t4 =t4 +_i;\
\
	_r =t19*__c + t20*__s;	t20=t20*__c - t19*__s;	t19=_r;	/* twiddle mpy by E^-1	*/\
	_r =t27*__s + t28*__c;	_i =t28*__s - t27*__c;	/* twiddle mpy by E^-3	*/\
	t27=t19-_r;				t19=t19+_r;\
	t28=t20-_i;				t20=t20+_i;\
\
	__B1r=t3+t19;			__B1i =t4+t20;\
	__B9r=t3-t19;			__B9i =t4-t20;\
\
	__B5r=t11+t28;			__B5i =t12-t27;			/* mpy by E^-4 =-I is inlined here...	*/\
	__BDr=t11-t28;			__BDi =t12+t27;\
\
	/*...Block 4: t7,15,23,31	*/\
	_r =(t15-t16)*ISRT2;	_i =(t15+t16)*ISRT2;	/* twiddle mpy by E^2 = -E^-6 is here...	*/\
	t15=t7 +_r;				t7 =t7 -_r;			/* ...and get E^6=(i-1)/sqrt2 by flipping signs here.	*/\
	t16=t8 +_i;				t8 =t8 -_i;\
\
	_r =t23*__s + t24*__c;	t24=t24*__s - t23*__c;	t23=_r;	/* twiddle mpy by E^-3	*/\
	_r =t31*__c + t32*__s;	_i =t32*__c - t31*__s;	/* twiddle mpy by E^-1 = -E^-9...	*/\
	t31=t23+_r;				t23=t23-_r;			/* ...and get E^9 by flipping signs here.	*/\
	t32=t24+_i;				t24=t24-_i;\
\
	__B3r=t7+t23;			__B3i =t8+t24;\
	__BBr=t7-t23;			__BBi =t8-t24;\
\
	__B7r=t15+t32;			__B7i =t16-t31;			/* mpy by E^-4 = -I is inlined here...	*/\
	__BFr=t15-t32;			__BFi =t16+t31;\
}

/* With-Twiddle radix-16 DFTs: */
/* Totals: 174 ADD, 84 MUL. */
#define RADIX_16_DIF_TWIDDLE(\
	__A0r,__A0i,__A1r,__A1i,__A2r,__A2i,__A3r,__A3i,__A4r,__A4i,__A5r,__A5i,__A6r,__A6i,__A7r,__A7i,__A8r,__A8i,__A9r,__A9i,__AAr,__AAi,__ABr,__ABi,__ACr,__ACi,__ADr,__ADi,__AEr,__AEi,__AFr,__AFi,\
	__B0r,__B0i,__B1r,__B1i,__B2r,__B2i,__B3r,__B3i,__B4r,__B4i,__B5r,__B5i,__B6r,__B6i,__B7r,__B7i,__B8r,__B8i,__B9r,__B9i,__BAr,__BAi,__BBr,__BBi,__BCr,__BCi,__BDr,__BDi,__BEr,__BEi,__BFr,__BFi,\
	__c1 ,__s1 ,__c2 ,__s2 ,__c3 ,__s3 ,__c4 ,__s4 ,__c5 ,__s5 ,__c6 ,__s6 ,__c7 ,__s7 ,__c8 ,__s8 ,__c9 ,__s9 ,__cA ,__sA ,__cB ,__sB ,__cC ,__sC ,__cD ,__sD ,__cE ,__sE ,__cF ,__sF ,\
	__c,__s)\
{\
	double t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32;\
	double _r,_i;\
	/* Gather the needed data (16 64-bit complex, i.e. 32 64-bit reals) and do first set of four length-4 transforms: */\
	/*...Block 1: */\
	t1 =__A0r;						t2 =__A0i;\
	rt =__A8r*__c8 -__A8i*__s8 ;	it =__A8i*__c8 +__A8r*__s8 ;\
	t3 =t1 -rt;						t1 =t1 +rt;\
	t4 =t2 -it;						t2 =t2 +it;\
\
	t5 =__A4r*__c4 -__A4i*__s4 ;	t6 =__A4i*__c4 +__A4r*__s4 ;\
	rt =__ACr*__cC -__ACi*__sC ;	it =__ACi*__cC +__ACr*__sC ;\
	t7 =t5 -rt;							t5 =t5 +rt;\
	t8 =t6 -it;							t6 =t6 +it;\
\
	rt =t5;	t5 =t1 -rt;					t1 =t1 +rt;\
	it =t6;	t6 =t2 -it;					t2 =t2 +it;\
\
	rt =t7;	t7 =t3 +t8;					t3 =t3 -t8;\
			t8 =t4 -rt;					t4 =t4 +rt;\
\
	/*...Block 2: */\
	t9 =__A2r*__c2 -__A2i*__s2 ;	t10=__A2i*__c2 +__A2r*__s2;\
	rt =__AAr*__cA -__AAi*__sA ;	it =__AAi*__cA +__AAr*__sA ;\
	t11=t9 -rt;						t9 =t9 +rt;\
	t12=t10-it;						t10=t10+it;\
\
	t13=__A6r*__c6 -__A6i*__s6 ;	t14=__A6i*__c6 +__A6r*__s6;\
	rt =__AEr*__cE -__AEi*__sE ;	it =__AEi*__cE +__AEr*__sE ;\
	t15=t13-rt;							t13=t13+rt;\
	t16=t14-it;							t14=t14+it;\
\
	rt =t13;	t13=t9 -rt;				t9 =t9 +rt;\
	it =t14;	t14=t10-it;				t10=t10+it;\
\
	rt =t15;	t15=t11+t16;			t11=t11-t16;\
				t16=t12-rt;				t12=t12+rt;\
\
	/*...Block 3: */\
	t17=__A1r*__c1 -__A1i*__s1 ;	t18=__A1i*__c1 +__A1r*__s1;\
	rt =__A9r*__c9 -__A9i*__s9 ;	it =__A9i*__c9 +__A9r*__s9;\
	t19=t17-rt;						t17=t17+rt;\
	t20=t18-it;						t18=t18+it;\
\
	t21=__A5r*__c5 -__A5i*__s5 ;	t22=__A5i*__c5 +__A5r*__s5;\
	rt =__ADr*__cD -__ADi*__sD ;	it =__ADi*__cD +__ADr*__sD ;\
	t23=t21-rt;							t21=t21+rt;\
	t24=t22-it;							t22=t22+it;\
\
	rt =t21;	t21=t17-rt;				t17=t17+rt;\
	it =t22;	t22=t18-it;				t18=t18+it;\
\
	rt =t23;	t23=t19+t24;			t19=t19-t24;\
				t24=t20-rt;				t20=t20+rt;\
\
	/*...Block 4: */\
	t25=__A3r*__c3 -__A3i*__s3 ;	t26=__A3i*__c3 +__A3r*__s3;\
	rt =__ABr*__cB -__ABi*__sB ;	it =__ABi*__cB +__ABr*__sB ;\
	t27=t25-rt;						t25=t25+rt;\
	t28=t26-it;						t26=t26+it;\
\
	t29=__A7r*__c7 -__A7i*__s7 ;	t30=__A7i*__c7 +__A7r*__s7;\
	rt =__AFr*__cF -__AFi*__sF ;	it =__AFi*__cF +__AFr*__sF ;\
	t31=t29-rt;							t29=t29+rt;\
	t32=t30-it;							t30=t30+it;\
\
	rt =t29;	t29=t25-rt;				t25=t25+rt;\
	it =t30;	t30=t26-it;				t26=t26+it;\
\
	rt =t31;	t31=t27+t32;			t27=t27-t32;\
				t32=t28-rt;				t28=t28+rt;\
\
/*...and now do four more radix-4 transforms, including the internal twiddle factors: */\
\
	/*...Block 1: t1,9,17,25 */\
	rt =t9;	t9 =t1 -rt;	t1 =t1 +rt;\
	it =t10;t10=t2 -it;	t2 =t2 +it;\
\
	rt =t25;t25=t17-rt;	t17=t17+rt;\
	it =t26;t26=t18-it;	t18=t18+it;\
\
	__A0r=t1+t17;	__A0i=t2+t18;\
	__A1r=t1-t17;	__A1i=t2-t18;\
\
	__A2r=t9 -t26;	__A2i=t10+t25;	/* mpy by E^4=i is inlined here... */\
	__A3r=t9 +t26;	__A3i=t10-t25;\
\
	/*...Block 3: t5,13,21,29 */\
	rt =t13;t13=t5 +t14;t5 =t5 -t14;		/* twiddle mpy by E^4 = I */\
			t14=t6 -rt;	t6 =t6 +rt;\
\
	rt =(t21-t22)*ISRT2;t22=(t21+t22)*ISRT2;	t21=rt;	/* twiddle mpy by E^2; could also do mpy by ISRT2 directly on t21,22,29,30 */\
	rt =(t30+t29)*ISRT2;it =(t30-t29)*ISRT2;		/* twiddle mpy by -E^6 is here... */\
	t29=t21+rt;			t21=t21-rt;			/* ...and get E^6=(i-1)/sqrt by flipping signs here. */\
	t30=t22+it;			t22=t22-it;\
\
	__A4r=t5+t21;	__A4i=t6+t22;\
	__A5r=t5-t21;	__A5i=t6-t22;\
\
	__A6r=t13-t30;	__A6i=t14+t29;	/* mpy by E^4=i is inlined here... */\
	__A7r=t13+t30;	__A7i=t14-t29;\
\
	/*...Block 2: t3,11,19,27 */\
	rt =(t11-t12)*ISRT2;it =(t11+t12)*ISRT2;		/* twiddle mpy by E^2; could also do mpy by ISRT2 directly on t11,12 */\
	t11=t3 -rt;			t3 =t3 +rt;\
	t12=t4 -it;			t4 =t4 +it;\
\
	rt =t19*__c - t20*__s;	t20=t20*__c + t19*__s;	t19=rt;	/* twiddle mpy by E^1 */\
	rt =t27*__s - t28*__c;	it =t28*__s + t27*__c;		/* twiddle mpy by E^3 */\
	t27=t19-rt;			t19=t19+rt;\
	t28=t20-it;			t20=t20+it;\
\
	__A8r=t3+t19;	__A8i=t4+t20;\
	__A9r=t3-t19;	__A9i=t4-t20;\
\
	__AAr=t11-t28;	__AAi=t12+t27;	/* mpy by E^4=i is inlined here... */\
	__ABr=t11+t28;	__ABi=t12-t27;\
\
	/*...Block 4: t7,15,23,31 */\
	rt =(t16+t15)*ISRT2;it =(t16-t15)*ISRT2;		/* twiddle mpy by -E^6 is here... */\
	t15=t7 +rt;			t7 =t7 -rt;			/* ...and get E^6=(i-1)/sqrt by flipping signs here. */\
	t16=t8 +it;			t8 =t8 -it;\
\
	rt =t23*__s - t24*__c;	t24=t24*__s + t23*__c;	t23=rt;	/* twiddle mpy by E^3 */\
	rt =t31*__c - t32*__s;	it =t32*__c + t31*__s;		/* twiddle mpy by E^1 = -E^9... */\
	t31=t23+rt;			t23=t23-rt;			/* ...and get E^9 by flipping signs here. */\
	t32=t24+it;			t24=t24-it;			/* Note: t23+rt = t23*(s+1) */\
\
	__ACr=t7+t23;	__ACi=t8+t24;\
	__ADr=t7-t23;	__ADi=t8-t24;\
\
	__AEr=t15-t32;	__AEi=t16+t31;	/* mpy by E^4=i is inlined here... */\
	__AFr=t15+t32;	__AFi=t16-t31;\
}

#define RADIX_16_DIT_TWIDDLE(\
	__A0r,__A0i,__A1r,__A1i,__A2r,__A2i,__A3r,__A3i,__A4r,__A4i,__A5r,__A5i,__A6r,__A6i,__A7r,__A7i,__A8r,__A8i,__A9r,__A9i,__AAr,__AAi,__ABr,__ABi,__ACr,__ACi,__ADr,__ADi,__AEr,__AEi,__AFr,__AFi,\
	__B0r,__B0i,__B1r,__B1i,__B2r,__B2i,__B3r,__B3i,__B4r,__B4i,__B5r,__B5i,__B6r,__B6i,__B7r,__B7i,__B8r,__B8i,__B9r,__B9i,__BAr,__BAi,__BBr,__BBi,__BCr,__BCi,__BDr,__BDi,__BEr,__BEi,__BFr,__BFi,\
	__c1 ,__s1 ,__c2 ,__s2 ,__c3 ,__s3 ,__c4 ,__s4 ,__c5 ,__s5 ,__c6 ,__s6 ,__c7 ,__s7 ,__c8 ,__s8 ,__c9 ,__s9 ,__cA ,__sA ,__cB ,__sB ,__cC ,__sC ,__cD ,__sD ,__cE ,__sE ,__cF ,__sF ,\
	__c,__s)\
{\
	double t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32;\
\
	/* Gather the needed data and do first set of four length-4 transforms: */\
	/*...Block 1: */\
	t1 =__A0r;	t2 =__A0i;\
	rt =__A1r;	it =__A1i;\
	t3 =t1 -rt;		t1 =t1 +rt;\
	t4 =t2 -it;		t2 =t2 +it;\
\
	t5 =__A2r;	t6 =__A2i;\
	rt =__A3r;	it =__A3i;\
	t7 =t5 -rt;		t5 =t5 +rt;\
	t8 =t6 -it;		t6 =t6 +it;\
\
	rt =t5;	t5 =t1 -rt ;	t1 =t1 +rt;\
	it =t6;	t6 =t2 -it ;	t2 =t2 +it;\
\
	rt =t7;	t7 =t3 -t8 ;	t3 =t3 +t8;\
			t8 =t4 +rt ;	t4 =t4 -rt;\
\
	/*...Block 2: */\
	rt =__A4r;	t10=__A4i;	t9 =rt;\
	rt =__A5r;	it =__A5i;\
	t11=t9 -rt;		t9 =t9 +rt;\
	t12=t10-it;		t10=t10+it;\
\
	rt =__A6r;	t14=__A6i;	t13=rt;\
	rt =__A7r;	it =__A7i;\
	t15=t13-rt;		t13=t13+rt;\
	t16=t14-it;		t14=t14+it;\
\
	rt =t13;	t13=t9 -rt ;	t9 =t9 +rt;\
	it =t14;	t14=t10-it ;	t10=t10+it;\
\
	rt =t15;	t15=t11-t16;	t11=t11+t16;\
				t16=t12+rt ;	t12=t12-rt;\
\
	/*...Block 3: */\
	rt =__A8r;	t18=__A8i;	t17=rt;\
	rt =__A9r;	it =__A9i;\
	t19=t17-rt;		t17=t17+rt;\
	t20=t18-it;		t18=t18+it;\
\
	rt =__AAr;	t22=__AAi;	t21=rt;\
	rt =__ABr;	it =__ABi;\
	t23=t21-rt;		t21=t21+rt;\
	t24=t22-it;		t22=t22+it;\
\
	rt =t21;	t21=t17-rt ;	t17=t17+rt;\
	it =t22;	t22=t18-it ;	t18=t18+it;\
\
	rt =t23;	t23=t19-t24;	t19=t19+t24;\
				t24=t20+rt ;	t20=t20-rt;\
\
	/*...Block 4: */\
	rt =__ACr;	t26=__ACi;	t25=rt;\
	rt =__ADr;	it =__ADi;\
	t27=t25-rt;		t25=t25+rt;\
	t28=t26-it;		t26=t26+it;\
\
	rt =__AEr;	t30=__AEi;	t29=rt;\
	rt =__AFr;	it =__AFi;\
	t31=t29-rt;		t29=t29+rt;\
	t32=t30-it;		t30=t30+it;\
\
	rt =t29;	t29=t25-rt ;	t25=t25+rt;\
	it =t30;	t30=t26-it ;	t26=t26+it;\
\
	rt =t31;	t31=t27-t32;	t27=t27+t32;\
				t32=t28+rt ;	t28=t28-rt;\
\
/*...and now do four more radix-4 transforms, including the internal twiddle factors: */\
\
	/*...Block 1: t1,9,17,25 */\
	rt =t9;	t9 =t1 -rt;	t1 =t1 +rt;\
	it =t10;	t10=t2 -it;	t2 =t2 +it;\
\
	rt =t25;	t25=t17-rt;	t17=t17+rt;\
	it =t26;	t26=t18-it;	t18=t18+it;\
\
	__B0r = t1+t17;			__B0i = t2+t18;\
	t1	     =t1-t17;			t2		 =t2-t18;\
	__B8r = t1 *__c8 +t2 *__s8;	__B8i = t2 *__c8 -t1 *__s8;\
\
	rt	   =t9 +t26;		it		 =t10-t25;	/* mpy by E^-4 = -I is inlined here... */\
	t9	   =t9 -t26;		t10		=t10+t25;\
	__B4r = rt *__c4 +it *__s4;	__B4i = it *__c4 -rt *__s4;\
	__BCr = t9 *__cC +t10*__sC;	__BCi = t10*__cC -t9 *__sC;\
\
	/*...Block 3: t5,13,21,29 */\
	rt =t13;	t13=t5 -t14;	t5 =t5 +t14;		/* twiddle mpy by E^4 =-I */\
		t14=t6 +rt;	t6 =t6 -rt;\
\
	rt =(t22+t21)*ISRT2;t22=(t22-t21)*ISRT2;	t21=rt;	/* twiddle mpy by E^-2 */\
	rt =(t29-t30)*ISRT2;it =(t29+t30)*ISRT2;		/* twiddle mpy by E^2 = -E^-6 is here... */\
	t29=t21+rt;		t21=t21-rt;			/* ...and get E^-6 by flipping signs here. */\
	t30=t22+it;		t22=t22-it;\
\
	rt	   =t5 +t21;		it		 =t6 +t22;\
	t5	   =t5 -t21;		t6		 =t6 -t22;\
	__B2r = rt *__c2 +it *__s2;	__B2i = it *__c2 -rt *__s2;\
	__BAr = t5 *__cA +t6 *__sA;	__BAi = t6 *__cA -t5 *__sA;\
\
	rt	   =t13+t30;		it		 =t14-t29;	/* mpy by E^-4 = -I is inlined here... */\
	t13	  =t13-t30;		t14		=t14+t29;\
	__B6r = rt *__c6 +it *__s6;	__B6i = it *__c6 -rt *__s6;\
	__BEr = t13*__cE +t14*__sE;	__BEi = t14*__cE -t13*__sE;\
\
	/*...Block 2: t3,11,19,27 */\
	rt =(t12+t11)*ISRT2;it =(t12-t11)*ISRT2;		/* twiddle mpy by E^-2 */\
	t11=t3 -rt;		t3 =t3 +rt;\
	t12=t4 -it;		t4 =t4 +it;\
\
	rt =t19*__c + t20*__s;	t20=t20*__c - t19*__s;	t19=rt;	/* twiddle mpy by E^-1 */\
	rt =t27*__s + t28*__c;	it =t28*__s - t27*__c;		/* twiddle mpy by E^-3 */\
	t27=t19-rt;		t19=t19+rt;\
	t28=t20-it;		t20=t20+it;\
\
	rt	   =t3 +t19;		it		 =t4 +t20;\
	t3	   =t3 -t19;		t4		 =t4 -t20;\
	__B1r = rt *__c1 +it *__s1;	__B1i = it *__c1 -rt *__s1;\
	__B9r = t3 *__c9 +t4 *__s9;	__B9i = t4 *__c9 -t3 *__s9;\
\
	rt	   =t11+t28;		it		 =t12-t27;	/* mpy by E^-4 = -I is inlined here... */\
	t11	  =t11-t28;		t12		=t12+t27;\
	__B5r = rt *__c5 +it *__s5;	__B5i = it *__c5 -rt *__s5;\
	__BDr = t11*__cD +t12*__sD;	__BDi = t12*__cD -t11*__sD;\
\
	/*...Block 4: t7,15,23,31 */\
	rt =(t15-t16)*ISRT2;it =(t15+t16)*ISRT2;		/* twiddle mpy by E^2 = -E^-6 is here... */\
	t15=t7 +rt;		t7 =t7 -rt;			/* ...and get E^6=(i-1)/sqrt by flipping signs here. */\
	t16=t8 +it;		t8 =t8 -it;\
\
	rt =t23*__s + t24*__c;	t24=t24*__s - t23*__c;	t23=rt;	/* twiddle mpy by E^-3 */\
	rt =t31*__c + t32*__s;	it =t32*__c - t31*__s;		/* twiddle mpy by E^-1 = -E^-9... */\
	t31=t23+rt;		t23=t23-rt;			/* ...and get E^9 by flipping signs here. */\
	t32=t24+it;		t24=t24-it;\
\
	rt	   =t7 +t23;		it		 =t8 +t24;\
	t7	   =t7 -t23;		t8		 =t8 -t24;\
	__B3r = rt *__c3 +it *__s3;	__B3i = it *__c3 -rt *__s3;\
	__BBr = t7 *__cB +t8 *__sB;	__BBi = t8 *__cB -t7 *__sB;\
\
	rt	   =t15+t32;		it		 =t16-t31;	/* mpy by E^-4 = -I is inlined here... */\
	t15	  =t15-t32;		t16		=t16+t31;\
	__B7r = rt *__c7 +it *__s7;	__B7i = it *__c7 -rt *__s7;\
	__BFr = t15*__cF +t16*__sF;	__BFi = t16*__cF -t15*__sF;\
}

/*================================*/

/*** To test the FMA-based macros below:
1. Change name of macro-to-test to remove trailing number and thus ending just in ...FMA;
2. cd ~/Mlucas/SRC/OBJ_SCALAR
3. gcc -c -O0 -DUSE_SCALAR_DFT_MACRO ../radix16*p*c && gcc -o Mlucas *o -lm && ./Mlucas -fftlen 128 -iters 100 -radset 1
***/

/* radix-16 FFT algorithms optimized for multiply-add instruction.
Floating-point Arithmetic Opcounts: If use FMA for all arithmetic including paired add/sub:

Algorithm includes 15 complex multiplications by twiddle factors omega_j = cr_j + I*ci_j, j = 1-15,
which can be computed in advance, in either raw (sincos) or "tangentized" form, computed below:

Total: 174 Fused mul/add [Breakdown: 87 FMA, 2 FMS, 85 FNMA, 0 FNMS.]
*/
// Version #1 uses 4-operand FMA4 syntax, similar to what we might use for SIMD assembler on SSE5-supporting AMD CPUs:
#define RADIX_16_DIF_FMA4(\
	__A0r,__A0i,__A1r,__A1i,__A2r,__A2i,__A3r,__A3i,__A4r,__A4i,__A5r,__A5i,__A6r,__A6i,__A7r,__A7i,__A8r,__A8i,__A9r,__A9i,__AAr,__AAi,__ABr,__ABi,__ACr,__ACi,__ADr,__ADi,__AEr,__AEi,__AFr,__AFi,\
	__B0r,__B0i,__B1r,__B1i,__B2r,__B2i,__B3r,__B3i,__B4r,__B4i,__B5r,__B5i,__B6r,__B6i,__B7r,__B7i,__B8r,__B8i,__B9r,__B9i,__BAr,__BAi,__BBr,__BBi,__BCr,__BCi,__BDr,__BDi,__BEr,__BEi,__BFr,__BFi,\
	__c1 ,__s1 ,__c2 ,__s2 ,__c3 ,__s3 ,__c4 ,__s4 ,__c5 ,__s5 ,__c6 ,__s6 ,__c7 ,__s7 ,__c8 ,__s8 ,__c9 ,__s9 ,__cA ,__sA ,__cB ,__sB ,__cC ,__sC ,__cD ,__sD ,__cE ,__sE ,__cF ,__sF ,\
	__c,__s)\
{\
	double t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32;\
	double _c,_d, _a,_b;\
	/* The divides in the following sincos data will be precomputed in an optimized version: */\
	double __t1 = __s1/__c1;\
	double __t2 = __s2/__c2;\
	double __t3 = __s3/__c3;\
	double __t4 = __s4/__c4;\
	double __t5 = __s5/__c5;\
	double __t6 = __s6/__c6;\
	double __t7 = __s7/__c7;\
	double __t8 = __s8/__c8;\
	double __t9 = __s9/__c9;\
	double __tA = __sA/__cA;\
	double __tB = __sB/__cB;\
	double __tC = __sC/__cC;\
	double __tD = __sD/__cD;\
	double __tE = __sE/__cE;\
	double __tF = __sF/__cF;\
	double __sc = __s/__c;\
	double __c1_c = __c1*__c;\
	double __c1i2 = __c1*ISRT2;\
	double __c2i2 = __c2*ISRT2;\
	double __c31 = __c3/__c1;\
	double __c51 = __c5/__c1;\
	double __c62 = __c6/__c2;\
	double __c73 = __c7/__c3;\
	double __c91 = __c9/__c1;\
	double __cA2 = __cA/__c2;\
	double __cB3 = __cB/__c3;\
	double __cC4 = __cC/__c4;\
	double __cD5 = __cD/__c5;\
	double __cE6 = __cE/__c6;\
	double __cF7 = __cF/__c7;\
\
	/* Gather the needed data and do first set of four length-4 transforms: */\
	/*...Block 1: */\
	t1 =__A0r;						t2 =__A0i;\
	_a = FNMA4(__A8i,__t8,__A8r);	_b = FMA4(__A8r,__t8,__A8i);\
	t3 = FNMA4(_a,__c8,t1);			t1 = FMA4(_a,__c8,t1);\
	t4 = FNMA4(_b,__c8,t2);			t2 = FMA4(_b,__c8,t2);\
\
	t5 = FNMA4(__A4i,__t4,__A4r);	t6 = FMA4(__A4r,__t4,__A4i);\
	_a = FNMA4(__ACi,__tC,__ACr);	_b = FMA4(__ACr,__tC,__ACi);\
	t7 = FNMA4(_a,__cC4,t5);			t5 = FMA4(_a,__cC4,t5);\
	t8 = FNMA4(_b,__cC4,t6);			t6 = FMA4(_b,__cC4,t6);\
	/* t5-8 need remultiply by c4: */\
	_a =t5;	t5 =FNMA4(_a,__c4,t1);	t1 = FMA4(_a,__c4,t1);\
	_b =t6;	t6 =FNMA4(_b,__c4,t2);	t2 = FMA4(_b,__c4,t2);\
	_a =t7;	t7 = FMA4(t8,__c4,t3);	t3 =FNMA4(t8,__c4,t3);\
			t8 =FNMA4(_a,__c4,t4);	t4 = FMA4(_a,__c4,t4);\
\
	/*...Block 2: */\
	t9 = FNMA4(__A2i,__t2,__A2r);	t10= FMA4(__A2r,__t2,__A2i);\
	_a = FNMA4(__AAi,__tA,__AAr);	_b = FMA4(__AAr,__tA,__AAi);\
	t11= FNMA4(_a,__cA2,t9 );		t9 = FMA4(_a,__cA2,t9 );\
	t12= FNMA4(_b,__cA2,t10);		t10= FMA4(_b,__cA2,t10);\
	/* t9-12 need remultiply by c2: */\
\
	t13= FNMA4(__A6i,__t6,__A6r);	t14= FMA4(__A6r,__t6,__A6i);\
	_a = FNMA4(__AEi,__tE,__AEr);	_b = FMA4(__AEr,__tE,__AEi);\
	t15= FNMA4(_a,__cE6,t13);		t13= FMA4(_a,__cE6,t13);\
	t16= FNMA4(_b,__cE6,t14);		t14= FMA4(_b,__cE6,t14);\
	/* t13-16 need remultiply by c6: */\
\
	_a =t13;t13=FNMA4( _a,__c62,t9 );t9 = FMA4( _a,__c62,t9 );\
	_b =t14;t14=FNMA4( _b,__c62,t10);t10= FMA4( _b,__c62,t10);\
	_a =t15;t15= FMA4(t16,__c62,t11);t11=FNMA4(t16,__c62,t11);\
			t16=FNMA4( _a,__c62,t12);t12= FMA4( _a,__c62,t12);\
	/* t9-16 need remultiply by c2: */\
\
	/*...Block 3: */\
	t17= FNMA4(__A1i,__t1,__A1r);	t18= FMA4(__A1r,__t1,__A1i);\
	_a = FNMA4(__A9i,__t9,__A9r);	_b = FMA4(__A9r,__t9,__A9i);\
	t19= FNMA4(_a,__c91,t17);		t17= FMA4(_a,__c91,t17);\
	t20= FNMA4(_b,__c91,t18);		t18= FMA4(_b,__c91,t18);\
	/* t17-20 need remultiply by c1: */\
\
	t21= FNMA4(__A5i,__t5,__A5r);	t22= FMA4(__A5r,__t5,__A5i);\
	_a = FNMA4(__ADi,__tD,__ADr);	_b = FMA4(__ADr,__tD,__ADi);\
	t23= FNMA4(_a,__cD5,t21);		t21= FMA4(_a,__cD5,t21);\
	t24= FNMA4(_b,__cD5,t22);		t22= FMA4(_b,__cD5,t22);\
	/* t21-24 need remultiply by c5: */\
\
	_a =t21;t21=FNMA4( _a,__c51,t17);t17= FMA4( _a,__c51,t17);\
	_b =t22;t22=FNMA4( _b,__c51,t18);t18= FMA4( _b,__c51,t18);\
	_a =t23;t23= FMA4(t24,__c51,t19);t19=FNMA4(t24,__c51,t19);\
			t24=FNMA4( _a,__c51,t20);t20= FMA4( _a,__c51,t20);\
	/* t17-24 need remultiply by c1: */\
\
	/*...Block 4: */\
	t25= FNMA4(__A3i,__t3,__A3r);	t26= FMA4(__A3r,__t3,__A3i);\
	_a = FNMA4(__ABi,__tB,__ABr);	_b = FMA4(__ABr,__tB,__ABi);\
	t27= FNMA4(_a,__cB3,t25);		t25= FMA4(_a,__cB3,t25);\
	t28= FNMA4(_b,__cB3,t26);		t26= FMA4(_b,__cB3,t26);\
	/* t25-28 need remultiply by c3: */\
\
	t29= FNMA4(__A7i,__t7,__A7r);	t30= FMA4(__A7r,__t7,__A7i);\
	_a = FNMA4(__AFi,__tF,__AFr);	_b = FMA4(__AFr,__tF,__AFi);\
	t31= FNMA4(_a,__cF7,t29);		t29= FMA4(_a,__cF7,t29);\
	t32= FNMA4(_b,__cF7,t30);		t30= FMA4(_b,__cF7,t30);\
	/* t29-32 need remultiply by c7: */\
\
	_a =t29;t29=FNMA4( _a,__c73,t25);t25= FMA4( _a,__c73,t25);\
	_b =t30;t30=FNMA4( _b,__c73,t26);t26= FMA4( _b,__c73,t26);\
	_a =t31;t31= FMA4(t32,__c73,t27);t27=FNMA4(t32,__c73,t27);\
			t32=FNMA4( _a,__c73,t28);t28= FMA4( _a,__c73,t28);\
	/* t25-32 need remultiply by c3: */\
\
/*...and now do four more radix-4 transforms, including the internal twiddle factors: */\
\
	/*...Block 1: t1,9,17,25 */\
	_a =t9;	t9 =FNMA4(_a,__c2,t1);	t1 = FMA4(_a,__c2,t1);\
	_b =t10;t10=FNMA4(_b,__c2,t2);	t2 = FMA4(_b,__c2,t2);\
\
	_a =t25;t25=FNMA4(_a,__c31,t17);	t17= FMA4(_a,__c31,t17);\
	_b =t26;t26=FNMA4(_b,__c31,t18);	t18= FMA4(_b,__c31,t18);\
\
	__B0r= FMA4(t17,__c1,t1);		__B0i= FMA4(t18,__c1,t2);\
	__B1r=FNMA4(t17,__c1,t1);		__B1i=FNMA4(t18,__c1,t2);\
\
	__B2r=FNMA4(t26,__c1,t9);		__B2i= FMA4(t25,__c1,t10);	/* mpy by E^4=i is inlined here... */\
	__B3r= FMA4(t26,__c1,t9);		__B3i=FNMA4(t25,__c1,t10);\
\
	/*...Block 3: t5,13,21,29 */\
	_a =t13;t13= FMA4(t14,__c2,t5 );	t5 =FNMA4(t14,__c2,t5 );		/* twiddle mpy by E^4 = I */\
			t14=FNMA4( _a,__c2,t6 );	t6 = FMA4( _a,__c2,t6 );\
\
	_c =FNMA4(t22,1.0,t21);			_d = FMA4(t22,1.0,t21);		/* twiddle mpy by E^2; ISRT2 absorbed into various multipliers */\
	_a = FMA4(t29,1.0,t30);			_b =FNMA4(t29,1.0,t30);		/* twiddle mpy by -E^6 is here... */\
	t29= FMA4(_a,__c31,_c);			t21=FNMA4(_a,__c31,_c);		/* ...and get E^6=[i-1]/sqrt2 by flipping signs here. */\
	t30= FMA4(_b,__c31,_d);			t22=FNMA4(_b,__c31,_d);\
\
	__B4r= FMA4(t21,__c1i2,t5 );		__B4i= FMA4(t22,__c1i2,t6 );\
	__B5r=FNMA4(t21,__c1i2,t5 );		__B5i=FNMA4(t22,__c1i2,t6 );\
\
	__B6r=FNMA4(t30,__c1i2,t13);		__B6i= FMA4(t29,__c1i2,t14);	/* mpy by E^4=i is inlined here... */\
	__B7r= FMA4(t30,__c1i2,t13);		__B7i=FNMA4(t29,__c1i2,t14);\
\
	/*...Block 2: t3,11,19,27 */\
	_a =FNMA4(t12,1.0,t11);			_b = FMA4(t12,1.0,t11);		/* twiddle mpy by -E^6 is here... */\
	t11=FNMA4(_a,__c2i2,t3);			t3 = FMA4(_a,__c2i2,t3);		/* ...and get E^6=[i-1]/sqrt2 by flipping signs here. */\
	t12=FNMA4(_b,__c2i2,t4);			t4 = FMA4(_b,__c2i2,t4);\
\
	_c =FNMA4(t20,__sc,t19);			_d = FMA4(t19,__sc,t20);		/* twiddle mpy by E^1 */\
	_a = FMS4(t27,__sc,t28);			_b = FMA4(t28,__sc,t27);		/* twiddle mpy by E^3 */\
	t27=FNMA4(_a,__c31,_c);			t19= FMA4(_a,__c31,_c);\
	t28=FNMA4(_b,__c31,_d);			t20= FMA4(_b,__c31,_d);\
\
	__B8r= FMA4(t19,__c1_c,t3 );		__B8i= FMA4(t20,__c1_c,t4 );\
	__B9r=FNMA4(t19,__c1_c,t3 );		__B9i=FNMA4(t20,__c1_c,t4 );\
\
	__BAr=FNMA4(t28,__c1_c,t11);		__BAi= FMA4(t27,__c1_c,t12);	/* mpy by E^4=i is inlined here... */\
	__BBr= FMA4(t28,__c1_c,t11);		__BBi=FNMA4(t27,__c1_c,t12);\
\
	/*...Block 4: t7,15,23,31 */\
	_a = FMA4(t15,1.0,t16);			_b =FNMA4(t15,1.0,t16);		/* twiddle mpy by -E^6 is here... */\
	t15= FMA4(_a,__c2i2,t7);			t7 =FNMA4(_a,__c2i2,t7);		/* ...and get E^6=[i-1]/sqrt2 by flipping signs here. */\
	t16= FMA4(_b,__c2i2,t8);			t8 =FNMA4(_b,__c2i2,t8);\
\
	_c = FMS4(t23,__sc,t24);			_d = FMA4(t24,__sc,t23);		/* twiddle mpy by E^3 */\
	_a =FNMA4(t32,__sc,t31);			_b = FMA4(t31,__sc,t32);		/* twiddle mpy by E^1 =__c*[ -E^9... */\
	t31= FMA4(_a,__c31,_c);			t23=FNMA4(_a,__c31,_c);		/* ...and get E^9 by flipping signs here. */\
	t32= FMA4(_b,__c31,_d);			t24=FNMA4(_b,__c31,_d);		/* Note: t23+_a = t23*[s+1] */\
\
	__BCr= FMA4(t23,__c1_c,t7 );		__BCi= FMA4(t24,__c1_c,t8 );\
	__BDr=FNMA4(t23,__c1_c,t7 );		__BDi=FNMA4(t24,__c1_c,t8 );\
\
	__BEr=FNMA4(t32,__c1_c,t15);		__BEi= FMA4(t31,__c1_c,t16);	/* mpy by E^4=i is inlined here... */\
	__BFr= FMA4(t32,__c1_c,t15);		__BFi=FNMA4(t31,__c1_c,t16);\
}

// Version #2 uses a 3-operand FMA syntax in which the RIGHTmost of the 3 fused mul-add inputs -
// that is, the addend - is overwritten with the result:
#define RADIX_16_DIF_FMA(\
	__A0r,__A0i,__A1r,__A1i,__A2r,__A2i,__A3r,__A3i,__A4r,__A4i,__A5r,__A5i,__A6r,__A6i,__A7r,__A7i,__A8r,__A8i,__A9r,__A9i,__AAr,__AAi,__ABr,__ABi,__ACr,__ACi,__ADr,__ADi,__AEr,__AEi,__AFr,__AFi,\
	__B0r,__B0i,__B1r,__B1i,__B2r,__B2i,__B3r,__B3i,__B4r,__B4i,__B5r,__B5i,__B6r,__B6i,__B7r,__B7i,__B8r,__B8i,__B9r,__B9i,__BAr,__BAi,__BBr,__BBi,__BCr,__BCi,__BDr,__BDi,__BEr,__BEi,__BFr,__BFi,\
	__c1 ,__s1 ,__c2 ,__s2 ,__c3 ,__s3 ,__c4 ,__s4 ,__c5 ,__s5 ,__c6 ,__s6 ,__c7 ,__s7 ,__c8 ,__s8 ,__c9 ,__s9 ,__cA ,__sA ,__cB ,__sB ,__cC ,__sC ,__cD ,__sD ,__cE ,__sE ,__cF ,__sF ,\
	__c,__s)\
{\
	double t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31;\
	double _a,_b,_c,_d,_e,_f;\
	/* The divides in the following sincos data will be precomputed in an optimized version. */\
	/* In addition to the 30 derived multipliers above, we need the unmodified __c1,2,4,8, for a total of 34: */\
	double __sc = __s/__c;		/* In SIMD code, this will take the place of __s */\
	double __c1_c = __c1*__c;	/* In SIMD code, this will take the place of __c */\
	double __r1 = __s1/__c1;	/* In SIMD code,  these 15 tangents will take the place of __s1-F: */\
	double __r2 = __s2/__c2;\
	double __r3 = __s3/__c3;\
	double __r4 = __s4/__c4;\
	double __r5 = __s5/__c5;\
	double __r6 = __s6/__c6;\
	double __r7 = __s7/__c7;\
	double __r8 = __s8/__c8;\
	double __r9 = __s9/__c9;\
	double __rA = __sA/__cA;\
	double __rB = __sB/__cB;\
	double __rC = __sC/__cC;\
	double __rD = __sD/__cD;\
	double __rE = __sE/__cE;\
	double __rF = __sF/__cF;\
	double __c1i2 = __c1*ISRT2;	/* In SIMD code, this will take the place of __c0 */\
	double __c2i2 = __c2*ISRT2;	/* In SIMD code, this will take the place of __s0 */\
/*	double __c1;*/	/* In SIMD code,  these 15 cosine-derived terms will take the place of __c1-F: */\
/*	double __c2;*/\
	double __c31 = __c3/__c1;\
/*	double __c4;*/\
	double __c51 = __c5/__c1;\
	double __c62 = __c6/__c2;\
	double __c73 = __c7/__c3;\
/*	double __c8;*/\
	double __c91 = __c9/__c1;\
	double __cA2 = __cA/__c2;\
	double __cB3 = __cB/__c3;\
	double __cC4 = __cC/__c4;\
	double __cD5 = __cD/__c5;\
	double __cE6 = __cE/__c6;\
	double __cF7 = __cF/__c7;\
\
/* Gather the needed data and do first set of four length-4 transforms: */\
/* Remember that for the FMA231 macro (defined in float_intrin.h), the result overwrites the rightmost input! */\
	/*...Block 1: */\
	t04 =t06 =__A8r;			t05 =__A8i;\
	t00 =__A0r;					t01 =__A0i;/* load __c8,__r8 into pair of regs */\
	FNMA231(   t05,__r8,t04);	 FMA231(   t06,__r8,t05);/* Try to queue up as many FMAs as registers allow */\
	_a =__A4r;					_b =__A4i;/* load __cC4,__r4,__rC into trio of regs */\
	FNMA231(__A4i,__r4,_a);		 FMA231(__A4r,__r4,_b);\
	t06 =__ACr;					t07 =__ACi;\
	FNMA231(__ACi,__rC,t06);	 FMA231(__ACr,__rC,t07);\
	/* Swapped upper left and lower right FMA-pairs of octet below since _a/_b/t0/1 needed before _c/_d/t2/3: */\
	_c = _a;	t02 = t00;\
	 FMA231(t06,__cC4,_a);		 FMA231(t04,__c8,t00);\
	_d = _b;	t03 = t01;\
	 FMA231(t07,__cC4,_b);		 FMA231(t05,__c8,t01);\
	FNMA231(t06,__cC4,_c);		FNMA231(t04,__c8,t02);\
	FNMA231(t07,__cC4,_d);		FNMA231(t05,__c8,t03);\
	/* This snip then needs 13 registers: Data in t0-7, tmps t4,b,c,x,y,u,v, sincos datum __c4. */\
	/* Arrange these register-copies in same order as outputs of preceding 8 FMAs: */\
	/* load __c4 into reg */\
	t04 =t00; t05 =t01;\
	FNMA231(_a ,__c4 ,t04);		 FMA231(_a ,__c4 ,t00);\
	FNMA231(_b ,__c4 ,t05);		 FMA231(_b ,__c4 ,t01);\
	t06 =t02;	t07 =t03;\
	 FMA231(_d ,__c4 ,t06);		FNMA231(_d ,__c4 ,t02);\
	FNMA231(_c ,__c4 ,t07);		 FMA231(_c ,__c4 ,t03);\
	/* [In ASM version, write outputs into local store here] */\
\
	/*...Block 2: */\
	t08 =__A2r;					t09=__A2i;\
	FNMA231(__A2i,__r2,t08 );	 FMA231(__A2r,__r2,t09);\
	t12 =__AAr;					t13 =__AAi;\
	FNMA231(__AAi,__rA,t12 );	 FMA231(__AAr,__rA,t13 );\
	_a=__A6r;					_b=__A6i;\
	FNMA231(__A6i,__r6,_a);		 FMA231(__A6r,__r6,_b);\
	t14 =__AEr;					t15 =__AEi;\
	FNMA231(__AEi,__rE,t14 );	 FMA231(__AEr,__rE,t15 );\
	/* Swapped upper left and lower right FMA-pairs of octet below to hide latency: */\
	_c= _a;	t10= t08;\
	 FMA231(t14,__cE6,_a);		 FMA231(t12,__cA2,t08 );\
	_d= _b;	t11= t09;\
	 FMA231(t15,__cE6,_b);		 FMA231(t13,__cA2,t09);\
	FNMA231(t14,__cE6,_c);		FNMA231(t12,__cA2,t10);\
	FNMA231(t15,__cE6,_d);		FNMA231(t13,__cA2,t11);\
	/* t8-12 need remultiply by c2; _a-16 by c6: */\
	t12 =t08 ;	t13 =t09;\
	FNMA231(_a,__c62,t12);		 FMA231(_a,__c62,t08 );\
	FNMA231(_b,__c62,t13);		 FMA231(_b,__c62,t09);\
	t14 =t10;	t15 =t11;\
	 FMA231(_d,__c62,t14);		FNMA231(_d,__c62,t10);\
	FNMA231(_c,__c62,t15);		 FMA231(_c,__c62,t11);\
	/* t8-16 need remultiply by c2: */\
	/* [In ASM version, write outputs into local store here] */\
\
	/*...Block 3: */\
	t16=__A1r;	t17=__A1i;\
	FNMA231(__A1i,__r1,t16);	 FMA231(__A1r,__r1,t17);\
	t20 =__A9r;	t21 =__A9i;\
	FNMA231(__A9i,__r9,t20 );	 FMA231(__A9r,__r9,t21 );\
	_a=__A5r;	_b=__A5i;\
	FNMA231(__A5i,__r5,_a);		 FMA231(__A5r,__r5,_b);\
	t22 =__ADr;	t23 =__ADi;\
	FNMA231(__ADi,__rD,t22 );	 FMA231(__ADr,__rD,t23 );\
	/* Swapped upper left and lower right FMA-pairs of octet below to hide latency: */\
	_c= _a;	t18= t16;\
	 FMA231(t22,__cD5,_a);		 FMA231(t20,__c91,t16);\
	_d= _b;	t19= t17;\
	 FMA231(t23,__cD5,_b);		 FMA231(t21,__c91,t17);\
	FNMA231(t22,__cD5,_c);		FNMA231(t20,__c91,t18);\
	FNMA231(t23,__cD5,_d);		FNMA231(t21,__c91,t19);\
	/* t16-19 need remultiply by c1; t20-23 by c5: */\
	t20 =t16;	t21 =t17;\
	FNMA231(_a,__c51,t20);		 FMA231(_a,__c51,t16);\
	FNMA231(_b,__c51,t21);		 FMA231(_b,__c51,t17);\
	t22 =t18;	t23 =t19;\
	 FMA231(_d,__c51,t22);		FNMA231(_d,__c51,t18);\
	FNMA231(_c,__c51,t23);		 FMA231(_c,__c51,t19);\
	/* t16-23 need remultiply by c1: */\
	/* [In ASM version, write outputs into local store here] */\
\
	/*...Block 4: */\
	t24=__A3r;	t25=__A3i;\
	FNMA231(__A3i,__r3,t24);	 FMA231(__A3r,__r3,t25);\
	t28 =__ABr;	t29 =__ABi;\
	FNMA231(__ABi,__rB,t28 );	 FMA231(__ABr,__rB,t29 );\
	_a=__A7r;	_b=__A7i;\
	FNMA231(__A7i,__r7,_a);		 FMA231(__A7r,__r7,_b);\
	t30 =__AFr;	t31 =__AFi;\
	FNMA231(__AFi,__rF,t30 );	 FMA231(__AFr,__rF,t31 );\
	/* Swapped upper left and lower right FMA-pairs of octet below to hide latency: */\
	_c= _a;	t26= t24;\
	 FMA231(t30,__cF7,_a);		 FMA231(t28,__cB3,t24);\
	_d= _b;	t27= t25;\
	 FMA231(t31,__cF7,_b);		 FMA231(t29,__cB3,t25);\
	FNMA231(t30,__cF7,_c);		FNMA231(t28,__cB3,t26);\
	FNMA231(t31,__cF7,_d);		FNMA231(t29,__cB3,t27);\
	/* t24-27 need remultiply by c3; t28-31 by c7: */\
	t28 =t24;	t29 =t25;\
	FNMA231(_a,__c73,t28);		 FMA231(_a,__c73,t24);\
	FNMA231(_b,__c73,t29);		 FMA231(_b,__c73,t25);\
	t30 =t26;	t31 =t27;\
	 FMA231(_d,__c73,t30);		FNMA231(_d,__c73,t26);\
	FNMA231(_c,__c73,t31);		 FMA231(_c,__c73,t27);\
	/* t24-31 need remultiply by c3: */\
	/* [In ASM version, write outputs into local store here] */\
\
/*...and now do four more radix-4 transforms, including the internal twiddle factors: */\
\
	/*...Block 1: t0,8,16,24 ... Do the 4 Im-part FMAs first, because their results needed 1st below */\
	_a=t00; _c=t16;\
	 FMA231(t08,__c2 ,t00);		 FMA231(t24,__c31,t16);\
	_b=t01; _d=t17;\
	 FMA231(t09,__c2 ,t01);		 FMA231(t25,__c31,t17);\
	FNMA231(t08,__c2 ,_a );		FNMA231(t24,__c31,_c );\
	FNMA231(t09,__c2 ,_b );		FNMA231(t25,__c31,_d );\
	/* What FMA4 calls t8/9/24/25 are output in _a/_b/_c/_d in FMA3 version: */\
/******* need to overlap 2 blocks at a time to hide FMA latency - consider starting\
    Block 2 right here, using the just-freed t8/9/24/25 as temps for that *******/\
	_e = t00; _f = t01;\
	 FMA231(t16,__c1 ,t00);		 FMA231(t17,__c1 ,t01);\
	FNMA231(t16,__c1 ,_e );		FNMA231(t17,__c1 ,_f );\
	t08 = _a ; t09 = _b;\
	FNMA231(_d ,__c1 ,t08);		 FMA231(_c ,__c1 ,t09);	/* mpy by E^4=i is inlined here... */\
	 FMA231(_d ,__c1 ,_a );		FNMA231(_c ,__c1 ,_b );\
	/* Write outputs back to main array: */\
	__B0r= t00;		__B0i= t01;\
	__B1r= _e ;		__B1i= _f ;\
	__B2r= t08;		__B2i= t09;\
	__B3r= _a ;		__B3i= _b ;\
\
	/*...Block 3: t4,12,20,28 */\
	_c = t20; _d = t29;\
	FNMA231(t21,1.0,_c );		 FMA231(t21,1.0,t20);		/* twiddle mpy by E^2; ISRT2 absorbed into various multipliers */\
	 FMA231(t28,1.0,_d );		FNMA231(t28,1.0,t29);		/* twiddle mpy by -E^6 is here... */\
	/* Swapped upper right and lower left FMA-pairs of octet below to hide latency: */\
	_a = t04; _b = t05;\
	t21 = _c; _e = t20;\
	 FMA231(t13,__c2 ,_a );		 FMA231(_d ,__c31,t21);		/* twiddle mpy by E^4 = I */\
	FNMA231(t12,__c2 ,_b );		 FMA231(t29,__c31,t20);\
	FNMA231(t13,__c2 ,t04);		FNMA231(_d ,__c31,_c );		/* ...and get E^6=[i-1]/sqrt2 by flipping signs here. */\
	 FMA231(t12,__c2 ,t05);		FNMA231(t29,__c31,_e );\
	t12 = _a; t13 = _b;\
	FNMA231(t20,__c1i2,t12);	 FMA231(t21,__c1i2,t13);	/* mpy by E^4=i is inlined here... */\
	 FMA231(t20,__c1i2,_a );	FNMA231(t21,__c1i2,_b );\
	_d = t04; t29 = t05;\
	 FMA231(_c ,__c1i2,t04);	 FMA231(_e ,__c1i2,t05);\
	FNMA231(_c ,__c1i2,_d );	FNMA231(_e ,__c1i2,t29);\
	 								_e = t10-t11;	_f = t10+t11;	/* [Block 2] twiddle mpy by -E^6 is here... */\
	/* Write outputs back to main array: */\
	__B6r= t12;		__B6i= t13;\
	__B7r= _a ;		__B7i= _b ;\
	__B4r= t04;		__B4i= t05;\
	__B5r= _d ;		__B5i= t29;\
\
	/*...Block 2: t2,10,18,26 */\
	/* Swapped upper right and lower left FMA-pairs of quartet below to hide latency: */\
	_c = t18; _a = t27;\
	FNMA231(t19,__sc,_c );		 FMS231(t26,__sc,_a );		/* twiddle mpy by E^1 */\
	_d = t19; _b = t26;\
	 FMA231(t18,__sc,_d );		 FMA231(t27,__sc,_b );		/* twiddle mpy by E^3 */\
	/* Swapped upper left and lower right FMA-pairs of octet below to hide latency: */\
	t18 = _c;	t10 = t02;\
	 FMA231(_a ,__c31,t18);		 FMA231(_e ,__c2i2,t02);		/* ...and get E^6=[i-1]/sqrt2 by flipping signs here. */\
	t19 = _d;	t11 = t03;\
	 FMA231(_b ,__c31,t19);		 FMA231(_f ,__c2i2,t03);\
	FNMA231(_a ,__c31,_c );		FNMA231(_e ,__c2i2,t10);\
	FNMA231(_b ,__c31,_d );		FNMA231(_f ,__c2i2,t11);\
	_a = t02; _b = t03;\
	 FMA231(t18,__c1_c,t02);	 FMA231(t19,__c1_c,t03);\
	FNMA231(t18,__c1_c,_a );	FNMA231(t19,__c1_c,_b );\
	_e = t10; _f = t11;\
	FNMA231(_d ,__c1_c,t10);	 FMA231(_c ,__c1_c,t11);	/* mpy by E^4=i is inlined here... */\
	 FMA231(_d ,__c1_c,_e );	FNMA231(_c ,__c1_c,_f );\
	 								_c = t14+t15;	_d = t15-t14;	/* [Block 4] twiddle mpy by -E^6 is here... */\
	/* Write outputs back to main array: */\
	__B8r= t02;		__B8i= t03;\
	__B9r= _a ;		__B9i= _b ;\
	__BAr= t10;		__BAi= t11;\
	__BBr= _e ;		__BBi= _f ;\
\
	/*...Block 4: t6,14,22,30 */\
	_e = t23; _f = t22;\
	 FMS231(t22,__sc,_e );		 FMA231(t23,__sc,_f );		/* twiddle mpy by E^3 */\
	_a = t30; _b = t31;\
	FNMA231(t31,__sc,_a );		 FMA231(t30,__sc,_b );		/* twiddle mpy by E^1 =__c*[ -E^9... */\
	/* Rowwise-shuffled FMA-pairs of octet below to hide latency: */\
	t15= t07;	t14= t06;\
	FNMA231(_d ,__c2i2,t07);	FNMA231(_c ,__c2i2,t06);	/* ...and get E^6=[i-1]/sqrt2 by flipping signs here. */\
	t22= _e; t23= _f;\
	FNMA231(_a ,__c31 ,t22);	FNMA231(_b ,__c31 ,t23);	/* ...and get E^9 by flipping signs here. */\
	 FMA231(_a ,__c31 ,_e );	 FMA231(_b ,__c31 ,_f );\
	 FMA231(_c ,__c2i2,t14);	 FMA231(_d ,__c2i2,t15);\
	_a = t06; _b = t07;\
	 FMA231(t22,__c1_c,t06);	 FMA231(t23,__c1_c,t07);\
	FNMA231(t22,__c1_c,_a );	FNMA231(t23,__c1_c,_b );\
	_c = t14; _d = t15;\
	FNMA231(_f ,__c1_c,t14);	 FMA231(_e ,__c1_c,t15);	/* mpy by E^4=i is inlined here... */\
	 FMA231(_f ,__c1_c,_c );	FNMA231(_e ,__c1_c,_d );\
	/* Write outputs back to main array: */\
	__BCr= t06;		__BCi= t07;\
	__BDr= _a ;		__BDi= _b ;\
	__BEr= t14;		__BEi= t15;\
	__BFr= _c ;		__BFi= _d ;\
}

/************************************/
/*** DIT analogs of above macros: ***/
/************************************/
/* Version #1 uses 4-operand FMA4 syntax, similar to what we might use for SIMD assembler on SSE5-supporting AMD CPUs.

Total: 174 FMA [a whopping 112 of which have a trivial unity multiplicand], 34 MUL.

Extra MUL cost vs DIF is due to the post-twiddling, which is less amenable to "multiplier absorption".
*/
#define RADIX_16_DIT_FMA4(\
	__A0r,__A0i,__A1r,__A1i,__A2r,__A2i,__A3r,__A3i,__A4r,__A4i,__A5r,__A5i,__A6r,__A6i,__A7r,__A7i,__A8r,__A8i,__A9r,__A9i,__AAr,__AAi,__ABr,__ABi,__ACr,__ACi,__ADr,__ADi,__AEr,__AEi,__AFr,__AFi,\
	__B0r,__B0i,__B1r,__B1i,__B2r,__B2i,__B3r,__B3i,__B4r,__B4i,__B5r,__B5i,__B6r,__B6i,__B7r,__B7i,__B8r,__B8i,__B9r,__B9i,__BAr,__BAi,__BBr,__BBi,__BCr,__BCi,__BDr,__BDi,__BEr,__BEi,__BFr,__BFi,\
	__c1 ,__s1 ,__c2 ,__s2 ,__c3 ,__s3 ,__c4 ,__s4 ,__c5 ,__s5 ,__c6 ,__s6 ,__c7 ,__s7 ,__c8 ,__s8 ,__c9 ,__s9 ,__cA ,__sA ,__cB ,__sB ,__cC ,__sC ,__cD ,__sD ,__cE ,__sE ,__cF ,__sF ,\
	__c,__s)\
{\
	double t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32;\
	/* Note no cosine ratios needed in DIT version. */\
	/* The divides in the following sincos data will be precomputed in an optimized version: */\
	double __t1 = __s1/__c1;\
	double __t2 = __s2/__c2;\
	double __t3 = __s3/__c3;\
	double __t4 = __s4/__c4;\
	double __t5 = __s5/__c5;\
	double __t6 = __s6/__c6;\
	double __t7 = __s7/__c7;\
	double __t8 = __s8/__c8;\
	double __t9 = __s9/__c9;\
	double __tA = __sA/__cA;\
	double __tB = __sB/__cB;\
	double __tC = __sC/__cC;\
	double __tD = __sD/__cD;\
	double __tE = __sE/__cE;\
	double __tF = __sF/__cF;\
	double __sc = __s/__c;\
	/* Gather the needed data and do first set of four length-4 transforms: */\
	/*...Block 1: */\
	t1 =__A0r;	t2 =__A0i;\
	rt =__A1r;	it =__A1i;\
	t3 =FNMA4(1.0,rt,t1);		t1 = FMA4(1.0,rt,t1);\
	t4 =FNMA4(1.0,it,t2);		t2 = FMA4(1.0,it,t2);\
\
	t5 =__A2r;	t6 =__A2i;\
	rt =__A3r;	it =__A3i;\
	t7 =FNMA4(1.0,rt,t5);		t5 = FMA4(1.0,rt,t5);\
	t8 =FNMA4(1.0,it,t6);		t6 = FMA4(1.0,it,t6);\
\
	rt =t5;	t5 =FNMA4(1.0,rt,t1);	t1 = FMA4(1.0,rt,t1);\
	it =t6;	t6 =FNMA4(1.0,it,t2);	t2 = FMA4(1.0,it,t2);\
\
	it =t7;	t7 =FNMA4(1.0,t8,t3);	t3 = FMA4(1.0,t8,t3);\
			t8 = FMA4(1.0,it,t4);	t4 =FNMA4(1.0,it,t4);\
\
	/*...Block 2: */\
	t9 =__A4r;	t10=__A4i;\
	rt =__A5r;	it =__A5i;\
	t11=FNMA4(1.0,rt,t9 );		t9 = FMA4(1.0,rt,t9 );\
	t12=FNMA4(1.0,it,t10);		t10= FMA4(1.0,it,t10);\
\
	t13=__A6r;	t14=__A6i;\
	rt =__A7r;	it =__A7i;\
	t15=FNMA4(1.0,rt,t13);		t13= FMA4(1.0,rt,t13);\
	t16=FNMA4(1.0,it,t14);		t14= FMA4(1.0,it,t14);\
\
	rt =t13;	t13=FNMA4(1.0,rt,t9 );	t9 = FMA4(1.0,rt,t9 );\
	it =t14;	t14=FNMA4(1.0,it,t10);	t10= FMA4(1.0,it,t10);\
\
	it =t15;	t15=FNMA4(1.0,t16,t11);	t11= FMA4(1.0,t16,t11);\
				t16= FMA4(1.0, it,t12);	t12=FNMA4(1.0, it,t12);\
\
	/*...Block 3: */\
	t17=__A8r;	t18=__A8i;\
	rt =__A9r;	it =__A9i;\
	t19=FNMA4(1.0,rt,t17);		t17= FMA4(1.0,rt,t17);\
	t20=FNMA4(1.0,it,t18);		t18= FMA4(1.0,it,t18);\
\
	t21=__AAr;	t22=__AAi;\
	rt =__ABr;	it =__ABi;\
	t23=FNMA4(1.0,rt,t21);		t21= FMA4(1.0,rt,t21);\
	t24=FNMA4(1.0,it,t22);		t22= FMA4(1.0,it,t22);\
\
	rt =t21;	t21=FNMA4(1.0,rt,t17);	t17= FMA4(1.0,rt,t17);\
	it =t22;	t22=FNMA4(1.0,it,t18);	t18= FMA4(1.0,it,t18);\
\
	it =t23;	t23=FNMA4(1.0,t24,t19);	t19= FMA4(1.0,t24,t19);\
				t24= FMA4(1.0, it,t20);	t20=FNMA4(1.0, it,t20);\
\
	/*...Block 4: */\
	t25=__ACr;	t26=__ACi;\
	rt =__ADr;	it =__ADi;\
	t27=FNMA4(1.0,rt,t25);		t25= FMA4(1.0,rt,t25);\
	t28=FNMA4(1.0,it,t26);		t26= FMA4(1.0,it,t26);\
\
	t29=__AEr;	t30=__AEi;\
	rt =__AFr;	it =__AFi;\
	t31=FNMA4(1.0,rt,t29);		t29= FMA4(1.0,rt,t29);\
	t32=FNMA4(1.0,it,t30);		t30= FMA4(1.0,it,t30);\
\
	rt =t29;	t29=FNMA4(1.0,rt,t25);	t25= FMA4(1.0,rt,t25);\
	it =t30;	t30=FNMA4(1.0,it,t26);	t26= FMA4(1.0,it,t26);\
\
	it =t31;	t31=FNMA4(1.0,t32,t27);	t27= FMA4(1.0,t32,t27);\
				t32= FMA4(1.0, it,t28);	t28=FNMA4(1.0, it,t28);\
\
/*...and now do four more radix-4 transforms, including the internal twiddle factors: */\
\
	/*...Block 1: t1,9,17,25 */\
	rt =t9 ;	t9 =FNMA4(1.0,rt,t1 );	t1 = FMA4(1.0,rt,t1 );\
	it =t10;	t10=FNMA4(1.0,it,t2 );	t2 = FMA4(1.0,it,t2 );\
\
	rt =t25;	t25=FNMA4(1.0,rt,t17);	t17= FMA4(1.0,rt,t17);\
	it =t26;	t26=FNMA4(1.0,it,t18);	t18= FMA4(1.0,it,t18);\
\
	__B0r = FMA4(1.0,t17,t1 );		__B0i = FMA4(1.0,t18,t2 );\
	t1    =FNMA4(1.0,t17,t1 );		t2    =FNMA4(1.0,t18,t2 );\
	rt    = FMA4(t2,__t8,t1 );		it    =FNMA4(t1,__t8,t2 );\
	__B8r = rt *__c8;	__B8i = it *__c8;/*** Are these pre-output MULs in the DIF/post-twiddle algo avoidable? ***/\
\
	rt    = FMA4(1.0,t26,t9 );		it	  =FNMA4(1.0,t25,t10);	/* mpy by E^-4 = -I is inlined here... */\
	t9    =FNMA4(1.0,t26,t9 );		t10	  = FMA4(1.0,t25,t10);\
	t25= FMA4(it ,__t4,rt );		t26=FNMA4(rt ,__t4,it );\
	rt = FMA4(t10,__tC,t9 );		it =FNMA4(t9 ,__tC,t10);\
	__B4r = t25*__c4;	__B4i = t26*__c4;\
	__BCr = rt *__cC;	__BCi = it *__cC;\
\
	/*...Block 3: t5,13,21,29 */\
	rt =t13;	t13=FNMA4(t14,1.0,t5 );	t5 = FMA4(t14,1.0,t5 );		/* twiddle mpy by E^4 =-I */\
				t14= FMA4( rt,1.0,t6 );	t6 =FNMA4( rt,1.0,t6 );\
\
	rt = FMA4(1.0,t21,t22);				t22=FNMA4(1.0,t21,t22);	t21=rt;	/* twiddle mpy by E^-2 */\
	rt = FMS4(1.0,t29,t30);				it = FMA4(1.0,t29,t30);		/* twiddle mpy by E^2 = -E^-6 is here... */\
	t29= FMA4(1.0,t21,rt );				t21= FMS4(1.0,t21,rt );			/* ...and get E^-6 by flipping signs here. */\
	t30= FMA4(1.0,t22,it );				t22= FMS4(1.0,t22,it );\
\
	rt  = FMA4(t21,ISRT2,t5 );			it = FMA4(t22,ISRT2,t6 );\
	t5  =FNMA4(t21,ISRT2,t5 );			t6 =FNMA4(t22,ISRT2,t6 );\
	__B2r =( FMA4(it ,__t2,rt ))*__c2;	__B2i =(FNMA4(rt ,__t2,it ))*__c2;\
	__BAr =( FMA4(t6 ,__tA,t5 ))*__cA;	__BAi =(FNMA4(t5 ,__tA,t6 ))*__cA;\
\
	rt  = FMA4(t30,ISRT2,t13);			it  =FNMA4(t29,ISRT2,t14);	/* mpy by E^-4 =-I is inlined here... */\
	t13	=FNMA4(t30,ISRT2,t13);			t14 = FMA4(t29,ISRT2,t14);\
	__B6r =( FMA4(it ,__t6,rt ))*__c6;	__B6i =(FNMA4(rt ,__t6,it ))*__c6;\
	__BEr =( FMA4(t14,__tE,t13))*__cE;	__BEi =(FNMA4(t13,__tE,t14))*__cE;\
\
	/*...Block 2: t3,11,19,27 */\
	rt = FMA4(1.0,t12,t11);				it = FMS4(1.0,t12,t11);		/* twiddle mpy by E^-2 */\
	t11=FNMA4(rt ,ISRT2,t3 );			t3 = FMA4(rt ,ISRT2,t3 );\
	t12=FNMA4(it ,ISRT2,t4 );			t4 = FMA4(it ,ISRT2,t4 );\
\
	rt = FMA4(t20,__sc,t19);			t20=FNMA4(t19,__sc,t20);	t19=rt;	/* twiddle mpy by E^-1 */\
	rt = FMA4(t27,__sc,t28);			it = FMS4(t28,__sc,t27);		/* twiddle mpy by E^-3 */\
	rt *= __c;	it *= __c;	/* Doing things this way incurs just 2 extra MULs, other 2 absorbed into the ensuing 4 FMAs */\
	t27= FMS4(__c,t19,rt );				t19= FMA4(__c,t19,rt );\
	t28= FMS4(__c,t20,it );				t20= FMA4(__c,t20,it );\
\
	rt  = FMA4(1.0,t3 ,t19 );			it  = FMA4(1.0,t4 ,t20 );\
	t3  = FMS4(1.0,t3 ,t19 );			t4  = FMS4(1.0,t4 ,t20 );\
	__B1r =( FMA4(it ,__t1,rt ))*__c1;	__B1i =(FNMA4(rt ,__t1,it ))*__c1;\
	__B9r =( FMA4(t4 ,__t9,t3 ))*__c9;	__B9i =(FNMA4(t3 ,__t9,t4 ))*__c9;\
\
	rt  = FMA4(1.0,t11,t28);			it  = FMS4(1.0,t12,t27);	/* mpy by E^-4 = -I is inlined here... */\
	t11 = FMS4(1.0,t11,t28);			t12 = FMA4(1.0,t12,t27);\
	__B5r =( FMA4(it ,__t5,rt ))*__c5;	__B5i =(FNMA4(rt ,__t5,it ))*__c5;\
	__BDr =( FMA4(t12,__tD,t11))*__cD;	__BDi =(FNMA4(t11,__tD,t12))*__cD;\
\
	/*...Block 4: t7,15,23,31 */\
	rt = FMS4(1.0,t15,t16);				it = FMA4(1.0,t15,t16);		/* twiddle mpy by E^2 = -E^-6 is here... */\
	t15= FMA4(rt ,ISRT2,t7 );			t7 =FNMA4(rt ,ISRT2,t7 );	/* ...and get E^6=(i-1)/sqrt by flipping signs here. */\
	t16= FMA4(it ,ISRT2,t8 );			t8 =FNMA4(it ,ISRT2,t8 );\
\
	rt = FMA4(t23,__sc,t24);			t24= FMS4(t24,__sc,t23);	t23=rt;	/* twiddle mpy by E^-3 */\
	rt = FMA4(t32,__sc,t31);			it =FNMA4(t31,__sc,t32);		/* twiddle mpy by E^-1 = -E^-9... */\
	rt *= __c;	it *= __c;	/* Doing things this way incurs just 2 extra MULs, other 2 absorbed into the ensuing 4 FMAs */\
	t31= FMA4(__c,t23,rt );				t23= FMS4(__c,t23,rt );			/* ...and get E^9 by flipping signs here. */\
	t32= FMA4(__c,t24,it );				t24= FMS4(__c,t24,it );\
\
	rt  = FMA4(1.0,t7 ,t23);			it  = FMA4(1.0,t8 ,t24);\
	t7  = FMS4(1.0,t7 ,t23);			t8  = FMS4(1.0,t8 ,t24);\
	__B3r =( FMA4(it ,__t3,rt ))*__c3;	__B3i =(FNMA4(rt ,__t3,it ))*__c3;\
	__BBr =( FMA4(t8 ,__tB,t7 ))*__cB;	__BBi =(FNMA4(t7 ,__tB,t8 ))*__cB;\
\
	rt  = FMA4(1.0,t15,t32);			it  = FMS4(1.0,t16,t31);	/* mpy by E^-4 = -I is inlined here... */\
	t15 = FMS4(1.0,t15,t32);			t16 = FMA4(1.0,t16,t31);\
	__B7r =( FMA4(it ,__t7,rt ))*__c7;	__B7i =(FNMA4(rt ,__t7,it ))*__c7;\
	__BFr =( FMA4(t16,__tF,t15))*__cF;	__BFi =(FNMA4(t15,__tF,t16))*__cF;\
}

// Version #2 uses a 3-operand FMA3 syntax in which the RIGHTmost of the 3 fused mul-add inputs -
// that is, the addend - is overwritten with the result:
#define RADIX_16_DIT_FMA_1(\
	__A0r,__A0i,__A1r,__A1i,__A2r,__A2i,__A3r,__A3i,__A4r,__A4i,__A5r,__A5i,__A6r,__A6i,__A7r,__A7i,__A8r,__A8i,__A9r,__A9i,__AAr,__AAi,__ABr,__ABi,__ACr,__ACi,__ADr,__ADi,__AEr,__AEi,__AFr,__AFi,\
	__B0r,__B0i,__B1r,__B1i,__B2r,__B2i,__B3r,__B3i,__B4r,__B4i,__B5r,__B5i,__B6r,__B6i,__B7r,__B7i,__B8r,__B8i,__B9r,__B9i,__BAr,__BAi,__BBr,__BBi,__BCr,__BCi,__BDr,__BDi,__BEr,__BEi,__BFr,__BFi,\
	__c1 ,__s1 ,__c2 ,__s2 ,__c3 ,__s3 ,__c4 ,__s4 ,__c5 ,__s5 ,__c6 ,__s6 ,__c7 ,__s7 ,__c8 ,__s8 ,__c9 ,__s9 ,__cA ,__sA ,__cB ,__sB ,__cC ,__sC ,__cD ,__sD ,__cE ,__sE ,__cF ,__sF ,\
	__c,__s)\
{\
	double t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32;\
	double rt,it;\
	/* Note no cosine ratios needed in DIT version. */\
	/* The divides in the following sincos data will be precomputed in an optimized version: */\
	double __t1 = __s1/__c1;\
	double __t2 = __s2/__c2;\
	double __t3 = __s3/__c3;\
	double __t4 = __s4/__c4;\
	double __t5 = __s5/__c5;\
	double __t6 = __s6/__c6;\
	double __t7 = __s7/__c7;\
	double __t8 = __s8/__c8;\
	double __t9 = __s9/__c9;\
	double __tA = __sA/__cA;\
	double __tB = __sB/__cB;\
	double __tC = __sC/__cC;\
	double __tD = __sD/__cD;\
	double __tE = __sE/__cE;\
	double __tF = __sF/__cF;\
	double __sc = __s/__c;\
	/* Gather the needed data and do first set of four length-4 transforms: */\
	/*...Block 1: Each complex radix-4 block needs 11 registers (8 t-temps plus rt,it, plus const 1.0) */\
	/* Restructure to do batches of 4 FMA side-by-side	Swap t7,t8 on input so they come out "right" w/o extra reg-copying */\
	/*...Block 1: */\
	t1 =__A0r;	t2 =__A0i;\
	rt =__A1r;	it =__A1i;\
	t3 = t1 ;	t4 = t2 ;\
	FNMA231(1.0,rt,t3 );		FMA231(1.0,rt,t1 );\
	FNMA231(1.0,it,t4 );		FMA231(1.0,it,t2 );\
\
	t5 =__A2r;	t6 =__A2i;\
	rt =__A3r;	it =__A3i;\
	t7 = t5 ;	t8 = t6 ;\
	FNMA231(1.0,rt,t7 );		FMA231(1.0,rt,t5 );\
	FNMA231(1.0,it,t8 );		FMA231(1.0,it,t6 );\
\
	rt = t5 ;	it = t6 ;\
	t5 = t1 ;	t6 = t2 ;\
	FNMA231(1.0,rt,t5 );		FMA231(1.0,rt,t1 );\
	FNMA231(1.0,it,t6 );		FMA231(1.0,it,t2 );\
\
	it = t7 ;	rt = t8 ;\
	t7 = t3 ;	t8 = t4 ;\
	FNMA231(1.0,rt,t7 );		 FMA231(1.0,rt,t3 );\
	 FMA231(1.0,it,t8 );		FNMA231(1.0,it,t4 );\
\
	/*...Block 2: */\
	t9 =__A4r;	t10=__A4i;\
	rt =__A5r;	it =__A5i;\
	t11= t9 ;	t12= t10;\
	FNMA231(1.0,rt,t11);		FMA231(1.0,rt,t9 );\
	FNMA231(1.0,it,t12);		FMA231(1.0,it,t10);\
\
	t13=__A6r;	t14=__A6i;\
	rt =__A7r;	it =__A7i;\
	t15= t13;	t16= t14;\
	FNMA231(1.0,rt,t15);		FMA231(1.0,rt,t13);\
	FNMA231(1.0,it,t16);		FMA231(1.0,it,t14);\
\
	rt = t13;	it = t14;\
	t13= t9 ;	t14= t10;\
	FNMA231(1.0,rt,t13);		FMA231(1.0,rt,t9 );\
	FNMA231(1.0,it,t14);		FMA231(1.0,it,t10);\
\
	it = t15;	rt = t16;\
	t15= t11;	t16= t12;\
	FNMA231(1.0,rt,t15);		 FMA231(1.0,rt,t11);\
	 FMA231(1.0,it,t16);		FNMA231(1.0,it,t12);\
\
	/*...Block 3: */\
	t17=__A8r;	t18=__A8i;\
	rt =__A9r;	it =__A9i;\
	t19= t17;	t20= t18;\
	FNMA231(1.0,rt,t19);		FMA231(1.0,rt,t17);\
	FNMA231(1.0,it,t20);		FMA231(1.0,it,t18);\
\
	t21=__AAr;	t22=__AAi;\
	rt =__ABr;	it =__ABi;\
	t23= t21;	t24= t22;\
	FNMA231(1.0,rt,t23);		FMA231(1.0,rt,t21);\
	FNMA231(1.0,it,t24);		FMA231(1.0,it,t22);\
\
	rt = t21;	it = t22;\
	t21= t17;	t22= t18;\
	FNMA231(1.0,rt,t21);		FMA231(1.0,rt,t17);\
	FNMA231(1.0,it,t22);		FMA231(1.0,it,t18);\
\
	it = t23;	rt = t24;\
	t23 = t19;	t24 = t20;\
	FNMA231(1.0,rt,t23);		 FMA231(1.0,rt,t19);\
	 FMA231(1.0,it,t24);		FNMA231(1.0,it,t20);\
\
	/*...Block 4: */\
	t25=__ACr;	t26=__ACi;\
	rt =__ADr;	it =__ADi;\
	t27= t25;	t28= t26;\
	FNMA231(1.0,rt,t27);		FMA231(1.0,rt,t25);\
	FNMA231(1.0,it,t28);		FMA231(1.0,it,t26);\
\
	t29=__AEr;	t30=__AEi;\
	rt =__AFr;	it =__AFi;\
	t31= t29;	t32= t30;\
	FNMA231(1.0,rt,t31);		FMA231(1.0,rt,t29);\
	FNMA231(1.0,it,t32);		FMA231(1.0,it,t30);\
\
	rt = t29;	it = t30;\
	t29= t25;	t30= t26;\
	FNMA231(1.0,rt,t29);		FMA231(1.0,rt,t25);\
	FNMA231(1.0,it,t30);		FMA231(1.0,it,t26);\
\
	it = t31;	rt = t32;\
	t31= t27;	t32= t28;\
	FNMA231(1.0,rt,t31);		 FMA231(1.0,rt,t27);\
	 FMA231(1.0,it,t32);		FNMA231(1.0,it,t28);\
\
/*...and now do four more radix-4 transforms, including the internal twiddle factors: */\
\
	/*...Block 1: t1,9,17,25 */\
	rt = t1;	it = t2;\
	 FMA231(1.0,t9 ,t1 );		FNMA132(1.0,t9 ,rt );\
	 FMA231(1.0,t10,t2 );		FNMA132(1.0,t10,it );\
\
	rt =t17;	it =t18;\
	 FMA231(1.0,t25,t17);		FNMA132(1.0,t25,rt );\
	 FMA231(1.0,t26,t18);		FNMA132(1.0,t26,it );\
\
	rt =t1;	it =t2;\
	 FMA231(1.0,t17,rt );		 FMA231(1.0,t18,it );\
	FNMA231(1.0,t17,t1 );		FNMA231(1.0,t18,t2 );\
	__B0r = rt;			__B0i = it;\
	rt =t1;\
	 FMA231(t2,__t8,t1 );		FNMA231(rt,__t8,t2 );\
	__B8r = t1 *__c8;	__B8i = t2 *__c8;\
\
	rt =t9 ;	it =t10;\
	FNMA231(1.0,t26,t9 );		 FMA132(1.0,t26,rt );	/* mpy by E^-4 = -I is inlined here... */\
	 FMA231(1.0,t25,t10);		FNMA132(1.0,t25,it );\
	rt =t26;	it =t9 ;\
	 FMA231(__t4,t25,t26);		FNMA231(__t4,rt ,t25);\
	 FMA231(__tC,t10,t9 );		FNMA231(__tC,it ,t10);\
	__B4r = t26*__c4;	__B4i = t25*__c4;\
	__BCr = t9 *__cC;	__BCi = t10*__cC;\
\
	/*...Block 3: t5,13,21,29 */\
	rt = t5;	it = t6;\
	 FMA231(1.0,t14,t5 );		FNMA132(1.0,t14,rt );	/* twiddle mpy by E^4 =-I */\
	FNMA231(1.0,t13,t6 );		 FMA132(1.0,t13,it );	/* t13,t14 swapped w.r.to FMA4 here */\
	rt = t22;	it = t30;\
	FNMA231(1.0,t21,t22);		 FMA132(1.0,t21,rt );	/* twiddle mpy by E^-2; note FMA132 for right MUL saves us the t21=rt copy */\
	 FMA231(1.0,t29,t30);		 FMS132(1.0,t29,it );	/* twiddle mpy by E^2 = -E^-6 is here... */\
	rt = t21;	it = t22;\
	FNMA231(1.0,t29,t21);		 FMA132(1.0,t29,rt );	/* ...and get E^-6 by flipping signs here. */\
	FNMA231(1.0,t30,t22);		 FMA132(1.0,t30,it );\
\
	rt = t5;	it = t6;\
	FNMA231(ISRT2,t21,t5 );		 FMA132(ISRT2,t21,rt );\
	FNMA231(ISRT2,t22,t6 );		 FMA132(ISRT2,t22,it );\
	rt =t21;	it = t5;\
	 FMA231(__t2,t22,t21);		FNMA231(__t2,rt ,t22);\
	 FMA231(__tA,t6 ,t5 );		FNMA231(__tA,it ,t6 );\
	__B2r = t21*__c2;	__B2i = t22*__c2;\
	__BAr = t5 *__cA;	__BAi = t6 *__cA;\
\
	rt =t13;	it =t14;\
	 FMA231(ISRT2,t29,t13);		FNMA132(ISRT2,t29,rt );	/* mpy by E^-4 =-I is inlined here... */\
	FNMA231(ISRT2,t30,t14);		 FMA132(ISRT2,t30,it );\
	rt =t30;	it =t14;\
	 FMA231(__t6,t29,t30);		FNMA231(__t6,rt ,t29);\
	 FMA231(__tE,t13,t14);		FNMA231(__tE,it ,t13);\
	__B6r = t30*__c6;	__B6i = t29*__c6;\
	__BEr = t14*__cE;	__BEi = t13*__cE;\
\
	/*...Block 2: t3,11,19,27 */\
	rt =t11;\
	 FMA231(1.0,t12,t11);		FNMA231(1.0,rt ,t12);	/* twiddle mpy by E^-2 */\
	rt =t3 ;	it =t4 ;\
	 FMA231(ISRT2,t11,t3 );		FNMA132(ISRT2,t11,rt );\
	 FMA231(ISRT2,t12,t4 );		FNMA132(ISRT2,t12,it );\
\
	rt = t19;	it = t27;\
	 FMA231(__sc,t20,t19);		FNMA231(__sc,rt ,t20);	/* twiddle mpy by E^-1 */\
	 FMA132(__sc,t27,t28);		 FMS132(__sc,t28,it );	/* twiddle mpy by E^-3; Use FMA132 for both because we need to swap t27/t28 on way to outputs */\
	t27 *= __c;	t28 *= __c;	/* Doing things this way incurs just 2 extra MULs, other 2 absorbed into the ensuing 4 FMAs */\
	rt = t27;	it = t28;\
	 FMS231(__c,t19,t27);		 FMA132(__c,t19,rt );	/* Use FMA132 for rhs MULs to put outputs into t19,t20 */\
	 FMS231(__c,t20,t28);		 FMA132(__c,t20,it );\
\
	rt = t3; it = t4;\
	 FMS132(1.0,t3 ,t19 );		 FMS132(1.0,t4 ,t20 );	/* outs in t3 ,t4  */\
	 FMA231(1.0,rt ,t19 );		 FMA231(1.0,it ,t20 );	/* outs in t19,t20 */\
	rt = t19;	it = t3;\
	 FMA231(__t1,t20,t19);		FNMA231(__t1,rt ,t20);\
	 FMA231(__t9,t4 ,t3 );		FNMA231(__t9,it ,t4 );\
	__B1r = t19*__c1;	__B1i = t20*__c1;\
	__B9r = t3 *__c9;	__B9i = t4 *__c9;\
\
	rt = t11; it = t12;\
	 FMS132(1.0,t11,t28);		 FMA132(1.0,t12,t27);	/* outs in t11,t12 */\
	 FMA231(1.0,rt ,t28);		 FMS231(1.0,it ,t27);	/* outs in t28,t27 */\
	rt = t28;	it = t11;\
	 FMA231(__t5,t27,t28);		FNMA231(__t5,rt ,t27);\
	 FMA231(__tD,t12,t11);		FNMA231(__tD,it ,t12);\
	__B5r = t28*__c5;	__B5i = t27*__c5;\
	__BDr = t11*__cD;	__BDi = t12*__cD;\
\
	/*...Block 4: t7,15,23,31 */\
	rt = t16;\
	 FMA231(1.0,t15,t16);		 FMS132(1.0,t15,rt );	/* twiddle mpy by E^2 = -E^-6 is here... */\
	rt =t7 ;	it =t8 ;\
	FNMA231(ISRT2,t15,t7 );		 FMA132(ISRT2,t15,rt );	/* ...and get E^6=(i-1)/sqrt by flipping signs here. */\
	FNMA231(ISRT2,t16,t8 );		 FMA132(ISRT2,t16,it );\
\
	rt = t23;\
	FMA132(__sc,t23,t24);		 FMS132(__sc,t24,rt );	/* twiddle mpy by E^-3 */\
	rt = t31;\
	FMA231(__sc,t32,t31);		FNMA231(__sc,rt ,t32);	/* twiddle mpy by E^-1 = -E^-9... */\
	t31 *= __c;	t32 *= __c;	/* Doing things this way incurs just 2 extra MULs, other 2 absorbed into the ensuing 4 FMAs */\
	rt = t31;	it = t32;\
	 FMA231(__c,t23,t31);		 FMS132(__c,t23,rt );	/* ...and get E^9 by flipping signs here. */\
	 FMA231(__c,t24,t32);		 FMS132(__c,t24,it );\
\
	rt = t7; it = t8;\
	 FMS132(1.0,t7 ,t23 );		 FMS132(1.0,t8 ,t24 );	/* outs in t7 ,t8  */\
	 FMA231(1.0,rt ,t23 );		 FMA231(1.0,it ,t24 );	/* outs in t23,t24 */\
	rt = t23;	it = t7;\
	 FMA231(__t3,t24,t23);		FNMA231(__t3,rt ,t24);\
	 FMA231(__tB,t8 ,t7 );		FNMA231(__tB,it ,t8 );\
	__B3r = t23*__c3;	__B3i = t24*__c3;\
	__BBr = t7 *__cB;	__BBi = t8 *__cB;\
\
	rt = t15; it = t16;\
	 FMS132(1.0,t15,t32);		 FMA132(1.0,t16,t31);	/* outs in t15,t16 */\
	 FMA231(1.0,rt ,t32);		 FMS231(1.0,it ,t31);	/* outs in t32,t31 */\
	rt = t32;	it = t15;\
	 FMA231(__t7,t31,t32);		FNMA231(__t7,rt ,t31);\
	 FMA231(__tF,t16,t15);		FNMA231(__tF,it ,t16);\
	__B7r = t32*__c7;	__B7i = t31*__c7;\
	__BFr = t15*__cF;	__BFi = t16*__cF;\
}

// Version #3 is v2, restructured to take account of 16 SIMd registers and 5-cycle FMA latency of Intel Haswell:
#define RADIX_16_DIT_FMA(\
	__A0r,__A0i,__A1r,__A1i,__A2r,__A2i,__A3r,__A3i,__A4r,__A4i,__A5r,__A5i,__A6r,__A6i,__A7r,__A7i,__A8r,__A8i,__A9r,__A9i,__AAr,__AAi,__ABr,__ABi,__ACr,__ACi,__ADr,__ADi,__AEr,__AEi,__AFr,__AFi,\
	__B0r,__B0i,__B1r,__B1i,__B2r,__B2i,__B3r,__B3i,__B4r,__B4i,__B5r,__B5i,__B6r,__B6i,__B7r,__B7i,__B8r,__B8i,__B9r,__B9i,__BAr,__BAi,__BBr,__BBi,__BCr,__BCi,__BDr,__BDi,__BEr,__BEi,__BFr,__BFi,\
	__c1 ,__s1 ,__c2 ,__s2 ,__c3 ,__s3 ,__c4 ,__s4 ,__c5 ,__s5 ,__c6 ,__s6 ,__c7 ,__s7 ,__c8 ,__s8 ,__c9 ,__s9 ,__cA ,__sA ,__cB ,__sB ,__cC ,__sC ,__cD ,__sD ,__cE ,__sE ,__cF ,__sF ,\
	__c,__s)\
{\
	double t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32;\
	double rt,it,re,im;\
	/* Note no cosine ratios needed in DIT version. */\
	/* The divides in the following sincos data will be precomputed in an optimized version: */\
	double __t1 = __s1;	/*/__c1;*;	/*/\
	double __t2 = __s2;	/*/__c2;*;	/*/\
	double __t3 = __s3;	/*/__c3;*;	/*/\
	double __t4 = __s4;	/*/__c4;*;	/*/\
	double __t5 = __s5;	/*/__c5;*;	/*/\
	double __t6 = __s6;	/*/__c6;*;	/*/\
	double __t7 = __s7;	/*/__c7;*;	/*/\
	double __t8 = __s8;	/*/__c8;*;	/*/\
	double __t9 = __s9;	/*/__c9;*;	/*/\
	double __tA = __sA;	/*/__cA;*;	/*/\
	double __tB = __sB;	/*/__cB;*;	/*/\
	double __tC = __sC;	/*/__cC;*;	/*/\
	double __tD = __sD;	/*/__cD;*;	/*/\
	double __tE = __sE;	/*/__cE;*;	/*/\
	double __tF = __sF;	/*/__cF;*;	/*/\
	double __sc = __s /__c;\
	/* Gather the needed data and do first set of four length-4 transforms: */\
	/* Each complex radix-4 block needs 11 registers (8 t-temps plus rt,it, plus const 1.0), though we add 2 */\
	/* more tmps re,im to allow independent Re/Im subsections to overlap w/o introducing false dependencies. */\
	/* Restructure original FMA3 DIT macro above to put batches of 4 FMA side-by-side; no clue whether the */\
	/* microcode execution engine cares about that, or recognizes that consecutuve pair of such 4-FMA-groups */\
	/* can be performed independtly (we hope so, though, and structure the code to make this clear to the coder. */\
	/*...Block 1:									Swap t7,t8 on input so they come out "right" w/o extra reg-copying */\
	t1  =__A0r;	rt  =__A1r;	t3  = t1 ;				t5  =__A2r;	it  =__A3r;	t8  = t5 ;	/* Re part 1: */\
	FNMA231(1.0,rt ,t3 );	 FMA231(1.0,rt ,t1 );	FNMA231(1.0,it ,t8 );	 FMA231(1.0,it ,t5 );\
	t2  =__A0i;	re  =__A1i;	t4  = t2 ;				t6  =__A2i;	im  =__A3i;	t7  = t6 ;	/* Im part 1: */\
	FNMA231(1.0,re ,t4 );	 FMA231(1.0,re ,t2 );	FNMA231(1.0,im ,t7 );	 FMA231(1.0,im ,t6 );\
\
	rt  = t5 ;										it  = t3 ;							/* Re part 2: */\
	FNMA132(1.0,t5 ,t1 );	 FMA231(1.0,rt ,t1 );	 FMA231(1.0,t7 ,t3 );	FNMA132(1.0,t7 ,it );\
	re  = t6 ;										im  = t4 ;							/* Im part 2: */\
	FNMA132(1.0,t6 ,t2 );	 FMA231(1.0,re ,t2 );	FNMA231(1.0,t8 ,t4 );	 FMA132(1.0,t8 ,im );\
\
	/*...Block 2: */\
	t9  =__A4r;	rt  =__A5r;	t11 = t9 ;				t13 =__A6r;	it  =__A7r;	t16 = t13;	/* Re part 1: */\
	FNMA231(1.0,rt ,t11);	 FMA231(1.0,rt ,t9 );	FNMA231(1.0,it ,t16);	 FMA231(1.0,it ,t13);\
	t10 =__A4i;	re  =__A5i;	t12 = t10;				t14 =__A6i;	im  =__A7i;	t15 = t14;	/* Im part 1: */\
	FNMA231(1.0,re ,t12);	 FMA231(1.0,re ,t10);	FNMA231(1.0,im ,t15);	 FMA231(1.0,im ,t14);\
\
	rt  = t13;										it  = t11;							/* Re part 2: */\
	FNMA132(1.0,t13,t9 );	 FMA231(1.0,rt ,t9 );	 FMA231(1.0,t15,t11);	FNMA132(1.0,t15,it );\
	re  = t14;										im  = t12;							/* Im part 2: */\
	FNMA132(1.0,t14,t10);	 FMA231(1.0,re ,t10);	FNMA231(1.0,t16,t12);	 FMA132(1.0,t16,im );\
\
	/*...Block 3: */\
	t17 =__A8r;	rt  =__A9r;	t19 = t17;				t21 =__AAr;	it  =__ABr;	t24 = t21;	/* Re part 1: */\
	FNMA231(1.0,rt ,t19);	 FMA231(1.0,rt ,t17);	FNMA231(1.0,it ,t24);	 FMA231(1.0,it ,t21);\
	t18 =__A8i;	re  =__A9i;	t20 = t18;				t22 =__AAi;	im  =__ABi;	t23 = t22;	/* Im part 1: */\
	FNMA231(1.0,re ,t20);	 FMA231(1.0,re ,t18);	FNMA231(1.0,im ,t23);	 FMA231(1.0,im ,t22);\
\
	rt  = t21;										it  = t19;							/* Re part 2: */\
	FNMA132(1.0,t21,t17);	 FMA231(1.0,rt ,t17);	 FMA231(1.0,t23,t19);	FNMA132(1.0,t23,it );\
	re  = t22;										im  = t20;							/* Im part 2: */\
	FNMA132(1.0,t22,t18);	 FMA231(1.0,re ,t18);	FNMA231(1.0,t24,t20);	 FMA132(1.0,t24,im );\
\
	/*...Block 4: */\
	t25 =__ACr;	rt  =__ADr;	t27 = t25;				t29 =__AEr;	it  =__AFr;	t32 = t29;	/* Re part 1: */\
	FNMA231(1.0,rt ,t27);	 FMA231(1.0,rt ,t25);	FNMA231(1.0,it ,t32);	 FMA231(1.0,it ,t29);\
	t26 =__ACi;	re  =__ADi;	t28 = t26;				t30 =__AEi;	im  =__AFi;	t31 = t30;	/* Im part 1: */\
	FNMA231(1.0,re ,t28);	 FMA231(1.0,re ,t26);	FNMA231(1.0,im ,t31);	 FMA231(1.0,im ,t30);\
\
	rt  = t29;										it  = t27;							/* Re part 2: */\
	FNMA132(1.0,t29,t25);	 FMA231(1.0,rt ,t25);	 FMA231(1.0,t31,t27);	FNMA132(1.0,t31,it );\
	re  = t30;										im  = t28;							/* Im part 2: */\
	FNMA132(1.0,t30,t26);	 FMA231(1.0,re ,t26);	FNMA231(1.0,t32,t28);	 FMA132(1.0,t32,im );\
\
/*...and now do four more radix-4 transforms, including the internal twiddle factors: */\
\
	/*...Block 1: t1,9,17,25 */\
	rt = t1 ;	it = t2 ;\
	 FMA231(1.0,t9 ,t1 );	FNMA132(1.0,t9 ,rt );	 FMA231(1.0,t10,t2 );	FNMA132(1.0,t10,it );\
	re =t17;	im =t18;\
	 FMA231(1.0,t25,t17);	FNMA132(1.0,t25,re );	 FMA231(1.0,t26,t18);	FNMA132(1.0,t26,im );\
\
	rt =t1 ;	it =t2 ;\
	 FMA231(1.0,t17,rt );	 FMA231(1.0,t18,it );	FNMA231(1.0,t17,t1 );	FNMA231(1.0,t18,t2 );\
	re =t9 ;	im =t10;	/* mpy by E^-4 = -I is inlined here... */\
	FNMA231(1.0,t26,t9 );	 FMA132(1.0,t26,re );	 FMA231(1.0,t25,t10);	FNMA132(1.0,t25,im );\
\
	__B0r = rt;	__B0i = it;	rt =t1;					re =t26;	im =t9 ;\
	 FMA231(t2 ,__t8,t1 );	FNMA231(rt ,__t8,t2 );	 FMA231(__t4,t25,t26);	FNMA231(__t4,re ,t25);	 FMA231(__tC,t10,t9 );	FNMA231(__tC,im ,t10);\
	t1 *= __c8;				t2 *= __c8;				t26 *= __c4;	t25 *= __c4;	t9  *= __cC;	t10 *= __cC;\
	__B8r = t1 ;	__B8i = t2 ;					__B4r = t26;	__B4i = t25;	__BCr = t9 ;	__BCi = t10;\
\
	/*...Block 3: t5,13,21,29 */\
	rt = t5;	it = t6;	/* twiddle mpy by E^4 =-I */\
	 FMA231(1.0,t14,t5 );	FNMA132(1.0,t14,rt );	FNMA231(1.0,t13,t6 );	 FMA132(1.0,t13,it );\
	re = t22;	im = t30;\
	FNMA231(1.0,t21,t22);	 FMA132(1.0,t21,re );	 FMA231(1.0,t29,t30);	 FMS132(1.0,t29,im );\
\
	rt = t21;	it = t22;\
	FNMA231(1.0,t29,t21);	 FMA132(1.0,t29,rt );	FNMA231(1.0,t30,t22);	 FMA132(1.0,t30,it );\
\
	rt = t5;	it = t6;\
	FNMA231(ISRT2,t21,t5 );	 FMA132(ISRT2,t21,rt );	FNMA231(ISRT2,t22,t6 );	 FMA132(ISRT2,t22,it );\
	re =t13;	im =t14;\
	 FMA231(ISRT2,t29,t13);	FNMA132(ISRT2,t29,re );	FNMA231(ISRT2,t30,t14);	 FMA132(ISRT2,t30,im );\
\
	rt =t21;	it = t5;\
	 FMA231(__t2,t22,t21);	FNMA231(__t2,rt ,t22);	 FMA231(__tA,t6 ,t5 );	FNMA231(__tA,it ,t6 );\
	re =t30;	im =t14;\
	 FMA231(__t6,t29,t30);	FNMA231(__t6,re ,t29);	 FMA231(__tE,t13,t14);	FNMA231(__tE,im ,t13);\
	t21 *= __c2;	t22 *= __c2;	t5  *= __cA;	t6  *= __cA;\
	t30 *= __c6;	t29 *= __c6;	t14 *= __cE;	t13 *= __cE;\
	__B2r = t21;	__B2i = t22;	__BAr = t5 ;	__BAi = t6 ;\
	__B6r = t30;	__B6i = t29;	__BEr = t14;	__BEi = t13;\
\
	/*...Block 2: t3,11,19,27 */\
	rt =t11;\
	 FMA231(1.0,t12,t11);	FNMA231(1.0,rt ,t12);\
	re = t19;	im = t27;\
	 FMA231(__sc,t20,t19);	FNMA231(__sc,re ,t20);	 FMA132(__sc,t27,t28);	 FMS132(__sc,t28,im );\
\
	rt =t3 ;	it =t4 ;\
	 FMA231(ISRT2,t11,t3 );	FNMA132(ISRT2,t11,rt );	 FMA231(ISRT2,t12,t4 );	FNMA132(ISRT2,t12,it );\
	t27 *= __c;	t28 *= __c;\
\
	re = t27;	im = t28;\
	 FMS231(__c,t19,t27);	 FMA132(__c,t19,re );	 FMS231(__c,t20,t28);	 FMA132(__c,t20,im );\
\
	rt = t3; it = t4;\
	 FMS132(1.0,t3 ,t19);	 FMS132(1.0,t4 ,t20);	 FMA231(1.0,rt ,t19);	 FMA231(1.0,it ,t20);\
	re = t11; im = t12;\
	 FMS132(1.0,t11,t28);	 FMA132(1.0,t12,t27);	 FMA231(1.0,re ,t28);	 FMS231(1.0,im ,t27);\
\
	rt = t19;	it = t3;\
	 FMA231(__t1,t20,t19);	FNMA231(__t1,rt ,t20);	 FMA231(__t9,t4 ,t3 );	FNMA231(__t9,it ,t4 );\
	re = t28;	im = t11;\
	 FMA231(__t5,t27,t28);	FNMA231(__t5,re ,t27);	 FMA231(__tD,t12,t11);	FNMA231(__tD,im ,t12);\
\
	t19 *= __c1;	t20 *= __c1;	t3  *= __c9;	t4  *= __c9;\
	t28 *= __c5;	t27 *= __c5;	t11 *= __cD;	t12 *= __cD;\
\
	__B1r = t19;	__B1i = t20;	__B9r = t3 ;	__B9i = t4 ;\
	__B5r = t28;	__B5i = t27;	__BDr = t11;	__BDi = t12;\
\
	/*...Block 4: t7,15,23,31 */\
	rt = t16;\
	 FMA231(1.0,t15,t16);	 FMS132(1.0,t15,rt );\
	re = t23;	im = t31;\
	FMA231(__sc,t32,t31);	FNMA231(__sc,im ,t32);	FMA132(__sc,t23,t24);	 FMS132(__sc,t24,re );\
\
	rt =t7 ;	it =t8 ;\
	FNMA231(ISRT2,t15,t7 );	 FMA132(ISRT2,t15,rt );	FNMA231(ISRT2,t16,t8 );	 FMA132(ISRT2,t16,it );\
	t31 *= __c;	t32 *= __c;	/* Option A */\
							/* Option B: Defer *= __c on t23/24,31/32 until they get add/subbed with t15/16,7/8 ***/\
\
	re = t31;	im = t32;\
	 FMA231(__c,t23,t31);	 FMS132(__c,t23,re );	 FMA231(__c,t24,t32);	 FMS132(__c,t24,im );/* A */\
/*	 FMA231(1.0,t23,t31);	 FMS132(1.0,t23,re );	 FMA231(1.0,t24,t32);	 FMS132(1.0,t24,im );** B */\
\
	rt = t15; it = t16;\
	 FMS132(1.0,t15,t32);	 FMA132(1.0,t16,t31);	 FMA231(1.0,rt ,t32);	 FMS231(1.0,it ,t31);/* A */\
/*	FNMA231(__c,t32,t15);	 FMA231(__c,t31,t16);	 FMA132(__c,t32,rt );	FNMA132(__c,t31,it );** B */\
	re = t7; im = t8;\
	 FMS132(1.0,t7 ,t23 );	 FMS132(1.0,t8 ,t24 );	 FMA231(1.0,re ,t23 );	 FMA231(1.0,im ,t24);/* A */\
/*	FNMA231(__c,t23 ,t7 );	FNMA231(__c,t24 ,t8 );	 FMA132(__c,t23 ,re );	 FMA132(__c,t24,im );** B */\
\
	rt = t23;	it = t7;\
	 FMA231(__t3,t24,t23);	FNMA231(__t3,rt ,t24);	 FMA231(__tB,t8 ,t7 );	FNMA231(__tB,it ,t8 );\
	re = t32;	im = t15;\
	 FMA231(__t7,t31,t32);	FNMA231(__t7,re ,t31);	 FMA231(__tF,t16,t15);	FNMA231(__tF,im ,t16);\
\
	t23 *= __c3;	t24 *= __c3;	t7  *= __cB;	t8  *= __cB;\
	t32 *= __c7;	t31 *= __c7;	t15 *= __cF;	t16 *= __cF;\
\
	__B3r = t23;	__B3i = t24;	__BBr = t7 ;	__BBi = t8 ;\
	__B7r = t32;	__B7i = t31;	__BFr = t15;	__BFi = t16;\
}


/************** RADIX-32 DIF/DIT: *****************************/
/* Totals: 376 ADD, 88 MUL	*/
/* Because of the way the original (non-macro-ized) code was written, it's convenient to index the A-inputs in terms
of 4 length-8 blocks with octal indices, and the B-outputs in terms of 2 length-16 blocks with hexadecimal indices.
MSVC allows a maximum of 'only' 127 macro args (which is probably a good thing), so unlike the smaller-radix DFT
macros which use actual array-indexed terms as args, here we use pointers to the real part of each complex arg:
*/
#define RADIX_32_DIF(\
	__A,__idx,	/*  Inputs: Base address plus 32 (index) offsets */\
	__B,__odx	/* Outputs: Base address plus 32 (index) offsets */\
)\
{\
	double __rt,__it\
		,__t00,__t01,__t02,__t03,__t04,__t05,__t06,__t07,__t08,__t09,__t0A,__t0B,__t0C,__t0D,__t0E,__t0F\
		,__t10,__t11,__t12,__t13,__t14,__t15,__t16,__t17,__t18,__t19,__t1A,__t1B,__t1C,__t1D,__t1E,__t1F\
		,__t20,__t21,__t22,__t23,__t24,__t25,__t26,__t27,__t28,__t29,__t2A,__t2B,__t2C,__t2D,__t2E,__t2F\
		,__t30,__t31,__t32,__t33,__t34,__t35,__t36,__t37,__t38,__t39,__t3A,__t3B,__t3C,__t3D,__t3E,__t3F;\
	double *Aim = __A + RE_IM_STRIDE, *Bim = __B + RE_IM_STRIDE;\
	/* Gather the needed data (32 64-bit complex, i.e. 64 64-bit reals) and do the first set of four length-8 transforms...	*/\
	/* Each complex radix-8 subtransform needs 52 ADD, 4 MUL (not counting load/store/address-compute): */\
	/*...Block 1: */\
		__t00=*(__A+__idx[0x00]);			__t01=*(Aim+__idx[0x00]);\
		__rt =*(__A+__idx[0x10]);			__it =*(Aim+__idx[0x10]);\
		__t02=__t00-__rt;					__t03=__t01-__it;\
		__t00=__t00+__rt;					__t01=__t01+__it;\
\
		__t04=*(__A+__idx[0x08]);			__t05=*(Aim+__idx[0x08]);\
		__rt =*(__A+__idx[0x18]);			__it =*(Aim+__idx[0x18]);\
		__t06=__t04-__rt;					__t07=__t05-__it;\
		__t04=__t04+__rt;					__t05=__t05+__it;\
\
		__rt =__t04;						__it =__t05;\
		__t04=__t00-__rt;					__t05=__t01-__it;\
		__t00=__t00+__rt;					__t01=__t01+__it;\
\
		__rt =__t06;						__it =__t07;\
		__t06=__t02+__it;					__t07=__t03-__rt;\
		__t02=__t02-__it;					__t03=__t03+__rt;\
\
		__t08=*(__A+__idx[0x04]);			__t09=*(Aim+__idx[0x04]);\
		__rt =*(__A+__idx[0x14]);			__it =*(Aim+__idx[0x14]);\
		__t0A=__t08-__rt;					__t0B=__t09-__it;\
		__t08=__t08+__rt;					__t09=__t09+__it;\
\
		__t0C=*(__A+__idx[0x0C]);			__t0D=*(Aim+__idx[0x0C]);\
		__rt =*(__A+__idx[0x1C]);			__it =*(Aim+__idx[0x1C]);\
		__t0E=__t0C-__rt;					__t0F=__t0D-__it;\
		__t0C=__t0C+__rt;					__t0D=__t0D+__it;\
\
		__rt =__t0C;						__it =__t0D;\
		__t0C=__t08-__rt;					__t0D=__t09-__it;\
		__t08=__t08+__rt;					__t09=__t09+__it;\
\
		__rt =__t0E;						__it =__t0F;\
		__t0E=__t0A+__it;					__t0F=__t0B-__rt;\
		__t0A=__t0A-__it;					__t0B=__t0B+__rt;\
\
		__rt =__t08;						__it =__t09;\
		__t08=__t00-__rt;					__t09=__t01-__it;\
		__t00=__t00+__rt;					__t01=__t01+__it;\
\
		__rt =__t0C;						__it =__t0D;\
		__t0C=__t04+__it;					__t0D=__t05-__rt;\
		__t04=__t04-__it;					__t05=__t05+__rt;\
\
		__rt =(__t0A-__t0B)*ISRT2;			__it =(__t0A+__t0B)*ISRT2;\
		__t0A=__t02-__rt;					__t0B=__t03-__it;\
		__t02=__t02+__rt;					__t03=__t03+__it;\
\
		__rt =(__t0E+__t0F)*ISRT2;			__it =(__t0F-__t0E)*ISRT2;\
		__t0E=__t06+__rt;					__t0F=__t07+__it;\
		__t06=__t06-__rt;					__t07=__t07-__it;\
\
	/*...Block 2:;*/\
		__t10=*(__A+__idx[0x02]);			__t11=*(Aim+__idx[0x02]);\
		__rt =*(__A+__idx[0x12]);			__it =*(Aim+__idx[0x12]);\
		__t12=__t10-__rt;					__t13=__t11-__it;\
		__t10=__t10+__rt;					__t11=__t11+__it;\
\
		__t14=*(__A+__idx[0x0A]);			__t15=*(Aim+__idx[0x0A]);\
		__rt =*(__A+__idx[0x1A]);			__it =*(Aim+__idx[0x1A]);\
		__t16=__t14-__rt;					__t17=__t15-__it;\
		__t14=__t14+__rt;					__t15=__t15+__it;\
\
		__rt =__t14;						__it =__t15;\
		__t14=__t10-__rt;					__t15=__t11-__it;\
		__t10=__t10+__rt;					__t11=__t11+__it;\
\
		__rt =__t16;						__it =__t17;\
		__t16=__t12+__it;					__t17=__t13-__rt;\
		__t12=__t12-__it;					__t13=__t13+__rt;\
\
		__t18=*(__A+__idx[0x06]);			__t19=*(Aim+__idx[0x06]);\
		__rt =*(__A+__idx[0x16]);			__it =*(Aim+__idx[0x16]);\
		__t1A=__t18-__rt;					__t1B=__t19-__it;\
		__t18=__t18+__rt;					__t19=__t19+__it;\
\
		__t1C=*(__A+__idx[0x0E]);			__t1D=*(Aim+__idx[0x0E]);\
		__rt =*(__A+__idx[0x1E]);			__it =*(Aim+__idx[0x1E]);\
		__t1E=__t1C-__rt;					__t1F=__t1D-__it;\
		__t1C=__t1C+__rt;					__t1D=__t1D+__it;\
\
		__rt =__t1C;						__it =__t1D;\
		__t1C=__t18-__rt;					__t1D=__t19-__it;\
		__t18=__t18+__rt;					__t19=__t19+__it;\
\
		__rt =__t1E;						__it =__t1F;\
		__t1E=__t1A+__it;					__t1F=__t1B-__rt;\
		__t1A=__t1A-__it;					__t1B=__t1B+__rt;\
\
		__rt =__t18;						__it =__t19;\
		__t18=__t10-__rt;					__t19=__t11-__it;\
		__t10=__t10+__rt;					__t11=__t11+__it;\
\
		__rt =__t1C;						__it =__t1D;\
		__t1C=__t14+__it;					__t1D=__t15-__rt;\
		__t14=__t14-__it;					__t15=__t15+__rt;\
\
		__rt =(__t1A-__t1B)*ISRT2;			__it =(__t1A+__t1B)*ISRT2;\
		__t1A=__t12-__rt;					__t1B=__t13-__it;\
		__t12=__t12+__rt;					__t13=__t13+__it;\
\
		__rt =(__t1E+__t1F)*ISRT2;			__it =(__t1F-__t1E)*ISRT2;\
		__t1E=__t16+__rt;					__t1F=__t17+__it;\
		__t16=__t16-__rt;					__t17=__t17-__it;\
\
	/*...Block 3: */\
		__t20=*(__A+__idx[0x01]);			__t21=*(Aim+__idx[0x01]);\
		__rt =*(__A+__idx[0x11]);			__it =*(Aim+__idx[0x11]);\
		__t22=__t20-__rt;					__t23=__t21-__it;\
		__t20=__t20+__rt;					__t21=__t21+__it;\
\
		__t24=*(__A+__idx[0x09]);			__t25=*(Aim+__idx[0x09]);\
		__rt =*(__A+__idx[0x19]);			__it =*(Aim+__idx[0x19]);\
		__t26=__t24-__rt;					__t27=__t25-__it;\
		__t24=__t24+__rt;					__t25=__t25+__it;\
\
		__rt =__t24;						__it =__t25;\
		__t24=__t20-__rt;					__t25=__t21-__it;\
		__t20=__t20+__rt;					__t21=__t21+__it;\
\
		__rt =__t26;						__it =__t27;\
		__t26=__t22+__it;					__t27=__t23-__rt;\
		__t22=__t22-__it;					__t23=__t23+__rt;\
\
		__t28=*(__A+__idx[0x05]);			__t29=*(Aim+__idx[0x05]);\
		__rt =*(__A+__idx[0x15]);			__it =*(Aim+__idx[0x15]);\
		__t2A=__t28-__rt;					__t2B=__t29-__it;\
		__t28=__t28+__rt;					__t29=__t29+__it;\
\
		__t2C=*(__A+__idx[0x0D]);			__t2D=*(Aim+__idx[0x0D]);\
		__rt =*(__A+__idx[0x1D]);			__it =*(Aim+__idx[0x1D]);\
		__t2E=__t2C-__rt;					__t2F=__t2D-__it;\
		__t2C=__t2C+__rt;					__t2D=__t2D+__it;\
\
		__rt =__t2C;						__it =__t2D;\
		__t2C=__t28-__rt;					__t2D=__t29-__it;\
		__t28=__t28+__rt;					__t29=__t29+__it;\
\
		__rt =__t2E;						__it =__t2F;\
		__t2E=__t2A+__it;					__t2F=__t2B-__rt;\
		__t2A=__t2A-__it;					__t2B=__t2B+__rt;\
\
		__rt =__t28;						__it =__t29;\
		__t28=__t20-__rt;					__t29=__t21-__it;\
		__t20=__t20+__rt;					__t21=__t21+__it;\
\
		__rt =__t2C;						__it =__t2D;\
		__t2C=__t24+__it;					__t2D=__t25-__rt;\
		__t24=__t24-__it;					__t25=__t25+__rt;\
\
		__rt =(__t2A-__t2B)*ISRT2;			__it =(__t2A+__t2B)*ISRT2;\
		__t2A=__t22-__rt;					__t2B=__t23-__it;\
		__t22=__t22+__rt;					__t23=__t23+__it;\
\
		__rt =(__t2E+__t2F)*ISRT2;			__it =(__t2F-__t2E)*ISRT2;\
		__t2E=__t26+__rt;					__t2F=__t27+__it;\
		__t26=__t26-__rt;					__t27=__t27-__it;\
\
	/*...Block 4: */\
		__t30=*(__A+__idx[0x03]);			__t31=*(Aim+__idx[0x03]);\
		__rt =*(__A+__idx[0x13]);			__it =*(Aim+__idx[0x13]);\
		__t32=__t30-__rt;					__t33=__t31-__it;\
		__t30=__t30+__rt;					__t31=__t31+__it;\
\
		__t34=*(__A+__idx[0x0B]);			__t35=*(Aim+__idx[0x0B]);\
		__rt =*(__A+__idx[0x1B]);			__it =*(Aim+__idx[0x1B]);\
		__t36=__t34-__rt;					__t37=__t35-__it;\
		__t34=__t34+__rt;					__t35=__t35+__it;\
\
		__rt =__t34;						__it =__t35;\
		__t34=__t30-__rt;					__t35=__t31-__it;\
		__t30=__t30+__rt;					__t31=__t31+__it;\
\
		__rt =__t36;						__it =__t37;\
		__t36=__t32+__it;					__t37=__t33-__rt;\
		__t32=__t32-__it;					__t33=__t33+__rt;\
\
		__t38=*(__A+__idx[0x07]);			__t39=*(Aim+__idx[0x07]);\
		__rt =*(__A+__idx[0x17]);			__it =*(Aim+__idx[0x17]);\
		__t3A=__t38-__rt;					__t3B=__t39-__it;\
		__t38=__t38+__rt;					__t39=__t39+__it;\
\
		__t3C=*(__A+__idx[0x0F]);			__t3D=*(Aim+__idx[0x0F]);\
		__rt =*(__A+__idx[0x1F]);			__it =*(Aim+__idx[0x1F]);\
		__t3E=__t3C-__rt;					__t3F=__t3D-__it;\
		__t3C=__t3C+__rt;					__t3D=__t3D+__it;\
\
		__rt =__t3C;						__it =__t3D;\
		__t3C=__t38-__rt;					__t3D=__t39-__it;\
		__t38=__t38+__rt;					__t39=__t39+__it;\
\
		__rt =__t3E;						__it =__t3F;\
		__t3E=__t3A+__it;					__t3F=__t3B-__rt;\
		__t3A=__t3A-__it;					__t3B=__t3B+__rt;\
\
		__rt =__t38;						__it =__t39;\
		__t38=__t30-__rt;					__t39=__t31-__it;\
		__t30=__t30+__rt;					__t31=__t31+__it;\
\
		__rt =__t3C;						__it =__t3D;\
		__t3C=__t34+__it;					__t3D=__t35-__rt;\
		__t34=__t34-__it;					__t35=__t35+__rt;\
\
		__rt =(__t3A-__t3B)*ISRT2;			__it =(__t3A+__t3B)*ISRT2;\
		__t3A=__t32-__rt;					__t3B=__t33-__it;\
		__t32=__t32+__rt;					__t33=__t33+__it;\
\
		__rt =(__t3E+__t3F)*ISRT2;			__it =(__t3F-__t3E)*ISRT2;\
		__t3E=__t36+__rt;					__t3F=__t37+__it;\
		__t36=__t36-__rt;					__t37=__t37-__it;\
\
	/*...and now do eight radix-4 transforms, including the internal twiddle factors: */\
	/* Totals for the eight radix-4: 168 ADD, 72 MUL: */\
	/*...Block 1: __t00,__t10,__t20,__t30	*/\
		__rt =__t10;	__t10=__t00-__rt;	__t00=__t00+__rt;\
		__it =__t11;	__t11=__t01-__it;	__t01=__t01+__it;\
\
		__rt =__t30;	__t30=__t20-__rt;	__t20=__t20+__rt;\
		__it =__t31;	__t31=__t21-__it;	__t21=__t21+__it;\
	/* 16 ADD, 0 MUL: */\
		*(__B+__odx[0x00])=__t00+__t20;		*(Bim+__odx[0x00])=__t01+__t21;\
		*(__B+__odx[0x01])=__t00-__t20;		*(Bim+__odx[0x01])=__t01-__t21;\
		*(__B+__odx[0x02])=__t10-__t31;		*(Bim+__odx[0x02])=__t11+__t30;\
		*(__B+__odx[0x03])=__t10+__t31;		*(Bim+__odx[0x03])=__t11-__t30;\
\
	/*...Block 5: __t08,__t18,__t28,__t38	*/\
		__rt =__t18;\
		__t18=__t08+__t19;					__t08=__t08-__t19;\
		__t19=__t09-__rt;					__t09=__t09+__rt;\
\
		__rt =(__t28-__t29)*ISRT2;			__t29=(__t28+__t29)*ISRT2;	__t28=__rt;\
		__rt =(__t39+__t38)*ISRT2;			__it =(__t39-__t38)*ISRT2;\
		__t38=__t28+__rt;					__t28=__t28-__rt;\
		__t39=__t29+__it;					__t29=__t29-__it;\
	/* 20 ADD, 4 MUL: */\
		*(__B+__odx[0x04])=__t08+__t28;		*(Bim+__odx[0x04])=__t09+__t29;\
		*(__B+__odx[0x05])=__t08-__t28;		*(Bim+__odx[0x05])=__t09-__t29;\
		*(__B+__odx[0x06])=__t18-__t39;		*(Bim+__odx[0x06])=__t19+__t38;\
		*(__B+__odx[0x07])=__t18+__t39;		*(Bim+__odx[0x07])=__t19-__t38;\
\
	/*...Block 3: __t04,__t14,__t24,__t34	*/\
		__rt =(__t14-__t15)*ISRT2;			__it =(__t14+__t15)*ISRT2;\
		__t14=__t04-__rt;					__t04=__t04+__rt;\
		__t15=__t05-__it;					__t05=__t05+__it;\
\
		__rt =__t24*c - __t25*s;			__t25=__t25*c + __t24*s;	__t24=__rt;\
		__rt =__t34*s - __t35*c;			__it =__t35*s + __t34*c;\
		__t34=__t24-__rt;					__t24=__t24+__rt;\
		__t35=__t25-__it;					__t25=__t25+__it;\
	/* 22 ADD, 10 MUL: */\
		*(__B+__odx[0x08])=__t04+__t24;		*(Bim+__odx[0x08])=__t05+__t25;\
		*(__B+__odx[0x09])=__t04-__t24;		*(Bim+__odx[0x09])=__t05-__t25;\
		*(__B+__odx[0x0A])=__t14-__t35;		*(Bim+__odx[0x0A])=__t15+__t34;\
		*(__B+__odx[0x0B])=__t14+__t35;		*(Bim+__odx[0x0B])=__t15-__t34;\
\
	/*...Block 7: __t0C,__t1C,__t2C,__t3C	*/\
		__rt =(__t1D+__t1C)*ISRT2;			__it =(__t1D-__t1C)*ISRT2;\
		__t1C=__t0C+__rt;					__t0C=__t0C-__rt;\
		__t1D=__t0D+__it;					__t0D=__t0D-__it;\
\
		__rt =__t2C*s - __t2D*c;			__t2D=__t2D*s + __t2C*c;	__t2C=__rt;\
		__rt =__t3C*c - __t3D*s;			__it =__t3D*c + __t3C*s;\
		__t3C=__t2C+__rt;					__t2C=__t2C-__rt;\
		__t3D=__t2D+__it;					__t2D=__t2D-__it;\
	/* 22 ADD, 10 MUL: */\
		*(__B+__odx[0x0C])=__t0C+__t2C;		*(Bim+__odx[0x0C])=__t0D+__t2D;\
		*(__B+__odx[0x0D])=__t0C-__t2C;		*(Bim+__odx[0x0D])=__t0D-__t2D;\
		*(__B+__odx[0x0E])=__t1C-__t3D;		*(Bim+__odx[0x0E])=__t1D+__t3C;\
		*(__B+__odx[0x0F])=__t1C+__t3D;		*(Bim+__odx[0x0F])=__t1D-__t3C;\
\
	/*...Block 2: __t02,__t12,__t22,__t32	*/\
		__rt =__t12*c - __t13*s;			__it =__t13*c + __t12*s;\
		__t12=__t02-__rt;					__t02=__t02+__rt;\
		__t13=__t03-__it;					__t03=__t03+__it;\
\
		__rt =__t22*c32_1 - __t23*s32_1;	__t23=__t23*c32_1 + __t22*s32_1;	__t22=__rt;\
		__rt =__t32*c32_3 - __t33*s32_3;	__it =__t33*c32_3 + __t32*s32_3;\
		__t32=__t22-__rt;					__t22=__t22+__rt;\
		__t33=__t23-__it;					__t23=__t23+__it;\
	/* 22 ADD, 12 MUL: */\
		*(__B+__odx[0x10])=__t02+__t22;		*(Bim+__odx[0x10])=__t03+__t23;\
		*(__B+__odx[0x11])=__t02-__t22;		*(Bim+__odx[0x11])=__t03-__t23;\
		*(__B+__odx[0x12])=__t12-__t33;		*(Bim+__odx[0x12])=__t13+__t32;\
		*(__B+__odx[0x13])=__t12+__t33;		*(Bim+__odx[0x13])=__t13-__t32;\
\
	/*...Block 6: __t0A,__t1A,__t2A,__t3A	*/\
		__rt =__t1A*s + __t1B*c;			__it =__t1B*s - __t1A*c;\
		__t1A=__t0A+__rt;					__t0A=__t0A-__rt;\
		__t1B=__t0B+__it;					__t0B=__t0B-__it;\
\
		__rt =__t2A*s32_3 - __t2B*c32_3;	__t2B=__t2B*s32_3 + __t2A*c32_3;	__t2A=__rt;\
		__rt =__t3A*c32_1 + __t3B*s32_1;	__it =__t3B*c32_1 - __t3A*s32_1;\
		__t3A=__t2A+__rt;					__t2A=__t2A-__rt;\
		__t3B=__t2B+__it;					__t2B=__t2B-__it;\
	/* 22 ADD, 12 MUL: */\
		*(__B+__odx[0x14])=__t0A+__t2A;		*(Bim+__odx[0x14])=__t0B+__t2B;\
		*(__B+__odx[0x15])=__t0A-__t2A;		*(Bim+__odx[0x15])=__t0B-__t2B;\
		*(__B+__odx[0x16])=__t1A-__t3B;		*(Bim+__odx[0x16])=__t1B+__t3A;\
		*(__B+__odx[0x17])=__t1A+__t3B;		*(Bim+__odx[0x17])=__t1B-__t3A;\
\
	/*...Block 4: __t06,__t16,__t26,__t36	*/\
		__rt =__t16*s - __t17*c;			__it =__t17*s + __t16*c;\
		__t16=__t06-__rt;					__t06=__t06+__rt;\
		__t17=__t07-__it;					__t07=__t07+__it;\
\
		__rt =__t26*c32_3 - __t27*s32_3;	__t27=__t27*c32_3 + __t26*s32_3;	__t26=__rt;\
		__rt =__t36*s32_1 + __t37*c32_1;	__it =__t37*s32_1 - __t36*c32_1;\
		__t36=__t26+__rt;					__t26=__t26-__rt;\
		__t37=__t27+__it;					__t27=__t27-__it;\
	/* 22 ADD, 12 MUL: */\
		*(__B+__odx[0x18])=__t06+__t26;		*(Bim+__odx[0x18])=__t07+__t27;\
		*(__B+__odx[0x19])=__t06-__t26;		*(Bim+__odx[0x19])=__t07-__t27;\
		*(__B+__odx[0x1A])=__t16-__t37;		*(Bim+__odx[0x1A])=__t17+__t36;\
		*(__B+__odx[0x1B])=__t16+__t37;		*(Bim+__odx[0x1B])=__t17-__t36;\
\
	/*...Block 8: __t0E,__t1E,__t2E,__t3E	*/\
		__rt =__t1E*c + __t1F*s;			__it =__t1F*c - __t1E*s;\
		__t1E=__t0E+__rt;					__t0E=__t0E-__rt;\
		__t1F=__t0F+__it;					__t0F=__t0F-__it;\
\
		__rt =__t2E*s32_1 - __t2F*c32_1;	__t2F=__t2F*s32_1 + __t2E*c32_1;	__t2E=__rt;\
		__rt =__t3E*s32_3 - __t3F*c32_3;	__it =__t3F*s32_3 + __t3E*c32_3;\
		__t3E=__t2E+__rt;					__t2E=__t2E-__rt;\
		__t3F=__t2F+__it;					__t2F=__t2F-__it;\
	/* 22 ADD, 12 MUL: */\
		*(__B+__odx[0x1C])=__t0E+__t2E;		*(Bim+__odx[0x1C])=__t0F+__t2F;\
		*(__B+__odx[0x1D])=__t0E-__t2E;		*(Bim+__odx[0x1D])=__t0F-__t2F;\
		*(__B+__odx[0x1E])=__t1E-__t3F;		*(Bim+__odx[0x1E])=__t1F+__t3E;\
		*(__B+__odx[0x1F])=__t1E+__t3F;		*(Bim+__odx[0x1F])=__t1F-__t3E;\
}

#define RADIX_32_DIT(\
	__A,__idx,	/*  Inputs: Base address plus 32 (index) offsets */\
	__B,__odx	/* Outputs: Base address plus 32 (index) offsets */\
)\
{\
	double __rt,__it\
		,__t00,__t01,__t02,__t03,__t04,__t05,__t06,__t07,__t08,__t09,__t0A,__t0B,__t0C,__t0D,__t0E,__t0F\
		,__t10,__t11,__t12,__t13,__t14,__t15,__t16,__t17,__t18,__t19,__t1A,__t1B,__t1C,__t1D,__t1E,__t1F\
		,__t20,__t21,__t22,__t23,__t24,__t25,__t26,__t27,__t28,__t29,__t2A,__t2B,__t2C,__t2D,__t2E,__t2F\
		,__t30,__t31,__t32,__t33,__t34,__t35,__t36,__t37,__t38,__t39,__t3A,__t3B,__t3C,__t3D,__t3E,__t3F;\
	double *Aim = __A + RE_IM_STRIDE, *Bim = __B + RE_IM_STRIDE;\
	/* Gather the needed data (32 64-bit complex, i.e. 64 64-bit reals) and do the first set of four length-8 transforms...	*/\
	/* Each complex radix-8 subtransform needs 52 ADD, 4 MUL (not counting load/store/address-compute): */\
	/*...Block 1: */\
		__t00=*(__A+__idx[0x00]);			__t01=*(Aim+__idx[0x00]);\
		__rt =*(__A+__idx[0x01]);			__it =*(Aim+__idx[0x01]);\
		__t02=__t00-__rt;					__t03=__t01-__it;\
		__t00=__t00+__rt;					__t01=__t01+__it;\
\
		__t04=*(__A+__idx[0x02]);			__t05=*(Aim+__idx[0x02]);\
		__rt =*(__A+__idx[0x03]);			__it =*(Aim+__idx[0x03]);\
		__t06=__t04-__rt;					__t07=__t05-__it;\
		__t04=__t04+__rt;					__t05=__t05+__it;\
\
		__rt =__t04;						__it =__t05;\
		__t04=__t00-__rt;					__t05=__t01-__it;\
		__t00=__t00+__rt;					__t01=__t01+__it;\
\
		__rt =__t06;						__it =__t07;\
		__t06=__t02-__it;					__t07=__t03+__rt;\
		__t02=__t02+__it;					__t03=__t03-__rt;\
\
		__t08=*(__A+__idx[0x04]);			__t09=*(Aim+__idx[0x04]);\
		__rt =*(__A+__idx[0x05]);			__it =*(Aim+__idx[0x05]);\
		__t0A=__t08-__rt;					__t0B=__t09-__it;\
		__t08=__t08+__rt;					__t09=__t09+__it;\
\
		__t0C=*(__A+__idx[0x06]);			__t0D=*(Aim+__idx[0x06]);\
		__rt =*(__A+__idx[0x07]);			__it =*(Aim+__idx[0x07]);\
		__t0E=__t0C-__rt;					__t0F=__t0D-__it;\
		__t0C=__t0C+__rt;					__t0D=__t0D+__it;\
\
		__rt =__t0C;						__it =__t0D;\
		__t0C=__t08-__rt;					__t0D=__t09-__it;\
		__t08=__t08+__rt;					__t09=__t09+__it;\
\
		__rt =__t0E;						__it =__t0F;\
		__t0E=__t0A-__it;					__t0F=__t0B+__rt;\
		__t0A=__t0A+__it;					__t0B=__t0B-__rt;\
\
		__rt =__t08;						__it =__t09;\
		__t08=__t00-__rt;					__t09=__t01-__it;\
		__t00=__t00+__rt;					__t01=__t01+__it;\
\
		__rt =__t0C;						__it =__t0D;\
		__t0C=__t04-__it;					__t0D=__t05+__rt;\
		__t04=__t04+__it;					__t05=__t05-__rt;\
\
		__rt =(__t0A+__t0B)*ISRT2;			__it =(__t0A-__t0B)*ISRT2;\
		__t0A=__t02-__rt;					__t0B=__t03+__it;\
		__t02=__t02+__rt;					__t03=__t03-__it;\
\
		__rt =(__t0E-__t0F)*ISRT2;			__it =(__t0F+__t0E)*ISRT2;\
		__t0E=__t06+__rt;					__t0F=__t07+__it;\
		__t06=__t06-__rt;					__t07=__t07-__it;\
\
	/*...Block 2:;*/\
		__t10=*(__A+__idx[0x08]);			__t11=*(Aim+__idx[0x08]);\
		__rt =*(__A+__idx[0x09]);			__it =*(Aim+__idx[0x09]);\
		__t12=__t10-__rt;					__t13=__t11-__it;\
		__t10=__t10+__rt;					__t11=__t11+__it;\
\
		__t14=*(__A+__idx[0x0A]);			__t15=*(Aim+__idx[0x0A]);\
		__rt =*(__A+__idx[0x0B]);			__it =*(Aim+__idx[0x0B]);\
		__t16=__t14-__rt;					__t17=__t15-__it;\
		__t14=__t14+__rt;					__t15=__t15+__it;\
\
		__rt =__t14;						__it =__t15;\
		__t14=__t10-__rt;					__t15=__t11-__it;\
		__t10=__t10+__rt;					__t11=__t11+__it;\
\
		__rt =__t16;						__it =__t17;\
		__t16=__t12-__it;					__t17=__t13+__rt;\
		__t12=__t12+__it;					__t13=__t13-__rt;\
\
		__t18=*(__A+__idx[0x0C]);			__t19=*(Aim+__idx[0x0C]);\
		__rt =*(__A+__idx[0x0D]);			__it =*(Aim+__idx[0x0D]);\
		__t1A=__t18-__rt;					__t1B=__t19-__it;\
		__t18=__t18+__rt;					__t19=__t19+__it;\
\
		__t1C=*(__A+__idx[0x0E]);			__t1D=*(Aim+__idx[0x0E]);\
		__rt =*(__A+__idx[0x0F]);			__it =*(Aim+__idx[0x0F]);\
		__t1E=__t1C-__rt;					__t1F=__t1D-__it;\
		__t1C=__t1C+__rt;					__t1D=__t1D+__it;\
\
		__rt =__t1C;						__it =__t1D;\
		__t1C=__t18-__rt;					__t1D=__t19-__it;\
		__t18=__t18+__rt;					__t19=__t19+__it;\
\
		__rt =__t1E;						__it =__t1F;\
		__t1E=__t1A-__it;					__t1F=__t1B+__rt;\
		__t1A=__t1A+__it;					__t1B=__t1B-__rt;\
\
		__rt =__t18;						__it =__t19;\
		__t18=__t10-__rt;					__t19=__t11-__it;\
		__t10=__t10+__rt;					__t11=__t11+__it;\
\
		__rt =__t1C;						__it =__t1D;\
		__t1C=__t14-__it;					__t1D=__t15+__rt;\
		__t14=__t14+__it;					__t15=__t15-__rt;\
\
		__rt =(__t1A+__t1B)*ISRT2;			__it =(__t1A-__t1B)*ISRT2;\
		__t1A=__t12-__rt;					__t1B=__t13+__it;\
		__t12=__t12+__rt;					__t13=__t13-__it;\
\
		__rt =(__t1E-__t1F)*ISRT2;			__it =(__t1F+__t1E)*ISRT2;\
		__t1E=__t16+__rt;					__t1F=__t17+__it;\
		__t16=__t16-__rt;					__t17=__t17-__it;\
\
	/*...Block 3: */\
		__t20=*(__A+__idx[0x10]);			__t21=*(Aim+__idx[0x10]);\
		__rt =*(__A+__idx[0x11]);			__it =*(Aim+__idx[0x11]);\
		__t22=__t20-__rt;					__t23=__t21-__it;\
		__t20=__t20+__rt;					__t21=__t21+__it;\
\
		__t24=*(__A+__idx[0x12]);			__t25=*(Aim+__idx[0x12]);\
		__rt =*(__A+__idx[0x13]);			__it =*(Aim+__idx[0x13]);\
		__t26=__t24-__rt;					__t27=__t25-__it;\
		__t24=__t24+__rt;					__t25=__t25+__it;\
\
		__rt =__t24;						__it =__t25;\
		__t24=__t20-__rt;					__t25=__t21-__it;\
		__t20=__t20+__rt;					__t21=__t21+__it;\
\
		__rt =__t26;						__it =__t27;\
		__t26=__t22-__it;					__t27=__t23+__rt;\
		__t22=__t22+__it;					__t23=__t23-__rt;\
\
		__t28=*(__A+__idx[0x14]);			__t29=*(Aim+__idx[0x14]);\
		__rt =*(__A+__idx[0x15]);			__it =*(Aim+__idx[0x15]);\
		__t2A=__t28-__rt;					__t2B=__t29-__it;\
		__t28=__t28+__rt;					__t29=__t29+__it;\
\
		__t2C=*(__A+__idx[0x16]);			__t2D=*(Aim+__idx[0x16]);\
		__rt =*(__A+__idx[0x17]);			__it =*(Aim+__idx[0x17]);\
		__t2E=__t2C-__rt;					__t2F=__t2D-__it;\
		__t2C=__t2C+__rt;					__t2D=__t2D+__it;\
\
		__rt =__t2C;						__it =__t2D;\
		__t2C=__t28-__rt;					__t2D=__t29-__it;\
		__t28=__t28+__rt;					__t29=__t29+__it;\
\
		__rt =__t2E;						__it =__t2F;\
		__t2E=__t2A-__it;					__t2F=__t2B+__rt;\
		__t2A=__t2A+__it;					__t2B=__t2B-__rt;\
\
		__rt =__t28;						__it =__t29;\
		__t28=__t20-__rt;					__t29=__t21-__it;\
		__t20=__t20+__rt;					__t21=__t21+__it;\
\
		__rt =__t2C;						__it =__t2D;\
		__t2C=__t24-__it;					__t2D=__t25+__rt;\
		__t24=__t24+__it;					__t25=__t25-__rt;\
\
		__rt =(__t2A+__t2B)*ISRT2;			__it =(__t2A-__t2B)*ISRT2;\
		__t2A=__t22-__rt;					__t2B=__t23+__it;\
		__t22=__t22+__rt;					__t23=__t23-__it;\
\
		__rt =(__t2E-__t2F)*ISRT2;			__it =(__t2F+__t2E)*ISRT2;\
		__t2E=__t26+__rt;					__t2F=__t27+__it;\
		__t26=__t26-__rt;					__t27=__t27-__it;\
\
	/*...Block 4: */\
		__t30=*(__A+__idx[0x18]);			__t31=*(Aim+__idx[0x18]);\
		__rt =*(__A+__idx[0x19]);			__it =*(Aim+__idx[0x19]);\
		__t32=__t30-__rt;					__t33=__t31-__it;\
		__t30=__t30+__rt;					__t31=__t31+__it;\
\
		__t34=*(__A+__idx[0x1A]);			__t35=*(Aim+__idx[0x1A]);\
		__rt =*(__A+__idx[0x1B]);			__it =*(Aim+__idx[0x1B]);\
		__t36=__t34-__rt;					__t37=__t35-__it;\
		__t34=__t34+__rt;					__t35=__t35+__it;\
\
		__rt =__t34;						__it =__t35;\
		__t34=__t30-__rt;					__t35=__t31-__it;\
		__t30=__t30+__rt;					__t31=__t31+__it;\
\
		__rt =__t36;						__it =__t37;\
		__t36=__t32-__it;					__t37=__t33+__rt;\
		__t32=__t32+__it;					__t33=__t33-__rt;\
\
		__t38=*(__A+__idx[0x1C]);			__t39=*(Aim+__idx[0x1C]);\
		__rt =*(__A+__idx[0x1D]);			__it =*(Aim+__idx[0x1D]);\
		__t3A=__t38-__rt;					__t3B=__t39-__it;\
		__t38=__t38+__rt;					__t39=__t39+__it;\
\
		__t3C=*(__A+__idx[0x1E]);			__t3D=*(Aim+__idx[0x1E]);\
		__rt =*(__A+__idx[0x1F]);			__it =*(Aim+__idx[0x1F]);\
		__t3E=__t3C-__rt;					__t3F=__t3D-__it;\
		__t3C=__t3C+__rt;					__t3D=__t3D+__it;\
\
		__rt =__t3C;						__it =__t3D;\
		__t3C=__t38-__rt;					__t3D=__t39-__it;\
		__t38=__t38+__rt;					__t39=__t39+__it;\
\
		__rt =__t3E;						__it =__t3F;\
		__t3E=__t3A-__it;					__t3F=__t3B+__rt;\
		__t3A=__t3A+__it;					__t3B=__t3B-__rt;\
\
		__rt =__t38;						__it =__t39;\
		__t38=__t30-__rt;					__t39=__t31-__it;\
		__t30=__t30+__rt;					__t31=__t31+__it;\
\
		__rt =__t3C;						__it =__t3D;\
		__t3C=__t34-__it;					__t3D=__t35+__rt;\
		__t34=__t34+__it;					__t35=__t35-__rt;\
\
		__rt =(__t3A+__t3B)*ISRT2;			__it =(__t3A-__t3B)*ISRT2;\
		__t3A=__t32-__rt;					__t3B=__t33+__it;\
		__t32=__t32+__rt;					__t33=__t33-__it;\
\
		__rt =(__t3E-__t3F)*ISRT2;			__it =(__t3F+__t3E)*ISRT2;\
		__t3E=__t36+__rt;					__t3F=__t37+__it;\
		__t36=__t36-__rt;					__t37=__t37-__it;\
\
	/*...and now do eight radix-4 transforms, including the internal twiddle factors: */\
	/* Totals for the eight radix-4: 168 ADD, 72 MUL: */\
	/*...Block 1: __t00,__t10,__t20,__t30	*/\
		__rt =__t10;	__t10=__t00-__rt;	__t00=__t00+__rt;\
		__it =__t11;	__t11=__t01-__it;	__t01=__t01+__it;\
\
		__rt =__t30;	__t30=__t20-__rt;	__t20=__t20+__rt;\
		__it =__t31;	__t31=__t21-__it;	__t21=__t21+__it;\
	/* 16 ADD, 0 MUL: */\
		*(__B+__odx[0x00])=__t00+__t20;			*(Bim+__odx[0x00])=__t01+__t21;\
		*(__B+__odx[0x10])=__t00-__t20;			*(Bim+__odx[0x10])=__t01-__t21;\
		*(__B+__odx[0x08])=__t10+__t31;			*(Bim+__odx[0x08])=__t11-__t30;\
		*(__B+__odx[0x18])=__t10-__t31;			*(Bim+__odx[0x18])=__t11+__t30;\
\
	/*...Block 5: __t08,__t18,__t28,__t38	*/\
		__rt =__t18;\
		__t18=__t08-__t19;					__t08=__t08+__t19;\
		__t19=__t09+__rt;					__t09=__t09-__rt;\
\
		__rt =(__t29+__t28)*ISRT2;			__t29=(__t29-__t28)*ISRT2;	__t28=__rt;\
		__rt =(__t38-__t39)*ISRT2;			__it =(__t38+__t39)*ISRT2;\
		__t38=__t28+__rt;					__t28=__t28-__rt;\
		__t39=__t29+__it;					__t29=__t29-__it;\
	/* 20 ADD, 4 MUL: */\
		*(__B+__odx[0x04])=__t08+__t28;			*(Bim+__odx[0x04])=__t09+__t29;\
		*(__B+__odx[0x14])=__t08-__t28;			*(Bim+__odx[0x14])=__t09-__t29;\
		*(__B+__odx[0x0C])=__t18+__t39;			*(Bim+__odx[0x0C])=__t19-__t38;\
		*(__B+__odx[0x1C])=__t18-__t39;			*(Bim+__odx[0x1C])=__t19+__t38;\
\
	/*...Block 3: __t04,__t14,__t24,__t34	*/\
		__rt =(__t15+__t14)*ISRT2;			__it =(__t15-__t14)*ISRT2;\
		__t14=__t04-__rt;					__t04=__t04+__rt;\
		__t15=__t05-__it;					__t05=__t05+__it;\
\
		__rt =__t24*c + __t25*s;			__t25=__t25*c - __t24*s;	__t24=__rt;\
		__rt =__t34*s + __t35*c;			__it =__t35*s - __t34*c;\
		__t34=__t24-__rt;					__t24=__t24+__rt;\
		__t35=__t25-__it;					__t25=__t25+__it;\
	/* 22 ADD, 10 MUL: */\
		*(__B+__odx[0x02])=__t04+__t24;			*(Bim+__odx[0x02])=__t05+__t25;\
		*(__B+__odx[0x12])=__t04-__t24;			*(Bim+__odx[0x12])=__t05-__t25;\
		*(__B+__odx[0x0A])=__t14+__t35;			*(Bim+__odx[0x0A])=__t15-__t34;\
		*(__B+__odx[0x1A])=__t14-__t35;			*(Bim+__odx[0x1A])=__t15+__t34;\
\
	/*...Block 7: __t0C,__t1C,__t2C,__t3C	*/\
		__rt =(__t1C-__t1D)*ISRT2;			__it =(__t1C+__t1D)*ISRT2;\
		__t1C=__t0C+__rt;					__t0C=__t0C-__rt;\
		__t1D=__t0D+__it;					__t0D=__t0D-__it;\
\
		__rt =__t2C*s + __t2D*c;			__t2D=__t2D*s - __t2C*c;	__t2C=__rt;\
		__rt =__t3C*c + __t3D*s;			__it =__t3D*c - __t3C*s;\
		__t3C=__t2C+__rt;					__t2C=__t2C-__rt;\
		__t3D=__t2D+__it;					__t2D=__t2D-__it;\
	/* 22 ADD, 10 MUL: */\
		*(__B+__odx[0x06])=__t0C+__t2C;			*(Bim+__odx[0x06])=__t0D+__t2D;\
		*(__B+__odx[0x16])=__t0C-__t2C;			*(Bim+__odx[0x16])=__t0D-__t2D;\
		*(__B+__odx[0x0E])=__t1C+__t3D;			*(Bim+__odx[0x0E])=__t1D-__t3C;\
		*(__B+__odx[0x1E])=__t1C-__t3D;			*(Bim+__odx[0x1E])=__t1D+__t3C;\
\
	/*...Block 2: __t02,__t12,__t22,__t32	*/\
		__rt =__t12*c + __t13*s;			__it =__t13*c - __t12*s;\
		__t12=__t02-__rt;					__t02=__t02+__rt;\
		__t13=__t03-__it;					__t03=__t03+__it;\
\
		__rt =__t22*c32_1 + __t23*s32_1;	__t23=__t23*c32_1 - __t22*s32_1;	__t22=__rt;\
		__rt =__t32*c32_3 + __t33*s32_3;	__it =__t33*c32_3 - __t32*s32_3;\
		__t32=__t22-__rt;					__t22=__t22+__rt;\
		__t33=__t23-__it;					__t23=__t23+__it;\
	/* 22 ADD, 12 MUL: */\
		*(__B+__odx[0x01])=__t02+__t22;			*(Bim+__odx[0x01])=__t03+__t23;\
		*(__B+__odx[0x11])=__t02-__t22;			*(Bim+__odx[0x11])=__t03-__t23;\
		*(__B+__odx[0x09])=__t12+__t33;			*(Bim+__odx[0x09])=__t13-__t32;\
		*(__B+__odx[0x19])=__t12-__t33;			*(Bim+__odx[0x19])=__t13+__t32;\
\
	/*...Block 6: __t0A,__t1A,__t2A,__t3A	*/\
		__rt =__t1A*s - __t1B*c;			__it =__t1B*s + __t1A*c;\
		__t1A=__t0A+__rt;					__t0A=__t0A-__rt;\
		__t1B=__t0B+__it;					__t0B=__t0B-__it;\
\
		__rt =__t2A*s32_3 + __t2B*c32_3;	__t2B=__t2B*s32_3 - __t2A*c32_3;	__t2A=__rt;\
		__rt =__t3A*c32_1 - __t3B*s32_1;	__it =__t3B*c32_1 + __t3A*s32_1;\
		__t3A=__t2A+__rt;					__t2A=__t2A-__rt;\
		__t3B=__t2B+__it;					__t2B=__t2B-__it;\
	/* 22 ADD, 12 MUL: */\
		*(__B+__odx[0x05])=__t0A+__t2A;			*(Bim+__odx[0x05])=__t0B+__t2B;\
		*(__B+__odx[0x15])=__t0A-__t2A;			*(Bim+__odx[0x15])=__t0B-__t2B;\
		*(__B+__odx[0x0D])=__t1A+__t3B;			*(Bim+__odx[0x0D])=__t1B-__t3A;\
		*(__B+__odx[0x1D])=__t1A-__t3B;			*(Bim+__odx[0x1D])=__t1B+__t3A;\
\
	/*...Block 4: __t06,__t16,__t26,__t36	*/\
		__rt =__t16*s + __t17*c;			__it =__t17*s - __t16*c;\
		__t16=__t06-__rt;					__t06=__t06+__rt;\
		__t17=__t07-__it;					__t07=__t07+__it;\
\
		__rt =__t26*c32_3 + __t27*s32_3;	__t27=__t27*c32_3 - __t26*s32_3;	__t26=__rt;\
		__rt =__t36*s32_1 - __t37*c32_1;	__it =__t37*s32_1 + __t36*c32_1;\
		__t36=__t26+__rt;					__t26=__t26-__rt;\
		__t37=__t27+__it;					__t27=__t27-__it;\
	/* 22 ADD, 12 MUL: */\
		*(__B+__odx[0x03])=__t06+__t26;			*(Bim+__odx[0x03])=__t07+__t27;\
		*(__B+__odx[0x13])=__t06-__t26;			*(Bim+__odx[0x13])=__t07-__t27;\
		*(__B+__odx[0x0B])=__t16+__t37;			*(Bim+__odx[0x0B])=__t17-__t36;\
		*(__B+__odx[0x1B])=__t16-__t37;			*(Bim+__odx[0x1B])=__t17+__t36;\
\
	/*...Block 8: __t0E,__t1E,__t2E,__t3E	*/\
		__rt =__t1E*c - __t1F*s;			__it =__t1F*c + __t1E*s;\
		__t1E=__t0E+__rt;					__t0E=__t0E-__rt;\
		__t1F=__t0F+__it;					__t0F=__t0F-__it;\
\
		__rt =__t2E*s32_1 + __t2F*c32_1;	__t2F=__t2F*s32_1 - __t2E*c32_1;	__t2E=__rt;\
		__rt =__t3E*s32_3 + __t3F*c32_3;	__it =__t3F*s32_3 - __t3E*c32_3;\
		__t3E=__t2E+__rt;					__t2E=__t2E-__rt;\
		__t3F=__t2F+__it;					__t2F=__t2F-__it;\
	/* 22 ADD, 12 MUL: */\
		*(__B+__odx[0x07])=__t0E+__t2E;			*(Bim+__odx[0x07])=__t0F+__t2F;\
		*(__B+__odx[0x17])=__t0E-__t2E;			*(Bim+__odx[0x17])=__t0F-__t2F;\
		*(__B+__odx[0x0F])=__t1E+__t3F;			*(Bim+__odx[0x0F])=__t1F-__t3E;\
		*(__B+__odx[0x1F])=__t1E-__t3F;			*(Bim+__odx[0x1F])=__t1F+__t3E;\
}

#endif	/* #ifndef dft_macro_included */
