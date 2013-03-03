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

/****************************************************************************
 * We now include this header file if it was not included before.
 ****************************************************************************/
#ifndef dft_macro_included
#define dft_macro_included

/* All these radix-N macros are structured so that any 2 of the 3 sets of arguments (A, t, B) may point to the same set of
   addresses, i.e. so the DFT may be done in place (A == B) or re-using a single set of temporaries (A == t or B == t).
   By "in place" we mean that A == B up to an arbitrary permutation of array indices (or element subscripts). */

/****** RADIX = 3: ALLOWS IN-PLACE ******/

/* Totals: 12 FADD, 4 FMUL	*/
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
	__tr0 =__Ar2;					__ti0 =__Ai2;			\
	__tr2 =__Ar1 - __tr0;			__ti2 =__Ai1 - __ti0;	\
	__tr1 =__Ar1 + __tr0;			__ti1 =__Ai1 + __ti0;	\
	__Br0 =__Ar0 + __tr1;			__Bi0 =__Ai0 + __ti1;	\
	__Br1 =__Br0 + __tr1*__c3m1;	__Bi1 =__Bi0 + __ti1*__c3m1;\
	__tr0 =__s*__tr2;				__ti0 =__s*__ti2;			\
	__Br2 =__Br1 + __ti0;			__Bi2 =__Bi1 - __tr0;	\
	__Br1 =__Br1 - __ti0;			__Bi1 =__Bi1 + __tr0;	\
}

/* Totals: 12 FADD, 4 FMUL	*/
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
	__tr0 =__Ar2;					__ti0 =__Ai2;			\
	__tr2 =__Ar1 - __tr0;			__ti2 =__Ai1 - __ti0;	\
	__tr1 =__Ar1 + __tr0;			__ti1 =__Ai1 + __ti0;	\
	__Br0 =__Ar0 + __tr1;			__Bi0 =__Ai0 + __ti1;	\
	__Br1 =__Br0 + __tr1*__c3m1;	__Bi1 =__Bi0 + __ti1*__c3m1;\
	\
	addr = add0 + __add_offset;	prefetch_p_doubles(addr);	\
	\
	__tr0 =__s*__tr2;				__ti0 =__s*__ti2;			\
	__Br2 =__Br1 + __ti0;			__Bi2 =__Bi1 - __tr0;	\
	__Br1 =__Br1 - __ti0;			__Bi1 =__Bi1 + __tr0;	\
}

/****** RADIX = 4: ******/

/* Totals: 16 FADD, 0 FMUL	*/
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
	__rt = __Ar2;	__Ar2 =__Ar0 - __rt;	__Ar0 = __Ar0 + __rt;	\
	__it = __Ai2;	__Ai2 =__Ai0 - __it;	__Ai0 = __Ai0 + __it;	\
	__rt = __Ar3;	__Ar3 =__Ar1 - __rt;	__Ar1 = __Ar1 + __rt;	\
	__it = __Ai3;	__Ai3 =__Ai1 - __it;	__Ai1 = __Ai1 + __it;	\
	\
	__rt  = __Ar1;	__Br1 = __Ar0 - __rt;	__Br0 = __Ar0 + __rt;	\
	__it  = __Ai1;	__Bi1 = __Ai0 - __it;	__Bi0 = __Ai0 + __it;	\
	/* mpy by I is inlined here...	*/\
	__rt  = __Ar3;	__Br3 = __Ar2 + __Ai3;	__Br2 = __Ar2 - __Ai3;	\
			__Bi3 = __Ai2 - __rt;			__Bi2 = __Ai2 + __rt;	\
}

/* Totals: 16 FADD, 0 FMUL	*/
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
	__rt = __Ar2;	__Ar2 =__Ar0 - __rt;	__Ar0 = __Ar0 + __rt;	\
	__it = __Ai2;	__Ai2 =__Ai0 - __it;	__Ai0 = __Ai0 + __it;	\
	__rt = __Ar3;	__Ar3 =__Ar1 - __rt;	__Ar1 = __Ar1 + __rt;	\
	__it = __Ai3;	__Ai3 =__Ai1 - __it;	__Ai1 = __Ai1 + __it;	\
	\
	__addp = __addr + __add_offset0;	prefetch_p_doubles(__addp);	\
	\
	__rt  = __Ar1;	__Br1 = __Ar0 - __rt;	__Br0 = __Ar0 + __rt;	\
	__it  = __Ai1;	__Bi1 = __Ai0 - __it;	__Bi0 = __Ai0 + __it;	\
	/* mpy by I is inlined here...	*/\
	__rt  = __Ar3;	__Br3 = __Ar2 + __Ai3;	__Br2 = __Ar2 - __Ai3;	\
			__Bi3 = __Ai2 - __rt;			__Bi2 = __Ai2 + __rt;	\
	\
	__addp = __addr + __add_offset1;	prefetch_p_doubles(__addp);	\
}

/* Totals: 16 FADD, 0 FMUL	*/
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
	__rt = __Ar1;	__Ar1 =__Ar0 - __rt;	__Ar0 = __Ar0 + __rt;	\
	__it = __Ai1;	__Ai1 =__Ai0 - __it;	__Ai0 = __Ai0 + __it;	\
	__rt = __Ar3;	__Ar3 =__Ar2 - __rt;	__Ar2 = __Ar2 + __rt;	\
	__it = __Ai3;	__Ai3 =__Ai2 - __it;	__Ai2 = __Ai2 + __it;	\
\
	__rt  = __Ar2;	__Br2 = __Ar0 - __rt;	__Br0 = __Ar0 + __rt;	\
	__it  = __Ai2;	__Bi2 = __Ai0 - __it;	__Bi0 = __Ai0 + __it;	\
	/* mpy by I is inlined here...	*/\
	__rt  = __Ar3;	__Br3 = __Ar1 - __Ai3;	__Br1 = __Ar1 + __Ai3;	\
					__Bi3 = __Ai1 + __rt;	__Bi1 = __Ai1 - __rt;	\
}

/****** RADIX = 5: ******/

/* Totals: 34 FADD, 10 FMUL	*/
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
	__Br0 = __Ar0;					__Bi0 = __Ai0;					\
	__rt  = __Ar4;					__it  = __Ai4;					\
	__Br4 = __Ar1 - __rt;			__Bi4 = __Ai1 - __it;			\
	__Br1 = __Ar1 + __rt;			__Bi1 = __Ai1 + __it;			\
	__rt  = __Ar3;					__it  = __Ai3;					\
	__Br3 = __Ar2 - __rt;			__Bi3 = __Ai2 - __it;			\
	__Br2 = __Ar2 + __rt;			__Bi2 = __Ai2 + __it;			\
	__rt  = __Br1 + __Br2;			__it  = __Bi1 + __Bi2;			\
	__Br0 = __Br0 + __rt;			__Bi0 = __Bi0 + __it;		/* y0 */\
	__rt  = __Br0 + __rt  *__cc1;	__it  = __Bi0 + __it  *__cc1;	\
	__Br2 =(__Br1 - __Br2)*__cc2;	__Bi2 =(__Bi1 - __Bi2)*__cc2;	\
	__Br1 = __rt  + __Br2;			__Bi1 = __it  + __Bi2;			\
	__Br2 = __rt  - __Br2;			__Bi2 = __it  - __Bi2;			\
	__rt  =(__Br4 - __Br3)* __s2;	__it  =(__Bi4 - __Bi3)* __s2;	\
	__Br3 = __rt  + __Br3 *__ss1;	__Bi3 = __it  + __Bi3 *__ss1;	\
	__Br4 = __rt  - __Br4 *__ss2;	__Bi4 = __it  - __Bi4 *__ss2;	\
	__rt  = __Br4;					__it  = __Bi4;					\
	__Br4 = __Br1 + __Bi3;			__Bi4 = __Bi1 - __Br3;		/* y4 - note swap with y3 below! */\
	__Br1 = __Br1 - __Bi3;			__Bi1 = __Bi1 + __Br3;		/* y1 */\
	__Br3 = __Br2 + __it;			__Bi3 = __Bi2 - __rt;		/* y3 */\
	__Br2 = __Br2 - __it;			__Bi2 = __Bi2 + __rt;		/* y2 */\
}

/* Totals: 34 FADD, 10 FMUL	*/
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
	__Br0 = __Ar0;					__Bi0 = __Ai0;					\
	__rt  = __Ar4;					__it  = __Ai4;					\
	__Br4 = __Ar1 - __rt;			__Bi4 = __Ai1 - __it;			\
	__Br1 = __Ar1 + __rt;			__Bi1 = __Ai1 + __it;			\
	__rt  = __Ar3;					__it  = __Ai3;					\
	__Br3 = __Ar2 - __rt;			__Bi3 = __Ai2 - __it;			\
	__Br2 = __Ar2 + __rt;			__Bi2 = __Ai2 + __it;			\
	__rt  = __Br1 + __Br2;			__it  = __Bi1 + __Bi2;			\
	__Br0 = __Br0 + __rt;			__Bi0 = __Bi0 + __it;		/* y0 */\
	__rt  = __Br0 + __rt  *__cc1;	__it  = __Bi0 + __it  *__cc1;	\
	__Br2 =(__Br1 - __Br2)*__cc2;	__Bi2 =(__Bi1 - __Bi2)*__cc2;	\
	\
	__addp = __addr + __add_offset0;	prefetch_p_doubles(__addp);	\
	\
	__Br1 = __rt  + __Br2;			__Bi1 = __it  + __Bi2;			\
	__Br2 = __rt  - __Br2;			__Bi2 = __it  - __Bi2;			\
	__rt  =(__Br4 - __Br3)* __s2;	__it  =(__Bi4 - __Bi3)* __s2;	\
	__Br3 = __rt  + __Br3 *__ss1;	__Bi3 = __it  + __Bi3 *__ss1;	\
	__Br4 = __rt  - __Br4 *__ss2;	__Bi4 = __it  - __Bi4 *__ss2;	\
	__rt  = __Br4;					__it  = __Bi4;					\
	__Br4 = __Br1 + __Bi3;			__Bi4 = __Bi1 - __Br3;		/* y4 - note swap with y3 below! */\
	__Br1 = __Br1 - __Bi3;			__Bi1 = __Bi1 + __Br3;		/* y1 */\
	__Br3 = __Br2 + __it;			__Bi3 = __Bi2 - __rt;		/* y3 */\
	__Br2 = __Br2 - __it;			__Bi2 = __Bi2 + __rt;		/* y2 */\
	\
	__addp = __addr + __add_offset1;	prefetch_p_doubles(__addp);	\
}

/****** RADIX = 7: ******/

/* Totals: 60 FADD, 36 FMUL	*/
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
	__tr0 = __Ar0;												__ti0 = __Ai0;				\
	__tr6 = __Ar1 - __Ar6;										__ti6 = __Ai1 - __Ai6;	/* x1 - x6 */	\
	__tr1 = __Ar1 + __Ar6;										__ti1 = __Ai1 + __Ai6;	/* x1 + x6 */	\
	\
	__tr5 = __Ar2 - __Ar5;										__ti5 = __Ai2 - __Ai5;	/* x2 - x5 */	\
	__tr2 = __Ar2 + __Ar5;										__ti2 = __Ai2 + __Ai5;	/* x2 + x5 */	\
	\
	__tr4 = __Ar3 - __Ar4;										__ti4 = __Ai3 - __Ai4;	/* x3 - x4 */	\
	__tr3 = __Ar3 + __Ar4;										__ti3 = __Ai3 + __Ai4;	/* x3 + x4 */	\
	\
	__Br0 = __tr0 + __tr1 + __tr2 + __tr3;						__Bi0 = __ti0 + __ti1 + __ti2 + __ti3;			\
	__rt  = __tr0 + __cc1*__tr1 + __cc2*__tr2 + __cc3*__tr3;	__it  = __ti0 + __cc1*__ti1 + __cc2*__ti2 + __cc3*__ti3;	\
	__re  = __tr0 + __cc2*__tr1 + __cc3*__tr2 + __cc1*__tr3;	__im  = __ti0 + __cc2*__ti1 + __cc3*__ti2 + __cc1*__ti3;	\
	__tr0 = __tr0 + __cc3*__tr1 + __cc1*__tr2 + __cc2*__tr3;	__ti0 = __ti0 + __cc3*__ti1 + __cc1*__ti2 + __cc2*__ti3;	\
	\
	__tr1 = __ss1*__tr6 + __ss2*__tr5 + __ss3*__tr4;			__ti1 = __ss1*__ti6 + __ss2*__ti5 + __ss3*__ti4;	\
	__tr2 = __ss2*__tr6 - __ss3*__tr5 - __ss1*__tr4;			__ti2 = __ss2*__ti6 - __ss3*__ti5 - __ss1*__ti4;	\
	__tr3 = __ss3*__tr6 - __ss1*__tr5 + __ss2*__tr4;			__ti3 = __ss3*__ti6 - __ss1*__ti5 + __ss2*__ti4;	\
	/* Output permutation causes signs to get flipped here: */\
	__Br1 = __rt  - __ti1;										__Bi1 = __it  + __tr1;	\
	__Br2 = __re  - __ti2;										__Bi2 = __im  + __tr2;	\
	__Br3 = __tr0 - __ti3;										__Bi3 = __ti0 + __tr3;	\
	__Br4 = __tr0 + __ti3;										__Bi4 = __ti0 - __tr3;	\
	__Br5 = __re  + __ti2;										__Bi5 = __im  - __tr2;	\
	__Br6 = __rt  + __ti1;										__Bi6 = __it  - __tr1;	\
}

/* Totals: 72 FADD, 16 FMUL	*/
#define RADIX_07_DFT_NUSSBAUMER(\
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
	__t0r = __A0r;						__t0i = __A0i;				\
	__t6r = __A1r - __A6r;				__t6i = __A1i - __A6i;	/* x1 - x6 */	\
	__t1r = __A1r + __A6r;				__t1i = __A1i + __A6i;	/* x1 + x6 */	\
	\
	__t5r = __A2r - __A5r;				__t5i = __A2i - __A5i;	/* x2 - x5 */	\
	__t2r = __A2r + __A5r;				__t2i = __A2i + __A5i;	/* x2 + x5 */	\
	\
	__t4r = __A3r - __A4r;				__t4i = __A3i - __A4i;	/* x3 - x4 */	\
	__t3r = __A3r + __A4r;				__t3i = __A3i + __A4i;	/* x3 + x4 */	\
	\
	__rt = __t1r+__t2r+__t3r;			__it = __t1i+__t2i+__t3i;		\
	__B0r= __rt+__t0r;					__B0i= __it+__t0i;				\
	__t0r= __rt*__cx0+__t0r;			__t0i= __it*__cx0+__t0i;		\
	__t1r= __t1r-__t2r;					__t1i= __t1i-__t2i;				\
	__t2r= __t3r-__t2r;					__t2i= __t3i-__t2i;				\
	__t3r=(__t1r+__t2r)*__cx3;			__t3i=(__t1i+__t2i)*__cx3;		\
	__t1r= __t1r*__cx1;					__t1i= __t1i*__cx1;				\
	__t2r= __t2r*__cx2;					__t2i= __t2i*__cx2;				\
	__rt = __t1r-__t3r;					__it = __t1i-__t3i;				\
	__t2r= __t2r-__t3r;					__t2i= __t2i-__t3i;				\
																		\
	__t1r= __t0r-__rt-__t2r;			__t1i= __t0i-__it-__t2i;		\
	__t2r= __t0r+__t2r;					__t2i= __t0i+__t2i;				\
	__t0r= __t0r+__rt;					__t0i= __t0i+__it;				\
																		\
	__t3r=(__t6r-__t4r+__t5r)*__sx0;	__t3i=(__t6i-__t4i+__t5i)*__sx0;\
	__t6r= __t6r-__t5r;					__t6i= __t6i-__t5i;				\
	__t5r= __t4r+__t5r;					__t5i= __t4i+__t5i;				\
	__t4r=(__t5r-__t6r)*__sx3;			__t4i=(__t5i-__t6i)*__sx3;		\
	__t6r= __t6r*__sx1;					__t6i= __t6i*__sx1;				\
	__t5r= __t5r*__sx2;					__t5i= __t5i*__sx2;				\
	__t6r= __t4r+__t6r;					__t6i= __t4i+__t6i;				\
	__t5r= __t4r-__t5r;					__t5i= __t4i-__t5i;				\
																		\
	__t4r= __t3r-__t6r-__t5r;			__t4i= __t3i-__t6i-__t5i;		\
	__t5r= __t3r+__t5r;					__t5i= __t3i+__t5i;				\
	__t3r= __t3r+__t6r;					__t3i= __t3i+__t6i;				\
																		\
	__B1r =__t0r-__t3i;					__B1i =__t0i+__t3r;				\
	__B2r =__t1r-__t4i;					__B2i =__t1i+__t4r;				\
	__B3r =__t2r+__t5i;					__B3i =__t2i-__t5r;				\
	__B4r =__t2r-__t5i;					__B4i =__t2i+__t5r;				\
	__B5r =__t1r+__t4i;					__B5i =__t1i-__t4r;				\
	__B6r =__t0r+__t3i;					__B6i =__t0i-__t3r;				\
}

/* Totals: 60 FADD, 36 FMUL	*/
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
	__tr0 = __Ar0;												__ti0 = __Ai0;				\
	__tr6 = __Ar1 - __Ar6;										__ti6 = __Ai1 - __Ai6;	/* x1 - x6 */	\
	__tr1 = __Ar1 + __Ar6;										__ti1 = __Ai1 + __Ai6;	/* x1 + x6 */	\
	\
	__tr5 = __Ar2 - __Ar5;										__ti5 = __Ai2 - __Ai5;	/* x2 - x5 */	\
	__tr2 = __Ar2 + __Ar5;										__ti2 = __Ai2 + __Ai5;	/* x2 + x5 */	\
	\
	__tr4 = __Ar3 - __Ar4;										__ti4 = __Ai3 - __Ai4;	/* x3 - x4 */	\
	__tr3 = __Ar3 + __Ar4;										__ti3 = __Ai3 + __Ai4;	/* x3 + x4 */	\
	\
	__addp = __addr + __add_offset0;	prefetch_p_doubles(__addp);	\
	\
	__Br0 = __tr0 + __tr1 + __tr2 + __tr3;						__Bi0 = __ti0 + __ti1 + __ti2 + __ti3;			\
	__rt  = __tr0 + __cc1*__tr1 + __cc2*__tr2 + __cc3*__tr3;	__it  = __ti0 + __cc1*__ti1 + __cc2*__ti2 + __cc3*__ti3;	\
	__re  = __tr0 + __cc2*__tr1 + __cc3*__tr2 + __cc1*__tr3;	__im  = __ti0 + __cc2*__ti1 + __cc3*__ti2 + __cc1*__ti3;	\
	__tr0 = __tr0 + __cc3*__tr1 + __cc1*__tr2 + __cc2*__tr3;	__ti0 = __ti0 + __cc3*__ti1 + __cc1*__ti2 + __cc2*__ti3;	\
	\
	__addp = __addr + __add_offset1;	prefetch_p_doubles(__addp);	\
	\
	__tr1 = __ss1*__tr6 + __ss2*__tr5 + __ss3*__tr4;			__ti1 = __ss1*__ti6 + __ss2*__ti5 + __ss3*__ti4;	\
	__tr2 = __ss2*__tr6 - __ss3*__tr5 - __ss1*__tr4;			__ti2 = __ss2*__ti6 - __ss3*__ti5 - __ss1*__ti4;	\
	__tr3 = __ss3*__tr6 - __ss1*__tr5 + __ss2*__tr4;			__ti3 = __ss3*__ti6 - __ss1*__ti5 + __ss2*__ti4;	\
	\
	__addp = __addr + __add_offset2;	prefetch_p_doubles(__addp);	\
	\
	/* Output permutation causes signs to get flipped here: */\
	__Br1 = __rt  - __ti1;										__Bi1 = __it  + __tr1;	\
	__Br2 = __re  - __ti2;										__Bi2 = __im  + __tr2;	\
	__Br3 = __tr0 - __ti3;										__Bi3 = __ti0 + __tr3;	\
	__Br4 = __tr0 + __ti3;										__Bi4 = __ti0 - __tr3;	\
	__Br5 = __re  + __ti2;										__Bi5 = __im  - __tr2;	\
	__Br6 = __rt  + __ti1;										__Bi6 = __it  - __tr1;	\
}

/****** RADIX = 8: ******/

/* Totals: 52 FADD, 4 FMUL	*/
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
	__tr0 = __Ar0;			__ti0 = __Ai0;			\
	__rt  = __Ar4;			__it  = __Ai4;			\
	__tr1 = __tr0 - __rt;		__ti1 = __ti0 - __it;		\
	__tr0 = __tr0 + __rt;		__ti0 = __ti0 + __it;		\
\
	__tr2 = __Ar2;			__ti2 = __Ai2;			\
	__rt  = __Ar6;			__it  = __Ai6;			\
	__tr3 = __tr2 - __rt;		__ti3 = __ti2 - __it;		\
	__tr2 = __tr2 + __rt;		__ti2 = __ti2 + __it;		\
\
	__tr4 = __Ar1;			__ti4 = __Ai1;			\
	__rt  = __Ar5;			__it  = __Ai5;			\
	__tr5 = __tr4 - __rt;		__ti5 = __ti4 - __it;		\
	__tr4 = __tr4 + __rt;		__ti4 = __ti4 + __it;		\
\
	__tr6 = __Ar3;			__ti6 = __Ai3;			\
	__rt  = __Ar7;			__it  = __Ai7;			\
	__tr7 = __tr6 - __rt;		__ti7 = __ti6 - __it;		\
	__tr6 = __tr6 + __rt;		__ti6 = __ti6 + __it;		\
/*       combine to get the 2 length-4 transforms...	*/		\
	__rt  = __tr2;			__it  = __ti2;			\
	__tr2 = __tr0 - __rt;		__ti2 = __ti0 - __it;		\
	__tr0 = __tr0 + __rt;		__ti0 = __ti0 + __it;		\
\
	__rt  = __tr3;			__it  = __ti3;			\
	__tr3 = __tr1 + __it;		__ti3 = __ti1 - __rt;		\
	__tr1 = __tr1 - __it;		__ti1 = __ti1 + __rt;		\
\
	__rt  = __tr6;			__it  = __ti6;			\
	__tr6 = __tr4 - __rt;		__ti6 = __ti4 - __it;		\
	__tr4 = __tr4 + __rt;		__ti4 = __ti4 + __it;		\
\
	__rt  = __tr7;			__it  = __ti7;			\
	__tr7 = __tr5 + __it;		__ti7 = __ti5 - __rt;		\
	__tr5 = __tr5 - __it;		__ti5 = __ti5 + __rt;		\
/*       now combine the two half-transforms	*/			\
	__Br0 = __tr0 + __tr4;		__Bi0 = __ti0 + __ti4;		\
	__Br1 = __tr0 - __tr4;		__Bi1 = __ti0 - __ti4;		\
\
	__Br2 = __tr2 - __ti6;		__Bi2 = __ti2 + __tr6;		\
	__Br3 = __tr2 + __ti6;		__Bi3 = __ti2 - __tr6;		\
\
	__rt = (__tr5 - __ti5)*ISRT2;	__it = (__tr5 + __ti5)*ISRT2;	\
	__Br4 = __tr1 + __rt;		__Bi4 = __ti1 + __it;		\
	__Br5 = __tr1 - __rt;		__Bi5 = __ti1 - __it;		\
\
	__rt = (__tr7 + __ti7)*ISRT2;	__it = (__ti7 - __tr7)*ISRT2;	\
	__Br6 = __tr3 - __rt;		__Bi6 = __ti3 - __it;		\
	__Br7 = __tr3 + __rt;		__Bi7 = __ti3 + __it;		\
}

/* Totals: 52 FADD, 4 FMUL	*/
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
	__tr0 = __Ar0;			__ti0 = __Ai0;			\
	__rt  = __Ar4;			__it  = __Ai4;			\
	__tr1 = __tr0 - __rt;		__ti1 = __ti0 - __it;		\
	__tr0 = __tr0 + __rt;		__ti0 = __ti0 + __it;		\
\
	__tr2 = __Ar2;			__ti2 = __Ai2;			\
	__rt  = __Ar6;			__it  = __Ai6;			\
	__tr3 = __tr2 - __rt;		__ti3 = __ti2 - __it;		\
	__tr2 = __tr2 + __rt;		__ti2 = __ti2 + __it;		\
\
addr = add0 + __add_offset0;	prefetch_p_doubles(addr);	\
\
	__tr4 = __Ar1;			__ti4 = __Ai1;			\
	__rt  = __Ar5;			__it  = __Ai5;			\
	__tr5 = __tr4 - __rt;		__ti5 = __ti4 - __it;		\
	__tr4 = __tr4 + __rt;		__ti4 = __ti4 + __it;		\
\
	__tr6 = __Ar3;			__ti6 = __Ai3;			\
	__rt  = __Ar7;			__it  = __Ai7;			\
	__tr7 = __tr6 - __rt;		__ti7 = __ti6 - __it;		\
	__tr6 = __tr6 + __rt;		__ti6 = __ti6 + __it;		\
\
addr = add0 + __add_offset1;	prefetch_p_doubles(addr);	\
\
/*       combine to get the 2 length-4 transforms...	*/		\
	__rt  = __tr2;			__it  = __ti2;			\
	__tr2 = __tr0 - __rt;		__ti2 = __ti0 - __it;		\
	__tr0 = __tr0 + __rt;		__ti0 = __ti0 + __it;		\
\
	__rt  = __tr3;			__it  = __ti3;			\
	__tr3 = __tr1 + __it;		__ti3 = __ti1 - __rt;		\
	__tr1 = __tr1 - __it;		__ti1 = __ti1 + __rt;		\
\
addr = add0 + __add_offset2;	prefetch_p_doubles(addr);	\
\
	__rt  = __tr6;			__it  = __ti6;			\
	__tr6 = __tr4 - __rt;		__ti6 = __ti4 - __it;		\
	__tr4 = __tr4 + __rt;		__ti4 = __ti4 + __it;		\
\
	__rt  = __tr7;			__it  = __ti7;			\
	__tr7 = __tr5 + __it;		__ti7 = __ti5 - __rt;		\
	__tr5 = __tr5 - __it;		__ti5 = __ti5 + __rt;		\
\
addr = add0 + __add_offset3;	prefetch_p_doubles(addr);	\
\
/*       now combine the two half-transforms	*/			\
	__Br0 = __tr0 + __tr4;		__Bi0 = __ti0 + __ti4;		\
	__Br1 = __tr0 - __tr4;		__Bi1 = __ti0 - __ti4;		\
\
	__Br2 = __tr2 - __ti6;		__Bi2 = __ti2 + __tr6;		\
	__Br3 = __tr2 + __ti6;		__Bi3 = __ti2 - __tr6;		\
\
addr = add0 + __add_offset4;	prefetch_p_doubles(addr);	\
\
	__rt = (__tr5 - __ti5)*ISRT2;	__it = (__tr5 + __ti5)*ISRT2;	\
	__Br4 = __tr1 + __rt;		__Bi4 = __ti1 + __it;		\
	__Br5 = __tr1 - __rt;		__Bi5 = __ti1 - __it;		\
\
	__rt = (__tr7 + __ti7)*ISRT2;	__it = (__ti7 - __tr7)*ISRT2;	\
	__Br6 = __tr3 - __rt;		__Bi6 = __ti3 - __it;		\
	__Br7 = __tr3 + __rt;		__Bi7 = __ti3 + __it;		\
}

/* Totals: 52 FADD, 4 FMUL	*/
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
	__tr0 = __Ar0;			__ti0 = __Ai0;			\
	__rt  = __Ar1;			__it  = __Ai1;			\
	__tr1 = __tr0 - __rt;		__ti1 = __ti0 - __it;		\
	__tr0 = __tr0 + __rt;		__ti0 = __ti0 + __it;		\
\
	__tr2 = __Ar2;			__ti2 = __Ai2;			\
	__rt  = __Ar3;			__it  = __Ai3;			\
	__tr3 = __tr2 - __rt;		__ti3 = __ti2 - __it;		\
	__tr2 = __tr2 + __rt;		__ti2 = __ti2 + __it;		\
\
	__tr4 = __Ar4;			__ti4 = __Ai4;			\
	__rt  = __Ar5;			__it  = __Ai5;			\
	__tr5 = __tr4 - __rt;		__ti5 = __ti4 - __it;		\
	__tr4 = __tr4 + __rt;		__ti4 = __ti4 + __it;		\
\
	__tr6 = __Ar6;			__ti6 = __Ai6;			\
	__rt  = __Ar7;			__it  = __Ai7;			\
	__tr7 = __tr6 - __rt;		__ti7 = __ti6 - __it;		\
	__tr6 = __tr6 + __rt;		__ti6 = __ti6 + __it;		\
/*       combine to get the 2 length-4 transforms...	*/		\
	__rt  = __tr2;			__it  = __ti2;			\
	__tr2 = __tr0 - __rt;		__ti2 = __ti0 - __it;		\
	__tr0 = __tr0 + __rt;		__ti0 = __ti0 + __it;		\
\
	__rt  = __tr3;			__it  = __ti3;			\
	__tr3 = __tr1 + __it;		__ti3 = __ti1 - __rt;		\
	__tr1 = __tr1 - __it;		__ti1 = __ti1 + __rt;		\
\
	__rt  = __tr6;			__it  = __ti6;			\
	__tr6 = __tr4 - __rt;		__ti6 = __ti4 - __it;		\
	__tr4 = __tr4 + __rt;		__ti4 = __ti4 + __it;		\
\
	__rt  = __tr7;			__it  = __ti7;			\
	__tr7 = __tr5 + __it;		__ti7 = __ti5 - __rt;		\
	__tr5 = __tr5 - __it;		__ti5 = __ti5 + __rt;		\
/*       now combine the two half-transforms	*/			\
	__Br0 = __tr0 + __tr4;		__Bi0 = __ti0 + __ti4;		\
	__Br4 = __tr0 - __tr4;		__Bi4 = __ti0 - __ti4;		\
\
	__Br2 = __tr2 - __ti6;		__Bi2 = __ti2 + __tr6;		\
	__Br6 = __tr2 + __ti6;		__Bi6 = __ti2 - __tr6;		\
\
	__rt = (__tr5 - __ti5)*ISRT2;	__it = (__tr5 + __ti5)*ISRT2;	\
	__Br1 = __tr1 + __rt;		__Bi1 = __ti1 + __it;		\
	__Br5 = __tr1 - __rt;		__Bi5 = __ti1 - __it;		\
\
	__rt = (__tr7 + __ti7)*ISRT2;	__it = (__ti7 - __tr7)*ISRT2;	\
	__Br3 = __tr3 - __rt;		__Bi3 = __ti3 - __it;		\
	__Br7 = __tr3 + __rt;		__Bi7 = __ti3 + __it;		\
}


/****** RADIX = 9: ******/

/*...Radix-9 DIF: A/Bs are in/outputs and t's are temporaries (all doubles):	*/
/* Totals: 80 FADD, 40 FMUL	*/
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
/*...and now do three more radix-3 transforms, including the twiddle factors:	*/\
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
/* Totals: 68 FADD, 40 FMUL	*/
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
/*...and now do three more radix-3 transforms, including the twiddle factors:	*/\
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

/*...Simple-algo Radix-11 DFT.  Totals: 140 FADD, 100 FMUL	*/
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

/*...Radix-11 DFT using length-5 cyclic convolution scheme for the 5 x 5 matrix submultiplies.  Totals: 160 FADD,  44 FMUL	*/
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
		xr7 = __A6r + __A7r;			xi7 = __A6i + __A7i;		\
		xr5 = __A5r + __A8r;			xi5 = __A5i + __A8i;		\
		xr4 = __A4r + __A9r;			xi4 = __A4i + __A9i;		\
		xr6 = __A3r + __Aar;			xi6 = __A3i + __Aai;		\
		xr3 = __A2r + __Abr;			xi3 = __A2i + __Abi;		\
		xr1 = __A1r + __Acr;			xi1 = __A1i + __Aci;		\
\
		yr1 = __A1r - __Acr;			yi1 = __A1i - __Aci;		\
		yr2 = __A2r - __Abr;			yi2 = __A2i - __Abi;		\
		yr5 = __A3r - __Aar;			yi5 = __A3i - __Aai;		\
		yr3 = __A4r - __A9r;			yi3 = __A4i - __A9i;		\
		yr4 = __A5r - __A8r;			yi4 = __A5i - __A8i;		\
		yr6 = __A6r - __A7r;			yi6 = __A6i - __A7i;		\
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
/* Totals: 164 FADD, 64 FMUL. */\
}

/****** RADIX = 15: ******/

/* Totals: 162 FADD, 50 FMUL	*/
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
	__tr00 = __Ar00;				__ti00 = __Ai00;		\
	__tr01 = __Ar12;				__ti01 = __Ai12;		\
	__rt   = __Ar03;				__it   = __Ai03;		\
	__tr03 = __tr01-__rt;			__ti03 = __ti01-__it;		\
	__tr01 = __tr01+__rt;			__ti01 = __ti01+__it;		\
	__tr02 = __Ar09;				__ti02 = __Ai09;		\
	__rt   = __Ar06;				__it   = __Ai06;		\
	__tr04 = __tr02-__rt;			__ti04 = __ti02-__it;		\
	__tr02 = __tr02+__rt;			__ti02 = __ti02+__it;		\
									\
	__rt   = __tr01+__tr02;			__it   = __ti01+__ti02;		\
	__tr00 = __tr00+__rt;			__ti00 = __ti00+__it;		\
	__rt   = __tr00+cn1*__rt;		__it   = __ti00+cn1*__it;	\
	__tr02 = cn2*(__tr01-__tr02);	__ti02 = cn2*(__ti01-__ti02);	\
	__tr01 = __rt+__tr02;			__ti01 = __it+__ti02;		\
	__tr02 = __rt-__tr02;			__ti02 = __it-__ti02;		\
	__rt   = ss3*(__tr03-__tr04);	__it   = ss3*(__ti03-__ti04);	\
	__tr04 = __rt+sn1*__tr04;		__ti04 = __it+sn1*__ti04;	\
	__tr03 = __rt-sn2*__tr03;		__ti03 = __it-sn2*__ti03;	\
	__rt   = __tr04;				__it   = __ti04;		\
	__tr04 = __tr01+__it;			__ti04 = __ti01-__rt;		/*<==prefer these to be stored in __tr,i04 */\
	__tr01 = __tr01-__it;			__ti01 = __ti01+__rt;		\
	__rt   = __tr03;				__it   = __ti03;		\
	__tr03 = __tr02+__it;			__ti03 = __ti02-__rt;		/*<==prefer these to be stored in __tr,i03 */\
	__tr02 = __tr02-__it;			__ti02 = __ti02+__rt;		\
									\
/*...Block 2: Swap inputs as 01 <-> 10, 04 <-> 07: */\
	__tr05 = __Ar10;				__ti05 = __Ai10;		\
	__tr06 = __Ar07;				__ti06 = __Ai07;		\
	__rt   = __Ar13;				__it   = __Ai13;		\
	__tr08 = __tr06 -__rt;			__ti08 = __ti06 -__it;		\
	__tr06 = __tr06 +__rt;			__ti06 = __ti06 +__it;		\
	__tr07 = __Ar04;				__ti07 = __Ai04;		\
	__rt   = __Ar01;				__it   = __Ai01;		\
	__tr09 = __tr07 -__rt;			__ti09 = __ti07 -__it;		\
	__tr07 = __tr07 +__rt;			__ti07 = __ti07 +__it;		\
									\
	__rt   = __tr06+__tr07;			__it   = __ti06+__ti07;		\
	__tr05 = __tr05+__rt;			__ti05 = __ti05+__it;		\
	__rt   = __tr05+cn1*__rt;		__it   = __ti05+cn1*__it;	\
	__tr07 = cn2*(__tr06-__tr07);	__ti07 = cn2*(__ti06-__ti07);	\
	__tr06 = __rt+__tr07;			__ti06 = __it+__ti07;		\
	__tr07 = __rt-__tr07;			__ti07 = __it-__ti07;		\
	__rt   = ss3*(__tr08-__tr09);	__it   = ss3*(__ti08-__ti09);	\
	__tr09 = __rt+sn1*__tr09;		__ti09 = __it+sn1*__ti09;	\
	__tr08 = __rt-sn2*__tr08;		__ti08 = __it-sn2*__ti08;	\
	__rt   = __tr09;				__it   = __ti09;		\
	__tr09 = __tr06+__it;			__ti09 = __ti06-__rt;		\
	__tr06 = __tr06-__it;			__ti06 = __ti06+__rt;		\
	__rt   = __tr08;				__it   = __ti08;		\
	__tr08 = __tr07+__it;			__ti08 = __ti07-__rt;		\
	__tr07 = __tr07-__it;			__ti07 = __ti07+__rt;		\
									\
/*...Block 3: Swap inputs as 02 <-> 05, 08 <-> 14: */\
	__tr10 = __Ar05;				__ti10 = __Ai05;		\
	__tr11 = __Ar02;				__ti11 = __Ai02;		\
	__rt   = __Ar08;				__it   = __Ai08;		\
	__tr13 = __tr11 -__rt;			__ti13 = __ti11 -__it;		\
	__tr11 = __tr11 +__rt;			__ti11 = __ti11 +__it;		\
	__tr12 = __Ar14;				__ti12 = __Ai14;		\
	__rt   = __Ar11;				__it   = __Ai11;		\
	__tr14 = __tr12 -__rt;			__ti14 = __ti12 -__it;		\
	__tr12 = __tr12 +__rt;			__ti12 = __ti12 +__it;		\
									\
	__rt   = __tr11+__tr12;			__it   = __ti11+__ti12;		\
	__tr10 = __tr10+__rt;			__ti10 = __ti10+__it;		\
	__rt   = __tr10+cn1*__rt;		__it   = __ti10+cn1*__it;	\
	__tr12 = cn2*(__tr11-__tr12);	__ti12 = cn2*(__ti11-__ti12);	\
	__tr11 = __rt+__tr12;			__ti11 = __it+__ti12;		\
	__tr12 = __rt-__tr12;			__ti12 = __it-__ti12;		\
	__rt   = ss3*(__tr13-__tr14);	__it   = ss3*(__ti13-__ti14);	\
	__tr14 = __rt+sn1*__tr14;		__ti14 = __it+sn1*__ti14;	\
	__tr13 = __rt-sn2*__tr13;		__ti13 = __it-sn2*__ti13;	\
	__rt   = __tr14;				__it   = __ti14;		\
	__tr14 = __tr11+__it;			__ti14 = __ti11-__rt;		\
	__tr11 = __tr11-__it;			__ti11 = __ti11+__rt;		\
	__rt   = __tr13;				__it   = __ti13;		\
	__tr13 = __tr12+__it;			__ti13 = __ti12-__rt;		\
	__tr12 = __tr12-__it;			__ti12 = __ti12+__rt;		\
									\
/*...and now do five radix-3 transforms.	\
	The required output permutation is [0,1,2,13,14,12,9,10,11,8,6,7,4,5,3]. */\
/*...Block 1:	*/					\
	__rt   = __tr10;				__it   = __ti10;		\
	__tr10 = __tr05-__rt;			__ti10 = __ti05-__it;		\
	__tr05 = __tr05+__rt;			__ti05 = __ti05+__it;		\
	__tr00 = __tr00+__tr05;			__ti00 = __ti00+__ti05;		\
	__Br00 = __tr00;				__Bi00 = __ti00;		\
	__tr05 = __tr00+c3m1*__tr05;	__ti05 = __ti00+c3m1*__ti05;	\
	__rt   = s*__tr10;				__it   = s*__ti10;		\
	__Br01 = __tr05-__it;			__Bi01 = __ti05+__rt;		\
	__Br02 = __tr05+__it;			__Bi02 = __ti05-__rt;		\
									\
/*...Block 2:	*/					\
	__rt   = __tr11;				__it   = __ti11;		\
	__tr11 = __tr06-__rt;			__ti11 = __ti06-__it;		\
	__tr06 = __tr06+__rt;			__ti06 = __ti06+__it;		\
	__tr01 = __tr01+__tr06;			__ti01 = __ti01+__ti06;		\
	__Br13 = __tr01;				__Bi13 = __ti01;		\
	__tr06 = __tr01+c3m1*__tr06;	__ti06 = __ti01+c3m1*__ti06;	\
	__rt   = s*__tr11;				__it   = s*__ti11;		\
	__Br14 = __tr06-__it;			__Bi14 = __ti06+__rt;		\
	__Br12 = __tr06+__it;			__Bi12 = __ti06-__rt;		\
									\
/*...Block 3:	*/					\
	__rt   = __tr12;				__it   = __ti12;		\
	__tr12 = __tr07-__rt;			__ti12 = __ti07-__it;		\
	__tr07 = __tr07+__rt;			__ti07 = __ti07+__it;		\
	__tr02 = __tr02+__tr07;			__ti02 = __ti02+__ti07;		\
	__Br09 = __tr02;				__Bi09 = __ti02;		\
	__tr07 = __tr02+c3m1*__tr07;	__ti07 = __ti02+c3m1*__ti07;	\
	__rt   = s*__tr12;				__it   = s*__ti12;		\
	__Br10 = __tr07-__it;			__Bi10 = __ti07+__rt;		\
	__Br11 = __tr07+__it;			__Bi11 = __ti07-__rt;		\
									\
/*...Block 4:	*/					\
	__rt   = __tr13;				__it   = __ti13;		\
	__tr13 = __tr08-__rt;			__ti13 = __ti08-__it;		\
	__tr08 = __tr08+__rt;			__ti08 = __ti08+__it;		\
	__tr03 = __tr03+__tr08;			__ti03 = __ti03+__ti08;		\
	__Br08 = __tr03;				__Bi08 = __ti03;		\
	__tr08 = __tr03+c3m1*__tr08;	__ti08 = __ti03+c3m1*__ti08;	\
	__rt   = s*__tr13;				__it   = s*__ti13;		\
	__Br06 = __tr08-__it;			__Bi06 = __ti08+__rt;		\
	__Br07 = __tr08+__it;			__Bi07 = __ti08-__rt;		\
									\
/*...Block 5:	*/					\
	__rt   = __tr14;				__it   = __ti14;		\
	__tr14 = __tr09-__rt;			__ti14 = __ti09-__it;		\
	__tr09 = __tr09+__rt;			__ti09 = __ti09+__it;		\
	__tr04 = __tr04+__tr09;			__ti04 = __ti04+__ti09;		\
	__Br04 = __tr04;				__Bi04 = __ti04;		\
	__tr09 = __tr04+c3m1*__tr09;	__ti09 = __ti04+c3m1*__ti09;	\
	__rt   = s*__tr14;				__it   = s*__ti14;		\
	__Br05 = __tr09-__it;			__Bi05 = __ti09+__rt;		\
	__Br03 = __tr09+__it;			__Bi03 = __ti09-__rt;		\
}

/* Totals: 162 FADD, 50 FMUL	*/
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
	    __tr00 = __Ar00;				__ti00 = __Ai00;		\
	    __tr01 = __Ar02;				__ti01 = __Ai02;		\
	    __rt   = __Ar01;				__it   = __Ai01;		\
	    __tr02 = __tr01-__rt;			__ti02 = __ti01-__it;		\
	    __tr01 = __tr01+__rt;			__ti01 = __ti01+__it;		\
	    __tr00 = __tr00+__tr01;			__ti00 = __ti00+__ti01;		\
	    __tr01 = __tr00+c3m1*__tr01;	__ti01 = __ti00+c3m1*__ti01;	\
	    __rt   = s*__tr02;				__it   = s*__ti02;		\
	    __tr02 = __tr01-__it;			__ti02 = __ti01+__rt;		\
	    __tr01 = __tr01+__it;			__ti01 = __ti01-__rt;		\
										\
/*...Block 2: Swap inputs 3,4,5 -> 8,7,6: */\
	    __tr03 = __Ar08;				__ti03 = __Ai08;		\
	    __tr04 = __Ar07;				__ti04 = __Ai07;		\
	    __rt   = __Ar06;				__it   = __Ai06;		\
	    __tr05 = __tr04-__rt;			__ti05 = __ti04-__it;		\
	    __tr04 = __tr04+__rt;			__ti04 = __ti04+__it;		\
	    __tr03 = __tr03+__tr04;			__ti03 = __ti03+__ti04;		\
	    __tr04 = __tr03+c3m1*__tr04;	__ti04 = __ti03+c3m1*__ti04;	\
	    __rt   = s*__tr05;				__it   = s*__ti05;		\
	    __tr05 = __tr04-__it;			__ti05 = __ti04+__rt;		\
	    __tr04 = __tr04+__it;			__ti04 = __ti04-__rt;		\
										\
/*...Block 3: Swap inputs 6,7,8 -> 13,12,14: */\
	    __tr06 = __Ar13;				__ti06 = __Ai13;		\
	    __tr07 = __Ar12;				__ti07 = __Ai12;		\
	    __rt   = __Ar14;				__it   = __Ai14;		\
	    __tr08 = __tr07-__rt;			__ti08 = __ti07-__it;		\
	    __tr07 = __tr07+__rt;			__ti07 = __ti07+__it;		\
	    __tr06 = __tr06+__tr07;			__ti06 = __ti06+__ti07;		\
	    __tr07 = __tr06+c3m1*__tr07;	__ti07 = __ti06+c3m1*__ti07;	\
	    __rt   = s*__tr08;				__it   = s*__ti08;		\
	    __tr08 = __tr07-__it;			__ti08 = __ti07+__rt;		\
	    __tr07 = __tr07+__it;			__ti07 = __ti07-__rt;		\
										\
/*...Block 4: Swap inputs 9,10,11 -> 4,3,5: */\
	    __tr09 = __Ar04;				__ti09 = __Ai04;		\
	    __tr10 = __Ar03;				__ti10 = __Ai03;		\
	    __rt   = __Ar05;				__it   = __Ai05;		\
	    __tr11 = __tr10-__rt;			__ti11 = __ti10-__it;		\
	    __tr10 = __tr10+__rt;			__ti10 = __ti10+__it;		\
	    __tr09 = __tr09+__tr10;			__ti09 = __ti09+__ti10;		\
	    __tr10 = __tr09+c3m1*__tr10;	__ti10 = __ti09+c3m1*__ti10;	\
	    __rt   = s*__tr11;				__it   = s*__ti11;		\
	    __tr11 = __tr10-__it;			__ti11 = __ti10+__rt;		\
	    __tr10 = __tr10+__it;			__ti10 = __ti10-__rt;		\
										\
/*...Block 2: Swap inputs 12,13,14 -> 9,11,10: */\
	    __tr12 = __Ar09;				__ti12 = __Ai09;		\
	    __tr13 = __Ar11;				__ti13 = __Ai11;		\
	    __rt   = __Ar10;				__it   = __Ai10;		\
	    __tr14 = __tr13-__rt;			__ti14 = __ti13-__it;		\
	    __tr13 = __tr13+__rt;			__ti13 = __ti13+__it;		\
	    __tr12 = __tr12+__tr13;			__ti12 = __ti12+__ti13;		\
	    __tr13 = __tr12+c3m1*__tr13;	__ti13 = __ti12+c3m1*__ti13;	\
	    __rt   = s*__tr14;				__it   = s*__ti14;		\
	    __tr14 = __tr13-__it;			__ti14 = __ti13+__rt;		\
	    __tr13 = __tr13+__it;			__ti13 = __ti13-__rt;		\
										\
/*...and now do three radix-5 transforms.	\
	The required output permutation is [0,5,10,9,14,4,3,8,13,12,2,7,6,11,1]. */\
/*...Block 1: output permutation is 0,9,3,12,6 */\
	    __rt   = __tr12;				__it   = __ti12;		\
	    __tr12 = __tr03-__rt;			__ti12 = __ti03-__it;		\
	    __tr03 = __tr03+__rt;			__ti03 = __ti03+__it;		\
	    __rt   = __tr09;				__it   = __ti09;		\
	    __tr09 = __tr06-__rt;			__ti09 = __ti06-__it;		\
	    __tr06 = __tr06+__rt;			__ti06 = __ti06+__it;		\
										\
	    __rt   = __tr03+__tr06;			__it   = __ti03+__ti06;		\
	    __tr00 = __tr00+__rt;			__ti00 = __ti00+__it;		\
	    __rt   = __tr00+cn1*__rt;		__it   = __ti00+cn1*__it;	\
	    __tr06 = cn2*(__tr03-__tr06);	__ti06 = cn2*(__ti03-__ti06);	\
	    __tr03 = __rt+__tr06;			__ti03 = __it+__ti06;		\
	    __tr06 = __rt-__tr06;			__ti06 = __it-__ti06;		\
	    __rt   = ss3*(__tr09-__tr12);	__it   = ss3*(__ti09-__ti12);	\
	    __tr09 = __rt-sn1*__tr09;		__ti09 = __it-sn1*__ti09;	\
	    __tr12 = __rt+sn2*__tr12;		__ti12 = __it+sn2*__ti12;	\
										\
	    __Br00 = __tr00;				__Bi00 = __ti00;		\
	    __Br09 = __tr03-__ti09;			__Bi09 = __ti03+__tr09;		\
	    __Br03 = __tr06-__ti12;			__Bi03 = __ti06+__tr12;		\
	    __Br12 = __tr06+__ti12;			__Bi12 = __ti06-__tr12;		\
	    __Br06 = __tr03+__ti09;			__Bi06 = __ti03-__tr09;		\
										\
/*...Block 2: output permutation is 5,14,8,2,11 */\
	    __rt   = __tr13;				__it   = __ti13;		\
	    __tr13 = __tr04-__rt;			__ti13 = __ti04-__it;		\
	    __tr04 = __tr04+__rt;			__ti04 = __ti04+__it;		\
	    __rt   = __tr10;				__it   = __ti10;		\
	    __tr10 = __tr07-__rt;			__ti10 = __ti07-__it;		\
	    __tr07 = __tr07+__rt;			__ti07 = __ti07+__it;		\
										\
	    __rt   = __tr04+__tr07;			__it   = __ti04+__ti07;		\
	    __tr01 = __tr01+__rt;			__ti01 = __ti01+__it;		\
	    __rt   = __tr01+cn1*__rt;		__it   = __ti01+cn1*__it;	\
	    __tr07 = cn2*(__tr04-__tr07);	__ti07 = cn2*(__ti04-__ti07);	\
	    __tr04 = __rt+__tr07;			__ti04 = __it+__ti07;		\
	    __tr07 = __rt-__tr07;			__ti07 = __it-__ti07;		\
	    __rt   = ss3*(__tr10-__tr13);	__it   = ss3*(__ti10-__ti13);	\
	    __tr10 = __rt-sn1*__tr10;		__ti10 = __it-sn1*__ti10;	\
	    __tr13 = __rt+sn2*__tr13;		__ti13 = __it+sn2*__ti13;	\
										\
	    __Br05 = __tr01;				__Bi05 = __ti01;		\
	    __Br14 = __tr04-__ti10;			__Bi14 = __ti04+__tr10;		\
	    __Br08 = __tr07-__ti13;			__Bi08 = __ti07+__tr13;		\
	    __Br02 = __tr07+__ti13;			__Bi02 = __ti07-__tr13;		\
	    __Br11 = __tr04+__ti10;			__Bi11 = __ti04-__tr10;		\
										\
/*...Block 3: output permutation is 10,4,13,7,1 */\
	    __rt   = __tr14;				__it   = __ti14;		\
	    __tr14 = __tr05-__rt;			__ti14 = __ti05-__it;		\
	    __tr05 = __tr05+__rt;			__ti05 = __ti05+__it;		\
	    __rt   = __tr11;				__it   = __ti11;		\
	    __tr11 = __tr08-__rt;			__ti11 = __ti08-__it;		\
	    __tr08 = __tr08+__rt;			__ti08 = __ti08+__it;		\
										\
	    __rt   = __tr05+__tr08;			__it   = __ti05+__ti08;		\
	    __tr02 = __tr02+__rt;			__ti02 = __ti02+__it;		\
	    __rt   = __tr02+cn1*__rt;		__it   = __ti02+cn1*__it;	\
	    __tr08 = cn2*(__tr05-__tr08);	__ti08 = cn2*(__ti05-__ti08);	\
	    __tr05 = __rt+__tr08;			__ti05 = __it+__ti08;		\
	    __tr08 = __rt-__tr08;			__ti08 = __it-__ti08;		\
	    __rt   = ss3*(__tr11-__tr14);	__it   = ss3*(__ti11-__ti14);	\
	    __tr11 = __rt-sn1*__tr11;		__ti11 = __it-sn1*__ti11;	\
	    __tr14 = __rt+sn2*__tr14;		__ti14 = __it+sn2*__ti14;	\
										\
	    __Br10 = __tr02;				__Bi10 = __ti02;		\
	    __Br04 = __tr05-__ti11;			__Bi04 = __ti05+__tr11;		\
	    __Br13 = __tr08-__ti14;			__Bi13 = __ti08+__tr14;		\
	    __Br07 = __tr08+__ti14;			__Bi07 = __ti08-__tr14;		\
	    __Br01 = __tr05+__ti11;			__Bi01 = __ti05-__tr11;		\
}

/* Totals: 144 FADD, 24 FMUL	*/
#define RADIX_16_DIF(\
	__A0r,__A0i,\
	__A1r,__A1i,\
	__A2r,__A2i,\
	__A3r,__A3i,\
	__A4r,__A4i,\
	__A5r,__A5i,\
	__A6r,__A6i,\
	__A7r,__A7i,\
	__A8r,__A8i,\
	__A9r,__A9i,\
	__AAr,__AAi,\
	__ABr,__ABi,\
	__ACr,__ACi,\
	__ADr,__ADi,\
	__AEr,__AEi,\
	__AFr,__AFi,\
	__t1,__t2,__t3,__t4,__t5,__t6,__t7,__t8,__t9,__t10,__t11,__t12,__t13,__t14,__t15,__t16,__t17,__t18,__t19,__t20,__t21,__t22,__t23,__t24,__t25,__t26,__t27,__t28,__t29,__t30,__t31,__t32,\
	__B0r,__B0i,\
	__B1r,__B1i,\
	__B2r,__B2i,\
	__B3r,__B3i,\
	__B4r,__B4i,\
	__B5r,__B5i,\
	__B6r,__B6i,\
	__B7r,__B7i,\
	__B8r,__B8i,\
	__B9r,__B9i,\
	__BAr,__BAi,\
	__BBr,__BBi,\
	__BCr,__BCi,\
	__BDr,__BDi,\
	__BEr,__BEi,\
	__BFr,__BFi,\
	__rt,__it,\
	__c,__s)\
{\
	/*       gather the needed data (16 64-bit complex, i.e. 32 64-bit reals) and do the first set of four length-4 transforms.	*/\
	/*...Block 1:	*/\
	__t3 =__A0r -__A8r;	__t1 =__A0r +__A8r;\
	__t4 =__A0i -__A8i;	__t2 =__A0i +__A8i;\
	\
	__t7 =__A4r -__ACr;	__t5 =__A4r +__ACr;\
	__t8 =__A4i -__ACi;	__t6 =__A4i +__ACi;\
	\
	__rt =__t5;		__t5 =__t1 -__rt;	__t1 =__t1 +__rt;\
	__it =__t6;		__t6 =__t2 -__it;	__t2 =__t2 +__it;\
	\
	__rt =__t7;		__t7 =__t3 +__t8;	__t3 =__t3 -__t8;\
					__t8 =__t4 -__rt;	__t4 =__t4 +__rt;\
	/*...Block 2:	*/\
	__t11=__A2r -__AAr;	__t9 =__A2r +__AAr;\
	__t12=__A2i -__AAi;	__t10=__A2i +__AAi;\
	\
	__t15=__A6r -__AEr;	__t13=__A6r +__AEr;\
	__t16=__A6i -__AEi;	__t14=__A6i +__AEi;\
	\
	__rt =__t13;	__t13=__t9 -__rt;	__t9 =__t9 +__rt;\
	__it =__t14;	__t14=__t10-__it;	__t10=__t10+__it;\
	\
	__rt =__t15;	__t15=__t11+__t16;__t11=__t11-__t16;\
					__t16=__t12-__rt;	__t12=__t12+__rt;\
	/*...Block 3:	*/\
	__t19=__A1r -__A9r; __t17=__A1r +__A9r;\
	__t20=__A1i -__A9i;	__t18=__A1i +__A9i;\
	\
	__t23=__A5r -__ADr;	__t21=__A5r +__ADr;\
	__t24=__A5i -__ADi;	__t22=__A5i +__ADi;\
	\
	__rt =__t21;	__t21=__t17-__rt;	__t17=__t17+__rt;\
	__it =__t22;	__t22=__t18-__it;	__t18=__t18+__it;\
	\
	__rt =__t23;	__t23=__t19+__t24;	__t19=__t19-__t24;\
					__t24=__t20-__rt;	__t20=__t20+__rt;\
	/*...Block 4:	*/\
	__t27=__A3r -__ABr;	__t25=__A3r +__ABr;\
	__t28=__A3i -__ABi;	__t26=__A3i +__ABi;\
	\
	__t31=__A7r -__AFr;	__t29=__A7r +__AFr;\
	__t32=__A7i -__AFi;	__t30=__A7i +__AFi;\
	\
	__rt =__t29;	__t29=__t25-__rt;	__t25=__t25+__rt;\
	__it =__t30;	__t30=__t26-__it;	__t26=__t26+__it;\
	\
	__rt =__t31;	__t31=__t27+__t32;__t27=__t27-__t32;\
					__t32=__t28-__rt;	__t28=__t28+__rt;\
	\
	/*...and now do four more radix-4 __transforms, including the internal twiddle factors: */\
	\
	/*...Block 1: __t1,9,17,25	*/\
	__rt =__t9 ;	__t9 =__t1 -__rt;	__t1 =__t1 +__rt;\
	__it =__t10;	__t10=__t2 -__it;	__t2 =__t2 +__it;\
	\
	__rt =__t25;	__t25=__t17-__rt;	__t17=__t17+__rt;\
	__it =__t26;	__t26=__t18-__it;	__t18=__t18+__it;\
	\
	__B0r =__t1+__t17;			__B0i =__t2+__t18;\
	__B1r =__t1-__t17;			__B1i =__t2-__t18;\
	\
	__B2r =__t9 -__t26;			__B2i =__t10+__t25;			/* mpy by E^4=i is inlined here...	*/\
	__B3r =__t9 +__t26;			__B3i =__t10-__t25;\
	\
	/*...Block 3: __t5,13,21,29	*/\
	__rt =__t13;	__t13=__t5 +__t14;	__t5 =__t5 -__t14;	/* twiddle mpy by E^4 = I	*/\
					__t14=__t6 -__rt;	__t6 =__t6 +__rt;\
	\
	__rt =(__t21-__t22)*ISRT2;	__t22=(__t21+__t22)*ISRT2;	__t21=__rt;	/* twiddle mpy by E^2; could also do mpy by ISRT2 directly on t21,22,29,30	*/\
	__rt =(__t30+__t29)*ISRT2;	__it =(__t30-__t29)*ISRT2;	/* twiddle mpy by -E^6 is here...	*/\
	__t29=__t21+__rt;			__t21=__t21-__rt;			/* ...and get E^6=(i-1)/sq__rt by flipping signs here.	*/\
	__t30=__t22+__it;			__t22=__t22-__it;\
	\
	__B4r =__t5+__t21;			__B4i =__t6+__t22;\
	__B5r =__t5-__t21;			__B5i =__t6-__t22;\
	\
	__B6r =__t13-__t30;			__B6i =__t14+__t29;			/* mpy by E^4=i is inlined here...	*/\
	__B7r =__t13+__t30;			__B7i =__t14-__t29;\
	\
	/*...Block 2: __t3,11,19,27	*/\
	__rt =(__t11-__t12)*ISRT2;	__it =(__t11+__t12)*ISRT2;	/* twiddle mpy by E^2; could also do mpy by ISRT2 directly on t11,12	*/\
	__t11=__t3 -__rt;			__t3 =__t3 +__rt;\
	__t12=__t4 -__it;			__t4 =__t4 +__it;\
	\
	__rt =__t19*__c - __t20*__s;	__t20=__t20*__c + __t19*__s;	__t19=__rt;	/* twiddle mpy by E^1	*/\
	__rt =__t27*__s - __t28*__c;	__it =__t28*__s + __t27*__c;	/* twiddle mpy by E^3	*/\
	__t27=__t19-__rt;			__t19=__t19+__rt;\
	__t28=__t20-__it;			__t20=__t20+__it;\
	\
	__B8r =__t3+__t19;			__B8i =__t4+__t20;\
	__B9r =__t3-__t19;			__B9i =__t4-__t20;\
	\
	__BAr =__t11-__t28;			__BAi =__t12+__t27;			/* mpy by E^4=i is inlined here...	*/\
	__BBr =__t11+__t28;			__BBi =__t12-__t27;\
	\
	/*...Block 4: __t7,15,23,31	*/\
	__rt =(__t16+__t15)*ISRT2;	__it =(__t16-__t15)*ISRT2;	/* twiddle mpy by -E^6 is here...	*/\
	__t15=__t7 +__rt;			__t7 =__t7 -__rt;			/* ...and get E^6=(i-1)/sq__rt by flipping signs here.	*/\
	__t16=__t8 +__it;			__t8 =__t8 -__it;\
	\
	__rt =__t23*__s - __t24*__c;	__t24=__t24*__s + __t23*__c;	__t23=__rt;	/* twiddle mpy by E^3	*/\
	__rt =__t31*__c - __t32*__s;	__it =__t32*__c + __t31*__s;	/* twiddle mpy by E^1 = -E^9...	*/\
	__t31=__t23+__rt;			__t23=__t23-__rt;			/* ...and get E^9 by flipping signs here.	*/\
	__t32=__t24+__it;			__t24=__t24-__it;			/* Note: t23+__rt = t23*(s+1)	*/\
	\
	__BCr =__t7+__t23;			__BCi =__t8+__t24;\
	__BDr =__t7-__t23;			__BDi =__t8-__t24;\
	\
	__BEr =__t15-__t32;			__BEi =__t16+__t31;			/* mpy by E^4=i is inlined here...	*/\
	__BFr =__t15+__t32;			__BFi =__t16-__t31;\
}

/* Totals: 144 FADD, 24 FMUL	*/
#define RADIX_16_DIT(\
	__A0r,__A0i,\
	__A1r,__A1i,\
	__A2r,__A2i,\
	__A3r,__A3i,\
	__A4r,__A4i,\
	__A5r,__A5i,\
	__A6r,__A6i,\
	__A7r,__A7i,\
	__A8r,__A8i,\
	__A9r,__A9i,\
	__AAr,__AAi,\
	__ABr,__ABi,\
	__ACr,__ACi,\
	__ADr,__ADi,\
	__AEr,__AEi,\
	__AFr,__AFi,\
	__t1,__t2,__t3,__t4,__t5,__t6,__t7,__t8,__t9,__t10,__t11,__t12,__t13,__t14,__t15,__t16,__t17,__t18,__t19,__t20,__t21,__t22,__t23,__t24,__t25,__t26,__t27,__t28,__t29,__t30,__t31,__t32,\
	__B0r,__B0i,\
	__B1r,__B1i,\
	__B2r,__B2i,\
	__B3r,__B3i,\
	__B4r,__B4i,\
	__B5r,__B5i,\
	__B6r,__B6i,\
	__B7r,__B7i,\
	__B8r,__B8i,\
	__B9r,__B9i,\
	__BAr,__BAi,\
	__BBr,__BBi,\
	__BCr,__BCi,\
	__BDr,__BDi,\
	__BEr,__BEi,\
	__BFr,__BFi,\
	__rt,__it,\
	__c,__s)\
{\
	/*...Block 1:	*/\
	__t1 =__A0r;					__t2 =__A0i ;\
	__rt =__A1r;					__it =__A1i ;\
	__t3 =__t1 -__rt;  				__t1 =__t1 +__rt;\
	__t4 =__t2 -__it;				__t2 =__t2 +__it;\
	\
	__t5 =__A2r;					__t6 =__A2i ;\
	__rt =__A3r;					__it =__A3i ;\
	__t7 =__t5 -__rt;  				__t5 =__t5 +__rt;\
	__t8 =__t6 -__it;  				__t6 =__t6 +__it;\
	\
	__rt =__t5;	__t5 =__t1 -__rt;	__t1 =__t1 +__rt;\
	__it =__t6;	__t6 =__t2 -__it;	__t2 =__t2 +__it;\
	\
	__rt =__t7;	__t7 =__t3 -__t8;	__t3 =__t3 +__t8;\
				__t8 =__t4 +__rt;	__t4 =__t4 -__rt;\
	\
	/*...Block 2:	*/\
	__t9 =__A4r;					__t10=__A4i ;\
	__rt =__A5r;					__it =__A5i ;\
	__t11=__t9 -__rt;  				__t9 =__t9 +__rt;\
	__t12=__t10-__it;				__t10=__t10+__it;\
	\
	__t13=__A6r;					__t14=__A6i ;\
	__rt =__A7r;					__it =__A7i ;\
	__t15=__t13-__rt;  				__t13=__t13+__rt;\
	__t16=__t14-__it;				__t14=__t14+__it;\
	\
	__rt =__t13;__t13=__t9 -__rt;	__t9 =__t9 +__rt;\
	__it =__t14;__t14=__t10-__it;	__t10=__t10+__it;\
	\
	__rt =__t15;__t15=__t11-__t16;	__t11=__t11+__t16;\
				__t16=__t12+__rt;	__t12=__t12-__rt;\
	\
	/*...Block 3:	*/\
	__t17=__A8r;					__t18=__A8i ;\
	__rt =__A9r;					__it =__A9i ;\
	__t19=__t17-__rt;  				__t17=__t17+__rt;\
	__t20=__t18-__it;				__t18=__t18+__it;\
	\
	__t21=__AAr;					__t22=__AAi ;\
	__rt =__ABr;					__it =__ABi ;\
	__t23=__t21-__rt;  				__t21=__t21+__rt;\
	__t24=__t22-__it;				__t22=__t22+__it;\
	\
	__rt =__t21;__t21=__t17-__rt;	__t17=__t17+__rt;\
	__it =__t22;__t22=__t18-__it;	__t18=__t18+__it;\
	\
	__rt =__t23;__t23=__t19-__t24;	__t19=__t19+__t24;\
				__t24=__t20+__rt;	__t20=__t20-__rt;\
	\
	/*...Block 4:	*/\
	__t25=__ACr;					__t26=__ACi ;\
	__rt =__ADr;					__it =__ADi ;\
	__t27=__t25-__rt;  				__t25=__t25+__rt; \
	__t28=__t26-__it;				__t26=__t26+__it;	\
	\
	__t29=__AEr;					__t30=__AEi ;\
	__rt =__AFr;					__it =__AFi ;\
	__t31=__t29-__rt;  				__t29=__t29+__rt;\
	__t32=__t30-__it;				__t30=__t30+__it;\
	\
	__rt =__t29;__t29=__t25-__rt;	__t25=__t25+__rt;\
	__it =__t30;__t30=__t26-__it;	__t26=__t26+__it;\
	\
	__rt =__t31;__t31=__t27-__t32;	__t27=__t27+__t32;\
				__t32=__t28+__rt;	__t28=__t28-__rt;\
	\
	/*...and now do four more radix-4 __transforms, including the internal twiddle factors: */\
	\
	/*...Block 1: __t1,9,17,25	*/\
	__rt =__t9;	__t9 =__t1 -__rt;	__t1 =__t1 +__rt;\
	__it =__t10;__t10=__t2 -__it;	__t2 =__t2 +__it;\
	\
	__rt =__t25;__t25=__t17-__rt;	__t17=__t17+__rt;\
	__it =__t26;__t26=__t18-__it;	__t18=__t18+__it;\
	\
	__B0r=__t1+__t17;			__B0i =__t2+__t18;\
	__B8r=__t1-__t17;			__B8i =__t2-__t18;\
	\
	__B4r=__t9 +__t26;			__B4i =__t10-__t25;	/* mpy by E^-4 = -I is inlined here...	*/\
	__BCr=__t9 -__t26;			__BCi =__t10+__t25;\
	\
	/*...Block 3: __t5,13,21,29	*/\
	__rt =__t13;__t13=__t5 -__t14;	__t5 =__t5 +__t14;	/* __twiddle mpy by E^4 =-I	*/\
				__t14=__t6 +__rt;	__t6 =__t6 -__rt;\
	\
	__rt =(__t22+__t21)*ISRT2;	__t22=(__t22-__t21)*ISRT2;	__t21=__rt;	/* twiddle mpy by E^-2	*/\
	__rt =(__t29-__t30)*ISRT2;	__it =(__t29+__t30)*ISRT2;		/* twiddle mpy by E^2 = -E^-6 is here...	*/\
	__t29=__t21+__rt;			__t21=__t21-__rt;				/* ...and get E^-6 by flipping signs here.	*/\
	__t30=__t22+__it;			__t22=__t22-__it;\
	\
	__B2r=__t5+__t21;			__B2i =__t6+__t22;\
	__BAr=__t5-__t21;			__BAi =__t6-__t22;\
	\
	__B6r=__t13+__t30;			__B6i =__t14-__t29;	/* mpy by E^-4 =-I is inlined here...	*/\
	__BEr=__t13-__t30;			__BEi =__t14+__t29;\
	\
	/*...Block 2: __t3,11,19,27	*/\
	__rt =(__t12+__t11)*ISRT2;	__it =(__t12-__t11)*ISRT2;	/* twiddle mpy by E^-2	*/\
	__t11=__t3 -__rt;			__t3 =__t3 +__rt;\
	__t12=__t4 -__it;			__t4 =__t4 +__it;\
	\
	__rt =__t19*__c + __t20*__s;	__t20=__t20*__c - __t19*__s;	__t19=__rt;	/* twiddle mpy by E^-1	*/\
	__rt =__t27*__s + __t28*__c;	__it =__t28*__s - __t27*__c;	/* twiddle mpy by E^-3	*/\
	__t27=__t19-__rt;			__t19=__t19+__rt;\
	__t28=__t20-__it;			__t20=__t20+__it;\
	\
	__B1r=__t3+__t19;			__B1i =__t4+__t20;\
	__B9r=__t3-__t19;			__B9i =__t4-__t20;\
	\
	__B5r=__t11+__t28;			__B5i =__t12-__t27;	/* mpy by E^-4 =-I is inlined here...	*/\
	__BDr=__t11-__t28;			__BDi =__t12+__t27;\
	\
	/*...Block 4: __t7,15,23,31	*/\
	__rt =(__t15-__t16)*ISRT2;	__it =(__t15+__t16)*ISRT2;	/* twiddle mpy by E^2 = -E^-6 is here...	*/\
	__t15=__t7 +__rt;			__t7 =__t7 -__rt;			/* ...and get E^6=(i-1)/sqrt by flipping signs here.	*/\
	__t16=__t8 +__it;			__t8 =__t8 -__it;\
	\
	__rt =__t23*__s + __t24*__c;	__t24=__t24*__s - __t23*__c;	__t23=__rt;	/* twiddle mpy by E^-3	*/\
	__rt =__t31*__c + __t32*__s;	__it =__t32*__c - __t31*__s;	/* twiddle mpy by E^-1 = -E^-9...	*/\
	__t31=__t23+__rt;			__t23=__t23-__rt;			/* ...and get E^9 by flipping signs here.	*/\
	__t32=__t24+__it;			__t24=__t24-__it;\
	\
	__B3r=__t7+__t23;			__B3i =__t8+__t24;\
	__BBr=__t7-__t23;			__BBi =__t8-__t24;\
	\
	__B7r=__t15+__t32;			__B7i =__t16-__t31;	/* mpy by E^-4 = -I is inlined here...	*/\
	__BFr=__t15-__t32;			__BFi =__t16+__t31;\
}

/************** RADIX-32 DIF/DIT: *****************************/
/* Totals: 376 FADD, 88 FMUL	*/
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
	/* Each complex radix-8 subtransform needs 52 ADD, 4 MUL (not counting load/store/address-compute):	*/\
	/*...Block 1:	*/\
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
	/*...Block 3:	*/\
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
	/*...Block 4:	*/\
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
	/*...and now do eight radix-4 transforms, including the internal twiddle factors:	*/\
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
	/* Each complex radix-8 subtransform needs 52 ADD, 4 MUL (not counting load/store/address-compute):	*/\
	/*...Block 1:	*/\
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
	/*...Block 3:	*/\
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
	/*...Block 4:	*/\
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
	/*...and now do eight radix-4 transforms, including the internal twiddle factors:	*/\
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
