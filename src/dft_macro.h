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
#ifndef dft_macro_included
#define dft_macro_included

#if 0
	NOTE: To make macros here easier to browse in a code-smart editor, temp-uncomment the placeholder-function lines:
	//void ___RADIX...
#endif

/* All these radix-N macros are structured so that any 2 of the 3 sets of arguments (A, t, B) may point to the same set of
   addresses, i.e. so the DFT may be done in place (A == B) or re-using a single set of temporaries (A == t or B == t).
   By "in place" we mean that A == B up to an arbitrary permutation of array indices (or element subscripts). */

/****** RADIX = 3: ALLOWS IN-PLACE ******/

/* Totals: 12 ADD, 4 MUL	*/
//void ___RADIX_03_DFT() {}	// placeholder to get macro-of-same-name locations to appear in my editor's function-display
#define RADIX_03_DFT(\
	__s,__c3m1,\
	__Ar0,__Ai0,__Ar1,__Ai1,__Ar2,__Ai2,\
	__tr0,__ti0,__tr1,__ti1,__tr2,__ti2,\
	__Br0,__Bi0,__Br1,__Bi1,__Br2,__Bi2)\
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
//void ___RADIX_03_DFT_PFETCH() {}
#define RADIX_03_DFT_PFETCH(\
	__s,__c3m1,\
	__Ar0,__Ai0,__Ar1,__Ai1,__Ar2,__Ai2,\
	__tr0,__ti0,__tr1,__ti1,__tr2,__ti2,\
	__Br0,__Bi0,__Br1,__Bi1,__Br2,__Bi2,\
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
//void ___RADIX_04_DIF() {}
#define RADIX_04_DIF(\
	__Ar0,__Ai0,__Ar1,__Ai1,__Ar2,__Ai2,__Ar3,__Ai3,\
	__Br0,__Bi0,__Br1,__Bi1,__Br2,__Bi2,__Br3,__Bi3,\
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
//void ___RADIX_04_DIF_PFETCH() {}
#define RADIX_04_DIF_PFETCH(\
	__Ar0,__Ai0,__Ar1,__Ai1,__Ar2,__Ai2,__Ar3,__Ai3,\
	__Br0,__Bi0,__Br1,__Bi1,__Br2,__Bi2,__Br3,__Bi3,\
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

/* Totals: 22 ADD, 12 MUL, which is the cost of the twiddle-less version plus 3 CMUL at [2 ADD, 4 MUL] each. Combine
8 such macro calls (2 passes of 4 calls each) to create an alternate uniform-macro version of RADIX_16_DIF_TWIDDLE,
with Cost [176 ADD, 96 MUL], 2 more ADD and 12 more MUL, but represents a smaller optimization target, and is thus
ideal for SIMD inline-ASM implementation.
*/
//void ___RADIX_04_DIF() {}
#define RADIX_04_DIF_3TWIDDLE(\
	__Ar0,__Ai0,__Ar1,__Ai1,__Ar2,__Ai2,__Ar3,__Ai3,\
	__Br0,__Bi0,__Br1,__Bi1,__Br2,__Bi2,__Br3,__Bi3,\
	__c1,__s1 ,__c2,__s2 ,__c3,__s3,\
	__r,__i)\
{\
	__Br0 = __Ar0;						__Bi0 = __Ai0;\
	__r   = __Ar2*__c2 - __Ai2*__s2;	__i   = __Ai2*__c2 + __Ar2*__s2 ;\
	__Br1 = __Br0 - __r;				__Br0 = __Br0 + __r;\
	__Bi1 = __Bi0 - __i;				__Bi0 = __Bi0 + __i;\
\
	__Br2 = __Ar1*__c1 - __Ai1*__s1;	__Bi2 = __Ai1*__c1 + __Ar1*__s1 ;\
	__r   = __Ar3*__c3 - __Ai3*__s3;	__i   = __Ai3*__c3 + __Ar3*__s3 ;\
	__Br3 = __Br2 - __r;				__Br2 = __Br2 + __r;\
	__Bi3 = __Bi2 - __i;				__Bi2 = __Bi2 + __i;\
\
	__r = __Br2;	__Br2 = __Br0 - __r;		__Br0 = __Br0 + __r;\
	__i = __Bi2;	__Bi2 = __Bi0 - __i;		__Bi0 = __Bi0 + __i;\
\
	__r = __Br3;	__Br3 = __Br1 + __Bi3;		__Br1 = __Br1 - __Bi3;\
					__Bi3 = __Bi1 - __r;		__Bi1 = __Bi1 + __r;\
}

/* Totals: 16 ADD, 0 MUL	*/
//void ___RADIX_04_DIT() {}
#define RADIX_04_DIT(\
	__Ar0,__Ai0,__Ar1,__Ai1,__Ar2,__Ai2,__Ar3,__Ai3,\
	__Br0,__Bi0,__Br1,__Bi1,__Br2,__Bi2,__Br3,__Bi3,\
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

/* Totals: 22 ADD, 12 MUL, which is the cost of the twiddle-less version plus 3 CMUL at [2 ADD, 4 MUL] each. Combine
8 such macro calls (2 passes of 4 calls each) to create an alternate uniform-macro version of RADIX_16_DIT_TWIDDLE,
with Cost [176 ADD, 96 MUL], 2 more ADD and 12 more MUL, but represents a smaller optimization target, and is thus
ideal for SIMD inline-ASM implementation.
Notes:
o DIF version of this uses pretwiddles; DIT uses posttwidles;
o DIT is geared toward inverse FFT but with roots possibly shared with DIF, thus assumes roots are
  those with a positive-signed complex exponential argument, and conjugates them in-place for DIT.
*/
//void ___RADIX_04_DIF() {}
#define RADIX_04_DIT_3TWIDDLE(\
	__Ar0,__Ai0,__Ar1,__Ai1,__Ar2,__Ai2,__Ar3,__Ai3,\
	__Br0,__Bi0,__Br1,__Bi1,__Br2,__Bi2,__Br3,__Bi3,\
	__c1,__s1 ,__c2,__s2 ,__c3,__s3,\
	__r,__i)\
{\
	double _tr0,_ti0,_tr1,_ti1,_tr2,_ti2,_tr3,_ti3;\
	_tr1 = __Ar0 - __Ar1;			_tr0 = __Ar0 + __Ar1;\
	_ti1 = __Ai0 - __Ai1;			_ti0 = __Ai0 + __Ai1;\
\
	_tr3 = __Ar2 - __Ar3;			_tr2 = __Ar2 + __Ar3;\
	_ti3 = __Ai2 - __Ai3;			_ti2 = __Ai2 + __Ai3;\
\
	__Br0 = _tr0 + _tr2;			__Bi0 = _ti0 + _ti2;\
	 _tr0 = _tr0 - _tr2;			 _ti0 = _ti0 - _ti2;\
	__Br2 = _tr0*__c2 + _ti0*__s2;	__Bi2 = _ti0*__c2 - _tr0*__s2;	/* twiddle = ~w2 = c2-I.s1 */\
	/* mpy by -I is inlined here...	*/\
	__r = _tr3;	_tr3 = _tr1 - _ti3;	_tr1 = _tr1 + _ti3;\
				_ti3 = _ti1 + __r;	_ti1 = _ti1 - __r;\
	__Br1 = _tr1*__c1 + _ti1*__s1;	__Bi1 = _ti1*__c1 - _tr1*__s1;	/* twiddle = ~w1 = c1-I.s1 */\
	__Br3 = _tr3*__c3 + _ti3*__s3;	__Bi3 = _ti3*__c3 - _tr3*__s3;	/* twiddle = ~w3 = c3-I.s3 */\
}

/****** RADIX = 5: ******/

/* Totals: 34 ADD, 10 MUL	*/
//void ___RADIX_05_DFT() {}
#define RADIX_05_DFT(\
	__cc1, __cc2, __s2, __ss1, __ss2,\
	__Ar0,__Ai0,__Ar1,__Ai1,__Ar2,__Ai2,__Ar3,__Ai3,__Ar4,__Ai4,\
	__Br0,__Bi0,__Br1,__Bi1,__Br2,__Bi2,__Br3,__Bi3,__Br4,__Bi4,\
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
//void ___RADIX_05_DFT_PFETCH() {}
#define RADIX_05_DFT_PFETCH(\
	__cc1, __cc2, __s2, __ss1, __ss2,\
	__Ar0,__Ai0,__Ar1,__Ai1,__Ar2,__Ai2,__Ar3,__Ai3,__Ar4,__Ai4,\
	__Br0,__Bi0,__Br1,__Bi1,__Br2,__Bi2,__Br3,__Bi3,__Br4,__Bi4,\
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

/* Simple version may be best if ADD throughput is limiting factor. Totals: 60 ADD, 36 MUL	*/
/*
Base-root is E := exp(I*2*Pi/7), with symmetries E^(n-j) = ~E^j and E^n = E^(n mod 7).
Thus our DFT matrix is

		1	1	1	1	|	1	1	1
		----------------|--------------
		1	E1	E2	E3	|	E3~	E2~	E1~
		1	E2	E3~	E1~	|	E1	E3	E2~
		1	E3	E1~	E2	|	E2~	E1	E3~
		----------------|--------------
		1	E3~	E1	E2~	|	E2	E1~	E3
		1	E2~	E3	E1	|	E1~	E3~	E2
		1	E1~	E2~	E3~	|	E3	E2	E1 ;

the upper left 3 x 3 E-submatrix gives the form of our mini-convolutions, with ~ terms having negative
sine multipliers, i.e. Ej~ ==> -sin(2*Pi*j/7) .
*/
//void ___RADIX_07_DFT() {}
#define RADIX_07_DFT(\
	__Ar0,__Ai0,__Ar1,__Ai1,__Ar2,__Ai2,__Ar3,__Ai3,__Ar4,__Ai4,__Ar5,__Ai5,__Ar6,__Ai6,\
	__tr0,__ti0,__tr1,__ti1,__tr2,__ti2,__tr3,__ti3,__tr4,__ti4,__tr5,__ti5,__tr6,__ti6,\
	__Br0,__Bi0,__Br1,__Bi1,__Br2,__Bi2,__Br3,__Bi3,__Br4,__Bi4,__Br5,__Bi5,__Br6,__Bi6,\
	__cc1,__ss1,__cc2,__ss2,__cc3,__ss3,\
	__rt,__it,__re,__im)\
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

/* FMA-oriented version - above simple-structure better for FMA than Nussbaumer:
Totals: 30 ADD, 6 MUL, 30 FMA, i.e. trade 30 ADD + 30 MUL for 30 FMA, very nice.
We can also use FMA to save registers in computing radix-2 butterflies x +- y,
at a cost of replacing 12 ADD ==> 12 FMA, for a total of [18 ADD, 6 MUL, 42 FMA].
(This is what we do in our AVX2 implementation in sse_macro_gcc64.h:SSE2_RADIX_07_DFT().)
*/
//void ___RADIX_07_DFT_FMA() {}
#define RADIX_07_DFT_FMA(\
	__Ar0,__Ai0,__Ar1,__Ai1,__Ar2,__Ai2,__Ar3,__Ai3,__Ar4,__Ai4,__Ar5,__Ai5,__Ar6,__Ai6,\
	__tr0,__ti0,__tr1,__ti1,__tr2,__ti2,__tr3,__ti3,__tr4,__ti4,__tr5,__ti5,__tr6,__ti6,\
	__Br0,__Bi0,__Br1,__Bi1,__Br2,__Bi2,__Br3,__Bi3,__Br4,__Bi4,__Br5,__Bi5,__Br6,__Bi6,\
	__cc1,__ss1,__cc2,__ss2,__cc3,__ss3,\
	__rt,__it,__re,__im)\
{\
	__tr0 = __Ar0;									__ti0 = __Ai0;\
	__tr6 = __Ar1 - __Ar6;							__ti6 = __Ai1 - __Ai6;	/* x1 - x6 */\
	__tr5 = __Ar2 - __Ar5;							__ti5 = __Ai2 - __Ai5;	/* x2 - x5 */\
	__tr4 = __Ar3 - __Ar4;							__ti4 = __Ai3 - __Ai4;	/* x3 - x4 */\
	__tr1 = __Ar1 + __Ar6;							__ti1 = __Ai1 + __Ai6;	/* x1 + x6 */\
	__tr2 = __Ar2 + __Ar5;							__ti2 = __Ai2 + __Ai5;	/* x2 + x5 */\
	__tr3 = __Ar3 + __Ar4;							__ti3 = __Ai3 + __Ai4;	/* x3 + x4 */\
\
/* Structure this next FMA-heavy sequence with a view toward the Intel-Haswell issue-rate/	\
latency model, i.e. attempting to mitigate a 5-cycle latency, with 2 FMAs issuing per cycle.\
The Broadwell reduced 3-cycle latency would make things easy, but 5-cycle means doing stuff	\
like interleaving the accumulate-ADDs for the DC outputs Br0,Bi0 with the FMAs: */\
	__Br0 = __tr0;									__Bi0 = __ti0;\
/*	__rt  = __tr0 + __cc1*__tr1 + __cc2*__tr2 + __cc3*__tr3;	__it  = __ti0 + __cc1*__ti1 + __cc2*__ti2 + __cc3*__ti3;*/\
/*	__re  = __tr0 + __cc2*__tr1 + __cc3*__tr2 + __cc1*__tr3;	__im  = __ti0 + __cc2*__ti1 + __cc3*__ti2 + __cc1*__ti3;*/\
/*	__tr0 = __tr0 + __cc3*__tr1 + __cc1*__tr2 + __cc2*__tr3;	__ti0 = __ti0 + __cc3*__ti1 + __cc1*__ti2 + __cc2*__ti3;*/\
	__rt  = __re = __tr0;							__it  = __im = __ti0;\
	__rt  = __FMADD(__cc1,__tr1, __rt );			__it  = __FMADD(__cc1,__ti1, __it );\
	__re  = __FMADD(__cc2,__tr1, __re );			__im  = __FMADD(__cc2,__ti1, __im );\
	__tr0 = __FMADD(__cc3,__tr1, __tr0);			__ti0 = __FMADD(__cc3,__ti1, __ti0);\
	__Br0 += __tr1;									__Bi0 += __ti1;\
	__rt  = __FMADD(__cc2,__tr2, __rt );			__it  = __FMADD(__cc2,__ti2, __it );\
	__re  = __FMADD(__cc3,__tr2, __re );			__im  = __FMADD(__cc3,__ti2, __im );\
	__tr0 = __FMADD(__cc1,__tr2, __tr0);			__ti0 = __FMADD(__cc1,__ti2, __ti0);\
	__Br0 += __tr2;									__Bi0 += __ti2;\
	__rt  = __FMADD(__cc3,__tr3, __rt );			__it  = __FMADD(__cc3,__ti3, __it );\
	__re  = __FMADD(__cc1,__tr3, __re );			__im  = __FMADD(__cc1,__ti3, __im );\
	__tr0 = __FMADD(__cc2,__tr3, __tr0);			__ti0 = __FMADD(__cc2,__ti3, __ti0);\
	__Br0 += __tr3;									__Bi0 += __ti3;\
\
/* For the Im-parts we have no DC accumulation o interleave with the FMAs, alas: */\
/*	__tr1 = __ss1*__tr6 + __ss2*__tr5 + __ss3*__tr4;		__ti1 = __ss1*__ti6 + __ss2*__ti5 + __ss3*__ti4;*/\
/*	__tr2 = __ss2*__tr6 - __ss3*__tr5 - __ss1*__tr4;		__ti2 = __ss2*__ti6 - __ss3*__ti5 - __ss1*__ti4;*/\
/*	__tr3 = __ss3*__tr6 - __ss1*__tr5 + __ss2*__tr4;		__ti3 = __ss3*__ti6 - __ss1*__ti5 + __ss2*__ti4;*/\
	__tr1 = __ss1*__tr6;							__ti1 = __ss1*__ti6;\
	__tr2 = __ss2*__tr6;							__ti2 = __ss2*__ti6;\
	__tr3 = __ss3*__tr6;							__ti3 = __ss3*__ti6;\
\
	__tr1 =  __FMADD(__ss2,__tr5, __tr1);			__ti1 =  __FMADD(__ss2,__ti5, __ti1);\
	__tr2 = __FNMADD(__ss3,__tr5, __tr2);			__ti2 = __FNMADD(__ss3,__ti5, __ti2);\
	__tr3 = __FNMADD(__ss1,__tr5, __tr3);			__ti3 = __FNMADD(__ss1,__ti5, __ti3);\
\
	__tr1 =  __FMADD(__ss3,__tr4, __tr1);			__ti1 =  __FMADD(__ss3,__ti4, __ti1);\
	__tr2 = __FNMADD(__ss1,__tr4, __tr2);			__ti2 = __FNMADD(__ss1,__ti4, __ti2);\
	__tr3 =  __FMADD(__ss2,__tr4, __tr3);			__ti3 =  __FMADD(__ss2,__ti4, __ti3);\
	/* Output permutation causes signs to get flipped here: */\
	__Br1 = __rt  - __ti1;							__Bi1 = __it  + __tr1;\
	__Br2 = __re  - __ti2;							__Bi2 = __im  + __tr2;\
	__Br3 = __tr0 - __ti3;							__Bi3 = __ti0 + __tr3;\
	__Br4 = __tr0 + __ti3;							__Bi4 = __ti0 - __tr3;\
	__Br5 = __re  + __ti2;							__Bi5 = __im  - __tr2;\
	__Br6 = __rt  + __ti1;							__Bi6 = __it  - __tr1;\
}

/* Low-mul Nussbaumer-style DFT implementation has 20 fewer MUL but at cost of 12 more ADD.
Totals: 72 ADD, 16 MUL	*/
//void ___RADIX_05_DFT_MUSS() {}
#define RADIX_07_DFT_NUSS(\
	__Ar0,__Ai0,__Ar1,__Ai1,__Ar2,__Ai2,__Ar3,__Ai3,__Ar4,__Ai4,__Ar5,__Ai5,__Ar6,__Ai6,\
	__tr0,__ti0,__tr1,__ti1,__tr2,__ti2,__tr3,__ti3,__tr4,__ti4,__tr5,__ti5,__tr6,__ti6,\
	__Br0,__Bi0,__Br1,__Bi1,__Br2,__Bi2,__Br3,__Bi3,__Br4,__Bi4,__Br5,__Bi5,__Br6,__Bi6,\
	__cx0,__sx0,__cx1,__sx1,__cx2,__sx2,__cx3,__sx3,\
	__rt,__it)\
{\
	__tr0 = __Ar0;						__ti0 = __Ai0;\
	__tr6 = __Ar1 - __Ar6;				__ti6 = __Ai1 - __Ai6;	/* x1 - x6 */\
	__tr1 = __Ar1 + __Ar6;				__ti1 = __Ai1 + __Ai6;	/* x1 + x6 */\
\
	__tr5 = __Ar2 - __Ar5;				__ti5 = __Ai2 - __Ai5;	/* x2 - x5 */\
	__tr2 = __Ar2 + __Ar5;				__ti2 = __Ai2 + __Ai5;	/* x2 + x5 */\
\
	__tr4 = __Ar3 - __Ar4;				__ti4 = __Ai3 - __Ai4;	/* x3 - x4 */\
	__tr3 = __Ar3 + __Ar4;				__ti3 = __Ai3 + __Ai4;	/* x3 + x4 */\
\
	__rt = __tr1+__tr2+__tr3;			__it = __ti1+__ti2+__ti3;\
	__Br0= __rt+__tr0;					__Bi0= __it+__ti0;\
	__tr0= __rt*__cx0+__tr0;			__ti0= __it*__cx0+__ti0;\
	__tr1= __tr1-__tr2;					__ti1= __ti1-__ti2;\
	__tr2= __tr3-__tr2;					__ti2= __ti3-__ti2;\
	__tr3=(__tr1+__tr2)*__cx3;			__ti3=(__ti1+__ti2)*__cx3;\
	__tr1= __tr1*__cx1;					__ti1= __ti1*__cx1;\
	__tr2= __tr2*__cx2;					__ti2= __ti2*__cx2;\
	__rt = __tr1-__tr3;					__it = __ti1-__ti3;\
	__tr2= __tr2-__tr3;					__ti2= __ti2-__ti3;\
\
	__tr1= __tr0-__rt-__tr2;			__ti1= __ti0-__it-__ti2;\
	__tr2= __tr0+__tr2;					__ti2= __ti0+__ti2;\
	__tr0= __tr0+__rt;					__ti0= __ti0+__it;\
\
	__tr3=(__tr6-__tr4+__tr5)*__sx0;	__ti3=(__ti6-__ti4+__ti5)*__sx0;\
	__tr6= __tr6-__tr5;					__ti6= __ti6-__ti5;\
	__tr5= __tr4+__tr5;					__ti5= __ti4+__ti5;\
	__tr4=(__tr5-__tr6)*__sx3;			__ti4=(__ti5-__ti6)*__sx3;\
	__tr6= __tr6*__sx1;					__ti6= __ti6*__sx1;\
	__tr5= __tr5*__sx2;					__ti5= __ti5*__sx2;\
	__tr6= __tr4+__tr6;					__ti6= __ti4+__ti6;\
	__tr5= __tr4-__tr5;					__ti5= __ti4-__ti5;\
\
	__tr4= __tr3-__tr6-__tr5;			__ti4= __ti3-__ti6-__ti5;\
	__tr5= __tr3+__tr5;					__ti5= __ti3+__ti5;\
	__tr3= __tr3+__tr6;					__ti3= __ti3+__ti6;\
\
	__Br1 =__tr0-__ti3;					__Bi1 =__ti0+__tr3;\
	__Br2 =__tr1-__ti4;					__Bi2 =__ti1+__tr4;\
	__Br3 =__tr2+__ti5;					__Bi3 =__ti2-__tr5;\
	__Br4 =__tr2-__ti5;					__Bi4 =__ti2+__tr5;\
	__Br5 =__tr1+__ti4;					__Bi5 =__ti1-__tr4;\
	__Br6 =__tr0+__ti3;					__Bi6 =__ti0-__tr3;\
}

/* Totals: 60 ADD, 36 MUL	*/
//void ___RADIX_07_DFT_PFETCH() {}
#define RADIX_07_DFT_PFETCH(\
	__Ar0,__Ai0,__Ar1,__Ai1,__Ar2,__Ai2,__Ar3,__Ai3,__Ar4,__Ai4,__Ar5,__Ai5,__Ar6,__Ai6,\
	__tr0,__ti0,__tr1,__ti1,__tr2,__ti2,__tr3,__ti3,__tr4,__ti4,__tr5,__ti5,__tr6,__ti6,\
	__Br0,__Bi0,__Br1,__Bi1,__Br2,__Bi2,__Br3,__Bi3,__Br4,__Bi4,__Br5,__Bi5,__Br6,__Bi6,\
	__cc1,__ss1,__cc2,__ss2,__cc3,__ss3,\
	__rt,__it,__re,__im,\
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
//void ___RADIX_08_DIF() {}
#define RADIX_08_DIF(\
	__Ar0,__Ai0,__Ar1,__Ai1,__Ar2,__Ai2,__Ar3,__Ai3,__Ar4,__Ai4,__Ar5,__Ai5,__Ar6,__Ai6,__Ar7,__Ai7,\
	__tr0,__ti0,__tr1,__ti1,__tr2,__ti2,__tr3,__ti3,__tr4,__ti4,__tr5,__ti5,__tr6,__ti6,__tr7,__ti7,\
	__Br0,__Bi0,__Br1,__Bi1,__Br2,__Bi2,__Br3,__Bi3,__Br4,__Bi4,__Br5,__Bi5,__Br6,__Bi6,__Br7,__Bi7,\
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
//void ___RADIX_08_DIF_OOP() {}
#define RADIX_08_DIF_OOP(\
	__Ar0,__Ai0,__Ar1,__Ai1,__Ar2,__Ai2,__Ar3,__Ai3,__Ar4,__Ai4,__Ar5,__Ai5,__Ar6,__Ai6,__Ar7,__Ai7,\
	__tr0,__ti0,__tr1,__ti1,__tr2,__ti2,__tr3,__ti3,__tr4,__ti4,__tr5,__ti5,__tr6,__ti6,__tr7,__ti7)\
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
//void ___RADIX_08_DIF_TWIDDLE_OOP() {}
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
//void ___RADIX_08_DIT_OOP() {}
#define RADIX_08_DIT_OOP(\
	__Ar0,__Ai0,__Ar1,__Ai1,__Ar2,__Ai2,__Ar3,__Ai3,__Ar4,__Ai4,__Ar5,__Ai5,__Ar6,__Ai6,__Ar7,__Ai7,\
	__tr0,__ti0,__tr1,__ti1,__tr2,__ti2,__tr3,__ti3,__tr4,__ti4,__tr5,__ti5,__tr6,__ti6,__tr7,__ti7)\
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
// The DIT version of this macro processes the twiddles in-order:
// For the scalar-data versino of this macro it doesn't matter whether the in/outputs are local
// scalars or array locations, but the __B-naming indicates that this macro is intended to serve
// as a prototype for a SIMD/ASM macro in which the __t-inputs are read from a block of contiguous
// local-allocated storage and the __B-outputs go into an array with arbitrary index stride:
//
#if 1	// Set = 1 to enable debug-prints in easier-to-follow 1st version
/* Original-version code: */\
//void ___RADIX_08_DIT_TWIDDLE_OOP() {}
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
	double __rt,__it;\
/*\
printf("DIT twiddles:\n");\
printf("1 = %20.10e %20.10e\n",__Wr1,__Wi1);\
printf("2 = %20.10e %20.10e\n",__Wr2,__Wi2);\
printf("3 = %20.10e %20.10e\n",__Wr3,__Wi3);\
printf("4 = %20.10e %20.10e\n",__Wr4,__Wi4);\
printf("5 = %20.10e %20.10e\n",__Wr5,__Wi5);\
printf("6 = %20.10e %20.10e\n",__Wr6,__Wi6);\
printf("7 = %20.10e %20.10e\n",__Wr7,__Wi7);\
*/\
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
/*\
printf("DIT CMUL Outputs:\n");\
printf("0 = %20.10e %20.10e\n",__tr0,__ti0);\
printf("1 = %20.10e %20.10e\n",__tr1,__ti1);\
printf("2 = %20.10e %20.10e\n",__tr2,__ti2);\
printf("3 = %20.10e %20.10e\n",__tr3,__ti3);\
printf("4 = %20.10e %20.10e\n",__tr4,__ti4);\
printf("5 = %20.10e %20.10e\n",__tr5,__ti5);\
printf("6 = %20.10e %20.10e\n",__tr6,__ti6);\
printf("7 = %20.10e %20.10e\n",__tr7,__ti7);\
*/\
/*** Combine to get the 2 length-4 transforms... ***/\
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
/*\
printf("DIT Step 2 Outputs:\n");\
printf("0 = %20.10e %20.10e\n",__tr0,__ti0);\
printf("1 = %20.10e %20.10e\n",__tr1,__ti1);\
printf("2 = %20.10e %20.10e\n",__tr2,__ti2);\
printf("3 = %20.10e %20.10e\n",__tr3,__ti3);\
printf("4 = %20.10e %20.10e\n",__tr4,__ti4);\
printf("5 = %20.10e %20.10e\n",__tr5,__ti5);\
printf("6 = %20.10e %20.10e\n",__tr6,__ti6);\
printf("7 = %20.10e %20.10e\n",__tr7,__ti7);\
*/\
/*** Now combine the two half-transforms: ***/\
	__rt  = __tr4;			__it  = __ti4;\
	__Br0 = __tr0 + __rt;	__Bi0 = __ti0 + __it;\
	__Br1 = __tr0 - __rt;	__Bi1 = __ti0 - __it;\
\
	__rt  = __tr6;			__it  = __ti6;\
	__Br2 = __tr2 + __it;	__Bi2 = __ti2 - __rt;\
	__Br3 = __tr2 - __it;	__Bi3 = __ti2 + __rt;\
\
	__rt  = __tr5 + __ti5;	__it  = __tr5 - __ti5;/* rt,it swapped vs DIF */\
	__rt *= ISRT2;			__it *= ISRT2;\
	__Br4 = __tr1 + __rt;	__Bi4 = __ti1 - __it;\
	__Br5 = __tr1 - __rt;	__Bi5 = __ti1 + __it;\
\
	__rt  = __tr7 - __ti7;	__it  = __tr7 + __ti7;/* rt,it = DIF(-it,+rt) */\
	__rt *= ISRT2;			__it *= ISRT2;\
	__Br6 = __tr3 - __rt;	__Bi6 = __ti3 - __it;\
	__Br7 = __tr3 + __rt;	__Bi7 = __ti3 + __it;\
}
#else
/* Rewrite with x86 ISA and SSE2/16-register target: Assume 16 registers _r0-f available, all arithmetic is on these: */
/* In 16-register mode Block pairs 0/1,2/3 and 4/5,6/7 can be done at same time; 8-register mode just replace _r8-f with _r0-7. */
//void ___RADIX_08_DIT_TWIDDLE_OOP() {}
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
#endif	// 0|1 toggle


/* Totals: 52 ADD, 4 MUL	*/
//void ___RADIX_08_DIF_PFETCH() {}
#define RADIX_08_DIF_PFETCH(\
	__Ar0,__Ai0,__Ar1,__Ai1,__Ar2,__Ai2,__Ar3,__Ai3,__Ar4,__Ai4,__Ar5,__Ai5,__Ar6,__Ai6,__Ar7,__Ai7,\
	__tr0,__ti0,__tr1,__ti1,__tr2,__ti2,__tr3,__ti3,__tr4,__ti4,__tr5,__ti5,__tr6,__ti6,__tr7,__ti7,\
	__Br0,__Bi0,__Br1,__Bi1,__Br2,__Bi2,__Br3,__Bi3,__Br4,__Bi4,__Br5,__Bi5,__Br6,__Bi6,__Br7,__Bi7,\
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
//void ___RADIX_08_DIT() {}
#define RADIX_08_DIT(\
	__Ar0,__Ai0,__Ar1,__Ai1,__Ar2,__Ai2,__Ar3,__Ai3,__Ar4,__Ai4,__Ar5,__Ai5,__Ar6,__Ai6,__Ar7,__Ai7,\
	__tr0,__ti0,__tr1,__ti1,__tr2,__ti2,__tr3,__ti3,__tr4,__ti4,__tr5,__ti5,__tr6,__ti6,__tr7,__ti7,\
	__Br0,__Bi0,__Br1,__Bi1,__Br2,__Bi2,__Br3,__Bi3,__Br4,__Bi4,__Br5,__Bi5,__Br6,__Bi6,__Br7,__Bi7,\
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

// Radix-9 DIF: A/Bs are in/outputs and t's are temporaries (all doubles). Totals: [80 ADD, 40 MUL]
/*
Apr 2016: On revisiting the radix-9 DFTs (cf. FMA-oriented macros below), my test_fft_radix code
says the below DIF outputs need to be permuted as [0,3,6,1,4,7,2,5,8], i.e. should occur in-order:
[0,1,2,3,4,5,6,7,8] ... above scrambling permits in-place-ness w.r.to the non-array temporaries t00-t0h.
In any event said higher radices are based on the above ordering, just something for future readers of
this code to be aware of. If you need a radix-9 DIF macro with the proper output ordering use the
RADIX_09_DIF_FMA-named analog of this macro below.
*/
//void ___RADIX_09_DIF() {}
#define RADIX_09_DIF(\
	__A0r,__A0i,__A1r,__A1i,__A2r,__A2i,__A3r,__A3i,__A4r,__A4i,__A5r,__A5i,__A6r,__A6i,__A7r,__A7i,__A8r,__A8i,\
	__t00,__t01,__t02,__t03,__t04,__t05,__t06,__t07,__t08,__t09,__t0a,__t0b,__t0c,__t0d,__t0e,__t0f,__t0g,__t0h,\
	__rt,__it,__tt)\
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

// Radix-9 DIT: t's are temporaries and As are outputs (all doubles). Totals: [80 ADD, 40 MUL]
/*
Apr 2016: On revisiting the radix-9 DFTs (cf. FMA-oriented macros below), my test_fft_radix code
says the above DIT outputs need to be permuted as [0,4,8,3,7,2,6,1,5], i.e. should occur in order
[0,3,6,1,4,7,2,5,8] ... not sure now why I use the above scrambling but suspect it has to do with
easing the unscrambling/in-place-ness of higher-radices based on radix-9 such as 36.144,288 ...
in the DIF version of radix-9 a similar output scrambling is justified by the outputs being written
what are assumed to be local scalars, but that is not the case for DIT.
In any event said higher radices are based on the above ordering, just something for futre readers of
this code to be aware of. If you need a radix-9 DIT macro with the proper output ordering use the
RADIX_09_DIT_FMA-named analog of this macro below.
*/
//void ___RADIX_09_DIT() {}
#define RADIX_09_DIT(\
	__t00,__t01,__t02,__t03,__t04,__t05,__t06,__t07,__t08,__t09,__t0a,__t0b,__t0c,__t0d,__t0e,__t0f,__t0g,__t0h,\
	__A0r,__A0i,__A1r,__A1i,__A2r,__A2i,__A3r,__A3i,__A4r,__A4i,__A5r,__A5i,__A6r,__A6i,__A7r,__A7i,__A8r,__A8i,\
	__rt,__it,__tt)\
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
/* If A and t-terms map to same memlocs (i.e. caller doing an in-place radix-9 DIT), above
radix-3 pass 2 output-writes just clobbered t00,1,6,7,c,d (via A0,3,6)...  */\
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


/* FMA-optimized Radix-9 DIF -- Totals: [84 ADD, 44 MUL], thus slightly worse than above
two-[radix-3-triplet]-passes version for non-FMA implementation. With FMA, the cost drops
to [44 ADD, 44 MUL/FMA], in FMA-heaviest version (but still one with no trivial FMA). The
latter is just a tad better than the [46 ADD, 44 MUL/FMA] cost of my AVX2 implementation
of the above two-[radix-3-triplet]-passes version.

Basic approach here same as for add-prime radices, plus some added opts resulting from the
special form of the E3-multiplied terms in rows 3,6, cols 3,6 of our matrix-multiply DFT.
Base-root is E := exp(I*2*Pi/9), here are needed powers and symmetries:

	E   =  0.7660444431189780352023926506 + 0.6427876096865393263226434099*I
	E^2 =  0.1736481776669303488517166268 + 0.9848077530122080593667430246*I
	E^3 = -0.5000000000000000000000000000 + 0.8660254037844386467637231707*I = [-1 + I*sqrt(3)]/2
	E^4 = -0.9396926207859083840541092773 + 0.3420201433256687330440996146*I
	E^5,6,78 = ~E4,3,2,1, and higher powers map to these simply via E^n = E^(n mod 9).

Thus our DFT matrix is

		1	1	1	1	1	|	1	1	1	1
		--------------------|------------------
		1	E1	E2	E3	E4	|	E4~	E3~	E2~	E1~
		1	E2	E4	E3~	E1~	|	E1	E3	E4~	E2~
		1	E3	E3~	1	E3	|	E3~	1	E3	E3~
		1	E4	E1~	E3	E2~	|	E2	E3~	E1	E4~
		--------------------|------------------
		1	E4~	E1	E3~	E2	|	E2~	E3	E1~	E4
		1	E3~	E3	1	E3~	|	E3	1	E3~	E3
		1	E2~	E4~	E3	E1	|	E1~	E3~	E4	E2
		1	E1~	E2~	E3~	E4~	|	E4	E3	E2	E1 ;

the upper left 3 x 3 E-submatrix gives the form of our mini-convolutions, with ~ terms having negative
sine multipliers, i.e. Ej~ ==> -sin(2*Pi*j/9) .
*/
//void ___RADIX_09_DIF_FMA() {}
#define RADIX_09_DIF_FMA(\
	__Ar0,__Ai0,__Ar1,__Ai1,__Ar2,__Ai2,__Ar3,__Ai3,__Ar4,__Ai4,__Ar5,__Ai5,__Ar6,__Ai6,__Ar7,__Ai7,__Ar8,__Ai8,\
	__Br0,__Bi0,__Br1,__Bi1,__Br2,__Bi2,__Br3,__Bi3,__Br4,__Bi4,__Br5,__Bi5,__Br6,__Bi6,__Br7,__Bi7,__Br8,__Bi8,\
	__cc1,__ss1, __cc2,__ss2, __cc3,__ss3, __cc4,__ss4)\
{\
double __tr0,__tr1,__tr2,__tr3,__tr4,__tr5,__tr6,__tr7, __ti0,__ti1,__ti2,__ti3,__ti4,__ti5,__ti6,__ti7, __rt,__it;\
	__rt  = __Ar8;								__it  = __Ai8;	/* Copies needed for in-place capability */\
	__Br8 = __Ar1 - __rt;						__Bi8 = __Ai1 - __it;	/* x1 - x8 */\
	__Br1 = __Ar1 + __rt;						__Bi1 = __Ai1 + __it;	/* x1 + x8 */\
	__rt  = __Ar7;								__it  = __Ai7;	/* Copies needed for in-place capability */\
	__Br7 = __Ar2 - __rt;						__Bi7 = __Ai2 - __it;	/* x2 - x7 */\
	__Br2 = __Ar2 + __rt;						__Bi2 = __Ai2 + __it;	/* x2 + x7 */\
	__rt  = __Ar6;								__it  = __Ai6;	/* Copies needed for in-place capability */\
	__Br6 = __Ar3 - __rt;						__Bi6 = __Ai3 - __it;	/* x3 - x6 */\
	__Br3 = __Ar3 + __rt;						__Bi3 = __Ai3 + __it;	/* x3 + x6 */\
	__rt  = __Ar5;								__it  = __Ai5;	/* Copies needed for in-place capability */\
	__Br5 = __Ar4 - __rt;						__Bi5 = __Ai4 - __it;	/* x4 - x5 */\
	__Br4 = __Ar4 + __rt;						__Bi4 = __Ai4 + __it;	/* x4 + x5 *//* 16 add */\
	/* Total cost of 'mini-convo' mid-sequence: [52 add, 44 mul] for non-fma; [12 add, 44 fma] with fma.
	Must compute A0 + cc3*B3 before B0 here since A0,B0 may coincide: */\
	__tr3 = __Ar0 + __cc3*__Br3;				__ti3 = __Ai0 + __cc3*__Bi3;	/* 2 fma */\
	__Br0 = __Ar0 + __Br3;						__Bi0 = __Ai0 + __Bi3;			/* 2 add */\
	__rt  = __Br1 + __Br2 + __Br4;				__it  = __Bi1 + __Bi2 + __Bi4;	/* 4 add */\
	__tr2 = __Br0 + __cc3*__rt;					__ti2 = __Bi0 + __cc3*__it;		/* 2 fma */\
	__Br0 += __rt;								__Bi0 += __it;	/* DC output */	/* 2 add */\
	__tr0 = __tr3 + __cc1*__Br1 + __cc2*__Br2 + __cc4*__Br4;	/* 6 fma (Combining counts of Re,Im-parts) */\
	__tr1 = __tr3 + __cc2*__Br1 + __cc4*__Br2 + __cc1*__Br4;	/* 6 fma */\
	__tr3 +=        __cc4*__Br1 + __cc1*__Br2 + __cc2*__Br4;	/* 6 fma */\
												__ti0 = __ti3 + __cc1*__Bi1 + __cc2*__Bi2 + __cc4*__Bi4;\
												__ti1 = __ti3 + __cc2*__Bi1 + __cc4*__Bi2 + __cc1*__Bi4;\
												__ti3 +=        __cc4*__Bi1 + __cc1*__Bi2 + __cc2*__Bi4;\
	__rt  = __ss3*__Br6;						__it  = __ss3*__Bi6;/* 2 mul */\
	__tr4 =  __rt + __ss1*__Br8 + __ss2*__Br7 + __ss4*__Br5;	/* 6 fma (Combining counts of Re,Im-parts) */\
	__tr5 = -__rt + __ss2*__Br8 + __ss4*__Br7 - __ss1*__Br5;	/* 6 fma */\
	__tr6 =        __ss3*(__Br8 -       __Br7 +       __Br5);	/* 2 mul, 4 add */\
	__tr7 =  __rt + __ss4*__Br8 - __ss1*__Br7 - __ss2*__Br5;	/* 6 fma */\
												__ti4 =  __it + __ss1*__Bi8 + __ss2*__Bi7 + __ss4*__Bi5;\
												__ti5 = -__it + __ss2*__Bi8 + __ss4*__Bi7 - __ss1*__Bi5;\
												__ti6 =        __ss3*(__Bi8 -       __Bi7 +       __Bi5);\
												__ti7 =  __it + __ss4*__Bi8 - __ss1*__Bi7 - __ss2*__Bi5;\
	/* Output permutation = [0,3,6,1,4,7,2,5,8]: */\
	__Br3 = __tr0 - __ti4;						__Bi3 = __ti0 + __tr4;\
	__Br6 = __tr1 - __ti5;						__Bi6 = __ti1 + __tr5;\
	__Br1 = __tr2 - __ti6;						__Bi1 = __ti2 + __tr6;\
	__Br4 = __tr3 - __ti7;						__Bi4 = __ti3 + __tr7;\
	__Br7 = __tr3 + __ti7;						__Bi7 = __ti3 - __tr7;\
	__Br2 = __tr2 + __ti6;						__Bi2 = __ti2 - __tr6;\
	__Br5 = __tr1 + __ti5;						__Bi5 = __ti1 - __tr5;\
	__Br8 = __tr0 + __ti4;						__Bi8 = __ti0 - __tr4;	/* 16 add */\
}

//void ___RADIX_09_DIT_FMA() {}
#define RADIX_09_DIT_FMA(\
	__Ar0,__Ai0,__Ar1,__Ai1,__Ar2,__Ai2,__Ar3,__Ai3,__Ar4,__Ai4,__Ar5,__Ai5,__Ar6,__Ai6,__Ar7,__Ai7,__Ar8,__Ai8,\
	__Br0,__Bi0,__Br1,__Bi1,__Br2,__Bi2,__Br3,__Bi3,__Br4,__Bi4,__Br5,__Bi5,__Br6,__Bi6,__Br7,__Bi7,__Br8,__Bi8,\
	__cc1,__ss1, __cc2,__ss2, __cc3,__ss3, __cc4,__ss4)\
{\
double __tr0,__tr1,__tr2,__tr3,__tr4,__tr5,__tr6,__tr7, __ti0,__ti1,__ti2,__ti3,__ti4,__ti5,__ti6,__ti7, __rt,__it;\
/* Input ordering [0,3,6,1,4,7,2,5,8], i.e. unpermuted index pairs 1/8,2/7,3/6,4/5 map to 3/8,6/5,1/2,4/7: */\
	__rt  = __Ar8;								__it  = __Ai8;	/* Copies needed for in-place capability */\
	__Br8 = __Ar3 - __rt;						__Bi8 = __Ai3 - __it;	/* x1 - x8 */\
	__Br3 = __Ar3 + __rt;						__Bi3 = __Ai3 + __it;	/* x1 + x8 */\
	__rt  = __Ar5;								__it  = __Ai5;	/* Copies needed for in-place capability */\
	__Br5 = __Ar6 - __rt;						__Bi5 = __Ai6 - __it;	/* x2 - x7 */\
	__Br6 = __Ar6 + __rt;						__Bi6 = __Ai6 + __it;	/* x2 + x7 */\
	__rt  = __Ar2;								__it  = __Ai2;	/* Copies needed for in-place capability */\
	__Br2 = __Ar1 - __rt;						__Bi2 = __Ai1 - __it;	/* x3 - x6 */\
	__Br1 = __Ar1 + __rt;						__Bi1 = __Ai1 + __it;	/* x3 + x6 */\
	__rt  = __Ar7;								__it  = __Ai7;	/* Copies needed for in-place capability */\
	__Br7 = __Ar4 - __rt;						__Bi7 = __Ai4 - __it;	/* x4 - x5 */\
	__Br4 = __Ar4 + __rt;						__Bi4 = __Ai4 + __it;	/* x4 + x5 */\
/* B-terms here subject to same permutation as inputs, [0,3,6,1,4,7,2,5,8]: */\
	/* Must compute A0 + cc3*B3 before B0 here since A0,B0 may coincide: */\
	__tr3 = __Ar0 + __cc3*__Br1;				__ti3 = __Ai0 + __cc3*__Bi1;	\
	__Br0 = __Ar0 + __Br1;						__Bi0 = __Ai0 + __Bi1;			\
	__rt  = __Br3 + __Br6 + __Br4;				__it  = __Bi3 + __Bi6 + __Bi4;	\
	__tr2 = __Br0 + __cc3*__rt;					__ti2 = __Bi0 + __cc3*__it;		\
	__Br0 += __rt;								__Bi0 += __it;	/* DB output */	\
	__tr0 = __tr3 + __cc1*__Br3 + __cc2*__Br6 + __cc4*__Br4;\
	__tr1 = __tr3 + __cc2*__Br3 + __cc4*__Br6 + __cc1*__Br4;\
	__tr3 +=        __cc4*__Br3 + __cc1*__Br6 + __cc2*__Br4;\
												__ti0 = __ti3 + __cc1*__Bi3 + __cc2*__Bi6 + __cc4*__Bi4;\
												__ti1 = __ti3 + __cc2*__Bi3 + __cc4*__Bi6 + __cc1*__Bi4;\
												__ti3 +=        __cc4*__Bi3 + __cc1*__Bi6 + __cc2*__Bi4;\
	__rt  = __ss3*__Br2;						__it  = __ss3*__Bi2;\
	__tr4 =  __rt + __ss1*__Br8 + __ss2*__Br5 + __ss4*__Br7;\
	__tr5 = -__rt + __ss2*__Br8 + __ss4*__Br5 - __ss1*__Br7;\
	__tr6 =        __ss3*(__Br8 -       __Br5 +       __Br7);\
	__tr7 =  __rt + __ss4*__Br8 - __ss1*__Br5 - __ss2*__Br7;\
												__ti4 =  __it + __ss1*__Bi8 + __ss2*__Bi5 + __ss4*__Bi7;\
												__ti5 = -__it + __ss2*__Bi8 + __ss4*__Bi5 - __ss1*__Bi7;\
												__ti6 =        __ss3*(__Bi8 -       __Bi5 +       __Bi7);\
												__ti7 =  __it + __ss4*__Bi8 - __ss1*__Bi5 - __ss2*__Bi7;\
	/* No output permutation in DIT case - note if wanted to share roots with DIF (i.e. use non-conjugated
	9th roots of unity here), would need to flip +- signs here: */\
	__Br1 = __tr0 + __ti4;						__Bi1 = __ti0 - __tr4;\
	__Br2 = __tr1 + __ti5;						__Bi2 = __ti1 - __tr5;\
	__Br3 = __tr2 + __ti6;						__Bi3 = __ti2 - __tr6;\
	__Br4 = __tr3 + __ti7;						__Bi4 = __ti3 - __tr7;\
	__Br5 = __tr3 - __ti7;						__Bi5 = __ti3 + __tr7;\
	__Br6 = __tr2 - __ti6;						__Bi6 = __ti2 + __tr6;\
	__Br7 = __tr1 - __ti5;						__Bi7 = __ti1 + __tr5;\
	__Br8 = __tr0 - __ti4;						__Bi8 = __ti0 + __tr4;\
}

/****** RADIX = 11: ******/

/*...Simple-algo Radix-11 DFT.  Totals: 140 ADD, 100 MUL	*/
//void ___RADIX_11_DFT_BASIC() {}
#define RADIX_11_DFT_BASIC(\
	__A0r,__A0i,__A1r,__A1i,__A2r,__A2i,__A3r,__A3i,__A4r,__A4i,__A5r,__A5i,__A6r,__A6i,__A7r,__A7i,__A8r,__A8i,__A9r,__A9i,__Aar,__Aai,\
	__B0r,__B0i,__B1r,__B1i,__B2r,__B2i,__B3r,__B3i,__B4r,__B4i,__B5r,__B5i,__B6r,__B6i,__B7r,__B7i,__B8r,__B8i,__B9r,__B9i,__Bar,__Bai,\
	_cc1,_cc2,_cc3,_cc4,_cc5,_ss1,_ss2,_ss3,_ss4,_ss5)\
{\
	double _t0r,_t0i,_t1r,_t1i,_t2r,_t2i,_t3r,_t3i,_t4r,_t4i,_t5r,_t5i,_t6r,_t6i,_t7r,_t7i,_t8r,_t8i,_t9r,_t9i,_tar,_tai;\
	double _cr1,_cr2,_cr3,_cr4,_cr5,_ci1,_ci2,_ci3,_ci4,_ci5,_sr1,_sr2,_sr3,_sr4,_sr5,_si1,_si2,_si3,_si4,_si5;\
\
	_t0r = __A0r       ;					_t0i = __A0i        ;	/* x0		*/\
	_t1r = __A1r +__Aar;					_t1i = __A1i + __Aai;	/* x1 + x10	*/\
	_t2r = __A2r +__A9r;					_t2i = __A2i + __A9i;	/* x2 + x9	*/\
	_t3r = __A3r +__A8r;					_t3i = __A3i + __A8i;	/* x3 + x8	*/\
	_t4r = __A4r +__A7r;					_t4i = __A4i + __A7i;	/* x4 + x7	*/\
	_t5r = __A5r +__A6r;					_t5i = __A5i + __A6i;	/* x5 + x6	*/\
	_t6r = __A5r -__A6r;					_t6i = __A5i - __A6i;	/* x5 - x6	*/\
	_t7r = __A4r -__A7r;					_t7i = __A4i - __A7i;	/* x4 - x7	*/\
	_t8r = __A3r -__A8r;					_t8i = __A3i - __A8i;	/* x3 - x8	*/\
	_t9r = __A2r -__A9r;					_t9i = __A2i - __A9i;	/* x2 - x9	*/\
	_tar = __A1r -__Aar;					_tai = __A1i - __Aai;	/* x1 - x10	*/\
/*\
	_cr1 = _t0r+_cc1*_t1r+_cc2*_t2r+_cc3*_t3r+_cc4*_t4r+_cc5*_t5r;	_ci1 = _t0i+_cc1*_t1i+_cc2*_t2i+_cc3*_t3i+_cc4*_t4i+_cc5*_t5i;	// C1 //\
	_cr2 = _t0r+_cc2*_t1r+_cc4*_t2r+_cc5*_t3r+_cc3*_t4r+_cc1*_t5r;	_ci2 = _t0i+_cc2*_t1i+_cc4*_t2i+_cc5*_t3i+_cc3*_t4i+_cc1*_t5i;	// C2 //\
	_cr3 = _t0r+_cc3*_t1r+_cc5*_t2r+_cc2*_t3r+_cc1*_t4r+_cc4*_t5r;	_ci3 = _t0i+_cc3*_t1i+_cc5*_t2i+_cc2*_t3i+_cc1*_t4i+_cc4*_t5i;	// C3 //\
	_cr4 = _t0r+_cc4*_t1r+_cc3*_t2r+_cc1*_t3r+_cc5*_t4r+_cc2*_t5r;	_ci4 = _t0i+_cc4*_t1i+_cc3*_t2i+_cc1*_t3i+_cc5*_t4i+_cc2*_t5i;	// C4 //\
	_cr5 = _t0r+_cc5*_t1r+_cc1*_t2r+_cc4*_t3r+_cc2*_t4r+_cc3*_t5r;	_ci5 = _t0i+_cc5*_t1i+_cc1*_t2i+_cc4*_t3i+_cc2*_t4i+_cc3*_t5i;	// C5 //\
\
	_sr1 =      _ss1*_tar+_ss2*_t9r+_ss3*_t8r+_ss4*_t7r+_ss5*_t6r;	_si1 =      _ss1*_tai+_ss2*_t9i+_ss3*_t8i+_ss4*_t7i+_ss5*_t6i;	// S1 //\
	_sr2 =      _ss2*_tar+_ss4*_t9r-_ss5*_t8r-_ss3*_t7r-_ss1*_t6r;	_si2 =      _ss2*_tai+_ss4*_t9i-_ss5*_t8i-_ss3*_t7i-_ss1*_t6i;	// S2 //\
	_sr3 =      _ss3*_tar-_ss5*_t9r-_ss2*_t8r+_ss1*_t7r+_ss4*_t6r;	_si3 =      _ss3*_tai-_ss5*_t9i-_ss2*_t8i+_ss1*_t7i+_ss4*_t6i;	// S3 //\
	_sr4 =      _ss4*_tar-_ss3*_t9r+_ss1*_t8r+_ss5*_t7r-_ss2*_t6r;	_si4 =      _ss4*_tai-_ss3*_t9i+_ss1*_t8i+_ss5*_t7i-_ss2*_t6i;	// S4 //\
	_sr5 =      _ss5*_tar-_ss1*_t9r+_ss4*_t8r-_ss2*_t7r+_ss3*_t6r;	_si5 =      _ss5*_tai-_ss1*_t9i+_ss4*_t8i-_ss2*_t7i+_ss3*_t6i;	// S5 //\
\
	__B0r = _t0r+_t1r+_t2r+_t3r+_t4r+_t5r;	__B0i = _t0i+_t1i+_t2i+_t3i+_t4i+_t5i;	// X0 //\
*/\
	/* Try a version which aims to make it esier for the compiler to do a good job: */\
	__B0r = _t0r;	__B0i = _t0i;	/* X0 accumulators */\
	_cr1 = _t0r+_cc1*_t1r;	_ci1 = _t0i+_cc1*_t1i;\
	_cr2 = _t0r+_cc2*_t1r;	_ci2 = _t0i+_cc2*_t1i;\
	_cr3 = _t0r+_cc3*_t1r;	_ci3 = _t0i+_cc3*_t1i;\
	_cr4 = _t0r+_cc4*_t1r;	_ci4 = _t0i+_cc4*_t1i;\
	_cr5 = _t0r+_cc5*_t1r;	_ci5 = _t0i+_cc5*_t1i;\
	__B0r += _t1r;	__B0i += _t1i;\
\
	_cr1 += _cc2*_t2r;	_ci1 += _cc2*_t2i;\
	_cr2 += _cc4*_t2r;	_ci2 += _cc4*_t2i;\
	_cr3 += _cc5*_t2r;	_ci3 += _cc5*_t2i;\
	_cr4 += _cc3*_t2r;	_ci4 += _cc3*_t2i;\
	_cr5 += _cc1*_t2r;	_ci5 += _cc1*_t2i;\
	__B0r += _t2r;	__B0i += _t2i;\
\
	_cr1 += _cc3*_t3r;	_ci1 += _cc3*_t3i;\
	_cr2 += _cc5*_t3r;	_ci2 += _cc5*_t3i;\
	_cr3 += _cc2*_t3r;	_ci3 += _cc2*_t3i;\
	_cr4 += _cc1*_t3r;	_ci4 += _cc1*_t3i;\
	_cr5 += _cc4*_t3r;	_ci5 += _cc4*_t3i;\
	__B0r += _t3r;	__B0i += _t3i;\
\
	_cr1 += _cc4*_t4r;	_ci1 += _cc4*_t4i;\
	_cr2 += _cc3*_t4r;	_ci2 += _cc3*_t4i;\
	_cr3 += _cc1*_t4r;	_ci3 += _cc1*_t4i;\
	_cr4 += _cc5*_t4r;	_ci4 += _cc5*_t4i;\
	_cr5 += _cc2*_t4r;	_ci5 += _cc2*_t4i;\
	__B0r += _t4r;	__B0i += _t4i;\
\
	_cr1 += _cc5*_t5r;	_ci1 += _cc5*_t5i;\
	_cr2 += _cc1*_t5r;	_ci2 += _cc1*_t5i;\
	_cr3 += _cc4*_t5r;	_ci3 += _cc4*_t5i;\
	_cr4 += _cc2*_t5r;	_ci4 += _cc2*_t5i;\
	_cr5 += _cc3*_t5r;	_ci5 += _cc3*_t5i;\
	__B0r += _t5r;	__B0i += _t5i;\
\
	_sr1 =      _ss1*_tar;	_si1 =      _ss1*_tai;\
	_sr2 =      _ss2*_tar;	_si2 =      _ss2*_tai;\
	_sr3 =      _ss3*_tar;	_si3 =      _ss3*_tai;\
	_sr4 =      _ss4*_tar;	_si4 =      _ss4*_tai;\
	_sr5 =      _ss5*_tar;	_si5 =      _ss5*_tai;\
\
	_sr1 += _ss2*_t9r;	_si1 += _ss2*_t9i;\
	_sr2 += _ss4*_t9r;	_si2 += _ss4*_t9i;\
	_sr3 -= _ss5*_t9r;	_si3 -= _ss5*_t9i;\
	_sr4 -= _ss3*_t9r;	_si4 -= _ss3*_t9i;\
	_sr5 -= _ss1*_t9r;	_si5 -= _ss1*_t9i;\
\
	_sr1 += _ss3*_t8r;	_si1 += _ss3*_t8i;\
	_sr2 -= _ss5*_t8r;	_si2 -= _ss5*_t8i;\
	_sr3 -= _ss2*_t8r;	_si3 -= _ss2*_t8i;\
	_sr4 += _ss1*_t8r;	_si4 += _ss1*_t8i;\
	_sr5 += _ss4*_t8r;	_si5 += _ss4*_t8i;\
\
	_sr1 += _ss4*_t7r;	_si1 += _ss4*_t7i;\
	_sr2 -= _ss3*_t7r;	_si2 -= _ss3*_t7i;\
	_sr3 += _ss1*_t7r;	_si3 += _ss1*_t7i;\
	_sr4 += _ss5*_t7r;	_si4 += _ss5*_t7i;\
	_sr5 -= _ss2*_t7r;	_si5 -= _ss2*_t7i;\
\
	_sr1 += _ss5*_t6r;	_si1 += _ss5*_t6i;\
	_sr2 -= _ss1*_t6r;	_si2 -= _ss1*_t6i;\
	_sr3 += _ss4*_t6r;	_si3 += _ss4*_t6i;\
	_sr4 -= _ss2*_t6r;	_si4 -= _ss2*_t6i;\
	_sr5 += _ss3*_t6r;	_si5 += _ss3*_t6i;\
/**/\
	__B1r = _cr1 - _si1;					__B1i = _ci1 + _sr1;	/* X1 = C1 + I*S1	*/\
	__B2r = _cr2 - _si2;					__B2i = _ci2 + _sr2;	/* X2 = C2 + I*S2	*/\
	__B3r = _cr3 - _si3;					__B3i = _ci3 + _sr3;	/* X3 = C3 + I*S3	*/\
	__B4r = _cr4 - _si4;					__B4i = _ci4 + _sr4;	/* X4 = C4 + I*S4	*/\
	__B5r = _cr5 - _si5;					__B5i = _ci5 + _sr5;	/* X5 = C5 + I*S5	*/\
	__B6r = _cr5 + _si5;					__B6i = _ci5 - _sr5;	/* X6 =	C5 - I*S5	*/\
	__B7r = _cr4 + _si4;					__B7i = _ci4 - _sr4;	/* X7 =	C4 - I*S4	*/\
	__B8r = _cr3 + _si3;					__B8i = _ci3 - _sr3;	/* X8 =	C3 - I*S3	*/\
	__B9r = _cr2 + _si2;					__B9i = _ci2 - _sr2;	/* X9 =	C2 - I*S2	*/\
	__Bar = _cr1 + _si1;					__Bai = _ci1 - _sr1;	/* X10=	C1 - I*S1	*/\
}

/* FMA-oriented version - above simple-structure better for FMA than Nussbaumer:
Totals: 50 ADD, 10 MUL, 90 FMA, i.e. trade 90 ADD + 90 MUL for 90 FMA, very nice.
*/
//void ___RADIX_11_DFT_FMA() {}
#define RADIX_11_DFT_FMA(\
	__A0r,__A0i,__A1r,__A1i,__A2r,__A2i,__A3r,__A3i,__A4r,__A4i,__A5r,__A5i,__A6r,__A6i,__A7r,__A7i,__A8r,__A8i,__A9r,__A9i,__Aar,__Aai,\
	__B0r,__B0i,__B1r,__B1i,__B2r,__B2i,__B3r,__B3i,__B4r,__B4i,__B5r,__B5i,__B6r,__B6i,__B7r,__B7i,__B8r,__B8i,__B9r,__B9i,__Bar,__Bai,\
	_cc1,_cc2,_cc3,_cc4,_cc5,_ss1,_ss2,_ss3,_ss4,_ss5)\
{\
	double _t0r,_t0i,_t1r,_t1i,_t2r,_t2i,_t3r,_t3i,_t4r,_t4i,_t5r,_t5i,_t6r,_t6i,_t7r,_t7i,_t8r,_t8i,_t9r,_t9i,_tar,_tai;\
	double _cr1,_cr2,_cr3,_cr4,_cr5,_ci1,_ci2,_ci3,_ci4,_ci5,_sr1,_sr2,_sr3,_sr4,_sr5,_si1,_si2,_si3,_si4,_si5;\
\
	_t0r = __A0r       ;					_t0i = __A0i        ;	/* x0		*/\
	_t1r = __A1r +__Aar;					_t1i = __A1i + __Aai;	/* x1 + x10	*/\
	_t2r = __A2r +__A9r;					_t2i = __A2i + __A9i;	/* x2 + x9	*/\
	_t3r = __A3r +__A8r;					_t3i = __A3i + __A8i;	/* x3 + x8	*/\
	_t4r = __A4r +__A7r;					_t4i = __A4i + __A7i;	/* x4 + x7	*/\
	_t5r = __A5r +__A6r;					_t5i = __A5i + __A6i;	/* x5 + x6	*/\
	_tar = __A1r -__Aar;					_tai = __A1i - __Aai;	/* x1 - x10	*/\
	_t9r = __A2r -__A9r;					_t9i = __A2i - __A9i;	/* x2 - x9	*/\
	_t8r = __A3r -__A8r;					_t8i = __A3i - __A8i;	/* x3 - x8	*/\
	_t7r = __A4r -__A7r;					_t7i = __A4i - __A7i;	/* x4 - x7	*/\
	_t6r = __A5r -__A6r;					_t6i = __A5i - __A6i;	/* x5 - x6	*/\
/*\
	_cr1 = _t0r+_cc1*_t1r+_cc2*_t2r+_cc3*_t3r+_cc4*_t4r+_cc5*_t5r;	_ci1 = _t0i+_cc1*_t1i+_cc2*_t2i+_cc3*_t3i+_cc4*_t4i+_cc5*_t5i;	// C1
	_cr2 = _t0r+_cc2*_t1r+_cc4*_t2r+_cc5*_t3r+_cc3*_t4r+_cc1*_t5r;	_ci2 = _t0i+_cc2*_t1i+_cc4*_t2i+_cc5*_t3i+_cc3*_t4i+_cc1*_t5i;	// C2
	_cr3 = _t0r+_cc3*_t1r+_cc5*_t2r+_cc2*_t3r+_cc1*_t4r+_cc4*_t5r;	_ci3 = _t0i+_cc3*_t1i+_cc5*_t2i+_cc2*_t3i+_cc1*_t4i+_cc4*_t5i;	// C3
	_cr4 = _t0r+_cc4*_t1r+_cc3*_t2r+_cc1*_t3r+_cc5*_t4r+_cc2*_t5r;	_ci4 = _t0i+_cc4*_t1i+_cc3*_t2i+_cc1*_t3i+_cc5*_t4i+_cc2*_t5i;	// C4
	_cr5 = _t0r+_cc5*_t1r+_cc1*_t2r+_cc4*_t3r+_cc2*_t4r+_cc3*_t5r;	_ci5 = _t0i+_cc5*_t1i+_cc1*_t2i+_cc4*_t3i+_cc2*_t4i+_cc3*_t5i;	// C5
	__B0r = _t0r+_t1r+_t2r+_t3r+_t4r+_t5r;	__B0i = _t0i+_t1i+_t2i+_t3i+_t4i+_t5i;	// X0\
*/\
	/* In a 16-reg FMA model, makes sense to init the c-terms to the t0-summand,
	and keep the 10 c-terms and 4 of the 5 shared sincos data in registers. That uses 14 regs,
	leaving 2 for the 2 t*[r,i] data shared by each such block. Those need to be in-reg (at least
	in the cosine-mults section) because of the need to accumulate the DC terms: E.g. at end of the
	first block below (after doing the 10 c-term-update FMAs) we add the current values of B0r,i
	(which we init to t0r,i in the ASM version) to t1r,i, then write the result to memory and
	load t2r,i into the same 2 regs in preparation for the next such block.
	*/\
	_cr1 = _t0r + _cc1*_t1r;	_ci1 = _t0i + _cc1*_t1i;\
	_cr2 = _t0r + _cc2*_t1r;	_ci2 = _t0i + _cc2*_t1i;\
	_cr3 = _t0r + _cc3*_t1r;	_ci3 = _t0i + _cc3*_t1i;\
	_cr4 = _t0r + _cc4*_t1r;	_ci4 = _t0i + _cc4*_t1i;\
	_cr5 = _t0r + _cc5*_t1r;	_ci5 = _t0i + _cc5*_t1i;\
	__B0r = _t0r + _t1r;		__B0i = _t0i + _t1i;	\
\
	_cr1 += _cc2*_t2r;	_ci1 += _cc2*_t2i;\
	_cr2 += _cc4*_t2r;	_ci2 += _cc4*_t2i;\
	_cr3 += _cc5*_t2r;	_ci3 += _cc5*_t2i;\
	_cr4 += _cc3*_t2r;	_ci4 += _cc3*_t2i;\
	_cr5 += _cc1*_t2r;	_ci5 += _cc1*_t2i;\
	__B0r += _t2r;	__B0i += _t2i;\
\
	_cr1 += _cc3*_t3r;	_ci1 += _cc3*_t3i;\
	_cr2 += _cc5*_t3r;	_ci2 += _cc5*_t3i;\
	_cr3 += _cc2*_t3r;	_ci3 += _cc2*_t3i;\
	_cr4 += _cc1*_t3r;	_ci4 += _cc1*_t3i;\
	_cr5 += _cc4*_t3r;	_ci5 += _cc4*_t3i;\
	__B0r += _t3r;	__B0i += _t3i;\
\
	_cr1 += _cc4*_t4r;	_ci1 += _cc4*_t4i;\
	_cr2 += _cc3*_t4r;	_ci2 += _cc3*_t4i;\
	_cr3 += _cc1*_t4r;	_ci3 += _cc1*_t4i;\
	_cr4 += _cc5*_t4r;	_ci4 += _cc5*_t4i;\
	_cr5 += _cc2*_t4r;	_ci5 += _cc2*_t4i;\
	__B0r += _t4r;	__B0i += _t4i;\
\
	_cr1 += _cc5*_t5r;	_ci1 += _cc5*_t5i;\
	_cr2 += _cc1*_t5r;	_ci2 += _cc1*_t5i;\
	_cr3 += _cc4*_t5r;	_ci3 += _cc4*_t5i;\
	_cr4 += _cc2*_t5r;	_ci4 += _cc2*_t5i;\
	_cr5 += _cc3*_t5r;	_ci5 += _cc3*_t5i;\
	__B0r += _t5r;	__B0i += _t5i;\
\
/*\
	In a 16-reg FMA model, can keep the 10 s-terms and all 5 shared sincos data in registers, but that
	requires us to reload one of the 2 repeated t*-mults (e.g. _tai in the 1st block) 5 times per block.
	OTOH if we keep both t*r,i-mults and 4 of the 5 ss-mults in-reg, we just need 2 loads-from-mem per block:

	_sr1 =      _ss1*_tar+_ss2*_t9r+_ss3*_t8r+_ss4*_t7r+_ss5*_t6r;	_si1 =      _ss1*_tai+_ss2*_t9i+_ss3*_t8i+_ss4*_t7i+_ss5*_t6i;	// S1
	_sr2 =      _ss2*_tar+_ss4*_t9r-_ss5*_t8r-_ss3*_t7r-_ss1*_t6r;	_si2 =      _ss2*_tai+_ss4*_t9i-_ss5*_t8i-_ss3*_t7i-_ss1*_t6i;	// S2
	_sr3 =      _ss3*_tar-_ss5*_t9r-_ss2*_t8r+_ss1*_t7r+_ss4*_t6r;	_si3 =      _ss3*_tai-_ss5*_t9i-_ss2*_t8i+_ss1*_t7i+_ss4*_t6i;	// S3
	_sr4 =      _ss4*_tar-_ss3*_t9r+_ss1*_t8r+_ss5*_t7r-_ss2*_t6r;	_si4 =      _ss4*_tai-_ss3*_t9i+_ss1*_t8i+_ss5*_t7i-_ss2*_t6i;	// S4
	_sr5 =      _ss5*_tar-_ss1*_t9r+_ss4*_t8r-_ss2*_t7r+_ss3*_t6r;	_si5 =      _ss5*_tai-_ss1*_t9i+_ss4*_t8i-_ss2*_t7i+_ss3*_t6i;	// S5
*/\
	_sr1 = _ss1*_tar;	_si1 = _ss1*_tai;\
	_sr2 = _ss2*_tar;	_si2 = _ss2*_tai;\
	_sr3 = _ss3*_tar;	_si3 = _ss3*_tai;\
	_sr4 = _ss4*_tar;	_si4 = _ss4*_tai;\
	_sr5 = _ss5*_tar;	_si5 = _ss5*_tai;\
\
	_sr1 += _ss2*_t9r;	_si1 += _ss2*_t9i;\
	_sr2 += _ss4*_t9r;	_si2 += _ss4*_t9i;\
	_sr3 -= _ss5*_t9r;	_si3 -= _ss5*_t9i;\
	_sr4 -= _ss3*_t9r;	_si4 -= _ss3*_t9i;\
	_sr5 -= _ss1*_t9r;	_si5 -= _ss1*_t9i;\
\
	_sr1 += _ss3*_t8r;	_si1 += _ss3*_t8i;\
	_sr2 -= _ss5*_t8r;	_si2 -= _ss5*_t8i;\
	_sr3 -= _ss2*_t8r;	_si3 -= _ss2*_t8i;\
	_sr4 += _ss1*_t8r;	_si4 += _ss1*_t8i;\
	_sr5 += _ss4*_t8r;	_si5 += _ss4*_t8i;\
\
	_sr1 += _ss4*_t7r;	_si1 += _ss4*_t7i;\
	_sr2 -= _ss3*_t7r;	_si2 -= _ss3*_t7i;\
	_sr3 += _ss1*_t7r;	_si3 += _ss1*_t7i;\
	_sr4 += _ss5*_t7r;	_si4 += _ss5*_t7i;\
	_sr5 -= _ss2*_t7r;	_si5 -= _ss2*_t7i;\
\
	_sr1 += _ss5*_t6r;	_si1 += _ss5*_t6i;\
	_sr2 -= _ss1*_t6r;	_si2 -= _ss1*_t6i;\
	_sr3 += _ss4*_t6r;	_si3 += _ss4*_t6i;\
	_sr4 -= _ss2*_t6r;	_si4 -= _ss2*_t6i;\
	_sr5 += _ss3*_t6r;	_si5 += _ss3*_t6i;\
\
	__B1r = _cr1 - _si1;					__B1i = _ci1 + _sr1;	/* X1 = C1 + I*S1	*/\
	__B2r = _cr2 - _si2;					__B2i = _ci2 + _sr2;	/* X2 = C2 + I*S2	*/\
	__B3r = _cr3 - _si3;					__B3i = _ci3 + _sr3;	/* X3 = C3 + I*S3	*/\
	__B4r = _cr4 - _si4;					__B4i = _ci4 + _sr4;	/* X4 = C4 + I*S4	*/\
	__B5r = _cr5 - _si5;					__B5i = _ci5 + _sr5;	/* X5 = C5 + I*S5	*/\
	__B6r = _cr5 + _si5;					__B6i = _ci5 - _sr5;	/* X6 =	C5 - I*S5	*/\
	__B7r = _cr4 + _si4;					__B7i = _ci4 - _sr4;	/* X7 =	C4 - I*S4	*/\
	__B8r = _cr3 + _si3;					__B8i = _ci3 - _sr3;	/* X8 =	C3 - I*S3	*/\
	__B9r = _cr2 + _si2;					__B9i = _ci2 - _sr2;	/* X9 =	C2 - I*S2	*/\
	__Bar = _cr1 + _si1;					__Bai = _ci1 - _sr1;	/* X10=	C1 - I*S1	*/\
}

/*...Radix-11 DFT using length-5 cyclic convolution scheme for the 5 x 5 matrix submultiplies.  Totals: 160 ADD,  44 MUL	*/
//void ___RADIX_11_DFT() {}
#define RADIX_11_DFT(\
	__A0r,__A0i,__A1r,__A1i,__A2r,__A2i,__A3r,__A3i,__A4r,__A4i,__A5r,__A5i,__A6r,__A6i,__A7r,__A7i,__A8r,__A8i,__A9r,__A9i,__Aar,__Aai,\
	__B0r,__B0i,__B1r,__B1i,__B2r,__B2i,__B3r,__B3i,__B4r,__B4i,__B5r,__B5i,__B6r,__B6i,__B7r,__B7i,__B8r,__B8i,__B9r,__B9i,__Bar,__Bai,\
	_a0,_a1,_a2,_a3,_a4,_a5,_a6,_a7,_a8,_a9, _b0,_b1,_b2,_b3,_b4,_b5,_b6,_b7,_b8,_b9)\
{\
	double _t3,_t4,_t5,_t6,_t7,_t8,_t9,_t10,_t11,_t12,_t13,_t14,_t15,_t16,_t17,_t18,_t19,_t20,_t21,_t22;\
	double _cr1,_cr2,_cr3,_cr4,_cr5,_ci1,_ci2,_ci3,_ci4,_ci5,_sr1,_sr2,_sr3,_sr4,_sr5,_si1,_si2,_si3,_si4,_si5;\
	double _c1,_c2,_c3,_c4,_c5,_c6,_c7,_c8,_c9,_c10,_s1,_s2,_s3,_s4,_s5,_s6,_s7,_s8,_s9,_s10;\
\
	__B0r = __A0r;					__B0i = __A0i;			/* x0		; store in t2 */\
	_t3  = __A1r +__Aar;			_t4  = __A1i + __Aai;	/* x1 + x10	; load a1-a4,a7-aa into xmm0-7, compute t4-t10 (save in xmm0-3) */\
	_t5  = __A2r +__A9r;			_t6  = __A2i + __A9i;	/* x2 + x9	; and t16-22 (dump to memory), load a0,5,6 into xmm4-6, compute t12,14, */\
	_t7  = __A3r +__A8r;			_t8  = __A3i + __A8i;	/* x3 + x8	; dump t14 to memory. */\
	_t9  = __A4r +__A7r;			_t10 = __A4i + __A7i;	/* x4 + x7	*/\
	_t11 = __A5r +__A6r;			_t12 = __A5i + __A6i;	/* x5 + x6	*/\
	_t13 = __A5r -__A6r;			_t14 = __A5i - __A6i;	/* x5 - x6	*/\
	_t15 = __A4r -__A7r;			_t16 = __A4i - __A7i;	/* x4 - x7	*/\
	_t17 = __A3r -__A8r;			_t18 = __A3i - __A8i;	/* x3 - x8	*/\
	_t19 = __A2r -__A9r;			_t20 = __A2i - __A9i;	/* x2 - x9	*/\
	_t21 = __A1r -__Aar;			_t22 = __A1i - __Aai;	/* x1 - x10	*/\
\
/*  Here are the 5 cosine terms:\
	let y0 = (x1+x10) = t3,4, y4 = (x2+x9) = t5,6, y2 = (x3+x8) = t7,8, y3 = (x4+x7) = t9,10, y1 = (x5+x6) = t11,12, then form\
*/\
	/* b0/t2,t4,t6,t8,t10,t12 in xmm0-5, xmm6,7 free */\
	_c1 = _t3 -_t5;					_s1 = _t4 -_t6;			/* store in s1/t4,  make 2 copies (call them t0,t14) */\
	_c2 = _t11-_t5;					_s2 = _t12-_t6;			/* store in s2/t12  */\
	_c4 = _t7 -_t5;					_s4 = _t8 -_t6;			/* store in s4/t8 */\
	_c5 = _t9 -_t5;					_s5 = _t10-_t6;			/* store in s5/t10, then mul t6*5 */\
	_c3 = _c1 +_c2;					_s3 = _s1 +_s2;			/* t4+t12-2*t6; store in s3/t0 */\
	_c7 = _c1 -_c4;					_s7 = _s1 -_s4;			/* t4-t8; store in t14, mul by a6 and store. */\
	_c6 = _c4 +_c5;					_s6 = _s4 +_s5;			/* t8+t10-2*t6; store in t8 (do after s1-s4, since this will overwrite s4) */\
	_c8 = _c2 -_c5;					_s8 = _s2 -_s5;			/* t12-t10; copy s2 to t14, mul s2-s5 by a7 and store. */\
	_c9 = _c3 -_c6;					_s9 = _s3 -_s6;			/* t4+t12-(t8+t10); copy s3 to t14, store result in t14. */\
	_c10= _c3 +_c6+5*_t5;			_s10= _s3 +_s6+5*_t6;	/* c10= t3+t5+t7+t9+t11;	 s10= t4+t6+t8+t10+t12; store in t6. */\
\
	__B0r += _c10;				__B0i += _s10;	/* X0/t2, store result to memory */\
\
	_c3 = _a2*_c3;					_s3 = _a2*_s3;	/* t0  */\
	_c6 = _a5*_c6;					_s6 = _a5*_s6;	/* t8  */\
	_c9 = _a8*_c9;					_s9 = _a8*_s9;	/* t14 */\
	_c10= _a9*_c10+__B0r;			_s10= _a9*_s10+__B0i;	/* t6  */\
\
	_c1 = _a0*_c1+_c3;				_s1 = _a0*_s1+_s3;\
	_c2 = _a1*_c2+_c3;				_s2 = _a1*_s2+_s3;\
	_c4 = _a3*_c4+_c6;				_s4 = _a3*_s4+_s6;\
	_c5 = _a4*_c5+_c6;				_s5 = _a4*_s5+_s6;\
	_c7 = _a6*_c7+_c9;				_s7 = _a6*_s7+_s9;\
	_c8 = _a7*_c8+_c9;				_s8 = _a7*_s8+_s9;\
\
	_cr1 = _c10+_c1-_c7;			_ci1 = _s10+_s1-_s7;\
	_cr2 = _c10-_c1-_c2-_c4-_c5;	_ci2 = _s10-_s1-_s2-_s4-_s5;\
	_cr3 = _c10+_c4+_c7;			_ci3 = _s10+_s4+_s7;\
	_cr4 = _c10+_c5+_c8;			_ci4 = _s10+_s5+_s8;\
	_cr5 = _c10+_c2-_c8;			_ci5 = _s10+_s2-_s8;\
\
/*  Here are the 5 sine terms:\
let y0 = (x1-x10) = t21,22, y4 = (x9-x2) = -t19,20, y2 = (x3-x8) = t17,18, y3 = (x4-x7) = t15,16, y1 = (x5-x6) = t13,t14, then form\
*/\
	_c1 = _t21+_t19;				_s1 = _t22+_t20;\
	_c2 = _t13+_t19;				_s2 = _t14+_t20;\
	_c3 = _c1+_c2;					_s3 = _s1+_s2;\
	_c4 = _t17+_t19;				_s4 = _t18+_t20;\
	_c5 = _t15+_t19;				_s5 = _t16+_t20;\
	_c6 = _c4+_c5;					_s6 = _s4+_s5;\
	_c7 = _c1-_c4;					_s7 = _s1-_s4;\
	_c8 = _c2-_c5;					_s8 = _s2-_s5;\
	_c9 = _c3-_c6;					_s9 = _s3-_s6;\
	_c10= _c3+_c6-5*_t19;			_s10= _s3+_s6-5*_t20;	/* _c10= _t21-_t19+_t17+_t15+_t13;	_s10= _t22-_t20+_t18+_t16+_t14; */\
\
	_c3 = _b2*_c3;					_s3 = _b2*_s3;\
	_c6 = _b5*_c6;					_s6 = _b5*_s6;\
	_c9 = _b8*_c9;					_s9 = _b8*_s9;\
	_c10= _b9*_c10;					_s10= _b9*_s10;\
\
	_c1 = _b0*_c1+_c3;				_s1 = _b0*_s1+_s3;\
	_c2 = _b1*_c2+_c3;				_s2 = _b1*_s2+_s3;\
	_c3 = _b3*_c4+_c6;				_s3 = _b3*_s4+_s6;\
	_c4 = _b4*_c5+_c6;				_s4 = _b4*_s5+_s6;\
	_c5 = _b6*_c7+_c9;				_s5 = _b6*_s7+_s9;\
	_c6 = _b7*_c8+_c9;				_s6 = _b7*_s8+_s9;\
\
	_sr1 = _c10+_c1-_c5;			_si1 = _s10+_s1-_s5;\
	_sr2 = _c1+_c2+_c3+_c4-_c10;	_si2 = _s1+_s2+_s3+_s4-_s10;\
	_sr3 = _c10+_c3+_c5;			_si3 = _s10+_s3+_s5;\
	_sr4 = _c10+_c4+_c6;			_si4 = _s10+_s4+_s6;\
	_sr5 = _c10+_c2-_c6;			_si5 = _s10+_s2-_s6;\
\
	__B1r = _cr1 - _si1;			__B1i = _ci1 + _sr1;	/* X1 = C1 + I*S1	*/\
	__B2r = _cr2 - _si2;			__B2i = _ci2 + _sr2;	/* X2 = C2 + I*S2	*/\
	__B3r = _cr3 - _si3;			__B3i = _ci3 + _sr3;	/* X3 = C3 + I*S3	*/\
	__B4r = _cr4 - _si4;			__B4i = _ci4 + _sr4;	/* X4 = C4 + I*S4	*/\
	__B5r = _cr5 - _si5;			__B5i = _ci5 + _sr5;	/* X5 = C5 + I*S5	*/\
	__B6r = _cr5 + _si5;			__B6i = _ci5 - _sr5;	/* X6 =	C5 - I*S5	*/\
	__B7r = _cr4 + _si4;			__B7i = _ci4 - _sr4;	/* X7 =	C4 - I*S4	*/\
	__B8r = _cr3 + _si3;			__B8i = _ci3 - _sr3;	/* X8 =	C3 - I*S3	*/\
	__B9r = _cr2 + _si2;			__B9i = _ci2 - _sr2;	/* X9 =	C2 - I*S2	*/\
	__Bar = _cr1 + _si1;			__Bai = _ci1 - _sr1;	/* X10=	C1 - I*S1	*/\
}

/****** RADIX = 13: ******/

//...Simple-algo Radix-13 DFT, which also works well in an FMA context.  Totals: 140 ADD, 100 MUL, or 40 ADD, 100 FMA.
//
//void ___RADIX_13_DFT_BASIC() {}
#define RADIX_13_DFT_BASIC(\
	__A0r,__A0i,__A1r,__A1i,__A2r,__A2i,__A3r,__A3i,__A4r,__A4i,__A5r,__A5i,__A6r,__A6i,__A7r,__A7i,__A8r,__A8i,__A9r,__A9i,__Aar,__Aai,__Abr,__Abi,__Acr,__Aci,\
	__B0r,__B0i,__B1r,__B1i,__B2r,__B2i,__B3r,__B3i,__B4r,__B4i,__B5r,__B5i,__B6r,__B6i,__B7r,__B7i,__B8r,__B8i,__B9r,__B9i,__Bar,__Bai,__Bbr,__Bbi,__Bcr,__Bci,\
	_cc1,_cc2,_cc3,_cc4,_cc5,_cc6,_ss1,_ss2,_ss3,_ss4,_ss5,_ss6)\
{\
	double _t0r,_t0i,_t1r,_t1i,_t2r,_t2i,_t3r,_t3i,_t4r,_t4i,_t5r,_t5i,_t6r,_t6i,_t7r,_t7i,_t8r,_t8i,_t9r,_t9i,_tar,_tai,_tbr,_tbi,_tcr,_tci;\
	double _cr1,_cr2,_cr3,_cr4,_cr5,_cr6,_ci1,_ci2,_ci3,_ci4,_ci5,_ci6,_sr1,_sr2,_sr3,_sr4,_sr5,_sr6,_si1,_si2,_si3,_si4,_si5,_si6;\
\
	_t0r = __A0r       ;					_t0i = __A0i        ;	/* x0		*/\
	_t1r = __A1r +__Acr;					_t1i = __A1i + __Aci;	/* x1 + x12	*/\
	_t2r = __A2r +__Abr;					_t2i = __A2i + __Abi;	/* x2 + x11	*/\
	_t3r = __A3r +__Aar;					_t3i = __A3i + __Aai;	/* x3 + x10	*/\
	_t4r = __A4r +__A9r;					_t4i = __A4i + __A9i;	/* x4 + x9 	*/\
	_t5r = __A5r +__A8r;					_t5i = __A5i + __A8i;	/* x5 + x8 	*/\
	_t6r = __A6r +__A7r;					_t6i = __A6i + __A7i;	/* x6 - x7 	*/\
\
	_tcr = __A1r -__Acr;					_tci = __A1i - __Aci;	/* x1 - x12	*/\
	_tbr = __A2r -__Abr;					_tbi = __A2i - __Abi;	/* x2 - x11	*/\
	_tar = __A3r -__Aar;					_tai = __A3i - __Aai;	/* x3 - x10	*/\
	_t9r = __A4r -__A9r;					_t9i = __A4i - __A9i;	/* x4 - x9 	*/\
	_t8r = __A5r -__A8r;					_t8i = __A5i - __A8i;	/* x5 - x8 	*/\
	_t7r = __A6r -__A7r;					_t7i = __A6i - __A7i;	/* x6 - x7 	*/\
\
	_cr1 = _t0r+_cc1*_t1r+_cc2*_t2r+_cc3*_t3r+_cc4*_t4r+_cc5*_t5r+_cc6*_t6r;	_ci1 = _t0i+_cc1*_t1i+_cc2*_t2i+_cc3*_t3i+_cc4*_t4i+_cc5*_t5i+_cc6*_t6i;\
	_cr2 = _t0r+_cc2*_t1r+_cc4*_t2r+_cc6*_t3r+_cc5*_t4r+_cc3*_t5r+_cc1*_t6r;	_ci2 = _t0i+_cc2*_t1i+_cc4*_t2i+_cc6*_t3i+_cc5*_t4i+_cc3*_t5i+_cc1*_t6i;\
	_cr3 = _t0r+_cc3*_t1r+_cc6*_t2r+_cc4*_t3r+_cc1*_t4r+_cc2*_t5r+_cc5*_t6r;	_ci3 = _t0i+_cc3*_t1i+_cc6*_t2i+_cc4*_t3i+_cc1*_t4i+_cc2*_t5i+_cc5*_t6i;\
	_cr4 = _t0r+_cc4*_t1r+_cc5*_t2r+_cc1*_t3r+_cc3*_t4r+_cc6*_t5r+_cc2*_t6r;	_ci4 = _t0i+_cc4*_t1i+_cc5*_t2i+_cc1*_t3i+_cc3*_t4i+_cc6*_t5i+_cc2*_t6i;\
	_cr5 = _t0r+_cc5*_t1r+_cc3*_t2r+_cc2*_t3r+_cc6*_t4r+_cc1*_t5r+_cc4*_t6r;	_ci5 = _t0i+_cc5*_t1i+_cc3*_t2i+_cc2*_t3i+_cc6*_t4i+_cc1*_t5i+_cc4*_t6i;\
	_cr6 = _t0r+_cc6*_t1r+_cc1*_t2r+_cc5*_t3r+_cc2*_t4r+_cc4*_t5r+_cc3*_t6r;	_ci6 = _t0i+_cc6*_t1i+_cc1*_t2i+_cc5*_t3i+_cc2*_t4i+_cc4*_t5i+_cc3*_t6i;\
\
	_sr1 =      _ss1*_tcr+_ss2*_tbr+_ss3*_tar+_ss4*_t9r+_ss5*_t8r+_ss6*_t7r;	_si1 =      _ss1*_tci+_ss2*_tbi+_ss3*_tai+_ss4*_t9i+_ss5*_t8i+_ss6*_t7i;	/* S1 */\
	_sr2 =      _ss2*_tcr+_ss4*_tbr+_ss6*_tar-_ss5*_t9r-_ss3*_t8r-_ss1*_t7r;	_si2 =      _ss2*_tci+_ss4*_tbi+_ss6*_tai-_ss5*_t9i-_ss3*_t8i-_ss1*_t7i;	/* S2 */\
	_sr3 =      _ss3*_tcr+_ss6*_tbr-_ss4*_tar-_ss1*_t9r+_ss2*_t8r+_ss5*_t7r;	_si3 =      _ss3*_tci+_ss6*_tbi-_ss4*_tai-_ss1*_t9i+_ss2*_t8i+_ss5*_t7i;	/* S3 */\
	_sr4 =      _ss4*_tcr-_ss5*_tbr-_ss1*_tar+_ss3*_t9r-_ss6*_t8r-_ss2*_t7r;	_si4 =      _ss4*_tci-_ss5*_tbi-_ss1*_tai+_ss3*_t9i-_ss6*_t8i-_ss2*_t7i;	/* S4 */\
	_sr5 =      _ss5*_tcr-_ss3*_tbr+_ss2*_tar-_ss6*_t9r-_ss1*_t8r+_ss4*_t7r;	_si5 =      _ss5*_tci-_ss3*_tbi+_ss2*_tai-_ss6*_t9i-_ss1*_t8i+_ss4*_t7i;	/* S5 */\
	_sr6 =      _ss6*_tcr-_ss1*_tbr+_ss5*_tar-_ss2*_t9r+_ss4*_t8r-_ss3*_t7r;	_si6 =      _ss6*_tci-_ss1*_tbi+_ss5*_tai-_ss2*_t9i+_ss4*_t8i-_ss3*_t7i;	/* S5 */\
\
	__B0r = _t0r+_t1r+_t2r+_t3r+_t4r+_t5r+_t6r;	__B0i = _t0i+_t1i+_t2i+_t3i+_t4i+_t5i+_t6i;	/* X0	*/\
	__B1r = _cr1 - _si1;						__B1i = _ci1 + _sr1;	/* X1 = C1 + I*S1	*/\
	__B2r = _cr2 - _si2;						__B2i = _ci2 + _sr2;	/* X2 = C2 + I*S2	*/\
	__B3r = _cr3 - _si3;						__B3i = _ci3 + _sr3;	/* X3 = C3 + I*S3	*/\
	__B4r = _cr4 - _si4;						__B4i = _ci4 + _sr4;	/* X4 = C4 + I*S4	*/\
	__B5r = _cr5 - _si5;						__B5i = _ci5 + _sr5;	/* X5 = C5 + I*S5	*/\
	__B6r = _cr6 - _si6;						__B6i = _ci6 + _sr6;	/* X6 = C6 + I*S6	*/\
\
	__Bcr = _cr1 + _si1;						__Bci = _ci1 - _sr1;	/* X12=	C1 - I*S1	*/\
	__Bbr = _cr2 + _si2;						__Bbi = _ci2 - _sr2;	/* X11=	C2 - I*S2	*/\
	__Bar = _cr3 + _si3;						__Bai = _ci3 - _sr3;	/* X10=	C3 - I*S3	*/\
	__B9r = _cr4 + _si4;						__B9i = _ci4 - _sr4;	/* X9 =	C4 - I*S4	*/\
	__B8r = _cr5 + _si5;						__B8i = _ci5 - _sr5;	/* X8 =	C5 - I*S5	*/\
	__B7r = _cr6 + _si6;						__B7i = _ci6 - _sr6;	/* X7 =	C6 - I*S6	*/\
}

/* This version implements a van Buskirk tangent-based radix-13 DFT and only tries to minimize the number
of temporaries used in each of the major phases, with view toward a 16-register-optimized version.
(That is why we do the x and y-pieces side-by-side).
For an 8-register-optimized version which moreover mimics x86-style destructive (2-input, one overwritten with output)
register arithmetic, see the RADIX_13_XYZ version of the same macro in radix13_ditN_cy_dif1.c .

Output ordering is as for DIF, so for DIT caller must reverse order of last 12 B-output args.

ASSUMES that the needed DC and DS-constants have been defined in the calling routine via
inclusion of the radix13.h header file..
*/
//void ___RADIX_13_DFT() {}
#define RADIX_13_DFT(\
	__A0r,__A0i,__A1r,__A1i,__A2r,__A2i,__A3r,__A3i,__A4r,__A4i,__A5r,__A5i,__A6r,__A6i,__A7r,__A7i,__A8r,__A8i,__A9r,__A9i,__Aar,__Aai,__Abr,__Abi,__Acr,__Aci,\
	__B0r,__B0i,__B1r,__B1i,__B2r,__B2i,__B3r,__B3i,__B4r,__B4i,__B5r,__B5i,__B6r,__B6i,__B7r,__B7i,__B8r,__B8i,__B9r,__B9i,__Bar,__Bai,__Bbr,__Bbi,__Bcr,__Bci)\
{\
	double _xr0, _xi0, _xr1, _xi1, _xr2, _xi2, _xr3, _xi3, _xr4, _xi4, _xr5, _xi5, _xr6, _xi6, _xr7, _xi7;\
	double _yr1, _yi1, _yr2, _yi2, _yr3, _yi3, _yr4, _yi4, _yr5, _yi5, _yr6, _yi6;\
\
		_yr1 = __A1r - __Acr;			_yi1 = __A1i - __Aci;\
		_yr2 = __A2r - __Abr;			_yi2 = __A2i - __Abi;\
		_yr5 = __A3r - __Aar;			_yi5 = __A3i - __Aai;\
		_yr3 = __A4r - __A9r;			_yi3 = __A4i - __A9i;\
		_yr4 = __A5r - __A8r;			_yi4 = __A5i - __A8i;\
		_yr6 = __A6r - __A7r;			_yi6 = __A6i - __A7i;\
\
		_xr1 = __A1r + __Acr;			_xi1 = __A1i + __Aci;\
		_xr3 = __A2r + __Abr;			_xi3 = __A2i + __Abi;\
		_xr4 = __A4r + __A9r;			_xi4 = __A4i + __A9i;\
		_xr5 = __A5r + __A8r;			_xi5 = __A5i + __A8i;\
		_xr6 = __A3r + __Aar;			_xi6 = __A3i + __Aai;\
		_xr7 = __A6r + __A7r;			_xi7 = __A6i + __A7i;\
\
/* x-terms need 8 registers for each side: */\
		_xr2 = _xr1+_xr5;					_xi2 = _xi1+_xi5;\
		_xr5 = _xr1-_xr5;					_xi5 = _xi1-_xi5;\
		_xr0 = _xr3+_xr6;					_xi0 = _xi3+_xi6;\
		_xr6 = _xr3-_xr6;					_xi6 = _xi3-_xi6;\
		_xr3 = _xr4+_xr7;					_xi3 = _xi4+_xi7;\
		_xr7 = _xr4-_xr7;					_xi7 = _xi4-_xi7;\
		_xr1 = _xr2+_xr0+_xr3;				_xi1 = _xi2+_xi0+_xi3;\
		_xr4 = __A0r+_xr1;					_xi4 = __A0i+_xi1;\
		__B0r = _xr4;						__B0i = _xi4;\
		_xr4 +=    DC1*_xr1;				_xi4 +=    DC1*_xi1;\
		_xr1 = _xr4+DC2*_xr2+DC3*_xr0;		_xi1 = _xi4+DC2*_xi2+DC3*_xi0;\
		_xr2 = _xr4+DC2*_xr3+DC3*_xr2;		_xi2 = _xi4+DC2*_xi3+DC3*_xi2;\
		_xr0 = _xr4+DC2*_xr0+DC3*_xr3;		_xi0 = _xi4+DC2*_xi0+DC3*_xi3;\
		_xr4 = DC4*_xr5+DC5*_xr6+DC6*_xr7;	_xi4 = DC4*_xi5+DC5*_xi6+DC6*_xi7;\
		_xr3 = DC4*_xr6+DC5*_xr7-DC6*_xr5;	_xi3 = DC4*_xi6+DC5*_xi7-DC6*_xi5;\
		_xr7 = DC4*_xr7-DC5*_xr5-DC6*_xr6;	_xi7 = DC4*_xi7-DC5*_xi5-DC6*_xi6;\
		_xr5 = _xr0+_xr4;					_xi5 = _xi0+_xi4;\
		_xr0 = _xr0-_xr4;					_xi0 = _xi0-_xi4;\
		_xr6 = _xr2+_xr3;					_xi6 = _xi2+_xi3;\
		_xr2 = _xr2-_xr3;					_xi2 = _xi2-_xi3;\
		_xr3 = _xr1+_xr7;					_xi3 = _xi1+_xi7;\
		_xr7 = _xr1-_xr7;					_xi7 = _xi1-_xi7;\
/* _x4,1 free...now do _y-terms: */\
		_xr4 = _yr1-_yr3+_yr5;				_xi4 = _yi1-_yi3+_yi5;\
		_yr1 = _yr1+_yr3;					_yi1 = _yi1+_yi3;\
		_yr5 = _yr3+_yr5;		 			_yi5 = _yi3+_yi5;\
		_yr3 = _yr2+_yr4+_yr6;				_yi3 = _yi2+_yi4+_yi6;\
		_yr2 = _yr2-_yr4;					_yi2 = _yi2-_yi4;\
		_yr4 = _yr6-_yr4;		 			_yi4 = _yi6-_yi4;\
		_xr1 = DS1*_xr4+DS2*_yr3;			_xi1 = DS1*_xi4+DS2*_yi3;\
/* Need to use one output as a temporary here if have just 8 registers: */\
		__B1i = DS1*_yr3-DS2*_xr4;			__B1r = DS1*_yi3-DS2*_xi4;\
		_yr3 = _yr1+_yr2;		 			_yi3 = _yi1+_yi2;\
		_yr6 = _yr5+_yr4;					_yi6 = _yi5+_yi4;\
		_xr4 = DS3*_yr3-DS6*_yr6;			_xi4 = DS3*_yi3-DS6*_yi6;\
		_yr3 = DS9*_yr3-DS3*_yr6;			_yi3 = DS9*_yi3-DS3*_yi6;\
		_yr6 = _xr4+DS4*_yr2-DS7*_yr4;		_yi6 = _xi4+DS4*_yi2-DS7*_yi4;\
		_yr2 = _yr3+DSa*_yr2-DS4*_yr4;		_yi2 = _yi3+DSa*_yi2-DS4*_yi4;\
		_yr4 = _xr4+DS5*_yr1-DS8*_yr5;		_yi4 = _xi4+DS5*_yi1-DS8*_yi5;\
		_yr3 = _yr3+DSb*_yr1-DS5*_yr5;		_yi3 = _yi3+DSb*_yi1-DS5*_yi5;\
		_xr4 = __B1i;						_xi4 = __B1r;\
		_yr5 = _xr1+_yr6;					_yi5 = _xi1+_yi6;\
		_yr6 = _yr6-_xr1-_yr2;				_yi6 = _yi6-_xi1-_yi2;\
		_yr2 = _xr1-_yr2;					_yi2 = _xi1-_yi2;\
		_yr1 = _xr4+_yr4;					_yi1 = _xi4+_yi4;\
		_yr4 = _yr4-_xr4-_yr3;				_yi4 = _yi4-_xi4-_yi3;\
		_xr4 = _xr4-_yr3;					_xi4 = _xi4-_yi3;\
/* In ASM, do xr and yi-terms and combine to get real parts of outputs,\
   then do xi and yr-terms and combine to get imaginary parts: */\
		__B1r = _xr7-_xi4;				__B1i = _xi7+_xr4;\
		__Bcr = _xr7+_xi4;	 			__Bci = _xi7-_xr4;\
		__B2r = _xr2-_yi2;	 			__B2i = _xi2+_yr2;\
		__Bbr = _xr2+_yi2;				__Bbi = _xi2-_yr2;\
		__B3r = _xr6-_yi1;	 			__B3i = _xi6+_yr1;\
		__Bar = _xr6+_yi1;				__Bai = _xi6-_yr1;\
		__B4r = _xr0-_yi4;	 			__B4i = _xi0+_yr4;\
		__B9r = _xr0+_yi4;				__B9i = _xi0-_yr4;\
		__B5r = _xr3+_yi6;				__B5i = _xi3-_yr6;\
		__B8r = _xr3-_yi6;				__B8i = _xi3+_yr6;\
		__B6r = _xr5-_yi5;				__B6i = _xi5+_yr5;\
		__B7r = _xr5+_yi5;				__B7i = _xi5-_yr5;\
/* Totals: 164 ADD, 64 MUL. */\
}

/****** RADIX = 15: ******/

/* Totals: 162 ADD, 50 MUL	*/
//void ___RADIX_15_DIF() {}
#define RADIX_15_DIF(\
	__Ar00,__Ai00,__Ar01,__Ai01,__Ar02,__Ai02,__Ar03,__Ai03,__Ar04,__Ai04,\
	__Ar05,__Ai05,__Ar06,__Ai06,__Ar07,__Ai07,__Ar08,__Ai08,__Ar09,__Ai09,\
	__Ar10,__Ai10,__Ar11,__Ai11,__Ar12,__Ai12,__Ar13,__Ai13,__Ar14,__Ai14,\
	__Br00,__Bi00,__Br01,__Bi01,__Br02,__Bi02,__Br03,__Bi03,__Br04,__Bi04,\
	__Br05,__Bi05,__Br06,__Bi06,__Br07,__Bi07,__Br08,__Bi08,__Br09,__Bi09,\
	__Br10,__Bi10,__Br11,__Bi11,__Br12,__Bi12,__Br13,__Bi13,__Br14,__Bi14)\
{\
	double _tr00,_ti00,_tr01,_ti01,_tr02,_ti02,_tr03,_ti03,_tr04,_ti04,_tr05,_ti05,_tr06,_ti06,_tr07,_ti07,_tr08,_ti08,_tr09,_ti09,_tr10,_ti10,_tr11,_ti11,_tr12,_ti12,_tr13,_ti13,_tr14,_ti14;\
	double _rt,_it;\
/* Default input order for a twiddleless radix-15 DIF (cf.radix15_dif_pass1()) is [0,10,5,12,7,2,9,4,14,6,1,11,3,13,8]. */\
/*...gather the needed data (15 64-bit complex, i.e. 30 64-bit reals) and do 3 radix-5 transforms...*/\
/*...Block 1: Swap inputs as 03 <-> 12, 06 <-> 09: */\
	_tr00 = __Ar00;				_ti00 = __Ai00;\
	_tr01 = __Ar12;				_ti01 = __Ai12;\
	_rt   = __Ar03;				_it   = __Ai03;\
	_tr03 = _tr01-_rt;			_ti03 = _ti01-_it;\
	_tr01 = _tr01+_rt;			_ti01 = _ti01+_it;\
	_tr02 = __Ar09;				_ti02 = __Ai09;\
	_rt   = __Ar06;				_it   = __Ai06;\
	_tr04 = _tr02-_rt;			_ti04 = _ti02-_it;\
	_tr02 = _tr02+_rt;			_ti02 = _ti02+_it;\
\
	_rt   = _tr01+_tr02;		_it   = _ti01+_ti02;\
	_tr00 = _tr00+_rt;			_ti00 = _ti00+_it;\
	_rt   = _tr00+cn1*_rt;		_it   = _ti00+cn1*_it;\
	_tr02 = cn2*(_tr01-_tr02);	_ti02 = cn2*(_ti01-_ti02);\
	_tr01 = _rt+_tr02;			_ti01 = _it+_ti02;\
	_tr02 = _rt-_tr02;			_ti02 = _it-_ti02;\
	_rt   = ss3*(_tr03-_tr04);	_it   = ss3*(_ti03-_ti04);\
	_tr04 = _rt+sn1*_tr04;		_ti04 = _it+sn1*_ti04;\
	_tr03 = _rt-sn2*_tr03;		_ti03 = _it-sn2*_ti03;\
	_rt   = _tr04;				_it   = _ti04;\
	_tr04 = _tr01+_it;			_ti04 = _ti01-_rt;		/*<==prefer these to be stored in _tr,i04 */\
	_tr01 = _tr01-_it;			_ti01 = _ti01+_rt;\
	_rt   = _tr03;				_it   = _ti03;\
	_tr03 = _tr02+_it;			_ti03 = _ti02-_rt;		/*<==prefer these to be stored in _tr,i03 */\
	_tr02 = _tr02-_it;			_ti02 = _ti02+_rt;\
\
/*...Block 2: Swap inputs as 01 <-> 10, 04 <-> 07: */\
	_tr05 = __Ar10;				_ti05 = __Ai10;\
	_tr06 = __Ar07;				_ti06 = __Ai07;\
	_rt   = __Ar13;				_it   = __Ai13;\
	_tr08 = _tr06 -_rt;			_ti08 = _ti06 -_it;\
	_tr06 = _tr06 +_rt;			_ti06 = _ti06 +_it;\
	_tr07 = __Ar04;				_ti07 = __Ai04;\
	_rt   = __Ar01;				_it   = __Ai01;\
	_tr09 = _tr07 -_rt;			_ti09 = _ti07 -_it;\
	_tr07 = _tr07 +_rt;			_ti07 = _ti07 +_it;\
\
	_rt   = _tr06+_tr07;		_it   = _ti06+_ti07;\
	_tr05 = _tr05+_rt;			_ti05 = _ti05+_it;\
	_rt   = _tr05+cn1*_rt;		_it   = _ti05+cn1*_it;\
	_tr07 = cn2*(_tr06-_tr07);	_ti07 = cn2*(_ti06-_ti07);\
	_tr06 = _rt+_tr07;			_ti06 = _it+_ti07;\
	_tr07 = _rt-_tr07;			_ti07 = _it-_ti07;\
	_rt   = ss3*(_tr08-_tr09);	_it   = ss3*(_ti08-_ti09);\
	_tr09 = _rt+sn1*_tr09;		_ti09 = _it+sn1*_ti09;\
	_tr08 = _rt-sn2*_tr08;		_ti08 = _it-sn2*_ti08;\
	_rt   = _tr09;				_it   = _ti09;\
	_tr09 = _tr06+_it;			_ti09 = _ti06-_rt;\
	_tr06 = _tr06-_it;			_ti06 = _ti06+_rt;\
	_rt   = _tr08;				_it   = _ti08;\
	_tr08 = _tr07+_it;			_ti08 = _ti07-_rt;\
	_tr07 = _tr07-_it;			_ti07 = _ti07+_rt;\
\
/*...Block 3: Swap inputs as 02 <-> 05, 08 <-> 14: */\
	_tr10 = __Ar05;				_ti10 = __Ai05;\
	_tr11 = __Ar02;				_ti11 = __Ai02;\
	_rt   = __Ar08;				_it   = __Ai08;\
	_tr13 = _tr11 -_rt;			_ti13 = _ti11 -_it;\
	_tr11 = _tr11 +_rt;			_ti11 = _ti11 +_it;\
	_tr12 = __Ar14;				_ti12 = __Ai14;\
	_rt   = __Ar11;				_it   = __Ai11;\
	_tr14 = _tr12 -_rt;			_ti14 = _ti12 -_it;\
	_tr12 = _tr12 +_rt;			_ti12 = _ti12 +_it;\
\
	_rt   = _tr11+_tr12;		_it   = _ti11+_ti12;\
	_tr10 = _tr10+_rt;			_ti10 = _ti10+_it;\
	_rt   = _tr10+cn1*_rt;		_it   = _ti10+cn1*_it;\
	_tr12 = cn2*(_tr11-_tr12);	_ti12 = cn2*(_ti11-_ti12);\
	_tr11 = _rt+_tr12;			_ti11 = _it+_ti12;\
	_tr12 = _rt-_tr12;			_ti12 = _it-_ti12;\
	_rt   = ss3*(_tr13-_tr14);	_it   = ss3*(_ti13-_ti14);\
	_tr14 = _rt+sn1*_tr14;		_ti14 = _it+sn1*_ti14;\
	_tr13 = _rt-sn2*_tr13;		_ti13 = _it-sn2*_ti13;\
	_rt   = _tr14;				_it   = _ti14;\
	_tr14 = _tr11+_it;			_ti14 = _ti11-_rt;\
	_tr11 = _tr11-_it;			_ti11 = _ti11+_rt;\
	_rt   = _tr13;				_it   = _ti13;\
	_tr13 = _tr12+_it;			_ti13 = _ti12-_rt;\
	_tr12 = _tr12-_it;			_ti12 = _ti12+_rt;\
\
/*...and now do five radix-3 transforms.\
	The required output permutation is [0,1,2,13,14,12,9,10,11,8,6,7,4,5,3]. */\
/*...Block 1: */\
	_rt   = _tr10;				_it   = _ti10;\
	_tr10 = _tr05-_rt;			_ti10 = _ti05-_it;\
	_tr05 = _tr05+_rt;			_ti05 = _ti05+_it;\
	_tr00 = _tr00+_tr05;		_ti00 = _ti00+_ti05;\
	__Br00 = _tr00;				__Bi00 = _ti00;\
	_tr05 = _tr00+c3m1*_tr05;	_ti05 = _ti00+c3m1*_ti05;\
	_rt   = s*_tr10;			_it   = s*_ti10;\
	__Br01 = _tr05-_it;			__Bi01 = _ti05+_rt;\
	__Br02 = _tr05+_it;			__Bi02 = _ti05-_rt;\
\
/*...Block 2: */\
	_rt   = _tr11;				_it   = _ti11;\
	_tr11 = _tr06-_rt;			_ti11 = _ti06-_it;\
	_tr06 = _tr06+_rt;			_ti06 = _ti06+_it;\
	_tr01 = _tr01+_tr06;		_ti01 = _ti01+_ti06;\
	__Br13 = _tr01;				__Bi13 = _ti01;\
	_tr06 = _tr01+c3m1*_tr06;	_ti06 = _ti01+c3m1*_ti06;\
	_rt   = s*_tr11;			_it   = s*_ti11;\
	__Br14 = _tr06-_it;			__Bi14 = _ti06+_rt;\
	__Br12 = _tr06+_it;			__Bi12 = _ti06-_rt;\
\
/*...Block 3: */\
	_rt   = _tr12;				_it   = _ti12;\
	_tr12 = _tr07-_rt;			_ti12 = _ti07-_it;\
	_tr07 = _tr07+_rt;			_ti07 = _ti07+_it;\
	_tr02 = _tr02+_tr07;		_ti02 = _ti02+_ti07;\
	__Br09 = _tr02;				__Bi09 = _ti02;\
	_tr07 = _tr02+c3m1*_tr07;	_ti07 = _ti02+c3m1*_ti07;\
	_rt   = s*_tr12;			_it   = s*_ti12;\
	__Br10 = _tr07-_it;			__Bi10 = _ti07+_rt;\
	__Br11 = _tr07+_it;			__Bi11 = _ti07-_rt;\
\
/*...Block 4: */\
	_rt   = _tr13;				_it   = _ti13;\
	_tr13 = _tr08-_rt;			_ti13 = _ti08-_it;\
	_tr08 = _tr08+_rt;			_ti08 = _ti08+_it;\
	_tr03 = _tr03+_tr08;		_ti03 = _ti03+_ti08;\
	__Br08 = _tr03;				__Bi08 = _ti03;\
	_tr08 = _tr03+c3m1*_tr08;	_ti08 = _ti03+c3m1*_ti08;\
	_rt   = s*_tr13;			_it   = s*_ti13;\
	__Br06 = _tr08-_it;			__Bi06 = _ti08+_rt;\
	__Br07 = _tr08+_it;			__Bi07 = _ti08-_rt;\
\
/*...Block 5: */\
	_rt   = _tr14;				_it   = _ti14;\
	_tr14 = _tr09-_rt;			_ti14 = _ti09-_it;\
	_tr09 = _tr09+_rt;			_ti09 = _ti09+_it;\
	_tr04 = _tr04+_tr09;		_ti04 = _ti04+_ti09;\
	__Br04 = _tr04;				__Bi04 = _ti04;\
	_tr09 = _tr04+c3m1*_tr09;	_ti09 = _ti04+c3m1*_ti09;\
	_rt   = s*_tr14;			_it   = s*_ti14;\
	__Br05 = _tr09-_it;			__Bi05 = _ti09+_rt;\
	__Br03 = _tr09+_it;			__Bi03 = _ti09-_rt;\
}

/* Totals: 162 ADD, 50 MUL	*/
//void ___RADIX_15_DIT() {}
#define RADIX_15_DIT(\
	__Ar00,__Ai00,__Ar01,__Ai01,__Ar02,__Ai02,__Ar03,__Ai03,__Ar04,__Ai04,\
	__Ar05,__Ai05,__Ar06,__Ai06,__Ar07,__Ai07,__Ar08,__Ai08,__Ar09,__Ai09,\
	__Ar10,__Ai10,__Ar11,__Ai11,__Ar12,__Ai12,__Ar13,__Ai13,__Ar14,__Ai14,\
	__Br00,__Bi00,__Br01,__Bi01,__Br02,__Bi02,__Br03,__Bi03,__Br04,__Bi04,\
	__Br05,__Bi05,__Br06,__Bi06,__Br07,__Bi07,__Br08,__Bi08,__Br09,__Bi09,\
	__Br10,__Bi10,__Br11,__Bi11,__Br12,__Bi12,__Br13,__Bi13,__Br14,__Bi14)\
{\
	double _tr00,_ti00,_tr01,_ti01,_tr02,_ti02,_tr03,_ti03,_tr04,_ti04,_tr05,_ti05,_tr06,_ti06,_tr07,_ti07,_tr08,_ti08,_tr09,_ti09,_tr10,_ti10,_tr11,_ti11,_tr12,_ti12,_tr13,_ti13,_tr14,_ti14;\
	double _rt,_it;\
/* Default input order for a twiddleless radix-15 DIT (cf.radix15_dit_pass1()) is [0,2,1,8,7,6,13,12,14,4,3,5,9,11,10]. */\
/*...gather the needed data (15 64-bit complex, i.e. 30 64-bit reals) and do 5 radix-3 transforms...*/\
/*...Block 1: Swap inputs 01 <-> 02: */\
	    _tr00 = __Ar00;				_ti00 = __Ai00;\
	    _tr01 = __Ar02;				_ti01 = __Ai02;\
	    _rt   = __Ar01;				_it   = __Ai01;\
	    _tr02 = _tr01-_rt;			_ti02 = _ti01-_it;\
	    _tr01 = _tr01+_rt;			_ti01 = _ti01+_it;\
	    _tr00 = _tr00+_tr01;		_ti00 = _ti00+_ti01;\
	    _tr01 = _tr00+c3m1*_tr01;	_ti01 = _ti00+c3m1*_ti01;\
	    _rt   = s*_tr02;			_it   = s*_ti02;\
	    _tr02 = _tr01-_it;			_ti02 = _ti01+_rt;\
	    _tr01 = _tr01+_it;			_ti01 = _ti01-_rt;\
\
/*...Block 2: Swap inputs 3,4,5 -> 8,7,6: */\
	    _tr03 = __Ar08;				_ti03 = __Ai08;\
	    _tr04 = __Ar07;				_ti04 = __Ai07;\
	    _rt   = __Ar06;				_it   = __Ai06;\
	    _tr05 = _tr04-_rt;			_ti05 = _ti04-_it;\
	    _tr04 = _tr04+_rt;			_ti04 = _ti04+_it;\
	    _tr03 = _tr03+_tr04;		_ti03 = _ti03+_ti04;\
	    _tr04 = _tr03+c3m1*_tr04;	_ti04 = _ti03+c3m1*_ti04;\
	    _rt   = s*_tr05;			_it   = s*_ti05;\
	    _tr05 = _tr04-_it;			_ti05 = _ti04+_rt;\
	    _tr04 = _tr04+_it;			_ti04 = _ti04-_rt;\
\
/*...Block 3: Swap inputs 6,7,8 -> 13,12,14: */\
	    _tr06 = __Ar13;				_ti06 = __Ai13;\
	    _tr07 = __Ar12;				_ti07 = __Ai12;\
	    _rt   = __Ar14;				_it   = __Ai14;\
	    _tr08 = _tr07-_rt;			_ti08 = _ti07-_it;\
	    _tr07 = _tr07+_rt;			_ti07 = _ti07+_it;\
	    _tr06 = _tr06+_tr07;		_ti06 = _ti06+_ti07;\
	    _tr07 = _tr06+c3m1*_tr07;	_ti07 = _ti06+c3m1*_ti07;\
	    _rt   = s*_tr08;			_it   = s*_ti08;\
	    _tr08 = _tr07-_it;			_ti08 = _ti07+_rt;\
	    _tr07 = _tr07+_it;			_ti07 = _ti07-_rt;\
\
/*...Block 4: Swap inputs 9,10,11 -> 4,3,5: */\
	    _tr09 = __Ar04;				_ti09 = __Ai04;\
	    _tr10 = __Ar03;				_ti10 = __Ai03;\
	    _rt   = __Ar05;				_it   = __Ai05;\
	    _tr11 = _tr10-_rt;			_ti11 = _ti10-_it;\
	    _tr10 = _tr10+_rt;			_ti10 = _ti10+_it;\
	    _tr09 = _tr09+_tr10;		_ti09 = _ti09+_ti10;\
	    _tr10 = _tr09+c3m1*_tr10;	_ti10 = _ti09+c3m1*_ti10;\
	    _rt   = s*_tr11;			_it   = s*_ti11;\
	    _tr11 = _tr10-_it;			_ti11 = _ti10+_rt;\
	    _tr10 = _tr10+_it;			_ti10 = _ti10-_rt;\
\
/*...Block 2: Swap inputs 12,13,14 -> 9,11,10: */\
	    _tr12 = __Ar09;				_ti12 = __Ai09;\
	    _tr13 = __Ar11;				_ti13 = __Ai11;\
	    _rt   = __Ar10;				_it   = __Ai10;\
	    _tr14 = _tr13-_rt;			_ti14 = _ti13-_it;\
	    _tr13 = _tr13+_rt;			_ti13 = _ti13+_it;\
	    _tr12 = _tr12+_tr13;		_ti12 = _ti12+_ti13;\
	    _tr13 = _tr12+c3m1*_tr13;	_ti13 = _ti12+c3m1*_ti13;\
	    _rt   = s*_tr14;			_it   = s*_ti14;\
	    _tr14 = _tr13-_it;			_ti14 = _ti13+_rt;\
	    _tr13 = _tr13+_it;			_ti13 = _ti13-_rt;\
\
/*...and now do three radix-5 transforms.\
	The required output permutation is [0,5,10,9,14,4,3,8,13,12,2,7,6,11,1]. */\
/*...Block 1: output permutation is 0,9,3,12,6 */\
	    _rt   = _tr12;				_it   = _ti12;\
	    _tr12 = _tr03-_rt;			_ti12 = _ti03-_it;\
	    _tr03 = _tr03+_rt;			_ti03 = _ti03+_it;\
	    _rt   = _tr09;				_it   = _ti09;\
	    _tr09 = _tr06-_rt;			_ti09 = _ti06-_it;\
	    _tr06 = _tr06+_rt;			_ti06 = _ti06+_it;\
\
	    _rt   = _tr03+_tr06;		_it   = _ti03+_ti06;\
	    _tr00 = _tr00+_rt;			_ti00 = _ti00+_it;\
	    _rt   = _tr00+cn1*_rt;		_it   = _ti00+cn1*_it;\
	    _tr06 = cn2*(_tr03-_tr06);	_ti06 = cn2*(_ti03-_ti06);\
	    _tr03 = _rt+_tr06;			_ti03 = _it+_ti06;\
	    _tr06 = _rt-_tr06;			_ti06 = _it-_ti06;\
	    _rt   = ss3*(_tr09-_tr12);	_it   = ss3*(_ti09-_ti12);\
	    _tr09 = _rt-sn1*_tr09;		_ti09 = _it-sn1*_ti09;\
	    _tr12 = _rt+sn2*_tr12;		_ti12 = _it+sn2*_ti12;\
\
	    __Br00 = _tr00;				__Bi00 = _ti00;\
	    __Br09 = _tr03-_ti09;		__Bi09 = _ti03+_tr09;\
	    __Br03 = _tr06-_ti12;		__Bi03 = _ti06+_tr12;\
	    __Br12 = _tr06+_ti12;		__Bi12 = _ti06-_tr12;\
	    __Br06 = _tr03+_ti09;		__Bi06 = _ti03-_tr09;\
\
/*...Block 2: output permutation is 5,14,8,2,11 */\
	    _rt   = _tr13;				_it   = _ti13;\
	    _tr13 = _tr04-_rt;			_ti13 = _ti04-_it;\
	    _tr04 = _tr04+_rt;			_ti04 = _ti04+_it;\
	    _rt   = _tr10;				_it   = _ti10;\
	    _tr10 = _tr07-_rt;			_ti10 = _ti07-_it;\
	    _tr07 = _tr07+_rt;			_ti07 = _ti07+_it;\
\
	    _rt   = _tr04+_tr07;		_it   = _ti04+_ti07;\
	    _tr01 = _tr01+_rt;			_ti01 = _ti01+_it;\
	    _rt   = _tr01+cn1*_rt;		_it   = _ti01+cn1*_it;\
	    _tr07 = cn2*(_tr04-_tr07);	_ti07 = cn2*(_ti04-_ti07);\
	    _tr04 = _rt+_tr07;			_ti04 = _it+_ti07;\
	    _tr07 = _rt-_tr07;			_ti07 = _it-_ti07;\
	    _rt   = ss3*(_tr10-_tr13);	_it   = ss3*(_ti10-_ti13);\
	    _tr10 = _rt-sn1*_tr10;		_ti10 = _it-sn1*_ti10;\
	    _tr13 = _rt+sn2*_tr13;		_ti13 = _it+sn2*_ti13;\
\
	    __Br05 = _tr01;				__Bi05 = _ti01;\
	    __Br14 = _tr04-_ti10;		__Bi14 = _ti04+_tr10;\
	    __Br08 = _tr07-_ti13;		__Bi08 = _ti07+_tr13;\
	    __Br02 = _tr07+_ti13;		__Bi02 = _ti07-_tr13;\
	    __Br11 = _tr04+_ti10;		__Bi11 = _ti04-_tr10;\
\
/*...Block 3: output permutation is 10,4,13,7,1 */\
	    _rt   = _tr14;				_it   = _ti14;\
	    _tr14 = _tr05-_rt;			_ti14 = _ti05-_it;\
	    _tr05 = _tr05+_rt;			_ti05 = _ti05+_it;\
	    _rt   = _tr11;				_it   = _ti11;\
	    _tr11 = _tr08-_rt;			_ti11 = _ti08-_it;\
	    _tr08 = _tr08+_rt;			_ti08 = _ti08+_it;\
\
	    _rt   = _tr05+_tr08;		_it   = _ti05+_ti08;\
	    _tr02 = _tr02+_rt;			_ti02 = _ti02+_it;\
	    _rt   = _tr02+cn1*_rt;		_it   = _ti02+cn1*_it;\
	    _tr08 = cn2*(_tr05-_tr08);	_ti08 = cn2*(_ti05-_ti08);\
	    _tr05 = _rt+_tr08;			_ti05 = _it+_ti08;\
	    _tr08 = _rt-_tr08;			_ti08 = _it-_ti08;\
	    _rt   = ss3*(_tr11-_tr14);	_it   = ss3*(_ti11-_ti14);\
	    _tr11 = _rt-sn1*_tr11;		_ti11 = _it-sn1*_ti11;\
	    _tr14 = _rt+sn2*_tr14;		_ti14 = _it+sn2*_ti14;\
\
	    __Br10 = _tr02;				__Bi10 = _ti02;\
	    __Br04 = _tr05-_ti11;		__Bi04 = _ti05+_tr11;\
	    __Br13 = _tr08-_ti14;		__Bi13 = _ti08+_tr14;\
	    __Br07 = _tr08+_ti14;		__Bi07 = _ti08-_tr14;\
	    __Br01 = _tr05+_ti11;		__Bi01 = _ti05-_tr11;\
}

/* Scalar-data test macros for comparing with the SSE2 radix-15 DFTs: */
//void ___RADIX_15_DIF_B() {}
#define RADIX_15_DIF_B(\
	__s,__c3m1,\
	__cn1, __cn2, __ss3, __sn1, __sn2,\
	__Ar00,__Ai00,__Ar01,__Ai01,__Ar02,__Ai02,__Ar03,__Ai03,__Ar04,__Ai04,__Ar05,__Ai05,__Ar06,__Ai06,__Ar07,__Ai07,__Ar08,__Ai08,__Ar09,__Ai09,__Ar10,__Ai10,__Ar11,__Ai11,__Ar12,__Ai12,__Ar13,__Ai13,__Ar14,__Ai14,\
	__Br00,__Bi00,__Br01,__Bi01,__Br02,__Bi02,__Br03,__Bi03,__Br04,__Bi04,__Br05,__Bi05,__Br06,__Bi06,__Br07,__Bi07,__Br08,__Bi08,__Br09,__Bi09,__Br10,__Bi10,__Br11,__Bi11,__Br12,__Bi12,__Br13,__Bi13,__Br14,__Bi14\
)\
{\
	double _tr00,_ti00,_tr01,_ti01,_tr02,_ti02,_tr03,_ti03,_tr04,_ti04,_tr05,_ti05,_tr06,_ti06,_tr07,_ti07,_tr08,_ti08,_tr09,_ti09,_tr10,_ti10,_tr11,_ti11,_tr12,_ti12,_tr13,_ti13,_tr14,_ti14;\
	double _x0,_x1,_x2,_x3,_x4,_x5, _rt,_it;\
\
	RADIX_05_DFT(__cn1,__cn2,__ss3,__sn1,__sn2,__Ar00,__Ai00,__Ar12,__Ai12,__Ar09,__Ai09,__Ar06,__Ai06,__Ar03,__Ai03, _tr00,_ti00,_tr01,_ti01,_tr02,_ti02,_tr03,_ti03,_tr04,_ti04, _rt,_it);\
	RADIX_05_DFT(__cn1,__cn2,__ss3,__sn1,__sn2,__Ar10,__Ai10,__Ar07,__Ai07,__Ar04,__Ai04,__Ar01,__Ai01,__Ar13,__Ai13, _tr05,_ti05,_tr06,_ti06,_tr07,_ti07,_tr08,_ti08,_tr09,_ti09, _rt,_it);\
	RADIX_05_DFT(__cn1,__cn2,__ss3,__sn1,__sn2,__Ar05,__Ai05,__Ar02,__Ai02,__Ar14,__Ai14,__Ar11,__Ai11,__Ar08,__Ai08, _tr10,_ti10,_tr11,_ti11,_tr12,_ti12,_tr13,_ti13,_tr14,_ti14, _rt,_it);\
\
	RADIX_03_DFT(__s,__c3m1,_tr00,_ti00,_tr05,_ti05,_tr10,_ti10, _x0,_x1,_x2,_x3,_x4,_x5, __Br00,__Bi00,__Br01,__Bi01,__Br02,__Bi02);\
	RADIX_03_DFT(__s,__c3m1,_tr01,_ti01,_tr06,_ti06,_tr11,_ti11, _x0,_x1,_x2,_x3,_x4,_x5, __Br13,__Bi13,__Br14,__Bi14,__Br12,__Bi12);\
	RADIX_03_DFT(__s,__c3m1,_tr02,_ti02,_tr07,_ti07,_tr12,_ti12, _x0,_x1,_x2,_x3,_x4,_x5, __Br09,__Bi09,__Br10,__Bi10,__Br11,__Bi11);\
	RADIX_03_DFT(__s,__c3m1,_tr03,_ti03,_tr08,_ti08,_tr13,_ti13, _x0,_x1,_x2,_x3,_x4,_x5, __Br08,__Bi08,__Br06,__Bi06,__Br07,__Bi07);\
	RADIX_03_DFT(__s,__c3m1,_tr04,_ti04,_tr09,_ti09,_tr14,_ti14, _x0,_x1,_x2,_x3,_x4,_x5, __Br04,__Bi04,__Br05,__Bi05,__Br03,__Bi03);\
}

//void ___RADIX_15_DIT_B() {}
#define RADIX_15_DIT_B(\
	__s,__c3m1,\
	__cn1, __cn2, __ss3, __sn1, __sn2,\
	__Ar00,__Ai00,__Ar01,__Ai01,__Ar02,__Ai02,__Ar03,__Ai03,__Ar04,__Ai04,__Ar05,__Ai05,__Ar06,__Ai06,__Ar07,__Ai07,__Ar08,__Ai08,__Ar09,__Ai09,__Ar10,__Ai10,__Ar11,__Ai11,__Ar12,__Ai12,__Ar13,__Ai13,__Ar14,__Ai14,\
	__Br00,__Bi00,__Br01,__Bi01,__Br02,__Bi02,__Br03,__Bi03,__Br04,__Bi04,__Br05,__Bi05,__Br06,__Bi06,__Br07,__Bi07,__Br08,__Bi08,__Br09,__Bi09,__Br10,__Bi10,__Br11,__Bi11,__Br12,__Bi12,__Br13,__Bi13,__Br14,__Bi14\
)\
{\
	double _tr00,_ti00,_tr01,_ti01,_tr02,_ti02,_tr03,_ti03,_tr04,_ti04,_tr05,_ti05,_tr06,_ti06,_tr07,_ti07,_tr08,_ti08,_tr09,_ti09,_tr10,_ti10,_tr11,_ti11,_tr12,_ti12,_tr13,_ti13,_tr14,_ti14;\
	double _x0,_x1,_x2,_x3,_x4,_x5, _rt,_it;\
	/* Swap the 2nd pair of each output triplet to effect iDFT: */\
	RADIX_03_DFT(__s,__c3m1,__Ar00,__Ai00,__Ar02,__Ai02,__Ar01,__Ai01,_x0, _x1,_x2,_x3,_x4,_x5, _tr00,_ti00,_tr02,_ti02,_tr01,_ti01);\
	RADIX_03_DFT(__s,__c3m1,__Ar08,__Ai08,__Ar07,__Ai07,__Ar06,__Ai06,_x0, _x1,_x2,_x3,_x4,_x5, _tr03,_ti03,_tr05,_ti05,_tr04,_ti04);\
	RADIX_03_DFT(__s,__c3m1,__Ar13,__Ai13,__Ar12,__Ai12,__Ar14,__Ai14,_x0, _x1,_x2,_x3,_x4,_x5, _tr06,_ti06,_tr08,_ti08,_tr07,_ti07);\
	RADIX_03_DFT(__s,__c3m1,__Ar04,__Ai04,__Ar03,__Ai03,__Ar05,__Ai05,_x0, _x1,_x2,_x3,_x4,_x5, _tr09,_ti09,_tr11,_ti11,_tr10,_ti10);\
	RADIX_03_DFT(__s,__c3m1,__Ar09,__Ai09,__Ar11,__Ai11,__Ar10,__Ai10,_x0, _x1,_x2,_x3,_x4,_x5, _tr12,_ti12,_tr14,_ti14,_tr13,_ti13);\
	/* Output perm here is 0123456789abcde --> 05a6b1c2738d9e4: */\
	RADIX_05_DFT(__cn1,__cn2,__ss3,__sn1,__sn2, _tr00,_ti00,_tr03,_ti03,_tr06,_ti06,_tr09,_ti09,_tr12,_ti12, __Br00,__Bi00,__Br06,__Bi06,__Br12,__Bi12,__Br03,__Bi03,__Br09,__Bi09,_rt,_it);\
	RADIX_05_DFT(__cn1,__cn2,__ss3,__sn1,__sn2, _tr01,_ti01,_tr04,_ti04,_tr07,_ti07,_tr10,_ti10,_tr13,_ti13, __Br05,__Bi05,__Br11,__Bi11,__Br02,__Bi02,__Br08,__Bi08,__Br14,__Bi14,_rt,_it);\
	RADIX_05_DFT(__cn1,__cn2,__ss3,__sn1,__sn2, _tr02,_ti02,_tr05,_ti05,_tr08,_ti08,_tr11,_ti11,_tr14,_ti14, __Br10,__Bi10,__Br01,__Bi01,__Br07,__Bi07,__Br13,__Bi13,__Br04,__Bi04,_rt,_it);\
}

/* Twiddleless radix-16 DFTs: */
/* Totals: 144 ADD, 24 MUL. With twiddles adds 30 CMUL @[4 MUL, 2 ADD] each and thus needs 174 ADD, 84 MUL. */
//void ___RADIX_16_DIF() {}
#define RADIX_16_DIF(\
	__A0r,__A0i,__A1r,__A1i,__A2r,__A2i,__A3r,__A3i,__A4r,__A4i,__A5r,__A5i,__A6r,__A6i,__A7r,__A7i,__A8r,__A8i,__A9r,__A9i,__AAr,__AAi,__ABr,__ABi,__ACr,__ACi,__ADr,__ADi,__AEr,__AEi,__AFr,__AFi,\
	__B0r,__B0i,__B1r,__B1i,__B2r,__B2i,__B3r,__B3i,__B4r,__B4i,__B5r,__B5i,__B6r,__B6i,__B7r,__B7i,__B8r,__B8i,__B9r,__B9i,__BAr,__BAi,__BBr,__BBi,__BCr,__BCi,__BDr,__BDi,__BEr,__BEi,__BFr,__BFi,\
	__c,__s)\
{\
	double _t1,_t2,_t3,_t4,_t5,_t6,_t7,_t8,_t9,_t10,_t11,_t12,_t13,_t14,_t15,_t16,_t17,_t18,_t19,_t20,_t21,_t22,_t23,_t24,_t25,_t26,_t27,_t28,_t29,_t30,_t31,_t32;\
	double _r,_i;\
	/*       gather the needed data (16 64-bit complex, i.e. 32 64-bit reals) and do the first set of four length-4 transforms.	*/\
	/*...Block 1: */\
	_t3 =__A0r -__A8r;	_t1 =__A0r +__A8r;\
	_t4 =__A0i -__A8i;	_t2 =__A0i +__A8i;\
\
	_t7 =__A4r -__ACr;	_t5 =__A4r +__ACr;\
	_t8 =__A4i -__ACi;	_t6 =__A4i +__ACi;\
\
	_r =_t5;	_t5 =_t1 -_r;	_t1 =_t1 +_r;\
	_i =_t6;	_t6 =_t2 -_i;	_t2 =_t2 +_i;\
\
	_r =_t7;	_t7 =_t3 +_t8;	_t3 =_t3 -_t8;\
				_t8 =_t4 -_r;	_t4 =_t4 +_r;\
\
	/*...Block 2: */\
	_t11=__A2r -__AAr;	_t9 =__A2r +__AAr;\
	_t12=__A2i -__AAi;	_t10=__A2i +__AAi;\
\
	_t15=__A6r -__AEr;	_t13=__A6r +__AEr;\
	_t16=__A6i -__AEi;	_t14=__A6i +__AEi;\
\
	_r =_t13;	_t13=_t9 -_r;		_t9 =_t9 +_r;\
	_i =_t14;	_t14=_t10-_i;		_t10=_t10+_i;\
\
	_r =_t15;	_t15=_t11+_t16;	_t11=_t11-_t16;\
				_t16=_t12-_r;		_t12=_t12+_r;\
\
	/*...Block 3: */\
	_t19=__A1r -__A9r;	_t17=__A1r +__A9r;\
	_t20=__A1i -__A9i;	_t18=__A1i +__A9i;\
\
	_t23=__A5r -__ADr;	_t21=__A5r +__ADr;\
	_t24=__A5i -__ADi;	_t22=__A5i +__ADi;\
\
	_r =_t21;	_t21=_t17-_r;		_t17=_t17+_r;\
	_i =_t22;	_t22=_t18-_i;		_t18=_t18+_i;\
\
	_r =_t23;	_t23=_t19+_t24;	_t19=_t19-_t24;\
				_t24=_t20-_r;		_t20=_t20+_r;\
\
	/*...Block 4: */\
	_t27=__A3r -__ABr;	_t25=__A3r +__ABr;\
	_t28=__A3i -__ABi;	_t26=__A3i +__ABi;\
\
	_t31=__A7r -__AFr;	_t29=__A7r +__AFr;\
	_t32=__A7i -__AFi;	_t30=__A7i +__AFi;\
\
	_r =_t29;	_t29=_t25-_r;		_t25=_t25+_r;\
	_i =_t30;	_t30=_t26-_i;		_t26=_t26+_i;\
\
	_r =_t31;	_t31=_t27+_t32;	_t27=_t27-_t32;\
				_t32=_t28-_r;		_t28=_t28+_r;\
\
/**************************************************************************************/\
/*...and now do four more radix-4 transforms, including the internal twiddle factors: */\
/**************************************************************************************/\
\
	/*...Block 1: t1,9,17,25	*/\
	_r =_t9 ;	_t9 =_t1 -_r;	_t1 =_t1 +_r;\
	_i =_t10;	_t10=_t2 -_i;	_t2 =_t2 +_i;\
\
	_r =_t25;	_t25=_t17-_r;	_t17=_t17+_r;\
	_i =_t26;	_t26=_t18-_i;	_t18=_t18+_i;\
\
	__B0r =_t1+_t17;			__B0i =_t2+_t18;\
	__B1r =_t1-_t17;			__B1i =_t2-_t18;\
\
	__B2r =_t9 -_t26;			__B2i =_t10+_t25;	/* mpy by E^4=i is inlined here...	*/\
	__B3r =_t9 +_t26;			__B3i =_t10-_t25;\
\
	/*...Block 3: _t5,13,21,29	*/\
	_r =_t13;	_t13=_t5 +_t14;	_t5 =_t5 -_t14;		/* Twiddle mpy by E^4 = I	*/\
				_t14=_t6 -_r;	_t6 =_t6 +_r;\
\
	_r =(_t21-_t22)*ISRT2;		_t22=(_t21+_t22)*ISRT2;	_t21=_r;	/* Twiddle mpy by E^2; could also do mpy by ISRT2 directly on t21,22,29,30	*/\
	_r =(_t30+_t29)*ISRT2;		_i =(_t30-_t29)*ISRT2;	/* Twiddle mpy by -E^6 is here...	*/\
	_t29=_t21+_r;				_t21=_t21-_r;				/* ...and get E^6=(i-1)/sq_r by flipping signs here.	*/\
	_t30=_t22+_i;				_t22=_t22-_i;\
\
	__B4r =_t5+_t21;			__B4i =_t6+_t22;\
	__B5r =_t5-_t21;			__B5i =_t6-_t22;\
\
	__B6r =_t13-_t30;			__B6i =_t14+_t29;			/* mpy by E^4=i is inlined here...	*/\
	__B7r =_t13+_t30;			__B7i =_t14-_t29;\
\
	/*...Block 2: _t3,11,19,27	*/\
	_r =(_t11-_t12)*ISRT2;		_i =(_t11+_t12)*ISRT2;	/* Twiddle mpy by E^2; could also do mpy by ISRT2 directly on t11,12	*/\
	_t11=_t3 -_r;				_t3 =_t3 +_r;\
	_t12=_t4 -_i;				_t4 =_t4 +_i;\
\
	_r =_t19*__c - _t20*__s;	_t20=_t20*__c + _t19*__s;	_t19=_r;	/* Twiddle mpy by E^1	*/\
	_r =_t27*__s - _t28*__c;	_i =_t28*__s + _t27*__c;	/* Twiddle mpy by E^3	*/\
	_t27=_t19-_r;				_t19=_t19+_r;\
	_t28=_t20-_i;				_t20=_t20+_i;\
\
	__B8r =_t3+_t19;			__B8i =_t4+_t20;\
	__B9r =_t3-_t19;			__B9i =_t4-_t20;\
\
	__BAr =_t11-_t28;			__BAi =_t12+_t27;			/* mpy by E^4=i is inlined here...	*/\
	__BBr =_t11+_t28;			__BBi =_t12-_t27;\
\
	/*...Block 4: _t7,15,23,31	*/\
	_r =(_t16+_t15)*ISRT2;		_i =(_t16-_t15)*ISRT2;	/* Twiddle mpy by -E^6 is here...	*/\
	_t15=_t7 +_r;				_t7 =_t7 -_r;			/* ...and get E^6=(i-1)/sqrt2 by flipping signs here.	*/\
	_t16=_t8 +_i;				_t8 =_t8 -_i;\
\
	_r =_t23*__s - _t24*__c;	_t24=_t24*__s + _t23*__c;	_t23=_r;	/* Twiddle mpy by E^3	*/\
	_r =_t31*__c - _t32*__s;	_i =_t32*__c + _t31*__s;	/* Twiddle mpy by E^1 = -E^9...	*/\
	_t31=_t23+_r;				_t23=_t23-_r;			/* ...and get E^9 by flipping signs here.	*/\
	_t32=_t24+_i;				_t24=_t24-_i;			/* Note: t23+_r = t23*(s+1)	*/\
\
	__BCr =_t7+_t23;			__BCi =_t8+_t24;\
	__BDr =_t7-_t23;			__BDi =_t8-_t24;\
\
	__BEr =_t15-_t32;			__BEi =_t16+_t31;			/* mpy by E^4=i is inlined here...	*/\
	__BFr =_t15+_t32;			__BFi =_t16-_t31;\
}

/* Totals: 144 ADD, 24 MUL	*/
//void ___RADIX_16_DIT() {}
#define RADIX_16_DIT(\
	__A0r,__A0i,__A1r,__A1i,__A2r,__A2i,__A3r,__A3i,__A4r,__A4i,__A5r,__A5i,__A6r,__A6i,__A7r,__A7i,__A8r,__A8i,__A9r,__A9i,__AAr,__AAi,__ABr,__ABi,__ACr,__ACi,__ADr,__ADi,__AEr,__AEi,__AFr,__AFi,\
	__B0r,__B0i,__B1r,__B1i,__B2r,__B2i,__B3r,__B3i,__B4r,__B4i,__B5r,__B5i,__B6r,__B6i,__B7r,__B7i,__B8r,__B8i,__B9r,__B9i,__BAr,__BAi,__BBr,__BBi,__BCr,__BCi,__BDr,__BDi,__BEr,__BEi,__BFr,__BFi,\
	__c,__s)\
{\
	double _t1,_t2,_t3,_t4,_t5,_t6,_t7,_t8,_t9,_t10,_t11,_t12,_t13,_t14,_t15,_t16,_t17,_t18,_t19,_t20,_t21,_t22,_t23,_t24,_t25,_t26,_t27,_t28,_t29,_t30,_t31,_t32;\
	double _r,_i;\
	/*...Block 1: */\
	_t3 =__A0r -__A1r;	_t1 =__A0r +__A1r;\
	_t4 =__A0i -__A1i;	_t2 =__A0i +__A1i;\
\
	_t7 =__A2r -__A3r;	_t5 =__A2r +__A3r;\
	_t8 =__A2i -__A3i;	_t6 =__A2i +__A3i;\
\
	_r =_t5;		_t5 =_t1 -_r;	_t1 =_t1 +_r;\
	_i =_t6;		_t6 =_t2 -_i;	_t2 =_t2 +_i;\
\
	_r =_t7;		_t7 =_t3 -_t8;	_t3 =_t3 +_t8;\
					_t8 =_t4 +_r;	_t4 =_t4 -_r;\
\
	/*...Block 2: */\
	_t11=__A4r -__A5r;	_t9 =__A4r +__A5r;\
	_t12=__A4i -__A5i;	_t10=__A4i +__A5i;\
\
	_t15=__A6r -__A7r;	_t13=__A6r +__A7r;\
	_t16=__A6i -__A7i;	_t14=__A6i +__A7i;\
\
	_r =_t13;	_t13=_t9 -_r;		_t9 =_t9 +_r;\
	_i =_t14;	_t14=_t10-_i;		_t10=_t10+_i;\
\
	_r =_t15;	_t15=_t11-_t16;	_t11=_t11+_t16;\
				_t16=_t12+_r;		_t12=_t12-_r;\
\
	/*...Block 3: */\
	_t19=__A8r -__A9r;	_t17=__A8r +__A9r;\
	_t20=__A8i -__A9i;	_t18=__A8i +__A9i;\
\
	_t23=__AAr -__ABr;	_t21=__AAr +__ABr;\
	_t24=__AAi -__ABi;	_t22=__AAi +__ABi;\
\
	_r =_t21;	_t21=_t17-_r;		_t17=_t17+_r;\
	_i =_t22;	_t22=_t18-_i;		_t18=_t18+_i;\
\
	_r =_t23;	_t23=_t19-_t24;	_t19=_t19+_t24;\
				_t24=_t20+_r;		_t20=_t20-_r;\
\
	/*...Block 4: */\
	_t27=__ACr -__ADr;	_t25=__ACr +__ADr; \
	_t28=__ACi -__ADi;	_t26=__ACi +__ADi;\
\
	_t31=__AEr -__AFr;	_t29=__AEr +__AFr;\
	_t32=__AEi -__AFi;	_t30=__AEi +__AFi;\
\
	_r =_t29;	_t29=_t25-_r;		_t25=_t25+_r;\
	_i =_t30;	_t30=_t26-_i;		_t26=_t26+_i;\
\
	_r =_t31;	_t31=_t27-_t32;	_t27=_t27+_t32;\
				_t32=_t28+_r;		_t28=_t28-_r;\
\
/**************************************************************************************/\
/*...and now do four more radix-4 transforms, including the internal twiddle factors: */\
/**************************************************************************************/\
\
	/*...Block 1: _t1,9,17,25	*/\
	_r =_t9 ;	_t9 =_t1 -_r;	_t1 =_t1 +_r;\
	_i =_t10;	_t10=_t2 -_i;	_t2 =_t2 +_i;\
\
	_r =_t25;	_t25=_t17-_r;	_t17=_t17+_r;\
	_i =_t26;	_t26=_t18-_i;	_t18=_t18+_i;\
\
	__B0r=_t1+_t17;				__B0i =_t2+_t18;\
	__B8r=_t1-_t17;				__B8i =_t2-_t18;\
\
	__B4r=_t9 +_t26;			__B4i =_t10-_t25;	/* mpy by E^-4 = -I is inlined here...	*/\
	__BCr=_t9 -_t26;			__BCi =_t10+_t25;\
\
	/*...Block 3: _t5,13,21,29	*/\
	_r =_t13;_t13=_t5 -_t14;	_t5 =_t5 +_t14;	/* Twiddle mpy by E^4 =-I	*/\
			_t14=_t6 +_r;		_t6 =_t6 -_r;\
\
	_r =(_t22+_t21)*ISRT2;		_t22=(_t22-_t21)*ISRT2;	_t21=_r;	/* Twiddle mpy by E^-2	*/\
	_r =(_t29-_t30)*ISRT2;		_i =(_t29+_t30)*ISRT2;	/* Twiddle mpy by E^2 = -E^-6 is here...	*/\
	_t29=_t21+_r;				_t21=_t21-_r;				/* ...and get E^-6 by flipping signs here.	*/\
	_t30=_t22+_i;				_t22=_t22-_i;\
\
	__B2r=_t5+_t21;				__B2i =_t6+_t22;\
	__BAr=_t5-_t21;				__BAi =_t6-_t22;\
\
	__B6r=_t13+_t30;			__B6i =_t14-_t29;			/* mpy by E^-4 =-I is inlined here...	*/\
	__BEr=_t13-_t30;			__BEi =_t14+_t29;\
\
	/*...Block 2: _t3,11,19,27	*/\
	_r =(_t12+_t11)*ISRT2;		_i =(_t12-_t11)*ISRT2;	/* Twiddle mpy by E^-2	*/\
	_t11=_t3 -_r;				_t3 =_t3 +_r;\
	_t12=_t4 -_i;				_t4 =_t4 +_i;\
\
	_r =_t19*__c + _t20*__s;	_t20=_t20*__c - _t19*__s;	_t19=_r;	/* Twiddle mpy by E^-1	*/\
	_r =_t27*__s + _t28*__c;	_i =_t28*__s - _t27*__c;	/* Twiddle mpy by E^-3	*/\
	_t27=_t19-_r;				_t19=_t19+_r;\
	_t28=_t20-_i;				_t20=_t20+_i;\
\
	__B1r=_t3+_t19;				__B1i =_t4+_t20;\
	__B9r=_t3-_t19;				__B9i =_t4-_t20;\
\
	__B5r=_t11+_t28;			__B5i =_t12-_t27;			/* mpy by E^-4 =-I is inlined here...	*/\
	__BDr=_t11-_t28;			__BDi =_t12+_t27;\
\
	/*...Block 4: _t7,15,23,31	*/\
	_r =(_t15-_t16)*ISRT2;		_i =(_t15+_t16)*ISRT2;	/* Twiddle mpy by E^2 = -E^-6 is here...	*/\
	_t15=_t7 +_r;				_t7 =_t7 -_r;			/* ...and get E^6=(i-1)/sqrt2 by flipping signs here.	*/\
	_t16=_t8 +_i;				_t8 =_t8 -_i;\
\
	_r =_t23*__s + _t24*__c;	_t24=_t24*__s - _t23*__c;	_t23=_r;	/* Twiddle mpy by E^-3	*/\
	_r =_t31*__c + _t32*__s;	_i =_t32*__c - _t31*__s;	/* Twiddle mpy by E^-1 = -E^-9...	*/\
	_t31=_t23+_r;				_t23=_t23-_r;			/* ...and get E^9 by flipping signs here.	*/\
	_t32=_t24+_i;				_t24=_t24-_i;\
\
	__B3r=_t7+_t23;				__B3i =_t8+_t24;\
	__BBr=_t7-_t23;				__BBi =_t8-_t24;\
\
	__B7r=_t15+_t32;			__B7i =_t16-_t31;			/* mpy by E^-4 = -I is inlined here...	*/\
	__BFr=_t15-_t32;			__BFi =_t16+_t31;\
}

/* With-Twiddle in=place radix-16 DFTs: */
/* Totals: 174 ADD, 84 MUL, which is the cost of the twiddle-less version plus 15 CMUL at [2 ADD, 4 MUL] each. */
//void ___RADIX_16_DIF_TWIDDLE() {}
#define RADIX_16_DIF_TWIDDLE(\
	__A0r,__A0i,__A1r,__A1i,__A2r,__A2i,__A3r,__A3i,__A4r,__A4i,__A5r,__A5i,__A6r,__A6i,__A7r,__A7i,__A8r,__A8i,__A9r,__A9i,__AAr,__AAi,__ABr,__ABi,__ACr,__ACi,__ADr,__ADi,__AEr,__AEi,__AFr,__AFi,\
	__c1 ,__s1 ,__c2 ,__s2 ,__c3 ,__s3 ,__c4 ,__s4 ,__c5 ,__s5 ,__c6 ,__s6 ,__c7 ,__s7 ,__c8 ,__s8 ,__c9 ,__s9 ,__cA ,__sA ,__cB ,__sB ,__cC ,__sC ,__cD ,__sD ,__cE ,__sE ,__cF ,__sF ,\
	__c,__s)\
{\
	double _t1,_t2,_t3,_t4,_t5,_t6,_t7,_t8,_t9,_t10,_t11,_t12,_t13,_t14,_t15,_t16,_t17,_t18,_t19,_t20,_t21,_t22,_t23,_t24,_t25,_t26,_t27,_t28,_t29,_t30,_t31,_t32;\
	double _r,_i;\
	/* Gather the needed data (16 64-bit complex, i.e. 32 64-bit reals) and do first set of four length-4 transforms: */\
/* Pass 1 - each sub-4-DFT costs [24 ADD, 16 MUL], except for the first, which has just 3 non-unity twiddles
and as a result costs 2 fewer ADD and 4 fewer MUL. Total Pass 1 cost: [94 ADD, 60 MUL]. Pass 2 cost: [80 ADD, 24 MUL]
*/\
	/*...Block 1: */\
	_t1 =__A0r;						_t2 =__A0i;\
	_r  =__A8r*__c8 -__A8i*__s8 ;	_i  =__A8i*__c8 +__A8r*__s8 ;\
	_t3 =_t1 -_r;					_t1 =_t1 +_r;	/* t 1, 2: w0.x0 + w8.x8 */\
	_t4 =_t2 -_i;					_t2 =_t2 +_i;	/* t 3, 4: w0.x0 - w8.x8 */\
\
	_t5 =__A4r*__c4 -__A4i*__s4 ;	_t6 =__A4i*__c4 +__A4r*__s4 ;\
	_r  =__ACr*__cC -__ACi*__sC ;	_i  =__ACi*__cC +__ACr*__sC ;\
	_t7 =_t5 -_r;					_t5 =_t5 +_r;	/* t 5, 6: w4.x4 + wc.xc */\
	_t8 =_t6 -_i;					_t6 =_t6 +_i;	/* t 7, 8: w4.x4 - wc.xc */\
\
	_r =_t5;	_t5 =_t1 -_r;		_t1 =_t1 +_r;	/* y0 = t 1, 2: [w0.x0 + w8.x8] +   [w4.x4 + wc.xc] */\
	_i =_t6;	_t6 =_t2 -_i;		_t2 =_t2 +_i;	/* y2 = t 5, 6: [w0.x0 + w8.x8] -   [w4.x4 + wc.xc] */\
\
	_r =_t7;	_t7 =_t3 +_t8;		_t3 =_t3 -_t8;	/* y1 = t 3, 4: [w0.x0 - w8.x8] + I.[w4.x4 - wc.xc] */\
				_t8 =_t4 -_r;		_t4 =_t4 +_r;	/* y3 = t 7, 8: [w0.x0 - w8.x8] - I.[w4.x4 - wc.xc] */\
\
	/*...Block 2: */\
	_t9 =__A2r*__c2 -__A2i*__s2 ;	_t10=__A2i*__c2 +__A2r*__s2;\
	_r  =__AAr*__cA -__AAi*__sA ;	_i  =__AAi*__cA +__AAr*__sA ;\
	_t11=_t9 -_r;					_t9 =_t9 +_r;\
	_t12=_t10-_i;					_t10=_t10+_i;	/* w2.x2 +- wa.xa */\
\
	_t13=__A6r*__c6 -__A6i*__s6 ;	_t14=__A6i*__c6 +__A6r*__s6;\
	_r  =__AEr*__cE -__AEi*__sE ;	_i  =__AEi*__cE +__AEr*__sE ;\
	_t15=_t13-_r;					_t13=_t13+_r;\
	_t16=_t14-_i;					_t14=_t14+_i;	/* w6.x6 +- we.xe */\
\
	_r =_t13;	_t13=_t9 -_r;		_t9 =_t9 +_r;\
	_i =_t14;	_t14=_t10-_i;		_t10=_t10+_i;\
\
	_r =_t15;	_t15=_t11+_t16;		_t11=_t11-_t16;\
				_t16=_t12-_r;		_t12=_t12+_r;\
\
	/*...Block 3: */\
	_t17=__A1r*__c1 -__A1i*__s1 ;	_t18=__A1i*__c1 +__A1r*__s1;\
	_r  =__A9r*__c9 -__A9i*__s9 ;	_i  =__A9i*__c9 +__A9r*__s9;\
	_t19=_t17-_r;					_t17=_t17+_r;\
	_t20=_t18-_i;					_t18=_t18+_i;	/* w1.x1 +- w9.x9 */\
\
	_t21=__A5r*__c5 -__A5i*__s5 ;	_t22=__A5i*__c5 +__A5r*__s5;\
	_r  =__ADr*__cD -__ADi*__sD ;	_i  =__ADi*__cD +__ADr*__sD ;\
	_t23=_t21-_r;					_t21=_t21+_r;\
	_t24=_t22-_i;					_t22=_t22+_i;	/* w5.x5 +- wd.xd */\
\
	_r =_t21;	_t21=_t17-_r;		_t17=_t17+_r;\
	_i =_t22;	_t22=_t18-_i;		_t18=_t18+_i;\
\
	_r =_t23;	_t23=_t19+_t24;		_t19=_t19-_t24;\
				_t24=_t20-_r;		_t20=_t20+_r;\
\
	/*...Block 4: */\
	_t25=__A3r*__c3 -__A3i*__s3 ;	_t26=__A3i*__c3 +__A3r*__s3;\
	_r  =__ABr*__cB -__ABi*__sB ;	_i  =__ABi*__cB +__ABr*__sB ;\
	_t27=_t25-_r;					_t25=_t25+_r;\
	_t28=_t26-_i;					_t26=_t26+_i;	/* w3.x3 +- wb.xb */\
\
	_t29=__A7r*__c7 -__A7i*__s7 ;	_t30=__A7i*__c7 +__A7r*__s7;\
	_r  =__AFr*__cF -__AFi*__sF ;	_i  =__AFi*__cF +__AFr*__sF ;\
	_t31=_t29-_r;					_t29=_t29+_r;\
	_t32=_t30-_i;					_t30=_t30+_i;	/* w7.x7 +- wf.xf */\
\
	_r =_t29;	_t29=_t25-_r;		_t25=_t25+_r;\
	_i =_t30;	_t30=_t26-_i;		_t26=_t26+_i;\
\
	_r =_t31;	_t31=_t27+_t32;		_t27=_t27-_t32;\
				_t32=_t28-_r;		_t28=_t28+_r;\
\
/************************************************************************************/\
/* Pass 2: Do four more radix-4 transforms, including the internal twiddle factors: */\
/************************************************************************************/\
\
/* Pass 1 - each sub-4-DFT costs [24 ADD, 16 MUL], except for the first, which has just 3 non-unity twiddles
and as a result costs 2 fewer ADD and 4 fewer MUL. Total Pass 1 cost: [94 ADD, 60 MUL].
*/\
	/*...Block 1: Inputs y0,4,8,c: */\
	_r =_t9;	_t9 =_t1 -_r;	_t1 =_t1 +_r;	/* t 1, 2: y0 + y4 */\
	_i =_t10;	_t10=_t2 -_i;	_t2 =_t2 +_i;	/* t 9,10: y0 - y4 */\
\
	_r =_t25;	_t25=_t17-_r;	_t17=_t17+_r;	/* t17,18: y8 + yc */\
	_i =_t26;	_t26=_t18-_i;	_t18=_t18+_i;	/* t25,26: y8 - yc */\
\
	__A0r=_t1 +_t17;	__A0i=_t2 +_t18;	/* A0: [y0 + y4] +   [y8 + yc] */\
	__A1r=_t1 -_t17;	__A1i=_t2 -_t18;	/* A1: [y0 + y4] -   [y8 + yc] */\
\
	__A2r=_t9 -_t26;	__A2i=_t10+_t25;	/* A2: [y0 - y4] + I.[y8 - yc] */\
	__A3r=_t9 +_t26;	__A3i=_t10-_t25;	/* A3: [y0 - y4] - I.[y8 - yc] */\
\
	/*...Block 3: Inputs y2,6,a,e: */\
				/* Twiddle mpy by E^4 = I: */\
	_r =_t13;	_t13=_t5 +_t14;	_t5 =_t5 -_t14;	/* t 5, 6: y2 + I.y6 */\
				_t14=_t6 -_r;	_t6 =_t6 +_r;	/* t13,14: y2 - I.y6 */\
\
	_r =(_t21-_t22)*ISRT2;	_t22=(_t21+_t22)*ISRT2;	_t21=_r;/* Twiddle mpy by E^2; could also do mpy by ISRT2 directly on t21,22,29,30: */\
	_r =(_t30+_t29)*ISRT2;	_i  =(_t30-_t29)*ISRT2;			/* Twiddle mpy by -E^6 is here... */\
	/* ...and get E^6 = I.E^2 = (I-1)/sqrt2 by flipping signs on _r,_i here: */\
	_t29=_t21+_r;			_t21=_t21-_r;	/* t21,22: E^2.ya + I.E^2.ye */\
	_t30=_t22+_i;			_t22=_t22-_i;	/* t29,30: E^2.ya - I.E^2.ye */\
\
	__A4r=_t5 +_t21;	__A4i=_t6 +_t22;	/* A4: [y2 + I.y6] +   [E^2.ya + I.E^2.ye] */\
	__A5r=_t5 -_t21;	__A5i=_t6 -_t22;	/* A5: [y2 + I.y6] -   [E^2.ya + I.E^2.ye] */\
\
	__A6r=_t13-_t30;	__A6i=_t14+_t29;	/* A6: [y2 - I.y6] + I.[E^2.ya - I.E^2.ye] */\
	__A7r=_t13+_t30;	__A7i=_t14-_t29;	/* A7: [y2 - I.y6] + I.[E^2.ya - I.E^2.ye] */\
\
	/*...Block 2: Inputs y1,5,9,d: */\
	_r =(_t11-_t12)*ISRT2;_i =(_t11+_t12)*ISRT2;		/* Twiddle mpy by E^2; could also do mpy by ISRT2 directly on t11,12 */\
	_t11=_t3 -_r;			_t3 =_t3 +_r;\
	_t12=_t4 -_i;			_t4 =_t4 +_i;\
\
	_r =_t19*__c - _t20*__s;	_t20=_t20*__c + _t19*__s;	_t19=_r;	/* Twiddle mpy by E^1 */\
	_r =_t27*__s - _t28*__c;	_i =_t28*__s + _t27*__c;		/* Twiddle mpy by E^3 */\
	_t27=_t19-_r;			_t19=_t19+_r;\
	_t28=_t20-_i;			_t20=_t20+_i;\
\
	__A8r=_t3+_t19;	__A8i=_t4+_t20;\
	__A9r=_t3-_t19;	__A9i=_t4-_t20;\
\
	__AAr=_t11-_t28;	__AAi=_t12+_t27;	/* mpy by E^4=i is inlined here... */\
	__ABr=_t11+_t28;	__ABi=_t12-_t27;\
\
	/*...Block 4: Inputs y3,7,b,f: */\
	_r =(_t16+_t15)*ISRT2;_i =(_t16-_t15)*ISRT2;		/* Twiddle mpy by -E^6 is here... */\
	_t15=_t7 +_r;			_t7 =_t7 -_r;			/* ...and get E^6=(i-1)/sqrt by flipping signs here. */\
	_t16=_t8 +_i;			_t8 =_t8 -_i;\
\
	_r =_t23*__s - _t24*__c;	_t24=_t24*__s + _t23*__c;	_t23=_r;	/* Twiddle mpy by E^3 */\
	_r =_t31*__c - _t32*__s;	_i =_t32*__c + _t31*__s;		/* Twiddle mpy by E^1 = -E^9... */\
	_t31=_t23+_r;			_t23=_t23-_r;			/* ...and get E^9 by flipping signs here. */\
	_t32=_t24+_i;			_t24=_t24-_i;			/* Note: t23+rt = t23*(s+1) */\
\
	__ACr=_t7+_t23;	__ACi=_t8+_t24;\
	__ADr=_t7-_t23;	__ADi=_t8-_t24;\
\
	__AEr=_t15-_t32;	__AEi=_t16+_t31;	/* mpy by E^4=i is inlined here... */\
	__AFr=_t15+_t32;	__AFi=_t16-_t31;\
}

/* Totals: 174 ADD, 84 MUL, which is the cost of the twiddle-less version plus 15 CMUL at [2 ADD, 4 MUL] each. */
//void ___RADIX_16_DIT_TWIDDLE() {}
#define RADIX_16_DIT_TWIDDLE(\
	__A0r,__A0i,__A1r,__A1i,__A2r,__A2i,__A3r,__A3i,__A4r,__A4i,__A5r,__A5i,__A6r,__A6i,__A7r,__A7i,__A8r,__A8i,__A9r,__A9i,__AAr,__AAi,__ABr,__ABi,__ACr,__ACi,__ADr,__ADi,__AEr,__AEi,__AFr,__AFi,\
	__c1 ,__s1 ,__c2 ,__s2 ,__c3 ,__s3 ,__c4 ,__s4 ,__c5 ,__s5 ,__c6 ,__s6 ,__c7 ,__s7 ,__c8 ,__s8 ,__c9 ,__s9 ,__cA ,__sA ,__cB ,__sB ,__cC ,__sC ,__cD ,__sD ,__cE ,__sE ,__cF ,__sF ,\
	__c,__s)\
{\
	double _t1,_t2,_t3,_t4,_t5,_t6,_t7,_t8,_t9,_t10,_t11,_t12,_t13,_t14,_t15,_t16,_t17,_t18,_t19,_t20,_t21,_t22,_t23,_t24,_t25,_t26,_t27,_t28,_t29,_t30,_t31,_t32;\
	double _r,_i;\
\
	/* Gather the needed data and do first set of four length-4 transforms: */\
	/*...Block 1: */\
	_t1 =__A0r;	_t2 =__A0i;\
	_r =__A1r;	_i =__A1i;\
	_t3 =_t1 -_r;		_t1 =_t1 +_r;	/* t1,2: x0 + x1 */\
	_t4 =_t2 -_i;		_t2 =_t2 +_i;	/* t3,4: x0 - x1 */\
\
	_t5 =__A2r;	_t6 =__A2i;\
	_r =__A3r;	_i =__A3i;\
	_t7 =_t5 -_r;		_t5 =_t5 +_r;	/* t5,6: x2 + x3 */\
	_t8 =_t6 -_i;		_t6 =_t6 +_i;	/* t7,8: x2 - x3 */\
\
	_r =_t5;	_t5 =_t1 -_r ;	_t1 =_t1 +_r;	/* t1,2: y0 = [x0 + x1] + [x2 + x3] */\
	_i =_t6;	_t6 =_t2 -_i ;	_t2 =_t2 +_i;	/* t5,6: y2 = [x0 + x1] - [x2 + x3] */\
\
	_r =_t7;	_t7 =_t3 -_t8;	_t3 =_t3 +_t8;	/* t3,4: y1 = [x0 - x1] - I.[x2 - x3] */\
				_t8 =_t4 +_r ;	_t4 =_t4 -_r ;	/* t7,8: y3 = [x0 - x1] + I.[x2 - x3] */\
										/* Remaining 3 Pass 1 DFTs same, just with x,y-indices +=4 for each */\
	/*...Block 2: */\
	_r =__A4r;	_t10=__A4i;	_t9 =_r;\
	_r =__A5r;	_i =__A5i;\
	_t11=_t9 -_r;		_t9 =_t9 +_r;\
	_t12=_t10-_i;		_t10=_t10+_i;\
\
	_r =__A6r;	_t14=__A6i;	_t13=_r;\
	_r =__A7r;	_i =__A7i;\
	_t15=_t13-_r;		_t13=_t13+_r;\
	_t16=_t14-_i;		_t14=_t14+_i;\
\
	_r =_t13;	_t13=_t9 -_r  ;	_t9 =_t9 +_r  ;	/* y4 = t9 ,10 */\
	_i =_t14;	_t14=_t10-_i  ;	_t10=_t10+_i  ;	/* y6 = t13,14 */\
\
	_r =_t15;	_t15=_t11-_t16;	_t11=_t11+_t16;	/* y5 = t11,12 */\
				_t16=_t12+_r  ;	_t12=_t12-_r  ;	/* y7 = t15,16 */\
\
	/*...Block 3: */\
	_r =__A8r;	_t18=__A8i;	_t17=_r;\
	_r =__A9r;	_i =__A9i;\
	_t19=_t17-_r;		_t17=_t17+_r;\
	_t20=_t18-_i;		_t18=_t18+_i;\
\
	_r =__AAr;	_t22=__AAi;	_t21=_r;\
	_r =__ABr;	_i =__ABi;\
	_t23=_t21-_r;		_t21=_t21+_r;\
	_t24=_t22-_i;		_t22=_t22+_i;\
\
	_r =_t21;	_t21=_t17-_r  ;	_t17=_t17+_r  ;	/* y8 = t17,18 */\
	_i =_t22;	_t22=_t18-_i  ;	_t18=_t18+_i  ;	/* y9 = t19,20 */\
\
	_r =_t23;	_t23=_t19-_t24;	_t19=_t19+_t24;	/* yA = t21,22 */\
				_t24=_t20+_r  ;	_t20=_t20-_r  ;	/* yB = t23,24 */\
\
	/*...Block 4: */\
	_r =__ACr;	_t26=__ACi;	_t25=_r;\
	_r =__ADr;	_i =__ADi;\
	_t27=_t25-_r;		_t25=_t25+_r;\
	_t28=_t26-_i;		_t26=_t26+_i;\
\
	_r =__AEr;	_t30=__AEi;	_t29=_r;\
	_r =__AFr;	_i =__AFi;\
	_t31=_t29-_r;		_t29=_t29+_r;\
	_t32=_t30-_i;		_t30=_t30+_i;\
\
	_r =_t29;	_t29=_t25-_r  ;	_t25=_t25+_r  ;	/* yC = t25,26 */\
	_i =_t30;	_t30=_t26-_i  ;	_t26=_t26+_i  ;	/* yD = t27,28 */\
\
	_r =_t31;	_t31=_t27-_t32;	_t27=_t27+_t32;	/* yE = t29,30 */\
				_t32=_t28+_r  ;	_t28=_t28-_r  ;	/* yF = t31,32 */\
\
/**************************************************************************************/\
/*...and now do four more radix-4 transforms, including the internal twiddle factors: */\
/**************************************************************************************/\
\
	/*...Block 1: _t1,9,17,25 */\
	_r =_t9;	_t9 =_t1 -_r;		_t1 =_t1 +_r;	/* y0 + y4 */\
	_i =_t10;	_t10=_t2 -_i;		_t2 =_t2 +_i;	/* y0 - y4 */\
\
	_r =_t25;	_t25=_t17-_r;		_t17=_t17+_r;	/* y8 + yC */\
	_i =_t26;	_t26=_t18-_i;		_t18=_t18+_i;	/* y8 - yC */\
\
	__A0r = _t1 +_t17;				__A0i = _t2 +_t18;				/* [(y0 + y4) + (y8 + yC)] */\
	_t1   = _t1 -_t17;				_t2   = _t2 -_t18;\
	__A8r = _t1 *__c8 +_t2 *__s8;	__A8i = _t2 *__c8 -_t1 *__s8;	/* [(y0 + y4) - (y8 + yC)].w8 */\
	/* mpy by E^-4 = -I is inlined here... */\
	_r    = _t9 +_t26;				_i    = _t10-_t25;	/* [(y0 - y4) - I.(y8 - yC)] */\
	_t9   = _t9 -_t26;				_t10  = _t10+_t25;	/* [(y0 - y4) + I.(y8 - yC)] */\
	__A4r = _r  *__c4 +_i  *__s4;	__A4i = _i  *__c4 -_r  *__s4;	/* [(y0 - y4) - I.(y8 - yC)].w4 */\
	__ACr = _t9 *__cC +_t10*__sC;	__ACi = _t10*__cC -_t9 *__sC;	/* [(y0 - y4) + I.(y8 - yC)].wC */\
\
	/*...Block 3: _t5,13,21,29 */\
	_r =_t13;	_t13=_t5 -_t14;		_t5 =_t5 +_t14;		/* Twiddle mpy by E^4 = -I */\
				_t14=_t6 +_r  ;		_t6 =_t6 -_r  ;		/* y2 +- E^4.y6 = y2 -+ I.y6 */\
\
	_r = (_t22+_t21)*ISRT2;	_t22 = (_t22-_t21)*ISRT2; _t21=_r;	/* E^2 = (1-I).ISRT2, thus (x+I.y).E^2 = [(y+x) + I.(y-x)].ISRT2 */\
	_r = (_t29-_t30)*ISRT2;	_i   = (_t29+_t30)*ISRT2;			/* Twiddle mpy by -E^6 = (1+I).ISRT2, thus (x+I.y).E^2 = [(x-y) + I.(x+y)].ISRT2 is here... */\
	_t29=_t21+_r;		_t21=_t21-_r;			/* ...and get E^6 by flipping signs here. */\
	_t30=_t22+_i;		_t22=_t22-_i;			/* yA.E^2 +- yE.E^6 = yA.E^2 -+ yE.E^2.I */\
\
	_r    =_t5 +_t21;				_i    =_t6 +_t22;	/* (y2 + E^4.y6) + (yA.E^2 + yE.E^6) */\
	_t5   =_t5 -_t21;				_t6   =_t6 -_t22;	/* (y2 + E^4.y6) - (yA.E^2 + yE.E^6) */\
	__A2r = _r  *__c2 +_i  *__s2;	__A2i = _i  *__c2 -_r  *__s2;	/* [(y2 + E^4.y6) +   (yA.E^2 + yE.E^6)].w2 */\
	__AAr = _t5 *__cA +_t6 *__sA;	__AAi = _t6 *__cA -_t5 *__sA;	/* [(y2 + E^4.y6) -   (yA.E^2 + yE.E^6)].wA */\
	/* mpy by E^4 = -I is inlined here: */\
	_r   =_t13+_t30;				_i   =_t14-_t29;	/* (y2 - E^4.y6) - I.(yA.E^2 - yE.E^6) */\
	_t13 =_t13-_t30;				_t14 =_t14+_t29;	/* (y2 - E^4.y6) + I.(yA.E^2 - yE.E^6) */\
	__A6r = _r  *__c6 +_i  *__s6;	__A6i = _i  *__c6 -_r  *__s6;	/* [(y2 - E^4.y6) - I.(yA.E^2 - yE.E^6)].w6 */\
	__AEr = _t13*__cE +_t14*__sE;	__AEi = _t14*__cE -_t13*__sE;	/* [(y2 - E^4.y6) + I.(yA.E^2 - yE.E^6)].wE */\
\
	/*...Block 2: _t3,11,19,27 */\
	_r =(_t12+_t11)*ISRT2;	_i =(_t12-_t11)*ISRT2;		/* Twiddle mpy by E^2 = (1-I)*ISRT2 */\
	_t11=_t3 -_r;					_t3 =_t3 +_r;\
	_t12=_t4 -_i;					_t4 =_t4 +_i;		/* y1 +- E^2.y5 */\
	/* Twiddle muls by E^1 = (c-I.s) and E^3 = (s-I.c): */\
	_r  =_t19*__c + _t20*__s;		_t20=_t20*__c - _t19*__s; _t19=_r;	/* y9.E^1 */\
	_r  =_t27*__s + _t28*__c;		_i  =_t28*__s - _t27*__c;			/* yD.E^3 */\
	_t27=_t19-_r;					_t19=_t19+_r;\
	_t28=_t20-_i;					_t20=_t20+_i;		/* y9.E^1 +- yD.E^3 */\
\
	_r    =_t3 +_t19;				_i    =_t4 +_t20;	/* (y1 + E^2.y5) +   (y9.E^1 + yD.E^3) */\
	_t3   =_t3 -_t19;				_t4   =_t4 -_t20;	/* (y1 + E^2.y5) -   (y9.E^1 + yD.E^3) */\
	__A1r = _r  *__c1 +_i  *__s1;	__A1i = _i  *__c1 -_r  *__s1;	/* [(y1 + E^2.y5) +   (y9.E^1 + yD.E^3)].w1 */\
	__A9r = _t3 *__c9 +_t4 *__s9;	__A9i = _t4 *__c9 -_t3 *__s9;	/* [(y1 + E^2.y5) -   (y9.E^1 + yD.E^3)].w9 */\
	/* mpy by E^-4 = -I is inlined here: */\
	_r    =_t11+_t28;				_i    =_t12-_t27;	/* (y1 - E^2.y5) - I.(y9.E^1 - yD.E^3) */\
	_t11  =_t11-_t28;				_t12  =_t12+_t27;	/* (y1 - E^2.y5) + I.(y9.E^1 - yD.E^3) */\
	__A5r = _r  *__c5 +_i  *__s5;	__A5i = _i  *__c5 -_r  *__s5;	/* [(y1 - E^2.y5) - I.(y9.E^1 - yD.E^3)].w5 */\
	__ADr = _t11*__cD +_t12*__sD;	__ADi = _t12*__cD -_t11*__sD;	/* [(y1 - E^2.y5) + I.(y9.E^1 - yD.E^3)].wD */\
\
	/*...Block 4: _t7,15,23,31 */\
	_r =(_t15-_t16)*ISRT2;	_i =(_t15+_t16)*ISRT2;		/* Twiddle mpy by -E^6 = (1+I)*ISRT2 is here... */\
	_t15=_t7 +_r;					_t7 =_t7 -_r;		/* ...and get E^6 = -(1+I)*ISRT2 by flipping signs here. */\
	_t16=_t8 +_i;					_t8 =_t8 -_i;		/* y3 +- E^6.y7 */\
	/* Twiddle muls by E^3 = (s-I.c) and E^9 = -E^1 = -(c-I.s): */\
	_r  =_t23*__s + _t24*__c;		_t24=_t24*__s - _t23*__c; _t23=_r;	/* yB.E^3 */\
	_r  =_t31*__c + _t32*__s;		_i  =_t32*__c - _t31*__s;			/* yF.E^1 = yF.(c-I.s)... */\
	_t31=_t23+_r;					_t23=_t23-_r;		/* ... and get yF.E^9 = -yF.E^1 by flipping signs here. */\
	_t32=_t24+_i;					_t24=_t24-_i;		/* yB.E^3 +- yF.E^9 */\
\
	_r    =_t7 +_t23;				_i    =_t8 +_t24;	/* (y3 + E^6.y7) +   (yB.E^3 + yF.E^9) */\
	_t7   =_t7 -_t23;				_t8   =_t8 -_t24;	/* (y3 + E^6.y7) -   (yB.E^3 + yF.E^9) */\
	__A3r = _r  *__c3 +_i  *__s3;	__A3i = _i  *__c3 -_r  *__s3;	/* [(y3 + E^6.y7) +   (yB.E^3 + yF.E^9)].w3 */\
	__ABr = _t7 *__cB +_t8 *__sB;	__ABi = _t8 *__cB -_t7 *__sB;	/* [(y3 + E^6.y7) +   (yB.E^3 + yF.E^9)].wB */\
	/* mpy by E^-4 = -I is inlined here: */\
	_r    =_t15+_t32;				_i    =_t16-_t31;\
	_t15  =_t15-_t32;				_t16  =_t16+_t31;\
	__A7r = _r  *__c7 +_i  *__s7;	__A7i = _i  *__c7 -_r  *__s7;	/* [(y3 - E^6.y7) - I.(yB.E^3 - yF.E^9)].w7 */\
	__AFr = _t15*__cF +_t16*__sF;	__AFi = _t16*__cF -_t15*__sF;	/* [(y3 - E^6.y7) + I.(yB.E^3 - yF.E^9)].wF */\
}


// With-twiddles out-of-place analog of above twiddleless DIF macro: 15 nontrivial complex input twiddles E1-f [E0 assumed = 1],
// The DIF version of this macro processes the twiddles in bit-reversed order: 0,8,4,c,2,a,6,e,1,9,5,d,3,b,7,f.
//void ___RADIX_16_DIF_TWIDDLE_OOP() {}
#define RADIX_16_DIF_TWIDDLE_OOP(\
	/* Inputs: */\
	__A0r,__A0i,__A1r,__A1i,__A2r,__A2i,__A3r,__A3i,__A4r,__A4i,__A5r,__A5i,__A6r,__A6i,__A7r,__A7i,__A8r,__A8i,__A9r,__A9i,__AAr,__AAi,__ABr,__ABi,__ACr,__ACi,__ADr,__ADi,__AEr,__AEi,__AFr,__AFi,\
	/* Outputs: */\
	__B0r,__B0i,__B1r,__B1i,__B2r,__B2i,__B3r,__B3i,__B4r,__B4i,__B5r,__B5i,__B6r,__B6i,__B7r,__B7i,__B8r,__B8i,__B9r,__B9i,__BAr,__BAi,__BBr,__BBi,__BCr,__BCi,__BDr,__BDi,__BEr,__BEi,__BFr,__BFi,\
	/* Twiddles: */\
	__c8 ,__s8 ,__c4 ,__s4 ,__cC ,__sC ,__c2 ,__s2 ,__cA ,__sA ,__c6 ,__s6 ,__cE ,__sE ,__c1 ,__s1 ,__c9 ,__s9 ,__c5 ,__s5 ,__cD ,__sD ,__c3 ,__s3 ,__cB ,__sB ,__c7 ,__s7 ,__cF ,__sF ,\
	__c,__s)\
{\
	double _t1,_t2,_t3,_t4,_t5,_t6,_t7,_t8,_t9,_t10,_t11,_t12,_t13,_t14,_t15,_t16,_t17,_t18,_t19,_t20,_t21,_t22,_t23,_t24,_t25,_t26,_t27,_t28,_t29,_t30,_t31,_t32;\
	double _r,_i;\
	/* Gather the needed data (16 64-bit complex, i.e. 32 64-bit reals) and do first set of four length-4 transforms: */\
	/*...Block 1: */\
	_t1 =__A0r;						_t2 =__A0i;\
	_r =__A8r*__c8 -__A8i*__s8 ;	_i =__A8i*__c8 +__A8r*__s8 ;\
	_t3 =_t1 -_r;						_t1 =_t1 +_r;\
	_t4 =_t2 -_i;						_t2 =_t2 +_i;\
\
	_t5 =__A4r*__c4 -__A4i*__s4 ;	_t6 =__A4i*__c4 +__A4r*__s4 ;\
	_r =__ACr*__cC -__ACi*__sC ;	_i =__ACi*__cC +__ACr*__sC ;\
	_t7 =_t5 -_r;							_t5 =_t5 +_r;\
	_t8 =_t6 -_i;							_t6 =_t6 +_i;\
\
	_r =_t5;	_t5 =_t1 -_r;					_t1 =_t1 +_r;\
	_i =_t6;	_t6 =_t2 -_i;					_t2 =_t2 +_i;\
\
	_r =_t7;	_t7 =_t3 +_t8;					_t3 =_t3 -_t8;\
			_t8 =_t4 -_r;					_t4 =_t4 +_r;\
\
	/*...Block 2: */\
	_t9 =__A2r*__c2 -__A2i*__s2;	_t10=__A2i*__c2 +__A2r*__s2;\
	_r =__AAr*__cA -__AAi*__sA ;	_i =__AAi*__cA +__AAr*__sA ;\
	_t11=_t9 -_r;						_t9 =_t9 +_r;\
	_t12=_t10-_i;						_t10=_t10+_i;\
\
	_t13=__A6r*__c6 -__A6i*__s6 ;	_t14=__A6i*__c6 +__A6r*__s6;\
	_r =__AEr*__cE -__AEi*__sE ;	_i =__AEi*__cE +__AEr*__sE ;\
	_t15=_t13-_r;							_t13=_t13+_r;\
	_t16=_t14-_i;							_t14=_t14+_i;\
\
	_r =_t13;	_t13=_t9 -_r;				_t9 =_t9 +_r;\
	_i =_t14;	_t14=_t10-_i;				_t10=_t10+_i;\
\
	_r =_t15;	_t15=_t11+_t16;			_t11=_t11-_t16;\
				_t16=_t12-_r;				_t12=_t12+_r;\
\
	/*...Block 3: */\
	_t17=__A1r*__c1 -__A1i*__s1 ;	_t18=__A1i*__c1 +__A1r*__s1;\
	_r =__A9r*__c9 -__A9i*__s9 ;	_i =__A9i*__c9 +__A9r*__s9;\
	_t19=_t17-_r;						_t17=_t17+_r;\
	_t20=_t18-_i;						_t18=_t18+_i;\
\
	_t21=__A5r*__c5 -__A5i*__s5 ;	_t22=__A5i*__c5 +__A5r*__s5;\
	_r =__ADr*__cD -__ADi*__sD ;	_i =__ADi*__cD +__ADr*__sD ;\
	_t23=_t21-_r;							_t21=_t21+_r;\
	_t24=_t22-_i;							_t22=_t22+_i;\
\
	_r =_t21;	_t21=_t17-_r;				_t17=_t17+_r;\
	_i =_t22;	_t22=_t18-_i;				_t18=_t18+_i;\
\
	_r =_t23;	_t23=_t19+_t24;			_t19=_t19-_t24;\
				_t24=_t20-_r;				_t20=_t20+_r;\
\
	/*...Block 4: */\
	_t25=__A3r*__c3 -__A3i*__s3 ;	_t26=__A3i*__c3 +__A3r*__s3;\
	_r =__ABr*__cB -__ABi*__sB ;	_i =__ABi*__cB +__ABr*__sB ;\
	_t27=_t25-_r;						_t25=_t25+_r;\
	_t28=_t26-_i;						_t26=_t26+_i;\
\
	_t29=__A7r*__c7 -__A7i*__s7 ;	_t30=__A7i*__c7 +__A7r*__s7;\
	_r =__AFr*__cF -__AFi*__sF ;	_i =__AFi*__cF +__AFr*__sF ;\
	_t31=_t29-_r;							_t29=_t29+_r;\
	_t32=_t30-_i;							_t30=_t30+_i;\
\
	_r =_t29;	_t29=_t25-_r;				_t25=_t25+_r;\
	_i =_t30;	_t30=_t26-_i;				_t26=_t26+_i;\
\
	_r =_t31;	_t31=_t27+_t32;			_t27=_t27-_t32;\
				_t32=_t28-_r;				_t28=_t28+_r;\
\
/**************************************************************************************/\
/*...and now do four more radix-4 transforms, including the internal twiddle factors: */\
/**************************************************************************************/\
\
	/*...Block 1: _t1,9,17,25 */\
	_r =_t9;	_t9 =_t1 -_r;	_t1 =_t1 +_r;\
	_i =_t10;_t10=_t2 -_i;	_t2 =_t2 +_i;\
\
	_r =_t25;_t25=_t17-_r;	_t17=_t17+_r;\
	_i =_t26;_t26=_t18-_i;	_t18=_t18+_i;\
\
	__B0r=_t1+_t17;	__B0i=_t2+_t18;\
	__B1r=_t1-_t17;	__B1i=_t2-_t18;\
\
	__B2r=_t9 -_t26;	__B2i=_t10+_t25;	/* mpy by E^4=i is inlined here... */\
	__B3r=_t9 +_t26;	__B3i=_t10-_t25;\
\
	/*...Block 3: _t5,13,21,29 */\
	_r =_t13;_t13=_t5 +_t14;_t5 =_t5 -_t14;		/* Twiddle mpy by E^4 = I */\
			_t14=_t6 -_r;	_t6 =_t6 +_r;\
\
	_r =(_t21-_t22)*ISRT2;_t22=(_t21+_t22)*ISRT2;	_t21=_r;	/* Twiddle mpy by E^2; could also do mpy by ISRT2 directly on t21,22,29,30 */\
	_r =(_t30+_t29)*ISRT2;_i =(_t30-_t29)*ISRT2;		/* Twiddle mpy by -E^6 is here... */\
	_t29=_t21+_r;			_t21=_t21-_r;			/* ...and get E^6=(i-1)/sqrt by flipping signs here. */\
	_t30=_t22+_i;			_t22=_t22-_i;\
\
	__B4r=_t5+_t21;	__B4i=_t6+_t22;\
	__B5r=_t5-_t21;	__B5i=_t6-_t22;\
\
	__B6r=_t13-_t30;	__B6i=_t14+_t29;	/* mpy by E^4=i is inlined here... */\
	__B7r=_t13+_t30;	__B7i=_t14-_t29;\
\
	/*...Block 2: _t3,11,19,27 */\
	_r =(_t11-_t12)*ISRT2;_i =(_t11+_t12)*ISRT2;		/* Twiddle mpy by E^2; could also do mpy by ISRT2 directly on t11,12 */\
	_t11=_t3 -_r;			_t3 =_t3 +_r;\
	_t12=_t4 -_i;			_t4 =_t4 +_i;\
\
	_r =_t19*__c - _t20*__s;	_t20=_t20*__c + _t19*__s;	_t19=_r;	/* Twiddle mpy by E^1 */\
	_r =_t27*__s - _t28*__c;	_i =_t28*__s + _t27*__c;		/* Twiddle mpy by E^3 */\
	_t27=_t19-_r;			_t19=_t19+_r;\
	_t28=_t20-_i;			_t20=_t20+_i;\
\
	__B8r=_t3+_t19;	__B8i=_t4+_t20;\
	__B9r=_t3-_t19;	__B9i=_t4-_t20;\
\
	__BAr=_t11-_t28;	__BAi=_t12+_t27;	/* mpy by E^4=i is inlined here... */\
	__BBr=_t11+_t28;	__BBi=_t12-_t27;\
\
	/*...Block 4: _t7,15,23,31 */\
	_r =(_t16+_t15)*ISRT2;_i =(_t16-_t15)*ISRT2;		/* Twiddle mpy by -E^6 is here... */\
	_t15=_t7 +_r;			_t7 =_t7 -_r;			/* ...and get E^6=(i-1)/sqrt by flipping signs here. */\
	_t16=_t8 +_i;			_t8 =_t8 -_i;\
\
	_r =_t23*__s - _t24*__c;	_t24=_t24*__s + _t23*__c;	_t23=_r;	/* Twiddle mpy by E^3 */\
	_r =_t31*__c - _t32*__s;	_i =_t32*__c + _t31*__s;		/* Twiddle mpy by E^1 = -E^9... */\
	_t31=_t23+_r;			_t23=_t23-_r;			/* ...and get E^9 by flipping signs here. */\
	_t32=_t24+_i;			_t24=_t24-_i;			/* Note: t23+rt = t23*(s+1) */\
\
	__BCr=_t7+_t23;	__BCi=_t8+_t24;\
	__BDr=_t7-_t23;	__BDi=_t8-_t24;\
\
	__BEr=_t15-_t32;	__BEi=_t16+_t31;	/* mpy by E^4=i is inlined here... */\
	__BFr=_t15+_t32;	__BFi=_t16-_t31;\
}


// With-twiddles out-of-place analog of above twiddleless DIT macro: 15 nontrivial complex input twiddles E1-f [E0 assumed = 1],
// The DIT version of this macro processes the twiddles in-order.
// NOTE: SINCE THIS MACRO IS SPECIFICALLY DESIGNED AS THE 2ND-PASS OF LARGE-POWER-OF-2-TWIDDLELESS DFT SYNTHESIS, THE
// "TWIDDLES" HERE ARE PURE OF THE DFT-INTERNAL VARIETY, AND THUS APPLIED TO THE INPUTS, JUST AS FOR THE ABOVE DIF COUNTERPART.
//void ___RADIX_16_DIT_TWIDDLE_OOP() {}
#define RADIX_16_DIT_TWIDDLE_OOP(\
	/* Inputs: */\
	__A0r,__A0i,__A1r,__A1i,__A2r,__A2i,__A3r,__A3i,__A4r,__A4i,__A5r,__A5i,__A6r,__A6i,__A7r,__A7i,__A8r,__A8i,__A9r,__A9i,__AAr,__AAi,__ABr,__ABi,__ACr,__ACi,__ADr,__ADi,__AEr,__AEi,__AFr,__AFi,\
	/* Outputs: */\
	__B0r,__B0i,__B1r,__B1i,__B2r,__B2i,__B3r,__B3i,__B4r,__B4i,__B5r,__B5i,__B6r,__B6i,__B7r,__B7i,__B8r,__B8i,__B9r,__B9i,__BAr,__BAi,__BBr,__BBi,__BCr,__BCi,__BDr,__BDi,__BEr,__BEi,__BFr,__BFi,\
	/* Twiddles: */\
	__c1 ,__s1 ,__c2 ,__s2 ,__c3 ,__s3 ,__c4 ,__s4 ,__c5 ,__s5 ,__c6 ,__s6 ,__c7 ,__s7 ,__c8 ,__s8 ,__c9 ,__s9 ,__cA ,__sA ,__cB ,__sB ,__cC ,__sC ,__cD ,__sD ,__cE ,__sE ,__cF ,__sF ,\
	__c,__s)\
{\
	double _t1,_t2,_t3,_t4,_t5,_t6,_t7,_t8,_t9,_t10,_t11,_t12,_t13,_t14,_t15,_t16,_t17,_t18,_t19,_t20,_t21,_t22,_t23,_t24,_t25,_t26,_t27,_t28,_t29,_t30,_t31,_t32;\
	double _r,_i;\
	/*...Block 1: */\
	_t1 =__A0r;						_t2 =__A0i;\
	_r = __A1r*__c1 +__A1i*__s1;	_i  =__A1i*__c1 -__A1r*__s1;\
	_t3 =_t1 -_r;	_t1 += _r;\
	_t4 =_t2 -_i;	_t2 += _i;\
\
	_t5 =__A2r*__c2 +__A2i*__s2;	_t6 =__A2i*__c2 -__A2r*__s2;\
	_r  =__A3r*__c3 +__A3i*__s3;	_i  =__A3i*__c3 -__A3r*__s3;\
	_t7 =_t5 -_r;	_t5 += _r;\
	_t8 =_t6 -_i;	_t6 += _i;\
\
	_r =_t5;		_t5 =_t1 -_r;	_t1 =_t1 +_r;\
	_i =_t6;		_t6 =_t2 -_i;	_t2 =_t2 +_i;\
\
	_r =_t7;		_t7 =_t3 -_t8;	_t3 =_t3 +_t8;\
				_t8 =_t4 +_r;	_t4 =_t4 -_r;\
\
	/*...Block 2: */\
	_t9 =__A4r*__c4 +__A4i*__s4;	_t10 =__A4i*__c4 -__A4r*__s4;\
	_r  =__A5r*__c5 +__A5i*__s5;	_i   =__A5i*__c5 -__A5r*__s5;\
	_t11=_t9  -_r;	_t9  += _r;\
	_t12=_t10 -_i;	_t10 += _i;\
\
	_t13 =__A6r*__c6 +__A6i*__s6;	_t14 =__A6i*__c6 -__A6r*__s6;\
	_r   =__A7r*__c7 +__A7i*__s7;	_i   =__A7i*__c7 -__A7r*__s7;\
	_t15=_t13 -_r;	_t13 += _r;\
	_t16=_t14 -_i;	_t14 += _i;\
\
	_r =_t13;	_t13=_t9 -_r;		_t9 =_t9 +_r;\
	_i =_t14;	_t14=_t10-_i;		_t10=_t10+_i;\
\
	_r =_t15;	_t15=_t11-_t16;	_t11=_t11+_t16;\
				_t16=_t12+_r;		_t12=_t12-_r;\
\
	/*...Block 3: */\
	_t17 =__A8r*__c8 +__A8i*__s8;	_t18 =__A8i*__c8 -__A8r*__s8;\
	_r   =__A9r*__c9 +__A9i*__s9;	_i   =__A9i*__c9 -__A9r*__s9;\
	_t19=_t17 -_r;	_t17 += _r;\
	_t20=_t18 -_i;	_t18 += _i;\
\
	_t21 =__AAr*__cA +__AAi*__sA;	_t22 =__AAi*__cA -__AAr*__sA;\
	_r   =__ABr*__cB +__ABi*__sB;	_i   =__ABi*__cB -__ABr*__sB;\
	_t23=_t21 -_r;	_t21 += _r;\
	_t24=_t22 -_i;	_t22 += _i;\
\
	_r =_t21;	_t21=_t17-_r;		_t17=_t17+_r;\
	_i =_t22;	_t22=_t18-_i;		_t18=_t18+_i;\
\
	_r =_t23;	_t23=_t19-_t24;	_t19=_t19+_t24;\
				_t24=_t20+_r;		_t20=_t20-_r;\
\
	/*...Block 4: */\
	_t25 =__ACr*__cC +__ACi*__sC;	_t26 =__ACi*__cC -__ACr*__sC;\
	_r   =__ADr*__cD +__ADi*__sD;	_i   =__ADi*__cD -__ADr*__sD;\
	_t27=_t25 -_r;	_t25 += _r; \
	_t28=_t26 -_i;	_t26 += _i;\
\
	_t29 =__AEr*__cE +__AEi*__sE;	_t30 =__AEi*__cE -__AEr*__sE;\
	_r   =__AFr*__cF +__AFi*__sF;	_i   =__AFi*__cF -__AFr*__sF;\
	_t31=_t29 -_r;	_t29 += _r;\
	_t32=_t30 -_i;	_t30 += _i;\
\
	_r =_t29;	_t29=_t25-_r;		_t25=_t25+_r;\
	_i =_t30;	_t30=_t26-_i;		_t26=_t26+_i;\
\
	_r =_t31;	_t31=_t27-_t32;	_t27=_t27+_t32;\
				_t32=_t28+_r;		_t28=_t28-_r;\
\
/**************************************************************************************/\
/*...and now do four more radix-4 transforms, including the internal twiddle factors: */\
/**************************************************************************************/\
\
	/*...Block 1: _t1,9,17,25	*/\
	_r =_t9 ;	_t9 =_t1 -_r;	_t1 =_t1 +_r;\
	_i =_t10;	_t10=_t2 -_i;	_t2 =_t2 +_i;\
\
	_r =_t25;	_t25=_t17-_r;	_t17=_t17+_r;\
	_i =_t26;	_t26=_t18-_i;	_t18=_t18+_i;\
\
	__B0r=_t1+_t17;			__B0i =_t2+_t18;\
	__B8r=_t1-_t17;			__B8i =_t2-_t18;\
\
	__B4r=_t9 +_t26;			__B4i =_t10-_t25;	/* mpy by E^-4 = -I is inlined here...	*/\
	__BCr=_t9 -_t26;			__BCi =_t10+_t25;\
\
	/*...Block 3: _t5,13,21,29	*/\
	_r =_t13;_t13=_t5 -_t14;	_t5 =_t5 +_t14;	/* Twiddle mpy by E^4 =-I	*/\
			_t14=_t6 +_r;		_t6 =_t6 -_r;\
\
	_r =(_t22+_t21)*ISRT2;	_t22=(_t22-_t21)*ISRT2;	_t21=_r;	/* Twiddle mpy by E^-2	*/\
	_r =(_t29-_t30)*ISRT2;	_i =(_t29+_t30)*ISRT2;	/* Twiddle mpy by E^2 = -E^-6 is here...	*/\
	_t29=_t21+_r;				_t21=_t21-_r;				/* ...and get E^-6 by flipping signs here.	*/\
	_t30=_t22+_i;				_t22=_t22-_i;\
\
	__B2r=_t5+_t21;			__B2i =_t6+_t22;\
	__BAr=_t5-_t21;			__BAi =_t6-_t22;\
\
	__B6r=_t13+_t30;			__B6i =_t14-_t29;			/* mpy by E^-4 =-I is inlined here...	*/\
	__BEr=_t13-_t30;			__BEi =_t14+_t29;\
\
	/*...Block 2: _t3,11,19,27	*/\
	_r =(_t12+_t11)*ISRT2;	_i =(_t12-_t11)*ISRT2;	/* Twiddle mpy by E^-2	*/\
	_t11=_t3 -_r;				_t3 =_t3 +_r;\
	_t12=_t4 -_i;				_t4 =_t4 +_i;\
\
	_r =_t19*__c + _t20*__s;	_t20=_t20*__c - _t19*__s;	_t19=_r;	/* Twiddle mpy by E^-1	*/\
	_r =_t27*__s + _t28*__c;	_i =_t28*__s - _t27*__c;	/* Twiddle mpy by E^-3	*/\
	_t27=_t19-_r;				_t19=_t19+_r;\
	_t28=_t20-_i;				_t20=_t20+_i;\
\
	__B1r=_t3+_t19;			__B1i =_t4+_t20;\
	__B9r=_t3-_t19;			__B9i =_t4-_t20;\
\
	__B5r=_t11+_t28;			__B5i =_t12-_t27;			/* mpy by E^-4 =-I is inlined here...	*/\
	__BDr=_t11-_t28;			__BDi =_t12+_t27;\
\
	/*...Block 4: _t7,15,23,31	*/\
	_r =(_t15-_t16)*ISRT2;	_i =(_t15+_t16)*ISRT2;	/* Twiddle mpy by E^2 = -E^-6 is here...	*/\
	_t15=_t7 +_r;				_t7 =_t7 -_r;			/* ...and get E^6=(i-1)/sqrt2 by flipping signs here.	*/\
	_t16=_t8 +_i;				_t8 =_t8 -_i;\
\
	_r =_t23*__s + _t24*__c;	_t24=_t24*__s - _t23*__c;	_t23=_r;	/* Twiddle mpy by E^-3	*/\
	_r =_t31*__c + _t32*__s;	_i =_t32*__c - _t31*__s;	/* Twiddle mpy by E^-1 = -E^-9...	*/\
	_t31=_t23+_r;				_t23=_t23-_r;			/* ...and get E^9 by flipping signs here.	*/\
	_t32=_t24+_i;				_t24=_t24-_i;\
\
	__B3r=_t7+_t23;			__B3i =_t8+_t24;\
	__BBr=_t7-_t23;			__BBi =_t8-_t24;\
\
	__B7r=_t15+_t32;			__B7i =_t16-_t31;			/* mpy by E^-4 = -I is inlined here...	*/\
	__BFr=_t15-_t32;			__BFi =_t16+_t31;\
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
//void ___RADIX_16_DIF_FMA4() {}
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
/**************************************************************************************/\
/*...and now do four more radix-4 transforms, including the internal twiddle factors: */\
/**************************************************************************************/\
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
//void ___RADIX_16_DIF_FMA() {}
#define RADIX_16_DIF_FMA(\
	__A0r,__A0i,__A1r,__A1i,__A2r,__A2i,__A3r,__A3i,__A4r,__A4i,__A5r,__A5i,__A6r,__A6i,__A7r,__A7i,__A8r,__A8i,__A9r,__A9i,__AAr,__AAi,__ABr,__ABi,__ACr,__ACi,__ADr,__ADi,__AEr,__AEi,__AFr,__AFi,\
	__B0r,__B0i,__B1r,__B1i,__B2r,__B2i,__B3r,__B3i,__B4r,__B4i,__B5r,__B5i,__B6r,__B6i,__B7r,__B7i,__B8r,__B8i,__B9r,__B9i,__BAr,__BAi,__BBr,__BBi,__BCr,__BCi,__BDr,__BDi,__BEr,__BEi,__BFr,__BFi,\
	__c1 ,__t1 ,__c2 ,__t2 ,__c31,__t3 ,__c4 ,__t4 ,__c51,__t5 ,__c62,__t6 ,__c73,__t7 ,__c8 ,__t8 ,__c91,__t9 ,__cA2,__tA ,__cB3,__tB ,__cC4,__tC ,__cD5,__tD ,__cE6,__tE ,__cF7,__tF ,\
	__c1_c,__sc,__c1i2,__c2i2)\
{\
	double t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31;\
	double _a,_b,_c,_d,_e,_f;\
	/* The divides in the following sincos data will be precomputed in an optimized version. */\
	/* Sine terms here replaced by (and assumed to hold) tangents; cosine terms defined like so: */\
	/*	double __c1_c = __c1*__c;	- In SIMD code, this will take the place of __c */\
	/*	double __sc   = __s/__c;	- In SIMD code, this will take the place of __s */\
	/*	double __c1i2 = __c1*ISRT2; */\
	/*	double __c2i2 = __c2*ISRT2; */\
	/*	double __c1 ;*/\
	/*	double __c2 ;*/\
	/*	double __c31 = __c3/__c1; */\
	/*	double __c4 ;*/\
	/*	double __c51 = __c5/__c1; */\
	/*	double __c62 = __c6/__c2; */\
	/*	double __c73 = __c7/__c3; */\
	/*	double __c8 ;*/\
	/*	double __c91 = __c9/__c1; */\
	/*	double __cA2 = __cA/__c2; */\
	/*	double __cB3 = __cB/__c3; */\
	/*	double __cC4 = __cC/__c4; */\
	/*	double __cD5 = __cD/__c5; */\
	/*	double __cE6 = __cE/__c6; */\
	/*	double __cF7 = __cF/__c7; */\
	/* In addition to the 30 derived multipliers, we need the unmodified __c1,2,4,8, for a total of 34. */\
\
/* Gather the needed data and do first set of four length-4 transforms: */\
/* Remember that for the FMA231 macro (defined in float_intrin.h), the result overwrites the rightmost input! */\
	/*...Block 1: */\
	t04 =t06 =__A8r;			t05 =__A8i;\
	t00 =__A0r;					t01 =__A0i;/* load __c8,__t8 into pair of regs */\
	FNMA231(   t05,__t8,t04);	 FMA231(   t06,__t8,t05);/* Try to queue up as many FMAs as registers allow */\
	_a =__A4r;					_b =__A4i;/* load __cC4,__t4,__tC into trio of regs */\
	FNMA231(__A4i,__t4,_a);		 FMA231(__A4r,__t4,_b);\
	t06 =__ACr;					t07 =__ACi;\
	FNMA231(__ACi,__tC,t06);	 FMA231(__ACr,__tC,t07);\
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
	FNMA231(__A2i,__t2,t08 );	 FMA231(__A2r,__t2,t09);\
	t12 =__AAr;					t13 =__AAi;\
	FNMA231(__AAi,__tA,t12 );	 FMA231(__AAr,__tA,t13 );\
	_a=__A6r;					_b=__A6i;\
	FNMA231(__A6i,__t6,_a);		 FMA231(__A6r,__t6,_b);\
	t14 =__AEr;					t15 =__AEi;\
	FNMA231(__AEi,__tE,t14 );	 FMA231(__AEr,__tE,t15 );\
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
	FNMA231(__A1i,__t1,t16);	 FMA231(__A1r,__t1,t17);\
	t20 =__A9r;	t21 =__A9i;\
	FNMA231(__A9i,__t9,t20 );	 FMA231(__A9r,__t9,t21 );\
	_a=__A5r;	_b=__A5i;\
	FNMA231(__A5i,__t5,_a);		 FMA231(__A5r,__t5,_b);\
	t22 =__ADr;	t23 =__ADi;\
	FNMA231(__ADi,__tD,t22 );	 FMA231(__ADr,__tD,t23 );\
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
	FNMA231(__A3i,__t3,t24);	 FMA231(__A3r,__t3,t25);\
	t28 =__ABr;	t29 =__ABi;\
	FNMA231(__ABi,__tB,t28 );	 FMA231(__ABr,__tB,t29 );\
	_a=__A7r;	_b=__A7i;\
	FNMA231(__A7i,__t7,_a);		 FMA231(__A7r,__t7,_b);\
	t30 =__AFr;	t31 =__AFi;\
	FNMA231(__AFi,__tF,t30 );	 FMA231(__AFr,__tF,t31 );\
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
/**************************************************************************************/\
/*...and now do four more radix-4 transforms, including the internal twiddle factors: */\
/**************************************************************************************/\
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

// This is a DIT analog of the above, specially designed for the 2nd-pass phase of large power-of-2 DIT DFTs.
// We lump it with the DIF macros because like those, it has the twiddles applied to the inputs, and in fact
// is based on the above DIF macro, with sign flips where needed to effect roots-of-unity conjugacy:
//void ___RADIX_16_DIT_FMA_PRETWIDDLE() {}
#define RADIX_16_DIT_FMA_PRETWIDDLE(\
	__A0r,__A0i,__A1r,__A1i,__A2r,__A2i,__A3r,__A3i,__A4r,__A4i,__A5r,__A5i,__A6r,__A6i,__A7r,__A7i,__A8r,__A8i,__A9r,__A9i,__AAr,__AAi,__ABr,__ABi,__ACr,__ACi,__ADr,__ADi,__AEr,__AEi,__AFr,__AFi,\
	__B0r,__B0i,__B1r,__B1i,__B2r,__B2i,__B3r,__B3i,__B4r,__B4i,__B5r,__B5i,__B6r,__B6i,__B7r,__B7i,__B8r,__B8i,__B9r,__B9i,__BAr,__BAi,__BBr,__BBi,__BCr,__BCi,__BDr,__BDi,__BEr,__BEi,__BFr,__BFi,\
	__c1 ,__t1 ,__c2 ,__t2 ,__c31,__t3 ,__c4 ,__t4 ,__c51,__t5 ,__c62,__t6 ,__c73,__t7 ,__c8 ,__t8 ,__c91,__t9 ,__cA2,__tA ,__cB3,__tB ,__cC4,__tC ,__cD5,__tD ,__cE6,__tE ,__cF7,__tF ,\
	__c1_c,__sc,__c1i2,__c2i2)\
{\
	double t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31;\
	double _a,_b,_c,_d,_e,_f;\
\
/* Gather the needed data and do first set of four length-4 transforms: */\
/* Remember that for the FMA231 macro (defined in float_intrin.h), the result overwrites the rightmost input! */\
	/*...Block 1: */\
	t04 =t06 =__A1r;			t05 =__A1i;\
	t00 =__A0r;					t01 =__A0i;/* load __c8,__t8 into pair of regs */\
	 FMA231(  t05,__t8,t04);	FNMA231(  t06,__t8,t05);/* Try to queue up as many FMAs as registers allow */\
	_a =__A2r;					_b =__A2i;/* load __cC4,__t4,__tC into trio of regs */\
	 FMA231(__A2i,__t4,_a);		FNMA231(__A2r,__t4,_b);\
	t06 =__A3r;					t07 =__A3i;\
	 FMA231(__A3i,__tC,t06);	FNMA231(__A3r,__tC,t07);\
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
	 FMA231(_d ,__c4 ,t02);		FNMA231(_d ,__c4 ,t06);\
	FNMA231(_c ,__c4 ,t03);		 FMA231(_c ,__c4 ,t07);\
	/* [In ASM version, write outputs into local store here] */\
\
	/*...Block 2: */\
	t08 =__A4r;					t09=__A4i;\
	 FMA231(__A4i,__t2,t08 );	FNMA231(__A4r,__t2,t09);\
	t12 =__A5r;					t13 =__A5i;\
	 FMA231(__A5i,__tA,t12 );	FNMA231(__A5r,__tA,t13 );\
	_a=__A6r;					_b=__A6i;\
	 FMA231(__A6i,__t6,_a);		FNMA231(__A6r,__t6,_b);\
	t14 =__A7r;					t15 =__A7i;\
	 FMA231(__A7i,__tE,t14 );	FNMA231(__A7r,__tE,t15 );\
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
	 FMA231(_d,__c62,t10);		FNMA231(_d,__c62,t14);\
	FNMA231(_c,__c62,t11);		 FMA231(_c,__c62,t15);\
	/* t8-16 need remultiply by c2: */\
	/* [In ASM version, write outputs into local store here] */\
\
	/*...Block 3: */\
	t16=__A8r;	t17=__A8i;\
	 FMA231(__A8i,__t1,t16);	FNMA231(__A8r,__t1,t17);\
	t20 =__A9r;	t21 =__A9i;\
	 FMA231(__A9i,__t9,t20 );	FNMA231(__A9r,__t9,t21 );\
	_a=__AAr;	_b=__AAi;\
	 FMA231(__AAi,__t5,_a);		FNMA231(__AAr,__t5,_b);\
	t22 =__ABr;	t23 =__ABi;\
	 FMA231(__ABi,__tD,t22 );	FNMA231(__ABr,__tD,t23 );\
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
	 FMA231(_d,__c51,t18);		FNMA231(_d,__c51,t22);\
	FNMA231(_c,__c51,t19);		 FMA231(_c,__c51,t23);\
	/* t16-23 need remultiply by c1: */\
	/* [In ASM version, write outputs into local store here] */\
\
	/*...Block 4: */\
	t24=__ACr;	t25=__ACi;\
	 FMA231(__ACi,__t3,t24);	FNMA231(__ACr,__t3,t25);\
	t28 =__ADr;	t29 =__ADi;\
	 FMA231(__ADi,__tB,t28 );	FNMA231(__ADr,__tB,t29 );\
	_a=__AEr;	_b=__AEi;\
	 FMA231(__AEi,__t7,_a);		FNMA231(__AEr,__t7,_b);\
	t30 =__AFr;	t31 =__AFi;\
	 FMA231(__AFi,__tF,t30 );	FNMA231(__AFr,__tF,t31 );\
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
	 FMA231(_d,__c73,t26);		FNMA231(_d,__c73,t30);\
	FNMA231(_c,__c73,t27);		 FMA231(_c,__c73,t31);\
	/* t24-31 need remultiply by c3: */\
	/* [In ASM version, write outputs into local store here] */\
\
/**************************************************************************************/\
/*...and now do four more radix-4 transforms, including the internal twiddle factors: */\
/**************************************************************************************/\
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
	__B8r= _e ;		__B8i= _f ;\
	__BCr= t08;		__BCi= t09;/* Swap 4/C outputs for DIT */\
	__B4r= _a ;		__B4i= _b ;\
\
	/*...Block 3: t4,12,20,28 */\
	_c = t20; _d = t29;\
	 FMA231(t21,1.0,_c );		 FMS231(t21,1.0,t20);		/* twiddle mpy by E^2; ISRT2 absorbed into various multipliers */\
	 FMS231(t28,1.0,_d );		 FMA231(t28,1.0,t29);		/* twiddle mpy by -E^6 is here... */\
	/* Swapped upper right and lower left FMA-pairs of octet below to hide latency: */\
	_a = t04; _b = t05;\
	t21 = _c; _e = t20;\
	FNMA231(t13,__c2 ,_a );		 FMA231(_d ,__c31,t21);		/* twiddle mpy by E^4 = I */\
	 FMA231(t12,__c2 ,_b );		 FMA231(t29,__c31,t20);\
	 FMA231(t13,__c2 ,t04);		FNMA231(_d ,__c31,_c );		/* ...and get E^6=[i-1]/sqrt2 by flipping signs here. */\
	FNMA231(t12,__c2 ,t05);		FNMA231(t29,__c31,_e );\
	t12 = _a; t13 = _b;\
	FNMA231(t20,__c1i2,t12);	 FMA231(t21,__c1i2,t13);	/* mpy by E^4=i is inlined here... */\
	 FMA231(t20,__c1i2,_a );	FNMA231(t21,__c1i2,_b );\
	_d = t04; t29 = t05;\
	 FMA231(_c ,__c1i2,t04);	 FMA231(_e ,__c1i2,t05);\
	FNMA231(_c ,__c1i2,_d );	FNMA231(_e ,__c1i2,t29);\
	 								_e = t11+t10;	_f = t11-t10;	/* [Block 2] twiddle mpy by -E^6 is here... */\
	/* Write outputs back to main array - Not sure why, but need apply 6/E swap here *and* then pairwise swap, i.e. 2A[6][E] => [2A][E6] => E62A: */\
	__BEr= t12;		__BEi= t13;\
	__B6r= _a ;		__B6i= _b ;\
	__B2r= t04;		__B2i= t05;/* Swap 6/E outputs for DIT */\
	__BAr= _d ;		__BAi= t29;\
\
	/*...Block 2: t2,10,18,26 */\
	/* Swapped upper right and lower left FMA-pairs of quartet below to hide latency: */\
	_c = t18; _a = t27;\
	 FMA231(t19,__sc,_c );		 FMA231(t26,__sc,_a );		/* twiddle mpy by E^1 */\
	_d = t19; _b = t26;\
	FNMA231(t18,__sc,_d );		 FMS231(t27,__sc,_b );		/* twiddle mpy by E^3 */\
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
	 								_c = t14-t15;	_d = t14+t15;	/* [Block 4] twiddle mpy by -E^6 is here... */\
	/* Write outputs back to main array: */\
	__B1r= t02;		__B1i= t03;\
	__B9r= _a ;		__B9i= _b ;\
	__BDr= t10;		__BDi= t11;/* Swap 5/D outputs for DIT */\
	__B5r= _e ;		__B5i= _f ;\
\
	/*...Block 4: t6,14,22,30 */\
	_e = t23; _f = t22;\
	 FMA231(t22,__sc,_e );		 FMS231(t23,__sc,_f );		/* twiddle mpy by E^3 */\
	_a = t30; _b = t31;\
	 FMA231(t31,__sc,_a );		FNMA231(t30,__sc,_b );		/* twiddle mpy by E^1 =__c*[ -E^9... */\
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
	__B3r= t06;		__B3i= t07;\
	__BBr= _a ;		__BBi= _b ;\
	__BFr= t14;		__BFi= t15;/* Swap 7/F outputs for DIT */\
	__B7r= _c ;		__B7i= _d ;\
}

/************************************/
/*** DIT analogs of above macros: ***/
/************************************/
/* Version #1 uses 4-operand FMA4 syntax, similar to what we might use for SIMD assembler on SSE5-supporting AMD CPUs.

Total: 174 FMA [a whopping 112 of which have a trivial unity multiplicand], 34 MUL.

Extra MUL cost vs DIF is due to the post-twiddling, which is less amenable to "multiplier absorption".
*/
//void ___RADIX_16_DIT_FMA4() {}
#define RADIX_16_DIT_FMA4(\
	__A0r,__A0i,__A1r,__A1i,__A2r,__A2i,__A3r,__A3i,__A4r,__A4i,__A5r,__A5i,__A6r,__A6i,__A7r,__A7i,__A8r,__A8i,__A9r,__A9i,__AAr,__AAi,__ABr,__ABi,__ACr,__ACi,__ADr,__ADi,__AEr,__AEi,__AFr,__AFi,\
	__B0r,__B0i,__B1r,__B1i,__B2r,__B2i,__B3r,__B3i,__B4r,__B4i,__B5r,__B5i,__B6r,__B6i,__B7r,__B7i,__B8r,__B8i,__B9r,__B9i,__BAr,__BAi,__BBr,__BBi,__BCr,__BCi,__BDr,__BDi,__BEr,__BEi,__BFr,__BFi,\
	__c1 ,__s1 ,__c2 ,__s2 ,__c3 ,__s3 ,__c4 ,__s4 ,__c5 ,__s5 ,__c6 ,__s6 ,__c7 ,__s7 ,__c8 ,__s8 ,__c9 ,__s9 ,__cA ,__sA ,__cB ,__sB ,__cC ,__sC ,__cD ,__sD ,__cE ,__sE ,__cF ,__sF ,\
	__c,__s)\
{\
	double t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32;\
	double _r,_i;\
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
	_r =__A1r;	_i =__A1i;\
	t3 =FNMA4(1.0,_r,t1);		t1 = FMA4(1.0,_r,t1);\
	t4 =FNMA4(1.0,_i,t2);		t2 = FMA4(1.0,_i,t2);\
\
	t5 =__A2r;	t6 =__A2i;\
	_r =__A3r;	_i =__A3i;\
	t7 =FNMA4(1.0,_r,t5);		t5 = FMA4(1.0,_r,t5);\
	t8 =FNMA4(1.0,_i,t6);		t6 = FMA4(1.0,_i,t6);\
\
	_r =t5;	t5 =FNMA4(1.0,_r,t1);	t1 = FMA4(1.0,_r,t1);\
	_i =t6;	t6 =FNMA4(1.0,_i,t2);	t2 = FMA4(1.0,_i,t2);\
\
	_i =t7;	t7 =FNMA4(1.0,t8,t3);	t3 = FMA4(1.0,t8,t3);\
			t8 = FMA4(1.0,_i,t4);	t4 =FNMA4(1.0,_i,t4);\
\
	/*...Block 2: */\
	t9 =__A4r;	t10=__A4i;\
	_r =__A5r;	_i =__A5i;\
	t11=FNMA4(1.0,_r,t9 );		t9 = FMA4(1.0,_r,t9 );\
	t12=FNMA4(1.0,_i,t10);		t10= FMA4(1.0,_i,t10);\
\
	t13=__A6r;	t14=__A6i;\
	_r =__A7r;	_i =__A7i;\
	t15=FNMA4(1.0,_r,t13);		t13= FMA4(1.0,_r,t13);\
	t16=FNMA4(1.0,_i,t14);		t14= FMA4(1.0,_i,t14);\
\
	_r =t13;	t13=FNMA4(1.0,_r,t9 );	t9 = FMA4(1.0,_r,t9 );\
	_i =t14;	t14=FNMA4(1.0,_i,t10);	t10= FMA4(1.0,_i,t10);\
\
	_i =t15;	t15=FNMA4(1.0,t16,t11);	t11= FMA4(1.0,t16,t11);\
				t16= FMA4(1.0, _i,t12);	t12=FNMA4(1.0, _i,t12);\
\
	/*...Block 3: */\
	t17=__A8r;	t18=__A8i;\
	_r =__A9r;	_i =__A9i;\
	t19=FNMA4(1.0,_r,t17);		t17= FMA4(1.0,_r,t17);\
	t20=FNMA4(1.0,_i,t18);		t18= FMA4(1.0,_i,t18);\
\
	t21=__AAr;	t22=__AAi;\
	_r =__ABr;	_i =__ABi;\
	t23=FNMA4(1.0,_r,t21);		t21= FMA4(1.0,_r,t21);\
	t24=FNMA4(1.0,_i,t22);		t22= FMA4(1.0,_i,t22);\
\
	_r =t21;	t21=FNMA4(1.0,_r,t17);	t17= FMA4(1.0,_r,t17);\
	_i =t22;	t22=FNMA4(1.0,_i,t18);	t18= FMA4(1.0,_i,t18);\
\
	_i =t23;	t23=FNMA4(1.0,t24,t19);	t19= FMA4(1.0,t24,t19);\
				t24= FMA4(1.0, _i,t20);	t20=FNMA4(1.0, _i,t20);\
\
	/*...Block 4: */\
	t25=__ACr;	t26=__ACi;\
	_r =__ADr;	_i =__ADi;\
	t27=FNMA4(1.0,_r,t25);		t25= FMA4(1.0,_r,t25);\
	t28=FNMA4(1.0,_i,t26);		t26= FMA4(1.0,_i,t26);\
\
	t29=__AEr;	t30=__AEi;\
	_r =__AFr;	_i =__AFi;\
	t31=FNMA4(1.0,_r,t29);		t29= FMA4(1.0,_r,t29);\
	t32=FNMA4(1.0,_i,t30);		t30= FMA4(1.0,_i,t30);\
\
	_r =t29;	t29=FNMA4(1.0,_r,t25);	t25= FMA4(1.0,_r,t25);\
	_i =t30;	t30=FNMA4(1.0,_i,t26);	t26= FMA4(1.0,_i,t26);\
\
	_i =t31;	t31=FNMA4(1.0,t32,t27);	t27= FMA4(1.0,t32,t27);\
				t32= FMA4(1.0, _i,t28);	t28=FNMA4(1.0, _i,t28);\
\
/**************************************************************************************/\
/*...and now do four more radix-4 transforms, including the internal twiddle factors: */\
/**************************************************************************************/\
\
	/*...Block 1: t1,9,17,25 */\
	_r =t9 ;	t9 =FNMA4(1.0,_r,t1 );	t1 = FMA4(1.0,_r,t1 );\
	_i =t10;	t10=FNMA4(1.0,_i,t2 );	t2 = FMA4(1.0,_i,t2 );\
\
	_r =t25;	t25=FNMA4(1.0,_r,t17);	t17= FMA4(1.0,_r,t17);\
	_i =t26;	t26=FNMA4(1.0,_i,t18);	t18= FMA4(1.0,_i,t18);\
\
	__B0r = FMA4(1.0,t17,t1 );		__B0i = FMA4(1.0,t18,t2 );\
	t1    =FNMA4(1.0,t17,t1 );		t2    =FNMA4(1.0,t18,t2 );\
	_r    = FMA4(t2,__t8,t1 );		_i    =FNMA4(t1,__t8,t2 );\
	__B8r = _r *__c8;	__B8i = _i *__c8;/*** Are these pre-output MULs in the DIF/post-twiddle algo avoidable? ***/\
\
	_r    = FMA4(1.0,t26,t9 );		_i	  =FNMA4(1.0,t25,t10);	/* mpy by E^-4 = -I is inlined here... */\
	t9    =FNMA4(1.0,t26,t9 );		t10	  = FMA4(1.0,t25,t10);\
	t25= FMA4(_i ,__t4,_r );		t26=FNMA4(_r ,__t4,_i );\
	_r = FMA4(t10,__tC,t9 );		_i =FNMA4(t9 ,__tC,t10);\
	__B4r = t25*__c4;	__B4i = t26*__c4;\
	__BCr = _r *__cC;	__BCi = _i *__cC;\
\
	/*...Block 3: t5,13,21,29 */\
	_r =t13;	t13=FNMA4(t14,1.0,t5 );	t5 = FMA4(t14,1.0,t5 );		/* twiddle mpy by E^4 =-I */\
				t14= FMA4( _r,1.0,t6 );	t6 =FNMA4( _r,1.0,t6 );\
\
	_r = FMA4(1.0,t21,t22);				t22=FNMA4(1.0,t21,t22);	t21=_r;	/* twiddle mpy by E^-2 */\
	_r = FMS4(1.0,t29,t30);				_i = FMA4(1.0,t29,t30);		/* twiddle mpy by E^2 = -E^-6 is here... */\
	t29= FMA4(1.0,t21,_r );				t21= FMS4(1.0,t21,_r );			/* ...and get E^-6 by flipping signs here. */\
	t30= FMA4(1.0,t22,_i );				t22= FMS4(1.0,t22,_i );\
\
	_r  = FMA4(t21,ISRT2,t5 );			_i = FMA4(t22,ISRT2,t6 );\
	t5  =FNMA4(t21,ISRT2,t5 );			t6 =FNMA4(t22,ISRT2,t6 );\
	__B2r =( FMA4(_i ,__t2,_r ))*__c2;	__B2i =(FNMA4(_r ,__t2,_i ))*__c2;\
	__BAr =( FMA4(t6 ,__tA,t5 ))*__cA;	__BAi =(FNMA4(t5 ,__tA,t6 ))*__cA;\
\
	_r  = FMA4(t30,ISRT2,t13);			_i  =FNMA4(t29,ISRT2,t14);	/* mpy by E^-4 =-I is inlined here... */\
	t13	=FNMA4(t30,ISRT2,t13);			t14 = FMA4(t29,ISRT2,t14);\
	__B6r =( FMA4(_i ,__t6,_r ))*__c6;	__B6i =(FNMA4(_r ,__t6,_i ))*__c6;\
	__BEr =( FMA4(t14,__tE,t13))*__cE;	__BEi =(FNMA4(t13,__tE,t14))*__cE;\
\
	/*...Block 2: t3,11,19,27 */\
	_r = FMA4(1.0,t12,t11);				_i = FMS4(1.0,t12,t11);		/* twiddle mpy by E^-2 */\
	t11=FNMA4(_r ,ISRT2,t3 );			t3 = FMA4(_r ,ISRT2,t3 );\
	t12=FNMA4(_i ,ISRT2,t4 );			t4 = FMA4(_i ,ISRT2,t4 );\
\
	_r = FMA4(t20,__sc,t19);			t20=FNMA4(t19,__sc,t20);	t19=_r;	/* twiddle mpy by E^-1 */\
	_r = FMA4(t27,__sc,t28);			_i = FMS4(t28,__sc,t27);		/* twiddle mpy by E^-3 */\
	_r *= __c;	_i *= __c;	/* Doing things this way incurs just 2 extra MULs, other 2 absorbed into the ensuing 4 FMAs */\
	t27= FMS4(__c,t19,_r );				t19= FMA4(__c,t19,_r );\
	t28= FMS4(__c,t20,_i );				t20= FMA4(__c,t20,_i );\
\
	_r  = FMA4(1.0,t3 ,t19 );			_i  = FMA4(1.0,t4 ,t20 );\
	t3  = FMS4(1.0,t3 ,t19 );			t4  = FMS4(1.0,t4 ,t20 );\
	__B1r =( FMA4(_i ,__t1,_r ))*__c1;	__B1i =(FNMA4(_r ,__t1,_i ))*__c1;\
	__B9r =( FMA4(t4 ,__t9,t3 ))*__c9;	__B9i =(FNMA4(t3 ,__t9,t4 ))*__c9;\
\
	_r  = FMA4(1.0,t11,t28);			_i  = FMS4(1.0,t12,t27);	/* mpy by E^-4 = -I is inlined here... */\
	t11 = FMS4(1.0,t11,t28);			t12 = FMA4(1.0,t12,t27);\
	__B5r =( FMA4(_i ,__t5,_r ))*__c5;	__B5i =(FNMA4(_r ,__t5,_i ))*__c5;\
	__BDr =( FMA4(t12,__tD,t11))*__cD;	__BDi =(FNMA4(t11,__tD,t12))*__cD;\
\
	/*...Block 4: t7,15,23,31 */\
	_r = FMS4(1.0,t15,t16);				_i = FMA4(1.0,t15,t16);		/* twiddle mpy by E^2 = -E^-6 is here... */\
	t15= FMA4(_r ,ISRT2,t7 );			t7 =FNMA4(_r ,ISRT2,t7 );	/* ...and get E^6=(i-1)/sq_r by flipping signs here. */\
	t16= FMA4(_i ,ISRT2,t8 );			t8 =FNMA4(_i ,ISRT2,t8 );\
\
	_r = FMA4(t23,__sc,t24);			t24= FMS4(t24,__sc,t23);	t23=_r;	/* twiddle mpy by E^-3 */\
	_r = FMA4(t32,__sc,t31);			_i =FNMA4(t31,__sc,t32);		/* twiddle mpy by E^-1 = -E^-9... */\
	_r *= __c;	_i *= __c;	/* Doing things this way incurs just 2 extra MULs, other 2 absorbed into the ensuing 4 FMAs */\
	t31= FMA4(__c,t23,_r );				t23= FMS4(__c,t23,_r );			/* ...and get E^9 by flipping signs here. */\
	t32= FMA4(__c,t24,_i );				t24= FMS4(__c,t24,_i );\
\
	_r  = FMA4(1.0,t7 ,t23);			_i  = FMA4(1.0,t8 ,t24);\
	t7  = FMS4(1.0,t7 ,t23);			t8  = FMS4(1.0,t8 ,t24);\
	__B3r =( FMA4(_i ,__t3,_r ))*__c3;	__B3i =(FNMA4(_r ,__t3,_i ))*__c3;\
	__BBr =( FMA4(t8 ,__tB,t7 ))*__cB;	__BBi =(FNMA4(t7 ,__tB,t8 ))*__cB;\
\
	_r  = FMA4(1.0,t15,t32);			_i  = FMS4(1.0,t16,t31);	/* mpy by E^-4 = -I is inlined here... */\
	t15 = FMS4(1.0,t15,t32);			t16 = FMA4(1.0,t16,t31);\
	__B7r =( FMA4(_i ,__t7,_r ))*__c7;	__B7i =(FNMA4(_r ,__t7,_i ))*__c7;\
	__BFr =( FMA4(t16,__tF,t15))*__cF;	__BFi =(FNMA4(t15,__tF,t16))*__cF;\
}

// Version #2 uses a 3-operand FMA3 syntax in which the RIGHTmost of the 3 fused mul-add inputs -
// that is, the addend - is overwritten with the result:
//void ___RADIX_16_DIT_FMA_1() {}
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
/**************************************************************************************/\
/*...and now do four more radix-4 transforms, including the internal twiddle factors: */\
/**************************************************************************************/\
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

// Version #3 is v2, restructured to take account of 16 SIMD registers and 5-cycle FMA latency of Intel Haswell:
/*** NOTE: Due the post-twiddles impl of my DIT FFTs, of the 174 FMAs here, a whopping 112 involve a unity multiplicand,
			i.e. we can - and definitely should attempt-and-compare-timings - replace these by ADD/SUBs.

Total: 174 Fused mul/add. Breakdown, and add/sub instruction each is equivalent to in the one-unity-multiplicand case:
				F*M*231(1.0,b,a)=								F*M*132(1.0,a,b)=
70  FMA231	a = + b.1.0 + a		ADD(b,a,a)			23  FMA132	a = + b.1.0 + a		ADD(b,a,a)
 4  FMS231	a = + b.1.0 - a		SUB(a,b,a)			12  FMS132	a = + b.1.0 - a		SUB(a,b,a)
43 FNMA231	a = - b.1.0 + a		SUB(b,a,a)			22 FNMA132	a = - b.1.0 + a		SUB(b,a,a)

***/

//void ___RADIX_16_DIT_FMA() {}
#define RADIX_16_DIT_FMA(\
	__A0r,__A0i,__A1r,__A1i,__A2r,__A2i,__A3r,__A3i,__A4r,__A4i,__A5r,__A5i,__A6r,__A6i,__A7r,__A7i,__A8r,__A8i,__A9r,__A9i,__AAr,__AAi,__ABr,__ABi,__ACr,__ACi,__ADr,__ADi,__AEr,__AEi,__AFr,__AFi,\
	__B0r,__B0i,__B1r,__B1i,__B2r,__B2i,__B3r,__B3i,__B4r,__B4i,__B5r,__B5i,__B6r,__B6i,__B7r,__B7i,__B8r,__B8i,__B9r,__B9i,__BAr,__BAi,__BBr,__BBi,__BCr,__BCi,__BDr,__BDi,__BEr,__BEi,__BFr,__BFi,\
	/* Sine terms here replaced by (and assumed to hold) tangents: */\
	__c1 ,__t1 ,__c2 ,__t2 ,__c3 ,__t3 ,__c4 ,__t4 ,__c5 ,__t5 ,__c6 ,__t6 ,__c7 ,__t7 ,__c8 ,__t8 ,__c9 ,__t9 ,__cA ,__tA ,__cB ,__tB ,__cC ,__tC ,__cD ,__tD ,__cE ,__tE ,__cF ,__tF ,\
	__c,__sc)\
{\
	/* Note no cosine ratios needed in DIT version: */\
	double t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32;\
	double rt,it,re,im;\
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
/**************************************************************************************/\
/*...and now do four more radix-4 transforms, including the internal twiddle factors: */\
/**************************************************************************************/\
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
/*	 FMA231(ONE,t23,t31);	 FMS132(ONE,t23,re );	 FMA231(ONE,t24,t32);	 FMS132(ONE,t24,im );** B */\
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

// Inlining these larger macros leads to very long compiles, so instead prototype here and define them as functions in dft_macro.c.

/* Totals: 376 ADD, 88 MUL	*/
/* Because of the way the original (non-macro-ized) code was written, it's convenient to index the A-inputs in terms
of 4 length-8 blocks with octal indices, and the B-outputs in terms of 2 length-16 blocks with hexadecimal indices.
MSVC allows a maximum of 'only' 127 macro args (which is probably a good thing), so unlike the smaller-radix DFT
macros which use actual array-indexed terms as args, here we use pointers to the real part of each complex arg:
*/
void RADIX_32_DIF(
	double *__A, const int *__idx, const int __re_im_stride_in,	/*  Inputs: Base address plus 32 (index) offsets */
	double *__B, const int *__odx, const int __re_im_stride_out	/* Outputs: Base address plus 32 (index) offsets */
);

// With-twiddles out-of-place analog of above twiddleless DIF macro: 31 nontrivial complex input twiddles E01-1f [E0 assumed = 1],
// The DIF version of this macro processes the twiddles in bit-reversed order: 0,8,4,c,2,a,6,e,1,9,5,d,3,b,7,f.
void RADIX_32_DIF_TWIDDLE_OOP(
	double *__A, const int *__idx,	/*  Inputs: Base address plus 32 (index) offsets */
	double *__B, const int *__odx,	/* Outputs: Base address plus 32 (index) offsets */
	/* Twiddles: */
	double __c10,double __s10,
	double __c08,double __s08,
	double __c18,double __s18,
	/**/
	double __c04,double __s04,
	double __c14,double __s14,
	double __c0C,double __s0C,
	double __c1C,double __s1C,
	/**/
	double __c02,double __s02,
	double __c12,double __s12,
	double __c0A,double __s0A,
	double __c1A,double __s1A,
	/**/
	double __c06,double __s06,
	double __c16,double __s16,
	double __c0E,double __s0E,
	double __c1E,double __s1E,
	/**/
	double __c01,double __s01,
	double __c11,double __s11,
	double __c09,double __s09,
	double __c19,double __s19,
	/**/
	double __c05,double __s05,
	double __c15,double __s15,
	double __c0D,double __s0D,
	double __c1D,double __s1D,
	/**/
	double __c03,double __s03,
	double __c13,double __s13,
	double __c0B,double __s0B,
	double __c1B,double __s1B,
	/**/
	double __c07,double __s07,
	double __c17,double __s17,
	double __c0F,double __s0F,
	double __c1F,double __s1F
);

void RADIX_32_DIT(
	double *__A, const int *__idx, const int __re_im_stride_in,	/*  Inputs: Base address plus 32 (index) offsets */
	double *__B, const int *__odx, const int __re_im_stride_out	/* Outputs: Base address plus 32 (index) offsets */
);

// Post-twiddles out-of-place-permitting analog of above twiddleless DIF macro: 31 nontrivial complex input twiddles E01-1f [E0 assumed = 1],
// The DIF version of this macro processes the twiddles in bit-reversed order: 0,8,4,c,2,a,6,e,1,9,5,d,3,b,7,f.
void RADIX_32_DIT_TWIDDLE(
	double *__A, const int *__idx,	/*  Inputs: Base address plus 32 (index) offsets */
	double *__B, const int *__odx,	/* Outputs: Base address plus 32 (index) offsets */
	/* Twiddles: */
	double __c10,double __s10,
	double __c08,double __s08,
	double __c18,double __s18,
	/**/
	double __c04,double __s04,
	double __c14,double __s14,
	double __c0C,double __s0C,
	double __c1C,double __s1C,
	/**/
	double __c02,double __s02,
	double __c12,double __s12,
	double __c0A,double __s0A,
	double __c1A,double __s1A,
	/**/
	double __c06,double __s06,
	double __c16,double __s16,
	double __c0E,double __s0E,
	double __c1E,double __s1E,
	/**/
	double __c01,double __s01,
	double __c11,double __s11,
	double __c09,double __s09,
	double __c19,double __s19,
	/**/
	double __c05,double __s05,
	double __c15,double __s15,
	double __c0D,double __s0D,
	double __c1D,double __s1D,
	/**/
	double __c03,double __s03,
	double __c13,double __s13,
	double __c0B,double __s0B,
	double __c1B,double __s1B,
	/**/
	double __c07,double __s07,
	double __c17,double __s17,
	double __c0F,double __s0F,
	double __c1F,double __s1F
);

// With-twiddles out-of-place analog of above twiddleless DIT macro: 31 nontrivial complex input twiddles E01-1f [E0 assumed = 1],
// The DIT version of this macro processes the twiddles in order.
void RADIX_32_DIT_TWIDDLE_OOP(
	double *__A, const int *__idx,	/*  Inputs: Base address plus 32 (index) offsets */
	double *__B, const int *__odx,	/* Outputs: Base address plus 32 (index) offsets */
	/* Twiddles: */
	double __c01,double __s01,
	double __c02,double __s02,
	double __c03,double __s03,
	double __c04,double __s04,
	double __c05,double __s05,
	double __c06,double __s06,
	double __c07,double __s07,
	double __c08,double __s08,
	double __c09,double __s09,
	double __c0A,double __s0A,
	double __c0B,double __s0B,
	double __c0C,double __s0C,
	double __c0D,double __s0D,
	double __c0E,double __s0E,
	double __c0F,double __s0F,
	double __c10,double __s10,
	double __c11,double __s11,
	double __c12,double __s12,
	double __c13,double __s13,
	double __c14,double __s14,
	double __c15,double __s15,
	double __c16,double __s16,
	double __c17,double __s17,
	double __c18,double __s18,
	double __c19,double __s19,
	double __c1A,double __s1A,
	double __c1B,double __s1B,
	double __c1C,double __s1C,
	double __c1D,double __s1D,
	double __c1E,double __s1E,
	double __c1F,double __s1F
);

// Twiddleless Radix-32 DIF subtransform macro for use by larger radix-32*k twiddleless-DFT macros.
// OOP = out of place, i.e. Assumes output locs != input locs.
// For the scalar-data versino of this macro it doesn't fundamentally matter whether the in/outputs are local
// scalars or array locations, but our implementation of this macro is driven by our twiddleless-DFT scheme
// in which the DIF step has the power-of-2 component macros following a set of odd-radix ones, thus with:
//
//	o __A-inputs read from a block of contiguous local-allocated storage with unit index stride;
//	o __B-outputs written to an array with arbitrary index strides encoded in the __odx auxiliary array.
//
void RADIX_32_DIF_OOP(
	double *__A,				/*  Inputs: Base address plus 32 (index) offsets */
	double *__B, const int *__odx	/* Outputs: Base address plus 32 (index) offsets */
);

// Twiddleless Radix-32 DIT subtransform macro for use by larger radix-8*k macros.
// OOP = out of place, i.e. Assumes output locs != input locs.
// For the scalar-data versino of this macro it doesn't matter whether the in/outputs are local
// scalars or array locations, but the __A-naming indicates that this macro is intended to serve
// as a prototype for a SIMD/ASM macro in which the __A-inputs are read from an array with arbitrary
// index stride and the __t-outputs go into a block of contiguous local-allocated storage:
void RADIX_32_DIT_OOP(
	double *__A, const int *__idx,/*  Inputs: Base address plus 32 (index) offsets */
	double *__B				/* Outputs: Base address plus 32 (index) offsets */
);

void RADIX_63_DIF(
	double *__A, const int *__idx, const int __re_im_stride_in,	/*  Inputs: Base address plus 63 (index) offsets */
	double *__B, const int *__odx, const int __re_im_stride_out	/* Outputs: Base address plus 63 (index) offsets */
);

void RADIX_63_DIT(
	double *__A, const int *__idx, const int __re_im_stride_in,	/*  Inputs: Base address plus 63 (index) offsets */
	double *__B, const int *__odx, const int __re_im_stride_out	/* Outputs: Base address plus 63 (index) offsets */
);

// For power-of-2 radix > 32 we use a common macro for both the first-pass-of-2-pass-pow2
// and twiddleless-subtransform-to-be-combined-with-odd-radix cases:

void RADIX_64_DIF(
	double *__A, const int *__idx, const int __re_im_stride_in,	/*  Inputs: Base address plus 64 (index) offsets */
	double *__B, const int *__odx, const int __re_im_stride_out	/* Outputs: Base address plus 64 (index) offsets */
);

void RADIX_64_DIT(
	double *__A, const int *__idx, const int __re_im_stride_in,	/*  Inputs: Base address plus 64 (index) offsets */
	double *__B, const int *__odx, const int __re_im_stride_out	/* Outputs: Base address plus 64 (index) offsets */
);

void RADIX_128_DIF(
	double *__A, const int *__idx, const int __re_im_stride_in,	/*  Inputs: Base address plus 128 (index) offsets */
	double *__B, const int *__odx, const int __re_im_stride_out	/* Outputs: Base address plus 128 (index) offsets */
);

void RADIX_128_DIT(
	double *__A, const int *__idx, const int __re_im_stride_in,	/*  Inputs: Base address plus 128 (index) offsets */
	double *__B, const int *__odx, const int __re_im_stride_out	/* Outputs: Base address plus 128 (index) offsets */
);

void RADIX_256_DIF(
	// inputs: Base address plus 32 index offsets:
	double *__A,
	// Since in-array may either coincide with oarray [if using it for a complete radix-256 pass]
	// or a local-complex-data scratch array [if using for phase of a larger-radix DFT algo] and we need
	// this routine to work correctly with both scalar-complex and SIMD main-array data layouts, need a param
	// telling routine whether to use scalar-complex re/im stride of 1 or SIMD [2,4,8... depending on type]:
	const int __re_im_stride_in,
	int *i_offsets_lo,	// Array storing  low parts of input index offsets in 16 slots
	int *i_offsets_hi,	// Array storing high parts of input index offsets in 16 slots
	// outputs: Base address plus index offsets:
	double *__B,
	const int __re_im_stride_out,	// Similar as for in-array
	int *o_offsets_lo,	// Array storing  low parts of output index offsets in 16 slots
	uint32 o_idx,	// Bitfield encoding the sequence of the o_offsets_lo sub-vectors to use for the radix-256 DFT's outputs
	int *o_offsets_hi	// Array storing high parts of output index offsets in 16 slots
);

void RADIX_256_DIT(
	// inputs: Base address plus 32 index offsets:
	double *__A,
	const int __re_im_stride_in,
	int *i_offsets_lo,	// Array storing  low parts of input index offsets in 16 slots
	uint32 i_idx,	// Bitfield encoding the sequence of the i_offsets_lo sub-vectors to use for the radix-256 DFT's inputs
	int *i_offsets_hi,	// Array storing high parts of input index offsets in 16 slots
	// outputs: Base address plus index offsets:
	double *__B,
	const int __re_im_stride_out,
	int *o_offsets_lo,	// Array storing  low parts of output index offsets in 16 slots
	int *o_offsets_hi	// Array storing high parts of output index offsets in 16 slots
);

#ifdef USE_SSE2

void SSE2_RADIX_63_DIF(
	const int init,	// Init consts (in 1-thread mode only!) and exit
	const int thr_id,
	vec_dbl *__A, const int *__idx,	/* Inputs : Base address plus 63 (index) offsets */
	vec_dbl *__B, const int *__odx	/* Outputs: Base address plus 63 (index) offsets */
);

void SSE2_RADIX_63_DIT(
	const int init,	// Init consts (in 1-thread mode only!) and exit
	const int thr_id,
	vec_dbl *__A, const int *__idx,	/* Inputs : Base address plus 63 (index) offsets */
	vec_dbl *__B, const int *__odx	/* Outputs: Base address plus 63 (index) offsets */
);

void SSE2_RADIX_64_DIF(
	const int init,	// Init consts (in 1-thread mode only!) and exit
	const int thr_id,
	// Need to know if the DIF-64 is standalone within a contig block of 64 vec_cmplx data
	// (as occurs if it's a standalone or part of a larger transform of length N = odd*64), or
	// part of a larger power-of-2 transform of length N = 2^k > 64, in which we need to adjust data strides:
	const int pow2_stride_shift,	// set = trailz(N) - trailz(64)
	// Inputs: Base address plus index offsets:
	double *__A, const int *i_offsets,
	// Intermediates-storage pointer:
	vec_dbl*r00,
	// Outputs: Base address plus index offsets:
	double *__B, const int *o_offsets
);

void SSE2_RADIX_64_DIT(
	const int init,	// Init consts (in 1-thread mode only!) and exit
	const int thr_id,
	// Inputs: Base address plus index offsets:
	double *__A, const int *i_offsets,
	// Intermediates-storage pointer:
	vec_dbl*r00,
	// Outputs: Base address plus index offsets:
	vec_dbl*__B, const int *o_offsets
);

void SSE2_RADIX128_DIF(
	// Input pointer:
	vec_dbl*__A,
	// Intermediates-storage pointer:
	vec_dbl*tmp,
	// Pointers to base-roots data and first of 16 twiddle vectors:
	vec_dbl*isrt2, vec_dbl*two, vec_dbl*twid0,
	// Outputs: Base address plus index offsets:
	double *__B,
	int *o_offsets	// Array storing output index offsets
);

void SSE2_RADIX128_DIT(
	// Inputs: Base address plus index offsets:
	double *__A,
	int *i_offsets,	// Array storing input index offsets
	// Intermediates-storage pointer:
	vec_dbl*tmp,
	// Pointers to base-roots data and first of 16 twiddle vectors:
	vec_dbl*isrt2, vec_dbl*two, vec_dbl*twid0,
	// Output pointer: Base ptr of 16 local-mem:
	vec_dbl*__B
);

void SSE2_RADIX256_DIF(
	// Input pointer:
	vec_dbl*__A,
	// Intermediates-storage pointer:
	vec_dbl*tmp,
	// Pointers to base-roots data and first of 16 twiddle vectors:
	vec_dbl*isrt2, vec_dbl*two, vec_dbl*twid0,
	// Outputs: Base address plus index offsets:
	double *__B,
	int *o_offsets_lo,	// Array storing  low parts of output index offsets in 16 slots
	uint32 o_idx,	// Bitfield encoding the sequence of the o_offsets_lo sub-vectors to use for the radix-256 DFT's outputs
	int *o_offsets_hi	// Array storing high parts of output index offsets in 16 slots
);

void SSE2_RADIX256_DIT(
	// Inputs: Base address plus index offsets:
	double *__A,
	int *i_offsets_lo,	// Array storing  low parts of input index offsets in 16 slots
	uint32 i_idx,	// Bitfield encoding the sequence of the i_offsets_lo sub-vectors to use for the radix-256 DFT's inputs
	int *i_offsets_hi,	// Array storing high parts of input index offsets in 16 slots
	// Intermediates-storage pointer:
	vec_dbl*tmp,
	// Pointers to base-roots data and first of 16 twiddle vectors:
	vec_dbl*isrt2, vec_dbl*two, vec_dbl*twid0,
	// Output pointer: Base ptr of 16 local-mem:
	vec_dbl*__B
);

#endif

#endif	/* #ifndef dft_macro_included */
