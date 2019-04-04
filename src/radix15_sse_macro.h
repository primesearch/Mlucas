/*******************************************************************************
*                                                                              *
*   (C) 1997-2017 by Ernst W. Mayer.                                           *
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
#ifndef radix15_sse_macro_h_included
#define radix15_sse_macro_h_included

#include "sse2_macro.h"

/* General indexing for twiddleless radix-15 done as 3*radix-5 followed by 5*radix-3 is as for the scalar macro above:
RADIX_15_DIF(00,01,02,03,04,05,06,07,08,09,0A,0B,0C,0D,0E)
->
	RADIX_05_DFT(i0,iC,i9,i6,i3, t0,t1,t2,t3,t4)
	RADIX_05_DFT(iA,i7,i4,i1,iD, t5,t6,t7,t8,t9)
	RADIX_05_DFT(i5,i2,iE,iB,i8, tA,tB,tC,tD,tE)

	RADIX_03_DFT(t0,t5,tA, o0,o1,o2,)
	RADIX_03_DFT(t1,t6,tB, oD,oE,oB,)
	RADIX_03_DFT(t2,t7,tC, o9,oA,oB,)
	RADIX_03_DFT(t3,t8,tD, o8,o6,o7,)
	RADIX_03_DFT(t4,t9,tE, o4,o5,o3,)

In our impl below, the __i are input pointers, which may overlap the __o outputs;
..cc0 and cc1 are ptrs to the radix-3 and radix-5 SSE2 sincos constants (c3m1 and cn1);
__t0-E are ptr to scratch local storage (i.e. the address block pointed to by r00-r3e).
*/
// Aug 2014: Need arbitrary-pointer-offsets to support I/O permutations needed by
// larger-radix DFTs of length 15 * 2^n

#define SSE2_RADIX_15_DIF(\
	__cc0, __cc1,\
	__i0,__i1,__i2,__i3,__i4,__i5,__i6,__i7,__i8,__i9,__iA,__iB,__iC,__iD,__iE,\
	__t0,__t1,__t2,__t3,__t4,__t5,__t6,__t7,__t8,__t9,__tA,__tB,__tC,__tD,__tE,\
	__o0,__o1,__o2,__o3,__o4,__o5,__o6,__o7,__o8,__o9,__oA,__oB,__oC,__oD,__oE)\
{\
	SSE2_RADIX_05_DFT_0TWIDDLE(__i0,__iC,__i9,__i6,__i3, __cc1, __t0,__t1,__t2,__t3,__t4);\
	SSE2_RADIX_05_DFT_0TWIDDLE(__iA,__i7,__i4,__i1,__iD, __cc1, __t5,__t6,__t7,__t8,__t9);\
	SSE2_RADIX_05_DFT_0TWIDDLE(__i5,__i2,__iE,__iB,__i8, __cc1, __tA,__tB,__tC,__tD,__tE);\
\
	SSE2_RADIX_03_DFT(__t0,__t5,__tA, __cc0, __o0,__o1,__o2);\
	SSE2_RADIX_03_DFT(__t1,__t6,__tB, __cc0, __oD,__oE,__oC);\
	SSE2_RADIX_03_DFT(__t2,__t7,__tC, __cc0, __o9,__oA,__oB);\
	SSE2_RADIX_03_DFT(__t3,__t8,__tD, __cc0, __o8,__o6,__o7);\
	SSE2_RADIX_03_DFT(__t4,__t9,__tE, __cc0, __o4,__o5,__o3);\
}

#define SSE2_RADIX_15_DIT(\
	__cc0, __cc1,\
	__i0,__i1,__i2,__i3,__i4,__i5,__i6,__i7,__i8,__i9,__iA,__iB,__iC,__iD,__iE,\
	__t0,__t1,__t2,__t3,__t4,__t5,__t6,__t7,__t8,__t9,__tA,__tB,__tC,__tD,__tE,\
	__o0,__o1,__o2,__o3,__o4,__o5,__o6,__o7,__o8,__o9,__oA,__oB,__oC,__oD,__oE)\
{\
/* Swap the 2nd pair of each output triplet to effect iDFT: */\
	SSE2_RADIX_03_DFT(__i0,__i2,__i1, __cc0, __t0,__t2,__t1);\
	SSE2_RADIX_03_DFT(__i8,__i7,__i6, __cc0, __t3,__t5,__t4);\
	SSE2_RADIX_03_DFT(__iD,__iC,__iE, __cc0, __t6,__t8,__t7);\
	SSE2_RADIX_03_DFT(__i4,__i3,__i5, __cc0, __t9,__tB,__tA);\
	SSE2_RADIX_03_DFT(__i9,__iB,__iA, __cc0, __tC,__tE,__tD);\
\
/* Output perm here is 0123456789abcde --> 05a6b1c2738d9e4: */\
	SSE2_RADIX_05_DFT_0TWIDDLE(__t0,__t3,__t6,__t9,__tC, __cc1, __o0,__o6,__oC,__o3,__o9);\
	SSE2_RADIX_05_DFT_0TWIDDLE(__t1,__t4,__t7,__tA,__tD, __cc1, __o5,__oB,__o2,__o8,__oE);\
	SSE2_RADIX_05_DFT_0TWIDDLE(__t2,__t5,__t8,__tB,__tE, __cc1, __oA,__o1,__o7,__oD,__o4);\
}

// Cost: 12 DP-math, 17 vector MOV for each of the two side-by-side 3-DFTs in SSE2_RADIX_03_DFT_X2
//       38 DP-math, 31 vector MOV for each of the two side-by-side 5-DFTs in SSE2_RADIX_05_DFT_0TWIDDLE_X2. Thus
// 150 DP-math, 144 vector MOV for each of the two side-by-side 15-DFTs in each of these two [DIF and DIT] 15-DFT macro-of-macros.
// Compare to van-Buskirk 13-DFT: 198 DP-math, 168 vector MOV.
#define SSE2_RADIX_15_DIF_X2(\
	__cc0, __cc1, __two,\
	__i0,__i1,__i2,__i3,__i4,__i5,__i6,__i7,__i8,__i9,__iA,__iB,__iC,__iD,__iE,\
	__s0,__s1,__s2,__s3,__s4,__s5,__s6,__s7,__s8,__s9,__sA,__sB,__sC,__sD,__sE,\
	__o0,__o1,__o2,__o3,__o4,__o5,__o6,__o7,__o8,__o9,__oA,__oB,__oC,__oD,__oE,\
			__j0,__j1,__j2,__j3,__j4,__j5,__j6,__j7,__j8,__j9,__jA,__jB,__jC,__jD,__jE,\
			__t0,__t1,__t2,__t3,__t4,__t5,__t6,__t7,__t8,__t9,__tA,__tB,__tC,__tD,__tE,\
			__u0,__u1,__u2,__u3,__u4,__u5,__u6,__u7,__u8,__u9,__uA,__uB,__uC,__uD,__uE)\
{\
	SSE2_RADIX_05_DFT_0TWIDDLE_X2(__cc1,__two, __i0,__iC,__i9,__i6,__i3, __s0,__s1,__s2,__s3,__s4,	__j0,__jC,__j9,__j6,__j3, __t0,__t1,__t2,__t3,__t4);\
	SSE2_RADIX_05_DFT_0TWIDDLE_X2(__cc1,__two, __iA,__i7,__i4,__i1,__iD, __s5,__s6,__s7,__s8,__s9,	__jA,__j7,__j4,__j1,__jD, __t5,__t6,__t7,__t8,__t9);\
	SSE2_RADIX_05_DFT_0TWIDDLE_X2(__cc1,__two, __i5,__i2,__iE,__iB,__i8, __sA,__sB,__sC,__sD,__sE,	__j5,__j2,__jE,__jB,__j8, __tA,__tB,__tC,__tD,__tE);\
\
	SSE2_RADIX_03_DFT_X2(__cc0, __s0,__s5,__sA, __o0,__o1,__o2,		__t0,__t5,__tA, __u0,__u1,__u2);\
	SSE2_RADIX_03_DFT_X2(__cc0, __s1,__s6,__sB, __oD,__oE,__oC,		__t1,__t6,__tB, __uD,__uE,__uC);\
	SSE2_RADIX_03_DFT_X2(__cc0, __s2,__s7,__sC, __o9,__oA,__oB,		__t2,__t7,__tC, __u9,__uA,__uB);\
	SSE2_RADIX_03_DFT_X2(__cc0, __s3,__s8,__sD, __o8,__o6,__o7,		__t3,__t8,__tD, __u8,__u6,__u7);\
	SSE2_RADIX_03_DFT_X2(__cc0, __s4,__s9,__sE, __o4,__o5,__o3,		__t4,__t9,__tE, __u4,__u5,__u3);\
}

#define SSE2_RADIX_15_DIT_X2(\
	__cc0, __cc1,__two,\
	__i0,__i1,__i2,__i3,__i4,__i5,__i6,__i7,__i8,__i9,__iA,__iB,__iC,__iD,__iE,\
	__s0,__s1,__s2,__s3,__s4,__s5,__s6,__s7,__s8,__s9,__sA,__sB,__sC,__sD,__sE,\
	__o0,__o1,__o2,__o3,__o4,__o5,__o6,__o7,__o8,__o9,__oA,__oB,__oC,__oD,__oE,\
			__j0,__j1,__j2,__j3,__j4,__j5,__j6,__j7,__j8,__j9,__jA,__jB,__jC,__jD,__jE,\
			__t0,__t1,__t2,__t3,__t4,__t5,__t6,__t7,__t8,__t9,__tA,__tB,__tC,__tD,__tE,\
			__u0,__u1,__u2,__u3,__u4,__u5,__u6,__u7,__u8,__u9,__uA,__uB,__uC,__uD,__uE)\
{\
/* Swap the 2nd pair of each output triplet to effect iDFT: */\
	SSE2_RADIX_03_DFT_X2(__cc0, __i0,__i2,__i1, __s0,__s2,__s1,		__j0,__j2,__j1, __t0,__t2,__t1);\
	SSE2_RADIX_03_DFT_X2(__cc0, __i8,__i7,__i6, __s3,__s5,__s4,		__j8,__j7,__j6, __t3,__t5,__t4);\
	SSE2_RADIX_03_DFT_X2(__cc0, __iD,__iC,__iE, __s6,__s8,__s7,		__jD,__jC,__jE, __t6,__t8,__t7);\
	SSE2_RADIX_03_DFT_X2(__cc0, __i4,__i3,__i5, __s9,__sB,__sA,		__j4,__j3,__j5, __t9,__tB,__tA);\
	SSE2_RADIX_03_DFT_X2(__cc0, __i9,__iB,__iA, __sC,__sE,__sD,		__j9,__jB,__jA, __tC,__tE,__tD);\
\
/* Output perm here is 0123456789abcde --> 05a6b1c2738d9e4: */\
	SSE2_RADIX_05_DFT_0TWIDDLE_X2(__cc1,__two, __s0,__s3,__s6,__s9,__sC, __o0,__o6,__oC,__o3,__o9,	__t0,__t3,__t6,__t9,__tC, __u0,__u6,__uC,__u3,__u9);\
	SSE2_RADIX_05_DFT_0TWIDDLE_X2(__cc1,__two, __s1,__s4,__s7,__sA,__sD, __o5,__oB,__o2,__o8,__oE,	__t1,__t4,__t7,__tA,__tD, __u5,__uB,__u2,__u8,__uE);\
	SSE2_RADIX_05_DFT_0TWIDDLE_X2(__cc1,__two, __s2,__s5,__s8,__sB,__sE, __oA,__o1,__o7,__oD,__o4,	__t2,__t5,__t8,__tB,__tE, __uA,__u1,__u7,__uD,__u4);\
}

#endif	/* radix15_sse_macro_h_included */

