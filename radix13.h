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
#ifndef radix13_included
#define radix13_included

#define DC1 ((double)-1.66277335223429382656331008444690382)
#define DC2 ((double) 0.73124599097534822519618254560377760)
#define DC3 ((double) 1.0070740657275332544937477077369340)
#define DC4 ((double)-0.30816846519175820059367219184688543)
#define DC5 ((double) 0.81698338691215549726750306085822509)
#define DC6 ((double) 0.22376403323791637458136993583748652)
#define DS1 ((double) 0.57514072947400312136838554745545335)
#define DS2 ((double)-0.17413860115213590500566079492926474)
#define DS3 ((double)-0.33582506518644535421963182119524142)
#define DS4 ((double) 4.5240494294812713569277280991401412E-0002)
#define DS5 ((double) 1.1543953381323634420147226757584967)
#define DS6 ((double)-8.7981928766792081008399945211312915E-0002)
#define DS7 ((double) 0.90655220171271016880349079977456841)
#define DS8 ((double) 1.1971367726043428094538453399784083)
#define DS9 ((double)-0.24784313641965327321123187598392850)
#define DSa ((double)-0.86131170741789745523421351878316690)
#define DSb ((double)-4.2741434471979367439122664219911502E-0002)

/*
The basic tangent-based radix-13 DFT macro is defined in radix13.h::RADIX_13_DFT(...).
Here is an 8-register-optimized version which moreover mimics x86-style destructive (2-input, one overwritten with output)
register arithmetic. This version uses the copies-into-the-__t-temps to allow for an in-place (inputs overlap outputs) calling
convention; if only an out-of-place version is needed (as in the carry routine proper, delete the __t*r = __A*r assignments at
the beginning, and replace the __t*r right-hand-side operands in the opening the yi-terms of the imaginary parts with __A*r.

This version increases the number of multiplies somewhat relative to the original temporary-heavy implemetation (from 68 to 88)
since it uses muls-by-2 to effect the in-place doubling needed by the very frequent radix-2 butterfly operation sequence which
takes inputs x,y and outputs x+y and x-y:

	x -= y
	y *= 2
	y += x

These extra constants are needed for this version:
	DC23 =  DC2/DC3	(Note that DC6/DC4 = -DC2/DC3)
	DC54 =  DC5/DC4
	DC65 =  DC6/DC5
	DS63 = -DS6/DS3
	DS74 = -DS7/DS4
	DS85 = -DS8/DS5
	DS93 = -DS9/DS3
	DSa4 = -DSa/DS4
	DSb5 = -DSb/DS5
Note we do not need the above consts DC2,DC5,DC6,DS6,DS7,DS8,DS9,DSa,DSb, i.e. this versions needs the same number of precomputed constants, 17:
*/
#define DC23 ((double) 0.7261094450357824054685101554)
#define DC54 ((double) -2.651093408937175306253240338)
#define DC65 ((double)  0.273890554964217594531489845)
#define DS63 ((double) -0.2619873794052440703141563891)
#define DS74 ((double)-20.03851230725071170999037075)
#define DS85 ((double) -1.037024954155764805303318302)
#define DS93 ((double) -0.7380126205947559296858436109)
#define DSa4 ((double) 19.03851230725071170999037075)
#define DSb5 ((double)  0.03702495415576480530331830225)

#define RADIX_13_XYZ(\
	__A0r,__A0i,__A1r,__A1i,__A2r,__A2i,__A3r,__A3i,__A4r,__A4i,__A5r,__A5i,__A6r,__A6i,__A7r,__A7i,__A8r,__A8i,__A9r,__A9i,__Aar,__Aai,__Abr,__Abi,__Acr,__Aci,\
	__B0r,__B0i,__B1r,__B1i,__B2r,__B2i,__B3r,__B3i,__B4r,__B4i,__B5r,__B5i,__B6r,__B6i,__B7r,__B7i,__B8r,__B8i,__B9r,__B9i,__Bar,__Bai,__Bbr,__Bbi,__Bcr,__Bci\
)\
{\
	double xr0, xr1, xr2, xr3, xr4, xr5, xr6, xr7;\
	double __t1r,__t2r,__t3r,__t4r,__t5r,__t6r,__t7r,__t8r,__t9r,__tar,__tbr,__tcr;\
		__t1r = __A1r;\
		__t2r = __A2r;\
		__t3r = __A3r;\
		__t4r = __A4r;\
		__t5r = __A5r;\
		__t6r = __A6r;\
		__t7r = __A7r;\
		__t8r = __A8r;\
		__t9r = __A9r;\
		__tar = __Aar;\
		__tbr = __Abr;\
		__tcr = __Acr;\
	/***************/\
	/* REAL PARTS: */\
	/***************/\
	/* xr-terms need 8 registers for each side: */\
		xr7 = __A6r;\
		xr5 = __A5r;\
		xr4 = __A4r;\
		xr6 = __A3r;\
		xr3 = __A2r;\
		xr1 = __A1r;\
		xr7 += __A7r;\
		xr5 += __A8r;\
		xr4 += __A9r;\
		xr6 += __Aar;\
		xr3 += __Abr;\
		xr1 += __Acr;\
		xr1 -= xr5;	xr5 *= 2.0;\
		xr3 -= xr6;	xr6 *= 2.0;\
		xr4 -= xr7;	xr7 *= 2.0;\
		xr5 += xr1;\
		xr6 += xr3;\
		xr7 += xr4;\
		xr2 = xr5;\
		xr2 += xr6;\
		xr2 += xr7;\
		xr0 = __A0r;\
		xr0 += xr2;\
		xr2 *= DC1;\
		__B0r = xr0;\
		xr0 += xr2;\
		xr2 = DC3;\
		xr6 *= xr2;\
		xr5 *= xr2;\
		xr7 *= xr2;\
		__B1r = xr6;\
		__B2r = xr5;\
		__B3r = xr7;\
		xr2 = DC23;\
		xr5 *= xr2;\
		xr7 *= xr2;\
		xr6 *= xr2;\
		xr5 += xr0;\
		xr7 += xr0;\
		xr6 += xr0;\
		xr5 += __B1r;\
		xr7 += __B2r;\
		xr6 += __B3r;\
		/* Rearrange the 3x3 multiply block to ease in-place computation: */\
		xr2 = DC65;\
		xr0 = DC54;\
		__B1r = xr4;\
		__B2r = xr1;\
		__B3r = xr3;\
		xr4 *= xr2;  /* *= DC65 */\
		xr1 *= xr2;  /* *= DC65 */\
		xr3 *= xr2;  /* *= DC65 */\
		xr4 += __B3r;\
		xr1 -= __B1r;\
		xr3 += __B2r;\
		xr4 *= xr0;  /* *= DC54 */\
		xr1 *= xr0;  /* *= DC54 */\
		xr3 *= xr0;  /* *= DC54 */\
		xr2 = DC4;\
		xr4 += __B2r;\
		xr1 -= __B3r;\
		xr3 -= __B1r;\
		xr4 *= xr2;  /* *= DC4 */\
		xr1 *= xr2;  /* *= DC4 */\
		xr3 *= xr2;  /* *= DC4 */\
	/* Spill into destination outputs: */\
		xr6 -= xr4;	xr4 *= 2.0;\
		xr7 -= xr1;	xr1 *= 2.0;\
		xr5 -= xr3;	xr3 *= 2.0;\
		__B3r = xr7;\
		__B8r = xr5;\
		xr4 += xr6;\
		xr7 += xr1;\
		xr5 += xr3;\
		__B1r = xr5;\
		__B2r = xr7;\
		__B4r = xr6;\
		__B6r = xr4;\
	/* yi-terms: */\
		xr1 = __A1i;\
		xr2 = __A2i;\
		xr5 = __A3i;\
		xr3 = __A4i;\
		xr4 = __A5i;\
		xr6 = __A6i;\
		xr1 -= __Aci;\
		xr2 -= __Abi;\
		xr5 -= __Aai;\
		xr3 -= __A9i;\
		xr4 -= __A8i;\
		xr6 -= __A7i;\
		xr7 = xr1;\
		xr0 = xr2;\
		xr7 -= xr3;\
		xr0 += xr4;\
		xr7 += xr5;\
		xr0 += xr6;\
		xr1 += xr3;\
		xr5 += xr3;\
		xr2 -= xr4;\
		xr6 -= xr4;\
		xr4 = xr0;\
		xr3 = xr7;\
		xr0 *= DS2;\
		xr7 *= DS2;\
		xr4 *= DS1;\
		xr3 *= DS1;\
		xr4 -= xr7;\
		xr3 += xr0;\
		__Bcr = xr4;	/* tmp-store in Bcr */\
		xr0 = xr1;\
		xr4 = xr5;\
		xr0 += xr2;\
		xr4 += xr6;\
		/*\
		xr7 = DS3*xr0-DS6*xr4;\
		xr0 = DS3*xr4-DS9*xr0;\
		*/\
		__Bbr = xr0;	/* tmp-store in Bbr */\
		  xr7 = xr4;\
		xr4 *= DS63;\
		xr0 *= DS93;\
		xr4 += __Bbr;\
		xr0 +=   xr7;\
		xr4 *= DS3;\
		xr0 *= DS3;\
		/*\
		xr7 = xr4+DS4*xr2-DS7*xr6;\
		xr2 = xr0+DS4*xr6-DSa*xr2;\
		*/\
		__Bbr = xr2;\
		  xr7 = xr6;\
		xr6 *= DS74;\
		xr2 *= DSa4;\
		xr6 += __Bbr;\
		xr2 +=   xr7;\
		xr6 *= DS4;\
		xr2 *= DS4;\
		xr6 += xr4;\
		xr2 += xr0;\
		/*\
		xr7 = xr4+DS5*xr1-DS8*xr5;\
		xr0 = xr0+DS5*xr5-DSb*xr1;\
		*/\
		__Bbr = xr1;\
		  xr7 = xr5;\
		xr5 *= DS85;\
		xr1 *= DSb5;\
		xr5 += __Bbr;\
		xr1 +=   xr7;\
		xr5 *= DS5;\
		xr1 *= DS5;\
		xr5 += xr4;\
		xr1 += xr0;\
		xr4 = __Bcr;\
		/*\
		xr7 = xr3+xr6;\
		xr6 = xr6-xr3+xr2;\
		xr2 = xr3+xr2;\
		xr0 = xr4+xr5;\
		xr5 = xr5-xr4+xr1;\
		xr4 = xr4+xr1;\
		*/\
		xr7 = xr3;\
		xr0 = xr4;\
		xr7+= xr6;\
		xr0+= xr5;\
		xr6-= xr3;\
		xr5-= xr4;\
		xr6+= xr2;\
		xr5+= xr1;\
		xr2+= xr3;\
		xr4+= xr1;\
	/* Combine xr and yi-terms to get real parts of outputs: */\
		xr1 = 2.0;\
		xr3 = __B6r;	xr3 -= xr7;	__B6r = xr3;	xr7 *= xr1;	xr7 += xr3;	__B7r = xr7;\
		xr3 = __B8r;	xr3 -= xr6;	__B8r = xr3;	xr6 *= xr1;	xr6 += xr3;	__B5r = xr6;\
		xr3 = __B2r;	xr3 -= xr2;	__B2r = xr3;	xr2 *= xr1;	xr2 += xr3;	__Bbr = xr2;\
		xr3 = __B3r;	xr3 -= xr0;	__B3r = xr3;	xr0 *= xr1;	xr0 += xr3;	__Bar = xr0;\
		xr3 = __B4r;	xr3 -= xr5;	__B4r = xr3;	xr5 *= xr1;	xr5 += xr3;	__B9r = xr5;\
		xr3 = __B1r;	xr3 -= xr4;	__B1r = xr3;	xr4 *= xr1;	xr4 += xr3;	__Bcr = xr4;\
		/*\
		__B6r -= xr7;	xr7 *= xr1;	xr7 += __B6r;	__B7r = xr7;\
		__B8r -= xr6;	xr6 *= xr1;	xr6 += __B8r;	__B5r = xr6;\
		__B2r -= xr2;	xr2 *= xr1;	xr2 += __B2r;	__Bbr = xr2;\
		__B3r -= xr0;	xr0 *= xr1;	xr0 += __B3r;	__Bar = xr0;\
		__B4r -= xr5;	xr5 *= xr1;	xr5 += __B4r;	__B9r = xr5;\
		__B1r -= xr4;	xr4 *= xr1;	xr4 += __B1r;	__Bcr = xr4;\
		*/\
	/***************/\
	/* IMAG PARTS: Swap __**r <--> __**i, Use __t*r (instead of __A*r) as inputs for yr-terms, Replace __B[j] with __B[13-j] for j > 0: */\
	/***************/\
	/* xi-terms need 8 registers for each side: */\
		xr7 = __A6i;\
		xr5 = __A5i;\
		xr4 = __A4i;\
		xr6 = __A3i;\
		xr3 = __A2i;\
		xr1 = __A1i;\
		xr7 += __A7i;\
		xr5 += __A8i;\
		xr4 += __A9i;\
		xr6 += __Aai;\
		xr3 += __Abi;\
		xr1 += __Aci;\
		xr1 -= xr5;	xr5 *= 2.0;\
		xr3 -= xr6;	xr6 *= 2.0;\
		xr4 -= xr7;	xr7 *= 2.0;\
		xr5 += xr1;\
		xr6 += xr3;\
		xr7 += xr4;\
		xr2 = xr5;\
		xr2 += xr6;\
		xr2 += xr7;\
		xr0 = __A0i;\
		xr0 += xr2;\
		xr2 *= DC1;\
		__B0i = xr0;\
		xr0 += xr2;\
		xr2 = DC3;\
		xr6 *= xr2;\
		xr5 *= xr2;\
		xr7 *= xr2;\
		__Bci = xr6;\
		__Bbi = xr5;\
		__Bai = xr7;\
		xr2 = DC23;\
		xr5 *= xr2;\
		xr7 *= xr2;\
		xr6 *= xr2;\
		xr5 += xr0;\
		xr7 += xr0;\
		xr6 += xr0;\
		xr5 += __Bci;\
		xr7 += __Bbi;\
		xr6 += __Bai;\
		/* Rearrange the 3x3 multiply block to ease in-place computation: */\
		xr2 = DC65;\
		xr0 = DC54;\
		__Bci = xr4;\
		__Bbi = xr1;\
		__Bai = xr3;\
		xr4 *= xr2;  /* *= DC65 */\
		xr1 *= xr2;  /* *= DC65 */\
		xr3 *= xr2;  /* *= DC65 */\
		xr4 += __Bai;\
		xr1 -= __Bci;\
		xr3 += __Bbi;\
		xr4 *= xr0;  /* *= DC54 */\
		xr1 *= xr0;  /* *= DC54 */\
		xr3 *= xr0;  /* *= DC54 */\
		xr2 = DC4;\
		xr4 += __Bbi;\
		xr1 -= __Bai;\
		xr3 -= __Bci;\
		xr4 *= xr2;  /* *= DC4 */\
		xr1 *= xr2;  /* *= DC4 */\
		xr3 *= xr2;  /* *= DC4 */\
	/* Spill into destination outputs: */\
		xr6 -= xr4;	xr4 *= 2.0;\
		xr7 -= xr1;	xr1 *= 2.0;\
		xr5 -= xr3;	xr3 *= 2.0;\
		__Bai = xr7;\
		__B5i = xr5;\
		xr4 += xr6;\
		xr7 += xr1;\
		xr5 += xr3;\
		__Bci = xr5;\
		__Bbi = xr7;\
		__B9i = xr6;\
		__B7i = xr4;\
	/* yr-terms: */\
		xr1 = __t1r;\
		xr2 = __t2r;\
		xr5 = __t3r;\
		xr3 = __t4r;\
		xr4 = __t5r;\
		xr6 = __t6r;\
		xr1 -= __tcr;\
		xr2 -= __tbr;\
		xr5 -= __tar;\
		xr3 -= __t9r;\
		xr4 -= __t8r;\
		xr6 -= __t7r;\
		xr7 = xr1;\
		xr0 = xr2;\
		xr7 -= xr3;\
		xr0 += xr4;\
		xr7 += xr5;\
		xr0 += xr6;\
		xr1 += xr3;\
		xr5 += xr3;\
		xr2 -= xr4;\
		xr6 -= xr4;\
		xr4 = xr0;\
		xr3   = xr7;\
		xr0 *= DS2;\
		xr7 *= DS2;\
		xr4 *= DS1;\
		xr3   *= DS1;\
		xr4 -= xr7;\
		xr3   += xr0;\
		__B1i = xr4;	/* tmp-store in B1i */\
		xr0 = xr1;\
		xr4 = xr5;\
		xr0 += xr2;\
		xr4 += xr6;\
		/*\
		xr7 = DS3*xr0-DS6*xr4;\
		xr0 = DS3*xr4-DS9*xr0;\
		*/\
		__B2i = xr0;	/* tmp-store in B2i */\
		  xr7 = xr4;\
		xr4 *= DS63;\
		xr0 *= DS93;\
		xr4 += __B2i;\
		xr0 +=   xr7;\
		xr4 *= DS3;\
		xr0 *= DS3;\
		/*\
		xr7 = xr4+DS4*xr2-DS7*xr6;\
		xr2 = xr0+DS4*xr6-DSa*xr2;\
		*/\
		__B2i = xr2;\
		  xr7 = xr6;\
		xr6 *= DS74;\
		xr2 *= DSa4;\
		xr6 += __B2i;\
		xr2 +=   xr7;\
		xr6 *= DS4;\
		xr2 *= DS4;\
		xr6 += xr4;\
		xr2 += xr0;\
		/*\
		xr7 = xr4+DS5*xr1-DS8*xr5;\
		xr0 = xr0+DS5*xr5-DSb*xr1;\
		*/\
		__B2i = xr1;\
		  xr7 = xr5;\
		xr5 *= DS85;\
		xr1 *= DSb5;\
		xr5 += __B2i;\
		xr1 +=   xr7;\
		xr5 *= DS5;\
		xr1 *= DS5;\
		xr5 += xr4;\
		xr1 += xr0;\
		xr4 = __B1i;\
		/*\
		xr7 = xr3+xr6;\
		xr6 = xr6-xr3+xr2;\
		xr2 = xr3+xr2;\
		xr0 = xr4+xr5;\
		xr5 = xr5-xr4+xr1;\
		xr4 = xr4+xr1;\
		*/\
		xr7 = xr3;\
		xr0 = xr4;\
		xr7+= xr6;\
		xr0+= xr5;\
		xr6-= xr3;\
		xr5-= xr4;\
		xr6+= xr2;\
		xr5+= xr1;\
		xr2+= xr3;\
		xr4+= xr1;\
	/* Combine xr and yi-terms to get real parts of outputs: */\
		xr1 = 2.0;\
		xr3 = __B7i;	xr3 -= xr7;	__B7i = xr3;	xr7 *= xr1;	xr7 += xr3;	__B6i = xr7;\
		xr3 = __B5i;	xr3 -= xr6;	__B5i = xr3;	xr6 *= xr1;	xr6 += xr3;	__B8i = xr6;\
		xr3 = __Bbi;	xr3 -= xr2;	__Bbi = xr3;	xr2 *= xr1;	xr2 += xr3;	__B2i = xr2;\
		xr3 = __Bai;	xr3 -= xr0;	__Bai = xr3;	xr0 *= xr1;	xr0 += xr3;	__B3i = xr0;\
		xr3 = __B9i;	xr3 -= xr5;	__B9i = xr3;	xr5 *= xr1;	xr5 += xr3;	__B4i = xr5;\
		xr3 = __Bci;	xr3 -= xr4;	__Bci = xr3;	xr4 *= xr1;	xr4 += xr3;	__B1i = xr4;\
/* Totals: 164 FADD, 88 FMUL. */\
}

#endif	/* #ifndef radix13_included */
