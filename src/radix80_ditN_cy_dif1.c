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

#include "Mlucas.h"

#define RADIX 80	// Use #define rather than const int to ensure it's really a compile-time const in the C sense

/***************

int radix80_ditN_cy_dif1(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], double base[], double baseinv[], int iter, double *fracmax, uint64 p)
{
}

****************/

void radix80_dif_pass1(double a[], int n)
{
/*
!...Acronym: DIF = Decimation In Frequency
!
!...Subroutine to perform an initial radix-80 complex DIF FFT pass on the data in the length-N real vector A.
!
!   See the documentation in radix16_dif_pass for further details on storage and indexing.
!
!   See the documentation in radix5,16_dif_pass for details on the radix-5,16 subtransforms.
*/
	int j,j1,j2,jp,jt;
	static int NDIVR,p01,p02,p03,p04,p05,p06,p07,p08,p09,p10,p11,p12,p13,p14,p15,p16,p32,p48,p64, first_entry=TRUE;
	const double c = 0.9238795325112867561281831, s = 0.3826834323650897717284599,
				uc1 = -1.25000000000000000000,	/* [cos(u)+cos(2u)]/2-1 = -5/4; u = 2*pi/5 */
				uc2 =  0.55901699437494742409,	/* [cos(u)-cos(2u)]/2 */
				us1 =  0.95105651629515357211,	/*  sin(u) */
				us2 =  1.53884176858762670130,	/* [sin(u)+sin(2u)] */
				us3 =  0.36327126400268044292;	/* [sin(u)-sin(2u)] */
	double rt,it
	,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09
	,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19
	,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29
	,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39
	,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49
	,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59
	,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69
	,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79
	,t80,t81,t82,t83,t84,t85,t86,t87,t88,t89
	,t90,t91,t92,t93,t94,t95,t96,t97,t98,t99
	,tA0,tA1,tA2,tA3,tA4,tA5,tA6,tA7,tA8,tA9
	,tB0,tB1,tB2,tB3,tB4,tB5,tB6,tB7,tB8,tB9
	,tC0,tC1,tC2,tC3,tC4,tC5,tC6,tC7,tC8,tC9
	,tD0,tD1,tD2,tD3,tD4,tD5,tD6,tD7,tD8,tD9
	,tE0,tE1,tE2,tE3,tE4,tE5,tE6,tE7,tE8,tE9
	,tF0,tF1,tF2,tF3,tF4,tF5,tF6,tF7,tF8,tF9;

	if(!first_entry && (n/RADIX) != NDIVR)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		first_entry=FALSE;
		NDIVR=n/RADIX;

		/*   constant index offsets for array load/stores are here.	*/
		p01 = NDIVR;
		p02 = p01 + p01;
		p03 = p02 + p01;
		p04 = p03 + p01;
		p05 = p04 + p01;
		p06 = p05 + p01;
		p07 = p06 + p01;
		p08 = p07 + p01;
		p09 = p08 + p01;
		p10 = p09 + p01;
		p11 = p10 + p01;
		p12 = p11 + p01;
		p13 = p12 + p01;
		p14 = p13 + p01;
		p15 = p14 + p01;
		p16 = p15 + p01;
		p32 = p16 + p16;
		p48 = p32 + p16;
		p64 = p48 + p16;

		p01 = p01 + ( (p01 >> DAT_BITS) << PAD_BITS );
		p02 = p02 + ( (p02 >> DAT_BITS) << PAD_BITS );
		p03 = p03 + ( (p03 >> DAT_BITS) << PAD_BITS );
		p04 = p04 + ( (p04 >> DAT_BITS) << PAD_BITS );
		p05 = p05 + ( (p05 >> DAT_BITS) << PAD_BITS );
		p06 = p06 + ( (p06 >> DAT_BITS) << PAD_BITS );
		p07 = p07 + ( (p07 >> DAT_BITS) << PAD_BITS );
		p08 = p08 + ( (p08 >> DAT_BITS) << PAD_BITS );
		p09 = p09 + ( (p09 >> DAT_BITS) << PAD_BITS );
		p10 = p10 + ( (p10 >> DAT_BITS) << PAD_BITS );
		p11 = p11 + ( (p11 >> DAT_BITS) << PAD_BITS );
		p12 = p12 + ( (p12 >> DAT_BITS) << PAD_BITS );
		p13 = p13 + ( (p13 >> DAT_BITS) << PAD_BITS );
		p14 = p14 + ( (p14 >> DAT_BITS) << PAD_BITS );
		p15 = p15 + ( (p15 >> DAT_BITS) << PAD_BITS );
		p16 = p16 + ( (p16 >> DAT_BITS) << PAD_BITS );
		p32 = p32 + ( (p32 >> DAT_BITS) << PAD_BITS );
		p48 = p48 + ( (p48 >> DAT_BITS) << PAD_BITS );
		p64 = p64 + ( (p64 >> DAT_BITS) << PAD_BITS );
	}

/*...The radix-80 pass is here.	*/

	for(j=0; j < NDIVR; j += 2)
	{
	#ifdef USE_AVX
		j1 = (j & mask02) + br8[j&7];
	#elif defined(USE_SSE2)
		j1 = (j & mask01) + br4[j&3];
	#else
		j1 = j;
	#endif
		j1 =j1 + ( (j1>> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
		j2 = j1+RE_IM_STRIDE;
	/*
	Twiddleless version arranges 8 sets of radix-5 DFT inputs as follows:
	0 in upper left corner, decrement 16 horizontally and 5 vertically:

		RADIX_05_DFT(00,64,48,32,16)
		RADIX_05_DFT(75,59,43,27,11)
		RADIX_05_DFT(70,54,38,22,06)
		RADIX_05_DFT(65,49,33,17,01)
		RADIX_05_DFT(60,44,28,12,76)
		RADIX_05_DFT(55,39,23,07,71)
		RADIX_05_DFT(50,34,18,02,66)
		RADIX_05_DFT(45,29,13,77,61)
		RADIX_05_DFT(40,24,08,72,56)
		RADIX_05_DFT(35,19,03,67,51)
		RADIX_05_DFT(30,14,78,62,46)
		RADIX_05_DFT(25,09,73,57,41)
		RADIX_05_DFT(20,04,68,52,36)
		RADIX_05_DFT(15,79,63,47,31)
		RADIX_05_DFT(10,74,58,42,26)
		RADIX_05_DFT(05,69,53,37,21)

	Use the supercalafragalistic Ancient Chinese Secret index-munging permutation formula [SACSIMPF]
	to properly permute the ordering of the outputs of the ensuing 5 radix-16 DFTs:

		00,01,03,02,06,07,04,05,12,13,15,14,09,08,10,11
	64+	03,02,01,00,04,05,07,06,15,14,13,12,10,11,08,09
	48+	09,08,10,11,15,14,13,12,03,02,01,00,04,05,07,06
	16+	12,13,15,14,09,08,10,11,06,07,04,05,03,02,01,00
	32+	06,07,04,05,03,02,01,00,09,08,10,11,15,14,13,12
	*/
	/*...gather the needed data (80 64-bit complex) and do 16 radix-5 transforms...*/
		/*[--y3-] [--y4-] <<<<< swap last 2 outputs to undo swap of these in macro!: */
					/*    sincos:     */ /*                                            inputs                                             */ /*          outputs                  */
		RADIX_05_DFT(uc1,uc2,us1,us2,us3,a[j1    ],a[j2    ],a[j1+p64],a[j2+p64],a[j1+p48],a[j2+p48],a[j1+p32],a[j2+p32],a[j1+p16],a[j2+p16],t00,t01,t02,t03,t04,t05,t08,t09,t06,t07,rt,it);	jt = j1+p11; jp = j2+p11;
		RADIX_05_DFT(uc1,uc2,us1,us2,us3,a[jt+p64],a[jp+p64],a[jt+p48],a[jp+p48],a[jt+p32],a[jp+p32],a[jt+p16],a[jp+p16],a[jt    ],a[jp    ],t10,t11,t12,t13,t14,t15,t18,t19,t16,t17,rt,it);	jt = j1+p06; jp = j2+p06;
		RADIX_05_DFT(uc1,uc2,us1,us2,us3,a[jt+p64],a[jp+p64],a[jt+p48],a[jp+p48],a[jt+p32],a[jp+p32],a[jt+p16],a[jp+p16],a[jt    ],a[jp    ],t20,t21,t22,t23,t24,t25,t28,t29,t26,t27,rt,it);	jt = j1+p01; jp = j2+p01;
		RADIX_05_DFT(uc1,uc2,us1,us2,us3,a[jt+p64],a[jp+p64],a[jt+p48],a[jp+p48],a[jt+p32],a[jp+p32],a[jt+p16],a[jp+p16],a[jt    ],a[jp    ],t30,t31,t32,t33,t34,t35,t38,t39,t36,t37,rt,it);	jt = j1+p12; jp = j2+p12;
		RADIX_05_DFT(uc1,uc2,us1,us2,us3,a[jt+p48],a[jp+p48],a[jt+p32],a[jp+p32],a[jt+p16],a[jp+p16],a[jt    ],a[jp    ],a[jt+p64],a[jp+p64],t40,t41,t42,t43,t44,t45,t48,t49,t46,t47,rt,it);	jt = j1+p07; jp = j2+p07;
		RADIX_05_DFT(uc1,uc2,us1,us2,us3,a[jt+p48],a[jp+p48],a[jt+p32],a[jp+p32],a[jt+p16],a[jp+p16],a[jt    ],a[jp    ],a[jt+p64],a[jp+p64],t50,t51,t52,t53,t54,t55,t58,t59,t56,t57,rt,it);	jt = j1+p02; jp = j2+p02;
		RADIX_05_DFT(uc1,uc2,us1,us2,us3,a[jt+p48],a[jp+p48],a[jt+p32],a[jp+p32],a[jt+p16],a[jp+p16],a[jt    ],a[jp    ],a[jt+p64],a[jp+p64],t60,t61,t62,t63,t64,t65,t68,t69,t66,t67,rt,it);	jt = j1+p13; jp = j2+p13;
		RADIX_05_DFT(uc1,uc2,us1,us2,us3,a[jt+p32],a[jp+p32],a[jt+p16],a[jp+p16],a[jt    ],a[jp    ],a[jt+p64],a[jp+p64],a[jt+p48],a[jp+p48],t70,t71,t72,t73,t74,t75,t78,t79,t76,t77,rt,it);	jt = j1+p08; jp = j2+p08;
		RADIX_05_DFT(uc1,uc2,us1,us2,us3,a[jt+p32],a[jp+p32],a[jt+p16],a[jp+p16],a[jt    ],a[jp    ],a[jt+p64],a[jp+p64],a[jt+p48],a[jp+p48],t80,t81,t82,t83,t84,t85,t88,t89,t86,t87,rt,it);	jt = j1+p03; jp = j2+p03;
		RADIX_05_DFT(uc1,uc2,us1,us2,us3,a[jt+p32],a[jp+p32],a[jt+p16],a[jp+p16],a[jt    ],a[jp    ],a[jt+p64],a[jp+p64],a[jt+p48],a[jp+p48],t90,t91,t92,t93,t94,t95,t98,t99,t96,t97,rt,it);	jt = j1+p14; jp = j2+p14;
		RADIX_05_DFT(uc1,uc2,us1,us2,us3,a[jt+p16],a[jp+p16],a[jt    ],a[jp    ],a[jt+p64],a[jp+p64],a[jt+p48],a[jp+p48],a[jt+p32],a[jp+p32],tA0,tA1,tA2,tA3,tA4,tA5,tA8,tA9,tA6,tA7,rt,it);	jt = j1+p09; jp = j2+p09;
		RADIX_05_DFT(uc1,uc2,us1,us2,us3,a[jt+p16],a[jp+p16],a[jt    ],a[jp    ],a[jt+p64],a[jp+p64],a[jt+p48],a[jp+p48],a[jt+p32],a[jp+p32],tB0,tB1,tB2,tB3,tB4,tB5,tB8,tB9,tB6,tB7,rt,it);	jt = j1+p04; jp = j2+p04;
		RADIX_05_DFT(uc1,uc2,us1,us2,us3,a[jt+p16],a[jp+p16],a[jt    ],a[jp    ],a[jt+p64],a[jp+p64],a[jt+p48],a[jp+p48],a[jt+p32],a[jp+p32],tC0,tC1,tC2,tC3,tC4,tC5,tC8,tC9,tC6,tC7,rt,it);	jt = j1+p15; jp = j2+p15;
		RADIX_05_DFT(uc1,uc2,us1,us2,us3,a[jt    ],a[jp    ],a[jt+p64],a[jp+p64],a[jt+p48],a[jp+p48],a[jt+p32],a[jp+p32],a[jt+p16],a[jp+p16],tD0,tD1,tD2,tD3,tD4,tD5,tD8,tD9,tD6,tD7,rt,it);	jt = j1+p10; jp = j2+p10;
		RADIX_05_DFT(uc1,uc2,us1,us2,us3,a[jt    ],a[jp    ],a[jt+p64],a[jp+p64],a[jt+p48],a[jp+p48],a[jt+p32],a[jp+p32],a[jt+p16],a[jp+p16],tE0,tE1,tE2,tE3,tE4,tE5,tE8,tE9,tE6,tE7,rt,it);	jt = j1+p05; jp = j2+p05;
		RADIX_05_DFT(uc1,uc2,us1,us2,us3,a[jt    ],a[jp    ],a[jt+p64],a[jp+p64],a[jt+p48],a[jp+p48],a[jt+p32],a[jp+p32],a[jt+p16],a[jp+p16],tF0,tF1,tF2,tF3,tF4,tF5,tF8,tF9,tF6,tF7,rt,it)

		/*...and now do 5 radix-16 transforms:	*/
					 /*                                                          inputs                                                           */  /*                                                                                           outputs                                                                                                                                                                                                                      */  /* consts */
		RADIX_16_DIF(t00,t01,t10,t11,t20,t21,t30,t31,t40,t41,t50,t51,t60,t61,t70,t71,t80,t81,t90,t91,tA0,tA1,tB0,tB1,tC0,tC1,tD0,tD1,tE0,tE1,tF0,tF1, a[j1    ],a[j2    ],a[j1+p01],a[j2+p01],a[j1+p03],a[j2+p03],a[j1+p02],a[j2+p02],a[j1+p06],a[j2+p06],a[j1+p07],a[j2+p07],a[j1+p04],a[j2+p04],a[j1+p05],a[j2+p05],a[j1+p12],a[j2+p12],a[j1+p13],a[j2+p13],a[j1+p15],a[j2+p15],a[j1+p14],a[j2+p14],a[j1+p09],a[j2+p09],a[j1+p08],a[j2+p08],a[j1+p10],a[j2+p10],a[j1+p11],a[j2+p11], c,s);	jt = j1+p64; jp = j2+p64;
		RADIX_16_DIF(t02,t03,t12,t13,t22,t23,t32,t33,t42,t43,t52,t53,t62,t63,t72,t73,t82,t83,t92,t93,tA2,tA3,tB2,tB3,tC2,tC3,tD2,tD3,tE2,tE3,tF2,tF3, a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],a[jt+p15],a[jp+p15],a[jt+p14],a[jp+p14],a[jt+p13],a[jp+p13],a[jt+p12],a[jp+p12],a[jt+p10],a[jp+p10],a[jt+p11],a[jp+p11],a[jt+p08],a[jp+p08],a[jt+p09],a[jp+p09], c,s);	jt = j1+p48; jp = j2+p48;
		RADIX_16_DIF(t04,t05,t14,t15,t24,t25,t34,t35,t44,t45,t54,t55,t64,t65,t74,t75,t84,t85,t94,t95,tA4,tA5,tB4,tB5,tC4,tC5,tD4,tD5,tE4,tE5,tF4,tF5, a[jt+p09],a[jp+p09],a[jt+p08],a[jp+p08],a[jt+p10],a[jp+p10],a[jt+p11],a[jp+p11],a[jt+p15],a[jp+p15],a[jt+p14],a[jp+p14],a[jt+p13],a[jp+p13],a[jt+p12],a[jp+p12],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06], c,s);	jt = j1+p16; jp = j2+p16;
		RADIX_16_DIF(t06,t07,t16,t17,t26,t27,t36,t37,t46,t47,t56,t57,t66,t67,t76,t77,t86,t87,t96,t97,tA6,tA7,tB6,tB7,tC6,tC7,tD6,tD7,tE6,tE7,tF6,tF7, a[jt+p12],a[jp+p12],a[jt+p13],a[jp+p13],a[jt+p15],a[jp+p15],a[jt+p14],a[jp+p14],a[jt+p09],a[jp+p09],a[jt+p08],a[jp+p08],a[jt+p10],a[jp+p10],a[jt+p11],a[jp+p11],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ], c,s);	jt = j1+p32; jp = j2+p32;
		RADIX_16_DIF(t08,t09,t18,t19,t28,t29,t38,t39,t48,t49,t58,t59,t68,t69,t78,t79,t88,t89,t98,t99,tA8,tA9,tB8,tB9,tC8,tC9,tD8,tD9,tE8,tE9,tF8,tF9, a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p09],a[jp+p09],a[jt+p08],a[jp+p08],a[jt+p10],a[jp+p10],a[jt+p11],a[jp+p11],a[jt+p15],a[jp+p15],a[jt+p14],a[jp+p14],a[jt+p13],a[jp+p13],a[jt+p12],a[jp+p12], c,s);
	}
}

/***************/

void radix80_dit_pass1(double a[], int n)
{
/*
!...Acronym: DIF = Decimation In Frequency
!
!...Subroutine to perform an initial radix-80 complex DIT FFT pass on the data in the length-N real vector A.
!
!   See the documentation in radix16_dif_pass for further details on storage and indexing.
!
!   See the documentation in radix7_dif_pass for details on the radix-15 subtransforms.
*/
	int j,j1,j2,jp,jt;
	static int NDIVR,p01,p02,p03,p04,p05,p06,p07,p08,p09,p10,p11,p12,p13,p14,p15,p16,p32,p48,p64, first_entry=TRUE;
	const double c = 0.9238795325112867561281831, s = 0.3826834323650897717284599,
				uc1 = -1.25000000000000000000,	/* [cos(u)+cos(2u)]/2-1 = -5/4; u = 2*pi/5 */
				uc2 =  0.55901699437494742409,	/* [cos(u)-cos(2u)]/2 */
				us1 =  0.95105651629515357211,	/*  sin(u) */
				us2 =  1.53884176858762670130,	/* [sin(u)+sin(2u)] */
				us3 =  0.36327126400268044292;	/* [sin(u)-sin(2u)] */
	double rt,it
	,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09
	,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19
	,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29
	,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39
	,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49
	,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59
	,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69
	,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79
	,t80,t81,t82,t83,t84,t85,t86,t87,t88,t89
	,t90,t91,t92,t93,t94,t95,t96,t97,t98,t99
	,tA0,tA1,tA2,tA3,tA4,tA5,tA6,tA7,tA8,tA9
	,tB0,tB1,tB2,tB3,tB4,tB5,tB6,tB7,tB8,tB9
	,tC0,tC1,tC2,tC3,tC4,tC5,tC6,tC7,tC8,tC9
	,tD0,tD1,tD2,tD3,tD4,tD5,tD6,tD7,tD8,tD9
	,tE0,tE1,tE2,tE3,tE4,tE5,tE6,tE7,tE8,tE9
	,tF0,tF1,tF2,tF3,tF4,tF5,tF6,tF7,tF8,tF9;

	if(!first_entry && (n/RADIX) != NDIVR)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		first_entry=FALSE;
		NDIVR=n/RADIX;

		/*   constant index offsets for array load/stores are here.	*/
		p01 = NDIVR;
		p02 = p01 + p01;
		p03 = p02 + p01;
		p04 = p03 + p01;
		p05 = p04 + p01;
		p06 = p05 + p01;
		p07 = p06 + p01;
		p08 = p07 + p01;
		p09 = p08 + p01;
		p10 = p09 + p01;
		p11 = p10 + p01;
		p12 = p11 + p01;
		p13 = p12 + p01;
		p14 = p13 + p01;
		p15 = p14 + p01;
		p16 = p15 + p01;
		p32 = p16 + p16;
		p48 = p32 + p16;
		p64 = p48 + p16;

		p01 = p01 + ( (p01 >> DAT_BITS) << PAD_BITS );
		p02 = p02 + ( (p02 >> DAT_BITS) << PAD_BITS );
		p03 = p03 + ( (p03 >> DAT_BITS) << PAD_BITS );
		p04 = p04 + ( (p04 >> DAT_BITS) << PAD_BITS );
		p05 = p05 + ( (p05 >> DAT_BITS) << PAD_BITS );
		p06 = p06 + ( (p06 >> DAT_BITS) << PAD_BITS );
		p07 = p07 + ( (p07 >> DAT_BITS) << PAD_BITS );
		p08 = p08 + ( (p08 >> DAT_BITS) << PAD_BITS );
		p09 = p09 + ( (p09 >> DAT_BITS) << PAD_BITS );
		p10 = p10 + ( (p10 >> DAT_BITS) << PAD_BITS );
		p11 = p11 + ( (p11 >> DAT_BITS) << PAD_BITS );
		p12 = p12 + ( (p12 >> DAT_BITS) << PAD_BITS );
		p13 = p13 + ( (p13 >> DAT_BITS) << PAD_BITS );
		p14 = p14 + ( (p14 >> DAT_BITS) << PAD_BITS );
		p15 = p15 + ( (p15 >> DAT_BITS) << PAD_BITS );
		p16 = p16 + ( (p16 >> DAT_BITS) << PAD_BITS );
		p32 = p32 + ( (p32 >> DAT_BITS) << PAD_BITS );
		p48 = p48 + ( (p48 >> DAT_BITS) << PAD_BITS );
		p64 = p64 + ( (p64 >> DAT_BITS) << PAD_BITS );
	}

/*...The radix-40 pass is here.	*/

	for(j=0; j < NDIVR; j += 2)
	{
	#ifdef USE_AVX
		j1 = (j & mask02) + br8[j&7];
	#elif defined(USE_SSE2)
		j1 = (j & mask01) + br4[j&3];
	#else
		j1 = j;
	#endif
		j1 =j1 + ( (j1>> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
		j2 = j1+RE_IM_STRIDE;
	/*
	Twiddleless version uses same linear-index-vector-form permutation as in DIF:

		00,01,02,03,04,05,06,07,08,09,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79

		00,64,48,32,16,75,59,43,27,11,70,54,38,22,06,65,49,33,17,01,60,44,28,12,76,55,39,23,07,71,50,34,18,02,66,45,29,13,77,61,40,24,08,72,56,35,19,03,67,51,30,14,78,62,46,25,09,73,57,41,20,04,68,52,36,15,79,63,47,31,10,74,58,42,26,05,69,53,37,21.	(*)

	Remember, inputs to DIT are bit-reversed, so

		a[00,01,02,03,04,05,06,07,08,09,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79] contain
		x[00,40,20,60,10,50,30,70,05,45,25,65,15,55,35,75,01,41,21,61,11,51,31,71,06,46,26,66,16,56,36,76,02,42,22,62,12,52,32,72,07,47,27,67,17,57,37,77,03,43,23,63,13,53,33,73,08,48,28,68,18,58,38,78,04,44,24,64,14,54,34,74,09,49,29,69,19,59,39,79], which get swapped [using the permutation (*) on the index *values*] to
		x[00,40,60,20,70,30,50,10,75,35,55,15,65,25,45,05,64,24,44,04,54,14,34,74,59,19,39,79,49,09,29,69,48,08,28,68,38,78,18,58,43,03,23,63,33,73,13,53,32,72,12,52,22,62,02,42,27,67,07,47,17,57,77,37,16,56,76,36,06,46,66,26,11,51,71,31,01,41,61,21], which means the a-indices get swapped as
		a[00,01,03,02,07,06,05,04,15,14,13,12,11,10,09,08|67,66,65,64,69,68,70,71,77,76,78,79,73,72,74,75|57,56,58,59,62,63,60,61,49,48,50,51,54,55,52,53|38,39,36,37,34,35,32,33,42,43,40,41,44,45,47,46|28,29,31,30,24,25,27,26,20,21,23,22,16,17,19,18]. These are the 5 16-tets going into the radix-16 DFTs.
	*/
	/*...gather the needed data (80 64-bit complex) and do 5 radix-16 transforms:	*/
					 /*                                                          inputs                                                                                                                                                                                                                                                           */  /*                                                     outputs                                                           */  /* consts */
		RADIX_16_DIT(a[j1    ],a[j2    ],a[j1+p01],a[j2+p01],a[j1+p03],a[j2+p03],a[j1+p02],a[j2+p02],a[j1+p07],a[j2+p07],a[j1+p06],a[j2+p06],a[j1+p05],a[j2+p05],a[j1+p04],a[j2+p04],a[j1+p15],a[j2+p15],a[j1+p14],a[j2+p14],a[j1+p13],a[j2+p13],a[j1+p12],a[j2+p12],a[j1+p11],a[j2+p11],a[j1+p10],a[j2+p10],a[j1+p09],a[j2+p09],a[j1+p08],a[j2+p08], t00,t01,t10,t11,t20,t21,t30,t31,t40,t41,t50,t51,t60,t61,t70,t71,t80,t81,t90,t91,tA0,tA1,tB0,tB1,tC0,tC1,tD0,tD1,tE0,tE1,tF0,tF1, c,s);	jt = j1+p64; jp = j2+p64;
		RADIX_16_DIT(a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p13],a[jp+p13],a[jt+p12],a[jp+p12],a[jt+p14],a[jp+p14],a[jt+p15],a[jp+p15],a[jt+p09],a[jp+p09],a[jt+p08],a[jp+p08],a[jt+p10],a[jp+p10],a[jt+p11],a[jp+p11], t02,t03,t12,t13,t22,t23,t32,t33,t42,t43,t52,t53,t62,t63,t72,t73,t82,t83,t92,t93,tA2,tA3,tB2,tB3,tC2,tC3,tD2,tD3,tE2,tE3,tF2,tF3, c,s);	jt = j1+p48; jp = j2+p48;
		RADIX_16_DIT(a[jt+p09],a[jp+p09],a[jt+p08],a[jp+p08],a[jt+p10],a[jp+p10],a[jt+p11],a[jp+p11],a[jt+p14],a[jp+p14],a[jt+p15],a[jp+p15],a[jt+p12],a[jp+p12],a[jt+p13],a[jp+p13],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05], t04,t05,t14,t15,t24,t25,t34,t35,t44,t45,t54,t55,t64,t65,t74,t75,t84,t85,t94,t95,tA4,tA5,tB4,tB5,tC4,tC5,tD4,tD5,tE4,tE5,tF4,tF5, c,s);	jt = j1+p32; jp = j2+p32;
		RADIX_16_DIT(a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p10],a[jp+p10],a[jt+p11],a[jp+p11],a[jt+p08],a[jp+p08],a[jt+p09],a[jp+p09],a[jt+p12],a[jp+p12],a[jt+p13],a[jp+p13],a[jt+p15],a[jp+p15],a[jt+p14],a[jp+p14], t06,t07,t16,t17,t26,t27,t36,t37,t46,t47,t56,t57,t66,t67,t76,t77,t86,t87,t96,t97,tA6,tA7,tB6,tB7,tC6,tC7,tD6,tD7,tE6,tE7,tF6,tF7, c,s);	jt = j1+p16; jp = j2+p16;
		RADIX_16_DIT(a[jt+p12],a[jp+p12],a[jt+p13],a[jp+p13],a[jt+p15],a[jp+p15],a[jt+p14],a[jp+p14],a[jt+p08],a[jp+p08],a[jt+p09],a[jp+p09],a[jt+p11],a[jp+p11],a[jt+p10],a[jp+p10],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02], t08,t09,t18,t19,t28,t29,t38,t39,t48,t49,t58,t59,t68,t69,t78,t79,t88,t89,t98,t99,tA8,tA9,tB8,tB9,tC8,tC9,tD8,tD9,tE8,tE9,tF8,tF9, c,s);

	/*...and now do 16 radix-5 transforms:
	Use the inverse supercalafragalistic Ancient Chinese Secret index-munging formula [iSACSIMPF]
	to properly permute the outputs of the 16 radix-5 DFTs:
	*/
            	     /*   sincos     */ /*                                            inputs                                             */ /*          outputs                  */
		RADIX_05_DFT(uc1,uc2,us1,us2,us3,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,a[j1    ],a[j2    ],a[j1+p16],a[j2+p16],a[j1+p32],a[j2+p32],a[j1+p48],a[j2+p48],a[j1+p64],a[j2+p64],rt,it);	jt = j1+p15; jp = j2+p15;
		RADIX_05_DFT(uc1,uc2,us1,us2,us3,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,a[jt    ],a[jp    ],a[jt+p16],a[jp+p16],a[jt+p32],a[jp+p32],a[jt+p48],a[jp+p48],a[jt+p64],a[jp+p64],rt,it);	jt = j1+p14; jp = j2+p14;
		RADIX_05_DFT(uc1,uc2,us1,us2,us3,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,a[jt+p16],a[jp+p16],a[jt+p32],a[jp+p32],a[jt+p48],a[jp+p48],a[jt+p64],a[jp+p64],a[jt    ],a[jp    ],rt,it);	jt = j1+p13; jp = j2+p13;
		RADIX_05_DFT(uc1,uc2,us1,us2,us3,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,a[jt+p32],a[jp+p32],a[jt+p48],a[jp+p48],a[jt+p64],a[jp+p64],a[jt    ],a[jp    ],a[jt+p16],a[jp+p16],rt,it);	jt = j1+p12; jp = j2+p12;
		RADIX_05_DFT(uc1,uc2,us1,us2,us3,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,a[jt+p48],a[jp+p48],a[jt+p64],a[jp+p64],a[jt    ],a[jp    ],a[jt+p16],a[jp+p16],a[jt+p32],a[jp+p32],rt,it);	jt = j1+p11; jp = j2+p11;
		RADIX_05_DFT(uc1,uc2,us1,us2,us3,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,a[jt+p64],a[jp+p64],a[jt    ],a[jp    ],a[jt+p16],a[jp+p16],a[jt+p32],a[jp+p32],a[jt+p48],a[jp+p48],rt,it);	jt = j1+p10; jp = j2+p10;
		RADIX_05_DFT(uc1,uc2,us1,us2,us3,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,a[jt    ],a[jp    ],a[jt+p16],a[jp+p16],a[jt+p32],a[jp+p32],a[jt+p48],a[jp+p48],a[jt+p64],a[jp+p64],rt,it);	jt = j1+p09; jp = j2+p09;
		RADIX_05_DFT(uc1,uc2,us1,us2,us3,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,a[jt+p16],a[jp+p16],a[jt+p32],a[jp+p32],a[jt+p48],a[jp+p48],a[jt+p64],a[jp+p64],a[jt    ],a[jp    ],rt,it);	jt = j1+p08; jp = j2+p08;
		RADIX_05_DFT(uc1,uc2,us1,us2,us3,t80,t81,t82,t83,t84,t85,t86,t87,t88,t89,a[jt+p32],a[jp+p32],a[jt+p48],a[jp+p48],a[jt+p64],a[jp+p64],a[jt    ],a[jp    ],a[jt+p16],a[jp+p16],rt,it);	jt = j1+p07; jp = j2+p07;
		RADIX_05_DFT(uc1,uc2,us1,us2,us3,t90,t91,t92,t93,t94,t95,t96,t97,t98,t99,a[jt+p48],a[jp+p48],a[jt+p64],a[jp+p64],a[jt    ],a[jp    ],a[jt+p16],a[jp+p16],a[jt+p32],a[jp+p32],rt,it);	jt = j1+p06; jp = j2+p06;
		RADIX_05_DFT(uc1,uc2,us1,us2,us3,tA0,tA1,tA2,tA3,tA4,tA5,tA6,tA7,tA8,tA9,a[jt+p64],a[jp+p64],a[jt    ],a[jp    ],a[jt+p16],a[jp+p16],a[jt+p32],a[jp+p32],a[jt+p48],a[jp+p48],rt,it);	jt = j1+p05; jp = j2+p05;
		RADIX_05_DFT(uc1,uc2,us1,us2,us3,tB0,tB1,tB2,tB3,tB4,tB5,tB6,tB7,tB8,tB9,a[jt    ],a[jp    ],a[jt+p16],a[jp+p16],a[jt+p32],a[jp+p32],a[jt+p48],a[jp+p48],a[jt+p64],a[jp+p64],rt,it);	jt = j1+p04; jp = j2+p04;
		RADIX_05_DFT(uc1,uc2,us1,us2,us3,tC0,tC1,tC2,tC3,tC4,tC5,tC6,tC7,tC8,tC9,a[jt+p16],a[jp+p16],a[jt+p32],a[jp+p32],a[jt+p48],a[jp+p48],a[jt+p64],a[jp+p64],a[jt    ],a[jp    ],rt,it);	jt = j1+p03; jp = j2+p03;
		RADIX_05_DFT(uc1,uc2,us1,us2,us3,tD0,tD1,tD2,tD3,tD4,tD5,tD6,tD7,tD8,tD9,a[jt+p32],a[jp+p32],a[jt+p48],a[jp+p48],a[jt+p64],a[jp+p64],a[jt    ],a[jp    ],a[jt+p16],a[jp+p16],rt,it);	jt = j1+p02; jp = j2+p02;
		RADIX_05_DFT(uc1,uc2,us1,us2,us3,tE0,tE1,tE2,tE3,tE4,tE5,tE6,tE7,tE8,tE9,a[jt+p48],a[jp+p48],a[jt+p64],a[jp+p64],a[jt    ],a[jp    ],a[jt+p16],a[jp+p16],a[jt+p32],a[jp+p32],rt,it);	jt = j1+p01; jp = j2+p01;
		RADIX_05_DFT(uc1,uc2,us1,us2,us3,tF0,tF1,tF2,tF3,tF4,tF5,tF6,tF7,tF8,tF9,a[jt+p64],a[jp+p64],a[jt    ],a[jp    ],a[jt+p16],a[jp+p16],a[jt+p32],a[jp+p32],a[jt+p48],a[jp+p48],rt,it)
	}
}

#undef RADIX
