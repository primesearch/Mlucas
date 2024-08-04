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

#include "Mlucas.h"

#define RADIX 72	// Use #define rather than const int to ensure it's really a compile-time const in the C sense

/***************

int radix72_ditN_cy_dif1(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], double base[], double baseinv[], int iter, double *fracmax, uint64 p)
{
}

****************/

void radix72_dif_pass1(double a[], int n)
{
/*
!...Acronym: DIF = Decimation In Frequency
!
!...Subroutine to perform an initial radix-72 complex DIF FFT pass on the data in the length-N real vector A.
!
!   See the documentation in radix16_dif_pass for further details on storage and indexing.
!
!   See the documentation in radix5,16_dif_pass for details on the radix-5,16 subtransforms.
*/
	int j,j1,j2,jp,jt;
	static int NDIVR,p01,p02,p03,p04,p05,p06,p07,p08,p16,p24,p32,p40,p48,p56,p64, first_entry=TRUE;
	const double	c   =  0.76604444311897803520,	/* cos(2*pi/9) */
					s   =  0.64278760968653932631,	/* sin(2*pi/9) */
					c2  =  0.17364817766693034887,	/* cos(2*u) */
					s2  =  0.98480775301220805936,	/* sin(2*u) */
					c3m1= -1.50000000000000000000,	/* cos(3*u)-1 */
					s3  =  0.86602540378443864677,	/* sin(3*u) */
					c4  = -0.93969262078590838404,	/* cos(4*u) */
					s4  =  0.34202014332566873307;	/* sin(4*u) */
	double re,rt,it, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15
		,r00r,r00i,r01r,r01i,r02r,r02i,r03r,r03i,r04r,r04i,r05r,r05i,r06r,r06i,r07r,r07i,r08r,r08i,r09r,r09i
		,r10r,r10i,r11r,r11i,r12r,r12i,r13r,r13i,r14r,r14i,r15r,r15i,r16r,r16i,r17r,r17i,r18r,r18i,r19r,r19i
		,r20r,r20i,r21r,r21i,r22r,r22i,r23r,r23i,r24r,r24i,r25r,r25i,r26r,r26i,r27r,r27i,r28r,r28i,r29r,r29i
		,r30r,r30i,r31r,r31i,r32r,r32i,r33r,r33i,r34r,r34i,r35r,r35i,r36r,r36i,r37r,r37i,r38r,r38i,r39r,r39i
		,r40r,r40i,r41r,r41i,r42r,r42i,r43r,r43i,r44r,r44i,r45r,r45i,r46r,r46i,r47r,r47i,r48r,r48i,r49r,r49i
		,r50r,r50i,r51r,r51i,r52r,r52i,r53r,r53i,r54r,r54i,r55r,r55i,r56r,r56i,r57r,r57i,r58r,r58i,r59r,r59i
		,r60r,r60i,r61r,r61i,r62r,r62i,r63r,r63i,r64r,r64i,r65r,r65i,r66r,r66i,r67r,r67i,r68r,r68i,r69r,r69i
		,r70r,r70i,r71r,r71i;

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
		p16 = p08 + p08;
		p24 = p16 + p08;
		p32 = p24 + p08;
		p40 = p32 + p08;
		p48 = p40 + p08;
		p56 = p48 + p08;
		p64 = p56 + p08;

		p01 = p01 + ( (p01 >> DAT_BITS) << PAD_BITS );
		p02 = p02 + ( (p02 >> DAT_BITS) << PAD_BITS );
		p03 = p03 + ( (p03 >> DAT_BITS) << PAD_BITS );
		p04 = p04 + ( (p04 >> DAT_BITS) << PAD_BITS );
		p05 = p05 + ( (p05 >> DAT_BITS) << PAD_BITS );
		p06 = p06 + ( (p06 >> DAT_BITS) << PAD_BITS );
		p07 = p07 + ( (p07 >> DAT_BITS) << PAD_BITS );
		p08 = p08 + ( (p08 >> DAT_BITS) << PAD_BITS );
		p16 = p16 + ( (p16 >> DAT_BITS) << PAD_BITS );
		p24 = p24 + ( (p24 >> DAT_BITS) << PAD_BITS );
		p32 = p32 + ( (p32 >> DAT_BITS) << PAD_BITS );
		p40 = p40 + ( (p40 >> DAT_BITS) << PAD_BITS );
		p48 = p48 + ( (p48 >> DAT_BITS) << PAD_BITS );
		p56 = p56 + ( (p56 >> DAT_BITS) << PAD_BITS );
		p64 = p64 + ( (p64 >> DAT_BITS) << PAD_BITS );
	}

/*...The radix-72 pass is here.	*/

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
	Twiddleless version arranges 8 sets of radix-9 DFT inputs as follows:
	0 in upper left corner, decrement 8 horizontally and 9 vertically:

		00,64,56,48,40,32,24,16,08   00,64,56,48,40,32,24,16,08 + p00
		63,55,47,39,31,23,15,07,71   56,48,40,32,24,16,08,00,64 + p07
		54,46,38,30,22,14,06,70,62   48,40,32,24,16,08,00,64,56 + p06
		45,37,29,21,13,05,69,61,53 = 40,32,24,16,08,00,64,56,48 + p05
		36,28,20,12,04,68,60,52,44   32,24,16,08,00,64,56,48,40 + p04
		27,19,11,03,67,59,51,43,35   24,16,08,00,64,56,48,40,32 + p03
		18,10,02,66,58,50,42,34,26   16,08,00,64,56,48,40,32,24 + p02
		09,01,65,57,49,41,33,25,17   08,00,64,56,48,40,32,24,16 + p01
	*/
	/*...gather the needed data (72 64-bit complex) and do 8 radix-9 transforms...*/
		RADIX_09_DIF(a[j1    ],a[j2    ],a[j1+p64],a[j2+p64],a[j1+p56],a[j2+p56],a[j1+p48],a[j2+p48],a[j1+p40],a[j2+p40],a[j1+p32],a[j2+p32],a[j1+p24],a[j2+p24],a[j1+p16],a[j2+p16],a[j1+p08],a[j2+p08], r00r,r00i,r08r,r08i,r16r,r16i,r24r,r24i,r32r,r32i,r40r,r40i,r48r,r48i,r56r,r56i,r64r,r64i, rt,it,re);	jt = j1+p07; jp = j2+p07;
		RADIX_09_DIF(a[jt+p56],a[jp+p56],a[jt+p48],a[jp+p48],a[jt+p40],a[jp+p40],a[jt+p32],a[jp+p32],a[jt+p24],a[jp+p24],a[jt+p16],a[jp+p16],a[jt+p08],a[jp+p08],a[jt    ],a[jp    ],a[jt+p64],a[jp+p64], r01r,r01i,r09r,r09i,r17r,r17i,r25r,r25i,r33r,r33i,r41r,r41i,r49r,r49i,r57r,r57i,r65r,r65i, rt,it,re);	jt = j1+p06; jp = j2+p06;
		RADIX_09_DIF(a[jt+p48],a[jp+p48],a[jt+p40],a[jp+p40],a[jt+p32],a[jp+p32],a[jt+p24],a[jp+p24],a[jt+p16],a[jp+p16],a[jt+p08],a[jp+p08],a[jt    ],a[jp    ],a[jt+p64],a[jp+p64],a[jt+p56],a[jp+p56], r02r,r02i,r10r,r10i,r18r,r18i,r26r,r26i,r34r,r34i,r42r,r42i,r50r,r50i,r58r,r58i,r66r,r66i, rt,it,re);	jt = j1+p05; jp = j2+p05;
		RADIX_09_DIF(a[jt+p40],a[jp+p40],a[jt+p32],a[jp+p32],a[jt+p24],a[jp+p24],a[jt+p16],a[jp+p16],a[jt+p08],a[jp+p08],a[jt    ],a[jp    ],a[jt+p64],a[jp+p64],a[jt+p56],a[jp+p56],a[jt+p48],a[jp+p48], r03r,r03i,r11r,r11i,r19r,r19i,r27r,r27i,r35r,r35i,r43r,r43i,r51r,r51i,r59r,r59i,r67r,r67i, rt,it,re);	jt = j1+p04; jp = j2+p04;
		RADIX_09_DIF(a[jt+p32],a[jp+p32],a[jt+p24],a[jp+p24],a[jt+p16],a[jp+p16],a[jt+p08],a[jp+p08],a[jt    ],a[jp    ],a[jt+p64],a[jp+p64],a[jt+p56],a[jp+p56],a[jt+p48],a[jp+p48],a[jt+p40],a[jp+p40], r04r,r04i,r12r,r12i,r20r,r20i,r28r,r28i,r36r,r36i,r44r,r44i,r52r,r52i,r60r,r60i,r68r,r68i, rt,it,re);	jt = j1+p03; jp = j2+p03;
		RADIX_09_DIF(a[jt+p24],a[jp+p24],a[jt+p16],a[jp+p16],a[jt+p08],a[jp+p08],a[jt    ],a[jp    ],a[jt+p64],a[jp+p64],a[jt+p56],a[jp+p56],a[jt+p48],a[jp+p48],a[jt+p40],a[jp+p40],a[jt+p32],a[jp+p32], r05r,r05i,r13r,r13i,r21r,r21i,r29r,r29i,r37r,r37i,r45r,r45i,r53r,r53i,r61r,r61i,r69r,r69i, rt,it,re);	jt = j1+p02; jp = j2+p02;
		RADIX_09_DIF(a[jt+p16],a[jp+p16],a[jt+p08],a[jp+p08],a[jt    ],a[jp    ],a[jt+p64],a[jp+p64],a[jt+p56],a[jp+p56],a[jt+p48],a[jp+p48],a[jt+p40],a[jp+p40],a[jt+p32],a[jp+p32],a[jt+p24],a[jp+p24], r06r,r06i,r14r,r14i,r22r,r22i,r30r,r30i,r38r,r38i,r46r,r46i,r54r,r54i,r62r,r62i,r70r,r70i, rt,it,re);	jt = j1+p01; jp = j2+p01;
		RADIX_09_DIF(a[jt+p08],a[jp+p08],a[jt    ],a[jp    ],a[jt+p64],a[jp+p64],a[jt+p56],a[jp+p56],a[jt+p48],a[jp+p48],a[jt+p40],a[jp+p40],a[jt+p32],a[jp+p32],a[jt+p24],a[jp+p24],a[jt+p16],a[jp+p16], r07r,r07i,r15r,r15i,r23r,r23i,r31r,r31i,r39r,r39i,r47r,r47i,r55r,r55i,r63r,r63i,r71r,r71i, rt,it,re);

	/*...and now do 9 radix-8 transforms.
	Use the supercalafragalistic Ancient Chinese Secret index-munging formula [SACSIMPF]
	to properly permute the outputs of the ensuing 9 radix-8 DFTs to the required ordering, which in terms of our 9 DFTs is

		 0, 1, 3, 2, 7, 6, 5, 4   01327654 + p00
		64,65,67,66,71,70,69,68   01327654 + p64
		44,45,47,46,40,41,43,42   45760132 + p40
		18,19,16,17,20,21,23,22   23014576 + p16
		62,63,60,61,58,59,56,57 = 67452301 + p56
		33,32,34,35,38,39,36,37   10236745 + p32
		13,12,14,15, 9, 8,10,11   54671023 + p08
		51,50,49,48,53,52,54,55   32105467 + p48
		31,30,29,28,27,26,25,24   76543210 + p24
	*/
					 /*                                    inputs                                 */ /*                       intermediates                       */ /*                 outputs                   */
		RADIX_08_DIF(r00r,r00i,r01r,r01i,r02r,r02i,r03r,r03i,r04r,r04i,r05r,r05i,r06r,r06i,r07r,r07i, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15, a[j1    ],a[j2    ],a[j1+p01],a[j2+p01],a[j1+p03],a[j2+p03],a[j1+p02],a[j2+p02],a[j1+p07],a[j2+p07],a[j1+p06],a[j2+p06],a[j1+p05],a[j2+p05],a[j1+p04],a[j2+p04], rt,it);	jt = j1+p64; jp = j2+p64;
		RADIX_08_DIF(r08r,r08i,r09r,r09i,r10r,r10i,r11r,r11i,r12r,r12i,r13r,r13i,r14r,r14i,r15r,r15i, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15, a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04], rt,it);	jt = j1+p40; jp = j2+p40;
		RADIX_08_DIF(r16r,r16i,r17r,r17i,r18r,r18i,r19r,r19i,r20r,r20i,r21r,r21i,r22r,r22i,r23r,r23i, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15, a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02], rt,it);	jt = j1+p16; jp = j2+p16;
		RADIX_08_DIF(r24r,r24i,r25r,r25i,r26r,r26i,r27r,r27i,r28r,r28i,r29r,r29i,r30r,r30i,r31r,r31i, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15, a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06], rt,it);	jt = j1+p56; jp = j2+p56;
		RADIX_08_DIF(r32r,r32i,r33r,r33i,r34r,r34i,r35r,r35i,r36r,r36i,r37r,r37i,r38r,r38i,r39r,r39i, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15, a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01], rt,it);	jt = j1+p32; jp = j2+p32;
		RADIX_08_DIF(r40r,r40i,r41r,r41i,r42r,r42i,r43r,r43i,r44r,r44i,r45r,r45i,r46r,r46i,r47r,r47i, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15, a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05], rt,it);	jt = j1+p08; jp = j2+p08;
		RADIX_08_DIF(r48r,r48i,r49r,r49i,r50r,r50i,r51r,r51i,r52r,r52i,r53r,r53i,r54r,r54i,r55r,r55i, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15, a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03], rt,it);	jt = j1+p48; jp = j2+p48;
		RADIX_08_DIF(r56r,r56i,r57r,r57i,r58r,r58i,r59r,r59i,r60r,r60i,r61r,r61i,r62r,r62i,r63r,r63i, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15, a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07], rt,it);	jt = j1+p24; jp = j2+p24;
		RADIX_08_DIF(r64r,r64i,r65r,r65i,r66r,r66i,r67r,r67i,r68r,r68i,r69r,r69i,r70r,r70i,r71r,r71i, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15, a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ], rt,it);
	}
}

/***************/

void radix72_dit_pass1(double a[], int n)
{
/*
!...Acronym: DIF = Decimation In Frequency
!
!...Subroutine to perform an initial radix-72 complex DIT FFT pass on the data in the length-N real vector A.
!
!   See the documentation in radix16_dif_pass for further details on storage and indexing.
!
!   See the documentation in radix7_dif_pass for details on the radix-15 subtransforms.
*/
	int j,j1,j2,jp,jt;
	static int NDIVR,p01,p02,p03,p04,p05,p06,p07,p08,p16,p24,p32,p40,p48,p56,p64, first_entry=TRUE;
	const double	c   =  0.76604444311897803520,	/* cos(2*pi/9) */
					s   =  0.64278760968653932631,	/* sin(2*pi/9) */
					c2  =  0.17364817766693034887,	/* cos(2*u) */
					s2  =  0.98480775301220805936,	/* sin(2*u) */
					c3m1= -1.50000000000000000000,	/* cos(3*u)-1 */
					s3  =  0.86602540378443864677,	/* sin(3*u) */
					c4  = -0.93969262078590838404,	/* cos(4*u) */
					s4  =  0.34202014332566873307;	/* sin(4*u) */
	double re,rt,it, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15
		,r00r,r00i,r01r,r01i,r02r,r02i,r03r,r03i,r04r,r04i,r05r,r05i,r06r,r06i,r07r,r07i,r08r,r08i,r09r,r09i
		,r10r,r10i,r11r,r11i,r12r,r12i,r13r,r13i,r14r,r14i,r15r,r15i,r16r,r16i,r17r,r17i,r18r,r18i,r19r,r19i
		,r20r,r20i,r21r,r21i,r22r,r22i,r23r,r23i,r24r,r24i,r25r,r25i,r26r,r26i,r27r,r27i,r28r,r28i,r29r,r29i
		,r30r,r30i,r31r,r31i,r32r,r32i,r33r,r33i,r34r,r34i,r35r,r35i,r36r,r36i,r37r,r37i,r38r,r38i,r39r,r39i
		,r40r,r40i,r41r,r41i,r42r,r42i,r43r,r43i,r44r,r44i,r45r,r45i,r46r,r46i,r47r,r47i,r48r,r48i,r49r,r49i
		,r50r,r50i,r51r,r51i,r52r,r52i,r53r,r53i,r54r,r54i,r55r,r55i,r56r,r56i,r57r,r57i,r58r,r58i,r59r,r59i
		,r60r,r60i,r61r,r61i,r62r,r62i,r63r,r63i,r64r,r64i,r65r,r65i,r66r,r66i,r67r,r67i,r68r,r68i,r69r,r69i
		,r70r,r70i,r71r,r71i;

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
		p16 = p08 + p08;
		p24 = p16 + p08;
		p32 = p24 + p08;
		p40 = p32 + p08;
		p48 = p40 + p08;
		p56 = p48 + p08;
		p64 = p56 + p08;

		p01 = p01 + ( (p01 >> DAT_BITS) << PAD_BITS );
		p02 = p02 + ( (p02 >> DAT_BITS) << PAD_BITS );
		p03 = p03 + ( (p03 >> DAT_BITS) << PAD_BITS );
		p04 = p04 + ( (p04 >> DAT_BITS) << PAD_BITS );
		p05 = p05 + ( (p05 >> DAT_BITS) << PAD_BITS );
		p06 = p06 + ( (p06 >> DAT_BITS) << PAD_BITS );
		p07 = p07 + ( (p07 >> DAT_BITS) << PAD_BITS );
		p08 = p08 + ( (p08 >> DAT_BITS) << PAD_BITS );
		p16 = p16 + ( (p16 >> DAT_BITS) << PAD_BITS );
		p24 = p24 + ( (p24 >> DAT_BITS) << PAD_BITS );
		p32 = p32 + ( (p32 >> DAT_BITS) << PAD_BITS );
		p40 = p40 + ( (p40 >> DAT_BITS) << PAD_BITS );
		p48 = p48 + ( (p48 >> DAT_BITS) << PAD_BITS );
		p56 = p56 + ( (p56 >> DAT_BITS) << PAD_BITS );
		p64 = p64 + ( (p64 >> DAT_BITS) << PAD_BITS );
	}

/*...The radix-72 pass is here.	*/

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

		00,01,02,03,04,05,06,07,08,09,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,
		36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71
	map index-by-index [e.g. the index 01 gets replaced by 64, not the index at *position*...] to
		00,64,56,48,40,32,24,16,08,63,55,47,39,31,23,15,07,71,54,46,38,30,22,14,06,70,62,45,37,29,21,13,05,69,61,53,
		36,28,20,12,04,68,60,52,44,27,19,11,03,67,59,51,43,35,18,10,02,66,58,50,42,34,26,09,01,65,57,49,41,33,25,17.	(*)

	(***NOTE*** The following set of permutations can be auto-generated by running test_dft_radix()
	for the radix in question in TEST_TYPE = 1 [DIT] mode, skipping the actual DIT-pass step initially.)

	Remember, inputs to DIT are bit-reversed, so

		a[00,01,02,03,04,05,06,07|08,09,10,11,12,13,14,15|16,17,18,19,20,21,22,23|24,25,26,27,28,29,30,31|32,33,34,35,36,37,38,39|40,41,42,43,44,45,46,47|48,49,50,51,52,53,54,55|56,57,58,59,60,61,62,63|64,65,66,67,68,69,70,71] contain [cf. bit-reversal index array output by test_fft_radix()]
		x[00,36,18,54,09,45,27,63|03,39,21,57,12,48,30,66|06,42,24,60,15,51,33,69|01,37,19,55,10,46,28,64|04,40,22,58,13,49,31,67|07,43,25,61,16,52,34,70|02,38,20,56,11,47,29,65|05,41,23,59,14,50,32,68|08,44,26,62,17,53,35,71], which get swapped [using the permutation (*)] to
		x[00,36,54,18,63,27,45,09|48,12,30,66,39,03,21,57|24,60,06,42,15,51,69,33|64,28,46,10,55,19,37,01|40,04,22,58,31,67,13,49|16,52,70,34,07,43,61,25|56,20,38,02,47,11,29,65|32,68,14,50,23,59,05,41|08,44,62,26,71,35,53,17], which means the a-indices get swapped as
		a[00,01,03,02,07,06,05,04|13,12,14,15,09,08,10,11|18,19,16,17,20,21,23,22|31,30,29,28,27,26,25,24|33,32,34,35,38,39,36,37|44,45,47,46,40,41,43,42|51,50,49,48,53,52,54,55|62,63,60,61,58,59,56,57|64,65,67,66,71,70,69,68].

	Breaking the final index-row into octets are the 9 octets going into the radix-8 DFTs:

		RADIX_08_DIT(00,01,03,02,07,06,05,04)   01327654
		RADIX_08_DIT(13,12,14,15,09,08,10,11)   54671023 + p08
		RADIX_08_DIT(18,19,16,17,20,21,23,22)   23014576 + p16
		RADIX_08_DIT(31,30,29,28,27,26,25,24)   76543210 + p24
		RADIX_08_DIT(33,32,34,35,38,39,36,37) = 10236745 + p32
		RADIX_08_DIT(44,45,47,46,40,41,43,42)   45760132 + p40
		RADIX_08_DIT(51,50,49,48,53,52,54,55)   32105467 + p48
		RADIX_08_DIT(62,63,60,61,58,59,56,57)   67452301 + p56
		RADIX_08_DIT(64,65,67,66,71,70,69,68)   01327654 + p64
*/
	/*...gather the needed data (120 64-bit complex, i.e. 240 64-bit reals) and do 9 radix-8 transforms...*/
					/*                                    inputs                                  */ /*                         outputs                   */
		RADIX_08_DIT(a[j1    ],a[j2    ],a[j1+p01],a[j2+p01],a[j1+p03],a[j2+p03],a[j1+p02],a[j2+p02],a[j1+p07],a[j2+p07],a[j1+p06],a[j2+p06],a[j1+p05],a[j2+p05],a[j1+p04],a[j2+p04], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15, r00r,r00i,r01r,r01i,r02r,r02i,r03r,r03i,r04r,r04i,r05r,r05i,r06r,r06i,r07r,r07i, rt,it);	jt = j1+p08; jp = j2+p08;
		RADIX_08_DIT(a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15, r08r,r08i,r09r,r09i,r10r,r10i,r11r,r11i,r12r,r12i,r13r,r13i,r14r,r14i,r15r,r15i, rt,it);	jt = j1+p16; jp = j2+p16;
		RADIX_08_DIT(a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15, r16r,r16i,r17r,r17i,r18r,r18i,r19r,r19i,r20r,r20i,r21r,r21i,r22r,r22i,r23r,r23i, rt,it);	jt = j1+p24; jp = j2+p24;
		RADIX_08_DIT(a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15, r24r,r24i,r25r,r25i,r26r,r26i,r27r,r27i,r28r,r28i,r29r,r29i,r30r,r30i,r31r,r31i, rt,it);	jt = j1+p32; jp = j2+p32;
		RADIX_08_DIT(a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15, r32r,r32i,r33r,r33i,r34r,r34i,r35r,r35i,r36r,r36i,r37r,r37i,r38r,r38i,r39r,r39i, rt,it);	jt = j1+p40; jp = j2+p40;
		RADIX_08_DIT(a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15, r40r,r40i,r41r,r41i,r42r,r42i,r43r,r43i,r44r,r44i,r45r,r45i,r46r,r46i,r47r,r47i, rt,it);	jt = j1+p48; jp = j2+p48;
		RADIX_08_DIT(a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15, r48r,r48i,r49r,r49i,r50r,r50i,r51r,r51i,r52r,r52i,r53r,r53i,r54r,r54i,r55r,r55i, rt,it);	jt = j1+p56; jp = j2+p56;
		RADIX_08_DIT(a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15, r56r,r56i,r57r,r57i,r58r,r58i,r59r,r59i,r60r,r60i,r61r,r61i,r62r,r62i,r63r,r63i, rt,it);	jt = j1+p64; jp = j2+p64;
		RADIX_08_DIT(a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15, r64r,r64i,r65r,r65i,r66r,r66i,r67r,r67i,r68r,r68i,r69r,r69i,r70r,r70i,r71r,r71i, rt,it);

	/*...and now do 8 radix-9 transforms...

	Use the inverse supercalafragalistic Ancient Chinese Secret index-munging formula [iSACSIMPF]
	to properly permute the outputs of the 8 radix-9 DFTs:

		00,32,64,24,56,16,48,08,40   00,32,64,24,56,16,48,08,40
		09,41,01,33,65,25,57,17,49   08,40,00,32,64,24,56,16,48 + p01
		18,50,10,42,02,34,66,26,58   16,48,08,40,00,32,64,24,56 + p02
		27,59,19,51,11,43,03,35,67 = 24,56,16,48,08,40,00,32,64 + p03
		36,68,28,60,20,52,12,44,04   32,64,24,56,16,48,08,40,00 + p04
		45,05,37,69,29,61,21,53,13   40,00,32,64,24,56,16,48,08 + p05
		54,14,46,06,38,70,30,62,22   48,08,40,00,32,64,24,56,16 + p06
		63,23,55,15,47,07,39,71,31   56,16,48,08,40,00,32,64,24 + p07
	*/
		RADIX_09_DIT(r00r,r00i, r08r,r08i, r16r,r16i, r24r,r24i, r32r,r32i, r40r,r40i, r48r,r48i, r56r,r56i, r64r,r64i ,a[j1    ],a[j2    ],a[j1+p32],a[j2+p32],a[j1+p64],a[j2+p64],a[j1+p24],a[j2+p24],a[j1+p56],a[j2+p56],a[j1+p16],a[j2+p16],a[j1+p48],a[j2+p48],a[j1+p08],a[j2+p08],a[j1+p40],a[j2+p40], rt,it,re);	jt = j1+p01; jp = j2+p01;
		RADIX_09_DIT(r01r,r01i, r09r,r09i, r17r,r17i, r25r,r25i, r33r,r33i, r41r,r41i, r49r,r49i, r57r,r57i, r65r,r65i ,a[jt+p08],a[jp+p08],a[jt+p40],a[jp+p40],a[jt    ],a[jp    ],a[jt+p32],a[jp+p32],a[jt+p64],a[jp+p64],a[jt+p24],a[jp+p24],a[jt+p56],a[jp+p56],a[jt+p16],a[jp+p16],a[jt+p48],a[jp+p48], rt,it,re);	jt = j1+p02; jp = j2+p02;
		RADIX_09_DIT(r02r,r02i, r10r,r10i, r18r,r18i, r26r,r26i, r34r,r34i, r42r,r42i, r50r,r50i, r58r,r58i, r66r,r66i ,a[jt+p16],a[jp+p16],a[jt+p48],a[jp+p48],a[jt+p08],a[jp+p08],a[jt+p40],a[jp+p40],a[jt    ],a[jp    ],a[jt+p32],a[jp+p32],a[jt+p64],a[jp+p64],a[jt+p24],a[jp+p24],a[jt+p56],a[jp+p56], rt,it,re);	jt = j1+p03; jp = j2+p03;
		RADIX_09_DIT(r03r,r03i, r11r,r11i, r19r,r19i, r27r,r27i, r35r,r35i, r43r,r43i, r51r,r51i, r59r,r59i, r67r,r67i ,a[jt+p24],a[jp+p24],a[jt+p56],a[jp+p56],a[jt+p16],a[jp+p16],a[jt+p48],a[jp+p48],a[jt+p08],a[jp+p08],a[jt+p40],a[jp+p40],a[jt    ],a[jp    ],a[jt+p32],a[jp+p32],a[jt+p64],a[jp+p64], rt,it,re);	jt = j1+p04; jp = j2+p04;
		RADIX_09_DIT(r04r,r04i, r12r,r12i, r20r,r20i, r28r,r28i, r36r,r36i, r44r,r44i, r52r,r52i, r60r,r60i, r68r,r68i ,a[jt+p32],a[jp+p32],a[jt+p64],a[jp+p64],a[jt+p24],a[jp+p24],a[jt+p56],a[jp+p56],a[jt+p16],a[jp+p16],a[jt+p48],a[jp+p48],a[jt+p08],a[jp+p08],a[jt+p40],a[jp+p40],a[jt    ],a[jp    ], rt,it,re);	jt = j1+p05; jp = j2+p05;
		RADIX_09_DIT(r05r,r05i, r13r,r13i, r21r,r21i, r29r,r29i, r37r,r37i, r45r,r45i, r53r,r53i, r61r,r61i, r69r,r69i ,a[jt+p40],a[jp+p40],a[jt    ],a[jp    ],a[jt+p32],a[jp+p32],a[jt+p64],a[jp+p64],a[jt+p24],a[jp+p24],a[jt+p56],a[jp+p56],a[jt+p16],a[jp+p16],a[jt+p48],a[jp+p48],a[jt+p08],a[jp+p08], rt,it,re);	jt = j1+p06; jp = j2+p06;
		RADIX_09_DIT(r06r,r06i, r14r,r14i, r22r,r22i, r30r,r30i, r38r,r38i, r46r,r46i, r54r,r54i, r62r,r62i, r70r,r70i ,a[jt+p48],a[jp+p48],a[jt+p08],a[jp+p08],a[jt+p40],a[jp+p40],a[jt    ],a[jp    ],a[jt+p32],a[jp+p32],a[jt+p64],a[jp+p64],a[jt+p24],a[jp+p24],a[jt+p56],a[jp+p56],a[jt+p16],a[jp+p16], rt,it,re);	jt = j1+p07; jp = j2+p07;
		RADIX_09_DIT(r07r,r07i, r15r,r15i, r23r,r23i, r31r,r31i, r39r,r39i, r47r,r47i, r55r,r55i, r63r,r63i, r71r,r71i ,a[jt+p56],a[jp+p56],a[jt+p16],a[jp+p16],a[jt+p48],a[jp+p48],a[jt+p08],a[jp+p08],a[jt+p40],a[jp+p40],a[jt    ],a[jp    ],a[jt+p32],a[jp+p32],a[jt+p64],a[jp+p64],a[jt+p24],a[jp+p24], rt,it,re);	jt = j1+p08; jp = j2+p08;
	}
}

#undef RADIX
