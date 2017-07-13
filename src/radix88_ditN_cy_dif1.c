/*******************************************************************************
*                                                                              *
*   (C) 1997-2015 by Ernst W. Mayer.                                           *
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

#define RADIX 88	// Use #define rather than const int to ensure it's really a compile-time const in the C sense

/***************

int radix88_ditN_cy_dif1(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], double base[], double baseinv[], int iter, double *fracmax, uint64 p)
{
}

****************/

void radix88_dif_pass1(double a[], int n)
{
/*
!...Acronym: DIF = Decimation In Frequency
!
!...Subroutine to perform an initial radix-88 complex DIF FFT pass on the data in the length-N real vector A.
!
!   See the documentation in radix16_dif_pass for further details on storage and indexing.
!
!   See the documentation in radix11,16_dif_pass for details on the radix-5,16 subtransforms.
*/
	int j,j1,j2,jp,jt;
	// In order to preserve length-2 numeric-index-offset property here, use hex for digit:
	static int NDIVR,p01,p02,p03,p04,p05,p06,p07,p08,p10,p18,p20,p28,p30,p38,p40,p48,p50, first_entry=TRUE;
	const double cc1 =  0.84125353283118116886,	/* Real part of exp(i*2*pi/11), the radix-11 fundamental sincos datum	*/
			ss1 =  0.54064081745559758210,	/* Imag part of exp(i*2*pi/11).	*/
			cc2 =  0.41541501300188642553,	/* cos(2u)	*/
			ss2 =  0.90963199535451837140,	/* sin(2u)	*/
			cc3 = -0.14231483827328514043,	/* cos(3u)	*/
			ss3 =  0.98982144188093273237,	/* sin(3u)	*/
			cc4 = -0.65486073394528506404,	/* cos(4u)	*/
			ss4 =  0.75574957435425828378,	/* sin(4u)	*/
			cc5 = -0.95949297361449738988,	/* cos(5u)	*/
			ss5 =  0.28173255684142969773;	/* sin(5u)	*/
	double rt,it, x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i
		,r00r,r00i,r01r,r01i,r02r,r02i,r03r,r03i,r04r,r04i,r05r,r05i,r06r,r06i,r07r,r07i,r08r,r08i,r09r,r09i
		,r10r,r10i,r11r,r11i,r12r,r12i,r13r,r13i,r14r,r14i,r15r,r15i,r16r,r16i,r17r,r17i,r18r,r18i,r19r,r19i
		,r20r,r20i,r21r,r21i,r22r,r22i,r23r,r23i,r24r,r24i,r25r,r25i,r26r,r26i,r27r,r27i,r28r,r28i,r29r,r29i
		,r30r,r30i,r31r,r31i,r32r,r32i,r33r,r33i,r34r,r34i,r35r,r35i,r36r,r36i,r37r,r37i,r38r,r38i,r39r,r39i
		,r40r,r40i,r41r,r41i,r42r,r42i,r43r,r43i,r44r,r44i,r45r,r45i,r46r,r46i,r47r,r47i,r48r,r48i,r49r,r49i
		,r50r,r50i,r51r,r51i,r52r,r52i,r53r,r53i,r54r,r54i,r55r,r55i,r56r,r56i,r57r,r57i,r58r,r58i,r59r,r59i
		,r60r,r60i,r61r,r61i,r62r,r62i,r63r,r63i,r64r,r64i,r65r,r65i,r66r,r66i,r67r,r67i,r68r,r68i,r69r,r69i
		,r70r,r70i,r71r,r71i,r72r,r72i,r73r,r73i,r74r,r74i,r75r,r75i,r76r,r76i,r77r,r77i,r78r,r78i,r79r,r79i
		,r80r,r80i,r81r,r81i,r82r,r82i,r83r,r83i,r84r,r84i,r85r,r85i,r86r,r86i,r87r,r87i;

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
		p10 = p08 + p08;
		p18 = p10 + p08;
		p20 = p18 + p08;
		p28 = p20 + p08;
		p30 = p28 + p08;
		p38 = p30 + p08;
		p40 = p38 + p08;
		p48 = p40 + p08;
		p50 = p48 + p08;

		p01 = p01 + ( (p01 >> DAT_BITS) << PAD_BITS );
		p02 = p02 + ( (p02 >> DAT_BITS) << PAD_BITS );
		p03 = p03 + ( (p03 >> DAT_BITS) << PAD_BITS );
		p04 = p04 + ( (p04 >> DAT_BITS) << PAD_BITS );
		p05 = p05 + ( (p05 >> DAT_BITS) << PAD_BITS );
		p06 = p06 + ( (p06 >> DAT_BITS) << PAD_BITS );
		p07 = p07 + ( (p07 >> DAT_BITS) << PAD_BITS );
		p08 = p08 + ( (p08 >> DAT_BITS) << PAD_BITS );
		p10 = p10 + ( (p10 >> DAT_BITS) << PAD_BITS );
		p18 = p18 + ( (p18 >> DAT_BITS) << PAD_BITS );
		p20 = p20 + ( (p20 >> DAT_BITS) << PAD_BITS );
		p28 = p28 + ( (p28 >> DAT_BITS) << PAD_BITS );
		p30 = p30 + ( (p30 >> DAT_BITS) << PAD_BITS );
		p38 = p38 + ( (p38 >> DAT_BITS) << PAD_BITS );
		p40 = p40 + ( (p40 >> DAT_BITS) << PAD_BITS );
		p48 = p48 + ( (p48 >> DAT_BITS) << PAD_BITS );
		p50 = p50 + ( (p50 >> DAT_BITS) << PAD_BITS );
	}

/*...The radix-88 pass is here.	*/

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
	Twiddleless version arranges 8 sets of radix-11 DFT inputs as follows:
	0 in upper left corner, decrement 8 horizontally and 11 vertically [rep'd in hex here to match p-indexing]:
	[Can get these by running test_fft_radix with TTYPE = 0]:

		00,50,48,40,38,30,28,20,18,10,08   00,50,48,40,38,30,28,20,18,10,08 + p00
		4d,45,3d,35,2d,25,1d,15,0d,05,55   48,40,38,30,28,20,18,10,08,00,50 + p05
		42,3a,32,2a,22,1a,12,0a,02,52,4a   40,38,30,28,20,18,10,08,00,50,48 + p02
		37,2f,27,1f,17,0f,07,57,4f,47,3f = 30,28,20,18,10,08,00,50,48,40,38 + p07
		2c,24,1c,14,0c,04,54,4c,44,3c,34   28,20,18,10,08,00,50,48,40,38,30 + p04
		21,19,11,09,01,51,49,41,39,31,29   20,18,10,08,00,50,48,40,38,30,28 + p01
		16,0e,06,56,4e,46,3e,36,2e,26,1e   10,08,00,50,48,40,38,30,28,20,18 + p06
		0b,03,53,4b,43,3b,33,2b,23,1b,13   08,00,50,48,40,38,30,28,20,18,10 + p03
	*/
	/*...gather the needed data (88 64-bit complex) and do 8 radix-9 transforms...*/
		RADIX_11_DFT_BASIC(a[j1    ],a[j2    ],a[j1+p50],a[j2+p50],a[j1+p48],a[j2+p48],a[j1+p40],a[j2+p40],a[j1+p38],a[j2+p38],a[j1+p30],a[j2+p30],a[j1+p28],a[j2+p28],a[j1+p20],a[j2+p20],a[j1+p18],a[j2+p18],a[j1+p10],a[j2+p10],a[j1+p08],a[j2+p08], r00r,r00i,r01r,r01i,r02r,r02i,r03r,r03i,r04r,r04i,r05r,r05i,r06r,r06i,r07r,r07i,r08r,r08i,r09r,r09i,r10r,r10i, cc1,cc2,cc3,cc4,cc5,ss1,ss2,ss3,ss4,ss5);	jt = j1+p05; jp = j2+p05;
		RADIX_11_DFT_BASIC(a[jt+p48],a[jp+p48],a[jt+p40],a[jp+p40],a[jt+p38],a[jp+p38],a[jt+p30],a[jp+p30],a[jt+p28],a[jp+p28],a[jt+p20],a[jp+p20],a[jt+p18],a[jp+p18],a[jt+p10],a[jp+p10],a[jt+p08],a[jp+p08],a[jt    ],a[jp    ],a[jt+p50],a[jp+p50], r11r,r11i,r12r,r12i,r13r,r13i,r14r,r14i,r15r,r15i,r16r,r16i,r17r,r17i,r18r,r18i,r19r,r19i,r20r,r20i,r21r,r21i, cc1,cc2,cc3,cc4,cc5,ss1,ss2,ss3,ss4,ss5);	jt = j1+p02; jp = j2+p02;
		RADIX_11_DFT_BASIC(a[jt+p40],a[jp+p40],a[jt+p38],a[jp+p38],a[jt+p30],a[jp+p30],a[jt+p28],a[jp+p28],a[jt+p20],a[jp+p20],a[jt+p18],a[jp+p18],a[jt+p10],a[jp+p10],a[jt+p08],a[jp+p08],a[jt    ],a[jp    ],a[jt+p50],a[jp+p50],a[jt+p48],a[jp+p48], r22r,r22i,r23r,r23i,r24r,r24i,r25r,r25i,r26r,r26i,r27r,r27i,r28r,r28i,r29r,r29i,r30r,r30i,r31r,r31i,r32r,r32i, cc1,cc2,cc3,cc4,cc5,ss1,ss2,ss3,ss4,ss5);	jt = j1+p07; jp = j2+p07;
		RADIX_11_DFT_BASIC(a[jt+p30],a[jp+p30],a[jt+p28],a[jp+p28],a[jt+p20],a[jp+p20],a[jt+p18],a[jp+p18],a[jt+p10],a[jp+p10],a[jt+p08],a[jp+p08],a[jt    ],a[jp    ],a[jt+p50],a[jp+p50],a[jt+p48],a[jp+p48],a[jt+p40],a[jp+p40],a[jt+p38],a[jp+p38], r33r,r33i,r34r,r34i,r35r,r35i,r36r,r36i,r37r,r37i,r38r,r38i,r39r,r39i,r40r,r40i,r41r,r41i,r42r,r42i,r43r,r43i, cc1,cc2,cc3,cc4,cc5,ss1,ss2,ss3,ss4,ss5);	jt = j1+p04; jp = j2+p04;
		RADIX_11_DFT_BASIC(a[jt+p28],a[jp+p28],a[jt+p20],a[jp+p20],a[jt+p18],a[jp+p18],a[jt+p10],a[jp+p10],a[jt+p08],a[jp+p08],a[jt    ],a[jp    ],a[jt+p50],a[jp+p50],a[jt+p48],a[jp+p48],a[jt+p40],a[jp+p40],a[jt+p38],a[jp+p38],a[jt+p30],a[jp+p30], r44r,r44i,r45r,r45i,r46r,r46i,r47r,r47i,r48r,r48i,r49r,r49i,r50r,r50i,r51r,r51i,r52r,r52i,r53r,r53i,r54r,r54i, cc1,cc2,cc3,cc4,cc5,ss1,ss2,ss3,ss4,ss5);	jt = j1+p01; jp = j2+p01;
		RADIX_11_DFT_BASIC(a[jt+p20],a[jp+p20],a[jt+p18],a[jp+p18],a[jt+p10],a[jp+p10],a[jt+p08],a[jp+p08],a[jt    ],a[jp    ],a[jt+p50],a[jp+p50],a[jt+p48],a[jp+p48],a[jt+p40],a[jp+p40],a[jt+p38],a[jp+p38],a[jt+p30],a[jp+p30],a[jt+p28],a[jp+p28], r55r,r55i,r56r,r56i,r57r,r57i,r58r,r58i,r59r,r59i,r60r,r60i,r61r,r61i,r62r,r62i,r63r,r63i,r64r,r64i,r65r,r65i, cc1,cc2,cc3,cc4,cc5,ss1,ss2,ss3,ss4,ss5);	jt = j1+p06; jp = j2+p06;
		RADIX_11_DFT_BASIC(a[jt+p10],a[jp+p10],a[jt+p08],a[jp+p08],a[jt    ],a[jp    ],a[jt+p50],a[jp+p50],a[jt+p48],a[jp+p48],a[jt+p40],a[jp+p40],a[jt+p38],a[jp+p38],a[jt+p30],a[jp+p30],a[jt+p28],a[jp+p28],a[jt+p20],a[jp+p20],a[jt+p18],a[jp+p18], r66r,r66i,r67r,r67i,r68r,r68i,r69r,r69i,r70r,r70i,r71r,r71i,r72r,r72i,r73r,r73i,r74r,r74i,r75r,r75i,r76r,r76i, cc1,cc2,cc3,cc4,cc5,ss1,ss2,ss3,ss4,ss5);	jt = j1+p03; jp = j2+p03;
		RADIX_11_DFT_BASIC(a[jt+p08],a[jp+p08],a[jt    ],a[jp    ],a[jt+p50],a[jp+p50],a[jt+p48],a[jp+p48],a[jt+p40],a[jp+p40],a[jt+p38],a[jp+p38],a[jt+p30],a[jp+p30],a[jt+p28],a[jp+p28],a[jt+p20],a[jp+p20],a[jt+p18],a[jp+p18],a[jt+p10],a[jp+p10], r77r,r77i,r78r,r78i,r79r,r79i,r80r,r80i,r81r,r81i,r82r,r82i,r83r,r83i,r84r,r84i,r85r,r85i,r86r,r86i,r87r,r87i, cc1,cc2,cc3,cc4,cc5,ss1,ss2,ss3,ss4,ss5);

	/*...and now do 13 radix-8 transforms.
	Use the supercalafragalistic Ancient Chinese Secret index-munging formula [SACSIMPF]
	to properly permute the outputs of the ensuing 13 radix-8 DFTs to the required ordering, which in terms of our 13 DFTs is

	[these indices are base-10]				[hex]
		00,01,02,03,05,04,07,06   01235476 + p00
		82,83,81,80,87,86,84,85   23107645 + p50
		77,76,79,78,74,75,73,72   54762310 + p48
		64,65,66,67,69,68,71,70   01235476 + p40
		62,63,61,60,56,57,58,59   67540123 + p38
		51,50,48,49,54,55,53,52   32016754 + p30
		44,45,46,47,43,42,40,41   45673201 + p28
		33,32,35,34,36,37,38,39 = 10324567 + p20
		31,30,28,29,25,24,27,26   76451032 + p18
		18,19,17,16,23,22,20,21   23107645 + p10
		13,12,15,14,10,11,09,08   54762310 + p08
	*/
					 /*                                    inputs                                 */ /*                       intermediates                       */ /*                 outputs                   */
		RADIX_08_DIF(r00r,r00i,r11r,r11i,r22r,r22i,r33r,r33i,r44r,r44i,r55r,r55i,r66r,r66i,r77r,r77i, x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i, a[j1    ],a[j2    ],a[j1+p01],a[j2+p01],a[j1+p02],a[j2+p02],a[j1+p03],a[j2+p03],a[j1+p05],a[j2+p05],a[j1+p04],a[j2+p04],a[j1+p07],a[j2+p07],a[j1+p06],a[j2+p06], rt,it);	jt = j1+p50; jp = j2+p50;
		RADIX_08_DIF(r01r,r01i,r12r,r12i,r23r,r23i,r34r,r34i,r45r,r45i,r56r,r56i,r67r,r67i,r78r,r78i, x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i, a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05], rt,it);	jt = j1+p48; jp = j2+p48;
		RADIX_08_DIF(r02r,r02i,r13r,r13i,r24r,r24i,r35r,r35i,r46r,r46i,r57r,r57i,r68r,r68i,r79r,r79i, x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i, a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ], rt,it);	jt = j1+p40; jp = j2+p40;
		RADIX_08_DIF(r03r,r03i,r14r,r14i,r25r,r25i,r36r,r36i,r47r,r47i,r58r,r58i,r69r,r69i,r80r,r80i, x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i, a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06], rt,it);	jt = j1+p38; jp = j2+p38;
		RADIX_08_DIF(r04r,r04i,r15r,r15i,r26r,r26i,r37r,r37i,r48r,r48i,r59r,r59i,r70r,r70i,r81r,r81i, x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i, a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03], rt,it);	jt = j1+p30; jp = j2+p30;
		RADIX_08_DIF(r05r,r05i,r16r,r16i,r27r,r27i,r38r,r38i,r49r,r49i,r60r,r60i,r71r,r71i,r82r,r82i, x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i, a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04], rt,it);	jt = j1+p28; jp = j2+p28;
		RADIX_08_DIF(r06r,r06i,r17r,r17i,r28r,r28i,r39r,r39i,r50r,r50i,r61r,r61i,r72r,r72i,r83r,r83i, x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i, a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01], rt,it);	jt = j1+p20; jp = j2+p20;
		RADIX_08_DIF(r07r,r07i,r18r,r18i,r29r,r29i,r40r,r40i,r51r,r51i,r62r,r62i,r73r,r73i,r84r,r84i, x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i, a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07], rt,it);	jt = j1+p18; jp = j2+p18;
		RADIX_08_DIF(r08r,r08i,r19r,r19i,r30r,r30i,r41r,r41i,r52r,r52i,r63r,r63i,r74r,r74i,r85r,r85i, x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i, a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02], rt,it);	jt = j1+p10; jp = j2+p10;
		RADIX_08_DIF(r09r,r09i,r20r,r20i,r31r,r31i,r42r,r42i,r53r,r53i,r64r,r64i,r75r,r75i,r86r,r86i, x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i, a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05], rt,it);	jt = j1+p08; jp = j2+p08;
		RADIX_08_DIF(r10r,r10i,r21r,r21i,r32r,r32i,r43r,r43i,r54r,r54i,r65r,r65i,r76r,r76i,r87r,r87i, x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i, a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ], rt,it);
	}
}

/***************/

void radix88_dit_pass1(double a[], int n)
{
/*
!...Acronym: DIF = Decimation In Frequency
!
!...Subroutine to perform an initial radix-88 complex DIT FFT pass on the data in the length-N real vector A.
!
!   See the documentation in radix16_dif_pass for further details on storage and indexing.
*/
	int j,j1,j2,jp,jt;
	// In order to preserve length-2 numeric-index-offset property here, use hex for digit:
	static int NDIVR,p01,p02,p03,p04,p05,p06,p07,p08,p10,p18,p20,p28,p30,p38,p40,p48,p50, first_entry=TRUE;
	const double cc1 =  0.84125353283118116886,	/* Real part of exp(i*2*pi/11), the radix-11 fundamental sincos datum	*/
			ss1 =  0.54064081745559758210,	/* Imag part of exp(i*2*pi/11).	*/
			cc2 =  0.41541501300188642553,	/* cos(2u)	*/
			ss2 =  0.90963199535451837140,	/* sin(2u)	*/
			cc3 = -0.14231483827328514043,	/* cos(3u)	*/
			ss3 =  0.98982144188093273237,	/* sin(3u)	*/
			cc4 = -0.65486073394528506404,	/* cos(4u)	*/
			ss4 =  0.75574957435425828378,	/* sin(4u)	*/
			cc5 = -0.95949297361449738988,	/* cos(5u)	*/
			ss5 =  0.28173255684142969773;	/* sin(5u)	*/
	double rt,it, x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i
		,r00r,r00i,r01r,r01i,r02r,r02i,r03r,r03i,r04r,r04i,r05r,r05i,r06r,r06i,r07r,r07i,r08r,r08i,r09r,r09i
		,r10r,r10i,r11r,r11i,r12r,r12i,r13r,r13i,r14r,r14i,r15r,r15i,r16r,r16i,r17r,r17i,r18r,r18i,r19r,r19i
		,r20r,r20i,r21r,r21i,r22r,r22i,r23r,r23i,r24r,r24i,r25r,r25i,r26r,r26i,r27r,r27i,r28r,r28i,r29r,r29i
		,r30r,r30i,r31r,r31i,r32r,r32i,r33r,r33i,r34r,r34i,r35r,r35i,r36r,r36i,r37r,r37i,r38r,r38i,r39r,r39i
		,r40r,r40i,r41r,r41i,r42r,r42i,r43r,r43i,r44r,r44i,r45r,r45i,r46r,r46i,r47r,r47i,r48r,r48i,r49r,r49i
		,r50r,r50i,r51r,r51i,r52r,r52i,r53r,r53i,r54r,r54i,r55r,r55i,r56r,r56i,r57r,r57i,r58r,r58i,r59r,r59i
		,r60r,r60i,r61r,r61i,r62r,r62i,r63r,r63i,r64r,r64i,r65r,r65i,r66r,r66i,r67r,r67i,r68r,r68i,r69r,r69i
		,r70r,r70i,r71r,r71i,r72r,r72i,r73r,r73i,r74r,r74i,r75r,r75i,r76r,r76i,r77r,r77i,r78r,r78i,r79r,r79i
		,r80r,r80i,r81r,r81i,r82r,r82i,r83r,r83i,r84r,r84i,r85r,r85i,r86r,r86i,r87r,r87i;

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
		p10 = p08 + p08;
		p18 = p10 + p08;
		p20 = p18 + p08;
		p28 = p20 + p08;
		p30 = p28 + p08;
		p38 = p30 + p08;
		p40 = p38 + p08;
		p48 = p40 + p08;
		p50 = p48 + p08;

		p01 = p01 + ( (p01 >> DAT_BITS) << PAD_BITS );
		p02 = p02 + ( (p02 >> DAT_BITS) << PAD_BITS );
		p03 = p03 + ( (p03 >> DAT_BITS) << PAD_BITS );
		p04 = p04 + ( (p04 >> DAT_BITS) << PAD_BITS );
		p05 = p05 + ( (p05 >> DAT_BITS) << PAD_BITS );
		p06 = p06 + ( (p06 >> DAT_BITS) << PAD_BITS );
		p07 = p07 + ( (p07 >> DAT_BITS) << PAD_BITS );
		p08 = p08 + ( (p08 >> DAT_BITS) << PAD_BITS );
		p10 = p10 + ( (p10 >> DAT_BITS) << PAD_BITS );
		p18 = p18 + ( (p18 >> DAT_BITS) << PAD_BITS );
		p20 = p20 + ( (p20 >> DAT_BITS) << PAD_BITS );
		p28 = p28 + ( (p28 >> DAT_BITS) << PAD_BITS );
		p30 = p30 + ( (p30 >> DAT_BITS) << PAD_BITS );
		p38 = p38 + ( (p38 >> DAT_BITS) << PAD_BITS );
		p40 = p40 + ( (p40 >> DAT_BITS) << PAD_BITS );
		p48 = p48 + ( (p48 >> DAT_BITS) << PAD_BITS );
		p50 = p50 + ( (p50 >> DAT_BITS) << PAD_BITS );
	}

/*...The radix-88 pass is here.	*/

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

		00,01,02,03,04,05,06,07,08,09,0a,0b,0c,0d,0e,0f,10,11,12,13,14,15,
		16,17,18,19,1a,1b,1c,1d,1e,1f,20,21,22,23,24,25,26,27,28,29,2a,2b,
		2c,2d,2e,2f,30,31,32,33,34,35,36,37,38,39,3a,3b,3c,3d,3e,3f,40,41,
		42,43,44,45,46,47,48,49,4a,4b,4c,4d,4e,4f,50,51,52,53,54,55,56,57.
	map index-by-index [e.g. the index 01 gets replaced by 0x60, not the index at *position*...] to
		00,50,48,40,38,30,28,20,18,10,08,4d,45,3d,35,2d,25,1d,15,0d,05,55,
		42,3a,32,2a,22,1a,12,0a,02,52,4a,37,2f,27,1f,17,0f,07,57,4f,47,3f,
		2c,24,1c,14,0c,04,54,4c,44,3c,34,21,19,11,09,01,51,49,41,39,31,29,
		16,0e,06,56,4e,46,3e,36,2e,26,1e,0b,03,53,4b,43,3b,33,2b,23,1b,13.	(*)

	(***NOTE*** The following set of permutations can be auto-generated by running test_dft_radix()
	for the radix in question in TEST_TYPE = 1 [DIT] mode, skipping the actual DIT-pass step initially.)

	Remember, inputs to DIT are bit-reversed, so [all indices hex here]

		a[00,01,02,03,04,05,06,07,08,09,0a,0b,0c,0d,0e,0f,10,11,12,13,14,15,16,17,18,19,1a,1b,1c,1d,1e,1f,20,21,22,23,24,25,26,27,28,29,2a,2b,2c,2d,2e,2f,30,31,32,33,34,35,36,37,38,39,3a,3b,3c,3d,3e,3f,40,41,42,43,44,45,46,47,48,49,4a,4b,4c,4d,4e,4f,50,51,52,53,54,55,56,57] contain [cf. bit-reversal index array output by test_fft_radix()]
		x[00,2c,16,42,0b,37,21,4d,01,2d,17,43,0c,38,22,4e,02,2e,18,44,0d,39,23,4f,03,2f,19,45,0e,3a,24,50,04,30,1a,46,0f,3b,25,51,05,31,1b,47,10,3c,26,52,06,32,1c,48,11,3d,27,53,07,33,1d,49,12,3e,28,54,08,34,1e,4a,13,3f,29,55,09,35,1f,4b,14,40,2a,56,0a,36,20,4c,15,41,2b,57], which get swapped [using the permutation (*)] to
		x[00,2c,42,16,4d,21,37,0b,50,24,3a,0e,45,19,2f,03,48,1c,32,06,3d,11,27,53,40,14,2a,56,35,09,1f,4b,38,0c,22,4e,2d,01,17,43,30,04,1a,46,25,51,0f,3b,28,54,12,3e,1d,49,07,33,20,4c,0a,36,15,41,57,2b,18,44,02,2e,0d,39,4f,23,10,3c,52,26,05,31,47,1b,08,34,4a,1e,55,29,3f,13], which means the a-indices get swapped as
		a[00,01,03,02,07,06,05,04,1f,1e,1d,1c,1b,1a,19,18,33,32,31,30,35,34,36,37,4d,4c,4e,4f,49,48,4a,4b,0d,0c,0e,0f,09,08,0a,0b,21,20,22,23,26,27,24,25,3e,3f,3c,3d,3a,3b,38,39,52,53,50,51,54,55,57,56,12,13,10,11,14,15,17,16,2c,2d,2f,2e,28,29,2b,2a,40,41,43,42,47,46,45,44].

	Breaking the final index-row into octets are the 11 octets going into the radix-8 DFTs:

					[all indices hex here]
		00,01,03,02,07,06,05,04   01327654 + p00
		1f,1e,1d,1c,1b,1a,19,18   76543210 + p18
		33,32,31,30,35,34,36,37   32105467 + p30
		4d,4c,4e,4f,49,48,4a,4b   54671023 + p48
		0d,0c,0e,0f,09,08,0a,0b   54671023 + p08
		21,20,22,23,26,27,24,25   10236745 + p20
		3e,3f,3c,3d,3a,3b,38,39   67452301 + p38
		52,53,50,51,54,55,57,56 = 23014576 + p50
		12,13,10,11,14,15,17,16   23014576 + p10
		2c,2d,2f,2e,28,29,2b,2a   45760132 + p28
		40,41,43,42,47,46,45,44   01327654 + p40
*/
	/*...gather the needed data (88 64-bit complex) and do 11 radix-8 transforms...*/
		RADIX_08_DIT(a[j1    ],a[j2    ],a[j1+p01],a[j2+p01],a[j1+p03],a[j2+p03],a[j1+p02],a[j2+p02],a[j1+p07],a[j2+p07],a[j1+p06],a[j2+p06],a[j1+p05],a[j2+p05],a[j1+p04],a[j2+p04], x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i, r00r,r00i,r11r,r11i,r22r,r22i,r33r,r33i,r44r,r44i,r55r,r55i,r66r,r66i,r77r,r77i, rt,it);	jt = j1+p18; jp = j2+p18;
		RADIX_08_DIT(a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ], x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i, r01r,r01i,r12r,r12i,r23r,r23i,r34r,r34i,r45r,r45i,r56r,r56i,r67r,r67i,r78r,r78i, rt,it);	jt = j1+p30; jp = j2+p30;
		RADIX_08_DIT(a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07], x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i, r02r,r02i,r13r,r13i,r24r,r24i,r35r,r35i,r46r,r46i,r57r,r57i,r68r,r68i,r79r,r79i, rt,it);	jt = j1+p48; jp = j2+p48;
		RADIX_08_DIT(a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03], x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i, r03r,r03i,r14r,r14i,r25r,r25i,r36r,r36i,r47r,r47i,r58r,r58i,r69r,r69i,r80r,r80i, rt,it);	jt = j1+p08; jp = j2+p08;
		RADIX_08_DIT(a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03], x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i, r04r,r04i,r15r,r15i,r26r,r26i,r37r,r37i,r48r,r48i,r59r,r59i,r70r,r70i,r81r,r81i, rt,it);	jt = j1+p20; jp = j2+p20;
		RADIX_08_DIT(a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05], x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i, r05r,r05i,r16r,r16i,r27r,r27i,r38r,r38i,r49r,r49i,r60r,r60i,r71r,r71i,r82r,r82i, rt,it);	jt = j1+p38; jp = j2+p38;
		RADIX_08_DIT(a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01], x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i, r06r,r06i,r17r,r17i,r28r,r28i,r39r,r39i,r50r,r50i,r61r,r61i,r72r,r72i,r83r,r83i, rt,it);	jt = j1+p50; jp = j2+p50;
		RADIX_08_DIT(a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06], x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i, r07r,r07i,r18r,r18i,r29r,r29i,r40r,r40i,r51r,r51i,r62r,r62i,r73r,r73i,r84r,r84i, rt,it);	jt = j1+p10; jp = j2+p10;
		RADIX_08_DIT(a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06], x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i, r08r,r08i,r19r,r19i,r30r,r30i,r41r,r41i,r52r,r52i,r63r,r63i,r74r,r74i,r85r,r85i, rt,it);	jt = j1+p28; jp = j2+p28;
		RADIX_08_DIT(a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02], x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i, r09r,r09i,r20r,r20i,r31r,r31i,r42r,r42i,r53r,r53i,r64r,r64i,r75r,r75i,r86r,r86i, rt,it);	jt = j1+p40; jp = j2+p40;
		RADIX_08_DIT(a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04], x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i, r10r,r10i,r21r,r21i,r32r,r32i,r43r,r43i,r54r,r54i,r65r,r65i,r76r,r76i,r87r,r87i, rt,it);

	/*...and now do 8 radix-11 transforms: 

	Use the inverse supercalafragalistic Ancient Chinese Secret index-munging formula [iSACSIMPF]
	to properly permute the outputs of the 8 radix-11 DFTs:

								[all indices hex here]
		00,38,18,50,30,10,48,28,08,40,20   00,38,18,50,30,10,48,28,08,40,20
		21,01,39,19,51,31,11,49,29,09,41   20,00,38,18,50,30,10,48,28,08,40 + p01
		42,22,02,3a,1a,52,32,12,4a,2a,0a   40,20,00,38,18,50,30,10,48,28,08 + p02
		0b,43,23,03,3b,1b,53,33,13,4b,2b = 08,40,20,00,38,18,50,30,10,48,28 + p03
		2c,0c,44,24,04,3c,1c,54,34,14,4c   28,08,40,20,00,38,18,50,30,10,48 + p04
		4d,2d,0d,45,25,05,3d,1d,55,35,15   48,28,08,40,20,00,38,18,50,30,10 + p05
		16,4e,2e,0e,46,26,06,3e,1e,56,36   10,48,28,08,40,20,00,38,18,50,30 + p06
		37,17,4f,2f,0f,47,27,07,3f,1f,57   30,10,48,28,08,40,20,00,38,18,50 + p07
	*/
		RADIX_11_DFT_BASIC(r00r,r00i,r01r,r01i,r02r,r02i,r03r,r03i,r04r,r04i,r05r,r05i,r06r,r06i,r07r,r07i,r08r,r08i,r09r,r09i,r10r,r10i, a[j1    ],a[j2    ],a[j1+p38],a[j2+p38],a[j1+p18],a[j2+p18],a[j1+p50],a[j2+p50],a[j1+p30],a[j2+p30],a[j1+p10],a[j2+p10],a[j1+p48],a[j2+p48],a[j1+p28],a[j2+p28],a[j1+p08],a[j2+p08],a[j1+p40],a[j2+p40],a[j1+p20],a[j2+p20], cc1,cc2,cc3,cc4,cc5,ss1,ss2,ss3,ss4,ss5);	jt = j1+p01; jp = j2+p01;
		RADIX_11_DFT_BASIC(r11r,r11i,r12r,r12i,r13r,r13i,r14r,r14i,r15r,r15i,r16r,r16i,r17r,r17i,r18r,r18i,r19r,r19i,r20r,r20i,r21r,r21i, a[jt+p20],a[jp+p20],a[jt    ],a[jp    ],a[jt+p38],a[jp+p38],a[jt+p18],a[jp+p18],a[jt+p50],a[jp+p50],a[jt+p30],a[jp+p30],a[jt+p10],a[jp+p10],a[jt+p48],a[jp+p48],a[jt+p28],a[jp+p28],a[jt+p08],a[jp+p08],a[jt+p40],a[jp+p40], cc1,cc2,cc3,cc4,cc5,ss1,ss2,ss3,ss4,ss5);	jt = j1+p02; jp = j2+p02;
		RADIX_11_DFT_BASIC(r22r,r22i,r23r,r23i,r24r,r24i,r25r,r25i,r26r,r26i,r27r,r27i,r28r,r28i,r29r,r29i,r30r,r30i,r31r,r31i,r32r,r32i, a[jt+p40],a[jp+p40],a[jt+p20],a[jp+p20],a[jt    ],a[jp    ],a[jt+p38],a[jp+p38],a[jt+p18],a[jp+p18],a[jt+p50],a[jp+p50],a[jt+p30],a[jp+p30],a[jt+p10],a[jp+p10],a[jt+p48],a[jp+p48],a[jt+p28],a[jp+p28],a[jt+p08],a[jp+p08], cc1,cc2,cc3,cc4,cc5,ss1,ss2,ss3,ss4,ss5);	jt = j1+p03; jp = j2+p03;
		RADIX_11_DFT_BASIC(r33r,r33i,r34r,r34i,r35r,r35i,r36r,r36i,r37r,r37i,r38r,r38i,r39r,r39i,r40r,r40i,r41r,r41i,r42r,r42i,r43r,r43i, a[jt+p08],a[jp+p08],a[jt+p40],a[jp+p40],a[jt+p20],a[jp+p20],a[jt    ],a[jp    ],a[jt+p38],a[jp+p38],a[jt+p18],a[jp+p18],a[jt+p50],a[jp+p50],a[jt+p30],a[jp+p30],a[jt+p10],a[jp+p10],a[jt+p48],a[jp+p48],a[jt+p28],a[jp+p28], cc1,cc2,cc3,cc4,cc5,ss1,ss2,ss3,ss4,ss5);	jt = j1+p04; jp = j2+p04;
		RADIX_11_DFT_BASIC(r44r,r44i,r45r,r45i,r46r,r46i,r47r,r47i,r48r,r48i,r49r,r49i,r50r,r50i,r51r,r51i,r52r,r52i,r53r,r53i,r54r,r54i, a[jt+p28],a[jp+p28],a[jt+p08],a[jp+p08],a[jt+p40],a[jp+p40],a[jt+p20],a[jp+p20],a[jt    ],a[jp    ],a[jt+p38],a[jp+p38],a[jt+p18],a[jp+p18],a[jt+p50],a[jp+p50],a[jt+p30],a[jp+p30],a[jt+p10],a[jp+p10],a[jt+p48],a[jp+p48], cc1,cc2,cc3,cc4,cc5,ss1,ss2,ss3,ss4,ss5);	jt = j1+p05; jp = j2+p05;
		RADIX_11_DFT_BASIC(r55r,r55i,r56r,r56i,r57r,r57i,r58r,r58i,r59r,r59i,r60r,r60i,r61r,r61i,r62r,r62i,r63r,r63i,r64r,r64i,r65r,r65i, a[jt+p48],a[jp+p48],a[jt+p28],a[jp+p28],a[jt+p08],a[jp+p08],a[jt+p40],a[jp+p40],a[jt+p20],a[jp+p20],a[jt    ],a[jp    ],a[jt+p38],a[jp+p38],a[jt+p18],a[jp+p18],a[jt+p50],a[jp+p50],a[jt+p30],a[jp+p30],a[jt+p10],a[jp+p10], cc1,cc2,cc3,cc4,cc5,ss1,ss2,ss3,ss4,ss5);	jt = j1+p06; jp = j2+p06;
		RADIX_11_DFT_BASIC(r66r,r66i,r67r,r67i,r68r,r68i,r69r,r69i,r70r,r70i,r71r,r71i,r72r,r72i,r73r,r73i,r74r,r74i,r75r,r75i,r76r,r76i, a[jt+p10],a[jp+p10],a[jt+p48],a[jp+p48],a[jt+p28],a[jp+p28],a[jt+p08],a[jp+p08],a[jt+p40],a[jp+p40],a[jt+p20],a[jp+p20],a[jt    ],a[jp    ],a[jt+p38],a[jp+p38],a[jt+p18],a[jp+p18],a[jt+p50],a[jp+p50],a[jt+p30],a[jp+p30], cc1,cc2,cc3,cc4,cc5,ss1,ss2,ss3,ss4,ss5);	jt = j1+p07; jp = j2+p07;
		RADIX_11_DFT_BASIC(r77r,r77i,r78r,r78i,r79r,r79i,r80r,r80i,r81r,r81i,r82r,r82i,r83r,r83i,r84r,r84i,r85r,r85i,r86r,r86i,r87r,r87i, a[jt+p30],a[jp+p30],a[jt+p10],a[jp+p10],a[jt+p48],a[jp+p48],a[jt+p28],a[jp+p28],a[jt+p08],a[jp+p08],a[jt+p40],a[jp+p40],a[jt+p20],a[jp+p20],a[jt    ],a[jp    ],a[jt+p38],a[jp+p38],a[jt+p18],a[jp+p18],a[jt+p50],a[jp+p50], cc1,cc2,cc3,cc4,cc5,ss1,ss2,ss3,ss4,ss5);	
	}
}

#undef RADIX
