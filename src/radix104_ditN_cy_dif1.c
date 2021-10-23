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
#include "radix13.h"

#define RADIX 104	// Use #define rather than const int to ensure it's really a compile-time const in the C sense

/***************

int radix104_ditN_cy_dif1(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], double base[], double baseinv[], int iter, double *fracmax, uint64 p)
{
}

****************/

void radix104_dif_pass1(double a[], int n)
{
/*
!...Acronym: DIF = Decimation In Frequency
!
!...Subroutine to perform an initial radix-104 complex DIF FFT pass on the data in the length-N real vector A.
!
!   See the documentation in radix16_dif_pass for further details on storage and indexing.
!
!   See the documentation in radix5,16_dif_pass for details on the radix-5,16 subtransforms.
*/
	int j,j1,j2,jp,jt;
	// In order to preserve length-2 numeric-index-offset property here, use hex for digit:
	static int NDIVR,p01,p02,p03,p04,p05,p06,p07,p08,p10,p18,p20,p28,p30,p38,p40,p48,p50,p58,p60, first_entry=TRUE;
	double rt,it, x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i
		,r00r,r00i,r01r,r01i,r02r,r02i,r03r,r03i,r04r,r04i,r05r,r05i,r06r,r06i,r07r,r07i,r08r,r08i,r09r,r09i
		,r10r,r10i,r11r,r11i,r12r,r12i,r13r,r13i,r14r,r14i,r15r,r15i,r16r,r16i,r17r,r17i,r18r,r18i,r19r,r19i
		,r20r,r20i,r21r,r21i,r22r,r22i,r23r,r23i,r24r,r24i,r25r,r25i,r26r,r26i,r27r,r27i,r28r,r28i,r29r,r29i
		,r30r,r30i,r31r,r31i,r32r,r32i,r33r,r33i,r34r,r34i,r35r,r35i,r36r,r36i,r37r,r37i,r38r,r38i,r39r,r39i
		,r40r,r40i,r41r,r41i,r42r,r42i,r43r,r43i,r44r,r44i,r45r,r45i,r46r,r46i,r47r,r47i,r48r,r48i,r49r,r49i
		,r50r,r50i,r51r,r51i,r52r,r52i,r53r,r53i,r54r,r54i,r55r,r55i,r56r,r56i,r57r,r57i,r58r,r58i,r59r,r59i
		,r60r,r60i,r61r,r61i,r62r,r62i,r63r,r63i,r64r,r64i,r65r,r65i,r66r,r66i,r67r,r67i,r68r,r68i,r69r,r69i
		,r70r,r70i,r71r,r71i,r72r,r72i,r73r,r73i,r74r,r74i,r75r,r75i,r76r,r76i,r77r,r77i,r78r,r78i,r79r,r79i
		,r80r,r80i,r81r,r81i,r82r,r82i,r83r,r83i,r84r,r84i,r85r,r85i,r86r,r86i,r87r,r87i,r88r,r88i,r89r,r89i
		,r90r,r90i,r91r,r91i,r92r,r92i,r93r,r93i,r94r,r94i,r95r,r95i,r96r,r96i,r97r,r97i,r98r,r98i,r99r,r99i
		,rA0r,rA0i,rA1r,rA1i,rA2r,rA2i,rA3r,rA3i;

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
		p58 = p50 + p08;
		p60 = p58 + p08;

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
		p58 = p58 + ( (p58 >> DAT_BITS) << PAD_BITS );
		p60 = p60 + ( (p60 >> DAT_BITS) << PAD_BITS );
	}

/*...The radix-104 pass is here.	*/

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
	Twiddleless version arranges 8 sets of radix-13 DFT inputs as follows:
	0 in upper left corner, decrement 8 horizontally and 13 vertically [rep'd in hex here to match p-indexing]:
	[Can get these by running test_fft_radix with TTYPE = 0]:
		00,60,58,50,48,40,38,30,28,20,18,10,08   00,60,58,50,48,40,38,30,28,20,18,10,08 + p00
		5b,53,4b,43,3b,33,2b,23,1b,13,0b,03,63   58,50,48,40,38,30,28,20,18,10,08,00,60 + p03
		4e,46,3e,36,2e,26,1e,16,0e,06,66,5e,56   48,40,38,30,28,20,18,10,08,00,60,58,50 + p06
		41,39,31,29,21,19,11,09,01,61,59,51,49 = 40,38,30,28,20,18,10,08,00,60,58,50,48 + p01
		34,2c,24,1c,14,0c,04,64,5c,54,4c,44,3c   30,28,20,18,10,08,00,60,58,50,48,40,38 + p04
		27,1f,17,0f,07,67,5f,57,4f,47,3f,37,2f   20,18,10,08,00,60,58,50,48,40,38,30,28 + p07
		1a,12,0a,02,62,5a,52,4a,42,3a,32,2a,22   18,10,08,00,60,58,50,48,40,38,30,28,20 + p02
		0d,05,65,5d,55,4d,45,3d,35,2d,25,1d,15   08,00,60,58,50,48,40,38,30,28,20,18,10 + p05
	*/
	/*...gather the needed data (104 64-bit complex) and do 8 radix-9 transforms...*/
		RADIX_13_DFT(a[j1    ],a[j2    ],a[j1+p60],a[j2+p60],a[j1+p58],a[j2+p58],a[j1+p50],a[j2+p50],a[j1+p48],a[j2+p48],a[j1+p40],a[j2+p40],a[j1+p38],a[j2+p38],a[j1+p30],a[j2+p30],a[j1+p28],a[j2+p28],a[j1+p20],a[j2+p20],a[j1+p18],a[j2+p18],a[j1+p10],a[j2+p10],a[j1+p08],a[j2+p08], r00r,r00i,r01r,r01i,r02r,r02i,r03r,r03i,r04r,r04i,r05r,r05i,r06r,r06i,r07r,r07i,r08r,r08i,r09r,r09i,r10r,r10i,r11r,r11i,r12r,r12i);	jt = j1+p03; jp = j2+p03;
		RADIX_13_DFT(a[jt+p58],a[jp+p58],a[jt+p50],a[jp+p50],a[jt+p48],a[jp+p48],a[jt+p40],a[jp+p40],a[jt+p38],a[jp+p38],a[jt+p30],a[jp+p30],a[jt+p28],a[jp+p28],a[jt+p20],a[jp+p20],a[jt+p18],a[jp+p18],a[jt+p10],a[jp+p10],a[jt+p08],a[jp+p08],a[jt    ],a[jp    ],a[jt+p60],a[jp+p60], r13r,r13i,r14r,r14i,r15r,r15i,r16r,r16i,r17r,r17i,r18r,r18i,r19r,r19i,r20r,r20i,r21r,r21i,r22r,r22i,r23r,r23i,r24r,r24i,r25r,r25i);	jt = j1+p06; jp = j2+p06;
		RADIX_13_DFT(a[jt+p48],a[jp+p48],a[jt+p40],a[jp+p40],a[jt+p38],a[jp+p38],a[jt+p30],a[jp+p30],a[jt+p28],a[jp+p28],a[jt+p20],a[jp+p20],a[jt+p18],a[jp+p18],a[jt+p10],a[jp+p10],a[jt+p08],a[jp+p08],a[jt    ],a[jp    ],a[jt+p60],a[jp+p60],a[jt+p58],a[jp+p58],a[jt+p50],a[jp+p50], r26r,r26i,r27r,r27i,r28r,r28i,r29r,r29i,r30r,r30i,r31r,r31i,r32r,r32i,r33r,r33i,r34r,r34i,r35r,r35i,r36r,r36i,r37r,r37i,r38r,r38i);	jt = j1+p01; jp = j2+p01;
		RADIX_13_DFT(a[jt+p40],a[jp+p40],a[jt+p38],a[jp+p38],a[jt+p30],a[jp+p30],a[jt+p28],a[jp+p28],a[jt+p20],a[jp+p20],a[jt+p18],a[jp+p18],a[jt+p10],a[jp+p10],a[jt+p08],a[jp+p08],a[jt    ],a[jp    ],a[jt+p60],a[jp+p60],a[jt+p58],a[jp+p58],a[jt+p50],a[jp+p50],a[jt+p48],a[jp+p48], r39r,r39i,r40r,r40i,r41r,r41i,r42r,r42i,r43r,r43i,r44r,r44i,r45r,r45i,r46r,r46i,r47r,r47i,r48r,r48i,r49r,r49i,r50r,r50i,r51r,r51i);	jt = j1+p04; jp = j2+p04;
		RADIX_13_DFT(a[jt+p30],a[jp+p30],a[jt+p28],a[jp+p28],a[jt+p20],a[jp+p20],a[jt+p18],a[jp+p18],a[jt+p10],a[jp+p10],a[jt+p08],a[jp+p08],a[jt    ],a[jp    ],a[jt+p60],a[jp+p60],a[jt+p58],a[jp+p58],a[jt+p50],a[jp+p50],a[jt+p48],a[jp+p48],a[jt+p40],a[jp+p40],a[jt+p38],a[jp+p38], r52r,r52i,r53r,r53i,r54r,r54i,r55r,r55i,r56r,r56i,r57r,r57i,r58r,r58i,r59r,r59i,r60r,r60i,r61r,r61i,r62r,r62i,r63r,r63i,r64r,r64i);	jt = j1+p07; jp = j2+p07;
		RADIX_13_DFT(a[jt+p20],a[jp+p20],a[jt+p18],a[jp+p18],a[jt+p10],a[jp+p10],a[jt+p08],a[jp+p08],a[jt    ],a[jp    ],a[jt+p60],a[jp+p60],a[jt+p58],a[jp+p58],a[jt+p50],a[jp+p50],a[jt+p48],a[jp+p48],a[jt+p40],a[jp+p40],a[jt+p38],a[jp+p38],a[jt+p30],a[jp+p30],a[jt+p28],a[jp+p28], r65r,r65i,r66r,r66i,r67r,r67i,r68r,r68i,r69r,r69i,r70r,r70i,r71r,r71i,r72r,r72i,r73r,r73i,r74r,r74i,r75r,r75i,r76r,r76i,r77r,r77i);	jt = j1+p02; jp = j2+p02;
		RADIX_13_DFT(a[jt+p18],a[jp+p18],a[jt+p10],a[jp+p10],a[jt+p08],a[jp+p08],a[jt    ],a[jp    ],a[jt+p60],a[jp+p60],a[jt+p58],a[jp+p58],a[jt+p50],a[jp+p50],a[jt+p48],a[jp+p48],a[jt+p40],a[jp+p40],a[jt+p38],a[jp+p38],a[jt+p30],a[jp+p30],a[jt+p28],a[jp+p28],a[jt+p20],a[jp+p20], r78r,r78i,r79r,r79i,r80r,r80i,r81r,r81i,r82r,r82i,r83r,r83i,r84r,r84i,r85r,r85i,r86r,r86i,r87r,r87i,r88r,r88i,r89r,r89i,r90r,r90i);	jt = j1+p05; jp = j2+p05;
		RADIX_13_DFT(a[jt+p08],a[jp+p08],a[jt    ],a[jp    ],a[jt+p60],a[jp+p60],a[jt+p58],a[jp+p58],a[jt+p50],a[jp+p50],a[jt+p48],a[jp+p48],a[jt+p40],a[jp+p40],a[jt+p38],a[jp+p38],a[jt+p30],a[jp+p30],a[jt+p28],a[jp+p28],a[jt+p20],a[jp+p20],a[jt+p18],a[jp+p18],a[jt+p10],a[jp+p10], r91r,r91i,r92r,r92i,r93r,r93i,r94r,r94i,r95r,r95i,r96r,r96i,r97r,r97i,r98r,r98i,r99r,r99i,rA0r,rA0i,rA1r,rA1i,rA2r,rA2i,rA3r,rA3i);

	/*...and now do 13 radix-8 transforms.
	Use the supercalafragalistic Ancient Chinese Secret index-munging formula [SACSIMPF]
	to properly permute the outputs of the ensuing 13 radix-8 DFTs to the required ordering, which in terms of our 13 DFTs is

			[these indices are base-10]				[hex]
		  0,  1,  3,  2,  6,  7,  4,  5   01326745 + p00
		 97, 96, 98, 99,103,102,101,100   10237654 + p60
		 92, 93, 95, 94, 89, 88, 90, 91   45761023 + p58
		 83, 82, 81, 80, 84, 85, 87, 86   32104576 + p50
		 78, 79, 76, 77, 75, 74, 73, 72   67453210 + p48
		 64, 65, 67, 66, 70, 71, 68, 69   01326745 + p40
		 61, 60, 62, 63, 56, 57, 59, 58   54670132 + p38
		 50, 51, 48, 49, 53, 52, 54, 55 = 23015467 + p30
		 47, 46, 45, 44, 42, 43, 40, 41   76542301 + p28
		 33, 32, 34, 35, 39, 38, 37, 36   10237654 + p20
		 28, 29, 31, 30, 25, 24, 26, 27   45761023 + p18
		 19, 18, 17, 16, 20, 21, 23, 22   32104576 + p10
		 14, 15, 12, 13, 11, 10,  9,  8   67453210 + p08
	*/
					 /*                                    inputs                                 */ /*                       intermediates                       */ /*                 outputs                   */
		RADIX_08_DIF(r00r,r00i,r13r,r13i,r26r,r26i,r39r,r39i,r52r,r52i,r65r,r65i,r78r,r78i,r91r,r91i, x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i, a[j1    ],a[j2    ],a[j1+p01],a[j2+p01],a[j1+p03],a[j2+p03],a[j1+p02],a[j2+p02],a[j1+p06],a[j2+p06],a[j1+p07],a[j2+p07],a[j1+p04],a[j2+p04],a[j1+p05],a[j2+p05], rt,it);	jt = j1+p60; jp = j2+p60;
		RADIX_08_DIF(r01r,r01i,r14r,r14i,r27r,r27i,r40r,r40i,r53r,r53i,r66r,r66i,r79r,r79i,r92r,r92i, x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i, a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04], rt,it);	jt = j1+p58; jp = j2+p58;
		RADIX_08_DIF(r02r,r02i,r15r,r15i,r28r,r28i,r41r,r41i,r54r,r54i,r67r,r67i,r80r,r80i,r93r,r93i, x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i, a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03], rt,it);	jt = j1+p50; jp = j2+p50;
		RADIX_08_DIF(r03r,r03i,r16r,r16i,r29r,r29i,r42r,r42i,r55r,r55i,r68r,r68i,r81r,r81i,r94r,r94i, x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i, a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06], rt,it);	jt = j1+p48; jp = j2+p48;
		RADIX_08_DIF(r04r,r04i,r17r,r17i,r30r,r30i,r43r,r43i,r56r,r56i,r69r,r69i,r82r,r82i,r95r,r95i, x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i, a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ], rt,it);	jt = j1+p40; jp = j2+p40;
		RADIX_08_DIF(r05r,r05i,r18r,r18i,r31r,r31i,r44r,r44i,r57r,r57i,r70r,r70i,r83r,r83i,r96r,r96i, x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i, a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05], rt,it);	jt = j1+p38; jp = j2+p38;
		RADIX_08_DIF(r06r,r06i,r19r,r19i,r32r,r32i,r45r,r45i,r58r,r58i,r71r,r71i,r84r,r84i,r97r,r97i, x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i, a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02], rt,it);	jt = j1+p30; jp = j2+p30;
		RADIX_08_DIF(r07r,r07i,r20r,r20i,r33r,r33i,r46r,r46i,r59r,r59i,r72r,r72i,r85r,r85i,r98r,r98i, x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i, a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07], rt,it);	jt = j1+p28; jp = j2+p28;
		RADIX_08_DIF(r08r,r08i,r21r,r21i,r34r,r34i,r47r,r47i,r60r,r60i,r73r,r73i,r86r,r86i,r99r,r99i, x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i, a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01], rt,it);	jt = j1+p20; jp = j2+p20;
		RADIX_08_DIF(r09r,r09i,r22r,r22i,r35r,r35i,r48r,r48i,r61r,r61i,r74r,r74i,r87r,r87i,rA0r,rA0i, x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i, a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04], rt,it);	jt = j1+p18; jp = j2+p18;
		RADIX_08_DIF(r10r,r10i,r23r,r23i,r36r,r36i,r49r,r49i,r62r,r62i,r75r,r75i,r88r,r88i,rA1r,rA1i, x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i, a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03], rt,it);	jt = j1+p10; jp = j2+p10;
		RADIX_08_DIF(r11r,r11i,r24r,r24i,r37r,r37i,r50r,r50i,r63r,r63i,r76r,r76i,r89r,r89i,rA2r,rA2i, x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i, a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06], rt,it);	jt = j1+p08; jp = j2+p08;
		RADIX_08_DIF(r12r,r12i,r25r,r25i,r38r,r38i,r51r,r51i,r64r,r64i,r77r,r77i,r90r,r90i,rA3r,rA3i, x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i, a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ], rt,it);
	}
}

/***************/

void radix104_dit_pass1(double a[], int n)
{
/*
!...Acronym: DIF = Decimation In Frequency
!
!...Subroutine to perform an initial radix-104 complex DIT FFT pass on the data in the length-N real vector A.
!
!   See the documentation in radix16_dif_pass for further details on storage and indexing.
!
!   See the documentation in radix7_dif_pass for details on the radix-13 subtransforms.
*/
	int j,j1,j2,jp,jt;
	// In order to preserve length-2 numeric-index-offset property here, use hex for digit:
	static int NDIVR,p01,p02,p03,p04,p05,p06,p07,p08,p10,p18,p20,p28,p30,p38,p40,p48,p50,p58,p60, first_entry=TRUE;
	double rt,it, x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i
		,r00r,r00i,r01r,r01i,r02r,r02i,r03r,r03i,r04r,r04i,r05r,r05i,r06r,r06i,r07r,r07i,r08r,r08i,r09r,r09i
		,r10r,r10i,r11r,r11i,r12r,r12i,r13r,r13i,r14r,r14i,r15r,r15i,r16r,r16i,r17r,r17i,r18r,r18i,r19r,r19i
		,r20r,r20i,r21r,r21i,r22r,r22i,r23r,r23i,r24r,r24i,r25r,r25i,r26r,r26i,r27r,r27i,r28r,r28i,r29r,r29i
		,r30r,r30i,r31r,r31i,r32r,r32i,r33r,r33i,r34r,r34i,r35r,r35i,r36r,r36i,r37r,r37i,r38r,r38i,r39r,r39i
		,r40r,r40i,r41r,r41i,r42r,r42i,r43r,r43i,r44r,r44i,r45r,r45i,r46r,r46i,r47r,r47i,r48r,r48i,r49r,r49i
		,r50r,r50i,r51r,r51i,r52r,r52i,r53r,r53i,r54r,r54i,r55r,r55i,r56r,r56i,r57r,r57i,r58r,r58i,r59r,r59i
		,r60r,r60i,r61r,r61i,r62r,r62i,r63r,r63i,r64r,r64i,r65r,r65i,r66r,r66i,r67r,r67i,r68r,r68i,r69r,r69i
		,r70r,r70i,r71r,r71i,r72r,r72i,r73r,r73i,r74r,r74i,r75r,r75i,r76r,r76i,r77r,r77i,r78r,r78i,r79r,r79i
		,r80r,r80i,r81r,r81i,r82r,r82i,r83r,r83i,r84r,r84i,r85r,r85i,r86r,r86i,r87r,r87i,r88r,r88i,r89r,r89i
		,r90r,r90i,r91r,r91i,r92r,r92i,r93r,r93i,r94r,r94i,r95r,r95i,r96r,r96i,r97r,r97i,r98r,r98i,r99r,r99i
		,rA0r,rA0i,rA1r,rA1i,rA2r,rA2i,rA3r,rA3i;

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
		p58 = p50 + p08;
		p60 = p58 + p08;

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
		p58 = p58 + ( (p58 >> DAT_BITS) << PAD_BITS );
		p60 = p60 + ( (p60 >> DAT_BITS) << PAD_BITS );
	}

/*...The radix-104 pass is here.	*/

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

		00,01,02,03,04,05,06,07,08,09,0a,0b,0c,0d,0e,0f,10,11,12,13,14,15,16,17,18,19,
		1a,1b,1c,1d,1e,1f,20,21,22,23,24,25,26,27,28,29,2a,2b,2c,2d,2e,2f,30,31,32,33,
		34,35,36,37,38,39,3a,3b,3c,3d,3e,3f,40,41,42,43,44,45,46,47,48,49,4a,4b,4c,4d,
		4e,4f,50,51,52,53,54,55,56,57,58,59,5a,5b,5c,5d,5e,5f,60,61,62,63,64,65,66,67.
	map index-by-index [e.g. the index 01 gets replaced by 0x60, not the index at *position*...] to
		00,60,58,50,48,40,38,30,28,20,18,10,08,5b,53,4b,43,3b,33,2b,23,1b,13,0b,03,63,
		4e,46,3e,36,2e,26,1e,16,0e,06,66,5e,56,41,39,31,29,21,19,11,09,01,61,59,51,49,
		34,2c,24,1c,14,0c,04,64,5c,54,4c,44,3c,27,1f,17,0f,07,67,5f,57,4f,47,3f,37,2f,
		1a,12,0a,02,62,5a,52,4a,42,3a,32,2a,22,0d,05,65,5d,55,4d,45,3d,35,2d,25,1d,15.	(*)

	(***NOTE*** The following set of permutations can be auto-generated by running test_dft_radix()
	for the radix in question in TEST_TYPE = 1 [DIT] mode, skipping the actual DIT-pass step initially.)

	Remember, inputs to DIT are bit-reversed, so [all indices hex here]

		a[00,01,02,03,04,05,06,07,08,09,0a,0b,0c,0d,0e,0f,10,11,12,13,14,15,16,17,18,19,1a,1b,1c,1d,1e,1f,20,21,22,23,24,25,26,27,28,29,2a,2b,2c,2d,2e,2f,30,31,32,33,34,35,36,37,38,39,3a,3b,3c,3d,3e,3f,40,41,42,43,44,45,46,47,48,49,4a,4b,4c,4d,4e,4f,50,51,52,53,54,55,56,57,58,59,5a,5b,5c,5d,5e,5f,60,61,62,63,64,65,66,67] contain [cf. bit-reversal index array output by test_fft_radix()]
		x[00,34,1a,4e,0d,41,27,5b,01,35,1b,4f,0e,42,28,5c,02,36,1c,50,0f,43,29,5d,03,37,1d,51,10,44,2a,5e,04,38,1e,52,11,45,2b,5f,05,39,1f,53,12,46,2c,60,06,3a,20,54,13,47,2d,61,07,3b,21,55,14,48,2e,62,08,3c,22,56,15,49,2f,63,09,3d,23,57,16,4a,30,64,0a,3e,24,58,17,4b,31,65,0b,3f,25,59,18,4c,32,66,0c,40,26,5a,19,4d,33,67], which get swapped [using the permutation (*)] to
		x[00,34,4e,1a,5b,27,41,0d,60,2c,46,12,53,1f,39,05,58,24,3e,0a,4b,17,31,65,50,1c,36,02,43,0f,29,5d,48,14,2e,62,3b,07,21,55,40,0c,26,5a,33,67,19,4d,38,04,1e,52,2b,5f,11,45,30,64,16,4a,23,57,09,3d,28,5c,0e,42,1b,4f,01,35,20,54,06,3a,13,47,61,2d,18,4c,66,32,0b,3f,59,25,10,44,5e,2a,03,37,51,1d,08,3c,56,22,63,2f,49,15], which means the a-indices get swapped as
		a[00,01,03,02,07,06,05,04,2f,2e,2d,2c,2b,2a,29,28,53,52,51,50,55,54,56,57,13,12,11,10,15,14,16,17,3d,3c,3e,3f,39,38,3a,3b,61,60,62,63,66,67,64,65,21,20,22,23,26,27,24,25,4e,4f,4c,4d,4a,4b,48,49,0e,0f,0c,0d,0a,0b,08,09,32,33,30,31,34,35,37,36,5c,5d,5f,5e,58,59,5b,5a,1c,1d,1f,1e,18,19,1b,1a,40,41,43,42,47,46,45,44].

	Breaking the final index-row into octets are the 13 octets going into the radix-8 DFTs:

					[all indices hex here]
		00,01,03,02,07,06,05,04   01327654 + p00
		2f,2e,2d,2c,2b,2a,29,28   76543210 + p28
		53,52,51,50,55,54,56,57   32105467 + p50
		13,12,11,10,15,14,16,17   32105467 + p10
		3d,3c,3e,3f,39,38,3a,3b   54671023 + p38
		61,60,62,63,66,67,64,65   10236745 + p60
		21,20,22,23,26,27,24,25   10236745 + p20
		4e,4f,4c,4d,4a,4b,48,49 = 67452301 + p48
		0e,0f,0c,0d,0a,0b,08,09   67452301 + p08
		32,33,30,31,34,35,37,36   23014576 + p30
		5c,5d,5f,5e,58,59,5b,5a   45760132 + p58
		1c,1d,1f,1e,18,19,1b,1a   45760132 + p18
		40,41,43,42,47,46,45,44   01327654 + p40
*/
	/*...gather the needed data (104 64-bit complex) and do 13 radix-8 transforms...*/
		RADIX_08_DIT(a[j1    ],a[j2    ],a[j1+p01],a[j2+p01],a[j1+p03],a[j2+p03],a[j1+p02],a[j2+p02],a[j1+p07],a[j2+p07],a[j1+p06],a[j2+p06],a[j1+p05],a[j2+p05],a[j1+p04],a[j2+p04], x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i, r00r,r00i,r13r,r13i,r26r,r26i,r39r,r39i,r52r,r52i,r65r,r65i,r78r,r78i,r91r,r91i, rt,it);	jt = j1+p28; jp = j2+p28;
		RADIX_08_DIT(a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ], x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i, r01r,r01i,r14r,r14i,r27r,r27i,r40r,r40i,r53r,r53i,r66r,r66i,r79r,r79i,r92r,r92i, rt,it);	jt = j1+p50; jp = j2+p50;
		RADIX_08_DIT(a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07], x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i, r02r,r02i,r15r,r15i,r28r,r28i,r41r,r41i,r54r,r54i,r67r,r67i,r80r,r80i,r93r,r93i, rt,it);	jt = j1+p10; jp = j2+p10;
		RADIX_08_DIT(a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07], x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i, r03r,r03i,r16r,r16i,r29r,r29i,r42r,r42i,r55r,r55i,r68r,r68i,r81r,r81i,r94r,r94i, rt,it);	jt = j1+p38; jp = j2+p38;
		RADIX_08_DIT(a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03], x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i, r04r,r04i,r17r,r17i,r30r,r30i,r43r,r43i,r56r,r56i,r69r,r69i,r82r,r82i,r95r,r95i, rt,it);	jt = j1+p60; jp = j2+p60;
		RADIX_08_DIT(a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05], x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i, r05r,r05i,r18r,r18i,r31r,r31i,r44r,r44i,r57r,r57i,r70r,r70i,r83r,r83i,r96r,r96i, rt,it);	jt = j1+p20; jp = j2+p20;
		RADIX_08_DIT(a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05], x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i, r06r,r06i,r19r,r19i,r32r,r32i,r45r,r45i,r58r,r58i,r71r,r71i,r84r,r84i,r97r,r97i, rt,it);	jt = j1+p48; jp = j2+p48;
		RADIX_08_DIT(a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01], x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i, r07r,r07i,r20r,r20i,r33r,r33i,r46r,r46i,r59r,r59i,r72r,r72i,r85r,r85i,r98r,r98i, rt,it);	jt = j1+p08; jp = j2+p08;
		RADIX_08_DIT(a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01], x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i, r08r,r08i,r21r,r21i,r34r,r34i,r47r,r47i,r60r,r60i,r73r,r73i,r86r,r86i,r99r,r99i, rt,it);	jt = j1+p30; jp = j2+p30;
		RADIX_08_DIT(a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06], x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i, r09r,r09i,r22r,r22i,r35r,r35i,r48r,r48i,r61r,r61i,r74r,r74i,r87r,r87i,rA0r,rA0i, rt,it);	jt = j1+p58; jp = j2+p58;
		RADIX_08_DIT(a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02], x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i, r10r,r10i,r23r,r23i,r36r,r36i,r49r,r49i,r62r,r62i,r75r,r75i,r88r,r88i,rA1r,rA1i, rt,it);	jt = j1+p18; jp = j2+p18;
		RADIX_08_DIT(a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02], x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i, r11r,r11i,r24r,r24i,r37r,r37i,r50r,r50i,r63r,r63i,r76r,r76i,r89r,r89i,rA2r,rA2i, rt,it);	jt = j1+p40; jp = j2+p40;
		RADIX_08_DIT(a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04], x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i, r12r,r12i,r25r,r25i,r38r,r38i,r51r,r51i,r64r,r64i,r77r,r77i,r90r,r90i,rA3r,rA3i, rt,it);

	/*...and now do 8 radix-13 transforms: 

	Use the inverse supercalafragalistic Ancient Chinese Secret index-munging formula [iSACSIMPF]
	to properly permute the outputs of the 8 radix-13 DFTs:

								[all indices hex here]
		00,28,50,10,38,60,20,48,08,30,58,18,40   00,28,50,10,38,60,20,48,08,30,58,18,40
		41,01,29,51,11,39,61,21,49,09,31,59,19   40,00,28,50,10,38,60,20,48,08,30,58,18 + p01
		1a,42,02,2a,52,12,3a,62,22,4a,0a,32,5a   18,40,00,28,50,10,38,60,20,48,08,30,58 + p02
		5b,1b,43,03,2b,53,13,3b,63,23,4b,0b,33 = 58,18,40,00,28,50,10,38,60,20,48,08,30 + p03
		34,5c,1c,44,04,2c,54,14,3c,64,24,4c,0c   30,58,18,40,00,28,50,10,38,60,20,48,08 + p04
		0d,35,5d,1d,45,05,2d,55,15,3d,65,25,4d   08,30,58,18,40,00,28,50,10,38,60,20,48 + p05
		4e,0e,36,5e,1e,46,06,2e,56,16,3e,66,26   48,08,30,58,18,40,00,28,50,10,38,60,20 + p06
		27,4f,0f,37,5f,1f,47,07,2f,57,17,3f,67   20,48,08,30,58,18,40,00,28,50,10,38,60 + p07
*/
		RADIX_13_DFT(r00r,r00i,r01r,r01i,r02r,r02i,r03r,r03i,r04r,r04i,r05r,r05i,r06r,r06i,r07r,r07i,r08r,r08i,r09r,r09i,r10r,r10i,r11r,r11i,r12r,r12i, a[j1    ],a[j2    ],a[j1+p28],a[j2+p28],a[j1+p50],a[j2+p50],a[j1+p10],a[j2+p10],a[j1+p38],a[j2+p38],a[j1+p60],a[j2+p60],a[j1+p20],a[j2+p20],a[j1+p48],a[j2+p48],a[j1+p08],a[j2+p08],a[j1+p30],a[j2+p30],a[j1+p58],a[j2+p58],a[j1+p18],a[j2+p18],a[j1+p40],a[j2+p40]);	jt = j1+p01; jp = j2+p01;
		RADIX_13_DFT(r13r,r13i,r14r,r14i,r15r,r15i,r16r,r16i,r17r,r17i,r18r,r18i,r19r,r19i,r20r,r20i,r21r,r21i,r22r,r22i,r23r,r23i,r24r,r24i,r25r,r25i, a[jt+p40],a[jp+p40],a[jt    ],a[jp    ],a[jt+p28],a[jp+p28],a[jt+p50],a[jp+p50],a[jt+p10],a[jp+p10],a[jt+p38],a[jp+p38],a[jt+p60],a[jp+p60],a[jt+p20],a[jp+p20],a[jt+p48],a[jp+p48],a[jt+p08],a[jp+p08],a[jt+p30],a[jp+p30],a[jt+p58],a[jp+p58],a[jt+p18],a[jp+p18]);	jt = j1+p02; jp = j2+p02;
		RADIX_13_DFT(r26r,r26i,r27r,r27i,r28r,r28i,r29r,r29i,r30r,r30i,r31r,r31i,r32r,r32i,r33r,r33i,r34r,r34i,r35r,r35i,r36r,r36i,r37r,r37i,r38r,r38i, a[jt+p18],a[jp+p18],a[jt+p40],a[jp+p40],a[jt    ],a[jp    ],a[jt+p28],a[jp+p28],a[jt+p50],a[jp+p50],a[jt+p10],a[jp+p10],a[jt+p38],a[jp+p38],a[jt+p60],a[jp+p60],a[jt+p20],a[jp+p20],a[jt+p48],a[jp+p48],a[jt+p08],a[jp+p08],a[jt+p30],a[jp+p30],a[jt+p58],a[jp+p58]);	jt = j1+p03; jp = j2+p03;
		RADIX_13_DFT(r39r,r39i,r40r,r40i,r41r,r41i,r42r,r42i,r43r,r43i,r44r,r44i,r45r,r45i,r46r,r46i,r47r,r47i,r48r,r48i,r49r,r49i,r50r,r50i,r51r,r51i, a[jt+p58],a[jp+p58],a[jt+p18],a[jp+p18],a[jt+p40],a[jp+p40],a[jt    ],a[jp    ],a[jt+p28],a[jp+p28],a[jt+p50],a[jp+p50],a[jt+p10],a[jp+p10],a[jt+p38],a[jp+p38],a[jt+p60],a[jp+p60],a[jt+p20],a[jp+p20],a[jt+p48],a[jp+p48],a[jt+p08],a[jp+p08],a[jt+p30],a[jp+p30]);	jt = j1+p04; jp = j2+p04;
		RADIX_13_DFT(r52r,r52i,r53r,r53i,r54r,r54i,r55r,r55i,r56r,r56i,r57r,r57i,r58r,r58i,r59r,r59i,r60r,r60i,r61r,r61i,r62r,r62i,r63r,r63i,r64r,r64i, a[jt+p30],a[jp+p30],a[jt+p58],a[jp+p58],a[jt+p18],a[jp+p18],a[jt+p40],a[jp+p40],a[jt    ],a[jp    ],a[jt+p28],a[jp+p28],a[jt+p50],a[jp+p50],a[jt+p10],a[jp+p10],a[jt+p38],a[jp+p38],a[jt+p60],a[jp+p60],a[jt+p20],a[jp+p20],a[jt+p48],a[jp+p48],a[jt+p08],a[jp+p08]);	jt = j1+p05; jp = j2+p05;
		RADIX_13_DFT(r65r,r65i,r66r,r66i,r67r,r67i,r68r,r68i,r69r,r69i,r70r,r70i,r71r,r71i,r72r,r72i,r73r,r73i,r74r,r74i,r75r,r75i,r76r,r76i,r77r,r77i, a[jt+p08],a[jp+p08],a[jt+p30],a[jp+p30],a[jt+p58],a[jp+p58],a[jt+p18],a[jp+p18],a[jt+p40],a[jp+p40],a[jt    ],a[jp    ],a[jt+p28],a[jp+p28],a[jt+p50],a[jp+p50],a[jt+p10],a[jp+p10],a[jt+p38],a[jp+p38],a[jt+p60],a[jp+p60],a[jt+p20],a[jp+p20],a[jt+p48],a[jp+p48]);	jt = j1+p06; jp = j2+p06;
		RADIX_13_DFT(r78r,r78i,r79r,r79i,r80r,r80i,r81r,r81i,r82r,r82i,r83r,r83i,r84r,r84i,r85r,r85i,r86r,r86i,r87r,r87i,r88r,r88i,r89r,r89i,r90r,r90i, a[jt+p48],a[jp+p48],a[jt+p08],a[jp+p08],a[jt+p30],a[jp+p30],a[jt+p58],a[jp+p58],a[jt+p18],a[jp+p18],a[jt+p40],a[jp+p40],a[jt    ],a[jp    ],a[jt+p28],a[jp+p28],a[jt+p50],a[jp+p50],a[jt+p10],a[jp+p10],a[jt+p38],a[jp+p38],a[jt+p60],a[jp+p60],a[jt+p20],a[jp+p20]);	jt = j1+p07; jp = j2+p07;
		RADIX_13_DFT(r91r,r91i,r92r,r92i,r93r,r93i,r94r,r94i,r95r,r95i,r96r,r96i,r97r,r97i,r98r,r98i,r99r,r99i,rA0r,rA0i,rA1r,rA1i,rA2r,rA2i,rA3r,rA3i, a[jt+p20],a[jp+p20],a[jt+p48],a[jp+p48],a[jt+p08],a[jp+p08],a[jt+p30],a[jp+p30],a[jt+p58],a[jp+p58],a[jt+p18],a[jp+p18],a[jt+p40],a[jp+p40],a[jt    ],a[jp    ],a[jt+p28],a[jp+p28],a[jt+p50],a[jp+p50],a[jt+p10],a[jp+p10],a[jt+p38],a[jp+p38],a[jt+p60],a[jp+p60]);
	}
}

#undef RADIX
