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

#define RADIX 120	// Use #define rather than const int to ensure it's really a compile-time const in the C sense

/***************

int radix120_ditN_cy_dif1(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], double base[], double baseinv[], int iter, double *fracmax, uint64 p)
{
}

****************/

void radix120_dif_pass1(double a[], int n)
{
/*
!...Acronym: DIF = Decimation In Frequency
!
!...Subroutine to perform an initial radix-120 complex DIF FFT pass on the data in the length-N real vector A.
!
!   See the documentation in radix16_dif_pass for further details on storage and indexing.
!
!   See the documentation in radix5,16_dif_pass for details on the radix-5,16 subtransforms.
*/
	int j,j1,j2,jp,jt;
	// In order to preserve length-2 numeric-index-offset property here, use hex for digit:
	static int NDIVR,p01,p02,p03,p04,p05,p06,p07,p08,p10,p18,p20,p28,p30,p38,p40,p48,p50,p58,p60,p68,p70, first_entry=TRUE;
	const double	c3m1= -1.50000000000000000000,	/* cos(twopi/3)-1	*/
					s   =  0.86602540378443864675,	/* sin(twopi/3)		*/
					cn1 = -1.25000000000000000000,	/* [cos(u)+cos(2u)]/2-1 = -5/4 */
					cn2 =  0.55901699437494742409,	/* [cos(u)-cos(2u)]/2 */
					ss3 =  0.95105651629515357210,	/*  sin(u) */
					sn1 =  1.53884176858762670130,	/* [sin(u)+sin(2u)] */
					sn2 =  0.36327126400268044292;	/* [sin(u)-sin(2u)] */
	double rt,it
	,x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i
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
		,rA0r,rA0i,rA1r,rA1i,rA2r,rA2i,rA3r,rA3i,rA4r,rA4i,rA5r,rA5i,rA6r,rA6i,rA7r,rA7i,rA8r,rA8i,rA9r,rA9i
		,rB0r,rB0i,rB1r,rB1i,rB2r,rB2i,rB3r,rB3i,rB4r,rB4i,rB5r,rB5i,rB6r,rB6i,rB7r,rB7i,rB8r,rB8i,rB9r,rB9i;

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
		p68 = p60 + p08;
		p70 = p68 + p08;

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
		p68 = p68 + ( (p68 >> DAT_BITS) << PAD_BITS );
		p70 = p70 + ( (p70 >> DAT_BITS) << PAD_BITS );
	}

/*...The radix-120 pass is here.	*/

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
	Twiddleless version arranges 8 sets of radix-15 DFT inputs as follows:
	0 in upper left corner, decrement 8 horizontally and 15 vertically [rep'd in hex here to match p-indexing]:
	[Can get these by running test_fft_radix with TTYPE = 0]:

		00,70,68,60,58,50,48,40,38,30,28,20,18,10,08   00,70,68,60,58,50,48,40,38,30,28,20,18,10,08 + p00
		69,61,59,51,49,41,39,31,29,21,19,11,09,01,71   68,60,58,50,48,40,38,30,28,20,18,10,08,00,70 + p01
		5a,52,4a,42,3a,32,2a,22,1a,12,0a,02,72,6a,62   58,50,48,40,38,30,28,20,18,10,08,00,70,68,60 + p02
		4b,43,3b,33,2b,23,1b,13,0b,03,73,6b,63,5b,53 = 48,40,38,30,28,20,18,10,08,00,70,68,60,58,50 + p03
		3c,34,2c,24,1c,14,0c,04,74,6c,64,5c,54,4c,44   38,30,28,20,18,10,08,00,70,68,60,58,50,48,40 + p04
		2d,25,1d,15,0d,05,75,6d,65,5d,55,4d,45,3d,35   28,20,18,10,08,00,70,68,60,58,50,48,40,38,30 + p05
		1e,16,0e,06,76,6e,66,5e,56,4e,46,3e,36,2e,26   18,10,08,00,70,68,60,58,50,48,40,38,30,28,20 + p06
		0f,07,77,6f,67,5f,57,4f,47,3f,37,2f,27,1f,17   08,00,70,68,60,58,50,48,40,38,30,28,20,18,10 + p07
	*/
	/*...gather the needed data (120 64-bit complex) and do 8 radix-9 transforms...*/
		RADIX_15_DIF_B(s,c3m1,cn1,cn2,ss3,sn1,sn2, a[j1    ],a[j2    ],a[j1+p70],a[j2+p70],a[j1+p68],a[j2+p68],a[j1+p60],a[j2+p60],a[j1+p58],a[j2+p58],a[j1+p50],a[j2+p50],a[j1+p48],a[j2+p48],a[j1+p40],a[j2+p40],a[j1+p38],a[j2+p38],a[j1+p30],a[j2+p30],a[j1+p28],a[j2+p28],a[j1+p20],a[j2+p20],a[j1+p18],a[j2+p18],a[j1+p10],a[j2+p10],a[j1+p08],a[j2+p08], r00r,r00i,r01r,r01i,r02r,r02i,r03r,r03i,r04r,r04i,r05r,r05i,r06r,r06i,r07r,r07i,r08r,r08i,r09r,r09i,r10r,r10i,r11r,r11i,r12r,r12i,r13r,r13i,r14r,r14i);	jt = j1+p01; jp = j2+p01;
		RADIX_15_DIF_B(s,c3m1,cn1,cn2,ss3,sn1,sn2, a[jt+p68],a[jp+p68],a[jt+p60],a[jp+p60],a[jt+p58],a[jp+p58],a[jt+p50],a[jp+p50],a[jt+p48],a[jp+p48],a[jt+p40],a[jp+p40],a[jt+p38],a[jp+p38],a[jt+p30],a[jp+p30],a[jt+p28],a[jp+p28],a[jt+p20],a[jp+p20],a[jt+p18],a[jp+p18],a[jt+p10],a[jp+p10],a[jt+p08],a[jp+p08],a[jt    ],a[jp    ],a[jt+p70],a[jp+p70], r15r,r15i,r16r,r16i,r17r,r17i,r18r,r18i,r19r,r19i,r20r,r20i,r21r,r21i,r22r,r22i,r23r,r23i,r24r,r24i,r25r,r25i,r26r,r26i,r27r,r27i,r28r,r28i,r29r,r29i);	jt = j1+p02; jp = j2+p02;
		RADIX_15_DIF_B(s,c3m1,cn1,cn2,ss3,sn1,sn2, a[jt+p58],a[jp+p58],a[jt+p50],a[jp+p50],a[jt+p48],a[jp+p48],a[jt+p40],a[jp+p40],a[jt+p38],a[jp+p38],a[jt+p30],a[jp+p30],a[jt+p28],a[jp+p28],a[jt+p20],a[jp+p20],a[jt+p18],a[jp+p18],a[jt+p10],a[jp+p10],a[jt+p08],a[jp+p08],a[jt    ],a[jp    ],a[jt+p70],a[jp+p70],a[jt+p68],a[jp+p68],a[jt+p60],a[jp+p60], r30r,r30i,r31r,r31i,r32r,r32i,r33r,r33i,r34r,r34i,r35r,r35i,r36r,r36i,r37r,r37i,r38r,r38i,r39r,r39i,r40r,r40i,r41r,r41i,r42r,r42i,r43r,r43i,r44r,r44i);	jt = j1+p03; jp = j2+p03;
		RADIX_15_DIF_B(s,c3m1,cn1,cn2,ss3,sn1,sn2, a[jt+p48],a[jp+p48],a[jt+p40],a[jp+p40],a[jt+p38],a[jp+p38],a[jt+p30],a[jp+p30],a[jt+p28],a[jp+p28],a[jt+p20],a[jp+p20],a[jt+p18],a[jp+p18],a[jt+p10],a[jp+p10],a[jt+p08],a[jp+p08],a[jt    ],a[jp    ],a[jt+p70],a[jp+p70],a[jt+p68],a[jp+p68],a[jt+p60],a[jp+p60],a[jt+p58],a[jp+p58],a[jt+p50],a[jp+p50], r45r,r45i,r46r,r46i,r47r,r47i,r48r,r48i,r49r,r49i,r50r,r50i,r51r,r51i,r52r,r52i,r53r,r53i,r54r,r54i,r55r,r55i,r56r,r56i,r57r,r57i,r58r,r58i,r59r,r59i);	jt = j1+p04; jp = j2+p04;
		RADIX_15_DIF_B(s,c3m1,cn1,cn2,ss3,sn1,sn2, a[jt+p38],a[jp+p38],a[jt+p30],a[jp+p30],a[jt+p28],a[jp+p28],a[jt+p20],a[jp+p20],a[jt+p18],a[jp+p18],a[jt+p10],a[jp+p10],a[jt+p08],a[jp+p08],a[jt    ],a[jp    ],a[jt+p70],a[jp+p70],a[jt+p68],a[jp+p68],a[jt+p60],a[jp+p60],a[jt+p58],a[jp+p58],a[jt+p50],a[jp+p50],a[jt+p48],a[jp+p48],a[jt+p40],a[jp+p40], r60r,r60i,r61r,r61i,r62r,r62i,r63r,r63i,r64r,r64i,r65r,r65i,r66r,r66i,r67r,r67i,r68r,r68i,r69r,r69i,r70r,r70i,r71r,r71i,r72r,r72i,r73r,r73i,r74r,r74i);	jt = j1+p05; jp = j2+p05;
		RADIX_15_DIF_B(s,c3m1,cn1,cn2,ss3,sn1,sn2, a[jt+p28],a[jp+p28],a[jt+p20],a[jp+p20],a[jt+p18],a[jp+p18],a[jt+p10],a[jp+p10],a[jt+p08],a[jp+p08],a[jt    ],a[jp    ],a[jt+p70],a[jp+p70],a[jt+p68],a[jp+p68],a[jt+p60],a[jp+p60],a[jt+p58],a[jp+p58],a[jt+p50],a[jp+p50],a[jt+p48],a[jp+p48],a[jt+p40],a[jp+p40],a[jt+p38],a[jp+p38],a[jt+p30],a[jp+p30], r75r,r75i,r76r,r76i,r77r,r77i,r78r,r78i,r79r,r79i,r80r,r80i,r81r,r81i,r82r,r82i,r83r,r83i,r84r,r84i,r85r,r85i,r86r,r86i,r87r,r87i,r88r,r88i,r89r,r89i);	jt = j1+p06; jp = j2+p06;
		RADIX_15_DIF_B(s,c3m1,cn1,cn2,ss3,sn1,sn2, a[jt+p18],a[jp+p18],a[jt+p10],a[jp+p10],a[jt+p08],a[jp+p08],a[jt    ],a[jp    ],a[jt+p70],a[jp+p70],a[jt+p68],a[jp+p68],a[jt+p60],a[jp+p60],a[jt+p58],a[jp+p58],a[jt+p50],a[jp+p50],a[jt+p48],a[jp+p48],a[jt+p40],a[jp+p40],a[jt+p38],a[jp+p38],a[jt+p30],a[jp+p30],a[jt+p28],a[jp+p28],a[jt+p20],a[jp+p20], r90r,r90i,r91r,r91i,r92r,r92i,r93r,r93i,r94r,r94i,r95r,r95i,r96r,r96i,r97r,r97i,r98r,r98i,r99r,r99i,rA0r,rA0i,rA1r,rA1i,rA2r,rA2i,rA3r,rA3i,rA4r,rA4i);	jt = j1+p07; jp = j2+p07;
		RADIX_15_DIF_B(s,c3m1,cn1,cn2,ss3,sn1,sn2, a[jt+p08],a[jp+p08],a[jt    ],a[jp    ],a[jt+p70],a[jp+p70],a[jt+p68],a[jp+p68],a[jt+p60],a[jp+p60],a[jt+p58],a[jp+p58],a[jt+p50],a[jp+p50],a[jt+p48],a[jp+p48],a[jt+p40],a[jp+p40],a[jt+p38],a[jp+p38],a[jt+p30],a[jp+p30],a[jt+p28],a[jp+p28],a[jt+p20],a[jp+p20],a[jt+p18],a[jp+p18],a[jt+p10],a[jp+p10], rA5r,rA5i,rA6r,rA6i,rA7r,rA7i,rA8r,rA8i,rA9r,rA9i,rB0r,rB0i,rB1r,rB1i,rB2r,rB2i,rB3r,rB3i,rB4r,rB4i,rB5r,rB5i,rB6r,rB6i,rB7r,rB7i,rB8r,rB8i,rB9r,rB9i);

	/*...and now do 15 radix-8 transforms.
	Use the supercalafragalistic Ancient Chinese Secret index-munging formula [SACSIMPF]
	to properly permute the outputs of the ensuing 15 radix-8 DFTs to the required ordering, which in terms of our 15 DFTs is

			[these indices are base-10]				[hex]
		  0,  1,  2,  3,  4,  5,  6,  7   01234567 + p00
		 18, 19, 17, 16, 22, 23, 21, 20   23106754 + p10
		 13, 12, 15, 14, 11, 10,  8,  9   54763201 + p08
		115,114,112,113,119,118,116,117   32017645 + p70
		108,109,110,111,106,107,105,104   45672310 + p68
		 97, 96, 99, 98,101,100,103,102   10325476 + p60
		 93, 92, 95, 94, 91, 90, 88, 89   54763201 + p58
		 80, 81, 82, 83, 84, 85, 86, 87 = 01234567 + p50
		 78, 79, 77, 76, 73, 72, 75, 74   67541032 + p48
		 65, 64, 67, 66, 69, 68, 71, 70   10325476 + p40
		 63, 62, 60, 61, 56, 57, 58, 59   76450123 + p38
		 50, 51, 49, 48, 54, 55, 53, 52   23106754 + p30
		 46, 47, 45, 44, 41, 40, 43, 42   67541032 + p28
		 35, 34, 32, 33, 39, 38, 36, 37   32017645 + p20
		 28, 29, 30, 31, 26, 27, 25, 24   45672310 + p18
	*/
					 /*                                    inputs                                 */ /*                       intermediates                       */ /*                 outputs                   */
		RADIX_08_DIF(r00r,r00i,r15r,r15i,r30r,r30i,r45r,r45i,r60r,r60i,r75r,r75i,r90r,r90i,rA5r,rA5i, x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i, a[j1    ],a[j2    ],a[j1+p01],a[j2+p01],a[j1+p02],a[j2+p02],a[j1+p03],a[j2+p03],a[j1+p04],a[j2+p04],a[j1+p05],a[j2+p05],a[j1+p06],a[j2+p06],a[j1+p07],a[j2+p07], rt,it);	jt = j1+p10; jp = j2+p10;
		RADIX_08_DIF(r01r,r01i,r16r,r16i,r31r,r31i,r46r,r46i,r61r,r61i,r76r,r76i,r91r,r91i,rA6r,rA6i, x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i, a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04], rt,it);	jt = j1+p08; jp = j2+p08;
		RADIX_08_DIF(r02r,r02i,r17r,r17i,r32r,r32i,r47r,r47i,r62r,r62i,r77r,r77i,r92r,r92i,rA7r,rA7i, x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i, a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01], rt,it);	jt = j1+p70; jp = j2+p70;
		RADIX_08_DIF(r03r,r03i,r18r,r18i,r33r,r33i,r48r,r48i,r63r,r63i,r78r,r78i,r93r,r93i,rA8r,rA8i, x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i, a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05], rt,it);	jt = j1+p68; jp = j2+p68;
		RADIX_08_DIF(r04r,r04i,r19r,r19i,r34r,r34i,r49r,r49i,r64r,r64i,r79r,r79i,r94r,r94i,rA9r,rA9i, x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i, a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ], rt,it);	jt = j1+p60; jp = j2+p60;
		RADIX_08_DIF(r05r,r05i,r20r,r20i,r35r,r35i,r50r,r50i,r65r,r65i,r80r,r80i,r95r,r95i,rB0r,rB0i, x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i, a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06], rt,it);	jt = j1+p58; jp = j2+p58;
		RADIX_08_DIF(r06r,r06i,r21r,r21i,r36r,r36i,r51r,r51i,r66r,r66i,r81r,r81i,r96r,r96i,rB1r,rB1i, x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i, a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01], rt,it);	jt = j1+p50; jp = j2+p50;
		RADIX_08_DIF(r07r,r07i,r22r,r22i,r37r,r37i,r52r,r52i,r67r,r67i,r82r,r82i,r97r,r97i,rB2r,rB2i, x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i, a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07], rt,it);	jt = j1+p48; jp = j2+p48;
		RADIX_08_DIF(r08r,r08i,r23r,r23i,r38r,r38i,r53r,r53i,r68r,r68i,r83r,r83i,r98r,r98i,rB3r,rB3i, x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i, a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02], rt,it);	jt = j1+p40; jp = j2+p40;
		RADIX_08_DIF(r09r,r09i,r24r,r24i,r39r,r39i,r54r,r54i,r69r,r69i,r84r,r84i,r99r,r99i,rB4r,rB4i, x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i, a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06], rt,it);	jt = j1+p38; jp = j2+p38;
		RADIX_08_DIF(r10r,r10i,r25r,r25i,r40r,r40i,r55r,r55i,r70r,r70i,r85r,r85i,rA0r,rA0i,rB5r,rB5i, x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i, a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03], rt,it);	jt = j1+p30; jp = j2+p30;
		RADIX_08_DIF(r11r,r11i,r26r,r26i,r41r,r41i,r56r,r56i,r71r,r71i,r86r,r86i,rA1r,rA1i,rB6r,rB6i, x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i, a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04], rt,it);	jt = j1+p28; jp = j2+p28;
		RADIX_08_DIF(r12r,r12i,r27r,r27i,r42r,r42i,r57r,r57i,r72r,r72i,r87r,r87i,rA2r,rA2i,rB7r,rB7i, x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i, a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02], rt,it);	jt = j1+p20; jp = j2+p20;
		RADIX_08_DIF(r13r,r13i,r28r,r28i,r43r,r43i,r58r,r58i,r73r,r73i,r88r,r88i,rA3r,rA3i,rB8r,rB8i, x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i, a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05], rt,it);	jt = j1+p18; jp = j2+p18;
		RADIX_08_DIF(r14r,r14i,r29r,r29i,r44r,r44i,r59r,r59i,r74r,r74i,r89r,r89i,rA4r,rA4i,rB9r,rB9i, x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i, a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ], rt,it);
	}
}

/***************/

void radix120_dit_pass1(double a[], int n)
{
/*
!...Acronym: DIF = Decimation In Frequency
!
!...Subroutine to perform an initial radix-120 complex DIT FFT pass on the data in the length-N real vector A.
!
!   See the documentation in radix16_dif_pass for further details on storage and indexing.
!
!   See the documentation in radix7_dif_pass for details on the radix-15 subtransforms.
*/
	int j,j1,j2,jp,jt;
	// In order to preserve length-2 numeric-index-offset property here, use hex for digit:
	static int NDIVR,p01,p02,p03,p04,p05,p06,p07,p08,p10,p18,p20,p28,p30,p38,p40,p48,p50,p58,p60,p68,p70, first_entry=TRUE;
	const double	c3m1= -1.50000000000000000000,	/* cos(twopi/3)-1	*/
					s   =  0.86602540378443864675,	/* sin(twopi/3)		*/
					cn1 = -1.25000000000000000000,	/* [cos(u)+cos(2u)]/2-1 = -5/4 */
					cn2 =  0.55901699437494742409,	/* [cos(u)-cos(2u)]/2 */
					ss3 =  0.95105651629515357210,	/*  sin(u) */
					sn1 =  1.53884176858762670130,	/* [sin(u)+sin(2u)] */
					sn2 =  0.36327126400268044292;	/* [sin(u)-sin(2u)] */
	double rt,it
	,x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i
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
		,rA0r,rA0i,rA1r,rA1i,rA2r,rA2i,rA3r,rA3i,rA4r,rA4i,rA5r,rA5i,rA6r,rA6i,rA7r,rA7i,rA8r,rA8i,rA9r,rA9i
		,rB0r,rB0i,rB1r,rB1i,rB2r,rB2i,rB3r,rB3i,rB4r,rB4i,rB5r,rB5i,rB6r,rB6i,rB7r,rB7i,rB8r,rB8i,rB9r,rB9i;

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
		p68 = p60 + p08;
		p70 = p68 + p08;

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
		p68 = p68 + ( (p68 >> DAT_BITS) << PAD_BITS );
		p70 = p70 + ( (p70 >> DAT_BITS) << PAD_BITS );
	}

/*...The radix-120 pass is here.	*/

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

		00,01,02,03,04,05,06,07,08,09,0a,0b,0c,0d,0e,0f,10,11,12,13,14,15,16,17,18,19,1a,1b,1c,1d,
		1e,1f,20,21,22,23,24,25,26,27,28,29,2a,2b,2c,2d,2e,2f,30,31,32,33,34,35,36,37,38,39,3a,3b,
		3c,3d,3e,3f,40,41,42,43,44,45,46,47,48,49,4a,4b,4c,4d,4e,4f,50,51,52,53,54,55,56,57,58,59,
		5a,5b,5c,5d,5e,5f,60,61,62,63,64,65,66,67,68,69,6a,6b,6c,6d,6e,6f,70,71,72,73,74,75,76,77,
	map index-by-index [e.g. the index 01 gets replaced by 0x70, not the index at *position*...] to
		00,70,68,60,58,50,48,40,38,30,28,20,18,10,08,69,61,59,51,49,41,39,31,29,21,19,11,09,01,71,
		5a,52,4a,42,3a,32,2a,22,1a,12,0a,02,72,6a,62,4b,43,3b,33,2b,23,1b,13,0b,03,73,6b,63,5b,53,
		3c,34,2c,24,1c,14,0c,04,74,6c,64,5c,54,4c,44,2d,25,1d,15,0d,05,75,6d,65,5d,55,4d,45,3d,35,
		1e,16,0e,06,76,6e,66,5e,56,4e,46,3e,36,2e,26,0f,07,77,6f,67,5f,57,4f,47,3f,37,2f,27,1f,17.	(*)

	(***NOTE*** The following set of permutations can be auto-generated by running test_dft_radix()
	for the radix in question in TEST_TYPE = 1 [DIT] mode, skipping the actual DIT-pass step initially.)

	Remember, inputs to DIT are bit-reversed, so [all indices hex here]

		a[00,01,02,03,04,05,06,07,08,09,0a,0b,0c,0d,0e,0f,10,11,12,13,14,15,16,17,18,19,1a,1b,1c,1d,1e,1f,20,21,22,23,24,25,26,27,28,29,2a,2b,2c,2d,2e,2f,30,31,32,33,34,35,36,37,38,39,3a,3b,3c,3d,3e,3f,40,41,42,43,44,45,46,47,48,49,4a,4b,4c,4d,4e,4f,50,51,52,53,54,55,56,57,58,59,5a,5b,5c,5d,5e,5f,60,61,62,63,64,65,66,67,68,69,6a,6b,6c,6d,6e,6f,70,71,72,73,74,75,76,77] contain [cf. bit-reversal index array output by test_fft_radix()]
		x[00,3c,1e,5a,0f,4b,2d,69,05,41,23,5f,14,50,32,6e,0a,46,28,64,19,55,37,73,01,3d,1f,5b,10,4c,2e,6a,06,42,24,60,15,51,33,6f,0b,47,29,65,1a,56,38,74,02,3e,20,5c,11,4d,2f,6b,07,43,25,61,16,52,34,70,0c,48,2a,66,1b,57,39,75,03,3f,21,5d,12,4e,30,6c,08,44,26,62,17,53,35,71,0d,49,2b,67,1c,58,3a,76,04,40,22,5e,13,4f,31,6d,09,45,27,63,18,54,36,72,0e,4a,2c,68,1d,59,3b,77], which get swapped [using the permutation (*)] to
		x[00,3c,5a,1e,69,2d,4b,0f,50,14,32,6e,41,05,23,5f,28,64,0a,46,19,55,73,37,70,34,52,16,61,25,43,07,48,0c,2a,66,39,75,1b,57,20,5c,02,3e,11,4d,6b,2f,68,2c,4a,0e,59,1d,3b,77,40,04,22,5e,31,6d,13,4f,18,54,72,36,09,45,63,27,60,24,42,06,51,15,33,6f,38,74,1a,56,29,65,0b,47,10,4c,6a,2e,01,3d,5b,1f,58,1c,3a,76,49,0d,2b,67,30,6c,12,4e,21,5d,03,3f,08,44,62,26,71,35,53,17], which means the a-indices get swapped as
		a[00,01,03,02,07,06,05,04,0d,0c,0e,0f,09,08,0a,0b,12,13,10,11,14,15,17,16,3f,3e,3d,3c,3b,3a,39,38,41,40,42,43,46,47,44,45,32,33,30,31,34,35,37,36,73,72,71,70,75,74,76,77,61,60,62,63,66,67,64,65,6c,6d,6f,6e,68,69,6b,6a,23,22,21,20,25,24,26,27,2e,2f,2c,2d,2a,2b,28,29,1c,1d,1f,1e,18,19,1b,1a,5d,5c,5e,5f,59,58,5a,5b,4e,4f,4c,4d,4a,4b,48,49,50,51,53,52,57,56,55,54].

	Breaking the final index-row into octets are the 15 octets going into the radix-8 DFTs:

			[these indices are base-10]				[hex]
		  0,  1,  3,  2,  7,  6,  5,  4   01327654 + p00
		 13, 12, 14, 15,  9,  8, 10, 11   54671023 + p08
		 18, 19, 16, 17, 20, 21, 23, 22   23014576 + p10
		 63, 62, 61, 60, 59, 58, 57, 56   76543210 + p38
		 65, 64, 66, 67, 70, 71, 68, 69   10236745 + p40
		 50, 51, 48, 49, 52, 53, 55, 54   23014576 + p30
		115,114,113,112,117,116,118,119   32105467 + p70
		 97, 96, 98, 99,102,103,100,101 = 10236745 + p60
		108,109,111,110,104,105,107,106   45760132 + p68
		 35, 34, 33, 32, 37, 36, 38, 39   32105467 + p20
		 46, 47, 44, 45, 42, 43, 40, 41   67452301 + p28
		 28, 29, 31, 30, 24, 25, 27, 26   45760132 + p18
		 93, 92, 94, 95, 89, 88, 90, 91   54671023 + p58
		 78, 79, 76, 77, 74, 75, 72, 73   67452301 + p48
		 80, 81, 83, 82, 87, 86, 85, 84   01327654 + p50
*/
	/*...gather the needed data (120 64-bit complex) and do 15 radix-8 transforms...*/
		RADIX_08_DIT(a[j1    ],a[j2    ],a[j1+p01],a[j2+p01],a[j1+p03],a[j2+p03],a[j1+p02],a[j2+p02],a[j1+p07],a[j2+p07],a[j1+p06],a[j2+p06],a[j1+p05],a[j2+p05],a[j1+p04],a[j2+p04], x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i, r00r,r00i,r15r,r15i,r30r,r30i,r45r,r45i,r60r,r60i,r75r,r75i,r90r,r90i,rA5r,rA5i, rt,it);	jt = j1+p08; jp = j2+p08;
		RADIX_08_DIT(a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03], x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i, r01r,r01i,r16r,r16i,r31r,r31i,r46r,r46i,r61r,r61i,r76r,r76i,r91r,r91i,rA6r,rA6i, rt,it);	jt = j1+p10; jp = j2+p10;
		RADIX_08_DIT(a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06], x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i, r02r,r02i,r17r,r17i,r32r,r32i,r47r,r47i,r62r,r62i,r77r,r77i,r92r,r92i,rA7r,rA7i, rt,it);	jt = j1+p38; jp = j2+p38;
		RADIX_08_DIT(a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ], x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i, r03r,r03i,r18r,r18i,r33r,r33i,r48r,r48i,r63r,r63i,r78r,r78i,r93r,r93i,rA8r,rA8i, rt,it);	jt = j1+p40; jp = j2+p40;
		RADIX_08_DIT(a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05], x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i, r04r,r04i,r19r,r19i,r34r,r34i,r49r,r49i,r64r,r64i,r79r,r79i,r94r,r94i,rA9r,rA9i, rt,it);	jt = j1+p30; jp = j2+p30;
		RADIX_08_DIT(a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06], x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i, r05r,r05i,r20r,r20i,r35r,r35i,r50r,r50i,r65r,r65i,r80r,r80i,r95r,r95i,rB0r,rB0i, rt,it);	jt = j1+p70; jp = j2+p70;
		RADIX_08_DIT(a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07], x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i, r06r,r06i,r21r,r21i,r36r,r36i,r51r,r51i,r66r,r66i,r81r,r81i,r96r,r96i,rB1r,rB1i, rt,it);	jt = j1+p60; jp = j2+p60;
		RADIX_08_DIT(a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05], x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i, r07r,r07i,r22r,r22i,r37r,r37i,r52r,r52i,r67r,r67i,r82r,r82i,r97r,r97i,rB2r,rB2i, rt,it);	jt = j1+p68; jp = j2+p68;
		RADIX_08_DIT(a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02], x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i, r08r,r08i,r23r,r23i,r38r,r38i,r53r,r53i,r68r,r68i,r83r,r83i,r98r,r98i,rB3r,rB3i, rt,it);	jt = j1+p20; jp = j2+p20;
		RADIX_08_DIT(a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07], x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i, r09r,r09i,r24r,r24i,r39r,r39i,r54r,r54i,r69r,r69i,r84r,r84i,r99r,r99i,rB4r,rB4i, rt,it);	jt = j1+p28; jp = j2+p28;
		RADIX_08_DIT(a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01], x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i, r10r,r10i,r25r,r25i,r40r,r40i,r55r,r55i,r70r,r70i,r85r,r85i,rA0r,rA0i,rB5r,rB5i, rt,it);	jt = j1+p18; jp = j2+p18;
		RADIX_08_DIT(a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02], x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i, r11r,r11i,r26r,r26i,r41r,r41i,r56r,r56i,r71r,r71i,r86r,r86i,rA1r,rA1i,rB6r,rB6i, rt,it);	jt = j1+p58; jp = j2+p58;
		RADIX_08_DIT(a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03], x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i, r12r,r12i,r27r,r27i,r42r,r42i,r57r,r57i,r72r,r72i,r87r,r87i,rA2r,rA2i,rB7r,rB7i, rt,it);	jt = j1+p48; jp = j2+p48;
		RADIX_08_DIT(a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01], x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i, r13r,r13i,r28r,r28i,r43r,r43i,r58r,r58i,r73r,r73i,r88r,r88i,rA3r,rA3i,rB8r,rB8i, rt,it);	jt = j1+p50; jp = j2+p50;
		RADIX_08_DIT(a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04], x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i, r14r,r14i,r29r,r29i,r44r,r44i,r59r,r59i,r74r,r74i,r89r,r89i,rA4r,rA4i,rB9r,rB9i, rt,it);

	/*...and now do 8 radix-15 transforms: 

	Use the inverse supercalafragalistic Ancient Chinese Secret index-munging formula [iSACSIMPF]
	to properly permute the outputs of the 8 radix-15 DFTs:
								[all indices hex here]
		00,68,58,48,38,28,18,08,70,60,50,40,30,20,10   00,68,58,48,38,28,18,08,70,60,50,40,30,20,10
		69,59,49,39,29,19,09,71,61,51,41,31,21,11,01   68,58,48,38,28,18,08,70,60,50,40,30,20,10,00 + p01
		5a,4a,3a,2a,1a,0a,72,62,52,42,32,22,12,02,6a   58,48,38,28,18,08,70,60,50,40,30,20,10,00,68 + p02
		4b,3b,2b,1b,0b,73,63,53,43,33,23,13,03,6b,5b = 48,38,28,18,08,70,60,50,40,30,20,10,00,68,58 + p03
		3c,2c,1c,0c,74,64,54,44,34,24,14,04,6c,5c,4c   38,28,18,08,70,60,50,40,30,20,10,00,68,58,48 + p04
		2d,1d,0d,75,65,55,45,35,25,15,05,6d,5d,4d,3d   28,18,08,70,60,50,40,30,20,10,00,68,58,48,38 + p05
		1e,0e,76,66,56,46,36,26,16,06,6e,5e,4e,3e,2e   18,08,70,60,50,40,30,20,10,00,68,58,48,38,28 + p06
		0f,77,67,57,47,37,27,17,07,6f,5f,4f,3f,2f,1f   08,70,60,50,40,30,20,10,00,68,58,48,38,28,18 + p07
*/
		RADIX_15_DIT_B(s,c3m1,cn1,cn2,ss3,sn1,sn2, r00r,r00i,r01r,r01i,r02r,r02i,r03r,r03i,r04r,r04i,r05r,r05i,r06r,r06i,r07r,r07i,r08r,r08i,r09r,r09i,r10r,r10i,r11r,r11i,r12r,r12i,r13r,r13i,r14r,r14i, a[j1    ],a[j2    ],a[j1+p68],a[j2+p68],a[j1+p58],a[j2+p58],a[j1+p48],a[j2+p48],a[j1+p38],a[j2+p38],a[j1+p28],a[j2+p28],a[j1+p18],a[j2+p18],a[j1+p08],a[j2+p08],a[j1+p70],a[j2+p70],a[j1+p60],a[j2+p60],a[j1+p50],a[j2+p50],a[j1+p40],a[j2+p40],a[j1+p30],a[j2+p30],a[j1+p20],a[j2+p20],a[j1+p10],a[j2+p10]);	jt = j1+p01; jp = j2+p01;
		RADIX_15_DIT_B(s,c3m1,cn1,cn2,ss3,sn1,sn2, r15r,r15i,r16r,r16i,r17r,r17i,r18r,r18i,r19r,r19i,r20r,r20i,r21r,r21i,r22r,r22i,r23r,r23i,r24r,r24i,r25r,r25i,r26r,r26i,r27r,r27i,r28r,r28i,r29r,r29i, a[jt+p68],a[jp+p68],a[jt+p58],a[jp+p58],a[jt+p48],a[jp+p48],a[jt+p38],a[jp+p38],a[jt+p28],a[jp+p28],a[jt+p18],a[jp+p18],a[jt+p08],a[jp+p08],a[jt+p70],a[jp+p70],a[jt+p60],a[jp+p60],a[jt+p50],a[jp+p50],a[jt+p40],a[jp+p40],a[jt+p30],a[jp+p30],a[jt+p20],a[jp+p20],a[jt+p10],a[jp+p10],a[jt    ],a[jp    ]);	jt = j1+p02; jp = j2+p02;
		RADIX_15_DIT_B(s,c3m1,cn1,cn2,ss3,sn1,sn2, r30r,r30i,r31r,r31i,r32r,r32i,r33r,r33i,r34r,r34i,r35r,r35i,r36r,r36i,r37r,r37i,r38r,r38i,r39r,r39i,r40r,r40i,r41r,r41i,r42r,r42i,r43r,r43i,r44r,r44i, a[jt+p58],a[jp+p58],a[jt+p48],a[jp+p48],a[jt+p38],a[jp+p38],a[jt+p28],a[jp+p28],a[jt+p18],a[jp+p18],a[jt+p08],a[jp+p08],a[jt+p70],a[jp+p70],a[jt+p60],a[jp+p60],a[jt+p50],a[jp+p50],a[jt+p40],a[jp+p40],a[jt+p30],a[jp+p30],a[jt+p20],a[jp+p20],a[jt+p10],a[jp+p10],a[jt    ],a[jp    ],a[jt+p68],a[jp+p68]);	jt = j1+p03; jp = j2+p03;
		RADIX_15_DIT_B(s,c3m1,cn1,cn2,ss3,sn1,sn2, r45r,r45i,r46r,r46i,r47r,r47i,r48r,r48i,r49r,r49i,r50r,r50i,r51r,r51i,r52r,r52i,r53r,r53i,r54r,r54i,r55r,r55i,r56r,r56i,r57r,r57i,r58r,r58i,r59r,r59i, a[jt+p48],a[jp+p48],a[jt+p38],a[jp+p38],a[jt+p28],a[jp+p28],a[jt+p18],a[jp+p18],a[jt+p08],a[jp+p08],a[jt+p70],a[jp+p70],a[jt+p60],a[jp+p60],a[jt+p50],a[jp+p50],a[jt+p40],a[jp+p40],a[jt+p30],a[jp+p30],a[jt+p20],a[jp+p20],a[jt+p10],a[jp+p10],a[jt    ],a[jp    ],a[jt+p68],a[jp+p68],a[jt+p58],a[jp+p58]);	jt = j1+p04; jp = j2+p04;
		RADIX_15_DIT_B(s,c3m1,cn1,cn2,ss3,sn1,sn2, r60r,r60i,r61r,r61i,r62r,r62i,r63r,r63i,r64r,r64i,r65r,r65i,r66r,r66i,r67r,r67i,r68r,r68i,r69r,r69i,r70r,r70i,r71r,r71i,r72r,r72i,r73r,r73i,r74r,r74i, a[jt+p38],a[jp+p38],a[jt+p28],a[jp+p28],a[jt+p18],a[jp+p18],a[jt+p08],a[jp+p08],a[jt+p70],a[jp+p70],a[jt+p60],a[jp+p60],a[jt+p50],a[jp+p50],a[jt+p40],a[jp+p40],a[jt+p30],a[jp+p30],a[jt+p20],a[jp+p20],a[jt+p10],a[jp+p10],a[jt    ],a[jp    ],a[jt+p68],a[jp+p68],a[jt+p58],a[jp+p58],a[jt+p48],a[jp+p48]);	jt = j1+p05; jp = j2+p05;
		RADIX_15_DIT_B(s,c3m1,cn1,cn2,ss3,sn1,sn2, r75r,r75i,r76r,r76i,r77r,r77i,r78r,r78i,r79r,r79i,r80r,r80i,r81r,r81i,r82r,r82i,r83r,r83i,r84r,r84i,r85r,r85i,r86r,r86i,r87r,r87i,r88r,r88i,r89r,r89i, a[jt+p28],a[jp+p28],a[jt+p18],a[jp+p18],a[jt+p08],a[jp+p08],a[jt+p70],a[jp+p70],a[jt+p60],a[jp+p60],a[jt+p50],a[jp+p50],a[jt+p40],a[jp+p40],a[jt+p30],a[jp+p30],a[jt+p20],a[jp+p20],a[jt+p10],a[jp+p10],a[jt    ],a[jp    ],a[jt+p68],a[jp+p68],a[jt+p58],a[jp+p58],a[jt+p48],a[jp+p48],a[jt+p38],a[jp+p38]);	jt = j1+p06; jp = j2+p06;
		RADIX_15_DIT_B(s,c3m1,cn1,cn2,ss3,sn1,sn2, r90r,r90i,r91r,r91i,r92r,r92i,r93r,r93i,r94r,r94i,r95r,r95i,r96r,r96i,r97r,r97i,r98r,r98i,r99r,r99i,rA0r,rA0i,rA1r,rA1i,rA2r,rA2i,rA3r,rA3i,rA4r,rA4i, a[jt+p18],a[jp+p18],a[jt+p08],a[jp+p08],a[jt+p70],a[jp+p70],a[jt+p60],a[jp+p60],a[jt+p50],a[jp+p50],a[jt+p40],a[jp+p40],a[jt+p30],a[jp+p30],a[jt+p20],a[jp+p20],a[jt+p10],a[jp+p10],a[jt    ],a[jp    ],a[jt+p68],a[jp+p68],a[jt+p58],a[jp+p58],a[jt+p48],a[jp+p48],a[jt+p38],a[jp+p38],a[jt+p28],a[jp+p28]);	jt = j1+p07; jp = j2+p07;
		RADIX_15_DIT_B(s,c3m1,cn1,cn2,ss3,sn1,sn2, rA5r,rA5i,rA6r,rA6i,rA7r,rA7i,rA8r,rA8i,rA9r,rA9i,rB0r,rB0i,rB1r,rB1i,rB2r,rB2i,rB3r,rB3i,rB4r,rB4i,rB5r,rB5i,rB6r,rB6i,rB7r,rB7i,rB8r,rB8i,rB9r,rB9i, a[jt+p08],a[jp+p08],a[jt+p70],a[jp+p70],a[jt+p60],a[jp+p60],a[jt+p50],a[jp+p50],a[jt+p40],a[jp+p40],a[jt+p30],a[jp+p30],a[jt+p20],a[jp+p20],a[jt+p10],a[jp+p10],a[jt    ],a[jp    ],a[jt+p68],a[jp+p68],a[jt+p58],a[jp+p58],a[jt+p48],a[jp+p48],a[jt+p38],a[jp+p38],a[jt+p28],a[jp+p28],a[jt+p18],a[jp+p18]);
	}
}

#undef RADIX
