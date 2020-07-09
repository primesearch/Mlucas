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

#define RADIX 112	// Use #define rather than const int to ensure it's really a compile-time const in the C sense

/***************

int radix112_ditN_cy_dif1(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], double base[], double baseinv[], int iter, double *fracmax, uint64 p)
{
}

****************/

void radix112_dif_pass1(double a[], int n)
{
/*
!...Acronym: DIF = Decimation In Frequency
!
!...Subroutine to perform an initial radix-112 complex DIF FFT pass on the data in the length-N real vector A.
!
!   See the documentation in radix16_dif_pass for further details on storage and indexing.
!
!   See the documentation in radix7,16_dif_pass for details on the radix-7,16 subtransforms.
*/
	int j,j1,j2,jp,jt;
	// In order to preserve length-2 numeric-index-offset property here, use hex for digit:
	static int NDIVR,p01,p02,p03,p04,p05,p06,p07,p08,p09,p0a,p0b,p0c,p0d,p0e,p0f,p10,p20,p30,p40,p50,p60, first_entry=TRUE;
	const double c = 0.9238795325112867561281831, s = 0.3826834323650897717284599,
					uc1 = .62348980185873353053,	 /* cos(u) = Real part of exp(i*2*pi/7), the radix-7 fundamental sincos datum	*/
					us1 = .78183148246802980870,	 /* sin(u) = Imag part of exp(i*2*pi/7).	*/
					uc2 =-.22252093395631440426,	 /* cos(2u)	*/
					us2 = .97492791218182360702,	 /* sin(2u)	*/
					uc3 =-.90096886790241912622,	 /* cos(3u)	*/
					us3 = .43388373911755812050;	 /* sin(3u)	*/
	double re,im,rt,it, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13
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
		,rB0r,rB0i,rB1r,rB1i;

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
		p0a = p09 + p01;
		p0b = p0a + p01;
		p0c = p0b + p01;
		p0d = p0c + p01;
		p0e = p0d + p01;
		p0f = p0e + p01;
		p10 = p0f + p01;
		p20 = p10 + p10;
		p30 = p20 + p10;
		p40 = p30 + p10;
		p50 = p40 + p10;
		p60 = p50 + p10;

		p01 = p01 + ( (p01 >> DAT_BITS) << PAD_BITS );
		p02 = p02 + ( (p02 >> DAT_BITS) << PAD_BITS );
		p03 = p03 + ( (p03 >> DAT_BITS) << PAD_BITS );
		p04 = p04 + ( (p04 >> DAT_BITS) << PAD_BITS );
		p05 = p05 + ( (p05 >> DAT_BITS) << PAD_BITS );
		p06 = p06 + ( (p06 >> DAT_BITS) << PAD_BITS );
		p07 = p07 + ( (p07 >> DAT_BITS) << PAD_BITS );
		p08 = p08 + ( (p08 >> DAT_BITS) << PAD_BITS );
		p09 = p09 + ( (p09 >> DAT_BITS) << PAD_BITS );
		p0a = p0a + ( (p0a >> DAT_BITS) << PAD_BITS );
		p0b = p0b + ( (p0b >> DAT_BITS) << PAD_BITS );
		p0c = p0c + ( (p0c >> DAT_BITS) << PAD_BITS );
		p0d = p0d + ( (p0d >> DAT_BITS) << PAD_BITS );
		p0e = p0e + ( (p0e >> DAT_BITS) << PAD_BITS );
		p0f = p0f + ( (p0f >> DAT_BITS) << PAD_BITS );
		p10 = p10 + ( (p10 >> DAT_BITS) << PAD_BITS );
		p20 = p20 + ( (p20 >> DAT_BITS) << PAD_BITS );
		p30 = p30 + ( (p30 >> DAT_BITS) << PAD_BITS );
		p40 = p40 + ( (p40 >> DAT_BITS) << PAD_BITS );
		p50 = p50 + ( (p50 >> DAT_BITS) << PAD_BITS );
		p60 = p60 + ( (p60 >> DAT_BITS) << PAD_BITS );
	}

/*...The radix-112 pass is here.	*/

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
	Twiddleless version arranges 16 sets of radix-7 DFT inputs as follows:
	0 in upper left corner, decrement 16 = 0x10 horizontally and 7 vertically [rep'd in hex here to match p-indexing]:
	[Can get these by running test_fft_radix with TTYPE = 0]:

		00,60,50,40,30,20,10   00,60,50,40,30,20,10 + p00
		69,59,49,39,29,19,09   60,50,40,30,20,10,00 + p09
		62,52,42,32,22,12,02   60,50,40,30,20,10,00 + p02
		5b,4b,3b,2b,1b,0b,6b   50,40,30,20,10,00,60 + p0b
		54,44,34,24,14,04,64   50,40,30,20,10,00,60 + p04
		4d,3d,2d,1d,0d,6d,5d   40,30,20,10,00,60,50 + p0d
		46,36,26,16,06,66,56   40,30,20,10,00,60,50 + p06
		3f,2f,1f,0f,6f,5f,4f = 30,20,10,00,60,50,40 + p0f
		38,28,18,08,68,58,48   30,20,10,00,60,50,40 + p08
		31,21,11,01,61,51,41   30,20,10,00,60,50,40 + p01
		2a,1a,0a,6a,5a,4a,3a   20,10,00,60,50,40,30 + p0a
		23,13,03,63,53,43,33   20,10,00,60,50,40,30 + p03
		1c,0c,6c,5c,4c,3c,2c   10,00,60,50,40,30,20 + p0c
		15,05,65,55,45,35,25   10,00,60,50,40,30,20 + p05
		0e,6e,5e,4e,3e,2e,1e   00,60,50,40,30,20,10 + p0e
		07,67,57,47,37,27,17   00,60,50,40,30,20,10 + p07
	*/
	/*...gather the needed data (112 64-bit complex) and do 8 radix-9 transforms...*/
		RADIX_07_DFT(a[j1    ],a[j2    ],a[j1+p60],a[j2+p60],a[j1+p50],a[j2+p50],a[j1+p40],a[j2+p40],a[j1+p30],a[j2+p30],a[j1+p20],a[j2+p20],a[j1+p10],a[j2+p10], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, r00r,r00i,r01r,r01i,r02r,r02i,r03r,r03i,r04r,r04i,r05r,r05i,r06r,r06i, uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);	jt = j1+p09; jp = j2+p09;
		RADIX_07_DFT(a[jt+p60],a[jp+p60],a[jt+p50],a[jp+p50],a[jt+p40],a[jp+p40],a[jt+p30],a[jp+p30],a[jt+p20],a[jp+p20],a[jt+p10],a[jp+p10],a[jt    ],a[jp    ], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, r07r,r07i,r08r,r08i,r09r,r09i,r10r,r10i,r11r,r11i,r12r,r12i,r13r,r13i, uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);	jt = j1+p02; jp = j2+p02;
		RADIX_07_DFT(a[jt+p60],a[jp+p60],a[jt+p50],a[jp+p50],a[jt+p40],a[jp+p40],a[jt+p30],a[jp+p30],a[jt+p20],a[jp+p20],a[jt+p10],a[jp+p10],a[jt    ],a[jp    ], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, r14r,r14i,r15r,r15i,r16r,r16i,r17r,r17i,r18r,r18i,r19r,r19i,r20r,r20i, uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);	jt = j1+p0b; jp = j2+p0b;
		RADIX_07_DFT(a[jt+p50],a[jp+p50],a[jt+p40],a[jp+p40],a[jt+p30],a[jp+p30],a[jt+p20],a[jp+p20],a[jt+p10],a[jp+p10],a[jt    ],a[jp    ],a[jt+p60],a[jp+p60], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, r21r,r21i,r22r,r22i,r23r,r23i,r24r,r24i,r25r,r25i,r26r,r26i,r27r,r27i, uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);	jt = j1+p04; jp = j2+p04;
		RADIX_07_DFT(a[jt+p50],a[jp+p50],a[jt+p40],a[jp+p40],a[jt+p30],a[jp+p30],a[jt+p20],a[jp+p20],a[jt+p10],a[jp+p10],a[jt    ],a[jp    ],a[jt+p60],a[jp+p60], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, r28r,r28i,r29r,r29i,r30r,r30i,r31r,r31i,r32r,r32i,r33r,r33i,r34r,r34i, uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);	jt = j1+p0d; jp = j2+p0d;
		RADIX_07_DFT(a[jt+p40],a[jp+p40],a[jt+p30],a[jp+p30],a[jt+p20],a[jp+p20],a[jt+p10],a[jp+p10],a[jt    ],a[jp    ],a[jt+p60],a[jp+p60],a[jt+p50],a[jp+p50], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, r35r,r35i,r36r,r36i,r37r,r37i,r38r,r38i,r39r,r39i,r40r,r40i,r41r,r41i, uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);	jt = j1+p06; jp = j2+p06;
		RADIX_07_DFT(a[jt+p40],a[jp+p40],a[jt+p30],a[jp+p30],a[jt+p20],a[jp+p20],a[jt+p10],a[jp+p10],a[jt    ],a[jp    ],a[jt+p60],a[jp+p60],a[jt+p50],a[jp+p50], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, r42r,r42i,r43r,r43i,r44r,r44i,r45r,r45i,r46r,r46i,r47r,r47i,r48r,r48i, uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);	jt = j1+p0f; jp = j2+p0f;
		RADIX_07_DFT(a[jt+p30],a[jp+p30],a[jt+p20],a[jp+p20],a[jt+p10],a[jp+p10],a[jt    ],a[jp    ],a[jt+p60],a[jp+p60],a[jt+p50],a[jp+p50],a[jt+p40],a[jp+p40], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, r49r,r49i,r50r,r50i,r51r,r51i,r52r,r52i,r53r,r53i,r54r,r54i,r55r,r55i, uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);	jt = j1+p08; jp = j2+p08;
		RADIX_07_DFT(a[jt+p30],a[jp+p30],a[jt+p20],a[jp+p20],a[jt+p10],a[jp+p10],a[jt    ],a[jp    ],a[jt+p60],a[jp+p60],a[jt+p50],a[jp+p50],a[jt+p40],a[jp+p40], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, r56r,r56i,r57r,r57i,r58r,r58i,r59r,r59i,r60r,r60i,r61r,r61i,r62r,r62i, uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);	jt = j1+p01; jp = j2+p01;
		RADIX_07_DFT(a[jt+p30],a[jp+p30],a[jt+p20],a[jp+p20],a[jt+p10],a[jp+p10],a[jt    ],a[jp    ],a[jt+p60],a[jp+p60],a[jt+p50],a[jp+p50],a[jt+p40],a[jp+p40], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, r63r,r63i,r64r,r64i,r65r,r65i,r66r,r66i,r67r,r67i,r68r,r68i,r69r,r69i, uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);	jt = j1+p0a; jp = j2+p0a;
		RADIX_07_DFT(a[jt+p20],a[jp+p20],a[jt+p10],a[jp+p10],a[jt    ],a[jp    ],a[jt+p60],a[jp+p60],a[jt+p50],a[jp+p50],a[jt+p40],a[jp+p40],a[jt+p30],a[jp+p30], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, r70r,r70i,r71r,r71i,r72r,r72i,r73r,r73i,r74r,r74i,r75r,r75i,r76r,r76i, uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);	jt = j1+p03; jp = j2+p03;
		RADIX_07_DFT(a[jt+p20],a[jp+p20],a[jt+p10],a[jp+p10],a[jt    ],a[jp    ],a[jt+p60],a[jp+p60],a[jt+p50],a[jp+p50],a[jt+p40],a[jp+p40],a[jt+p30],a[jp+p30], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, r77r,r77i,r78r,r78i,r79r,r79i,r80r,r80i,r81r,r81i,r82r,r82i,r83r,r83i, uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);	jt = j1+p0c; jp = j2+p0c;
		RADIX_07_DFT(a[jt+p10],a[jp+p10],a[jt    ],a[jp    ],a[jt+p60],a[jp+p60],a[jt+p50],a[jp+p50],a[jt+p40],a[jp+p40],a[jt+p30],a[jp+p30],a[jt+p20],a[jp+p20], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, r84r,r84i,r85r,r85i,r86r,r86i,r87r,r87i,r88r,r88i,r89r,r89i,r90r,r90i, uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);	jt = j1+p05; jp = j2+p05;
		RADIX_07_DFT(a[jt+p10],a[jp+p10],a[jt    ],a[jp    ],a[jt+p60],a[jp+p60],a[jt+p50],a[jp+p50],a[jt+p40],a[jp+p40],a[jt+p30],a[jp+p30],a[jt+p20],a[jp+p20], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, r91r,r91i,r92r,r92i,r93r,r93i,r94r,r94i,r95r,r95i,r96r,r96i,r97r,r97i, uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);	jt = j1+p0e; jp = j2+p0e;
		RADIX_07_DFT(a[jt    ],a[jp    ],a[jt+p60],a[jp+p60],a[jt+p50],a[jp+p50],a[jt+p40],a[jp+p40],a[jt+p30],a[jp+p30],a[jt+p20],a[jp+p20],a[jt+p10],a[jp+p10], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, r98r,r98i,r99r,r99i,rA0r,rA0i,rA1r,rA1i,rA2r,rA2i,rA3r,rA3i,rA4r,rA4i, uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);	jt = j1+p07; jp = j2+p07;
		RADIX_07_DFT(a[jt    ],a[jp    ],a[jt+p60],a[jp+p60],a[jt+p50],a[jp+p50],a[jt+p40],a[jp+p40],a[jt+p30],a[jp+p30],a[jt+p20],a[jp+p20],a[jt+p10],a[jp+p10], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, rA5r,rA5i,rA6r,rA6i,rA7r,rA7i,rA8r,rA8i,rA9r,rA9i,rB0r,rB0i,rB1r,rB1i, uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);

	/*...and now do 7 radix-16 transforms.
	Use the supercalafragalistic Ancient Chinese Secret index-munging formula [SACSIMPF]
	to properly permute the outputs of the ensuing 15 radix-8 DFTs to the required ordering, which in terms of our 15 DFTs is

								[all indices hex here]
		00,01,02,03,04,05,06,07,09,08,0b,0a,0d,0c,0f,0e   0123456798badcfe + p00
		66,67,65,64,61,60,63,62,6f,6e,6c,6d,68,69,6a,6b   67541032fecd89ab + p60
		5b,5a,58,59,5f,5e,5c,5d,56,57,55,54,51,50,53,52   ba89fecd67541032 + p50
		42,43,41,40,46,47,45,44,4b,4a,48,49,4f,4e,4c,4d = 23106754ba89fecd + p40
		3d,3c,3f,3e,3b,3a,38,39,32,33,31,30,36,37,35,34   dcfeba8923106754 + p30
		24,25,26,27,22,23,21,20,2d,2c,2f,2e,2b,2a,28,29   45672310dcfeba89 + p20
		19,18,1b,1a,1d,1c,1f,1e,14,15,16,17,12,13,11,10   98badcfe45672310 + p10
	*/
					 /*                                    inputs                                 */ /*                       intermediates                       */ /*                 outputs                   */
		RADIX_16_DIF(r00r,r00i,r07r,r07i,r14r,r14i,r21r,r21i,r28r,r28i,r35r,r35i,r42r,r42i,r49r,r49i,r56r,r56i,r63r,r63i,r70r,r70i,r77r,r77i,r84r,r84i,r91r,r91i,r98r,r98i,rA5r,rA5i, a[j1    ],a[j2    ],a[j1+p01],a[j2+p01],a[j1+p02],a[j2+p02],a[j1+p03],a[j2+p03],a[j1+p04],a[j2+p04],a[j1+p05],a[j2+p05],a[j1+p06],a[j2+p06],a[j1+p07],a[j2+p07],a[j1+p09],a[j2+p09],a[j1+p08],a[j2+p08],a[j1+p0b],a[j2+p0b],a[j1+p0a],a[j2+p0a],a[j1+p0d],a[j2+p0d],a[j1+p0c],a[j2+p0c],a[j1+p0f],a[j2+p0f],a[j1+p0e],a[j2+p0e], c,s);	jt = j1+p60; jp = j2+p60;
		RADIX_16_DIF(r01r,r01i,r08r,r08i,r15r,r15i,r22r,r22i,r29r,r29i,r36r,r36i,r43r,r43i,r50r,r50i,r57r,r57i,r64r,r64i,r71r,r71i,r78r,r78i,r85r,r85i,r92r,r92i,r99r,r99i,rA6r,rA6i, a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p0f],a[jp+p0f],a[jt+p0e],a[jp+p0e],a[jt+p0c],a[jp+p0c],a[jt+p0d],a[jp+p0d],a[jt+p08],a[jp+p08],a[jt+p09],a[jp+p09],a[jt+p0a],a[jp+p0a],a[jt+p0b],a[jp+p0b], c,s);	jt = j1+p50; jp = j2+p50;
		RADIX_16_DIF(r02r,r02i,r09r,r09i,r16r,r16i,r23r,r23i,r30r,r30i,r37r,r37i,r44r,r44i,r51r,r51i,r58r,r58i,r65r,r65i,r72r,r72i,r79r,r79i,r86r,r86i,r93r,r93i,rA0r,rA0i,rA7r,rA7i, a[jt+p0b],a[jp+p0b],a[jt+p0a],a[jp+p0a],a[jt+p08],a[jp+p08],a[jt+p09],a[jp+p09],a[jt+p0f],a[jp+p0f],a[jt+p0e],a[jp+p0e],a[jt+p0c],a[jp+p0c],a[jt+p0d],a[jp+p0d],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02], c,s);	jt = j1+p40; jp = j2+p40;
		RADIX_16_DIF(r03r,r03i,r10r,r10i,r17r,r17i,r24r,r24i,r31r,r31i,r38r,r38i,r45r,r45i,r52r,r52i,r59r,r59i,r66r,r66i,r73r,r73i,r80r,r80i,r87r,r87i,r94r,r94i,rA1r,rA1i,rA8r,rA8i, a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04],a[jt+p0b],a[jp+p0b],a[jt+p0a],a[jp+p0a],a[jt+p08],a[jp+p08],a[jt+p09],a[jp+p09],a[jt+p0f],a[jp+p0f],a[jt+p0e],a[jp+p0e],a[jt+p0c],a[jp+p0c],a[jt+p0d],a[jp+p0d], c,s);	jt = j1+p30; jp = j2+p30;
		RADIX_16_DIF(r04r,r04i,r11r,r11i,r18r,r18i,r25r,r25i,r32r,r32i,r39r,r39i,r46r,r46i,r53r,r53i,r60r,r60i,r67r,r67i,r74r,r74i,r81r,r81i,r88r,r88i,r95r,r95i,rA2r,rA2i,rA9r,rA9i, a[jt+p0d],a[jp+p0d],a[jt+p0c],a[jp+p0c],a[jt+p0f],a[jp+p0f],a[jt+p0e],a[jp+p0e],a[jt+p0b],a[jp+p0b],a[jt+p0a],a[jp+p0a],a[jt+p08],a[jp+p08],a[jt+p09],a[jp+p09],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04], c,s);	jt = j1+p20; jp = j2+p20;
		RADIX_16_DIF(r05r,r05i,r12r,r12i,r19r,r19i,r26r,r26i,r33r,r33i,r40r,r40i,r47r,r47i,r54r,r54i,r61r,r61i,r68r,r68i,r75r,r75i,r82r,r82i,r89r,r89i,r96r,r96i,rA3r,rA3i,rB0r,rB0i, a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p0d],a[jp+p0d],a[jt+p0c],a[jp+p0c],a[jt+p0f],a[jp+p0f],a[jt+p0e],a[jp+p0e],a[jt+p0b],a[jp+p0b],a[jt+p0a],a[jp+p0a],a[jt+p08],a[jp+p08],a[jt+p09],a[jp+p09], c,s);	jt = j1+p10; jp = j2+p10;
		RADIX_16_DIF(r06r,r06i,r13r,r13i,r20r,r20i,r27r,r27i,r34r,r34i,r41r,r41i,r48r,r48i,r55r,r55i,r62r,r62i,r69r,r69i,r76r,r76i,r83r,r83i,r90r,r90i,r97r,r97i,rA4r,rA4i,rB1r,rB1i, a[jt+p09],a[jp+p09],a[jt+p08],a[jp+p08],a[jt+p0b],a[jp+p0b],a[jt+p0a],a[jp+p0a],a[jt+p0d],a[jp+p0d],a[jt+p0c],a[jp+p0c],a[jt+p0f],a[jp+p0f],a[jt+p0e],a[jp+p0e],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ], c,s);
	}
}

/***************/

void radix112_dit_pass1(double a[], int n)
{
/*
!...Acronym: DIF = Decimation In Frequency
!
!...Subroutine to perform an initial radix-112 complex DIT FFT pass on the data in the length-N real vector A.
!
!   See the documentation in radix16_dif_pass for further details on storage and indexing.
!
!   See the documentation in radix7_dif_pass for details on the radix-7 subtransforms.
*/
	int j,j1,j2,jp,jt;
	// In order to preserve length-2 numeric-index-offset property here, use hex for digit:
	static int NDIVR,p01,p02,p03,p04,p05,p06,p07,p08,p09,p0a,p0b,p0c,p0d,p0e,p0f,p10,p20,p30,p40,p50,p60, first_entry=TRUE;
	const double c = 0.9238795325112867561281831, s = 0.3826834323650897717284599,
					uc1 = .62348980185873353053,	 /* cos(u) = Real part of exp(i*2*pi/7), the radix-7 fundamental sincos datum	*/
					us1 = .78183148246802980870,	 /* sin(u) = Imag part of exp(i*2*pi/7).	*/
					uc2 =-.22252093395631440426,	 /* cos(2u)	*/
					us2 = .97492791218182360702,	 /* sin(2u)	*/
					uc3 =-.90096886790241912622,	 /* cos(3u)	*/
					us3 = .43388373911755812050;	 /* sin(3u)	*/
	double re,im,rt,it, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13
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
		,rB0r,rB0i,rB1r,rB1i;

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
		p0a = p09 + p01;
		p0b = p0a + p01;
		p0c = p0b + p01;
		p0d = p0c + p01;
		p0e = p0d + p01;
		p0f = p0e + p01;
		p10 = p0f + p01;
		p20 = p10 + p10;
		p30 = p20 + p10;
		p40 = p30 + p10;
		p50 = p40 + p10;
		p60 = p50 + p10;

		p01 = p01 + ( (p01 >> DAT_BITS) << PAD_BITS );
		p02 = p02 + ( (p02 >> DAT_BITS) << PAD_BITS );
		p03 = p03 + ( (p03 >> DAT_BITS) << PAD_BITS );
		p04 = p04 + ( (p04 >> DAT_BITS) << PAD_BITS );
		p05 = p05 + ( (p05 >> DAT_BITS) << PAD_BITS );
		p06 = p06 + ( (p06 >> DAT_BITS) << PAD_BITS );
		p07 = p07 + ( (p07 >> DAT_BITS) << PAD_BITS );
		p08 = p08 + ( (p08 >> DAT_BITS) << PAD_BITS );
		p09 = p09 + ( (p09 >> DAT_BITS) << PAD_BITS );
		p0a = p0a + ( (p0a >> DAT_BITS) << PAD_BITS );
		p0b = p0b + ( (p0b >> DAT_BITS) << PAD_BITS );
		p0c = p0c + ( (p0c >> DAT_BITS) << PAD_BITS );
		p0d = p0d + ( (p0d >> DAT_BITS) << PAD_BITS );
		p0e = p0e + ( (p0e >> DAT_BITS) << PAD_BITS );
		p0f = p0f + ( (p0f >> DAT_BITS) << PAD_BITS );
		p10 = p10 + ( (p10 >> DAT_BITS) << PAD_BITS );
		p20 = p20 + ( (p20 >> DAT_BITS) << PAD_BITS );
		p30 = p30 + ( (p30 >> DAT_BITS) << PAD_BITS );
		p40 = p40 + ( (p40 >> DAT_BITS) << PAD_BITS );
		p50 = p50 + ( (p50 >> DAT_BITS) << PAD_BITS );
		p60 = p60 + ( (p60 >> DAT_BITS) << PAD_BITS );
	}

/*...The radix-112 pass is here.	*/

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
	Twiddleless version uses same linear-index-vector-form permutation as in DIF -
	Remember, inputs to DIT are bit-reversed, so using output of test_fft_radix(),
	here are the 7 index-offset 16-tets going into the radix-16 DFTs:

	Combined DIT input-scramble array =

										[all indices hex here]
	00,01,03,02,07,06,05,04,0f,0e,0d,0c,0b,0a,09,08   01327654fedcba98 + p00
	5b,5a,59,58,5d,5c,5e,5f,53,52,51,50,55,54,56,57   ba98dcef32105467 + p50
	3d,3c,3e,3f,39,38,3a,3b,35,34,36,37,31,30,32,33   dcef98ab54671023 + p30
	19,18,1a,1b,1e,1f,1c,1d,11,10,12,13,16,17,14,15   98abefcd10236745 + p10
	66,67,64,65,62,63,60,61,6a,6b,68,69,6c,6d,6f,6e   67452301ab89cdfe + p60
	42,43,40,41,44,45,47,46,4c,4d,4f,4e,48,49,4b,4a = 23014576cdfe89ba + p40
	24,25,27,26,20,21,23,22,28,29,2b,2a,2f,2e,2d,2c   4576013289bafedc + p20
	*/
	/*...gather the needed data (112 64-bit complex) and do 7 radix-16 transforms...*/
		RADIX_16_DIT(a[j1    ],a[j2    ],a[j1+p01],a[j2+p01],a[j1+p03],a[j2+p03],a[j1+p02],a[j2+p02],a[j1+p07],a[j2+p07],a[j1+p06],a[j2+p06],a[j1+p05],a[j2+p05],a[j1+p04],a[j2+p04],a[j1+p0f],a[j2+p0f],a[j1+p0e],a[j2+p0e],a[j1+p0d],a[j2+p0d],a[j1+p0c],a[j2+p0c],a[j1+p0b],a[j2+p0b],a[j1+p0a],a[j2+p0a],a[j1+p09],a[j2+p09],a[j1+p08],a[j2+p08], r00r,r00i,r07r,r07i,r14r,r14i,r21r,r21i,r28r,r28i,r35r,r35i,r42r,r42i,r49r,r49i,r56r,r56i,r63r,r63i,r70r,r70i,r77r,r77i,r84r,r84i,r91r,r91i,r98r,r98i,rA5r,rA5i, c,s);	jt = j1+p50; jp = j2+p50;
		RADIX_16_DIT(a[jt+p0b],a[jp+p0b],a[jt+p0a],a[jp+p0a],a[jt+p09],a[jp+p09],a[jt+p08],a[jp+p08],a[jt+p0d],a[jp+p0d],a[jt+p0c],a[jp+p0c],a[jt+p0e],a[jp+p0e],a[jt+p0f],a[jp+p0f],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07], r01r,r01i,r08r,r08i,r15r,r15i,r22r,r22i,r29r,r29i,r36r,r36i,r43r,r43i,r50r,r50i,r57r,r57i,r64r,r64i,r71r,r71i,r78r,r78i,r85r,r85i,r92r,r92i,r99r,r99i,rA6r,rA6i, c,s);	jt = j1+p30; jp = j2+p30;
		RADIX_16_DIT(a[jt+p0d],a[jp+p0d],a[jt+p0c],a[jp+p0c],a[jt+p0e],a[jp+p0e],a[jt+p0f],a[jp+p0f],a[jt+p09],a[jp+p09],a[jt+p08],a[jp+p08],a[jt+p0a],a[jp+p0a],a[jt+p0b],a[jp+p0b],a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03], r02r,r02i,r09r,r09i,r16r,r16i,r23r,r23i,r30r,r30i,r37r,r37i,r44r,r44i,r51r,r51i,r58r,r58i,r65r,r65i,r72r,r72i,r79r,r79i,r86r,r86i,r93r,r93i,rA0r,rA0i,rA7r,rA7i, c,s);	jt = j1+p10; jp = j2+p10;
		RADIX_16_DIT(a[jt+p09],a[jp+p09],a[jt+p08],a[jp+p08],a[jt+p0a],a[jp+p0a],a[jt+p0b],a[jp+p0b],a[jt+p0e],a[jp+p0e],a[jt+p0f],a[jp+p0f],a[jt+p0c],a[jp+p0c],a[jt+p0d],a[jp+p0d],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05], r03r,r03i,r10r,r10i,r17r,r17i,r24r,r24i,r31r,r31i,r38r,r38i,r45r,r45i,r52r,r52i,r59r,r59i,r66r,r66i,r73r,r73i,r80r,r80i,r87r,r87i,r94r,r94i,rA1r,rA1i,rA8r,rA8i, c,s);	jt = j1+p60; jp = j2+p60;
		RADIX_16_DIT(a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p0a],a[jp+p0a],a[jt+p0b],a[jp+p0b],a[jt+p08],a[jp+p08],a[jt+p09],a[jp+p09],a[jt+p0c],a[jp+p0c],a[jt+p0d],a[jp+p0d],a[jt+p0f],a[jp+p0f],a[jt+p0e],a[jp+p0e], r04r,r04i,r11r,r11i,r18r,r18i,r25r,r25i,r32r,r32i,r39r,r39i,r46r,r46i,r53r,r53i,r60r,r60i,r67r,r67i,r74r,r74i,r81r,r81i,r88r,r88i,r95r,r95i,rA2r,rA2i,rA9r,rA9i, c,s);	jt = j1+p40; jp = j2+p40;
		RADIX_16_DIT(a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],a[jt+p0c],a[jp+p0c],a[jt+p0d],a[jp+p0d],a[jt+p0f],a[jp+p0f],a[jt+p0e],a[jp+p0e],a[jt+p08],a[jp+p08],a[jt+p09],a[jp+p09],a[jt+p0b],a[jp+p0b],a[jt+p0a],a[jp+p0a], r05r,r05i,r12r,r12i,r19r,r19i,r26r,r26i,r33r,r33i,r40r,r40i,r47r,r47i,r54r,r54i,r61r,r61i,r68r,r68i,r75r,r75i,r82r,r82i,r89r,r89i,r96r,r96i,rA3r,rA3i,rB0r,rB0i, c,s);	jt = j1+p20; jp = j2+p20;
		RADIX_16_DIT(a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p08],a[jp+p08],a[jt+p09],a[jp+p09],a[jt+p0b],a[jp+p0b],a[jt+p0a],a[jp+p0a],a[jt+p0f],a[jp+p0f],a[jt+p0e],a[jp+p0e],a[jt+p0d],a[jp+p0d],a[jt+p0c],a[jp+p0c], r06r,r06i,r13r,r13i,r20r,r20i,r27r,r27i,r34r,r34i,r41r,r41i,r48r,r48i,r55r,r55i,r62r,r62i,r69r,r69i,r76r,r76i,r83r,r83i,r90r,r90i,r97r,r97i,rA4r,rA4i,rB1r,rB1i, c,s);

	/*...and now do 16 radix-7 transforms: 

	Use the inverse supercalafragalistic Ancient Chinese Secret index-munging formula [iSACSIMPF]
	to properly permute the outputs of the 8 radix-15 DFTs:
								[all indices hex here]
		00,40,10,50,20,60,30   38,28,18,08,70,60,50,40,30,20,10,00,68,58,48 + p00
		3f,0f,4f,1f,5f,2f,6f   28,18,08,70,60,50,40,30,20,10,00,68,58,48,38 + p0f
		0e,4e,1e,5e,2e,6e,3e   18,08,70,60,50,40,30,20,10,00,68,58,48,38,28 + p0e
		4d,1d,5d,2d,6d,3d,0d   08,70,60,50,40,30,20,10,00,68,58,48,38,28,18 + p0d
		1c,5c,2c,6c,3c,0c,4c   00,68,58,48,38,28,18,08,70,60,50,40,30,20,10 + p0c
		5b,2b,6b,3b,0b,4b,1b   68,58,48,38,28,18,08,70,60,50,40,30,20,10,00 + p0b
		2a,6a,3a,0a,4a,1a,5a   58,48,38,28,18,08,70,60,50,40,30,20,10,00,68 + p0a
		69,39,09,49,19,59,29 = 48,38,28,18,08,70,60,50,40,30,20,10,00,68,58 + p09
		38,08,48,18,58,28,68   38,28,18,08,70,60,50,40,30,20,10,00,68,58,48 + p08
		07,47,17,57,27,67,37   28,18,08,70,60,50,40,30,20,10,00,68,58,48,38 + p07
		46,16,56,26,66,36,06   18,08,70,60,50,40,30,20,10,00,68,58,48,38,28 + p06
		15,55,25,65,35,05,45   08,70,60,50,40,30,20,10,00,68,58,48,38,28,18 + p05
		54,24,64,34,04,44,14   38,28,18,08,70,60,50,40,30,20,10,00,68,58,48 + p04
		23,63,33,03,43,13,53   28,18,08,70,60,50,40,30,20,10,00,68,58,48,38 + p03
		62,32,02,42,12,52,22   18,08,70,60,50,40,30,20,10,00,68,58,48,38,28 + p02
		31,01,41,11,51,21,61   08,70,60,50,40,30,20,10,00,68,58,48,38,28,18 + p01
	*/
		RADIX_07_DFT(r00r,r00i,r01r,r01i,r02r,r02i,r03r,r03i,r04r,r04i,r05r,r05i,r06r,r06i, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[j1    ],a[j2    ],a[j1+p40],a[j2+p40],a[j1+p10],a[j2+p10],a[j1+p50],a[j2+p50],a[j1+p20],a[j2+p20],a[j1+p60],a[j2+p60],a[j1+p30],a[j2+p30], uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);	jt = j1+p0f; jp = j2+p0f;
		RADIX_07_DFT(r07r,r07i,r08r,r08i,r09r,r09i,r10r,r10i,r11r,r11i,r12r,r12i,r13r,r13i, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[jt+p30],a[jp+p30],a[jt    ],a[jp    ],a[jt+p40],a[jp+p40],a[jt+p10],a[jp+p10],a[jt+p50],a[jp+p50],a[jt+p20],a[jp+p20],a[jt+p60],a[jp+p60], uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);	jt = j1+p0e; jp = j2+p0e;
		RADIX_07_DFT(r14r,r14i,r15r,r15i,r16r,r16i,r17r,r17i,r18r,r18i,r19r,r19i,r20r,r20i, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[jt    ],a[jp    ],a[jt+p40],a[jp+p40],a[jt+p10],a[jp+p10],a[jt+p50],a[jp+p50],a[jt+p20],a[jp+p20],a[jt+p60],a[jp+p60],a[jt+p30],a[jp+p30], uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);	jt = j1+p0d; jp = j2+p0d;
		RADIX_07_DFT(r21r,r21i,r22r,r22i,r23r,r23i,r24r,r24i,r25r,r25i,r26r,r26i,r27r,r27i, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[jt+p40],a[jp+p40],a[jt+p10],a[jp+p10],a[jt+p50],a[jp+p50],a[jt+p20],a[jp+p20],a[jt+p60],a[jp+p60],a[jt+p30],a[jp+p30],a[jt    ],a[jp    ], uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);	jt = j1+p0c; jp = j2+p0c;
		RADIX_07_DFT(r28r,r28i,r29r,r29i,r30r,r30i,r31r,r31i,r32r,r32i,r33r,r33i,r34r,r34i, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[jt+p10],a[jp+p10],a[jt+p50],a[jp+p50],a[jt+p20],a[jp+p20],a[jt+p60],a[jp+p60],a[jt+p30],a[jp+p30],a[jt    ],a[jp    ],a[jt+p40],a[jp+p40], uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);	jt = j1+p0b; jp = j2+p0b;
		RADIX_07_DFT(r35r,r35i,r36r,r36i,r37r,r37i,r38r,r38i,r39r,r39i,r40r,r40i,r41r,r41i, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[jt+p50],a[jp+p50],a[jt+p20],a[jp+p20],a[jt+p60],a[jp+p60],a[jt+p30],a[jp+p30],a[jt    ],a[jp    ],a[jt+p40],a[jp+p40],a[jt+p10],a[jp+p10], uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);	jt = j1+p0a; jp = j2+p0a;
		RADIX_07_DFT(r42r,r42i,r43r,r43i,r44r,r44i,r45r,r45i,r46r,r46i,r47r,r47i,r48r,r48i, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[jt+p20],a[jp+p20],a[jt+p60],a[jp+p60],a[jt+p30],a[jp+p30],a[jt    ],a[jp    ],a[jt+p40],a[jp+p40],a[jt+p10],a[jp+p10],a[jt+p50],a[jp+p50], uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);	jt = j1+p09; jp = j2+p09;
		RADIX_07_DFT(r49r,r49i,r50r,r50i,r51r,r51i,r52r,r52i,r53r,r53i,r54r,r54i,r55r,r55i, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[jt+p60],a[jp+p60],a[jt+p30],a[jp+p30],a[jt    ],a[jp    ],a[jt+p40],a[jp+p40],a[jt+p10],a[jp+p10],a[jt+p50],a[jp+p50],a[jt+p20],a[jp+p20], uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);	jt = j1+p08; jp = j2+p08;
		RADIX_07_DFT(r56r,r56i,r57r,r57i,r58r,r58i,r59r,r59i,r60r,r60i,r61r,r61i,r62r,r62i, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[jt+p30],a[jp+p30],a[jt    ],a[jp    ],a[jt+p40],a[jp+p40],a[jt+p10],a[jp+p10],a[jt+p50],a[jp+p50],a[jt+p20],a[jp+p20],a[jt+p60],a[jp+p60], uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);	jt = j1+p07; jp = j2+p07;
		RADIX_07_DFT(r63r,r63i,r64r,r64i,r65r,r65i,r66r,r66i,r67r,r67i,r68r,r68i,r69r,r69i, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[jt    ],a[jp    ],a[jt+p40],a[jp+p40],a[jt+p10],a[jp+p10],a[jt+p50],a[jp+p50],a[jt+p20],a[jp+p20],a[jt+p60],a[jp+p60],a[jt+p30],a[jp+p30], uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);	jt = j1+p06; jp = j2+p06;
		RADIX_07_DFT(r70r,r70i,r71r,r71i,r72r,r72i,r73r,r73i,r74r,r74i,r75r,r75i,r76r,r76i, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[jt+p40],a[jp+p40],a[jt+p10],a[jp+p10],a[jt+p50],a[jp+p50],a[jt+p20],a[jp+p20],a[jt+p60],a[jp+p60],a[jt+p30],a[jp+p30],a[jt    ],a[jp    ], uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);	jt = j1+p05; jp = j2+p05;
		RADIX_07_DFT(r77r,r77i,r78r,r78i,r79r,r79i,r80r,r80i,r81r,r81i,r82r,r82i,r83r,r83i, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[jt+p10],a[jp+p10],a[jt+p50],a[jp+p50],a[jt+p20],a[jp+p20],a[jt+p60],a[jp+p60],a[jt+p30],a[jp+p30],a[jt    ],a[jp    ],a[jt+p40],a[jp+p40], uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);	jt = j1+p04; jp = j2+p04;
		RADIX_07_DFT(r84r,r84i,r85r,r85i,r86r,r86i,r87r,r87i,r88r,r88i,r89r,r89i,r90r,r90i, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[jt+p50],a[jp+p50],a[jt+p20],a[jp+p20],a[jt+p60],a[jp+p60],a[jt+p30],a[jp+p30],a[jt    ],a[jp    ],a[jt+p40],a[jp+p40],a[jt+p10],a[jp+p10], uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);	jt = j1+p03; jp = j2+p03;
		RADIX_07_DFT(r91r,r91i,r92r,r92i,r93r,r93i,r94r,r94i,r95r,r95i,r96r,r96i,r97r,r97i, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[jt+p20],a[jp+p20],a[jt+p60],a[jp+p60],a[jt+p30],a[jp+p30],a[jt    ],a[jp    ],a[jt+p40],a[jp+p40],a[jt+p10],a[jp+p10],a[jt+p50],a[jp+p50], uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);	jt = j1+p02; jp = j2+p02;
		RADIX_07_DFT(r98r,r98i,r99r,r99i,rA0r,rA0i,rA1r,rA1i,rA2r,rA2i,rA3r,rA3i,rA4r,rA4i, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[jt+p60],a[jp+p60],a[jt+p30],a[jp+p30],a[jt    ],a[jp    ],a[jt+p40],a[jp+p40],a[jt+p10],a[jp+p10],a[jt+p50],a[jp+p50],a[jt+p20],a[jp+p20], uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);	jt = j1+p01; jp = j2+p01;
		RADIX_07_DFT(rA5r,rA5i,rA6r,rA6i,rA7r,rA7i,rA8r,rA8i,rA9r,rA9i,rB0r,rB0i,rB1r,rB1i, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[jt+p30],a[jp+p30],a[jt    ],a[jp    ],a[jt+p40],a[jp+p40],a[jt+p10],a[jp+p10],a[jt+p50],a[jp+p50],a[jt+p20],a[jp+p20],a[jt+p60],a[jp+p60], uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);
	}
}

#undef RADIX
