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

#define RADIX 96	// Use #define rather than const int to ensure it's really a compile-time const in the C sense

/***************

int radix96_ditN_cy_dif1(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], double base[], double baseinv[], int iter, double *fracmax, uint64 p)
{
}

****************/

void radix96_dif_pass1(double a[], int n)
{
/*
!...Acronym: DIF = Decimation In Frequency
!
!...Subroutine to perform an initial radix-96 complex DIF FFT pass on the data in the length-N real vector A.
!
!   See the documentation in radix16_dif_pass for further details on storage and indexing.
!
!   See the documentation in radix3_dif_pass for details on the radix-3 subtransforms.
*/
	int j,j1,j2,jp,jt;
	// In order to preserve length-2 numeric-index-offset property here, use hex for digit:
	static int NDIVR,p01,p02,p03,p04,p05,p06,p07,p08,p09,p0a,p0b,p0c,p0d,p0e,p0f,p10,p20,p30,p40,p50, first_entry=TRUE;
	const double c3m1= -1.50000000000000000000, s     = 0.86602540378443864675;	/* cos(twopi/3)-1, sin(twopi/3)		*/
	static int arr_offsets[32];
	// Local storage: We must use an array here because scalars have no guarantees about relative address offsets
	// [and even if those are contiguous-as-hoped-for, they may run in reverse]; Make array type (struct complex)
	// to allow us to use the same offset-indexing as in the original radix-32 in-place DFT macros:
	struct complex r[RADIX];
	double t00,t01,t02,t03,t04,t05;

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
	}

/*...The radix-96 pass is here.	*/

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
	Twiddleless version arranges 32 sets of radix-3 DFT inputs as follows: 0 in upper left corner,
	decrement 32 horizontally and 3 vertically, all (mod 96):

		00,40,20   00,40,20 + p00		30,10,50   30,10,50 + p00
		5d,3d,1d   50,30,10 + p0d		2d,0d,4d   20,00,40 + p0d
		5a,3a,1a   50,30,10 + p0a		2a,0a,4a   20,00,40 + p0a
		57,37,17   50,30,10 + p07		27,07,47   20,00,40 + p07
		54,34,14   50,30,10 + p04		24,04,44   20,00,40 + p04
		51,31,11   50,30,10 + p01		21,01,41   20,00,40 + p01
		4e,2e,0e   40,20,00 + p0e		1e,5e,3e   10,50,30 + p0e
		4b,2b,0b   40,20,00 + p0b		1b,5b,3b   10,50,30 + p0b
		48,28,08   40,20,00 + p08		18,58,38   10,50,30 + p08
		45,25,05   40,20,00 + p05		15,55,35   10,50,30 + p05
		42,22,02   40,20,00 + p02		12,52,32   10,50,30 + p02
		3f,1f,5f   30,10,50 + p0f		0f,4f,2f   00,40,20 + p0f
		3c,1c,5c   30,10,50 + p0c		0c,4c,2c   00,40,20 + p0c
		39,19,59   30,10,50 + p09		09,49,29   00,40,20 + p09
		36,16,56   30,10,50 + p06		06,46,26   00,40,20 + p06
		33,13,53 = 30,10,50 + p03		03,43,23   00,40,20 + p03
		[cont. in rcol]
	*/
	/*...gather the needed data (96 64-bit complex) and do 16 radix-3 transforms - We want unit-strides in the radix32-DFT macros, so use large output strides here: */
					 /*                        inputs                         */ /*             intermediates                 */ /*                 outputs            */
		RADIX_03_DFT(s,c3m1, a[j1    ],a[j2    ],a[j1+p40],a[j2+p40],a[j1+p20],a[j2+p20], t00,t01,t02,t03,t04,t05, r[ 0].re,r[ 0].im,r[32].re,r[32].im,r[64].re,r[64].im);	jt = j1 + p0d; jp = j2 + p0d;
		RADIX_03_DFT(s,c3m1, a[jt+p50],a[jp+p50],a[jt+p30],a[jp+p30],a[jt+p10],a[jp+p10], t00,t01,t02,t03,t04,t05, r[ 1].re,r[ 1].im,r[33].re,r[33].im,r[65].re,r[65].im);	jt = j1 + p0a; jp = j2 + p0a;
		RADIX_03_DFT(s,c3m1, a[jt+p50],a[jp+p50],a[jt+p30],a[jp+p30],a[jt+p10],a[jp+p10], t00,t01,t02,t03,t04,t05, r[ 2].re,r[ 2].im,r[34].re,r[34].im,r[66].re,r[66].im);	jt = j1 + p07; jp = j2 + p07;
		RADIX_03_DFT(s,c3m1, a[jt+p50],a[jp+p50],a[jt+p30],a[jp+p30],a[jt+p10],a[jp+p10], t00,t01,t02,t03,t04,t05, r[ 3].re,r[ 3].im,r[35].re,r[35].im,r[67].re,r[67].im);	jt = j1 + p04; jp = j2 + p04;
		RADIX_03_DFT(s,c3m1, a[jt+p50],a[jp+p50],a[jt+p30],a[jp+p30],a[jt+p10],a[jp+p10], t00,t01,t02,t03,t04,t05, r[ 4].re,r[ 4].im,r[36].re,r[36].im,r[68].re,r[68].im);	jt = j1 + p01; jp = j2 + p01;
		RADIX_03_DFT(s,c3m1, a[jt+p50],a[jp+p50],a[jt+p30],a[jp+p30],a[jt+p10],a[jp+p10], t00,t01,t02,t03,t04,t05, r[ 5].re,r[ 5].im,r[37].re,r[37].im,r[69].re,r[69].im);	jt = j1 + p0e; jp = j2 + p0e;
		RADIX_03_DFT(s,c3m1, a[jt+p40],a[jp+p40],a[jt+p20],a[jp+p20],a[jt    ],a[jp    ], t00,t01,t02,t03,t04,t05, r[ 6].re,r[ 6].im,r[38].re,r[38].im,r[70].re,r[70].im);	jt = j1 + p0b; jp = j2 + p0b;
		RADIX_03_DFT(s,c3m1, a[jt+p40],a[jp+p40],a[jt+p20],a[jp+p20],a[jt    ],a[jp    ], t00,t01,t02,t03,t04,t05, r[ 7].re,r[ 7].im,r[39].re,r[39].im,r[71].re,r[71].im);	jt = j1 + p08; jp = j2 + p08;
		RADIX_03_DFT(s,c3m1, a[jt+p40],a[jp+p40],a[jt+p20],a[jp+p20],a[jt    ],a[jp    ], t00,t01,t02,t03,t04,t05, r[ 8].re,r[ 8].im,r[40].re,r[40].im,r[72].re,r[72].im);	jt = j1 + p05; jp = j2 + p05;
		RADIX_03_DFT(s,c3m1, a[jt+p40],a[jp+p40],a[jt+p20],a[jp+p20],a[jt    ],a[jp    ], t00,t01,t02,t03,t04,t05, r[ 9].re,r[ 9].im,r[41].re,r[41].im,r[73].re,r[73].im);	jt = j1 + p02; jp = j2 + p02;
		RADIX_03_DFT(s,c3m1, a[jt+p40],a[jp+p40],a[jt+p20],a[jp+p20],a[jt    ],a[jp    ], t00,t01,t02,t03,t04,t05, r[10].re,r[10].im,r[42].re,r[42].im,r[74].re,r[74].im);	jt = j1 + p0f; jp = j2 + p0f;
		RADIX_03_DFT(s,c3m1, a[jt+p30],a[jp+p30],a[jt+p10],a[jp+p10],a[jt+p50],a[jp+p50], t00,t01,t02,t03,t04,t05, r[11].re,r[11].im,r[43].re,r[43].im,r[75].re,r[75].im);	jt = j1 + p0c; jp = j2 + p0c;
		RADIX_03_DFT(s,c3m1, a[jt+p30],a[jp+p30],a[jt+p10],a[jp+p10],a[jt+p50],a[jp+p50], t00,t01,t02,t03,t04,t05, r[12].re,r[12].im,r[44].re,r[44].im,r[76].re,r[76].im);	jt = j1 + p09; jp = j2 + p09;
		RADIX_03_DFT(s,c3m1, a[jt+p30],a[jp+p30],a[jt+p10],a[jp+p10],a[jt+p50],a[jp+p50], t00,t01,t02,t03,t04,t05, r[13].re,r[13].im,r[45].re,r[45].im,r[77].re,r[77].im);	jt = j1 + p06; jp = j2 + p06;
		RADIX_03_DFT(s,c3m1, a[jt+p30],a[jp+p30],a[jt+p10],a[jp+p10],a[jt+p50],a[jp+p50], t00,t01,t02,t03,t04,t05, r[14].re,r[14].im,r[46].re,r[46].im,r[78].re,r[78].im);	jt = j1 + p03; jp = j2 + p03;
		RADIX_03_DFT(s,c3m1, a[jt+p30],a[jp+p30],a[jt+p10],a[jp+p10],a[jt+p50],a[jp+p50], t00,t01,t02,t03,t04,t05, r[15].re,r[15].im,r[47].re,r[47].im,r[79].re,r[79].im);
		RADIX_03_DFT(s,c3m1, a[j1+p30],a[j2+p30],a[j1+p10],a[j2+p10],a[j1+p50],a[j2+p50], t00,t01,t02,t03,t04,t05, r[16].re,r[16].im,r[48].re,r[48].im,r[80].re,r[80].im);	jt = j1 + p0d; jp = j2 + p0d;
		RADIX_03_DFT(s,c3m1, a[jt+p20],a[jp+p20],a[jt    ],a[jp    ],a[jt+p40],a[jp+p40], t00,t01,t02,t03,t04,t05, r[17].re,r[17].im,r[49].re,r[49].im,r[81].re,r[81].im);	jt = j1 + p0a; jp = j2 + p0a;
		RADIX_03_DFT(s,c3m1, a[jt+p20],a[jp+p20],a[jt    ],a[jp    ],a[jt+p40],a[jp+p40], t00,t01,t02,t03,t04,t05, r[18].re,r[18].im,r[50].re,r[50].im,r[82].re,r[82].im);	jt = j1 + p07; jp = j2 + p07;
		RADIX_03_DFT(s,c3m1, a[jt+p20],a[jp+p20],a[jt    ],a[jp    ],a[jt+p40],a[jp+p40], t00,t01,t02,t03,t04,t05, r[19].re,r[19].im,r[51].re,r[51].im,r[83].re,r[83].im);	jt = j1 + p04; jp = j2 + p04;
		RADIX_03_DFT(s,c3m1, a[jt+p20],a[jp+p20],a[jt    ],a[jp    ],a[jt+p40],a[jp+p40], t00,t01,t02,t03,t04,t05, r[20].re,r[20].im,r[52].re,r[52].im,r[84].re,r[84].im);	jt = j1 + p01; jp = j2 + p01;
		RADIX_03_DFT(s,c3m1, a[jt+p20],a[jp+p20],a[jt    ],a[jp    ],a[jt+p40],a[jp+p40], t00,t01,t02,t03,t04,t05, r[21].re,r[21].im,r[53].re,r[53].im,r[85].re,r[85].im);	jt = j1 + p0e; jp = j2 + p0e;
		RADIX_03_DFT(s,c3m1, a[jt+p10],a[jp+p10],a[jt+p50],a[jp+p50],a[jt+p30],a[jp+p30], t00,t01,t02,t03,t04,t05, r[22].re,r[22].im,r[54].re,r[54].im,r[86].re,r[86].im);	jt = j1 + p0b; jp = j2 + p0b;
		RADIX_03_DFT(s,c3m1, a[jt+p10],a[jp+p10],a[jt+p50],a[jp+p50],a[jt+p30],a[jp+p30], t00,t01,t02,t03,t04,t05, r[23].re,r[23].im,r[55].re,r[55].im,r[87].re,r[87].im);	jt = j1 + p08; jp = j2 + p08;
		RADIX_03_DFT(s,c3m1, a[jt+p10],a[jp+p10],a[jt+p50],a[jp+p50],a[jt+p30],a[jp+p30], t00,t01,t02,t03,t04,t05, r[24].re,r[24].im,r[56].re,r[56].im,r[88].re,r[88].im);	jt = j1 + p05; jp = j2 + p05;
		RADIX_03_DFT(s,c3m1, a[jt+p10],a[jp+p10],a[jt+p50],a[jp+p50],a[jt+p30],a[jp+p30], t00,t01,t02,t03,t04,t05, r[25].re,r[25].im,r[57].re,r[57].im,r[89].re,r[89].im);	jt = j1 + p02; jp = j2 + p02;
		RADIX_03_DFT(s,c3m1, a[jt+p10],a[jp+p10],a[jt+p50],a[jp+p50],a[jt+p30],a[jp+p30], t00,t01,t02,t03,t04,t05, r[26].re,r[26].im,r[58].re,r[58].im,r[90].re,r[90].im);	jt = j1 + p0f; jp = j2 + p0f;
		RADIX_03_DFT(s,c3m1, a[jt    ],a[jp    ],a[jt+p40],a[jp+p40],a[jt+p20],a[jp+p20], t00,t01,t02,t03,t04,t05, r[27].re,r[27].im,r[59].re,r[59].im,r[91].re,r[91].im);	jt = j1 + p0c; jp = j2 + p0c;
		RADIX_03_DFT(s,c3m1, a[jt    ],a[jp    ],a[jt+p40],a[jp+p40],a[jt+p20],a[jp+p20], t00,t01,t02,t03,t04,t05, r[28].re,r[28].im,r[60].re,r[60].im,r[92].re,r[92].im);	jt = j1 + p09; jp = j2 + p09;
		RADIX_03_DFT(s,c3m1, a[jt    ],a[jp    ],a[jt+p40],a[jp+p40],a[jt+p20],a[jp+p20], t00,t01,t02,t03,t04,t05, r[29].re,r[29].im,r[61].re,r[61].im,r[93].re,r[93].im);	jt = j1 + p06; jp = j2 + p06;
		RADIX_03_DFT(s,c3m1, a[jt    ],a[jp    ],a[jt+p40],a[jp+p40],a[jt+p20],a[jp+p20], t00,t01,t02,t03,t04,t05, r[30].re,r[30].im,r[62].re,r[62].im,r[94].re,r[94].im);	jt = j1 + p03; jp = j2 + p03;
		RADIX_03_DFT(s,c3m1, a[jt    ],a[jp    ],a[jt+p40],a[jp+p40],a[jt+p20],a[jp+p20], t00,t01,t02,t03,t04,t05, r[31].re,r[31].im,r[63].re,r[63].im,r[95].re,r[95].im);

	/*...and now do 3 radix-32 transforms.
	Use the supercalafragalistic Ancient Chinese Secret index-munging formula [SACSIMPF]
	to properly permute the outputs to the required ordering, which in terms of our 3 radix-32 DFTs is

								[all indices hex here]
	 0, 1, 2, 3, 5, 4, 7, 6, a, b, 9, 8, f, e, c, d,15,14,17,16,12,13,11,10,1f,1e,1c,1d,19,18,1b,1a   [01235476ab98fecd + p00],[54762310fecd98ba + p10]
	4a,4b,49,48,4f,4e,4c,4d,45,44,47,46,42,43,41,40,5f,5e,5c,5d,59,58,5b,5a,52,53,51,50,57,56,54,55 = [ab98fecd54762310 + p40],[fecd98ba23107645 + p50]
	35,34,37,36,32,33,31,30,3f,3e,3c,3d,39,38,3b,3a,2a,2b,29,28,2f,2e,2c,2d,25,24,27,26,22,23,21,20   [54762310fecd98ba + p30],[ab98fecd54762310 + p20]
	*/
		// NOTE: Due to casting-to-double of inputs temp-array pointer Radix-32 macro doesn't play nice with r-array
		// offsets of form "r+32", so make address-taking explicit and wrapping in () prior to macro code execution:
		//
		// Set arr_offsets for radix-32 DFT outputs:	[01235476ab98fecd + p00],[54762310fecd98ba + p10]
		arr_offsets[0x00] = 0  ;		arr_offsets[0x10] = p05+p10;
		arr_offsets[0x01] = p01;		arr_offsets[0x11] = p04+p10;
		arr_offsets[0x02] = p02;		arr_offsets[0x12] = p07+p10;
		arr_offsets[0x03] = p03;		arr_offsets[0x13] = p06+p10;
		arr_offsets[0x04] = p05;		arr_offsets[0x14] = p02+p10;
		arr_offsets[0x05] = p04;		arr_offsets[0x15] = p03+p10;
		arr_offsets[0x06] = p07;		arr_offsets[0x16] = p01+p10;
		arr_offsets[0x07] = p06;		arr_offsets[0x17] =     p10;
		arr_offsets[0x08] = p0a;		arr_offsets[0x18] = p0f+p10;
		arr_offsets[0x09] = p0b;		arr_offsets[0x19] = p0e+p10;
		arr_offsets[0x0A] = p09;		arr_offsets[0x1A] = p0c+p10;
		arr_offsets[0x0B] = p08;		arr_offsets[0x1B] = p0d+p10;
		arr_offsets[0x0C] = p0f;		arr_offsets[0x1C] = p09+p10;
		arr_offsets[0x0D] = p0e;		arr_offsets[0x1D] = p08+p10;
		arr_offsets[0x0E] = p0c;		arr_offsets[0x1E] = p0b+p10;
		arr_offsets[0x0F] = p0d;		arr_offsets[0x1F] = p0a+p10;
		RADIX_32_DIF_OOP((double *)(r+00), a    ,arr_offsets);	/* Inputs in r[00-31] */		// Set arr_offsets for radix-32 DFT outputs:

		// Set arr_offsets for radix-32 DFT outputs:	[ab98fecd54762310 + p40],[fecd98ba23107645 + p50]
		arr_offsets[0x00] = p0a;		arr_offsets[0x10] = p0f+p10;
		arr_offsets[0x01] = p0b;		arr_offsets[0x11] = p0e+p10;
		arr_offsets[0x02] = p09;		arr_offsets[0x12] = p0c+p10;
		arr_offsets[0x03] = p08;		arr_offsets[0x13] = p0d+p10;
		arr_offsets[0x04] = p0f;		arr_offsets[0x14] = p09+p10;
		arr_offsets[0x05] = p0e;		arr_offsets[0x15] = p08+p10;
		arr_offsets[0x06] = p0c;		arr_offsets[0x16] = p0b+p10;
		arr_offsets[0x07] = p0d;		arr_offsets[0x17] = p0a+p10;
		arr_offsets[0x08] = p05;		arr_offsets[0x18] = p02+p10;
		arr_offsets[0x09] = p04;		arr_offsets[0x19] = p03+p10;
		arr_offsets[0x0A] = p07;		arr_offsets[0x1A] = p01+p10;
		arr_offsets[0x0B] = p06;		arr_offsets[0x1B] =     p10;
		arr_offsets[0x0C] = p02;		arr_offsets[0x1C] = p07+p10;
		arr_offsets[0x0D] = p03;		arr_offsets[0x1D] = p06+p10;
		arr_offsets[0x0E] = p01;		arr_offsets[0x1E] = p04+p10;
		arr_offsets[0x0F] = 0  ;		arr_offsets[0x1F] = p05+p10;
		RADIX_32_DIF_OOP((double *)(r+32), a+p40,arr_offsets);	/* Inputs in r[32-63] */		// Set arr_offsets for radix-32 DFT outputs:

		// Set arr_offsets for radix-32 DFT outputs:	[54762310fecd98ba + p30],[ab98fecd54762310 + p20]
		arr_offsets[0x00] = p05+p10;		arr_offsets[0x10] = p0a;
		arr_offsets[0x01] = p04+p10;		arr_offsets[0x11] = p0b;
		arr_offsets[0x02] = p07+p10;		arr_offsets[0x12] = p09;
		arr_offsets[0x03] = p06+p10;		arr_offsets[0x13] = p08;
		arr_offsets[0x04] = p02+p10;		arr_offsets[0x14] = p0f;
		arr_offsets[0x05] = p03+p10;		arr_offsets[0x15] = p0e;
		arr_offsets[0x06] = p01+p10;		arr_offsets[0x16] = p0c;
		arr_offsets[0x07] =     p10;		arr_offsets[0x17] = p0d;
		arr_offsets[0x08] = p0f+p10;		arr_offsets[0x18] = p05;
		arr_offsets[0x09] = p0e+p10;		arr_offsets[0x19] = p04;
		arr_offsets[0x0A] = p0c+p10;		arr_offsets[0x1A] = p07;
		arr_offsets[0x0B] = p0d+p10;		arr_offsets[0x1B] = p06;
		arr_offsets[0x0C] = p09+p10;		arr_offsets[0x1C] = p02;
		arr_offsets[0x0D] = p08+p10;		arr_offsets[0x1D] = p03;
		arr_offsets[0x0E] = p0b+p10;		arr_offsets[0x1E] = p01;
		arr_offsets[0x0F] = p0a+p10;		arr_offsets[0x1F] = 0  ;
		RADIX_32_DIF_OOP((double *)(r+64), a+p20,arr_offsets);	/* Inputs in r[64-95] */
	}
}

/***************/

void radix96_dit_pass1(double a[], int n)
{
/*
!...Acronym: DIT = Decimation In Time
!
!...Subroutine to perform an initial radix-96 complex DIT FFT pass on the data in the length-N real vector A.
!
!   See the documentation in radix16_dif_pass for further details on storage and indexing.
*/
	int j,j1,j2,jp,jt;
	// In order to preserve length-2 numeric-index-offset property here, use hex for digit:
	static int NDIVR,p01,p02,p03,p04,p05,p06,p07,p08,p09,p0a,p0b,p0c,p0d,p0e,p0f,p10,p20,p30,p40,p50, first_entry=TRUE;
	const double c3m1= -1.50000000000000000000, s     = 0.86602540378443864675;	/* cos(twopi/3)-1, sin(twopi/3)		*/
	static int arr_offsets[32];
	// Local storage: We must use an array here because scalars have no guarantees about relative address offsets
	// [and even if those are contiguous-as-hoped-for, they may run in reverse]; Make array type (struct complex)
	// to allow us to use the same offset-indexing as in the original radix-32 in-place DFT macros:
	struct complex r[RADIX];
	double t00,t01,t02,t03,t04,t05;

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
	}

/*...The radix-96 pass is here.	*/

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
	here are the 3 index-offset 16-tets going into the radix-32 DFTs:

	Combined DIT input-scramble array =

										[all indices hex here]
	00,01,03,02,07,06,05,04,0f,0e,0d,0c,0b,0a,09,08,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10   [01327654fedcba98 + p00],[fedcba9876543210 + p10]
	35,34,36,37,31,30,32,33,39,38,3a,3b,3e,3f,3c,3d,25,24,26,27,21,20,22,23,29,28,2a,2b,2e,2f,2c,2d = [5467102398abefcd + p30],[5467102398abefcd + p20]
	4a,4b,48,49,4c,4d,4f,4e,42,43,40,41,44,45,47,46,52,53,50,51,54,55,57,56,5c,5d,5f,5e,58,59,5b,5a   [ab89cdfe23014576 + p40],[23014576cdfe89ba + p50]
	*/
	/*...gather the needed data (96 64-bit complex) and do 3 radix-32 transforms,	*/

		// NOTE: Due to casting-to-double of inputs temp-array pointer Radix-32 macro doesn't play nice with r-array
		// offsets of form "r+32", so make address-taking explicit and wrapping in () prior to macro code execution:
		//
		// Set arr_offsets for radix-32 DFT inputs:	[01327654fedcba98 + p00],[fedcba9876543210 + p10]
		arr_offsets[0x00] =   0;		arr_offsets[0x10] = p0f+p10;
		arr_offsets[0x01] = p01;		arr_offsets[0x11] = p0e+p10;
		arr_offsets[0x02] = p03;		arr_offsets[0x12] = p0d+p10;
		arr_offsets[0x03] = p02;		arr_offsets[0x13] = p0c+p10;
		arr_offsets[0x04] = p07;		arr_offsets[0x14] = p0b+p10;
		arr_offsets[0x05] = p06;		arr_offsets[0x15] = p0a+p10;
		arr_offsets[0x06] = p05;		arr_offsets[0x16] = p09+p10;
		arr_offsets[0x07] = p04;		arr_offsets[0x17] = p08+p10;
		arr_offsets[0x08] = p0f;		arr_offsets[0x18] = p07+p10;
		arr_offsets[0x09] = p0e;		arr_offsets[0x19] = p06+p10;
		arr_offsets[0x0A] = p0d;		arr_offsets[0x1A] = p05+p10;
		arr_offsets[0x0B] = p0c;		arr_offsets[0x1B] = p04+p10;
		arr_offsets[0x0C] = p0b;		arr_offsets[0x1C] = p03+p10;
		arr_offsets[0x0D] = p0a;		arr_offsets[0x1D] = p02+p10;
		arr_offsets[0x0E] = p09;		arr_offsets[0x1E] = p01+p10;
		arr_offsets[0x0F] = p08;		arr_offsets[0x1F] =     p10;
		RADIX_32_DIT_OOP(a    ,arr_offsets, (double *)(r+00));	/* Outputs in r[00-31] */

		// Set arr_offsets for radix-32 DFT inputs:	[5467102398abefcd + p30],[5467102398abefcd + p20]
		arr_offsets[0x00] = p05+p10;		arr_offsets[0x10] = p05;
		arr_offsets[0x01] = p04+p10;		arr_offsets[0x11] = p04;
		arr_offsets[0x02] = p06+p10;		arr_offsets[0x12] = p06;
		arr_offsets[0x03] = p07+p10;		arr_offsets[0x13] = p07;
		arr_offsets[0x04] = p01+p10;		arr_offsets[0x14] = p01;
		arr_offsets[0x05] =     p10;		arr_offsets[0x15] =   0;
		arr_offsets[0x06] = p02+p10;		arr_offsets[0x16] = p02;
		arr_offsets[0x07] = p03+p10;		arr_offsets[0x17] = p03;
		arr_offsets[0x08] = p09+p10;		arr_offsets[0x18] = p09;
		arr_offsets[0x09] = p08+p10;		arr_offsets[0x19] = p08;
		arr_offsets[0x0A] = p0a+p10;		arr_offsets[0x1A] = p0a;
		arr_offsets[0x0B] = p0b+p10;		arr_offsets[0x1B] = p0b;
		arr_offsets[0x0C] = p0e+p10;		arr_offsets[0x1C] = p0e;
		arr_offsets[0x0D] = p0f+p10;		arr_offsets[0x1D] = p0f;
		arr_offsets[0x0E] = p0c+p10;		arr_offsets[0x1E] = p0c;
		arr_offsets[0x0F] = p0d+p10;		arr_offsets[0x1F] = p0d;
		RADIX_32_DIT_OOP(a+p20,arr_offsets, (double *)(r+32));	/* Outputs in r[32-63] */

		// Set arr_offsets for radix-32 DFT inputs:	[ab89cdfe23014576 + p40],[23014576cdfe89ba + p50]
		arr_offsets[0x00] = p0a;		arr_offsets[0x10] = p02+p10;
		arr_offsets[0x01] = p0b;		arr_offsets[0x11] = p03+p10;
		arr_offsets[0x02] = p08;		arr_offsets[0x12] =     p10;
		arr_offsets[0x03] = p09;		arr_offsets[0x13] = p01+p10;
		arr_offsets[0x04] = p0c;		arr_offsets[0x14] = p04+p10;
		arr_offsets[0x05] = p0d;		arr_offsets[0x15] = p05+p10;
		arr_offsets[0x06] = p0f;		arr_offsets[0x16] = p07+p10;
		arr_offsets[0x07] = p0e;		arr_offsets[0x17] = p06+p10;
		arr_offsets[0x08] = p02;		arr_offsets[0x18] = p0c+p10;
		arr_offsets[0x09] = p03;		arr_offsets[0x19] = p0d+p10;
		arr_offsets[0x0A] =   0;		arr_offsets[0x1A] = p0f+p10;
		arr_offsets[0x0B] = p01;		arr_offsets[0x1B] = p0e+p10;
		arr_offsets[0x0C] = p04;		arr_offsets[0x1C] = p08+p10;
		arr_offsets[0x0D] = p05;		arr_offsets[0x1D] = p09+p10;
		arr_offsets[0x0E] = p07;		arr_offsets[0x1E] = p0b+p10;
		arr_offsets[0x0F] = p06;		arr_offsets[0x1F] = p0a+p10;
		RADIX_32_DIT_OOP(a+p40,arr_offsets, (double *)(r+64));	/* Outputs in r[64-95] */

	/*...and now do 32 radix-3 transforms:

	Use the inverse supercalafragalistic Ancient Chinese Secret index-munging formula [iSACSIMPF]
	to properly permute the outputs of the 32 radix-3 DFTs:
		[all indices hex here]
		00,40,20   00,40,20 + p00		30,10,50   30,10,50 + p00
		3f,1f,5f   30,10,50 + p0f		0f,4f,2f   00,40,20 + p0f
		1e,5e,3e   10,50,30 + p0e		4e,2e,0e   40,20,00 + p0e
		5d,3d,1d   50,30,10 + p0d		2d,0d,4d   20,00,40 + p0d
		3c,1c,5c   30,10,50 + p0c		0c,4c,2c   00,40,20 + p0c
		1b,5b,3b   10,50,30 + p0b		4b,2b,0b   40,20,00 + p0b
		5a,3a,1a   50,30,10 + p0a		2a,0a,4a   20,00,40 + p0a
		39,19,59 = 30,10,50 + p09		09,49,29 = 00,40,20 + p09
		18,58,38   10,50,30 + p08		48,28,08   40,20,00 + p08
		57,37,17   50,30,10 + p07		27,07,47   20,00,40 + p07
		36,16,56   30,10,50 + p06		06,46,26   00,40,20 + p06
		15,55,35   10,50,30 + p05		45,25,05   40,20,00 + p05
		54,34,14   50,30,10 + p04		24,04,44   20,00,40 + p04
		33,13,53   30,10,50 + p03		03,43,23   00,40,20 + p03
		12,52,32   10,50,30 + p02		42,22,02   40,20,00 + p02
		51,31,11   50,30,10 + p01		21,01,41   20,00,40 + p01
		[cont. in rcol]
	*/
		RADIX_03_DFT(s,c3m1, r[ 0].re,r[ 0].im,r[32].re,r[32].im,r[64].re,r[64].im, t00,t01,t02,t03,t04,t05, a[j1    ],a[j2    ],a[j1+p40],a[j2+p40],a[j1+p20],a[j2+p20]);	jt = j1 + p0f; jp = j2 + p0f;
		RADIX_03_DFT(s,c3m1, r[ 1].re,r[ 1].im,r[33].re,r[33].im,r[65].re,r[65].im, t00,t01,t02,t03,t04,t05, a[jt+p30],a[jp+p30],a[jt+p10],a[jp+p10],a[jt+p50],a[jp+p50]);	jt = j1 + p0e; jp = j2 + p0e;
		RADIX_03_DFT(s,c3m1, r[ 2].re,r[ 2].im,r[34].re,r[34].im,r[66].re,r[66].im, t00,t01,t02,t03,t04,t05, a[jt+p10],a[jp+p10],a[jt+p50],a[jp+p50],a[jt+p30],a[jp+p30]);	jt = j1 + p0d; jp = j2 + p0d;
		RADIX_03_DFT(s,c3m1, r[ 3].re,r[ 3].im,r[35].re,r[35].im,r[67].re,r[67].im, t00,t01,t02,t03,t04,t05, a[jt+p50],a[jp+p50],a[jt+p30],a[jp+p30],a[jt+p10],a[jp+p10]);	jt = j1 + p0c; jp = j2 + p0c;
		RADIX_03_DFT(s,c3m1, r[ 4].re,r[ 4].im,r[36].re,r[36].im,r[68].re,r[68].im, t00,t01,t02,t03,t04,t05, a[jt+p30],a[jp+p30],a[jt+p10],a[jp+p10],a[jt+p50],a[jp+p50]);	jt = j1 + p0b; jp = j2 + p0b;
		RADIX_03_DFT(s,c3m1, r[ 5].re,r[ 5].im,r[37].re,r[37].im,r[69].re,r[69].im, t00,t01,t02,t03,t04,t05, a[jt+p10],a[jp+p10],a[jt+p50],a[jp+p50],a[jt+p30],a[jp+p30]);	jt = j1 + p0a; jp = j2 + p0a;
		RADIX_03_DFT(s,c3m1, r[ 6].re,r[ 6].im,r[38].re,r[38].im,r[70].re,r[70].im, t00,t01,t02,t03,t04,t05, a[jt+p50],a[jp+p50],a[jt+p30],a[jp+p30],a[jt+p10],a[jp+p10]);	jt = j1 + p09; jp = j2 + p09;
		RADIX_03_DFT(s,c3m1, r[ 7].re,r[ 7].im,r[39].re,r[39].im,r[71].re,r[71].im, t00,t01,t02,t03,t04,t05, a[jt+p30],a[jp+p30],a[jt+p10],a[jp+p10],a[jt+p50],a[jp+p50]);	jt = j1 + p08; jp = j2 + p08;
		RADIX_03_DFT(s,c3m1, r[ 8].re,r[ 8].im,r[40].re,r[40].im,r[72].re,r[72].im, t00,t01,t02,t03,t04,t05, a[jt+p10],a[jp+p10],a[jt+p50],a[jp+p50],a[jt+p30],a[jp+p30]);	jt = j1 + p07; jp = j2 + p07;
		RADIX_03_DFT(s,c3m1, r[ 9].re,r[ 9].im,r[41].re,r[41].im,r[73].re,r[73].im, t00,t01,t02,t03,t04,t05, a[jt+p50],a[jp+p50],a[jt+p30],a[jp+p30],a[jt+p10],a[jp+p10]);	jt = j1 + p06; jp = j2 + p06;
		RADIX_03_DFT(s,c3m1, r[10].re,r[10].im,r[42].re,r[42].im,r[74].re,r[74].im, t00,t01,t02,t03,t04,t05, a[jt+p30],a[jp+p30],a[jt+p10],a[jp+p10],a[jt+p50],a[jp+p50]);	jt = j1 + p05; jp = j2 + p05;
		RADIX_03_DFT(s,c3m1, r[11].re,r[11].im,r[43].re,r[43].im,r[75].re,r[75].im, t00,t01,t02,t03,t04,t05, a[jt+p10],a[jp+p10],a[jt+p50],a[jp+p50],a[jt+p30],a[jp+p30]);	jt = j1 + p04; jp = j2 + p04;
		RADIX_03_DFT(s,c3m1, r[12].re,r[12].im,r[44].re,r[44].im,r[76].re,r[76].im, t00,t01,t02,t03,t04,t05, a[jt+p50],a[jp+p50],a[jt+p30],a[jp+p30],a[jt+p10],a[jp+p10]);	jt = j1 + p03; jp = j2 + p03;
		RADIX_03_DFT(s,c3m1, r[13].re,r[13].im,r[45].re,r[45].im,r[77].re,r[77].im, t00,t01,t02,t03,t04,t05, a[jt+p30],a[jp+p30],a[jt+p10],a[jp+p10],a[jt+p50],a[jp+p50]);	jt = j1 + p02; jp = j2 + p02;
		RADIX_03_DFT(s,c3m1, r[14].re,r[14].im,r[46].re,r[46].im,r[78].re,r[78].im, t00,t01,t02,t03,t04,t05, a[jt+p10],a[jp+p10],a[jt+p50],a[jp+p50],a[jt+p30],a[jp+p30]);	jt = j1 + p01; jp = j2 + p01;
		RADIX_03_DFT(s,c3m1, r[15].re,r[15].im,r[47].re,r[47].im,r[79].re,r[79].im, t00,t01,t02,t03,t04,t05, a[jt+p50],a[jp+p50],a[jt+p30],a[jp+p30],a[jt+p10],a[jp+p10]);
		RADIX_03_DFT(s,c3m1, r[16].re,r[16].im,r[48].re,r[48].im,r[80].re,r[80].im, t00,t01,t02,t03,t04,t05, a[j1+p30],a[j2+p30],a[j1+p10],a[j2+p10],a[j1+p50],a[j2+p50]);	jt = j1 + p0f; jp = j2 + p0f;
		RADIX_03_DFT(s,c3m1, r[17].re,r[17].im,r[49].re,r[49].im,r[81].re,r[81].im, t00,t01,t02,t03,t04,t05, a[jt    ],a[jp    ],a[jt+p40],a[jp+p40],a[jt+p20],a[jp+p20]);	jt = j1 + p0e; jp = j2 + p0e;
		RADIX_03_DFT(s,c3m1, r[18].re,r[18].im,r[50].re,r[50].im,r[82].re,r[82].im, t00,t01,t02,t03,t04,t05, a[jt+p40],a[jp+p40],a[jt+p20],a[jp+p20],a[jt    ],a[jp    ]);	jt = j1 + p0d; jp = j2 + p0d;
		RADIX_03_DFT(s,c3m1, r[19].re,r[19].im,r[51].re,r[51].im,r[83].re,r[83].im, t00,t01,t02,t03,t04,t05, a[jt+p20],a[jp+p20],a[jt    ],a[jp    ],a[jt+p40],a[jp+p40]);	jt = j1 + p0c; jp = j2 + p0c;
		RADIX_03_DFT(s,c3m1, r[20].re,r[20].im,r[52].re,r[52].im,r[84].re,r[84].im, t00,t01,t02,t03,t04,t05, a[jt    ],a[jp    ],a[jt+p40],a[jp+p40],a[jt+p20],a[jp+p20]);	jt = j1 + p0b; jp = j2 + p0b;
		RADIX_03_DFT(s,c3m1, r[21].re,r[21].im,r[53].re,r[53].im,r[85].re,r[85].im, t00,t01,t02,t03,t04,t05, a[jt+p40],a[jp+p40],a[jt+p20],a[jp+p20],a[jt    ],a[jp    ]);	jt = j1 + p0a; jp = j2 + p0a;
		RADIX_03_DFT(s,c3m1, r[22].re,r[22].im,r[54].re,r[54].im,r[86].re,r[86].im, t00,t01,t02,t03,t04,t05, a[jt+p20],a[jp+p20],a[jt    ],a[jp    ],a[jt+p40],a[jp+p40]);	jt = j1 + p09; jp = j2 + p09;
		RADIX_03_DFT(s,c3m1, r[23].re,r[23].im,r[55].re,r[55].im,r[87].re,r[87].im, t00,t01,t02,t03,t04,t05, a[jt    ],a[jp    ],a[jt+p40],a[jp+p40],a[jt+p20],a[jp+p20]);	jt = j1 + p08; jp = j2 + p08;
		RADIX_03_DFT(s,c3m1, r[24].re,r[24].im,r[56].re,r[56].im,r[88].re,r[88].im, t00,t01,t02,t03,t04,t05, a[jt+p40],a[jp+p40],a[jt+p20],a[jp+p20],a[jt    ],a[jp    ]);	jt = j1 + p07; jp = j2 + p07;
		RADIX_03_DFT(s,c3m1, r[25].re,r[25].im,r[57].re,r[57].im,r[89].re,r[89].im, t00,t01,t02,t03,t04,t05, a[jt+p20],a[jp+p20],a[jt    ],a[jp    ],a[jt+p40],a[jp+p40]);	jt = j1 + p06; jp = j2 + p06;
		RADIX_03_DFT(s,c3m1, r[26].re,r[26].im,r[58].re,r[58].im,r[90].re,r[90].im, t00,t01,t02,t03,t04,t05, a[jt    ],a[jp    ],a[jt+p40],a[jp+p40],a[jt+p20],a[jp+p20]);	jt = j1 + p05; jp = j2 + p05;
		RADIX_03_DFT(s,c3m1, r[27].re,r[27].im,r[59].re,r[59].im,r[91].re,r[91].im, t00,t01,t02,t03,t04,t05, a[jt+p40],a[jp+p40],a[jt+p20],a[jp+p20],a[jt    ],a[jp    ]);	jt = j1 + p04; jp = j2 + p04;
		RADIX_03_DFT(s,c3m1, r[28].re,r[28].im,r[60].re,r[60].im,r[92].re,r[92].im, t00,t01,t02,t03,t04,t05, a[jt+p20],a[jp+p20],a[jt    ],a[jp    ],a[jt+p40],a[jp+p40]);	jt = j1 + p03; jp = j2 + p03;
		RADIX_03_DFT(s,c3m1, r[29].re,r[29].im,r[61].re,r[61].im,r[93].re,r[93].im, t00,t01,t02,t03,t04,t05, a[jt    ],a[jp    ],a[jt+p40],a[jp+p40],a[jt+p20],a[jp+p20]);	jt = j1 + p02; jp = j2 + p02;
		RADIX_03_DFT(s,c3m1, r[30].re,r[30].im,r[62].re,r[62].im,r[94].re,r[94].im, t00,t01,t02,t03,t04,t05, a[jt+p40],a[jp+p40],a[jt+p20],a[jp+p20],a[jt    ],a[jp    ]);	jt = j1 + p01; jp = j2 + p01;
		RADIX_03_DFT(s,c3m1, r[31].re,r[31].im,r[63].re,r[63].im,r[95].re,r[95].im, t00,t01,t02,t03,t04,t05, a[jt+p20],a[jp+p20],a[jt    ],a[jp    ],a[jt+p40],a[jp+p40]);
	}
}

#undef RADIX
