/*******************************************************************************
*                                                                              *
*   (C) 1997-2021  by Ernst W. Mayer.                                           *
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
#include "radix31.h"

/***************/

void radix31_dif_pass1(double a[], int n)
{
/*
!
!...Acronym: DIF = Decimation In Frequency
!
!...Subroutine to perform an initial radix-31 complex DIF FFT pass on the data in the length-N real vector A.
!
!   See the documentation in radix16_dif_pass for further details on storage and indexing.
*/
	int j,j1,j2;
	static int arr_offsets[31],
		n31,p01,p02,p03,p04,p05,p06,p07,p08,p09,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,p23,p24,p25,p26,p27,p28,p29,p30,
		first_entry=TRUE;
	double tmp[62];

/*...initialize things upon first entry	*/

	if(!first_entry && (n/31) != n31)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}
	if(first_entry)
	{
		first_entry=FALSE;
		n31=n/31;

/*   constant index offsets for array load/stores are here.	*/

		p01= n31;
		p02= p01+p01;
		p03= p02+p01;
		p04= p03+p01;
		p05= p04+p01;
		p06= p05+p01;
		p07= p06+p01;
		p08= p07+p01;
		p09= p08+p01;
		p10= p09+p01;
		p11= p10+p01;
		p12= p11+p01;
		p13= p12+p01;
		p14= p13+p01;
		p15= p14+p01;
		p16= p15+p01;
		p17= p16+p01;
		p18= p17+p01;
		p19= p18+p01;
		p20= p19+p01;
		p21= p20+p01;
		p22= p21+p01;
		p23= p22+p01;
		p24= p23+p01;
		p25= p24+p01;
		p26= p25+p01;
		p27= p26+p01;
		p28= p27+p01;
		p29= p28+p01;
		p30= p29+p01;

		p01= p01+ ( (p01>> DAT_BITS) << PAD_BITS );
		p02= p02+ ( (p02>> DAT_BITS) << PAD_BITS );
		p03= p03+ ( (p03>> DAT_BITS) << PAD_BITS );
		p04= p04+ ( (p04>> DAT_BITS) << PAD_BITS );
		p05= p05+ ( (p05>> DAT_BITS) << PAD_BITS );
		p06= p06+ ( (p06>> DAT_BITS) << PAD_BITS );
		p07= p07+ ( (p07>> DAT_BITS) << PAD_BITS );
		p08= p08+ ( (p08>> DAT_BITS) << PAD_BITS );
		p09= p09+ ( (p09>> DAT_BITS) << PAD_BITS );
		p10= p10+ ( (p10>> DAT_BITS) << PAD_BITS );
		p11= p11+ ( (p11>> DAT_BITS) << PAD_BITS );
		p12= p12+ ( (p12>> DAT_BITS) << PAD_BITS );
		p13= p13+ ( (p13>> DAT_BITS) << PAD_BITS );
		p14= p14+ ( (p14>> DAT_BITS) << PAD_BITS );
		p15= p15+ ( (p15>> DAT_BITS) << PAD_BITS );
		p16= p16+ ( (p16>> DAT_BITS) << PAD_BITS );
		p17= p17+ ( (p17>> DAT_BITS) << PAD_BITS );
		p18= p18+ ( (p18>> DAT_BITS) << PAD_BITS );
		p19= p19+ ( (p19>> DAT_BITS) << PAD_BITS );
		p20= p20+ ( (p20>> DAT_BITS) << PAD_BITS );
		p21= p21+ ( (p21>> DAT_BITS) << PAD_BITS );
		p22= p22+ ( (p22>> DAT_BITS) << PAD_BITS );
		p23= p23+ ( (p23>> DAT_BITS) << PAD_BITS );
		p24= p24+ ( (p24>> DAT_BITS) << PAD_BITS );
		p25= p25+ ( (p25>> DAT_BITS) << PAD_BITS );
		p26= p26+ ( (p26>> DAT_BITS) << PAD_BITS );
		p27= p27+ ( (p27>> DAT_BITS) << PAD_BITS );
		p28= p28+ ( (p28>> DAT_BITS) << PAD_BITS );
		p29= p29+ ( (p29>> DAT_BITS) << PAD_BITS );
		p30= p30+ ( (p30>> DAT_BITS) << PAD_BITS );

		arr_offsets[0x00] = 0;
		arr_offsets[0x01] = p01;
		arr_offsets[0x02] = p02;
		arr_offsets[0x03] = p03;
		arr_offsets[0x04] = p04;
		arr_offsets[0x05] = p05;
		arr_offsets[0x06] = p06;
		arr_offsets[0x07] = p07;
		arr_offsets[0x08] = p08;
		arr_offsets[0x09] = p09;
		arr_offsets[0x0A] = p10;
		arr_offsets[0x0B] = p11;
		arr_offsets[0x0C] = p12;
		arr_offsets[0x0D] = p13;
		arr_offsets[0x0E] = p14;
		arr_offsets[0x0F] = p15;
		arr_offsets[0x10] = p16;
		arr_offsets[0x11] = p17;
		arr_offsets[0x12] = p18;
		arr_offsets[0x13] = p19;
		arr_offsets[0x14] = p20;
		arr_offsets[0x15] = p21;
		arr_offsets[0x16] = p22;
		arr_offsets[0x17] = p23;
		arr_offsets[0x18] = p24;
		arr_offsets[0x19] = p25;
		arr_offsets[0x1A] = p26;
		arr_offsets[0x1B] = p27;
		arr_offsets[0x1C] = p28;
		arr_offsets[0x1D] = p29;
		arr_offsets[0x1E] = p30;
	}

/*...The radix-31 pass is here.	*/

	for(j=0; j < n31; j += 2)
	{
	#ifdef USE_SSE2
		j1 = (j & mask01) + br4[j&3];
		j1 =j1 + ( (j1>> DAT_BITS) << PAD_BITS );
	#else
		j1 = j + ( (j >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
	#endif
		j2 = j1+RE_IM_STRIDE;

	#if 0
		RADIX_31_DFT(a+j1,0,p01,p02,p03,p04,p05,p06,p07,p08,p09,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,p23,p24,p25,p26,p27,p28,p29,p30
	/* outputs: */	,a+j1,0,p01,p02,p03,p04,p05,p06,p07,p08,p09,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,p23,p24,p25,p26,p27,p28,p29,p30);
	#else	/* Out-of-place version, needed for combining radix-31 with powers-of-2: */
		RADIX_31_DIF(\
			a+j1,arr_offsets,\
			tmp\
		);
		a[j1    ] = tmp[ 0];	a[j2    ] = tmp[ 1];
		a[j1+p01] = tmp[ 2];	a[j2+p01] = tmp[ 3];
		a[j1+p02] = tmp[ 4];	a[j2+p02] = tmp[ 5];
		a[j1+p03] = tmp[ 6];	a[j2+p03] = tmp[ 7];
		a[j1+p04] = tmp[ 8];	a[j2+p04] = tmp[ 9];
		a[j1+p05] = tmp[10];	a[j2+p05] = tmp[11];
		a[j1+p06] = tmp[12];	a[j2+p06] = tmp[13];
		a[j1+p07] = tmp[14];	a[j2+p07] = tmp[15];
		a[j1+p08] = tmp[16];	a[j2+p08] = tmp[17];
		a[j1+p09] = tmp[18];	a[j2+p09] = tmp[19];
		a[j1+p10] = tmp[20];	a[j2+p10] = tmp[21];
		a[j1+p11] = tmp[22];	a[j2+p11] = tmp[23];
		a[j1+p12] = tmp[24];	a[j2+p12] = tmp[25];
		a[j1+p13] = tmp[26];	a[j2+p13] = tmp[27];
		a[j1+p14] = tmp[28];	a[j2+p14] = tmp[29];
		a[j1+p15] = tmp[30];	a[j2+p15] = tmp[31];
		a[j1+p16] = tmp[32];	a[j2+p16] = tmp[33];
		a[j1+p17] = tmp[34];	a[j2+p17] = tmp[35];
		a[j1+p18] = tmp[36];	a[j2+p18] = tmp[37];
		a[j1+p19] = tmp[38];	a[j2+p19] = tmp[39];
		a[j1+p20] = tmp[40];	a[j2+p20] = tmp[41];
		a[j1+p21] = tmp[42];	a[j2+p21] = tmp[43];
		a[j1+p22] = tmp[44];	a[j2+p22] = tmp[45];
		a[j1+p23] = tmp[46];	a[j2+p23] = tmp[47];
		a[j1+p24] = tmp[48];	a[j2+p24] = tmp[49];
		a[j1+p25] = tmp[50];	a[j2+p25] = tmp[51];
		a[j1+p26] = tmp[52];	a[j2+p26] = tmp[53];
		a[j1+p27] = tmp[54];	a[j2+p27] = tmp[55];
		a[j1+p28] = tmp[56];	a[j2+p28] = tmp[57];
		a[j1+p29] = tmp[58];	a[j2+p29] = tmp[59];
		a[j1+p30] = tmp[60];	a[j2+p30] = tmp[61];
	#endif
	}
}

/*
DIT - Since the radix is prime - differs from DIF only in having the + and - of the final output-combinations swapped.
*/
void radix31_dit_pass1(double a[], int n)
{
	int j,j1;
	static int n31,p01,p02,p03,p04,p05,p06,p07,p08,p09,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,p23,p24,p25,p26,p27,p28,p29,p30,
		first_entry=TRUE;

	if(!first_entry && (n/31) != n31)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		first_entry=FALSE;
		n31=n/31;

/*   constant index offsets for array load/stores are here.	*/

		p01= n31;
		p02= p01+p01;
		p03= p02+p01;
		p04= p03+p01;
		p05= p04+p01;
		p06= p05+p01;
		p07= p06+p01;
		p08= p07+p01;
		p09= p08+p01;
		p10= p09+p01;
		p11= p10+p01;
		p12= p11+p01;
		p13= p12+p01;
		p14= p13+p01;
		p15= p14+p01;
		p16= p15+p01;
		p17= p16+p01;
		p18= p17+p01;
		p19= p18+p01;
		p20= p19+p01;
		p21= p20+p01;
		p22= p21+p01;
		p23= p22+p01;
		p24= p23+p01;
		p25= p24+p01;
		p26= p25+p01;
		p27= p26+p01;
		p28= p27+p01;
		p29= p28+p01;
		p30= p29+p01;

		p01= p01+ ( (p01>> DAT_BITS) << PAD_BITS );
		p02= p02+ ( (p02>> DAT_BITS) << PAD_BITS );
		p03= p03+ ( (p03>> DAT_BITS) << PAD_BITS );
		p04= p04+ ( (p04>> DAT_BITS) << PAD_BITS );
		p05= p05+ ( (p05>> DAT_BITS) << PAD_BITS );
		p06= p06+ ( (p06>> DAT_BITS) << PAD_BITS );
		p07= p07+ ( (p07>> DAT_BITS) << PAD_BITS );
		p08= p08+ ( (p08>> DAT_BITS) << PAD_BITS );
		p09= p09+ ( (p09>> DAT_BITS) << PAD_BITS );
		p10= p10+ ( (p10>> DAT_BITS) << PAD_BITS );
		p11= p11+ ( (p11>> DAT_BITS) << PAD_BITS );
		p12= p12+ ( (p12>> DAT_BITS) << PAD_BITS );
		p13= p13+ ( (p13>> DAT_BITS) << PAD_BITS );
		p14= p14+ ( (p14>> DAT_BITS) << PAD_BITS );
		p15= p15+ ( (p15>> DAT_BITS) << PAD_BITS );
		p16= p16+ ( (p16>> DAT_BITS) << PAD_BITS );
		p17= p17+ ( (p17>> DAT_BITS) << PAD_BITS );
		p18= p18+ ( (p18>> DAT_BITS) << PAD_BITS );
		p19= p19+ ( (p19>> DAT_BITS) << PAD_BITS );
		p20= p20+ ( (p20>> DAT_BITS) << PAD_BITS );
		p21= p21+ ( (p21>> DAT_BITS) << PAD_BITS );
		p22= p22+ ( (p22>> DAT_BITS) << PAD_BITS );
		p23= p23+ ( (p23>> DAT_BITS) << PAD_BITS );
		p24= p24+ ( (p24>> DAT_BITS) << PAD_BITS );
		p25= p25+ ( (p25>> DAT_BITS) << PAD_BITS );
		p26= p26+ ( (p26>> DAT_BITS) << PAD_BITS );
		p27= p27+ ( (p27>> DAT_BITS) << PAD_BITS );
		p28= p28+ ( (p28>> DAT_BITS) << PAD_BITS );
		p29= p29+ ( (p29>> DAT_BITS) << PAD_BITS );
		p30= p30+ ( (p30>> DAT_BITS) << PAD_BITS );
	}

/*...The radix-31 pass is here.	*/

	for(j=0; j < n31; j += 2)
	{
	#ifdef USE_SSE2
		j1 = (j & mask01) + br4[j&3];
		j1 =j1 + ( (j1>> DAT_BITS) << PAD_BITS );
	#else
		j1 = j + ( (j >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
	#endif

		RADIX_31_DFT(a+j1,0,p01,p02,p03,p04,p05,p06,p07,p08,p09,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,p23,p24,p25,p26,p27,p28,p29,p30
	/* outputs: */	,a+j1,0,p30,p29,p28,p27,p26,p25,p24,p23,p22,p21,p20,p19,p18,p17,p16,p15,p14,p13,p12,p11,p10,p09,p08,p07,p06,p05,p04,p03,p02,p01);
	}
}

