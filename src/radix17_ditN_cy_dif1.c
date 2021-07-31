/*******************************************************************************
*                                                                              *
*   (C) 1997-2020 by Ernst W. Mayer.                                           *
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
#include "radix17_dft.h"

/***************/

int radix17_ditN_cy_dif1(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], double base[], double baseinv[], int iter, double *fracmax, uint64 p)
{
	return 1;
}

/***************/

void radix17_dif_pass1(double a[], int n)
{
/*
!...Acronym: DIF = Decimation In Frequency
!
!...Subroutine to perform an initial radix-17 complex DIF FFT pass on the data in the length-N real vector A.
!
!   See the documentation in radix16_dif_pass for further details on storage and indexing.
!
!   Given complex inputs (x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,xA,xB,xC,xD,xE,xF,xG), we need the following outputs
!   (here cJ = cos(2*J*pi/17), sJ = sin(2*J*pi/17)):
!
!	X0 = C0,      where C0 = x0+   (x1+xG)+   (x2+xF)+   (x3+xE)+   (x4+xD)+   (x5+xC)+   (x6+xB)+   (x7+xA)+   (x6+x9),
!					the cosine terms below get massaged into the form of a length-8 cyclic convolution:
!	X1 = C1 + I*S1		C1 =
!	X2 = C2 + I*S2
!	X3 = C3 + I*S3
!	X4 = C4 + I*S4
!	X5 = C5 + I*S5
!	X6 = C6 + I*S6
!	X7 = C7 + I*S7
!	X8 = C8 + I*S8
!					and the sine terms get massaged into the form of a length-8 acyclic convolution:
!	X9 = C8 - I*S8
!	XA = C7 - I*S7
!	XB = C6 - I*S6
!	XC = C5 - I*S5
!	XD = C4 - I*S4
!	XE = C3 - I*S3
!	XF = C2 - I*S2
!	XG = C1 - I*S1
!
!   We refer to the terms C1-8 (which do not explicitly involving the imaginary constant I)
!   as the "cosine part" of the output, and S1-8 (those multiplied by I) as the "sine part."
!	Opcount for general odd-prime radix R:
!   Totals :                                                        100 FMUL, 140 FADD,		(R-1)^2 fmul	(R+3)*(R-1) fadd
!                                                        compared to 16 FMUL,  96 FADD for radix-12. (Ouch!)
!
!   Relative cost := #FADD/(radix*lg2(radix)) = 3.679 .
*/
	int j,j1,j2;
	static int n17,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16, first_entry=TRUE;

	if(!first_entry && (n/17) != n17)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		first_entry=FALSE;
		n17 = n/17;
	// Constant index offsets for array load/stores are here:
		p1  = n17;
		p2  = p1 +p1;
		p3  = p2 +p1;
		p4  = p3 +p1;
		p5  = p4 +p1;
		p6  = p5 +p1;
		p7  = p6 +p1;
		p8  = p7 +p1;
		p9  = p8 +p1;
		p10 = p9 +p1;
		p11 = p10+p1;
		p12 = p11+p1;
		p13 = p12+p1;
		p14 = p13+p1;
		p15 = p14+p1;
		p16 = p15+p1;

		p1  += ( (p1 >> DAT_BITS) << PAD_BITS );
		p2  += ( (p2 >> DAT_BITS) << PAD_BITS );
		p3  += ( (p3 >> DAT_BITS) << PAD_BITS );
		p4  += ( (p4 >> DAT_BITS) << PAD_BITS );
		p5  += ( (p5 >> DAT_BITS) << PAD_BITS );
		p6  += ( (p6 >> DAT_BITS) << PAD_BITS );
		p7  += ( (p7 >> DAT_BITS) << PAD_BITS );
		p8  += ( (p8 >> DAT_BITS) << PAD_BITS );
		p9  += ( (p9 >> DAT_BITS) << PAD_BITS );
		p10 += ( (p10>> DAT_BITS) << PAD_BITS );
		p11 += ( (p11>> DAT_BITS) << PAD_BITS );
		p12 += ( (p12>> DAT_BITS) << PAD_BITS );
		p13 += ( (p13>> DAT_BITS) << PAD_BITS );
		p14 += ( (p14>> DAT_BITS) << PAD_BITS );
		p15 += ( (p15>> DAT_BITS) << PAD_BITS );
		p16 += ( (p16>> DAT_BITS) << PAD_BITS );
	}

/*...The radix-17 pass is here.	*/

	for(j=0; j < n17; j += 2)
	{
	#ifdef USE_SSE2
		j1 = (j & mask01) + br4[j&3];
			j1 =j1 + ( (j1>> DAT_BITS) << PAD_BITS );
	#else
			j1 = j + ( (j >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
	#endif
		j2 = j1+RE_IM_STRIDE;
		/* Call same radix-11 DFT macro as for DIF, but replace indices [0,1,2,3,4,5,6,7,8,9,10] with j*10%11, j = 0, ..., 10: */
		RADIX_17_DFT(a[j1],a[j2],a[j1+p1],a[j2+p1],a[j1+p2],a[j2+p2],a[j1+p3],a[j2+p3],a[j1+p4],a[j2+p4],a[j1+p5],a[j2+p5],a[j1+p6],a[j2+p6],a[j1+p7],a[j2+p7],a[j1+p8],a[j2+p8],a[j1+p9],a[j2+p9],a[j1+p10],a[j2+p10],a[j1+p11],a[j2+p11],a[j1+p12],a[j2+p12],a[j1+p13],a[j2+p13],a[j1+p14],a[j2+p14],a[j1+p15],a[j2+p15],a[j1+p16],a[j2+p16]
					,a+j1 ,a+j2 ,a+j1+p1 ,a+j2+p1 ,a+j1+p2 ,a+j2+p2 ,a+j1+p3 ,a+j2+p3 ,a+j1+p4 ,a+j2+p4 ,a+j1+p5 ,a+j2+p5 ,a+j1+p6 ,a+j2+p6 ,a+j1+p7 ,a+j2+p7 ,a+j1+p8 ,a+j2+p8 ,a+j1+p9 ,a+j2+p9 ,a+j1+p10 ,a+j2+p10 ,a+j1+p11 ,a+j2+p11 ,a+j1+p12 ,a+j2+p12 ,a+j1+p13 ,a+j2+p13 ,a+j1+p14 ,a+j2+p14 ,a+j1+p15 ,a+j2+p15 ,a+j1+p16 ,a+j2+p16 );
	}
}

/***************/

void radix17_dit_pass1(double a[], int n)
{
/*
!...Acronym: DIT = Decimation In Time
!
!...Subroutine to perform a final radix-17 complex DIT FFT pass on the data in the length-N real vector A.
!
!   See the documentation in radix16_dif_pass for further details on storage and indexing.
*/
	int j,j1,j2;
	static int n17,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16, first_entry=TRUE;

	if(!first_entry && (n/17) != n17)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		first_entry=FALSE;
		n17 = n/17;
	// Constant index offsets for array load/stores are here:
		p1  = n17;
		p2  = p1 +p1;
		p3  = p2 +p1;
		p4  = p3 +p1;
		p5  = p4 +p1;
		p6  = p5 +p1;
		p7  = p6 +p1;
		p8  = p7 +p1;
		p9  = p8 +p1;
		p10 = p9 +p1;
		p11 = p10+p1;
		p12 = p11+p1;
		p13 = p12+p1;
		p14 = p13+p1;
		p15 = p14+p1;
		p16 = p15+p1;

		p1  += ( (p1 >> DAT_BITS) << PAD_BITS );
		p2  += ( (p2 >> DAT_BITS) << PAD_BITS );
		p3  += ( (p3 >> DAT_BITS) << PAD_BITS );
		p4  += ( (p4 >> DAT_BITS) << PAD_BITS );
		p5  += ( (p5 >> DAT_BITS) << PAD_BITS );
		p6  += ( (p6 >> DAT_BITS) << PAD_BITS );
		p7  += ( (p7 >> DAT_BITS) << PAD_BITS );
		p8  += ( (p8 >> DAT_BITS) << PAD_BITS );
		p9  += ( (p9 >> DAT_BITS) << PAD_BITS );
		p10 += ( (p10>> DAT_BITS) << PAD_BITS );
		p11 += ( (p11>> DAT_BITS) << PAD_BITS );
		p12 += ( (p12>> DAT_BITS) << PAD_BITS );
		p13 += ( (p13>> DAT_BITS) << PAD_BITS );
		p14 += ( (p14>> DAT_BITS) << PAD_BITS );
		p15 += ( (p15>> DAT_BITS) << PAD_BITS );
		p16 += ( (p16>> DAT_BITS) << PAD_BITS );
	}

/*...The radix-17 pass is here.	*/

	for(j=0; j < n17; j += 2)
	{
	#ifdef USE_SSE2
		j1 = (j & mask01) + br4[j&3];
			j1 =j1 + ( (j1>> DAT_BITS) << PAD_BITS );
	#else
			j1 = j + ( (j >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
	#endif
		j2 = j1+RE_IM_STRIDE;
		// Call same radix-17 DFT macro as for DIF, but replace indices j = 1-16 with j*16%17, i.e. run in reverse order:
		RADIX_17_DFT(a[j1],a[j2],a[j1+p1],a[j2+p1],a[j1+p2],a[j2+p2],a[j1+p3],a[j2+p3],a[j1+p4],a[j2+p4],a[j1+p5],a[j2+p5],a[j1+p6],a[j2+p6],a[j1+p7],a[j2+p7],a[j1+p8],a[j2+p8],a[j1+p9],a[j2+p9],a[j1+p10],a[j2+p10],a[j1+p11],a[j2+p11],a[j1+p12],a[j2+p12],a[j1+p13],a[j2+p13],a[j1+p14],a[j2+p14],a[j1+p15],a[j2+p15],a[j1+p16],a[j2+p16]
					,a+j1 ,a+j2 ,a+j1+p16 ,a+j2+p16 ,a+j1+p15 ,a+j2+p15 ,a+j1+p14 ,a+j2+p14 ,a+j1+p13 ,a+j2+p13 ,a+j1+p12 ,a+j2+p12 ,a+j1+p11 ,a+j2+p11 ,a+j1+p10 ,a+j2+p10 ,a+j1+p9 ,a+j2+p9 ,a+j1+p8 ,a+j2+p8 ,a+j1+p7 ,a+j2+p7 ,a+j1+p6 ,a+j2+p6 ,a+j1+p5 ,a+j2+p5 ,a+j1+p4 ,a+j2+p4 ,a+j1+p3 ,a+j2+p3 ,a+j1+p2 ,a+j2+p2 ,a+j1+p1 ,a+j2+p1 );
	}
}
