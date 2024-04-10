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
#include "radix512.h"

#define RADIX 512	// Use #define rather than const int to ensure it's really a compile-time const in the C sense

/***************/

int radix512_ditN_cy_dif1(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], struct complex rn0[], struct complex rn1[], double base[], double baseinv[], int iter, double *fracmax, uint64 p)
{
	return 0;
}

/****************/

void radix512_dif_pass1(double a[], int n)
{
/*
!...Acronym: DIF = Decimation In Frequency
!
!...Subroutine to perform an initial radix-512 complex DIF FFT pass on the data in the length-N real vector A.
!
!   The data are stored in a 1-D zero-offset array, with 2^PAD_BITS 8-byte padding elements inserted
!   between every block of 2^DAT_BITS contiguous data. The array padding is to prevent data being accessed
!   in strides that are large powers of two and thus to minimize cache thrashing
!   (in cache-based microprocessor architectures) or bank conflicts (in supercomputers.)
!
!   See the documentation in radix64_dif_pass for further details on storage, indexing and DFT decomposition.
*/
	int i,j,j1,j2,jt,jp;
	static int NDIVR,first_entry=TRUE,
			p01,p02,p03,p04,p05,p06,p07,p08,p09,p0a,p0b,p0c,p0d,p0e,p0f,
			p10,p20,p30,p40,p50,p60,p70,p80,p90,pa0,pb0,pc0,pd0,pe0,pf0,p100;
	static int i_offsets[32], o_offsets[32], po_br[16];
	// Local storage: We must use an array here because scalars have no guarantees about relative address offsets
	// [and even if those are contiguous-as-hoped-for, they may run in reverse]; Make array type (struct complex)
	// to allow us to use the same offset-indexing as in the original radix-32 in-place DFT macros:
	double *addr,*addi;
	#include "radix1024_twiddles.h"	// Can share radix-1024 table, just use first 31 of 63 rows here
	struct complex t[RADIX], *tptr;

	if(!first_entry && (n >> 9) != NDIVR)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		ASSERT((double *)t == &(t[0x00].re), "Unexpected value for Tmp-array-start pointer!");
		first_entry=FALSE;
		NDIVR = n >> 9;

		p01 = NDIVR;
		p02 = p01 + NDIVR;		p01 = p01 + ( (p01 >> DAT_BITS) << PAD_BITS );
		p03 = p02 + NDIVR;		p02 = p02 + ( (p02 >> DAT_BITS) << PAD_BITS );
		p04 = p03 + NDIVR;		p03 = p03 + ( (p03 >> DAT_BITS) << PAD_BITS );
		p05 = p04 + NDIVR;		p04 = p04 + ( (p04 >> DAT_BITS) << PAD_BITS );
		p06 = p05 + NDIVR;		p05 = p05 + ( (p05 >> DAT_BITS) << PAD_BITS );
		p07 = p06 + NDIVR;		p06 = p06 + ( (p06 >> DAT_BITS) << PAD_BITS );
		p08 = p07 + NDIVR;		p07 = p07 + ( (p07 >> DAT_BITS) << PAD_BITS );
		p09 = p08 + NDIVR;		p08 = p08 + ( (p08 >> DAT_BITS) << PAD_BITS );
		p0a = p09 + NDIVR;		p09 = p09 + ( (p09 >> DAT_BITS) << PAD_BITS );
		p0b = p0a + NDIVR;		p0a = p0a + ( (p0a >> DAT_BITS) << PAD_BITS );
		p0c = p0b + NDIVR;		p0b = p0b + ( (p0b >> DAT_BITS) << PAD_BITS );
		p0d = p0c + NDIVR;		p0c = p0c + ( (p0c >> DAT_BITS) << PAD_BITS );
		p0e = p0d + NDIVR;		p0d = p0d + ( (p0d >> DAT_BITS) << PAD_BITS );
		p0f = p0e + NDIVR;		p0e = p0e + ( (p0e >> DAT_BITS) << PAD_BITS );
		NDIVR <<= 4;			p0f = p0f + ( (p0f >> DAT_BITS) << PAD_BITS );
		p10 = NDIVR;
		p20 = p10 + NDIVR;		p10 = p10 + ( (p10 >> DAT_BITS) << PAD_BITS );
		p30 = p20 + NDIVR;		p20 = p20 + ( (p20 >> DAT_BITS) << PAD_BITS );
		p40 = p30 + NDIVR;		p30 = p30 + ( (p30 >> DAT_BITS) << PAD_BITS );
		p50 = p40 + NDIVR;		p40 = p40 + ( (p40 >> DAT_BITS) << PAD_BITS );
		p60 = p50 + NDIVR;		p50 = p50 + ( (p50 >> DAT_BITS) << PAD_BITS );
		p70 = p60 + NDIVR;		p60 = p60 + ( (p60 >> DAT_BITS) << PAD_BITS );
		p80 = p70 + NDIVR;		p70 = p70 + ( (p70 >> DAT_BITS) << PAD_BITS );
		p90 = p80 + NDIVR;		p80 = p80 + ( (p80 >> DAT_BITS) << PAD_BITS );
		pa0 = p90 + NDIVR;		p90 = p90 + ( (p90 >> DAT_BITS) << PAD_BITS );
		pb0 = pa0 + NDIVR;		pa0 = pa0 + ( (pa0 >> DAT_BITS) << PAD_BITS );
		pc0 = pb0 + NDIVR;		pb0 = pb0 + ( (pb0 >> DAT_BITS) << PAD_BITS );
		pd0 = pc0 + NDIVR;		pc0 = pc0 + ( (pc0 >> DAT_BITS) << PAD_BITS );
		pe0 = pd0 + NDIVR;		pd0 = pd0 + ( (pd0 >> DAT_BITS) << PAD_BITS );
		pf0 = pe0 + NDIVR;		pe0 = pe0 + ( (pe0 >> DAT_BITS) << PAD_BITS );
		p100= pf0 + NDIVR;		pf0 = pf0 + ( (pf0 >> DAT_BITS) << PAD_BITS );
		NDIVR >>= 4;			p100= p100+ ( (p100>> DAT_BITS) << PAD_BITS );

		// Set array offsets for radix-32 DFT in/outputs:
		i_offsets[0x00] = 0  ;		i_offsets[0x10] = 0  +p100;
		i_offsets[0x01] = p10;		i_offsets[0x11] = p10+p100;
		i_offsets[0x02] = p20;		i_offsets[0x12] = p20+p100;
		i_offsets[0x03] = p30;		i_offsets[0x13] = p30+p100;
		i_offsets[0x04] = p40;		i_offsets[0x14] = p40+p100;
		i_offsets[0x05] = p50;		i_offsets[0x15] = p50+p100;
		i_offsets[0x06] = p60;		i_offsets[0x16] = p60+p100;
		i_offsets[0x07] = p70;		i_offsets[0x17] = p70+p100;
		i_offsets[0x08] = p80;		i_offsets[0x18] = p80+p100;
		i_offsets[0x09] = p90;		i_offsets[0x19] = p90+p100;
		i_offsets[0x0a] = pa0;		i_offsets[0x1a] = pa0+p100;
		i_offsets[0x0b] = pb0;		i_offsets[0x1b] = pb0+p100;
		i_offsets[0x0c] = pc0;		i_offsets[0x1c] = pc0+p100;
		i_offsets[0x0d] = pd0;		i_offsets[0x1d] = pd0+p100;
		i_offsets[0x0e] = pe0;		i_offsets[0x1e] = pe0+p100;
		i_offsets[0x0f] = pf0;		i_offsets[0x1f] = pf0+p100;

		o_offsets[0x00] = 0x00<<1;	o_offsets[0x10] = 0x10<<1;
		o_offsets[0x01] = 0x01<<1;	o_offsets[0x11] = 0x11<<1;
		o_offsets[0x02] = 0x02<<1;	o_offsets[0x12] = 0x12<<1;
		o_offsets[0x03] = 0x03<<1;	o_offsets[0x13] = 0x13<<1;
		o_offsets[0x04] = 0x04<<1;	o_offsets[0x14] = 0x14<<1;
		o_offsets[0x05] = 0x05<<1;	o_offsets[0x15] = 0x15<<1;
		o_offsets[0x06] = 0x06<<1;	o_offsets[0x16] = 0x16<<1;
		o_offsets[0x07] = 0x07<<1;	o_offsets[0x17] = 0x17<<1;
		o_offsets[0x08] = 0x08<<1;	o_offsets[0x18] = 0x18<<1;
		o_offsets[0x09] = 0x09<<1;	o_offsets[0x19] = 0x19<<1;
		o_offsets[0x0a] = 0x0a<<1;	o_offsets[0x1a] = 0x1a<<1;
		o_offsets[0x0b] = 0x0b<<1;	o_offsets[0x1b] = 0x1b<<1;
		o_offsets[0x0c] = 0x0c<<1;	o_offsets[0x1c] = 0x1c<<1;
		o_offsets[0x0d] = 0x0d<<1;	o_offsets[0x1d] = 0x1d<<1;
		o_offsets[0x0e] = 0x0e<<1;	o_offsets[0x1e] = 0x1e<<1;
		o_offsets[0x0f] = 0x0f<<1;	o_offsets[0x1f] = 0x1f<<1;

		po_br[0x0] = 0; po_br[0x1] = p08; po_br[0x2] = p04; po_br[0x3] = p0c; po_br[0x4] = p02; po_br[0x5] = p0a; po_br[0x6] = p06; po_br[0x7] = p0e; po_br[0x8] = p01; po_br[0x9] = p09; po_br[0xa] = p05; po_br[0xb] = p0d; po_br[0xc] = p03; po_br[0xd] = p0b; po_br[0xe] = p07; po_br[0xf] = p0f;
	}

/*...The radix-512 pass is here.	*/

	for(j = 0; j < NDIVR; j += 2)
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

	// Gather the needed data and do 16 twiddleless length-32 subtransforms, with p-offsets in br16 order: 084c2a6e195d3b7f:
	// NOTE that RADIX_32_DIF outputs are IN-ORDER rather than BR:

		for(i = 0, jp = 0; i < 16; i++, jp += 32) {
			jt = j1 + po_br[i];	// po_br[] = p[084c2a6e195d3b7f]
			//	NOTE that RADIX_64_DIF outputs are IN-ORDER rather than BR:
			RADIX_32_DIF((a+jt),i_offsets,RE_IM_STRIDE, (double *)(t+jp),o_offsets,1);
		}

	/*...and now do 32 radix-16 subtransforms, including the internal twiddles (here use E^n = -E^(n-0x100) for "easy mod" of larger powers):

		Block  0: twiddles = {   1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1}
		Block  1: twiddles = {   1, E^ 1, E^ 2, E^ 3, E^ 4, E^ 5, E^ 6, E^ 7, E^ 8, E^ 9, E^ a, E^ b, E^ c, E^ d, E^ e, E^ f}
		Block  2: twiddles = {   1, E^ 2, E^ 4, E^ 6, E^ 8, E^ a, E^ c, E^ e, E^10, E^12, E^14, E^16, E^18, E^1a, E^1c, E^1e}
		Block  3: twiddles = {   1, E^ 3, E^ 6, E^ 9, E^ c, E^ f, E^12, E^15, E^18, E^1b, E^1e, E^21, E^24, E^27, E^2a, E^2d}
		Block  4: twiddles = {   1, E^ 4, E^ 8, E^ c, E^10, E^14, E^18, E^1c, E^20, E^24, E^28, E^2c, E^30, E^34, E^38, E^3c}
		Block  5: twiddles = {   1, E^ 5, E^ a, E^ f, E^14, E^19, E^1e, E^23, E^28, E^2d, E^32, E^37, E^3c, E^41, E^46, E^4b}
		Block  6: twiddles = {   1, E^ 6, E^ c, E^12, E^18, E^1e, E^24, E^2a, E^30, E^36, E^3c, E^42, E^48, E^4e, E^54, E^5a}
		Block  7: twiddles = {   1, E^ 7, E^ e, E^15, E^1c, E^23, E^2a, E^31, E^38, E^3f, E^46, E^4d, E^54, E^5b, E^62, E^69}
		Block  8: twiddles = {   1, E^ 8, E^10, E^18, E^20, E^28, E^30, E^38, E^40, E^48, E^50, E^58, E^60, E^68, E^70, E^78}
		Block  9: twiddles = {   1, E^ 9, E^12, E^1b, E^24, E^2d, E^36, E^3f, E^48, E^51, E^5a, E^63, E^6c, E^75, E^7e, E^87}
		Block  a: twiddles = {   1, E^ a, E^14, E^1e, E^28, E^32, E^3c, E^46, E^50, E^5a, E^64, E^6e, E^78, E^82, E^8c, E^96}
		Block  b: twiddles = {   1, E^ b, E^16, E^21, E^2c, E^37, E^42, E^4d, E^58, E^63, E^6e, E^79, E^84, E^8f, E^9a, E^a5}
		Block  c: twiddles = {   1, E^ c, E^18, E^24, E^30, E^3c, E^48, E^54, E^60, E^6c, E^78, E^84, E^90, E^9c, E^a8, E^b4}
		Block  d: twiddles = {   1, E^ d, E^1a, E^27, E^34, E^41, E^4e, E^5b, E^68, E^75, E^82, E^8f, E^9c, E^a9, E^b6, E^c3}
		Block  e: twiddles = {   1, E^ e, E^1c, E^2a, E^38, E^46, E^54, E^62, E^70, E^7e, E^8c, E^9a, E^a8, E^b6, E^c4, E^d2}
		Block  f: twiddles = {   1, E^ f, E^1e, E^2d, E^3c, E^4b, E^5a, E^69, E^78, E^87, E^96, E^a5, E^b4, E^c3, E^d2, E^e1}
		Block 10: twiddles = {   1, E^10, E^20, E^30, E^40, E^50, E^60, E^70, E^80, E^90, E^a0, E^b0, E^c0, E^d0, E^e0, E^f0}
		Block 11: twiddles = {   1, E^11, E^22, E^33, E^44, E^55, E^66, E^77, E^88, E^99, E^aa, E^bb, E^cc, E^dd, E^ee, E^ff}
		Block 12: twiddles = {   1, E^12, E^24, E^36, E^48, E^5a, E^6c, E^7e, E^90, E^a2, E^b4, E^c6, E^d8, E^ea, E^fc,-E^ e}
		Block 13: twiddles = {   1, E^13, E^26, E^39, E^4c, E^5f, E^72, E^85, E^98, E^ab, E^be, E^d1, E^e4, E^f7,-E^ a,-E^1d}
		Block 14: twiddles = {   1, E^14, E^28, E^3c, E^50, E^64, E^78, E^8c, E^a0, E^b4, E^c8, E^dc, E^f0,-E^ 4,-E^18,-E^2c}
		Block 15: twiddles = {   1, E^15, E^2a, E^3f, E^54, E^69, E^7e, E^93, E^a8, E^bd, E^d2, E^e7, E^fc,-E^11,-E^26,-E^3b}
		Block 16: twiddles = {   1, E^16, E^2c, E^42, E^58, E^6e, E^84, E^9a, E^b0, E^c6, E^dc, E^f2,-E^ 8,-E^1e,-E^34,-E^4a}
		Block 17: twiddles = {   1, E^17, E^2e, E^45, E^5c, E^73, E^8a, E^a1, E^b8, E^cf, E^e6, E^fd,-E^14,-E^2b,-E^42,-E^59}
		Block 18: twiddles = {   1, E^18, E^30, E^48, E^60, E^78, E^90, E^a8, E^c0, E^d8, E^f0,-E^ 8,-E^20,-E^38,-E^50,-E^68}
		Block 19: twiddles = {   1, E^19, E^32, E^4b, E^64, E^7d, E^96, E^af, E^c8, E^e1, E^fa,-E^13,-E^2c,-E^45,-E^5e,-E^77}
		Block 1a: twiddles = {   1, E^1a, E^34, E^4e, E^68, E^82, E^9c, E^b6, E^d0, E^ea,-E^ 4,-E^1e,-E^38,-E^52,-E^6c,-E^86}
		Block 1b: twiddles = {   1, E^1b, E^36, E^51, E^6c, E^87, E^a2, E^bd, E^d8, E^f3,-E^ e,-E^29,-E^44,-E^5f,-E^7a,-E^95}
		Block 1c: twiddles = {   1, E^1c, E^38, E^54, E^70, E^8c, E^a8, E^c4, E^e0, E^fc,-E^18,-E^34,-E^50,-E^6c,-E^88,-E^a4}
		Block 1d: twiddles = {   1, E^1d, E^3a, E^57, E^74, E^91, E^ae, E^cb, E^e8,-E^ 5,-E^22,-E^3f,-E^5c,-E^79,-E^96,-E^b3}
		Block 1e: twiddles = {   1, E^1e, E^3c, E^5a, E^78, E^96, E^b4, E^d2, E^f0,-E^ e,-E^2c,-E^4a,-E^68,-E^86,-E^a4,-E^c2}
		Block 1f: twiddles = {   1, E^1f, E^3e, E^5d, E^7c, E^9b, E^ba, E^d9, E^f8,-E^17,-E^36,-E^55,-E^74,-E^93,-E^b2,-E^d1} .

	Highest power in lower right slot: 31*15 = 0x1d1 = -0xd1, checks.

	Now further reduce remaining powers:
		~ denotes complex conjugation, * denotes interchange of real and imaginary part, - denotes negation
	(i.e. ~E^n := ( Re(E),-Im(E)), *E^n := ( Im(E), Re(E)), -E^n := (-Re(E),-Im(E)), and any sequence of these
	operators is evaluated right-to-left, e.g. -~*E^n = -~(*E^n) = -~( Im(E), Re(E)) = -( Im(E),-Re(E)) = (-Im(E),+Re(E));
	note that the - and ~ operators commute, as do the - and *, but ~ and * anticommute, i.e. ~*E^n = -*~E^n) ,
	and the values of the individual exponentials are, in terms of the sincos parameters defined in this module:

		E^40 = exp(i* 1*twopi/8) = isrt2*( 1    , 1    )
		E^80 = exp(i* 2*twopi/8) =       ( 0    , 1    ) = I
		E^c0 = exp(i* 3*twopi/8) = isrt2*(-1    , 1    ) = *~E^40

		E^{41,42, ... ,7e,7f} =  *E^{3f,3e, ... , 2, 1}

		E^{81,82, ... ,be,bf} = I.E^{ 1, 2, ... ,3e,3f} = *~E^{ 1, 2, ... ,3e,3f}	(I.{} denotes complex multiply by I)

		E^{c1,c2, ... ,fe,ff} = -~E^{3f,3e, ... , 2, 1} .

	Exploiting these symmetries allows our 32x16 twiddles matrix to be expressed in terms of the powers E^1-3f, 40,
	the imaginary constant I and isrt2 (via E^40 = isrt2*(1+I)) as follows - Note the even rows (counting the top row
	as "row 0") are identical to the twiddles-sets of the radix-256 DFT:

		    0   1     2     3     4      5      6      7      8      9      a      b      c      d      e      f
	 0	{   1,    1,    1,    1,    1,     1,     1,     1,     1,     1,     1,     1,     1,     1,     1,     1}
	 1	{   1, E^ 1, E^ 2, E^ 3, E^ 4,  E^ 5,  E^ 6,  E^ 7,  E^ 8,  E^ 9,  E^ a,  E^ b,  E^ c,  E^ d,  E^ e,  E^ f}
	 2	{   1, E^ 2, E^ 4, E^ 6, E^ 8,  E^ a,  E^ c,  E^ e,  E^10,  E^12,  E^14,  E^16,  E^18,  E^1a,  E^1c,  E^1e}
	 3	{   1, E^ 3, E^ 6, E^ 9, E^ c,  E^ f,  E^12,  E^15,  E^18,  E^1b,  E^1e,  E^21,  E^24,  E^27,  E^2a,  E^2d}
	 4	{   1, E^ 4, E^ 8, E^ c, E^10,  E^14,  E^18,  E^1c,  E^20,  E^24,  E^28,  E^2c,  E^30,  E^34,  E^38,  E^3c}
	 5	{   1, E^ 5, E^ a, E^ f, E^14,  E^19,  E^1e,  E^23,  E^28,  E^2d,  E^32,  E^37,  E^3c, *E^3f, *E^3a, *E^35}
	 6	{   1, E^ 6, E^ c, E^12, E^18,  E^1e,  E^24,  E^2a,  E^30,  E^36,  E^3c, *E^3e, *E^38, *E^32, *E^2c, *E^26}
	 7	{   1, E^ 7, E^ e, E^15, E^1c,  E^23,  E^2a,  E^31,  E^38,  E^3f, *E^3a, *E^33, *E^2c, *E^25, *E^1e, *E^17}
	 8	{   1, E^ 8, E^10, E^18, E^20,  E^28,  E^30,  E^38,  E^40, *E^38, *E^30, *E^28, *E^20, *E^18, *E^10, *E^ 8}
	 9	{   1, E^ 9, E^12, E^1b, E^24,  E^2d,  E^36,  E^3f, *E^38, *E^2f, *E^26, *E^1d, *E^14, *E^ b, *E^ 2,*~E^ 7}
	 a	{   1, E^ a, E^14, E^1e, E^28,  E^32,  E^3c, *E^3a, *E^30, *E^26, *E^1c, *E^12, *E^ 8,*~E^ 2,*~E^ c,*~E^16}
	 b	{   1, E^ b, E^16, E^21, E^2c,  E^37, *E^3e, *E^33, *E^28, *E^1d, *E^12, *E^ 7,*~E^ 4,*~E^ f,*~E^1a,*~E^25}
	 c	{   1, E^ c, E^18, E^24, E^30,  E^3c, *E^38, *E^2c, *E^20, *E^14, *E^ 8,*~E^ 4,*~E^10,*~E^1c,*~E^28,*~E^34}
	 d	{   1, E^ d, E^1a, E^27, E^34, *E^3f, *E^32, *E^25, *E^18, *E^ b,*~E^ 2,*~E^ f,*~E^1c,*~E^29,*~E^36,-~E^3d}
	 e	{   1, E^ e, E^1c, E^2a, E^38, *E^3a, *E^2c, *E^1e, *E^10, *E^ 2,*~E^ c,*~E^1a,*~E^28,*~E^36,-~E^3c,-~E^2e}
	 f	{   1, E^ f, E^1e, E^2d, E^3c, *E^35, *E^26, *E^17, *E^ 8,*~E^ 7,*~E^16,*~E^25,*~E^34,-~E^3d,-~E^2e,-~E^1f}
	10	{   1, E^10, E^20, E^30, E^40, *E^30, *E^20, *E^10,  I.{},*~E^10,*~E^20,*~E^30,*~E^40,-~E^30,-~E^20,-~E^10}
	11	{   1, E^11, E^22, E^33,*E^3c, *E^2b, *E^1a, *E^ 9,*~E^ 8,*~E^19,*~E^2a,*~E^3b,-~E^34,-~E^23,-~E^12,-~E^ 1}
	12	{   1, E^12, E^24, E^36,*E^38, *E^26, *E^14, *E^ 2,*~E^10,*~E^22,*~E^34,-~E^3a,-~E^28,-~E^16,-~E^ 4, -E^ e}
	13	{   1, E^13, E^26, E^39,*E^34, *E^21, *E^ e,*~E^ 5,*~E^18,*~E^2b,*~E^3e,-~E^2f,-~E^1c,-~E^ 9, -E^ a, -E^1d}
	14	{   1, E^14, E^28, E^3c,*E^30, *E^1c, *E^ 8,*~E^ c,*~E^20,*~E^34,-~E^38,-~E^24,-~E^10, -E^ 4, -E^18, -E^2c}
	15	{   1, E^15, E^2a, E^3f,*E^2c, *E^17, *E^ 2,*~E^13,*~E^28,*~E^3d,-~E^2e,-~E^19,-~E^ 4, -E^11, -E^26, -E^3b}
	16	{   1, E^16, E^2c,*E^3e,*E^28, *E^12,*~E^ 4,*~E^1a,*~E^30,-~E^3a,-~E^24,-~E^ e, -E^ 8, -E^1e, -E^34,-*E^36}
	17	{   1, E^17, E^2e,*E^3b,*E^24, *E^ d,*~E^ a,*~E^21,*~E^38,-~E^31,-~E^1a,-~E^ 3, -E^14, -E^2b,-*E^3e,-*E^27}
	18	{   1, E^18, E^30,*E^38,*E^20, *E^ 8,*~E^10,*~E^28,*~E^40,-~E^28,-~E^10, -E^ 8, -E^20, -E^38,-*E^30,-*E^18}
	19	{   1, E^19, E^32,*E^35,*E^1c, *E^ 3,*~E^16,*~E^2f,-~E^38,-~E^1f,-~E^ 6, -E^13, -E^2c,-*E^3b,-*E^22,-*E^ 9}
	1a	{   1, E^1a, E^34,*E^32,*E^18,*~E^ 2,*~E^1c,*~E^36,-~E^30,-~E^16, -E^ 4, -E^1e, -E^38,-*E^2e,-*E^14,~*E^ 6}
	1b	{   1, E^1b, E^36,*E^2f,*E^14,*~E^ 7,*~E^22,*~E^3d,-~E^28,-~E^ d, -E^ e, -E^29,-*E^3c,-*E^21,-*E^ 6,~*E^15}
	1c	{   1, E^1c, E^38,*E^2c,*E^10,*~E^ c,*~E^28,-~E^3c,-~E^20,-~E^ 4, -E^18, -E^34,-*E^30,-*E^14,~*E^ 8,~*E^24}
	1d	{   1, E^1d, E^3a,*E^29,*E^ c,*~E^11,*~E^2e,-~E^35,-~E^18, -E^ 5, -E^22, -E^3f,-*E^24,-*E^ 7,~*E^16,~*E^33}
	1e	{   1, E^1e, E^3c,*E^26,*E^ 8,*~E^16,*~E^34,-~E^2e,-~E^10, -E^ e, -E^2c,-*E^36,-*E^18,~*E^ 6,~*E^24, ~E^3e}
	1f	{   1, E^1f, E^3e,*E^23,*E^ 4,*~E^1b,*~E^3a,-~E^27,-~E^ 8, -E^17, -E^36,-*E^2b,-*E^ c,~*E^13,~*E^32, ~E^2f} .

	Or with both rows rearranged in BR32 {0+8+4+c+2+a+6+e+1+9+5+d+3+b+7+f+} (where + denotes 'preceding index + 0x10') and
	columns in BR16 {084c2a6e195d3b7f} (pairwise-swap cols 1/8,2/4,3/c,5/a,7/e,b/d) order:

		    0    8     4      c     2      a      6      e     1      9      5      d     3      b      7      f
	 0	{   1,     1,    1,     1,    1,     1,     1,     1,    1,     1,     1,     1,    1,     1,     1,     1}
	10	{   1,  I.{}, E^40,*~E^40, E^20,*~E^20, *E^20,-~E^20, E^10,*~E^10, *E^30,-~E^30, E^30,*~E^30, *E^10,-~E^10}
	 8	{   1,  E^40, E^20, *E^20, E^10, *E^30,  E^30, *E^10, E^ 8, *E^38,  E^28, *E^18, E^18, *E^28,  E^38, *E^ 8}
	18	{   1,*~E^40,*E^20, -E^20, E^30,-~E^10,*~E^10,-*E^30, E^18,-~E^28, *E^ 8, -E^38,*E^38, -E^ 8,*~E^28,-*E^18}
	 4	{   1,  E^20, E^10,  E^30, E^ 8,  E^28,  E^18,  E^38, E^ 4,  E^24,  E^14,  E^34, E^ c,  E^2c,  E^1c,  E^3c}
	14	{   1,*~E^20,*E^30,-~E^10, E^28,-~E^38, *E^ 8, -E^18, E^14,*~E^34, *E^1c, -E^ 4, E^3c,-~E^24,*~E^ c, -E^2c}
	 c	{   1, *E^20, E^30,*~E^10, E^18, *E^ 8, *E^38,*~E^28, E^ c, *E^14,  E^3c,*~E^1c, E^24,*~E^ 4, *E^2c,*~E^34}
	1c	{   1,-~E^20,*E^10,-*E^30, E^38, -E^18,*~E^28,~*E^ 8, E^1c,-~E^ 4,*~E^ c,-*E^14,*E^2c, -E^34,-~E^3c,~*E^24}
	 2	{   1,  E^10, E^ 8,  E^18, E^ 4,  E^14,  E^ c,  E^1c, E^ 2,  E^12,  E^ a,  E^1a, E^ 6,  E^16,  E^ e,  E^1e}
	12	{   1,*~E^10,*E^38,-~E^28, E^24,*~E^34, *E^14,-~E^ 4, E^12,*~E^22, *E^26,-~E^16, E^36,-~E^3a, *E^ 2, -E^ e}
	 a	{   1, *E^30, E^28, *E^ 8, E^14, *E^1c,  E^3c,*~E^ c, E^ a, *E^26,  E^32,*~E^ 2, E^1e, *E^12, *E^3a,*~E^16}
	1a	{   1,-~E^30,*E^18, -E^38, E^34, -E^ 4,*~E^1c,-*E^14, E^1a,-~E^16,*~E^ 2,-*E^2e,*E^32, -E^1e,*~E^36,~*E^ 6}
	 6	{   1,  E^30, E^18, *E^38, E^ c,  E^3c,  E^24, *E^2c, E^ 6,  E^36,  E^1e, *E^32, E^12, *E^3e,  E^2a, *E^26}
	16	{   1,*~E^30,*E^28, -E^ 8, E^2c,-~E^24,*~E^ 4, -E^34, E^16,-~E^3a, *E^12, -E^1e,*E^3e,-~E^ e,*~E^1a,-*E^36}
	 e	{   1, *E^10, E^38,*~E^28, E^1c,*~E^ c, *E^2c,-~E^3c, E^ e, *E^ 2, *E^3a,*~E^36, E^2a,*~E^1a, *E^1e,-~E^2e}
	1e	{   1,-~E^10,*E^ 8,-*E^18, E^3c, -E^2c,*~E^34,~*E^24, E^1e, -E^ e,*~E^16,~*E^ 6,*E^26,-*E^36,-~E^2e, ~E^3e}
	 1	{   1,  E^ 8, E^ 4,  E^ c, E^ 2,  E^ a,  E^ 6,  E^ e, E^ 1,  E^ 9,  E^ 5,  E^ d, E^ 3,  E^ b,  E^ 7,  E^ f}
	11	{   1,*~E^ 8,*E^3c,-~E^34, E^22,*~E^2a, *E^1a,-~E^12, E^11,*~E^19, *E^2b,-~E^23, E^33,*~E^3b, *E^ 9,-~E^ 1}
	 9	{   1, *E^38, E^24, *E^14, E^12, *E^26,  E^36, *E^ 2, E^ 9, *E^2f,  E^2d, *E^ b, E^1b, *E^1d,  E^3f,*~E^ 7}
	19	{   1,-~E^38,*E^1c, -E^2c, E^32,-~E^ 6,*~E^16,-*E^22, E^19,-~E^1f, *E^ 3,-*E^3b,*E^35, -E^13,*~E^2f,-*E^ 9}
	 5	{   1,  E^28, E^14,  E^3c, E^ a,  E^32,  E^1e, *E^3a, E^ 5,  E^2d,  E^19, *E^3f, E^ f,  E^37,  E^23, *E^35}
	15	{   1,*~E^28,*E^2c,-~E^ 4, E^2a,-~E^2e, *E^ 2, -E^26, E^15,*~E^3d, *E^17, -E^11, E^3f,-~E^19,*~E^13, -E^3b}
	 d	{   1, *E^18, E^34,*~E^1c, E^1a,*~E^ 2, *E^32,*~E^36, E^ d, *E^ b, *E^3f,*~E^29, E^27,*~E^ f, *E^25,-~E^3d}
	1d	{   1,-~E^18,*E^ c,-*E^24, E^3a, -E^22,*~E^2e,~*E^16, E^1d, -E^ 5,*~E^11,-*E^ 7,*E^29, -E^3f,-~E^35,~*E^33}
	 3	{   1,  E^18, E^ c,  E^24, E^ 6,  E^1e,  E^12,  E^2a, E^ 3,  E^1b,  E^ f,  E^27, E^ 9,  E^21,  E^15,  E^2d}
	13	{   1,*~E^18,*E^34,-~E^1c, E^26,*~E^3e, *E^ e, -E^ a, E^13,*~E^2b, *E^21,-~E^ 9, E^39,-~E^2f,*~E^ 5, -E^1d}
	 b	{   1, *E^28, E^2c,*~E^ 4, E^16, *E^12, *E^3e,*~E^1a, E^ b, *E^1d,  E^37,*~E^ f, E^21, *E^ 7, *E^33,*~E^25}
	1b	{   1,-~E^28,*E^14,-*E^3c, E^36, -E^ e,*~E^22,-*E^ 6, E^1b,-~E^ d,*~E^ 7,-*E^21,*E^2f, -E^29,*~E^3d,~*E^15}
	 7	{   1,  E^38, E^1c, *E^2c, E^ e, *E^3a,  E^2a, *E^1e, E^ 7,  E^3f,  E^23, *E^25, E^15, *E^33,  E^31, *E^17}
	17	{   1,*~E^38,*E^24, -E^14, E^2e,-~E^1a,*~E^ a,-*E^3e, E^17,-~E^31, *E^ d, -E^2b,*E^3b,-~E^ 3,*~E^21,-*E^27}
	 f	{   1, *E^ 8, E^3c,*~E^34, E^1e,*~E^16, *E^26,-~E^2e, E^ f,*~E^ 7, *E^35,-~E^3d, E^2d,*~E^25, *E^17,-~E^1f}
	1f	{   1,-~E^ 8,*E^ 4,-*E^ c, E^3e, -E^36,*~E^3a,~*E^32, E^1f, -E^17,*~E^1b,~*E^13,*E^23,-*E^2b,-~E^27, ~E^2f} .

	Only the last 15 inputs to each of the radix-16 transforms ("Blocks") 1 through 15 are multiplied by non-unity twiddles.
	For DIF we process both the blocks, and the twiddles within each block, in bit-reversed order.
	One can see from the data below that aside from using a twiddleless DIF for Block 0 there is
	little to be gained from trying to exploit other "special form" twiddles such as I and isrt2*[+-1,+-1].
	Thus our radix-16-DIF-with-twiddles macro uses generic complex MUL for the 15 non-unity twiddles of each invocation.
	*/
/**** EXPRESS POWERS WHICH ARE MULTS OF 2,4,8,16,32,64 RESP. IN TERMS OF D,C,B,A := RADIX-256,128,64,32 FUNDAMENTAL ROOTS ****/
// First 7 (of 15 nontrivial) of each twiddles set must match those of the analogous radix-16 DFT call in the radix256 DIF.

	// Block 0: has all-unity twiddles
	tptr = t;
	jt = j1;	jp = j2;
	// Twiddleless DIF bit-reverses its outputs, so a_p* terms appear in BR-order [swap index pairs 1/8,2/4,3/c,5/a,7/e.b/d]:
	RADIX_16_DIF(
		tptr->re,tptr->im,(tptr+0x100)->re,(tptr+0x100)->im,(tptr+0x80)->re,(tptr+0x80)->im,(tptr+0x180)->re,(tptr+0x180)->im,(tptr+0x40)->re,(tptr+0x40)->im,(tptr+0x140)->re,(tptr+0x140)->im,(tptr+0xc0)->re,(tptr+0xc0)->im,(tptr+0x1c0)->re,(tptr+0x1c0)->im,(tptr+0x20)->re,(tptr+0x20)->im,(tptr+0x120)->re,(tptr+0x120)->im,(tptr+0xa0)->re,(tptr+0xa0)->im,(tptr+0x1a0)->re,(tptr+0x1a0)->im,(tptr+0x60)->re,(tptr+0x60)->im,(tptr+0x160)->re,(tptr+0x160)->im,(tptr+0xe0)->re,(tptr+0xe0)->im,(tptr+0x1e0)->re,(tptr+0x1e0)->im,
		a[jt],a[jp],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p08],a[jp+p08],a[jt+p09],a[jp+p09],a[jt+p0a],a[jp+p0a],a[jt+p0b],a[jp+p0b],a[jt+p0c],a[jp+p0c],a[jt+p0d],a[jp+p0d],a[jt+p0e],a[jp+p0e],a[jt+p0f],a[jp+p0f],
		c16,s16
	);	tptr++;

	// Remaining 63 sets of macro calls done in loop:
	for(i = 1; i < 32; i++) {
		jt = j1 + i_offsets[i]; jp = j2 + i_offsets[i];	// poffs[] = p10,p20,...,p1f0
		addr = DFT1024_TWIDDLES[i]; addi = addr+1;	// Pointer to required row of 2-D twiddles array
		RADIX_16_DIF_TWIDDLE_OOP(
			tptr->re,tptr->im,(tptr+0x100)->re,(tptr+0x100)->im,(tptr+0x80)->re,(tptr+0x80)->im,(tptr+0x180)->re,(tptr+0x180)->im,(tptr+0x40)->re,(tptr+0x40)->im,(tptr+0x140)->re,(tptr+0x140)->im,(tptr+0xc0)->re,(tptr+0xc0)->im,(tptr+0x1c0)->re,(tptr+0x1c0)->im,(tptr+0x20)->re,(tptr+0x20)->im,(tptr+0x120)->re,(tptr+0x120)->im,(tptr+0xa0)->re,(tptr+0xa0)->im,(tptr+0x1a0)->re,(tptr+0x1a0)->im,(tptr+0x60)->re,(tptr+0x60)->im,(tptr+0x160)->re,(tptr+0x160)->im,(tptr+0xe0)->re,(tptr+0xe0)->im,(tptr+0x1e0)->re,(tptr+0x1e0)->im,
			a[jt],a[jp],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p08],a[jp+p08],a[jt+p09],a[jp+p09],a[jt+p0a],a[jp+p0a],a[jt+p0b],a[jp+p0b],a[jt+p0c],a[jp+p0c],a[jt+p0d],a[jp+p0d],a[jt+p0e],a[jp+p0e],a[jt+p0f],a[jp+p0f],
			*(addr+0x00),*(addi+0x00), *(addr+0x02),*(addi+0x02), *(addr+0x04),*(addi+0x04), *(addr+0x06),*(addi+0x06), *(addr+0x08),*(addi+0x08), *(addr+0x0a),*(addi+0x0a), *(addr+0x0c),*(addi+0x0c), *(addr+0x0e),*(addi+0x0e), *(addr+0x10),*(addi+0x10), *(addr+0x12),*(addi+0x12), *(addr+0x14),*(addi+0x14), *(addr+0x16),*(addi+0x16), *(addr+0x18),*(addi+0x18), *(addr+0x1a),*(addi+0x1a), *(addr+0x1c),*(addi+0x1c),
			c16,s16
		);	tptr++;
	}

	}
}

/**************/

void radix512_dit_pass1(double a[], int n)
{
/*
!...Acronym: DIT = Decimation In Time
!
!...Subroutine to perform an initial radix-256 complex inverse DIT FFT pass on the data in the length-N real vector A.
!
!   See the documentation in radix256_0dif_pass for further details on storage, indexing and DFT decomposition.
*/
	int i,j,j1,j2,jt,jp,ju,jv;
	static int NDIVR,first_entry=TRUE,
			p01,p02,p03,p04,p05,p06,p07,p08,p09,p0a,p0b,p0c,p0d,p0e,p0f,
			p10,p20,p30,p40,p50,p60,p70,p80,p90,pa0,pb0,pc0,pd0,pe0,pf0,p100;
	static int i_offsets[32], o_offsets[32];
	static int poffs[16],po_br[32];
	// We prefer pointer-based array-element access, because that allows our radix16 DFT-with-twiddles
	// to look the same in terms of array-element arglists:
	double *addr,*addi;
	struct complex *tptr;
	#include "radix1024_twiddles.h"
	// Local storage: We must use an array here because scalars have no guarantees about relative address offsets
	// [and even if those are contiguous-as-hoped-for, they may run in reverse]; Make array type (struct complex)
	// to allow us to use the same offset-indexing as in the original radix-32 in-place DFT macros:
	struct complex t[RADIX];

	if(!first_entry && (n >> 9) != NDIVR)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		first_entry=FALSE;
		NDIVR = n >> 9;

		p01 = NDIVR;
		p02 = p01 + NDIVR;		p01 = p01 + ( (p01 >> DAT_BITS) << PAD_BITS );
		p03 = p02 + NDIVR;		p02 = p02 + ( (p02 >> DAT_BITS) << PAD_BITS );
		p04 = p03 + NDIVR;		p03 = p03 + ( (p03 >> DAT_BITS) << PAD_BITS );
		p05 = p04 + NDIVR;		p04 = p04 + ( (p04 >> DAT_BITS) << PAD_BITS );
		p06 = p05 + NDIVR;		p05 = p05 + ( (p05 >> DAT_BITS) << PAD_BITS );
		p07 = p06 + NDIVR;		p06 = p06 + ( (p06 >> DAT_BITS) << PAD_BITS );
		p08 = p07 + NDIVR;		p07 = p07 + ( (p07 >> DAT_BITS) << PAD_BITS );
		p09 = p08 + NDIVR;		p08 = p08 + ( (p08 >> DAT_BITS) << PAD_BITS );
		p0a = p09 + NDIVR;		p09 = p09 + ( (p09 >> DAT_BITS) << PAD_BITS );
		p0b = p0a + NDIVR;		p0a = p0a + ( (p0a >> DAT_BITS) << PAD_BITS );
		p0c = p0b + NDIVR;		p0b = p0b + ( (p0b >> DAT_BITS) << PAD_BITS );
		p0d = p0c + NDIVR;		p0c = p0c + ( (p0c >> DAT_BITS) << PAD_BITS );
		p0e = p0d + NDIVR;		p0d = p0d + ( (p0d >> DAT_BITS) << PAD_BITS );
		p0f = p0e + NDIVR;		p0e = p0e + ( (p0e >> DAT_BITS) << PAD_BITS );
		NDIVR <<= 4;			p0f = p0f + ( (p0f >> DAT_BITS) << PAD_BITS );
		p10 = NDIVR;
		p20 = p10 + NDIVR;		p10 = p10 + ( (p10 >> DAT_BITS) << PAD_BITS );
		p30 = p20 + NDIVR;		p20 = p20 + ( (p20 >> DAT_BITS) << PAD_BITS );
		p40 = p30 + NDIVR;		p30 = p30 + ( (p30 >> DAT_BITS) << PAD_BITS );
		p50 = p40 + NDIVR;		p40 = p40 + ( (p40 >> DAT_BITS) << PAD_BITS );
		p60 = p50 + NDIVR;		p50 = p50 + ( (p50 >> DAT_BITS) << PAD_BITS );
		p70 = p60 + NDIVR;		p60 = p60 + ( (p60 >> DAT_BITS) << PAD_BITS );
		p80 = p70 + NDIVR;		p70 = p70 + ( (p70 >> DAT_BITS) << PAD_BITS );
		p90 = p80 + NDIVR;		p80 = p80 + ( (p80 >> DAT_BITS) << PAD_BITS );
		pa0 = p90 + NDIVR;		p90 = p90 + ( (p90 >> DAT_BITS) << PAD_BITS );
		pb0 = pa0 + NDIVR;		pa0 = pa0 + ( (pa0 >> DAT_BITS) << PAD_BITS );
		pc0 = pb0 + NDIVR;		pb0 = pb0 + ( (pb0 >> DAT_BITS) << PAD_BITS );
		pd0 = pc0 + NDIVR;		pc0 = pc0 + ( (pc0 >> DAT_BITS) << PAD_BITS );
		pe0 = pd0 + NDIVR;		pd0 = pd0 + ( (pd0 >> DAT_BITS) << PAD_BITS );
		pf0 = pe0 + NDIVR;		pe0 = pe0 + ( (pe0 >> DAT_BITS) << PAD_BITS );
		p100= pf0 + NDIVR;		pf0 = pf0 + ( (pf0 >> DAT_BITS) << PAD_BITS );
		NDIVR >>= 4;			p100= p100+ ( (p100>> DAT_BITS) << PAD_BITS );

		// Set array offsets for radix-32 DFT in/outputs:
		i_offsets[0x00] = 0  ;		i_offsets[0x10] = 0  +p10;
		i_offsets[0x01] = p01;		i_offsets[0x11] = p01+p10;
		i_offsets[0x02] = p02;		i_offsets[0x12] = p02+p10;
		i_offsets[0x03] = p03;		i_offsets[0x13] = p03+p10;
		i_offsets[0x04] = p04;		i_offsets[0x14] = p04+p10;
		i_offsets[0x05] = p05;		i_offsets[0x15] = p05+p10;
		i_offsets[0x06] = p06;		i_offsets[0x16] = p06+p10;
		i_offsets[0x07] = p07;		i_offsets[0x17] = p07+p10;
		i_offsets[0x08] = p08;		i_offsets[0x18] = p08+p10;
		i_offsets[0x09] = p09;		i_offsets[0x19] = p09+p10;
		i_offsets[0x0a] = p0a;		i_offsets[0x1a] = p0a+p10;
		i_offsets[0x0b] = p0b;		i_offsets[0x1b] = p0b+p10;
		i_offsets[0x0c] = p0c;		i_offsets[0x1c] = p0c+p10;
		i_offsets[0x0d] = p0d;		i_offsets[0x1d] = p0d+p10;
		i_offsets[0x0e] = p0e;		i_offsets[0x1e] = p0e+p10;
		i_offsets[0x0f] = p0f;		i_offsets[0x1f] = p0f+p10;

		o_offsets[0x00] = 0x00<<1;	o_offsets[0x10] = 0x10<<1;
		o_offsets[0x01] = 0x01<<1;	o_offsets[0x11] = 0x11<<1;
		o_offsets[0x02] = 0x02<<1;	o_offsets[0x12] = 0x12<<1;
		o_offsets[0x03] = 0x03<<1;	o_offsets[0x13] = 0x13<<1;
		o_offsets[0x04] = 0x04<<1;	o_offsets[0x14] = 0x14<<1;
		o_offsets[0x05] = 0x05<<1;	o_offsets[0x15] = 0x15<<1;
		o_offsets[0x06] = 0x06<<1;	o_offsets[0x16] = 0x16<<1;
		o_offsets[0x07] = 0x07<<1;	o_offsets[0x17] = 0x17<<1;
		o_offsets[0x08] = 0x08<<1;	o_offsets[0x18] = 0x18<<1;
		o_offsets[0x09] = 0x09<<1;	o_offsets[0x19] = 0x19<<1;
		o_offsets[0x0a] = 0x0a<<1;	o_offsets[0x1a] = 0x1a<<1;
		o_offsets[0x0b] = 0x0b<<1;	o_offsets[0x1b] = 0x1b<<1;
		o_offsets[0x0c] = 0x0c<<1;	o_offsets[0x1c] = 0x1c<<1;
		o_offsets[0x0d] = 0x0d<<1;	o_offsets[0x1d] = 0x1d<<1;
		o_offsets[0x0e] = 0x0e<<1;	o_offsets[0x1e] = 0x1e<<1;
		o_offsets[0x0f] = 0x0f<<1;	o_offsets[0x1f] = 0x1f<<1;

		po_br[2*0x0] =   0; po_br[2*0x1] = p08; po_br[2*0x2] = p04; po_br[2*0x3] = p0c; po_br[2*0x4] = p02; po_br[2*0x5] = p0a; po_br[2*0x6] = p06; po_br[2*0x7] = p0e; po_br[2*0x8] = p01; po_br[2*0x9] = p09; po_br[2*0xa] = p05; po_br[2*0xb] = p0d; po_br[2*0xc] = p03; po_br[2*0xd] = p0b; po_br[2*0xe] = p07; po_br[2*0xf] = p0f;
		// Each of the foregoing 16 indices is head of a (i0,i0+p10) pair:
		for(i = 0; i < 32; i += 2) {
			po_br[i+1] = po_br[i] + p10;
		}
		poffs[0x0] = 0; poffs[0x1] = p20; poffs[0x2] = p40; poffs[0x3] = p60; poffs[0x4] = p80; poffs[0x5] = pa0; poffs[0x6] = pc0; poffs[0x7] = pe0; poffs[0x8] = p100; poffs[0x9] = p20 + p100; poffs[0xa] = p40 + p100; poffs[0xb] = p60 + p100; poffs[0xc] = p80 + p100; poffs[0xd] = pa0 + p100; poffs[0xe] = pc0 + p100; poffs[0xf] = pe0 + p100;
	}

/*...The radix-512 pass is here.	*/

	for(j = 0; j < NDIVR; j += 2)
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

	// Gather the needed data and do 16 twiddleless length-32 subtransforms, with p-offsets in-order:

		for(i = 0, jp = 0; i < 16; i++, jp += 32) {
			jt = j1 + poffs[i];	// poffs[] = p20,p40,...,p1e0
			RADIX_32_DIT((a+jt),i_offsets,RE_IM_STRIDE, (double *)(t+jp),o_offsets,1);
		}

	/*...and now do 32 radix-16 subtransforms, including the internal twiddle factors - we use the same positive-power
	roots as in the DIF here, just fiddle with signs within the macro to effect the conjugate-multiplies. Twiddles occur
	in the same order here as DIF, but the in-and-output-index offsets are BRed: j1 + p[0,8,4,c,2,a,6,e,1,9,5,d,3,b,7,f].
	*/
		// Block 0: has all-unity twiddles
		tptr = t;
		jt = j1;	jp = j2;
		ju = jt+p100;	jv = jp+p100;
		// Twiddleless DIF bit-reverses its outputs, so a_p* terms appear in BR-order [swap index pairs 1/8,2/4,3/c,5/a,7/e.b/d]:
		RADIX_16_DIT(
			tptr->re,tptr->im,(tptr+0x20)->re,(tptr+0x20)->im,(tptr+0x40)->re,(tptr+0x40)->im,(tptr+0x60)->re,(tptr+0x60)->im,(tptr+0x80)->re,(tptr+0x80)->im,(tptr+0xa0)->re,(tptr+0xa0)->im,(tptr+0xc0)->re,(tptr+0xc0)->im,(tptr+0xe0)->re,(tptr+0xe0)->im,(tptr+0x100)->re,(tptr+0x100)->im,(tptr+0x120)->re,(tptr+0x120)->im,(tptr+0x140)->re,(tptr+0x140)->im,(tptr+0x160)->re,(tptr+0x160)->im,(tptr+0x180)->re,(tptr+0x180)->im,(tptr+0x1a0)->re,(tptr+0x1a0)->im,(tptr+0x1c0)->re,(tptr+0x1c0)->im,(tptr+0x1e0)->re,(tptr+0x1e0)->im,
			a[jt    ],a[jp    ],a[jt+p20],a[jp+p20],a[jt+p40],a[jp+p40],a[jt+p60],a[jp+p60],a[jt+p80],a[jp+p80],a[jt+pa0],a[jp+pa0],a[jt+pc0],a[jp+pc0],a[jt+pe0],a[jp+pe0],a[ju],a[jv],a[ju+p20],a[jv+p20],a[ju+p40],a[jv+p40],a[ju+p60],a[jv+p60],a[ju+p80],a[jv+p80],a[ju+pa0],a[jv+pa0],a[ju+pc0],a[jv+pc0],a[ju+pe0],a[jv+pe0],
			c16,s16
		);

		// Remaining 31 sets of macro calls done in loop:
		for(i = 1; i < 32; i++) {
			tptr = t + reverse(i,5);
			jt = j1 + po_br[i]; jp = j2 + po_br[i];
			ju = jt+p100;	jv = jp+p100;
			addr = DFT1024_TWIDDLES[i]; addi = addr+1;	// Pointer to required row of 2-D twiddles array
			RADIX_16_DIT_TWIDDLE_OOP(
				tptr->re,tptr->im,(tptr+0x20)->re,(tptr+0x20)->im,(tptr+0x40)->re,(tptr+0x40)->im,(tptr+0x60)->re,(tptr+0x60)->im,(tptr+0x80)->re,(tptr+0x80)->im,(tptr+0xa0)->re,(tptr+0xa0)->im,(tptr+0xc0)->re,(tptr+0xc0)->im,(tptr+0xe0)->re,(tptr+0xe0)->im,(tptr+0x100)->re,(tptr+0x100)->im,(tptr+0x120)->re,(tptr+0x120)->im,(tptr+0x140)->re,(tptr+0x140)->im,(tptr+0x160)->re,(tptr+0x160)->im,(tptr+0x180)->re,(tptr+0x180)->im,(tptr+0x1a0)->re,(tptr+0x1a0)->im,(tptr+0x1c0)->re,(tptr+0x1c0)->im,(tptr+0x1e0)->re,(tptr+0x1e0)->im,
				a[jt    ],a[jp    ],a[jt+p20],a[jp+p20],a[jt+p40],a[jp+p40],a[jt+p60],a[jp+p60],a[jt+p80],a[jp+p80],a[jt+pa0],a[jp+pa0],a[jt+pc0],a[jp+pc0],a[jt+pe0],a[jp+pe0],a[ju],a[jv],a[ju+p20],a[jv+p20],a[ju+p40],a[jv+p40],a[ju+p60],a[jv+p60],a[ju+p80],a[jv+p80],a[ju+pa0],a[jv+pa0],a[ju+pc0],a[jv+pc0],a[ju+pe0],a[jv+pe0],
				*(addr+0x00),*(addi+0x00), *(addr+0x02),*(addi+0x02), *(addr+0x04),*(addi+0x04), *(addr+0x06),*(addi+0x06), *(addr+0x08),*(addi+0x08), *(addr+0x0a),*(addi+0x0a), *(addr+0x0c),*(addi+0x0c), *(addr+0x0e),*(addi+0x0e), *(addr+0x10),*(addi+0x10), *(addr+0x12),*(addi+0x12), *(addr+0x14),*(addi+0x14), *(addr+0x16),*(addi+0x16), *(addr+0x18),*(addi+0x18), *(addr+0x1a),*(addi+0x1a), *(addr+0x1c),*(addi+0x1c),
				c16,s16
			);
		}

	}
}

#undef RADIX
