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

/*
	Given an FFT length (in terms of 1K doubles, i.e. there are N = kblocks*2^10 doubles)
	and (optionally) a preferred radix set for the given FFT length (typically based on
	actual timing data for the target platform), returns the set of complex radices to be
	used in performing the FFT in the input rvec[] vector.

	*** NOTE ***: Every legal FFT length, i.e. one with a valid entry in the big table below,
		*must* have a corresponding set of pretested self-test residues in the MersVec table in Mlucas.c!
		[FFT lengths for which Fermat-mod is also supported should similarly have a corr. entry in Fermvec.]

	There are 3 different return values:

		ERR_FFTLENGTH_ILLEGAL	- There is no support for an FFT of length (kblocks * 2^10).
							This is either because kblocks has an odd factor which exceeds
							the largest supported odd radix (currently 15), or because the
							FFT length in question is smaller or larger than the range of
							the lookup table in this function, or some combination of the
							two (e.g. FFT length 2.5K or 7.5K).

		ERR_RADIXSET_UNAVAILABLE - The specified (by way of on radix-set index) radix set
							is not available. In this case the return value of *nradices if the pointer
							is non-null) is set to the number of distinct radix sets available for the
							FFT length in question.

		return 0;			- SUCCESS: both requested FFT length and radix_set supported.
							If the nradices and rvec[] pointers are non-null, the radices
							corresponding to the requested radix_set are written to the latter,
 							and the number of radices in the radix_set is returned via nradices.

	For any supported FFT length (i.e. return value != ERR_FFTLENGTH_ILLEGAL), there must
	be at least one supported FFT radix set, i.e. radix_set = 0 is a valid choice for such.

	The radices are intended be performed in the order given during a forward FFT, and in
	reverse order during an inverse FFT. For the case N = p*2^k with p an odd prime,
	we assume the odd radix p is done first during the forward FFT (possibly as part of
	a single radix-R pass with R = p*(small power of 2).

	Especially for large lengths we can gain a small amount of speed by restricting
	ourselves to radices >= 8 in the power-of-2 part.
*/

int	get_fft_radices(uint32 kblocks, int radix_set, uint32 *nradices, uint32 radix_vec[], int radix_vec_dim)
{
	uint32 rvec[10] = {0,0,0,0,0,0,0,0,0,0};	// Temporary storage for FFT radices
	uint32 i, n, numrad, rad_prod;
	if(!nradices) nradices = &numrad;	// If ptr not supplied, dump datum into local-var to ease code logic

	// 2016: We arrange the allowed radix sets with the larger-leading-radices-first, as SIMD builds typically
	// have restrictions on the power-of-2 component of the leading radix, e.g. avx512 requires this to be
	// at least 8 and AVX/SSE2 require the leading radix to be divisible by 4. These requirements are
	// encoded via #ifdefs within the selection table for each FFT length.
	/*
	[older commentary, subject to modification by the above notes]:
	Case selection based on real vector length - product of complex radices must equal complex vector length = N2.
	Note that the first radix must be one of (5,6,7,8,9,10,11,12,13,14,15,16) and the last must be 16 or 32.
	Intermediate radices must be a combination of 8, 16 and 32. For technical reasons related to the implementation
	of data-block-processing in the Mersenne-mod radix{16|32}_wrapper_square routines,
	THERE MUST AT LEAST ONE INTERMEDIATE FFT RADIX, i.e. there must be 3 or more total radices - that's why we e.g.
	can't do 1K - that would require either a first-radix-4 or last-radix-8 capability, and at this point it's just
	not been worth my time to implement either - for numbers that small, I suggest you use PARI or Mathematica.

	In general we should attempt to order the radices in ascending order, i.e. have the larger radices
	get done at the end of the DIF FFT and beginning of the DIT IFFT, since data accesses become increasingly
	local as one gets closer to the wrapper/square step. In fact, data accesses are contiguous in the DFT/IDFT
	blocks sandwiching the wrapper/square step, so we can afford to use a large-radix (say 32) DFT
	there if needed to minimize the number of radices to be processed. The one common exception to the
	ascending-radix rule is the final pass of the forward FFT (and hence initial pass of the inverse FFT) -
	on some architectures it is better to put the largest radix into the penultimate FFT pass and use a smaller
	one for the wrapper/square step. Running the standard self-test(s) will indicate the best choice here.
	*/

	if(kblocks == 0)
	{
		fprintf(stderr,"ERROR: %d K must be > 0.\n",kblocks);
		return ERR_FFTLENGTH_ILLEGAL;
	}
	n = kblocks << 10;
	if(kblocks != (n>>10))
	{
		fprintf(stderr,"ERROR: length %d K overflows N = 1024*K.\n",kblocks);
		return ERR_FFTLENGTH_ILLEGAL;
	}
	i = kblocks >> trailz32(kblocks);	// odd component of FFT length
	if(i > 15 && i != 31 && i != 63)	// Apr 2014: Add radix-31,63 support for purpose of ROE testing of F33-mod arithmetic
	{
		fprintf(stderr,"ERROR: specified FFT length %d K has non-power-of-2 component %d,\n", kblocks,i);
		fprintf(stderr,"       which exceeds the maximum allowed odd radix of 15.\n");
		return ERR_FFTLENGTH_ILLEGAL;
	}

	/* Real-array length, in CS notation (1K=1024 8-byte elements):	*/
	switch(kblocks)
	{
	case 1 :						/* 1K */
		switch(radix_set) {
		case 0 :
			numrad = 2; rvec[0] = 32; rvec[1] = 16; break;
		case 1 :
			numrad = 2; rvec[0] = 16; rvec[1] = 32; break;
		default :
			*nradices = 2;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 2 :						/* 2K */
		switch(radix_set) {
		case 0 :
			numrad = 3; rvec[0] =  8; rvec[1] =  8; rvec[2] = 16; break;
		case 1 :
			numrad = 2; rvec[0] = 32; rvec[1] = 32; break;
		default :
			*nradices = 2;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 3 :						/* 3K */
		switch(radix_set) {
		case 0 :
			numrad = 3; rvec[0] = 12; rvec[1] =  8; rvec[2] = 16; break;	// Note: Radix-12 has no SSE2-optimized DFT, but must support at least one radix set in vector-build mode
	  #ifdef USE_SSE2
		default :
			*nradices = 1;	return ERR_RADIXSET_UNAVAILABLE;
	  #else
		case 1 :
			numrad = 3; rvec[0] =  6; rvec[1] = 16; rvec[2] = 16; break;
		default :
			*nradices = 2;	return ERR_RADIXSET_UNAVAILABLE;
	  #endif
		}; break;
	case 4 :						/* 4K */
		switch(radix_set) {
		case 0 :
			numrad = 3; rvec[0] = 16; rvec[1] =  8; rvec[2] = 16; break;
		case 1 :
			numrad = 3; rvec[0] =  8; rvec[1] = 16; rvec[2] = 16; break;
	  #ifdef USE_SSE2
		default :
			*nradices = 2;	return ERR_RADIXSET_UNAVAILABLE;
	  #else
		case 2 :
			numrad = 3; rvec[0] =  8; rvec[1] =  8; rvec[2] = 32; break;
		default :
			*nradices = 3;	return ERR_RADIXSET_UNAVAILABLE;
	  #endif
		}; break;
	case 5 :						/* 5K */
		switch(radix_set) {
		case 0 :
			numrad = 3; rvec[0] = 20; rvec[1] =  8; rvec[2] = 16; break;
	  #ifdef USE_SSE2
		default :
			*nradices = 1;	return ERR_RADIXSET_UNAVAILABLE;
	  #else
		case 1 :
			numrad = 3; rvec[0] = 10; rvec[1] = 16; rvec[2] = 16; break;
		case 2 :
			numrad = 3; rvec[0] =  5; rvec[1] = 16; rvec[2] = 32; break;
		default :
			*nradices = 3;	return ERR_RADIXSET_UNAVAILABLE;
	  #endif
		}; break;
	case 6 :						/* 6K */
		switch(radix_set) {
		case 0 :
			numrad = 3; rvec[0] = 24; rvec[1] =  8; rvec[2] = 16; break;
	  #ifdef USE_SSE2
		default :
			*nradices = 1;	return ERR_RADIXSET_UNAVAILABLE;
	  #else
		case 1 :
			numrad = 3; rvec[0] = 12; rvec[1] = 16; rvec[2] = 16; break;
		case 2 :
			numrad = 3; rvec[0] =  6; rvec[1] = 16; rvec[2] = 32; break;
		default :
			*nradices = 3;	return ERR_RADIXSET_UNAVAILABLE;
	  #endif
		}; break;
	case 7 :						/* 7K */
		switch(radix_set) {
		case 0 :
			numrad = 3; rvec[0] = 28; rvec[1] =  8; rvec[2] = 16; break;
	  #ifdef USE_SSE2
		default :
			*nradices = 1;	return ERR_RADIXSET_UNAVAILABLE;
	  #else
		case 1 :
			numrad = 3; rvec[0] = 14; rvec[1] = 16; rvec[2] = 16; break;
		case 2 :
			numrad = 3; rvec[0] =  7; rvec[1] = 16; rvec[2] = 32; break;
		default :
			*nradices = 3;	return ERR_RADIXSET_UNAVAILABLE;
	  #endif
		}; break;
	case 8 :						/* 8K */
		switch(radix_set) {
		case 0 :
			numrad = 3; rvec[0] = 32; rvec[1] =  8; rvec[2] = 16; break;
		case 1 :
			numrad = 3; rvec[0] = 16; rvec[1] = 16; rvec[2] = 16; break;
		case 2 :
			numrad = 3; rvec[0] =  8; rvec[1] = 32; rvec[2] = 16; break;
		case 3 :
			numrad = 3; rvec[0] =  8; rvec[1] = 16; rvec[2] = 32; break;
		default :
			*nradices = 4;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 9 :						/* 9K */
		switch(radix_set) {
		case 0 :
			numrad = 3; rvec[0] = 36; rvec[1] =  8; rvec[2] = 16; break;
	  #ifdef USE_SSE2
		default :
			*nradices = 1;	return ERR_RADIXSET_UNAVAILABLE;
	  #else
		case 1 :
			numrad = 3; rvec[0] = 18; rvec[1] = 16; rvec[2] = 16; break;
		case 2 :
			numrad = 3; rvec[0] =  9; rvec[1] = 32; rvec[2] = 16; break;
		case 3 :
			numrad = 3; rvec[0] =  9; rvec[1] = 16; rvec[2] = 32; break;
		default :
			*nradices = 4;	return ERR_RADIXSET_UNAVAILABLE;
	  #endif
		}; break;
	case 10 :						/* 10K */
		switch(radix_set) {
		case 0 :
			numrad = 3; rvec[0] = 40; rvec[1] =  8; rvec[2] = 16; break;
		case 1 :
			numrad = 3; rvec[0] = 20; rvec[1] = 16; rvec[2] = 16; break;
	  #ifdef USE_SSE2
		default :
			*nradices = 2;	return ERR_RADIXSET_UNAVAILABLE;
	  #else
		case 2 :
			numrad = 3; rvec[0] = 10; rvec[1] = 32; rvec[2] = 16; break;
		case 3 :
			numrad = 3; rvec[0] = 10; rvec[1] = 16; rvec[2] = 32; break;
		case 4 :
			numrad = 3; rvec[0] =  5; rvec[1] = 32; rvec[2] = 32; break;
		case 5 :
			numrad = 4; rvec[0] =  5; rvec[1] =  8; rvec[2] =  8; rvec[3] = 16; break;
		default :
			*nradices = 6;	return ERR_RADIXSET_UNAVAILABLE;
	  #endif
		}; break;
	case 11 :						/* 11K */
		switch(radix_set) {
		case 0 :
			numrad = 3; rvec[0] = 44; rvec[1] =  8; rvec[2] = 16; break;
	  #ifdef USE_SSE2
		default :
			*nradices = 1;	return ERR_RADIXSET_UNAVAILABLE;
	  #else
		case 1 :
			numrad = 3; rvec[0] = 22; rvec[1] = 16; rvec[2] = 16; break;
		case 2 :
			numrad = 3; rvec[0] = 11; rvec[1] = 32; rvec[2] = 16; break;
		case 3 :
			numrad = 3; rvec[0] = 11; rvec[1] = 16; rvec[2] = 32; break;
		default :
			*nradices = 4;	return ERR_RADIXSET_UNAVAILABLE;
	  #endif
		}; break;
	case 12 :						/* 12K */
		switch(radix_set) {
		case 0 :
			numrad = 3; rvec[0] = 24; rvec[1] = 16; rvec[2] = 16; break;
		case 1 :
			numrad = 3; rvec[0] = 24; rvec[1] =  8; rvec[2] = 32; break;
		case 2 :
			numrad = 3; rvec[0] = 12; rvec[1] = 32; rvec[2] = 16; break;
		case 3 :
			numrad = 3; rvec[0] = 12; rvec[1] = 16; rvec[2] = 32; break;
	  #ifdef USE_SSE2
		default :
			*nradices = 4;	return ERR_RADIXSET_UNAVAILABLE;
	  #else
		case 4 :
			numrad = 3; rvec[0] =  6; rvec[1] = 32; rvec[2] = 32; break;
		case 5 :
			numrad = 4; rvec[0] =  6; rvec[1] =  8; rvec[2] =  8; rvec[3] = 16; break;
		default :
			*nradices = 6;	return ERR_RADIXSET_UNAVAILABLE;
	  #endif
		}; break;
	case 13 :						/* 13K */
		switch(radix_set) {
		case 0 :
			numrad = 3; rvec[0] = 52; rvec[1] =  8; rvec[2] = 16; break;
	  #ifdef USE_SSE2
		default :
			*nradices = 1;	return ERR_RADIXSET_UNAVAILABLE;
	  #else
		case 1 :
			numrad = 3; rvec[0] = 26; rvec[1] = 16; rvec[2] = 16; break;
		case 2 :
			numrad = 3; rvec[0] = 13; rvec[1] = 32; rvec[2] = 16; break;
		case 3 :
			numrad = 3; rvec[0] = 13; rvec[1] = 16; rvec[2] = 32; break;
		default :
			*nradices = 4;	return ERR_RADIXSET_UNAVAILABLE;
	  #endif
		}; break;
	case 14 :						/* 14K */
		switch(radix_set) {
		case 0 :
			numrad = 3; rvec[0] = 56; rvec[1] =  8; rvec[2] = 16; break;
		case 1 :
			numrad = 3; rvec[0] = 28; rvec[1] = 16; rvec[2] = 16; break;
	  #ifdef USE_SSE2
		default :
			*nradices = 2;	return ERR_RADIXSET_UNAVAILABLE;
	  #else
		case 2 :
			numrad = 3; rvec[0] = 14; rvec[1] = 32; rvec[2] = 16; break;
		case 3 :
			numrad = 3; rvec[0] = 14; rvec[1] = 16; rvec[2] = 32; break;
		case 4 :
			numrad = 3; rvec[0] =  7; rvec[1] = 32; rvec[2] = 32; break;
		case 5 :
			numrad = 4; rvec[0] =  7; rvec[1] =  8; rvec[2] =  8; rvec[3] = 16; break;
		default :
			*nradices = 6;	return ERR_RADIXSET_UNAVAILABLE;
	  #endif
		}; break;
	case 15 :						/* 15K */
		switch(radix_set) {
		case 0 :
			numrad = 3; rvec[0] = 60; rvec[1] =  8; rvec[2] = 16; break;
	  #ifdef USE_SSE2
		default :
			*nradices = 1;	return ERR_RADIXSET_UNAVAILABLE;
	  #else
		case 1 :
			numrad = 3; rvec[0] = 30; rvec[1] = 16; rvec[2] = 16; break;
		case 2 :
			numrad = 3; rvec[0] = 15; rvec[1] = 32; rvec[2] = 16; break;
		case 3 :
			numrad = 3; rvec[0] = 15; rvec[1] = 16; rvec[2] = 32; break;
		default :
			*nradices = 4;	return ERR_RADIXSET_UNAVAILABLE;
	  #endif
		}; break;
	case 16 :						/* 16K */
		switch(radix_set) {
		case 0 :
			numrad = 3; rvec[0] = 32; rvec[1] = 16; rvec[2] = 16; break;
		case 1 :
			numrad = 3; rvec[0] = 16; rvec[1] = 32; rvec[2] = 16; break;
		case 2 :
			numrad = 3; rvec[0] = 16; rvec[1] = 16; rvec[2] = 32; break;
	  #ifdef USE_SSE2
		default :
			*nradices = 3;	return ERR_RADIXSET_UNAVAILABLE;
	  #else
		case 3 :
			numrad = 3; rvec[0] =  8; rvec[1] = 32; rvec[2] = 32; break;
		case 4 :
			numrad = 4; rvec[0] =  8; rvec[1] =  8; rvec[2] =  8; rvec[3] = 16; break;
		default :
			*nradices = 5;	return ERR_RADIXSET_UNAVAILABLE;
	  #endif
		}; break;
	case 18 :						/* 18K */
		switch(radix_set) {
		case 0 :
			numrad = 3; rvec[0] = 36; rvec[1] = 16; rvec[2] = 16; break;
		case 1 :
			numrad = 3; rvec[0] = 36; rvec[1] =  8; rvec[2] = 32; break;
	  #ifdef USE_SSE2
		default :
			*nradices = 2;	return ERR_RADIXSET_UNAVAILABLE;
	  #else
		case 2 :
			numrad = 3; rvec[0] = 18; rvec[1] = 32; rvec[2] = 16; break;
		case 3 :
			numrad = 3; rvec[0] = 18; rvec[1] = 16; rvec[2] = 32; break;
		case 4 :
			numrad = 3; rvec[0] =  9; rvec[1] = 32; rvec[2] = 32; break;
		case 5 :
			numrad = 4; rvec[0] =  9; rvec[1] =  8; rvec[2] =  8; rvec[3] = 16; break;
		default :
			*nradices = 6;	return ERR_RADIXSET_UNAVAILABLE;
	  #endif
		}; break;
	case 20 :						/* 20K */
		switch(radix_set) {
		case 0 :
			numrad = 3; rvec[0] = 40; rvec[1] = 16; rvec[2] = 16; break;
		case 1 :
			numrad = 3; rvec[0] = 40; rvec[1] =  8; rvec[2] = 32; break;
		case 2 :
			numrad = 3; rvec[0] = 20; rvec[1] = 16; rvec[2] = 32; break;
		case 3 :
			numrad = 3; rvec[0] = 20; rvec[1] = 32; rvec[2] = 16; break;
	  #ifdef USE_SSE2
		default :
			*nradices = 4;	return ERR_RADIXSET_UNAVAILABLE;
	  #else
		case 4 :
			numrad = 3; rvec[0] = 10; rvec[1] = 32; rvec[2] = 32; break;
		case 5 :
			numrad = 4; rvec[0] =  5; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; break;
		default :
			*nradices = 6;	return ERR_RADIXSET_UNAVAILABLE;
	  #endif
		}; break;
	case 22 :						/* 22K */
		switch(radix_set) {
		case 0 :
			numrad = 3; rvec[0] = 44; rvec[1] = 16; rvec[2] = 16; break;
		case 1 :
			numrad = 3; rvec[0] = 44; rvec[1] =  8; rvec[2] = 32; break;
	  #ifdef USE_SSE2
		default :
			*nradices = 2;	return ERR_RADIXSET_UNAVAILABLE;
	  #else
		case 2 :
			numrad = 3; rvec[0] = 22; rvec[1] = 32; rvec[2] = 16; break;
		case 3 :
			numrad = 3; rvec[0] = 22; rvec[1] = 16; rvec[2] = 32; break;
		case 4 :
			numrad = 3; rvec[0] = 11; rvec[1] = 32; rvec[2] = 32; break;
		case 5 :
			numrad = 4; rvec[0] = 11; rvec[1] =  8; rvec[2] =  8; rvec[3] = 16; break;
		default :
			*nradices = 6;	return ERR_RADIXSET_UNAVAILABLE;
	  #endif
		}; break;
	case 24 :						/* 24K */
		switch(radix_set) {
		case 0 :
			numrad = 3; rvec[0] = 48; rvec[1] = 16; rvec[2] = 16; break;
		case 1 :
			numrad = 3; rvec[0] = 48; rvec[1] =  8; rvec[2] = 32; break;
		case 2 :
			numrad = 3; rvec[0] = 24; rvec[1] = 16; rvec[2] = 32; break;
		case 3 :
			numrad = 3; rvec[0] = 24; rvec[1] = 32; rvec[2] = 16; break;
	  #ifdef USE_SSE2
		default :
			*nradices = 4;	return ERR_RADIXSET_UNAVAILABLE;
	  #else
		case 4 :
			numrad = 3; rvec[0] = 12; rvec[1] = 32; rvec[2] = 32; break;
		case 5 :
			numrad = 4; rvec[0] =  6; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; break;
		default :
			*nradices = 6;	return ERR_RADIXSET_UNAVAILABLE;
	  #endif
		}; break;
	case 26 :						/* 26K */
		switch(radix_set) {
		case 0 :
			numrad = 3; rvec[0] = 52; rvec[1] = 16; rvec[2] = 16; break;
		case 1 :
			numrad = 3; rvec[0] = 52; rvec[1] =  8; rvec[2] = 32; break;
	  #ifdef USE_SSE2
		default :
			*nradices = 2;	return ERR_RADIXSET_UNAVAILABLE;
	  #else
		case 2 :
			numrad = 3; rvec[0] = 26; rvec[1] = 32; rvec[2] = 16; break;
		case 3 :
			numrad = 3; rvec[0] = 26; rvec[1] = 16; rvec[2] = 32; break;
		case 4 :
			numrad = 3; rvec[0] = 13; rvec[1] = 32; rvec[2] = 32; break;
		case 5 :
			numrad = 4; rvec[0] = 13; rvec[1] =  8; rvec[2] =  8; rvec[3] = 16; break;
		default :
			*nradices = 6;	return ERR_RADIXSET_UNAVAILABLE;
	  #endif
		}; break;
	case 28 :						/* 28K */
		switch(radix_set) {
		case 0 :
			numrad = 3; rvec[0] = 56; rvec[1] = 16; rvec[2] = 16; break;
		case 1 :
			numrad = 3; rvec[0] = 56; rvec[1] =  8; rvec[2] = 32; break;
		case 2 :
			numrad = 3; rvec[0] = 28; rvec[1] = 16; rvec[2] = 32; break;
		case 3 :
			numrad = 3; rvec[0] = 28; rvec[1] = 32; rvec[2] = 16; break;
	  #ifdef USE_SSE2
		default :
			*nradices = 4;	return ERR_RADIXSET_UNAVAILABLE;
	  #else
		case 4 :
			numrad = 3; rvec[0] = 14; rvec[1] = 32; rvec[2] = 32; break;
		case 5 :
			numrad = 4; rvec[0] =  7; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; break;
		default :
			*nradices = 6;	return ERR_RADIXSET_UNAVAILABLE;
	  #endif
		}; break;
	case 30 :						/* 30K */
		switch(radix_set) {
		case 0 :
			numrad = 3; rvec[0] = 60; rvec[1] = 16; rvec[2] = 16; break;
		case 1 :
			numrad = 3; rvec[0] = 60; rvec[1] =  8; rvec[2] = 32; break;
	  #ifdef USE_SSE2
		default :
			*nradices = 2;	return ERR_RADIXSET_UNAVAILABLE;
	  #else
		case 2 :
			numrad = 3; rvec[0] = 30; rvec[1] = 32; rvec[2] = 16; break;
		case 3 :
			numrad = 3; rvec[0] = 30; rvec[1] = 16; rvec[2] = 32; break;
		case 4 :
			numrad = 3; rvec[0] = 15; rvec[1] = 32; rvec[2] = 32; break;
		case 5 :
			numrad = 4; rvec[0] = 15; rvec[1] =  8; rvec[2] =  8; rvec[3] = 16; break;
		default :
			*nradices = 6;	return ERR_RADIXSET_UNAVAILABLE;
	  #endif
		}; break;
	case 32 :						/* 32K */
		switch(radix_set) {
		case 0 :
			numrad = 3; rvec[0] = 32; rvec[1] = 32; rvec[2] = 16; break;
		case 1 :
			numrad = 3; rvec[0] = 32; rvec[1] = 16; rvec[2] = 32; break;
		case 2 :
			numrad = 3; rvec[0] = 16; rvec[1] = 32; rvec[2] = 32; break;
		case 3 :
			numrad = 4; rvec[0] = 16; rvec[1] =  8; rvec[2] =  8; rvec[3] = 16; break;
	  #ifdef USE_SSE2
		default :
			*nradices = 4;	return ERR_RADIXSET_UNAVAILABLE;
	  #else
		case 4 :
			numrad = 4; rvec[0] =  8; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; break;
		default :
			*nradices = 5;	return ERR_RADIXSET_UNAVAILABLE;
	  #endif
		}; break;
	case 36 :						/* 36K */
		switch(radix_set) {
		case 0 :
			numrad = 3; rvec[0] = 36; rvec[1] = 32; rvec[2] = 16; break;
		case 1 :
			numrad = 3; rvec[0] = 36; rvec[1] = 16; rvec[2] = 32; break;
	  #ifdef USE_SSE2
		default :
			*nradices = 2;	return ERR_RADIXSET_UNAVAILABLE;
	  #else
		case 2 :
			numrad = 3; rvec[0] = 18; rvec[1] = 32; rvec[2] = 32; break;
		case 3 :
			numrad = 4; rvec[0] =  9; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; break;
		default :
			*nradices = 4;	return ERR_RADIXSET_UNAVAILABLE;
	  #endif
		}; break;
	case 40 :						/* 40K */
		switch(radix_set) {
		case 0 :
			numrad = 3; rvec[0] = 40; rvec[1] = 32; rvec[2] = 16; break;
		case 1 :
			numrad = 3; rvec[0] = 40; rvec[1] = 16; rvec[2] = 32; break;
		case 2 :
			numrad = 3; rvec[0] = 20; rvec[1] = 32; rvec[2] = 32; break;
	  #ifdef USE_SSE2
		default :
			*nradices = 3;	return ERR_RADIXSET_UNAVAILABLE;
	  #else
		case 3 :
			numrad = 4; rvec[0] = 10; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; break;
		default :
			*nradices = 4;	return ERR_RADIXSET_UNAVAILABLE;
	  #endif
		}; break;
	case 44 :						/* 44K */
		switch(radix_set) {
		case 0 :
			numrad = 3; rvec[0] = 44; rvec[1] = 32; rvec[2] = 16; break;
		case 1 :
			numrad = 3; rvec[0] = 44; rvec[1] = 16; rvec[2] = 32; break;
	  #ifdef USE_SSE2
		default :
			*nradices = 2;	return ERR_RADIXSET_UNAVAILABLE;
	  #else
		case 2 :
			numrad = 3; rvec[0] = 22; rvec[1] = 32; rvec[2] = 32; break;
		case 3 :
			numrad = 4; rvec[0] = 11; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; break;
		default :
			*nradices = 4;	return ERR_RADIXSET_UNAVAILABLE;
	  #endif
		}; break;
	case 48 :						/* 48K */
		switch(radix_set) {
		case 0 :
			numrad = 3; rvec[0] = 48; rvec[1] = 32; rvec[2] = 16; break;
		case 1 :
			numrad = 3; rvec[0] = 48; rvec[1] = 16; rvec[2] = 32; break;
		case 2 :
			numrad = 3; rvec[0] = 24; rvec[1] = 32; rvec[2] = 32; break;
	  #ifdef USE_SSE2
		default :
			*nradices = 3;	return ERR_RADIXSET_UNAVAILABLE;
	  #else
		case 3 :
			numrad = 4; rvec[0] = 12; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; break;
		default :
			*nradices = 4;	return ERR_RADIXSET_UNAVAILABLE;
	  #endif
		}; break;
	case 52 :						/* 52K */
		switch(radix_set) {
		case 0 :
			numrad = 3; rvec[0] = 52; rvec[1] = 32; rvec[2] = 16; break;
		case 1 :
			numrad = 3; rvec[0] = 52; rvec[1] = 16; rvec[2] = 32; break;
	  #ifdef USE_SSE2
		default :
			*nradices = 2;	return ERR_RADIXSET_UNAVAILABLE;
	  #else
		case 2 :
			numrad = 3; rvec[0] = 26; rvec[1] = 32; rvec[2] = 32; break;
		case 3 :
			numrad = 4; rvec[0] = 13; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; break;
		default :
			*nradices = 4;	return ERR_RADIXSET_UNAVAILABLE;
	  #endif
		}; break;
	case 56 :						/* 56K */
		switch(radix_set) {
		case 0 :
			numrad = 3; rvec[0] = 56; rvec[1] = 32; rvec[2] = 16; break;
		case 1 :
			numrad = 3; rvec[0] = 56; rvec[1] = 16; rvec[2] = 32; break;
		case 2 :
			numrad = 3; rvec[0] = 28; rvec[1] = 32; rvec[2] = 32; break;
	  #ifdef USE_SSE2
		default :
			*nradices = 3;	return ERR_RADIXSET_UNAVAILABLE;
	  #else
		case 3 :
			numrad = 4; rvec[0] = 14; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; break;
		default :
			*nradices = 4;	return ERR_RADIXSET_UNAVAILABLE;
	  #endif
		}; break;
	case 60 :						/* 60K */
		switch(radix_set) {
		case 0 :
			numrad = 3; rvec[0] = 60; rvec[1] = 32; rvec[2] = 16; break;
		case 1 :
			numrad = 3; rvec[0] = 60; rvec[1] = 16; rvec[2] = 32; break;
	  #ifdef USE_SSE2
		default :
			*nradices = 2;	return ERR_RADIXSET_UNAVAILABLE;
	  #else
		case 2 :
			numrad = 3; rvec[0] = 30; rvec[1] = 32; rvec[2] = 32; break;
		case 3 :
			numrad = 4; rvec[0] = 15; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; break;
		default :
			*nradices = 4;	return ERR_RADIXSET_UNAVAILABLE;
	  #endif
		}; break;
	case 64 :						/* 64K */
		switch(radix_set) {
		case 0 :
			numrad = 3; rvec[0] = 64; rvec[1] = 32; rvec[2] = 16; break;
		case 1 :
			numrad = 3; rvec[0] = 64; rvec[1] = 16; rvec[2] = 32; break;
		case 2 :
			numrad = 3; rvec[0] = 32; rvec[1] = 32; rvec[2] = 32; break;
		case 3 :
			numrad = 4; rvec[0] = 32; rvec[1] =  8; rvec[2] =  8; rvec[3] = 16; break;
		case 4 :
			numrad = 4; rvec[0] = 16; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; break;
	  #ifdef USE_SSE2
		default :
			*nradices = 5;	return ERR_RADIXSET_UNAVAILABLE;
	  #else
		case 5 :
			numrad = 4; rvec[0] =  8; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; break;
		default :
			*nradices = 6;	return ERR_RADIXSET_UNAVAILABLE;
	  #endif
		}; break;
	case 72 :						/* 72K */
		switch(radix_set) {
		case 0 :
			numrad = 3; rvec[0] = 36; rvec[1] = 32; rvec[2] = 32; break;
		case 1 :
			numrad = 4; rvec[0] = 36; rvec[1] =  8; rvec[2] =  8; rvec[3] = 16; break;
	  #ifdef USE_SSE2
		default :
			*nradices = 2;	return ERR_RADIXSET_UNAVAILABLE;
	  #else
		case 2 :
			numrad = 4; rvec[0] = 18; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; break;
		case 3 :
			numrad = 4; rvec[0] =  9; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; break;
		default :
			*nradices = 4;	return ERR_RADIXSET_UNAVAILABLE;
	  #endif
		}; break;
	case 80 :						/* 80K */
		switch(radix_set) {
		case 0 :
			numrad = 3; rvec[0] = 40; rvec[1] = 32; rvec[2] = 32; break;
		case 1 :
			numrad = 4; rvec[0] = 40; rvec[1] =  8; rvec[2] =  8; rvec[3] = 16; break;
		case 2 :
			numrad = 4; rvec[0] = 20; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; break;
	  #ifdef USE_SSE2
		default :
			*nradices = 3;	return ERR_RADIXSET_UNAVAILABLE;
	  #else
		case 3 :
			numrad = 4; rvec[0] = 10; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; break;
		default :
			*nradices = 4;	return ERR_RADIXSET_UNAVAILABLE;
	  #endif
		}; break;
	case 88 :						/* 88K */
		switch(radix_set) {
		case 0 :
			numrad = 3; rvec[0] = 44; rvec[1] = 32; rvec[2] = 32; break;
		case 1 :
			numrad = 4; rvec[0] = 44; rvec[1] =  8; rvec[2] =  8; rvec[3] = 16; break;
	  #ifdef USE_SSE2
		default :
			*nradices = 2;	return ERR_RADIXSET_UNAVAILABLE;
	  #else
		case 2 :
			numrad = 4; rvec[0] = 22; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; break;
		case 3 :
			numrad = 4; rvec[0] = 11; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; break;
		default :
			*nradices = 4;	return ERR_RADIXSET_UNAVAILABLE;
	  #endif
		}; break;
	case 96 :						/* 96K */
		switch(radix_set) {
		case 0 :
			numrad = 3; rvec[0] = 48; rvec[1] = 32; rvec[2] = 32; break;
		case 1 :
			numrad = 4; rvec[0] = 24; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; break;
	  #ifdef USE_SSE2
		default :
			*nradices = 2;	return ERR_RADIXSET_UNAVAILABLE;
	  #else
		case 2 :
			numrad = 4; rvec[0] = 12; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; break;
		case 3 :
			numrad = 4; rvec[0] =  6; rvec[1] = 16; rvec[2] = 32; rvec[3] = 16; break;
		case 4 :
			numrad = 4; rvec[0] =  6; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; break;
		default :
			*nradices = 5;	return ERR_RADIXSET_UNAVAILABLE;
	  #endif
		}; break;
	case 104 :						/* 104K */
		switch(radix_set) {
		case 0 :
			numrad = 3; rvec[0] = 52; rvec[1] = 32; rvec[2] = 32; break;
		case 1 :
			numrad = 4; rvec[0] = 52; rvec[1] =  8; rvec[2] =  8; rvec[3] = 16; break;
	  #ifdef USE_SSE2
		default :
			*nradices = 2;	return ERR_RADIXSET_UNAVAILABLE;
	  #else
		case 2 :
			numrad = 4; rvec[0] = 26; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; break;
		case 3 :
			numrad = 4; rvec[0] = 13; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; break;
		default :
			*nradices = 4;	return ERR_RADIXSET_UNAVAILABLE;
	  #endif
		}; break;
	case 112 :						/* 112K */
		switch(radix_set) {
		case 0 :
			numrad = 3; rvec[0] = 56; rvec[1] = 32; rvec[2] = 32; break;
		case 1 :
			numrad = 4; rvec[0] = 28; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; break;
	  #ifdef USE_SSE2
		default :
			*nradices = 2;	return ERR_RADIXSET_UNAVAILABLE;
	  #else
		case 2 :
			numrad = 4; rvec[0] = 14; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; break;
		case 3 :
			numrad = 4; rvec[0] =  7; rvec[1] = 16; rvec[2] = 32; rvec[3] = 16; break;
		case 4 :
			numrad = 4; rvec[0] =  7; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; break;
		default :
			*nradices = 5;	return ERR_RADIXSET_UNAVAILABLE;
	  #endif
		}; break;
	case 120 :						/* 120K */
		switch(radix_set) {
		case 0 :
			numrad = 3; rvec[0] = 60; rvec[1] = 32; rvec[2] = 32; break;
		case 1 :
			numrad = 4; rvec[0] = 60; rvec[1] =  8; rvec[2] =  8; rvec[3] = 16; break;
	  #ifdef USE_SSE2
		default :
			*nradices = 2;	return ERR_RADIXSET_UNAVAILABLE;
	  #else
		case 2 :
			numrad = 4; rvec[0] = 30; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; break;
		case 3 :
			numrad = 4; rvec[0] = 15; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; break;
		default :
			*nradices = 4;	return ERR_RADIXSET_UNAVAILABLE;
	  #endif
		}; break;
	case 128 :						/* 128K */
		switch(radix_set) {
		case 0 :
			numrad = 3; rvec[0] =128; rvec[1] = 16; rvec[2] = 32; break;
		case 1 :
			numrad = 3; rvec[0] =128; rvec[1] = 32; rvec[2] = 16; break;
		case 2 :
			numrad = 3; rvec[0] = 64; rvec[1] = 32; rvec[2] = 32; break;
		case 3 :
			numrad = 4; rvec[0] = 32; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; break;
		case 4 :
			numrad = 4; rvec[0] = 16; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; break;
	  #ifdef USE_SSE2
		default :
			*nradices = 5;	return ERR_RADIXSET_UNAVAILABLE;
	  #else
		case 5 :
			numrad = 4; rvec[0] =  8; rvec[1] = 16; rvec[2] = 32; rvec[3] = 16; break;
		case 6 :
			numrad = 4; rvec[0] =  8; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; break;
		default :
			*nradices = 7;	return ERR_RADIXSET_UNAVAILABLE;
	  #endif
		}; break;
	case 144 :						/* 144K */
		switch(radix_set) {
		case 0 :
			numrad = 3; rvec[0] =144; rvec[1] = 16; rvec[2] = 32; break;
		case 1 :
			numrad = 3; rvec[0] =144; rvec[1] = 32; rvec[2] = 16; break;
		case 2 :
			numrad = 4; rvec[0] = 36; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; break;
	  #ifdef USE_SSE2
		default :
			*nradices = 3;	return ERR_RADIXSET_UNAVAILABLE;
	  #else
		case 3 :
			numrad = 4; rvec[0] = 18; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; break;
		case 4 :
			numrad = 4; rvec[0] =  9; rvec[1] = 16; rvec[2] = 32; rvec[3] = 16; break;
		case 5 :
			numrad = 4; rvec[0] =  9; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; break;
		default :
			*nradices = 6;	return ERR_RADIXSET_UNAVAILABLE;
	  #endif
		}; break;
	case 160 :						/* 160K */
		switch(radix_set) {
		case 0 :
			numrad = 3; rvec[0] =160; rvec[1] = 16; rvec[2] = 32; break;
		case 1 :
			numrad = 3; rvec[0] =160; rvec[1] = 32; rvec[2] = 16; break;
		case 2 :
			numrad = 4; rvec[0] = 40; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; break;
		case 3 :
			numrad = 4; rvec[0] = 20; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; break;
	  #ifdef USE_SSE2
		default :
			*nradices = 4;	return ERR_RADIXSET_UNAVAILABLE;
	  #else
		case 4 :
			numrad = 4; rvec[0] = 10; rvec[1] = 16; rvec[2] = 32; rvec[3] = 16; break;
		case 5 :
			numrad = 4; rvec[0] = 10; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; break;
		default :
			*nradices = 6;	return ERR_RADIXSET_UNAVAILABLE;
	  #endif
		}; break;
	case 176 :						/* 176K */
		switch(radix_set) {
		case 0 :
			numrad = 3; rvec[0] =176; rvec[1] = 16; rvec[2] = 32; break;
		case 1 :
			numrad = 3; rvec[0] =176; rvec[1] = 32; rvec[2] = 16; break;
		case 2 :
			numrad = 4; rvec[0] = 44; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; break;
	  #ifdef USE_SSE2
		default :
			*nradices = 3;	return ERR_RADIXSET_UNAVAILABLE;
	  #else
		case 3 :
			numrad = 4; rvec[0] = 22; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; break;
		case 4 :
			numrad = 4; rvec[0] = 11; rvec[1] = 16; rvec[2] = 32; rvec[3] = 16; break;
		case 5 :
			numrad = 4; rvec[0] = 11; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; break;
		default :
			*nradices = 6;	return ERR_RADIXSET_UNAVAILABLE;
	  #endif
		}; break;
	case 192 :						/* 192K */
		switch(radix_set) {
		case 0 :
			numrad = 3; rvec[0] =192; rvec[1] = 16; rvec[2] = 32; break;
		case 1 :
			numrad = 3; rvec[0] =192; rvec[1] = 32; rvec[2] = 16; break;
		case 2 :
			numrad = 4; rvec[0] = 48; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; break;
		case 3 :
			numrad = 4; rvec[0] = 24; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; break;
	  #ifdef USE_SSE2
		default :
			*nradices = 4;	return ERR_RADIXSET_UNAVAILABLE;
	  #else
		case 4 :
			numrad = 4; rvec[0] = 12; rvec[1] = 16; rvec[2] = 32; rvec[3] = 16; break;
		case 5 :
			numrad = 4; rvec[0] = 12; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; break;
		default :
			*nradices = 6;	return ERR_RADIXSET_UNAVAILABLE;
	  #endif
		}; break;
	case 208 :						/* 208K */
		switch(radix_set) {
		case 0 :
			numrad = 3; rvec[0] =208; rvec[1] = 16; rvec[2] = 32; break;
		case 1 :
			numrad = 3; rvec[0] =208; rvec[1] = 32; rvec[2] = 16; break;
		case 2 :
			numrad = 4; rvec[0] = 52; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; break;
	  #ifdef USE_SSE2
		default :
			*nradices = 3;	return ERR_RADIXSET_UNAVAILABLE;
	  #else
		case 3 :
			numrad = 4; rvec[0] = 26; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; break;
		case 4 :
			numrad = 4; rvec[0] = 13; rvec[1] = 16; rvec[2] = 32; rvec[3] = 16; break;
		case 5 :
			numrad = 4; rvec[0] = 13; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; break;
		default :
			*nradices = 6;	return ERR_RADIXSET_UNAVAILABLE;
	  #endif
		}; break;
	case 224 :						/* 224K */
		switch(radix_set) {
		case 0 :
			numrad = 3; rvec[0] =224; rvec[1] = 16; rvec[2] = 32; break;
		case 1 :
			numrad = 3; rvec[0] =224; rvec[1] = 32; rvec[2] = 16; break;
		case 2 :
			numrad = 4; rvec[0] = 56; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; break;
		case 3 :
			numrad = 4; rvec[0] = 28; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; break;
	  #ifdef USE_SSE2
		default :
			*nradices = 4;	return ERR_RADIXSET_UNAVAILABLE;
	  #else
		case 4 :
			numrad = 4; rvec[0] = 14; rvec[1] = 16; rvec[2] = 32; rvec[3] = 16; break;
		case 5 :
			numrad = 4; rvec[0] = 14; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; break;
		default :
			*nradices = 6;	return ERR_RADIXSET_UNAVAILABLE;
	  #endif
		}; break;
	case 240 :						/* 240K */
		switch(radix_set) {
		case 0 :
			numrad = 3; rvec[0] =240; rvec[1] = 16; rvec[2] = 32; break;
		case 1 :
			numrad = 3; rvec[0] =240; rvec[1] = 32; rvec[2] = 16; break;
		case 2 :
			numrad = 4; rvec[0] = 60; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; break;
	  #ifdef USE_SSE2
		default :
			*nradices = 3;	return ERR_RADIXSET_UNAVAILABLE;
	  #else
		case 3 :
			numrad = 4; rvec[0] = 30; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; break;
		case 4 :
			numrad = 4; rvec[0] = 15; rvec[1] = 16; rvec[2] = 32; rvec[3] = 16; break;
		case 5 :
			numrad = 4; rvec[0] = 15; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; break;
		default :
			*nradices = 6;	return ERR_RADIXSET_UNAVAILABLE;
	  #endif
		}; break;
	case 256 :						/* 256K */
		switch(radix_set) {
		case 0 :
			numrad = 3; rvec[0] =256; rvec[1] = 16; rvec[2] = 32; break;
		case 1 :
			numrad = 3; rvec[0] =128; rvec[1] = 32; rvec[2] = 32; break;
		case 2 :
			numrad = 4; rvec[0] = 64; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; break;
		case 3 :
			numrad = 4; rvec[0] = 32; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; break;
		case 4 :
			numrad = 4; rvec[0] = 16; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; break;
		case 5 :
			numrad = 4; rvec[0] = 16; rvec[1] = 32; rvec[2] = 16; rvec[3] = 16; break;
		case 6 :
			numrad = 5; rvec[0] = 16; rvec[1] =  8; rvec[2] =  8; rvec[3] =  8; rvec[4] = 16; break;
	  #ifdef USE_SSE2
		default :
			*nradices = 7;	return ERR_RADIXSET_UNAVAILABLE;
	  #else
		case 7 :
			numrad = 4; rvec[0] =  8; rvec[1] = 32; rvec[2] = 32; rvec[3] = 16; break;
		case 8 :
			numrad = 4; rvec[0] =  8; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; break;
		case 9 :
			numrad = 5; rvec[0] =  8; rvec[1] =  8; rvec[2] =  8; rvec[3] = 16; rvec[4] = 16; break;
		default :
			*nradices = 10;	return ERR_RADIXSET_UNAVAILABLE;
	  #endif
		}; break;
	case 288 :						/* 288K */
		switch(radix_set) {
		case 0 :
			numrad = 3; rvec[0] =144; rvec[1] = 32; rvec[2] = 32; break;
		case 1 :
			numrad = 4; rvec[0] =144; rvec[1] =  8; rvec[2] =  8; rvec[3] = 16; break;
		case 2 :
			numrad = 4; rvec[0] = 36; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; break;
		default :
			*nradices = 3;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 320 :						/* 320K */
		switch(radix_set) {
		case 0 :
			numrad = 3; rvec[0] =160; rvec[1] = 32; rvec[2] = 32; break;
		case 1 :
			numrad = 4; rvec[0] =160; rvec[1] =  8; rvec[2] =  8; rvec[3] = 16; break;
		case 2 :
			numrad = 4; rvec[0] = 40; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; break;
		case 3 :
			numrad = 4; rvec[0] = 40; rvec[1] =  8; rvec[2] = 16; rvec[3] = 32; break;
		case 4 :
			numrad = 4; rvec[0] = 20; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; break;
		case 5 :
			numrad = 4; rvec[0] = 20; rvec[1] = 32; rvec[2] = 16; rvec[3] = 16; break;
		default :
			*nradices = 6;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 352 :						/* 352K */
		switch(radix_set) {
		case 0 :
			numrad = 3; rvec[0] =176; rvec[1] = 32; rvec[2] = 32; break;
		case 1 :
			numrad = 4; rvec[0] =176; rvec[1] =  8; rvec[2] =  8; rvec[3] = 16; break;
		case 2 :
			numrad = 4; rvec[0] = 44; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; break;
		default :
			*nradices = 3;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 384 :						/* 384K */
		switch(radix_set) {
		case 0 :
			numrad = 3; rvec[0] =192; rvec[1] = 32; rvec[2] = 32; break;
		case 1 :
			numrad = 4; rvec[0] = 48; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; break;
		case 2 :
			numrad = 4; rvec[0] = 48; rvec[1] =  8; rvec[2] = 16; rvec[3] = 32; break;
		case 3 :
			numrad = 4; rvec[0] = 24; rvec[1] = 32; rvec[2] = 16; rvec[3] = 16; break;
		case 4 :
			numrad = 5; rvec[0] = 24; rvec[1] =  8; rvec[2] =  8; rvec[3] =  8; rvec[4] = 16; break;
		default :
			*nradices = 5;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
  	case 416 :						/* 416K */
		switch(radix_set) {
		case 0 :
			numrad = 3; rvec[0] =208; rvec[1] = 32; rvec[2] = 32; break;
		case 1 :
			numrad = 4; rvec[0] =208; rvec[1] =  8; rvec[2] =  8; rvec[3] = 16; break;
		case 2 :
			numrad = 4; rvec[0] = 52; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; break;
		default :
			*nradices = 3;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 448 :						/* 448K */
		switch(radix_set) {
		case 0 :
			numrad = 3; rvec[0] =224; rvec[1] = 32; rvec[2] = 32; break;
		case 1 :
			numrad = 4; rvec[0] =224; rvec[1] =  8; rvec[2] =  8; rvec[3] = 16; break;
		case 2 :
			numrad = 4; rvec[0] = 56; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; break;
		case 3 :
			numrad = 4; rvec[0] = 28; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; break;
		case 4 :
			numrad = 4; rvec[0] = 28; rvec[1] = 32; rvec[2] = 16; rvec[3] = 16; break;
		case 5 :
			numrad = 5; rvec[0] = 28; rvec[1] =  8; rvec[2] =  8; rvec[3] =  8; rvec[4] = 16; break;
		default :
			*nradices = 6;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 480 :						/* 480K */
		switch(radix_set) {
		case 0 :
			numrad = 3; rvec[0] =240; rvec[1] = 32; rvec[2] = 32; break;
		case 1 :
			numrad = 4; rvec[0] =240; rvec[1] =  8; rvec[2] =  8; rvec[3] = 16; break;
		case 2 :
			numrad = 4; rvec[0] = 60; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; break;
		default :
			*nradices = 3;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 512 :						/* 512K */
		switch(radix_set) {
		case 0 :
			numrad = 3; rvec[0] =256; rvec[1] = 32; rvec[2] = 32; break;
		case 1 :
			numrad = 4; rvec[0] =128; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; break;
		case 2 :
			numrad = 4; rvec[0] = 64; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; break;
		case 3 :
			numrad = 4; rvec[0] = 32; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; break;
		case 4 :
			numrad = 4; rvec[0] = 32; rvec[1] = 32; rvec[2] = 16; rvec[3] = 16; break;
		case 5 :
			numrad = 5; rvec[0] = 32; rvec[1] =  8; rvec[2] =  8; rvec[3] =  8; rvec[4] = 16; break;
		case 6 :
			numrad = 4; rvec[0] = 16; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; break;
		case 7 :
			numrad = 4; rvec[0] = 16; rvec[1] = 32; rvec[2] = 32; rvec[3] = 16; break;
		case 8 :
			numrad = 5; rvec[0] = 16; rvec[1] =  8; rvec[2] =  8; rvec[3] = 16; rvec[4] = 16; break;
		default :
			*nradices = 9;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 576 :						/* 576K */
		switch(radix_set) {
		case 0 :
			numrad = 3; rvec[0] =288; rvec[1] = 32; rvec[2] = 32; break;
		case 1 :
			numrad = 4; rvec[0] =144; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; break;
		case 2 :
			numrad = 4; rvec[0] = 36; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; break;
		case 3 :
			numrad = 4; rvec[0] = 36; rvec[1] = 32; rvec[2] = 16; rvec[3] = 16; break;
		case 4 :
			numrad = 5; rvec[0] = 36; rvec[1] =  8; rvec[2] =  8; rvec[3] =  8; rvec[4] = 16; break;
		default :
			*nradices = 5;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 640 :						/* 640K */
		switch(radix_set) {
		case 0 :
			numrad = 3; rvec[0] =320; rvec[1] = 32; rvec[2] = 32; break;
		case 1 :
			numrad = 4; rvec[0] =160; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; break;
		case 2 :
			numrad = 4; rvec[0] = 40; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; break;
		case 3 :
			numrad = 4; rvec[0] = 40; rvec[1] = 32; rvec[2] = 16; rvec[3] = 16; break;
		case 4 :
			numrad = 5; rvec[0] = 40; rvec[1] =  8; rvec[2] =  8; rvec[3] =  8; rvec[4] = 16; break;
		case 5 :
			numrad = 4; rvec[0] = 20; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; break;
		case 6 :
			numrad = 4; rvec[0] = 20; rvec[1] = 32; rvec[2] = 32; rvec[3] = 16; break;
		default :
			*nradices = 7;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 704 :						/* 704K */
		switch(radix_set) {
		case 0 :
			numrad = 3; rvec[0] =352; rvec[1] = 32; rvec[2] = 32; break;
		case 1 :
			numrad = 4; rvec[0] =176; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; break;
		case 2 :
			numrad = 4; rvec[0] = 44; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; break;
		case 3 :
			numrad = 4; rvec[0] = 44; rvec[1] = 32; rvec[2] = 16; rvec[3] = 16; break;
		case 4 :
			numrad = 5; rvec[0] = 44; rvec[1] =  8; rvec[2] =  8; rvec[3] =  8; rvec[4] = 16; break;
		default :
			*nradices = 5;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 768 :						/* 768K */
		switch(radix_set) {	// 768K too small to use rvec[0] = 768 in || mode - CY_THREADS > 1 fails the (n_div_nwt%CY_THREADS == 0) assertion check.
		case 0 :
			numrad = 3; rvec[0] =768; rvec[1] = 16; rvec[2] = 32; break;
		case 1 :
			numrad = 4; rvec[0] =192; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; break;
		case 2 :
			numrad = 4; rvec[0] = 48; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; break;
		case 3 :
			numrad = 4; rvec[0] = 48; rvec[1] = 32; rvec[2] = 16; rvec[3] = 16; break;
		case 4 :
			numrad = 5; rvec[0] = 48; rvec[1] =  8; rvec[2] =  8; rvec[3] =  8; rvec[4] = 16; break;
		case 5 :
			numrad = 4; rvec[0] = 24; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; break;
		case 6 :
			numrad = 4; rvec[0] = 24; rvec[1] = 32; rvec[2] = 32; rvec[3] = 16; break;
		default :
			*nradices = 7;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 832 :						/* 832K */
		switch(radix_set) {
		case 0 :
			numrad = 4; rvec[0] =208; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; break;
		case 1 :
			numrad = 4; rvec[0] = 52; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; break;
		case 2 :
			numrad = 4; rvec[0] = 52; rvec[1] = 32; rvec[2] = 16; rvec[3] = 16; break;
		case 3 :
			numrad = 5; rvec[0] = 52; rvec[1] =  8; rvec[2] =  8; rvec[3] =  8; rvec[4] = 16; break;
		default :
			*nradices = 4;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 896 :						/* 896K */
		switch(radix_set) {
		case 0 :
			numrad = 4; rvec[0] =224; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; break;
		case 1 :
			numrad = 4; rvec[0] = 56; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; break;
		case 2 :
			numrad = 4; rvec[0] = 56; rvec[1] = 32; rvec[2] = 16; rvec[3] = 16; break;
		case 3 :
			numrad = 5; rvec[0] = 56; rvec[1] =  8; rvec[2] =  8; rvec[3] =  8; rvec[4] = 16; break;
		case 4 :
			numrad = 4; rvec[0] = 28; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; break;
		case 5 :
			numrad = 4; rvec[0] = 28; rvec[1] = 32; rvec[2] = 32; rvec[3] = 16; break;
		default :
			*nradices = 6;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 960 :						/* 960K */
		switch(radix_set) {
		case 0 :
			numrad = 3; rvec[0] =960; rvec[1] = 16; rvec[2] = 32; break;
		case 1 :
			numrad = 4; rvec[0] =240; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; break;
		case 2 :
			numrad = 4; rvec[0] = 60; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; break;
		case 3 :
			numrad = 4; rvec[0] = 60; rvec[1] = 32; rvec[2] = 16; rvec[3] = 16; break;
		case 4 :
			numrad = 5; rvec[0] = 60; rvec[1] =  8; rvec[2] =  8; rvec[3] =  8; rvec[4] = 16; break;
		default :
			*nradices = 5;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 992 :
		switch(radix_set) {
		case 0 :
			numrad = 3; rvec[0] =992; rvec[1] = 16; rvec[2] = 32; break;
		case 1 :
			numrad = 3; rvec[0] =992; rvec[1] = 32; rvec[2] = 16; break;
		default :
			*nradices = 2;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 1008 :
		switch(radix_set) {
		case 0 :
			numrad = 3; rvec[0]=1008; rvec[1] = 16; rvec[2] = 32; break;
		case 1 :
			numrad = 3; rvec[0]=1008; rvec[1] = 32; rvec[2] = 16; break;
		case 2 :
			numrad = 4; rvec[0] = 63; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; break;
		default :
			*nradices = 3;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 1024 :					/* 1M = 1024K */
		switch(radix_set) {
		case 0 :
			numrad = 3; rvec[0]=1024; rvec[1] = 16; rvec[2] = 32; break;
		case 1 :
			numrad = 3; rvec[0]=1024; rvec[1] = 32; rvec[2] = 16; break;
		case 2 :
			numrad = 4; rvec[0] =256; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; break;
		case 3 :
			numrad = 4; rvec[0] =128; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; break;
		case 4 :
			numrad = 4; rvec[0] = 64; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; break;
		case 5 :
			numrad = 4; rvec[0] = 64; rvec[1] = 32; rvec[2] = 16; rvec[3] = 16; break;
		case 6 :
			numrad = 5; rvec[0] = 64; rvec[1] =  8; rvec[2] =  8; rvec[3] =  8; rvec[4] = 16; break;
		case 7 :
			numrad = 4; rvec[0] = 32; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; break;
		case 8 :
			numrad = 4; rvec[0] = 32; rvec[1] = 32; rvec[2] = 32; rvec[3] = 16; break;
		case 9 :
			numrad = 5; rvec[0] = 32; rvec[1] =  8; rvec[2] =  8; rvec[3] = 16; rvec[4] = 16; break;
		case 10:
			numrad = 4; rvec[0] = 16; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; break;
		case 11:
			numrad = 5; rvec[0] = 16; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
		default :
			*nradices = 12;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 1152 :					/* 1.125M = 1152K */
		switch(radix_set) {
		case 0 :
			numrad = 4; rvec[0] =288; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; break;
		case 1 :
			numrad = 4; rvec[0] =144; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; break;
		case 2 :
			numrad = 4; rvec[0] = 36; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; break;
		case 3 :
			numrad = 4; rvec[0] = 36; rvec[1] = 32; rvec[2] = 32; rvec[3] = 16; break;
		default :
			*nradices = 4;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 1280 :					/* 1.25M = 1280K */
		switch(radix_set) {
		case 0 :
			numrad = 4; rvec[0] =320; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; break;
		case 1 :
			numrad = 4; rvec[0] =160; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; break;
		case 2 :
			numrad = 4; rvec[0] = 40; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; break;
		case 3 :
			numrad = 4; rvec[0] = 40; rvec[1] = 32; rvec[2] = 32; rvec[3] = 16; break;
		case 4 :
			numrad = 4; rvec[0] = 20; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; break;
		case 5 :
			numrad = 5; rvec[0] = 20; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
		default :
			*nradices = 6;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 1408 :					/* 1.375M = 1408K */
		switch(radix_set) {
		case 0 :
			numrad = 4; rvec[0] =352; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; break;
		case 1 :
			numrad = 4; rvec[0] =176; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; break;
		case 2 :
			numrad = 4; rvec[0] = 44; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; break;
		case 3 :
			numrad = 4; rvec[0] = 44; rvec[1] = 32; rvec[2] = 32; rvec[3] = 16; break;
		default :
			*nradices = 4;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 1536 :					/* 1.5M = 1536K */
		switch(radix_set) {
		case 0 :
			numrad = 3; rvec[0] =768; rvec[1] = 32; rvec[2] = 32; break;
		case 1 :
			numrad = 4; rvec[0] =192; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; break;
		case 2 :
			numrad = 4; rvec[0] = 48; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; break;
		case 3 :
			numrad = 4; rvec[0] = 48; rvec[1] = 32; rvec[2] = 32; rvec[3] = 16; break;
		case 4 :
			numrad = 4; rvec[0] = 24; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; break;
		default :
			*nradices = 5;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 1664 :					/* 1.625M = 1664K */
		switch(radix_set) {
		case 0 :
			numrad = 4; rvec[0] =208; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; break;
		case 1 :
			numrad = 4; rvec[0] = 52; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; break;
		case 2 :
			numrad = 4; rvec[0] = 52; rvec[1] = 32; rvec[2] = 32; rvec[3] = 16; break;
		default :
			*nradices = 3;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 1792 :					/* 1.75M = 1792K */
		switch(radix_set) {
		case 0 :
			numrad = 4; rvec[0] =224; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; break;
		case 1 :
			numrad = 4; rvec[0] = 56; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; break;
		case 2 :
			numrad = 4; rvec[0] = 56; rvec[1] = 32; rvec[2] = 32; rvec[3] = 16; break;
		case 3 :
			numrad = 4; rvec[0] = 28; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; break;
		case 4 :
			numrad = 5; rvec[0] = 28; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
		default :
			*nradices = 5;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 1920 :					/* 1.875M = 1920K */
		switch(radix_set) {
		case 0 :
			numrad = 3; rvec[0] =960; rvec[1] = 32; rvec[2] = 32; break;
		case 1 :
			numrad = 4; rvec[0] =240; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; break;
		case 2 :
			numrad = 4; rvec[0] = 60; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; break;
		case 3 :
			numrad = 4; rvec[0] = 60; rvec[1] = 32; rvec[2] = 32; rvec[3] = 16; break;
		default :
			*nradices = 4;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 1984 :
		switch(radix_set) {
		case 0 :
			numrad = 3; rvec[0] =992; rvec[1] = 32; rvec[2] = 32; break;
		default :
			*nradices = 1;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 2016 :
		switch(radix_set) {
		case 0 :
			numrad = 3; rvec[0]=1008; rvec[1] = 32; rvec[2] = 32; break;
		case 1 :
			numrad = 4; rvec[0] = 63; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; break;
		case 2 :
			numrad = 4; rvec[0] = 63; rvec[1] = 32; rvec[2] = 32; rvec[3] = 16; break;
		default :
			*nradices = 3;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 2048 :					/* 2M = 2048K */
		switch(radix_set) {
		case 0 :
			numrad = 3; rvec[0]=1024; rvec[1] = 32; rvec[2] = 32; break;
		case 1 :
			numrad = 4; rvec[0] =256; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; break;
		case 2 :
			numrad = 4; rvec[0] =128; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; break;
		case 3 :
			numrad = 4; rvec[0] =128; rvec[1] = 32; rvec[2] = 16; rvec[3] = 16; break;
		case 4 :
			numrad = 4; rvec[0] = 64; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; break;
		case 5 :
			numrad = 4; rvec[0] = 64; rvec[1] = 32; rvec[2] = 32; rvec[3] = 16; break;
		case 6 :
			numrad = 4; rvec[0] = 32; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; break;
		case 7 :
			numrad = 5; rvec[0] = 32; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
		case 8 :
			numrad = 5; rvec[0] = 16; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
		case 9 :
			numrad = 6; rvec[0] = 16; rvec[1] =  8; rvec[2] =  8; rvec[3] =  8; rvec[4] =  8; rvec[5] = 16; break;
		default :
			*nradices = 10;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 2304 :					/* 2.25M = 2304K */
		switch(radix_set) {
		case 0 :
			numrad = 4; rvec[0] =288; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; break;
		case 1 :
			numrad = 4; rvec[0] =144; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; break;
		case 2 :
			numrad = 4; rvec[0] =144; rvec[1] = 32; rvec[2] = 16; rvec[3] = 16; break;
		case 3 :
			numrad = 5; rvec[0] =144; rvec[1] =  8; rvec[2] =  8; rvec[3] =  8; rvec[4] = 16; break;
		case 4 :
			numrad = 4; rvec[0] = 36; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; break;
		default :
			*nradices = 5;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 2560 :					/* 2.5M = 2560K */
		switch(radix_set) {
		case 0 :
			numrad = 4; rvec[0] =320; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; break;
		case 1 :
			numrad = 4; rvec[0] =160; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; break;
		case 2 :
			numrad = 4; rvec[0] =160; rvec[1] = 32; rvec[2] = 16; rvec[3] = 16; break;
		case 3 :
			numrad = 5; rvec[0] =160; rvec[1] =  8; rvec[2] =  8; rvec[3] =  8; rvec[4] = 16; break;
		case 4 :
			numrad = 4; rvec[0] = 40; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; break;
		case 5 :
			numrad = 5; rvec[0] = 20; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
		default :
			*nradices = 6;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 2816 :					/* 2.75M = 2816K */
		switch(radix_set) {
		case 0 :
			numrad = 4; rvec[0] =352; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; break;
		case 1 :
			numrad = 4; rvec[0] =176; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; break;
		case 2 :
			numrad = 4; rvec[0] =176; rvec[1] = 32; rvec[2] = 16; rvec[3] = 16; break;
		case 3 :
			numrad = 5; rvec[0] =176; rvec[1] =  8; rvec[2] =  8; rvec[3] =  8; rvec[4] = 16; break;
		case 4 :
			numrad = 4; rvec[0] = 44; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; break;
		case 5 :
			numrad = 5; rvec[0] = 44; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
		default :
			*nradices = 6;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 3072 :					/* 3M = 3072K */
		switch(radix_set) {
		case 0 :
			numrad = 4; rvec[0] =768; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; break;
		case 1 :
			numrad = 4; rvec[0] =192; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; break;
		case 2 :
			numrad = 4; rvec[0] =192; rvec[1] = 32; rvec[2] = 16; rvec[3] = 16; break;
		case 3 :
			numrad = 5; rvec[0] =192; rvec[1] =  8; rvec[2] =  8; rvec[3] =  8; rvec[4] = 16; break;
		case 4 :
			numrad = 4; rvec[0] = 48; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; break;
		case 5 :
			numrad = 5; rvec[0] = 24; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
		default :
			*nradices = 6;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 3328 :					/* 3.25M = 3328K */
		switch(radix_set) {
		case 0 :
			numrad = 4; rvec[0] =208; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; break;
		case 1 :
			numrad = 4; rvec[0] =208; rvec[1] = 32; rvec[2] = 16; rvec[3] = 16; break;
		case 2 :
			numrad = 5; rvec[0] =208; rvec[1] =  8; rvec[2] =  8; rvec[3] =  8; rvec[4] = 16; break;
		case 3 :
			numrad = 4; rvec[0] = 52; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; break;
		default :
			*nradices = 4;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 3584 :					/* 3.5M = 3584K */
		switch(radix_set) {
		case 0 :
			numrad = 4; rvec[0] =224; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; break;
		case 1 :
			numrad = 4; rvec[0] =224; rvec[1] = 32; rvec[2] = 16; rvec[3] = 16; break;
		case 2 :
			numrad = 5; rvec[0] =224; rvec[1] =  8; rvec[2] =  8; rvec[3] =  8; rvec[4] = 16; break;
		case 3 :
			numrad = 4; rvec[0] = 56; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; break;
		case 4 :
			numrad = 5; rvec[0] = 28; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
		default :
			*nradices = 5;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 3840 :					/* 3.75M = 3840K */
		switch(radix_set) {
		case 0 :
			numrad = 4; rvec[0]= 960; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; break;
		case 1 :
			numrad = 4; rvec[0] =240; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; break;
		case 2 :
			numrad = 4; rvec[0] =240; rvec[1] = 32; rvec[2] = 16; rvec[3] = 16; break;
		case 3 :
			numrad = 5; rvec[0] =240; rvec[1] =  8; rvec[2] =  8; rvec[3] =  8; rvec[4] = 16; break;
		case 4 :
			numrad = 4; rvec[0] = 60; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; break;
		default :
			*nradices = 5;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 3968 :
		switch(radix_set) {
		case 0 :
			numrad = 4; rvec[0]= 992; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; break;
		default :
			*nradices = 1;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 4032 :
		switch(radix_set) {
		case 0 :
			numrad = 3; rvec[0]=4032; rvec[1] = 16; rvec[2] = 32; break;
		case 1 :
			numrad = 3; rvec[0]=4032; rvec[1] = 32; rvec[2] = 16; break;
		case 2 :
			numrad = 4; rvec[0]=1008; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; break;
		case 3 :
			numrad = 4; rvec[0] = 63; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; break;
		default :
			*nradices = 4;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 4096 :						/* 4M = 4096K */
		switch(radix_set) {
		case 0 :
			numrad = 4; rvec[0]=1024; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; break;
		case 1 :
			numrad = 4; rvec[0] =256; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; break;
		case 2 :
			numrad = 4; rvec[0] =256; rvec[1] = 32; rvec[2] = 16; rvec[3] = 16; break;
		case 3 :
			numrad = 5; rvec[0] =256; rvec[1] =  8; rvec[2] =  8; rvec[3] =  8; rvec[4] = 16; break;
		case 4 :
			numrad = 4; rvec[0] =128; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; break;
		case 5 :
			numrad = 4; rvec[0] =128; rvec[1] = 32; rvec[2] = 32; rvec[3] = 16; break;
		case 6 :
			numrad = 5; rvec[0] =128; rvec[1] =  8; rvec[2] =  8; rvec[3] = 16; rvec[4] = 16; break;
		case 7 :
			numrad = 4; rvec[0] = 64; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; break;
		case 8 :
			numrad = 5; rvec[0] = 64; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
		case 9 :
			numrad = 5; rvec[0] = 32; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
		case 10:
			numrad = 5; rvec[0] = 16; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; break;
		default :
			*nradices = 11;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 4608 :					/* 4.5M = 4608K */
		switch(radix_set) {
		case 0 :
			numrad = 4; rvec[0] =288; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; break;
		case 1 :
			numrad = 4; rvec[0] =288; rvec[1] = 32; rvec[2] = 16; rvec[3] = 16; break;
		case 2 :
			numrad = 5; rvec[0] =288; rvec[1] =  8; rvec[2] =  8; rvec[3] =  8; rvec[4] = 16; break;
		case 3 :
			numrad = 4; rvec[0] =144; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; break;
		case 4 :
			numrad = 4; rvec[0] =144; rvec[1] = 32; rvec[2] = 32; rvec[3] = 16; break;
		case 5 :
			numrad = 5; rvec[0] =144; rvec[1] =  8; rvec[2] =  8; rvec[3] = 16; rvec[4] = 16; break;
		case 6 :
			numrad = 5; rvec[0] = 36; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
		default :
			*nradices = 7;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 5120 :					/* 5M = 5120K */
		switch(radix_set) {
		case 0 :
			numrad = 4; rvec[0] =320; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; break;
		case 1 :
			numrad = 4; rvec[0] =320; rvec[1] = 32; rvec[2] = 16; rvec[3] = 16; break;
		case 2 :
			numrad = 4; rvec[0] =160; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; break;
		case 3 :
			numrad = 4; rvec[0] =160; rvec[1] = 32; rvec[2] = 32; rvec[3] = 16; break;
		case 4 :
			numrad = 5; rvec[0] =160; rvec[1] =  8; rvec[2] =  8; rvec[3] = 16; rvec[4] = 16; break;
		case 5 :
			numrad = 5; rvec[0] = 40; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
		case 6 :
			numrad = 5; rvec[0] = 20; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; break;
		case 7 :
			numrad = 5; rvec[0] = 20; rvec[1] = 32; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
		default :
			*nradices = 8;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 5632 :					/* 5.5M = 5632K */
		switch(radix_set) {
		case 0 :
			numrad = 4; rvec[0] =352; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; break;
		case 1 :
			numrad = 4; rvec[0] =352; rvec[1] = 32; rvec[2] = 16; rvec[3] = 16; break;
		case 2 :
			numrad = 4; rvec[0] =176; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; break;
		case 3 :
			numrad = 4; rvec[0] =176; rvec[1] = 32; rvec[2] = 32; rvec[3] = 16; break;
		case 4 :
			numrad = 5; rvec[0] =176; rvec[1] =  8; rvec[2] =  8; rvec[3] = 16; rvec[4] = 16; break;
		case 5 :
			numrad = 5; rvec[0] = 44; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
		default :
			*nradices = 6;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 6144 :					/* 6M = 6144K */
		switch(radix_set) {
		case 0 :
			numrad = 4; rvec[0] =768; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; break;
		case 1 :
			numrad = 4; rvec[0] =192; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; break;
		case 2 :
			numrad = 4; rvec[0] =192; rvec[1] = 32; rvec[2] = 32; rvec[3] = 16; break;
		case 3 :
			numrad = 5; rvec[0] =192; rvec[1] =  8; rvec[2] =  8; rvec[3] = 16; rvec[4] = 16; break;
		case 4 :
			numrad = 5; rvec[0] = 48; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
		case 5 :
			numrad = 5; rvec[0] = 24; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; break;
		case 6 :
			numrad = 5; rvec[0] = 24; rvec[1] = 32; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
		default :
			*nradices = 7;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 6656 :					/* 6.5M = 6656K */
		switch(radix_set) {
		case 0 :
			numrad = 4; rvec[0] =208; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; break;
		case 1 :
			numrad = 4; rvec[0] =208; rvec[1] = 32; rvec[2] = 32; rvec[3] = 16; break;
		case 2 :
			numrad = 5; rvec[0] =208; rvec[1] =  8; rvec[2] =  8; rvec[3] = 16; rvec[4] = 16; break;
		case 3 :
			numrad = 5; rvec[0] = 52; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
		default :
			*nradices = 4;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 7168 :					/* 7M = 7168K */
		switch(radix_set) {
		case 0 :
			numrad = 4; rvec[0] =224; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; break;
		case 1 :
			numrad = 4; rvec[0] =224; rvec[1] = 32; rvec[2] = 32; rvec[3] = 16; break;
		case 2 :
			numrad = 5; rvec[0] =224; rvec[1] =  8; rvec[2] =  8; rvec[3] = 16; rvec[4] = 16; break;
		case 3 :
			numrad = 5; rvec[0] = 56; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
		case 4 :
			numrad = 5; rvec[0] = 28; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; break;
		case 5 :
			numrad = 5; rvec[0] = 28; rvec[1] = 32; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
		default :
			*nradices = 6;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 7680 :					/* 7.5M = 7680K */
		switch(radix_set) {
		case 0 :
			numrad = 4; rvec[0]= 960; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; break;
		case 1 :
			numrad = 4; rvec[0] =240; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; break;
		case 2 :
			numrad = 4; rvec[0] =240; rvec[1] = 32; rvec[2] = 32; rvec[3] = 16; break;
		case 3 :
			numrad = 5; rvec[0] =240; rvec[1] =  8; rvec[2] =  8; rvec[3] = 16; rvec[4] = 16; break;
		case 4 :
			numrad = 5; rvec[0] = 60; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
		default :
			*nradices = 5;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 7936 :
		switch(radix_set) {
		case 0 :
			numrad = 4; rvec[0]= 992; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; break;
		default :
			*nradices = 1;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 8064 :
		switch(radix_set) {
		case 0 :
			numrad = 3; rvec[0]=4032; rvec[1] = 32; rvec[2] = 32; break;
		case 1 :
			numrad = 4; rvec[0]=1008; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; break;
		case 2 :
			numrad = 5; rvec[0] = 63; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
		default :
			*nradices = 3;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 8192 :						/* 8M = 8192K */
		switch(radix_set) {
		case 0 :
			numrad = 4; rvec[0]=1024; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; break;
		case 1 :
			numrad = 4; rvec[0] =256; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; break;
		case 2 :
			numrad = 4; rvec[0] =256; rvec[1] = 32; rvec[2] = 32; rvec[3] = 16; break;
		case 3 :
			numrad = 5; rvec[0] =256; rvec[1] =  8; rvec[2] =  8; rvec[3] = 16; rvec[4] = 16; break;
		case 4 :
			numrad = 4; rvec[0] =128; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; break;
		case 5 :
			numrad = 5; rvec[0] = 64; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
		case 6 :
			numrad = 5; rvec[0] = 32; rvec[1] = 32; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
		case 7 :
			numrad = 5; rvec[0] = 16; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; break;
		default :
			*nradices = 8;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 9216 :					/* 9M = 9216K */
		switch(radix_set) {
		case 0 :
			numrad = 4; rvec[0] =288; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; break;
		case 1 :
			numrad = 4; rvec[0] =288; rvec[1] = 32; rvec[2] = 32; rvec[3] = 16; break;
		case 2 :
			numrad = 5; rvec[0] =288; rvec[1] =  8; rvec[2] =  8; rvec[3] = 16; rvec[4] = 16; break;
		case 3 :
			numrad = 4; rvec[0] =144; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; break;
		case 4 :
			numrad = 5; rvec[0] =144; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
		case 5 :
			numrad = 5; rvec[0] = 36; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; break;
		case 6 :
			numrad = 5; rvec[0] = 36; rvec[1] = 32; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
		default :
			*nradices = 7;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 10240 :				/* 10M = 10240K */
		switch(radix_set) {
		case 0 :
			numrad = 4; rvec[0] =320; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; break;
		case 1 :
			numrad = 4; rvec[0] =320; rvec[1] = 32; rvec[2] = 32; rvec[3] = 16; break;
		case 2 :
			numrad = 4; rvec[0] =160; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; break;
		case 3 :
			numrad = 5; rvec[0] =160; rvec[1] = 16; rvec[2] =  8; rvec[3] = 16; rvec[4] = 16; break;
		case 4 :
			numrad = 5; rvec[0] = 40; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; break;
		case 5 :
			numrad = 5; rvec[0] = 40; rvec[1] = 32; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
		case 6 :
			numrad = 5; rvec[0] = 20; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; break;
		case 7 :
			numrad = 5; rvec[0] = 20; rvec[1] = 32; rvec[2] = 32; rvec[3] = 16; rvec[4] = 16; break;
		default :
			*nradices = 8;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 11264 :				/* 11M = 11264K */
		switch(radix_set) {
		case 0 :
			numrad = 4; rvec[0] =352; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; break;
		case 1 :
			numrad = 4; rvec[0] =352; rvec[1] = 32; rvec[2] = 32; rvec[3] = 16; break;
		case 2 :
			numrad = 4; rvec[0] =176; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; break;
		case 3 :
			numrad = 5; rvec[0] = 44; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; break;
		case 4 :
			numrad = 5; rvec[0] = 44; rvec[1] = 32; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
		default :
			*nradices = 5;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 12288 :				/* 12M = 12288K */
		switch(radix_set) {
		case 0 :
			numrad = 4; rvec[0] =768; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; break;
		case 1 :
			numrad = 4; rvec[0] =768; rvec[1] = 32; rvec[2] = 16; rvec[3] = 16; break;
		case 2 :
			numrad = 4; rvec[0] =192; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; break;
		case 3 :
			numrad = 5; rvec[0] = 48; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; break;
		case 4 :
			numrad = 5; rvec[0] = 48; rvec[1] = 32; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
		case 5 :
			numrad = 5; rvec[0] = 24; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; break;
		case 6 :
			numrad = 5; rvec[0] = 24; rvec[1] = 32; rvec[2] = 32; rvec[3] = 16; rvec[4] = 16; break;
		default :
			*nradices = 7;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 13312 :				/* 13M = 13312K */
		switch(radix_set) {
		case 0 :
			numrad = 4; rvec[0] =208; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; break;
		case 1 :
			numrad = 5; rvec[0] = 52; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; break;
		case 2 :
			numrad = 5; rvec[0] = 52; rvec[1] = 32; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
		default :
			*nradices = 3;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 14336 :				/* 14M = 14336K */
		switch(radix_set) {
		case 0 :
			numrad = 4; rvec[0] =224; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; break;
		case 1 :
			numrad = 5; rvec[0] = 56; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; break;
		case 2 :
			numrad = 5; rvec[0] = 28; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; break;
		default :
			*nradices = 3;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 15360 :				/* 15M = 15360K */
		switch(radix_set) {
		case 0 :
			numrad = 4; rvec[0]= 960; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; break;
		case 1 :
			numrad = 4; rvec[0] =240; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; break;
		case 2 :
			numrad = 5; rvec[0] = 60; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; break;
		default :
			*nradices = 3;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 15872 :
		switch(radix_set) {
		case 0 :
			numrad = 4; rvec[0]= 992; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; break;
		default :
			*nradices = 1;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 16128 :
		switch(radix_set) {
		case 0 :
			numrad = 4; rvec[0]=4032; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; break;
		case 1 :
			numrad = 4; rvec[0]=1008; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; break;
		case 2 :
			numrad = 5; rvec[0] = 63; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; break;
		default :
			*nradices = 3;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 16384 :				/* 16M = 16384K */
		switch(radix_set) {
		case 0 :
			numrad = 4; rvec[0]=1024; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; break;
		case 1 :
			numrad = 4; rvec[0] =256; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; break;
		case 2 :
			numrad = 5; rvec[0] =256; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
		case 3 :
			numrad = 5; rvec[0] =128; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
		case 4 :
			numrad = 5; rvec[0] = 64; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; break;
		case 5 :
			numrad = 5; rvec[0] = 32; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; break;
		case 6 :
			numrad = 5; rvec[0] = 16; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		default :
			*nradices = 7;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 18432 :				/* 18M = 18432K */
		switch(radix_set) {
		case 0 :
			numrad = 4; rvec[0] =288; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; break;
		case 1 :
			numrad = 5; rvec[0] =144; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
		case 2 :
			numrad = 5; rvec[0] = 36; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; break;
		default :
			*nradices = 3;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 20480 :				/* 20M = 20480K */
		switch(radix_set) {
		case 0 :
			numrad = 4; rvec[0] =320; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; break;
		case 1 :
			numrad = 5; rvec[0] =160; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
		case 2 :
			numrad = 5; rvec[0] = 40; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; break;
		case 3 :
			numrad = 5; rvec[0] = 20; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		default :
			*nradices = 4;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 22528 :				/* 22M = 22528K */
		switch(radix_set) {
		case 0 :
			numrad = 4; rvec[0] =352; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; break;
		case 1 :
			numrad = 5; rvec[0] =176; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
		case 2 :
			numrad = 5; rvec[0] = 44; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; break;
		default :
			*nradices = 3;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 24576 :				/* 24M = 24576K */
		switch(radix_set) {
		case 0 :
			numrad = 4; rvec[0] =768; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; break;
		case 1 :
			numrad = 5; rvec[0] =192; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
		case 2 :
			numrad = 5; rvec[0] = 48; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; break;
		case 3 :
			numrad = 5; rvec[0] = 24; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		default :
			*nradices = 4;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 26624 :				/* 26M = 26624K */
		switch(radix_set) {
		case 0 :
			numrad = 5; rvec[0] =208; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
		case 1 :
			numrad = 5; rvec[0] = 52; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; break;
		default :
			*nradices = 2;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 28672 :				/* 28M = 28672K */
		switch(radix_set) {
		case 0 :
			numrad = 5; rvec[0] =224; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
		case 1 :
			numrad = 5; rvec[0] = 56; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; break;
		case 2 :
			numrad = 5; rvec[0] = 28; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		default :
			*nradices = 3;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 30720 :				/* 30M = 30720K */
		switch(radix_set) {
		case 0 :
			numrad = 4; rvec[0]= 960; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; break;
		case 1 :
			numrad = 5; rvec[0] =240; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
		case 2 :
			numrad = 5; rvec[0] = 60; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; break;
		default :
			*nradices = 3;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 31744 :
		switch(radix_set) {
		case 0 :
			numrad = 4; rvec[0]= 992; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; break;
		case 1 :
			numrad = 4; rvec[0]= 992; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; break;
		default :
			*nradices = 2;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 32256 :
		switch(radix_set) {
		case 0 :
			numrad = 4; rvec[0]=4032; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; break;
		case 1 :
			numrad = 4; rvec[0]=1008; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; break;
		case 2 :
			numrad = 5; rvec[0] = 63; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; break;
		default :
			*nradices = 3;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 32768 :				/* 32M = 32768K */
		switch(radix_set) {
		case 0 :
			numrad = 4; rvec[0]=1024; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; break;
		case 1 :
			numrad = 5; rvec[0] =256; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
		case 2 :
			numrad = 5; rvec[0] =128; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; break;
		case 3 :
			numrad = 5; rvec[0] = 64; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; break;
		case 4 :
			numrad = 5; rvec[0] = 32; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		case 5 :
			numrad = 5; rvec[0] = 32; rvec[1] = 32; rvec[2] = 32; rvec[3] = 16; rvec[4] = 32; break;
		case 6 :
			numrad = 5; rvec[0] = 16; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		default :
			*nradices = 7;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 36864 :				/* 36M = 36864K */
		switch(radix_set) {
		case 0 :
			numrad = 5; rvec[0] =288; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
		case 1 :
			numrad = 5; rvec[0] =144; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; break;
		case 2 :
			numrad = 5; rvec[0] = 36; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		default :
			*nradices = 3;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 40960 :				/* 40M = 40960K */
		switch(radix_set) {
		case 0 :
			numrad = 5; rvec[0] =320; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
		case 1 :
			numrad = 5; rvec[0] =160; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; break;
		case 2 :
			numrad = 5; rvec[0] = 40; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		case 3 :
			numrad = 5; rvec[0] = 20; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		case 4 :
			numrad = 6; rvec[0] = 20; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 16; break;
		default :
			*nradices = 5;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 45056 :				/* 44M = 45056K */
		switch(radix_set) {
		case 0 :
			numrad = 5; rvec[0] =352; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
		case 1 :
			numrad = 5; rvec[0] =176; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; break;
		case 2 :
			numrad = 5; rvec[0] = 44; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		default :
			*nradices = 3;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 49152 :				/* 48M = 49152K */
		switch(radix_set) {
		case 0 :
			numrad = 4; rvec[0] =768; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; break;
		case 1 :
			numrad = 5; rvec[0] =192; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; break;
		case 2 :
			numrad = 5; rvec[0] = 48; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		case 3 :
			numrad = 5; rvec[0] = 24; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		case 4 :
			numrad = 6; rvec[0] = 24; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 16; break;
		default :
			*nradices = 5;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 53248 :				/* 52M = 53248K */
		switch(radix_set) {
		case 0 :
			numrad = 5; rvec[0] =208; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; break;
		case 1 :
			numrad = 5; rvec[0] = 52; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		default :
			*nradices = 2;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 57344 :				/* 56M = 57344K */
		switch(radix_set) {
		case 0 :
			numrad = 5; rvec[0] =224; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; break;
		case 1 :
			numrad = 5; rvec[0] = 56; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		case 2 :
			numrad = 5; rvec[0] = 28; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		case 3 :
			numrad = 6; rvec[0] = 28; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 16; break;
		default :
			*nradices = 4;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 61440 :				/* 60M = 61440K */
		switch(radix_set) {
		case 0 :
			numrad = 4; rvec[0]= 960; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; break;
		case 1 :
			numrad = 5; rvec[0] =240; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; break;
		case 2 :
			numrad = 5; rvec[0] =240; rvec[1] = 32; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
		case 3 :
			numrad = 5; rvec[0] = 60; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		case 4 :
			numrad = 5; rvec[0] = 60; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; rvec[4] = 16; break;
		default :
			*nradices = 5;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 63488 :
		switch(radix_set) {
		case 0 :
			numrad = 4; rvec[0]= 992; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; break;
		default :
			*nradices = 1;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 64512 :
		switch(radix_set) {
		case 0 :
			numrad = 4; rvec[0]=4032; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; break;
		case 1 :
			numrad = 4; rvec[0]=1008; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; break;
		case 2 :
			numrad = 5; rvec[0] = 63; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		default :
			*nradices = 3;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 65536 :				/* 64M = 65536K */
		switch(radix_set) {
		case 0 :
			numrad = 4; rvec[0]=1024; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; break;
		case 1 :
			numrad = 5; rvec[0] =256; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; break;
		case 2 :
			numrad = 5; rvec[0] =256; rvec[1] = 32; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
		case 3 :
			numrad = 5; rvec[0] =128; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; break;
		case 4 :
			numrad = 5; rvec[0] =128; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; break;
		case 5 :
			numrad = 5; rvec[0] = 64; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		case 6 :
			numrad = 5; rvec[0] = 64; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		case 7 :
			numrad = 5; rvec[0] = 32; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		case 8 :
			numrad = 6; rvec[0] = 32; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 16; break;
		case 9 :
			numrad = 6; rvec[0] = 16; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 32; break;
		default :
			*nradices = 10;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 73728 :				/* 72M = 73728K */
		switch(radix_set) {
		case 0 :
			numrad = 5; rvec[0] =288; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; break;
		case 1 :
			numrad = 5; rvec[0] =144; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; break;
		case 2 :
			numrad = 5; rvec[0] = 36; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		case 3 :
			numrad = 6; rvec[0] = 36; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 16; break;
		default :
			*nradices = 4;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 81920 :				/* 80M = 81920K */
		switch(radix_set) {
		case 0 :
			numrad = 5; rvec[0] =320; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; break;
		case 1 :
			numrad = 5; rvec[0] =160; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; break;
		case 2 :
			numrad = 5; rvec[0] = 40; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		case 3 :
			numrad = 6; rvec[0] = 40; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 16; break;
		case 4 :
			numrad = 6; rvec[0] = 20; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 32; break;
		default :
			*nradices = 5;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 90112 :				/* 88M = 90112K */
		switch(radix_set) {
		case 0 :
			numrad = 5; rvec[0] =352; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; break;
		case 1 :
			numrad = 5; rvec[0] =176; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; break;
		case 2 :
			numrad = 5; rvec[0] = 44; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		case 3 :
			numrad = 6; rvec[0] = 44; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 16; break;
		default :
			*nradices = 4;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 98304 :				/* 96M = 98304K */
		switch(radix_set) {
		case 0 :
			numrad = 5; rvec[0] =768; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
		case 1 :
			numrad = 5; rvec[0] =192; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; break;
		case 2 :
			numrad = 5; rvec[0] = 48; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		case 3 :
			numrad = 6; rvec[0] = 48; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 16; break;
		case 4 :
			numrad = 6; rvec[0] = 24; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 32; break;
		default :
			*nradices = 5;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 106496 :				/* 104M = 106496K */
		switch(radix_set) {
		case 0 :
			numrad = 5; rvec[0] =208; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; break;
		case 1 :
			numrad = 5; rvec[0] = 52; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		case 2 :
			numrad = 6; rvec[0] = 52; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 16; break;
		default :
			*nradices = 3;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 114688 :				/* 112M = 114688K */
		switch(radix_set) {
		case 0 :
			numrad = 5; rvec[0] =224; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; break;
		case 1 :
			numrad = 5; rvec[0] = 56; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		case 2 :
			numrad = 6; rvec[0] = 56; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 16; break;
		case 3 :
			numrad = 6; rvec[0] = 28; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 32; break;
		default :
			*nradices = 4;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 122880 :				/* 120M = 122880K */
		switch(radix_set) {
		case 0 :
			numrad = 5; rvec[0] =960; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
		case 1 :
			numrad = 5; rvec[0] =240; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; break;
		case 2 :
			numrad = 5; rvec[0] = 60; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		default :
			*nradices = 3;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 129024 :
		switch(radix_set) {
		case 0 :
			numrad = 4; rvec[0]=4032; rvec[1] = 32; rvec[2] = 32; rvec[3] = 16; break;
		case 1 :
			numrad = 5; rvec[0]=1008; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
		case 2 :
			numrad = 5; rvec[0] = 63; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		default :
			*nradices = 3;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 131072 :				/* 128M = 131072K */
		switch(radix_set) {
		case 0 :
			numrad = 5; rvec[0]=1024; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
		case 1 :
			numrad = 5; rvec[0] =256; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; break;
		case 2 :
			numrad = 5; rvec[0] =128; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		case 3 :
			numrad = 5; rvec[0] = 64; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		case 4 :
			numrad = 6; rvec[0] = 32; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 32; break;
		case 5 :
			numrad = 6; rvec[0] = 32; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; rvec[5] = 16; break;
		case 6 :
			numrad = 6; rvec[0] = 16; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; rvec[5] = 32; break;
		default :
			*nradices = 7;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 147456 :				/* 144M = 147456K */
		switch(radix_set) {
		case 0 :
			numrad = 5; rvec[0] =288; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; break;
		case 1 :
			numrad = 5; rvec[0] =144; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		case 2 :
			numrad = 6; rvec[0] = 36; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 32; break;
		default :
			*nradices = 3;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 163840 :				/* 160M = 163840K */
		switch(radix_set) {
		case 0 :
			numrad = 5; rvec[0] =320; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; break;
		case 1 :
			numrad = 5; rvec[0] =160; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		case 2 :
			numrad = 6; rvec[0] = 40; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 32; break;
		case 3 :
			numrad = 6; rvec[0] = 20; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; rvec[5] = 32; break;
		default :
			*nradices = 4;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 180224 :				/* 176M = 180224K */
		switch(radix_set) {
		case 0 :
			numrad = 5; rvec[0] =352; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; break;
		case 1 :
			numrad = 5; rvec[0] =176; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		case 2 :
			numrad = 6; rvec[0] = 44; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 32; break;
		default :
			*nradices = 3;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 196608 :				/* 192M = 196608K */
		switch(radix_set) {
		case 0 :
			numrad = 5; rvec[0] =768; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; break;
		case 1 :
			numrad = 5; rvec[0] =192; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		case 2 :
			numrad = 6; rvec[0] = 48; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 32; break;
		case 3 :
			numrad = 6; rvec[0] = 24; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; rvec[5] = 32; break;
		default :
			*nradices = 4;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 212992 :				/* 208M = 212992K */
		switch(radix_set) {
		case 0 :
			numrad = 5; rvec[0] =208; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		case 1 :
			numrad = 6; rvec[0] = 52; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 32; break;
		default :
			*nradices = 2;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 229376 :				/* 224M = 229376K */
		switch(radix_set) {
		case 0 :
			numrad = 5; rvec[0] =224; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		case 1 :
			numrad = 6; rvec[0] = 56; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 32; break;
		case 2 :
			numrad = 6; rvec[0] = 28; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; rvec[5] = 32; break;
		default :
			*nradices = 3;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 245760 :				/* 240M = 245760K */
		switch(radix_set) {
		case 0 :
			numrad = 5; rvec[0] =960; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; break;
		case 1 :
			numrad = 5; rvec[0] =240; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		case 2 :
			numrad = 6; rvec[0] = 60; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 32; break;
		default :
			*nradices = 3;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 258048 :
		switch(radix_set) {
		case 0 :
			numrad = 4; rvec[0]=4032; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; break;
		case 1 :
			numrad = 5; rvec[0]=1008; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; break;
		case 2 :
			numrad = 6; rvec[0] = 63; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 32; break;
		default :
			*nradices = 3;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 262144 :				/* 256M = 262144K */
		switch(radix_set) {
		case 0 :
			numrad = 5; rvec[0]=1024; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; break;
		case 1 :
			numrad = 5; rvec[0] =256; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		case 2 :
			numrad = 5; rvec[0] =128; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		case 3 :
			numrad = 6; rvec[0] =128; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 16; break;
		case 4 :
			numrad = 6; rvec[0] = 64; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 32; break;
		case 5 :
			numrad = 6; rvec[0] = 32; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; rvec[5] = 32; break;
		case 6 :
			numrad = 6; rvec[0] = 32; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; rvec[5] = 16; break;
		case 7 :
			numrad = 6; rvec[0] = 16; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; rvec[5] = 32; break;
		default :
			*nradices = 8;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 294912 :				/* 288M = 294912K */
		switch(radix_set) {
		case 0 :
			numrad = 5; rvec[0] =288; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		case 1 :
			numrad = 5; rvec[0] =144; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		case 2 :
			numrad = 6; rvec[0] =144; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 16; break;
		case 3 :
			numrad = 6; rvec[0] = 36; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; rvec[5] = 32; break;
		default :
			*nradices = 4;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 327680 :				/* 320M = 32768K */
		switch(radix_set) {
		case 0 :
			numrad = 5; rvec[0] =320; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		case 1 :
			numrad = 5; rvec[0] =160; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		case 2 :
			numrad = 6; rvec[0] =160; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 16; break;
		case 3 :
			numrad = 6; rvec[0] = 40; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; rvec[5] = 32; break;
		case 4 :
			numrad = 6; rvec[0] = 20; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; rvec[5] = 32; break;
		default :
			*nradices = 5;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 360448 :				/* 352M = 360448K */
		switch(radix_set) {
		case 0 :
			numrad = 5; rvec[0] =352; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		case 1 :
			numrad = 5; rvec[0] =176; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		case 2 :
			numrad = 6; rvec[0] =176; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 16; break;
		case 3 :
			numrad = 6; rvec[0] = 44; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; rvec[5] = 32; break;
		default :
			*nradices = 4;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 393216 :				/* 384M = 393216K */
		switch(radix_set) {
		case 0 :
			numrad = 5; rvec[0] =768; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; break;
		case 1 :
			numrad = 5; rvec[0] =192; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		case 2 :
			numrad = 6; rvec[0] =192; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 16; break;
		case 3 :
			numrad = 6; rvec[0] = 48; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; rvec[5] = 32; break;
		case 4 :
			numrad = 6; rvec[0] = 24; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; rvec[5] = 32; break;
		default :
			*nradices = 5;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 425984 :				/* 416M = 425984K */
		switch(radix_set) {
		case 0 :
			numrad = 5; rvec[0] =208; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		case 1 :
			numrad = 6; rvec[0] =208; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 16; break;
		case 2 :
			numrad = 6; rvec[0] = 52; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; rvec[5] = 32; break;
		default :
			*nradices = 3;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 458752 :				/* 448M = 458752K */
		switch(radix_set) {
		case 0 :
			numrad = 5; rvec[0] =224; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		case 1 :
			numrad = 6; rvec[0] =224; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 16; break;
		case 2 :
			numrad = 6; rvec[0] = 56; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; rvec[5] = 32; break;
		case 3 :
			numrad = 6; rvec[0] = 28; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; rvec[5] = 32; break;
		default :
			*nradices = 4;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 491520 :				/* 480M = 491520K */
		switch(radix_set) {
		case 0 :
			numrad = 5; rvec[0] =960; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; break;
		case 1 :
			numrad = 5; rvec[0] =240; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		case 2 :
			numrad = 6; rvec[0] =240; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 16; break;
		case 3 :
			numrad = 6; rvec[0] = 60; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; rvec[5] = 32; break;
		default :
			*nradices = 4;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 516096 :				/* 504M = 516096K */
		switch(radix_set) {
		case 0 :
			numrad = 5; rvec[0]=4032; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
		case 1 :
			numrad = 5; rvec[0]=1008; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; break;
		default :
			*nradices = 2;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 524288 :				/* 512M = 524288K */
		switch(radix_set) {
		case 0 :
			numrad = 5; rvec[0]=1024; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; break;
		case 1 :
			numrad = 5; rvec[0] =256; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		case 2 :
			numrad = 6; rvec[0] =256; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 16; break;
		case 3 :
			numrad = 6; rvec[0] =128; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 32; break;
		case 4 :
			numrad = 6; rvec[0] = 64; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; rvec[5] = 32; break;
		case 5 :
			numrad = 6; rvec[0] = 32; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; rvec[5] = 32; break;
		default :
			*nradices = 6;	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	default :
		 return ERR_FFTLENGTH_ILLEGAL;
	}

	ASSERT(rvec[0] <= MAX_RADIX, "Leading radix exceeds value of MAX_RADIX set in Mdata.h file!");

	// Check that there are at least 2 radices:
	if(numrad < 2) {
		fprintf(stderr,"ERROR: get_fft_radices: Specified %d radices; this must be >= 2.\n",numrad);
		return ERR_FFTLENGTH_ILLEGAL;
	}

	/* If user provided a radix array, make sure they gave a valid dimension: */
	if(radix_vec)
		ASSERT(radix_vec_dim >=	numrad,"get_fft_radices: radix_vec_dim has illegal value!");

	/* Check that N/2 = {product of the radices}, and if valid nradices and radix_vec pointers supplied,
	copy radices to the latter and	numrad to the former: */
	rad_prod = 1;

	if(nradices)
		*nradices =	numrad;

	for(i = 0; i <	numrad; i++)
	{
		rad_prod *= rvec[i];
		if(radix_vec)
			radix_vec[i] = rvec[i];
	}
	/* Zero-pad: */
	for(i =	numrad; i < radix_vec_dim; i++)
	{
		if(radix_vec)
			radix_vec[i] = 0;
	}

	if(rad_prod != n/2)
	{
		fprintf(stderr,"N = %u, radix_set = %u : product of complex radices %u != (FFT length/2)\n", n, radix_set, rad_prod);
		ASSERT(0,"0");
	}

	return 0;
}

/*
Simple self-tester for get_fft_radices()
*/
void	test_fft_radixtables()
{
	uint32 kblocks;
	uint32 radix0;
	int i,radset;
	int retval;
	uint32 nradices;

	i = 0;	/* Shift count for leading radix, i.e. length being tried = radix0*(1024)*(2^i). */
	while(1)
	{
		radset = 0;
		for(radix0 = 9; radix0 <= 16; ++radix0)
		{
			kblocks = (radix0 << i);

			if(kblocks > MAX_FFT_LENGTH_IN_K)
				return;

			retval = get_fft_radices(kblocks, radset, &nradices, 0x0, 0);
			if(retval == 0)
			{
				++radset;
				continue;
			}
			else if(retval == ERR_RADIXSET_UNAVAILABLE)
			{
				ASSERT(radset != 0, "test_fft_radixtables: Should only see ERR_RADIXSET_UNAVAILABLE for nonzero radix set index!");
				break;
			}
			else if(retval == ERR_FFTLENGTH_ILLEGAL)
			{
				fprintf(stderr,"ERROR: illegal FFT length %u K in test_fft_radixtables self-test!\n",kblocks);
				ASSERT(0,"0");
			}
			else
			{
				fprintf(stderr,"ERROR: unknown return value %d in test_fft_radixtables self-test; i = %d, kblocks = %u, radset = %u.\n", retval, i, kblocks, radset);
				ASSERT(0,"0");
			}
		}
		++i;
	}
}

/*
!...Set vector length, based on number of bits (p) in numbers to be FFT-multiplied.
Returns: FFT length in units of kdoubles, i.e. raw FFT length = (return value << 10) doubles:
*/
uint32 get_default_fft_length(uint64 p)
{
	uint32 nradices;
	uint32 leadingRadixVec[N_LEADING_RADICES] = {8,9,10,11,12,13,14,15};
	uint32 i, twoK, fftLen;

	ASSERT(PMAX > PMIN,"get_default_fft_length: PMAX > PMIN");
	if(p < PMIN || p > PMAX)
	{
		fprintf(stderr,"get_default_fft_length: invalid value for exponent %" PRIu64 "\n",p);
		ASSERT(0,"0");
		return 0;
	}

	/* Starting with N = 1K, Loop over all FFT lengths of form {8,9,10,11,12,13,14,15}*2^m,
	and return the smallest one for which maxP >= p: */
	i = 0;
	ASSERT(1024%leadingRadixVec[i] == 0,"get_default_fft_length: 1024%leadingRadixVec[0] == 0");
	twoK = 1024/leadingRadixVec[i];
	fftLen = leadingRadixVec[i]*twoK;
	for(;;)
	{
		if((fftLen >> 10) > MAX_FFT_LENGTH_IN_K)
			break;

		/* Only consider FFT lengths that are integer multiples of 1K: */
		if((fftLen & 1023) != 0)
			goto CYCLE;

		/* Make sure the FFT length is supported: */
		if(get_fft_radices((fftLen >> 10), 0, &nradices, 0x0, 0) != 0)
			goto CYCLE;

		/* As soon as the maxP >= p, we have our desired FFT length: */
		if(given_N_get_maxP(fftLen) >= p)
			return (fftLen >> 10);

	CYCLE:
		i = (i+1) % N_LEADING_RADICES;
		if(i == 0)
			twoK *= 2;
		fftLen = leadingRadixVec[i]*twoK;
	}
	if((fftLen >> 10) == 589824)
		fprintf(stderr,"get_default_fft_length: Allowing fftLen 576M just for informational purposes ... note this length is not supported.\n");
	else
		ASSERT(0,"get_default_fft_length: fftLen > MAX_FFT_LENGTH_IN_K!");
	return 0;
}

/*
!...Given a legal FFT length, return the next-larger supported length.
Returns: raw FFT length in units of doubles (contrast with get_default_fft_length(), which returns Kdoubles):
*/
uint32 get_nextlarger_fft_length(uint32 n)
{
	uint32 fftLen, lead4, rem2;	/* Leading 4 bits of input FFT length and number of bits in the remaining power of 2 */

	if(get_fft_radices((n >> 10), 0, 0x0, 0x0, 0) != 0)
	{
		sprintf(cbuf, "get_nextlarger_fft_length: Illegal or Unsupported input FFT length %u\n", n);
		ASSERT(0, cbuf);
	}

	/* Extract leading 4 bits of input FFT lengths, thus decomposing it into the form {8,9,10,11,12,13,14,15}*2^m,
	(m stored in rem2 variable here, decompose as n = lead4*2^rem2) then increment the leading portion -
	don't need to maintain leading part normalized in [8,15] here:
	*/
	rem2 = 32 - leadz32(n) - 4;
	lead4 = n >> rem2;
	ASSERT(lead4 > 7 && lead4 < 16,"get_nextlarger_fft_length: leading 4 bits of input FFT length out of range!");

	/* Make sure next-larger FFT length is supported: */
	++lead4;
	fftLen = lead4 << rem2;
	if(get_fft_radices((fftLen >> 10), 0, 0x0, 0x0, 0) != 0)
	{
		sprintf(cbuf, "get_nextlarger_fft_length: Next-larger FFT length %u not supported!\n", fftLen);
		DBG_WARN(HERE, cbuf, "", 0);
		return 0;
	}
	else
		return fftLen;
}

/*
For a given FFT length, estimate maximum exponent that can be tested.

This implements formula (8) in the F24 paper (Math Comp. 72 (243), pp.1555-1572,
December 2002) in order to estimate the maximum average wordsize for a given FFT length.
For roughly IEEE64-compliant arithmetic, an asymptotic constant of 0.6 (log2(C) in the
the paper, which recommends something around unity) seems to fit the observed data best.
*/
uint64 given_N_get_maxP(uint32 N)
{
	const double Bmant = 53;
#ifdef USE_FMADD
	const double AsympConst = 0.4;	// Allow slightly larger maxp if FMA used for floating-point arithmetic
#else
	const double AsympConst = 0.6;
#endif
	const double ln2inv = 1.0/log(2.0);
	double ln_N, lnln_N, l2_N, lnl2_N, l2l2_N, lnlnln_N, l2lnln_N;
	double Wbits, maxExp2;

	ln_N     = log(1.0*N);
	lnln_N   = log(ln_N);
	l2_N     = ln2inv*ln_N;
	lnl2_N   = log(l2_N);
	l2l2_N   = ln2inv*lnl2_N;
	lnlnln_N = log(lnln_N);
	l2lnln_N = ln2inv*lnlnln_N;

	Wbits = 0.5*( Bmant - AsympConst - 0.5*(l2_N + l2l2_N) - 1.5*(l2lnln_N) );
	maxExp2 = Wbits*N;
/*
	fprintf(stderr,"N = %8u K  maxP = %10u\n", N>>10, (uint32)maxExp2);
*/
	return (uint64)maxExp2;
}
/* Here a simple PARI 'script' to return maxExp for any desired N [FFT length in #doubles] and Bmant [#mantissa bits in float type]:
	maxp(Bmant, N, AsympConst) = { \
		ln2inv = 1.0/log(2.0); \
		ln_N = log(1.0*N); lnln_N = log(ln_N); l2_N = ln2inv*ln_N; lnl2_N = log(l2_N); l2l2_N = ln2inv*lnl2_N; lnlnln_N = log(lnln_N); l2lnln_N = ln2inv*lnlnln_N; \
		Wbits = 0.5*( Bmant - AsympConst - 0.5*(l2_N + l2l2_N) - 1.5*(l2lnln_N) ); \
		return(Wbits*N); \
	}
Ex: maxp(2240<<10, 0.4) = 43235170

And here a *nix bc function (invoke bc in floating-point mode, as 'bc -l'):
	define maxp(bmant, n, asympconst) {
		auto maxp, wbits, ln2inv, ln_n, lnln_n, l2_n, lnl2_n, l2l2_n, lnlnln_n, l2lnln_n;
		ln2inv = 1.0/l(2.0);
		ln_n = l(1.0*n); lnln_n = l(ln_n); l2_n = ln2inv*ln_n; lnl2_n = l(l2_n); l2l2_n = ln2inv*lnl2_n; lnlnln_n = l(lnln_n); l2lnln_n = ln2inv*lnlnln_n;
		wbits = 0.5*( bmant - asympconst - 0.5*(l2_n + l2l2_n) - 1.5*(l2lnln_n) );
		maxp = wbits*n;
		return maxp - maxp % 10^scale;
	}
Examples:
	Using DP-float: maxp(53,2048*2^10,0.4) = 39606917
					maxp(53,2048*2^10,0.6) = 39397201
	Using SP-float: maxp(24,2048*2^10,0.4) =  9198213
					maxp(24,2048*2^10,0.6) =  8988497
Here are maxp values for DP and SP-float per the above heuristic, both using the more-aggressive length-setting
given by using AsympConst = 0.6, for a range of FFT lengths covering current GIMPS DCs and testing wavefront (M = 2^20):
[code]FFT	DP maxp		SP maxp
	---	--------	--------
	 2M	 39397201	  8988497
	 3M	 58569855	 12956799
	 4M	 77597294	 16779886
	 5M	 96517023	 20495263
	 6M 115351074	 24124962
	...
	10M 190066770	 38023250
	20M				 70146273
	30M				100064264
	40M				128554502[/code]
I retested my several-years-ago SP-only hack of Mlucas side-by-side vs the DP version at 2M using a prime exponent ~9M - the SP version looks like it needs further debug of the DWT-weight/carry macro, so have not been able to confirm if the above SP numbers are realistic - but they're certain to be closer to realistic than my previous backwards-arguments numbers were.

With AsympConst = 0.6, this gives the following maxP values for various FFT lengths:
				maxn(N) as a function of AsympConst:
             N       AC = 0.6    AC = 0.4
        --------	----------	----------
             1 K		 22686		 22788	[0.45% larger]
             2 K		 44683		 44888
             3 K		 66435		 66742
             4 K		 88029		 88438
             5 K		109506		110018
             6 K		130892		131506
             7 K		152201		152918
             8 K		173445		174264
             9 K		194632		195554
            10 K		215769		216793
            12 K		257912		259141
            14 K		299904		301338
            16 K		341769		343407
            18 K		383521		385364
            20 K		425174		427222
            24 K		508222		510679
            28 K		590972		593840
            32 K		673470		676747
            36 K		755746		759433
            40 K		837827		841923
            48 K	   1001477	   1006392
            56 K	   1164540	   1170275
            64 K	   1327103	   1333656
            72 K	   1489228	   1496601
            80 K	   1650966	   1659158
            96 K	   1973430	   1983260
           112 K	   2294732	   2306201
           128 K	   2615043	   2628150
           144 K	   2934488	   2949234
           160 K	   3253166	   3269550
           176 K	   3571154	   3589176
           192 K	   3888516	   3908176
           208 K	   4205305	   4226604
           224 K	   4521565	   4544502
           240 K	   4837335	   4861911
           256 K	   5152648	   5178863
           288 K	   5782016	   5811507
           320 K	   6409862	   6442630
           352 K	   7036339	   7072384
           384 K	   7661575	   7700897
           416 K	   8285675	   8328273
           448 K	   8908726	   8954601
           480 K	   9530805	   9579957
           512 K	  10151977	  10204406
           576 K	  11391823	  11450805
           640 K	  12628648	  12694184
           704 K	  13862759	  13934849
           768 K	  15094405	  15173048
           832 K	  16323795	  16408992
           896 K	  17551103	  17642854
           960 K	  18776481	  18874785
          1024 K	  20000058	  20104916	[0.52% larger]
          1152 K	  22442252	  22560217
          1280 K	  24878447	  25009519
          1408 K	  27309250	  27453429
          1536 K	  29735157	  29892444
          1664 K	  32156582	  32326975
          1792 K	  34573872	  34757372
          1920 K	  36987325	  37183933
          2048 K	  39397201	  39606917
          2304 K	  44207097	  44443027
          2560 K	  49005071	  49267215
          2816 K	  53792328	  54080687
          3072 K	  58569855	  58884428
          3328 K	  63338470	  63679258
          3584 K	  68098867	  68465868
          3840 K	  72851637	  73244853
          4096 K	  77597294	  78016725
          4608 K	  87069012	  87540871
          5120 K	  96517023	  97041311
          5632 K	 105943724	 106520441
          6144 K	 115351074	 115980220
          6656 K	 124740700	 125422275
          7168 K	 134113980	 134847983	<*** F27 = 2^134217728 + 1

          7680 K	 143472090	 144258522
  8 M =   8192 K	 152816052	 153654913
  9 M =   9216 K	 171464992	 172408710
 10 M =  10240 K	 190066770	 191115346
 11 M =  11264 K	 208626152	 209779586
 12 M =  12288 K	 227147031	 228405322
 13 M =  13312 K	 245632644	 246995793
 14 M =  14336 K	 264085729	 265553736
 15 M =  15360 K	 282508628	 284081492	<*** F28 too large for 14M, needs >= 15M
 16 M =  16384 K	 300903371	 302581093
 18 M =  18432 K	 337615274	 339502711	<*** smallest 100-Mdigit moduli ***
 20 M =  20480 K	 374233313	 376330465
 22 M =  22528 K	 410766968	 413073835
 24 M =  24576 K	 447223981	 449740563
 26 M =  26624 K	 483610796	 486337093
 28 M =  28672 K	 519932856	 522868869
 30 M =  30720 K	 556194824	 559340552	<*** F29
 32 M =  32768 K	 592400738	 595756181	<*** Nov 2015: No ROE issues with a run of p = 595799947 [maxErr = 0.375], corr. to AsympConst ~= 0.4
 36 M =  36864 K	 664658102	 668432976
 40 M =  40960 K	 736728582	 740922886
 44 M =  45056 K	 808631042	 813244776
 48 M =  49152 K	 880380890	 885414055
 52 M =  53248 K	 951990950	 957443546
 56 M =  57344 K	1023472059	1029344085
 60 M =  61440 K	1094833496	1101124952	<*** F30
 64 M =  65536 K	1166083299	1172794185
 72 M =  73728 K	1308275271	1315825018
 80 M =  81920 K	1450095024	1458483632
 88 M =  90112 K	1591580114	1600807583
 96 M =  98304 K	1732761219	1742827549
104 M = 106496 K	1873663870	1884569060
112 M = 114688 K	2014309644	2026053695
120 M = 122880 K	2154717020	2167299932
128 M = 131072 K	2294902000	2308323773
144 M = 147456 K	2574659086	2589758580
160 M = 163840 K	2853674592	2870451808
176 M = 180224 K	3132023315	3150478252
192 M = 196608 K	3409766353	3429899013	<<<<*** Allows numbers slightly > 1Gdigit
208 M = 212992 K	3686954556	3708764937
224 M = 229376 K	3963630903	3987119006
240 M = 245760 K	4239832202	4264998026	<*** largest FFT length for 32-bit exponents
256 M = 262144 K	4515590327	4542433872
288 M = 294912 K	5065885246	5096084235
320 M = 327680 K	5614702299	5648256731
352 M = 360448 K	6162190494	6199100369
384 M = 393216 K	6708471554	6748736872
416 M = 425984 K	7253646785	7297267547
448 M = 458752 K	7797801823	7844778027
480 M = 491520 K	8341010002	8391341650
512 M = 524288 K	8883334834	8937021925	[0.60% larger]

For F33 have 2^33 = 8589934592, so could try a special FFT length based on leading radix = 31,
which should be just large enough to run without overly frequent fatal ROE > 0.4:

496 M = 507904 K	8612279197

But as radix-31 DFT is relatively inefficient in terms of flops/input, better to use the slightly
larger but much smoother odd radix 63, via 504 M = 516096 K . Here the next few Fermat numbers:
		AsympConst = 0.6:	AsympConst = 0.4:
  1 G		17472053085			17579427268		2^34 = 17179869184
  2 G		34356866077			34571614442		2^35 = 34359738368
  4 G		67542956696			67972453426		2^36 = 68719476736

So the transition to the power of 2 FFT length no longer being sufficient is very gradual, should be able to
handle F35 just fine, and even F36 should be OK, if we carefully maximize the ROE tolerance of our build,
say by improving the accuracy of the FFT-twiddles computation.
*/
