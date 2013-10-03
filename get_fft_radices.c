/*******************************************************************************
*                                                                              *
*   (C) 1997-2013 by Ernst W. Mayer.                                           *
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

#undef USE_ONLY_LARGE_LEAD_RADICES
#if defined(USE_SSE2) || defined(MULTITHREAD)
	#define USE_ONLY_LARGE_LEAD_RADICES
#endif

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

int	get_fft_radices(uint32 kblocks, int radix_set, int *nradices, int radix_vec[], int radix_vec_dim)
{
	int rvec[10];	/* Temporary storage for FFT radices */
	uint32 i, n,	numrad, rad_prod;
	// We arrange the allowed radix sets with the larger-leading-radices-first, as these typically (except for
	// very small FFT lengths) have an SSE2/AVX-specialized leading-radix DFT surrounding the carry step.
	// This mask permits only such radix sets in SSE2/AVX-enabled builds: In the index-out-of-range cases, we see
	//
	//      *nradices = A + (B & nrad_scalar_mask);
	//
	// Here, A is the number of SSE2-permitted radix sets (which must be listed first in the case set),
	// and B is the number of additional radix sets permitted for scalar (non-SSE2/AVX) builds.
	// The sum (A+B) must equal the total number of distinct radix-set cases listed for the given FFT length:
#ifdef USE_ONLY_LARGE_LEAD_RADICES
	int nrad_scalar_mask = 0;
#else
	int nrad_scalar_mask = -1;
#endif
	/*
	case selection based on real vector length - product of complex radices must equal complex vector length = N2.
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
	if((kblocks >> trailz32(kblocks)) > 15)
	{
		fprintf(stderr,"ERROR: specified FFT length %d K has non-power-of-2 component %d,\n", kblocks, (kblocks >> trailz32(kblocks)));
		fprintf(stderr,"       which exceeds the maximum allowed odd radix of 15.\n");
		return ERR_FFTLENGTH_ILLEGAL;
	}

	/* Real-array length, in CS notation (1K=1024 8-byte elements):	*/
	switch(kblocks)
	{
	case 2 :						/* 2K */
		switch(radix_set)
		{
		case 0 :
			numrad = 3; rvec[0] =  8; rvec[1] =  8; rvec[2] = 16; break;
		default :
			if(nradices){*nradices = 1 + (0 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 3 :						/* 3K */
		switch(radix_set)
		{
		case 0 :
			numrad = 3; rvec[0] = 12; rvec[1] =  8; rvec[2] = 16; break;	// Note: Radix-12 has no SSE2-optimized DFT, bu must support at least one radix set in vector-build mode
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 1 :
			numrad = 3; rvec[0] =  6; rvec[1] = 16; rvec[2] = 16; break;
	  #endif
		default :
			if(nradices){*nradices = 1 + (1 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 4 :						/* 4K */
		switch(radix_set)
		{
		case 0 :
			numrad = 3; rvec[0] =  8; rvec[1] = 16; rvec[2] = 16; break;
		case 1 :
			numrad = 3; rvec[0] = 16; rvec[1] =  8; rvec[2] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 2 :
			numrad = 3; rvec[0] =  8; rvec[1] =  8; rvec[2] = 32; break;
	  #endif
		default :
			if(nradices){*nradices = 2 + (1 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 5 :						/* 5K */
		switch(radix_set)
		{
		case 0 :
			numrad = 3; rvec[0] = 20; rvec[1] =  8; rvec[2] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 1 :
			numrad = 3; rvec[0] = 10; rvec[1] = 16; rvec[2] = 16; break;
		case 2 :
			numrad = 3; rvec[0] =  5; rvec[1] = 16; rvec[2] = 32; break;
	  #endif
		default :
			if(nradices){*nradices = 1 + (2 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 6 :						/* 6K */
		switch(radix_set)
		{
		case 0 :
			numrad = 3; rvec[0] = 24; rvec[1] =  8; rvec[2] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 1 :
			numrad = 3; rvec[0] = 12; rvec[1] = 16; rvec[2] = 16; break;
		case 2 :
			numrad = 3; rvec[0] =  6; rvec[1] = 16; rvec[2] = 32; break;
	  #endif
		default :
			if(nradices){*nradices = 1 + (2 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 7 :						/* 7K */
		switch(radix_set)
		{
		case 0 :
			numrad = 3; rvec[0] = 28; rvec[1] =  8; rvec[2] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 1 :
			numrad = 3; rvec[0] = 14; rvec[1] = 16; rvec[2] = 16; break;
		case 2 :
			numrad = 3; rvec[0] =  7; rvec[1] = 16; rvec[2] = 32; break;
	  #endif
		default :
			if(nradices){*nradices = 1 + (2 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 8 :						/* 8K */
		switch(radix_set)
		{
		case 0 :
			numrad = 3; rvec[0] = 32; rvec[1] =  8; rvec[2] = 16; break;
		case 1 :
			numrad = 3; rvec[0] = 16; rvec[1] = 16; rvec[2] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 2 :
			numrad = 3; rvec[0] =  8; rvec[1] = 16; rvec[2] = 32; break;
	  #endif
		default :
			if(nradices){*nradices = 2 + (1 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 9 :						/* 9K */
		switch(radix_set)
		{
		case 0 :
			numrad = 3; rvec[0] = 36; rvec[1] =  8; rvec[2] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 1 :
			numrad = 3; rvec[0] = 18; rvec[1] = 16; rvec[2] = 16; break;
		case 2 :
			numrad = 3; rvec[0] =  9; rvec[1] = 32; rvec[2] = 16; break;
		case 3 :
			numrad = 3; rvec[0] =  9; rvec[1] = 16; rvec[2] = 32; break;
	  #endif
		default :
			if(nradices){*nradices = 1 + (3 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 10 :						/* 10K */
		switch(radix_set)
		{
		case 0 :
			numrad = 3; rvec[0] = 40; rvec[1] =  8; rvec[2] = 16; break;
		case 1 :
			numrad = 3; rvec[0] = 20; rvec[1] = 16; rvec[2] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 2 :
			numrad = 3; rvec[0] = 10; rvec[1] = 32; rvec[2] = 16; break;
		case 3 :
			numrad = 3; rvec[0] = 10; rvec[1] = 16; rvec[2] = 32; break;
		case 4 :
			numrad = 3; rvec[0] =  5; rvec[1] = 32; rvec[2] = 32; break;
		case 5 :
			numrad = 4; rvec[0] =  5; rvec[1] =  8; rvec[2] =  8; rvec[3] = 16; break;
	  #endif
		default :
			if(nradices){*nradices = 2 + (4 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 11 :						/* 11K */
		switch(radix_set)
		{
		case 0 :
			numrad = 3; rvec[0] = 44; rvec[1] =  8; rvec[2] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 1 :
			numrad = 3; rvec[0] = 22; rvec[1] = 16; rvec[2] = 16; break;
		case 2 :
			numrad = 3; rvec[0] = 11; rvec[1] = 32; rvec[2] = 16; break;
		case 3 :
			numrad = 3; rvec[0] = 11; rvec[1] = 16; rvec[2] = 32; break;
	  #endif
		default :
			if(nradices){*nradices = 1 + (3 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 12 :						/* 12K */
		switch(radix_set)
		{
		case 0 :
			numrad = 3; rvec[0] = 24; rvec[1] = 16; rvec[2] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 1 :
			numrad = 3; rvec[0] = 12; rvec[1] = 32; rvec[2] = 16; break;
		case 2 :
			numrad = 3; rvec[0] = 12; rvec[1] = 16; rvec[2] = 32; break;
		case 3 :
			numrad = 3; rvec[0] =  6; rvec[1] = 32; rvec[2] = 32; break;
		case 4 :
			numrad = 4; rvec[0] =  6; rvec[1] =  8; rvec[2] =  8; rvec[3] = 16; break;
	  #endif
		default :
			if(nradices){*nradices = 1 + (4 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 13 :						/* 13K */
		switch(radix_set)
		{
		case 0 :
			numrad = 3; rvec[0] = 52; rvec[1] =  8; rvec[2] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 1 :
			numrad = 3; rvec[0] = 26; rvec[1] = 16; rvec[2] = 16; break;
		case 2 :
			numrad = 3; rvec[0] = 13; rvec[1] = 32; rvec[2] = 16; break;
		case 3 :
			numrad = 3; rvec[0] = 13; rvec[1] = 16; rvec[2] = 32; break;
	  #endif
		default :
			if(nradices){*nradices = 1 + (3 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 14 :						/* 14K */
		switch(radix_set)
		{
		case 0 :
			numrad = 3; rvec[0] = 28; rvec[1] = 16; rvec[2] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 1 :
			numrad = 3; rvec[0] = 14; rvec[1] = 32; rvec[2] = 16; break;
		case 2 :
			numrad = 3; rvec[0] = 14; rvec[1] = 16; rvec[2] = 32; break;
		case 3 :
			numrad = 3; rvec[0] =  7; rvec[1] = 32; rvec[2] = 32; break;
		case 4 :
			numrad = 4; rvec[0] =  7; rvec[1] =  8; rvec[2] =  8; rvec[3] = 16; break;
	  #endif
		default :
			if(nradices){*nradices = 1 + (4 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 15 :						/* 15K */
		switch(radix_set)
		{
		case 0 :
			numrad = 3; rvec[0] = 60; rvec[1] =  8; rvec[2] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 1 :
			numrad = 3; rvec[0] = 30; rvec[1] = 16; rvec[2] = 16; break;
		case 2 :
			numrad = 3; rvec[0] = 15; rvec[1] = 32; rvec[2] = 16; break;
		case 3 :
			numrad = 3; rvec[0] = 15; rvec[1] = 16; rvec[2] = 32; break;
	  #endif
		default :
			if(nradices){*nradices = 1 + (3 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 16 :						/* 16K */
		switch(radix_set)
		{
		case 0 :
			numrad = 3; rvec[0] = 32; rvec[1] = 16; rvec[2] = 16; break;
		case 1 :
			numrad = 3; rvec[0] = 16; rvec[1] = 32; rvec[2] = 16; break;
		case 2 :
			numrad = 3; rvec[0] = 16; rvec[1] = 16; rvec[2] = 32; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 3 :
			numrad = 3; rvec[0] =  8; rvec[1] = 32; rvec[2] = 32; break;
		case 4 :
			numrad = 4; rvec[0] =  8; rvec[1] =  8; rvec[2] =  8; rvec[3] = 16; break;
	  #endif
		default :
			if(nradices){*nradices = 3 + (2 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 18 :						/* 18K */
		switch(radix_set)
		{
		case 0 :
			numrad = 3; rvec[0] = 36; rvec[1] = 16; rvec[2] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 1 :
			numrad = 3; rvec[0] = 18; rvec[1] = 32; rvec[2] = 16; break;
		case 2 :
			numrad = 3; rvec[0] = 18; rvec[1] = 16; rvec[2] = 32; break;
		case 3 :
			numrad = 3; rvec[0] =  9; rvec[1] = 32; rvec[2] = 32; break;
		case 4 :
			numrad = 4; rvec[0] =  9; rvec[1] =  8; rvec[2] =  8; rvec[3] = 16; break;
	  #endif
		default :
			if(nradices){*nradices = 1 + (4 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 20 :						/* 20K */
		switch(radix_set)
		{
		case 0 :
			numrad = 3; rvec[0] = 40; rvec[1] = 16; rvec[2] = 16; break;
		case 1 :
			numrad = 3; rvec[0] = 20; rvec[1] = 32; rvec[2] = 16; break;
		case 2 :
			numrad = 3; rvec[0] = 20; rvec[1] = 16; rvec[2] = 32; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 3 :
			numrad = 3; rvec[0] = 10; rvec[1] = 32; rvec[2] = 32; break;
		case 4 :
			numrad = 4; rvec[0] =  5; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; break;
	  #endif
		default :
			if(nradices){*nradices = 3 + (2 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 22 :						/* 22K */
		switch(radix_set)
		{
		case 0 :
			numrad = 3; rvec[0] = 44; rvec[1] = 16; rvec[2] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 1 :
			numrad = 3; rvec[0] = 22; rvec[1] = 32; rvec[2] = 16; break;
		case 2 :
			numrad = 3; rvec[0] = 22; rvec[1] = 16; rvec[2] = 32; break;
		case 3 :
			numrad = 3; rvec[0] = 11; rvec[1] = 32; rvec[2] = 32; break;
		case 4 :
			numrad = 4; rvec[0] = 11; rvec[1] =  8; rvec[2] =  8; rvec[3] = 16; break;
	  #endif
		default :
			if(nradices){*nradices = 1 + (4 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 24 :						/* 24K */
		switch(radix_set)
		{
		case 0 :
			numrad = 3; rvec[0] = 24; rvec[1] = 32; rvec[2] = 16; break;
		case 1 :
			numrad = 3; rvec[0] = 24; rvec[1] = 16; rvec[2] = 32; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 2 :
			numrad = 3; rvec[0] = 12; rvec[1] = 32; rvec[2] = 32; break;
		case 3 :
			numrad = 4; rvec[0] =  6; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; break;
	  #endif
		default :
			if(nradices){*nradices = 2 + (2 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 26 :						/* 26K */
		switch(radix_set)
		{
		case 0 :
			numrad = 3; rvec[0] = 52; rvec[1] = 16; rvec[2] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 1 :
			numrad = 3; rvec[0] = 26; rvec[1] = 32; rvec[2] = 16; break;
		case 2 :
			numrad = 3; rvec[0] = 26; rvec[1] = 16; rvec[2] = 32; break;
		case 3 :
			numrad = 3; rvec[0] = 13; rvec[1] = 32; rvec[2] = 32; break;
		case 4 :
			numrad = 4; rvec[0] = 13; rvec[1] =  8; rvec[2] =  8; rvec[3] = 16; break;
	  #endif
		default :
			if(nradices){*nradices = 1 + (4 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 28 :						/* 28K */
		switch(radix_set)
		{
		case 0 :
			numrad = 3; rvec[0] = 28; rvec[1] = 32; rvec[2] = 16; break;
		case 1 :
			numrad = 3; rvec[0] = 28; rvec[1] = 16; rvec[2] = 32; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 2 :
			numrad = 3; rvec[0] = 14; rvec[1] = 32; rvec[2] = 32; break;
		case 3 :
			numrad = 4; rvec[0] =  7; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; break;
	  #endif
		default :
			if(nradices){*nradices = 2 + (2 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 30 :						/* 30K */
		switch(radix_set)
		{
		case 0 :
			numrad = 3; rvec[0] = 60; rvec[1] = 16; rvec[2] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 1 :
			numrad = 3; rvec[0] = 30; rvec[1] = 32; rvec[2] = 16; break;
		case 2 :
			numrad = 3; rvec[0] = 30; rvec[1] = 16; rvec[2] = 32; break;
		case 3 :
			numrad = 3; rvec[0] = 15; rvec[1] = 32; rvec[2] = 32; break;
		case 4 :
			numrad = 4; rvec[0] = 15; rvec[1] =  8; rvec[2] =  8; rvec[3] = 16; break;
	  #endif
		default :
			if(nradices){*nradices = 1 + (4 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 32 :						/* 32K */
		switch(radix_set)
		{
		case 0 :
			numrad = 3; rvec[0] = 32; rvec[1] = 32; rvec[2] = 16; break;
		case 1 :
			numrad = 3; rvec[0] = 32; rvec[1] = 16; rvec[2] = 32; break;
		case 2 :
			numrad = 3; rvec[0] = 16; rvec[1] = 32; rvec[2] = 32; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 3 :
			numrad = 4; rvec[0] =  8; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; break;
	  #endif
		default :
			if(nradices){*nradices = 3 + (1 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 36 :						/* 36K */
		switch(radix_set)
		{
		case 0 :
			numrad = 3; rvec[0] = 36; rvec[1] = 32; rvec[2] = 16; break;
		case 1 :
			numrad = 3; rvec[0] = 36; rvec[1] = 16; rvec[2] = 32; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 2 :
			numrad = 3; rvec[0] = 18; rvec[1] = 32; rvec[2] = 32; break;
		case 3 :
			numrad = 4; rvec[0] =  9; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; break;
	  #endif
		default :
			if(nradices){*nradices = 2 + (2 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 40 :						/* 40K */
		switch(radix_set)
		{
		case 0 :
			numrad = 3; rvec[0] = 40; rvec[1] = 32; rvec[2] = 16; break;
		case 1 :
			numrad = 3; rvec[0] = 40; rvec[1] = 16; rvec[2] = 32; break;
		case 2 :
			numrad = 3; rvec[0] = 20; rvec[1] = 32; rvec[2] = 32; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 3 :
			numrad = 4; rvec[0] = 10; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; break;
	  #endif
		default :
			if(nradices){*nradices = 3 + (1 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 44 :						/* 44K */
		switch(radix_set)
		{
		case 0 :
			numrad = 3; rvec[0] = 44; rvec[1] = 32; rvec[2] = 16; break;
		case 1 :
			numrad = 3; rvec[0] = 44; rvec[1] = 16; rvec[2] = 32; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 2 :
			numrad = 3; rvec[0] = 22; rvec[1] = 32; rvec[2] = 32; break;
		case 3 :
			numrad = 4; rvec[0] = 11; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; break;
	  #endif
		default :
			if(nradices){*nradices = 2 + (2 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 48 :						/* 48K */
		switch(radix_set)
		{
		case 0 :
			numrad = 3; rvec[0] = 24; rvec[1] = 32; rvec[2] = 32; break;
		case 1 :
			numrad = 4; rvec[0] = 24; rvec[1] =  8; rvec[2] =  8; rvec[3] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 2 :
			numrad = 4; rvec[0] = 12; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; break;
		case 3 :
			numrad = 4; rvec[0] =  6; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; break;
	  #endif
		default :
			if(nradices){*nradices = 2 + (2 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 52 :						/* 52K */
		switch(radix_set)
		{
		case 0 :
			numrad = 3; rvec[0] = 52; rvec[1] = 32; rvec[2] = 16; break;
		case 1 :
			numrad = 3; rvec[0] = 52; rvec[1] = 16; rvec[2] = 32; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 2 :
			numrad = 3; rvec[0] = 26; rvec[1] = 32; rvec[2] = 32; break;
		case 3 :
			numrad = 4; rvec[0] = 13; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; break;
	  #endif
		default :
			if(nradices){*nradices = 2 + (2 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 56 :						/* 56K */
		switch(radix_set)
		{
		case 0 :
			numrad = 3; rvec[0] = 28; rvec[1] = 32; rvec[2] = 32; break;
		case 1 :
			numrad = 4; rvec[0] = 28; rvec[1] =  8; rvec[2] =  8; rvec[3] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 2 :
			numrad = 4; rvec[0] = 14; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; break;
		case 3 :
			numrad = 4; rvec[0] =  7; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; break;
	  #endif
		default :
			if(nradices){*nradices = 2 + (2 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 60 :						/* 60K */
		switch(radix_set)
		{
		case 0 :
			numrad = 3; rvec[0] = 60; rvec[1] = 32; rvec[2] = 16; break;
		case 1 :
			numrad = 3; rvec[0] = 60; rvec[1] = 16; rvec[2] = 32; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 2 :
			numrad = 3; rvec[0] = 30; rvec[1] = 32; rvec[2] = 32; break;
		case 3 :
			numrad = 4; rvec[0] = 15; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; break;
	  #endif
		default :
			if(nradices){*nradices = 2 + (2 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 64 :						/* 64K */
		switch(radix_set)
		{
		case 0 :
			numrad = 3; rvec[0] = 32; rvec[1] = 32; rvec[2] = 32; break;
		case 1 :
			numrad = 4; rvec[0] = 32; rvec[1] =  8; rvec[2] =  8; rvec[3] = 16; break;
		case 2 :
			numrad = 4; rvec[0] = 16; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 3 :
			numrad = 4; rvec[0] =  8; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; break;
	  #endif
		default :
			if(nradices){*nradices = 3 + (1 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 72 :						/* 72K */
		switch(radix_set)
		{
		case 0 :
			numrad = 3; rvec[0] = 36; rvec[1] = 32; rvec[2] = 32; break;
		case 1 :
			numrad = 4; rvec[0] = 36; rvec[1] =  8; rvec[2] =  8; rvec[3] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 2 :
			numrad = 4; rvec[0] = 18; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; break;
		case 3 :
			numrad = 4; rvec[0] =  9; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; break;
	  #endif
		default :
			if(nradices){*nradices = 2 + (2 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 80 :						/* 80K */
		switch(radix_set)
		{
		case 0 :
			numrad = 3; rvec[0] = 40; rvec[1] = 32; rvec[2] = 32; break;
		case 1 :
			numrad = 4; rvec[0] = 40; rvec[1] =  8; rvec[2] =  8; rvec[3] = 16; break;
		case 2 :
			numrad = 4; rvec[0] = 20; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 3 :
			numrad = 4; rvec[0] = 10; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; break;
	  #endif
		default :
			if(nradices){*nradices = 3 + (1 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 88 :						/* 88K */
		switch(radix_set)
		{
		case 0 :
			numrad = 3; rvec[0] = 44; rvec[1] = 32; rvec[2] = 32; break;
		case 1 :
			numrad = 4; rvec[0] = 44; rvec[1] =  8; rvec[2] =  8; rvec[3] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 2 :
			numrad = 4; rvec[0] = 22; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; break;
		case 3 :
			numrad = 4; rvec[0] = 11; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; break;
	  #endif
		default :
			if(nradices){*nradices = 2 + (2 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 96 :						/* 96K */
		switch(radix_set)
		{
		case 0 :
			numrad = 3; rvec[0] = 48; rvec[1] = 32; rvec[2] = 32; break;
		case 1 :
			numrad = 4; rvec[0] = 24; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 2 :
			numrad = 4; rvec[0] = 12; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; break;
		case 3 :
			numrad = 4; rvec[0] =  6; rvec[1] = 16; rvec[2] = 32; rvec[3] = 16; break;
		case 4 :
			numrad = 4; rvec[0] =  6; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; break;
	  #endif
		default :
			if(nradices){*nradices = 2 + (3 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 104 :						/* 104K */
		switch(radix_set)
		{
		case 0 :
			numrad = 3; rvec[0] = 52; rvec[1] = 32; rvec[2] = 32; break;
		case 1 :
			numrad = 4; rvec[0] = 52; rvec[1] =  8; rvec[2] =  8; rvec[3] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 2 :
			numrad = 4; rvec[0] = 26; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; break;
		case 3 :
			numrad = 4; rvec[0] = 13; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; break;
	  #endif
		default :
			if(nradices){*nradices = 2 + (2 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 112 :						/* 112K */
		switch(radix_set)
		{
		case 0 :
			numrad = 3; rvec[0] = 56; rvec[1] = 32; rvec[2] = 32; break;
		case 1 :
			numrad = 4; rvec[0] = 28; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 2 :
			numrad = 4; rvec[0] = 14; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; break;
		case 3 :
			numrad = 4; rvec[0] =  7; rvec[1] = 16; rvec[2] = 32; rvec[3] = 16; break;
		case 4 :
			numrad = 4; rvec[0] =  7; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; break;
	  #endif
		default :
			if(nradices){*nradices = 2 + (3 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 120 :						/* 120K */
		switch(radix_set)
		{
		case 0 :
			numrad = 3; rvec[0] = 60; rvec[1] = 32; rvec[2] = 32; break;
		case 1 :
			numrad = 4; rvec[0] = 60; rvec[1] =  8; rvec[2] =  8; rvec[3] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 2 :
			numrad = 4; rvec[0] = 30; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; break;
		case 3 :
			numrad = 4; rvec[0] = 15; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; break;
	  #endif
		default :
			if(nradices){*nradices = 2 + (2 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 128 :						/* 128K */
		switch(radix_set)
		{
		case 0 :
			numrad = 3; rvec[0] = 64; rvec[1] = 32; rvec[2] = 32; break;
		case 1 :
			numrad = 4; rvec[0] = 32; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; break;
		case 2 :
			numrad = 4; rvec[0] = 16; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 3 :
			numrad = 4; rvec[0] =  8; rvec[1] = 16; rvec[2] = 32; rvec[3] = 16; break;
		case 4 :
			numrad = 4; rvec[0] =  8; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; break;
	  #endif
		default :
			if(nradices){*nradices = 3 + (2 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 144 :						/* 144K */
		switch(radix_set)
		{
		case 0 :
			numrad = 4; rvec[0] = 36; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 1 :
			numrad = 4; rvec[0] = 18; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; break;
		case 2 :
			numrad = 4; rvec[0] =  9; rvec[1] = 16; rvec[2] = 32; rvec[3] = 16; break;
		case 3 :
			numrad = 4; rvec[0] =  9; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; break;
	  #endif
		default :
			if(nradices){*nradices = 1 + (3 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 160 :						/* 160K */
		switch(radix_set)
		{
		case 0 :
			numrad = 4; rvec[0] = 40; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; break;
		case 1 :
			numrad = 4; rvec[0] = 20; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 2 :
			numrad = 4; rvec[0] = 10; rvec[1] = 16; rvec[2] = 32; rvec[3] = 16; break;
		case 3 :
			numrad = 4; rvec[0] = 10; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; break;
		case 4 :
			numrad = 4; rvec[0] =  5; rvec[1] = 32; rvec[2] = 32; rvec[3] = 16; break;
	  #endif
		default :
			if(nradices){*nradices = 2 + (3 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 176 :						/* 176K */
		switch(radix_set)
		{
		case 0 :
			numrad = 4; rvec[0] = 44; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 1 :
			numrad = 4; rvec[0] = 22; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; break;
		case 2 :
			numrad = 4; rvec[0] = 11; rvec[1] = 16; rvec[2] = 32; rvec[3] = 16; break;
		case 3 :
			numrad = 4; rvec[0] = 11; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; break;
	  #endif
		default :
			if(nradices){*nradices = 1 + (3 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 192 :						/* 192K */
		switch(radix_set)
		{
		case 0 :
			numrad = 4; rvec[0] = 48; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; break;
		case 1 :
			numrad = 4; rvec[0] = 24; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 2 :
			numrad = 4; rvec[0] = 12; rvec[1] = 16; rvec[2] = 32; rvec[3] = 16; break;
		case 3 :
			numrad = 4; rvec[0] = 12; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; break;
		case 4 :
			numrad = 4; rvec[0] =  6; rvec[1] = 32; rvec[2] = 32; rvec[3] = 16; break;
		case 5 :
			numrad = 4; rvec[0] =  6; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; break;
		case 6 :
			numrad = 5; rvec[0] =  6; rvec[1] =  8; rvec[2] =  8; rvec[3] = 16; rvec[4] = 16; break;
	  #endif
		default :
			if(nradices){*nradices = 2 + (5 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 208 :						/* 208K */
		switch(radix_set)
		{
		case 0 :
			numrad = 4; rvec[0] = 52; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 1 :
			numrad = 4; rvec[0] = 26; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; break;
		case 2 :
			numrad = 4; rvec[0] = 13; rvec[1] = 16; rvec[2] = 32; rvec[3] = 16; break;
		case 3 :
			numrad = 4; rvec[0] = 13; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; break;
	  #endif
		default :
			if(nradices){*nradices = 1 + (3 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 224 :						/* 224K */
		switch(radix_set)
		{
		case 0 :
			numrad = 4; rvec[0] = 56; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; break;
		case 1 :
			numrad = 4; rvec[0] = 28; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 2 :
			numrad = 4; rvec[0] = 14; rvec[1] = 16; rvec[2] = 32; rvec[3] = 16; break;
		case 3 :
			numrad = 4; rvec[0] = 14; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; break;
		case 4 :
			numrad = 4; rvec[0] =  7; rvec[1] = 32; rvec[2] = 32; rvec[3] = 16; break;
		case 5 :
			numrad = 4; rvec[0] =  7; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; break;
		case 6 :
			numrad = 5; rvec[0] =  7; rvec[1] =  8; rvec[2] =  8; rvec[3] = 16; rvec[4] = 16; break;
	  #endif
		default :
			if(nradices){*nradices = 2 + (5 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 240 :						/* 240K */
		switch(radix_set)
		{
		case 0 :
			numrad = 4; rvec[0] = 60; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 1 :
			numrad = 4; rvec[0] = 30; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; break;
		case 2 :
			numrad = 4; rvec[0] = 15; rvec[1] = 16; rvec[2] = 32; rvec[3] = 16; break;
		case 3 :
			numrad = 4; rvec[0] = 15; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; break;
	  #endif
		default :
			if(nradices){*nradices = 1 + (3 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 256 :						/* 256K */
		switch(radix_set)
		{
		case 0 :
			numrad = 4; rvec[0] = 64; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; break;
		case 1 :
			numrad = 4; rvec[0] = 32; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; break;
		case 2 :
			numrad = 4; rvec[0] = 16; rvec[1] = 16; rvec[2] = 32; rvec[3] = 16; break;
		case 3 :
			numrad = 4; rvec[0] = 16; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 4 :
			numrad = 4; rvec[0] =  8; rvec[1] = 32; rvec[2] = 32; rvec[3] = 16; break;
		case 5 :
			numrad = 4; rvec[0] =  8; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; break;
		case 6 :
			numrad = 5; rvec[0] =  8; rvec[1] =  8; rvec[2] =  8; rvec[3] = 16; rvec[4] = 16; break;
	  #endif
		default :
			if(nradices){*nradices = 4 + (3 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 288 :						/* 288K */
		switch(radix_set)
		{
		case 0 :
			numrad = 4; rvec[0] = 36; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 1 :
			numrad = 4; rvec[0] = 18; rvec[1] = 16; rvec[2] = 32; rvec[3] = 16; break;
		case 2 :
			numrad = 4; rvec[0] = 18; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; break;
		case 3 :
			numrad = 4; rvec[0] =  9; rvec[1] = 32; rvec[2] = 32; rvec[3] = 16; break;
		case 4 :
			numrad = 4; rvec[0] =  9; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; break;
		case 5 :
			numrad = 5; rvec[0] =  9; rvec[1] =  8; rvec[2] =  8; rvec[3] = 16; rvec[4] = 16; break;
	  #endif
		default :
			if(nradices){*nradices = 1 + (5 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 320 :						/* 320K */
		switch(radix_set)
		{
		case 0 :
			numrad = 4; rvec[0] = 40; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; break;
		case 1 :
			numrad = 4; rvec[0] = 20; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; break;
		case 2 :
			numrad = 4; rvec[0] = 20; rvec[1] = 16; rvec[2] = 32; rvec[3] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 3 :
			numrad = 4; rvec[0] = 10; rvec[1] = 32; rvec[2] = 32; rvec[3] = 16; break;
		case 4 :
			numrad = 4; rvec[0] = 10; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; break;
		case 5 :
			numrad = 5; rvec[0] =  5; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
	  #endif
		default :
			if(nradices){*nradices = 3 + (3 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 352 :						/* 352K */
		switch(radix_set)
		{
		case 0 :
			numrad = 4; rvec[0] = 44; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 1 :
			numrad = 4; rvec[0] = 22; rvec[1] = 16; rvec[2] = 32; rvec[3] = 16; break;
		case 2 :
			numrad = 4; rvec[0] = 22; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; break;
		case 3 :
			numrad = 4; rvec[0] = 11; rvec[1] = 32; rvec[2] = 32; rvec[3] = 16; break;
		case 4 :
			numrad = 4; rvec[0] = 11; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; break;
	  #endif
		default :
			if(nradices){*nradices = 1 + (4 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 384 :						/* 384K */
		switch(radix_set)
		{
		case 0 :
			numrad = 4; rvec[0] = 48; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; break;
		case 1 :
			numrad = 4; rvec[0] = 24; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; break;
		case 2 :
			numrad = 4; rvec[0] = 24; rvec[1] = 16; rvec[2] = 32; rvec[3] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 3 :
			numrad = 4; rvec[0] = 12; rvec[1] = 32; rvec[2] = 32; rvec[3] = 16; break;
		case 4 :
			numrad = 4; rvec[0] = 12; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; break;
		case 5 :
			numrad = 5; rvec[0] =  6; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
	  #endif
		default :
			if(nradices){*nradices = 3 + (3 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
  	case 416 :						/* 416K */
		switch(radix_set)
		{
		case 0 :
			numrad = 4; rvec[0] = 52; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 1 :
			numrad = 4; rvec[0] = 26; rvec[1] = 16; rvec[2] = 32; rvec[3] = 16; break;
		case 2 :
			numrad = 4; rvec[0] = 26; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; break;
		case 3 :
			numrad = 4; rvec[0] = 13; rvec[1] = 32; rvec[2] = 32; rvec[3] = 16; break;
		case 4 :
			numrad = 4; rvec[0] = 13; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; break;
	  #endif
		default :
			if(nradices){*nradices = 1 + (4 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 448 :						/* 448K */
		switch(radix_set)
		{
		case 0 :
			numrad = 4; rvec[0] = 56; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; break;
		case 1 :
			numrad = 4; rvec[0] = 28; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; break;
		case 2 :
			numrad = 4; rvec[0] = 28; rvec[1] = 16; rvec[2] = 32; rvec[3] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 3 :
			numrad = 4; rvec[0] = 14; rvec[1] = 32; rvec[2] = 32; rvec[3] = 16; break;
		case 4 :
			numrad = 4; rvec[0] = 14; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; break;
		case 5 :
			numrad = 5; rvec[0] =  7; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
	  #endif
		default :
			if(nradices){*nradices = 3 + (3 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 480 :						/* 480K */
		switch(radix_set)
		{
		case 0 :
			numrad = 4; rvec[0] = 60; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 1 :
			numrad = 4; rvec[0] = 30; rvec[1] = 16; rvec[2] = 32; rvec[3] = 16; break;
		case 2 :
			numrad = 4; rvec[0] = 30; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; break;
		case 3 :
			numrad = 4; rvec[0] = 15; rvec[1] = 32; rvec[2] = 32; rvec[3] = 16; break;
		case 4 :
			numrad = 4; rvec[0] = 15; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; break;
	  #endif
		default :
			if(nradices){*nradices = 1 + (4 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 512 :						/* 512K */
		switch(radix_set)
		{
		case 0 :
			numrad = 4; rvec[0] = 64; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; break;
		case 1 :
			numrad = 4; rvec[0] = 32; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; break;
		case 2 :
			numrad = 4; rvec[0] = 32; rvec[1] = 16; rvec[2] = 32; rvec[3] = 16; break;
		case 3 :
			numrad = 4; rvec[0] = 16; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; break;
		case 4 :
			numrad = 4; rvec[0] = 16; rvec[1] = 32; rvec[2] = 32; rvec[3] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 5 :
			numrad = 5; rvec[0] =  8; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
	  #endif
		default :
			if(nradices){*nradices = 5 + (1 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 576 :						/* 576K */
		switch(radix_set)
		{
		case 0 :
			numrad = 4; rvec[0] = 36; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; break;
		case 1 :
			numrad = 4; rvec[0] = 36; rvec[1] = 16; rvec[2] = 32; rvec[3] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 2 :
			numrad = 4; rvec[0] = 18; rvec[1] = 32; rvec[2] = 32; rvec[3] = 16; break;
		case 3 :
			numrad = 4; rvec[0] = 18; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; break;
		case 4 :
			numrad = 5; rvec[0] =  9; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
	  #endif
		default :
			if(nradices){*nradices = 2 + (3 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 640 :						/* 640K */
		switch(radix_set)
		{
		case 0 :
			numrad = 4; rvec[0] = 40; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; break;
		case 1 :
			numrad = 4; rvec[0] = 40; rvec[1] = 16; rvec[2] = 32; rvec[3] = 16; break;
		case 2 :
			numrad = 4; rvec[0] = 20; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; break;
		case 3 :
			numrad = 4; rvec[0] = 20; rvec[1] = 32; rvec[2] = 32; rvec[3] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 4 :
			numrad = 5; rvec[0] = 10; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
	  #endif
		default :
			if(nradices){*nradices = 4 + (1 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 704 :						/* 704K */
		switch(radix_set)
		{
		case 0 :
			numrad = 4; rvec[0] = 44; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; break;
		case 1 :
			numrad = 4; rvec[0] = 44; rvec[1] = 16; rvec[2] = 32; rvec[3] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 2 :
			numrad = 4; rvec[0] = 22; rvec[1] = 32; rvec[2] = 32; rvec[3] = 16; break;
		case 3 :
			numrad = 4; rvec[0] = 22; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; break;
		case 4 :
			numrad = 5; rvec[0] = 11; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
	  #endif
		default :
			if(nradices){*nradices = 2 + (3 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 768 :						/* 768K */
		switch(radix_set)
		{
		case 0 :
			numrad = 4; rvec[0] = 48; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; break;
		case 1 :
			numrad = 4; rvec[0] = 48; rvec[1] = 16; rvec[2] = 32; rvec[3] = 16; break;
		case 2 :
			numrad = 4; rvec[0] = 24; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; break;
		case 3 :
			numrad = 4; rvec[0] = 24; rvec[1] = 32; rvec[2] = 32; rvec[3] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 4 :
			numrad = 4; rvec[0] = 12; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; break;
		case 5 :
			numrad = 5; rvec[0] = 12; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
		case 6 :
			numrad = 5; rvec[0] =  6; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
	  #endif
		default :
			if(nradices){*nradices = 4 + (3 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 832 :						/* 832K */
		switch(radix_set)
		{
		case 0 :
			numrad = 4; rvec[0] = 52; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; break;
		case 1 :
			numrad = 4; rvec[0] = 52; rvec[1] = 16; rvec[2] = 32; rvec[3] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 2 :
			numrad = 4; rvec[0] = 26; rvec[1] = 32; rvec[2] = 32; rvec[3] = 16; break;
		case 3 :
			numrad = 4; rvec[0] = 26; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; break;
		case 4 :
			numrad = 5; rvec[0] = 13; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
	  #endif
		default :
			if(nradices){*nradices = 2 + (3 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 896 :						/* 896K */
		switch(radix_set)
		{
		case 0 :
			numrad = 4; rvec[0] = 56; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; break;
		case 1 :
			numrad = 4; rvec[0] = 56; rvec[1] = 16; rvec[2] = 32; rvec[3] = 16; break;
		case 2 :
			numrad = 4; rvec[0] = 28; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; break;
		case 3 :
			numrad = 4; rvec[0] = 28; rvec[1] = 32; rvec[2] = 32; rvec[3] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 4 :
			numrad = 4; rvec[0] = 14; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; break;
		case 5 :
			numrad = 5; rvec[0] = 14; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
		case 6 :
			numrad = 5; rvec[0] =  7; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
	  #endif
		default :
			if(nradices){*nradices = 4 + (3 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 960 :						/* 960K */
		switch(radix_set)
		{
		case 0 :
			numrad = 4; rvec[0] = 60; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; break;
		case 1 :
			numrad = 4; rvec[0] = 60; rvec[1] = 16; rvec[2] = 32; rvec[3] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 2 :
			numrad = 4; rvec[0] = 30; rvec[1] = 32; rvec[2] = 32; rvec[3] = 16; break;
		case 3 :
			numrad = 4; rvec[0] = 30; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; break;
		case 4 :
			numrad = 5; rvec[0] = 15; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
	  #endif
		default :
			if(nradices){*nradices = 2 + (3 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 1024 :					/* 1M = 1024K */
		switch(radix_set)
		{
		case 0 :
			numrad = 4; rvec[0] = 64; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; break;
		case 1 :
			numrad = 4; rvec[0] = 32; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; break;
		case 2 :
			numrad = 4; rvec[0] = 32; rvec[1] = 32; rvec[2] = 32; rvec[3] = 16; break;
		case 3 :
			numrad = 4; rvec[0] = 16; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; break;
		case 4 :
			numrad = 5; rvec[0] = 16; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 5 :
			numrad = 5; rvec[0] =  8; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
	  #endif
		default :
			if(nradices){*nradices = 5 + (1 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 1152 :					/* 1.125M = 1152K */
		switch(radix_set)
		{
		case 0 :
			numrad = 4; rvec[0] = 36; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; break;
		case 1 :
			numrad = 4; rvec[0] = 36; rvec[1] = 32; rvec[2] = 32; rvec[3] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 2 :
			numrad = 4; rvec[0] = 18; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; break;
		case 3 :
			numrad = 5; rvec[0] = 18; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
		case 4 :
			numrad = 5; rvec[0] =  9; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
	  #endif
		default :
			if(nradices){*nradices = 2 + (3 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 1280 :					/* 1.25M = 1280K */
		switch(radix_set)
		{
		case 0 :
			numrad = 4; rvec[0] = 40; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; break;
		case 1 :
			numrad = 4; rvec[0] = 40; rvec[1] = 32; rvec[2] = 32; rvec[3] = 16; break;
		case 2 :
			numrad = 4; rvec[0] = 20; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; break;
		case 3 :
			numrad = 5; rvec[0] = 20; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 4 :
			numrad = 5; rvec[0] = 10; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
	  #endif
		default :
			if(nradices){*nradices = 4 + (1 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 1408 :					/* 1.375M = 1408K */
		switch(radix_set)
		{
		case 0 :
			numrad = 4; rvec[0] = 44; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; break;
		case 1 :
			numrad = 4; rvec[0] = 44; rvec[1] = 32; rvec[2] = 32; rvec[3] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 2 :
			numrad = 4; rvec[0] = 22; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; break;
		case 3 :
			numrad = 5; rvec[0] = 22; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
		case 4 :
			numrad = 5; rvec[0] = 11; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
	  #endif
		default :
			if(nradices){*nradices = 2 + (3 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 1536 :					/* 1.5M = 1536K */
		switch(radix_set)
		{
		case 0 :
			numrad = 4; rvec[0] = 48; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; break;
		case 1 :
			numrad = 4; rvec[0] = 48; rvec[1] = 32; rvec[2] = 32; rvec[3] = 16; break;
		case 2 :
			numrad = 4; rvec[0] = 24; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; break;
		case 3 :
			numrad = 5; rvec[0] = 24; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 4 :
			numrad = 5; rvec[0] = 12; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
		case 5 :
			numrad = 5; rvec[0] =  6; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 16; break;
		case 6 :
			numrad = 5; rvec[0] =  6; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; break;
	  #endif
		default :
			if(nradices){*nradices = 4 + (3 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 1664 :					/* 1.625M = 1664K */
		switch(radix_set)
		{
		case 0 :
			numrad = 4; rvec[0] = 52; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; break;
		case 1 :
			numrad = 4; rvec[0] = 52; rvec[1] = 32; rvec[2] = 32; rvec[3] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 2 :
			numrad = 4; rvec[0] = 26; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; break;
		case 3 :
			numrad = 5; rvec[0] = 26; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
		case 4 :
			numrad = 5; rvec[0] = 13; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
	  #endif
		default :
			if(nradices){*nradices = 2 + (3 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 1792 :					/* 1.75M = 1792K */
		switch(radix_set)
		{
		case 0 :
			numrad = 4; rvec[0] = 56; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; break;
		case 1 :
			numrad = 4; rvec[0] = 56; rvec[1] = 32; rvec[2] = 32; rvec[3] = 16; break;
		case 2 :
			numrad = 4; rvec[0] = 28; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; break;
		case 3 :
			numrad = 5; rvec[0] = 28; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 4 :
			numrad = 5; rvec[0] = 14; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
		case 5 :
			numrad = 5; rvec[0] =  7; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 16; break;
		case 6 :
			numrad = 5; rvec[0] =  7; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; break;
	  #endif
		default :
			if(nradices){*nradices = 4 + (3 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 1920 :					/* 1.875M = 1920K */
		switch(radix_set)
		{
		case 0 :
			numrad = 4; rvec[0] = 60; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; break;
		case 1 :
			numrad = 4; rvec[0] = 60; rvec[1] = 32; rvec[2] = 32; rvec[3] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 2 :
			numrad = 4; rvec[0] = 30; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; break;
		case 3 :
			numrad = 5; rvec[0] = 30; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
		case 4 :
			numrad = 5; rvec[0] = 15; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
	  #endif
		default :
			if(nradices){*nradices = 2 + (3 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 2048 :					/* 2M = 2048K */
		switch(radix_set)
		{
		case 0 :
			numrad = 4; rvec[0] = 64; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; break;
		case 1 :
			numrad = 4; rvec[0] = 32; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; break;
		case 2 :
			numrad = 5; rvec[0] = 32; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
		case 3 :
			numrad = 5; rvec[0] = 16; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 4 :
			numrad = 5; rvec[0] =  8; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 16; break;
		case 5 :
			numrad = 5; rvec[0] =  8; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; break;
	  #endif
		default :
			if(nradices){*nradices = 4 + (2 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 2304 :					/* 2.25M = 2304K */
		switch(radix_set)
		{
		case 0 :
			numrad = 4; rvec[0] = 36; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; break;
		case 1 :
			numrad = 5; rvec[0] = 36; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 2 :
			numrad = 5; rvec[0] = 18; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
		case 3 :
			numrad = 5; rvec[0] =  9; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 16; break;
		case 4 :
			numrad = 5; rvec[0] =  9; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; break;
	  #endif
		default :
			if(nradices){*nradices = 2 + (3 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 2560 :					/* 2.5M = 2560K */
		switch(radix_set)
		{
		case 0 :
			numrad = 4; rvec[0] = 40; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; break;
		case 1 :
			numrad = 5; rvec[0] = 40; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
		case 2 :
			numrad = 5; rvec[0] = 20; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 3 :
			numrad = 5; rvec[0] = 10; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 16; break;
		case 4 :
			numrad = 5; rvec[0] = 10; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; break;
	  #endif
		default :
			if(nradices){*nradices = 3 + (2 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
/*
		case 5 :
			numrad = 4; rvec[0] = 80; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; break;
		case 6 :
			numrad = 4; rvec[0] = 80; rvec[1] = 32; rvec[2] = 16; rvec[3] = 32; break;
		case 7 :
			numrad = 4; rvec[0] = 80; rvec[1] = 32; rvec[2] = 32; rvec[3] = 16; break;
	  #endif
		default :
			if(nradices){*nradices = 8 + ( & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
*/
		}; break;
	case 2816 :					/* 2.75M = 2816K */
		switch(radix_set)
		{
		case 0 :
			numrad = 4; rvec[0] = 44; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; break;
		case 1 :
			numrad = 5; rvec[0] = 44; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 2 :
			numrad = 5; rvec[0] = 22; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
		case 3 :
			numrad = 5; rvec[0] = 11; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 16; break;
		case 4 :
			numrad = 5; rvec[0] = 11; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; break;
	  #endif
		default :
			if(nradices){*nradices = 2 + (3 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 3072 :					/* 3M = 3072K */
		switch(radix_set)
		{
		case 0 :
			numrad = 4; rvec[0] = 48; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; break;
		case 1 :
			numrad = 5; rvec[0] = 48; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
		case 2 :
			numrad = 5; rvec[0] = 24; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 3 :
			numrad = 5; rvec[0] = 12; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 16; break;
		case 4 :
			numrad = 5; rvec[0] = 12; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; break;
		case 5 :
			numrad = 6; rvec[0] =  6; rvec[1] =  8; rvec[2] =  8; rvec[3] = 16; rvec[4] = 16; rvec[5] = 16; break;
	  #endif
		default :
			if(nradices){*nradices = 3 + (3 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 3328 :					/* 3.25M = 3328K */
		switch(radix_set)
		{
		case 0 :
			numrad = 4; rvec[0] = 52; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; break;
		case 1 :
			numrad = 5; rvec[0] = 52; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 2 :
			numrad = 5; rvec[0] = 26; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
		case 3 :
			numrad = 5; rvec[0] = 13; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 16; break;
		case 4 :
			numrad = 5; rvec[0] = 13; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; break;
	  #endif
		default :
			if(nradices){*nradices = 2 + (3 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 3584 :					/* 3.5M = 3584K */
		switch(radix_set)
		{
		case 0 :
			numrad = 4; rvec[0] = 56; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; break;
		case 1 :
			numrad = 5; rvec[0] = 56; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
		case 2 :
			numrad = 5; rvec[0] = 28; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 3 :
			numrad = 5; rvec[0] = 14; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 16; break;
		case 4 :
			numrad = 5; rvec[0] = 14; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; break;
		case 5 :
			numrad = 6; rvec[0] =  7; rvec[1] =  8; rvec[2] =  8; rvec[3] = 16; rvec[4] = 16; rvec[5] = 16; break;
	  #endif
		default :
			if(nradices){*nradices = 3 + (3 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 3840 :					/* 3.75M = 3840K */
		switch(radix_set)
		{
		case 0 :
			numrad = 4; rvec[0] = 60; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; break;
		case 1 :
			numrad = 5; rvec[0] = 60; rvec[1] =  8; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 2 :
			numrad = 5; rvec[0] = 30; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
		case 3 :
			numrad = 5; rvec[0] = 15; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 16; break;
		case 4 :
			numrad = 5; rvec[0] = 15; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; break;
	  #endif
		default :
			if(nradices){*nradices = 2 + (3 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 4096 :						/* 4M = 4096K */
		switch(radix_set)
		{
		case 0 :
			numrad = 4; rvec[0] = 64; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; break;
		case 1 :
			numrad = 5; rvec[0] = 32; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
		case 2 :
			numrad = 5; rvec[0] = 16; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; break;
		case 3 :
			numrad = 5; rvec[0] = 16; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 4 :
			numrad = 5; rvec[0] =  8; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 16; break;
		case 5 :
			numrad = 5; rvec[0] =  8; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; break;
	  #endif
		default :
			if(nradices){*nradices = 4 + (2 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 4608 :					/* 4.5M = 4608K */
		switch(radix_set)
		{
		case 0 :
			numrad = 5; rvec[0] = 36; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 1 :
			numrad = 5; rvec[0] = 18; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 16; break;
		case 2 :
			numrad = 5; rvec[0] = 18; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; break;
		case 3 :
			numrad = 5; rvec[0] =  9; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 16; break;
		case 4 :
			numrad = 5; rvec[0] =  9; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; break;
	  #endif
		default :
			if(nradices){*nradices = 1 + (4 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 5120 :					/* 5M = 5120K */
		switch(radix_set)
		{
		case 0 :
			numrad = 5; rvec[0] = 40; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
		case 1 :
			numrad = 5; rvec[0] = 20; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; break;
		case 2 :
			numrad = 5; rvec[0] = 20; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 3 :
			numrad = 5; rvec[0] = 10; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 16; break;
		case 4 :
			numrad = 5; rvec[0] = 10; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; break;
	  #endif
		default :
			if(nradices){*nradices = 3 + (2 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
/*
		case 7 :
			numrad = 4; rvec[0] = 80; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; break;
	  #endif
		default :
			if(nradices){*nradices = 8 + ( & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
*/
		}; break;
	case 5632 :					/* 5.5M = 5632K */
		switch(radix_set)
		{
		case 0 :
			numrad = 5; rvec[0] = 44; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 1 :
			numrad = 5; rvec[0] = 22; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 16; break;
		case 2 :
			numrad = 5; rvec[0] = 22; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; break;
		case 3 :
			numrad = 5; rvec[0] = 11; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 16; break;
		case 4 :
			numrad = 5; rvec[0] = 11; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; break;
	  #endif
		default :
			if(nradices){*nradices = 1 + (4 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 6144 :					/* 6M = 6144K */
		switch(radix_set)
		{
		case 0 :
			numrad = 5; rvec[0] = 48; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
		case 1 :
			numrad = 5; rvec[0] = 24; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; break;
		case 2 :
			numrad = 5; rvec[0] = 24; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 3 :
			numrad = 5; rvec[0] = 12; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 16; break;
		case 4 :
			numrad = 5; rvec[0] = 12; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; break;
		case 5 :
			numrad = 5; rvec[0] =  6; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; rvec[4] = 16; break;
		case 6 :
			numrad = 5; rvec[0] =  6; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
	  #endif
		default :
			if(nradices){*nradices = 3 + (4 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 6656 :					/* 6.5M = 6656K */
		switch(radix_set)
		{
		case 0 :
			numrad = 5; rvec[0] = 52; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 1 :
			numrad = 5; rvec[0] = 26; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 16; break;
		case 2 :
			numrad = 5; rvec[0] = 26; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; break;
		case 3 :
			numrad = 5; rvec[0] = 13; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 16; break;
		case 4 :
			numrad = 5; rvec[0] = 13; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; break;
	  #endif
		default :
			if(nradices){*nradices = 1 + (4 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 7168 :					/* 7M = 7168K */
		switch(radix_set)
		{
		case 0 :
			numrad = 5; rvec[0] = 56; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
		case 1 :
			numrad = 5; rvec[0] = 28; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; break;
		case 2 :
			numrad = 5; rvec[0] = 28; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 3 :
			numrad = 5; rvec[0] = 14; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 16; break;
		case 4 :
			numrad = 5; rvec[0] = 14; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; break;
		case 5 :
			numrad = 5; rvec[0] =  7; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; rvec[4] = 16; break;
		case 6 :
			numrad = 5; rvec[0] =  7; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
	  #endif
		default :
			if(nradices){*nradices = 3 + (4 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 7680 :					/* 7.5M = 7680K */
		switch(radix_set)
		{
		case 0 :
			numrad = 5; rvec[0] = 60; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 1 :
			numrad = 5; rvec[0] = 30; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 16; break;
		case 2 :
			numrad = 5; rvec[0] = 30; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; break;
		case 3 :
			numrad = 5; rvec[0] = 15; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 16; break;
		case 4 :
			numrad = 5; rvec[0] = 15; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; break;
	  #endif
		default :
			if(nradices){*nradices = 1 + (4 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 8192 :						/* 8M = 8192K */
		switch(radix_set)
		{
		case 0 :
			numrad = 5; rvec[0] = 64; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
		case 1 :
			numrad = 5; rvec[0] = 32; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; break;
		case 2 :
			numrad = 5; rvec[0] = 32; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 16; break;
		case 3 :
			numrad = 5; rvec[0] = 16; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; break;
		case 4 :
			numrad = 5; rvec[0] = 16; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 5 :
			numrad = 5; rvec[0] =  8; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; rvec[4] = 16; break;
		case 6 :
			numrad = 5; rvec[0] =  8; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
	  #endif
		default :
			if(nradices){*nradices = 5 + (2 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 9216 :					/* 9M = 9216K */
		switch(radix_set)
		{
		case 0 :
			numrad = 5; rvec[0] = 36; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; break;
		case 1 :
			numrad = 5; rvec[0] = 36; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 2 :
			numrad = 5; rvec[0] = 18; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 16; break;
		case 3 :
			numrad = 5; rvec[0] = 18; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; break;
		case 4 :
			numrad = 5; rvec[0] = 10; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; rvec[4] = 16; break;
		case 5 :
			numrad = 5; rvec[0] = 10; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
	  #endif
		default :
			if(nradices){*nradices = 2 + (4 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 10240 :				/* 10M = 10240K */
		switch(radix_set)
		{
		case 0 :
			numrad = 5; rvec[0] = 40; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; break;
		case 1 :
			numrad = 5; rvec[0] = 40; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 16; break;
		case 2 :
			numrad = 5; rvec[0] = 20; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; break;
		case 3 :
			numrad = 5; rvec[0] = 20; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 4 :
			numrad = 5; rvec[0] = 10; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; rvec[4] = 16; break;
		case 5 :
			numrad = 5; rvec[0] = 10; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
	  #endif
		default :
			if(nradices){*nradices = 4 + (2 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
/*
		case 7 :
			numrad = 5; rvec[0] = 80; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; break;
	  #endif
		default :
			if(nradices){*nradices =10 + ( & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
*/
		}; break;
	case 11264 :				/* 11M = 11264K */
		switch(radix_set)
		{
		case 0 :
			numrad = 5; rvec[0] = 44; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; break;
		case 1 :
			numrad = 5; rvec[0] = 44; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 2 :
			numrad = 5; rvec[0] = 22; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 16; break;
		case 3 :
			numrad = 5; rvec[0] = 22; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; break;
		case 4 :
			numrad = 5; rvec[0] = 11; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; rvec[4] = 16; break;
		case 5 :
			numrad = 5; rvec[0] = 11; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
	  #endif
		default :
			if(nradices){*nradices = 2 + (4 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 12288 :				/* 12M = 12288K */
		switch(radix_set)
		{
		case 0 :
			numrad = 5; rvec[0] = 48; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; break;
		case 1 :
			numrad = 5; rvec[0] = 48; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 16; break;
		case 2 :
			numrad = 5; rvec[0] = 24; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; break;
		case 3 :
			numrad = 5; rvec[0] = 24; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 4 :
			numrad = 5; rvec[0] = 12; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; rvec[4] = 16; break;
		case 5 :
			numrad = 5; rvec[0] = 12; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		case 6 :
			numrad = 5; rvec[0] =  6; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		case 7 :
			numrad = 6; rvec[0] =  6; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 16; break;
	  #endif
		default :
			if(nradices){*nradices = 4 + (4 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 13312 :				/* 13M = 13312K */
		switch(radix_set)
		{
		case 0 :
			numrad = 5; rvec[0] = 52; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; break;
		case 1 :
			numrad = 5; rvec[0] = 52; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 2 :
			numrad = 5; rvec[0] = 26; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 16; break;
		case 3 :
			numrad = 5; rvec[0] = 26; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; break;
		case 4 :
			numrad = 5; rvec[0] = 13; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; rvec[4] = 16; break;
		case 5 :
			numrad = 5; rvec[0] = 13; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
	  #endif
		default :
			if(nradices){*nradices = 2 + (4 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 14336 :				/* 14M = 14336K */
		switch(radix_set)
		{
		case 0 :
			numrad = 5; rvec[0] = 56; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; break;
		case 1 :
			numrad = 5; rvec[0] = 56; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 16; break;
		case 2 :
			numrad = 5; rvec[0] = 28; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; break;
		case 3 :
			numrad = 5; rvec[0] = 28; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 4 :
			numrad = 5; rvec[0] = 14; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; rvec[4] = 16; break;
		case 5 :
			numrad = 5; rvec[0] = 14; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		case 6 :
			numrad = 5; rvec[0] =  7; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		case 7 :
			numrad = 6; rvec[0] =  7; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 16; break;
	  #endif
		default :
			if(nradices){*nradices = 4 + (4 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 15360 :				/* 15M = 15360K */
		switch(radix_set)
		{
		case 0 :
			numrad = 5; rvec[0] = 60; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; break;
		case 1 :
			numrad = 5; rvec[0] = 60; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 2 :
			numrad = 5; rvec[0] = 30; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 16; break;
		case 3 :
			numrad = 5; rvec[0] = 30; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; break;
		case 4 :
			numrad = 5; rvec[0] = 15; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; rvec[4] = 16; break;
		case 5 :
			numrad = 5; rvec[0] = 15; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
	  #endif
		default :
			if(nradices){*nradices = 2 + (4 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 16384 :				/* 16M = 16384K */
		switch(radix_set)
		{
		case 0 :
			numrad = 5; rvec[0] = 64; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; break;
		case 1 :
			numrad = 5; rvec[0] = 32; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; break;
		case 2 :
			numrad = 5; rvec[0] = 32; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 16; break;
		case 3 :
			numrad = 5; rvec[0] = 16; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		case 4 :
			numrad = 5; rvec[0] = 16; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; rvec[4] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 5 :
			numrad = 5; rvec[0] =  8; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		case 6 :
			numrad = 6; rvec[0] =  8; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 16; break;
	  #endif
		default :
			if(nradices){*nradices = 5 + (2 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 18432 :				/* 18M = 18432K */
		switch(radix_set)
		{
		case 0 :
			numrad = 5; rvec[0] = 36; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; break;
		case 1 :
			numrad = 5; rvec[0] = 36; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 2 :
			numrad = 5; rvec[0] = 18; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; rvec[4] = 16; break;
		case 3 :
			numrad = 5; rvec[0] = 18; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		case 4 :
			numrad = 5; rvec[0] =  9; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		case 5 :
			numrad = 6; rvec[0] =  9; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 16; break;
	  #endif
		default :
			if(nradices){*nradices = 2 + (4 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 20480 :				/* 20M = 20480K */
		switch(radix_set)
		{
		case 0 :
			numrad = 5; rvec[0] = 40; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; break;
		case 1 :
			numrad = 5; rvec[0] = 40; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 16; break;
		case 2 :
			numrad = 5; rvec[0] = 20; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; rvec[4] = 16; break;
		case 3 :
			numrad = 5; rvec[0] = 20; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 4 :
			numrad = 5; rvec[0] = 10; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		case 5 :
			numrad = 6; rvec[0] = 10; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 16; break;
	  #endif
		default :
			if(nradices){*nradices = 4 + (2 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 22528 :				/* 22M = 22528K */
		switch(radix_set)
		{
		case 0 :
			numrad = 5; rvec[0] = 44; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; break;
		case 1 :
			numrad = 5; rvec[0] = 44; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 2 :
			numrad = 5; rvec[0] = 22; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; rvec[4] = 16; break;
		case 3 :
			numrad = 5; rvec[0] = 22; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		case 4 :
			numrad = 5; rvec[0] = 11; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		case 5 :
			numrad = 6; rvec[0] = 11; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 16; break;
	  #endif
		default :
			if(nradices){*nradices = 2 + (4 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 24576 :				/* 24M = 24576K */
		switch(radix_set)
		{
		case 0 :
			numrad = 5; rvec[0] = 48; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; break;
		case 1 :
			numrad = 5; rvec[0] = 48; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 16; break;
		case 2 :
			numrad = 5; rvec[0] = 24; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		case 3 :
			numrad = 5; rvec[0] = 24; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; rvec[4] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 4 :
			numrad = 5; rvec[0] = 12; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		case 5 :
			numrad = 6; rvec[0] = 12; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 16; break;
		case 6 :
			numrad = 6; rvec[0] =  6; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; rvec[5] = 16; break;
		case 7 :
			numrad = 6; rvec[0] =  6; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 32; break;
	  #endif
		default :
			if(nradices){*nradices = 4 + (4 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 26624 :				/* 26M = 26624K */
		switch(radix_set)
		{
		case 0 :
			numrad = 5; rvec[0] = 52; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; break;
		case 1 :
			numrad = 5; rvec[0] = 52; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 2 :
			numrad = 5; rvec[0] = 26; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; rvec[4] = 16; break;
		case 3 :
			numrad = 5; rvec[0] = 26; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		case 4 :
			numrad = 5; rvec[0] = 13; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		case 5 :
			numrad = 6; rvec[0] = 13; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 16; break;
	  #endif
		default :
			if(nradices){*nradices = 2 + (4 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 28672 :				/* 28M = 28672K */
		switch(radix_set)
		{
		case 0 :
			numrad = 5; rvec[0] = 56; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; break;
		case 1 :
			numrad = 5; rvec[0] = 56; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 16; break;
		case 2 :
			numrad = 5; rvec[0] = 28; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		case 3 :
			numrad = 5; rvec[0] = 28; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; rvec[4] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 4 :
			numrad = 5; rvec[0] = 14; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		case 5 :
			numrad = 6; rvec[0] = 14; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 16; break;
		case 6 :
			numrad = 6; rvec[0] =  7; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; rvec[5] = 16; break;
		case 7 :
			numrad = 6; rvec[0] =  7; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 32; break;
	  #endif
		default :
			if(nradices){*nradices = 4 + (4 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 30720 :				/* 30M = 30720K */
		switch(radix_set)
		{
		case 0 :
			numrad = 5; rvec[0] = 60; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; break;
		case 1 :
			numrad = 5; rvec[0] = 60; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 2 :
			numrad = 5; rvec[0] = 30; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; rvec[4] = 16; break;
		case 3 :
			numrad = 5; rvec[0] = 30; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		case 4 :
			numrad = 5; rvec[0] = 15; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		case 5 :
			numrad = 6; rvec[0] = 15; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 16; break;
	  #endif
		default :
			if(nradices){*nradices = 2 + (4 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 32768 :				/* 32M = 32768K */
		switch(radix_set)
		{
		case 0 :
			numrad = 5; rvec[0] = 64; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; break;
		case 1 :
			numrad = 5; rvec[0] = 32; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		case 2 :
			numrad = 5; rvec[0] = 32; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; rvec[4] = 16; break;
		case 3 :
			numrad = 5; rvec[0] = 16; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		case 4 :
			numrad = 6; rvec[0] = 16; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 5 :
			numrad = 6; rvec[0] =  8; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; rvec[5] = 16; break;
		case 6 :
			numrad = 6; rvec[0] =  8; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 32; break;
	  #endif
		default :
			if(nradices){*nradices = 5 + (2 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 36864 :				/* 36M = 36864K */
		switch(radix_set)
		{
		case 0 :
			numrad = 5; rvec[0] = 36; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		case 1 :
			numrad = 5; rvec[0] = 36; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; rvec[4] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 2 :
			numrad = 5; rvec[0] = 18; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		case 3 :
			numrad = 6; rvec[0] = 18; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 16; break;
		case 4 :
			numrad = 6; rvec[0] =  9; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; rvec[5] = 16; break;
		case 5 :
			numrad = 6; rvec[0] =  9; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 32; break;
	  #endif
		default :
			if(nradices){*nradices = 2 + (4 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 40960 :				/* 40M = 40960K */
		switch(radix_set)
		{
		case 0 :
			numrad = 5; rvec[0] = 40; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		case 1 :
			numrad = 5; rvec[0] = 40; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; rvec[4] = 16; break;
		case 2 :
			numrad = 5; rvec[0] = 20; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		case 3 :
			numrad = 6; rvec[0] = 20; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 4 :
			numrad = 6; rvec[0] = 10; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; rvec[5] = 16; break;
		case 5 :
			numrad = 6; rvec[0] = 10; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 32; break;
	  #endif
		default :
			if(nradices){*nradices = 4 + (2 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 45056 :				/* 44M = 45056K */
		switch(radix_set)
		{
		case 0 :
			numrad = 5; rvec[0] = 44; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		case 1 :
			numrad = 5; rvec[0] = 44; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; rvec[4] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 2 :
			numrad = 5; rvec[0] = 22; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		case 3 :
			numrad = 6; rvec[0] = 22; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 16; break;
		case 4 :
			numrad = 6; rvec[0] = 11; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; rvec[5] = 16; break;
		case 5 :
			numrad = 6; rvec[0] = 11; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 32; break;
	  #endif
		default :
			if(nradices){*nradices = 2 + (4 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 49152 :				/* 48M = 49152K */
		switch(radix_set)
		{
		case 0 :
			numrad = 5; rvec[0] = 48; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		case 1 :
			numrad = 5; rvec[0] = 48; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; rvec[4] = 16; break;
		case 2 :
			numrad = 5; rvec[0] = 24; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		case 3 :
			numrad = 6; rvec[0] = 24; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 4 :
			numrad = 6; rvec[0] = 12; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; rvec[5] = 16; break;
		case 5 :
			numrad = 6; rvec[0] = 12; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 32; break;
		case 6 :
			numrad = 6; rvec[0] =  6; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; rvec[5] = 16; break;
		case 7 :
			numrad = 6; rvec[0] =  6; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; rvec[5] = 32; break;
	  #endif
		default :
			if(nradices){*nradices = 4 + (4 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 53248 :				/* 52M = 53248K */
		switch(radix_set)
		{
		case 0 :
			numrad = 5; rvec[0] = 52; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		case 1 :
			numrad = 5; rvec[0] = 52; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; rvec[4] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 2 :
			numrad = 5; rvec[0] = 26; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		case 3 :
			numrad = 6; rvec[0] = 26; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 16; break;
		case 4 :
			numrad = 6; rvec[0] = 13; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; rvec[5] = 16; break;
		case 5 :
			numrad = 6; rvec[0] = 13; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 32; break;
	  #endif
		default :
			if(nradices){*nradices = 2 + (4 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 57344 :				/* 56M = 57344K */
		switch(radix_set)
		{
		case 0 :
			numrad = 5; rvec[0] = 56; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		case 1 :
			numrad = 5; rvec[0] = 56; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; rvec[4] = 16; break;
		case 2 :
			numrad = 5; rvec[0] = 28; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		case 3 :
			numrad = 6; rvec[0] = 28; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 4 :
			numrad = 6; rvec[0] = 14; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; rvec[5] = 16; break;
		case 5 :
			numrad = 6; rvec[0] = 14; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 32; break;
		case 6 :
			numrad = 6; rvec[0] =  7; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; rvec[5] = 16; break;
		case 7 :
			numrad = 6; rvec[0] =  7; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; rvec[5] = 32; break;
	  #endif
		default :
			if(nradices){*nradices = 4 + (4 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 61440 :				/* 60M = 61440K */
		switch(radix_set)
		{
		case 0 :
			numrad = 5; rvec[0] = 60; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		case 1 :
			numrad = 5; rvec[0] = 60; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; rvec[4] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 2 :
			numrad = 5; rvec[0] = 30; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		case 3 :
			numrad = 6; rvec[0] = 30; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 16; break;
		case 4 :
			numrad = 6; rvec[0] = 15; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; rvec[5] = 16; break;
		case 5 :
			numrad = 6; rvec[0] = 15; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 32; break;
	  #endif
		default :
			if(nradices){*nradices = 2 + (4 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 65536 :				/* 64M = 65536K */
		switch(radix_set)
		{
		case 0 :
			numrad = 5; rvec[0] = 64; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		case 1 :
			numrad = 5; rvec[0] = 32; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		case 2 :
			numrad = 6; rvec[0] = 32; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 16; break;
		case 3 :
			numrad = 6; rvec[0] = 16; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 32; break;
		case 4 :
			numrad = 6; rvec[0] = 16; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; rvec[5] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 5 :
			numrad = 6; rvec[0] =  8; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; rvec[5] = 16; break;
		case 6 :
			numrad = 6; rvec[0] =  8; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; rvec[5] = 32; break;
	  #endif
		default :
			if(nradices){*nradices = 5 + (2 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 73728 :				/* 72M = 73728K */
		switch(radix_set)
		{
		case 0 :
			numrad = 5; rvec[0] = 36; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		case 1 :
			numrad = 6; rvec[0] = 36; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 2 :
			numrad = 6; rvec[0] = 18; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; rvec[5] = 16; break;
		case 3 :
			numrad = 6; rvec[0] = 18; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 32; break;
		case 4 :
			numrad = 6; rvec[0] =  9; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; rvec[5] = 16; break;
		case 5 :
			numrad = 6; rvec[0] =  9; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; rvec[5] = 32; break;
	  #endif
		default :
			if(nradices){*nradices = 2 + (4 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 81920 :				/* 80M = 81920K */
		switch(radix_set)
		{
		case 0 :
			numrad = 5; rvec[0] = 40; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		case 1 :
			numrad = 6; rvec[0] = 40; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 16; break;
		case 2 :
			numrad = 6; rvec[0] = 20; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 32; break;
		case 3 :
			numrad = 6; rvec[0] = 20; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; rvec[5] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 4 :
			numrad = 6; rvec[0] = 10; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; rvec[5] = 16; break;
		case 5 :
			numrad = 6; rvec[0] = 10; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; rvec[5] = 32; break;
	  #endif
		default :
			if(nradices){*nradices = 4 + (2 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 90112 :				/* 88M = 90112K */
		switch(radix_set)
		{
		case 0 :
			numrad = 5; rvec[0] = 44; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		case 1 :
			numrad = 6; rvec[0] = 44; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 2 :
			numrad = 6; rvec[0] = 22; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; rvec[5] = 16; break;
		case 3 :
			numrad = 6; rvec[0] = 22; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 32; break;
		case 4 :
			numrad = 6; rvec[0] = 11; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; rvec[5] = 16; break;
		case 5 :
			numrad = 6; rvec[0] = 11; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; rvec[5] = 32; break;
	  #endif
		default :
			if(nradices){*nradices = 2 + (4 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 98304 :				/* 96M = 98304K */
		switch(radix_set)
		{
		case 0 :
			numrad = 5; rvec[0] = 48; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		case 1 :
			numrad = 6; rvec[0] = 48; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 16; break;
		case 2 :
			numrad = 6; rvec[0] = 24; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 32; break;
		case 3 :
			numrad = 6; rvec[0] = 24; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; rvec[5] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 4 :
			numrad = 6; rvec[0] = 12; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; rvec[5] = 16; break;
		case 5 :
			numrad = 6; rvec[0] = 12; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; rvec[5] = 32; break;
		case 6 :
			numrad = 6; rvec[0] =  6; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; rvec[5] = 16; break;
		case 7 :
			numrad = 6; rvec[0] =  6; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; rvec[5] = 32; break;
	  #endif
		default :
			if(nradices){*nradices = 4 + (4 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 106496 :				/* 104M = 106496K */
		switch(radix_set)
		{
		case 0 :
			numrad = 5; rvec[0] = 52; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		case 1 :
			numrad = 6; rvec[0] = 52; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 2 :
			numrad = 6; rvec[0] = 26; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; rvec[5] = 16; break;
		case 3 :
			numrad = 6; rvec[0] = 26; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 32; break;
		case 4 :
			numrad = 6; rvec[0] = 13; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; rvec[5] = 16; break;
		case 5 :
			numrad = 6; rvec[0] = 13; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; rvec[5] = 32; break;
	  #endif
		default :
			if(nradices){*nradices = 2 + (4 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 114688 :				/* 112M = 114688K */
		switch(radix_set)
		{
		case 0 :
			numrad = 5; rvec[0] = 56; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		case 1 :
			numrad = 6; rvec[0] = 56; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 16; break;
		case 2 :
			numrad = 6; rvec[0] = 28; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 32; break;
		case 3 :
			numrad = 6; rvec[0] = 28; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; rvec[5] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 4 :
			numrad = 6; rvec[0] = 14; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; rvec[5] = 16; break;
		case 5 :
			numrad = 6; rvec[0] = 14; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; rvec[5] = 32; break;
		case 6 :
			numrad = 6; rvec[0] =  7; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; rvec[5] = 16; break;
		case 7 :
			numrad = 6; rvec[0] =  7; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; rvec[5] = 32; break;
	  #endif
		default :
			if(nradices){*nradices = 4 + (4 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 122880 :				/* 120M = 122880K */
		switch(radix_set)
		{
		case 0 :
			numrad = 5; rvec[0] = 60; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		case 1 :
			numrad = 6; rvec[0] = 60; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 2 :
			numrad = 6; rvec[0] = 30; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; rvec[5] = 16; break;
		case 3 :
			numrad = 6; rvec[0] = 30; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 32; break;
		case 4 :
			numrad = 6; rvec[0] = 15; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; rvec[5] = 16; break;
		case 5 :
			numrad = 6; rvec[0] = 15; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; rvec[5] = 32; break;
	  #endif
		default :
			if(nradices){*nradices = 2 + (4 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 131072 :				/* 128M = 131072K */
		switch(radix_set)
		{
		case 0 :
			numrad = 5; rvec[0] = 64; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; break;
		case 1 :
			numrad = 6; rvec[0] = 32; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 32; break;
		case 2 :
			numrad = 6; rvec[0] = 32; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; rvec[5] = 16; break;
		case 3 :
			numrad = 6; rvec[0] = 16; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; rvec[5] = 32; break;
		case 4 :
			numrad = 6; rvec[0] = 16; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; rvec[5] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 5 :
			numrad = 6; rvec[0] =  8; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; rvec[5] = 16; break;
		case 6 :
			numrad = 6; rvec[0] =  8; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; rvec[5] = 32; break;
	  #endif
		default :
			if(nradices){*nradices = 5 + (2 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 147456 :				/* 144M = 147456K */
		switch(radix_set)
		{
		case 0 :
			numrad = 6; rvec[0] = 36; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 32; break;
		case 1 :
			numrad = 6; rvec[0] = 36; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; rvec[5] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 2 :
			numrad = 6; rvec[0] = 18; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; rvec[5] = 16; break;
		case 3 :
			numrad = 6; rvec[0] = 18; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; rvec[5] = 32; break;
		case 4 :
			numrad = 6; rvec[0] =  9; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; rvec[5] = 16; break;
		case 5 :
			numrad = 6; rvec[0] =  9; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; rvec[5] = 32; break;
	  #endif
		default :
			if(nradices){*nradices = 2 + (4 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 163840 :				/* 160M = 163840K */
		switch(radix_set)
		{
		case 0 :
			numrad = 6; rvec[0] = 40; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 32; break;
		case 1 :
			numrad = 6; rvec[0] = 40; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; rvec[5] = 16; break;
		case 2 :
			numrad = 6; rvec[0] = 20; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; rvec[5] = 32; break;
		case 3 :
			numrad = 6; rvec[0] = 20; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; rvec[5] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 4 :
			numrad = 6; rvec[0] = 10; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; rvec[5] = 16; break;
		case 5 :
			numrad = 6; rvec[0] = 10; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; rvec[5] = 32; break;
	  #endif
		default :
			if(nradices){*nradices = 4 + (2 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 180224 :				/* 176M = 180224K */
		switch(radix_set)
		{
		case 0 :
			numrad = 6; rvec[0] = 44; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 32; break;
		case 1 :
			numrad = 6; rvec[0] = 44; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; rvec[5] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 2 :
			numrad = 6; rvec[0] = 22; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; rvec[5] = 16; break;
		case 3 :
			numrad = 6; rvec[0] = 22; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; rvec[5] = 32; break;
		case 4 :
			numrad = 6; rvec[0] = 11; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; rvec[5] = 16; break;
		case 5 :
			numrad = 6; rvec[0] = 11; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; rvec[5] = 32; break;
	  #endif
		default :
			if(nradices){*nradices = 2 + (4 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 196608 :				/* 192M = 196608K */
		switch(radix_set)
		{
		case 0 :
			numrad = 6; rvec[0] = 48; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 32; break;
		case 1 :
			numrad = 6; rvec[0] = 48; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; rvec[5] = 16; break;
		case 2 :
			numrad = 6; rvec[0] = 24; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; rvec[5] = 32; break;
		case 3 :
			numrad = 6; rvec[0] = 24; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; rvec[5] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 4 :
			numrad = 6; rvec[0] = 12; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; rvec[5] = 16; break;
		case 5 :
			numrad = 6; rvec[0] = 12; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; rvec[5] = 32; break;
		case 6 :
			numrad = 6; rvec[0] =  6; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; rvec[5] = 16; break;
		case 7 :
			numrad = 6; rvec[0] =  6; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; rvec[5] = 32; break;
	  #endif
		default :
			if(nradices){*nradices = 4 + (4 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 212992 :				/* 208M = 212992K */
		switch(radix_set)
		{
		case 0 :
			numrad = 6; rvec[0] = 52; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 32; break;
		case 1 :
			numrad = 6; rvec[0] = 52; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; rvec[5] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 2 :
			numrad = 6; rvec[0] = 26; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; rvec[5] = 16; break;
		case 3 :
			numrad = 6; rvec[0] = 26; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; rvec[5] = 32; break;
		case 4 :
			numrad = 6; rvec[0] = 13; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; rvec[5] = 16; break;
		case 5 :
			numrad = 6; rvec[0] = 13; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; rvec[5] = 32; break;
	  #endif
		default :
			if(nradices){*nradices = 2 + (4 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 229376 :				/* 224M = 229376K */
		switch(radix_set)
		{
		case 0 :
			numrad = 6; rvec[0] = 56; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 32; break;
		case 1 :
			numrad = 6; rvec[0] = 56; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; rvec[5] = 16; break;
		case 2 :
			numrad = 6; rvec[0] = 28; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; rvec[5] = 32; break;
		case 3 :
			numrad = 6; rvec[0] = 28; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; rvec[5] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 4 :
			numrad = 6; rvec[0] = 14; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; rvec[5] = 16; break;
		case 5 :
			numrad = 6; rvec[0] = 14; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; rvec[5] = 32; break;
		case 6 :
			numrad = 6; rvec[0] =  7; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; rvec[5] = 16; break;
		case 7 :
			numrad = 6; rvec[0] =  7; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; rvec[5] = 32; break;
	  #endif
		default :
			if(nradices){*nradices = 4 + (4 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 245760 :				/* 240M = 245760K */
		switch(radix_set)
		{
		case 0 :
			numrad = 6; rvec[0] = 60; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 32; break;
		case 1 :
			numrad = 6; rvec[0] = 60; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; rvec[5] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 2 :
			numrad = 6; rvec[0] = 30; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; rvec[5] = 16; break;
		case 3 :
			numrad = 6; rvec[0] = 30; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; rvec[5] = 32; break;
		case 4 :
			numrad = 6; rvec[0] = 15; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; rvec[5] = 16; break;
		case 5 :
			numrad = 6; rvec[0] = 15; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; rvec[5] = 32; break;
	  #endif
		default :
			if(nradices){*nradices = 2 + (4 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	case 262144 :				/* 256M = 262144K */
		switch(radix_set)
		{
		case 0 :
			numrad = 6; rvec[0] = 64; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 16; rvec[5] = 32; break;
		case 1 :
			numrad = 6; rvec[0] = 32; rvec[1] = 16; rvec[2] = 16; rvec[3] = 16; rvec[4] = 32; rvec[5] = 32; break;
		case 2 :
			numrad = 6; rvec[0] = 32; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; rvec[5] = 16; break;
		case 3 :
			numrad = 6; rvec[0] = 16; rvec[1] = 16; rvec[2] = 16; rvec[3] = 32; rvec[4] = 32; rvec[5] = 32; break;
		case 4 :
			numrad = 6; rvec[0] = 16; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; rvec[5] = 16; break;
	  #ifndef USE_ONLY_LARGE_LEAD_RADICES
		case 5 :
			numrad = 6; rvec[0] =  8; rvec[1] = 32; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; rvec[5] = 16; break;
		case 6 :
			numrad = 6; rvec[0] =  8; rvec[1] = 16; rvec[2] = 32; rvec[3] = 32; rvec[4] = 32; rvec[5] = 32; break;
	  #endif
		default :
			if(nradices){*nradices = 5 + (2 & nrad_scalar_mask);}	return ERR_RADIXSET_UNAVAILABLE;
		}; break;
	default :
		 return ERR_FFTLENGTH_ILLEGAL;
	}

	/* Check that there are at least 3 radices: */
	if(numrad < 3)
	{
		fprintf(stderr,"ERROR: get_fft_radices: %d #radices must be >= 3.\n",numrad);
		return ERR_FFTLENGTH_ILLEGAL;
	}

	/* If user provided a radix array, make sure they gave a valid dimension: */
	if(radix_vec)
		ASSERT(HERE, radix_vec_dim >=	numrad,"get_fft_radices: radix_vec_dim has illegal value!");

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
		ASSERT(HERE, 0,"0");
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
	int nradices;

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
				ASSERT(HERE, radset != 0, "test_fft_radixtables: Should only see ERR_RADIXSET_UNAVAILABLE for nonzero radix set index!");
				break;
			}
			else if(retval == ERR_FFTLENGTH_ILLEGAL)
			{
				fprintf(stderr,"ERROR: illegal FFT length %u K in test_fft_radixtables self-test!\n",kblocks);
				ASSERT(HERE, 0,"0");
			}
			else
			{
				fprintf(stderr,"ERROR: unknown return value %d in test_fft_radixtables self-test; i = %d, kblocks = %u, radset = %u.\n", retval, i, kblocks, radset);
				ASSERT(HERE, 0,"0");
			}
		}
		++i;
	}
}

/*
!...Set vector length, based on number of bits (p) in numbers to be FFT-multiplied.
*/
uint32 get_default_fft_length(uint32 p)
{
	int nradices;
	uint32 leadingRadixVec[N_LEADING_RADICES] = {8,9,10,11,12,13,14,15};
	uint32 i, twoK, fftLen;

	ASSERT(HERE, PMAX > PMIN,"get_default_fft_length: PMAX > PMIN");
	if(p < PMIN || p > PMAX)
	{
		fprintf(stderr,"get_default_fft_length: invalid value for exponent %u\n",p);
		ASSERT(HERE, 0,"0");
		return 0;
	}

	/* Starting with N = 1K, Loop over all FFT lengths of form {8,9,10,11,12,13,14,15}*2^m,
	and return the smallest one for which maxP >= p: */
	i = 0;
	ASSERT(HERE, 1024%leadingRadixVec[i] == 0,"get_default_fft_length: 1024%leadingRadixVec[0] == 0");
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
	ASSERT(HERE, 0,"0");
	return 0;
}

/*
!...Given a legal FFT length, return the next-larger supported length.
*/
uint32 get_nextlarger_fft_length(uint32 n)
{
	uint32 fftLen, lead4, rem2;	/* Leading 4 bits of input FFT length and number of bits in the remaining power of 2 */

	if(get_fft_radices((n >> 10), 0, 0x0, 0x0, 0) != 0)
	{
		sprintf(cbuf, "get_nextlarger_fft_length: Illegal or Unsupported input FFT length %u\n", n);
		ASSERT(HERE, 0, cbuf);
	}

	/* Extract leading 4 bits of input FFT lengths, thus decomposing it into the form {8,9,10,11,12,13,14,15}*2^m,
	(m stored in rem2 variable here, decompose as n = lead4*2^rem2) then increment the leading portion -
	don't need to maintain leading part normalized in [8,15] here:
	*/
	rem2 = 32 - leadz32(n) - 4;
	lead4 = n >> rem2;
	ASSERT(HERE, lead4 > 7 && lead4 < 16,"get_nextlarger_fft_length: leading 4 bits of input FFT length out of range!");

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
uint32 given_N_get_maxP(uint32 N)
{
	const double Bmant = 53;
	const double AsympConst = 0.6;
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

	/* 3/10/05: Future versions will need to loosen this p < 2^32 restriction: */
	ASSERT(HERE, maxExp2 <= 1.0*0xffffffff,"given_N_get_maxP: maxExp2 <= 1.0*0xffffffff");
/*
	fprintf(stderr,"N = %8u K  maxP = %10u\n", N>>10, (uint32)maxExp2);
*/
	return (uint32)maxExp2;
}
/* Here a simple PARI 'script' to return maxExp for any desired

N = [enter FFT length in #doubles], can use [FFTlen in K]<<10 here
Bmant = 53.; AsympConst = 0.6; ln2inv = 1.0/log(2.0);
ln_N = log(1.0*N); lnln_N = log(ln_N); l2_N = ln2inv*ln_N; lnl2_N = log(l2_N); l2l2_N = ln2inv*lnl2_N; lnlnln_N = log(lnln_N); l2lnln_N = ln2inv*lnlnln_N;
Wbits = 0.5*( Bmant - AsympConst - 0.5*(l2_N + l2l2_N) - 1.5*(l2lnln_N) )
maxExp2 = Wbits*N

With AsympConst = 0.6, this gives the following maxP values for various FFT lengths:

             N          maxP
        --------  ----------
             1 K       22686
             2 K       44683
             3 K       66435
             4 K       88029
             5 K      109506
             6 K      130892
             7 K      152201
             8 K      173445
             9 K      194632
            10 K      215769
            12 K      257912
            14 K      299904
            16 K      341769
            18 K      383521
            20 K      425174
            24 K      508222
            28 K      590972
            32 K      673470
            36 K      755746
            40 K      837827
            48 K     1001477
            56 K     1164540
            64 K     1327103
            72 K     1489228
            80 K     1650966
            96 K     1973430
           112 K     2294732
           128 K     2615043
           144 K     2934488
           160 K     3253166
           176 K     3571154
           192 K     3888516
           208 K     4205305
           224 K     4521565
           240 K     4837335
           256 K     5152648
           288 K     5782016
           320 K     6409862
           352 K     7036339
           384 K     7661575
           416 K     8285675
           448 K     8908726
           480 K     9530805
           512 K    10151977
           576 K    11391823
           640 K    12628648
           704 K    13862759
           768 K    15094405
           832 K    16323795
           896 K    17551103
           960 K    18776481
          1024 K    20000058
          1152 K    22442252
          1280 K    24878447
          1408 K    27309250
          1536 K    29735157
          1664 K    32156582
          1792 K    34573872
          1920 K    36987325
          2048 K    39397201
          2304 K    44207097
          2560 K    49005071
          2816 K    53792328
          3072 K    58569855
          3328 K    63338470
          3584 K    68098867
          3840 K    72851637
          4096 K    77597294
          4608 K    87069012
          5120 K    96517023
          5632 K   105943724
          6144 K   115351074
          6656 K   124740700
          7168 K   134113980
          7680 K   143472090
  8 M =   8192 K   152816052
  9 M =   9216 K   171464992
 10 M =  10240 K   190066770
 11 M =  11264 K   208626152
 12 M =  12288 K   227147031
 13 M =  13312 K   245632644
 14 M =  14336 K   264085729
 15 M =  15360 K   282508628
 16 M =  16384 K   300903371
 18 M =  18432 K   337615274
 20 M =  20480 K   374233313
 22 M =  22528 K   410766968
 24 M =  24576 K   447223981
 26 M =  26624 K   483610796
 28 M =  28672 K   519932856
 30 M =  30720 K   556194824
 32 M =  32768 K   592400738
 36 M =  36864 K   664658102
 40 M =  40960 K   736728582
 44 M =  45056 K   808631042
 48 M =  49152 K   880380890
 52 M =  53248 K   951990950
 56 M =  57344 K  1023472059
 60 M =  61440 K  1094833496
 64 M =  65536 K  1166083299
 72 M =  73728 K  1308275271
 80 M =  81920 K  1450095024
 88 M =  90112 K  1591580114
 96 M =  98304 K  1732761219
104 M = 106496 K  1873663870
112 M = 114688 K  2014309644
120 M = 122880 K  2154717020
128 M = 131072 K  2294902000
144 M = 147456 K  2574659086
160 M = 163840 K  2853674592
176 M = 180224 K  3132023315
192 M = 196608 K  3409766353	<<<<*** Allows numbers slightly > 1Gdigit
208 M = 212992 K  3686954556
224 M = 229376 K  3963630903
240 M = 245760 K  4239832202	<*** largest FFT length for 32-bit exponents
256 M = 262144 K  4515590327
288 M = 294912 K  5065885246
320 M = 327680 K  5614702299
352 M = 360448 K  6162190494
384 M = 393216 K  6708471554
416 M = 425984 K  7253646785
448 M = 458752 K  7797801823
480 M = 491520 K  8341010002
512 M = 524288 K  8883334834

For F33 have 2^33 = 8589934592, so could try a special FFT length based on leading radix = 31,
which should be just large enough to run without overly frequent fatal ROE > 0.4:

496 M = 507904 K	8612279197

But as radix-31 DFT is relatively inefficient in terms of flops/input, better to use the slightly
larger but much soother odd radix 63:

504 M = 516096 K
*/

