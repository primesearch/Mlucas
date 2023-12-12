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

//================================ GPU stuff: ==================================
#ifdef USE_GPU

#include "factor.h"
#include "align.h"

#define BIT_SET_8(_arr,_bit) {\
	uint32 n = _bit>>3;\
	*((uint8*)_arr+n) |=  ((uint8)1 << (_bit&7));\
}

#define BIT_CLR_8(_arr,_bit) {\
	uint32 n = _bit>>3;\
	*((uint8*)_arr+n) &= ~((uint8)1 << (_bit&7));\
}

#define BIT_CLR32(_arr,_bit) {\
	uint32 n = _bit>>5;\
	_arr[n] &= ~((uint32)1 << (_bit&31));\
}

// Inline to sub a nonnegative constant (del) mod p.  That is given i and del in [0, p-1] , return ((i - del) % p)
__device__ __inline static int decr_mod_p(int i, int del, int p)
{
	int	x, j;
	i = i - del; j = i + p;
	asm("slct.s32.s32 %0, %1, %2, %1;" : "=r" (x) : "r" (i), "r" (j));
	return x;
}

// Double-word left-shift of a pair of 32-bit ints x,y; x has high 32 bits of the pair. No arg-checking done; assumes n in (0,31]:
// 	NOTE! User must avoid calling with n = 0 due to undef'd behavior of >> (32-0)
#define DSHIFT32_L(x,y,n)	((uint32)(x) << n) + ((uint32)(y) >> (32-n))

// ...And the 64-bit version:
#define DSHIFT64_L(x,y,n)	((uint64)(x) << n) + ((uint64)(y) >> (64-n))

// Init 'live'-mem bitmap using (ncopies_of_master) master-copies, which is assumed 32-bit aligned with no padding of any kind:
__global__ void refreshBitmap(const uint32*bit_map1, uint32*bit_map2, const uint32 nthr, const uint32 nwords, const uint32 ncopies_of_master)
{
	uint32 idx = blockDim.x * blockIdx.x + threadIdx.x;	// Hardware thread ID:
	if(idx < nthr) {
		uint32 i, nc;
		const uint32 *src_ptr  = bit_map1 + idx*nwords;	// Starting word of current thread's source-bitmap chunk
		      uint32 *dest_ptr = bit_map2 + idx*nwords;	// Starting word of current thread's destination-bitmap chunk
		for(nc = 0; nc < ncopies_of_master; nc++) {
			for(i = 0; i < nwords; i++) {
				*(dest_ptr + i) = *(src_ptr + i);
			}
			dest_ptr += nthr*nwords;	// There are (nthr*nwords) words in master (bit_map1), so incr dest-idx by that much on each copy-pass
		}
	}
}

/******************************************************************************************/
/* 2-pass (plus small serial intermediate partial-sums-propagation step) procedure for    */
/* parsing small-primes-cleared bitmap and converting each 1-bit to a factor-candidate k: */
/******************************************************************************************/

#define SERIAL	// This is the default, at least until I impl a 1024-thread version of the ||-accum algo below

// All-in-one version of the 3-function sequence below:
__global__ void parseBitmap_InitKoff3in1(uint32*bit_map, uint32*kvec, const uint32 nthr, const uint32 nwords)
{
	__shared__ uint32 popc_vec[256];
	uint32 idx = blockDim.x * blockIdx.x + threadIdx.x;	// Hardware thread ID:
	// Pass 1: Treating bitmap as (say) 256 equal-size chunks of [nwords] 32-bit words each,
	// accumulate #-set-bits in each chunk and deposit the resulting chunk-sum in an array:
	uint32 i, numk = 0;
	uint32 ioff = idx*nwords;
	uint32 *bm_ptr = bit_map + ioff;
	for(i = 0; i < nwords; i++) {
		numk += __popc(*(bm_ptr+i));
	}
	popc_vec[idx] = numk;

	// Intermediate partial-sums-propagation step. Two ways of thinking about this:
	// First five tallies remain within one warp.  Should be in lock-step.
	if (threadIdx.x & 1)        // If we are running on any thread 0bxxxxxxx1, tally neighbor's count.
	popc_vec[threadIdx.x] += popc_vec[threadIdx.x - 1];

	if (threadIdx.x & 2)        // If we are running on any thread 0bxxxxxx1x, tally neighbor's count.
	popc_vec[threadIdx.x] += popc_vec[(threadIdx.x - 2) | 1];

	if (threadIdx.x & 4)        // If we are running on any thread 0bxxxxx1xx, tally neighbor's count.
	popc_vec[threadIdx.x] += popc_vec[(threadIdx.x - 4) | 3];

	if (threadIdx.x & 8)        // If we are running on any thread 0bxxxx1xxx, tally neighbor's count.
	popc_vec[threadIdx.x] += popc_vec[(threadIdx.x - 8) | 7];

	if (threadIdx.x & 16)       // If we are running on any thread 0bxxx1xxxx, tally neighbor's count.
	popc_vec[threadIdx.x] += popc_vec[(threadIdx.x - 16) | 15];

	// Further tallies are across warps.  Must synchronize
	__syncthreads();
	if (threadIdx.x  & 32)      // If we are running on any thread 0bxx1xxxxx, tally neighbor's count.
	popc_vec[threadIdx.x] += popc_vec[(threadIdx.x - 32) | 31];

	__syncthreads();
	if (threadIdx.x & 64)       // If we are running on any thread 0bx1xxxxxx, tally neighbor's count.
	popc_vec[threadIdx.x] += popc_vec[(threadIdx.x - 64) | 63];

	__syncthreads();
	if (threadIdx.x & 128)       // If we are running on any thread 0b1xxxxxxx, tally neighbor's count.
	popc_vec[threadIdx.x] += popc_vec[127];

	// At this point, popc_vec[...] contains the total number of bits for the indexed
	// thread plus all lower-numbered threads.  I.e., popc_vec[255] is the total count.
	__syncthreads();

	numk = popc_vec[255];	// This also serves as the return mechanism for numk to the calling __host__ routine
//if(idx == 0)printf("parseBitmap_InitKoff3in1 gives bitmapPop = %u\n",numk);
	bit_map[0] = numk;		// Since done with current copy of bit_map, use0-elt to stash numk
	// Zero the 'one-beyond' elt of the bitmap - on return from the ||-modpow,
	// this will hold any factor-k found in the current interval:
	kvec[numk] = kvec[numk+1] = 0x0;	// Need 2 zero-elts in 32-bit-offset mode

	// Pass 2: Treating bitmap in same decomposed manner as Pass 1, compute the k-offsets for the set bits in each chunk and
	// store them in an array, using the array-index offsets computed in the above intermediate partial-sums-propagation step:
	uint32 j, wd, idx_mask = -(int32)(idx > 0);
	numk = 0;
	// Each elt of popc_vec[] store #set-bits of the current chunk and the
	// ones below it, thus serves as the kvec-offset for the NEXT-HIGHER block.
	// First AND here avoids out-of-range -1 array index; 2nd makes offset = 0 for idx = 0.
	uint32 *kv_ptr = kvec + ( (popc_vec[(idx-1) & idx_mask]) & idx_mask );
	ioff = 60*(ioff << 5);
	for(i = 0; i < nwords; i++) {
		wd = *(bm_ptr+i);
	#if 0
		for(j = 0; j < 32; j++, ioff += 60) {
			if(((wd >> j) & 1) != 0) {
				*(kv_ptr + numk++) = ioff;
			}
		}
	#elif 1	// 2-bits-at-a-time ~10% faster than 1-bit (4-bits is much slower than either 1 or 2 due to explosion of the logic tree):
		for(j = 0; j < 32; j += 2, ioff += 120) {
			uint32 two_bits = (wd >> j) & 3;
			if(two_bits == 1) {
				*(kv_ptr + numk++) = ioff   ;
			} else if(two_bits == 2) {
				*(kv_ptr + numk++) = ioff+60;
			} else if(two_bits == 3) {
				*(kv_ptr + numk++) = ioff   ;
				*(kv_ptr + numk++) = ioff+60;
			}
		}
	#endif
	}
}


// Pass 1: Treating bitmap as (say) 256 equal-size chunks of [nwords] 32-bit words each,
// accumulate #-set-bits in each chunk and deposit the resulting chunk-sum in an array:
__global__ void parseBitmap_CountSetBits(const uint32*bit_map, uint32*popc_vec, const uint32 nthr, const uint32 nwords)
{
	uint32 idx = blockDim.x * blockIdx.x + threadIdx.x;	// Hardware thread ID:
/*	if(idx == 0) {
		popc_vec[0] = 0;
	} else */
	if(idx < nthr) {
		uint32 i, numk = 0;
//		const uint32 *bm_ptr = bit_map + (idx-1)*nwords;
		const uint32 *bm_ptr = bit_map + idx*nwords;
		for(i = 0; i < nwords; i++) {
			numk += __popc(*(bm_ptr+i));
		}
		popc_vec[idx] = numk;
	}
}

/* Intermediate partial-sums-propagation step. Two ways of thinking about this:

1. [lazy bastard serialist] "This step is alas strictly serial, but consists of just (nthr) integer adds and array writes";

2. [bright-eyed and bushy-tailed parallelist] Parallel partial-sums algo -
	We desire a partial-sum of an n-vec, that is given input vector v[0],...,v[n-1],
	to compute an ovec w[] such that w[i] = Sum_0^i ( v[j] ).
	Example:
	In:
		v[8] = 3, 1, 4, 5, 9, 2, 6, 2
	Out
		w[8] = 3, 4, 8,13,22,24,30,32

	How to parallelize this? Here is a simple divide-and-conquer scheme:

	Step 1: group into pairs, replace hi-idx elt of each pair by the sum of te 2 elts in the pair:

		| 3, 1|,| 4, 5|,| 9, 2|,| 6, 2| ==> | 3, 4|,| 4, 9|,| 9,11|,| 6, 8|

	Step 2: group into quartets (viewed as pairs of pairs), add hi-idx elt of the lo pair of each quartet to each elt in the hi pair:

		| 3, 4|,| 4, 9|	| 9,11|,| 6, 8| ==> | 3, 4|,| 8,13|	| 9,11|,|17,19|

	Step 2: group into octets (viewed as pairs of quartets), add hi-idx elt of the lo quartet of each octet to each elt in the hi quartet:

		| 3, 4, 8,13|,| 9,11,17,19| ==> | 3, 4, 8,13|,|22,24,30,32| .

	Et voila! The generalization to larger power-of-2 input lenghts is obvious. Only downsides are

		(a) Requirement of power-of-2 vector length (we currently require length 256,
			and enforce this by suitably adjusting the contiguous-chunk wordcounts in the 2 bookending routines), and
		(b) Halving of the potential parallelism at each step.

Basic algo was easy enough for me to reinvent on my own, but the sweet nVidia-optimized implementation below is due to Rocke Verser:

In practice there is very little speedup (< 1%) from using the || algo, but it's a good small example of "thinking parallelly" for the nVidia arch.
*/
__global__ void parseBitmap_PopcountSums(                      uint32*kvec, uint32*popc_vec, const uint32 nthr)
{
	uint32 numk;
#ifdef SERIAL
//	#error Must use parallel algo!
	uint32 idx = blockDim.x * blockIdx.x + threadIdx.x;	// Hardware thread ID:
	if(idx == 0) {
		uint32 i = 0;
//printf("parseBitmap_PopcountSums: bitmap sub-blocks have popcounts = \n%u [%u]",popc_vec[i],i);
		for(i = 1; i < nthr; i++) {	// Skip 0-elt
//printf(", %u [%u]",popc_vec[i],i);
			popc_vec[i] += popc_vec[i-1];
		}
		numk = popc_vec[nthr-1];	// This also serves as the return mechanism for numk to the calling __host__ routine
//printf("parseBitmap_PopcountSums [SERIAL] gives bitmapPop = %u\n",numk);
	}

#else
	// First five tallies remain within one warp.  Should be in lock-step.
	if (threadIdx.x & 1)        // If we are running on any thread 0bxxxxxxx1, tally neighbor's count.
	popc_vec[threadIdx.x] += popc_vec[threadIdx.x - 1];

	if (threadIdx.x & 2)        // If we are running on any thread 0bxxxxxx1x, tally neighbor's count.
	popc_vec[threadIdx.x] += popc_vec[(threadIdx.x - 2) | 1];

	if (threadIdx.x & 4)        // If we are running on any thread 0bxxxxx1xx, tally neighbor's count.
	popc_vec[threadIdx.x] += popc_vec[(threadIdx.x - 4) | 3];

	if (threadIdx.x & 8)        // If we are running on any thread 0bxxxx1xxx, tally neighbor's count.
	popc_vec[threadIdx.x] += popc_vec[(threadIdx.x - 8) | 7];

	if (threadIdx.x & 16)       // If we are running on any thread 0bxxx1xxxx, tally neighbor's count.
	popc_vec[threadIdx.x] += popc_vec[(threadIdx.x - 16) | 15];

	// Further tallies are across warps.  Must synchronize
	__syncthreads();
	if (threadIdx.x  & 32)      // If we are running on any thread 0bxx1xxxxx, tally neighbor's count.
	popc_vec[threadIdx.x] += popc_vec[(threadIdx.x - 32) | 31];

	__syncthreads();
	if (threadIdx.x & 64)       // If we are running on any thread 0bx1xxxxxx, tally neighbor's count.
	popc_vec[threadIdx.x] += popc_vec[(threadIdx.x - 64) | 63];

	__syncthreads();
	if (threadIdx.x & 128)       // If we are running on any thread 0b1xxxxxxx, tally neighbor's count.
	popc_vec[threadIdx.x] += popc_vec[127];

	// At this point, popc_vec[...] contains the total number of bits for the indexed
	// thread plus all lower-numbered threads.  I.e., popc_vec[255] is the total count.
	__syncthreads();

	numk = popc_vec[255];	// This also serves as the return mechanism for numk to the calling __host__ routine
//(idx == 0)printf("parseBitmap_PopcountSums [PARALLEL] gives bitmapPop = %u\n",numk);

#endif

	// Zero the 'one-beyond' elt of the bitmap - on return from the ||-modpow,
	// this will hold any factor-k found in the current interval:
	kvec[numk] = kvec[numk+1] = 0x0;	// Need 2 zero-elts in 32-bit-offset mode

/* Debug:
	if(idx == 0) {
		uint32 i;
		printf("partial sums:");
		for(i = 0; i < nthr; i++) {
			printf("%5u..",popc_vec[i]);
			if((i&15) == 0)printf("\n");
		}
		printf("final sum = %u\n",numk);
	}
*/
}

// Pass 2: Treating bitmap in same decomposed manner as Pass 1, compute the k-offsets for the set bits in each chunk and
// store them in an array, using the array-index offsets computed in the above intermediate partial-sums-propagation step:
__global__ void parseBitmap_InitKoffsets(const uint32*bit_map, uint32*kvec, const uint32*popc_vec, const uint32 nthr, const uint32 nwords)
{
	uint32 idx = blockDim.x * blockIdx.x + threadIdx.x;	// Hardware thread ID:
	if(idx < nthr) {
		uint32 ioff = idx*nwords;
		uint32 i, j, wd, numk = 0, idx_mask = -(int32)(idx > 0);
		// Each elt of popc_vec[] store #set-bits of the current chunk and the
		// ones below it, thus serves as the kvec-offset for the NEXT-HIGHER block.
		// First AND here avoids out-of-range -1 array index; 2nd makes offset = 0 for idx = 0.
		uint32 *kv_ptr = kvec + ( (popc_vec[(idx-1) & idx_mask]) & idx_mask );
		const uint32 *bm_ptr = bit_map + ioff;
		ioff = 60*(ioff << 5);
		for(i = 0; i < nwords; i++) {
			wd = *(bm_ptr+i);
		#if 0
			for(j = 0; j < 32; j++, ioff += 60) {
				if(((wd >> j) & 1) != 0) {
					*(kv_ptr + numk++) = ioff;
				}
			}
		#elif 1	// 2-bits-at-a-time ~10% faster than 1-bit:
			for(j = 0; j < 32; j += 2, ioff += 120) {
				uint32 two_bits = (wd >> j) & 3;
				if(two_bits == 1) {
					*(kv_ptr + numk++) = ioff   ;
				} else if(two_bits == 2) {
					*(kv_ptr + numk++) = ioff+60;
				} else if(two_bits == 3) {
					*(kv_ptr + numk++) = ioff   ;
					*(kv_ptr + numk++) = ioff+60;
				}
			}
		#endif
		}
	}
}

/*
__global__ void retrieveBitmap(uint32*bit_map1, uint32*bit_map2, const uint32 nthr, const uint32 nwords)
{
	uint32 idx = blockDim.x * blockIdx.x + threadIdx.x;	// Hardware thread ID:
	if(idx < nthr) {
		uint32 i;
		uint32 *ptr_out = bit_map1 + idx*nwords, *ptr_in = bit_map2 + idx*nwords;	// Starting word or current thread's bit_map chunk
		for(i = 0; i < nwords; i++) {
			*(ptr_out + i) = *(ptr_in + i);
		}
	}
}
*/

// GPU-|| function to update startval array storing first bit 'hit' by each odd prime in bitmap encoding the current k-interval.
// Caller must precompute array of every-100th-prime (starting with p_last_small, presumed to be in 0-slot of said array)
// and thus having the 0-thread only process from that to p(100)) and supply via the p100 pointer.
//
__global__ void updateStartval(
	const uint32 bit_len,
	const uint32 nthr,
	const uint32 nclear,
	const uint32*p100,
	const uint8 *pdiff,
	      uint32*startval)
{
	uint32 idx = blockDim.x * blockIdx.x + threadIdx.x;	// Hardware thread ID:
	uint32 curr_p,i,ilo,ihi,l,m;
	if(idx < nthr)
	{
		ilo = MAX(nclear, idx*100);	// 0-thread uses ilo = nclear rather than 0; all others use ilo = idx*100
		ihi = (idx+1)*100;		// NB: For some reason, using ilo = MAX(nclear, idx*100) caused NVCC to make the if() body here into a no-op.
		curr_p = p100[idx];
//(idx == 0)printf("Idx = %2u: ilo,ihi,curr_p = %u,%u,%u, sval0 = %u\n",idx,ilo,ihi,curr_p,startval[ilo]);
		for(m = ilo; m < ihi; m++)
		{
			curr_p += (pdiff[m] << 1);
			// Here is a quick way - most expensive op is the modulo - to compute the update to startval[m] each time through the sievelet:
			l = startval[m];								// l = old startval
			i = bit_len % curr_p;
			i = (l - i) + (curr_p & -((int32)(l - i) < 0));	// i = new startval
			startval[m] = i;
//(m == ilo)printf("curr_p[%u] = %u, sval0 = %u, sval1 = %u\n",m,curr_p,l,startval[m]);
		}
	}
}

// Device-parallel function in which each GPU thread takes 64-bits of our current GPU sieve (k-interval) and clears bits
// corr. to as-yet-uncleared (that is, not built-into-sieve) primes < 64, i.e. primes 19,23,29,31,37,41,43,47,53,59,61.
// Assumes caller provides pointer to proper element of startval array, containing "first hit" bit offsets in the
// sieve of each odd prime greater than the 'built-in' set [3,5,7,11,13,17], i.e. *(startval+0) gives to first bit
// corresponding to a factor candidate divisible by 19:
//
#if 0	// Original version with const 0-shift initial bitmasks:

	__global__ void clrPrimes_Lt64(uint64*bit_map, const uint32 nthr, const pinv_t*pinv, const uint32*startval)
	{
		const uint8 curr_p[11] = {19,23,29,31,37,41,43,47,53,59,61};
		const uint64 maskL[11] = {
			~0x0200004000080001ull,	// 1+2^19+2^38+2^57
			~0x0000400000800001ull,	// 1+2^23+2^46
			~0x0400000020000001ull,	// 1+2^29+2^58
			~0x4000000080000001ull,	// 1+2^31+2^62
			~0x0000002000000001ull,	// 1+2^37
			~0x0000020000000001ull,	// 1+2^41
			~0x0000080000000001ull,	// 1+2^43
			~0x0000800000000001ull,	// 1+2^47
			~0x0020000000000001ull,	// 1+2^53
			~0x0800000000000001ull,	// 1+2^59
			~0x2000000000000001ull	// 1+2^61
		};
		// For each left-shift of a small-prime-associated mask, the elements of this array provide the bits to shift on from the right.
		const uint64 maskR[11] = {	// Negative exponents alias to their (mod 64) values:
			~0x0000200004000080ull,	// 2^-19 + 2^-38 + 2^-57 = 2^45 + 2^26 + 2^7
			~0x0000020000040000ull,	// 2^-23 + 2^-46 = 2^41 + 2^18
			~0x0000000800000040ull,	// 2^-29 + 2^-58 = 2^35 + 2^6
			~0x0000000200000004ull,	// 2^-31 + 2^-62 = 2^33 + 2^2
			~0x0000000008000000ull,	// 2^-37 = 2^27
			~0x0000000000800000ull,	// 2^-41 = 2^23
			~0x0000000000200000ull,	// 2^-43 = 2^21
			~0x0000000000020000ull,	// 2^-47 = 2^17
			~0x0000000000000800ull,	// 2^-53 = 2^11
			~0x0000000000000020ull,	// 2^-59 = 2^ 5
			~0x0000000000000008ull	// 2^-61 = 2^ 3
		};

		uint32 idx = blockDim.x * blockIdx.x + threadIdx.x;	// Hardware thread ID:

		// Since this function deals with non-power-of-2 array sizes (and hence threadcounts), it gets fed out-of-range thread indices by
		// the warp-scheduler's insistence on 'filling in' #threads to the next-higher power of 2. Thus must wrap body in explicit range check:
		if(idx < nthr)
		{
			/* Here is a quick way - most expensive op is the modulo - to compute the update to startval[m] each time through the sievelet:
				i = bit_len % curr_p;							// l = old startval
				i = (l - i) + (curr_p & -((int32)(l - i) < 0));	// i = new startval
			In the current instance, bit_len is the number of bits from bit 0 of the full-length bitmap
			for the 64-bit word handled by the current thread, thus 64*idx:
			*/
			uint32 i,j,l;
			uint64 mask = -1;
			for(i = 0; i < 11; i++) {
				// Replacing library mod here (and below) with float-based alternative cuts overall time by ~10%, sieve-time by ~30%:
			#if 0
				j = (idx<<6) % curr_p[i];
			#else
				j = (idx<<6);
				j -= (uint32)(j*pinv[i])*curr_p[i];
			#endif
				l = startval[i];									// l = startval for idx = 0
				j = (l - j) + (curr_p[i] & -((int32)(l - j) < 0));	// j = startval for this thread's idx
			//	ASSERT(HERE, maskL[i] == DSHIFT64_L(maskL[i], maskR[i], curr_p[i]), "MASK bad!");
				// Must avoid 0-count shifts due to undefined behavior of >> (64-0)
				if(j) {
					mask &= DSHIFT64_L(maskL[i], maskR[i], j);
				} else {
					mask &=            maskL[i]              ;
				}
		//		if(idx ==3) { printf("Idx = %2u: p = %2u, sval0 = %2u, sval1 = %2u, mask = %llX, cshft = %llX\n",idx,curr_p[i],l,j,maskL[i],mask); }
			}
			bit_map[idx] &= mask;
		}
	}

#else	// Version in which we compute the 64-bit masks on the fly - only slightly faster, but code much more compact:

	__global__ void clrPrimes_Lt64(uint64*bit_map, const uint32 nthr, const pinv_t*pinv, const uint32*startval)
	{
		const uint8 curr_p[11] = {19,23,29,31,37,41,43,47,53,59,61};
		uint32 idx = blockDim.x * blockIdx.x + threadIdx.x;	// Hardware thread ID:
		if(idx < nthr)
		{
			uint32 i,j,l,p;
			uint64 mask = 0;	// Use a 'negative mask' and invert in final global-mem-wordmask step
			for(i = 0; i < 11; i++) {
				p = curr_p[i];
				// Replacing library mod here (and below) with float-based alternative cuts overall time by ~10%, sieve-time by ~30%:
				j = (idx<<6) - (uint32)((idx*64)*pinv[i])*p;
				l = startval[i];							// l = startval for idx = 0
			#if 0
				j = (l - j) + (p & -((int32)(l - j) < 0));	// j = startval for this thread's idx
			#else
				j = decr_mod_p(l, j, p);
			#endif
				// Loop over occurrences of curr_p-hits in current 64-bit word and OR them in:
				while(j < 64) {
					mask |= (uint64)1 << j;
					j += p;
				}
			}
			bit_map[idx] &= ~mask;
		}
	}

#endif

// Each thread handles four 32-bit bitmap words:
__global__ void clrPrimes_Lt128(uint32*bit_map, const uint32 nthr, const pinv_t*pinv, const uint32*startval)
{
	const uint8 curr_p[24] = {19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,127};
	uint32 idx = blockDim.x * blockIdx.x + threadIdx.x;	// Hardware thread ID:
	if(idx < nthr)
	{
		uint32 i,j,l,p;
		uint32 mask[4] = {0,0,0,0};	// Use a 'negative mask' and invert in final global-mem-wordmask step
		for(i = 0; i < 24; i++) {
			p = curr_p[i];
			// Replacing library mod here (and below) with float-based alternative is faster:
			j = (idx<<7) - (uint32)((idx*128)*pinv[i])*p;
			l = startval[i];							// l = startval for idx = 0
			j = decr_mod_p(l, j, p);
			// OR occurrences of curr_p-hits in current 128-bit mask. Each p hits 1 or 2 bits:
			while(j < 128) {
				BIT_SET_8(mask,j);
				j += p;
			}
		}
		i = idx << 2;
		bit_map[i  ] &= ~mask[0];
		bit_map[i+1] &= ~mask[1];
		bit_map[i+2] &= ~mask[2];
		bit_map[i+3] &= ~mask[3];
	}
}

// Each thread handles four 32-bit bitmap words:
__global__ void clrPrimes_Gt64_Lt128(uint32*bit_map, const uint32 nthr, const pinv_t*pinv, const uint32*startval)
{
	const uint8 curr_p[13] = {67,71,73,79,83,89,97,101,103,107,109,113,127};
	uint32 idx = blockDim.x * blockIdx.x + threadIdx.x;	// Hardware thread ID:
	if(idx < nthr)
	{
		uint32 i,j,l,p;
		uint32 mask[4] = {0,0,0,0};	// Use a 'negative mask' and invert in final global-mem-wordmask step
		for(i = 0; i < 13; i++) {
			p = curr_p[i];
			// Replacing library mod here (and below) with float-based alternative is faster:
			j = (idx<<7) - (uint32)((idx*128)*pinv[i])*p;
			l = startval[i];							// l = startval for idx = 0
			j = decr_mod_p(l, j, p);
			// OR occurrences of curr_p-hits in current 128-bit mask. Each p hits 1 or 2 bits:
			while(j < 128) {
				BIT_SET_8(mask,j);
				j += p;
			}
		}
		i = idx << 2;
		bit_map[i  ] &= ~mask[0];
		bit_map[i+1] &= ~mask[1];
		bit_map[i+2] &= ~mask[2];
		bit_map[i+3] &= ~mask[3];
	}
}
/*
	Restructuring the while(j < 128){} above as a paired if(j < 64) / if(j < 128) burns lotsa regs:
	Before (i.e. as above):
		Used 16 registers, 64 bytes cmem[0], 4 bytes cmem[16]	<*** Same cmem irresp. of whether const curr_p decl uint8 or uint32, but stack frame = 72 bytes in the latter case ***
		32 bytes stack frame, 0 bytes spill stores, 0 bytes spill loads
	Using (paired if()):
			// OR occurrences of curr_p-hits in current 128-bit mask. Each p hits 1 or 2 bits:
			if(j < 64) {
				mask[0] |= (uint64)1 << j;
				j += p;
			}
			if(j < 128) {
				mask[1] |= (uint64)1 << (j-64);
			}
		Used 8 registers, 52 bytes cmem[0]
		40 bytes stack frame, 44 bytes spill stores, 40 bytes spill loads
*/

#if 0	// Next 70 primes are below - we can compactly store these (in fact primes <= 641) in uint8 form as (p - 131)/2:
curr_p[ 6-29] = 19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,127;
curr_p[30-99] =
131,137,139,149,151,157,163,167,173,179,
181,191,193,197,199,211,223,227,229,233,
239,241,251,257,263,269,271,277,281,283,
293,307,311,313,317,331,337,347,349,353,
359,367,373,379,383,389,397,401,409,419,
421,431,433,439,443,449,457,461,463,467,
479,487,491,499,503,509,521,523		// ...,541,547,...

// Each thread handles four 32-bit bitmap words:
__global__ void clrPrimes_Gt128_Lt540(uint32*bit_map, const uint32 nthr, const pinv_t*pinv, const uint32*startval)
{
	const uint8 curr_p[68] = {0,3,4,9,10,13,16,18,21,24,25,30,31,33,34,40,46,48,49,51,54,55,60,63,66,69,70,73,75,76,81,88,90,91,93,100,103,108,109,111,114,118,121,124,126,129,133,135,139,144,145,150,151,154,156,159,163,165,166,168,174,178,180,184,186,189,195,196};
	uint32 idx = blockDim.x * blockIdx.x + threadIdx.x;	// Hardware thread ID:
	if(idx < nthr)
	{
		uint32 i,j,l,p;
		uint32 mask[8] = {0,0,0,0,0,0,0,0};	// Use a 'negative mask' and invert in final global-mem-wordmask step
		for(i = 0; i < 68; i++) {
			p = 131 + ((uint32)curr_p[i]<<1);
			// Replacing library mod here (and below) with float-based alternative is faster:
			j = (idx<<8) - (uint32)((idx*256)*pinv[i])*p;
			l = startval[i];							// l = startval for idx = 0
			j = decr_mod_p(l, j, p);
			// OR occurrences of curr_p-hits in current 128-bit mask. Each p hits 1 or 2 bits:
			while(j < 256) {
				BIT_SET_8(mask,j);
				j += p;
			}
		}
		i = idx << 3;
		bit_map[i  ] &= ~mask[0];
		bit_map[i+1] &= ~mask[1];
		bit_map[i+2] &= ~mask[2];
		bit_map[i+3] &= ~mask[3];
		bit_map[i+4] &= ~mask[4];
		bit_map[i+5] &= ~mask[5];
		bit_map[i+6] &= ~mask[6];
		bit_map[i+7] &= ~mask[7];
	}
}

#endif

// Name here is misleading - this is variant of above which handles all primes >= 19 and < 540:
__global__ void clrPrimes_Gt128_Lt540(uint32*bit_map, const uint32 nthr, const pinv_t*pinv, const uint32*startval)
{
	const uint8 curr_p[92] = {0,2,5,6,9,11,12,14,17,20,21,24,26,27,30,32,35,39,41,42,44,45,47,54,56,59,60,65,66,69,72,74,77,80,81,86,87,89,90,96,102,104,105,107,110,111,116,119,122,125,126,129,131,132,137,144,146,147,149,156,159,164,165,167,170,174,177,180,182,185,189,191,195,200,201,206,207,210,212,215,219,221,222,224,230,234,236,240,242,245,251,252};
	uint32 idx = blockDim.x * blockIdx.x + threadIdx.x;	// Hardware thread ID:
	if(idx < nthr)
	{
		uint32 i,j,l,p;
		uint32 mask[8] = {0,0,0,0,0,0,0,0};	// Use a 'negative mask' and invert in final global-mem-wordmask step
		for(i = 0; i < 92; i++) {
			p = 19 + ((uint32)curr_p[i]<<1);
			// Replacing library mod here (and below) with float-based alternative is faster:
			j = (idx<<8) - (uint32)((idx*256)*pinv[i])*p;
			l = startval[i];							// l = startval for idx = 0
			j = decr_mod_p(l, j, p);
			// OR occurrences of curr_p-hits in current 128-bit mask. Each p hits 1 or 2 bits:
			while(j < 256) {
				BIT_SET_8(mask,j);
				j += p;
			}
		}
		i = idx << 3;
		bit_map[i  ] &= ~mask[0];
		bit_map[i+1] &= ~mask[1];
		bit_map[i+2] &= ~mask[2];
		bit_map[i+3] &= ~mask[3];
		bit_map[i+4] &= ~mask[4];
		bit_map[i+5] &= ~mask[5];
		bit_map[i+6] &= ~mask[6];
		bit_map[i+7] &= ~mask[7];
	}
}

/********** Code ripped off from GW's gpusieve.cu: My added notes marked with 'EWM' ***********/

// Inline to calculate x mod p using 2^32 / p, for unsigned inputs with (obviously) p != 0.
// This version is for x,p < 2^31. For a general 32-bit version see util.c : check_nbits_in_types().

#define gen_pinv(p)	(0xFFFFFFFF / (p))

__device__ __inline static uint32 mod_p (const uint32 x, const uint32 p, const uint32 pinv)
{
	uint32 r;

#if 0	// Generic-C version (with the CUDA __umulhi() intrinsic inlined for the MULH32):

	uint32 y = x;	// Use x-copy in place of x below, since must leave inputs as-is
	r = __MULH32(y, pinv);	//	r = __mulhi (y, pinv);
	r = __MULL32(r, p   );	//	r = r * p;
	y = y - r;
	r = y - p;
	r = (r >= p) ? y : r ;

#elif 0	// For the number ranges of interest to us, the extra subtract-correction is almost never needed,
		// and skipping it at worst causes an occasional bit-which-should-get-cleared to remain lit:
	r = __MULH32(x, pinv);
	r = x + __MULL32(r,-p);	// This needs just a single FMA in an inline-PTX-asm implementation

#else

	// EWM: Inline-asm of above streamlined code sequence:
	asm ("mul.hi.u32	%0, %1, %2;		\n\t"		//	r = __MULH32(x, pinv);
	     "mad.lo.u32	%0, %0, %3, %1;	\n\t"		//	r = x - r * p, restructured as x + r*(-p) to permit FMA use
	     : "=r" (r), "+r" (x) : "r" (pinv), "r" (-p));

#endif
	return r;
}

// Device-parallel function in which each GPU thread takes an (nword)-long contiguous segment of 32-bit words of our sieve
// (each thread handles a disjoint segment such that a covering-set union of the segments, blah blah) and clears bits corr.
// to as-yet-uncleared primes > (A) and <= (B), where N is the (np)th odd prime > curr_p. The input array pdiff contains
// the half-diffs needed to generate the primes in question, starting from the initial prime curr_p, also supplied via arglist.
//
// The caller must pad out said array to the next-larger multiple of nwords and fill in the pads
// with 0s, if its unpadded size is not naturally a multiple of nwords.
//
// Assumes caller provides pointer to proper elt of startval array (containing "first hit" bit offsets in the
// sievelet of each prime, i.e. *(startval+0) corr. to first bit corresponding to a factor candidate divisible by 67.
//
__global__ void clrPrimes_GtA_LtB(uint32*bit_map, const uint32 nthr, const uint32 nwords, const uint32 np, const uint8*pdiff, const pinv_t*pinv, const uint32*startval, const uint32 curr_p)
{
	const uint32 WMAX = 16;
	const uint32 nbits = (nwords << 5);
	uint32 idx = blockDim.x * blockIdx.x + threadIdx.x;	// Hardware thread ID:
	if(idx < nthr) {
		uint32 i,j,l,p = curr_p;	// largest prime pvsly cleared
		uint32 mask[WMAX];	// nwords may vary between calls, so use largest-allowable dimension
		uint32*ptr = bit_map + idx*nwords;	// Starting word or current thread's bit_map chunk

		for(i = 0; i < nwords; i++) {
			mask[i] = 0;
		}
		for(i = 0; i < np; i++)
		{
			p += (pdiff[i] << 1);
			// Replacing library mod here (and below) with float-based alternative cuts overall time by ~25%, sieve-time by ~50%:
			// But, trying other 'obvious optimimizations' like precomputing loop const (idx<<5)*nwords and storing
			// same as float in order to save a type conversion each loop pass actually slowed things down ~20%! Bizarre...
		#ifdef FLOAT_PINV	// Float-inverse version:
			j = (idx<<5)*nwords;		// Fastest float-based option ... note moving (idx<<5)*nwords outside loop and computing just once per thread kills timing(!?)
			j -= (uint32)(j*pinv[i])*p;	// Cast-to-int gives same timing as truncf()
		//	if( j != (((idx<<5)*nwords) % p) ) {
		//		printf("ERROR: Idx = %d: j,p[%u],fpinv = %u,%u,%e: Float-pinv gives j (mod p) = %u, exact = %u\n",idx,i,j,p,pinv[i], j,((idx<<5)*nwords)%p);
		//	}
		#else
			j = (idx<<5)*nwords;
			j = mod_p(j,p,pinv[i]);	// Int-inverse here
		//	if( j != (((idx<<5)*nwords) % p) ) {
		//		printf("ERROR: Idx = %d: j,p[%u],ipinv = %u,%u,%u: Int ipinv gives j (mod p) = %u, exact = %u\n",idx,i,j,p,pinv[i], j,((idx<<5)*nwords)%p);
		//	}
		#endif
			l = startval[i];							// l = startval for idx = 0
			j = decr_mod_p(l, j, p);
			// Mask off the appropriate bit/word of the mask array:
			while(j < nbits) {
				BIT_SET_8(mask,j);
//(idx==0)printf("i = %3u, curr_p = %3u, clear bit %2u, mask = %16llX\n",i,p,j,mask[j >> 5]);
				j += p;
			}
		}
		for(i = 0; i < nwords; i++) {
			*(ptr + i) &= ~mask[i];
		}
	}
__syncthreads();	// 3/23/15: Not necessary, but gives a 2-3% speedup (!?)
}

__global__ void popcountBitmap(uint32*bit_map, const uint32 nwords)
{
	uint32 i, idx = blockDim.x * blockIdx.x + threadIdx.x;	// Hardware thread ID:
	if(idx == 0) {	// This function is debug-diagnostic, hence serial only
		uint32 numk = 0;
//printf("bitmap per-word popcounts = \n");
		for(i = 0; i < nwords; ++i) {
//printf(", %u [%u]",__popc(bit_map[i]),i);
//printf("%u+",__popc(bit_map[i]));
			numk += __popc(bit_map[i]);
//if(i%100==0)printf("\ti = %u: currPop = %u\n",i,numk);
		}
//printf("bitmapPop = %u\n",numk);
	}
}

#define USE_SHARED_MEM	0
#if !USE_SHARED_MEM	// 'else' has shared-mem version (limit 256 threads) which is FUBAR, NEEDS DEBUG.

	__global__ void clrPrimes_LuchaLibre(uint32*bit_map, const uint32 nthr, const uint32 bit_len, const uint32*prim, uint32*startval)
	{
		uint32 l,p, idx = blockDim.x * blockIdx.x + threadIdx.x;	// Hardware thread ID:
		// Each thread grabs one of the first 256 primes above our start threshold, then processes every 256th prime above:
		if(idx < nthr) {
			l = startval[idx];
			p = prim[idx];
		/* Special-handling code for p == curr_p case:
			if(l == 0xffffffff) {
				return;
			}
		*/
			// Need to account for the fact that primes greater than the # of bits in the sieve
			// may not hit *any* sieve bit on the current pass through the sieve:
			while(l < bit_len)
			{
			// Atomic AND avoids thread collisions and the resulting survival of
			// bits-which-should-have-been-cleared, but result is 10-20% slower overall:
			#if 0
				uint32 n = l>>5;	// Atomics require 32-bit ints (or 64-bit)
				atomicAnd(bit_map + n, ~((uint32)1 << (l&31)));
			#else
				BIT_CLR_8(bit_map, l);
			#endif
				l += p;
			}
			// save new startvalue:
			startval[idx] = l-bit_len;
		}
	}

#else

	// Shared-mem requires Larger-array-chunk-per-thread version, as we are limited to 256 threads:
	__global__ void clrPrimes_LuchaLibre(uint32*bit_map, const uint32 np_per_thr, const uint32 bit_len, const uint32*prim, uint32*startval)
	{
		uint32 i,j,l,curr_p, idx = blockDim.x * blockIdx.x + threadIdx.x;	// Hardware thread ID:
	#define LOCAL__	1
	#if LOCAL__
		const uint32 nwords = 8704;	// 2*4352
		uint32 ioff = 34*idx;		// = idx*(nwords/256)
		__shared__ uint32 local_map[nwords];
		uint8 *map_ptr = (uint8*)local_map;
	  #if 0
		// Use first (nwords) threads to copy input global-array bitmap to fast shared memory, 34 x 32-bit words per thread:
		if(idx < nwords) {
			local_map[idx] = 0;
		}
	  #else
		// Use first 256 threads to copy input global-array bitmap to fast shared memory, 34 x 32-bit words per thread:
		if(idx < 256) {
			for(i = 0; i < 34; i++) {
				local_map[i+ioff] = 0;
			}
		}
	  #endif
		__syncthreads();
	#else
		uint32*map_ptr =   bit_map;
	#endif
		// Each thread grabs one of the first 256 primes above our start threshold, then processes every 256th prime above:
		for(i = 0; i < np_per_thr; i++) {
			j = idx + (i << 8);
			curr_p = prim[j];
			l = startval[j];
			// Special-handling code for p == curr_p case:
			if(l == 0xffffffff) {
				return;
			}
			// Need to account for the fact that primes greater than the # of bits in the sieve
			// may not hit *any* sieve bit on the current pass through the sieve:
			if(l >= bit_len) {
				startval[j] -= bit_len;
			} else {
				while(l < bit_len)
				{
				// Atomic AND avoids thread collisions and the resulting survival of
				// bits-which-should-have-been-cleared, but result is 10-20% slower overall:
				#if LOCAL__
					BIT_SET_8(map_ptr,l);
				#else
				  #if 0
					uint32 n = l>>5;	// Atomics require 32-bit ints (or 64-bit)
					atomicAnd((uint32*)map_ptr+n, ~((uint32)1 << (l&31)));
				  #else
					BIT_CLR_8(map_ptr,l);
				  #endif
				#endif
					l += curr_p;
				}
				// save new startvalue:
				startval[j] = l-bit_len;
			}
		}
		__syncthreads();

	// Use first 256 threads to copy local shared-mem bitmap copy back to global memory, 34 x 32-bit words per thread:
	#if LOCAL__
	  #if 0
		// Use first (nwords) threads to transfer mask bits from shared memory to global, 34 x 32-bit words per thread:
		if(idx < nwords) {
			bit_map[idx] &= ~local_map[idx];
		}
	  #else
		// Use first 256 threads to transfer mask bits from shared memory to global, 34 x 32-bit words per thread:
		if(idx < 256) {
			for(i = 0; i < 34; i++) {
				bit_map[i+ioff] &= ~local_map[i+ioff];
			}
		}
	  #endif
	__syncthreads();
	#endif
	}

#endif

__host__ uint64 PerPass_tfSieve(
	const uint64 interval_lo, const uint64 interval_hi,
	const double fbits_in_2p,
	const uint32 nclear,
	const uint32 sieve_len,	// #64-bit words in the full-length sieving bitmap setup by the (CPU-side) caller.
							// This routine uses only a specific (len/60)-sized "sievelet" excerpt.
							// sieve_len = 64*(product of first NCLEAR odd primes)
	const uint32 p_last_small,	//largest odd prime appearing in the product; that is, the (nclear)th odd prime.
	const uint32 NUM_SIEVING_PRIME,	// #sieving primes (counting from 3)
	const uint32 MAX_SIEVING_PRIME,
	const uint8 *pdiff,
	      uint32*startval,	// This gets updated within
		  uint64*k_to_try,	// Unused by GPU code; include to yield a (mostly) uniform API
	      uint64*factor_k,	// List of found factors for each p gets updated (we hope) within
	      uint32*nfactor,	// Here the '*' is to denote a writeable scalar
	const uint32 findex,
	const uint64 p,		// GPU version currently only supports 32-bit p and 64-bit k
	const uint32 incr,
	const uint64 kstart,
	const uint64*bit_map, uint64*bit_map2,	// GPU version uses only 1st of these
	double *tdiff,
	const int MODULUS_TYPE,
	const char*VERSION,
	const char*OFILE
)
{
	static int first_entry = 1;		// For managing one-time allocs; factor.c calls this function once for each (mod 60) k-class
	cudaError_t cudaError = cudaGetLastError();	// Call this to reset error flag to 0
	if(cudaError != cudaSuccess)
	{
		printf("ERROR: cudaGetLastError() returned %d: %s\n", cudaError, cudaGetErrorString(cudaError));
		ASSERT(HERE, 0, "gpu_sieve.cu : Error on entry to PerPass_tfSieve - check upstream code!");
	}
	const uint32 nprime = NUM_SIEVING_PRIME;	// #sieving primes (counting from 3)
	uint32 nthr;		// Encodes # of GPU threads for each device-function call
	// sieve_len = 3*5*7*11*13*17 = 255255, thus bit_len = 16*7*11*13*17 = 272272 bits = 34 kB; ONLY use low 16 bits of last 32-bit word.
//	const uint32 bit_map_len = (sieve_len/60+1);// Number of 64-bit words in each of the 16 mod-60 sievelets = 4255
	const uint32 bit_len = 64*sieve_len/60; 	// Number of         bits in each of the 16 mod-60 sievelets = 272272.
//printf("sieve_len = %u\n",sieve_len);
//printf("bit_map_len,bit_len = %u, %u\n",bit_map_len,bit_len);	exit(0);

	// Do-sieving-on-GPU needs these:
	const uint32 NSIEVE_COPIES_IN_GPU_SIEVE = 4;	// Number of compact (~34kB) sievelet copies copied (sans padding) into the GPU-master.
							// Recall that the CPU-generated sievelet == 16 (mod 64), so take = 4 to ensure 64-bit alignment of the result.
	const uint32 gpu_bit_len = bit_len*NSIEVE_COPIES_IN_GPU_SIEVE;	// #bits in GPU-master. Basic GPU-side bitmap copy is (4 x bit_len) bits,
																	// which is 64-bit aligned, but only require 32-bit alignment.
//	ASSERT(HERE,(gpu_bit_len & 31) == 0, "GPU-master bitmap must have #bits divisible by 32!");
	const uint32 gpu_map_len = gpu_bit_len >> 5;	// GPU-side analog of bit_map_len: Number of 32-bit words in the master GPU-bitmap.
//****** start here, NMASTER_COPIES_IN_LIVE_SIEVE support *************
	const uint32 NMASTER_COPIES_IN_LIVE_SIEVE = 16;	// #Copies of GPU-sieve master (gpu_bit_map) to place into working sieve (gpu_bit_map2.)
/************** To-Do: This is currently limited by the 32-bit k-offsets scheme, since 2^32/60 = 71582788 bits ... 64 here gives nbits = 34034*64*32 = 69701632, which is just under the limit. ***************/
		// The 34kB CPU sieve was designed with a CPU L1 cache in mind, but on GPU, want to minimize frequency of kernel calls, within reason.
	const uint32 pad_map_len = ((gpu_map_len*NMASTER_COPIES_IN_LIVE_SIEVE + 255) & 0xffffff00);
	// get_startval() routine is generalized to multiword-p/q, thus expects a ptr to a two_p array as one of its args:
	static uint64 two_p[1] = { (uint64)p+p };
	const uint32 kvsz = 100000*NSIEVE_COPIES_IN_GPU_SIEVE*NMASTER_COPIES_IN_LIVE_SIEVE;
	// Using first 10k primes get ~50k survivors per CPU-size bitmap, so a multiplier of 100k gives ~2x the expected required-allocation here:
	uint32 curr_p, m;
	uint64 count = 0, k = kstart, kfac = 0, sweep;
	const double ILG2 = 1.0/log(2.0);
	double fbits_in_k = 0, fbits_in_q = 0;
	int threadsPerBlock = 256, blocksPerGrid = 0;	// Fix first of these and adjust 2nd based on needed threadcount for the procedure in question

/***** These pointers are static due to 1-time-init aspect: *****/
// CPU-side malloc-arrays:
	static uint32 *cpu_p100 = 0x0;	// Array of every-100th-prime, starting with p_last_small, presumed to be in 0-slot of said array.
	static uint32 *gpu_p100 = 0x0;
	static pinv_t *cpu_pinv = 0x0;
	static pinv_t *gpu_pinv = 0x0;
	// GPU-side malloc-arrays:
	static uint32 *gpu_kvec = 0x0;	// GPU device copy of kvec.
	static uint32 *gpu_bit_map = 0x0, *gpu_bit_map2 = 0x0;
	static uint8  *gpu_pdiff   = 0x0;
	static uint32 *gpu_startval= 0x0;
	static uint32 *gpu_popc_vec= 0x0;
// 2/28/2015: Larger-primes "lucha libre" step needs explicit primes arrays:
	static uint32 *cpu_prim    = 0x0;
	static uint32 *gpu_prim    = 0x0;

if(first_entry) {
	first_entry = 0;	// FALSE

	// Set Prefs to prefer L1 cache over shared-mem via the Runtime API: For cc2.0,
	// cudaFuncCachePreferShared: shared memory is 48 KB, L1 cache is 16 KB
	// cudaFuncCachePreferL1:     shared memory is 16 KB, L1 cache is 48 KB
	cudaFuncSetCacheConfig(refreshBitmap, cudaFuncCachePreferL1);	// Doing it for the first-called kernel makes it the default for the rest of the run
	cudaError = cudaGetLastError();
	if(cudaError != cudaSuccess)
	{
		printf("ERROR: cudaGetLastError() returned %d: %s\n", cudaError, cudaGetErrorString(cudaError));
		ASSERT(HERE, 0, "PerPass_tfSieve : Call to cudaFuncSetCacheConfig failed!");
	}

	cpu_prim = ALLOC_UINT(cpu_prim, nprime);
	cudaMalloc(&gpu_prim    , (nprime)<<2 );

	// GPU device copy of kvec:
	cudaMalloc(&gpu_kvec, (kvsz << 2));

	// Array of every-100th-prime, starting with p_last_small, presumed to be in 0-slot of said array.
	cpu_p100 = ALLOC_UINT(cpu_p100, nprime/100);
	cudaMalloc(&gpu_p100    , (nprime/100)<<2 );
	// GPU-master bitmap is treated as const once inited from the CPU-generated master:
	cudaMalloc(&gpu_bit_map , gpu_map_len << 2 );
//intf("Dim of GPU master bitmap: %u words\n",gpu_map_len);
	// #words in live bitmap must be multiple of threadsPerBlock (256)
  	cudaMalloc(&gpu_bit_map2, pad_map_len << 2 );
	// Zero the padding words in the live bitmap:
	nthr = pad_map_len - gpu_map_len*NMASTER_COPIES_IN_LIVE_SIEVE;	// Not a thread count here - just use nthr as a convenient place to stash #pads
	for(m = 0; m < nthr; m++) {
		cpu_prim[m] = 0;	// And cpu_prim[] is a convenient place to stash (npad) 0s
	}
//intf("Dim of GPU   live bitmap: %u words, %u of which are 0-pads\n",pad_map_len, nthr);
	cudaMemcpy(gpu_bit_map2 + gpu_map_len*NMASTER_COPIES_IN_LIVE_SIEVE , cpu_prim, nthr << 2, cudaMemcpyHostToDevice);

  	cudaMalloc(&gpu_pdiff   ,  nprime );
	cudaMalloc(&gpu_startval, (nprime)<<2);
#ifdef FLOAT_PINV
	cpu_pinv = ALLOC_FLOAT(cpu_pinv, nprime);	cpu_pinv = ALIGN_FLOAT(cpu_pinv);
	cudaMalloc(&gpu_pinv    , (nprime)<<2);
#else
	cpu_pinv = ALLOC_UINT(cpu_pinv, nprime);	cpu_pinv = ALIGN_UINT(cpu_pinv);
	cudaMalloc(&gpu_pinv    , (nprime)<<2);
#endif
	cudaMalloc(&gpu_popc_vec, ( 1024 )<<2);

  	cudaMemcpy(gpu_pdiff   , pdiff   ,  nprime         , cudaMemcpyHostToDevice);

	// Array of every-100th-prime, starting with p_last_small, presumed to be in 0-slot of said array.
	// We use ((m+1) % 100) == 0 as the write-datum condition, because because we really need the [99,199,299,...]th primes so
	// that the curr_p update which opens the loop in updateStartval like the following one yields the ( % 100) == 0 primes:
	cpu_p100[0] = curr_p = p_last_small;
	for(m = nclear; m < nprime; m++)
	{
		curr_p += (pdiff[m] << 1);
	//	if(m < 100)printf("curr_p[%u] = %u\n",m,curr_p);
	//	if(m < 100)printf("curr_p[%u] = %u; (curr_p-19)/2 = %u\n",m,curr_p,(curr_p-19)/2);
		cpu_prim[m] = curr_p;
	#ifdef FLOAT_PINV	// Float-inverse version:
		cpu_pinv[m] = 1.0f/curr_p;
	#else
		cpu_pinv[m] = gen_pinv(curr_p);
	//	if(m < 100)printf("[%u]:%u,%u : p*pinv = %u\n",m,curr_p,cpu_pinv[m], curr_p*cpu_pinv[m]);
	#endif
		if(((m+1) % 100) == 0) {
			cpu_p100[(m+1)/100] = curr_p;
	//	if(m < 1000)printf("Idx = %2u: p100 = %u\n",m/100,curr_p);
		}
	}
	cudaMemcpy(gpu_prim    , cpu_prim, (nprime    )<<2 , cudaMemcpyHostToDevice);
	cudaMemcpy(gpu_p100    , cpu_p100, (nprime/100)<<2 , cudaMemcpyHostToDevice);
	cudaMemcpy(gpu_pinv    , cpu_pinv, (nprime     <<2), cudaMemcpyHostToDevice);
}	// End of 1-time-only inits

	// The CPU-generated master bitmap and startval data vary with (k mod 60) pass, thus redo these each time this function is called:
//intf("k0 = %llu: Copying %u-byte CPU_bitmap in [%u]x fashion to GPU\n", kstart, bit_len >> 3, NSIEVE_COPIES_IN_GPU_SIEVE);
	for(m = 0; m < NSIEVE_COPIES_IN_GPU_SIEVE; m++) {
		//         vvvvvvvv Cast to byte here to ease offset math (offset = m x [same byte count used to control size of xfer])
		cudaMemcpy((uint8*)gpu_bit_map + m*( bit_len >> 3 ), bit_map, ( bit_len >> 3 ), cudaMemcpyHostToDevice);
	}
	#if 0
		// Debug: compare lead/trail words of CPU bitmap to first 64-bit 'catgap' word of clone-copy:
		//	bit_map:     [bit_map_len-1:0x000000000000|WXYZ] ... [0:0xFEDCBA9876543210] (bit & word significance increasing leftward)
		//	gpu_map: ... [bit_map_len-1:0xBA9876543210|WXYZ] ... [0:0xFEDCBA9876543210]
		cudaMemcpy(&kfac, (uint64*)gpu_bit_map + bit_map_len-1, 1<<3, cudaMemcpyDeviceToHost);
		printf("Lead48/Trail16 bits of source = 0x[%12llX][%4X]\n",bit_map[0] & 0x0000ffffffffffff,bit_map[bit_map_len-1]);
		printf("First Catgap word in GPU copy = 0x[%12llX][%4X]\n",kfac >> 16                     ,kfac & 0xffff         );
		exit(0);
	#endif

	// Compute startbit k (occurrence of first multiple of prime curr_p in first pass through the relevant sievelet:
	if(p <= MAX_SIEVING_PRIME)							//   v----- Restrict to 1-rod p (and 1-word 2*p) for now
		get_startval(MODULUS_TYPE, (uint64)p, findex, two_p, 1, bit_len, interval_lo, incr, nclear, nprime, p_last_small, pdiff, startval);
	else
		get_startval(MODULUS_TYPE,      0ull, findex, two_p, 1, bit_len, interval_lo, incr, nclear, nprime, p_last_small, pdiff, startval);

//printf("k0 = %llu, bit_len = %u, incr = %u, nclear = %u, nprime = %u, p_last_small = %u, pdiff[10] = %u\n",kstart, bit_len, incr, nclear, nprime, p_last_small, pdiff[10]);
//printf("sweep %llu, startval_sample = %u,%u,%u,%u\n",interval_lo,startval[1000],startval[2000],startval[3000],startval[4000]);

	cudaMemcpy(gpu_startval, startval, (nprime << 2), cudaMemcpyHostToDevice);

	for(sweep = interval_lo; sweep < interval_hi; sweep += NSIEVE_COPIES_IN_GPU_SIEVE*NMASTER_COPIES_IN_LIVE_SIEVE)
	{
		// Load fresh copy of master sievelet:
		uint32 wds_per_thr = 34;	// Map contains some integer multiple of 34034 = 2.7.11.13.17 32-bit words, so copy 2.17 = 34 words per thread per master copy
		nthr = gpu_map_len / wds_per_thr;
//printf("sweep = %u [incr = %u]\n",sweep, NSIEVE_COPIES_IN_GPU_SIEVE*NMASTER_COPIES_IN_LIVE_SIEVE);
//if(sweep == interval_lo)printf("#32-bit words in GPU-map = %u: Using %u threads, %u wds/thread in refreshBitmap\n",gpu_map_len,nthr,wds_per_thr*NMASTER_COPIES_IN_LIVE_SIEVE);
//if(sweep > interval_lo)
//exit(0);
		blocksPerGrid = ceil((float)nthr / threadsPerBlock);
		// GPU version: master is arg1, live is arg2. This places NMASTER_COPIES_IN_LIVE_SIEVE copies of GPU-master () in 'live' bitmap:
		refreshBitmap<<<blocksPerGrid, threadsPerBlock>>>(gpu_bit_map, gpu_bit_map2, nthr, wds_per_thr, NMASTER_COPIES_IN_LIVE_SIEVE);
	//
		if((sweep & 127) == 0)
		{
			printf(".");
			fflush(stdout);
		}
	//
		// To track lg(q) = lg(2.k.p+1), use approximation q ~= 2.k.p, thus lg(q) ~= lg(2.p) + lg(k).
		// At start of each pass through the k-based sieve, use 2.k.p with k = [starting k + sievebits]
		// to bound qmax from above, and compute lg(qmax) using the above logarithmic sum.
		fbits_in_k = log((double)k + 60*gpu_bit_len*NMASTER_COPIES_IN_LIVE_SIEVE)*ILG2;	// Use k-value at end of upcoming pass thru sieve as upper-bound
		fbits_in_q = fbits_in_2p + fbits_in_k;
//	if(fbits_in_q > 64)
//		printf("sweep = %d, k0 = %llu: fbits_in_q = fbits_in_2p [%10.4f] + fbits_in_k [%10.4f] = %10.4f\n",sweep,k,fbits_in_2p,fbits_in_k,fbits_in_q);

		int last_batch = 1;	// Fix = TRUE for now since testing each bitmap interval as it comes seems a tad faster
	//	int last_batch = (sweep == (interval_hi-1));	// TRUE triggers GPU_tfBatch to try the current set of factor candidates rather than waiting to accumulate more
		uint32 ntry = GPU_tfBatch(	// Return value = #candidates tried (set bits) in current sieve interval
			p, k, nclear, nprime,	// nprime = #sieving primes (counting from 3)
			gpu_prim,
			gpu_p100,
			gpu_pdiff,
			gpu_pinv,
			fbits_in_q,
			gpu_startval,
			gpu_bit_map2,
			gpu_kvec, kvsz,
			gpu_popc_vec,
			gpu_map_len*NMASTER_COPIES_IN_LIVE_SIEVE,	// #32-bit words in GPU live bitmap
			&kfac,
			last_batch
		);

		// Accumulate #candidates tried (set bits) in current sieve interval:
		count += ntry;
		// Return value of 0 from above function means it's still accumulating current batch of factor candidates, no factor-k in play yet:
		if((ntry > 0) && kfac) {
			factor_k[(*nfactor)++] = kfac;
			printf("\n\tFound factor with k = %llu.\n",kfac);
			kfac = 0;
		}
		// Update k to reflect another completed pass thru the bitmap:
		k += 60*gpu_bit_len*NMASTER_COPIES_IN_LIVE_SIEVE;

	} /* end of sweep loop	*/
	return count;
}

// TRYQ treated as 1 in this version. Return value = #candidates tried (set bits) in current sieve interval.
// Returns k-value of any factor found in current sieve interval in the input-pointer kfac arg.
__host__ uint32 GPU_tfBatch(
	const uint64 p,
	const uint64 k0,
	const uint32 nclear,
	const uint32 nprime,		// #sieving primes (counting from 3)
	const uint32*gpu_prim,
	const uint32*gpu_p100,
	const uint8 *gpu_pdiff,
	const pinv_t*gpu_pinv,
	const double fbits_in_q,	// Bit-ness of largest candidate-q in each batch
	uint32*gpu_startval,
	uint32*gpu_bit_map,
	uint32*gpu_kvec, const uint32 kvsz,
	uint32*gpu_popc_vec,
	const uint32 bit_map_len,	// #32-bit words in GPU bitmap
	uint64*kfac,
	int last_batch
)
{
	const uint32 bit_len = bit_map_len << 5;	// #bits in GPU bitmap
	ASSERT(HERE, kfac != 0, "Null kfac pointer!");
	static uint32 ntry = 0;	// Use this to accumulate batches of candidates from multiple sieve intervals,
							// and only test current batch when get close to the kvec-dimensioned limit.
	// Fix first of these and adjust 2nd based on needed threadcount for the procedure in question
	int threadsPerBlock = 256, blocksPerGrid = 0;	ASSERT(HERE, threadsPerBlock == 256, "threadsPerBlock must have fixed value = 256!");
	uint32 i, m, nthr, curr_p;
	static uint64 pshift = 0;
	static uint32 jshift, zshift, start_index, leadb;
	uint32 FERMAT = isPow2_64(p)<<1;	// *2 is b/c need to add 2 to the usual Mers-mod residue in the Fermat case

/********************************************************************************************/
/***** Clear the bits corresponding to the small primes from the current-pass sievelet: *****/
/********************************************************************************************/

	// Primes > p_last_small but < 64:
	m = nclear; curr_p = 17;
#if 0
	nthr = bit_map_len/2;	// Break bitmap into 64-bit chunk-per-thread in first phase of small-p clearing
	blocksPerGrid = ceil((float)nthr / threadsPerBlock);
//if(k0==5){ cudaThreadSynchronize();	printf("Input bitmap: ");	popcountBitmap<<<1,1>>>(gpu_bit_map, bit_map_len); }
	clrPrimes_Lt64      <<<blocksPerGrid, threadsPerBlock>>>((uint64 *)gpu_bit_map, nthr, gpu_pinv+m, gpu_startval+m);
//cudaThreadSynchronize();	printf("After clearing p < 64: ");	popcountBitmap<<<1,1>>>(gpu_bit_map, bit_map_len);
	m = nclear + 11;	curr_p = 61;
	// Primes > 64 but < 128:
	nthr /= 2;
	blocksPerGrid = ceil((float)nthr / threadsPerBlock);
	clrPrimes_Gt64_Lt128<<<blocksPerGrid, threadsPerBlock>>>((uint32 *)gpu_bit_map, nthr, gpu_pinv+m, gpu_startval+m);
#elif 0
	// Fusing above 2 kernel calls into 1 gives a nice speedup:
	nthr = bit_map_len/4;	// Break bitmap into 128-bit (4-word) chunks-per-thread in this phase of small-p clearing
	blocksPerGrid = ceil((float)nthr / threadsPerBlock);
	clrPrimes_Lt128     <<<blocksPerGrid, threadsPerBlock>>>((uint32 *)gpu_bit_map, nthr, gpu_pinv+m, gpu_startval+m);
	m = nclear + 24;	curr_p = 127;
#elif 0	// On my GTX430, this option is best (but only a smidge faster than the clrPrimes_Lt128 option)
	// Primes: > p_last_small but < 540:
	nthr = bit_map_len/8;	// Break bitmap into 256-bit (8-word) chunks-per-thread in this phase of small-p clearing
	blocksPerGrid = ceil((float)nthr / threadsPerBlock);
	clrPrimes_Gt128_Lt540<<<blocksPerGrid, threadsPerBlock>>>((uint32 *)gpu_bit_map, nthr, gpu_pinv+m, gpu_startval+m);
	m = 98;	curr_p = 523;
//cudaThreadSynchronize();	printf("After clearing p < 540: ");	popcountBitmap<<<1,1>>>(gpu_bit_map, bit_map_len);
#endif

//	uint32 num_p = nprime;	*** Until I de-crapify my GPU sieve no point in using more than 1000 primes ***
	uint32 num_p = 1000;	// Must be <= nprime, and corr. to #primes smaller than those done in Lucha-Libre bit-clearing phase
	uint32 np = num_p - m;	// Subtract #primes built into the sieve and alredy handled above,
							// so np reflects true #primes cleared after this next call.

	if(pshift == 0ull && num_p < nprime)	// (pshift == 0) in this function is equivalent to first_entry = TRUE
	printf("INFO: GPU sieve overriding factor.c setting [%u] and using only the first %u odd primes\n",nprime,num_p);

	uint32 lg_nwd = 3, nwords = 1<<lg_nwd;	// #64-bit words per thread ... 8 appears optimal in my early tests.
	nthr = (bit_map_len + (nwords-1))>>lg_nwd;
//printf("#p < 1000: nthr = %u\n",nthr);
	blocksPerGrid = ceil((float)nthr / threadsPerBlock);
	clrPrimes_GtA_LtB<<<blocksPerGrid, threadsPerBlock>>>(gpu_bit_map, nthr, nwords, np, gpu_pdiff+m, gpu_pinv+m, gpu_startval+m, curr_p);
//	cudaThreadSynchronize();	printf("After clearing first %u primes: ",num_p);	popcountBitmap<<<1,1>>>(gpu_bit_map, bit_map_len);

#if 0
	1000: TD = 253313872	NOTE: The following clrPrimes_LuchaLibre-based popcounts are not precisely reproducible from run-to-run due to use of non-atomic ANDing:
	10^4: TD = 196744856	<*** LL code below finds factor, but gives 209459650 ... likely write-collisions between threads
	10^5: TD = 162878654	<*** LL code below finds factor, but gives 175764726 ... runtime ~same as above 2, but that is *way* better than using above code for all primes
																	166785672 using byte-sized masking, which cuts #write-collisions by ~4x as expected
	Input bitmap:
	bitmapPop = 184324
	After clearing p < 64:
	bitmapPop = 134342
	After clearing first 1000 primes:
	bitmapPop = 63619
	After clearing first 10000 primes:
	bitmapPop = 49967	(LL: 50940; shared-mem: 60503)
	After clearing first 100000 primes:
	bitmapPop = 43697	<*** This LL-done one will vary slightly from run to run *** 41003 using byte-sized masking, which cuts #write-collisions as expected

	Debug: Doing just first 256 primes > #1000:
	bitmapPop = 61905	(run 2: 61915)
	shared gives: 63277	(run 2: 63276)
#endif

#if 0
	// 2/28: Larger-primes "Lucha Libre" step, in which each remaining prime-to-be-cleared gets its own thread:
  #if !USE_SHARED_MEM	// 'else' has shared-mem version (limit 256 threads) which is FUBAR, NEEDS DEBUG.
	nthr = nprime - num_p;
	blocksPerGrid = ceil((float)nthr / threadsPerBlock);
	clrPrimes_LuchaLibre<<<blocksPerGrid, threadsPerBlock>>>(gpu_bit_map, nthr, bit_len, gpu_prim+num_p, gpu_startval+num_p);
  #else
	#error Shared-mem usage in bit-clearing not currently supported!
	m = (nprime - num_p)/256;
	blocksPerGrid = 1;	// Limit total thread count to threadsPerBlock due to shared-mem usage in the function
	clrPrimes_LuchaLibre<<<blocksPerGrid, threadsPerBlock>>>(gpu_bit_map, m   , bit_len, gpu_prim+num_p, gpu_startval+num_p);
  #endif
//if(k0==5){ cudaThreadSynchronize();	printf("After clearing first %u primes: ",nprime);	popcountBitmap<<<1,1>>>(gpu_bit_map, bit_map_len); }
//exit(0);	// May need to move this down ~20 lines to capture above printf

#endif

	// For the exact-clear-stage primes, update startval array to reflect a full pass through the sievelet:
	np = 100;	// #small primes per thread
	nthr = num_p/np;
	blocksPerGrid = ceil((float)nthr / threadsPerBlock);
	updateStartval<<<blocksPerGrid, threadsPerBlock>>>(bit_len, nthr, nclear, gpu_p100, gpu_pdiff, gpu_startval);

/******************************************************************/
/***** Count the surviving bits in the small-p-cleared sieve: *****/
/******************************************************************/

	blocksPerGrid = 1;
	// #words in live bitmap must be multiple of threadsPerBlock (256)
	const uint32 pad_map_len = ((bit_map_len + 255) & 0xffffff00);
//if(k0==5)printf("pad_map_len = %u\n",pad_map_len);
	ASSERT(HERE, pad_map_len%threadsPerBlock == 0, "pad_map_len not divisible by threadsPerBlock!");
	nwords = pad_map_len/threadsPerBlock;
//printf("nwds = %u\n",nwords);

#if 0	// All-in-one version of the 3-function sequence below:

	parseBitmap_InitKoff3in1<<<blocksPerGrid, threadsPerBlock>>>(gpu_bit_map, gpu_kvec,               threadsPerBlock, nwords);
	// The 0-element of gpu_bit_map stores ntry:
	cudaMemcpy(&ntry, gpu_bit_map, 1<<2, cudaMemcpyDeviceToHost);

#else	// Separate-functions has advantage that we can time/[;ay-with subcomponents:

  #ifdef SERIAL
  	nthr = 1001;
  	nwords = bit_map_len/nthr;
	ASSERT(HERE, bit_map_len%nthr == 0, "bit_map_len not divisible by nthr!");
//printf("nwords = %u\n",nwords);
	blocksPerGrid = ceil((float)nthr / threadsPerBlock);
  #else
  	nthr = threadsPerBlock;
	blocksPerGrid = 1;
	nwords = pad_map_len/threadsPerBlock;
  #endif
	parseBitmap_CountSetBits<<<blocksPerGrid, threadsPerBlock>>>(gpu_bit_map,           gpu_popc_vec, nthr, nwords);
	parseBitmap_PopcountSums<<<blocksPerGrid, threadsPerBlock>>>(             gpu_kvec, gpu_popc_vec, nthr);
	parseBitmap_InitKoffsets<<<blocksPerGrid, threadsPerBlock>>>(gpu_bit_map, gpu_kvec, gpu_popc_vec, nthr, nwords);
	// In this version, the highest [255] element of gpu_popc_vec stores ntry:
	cudaMemcpy(&ntry, gpu_popc_vec + nthr-1, 1<<2, cudaMemcpyDeviceToHost);
#endif
//if(k0==5)printf("parseBitmap_PopcountSums gives bitmapPop = %u\n",ntry);
//exit(0);

	// Still below our accumulation threshold?
	if(!last_batch && (ntry < 1000000))
		return 0;

	blocksPerGrid = ceil((float)ntry / threadsPerBlock);
//(k0 < 60)printf("In GPU_tfBatch with ntry = %u, [kbeg,kend] = [%llu,%llu] ==> <blocksPerGrid = %d, threadsPerBlock = %d>...\n",ntry, k0,k0+60*bit_len, blocksPerGrid, threadsPerBlock);
	if(ntry > kvsz) {
		printf("In GPU_tfBatch with ntry = %u, kvsz = %u\n",ntry, kvsz);
		ASSERT(HERE, ntry <= kvsz, "ntry exceeds allocated memory!");
	}

	if(fbits_in_q >= 96) {		// q > 2^96
		ASSERT(HERE, 0, "q's > 2^96 not yet handled by GPU code!");
  #ifdef USE_FLOAT
	} else if((fbits_in_q >= 64) && (fbits_in_q < 78) && (pshift != p + 78)) {			// q in (2^64, 2^78)
		// Compute auxiliary TF data needed for 78-bit moduli:
	//	printf("Setting modpow params for 78-bits.\n");
		pshift = p + 78;
		jshift = leadz64(pshift);
		// Extract leftmost 7 bits of pshift (if > 77, use the leftmost 6) and subtract from 78:
		leadb = ((pshift<<jshift) >> 57);
		if(leadb > 77) {
			leadb >>= 1;
			start_index =  64-jshift-6;	// Use only the leftmost 6 bits
		} else {
			start_index =  64-jshift-7;
		}
		zshift = 77 - leadb;
		zshift <<= 1;				// Doubling the shift count here takes cares of the first SQR_LOHI
		pshift = ~pshift;
  #endif
	} else if((fbits_in_q >= 64) && (pshift != p + 96)) {	// q in (2^64, 2^96) or (2^78, 2^96) if USE_FLOAT defined
		// Compute auxiliary TF data needed for 96-bit moduli:
	//	printf("Setting modpow params for 96-bits.\n");
		pshift = p + 96;
		jshift = leadz64(pshift);
		// Extract leftmost 7 bits of pshift (if > 95, use the leftmost 6) and subtract from 96: */
		leadb = ((pshift<<jshift) >> 57);	// In [64,127]
		if(leadb > 95) {
			leadb >>= 1;	// Guarantees that leadb in [48,95]
			start_index =  64-jshift-6;	// Use only the leftmost 6 bits
		} else {
			start_index =  64-jshift-7;
		}
		zshift = 95 - leadb;
		zshift <<= 1;				// Doubling the shift count here takes cares of the first SQR_LOHI
		pshift = ~pshift;
	} else if (pshift != p + 64) {							// q < 2^64
		// Compute auxiliary TF data needed for 64-bit moduli:
	//	printf("Setting modpow params for 64-bits.\n");
		pshift = p + 64;
		jshift = leadz64(pshift);
		// Extract leftmost 6 bits of pshift and subtract from 64:
		leadb = ((pshift<<jshift) >> 58);
		start_index = 64-jshift-6;
		zshift = 63 - leadb;
		zshift <<= 1;				// Doubling the shift count here takes cares of the first SQR_LOHI
		pshift = ~pshift;
	}

	blocksPerGrid = ceil((float)ntry / threadsPerBlock);
//	printf("fbits_in_q = %10.4f, ntry = %u: Sending range starting at k0 = %llu to GPU for testing\n",fbits_in_q,ntry,k0);

	// Call appropriate batch-TF function, based on bitness of max. candidate in each batch:
	if(fbits_in_q >= 96) {			// q > 2^96
		ASSERT(HERE, 0, "q's > 2^96 not yet handled by GPU code!");
  #ifndef USE_FLOAT
	} else if(fbits_in_q >= 64) {				// q in (2^64, 2^96)
	//	printf("96-bit batch-modpow...\n");
		GPU_TF96_pop64<<<blocksPerGrid, threadsPerBlock>>>(p, pshift, FERMAT, zshift, start_index, k0, (uint32*)gpu_kvec, ntry);
  #else	// USE_FLOAT = true:
	} else if(fbits_in_q >= 78) {	// q in (2^78, 2^96)
	//	printf("96-bit batch-modpow...\n");
		GPU_TF96_pop64<<<blocksPerGrid, threadsPerBlock>>>(p, pshift, FERMAT, zshift, start_index, k0, (uint32*)gpu_kvec, ntry);
	} else if(fbits_in_q >= 64) {				// q in (2^64, 2^78) - Use float-based modexp:
	//	printf("78-bit batch-modpow...\n");
		GPU_TF78_pop64<<<blocksPerGrid, threadsPerBlock>>>(p, pshift, FERMAT, zshift, start_index, k0, (uint32*)gpu_kvec, ntry);
  #endif
	} else {					// q < 2^64
	//	printf("64-bit batch-modpow...\n");
		GPU_TF64_pop64<<<blocksPerGrid, threadsPerBlock>>>(p, pshift, FERMAT, zshift, start_index, k0, (uint32*)gpu_kvec, ntry);
	}

	// Retrieve 'one-beyond' elt of kvec which stores any factor-k:
	cudaMemcpy(kfac, gpu_kvec+ntry, 1<<3, cudaMemcpyDeviceToHost);
	// Rezero the statics in preparation for start of next-batch accumulation (but first must make a copy of ntry, which is our return value):
	i = ntry;
	ntry = 0;
	return i;
}

#endif	// USE_GPU ?

