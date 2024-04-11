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

/***************/

/* Initialize the various arrays of indices used in radix16_wrapper_square, so we can execute
   the processing of the [radix0] disjoint data blocks by that routine in parallel, if desired.
*/
void radix16_wrapper_ini(int n, int radix0, int iblock, int nradices_prim, int radix_prim[], int ws_i[], int ws_j1[], int ws_j2[], int ws_j2_start[], int ws_k[], int ws_m[], int ws_blocklen[], int ws_blocklen_sum[])
{
	static int i,j1,j2,j2_start,k,m,blocklen,blocklen_sum;
	int iblock_next;

	if(iblock <= 1 && !(radix0 & 1))
	  	iblock_next = iblock + 1;
	else
		iblock_next = iblock + 2;

	if(iblock == 0)	// j1 = real-array index (double the complex-array index) of the 1st element of each floating pair.
	{
		// No need to init I and M here, since they are set by entry into the nested I/M loop in radix16_pairFFT_mul_square:
		j1           =  0;
		j2           = 32;
		j2_start     = j2;	// j2 = real-array index (double the complex-array index) of 2nd element of each floating pair.
		k            =  0;
		blocklen     = 16;	// = half of complex blocklength, since process 2 complex data for each value of loop index L.
		blocklen_sum =  0;

		ws_i           [iblock] = i           ;
		ws_j1          [iblock] = j1          ;
		ws_j2          [iblock] = j2          ;
		ws_j2_start    [iblock] = j2_start    ;
		ws_k           [iblock] = k           ;
		ws_m           [iblock] = m           ;
		ws_blocklen    [iblock] = blocklen    ;
		ws_blocklen_sum[iblock] = blocklen_sum;
	} else {
		goto jump_in;
	}

	for(i = nradices_prim-5; i >= 0; i-- )	// Main loop: lower bound = nradices_prim - radix_now.
	{										// Remember, radices get processed in reverse order here as in forward FFT.
		for(m = 0; m < (blocklen-1)>>1; m += 8) // Do two 16-element sets per loop, so only execute loop half as many times as before.
		{
			// This tells us when we've reached the end of the current data block:
			// Apr 2014: Must store intermediate product j1*radix0 in a 64-bit int to prevent overflow!
			if(j1 && ((uint64)j1*radix0)%n == 0)
			{
				ws_i           [iblock_next] = i           ;
				ws_j1          [iblock_next] = j1          ;
				ws_j2          [iblock_next] = j2          ;
				ws_j2_start    [iblock_next] = j2_start    ;
				ws_k           [iblock_next] = k           ;
				ws_m           [iblock_next] = m           ;
				ws_blocklen    [iblock_next] = blocklen    ;
				ws_blocklen_sum[iblock_next] = blocklen_sum;
			//	printf("%8" PRIu64 "  %20" PRIu64 "  %8" PRIu64 ": init ws_k[%3d] = %10d\n",j1,((uint64)j1*radix0),j2,iblock_next,k);
				return;
			}
	jump_in:	// Entry point for all blocks but the first.
			k += 2;	// increment sincos array index
			// And update the data (j1 and j2) array indices:
			j1 += 32;
			j2 -= 32;
		}
	/*
	!...Since the foregoing loop only gets executed half as many times as in the simple version, to properly position
	!   ourselves in the data array for the start of the next block, need to bump up j1 by as much as would occur in a
	!   second execution of the above loop. The exception is the first loop execution, where j1 needs to be doubled (32 x 2).
	*/
		j1 += (blocklen << 1);
		if(j2_start == n-32) {
		//	printf("(j2_start == n-32) return with j2_start = %d\n",j2_start);
			return;
		}

	/*...Reset half-complex-blocklength for next pass. If K >> 1 has a zero trailing bit,
		 we multiply the blocklength by K >> 1 in preparation for the final block.	*/

		blocklen_sum += blocklen;
		blocklen = (blocklen_sum) * (radix_prim[i-1]-1);

	/*...Next j2_start is previous one plus the (real) length of the current block = 4*(half-complex-blocklength) */

		j2_start += (blocklen<<2);
		j2 = j2_start;			/* Reset j2 for start of the next block. */
	//	printf("newblock: blocklen = %8d blocklen_sum = %8d j2 = %8d\n",blocklen,blocklen_sum,j2);
	}	 /* End of Main loop */
}

/*
Jun 2014: Possible UMR bug? Note the aside from i=1, only even-idx elts of the ws-arrays get inited ...
so how do the odd-index reads not hose the result?

Using complex FFT radices        16         8        16        16        16
                                          init ws_k[  0] =          0
   65536               1048576    131040: init ws_k[  1] =       2048
  131072               2097152    262112: init ws_k[  2] =       4096
  262144               4194304    524256: init ws_k[  4] =       8192
  327680               5242880    458720: init ws_k[  6] =      12288
  524288               8388608   1048544: init ws_k[  8] =      16384
  589824               9437184    983008: init ws_k[ 10] =      20480
  655360              10485760    917472: init ws_k[ 12] =      24576
  720896              11534336    851936: init ws_k[ 14] =      28672
Mers_mod_square: Init threadpool of 1 threads
Setting CPU = 0 affinity of worker thread id 0, mach_id = 3843
radix16_wrapper_square with ws[]-index = 0
stride = 32
On entry: i = 0, j1,j2,j2_start = 0, 32, 32, k,m = 0, 0, nrad_prim = 19, blocklen,sum = 16
radix16_wrapper_square with ws[]-index = 1
stride = 32
On entry: i = 3, j1,j2,j2_start = 65536, 131040, 131040, k,m = 2048, 0, nrad_prim = 19, blocklen,sum = 16384
radix16_wrapper_square with ws[]-index = 2
stride = 32
On entry: i = 2, j1,j2,j2_start = 131072, 262112, 262112, k,m = 4096, 0, nrad_prim = 19, blocklen,sum = 32768
radix16_wrapper_square with ws[]-index = 3	<*** inited where? ***
stride = 32						vvvv <*** j1 = 0, so no "jump_in": ***
On entry: i = 0, j1,j2,j2_start = 0, 0, 0, k,m = 0, 0, nrad_prim = 19, blocklen,sum = 0
========================
So whn j1 = 0 on entry we exit immediately via:
		if(j1 && ((uint64)j1*radix0)%n == 0)
		{
		//	fprintf(stderr,"(j1 && j1*radix0 == 0 (mod n)) check hit: returning\n");
			return;
		}
========================
radix16_wrapper_square with ws[]-index = 4
stride = 32
On entry: i = 1, j1,j2,j2_start = 262144, 524256, 524256, k,m = 8192, 0, nrad_prim = 19, blocklen,sum = 65536
radix16_wrapper_square with ws[]-index = 5	<*** j1 = 0
stride = 32
On entry: i = 0, j1,j2,j2_start = 0, 0, 0, k,m = 0, 0, nrad_prim = 19, blocklen,sum = 0
radix16_wrapper_square with ws[]-index = 6
stride = 32
On entry: i = 1, j1,j2,j2_start = 327680, 458720, 524256, k,m = 12288, 16384, nrad_prim = 19, blocklen,sum = 65536
radix16_wrapper_square with ws[]-index = 7	<*** j1 = 0
stride = 32
On entry: i = 0, j1,j2,j2_start = 0, 0, 0, k,m = 0, 0, nrad_prim = 19, blocklen,sum = 0
radix16_wrapper_square with ws[]-index = 8
stride = 32
On entry: i = 0, j1,j2,j2_start = 524288, 1048544, 1048544, k,m = 16384, 0, nrad_prim = 19, blocklen,sum = 131072
radix16_wrapper_square with ws[]-index = 9	<*** j1 = 0
stride = 32
On entry: i = 0, j1,j2,j2_start = 0, 0, 0, k,m = 0, 0, nrad_prim = 19, blocklen,sum = 0
radix16_wrapper_square with ws[]-index = 10
stride = 32
On entry: i = 0, j1,j2,j2_start = 589824, 983008, 1048544, k,m = 20480, 16384, nrad_prim = 19, blocklen,sum = 131072
radix16_wrapper_square with ws[]-index = 11	<*** j1 = 0
stride = 32
On entry: i = 0, j1,j2,j2_start = 0, 0, 0, k,m = 0, 0, nrad_prim = 19, blocklen,sum = 0
radix16_wrapper_square with ws[]-index = 12
stride = 32
On entry: i = 0, j1,j2,j2_start = 655360, 917472, 1048544, k,m = 24576, 32768, nrad_prim = 19, blocklen,sum = 131072
radix16_wrapper_square with ws[]-index = 13	<*** j1 = 0
stride = 32
On entry: i = 0, j1,j2,j2_start = 0, 0, 0, k,m = 0, 0, nrad_prim = 19, blocklen,sum = 0
radix16_wrapper_square with ws[]-index = 14
stride = 32
On entry: i = 0, j1,j2,j2_start = 720896, 851936, 1048544, k,m = 28672, 49152, nrad_prim = 19, blocklen,sum = 131072
radix16_wrapper_square with ws[]-index = 15	<*** j1 = 0
stride = 32
On entry: i = 0, j1,j2,j2_start = 0, 0, 0, k,m = 0, 0, nrad_prim = 19, blocklen,sum = 0

Thus, j1 = 0 is how the odd-idx uninit is handled in practice - BUT NEED TO ENSURE ALL THE J1-DATA ARE INITED = 0 AT OUTSET

Thus, switch ws_* allocs in mers_mod_square from malloc to calloc.

(Surprised this issue took so long to manifest...)
*/
