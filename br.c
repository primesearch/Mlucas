/*******************************************************************************
*                                                                              *
*   (C) 1997-2009 by Ernst W. Mayer.                                           *
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

/****************/

int reverse(uint32 i, uint32 n)
{
/*...take index i of a set of N = 2^k and return the bit-reversed complement integer.
     Since Mlucas isn't restricted to power-of-2 FFT lengths, we don't actually use
     this function, but include it for purposes of reference. */

	uint32 tmp = 0;

/*...Make sure each new N is a power of 2. For high-performance implementations
     (i.e. where one does tons of these,) one could skip this part after verifying
     it on an initial test run. */

	if((n >> trailz32(n)) != 1)
	{
		sprintf(cbuf,"FATAL: non-power-of-2 length encountered in REVERSE.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
	}

	n >>= 1;

	while(n)
	{
		tmp <<= 1;
		tmp += (i & 1);
		i >>= 1;
		n >>= 1;
	}

	return(tmp);
}

/****************/

/*   Given an integer array of length N and a general set of radices
     r1,...,rJ, with r1*...*rJ = N, performs a bit-reversal reordering of the data in
     the array with respect to the first J-1 radices, r1,...,r{J-1}. Note that to obtain
     a reordering suitable for a decimation-in-time FFT which consists of the above
     sequence of radices (or which would result from a decimation-in-frequency transform
     via the sequence of radices,) one must perform a bit-reversal which processes the radices
     in REVERSE (the direction, not the above function :) order, i.e. set the radix array to
     RADIX(1:J)=RADIX(J:1:-1) prior to calling the BIT_REVERSE_INT subroutine.

     Result is returned in the calling array.

     Some small examples that illustrate bit reversals with respect to a general sequence of radices
     and the differences between the power-of-2 and non-power-of-2 cases follow.

  ...Example 1: n=16, FFT having radices {2,2,2,2}: we do three BR passes, corresponding to radices {2,2,2}.
  	Initial vector:
  	X   =    1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16
  		-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --	<---16 length-1 blocks
  	Pass 1: stride = 16/2 = 8, blocklen=1 gives
  -->	X'  =    1, 9, 2,10, 3,11, 4,12, 5,13, 6,14, 7,15, 8,16
  		----- ----- ----- ----- ----- ----- ----- ----- <--- 8 length-2 blocks
  	Pass 2: stride = 16/2 = 8, blocklen=2
  -->	X'' =    1, 9, 5,13, 2,10, 6,14, 3,11, 7,15, 4,12, 8,16
  		----------- ----------- ----------- ----------- <--- 4 length-4 blocks
  	Pass 3: stride = 16/2 = 8, blocklen=4
  -->	X'' =    1, 9, 5,13, 3,11, 7,15, 2,10, 6,14, 4,12, 8,16
                 ----------------------- ----------------------- <--- 2 length-8 blocks
  	Done:   stride = 16/2 = 8, blocklen=8, i.e. stride=blocklen.

  	NOTE that for n a power of 2, the order of the radices is moot. Also, we can
  	achieve the same reordering as illustrated in pass-by-pass fashion above by
  	redefining X as a zero-offset vector of 4-bit elements,

  	X =    0,   1,   2,   3,   4,   5,   6,   7,   9,  10,  11,  12,  13,  14,  15,

  	or, in binary,

  	X = 0000,0001,0010,0011,0100,0101,0110,0111,1000,1001,1010,1011,1100,1101,1110,1111.

  	Reversing the order of the bits in each element gives

  	X = 0000,1000,0100,1100,0010,1010,0110,1110,0001,1001,0101,1101,0011,1011,0111,1111,

  	and rewriting these in decimal and adding one to each gives the same result as above:

  	X =    1,   9,   5,  13,   3,  11,   7,  15,   2,  10,   6,  14,   4,  12,  16,

  	hence the name BIT-REVERSAL REORDERING. This is a misnomer when the number of elements
  	is not a power of 2, but the term has become so ingrained in the FFT world that we'll
  	use it even in this case, with the understanding that we really mean a reordering
  	corresponding to a specific sequence of radices and having nothing to do with the
  	actual bit patterns of the element indices in question.

  ...Example 1: n=24, FFT having radices {2,2,2,3}: we do three BR passes, corresponding to radices {3,2,2} (Note the
     reversal of the order of the FFT radices in performing the BR - for power-of-2 it's all the same, but not in general  )
  	Initial vector:
  	X   =    0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23
  		-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --	-- -- -- -- -- -- -- -- <---24 length-1 blocks
  	Pass 1: stride = 24/3 =  8, blocklen=1
  -->	X'  =    0, 8,16, 1, 9,17, 2,10,18, 3,11,19, 4,12,20, 5,13,21, 6,14,22, 7,15,23
  		-------- -------- -------- -------- -------- -------- -------- -------- <--- 8 length-3 blocks
  	Pass 2: stride = 24/2 = 12, blocklen=3
  -->	X'' =    0, 8,16, 4,12,20, 1, 9,17, 5,13,21, 2,10,18, 6,14,22, 3,11,19, 7,15,23
  		----------------- ----------------- ----------------- ----------------- <--- 4 length-6 blocks
  	Pass 3: stride = 24/2 = 12, blocklen=6
  -->	X'''=    0, 8,16, 4,12,20, 2,10,18, 6,14,22, 1, 9,17, 5,13,21, 3,11,19, 7,15,23
  		----------------------------------- ----------------------------------- <--- 2 length-12 blocks
  	Done:   stride = 24/2 = 12, blocklen=12, i.e. stride=blocklen.

  ...Example 2: n=24, FFT having radices {3,2,2,2}: we do three BR passes, corresponding to radices {2,2,2}.
  	Initial vector:
  	X   =    0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23
  		-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --	-- -- -- -- -- -- -- -- <---24 length-1 blocks
  	Pass 1: stride = 24/2 = 12, blocklen=1
  -->	X'  =    0,12, 1,13, 2,14, 3,15, 4,16, 5,17, 6,18, 7,19, 8,20, 9,21,10,22,11,23
  		----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- <---12 length-2 blocks
  	Pass 2: stride = 24/2 = 12, blocklen=2
  -->	X'' =    0,12, 6,18, 1,13, 7,19, 2,14, 8,20, 3,15, 9,21, 4,16,10,22, 5,17,11,23
  		----------- ----------- ------------ ---------- ----------- ----------- <--- 6 length-4 blocks
  	Pass 3: stride = 24/2 = 12, blocklen=4
  -->	X'''=    0,12, 6,18, 3,15, 9,21, 1,13, 7,19, 4,16,10,22, 2,14, 8,20, 5,17,11,23
  		----------------------- ----------------------- ----------------------- <--- 3 length-8 blocks
  	Done:   stride = 24/3 =  8, blocklen=8, i.e. stride=blocklen.

  	NOTE how in the non-power-of-2 case, the specific order in which the radices are processed matters very much.

  ...Example 3: n=24, FFT having radices {2,3,2,2}: we do three BR passes, corresponding to radices {2,2,3}.
  	Initial vector:
  	X   =    0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23
  		-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --	-- -- -- -- -- -- -- -- <---24 length-1 blocks
  	Pass 1: stride = 24/2 = 12, blocklen=1
  -->	X'  =    0,12, 1,13, 2,14, 3,15, 4,16, 5,17, 6,18, 7,19, 8,20, 9,21,10,22,11,23
  		----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- <---12 length-2 blocks
  	Pass 2: stride = 24/2 = 12, blocklen=2
  -->	X'' =    0,12, 6,18, 1,13, 7,19, 2,14, 8,20, 3,15, 9,21, 4,16,10,22, 5,17,11,23
  		----------- ----------- ------------ ---------- ----------- ----------- <--- 6 length-4 blocks
  	Pass 3: stride = 24/3 =  8, blocklen=4
  -->	X'''=    0,12, 6,18, 2,14, 8,20, 4,16,10,22, 1,13, 7,19, 3,15, 9,21, 5,17,11,23
  		----------------------------------- ----------------------------------- <--- 2 length-12 blocks
  	Done:   stride = 24/2 = 12, blocklen=12, i.e. stride=blocklen.

     Next, look at how the (j,N-j) real/complex wrapper index correlations fare under the various reorderings.
     Since the 0 and N/2 = 12 elements get special treatment, we prefer them to wind up next to each other, i.e.
     the last pass of the forward transform must be a radix-2. This rules out example 1 above. For examples 2 and 3,
     the wrapper_square routine will treat the data as a set of blocks of ever-increasing size (the blocklength
     for the current pass is (radix_now - 1)*Blocklen_sum), and working through each block from both ends
     gives two sets of (j,N-j) correlations:

     Example 2: FFT having radices {3,2,2,2}:
     Blocklen	Blocklen_sum	radix_now-1	elements
     2		2		1		 0,12
     2		4		1		 6,18
     4		8		1		 3,15, 9,21
     16				2		 1,13, 7,19, 4,16,10,22, 2,14, 8,20, 5,17,11,23

     Example 3: FFT having radices {2,3,2,2}:
     Blocklen	Blocklen_sum	radix_now	elements
     2		2		1		 0,12
     2		4		1		 6,18
     8		12		2		 2,14, 8,20, 4,16,10,22
     12				1		 1,13, 7,19, 3,15, 9,21, 5,17,11,23
*/
void bit_reverse_int(int vec[], int n, int nradices, int radix[], int incr, int*arr_scratch)
{
	int blocklen,count,i,iput,j,jstart,k,l,nblock,stride;
	static int *tmp = 0x0;

	/* If no scratch-space array provided, create one locally: */
	if(arr_scratch)
	{
		/* Don't allow reuse of main array for inits at this time: */
		ASSERT(HERE, &vec[0] != &arr_scratch[0], "Array re-use not currently supported!");

		tmp = arr_scratch;
	}
	else
	{
		tmp = (int *)malloc(n*sizeof(int));
	}
	memset(tmp, 0, n*sizeof(int));

/*...Check that product of radices = vector length */

	i = 0;
	k = radix[i];
	for(count = 1; count < nradices; count++)
	{
	  i += incr;
	  k *= radix[i];
	}

	if(k != n)
	{
	  printf("FATAL: product of radices != vector length in BIT_REVERSE_INT.\n");
	  printf("Radices = ");
	  for(count = 0, i = 0; count < nradices; count++)
	  {
	    printf("%4d",radix[i]);
	    i += incr;
	  }
	  printf("\n");
	  printf("Length = %10d\n",n);
	  fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
	}

/*...We don't use the final radix for the bit reversal, we simply need it for array bounds checking. */

	nblock	=n;
	stride  =1;
	blocklen=1;

	for(count = 0, i = 0; count < nradices-1; count++)
	{
	  iput	=0;
	  stride=     n/radix[i];	/* STRIDE = number of data elements between adjacent fetch locations.	*/
	  nblock=nblock/radix[i];	/* NBLOCK = number of data blocks   between adjacent fetch locations.	*/

	  for(j = 0; j < stride; j += blocklen)		/* Index of first element fetched on each pass runs through a full stride length.	*/
	  {
	    jstart=j;

	    for(k = 0; k < radix[i]; k++)		/* On the current sweep, we fetch elements from RADIX[I] equally-spaced data blocks...	*/
	    {

	      for(l = jstart; l < jstart+blocklen; l++)	/* ...each of length BLOCKLEN.	*/
	      {
		tmp[iput++]=vec[l];	/* This is if we want VEC to contain the fetch locations of BR data, i.e. do BR via A=A(VEC).	*/
/*		tmp[l]=vec[iput++];*/	/* This is if we want VEC to contain the  put  locations of BR data, i.e. do BR via A(VEC)=A.	*/
	      }

	      jstart += stride;	/* Start of next fetch block.	*/

	    }

	  }

	  for(j = 0; j < n; j++)
	  {
	    vec[j] = tmp[j];	/* Put result of current pass back into VEC.	*/
	  }

	  blocklen *= radix[i];
	  i += incr;
	}

	if(!arr_scratch)
	{
		free(tmp); tmp = 0x0;
	}
}
