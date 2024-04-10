/*******************************************************************************
*                                                                              *
*   (C) 1997-2016 by Ernst W. Mayer.                                           *
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

/* To test a particular DFT routine:

First, for the DFT in question, assuming it is a new under-development one,
MAKE SURE THE OUTPUT INDEXING IS STRICLY ASCENDING-ORDERED! - The test_fft_radix() index-display assumes this.

1. To set desired radix for testing and test-type (see comments around TTYPE below), compile this file
using (note the [] are just for emphasis):
	gcc/clang -c -DRADIX=[desired radix] -DTTYPE=[desired test type] test_fft_radix.c

2. Compile util.c via 'gcc/clang -c -DTEST_FFT_RADIX util.c' [Note: If -DUSE_THREADS was|was-not used for other
	sourcefiles in the build, need to use|not-use it here, too];

3. Link into an Mlucas binary and run-as-if-timing-test at any desired length to test the DFT radix.

To deduce the required output-idx ordering, sort Actual-outputs data by left col-of-real-parts, move 
resulting [re,im,i] block to file, repeat procedure for Expected-outputs, then compare the two files.
If they differ only in their respective rcol indices, those two cols give the required index mapping,
i.e. the desired output-index mapping.
*/

#include "Mlucas.h"

#ifndef RADIX
	#define RADIX	1024
#endif

// May 2014: In preparation for the push toward DFT radices ~4096 (which allow mod-F33 convolution
// to be done with 1 MB chunksizes, a crucial performance-related parameter on x86), up this from 1024
// to 4096 and quadruple the digits-of-Pi used as entries in our test-inputs const array:
#if RADIX > 4096
	#error ref-digits const array needs expanding, and radix_prim array needs dimension > 12!
#endif

#ifndef TTYPE
	#define TTYPE	3	// 0 = Show input-index-scramblings-needed for DIF and DIT, 1 = DIF, 2 = DIT,
						// 3 = DIF+DIT (which should return the original inputs after dividing by N)
#endif

// Structs used in index-permutation extraction:
struct int_pair{
	int a;
	int b;
};
struct idx_cmplx{
	int idx;
	struct complex cdat;
};

// Sorting predicate for the int_pair struct:
int int_pair_cmp(const void *va, const void *vb)	// 'v' prefix stands for 'void' in re, the pointer type required by C qsort
{
	/* Need to store arguments in appropriate type before using */
	const int a = ((const struct int_pair *)va)->a;	// <*** Use this predicate for sorting int-pair lists by *first* index,
	const int b = ((const struct int_pair *)vb)->a;	// that's why we compare the a-data of the 2 inputs
	/* Return -1 or +1 if a < b or a > b, resp: */
	if (a > b) return +1;
	else if (a < b) return -1;
	/* Return 0 if equal - Pull out of the conditional so this can also serve as a catchall return value,
	which avoids the usual compiler warning about missing return-val due to most compilers not being able
	to discern that if(>)/elif(<)/else covers all possibilities: */
	return 0;
}

// Need 2 flavors of Sorting predicates for the idx_cmplx struct:
// [1] Sort-by-index:
int cmp_by_idx(const void *va, const void *vb)
{
	const int a = ((const struct idx_cmplx *)va)->idx;
	const int b = ((const struct idx_cmplx *)vb)->idx;

	if (a > b) return +1;
	else if (a < b) return -1;
	return 0;
}

// [2] Sort-by-value: Treating the re/im parts of the complex datum
// like high and low words of a 2-word multiprecision datum suffices for our present purposes:
int cmp_by_val(const void *va, const void *vb)
{
	const double EPS = 1e-10;	double adiff;
	int retval = 0;
	const double ar = (((const struct idx_cmplx *)va)->cdat).re, ai = (((const struct idx_cmplx *)va)->cdat).im;
	const double br = (((const struct idx_cmplx *)vb)->cdat).re, bi = (((const struct idx_cmplx *)vb)->cdat).im;
	//printf("cmp_by_val: [ar,ai] = %20.10f,%20.10f; [br,bi] = %20.10f,%20.10f\n",ar,ai,br,bi);
	adiff = fabs(ar - br);
	if(adiff > EPS) {
		if (ar > br) {
			retval = +1;	//printf("cmp_by_val: ar > br: retval = %2d\n",retval);
		} else if (ar < br) {
			retval = -1;	//printf("cmp_by_val: ar < br: retval = %2d\n",retval);
		}
		//printf("cmp_by_val: [ar-br] = %20.10e: retval = %2d\n",adiff,retval);
		return retval;
	}
	/* If real parts equal, proceed to imaginary parts: */
	adiff = fabs(ai - bi);
	if(adiff > EPS) {
		if (ai > bi) {
			retval = +1;	//printf("cmp_by_val: ai > bi: retval = %2d\n",retval);
		} else if (ai < bi) {
			retval = -1;	//printf("cmp_by_val: ai < bi: retval = %2d\n",retval);
		}
		//printf("cmp_by_val: [ai-bi] = %20.10e: retval = %2d\n",adiff,retval);
		return retval;
	}
	//printf("cmp_by_val: complex data equal within roundoff tolerance, retval = %2d\n",retval);
	return retval;
}
/*
30 Nov 2015: Bizarre result was hosing a small-DFT test (radix-8 DIF) with 2 identical Re-part outputs, [-6,-14] and [-6,-10]:
This one OK:
	cmp_by_val: [ar,ai] =        -6.0000000000,      -14.0000000000; [br,bi] =        -6.0000000000,      -10.0000000000
	cmp_by_val: ar == ai: Im-parts retval = -1
This is BAD:
	cmp_by_val: [ar,ai] =        -6.0000000000,      -14.0000000000; [br,bi] =        -6.0000000000,      -10.0000000000
	cmp_by_val: ar > br: retval =  1
Had to add ROE-tolerance 'fuzzy compare' to cure that.
*/

void matmul_double (double **, double *, double *, int, int);
void matmul_complex(struct complex **, struct complex *, struct complex *, int, int);
#ifdef USE_FGT61
// Modular analog of matmul_complex, over FGT(M61^2):
void matmul_fgtmod (uint128 **, uint128 *, uint128 *, int, int);
#endif

#define CABS(a,b)	sqrt((a)*(a) + (b)*(b))

void test_fft_radix(void)
{
	struct int_pair ipair[RADIX];
	struct idx_cmplx sort_arr1[RADIX], sort_arr2[RADIX];
	const int sz_int_pair = sizeof(struct int_pair), sz_idx_cmplx = sizeof(struct idx_cmplx);
	const double err_theshold = 1e-10, sqrt_radix = sqrt(RADIX);
	double iradix = 1.0/RADIX;
	int rmul = 2*RE_IM_STRIDE;
	int i,j,j1,j2,k,l,nradices, *index = 0x0, nerr,idiff, print_pass, pow2, podd;
	// aadix_prim stores factorization of RADIX in order reflecting order of sub-radix processing;
	int radix_prim[12];
	int *dit_scramble = 0x0;	/* This holds the input-perm for the DIF of the given length ... compute this using the primitive radices */
	double err_r, err_i, abserr, maxerr, avgerr;
	/* "Random" inputs are just decimal digits of Pi ... make this as big as needed, currently support
	up to complex length = 4096. NB: Do NOT use Unix bc for this! (Way too slow). use pari, e.g.:

		? default(realprecision,8200)
		? Pi

	For the curious, the number of occurrences of each decimal digit [0-9] in these first 8192 is as follows:

		775,850,825,801,825,854,841,805,785,831

	The latest occurrence of any digit 0-9 is 0, which does not occur until the 33rd digit.
	*/
	const char ref[8192] = {	// RADIX/32 rows of 64 real (32 complex-pair) data each:
	3,1,4,1,5,9,2,6,5,3,5,8,9,7,9,3,2,3,8,4,6,2,6,4,3,3,8,3,2,7,9,5,0,2,8,8,4,1,9,7,1,6,9,3,9,9,3,7,5,1,0,5,8,2,0,9,7,4,9,4,4,5,9,2,
	3,0,7,8,1,6,4,0,6,2,8,6,2,0,8,9,9,8,6,2,8,0,3,4,8,2,5,3,4,2,1,1,7,0,6,7,9,8,2,1,4,8,0,8,6,5,1,3,2,8,2,3,0,6,6,4,7,0,9,3,8,4,4,6,
	0,9,5,5,0,5,8,2,2,3,1,7,2,5,3,5,9,4,0,8,1,2,8,4,8,1,1,1,7,4,5,0,2,8,4,1,0,2,7,0,1,9,3,8,5,2,1,1,0,5,5,5,9,6,4,4,6,2,2,9,4,8,9,5,
	4,9,3,0,3,8,1,9,6,4,4,2,8,8,1,0,9,7,5,6,6,5,9,3,3,4,4,6,1,2,8,4,7,5,6,4,8,2,3,3,7,8,6,7,8,3,1,6,5,2,7,1,2,0,1,9,0,9,1,4,5,6,4,8,
	5,6,6,9,2,3,4,6,0,3,4,8,6,1,0,4,5,4,3,2,6,6,4,8,2,1,3,3,9,3,6,0,7,2,6,0,2,4,9,1,4,1,2,7,3,7,2,4,5,8,7,0,0,6,6,0,6,3,1,5,5,8,8,1,
	7,4,8,8,1,5,2,0,9,2,0,9,6,2,8,2,9,2,5,4,0,9,1,7,1,5,3,6,4,3,6,7,8,9,2,5,9,0,3,6,0,0,1,1,3,3,0,5,3,0,5,4,8,8,2,0,4,6,6,5,2,1,3,8,
	4,1,4,6,9,5,1,9,4,1,5,1,1,6,0,9,4,3,3,0,5,7,2,7,0,3,6,5,7,5,9,5,9,1,9,5,3,0,9,2,1,8,6,1,1,7,3,8,1,9,3,2,6,1,1,7,9,3,1,0,5,1,1,8,
	5,4,8,0,7,4,4,6,2,3,7,9,9,6,2,7,4,9,5,6,7,3,5,1,8,8,5,7,5,2,7,2,4,8,9,1,2,2,7,9,3,8,1,8,3,0,1,1,9,4,9,1,2,9,8,3,3,6,7,3,3,6,2,4,
	4,0,6,5,6,6,4,3,0,8,6,0,2,1,3,9,4,9,4,6,3,9,5,2,2,4,7,3,7,1,9,0,7,0,2,1,7,9,8,6,0,9,4,3,7,0,2,7,7,0,5,3,9,2,1,7,1,7,6,2,9,3,1,7,
	6,7,5,2,3,8,4,6,7,4,8,1,8,4,6,7,6,6,9,4,0,5,1,3,2,0,0,0,5,6,8,1,2,7,1,4,5,2,6,3,5,6,0,8,2,7,7,8,5,7,7,1,3,4,2,7,5,7,7,8,9,6,0,9,
	1,7,3,6,3,7,1,7,8,7,2,1,4,6,8,4,4,0,9,0,1,2,2,4,9,5,3,4,3,0,1,4,6,5,4,9,5,8,5,3,7,1,0,5,0,7,9,2,2,7,9,6,8,9,2,5,8,9,2,3,5,4,2,0,
	1,9,9,5,6,1,1,2,1,2,9,0,2,1,9,6,0,8,6,4,0,3,4,4,1,8,1,5,9,8,1,3,6,2,9,7,7,4,7,7,1,3,0,9,9,6,0,5,1,8,7,0,7,2,1,1,3,4,9,9,9,9,9,9,
	8,3,7,2,9,7,8,0,4,9,9,5,1,0,5,9,7,3,1,7,3,2,8,1,6,0,9,6,3,1,8,5,9,5,0,2,4,4,5,9,4,5,5,3,4,6,9,0,8,3,0,2,6,4,2,5,2,2,3,0,8,2,5,3,
	3,4,4,6,8,5,0,3,5,2,6,1,9,3,1,1,8,8,1,7,1,0,1,0,0,0,3,1,3,7,8,3,8,7,5,2,8,8,6,5,8,7,5,3,3,2,0,8,3,8,1,4,2,0,6,1,7,1,7,7,6,6,9,1,
	4,7,3,0,3,5,9,8,2,5,3,4,9,0,4,2,8,7,5,5,4,6,8,7,3,1,1,5,9,5,6,2,8,6,3,8,8,2,3,5,3,7,8,7,5,9,3,7,5,1,9,5,7,7,8,1,8,5,7,7,8,0,5,3,
	2,1,7,1,2,2,6,8,0,6,6,1,3,0,0,1,9,2,7,8,7,6,6,1,1,1,9,5,9,0,9,2,1,6,4,2,0,1,9,8,9,3,8,0,9,5,2,5,7,2,0,1,0,6,5,4,8,5,8,6,3,2,7,8,
	8,6,5,9,3,6,1,5,3,3,8,1,8,2,7,9,6,8,2,3,0,3,0,1,9,5,2,0,3,5,3,0,1,8,5,2,9,6,8,9,9,5,7,7,3,6,2,2,5,9,9,4,1,3,8,9,1,2,4,9,7,2,1,7,
	7,5,2,8,3,4,7,9,1,3,1,5,1,5,5,7,4,8,5,7,2,4,2,4,5,4,1,5,0,6,9,5,9,5,0,8,2,9,5,3,3,1,1,6,8,6,1,7,2,7,8,5,5,8,8,9,0,7,5,0,9,8,3,8,
	1,7,5,4,6,3,7,4,6,4,9,3,9,3,1,9,2,5,5,0,6,0,4,0,0,9,2,7,7,0,1,6,7,1,1,3,9,0,0,9,8,4,8,8,2,4,0,1,2,8,5,8,3,6,1,6,0,3,5,6,3,7,0,7,
	6,6,0,1,0,4,7,1,0,1,8,1,9,4,2,9,5,5,5,9,6,1,9,8,9,4,6,7,6,7,8,3,7,4,4,9,4,4,8,2,5,5,3,7,9,7,7,4,7,2,6,8,4,7,1,0,4,0,4,7,5,3,4,6,
	4,6,2,0,8,0,4,6,6,8,4,2,5,9,0,6,9,4,9,1,2,9,3,3,1,3,6,7,7,0,2,8,9,8,9,1,5,2,1,0,4,7,5,2,1,6,2,0,5,6,9,6,6,0,2,4,0,5,8,0,3,8,1,5,
	0,1,9,3,5,1,1,2,5,3,3,8,2,4,3,0,0,3,5,5,8,7,6,4,0,2,4,7,4,9,6,4,7,3,2,6,3,9,1,4,1,9,9,2,7,2,6,0,4,2,6,9,9,2,2,7,9,6,7,8,2,3,5,4,
	7,8,1,6,3,6,0,0,9,3,4,1,7,2,1,6,4,1,2,1,9,9,2,4,5,8,6,3,1,5,0,3,0,2,8,6,1,8,2,9,7,4,5,5,5,7,0,6,7,4,9,8,3,8,5,0,5,4,9,4,5,8,8,5,
	8,6,9,2,6,9,9,5,6,9,0,9,2,7,2,1,0,7,9,7,5,0,9,3,0,2,9,5,5,3,2,1,1,6,5,3,4,4,9,8,7,2,0,2,7,5,5,9,6,0,2,3,6,4,8,0,6,6,5,4,9,9,1,1,
	9,8,8,1,8,3,4,7,9,7,7,5,3,5,6,6,3,6,9,8,0,7,4,2,6,5,4,2,5,2,7,8,6,2,5,5,1,8,1,8,4,1,7,5,7,4,6,7,2,8,9,0,9,7,7,7,7,2,7,9,3,8,0,0,
	0,8,1,6,4,7,0,6,0,0,1,6,1,4,5,2,4,9,1,9,2,1,7,3,2,1,7,2,1,4,7,7,2,3,5,0,1,4,1,4,4,1,9,7,3,5,6,8,5,4,8,1,6,1,3,6,1,1,5,7,3,5,2,5,
	5,2,1,3,3,4,7,5,7,4,1,8,4,9,4,6,8,4,3,8,5,2,3,3,2,3,9,0,7,3,9,4,1,4,3,3,3,4,5,4,7,7,6,2,4,1,6,8,6,2,5,1,8,9,8,3,5,6,9,4,8,5,5,6,
	2,0,9,9,2,1,9,2,2,2,1,8,4,2,7,2,5,5,0,2,5,4,2,5,6,8,8,7,6,7,1,7,9,0,4,9,4,6,0,1,6,5,3,4,6,6,8,0,4,9,8,8,6,2,7,2,3,2,7,9,1,7,8,6,
	0,8,5,7,8,4,3,8,3,8,2,7,9,6,7,9,7,6,6,8,1,4,5,4,1,0,0,9,5,3,8,8,3,7,8,6,3,6,0,9,5,0,6,8,0,0,6,4,2,2,5,1,2,5,2,0,5,1,1,7,3,9,2,9,
	8,4,8,9,6,0,8,4,1,2,8,4,8,8,6,2,6,9,4,5,6,0,4,2,4,1,9,6,5,2,8,5,0,2,2,2,1,0,6,6,1,1,8,6,3,0,6,7,4,4,2,7,8,6,2,2,0,3,9,1,9,4,9,4,
	5,0,4,7,1,2,3,7,1,3,7,8,6,9,6,0,9,5,6,3,6,4,3,7,1,9,1,7,2,8,7,4,6,7,7,6,4,6,5,7,5,7,3,9,6,2,4,1,3,8,9,0,8,6,5,8,3,2,6,4,5,9,9,5,
	8,1,3,3,9,0,4,7,8,0,2,7,5,9,0,0,9,9,4,6,5,7,6,4,0,7,8,9,5,1,2,6,9,4,6,8,3,9,8,3,5,2,5,9,5,7,0,9,8,2,5,8,2,2,6,2,0,5,2,2,4,8,9,4,
	0,7,7,2,6,7,1,9,4,7,8,2,6,8,4,8,2,6,0,1,4,7,6,9,9,0,9,0,2,6,4,0,1,3,6,3,9,4,4,3,7,4,5,5,3,0,5,0,6,8,2,0,3,4,9,6,2,5,2,4,5,1,7,4,
	9,3,9,9,6,5,1,4,3,1,4,2,9,8,0,9,1,9,0,6,5,9,2,5,0,9,3,7,2,2,1,6,9,6,4,6,1,5,1,5,7,0,9,8,5,8,3,8,7,4,1,0,5,9,7,8,8,5,9,5,9,7,7,2,
	9,7,5,4,9,8,9,3,0,1,6,1,7,5,3,9,2,8,4,6,8,1,3,8,2,6,8,6,8,3,8,6,8,9,4,2,7,7,4,1,5,5,9,9,1,8,5,5,9,2,5,2,4,5,9,5,3,9,5,9,4,3,1,0,
	4,9,9,7,2,5,2,4,6,8,0,8,4,5,9,8,7,2,7,3,6,4,4,6,9,5,8,4,8,6,5,3,8,3,6,7,3,6,2,2,2,6,2,6,0,9,9,1,2,4,6,0,8,0,5,1,2,4,3,8,8,4,3,9,
	0,4,5,1,2,4,4,1,3,6,5,4,9,7,6,2,7,8,0,7,9,7,7,1,5,6,9,1,4,3,5,9,9,7,7,0,0,1,2,9,6,1,6,0,8,9,4,4,1,6,9,4,8,6,8,5,5,5,8,4,8,4,0,6,
	3,5,3,4,2,2,0,7,2,2,2,5,8,2,8,4,8,8,6,4,8,1,5,8,4,5,6,0,2,8,5,0,6,0,1,6,8,4,2,7,3,9,4,5,2,2,6,7,4,6,7,6,7,8,8,9,5,2,5,2,1,3,8,5,
	2,2,5,4,9,9,5,4,6,6,6,7,2,7,8,2,3,9,8,6,4,5,6,5,9,6,1,1,6,3,5,4,8,8,6,2,3,0,5,7,7,4,5,6,4,9,8,0,3,5,5,9,3,6,3,4,5,6,8,1,7,4,3,2,
	4,1,1,2,5,1,5,0,7,6,0,6,9,4,7,9,4,5,1,0,9,6,5,9,6,0,9,4,0,2,5,2,2,8,8,7,9,7,1,0,8,9,3,1,4,5,6,6,9,1,3,6,8,6,7,2,2,8,7,4,8,9,4,0,
	5,6,0,1,0,1,5,0,3,3,0,8,6,1,7,9,2,8,6,8,0,9,2,0,8,7,4,7,6,0,9,1,7,8,2,4,9,3,8,5,8,9,0,0,9,7,1,4,9,0,9,6,7,5,9,8,5,2,6,1,3,6,5,5,
	4,9,7,8,1,8,9,3,1,2,9,7,8,4,8,2,1,6,8,2,9,9,8,9,4,8,7,2,2,6,5,8,8,0,4,8,5,7,5,6,4,0,1,4,2,7,0,4,7,7,5,5,5,1,3,2,3,7,9,6,4,1,4,5,
	1,5,2,3,7,4,6,2,3,4,3,6,4,5,4,2,8,5,8,4,4,4,7,9,5,2,6,5,8,6,7,8,2,1,0,5,1,1,4,1,3,5,4,7,3,5,7,3,9,5,2,3,1,1,3,4,2,7,1,6,6,1,0,2,
	1,3,5,9,6,9,5,3,6,2,3,1,4,4,2,9,5,2,4,8,4,9,3,7,1,8,7,1,1,0,1,4,5,7,6,5,4,0,3,5,9,0,2,7,9,9,3,4,4,0,3,7,4,2,0,0,7,3,1,0,5,7,8,5,
	3,9,0,6,2,1,9,8,3,8,7,4,4,7,8,0,8,4,7,8,4,8,9,6,8,3,3,2,1,4,4,5,7,1,3,8,6,8,7,5,1,9,4,3,5,0,6,4,3,0,2,1,8,4,5,3,1,9,1,0,4,8,4,8,
	1,0,0,5,3,7,0,6,1,4,6,8,0,6,7,4,9,1,9,2,7,8,1,9,1,1,9,7,9,3,9,9,5,2,0,6,1,4,1,9,6,6,3,4,2,8,7,5,4,4,4,0,6,4,3,7,4,5,1,2,3,7,1,8,
	1,9,2,1,7,9,9,9,8,3,9,1,0,1,5,9,1,9,5,6,1,8,1,4,6,7,5,1,4,2,6,9,1,2,3,9,7,4,8,9,4,0,9,0,7,1,8,6,4,9,4,2,3,1,9,6,1,5,6,7,9,4,5,2,
	0,8,0,9,5,1,4,6,5,5,0,2,2,5,2,3,1,6,0,3,8,8,1,9,3,0,1,4,2,0,9,3,7,6,2,1,3,7,8,5,5,9,5,6,6,3,8,9,3,7,7,8,7,0,8,3,0,3,9,0,6,9,7,9,
	2,0,7,7,3,4,6,7,2,2,1,8,2,5,6,2,5,9,9,6,6,1,5,0,1,4,2,1,5,0,3,0,6,8,0,3,8,4,4,7,7,3,4,5,4,9,2,0,2,6,0,5,4,1,4,6,6,5,9,2,5,2,0,1,
	4,9,7,4,4,2,8,5,0,7,3,2,5,1,8,6,6,6,0,0,2,1,3,2,4,3,4,0,8,8,1,9,0,7,1,0,4,8,6,3,3,1,7,3,4,6,4,9,6,5,1,4,5,3,9,0,5,7,9,6,2,6,8,5,
	6,1,0,0,5,5,0,8,1,0,6,6,5,8,7,9,6,9,9,8,1,6,3,5,7,4,7,3,6,3,8,4,0,5,2,5,7,1,4,5,9,1,0,2,8,9,7,0,6,4,1,4,0,1,1,0,9,7,1,2,0,6,2,8,
	0,4,3,9,0,3,9,7,5,9,5,1,5,6,7,7,1,5,7,7,0,0,4,2,0,3,3,7,8,6,9,9,3,6,0,0,7,2,3,0,5,5,8,7,6,3,1,7,6,3,5,9,4,2,1,8,7,3,1,2,5,1,4,7,
	1,2,0,5,3,2,9,2,8,1,9,1,8,2,6,1,8,6,1,2,5,8,6,7,3,2,1,5,7,9,1,9,8,4,1,4,8,4,8,8,2,9,1,6,4,4,7,0,6,0,9,5,7,5,2,7,0,6,9,5,7,2,2,0,
	9,1,7,5,6,7,1,1,6,7,2,2,9,1,0,9,8,1,6,9,0,9,1,5,2,8,0,1,7,3,5,0,6,7,1,2,7,4,8,5,8,3,2,2,2,8,7,1,8,3,5,2,0,9,3,5,3,9,6,5,7,2,5,1,
	2,1,0,8,3,5,7,9,1,5,1,3,6,9,8,8,2,0,9,1,4,4,4,2,1,0,0,6,7,5,1,0,3,3,4,6,7,1,1,0,3,1,4,1,2,6,7,1,1,1,3,6,9,9,0,8,6,5,8,5,1,6,3,9,
	8,3,1,5,0,1,9,7,0,1,6,5,1,5,1,1,6,8,5,1,7,1,4,3,7,6,5,7,6,1,8,3,5,1,5,5,6,5,0,8,8,4,9,0,9,9,8,9,8,5,9,9,8,2,3,8,7,3,4,5,5,2,8,3,
	3,1,6,3,5,5,0,7,6,4,7,9,1,8,5,3,5,8,9,3,2,2,6,1,8,5,4,8,9,6,3,2,1,3,2,9,3,3,0,8,9,8,5,7,0,6,4,2,0,4,6,7,5,2,5,9,0,7,0,9,1,5,4,8,
	1,4,1,6,5,4,9,8,5,9,4,6,1,6,3,7,1,8,0,2,7,0,9,8,1,9,9,4,3,0,9,9,2,4,4,8,8,9,5,7,5,7,1,2,8,2,8,9,0,5,9,2,3,2,3,3,2,6,0,9,7,2,9,9,
	7,1,2,0,8,4,4,3,3,5,7,3,2,6,5,4,8,9,3,8,2,3,9,1,1,9,3,2,5,9,7,4,6,3,6,6,7,3,0,5,8,3,6,0,4,1,4,2,8,1,3,8,8,3,0,3,2,0,3,8,2,4,9,0,
	3,7,5,8,9,8,5,2,4,3,7,4,4,1,7,0,2,9,1,3,2,7,6,5,6,1,8,0,9,3,7,7,3,4,4,4,0,3,0,7,0,7,4,6,9,2,1,1,2,0,1,9,1,3,0,2,0,3,3,0,3,8,0,1,
	9,7,6,2,1,1,0,1,1,0,0,4,4,9,2,9,3,2,1,5,1,6,0,8,4,2,4,4,4,8,5,9,6,3,7,6,6,9,8,3,8,9,5,2,2,8,6,8,4,7,8,3,1,2,3,5,5,2,6,5,8,2,1,3,
	1,4,4,9,5,7,6,8,5,7,2,6,2,4,3,3,4,4,1,8,9,3,0,3,9,6,8,6,4,2,6,2,4,3,4,1,0,7,7,3,2,2,6,9,7,8,0,2,8,0,7,3,1,8,9,1,5,4,4,1,1,0,1,0,
	4,4,6,8,2,3,2,5,2,7,1,6,2,0,1,0,5,2,6,5,2,2,7,2,1,1,1,6,6,0,3,9,6,6,6,5,5,7,3,0,9,2,5,4,7,1,1,0,5,5,7,8,5,3,7,6,3,4,6,6,8,2,0,6,
	5,3,1,0,9,8,9,6,5,2,6,9,1,8,6,2,0,5,6,4,7,6,9,3,1,2,5,7,0,5,8,6,3,5,6,6,2,0,1,8,5,5,8,1,0,0,7,2,9,3,6,0,6,5,9,8,7,6,4,8,6,1,1,7,
	9,1,0,4,5,3,3,4,8,8,5,0,3,4,6,1,1,3,6,5,7,6,8,6,7,5,3,2,4,9,4,4,1,6,6,8,0,3,9,6,2,6,5,7,9,7,8,7,7,1,8,5,5,6,0,8,4,5,5,2,9,6,5,4,
	1,2,6,6,5,4,0,8,5,3,0,6,1,4,3,4,4,4,3,1,8,5,8,6,7,6,9,7,5,1,4,5,6,6,1,4,0,6,8,0,0,7,0,0,2,3,7,8,7,7,6,5,9,1,3,4,4,0,1,7,1,2,7,4,
	9,4,7,0,4,2,0,5,6,2,2,3,0,5,3,8,9,9,4,5,6,1,3,1,4,0,7,1,1,2,7,0,0,0,4,0,7,8,5,4,7,3,3,2,6,9,9,3,9,0,8,1,4,5,4,6,6,4,6,4,5,8,8,0,
	7,9,7,2,7,0,8,2,6,6,8,3,0,6,3,4,3,2,8,5,8,7,8,5,6,9,8,3,0,5,2,3,5,8,0,8,9,3,3,0,6,5,7,5,7,4,0,6,7,9,5,4,5,7,1,6,3,7,7,5,2,5,4,2,
	0,2,1,1,4,9,5,5,7,6,1,5,8,1,4,0,0,2,5,0,1,2,6,2,2,8,5,9,4,1,3,0,2,1,6,4,7,1,5,5,0,9,7,9,2,5,9,2,3,0,9,9,0,7,9,6,5,4,7,3,7,6,1,2,
	5,5,1,7,6,5,6,7,5,1,3,5,7,5,1,7,8,2,9,6,6,6,4,5,4,7,7,9,1,7,4,5,0,1,1,2,9,9,6,1,4,8,9,0,3,0,4,6,3,9,9,4,7,1,3,2,9,6,2,1,0,7,3,4,
	0,4,3,7,5,1,8,9,5,7,3,5,9,6,1,4,5,8,9,0,1,9,3,8,9,7,1,3,1,1,1,7,9,0,4,2,9,7,8,2,8,5,6,4,7,5,0,3,2,0,3,1,9,8,6,9,1,5,1,4,0,2,8,7,
	0,8,0,8,5,9,9,0,4,8,0,1,0,9,4,1,2,1,4,7,2,2,1,3,1,7,9,4,7,6,4,7,7,7,2,6,2,2,4,1,4,2,5,4,8,5,4,5,4,0,3,3,2,1,5,7,1,8,5,3,0,6,1,4,
	2,2,8,8,1,3,7,5,8,5,0,4,3,0,6,3,3,2,1,7,5,1,8,2,9,7,9,8,6,6,2,2,3,7,1,7,2,1,5,9,1,6,0,7,7,1,6,6,9,2,5,4,7,4,8,7,3,8,9,8,6,6,5,4,
	9,4,9,4,5,0,1,1,4,6,5,4,0,6,2,8,4,3,3,6,6,3,9,3,7,9,0,0,3,9,7,6,9,2,6,5,6,7,2,1,4,6,3,8,5,3,0,6,7,3,6,0,9,6,5,7,1,2,0,9,1,8,0,7,
	6,3,8,3,2,7,1,6,6,4,1,6,2,7,4,8,8,8,8,0,0,7,8,6,9,2,5,6,0,2,9,0,2,2,8,4,7,2,1,0,4,0,3,1,7,2,1,1,8,6,0,8,2,0,4,1,9,0,0,0,4,2,2,9,
	6,6,1,7,1,1,9,6,3,7,7,9,2,1,3,3,7,5,7,5,1,1,4,9,5,9,5,0,1,5,6,6,0,4,9,6,3,1,8,6,2,9,4,7,2,6,5,4,7,3,6,4,2,5,2,3,0,8,1,7,7,0,3,6,
	7,5,1,5,9,0,6,7,3,5,0,2,3,5,0,7,2,8,3,5,4,0,5,6,7,0,4,0,3,8,6,7,4,3,5,1,3,6,2,2,2,2,4,7,7,1,5,8,9,1,5,0,4,9,5,3,0,9,8,4,4,4,8,9,
	3,3,3,0,9,6,3,4,0,8,7,8,0,7,6,9,3,2,5,9,9,3,9,7,8,0,5,4,1,9,3,4,1,4,4,7,3,7,7,4,4,1,8,4,2,6,3,1,2,9,8,6,0,8,0,9,9,8,8,8,6,8,7,4,
	1,3,2,6,0,4,7,2,1,5,6,9,5,1,6,2,3,9,6,5,8,6,4,5,7,3,0,2,1,6,3,1,5,9,8,1,9,3,1,9,5,1,6,7,3,5,3,8,1,2,9,7,4,1,6,7,7,2,9,4,7,8,6,7,
	2,4,2,2,9,2,4,6,5,4,3,6,6,8,0,0,9,8,0,6,7,6,9,2,8,2,3,8,2,8,0,6,8,9,9,6,4,0,0,4,8,2,4,3,5,4,0,3,7,0,1,4,1,6,3,1,4,9,6,5,8,9,7,9,
	4,0,9,2,4,3,2,3,7,8,9,6,9,0,7,0,6,9,7,7,9,4,2,2,3,6,2,5,0,8,2,2,1,6,8,8,9,5,7,3,8,3,7,9,8,6,2,3,0,0,1,5,9,3,7,7,6,4,7,1,6,5,1,2,
	2,8,9,3,5,7,8,6,0,1,5,8,8,1,6,1,7,5,5,7,8,2,9,7,3,5,2,3,3,4,4,6,0,4,2,8,1,5,1,2,6,2,7,2,0,3,7,3,4,3,1,4,6,5,3,1,9,7,7,7,7,4,1,6,
	0,3,1,9,9,0,6,6,5,5,4,1,8,7,6,3,9,7,9,2,9,3,3,4,4,1,9,5,2,1,5,4,1,3,4,1,8,9,9,4,8,5,4,4,4,7,3,4,5,6,7,3,8,3,1,6,2,4,9,9,3,4,1,9,
	1,3,1,8,1,4,8,0,9,2,7,7,7,7,1,0,3,8,6,3,8,7,7,3,4,3,1,7,7,2,0,7,5,4,5,6,5,4,5,3,2,2,0,7,7,7,0,9,2,1,2,0,1,9,0,5,1,6,6,0,9,6,2,8,
	0,4,9,0,9,2,6,3,6,0,1,9,7,5,9,8,8,2,8,1,6,1,3,3,2,3,1,6,6,6,3,6,5,2,8,6,1,9,3,2,6,6,8,6,3,3,6,0,6,2,7,3,5,6,7,6,3,0,3,5,4,4,7,7,
	6,2,8,0,3,5,0,4,5,0,7,7,7,2,3,5,5,4,7,1,0,5,8,5,9,5,4,8,7,0,2,7,9,0,8,1,4,3,5,6,2,4,0,1,4,5,1,7,1,8,0,6,2,4,6,4,3,6,2,6,7,9,4,5,
	6,1,2,7,5,3,1,8,1,3,4,0,7,8,3,3,0,3,3,6,2,5,4,2,3,2,7,8,3,9,4,4,9,7,5,3,8,2,4,3,7,2,0,5,8,3,5,3,1,1,4,7,7,1,1,9,9,2,6,0,6,3,8,1,
	3,3,4,6,7,7,6,8,7,9,6,9,5,9,7,0,3,0,9,8,3,3,9,1,3,0,7,7,1,0,9,8,7,0,4,0,8,5,9,1,3,3,7,4,6,4,1,4,4,2,8,2,2,7,7,2,6,3,4,6,5,9,4,7,
	0,4,7,4,5,8,7,8,4,7,7,8,7,2,0,1,9,2,7,7,1,5,2,8,0,7,3,1,7,6,7,9,0,7,7,0,7,1,5,7,2,1,3,4,4,4,7,3,0,6,0,5,7,0,0,7,3,3,4,9,2,4,3,6,
	9,3,1,1,3,8,3,5,0,4,9,3,1,6,3,1,2,8,4,0,4,2,5,1,2,1,9,2,5,6,5,1,7,9,8,0,6,9,4,1,1,3,5,2,8,0,1,3,1,4,7,0,1,3,0,4,7,8,1,6,4,3,7,8,
	8,5,1,8,5,2,9,0,9,2,8,5,4,5,2,0,1,1,6,5,8,3,9,3,4,1,9,6,5,6,2,1,3,4,9,1,4,3,4,1,5,9,5,6,2,5,8,6,5,8,6,5,5,7,0,5,5,2,6,9,0,4,9,6,
	5,2,0,9,8,5,8,0,3,3,8,5,0,7,2,2,4,2,6,4,8,2,9,3,9,7,2,8,5,8,4,7,8,3,1,6,3,0,5,7,7,7,7,5,6,0,6,8,8,8,7,6,4,4,6,2,4,8,2,4,6,8,5,7,
	9,2,6,0,3,9,5,3,5,2,7,7,3,4,8,0,3,0,4,8,0,2,9,0,0,5,8,7,6,0,7,5,8,2,5,1,0,4,7,4,7,0,9,1,6,4,3,9,6,1,3,6,2,6,7,6,0,4,4,9,2,5,6,2,
	7,4,2,0,4,2,0,8,3,2,0,8,5,6,6,1,1,9,0,6,2,5,4,5,4,3,3,7,2,1,3,1,5,3,5,9,5,8,4,5,0,6,8,7,7,2,4,6,0,2,9,0,1,6,1,8,7,6,6,7,9,5,2,4,
	0,6,1,6,3,4,2,5,2,2,5,7,7,1,9,5,4,2,9,1,6,2,9,9,1,9,3,0,6,4,5,5,3,7,7,9,9,1,4,0,3,7,3,4,0,4,3,2,8,7,5,2,6,2,8,8,8,9,6,3,9,9,5,8,
	7,9,4,7,5,7,2,9,1,7,4,6,4,2,6,3,5,7,4,5,5,2,5,4,0,7,9,0,9,1,4,5,1,3,5,7,1,1,1,3,6,9,4,1,0,9,1,1,9,3,9,3,2,5,1,9,1,0,7,6,0,2,0,8,
	2,5,2,0,2,6,1,8,7,9,8,5,3,1,8,8,7,7,0,5,8,4,2,9,7,2,5,9,1,6,7,7,8,1,3,1,4,9,6,9,9,0,0,9,0,1,9,2,1,1,6,9,7,1,7,3,7,2,7,8,4,7,6,8,
	4,7,2,6,8,6,0,8,4,9,0,0,3,3,7,7,0,2,4,2,4,2,9,1,6,5,1,3,0,0,5,0,0,5,1,6,8,3,2,3,3,6,4,3,5,0,3,8,9,5,1,7,0,2,9,8,9,3,9,2,2,3,3,4,
	5,1,7,2,2,0,1,3,8,1,2,8,0,6,9,6,5,0,1,1,7,8,4,4,0,8,7,4,5,1,9,6,0,1,2,1,2,2,8,5,9,9,3,7,1,6,2,3,1,3,0,1,7,1,1,4,4,4,8,4,6,4,0,9,
	0,3,8,9,0,6,4,4,9,5,4,4,4,0,0,6,1,9,8,6,9,0,7,5,4,8,5,1,6,0,2,6,3,2,7,5,0,5,2,9,8,3,4,9,1,8,7,4,0,7,8,6,6,8,0,8,8,1,8,3,3,8,5,1,
	0,2,2,8,3,3,4,5,0,8,5,0,4,8,6,0,8,2,5,0,3,9,3,0,2,1,3,3,2,1,9,7,1,5,5,1,8,4,3,0,6,3,5,4,5,5,0,0,7,6,6,8,2,8,2,9,4,9,3,0,4,1,3,7,
	7,6,5,5,2,7,9,3,9,7,5,1,7,5,4,6,1,3,9,5,3,9,8,4,6,8,3,3,9,3,6,3,8,3,0,4,7,4,6,1,1,9,9,6,6,5,3,8,5,8,1,5,3,8,4,2,0,5,6,8,5,3,3,8,
	6,2,1,8,6,7,2,5,2,3,3,4,0,2,8,3,0,8,7,1,1,2,3,2,8,2,7,8,9,2,1,2,5,0,7,7,1,2,6,2,9,4,6,3,2,2,9,5,6,3,9,8,9,8,9,8,9,3,5,8,2,1,1,6,
	7,4,5,6,2,7,0,1,0,2,1,8,3,5,6,4,6,2,2,0,1,3,4,9,6,7,1,5,1,8,8,1,9,0,9,7,3,0,3,8,1,1,9,8,0,0,4,9,7,3,4,0,7,2,3,9,6,1,0,3,6,8,5,4,
	0,6,6,4,3,1,9,3,9,5,0,9,7,9,0,1,9,0,6,9,9,6,3,9,5,5,2,4,5,3,0,0,5,4,5,0,5,8,0,6,8,5,5,0,1,9,5,6,7,3,0,2,2,9,2,1,9,1,3,9,3,3,9,1,
	8,5,6,8,0,3,4,4,9,0,3,9,8,2,0,5,9,5,5,1,0,0,2,2,6,3,5,3,5,3,6,1,9,2,0,4,1,9,9,4,7,4,5,5,3,8,5,9,3,8,1,0,2,3,4,3,9,5,5,4,4,9,5,9,
	7,7,8,3,7,7,9,0,2,3,7,4,2,1,6,1,7,2,7,1,1,1,7,2,3,6,4,3,4,3,5,4,3,9,4,7,8,2,2,1,8,1,8,5,2,8,6,2,4,0,8,5,1,4,0,0,6,6,6,0,4,4,3,3,
	2,5,8,8,8,5,6,9,8,6,7,0,5,4,3,1,5,4,7,0,6,9,6,5,7,4,7,4,5,8,5,5,0,3,3,2,3,2,3,3,4,2,1,0,7,3,0,1,5,4,5,9,4,0,5,1,6,5,5,3,7,9,0,6,
	8,6,6,2,7,3,3,3,7,9,9,5,8,5,1,1,5,6,2,5,7,8,4,3,2,2,9,8,8,2,7,3,7,2,3,1,9,8,9,8,7,5,7,1,4,1,5,9,5,7,8,1,1,1,9,6,3,5,8,3,3,0,0,5,
	9,4,0,8,7,3,0,6,8,1,2,1,6,0,2,8,7,6,4,9,6,2,8,6,7,4,4,6,0,4,7,7,4,6,4,9,1,5,9,9,5,0,5,4,9,7,3,7,4,2,5,6,2,6,9,0,1,0,4,9,0,3,7,7,
	8,1,9,8,6,8,3,5,9,3,8,1,4,6,5,7,4,1,2,6,8,0,4,9,2,5,6,4,8,7,9,8,5,5,6,1,4,5,3,7,2,3,4,7,8,6,7,3,3,0,3,9,0,4,6,8,8,3,8,3,4,3,6,3,
	4,6,5,5,3,7,9,4,9,8,6,4,1,9,2,7,0,5,6,3,8,7,2,9,3,1,7,4,8,7,2,3,3,2,0,8,3,7,6,0,1,1,2,3,0,2,9,9,1,1,3,6,7,9,3,8,6,2,7,0,8,9,4,3,
	8,7,9,9,3,6,2,0,1,6,2,9,5,1,5,4,1,3,3,7,1,4,2,4,8,9,2,8,3,0,7,2,2,0,1,2,6,9,0,1,4,7,5,4,6,6,8,4,7,6,5,3,5,7,6,1,6,4,7,7,3,7,9,4,
	6,7,5,2,0,0,4,9,0,7,5,7,1,5,5,5,2,7,8,1,9,6,5,3,6,2,1,3,2,3,9,2,6,4,0,6,1,6,0,1,3,6,3,5,8,1,5,5,9,0,7,4,2,2,0,2,0,2,0,3,1,8,7,2,
	7,7,6,0,5,2,7,7,2,1,9,0,0,5,5,6,1,4,8,4,2,5,5,5,1,8,7,9,2,5,3,0,3,4,3,5,1,3,9,8,4,4,2,5,3,2,2,3,4,1,5,7,6,2,3,3,6,1,0,6,4,2,5,0,
	6,3,9,0,4,9,7,5,0,0,8,6,5,6,2,7,1,0,9,5,3,5,9,1,9,4,6,5,8,9,7,5,1,4,1,3,1,0,3,4,8,2,2,7,6,9,3,0,6,2,4,7,4,3,5,3,6,3,2,5,6,9,1,6,
	0,7,8,1,5,4,7,8,1,8,1,1,5,2,8,4,3,6,6,7,9,5,7,0,6,1,1,0,8,6,1,5,3,3,1,5,0,4,4,5,2,1,2,7,4,7,3,9,2,4,5,4,4,9,4,5,4,2,3,6,8,2,8,8,
	6,0,6,1,3,4,0,8,4,1,4,8,6,3,7,7,6,7,0,0,9,6,1,2,0,7,1,5,1,2,4,9,1,4,0,4,3,0,2,7,2,5,3,8,6,0,7,6,4,8,2,3,6,3,4,1,4,3,3,4,6,2,3,5,
	1,8,9,7,5,7,6,6,4,5,2,1,6,4,1,3,7,6,7,9,6,9,0,3,1,4,9,5,0,1,9,1,0,8,5,7,5,9,8,4,4,2,3,9,1,9,8,6,2,9,1,6,4,2,1,9,3,9,9,4,9,0,7,2,
	3,6,2,3,4,6,4,6,8,4,4,1,1,7,3,9,4,0,3,2,6,5,9,1,8,4,0,4,4,3,7,8,0,5,1,3,3,3,8,9,4,5,2,5,7,4,2,3,9,9,5,0,8,2,9,6,5,9,1,2,2,8,5,0,
	8,5,5,5,8,2,1,5,7,2,5,0,3,1,0,7,1,2,5,7,0,1,2,6,6,8,3,0,2,4,0,2,9,2,9,5,2,5,2,2,0,1,1,8,7,2,6,7,6,7,5,6,2,2,0,4,1,5,4,2,0,5,1,6,
	1,8,4,1,6,3,4,8,4,7,5,6,5,1,6,9,9,9,8,1,1,6,1,4,1,0,1,0,0,2,9,9,6,0,7,8,3,8,6,9,0,9,2,9,1,6,0,3,0,2,8,8,4,0,0,2,6,9,1,0,4,1,4,0,
	7,9,2,8,8,6,2,1,5,0,7,8,4,2,4,5,1,6,7,0,9,0,8,7,0,0,0,6,9,9,2,8,2,1,2,0,6,6,0,4,1,8,3,7,1,8,0,6,5,3,5,5,6,7,2,5,2,5,3,2,5,6,7,5,
	3,2,8,6,1,2,9,1,0,4,2,4,8,7,7,6,1,8,2,5,8,2,9,7,6,5,1,5,7,9,5,9,8,4,7,0,3,5,6,2,2,2,6,2,9,3,4,8,6,0,0,3,4,1,5,8,7,2,2,9,8,0,5,3,
	4,9,8,9,6,5,0,2,2,6,2,9,1,7,4,8,7,8,8,2,0,2,7,3,4,2,0,9,2,2,2,2,4,5,3,3,9,8,5,6,2,6,4,7,6,6,9,1,4,9,0,5,5,6,2,8,4,2,5,0,3,9,1,2,
	7,5,7,7,1,0,2,8,4,0,2,7,9,9,8,0,6,6,3,6,5,8,2,5,4,8,8,9,2,6,4,8,8,0,2,5,4,5,6,6,1,0,1,7,2,9,6,7,0,2,6,6,4,0,7,6,5,5,9,0,4,2,9,0,
	9,9,4,5,6,8,1,5,0,6,5,2,6,5,3,0,5,3,7,1,8,2,9,4,1,2,7,0,3,3,6,9,3,1,3,7,8,5,1,7,8,6,0,9,0,4,0,7,0,8,6,6,7,1,1,4,9,6,5,5,8,3,4,3,
	4,3,4,7,6,9,3,3,8,5,7,8,1,7,1,1,3,8,6,4,5,5,8,7,3,6,7,8,1,2,3,0,1,4,5,8,7,6,8,7,1,2,6,6,0,3,4,8,9,1,3,9,0,9,5,6,2,0,0,9,9,3,9,3,
	};
	const char* test_info_str[] = {"Show input-index-scramblings-needed for DIF and DIT","DIF","DIT","Combined DIF+DIT"};
	double *a = 0x0, *b = 0x0, *arrtmp = 0x0, *ptmp = 0x0;
	struct complex *ac, *bc;
	struct complex **mat = 0x0, **matp = 0x0, **ctmpp = 0x0, *ctmp = 0x0;
	double t0,t1,t2,t3;
	double theta, twopi = 6.2831853071795864769;
  #ifdef USE_FGT61
	#if defined(USE_SSE2) || (RADIX != 16)
	  #error Currently only Scalar-mode Radix-16 supported for FGT-test mode!
	#endif
	const uint64 q  = 0x1FFFFFFFFFFFFFFFull;	// q = 2^61 - 1
	uint64 *amod = 0x0, *bmod = 0x0, *iptr = 0x0;
	uint128 *am, *bm;
	uint128 **matmod = 0x0, **matmodp = 0x0, **itmpp = 0x0, *itmp = 0x0;
	uint64 order,root_re,root_im,rm,im,m0,m1,m2,m3;
  #endif

	/********* allocate all radix-dependent arrays dynamically: ********/
	index        = ALLOC_INT(index       , RADIX);
	dit_scramble = ALLOC_INT(dit_scramble, RADIX);
	/* double a[rmul*RADIX], b[rmul*RADIX], arrtmp[rmul*RADIX]: */
	ptmp = ALLOC_DOUBLE(ptmp, rmul*RADIX);	ASSERT((ptmp != 0x0), "FATAL: unable to allocate array A in test_fft_radix.\n");
	a    = ALIGN_DOUBLE(ptmp);	ptmp = 0x0;
	ac = (struct complex *)a;
	ASSERT(((long)((void *)a) & 63) == 0x0,"test_fft_radix: A[] not aligned on 64-byte boundary!");
	ptmp = ALLOC_DOUBLE(ptmp, rmul*RADIX);	ASSERT((ptmp != 0x0), "FATAL: unable to allocate array B in test_fft_radix.\n");
	b    = ALIGN_DOUBLE(ptmp);	ptmp = 0x0;
	ASSERT(((long)((void *)b) & 63) == 0x0,"test_fft_radix: B[] not aligned on 64-byte boundary!");
	bc = (struct complex *)b;
	ptmp = ALLOC_DOUBLE(ptmp, rmul*RADIX);	ASSERT((ptmp != 0x0), "FATAL: unable to allocate array A_ptmp in test_fft_radix.\n");
	arrtmp = ALIGN_DOUBLE(ptmp);	ptmp = 0x0;
	ASSERT(((long)((void *)arrtmp) & 63) == 0x0,"test_fft_radix: arrtmp[] not aligned on 64-byte boundary!");
	/* struct complex mat[radix][RADIX], *matp[RADIX]: */
	ctmpp = ALLOC_POINTER(ctmpp,struct complex*, RADIX);	ASSERT((ctmpp != 0x0), "FATAL: unable to allocate array MATP in test_fft_radix.\n");
	matp  = ALIGN_POINTER(ctmpp,struct complex*);
	ctmpp = ALLOC_POINTER(ctmpp,struct complex*, RADIX);	ASSERT((ctmpp != 0x0), "FATAL: unable to allocate array MAT[][] in test_fft_radix.\n");
	mat   = ALIGN_POINTER(ctmpp,struct complex*);
	for(i = 0; i < RADIX; ++i) {
		ctmp = ALLOC_COMPLEX(ctmp, RADIX);	ASSERT((ctmp != 0x0), "FATAL: unable to allocate array Ctmp in test_fft_radix.\n");
		mat[i] = ALIGN_COMPLEX(ctmp);
		ctmp = 0x0;	/* Must re-init pointer so the realloc used by the ALLOC macro allocates new fresh memory for each row */
	}
  #ifdef USE_FGT61
	iptr = ALLOC_UINT64(iptr, rmul*RADIX);	ASSERT((iptr != 0x0), "FATAL: unable to allocate array AMOD in test_fft_radix.\n");
	amod = ALIGN_UINT64(iptr);	iptr = 0x0;
	am = (uint128 *)amod;
	ASSERT(((long)((void *)amod) & 63) == 0x0,"test_fft_radix: AMOD[] not aligned on 64-byte boundary!");
	iptr = ALLOC_UINT64(iptr, rmul*RADIX);	ASSERT((iptr != 0x0), "FATAL: unable to allocate array BMOD in test_fft_radix.\n");
	bmod = ALIGN_UINT64(iptr);	iptr = 0x0;
	ASSERT(((long)((void *)bmod) & 63) == 0x0,"test_fft_radix: BMOD[] not aligned on 64-byte boundary!");
	bm = (uint128 *)bmod;
	iptr = ALLOC_UINT64(iptr, rmul*RADIX);	ASSERT((iptr != 0x0), "FATAL: unable to allocate array A_iptr in test_fft_radix.\n");
	itmpp = ALLOC_POINTER(itmpp,uint128*, RADIX);	ASSERT((itmpp != 0x0), "FATAL: unable to allocate array MATP in test_fft_radix.\n");
	matmodp  = ALIGN_POINTER(itmpp,uint128*);
	itmpp = ALLOC_POINTER(itmpp,uint128*, RADIX);	ASSERT((itmpp != 0x0), "FATAL: unable to allocate array MAT[][] in test_fft_radix.\n");
	matmod   = ALIGN_POINTER(itmpp,uint128*);
	for(i = 0; i < RADIX; ++i) {
		itmp = ALLOC_UINT128(itmp, RADIX);	ASSERT((itmp != 0x0), "FATAL: unable to allocate array Ctmp in test_fft_radix.\n");
		matmod[i] = ALIGN_UINT128(itmp);
		itmp = 0x0;	/* Must re-init pointer so the realloc used by the ALLOC macro allocates new fresh memory for each row */
	}
  #endif

	fprintf(stderr, "test_fft_radix: Testing radix-%d %s dft:\n", RADIX, test_info_str[TTYPE]);

	/* Power-of-2 component of the DFT length: */
	pow2 = 1 << trailz32(RADIX);
	podd = RADIX >> trailz32(RADIX);
	ASSERT(RADIX == pow2*podd, "Radix decomposition failed!");
	ASSERT((podd < 16 || podd == 31 || podd == 63), "test_fft_radix: Illegal radix; must be odd*2^n with odd = [3,5,7,9,11,13,15,31,63]");
	/* These may not have been init'ed yet, so do it here: */
	DAT_BITS = DAT_BITS_DEF;
	PAD_BITS = PAD_BITS_DEF;

	/* Init data array in scalar-layout mode for reference matrix-multiply-DFT computation: */
#if RADIX == 1024
	t0 = t1 = 0.;
	for(j = 0; j < 32 ; j++)
	{
		t2 = t3 = 0.;
		for(i = (j << 5); i < (j << 5) + 32 ; i++)
		{
			a[2*i  ] = ref[2*i  ];	t2 += ref[2*i  ];
			a[2*i+1] = ref[2*i+1];	t3 += ref[2*i+1];
		//	printf("Input %03x = %15.5f  %15.5f\n",i, a[2*i  ],a[2*i+1]);
		}
	//	printf("Block %03x sum[Re,Im] = %15.5f  %15.5f\n",j, t2,t3);
		t0 += t2; t1 += t3;
	}
	printf("DC signal components: sum[Re,Im] = %15.5f  %15.5f\n",t0,t1);
#else
	t0 = t1 = 0.;
	for(i = 0; i < RADIX ; i++)
	{
		a[2*i  ] = ref[2*i  ];	t0 += ref[2*i  ];
		a[2*i+1] = ref[2*i+1];	t1 += ref[2*i+1];
	#ifdef USE_FGT61
		amod[2*i  ] = ref[2*i  ];
		amod[2*i+1] = ref[2*i+1];
	#endif
	}
	printf("DC signal components: sum[Re,Im] = %15.5f  %15.5f\n",t0,t1);
#endif
	/* Init DFT matrix */
#ifdef USE_FGT61
	order = RADIX;	prim_root_q(order, &root_re,&root_im);	// RADIXth primitive root of unity
	// primitive 16th root of unity, scaled by *8:
	ASSERT(root_re == 1693317751237720973ull && root_im == 2283815672160731785ull,"Bad prim-root[16]!");;
#endif
	for(i = 0; i < RADIX; i++)
	{
		theta = i*twopi/RADIX;
	#ifdef USE_FGT61
		pow_modq((uint64)i, root_re,root_im, &m0,&m1);	// m0,m1 = Ith power of prim-root
		if(i == 0) ASSERT(m0 == 1ull && m1 == 0ull, "Bad 0th power of prim-root!");
		rm = 1ull;	im = 0ull;	// leftmost col has [m0,m1]^0 = [1,0]...
	//	printf("DFT-int matrix row %d:\n",i);
	#endif
		for(j = 0; j < RADIX; j++)
		{
			mat[i][j].re = cos(j*theta);
			mat[i][j].im = sin(j*theta);
		/*	printf("mat[%4d][%4d] = %15.10f  %15.10f\n",i,j, mat[i][j].re, mat[i][j].im);*/
		#ifdef USE_FGT61
			matmod[i][j].d0 = rm;
			matmod[i][j].d1 = im;
	//	printf("\t[%2d] = %20llu, %20llu\n",j, rm,im);
			cmul_modq(m0,m1, rm,im, &rm,&im);	// ... [j]col has [m0,m1]^j
			rm = qreduce_full(rm);	im = qreduce_full(im);
		#endif
		}
	}
	matmul_complex(mat,ac,bc,RADIX,RADIX);
#ifdef USE_FGT61
	matmul_fgtmod(matmod,am,bm,RADIX,RADIX);
#endif

#ifdef USE_SSE2
	/* In SSE2 mode re-Init data array, using [re,re,im,im] data layout: */
	i = RADIX-1;	i = i + ( (i >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
	ASSERT(i == RADIX-1, "for large radix, need to enable array padding in indexing here!!");
	ASSERT(rmul == 4, "!");
	for(i = 0; i < RADIX ; i++)
	{
		a[ 2*i   *RE_IM_STRIDE  ] = ref[2*i  ];
		a[ 2*i   *RE_IM_STRIDE+1] = ref[2*i  ];
		a[(2*i+1)*RE_IM_STRIDE  ] = ref[2*i+1];
		a[(2*i+1)*RE_IM_STRIDE+1] = ref[2*i+1];
	}
#else
	// Reindex a[]-array data using padded-index scheme prior to calling FFT code:
	for(i = 0; i < RADIX+RADIX ; i += 2)
	{
		j1 = i + ( (i >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
		a[j1  ] = ref[i  ];
		a[j1+1] = ref[i+1];
	#ifdef USE_FGT61
		amod[j1  ] = ref[i  ];
		amod[j1+1] = ref[i+1];
	#endif
	}
#endif

	// If doing a DIF followed by a DIT, save a copy of the original data:
#if TTYPE == 3
	for(i = 0; i < RADIX ; i++)
	{
		arrtmp[2*i  ] = ref[2*i  ];	// These data are whole-number, thus no need
		arrtmp[2*i+1] = ref[2*i+1];	// for a separate brrtmp-array in the FGT case
	}
#endif

/*...Forward (DIF) FFT sincos data are in bit-reversed order.	*/

	l =0;

	switch(RADIX){
	case 2 :
		nradices = 1;
		radix_prim[l++] = 2; break;
	case 3 :
		nradices = 1;
		radix_prim[l++] = 3; break;
	case 4 :
		nradices = 2;
		radix_prim[l++] = 2; radix_prim[l++] = 2; break;
	case 5 :
		nradices = 1;
		radix_prim[l++] = 5; break;
	case 6 :
		nradices = 2;
		radix_prim[l++] = 3; radix_prim[l++] = 2; break;
	case 7 :
		nradices = 1;
		radix_prim[l++] = 7; break;
	case 8 :
		nradices = 3;
		radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
	case 9 :
		nradices = 2;
		radix_prim[l++] = 3; radix_prim[l++] = 3; break;
	case 10 :
		nradices = 2;
		radix_prim[l++] = 5; radix_prim[l++] = 2; break;
	case 11 :
		nradices = 1;
		radix_prim[l++] =11; break;
	case 12 :
		nradices = 3;
		radix_prim[l++] = 3; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
	case 13 :
		nradices = 1;
		radix_prim[l++] =13; break;
	case 14 :
		nradices = 2;
		radix_prim[l++] = 7; radix_prim[l++] = 2; break;
	case 15 :
		nradices = 2;
		radix_prim[l++] = 5; radix_prim[l++] = 3; break;
	case 16 :
		nradices = 4;
		radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
	case 18 :
		nradices = 3;
		radix_prim[l++] = 3; radix_prim[l++] = 3; radix_prim[l++] = 2; break;
	case 20 :
		nradices = 3;
		radix_prim[l++] = 5; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
	case 22 :
		nradices = 2;
		radix_prim[l++] =11; radix_prim[l++] = 2; break;
	case 24 :
		nradices = 4;
		radix_prim[l++] = 3; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
	case 26 :
		nradices = 2;
		radix_prim[l++] =13; radix_prim[l++] = 2; break;
	case 28 :
		nradices = 3;
		radix_prim[l++] = 7; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
	case 30 :
		nradices = 3;
		radix_prim[l++] = 5; radix_prim[l++] = 3; radix_prim[l++] = 2; break;
	case 31 :
		nradices = 1;
		radix_prim[l++] =31; break;
	case 32 :
		nradices = 5;
		radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
	case 36 :
		nradices = 4;
		radix_prim[l++] = 3; radix_prim[l++] = 3; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
	case 40:
		nradices = 4;
		radix_prim[l++] = 5; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
	case 44 :
		nradices = 3;
		radix_prim[l++] =11; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
	case 48 :
		nradices = 5;
		radix_prim[l++] = 3; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
	case 52 :
		nradices = 3;
		radix_prim[l++] =13; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
	case 56 :
		nradices = 4;
		radix_prim[l++] = 7; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
	case 60 :
		nradices = 4;
		radix_prim[l++] = 5; radix_prim[l++] = 3; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
	case 63 :
		nradices = 3;
		radix_prim[l++] = 7; radix_prim[l++] = 3; radix_prim[l++] = 3; break;
//		radix_prim[l++] = 3; radix_prim[l++] = 3; radix_prim[l++] = 7; break;
	case 64 :
		nradices = 6;
		radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
	case 72 :
		nradices = 5;
		radix_prim[l++] = 3; radix_prim[l++] = 3; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
	case 80 :
		nradices = 5;
		radix_prim[l++] = 5; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
	case 88 :
		nradices = 4;
		radix_prim[l++] =11; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
	case 96 :
		nradices = 6;
		radix_prim[l++] = 3; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
	case 104:
		nradices = 4;
		radix_prim[l++] =13; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
	case 112:
		nradices = 5;
		radix_prim[l++] = 7; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
	case 120:
		nradices = 5;
		radix_prim[l++] = 5; radix_prim[l++] = 3; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
	case 128:
		nradices = 7;
		radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
	case 144:
		nradices = 6;
		radix_prim[l++] = 3; radix_prim[l++] = 3; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
	case 160:
		nradices = 6;
		radix_prim[l++] = 5; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
	case 176:
		nradices = 5;
		radix_prim[l++] =11; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
	case 192:
		nradices = 7;
		radix_prim[l++] = 3; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
	case 208:
		nradices = 5;
		radix_prim[l++] =13; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
	case 224:
		nradices = 6;
		radix_prim[l++] = 7; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
	case 240:
		nradices = 6;
		radix_prim[l++] = 5; radix_prim[l++] = 3; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
	case 256:
		nradices = 8;
		radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
	case 288:
		nradices = 7;
		radix_prim[l++] = 3; radix_prim[l++] = 3; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
	case 320:
		nradices = 7;
		radix_prim[l++] = 5; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
	case 352:
		nradices = 6;
		radix_prim[l++] =11; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
	case 384:
		nradices = 8;
		radix_prim[l++] = 3; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
	case 416:
		nradices = 6;
		radix_prim[l++] =13; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
	case 448:
		nradices = 7;
		radix_prim[l++] = 7; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
	case 480:
		nradices = 7;
		radix_prim[l++] = 5; radix_prim[l++] = 3; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
	case 512:
		nradices = 9;
		radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
	case 576:
		nradices = 8;
		radix_prim[l++] = 3; radix_prim[l++] = 3; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
	case 640:
		nradices = 8;
		radix_prim[l++] = 5; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
	case 704:
		nradices = 7;
		radix_prim[l++] =11; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
	case 768:
		nradices = 9;
		radix_prim[l++] = 3; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
	case 832:
		nradices = 7;
		radix_prim[l++] =13; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
	case 896:
		nradices = 8;
		radix_prim[l++] = 7; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
	case 960:
		nradices = 8;
		radix_prim[l++] = 5; radix_prim[l++] = 3; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
	case 992:
		nradices = 6;
		radix_prim[l++] =31; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
	case 1008:
		nradices = 7;
		radix_prim[l++] = 7; radix_prim[l++] = 3; radix_prim[l++] = 3; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
	case 1024:
		nradices = 10;
		radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
	// Next 2 are the leading radices R0 we intend to use for the 2 independent-FFT-length runs of F33 (@overall input-vector lengths of R0*2^17 doubles):
	case 4032:
		nradices = 9;
		radix_prim[l++] = 7; radix_prim[l++] = 3; radix_prim[l++] = 3; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
	case 4096:
		nradices = 12;
		radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
	default :
		printf("FATAL: radix %d not available. Halting...\n",RADIX); exit(EXIT_FAILURE);
	}

	/* Init dit_scramble array: This is the input permutation applied to DIT inputs, in top of the bit-reversal reordering: */
	k = 0;	/* Current perm index */
	l = 0;	/* Current perm value */
	if(pow2 > 1) {	// Composite = odd*pow2, assume odd-subradix DFTs get done first:
		for(i = 0; i < pow2; i++)
		{
			for(j = 0; j < podd; ++j)
			{
				dit_scramble[k++] = l;
				l = (l - pow2)%RADIX;	if(l < 0) { l += RADIX; }
			}
			l = (l - podd)%RADIX;	if(l < 0) { l += RADIX; }
		}
	} else if(nradices > 1 && radix_prim[0] != radix_prim[1]) {	// Composite = odd1*odd2, assume radix_prim[0] contains odd-subradix of DFTs which get done first:
		pow2 = RADIX/radix_prim[0];	// Not a power of 2, but acts same way as decrement
		podd = RADIX/pow2;
		for(i = 0; i < pow2; i++)
		{
			for(j = 0; j < podd; ++j)
			{
				dit_scramble[k++] = l;
				l = (l - pow2)%RADIX;	if(l < 0) { l += RADIX; }
			}
			l = (l - podd)%RADIX;	if(l < 0) { l += RADIX; }
		}
	} else {	// Odd-prime or odd-prime-power radix
		ASSERT((nradices == 1 && (radix_prim[0]&1))
				|| (nradices > 1 && (radix_prim[0] == radix_prim[1])), "Unexpected radix-decomposition!");
	}

/*...Allocate and initialize an index array containing (RADIX) indices...	*/

	for(i=0; i < RADIX; i++)
	{
		index[i]=i;
	}

/*...then bit-reverse INDEX with respect to the [nradices] primitive radices.
		 The order of radices sent to bit_reverse_int is the reverse of that in which these radices are processed
		 in the forward (decimation in frequency) FFT. This is moot for a power-of-2 FFT (or any FFT whose length
		 is a prime power), but necessary for general vector lengths which are a product of 2 or more distinct primes.

		 If the current (Ith) radix is composite with distinct prime factors (e.g. 15 = 3*5), we must specify these
		 factors here in the opposite order from that which is used in the actual FFT-pass routine. For example,
		 if the radix-15 pass implementation does 5 radix-3 DFTs, followed by 3 radix-5 DFTs, then we send (3,5)
		 as the corresponding reverse-ordered prime radices to the bit-reversal routine, not (5,3).	*/

	bit_reverse_int(&index[0],RADIX,nradices,&radix_prim[nradices-1],-1,(int *)0x0);
/*
	printf("bit-reversal index array = [");
	for(i=0; i < RADIX; i++)
	{
		printf(" %d",index[i]);
	}
	printf("]\n\n");
*/
	/* Input-index-permute diagnostics ... the final 'Combined' permutation is that needed for DIT inputs: */
#if TTYPE == 0

	printf("DIF/DIT input-scramble array = [");
	for(i=0; i < RADIX; i++)
	{
		printf("%3x,",dit_scramble[i]);
	}
	printf("]\n\n");
/*
	printf("Bit-reversal array = [");
	for(i=0; i < RADIX; i++)
	{
		j = index[i];
		printf("%3x,",j);
	}
	printf("]\n\n");

	printf("DIT input-scramble + bit-reversal array = [");
	for(i=0; i < RADIX; i++)
	{
		j = dit_scramble[index[i]];
		printf("%3x,",j);
	}
	printf("]\n\n");
*/
	/* Now find the location of each j-valued element above in the bit-reversal vector ... the *location* in that vector equals the final scrambled-DIT-input index: */
	printf("Combined DIT input-scramble array = [");
	for(i=0; i < RADIX; i++)
	{
		j = dit_scramble[index[i]];
		for(k = 0; k < RADIX; ++k)	// Really need to do this via sort and 1-loop rather than O(n^2) nested loop pair...
		{
			if(j == index[k])
			{
				printf("%3x,",k);
				break;
			}
		}
	}
	printf("]\n\n");
	exit(0);
#endif

#if TTYPE == 1 || TTYPE == 3

	#if RADIX == 2
		radix2_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 3
		radix3_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 4
		radix4_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 5
		radix5_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 6
		radix6_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 7
		radix7_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 8
		radix8_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 9
		radix9_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 10
		radix10_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 11
		radix11_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 12
		radix12_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 13
		radix13_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 14
		radix14_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 15
		radix15_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 16
	  #ifdef USE_FGT61
		radix16_dif_pass1 (a,amod,rmul*RADIX);
	  #else
		radix16_dif_pass1 (a,     rmul*RADIX);
	  #endif
	#elif RADIX == 18
		radix18_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 20
		radix20_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 22
		radix22_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 24
		radix24_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 26
		radix26_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 28
		radix28_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 30
		radix30_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 31
		radix31_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 32
		radix32_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 36
		radix36_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 40
		radix40_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 44
		radix44_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 48
		radix48_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 52
		radix52_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 56
		radix56_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 60
		radix60_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 63
		radix63_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 64
		radix64_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 72
		radix72_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 80
		radix80_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 88
		radix88_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 96
		radix96_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 104
		radix104_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 112
		radix112_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 120
		radix120_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 128
		radix128_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 144
		radix144_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 160
		radix160_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 176
		radix176_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 192
		radix192_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 208
		radix208_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 224
		radix224_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 240
		radix240_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 256
		radix256_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 288
		radix288_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 320
		radix320_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 352
		radix352_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 384
		radix384_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 416
		radix416_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 448
		radix448_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 480
		radix480_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 512
		radix512_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 576
		radix576_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 640
		radix640_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 704
		radix704_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 768
		radix768_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 832
		radix832_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 896
		radix896_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 960
		radix960_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 992
		radix992_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 1008
		radix1008_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 1024
		radix1024_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 4032
		radix4032_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 4096
		radix4096_dif_pass1 (a,rmul*RADIX);
	#else
		#error This DIF radix not yet implemented!
	#endif

#endif

#if TTYPE == 1

	// To deduce the required output-idx ordering, sort Actual-outputs data by left col-of-real-parts, move
	// resulting [re,im,i] block to file, repeat procedure for Expected-outputs, then compare the two files.
	// If they differ only in their respective rcol indices, those two cols give the required index mapping.
	// Apr 2014:
	// Fully automated the above procedure - matching problem reduced to an exercise in predicated sorting.
	nerr = idiff = 0;
	avgerr = maxerr = 0.0;
	// Save copy of the 2 resp. datasets in the sort-arrays:
	for(i = 0; i < RADIX ; i++)
	{
		j = index[i];
	#ifdef USE_SSE2
		ASSERT(a[2*i*RE_IM_STRIDE] == a[2*i*RE_IM_STRIDE+1] && a[(2*i+1)*RE_IM_STRIDE] == a[(2*i+1)*RE_IM_STRIDE+1], "1/2 components of SSE2-pack mismatch!");
	#endif
		j1 =  2*i   *RE_IM_STRIDE;	// Real part
		j2 = (2*i+1)*RE_IM_STRIDE;	// Imag part
		j1 += ( (j1 >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
		j2 += ( (j2 >> DAT_BITS) << PAD_BITS );

		sort_arr1[i].idx = i;			sort_arr2[i].idx = i;
		sort_arr1[i].cdat.re = a[j1];	sort_arr2[i].cdat.re = b[2*j];
		sort_arr1[i].cdat.im = a[j2];	sort_arr2[i].cdat.im = b[2*j+1];
	//	printf("I = %3u: DIF-ref: %20.15f  %20.15f,  Actual: %20.15f  %20.15f\n",i, b[2*j],b[2*j+1], a[j1],a[j2]);
		// We only deploy FGT-based DFTs once the floating version has been tested, so no point doing the sorting
		// here, just compare using the (presumably correct) output-index permutations derived for the float code:
	#ifdef USE_FGT61
		printf("I = %3u: DIF-ref: %20llu  %20llu,  FGT: %20llu  %20llu",i, bmod[2*j],bmod[2*j+1], amod[j1],amod[j2]);
		if(bmod[2*j] != amod[j1] || bmod[2*j+1] != amod[j2]) {
			if(bmod[2*j] != qreduce_full(amod[j1]) || bmod[2*j+1] != qreduce_full(amod[j2])) {
				printf("\tDiff = %20lld  %20lld\n",bmod[2*j]-amod[j1], bmod[2*j+1]-amod[j2]);
			} else {
				printf("\tMatch (mod q)\n");
			}
		} else {
			printf("\tMatch\n");
		}
	#endif
	}
	// Sort the 2 structified data arrays by complex value and then compare both complex data and their resp. indices...
	qsort((void *)sort_arr1, RADIX, sz_idx_cmplx, cmp_by_val);
	qsort((void *)sort_arr2, RADIX, sz_idx_cmplx, cmp_by_val);
	for(i = 0; i < RADIX ; i++)
	{
		err_r = fabs(sort_arr1[i].cdat.re - sort_arr2[i].cdat.re);
		err_i = fabs(sort_arr1[i].cdat.im - sort_arr2[i].cdat.im);
		abserr = CABS(err_r, err_i);
		avgerr += abserr;
		maxerr = MAX(maxerr, abserr);
	/*
		printf("%15.10f  %15.10f  %4d  %15.10f  %15.10f  %4d, ERR= %15.10e\n",
			sort_arr1[i].cdat.re,sort_arr1[i].cdat.im,sort_arr1[i].idx,
			sort_arr2[i].cdat.re,sort_arr2[i].cdat.im,sort_arr2[i].idx,
			abserr
		);
	*/
		// DIF and DIT get a more-lenient error level than the combined:
		if(abserr >= err_theshold*sqrt_radix) {
			++nerr;
		} else {
			// If complex data match, compare resp. list indices:
			idiff += (sort_arr1[i].idx != sort_arr2[i].idx);
		}
	}

	// If sorted lists match, extract any needed index-perm:
	if(!nerr && !idiff) {
		// All is well
	} else if(!nerr && idiff) {
		// Data match but need to apply an output iperm in the DFT:
		for(i = 0; i < RADIX ; i++)
		{
			ipair[i].a = sort_arr1[i].idx;		ipair[i].b = sort_arr2[i].idx;
		}
		qsort((void *)ipair, RADIX, sz_int_pair, int_pair_cmp);
		printf("Required output index permutation = [");
		for(i = 0; i < RADIX ; i++)
		{
			printf("%3x,",ipair[i].b);
		}
		printf("]\n");
	} else {
		// There are data mismatches - print the original output lists for the user to inspect prior to code debug:
		printf("\n");
		printf("         Actual outputs:             i          Expected outputs:            i BR oidx:\n");
		printf(" -------------------------------   ---   -------------------------------   --- -------\n");
		for(i = 0; i < RADIX ; i++)
		{
			j = index[i];
			j1 =  2*i   *RE_IM_STRIDE;
			j2 = (2*i+1)*RE_IM_STRIDE;
			j1 += ( (j1 >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
			j2 += ( (j2 >> DAT_BITS) << PAD_BITS );
	
			err_r = fabs(a[j1]-b[2*j]);  err_i = fabs(a[j2]-b[2*j+1]);
			abserr = CABS(err_r, err_i);

			// Unlike the combined fwd/inv test (TTYPE = 3), for solo DIF and DIT if even 1 mismatch we print all
			// compares (even if a match) because may need the full list for debug purposes:
			if(abserr < err_theshold*sqrt_radix) {
				printf("%15.10f  %15.10f  %4d  %15.10f  %15.10f  %4d  (%4d)\n",a[j1],a[j2],i, b[2*j],b[2*j+1],i,j);
			} else {
				printf("%15.10f  %15.10f  %4d  %15.10f  %15.10f  %4d  (%4d), ERR= %15.10e\n",a[j1],a[j2],i, b[2*j],b[2*j+1],i,j, abserr);
				// Early-abort if DC components fail to match:
				if(i == 0) {
					printf("DC components fail to match ... aborting.");
					exit(0);
				}
			}
		}
	}
	avgerr *= iradix;
	printf("test_fft_radix: %d Mismatches detected in DIF DFT; maxerr = %15.10e, avgerr = %15.10e\n", nerr, maxerr, avgerr);
	printf("\n");
	ASSERT(nerr == 0, "test_fft_radix: Mismatches detected in DIF transform!");

#endif

#if TTYPE == 2	// DIT-only, as opposed to combined DIF+DIT [TTYPE = 3]

	// Bit-reverse the inputs to the transform...
  #ifdef USE_SSE2

	t0 = t1 = t2 = t3 = 0.;
	for(i = 0; i < RADIX ; i++)
	{
		j = index[i];
		arrtmp[ 2*i   *RE_IM_STRIDE  ] = a[2*j*RE_IM_STRIDE  ];
		arrtmp[ 2*i   *RE_IM_STRIDE+1] = a[2*j*RE_IM_STRIDE+1];
		arrtmp[(2*i+1)*RE_IM_STRIDE  ] = -a[(2*j+1)*RE_IM_STRIDE  ];
		arrtmp[(2*i+1)*RE_IM_STRIDE+1] = -a[(2*j+1)*RE_IM_STRIDE+1];
		t0 += arrtmp[ 2*i   *RE_IM_STRIDE  ];
		t1 += arrtmp[ 2*i   *RE_IM_STRIDE+1];
		t2 += arrtmp[(2*i+1)*RE_IM_STRIDE  ];
		t3 += arrtmp[(2*i+1)*RE_IM_STRIDE+1];
/*printf("J = [%3d]: add %6d, %6d\n",j,(int)a[2*j  ],(int)a[2*j+1]);*/
	}
/*printf("sum[Re,Im] = %15.5f  %15.5f\n",t0,t2);*/
	ASSERT(t0==t1 && t2==t3, "!");
	for(i = 0; i < rmul*RADIX ; i+=2)
	{
		a[i  ] = arrtmp[i  ];
		a[i+1] = arrtmp[i+1];
		ASSERT(a[i  ] == a[i+1], "!");
	}

  #else

	for(i = 0; i < RADIX ; i++)
	{
		j = index[i];
		j1 =  2*j   ;
		j2 = (2*j+1);
		j1 += ( (j1 >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
		j2 += ( (j2 >> DAT_BITS) << PAD_BITS );
		arrtmp[2*i  ] = a[j1];
		arrtmp[2*i+1] =-a[j2];
	}
	// Now copy the DIT-index-scrambled data back to the A-array:
	for(i = 0; i < rmul*RADIX ; i += 2)
	{
		j1 = i + ( (i >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
		j2 = j1+1;
		a[j1] = arrtmp[i  ];
		a[j2] = arrtmp[i+1];
	#ifdef USE_FGT61
		// Since we negated Im-part above, must analogize to q - (pre-negation)a[j2] here:
		amod[j1] =             a[j1];	// Here, the cast-to-uint64 is implied by the assignment ...
		amod[j2] = q + (uint64)a[j2];	// ...but here need explicit cast to ensure integer addition.
		printf("DIT-in[%2u]: float = [%10.5f,%10.5f]; int = [ %llu, q - %llu]\n",i, a[j1],a[j2] ,amod[j1],q - amod[j2]);
	#endif
	}

  #endif

#endif

#if TTYPE > 1

	#if RADIX == 2
		radix2_dit_pass1 (a,rmul*RADIX);
	#elif RADIX == 3
		radix3_dit_pass1 (a,rmul*RADIX);
	#elif RADIX == 4
		radix4_dit_pass1 (a,rmul*RADIX);
	#elif RADIX == 5
		radix5_dit_pass1 (a,rmul*RADIX);
	#elif RADIX == 6
		radix6_dit_pass1 (a,rmul*RADIX);
	#elif RADIX == 7
		radix7_dit_pass1 (a,rmul*RADIX);
	#elif RADIX == 8
		radix8_dit_pass1 (a,rmul*RADIX);
	#elif RADIX == 9
		radix9_dit_pass1 (a,rmul*RADIX);
	#elif RADIX == 10
		radix10_dit_pass1 (a,rmul*RADIX);
	#elif RADIX == 11
		radix11_dit_pass1 (a,rmul*RADIX);
	#elif RADIX == 12
		radix12_dit_pass1 (a,rmul*RADIX);
	#elif RADIX == 13
		radix13_dit_pass1 (a,rmul*RADIX);
	#elif RADIX == 14
		radix14_dit_pass1 (a,rmul*RADIX);
	#elif RADIX == 15
		radix15_dit_pass1 (a,rmul*RADIX);
	#elif RADIX == 16
	  #ifdef USE_FGT61
		radix16_dit_pass1 (a,amod,rmul*RADIX);
	  #else
		radix16_dit_pass1 (a,     rmul*RADIX);
	  #endif
	#elif RADIX == 18
		radix18_dit_pass1 (a,rmul*RADIX);
	#elif RADIX == 20
		radix20_dit_pass1 (a,rmul*RADIX);
	#elif RADIX == 22
		radix22_dit_pass1 (a,rmul*RADIX);
	#elif RADIX == 24
		radix24_dit_pass1 (a,rmul*RADIX);
	#elif RADIX == 26
		radix26_dit_pass1 (a,rmul*RADIX);
	#elif RADIX == 28
		radix28_dit_pass1 (a,rmul*RADIX);
	#elif RADIX == 30
		radix30_dit_pass1 (a,rmul*RADIX);
	#elif RADIX == 31
		radix31_dit_pass1 (a,rmul*RADIX);
	#elif RADIX == 32
		radix32_dit_pass1 (a,rmul*RADIX);
	#elif RADIX == 36
		radix36_dit_pass1 (a,rmul*RADIX);
	#elif RADIX == 40
		radix40_dit_pass1 (a,rmul*RADIX);
	#elif RADIX == 44
		radix44_dit_pass1 (a,rmul*RADIX);
	#elif RADIX == 48
		radix48_dit_pass1 (a,rmul*RADIX);
	#elif RADIX == 52
		radix52_dit_pass1 (a,rmul*RADIX);
	#elif RADIX == 56
		radix56_dit_pass1 (a,rmul*RADIX);
	#elif RADIX == 60
		radix60_dit_pass1 (a,rmul*RADIX);
	#elif RADIX == 63
		radix63_dit_pass1 (a,rmul*RADIX);
	#elif RADIX == 64
		radix64_dit_pass1 (a,rmul*RADIX);
	#elif RADIX == 72
		radix72_dit_pass1 (a,rmul*RADIX);
	#elif RADIX == 80
		radix80_dit_pass1 (a,rmul*RADIX);
	#elif RADIX == 88
		radix88_dit_pass1 (a,rmul*RADIX);
	#elif RADIX == 96
		radix96_dit_pass1 (a,rmul*RADIX);
	#elif RADIX ==104
		radix104_dit_pass1(a,rmul*RADIX);
	#elif RADIX ==112
		radix112_dit_pass1(a,rmul*RADIX);
	#elif RADIX ==120
		radix120_dit_pass1(a,rmul*RADIX);
	#elif RADIX ==128
		radix128_dit_pass1(a,rmul*RADIX);
	#elif RADIX == 144
		radix144_dit_pass1 (a,rmul*RADIX);
	#elif RADIX == 160
		radix160_dit_pass1 (a,rmul*RADIX);
	#elif RADIX == 176
		radix176_dit_pass1 (a,rmul*RADIX);
	#elif RADIX == 192
		radix192_dit_pass1 (a,rmul*RADIX);
	#elif RADIX == 208
		radix208_dit_pass1(a,rmul*RADIX);
	#elif RADIX == 224
		radix224_dit_pass1(a,rmul*RADIX);
	#elif RADIX == 240
		radix240_dit_pass1(a,rmul*RADIX);
	#elif RADIX == 256
		radix256_dit_pass1(a,rmul*RADIX);
	#elif RADIX == 288
		radix288_dit_pass1 (a,rmul*RADIX);
	#elif RADIX == 320
		radix320_dit_pass1 (a,rmul*RADIX);
	#elif RADIX == 352
		radix352_dit_pass1 (a,rmul*RADIX);
	#elif RADIX == 384
		radix384_dit_pass1 (a,rmul*RADIX);
	#elif RADIX == 416
		radix416_dit_pass1(a,rmul*RADIX);
	#elif RADIX == 448
		radix448_dit_pass1(a,rmul*RADIX);
	#elif RADIX == 480
		radix480_dit_pass1(a,rmul*RADIX);
	#elif RADIX == 512
		radix512_dit_pass1(a,rmul*RADIX);
	#elif RADIX == 576
		radix576_dit_pass1 (a,rmul*RADIX);
	#elif RADIX == 640
		radix640_dit_pass1 (a,rmul*RADIX);
	#elif RADIX == 704
		radix704_dit_pass1 (a,rmul*RADIX);
	#elif RADIX == 768
		radix768_dit_pass1 (a,rmul*RADIX);
	#elif RADIX == 832
		radix832_dit_pass1(a,rmul*RADIX);
	#elif RADIX == 896
		radix896_dit_pass1(a,rmul*RADIX);
	#elif RADIX == 960
		radix960_dit_pass1(a,rmul*RADIX);
	#elif RADIX == 992
		radix992_dit_pass1(a,rmul*RADIX);
	#elif RADIX == 1008
		radix1008_dit_pass1 (a,rmul*RADIX);
	#elif RADIX == 1024
		radix1024_dit_pass1(a,rmul*RADIX);
	#elif RADIX == 4032
		radix4032_dit_pass1 (a,rmul*RADIX);
	#elif RADIX == 4096
		radix4096_dit_pass1(a,rmul*RADIX);
	#else
		#error This DIT radix not yet implemented!
	#endif

#endif

#if TTYPE == 2

	nerr = idiff = 0;
	avgerr = maxerr = 0.0;
	// Save copy of the 2 resp. datasets in the sort-arrays:
	for(i = 0; i < RADIX ; i++)
	{
		j1 =  2*i   *RE_IM_STRIDE;	// Real part
		j2 = (2*i+1)*RE_IM_STRIDE;	// Imag part
		j1 += ( (j1 >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
		j2 += ( (j2 >> DAT_BITS) << PAD_BITS );

		sort_arr1[i].idx = i;			sort_arr2[i].idx = i;
		sort_arr1[i].cdat.re = a[j1];	sort_arr2[i].cdat.re = b[2*i];
		sort_arr1[i].cdat.im = a[j2];	sort_arr2[i].cdat.im =-b[2*i+1];	// Flip sign on Im-part of ref-output
		// We only deploy FGT-based DFTs once the floating version has been tested, so no point doing the sorting
		// here, just compare using the (presumably correct) output-index permutations derived for the float code:
	#ifdef USE_FGT61
		// Flip sign on Im-part of ref-outputs:
		printf("I = %3u: DIT-ref: %20llu  %20llu,  FGT: %20llu  %20llu",i, bmod[2*i],q-bmod[2*i+1], amod[j1],amod[j2]);
		if(bmod[2*i] != amod[j1] || q-bmod[2*i+1] != amod[j2]) {
			if(bmod[2*i] != qreduce_full(amod[j1]) || q-bmod[2*i+1] != qreduce_full(amod[j2])) {
				printf("\tDiff = %20lld  %20lld\n",bmod[2*i]-amod[j1], (q-bmod[2*i+1])-amod[j2]);
			} else {
				printf("\tMatch (mod q)\n");
			}
		} else {
			printf("\tMatch\n");
		}
	#endif
	}
	// Sort the 2 structified data arrays by complex value and then compare both complex data and their resp. indices...
	qsort((void *)sort_arr1, RADIX, sz_idx_cmplx, cmp_by_val);
	qsort((void *)sort_arr2, RADIX, sz_idx_cmplx, cmp_by_val);
	for(i = 0; i < RADIX ; i++)
	{
		err_r = fabs(sort_arr1[i].cdat.re - sort_arr2[i].cdat.re);
		err_i = fabs(sort_arr1[i].cdat.im - sort_arr2[i].cdat.im);
		abserr = CABS(err_r, err_i);
		avgerr += abserr;
		maxerr = MAX(maxerr, abserr);
		// DIF and DIT get a more-lenient error level than the combined:
		if(abserr >= err_theshold*sqrt_radix) {
			++nerr;
		} else {
			// If complex data match, compare resp. list indices:
			idiff += (sort_arr1[i].idx != sort_arr2[i].idx);
		}
	}
	// If sorted lists match, extract any needed index-perm:
	if(!nerr && !idiff) {
		// All is well
	} else if(!nerr && idiff) {
		// Data match but need to apply an output iperm in the DFT:
		for(i = 0; i < RADIX ; i++)
		{
			ipair[i].a = sort_arr1[i].idx;		ipair[i].b = sort_arr2[i].idx;
		}
		qsort((void *)ipair, RADIX, sz_int_pair, int_pair_cmp);
		printf("Required output index permutation = [");
		for(i = 0; i < RADIX ; i++)
		{
			printf("%3x,",ipair[i].b);
		}
		printf("]\n");
	} else {
		// There are data mismatches - print the original output lists for the user to inspect prior to code debug:
		printf("\n");
		printf("         Actual outputs:             i          Expected outputs:            i BR oidx:\n");
		printf(" -------------------------------   ---   -------------------------------   --- -------\n");
		for(i = 0; i < RADIX ; i++)
		{
			j1 =  2*i   *RE_IM_STRIDE;
			j2 = (2*i+1)*RE_IM_STRIDE;
			j1 += ( (j1 >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
			j2 += ( (j2 >> DAT_BITS) << PAD_BITS );

			// Flip sign on Im-part of ref-output
			err_r = fabs(a[j1]-b[2*i]);  err_i = fabs(a[j2]+b[2*i+1]);
			abserr = CABS(err_r, err_i);

			// Unlike the combined fwd/inv test (TTYPE = 3), for solo DIF and DIT if even 1 mismatch we print all
			// compares (even if a match) because may need the full list for debug purposes:
			if(abserr < err_theshold*sqrt_radix) {
				printf("%15.10f  %15.10f  %4d  %15.10f  %15.10f  %4d\n",a[j1],a[j2],i, b[2*i],-b[2*i+1],i);
			} else {
				printf("%15.10f  %15.10f  %4d  %15.10f  %15.10f  %4d, ERR= %15.10e\n",a[j1],a[j2],i, b[2*i],-b[2*i+1],i, abserr);
				// Early-abort if DC components fail to match:
				if(i == 0) {
					printf("DC components fail to match ... aborting.");
					exit(0);
				}
			}
		}
	}
	avgerr *= iradix;
	printf("test_fft_radix: %d Mismatches detected in DIT DFT; maxerr = %15.10e, avgerr = %15.10e\n", nerr, maxerr, avgerr);
	printf("\n");
	ASSERT(nerr == 0, "test_fft_radix: Mismatches detected in DIT transform!");

#endif

#if TTYPE == 3

	nerr = 0;
	avgerr = maxerr = 0.0;
	for(i = 0; i < RADIX ; i++)
	{
	#ifdef USE_SSE2
		ASSERT(a[2*i*RE_IM_STRIDE] == a[2*i*RE_IM_STRIDE+1] && a[(2*i+1)*RE_IM_STRIDE] == a[(2*i+1)*RE_IM_STRIDE+1], "1/2 components of SSE2-pack mismatch!");
	#endif
		j1 =  2*i   *RE_IM_STRIDE;
		j2 = (2*i+1)*RE_IM_STRIDE;
		j1 += ( (j1 >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
		j2 += ( (j2 >> DAT_BITS) << PAD_BITS );
		a[j1] *= iradix;	a[j2] *= iradix;
		err_r = fabs(a[j1]-arrtmp[2*i]);  err_i = fabs(a[j2]-arrtmp[2*i+1]);
		abserr = CABS(err_r, err_i);
		avgerr += abserr;
		maxerr = MAX(maxerr, abserr);
		if(abserr >= err_theshold) {
			++nerr;
			printf("%4d  %25.15f  %25.15f, ERR= %15.10e\n",i,a[j1], a[j2], CABS(err_r, err_i));
		}
	#ifdef USE_FGT61
		printf("I = %3u: DIF+DIT ref: [%lld,%lld],  FGT: [%20llu,%20llu]",i, (uint64)arrtmp[2*i],(uint64)arrtmp[2*i+1], amod[j1]/RADIX,amod[j2]/RADIX);
		if((uint64)arrtmp[2*i] != amod[j1]/RADIX || (uint64)arrtmp[2*i+1] != amod[j2]/RADIX) {
			if((uint64)arrtmp[2*i] != qreduce_full(amod[j1])/RADIX || (uint64)arrtmp[2*i+1] != qreduce_full(amod[j2])/RADIX) {
				printf("\tMismatch! mod-outputs (mod RADIX) = [%20llu,%20llu]\n",amod[j1]%RADIX, amod[j2]%RADIX);
			} else {
				printf("\tMatch (mod q)\n");
			}
		} else {
			printf("\tMatch\n");
		}
	#endif
	}
	avgerr *= iradix;
	printf("test_fft_radix: %d Mismatches detected in DIF/DIT combo; maxerr = %15.10e, avgerr = %15.10e\n", nerr, maxerr, avgerr);
	printf("\n");
	ASSERT(nerr == 0, "test_fft_radix: Mismatches detected in DIF/DIT combo!");

#endif
	printf("");
}

void matmul_double(double **mat, double vec_in[], double vec_out[], int nrow, int ncol)
{
	int i,j;
	double tmp;

	/*	Loop over [nrow] rows of matrix; inner product of each row with length-[ncol] vector is stored in a second vector of length [nrow]. */
	for(i = 0; i < nrow; i++)
	{
		tmp = 0.0;
		for(j = 0; j < ncol; j++)
		{
			tmp += mat[i][j]*vec_in[j];
		}
		vec_out[i] = tmp;
	}

	return;
}

void matmul_complex(struct complex **mat, struct complex vec_in[], struct complex vec_out[], int nrow, int ncol)
{
	int i,j;
	struct complex tmp;

	/*	Loop over [nrow] rows of matrix; inner product of each row with length-[ncol] vector is stored in a second vector of length [nrow]. */
	for(i = 0; i < nrow; i++)
	{
		vec_out[i].re = 0.0;
		vec_out[i].im = 0.0;
		for(j = 0; j < ncol; j++)
		{
			tmp = cmul(&mat[i][j], &vec_in[j]);
			vec_out[i].re += tmp.re;
			vec_out[i].im += tmp.im;
		/*
			printf("mat[%4d][%4d] = %15.10f  %15.10f\n",i,j, mat[i][j].re, mat[i][j].im);
			printf("cmul * vec_in = %15.10f  %15.10f\n",     vec_in[j].re, vec_in[j].im);
		*/
		}
	/*
		printf("--------------------------------\n");
		printf("vec_out[%4d] = %15.10f  %15.10f\n",i, vec_out[i].re,vec_out[i].im);
	*/
	}

	return;
}

#ifdef USE_FGT61
// Modular analog of matmul_complex, over FGT(M61^2):
void matmul_fgtmod(uint128 **mat, uint128 vec_in[], uint128 vec_out[], int nrow, int ncol)
{
	const uint64 q  = 0x1FFFFFFFFFFFFFFFull;	// q = 2^61 - 1
	int i,j;
	uint64 rm,im;
	for(i = 0; i < nrow; i++)
	{
	//if(!i) printf("matmul_fgtmod: DC-component inner product terms:\n");
		vec_out[i].d0 = vec_out[i].d1 = 0ull;
		for(j = 0; j < ncol; j++)
		{
			cmul_modq(mat[i][j].d0,mat[i][j].d1, vec_in[j].d0,vec_in[j].d1, &rm,&im);
			// CMUL_MODQ outputs in 0,4b - must feed to qreduce() prior to accumulating:
			rm = qreduce(rm);	im = qreduce(im);
	//	if(!i) printf("\t[%2d] = [%llu,%llu] * [%llu,%llu] = [%llu,%llu]\n",j, mat[i][j].d0,mat[i][j].d1, vec_in[j].d0,vec_in[j].d1, rm,im);
			rm += vec_out[i].d0;
			im += vec_out[i].d1;
			// Normalize to ensure accumulated sum in [0,q-1]:
			vec_out[i].d0 = rm - (-(uint64)(rm >= q)) & q;
			vec_out[i].d1 = im - (-(uint64)(im >= q)) & q;
		}
	}

	return;
}
#endif
