/*******************************************************************************
*                                                                              *
*   (C) 1997-2012 by Ernst W. Mayer.                                           *
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

#include "gcd_lehmer.h"

/*
   Set == 1 to enable basic per-pass diagnostic printing.
   Set >= 2 to enable detailed per-pass diagnostic printing.
            WARNING: level-2 diagnostics not recommended for large vectors!

If the FFT-mul enabled, need all the stuff for factor.c + int-gcd,
including factor.c built with -DINCLUDE_PM1, plus the following FFT stuff:

g42 -c -Wall -g3 -ggdb -DFACTOR_STANDALONE -DTRYQ=4 -DNWORD -DINCLUDE_PM1 factor.c
g42 -c -Wall -O3 br.c get_fft*c *pairFFT*.c
g42 -c -Wall -O3 radix32_dif_dit_pass.c
g42 -c -Wall -O3 radix16_dif_dit_pass.c
g42 -c -Wall -O3 radix8_dif_dit_pass.c
g42 -c -Wall -O3 -DGCD_STANDALONE radix32_*cy*.c
g42 -c -Wall -O3 -DGCD_STANDALONE radix16_*cy*.c
g42 -c -Wall -O3 -DGCD_STANDALONE radix8_*cy*.c

*/
#define GCD_DEBUG	0

char string0[STR_MAX_LEN];
int gcd_debug = 0;	/* Set GCD_DEBUG >= 2 and this != 0 to selectively enable debug
				printing for some particular self-test or phase of the computation. */

#if GCD_DEBUG >= 1
	char string1[STR_MAX_LEN];
	char str_10k[10*1024];
#endif
#if GCD_DEBUG >= 2
	char string2[STR_MAX_LEN];
#endif

#ifdef GCD_STANDALONE
	/* In standalone, need to locally define some Mlucas.h-declared globals
	which would normally be defined in Mlucas.c. Init these via call to */

	double ISRT2;

	int32 DAT_BITS, PAD_BITS;	/* Array padding parameters */

	uint32 PMIN;	/* minimum #bits allowed for FFT-based mul */
	uint32 PMAX;	/* maximum #bits allowed depends on max. FFT length allowed
					  and will be determined at runtime, via call to given_N_get_maxP(). */
#endif

#define FFTMUL_THRESHOLD_BITS	16384	/* # of bits in a multiword base-2^64 integer
										at which to switch from pure-int grammar-school
										to floating-FFT-based multiply algorithm */

/***********************/

/* Do some basic inits and self-tests: */
void gcd_init(void)
{
	uint32 n;
	int nradices, radix_vec[10];	/* A pair of locals for storing these data
									returned by get_ftt_radices() can be handy... */

	host_init();

	/* Simple self-tester for get_fft_radices(): */
	test_fft_radixtables();

	/*...Set the array padding parameters to disable padding for now: */
	if(DAT_BITS < 0 || PAD_BITS < 0)
	{
	  #if 0
		if(kblocks > 32)
		{
			DAT_BITS = DAT_BITS_DEF;
			PAD_BITS = PAD_BITS_DEF;
		}
		else
		{
			DAT_BITS =31;	/* This causes the padding to go away */
			PAD_BITS = 0;
		}
	  #else
		DAT_BITS =31;	/* This causes the padding to go away */
		PAD_BITS = 0;
	  #endif
	}
	else
	{
		fprintf(stderr,"INFO: Array padding params previously set: DAT_BITS = %d, PAD_BITS = %d\n", DAT_BITS, PAD_BITS);
	}

	/* Set min. exponent (in terms of power of 2) that can be tested: */
	/* Check that the purported min. FFT length is actually supported: */
	ASSERT(HERE, get_fft_radices(MIN_FFT_LENGTH_IN_K, 0, &nradices, radix_vec, 10) == 0,"gcd_lehmer.c: get_fft_radices(MIN_FFT_LENGTH_IN_K, 0) == 0");
	n = (MIN_FFT_LENGTH_IN_K << 10);
	/* Make sure N didn't overflow */
	ASSERT(HERE, (n >> 10) == MIN_FFT_LENGTH_IN_K,"gcd_lehmer.c: (n >> 10) == MIN_FFT_LENGTH_IN_K");
	PMIN = 2*n;	/* 2 bits per input is about the smallest we can test without getting nonzero-carry errors */

	/* Set max. exponent (in terms of power of 2) that can be tested: */
	/* Check that the purported max. FFT length is actually supported: */
	ASSERT(HERE, get_fft_radices(MAX_FFT_LENGTH_IN_K, 0, &nradices, radix_vec, 10) == 0,"gcd_lehmer.c: get_fft_radices(MAX_FFT_LENGTH_IN_K, 0) == 0");
	n = (MAX_FFT_LENGTH_IN_K << 10);
	/* Make sure N didn't overflow */
	ASSERT(HERE, (n >> 10) == MAX_FFT_LENGTH_IN_K,"gcd_lehmer.c: (n >> 10) == MAX_FFT_LENGTH_IN_K");
	PMAX = given_N_get_maxP(n);

	ASSERT(HERE, PMAX > PMIN,"gcd_lehmer.c: PMAX > PMIN");
}

/***********************/

/*
Lehmer-type Euclidean GCD(U, V), where U and V are two multiword integers
stored in base-2^64 positive-digit form; U and V are assumed to both have
at least (ndim) allocated 64-bit elements, with U,V[0] being the least-significant.

RETURNS: The GCD in base-2^64 form in both U and V, and the length of the GCD
		(in terms of # of significant 64-bit words) in the function value.
		If the boolean EGCD is set on entry, an extended GCD is performed and
		pointers to the (absolute values of the) A and B-array multipliers are
		returned in the array pointers Ap, Bp, with corresponding signs indicated
		via the low 2 bits of the returned sign_AB variable.
		These array multipliers are defined such that A*U + B*V = GCD(U, V).

EFFECTS: If both inputs are nonzero, the ORIGINAL INPUT VECTOR CONTENTS ARE DESTROYED
		during the computation, and a copy of the GCD returned in each vector.

		If one or both inputs = 0, neither input is altered, and the function
		return value is = 0 by convention to signal that this is the case.
		In order to match mathematical convention, the user should treat the
		GCD as 0 only if both U and V = 0, and equal to max(U, V)
		(i.e. the nonzero vector among U and V) if only one of the pair = 0.
*/
uint32	mi64_gcd(uint64 u[], uint64 v[], uint32 ndim, uint32 EGCD, uint64 **Ap, uint64 **Bp, uint32 *len_AB, uint32 *sign_AB)
{
#if GCD_DEBUG>=1
	/* stuff needed for known-divisor check: */
	static uint64 *div = 0x0, *rem = 0x0;
#endif
	static uint64 *quo = 0x0;	// Use this in non-debug mode, as well
/* Set max. scalar multiplier bits = 63 if 128-bit integer multiply available; otherwise set = 52. */
#if USE_FLOATING_MULH64
	const int max_bits = 52;
#else
	const int max_bits = 63;
#endif

	int i,istart,j,k,len_diff;
	uint32 j1,j2,k1 = 0,k2 = 1,lz0,lz1,tz0,tz1,min_trailing_zeros=0,nshift,egcd_sign_swap = 0;
	uint64 *uv_ptr[2];
	uint64 *x;
	uint64 iop,lead_u[2],lead_v[2],coeff[2][2],diff[2],range[2][2];
	uint64 abmul[2],cdmul[2],max_abcd,t1,t2,t3,t4,lo,hi;
	const double rbase = pow(2.0, 64), log10_base = 19.26591972249479649367928926;	/* 2^64, log10(2^64) */

	double num,den,frac;
	uint32 pass=0, lenU, lenV;	/* Note that lenT stores the NUMBER OF SIGNIFICANT ELEMENTS of vector T
								(i.e. leading element is in T[lenT-1]), not the # of allocated elements! */
	uint32 tmp32;
	static uint32 ndim_save = 0;/* Stores the current allocation for the multiplier arrays
								A,B,C,D used for extended GCD. Only used if EGCD = TRUE. */
	uint32 len_abcd;
	static uint64 *A = 0x0,*B = 0x0,*C = 0x0,*D = 0x0;	/* Multiplier arrays for extended GCD. For simplicity
								(e.g. so we can use the mi64 functions to manipulate) we use unsigned ints to
								store the magnitudes, but separately track the corresponding sign bits... */
	uint32 sign_abcd;	/* ...in the low 2 bits of this variable (i.e. the 0-bit of sign_abcd
						stores sign(A) and the 1-bit of sign_abcd stores sign(C). Since each of
						B and D must have the opposite sign of A and C, respectively, we need not
						explicitly store the latter two signs. */
	/* Set up for easily-permutable access to the (vector) elements
	of the 2x2 eGCD multiplier matrix, in columnwise fashion: */
	uint64 *ac_ptr[2];
	uint64 *bd_ptr[2];

	ASSERT(HERE, ndim != 0,"ndim == 0");
	ASSERT(HERE,!mi64_iszero(u, ndim),"U-input 0");
	ASSERT(HERE,!mi64_iszero(v, ndim),"V-input 0");

	quo = (uint64 *)calloc(ndim, sizeof(uint64));
#if GCD_DEBUG>=1
	/* known divisor, quotient, remainder vectors needed for mi64_div check: */
	div = (uint64 *)calloc(ndim, sizeof(uint64));
	rem = (uint64 *)calloc(ndim, sizeof(uint64));
	// To check divisibility-by-known-divisor, set divisor above, then insert the following code (with the #if restored) as desired:
	#if 0//GCD_DEBUG>=1
	div[0] =  4235226679561594903ull;
	div[1] =    94004235929829273ull;
	if(gcd_debug) {
		mi64_div(uv_ptr[k2],div,quo,rem,lenU);
		ASSERT(HERE, mi64_iszero(rem, lenU), "divisor check fails!");
	}
	#endif
#endif
#if GCD_DEBUG>=2
	if(1)
	{
		printf("Inputs to GCD are:\n");
		printf("                         U                   V  \n");
		printf("              base=2^64;x=1;u=v=0;\n");
		for(i=0; i<ndim; i++)
		{
			convert_uint64_base10_char(string1, u[i]);
			convert_uint64_base10_char(string2, v[i]);
			printf("i = %8u:  u+=%s*x; v+=%s*x; x*=base;\n",i,string1,string2);
		}
		printf("\n");
	}
#endif

	/* Test for EGCD: */
	if(EGCD)
	{
		if(ndim_save < ndim)
		{
			/* If mul-array allocation not yet performed or insufficient, allocate
			to hold (ndim) uint64s: */
			free((void *)A); A = (uint64 *)calloc(ndim, sizeof(uint64));
			free((void *)B); B = (uint64 *)calloc(ndim, sizeof(uint64));
			free((void *)C); C = (uint64 *)calloc(ndim, sizeof(uint64));
			free((void *)D); D = (uint64 *)calloc(ndim, sizeof(uint64));
			ndim_save = ndim;
		}

		/* Initialize:

			(  A  -B )   ( 1  0 )
			(        ) = (      )
			( -C   D )   ( 0  1 ) .

		We update the vectors forming the elements of this 2x2 multiplier matrix
		such that any stage of the computation, we have

			( u )   (  A  -B ) ( U )
			(   ) = (        )*(   )
			( v )   ( -C   D ) ( V ) , where U and V are the original input vectors (not stored).
		*/
		len_abcd  = 1;
		/* We separately store the magnitudes and signs of the elements of the A/B/C/D matrix -
		the initial magnitudes are those of the 2x2 identity-matrix elements, i.e. 1/0/0/1;
		the initial signs are A/B = +/- and C/D = -/+, which maps to bit pattern 10
		in the low 2 bits of sign_abcd: */
		sign_abcd = 0x2;
		A[0] = D[0] = 1;
		B[0] = C[0] = 0;
	}
	else
	{
		ASSERT(HERE, (A==0x0 && B==0x0 && C==0x0 && D==0x0),"(A==0x0 && B==0x0 && C==0x0 && D==0x0)");
	}


	/* Init the pointers to the U/V arrays: */
	uv_ptr[0] = u;	uv_ptr[1] = v;

	/* Init the pointers to the A/B/C/D arrays: */
	ac_ptr[0] = A;	bd_ptr[0] = B;
	ac_ptr[1] = C;	bd_ptr[1] = D;

FIND_LARGER:

	/* k1 and k2 index into the uv_ptr array - uv_ptr[k1] points to larger of u and v: */
	istart = ndim-1;
	lenU = 0;	/* This allows us to determine whether lenU has already been set on a previous pass in the ensuing loop. */
	for(i=istart; i>=0; i--)
	{
		if((u[i] != 0 || v[i] != 0) && lenU == 0)
			lenU = i+1;
		/* Even if we've already found the leading nonzero element,
		in order to properly set k1 and k2 continue until find a differing element:
		*/
		if(u[i] == v[i])
			continue;
		else
		{
			k1 = (u[i] < v[i]);
			k2 = k1 ^ 0x00000001;
			break;
		}
	}
	/*	If U = V, both inputs identical:
	*/
	if(i == -1)
	{
		printf("INFO: Identical inputs to MI64_GCD.\n");
		lenV = lenU;	k1 = 0; k2 = 1;
		goto GCD_DONE;
	}

	/* Otherwise, parse the smaller array until find leading nonzero element: */
	x = uv_ptr[k2];
	for(i=lenU-1; i>=0; i--)
	{
		if(x[i] != 0)
			break;
	}
	lenV = i + 1;

	/*	If lenV = 0, one of the inputs = 0, which is an error: */
	if(lenV == 0)
	{
		fprintf(stderr, "ERROR: Zero input vector to MI64_GCD!\n");
		ASSERT(HERE, 0,"0");
	}

	/* A small additional bit of preprocessing: count the trailing zeros in each vector and
	right-shift both by the respective amounts, saving the smaller of the 2 shift counts
	for output of GCD result as GCD(right-justified inputs) * 2^(initial right-shift count).
	We need to make both inputs into the main loop odd here, since in the loop there's a
	left-justify-the-smaller-of-the-two-vectors step which would otherwise lead to spurious
	powers of 2 appearing in the final result.
	*/

#if 1
	tz0 = mi64_trailz(uv_ptr[k1], lenU);
	tz1 = mi64_trailz(uv_ptr[k2], lenV);
	min_trailing_zeros = MIN(tz0,tz1);
#else
/**** 8/23/2005: disable for now, since this really complicates things when doing eGCD. ****/
	/* Right-shift U: */
	tz0 = mi64_trailz(uv_ptr[k1], lenU);
	mi64_shrl(uv_ptr[k1], uv_ptr[k1], tz0, lenU);	lenU = mi64_getlen(uv_ptr[k1], lenU);

	/* Right-shift V: */
	tz1 = mi64_trailz(uv_ptr[k2], lenV);
	mi64_shrl(uv_ptr[k2], uv_ptr[k2], tz1, lenV);	lenV = mi64_getlen(uv_ptr[k2], lenV);

	/* 21 Aug 2005:
	How does such a right-shift affect the eGCD multiplier arrays?
	In doing an eGCD, at any stage of the computation we are assumed to have
	(where U and V are the original input vectors, not stored)

		( u )   (  A  -B ) ( U )
		(   ) = (        )*(   )
		( v )   ( -C   D ) ( V ) ,

	So multiplying u or v by some constant is equivalent to multiplying both of A/B or C/D, respectively,
	by the same constant. In the case of a right-shift of u or v, since we can't assume that each of
	A/B or C/D will be divisible by the same power of 2, we instead LEFT-SHIFT the complement vector (v or u)
	by the same number of places, accumulating these u/v shifts in separate counters so we can undo them
	at the conclusion of the computation.
	*/

	if(EGCD)
	{
		if     (tz0 > tz1)
		{
			/* Check left-shift return value to see if len_ancd needs to be updated: */
			******************** start here - need to pad length by abs(tz0-tz1)/64 + 1 prior to shifting *********************************
			mi64_shl(ac_ptr[k2], ac_ptr[k2], tz0-tz1, len_abcd);
			mi64_shl(bd_ptr[k2], bd_ptr[k2], tz0-tz1, len_abcd);	len_abcd = MAX(mi64_getlen(ac_ptr[k2], len_abcd) , mi64_getlen(bd_ptr[k2], len_abcd));
		}
		else if(tz0 < tz1)
		{
			/* len_abcd will be automatically updated during the 2 shifts
			so the final value corresponds to the larger of the 2 result vectors: */
			mi64_shl(ac_ptr[k1], ac_ptr[k1], tz1-tz0, len_abcd);
			mi64_shl(bd_ptr[k1], bd_ptr[k1], tz1-tz0, len_abcd);	len_abcd = MAX(mi64_getlen(ac_ptr[k1], len_abcd) , mi64_getlen(bd_ptr[k1], len_abcd));
		}
	}

	/* If this is our second and final pass this way, make sure both operands are odd, i.e. tz0 = tz1 = 0: */
	if(min_trailing_zeros != 0)
	{
		ASSERT(HERE, MAX(tz0,tz1) == 0,"gcd_lehmer.c: MAX(tz0,tz1) == 0");
	}
	else if(MIN(tz0,tz1) > 0)
	{
		ASSERT(HERE, min_trailing_zeros == 0,"gcd_lehmer.c: min_trailing_zeros == 0");	/* Should only have to do this once. */
		min_trailing_zeros = MIN(tz0,tz1);
		goto FIND_LARGER;
	}

  #if GCD_DEBUG>=2
	if(gcd_debug)
	{
		printf("Right-justified Inputs to GCD are:\n");
		printf("                         U                   V  \n");
		for(i=0; i<lenU; i++)
		{
			convert_uint64_base10_char(string1, u[i]);
			convert_uint64_base10_char(string2, v[i]);
			printf("i = %8u:  %s   %s\n",i,string1,string2);
		}
		printf("\n");
	}
  #endif

#endif

#if GCD_DEBUG>=2
if(gcd_debug)
{
	printf("Indices into U/V = %u, %u\n",k1,k2);
	printf("\n");
}
#endif

/*	Main loop is here.	*/

#if GCD_DEBUG>=1
if(gcd_debug)
	printf("GCD pass = %u, lenU = %u\n",pass,lenU);
#endif

for(;;)	/* MAIN */
{
	pass++;

#if 1
	lz0 = mi64_leadz(uv_ptr[k1], lenU);
	lz1 = mi64_leadz(uv_ptr[k2], lenV) + ((lenU - lenV)<<6);
	ASSERT(HERE, lz0 <= lz1,"gcd_lehmer.c: lz0 > lz1!");
	#if GCD_DEBUG>=1
	if(gcd_debug)
		printf("lz_diff = %u bits\n",lz1-lz0);
	#endif
	if(lz1 - lz0 > 64)
	{
		sprintf(string0,"U (lz0 = %u bits) more than 64 bits larger than V (lz1 = %u bits) in MI64_GCD!\n",lz0,lz1);
		WARN(HERE, string0, "", TRUE);
	}
#endif
/**** 8/23/2005: disable for now, since this really complicates things when doing eGCD. ****/
	/* In order to gurantee a nontrivial 2x2 scalar multiplier matrix,
	make sure that U is no more than max_bits bits greater than V. in
	fact, while we're at it, let's make the difference < max_bits/2: */
#if 1
	/* First, do any needed full-word-sized shifts: */
	len_diff = lenU - lenV;
	ASSERT(HERE, len_diff >= 0,"gcd_lehmer.c: len_diff >= 0");

	if(len_diff > (lenU>>1))	// V less than 1/2 length of U
	{
		// Try division-with-remainder in large-size-disparity cases - this is very slow in terms of timing,
		// but only because the current div implementation is bitwise-quadratic:
		mi64_div(uv_ptr[k1],uv_ptr[k2],lenU,lenU,quo,uv_ptr[k1]);	// To-Do: Use quotient to adjuct the eGCD coeffs to reflect the div.
		ASSERT(HERE, mi64_cmpult(uv_ptr[k1],uv_ptr[k2], lenU), "divisor check fails!");
		j = mi64_getlen(uv_ptr[k1], lenU);
		ASSERT(HERE, (j <= lenV) && mi64_cmpult(uv_ptr[k1],uv_ptr[k2],lenV),"gcd_lehmer.c: lenV == lenU");
		lenU = lenV; lenV = j;
		k1 ^= 0x00000001;
		k2 ^= 0x00000001;
		/* If lenV = 0, we're done. */
		if(lenV == 0)
			goto GCD_DONE;
	}
	else if(len_diff > 0)
	{
		/*	if lenU > lenV and leading digit(v) < leading digit(u),
		left-shift V by len_diff places	*/
		if(uv_ptr[k2][lenV-1] < uv_ptr[k1][lenU-1])
		{
		#if GCD_DEBUG>=2
			if(gcd_debug)
			{
				printf("pass %u , A: left-shifting V-array by %d words\n",pass,len_diff);
			}
		#endif
			mi64_shl(uv_ptr[k2], uv_ptr[k2],  len_diff<<6, lenU);
			lenV = mi64_getlen(uv_ptr[k2], lenU);
			ASSERT(HERE, lenV == lenU,"gcd_lehmer.c: lenV == lenU");
		}
		/*	otherwise left-shift V by (len_diff - 1) places.	*/
		else if(len_diff > 1)
		{
		#if GCD_DEBUG>=2
			if(gcd_debug)
				printf("pass %u , B: left-shifting array %u by %d words\n",pass,k2,len_diff-1);
		#endif
			mi64_shl(uv_ptr[k2], uv_ptr[k2], (len_diff-1)<<6, lenU);
			lenV = mi64_getlen(uv_ptr[k2], lenU);
			ASSERT(HERE, lenV == lenU-1,"gcd_lehmer.c: lenV == lenU-1");
		}

		/*	Now that lenV is no less than lenU - 1, count leading zeros and
		see if any less-than-full-word-sized supplemental shift is warranted:
		*/
		lz0 = mi64_leadz(uv_ptr[k1], lenU);
		lz1 = mi64_leadz(uv_ptr[k2], lenV) + ((lenU - lenV)<<6);
		j = lz1 - lz0;
		ASSERT(HERE, j >= 0,"ERROR: lz0 > lz1 in MI64_GCD!");
		if(j >= (max_bits>>1))
		{
			// If larger vector has any trailing zeros resulting from a previous such shift, right-justify it first:
			i = trailz64(uv_ptr[k1][0]);
			if(i)
			{
				ASSERT(HERE, i < j, "Excessive trailing zeros in MI64_GCD!");
				mi64_shrl(uv_ptr[k1], uv_ptr[k1], i, lenU);	lenU = mi64_getlen(uv_ptr[k1], lenU);
			}
		#if GCD_DEBUG>=1
			if(gcd_debug)
			{
				printf("pass %u, C: left-shifting vector %u by %u-%u-1 = %u bits;\n",pass,k2,j,i,j-i-1);
			//	printf("Input  = %s\n",&str_10k[__convert_mi64_base10_char(str_10k, 10<<10, uv_ptr[k2], lenV, 0)]);
			}
		#endif
			// Must use lenU as shift-operand length, since result may have a carry into next-higher word, which is discarded otherwise.
			mi64_shl(uv_ptr[k2], uv_ptr[k2], j-i-1, lenU);	lenV = mi64_getlen(uv_ptr[k2], lenU);
		}
	}

  #if GCD_DEBUG>=2
	if(gcd_debug && pass>=0)
	{
		printf("Left-justified Inputs to GCD are:\n");
		printf("                         U                   V  \n");
		for(i=0; i<lenU; i++)
		{
			if((u[i] != 0) || (v[i] != 0))
			{
				convert_uint64_base10_char(string1, u[i]);
				convert_uint64_base10_char(string2, v[i]);
				printf("i = %8u:  %s   %s\n",i,string1,string2);
			}
		}
		printf("\n");
	}
  #endif

#endif	/* #if(1) */

	/*
	Get leading 128 bits of U and V (starting with the MSB of max(U, V))
	to form the startvalues for the Lehmer GCD:
	*/
	nshift = leadz64(uv_ptr[k1][lenU-1]);
	ASSERT(HERE, nshift < 64,"gcd_lehmer.c: nshift < 64");

	mi64_extract_lead128(uv_ptr[k1], lenU, nshift, lead_u);
	mi64_extract_lead128(uv_ptr[k2], lenU, nshift, lead_v);

#if 0	/****** OBSOLETE ******/
	#error Do not use this code!

	mvbits64(uv_ptr[k1][lenU-1],0        ,64-nshift,&range[0][1],nshift);	/* move nonzero bits of leading digit into high word (2) of 128-bit slot, left-justifying...*/
	mvbits64(uv_ptr[k2][lenU-1],0        ,64-nshift,&range[1][1],nshift);

  if(lenU > 1)
  {
	mvbits64(uv_ptr[k1][lenU-2],64-nshift,nshift   ,&range[0][1],0     );	/* if leading digit had any leftmost zero bits, fill low-order bits of high word with high bits of next digit...*/
	mvbits64(uv_ptr[k2][lenU-2],64-nshift,nshift   ,&range[1][1],0     );

	mvbits64(uv_ptr[k1][lenU-2],0        ,64-nshift,&range[0][0],nshift);	/* move leading bits of (lenU-1)st digit into low word (1) of 128-bit slot, left-justifying...*/
	mvbits64(uv_ptr[k2][lenU-2],0        ,64-nshift,&range[1][0],nshift);
  }

  if(lenU > 2)
  {
	mvbits64(uv_ptr[k1][lenU-3],64-nshift,nshift   ,&range[0][0],0     );	/* if necessary, fill low-order bits of low word with high bits of (lenU-2)nd digit. */
	mvbits64(uv_ptr[k2][lenU-3],64-nshift,nshift   ,&range[1][0],0     );
  }
#endif

#if GCD_DEBUG>=2
if(gcd_debug)
{
	printf("Indices into U/V = %u, %u\n",k1,k2);
	printf("\n");

	printf("NSHIFT = %u:\n", nshift);
	printf("\n");
	convert_uint64_base10_char(string1, lead_u[0]);
	convert_uint64_base10_char(string2, lead_u[1]);
	printf("leadU[1,0] = %s   %s\n",string2,string1);

	convert_uint64_base10_char(string1, lead_v[0]);
	convert_uint64_base10_char(string2, lead_v[1]);
	printf("leadV[1,0] = %s   %s\n",string2,string1);
	printf("\n");
}
else if(0)
	printf("Indices into U/V = %u, %u\n",k1,k2);
#endif

/*	This is the loop that calculates the scalar multipliers for the Lehmer scheme: */

	/*	Initialize 2 x 2 matrix coefficients, difference and range arrays:	*/

	coeff[0][0] = coeff[1][1] = 1;		/* x,y multipliers here;		*/
	coeff[0][1] = coeff[1][0] = 0;		/* init to 2x2 identity matrix. */
	diff [0] = diff [1] = 1;	/* accumulated error here. */

	range[0][0] = lead_u[0];	range[0][1] = lead_u[1];
	range[1][0] = lead_v[0];	range[1][1] = lead_v[1];

	i = 0;	/* initial row A index */
	j = 1;	/* initial row B index */

	for(;;)	/* CALC_MULTS */
	{
	/*	use FP approximation to diagonal terms and FDIV to form integer multiplier. */

		/* calculate range + diff (128-bit sum), store in (t1,t2). */

		t1 = range[j][0] + diff[j];
		t2 = range[j][1] + (t1 < range[j][0]);	/* if (range[j,0] + diff[j] < range[j,0]), had a carry out of the 64-bit add. */

		num = (double)range[i][1]*rbase + (double)range[i][0];
		den = (double)t2         *rbase + (double)t1;

		/* If zero denominator, break: */
		if(den == 0)
		{
		#if GCD_DEBUG>=2
			if(gcd_debug)
			{
				printf("EXIT oo\n");
			}
		#endif
			break;	/* EXIT CALC_MULTS */
		}

		frac = num/den;
		iop = (uint64)frac;	// Use FDIV to calculate FP approximate to ratio, then integer truncate upon store.
	#if GCD_DEBUG>=1
		if(gcd_debug) {
			printf("frac = %20.15e, iop = %20llu\n",frac,iop);
		}
	#endif

		// If FP result is just a smidge < 1, force = 1 to prevent premature loop exit for nearly-indeitcal 128-bit num,den values:
		if(iop == 0 && frac > 0.999999999999999) {
			iop = 1;
			#if GCD_DEBUG>=1
				if(gcd_debug) {
				//	printf("frac = %20.15e, iop = %20llu\n",frac,iop);
					printf("n1 = %20llu, n0 = %20llu\n",range[i][1],range[i][0]);
					printf("d1 = %20llu, d0 = %20llu\n",t2,t1);
				}
			#endif
		}

		/* If ratio = 0 or is out of range of uint64, break: */
		if(iop == 0)
		{
		#if GCD_DEBUG>=1
			if(gcd_debug) printf("EXIT 0\n");
		#endif
			break;	/* EXIT CALC_MULTS */
		}
		else if(frac >= TWO63FLOAT)	// Q: Should we require ratio to be in range of double-float mantissa instead?
		{
		#if GCD_DEBUG>=1
			if(gcd_debug) printf("EXIT oo\n");
		#endif
			break;	/* EXIT CALC_MULTS */
		}

	/*	  Calculate diff[i] =  diff[i] + iop* diff(j), checking for integer overflow: */

	#ifdef MUL_LOHI64_SUBROUTINE
		MUL_LOHI64(iop,diff[j],&lo,&hi);
	#else
		MUL_LOHI64(iop,diff[j], lo, hi);
	#endif
		t3 = lo + diff[i];
		/* If low part carries on add or high part>0, diff[i] + iop* diff(j) > 2^64 */
		if((t3 < lo) || hi!=0)
		{
		#if GCD_DEBUG>=2
			if(gcd_debug)
			{
				printf("EXIT A\n");
			}
		#endif
			break;	/* EXIT CALC_MULTS */
		}

		diff[i] = t3;

	/*
	Similarly calculate coeff[i] = coeff[i] + iop*coeff(j). Since carries
	in matrix-vector multiply routine are signed, must limit coefficients to be < 2^63.
	For floating emulation of MULH64, limit coefficients to be < 2^52.
	*/
	#ifdef MUL_LOHI64_SUBROUTINE
		MUL_LOHI64(iop,coeff[j][0],&lo,&hi);
	#else
		MUL_LOHI64(iop,coeff[j][0], lo, hi);
	#endif
		t3 = lo + coeff[i][0];
		hi = hi + (t3<lo);
		if(hi!=0 || (t3>>max_bits) != 0)
		{
		#if GCD_DEBUG>=2
			if(gcd_debug)
			{
				printf("EXIT B\n");
			}
		#endif
			break;	/* EXIT CALC_MULTS */
		}

	#ifdef MUL_LOHI64_SUBROUTINE
		MUL_LOHI64(iop,coeff[j][1],&lo,&hi);
	#else
		MUL_LOHI64(iop,coeff[j][1], lo, hi);
	#endif
		t4 = lo + coeff[i][1];
		hi = hi + (t4<lo);
		if(hi!=0 || (t4>>max_bits) != 0)
		{
		#if GCD_DEBUG>=2
			if(gcd_debug)
			{
				printf("EXIT C\n");
			}
		#endif
			break;	/* EXIT CALC_MULTS */
		}

		coeff[i][0] = t3;
		coeff[i][1] = t4;

	/*	update row A of range:	*/

		/* calculate iop*(range(j) + diff(j)). This is guaranteed not to exceed 128 bits. */

	#ifdef MUL_LOHI64_SUBROUTINE
		MUL_LOHI64(iop,t1,&t3,&t4);
	#else
		MUL_LOHI64(iop,t1, t3, t4);
	#endif
		t1 = t3;			/*  low part = MULL64(iop,lo)	*/
		t2 = iop*t2 + t4;	/* high part = MULL64(iop,hi) + MULH64(iop,lo) */

		/* range[i] = range[i] - iop*(range(j) + diff(j)). */

		t1 = range[i][0] - t1;
		range[i][1] = range[i][1] - t2 - (t1 > range[i][0]);	/* see if had a borrow. */
		range[i][0] = t1;

	/*	interchange row indices: */

		k = i; i = j; j = k;

#if GCD_DEBUG>=3
if(gcd_debug)
{
	convert_uint64_base10_char(string1, coeff[0][0]);
	convert_uint64_base10_char(string2, coeff[0][1]);
	printf("\n");
	printf("ABmul=%s   %s\n",string1,string2);

	convert_uint64_base10_char(string1, coeff[1][0]);
	convert_uint64_base10_char(string2, coeff[1][1]);
	printf("CDmul=%s   %s\n",string1,string2);

	printf("i,j = %d, %d\n",i,j);
	printf("\n");
}
#endif
	}	/* enddo CALC_MULTS */

	ASSERT(HERE, i == (j^0x00000001),"gcd_lehmer.c: i == (j^0x00000001)");

	max_abcd = 0;	/* Allows us to verify that our max. multiplier >= 1 - we allow equality
					e.g. if our scalar-multiplier matrix has form

							( a  -b )   ( 1   0 )
							(       ) = (       )
							( c  -d )   ( 1  -1 ) .

					i.e. mainly ensure that we're not doing a trivial multiply by the identity matrix. */

	abmul[k1] = coeff[i][0];	max_abcd = MAX(abmul[k1], max_abcd);	/* a is here */
	abmul[k2] = coeff[i][1];	max_abcd = MAX(abmul[k2], max_abcd);	/* b is here */

	cdmul[k1] = coeff[j][0];	max_abcd = MAX(cdmul[k1], max_abcd);	/* c is here */
	cdmul[k2] = coeff[j][1];	max_abcd = MAX(cdmul[k2], max_abcd);	/* d is here */

#if GCD_DEBUG >= 1
if(gcd_debug)
{
	printf("a = %20llu; b = %20llu;\n", abmul[ 0], abmul[ 1]);
	printf("c = %20llu; d = %20llu;\n", cdmul[ 0], cdmul[ 1]);
}
#endif

	if(max_abcd <= 1)
	{
		/* Make sure there are at least 3 nonzero matrix coefficients: */
		if((abmul[k1]+abmul[k2] + cdmul[k1]+cdmul[k2]) < 3)
		{
			sprintf(string0, "2 or more zero matrix coefficients: %u, %u, %u, %u.\n", (uint32)abmul[k1], (uint32)abmul[k2], (uint32)cdmul[k1], (uint32)cdmul[k2]);
			ASSERT(HERE, 0, string0);
		}
	}

	/*	Now calculate the matrix-vector product

	( u')	   ( a  -b ) ( u )	   ( |a*u - b*v| )
	(   ) =	ABS(       )*(   )	== (             )
	( v')	   ( c  -d ) ( v )	   ( |c*u - d*v| )
	*/
	tmp32 = matrix_vector_product_sub(abmul, cdmul, uv_ptr, lenU);

	/* Bottom bit of return value indicates whether either result vector
	had greater magnitude than the input vectors - should be 0 in this case: */
	ASSERT(HERE, (tmp32&1) == 0, "gcd_lehmer: unexpected exit carry from matrix_vector_product_sub!");

	/*
	Store any sign-flip that may have occurred (e.g. the routine detected that
	the full-length scalar*vector product a*u - b*v was actually < 0, and as a
	result returned the inverted b*v - a*u in the u-vector slot) in j1 and j2:
	*/
	j1 = (tmp32>>(k1+1))&1;
	j2 = (tmp32>>(k2+1))&1;

#if GCD_DEBUG>=2
/* Insert a nonzero pass number here to enable pass-specific debug: */
if(gcd_debug && pass==0)
{
	printf("Result of matrix_vector_product:\n");
	printf("                         U                   V  \n");
	for(i=0; i<lenU; i++)
	{
		convert_uint64_base10_char(string1, u[i]);
		convert_uint64_base10_char(string2, v[i]);
		printf("i = %8u:  %s   %s\n",i,string1,string2);
	}
	printf("\n");
}
#endif

	if(EGCD)
	{
	#if 0/*GCD_DEBUG>=2*/
		/* Compute the leading 128 bits of the
		A' = (a*A + b*C), B' = (a*B + b*D), C' = (c*A + d*C), D' = (c*B + d*D)
		products and the leading 256 bits of the resulting A'*u,B'*v,C'*u,D'*v terms.
		Then repeat, but for the subtractive
		A" = (a*A - b*C), B" = (a*B - b*D), C" = (c*A - d*C), D" = (c*B - d*D) terms.
		This is to see if we want to + or - in updating the eGCD multipliers.
		*/
		if(gcd_debug)
		{
			/* Leading 128 bits of A/B: */
			if(mi64_cmpugt(A, B, len_abcd)
			{
				tmp32 = mi64_getlen(A, len_abcd);
				i = leadz64(A, tmp32);
			}
			else
			{
				tmp32 = mi64_getlen(B, len_abcd);
				i = leadz64(B, tmp32);
			}
			mi64_extract_lead128(A, tmp32, i, a128);
			mi64_extract_lead128(B, tmp32, i, b128);

			/* Leading 128 bits of C/D: */
			if(mi64_cmpugt(C, D, len_abcd)
			{
				tmp32 = mi64_getlen(C, len_abcd);
				i = leadz64(C, tmp32);
			}
			else
			{
				tmp32 = mi64_getlen(D, len_abcd);
				i = leadz64(D, tmp32);
			}
			mi64_extract_lead128(C, tmp32, i, c128);
			mi64_extract_lead128(D, tmp32, i, d128);

			/* A' = (a*A + b*C)*/
			mi64_mul_scalar(a128, ab_mul[0], x192, 2);
			mi64_mul_scalar(c128, ab_mul[1], y192, 2);
			mi64_add(x192,y192,x192);
			ap128[0]=x192[1];	ap128[1]=x192[2];	/* extract leading 128 bits of sum */
			/* B' = (a*B + b*D)*/
			mi64_mul_scalar(b128, ab_mul[0], x192, 2);
			mi64_mul_scalar(d128, ab_mul[1], y192, 2);
			mi64_add(x192,y192,x192);
			bp128[0]=x192[1];	bp128[1]=x192[2];	/* extract leading 128 bits of sum */
			/* C' = (c*A + d*C)*/
			mi64_mul_scalar(a128, cd_mul[0], x192, 2);
			mi64_mul_scalar(c128, cd_mul[1], y192, 2);
			mi64_add(x192,y192,x192);
			cp128[0]=x192[1];	cp128[1]=x192[2];	/* extract leading 128 bits of sum */
			/* D' = (c*B + d*D)*/
			mi64_mul_scalar(b128, cd_mul[0], x192, 2);
			mi64_mul_scalar(d128, cd_mul[1], y192, 2);
			mi64_add(x192,y192,x192);
			dp128[0]=x192[1];	dp128[1]=x192[2];	/* extract leading 128 bits of sum */

			/* Now calculate & print the leading terms of |A'*u-B'*v|, |C'*u-D'*v|: */

***************************************************************************
/***No - this won't work, since we no longer have the original U and V!***/
		}
	#endif
		/*
		If doing an eGCD, update the 2 columns of the eGCD multiplier array similarly.
		The basic computation is of the form

			(  A' -B')   ( a  -b ) (  A  -B )   ( (a*A + b*C)  -(a*B + b*D) )
			(        ) = (       )*(        ) = (                           )
			( -C'  D')   ( c  -d ) ( -C   D )   ( (c*A + d*C)  -(c*B + d*D) )

		The indices j1 and j2 tell us whether the elements of the (a, -b) and (c, -d)
		pairs have been flipped (e.g. j1 = 1 means that the first row of our scalar-multipliers
		matrix was ( -b  a ) rather than ( a  -b )), i.e. whether we need to flip the corresponding
		sign bit of sign_abcd. Since the sign-flips indicated by j1 and j2 can be overridden post hoc
		during the execution of matrix_vector_product_sub, we XOR the a priori sign-flip bits
		(j1 +2*j2) with the post-hoc-override bits (stored in tmp32) to get the final flip-values:
		*/
		sign_abcd ^= (j1 + (j2 << 1));

		/* The u/v-update call to matrix_vector_product_sub normally expects
		a single sign flip in computing |a*u-b*v| and |c*u-d*v|, in which case
		we call matrix_vector_product_add to update the eGCD multipliers. OTOH
		if there were 0 or 2 sign flips in the u/v-update, call matrix_vector_product_sub
		to update eGCD mutlipliers: */

/*****HOW TO PROPERLY DEAL WITH CASE WHERE A*u-B*v or C*u-D*v = 0?******/

		if(egcd_sign_swap)
		{
		#if GCD_DEBUG>=2
			if(gcd_debug)
			{
				printf("Pass = %u\n", pass);
				printf("A/B eGCD multiplier vectors:\n");
				if(sign_abcd & 0x1)
				printf("                         -A                  +B  \n");
				else
				printf("                         +A                  -B  \n");

				printf("              base=2^64;x=1;A=B=0;\n");
				for(i=0; i<len_abcd; i++)
				{
					convert_uint64_base10_char(string1, A[i]);
					convert_uint64_base10_char(string2, B[i]);
					printf("i = %8u:  A+=%s*x; B+=%s*x; x*=base;\n",i,string1,string2);
				}
				printf("\n");

				printf("C/D eGCD multiplier vectors:\n");
				if(sign_abcd & 0x2)
				printf("                         -C                  +D  \n");
				else
				printf("                         +C                  -D  \n");

				printf("              base=2^64;x=1;C=D=0;\n");
				for(i=0; i<len_abcd; i++)
				{
					convert_uint64_base10_char(string1, C[i]);
					convert_uint64_base10_char(string2, D[i]);
					printf("i = %8u:  C+=%s*x; D+=%s*x; x*=base;\n",i,string1,string2);
				}
				printf("\n");
			}	/* endif(!gcd_debug) */
		#endif

			/* Don't combine the 2 calls into an if(call1() || call2()),
			since we need to be sure both function calls always occur: */
			tmp32  = matrix_vector_product_sub(abmul, cdmul, ac_ptr, len_abcd);
			tmp32 |= matrix_vector_product_sub(abmul, cdmul, bd_ptr, len_abcd);
			len_abcd += (tmp32&1);
		}
		else
		{
			tmp32  = matrix_vector_product_add(abmul, cdmul, ac_ptr, len_abcd);
			tmp32 |= matrix_vector_product_add(abmul, cdmul, bd_ptr, len_abcd);
			len_abcd += (tmp32&1);
		}

		/* Lookahead scheme? */
		if(j1 ^ j2)
			egcd_sign_swap = FALSE;
		else
			egcd_sign_swap = TRUE;

#if GCD_DEBUG>=2
		/* Insert a nonzero pass number here to enable pass-specific debug: */
		if(gcd_debug)
		{
			printf("Pass = %u\n", pass);
			printf("A/B eGCD multiplier vectors:\n");
			if(sign_abcd & 0x1)
			printf("                         -A                  +B  \n");
			else
			printf("                         +A                  -B  \n");

			printf("              base=2^64;x=1;A=B=0;\n");
			for(i=0; i<len_abcd; i++)
			{
				convert_uint64_base10_char(string1, A[i]);
				convert_uint64_base10_char(string2, B[i]);
				printf("i = %8u:  A+=%s*x; B+=%s*x; x*=base;\n",i,string1,string2);
			}
			printf("\n");

			printf("C/D eGCD multiplier vectors:\n");
			if(sign_abcd & 0x2)
			printf("                         -C                  +D  \n");
			else
			printf("                         +C                  -D  \n");

			printf("              base=2^64;x=1;C=D=0;\n");
			for(i=0; i<len_abcd; i++)
			{
				convert_uint64_base10_char(string1, C[i]);
				convert_uint64_base10_char(string2, D[i]);
				printf("i = %8u:  C+=%s*x; D+=%s*x; x*=base;\n",i,string1,string2);
			}
			printf("\n");
		}
	#endif

	}

/**** 8/23/2005: disable for now, since this really complicates things when doing eGCD. ****/
#if 0
	uint32 shift_ab = 0, shift_cd = 0;
	/* Right-justify U and V, one of which (but not both)
	might be even as a result of the matrix_vector_product call.
	Must do this to prevent accumulation of spurious powers of 2
	in the gcd, i.e. if we see at the start of the next pass through
	the main for-loop that V needs left-justifying, we need to ensure
	that U is odd, otherwise the gcd will accumulate a spurious 2:
	*/
	/* lenU,B remain unchanged: */
	tz0 = mi64_trailz(u, lenU);
	if(tz0)
	{
		mi64_shrl(u, u, tz0, lenU);
shift_cd += tz0;
		if(EGCD)
		{
			/* Left-shift the complement-vector eGCD multipliers (C/D) by tz0 places: */
			mi64_shl(ac_ptr[1], ac_ptr[1], tz0, len_abcd);
			mi64_shl(bd_ptr[1], bd_ptr[1], tz0, len_abcd);	len_abcd = MAX(mi64_getlen(ac_ptr[1], len_abcd) , mi64_getlen(bd_ptr[1], len_abcd));
		}
	}

	tz1 = mi64_trailz(v, lenU);
	if(tz1)
	{
		mi64_shrl(v, v, tz1, lenU);
shift_ab += tz1;
		if(EGCD)
		{
			/* Left-shift the complement-vector eGCD multipliers by (A/B) tz1 places: */
			mi64_shl(ac_ptr[0], ac_ptr[0], tz1, len_abcd);
			mi64_shl(bd_ptr[0], bd_ptr[0], tz1, len_abcd);	len_abcd = MAX(mi64_getlen(ac_ptr[0], len_abcd) , mi64_getlen(bd_ptr[0], len_abcd));
		}
	}
#endif

	/*	set k1 to row index of larger element: */
	istart = lenU-1;
	lenU = 0;	/* This allows us to determine whether lenU has already been set on a previous pass in the ensuing loop. */
	for(i=istart; i>=0; i--)
	{
		if((u[i] != 0 || v[i] != 0) && lenU == 0)
			lenU = i+1;
		/* Even if we've already found the leading nonzero element,
		in order to properly set k1 and k2 continue until find a differing element:
		*/
		if(u[i] == v[i])	//<*** On pass 2, get upper 2 slots = 0 as expected, but next 350 or so all have u[i] = v[i] = 14757395258967641292 ***
			continue;		//That leads to u/v having leading 350 64-bit words identical, i.e. on next cal-mults pass expect a linear combo = u-v
		else				//in on of the slots, other should be u or v unchanged, i.e. scalar-mults matrix = 1  -1 / 1   0 .
		{
			k1 = (u[i] < v[i]);
			k2 = k1 ^ 0x00000001;
			break;
		}
	}
	/*	If U = V, we're done, with a copy of a nontrivial GCD in each of the 2 vectors. */
	if(i == -1)
		goto GCD_DONE;

	/* Otherwise, parse the smaller array until find leading nonzero element: */
	x = uv_ptr[k2];
	for(i=lenU-1; i>=0; i--)
	{
		if(x[i] != 0)
			break;
	}
	lenV = i + 1;
	/* If lenV = 0, we're done. */
	if(lenV == 0)
		goto GCD_DONE;

#if GCD_DEBUG>=1
if(gcd_debug || (pass&1023)==0) {
	printf("PASS %d: lenU = %u, lenV = %u, k1,k2 = %d,%d\n",pass,lenU,lenV,k1,k2);
//	if(gcd_debug && (pass%5)==0) {
//		printf("");
//	}
}
#endif
}	/* enddo MAIN	*/

GCD_DONE:

	/*	Print GCD: */
#if GCD_DEBUG>=0
	printf("GCD finished in %u passes.\n", pass);
#endif
	ASSERT(HERE, (lenU+lenV != 0), "mi_gcd: final U & V-vectors both zero-length!");

	/* Put copies of the GCD into both uv_ptr[k1] and uv_ptr[k2]: */
	for(i = 0; i < lenU; i++)
	{
		uv_ptr[k2][i] = uv_ptr[k1][i];
	}

#if GCD_DEBUG>=2
	if(EGCD && gcd_debug)
	{
		ASSERT(HERE, min_trailing_zeros == 0, "min_trailing_zeros != 0");

		printf("A/B eGCD multiplier vectors:\n");
		if(sign_abcd & 0x1)
		printf("                         -A                  +B  \n");
		else
		printf("                         +A                  -B  \n");

		printf("              base=2^64;x=1;A=B=0;\n");
		for(i=0; i<len_abcd; i++)
		{
			convert_uint64_base10_char(string1, A[i]);
			convert_uint64_base10_char(string2, B[i]);
			printf("i = %8u:  A+=%s*x; B+=%s*x; x*=base;\n",i,string1,string2);
		}
		printf("\n");

		printf("C/D eGCD multiplier vectors:\n");
		if(sign_abcd & 0x2)
		printf("                         -C                  +D  \n");
		else
		printf("                         +C                  -D  \n");

		printf("              base=2^64;x=1;C=D=0;\n");
		for(i=0; i<len_abcd; i++)
		{
			convert_uint64_base10_char(string1, C[i]);
			convert_uint64_base10_char(string2, D[i]);
			printf("i = %8u:  C+=%s*x; D+=%s*x; x*=base;\n",i,string1,string2);
		}
		printf("\n");
	}
#endif

	/* Restore any common power of 2: */
	if(min_trailing_zeros > 0)
	{
		lenU += (min_trailing_zeros>>6)+1;
		mi64_shl(uv_ptr[k1], uv_ptr[k1], min_trailing_zeros, lenU);
		lenU = mi64_getlen(uv_ptr[k1], lenU);	/* 12/16/08: getlen() had lenV here. (?) */
	}

	/* Currently printing functions can only handle up to STR_MAX_LEN digits: */
	if(lenU*log10_base >= STR_MAX_LEN)
	{
		printf("GCD exceeds printable string length - printing in raw base-2^64 form, least-significant first:\n");
		printf("GCD = 2^%u * [\n", min_trailing_zeros);
		for(i = 0; i < lenU; i++)
		{
			convert_uint64_base10_char(string0, uv_ptr[k1][i]);
			printf("i = %8u:  %s\n",i,string0);
		}
		printf("]\n");
	}
	else
	{
		printf("GCD = %s\n", &string0[convert_mi64_base10_char(string0, uv_ptr[k1], lenU, 0)]);
	}

	/*	********TODO:********** base-3 PRP test	*/

	return lenU;
}

/************* Helper functions: *************/

/*	Calculates the matrix-vector product

		( u' )      ( a  -b ) ( u )
		(    ) = ABS(       )*(   )
		( v' )      ( c  -d ) ( v ) ,

	where the absolute-value applies to the separate linear combinations resulting
	from the matrix-vector multiply.

	Here, a, b, c and d are unsigned 64-bit integers having magnitude <= 2^63 and
	u and v are pointers to base-2^64 multiword integer arrays assumed to contain at most
	(len) nonzero elements. The routine examines the inputs and decides how to order
	the a*u/b*v and c*u/d*v scalr*vector multiplies so as to make the result of each
	nonnegative - since the criterion it uses to do so examines only the leading 192
	bits of each of the 4 aforementioned products it is not completely deterministic,
	and there is postprocessing code to fix things up if the 192-bit ordering-guess
	proves incorrect. (But in general we hope the 192-bit guess is correct most of the time.)

	For each row-times-vector inner product defined by the above matrix-vector multiply, the
	u'/v'-pointers point (in terms of the main fixed-address u/v data arrays defined in the
	calling function) either to (u, v) or (v, u), as determined by the arguments ii1 for the
	(a, -b) row and ii2 for the (c, -d) row multiplies, respectively. This mechanism is designed
	so that via the above-described small amount of preprocessing, ii1 and ii2 can be set to
	greatly enhance the odds) that each of the resulting row-times-vector inner products
	will be nonnegative. In terms of the fixed-address u/v data arrays, the 4 possible values
	of the (ii1, ii2) bit pair define 4 possibilities, as follows:

							(  a  -b ) ( u )
	(ii1, ii2) = (0, 0):	(        )*(   )
							(  c  -d ) ( v )

							(  a  -b ) ( u )
	(ii1, ii2) = (0, 1):	(        )*(   )
							( -d   c ) ( v )

							( -b   a ) ( u )
	(ii1, ii2) = (1, 0):	(        )*(   )
							(  c  -d ) ( v )

							( -b   a ) ( u )
	(ii1, ii2) = (1, 1):	(        )*(   )
							( -d   c ) ( v ) .

	In our implementation below we define indices j1 = ii1, j2 = 1 - j1, j3 = ii2, j4 = 1 - j3
	and access u and v via uv_ptr[0] = &u and uv_ptr[1] = &v, in terms of which our 2 output
	vectors (which overwrite the inputs) are:

		u" = abmul[j1]*( uv_ptr[j1] ) - abmul[j2]*( uv_ptr[j2] )

		v" = cdmul[j3]*( uv_ptr[j3] ) - cdmul[j4]*( uv_ptr[j4] )

	(the aforementioned calling-function preprocessing comes into play here, in that it is
	desirable that the user has arranged the scalar multipliers a,b,c,d such that u" and v" >= 0.)

	EFFECTS:
			- Overwrites the input vectors with the result vectors.

	RETURNS: the following result-value bits are meaningful:

			<0:0> = 1 if either of the |a*u-b*v| or |c*u-d*v| linear combination
			  had a nonzero output carry (which is placed in the respective u,v[len]
			  vector slot), which should only occur when the routine is being called
			  to update eGCD multipliers - in this case caller should increment (len).

			  In the case where the routine is called to update GCD vectors (i.e.
			  the outputs have size less than the inputs, leaves any needed zero-leading-
			  element detection and subsequent vector-length-decrementing to the caller.

			<1:2> = 1 if the respective a*u-b*v or c*u-d*v linear combination
			  needed a negation in order to make the result nonnegative.
*/
uint32 matrix_vector_product_sub
(
	uint64c abmul[],	/* It is assumed here that abmul[0,1] multiply uv_ptr[0,1]! */
	uint64c cdmul[],	/* It is assumed here that cdmul[0,1] multiply uv_ptr[0,1]! */
	uint64 *uv_ptr[],
	uint32  len
)
{
    uint32 i, j1,j2,j3,j4, nshift, retval = 0;
    sint64 cy1 = 0, cy2 = 0;	/* Carries */
	uint64 cylo, hi1, hi2, lo1, lo2, lo;
	int cyhi;
    sint64 a, b, c, d;
	uint64 *u, *v, *A, *B, *C, *D;
	uint64 lead_u[2], lead_v[2];
	uint64 val1, val2, val3, val4;

	/* Using absolute indices here (rather than k1/k2) during u/v-update calls
	eases later eGCD bookkeeping; */
	u = uv_ptr[0];
	v = uv_ptr[1];

	/* Find larger of the 2 vectors in order to properly set shift
	count prior to leading-bits extraction: */
	if(mi64_cmpult(v, u, len))
	{
		nshift = leadz64(u[len-1]);	/* If U > V, use it to set shift count... */
	}
	else
	{
		nshift = leadz64(v[len-1]);	/* ...otherwise use V to set shift count. */
	}

	/*	Figure out which order to do the multiplies in (i.e. a*u - b*v or b*v - a*u)
	to ensure (or at least maximize the chance, as O(len)-costly postprocessing is
	needed otherwise) that result >= 0:
	*/

	/*	If leading 192 bits of a*u < b*v, swap a,b and the ptrs to u,v: */
	mi64_extract_lead128(u,len, nshift, lead_u);
/*	ASSERT(HERE, (lead_u[0]+lead_u[1]) != 0, "(lead_u[0]+lead_u[1]) == 0!");	*/

	mi64_extract_lead128(v,len, nshift, lead_v);
/*	ASSERT(HERE, (lead_v[0]+lead_v[1]) != 0, "(lead_v[0]+lead_v[1]) == 0!");	*/

	j1 = 0;
	if(CMP_LT_PROD192(abmul[0], lead_u[0], lead_u[1], abmul[1], lead_v[0], lead_v[1]))
	{
	#if GCD_DEBUG>=2
		if(gcd_debug)
		{
			printf("(a*u) < (b*v) : swapping these\n");
		}
	#endif
		j1 ^= 0x00000001;
		retval ^= 0x00000002;
	}
	j2 = j1^0x00000001;

/*	Now do same for c and d: */

	j3 = 0;
	if(CMP_LT_PROD192(cdmul[0], lead_u[0], lead_u[1], cdmul[1], lead_v[0], lead_v[1]))
	{
	#if GCD_DEBUG>=2
		if(gcd_debug)
		{
			printf("(c*u) < (d*v) : swapping these\n");
		}
	#endif
		j3 ^= 0x00000001;
		retval ^= 0x00000004;
	}
	j4 = j3^0x00000001;

	/*
	Proceed with the pairwise-linear-combination step, writing results back into u/v:
	*/
    a = abmul[j1];
    b = abmul[j2];
    c = cdmul[j3];
    d = cdmul[j4];

	/* Check magnitudes of a,b,c,d < 2^63: */
	ASSERT(HERE, a >= 0 && b >= 0 && c >= 0 && d >= 0,"matrix_vector_product_sub: a >= 0 && b >= 0 && c >= 0 && d >= 0");

	A = uv_ptr[j1];
	B = uv_ptr[j2];
	C = uv_ptr[j3];
	D = uv_ptr[j4];

	// 04 May 2012 x86_64 timing experiments: 2x loop-unroll was worse, rejiggering MULs as noted below helped:
    for (i = 0; i < len; i++)
    {
    	/* Make copies of the current 4 vector elts since some of these will get clobbered: */
		val1 = A[i];
		val2 = B[i];
	#ifdef MUL_LOHI64_SUBROUTINE	// Moving the first pair of MULs from just below the val4 = to here gave a 7% speedup
		MUL_LOHI64(a, val1,&lo1,&hi1);
		MUL_LOHI64(b, val2,&lo2,&hi2);
	#else
		MUL_LOHI64(a, val1, lo1, hi1);
		MUL_LOHI64(b, val2, lo2, hi2);
	#endif
		val3 = C[i];
		val4 = D[i];

	/* a*u[i] - b*v[i] */
		/*
		Compute cy1 + (hi1, lo1) - (hi2, lo2)
			= (cyhi, cylo) - [(hi2, lo2) - (hi1, lo1)] .
		Note cy1 is signed.
		*/
		cylo = (uint64)cy1;
		cyhi = -(cy1 < 0);	/* Sign extension of cy1 */

		lo = lo2 - lo1;
		cyhi = cyhi - (cylo < lo) + (lo2 < lo1);
		cy1 = (sint64)cyhi + (sint64)(hi1 - hi2);

		u[i] = cylo - lo;

/*	printf("difference = %lld + %lld *2^64\n", u[i], cy1);	*/

	/* c*u[i] - d*v[i] */
	#ifdef MUL_LOHI64_SUBROUTINE
		MUL_LOHI64(c, val3,&lo1,&hi1);
		MUL_LOHI64(d, val4,&lo2,&hi2);
	#else
		MUL_LOHI64(c, val3, lo1, hi1);
		MUL_LOHI64(d, val4, lo2, hi2);
	#endif
		/*
		Compute cy2 in same fashion as cy1 above:
		*/
		cylo = (uint64)cy2;
		cyhi = -(cy2 < 0);	/* Sign extension of cy2 */

		lo = lo2 - lo1;
		cyhi = cyhi - (cylo < lo) + (lo2 < lo1);
		cy2 = (sint64)cyhi + (sint64)(hi1 - hi2);
		v[i] = cylo - lo;

/*	printf("difference = %lld + %lld *2^64\n", v[i], cy2);	*/
    } /* for i */

#if GCD_DEBUG>=2
    printf("cy1 = %lld, cy2 = %lld, cy1|cy2 = %lld\n", cy1, cy2, cy1|cy2);
#endif

	/*	For GCD main-vector updates we know exit carries are 0 or -1,
	but to handle eGCD output carries (which may be as large as the difference of
	two 63-bit ints) need to generalize this to detect any negative 64-bit int: */
    if(cy1 < 0)
    {
#if GCD_DEBUG
    printf("cy1 < 0...complementing.\n");
#endif
		for (i = 0; i < len; i++)
		{
			u[i] = ~u[i];	     /* Ones' complement */
		}
		cy1 = ~cy1;

		for (i = 0; i < len; i++)
		{    /* Add -cy1 */
			u[i]++;
			if (u[i] != 0) break;
		}

		/* Post-processing negation flips any corresponding sign bit in the return field: */
		retval ^= 0x00000002;
    }

	if(cy2 < 0)
    {
#if GCD_DEBUG
    printf("cy2 < 0...complementing.\n");
#endif
		for (i = 0; i < len; i++)
		{
			v[i] = ~v[i];
		}
		cy2 = ~cy2;

		for (i = 0; i < len; i++)
		{    /* Add -cy2 */
			v[i]++;
			if (v[i] != 0) break;
		}

		retval ^= 0x00000004;
    }

	/* Does vector length need incrementing? */
    if(cy1 || cy2)
    {
		u[len] = cy1;
		v[len] = cy2;
		retval |= 1;
	}

	return retval;

} /* matrix_vector_product_sub */

/*	Calculates the matrix-vector product

		( u' )   ( a  b ) ( u )
		(    ) = (      )*(   )
		( v' )   ( c  d ) ( v ) .

	Here, a, b, c and d are unsigned 64-bit integers having magnitude <= 2^63 and
	u and v are pointers to base-2^64 multiword integer arrays assumed to contain at most
	(len) nonzero elements.

	EFFECTS:
			- Overwrites the input vectors with the result vectors.

	RETURNS:
			- 1 if either of the linear combinations had a nonzero output carry,
			  indicating to the caller that the vector length field needs incrementing
*/
uint32 matrix_vector_product_add
(
	uint64c abmul[],
	uint64c cdmul[],
	uint64 *uv_ptr[],
	uint32  len
)
{
    uint32 i;
    uint64 cyloA = 0, cyhiA = 0, cyloC = 0, cyhiC = 0;
    uint64 hi1, hi2, lo1, lo2;
	sint64c a = abmul[0], b = abmul[1], c = cdmul[0], d = cdmul[1];
	uint64 Ai, Ci;
	uint64 *u = uv_ptr[0], *v = uv_ptr[1];

	/* Check magnitudes of a,b,c,d < 2^63: */
	ASSERT(HERE, a >= 0 && b >= 0 && c >= 0 && d >= 0,"matrix_vector_product_add: a >= 0 && b >= 0 && c >= 0 && d >= 0");

    for (i = 0; i < len; i++)
    {
		Ai = u[i];
		Ci = v[i];

	/* a*u[i] + b*v[i] */
	#ifdef MUL_LOHI64_SUBROUTINE
		MUL_LOHI64(a, Ai,&lo1,&hi1);
		MUL_LOHI64(b, Ci,&lo2,&hi2);
	#else
		MUL_LOHI64(a, Ai, lo1, hi1);
		MUL_LOHI64(b, Ci, lo2, hi2);
	#endif
		/*
			(u[i], cylo, cyhi) <-- (lo1 + lo2 + cylo, hi1 + hi2 + cyhi, 0)	(3-word vector add-with carry)
		*/
		/* Low part, with carryin from previous loop execution: */
		lo1 += cyloA;	cyhiA += lo1 < cyloA;
		lo1 +=   lo2;	cyhiA += lo1 < lo2;
		u[i] = lo1;

		/* High part, with carrying from low part: */
		cyloA  = hi1 + hi2;		/* hi1, hi2 < 2^63, so this can't overflow */
		cyloA += cyhiA;	cyhiA = cyloA < cyhiA;	/* Carries into the next loop execution */

	/* c*u[i] + d*v[i] */
	#ifdef MUL_LOHI64_SUBROUTINE
		MUL_LOHI64(c, Ai,&lo1,&hi1);
		MUL_LOHI64(d, Ci,&lo2,&hi2);
	#else
		MUL_LOHI64(c, Ai, lo1, hi1);
		MUL_LOHI64(d, Ci, lo2, hi2);
	#endif
		/* Low part, with carryin from previous loop execution: */
		lo1 += cyloC;	cyhiC += lo1 < cyloC;
		lo1 +=   lo2;	cyhiC += lo1 < lo2;
		v[i] = lo1;

		/* High part, with carrying from low part: */
		cyloC  = hi1 + hi2;
		cyloC += cyhiC;	cyhiC = cyloC < cyhiC;	/* Carries into the next loop execution */

    } /* for i */

	/*	Make sure exit carries are at most one word long: */
	ASSERT(HERE, (cyhiA+cyhiC == 0),"Carryouts from matrix_vector_product_add(abmul, cdmul, uv_ptr...) out of range!");
	if(cyloA+cyloC != 0)
	{
		u[len] = cyloA;
		v[len] = cyloC;
		return 1;
	}
	return 0;

} /* matrix_vector_product_add */


/*******************/

/*
!   Unsigned compare of two 192-bit products: a*(xhi*2^64 + xlo) < b*(yhi*2^64 + ylo).
!   where all operands are treated as unsigned 64-bit integers.
!   If true, returns 1; otherwise returns zero.
*/
int	CMP_LT_PROD192(uint64 a, uint64 xlo, uint64 xhi, uint64 b, uint64 ylo, uint64 yhi)
{
	uint192 ax, by;
	uint64 tmp;

/*	calculate the two products. */
/*
!   a*x = a*(xhi*2^64 + xlo) = MULT_HIGH(a*xhi)*2^128 + [ a*xhi + MULT_HIGH(a*xlo) ]*2^64 + a*xlo
!                            = ax.d2*2^128 + ax.d1*2^64 + ax.d0 is here...
*/
	#ifdef MUL_LOHI64_SUBROUTINE
		MUL_LOHI64(a,xlo,&ax.d0,&  tmp);
		MUL_LOHI64(a,xhi,&ax.d1,&ax.d2);
	#else
		MUL_LOHI64(a,xlo, ax.d0,   tmp);
		MUL_LOHI64(a,xhi, ax.d1, ax.d2);
	#endif
	ax.d1 = ax.d1 + tmp;				/* middle part is here... */
	ax.d2 = ax.d2 + (ax.d1<tmp);	/* had a carry into high part. */
/*
!   b*y = b*(yhi*2^64 + ylo) = MULT_HIGH(b*yhi)*2^128 + [ b*yhi + MULT_HIGH(b*ylo) ]*2^64 + b*ylo
!                            = bx.d2*2^128 + bx.d1*2^64 + bx.d0 is here...
*/
	#ifdef MUL_LOHI64_SUBROUTINE
		MUL_LOHI64(b,ylo,&by.d0,&  tmp);
		MUL_LOHI64(b,yhi,&by.d1,&by.d2);
	#else
		MUL_LOHI64(b,ylo, by.d0,   tmp);
		MUL_LOHI64(b,yhi, by.d1, by.d2);
	#endif
	by.d1 = by.d1 + tmp;				/* middle part is here... */
	by.d2 = by.d2 + (by.d1<tmp);	/* had a carry into high part. */

/*	Unsigned compare of products. */
	return (int)CMPULT192(ax,by);
}

/*
!***************

	subroutine mv_dwtvarbase_to_int64(x,p,m,u,ndim)
!...Move an arbitrary-length integer, stored in balanced-digit, variable-base (DWT) form in an integer*4 array X(1:m),
!   into row 2 of the integer*8 array U in base-2^64 positive-digit form. The BASE_INDEX array tells us whether a given
!   digit of X is a DWT bigword or littleword during the conversion.
!
	implicit none
	integer, intent(in) :: ndim,m,p
	integer :: x(m)
!...
	integer :: bits0,i,j,k,slot,word,bw,sw,bjmodn
	integer :: bit_end,bit_start,bit_sum,excess_bits
	integer*8 :: u(ndim)
	integer*8 :: cy1,t1

	bw    = mod(p,m)	!* Number of bigwords in the Crandall/Fagin mixed-radix representation = (Mersenne exponent) mod (vector length).
	sw    = m - bw		!* Number of smallwords.
	bits0 = p/m		!* Number of bits in a smallword.

!...make sure the V vector has positive leading digit...

	do j = m,1,-1
	  if(x(j) /= 0) then
	    if(x(j) < 0) x = -x	!* flips the sign of the entire vector
	    EXIT
	  endif
	enddo

!...Now move the V vector into row 2 of the 64-bit integer U(:,:) array, converting to a fixed base along the way...

	bit_start = 0
	i = 1		!* put index is here
	k = bits0 + 1	!* lowest-order digit is always a bigword.

	do j = 1,m-1

	  cy1 = x(j)			!* this converts to positive-digit form sans branches:
	  t1 = ISHFT(cy1,-63)		!* If the sign bit = 1, t1 = 1...
	  x(j+1) = x(j+1) - t1		!* subtract one from next-higher digit...
	  cy1 = cy1 + ISHFT(t1,k)	!* and add the appropriate base to the current digit.

	  bit_end = bit_start + k
	  excess_bits = bit_end - 64
	  if(excess_bits > 0)then
	    mvbits64(cy1,0,k - excess_bits,u(i),bit_start)
	    i = i + 1
	    mvbits64(cy1,k - excess_bits,excess_bits,u(i),0)
	    bit_start = excess_bits
	  elseif(excess_bits == 0)then
	    mvbits64(cy1,0,k,u(i),bit_start)
	    i = i + 1
	    bit_start = 0
	  else
	    mvbits64(cy1,0,k,u(i),bit_start)
	    bit_start = bit_end
	  endif

	  bjmodn = iand(bjmodn + bw,m-1)	!* Can only use this speedy mod for n a power of 2...
	  k = bits0 + ishft(sw - bjmodn, -31)

	enddo

!...Do the leading digit separately...

	k = bits0		!* leading digit is always a smallword.

	cy1 = x(m)		!* this converts to positive-digit form sans branches:
	if(btest(cy1,63))then; print*,'FATAL: V-vector has negative leading digit.'; STOP; endif

	bit_end = bit_start + k
	excess_bits = bit_end - 64
	if(excess_bits > 0)then
	  mvbits64(cy1,0,k - excess_bits,u(i),bit_start)
	  i = i + 1
	  mvbits64(cy1,k - excess_bits,excess_bits,u(i),0)
	elseif(excess_bits == 0)then
	  mvbits64(cy1,0,k,u(i),bit_start)
	else
	  mvbits64(cy1,0,k,u(i),bit_start)
	endif

!...And zero any remaining slots.

	u(i+1:ndim)=0

	end subroutine mv_dwtvarbase_to_int64
*/

/*
Various small GCD self-tests.
*/
#ifdef GCD_STANDALONE
  int main(int argc, char *argv[])
  {
#else
  int test_gcd()
  {
#endif

	/* Adjust this value to reflect max. vector length desired for self-tests: */
	#define MAX_ARRAY_DIM	34420
	uint32 i, j, p, lenX, lenY, lenU, lenV, lenZ, vec_len, fft_len;
	int fail = 0;
	uint64 q,qinv,tmp, mod, rem, tmp64 = 0;
	uint64 *u, *v, *w, *x, *y, *z;
	double *a, *b, *c;

	/* Extra stuff needed for eGCD: */
	uint32 eGCD = FALSE;	// This can be toggled here for the entire self-test suite, or set/unset on a case-by-case basis.
	uint64 *Ap, *Bp;		// *To-Do*: Currently, none of the tests actually makes use of the eGCD multipliers - Need to remedy that.
	uint32 len_AB, sign_AB;

	/* Multiwird primes of various lengths: */
	const uint32 num_mword_prime = 10;	/* Define precomputed primes of lengths 1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100 words */
	const uint32 len_mword_prime[] = {1,2,3,4,5,6,7,8,9,10};

	const uint64 y1[1] = {
		16357897499336320021ull
	};

	const uint64 y2[2] = {
		14921004030550884877ull, 2503926509447010181ull
	};

	const uint64 y3[3] = {
		14096952260802720589ull,17612861172509999137ull, 3613444211655844428ull
	};

	const uint64 y4[4] = {
		15441773537038068193ull, 3295852079737058441ull, 4277869928288117018ull, 7590955667915091630ull
	};

	const uint64 y5[5] = {
		10957599942581074327ull, 2992213073808157724ull,16384404950923186766ull, 9717091901949065414ull, 8787226726310349388ull
	};

	const uint64 y6[6] = {
		 3363579288736505507ull,18296582175061969760ull,11592634605935338547ull,  850486128012099330ull, 3591220199434432035ull, 8391710716962706019ull
	};

	const uint64 y7[7] = {
		 6502386164572629557ull, 2945843944052902425ull,13516514962117353801ull,15895036866380037651ull,17457756735211593046ull,15977169653466540555ull,13709678387685797452ull
	};

	const uint64 y8[8] = {
		 3786291846004455523ull, 7689517826572499780ull, 2216874701508481392ull, 2224995263328296486ull,16239032412117211245ull, 2072578804862300364ull, 5338315421425562581ull,16881966094930949545ull
	};

	const uint64 y9[9] = {
		 3396535911186315951ull,17961240762910685599ull,11234599856160171073ull,14800252658289931881ull,  421966617845574006ull, 7709528122749317599ull,14770838913912070975ull,15919994111674772454ull, 1022190188411309751ull
	};

	const uint64 y10[10] = {
		 6154942771176964881ull,14657802710986676243ull, 6367391185176237819ull, 7220279878409007123ull,17437454528603741845ull,11135248797733030345ull, 9466072315348543570ull,16100604701688312144ull,10143506867524673304ull,	 174955034380965ull
	};

	#define w431	4757395258967641292ull
	const uint64 y431_test[431] = {
		15441773537038068193ull, 3295852079737058441ull, 4277869928288117018ull, 7590955667915091630ull,w431,14096952260802720589ull,
		17612861172509999137ull, 3613444211655844428ull,w431,14921004030550884877ull, 2503926509447010181ull,w431,16357897499336320021ull,
		w431, 8589934593ull, 17179869187ull, 25769803781ull, 34359738375ull, 42949672969ull,w431,14757395255531667466ull,w431,w431,w431,
		w431,w431,w431, 3435973836ull,w431,14757395255547920448ull, 59954448353250508ull,w431,14757395255543332928ull, 47569549377981644ull,
		w431,14757395255535140928ull, 35261719295872204ull,w431,14757395255539597416ull, 38562521921932492ull,w431,14757395255540157512ull,
		7094523705632279756ull,14757395256418681544ull, 3435973836ull,w431,w431, 1863156813004ull,w431,14757395255531669504ull,
		243954142412ull,w431,14757395255531667470ull, 304083684556ull,w431,14757395255531667470ull, 46385646796ull,w431,14757395255531667514ull,
		5347289241946570ull, 817021423881133759ull, 9222773902559490395ull,w431,w431, 7730941132ull, 0ull,w431,w431,w431,w431,w431,w431,w431,
		w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431, 7730941132ull,w431,w431,w431,w431,w431,w431,
		w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,
		w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,
		w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,
		w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,
		w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,
		w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,
		w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,
		w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,
		w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,
		w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,
		w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,
		w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,
		w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431,w431
	};

	/* Note the revsre digit order here w.r.to the factor.c testFac* inits - that's because for the
	latter we declare an explicit digit ordering, here we have to accept the default, which is x[0], x[1],...
	*/
	const uint128 testp128[] =
	{
		{ 1654746039858251761ull,   12240518780192025ull},
		{  353773459776294223ull, 2923447923687422893ull},
		{ 7733695761441692271ull,     128099917305337ull},
		{ 5349936413307099433ull,   10811609908058563ull},
		{  701890237964104231ull,        700245770430ull},
		{ 7804835620876695225ull,  756146046438660814ull},
		{17050154159176967743ull,  106450884062962221ull},
		{15223106393317212577ull, 7448657723978021346ull},
		{16055621295463638505ull,     644719741813452ull},
		{ 4961431981940124743ull,      44696312570505ull},
		{ 1738540175825943815ull,   99970972632587991ull},
		{13887741953162944095ull,         67677680549ull},
		{ 2894679571106043497ull,    5287390011750720ull},
		{10375766809019373543ull,      11129117045170ull},
		{ 8321741389535251703ull,       1551337752834ull},
		{ 6586095673132787791ull,  133834206206032981ull},
		{ 2460710484528304153ull,       5747037125100ull},
		{ 7361144750966677159ull,      10824073357153ull},
		{ 8212830436061989903ull,   32559650929209964ull},
		{ 4235226679561594903ull,   94004235929829273ull},
		{ 8676403300852410079ull,        103337218078ull},
		{14211535226588354713ull,      62897895526806ull},
		{14854696485656401105ull,          5492917609ull},
		{10285356664749312993ull,          5799457823ull},
		{16578386512849109713ull,          4303087381ull},
		{18263664019678288919ull,          5202063708ull},
		{ 7607475409566672241ull,          7785579841ull},
		{ 5630449305759171207ull,          7593776864ull},
		{ 9058738039473012457ull,         20308449831ull},
		{ 5034609389988515233ull,         15531431134ull},
		{ 8291172543411688687ull,         18216394609ull},
		{15859685870849762975ull,         16259503442ull},
		{11995354231649723881ull,         20551559047ull},
		{15122069645900159367ull,         28364424832ull},
		{ 9636073161914837921ull,         34441477586ull},
		{ 1304857345634219175ull,         30977655046ull},
		{ 4788416424359163737ull,         43144178324ull},
		{ 2258783450670948535ull,         45963786472ull},
		{15262466751214122975ull,         66325700032ull},
		{ 3082735332820781609ull,         57210216387ull},
		{15689518004012743009ull,         80355238912ull},
		{10819675441336938065ull,        109346652057ull},
		{10329047311584913071ull,      17534809723250ull},
		{ 2595613477835803991ull,       1221873279710ull},
		{ 8612165677489771129ull,      12549422209078ull},
		{ 9015544550402598895ull,        112416184026ull},
		{11385509628023387489ull,        142220976614ull},
		{ 2773697912020767049ull,      14320762091913ull},
		{ 7515490546312285159ull,       1996508583829ull},
		{ 2388855518206098663ull,        365842230851ull},
		{  465403687781705377ull,        261078947686ull},
		{17079803649869853575ull,     199835753775288ull},
		{15105677628487752455ull,        202355339943ull},
		{16692905930976531153ull,      18738454648009ull},
		{18170931828058363183ull,        412571049040ull},
		{ 2216600112648316881ull,     534505286298455ull},
		{0ull,0ull}
	};

	const uint192 testp192[] =
	{
		{15875370168207932041ull,11545660419510266595ull,                  133ull},
		{  509892144742137431ull,15571349859840161706ull,                 1394ull},
		{14226674137430228263ull, 4492854135134704005ull,               121320ull},
		{10174116463236461383ull,14842833464112563611ull,           2649519282ull},
		{ 7660429456444636239ull,17652352551621896287ull,            655903171ull},
		{ 2219421057460140527ull,18314245293386716597ull,           1083827012ull},
		{12343089078196252631ull,18225246095436784582ull,             13161208ull},
		{ 8097149896429635207ull,14663183769241509326ull,                 4730ull},
		{17184239148975426263ull,  881920578744577810ull,               215159ull},
		{17733134473107607967ull, 9900144438119899815ull,            212724356ull},
		{ 2803405107698253561ull, 5238930328752646394ull,                  261ull},
		{16346425147370540471ull, 4415476118538293365ull,                    1ull},
		{ 6167785434693019223ull,11905462972019801043ull,                70130ull},
		{17951008765075981215ull,18429773635221665090ull,              5800574ull},
		{15903397166638806257ull,14500669099417213747ull,             22381525ull},
		{ 3893270457587058239ull, 3291757557782450881ull,                   14ull},
		{14288981644299514807ull, 1390029428449091172ull,                 1552ull},
		{ 5085420234315110585ull,14802171160149427175ull,                 2674ull},
		{ 4949688733053552967ull,14291576310931480037ull,                  664ull},
		{11405337619840706193ull, 6334326874596939334ull,               617742ull},
		{  329809049266961143ull,10558642444782195772ull,      157590042578912ull},
		{17814616685598394119ull, 1933308633079010416ull,        9118322195022ull},
		{ 3547755741880899889ull,17012949627558354271ull,       70286054459973ull},
		{16007877010440112335ull, 8040689323464953445ull,   492416983078691417ull},
		{10050950882119470361ull, 9565712986615012496ull,          59364131986ull},
		{0ull,0ull,0ull}
	};

	char char_buf[STR_MAX_LEN];
	double dbl;
	const uint64 *ptr_mword_prime[10], *curr_mword_prime;
	uint64 *out_array = 0x0;
	uint64 ONES64 = 0xFFFFFFFFFFFFFFFFull;	// In GCC, making this 'const' gives "warning: overflow in implicit constant conversion" wherever it is used.
uint64 k;
uint192 x192;
	/*...time-related stuff	*/
	clock_t clock1, clock2;
	double tdiff;

	ptr_mword_prime[0] = y1;
	ptr_mword_prime[1] = y2;
	ptr_mword_prime[2] = y3;
	ptr_mword_prime[3] = y4;
	ptr_mword_prime[4] = y5;
	ptr_mword_prime[5] = y6;
	ptr_mword_prime[6] = y7;
	ptr_mword_prime[7] = y8;
	ptr_mword_prime[8] = y9;
	ptr_mword_prime[9] = y10;

#ifdef GCD_STANDALONE
	gcd_init();
#endif

	/* Init the RNG: */
	rng_isaac_init(TRUE);

	/* Allocate the main data arrays: */
	u = (uint64 *)calloc(  MAX_ARRAY_DIM, sizeof(uint64));
	v = (uint64 *)calloc(  MAX_ARRAY_DIM, sizeof(uint64));
	w = (uint64 *)calloc(2*MAX_ARRAY_DIM, sizeof(uint64));
	x = (uint64 *)calloc(  MAX_ARRAY_DIM, sizeof(uint64));
	y = (uint64 *)calloc(  MAX_ARRAY_DIM, sizeof(uint64));
	z = (uint64 *)calloc(2*MAX_ARRAY_DIM, sizeof(uint64));
	a = (double *)calloc(8*MAX_ARRAY_DIM, sizeof(double));
	b = (double *)calloc(8*MAX_ARRAY_DIM, sizeof(double));
	c = (double *)calloc(8*MAX_ARRAY_DIM, sizeof(double));

	// Basic tests of mi64_add/sub:
	ASSERT(HERE, 1000 < MAX_ARRAY_DIM-2, "MAX_ARRAY_DIM too small!");
	for(i = 0; i < 1000; ++i)
	{
		for(j = 0; j <= i; ++j)
		{
			x[j] = u[j] = rng_isaac_rand();
			y[j] = v[j] = rng_isaac_rand();
		}
		z[j] = mi64_add_ref(x,y,z, i+1);
		w[j] = mi64_add    (u,v,w, i+1);
		if(!mi64_cmp_eq(w,z, i+2)) {
			for(j = 0; j <= i+1; ++j) {
				if(w[j] != z[j]) {
					printf("w[%3d] = %20llu != z[%3d] = %20llu\n",j,w[j],j,z[j]);
				}
			}
			ASSERT(HERE, 0, "mi64_add results differ!");
		}
	}

	// mi64_add timing test
	vec_len = 1000;
	for(j = 0; j <= i; ++j)
	{
		x[j] = u[j] = rng_isaac_rand();
		y[j] = v[j] = rng_isaac_rand();
	}
	clock1 = clock();
	for(i = 0; i < 1000000; i++) {
		ASSERT(HERE, mi64_add(x,y,z, vec_len) < 2, "GCD self-test mi64_add!");
	}
	clock2 = clock();
	tdiff = (double)(clock2 - clock1);
	printf	("GCD: mi64_add() self-test passed: %u*%u words processed, Time =%s\n",i,vec_len,get_time_str(tdiff));

// 8/08/2012: Track down mi64_div bug which leaves remainder with 1*divisor needing subtraction:
	// x= 184027369222821018387581692039:
	x[1] = 9976143675ull;	x[0] =  7539741246537263239ull;
	y[1] =         10ull;	y[0] = 10733642766940458262ull;
	mi64_div(x,y,2,2,u,v);	// quo in u[], rem in v[]
	ASSERT(HERE, u[1]==0 && u[0]==942757929 && v[1]==0 && v[0]==1, "divisor check fails!");

// 6/08/2012: Trial-factor of mm127 residue check for Paul L:
#if 0
	p = 127;	lenX = (p>>6);
	memset(x,ONES64,(lenX<<3));	x[lenX++] = (1ull << (p&63)) - 1;	// x has p
	k = 20;	// Really testing k = 10 ... absorb the 2* part of q = 2.k.p+1 into this constant
	y[lenX] = mi64_mul_scalar(x, k, y, lenX); ++y[0];	++lenX;		// y = q; Can add 1 sans carry check since 2.k.p even
	x192 = twopmodq192(*(uint192*)x,*(uint192*)y);
	printf("Res127 = %s\n", &string0[convert_mi64_base10_char(string0, (uint64*)&x192, lenX, 0)]);
exit(0);
#endif

/**** 7/24/2012: Collect stats on the 64-bit mod-inverses: Results summarized in http://www.mersenneforum.org/showpost.php?p=305815&postcount=35 ****/
#if 0
	for(j = 0; (q = testp128[j].d0) != 0; ++j)
	{
		q = (q & 0xFFFFFFFFFFFFFFF0) + 0xf;	// Control q (mod 16) here; 1/3/5/7/9/b/d/f give >= 6/5/6/5/5/6/5/6 good bits in qinv_0, resp.
		qinv = (q+q+q) ^ (uint64)2;
		printf("q = %16llX, qinv_0 = %16llX\n",q,qinv);
		for(i = 0; i < 4; i++)
		{
			tmp = q*qinv;
			qinv = qinv*((uint64)2 - tmp);
			printf("i = %1u, tmp = %16llX, qinv_%1u = %16llX\n",i+1,tmp,i+1,qinv);
		}
		printf("\n");
	}
	exit(0);
#elif 0
Outputs:
q = 16F6D6C18B3C47F1, qinv_0 = 44E48444A1B4D7D1	i = 1, tmp = EE480F97731622C1, qinv_1 = 124D858CE6733911	i = 2, tmp = E4139C5D82487001, qinv_2 = BC53762ACEB3C911	i = 3, tmp = 2C90F9C0CF000001, qinv_3 = 3437C0D60FB3C911	i = 4, tmp = D89F000000000001, qinv_4 = FAA8C0D60FB3C911
q =  4E8DB1658ED4D4F, qinv_0 =  EBA91430AC7E7EF	i = 1, tmp =  D39144454B675C1, qinv_1 = 36C7715847EFB9AF	i = 2, tmp = 712332AE5CD6F001, qinv_2 = 69BD01D56D91A9AF	i = 3, tmp = 14668F09DF000001, qinv_3 = FFA966DEFC91A9AF	i = 4, tmp = 8FBF000000000001, qinv_4 = A51866DEFC91A9AF
q = 6B539B51F5AD1A6F, qinv_0 = 41FAD1F5E1074F4F	i = 1, tmp = A48B004D839C6941, qinv_1 = C414C264DE88148F	i = 2, tmp = 44BF3C380EBA7001, qinv_2 =  84C7DC293A3848F	i = 3, tmp =  91613F90F000001, qinv_3 = 79795CE732A3848F	i = 4, tmp = D11F000000000001, qinv_4 = AD285CE732A3848F
q = 4A3ECAD29E18D929, qinv_0 = DEBC6077DA4A8B79	i = 1, tmp =  269678EE781E761, qinv_1 = 6095BE0443830F19	i = 2, tmp = 85986005C0219C01, qinv_2 =  48B7A2BB416D319	i = 3, tmp = 797D799668F00001, qinv_3 = 870DDD5DA4A6D319	i = 4, tmp = FBBC1F0000000001, qinv_4 = CA21D65DA4A6D319
q =  9BD9D73E12A3227, qinv_0 = 1D38D85BA37E9677	i = 1, tmp = F0A5835CE8322A21, qinv_1 = C15E253B33BE4197	i = 2, tmp = A4FA02D184917C01, qinv_2 =  5D2554F09721D97	i = 3, tmp = CFDD2D724BF00001, qinv_3 = FB95536A0EE21D97	i = 4, tmp = E7B97F0000000001, qinv_4 = 2BC86A6A0EE21D97
q = 6C5058A538ABB6B9, qinv_0 = 44F109EFAA032429	i = 1, tmp = 80CC377F375D47A1, qinv_1 = E5F22CBE97202B89	i = 2, tmp = C9B2F87C21B5DC01, qinv_2 = 5F3E952780D96F89	i = 3, tmp =  9227C172AF00001, qinv_3 = 8850BBB376696F89	i = 4, tmp = C1AC5F0000000001, qinv_4 = 2EE0E4B376696F89
q = EC9E50DF4764663F, qinv_0 = C5DAF29DD62D32BF	i = 1, tmp = 90B7AF4982F39701, qinv_1 = 7F0D8F2603F189BF	i = 2, tmp = 23A594DBFCEF0001, qinv_2 = 864FDAA966A089BF	i = 3, tmp = ED2E98DF00000001, qinv_3 = 3B5C754866A089BF	i = 4, tmp =                1, qinv_4 = 3B5C754866A089BF
q = D343569BF7A681A1, qinv_0 = 79CA03D3E6F384E1	i = 1, tmp = 312EF0FD3601F281, qinv_1 = 69F342363EB36261	i = 2, tmp = 6811DCAE3549C001, qinv_2 = 24EE54BC9241A261	i = 3, tmp = 42A76140F0000001, qinv_3 = FD9A7341A241A261	i = 4, tmp = 1F00000000000001, qinv_4 = 3E9A7341A241A261
q = DED106807C4EA5E9, qinv_0 = 9C73138174EBF1B9	i = 1, tmp = FA0F4CA8BFE93E61, qinv_1 = 211267246F857E59	i = 2, tmp = 41B7426A260D5C01, qinv_2 = 5C2191ADD5988259	i = 3, tmp =  BD7D7FD86F00001, qinv_3 = 526AED840C288259	i = 4, tmp = E27FDF0000000001, qinv_4 = CCB866840C288259
q = 44DA8C19CCCDAC47, qinv_0 = CE8FA44D666904D7	i = 1, tmp = B5CE51392A8BCBA1, qinv_1 = F4BFA30837328177	i = 2, tmp = 9A51E798BE48DC01, qinv_2 = 6BE751A614783D77	i = 3, tmp = C70284B382F00001, qinv_3 = 3E1CBB8106E83D77	i = 4, tmp = 21675F0000000001, qinv_4 = FD6C928106E83D77
q = 1820891427D9BD07, qinv_0 = 48619B3C778D3717	i = 1, tmp = CB2259CEBA077CA1, qinv_1 = BCD4192420FAA4B7	i = 2, tmp =  59D4D2F72949C01, qinv_2 = D5251E31ACCF20B7	i = 3, tmp = FC4E06CB40F00001, qinv_3 = 1C9A18B8413F20B7	i = 4, tmp = B2E71F0000000001, qinv_4 = B482EFB8413F20B7
q = C0BB2BB9DA911E5F, qinv_0 = 4231832D8FB35B1F	i = 1, tmp = 3F7C33555ACB7281, qinv_1 = DDA732DAAFDCFD9F	i = 2, tmp = 991348CF35C9C001, qinv_2 =  E938E5EE4CEBD9F	i = 3, tmp = FC6E5780F0000001, qinv_3 = 6FEFE419D4CEBD9F	i = 4, tmp = 1F00000000000001, qinv_4 = 2EEFE419D4CEBD9F
q = 282BF7BBB65B1669, qinv_0 = 7883E73323114339	i = 1, tmp = C61F0F8D491E7861, qinv_1 = 3855D5FFC36755D9	i = 2, tmp = 2631C1F8D0E5DC01, qinv_2 = 5F8385AADA83D9D9	i = 3, tmp = CA776C1CAAF00001, qinv_3 = 9608BBA88513D9D9	i = 4, tmp = A95C5F0000000001, qinv_4 = AD3534A88513D9D9
q = 8FFE20F08BE357E7, qinv_0 = AFFA62D1A3AA07B7	i = 1, tmp = E735F74EAB512721, qinv_1 = 99D5D1F899552FD7	i = 2, tmp = 9B26E48177C53C01, qinv_2 = EF985EE6E0ABCBD7	i = 3, tmp = F050924299F00001, qinv_3 = D123E596481BCBD7	i = 4, tmp = 83AF3F0000000001, qinv_4 = 9D00FC96481BCBD7
q = 737CC3B40BEA0CF7, qinv_0 = 5A764B1C23BE26E7	i = 1, tmp = 4722F0BAE2705CE1, qinv_1 = 83CA2045945118C7	i = 2, tmp = 92CA75331A4E3C01, qinv_2 = 967581A3CFE074C7	i = 3, tmp = 4F481BE761F00001, qinv_3 = 7A4C3665EE5074C7	i = 4, tmp = 7B683F0000000001, qinv_4 = 9FB73D65EE5074C7
q = 5B66831EBDD1284F, qinv_0 = 1233895C397378EF	i = 1, tmp = DCB76CFFF3A6A9C1, qinv_1 =  1A57B87A008FEAF	i = 2, tmp = 625A00BCEA70F001, qinv_2 = CC17103FCEB4EEAF	i = 3, tmp = 136CD76D1F000001, qinv_3 = 6BBC8DD59DB4EEAF	i = 4, tmp = 963F000000000001, qinv_4 = 24AB8DD59DB4EEAF
q = 22263325F5639019, qinv_0 = 66729971E02AB049	i = 1, tmp = 8E339AFBF08F4721, qinv_1 = 425AD75A1B296829	i = 2, tmp = 2CAB11B3767D3C01, qinv_2 = EB1BD216A4BACC29	i = 3, tmp = 4C0E676C59F00001, qinv_3 = 196BD1B0FD4ACC29	i = 4, tmp = DDE73F0000000001, qinv_4 = 132EBAB0FD4ACC29
q = 66280A2665999EA7, qinv_0 = 32781E7330CCDBF7	i = 1, tmp = BBD83F252F04F021, qinv_1 = 2A005F04401DCD17	i = 2, tmp = FBC0067BDDC3FC01, qinv_2 = E1818B19D6B62917	i = 3, tmp = F5BFE0DE1FF00001, qinv_3 = 17336F6788262917	i = 4, tmp = 97C3FF0000000001, qinv_4 = 36C0866788262917
q = 71F9D5C4A5FDCC0F, qinv_0 = 55ED814DF1F9642F	i = 1, tmp = 46F5606A36E552C1, qinv_1 = 42CD64A80CCC32EF	i = 2, tmp = 1CC1ECA9D9C07001, qinv_2 = 39F03AF4EC43A2EF	i = 3, tmp = 18345177CF000001, qinv_3 = 745B2E1CAB43A2EF	i = 4, tmp = E69F000000000001, qinv_4 = 87EA2E1CAB43A2EF
q = 3AC68C39D286A817, qinv_0 = B053A4AD7793F847	i = 1, tmp = E31CB6FC0C63E661, qinv_1 = 9562FC71579213A7	i = 2, tmp = EA286E6B026F5C01, qinv_2 = 943E2BAE05190FA7	i = 3, tmp = 83D32A1F16F00001, qinv_3 = 9669F19DFE890FA7	i = 4, tmp = AFD1DF0000000001, qinv_4 = 4170789DFE890FA7
q = 7868C6D52355A6DF, qinv_0 = 693A547F6A00F49F	i = 1, tmp = CEFD0FC7EB3F3081, qinv_1 = F257C1871387D51F	i = 2, tmp = 26E6865C17CFC001, qinv_2 = 4FFFA2C2169F951F	i = 3, tmp = 637702E7F0000001, qinv_3 = 633739FC069F951F	i = 4, tmp = FF00000000000001, qinv_4 = 823739FC069F951F
q = C53983FE1DBFD899, qinv_0 = 4FAC8BFA593F89C9	i = 1, tmp = 6B6D3D5E8531F121, qinv_1 = 784CCA1D941F17A9	i = 2, tmp =  7DCBE368EA2BC01, qinv_2 = 7EA3AAE407CCFBA9	i = 3, tmp = D4CFBFFD85F00001, qinv_3 = B9D00EF44C5CFBA9	i = 4, tmp = DD8CBF0000000001, qinv_4 = D4A0F7F44C5CFBA9
q = CE267BC009CB14D1, qinv_0 = 6A7373401D613E71	i = 1, tmp = 86D6F0B211DFCE41, qinv_1 = 4F6EA4279F63B431	i = 2, tmp = ADCDB1F28654F001, qinv_2 = F707CA5F6061C431	i = 3, tmp = ED0D38919F000001, qinv_3 = 6C6339C3F161C431	i = 4, tmp = 7F3F000000000001, qinv_4 = D55439C3F161C431
q = 8EBCED6476944BE1, qinv_0 = AC36C82D63BCE3A1	i = 1, tmp = 5BC2EE20CEC83B81, qinv_1 = 38BA2EA1990CF821	i = 2, tmp = 2A4843A0FA2BC001, qinv_2 = B0D1E4B9F7693821	i = 3, tmp = 6E568585F0000001, qinv_3 = 1BF5F1F607693821	i = 4, tmp = BF00000000000001, qinv_4 = 7CF5F1F607693821
q = E61242B2877F82D1, qinv_0 = B236C817967E8871	i = 1, tmp = 657FB8F965A5C641, qinv_1 = E7FBAA8093800631	i = 2, tmp = 9AAC1DB857F8F001, qinv_2 = F1CFB9DFE73A1631	i = 3, tmp = 858BDACE1F000001, qinv_3 = CB911BC1F83A1631	i = 4, tmp = 183F000000000001, qinv_4 = BD821BC1F83A1631
q = FD7591A92E01DC17, qinv_0 = F860B4FB8A059447	i = 1, tmp = 8A1F11BF88345661, qinv_1 = 5F49FE6911B21FA7	i = 2, tmp = 1AB6B646CBDB5C01, qinv_2 = E30277E88D751BA7	i = 3, tmp = 2AE4102176F00001, qinv_3 = B33ADED8A6E51BA7	i = 4, tmp = 1EDDDF0000000001, qinv_4 = AEF965D8A6E51BA7
q = 69932E9430BDA571, qinv_0 = 3CB98BBC9238F051	i = 1, tmp = E6365667D9D348C1, qinv_1 = 112B05F5282AEB91	i = 2, tmp = 58C295C97ED37001, qinv_2 = CB0A2F37DA987B91	i = 3, tmp = B8C54F1E2F000001, qinv_3 = 8F52F68A3B987B91	i = 4, tmp = F35F000000000001, qinv_4 = 1183F68A3B987B91
q = 4E235FCFEB9EC287, qinv_0 = EA6A1F6FC2DC4797	i = 1, tmp = 22F7AB5F7C9C2EA1, qinv_1 = 77B185A0A04E6737	i = 2, tmp = 2863F8CC24821C01, qinv_2 = 663AFDEFE7166337	i = 3, tmp =  D2E85FF8CF00001, qinv_3 = B337B3E7CF866337	i = 4, tmp = 8C489F0000000001, qinv_4 = 40208AE7CF866337
q = 7DB71A406C1B7EE9, qinv_0 = 79254EC144527CB9	i = 1, tmp = 6975C718FBF99261, qinv_1 = 4408B873E7513559	i = 2, tmp = 3BDFF5A9AD8E5C01, qinv_2 = 15F17C78FBC73959	i = 3, tmp = A048D07DCEF00001, qinv_3 = 7A03D1189A573959	i = 4, tmp = 4258DF0000000001, qinv_4 = D6774A189A573959
q = 45DE8690D11C8DA1, qinv_0 = D19B93B27355A8E1	i = 1, tmp = 29917D8EE07F2281, qinv_1 = 72575215C0745661	i = 2, tmp = F165FE7CC059C001, qinv_2 =  F8CC387B7F29661	i = 3, tmp = 7F874088F0000001, qinv_3 =  87A3304C7F29661	i = 4, tmp = 1F00000000000001, qinv_4 = 497A3304C7F29661
q = 7310297FED8694EF, qinv_0 = 59307C7FC893BECF	i = 1, tmp = D89CAFEE2F98CF41, qinv_1 =  FC4EC2AE072AA0F	i = 2, tmp = F5ABADAFBC377001, qinv_2 = D014593228D31A0F	i = 3, tmp = 9473677EAF000001, qinv_3 = 1E0C2DFFE7D31A0F	i = 4, tmp = 445F000000000001, qinv_4 = 767B2DFFE7D31A0F
q = DC18EC457548269F, qinv_0 = 944AC4D05FD873DF	i = 1, tmp = AB20DC962A5B1181, qinv_1 = B4C368A54DA7B55F	i = 2, tmp = 70317BDC8DCDC001, qinv_2 = B6AF5D91F58D755F	i = 3, tmp = C151B522F0000001, qinv_3 = 7DA5FDEAE58D755F	i = 4, tmp = 5F00000000000001, qinv_4 = 3CA5FDEAE58D755F
q = A6780F1E9C75E5E9, qinv_0 = F3682D5BD561B1B9	i = 1, tmp = 6AF927168F723E61, qinv_1 = 42DAE090AC523E59	i = 2, tmp = 623B32D5234D5C01, qinv_2 = 65F58D5332254259	i = 3, tmp = 9C76A97786F00001, qinv_3 = A292714B68B54259	i = 4, tmp = 65BFDF0000000001, qinv_4 = 205FEA4B68B54259
q = D1DC623908989987, qinv_0 = 759526AB19C9CC97	i = 1, tmp = 8B2DA079665922A1, qinv_1 = 47ACD042E3FEE037	i = 2, tmp = BF4A1E1F68111C01, qinv_2 =  CCB37A6C7D1DC37	i = 3, tmp = 6EFD501B44F00001, qinv_3 = A70DBA9CB841DC37	i = 4, tmp = 62CF9F0000000001, qinv_4 = A0CE919CB841DC37
q = 85BA357C959BEBA1, qinv_0 = 912EA075C0D3C2E1	i = 1, tmp = 868AFCB09C4D1A81, qinv_1 = 5FF0EEDD1B7A7861	i = 2, tmp = 5861807B0C41C001, qinv_2 = 44940BE92490B861	i = 3, tmp = F83B451CF0000001, qinv_3 = 3A6B0F723490B861	i = 4, tmp = 9F00000000000001, qinv_4 = FB6B0F723490B861
q = 121BC8D7A91634A7, qinv_0 = 36535A86FB429DF7	i = 1, tmp = 59A4AD8A12C53821, qinv_1 = 80AB672C134DD717	i = 2, tmp =  B23E9AC1271FC01, qinv_2 = ED16933BE46C3317	i = 3, tmp = D0C925CF8FF00001, qinv_3 = 4C695AA925DC3317	i = 4, tmp = 88F1FF0000000001, qinv_4 = A4DE71A925DC3317
q = 4273DF5E2A92F759, qinv_0 = C75B9E1A7FB8E609	i = 1, tmp = D6C803F0785CA821, qinv_1 = 59B021424E6A3CE9	i = 2, tmp = BBBD04D8BA95FC01, qinv_2 = 72C5FF62FCD7E0E9	i = 3, tmp = D30875F0AFF00001, qinv_3 = 94AC0AD0DB67E0E9	i = 4, tmp = 4515FF0000000001, qinv_4 = DB86F3D0DB67E0E9
q = 1F58CF94B1C5E0B7, qinv_0 = 5E0A6EBE1551A227	i = 1, tmp = CB169D036A4009E1, qinv_1 = 25EAE6A9DF706107	i = 2, tmp = D7C9EC8D8F9E7C01, qinv_2 = 63B20C05A51EFD07	i = 3, tmp = 8E38F75AB3F00001, qinv_3 = FD63DE96898EFD07	i = 4, tmp = DCC67F0000000001, qinv_4 = 56736596898EFD07
q = D3CF2CA56E53F7DF, qinv_0 = 7B6D85F04AFBE79F	i = 1, tmp = D47FAF57C8762C81, qinv_1 = 377B7BCD85EEC41F	i = 2, tmp = 6ADDE4EAF243C001, qinv_2 = B38DB421D0BA841F	i = 3, tmp =  8E0D711F0000001, qinv_3 = EB660935C0BA841F	i = 4, tmp = 3F00000000000001, qinv_4 = 4A660935C0BA841F
q = 2AC81373C1540A29, qinv_0 = 80583A5B43FC1E79	i = 1, tmp = DDA8F653F5459B61, qinv_1 = 1CD5529F1B606E19	i = 2, tmp = 9542716D1FF29C01, qinv_2 = 97A9E5CE3DA73219	i = 3, tmp = 29B2804CB0F00001, qinv_3 = 9EDBC73216373219	i = 4, tmp = 77351F0000000001, qinv_4 = F19DC03216373219
q = D9BC5D7F7E467961, qinv_0 = 8D35187E7AD36C21	i = 1, tmp = AB225661473D9181, qinv_1 = FDC43FD06601AAA1	i = 2, tmp = EBC2877C564DC001, qinv_2 = 37292880BD9BEAA1	i = 3, tmp = 4715AB62F0000001, qinv_3 = E3C89FE7CD9BEAA1	i = 4, tmp = 5F00000000000001, qinv_4 = 24C89FE7CD9BEAA1
q = 9627357D21E7EE51, qinv_0 = C275A07765B7CAF1	i = 1, tmp = EAB93D063E4A4441, qinv_1 = 7974B2AC1E330AB1	i = 2, tmp = FC563E6578CDF001, qinv_2 = 35C724A4BC701AB1	i = 3, tmp = 743EE955BF000001, qinv_3 =  8CD8AF5AD701AB1	i = 4, tmp = 9B7F000000000001, qinv_4 = 9FFE8AF5AD701AB1
q = 8F5825CDE338E6AF, qinv_0 = AE087169A9AAB40F	i = 1, tmp = DE1A4BAE09BE9041, qinv_1 = E870F3570F93404F	i = 2, tmp = CB0F75010FB7F001, qinv_2 = 3191366669D0304F	i = 3, tmp = 9E390DB6FF000001, qinv_3 = DE1D7B1DB8D0304F	i = 4, tmp = 6DFF000000000001, qinv_4 = 1C6C7B1DB8D0304F
q = 240578B4B89E8D57, qinv_0 = 6C106A1E29DBA807	i = 1, tmp = 9A12D090F183F561, qinv_1 = 72F05839AD18F267	i = 2, tmp = F7763BFC348F1C01, qinv_2 = E902ECBA9F0CAE67	i = 3, tmp = B48E6C9FB4F00001, qinv_3 = 76A3203DB27CAE67	i = 4, tmp = 3DFD9F0000000001, qinv_4 = 1F86273DB27CAE67
q = 77848F0DF196D279, qinv_0 = 668DAD29D4C47769	i = 1, tmp = 76A4802B8B5692A1, qinv_1 = C0B2E233AD99F3C9	i = 2, tmp = F5739DBF68851C01, qinv_2 = DBEB573A8782F7C9	i = 3, tmp = 725BCC09E4F00001, qinv_3 = 3CD9C5B23712F7C9	i = 4, tmp = 9A639F0000000001, qinv_4 = B838EEB23712F7C9
q = 7D1DA6008F74A7EF, qinv_0 = 7758F201AE5DF7CF	i = 1, tmp = 1F35C252012E6341, qinv_1 = 4ECAD9FF69D8F70F	i = 2, tmp = 575A05D1AE857001, qinv_2 = 188CA7CCF0F7670F	i = 3, tmp = 9CC576326F000001, qinv_3 = AFB3572F6FF7670F	i = 4, tmp = 73DF000000000001, qinv_4 = 2CA2572F6FF7670F
q = 9E0174AF051FB561, qinv_0 = DA045E0D0F5F2021	i = 1, tmp = 71BF81859DC18181, qinv_1 = 44EF8335FABD6EA1	i = 2, tmp = 723585E47B7DC001, qinv_2 = 4D2D681D8827AEA1	i = 3, tmp = A2EFEBBAF0000001, qinv_3 = 5916886C9827AEA1	i = 4, tmp = 5F00000000000001, qinv_4 = 9A16886C9827AEA1
q = 267E278E283D9149, qinv_0 = 737A76AA78B8B3D9	i = 1, tmp =  945D677CD3E31E1, qinv_1 = 18F6DBC62920CCF9	i = 2, tmp = ED7E76C30DC87C01, qinv_2 = E9972C267E5030F9	i = 3, tmp = 1D851E6603F00001, qinv_3 = B1F81F2FA9E030F9	i = 4, tmp = 9CB07F0000000001, qinv_4 = 127C982FA9E030F9
q = 684C62D545C00BE7, qinv_0 = 38E5287FD14023B7	i = 1, tmp = 9D260851AEA91721, qinv_1 = 2A6EDC33F0373BD7	i = 2, tmp = 52E312A0F5A93C01, qinv_2 = 2944BACBBA41D7D7	i = 3, tmp = CC47A347B9F00001, qinv_3 =  B5AC67601B1D7D7	i = 4, tmp = B5D33F0000000001, qinv_4 = 5E07DD7601B1D7D7
q = 2126EB6FE667FCE7, qinv_0 = 6374C24FB337F6B7	i = 1, tmp = C19F23DB7BFCC321, qinv_1 = 490560A7E227BAD7	i = 2, tmp = E2628A2E84463C01, qinv_2 = A5B03BD09F9356D7	i = 3, tmp = 89B9EEDB21F00001, qinv_3 = 5F6277907F0356D7	i = 4, tmp = AAE03F0000000001, qinv_4 = CAE38E907F0356D7
q =  67572302F62A6A1, qinv_0 = 136056908E27F3E1	i = 1, tmp = 8314DA58FA664681, qinv_1 = 4F206FFDA7D87D61	i = 2, tmp = B0B19D71BE95C001, qinv_2 = 4E67FF3B525ABD61	i = 3, tmp =  469DF66F0000001, qinv_3 = A54AFA0A625ABD61	i = 4, tmp = DF00000000000001, qinv_4 = 264AFA0A625ABD61
q = ED07A6ED47E9B387, qinv_0 = C716F4C7D7BD1A97	i = 1, tmp = 1F6C268412BF9AA1, qinv_1 =  B105E47CAE4A637	i = 2, tmp = 902189D967DB1C01, qinv_2 =   D89FFD76A9A237	i = 3, tmp =  8707BEF14F00001, qinv_3 = 3C9541F01719A237	i = 4, tmp = 46299F0000000001, qinv_4 = 4C0618F01719A237
q = D1A225C51193DB07, qinv_0 = 74E6714F34BB9117	i = 1, tmp = F0149150BD74A4A1, qinv_1 = 54DF24BA5D2226B7	i = 2, tmp = 2F666BC625229C01, qinv_2 = 35C4D74DF63CA2B7	i = 3, tmp = 4AABAA3A30F00001, qinv_3 = 1A7B3C7D1AACA2B7	i = 4, tmp =  9E51F0000000001, qinv_4 = 3614137D1AACA2B7
q = E7A91D805CDE12D1, qinv_0 = B6FB5881169A3871	i = 1, tmp = 3297607FD9DE0641, qinv_1 = 7FA47E1017BB7631	i = 2, tmp = 381BA558A8D8F001, qinv_2 = 3185263C77958631	i = 3, tmp = 361D8D2A1F000001, qinv_3 = 5C0C09F288958631	i = 4, tmp = D03F000000000001, qinv_4 = 85FD09F288958631
q = FC2C1E395142152F, qinv_0 = F4845AABF3C63F8F	i = 1, tmp = BB5E3CA0A67A6641, qinv_1 = C2F9357D00BD61CF	i = 2, tmp = 83CF38DF6228F001, qinv_2 = D544EC71EEB351CF	i = 3, tmp = CE0A81B41F000001, qinv_3 = 12F860FDDDB351CF	i = 4, tmp = 643F000000000001, qinv_4 = 150760FDDDB351CF
q = 1EC2F20EF37E67D1, qinv_0 = 5C48D62CDA7B3771	i = 1, tmp = 51E13EF97E84BA41, qinv_1 =  7361C7D81214131	i = 2, tmp = 748FC358667EF001, qinv_2 = EED1E93EB7E55131	i = 3, tmp = DE1199CEDF000001, qinv_3 = 6252891708E55131	i = 4, tmp = 59BF000000000001, qinv_4 = C5C3891708E55131

Now focus just on the final (4th) iteration, sort by #good bits in tmp:
q = 66280A2665999EA7	i = 4, tmp = 97C3FF0000000001, qinv_4 = 36C0866788262917
q = 4273DF5E2A92F759	i = 4, tmp = 4515FF0000000001, qinv_4 = DB86F3D0DB67E0E9
q = 121BC8D7A91634A7	i = 4, tmp = 88F1FF0000000001, qinv_4 = A4DE71A925DC3317
q = 3AC68C39D286A817	i = 4, tmp = AFD1DF0000000001, qinv_4 = 4170789DFE890FA7
q = 7DB71A406C1B7EE9	i = 4, tmp = 4258DF0000000001, qinv_4 = D6774A189A573959
q = A6780F1E9C75E5E9	i = 4, tmp = 65BFDF0000000001, qinv_4 = 205FEA4B68B54259
q = DED106807C4EA5E9	i = 4, tmp = E27FDF0000000001, qinv_4 = CCB866840C288259
q = FD7591A92E01DC17	i = 4, tmp = 1EDDDF0000000001, qinv_4 = AEF965D8A6E51BA7
q = C53983FE1DBFD899	i = 4, tmp = DD8CBF0000000001, qinv_4 = D4A0F7F44C5CFBA9
q = 4E235FCFEB9EC287	i = 4, tmp = 8C489F0000000001, qinv_4 = 40208AE7CF866337
q = D1DC623908989987	i = 4, tmp = 62CF9F0000000001, qinv_4 = A0CE919CB841DC37
q = ED07A6ED47E9B387	i = 4, tmp = 46299F0000000001, qinv_4 = 4C0618F01719A237
q = 77848F0DF196D279	i = 4, tmp = 9A639F0000000001, qinv_4 = B838EEB23712F7C9
q = 240578B4B89E8D57	i = 4, tmp = 3DFD9F0000000001, qinv_4 = 1F86273DB27CAE67
q = 1F58CF94B1C5E0B7	i = 4, tmp = DCC67F0000000001, qinv_4 = 56736596898EFD07
q =  9BD9D73E12A3227	i = 4, tmp = E7B97F0000000001, qinv_4 = 2BC86A6A0EE21D97
q = 267E278E283D9149	i = 4, tmp = 9CB07F0000000001, qinv_4 = 127C982FA9E030F9
q = 6C5058A538ABB6B9	i = 4, tmp = C1AC5F0000000001, qinv_4 = 2EE0E4B376696F89
q = 44DA8C19CCCDAC47	i = 4, tmp = 21675F0000000001, qinv_4 = FD6C928106E83D77
q = 282BF7BBB65B1669	i = 4, tmp = A95C5F0000000001, qinv_4 = AD3534A88513D9D9
q = 8FFE20F08BE357E7	i = 4, tmp = 83AF3F0000000001, qinv_4 = 9D00FC96481BCBD7
q = 22263325F5639019	i = 4, tmp = DDE73F0000000001, qinv_4 = 132EBAB0FD4ACC29
q = 737CC3B40BEA0CF7	i = 4, tmp = 7B683F0000000001, qinv_4 = 9FB73D65EE5074C7
q = 684C62D545C00BE7	i = 4, tmp = B5D33F0000000001, qinv_4 = 5E07DD7601B1D7D7
q = 2126EB6FE667FCE7	i = 4, tmp = AAE03F0000000001, qinv_4 = CAE38E907F0356D7
q = 1820891427D9BD07	i = 4, tmp = B2E71F0000000001, qinv_4 = B482EFB8413F20B7
q = D1A225C51193DB07	i = 4, tmp = 09E51F0000000001, qinv_4 = 3614137D1AACA2B7
q = 4A3ECAD29E18D929	i = 4, tmp = FBBC1F0000000001, qinv_4 = CA21D65DA4A6D319
q = 2AC81373C1540A29	i = 4, tmp = 77351F0000000001, qinv_4 = F19DC03216373219
q = 8F5825CDE338E6AF	i = 4, tmp = 6DFF000000000001, qinv_4 = 1C6C7B1DB8D0304F
q = 7D1DA6008F74A7EF	i = 4, tmp = 73DF000000000001, qinv_4 = 2CA2572F6FF7670F
q =  4E8DB1658ED4D4F	i = 4, tmp = 8FBF000000000001, qinv_4 = A51866DEFC91A9AF
q = 1EC2F20EF37E67D1	i = 4, tmp = 59BF000000000001, qinv_4 = C5C3891708E55131
q = 71F9D5C4A5FDCC0F	i = 4, tmp = E69F000000000001, qinv_4 = 87EA2E1CAB43A2EF
q = 16F6D6C18B3C47F1	i = 4, tmp = D89F000000000001, qinv_4 = FAA8C0D60FB3C911
q = 9627357D21E7EE51	i = 4, tmp = 9B7F000000000001, qinv_4 = 9FFE8AF5AD701AB1
q = 7310297FED8694EF	i = 4, tmp = 445F000000000001, qinv_4 = 767B2DFFE7D31A0F
q = 69932E9430BDA571	i = 4, tmp = F35F000000000001, qinv_4 = 1183F68A3B987B91
q = FC2C1E395142152F	i = 4, tmp = 643F000000000001, qinv_4 = 150760FDDDB351CF
q = E7A91D805CDE12D1	i = 4, tmp = D03F000000000001, qinv_4 = 85FD09F288958631
q = 5B66831EBDD1284F	i = 4, tmp = 963F000000000001, qinv_4 = 24AB8DD59DB4EEAF
q = CE267BC009CB14D1	i = 4, tmp = 7F3F000000000001, qinv_4 = D55439C3F161C431
q = E61242B2877F82D1	i = 4, tmp = 183F000000000001, qinv_4 = BD821BC1F83A1631
q = 6B539B51F5AD1A6F	i = 4, tmp = D11F000000000001, qinv_4 = AD285CE732A3848F
q = 7868C6D52355A6DF	i = 4, tmp = FF00000000000001, qinv_4 = 823739FC069F951F
q =  67572302F62A6A1	i = 4, tmp = DF00000000000001, qinv_4 = 264AFA0A625ABD61
q = 8EBCED6476944BE1	i = 4, tmp = BF00000000000001, qinv_4 = 7CF5F1F607693821
q = 85BA357C959BEBA1	i = 4, tmp = 9F00000000000001, qinv_4 = FB6B0F723490B861
q = DC18EC457548269F	i = 4, tmp = 5F00000000000001, qinv_4 = 3CA5FDEAE58D755F
q = 9E0174AF051FB561	i = 4, tmp = 5F00000000000001, qinv_4 = 9A16886C9827AEA1
q = D9BC5D7F7E467961	i = 4, tmp = 5F00000000000001, qinv_4 = 24C89FE7CD9BEAA1
q = D3CF2CA56E53F7DF	i = 4, tmp = 3F00000000000001, qinv_4 = 4A660935C0BA841F
q = D343569BF7A681A1	i = 4, tmp = 1F00000000000001, qinv_4 = 3E9A7341A241A261
q = C0BB2BB9DA911E5F	i = 4, tmp = 1F00000000000001, qinv_4 = 2EEFE419D4CEBD9F
q = 45DE8690D11C8DA1	i = 4, tmp = 1F00000000000001, qinv_4 = 497A3304C7F29661
q = EC9E50DF4764663F	i = 4, tmp = 0000000000000001, qinv_4 = 3B5C754866A089BF

So q == 7 or 9 (mod 16) gives 5 good bits in qinv_0, translating to 40 good bits in the q*qinv_3 product computed on the final iteration.
   q == 1 of F (mod 16) gives 6 or more good bits in qinv_0, with the odds of x good bits roughly halving with each increment in x.

But what about the "missing moduli", i.e. q == 3,5,B,D (mod 16)?

Tried those by twiddling the bottom 4 bits of the same dataset, all give >= 5 bits.
But interestingly, the mod-16 values which give "at least 5" (3,7,9,13) *never* seem to give more than 5,
whereas the ones which give "at least 6" (1,5,11,15) all give 6,7,8,... good bits with probability hakving with each additional "lucky" bit of precision.

#endif

// 07/31/2012: Test code using random inputs to test why "missing borrows" do not hose
// the loop-folded quotient algorithm (M(p) bad in this context since all the digits = B-1):
vec_len = 20;
for(i = 0; i < vec_len; i++)
{
	x[i]  = rng_isaac_rand();
}
mod = 16357897499336320049ull;
rem = mi64_div_by_scalar64(x, mod, vec_len, u);
ASSERT(HERE, rem == mi64_div_by_scalar64_u2(x, mod, vec_len, v) && mi64_cmp_eq(u,v,vec_len), "mi64_div_by_scalar64 Test #0.a!");
ASSERT(HERE, rem == mi64_div_by_scalar64_u4(x, mod, vec_len, v) && mi64_cmp_eq(u,v,vec_len), "mi64_div_by_scalar64 Test #0.b!");
exit(0);

/**** 5/16/2012: Tmp-code to test Montgomery-mul remainder for M977, see if can efficiently convert
carryout of my RL divisibility algo to true remainder:
*/
p = 977;	lenX = (p>>6);
memset(x,ONES64,(lenX<<3));	x[lenX++] = (1ull << (p&63)) - 1;
// Try an odd modulus close to 2^64:
mod = 16357897499336320049ull;
rem =  8623243291871090711ull;
ASSERT(HERE, !mi64_is_div_by_scalar64   (x, mod, lenX), "mi64_is_div_by_scalar64 Test #0.1!");
ASSERT(HERE, !mi64_is_div_by_scalar64_u2(x, mod, lenX), "mi64_is_div_by_scalar64 Test #0.2!");
ASSERT(HERE, !mi64_is_div_by_scalar64_u4(x, mod, lenX), "mi64_is_div_by_scalar64 Test #0.4!");
ASSERT(HERE, rem == mi64_div_by_scalar64   (x, mod, lenX, y), "mi64_div_by_scalar64 Test #0.1!");
ASSERT(HERE, rem == mi64_div_by_scalar64_u2(x, mod, lenX, y), "mi64_div_by_scalar64 Test #0.2!");
ASSERT(HERE, rem == mi64_div_by_scalar64_u4(x, mod, lenX, y), "mi64_div_by_scalar64 Test #0.4!");

// A second odd modulus: Use x/2 for modding here to help debug the even-modulus code tested by the next case:
mod = 8040689323464953445ull;
rem =    2985496175289855ull;
mi64_div2(x, y, lenX);
ASSERT(HERE, mi64_div_by_scalar64(y, mod, lenX, y) == rem, "mi64_div_by_scalar64 Test #1!");
// Use that to create an even modulus
mod=2*8040689323464953445ull;
rem =    5970992350579711ull;
ASSERT(HERE, mi64_div_by_scalar64(x, mod, lenX, y) == rem
			&& !mi64_mul_scalar(y, mod, y, lenX)
			&& (x[0] == y[0] + rem) && mi64_cmp_eq(x+1,y+1,lenX-1), "mi64_div_by_scalar64 Test #2!");
// And a modulus divisible by 4:
mod = 1635789749933632004ull;
rem =  816586808693101463ull;
ASSERT(HERE, mi64_div_by_scalar64(x, mod, lenX, y) == rem
			&& !mi64_mul_scalar(y, mod, y, lenX)
			&& (x[0] == y[0] + rem) && mi64_cmp_eq(x+1,y+1,lenX-1), "mi64_div_by_scalar64 Test #3!");
// And a modulus divisible by 8:
mod=2*1635789749933632004ull;
rem =  816586808693101463ull;
ASSERT(HERE, mi64_div_by_scalar64(x, mod, lenX, y) == rem
			&& !mi64_mul_scalar(y, mod, y, lenX)
			&& (x[0] == y[0] + rem) && mi64_cmp_eq(x+1,y+1,lenX-1), "mi64_div_by_scalar64 Test #4!");
// And a modulus divisible by 16:
mod=4*1635789749933632004ull;
rem =  4088166308560365471ull;
ASSERT(HERE, mi64_div_by_scalar64(x, mod, lenX, y) == rem
			&& !mi64_mul_scalar(y, mod, y, lenX)
			&& (x[0] == y[0] + rem) && mi64_cmp_eq(x+1,y+1,lenX-1), "mi64_div_by_scalar64 Test #5!");
// And a modulus divisible by 32:
mod=8*1635789749933632004ull;
rem =  4088166308560365471ull;
ASSERT(HERE, mi64_div_by_scalar64(x, mod, lenX, y) == rem
			&& !mi64_mul_scalar(y, mod, y, lenX)
			&& (x[0] == y[0] + rem) && mi64_cmp_eq(x+1,y+1,lenX-1), "mi64_div_by_scalar64 Test #6!");

	/* Test shift routines: */
	vec_len = MAX_ARRAY_DIM;
	for(i = 0; i < vec_len; i++)
	{
		x[i]  = rng_isaac_rand();
	}
	for(i = 0; i < 100; ++i)
	{
		j = (vec_len * 65)*rng_isaac_rand_double_norm_pos();
		if(j > (vec_len << 6))
		{
			mi64_shl (x, u, j, vec_len);	ASSERT(HERE, mi64_iszero(u, vec_len), "u != 0");
			mi64_shrl(x, v, j, vec_len);	ASSERT(HERE, mi64_iszero(v, vec_len), "v != 0");
		}
		else
		{
			mi64_shrl(x, v, j, vec_len);
			mi64_shl (v, v, j, vec_len);
											j = (vec_len << 6) - j;
			mi64_shl (x, u, j, vec_len);
			mi64_shrl(u, u, j, vec_len);
			mi64_add (u, v, v, vec_len);
			ASSERT(HERE, mi64_cmp_eq(x, v, vec_len), "x != v");
		}
	}

	/*
	Test #0: check the uint64[] <==> scalar double interconversion routine:
	*/
	vec_len = 1;
	x[0] = 0ull;
	dbl = mi64_cvt_double(x, vec_len);
	ASSERT(HERE, fabs(dbl - 0.0) < 1e-10, "GCD self-test 0a");

	x[0] = 1ull;
	dbl = mi64_cvt_double(x, vec_len);
	ASSERT(HERE, fabs(dbl - 1.0) < 1e-10, "GCD self-test 0b");

	x[0] = 2ull;
	dbl = mi64_cvt_double(x, vec_len);
	ASSERT(HERE, fabs(dbl - 2.0) < 1e-10, "GCD self-test 0c");

	x[0] = 16357897499336320021ull;
	dbl = mi64_cvt_double(x, vec_len);
	ASSERT(HERE, fabs(dbl - 16357897499336320021.0)/16357897499336320021.0 < 1e-10, "GCD self-test 0d");

	/* x = 46189291499145878679976776583847887373 */
	vec_len = 2;
	x[0] = 14921004030550884877ull;
	x[1] =  2503926509447010181ull;
	dbl = mi64_cvt_double(x, vec_len);
	ASSERT(HERE, fabs(dbl - 46189291499145878679976776583847887373.0)/46189291499145878679976776583847887373.0 < 1e-10, "GCD self-test 0e");

	mi64_set_eq(y, y10, 10);
	mi64_add_scalar(y, 1ull, y, 10);	/* p+1 */
	ASSERT(HERE, mi64_is_div_by_scalar64(y,     30ull, 10) == TRUE, "GCD self-test #0g");
	ASSERT(HERE, mi64_is_div_by_scalar64(y,    997ull, 10) == TRUE, "GCD self-test #0h");
	ASSERT(HERE, mi64_is_div_by_scalar64(y,   2113ull, 10) == TRUE, "GCD self-test #0i");
	ASSERT(HERE, mi64_is_div_by_scalar64(y,  87643ull, 10) == TRUE, "GCD self-test #0j");
	ASSERT(HERE, mi64_is_div_by_scalar64(y, 219607ull, 10) == TRUE, "GCD self-test #0k");	i = (uint32)6; i *= 2113; i *= 219607;	/* Product of selected small factors slightly below 2^64 */
	ASSERT(HERE, mi64_is_div_by_scalar64(y, (uint64)i, 10) == TRUE, "GCD self-test #0l");	i = (uint32)87643;
	ASSERT(HERE, mi64_is_div_by_scalar64(y,5ull*997*i, 10) == TRUE, "GCD self-test #0m");

	ASSERT(HERE, mi64_is_div_by_scalar32((uint32 *)y,     30, 10) == TRUE, "GCD self-test #0n");
	ASSERT(HERE, mi64_is_div_by_scalar32((uint32 *)y,    997, 10) == TRUE, "GCD self-test #0o");
	ASSERT(HERE, mi64_is_div_by_scalar32((uint32 *)y,   2113, 10) == TRUE, "GCD self-test #0p");
	ASSERT(HERE, mi64_is_div_by_scalar32((uint32 *)y,  87643, 10) == TRUE, "GCD self-test #0q");
	ASSERT(HERE, mi64_is_div_by_scalar32((uint32 *)y, 219607, 10) == TRUE, "GCD self-test #0r");	i = (uint32)6; i *= 2113; i *= 219607;	/* Product of selected small factors slightly below 2^32 */
	ASSERT(HERE, mi64_is_div_by_scalar32((uint32 *)y,      i, 10) == TRUE, "GCD self-test #0s");	i = (uint32)87643;
	ASSERT(HERE, mi64_is_div_by_scalar32((uint32 *)y,5*997*i, 10) == TRUE, "GCD self-test #0t");

	vec_len = 1;
	while(vec_len <= 1024)
	{
		for(j=0;j<100;++j)
		{
			/* Init an array of (vec_len) quasirandom 64-bit words: */
			for(i = 0; i < vec_len; i++)
			{
				x[i]  = rng_isaac_rand();
				tmp64 = rng_isaac_rand();
			}
			x[vec_len] = mi64_mul_scalar(x, tmp64, x, vec_len);
			ASSERT(HERE, mi64_is_div_by_scalar64(x, tmp64, vec_len+1) == TRUE, "mi64_mod_scalar test!");
		}
		vec_len *= 2;
		/* Want to do at least one test right at the max. length: */
		if(vec_len > MAX_ARRAY_DIM)
		{
			if(vec_len >= 2*MAX_ARRAY_DIM)
				break;
			else
				vec_len = MAX_ARRAY_DIM;
		}
	}

	/************* Test #??: shakedown tests of loop-folded DIV algo: ****************/
	printf("\nPerforming GCD self-test DIV...\n");
	vec_len = MAX_ARRAY_DIM;
	for(i = 0; i < vec_len; i++)
	{
		x[i]  = rng_isaac_rand();
	}
	for(i = 4; i < MAX_ARRAY_DIM; i += 4)
	{
		mod = rng_isaac_rand() | 1;	// Current loop-folded DIV only allows odd modulus
		rem   = mi64_div_by_scalar64   (x, mod, i, y);
		tmp64 = mi64_div_by_scalar64_u2(x, mod, i, z);	ASSERT(HERE, rem == tmp64 && mi64_cmp_eq(y,z,i), "GCD self-test #DIV");
		tmp64 = mi64_div_by_scalar64_u4(x, mod, i, z);	ASSERT(HERE, rem == tmp64 && mi64_cmp_eq(y,z,i), "GCD self-test #DIV");
	}
	/********************************************************************************/

	/*
	Test #1: check the uint64[] pure-integer multiply routines - these use grammar-school algorithm:
	*/
	vec_len = 1;
	while(vec_len <= 1024)
	{
		/* Init an array of (vec_len) quasirandom 64-bit words: */
		for(i = 0; i < vec_len; i++)
		{
			x[i] = rng_isaac_rand();
			y[i] = rng_isaac_rand();
		}

		mi64_mul_vector(x, vec_len, y, vec_len, z, &lenZ);
		mi64_mul_vector_lo_half(x, y, u, vec_len);
		ASSERT(HERE, mi64_cmp_eq(u, z, vec_len), "GCD self-test 1: mi64_mul_vector_lo_half");
		mi64_mul_vector_hi_half(x, y, v, vec_len);
		ASSERT(HERE, mi64_cmp_eq(v, &z[vec_len], vec_len), "GCD self-test 1: mi64_mul_vector_hi_half");

		/******************************************************************/
		/* ****TODO**** Compare pure-integer x*y with FFT-based multiply: */
		/******************************************************************/

		vec_len *= 2;
		/* Want to do at least one test right at the max. length: */
		if(vec_len > MAX_ARRAY_DIM)
		{
			if(vec_len >= 2*MAX_ARRAY_DIM)
				break;
			else
				vec_len = MAX_ARRAY_DIM;
		}
	}

	/*
	Test #2: let p = 16357897499336320021 (64-bit test prime), check that mi64_pprimeF to various small-prime bases returns TRUE:
	*/
	printf("\nPerforming GCD self-test #2...\n");
	vec_len = 1;
	x[0] = 43112609ull;
	ASSERT(HERE, pprimeF   ((uint32)x[0], 2ull), "GCD self-test #2,02");
	ASSERT(HERE, mi64_pprimeF(x, 2ull, vec_len), "GCD self-test #2, 2");
	y[0] = 16357897499336320021ull;
	ASSERT(HERE, mi64_div_y32(y, (uint32)x[0], u, 1) == 7915399 && u[0] == 379422583758ull, "GCD self-test #3f");
	ASSERT(HERE, mi64_pprimeF(y, 2ull, vec_len), "GCD self-test #2, 2");
	ASSERT(HERE, mi64_pprimeF(y, 3ull, vec_len), "GCD self-test #2, 3");
	ASSERT(HERE, mi64_pprimeF(y, 5ull, vec_len), "GCD self-test #2, 5");
	ASSERT(HERE, mi64_pprimeF(y, 7ull, vec_len), "GCD self-test #2, 7");
	ASSERT(HERE, mi64_pprimeF(y,11ull, vec_len), "GCD self-test #2,11");
	ASSERT(HERE, mi64_pprimeF(y,13ull, vec_len), "GCD self-test #2,13");
	ASSERT(HERE, mi64_pprimeF(y,17ull, vec_len), "GCD self-test #2,17");
	ASSERT(HERE, mi64_pprimeF(y,19ull, vec_len), "GCD self-test #2,19");
	ASSERT(HERE, mi64_pprimeF(y,23ull, vec_len), "GCD self-test #2,23");
	ASSERT(HERE, mi64_pprimeF(y,29ull, vec_len), "GCD self-test #2,29");

	/*
	Test #2a: 128-bit test primes:
	*/
	vec_len = 2;
	for(j = 0; testp128[j].d0 != 0; ++j)
	{
		x[0] = testp128[j].d0;	x[1] =  testp128[j].d1;
		/*
		fprintf(stderr, "mi64_pprimeF: testing p = %s\n", &char_buf[convert_uint128_base10_char(char_buf, testp128[j], 0)]);
		*/
		ASSERT(HERE, mi64_pprimeF(x, 3ull, vec_len) == 1, "GCD self-test #2a");
	}

	/*
	Test #2b: 192-bit test primes:
	*/
	vec_len = 3;
	for(j = 0; testp192[j].d0 != 0; ++j)
	{
		x[0] = testp192[j].d0;	x[1] =  testp192[j].d1;	x[2] =  testp192[j].d2;
		/*
		fprintf(stderr, "mi64_pprimeF: testing p = %s\n", &char_buf[convert_uint192_base10_char(char_buf, testp192[j], 0)]);
		*/
		ASSERT(HERE, mi64_pprimeF(x, 3ull, vec_len) == 1, "GCD self-test #2b");
	}

	ASSERT(HERE, mi64_pprimeF(y1 , 3ull, 1 ) == 1, "GCD self-test #2.1 ");
	ASSERT(HERE, mi64_pprimeF(y2 , 3ull, 2 ) == 1, "GCD self-test #2.2 ");
	ASSERT(HERE, mi64_pprimeF(y3 , 3ull, 3 ) == 1, "GCD self-test #2.3 ");
	ASSERT(HERE, mi64_pprimeF(y4 , 3ull, 4 ) == 1, "GCD self-test #2.4 ");
	ASSERT(HERE, mi64_pprimeF(y5 , 3ull, 5 ) == 1, "GCD self-test #2.5 ");
	ASSERT(HERE, mi64_pprimeF(y6 , 3ull, 6 ) == 1, "GCD self-test #2.6 ");
	ASSERT(HERE, mi64_pprimeF(y7 , 3ull, 7 ) == 1, "GCD self-test #2.7 ");
	ASSERT(HERE, mi64_pprimeF(y8 , 3ull, 8 ) == 1, "GCD self-test #2.8 ");
	ASSERT(HERE, mi64_pprimeF(y9 , 3ull, 9 ) == 1, "GCD self-test #2.9 ");
	ASSERT(HERE, mi64_pprimeF(y10, 3ull, 10) == 1, "GCD self-test #2.10");

	/*
	Test #3: Test mi64 scalar-mul routines:
	let q = 16357897499336320021;
	x = q * {some quasi-random multiword int D}
	u,v,x,y available for storage.
	*/
	printf("\nPerforming GCD self-test #3...\n");
	vec_len = 1;
	for(j = 0; j < 10; ++j)
	{
		vec_len = 100*j + 1;
		for(i = 0; i < vec_len; i++)
		{
			x[i] = rng_isaac_rand();
		}
		mi64_set_eq(y, x, vec_len);									/* y = D */
		x[vec_len] = mi64_mul_scalar(x, 16357897499336320021ull, x, vec_len);	/* x = q*D */
		++vec_len;	/* Changed x[vec_len++] = ... to a separate increment step to avoid order-of-eval ambiguity. */
		ASSERT(HERE, mi64_getlen(x, vec_len) == vec_len , "GCD self-test #3: expected nonzero carryout of mi64_mul_scalar!");
		mi64_set_eq_scalar(u, 16357897499336320021ull, vec_len);	/* u = q */
		mi64_div(x, u, vec_len, vec_len, v, x);	/* In-place x/q: dividend put into v, remainder back into x */
		ASSERT(HERE, mi64_cmp_eq       (y,  v, vec_len-1), "GCD self-test #3: (x*q)/q != x");
		ASSERT(HERE, mi64_cmp_eq_scalar(x, 0ull, vec_len), "GCD self-test #3: (x*q)%q != 0");
	}

	/*
	Test #4: Test mi64 scalar-div and remainder routines:
	x = p * {some quasi-random (N-1)-word int}
	y =     {some quasi-random  N   -word prime}, N of various lengths for which we have prestored multiword primes.
	*/
	printf("\nPerforming GCD self-test #4...\n");
	mi64_set_eq(y, y10, 10);
	mi64_add_scalar(y, 1ull, y, 10);	/* p+1 */
	ASSERT(HERE, mi64_div_y32(y,    30, 0x0, 10) == 0, "GCD self-test #4a");
	ASSERT(HERE, mi64_div_y32(y,   997, 0x0, 10) == 0, "GCD self-test #4b");
	ASSERT(HERE, mi64_div_y32(y,  2113, 0x0, 10) == 0, "GCD self-test #4c");
	ASSERT(HERE, mi64_div_y32(y, 87643, 0x0, 10) == 0, "GCD self-test #4d");
	ASSERT(HERE, mi64_div_y32(y,219607, 0x0, 10) == 0, "GCD self-test #4e");
	i = (uint32)6; i *= 2113; i *= 219607;	/* Product of selected small factors slightly below 2^32 */
	ASSERT(HERE, mi64_div_y32(y,     i,  u , 10) == 0, "GCD self-test #4f");
	i = (uint32)87643;
	ASSERT(HERE, mi64_div_y32(u,5*997*i, u , 10) == 0, "GCD self-test #4f");

	mi64_sub_scalar(y, 1ull, y, 10);	/* p */
	ASSERT(HERE, mi64_div_y32(y, (uint32)3578974993u, u, 10) == (uint32)2756613022u, "GCD self-test #4g");
	ASSERT(HERE, mi64_mul_scalar(u, (uint32)3578974993u, u, 10) == 0, "GCD self-test #4h");
	ASSERT(HERE, mi64_add_scalar(u, (uint32)2756613022u, u, 10) == 0, "GCD self-test #4i");
	ASSERT(HERE, mi64_cmp_eq(u, y10,10), "GCD self-test #4j");

	/*
	Test #5: let p = 16357897499336320049 [next-higher 64-bit prime above our 1-word reference prime]
	x = p * {some quasi-random (N-1)-word int}
	y =     {some quasi-random  N   -word prime}, N of various lengths for which we have prestored multiword primes.
	*/
	printf("\nPerforming GCD self-test #5...\n");
	for(j = 1; j < num_mword_prime; j++)	/* Start multipliers with 2-word primes and go up from there */
	{
		gcd_debug=0;
		vec_len = len_mword_prime[j]-1;
		curr_mword_prime = ptr_mword_prime[j];
		for(i = 0; i < vec_len; i++)
		{
			x[i] = rng_isaac_rand();
		}
		x[vec_len] = mi64_mul_scalar(x, 16357897499336320049ull, x, vec_len);
		ASSERT(HERE, x[vec_len] != 0, "GCD self-test #5: expected nonzero carryout of mi64_mul_scalar!");
		vec_len += 1;
		for(i = 0; i < vec_len; i++)
		{
			y[i] = curr_mword_prime[i];
		}
		fail += (mi64_gcd(x, y, vec_len, eGCD, &Ap, &Bp, &len_AB, &sign_AB) != 1);
		fail += (x[0] != 1ull);	/* Put this on a separate line because some compilers will return 1 as the function result of evaluating (a+b), i.e. fail += a+b increments fail by 1 even if a = b = 0. Strange but true... */
		ASSERT(HERE, fail==0,"GCD self-test #5");
	}

	/*
	Test #6: let p = 2-word prime
	y = p * {some quasi-random N-word prime, stored in y}
	x = p * {some quasi-random N-word int coprime to y[]}
	*/
	printf("\nPerforming GCD self-test #6...\n");
	for(j = 2; j < num_mword_prime; j++)	/* Start multipliers with 3-word primes and go up from there */
	{
		gcd_debug=0;
		vec_len = len_mword_prime[j];
		curr_mword_prime = ptr_mword_prime[j];
		for(i = 0; i < vec_len; i++)
		{
			v[i] = curr_mword_prime[i];
			u[i] = rng_isaac_rand();
		}
		ASSERT(HERE, mi64_pprimeF(v, 3ull, vec_len) == 1, "GCD self-test #6.prp");
		/* vector*vector MULs to get x and y are here: */
		mi64_mul_vector(y2, 2, v, vec_len, y, &i      );
		mi64_mul_vector(y2, 2, u, vec_len, x, &vec_len);
		vec_len = MAX(i, vec_len);
		fail += (mi64_gcd(x, y, vec_len, eGCD, &Ap, &Bp, &len_AB, &sign_AB) != 2);
		fail += (mi64_cmp_eq(x, y2, 2) != TRUE);
		ASSERT(HERE, fail==0,"GCD self-test #6");
	}

	/* Test #7: Use various small primes to test out the mi64_twopmodq routine: */
	printf("\nPerforming GCD self-test #7...\n");
	for(j = 0; j < num_mword_prime; j++)
	{
		vec_len = len_mword_prime[j];
		curr_mword_prime = ptr_mword_prime[j];
		mi64_set_eq(x, curr_mword_prime, vec_len);	/* x = p */
		mi64_sub_scalar(x, 1ull, y, vec_len);			/* y = p-1 */
		ASSERT(HERE, mi64_twopmodq(y, vec_len, 0, x, vec_len, 0x0) == 1, "GCD self-test #7: mi64_twopmodq != 1");	/* 2^(p-1) ?= 1 (mod p) */
	}

	/*
	Test #8: Multiply the following factors of M-numbers together (where w:= 2^64):

		M(    677):    157590042578912*w^2 + 10558642444782195772*w +   329809049266961143;
		M(    773):      9118322195022*w^2 +  1933308633079010416*w + 17814616685598394119;
		M(    971):     70286054459973*w^2 + 17012949627558354271*w +  3547755741880899889;
		M(    997): 492416983078691417*w^2 +  8040689323464953445*w + 16007877010440112335;
		M(   1001):        59364131986*w^2 +  9565712986615012496*w + 10050950882119470361;

	The product of the 5 M-numbers having these small factors needs P uint64 slots to store in compact form.
	Create a nontrivial y-vector by multiplying the F 64-bit words of the factor product by (P-F) random 64-bit ints.

	We also use this section to test the mi64 string I/O functionality.
	*/
	printf("\nPerforming GCD self-test #8...\n");

	/*** First multiply together all the M-numbers, storing the result in x[]: ***/
	/* M677: */
	p = 677;	lenX = (p>>6)+1;
	memset(x,ONES64,(lenX<<3));	x[lenX-1] = (1ull << (p&63)) - 1;
	ASSERT(HERE, STREQ(&char_buf[convert_mi64_base10_char(char_buf, x, lenX, 0)], "627057063764139831929324851379409869378845668175598843037877190478889006888518431438644711527536922839520331484815861906173161536477065546885468336421475511783984145060592245840032548652210559519683510271"), "GCD self-test #8");
	/* q677: */
	u[2] =    157590042578912ull; u[1] = 10558642444782195772ull; u[0] =   329809049266961143ull;	lenU = 3;
	ASSERT(HERE, STREQ(&char_buf[convert_mi64_base10_char(char_buf, u, lenU, 0)], "53625112691923843508117942311516428173021903300344567"), "GCD self-test #8");
	/* Test the back-conversion (string -> mi64): */
	ASSERT(HERE, 0x0 != (out_array = convert_base10_char_mi64("53625112691923843508117942311516428173021903300344567", &i)) && (i == 3), "0");
	ASSERT(HERE, mi64_cmp_eq(out_array, u, 3), "0");	free((void *)out_array);	out_array = 0x0;

	ASSERT(HERE, mi64_pprimeF(u, 3ull, lenU) == 1, "GCD self-test #8.pr1");
//	mi64_setlen(u, lenU, lenX);	only needed for legacy mi64_div_binary function
	mi64_div(x, u, lenX, lenU, w, z);	/* Check divisibility directly - quotient returned in w, remainder in z */
	ASSERT(HERE, STREQ(&char_buf[convert_mi64_base10_char(char_buf, w, lenX, 0)], "11693347245097298120823316616845179432964234640294687969094035694006198722115249470397309583454281704888388908348322516169275641497746779501202976102713"), "GCD self-test #8");
	ASSERT(HERE, mi64_iszero(z, lenU) == 1, "GCD div-test #8.pr1");
	ASSERT(HERE, STREQ(&char_buf[convert_mi64_base10_char(char_buf, z, lenU, 0)], "0"), "GCD self-test #8");

	/* ... *= M773: */
	p = 773;	lenY = (p>>6)+1;
	memset(y,ONES64,(lenY<<3));	y[lenY-1] = (1ull << (p&63)) - 1;
	ASSERT(HERE, STREQ(&char_buf[convert_mi64_base10_char(char_buf, y, lenY, 0)], "49680578953622685924767343630800081768220352547734291556449665216833630485964060362588109082516687294415607382308194342597490561411674060526217192801317796454542559232667196977608489140211150234408415974198927000028571099322113851391"), "GCD self-test #8");
	mi64_set_eq(w, x, lenX);	/* Use copy of X as multiplicand 1 so we can store the result back into X */
	mi64_mul_vector(w, lenX, y, lenY, x, &lenX);
	ASSERT(HERE, STREQ(&char_buf[convert_mi64_base10_char(char_buf, x, lenX, 0)], "31152557964761163897928842168277492768812223183099899891668902917976925560029252973898899702640514119214599026359946447505761700950328538576680016552349084733366485698028733665842440380429912269143850331229298040934224944520428395437344934546472616298890629611082440248854931228400162480344449363195362608681181976242618371232062995112351931255191940463434239006436301912973135868237655715193561884369357401722616533384552395786116136961"), "GCD self-test #8");
	/* ... *= q773: */
	v[2] =      9118322195022ull; v[1] =  1933308633079010416ull; v[0] = 17814616685598394119ull; lenV = 3;
	ASSERT(HERE, mi64_pprimeF(v, 3ull, lenV) == 1, "GCD self-test #8.pr2");
	mi64_set_eq(z, u, lenU);
	mi64_mul_vector(z, lenU, v, lenV, u, &lenU);
	ASSERT(HERE, STREQ(&char_buf[convert_mi64_base10_char(char_buf, u, lenU, 0)], "166388228042876887900053740282475640168826969539724862936313620741285381200003507633306121185343797304769"), "GCD self-test #8");
//	mi64_setlen(u, lenU, lenX);	only needed for legacy mi64_div_binary function
	mi64_div(x, u, lenX, lenU, w, z);	/* Check divisibility directly - quotient returned in w, remainder in z */
	ASSERT(HERE, mi64_iszero(z, lenU) == 1, "GCD div-test #8.pr2");

	/* ... *= M971: */
	p = 971;	lenY = (p>>6)+1;
	memset(y,ONES64,(lenY<<3));	y[lenY-1] = (1ull << (p&63)) - 1;
	mi64_set_eq(w, x, lenX);	/* Use copy of X as multiplicand 1 so we can store the result back into X */
	mi64_mul_vector(w, lenX, y, lenY, x, &lenX);
	/* ... *= q971: */
	v[2] =     70286054459973ull; v[1] = 17012949627558354271ull; v[0] =  3547755741880899889ull; lenV = 3;
	ASSERT(HERE, mi64_pprimeF(v, 3ull, lenV) == 1, "GCD self-test #8.pr3");
	mi64_set_eq(z, u, lenU);
	mi64_mul_vector(z, lenU, v, lenV, u, &lenU);
//	mi64_setlen(u, lenU, lenX);	only needed for legacy mi64_div_binary function
	mi64_div(x, u, lenX, lenU, w, z);	/* Check divisibility directly - quotient returned in w, remainder in z */
	ASSERT(HERE, mi64_iszero(z, lenU) == 1, "GCD div-test #8.pr3");

	/* ... *= M997: */
	p = 997;	lenY = (p>>6)+1;
	memset(y,ONES64,(lenY<<3));	y[lenY-1] = (1ull << (p&63)) - 1;
	mi64_set_eq(w, x, lenX);	/* Use copy of X as multiplicand 1 so we can store the result back into X */
	mi64_mul_vector(w, lenX, y, lenY, x, &lenX);
	/* ... *= q997: */
	v[2] = 492416983078691417ull; v[1] =  8040689323464953445ull; v[0] = 16007877010440112335ull; lenV = 3;
	ASSERT(HERE, mi64_pprimeF(v, 3ull, lenV) == 1, "GCD self-test #8.pr3");
	mi64_set_eq(z, u, lenU);
	mi64_mul_vector(z, lenU, v, lenV, u, &lenU);
//	mi64_setlen(u, lenU, lenX);	only needed for legacy mi64_div_binary function
	mi64_div(x, u, lenX, lenU, w, z);	/* Check divisibility directly - quotient returned in w, remainder in z */
	ASSERT(HERE, mi64_iszero(z, lenU) == 1, "GCD div-test #8.pr3");

	/* ... *= M1001: */
	p =1001;	lenY = (p>>6)+1;
	memset(y,ONES64,(lenY<<3));	y[lenY-1] = (1ull << (p&63)) - 1;
	mi64_set_eq(w, x, lenX);	/* Use copy of X as multiplicand 1 so we can store the result back into X */
	mi64_mul_vector(w, lenX, y, lenY, x, &lenX);
	/* ... *= q1001: */
	v[2] =        59364131986ull; v[1] =  9565712986615012496ull; v[0] = 10050950882119470361ull; lenV = 3;
	ASSERT(HERE, mi64_pprimeF(v, 3ull, lenV) == 1, "GCD self-test #8.pr4");
	mi64_set_eq(z, u, lenU);
	mi64_mul_vector(z, lenU, v, lenV, u, &lenU);
//	mi64_setlen(u, lenU, lenX);	only needed for legacy mi64_div_binary function
	mi64_div(x, u, lenX, lenU, w, z);	/* Check divisibility directly - quotient returned in w, remainder in z */
	ASSERT(HERE, mi64_iszero(z, lenU) == 1, "GCD div-test #8.pr4");

	/*** Lastly, pad the q-product stored in U] to be the same length as the M-product [stored in X]
	by multiplying the q-product by as many random 64-bit digits as are needed to do so: ***/
	lenV = lenX - lenU;
	ASSERT(HERE, lenV <= MAX_ARRAY_DIM, "GCD self-test #8: lenV > MAX_ARRAY_DIM!");
	for(i = 0; i < lenV; i++)
	{
		v[i] = rng_isaac_rand();
	}
	mi64_mul_vector(u, lenU, v, lenV, y, &lenY);
	ASSERT(HERE, lenY >= (lenX - 1), "GCD self-test #8: mi64_mul_vector(u, lenU, v, lenV, y, &lenV) < (vec_len - 1)");
	mi64_setlen(y, lenY, lenX);
	gcd_debug=0;

	/* First use the data to test mi64_div - add a 10-word const, quotient returned in w, remainder in z, rem = the const we added */
	tmp64 = mi64_add(x, y10, x, 10);	mi64_add_scalar(x+10, tmp64, x+10, lenX-10);
	mi64_div(x, u, lenX, lenU, w, z);
	ASSERT(HERE, mi64_cmp_eq(z, y10, 10), "GCD div-test #8.pr5: Remainder check failed");
	// Subtract the 10-word const back off:
	tmp64 = mi64_sub(x, y10, x, 10);	mi64_sub_scalar(x+10, tmp64, x+10, lenX-10);
	// ...And check that quotient*divisor == x:
	mi64_mul_vector(w,lenX,u,lenU,z,&i);
	ASSERT(HERE, i == lenX, "GCD div-test #8.pr5: Remultiply length check failed");
	ASSERT(HERE, mi64_cmp_eq(x, z, lenX), "GCD div-test #8.pr5: Remultiply product check failed");

	/* Get the GCD: */
	vec_len = mi64_gcd(x, y, lenX, eGCD, &Ap, &Bp, &len_AB, &sign_AB);

	/* Note that 1001 is in fact composite (7*11*13) so M1001 has many small factors. Thus, since the
	random multiplier we used to pad the large-prime-factor product ma well share some of those small multipliers,
	we take a further GCD with the unpadded factor product:
	*/
	lenX = vec_len;
	mi64_setlen(u, lenU, lenX);
	vec_len = mi64_gcd(x, u, lenX, eGCD, &Ap, &Bp, &len_AB, &sign_AB);
	ASSERT(HERE, vec_len == 14,"GCD self-test #8");
	fail = 0;
	fail += (x[ 0] !=  3593788887684319431ull);
	fail += (x[ 1] != 10937322091449271811ull);
	fail += (x[ 2] != 12946497122050840313ull);
	fail += (x[ 3] !=   647314691296882413ull);
	fail += (x[ 4] != 16358554045068224992ull);
	fail += (x[ 5] != 11116777718169420623ull);
	fail += (x[ 6] !=  1469016928937333456ull);
	fail += (x[ 7] !=  9391229552579049908ull);
	fail += (x[ 8] !=  5839179435288346904ull);
	fail += (x[ 9] !=  8876202416807555851ull);
	fail += (x[10] !=  5810184814721434260ull);
	fail += (x[11] !=  4707921738239363748ull);
	fail += (x[12] != 13066869379195810396ull);
	fail += (x[13] !=         470338845646ull);
	ASSERT(HERE, fail==0,"GCD self-test #8");
	ASSERT(HERE, STREQ(&char_buf[convert_mi64_base10_char(char_buf, x, vec_len, 0)], "13469989009602505896761896298151251254620625588479300721021471299084968313453138000611404859588356783789058816019125967988627321691441338348931654165690484559185870012613932145949653360324772652694002027486449514263873807210595986656151451735452053543049459632327"), "GCD self-test #8");

	/*
	Test #9: M27691 has the 121-bit factor 1734072082042172647364731231822850071.
	This M-number needs 433 uint64 slots to store in compact form. Create a nontrivial
	y-vector by multiplying the two 64-bit words of the factor by the 431 LSWs of y4[].
	*/
if(MAX_ARRAY_DIM >= 433)
{
	printf("\nPerforming GCD self-test #9...\n");

	p = 27691;	vec_len = (p>>6)+1;

	// First try a gcd of 2^p-1 (huge) with just the small factor:
	memset(x,ONES64,(vec_len<<3));	x[vec_len-1] = (1ull << (p&63)) - 1;
	memset(u,(uint64) 0,(vec_len<<3));
	/* Store q in u[]: */
	u[0] =  4235226679561594903ull;
	u[1] =    94004235929829273ull;
	ASSERT(HERE, mi64_pprimeF(u, 3ull, 2) == 1, "GCD self-test #9.prp");
	clock1 = clock();
	gcd_debug=0;
	fail += (mi64_gcd(x, u, vec_len, eGCD, &Ap, &Bp, &len_AB, &sign_AB) != 2);
	fail += (x[0] !=  4235226679561594903ull);
	fail += (x[1] !=    94004235929829273ull);
	ASSERT(HERE, fail==0,"GCD self-test #9a");

	clock2 = clock();
	tdiff = (double)(clock2 - clock1);
	printf	("GCD self-test #9a passed: GCD Time =%s\n",get_time_str(tdiff));

	// Now re-init x and u and do full-length gcd:
	mi64_set_eq(u, x, vec_len);
	memset(x,ONES64,(vec_len<<3));	x[vec_len-1] = (1ull << (p&63)) - 1;
	// Replace RHS with y431_test[i] here to recreate garbage-repeat-high-words testcase of 5/13/2012
	for(i = 0; i < 431; i++)
	{
		v[i] = y431_test[i];//rng_isaac_rand();
	//	printf("v[%3u] = %20llu\n",i,v[i]);
	}
	/* vector*vector MUL to get y is here: */
	mi64_mul_vector(u, 2, v, 431, y, &vec_len);
	ASSERT(HERE, vec_len == 433, "GCD self-test #9: vec_len == 433");
	clock1 = clock();
	gcd_debug=0;
	fail += (mi64_gcd(x, y, vec_len, eGCD, &Ap, &Bp, &len_AB, &sign_AB) != 2);
	fail += (x[0] !=  4235226679561594903ull);
	fail += (x[1] !=    94004235929829273ull);
	ASSERT(HERE, fail==0,"GCD self-test #9b");

	clock2 = clock();
	tdiff = (double)(clock2 - clock1);
	printf	("GCD self-test #9b passed: GCD Time =%s\n",get_time_str(tdiff));
}

	/*
	Test #10: M80239 has the factor q = 39654784768949071.
	This M-number needs 1254 uint64 slots to store in compact form. Create a nontrivial
	y-vector by multiplying the factor by a number whose base-2^64 digits consist of
	1253 randomly selected words of y4[]. Owing to
	the special form of M-number divisors it's exceedingly unlikely that the resulting y[]
	will have any factors in common with x[] other than q.
	*/
if(MAX_ARRAY_DIM >= 1254)
{
	printf("\nPerforming GCD self-test #10...\n");

	p = 80239;	vec_len = (p>>6)+1;

	// First try a gcd of 2^p-1 (huge) with just the small factor:
	memset(x,ONES64,(vec_len<<3));	x[vec_len-1] = (1ull << (p&63)) - 1;
	memset(u,(uint64) 0,(vec_len<<3));	u[0] = 39654784768949071ull;
	clock1 = clock();
	gcd_debug=0;
	fail += (mi64_gcd(x, u, vec_len, eGCD, &Ap, &Bp, &len_AB, &sign_AB) != 1);
	fail += (x[0] !=    39654784768949071ull);
	ASSERT(HERE, fail==0,"GCD self-test #10a");

	clock2 = clock();
	tdiff = (double)(clock2 - clock1);
	printf	("GCD self-test #10a passed: GCD Time =%s\n",get_time_str(tdiff));

	// Now re-init x and u and do full-length gcd:
	memset(x,ONES64,(vec_len<<3));	x[vec_len-1] = (1ull << (p&63)) - 1;
	for(i = 0; i < vec_len-1; i++)
	{
		y[i] = rng_isaac_rand();
	}
	y[vec_len-1] = mi64_mul_scalar(y, 39654784768949071ull, y, vec_len-1);
	ASSERT(HERE, y[vec_len-1] != 0, "GCD self-test #10: expected nonzero carryout of q0* y4[0:vec_len-1]!");

	clock1 = clock();
	gcd_debug=0;
	fail += (mi64_gcd(x, y, vec_len, eGCD, &Ap, &Bp, &len_AB, &sign_AB) != 1);
	fail += (x[0] !=    39654784768949071ull);
	ASSERT(HERE, fail==0,"GCD self-test #10");

	clock2 = clock();
	tdiff = (double)(clock2 - clock1);
	printf	("GCD self-test #10b passed: GCD Time =%s\n",get_time_str(tdiff));
}

	/*
	Test #11: M80239 has the factor q = 39654784768949071.
	Same as above, but in this case give the mi64_mul_vector routine a workout by
	multiplying M80239 by a quasirandom 747-word input to get a 2001-word x-input
	and multiplying q by the product of 2 quasirandom 1000-word inputs to get a 2001-word y.
	*/
gcd_debug=0;

if(MAX_ARRAY_DIM >= 1254)
{
	printf("\nPerforming GCD self-test #11...\n");

	p = 80239;	vec_len = (p>>6)+1;
	memset(u,ONES64,(vec_len<<3));	u[vec_len-1] = (1ull << (p&63)) - 1;
	/* Multiply u by random multiword integer to get 2001-word product: */
	for(i = 0; i < 2001-vec_len; i++)
	{
		v[i] = rng_isaac_rand();
	}

	/* vector*vector MUL to get x is here: */
	mi64_mul_vector(u, vec_len, v, 2001-vec_len, x, &vec_len);
	ASSERT(HERE, vec_len == 2001, "GCD self-test #11: vec_len == 2001");

	/* y: */
	for(i = 0; i < 1000; i++)
	{
		u[i] = rng_isaac_rand();
		v[i] = rng_isaac_rand();
	}

	/* vector*vector MUL to get y is here: */
	mi64_mul_vector(u, 1000, v, 1000, y, &vec_len);
	ASSERT(HERE, vec_len == 2000, "GCD self-test #11: vec_len == 2000");
	/* ...and in-place-multiply the result by q: */
	y[vec_len] = mi64_mul_scalar(y, 39654784768949071ull, y, vec_len);
	++vec_len;

	/* Because of the way we constructed the test vectors they may have other small-prime
	factors in common, so weed those out by taking a second short-length GCD, if necessary: */
	clock1 = clock();
	vec_len = mi64_gcd(x, y, vec_len, eGCD, &Ap, &Bp, &len_AB, &sign_AB);
	if(vec_len != 1 || x[0] != 39654784768949071ull)
	{
		mi64_clear(y, vec_len);
		y[0] = 39654784768949071ull;
		fail += (mi64_gcd(x, y, vec_len, eGCD, &Ap, &Bp, &len_AB, &sign_AB) != 1);
		fail += (x[0] != 39654784768949071ull);
	}
	ASSERT(HERE, fail==0,"GCD self-test #11");

	clock2 = clock();
	tdiff = (double)(clock2 - clock1);
	printf	("GCD self-test #11 passed: GCD Time =%s\n",get_time_str(tdiff));
}

	/*
	Test #12: M2202817 has the factor 87722769297534671.
	This M-number needs 34420 uint64 slots to store in compact form. Create a nontrivial
	y-vector by multiplying the factor by 1253 randomly selected words of y4[]. Owing to
	the special form of M-number divisors it's exceedingly unlikely that the resulting y[]
	will have any factors in common with x[] other than q.

	04 May 2012: ~9 sec on a single core of a 2 GHz Core2 Duo.
	*/
if(MAX_ARRAY_DIM >= 34420)
{
	printf("\nPerforming GCD self-test #12...\n");
	p = 2202817;	vec_len = (p>>6)+1;
	/*
	2202817 = 1000011001110011000001_2
	So the LR binary powering sequence is:
	bit	power	action
	---	-------	---------
	1	1		x=2;
	0	2		x=x^2%q;
	0	4		x=x^2%q;
	0	8		x=x^2%q;
	0	16		x=x^2%q;
	1	33		x=x^2%q; x=2*x%q;
	1	67		x=x^2%q; x=2*x%q;
	0	134		x=x^2%q;
	0	268		x=x^2%q;
	1	537		x=x^2%q; x=2*x%q;
	1	1075	x=x^2%q; x=2*x%q;
	1	2151	x=x^2%q; x=2*x%q;
	0	4302	x=x^2%q;
	0	8604	x=x^2%q;
	1	17209	x=x^2%q; x=2*x%q;
	1	34419	x=x^2%q; x=2*x%q;
	0	68838	x=x^2%q;
	0	137676	x=x^2%q;
	0	275352	x=x^2%q;
	0	550704	x=x^2%q;
	0	1101408	x=x^2%q;
	1	2202817	x=x^2%q; x=2*x%q;
	Which allows the quotient ( = x-1 at the end of the above) to easily be computed using e.g. Pari/GP or bc.
	*/
	// First use this M(p) for a large-vector timing test of the mi64_div_by_scalar64 fucntion:
	memset(x,ONES64,(vec_len<<3));	x[vec_len-1] = (1ull << (p&63)) - 1;
#if 1	// Set = 1 to do 64-bit is-divisible timing test
	clock1 = clock();
	for(i = 0; i < 10000; i++) {
		ASSERT(HERE, mi64_is_div_by_scalar64_u4(x, 87722769297534671ull, vec_len), "GCD self-test #12.0");
	}
	clock2 = clock();
	tdiff = (double)(clock2 - clock1);
	printf	("GCD self-test #12.0 passed: mi64_is_div_by_scalar64(), %u*%u words processed, Time =%s\n",i,vec_len,get_time_str(tdiff));
//	exit(0);
#endif
#if 1	// Set = 1 to do 64-bit div-with-remainder timing test
	clock1 = clock();
	for(i = 0; i < 10000; i++) {
	//	ASSERT(HERE, mi64_div_by_scalar64   (x, 16357897499336320049ull, vec_len, y) == 14995895313315469881ull, "GCD self-test #12.0");
	//	ASSERT(HERE, mi64_div_by_scalar64_u2(x, 16357897499336320049ull, vec_len, y) == 14995895313315469881ull, "GCD self-test #12.0");
		ASSERT(HERE, mi64_div_by_scalar64_u4(x, 16357897499336320049ull, vec_len, y) == 14995895313315469881ull, "GCD self-test #12.0");
	}
	clock2 = clock();
	ASSERT(HERE,!mi64_mul_scalar(y, 16357897499336320049ull, y, vec_len), "GCD self-test #12.1");
	ASSERT(HERE, (y[0] += 14995895313315469881ull) && mi64_cmp_eq(x, y, vec_len), "GCD self-test #12.1");
	tdiff = (double)(clock2 - clock1);
	printf	("GCD self-test #12.1 passed: mi64_div_by_scalar64(), %u*%u words processed, Time =%s\n",i,vec_len,get_time_str(tdiff));
exit(0);
#endif
	// Next try a gcd of 2^p-1 (huge) with just the small factor:
	memset(y,(uint64) 0,(vec_len<<3));	y[0] = 87722769297534671ull;
	clock1 = clock();
	gcd_debug=0;
	fail += (mi64_gcd(x, y, vec_len, eGCD, &Ap, &Bp, &len_AB, &sign_AB) != 1);
	fail += (x[0] !=    87722769297534671ull);
	ASSERT(HERE, fail==0,"GCD self-test #12a");

	clock2 = clock();
	tdiff = (double)(clock2 - clock1);
	printf	("GCD self-test #12a passed: GCD Time =%s\n",get_time_str(tdiff));

	/* Now create a full-length y-vector = [large random padding integer]*[known small factor]: */
	memset(x,ONES64,(vec_len<<3));	x[vec_len-1] = (1ull << (p&63)) - 1;
	for(i = 0; i < vec_len-1; i++)
	{
		y[i] = rng_isaac_rand();
	}
	y[vec_len-1] = 0;

	y[vec_len-1] = mi64_mul_scalar(y, 87722769297534671ull, y, vec_len-1);
	ASSERT(HERE, y[vec_len-1] != 0, "GCD self-test #12: expected nonzero carryout of q0* y4[0:vec_len-1]!");

	clock1 = clock();
	gcd_debug=0;
	fail += (mi64_gcd(x, y, vec_len, eGCD, &Ap, &Bp, &len_AB, &sign_AB) != 1);
	fail += (x[0] !=    87722769297534671ull);
	ASSERT(HERE, fail==0,"GCD self-test #12b");

	clock2 = clock();
	tdiff = (double)(clock2 - clock1);
	printf	("GCD self-test #12 passed: GCD Time =%s\n",get_time_str(tdiff));
}
exit(0);
	/*
	Test #13: check the uint64[] <==> double[]  interconversion routines:
	on arrays of increasing size, doubling the size until reach MAX_ARRAY_DIM:
	*/
	vec_len = 1;
	while(vec_len <= MAX_ARRAY_DIM)
	{
		/* Init an array of (vec_len) uint64's using quasirandom selections from y4[]: */
		for(i = 0; i < vec_len; i++)
		{
			x[i] = rng_isaac_rand();
			y[i] = rng_isaac_rand();
			u[i] = x[i];
			v[i] = y[i];
		}

		ASSERT(HERE, mi64_cvt_uint64_double(x,y,   vec_len, a) == 4*vec_len, "GCD self-test #13: unexpected return value of mi64_cvt_uint64_double!");
		ASSERT(HERE, mi64_cvt_double_uint64(a, 4*vec_len, x,y) ==   vec_len, "GCD self-test #13: nonzero carryout from mi64_cvt_double_uint64!");
		for(i = 0; i < vec_len; i++)
		{
			ASSERT(HERE, x[i] == u[i], "GCD self-test #13: x[i] != u[i]");
			ASSERT(HERE, y[i] == v[i], "GCD self-test #13: y[i] != v[i]");
		}
		vec_len *= 2;
		/* Want to do at least one test right at the max. length: */
		if(vec_len > MAX_ARRAY_DIM)
		{
			if(vec_len >= 2*MAX_ARRAY_DIM)
				break;
			else
				vec_len = MAX_ARRAY_DIM;
		}
	}

#if 0	/****************************************************************************************/
	/*
	Test #14: test the FFT-mul routines with vectors of various sizes,
	starting with 4K 16-bit reals (= 2K 16-bit complex = 512 64-bit complex):
	*/
	vec_len = 512;
	ASSERT(HERE, MAX_ARRAY_DIM >= 2*vec_len, "GCD self-test #14: MAX_ARRAY_DIM too small!");
	ASSERT(HERE, !ARRAYS_OVERLAP(u,vec_len,v,vec_len), "GCD self-test #14: u/v-arrays overlap!");
	ASSERT(HERE, !ARRAYS_OVERLAP(x,vec_len,y,vec_len), "GCD self-test #14: x/y-arrays overlap!");
	while(vec_len <= MAX_ARRAY_DIM)
	{
		/* Init an array of (vec_len) uint64's using quasirandom selections from y4[]: */
		for(i = 0; i < vec_len; i++)
		{
			x[i] = rng_isaac_rand();
			y[i] = rng_isaac_rand();
			u[i] = x[i];
			v[i] = y[i];
		}
		mi64_mul_vector(u, vec_len, u, vec_len, w, &lenZ);
		ASSERT(HERE, lenZ == 2*vec_len, "GCD self-test #14: lenZ != 2*vec_len");
		mi64_mul_vector(v, vec_len, v, vec_len, z, &lenZ);
		ASSERT(HERE, lenZ == 2*vec_len, "GCD self-test #14: lenZ != 2*vec_len");

		/* Convert x+I*y to packed-double form, store outputs in A: */
		fft_len = mi64_cvt_uint64_double(x,y,   vec_len, a);
		ASSERT(HERE, fft_len == 4*vec_len, "GCD self-test #14: unexpected return value of mi64_cvt_uint64_double!");
		fft_len *= 2;	/* Multiply by 2 here because mi64_cvt_uint64_double returns *complex* FFT length
						that results from odd/even interleaving and int64-to-real conversion of the 2 input arrays
						*/
		ASSERT(HERE, !ARRAYS_OVERLAP(a,2*fft_len,b,2*fft_len), "GCD self-test #14: a/b-arrays overlap!");
		/* Zero-pad the input vector in preparation for FFT-based multiply: */
		memset(a+fft_len, 0, (fft_len<<3));
		/* Call the init-FFT routines, using the (as yet uninited) b-array for the needed scratch space: */
		pairFFT_mul(  b, 0x0, 2*fft_len, TRUE , FALSE);	/* INIT_ARRAYS = TRUE , FORWARD_FFT_ONLY = FALSE */
		/* Save a copy of the input vector: */
		memcpy(b,a,(fft_len<<4));	// double the length of memcpy, thus no 2nd memset(0) step needed to 0-pad b.
		/* Forward FFT of A-vector: */
		pairFFT_mul(  a,   0, 2*fft_len, FALSE, TRUE );	/* INIT_ARRAYS = FALSE, FORWARD_FFT_ONLY = TRUE  */

		/* Now do the FFT-based squaring: */
		pairFFT_mul(  b,   0, 2*fft_len, FALSE, FALSE);	/* INIT_ARRAYS = FALSE, FORWARD_FFT_ONLY = FALSE */

		/* Convert result back to uint64[] form and compare 2 outputs to all-integer squarings: */
		ASSERT(HERE, (j = mi64_cvt_double_uint64(b, fft_len, x,y)) <= 2*vec_len, "GCD self-test #14: nonzero carryout from mi64_cvt_double_uint64!");

/*		ASSERT(HERE, mi64_cmp_eq(x,w, 2*vec_len), "GCD self-test #14: x*x != u*u");	*/
		for(i = 0; i < 2*vec_len; i++)
		{
			ASSERT(HERE, x[i] == w[i], "GCD self-test #14: x*x[i] != u*u[i]");
		}

/*		ASSERT(HERE, mi64_cmp_eq(x,w, 2*vec_len), "GCD self-test #14: x*x != u*u");	*/
		for(i = 0; i < 2*vec_len; i++)
		{
			ASSERT(HERE, y[i] == z[i], "GCD self-test #14: x*x[i] != u*u[i]");
		}

		vec_len *= 2;

		if(vec_len > MAX_ARRAY_DIM)
		{
			break;
		}
	}
#endif/****************************************************************************************/

	return fail;
}

