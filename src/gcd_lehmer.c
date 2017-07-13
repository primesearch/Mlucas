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

/*
If the FFT-mul enabled, need all the stuff for factor.c + int-gcd,
including factor.c built with -DINCLUDE_PM1, plus the following FFT stuff:

g42 -c -Wall -g3 -ggdb -DFACTOR_STANDALONE -DTRYQ=4 -DNWORD -DINCLUDE_PM1 factor.c
g42 -c -Wall -O3 br.c dft*.c get_fft*c *pairFFT*.c
g42 -c -Wall -O3 radix32_dif_dit_pass.c
g42 -c -Wall -O3 radix16_dif_dit_pass.c
g42 -c -Wall -O3 radix8_dif_dit_pass.c
g42 -c -Wall -O3 radix32_*cy*.c
g42 -c -Wall -O3 radix16_*cy*.c
g42 -c -Wall -O3 radix8_*cy*.c
*/

#include "gcd_lehmer.h"

// Use x86_64 inline-asm?
#undef YES_ASM
#if(defined(CPU_IS_X86_64) && defined(COMPILER_TYPE_GCC) && (OS_BITS == 64))
	#define YES_ASM
#endif

#define GCD_DEBUG	1
/*
   Set == 1 to enable basic per-pass diagnostic printing.
   Set >= 2 to enable detailed per-pass diagnostic printing.
			WARNING: level-2 diagnostics not recommended for large vectors!
*/
	int fft_gcd_debug = 0;
	FILE *fp;
	static char *file_access_mode[2] = {"a","w"};
	char string0[STR_MAX_LEN];
#if GCD_DEBUG >= 1
	// Set GCD_DEBUG >= 1 and this != 0 to selectively enable debug
	// printing for some particular self-test or phase of the computation:
	int gcd_debug = 0;
	char string1[STR_MAX_LEN];
	char str_10k[10*1024];
#endif
#if GCD_DEBUG >= 2
	char string2[STR_MAX_LEN];
#endif

#define FFTMUL_THRESHOLD_BITS	16384	/* # of bits in a multiword base-2^64 integer
										at which to switch from pure-int grammar-school
										to floating-FFT-based multiply algorithm */
/***********************/

// Returns: number of 64-bit words in computed abcd-multipliers.
uint32 fft_gcd_get_mults(
	const uint64 u[], const uint64 v[], const uint32 vec_len, const uint32 targ_len, const uint32 sec_len,
	uint64 ai[], uint64 bi[], uint64 ci[], uint64 di[], uint64 w[], uint64 x[], uint64 y[], uint64 z[],
	double a[], double b[], double c[], double d[], uint32*len_ab, uint32*len_cd, uint32*sign_ab, uint32*sign_cd)
{
	uint32 retval = 0;
	uint32 i,j,k,kk,cy,cy_fwd, eGCD = FALSE, HALF = FALSE;
	uint32 fft_len = sec_len<<4, pad_len = fft_len + ( (fft_len >> DAT_BITS) << PAD_BITS );	// Padded-array length
	 int64 cy_re,cy_im;	// Prefer these to be signed
	uint64 tmp64 = 0,cy1,cy2;
	const double sign_mult[2] = {1.0,-1.0};
	const int check_len = 4;	// #high 64-bit words to compute for hi-block cancellation check
	uint64 tvec0[check_len],tvec1[check_len],tvec2[check_len],tvec3[check_len];

		j = 2*sec_len+10;	// Adding 10 guard words to half-exit EGCD inputs to ensure accuracy of au-bv,cu-dv signs
		// Copy leading ? bits of input-vecs into a pair of uint64-arrays. 16 bits per double ==> 64-bit vec,
		// but need 2x as many bits in leadU,V as we will use in ensuing FFT-muls, hence the 2*sec_len:
		mi64_set_eq(w,u+vec_len-j,j);	// w = high j words of u
		mi64_set_eq(z,v+vec_len-j,j);	// z = high j words of v
		// Return value is length of abcd-multipliers:
		retval = mi64_gcd(w,z, j, eGCD=TRUE,ai,bi,len_ab,sign_ab, HALF=TRUE,ci,di,len_cd,sign_cd, sec_len);
		ASSERT(HERE, retval == sec_len, "Unexpected retval from mi64_gcd call to get abcd-multipliers!");
	#if 1
		if(!*sign_ab) { printf("a*u-b*v > 0\n"); } else if(*sign_ab == 1) { printf("a*u-b*v < 0\n"); } else { ASSERT(HERE, 0, "foo!"); }
		if(!*sign_cd) { printf("c*u-d*v > 0\n"); } else if(*sign_cd == 1) { printf("c*u-d*v < 0\n"); } else { ASSERT(HERE, 0, "bar!"); }
	#endif
		// Verify that resulting abcd-multipliers have at most (sec_len) set words:
		i = sec_len;
		ASSERT(HERE, ai[i] == 0ull && bi[i] == 0ull && ci[i] == 0ull && di[i] == 0ull, "abcd-multipliers exceed expected length!");
		// Convert a+I*b and c+I*d to packed-double form, store outputs in double-arrays A,B, resp.:
		cy = mi64_cvt_uint64_double(ai,bi,0,sec_len, a);
		ASSERT(HERE, cy == 0, "FFT_gcd: unexpected carryout of mi64_cvt_uint64_double!");
		cy = mi64_cvt_uint64_double(ci,di,0,sec_len, b);
		ASSERT(HERE, cy == 0, "FFT_gcd: unexpected carryout of mi64_cvt_uint64_double!");
		// Zero upper halves of FFT-multiplier vectors, that is the upper pad_len/2 doubles:
		memset(a+(pad_len>>1), 0, (pad_len<<2));
		memset(b+(pad_len>>1), 0, (pad_len<<2));

		// Forward FFT of A/B-vector [in a] and C/D vector [in b]:
		pairFFT_mul(  a,b,0, fft_len, FALSE, TRUE );	// INIT_ARRAYS = FALSE, FORWARD_FFT_ONLY = TRUE

		/************ debug-dump egcd-mults: *****************/
		if(0) {
			fp = fopen("c.txt",file_access_mode[vec_len > 16000]);	ASSERT(HERE, fp != 0x0, "Null file pointer!");
			fprintf(fp,"vec_len = %u\n",vec_len);
			fprintf(fp,"Leading u,v-elts:\n");
			for(i = vec_len; i > vec_len-10; i--) {
				fprintf(fp,"%5u\t%20llu\t%20llu\n",i,u[i],v[i]);
			}
			fprintf(fp,"egcd-mults:\n");
			for(i = 0; i < sec_len; i++) {
				fprintf(fp,"%5u\t%20llu\t%20llu\t%20llu\t%20llu\n",i,ai[i],bi[i],ci[i],di[i]);
			}
			fprintf(fp,"fFFt[a,b]:\n");
			for(i = 0; i < fft_len; i+=2) {
				j = i + ( (i >> DAT_BITS) << PAD_BITS );
				fprintf(fp,"%5u\t%12.3f\t%12.3f\t%12.3f\t%12.3f\n",i,a[j],a[j+1],b[j],b[j+1]);
			}
			fclose(fp); fp = 0x0;
			exit(0);
		}
		/*****************************************************/

		// [a*u, b*v], [c*u, d*v] separately computed for now:

		// Input vecs to FFT-gcd must be suitably 0-padded at top:
		i = vec_len;	j = sec_len - (i & (sec_len-1));	// mod-via-AND here assumes sec_len a power of 2
		ASSERT(HERE, mi64_iszero(u+i,j) && mi64_iszero(v+i,j), "Input vecs to FFT-gcd must be suitably 0-padded at top!");

	/********* 3rd-lower block beneath the highest one - this has the vector-carryins in the hi result half: **********/
		kk = vec_len-4*sec_len;	// Can try adding a few guard words to our block-start indices
		// Only really need to clear hi half since lo half overwritten by mi64_cvt_uint64_double call:
		memset(c, 0, (pad_len<<3));
		// Process length = sec_len section starting at index [kk]:
		printf("Processing length-%u section starting at index %u:\n",sec_len,kk);
		cy_fwd = mi64_cvt_uint64_double(u+kk,v+kk, 0, sec_len, c);	printf("cy_fwd = %d\n",cy_fwd);
		// Forward FFT of U/V-vector section [in c]:
		pairFFT_mul(  c,0,0, fft_len, FALSE, TRUE );	// INIT_ARRAYS = FALSE, FORWARD_FFT_ONLY = TRUE
		// d = copy of 0-padded u/v data stored in fFFTed form in c:
		memcpy(d, c, (pad_len<<3));
		// Do the [a*u, b*v] FFT-mul; the subtract-and-store-result-in-Re-part occurs in the carry step:
		pairFFT_mul(  c,a,0, fft_len, FALSE, 2*TRUE);	// INIT_ARRAYS = FALSE, FORWARD_FFT_ONLY = 2 (indicates both inputs already fFFTed)
		// Do the [c*u, d*v] FFT-mul:
		pairFFT_mul(  d,b,0, fft_len, FALSE, 2*TRUE);	// INIT_ARRAYS = FALSE, FORWARD_FFT_ONLY = 2
		// Interleave the [a*u-b*v, c*u-d*v] outputs; c*u-d*v go into odd slots of C-array:
		for(i = 0; i < fft_len; i+=2) {
			j = i + ( (i >> DAT_BITS) << PAD_BITS );
			c[j] *= sign_mult[*sign_ab];	c[j+1] = sign_mult[*sign_cd]*d[j];
		}
		cy = mi64_cvt_double_uint64(c,fft_len, x,y);

	/********* 2nd-lower block beneath the highest one - this has the vector-carryins in the hi result half: **********/
		kk = vec_len-3*sec_len;	// Can try adding a few guard words to our block-start indices
		// Only really need to clear hi half since lo half overwritten by mi64_cvt_uint64_double call:
		memset(c, 0, (pad_len<<3));
		// Process length = sec_len section starting at index [kk]:
		printf("Processing length-%u section starting at index %u:\n",sec_len,kk);
		cy_fwd = mi64_cvt_uint64_double(u+kk,v+kk, cy_fwd, sec_len, c);	printf("cy_fwd = %d\n",cy_fwd);
		// Forward FFT of U/V-vector section [in c]:
		pairFFT_mul(  c,0,0, fft_len, FALSE, TRUE );	// INIT_ARRAYS = FALSE, FORWARD_FFT_ONLY = TRUE
		// d = copy of 0-padded u/v data stored in fFFTed form in c:
		memcpy(d, c, (pad_len<<3));
		// Do the [a*u, b*v] FFT-mul; the subtract-and-store-result-in-Re-part occurs in the carry step:
		pairFFT_mul(  c,a,0, fft_len, FALSE, 2*TRUE);	// INIT_ARRAYS = FALSE, FORWARD_FFT_ONLY = 2 (indicates both inputs already fFFTed)
		// Do the [c*u, d*v] FFT-mul:
		pairFFT_mul(  d,b,0, fft_len, FALSE, 2*TRUE);	// INIT_ARRAYS = FALSE, FORWARD_FFT_ONLY = 2
		// Interleave the [a*u-b*v, c*u-d*v] outputs; c*u-d*v go into odd slots of C-array:
		for(i = 0; i < fft_len; i+=2) {
			j = i + ( (i >> DAT_BITS) << PAD_BITS );
			c[j] *= sign_mult[*sign_ab];	c[j+1] = sign_mult[*sign_cd]*d[j];
		}
		cy = mi64_cvt_double_uint64(c,fft_len, w,z);
		cy_re = mi64_add(w,x+sec_len,w, sec_len);
		cy_im = mi64_add(z,y+sec_len,z, sec_len);
		printf("\tAfter previous-block vector-carry add: cy_re,im = %lld,%lld ... ignoring.\n",cy_re,cy_im);

	/********* Penultimate block: *************************************/
		kk = vec_len-2*sec_len;
		// Only really need to clear hi half since lo half overwritten by mi64_cvt_uint64_double call:
		memset(c, 0, (pad_len<<3));
		// Process length = sec_len leading section:
		printf("Processing length-%u section starting at index %u, cy_fwd = %u:\n",sec_len,kk,cy_fwd);
		cy_fwd = mi64_cvt_uint64_double(u+kk,v+kk, cy_fwd, sec_len, c);	//ASSERT(HERE, cy_fwd == 0, "Unexpected carryout!");***if this were uppermost block
		cy1 = cy_fwd & 1ull; cy2 = (cy_fwd>>1) & 1ull;	// Save these for later low-64-bit check
		// Forward FFT of U/V-vector section [in c]:
		pairFFT_mul(  c,0,0, fft_len, FALSE, TRUE );	// INIT_ARRAYS = FALSE, FORWARD_FFT_ONLY = TRUE
		// d = copy of 0-padded u/v data stored in fFFTed form in c:
		memcpy(d, c, (pad_len<<3));
		// Do the [a*u, b*v] FFT-mul; the subtract-and-store-result-in-Re-part occurs in the carry step:
		pairFFT_mul(  c,a,0, fft_len, FALSE, 2*TRUE);	// INIT_ARRAYS = FALSE, FORWARD_FFT_ONLY = 2 (indicates both inputs already fFFTed)
		// Do the [c*u, d*v] FFT-mul:
		pairFFT_mul(  d,b,0, fft_len, FALSE, 2*TRUE);	// INIT_ARRAYS = FALSE, FORWARD_FFT_ONLY = 2
		// Interleave the [a*u-b*v, c*u-d*v] outputs; c*u-d*v go into odd slots of C-array:
		for(i = 0; i < fft_len; i+=2) {
			j = i + ( (i >> DAT_BITS) << PAD_BITS );
			c[j] *= sign_mult[*sign_ab];	c[j+1] = sign_mult[*sign_cd]*d[j];
		}
		// Convert c/d to uint64, yielding a pair of length-(2*sec_len) vectors x,y plus borrow-outs, expected = 0 here:
		cy = mi64_cvt_double_uint64(c,fft_len, x,y);	//ASSERT(HERE, cy==0, "Unexpected carryout!");***if this were uppermost block
		cy_re = mi64_add(x,w+sec_len,w+sec_len, sec_len);
		cy_im = mi64_add(y,z+sec_len,z+sec_len, sec_len);
		printf("\tAfter previous-block vector-carry add: cy_re,im = %lld,%lld ... ignoring.\n",cy_re,cy_im);

	/********* Uppermost block: *************************************/
		// Expect result of same FFT-mul-computed [a*u-b*v, c*u-d*v] for uppermost block to exactly
		// cancel upper half of penultimate-block lin-combo output, except possibly for low 64 bits,
		// which we cheaply compute via integer arithmetic:
		kk = vec_len-sec_len;
		mi64_add_scalar(u+kk,cy1, tvec0,check_len);
		mi64_add_scalar(v+kk,cy2, tvec1,check_len);
	// Compute w[2*sec_len+i] = x[sec_len] + SGN(ai[i]*(u[kk+i]+cy1) - bi[i]*(v[kk+i]+cy2), sign_ab) :
		mi64_mul_vector_lo_half(ai,tvec0, tvec2,check_len);	// ai[i]*(u[kk+i]+cy1)
		mi64_mul_vector_lo_half(bi,tvec1, tvec3,check_len);	// bi[i]*(v[kk+i]+cy2)
		if(*sign_ab)
			mi64_sub(tvec3,tvec2, tvec2,check_len);	// SGN(ai[i]*(u[kk+i]+cy1) - bi[i]*(v[kk+i]+cy2), sgn_ab)
		else
			mi64_sub(tvec2,tvec3, tvec2,check_len);
		mi64_add(x+sec_len,tvec2, w+2*sec_len,check_len);
	// Compute z[2*sec_len+i] = y[sec_len] + SGN(ci[i]*(u[kk+i]+cy1) - di[i]*(v[kk+i]+cy2), sign_cd) :
		mi64_mul_vector_lo_half(ci,tvec0, tvec2,check_len);	// ci[i]*(u[kk+i]+cy1)
		mi64_mul_vector_lo_half(di,tvec1, tvec3,check_len);	// di[i]*(v[kk+i]+cy2)
		if(*sign_cd)
			mi64_sub(tvec3,tvec2, tvec2,check_len);	// SGN(ci[i]*(u[kk+i]+cy1) - di[i]*(v[kk+i]+cy2), sgn_cd)
		else
			mi64_sub(tvec2,tvec3, tvec2,check_len);
		mi64_add(y+sec_len,tvec2, z+2*sec_len,check_len);
	/* Single-low-word version of above:
		w[2*sec_len] = x[sec_len] + SGN(ai[0]*(u[kk]+cy1) - bi[0]*(v[kk]+cy2), sign_ab);
		z[2*sec_len] = y[sec_len] + SGN(ci[0]*(u[kk]+cy1) - di[0]*(v[kk]+cy2), sign_cd);
	*/
		// At most the lowest 'one-beyond' result word should be nonzero:
		ASSERT(HERE, mi64_iszero(w+2*sec_len+1,check_len-1) && mi64_iszero(z+2*sec_len+1,check_len-1),"At most the lowest 'one-beyond' result word should be nonzero!");
	if(0) {
		fp = fopen("c.txt","w");	ASSERT(HERE, fp != 0x0, "Null file pointer!");
		fprintf(fp,"Uppermost block of %u FFT data:\n",fft_len);
		fprintf(fp,"\t\ta*u - b*v\t\tc*u - d*v:\n");
	//	for(i = 0; i < fft_len; i+=2) {
		for(i = 0; i < 2*sec_len+1; i++) {
		//	j = i + ( (i >> DAT_BITS) << PAD_BITS );
		//	fprintf(fp,"%5u\t%12.3f\t%12.3f\n",i,c[j],c[j+1]);
			fprintf(fp,"%5u\t%20llu\t%20llu\n",i+vec_len-3*sec_len,w[i],z[i]);
		}
		fclose(fp); fp = 0x0;
		exit(0);
	}

/*
Now that have 2*sec_len[+1] leading good words of the full-length u',v' result vectors, compute a fresh
length = sec_len set of abcd-linear-combination-reduction multipliers for those:

...And now compute the product of the two 2x2 multiplier arrays:

	/aa bb\     /a0 b0\     /a1 b1\     /a0*a1 + b0*c1  a0*b1 + b0*d1\
   |       | = |       | * |       | = |                              |
	\cc dd/     \c0 d0/     \c1 d1/     \c0*a1 + d0*c1  c0*b1 + d0*d1/ .

To do this, we need fwd-FFTs of the vector pairs [a0,b0], [c0,d0], [a1,c1], [b1,d1].
We already have the first 2 pairwise-fwd-FFTs of [a0,b0] and [c0,d0] in our floating-point a,b-vectors.
*/
		j = 2*sec_len+1;
		i = mi64_getlen(w, j);	j = mi64_getlen(v, j);	j = MAX(i,j);
		// Return value is length of abcd-multipliers:
		j = mi64_gcd(w,z, j, eGCD=TRUE,ai,bi,len_ab,sign_ab, HALF=TRUE,ci,di,len_cd,sign_cd, sec_len);
		ASSERT(HERE, j == sec_len, "Unexpected retval from mi64_gcd call to get abcd-multipliers!");
	#if 1
		if(!*sign_ab) { printf("a*u-b*v > 0\n"); } else if(*sign_ab == 1) { printf("a*u-b*v < 0\n"); } else { ASSERT(HERE, 0, "foo!"); }
		if(!*sign_cd) { printf("c*u-d*v > 0\n"); } else if(*sign_cd == 1) { printf("c*u-d*v < 0\n"); } else { ASSERT(HERE, 0, "bar!"); }
	#endif
		// Verify that resulting abcd-multipliers have at most (sec_len) set words:
		i = sec_len;
		ASSERT(HERE, ai[i] == 0ull && bi[i] == 0ull && ci[i] == 0ull && di[i] == 0ull, "abcd-multipliers exceed expected length!");
		// Convert a1+I*c1 and b1+I*d1 - note b/c swap vs abcd0-computation - to packed-double form:
		cy = mi64_cvt_uint64_double(ai,ci,0,sec_len, c);
		ASSERT(HERE, cy == 0, "FFT_gcd: unexpected carryout of mi64_cvt_uint64_double!");
		cy = mi64_cvt_uint64_double(bi,di,0,sec_len, d);
		ASSERT(HERE, cy == 0, "FFT_gcd: unexpected carryout of mi64_cvt_uint64_double!");
		// Zero upper halves of FFT-multiplier vectors, that is the upper pad_len/2 doubles:
		memset(c+(pad_len>>1), 0, (pad_len<<2));
		memset(d+(pad_len>>1), 0, (pad_len<<2));
		// Because the genfftmul_carry_norm_pow2_errcheck() macro is set up to compute difference of Re,Im=parts
		// of FFT=mul outputs (as needed for a*u-b*v, c*u-d*v GCD-reduction muls), need to negate Im-parts here
		// to -c1 and -d1 to yield the desired RHS sums in the above-described 2x2 matrix product:
		for(i = 0; i < pad_len; i+=2) {	// Save ourselves indexing work by negating both data & pads here
			c[i+1] = -c[i+1];	d[i+1] = -d[i+1];
		}
		// Forward FFT of a1/-c1-vector [in c] and b1/-d1 vector [in d]:
		pairFFT_mul(  c,d,0, fft_len, FALSE, TRUE );	// INIT_ARRAYS = FALSE, FORWARD_FFT_ONLY = TRUE
		// Stash a copy of fFFT(c) in next-higher [pad_len]-sized section:
		memcpy(c+pad_len, c, (pad_len<<3));
		// Do the [a0*a1,-b0*c1] FFT-mul; the subtract-and-store-[a0*a1 + b0*c1]-result-in-Re-part occurs in the carry step:
		pairFFT_mul(  c,a,0, fft_len, FALSE, 2*TRUE);	// INIT_ARRAYS = FALSE, FORWARD_FFT_ONLY = 2 (indicates both inputs already fFFTed)
		// Do the [a0*b1,-b0*d1] FFT-mul, yielding bb = [a0*b1 + b0*d1]; DONE WITH A-ARRAY AFTERWARD, thus store result in it:
		pairFFT_mul(  a,d,0, fft_len, FALSE, 2*TRUE);	// INIT_ARRAYS = FALSE, FORWARD_FFT_ONLY = 2
		// Interleave the [aa,bb] = [a0*a1 + b0*c1, a0*b1 + b0*d1] outputs, put into A-array:
		for(i = 0; i < pad_len; i+=2) {
			// Again do pads and data both; Must move the bb-data from odd to even A-slots first:
			a[i] = a[i+1];	a[i+1] = c[i];
		}
		// C-array now free again, but already stashed copy of fFFT(c) in next-higher [pad_len]-sized section,
		// so just do FFT-muls directly with that:
		// Do the [b0*a1,-d0*c1] FFT-mul, yielding cc = [b0*a1 + d0*c1]:
		pairFFT_mul(c+pad_len,b,0, fft_len, FALSE, 2*TRUE);	// INIT_ARRAYS = FALSE, FORWARD_FFT_ONLY = 2 (indicates both inputs already fFFTed)
		// Do the [b0*b1,-d0*d1] FFT-mul, yielding dd = [a0*b1 + b0*d1]; DONE WITH B-ARRAY AFTERWARD, thus store result in it:
		pairFFT_mul(  b,d,0, fft_len, FALSE, 2*TRUE);	// INIT_ARRAYS = FALSE, FORWARD_FFT_ONLY = 2
		// Interleave the [cc,dd] = [b0*a1 + d0*c1, a0*b1 + b0*d1] outputs, put into B-array:
		for(i = 0; i < pad_len; i+=2) {
			// Again do pads and data both; Must move the dd-data from odd to even B-slots first:
			b[i] = b[i+1];	b[i+1] = d[i];
		}
#warning to-do! *** cvt a,b back to int and use those 2*sec_len - long ints to do int-mul to test resulting gcd-reduction ***
		// Lastly, need to 0-pad the now-doubled-length mult-arrays and compute fFFTs of those:
		memset(a+pad_len, 0, (pad_len<<3));
		memset(b+pad_len, 0, (pad_len<<3));
		fft_len <<= 1;	pad_len <<= 1;	// Now double both FFT-len params
		pairFFT_mul(  a,b,0, fft_len, FALSE, TRUE );	// INIT_ARRAYS = FALSE, FORWARD_FFT_ONLY = TRUE

	retval += sec_len;	// = 2*sec_len
	return retval;	// Will return actual multiplier length, once have code to compute it
}

/***********************/

/*
FFT-mul-accelerated Lehmer-type Euclidean GCD(U, V), where U,V are multiword ints stored
in base-2^64 positive-digit form; U and V are assumed to both have at least (ndim)
allocated 64-bit words, with U[0],V[0] being the least-significant words.

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

CALLS:	mi64_gcd, both [1] in EGCD+HALF mode to generate the smallest FFT-sized
		multiword ABCD multipliers, and [2] for the final GCD computation once the operands
		have been reduced in size sufficiently for FFT-mul to lose its effectiveness.
*/
uint32	fft_gcd(
	uint64 u[], uint64 v[], const uint32 ndim,
	const uint32 EGCD, uint64 Ap[], uint64 Bp[], uint32 *len_AB, uint32 *sign_AB)
{
	uint32 i,j,kk,bw,cy,cy_fwd, vec_len, eGCD = FALSE, HALF = FALSE, len_ab,sign_ab, len_cd,sign_cd;
	static uint32 sec_len, fft_len,pad_len, first_entry=TRUE;
	 int64 cy_re,cy_im;	// Prefer these to be signed
	uint64 tmp64 = 0,cy1,cy2;
	static uint64 *w = 0x0, *x = 0x0, *y = 0x0, *z = 0x0, *ai = 0x0, *bi = 0x0, *ci = 0x0, *di = 0x0;
	static double *a = 0x0, *b = 0x0, *c = 0x0, *d = 0x0;
	const double sign_mult[2] = {1.0,-1.0};
	ASSERT(HERE, IS_ODD(u[0]) && IS_ODD(v[0]), "Even inputs not yet supported!");

	if(first_entry) {
		fprintf(stderr,"fft_gcd alloc with ndim = %u\n",ndim);
		ASSERT(HERE, w == 0x0, "alloc-related sentinel pointer non-null!");
		first_entry=FALSE;
		/* Allocate the main data arrays: */
		sec_len = FFTMUL_THRESHOLD_BITS>>2;
	i = ndim + 2*sec_len;	//*** make large enough to hold mi64_mul_vector results if debugging-to-file
//	i = sec_len<<1;
		x  = (uint64 *)calloc(i, sizeof(uint64));
		y  = (uint64 *)calloc(i, sizeof(uint64));
		w  = (uint64 *)calloc(i, sizeof(uint64));
		z  = (uint64 *)calloc(i, sizeof(uint64));
		ai = (uint64 *)calloc(i, sizeof(uint64));
		bi = (uint64 *)calloc(i, sizeof(uint64));
		ci = (uint64 *)calloc(i, sizeof(uint64));
		di = (uint64 *)calloc(i, sizeof(uint64));
		// For now just use 16 bits per double:
		fft_len = 4*sec_len;
		fft_len *= 2;	/* Extra multiply by 2 here because mi64_cvt_uint64_double returns *complex* FFT length
						that results from odd/even interleaving and int64-to-real conversion of the 2 input arrays.
						Thus, we now have the real-vector FFT length, BUT NOT YET THE ZERO-PADDED LENGTH. */
		pad_len = fft_len + ( (fft_len >> DAT_BITS) << PAD_BITS );	/* Padded-array length */
		fft_len *= 2;	pad_len *= 2;	// Both of these now reflect real-vector FFT length *with* 0-padding
		fprintf(stderr,"fft_gcd double-array allocs with pad_len = %u\n",pad_len);
		a = (double *)calloc(  pad_len, sizeof(double));
		b = (double *)calloc(  pad_len, sizeof(double));
		c = (double *)calloc(2*pad_len, sizeof(double));	fprintf(stderr,"Temporary: C-array gets 2x alloc to allow high-half to be used for scratch storage.\n");
		d = (double *)calloc(  pad_len, sizeof(double));
		ASSERT(HERE, !ARRAYS_OVERLAP(a,pad_len, b,pad_len), "FFT_gcd: a/b-arrays overlap!");
		// Call the init-FFT routines, using the (currently unused) c,d-arrays for the needed scratch space:
		pairFFT_mul(  c,d,0, fft_len, TRUE , FALSE);	// INIT_ARRAYS = TRUE , FORWARD_FFT_ONLY = FALSE
	}

	/* Test for EGCD: */
	if(EGCD) ASSERT(HERE, Ap != 0ull && Bp != 0ull && len_AB != 0ull && sign_AB != 0ull, "One or more of the EGCD-required return-pointers null!");

	i = mi64_getlen(u, ndim);	j = mi64_getlen(v, ndim);	vec_len = MAX(i,j);
	while(vec_len >= 2*sec_len)
	{
		// Input vecs to FFT-gcd must be suitably 0-padded at top:
		i = vec_len;	j = sec_len - (i & (sec_len-1));	// mod-via-AND here assumes sec_len a power of 2
		mi64_clear(u+i,j);	mi64_clear(v+i,j);
#warning May need to clear more words for half-length reduction!

		j = 2*sec_len+10;	// Adding 10 guard words to half-exit EGCD inputs to ensure accuracy of au-bv,cu-dv signs
		// Copy leading ? bits of input-vecs into a pair of uint64-arrays. 16 bits per double ==> 64-bit vec,
		// but need 2x as many bits in leadU,V as we will use in ensuing FFT-muls, hence the 2*sec_len:
		mi64_set_eq(w,u+vec_len-j,j);	// w = high j words of u
		mi64_set_eq(z,v+vec_len-j,j);	// z = high j words of v
		// Call recursive routine to obtain vector FFT-multipliers for up-to-half-length reduction of [u,v] vectors:
		// Return value is length of abcd-multipliers:
		uint32 targ_len = vec_len>>1;
		fft_gcd_get_mults(u,v,vec_len, targ_len,sec_len, ai,bi,ci,di, w,x,y,z, a,b,c,d, &len_ab,&len_cd,&sign_ab,&sign_cd);

		/************ debug-dump egcd-mults: *****************/
		if(0) {
			fp = fopen("a.txt",file_access_mode[vec_len > 16000]);	ASSERT(HERE, fp != 0x0, "Null file pointer!");
			fprintf(fp,"vec_len = %u\n",vec_len);
			fprintf(fp,"Leading u,v-elts:\n");
			for(i = vec_len; i > vec_len-10; i--) {
				fprintf(fp,"%5u\t%20llu\t%20llu\n",i,u[i],v[i]);
			}
			fprintf(fp,"egcd-mults:\n");
			for(i = 0; i < sec_len; i++) {
				fprintf(fp,"%5u\t%20llu\t%20llu\t%20llu\t%20llu\n",i,ai[i],bi[i],ci[i],di[i]);
			}
			fprintf(fp,"fFFt[a,b]:\n");
			for(i = 0; i < fft_len; i+=2) {
				j = i + ( (i >> DAT_BITS) << PAD_BITS );
				fprintf(fp,"%5u\t%12.3f\t%12.3f\t%12.3f\t%12.3f\n",i,a[j],a[j+1],b[j],b[j+1]);
			}
			fclose(fp); fp = 0x0;
			exit(0);
		}
		/*****************************************************/

		// [a*u, b*v], [c*u, d*v] separately computed for now:

		// Convert abcd-multiplier-sized sections of u+I*v to packed-double form, result ==> double-array C.
		// Unlike the above two mi64_cvt_uint64_double calls, we expect the possibility of
		// carryouts here because we are [uint64 -> balanced-digit-double] converting sections
		// of a pair of larger vectors, of which only the most-significant such section must have MSW > 0:

		/**************** pure-int a*u - b*v: ****************/
		if(fft_gcd_debug) {
			fp = fopen("a.txt",file_access_mode[vec_len > 16000]);	ASSERT(HERE, fp != 0x0, "Null file pointer!");
		// Compute a*u - b*v:
			i = vec_len;	j = sec_len;
			mi64_clear(x,i+j);	mi64_clear(w,i+j);	mi64_clear(z,i+j);	// w,z will hold intermediate products a*u, b*v
			mi64_mul_vector(u,i, ai,j, w,&kk);	fprintf(fp,"A*U has %u words = %20llu,%20llu,...\n",kk,w[kk-1],w[kk-2]);
			mi64_mul_vector(v,i, bi,j, z,&kk);	fprintf(fp,"B*V has %u words = %20llu,%20llu,...\n",kk,z[kk-1],z[kk-2]);
			for(i = kk-1; (int)i > 0; i--)
				if(w[i] != z[i]) break;
			for(j = i-1; (int)j > 0; j--)
				if(w[j] != (uint64)-1 || z[j] != 0ull) break;
			fprintf(fp,"First diff at word %5u, next at %5u:\n",i,j);
			fprintf(fp,"A*U[diff1-2] = %20llu,%20llu,%20llu,...,%20llu\n",w[i],w[i-1],w[i-2], w[j]);
			fprintf(fp,"B*V[diff1-2] = %20llu,%20llu,%20llu,...,%20llu\n",z[i],z[i-1],z[i-2], z[j]);
			// |a*u - b*v| result into x:
		if(!sign_ab) { tmp64 = mi64_sub(w,z, x, vec_len+sec_len); } else { tmp64 = mi64_sub(z,w, x, vec_len+sec_len); }
			fprintf(fp,"|A*U - B*V| has borrow = %lld:\n",tmp64);	ASSERT(HERE, tmp64 == 0ull,"|A*U - B*V| has nonzero borrow!");
		// Compute c*u - d*v:
			i = vec_len;	j = sec_len;
			mi64_clear(y,i+j);	mi64_clear(w,i+j);	mi64_clear(z,i+j);	// w,z will hold intermediate products d*u, d*v
			mi64_mul_vector(u,i, ci,j, w,&kk);	fprintf(fp,"C*U has %u words = %20llu,%20llu,...\n",kk,w[kk-1],w[kk-2]);
			mi64_mul_vector(v,i, di,j, z,&kk);	fprintf(fp,"D*V has %u words = %20llu,%20llu,...\n",kk,z[kk-1],z[kk-2]);
			for(i = kk-1; (int)i > 0; i--)
				if(w[i] != z[i]) break;
			for(j = i-1; (int)j > 0; j--)
				if(w[j] != (uint64)-1 || z[j] != 0ull) break;
			fprintf(fp,"First diff at word %5u, next at %5u:\n",i,j);
			fprintf(fp,"C*U[diff1-2] = %20llu,%20llu,%20llu,...,%20llu\n",w[i],w[i-1],w[i-2], w[j]);
			fprintf(fp,"D*V[diff1-2] = %20llu,%20llu,%20llu,...,%20llu\n",z[i],z[i-1],z[i-2], z[j]);
			// c*u - d*v result into y:
		if(!sign_cd) { tmp64 = mi64_sub(w,z, y, vec_len+sec_len); } else { tmp64 = mi64_sub(z,w, y, vec_len+sec_len); }
			fprintf(fp,"|C*U - D*V| has borrow = %lld:\n",tmp64);	ASSERT(HERE, tmp64 == 0ull,"|C*U - D*V| has nonzero borrow!");
		
			i = mi64_getlen(x, vec_len+sec_len);	fprintf(fp,"A*U - B*V has %u nonzero words\n",i);
			j = mi64_getlen(y, vec_len+sec_len);	fprintf(fp,"C*U - D*V has %u nonzero words\n",j);
			kk = MAX(i,j);
			fprintf(fp,"INT-mul [a*u - b*v, c*u - d*v] for vec_len = %u:\n",vec_len);
			for(i = 0; i < kk; i++) {
				fprintf(fp,"%5u\t%20llu\t%20llu\n",i,x[i],y[i]);
			}
			fclose(fp); fp = 0x0;
		exit(0);
		}
		/*****************************************************/

		cy_fwd = 0;		// This cy-out-of-the-untransformed-vector-section must be preserved for next block
		bw = cy = 0;	// These are to manage the 2-blocks-down forwarding of borrows out of the calls to
						// mi64_cvt_double_uint64 which convert signed FFT-mul-linear-combo data back to uint64[].
		cy_re = cy_im = 0ull;
		for(kk = 0; kk < vec_len; kk += sec_len) {
			// Only really need to clear hi half since lo half overwritten by mi64_cvt_uint64_double call:
			memset(c, 0, (pad_len<<3));
			// Process length = sec_len section starting at index [kk]:
			cy1 = cy_fwd & 1ull; cy2 = (cy_fwd>>1) & 1ull;//	printf("kk = %u: cy_re,im = %lld,%lld, cy_fwd = %u\n",kk,cy_re,cy_im,cy_fwd);
			cy_fwd = mi64_cvt_uint64_double(u+kk,v+kk,cy_fwd, sec_len, c);
		//	printf("\tlo64 check1: au = %llu; bv = %llu; cu = %llu; dv = %llu\n",ai[0]*(u[kk]+cy1),bi[0]*(v[kk]+cy2),ci[0]*(u[kk]+cy1),di[0]*(v[kk]+cy2));
			// Forward FFT of U/V-vector section [in c]:
			pairFFT_mul(  c,0,0, fft_len, FALSE, TRUE );	// INIT_ARRAYS = FALSE, FORWARD_FFT_ONLY = TRUE
			// d = copy of 0-padded u/v data stored in fFFTed form in c:
			memcpy(d, c, (pad_len<<3));
			// Do the [a*u, b*v] FFT-mul; the subtract-and-store-result-in-Re-part occurs in the carry step:
			pairFFT_mul(  c,a,0, fft_len, FALSE, 2*TRUE);	// INIT_ARRAYS = FALSE, FORWARD_FFT_ONLY = 2 (indicates both inputs already fFFTed)
			// Do the [c*u, d*v] FFT-mul:
			pairFFT_mul(  d,b,0, fft_len, FALSE, 2*TRUE);	// INIT_ARRAYS = FALSE, FORWARD_FFT_ONLY = 2
			// Interleave the [a*u-b*v, c*u-d*v] outputs; c*u-d*v go into odd slots of C-array:
			for(i = 0; i < fft_len; i+=2) {
				j = i + ( (i >> DAT_BITS) << PAD_BITS );
				c[j] *= sign_mult[sign_ab];	c[j+1] = sign_mult[sign_cd]*d[j];
			}
			// Convert c/d to uint64, yielding a pair of length-(2*sec_len) vectors x,y plus borrow-outs.
			// This conversion should only ever yield borrows (if anything), never carries, i.e. each of
			// the 2-bit fields encoding cy/bw-outs should == 0 or 3, hence acceptable values for cy = 0,3,12,15:
			bw = cy;	// Save pvs pass value of cy (really a borrow), since we need to preserve these across 2 blocks
						// After 1 loop pass, may have a nonzero cy ... thus after 2 passes we may have a nonzero bw.
			cy = mi64_cvt_double_uint64(c,fft_len, x,y);//	printf("\tmi64_cvt_double_uint64: cy_lo,hi = %u,%u\n",cy&3,cy>>2);
			if(cy)
				ASSERT(HERE, cy==3 || cy==12 || cy==15, "Unexpected carryout!");
		//	printf("\tlo64 check2: x0 = %llu; y0 = %llu\n",x[0],y[0]);
			ASSERT(HERE, SGN(x[0],sign_ab) == ai[0]*(u[kk]+cy1) - bi[0]*(v[kk]+cy2) && SGN(y[0],sign_cd) == ci[0]*(u[kk]+cy1) - di[0]*(v[kk]+cy2), "lo64 check fails!");
			// Lower halves of x,y directly overwrite the current u,v-sections:
			mi64_set_eq(u+kk,x,sec_len);	mi64_set_eq(v+kk,y,sec_len);
			// Upper halves of x,y from previous loop execution (stored in w,z) plus borrowouts from
			// mi64_cvt_double_uint64 conversion of those pvs x,y-vectors added to post-multiply sections of u,v:
			if(kk > 0) {
			#if 1	// Dumb no-op approach:
				// Will propagate cyouts into low words of upper halves of x,y when we copy latter into w,z below:
				cy_re = mi64_add(w,u+kk,u+kk, sec_len);
				cy_im = mi64_add(z,v+kk,v+kk, sec_len);
			//	printf("\tAfter previous-block vector-carry add: cy_re,im = %lld,%lld ... ignoring.\n",cy_re,cy_im);
			#else
				if(cy_re < 0) {
					cy_re = mi64_add(w,u+kk,u+kk, sec_len);				// I think we can ignore these if past the new reduced u,v lengths:
					tmp64 = mi64_sub_scalar(u+kk, 1ull, u+kk, sec_len);	if(tmp64 != 0ull) printf("INFO: Nonzero borrow = -%lld out of u--\n",tmp64);
				} else
					cy_re = mi64_add_cyin(w,u+kk,u+kk, sec_len, cy_re);

				if(cy_im < 0) {
					cy_im = mi64_add(z,v+kk,v+kk, sec_len);				// ... but may need handling otherwise. For now, just report:
					tmp64 = mi64_sub_scalar(v+kk, 1ull, v+kk, sec_len);	if(tmp64 != 0ull) printf("INFO: Nonzero borrow = -%lld out of u--\n",tmp64);
				} else
					cy_im = mi64_add_cyin(z,v+kk,v+kk, sec_len, cy_im);
				printf("\tAfter previous-block vector-carry add: cy_re,im = %lld,%lld\n",cy_re,cy_im);
			#endif
			}
			// Upper halves of x,y will get added to low halves of [a*u-bv, c*u-b*v] output for *NEXT* pair of
			// u,v-sections (which have not yet been linearly combined, i.e. must wait until next loop pass),
			// so now that have processed previous pair of vector-addends, save those in w,z for later addition;
			mi64_set_eq(w,x+sec_len,sec_len);
			mi64_set_eq(z,y+sec_len,sec_len);
		#if 1	// Dumb no-op approach:
			cy_re = cy_im = 0ull;
		#else
			if(cy_re < 0) {
				cy_re = mi64_sub_scalar(w,cy_re, w, sec_len);
				if(cy_re != 0ull) printf("\tNonzero bw = %lld out of w--!\n", ABS(cy_re));
			} else {
				cy_re = mi64_add_scalar(w,cy_re, w, sec_len);
				if(cy_re != 0ull) printf("\tNonzero cy = %lld out of w++!\n", ABS(cy_re));
			}
			if(cy_im < 0) {
				cy_im = mi64_sub_scalar(z,cy_im, z, sec_len);
				if(cy_im != 0ull) printf("\tNonzero bw = %lld out of z--!\n", ABS(cy_im));
			} else {
				cy_im = mi64_add_scalar(z,cy_im, z, sec_len);
				if(cy_im != 0ull) printf("\tNonzero cy = %lld out of z++!\n", ABS(cy_im));
			}
			printf("\tAfter w,z incr/decr: cy_re,im = %lld,%lld\n",cy_re,cy_im);
			if(bw) {
				if(sign_ab == 0) cy_re -= ((bw & 3) != 0);	// Only apply borrow if a*u-b*v > 0
				if(sign_cd == 0) cy_im -= ((bw >> 2) != 0);	// Only apply borrow if c*u-d*v > 0
				printf("\tsign_mask*bw = %d*%d,%d*%d ==> cy_re,im = %lld,%lld\n",(sign_ab == 0),-((bw & 3) != 0),(sign_cd == 0),-((bw >> 2) != 0),cy_re,cy_im);
			}
		#endif
		}	// endfor(kk = 0; kk < vec_len; kk += sec_len)
		/**************** FFT-mul a*u - b*v: ****************/
		if(fft_gcd_debug) {
			fp = fopen("c.txt",file_access_mode[vec_len > 16000]);	ASSERT(HERE, fp != 0x0, "Null file pointer!");
			fprintf(fp,"FFT-mul [a*u - b*v, c*u - d*v] for vec_len = %u:\n",vec_len);
			for(i = 0; i < vec_len; i++) {
				fprintf(fp,"%5u\t%20llu\t%20llu\n",i,u[i],v[i]);
			}
			fclose(fp); fp = 0x0;
		//	exit(0);
		}
		/*****************************************************/
		ASSERT(HERE, cy_fwd == 0, "Overall u,v-vectors must be > 0!");
		// For the final block, we *expect* both carries (and all elts of the w,z-vectors) to == -1;
		// the real sign of interest in this case is encoded in the high-order elements of u and v:
		kk -= sec_len;
		cy_re = u[kk-1];	ASSERT(HERE, ABS(cy_re) < 2, "Unexpected Re-carryout!");
		cy_im = v[kk-1];	ASSERT(HERE, ABS(cy_im) < 2, "Unexpected Im-carryout!");
	//	printf("\tOutput signs: cy_re,im = %lld,%lld\n",cy_re,cy_im);
		// If either (a*u - b*v) or (c*u - d*v) result < 0, negate it - to account for multiword uint64[] carry
		// chain here do as in matrix_vector_product_sub - logical negation and reinjection of -cy = +1 at bottom:
		if(cy_re == -1) {
			ASSERT(HERE, 0, "Post-mul negation should not be needed!");
			for(j = 0; j < kk; j++)
				u[j] = ~u[j];	/* Ones' complement */
			// Reinject -cy at bottom, exit as soon as ripple carry peters out:
			for(j = 0; j < kk; j++)
				if(++u[j] != 0) break;
		}
		// Im-part carryout in low 2 bits of cy:
		if(cy_im == -1) {
			ASSERT(HERE, 0, "Post-mul negation should not be needed!");
			for(j = 0; j < kk; j++)
				v[j] = ~v[j];	/* Ones' complement */
			// Reinject -cy at bottom, exit as soon as ripple carry peters out:
			for(j = 0; j < kk; j++)
				if(++v[j] != 0) break;
		}
		// Now compute updated vector lengths:
		i = mi64_getlen(u, kk);		j = mi64_getlen(v, kk);		vec_len = MAX(i,j);
	//	printf("\tFinal words: u[%u],v[%u] = %llu,%llu\n",i,j,u[i-1],v[j-1]);
	printf("New vec_len = %u",vec_len);
//	printf("\n");
	}	// endwhile(vec_len >= 2*sec_len)

	j = mi64_gcd(u, v, vec_len, eGCD = FALSE, 0x0,0x0,0x0,0x0, 0,0x0,0x0,0x0,0x0,0);
	if(EGCD) {
		ASSERT(HERE, 0, "EGCD not yet supported by FFT-gcd!");
	//**** Do final FFTmul of cleanup-mi64_gcd-produced AB multipliers with ones accumulated in FFT-gcd passes ****
	}
	return j;
}

/*
Lehmer-type Euclidean GCD(U, V), where U and V are two multiword integers stored
in base-2^64 positive-digit form; U and V are assumed to both have at least (ndim)
allocated 64-bit words, with U[0],V[0] being the least-significant words.

RETURNS: The GCD in base-2^64 form in both U and V, and the length of the GCD
		(in terms of # of significant 64-bit words) in the function value.

		If the boolean EGCD is set on entry, an extended GCD is performed.
		For simplicity (e.g. so we can use the mi64 functions to manipulate) use uint64 arrays to
		store magnitudes of EGCD multiplier-vectors, separately track the corresponding sign bits.
		Upon completion, pointers to the (absolute values of the) A and B-array multipliers are
		returned in the arglist arrays Ap, Bp, with corresponding signs indicated
		via the low 2 bits of the returned sign_AB variable.
		These array multipliers are defined such that A*U + B*V = GCD(U, V).

		If EGCD is set and the boolean HALF is also set on entry, an extended GCD
		is performed but 'only halfway', in the sense that we return as soon as
		the four multipliers A,B,C,D are all >= len_targ words long ('Half' refers to
		my original thought that this should be >= half the input length; somewhat of
		a misnomer for the final version, in which the user specifies the targeted
		length reduction of the partial-GCD).
			In this case, intended to generate multiword intermediate multipliers for
		use by an FFT-mul-accelerated reduction-via-linear-combination scheme (which
		in a further recursive framework yields a true subquadratic GCD algorithm), we
		also return the C and D-array multipliers in the arglist arrays Cp, Dp, with
		corresponding signs in the low 2 bits of sign_CD. If EGCD = HALF = true, the
		return value is the length of the largest ABCD-multiplier vectors.
		If EGCD = 0, HALF is ignored.

EFFECTS: If both inputs are nonzero, the ORIGINAL INPUT VECTOR CONTENTS ARE DESTROYED
		during the computation, and a copy of the GCD returned in each vector.

		If one or both inputs = 0, neither input is altered, and the function
		return value is = 0 by convention to signal that this is the case.
		In order to match mathematical convention, the user should treat the
		GCD as 0 only if both U and V = 0, and equal to max(U, V)
		(i.e. the nonzero vector among U and V) if only one of the pair = 0.
*/
uint32	mi64_gcd(
	uint64 u[], uint64 v[], uint32 const ndim,
	const uint32 EGCD, uint64 Ap[], uint64 Bp[], uint32 *len_AB, uint32 *sign_AB,
	const uint32 HALF, uint64 Cp[], uint64 Dp[], uint32 *len_CD, uint32 *sign_CD, const uint32 len_targ)
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
	uint64 abmul[2],cdmul[2],max_abcd,t1,t2,t3,t4,lo,hi, uv0[2];
	const double rbase = pow(2.0, 64), log10_base = 19.26591972249479649367928926;	/* 2^64, log10(2^64) */

	double num,den,frac;
	uint32 pass=0, lenU, lenV;	/* Note that lenFoo stores the NUMBER OF SIGNIFICANT ELEMENTS of vector Foo
								(i.e. leading element is in Foo[lenFoo-1]), not the # of allocated elements! */
	uint32 tmp32;
	static uint32 ndim_save = 0;/* Stores the current allocation for the multiplier arrays
								A,B,C,D used for extended GCD. Only used if EGCD = TRUE. */
	uint32 len_abcd;
	static uint64 *A = 0x0,*B = 0x0,*C = 0x0,*D = 0x0;	/* Multiplier arrays for extended GCD.
								For simplicity (e.g. so we can use the mi64 functions to manipulate) use uint64s
								to store the magnitudes, and separately track the corresponding sign bits... */
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
#if GCD_DEBUG>=1
	if(gcd_debug) {
		printf("Inputs to GCD are:\n");
		printf("                         U                   V  \n");

		// Print format uses the PARI/GP arithmetic operator precedence rules, under the
		// assumption that we may be checking the code outputs using that program:
		printf("u=v=0;\n");
		for(i=0; i<ndim; i++)
			printf("i = %u; u+=%20llu<<(i<<6); v+=%20llu<<(i<<6);\n",i,u[i],v[i]);	// For every i++, shift count += 64
		printf("\n");
	}
#endif

	/* Test for EGCD: */
	if(EGCD) {
		ASSERT(HERE, Ap != 0ull && Bp != 0ull && len_AB != 0ull && sign_AB != 0ull, "One or more of the EGCD-required return-pointers null!");
	  if(HALF)
		ASSERT(HERE, Cp != 0ull && Dp != 0ull && len_CD != 0ull && sign_CD != 0ull, "One or more of the HALF-required return-pointers null!");

		if(ndim_save < ndim) {
			// If mul-array allocation not yet performed or insufficient, allocate to hold (ndim) uint64s:
			free((void *)A); A = (uint64 *)calloc(ndim, sizeof(uint64));
			free((void *)B); B = (uint64 *)calloc(ndim, sizeof(uint64));
			free((void *)C); C = (uint64 *)calloc(ndim, sizeof(uint64));
			free((void *)D); D = (uint64 *)calloc(ndim, sizeof(uint64));
			ndim_save = ndim;
		//	fprintf(stderr,"mi64_gcd alloc with ndim = %u\n",ndim);
		}
		memset(A,0,ndim<<3);	memset(B,0,ndim<<3);	memset(C,0,ndim<<3);	memset(D,0,ndim<<3);

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
		the initial signs are A/B = +/- and C/D = -/+, we store sgn(A) in bit 0, sgn(B) in bit 1,
		thus A/C = +- maps to bit pattern 10 in the low 2 bits of sign_abcd: */
		sign_abcd = 0x2;
		A[0] = D[0] = 1;
		B[0] = C[0] = 0;
	}

	/* Init the pointers to the U/V arrays: */
	uv_ptr[0] = u;	uv_ptr[1] = v;

	/* Init the pointers to the A/B/C/D arrays: */
	ac_ptr[0] = A;	bd_ptr[0] = B;
	ac_ptr[1] = C;	bd_ptr[1] = D;
	uv0[0] = u[0];	uv0[1] = v[0];	// Save low limb of each input for EGCD sign determination and sanity check

FIND_LARGER:

	/* k1 and k2 index into the uv_ptr array - uv_ptr[k1] points to larger of u and v: */
	istart = ndim-1;
	lenU = 0;	/* This allows us to determine whether lenU has already been set on a previous pass in the ensuing loop. */
	for(i=istart; i>=0; i--) {
		if((u[i] != 0 || v[i] != 0) && lenU == 0)
			lenU = i+1;
		/* Even if we've already found the leading nonzero element,
		in order to properly set k1 and k2 continue until find a differing element:
		*/
		if(u[i] == v[i])
			continue;
		else {
			k1 = (u[i] < v[i]);
			k2 = k1 ^ 0x00000001;
			break;
		}
	}
	/*	If U = V, both inputs identical:
	*/
	if(i == -1) {
		printf("INFO: Identical inputs to MI64_GCD.\n");
		lenV = lenU;	k1 = 0; k2 = 1;
		goto GCD_DONE;
	}

	/* Otherwise, parse the smaller array until find leading nonzero element: */
	x = uv_ptr[k2];
	for(i=lenU-1; i>=0; i--) {
		if(x[i] != 0) break;
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
	if(EGCD) {
		if     (tz0 > tz1) {
			/* Check left-shift return value to see if len_abcd needs to be updated: */
// 2005: ******************** start here - need to pad length by abs(tz0-tz1)/64 + 1 prior to shifting *********************************
			mi64_shl(ac_ptr[k2], ac_ptr[k2], tz0-tz1, len_abcd);
			mi64_shl(bd_ptr[k2], bd_ptr[k2], tz0-tz1, len_abcd);	len_abcd = MAX(mi64_getlen(ac_ptr[k2], len_abcd) , mi64_getlen(bd_ptr[k2], len_abcd));
		} else if(tz0 < tz1) {
			/* len_abcd will be automatically updated during the 2 shifts
			so the final value corresponds to the larger of the 2 result vectors: */
			mi64_shl(ac_ptr[k1], ac_ptr[k1], tz1-tz0, len_abcd);
			mi64_shl(bd_ptr[k1], bd_ptr[k1], tz1-tz0, len_abcd);	len_abcd = MAX(mi64_getlen(ac_ptr[k1], len_abcd) , mi64_getlen(bd_ptr[k1], len_abcd));
		}
	}

	/* If this is our second and final pass this way, make sure both operands are odd, i.e. tz0 = tz1 = 0: */
	if(min_trailing_zeros != 0)
		ASSERT(HERE, MAX(tz0,tz1) == 0,"gcd_lehmer.c: MAX(tz0,tz1) == 0");
	else if(MIN(tz0,tz1) > 0) {
		ASSERT(HERE, min_trailing_zeros == 0,"gcd_lehmer.c: min_trailing_zeros == 0");	/* Should only have to do this once. */
		min_trailing_zeros = MIN(tz0,tz1);
		goto FIND_LARGER;
	}

  #if GCD_DEBUG>=2
	if(gcd_debug && (tz0 || tz1)) {
		printf("Right-justified Inputs to GCD are:\n");
		printf("                         U                   V  \n");
		printf("u=v=0;\n");
		for(i=0; i<lenU; i++)
			printf("i = %u; u+=%20llu<<(i<<6); v+=%20llu<<(i<<6);\n",i,u[i],v[i]);	// For every i++, shift count += 64
		printf("\n");
	}
  #endif

#endif

/*	Main loop is here.	*/

#if GCD_DEBUG>=2
if(gcd_debug) printf("PASS %d: lenU = %u, lenV = %u, k1,k2 = %d,%d\n",pass,lenU,lenV,k1,k2);
#endif

for(;;)	/* MAIN */
{
	pass++;

#if 1
	lz0 = mi64_leadz(uv_ptr[k1], lenU);
	lz1 = mi64_leadz(uv_ptr[k2], lenV) + ((lenU - lenV)<<6);
	ASSERT(HERE, lz0 <= lz1,"gcd_lehmer.c: lz0 > lz1!");
	#if GCD_DEBUG>=2
		if(gcd_debug) printf("lz_diff = %u bits\n",lz1-lz0);
	#endif
	if(lz1 - lz0 > 64)
	{
		sprintf(string0,"U (lz0 = %u bits) more than 64 bits larger than V (lz1 = %u bits) in MI64_GCD!\n",lz0,lz1);
		WARN(HERE, string0, "", TRUE);
	}
#endif
/**** 8/23/2005: disable for now, since this really complicates things when doing eGCD. ****/
	/* In order to guarantee a nontrivial 2x2 scalar multiplier matrix,
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
		#if GCD_DEBUG>=1
			if(gcd_debug) printf("pass %u , A: left-shifting V-array by %d words\n",pass,len_diff);
		#endif
			mi64_shl(uv_ptr[k2], uv_ptr[k2],  len_diff<<6, lenU);
			lenV = mi64_getlen(uv_ptr[k2], lenU);
			ASSERT(HERE, lenV == lenU,"gcd_lehmer.c: lenV == lenU");
		}
		/*	otherwise left-shift V by (len_diff - 1) places.	*/
		else if(len_diff > 1)
		{
		#if GCD_DEBUG>=1
			if(gcd_debug) printf("pass %u , B: left-shifting array %u by %d words\n",pass,k2,len_diff-1);
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
		j = lz1 - lz0;	ASSERT(HERE, j >= 0,"ERROR: lz0 > lz1 in MI64_GCD!");
		if(j >= (max_bits>>1))
		{
			// If larger vector has any trailing zeros resulting from a previous such shift, right-justify it first:
			i = trailz64(uv_ptr[k1][0]);
			if(i) {
				ASSERT(HERE, i < j, "Excessive trailing zeros in MI64_GCD!");
				mi64_shrl(uv_ptr[k1], uv_ptr[k1], i, lenU);	lenU = mi64_getlen(uv_ptr[k1], lenU);
			}
		#if GCD_DEBUG>=1
			if(gcd_debug) {
				printf("pass %u, C: left-shifting vector %u by %u-%u-1 = %u bits;\n",pass,k2,j,i,j-i-1);
			//	printf("Input  = %s\n",&str_10k[__convert_mi64_base10_char(str_10k, 10<<10, uv_ptr[k2], lenV, 0)]);
			}
		#endif
			// Must use lenU as shift-operand length, since result may have a carry into next-higher word, which is discarded otherwise.
			mi64_shl(uv_ptr[k2], uv_ptr[k2], j-i-1, lenU);	lenV = mi64_getlen(uv_ptr[k2], lenU);
		}
	}

  #if GCD_DEBUG>=2
	if(gcd_debug && pass>=0) {
		printf("Left-justified Inputs to GCD are:\n");
		printf("                         U                   V  \n");
		printf("u=v=0;\n");
		for(i=0; i<lenU; i++) {
			printf("i = %u; u+=%20llu<<(i<<6); v+=%20llu<<(i<<6);\n",i,u[i],v[i]);	// For every i++, shift count += 64
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
if(gcd_debug) {
	printf("Indices into U/V = %u, %u\n",k1,k2);
	printf("\n");
	printf("NSHIFT = %u:\n", nshift);
	printf("\n");
	printf("leadU[1,0] = %20llu   %20llu\n",lead_u[1],lead_u[0]);
	printf("leadV[1,0] = %20llu   %20llu\n",lead_v[1],lead_v[0]);
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
		if(den == 0) {
		#if GCD_DEBUG>=2
			if(gcd_debug) printf("EXIT oo\n");
		#endif
			break;	/* EXIT CALC_MULTS */
		}

		frac = num/den;
		iop = (uint64)frac;	// Use FDIV to calculate FP approximate to ratio, then integer truncate upon store.
	#if GCD_DEBUG>=2
		if(gcd_debug) printf("frac = %20.15e, iop = %20llu\n",frac,iop);
	#endif

		// If FP result is just a smidge < 1, force = 1 to prevent premature loop exit for nearly-indeitcal 128-bit num,den values:
		if(iop == 0 && frac > 0.999999999999999) {
			iop = 1;
			#if GCD_DEBUG>=2
				if(gcd_debug) {
				//	printf("frac = %20.15e, iop = %20llu\n",frac,iop);
					printf("n1 = %20llu, n0 = %20llu\n",range[i][1],range[i][0]);
					printf("d1 = %20llu, d0 = %20llu\n",t2,t1);
				}
			#endif
		}

		/* If ratio = 0 or is out of range of uint64, break: */
		if(iop == 0) {
		#if GCD_DEBUG>=2
			if(gcd_debug) printf("EXIT 0\n");
		#endif
			break;	/* EXIT CALC_MULTS */
		} else if(frac >= TWO63FLOAT) {	// Q: Should we require ratio to be in range of double-float mantissa instead?
		#if GCD_DEBUG>=2
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
		if((t3 < lo) || hi!=0) {
		#if GCD_DEBUG>=2
			if(gcd_debug) printf("EXIT A\n");
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
		if(hi!=0 || (t3>>max_bits) != 0) {
		#if GCD_DEBUG>=2
			if(gcd_debug) printf("EXIT B\n");
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
		if(hi!=0 || (t4>>max_bits) != 0) {
		#if GCD_DEBUG>=2
			if(gcd_debug) printf("EXIT C\n");
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
if(gcd_debug) {
	printf("\n");
	printf("ABmul=%20llu   %20llu\n",coeff[0][0],coeff[0][1]);
	printf("CDmul=%20llu   %20llu\n",coeff[1][0],coeff[1][1]);
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

#if GCD_DEBUG >= 2
if(gcd_debug) {
	printf("a = %20llu; b = %20llu;\n", abmul[ 0], abmul[ 1]);
	printf("c = %20llu; d = %20llu;\n", cdmul[ 0], cdmul[ 1]);
}
#endif

	if(max_abcd <= 1) {
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
if(gcd_debug && pass==0) {
	printf("Result of matrix_vector_product:\n");
	printf("                         U                   V  \n");
	printf("u=v=0;\n");
	for(i=0; i<lenU; i++)
		printf("i = %u; u+=%20llu<<(i<<6); v+=%20llu<<(i<<6);\n",i,u[i],v[i]);	// For every i++, shift count += 64
	printf("\n");
}
#endif

	if(EGCD) {
	#if 0/*GCD_DEBUG>=2*/
		/* Compute the leading 128 bits of the
		A' = (a*A + b*C), B' = (a*B + b*D), C' = (c*A + d*C), D' = (c*B + d*D)
		products and the leading 256 bits of the resulting A'*u,B'*v,C'*u,D'*v terms.
		Then repeat, but for the subtractive
		A" = (a*A - b*C), B" = (a*B - b*D), C" = (c*A - d*C), D" = (c*B - d*D) terms.
		This is to see if we want to + or - in updating the eGCD multipliers.
		*/
		if(gcd_debug) {
			/* Leading 128 bits of A/B: */
			if(mi64_cmpugt(A, B, len_abcd) {
				tmp32 = mi64_getlen(A, len_abcd);
				i = leadz64(A, tmp32);
			} else {
				tmp32 = mi64_getlen(B, len_abcd);
				i = leadz64(B, tmp32);
			}
			mi64_extract_lead128(A, tmp32, i, a128);
			mi64_extract_lead128(B, tmp32, i, b128);

			/* Leading 128 bits of C/D: */
			if(mi64_cmpugt(C, D, len_abcd) {
				tmp32 = mi64_getlen(C, len_abcd);
				i = leadz64(C, tmp32);
			} else {
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
		to update eGCD multipliers: */

/*****HOW TO PROPERLY DEAL WITH CASE WHERE A*u-B*v or C*u-D*v = 0?******/

		// No more need for this - found it easier to simply fix the needed signs at very end (cf. 'Nov 2015' note below):
		ASSERT(HERE, egcd_sign_swap == 0, "egcd_sign_swap unexpectedly set!");
		if(egcd_sign_swap) {
			ASSERT(HERE, 0, "Bad code!");
		#if GCD_DEBUG>=2
			if(gcd_debug) {
				printf("Pass = %u\n", pass);
				printf("A/B eGCD multiplier vectors:\n");
				if(sign_abcd & 0x1)
					printf("                         -A                  +B  \n");
				else
					printf("                         +A                  -B  \n");

				printf("a=b=0;\n");
				for(i=0; i<len_abcd; i++)
					printf("i = %u; a+=%20llu<<(i<<6); b+=%20llu<<(i<<6);\n",i,A[i],B[i]);
				printf("\n");

				printf("C/D eGCD multiplier vectors:\n");
				if(sign_abcd & 0x2)
					printf("                         -C                  +D  \n");
				else
					printf("                         +C                  -D  \n");

				printf("c=d=0;\n");
				for(i=0; i<len_abcd; i++)
					printf("i = %u; c+=%20llu<<(i<<6); d+=%20llu<<(i<<6);\n",i,C[i],D[i]);
				printf("\n");
			}
		#endif

			/* Don't combine the 2 calls into an if(call1() || call2()),
			since we need to be sure both function calls always occur: */
			tmp32  = matrix_vector_product_sub(abmul, cdmul, ac_ptr, len_abcd);
			tmp32 |= matrix_vector_product_sub(abmul, cdmul, bd_ptr, len_abcd);
			len_abcd += (tmp32&1);
		} else {
			tmp32  = matrix_vector_product_add(abmul, cdmul, ac_ptr, len_abcd);
			tmp32 |= matrix_vector_product_add(abmul, cdmul, bd_ptr, len_abcd);
			len_abcd += (tmp32&1);
		}

		// Early-abort in case using eGCD to accumulate "chunks of yardage" for subsequent FFT-mul reduction:
		if(HALF && len_abcd >= len_targ) goto HALF_EXIT;

		/* Lookahead scheme? */
		if(j1 ^ j2)
			egcd_sign_swap = FALSE;
		else
			egcd_sign_swap = TRUE;

	#if GCD_DEBUG>=2
		if(gcd_debug) {
			printf("A/B eGCD multiplier vectors, sign_abcd = %u:\n",sign_abcd);
			printf("a=b=0;\n");
			for(i=0; i<len_abcd; i++)
				printf("i = %u; a+=%20llu<<(i<<6); b+=%20llu<<(i<<6);\n",i,A[i],B[i]);
			printf("\n");

			printf("C/D eGCD multiplier vectors:\n");
			printf("c=d=0;\n");
			for(i=0; i<len_abcd; i++)
				printf("i = %u; c+=%20llu<<(i<<6); d+=%20llu<<(i<<6);\n",i,C[i],D[i]);
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
	if(tz0) {
		mi64_shrl(u, u, tz0, lenU);
shift_cd += tz0;
		if(EGCD) {
			/* Left-shift the complement-vector eGCD multipliers (C/D) by tz0 places: */
			mi64_shl(ac_ptr[1], ac_ptr[1], tz0, len_abcd);
			mi64_shl(bd_ptr[1], bd_ptr[1], tz0, len_abcd);
			len_abcd = MAX(mi64_getlen(ac_ptr[1], len_abcd) , mi64_getlen(bd_ptr[1], len_abcd));
		}
	}

	tz1 = mi64_trailz(v, lenU);
	if(tz1) {
		mi64_shrl(v, v, tz1, lenU);
shift_ab += tz1;
		if(EGCD) {
			/* Left-shift the complement-vector eGCD multipliers by (A/B) tz1 places: */
			mi64_shl(ac_ptr[0], ac_ptr[0], tz1, len_abcd);
			mi64_shl(bd_ptr[0], bd_ptr[0], tz1, len_abcd);
			len_abcd = MAX(mi64_getlen(ac_ptr[0], len_abcd) , mi64_getlen(bd_ptr[0], len_abcd));
		}
	}
#endif

	/*	set k1 to row index of larger element: */
	istart = lenU-1;
	lenU = 0;	/* This allows us to determine whether lenU has already been set on a previous pass in the ensuing loop. */
	for(i=istart; i>=0; i--) {
		if((u[i] != 0 || v[i] != 0) && lenU == 0)
			lenU = i+1;
		/* Even if we've already found the leading nonzero element,
		in order to properly set k1 and k2 continue until find a differing element:
		*/
		if(u[i] == v[i])	//<*** On pass 2, get upper 2 slots = 0 as expected, but next 350 or so all have u[i] = v[i] = 14757395258967641292 ***
			continue;		//That leads to u/v having leading 350 64-bit words identical, i.e. on next cal-mults pass expect a linear combo = u-v
		else {				//in on of the slots, other should be u or v unchanged, i.e. scalar-mults matrix = 1  -1 / 1   0 .
			k1 = (u[i] < v[i]);
			k2 = k1 ^ 0x00000001;
			break;
		}
	}
	/*	If U = V, we're done, with a copy of a nontrivial GCD in each of the 2 vectors. */
	if(i == -1) goto GCD_DONE;

	/* Otherwise, parse the smaller array until find leading nonzero element: */
	x = uv_ptr[k2];
	for(i=lenU-1; i>=0; i--) {
		if(x[i] != 0) break;
	}
	lenV = i + 1;
	/* If lenV = 0, we're done. */
	if(lenV == 0) goto GCD_DONE;

#if GCD_DEBUG>=2
if(gcd_debug || (pass&1023)==0) printf("PASS %d: lenU = %u, lenV = %u, k1,k2 = %d,%d\n",pass,lenU,lenV,k1,k2);
#endif
}	/* enddo MAIN	*/

GCD_DONE:

	/*	Print GCD: */
#if GCD_DEBUG>=0
	printf("GCD finished in %u passes.\n", pass);
#endif
	ASSERT(HERE, (lenU+lenV != 0), "mi_gcd: final U & V-vectors both zero-length!");

	/* Put copies of the GCD into both uv_ptr[k1] and uv_ptr[k2]: */
	for(i = 0; i < lenU; i++) {
		uv_ptr[k2][i] = uv_ptr[k1][i];
	}

	/* Nov 2015:
	Just compute the EGCD multiplies A,B,C,D sans signs through course of GCD, then at end simply
	do a 64-bit low-half multiply of 0-words of A*u and B*v, whichever of these 2 64-bit products
	is 1 greater than the other gets the + sign:
	*/
	if(EGCD) {
		uv0[0] *= A[0];	uv0[1] *= B[0];
		sign_abcd = (uv0[1] > uv0[0]);	// Result = index of whichever of A/B gets the + sign in GCD = +-A*u-+B*v :
		ASSERT(HERE, uv0[sign_abcd] - uv0[sign_abcd ^ 1] == 1ull, "EGCD multiplier sanity check failed!");
	}

HALF_EXIT:

	if(EGCD) {
		mi64_set_eq(Ap,A,len_abcd); mi64_set_eq(Bp,B,len_abcd); *len_AB = len_abcd; *sign_AB = sign_abcd & 0x1;
	  if(HALF) {
		mi64_set_eq(Cp,C,len_abcd); mi64_set_eq(Dp,D,len_abcd); *len_CD = len_abcd; *sign_CD = sign_abcd >> 1;
		ASSERT(HERE, len_abcd == len_targ, "HALF_EXIT: len_abcd != len_targ");
		// The sign_abcd-setting in the code above doesn't work in the half-gcd case, so instead
		// use signs of LSWs of inputs and outputs to deduce signs of a*u-b*v, c*u-d*v:
		lo = Ap[0]*uv0[0] - Bp[0]*uv0[1]; if(lo == u[0]) { *sign_AB = 0; } else if(lo == -u[0]) { *sign_AB = 1; } else { ASSERT(HERE, 0, "foo!"); }
		lo = Cp[0]*uv0[0] - Dp[0]*uv0[1]; if(lo == v[0]) { *sign_CD = 0; } else if(lo == -v[0]) { *sign_CD = 1; } else { ASSERT(HERE, 0, "bar!"); }
		// If HALF = true, return value is the length of the largest of the ABCD-multiplier vectors.
		return len_abcd;
	  }
	}

#if GCD_DEBUG>=1
	if(EGCD && gcd_debug) {
		ASSERT(HERE, min_trailing_zeros == 0, "min_trailing_zeros != 0");

		printf("A/B eGCD multiplier vectors:\n");
		if(sign_abcd & 0x1)
			printf("                         -A                  +B  \n");
		else
			printf("                         +A                  -B  \n");

		printf("a=b=0;\n");
		for(i=0; i<len_abcd; i++)
			printf("i = %u; a+=%20llu<<(i<<6); b+=%20llu<<(i<<6);\n",i,A[i],B[i]);
		printf("\n");

		printf("C/D eGCD multiplier vectors:\n");
		if(sign_abcd & 0x2)
			printf("                         -C                  +D  \n");
		else
			printf("                         +C                  -D  \n");

		printf("c=d=0;\n");
		for(i=0; i<len_abcd; i++)
			printf("i = %u; c+=%20llu<<(i<<6); d+=%20llu<<(i<<6);\n",i,C[i],D[i]);
		printf("\n");

		printf("GCD vectors are:\n");
		printf("                         U                   V  \n");
		printf("u=v=0;\n");
		for(i=0; i<lenU; i++)
			printf("i = %u; u+=%20llu<<(i<<6); v+=%20llu<<(i<<6);\n",i,u[i],v[i]);	// For every i++, shift count += 64
		printf("\n");
	}
#endif

	// Restore any common power of 2: *** Still need to generalize this to EGCD case ***
	if(!EGCD && !HALF && min_trailing_zeros > 0) {
		lenU += (min_trailing_zeros>>6)+1;
		mi64_shl(uv_ptr[k1], uv_ptr[k1], min_trailing_zeros, lenU);
		lenU = mi64_getlen(uv_ptr[k1], lenU);	/* 12/16/08: getlen() had lenV here. (?) */
	}

	// Skip GCD-printing if sending leading chunks of large vectors to eGCD in order
	// to generate multiword abcd-multipliers for FFT-mul-speeded length reduction:
	if(!HALF) {
		/* Currently printing functions can only handle up to STR_MAX_LEN digits: */
		if(lenU*log10_base >= STR_MAX_LEN) {
			printf("GCD exceeds printable string length - printing in raw base-2^64 form, least-significant first:\n");
			printf("GCD = 2^%u * [\n", min_trailing_zeros);
			for(i = 0; i < lenU; i++)
				printf("i = %u:  %20llu\n",i,uv_ptr[k1][i]);
			printf("]\n");
		} else {
			printf("GCD = %s\n", &string0[convert_mi64_base10_char(string0, uv_ptr[k1], lenU, 0)]);
		}
	}
	// Return value = number of 64-bit words needed to hold the GCD:
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
		nshift = leadz64(u[len-1]);	/* If U > V, use it to set shift count... */
	else
		nshift = leadz64(v[len-1]);	/* ...otherwise use V to set shift count. */

	/*	Figure out which order to do the multiplies in (i.e. a*u - b*v or b*v - a*u)
	to ensure (or at least maximize the chance, as O(len)-costly postprocessing is
	needed otherwise) that result >= 0:
	*/

	/*	If leading 192 bits of a*u < b*v, swap a,b and the ptrs to u,v: */
	mi64_extract_lead128(u,len, nshift, lead_u);	//	ASSERT(HERE, (lead_u[0]+lead_u[1]) != 0, "(lead_u[0]+lead_u[1]) == 0!");
	mi64_extract_lead128(v,len, nshift, lead_v);	//	ASSERT(HERE, (lead_v[0]+lead_v[1]) != 0, "(lead_v[0]+lead_v[1]) == 0!");

	j1 = 0;
	if(CMP_LT_PROD192(abmul[0], lead_u[0], lead_u[1], abmul[1], lead_v[0], lead_v[1]))
	{
	#if GCD_DEBUG>=2
		if(gcd_debug) printf("(a*u) < (b*v) : swapping these\n");
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
		if(gcd_debug) printf("(c*u) < (d*v) : swapping these\n");
	#endif
		j3 ^= 0x00000001;
		retval ^= 0x00000004;
	}
	j4 = j3^0x00000001;

	/*
	Proceed with the pairwise-linear-combination step, writing results back into u/v:
	*/
	/* Check magnitudes of a,b,c,d < 2^63: */
	a = abmul[j1];	b = abmul[j2];	c = cdmul[j3];	d = cdmul[j4];
	ASSERT(HERE, a >= 0 && b >= 0 && c >= 0 && d >= 0,"matrix_vector_product_sub: a >= 0 && b >= 0 && c >= 0 && d >= 0");

	A = uv_ptr[j1];	B = uv_ptr[j2];	C = uv_ptr[j3];	D = uv_ptr[j4];

#ifndef YES_ASM	// Std-C version of the critical loop:

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
		// These 2 are possibly-swapped copies to same data as val1,2, must init before overwriting u,v!
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

#else	// YES_ASM = true: x86_64 inline-asm version of the critical loop:
		// We get rid of most refs to var-signedness here, since at the machine-instruction level
		// that comes about via choice of instructions (e.g. arithmetic vs logical right-shift),
		// and only matters for a very few instructions.

  #if 1
	// Oct 2015: This inline-ASM version of the loop runs 25% faster than the C version on my Core2;
	// The main gain came not from the loop-body (which I ASM-ized) first, but from the finishing
	// touch of full-ASM of the loop control.
	__asm__ volatile (\
		/* 2 of our 4 vector-int addresses can be load-once, as these are not overwritten in the loop: */\
		"movq	%[__A],%%r12	\n\t	movq	%[__B],%%r13	\n\t"\

	"movslq	%[__len], %%rcx	\n\t"/* ASM loop structured as for(j = len; j != 0; --j){...} */\
	"xorq	%%rsi, %%rsi	\n\t"/* Index into the 3 arrays (really a uint64-array pointer offset), init = 0 */\
	"movq	%%rcx,%%rdi	\n\t"/* Copy of array-length in rdi */\
	"1:					\n\t"\

		/* Remaining 2 of our 4 vector-int addresses must be reloaded each loop pass: */\
		"movq	%[__C],%%r14	\n\t	movq	%[__D],%%r15	\n\t"\
		"movq	%[__cy1],%%r8 	\n\t	movq	%[__cy2],%%r9 	\n\t"/* Load cy1,2 */\
		"movq	(%%r12,%%rsi),%%rax 	\n\t"/* val1 = A[i] */\
		"mulq	%[__a]					\n\t"/* MUL_LOHI64(a, val1); lo1:hi1 in rax:rdx; CF (unused) set according to (hi == 0) */\
		"movq	%%rax,%%rbx				\n\t"/* lo1 */\
		"movq	%%rdx,%%rdi				\n\t"/* hi1 */\
		"movq	(%%r13,%%rsi),%%rax 	\n\t"/* val2 = B[i] */\
		"mulq	%[__b]					\n\t"/* MUL_LOHI64(b, val2); lo2:hi2 in rax:rdx */\

		/* Load val3,4 before overwriting u,v[i] below */\
		"movq	(%%r14,%%rsi),%%r14 	\n\t"/* val3 = C[i] */\
		"movq	(%%r15,%%rsi),%%r15 	\n\t"/* val4 = D[i] */\

		"movq	%%r8 ,%%r10				\n\t"/* cylo = cy1, cylo in r10 */\
		"sarq	$63,%%r8 				\n\t"/* cy1 = -(cy1 < 0); Sign extension of cy1 */\
		"subq	%%rbx,%%rax				\n\t"/* lo[rax] = lo2 - lo1 */\
		"adcq	$0 ,%%r8 				\n\t"/* cy1 += Borrow out of (lo2 - lo1) */\
		"subq	%%rax,%%r10				\n\t"/*        cylo - lo */\
		"sbbq	$0 ,%%r8 				\n\t"/* cy1 -= Borrow out of (cylo - lo) */\
		"movq	%[__u],%%r11			\n\t"\
		"movq	%%r10,(%%r11,%%rsi)	\n\t"/* u[i] = cylo - lo */\
		"addq	%%rdi,%%r8 				\n\t"/* cy1 += hi1 */\
		"subq	%%rdx,%%r8 				\n\t"/* cy1 -= hi2 */\
		"movq	%%r8 ,%[__cy1]		\n\t"/* Store cy1 */\

		"movq	%%r14,%%rax 	\n\t"/* val3 = C[i] */\
		"mulq	%[__c]					\n\t"/* MUL_LOHI64(c, val1); lo1:hi1 in rax:rdx */\
		"movq	%%rax,%%rbx				\n\t"/* lo1 */\
		"movq	%%rdx,%%rdi				\n\t"/* hi1 */\
		"movq	%%r15,%%rax 	\n\t"/* val4 = D[i] */\
		"mulq	%[__d]					\n\t"/* MUL_LOHI64(d, val2); lo2:hi2 in rax:rdx */\
		"movq	%%r9 ,%%r10				\n\t"/* cylo = cy2, cylo in r10 */\
		"sarq	$63,%%r9 				\n\t"/* cy2 = -(cy2 < 0); Sign extension of cy2 */\
		"subq	%%rbx,%%rax				\n\t"/* lo[rax] = lo2 - lo1 */\
		"adcq	$0 ,%%r9 				\n\t"/* cy2 += Borrow out of (lo2 - lo1) */\
		"subq	%%rax,%%r10				\n\t"/*        cylo - lo */\
		"sbbq	$0 ,%%r9 				\n\t"/* cy2 -= Borrow out of (cylo - lo) */\
		"movq	%[__v],%%r11			\n\t"\
		"movq	%%r10,(%%r11,%%rsi)	\n\t"/* v[i] = cylo - lo */\
		"addq	%%rdi,%%r9 				\n\t"/* cy2 += hi1 */\
		"subq	%%rdx,%%r9 				\n\t"/* cy2 -= hi2 */\
		"movq	%%r9 ,%[__cy2]		\n\t"/* Store cy2 */\

	"leaq	0x8(%%rsi), %%rsi	\n\t"/* Incr shared mem-offset */\
	"decq	%%rcx \n\t"\
	"jnz 1b 	\n\t"/* loop1 end; continue is via jump-back if rcx != 0 */\
		: /* outputs: none */\
		: /* All inputs from memory here (Trying to use mem/reg "g" gives 'invalid operand for instruction' errors) */\
		  [__u] "m" (u), [__v] "m" (v)	\
		 ,[__a] "m" (a), [__b] "m" (b), [__c] "m" (c), [__d] "m" (d)	\
		 ,[__A] "m" (A), [__B] "m" (B), [__C] "m" (C), [__D] "m" (D)	\
		 ,[__cy1] "m" (cy1), [__cy2] "m" (cy2)	\
		 ,[__len] "m" (len), [__i] "m" (i)	\
		: "cc","memory","rax","rbx","rdx","rsi","rdi","r8","r9","r10","r11","r12","r13","r14","r15"	/* Clobbered registers */\
	);

  #else

	for (i = 0; i < len; i++)
	{
		val1 = A[i];
		val2 = B[i];
		MUL_LOHI64(a, val1, lo1, hi1);
		MUL_LOHI64(b, val2, lo2, hi2);
		// These 2 are possibly-swapped copies to same data as val1,2, must init before overwriting u,v!
		val3 = C[i];
		val4 = D[i];

	/* a*u[i] - b*v[i] */
		cylo = cy1;
		cy1  = -(cy1 < 0);	/* Sign extension of cy1 */
		lo   = lo2 - lo1;
		cy1 += (lo2 < lo1);	// cy1 += Borrow out of (lo2 - lo1)
		u[i] = cylo - lo;
		cy1 -= (cylo < lo);	// cy1 -= Borrow out of (cylo - lo)
		cy1 += hi1;
		cy1 -= hi2;

	/* c*u[i] - d*v[i] */
		MUL_LOHI64(c, val3, lo1, hi1);
		MUL_LOHI64(d, val4, lo2, hi2);
		cylo = cy2;
		cy2  = -(cy2 < 0);	/* Sign extension of cy2 */
		lo   = lo2 - lo1;
		cy2 += (lo2 < lo1);	// Borrow out of (lo2 - lo1)?
		v[i] = cylo - lo;
		cy2 -= (cylo < lo);
		cy2 += hi1;
		cy2 -= hi2;
	} /* for i */

  #endif

#endif	// YES_ASM ?

#if GCD_DEBUG>=2
//	printf("cy1 = %lld, cy2 = %lld, cy1|cy2 = %lld\n", cy1, cy2, cy1|cy2);
#endif

	/*	For GCD main-vector updates we know exit carries are 0 or -1,
	but to handle eGCD output carries (which may be as large as the difference of
	two 63-bit ints) need to generalize this to detect any negative 64-bit int: */
	if(cy1 < 0) {
	#if GCD_DEBUG
		printf("cy1 < 0...complementing.\n");
	#endif
		for (i = 0; i < len; i++)
			u[i] = ~u[i];	/* Ones' complement */
		cy1 = ~cy1;
	    /* Add -cy1, exit as soon as result is nonzero, indicating that the carryout of the incrementing is 0: */
		for (i = 0; i < len; i++) {
			u[i]++;
			if (u[i] != 0) break;
		}
		/* Post-processing negation flips any corresponding sign bit in the return field: */
		retval ^= 0x00000002;
	}

	if(cy2 < 0) {
	#if GCD_DEBUG
		printf("cy2 < 0...complementing.\n");
	#endif
		for (i = 0; i < len; i++)
			v[i] = ~v[i];
		cy2 = ~cy2;
	    /* Add -cy2, exit as soon as result is nonzero, indicating that the carryout of the incrementing is 0: */
		for (i = 0; i < len; i++) {
			v[i]++;
			if (v[i] != 0) break;
		}
		retval ^= 0x00000004;
	}

	/* Does vector length need incrementing? */
	if(cy1 || cy2) {
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
  int test_gcd()
  {
	/* Adjust this value to reflect max. vector length desired for self-tests: */
	#define MAX_ARRAY_DIM	344222	// Set >= 34420 to enable the 2-Mbit test, >= 344222 for the 22-Mbit test.
	uint32 i,j,kk,cy,cy_fwd, p, lenX, lenY, lenU, lenV, lenZ, sec_len,vec_len, fft_len,pad_len;
	int fail = 0;
	uint64 q,qinv,tmp, mod, rem, tmp64 = 0;
	uint64 *u, *v, *w, *x, *y, *z;
	double *a, *b;
	/* Extra stuff needed for eGCD: */
	uint32 eGCD = FALSE;	// This can be toggled here for the entire self-test suite, or set/unset on a case-by-case basis.
	uint32 HALF = FALSE;	// Ditto, though note that this one only has an effect if eGCD = TRUE
	uint64 *Ap = 0x0, *Bp = 0x0, *Cp = 0x0, *Dp = 0x0;
	uint32 len_AB,sign_AB, len_CD,sign_CD;

	const uint32 num_mword_prime = 20;	// This value must match the # of const uint64 y**[] entries below
	/* Define precomputed primes of lengths 1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100... words */
	const uint32 len_mword_prime[] = {1,2,3,4,5,6,7,8,9,10,15,20,30,40,50,60,70,80,90,100,200,300,400,500,600,700,800,900,1000};

	/* Multiword primes of various lengths, least-significant [0]-word leftmost: */
	const uint64 y1 [ 1] = { 16357897499336320021ull };
	const uint64 y2 [ 2] = { 6288613973463385207ull, 17450374062388873929ull };
	const uint64 y3 [ 3] = { 10949365252588722703ull, 1187612066073480231ull, 13552797462131857749ull };
	const uint64 y4 [ 4] = { 2290166132859908083ull, 15141729687261062388ull, 3353845001322014797ull, 10668737080789349223ull };
	const uint64 y5 [ 5] = { 1909815749502324189ull, 3409223353484448796ull, 8041232679190714362ull, 18154516404802398221ull, 4919133996046537621ull };
	const uint64 y6 [ 6] = { 13090340153895430355ull, 13316003867294904363ull, 4086877960619352751ull, 4291335014170841853ull, 11238552196362275798ull, 16072837637891425689ull };
	const uint64 y7 [ 7] = { 6436712882813530525ull, 3090069294702847100ull, 14604617259536040931ull, 3881736123518103891ull, 12033780901417973631ull, 2121658959755889151ull, 1015222311804066800ull };
	const uint64 y8 [ 8] = { 14218778164426756345ull, 3881548058709958902ull, 9435344667623458196ull, 15825266020521010528ull, 17646912765656454557ull, 9792466261391312854ull, 13836951808143056624ull, 16061163011110659553ull };
	const uint64 y9 [ 9] = { 15531849983479353575ull, 13244756845085485741ull, 12235243054937127412ull, 4709695542336918044ull, 11468374613656802212ull, 10406739941509043845ull, 5882424034228837618ull, 1440633758287919928ull, 16966630538648771321ull };
	const uint64 y10[10] = { 6154942771176964881ull, 14657802710986676243ull, 6367391185176237819ull, 7220279878409007123ull, 17437454528603741845ull, 11135248797733030345ull, 9466072315348543570ull, 16100604701688312144ull, 10143506867524673304ull, 174955034380965ull };
	const uint64 y15[15] = { 8819938628747633077ull, 5942790028746468569ull, 4937492223042946033ull, 8730501809086437834ull, 1243723839784251790ull, 13716683457145268153ull, 488754166078206799ull, 1949794776773723448ull, 16629452879624053947ull, 1306630474066381588ull, 712163870246867324ull, 4434513443731937057ull, 9865801115224962539ull, 10472193284737100315ull, 17548140925265117600ull };
	const uint64 y20[20] = { 13716683457145268063ull, 488754166078206799ull, 1949794776773723448ull, 16629452879624053947ull, 1306630474066381588ull, 712163870246867324ull, 4434513443731937057ull, 9865801115224962539ull, 10472193284737100315ull, 17548140925265117600ull, 13850887988965333521ull, 15740887654137791778ull, 12429824020706400576ull, 1975658775442752547ull, 13083483410316005473ull, 2987141948949602865ull, 18052154556276307484ull, 15115700697529442503ull, 7255095742800246455ull, 15953457815171762164ull };
	const uint64 y30[30] = { 15546513674694818401ull, 6094687342574418678ull, 11123946222831718112ull, 12217883352413216347ull, 15550874203307310822ull, 7206187696399018449ull, 2089253583232126963ull, 11405684275551645807ull, 6595740366593535533ull, 513947309979730528ull, 623238905451609664ull, 14580584265174886892ull, 10164089984630032021ull, 6295726803342317635ull, 4169976721976571503ull, 11210444189919629661ull, 12253764729310025870ull, 15982375818961640458ull, 4835638556759726427ull, 4550010953562141766ull, 4708309259113265585ull, 7528936866787101258ull, 2993339143445663357ull, 4947092101451202277ull, 3470494295975615532ull, 7612382765688613360ull, 1218883132193815401ull, 14455287304819392609ull, 6565258255864629018ull, 446327481124760089ull };
	const uint64 y40[40] = { 5095253233679814563ull, 9094873052810277402ull, 5502432594497619071ull, 10253180810779792312ull, 11476072751332469357ull, 15370988773499066191ull, 17581225978208792286ull, 526277277651859297ull, 1885992293663728627ull, 10153305797529801974ull, 5464673869670139137ull, 17401163700214643993ull, 17952651913326118589ull, 17742345220822852562ull, 4303343635453243709ull, 4013620749860433596ull, 6677371295215675137ull, 8193725131766686871ull, 7659816303732807063ull, 2282385124936021778ull, 9113451334539710638ull, 10881171535438876472ull, 5476412011083724153ull, 4864713286685264433ull, 5628110023653081947ull, 6926699408022157704ull, 6403858520342783331ull, 6696289090173746018ull, 13288386934025662659ull, 4362860849045370749ull, 16932024427685935439ull, 10744939144246263002ull, 11156351709119530729ull, 5538822441276190107ull, 593146585947091059ull, 16287977467190924709ull, 7931952495909425541ull, 1323052627185110710ull, 13855398606528214565ull, 4866516855423411741ull };
	const uint64 y50[50] = { 803525447701230777ull, 11995973163883103671ull, 12513807631903294061ull, 13206931885979440508ull, 5781832881550132408ull, 2760760068053793998ull, 7279871441580359859ull, 8868806671540311732ull, 5211038110177683731ull, 12119673326921841832ull, 1820023893508557590ull, 15459246329388416701ull, 12576362623841191959ull, 5774269049365068898ull, 4014785386331757453ull, 16752022808542583240ull, 12646727565723050809ull, 17609658751476230416ull, 4567618588919495722ull, 18229471973955028755ull, 16365526703331694869ull, 6014402465390895934ull, 14379352264227400336ull, 7876171534707826122ull, 1228525137268097609ull, 4019073333129686499ull, 9768820187581486855ull, 15898165664400102772ull, 7343357198210053564ull, 17940883182259865333ull, 2218383524838541976ull, 839057351822990983ull, 15720245341684857544ull, 13837110097306272389ull, 18090746938060205953ull, 17626089812589788026ull, 17795072819199975687ull, 9122504416397444822ull, 7804210406658853893ull, 9085825988092390607ull, 38573679668204362ull, 6489947387271485874ull, 1580408331246418996ull, 843254109280709615ull, 9649575090905754564ull, 2194379072572275786ull, 1300697824291809687ull, 6758258572906841718ull, 6618325607590950226ull, 3084702059426541837ull };
	const uint64 y60[60] = { 7402034356984044041ull, 18246800267094834509ull, 9355646276203941041ull, 17708280664550799098ull, 9884444230844509075ull, 2445586268038335580ull, 15849763400441455983ull, 16576388534672358617ull, 827691637990772710ull, 11785188128139116957ull, 9456967256883473605ull, 8353050390135495047ull, 9898448512339846253ull, 3671369225156337572ull, 14994373522692786547ull, 16859594160103298137ull, 14510522032204748689ull, 2512816302537297738ull, 15923263995753141351ull, 7428472312605469391ull, 14733812811272604664ull, 2503760744579306444ull, 1243518692364200645ull, 17335828848204797042ull, 15001078003915713534ull, 16234865648979114395ull, 35478654369178396ull, 16879682265742406327ull, 988822640532990652ull, 6441252208187343363ull, 7471324117207928020ull, 9321993331818520759ull, 13530855453747965965ull, 16878852095836032681ull, 999692530134946768ull, 12507060968559488909ull, 1646244627744266941ull, 2863247263353257887ull, 17953355115889683125ull, 10752587360488671818ull, 10532860777350021446ull, 2364427453494076967ull, 7455771926569400730ull, 5195467989032989018ull, 1834988649296150770ull, 12020496442723221343ull, 13846846029514221882ull, 14513446473104247648ull, 17178375772904719250ull, 7603705109499324354ull, 4478326676333763172ull, 11816598820375785035ull, 15738093394210874391ull, 6652442946475616050ull, 11983808028397646472ull, 2973536586709636463ull, 12267436926127621622ull, 12148982840339802509ull, 15447503755701382904ull, 9170693773425904909ull };
	const uint64 y70[70] = { 13035174564231269173ull, 2031817238368888345ull, 14227779053044384239ull, 15697031109213648585ull, 16246797397621517264ull, 6332083408892383727ull, 9388765936565627601ull, 3055339951133711026ull, 17382999924880859970ull, 12771672880377584848ull, 8364171159194281015ull, 7076312770171612795ull, 2775074124666994714ull, 541377828095402822ull, 10583857457855614412ull, 2657209127254259452ull, 7971412721646259126ull, 15646619613237073924ull, 7193777446646360147ull, 13685152369053955928ull, 3288403155551002076ull, 16539390239681193976ull, 5892710676204249261ull, 7855521822881034680ull, 2976852341030657572ull, 7961678784864755097ull, 3979145662020335758ull, 16271488984939745189ull, 13690797836746321285ull, 11196731569817184326ull, 11885717758943982388ull, 17357160804673989517ull, 13804529479597294534ull, 16043721574971970400ull, 2096844946123222333ull, 9921975875987165637ull, 2768798678716382281ull, 12947331320923744074ull, 1084691708838543443ull, 10749401893225727149ull, 10768375991636447559ull, 3260458100291155455ull, 14413637195894780110ull, 15652838634222062970ull, 4048582062151630608ull, 8946406208257462796ull, 17853047967386028589ull, 11641662773640704651ull, 11303917129069076390ull, 17063452783458751347ull, 2049874958995099249ull, 15546160938745120266ull, 9538043553152286716ull, 18162221324682376697ull, 9628929869110154530ull, 6784665020163418060ull, 1378708342066296474ull, 11536753251347373535ull, 4534979510131407149ull, 17241580948451511362ull, 11374788804146286863ull, 6680542423372366105ull, 18262699497438852104ull, 2580645069957584047ull, 17399724939617397826ull, 13729671085246255165ull, 1262258216099317763ull, 8161435500821964714ull, 2568509089485667314ull, 13586984782852770896ull };
	const uint64 y80[80] = { 3686879699877952657ull, 453741885935125956ull, 3938341138369497622ull, 16419599976023570750ull, 17910950143366487185ull, 5210269363513904009ull, 6365406617924173421ull, 1137288991988458049ull, 7514907328692273864ull, 6142628874908786840ull, 11723456648662884184ull, 15315897270966342086ull, 12657568862368998877ull, 12168618695858229831ull, 7894411684859967070ull, 15459088869505907538ull, 874797005966515536ull, 8552702139089316709ull, 14992398079764467869ull, 7490232213726620779ull, 2531746272349789760ull, 12917674233239253212ull, 15625330359100556049ull, 9093415458570795535ull, 13538593178636244204ull, 9351204785170562746ull, 14303148925206010133ull, 8876380558379849305ull, 1375113876046730803ull, 6246506863689184382ull, 10848669702049448500ull, 5254421183383954712ull, 12524232939757077924ull, 14519354321980521579ull, 13978127845095476853ull, 2819464761450493294ull, 5706007102882530005ull, 4551246688993488148ull, 11434788483145254823ull, 16173841354543856846ull, 12509315208860406023ull, 4172036203402343042ull, 11799184975285764057ull, 12925189668663022489ull, 5747199812525005113ull, 8764377299440814754ull, 1428887035158948118ull, 2233649227163155683ull, 9403449734926482041ull, 14757393839014274516ull, 8304251627243296163ull, 1937119483778953015ull, 3038821567260813359ull, 12717451046176433525ull, 8660458500420561522ull, 16144546352425647343ull, 14737742462608257290ull, 8697636606013577263ull, 17293841164144319976ull, 15928064021931218237ull, 12115048999295683122ull, 3255001946126728228ull, 8723994184679982995ull, 13439587285723996107ull, 7669185267979994838ull, 17379309529588673134ull, 6846954217803119602ull, 12356564666511368161ull, 9997388457301655286ull, 8583845996545228489ull, 8699953439928823890ull, 15611507079256356466ull, 9868697811108222063ull, 2375318481144124919ull, 12135791452161853891ull, 17277216076912414267ull, 6137863990657693900ull, 6415440756785175233ull, 10505426288612333674ull, 96627997299058663ull };
	const uint64 y90[90] = { 14327012453873980177ull, 7896631668913294125ull, 10743171059613661959ull, 9254020797600318726ull, 7821195462223301752ull, 12067004296059438889ull, 11451145975088366539ull, 3356748286194447043ull, 5756955915065573524ull, 4967394651906704104ull, 16200443127959440649ull, 748986858458749867ull, 5935610987814150116ull, 7257751516081777223ull, 14632677804351388845ull, 10673271736460498614ull, 10598143576882860980ull, 5342508627593180883ull, 2289255404715875601ull, 18353036279070481307ull, 1462476853186759077ull, 3808265043807978242ull, 9872366095748063543ull, 10012349086639405479ull, 18005086452581661405ull, 16861687126321445049ull, 3294135959910501121ull, 5229794723077146894ull, 3794554444783317603ull, 12152833457765224981ull, 6372607249820535687ull, 3473297593761212178ull, 13771677630026388707ull, 345555404739096981ull, 12722074368935490792ull, 8379379876849799073ull, 10325518729628012158ull, 8698743719780584741ull, 1734303723533305594ull, 10717879912171351068ull, 18242600550810530527ull, 17323713309028296949ull, 552488179587813268ull, 3959798390820063107ull, 11914387527272528608ull, 8992644645129307207ull, 8962962353479785492ull, 3296326283812962026ull, 9706029827934915586ull, 6847717899506044476ull, 9393747518899414823ull, 15192326226324498214ull, 3713534964763546399ull, 9331349013517112684ull, 5279581406580265303ull, 15033055467318031785ull, 18014970481751534816ull, 16545256112998204669ull, 15574852613148839853ull, 3905640428107038313ull, 5804591361607298468ull, 14834959420280064335ull, 6522829411953094466ull, 9334674759583114244ull, 9293846761775796332ull, 16478175502621994545ull, 155420708107757789ull, 10526865276925211750ull, 77227172960405923ull, 15962009290802429017ull, 292405915149433336ull, 10802968607058396556ull, 17598447061487595913ull, 2597092985838943070ull, 5403151017845588797ull, 793021800535760345ull, 12220345118389884222ull, 11504865507696002441ull, 13842125103846679754ull, 14288554951294647296ull, 6646685620512546564ull, 11456722540097699043ull, 6371962759385722204ull, 18037500756150660222ull, 12372696248527741680ull, 17492587105238146773ull, 7029004413759659380ull, 9120104346702250711ull, 6438074115073266034ull, 8763185954741157488ull };
	const uint64 y100[100] = { 8651761500958765785ull, 18144327741091934514ull, 12270157337756595666ull, 3156099161148169142ull, 4279338241002989180ull, 12083717968964082031ull, 9633317538807782904ull, 4870798020911124321ull, 2981175300705172332ull, 2151188483678590934ull, 10874467985769776816ull, 9382306411753430533ull, 1036028362070479645ull, 715225889416666486ull, 17804187131157907299ull, 13617865426029361250ull, 14751831722551499909ull, 13260903344971084065ull, 1385748392519462027ull, 15263475433688704658ull, 2053361431173867803ull, 8865637984170331070ull, 14288141023742298787ull, 11164273184837158983ull, 7985730982309512025ull, 5702456664721074475ull, 6166119136247968676ull, 11050566233103716601ull, 17496346217083451029ull, 18255904514403111601ull, 12216469643673765517ull, 3204481586097061852ull, 11547311222167575704ull, 2092570891781513960ull, 10130609998198227026ull, 1334745727499238271ull, 11318496886609830379ull, 3746857903549267585ull, 18313219339809391475ull, 17625930300511294323ull, 9367046833890820376ull, 15656591948583642216ull, 14940178300211116794ull, 2126476221790779880ull, 5974469236764481950ull, 6112732136507428712ull, 13721023765697723568ull, 2058945967400097341ull, 71366495660924301ull, 6152375738491761391ull, 10890210021618721414ull, 14177711553690602706ull, 10008455012174366044ull, 6960104434898147378ull, 7195027169060283539ull, 12345827306126428934ull, 10177731113695244691ull, 6733615417560036667ull, 13374622979515710970ull, 12005838724605074293ull, 6824930360888017736ull, 6164640720528174496ull, 1722571809865283122ull, 8542194294112176484ull, 14588769469551781698ull, 13470679260057915358ull, 3544840540517842490ull, 45773415225067640ull, 14389667535480509479ull, 5665297163344087138ull, 11200764405566174782ull, 10255336710292127490ull, 9646756178641762584ull, 1816340593859371465ull, 11953993155776741801ull, 1311299289600998907ull, 17847651011707485923ull, 4862792314258313882ull, 251379236753398197ull, 16605144173203011556ull, 10895953446093503786ull, 9391485119640543730ull, 6675880166957346028ull, 5920536850090467440ull, 8456458013384015881ull, 1313125224932561774ull, 9654843907277435837ull, 3005807009564180859ull, 2367918741490012605ull, 13282919683781319355ull, 3686131585151930643ull, 14798551216059837626ull, 3718487697088598262ull, 1213842626237841938ull, 354332956691164894ull, 2416810242288740439ull, 18361392855425522495ull, 17083443342269785852ull, 17392421406996216228ull, 13453818413721276128ull };

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
	const uint64 *ptr_mword_prime[20], *curr_mword_prime;
	uint64 *out_array = 0x0;
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
	ptr_mword_prime[10] = y15;
	ptr_mword_prime[11] = y20;
	ptr_mword_prime[12] = y30;
	ptr_mword_prime[13] = y40;
	ptr_mword_prime[14] = y50;
	ptr_mword_prime[15] = y60;
	ptr_mword_prime[16] = y70;
	ptr_mword_prime[17] = y80;
	ptr_mword_prime[18] = y90;
	ptr_mword_prime[19] = y100;

	fprintf(stderr, "INFO: testing GCD routines...\n");

	/* Init the RNG: */
	rng_isaac_init(TRUE);

	/* Allocate the main data arrays: */
	sec_len = FFTMUL_THRESHOLD_BITS>>2;		// Input vecs to FFT-gcd must be suitably 0-padded at top:
//	j = sec_len - (MAX_ARRAY_DIM & (sec_len-1));	// mod-via-AND here assumes sec_len a power of 2
	j = sec_len;	// Use this simpler form to allow for pure-integer computation of a*u-b*v for check vs FFT result
	u = (uint64 *)calloc(j+MAX_ARRAY_DIM, sizeof(uint64));
	v = (uint64 *)calloc(j+MAX_ARRAY_DIM, sizeof(uint64));
	w = (uint64 *)calloc(2*MAX_ARRAY_DIM, sizeof(uint64));
	x = (uint64 *)calloc(j+MAX_ARRAY_DIM, sizeof(uint64));
	y = (uint64 *)calloc(j+MAX_ARRAY_DIM, sizeof(uint64));
	z = (uint64 *)calloc(2*MAX_ARRAY_DIM, sizeof(uint64));
	a = (double *)calloc(16*MAX_ARRAY_DIM, sizeof(double));
	b = (double *)calloc(16*MAX_ARRAY_DIM, sizeof(double));
//=====================================
#if 0
#warning Enable me to resume FFT-GCD debug!
	vec_len = 20000;	tmp64 = 4788416424359163823ull;// For len 10000 try j = 523 to expose cvt-cy bug
	for(j = 0; j < 1000; j++) {
		for(i = 0; i < vec_len-1; i++) {
			u[i] = rng_isaac_rand();	v[i] = rng_isaac_rand();
		}
	if(j < 0) continue;		// ...then stick the offending j into the (j < **) to skip GCD until hit the case in question.
	fft_gcd_debug = (j == -1);	// Set -1 to specific index you want to do detailed debug on
	printf("****** Trial %u ******\n",j);	// Initial trial-loop to find smallest j for which test fails...
	// Fix LSB = 1 to eliminate most common cause of the as-yet-unmultiplied-by-p vecs having GCD > 1, a common power of 2:
	u[0] |= 1ull;	v[0] |= 1ull;
	// Now mpy quasirandom 8999-word vecs by small prime GCD:
	u[vec_len-1] = mi64_mul_scalar(u, tmp64, u, vec_len-1);
	v[vec_len-1] = mi64_mul_scalar(v, tmp64, v, vec_len-1);
	// Do the FFT-GCD ... result may be a multiple of the shared small prime, need to account for that:
	lenX = fft_gcd(u,v, vec_len, eGCD = FALSE, 0x0,0x0,0x0,0x0);
	fail += !mi64_is_div_by_scalar64(u, tmp64, lenX);
	fail += !mi64_is_div_by_scalar64(v, tmp64, lenX);
	ASSERT(HERE, fail==0,"Failed testcase in FFT-GCD loop!");
}
exit(0);
#endif
//=====================================
	// Basic tests of mi64_add/sub:
	ASSERT(HERE, 1000 < MAX_ARRAY_DIM-2, "MAX_ARRAY_DIM too small!");
	for(i = 0; i < 1000; ++i) {
		for(j = 0; j <= i; ++j) {
			x[j] = u[j] = rng_isaac_rand();
			y[j] = v[j] = rng_isaac_rand();
		}
		z[j] = mi64_add_ref(x,y,z, i+1);
		w[j] = mi64_add    (u,v,w, i+1);
		if(!mi64_cmp_eq(w,z, i+2)) {
			for(j = 0; j <= i+1; ++j) {
				if(w[j] != z[j]) printf("w[%3d] = %20llu != z[%3d] = %20llu\n",j,w[j],j,z[j]);
			}
			ASSERT(HERE, 0, "mi64_add results differ!");
		}
	}

	/*
	Check the uint64[] <==> double[]  interconversion routines:
	on arrays of increasing size, doubling the size until reach MAX_ARRAY_DIM:
	*/
	vec_len = 1;
	while(vec_len <= MAX_ARRAY_DIM)
	{
		/* Init an array of (vec_len) uint64's using quasirandom selections from y4[]: */
		for(i = 0; i < vec_len; i++) {
			x[i] = rng_isaac_rand();
			y[i] = rng_isaac_rand();
			u[i] = x[i];
			v[i] = y[i];
		}
		ASSERT(HERE, (cy = mi64_cvt_uint64_double(x,y, 0, vec_len, a)) <= 3, "GCD self-test CVT: unexpected carryout of mi64_cvt_uint64_double!");
		// vec_len is for *each* 64-bit input vec, thus interleaved floating-double ovec has 8*vec_len terms.
		// Any carryouts of the mi64_cvt_uint64_double() call must be folded back into high words of a[]
		// prior to back-conversion:
		if((cy & 1) == 1) {
			i = 8*vec_len - 2;	i += ( (i >> DAT_BITS) << PAD_BITS );
			ASSERT(HERE, a[i] < 0, "Unexpected sign of a[i]!");
			a[i] += FFT_MUL_BASE;
		}
		if((cy >> 1) == 1) {
			i = 8*vec_len - 2;	i += ( (i >> DAT_BITS) << PAD_BITS );
			ASSERT(HERE, a[i+1] < 0, "Unexpected sign of a[i]!");
			a[i+1] += FFT_MUL_BASE;
		}
		ASSERT(HERE, mi64_cvt_double_uint64(a, 4*vec_len, x,y) == 0, "GCD self-test CVT: unexpected carryout of mi64_cvt_double_uint64!");
		for(i = 0; i < vec_len; i++) {
			ASSERT(HERE, x[i] == u[i], "GCD self-test CVT: x[i] != u[i]");
			ASSERT(HERE, y[i] == v[i], "GCD self-test CVT: y[i] != v[i]");
		}
		vec_len *= 2;
		/* Want to do at least one test right at the max. length: */
		if(vec_len > MAX_ARRAY_DIM) {
			if(vec_len >= 2*MAX_ARRAY_DIM)
				break;
			else
				vec_len = MAX_ARRAY_DIM;
		}
	}

#if 0	// Set == 1 to enable mi64_add timing test
	// mi64_add timing test
	printf	("GCD: Starting mi64_add() timing test loop:\n");
	vec_len = 1000;
	for(j = 0; j <= i; ++j) {
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
	exit(0);
#endif

// 8/08/2012: Track down mi64_div bug which leaves remainder with 1*divisor needing subtraction:
	// x= 184027369222821018387581692039:
	x[1] = 9976143675ull;	x[0] =  7539741246537263239ull;
	y[1] =         10ull;	y[0] = 10733642766940458262ull;
	mi64_div(x,y,2,2,u,v);	// quo in u[], rem in v[]
	ASSERT(HERE, u[1]==0 && u[0]==942757929 && v[1]==0 && v[0]==1, "divisor check fails!");

// 6/08/2012: Trial-factor of mm127 residue check for Paul L:
#if 0
	p = 127;	lenX = (p>>6);
	memset(x,0xff,(lenX<<3));	x[lenX++] = (1ull << (p&63)) - 1;	// x has p
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
	/*
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
	
	So q == 7 or 9 (mod 16) gives 5 good bits in qinv_0 ==> 40 good bits in q*qinv_3 product computed on final iteration.
	   q == 1 of F (mod 16) gives >= 6 bits in qinv_0, with odds of x good bits roughly halving with each increment in x.
	
	But what about the "missing moduli", i.e. q == 3,5,B,D (mod 16)?
	
	Tried those by twiddling the bottom 4 bits of the same dataset, all give >= 5 bits.
	But interestingly, the mod-16 values which give "at least 5" (3,7,9,13) *never* seem to give more than 5,
	whereas the ones which give "at least 6" (1,5,11,15) all give 6,7,8,... good bits with probability
	halving with each additional "lucky" bit of precision.
	*/
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
	
	/**** 5/16/2012: Tmp-code to test Montgomery-mul remainder for M977, see if can efficiently convert
	carryout of my RL divisibility algo to true remainder:
	*/
	p = 977;	lenX = (p>>6);
	memset(x,0xff,(lenX<<3));	x[lenX++] = (1ull << (p&63)) - 1;
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
	mi64_shrl_short(x, y, 1, lenX);	// y = x>>1
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
/*
b=2^64;t=1;x=2421895350930534175*t;t*=b;x+=2416981188099281967*t;t*=b;x+=17119718852038027230*t;t*=b;x+=792133985803089581*t;t*=b;x+=2967675029189239758*t;t*=b;x+=3104517425838997392*t;t*=b;x+=16792572262743320*t;t*=b;x+=1917901379003983519*t;t*=b;x+=15810729283012914840*t;t*=b;x+=10466991815671365502
*/
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

	/************* Test DIV1: shakedown tests of non-folded DIV algo: ****************/
	printf("Performing GCD self-test DIV1...\n");
	vec_len = 1;
	while(vec_len <= 1024)
	{
		for(j=0;j<100;++j)
		{
			/* Init an array of (vec_len) quasirandom 64-bit words: */
			for(i = 0; i < vec_len; i++)
			{
				x[i] = y[i] = rng_isaac_rand();
				tmp64 = rng_isaac_rand();
			}
			x[vec_len] = mi64_mul_scalar(x, tmp64, x, vec_len);
			ASSERT(HERE, mi64_is_div_by_scalar64(x, tmp64, vec_len+1) == TRUE, "GCD self-test #DIV1a!");
			rem = mi64_div_by_scalar64(x, tmp64, vec_len+1, x);
			ASSERT(HERE, rem == 0ull && mi64_cmp_eq(x,y,vec_len) == TRUE, "GCD self-test #DIV1b!");
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

	/************* Test DIV2: shakedown tests of loop-folded DIV algo: ****************/
	printf("Performing GCD self-test DIV2...\n");	fflush(stdout);
	vec_len = MAX_ARRAY_DIM;
	for(i = 0; i < vec_len; i++) {
		x[i] = rng_isaac_rand();
	}

	// Correctness tests of short (< 1-word) shifts:
	// rshift:
	for(vec_len = 1; vec_len < 1000; vec_len++) 	{
		tmp64 = x[vec_len];
		j =  tmp64     & 0x3f;	// 63 bits = max shift count for short-shift routines
		kk= (tmp64>>8) & 0x3f;	// Also randomly move start-word around between x[0-63] to test alignment-based code
		// We do in-place shifts below, so copy x-vectar data into working arrays y and z:
		memcpy(y+kk,x,vec_len<<3);	memcpy(z+kk,x,vec_len<<3);	
		tmp64 = mi64_shrl_short_ref(y+kk,y+kk,j,vec_len);
		tmp   = mi64_shrl_short    (z+kk,z+kk,j,vec_len);
		ASSERT(HERE, tmp == tmp64 && mi64_cmp_eq(z+kk,y+kk,vec_len), "short-rshift test fails!");
	}
	// lshift:
	for(vec_len = 1; vec_len < 1000; vec_len++) 	{
		tmp64 = x[vec_len];
		j =  tmp64     & 0x3f;	// 63 bits = max shift count for short-shift routines
		kk= (tmp64>>8) & 0x3f;	// Also randomly move start-word around between x[0-63] to test alignment-based code
		// We do in-place shifts below, so copy x-vectar data into working arrays y and z:
		memcpy(y+kk,x,vec_len<<3);	memcpy(z+kk,x,vec_len<<3);	
		tmp64 = mi64_shl_short_ref(y+kk,y+kk,j,vec_len);
		tmp   = mi64_shl_short    (z+kk,z+kk,j,vec_len);
		ASSERT(HERE, tmp == tmp64 && mi64_cmp_eq(z+kk,y+kk,vec_len), "short-lshift test fails!");
	}

	// Performance tests of short (< 1-word) shifts:
	vec_len = 1000;
	memcpy(y,x,vec_len<<3);	ASSERT(HERE, y[0] != 0 && y[0] == x[0], "init error!");
	j = 1000000;
	// I/0 vectors have same alignment:
	clock1 = clock();
	tmp = 0ull;
	for(i = 0; i < j; i++) 	{
		tmp64 = mi64_shrl_short(y,z,63,vec_len);	// 63 bits = max shift count for short-shift routines
		tmp  += mi64_shl_short (z,y,63,vec_len);	y[0] += (tmp64 >> 1);	// Add rshift off-shift bits back into low word, cy = 0 assured
	}
	ASSERT(HERE, mi64_cmp_eq(x,y,vec_len), "short-shift test fails!");
	clock2 = clock();
	tdiff = (double)(clock2 - clock1);
	printf("64-bit short-shift, #words = %llu, Time =%s, (sum = %llu)\n",(uint64)j*vec_len*2,get_time_str(tdiff),tmp);

	// I/0 vectors have 8-byte-offset alignment:
	clock1 = clock();
	tmp = 0ull;
	for(i = 0; i < j; i++) 	{
		tmp64 = mi64_shrl_short(y,z+1,63,vec_len);	// 63 bits = max shift count for short-shift routines
		tmp  += mi64_shl_short (z+1,y,63,vec_len);	y[0] += (tmp64 >> 1);	// Add rshift off-shift bits back into low word, cy = 0 assured
	}
	ASSERT(HERE, mi64_cmp_eq(x,y,vec_len), "short-shift test fails!");
	clock2 = clock();
	tdiff = (double)(clock2 - clock1);
	printf("64-bit short-shift, 8-byte-offset, Time =%s, (sum = %llu)\n",get_time_str(tdiff),tmp);

	// 64-bit Montgomery inverse, several versions:
	q = x[0];	q |= 1;
	clock1 = clock();
	tmp = 0ull;
	j = 100000000;
	while(j--) {
		qinv = (q+q+q) ^ (uint64)2;
		for(i = 0; i < 4; i++) {
			qinv = qinv*((uint64)2 - q*qinv);
		}
		q += 2;	tmp += qinv;
	}	// Core2: ~19.8 cycles; no faster unrolling the inner 4-pass loop
	clock2 = clock();
	tdiff = (double)(clock2 - clock1);
	printf("64-bit Mont-inverse v1, #calls = %u, Time =%s, (sum = %llu)\n",j,get_time_str(tdiff),tmp);

	q = x[0];	q |= 1;
	clock1 = clock();
	tmp = 0ull;
	j = 100000000;
	while(j--) {
		qinv = minv8[(q&0xff)>>1];
		for(i = 0; i < 3; i++) {
			qinv = qinv*((uint64)2 - q*qinv);
		}
		q += 2;	tmp += qinv;
	}	// Core2: ~14.3 cycles, 5.5 cycles faster than above v1
	clock2 = clock();
	tdiff = (double)(clock2 - clock1);
	printf("64-bit Mont-inverse v2, #calls = %u, Time =%s, (sum = %llu)\n",j,get_time_str(tdiff),tmp);

	uint32 q32,qi32;
	q = x[0];	q |= 1;
	clock1 = clock();
	tmp = 0ull;
	j = 100000000;
	while(j--) {
		q32  = q; qi32 = minv8[(q&0xff)>>1];
		qi32 = qi32*((uint32)2 - q32*qi32);
		qi32 = qi32*((uint32)2 - q32*qi32);	qinv = qi32;
		qinv = qinv*((uint64)2 - q*qinv);
		q += 2;	tmp += qinv;
	}	// Core2: ~12 cycles, ~8 cycles faster than v1
	clock2 = clock();
	tdiff = (double)(clock2 - clock1);
	printf("64-bit Mont-inverse v3, #calls = %u, Time =%s, (sum = %llu)\n",j,get_time_str(tdiff),tmp);

	// Correctness (w.r.to the internal fault handlers) tests of the radix_power64() routine,
	// using random 32/48/64-bit moduli (by which we mean < 2^32/48/64, respectively):
	kk = 10000000;
	// q < 2^32:
	clock1 = clock();
	tmp = 0ull;
	j = kk;
	while(j--) {
		q = rng_isaac_rand() >> 32; q |= 1;	// Force odd
		qinv = (q+q+q) ^ (uint64)2;
		for(i = 0; i < 4; i++) {
			qinv = qinv*((uint64)2 - q*qinv);
		}
		tmp += radix_power64(q,qinv,2);
	}
	clock2 = clock();
	tdiff = (double)(clock2 - clock1);
	printf("Shakedown-test of RADIX_POWER64, q < 2^32, #calls = %u (sum = %llu), Time =%s\n",kk,tmp,get_time_str(tdiff));

	// q < 2^48:
	clock1 = clock();
	tmp = 0ull;
	j = kk;
	while(j--) {
		q = rng_isaac_rand() >> 16; q |= 1;	// Force odd
		qinv = (q+q+q) ^ (uint64)2;
		for(i = 0; i < 4; i++) {
			qinv = qinv*((uint64)2 - q*qinv);
		}
		tmp += radix_power64(q,qinv,2);
	}
	clock2 = clock();
	tdiff = (double)(clock2 - clock1);
	printf("Shakedown-test of RADIX_POWER64, q < 2^48, #calls = %u (sum = %llu), Time =%s\n",kk,tmp,get_time_str(tdiff));

	// q < 2^64:
	clock1 = clock();
	tmp = 0ull;
	j = kk;
	while(j--) {
		q = rng_isaac_rand(); q |= 1;	// Force odd
		qinv = (q+q+q) ^ (uint64)2;
		for(i = 0; i < 4; i++) {
			qinv = qinv*((uint64)2 - q*qinv);
		}
		tmp += radix_power64(q,qinv,2);
	}
	clock2 = clock();
	tdiff = (double)(clock2 - clock1);
	printf("Shakedown-test of RADIX_POWER64, q < 2^64, #calls = %u (sum = %llu), Time =%s\n",kk,tmp,get_time_str(tdiff));

	// Timing tests of key parts of the radix_power64() routine:
	q = x[0];	q |= 1;
	qinv = (q+q+q) ^ (uint64)2;
	for(i = 0; i < 4; i++) {
		qinv = qinv*((uint64)2 - q*qinv);
	}

	// 64-bit stdlib mod:
	j = 100000000;
	clock1 = clock();
	tmp = 0ull;
	for(i = 0; i < j; i++) 	{
		tmp += x[i&0xffff] % q;
	}
	clock2 = clock();
	tdiff = (double)(clock2 - clock1);
	printf("64-bit stdlib mod, #calls = %u, Time =%s, (sum = %llu)\n",j,get_time_str(tdiff),tmp);

	// MONT_SQR64() macro:
	clock1 = clock();
	tmp = 0ull;	// Accumulate sum of outputs in loop to prevent compiler from optimizing away the macro as a no-op:
	for(i = 0; i < j; i++) 	{
		MONT_SQR64(x[1],q,qinv,rem);	tmp += rem;
	}
	clock2 = clock();
	tdiff = (double)(clock2 - clock1);
	printf("MONT_SQR64, #calls = %u, Time =%s, (sum = %llu)\n",j,get_time_str(tdiff),tmp);

	// 64-bit odd modulus q, B^2 mod q:
	clock1 = clock();
	for(i = 0; i < j; i++) 	{
		tmp = radix_power64(q,qinv,2);
	}
	clock2 = clock();
	tdiff = (double)(clock2 - clock1);
	printf("RADIX_POWER64, B^2 mod q, 64-bit modulus, #calls = %u, Time =%s\n",j,get_time_str(tdiff));

	// 64-bit odd modulus q, B^1000 mod q:
	j = 10000000;	// ~10x slower than B^2, so cut loop count by 10
	clock1 = clock();
	for(i = 0; i < j; i++) 	{
		tmp = radix_power64(q,qinv,1000);
	}
	clock2 = clock();
	tdiff = (double)(clock2 - clock1);
	printf("RADIX_POWER64, B^1000 mod q, 64-bit modulus, #calls = %u, Time =%s\n",j,get_time_str(tdiff));

	// 64-bit odd modulus q, B^1000000 mod q:
	j = 5000000;	// ~2x slower than B^1000, so cut loop count by 2
	clock1 = clock();
	for(i = 0; i < j; i++) 	{
		tmp = radix_power64(q,qinv,1000000);
	}
	clock2 = clock();
	tdiff = (double)(clock2 - clock1);
	printf("RADIX_POWER64, B^1000000 mod q, 64-bit modulus, #calls = %u, Time =%s\n",j,get_time_str(tdiff));

	tmp = 0ull;	// array-length accumulator
	// Use a smaller max. length for this time-consuming series of correctness tests:
	for(i = 1; i < 1000; i++)
	{
		tmp += i;
		mod = rng_isaac_rand();
	//	printf("Len = %u, q = %16llu...\n",i,mod);
		rem   = mi64_div_by_scalar64   (x, mod, i, y);	mi64_mul_scalar(y, mod, z, i);	mi64_add_scalar(z, rem, z, i);	ASSERT(HERE, mi64_cmp_eq(x, z, i), "GCD self-test #DIV2_u1");
		tmp64 = mi64_div_by_scalar64_u2(x, mod, i, z);	ASSERT(HERE, rem == tmp64 && mi64_cmp_eq(y,z,i), "GCD self-test #DIV2_u2");
		tmp64 = mi64_div_by_scalar64_u4(x, mod, i, z);	ASSERT(HERE, rem == tmp64 && mi64_cmp_eq(y,z,i), "GCD self-test #DIV2_u4");
	}
	printf("DIVREM u1|u2|u4 correctness self-test passed. Total vector Len = %llu\n",tmp);

	// Remainder-only, odd moduli:
	j = 300000;	//*** Note that the following 4 tests have work increasing with the *square* of the loop bound
	clock1 = clock();
	tmp = 0ull;	// array-length accumulator
	for(i = 4; i < j; i += 40 /*(((i>>2)+3) & 0xfffffffc)*/)	// Increment must keep i a multiple of 4
	{
		tmp += i;
		mod = rng_isaac_rand();
		mod |= 1;	// Force odd  modulus
	//	printf("Len = %u, q = %16llu...\n",i,mod);
		mi64_div_by_scalar64_u4(x, mod, i, 0x0);
	}
	clock2 = clock();
	tdiff = (double)(clock2 - clock1);
	printf("\tREM 4-fold performance test, odd -moduli: Total vector Len = %llu, Time =%s\n",tmp,get_time_str(tdiff));

	// Remainder-only, even moduli:
	clock1 = clock();
	tmp = 0ull;	// array-length accumulator
	for(i = 4; i < j; i += 40 /*(((i>>2)+3) & 0xfffffffc)*/)	// Increment must keep i a multiple of 4
	{
		tmp += i;
		mod = rng_isaac_rand();
		mod &= 0xFFFFFFFFFFFFFFFE;	// Force even modulus
	//	printf("Len = %u, q = %16llu...\n",i,mod);
		mi64_div_by_scalar64_u4(x, mod, i, 0x0);
	}
	clock2 = clock();
	tdiff = (double)(clock2 - clock1);
	printf("\tREM 4-fold performance test, even-moduli: Total vector Len = %llu, Time =%s\n",tmp,get_time_str(tdiff));

	// Remainder and quotient, odd moduli:
	clock1 = clock();
	tmp = 0ull;	// array-length accumulator
	for(i = 4; i < j; i += 40 /*(((i>>2)+3) & 0xfffffffc)*/)	// Increment must keep i a multiple of 4
	{
		tmp += i;
		mod = rng_isaac_rand();
		mod |= 1;	// Force odd  modulus
	//	printf("Len = %u, q = %16llu...\n",i,mod);
		mi64_div_by_scalar64_u4(x, mod, i, z);
	}
	clock2 = clock();
	tdiff = (double)(clock2 - clock1);
	printf("DIV+REM 4-fold performance test, odd -moduli: Total vector Len = %llu, Time =%s\n",tmp,get_time_str(tdiff));

	// Remainder and quotient, even moduli:
	clock1 = clock();
	tmp = 0ull;	// array-length accumulator
	for(i = 4; i < j; i += 40 /*(((i>>2)+3) & 0xfffffffc)*/)	// Increment must keep i a multiple of 4
	{
		tmp += i;
		mod = rng_isaac_rand();
		mod &= 0xFFFFFFFFFFFFFFFE;	// Force even modulus
	//	printf("Len = %u, q = %16llu...\n",i,mod);
		mi64_div_by_scalar64_u4(x, mod, i, z);
	}
	clock2 = clock();
	tdiff = (double)(clock2 - clock1);
	printf("DIV+REM 4-fold performance test, even-moduli: Total vector Len = %llu, Time =%s\n",tmp,get_time_str(tdiff));
	/********************************************************************************/
exit(0);
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
	Test #2: let p = 16357897499336320021 (64-bit test prime), check that mi64_pprimeF
	to various small-prime bases returns TRUE:
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
#if 0
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
#endif
	/*
	Test #5: let p = 16357897499336320049 [next-higher 64-bit prime above our 1-word reference prime]
	x = p * {some quasi-random (N-1)-word int}
	y =     {some quasi-random  N   -word prime}, N of various lengths for which we have prestored multiword primes.
	*/
	printf("\nPerforming GCD self-test #5...\n");
	for(j = 1; j < num_mword_prime; j++)	/* Start multipliers with 2-word primes and go up from there */
	{
	#if GCD_DEBUG >= 1
		if(gcd_debug) printf("\n********************************\n      Length-%u inputs:\n********************************\n",j+1);
	#endif
		vec_len = len_mword_prime[j]-1;
		curr_mword_prime = ptr_mword_prime[j];
		for(i = 0; i < vec_len; i++) {
			x[i] = rng_isaac_rand();
		}
		x[vec_len] = mi64_mul_scalar(x, 16357897499336320049ull, x, vec_len);
		ASSERT(HERE, x[vec_len] != 0, "GCD self-test #5: expected nonzero carryout of mi64_mul_scalar!");
		vec_len += 1;
		// Toggle this #if value from 1 to 0 to enable code to auto-generate more fixed-length base-2 PRPs:
	#if 1
		for(i = 0; i < vec_len; i++) {
			y[i] = curr_mword_prime[i];
		}
	#else
		// Start with a random n-word int, fiddle LSB to make odd, then sub-2 until find a 2-PRP:
		for(i = 0; i < vec_len; i++) {
			y[i] = rng_isaac_rand();
		}	y[0] |= 1ull;
		uint32 max_try = 10000, curr_p, small_div, ntry = 0;
		while(ntry < max_try) {
			mi64_sub_scalar(y, 1ull, z, vec_len);	// z = y-1
			small_div = 0;	// First trial-divide by the small primes (2-PRPs, actually):
			for(curr_p = 3; curr_p < 10000*vec_len; curr_p += 2) {
				if(!pprimeF(curr_p,2)) continue;
				if(mi64_is_div_by_scalar64(y, curr_p, vec_len)) {
					small_div = 1;	//	printf("Found small divisor %u\n",curr_p);
					break;
				}
			}
			if(!small_div && ntry++ && mi64_twopmodq(z, vec_len, 0, y, vec_len, 0x0) == 1) {
				printf("[%u]th try: Found %u-word 2-PRP = { ",ntry,vec_len);	for(i = 0; i < vec_len-1; i++) { printf("%lluull, ",y[i]); }	printf("%lluull };\n",y[i]);
				break;
				// This requires y to have no more digits than STR_MAX_LEN = 1024, i.e. ~50 words:
			//	printf("[%u]th try: Found %u-word 2-PRP = %s\n",i,vec_len,&string0[convert_mi64_base10_char(string0, y,vec_len, 0)]);
			}
			mi64_sub_scalar(y, 2ull, y, vec_len);	// y -= 2
		}
		ASSERT(HERE, i < max_try, "No 2-PRP found!");
	#endif
//	gcd_debug = (vec_len == 100);
	if(gcd_debug) {
		// Used this testcase to convince myself that really do need as many trailing 'guard bits' in the leading-order
		// arrays sections being fed to the eGCD routine as bits in the resulting abcd multipliers to get a resulting
		// array-length reduction matching the length of the multipliers:
		// The abcd-mults resulting from this reduced length from lg(u,v) = 6398 to lg(|a*u-b*v|,|c*u-d*v|) = 3233:
		fail += (mi64_gcd(x, y, vec_len, eGCD = TRUE, Ap,Bp,&len_AB,&sign_AB, HALF = TRUE,Cp,Dp,&len_CD,&sign_CD,vec_len/2) != 1);
		// The abcd-mults resulting from this reduced length from log(u,v) = 6398 to log() = 5712, -686 bits vs -3165:
		fail += (mi64_gcd(x+40, y+40, 60, eGCD = TRUE, Ap,Bp,&len_AB,&sign_AB, HALF = TRUE,Cp,Dp,&len_CD,&sign_CD,50) != 1);
		exit(0);
	} else
		fail += (mi64_gcd(x, y, vec_len, eGCD = FALSE, 0x0,0x0,0x0,0x0, 0,0x0,0x0,0x0,0x0,0) != 1);
		fail += (x[0] != 1ull);	/* Put this on a separate line because some compilers will return 1 as the function result of evaluating (a+b), i.e. fail += a+b increments fail by 1 even if a = b = 0. Strange but true... */
		ASSERT(HERE, fail==0,"GCD self-test #5");
	}

	/*
	Test #6: let p = 2-word prime
	y = p * {some quasi-random N-word prime, stored in y}
	x = p * {some quasi-random N-word int coprime to y[]}
	*/
	printf("\nPerforming GCD self-test #6...\n");
	for(j = 2; j < 10; j++)	/* Start multipliers with 3-word primes and go up from there */
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
		fail += (mi64_gcd(x, y, vec_len, eGCD = FALSE, 0x0,0x0,0x0,0x0, 0,0x0,0x0,0x0,0x0,0) != 2);
		fail += (mi64_cmp_eq(x, y2, 2) != TRUE);
		ASSERT(HERE, fail==0,"GCD self-test #6");
	}

	/* Test #7: Use various small primes to test out the mi64_twopmodq routine: */
	printf("\nPerforming GCD self-test #7...\n");
	for(j = 0; j < 10; j++) {
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
	memset(x,0xff,(lenX<<3));	x[lenX-1] = (1ull << (p&63)) - 1;
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
	memset(y,0xff,(lenY<<3));	y[lenY-1] = (1ull << (p&63)) - 1;
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
	memset(y,0xff,(lenY<<3));	y[lenY-1] = (1ull << (p&63)) - 1;
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
	memset(y,0xff,(lenY<<3));	y[lenY-1] = (1ull << (p&63)) - 1;
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
	memset(y,0xff,(lenY<<3));	y[lenY-1] = (1ull << (p&63)) - 1;
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
	for(i = 0; i < lenV; i++) {
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
	vec_len = mi64_gcd(x, y, lenX, eGCD = FALSE, 0x0,0x0,0x0,0x0, 0,0x0,0x0,0x0,0x0,0);

	/* Note that 1001 is in fact composite (7*11*13) so M1001 has many small factors. Thus, since the
	random multiplier we used to pad the large-prime-factor product ma well share some of those small multipliers,
	we take a further GCD with the unpadded factor product:
	*/
	lenX = vec_len;
	mi64_setlen(u, lenU, lenX);
	vec_len = mi64_gcd(x, u, lenX, eGCD = FALSE, 0x0,0x0,0x0,0x0, 0,0x0,0x0,0x0,0x0,0);
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
	memset(x,0xff,(vec_len<<3));	x[vec_len-1] = (1ull << (p&63)) - 1;
	memset(u,(uint64) 0,(vec_len<<3));
	/* Store q in u[]: */
	u[0] =  4235226679561594903ull;
	u[1] =    94004235929829273ull;
	ASSERT(HERE, mi64_pprimeF(u, 3ull, 2) == 1, "GCD self-test #9.prp");
	clock1 = clock();
	gcd_debug=0;
	fail += (mi64_gcd(x, u, vec_len, eGCD = FALSE, 0x0,0x0,0x0,0x0, 0,0x0,0x0,0x0,0x0,0) != 2);
	fail += (x[0] !=  4235226679561594903ull);
	fail += (x[1] !=    94004235929829273ull);
	ASSERT(HERE, fail==0,"GCD self-test #9a");

	clock2 = clock();
	tdiff = (double)(clock2 - clock1);
	printf	("GCD self-test #9a passed: GCD Time =%s\n",get_time_str(tdiff));

	// Now re-init x and u and do full-length gcd:
	mi64_set_eq(u, x, vec_len);
	memset(x,0xff,(vec_len<<3));	x[vec_len-1] = (1ull << (p&63)) - 1;
	// Replace RHS with y431_test[i] here to recreate garbage-repeat-high-words testcase of 5/13/2012
	for(i = 0; i < 431; i++) {
		v[i] = y431_test[i];//rng_isaac_rand();
	//	printf("v[%3u] = %20llu\n",i,v[i]);
	}
	/* vector*vector MUL to get y is here: */
	mi64_mul_vector(u, 2, v, 431, y, &vec_len);
	ASSERT(HERE, vec_len == 433, "GCD self-test #9: vec_len == 433");
	clock1 = clock();
	gcd_debug=0;
	fail += (mi64_gcd(x, y, vec_len, eGCD = FALSE, 0x0,0x0,0x0,0x0, 0,0x0,0x0,0x0,0x0,0) != 2);
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
	memset(x,0xff,(vec_len<<3));	x[vec_len-1] = (1ull << (p&63)) - 1;
	memset(u,(uint64) 0,(vec_len<<3));	u[0] = 39654784768949071ull;
	clock1 = clock();
	gcd_debug=0;
	fail += (mi64_gcd(x, u, vec_len, eGCD = FALSE, 0x0,0x0,0x0,0x0, 0,0x0,0x0,0x0,0x0,0) != 1);
	fail += (x[0] !=    39654784768949071ull);
	ASSERT(HERE, fail==0,"GCD self-test #10a");

	clock2 = clock();
	tdiff = (double)(clock2 - clock1);
	printf	("GCD self-test #10a passed: GCD Time =%s\n",get_time_str(tdiff));

	// Now re-init x and u and do full-length gcd:
	memset(x,0xff,(vec_len<<3));	x[vec_len-1] = (1ull << (p&63)) - 1;
	for(i = 0; i < vec_len-1; i++) {
		y[i] = rng_isaac_rand();
	}
	y[vec_len-1] = mi64_mul_scalar(y, 39654784768949071ull, y, vec_len-1);
	ASSERT(HERE, y[vec_len-1] != 0, "GCD self-test #10: expected nonzero carryout of q0* y4[0:vec_len-1]!");

	clock1 = clock();
	gcd_debug=0;
	fail += (mi64_gcd(x, y, vec_len, eGCD = FALSE, 0x0,0x0,0x0,0x0, 0,0x0,0x0,0x0,0x0,0) != 1);
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
	memset(u,0xff,(vec_len<<3));	u[vec_len-1] = (1ull << (p&63)) - 1;
	/* Multiply u by random multiword integer to get 2001-word product: */
	for(i = 0; i < 2001-vec_len; i++) {
		v[i] = rng_isaac_rand();
	}

	/* vector*vector MUL to get x is here: */
	mi64_mul_vector(u, vec_len, v, 2001-vec_len, x, &vec_len);
	ASSERT(HERE, vec_len == 2001, "GCD self-test #11: vec_len == 2001");

	/* y: */
	for(i = 0; i < 1000; i++) {
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
	vec_len = mi64_gcd(x, y, vec_len, eGCD = FALSE, 0x0,0x0,0x0,0x0, 0,0x0,0x0,0x0,0x0,0);
	if(vec_len != 1 || x[0] != 39654784768949071ull) {
		mi64_clear(y, vec_len);
		y[0] = 39654784768949071ull;
		fail += (mi64_gcd(x, y, vec_len, eGCD = FALSE, 0x0,0x0,0x0,0x0, 0,0x0,0x0,0x0,0x0,0) != 1);
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
	y-vector by multiplying the factor by 34419 randomly selected words of y4[]. Owing to
	the special form of M-number divisors it's exceedingly unlikely that the resulting y[]
	will have any factors in common with x[] other than q.

	04 May 2012: ~9 sec on a single core of a 2 GHz Core2 Duo.
	23 Dec 2015: 6.2 sec	 "			"				" using new FFT-enhanced GCD.
	*/
#if OS_BITS == 64
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
		memset(u,0xff,(vec_len<<3));	u[vec_len-1] = (1ull << (p&63)) - 1;
	
	#if 0	// Set = 1 to do 64-bit is-divisible timing test
		clock1 = clock();
		for(i = 0; i < 10000; i++) {
			ASSERT(HERE, mi64_is_div_by_scalar64_u4(u, 87722769297534671ull, vec_len), "GCD self-test #12.0");
		}
		clock2 = clock();
		tdiff = (double)(clock2 - clock1);
		printf	("GCD self-test #12.0 passed: mi64_is_div_by_scalar64(), %u*%u words processed, Time =%s\n",i,vec_len,get_time_str(tdiff));
	//	exit(0);
	#endif
	#if 0	// Set = 1 to do 64-bit div-with-remainder timing test
		clock1 = clock();
		for(i = 0; i < 10000; i++) {
		//	ASSERT(HERE, mi64_div_by_scalar64   (u, 16357897499336320049ull, vec_len, v) == 14995895313315469881ull, "GCD self-test #12.0");
		//	ASSERT(HERE, mi64_div_by_scalar64_u2(u, 16357897499336320049ull, vec_len, v) == 14995895313315469881ull, "GCD self-test #12.0");
			ASSERT(HERE, mi64_div_by_scalar64_u4(u, 16357897499336320049ull, vec_len, v) == 14995895313315469881ull, "GCD self-test #12.0");
		}
		clock2 = clock();
		ASSERT(HERE,!mi64_mul_scalar(v, 16357897499336320049ull, v, vec_len), "GCD self-test #12.1");
		ASSERT(HERE, (v[0] += 14995895313315469881ull) && mi64_cmp_eq(u, v, vec_len), "GCD self-test #12.1");
		tdiff = (double)(clock2 - clock1);
		printf	("GCD self-test #12.1 passed: mi64_div_by_scalar64(), %u*%u words processed, Time =%s\n",i,vec_len,get_time_str(tdiff));
		exit(0);
	#endif
	
		// Next try a gcd of 2^p-1 (huge) with just the small factor:
		memset(v,(uint64) 0,(vec_len<<3));	v[0] = 87722769297534671ull;
		clock1 = clock();
		gcd_debug=0;
		fail += (mi64_gcd(u, v, vec_len, eGCD = FALSE, 0x0,0x0,0x0,0x0, 0,0x0,0x0,0x0,0x0,0) != 1);
		fail += (u[0] !=    87722769297534671ull);
		ASSERT(HERE, fail==0,"GCD self-test #12a");
	
		clock2 = clock();
		tdiff = (double)(clock2 - clock1);
		printf	("GCD self-test #12a passed: GCD Time =%s\n",get_time_str(tdiff));
	
		/* Now create a full-length v-vector = [large random padding integer]*[known small factor]: */
		memset(u,0xff,(vec_len<<3));	u[vec_len-1] = (1ull << (p&63)) - 1;
		for(i = 0; i < vec_len-1; i++) {
			v[i] = rng_isaac_rand();
		}
		v[vec_len-1] = mi64_mul_scalar(v, 87722769297534671ull, v, vec_len-1);
		ASSERT(HERE, v[vec_len-1] != 0, "GCD self-test #12: expected nonzero carryout of q0* y4[0:vec_len-1]!");

	//******************** Dec 2015: Use this testcase to prototype the FFT-mul-enhanced GCD: *************************
		clock1 = clock();
		gcd_debug=0;
	// Do the FFT-GCD:
		fail  = (fft_gcd(u, v, vec_len, eGCD = FALSE, 0x0,0x0,0x0,0x0) != 1);
		fail += (u[0] !=    87722769297534671ull);
		ASSERT(HERE, fail==0,"GCD self-test #12b");
		clock2 = clock();
		tdiff = (double)(clock2 - clock1);
		printf	("GCD self-test #12fft passed: FFT-GCD Time =%s\n",get_time_str(tdiff));

	// Pure-int GCD of same inputs:
		p = 2202817;	vec_len = (p>>6)+1;
		memset(u,0xff,(vec_len<<3));	u[vec_len-1] = (1ull << (p&63)) - 1;
		for(i = 0; i < vec_len-1; i++) {
			v[i] = rng_isaac_rand();
		}
		v[vec_len-1] = mi64_mul_scalar(v, 87722769297534671ull, v, vec_len-1);
		ASSERT(HERE, v[vec_len-1] != 0, "GCD self-test #12: expected nonzero carryout of q0* y4[0:vec_len-1]!");
		clock1 = clock();
		gcd_debug=0;
		fail  = (mi64_gcd(u, v, vec_len, eGCD = FALSE, 0x0,0x0,0x0,0x0, 0,0x0,0x0,0x0,0x0,0) != 1);
		fail += (u[0] !=    87722769297534671ull);
		ASSERT(HERE, fail==0,"GCD self-test #12c");
		clock2 = clock();
		tdiff = (double)(clock2 - clock1);
		printf	("GCD self-test #12int passed: INT-GCD Time =%s\n",get_time_str(tdiff));

	}

#endif	// OS_BITS == 64 ?

	/*
	Test #13: M22030163 has factor 795869615818103511527145535321.
	This M-number needs 344222 uint64 slots to store in compact form. Create a nontrivial
	y-vector by multiplying the factor by 344220 random 64-bit ints. Owing to
	the special form of M-number divisors it's exceedingly unlikely that the resulting y[]
	will have any factors in common with x[] other than q.

	26 Dec 2015, on a single core of a 2 GHz Core2 Duo: ? sec using pure-int code, ? sec using new FFT-gcd.
	*/
#if OS_BITS == 64
	if(MAX_ARRAY_DIM >= 344222)
	{
		printf("\nPerforming GCD self-test #13...\n");
		p = 22030163;	vec_len = (p>>6)+1;
		// First use this M(p) for a large-vector timing test of the mi64_div_by_scalar64 fucntion:
		memset(x,0xff,(vec_len<<3));	x[vec_len-1] = (1ull << (p&63)) - 1;
	
		// Next try a gcd of 2^p-1 (huge) with just the small factor:
		memset(y,(uint64) 0,(vec_len<<3));	y[1] = 43144178324ull; y[0] = 4788416424359163737ull;
		clock1 = clock();
		gcd_debug=0;
		fail += (mi64_gcd(x, y, vec_len, eGCD = FALSE, 0x0,0x0,0x0,0x0, 0,0x0,0x0,0x0,0x0,0) != 2);
		fail += (x[1] != 43144178324ull || x[0] != 4788416424359163737ull);
		ASSERT(HERE, fail==0,"GCD self-test #13a");
		clock2 = clock();
		tdiff = (double)(clock2 - clock1);
		printf	("GCD self-test #13a passed: GCD Time =%s\n",get_time_str(tdiff));
	
		/* Now create a full-length y-vector = [large random padding integer]*[known small factor]: */
		memset(x,0xff,(vec_len<<3));	x[vec_len-1] = (1ull << (p&63)) - 1;
		y[0] = 0;
		for(i = 0; i < vec_len-2; i++)
		{
			x[i] = y[i+1] = rng_isaac_rand();
		}
		x[vec_len-1] = x[vec_len-2] = y[vec_len-1] = 0;
		// Now y = 2^64 * x:
		x[vec_len-2] = mi64_mul_scalar(x, 4788416424359163737ull, x, vec_len-2);
		ASSERT(HERE, x[vec_len-2] != 0, "GCD self-test #13b: expected nonzero carryout of q0* y4[0:vec_len-2]!");
		y[vec_len-1] = mi64_mul_scalar_add_vec2(y, 43144178324ull, x, y, vec_len-1);
		ASSERT(HERE, y[vec_len-1] != 0, "GCD self-test #13b: expected nonzero carryout of q1* y4[0:vec_len-1]!");
		// Check that the 2-word factor actually divides y:
		x[1] = 43144178324ull; x[0] = 4788416424359163737ull;
		ASSERT(HERE, mi64_div(y, x, vec_len, 2, 0x0, 0x0), "GCD self-test #13b: divisibility test fails!");
		// ...And re-init x to hold 2^p-1:
		memset(x,0xff,(vec_len<<3));	x[vec_len-1] = (1ull << (p&63)) - 1;

	// Do the FFT-GCD:
		memcpy(u,x,(vec_len<<3));	memcpy(v,y,(vec_len<<3));
		fail  = (fft_gcd(u, v, vec_len, eGCD = FALSE, 0x0,0x0,0x0,0x0) != 1);
		fail += (x[1] != 43144178324ull || x[0] != 4788416424359163737ull);
		ASSERT(HERE, fail==0,"GCD self-test #12b");
		clock2 = clock();
		tdiff = (double)(clock2 - clock1);
		printf	("GCD self-test #13fft passed: FFT-GCD Time =%s\n",get_time_str(tdiff));

	// Pure-int GCD of same inputs:
		clock1 = clock();
		gcd_debug=0;
		fail += (mi64_gcd(x, y, vec_len, eGCD = FALSE, 0x0,0x0,0x0,0x0, 0,0x0,0x0,0x0,0x0,0) != 2);
		fail += (x[1] != 43144178324ull || x[0] != 4788416424359163737ull);
		ASSERT(HERE, fail==0,"GCD self-test #13b");
		clock2 = clock();
		tdiff = (double)(clock2 - clock1);
		printf	("GCD self-test #13int passed: GCD Time =%s\n",get_time_str(tdiff));
	}
exit(0);
#endif	// OS_BITS == 64 ?

#if 1	/****************************************************************************************/
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
		/* Init an array of (vec_len) random uint64's: */
		mi64_rand(x,vec_len);
		mi64_rand(y,vec_len);
		mi64_set_eq(u,x,vec_len);
		mi64_set_eq(v,y,vec_len);
	printf("GCD self-test #14: vec_len = %u: x^2.lo64 = %20llu, y^2.lo64 = %20llu\n",vec_len,x[0]*x[0],y[0]*y[0]);
		mi64_mul_vector(u, vec_len, u, vec_len, w, &lenZ);
		ASSERT(HERE, lenZ == 2*vec_len, "GCD self-test #14: lenZ != 2*vec_len");
		mi64_mul_vector(v, vec_len, v, vec_len, z, &lenZ);
		ASSERT(HERE, lenZ == 2*vec_len, "GCD self-test #14: lenZ != 2*vec_len");
	/* This was for testing repeated auto-square within an iterative loop:
		for(i = 0; i < vec_len; i++) {
			x[i] = y[i] = 0ull;
		}
		x[0] = 3; y[0] = 7;
	*/
		/* Init A: Convert x+I*y to packed-double form, store outputs in A: */
		fft_len = mi64_cvt_uint64_double(x,y, 0,vec_len, a);
		// For now just use 16 bits per double:
		ASSERT(HERE, fft_len == 4*vec_len, "GCD self-test #14: unexpected return value of mi64_cvt_uint64_double!");
		fft_len *= 2;	/* Multiply by 2 here because mi64_cvt_uint64_double returns *complex* FFT length
						that results from odd/even interleaving and int64-to-real conversion of the 2 input arrays
						*/
		pad_len = fft_len + ( (fft_len >> DAT_BITS) << PAD_BITS );	/* Padded-array length */
		ASSERT(HERE, !ARRAYS_OVERLAP(a,2*pad_len,b,2*pad_len), "GCD self-test #14: a/b-arrays overlap!");
		/* Zero-pad the input vector in preparation for FFT-based multiply: */
		memset(a+pad_len, 0, (pad_len<<3));
	/*
		printf("FFTmul inputs:\n");
		for(i = 0; i < 2*fft_len; i+=2)
		{	*** need padded-index here! ***
			printf("[%3u]: %20.5f %20.5f\n",i>>1,a[i],a[i+1]);
		}
	*/
		/* Call the init-FFT routines, using the (as yet uninited) b-array for the needed scratch space: */
		pairFFT_mul(  b,a,0, 2*fft_len, TRUE , FALSE);	/* INIT_ARRAYS = TRUE , FORWARD_FFT_ONLY = FALSE */
		/* Save a copy of the input vector: */
		memcpy(b,a,(pad_len<<4));	// double the length of memcpy, thus no 2nd memset(0) step needed to 0-pad b.
		/* Forward FFT of A-vector: */
		pairFFT_mul(  a,0,0, 2*fft_len, FALSE, TRUE );	/* INIT_ARRAYS = FALSE, FORWARD_FFT_ONLY = TRUE  */
	/*
		printf("FFTmul outputs:\n");
		for(i = 0; i < 2*fft_len; i+=2)
		{	*** need padded-index here! ***
			printf("[%3u]: %20.5f %20.5f\n",i>>1,a[i],a[i+1]);
		}
	*/
		/* Now do the FFT-based squaring: */
		pairFFT_mul(  b,a,0, 2*fft_len, FALSE, FALSE);	/* INIT_ARRAYS = FALSE, FORWARD_FFT_ONLY = FALSE */

		/* Convert result back to uint64[] form and compare 2 outputs to all-integer squarings.

			NOTE: u,v,x,y each only alloc'ed to MAX_ARRAY_DIM elts, thus not large enough to store
		a full double-width product. But since we only need two pointers for the back-converted FFT outputs,
		we use the 2*MAX_ARRAY_DIM-long alloc for u+v (accessed via u-ptr) for Re part of the FFT output,
		and similarly use x+y (accessed via x-ptr) for the Im part of the FFT output:
		*/
		ASSERT(HERE, mi64_cvt_double_uint64(b, fft_len, u,x) == 0x0, "GCD self-test #14: nonzero carryout from mi64_cvt_double_uint64!");

	//	ASSERT(HERE, mi64_cmp_eq(u,w, 2*vec_len), "GCD self-test #14: Re part of the FFT output mismatches reference!");
		for(i = 0; i < 2*vec_len; i++) {
			ASSERT(HERE, u[i] == w[i], "GCD self-test #14: Re part of the FFT output mismatches reference!");
		}

	//	ASSERT(HERE, mi64_cmp_eq(x,z, 2*vec_len), "GCD self-test #14: Im part of the FFT output mismatches reference!");
		for(i = 0; i < 2*vec_len; i++) {
			ASSERT(HERE, x[i] == z[i], "GCD self-test #14: Im part of the FFT output mismatches reference!");
		}

		vec_len *= 2;

		if(vec_len > MAX_ARRAY_DIM)
			break;
	}

	/* Now turn the pairwise-linear-combination version of the genFFT_mul algo loose on some actual GCD vector pairs.
	Here is the basic idea:

		Using the old quadratic 64-bit-int version of the code (as captured in the mi64_gcd function), we compute
	GCD(u,v) by reducing both vectors u,v ~64 bits per iteration. To do this we use the leading 128 bits of u,v to
	construct 2 pairs of 64-bit multipliers [a,b] and [c,d] such that u' = a*u - b*v and v' = c*u - d*v are just
	64 bits each. Each of the 4 subproducts a*u, b*v, c*u and d*v is 192 bits, thus the multipliers are chosen such
	that the leading 128 bits of the pairs [a*u, b*v] and [c*u, d*v] match and cancel upon subtraction.

	Example: base = 2^64, 128-bit inputs:
		u = base*15243067522374418434 +  8684254506156651683;
		v = base* 2503926509447010181 + 14921004030550884877;

		a =    29965985744110959; b =   182422903527176837;
		c =   900835843694674153; d =  5483987465369195218;

a*u = base^2* 24761743440837882 + base*18356059860786658291 + 4126082876808790445
b*v = base^2* 24761743440837882 + base*18356059860786658288 + 8228777268749009601
a*u-b*v = 51237537829188435692, 65 bits.

c*u = base^2* 744386192877435340 + base*12608466216960890642 + 10700841904888342619
d*v = base^2* 744386192877435340 + base*12608466216960890642 + 11791441100596549802
c*u-d*v = -1090599195708207183, 64-bit.

But note that if we use only the leading 64-bit word of u,v in the 2 linear combinations, we get

a*u.hi64 = base^2*  24761743440837882 + base*18341952642423900894
b*v.hi64 = base^2*  24761743440837882 + base*18208503583759260185

c*u.hi64 = base^2* 744386192877435340 + base*12184375781074026962
d*v.hi64 = base^2* 744386192877435340 + base* 8172637392790005018 ,

and only the leading ~64 bits of the resulting subproducts cancel upon subtraction, leaving 64 bits
plus the unprocessed low words. This would indicate that we must apply the N-bit [a,b,c,d] multipliers
to the same leading 2N bits of the inputs u,v which were used to compute [a,b,c,d], but note what we
get by applying the same multipliers to low halves of u,v:

a*u.lo64 =                              base*   14107218362757397 +  4126082876808790445
b*v.lo64 =                              base*  147556277027398103 +  8228777268749009601

c*u.lo64 =                              base*  424090435886863680 + 10700841904888342619
d*v.lo64 =                              base* 4435828824170885624 + 11791441100596549802 .

Now add the base* coefficients in the hi64-products to those in the lo64-products:

18341952642423900894 +   14107218362757397 = 18356059860786658291	<*** low 2 bits
18208503583759260185 +  147556277027398103 = 18356059860786658288	<*** mismatch!

12184375781074026962 +  424090435886863680 = 12608466216960890642
 8172637392790005018 + 4435828824170885624 = 12608466216960890642 ,

and we see that except for the low 2 bits of the first case, the results match. (And if we had simply
included a few 'guard bits' into our computation, we could guarantee with a high probability that all
the results would match perfectly.) Thus, with appropriate bitcounts in our [a,b,c,d] multipliers and
the chunks of u,v to which they are applied, we *know* that the high 2N bits of the full length-3N
[a,b,c,d]*[lead 2N bits of u,v] linear combos will match, thus we only need to explicitly compute the
low N bits, via a single length-2N FFT-mul to get the length-2N linear combos [a*u.loN - b*v.loN] and
[c*u.loN - d*v.loN] and discarding of the high N-bit halves of the results.

x=d*v; hi = x>>64
lo = x-(hi<<64)
******* See if can actually get low-halves-only-FFT-mul scheme working, then modify appropriately: ********
For large u,v, inputs in order to take advantage of fast FFT-mul, we run mi64_gcd in eGCD mode to compute
(say) 8192-bit [a,b,c,d] multipliers at a time, using the leading 16384 bits of u,v. We then use FFT-mul
to compute our linear combos, in 2 stages:

[1] FFT-mul to compute [a*u.lo8192 - b*v.lo8192] and [c*u.lo8192 - d*v.lo8192], each 16384 bits,
where lo8192 refers to the lower halves of the leading 16384-bits of u,v used to compute the multipliers.

[2] FFT-mul to compute [a*u.hi8192 - b*v.hi8192] and [c*u.hi8192 - d*v.hi8192], each 16384 bits,
where hi8192 refer to the upper halves of the leading 16384-bits of u,v used to compute the multipliers.
We then add the upper 8192 bits of the 2 lincombos computed in [1] to the low 8192 bits of the high-part
lincombos and propagate carries, at which point we should find all but perhaps the lowest few of the high
16384 bits of the full length-3*N = 3*8192 = 24576-bit lincombos cancelling to 0.

Cost of the FFT-muls in [1],[2]:

[A] Cost of generating the 8192 mults [a,b,c,d] using pure-int eGCD and 64-bit-at-a-time reduction: O((N/64)^2);
[B] Length-2*N complex FFT to compute length-2*N fFFT(u.lo8192, v.lo8192);
[C] Length-2*N complex FFT to compute length-2*N fFFT(u.hi8192, v.hi8192);
[D] Length-2*N complex FFT to compute length-2*N fFFT(a,b);
[E] Length-2*N complex FFT to compute length-2*N fFFT(c,d);
[F] O(N) dyad-mul and Length-2*N complex iFFT to compute length-2*N lincombos in [1];
[G] O(N) dyad-mul and Length-2*N complex iFFT to compute length-2*N lincombos in [2];
[H} O(N) add/carry step to combine the 2 length-2*N vector pairs (via overlap of middle N bits).

Thus asymptotically the cost is that of 6 length-2*N complex FFTs.
	*/
// FFTMUL_THRESHOLD_BITS = #leading bits of u,v to use to generate FFT-multipliers.
// This corresponds to 2N in above commentary, thus said multipliers are half as long.

#endif	/****************************************************************************************/

	return fail;
}

#undef YES_ASM
