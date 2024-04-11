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

/* To build inside Mlucas/src on Mac:
clang -c -O3 -DINCLUDE_GMP=0 -DINCLUDE_HWLOC=0 get_cpuid.c get_fp_rnd_const.c imul_macro.c mi64.c qfloat.c rng_isaac.c twopmodq.c types.c util.c
clang -c -O3 -g3 -ggdb -DINCLUDE_GMP=0 -DINCLUDE_HWLOC=0 -DPM1_STANDALONE [-DPM1_DEBUG] -O3 pm1.c
clang -o pm1 *.o -Xlinker --no-demangle
Then to run, e.g.
./pm1 -bigstep [210|330|420|660|840] -b1 5000000 -b2 150000000 -m 1
*/
#include "Mlucas.h"
#define STR_MAX_LEN 1024

#ifdef PM1_STANDALONE
	#warning Building pm1.c in PM1_STANDALONE mode.
	char STATFILE[] = "pm1_debug.txt";
	uint32 pm1_standlone = 1;
	// Externs normally def'd in Mlucas.c:
	uint64 RES_SHIFT = 0ull;
	uint32 PRP_BASE = 0;
	uint64 *PM1_S1_PRODUCT = 0x0, PM1_S1_PROD_RES64 = 0ull;	// Vector to hold Stage 1 prime-powers product product in
					// most-significant-bit-deleted-and-result-bit-reversed form, and (mod 2^64) checksum on same.
	uint32 PM1_S1_PROD_B1 = 0, PM1_S1_PROD_BITS = 0;	// Stage 1 bound to which the current value of PM1_S1_PRODUCT corresponds, and its #bits
	uint32 PM1_S2_NBUF = 0;	// # of floating-double residue-length memblocks available for Stage 2
	uint32 B1 = 0;
	uint64 B2 = 0ull, B2_start = 0ull;
	char cbuf[STR_MAX_LEN*2];
	uint32 SYSTEM_RAM, MAX_RAM_USE;	// Total usable main memory size, and max. amount of that to use per instance, in MB
	double MME;
#else
	uint32 pm1_standlone = 0;
#endif

#undef RTIME
#undef CTIME

#ifndef PM1_STANDALONE
  #ifdef MULTITHREAD
	#define RTIME	// In multithreaded mode, need to use real (wall-clock) time

	#include "threadpool.h"
	// Simple threadpool-struct for parallel|SIMD vector-double subtract c[] = a[] - b[]:
	struct pm1_thread_data_t{
		int tid;
		int*retval;
		double*arr0;			// Output array0 = c[], including address-offset into a given thread-processed chunk
		double*arr1;			// Input  array1 = a[], ditto
		double*arr2;			// Input  array2 = b[], ditto
		int n;					// Chunksize
	};
	// Stick protos fo these SIMD utility functions used by the stage 2 loop here:
	void vec_double_sub(struct threadpool *tpool, struct pm1_thread_data_t *tdat, double c[]);
	void*vec_double_sub_loop(void*targ);
  #else
	#define CTIME	// In single-thread mode, prefer cycle-based time because of its finer granularity
	void vec_double_sub(double a[], double b[], double c[], uint32 n);
	void vec_double_sub_loop(double a[], double b[], double c[], uint32 n);
  #endif
#endif

/*************** Bytewise utility routines needed by prime-pairing algorithm ***************/
void bytevec_bitstr(uint8*x, int nbytes, char*ostr)
{
	int i;
	for(i = 0; i < nbytes; i++) {
		// High byte ==> leftmost 8 chars of output string, thus the (nbytes-1-i)
		byte_bitstr(x[nbytes-1-i], ostr + (i<<3));
	}
}

// Copy: xout = xin:
void bytevec_set_eq(uint8*xout, const uint8*xin, int nbytes) {
	int i;
	for(i = 0; i < nbytes; ++i) {
		xout[i] = xin[i];
	}
}

// Bit-reversal: brev8[] is lookup-table of bit-reversed bytes def'd in mi64.h. Allows in-place:
void bytevec_brev(uint8*xout, uint8*xin, int nbytes) {
	int i,j;
	uint8 tmp1,tmp2;
	// To support in-place bit-reversal (xin == xout), exchange bit-reversed pairs of bytes starting
	// with bookending byte-pair and working inward toward middle; 'middle byte' for odd #bytes means
	// redundancy in computing xi,xj, but eschew special post-loop code just to save a few ops in such cases:
	for(i = 0; i < (nbytes>>1); ++i) {
		j = nbytes-i-1;	// i,j = indices of the 2 bytes containing the bits to be swapped on this loop pass
		// This works even if i == j, i.e. final loop-exec when nbytes odd:
		tmp1 = brev8[xin[i]];	tmp2 = brev8[xin[j]];
		xout[i] =       tmp2;	xout[j] =       tmp1;
	}
}

// xout = xin1 [&,|] xin2:
void bytevec_and(uint8*xout, uint8*xin1, uint8*xin2, int nbytes) {
	int i;	for(i = 0; i < nbytes; ++i) { xout[i] = xin1[i] & xin2[i]; }
}
void bytevec_or(uint8*xout, uint8*xin1, uint8*xin2, int nbytes) {
	int i;	for(i = 0; i < nbytes; ++i) { xout[i] = xin1[i] | xin2[i]; }
}

// Test selected bit - we allow bidirectional relative offsets:
int bytevec_test_bit(uint8*x, int bit) {
	int i = bit>>3, j = bit - (i<<3);	// i = byte containing bit-to-be-tested, j = bit within said byte
	uint8 mask = 1 << j;
	return (x[i] & mask) != 0;
}

void bytevec_clear(uint8*x, int nbytes) {
	int i;
	for(i = 0; i < nbytes; ++i) {
		x[i] = 0;
	}
}

int bytevec_iszero(uint8*x, int nbytes) {
	int i;
	for(i = 0; i < nbytes; ++i) {
		if(x[i]) return 0;
	}
	return 1;
}

// Set specified bit in a byte-vec:
void bytevec_bset(uint8*x, int bit) {
	int i = bit>>3, j = bit - (i<<3);
	uint8 mask = 1 << j;
	x[i] |= mask;
}

// Clear specified bit in a byte-vec:
void bytevec_bclr(uint8*x, int bit) {
	int i = bit>>3, j = bit - (i<<3);
	uint8 mask = ~(1 << j);
	x[i] &= mask;
}

// Same-index bit-clear in 2 byte-vecs:
void bytevec_bclr_pair(uint8*x, uint8*y, int bit) {
	int i = bit>>3, j = bit - (i<<3);
	uint8 mask = ~(1 << j);
	x[i] &= mask;	y[i] &= mask;
}

/* Set default p-1 stage bounds for M(p), in the context of length-n-doubles FFT-modmul.
This needs to be enhanced to pick B1 and B2 based on the following nontrivial optimization problem:
	For given (TF done to 2^tf_bits with no factor found)
	and (p-1 attempted to stage bounds b1_prev and b2_prev with no factor found),
we seek to minimize (S*(1-P) - W), where
	S = (#primality tests saved if p-1 factor found at computed bounds),
	P = (probability of p-1 factor being found at computed bounds),
	W = (computational effort of p-1 to computed bounds, relative to a unit based on the computational effort of a primality test) .
The idea is that
	S = total work of doing the primality test(s) in the absence of further factoring being attempted;
	S*P = (probable, i.e. on-average) primality-test work saved via further factoring being attempted to the bounds-to-be-determined;
	W = work expended via further factoring being attempted to the bounds-to-be-determined.
In other words, (S*(1-P) - W) is an average-total-work function, which we attempt to minimize.

The resulting bounds are set in form of the globals B1 and B2.
Return value != 0 indicates success; return 0 indicates an error of some kind, in which case B1 and B2 are set = 0.
On successful execution, B2 > B1 > 0 and PM1_S2_NBUF will have been set based on available RAM and FFT-length-in-doubles, n.
*/
uint32 pm1_set_bounds(const uint64 p, const uint32 n, const uint32 tf_bits, const double tests_saved)
{
	const double inv100k = 1./100000, inv1m = 1./1000000;
/* NB: Since GIMPS testing will all soon be 1-shot PRP-with-proof-of-correctness, no need to make use of tests_saved. */
	// Get PM1_S2_NBUF based on available RAM; n doubles need n/2^17 MB; use KB and double-math for the intermediates here.
	// First compute FP-form #bufs if < 5 set PM1_S2_NBUF = 0, to prevent (dtmp - 5) from giving nonsense after cast-to-uint32:
	uint64 i64;
	double dtmp = (SYSTEM_RAM*(double)MAX_RAM_USE*0.01)*1024./(n>>7);
	sprintf(cbuf,"pm1_set_bounds: Stage 2 needs at least 24+5 buffers ... each buffer needs %u MB; avail-RAM allows %u such.\n",(n>>17),(uint32)dtmp);
	mlucas_fprint(cbuf,pm1_standlone+1);
	if(PM1_S2_NBUF && PM1_S2_NBUF > ((uint32)dtmp - 5)) {
		sprintf(cbuf,"WARNING: User-specified Stage 2 buffer count %u exceeds avail-RAM ... each buffer needs %u MB; avail-RAM allows %u such.\n",PM1_S2_NBUF,(n>>17),(uint32)dtmp);
		mlucas_fprint(cbuf,pm1_standlone+1);
	} else if(!PM1_S2_NBUF) {	// Respect any user-set value here, so long as it's >= minimum-buffer-count
		if(dtmp < 5.) {
			PM1_S2_NBUF = 0;
		} else {	// Assume 5 n-double chunks already alloc'ed for Stage 1 residue and other stuff:
			PM1_S2_NBUF = (uint32)dtmp - 5;
		}
	}
	// Force B1 >= 10^4 to avoid possible large-buffer-count underflow of qlo in stage 2.
	// Conservatively use (#bits in Stage 1 prime-powers product ~= 1.5*B1), must fit into a uint32, thus B1_max = 2^33/3 = 2863311530:
	i64 = p>>7;
	ASSERT(i64 <= 2863311530ull, "Stage 1 prime-powers product must fit into a uint32; default B1 for your exponent is too large!");
	B1 = MAX((uint32)i64,10000);	// #bits in Stage 1 prime-powers product ~= 1.4*B1, so e.g. B1 = p/128 gives a ~= 1.1*p/100 bits
	B1 = (B1 + 99999)*inv100k;	B1 *= 100000;	ASSERT(B1 >= 100000, "B1 unacceptably small!");	// Round up to nearest 100k:
	if(PM1_S2_NBUF < 24) {
		sprintf(cbuf,"pm1_set_bounds: Insufficient free memory for Stage 2 ... will run only Stage 1.\n");
		mlucas_fprint(cbuf,pm1_standlone+1);
		B2_start = B2 = (uint64)0;
		// If no stage 2 being done, run slightly deeper stage 1:
		B1 = (B1 * 5)>>2;
		B1 = (B1 + 99999)*inv100k;	B1 *= 100000;
	} else {
		/* Relative cost of Stage 2 for bigstep D = 210 and the minimal-memory scheme, stage2_mem_multiple := M = 1:
		Each D-sized bigstep yields ~ 13 modmul (roughly 1/5th in form of prime-pairs, ~4/5 as/2+)
		and needs 2 added modmuls for end-of-loop powering-up of the base multiplier, thus ~15 modmuls per D-step.
		Thus for given B1 and B2, min-memory Stage 2 needs ~(B2-B1)/14 modmuls, vs ~1.4*B1 modmuls for Stage 1.
		To roughly equalize the work-per-stage we have (B2-B1)/14 ~= 1.4*B1, or roughly B2 ~= 30*B1.

		For M --> oo, the minimum achievable cost ~= 0.5 * this. Cf. my pm1_compare.png plot for a graph of the data.
		Playing with various curve fits using Mac Grapher (whose 'log' = log10, hence the log10-ness of our fit) gives
		ratio ~= [1 - 0.93 log10(log10(#buffers/2.7))] as a decent approximation; multiply our min-memory B2 by 1/ratio:

		BC code - In a *nix terminal, type 'bc -l' and then paste the code into the resulting bc session::
			define log10(x) {
				return l(x)/l(10);
			}
			define b2_b1_ratio(nbuf) {
				x = 1 - 0.93*log10(log10(nbuf/2.7));
				r = 30/x;
				return r;
			}

		Flipping the curve-fit formula around, given some desired b2/b1 value, we want to find the minimal buffer count needed to achieve this. Given r = b2/b1, we have
			r = 30/(1 - 0.93*log10(log10(nbuf/2.7))), which we wish to solve for nbuf.
		->	30 = r*(1 - 0.93*log10(log10(nbuf/2.7)))
		->	(30-r) = -0.93*r*log10(log10(nbuf/2.7)))
		->	log10(log10(nbuf/2.7))) = (r-30)/(0.93*r); let x = (r-30)/(0.93*r).
		->	log10(nbuf/2.7)) = 10^x
		->	nbuf = 2.7*10^(10^x)

		Implementing the atter expression in bc-code, BC alas doesn't support exponentiation with fractional powers, so we roll our own - to compute x = a^b, we first compute log(x) = b*log(a) then exponentiate:

			define power(a,b) {
				return e(b*l(a));
			}
			define pow10(x) {
				return power(10,x);
			}
			define nbuf(b2_b1_ratio) {
				x = (b2_b1_ratio-30)/(0.93*b2_b1_ratio);
				nbuf = 2.7*pow10(pow10(x));
				return nbuf;
			}
		*/
		uint32 bigstep = 0, stage2_mem_multiple = 0, psmall = 0;	// For initial-bounds setting, no S2 relocation-prime set yet
		pm1_bigstep_size(&PM1_S2_NBUF, &bigstep, &stage2_mem_multiple,psmall);
		if(bigstep != 210 && bigstep != 330 && bigstep != 420 && bigstep != 660 && bigstep != 840) {
			sprintf(cbuf,"%u is unsupported value of bigstep!",bigstep);
			mlucas_fprint(cbuf,pm1_standlone+1);	ASSERT(0,cbuf);
		}
		double f2 = 30.0/(1 - 0.93*log10(log10(0.37037037037037037037*(double)PM1_S2_NBUF)));
		B2 = (uint64)(f2*(double)B1);
		// Round to nearest 1m:
		B2 = (B2 + 999999)*inv1m;	B2 *= 1000000;	ASSERT(B2 >= 1000000, "B2 unacceptably small!");
	}
	pm1_check_bounds();	// This sanity-checks the bounds and sets B2_start = B1 if unset.
	sprintf(cbuf,"Setting default p-1 stage bounds b1 = %u, b2_start = %" PRIu64 ", b2 = %" PRIu64 ".\n",B1,B2_start,B2);
	mlucas_fprint(cbuf,pm1_standlone+1);
	return 1;
}

/* Check the p-1 stage bounds stored in the globals B1 and B2. Return value != 0 indicates success;
return 0 indicates an error of some kind; in this case B1 and B2 are left as-is to aid in debugging.
*/
uint32 pm1_check_bounds() {
	while(1) {
		if(B1 == 0) { sprintf(cbuf,"P-1 requires at least a nonzero Stage 1 bound to be specified via the -b1 flag.\n"); break; }
		// Force B1 >= 10^4 to avoid possible large-buffer-count underflow of qlo in stage 2:
	#ifndef PM1_DEBUG	// Allow any B1 in debug-build mode
		if(B1 < 10000) { sprintf(cbuf,"The minimum P-1 Stage 1 bound = 10000; resetting to that.\n"); mlucas_fprint(cbuf,pm1_standlone+1); B1 = 10000; }
	#endif
		if(B2_start) {
			if(B2_start > B2) { sprintf(cbuf,"P-1 Stage 2 starting bound [= %" PRIu64 "] must be less than or equal to Stage 2 bound [= %" PRIu64 "].\n",B2_start,B2); break; }
			if(B1 > B2) { sprintf(cbuf,"P-1 Stage 2 bound [= %" PRIu64 "] must be greater than or equal to that of Stage 1 [= %u].\n",B2,B1); break; }
		} else if(B2) {	// Stage 2 takes off where Stage 1 left off
			if(B1 > B2) { sprintf(cbuf,"P-1 Stage 2 bound [= %" PRIu64 "] set nonzero but < Stage 1 bound [= %u] ... no Stage 2 will be run.\n",B2,B1); }
			B2_start = B1;
		} else {	// No Stage 2 - Can set both of these to 0 or B1 in this case
			B2_start = B2 = (uint64)0;
		}
		return 1;	// B1 and B2 legal.
	}
	mlucas_fprint(cbuf,pm1_standlone+1);	ASSERT(0,cbuf);
	return 0;	// Bzzt!
}

/* Re. mem-alloc for Stage 1 prime-powers product:
For a given n, there are Pi(n) ~ n/ln(n) primes, but this undercounts them. We seek a similarly-simple estimate
which slightly overcounts the primes up to 2^32. We approximate Pi(n) ~= n/(ln(n)-A), where A is an adjustable constant.
Using that Pi(10^9) = 50847534 gives A = ln(n) - n/Pi(n) = 1.0566287... for this choice of n.
For n = 10^6 this choice of A gives 78376, just less than the actual Pi(n) = 78498. So take A = 1.1 for safety's sake.
No point trying to be finer, since we further estimate the average bitness of primes < n as lg(n)-1. Thus our estimate
of #bits of the product of primes <= n is B(n) ~= (lg(n)-1)*n/(ln(n)-A). This does not even try to account for the fact
that our product of Stage 1 prime powers uses powers > 1 for the primes < sqrt(B1), but the contribution of such is
relatively negligible, as long as we choose our tuning parameter A loosely enough. Here are some computed-vs-estimated:
			#bits in S1 prime-powers product:
	b1		computed	estimated B(n) ~= (lg(n)-1)*n/(ln(n)-A), A = 1.1
	----	--------	--------
	10^3	1438		1526
	10^4	14460		15027
	10^5	144344		148946
	10^6	1442099		1480991
	10^7	14424867	14751202	Time on Core2: 140s
						sans mi64_mul_scalar call:   5s	Consider using FFT-mul to accumulate large sub-chunks.
						Quick-optimization to minimize number of mi64_mul_scalar calls, while keeping scalar-mult < 2^64
						makes little difference here, because e.g. for b1 = 10^7 the accumulation is dominated by primes
						22-23 bits in length, i.e. our first-cut strategy of combining 2 consecutive prime-power subproducts
						into each scalar-multiplier resulted in nearly the same average scalar-mult size, the savings is all
						for the small primes, which are relatively few in number. But the new routine still looks cleaner.
*/
/* Function allocs mem via global PM1_S1_PRODUCT array, wraps call to pm1_s1_ppow_prod(), etc.
Inputs: A base-2 exponent. Modulus can be a Mersenne M(p) or Fermat F(m); in the latter case our input p contains 2^m.
Assumes: Stage 1 bound has been stored in global B1 prior to call.
Effects: Stores computed stage 1 prime-powers product in PM1_S1_PRODUCT and #bits of same in PM1_S1_PROD_BITS; sets PRP_BASE = 3.
Returns: #limbs alloc'ed.

NOTE: If PM1_S1_PROD_B1 == B1, current Stage 1 prime-powers-product allocation is the correct one for this B1,
but in order to re-use it would need to divide out the Mersenne|Fermat-exponent-specific seed - which means another
global would be needed to store that - and remultiply by the appropriate one for the current p, not worth the effort.
*/
uint32 compute_pm1_s1_product(const uint64 p) {
	const double A = 1.1;
	ASSERT(B1 > 0, "Call to compute_pm1_s1_product needs Stage 1 bound global B1 to be set!");
	double ln = log(B1), lg = ln*ILG2;
	uint32 i,len = 0,nmul,nbits,ebits = (uint32)((lg-A)*B1/(ln-A));
	uint64 iseed,maxmult;
	char savefile[STR_MAX_LEN];

	// Compute Stage 1 prime-powers product, starting with alloc of needed memory:
	uint32 s1p_alloc = ((ebits + 63)>>6) + 1;	// Add 1 to account for seeding-by-binary-exponent described below
	PM1_S1_PRODUCT = ALLOC_UINT64(PM1_S1_PRODUCT, s1p_alloc);
	if(!PM1_S1_PRODUCT ){
		sprintf(cbuf, "ERROR: unable to allocate array PM1_S1_PRODUCT with %u linbs in main.\n",s1p_alloc);
		mlucas_fprint(cbuf,pm1_standlone+1);	ASSERT(0,cbuf);
	}

	// (E.g. on restart) First see if a savefile holding the precomputed/bit-reversed product for this p and B1 exists:
  #ifndef PM1_STANDALONE
	strcpy(savefile, RESTARTFILE);
	savefile[0] = ((MODULUS_TYPE == MODULUS_TYPE_MERSENNE) ? 'p' : 'f');
	strcat(savefile, ".s1_prod");
	if((len = read_pm1_s1_prod(savefile, p, &PM1_S1_PROD_BITS, PM1_S1_PRODUCT, &PM1_S1_PROD_RES64)) != 0) {
		sprintf(cbuf, "INFO: Successfully read precomputed/bit-reversed Stage 1 prime-powers product savefile for this modulus and B1 = %u.\n",B1);
		mlucas_fprint(cbuf,pm1_standlone+1);
	} else {	// Compute product from scratch:
  #endif
		// For M(p) want to seed the S1 prime-powers product with 2*p; for F(m) we want seed = 2^(m+2). Since in the latter
		// case our input p contains 2^m, can handle both cases via iseed = 4*p, giving an extra *2 in the Mersenne case:
		iseed = p<<2;	ASSERT((iseed>>2) == p,"Binary exponent overflows (uint64)4*p in compute_pm1_s1_product!");
		len = pm1_s1_ppow_prod(iseed, B1, PM1_S1_PRODUCT, &nmul, &maxmult);	PM1_S1_PROD_B1 = B1;
		nbits = (len<<6)-mi64_leadz(PM1_S1_PRODUCT,len);
		if(len > s1p_alloc) {
			sprintf(cbuf,"Size of S1 prime-powers product exceeds alloc of PM1_S1_PRODUCT[]!");
			mlucas_fprint(cbuf,pm1_standlone+1);	ASSERT(0,cbuf);
		}
	/*
		fprintf(stderr,"Product of Stage 1 prime powers used %u mi64_mul_scalar() calls; max-multiplier %u bits\n",nmul, 64-leadz64(maxmult));
		fprintf(stderr,"Limbs of PM1_S1_PRODUCT, low to high:\n");
		for(i = 0; i < len; i+=8) {
			fprintf(stderr,"%" PRIx64 ",%" PRIx64 ",%" PRIx64 ",%" PRIx64 ",%" PRIx64 ",%" PRIx64 ",%" PRIx64 ",%" PRIx64 "\n",PM1_S1_PRODUCT[i],PM1_S1_PRODUCT[i+1],PM1_S1_PRODUCT[i+2],PM1_S1_PRODUCT[i+3],PM1_S1_PRODUCT[i+4],PM1_S1_PRODUCT[i+5],PM1_S1_PRODUCT[i+6],PM1_S1_PRODUCT[i+7]);
		}
		exit(0);
	*/
	//	fprintf(stderr,"PM1_S1_PRODUCT limbs[%u,%u,...,1,0] = %016" PRIX64 ",%016" PRIX64 ",...,%016" PRIX64 ",%016" PRIX64 "\n",len-1,len-2,PM1_S1_PRODUCT[len-1],PM1_S1_PRODUCT[len-2],PM1_S1_PRODUCT[1],PM1_S1_PRODUCT[0]);
		// Ignore the #iters != 0 user needed to set to invoke selfTest mode, replace with nbits in S1 prime-powers product:
		PM1_S1_PROD_BITS = nbits-1;	// Leftmost bit accounted for by setting initial seed in the LR-modular binary powering
		// Bit-reverse s1 product, leaving leftmost 1-bit off. REMEMBER, this puts the 0-bits corresponding to the
		// power of 2 in the product leftmost, thus PM1_S1_PRODUCT may end up with < len nonzero words, we have some
		// number of 0-but-significant leftmost bits in the bit-reversed prime-powers product. Don't need to explicitly zero
		// the (64 - PM1_S1_PROD_BITS%64) upper bits of high word since mi64_brev's >>= pad_bits does it for us, but better safe than sorry:
		uint32 pad_bits = 64 - (PM1_S1_PROD_BITS & 63);
		if(pad_bits != 0)
			PM1_S1_PRODUCT[len-1] &= (-1ull >> pad_bits);
		else	// If leading-bit-ignored product a multiple of 64 bits long, our working length is 1 limb shorter than otherwise:
			--len;
		mi64_brev(PM1_S1_PRODUCT, PM1_S1_PROD_BITS);
		// Compute a simple (mod 2^64) checksum, which can be rechecked at start of every savefile-update interval for data integrity:
		PM1_S1_PROD_RES64 = 0ull;
		for(i = 0; i < len; i++) { PM1_S1_PROD_RES64 += PM1_S1_PRODUCT[i]; }
  #ifndef PM1_STANDALONE
		// Write result to savefile:
		if(!write_pm1_s1_prod(savefile, p, PM1_S1_PROD_BITS, PM1_S1_PRODUCT, PM1_S1_PROD_RES64)) {
			snprintf(cbuf,STR_MAX_LEN*2,"WARN: Unable to write precomputed/bit-reversed Stage 1 prime-powers product to savefile %s.\n",savefile);
			mlucas_fprint(cbuf,pm1_standlone+1);
		}
	} 	// endif(read_pm1_s1_prod)
  #endif
	sprintf(cbuf,"Product of Stage 1 prime powers with b1 = %u is %u bits (%u limbs), vs estimated %u. Setting PRP_BASE = 3.\n",B1,PM1_S1_PROD_BITS+1,len,ebits);
	mlucas_fprint(cbuf,pm1_standlone+1);
	PRP_BASE = 3;
	sprintf(cbuf,"BRed (PM1_S1_PRODUCT sans leading bit) has %u limbs, Res64 = %" PRIu64 "\n",len,PM1_S1_PROD_RES64);
	mlucas_fprint(cbuf,pm1_standlone+0);
	return len;	// return actual #limbs of product, not initial overestimate
}

// Compute product of Stage 1 prime powers and store in a uint64[] accumulator.
// Pointer-args nmul and maxmult return #mi64_mul_scalar calls and max value of the scalar multiplier for same:
uint32 pm1_s1_ppow_prod(const uint64 iseed, const uint32 b1, uint64 accum[], uint32 *nmul, uint64 *maxmult) {
	uint32 p = 2,i,j,len,maxbits = 64-leadz64(b1);
	uint32 loop = 64/maxbits;	// Number of prime-powers we can accumulate inside inner loop while remaining < 2^64
	uint64 tmp,prod,mult,cy = 0ull;
	ASSERT(accum != 0x0, "Null accum[] pointer in s1_ppow_prod()");
	ASSERT(accum != 0x0, "Zero initial seed in s1_ppow_prod()");
	accum[0] = iseed; len = 1; *nmul = 0; *maxmult = 0ull;
// Debug-only - allows testing of S1 on known-factor case without actually running S2:
#if 0
  if(MODULUS_TYPE == MODULUS_TYPE_FERMAT) {
	#warning DEBUG-ONLY! Including the known outlier-prime in q-1 for the testcase F31, q = 46931635677864055013377:
	mult = 140091319777ull;
	cy = mi64_mul_scalar(accum, mult, accum, len);	++*nmul;
	accum[len] = cy; len += (cy != 0ull);
	fprintf(stderr,"Pre-loop accumulator = %" PRIu64 " + 2^64*%" PRIu64,accum[0],accum[1]);
  }
#endif
//	fprintf(stderr,"Stage 1 exponent = %" PRIu64 ".",accum[0]);
	while(p < b1) {
		mult = 1ull;
		for(i = 0; i < loop; i++) {
			prod = p; tmp = prod*p; j = 1;
			// if() uses prod*p here so as to include only the 1st power of the primes >= sqrt(b1):
			while(tmp <= b1) {
				prod = tmp; tmp *= p; j++;
			}
			mult *= prod;
		/*	if(j > 1)
				fprintf(stderr,"%u^%u.",p,j);
			else
				fprintf(stderr,"%u.",p);	*/
			p = next_prime(p,1);
		}
		*maxmult = MAX(mult,*maxmult);
		cy = mi64_mul_scalar(accum, mult, accum, len);	++*nmul;
		accum[len] = cy; len += (cy != 0ull);
	}
//	fprintf(stderr,"\n");
	return len;
}

// Returns 1 on successful read, 0 otherwise:
int read_pm1_s1_prod(const char*fname, uint64 p, uint32*nbits, uint64 arr[], uint64*sum64)
{
	const char func[] = "read_pm1_s1_prod";
	int retval = 0;
	uint8 c;
	uint32 i,j,b1 = 0,nbytes,nlimbs;
	uint64 itmp64 = 0ull,isum64 = 0ull;
	ASSERT(arr != 0x0, "Null arr pointer!");
	ASSERT(strlen(fname) != 0, "Empty filename!");
  #ifdef PM1_STANDALONE
	FILE*fptr = 0x0;
	goto PM1_S1P_READ_RETURN;
  #else
	FILE*fptr = mlucas_fopen(fname, "rb");
	if(!fptr) {
		fprintf(stderr,"INFO: precomputed p-1 stage 1 primes-product file %s not found...computing from scratch.\n",fname);
		goto PM1_S1P_READ_RETURN;
	} else
		fprintf(stderr,"INFO: precomputed p-1 stage 1 primes-product file %s found...reading...\n",fname);

	if(!test_types_compatible((i = fgetc(fptr)), TEST_TYPE)) {
		sprintf(cbuf, "ERROR: %s: TEST_TYPE[%u] mismatches one[%u] read from savefile!\n",func,TEST_TYPE,i);
		goto PM1_S1P_READ_RETURN;
	}
	if((i = fgetc(fptr)) != MODULUS_TYPE) {
		sprintf(cbuf, "ERROR: %s: MODULUS_TYPE[%u] mismatches one[%u] read from savefile!\n",func,MODULUS_TYPE,i);
		goto PM1_S1P_READ_RETURN;
	}
	// compare b1 of savefile data vs global B1 presumed to have been set by caller:
 	for(j = 0; j < 32; j += 8) {
		i = fgetc(fptr);	b1 += (uint64)i << j;
	}
	if(B1 != b1) {
		sprintf(cbuf, "INFO: %s: B1 of current run[%u] mismatches one[%u] of savefile data.\n",func,B1,b1);
		goto PM1_S1P_READ_RETURN;
	}
	// Read bitlength of precomputed/bit-reversed product:
	*nbits = 0;
 	for(j = 0; j < 32; j += 8) {
		i = fgetc(fptr);	*nbits += i << j;
	}

	// Set the number of product bytes and zero the corr. target-array limbs:
	nbytes = (*nbits + 7)/8; nlimbs = (nbytes + 7)/8;
	for(i = 0; i < nlimbs; i++) { arr[i] = 0ull; }

	// Read the bytewise product into our array of 64-bit limbs:
	// j holds index of current byte of limb, (j>>3) = index of current limb of target
	for(j = 0; j < nbytes; j++) {				//vvvvvvvvv = 8*j (mod 64)
		c = fgetc(fptr);	arr[j>>3] += ((uint64)c << ((j<<3)&63));
	}
	// Read 8 bytes of simple (sum of limbs, mod 2^64) checksum, compare to one computed from read data:
	for(j = 0; j < 64; j += 8) {
		i = fgetc(fptr);	isum64 += (uint64)i << j;
	}
	for(i = 0; i < nlimbs; i++) { itmp64 += arr[i]; }
	if(itmp64 != isum64) {
		sprintf(cbuf, "INFO: %s: Computed checksum[%" PRIX64 "] mismatches one[%" PRIX64 "] appended to savefile data.\n",func,itmp64,isum64);
		*sum64 = 0ull;
		goto PM1_S1P_READ_RETURN;
	} else {
		*sum64 = isum64;
		retval = nlimbs;
	}
  #endif
PM1_S1P_READ_RETURN:
	if(fptr) { fclose(fptr); fptr = 0x0; }
	return retval;
}

#ifndef PM1_STANDALONE
	// Returns 1 on successful write, 0 otherwise:
	int write_pm1_s1_prod(const char*fname, uint64 p, uint32 nbits, uint64 arr[], uint64 sum64)
	{
		const char func[] = "write_pm1_s1_prod";
		int retval = 0;
		uint8 c;
		uint32 i,j,b1 = 0,nbytes,nlimbs;
		uint64 itmp64 = 0ull;
		ASSERT(arr != 0x0, "Null arr pointer!");
		ASSERT(strlen(fname) != 0, "Empty filename!");

		FILE*fptr = mlucas_fopen(fname, "wb");
		if(!fptr) {
			sprintf(cbuf,"ERROR: Unable to open precomputed p-1 stage 1 primes-product file %s for writing.\n",fname);
			mlucas_fprint(cbuf,pm1_standlone+1);	ASSERT(0, cbuf);
		}
		fprintf(stderr,"INFO: Opened precomputed p-1 stage 1 primes-product file %s for writing...\n",fname);

		fputc(TEST_TYPE,fptr);
		fputc(MODULUS_TYPE,fptr);
		// B1:
		for(j = 0; j < 32; j += 8) {
			fputc(B1>>j,fptr);
		}
		// Write bitlength of precomputed/bit-reversed product:
		for(j = 0; j < 32; j += 8) {
			fputc(nbits>>j,fptr);
		}

		// Set the number of product bytes and limbs:
		nbytes = (nbits + 7)/8; nlimbs = (nbytes + 7)/8;

		// Write the bytewise product stored in our array of 64-bit limbs:
		// j holds index of current byte of limb, (j>>3) = index of current limb of source
		for(j = 0; j < nbytes; j++) {				//vvvvvvvvv = 8*j (mod 64)
			c = arr[j>>3] >> ((j<<3)&63); fputc(c,fptr);
		}
		// Write 8 bytes of simple (sum of limbs, mod 2^64) checksum, after comparing arglist version to one computed from actual data:
		for(i = 0; i < nlimbs; i++) { itmp64 += arr[i]; }
		if(itmp64 != sum64) {
			sprintf(cbuf, "INFO: %s: Computed checksum[%" PRIX64 "] mismatches one[%" PRIX64 "] in arglist.\n",func,itmp64,sum64);
			goto PM1_S1P_WRITE_RETURN;
		}
		for(j = 0; j < 64; j += 8) {
			c = sum64 >> j; fputc(c,fptr);
		}
		retval = 1;

	PM1_S1P_WRITE_RETURN:
		if(fptr) { fclose(fptr); fptr = 0x0; }
		return retval;
	}
#endif

/* Given stage 2 #buffers determined from available RAM and modulus size, returns preferred bigstep size
(available values 210,330 and 420) and extended-prime-pairing-window multiplicity m for stage 2. Data set1 and set2
contain #modmul for simulated stage 2 with B1 = 5e6 and B2 = 150e6 for bigstep = 210 and 330, respectively.
There are 8095883 primes in [5000000,150000000], so as a simple sanity check, the lowest possible #modmul -
neglecting 2-modmul-per-bigstep-loop overhead associated with bigstep-up-multiplies - is 1/2 this = 4047941.
Taking into account said overhead & assuming 100% prime-pairing, lowest #modmul for various bigstep values D:
	D	#loops	min #modmul for [b1,b2] = [5e6,150e6]
	---	------	-------
	210	690476	5429010
	330	439394	4926846
	420	345238	4738534
	660	219697	4487452
	840	172620	4393298
and we see the predictably rapidly vanishing returns for ever-larger D.

Here is the procedure for creating our simple lookup table consisting of [#bufs,bigstep-to-use] pairs:
	1. Combine data for #bufs <= 2400 for bigstep D = [210|330|420|660|840], sort by ascending #bufs;
	2. For each distinct #bufs value, if multiple D hit it exactly, pick the one with fewest #modmul;
	3. For each distinct #bufs value, discard entry if some lower-#bufs value for some other D-value has fewer #modmul:
	in simpler terms, we are ensuring that #modmul is monotone decreasing with increasing #bufs;
	4. For input nbuf (maximum #bufs given available RAM and user's maxalloc setting), pick the entry having largest #buf <= nbuf,
	which given our above table-preprocessing is guaranteed to have the fewest #modmul at the given allowable-memory limit.
	5. For #bufs > 1392, D = 840 always wins, so check if that "large-memory default" holds first.
Here the data resulting from our prime-pairing algo, sorted as per above, in human-readable form - M is the
"extended prime-pairing window multiplicity" associated with each [#buf,D] pair - larger M yields larger pairing %,
but smaller D incur more per-loop up-multiply cost, so we have a complicated optimization problem best addressed
by crunching the actual numbers and letting the chips (modmuls) fall as they may. Dividing 8095883 (#primes in sample
S2 interval [5000000,150000000]) by #loops for each distinct D-value gives average number of stage 2 primes processed
per D-width loop pass as [11.7|18.4|23.5|46.9] for D = [210|330|420|840] (D = 660 was not best at any #buf).
%modmul is relative to min-#buf-24 baseline:
	#buf #modmul %modmul D  nloop  M	%pair		#buf #modmul %modmul D  nloop  M	%pair[cont. from bottom left:]

	  24 8482142 100.0% 210 690476  1	24.61		 672 5388440 63.53% 420 345238 14	83.96
	  40 8047788 94.88% 330 439394  1	22.93		 720 5347048 63.04% 420 345238 15	84.99
	  48 7799105 91.95% 420 345238  1	24.42		 768 5310865 62.61% 420 345238 16	85.88
	  72 7462541 87.98% 210 690476  3	49.79		 816 5278205 62.23% 420 345238 17	86.69
	  80 7440240 87.72% 330 439394  2	37.94		 864 5248799 61.88% 420 345238 18	87.41
	  96 7156961 84.38% 210 690476  4	57.34		 912 5222341 61.57% 420 345238 19	88.07
	 120 6929970 81.70% 210 690476  5	62.94		 960 5197995 61.28% 420 345238 20	88.67
	 144 6748924 79.57% 210 690476  6	67.42		1008 5176019 61.02% 420 345238 21	89.21
	 160 6720220 79.23% 330 439394  4	55.72		1056 5155931 60.79% 420 345238 22	89.71
	 168 6604449 77.86% 210 690476  7	70.98		1104 5137289 60.57% 420 345238 23	90.17
	 192 6487084 76.48% 210 690476  8	73.88		1152 5120242 60.36% 420 345238 24	90.59
	 200 6484045 76.44% 330 439394  5	61.55		1200 5104499 60.18% 420 345238 25	90.98
	 216 6389511 75.33% 210 690476  9	76.29		1248 5088667 59.99% 840 172620 13	82.84
	 240 6238029 73.54% 420 345238  5	62.98		1296 5076212 59.85% 420 345238 27	91.68
	 264 6236660 73.53% 210 690476 11	80.07		1344 5042544 59.45% 840 172620 14	83.98
	 280 6154343 72.56% 330 439394  7	69.70		1440 5001766 58.97% 840 172620 15	84.99
	 288 6058915 71.43% 420 345238  6	67.40		1536 4965014 58.53% 840 172620 16	85.90
	 320 6038492 71.19% 330 439394  8	72.56		1632 4932511 58.15% 840 172620 17	86.70
	 336 5914048 69.72% 420 345238  7	70.98		1728 4903454 57.81% 840 172620 18	87.42
	 384 5796905 68.34% 420 345238  8	73.87		1824 4876997 57.50% 840 172620 19	88.07
	 432 5699150 67.19% 420 345238  9	76.29		1920 4853017 57.21% 840 172620 20	88.66
	 480 5617042 66.22% 420 345238 10	78.32		2016 4830974 56.95% 840 172620 21	89.21
	 528 5547272 65.40% 420 345238 11	80.04		2112 4810589 56.71% 840 172620 22	89.71
	 576 5486926 64.69% 420 345238 12	81.53		2208 4791914 56.49% 840 172620 23	90.17
	 624 5434492 64.07% 420 345238 13	82.83		2304 4775023 56.30% 840 172620 24	90.59
			[continued at upper right:]				2400 4759140 56.11% 840 172620 25	90.98
6/09/21:
Here the above table redone with small-prime relocation enabled - %modmul still relative to above 24-buf count, 8482142.
Note how D = 660 is now competitive with D = 840 at larger #bufs values, since it has psmall = 7 (vs 840 with psmall = 11)
and thus benefits proportionally more from relocation:
	#buf #modmul %modmul D  nloop  M	%pair		#buf #modmul %modmul D  nloop  M	%pair[cont. from bottom left:]
	  24 8345678 98.39% 210 649351  1	25.94		 720 5268017 62.11% 420 324675 15	85.92
	  40 7824013 92.24% 330 389611  1	26.00		 768 5233484 61.70% 420 324675 16	86.77
	  48 7697860 90.75% 420 324675  1	25.91		 816 5202458 61.33% 420 324675 17	87.54
	  72 7300446 86.07% 210 649351  3	51.76		 864 5174553 61.01% 420 324675 18	88.23
	  80 7195798 84.83% 330 389611  2	41.52		 912 5149358 60.71% 420 324675 19	88.85
	  96 6999318 82.52% 210 649351  4	59.20		 960 5126327 60.44% 420 324675 20	89.42
	 120 6776892 79.90% 210 649351  5	64.69		1008 5105708 60.19% 420 324675 21	89.93
	 144 6601588 77.83% 210 649351  6	69.02		1040 5092650 60.04% 660 194805 13	83.84
	 160 6482198 76.42% 330 389611  4	59.14		1056 5086617 59.97% 420 324675 22	90.40
	 168 6461791 76.18% 210 649351  7	72.47		1104 5068966 59.76% 420 324675 23	90.84
	 192 6349082 74.85% 210 649351  8	75.26		1120 5048281 59.52% 660 194805 14	84.93
	 200 6258403 73.78% 330 389611  5	64.67		1200 5008982 59.05% 660 194805 15	85.90
	 216 6254572 73.74% 210 649351  9	77.59		1280 4974295 58.64% 660 194805 16	86.76
	 240 6083777 71.72% 330 389611  6	68.98		1360 4943217 58.28% 660 194805 17	87.53
	 280 5945213 70.09% 330 389611  7	72.40		1440 4915478 57.95% 660 194805 18	88.21
	 320 5831885 68.75% 330 389611  8	75.20		1520 4889742 57.65% 660 194805 19	88.85
	 336 5812502 68.53% 420 324675  7	72.47		1600 4866901 57.38% 660 194805 20	89.41
	 360 5737215 67.64% 330 389611  9	77.54		1680 4845955 57.13% 660 194805 21	89.93
	 384 5699633 67.20% 420 324675  8	75.26		1760 4826928 56.91% 660 194805 22	90.40
	 400 5658388 66.71% 330 389611 10	79.49		1824 4824707 56.88% 840 162338 19	88.85
	 432 5605740 66.09% 420 324675  9	77.58		1840 4809442 56.70% 660 194805 23	90.83
	 440 5590849 65.91% 330 389611 11	81.16		1920 4793157 56.51% 660 194805 24	91.23
	 480 5526455 65.15% 420 324675 10	79.54		2000 4778202 56.33% 660 194805 25	91.60
	 520 5481846 64.63% 330 389611 13	83.85		2080 4764475 56.17% 660 194805 26	91.94
	 528 5459529 64.36% 420 324675 11	81.19		2112 4761846 56.14% 840 162338 22	90.40
	 560 5437458 64.10% 330 389611 14	84.94		2160 4751637 56.02% 660 194805 27	92.26
	 576 5401805 63.68% 420 324675 12	82.62		2208 4744211 55.93% 840 162338 23	90.84
	 600 5398691 63.65% 330 389611 15	85.90		2240 4739863 55.88% 660 194805 28	92.55
	 624 5351372 63.09% 420 324675 13	83.86		2304 4728299 55.74% 840 162338 24	91.23
	 672 5307229 62.57% 420 324675 14	84.95		2400 4713322 55.57% 840 162338 25	91.60
			[continued at upper right:]

And here the data for the smaller stage 2 interval B1 = 1m, B2 = 30m, which contains 1779361 primes:
No small-prime relocation:
	#buf #modmul %modmul D  nloop  M	%pair		#buf #modmul %modmul D  nloop  M	%pair[cont. from bottom left:]
	  24 1818134 100.0% 210 138096  1	26.76		 720 1145216 62.98% 420  69048 15	86.85
	  40 1729614 95.13% 330  87880  1	25.42		 768 1138021 62.59% 420  69048 16	87.65
	  48 1679385 92.37% 420  69048  1	26.83		 816 1131353 62.23% 420  69048 17	88.40
	  72 1583838 87.11% 210 138096  3	53.08		 864 1125507 61.90% 420  69048 18	89.06
	  80 1591151 87.52% 330  87880  2	40.98		 912 1120118 61.61% 420  69048 19	89.67
	  96 1517322 83.45% 210 138096  4	60.56		 960 1115306 61.34% 420  69048 20	90.21
	 120 1468182 80.75% 210 138096  5	66.08		1008 1111002 61.11% 420  69048 21	90.69
	 144 1430020 78.65% 210 138096  6	70.36		1040 1120045 61.60% 660  43940 13	84.03
	 160 1430793 78.70% 330  87880  4	59.00		1056 1107002 60.89% 420  69048 22	91.14
	 168 1399838 76.99% 210 138096  7	73.75		1104 1103278 60.68% 420  69048 23	91.56
	 192 1375065 75.63% 210 138096  8	76.54		1120 1110435 61.08% 660  43940 14	85.12
	 200 1380361 75.92% 330  87880  5	64.66		1200 1096804 60.32% 420  69048 25	92.28
	 216 1354857 74.52% 210 138096  9	78.81		1280 1094055 60.17% 660  43940 16	86.96
	 240 1341246 73.77% 330  87880  6	69.06		1360 1087212 59.80% 660  43940 17	87.72
	 280 1310632 72.09% 330  87880  7	72.50		1440 1081064 59.46% 660  43940 18	88.41
	 320 1285829 70.72% 330  87880  8	75.28		1520 1075396 59.15% 660  43940 19	89.05
	 336 1261647 69.39% 420  69048  7	73.76		1600 1070283 58.87% 660  43940 20	89.63
	 360 1264910 69.57% 330  87880  9	77.63		1680 1065757 58.62% 660  43940 21	90.13
	 384 1237000 68.04% 420  69048  8	76.53		1760 1061719 58.40% 660  43940 22	90.59
	 400 1247389 68.61% 330  87880 10	79.60		1824 1050962 57.80% 840  34525 19	89.68
	 432 1216682 66.92% 420  69048  9	78.82		1840 1057952 58.19% 660  43940 23	91.01
	 440 1232358 67.78% 330  87880 11	81.29		1920 1054345 58.00% 660  43940 24	91.42
	 480 1199748 65.99% 420  69048 10	80.72		2000 1051033 57.81% 660  43940 25	91.79
	 520 1208462 66.47% 330  87880 13	83.97		2080 1048018 57.64% 660  43940 26	92.13
	 528 1185417 65.20% 420  69048 11	82.33		2112 1037777 57.08% 840  34525 22	91.16
	 560 1198362 65.91% 330  87880 14	85.11		2160 1045178 57.49% 660  43940 27	92.44
	 576 1173301 64.53% 420  69048 12	83.69		2208 1034068 56.88% 840  34525 23	91.58
	 600 1189946 65.45% 330  87880 15	86.06		2240 1042536 57.34% 660  43940 28	92.74
	 624 1162754 63.95% 420  69048 13	84.87		2304 1030679 56.69% 840  34525 24	91.96
	 672 1153389 63.44% 420  69048 14	85.93		2400 1027564 56.52% 840  34525 25	92.31
			[continued at upper right:]
With small-prime relocation:
	#buf #modmul %modmul D  nloop  M	%pair		#buf #modmul %modmul D  nloop  M	%pair[cont. from bottom left:]
	  24 1786179 98.24% 210 129871  1	28.50		 720 1129178 62.11% 420  64935 15	87.72
	  40 1681487 92.48% 330  77923  1	28.59		 768 1122364 61.73% 420  64935 16	88.49
	  48 1656511 91.11% 420  64935  1	28.47		 816 1116072 61.39% 420  64935 17	89.19
	  72 1549829 85.24% 210 129871  3	55.05		 864 1110516 61.08% 420  64935 18	89.82
	  80 1537437 84.56% 330  77923  2	44.77		 912 1105479 60.80% 420  64935 19	90.38
	  96 1484343 81.64% 210 129871  4	62.41		 960 1100884 60.55% 420  64935 20	90.90
	 120 1436840 79.03% 210 129871  5	67.75		1008 1096861 60.33% 420  64935 21	91.35
	 144 1399917 77.00% 210 129871  6	71.90		1040 1093873 60.16% 660  38961 13	85.85
	 160 1380408 75.92% 330  77923  4	62.42		1056 1093087 60.12% 420  64935 22	91.78
	 168 1370744 75.39% 210 129871  7	75.17		1104 1089522 59.93% 420  64935 23	92.17
	 192 1347022 74.09% 210 129871  8	77.84		1120 1085080 59.68% 660  38961 14	86.84
	 200 1332531 73.29% 330  77923  5	67.79		1200 1077231 59.25% 660  38961 15	87.72
	 216 1327505 73.01% 210 129871  9	80.03		1280 1070249 58.87% 660  38961 16	88.50
	 240 1295895 71.28% 330  77923  6	71.91		1360 1064021 58.52% 660  38961 17	89.20
	 280 1266916 69.68% 330  77923  7	75.17		1440 1058455 58.22% 660  38961 18	89.83
	 320 1243416 68.39% 330  77923  8	77.81		1520 1053409 57.94% 660  38961 19	90.39
	 336 1240671 68.24% 420  64935  7	75.19		1600 1048854 57.69% 660  38961 20	90.91
	 360 1223952 67.32% 330  77923  9	79.99		1680 1044726 57.46% 660  38961 21	91.37
	 384 1217030 66.94% 420  64935  8	77.85		1760 1041054 57.26% 660  38961 22	91.78
	 400 1207768 66.43% 330  77923 10	81.81		1824 1040554 57.23% 840  32468 19	90.38
	 432 1197594 65.87% 420  64935  9	80.03		1840 1037641 57.07% 660  38961 23	92.16
	 440 1193893 65.67% 330  77923 11	83.37		1920 1034445 56.90% 660  38961 24	92.52
	 480 1181404 64.98% 420  64935 10	81.85		2000 1031511 56.73% 660  38961 25	92.85
	 520 1172081 64.47% 330  77923 13	85.82		2080 1028865 56.59% 660  38961 26	93.15
	 528 1167552 64.22% 420  64935 11	83.41		2112 1028148 56.55% 840  32468 22	91.78
	 560 1162990 63.97% 330  77923 14	86.84		2160 1026396 56.45% 660  38961 27	93.43
	 576 1155956 63.58% 420  64935 12	84.71		2208 1024651 56.36% 840  32468 23	92.17
	 600 1155257 63.54% 330  77923 15	87.71		2240 1024091 56.32% 660  38961 28	93.69
	 624 1145852 63.02% 420  64935 13	85.85		2304 1021474 56.18% 840  32468 24	92.53
	 672 1136955 62.53% 420  64935 14	86.85		2400 1018574 56.02% 840  32468 25	92.85
			[continued at upper right:]
*/
// psmall stores any relocation-prime used for a previously-started but interrupted stage 2. On restart seed the call
// to pm1_bigstep_size() with that to ensure our restart-run bigstep shares the same relocation-prime:
void pm1_bigstep_size(uint32*nbuf, uint32*bigstep, uint32*m, uint32 psmall)
{
	// If psmall = 0, any of the above bigsteps D is in play - Use the last of the above tables, laid out in [#bufs,bigstep] pairs:
	const uint32 lut_all[120] = {
	  24,210,   40,330,   48,420,   72,210,   80,330,   96,210,  120,210,  144,210,  160,330,  168,210,  192,210,  200,330,
	 216,210,  240,330,  280,330,  320,330,  336,420,  360,330,  384,420,  400,330,  432,420,  440,330,  480,420,  520,330,
	 528,420,  560,330,  576,420,  600,330,  624,420,  672,420,  720,420,  768,420,  816,420,  864,420,  912,420,  960,420,
	1008,420, 1040,660, 1056,420, 1104,420, 1120,660, 1200,660, 1280,660, 1360,660, 1440,660, 1520,660, 1600,660, 1680,660,
	1760,660, 1824,840, 1840,660, 1920,660, 2000,660, 2080,660, 2112,840, 2160,660, 2208,840, 2240,660, 2304,840, 2400,840};
	// If psmall = 7, restrict ourselves to D = 330,660:
	const uint32 lut_psmall7[84] = {
	  40,330,   80,330,  120,330,  160,330,  200,330,  240,330,  280,330,  320,330,  360,330,  400,330,  440,330,  480,330,
	 520,330,  560,330,  600,330,  640,330,  680,330,  720,330,  760,330,  800,330,  840,330,  880,660,  920,330,  960,660,
	1040,660, 1120,660, 1200,660, 1280,660, 1360,660, 1440,660, 1520,660, 1600,660, 1680,660, 1760,660, 1840,660, 1920,660,
	2000,660, 2080,660, 2160,660, 2240,660, 2320,660, 2400,660};
	// If psmall = 11, restrict ourselves to D = 210,420,840:
	const uint32 lut_psmall11[88] = {
	   24,210,   48,420,   72,210,   96,210,  120,210,  144,210,  168,210,  192,210,  216,210,  240,420,  264,210,  288,420,
	  336,420,  384,420,  432,420,  480,420,  528,420,  576,420,  624,420,  672,420,  720,420,  768,420,  864,420,  912,420,
	  960,420, 1008,420, 1056,420, 1104,420, 1152,420, 1200,420, 1248,420, 1296,420, 1344,840, 1440,840, 1536,840, 1632,840,
	 1728,840, 1824,840, 1920,840, 2016,840, 2112,840, 2208,840, 2304,840, 2400,840,};
	const uint32*lut = 0x0;
	if(!psmall)
		lut = lut_all;
	else if(psmall == 7)
		lut = lut_psmall7;
	else if(psmall == 11)
		lut = lut_psmall11;
	else
		ASSERT(0, "pm1_bigstep_size: Bad input value of relocation-prime!");
	// High-RAM case - For given D and associated num_b, M = floor(nbuf/num_b), where num_b = 24|40|48|80|96
	// for D = 210|330|420|660|840. Only need to special-case psmall = 7 here, all others use D = 840:
	if(*nbuf >= 10000) {
		sprintf(cbuf,"WARNING: %u buffers requested; Stage 2 allows a maximum of 10000.\n",*nbuf);
		mlucas_fprint(cbuf,pm1_standlone+1);
		*nbuf = 10000;
	}
	if(*nbuf >= 2400) {
		if(psmall == 7) {
			*bigstep = 660;	*m = *nbuf/80; *nbuf = *m * 80;
		} else {
			*bigstep = 840;	*m = *nbuf/96; *nbuf = *m * 96;
		}
		return;
	}
	int i, num_b = 0;
	for(i = 0; ; i += 2) {	// Foregoing assert ensures i gets incremented at least once
		if(lut[i] > *nbuf) break;
	}
	if(!i)
		ASSERT(0, "P-1 stage 2 with relocation prime psmall = 7|11 needs at least 40|24 buffers of available RAM, respectively!");
	if(psmall) {
		sprintf(cbuf,"Previous Stage 2 work used relocation-prime %u ... enforcing compatibility with this: bigstep must be a multiple of %u.\n",psmall,18-psmall);
		mlucas_fprint(cbuf,pm1_standlone+1);
		// Here's why we don't declare psmall const in the arglist - it stores the smallest prime which does
		// not divide the bigstep value, in order to check divisibility replace it by its complement here:
		psmall = 18-psmall;
		ASSERT((lut[i-1]%psmall == 0), "P-1 stage 2 needs at least 24 buffers of available RAM!");
		/* First-go-round of this used just a single unified lut[] array and worked backward to the largest nbuf whoe D is compatible:
		for( ; ; i -= 2) {
			if(lut[i-1]%psmall == 0) break;
		}*/
	}
	// For each lut[] pair, lut[i] stores nbuf, lut[i+1] stores bigstep, but above loops' structure ensures
	// that we overshoot nbuf, i.e. on exit want to set nbuf = lut[i-2], bigstep = lut[i-1]:
	*bigstep = lut[i-1];	*nbuf = lut[i-2];
	// For a given #bufs value, the paired bigstep D yields an extended-prime-pairing-window multiplicity #bufs/num_b,
	// where num_b = #integers in [1,D/2-1] which are coprime to D and = [24|40|48|80|96] for D = [210|330|420|660|840]:
	for(i = 1; i < (*bigstep>>1); i++) {
		num_b += (gcd32(*bigstep,i) == 1);
	}
	*m = *nbuf/num_b;
	return;
}

/* In-place binary modpow, a^pow (mod [M(p) or F(lg2(p))]), using length-n FFT. Requires pow != 0.
Whether Input|Output vector a[] needs fwd-DWT weighting and pass 1 of fwd-FFT done on entry
and undone rior to return are controlled by low 2 bits of mode_flag arg, as detailed in the
long comment about this flag at top of [mers|ferm]_mod_square.c.

Assumes function-pointer arg has been set to point to either [mers|ferm]_mod_square by caller.
Our crappy modmul-API forces caller to specify an "input is pure-int?" flag.
Needs pointer to scratch storage vector b[] to be provided.

Result is left in fwd-FFT-pass1-done form, and overwrites input a[].
Returns nonzero if any of the modmuls returned an error code.
*/
int modpow(double a[], double b[], uint32 input_is_int, uint64 pow,
	int	(*func_mod_square)(double [], int [], int, int, int, uint64, uint64, int, double *, int, double *),
	uint64 p, int n, int scrnFlag, double *tdiff)
{
	*tdiff = 0.0;
	double tdif2 = 0.0;	// Need 2 timers here - tdif2 for the individual func_mod_square calls, accumulate in tdiff
	uint64 npad = n + ( (n >> DAT_BITS) << PAD_BITS ), nbytes = npad<<3;	// npad = length of padded data array
	// v20.1.1: Made npad, nbytes 64-bit here and in pm1_stage2() to allow FFT >= 512M. (Could also do nbytes = (uint64)npad<<3)
	uint32 i,ierr = 0,nerr = 0, mode_flag;
#ifdef PM1_DEBUG
	uint32 j; double dsum;
#endif
	ASSERT(a && b && n && func_mod_square,"Null input pointer or vector length in pm1.c::modpow!");
	// pow = 0: , b[1:n-1] = 0:
	if(!pow) {
		b[0] = 1.0;
		for(i = 1; i < npad; i++) { b[i] = 0.0; }
		return 0;
	}
	// pow = 1 is mathematically a no-op, but if a != b need to copy a into b:
	if(pow == 1 && a != b) {
		memcpy(b,a,nbytes);	// Assume disjoint memory footprints here
		return 0;
	}
	// Init b = fwdFFT(a); only need this if power != 2^k, in which case we only need autosquarings:
#ifdef PM1_DEBUG
	fprintf(stderr,"MODPOW: pow = %" PRIu64 "\n",pow);
#endif
	if(!isPow2_64(pow)) {
		memcpy(b,a,nbytes);	// b = a           vvvv + 4 to effect "Do in-place forward FFT only; low bit = 0 here implies pure-int input"
		ierr = func_mod_square(b, 0x0, n, 0,1, 4ull + (uint64)(!input_is_int), p, scrnFlag,&tdif2, FALSE, 0x0); *tdiff += tdif2; nerr += ierr;
		if(ierr == ERR_INTERRUPT) {
			return ierr;
		}
	#ifdef PM1_DEBUG
		dsum = 0; for(j = 0; j < npad; j++) { dsum += fabs(b[j]); }; fprintf(stderr,"b = fwdFFT(a) gives b[0] = %20.8f, b[1] = %20.8f, L1(b) = %20.8f\n",b[0],b[1],dsum/n); MME = 0;
	#endif
	}	ASSERT(nerr == 0, "func_mod_square returns error!");
	// Use LR binary modpow algorithm, though it's no faster in this general case than RL:
	uint32 len = nbits64(pow);
	pow = reverse64(pow,len)>>1;	// Leftmost bit of input power accounted for by implied (result = a[]) init;
									// This is why we start following loop at i = 1:
	/*
	Ex: a^35: len = 6, pow = 110001>>1, y=a, then loop:
	i	pow		len	action		y-value
	1	11000	5	y*=y		a^2
	2	1100	4	y*=y		a^4
	3	110		3	y*=y		a^8
	4	11		2	y*=y,y*=a	a^17
	5	1		1	y*=y,y*=a	a^35
	*/
	for(i = 1; i < len; i++, pow>>=1) {
		// input_is_int specifies whether a[] holds pure-int values or fwd-FFT-pas1-already-done, but only use on first call (i == 1):
		// bit 1 of mode_flag always = 1 since all FFT-mul outputs assumed to be getting re-used as inputs
		if(i == 1)	// i == 1: bit 0 = !input_is_int
			mode_flag = 2 + !input_is_int;
		else		// i != 1: bit 0 = 1
			mode_flag = 3;
		// y *= y:
		ierr = func_mod_square(a, 0x0, n, i,i+1, (uint64)mode_flag, p, scrnFlag,&tdif2, FALSE, 0x0); *tdiff += tdif2; nerr += ierr;
		if(ierr == ERR_INTERRUPT) {
			return ierr;
		}
	#ifdef PM1_DEBUG
		dsum = 0; for(j = 0; j < npad; j++) { dsum += fabs(a[j]); }; fprintf(stderr,"a^2: MME = %8.6f, a[0] = %20.8f, a[1] = %20.8f, L1(a) = %20.8f\n",MME,a[0],a[1],dsum/n); MME = 0;
	#endif
	  if(pow&1)	{	// y *= a; mode_flag for this fixed = 3:
		ierr = func_mod_square(a, 0x0, n, i,i+1, (uint64)b +  3ull, p, scrnFlag,&tdif2, FALSE, 0x0); *tdiff += tdif2; nerr += ierr;
		if(ierr == ERR_INTERRUPT) {
			return ierr;
		}
	  #ifdef PM1_DEBUG
		dsum = 0; for(j = 0; j < npad; j++) { dsum += fabs(a[j]); }; fprintf(stderr,"a*b: MME = %8.6f, a[0] = %20.8f, a[1] = %20.8f, L1(a) = %20.8f\n",MME,a[0],a[1],dsum/n); MME = 0;
	  #endif
	  }
	}
	// For initial release, no error handling - note we do have ROE handling in the main Stage 2 loop:
	if(nerr != 0) { sprintf(cbuf,"modpow hit one or more errors! Aborting."); ASSERT(0,cbuf); }
	// Result returned in a[]:
	return (nerr != 0);
}

/*
Savefile scheme for p-1 S1 and S2, considerations:

o To support distributed S2 runs we need to save S1 residue.
o S2 needs its own savefile, keep it simple for now, just do 1 gcd on completion of each stage.
o S1 uses regular-named savefile pair, with TEST_TYPE byte set to TEST_TYPE_PM1.
o S2 code appends ".s2" to p-or-f-prefix primary-savefile name & does its own savefile r+w, also with TEST_TYPE_PM1.

Here the overall p-1 savefile schema:

1. on (Re)start, we know it's p-1 via the worktodo assignment: read S1 savefile, see if S1 done
based on iteration count versus PM1_S1_PROD_BITS as computed from the B1 bound, then either:
	[a] If S1 incomplete (based on maxiter < PM1_S1_PROD_BITS), resume;
	[b] If S1 completed, including GCD, feed residue to S2.
2. S2 needs to do S1-res powering in both from-scratch and S2-restart-from-interrupt cases.
3. On entry, S2 will check for existence of ".s2" savefile; If exists, read nsquares, resume S2 there.

	6/08/21: Added small-prime relocation for the most-common use case of a contiguous (B2_start = B1 at
	start of stage 2) S2 continuation, in which we run S2 not from (B2_start = B1) to B2 but rather from
	B1*psmall to B2, where psmall = smallest prime which is not a factor of our bigstep D. In our case,
	that means psmall = 7 (for D = 330|660) or psmall = 11 (for D = 210|420|840). In the context of this
	relocation scheme, any restart - which may encounter differing avail-mem than initial S2 start and
	thus want to choose different D (and if a Pfactor= assignment, different B2) needs to be constrained
	to use a D with the same psmall as the preceding S2 work, and must also use the same B2.
4. Every so often, S2 writes savefile using end-of-k-loop value of residue 'pow';
5. On S2 completion, returns S2 bytewise residue in arrtmp.
*/
#ifdef PM1_DEBUG
	#warning Building p-1 code in PM1_DEBUG mode.
	const char asterisk_if_false[2] = {'*',' '};
#endif
#ifdef PM1_STANDALONE
  #warning Building pm1_stage2() in standalone (modmul-counting) mode!
  // Standalone mode proof-of-principle code is based on my original bc scripts, and used to test the q-pairing algo for large B2.
  int main(int argc, char *argv[])
  {
	uint32 bigstep = 0, m = 0, n = -1;	// Set n != 0 to trigger mem-alloc/init block
	uint32 interim_gcd = 0, q_old_10M = 0;
	const double inv10m = 1./10000000;
#else
	/* Production-mode stage 2 subroutines - assume stage bounds stored in externs B1,B2_start,B2 on entry,
	stage 1 powering residue in pure-integer-valued-double[] form supplied in pow[], and mem for needed
	m*[24 or 40] stage 2 buffers assumed alloc'ed via buf[], each element of which is assumed to point
	to a residue-length subarray. n = unpadded FFT length in #doubles.
	Assumes function-pointer arg has been set to point to either [mers|ferm]_mod_square by caller:
	*/
  int pm1_stage2(uint64 p, uint32 bigstep, uint32 m,
	double pow[],	// Stage 1 residue, overwritable after needed inits if needed
	double*mult[],	// scratch storage, in form of ptrs to 4 residue-length subarrays
	uint64 arrtmp[],	// scratch storage for bytewise savefile residues
	int	(*func_mod_square)(double [], int [], int, int, int, uint64, uint64, int, double *, int, double *),
	int n, int scrnFlag, double *tdiff, char*const gcd_str)
  {
	const double inv10m = 1./10000000; uint32 q_div_10M, q_old_10M = 0;	// For intermediate GCD-scheduling
	gcd_str[0] = '\0';	// should be null on entry, but better safe than sorry
	char savefile[STR_MAX_LEN];	// S1 residue-inverse file
	uint32 restart = 0, input_is_int = 0, mode_flag = 0, kblocks = (n>>10), interim_gcd = 1;
	// npad = length of padded data array:
	uint64 Res64,Res35m1,Res36m1, nalloc, npad = n + ( (n >> DAT_BITS) << PAD_BITS ), nbytes = npad<<3;
	double dtmp;
#if USE_PP1_MULTS		// (p+1)-style version of stage 2 loop-multipliers
	char inv_file[STR_MAX_LEN];	// S1 residue-inverse file
	uint32 s1_inverse = 0;
	const uint64 two35m1 = (uint64)0x00000007FFFFFFFFull, two36m1 = (uint64)0x0000000FFFFFFFFFull;	/* 2^35,36-1 */
#endif
  #define SIZE 256
	time_t calendar_time;
	struct tm *local_time;
	char timebuffer[SIZE];
  #ifdef CTIME
	clock_t clock1, clock2;
  #else
	double clock1, clock2;
  #endif
	double tdif2 = 0.0;	// tdif2 is dummy-arg for the time-ptr field of the modmul calls; as that neglects time used
						// for memcpy and FFT(a-b) calls, use arglist-ptr tdiff to accumulate overall elapsed-time.
	// Simple threadpool-struct for parallel FFT(a-b) = FFT(a) - FFT(b) used in Stage 2 loop:
  #ifdef MULTITHREAD
	static int *thr_ret = 0x0;
	static pthread_t *thread = 0x0;
	static pthread_attr_t attr;
	static struct pm1_thread_data_t *tdat = 0x0;
	// Threadpool-based dispatch:
	static struct threadpool *tpool = 0x0;
	static thread_control_t thread_control = {0,0,0};
  #endif
#endif	// #ifdef PM1_STANDALONE
	const char func[] = "pm1_stage2";
	char stFlag[STR_MAX_LEN];
	int retval = 0, i,j,jhi,jeven;	// Make i,j signed to allow for downward-running loop indices & loop control
	// num_b is #buffers per unit of extended-pairing-window size M; wsize is #bytes needed per 'word' of the associated bitmap
	uint32 bigstep_pow2,rsize, nq,np=0,ns=0,ierr,nerr,m2,m_is_odd,m_is_even,num_b,psmall,wsize, k,k0=0, nmodmul = 0,nmodmul_save = 0, p1,p2;
	uint32 word,bit;
	uint64 tmp,q,q0,q1,q2, qlo = 0ull,qhi, reloc_start, pinv64 = 0ull;
	// map_lo|hi intended as variable ptrs to various parts of map[], lo|hi as const ptrs to words beyond end of
	// "working map". Alas, since we alloc map[] at runtime, we can't actually declare hi|lo as const-ptrs. File under
	// "stupid C tricks" - lack of a "set once at runtime but compiler/OS treat as const subsequently" declaration option:
	static uint8 *map = 0x0,*lo = 0x0,*hi = 0x0; uint8 *map_lo,*map_hi, *rmap;
	static uint64 *vec1 = 0ull, *vec2 = 0ull;
	static uint32 nlimb, *b = 0x0;
	static double *a_ptmp = 0x0, *a = 0x0;	// Storage for stage 2 in form of ptrs to m*[24,40 or 48] residue-length subarrays, depending
	static double **buf = 0x0;	// on whether bigstep = [210,330,420 or 840]. a[] is simply a tmp-ptr used to init buf[]
								// and then as an alias for mult[3]. 	<**** NOTE! ****
	static double *vone = 0x0;	// Holds fwdFFT(1)
  #ifdef macintosh
	argc = ccommand(&argv);			/* Macintosh CW */
  #endif
  #ifdef PM1_STANDALONE
	/******** command-line-argument processing while() loop: ********/
	uint32 nargs = 1;
	while(argv[nargs]) {
		strncpy(stFlag, argv[nargs++], STR_MAX_LEN);
		if(stFlag[0] != '-') {
			fprintf(stderr, "*** ERROR: Illegal command-line option %s\n", stFlag);
			fprintf(stderr, "*** The required command-line options are -bigstep [210|330|420|660|840] -b1 [int > 0] -b2 [int > 0] -m [int > 0]\n");	return 1;
		}

		if(STREQ(stFlag, "-bigstep")) {
			strncpy(stFlag, argv[nargs++], STR_MAX_LEN);	bigstep = atoi(stFlag);
		} else if(STREQ(stFlag, "-b1")) {
			strncpy(stFlag, argv[nargs++], STR_MAX_LEN);	B1 = atoi(stFlag);
		} else if(STREQ(stFlag, "-b2")) {
			strncpy(stFlag, argv[nargs++], STR_MAX_LEN);	B2 = atoi(stFlag);
		} else if(STREQ(stFlag, "-m")) {
			strncpy(stFlag, argv[nargs++], STR_MAX_LEN);	m = atoi(stFlag);
		} else {
			fprintf(stderr, "*** ERROR: Unrecognized flag %s.\n", stFlag);	return 1;
		}
	}
	ASSERT(bigstep && B1 && B2 && m, "All 4 args bigstep,b1,b2,m must be set > 0!");
	B2_start = (uint64)B1;
  #else
	// Check function pointer to [mers|fermat]_mod_square based on modulus type:
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		ASSERT(func_mod_square == mers_mod_square  , "Mod-square function pointer incorrectly set in pm1_stage2!");
	else if(MODULUS_TYPE == MODULUS_TYPE_FERMAT)
		ASSERT(func_mod_square == fermat_mod_square, "Mod-square function pointer incorrectly set in pm1_stage2!");
	else
		ASSERT(0,"Modulus type not set in pm1_stage2!");
  #endif

  #ifndef PM1_STANDALONE
	// See if user opted out of doing interim GCDs (the default is interim_gcd = 1, i.e. do them):
	dtmp = mlucas_getOptVal(MLUCAS_INI_FILE,"InterimGCD");	// Any failure-to-find-or-parse can be checked for via isNaN(dtmp)
	if(dtmp == 0) {
		sprintf(cbuf,"User set InterimGCD = 0 in %s ... will only take GCD after completion of stage 2.\n",MLUCAS_INI_FILE);
		mlucas_fprint(cbuf,pm1_standlone+1);
		interim_gcd = 0;
	}
  #endif
	// Whichever bigstep values we want to allow here, the resulting num_b must be a multiple of 8:
	if(bigstep == 210) {		// 2.3.5.7
		num_b = 24;	bigstep_pow2 = 1; rsize = 3; psmall = 11;	// psmall = smallest prime which is not a factor of D
	} else if(bigstep == 330) {	// 2.3.5.11						// bigstep_pow2 = power of 2 in D; rsize = #bytes in relocation-maplets
		num_b = 40;	bigstep_pow2 = 1; rsize = 5; psmall =  7;
	} else if(bigstep == 420) {	// 2^2.3.5.7
		num_b = 48;	bigstep_pow2 = 2; rsize = 3; psmall = 11;
	} else if(bigstep == 660) {	// 2^2.3.5.11
		num_b = 80;	bigstep_pow2 = 2; rsize = 5; psmall =  7;
	} else if(bigstep == 840) {	// 2^3.3.5.7
		num_b = 96;	bigstep_pow2 = 3; rsize = 3; psmall = 11;
	} else {
		fprintf(stderr, "*** ERROR: -bigstep arg must be one of [210|330|420|660|840]; user entered %s.\n", stFlag);	return 1;
	}
	wsize = num_b>>2;	// bitmap word has 2*num_b bits and num_b>>2 bytes
	/*
	Small-prime relocation: If B2/psmall >= B1, we map the primes q in [B1,B2/psmall] to [B1*psmall,B2], run
	stage 2 not from B1 to B2 but from B1*psmall to B2, and incorporate them into the stage 2 accumulator when
	each corresponding semiprime q*psmall is hit in the remapped stage 2 interval:
	*/
	pm1_check_bounds();	// This sanity-checks the bounds and sets B2_start = B1 if unset.
	if(B2_start == B1 && B2 >= B1*psmall) {	// It's possible user running a S2 with B2/psmall < B1, hence the 2nd clause
		reloc_start = psmall*B1;				// Start including relocation-semiprimes once S2 passes this point, and
		B2_start = MIN(reloc_start, B2/psmall);	// shift B2_start upward to reflect the fact that primes in [B1,B2/psmall]
												// will be relocated to [B1*psmall,B2].
	} else {	// In the case of a standalone S2 interval (B2_small > B1), set psmall = 0 and reloc_start = UINT64_MAX:
		psmall = 0; reloc_start = -1ull;
	}
	sprintf(cbuf,"Using B2_start = %" PRIu64 ", B2 = %" PRIu64 ", Bigstep = %u, M = %u\n",B2_start,B2,bigstep,m);
	mlucas_fprint(cbuf,pm1_standlone+1);
	uint32 reloc_on = FALSE;	// Gets switched to TRUE (= start using semiprimes which are multiples of psmall) when q > reloc_start

	// Oct 2021: For small q0 and large #bufs, qlo can underflow, so check!
	q0 = B2_start + bigstep - B2_start%bigstep;
	if(q0 <= (m/2+1)*(uint64)bigstep) {	// Don't have m2 yet
		// Max #bufs computed by recasting if() expression as equality and solving for m, then decrementing m to
		// ensure result m satisfies q0 < (m/2+1)*bigstep (note strictly less-than!) and using PM1_S2_NBUF = m*num_b:
		tmp = (q0/(uint64)bigstep - 1)<<1;	// tmp holds m_max as a uint64
		// For this condition to be hit implies q0 quite small, but makes sure resulting m is 32-bit anyway:
		if((uint64)m < tmp) {
			sprintf(cbuf, "Nonsensical value of M_max = %" PRIu64 " in qlo-underflow check ... aborting.",tmp);
			mlucas_fprint(cbuf,pm1_standlone+1);	ASSERT(0,cbuf);
		}
		m = tmp-1;
		PM1_S2_NBUF = m*num_b;	// Don't use PM1_S2_NBUF per se in code below, but reset for consistency
		sprintf(cbuf, "WARNING: The stage 1 bound %u allows a maximum #buffers = %u; using that.\n",B1,PM1_S2_NBUF);
		mlucas_fprint(cbuf,pm1_standlone+1);
	}

	// May 2021: Added support for M even:
	m_is_odd = IS_ODD(m);
	m_is_even = !m_is_odd;
	ASSERT(RES_SHIFT == 0ull, "Shifted residues unsupported for p-1!\n");	// Need BASE_MULTIPLIER_BITS array = 0 for modmuls below!
	// Alloc the needed memory:
  #ifndef PM1_STANDALONE
	nlimb = (p+63+(MODULUS_TYPE == MODULUS_TYPE_FERMAT))>>6;	// # of 64-bit limbs in p-bit vector, alloc 2 of these for debug:
  #endif
/******************************************************************************************************/
/*** EVERY ERROR-RETURN FROM HERE ONWARD MUST DEALLOCATE THE ALLOCATED MEMORY, I.E. GOTO ERR_RETURN ***/
/******************************************************************************************************/
	// Init the needed array storage:
  #ifndef PM1_STANDALONE
	// Double arrays - need (num_b*m + [defined(USE_PP1_MULTS)]) arrays, each holding npad doubles; the +1 entry
	// is for the (p+1)-style loop-multiplier scheme, to store fwdFFT(1), which it uses for carry-normalization:
   #ifdef USE_PP1_MULTS
	const int use_pp1 = 1;
   #else
	const int use_pp1 = 0;
   #endif
	nalloc = (m*num_b + use_pp1)*npad;
	j = 0;
	if(nalloc & 7)
		j = 8 - (nalloc & 7);
	nalloc += j;	ASSERT((nalloc & 7) == 0,"nalloc must be a multiple of 8!");	// Ensure 64-byte alignment of a[]
	// double*a holds ptr to 1 scratch vector, double**buf holds ptrs to num_b*m double-vecs of same length npad:
	a_ptmp = ALLOC_DOUBLE(a_ptmp, nalloc);
	if(!a_ptmp){
		sprintf(cbuf, "ERROR: unable to allocate the needed %u buffers of p-1 Stage 2 storage.\n",num_b*m + use_pp1);
		mlucas_fprint(cbuf,pm1_standlone+1);	ASSERT(0,cbuf);
	}
	a      = ALIGN_DOUBLE(a_ptmp);	ASSERT(((intptr_t)a & 63) == 0x0,"a[] not aligned on 64-byte boundary!");
	buf = (double **)calloc(num_b*m,sizeof(double *));
	// ...and num_b*m "buffers" for precomputed bigstep-coprime odd-square powers of the stage 1 residue:
	for(i = 0; i < num_b*m; i++) {
		buf[i] = a + i*npad;
//		fprintf(stderr,"buf[%3d] = %#" PRIX64 "\n",i,(uint64)buf[i]);
		ASSERT(((intptr_t)(buf[i]) & 63) == 0x0,"buf[i] not aligned on 64-byte boundary!");
	}
	// Still do fwdFFT(1) as init-FFT step in non-(p+1) build, but use uppermost buf[] entry to hold as throwaway result:
	vone = a + (i - 1 + use_pp1)*npad;
	a = mult[3];	// Rest of the way, a[] serves as a handy alias for mult[3]

   #ifdef MULTITHREAD
	/* Threadpool for parallel FFT(a-b) = FFT(a) - FFT(b) used in Stage 2 loop: Ideally, each thread
	processes npad/NTHREADS (a-b) double-diffs, thus npad should be divisible by NTHREADS. Further, the
	vec_double_sub inner loop - which gets run on each thread's data chunk - is set up such that it needs
	said data chunk to be aligned properly for the applicable SIMD vector type, i.e. on a multiple of
	8*RE_IM_STRIDE bytes, where RE_IM_STRIDE = [1,2,4,8] for [no-simd,sse2,avx,avx512]. As these dual
	constraints may be impossible to satisfy for a given runlength n, assign each thread a data chunk
	containing as close to npad/NTHREADS doubles as possible, but satisfying the aforementioned starting-
	address alignment criterion:
	*/
	thr_ret  = (int *)calloc(NTHREADS, sizeof(int));
	thread   = (pthread_t *)calloc(NTHREADS, sizeof(pthread_t));
	tdat     = (struct pm1_thread_data_t *)calloc(NTHREADS, sizeof(struct pm1_thread_data_t));
	// Initialize and set thread detached attribute:
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	const int nbytes_simd_align = (RE_IM_STRIDE*8) - 1;	// And per-thread data chunk addresses with this to check SIMD alignment
	ASSERT(((intptr_t)mult[0] & nbytes_simd_align) == 0x0,"mult[0] not aligned on 64-byte boundary!");
	ASSERT(((intptr_t)buf [0] & nbytes_simd_align) == 0x0,"buf [0] not aligned on 64-byte boundary!");	// Since npad a multiple of RE_IM_STRIDE, only need to check buf[0] alignment
	j = npad / NTHREADS;	// j = #doubles in each thread-processed chunk
	/* Fiddle up-or-downward to make it a multiple of RE_IM_STRIDE; say this == 8. Since j == (npad/NTHREADS) - [0 or 1]
	due to truncation-on-integer-div, if jmod := (j % RE_IM_STRIDE) < RE_IM_STRIDE/2, subtract jmod from j, otherwise
	add (RE_IM_STRIDE - jmod) to j: */
	k = j & (RE_IM_STRIDE - 1);	// k holds jmod
	if(k < (RE_IM_STRIDE>>1))
		j -= k;
	else
		j += (RE_IM_STRIDE - k);
	// Populate the parts of the thread-specific data structs which remain fixed from one call to the next:
	for(i = 0, k = 0; i < NTHREADS; ++i, k += j) {
		tdat[i].tid = i;
		tdat[i].retval = &thr_ret[i];
		tdat[i].arr0 =       a + k;	// a[]-array (alias for mult[3]) takes output, a[] = (mult[0][] - buf[i][])
		tdat[i].arr1 = mult[0] + k;	// k = i*j = doubles-offset for this thread's pair of array pointers
		tdat[i].arr2 = (double *)(intptr_t)k;	// For array-pointer 2, init the fixed offsets, then add fixed base-pointer offset buf[i] to
									// each k-index offset at thread-dispatch time, re-subtract buf[i] after pool work completion
		tdat[i].n = j;	// Chunksize
	}
	tdat[NTHREADS-1].n = npad - (NTHREADS-1)*j;	// Fiddle the last thread's chunksize so the sum == npad
	ASSERT(0x0 != (tpool = threadpool_init(NTHREADS, MAX_THREADS, NTHREADS, &thread_control)), "threadpool_init failed!");
	printf("%s: Init threadpool of %d threads\n",func,NTHREADS);
   #endif	// PM1_STANDALONE?
  #endif	// MULTITHREAD?
	// Integer arrays:
	b = malloc(m*(bigstep>>1)*sizeof(uint32));	ASSERT(b != NULL, "B[]-array alloc failed!");
	/* Jun 2021: added (psmall) map words for psmall = (mod 7|11) bitmap needed to support small-prime relocation -
	optimization - This needs wsize bytes, hence the (...+1)*wsize: */
	map = calloc((m+2+1)*wsize,sizeof(uint8));	ASSERT(map != NULL, "map[]-array alloc failed!");
	// 2 extra word-slots at high end of map used for these temps - can't declare as const pointers,0x but treat as such below:
	lo = map + m*wsize; hi = lo + wsize;
	rmap = hi + wsize;
#ifndef PM1_STANDALONE
	// We know that the arrtmp vector passed in by caller is large enough to hold
	// at least 2 nlimb-sized arrays,0x so simply point vec1 and vec2 at successive nlimb-sections of it:
	vec1 = arrtmp; vec2 = arrtmp + nlimb;
#endif

	/* Init the psmall = 7 relocation bitmap (cf. pm1.txt for derivation) - this is
	computed in terms of [psmall] 2x40-bit [lo40,hi40] pairs, thus psmall*10 = 70 bytes:
	*/
	const uint8 reloc_mod7_bytemap[70] = {
		0x82,0x04,0x02,0x08,0x40,0x02,0x10,0x40,0x20,0x41,
		0x00,0x0A,0x40,0x20,0x82,0x04,0x02,0x90,0x00,0x04,
		0x11,0x20,0x10,0x80,0x04,0x10,0x88,0x20,0x82,0x00,
		0x04,0x80,0x20,0x02,0x11,0x00,0x20,0x01,0x08,0x12,
		0x48,0x10,0x80,0x04,0x00,0x88,0x40,0x04,0x01,0x20,
		0x00,0x41,0x04,0x11,0x08,0x20,0x01,0x08,0x04,0x88,
		0x20,0x00,0x09,0x40,0x20,0x41,0x04,0x02,0x50,0x00
	};
	/* Init the psmall = 11 relocation bitmap (cf. pm1.txt for derivation) - this is
	computed in terms of [psmall] 2x24-bit [lo24,hi24] pairs, thus psmall*6 = 66 bytes:
	q0%11:	full 48-bit:	lo24	hi24	bytewise,0x lo8 to hi8:
		0	000002400000	400000	000002	0x00,0x00,0x40,0x02,0x00,0x00,
		1	000400840108	840108	000400	0x08,0x01,0x84,0x00,0x04,0x00,
		2	201080200022	200022	201080	0x22,0x00,0x20,0x80,0x10,0x20,
		3	000210001000	001000	000210	0x00,0x10,0x00,0x10,0x02,0x00,
		4	020040008401	008401	020040	0x01,0x84,0x00,0x40,0x00,0x02,
		5	084008000080	000080	084008	0x80,0x00,0x00,0x08,0x40,0x08,
		6	010000100210	100210	010000	0x10,0x02,0x10,0x00,0x00,0x01,
		7	802100020040	020040	802100	0x40,0x00,0x02,0x00,0x21,0x80,
		8	000800084000	084000	000800	0x00,0x40,0x08,0x00,0x08,0x00,
		9	440004010804	010804	440004	0x04,0x08,0x01,0x04,0x00,0x44,
		10	108021002000	002000	108021	0x00,0x20,0x00,0x21,0x80,0x10.
	*/
	const uint8 reloc_mod11_bytemap[66] = {
		0x00,0x00,0x40,0x02,0x00,0x00,
		0x08,0x01,0x84,0x00,0x04,0x00,
		0x22,0x00,0x20,0x80,0x10,0x20,
		0x00,0x10,0x00,0x10,0x02,0x00,
		0x01,0x84,0x00,0x40,0x00,0x02,
		0x80,0x00,0x00,0x08,0x40,0x08,
		0x10,0x02,0x10,0x00,0x00,0x01,
		0x40,0x00,0x02,0x00,0x21,0x80,
		0x00,0x40,0x08,0x00,0x08,0x00,
		0x04,0x08,0x01,0x04,0x00,0x44,
		0x00,0x20,0x00,0x21,0x80,0x10
	};
	const uint8*reloc_mod_psmall_bytemap = 0x0;	// Clang disallowed separate 'uint8*' type-decl in each if/else branch
	if(psmall == 7) {	// pinv64 = 64-bit Montgomery inverse of psmall: psmall*pinv64 == 1 (mod 2^64)
		reloc_mod_psmall_bytemap = reloc_mod7_bytemap;	pinv64 = 0x6DB6DB6DB6DB6DB7ull;
	} else if(psmall == 11) {
		reloc_mod_psmall_bytemap = reloc_mod11_bytemap;	pinv64 = 0x2E8BA2E8BA2E8BA3ull;
	}
	// For our product-of-small-primes D [a.k.a. bigstep] giant-step, find all coprime b in [1,M*D/2-1]:
	j = 0;
	/*fprintf(stderr,"(b[i+1]^2-b[i]^2)/24 = ");*/
	for(i = 1; i < m*(bigstep>>1); i++) {
		if(gcd32(bigstep,i) == 1) { b[j++] = i;	/*if(j > 1 && i < bigstep) fprintf(stderr,"%u,",((b[j-1]-b[j-2])*(b[j-1]+b[j-2]))/24);*/ }
	}
	if(j != m*num_b) {
		sprintf(cbuf,"Error: #ints coprime to D should == M*%u, instead got %u\n",num_b,j);
		mlucas_fprint(cbuf,pm1_standlone+1);
		retval = 1; goto ERR_RETURN;
	}
	for(j = 2*num_b; j < m*num_b; j++) {
		ASSERT(b[j] == bigstep + b[j-2*num_b],"Bigstep-power-offset check fails!");
	}

#if !defined(PM1_STANDALONE) && defined(PM1_DEBUG)
  // Test the 2-input modmul:
	// Do a throwaway init-mode modmul using a[] as a scratch array for the init section of func_mod_square():
	memcpy(mult[0],pow,nbytes);
	mode_flag = 0;	// Let's generate pure-int output
	ierr  = func_mod_square(mult[0], (void *)a, n, 0,1,(uint64)mode_flag, p, scrnFlag,&tdif2, FALSE, 0x0);
	// Actually, instead of doing a throwaway-to-init, Test the 2-input modmul by using it to effect a modsqr:
	memcpy(mult[1],pow,nbytes);	// 1 copy in mult[1]...
	memcpy(      a,pow,nbytes);	// another in a[]...
	ierr += func_mod_square(      a, 0x0, n, 0,1,     4ull + mode_flag, p, scrnFlag,&tdif2, FALSE, 0x0);/* fwdFFT(a) */
	ierr += func_mod_square(mult[1], 0x0, n, 0,1,(uint64)a + mode_flag, p, scrnFlag,&tdif2, FALSE, 0x0);/* and done. */
	if(ierr != 0) {
		sprintf(cbuf,"Modmul test hit an error of type = %u! Aborting.",ierr);
		ASSERT(0,cbuf);
	}
	convert_res_FP_bytewise(mult[0],(uint8*)vec1, n, p, 0x0,0x0,0x0);
	convert_res_FP_bytewise(mult[1],(uint8*)vec2, n, p, 0x0,0x0,0x0);
	ASSERT(mi64_cmp_eq(vec1,vec2,nlimb), "Modmul-test results mismatch!");
  /********************************************************************************/
  /********* Known-stage-2-factor tests, starting with a stage 1 residue: *********/
  /********************************************************************************/
  // F31: Do a single stage-1-result-powering (pow^140091319777 - 1) and make sure the known factor divides the result:
  if(p == 2147483648) {
	ASSERT(MODULUS_TYPE == MODULUS_TYPE_FERMAT, "This p-1 self-test requires Fermat-mod mode!");
	input_is_int = TRUE;
	memcpy(a,pow,nbytes);
	modpow(a, mult[0], input_is_int, 140091319777ull, func_mod_square, p, n, scrnFlag,&tdif2);
	ierr = func_mod_square(       a, 0x0, n, 0,1, 8ull, p, scrnFlag,&tdif2, FALSE, 0x0);	/* undo pass 1 of fwd-FFT */
	// subtract 1:
	a[0] -= 1;
	convert_res_FP_bytewise(a,(uint8*)vec1, n, p, 0x0,0x0,0x0);
	uint64 rem[2] = {0ull,0ull}, q[2] = {3118754346955702273ull,2544ull};	// k = 3.13.140091319777; q = k.2^(m+2) + 1
	// In fact, F31 has nlimb+1 words, but the only way a p-1 residue R has the same high bit
	// set as F31 iif R == F31 (uninteresting) or R == 2^2^31, which implies GCD == 1:
	int isfact = mi64_div(vec1,q, nlimb,2, 0x0, rem);
	ASSERT(isfact != 0, "Failed to find known stage 2 factor!");
	fprintf(stderr,"%s p-1 known-stage-2 prime stage 1 powering success!\n",PSTRING);
  }
  // M(139788679): Do a stage-1-result-powering (pow^a - 1) with a = 9952471 and make sure the corresponding
  // known factor, q = 842944537391616 = 2.k.p+1 with k = 2^9.3^2.11^2.29.37.1187^2, divides the result.
  // With B1 < 1187^2 = 1408969 this factor is not found after stage 1 since this prime appears only as a single-power:
  if(p == 139788679) {
	ASSERT(MODULUS_TYPE == MODULUS_TYPE_MERSENNE, "This p-1 self-test requires Mersenne-mod mode!");
	// A^4002923: Use mult[0] as scratch array for modpow():
	input_is_int = TRUE;
	memcpy(a,pow,nbytes);
	modpow(a, mult[0], input_is_int, 1187ull, func_mod_square, p, n, scrnFlag,&tdif2);
	ierr = func_mod_square(       a, 0x0, n, 0,1, 8ull, p, scrnFlag,&tdif2, FALSE, 0x0);	// undo pass 1 of fwd-FFT
	// subtract 1:
	a[0] -= 1;
	convert_res_FP_bytewise(a,(uint8*)vec1, n, p, 0x0,0x0,0x0);
	uint64 rem[2] = {0ull,0ull}, q[2] = {11051162840690736129ull,12775ull};	// q = 1314651028704963254300497
	int isfact = mi64_div(vec1,q, nlimb,2, 0x0, rem);
	ASSERT(isfact != 0, "Failed to find known stage 2 factor!");
	fprintf(stderr,"%s p-1 known-stage-2 prime self-test success!\n",PSTRING);
	exit(0);
  }
#endif
	/* Here are the first 2*num_b entries for our 2 supported bigstep values - each sequence starts with 1, a copy of the stage 1
	residue, followed by ascending odd powers of same. Thus we also give the associated 2*incr sequences needed for an incremental-add loop:
	210: b[i] = 1,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,121,127,131,137,139,143,149,151,157,163,167,169,173,179,181,187,191,193,197,199,209
		*2-diffs: -,5,1,2,1,2,3,1,3,2,1,2,3,3,1,3,2,1,3,2,3,4,2,1,2,1,2,4,3,2,3,1,2,3,1,3,3,2,1,2,3,1,3,2,1,2,1,5
	330: b[i] = 1,7,13,17,19,23,29,31,37,41,43,47,49,53,59,61,67,71,73,79,83,89,91,97,101,103,107,109,113,119,127,131,133,137,139,149,151,157,161,163,167,169,173,179,181,191,193,197,199,203,211,217,221,223,227,229,233,239,241,247,251,257,259,263,269,271,277,281,283,287,289,293,299,301,307,311,313,317,323,32
		*2-diffs: -,3,3,2,1,2,3,1,3,2,1,2,1,2,3,1,3,2,1,3,2,3,1,3,2,1,2,1,2,3,4,2,1,2,1,5,1,3,2,1,2,1,2,3,1,5,1,2,1,2,4,3,2,1,2,1,2,3,1,3,2,3,1,2,3,1,3,2,1,2,1,2,3,1,3,2,1,2,3,3
	Once we have the first num_b powers, the next num_b satisfy the mirror-symmetry b[i] = bigstep - b[2*num_b-1-i],
		e.g. for 210, num_b = 24 and b[24] = 107 = 210-103 = bigstep - b[23].
	Once we have the first 2*num_b powers, each remaining one follows the simple formula b[i] = bigstep + b[i-2*num_b]
	Problem:
		We actually need the stage 1 residue (A) raised to the above sequence values, SQUARED. Thus we need the diffs between
	successive squares for our up-multiply loop. Here are those diffs for the first 2*num_b entries for bigstep value 210, for example:
	210: square-diffs, all are multiples of 24. For each distinct multiple of 24, @ means first appearance in the sequence,
			* means last; entries lacking these annotators appear just once in this low-2*num_b-terms square-diffs sequence:
	b[i] =     1,11,13,17,19 ,23,29 ,31,37 ,41 ,43,47,53 ,59 ,61,67 ,71 ,73,79,83,89 ,97 ,101,103,107,109,113,121,127,131,137,139,143,149,151,157,163,167,169,173,179,181,187,191,193,197,199,209
	s[i] = 24*[-,5@,2@*,5,3@*,7@,13@,5*,17@,13*,7*,15,25@,28@,10,32@,23@,12,38,27,43@,62@,33@,17*,35 ,18 ,37 ,78 ,62*,43*,67 ,23*,47 ,73 ,25*,77 ,80 ,55 ,28*,57 ,88 ,30 ,92 ,63 ,32*,65 ,33*,170]
	Simpler to use that diffs between successive odd squares follow a simple +=8 rule: 3^2-1^2 = 8, 5^2-3^2 = 16, 7^2-5^2 = 24, etc.
	Set 0-term of our A^(b[i]^2) sequence = A^1, fwd-FFT that. 3 modsquares of a 2nd copy of A^1 give A^8, which we also
	fwd-FFT and use as a fixed multiplier, with which we will generate ascending multiple-of-8 powers, A8,16,24,... . We then use
	these to run through ascending odd-square powers, and each time our loop index curr_i = 1,3,5,... hits the next entry in our b[i]
	sequence, we deposit the associated A^(b[i]^2) term in the corresponding buf[i] entry:
	i	j	mult[0]	mult[1]	b[i] ?=j	Action
	0	1	A^1 	A^8		1	1	buf[i++] = mult[0]
	1	3	A^9		A^16	11	0
	1	5	A^25	A^24	11	0
	1	7	A^49	A^32	11	0
	1	9	A^81	A^40	11	0
	1	11	A^121	A^48	11	1	buf[i++] = mult[0], etc.
	*/
#ifndef PM1_STANDALONE
  #ifdef CTIME
	clock1 = clock();
  #else
	clock1 = getRealTime();
  #endif
	nerr = 0;
	MME = 0.0;	// Reset (extern) maxROE
	/* Init a const-vector fwdFFT(1) to and store in highest-index buf[] slot; use for normalization-multiply
	prior to subtract-step in (p+1)-style Lucas-sequence multiplier computation. Since this may be the first modmul we do,
	use a[] as a scratch array for the init section of func_mod_square(). No residue-shift here, so rightmost arg = False.
	Alternatively, could use a simple loop-based init to set all the Re() entries of vone[] to 1, but that needs more code.
	Still compute fwdFFT(1) by way of an init-FFT in non-(p+1) build, but map vone[] to uppermost buf[] entry, i.e. throw away:
	*/
	vone[0] = 1.0; for(j = 1; j < n; j++) { vone[j] = 0.0; }
													//vvv-- Pure-int inputs, so mode_flag = 0
	ierr = func_mod_square(vone, (void *)a, n, 0,1, 4ull, p, scrnFlag,&tdif2, FALSE, 0x0);	ASSERT(ierr == 0,"fwdFFT(1) hit error!");

  #if !USE_PP1_MULTS		// Basic version:

	sprintf(cbuf,"Using Bigstep %u, pairing-window multiplicity M = %u: Init M*%u = %u [base^(A^(b^2)) %% n] buffers for Stage 2...\n",bigstep,m,num_b,m*num_b);
	mlucas_fprint(cbuf,pm1_standlone+1);
	// [a] Generate set of precomputed buffers A^(b^2) (A = s1 residue stored in pow[]) for b-values corr. to our choice of D:
	memcpy(buf[0] ,pow,nbytes);	// b[0] = 1 --> Copy of A^1 into buf[0]
	memcpy(mult[0],pow,nbytes);	// Another copy of A^1 into mult[0][] - this will hold ascending odd-square powers A^1,9,25,...
	memcpy(mult[1],pow,nbytes);	// A third copy of A^1 into mult[1][] - this will end up holding A^8 in fwd-FFTed form:
	mode_flag = 2;	// bit 1 of mode_flag = 1 since all FFT-mul outputs will be getting re-used as inputs
	// 3 mod-squares in-place to get A^8:
	ierr = func_mod_square(mult[1], 0x0, n, 0,3,        (uint64)mode_flag, p, scrnFlag,&tdif2, FALSE, 0x0); nerr += ierr;
	if(ierr == ERR_INTERRUPT) {
		return ierr;
	}
	mode_flag |= 1;	// Only 1st fft-mul of A and its copies need low bit of mode_flag = 0!
	memcpy(mult[2] ,mult[1],nbytes);	// mult[2][] holds ascending A^8,16,24,..., in *non*-fwd-FFTed form (that is, fwd-FFT-pass-1-done form)
	// mult[1] = fwdFFT(A^8):					  vvvv + 4 to effect "Do in-place forward FFT only" of pow^8:
	ierr = func_mod_square(mult[1], 0x0, n, 0,1, 4ull + (uint64)mode_flag, p, scrnFlag,&tdif2, FALSE, 0x0); nerr += ierr;
	// And now the big loop to compute the remaining pow^(b[i]^2) terms, each of which goes into buf[i++] in fwd-FFTed form.
	// Each pass through the loop costs 2.5 modsqr and 1 memcpy, with a 2nd memcpy whenever we hit a b[i] and write a buf[] entry:
	i = 1;
   #ifdef PM1_DEBUG
	fprintf(stderr,"Init buf[] = A^");
   #endif
	for(j = 3; j < m*(bigstep>>1); j += 2) {
		memcpy(a,mult[2],nbytes);	// a[] = Copy of fwd-FFT-pass-1-done(A^8,16,24,...) to be fwd-FFTed
		// a[] = FFT(A^8,16,24,...):
		ierr = func_mod_square(      a, 0x0, n, 0,1,         4ull + (uint64)mode_flag, p, scrnFlag,&tdif2, FALSE, 0x0); nerr += ierr;
		// mult[0] *= a[]: mult[0] holds result, thus is not fwd-FFTed (that is, in fwd-FFT-pass-1-done form) on entry;
		// Since mult[0] holds pure-int copy of stage 1 residue A on loop entry, bit 0 of mode_flag = 0 for just its first use:
		//                                                                                vvvvvvvv
		ierr = func_mod_square(mult[0], 0x0, n, 0,1, (uint64)a + (uint64)(mode_flag - (j==3)), p, scrnFlag,&tdif2, FALSE, 0x0); nerr += ierr;
		if(ierr == ERR_INTERRUPT) {
			return ierr;
		}
		if(b[i] == j) {
		#ifdef PM1_DEBUG
			fprintf(stderr,"%u^2.",j);
		#endif
//			fprintf(stderr,"buf[%3d] = %#" PRIX64 "\n",i,(uint64)buf[i]);
			ASSERT(((intptr_t)(buf[i]) & 63) == 0x0,"buf[i] not aligned on 64-byte boundary!");
			memcpy(buf[i++],mult[0],nbytes);	// buf[i++] = mult[0] = fwd-FFT-pass-1-done(A^1,9,25,...)
		}
		// Up-multiply the fwd-FFT-pass-1-done(A^8,16,24,...) by fixed multiplier fwd-FFT(A^8):
		// mult[2] = A^16,24,... :
		ierr = func_mod_square(mult[2], 0x0, n, 0,1, (uint64)mult[1] + (uint64)mode_flag, p, scrnFlag,&tdif2, FALSE, 0x0); nerr += ierr;
		if(ierr == ERR_INTERRUPT) {
			return ierr;
		}
	}

  #else	// USE_PP1_MULTS = True: (p+1)-style version of stage 2 loop-multipliers, needs just 1 modmul per bigstep-loop
	#warning Add signal-handling for non-fwd-only-FFTs!
	// vec1,vec2 will hold stage 1 residue A (stored in FP form in pow[]) and its mod-inverse in packed-bit form
	vec1[nlimb-1] = 0ull;
	convert_res_FP_bytewise(pow,(uint8*)vec1, n, p, &Res64,&Res35m1,&Res36m1);
fprintf(stderr,"#1: vec1 = A^+1 checksums = %" PRIu64 ",%" PRIu64 ",%" PRIu64 "; FP(A)[0:1] = %10.2f,%10.2f\n",Res64,Res35m1,Res36m1, pow[0],pow[1]);
	// First see if there's a savefile-copy of the s1 residue-inverse:
	strcpy(inv_file, RESTARTFILE);
	inv_file[0] = ((MODULUS_TYPE == MODULUS_TYPE_MERSENNE) ? 'p' : 'f');
	strcat(inv_file, ".inv");
	fp = mlucas_fopen(inv_file,"rb");
	// o If exists, read nsquares field (expected == 0 for s1-inv) into uint64 tmp:
	if(fp) {												// G-check residue fields all set NULL in this call:
		i = read_ppm1_savefiles(inv_file, p, &k, fp, &tmp, (uint8*)vec2, &Res64,&Res35m1,&Res36m1, 0x0,0x0,0x0,0x0);
		fclose(fp); fp = 0x0;
		ASSERT(tmp == 0ull, "Stage 1 residue-inverse savefile should have nsquares == 0!");
		if(!i) {
			/* First print any error message that may have been issued during the above function call: */
			if(strstr(cbuf, "read_ppm1_savefiles"))
				mlucas_fprint(cbuf,pm1_standlone+1);
			// And now for the official spokesmessage:
			snprintf(cbuf,STR_MAX_LEN*2, "Read of stage 1 residue-inverse savefile %s failed for reasons unknown. Computing inverse...\n",inv_file);
			mlucas_fprint(cbuf,pm1_standlone+1);
		} else {
			s1_inverse = TRUE;
		}
	}
	if(!s1_inverse) {
		snprintf(cbuf,STR_MAX_LEN*2, "Stage 2: Computing mod-inverse of Stage 1 residue...\n");	mlucas_fprint(cbuf,pm1_standlone+1);
		modinv(p,vec1,vec2,nlimb);	// Result in vec2
		Res64 = vec2[0];
		Res35m1 = mi64_div_by_scalar64(vec2,two35m1,nlimb,0x0);
		Res36m1 = mi64_div_by_scalar64(vec2,two36m1,nlimb,0x0);
		// Write inverse to savefile:
		fp = mlucas_fopen(inv_file, "wb");
		if(fp) {
			write_ppm1_savefiles(inv_file,p,n,fp, 0ull, (uint8*)vec2,Res64,Res35m1,Res36m1, 0x0,0x0,0x0,0x0);
			fclose(fp);	fp = 0x0;
		} else {
			snprintf(cbuf,STR_MAX_LEN*2, "ERROR: unable to open restart file %s for write of checkpoint data.\n",inv_file);
			mlucas_fprint(cbuf,pm1_standlone+1);	ASSERT(0,cbuf);
		}
	}
	convert_res_bytewise_FP((uint8*)vec2, a, n, p);	// Use a[] to hold inverse A^-1 until done with it
fprintf(stderr,"#1: vec2 = A^-1 checksums = %" PRIu64 ",%" PRIu64 ",%" PRIu64 "; FP(A^-1)[0:1] = %10.2f,%10.2f\n",Res64,Res35m1,Res36m1, a[0],a[1]);
   #ifdef PM1_DEBUG
	fprintf(stderr,"Checking mod-inverse...\n");
	// Debug: check inverse ... start with copies of A (pow[]) and A^-1 (a[]) into mult0-1:
	memcpy(mult[0],pow,nbytes);	memcpy(mult[1],a,nbytes);
	// mult[1] = fwdFFT(A^-1):					       vvvv + 4 to effect "Do in-place forward FFT only"; bit 0 of mode_flag (the only one fwd-FFT cares about) = 0, since mult[0] enters in pure-int form
	ierr += func_mod_square(mult[1],      0x0, n, 0,1, 4ull, p, scrnFlag,&tdif2, FALSE, 0x0);
	// mult[0] = A * A^-1, check that result = 1 as expected:
	mode_flag = 0;	// bits 0:1 of mode_flag = 0, since mult[0] enters in pure-int form and want output the same way
	ierr += func_mod_square(mult[0], 0x0, n, 0,1, (uint64)mult[1] + (uint64)mode_flag, p, scrnFlag,&tdif2, FALSE, 0x0);
	ASSERT(mult[0][0] == 1.0, "inverse-check fails!");
	for(i = 1; i < npad; i++) { ASSERT(mult[0][i] == 0.0, "inverse-check fails!"); }
   #endif
	/*** NOTE: vec1,vec2 hold stage 1 residue A and its mod-inverse in packed-bit form, can use any of
	pow[], mult[0-2][], a[] in our ensuing stage 2 inits and still re-obtain b and b^-1 from vec1,2 at any time. ***/
	/*
	In order to use the Lucas-sequence index-incrementing identity V[n+j] = V[n]*V[j] - V[n-j] to get e.g.
	V[3] = (A^(3*b) + A^-(3*b)) from (A^b + A^-b), first need V[2] = V[1]^2 - 2 (hello, Lucas-Lehmer!),
	then V[3] = V[2]*V[1] - V[1] = V[1]*(V[2] - 1). Subsequent odd-index Vs then can be computed as
	V[2*i+1] = V[2*i-1]*V[2] - V[2*i-3], e.g. 	V[5] = V[3]*V[2] - V[1] .
	*/
	sprintf(cbuf,"Using Bigstep %u, pairing-window multiplicity M = %u: Init M*%u = %u [base^(A^b + A^-b) %% n] buffers for Stage 2...\n",bigstep,m,num_b,m*num_b);
	mlucas_fprint(cbuf,pm1_standlone+1);
	// [a] Generate a set of precomputed buffers (A^b + A^-b) for the set of b-values corr. to our choice of D:
   #ifdef PM1_DEBUG
	fprintf(stderr,"Init buf[] = (A^j + A^-j) with j = 1.");
   #endif
	// Need to renormalize A + A^-1 sum prior to modmul - this is easier to do in packed-bit form, vec2 = vec1+vec2:
   #if 0
	tmp = mi64_add(vec1,vec2, vec2,nlimb-(MODULUS_TYPE == MODULUS_TYPE_FERMAT));	// Note, for Fermat-mod use 1 more
			// limb than p>>6, so shorten sub-vec-length by 1 to make sure any carry ends up in tmp, not vec2[nlimb-1]
	// Sign of wraparound carry depends on modulus type:
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE) {
		ASSERT(tmp == 0ull,"Mersenne-mod vec1+vec2 should never spill over into next word!");
		// Now get the really carry bit, bit p. Normalized residue only uses bits <0:p-1>:
		bit = p&63; word = p>>6; tmp = vec2[word]>>bit;
		q = mi64_add_scalar(vec2, tmp, vec2, nlimb);
		ASSERT(q == 0ull,"Mersenne-mod vec1+vec2 wraparound carry should never have carry-out!");
	} else {
		q = mi64_sub_scalar(vec2, tmp, vec2, nlimb-1);	// Again shorten vec-sub-length by 1 to properly check for borrow-out
		ASSERT(q == 0ull,"Fermat-mod vec1+vec2 wraparound carry should never have borrow-out!");
	}
	convert_res_bytewise_FP((uint8*)vec2,buf[0], n, p);	// buf[0] = V[1]
	// Now recover original vec2 = A^-1 from the FP version in a[]:
	convert_res_FP_bytewise(a,(uint8*)vec2, n, p, &Res64, &Res35m1, &Res36m1);
fprintf(stderr,"#2: vec2 = A^-1 checksums = %" PRIu64 ",%" PRIu64 ",%" PRIu64 "\n",Res64,Res35m1,Res36m1);
   #else
	// Original FP code for V[1] lacks post-add normalization:
	for(i = 0; i < npad; i++) { buf[0][i] = pow[i] + a[i]; }	// V[1] = A^1 + A^-1
   #endif
dtmp = 0; for(k = 0; k < npad; k++) { dtmp += fabs(buf[0][k]); } dtmp /= n;
fprintf(stderr,"V[1] has L1 = %20.10f\n",dtmp);

	memcpy(mult[0], buf[0],nbytes);	// Copy of V[1] into mult[0]
	// Put mult[0] = V[1] into fwd-FFT-pass-1-done form:
	ierr += func_mod_square(mult[0], 0x0, n, 0,1, -4ull, p, scrnFlag,&tdif2, FALSE, 0x0);
	memcpy(mult[1],mult[0],nbytes);	// Copy of fwd-FFT-pass-1-done(V[1]) into mult[1]
	memcpy(mult[2],mult[0],nbytes);	// Copy of fwd-FFT-pass-1-done(V[1]) into mult[2]
	mode_flag = 3;	// bit 1 of mode_flag = 1 since all FFT-mul outputs will be getting re-used as inputs
					// bit 0 also = 1 since mult[0-2] holds copies of stage 1 residue A in fwd-FFT-pass-1-done form.
	// 1 mod-square in-place and -2 to get V[2] - No residue-shift here, so rightmost arg = False:
	ierr += func_mod_square(mult[1], 0x0, n, 0,1, (uint64)mode_flag, p, scrnFlag,&tdif2, FALSE, 0x0);
	// Need to undo pass 1 of fwd-FFT prior to -= 2; do this just as with fwd-FFT-only, but with flag = 8 instead of 4:
	ierr += func_mod_square(mult[1], 0x0, n, 0,1, 8ull, p, scrnFlag,&tdif2, FALSE, 0x0);
	mult[1][0] -= 2;	// mult[1] = V[2]
#if 0	/******************************/
// Debug: Also compute V[2] via separate mod-squarings of A and A^-1, followed by summation and normalization-modmul:
// Have pure-int inputs and want fwd-FFT-pass-1-done outputs, so force 2 in place of mode_flag:
ierr  = func_mod_square(pow, 0x0, n, 0,1, 2ull, p, scrnFlag,&tdif2, FALSE, 0x0);	// A^+2
ierr  = func_mod_square(  a, 0x0, n, 0,1, 2ull, p, scrnFlag,&tdif2, FALSE, 0x0);	// A^-2
for(k = 0; k < npad; k++) { a[k] += pow[k]; }	// V[2] = A^2 + A^-2
// Normalization-modmul - Have fwd-FFT-pass-1-done inputs and want pure-int outputs, so force 1 in place of mode_flag:
ierr += func_mod_square(  a, 0x0, n, 0,1, (uint64)vone + 1ull, p, scrnFlag,&tdif2, FALSE, 0x0);
j = 0; for(k = 0; k < npad; k++) {	// j stores #mismatches
	if(a[k] != mult[1][k]) {
		fprintf(stderr,"V[2] check: a[%u] (%10.2f) != mult[1][%u] (%10.2f)\n",j,a[k],j,mult[1][k]);
		j++;
	}
	if(j > 100) break;
}
exit(0);
#endif	/******************************/
	// Ensuing do-full-fwd-FFT on V[2] again has pure-int inputs, so use 0 in place of mode_flag:
	// mult[1] = fwdFFT(V[2]):				      vvvv + 4 to effect "Do in-place forward FFT only":
	ierr += func_mod_square(mult[1], 0x0, n, 0,1, 4ull, p, scrnFlag,&tdif2, FALSE, 0x0);
	// V[3] = V[1]*(V[2] - 1) ... we compute as V[2]*V[1] - V[1] since it makes the logical sequencing less awkward:
MME = 0;
#warning *************************** MULSUB path here gives ROE = 0.5! ***************************
  #if 1//def USE_VEC_DBL_SUB
	// mult[0] = V[2]*V[1]:
	ierr += func_mod_square(mult[0], 0x0, n, 0,1, (uint64)mult[1] + (uint64)mode_flag, p, scrnFlag,&tdif2, FALSE, 0x0);
	// V[3] = V[2]*V[1] - V[1]; both sub-inputs are fwd-FFT-pass-1-done. The parallel-vector-sub leaves result in a[]:
   #ifdef MULTITHREAD
	vec_double_sub(tpool,tdat,mult[2]);
   #else
	vec_double_sub(mult[0],mult[2],a,npad);
   #endif
	memcpy(mult[0],a,nbytes);
  #else
	// V[3] = V[2]*V[1] - V[1]:
	ierr += func_mod_square(mult[0], 0x0, n, 0,1, (uint64)mult[1] + (uint64)mode_flag, p, scrnFlag,&tdif2, FALSE, (double*)(~((uint64)mult[2])));
/*	ierr += func_mod_square(mult[0], 0x0, n, 0,1, (uint64)vone + (uint64)mode_flag, p, scrnFlag,&tdif2, FALSE, 0x0);
   #ifdef MULTITHREAD
	vec_double_sub(tpool,tdat,mult[2]);
   #else
	vec_double_sub(mult[0],mult[2],a,npad);
   #endif
	memcpy(mult[0],a,nbytes);
*/
  #endif
dtmp = 0; for(k = 0; k < npad; k++) { dtmp += fabs(mult[0][k]); } dtmp /= n;
fprintf(stderr,"V[3] has L1 = %20.10f, gives IERR = %u, ROE = %10.6f\n",dtmp,ierr,MME);
	/* At this point:
		buf[0]  holds V[1] in pure-int form; mult[2] holds V[1] in fwd-FFT-pass-1-done form
		mult[0] holds V[3] in fwd-FFT-pass-1-done form
		mult[1] = fwdFFT(V[2])
	*/
	i = 1;
	if(b[i] == 3) {
	#ifdef PM1_DEBUG
		fprintf(stderr,"%u.",j);
	#endif
		memcpy(buf[i++],mult[0],nbytes);	// buf[i++] = mult[0] = (A^j + A^-j) in fwd-FFT-pass-1-done form
	}

	/* Loop to compute the remaining (A^b[i] + A^-b[i]) terms, each of which goes into buf[i++] in fwd-FFT-pass-1-done form.
	Note that odd-index V[] then can be computed as V[2*i+1] = V[2*i-1]*V[2] - V[2*i-3], e.g. 	V[5] = V[3]*V[2] - V[1] .
	On entry, mult[0-2] hold V[3,2,1], with the odd-index V[1,3]-terms in fwd-FFT-pass-1-done form and V[2] in fwd-FFTed form:
	Each loop costs 1 mulsub (modsqr & vector-subtract) and 3 memcpy, with a 3rd memcpy whenever we hit a b[i] and write a buf[] entry:
	*/
	jhi = m*(bigstep>>1);
	// Allow equality at upper end here so mult[0] = V[D/2] on loop exit, can use that to compute V[D]
	// and compare the latter to separate-path way of getting V[D] in the loop-multiplier-computation section:
	for(j = 5; j <= jhi; j += 2) {
		if(j < jhi) memcpy(buf[i],mult[0],nbytes);	// Stash a copy of V[2*i-1] in buf[i]
MME = 0;
	#if 1//def USE_VEC_DBL_SUB
	  /*** The problem with this build mode here is that while the V[odd]*V[2] intermediate
	  is normalized, we follow with an unrenormalized subtract, leading to the norm of the subtrahend
	  increasing monotonically: e.g. for p = 108268067, n = 6M have normalized-balanced-digit average
	  of 17.209 bpd, corresponding to a fractional average-base of ~151478, which means a balanced-digit
	  normalization should have L1 half that, ~75739. (But that crude average does not take into account
	  #bw-vs-sw: bw = p%n = 1313315 ~= 20% bigwords, thus L1 will be < 75739.) L1 data below based on FFT-pass-1-done outputs:
		V[2]*V[1]             has L1 =  705912.1840996306, ROE = 0.054688
	  Up to this point things are fine - V[1]^2 is properly normalized, and the unnormalized -= 2 does not signify.
	  But now we start using unnormalized whole-vector subtracts, e.g. V[3] = V[2]*V[1] - V[1]. So V[3] L1 jumps:
		V[3]                  has L1 = 1222716.6898910117, ROE = 0.054688
		V[3]*V[2]             has L1 =  705994.3160107044, ROE = 0.078125
	  Similarly for rest - the V[odd]*V[2] intermediate is normalized, but L1 of the subtrahend increases monotonically:
		V[ 5] = V[ 3]*V[2] - V[ 1] has L1 = 1222786.9329456112	V[ 5]*V[2] has L1 =  705730.2279694607, ROE = 0.070312
		V[ 7] = V[ 5]*V[2] - V[ 3] has L1 = 1411983.6492863866	V[ 7]*V[2] has L1 =  705652.3466845521, ROE = 0.085938
		V[ 9] = V[ 7]*V[2] - V[ 5] has L1 = 1411222.1661252610	V[ 9]*V[2] has L1 =  705793.4601404053, ROE = 0.078125
		V[11] = V[ 9]*V[2] - V[ 7] has L1 = 1578687.5101256447	V[11]*V[2] has L1 =  705779.8257252158, ROE = 0.093750
															...
		V[101]= V[99]*V[2] - V[97] has L1 = 3669572.1187865455	V[101]*V[2] has L1 = 705953.0830486178, ROE = 0.218750
		V[103]=V[101]*V[2] - V[99] has L1 =   3735250.9816183313
	  ***/
		// V[2*i-1] *= V[2]:
		ierr += func_mod_square(mult[0], 0x0, n, 0,1, (uint64)mult[1] + (uint64)mode_flag, p, scrnFlag,&tdif2, FALSE, 0x0);
		// V[2*i+1] = V[2*i-1]*V[2] - V[2*i-3]; the parallel-vector-sub leaves result in a[]:
	  #ifdef MULTITHREAD
		vec_double_sub(tpool,tdat,mult[2]);
	  #else
		vec_double_sub(mult[0],mult[2],a,npad);
	  #endif
		// Normalization-multiply by fwdFFT(1):
		ierr += func_mod_square(a, 0x0, n, 0,1, (uint64)vone + (uint64)mode_flag, p, scrnFlag,&tdif2, FALSE, 0x0);
		memcpy(mult[0],     a,nbytes);	// V[2*i+1] replaces V[2*i-1]:
		if(j < jhi) memcpy(mult[2],buf[i],nbytes);	// V[2*i-1] replaces V[2*i-3]
	#else
	  #error This path gives ROE = 0.5!
		// V[2*i+1] replaces V[2*i-1]:
		ierr += func_mod_square(mult[0], 0x0, n, 0,1, (uint64)mult[1] + (uint64)mode_flag, p, scrnFlag,&tdif2, FALSE, (double*)(~((uint64)mult[2])));
		memcpy(mult[2],buf[i],nbytes);	// V[2*i-1] replaces V[2*i-3]
		// MULSUB version needs subtrahend fully-fwdFFTed:
		ierr += func_mod_square(mult[2], 0x0, n, 0,1, 4ull + (uint64)mode_flag, p, scrnFlag,&tdif2, FALSE, 0x0);
	#endif
//	dtmp = 0; for(k = 0; k < npad; k++) { dtmp += fabs(mult[0][k]); } dtmp /= n;
//	fprintf(stderr,"V[%u] = V[%u]*V[2] - V[%u] has L1 = %20.10f; IERR = %u, ROE = %10.6f\n",j,j-2,j-4,dtmp,ierr,MME);
		if(b[i] == j) {
		#ifdef PM1_DEBUG
			fprintf(stderr,"%u.",j);
		#endif
			memcpy(buf[i++],mult[0],nbytes);	// buf[i++] = mult[0] = (A^j + A^-j) in fwd-FFT-pass-1-done form
		}
	}

   #if 0
	#define pp1_vd2_cross_check
	// PP1 Debug: copy V[D/2] into vector :
	memcpy(a,mult[0],nbytes);
	dtmp = 0;
	for(k = 0; k < npad; k++) { dtmp += fabs(a[k]); }
	dtmp /= n;
	fprintf(stderr,"Autosquare: V[D/2] with L1 = %20.10f",dtmp);
	// 1 mod-square in-place -2 to get V[D]:
	ierr += func_mod_square(a, 0x0, n, 0,1, (uint64)mode_flag, p, scrnFlag,&tdif2, FALSE, 0x0);
	fprintf(stderr," gives IERR = %u, ROE = %10.6f\n",ierr,MME);
	// Need to undo pass 1 of fwd-FFT prior to -= 2; do this just as with fwd-FFT-only, but with flag = 8 instead of 4:
	ierr += func_mod_square(a, 0x0, n, 0,1,              8ull, p, scrnFlag,&tdif2, FALSE, 0x0);
	a[0] -= 2;	// a = V[D]
	convert_res_FP_bytewise(a,(uint8*)vec1, n, p, &Res64, &Res35m1, &Res36m1);
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE) {
		// Get the carry bit, bit p. Normalized residue only uses bits <0:p-1>:
		bit = p&63; word = p>>6; tmp = vec1[word]>>bit;
	} else {
		tmp = vec1[nlimb-1];
	}
	ASSERT(tmp == 0ull,"Properly normalized residue should never spill over into next word!");
   #endif	// PM1_DEBUG?

  #endif	// USE_PP1_MULTS?

  #ifdef PM1_DEBUG
	fprintf(stderr,"\n");
  #endif
	if(nerr != 0) {
		sprintf(cbuf,"Stage 2 buffer-init hit 1 or more fatal errors! Aborting.");
		mlucas_fprint(cbuf,pm1_standlone+0);	ASSERT(0,cbuf);
	}
	if(i != m*num_b) {
		sprintf(cbuf,"Stage 2: Incorrect loop-exit value of buffer-index!");
		mlucas_fprint(cbuf,pm1_standlone+0);	ASSERT(0,cbuf);
	}
	// buf[] entries all need to be rest-of-fwd-FFTed;
	for(i = 0; i < m*num_b; i++) {
		// Since buf[0] holds pure-int copy of stage 1 residue A on loop entry, bit 0 of mode_flag = 0 for just it:
		//                                                                    vvvvvvvv
		ierr = func_mod_square(buf[i], 0x0, n, 0,1, 4ull + (uint64)(mode_flag - (i==0)), p, scrnFlag,&tdif2, FALSE, 0x0); nerr += ierr;
	}	ASSERT(nerr == 0, "fwdFFT of buf[] entries returns error!");

	// Accumulate the cycle count in a floating double on each pass to avoid problems
	// with integer overflow of the clock() result, if clock_t happens to be 32-bit int on the host platform:
  #ifdef CTIME
	clock2 = clock();
  #else
	clock2 = getRealTime();
  #endif
	*tdiff = clock2 - clock1; clock1 = clock2;
	snprintf(cbuf,STR_MAX_LEN*2, "Buffer-init done; clocks =%s, MaxErr = %10.9f.\n",get_time_str(*tdiff), MME);
	mlucas_fprint(cbuf,pm1_standlone+1);

	/********************* RESTART FILE STUFF: **********************/
	// Make sure we append '.s2' to the primary restart file name:
	strcpy(savefile, RESTARTFILE);
	savefile[0] = ((MODULUS_TYPE == MODULUS_TYPE_MERSENNE) ? 'p' : 'f');
	strcat(savefile, ".s2");
	// [From the above p-1 savefile schema] 3. On entry, S2 checks for existence of ".s2" savefile:
	fp = mlucas_fopen(savefile,"r");
	// o If exists, read nsquares field into uint64 qlo, mask off high byte (which stores the value of any relocation-prime
	// psmall used for stage 2), compare vs original-assignment B2_start read (or inferred, as B2_start = B1) from worktodo entry:
	if(fp) {												// G-check residue fields all set NULL in this call:
		i = read_ppm1_savefiles(savefile, p, &k, fp, &qlo, (uint8*)arrtmp, &Res64,&Res35m1,&Res36m1, 0x0,0x0,0x0,0x0);
		fclose(fp); fp = 0x0;
		if(i && psmall)	{
			// We expect the main-program S2-invocation code to have resolved this kind of mismatch via bigstep selection:
			if(qlo >> 56)
				ASSERT((uint32)(qlo >> 56) == psmall, "Mismatch between relocation-prime set for stage 2 restart and one read from S2 savefile!");
		}
		qlo &= 0x00ffffffffffffffull;	// Mask off high byte storing psmall
		// If savefile-read fails, start stage 2 from B2_start:
		if(!i) {
			qlo = B2_start;
		} else {
			if(k != kblocks) {
				sprintf(cbuf,"INFO: %s savefile has kblocks = (n>>10) = %u; Instead using kblocks %u for stage 2.\n",func,k,kblocks);
				mlucas_fprint(cbuf,pm1_standlone+1);
			}
			// If nsquares > B2_start, arrtmp holds the S2 interim residue for q = nsquares; set up to restart S2 at that point.
			if(qlo >= B2_start) {
				snprintf(cbuf,STR_MAX_LEN*2, "Read stage 2 savefile %s ... restarting stage 2 from q = %" PRIu64 ".\n",savefile,qlo);
			} else {	// If user running a new partial S2 interval with bounds larger than a previous S2 run, allow but info-print to that effect:
				snprintf(cbuf,STR_MAX_LEN*2, "INFO: %s savefile has qlo[%" PRIu64 "] <= B2_start[%" PRIu64 "] ... Stage 2 interval will skip intervening primes.\n",func,qlo,B2_start);
			}
			mlucas_fprint(cbuf,pm1_standlone+1);
			restart = TRUE;
			if(qlo >= B2)	// qlo >= S2 upper limit - nothing to do but proceed to gcd
				goto S2_RETURN;
		}
	}

#endif	// #ifndef PM1_STANDALONE

	/***** Small-prime-relocation: for the most-common use case of a contiguous (B2_start == B1 at
	start of stage 2) S2 with bounds [B1,B2], we run S2 not from (B2_start = B1) to B2 but rather from
	B2/psmall to B2, where psmall = smallest prime which is not a factor of our bigstep D: *****/
	if(!qlo) {			// If qlo unset, set = default stage 2 starting point ... if qlo already set via restart-file
		qlo = B2_start;	// read, it will automatically be > our small-prime-relocation-reflecting value of B2_start.
		if(psmall && B2_start > B1) {	// If psmall = 0, it's an S2 continuation run, no relocation done
			sprintf(cbuf,"Small-prime[%u] relocation: will start Stage 2 at bound %" PRIu64 "\n",psmall,qlo);
			mlucas_fprint(cbuf,pm1_standlone+1);
		}
	}
	/* [b] Stage 2 starts at q0, the smallest multiple of D nearest but not exceeding qlo + D/2;
	once have that, compute k0 = q0/D: */
	q0 = qlo + bigstep - qlo%bigstep;
	q_old_10M = (uint32)(q0 * inv10m);
	// The 1st clause is is to ensure that Stage 2 lower bound > D/2, otherwise we get some negative q's on the q - b[i] side:
	if(qlo > bigstep/2 && q0 > (qlo + bigstep/2)) {
		q0 -= bigstep;
	}
	if(q0 > qlo + bigstep/2) q0 -= bigstep;
	k = k0 = q0/bigstep;	// Now set k to its 'real' value
	if((uint64)k*bigstep != q0) {
		sprintf(cbuf,"k must be 32-bit!");
		mlucas_fprint(cbuf,pm1_standlone+1);	ASSERT(0,cbuf);
	}
	sprintf(cbuf,"Stage 2 q0 = %" PRIu64 ", k0 = %u\n",q0,k0);
	mlucas_fprint(cbuf,pm1_standlone+1);
	/*
	Expanded-match-window scheme needs us to precompute singleton-prime-q's bitmap corresponding to M intervals
	bigstep = D wide centered on q0, with a few key distinctions between M odd and even:

	o For M odd, the midpoint of each map word is a q-value which is a multiple of the bigstep D;
		the minimal-memory (M = 1) prime-pairing algo consists of matching primes symmetric about the midpoint
		of the 0-interval, which is the only interval for M = 1. For M odd > 1, the extended-matching consists
		of further matching prime pairs, still symmetrically about the midpoint of the 0-interval, but in the
		M-1 D-wide interval *pairs* bracketing the 0-interval.

	o For M even, there is no 0-interval. Instead, our map-midpoint, still a multiple of D, is the boundary between
		the middle 2 intervals and the extended-match prime-pairing algo again consists of matching primes symmetric
		about this midpoint. This means we init the primes in fresh intervals entering from the right end of the
		current map differently than for M odd.

	Init M 2*num_b-bit words of pairing bitmap = 0; e.g. for D = 210 we have M 48-bit words. The low num_b bits
	of each word correspond to singleton (in the sense that the corr. q2 is composite) q1 = (q - b[i])'s relative
	to the center of each of the M D-intervals, the hi bits to the q2 = (q + b[i])'s.
	*/
	m2 = m/2;	// If m odd, m2 = (m-1)/2 = # of D-interval on either side of the central 0-interval
	/*
	We save ourselves some awkward preprocessing by starting the stage 2 loop such that at end of the first loop pass,
	the D-interval centered on q0 (M odd) or just left of q0 (M even) enters the as the new upper-interval.
	We only start actual 0-interval and extended-window pairing when said interval has shifted to the middle
	of the extended pairing window, i.e. is the 0-interval (M odd), or shifted just left of the map midpoint (M even):
	*/
	ASSERT(q0 > (m2+1)*(uint64)bigstep, "ERROR: qlo underflows in p-1 stage 2.");
	qlo = q0 - (m2+1)*(uint64)bigstep;
	/*
	[c] Our A^(a^2) values = A^((k*D)^2) and we'll be incrementing k between sweeps over the set of b's.
	(k+1)^2 = k^2 + 2*k + 1, thus A^(((k+1)*D)^2) = A^(D^2*(k^2 + 2*k + 1)) = A^(D^2*k^2) * A^(D^2*(2*k + 1)),
	so for a given starting value k0 precompute the following set of multipliers (all mod n):
	*/
	// At this point pow = A[stage 1 residue]; need either A^(D^2) or (A^D + A^-D), where D = bigstep:
#ifndef PM1_STANDALONE
	snprintf(cbuf,STR_MAX_LEN*2, "Computing Stage 2 loop-multipliers...\n");	mlucas_fprint(cbuf,pm1_standlone+1);
	MME = 0.0;	// Reset maxROE
	// Raise A to power D^2, using mult[0] as a scratch array; again crap-API forces us to specify an "input is pure-int?" flag:
	input_is_int = TRUE;

  #if !USE_PP1_MULTS		// Basic version, needs 2 modmul per D-loop:

	// pow = A^(D^2), in-place, using mult[0] as a scratch array:
	modpow(pow, mult[0], input_is_int, (uint64)bigstep*bigstep, func_mod_square, p, n, scrnFlag,&tdif2);
	// mult[0] = A^((k0*D)^2) = (A^D^2)^(k0^2) = pow^(k0^2), using mult[2] as a scratch array:
	memcpy(mult[0],pow,nbytes);
	input_is_int = FALSE;	// modpow() leaves result in fwdFFT-pass1-done form, so negate this flag for further raising steps
	modpow(mult[0], mult[2], input_is_int, (uint64)k0*k0, func_mod_square, p, n, scrnFlag,&tdif2);
	// mult[1] = A^(D^2*(2*k0 + 1)) = pow^(2*k0+1), using mult[2] as a scratch array:
	memcpy(mult[1],pow,nbytes);
	modpow(mult[1], mult[2], input_is_int, (uint64)2*k0+1, func_mod_square, p, n, scrnFlag,&tdif2);
	// mult[2] = A^(2*D^2) = pow^2:
	memcpy(mult[2],pow,nbytes);
	ierr = func_mod_square(mult[2], 0x0, n, 0,1,        (uint64)mode_flag, p, scrnFlag,&tdif2, FALSE, 0x0);
	if(ierr == ERR_INTERRUPT) {
		return ierr;
	}

  #else	// USE_PP1_MULTS = True: Inits here more involved but then needs just 1 Lucas-sequence-modmul per loop:
	/*
	Montgomery ["Speeding the Pollard and Elliptic Curve Methods of Factorization", Math. Comp. 48(177), pp 243-264, 1987,
	https://doi.org/10.1090/S0025-5718-1987-0866113-7 ], in the p+1 section of the paper suggests, with variable names
	replaced to match mine above:

		"In Subsection 4.1, we used f(D) = A^(D^2) mod N. Another acceptable selection is
			f(D) = (A^D + A^-D) = V[D](A^1 + A^-1) mod N.
		This is well defined since GCD(A,N) = 1."

	Here, V[n] is the (n)th element of the Lucas sequence of the 2nd kind.
	Q: Does this preclude p-1 on known-stage-1-factor runs, i.e. ones where we seek additional factors using stage 2?
	A: No, because have a stage 1 factor if GCD(A-1,N) > 1, not GCD(A,N) > 1.

	The identity for index-incrementing we need is

		V[n+j] = V[n]*V[j] - V[n-j] . [*]

	Implementation of [*] ideally needs an FFT-based fused multiply-subtract.
	For first-look impl just compute V[n]*V[j] via FFT-mul, fwd-FFT the reult, and subtract V[n-j], which is stored in fwd-FFTed form:
	*/
	// vec1,vec2 hold stage 1 residue A and its mod-inverse in packed-bit form,
	// but for the (p+1)-style powering pow[] also still holds A in balanced-digit FP form:
	memcpy(mult[2],pow,nbytes);
	input_is_int = TRUE;
  /******************************************************************************************************/
  /********** Compute A^+[D,(k0-1)*D,k0*D], then add the ± components of each same-power pair: **********/
  /******************************************************************************************************/
	// mult[2] = A^D, in-place, using pow as a scratch array:
	modpow(mult[2], pow, input_is_int, (uint64)bigstep, func_mod_square, p, n, scrnFlag,&tdif2);
	input_is_int = FALSE;	// modpow() leaves result in fwdFFT-pass1-done form, so negate this flag for further raising steps
	// mult[1] = A^((k0-1)*D), using pow as a scratch array:
	memcpy(mult[1],mult[2],nbytes);
	modpow(mult[1], pow, input_is_int, (uint64)k0-1, func_mod_square, p, n, scrnFlag,&tdif2);
	// mult[1] gets fwdFFTed:
	mode_flag = 3;	// bit 0 of mode_flag = 1 since modpow() leaves result in fwdFFT-pass1-done form;
					// bit 1 of mode_flag = 1 since all FFT-mul outputs will be getting re-used as inputs
	ierr += func_mod_square(mult[1], 0x0, n, 0,1, 4ull + (uint64)mode_flag, p, scrnFlag,&tdif2, FALSE, 0x0);
	// mult[0] = A^D * A^((k0-1)*D) = A^(k0*D):
	memcpy(mult[0],mult[2],nbytes);
	ierr += func_mod_square(mult[0], 0x0, n, 0,1, (uint64)mult[1] + (uint64)mode_flag, p, scrnFlag,&tdif2, FALSE, 0x0);
  /******************************************************************************************************/
  /********** Repeat above sequence for A^-1, adding resulting 3 vectors to mult[0-2] as we go: *********/
  /******************************************************************************************************/
	convert_res_bytewise_FP((uint8*)vec2,  a, n, p);
	input_is_int = TRUE;
	// a[] = A^-D, in-place, using pow as a scratch array:
	modpow(a, pow, input_is_int, (uint64)bigstep, func_mod_square, p, n, scrnFlag,&tdif2);

	// Debug: check inverse A^-D by mpying A^+D (mult[2]) and A^-D (a[]):
   #if 0
	// a[] = fwdFFT(A^-D):							   vvvv + 4 to effect "Do in-place forward FFT only"; bit 0 of mode_flag (the only one fwd-FFT cares about) = 1, since mult[0] enters in pure-int form
	ierr += func_mod_square(      a, 0x0, n, 0,1,      4ull + (uint64)mode_flag, p, scrnFlag,&tdif2, FALSE, 0x0);
	double dsum = 0.0; for(i = 0; i < npad; i++) { dsum += a[i]; }
	fprintf(stderr,"fwdFFT(A^-D) has element 0 = %20.5f, sum = %20.5f\n",a[0],dsum);	// [0] = 66475254.49010, sum = 69950177280.00076
	// mult[0] = A^+D * A^-D, check that result = 1 as expected:
	mode_flag = 1;	// bits 0:1 of mode_flag = 1,0, since mult[2] enters in fwd-FFT-pass-1-done form and want output in pure-int form
	ierr += func_mod_square(mult[2], 0x0, n, 0,1, (uint64)a + (uint64)mode_flag, p, scrnFlag,&tdif2, FALSE, 0x0);
	ASSERT(mult[2][0] == 1.0, "inverse-check fails!");
	for(i = 1; i < npad; i++) {
		ASSERT(mult[2][i] == 0.0, "inverse-check fails!");
	}
	fprintf(stderr,"A^-D inverse check passed ... exiting.\n");
	exit(0);	// Since above debug overwrites mult[2], must quit.
   #endif

	// mult[2] = V[D] = A^D + A^-D:
	for(i = 0; i < npad; i++) { mult[2][i] += a[i]; }
	// Normalization-multiply by fwdFFT(1):
	ierr += func_mod_square(mult[2], 0x0, n, 0,1, (uint64)vone + (uint64)mode_flag, p, scrnFlag,&tdif2, FALSE, 0x0);

   #if 0
	#ifndef pp1_vd2_cross_check
		#error The V[D] cross-check code needs the above pp1_vd2_cross_check-defining snippet to be enabled!
	#endif
	// PP1 Debug: compare V[D] computed this way vs alternate-path V[D] computed in above buffer-init code:
	memcpy(a,mult[2],nbytes);
	// Need to undo pass 1 of fwd-FFT; do this just as with fwd-FFT-only, but with flag = 8 instead of 4:
	ierr += func_mod_square(a, 0x0, n, 0,1, 8ull, p, scrnFlag,&tdif2, FALSE, 0x0);
	convert_res_FP_bytewise(a,(uint8*)vec2, n, p, &Res64, &Res35m1, &Res36m1);
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE) {
		// Get the carry bit, bit p. Normalized residue only uses bits <0:p-1>:
		bit = p&63; word = p>>6; tmp = vec2[word]>>bit;
		q = mi64_add_scalar(vec2, tmp, vec2, nlimb);
	} else {
		tmp = vec2[nlimb-1];
		q = mi64_sub_scalar(vec2, tmp, vec2, nlimb-1);	// Again shorten vec-sub-length by 1 to properly check for borrow-out
	}
	ASSERT(q == 0ull,"Properly normalized vec1+vec2 wraparound carry should never have borrow-out!");
	// Now compare to alternate-path V[D] computed in above buffer-init code:
	ASSERT(mi64_cmp_eq(vec1,vec2,nlimb), "V[D] results mismatch!");
	fprintf(stderr,"V[D] cross-check passed ... exiting.\n");
	exit(0);	// Since above debug overwrites a[], must quit.
   #endif	// PM1_DEBUG?

	input_is_int = FALSE;	// modpow() leaves result in fwdFFT-pass1-done form, so negate this flag for further raising steps
	// pow[] = (A^-D)^(k0-1) = A^-((k0-1)*D), in-place, using a[] as a scratch array:
	memcpy(pow,a,nbytes);	// See "Problem" comment below for why we use pow[] as modpow input
							// and a[] as scratch array here rather than simply doing in-place with a[] as input.
	modpow(pow, a, input_is_int, (uint64)k0-1, func_mod_square, p, n, scrnFlag,&tdif2);

	// Debug: check inverse A^-((k0-1)*D) by mpying A^+((k0-1)*D) (mult[1], already fully fwdFFTed) and A^-((k0-1)*D) (pow[]):
   #if 0
	// pow[] = A^+((k0-1)*D) * A^-((k0-1)*D), check that result = 1 as expected:
	mode_flag = 1;	// bits 0:1 of mode_flag = 1,0, since pow[] enters in fwd-FFT-pass-1-done form and want output in pure-int form
	ierr += func_mod_square(pow, 0x0, n, 0,1, (uint64)mult[1] + (uint64)mode_flag, p, scrnFlag,&tdif2, FALSE, 0x0);
	ASSERT(pow[0] == 1.0, "A^-(k0-1)*D inverse-check fails!");
	for(i = 1; i < npad; i++) {
		ASSERT(pow[i] == 0.0, "A^-(k0-1)*D inverse-check fails!");
	}
	fprintf(stderr,"A^-(k0-1)*D inverse check passed ... exiting.\n");
	exit(0);	// Since above debug overwrites pow[], must quit.
   #endif

	// mult[1] is fully fwdFFTed, so pow[] needs the same treatment prior to add:
	ierr += func_mod_square(pow, 0x0, n, 0,1, 4ull + (uint64)mode_flag, p, scrnFlag,&tdif2, FALSE, 0x0);
	// mult[1] = A^((k0-1)*D) + A^-((k0-1)*D):
	for(i = 0; i < npad; i++) { mult[1][i] += pow[i]; }
	// Normalization-multiply by fwdFFT(1):
	ierr += func_mod_square(mult[1], 0x0, n, 0,1, (uint64)vone + 0xC + mode_flag, p, scrnFlag,&tdif2, FALSE, 0x0);

	/* Problem: destroyed A^-D in computing A^-((k0-1)*D) ... how can we compute A^-D * A^-((k0-1)*D) = A^-(k0*D)
	with just those 2 double-float residue storage arrays and A^-D destroyed? Here's where we get absolutely filthy:
	We know that as a side effect, modpow() uses the scratch-array to store the fwd-FFT of the input when the exponent
	is not a power of 2 (in the power of 2 case, the scratch-array is unused). So if (k0-1) != 2^n, a[] holds FFT(A^-D)
	at this point; if (k0-1) == 2^n, a[] holds A^-D and still needs fwd-FFTing: */
	if(isPow2(k0-1)) {
		ierr += func_mod_square(a, 0x0, n, 0,1, 4ull + (uint64)mode_flag, p, scrnFlag,&tdif2, FALSE, 0x0);
	}
	// pow = A^-(k0*D) = A^-D * A^-((k0-1)*D) ... unlike positive-powers case, both inputs here already fwd-FFTed:
	ierr += func_mod_square(pow, 0x0, n, 0,1, (uint64)a + 0xC + mode_flag, p, scrnFlag,&tdif2, FALSE, 0x0);

	// Debug: check inverse A^-(k0*D) by mpying A^+(k0*D) (mult[0]) and A^-(k0*D) (pow[]):
   #if 1
	ierr += func_mod_square(mult[0], 0x0, n, 0,1, 4ull + (uint64)mode_flag, p, scrnFlag,&tdif2, FALSE, 0x0);
	mode_flag = 1;	// bits 0:1 of mode_flag = 1,0, since pow[] enters in fwd-FFT-pass-1-done form and want output in pure-int form
	ierr += func_mod_square(pow, 0x0, n, 0,1, (uint64)mult[0] + (uint64)mode_flag, p, scrnFlag,&tdif2, FALSE, 0x0);
	ASSERT(pow[0] == 1.0, "A^-k0*D inverse-check fails!");
	for(i = 1; i < npad; i++) {
		ASSERT(pow[i] == 0.0, "A^-k0*D inverse-check fails!");
	}
	fprintf(stderr,"A^-k0*D inverse check passed ... exiting.\n");
	exit(0);	// Since above debug overwrites pow[], must quit.
   #endif

	// mult[0] = A^(k0*D) + A^-(k0*D):
	for(i = 0; i < npad; i++) { mult[0][i] += pow[i]; }
	// Normalization-multiply by fwdFFT(1):
	ierr += func_mod_square(mult[0], 0x0, n, 0,1, (uint64)vone + (uint64)mode_flag, p, scrnFlag,&tdif2, FALSE, 0x0);

  #endif	// ifdef USE_PP1_MULTS

	if(restart) {	// If restart, convert bytewise-residue S2 accumulator read from file to floating-point form:
		if(!convert_res_bytewise_FP((uint8*)arrtmp, pow, n, p)) {
			snprintf(cbuf,STR_MAX_LEN*2, "ERROR: convert_res_bytewise_FP Failed on primality-test residue read from savefile %s!\n",savefile);
			mlucas_fprint(cbuf,pm1_standlone+1);	ASSERT(0,cbuf);
		}
		// Restart-file-read S2 interim residue in pow[] needs fwd-weight and FFT-pass1-done:
		ierr = func_mod_square(pow, 0x0, n, 0,1, -4ull, p, scrnFlag,&tdif2, FALSE, 0x0);
	} else {	// Init stage 2 accumulator pow = mult[0]; re-use same var as for stage 1 accum, since done with that):
	  #if !USE_PP1_MULTS
		memcpy(pow,mult[0],nbytes);	// pow = A^((k0*D)^2) for regular (p-1); pow = A^(k0*D) + A^-(k0*D) for (p+1)-style stage 2
	  #else
	/************ Q: Does (p+1) code need to re-init pow from S1 residue here? **************/
	   #if 0	// A: No, because pow = A^(k0*D) + A^-(k0*D) is perfectly fine as S2 init-accumulator
		vec1[nlimb-1] = 0ull;
		if(!convert_res_bytewise_FP((uint8*)vec1, pow, n, p)) {
			snprintf(cbuf,STR_MAX_LEN*2, "ERROR: convert_res_bytewise_FP Failed on S1 residue in vec1!\n");
			mlucas_fprint(cbuf,pm1_standlone+1);	ASSERT(0,cbuf);
		}
		// Pure-int S1 residue in pow[] needs fwd-weight and FFT-pass1-done:
		ierr = func_mod_square(pow, 0x0, n, 0,1, -4ull, p, scrnFlag,&tdif2, FALSE, 0x0);
	   #endif
	/*************************************************************************************/
	  #endif
	}
	nerr = 0;
	// mult0-2 need to be fwdFFTed ... if using (p+1)-style multiplier scheme, mult[1] is already so:
	ierr = func_mod_square(mult[0], 0x0, n, 0,1, 4ull + (uint64)mode_flag, p, scrnFlag,&tdif2, FALSE, 0x0);	nerr += ierr;
  #if !USE_PP1_MULTS
	ierr = func_mod_square(mult[1], 0x0, n, 0,1, 4ull + (uint64)mode_flag, p, scrnFlag,&tdif2, FALSE, 0x0); nerr += ierr;
  #endif
	ierr = func_mod_square(mult[2], 0x0, n, 0,1, 4ull + (uint64)mode_flag, p, scrnFlag,&tdif2, FALSE, 0x0); nerr += ierr;
	if(nerr != 0) {
		sprintf(cbuf,"Stage 2 loop-multipliers computation hit one or more fatal errors! Aborting.");
		mlucas_fprint(cbuf,pm1_standlone+1);	ASSERT(0,cbuf);
	}
  #ifdef CTIME
	clock2 = clock();
  #else
	clock2 = getRealTime();
  #endif
	*tdiff = clock2 - clock1; clock1 = clock2;
	snprintf(cbuf,STR_MAX_LEN*2, "Stage 2 loop-multipliers: clocks =%s, MaxErr = %10.9f.\n",get_time_str(*tdiff), MME);
	mlucas_fprint(cbuf,pm1_standlone+1);
	*tdiff = AME = MME = 0.0;	// Reset timer and maxROE, now also init AvgROE
	AME_ITER_START = 0;	// For p-1 stage 2, start collecting AvgROE data immediately, no need t wait for residue to "fill in"
#endif
	nq = np = 0;	// nq = # of stage 2 q's processed; np = #prime-pairs
	ns = 0;			// ns = # singleton-prime q's, i.e. prime q which end up getting paired with a composite
	/*
	There are 2 kinds of modmul in the stage 2 accumulation loop: the first is by an element of the precomputed
	table of the stage 1 residue raised to successive even powers; in an actual bignum-code implementation these
	are stored in forward-transformed form, thus the mult-update modmul costs the same as a mod-autosquare.

	The second kind is to update the stage 2 product accumulator, which also costs the same as a modsquare.
	*/
	nerr = 0;
	// Oct 2021: extend qhi by several times D to make sure we don't orphan any unpaired-primes just below B2:
	qhi = (B2 + m*bigstep);
	for(q = qlo; q < qhi; q += bigstep)
	{
		if(!reloc_on && q >= reloc_start) {	// Start including relocation-semiprimes once S@ passes this point
			reloc_on = TRUE;
			sprintf(cbuf,"Hit q = %" PRIu64 " >= reloc_start[%" PRIu64 "] ... enabling small-prime relocation.\n",q,reloc_start);
			mlucas_fprint(cbuf,pm1_standlone+1);
		}
		// Only start actual 0-interval and extended-window pairing when q hits q0:
		if(q >= q0) {
		  if(m_is_odd) {	// prime-pairing between lo|hi halves of 0-interval only done if M odd
			/* [d] 0-interval q's were already computed once when said interval entered at upper end of expanded-match
			window, but that pass simply set bitmap bits corr. to singleton-prime-q's = 1. Once interval has shifted to
			center of the expanded-match window, recompute same batch of stage 2 powering pairs q1,q2[i] = k*D +- b[i]
			and process the both-q1-and-q2-prime pairs: */
		#ifdef PM1_DEBUG
			fprintf(stderr,"k = %u: q = %" PRIu64 "\n",k,q);
			fprintf(stderr,"Processing 0-interval prime pairs:\n");
		#endif
			map_lo = map + m2*wsize + (wsize>>1);	// ptr to midpoint of 0-interval map word
			for(i = 0; i < num_b; i++) {
				nq += 2;
				p1 = bytevec_test_bit(map_lo,-i-1); p2 = bytevec_test_bit(map_lo,+i);
			#ifdef PM1_DEBUG
				q1 = q - b[i]; q2 = q + b[i];
				bit = pprimeF64(q1,2ull); if(p1 != bit)
					fprintf(stderr,"Mismatch: q1 = %" PRIu64 "[%u], bytevec_test_bit returns %u\n",q1,bit,p1);
				bit = pprimeF64(q2,2ull); if(p2 != bit)
					fprintf(stderr,"Mismatch: q2 = %" PRIu64 "[%u], bytevec_test_bit returns %u\n",q2,bit,p2);
			#endif
				// Skip a given value of i if one or both of q1,q2[i] are composite according to a 2-prp test:
				j = p1+p2;
				if(j < 2) {
					// For M = 1, must process 1-primes here [no need to clear map bit], otherwise we can wait
					// and hope for singletons to get extended-window-matched; this == (!j || (m == 1)):
					if(j < m)
						continue;
				#ifdef PM1_DEBUG
					fprintf(stderr,"\tq1 = %" PRIu64 "[%u], q2 = %" PRIu64 "[%u], 1-prime\n",q1,p1,q2,p2);
				#endif
					ns++;
				} else {
				#ifdef PM1_DEBUG
					fprintf(stderr,"\tq1 = %" PRIu64 "[%u], q2 = %" PRIu64 "[%u], both prime\n",q1,p1,q2,p2);
				#endif
					np++;
				}
				// Clear both paired map bits, irrespective of whether 1 or both = 1:
				bytevec_bclr(map_lo,-i-1); bytevec_bclr(map_lo,+i);
				// [e] For each prime-pair q1,q2 = -+ b[i], update stage 2 accumulator as *= A^((k*D)^2) - buf[i]:
			#ifndef PM1_STANDALONE
				/* Thanks to linearity of FFT, fwdFFT((mult[0] - buf[...]) = fwdFFT(mult[0]) - fwdFFT(buf[...]);
				we don't need to worry about excess ROE here because the 2 inputs are effectively uncorrelated,
				i.e. the vector-subtract acts much the same as a single added radix-2 FFT pass sans twiddles.
				Use scratch array a[] to hold (mult[0] - buf[...]) vectors:
				*/
				// In ||-exec case, thread-data struct prepopulated with mult[0], a, npad, only need fixed-offset buf[i] at runtime:
			 #ifdef USE_VEC_DBL_SUB
			  #ifdef MULTITHREAD
				vec_double_sub(tpool,tdat,buf[i]);
			  #else
				vec_double_sub(mult[0],buf[i],a,npad);	// a[] = (mult[0][] - buf[i][])
			  #endif
				// pow = pow*(mult[0] - buf[i]) % n: Don't increment nmodmul until after call due to ambiguity in eval order of func(i,++i):
				ierr = func_mod_square(pow, 0x0, n, nmodmul,nmodmul+1, (uint64)a       + (uint64)mode_flag, p, scrnFlag,&tdif2, FALSE,    0x0);	nerr += ierr;
			 #else
				ierr = func_mod_square(pow, 0x0, n, nmodmul,nmodmul+1, (uint64)mult[0] + (uint64)mode_flag, p, scrnFlag,&tdif2, FALSE, buf[i]); nerr += ierr;
			 #endif
				if(ierr == ERR_INTERRUPT) {
					return ierr;
				}
				++nmodmul;
			#endif
			}
		  }	// m_is_odd?
			/* Loop over the m2 interval-pairs symmetric about the 0-interval (M odd) or midpoint of the
			current set of M extended-pairing windows (M even) and process any resulting pairings.
			Example:
				m = 3: m2 = 1, i = 1, pair map words (m2-+i) = [0,2]
				m = 4: m2 = 2, pair map words [1,2],[0,3] for i = 1,2. Thus generalized word pair is [m2-i,m2+i-m_is_even]:
			*/
			for(i = 1; i <= m2; i++) {	// start with pair bracketing 0-interval and work outward
			#ifdef PM1_DEBUG
				fprintf(stderr,"Extended-window pair [%u]\n",i);
			#endif
				// In b[]-index terms, LSB of low word [i]th pair this far left of 0-interval ctr:
				tmp = (i+i-1-m_is_even)*num_b;
				// In context of extended-window bitmap word-pairs, lo,hi refer to entire [2*num_b]-bit map words:
				map_lo = map + (m2-i          )*wsize;	// lo = map[m2-i]
				map_hi = map + (m2+i-m_is_even)*wsize;	// hi = map[m2+i-m_is_even];
				// Testing bits of lo,hi starting with 0-bits corr. to working outward through paired D-intervals:
				bytevec_and(lo,map_lo,map_hi,wsize);	// Use lo as temp to hold AND result
				// Process and clear selected paired 1-bits of the 2 map words:
				for(j = 0; j < (num_b<<1); j++) {
				#ifdef PM1_DEBUG
					p1 = bytevec_test_bit(map_lo,j); p2 = bytevec_test_bit(map_hi,j);
					q1 = q - b[tmp+j]; q2 = q + b[tmp+j];
					bit = pprimeF64(q1,2ull); if(p1 != bit)
						fprintf(stderr,"Mismatch: q1 = %" PRIu64 "[%u], bytevec_test_bit returns %u\n",q1,bit,p1);
					bit = pprimeF64(q2,2ull); if(p2 != bit)
						fprintf(stderr,"Mismatch: q2 = %" PRIu64 "[%u], bytevec_test_bit returns %u\n",q2,p2,bit);
				#endif
					// For thus-paired prime-q's, update stage 2 accumulator:
					p1 = bytevec_test_bit(lo,j);
					if(p1) {
						np++;
					#ifdef PM1_DEBUG
						fprintf(stderr,"\tq = %" PRIu64 " -+ %u: q1 = %" PRIu64 "[%u], q2 = %" PRIu64 "[%u], paired singles\n",q,b[tmp+j],q-b[tmp+j],p1,q+b[tmp+j],p2);
					#endif
					#ifndef PM1_STANDALONE
					 #ifdef USE_VEC_DBL_SUB
					  #ifdef MULTITHREAD
						vec_double_sub(tpool,tdat,buf[tmp+j]);
					  #else
						vec_double_sub(mult[0],buf[tmp+j],a,npad);
					  #endif
						// pow = pow*(mult[0] - buf[tmp+j]) % n;
						ierr = func_mod_square(pow, 0x0, n, nmodmul,nmodmul+1, (uint64)a       + (uint64)mode_flag, p, scrnFlag,&tdif2, FALSE,        0x0); nerr += ierr;
					 #else
						ierr = func_mod_square(pow, 0x0, n, nmodmul,nmodmul+1, (uint64)mult[0] + (uint64)mode_flag, p, scrnFlag,&tdif2, FALSE, buf[tmp+j]); nerr += ierr;
					 #endif
						if(ierr == ERR_INTERRUPT) {
							return ierr;
						}
						++nmodmul;
					#endif
						bytevec_bclr_pair(map_lo,map_hi,j);	// Clear the paired 1-bits from the map words
					}
				}
			}
		}	// endif(q == q0)
		/* Process remaining q-singles in map[0], which is about to get shifted off. Remember, since we are
		left of 0-interval, 0-bit of bitmap word corr. to rightmost q of interval (closest to 0-interval): */
		if(m > 1) {		// code below assumes an expanded-match window, i.e. m >= 3
			if(!bytevec_iszero(map+0,wsize)) {
				bytevec_set_eq(lo,map+0,wsize);	// lo = map[0];
				// In b[]-index terms, LSB of low word [m2]th pair this far left of 0-interval ctr:
				tmp = (m-2)*num_b;
				// Testing bits of lo,hi starting with 0-bits corr. to working outward through paired D-intervals:
			#ifdef PM1_DEBUG
				fprintf(stderr,"Word 0, processing remaining p1-prime singletons:\n");
			#endif
				for(i = 0; i < (num_b<<1); i++) {
					// Pair remaining prime-singles in map[0] with corr. high-interval map[M-1] composites:
					if(bytevec_test_bit(lo,i)) {
						ns++;
					#ifdef PM1_DEBUG
						q1 = q-b[tmp+i]; q2 = q+b[tmp+i];
						p1 = pprimeF64(q1,2ull); p2 = pprimeF64(q2,2ull);	// Run q1,q2 through a base-2 Fermat-composite test
						fprintf(stderr,"\tq = %" PRIu64 " -+ %u: q1 = %" PRIu64 "[%u], q2 = %" PRIu64 "[%u], 1-prime\n",q,b[tmp+i],q1,p1,q2,p2);
					#endif
					#ifndef PM1_STANDALONE
					 #ifdef USE_VEC_DBL_SUB
					  #ifdef MULTITHREAD
						vec_double_sub(tpool,tdat,buf[tmp+i]);
					  #else
						vec_double_sub(mult[0],buf[tmp+i],a,npad);
					  #endif
						// pow = pow*(mult[0] - buf[tmp+i]) % n:
						ierr += func_mod_square(pow, 0x0, n, nmodmul,nmodmul+1, (uint64)a       + (uint64)mode_flag, p, scrnFlag,&tdif2, FALSE,        0x0); nerr += ierr;
					 #else
						ierr += func_mod_square(pow, 0x0, n, nmodmul,nmodmul+1, (uint64)mult[0] + (uint64)mode_flag, p, scrnFlag,&tdif2, FALSE, buf[tmp+i]); nerr += ierr;
					 #endif
						if(ierr == ERR_INTERRUPT) {
							return ierr;
						}
						++nmodmul;
					#endif
					}
				}
			}
			// Since current 0-interval will end up leftward of new one after shift, bit-reverse its map word:
			map_hi = map + m2*wsize;
			bytevec_brev(map_hi,map_hi,wsize);
			// shift words in map down one slot, vacating the upper-interval slot:
			map_lo = map + 0; map_hi = map_lo + wsize;
			for(i = 0; i < m-1; i++, map_lo += wsize, map_hi += wsize) {
				bytevec_set_eq(map_lo,map_hi,wsize);	// map[i] = map[i+1]
			}
		}
		/* shift in new upper-interval and tag its within-interval singletons - remember, our bitmap consists of m
		bigstep-D-representing 'words' of 2*num_b bits each, laid out as follows, depending on if M odd or even -
		in either case the midpoint of each map word must be a multiple of the bigstep D:
		M odd:
			Middle word of map (a.k.a. 0-interval) corr. to interval (q-+D/2); Midpoints of our M intervals map[0,...,m-1] are
			[q - m2*D, q - (m2-1)*D, ... , q, ... , q + (m2-1)*D, q + m2*D].
		M even:
			There is no 0-interval; our map-midpoint, still a multiple of D, is the boundary between the middle
			2 intervals. This means we init the primes in fresh intervals entering from the right end of the
			current map differently than for M odd.
		Here we are processing the interval one-beyond the right end of the map, whose midpoint is q + (m2+1)*D
		for both M odd and even:
		*/
		tmp = q + (m2+1)*bigstep;
		/* tmp = midpoint (M odd) or right endpoint (M even) of new upper-interval
		M odd:
			Working outward from middle of interval, init separate num_b-bit accumulators
			for lo,hi halves of bitmap word and then combine post-loop.
		M even:
			Working leftward from right endpoint of interval, init single 2*num_b-bit accumulator.
		*/
	#ifdef PM1_DEBUG
		fprintf(stderr,"New upper-interval with q0 = %" PRIu64 ", tagging its primes:\n",tmp);
	#endif
	/*
		Prime relocation: Illustrate using psmall = 11, but analogous pattering holds for psmall = 7 (D = 330|660):
		Table-lookup scheme for 24-bit (mod psmall) maplets - Note that we don't need to consider M > 2 because the rules apply
		per-width-D-interval.
		psmall = 11: Defining q0mod := q0%11 and computing q0mod for initial interval on loop entry:
		D = 210:
			M odd:	Each new upper interval gets mod-11 map [q0mod.hi,q0mod.lo], q0mod += 1 (mod 11) for next interval
			M even:	Each new upper interval gets mod-11 map [q0mod.lo,(q0mod-1).hi], q0mod += 1 (mod 11) for next interval
		D = 420:
			M odd:	Each new upper interval gets mod-11 map [(q0mod+1).lo,q0mod.hi,q0mod.lo,(q0mod-1).hi], q0mod += 2 (mod 11) for next interval
			M even:	Each new upper interval gets mod-11 map [q0mod.lo,(q0mod-1).hi,(q0mod-1).lo,(q0mod-2).hi], q0mod += 2 (mod 11) for next interval
		D = 840:
			M odd:	Each new upper interval gets mod-11 map [q0mod+(2,...,-2).(lo,hi)], q0mod += 4 (mod 11) for next interval
			M even:	Each new upper interval gets mod-11 map [q0mod+(0,...,-4).(lo,hi)], q0mod += 4 (mod 11) for next interval
	*/
		/* Once small-prime-relocation has begin:
		Copy the needed wsize bytes in (2^bigstep_pow2) rsize-byte 'maplet' chunks from the appropriate one of the
		reloc_mod_psmall_bytemap const init-vectors to rmap in rsize-byte chunks. We first compute ptr to the low (rightmost)
		maplet of the above-described rmap, which may be [q0mod-(something)].lo or [q0mod-(something)].hi, and work our way leftward:
		*/
	  if(reloc_on) {	// Otherwise rmap = 0, set by the initial calloc()
		// nmaplet = 2|4|8 for D = 210|420|840:
		uint32 nmaplet = 1<<bigstep_pow2, q0mod = tmp % psmall;	// Here, tmp holds the current q0
		bytevec_clear(rmap,wsize);
		if(m_is_odd) {
		  if(bigstep_pow2 == 1) {	// Special case - just copy the 2 [q0mod.hi,q0mod.lo] maplets in-order:
			i = q0mod*rsize; bytevec_set_eq(rmap, reloc_mod_psmall_bytemap + (i<<1), 2*rsize);
		  } else {
			i = nmaplet>>2;	// Rightmost-maplet negative-offset-index i = 1|2 for D = 420|840
			MOD_SUB32(q0mod,i,psmall, i);	// i = (q0mod - i)%psmall
			for(j = 0; j < nmaplet; j++) {	// 2*i gets us to the (q0mod-1).[hi,lo] pair, +1 to the hi maplet of the pair
				jeven = !(j & 1);
				bytevec_set_eq(rmap + j*rsize, reloc_mod_psmall_bytemap + (2*i+jeven)*rsize, rsize);	// j even: copy .hi; j odd: copy .lo
				if(jeven)	// Only mod-increment q0mod-offset every other maplet:
					MOD_ADD32(i,1,psmall, i);
			}
		  }
		} else {
			i = nmaplet>>1;	// Rightmost-maplet negative-offset-index i = 1|2|4 for D = 210|420|840
			MOD_SUB32(q0mod,i,psmall, i);	// i = (q0mod - i)%psmall
			for(j = 0; j < nmaplet; j++) {
				jeven = !(j & 1);
				bytevec_set_eq(rmap + j*rsize, reloc_mod_psmall_bytemap + (2*i+jeven)*rsize, rsize);	// j even: copy .hi; j odd: copy .lo
				if(jeven)	// Only mod-increment q0mod-offset every other maplet:
					MOD_ADD32(i,1,psmall, i);
			}
		}
	  }
		 map_hi =  map + (m-1)*wsize;	// ptr to high map word
		bytevec_clear( map_hi,wsize);
		if(m_is_odd) {
			for(i = 0,j = num_b-1; i < num_b; i++,j--) {
				q1 = tmp - b[i]; q2 = tmp + b[i];
				// If q divisible by our relocation prime psmall, first compute quotient q/psmall via
				// mul-by-precomputed-Montgomery inverse (mod 2^64), then check that for 2-PRP-ness.
				// If no relocation being done, rmap = 0 and the test_bit() calls will return FALSE, but no point even doing them:
				if(psmall) {
					if(bytevec_test_bit(rmap,j      )) {
						q1 *= pinv64;
					#ifdef PM1_DEBUG
						fprintf(stderr,"reloc q1: %" PRIu64 " => %" PRIu64 "\n",q1,q1*psmall);
					#endif
					}
					if(bytevec_test_bit(rmap,i+num_b)) {
						q2 *= pinv64;
					#ifdef PM1_DEBUG
						fprintf(stderr,"reloc q2: %" PRIu64 " => %" PRIu64 "\n",q2,q2*psmall);
					#endif
					}
				}
				p1 = pprimeF64(q1,2ull); p2 = pprimeF64(q2,2ull);	// Run q1,q2 through a base-2 Fermat-composite test
				if(p1) bytevec_bset(map_hi,j);
				if(p2) bytevec_bset(map_hi,i+num_b);	// High half of a 2*num_b-bit map word
			}
			/* hi ends up with q nearest the 0-interval (i = 0, the leftmost q in the hi[] half of the current interval)
			in 0-bit, just as we want; lo's leftmost q corr. to i = num_b-1, i.e. lo[] needs bit-reversing.
			*/
		//	bytevec_brev(lo,lo,wsize>>1);	Jun 2021: inlined the BR vai downward-running 2nd index j, and the below bytevec-OR:
		//	bytevec_or(map_hi,lo,hi,wsize);// Since hi already left-shifted num_b bits, sum is just OR of 2 disjoint num_b-bit pieces
		} else {	// M even:
			for(i = 0,j = 2*num_b-1; j >= 0; i++,j--) {	// Note: we declared i,j,k,l as signed
				q1 = tmp - b[i];
				// If q divisible by our relocation prime psmall, first compute quotient q/psmall via
				// mul-by-precomputed-Montgomery inverse (mod 2^64), then check that for 2-PRP-ness:
				if(psmall && bytevec_test_bit(rmap,j)) {
					q1 *= pinv64;
				#ifdef PM1_DEBUG
					fprintf(stderr,"reloc q: %" PRIu64 " => %" PRIu64 "\n",q1,q1*psmall);
				#endif
				}
				p1 = pprimeF64(q1,2ull);
				if(p1) bytevec_bset(map_hi,j);
			}
			// map_hi ends up with q nearest the map midpoint (i = 0, the leftmost q in the current interval) in 0-bit, just as we want.
		}
	#ifdef PM1_DEBUG
		// Print the string-form q%psmall bitmap:
		if(psmall) {
			bytevec_bitstr(rmap,wsize,cbuf);	//********************** Once debugged, mirror map-element BR/shifts in above
			fprintf(stderr,"\tq == 0 (mod %u) bitmap: %s\n",psmall,cbuf);
		}
	#endif

		// Only update mults once start actual 0-interval and extended-window pairing above:
		if(q >= q0) {
			// [f] Each ensuing k++ needs 2 modmuls to update A^((k*D)^2) for the incremented k, again all (mod n):
		#ifndef PM1_STANDALONE
//fprintf(stderr,"MME = %10.8f\n",MME);
		  #if !USE_PP1_MULTS		// Basic version, needs 2 modmul per D-loop:
		/*	mult[0] = mult[0]*mult[1] % n;	// 1. k^2 += 2*k + 1
			mult[1] = mult[1]*mult[2] % n;	// 2. 2*k++
		Since mult[0-2] all fwd-FFTed, this costs 2 x [dyadic-mul, inv-FFT, carry, fwd-FFT] = equivalent of 2 mod-squares.
		*/
			// Only increment nmodmul every 2nd call here, since each call is 1-FFT:
/* [1a]: */	ierr = func_mod_square(mult[0], 0x0, n, nmodmul,nmodmul+1, (uint64)mult[1] + 0xC + mode_flag, p, scrnFlag,&tdif2, FALSE, 0x0); nerr += ierr;
			if(ierr == ERR_INTERRUPT) {
				return ierr;
			}
/* [1b]: */	ierr = func_mod_square(mult[0], 0x0, n, nmodmul,nmodmul+1,            4ull       + mode_flag, p, scrnFlag,&tdif2, FALSE, 0x0); nerr += ierr;
/* [2a]: */	ierr = func_mod_square(mult[1], 0x0, n, nmodmul,nmodmul+1, (uint64)mult[2] + 0xC + mode_flag, p, scrnFlag,&tdif2, FALSE, 0x0); nerr += ierr;
			if(ierr == ERR_INTERRUPT) {
				return ierr;
			}
/* [2b]: */	ierr = func_mod_square(mult[1], 0x0, n, nmodmul,nmodmul+1,            4ull       + mode_flag, p, scrnFlag,&tdif2, FALSE, 0x0); nerr += ierr;
			nmodmul += 2;

		  #else	// USE_PP1_MULTS = True, needs just 1 modmul per D-loop:

			/* Up-multiply uses Lucas sequence of 2nd kind:
				V[n+j] = V[n]*V[j] - V[n-j] . [*]
			Here - all 3 terms stored in fwdFFTed form:
				V[n  ] = A^(k*D) + A^-(k*D), stored in mult[0];
				V[n-j] = A^((k-1)*D) + A^-((k-1)*D), stored in mult[1].
				V[  j] = A^D + A^-D, stored in mult[2];
			The result V[n+j] = A^((k+1)*D) + A^-((k+1)*D) overwrites mult[0], so must save copy of it first,
			since the old value must end up in mult[1]. This would need multiple memcpy to support, which would
			more than offset gain from saving 1 modmul, so instead directly support it via one more modification
			to the dyadic-mul routines, this one to add support for a fused mulsub operation, d = a*b - c,
			with b[] treated as const (our V[j] above), the input a[] (our V[n]) overwriting c[] (our V[n-j])
			after the subtract, and the result d[] (our V[n+j]) overwriting the input a[].
			*/
			// Only increment nmodmul once here, since each call is 1-FFT;
			// The ~() on the 64-bit subtrahend address tells the modmul routine that it's a mulsub rather than a submul:
			/******** For initial impl, use the sequence: ********/
			memcpy(a,mult[0],nbytes);	// Save copy of V[n] = A^(k*D) + A^-(k*D) in a[];
				// V[n] *= V[j] overwrites mult[0]:
			ierr = func_mod_square(mult[0], 0x0, n, nmodmul,nmodmul+1, (uint64)mult[2] + 0xC + mode_flag, p, scrnFlag,&tdif2, FALSE, 0x0); nerr += ierr;
			if(ierr == ERR_INTERRUPT) {
				return ierr;
			}
			ierr = func_mod_square(mult[0], 0x0, n, nmodmul,nmodmul+1,            4ull       + mode_flag, p, scrnFlag,&tdif2, FALSE, 0x0); nerr += ierr;
			++nmodmul;
				// Subtract mult[0] -= V[n-1], yielding V[n+1] = V[n]*V[j] - V[n-1];
		  #ifdef MULTITHREAD
			vec_double_sub(tpool,tdat,mult[1]);
		  #else
			vec_double_sub(mult[0],mult[1],a,npad);
		  #endif
			memcpy(mult[1],a,nbytes);	// Copy of V[n] = A^(k*D) + A^-(k*D) stored in a[] overwrites mult[1].
			/*****************************************************/
		//	ierr += func_mod_square(mult[0], 0x0, n, nmodmul,nmodmul+1, (uint64)mult[2] + 0xC + mode_flag, p, scrnFlag,&tdif2, FALSE, (double*)(~(uint64)mult[1]));
		//	ierr += func_mod_square(mult[0], 0x0, n, nmodmul,nmodmul+1,            4ull       + mode_flag, p, scrnFlag,&tdif2, FALSE,                         0x0);	++nmodmul;
		  #endif	// USE_PP1_MULTS ?
		#endif	// PM1_STANDALONE ?
			k++;
			// (k - k0) = #bigstep-blocks (passes thru above loop) used in stage 2; #modmul = np + ns + 2*(k - k0).
		}	// endif(q >= q0)
	#ifndef PM1_STANDALONE
		// In single-threaded mode, accumulate the cycle count in a floating double on each pass to avoid problems
		// with integer overflow of the clock() result, if clock_t happens to be 32-bit int on the host platform:
	  #ifdef CTIME
		clock2 = clock();	*tdiff += (double)(clock2 - clock1);	clock1 = clock2;
	  #endif
		// Only handle errs of type ROE in p-1 stage 2 - we prefer to handle such before doing any savefile-updating.
		/*** If multiple errs/types per bigstep-loop ever become an issue, change modmul-retval handling to 'nerr |= 1<<ierr;'
		and in stage-2-return-value-handling code check for various errtypes via e.g. 'if(ierr | (1<<ERR_ROUNDOFF))' ***/
		if(nerr) {
			retval = nerr; goto ERR_RETURN;
		}
		/*...Every (ITERS_BETWEEN_CHECKPOINTS)th modmuls, print timings to stdout or STATFILE.
		If it's a regular (i.e. non-timing) test, also write current residue to restart files.
		*/
		if((nmodmul - nmodmul_save) >= ITERS_BETWEEN_CHECKPOINTS || (q+bigstep) >= qhi) {
			// Copy current S2 residue into a[] and undo pass 1 of fwd-FFT:
			memcpy(a,pow,nbytes);
			ierr = func_mod_square(a, 0x0, n, 0,1, 8ull, p, scrnFlag,&tdif2, FALSE, 0x0);
			arrtmp[nlimb-1] = 0ull;
			convert_res_FP_bytewise(a, (uint8*)arrtmp, n, p, &Res64, &Res35m1, &Res36m1);
		  #ifdef RTIME
			clock2 = getRealTime();	*tdiff = clock2 - clock1;	clock1 = clock2;
		  #endif
			// Get a timestamp:
			calendar_time = time(NULL);
			local_time = localtime(&calendar_time);
			strftime(timebuffer,SIZE,"%Y-%m-%d %H:%M:%S",local_time);
			AME /= (nmodmul - nmodmul_save);
			// Print [date in hh:mm:ss | p | stage progress | %-complete | time | per-iter time | Res64 | max ROE:
			snprintf(cbuf,STR_MAX_LEN*2, "[%s] %s %s = %" PRIu64 " [%5.2f%% complete] clocks =%s [%8.4f msec/iter] Res64: %016" PRIX64 ". AvgMaxErr = %10.9f. MaxErr = %10.9f.\n"
				, timebuffer, PSTRING, "S2 at q", q+bigstep, (float)(q-B2_start)/(float)(B2-B2_start) * 100,get_time_str(*tdiff)
				, 1000*get_time(*tdiff)/(nmodmul - nmodmul_save), Res64, AME, MME);
			mlucas_fprint(cbuf,pm1_standlone+scrnFlag);
			*tdiff = MME = 0.0;	// Reset timer and maxerr at end of each iteration interval
			fp = mlucas_fopen(savefile, "wb");
			if(fp) {
				// q won't get += bigstep until we loop, so here, (q + bigstep) is the q-value corr. to just-incremented k.
				// Also write any relocation-prime psmall into high bit of the resulting nsquares field:
				write_ppm1_savefiles(savefile,p,n,fp, ((uint64)psmall<<56) + q + bigstep, (uint8*)arrtmp,Res64,Res35m1,Res36m1, 0x0,0x0,0x0,0x0);
				fclose(fp);	fp = 0x0;
			} else {
				snprintf(cbuf,STR_MAX_LEN*2, "ERROR: unable to open restart file %s for write of checkpoint data.\n",savefile);
				mlucas_fprint(cbuf,pm1_standlone+1);	ASSERT(0,cbuf);
			}
			// If interim-GCDs enabled (default) and latest S2 interval crossed a 10M mark, take a GCD; if factor found, early-return;
			if(interim_gcd) {
				if((q_div_10M = (uint32)(q * inv10m)) > q_old_10M) {
					q_old_10M = q_div_10M;
					i = gcd(2,p,arrtmp,0x0,nlimb,gcd_str);	// 1st arg = stage just completed
					if(i) {
						B2 = q;	// Reset B2 to reflect the actual interval run
						goto S2_RETURN;
					}
				#ifdef CTIME
					clock2 = clock();	*tdiff += (double)(clock2 - clock1);	clock1 = clock2;
				#else
					clock2 = getRealTime();	*tdiff = clock2 - clock1;	clock1 = clock2;
				#endif
				}
			}
			nmodmul_save = nmodmul;
		}
	#endif	// #ifndef PM1_STANDALONE
	}	// endfor(q = qlo; q < qhi; q += bigstep)
	ASSERT(nerr == 0, "Stage 2 loop hit a modmul error!");
#ifndef PM1_STANDALONE
	// Need to undo pass 1 of fwd-FFT on loop-exit; do this just as with fwd-FFT-only, but with flag = 8 instead of 4:
	ierr = func_mod_square(pow, 0x0, n, 0,1, 8ull, p, scrnFlag,&tdif2, FALSE, 0x0);
	arrtmp[nlimb-1] = 0ull;
	convert_res_FP_bytewise(pow, (uint8*)arrtmp, n, p, &Res64, &Res35m1, &Res36m1);
S2_RETURN:
#endif
	// (k - k0) = #bigstep-blocks (passes thru above loop) used in stage 2; np + ns + 2*(k - k0) = #modmul:
	nmodmul = np + ns + 2*(k - k0);	// This is actually redundant, but just to spell it out
	snprintf(cbuf,STR_MAX_LEN*2,"M = %2u: #buf = %4u, #pairs: %u, #single: %u (%5.2f%% paired), #blocks: %u, #modmul: %u\n",m,m*num_b,np,ns,100.0*2*np/(2*np+ns),k-k0,nmodmul);
	mlucas_fprint(cbuf,pm1_standlone+1);
#ifndef PM1_STANDALONE

  #ifdef PM1_DEBUG
  #warning Revert this preprocessor flag!
	fprintf(stderr,"Res64 = %#016" PRIX64 "; clocks =%s, MaxErr = %10.9f\n",arrtmp[0],get_time_str(*tdiff),MME);
  if(p == 33554432) {  // F25: check if the known factor divides the S2 result:
	ASSERT(MODULUS_TYPE == MODULUS_TYPE_FERMAT, "This p-1 self-test requires Fermat-mod mode!");
	int isfact = mi64_is_div_by_scalar64(arrtmp,2170072644496392193ull,nlimb);	// k = 2^5.3^2.37.997.11066599
	ASSERT(isfact != 0, "Failed to find known stage 2 factor!");
	fprintf(stderr,"%s p-1 known-stage-2 prime self-test success!\n",PSTRING);
  }
  if(p == 108268067) {
	ASSERT(MODULUS_TYPE == MODULUS_TYPE_MERSENNE, "This p-1 self-test requires Mersenne-mod mode!");
	uint64 rem[2] = {0ull,0ull}, q[2] = {11943519037290122063ull,18561975ull};	// k = 7.17.19.61.294313.38955941; q = k.2^p + 1
	int isfact = mi64_div(arrtmp,q, nlimb,2, 0x0, rem);
	ASSERT(isfact != 0, "Failed to find known stage 2 factor!");
	fprintf(stderr,"%s p-1 known-stage-2 prime self-test success!\n",PSTRING);
  }
  if(p == 2147483648) {  // F31: check if the known factor divides the S2 result:
	ASSERT(MODULUS_TYPE == MODULUS_TYPE_FERMAT, "This p-1 self-test requires Fermat-mod mode!");
	uint64 rem[2] = {0ull,0ull}, q[2] = {3118754346955702273ull,2544ull};	// k = 3.13.140091319777; q = k.2^(m+2) + 1
	int isfact = mi64_div(arrtmp,q, nlimb,2, 0x0, rem);
	ASSERT(isfact != 0, "Failed to find known stage 2 factor!");
	fprintf(stderr,"%s p-1 known-stage-2 prime self-test success!\n",PSTRING);
  }
  #endif	// PM1_DEBUG

	// In case of normal (non-early) return, caller will handle the GCD:
	if(strlen(gcd_str)) {
		snprintf(cbuf,STR_MAX_LEN*2, "Stage 2 early-return due to factor found; MaxErr = %10.9f.\n",MME);
	} else {
		snprintf(cbuf,STR_MAX_LEN*2, "Stage 2 done; MaxErr = %10.9f. Taking GCD...\n",MME);
	}
	mlucas_fprint(cbuf,pm1_standlone+scrnFlag);
#endif
ERR_RETURN:
	// Free the memory:
	free((void *)a_ptmp); a_ptmp = a = 0x0; buf = 0x0;
	free((void *)b); b = 0x0;
	free((void *)map); map = 0x0;
  #ifdef MULTITHREAD
	free((void *)thr_ret ); thr_ret  = 0x0;
	free((void *)thread  ); thread   = 0x0;
	free((void *)tdat    ); tdat     = 0x0;
  #endif
	return retval;
}

/************ Parallel|SIMD utility functions for p-1 Stage 2, a.k.a. "Amdahl's Law section": ************/

// SIMD n-double vector subtract: c[] = a[] - b[]. Assumes inputs properly aligned
// on any applicable SIMD-related memory addresses. Allows in-place:
#ifndef PM1_STANDALONE
  #ifdef MULTITHREAD
	void vec_double_sub(struct threadpool *tpool, struct pm1_thread_data_t *tdat, double c[])
	{
		// Threadpool-based dispatch:
		// First 3 subfields same for all threads, 4th provides thread-specifc data, will be inited with fixed data
		// (mult[0] + offset, per-thread chunksize) here, variable ones (subtrahend buf[i] + offset) at thread dispatch:
		static task_control_t task_control = {NULL, (void*)vec_double_sub_loop, NULL, 0x0};
		static int task_is_blocking = TRUE;
		uint32 i;
		for(i = 0; i < NTHREADS; ++i) {
			// Add fixed-offset represented by the address of the subtrahend-array c[] to each
			// precomputed datachunk offset index. Pointer arithmetic takes case of the *= 8 scaling,
			// but first cast index stored in tdat[i].arr2 to int to avoid illegal operation addition of pointers:
			tdat[i].arr2 = c + (intptr_t)(tdat[i].arr2);
			task_control.data = (void*)(&tdat[i]);
		//	printf("adding pool task %d\n",i);
			threadpool_add_task(tpool, &task_control, task_is_blocking);
		//	printf("; #tasks = %d, #free_tasks = %d\n", tpool->tasks_queue.num_tasks, tpool->free_tasks_queue.num_tasks);
		}
	//	printf("start; #tasks = %d, #free_tasks = %d\n", tpool->tasks_queue.num_tasks, tpool->free_tasks_queue.num_tasks);
		struct timespec ns_time;	// We want a sleep interval of 10 nSec here...
		ns_time.tv_sec  =  0;	// (time_t)seconds - Don't use this because under OS X it's of type __darwin_time_t, which is long rather than double as under most linux distros
		ns_time.tv_nsec = 10;	// (long)nanoseconds - Get our desired 0.1 mSec as 10^5 nSec here

	//	while(tpool->tasks_queue.num_tasks != 0) {	//*** not safe, since can have #tasks == 0 with some tasks still in flight ***
		while(tpool->free_tasks_queue.num_tasks != NTHREADS) {
			ASSERT(0 == mlucas_nanosleep(&ns_time), "nanosleep re-call-on-signal fail!");
		//	printf("sleep; #tasks = %d, #free_tasks = %d\n", tpool->tasks_queue.num_tasks, tpool->free_tasks_queue.num_tasks);
		}
	//	printf("end  ; #tasks = %d, #free_tasks = %d\n", tpool->tasks_queue.num_tasks, tpool->free_tasks_queue.num_tasks);
		for(i = 0; i < NTHREADS; ++i) {
			// Pointer subtraction is legal, and automatically undoes the earlier *8 pointer-arithmetic scaling, but result
			tdat[i].arr2 = (double *)(tdat[i].arr2 - c);	// treated as integer, so cast back to double* to avoid warning.
		}
	}
  #else
	// In unthreaded version, vec_double_sub() is just a wrapper for a single loop-exec:
	void vec_double_sub(double a[], double b[], double c[], uint32 n)
	{
		vec_double_sub_loop(a,b,c,n);
	}
  #endif

  #ifdef MULTITHREAD
	void*vec_double_sub_loop(void*targ)	// Thread-arg pointer *must* be cast to void and specialized inside the function
	{
		struct pm1_thread_data_t* thread_arg = targ;
		int thr_id = thread_arg->tid;
		double*c = thread_arg->arr0;	// Output array0 = c[], including address-offset into a given thread-processed chunk
		double*a = thread_arg->arr1;	// Input  array1 = a[], ditto
		double*b = thread_arg->arr2;	// Input  array2 = b[], ditto
		int n    = thread_arg->n;		// Chunksize
  #else
	void vec_double_sub_loop(double a[], double b[], double c[], uint32 n)
	{
		int thr_id = 0;	/* In unthreaded mode this must always = 0 */
  #endif
		// Inner asm-macro processes 8 RE_IM_STRIDE-double vec section per pass, #passes = n >> (3+(L2_SZ_VD-3)) = n >> L2_SZ_VD:
		uint32 i, nloop = n >> L2_SZ_VD;
	//	fprintf(stderr,"Thread %d: vec_double_sub_loop with nloop = %u...\n",thr_id,nloop);
	#ifndef USE_SSE2
		// Cleanup loop shared by SIMD and non-SIMD; in the latter case it does all the work:
		nloop = 0;
	#else
		double *aptr = a, *bptr = b, *cptr = c;
	//	for(i = 0; i < nloop; i++) {
	  #ifdef USE_ARM_V8_SIMD
		__asm__ volatile (\
		"ldr w0,%[__nloop]	\n\t"\
			"ldr	x1,%[__aptr]	\n\t	ldr	x2,%[__bptr]	\n\t	ldr	x3,%[__cptr]	\n\t"\
		"loop0:	\n\t"\
			"ldp	q0,q1,[x1]				\n\t	ldp	q8 ,q9 ,[x2      ]		\n\t"\
			"ldp	q2,q3,[x1,#0x20]		\n\t	ldp	q10,q11,[x2,#0x20]		\n\t"\
			"ldp	q4,q5,[x1,#0x40]		\n\t	ldp	q12,q13,[x2,#0x40]		\n\t"\
			"ldp	q6,q7,[x1,#0x60]		\n\t	ldp	q14,q15,[x2,#0x60]		\n\t"\
			"add	x1,x1,#0x80				\n\t	add	x2,x2,#0x80				\n\t"\
			"fsub	v0.2d,v0.2d,v8.2d		\n\t	fsub	v1.2d,v1.2d,v9.2d	\n\t"\
			"fsub	v2.2d,v2.2d,v10.2d		\n\t	fsub	v3.2d,v3.2d,v11.2d	\n\t"\
			"fsub	v4.2d,v4.2d,v12.2d		\n\t	fsub	v5.2d,v5.2d,v13.2d	\n\t"\
			"fsub	v6.2d,v6.2d,v14.2d		\n\t	fsub	v7.2d,v7.2d,v15.2d	\n\t"\
			"stp	q0,q1,[x3]				\n\t"\
			"stp	q2,q3,[x3,#0x20]		\n\t"\
			"stp	q4,q5,[x3,#0x40]		\n\t"\
			"stp	q6,q7,[x3,#0x60]		\n\t	add	x3,x3,#0x80				\n\t"\
		"sub	w0,w0,#1	\n\t"/* decrement loop counter */\
		"cmp	w0,0		\n\t"/* loop end; continue is via jump-back if w0 != 0 */\
		"bgt	loop0		\n\t"\
			:					// outputs: none
			: [__nloop] "m" (nloop)	/* All inputs from memory addresses here */\
			 ,[__aptr] "m" (aptr)	\
			 ,[__bptr] "m" (bptr)	\
			 ,[__cptr] "m" (cptr)	\
			: "cc","memory","x0","x1","x2","x3",\
			"v0","v1","v2","v3","v4","v5","v6","v7","v8","v9","v10","v11","v12","v13","v14","v15"	/* Clobbered registers */\
		);
	  #elif defined(USE_AVX512)
		__asm__ volatile (\
		"movslq	%[__nloop],%%rsi	\n\t"/* ASM loop structured as for(i = nloop; i != 0; i--){...} */\
			"movq %[__aptr],%%rax	\n\t	movq %[__bptr],%%rbx	\n\t	movq %[__cptr],%%rcx	\n\t"\
		"0:	\n\t"/* Note: Clang did not like longer label here, e.g. vdsub_loop in place of 0 gave "invalid symbol redefinition" error on __nloop */\
			"vmovaps	     (%%rax),%%zmm0	\n\t	vmovaps	0x040(%%rax),%%zmm1	\n\t"\
			"vmovaps	0x080(%%rax),%%zmm2	\n\t	vmovaps	0x0c0(%%rax),%%zmm3	\n\t"\
			"vmovaps	0x100(%%rax),%%zmm4	\n\t	vmovaps	0x140(%%rax),%%zmm5	\n\t"\
			"vmovaps	0x180(%%rax),%%zmm6	\n\t	vmovaps	0x1c0(%%rax),%%zmm7	\n\t"\
			"vsubpd	     (%%rbx),%%zmm0,%%zmm0	\n\t	vsubpd	0x040(%%rbx),%%zmm1,%%zmm1	\n\t"\
			"vsubpd	0x080(%%rbx),%%zmm2,%%zmm2	\n\t	vsubpd	0x0c0(%%rbx),%%zmm3,%%zmm3	\n\t"\
			"vsubpd	0x100(%%rbx),%%zmm4,%%zmm4	\n\t	vsubpd	0x140(%%rbx),%%zmm5,%%zmm5	\n\t"\
			"vsubpd	0x180(%%rbx),%%zmm6,%%zmm6	\n\t	vsubpd	0x1c0(%%rbx),%%zmm7,%%zmm7	\n\t"\
			"vmovaps	%%zmm0,     (%%rcx)	\n\t	vmovaps	%%zmm1,0x040(%%rcx)	\n\t"\
			"vmovaps	%%zmm2,0x080(%%rcx)	\n\t	vmovaps	%%zmm3,0x0c0(%%rcx)	\n\t"\
			"vmovaps	%%zmm4,0x100(%%rcx)	\n\t	vmovaps	%%zmm5,0x140(%%rcx)	\n\t"\
			"vmovaps	%%zmm6,0x180(%%rcx)	\n\t	vmovaps	%%zmm7,0x1c0(%%rcx)	\n\t"\
			/* [a|b|c]ptr += SZ_VD: Double-pointer arithmetic, thus += SZ_VD*8 in raw bytewise form. AVX-512 has SZ_VD = 64: */\
			"addq $0x200,%%rax	\n\t	addq $0x200,%%rbx	\n\t	addq $0x200,%%rcx	\n\t"\
		"decq	%%rsi			\n\t"\
		"jnz 0b	\n\t"/* loop end; continue is via jump-back if rsi != 0 */\
			:	/* outputs: none */\
			: [__nloop] "m" (nloop)	/* All inputs from memory addresses here */\
			 ,[__aptr] "m" (aptr)	\
			 ,[__bptr] "m" (bptr)	\
			 ,[__cptr] "m" (cptr)	\
			: "cc","memory","rax","rbx","rcx","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"	/* Clobbered registers */\
			);
		#elif defined(USE_AVX)	// AVX|AVX2
		__asm__ volatile (\
		"movslq	%[__nloop],%%rsi	\n\t"/* ASM loop structured as for(i = nloop; i != 0; i--){...} */\
			"movq %[__aptr],%%rax	\n\t	movq %[__bptr],%%rbx	\n\t	movq %[__cptr],%%rcx	\n\t"\
		"0:	\n\t"/* Note: Clang did not like longer label here, e.g. vdsub_loop in place of 0 gave "invalid symbol redefinition" error on __nloop */\
			"vmovaps	    (%%rax),%%ymm0	\n\t	vmovaps	0x20(%%rax),%%ymm1	\n\t"\
			"vmovaps	0x40(%%rax),%%ymm2	\n\t	vmovaps	0x60(%%rax),%%ymm3	\n\t"\
			"vmovaps	0x80(%%rax),%%ymm4	\n\t	vmovaps	0xa0(%%rax),%%ymm5	\n\t"\
			"vmovaps	0xc0(%%rax),%%ymm6	\n\t	vmovaps	0xe0(%%rax),%%ymm7	\n\t"\
			"vsubpd	    (%%rbx),%%ymm0,%%ymm0	\n\t	vsubpd	0x20(%%rbx),%%ymm1,%%ymm1	\n\t"\
			"vsubpd	0x40(%%rbx),%%ymm2,%%ymm2	\n\t	vsubpd	0x60(%%rbx),%%ymm3,%%ymm3	\n\t"\
			"vsubpd	0x80(%%rbx),%%ymm4,%%ymm4	\n\t	vsubpd	0xa0(%%rbx),%%ymm5,%%ymm5	\n\t"\
			"vsubpd	0xc0(%%rbx),%%ymm6,%%ymm6	\n\t	vsubpd	0xe0(%%rbx),%%ymm7,%%ymm7	\n\t"\
			"vmovaps	%%ymm0,    (%%rcx)	\n\t	vmovaps	%%ymm1,0x20(%%rcx)	\n\t"\
			"vmovaps	%%ymm2,0x40(%%rcx)	\n\t	vmovaps	%%ymm3,0x60(%%rcx)	\n\t"\
			"vmovaps	%%ymm4,0x80(%%rcx)	\n\t	vmovaps	%%ymm5,0xa0(%%rcx)	\n\t"\
			"vmovaps	%%ymm6,0xc0(%%rcx)	\n\t	vmovaps	%%ymm7,0xe0(%%rcx)	\n\t"\
			/* [a|b|c]ptr += SZ_VD: Double-pointer arithmetic, thus += SZ_VD*8 in raw bytewise form. AVX has SZ_VD = 32: */\
			"addq $0x100,%%rax	\n\t	addq $0x100,%%rbx	\n\t	addq $0x100,%%rcx	\n\t"\
		"decq	%%rsi			\n\t"\
		"jnz 0b	\n\t"/* loop end; continue is via jump-back if rsi != 0 */\
			:	/* outputs: none */\
			: [__nloop] "m" (nloop)	/* All inputs from memory addresses here */\
			 ,[__aptr] "m" (aptr)	\
			 ,[__bptr] "m" (bptr)	\
			 ,[__cptr] "m" (cptr)	\
			: "cc","memory","rax","rbx","rcx","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"	/* Clobbered registers */\
			);
		#else				// SSE2:
		__asm__ volatile (\
		"movslq	%[__nloop],%%rsi	\n\t"/* ASM loop structured as for(i = nloop; i != 0; i--){...} */\
			"movq %[__aptr],%%rax	\n\t	movq %[__bptr],%%rbx	\n\t	movq %[__cptr],%%rcx	\n\t"\
		"0:	\n\t"/* Note: Clang did not like longer label here, e.g. vdsub_loop in place of 0 gave "invalid symbol redefinition" error on __nloop */\
			"movaps	    (%%rax),%%xmm0	\n\t	movaps	0x10(%%rax),%%xmm1	\n\t"\
			"movaps	0x20(%%rax),%%xmm2	\n\t	movaps	0x30(%%rax),%%xmm3	\n\t"\
			"movaps	0x40(%%rax),%%xmm4	\n\t	movaps	0x50(%%rax),%%xmm5	\n\t"\
			"movaps	0x60(%%rax),%%xmm6	\n\t	movaps	0x70(%%rax),%%xmm7	\n\t"\
			"subpd	    (%%rbx),%%xmm0	\n\t	subpd	0x10(%%rbx),%%xmm1	\n\t"\
			"subpd	0x20(%%rbx),%%xmm2	\n\t	subpd	0x30(%%rbx),%%xmm3	\n\t"\
			"subpd	0x40(%%rbx),%%xmm4	\n\t	subpd	0x50(%%rbx),%%xmm5	\n\t"\
			"subpd	0x60(%%rbx),%%xmm6	\n\t	subpd	0x70(%%rbx),%%xmm7	\n\t"\
			"movaps	%%xmm0,    (%%rcx)	\n\t	movaps	%%xmm1,0x10(%%rcx)	\n\t"\
			"movaps	%%xmm2,0x20(%%rcx)	\n\t	movaps	%%xmm3,0x30(%%rcx)	\n\t"\
			"movaps	%%xmm4,0x40(%%rcx)	\n\t	movaps	%%xmm5,0x50(%%rcx)	\n\t"\
			"movaps	%%xmm6,0x60(%%rcx)	\n\t	movaps	%%xmm7,0x70(%%rcx)	\n\t"\
			/* [a|b|c]ptr += SZ_VD: Double-pointer arithmetic, thus += SZ_VD*8 in raw bytewise form. SSE2 has SZ_VD = 16: */\
			"addq $0x80,%%rax	\n\t	addq $0x80,%%rbx	\n\t	addq $0x80,%%rcx	\n\t"\
		"decq	%%rsi			\n\t"\
		"jnz 0b	\n\t"/* loop end; continue is via jump-back if rsi != 0 */\
			:	/* outputs: none */\
			: [__nloop] "m" (nloop)	/* All inputs from memory addresses here */\
			 ,[__aptr] "m" (aptr)	\
			 ,[__bptr] "m" (bptr)	\
			 ,[__cptr] "m" (cptr)	\
			: "cc","memory","rax","rbx","rcx","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"	/* Clobbered registers */\
			);
		#endif
	//		aptr += SZ_VD; bptr += SZ_VD; cptr += SZ_VD;
	//	}
	#endif
		for(i = (nloop << L2_SZ_VD); i < n; i++) {
			c[i] = a[i] - b[i];
		}
	#ifdef MULTITHREAD
		*(thread_arg->retval) = 0;	// 0 indicates successful return of current thread
	//	printf("Return from Thread %d ... ",thr_id);
		return 0x0;
	#endif
	}
#endif	// defined(PM1_STANDALONE)?

