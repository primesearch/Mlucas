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

#include "util.h"
#include "imul_macro.h"

/* Uncomment to enable Mp p-1 smoothness code - Note this is incompatible with standalone factor.c builds! */
//#define ENABLE_MPRIME_PM1_SMOOTH

#ifdef ENABLE_MPRIME_PM1_SMOOTH

	/* Tables of Fermat-base-2 pseudoprimes needed by the small-primes sieve code: */
	#include "f2psp_3_5.h"

	#undef psmooth
	struct psmooth
	{
		uint32 p;
		uint32 b;	// Standard B-smooth measure based on largest prime factor
		double r;	// L2 "roughness" metric in (0,1] defined by L2 norm of log factor sizes
	};

	// Decimal-print m(p) to a file in 100-digit chunks:
	void print_mp_dec()
	{
		const uint32 p = 13466917;
		const char fname[] = "p13466917_decimal.txt";
		FILE*fp = 0x0;
		uint32 i, lenX, lenD, nchars,nc, wrap_every = 100;	// Insert a LF every 100 digits
		uint64 *x,*y,*d,*r;
		uint64 ONES64 = 0xFFFFFFFFFFFFFFFFull;	// In GCC, making this 'const' gives "warning: overflow in implicit constant conversion" wherever it is used.
		char *str;

		// Allocate the array containing first M(p) and then subsequent divide-by-10^100 results.
		// Due to the requirement in mi64_div() that dividend and quotient arrays may not point
		// to the same memory, bounce successive-divide results between 2 arrays, x and y:
		lenX = (p>>6);
	//	x = (uint64 *)calloc(lenX + 1, sizeof(uint64));
		x = (uint64 *)calloc(((lenX + 3) & ~3), sizeof(uint64));	// Zero-pad to make multiple of 4, allowing 64-bit DIV algo to use 4-way-folded loops
		memset(x,ONES64,(lenX<<3));	x[lenX++] = (1ull << (p&63)) - 1;
		nchars = ceil(p * log(2.0)/log(10.));
		fprintf(stderr,"Generating decimal printout of M(%u), which has [%u] decimal digits; will write results to file '%s'...",p,nchars,fname);

	#if 0	// Hideously slow-but-reliable one-digit-at-a-time way:

		nchars += nchars/wrap_every + 1;
		str = (char *)calloc(nchars, sizeof(char));
		__convert_mi64_base10_char(str, nchars, x, lenX, wrap_every);

	#elif 0	// Try using fast right-to-left div-and-mod algo, modulo 10^100:

		nchars += nchars/100 + 1;
		str = (char *)calloc(nchars, sizeof(char));
		y = (uint64 *)calloc(lenX + 1, sizeof(uint64));
		// 10^100 has 333 bits, thus needs 6 uint64s, as do the mod-10^100 remainders,
		// but we allow the convert_base10_char_mi64() utility to do the allocation of the former for us:
		ASSERT(HERE, 0x0 != (d = convert_base10_char_mi64("10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000", &lenD)) && (lenD == 6), "0");
		r = (uint64 *)calloc(lenD, sizeof(uint64));
		nc = nchars-101;	// starting char of first 100 digit chunk
		for(i = 0; ; i+=2)	// i = #divides counter
		{
			mi64_div(x, d, lenX, lenD, y, r);	// dividend in y, remainder in r
			convert_mi64_base10_char_print_lead0(str + nc, r, lenD, 100,0);	nc -= 101;
			lenX = mi64_getlen(y, lenX);
			if( (lenX < lenD) || ((lenX == lenD) && mi64_cmpult(y,d,lenX)) ) {
				convert_mi64_base10_char(str, y, lenX, 100);
				break;
			}
			mi64_div(y, d, lenX, lenD, x, r);	// dividend in y, remainder in r
			convert_mi64_base10_char_print_lead0(str + nc, r, lenD, 100,0);	nc -= 101;
			lenX = mi64_getlen(x, lenX);
			if( (lenX < lenD) || ((lenX == lenD) && mi64_cmpult(x,d,lenX)) ) {
				convert_mi64_base10_char(str, x, lenX, 100);
				break;
			}
		}

	#else	// Same ideas as above but using modulo 10^27, the largest power of 10 whose odd factor (5^27) fits in a uint64,
			// thus allowing the core div-and-mod loops to use 1-word arguments:

		nchars += nchars/27 + 1;
		str = (char *)calloc(nchars, sizeof(char));
		y = (uint64 *)calloc(lenX + 1, sizeof(uint64));
		// 10^100 has 333 bits, thus needs 6 uint64s, as do the mod-10^100 remainders,
		// but we allow the convert_base10_char_mi64() utility to do the allocation of the former for us:
		ASSERT(HERE, 0x0 != (d = convert_base10_char_mi64("1000000000000000000000000000", &lenD)) && (lenD == 2), "0");
		r = (uint64 *)calloc(lenD, sizeof(uint64));
		nc = nchars- 28;	// starting char of first 27-digit chunk
		for(i = 0; ; i+=2)	// i = #divides counter
		{
			mi64_div(x, d, lenX, lenD, y, r);	// dividend in y, remainder in r
			convert_mi64_base10_char_print_lead0(str + nc, r, lenD, 27,0);	nc -= 28;	mi64_clear(r, lenD);
			lenX = mi64_getlen(y, lenX);
			if( (lenX < lenD) || ((lenX == lenD) && mi64_cmpult(y,d,lenX)) ) {
				convert_mi64_base10_char(str, y, lenX, 27);
				break;
			}
			mi64_div(y, d, lenX, lenD, x, r);	// dividend in y, remainder in r
			convert_mi64_base10_char_print_lead0(str + nc, r, lenD, 27,0);	nc -= 28;	mi64_clear(r, lenD);
			lenX = mi64_getlen(x, lenX);
			if( (lenX < lenD) || ((lenX == lenD) && mi64_cmpult(x,d,lenX)) ) {
				convert_mi64_base10_char(str, x, lenX, 27);
				break;
			}
		}

	#endif
		str[nchars-1] = '\0';

		fp = fopen(fname, "w");
		ASSERT(HERE, fp != 0x0, "Null file pointer!");
		fprintf(fp,"%s\n", str);
		fclose(fp);	fp = 0x0;
		exit(0);
	}

	// Binary predicates for use of stdlib qsort() on the b-subfield of the above psmooth struct:
	int psmooth_cmp_b(const void *x, const void *y)	// Default-int compare predicate
	{
		uint32 a = ((struct psmooth*)x)->b, b = ((struct psmooth*)y)->b;
		return ncmp_uint32( (void*)&a, (void*)&b );
	}

	// Binary predicates for use of stdlib qsort() on the r-subfield of the above psmooth struct:
	int psmooth_cmp_r(const void *x, const void *y)	// Default-int compare predicate
	{
		double two53float = (double)1.0*0x08000000*0x04000000;
		uint64 a = two53float*((struct psmooth*)x)->r, b = two53float*((struct psmooth*)y)->r;
		return ncmp_uint64( (void*)&a, (void*)&b );
	}

	void test_mp_pm1_smooth(uint32 p)
	{
		double u_so_smoove, logf, ilogn, dtmp;
		const double ln2 = log(2.0);
		uint32 nprime = 1000, pm_gap = 10000, thresh = 100000;
		uint32 curr_p,f2psp_idx,i,ihi,itmp32,j,jlo,jhi,k,max_diff,m,nfac,np,pm1;
		const uint32 pdiff_8[8] = {2,1,2,1,2,3,1,3}, pdsum_8[8] = { 0, 2, 6, 8,12,18,20,26};
		// Compact table storing the (difference/2) between adjacent odd primes.
		unsigned char *pdiff = (unsigned char *)calloc(nprime, sizeof(unsigned char));	// 1000 primes is plenty for this task
		// Struct used for storing smoothness data ... make big enough to store all primes in [p - pm_gap, p + pm_gap] with a safety factor
		struct psmooth sdat;
		// .../10 here is an approximation based on prime density for primes > 100000;
		// note the code uses an interval [p-pm_gap, p+pm_gap], i.e. of length 2*pm_gap, so the calloc needs to be twice pm_gap/10:
		struct psmooth*psmooth_vec = (struct psmooth *)calloc(2*pm_gap/10, sizeof(struct psmooth));

		/* Init first few diffs between 3/5, 5/7, 7/11, so can start loop with curr_p = 11 == 1 (mod 10), as required by twopmodq32_x8(): */
		pdiff[1] = 1;
		pdiff[2] = 1;
		ihi = curr_p = 11;
		/* Process chunks of length 30, starting with curr_p == 11 (mod 30). Applying the obvious divide-by-3,5 mini-sieve,
		we have 8 candidates in each interval: curr_p + [ 0, 2, 6, 8,12,18,20,26].
		For example: curr_p = 11 gives the 8 candidates: 11,13,17,19,23,29,31,37.
		*/
		f2psp_idx = 0;	// Index to next-expected Fermat base-2 pseudoprime in the precomputed table
		for(i = 3; i < nprime; curr_p += 30)
		{
			/* Make sure (curr_p + 29) < 2^32: */
			if(curr_p > 0xffffffe3)
			{
				fprintf(stderr,"curr_p overflows 32 bits!");
				nprime = i;
				break;
			}
			/* Do a quick Fermat base-2 compositeness test before invoking the more expensive mod operations: */
			itmp32 = twopmodq32_x8(curr_p, curr_p+ 2, curr_p+ 6, curr_p+ 8, curr_p+12, curr_p+18, curr_p+20, curr_p+26);
			for(j = 0; j < 8; ++j)
			{
				if((itmp32 >> j)&0x1)	// It's a PRP, so check against the table of known pseudoprimes and
				{						// (if it's not a PSP) init for the next gap
					ASSERT(HERE, curr_p <= f2psp[f2psp_idx],"Error in pseudoprime sieve");
					if((curr_p + pdsum_8[j]) == f2psp[f2psp_idx])	/* It's a base-2 pseudoprime */
					{
						++f2psp_idx;
						pdiff[i] += pdiff_8[j];
						continue;
					}
					else	/* It's prime - add final increment to current pdiff[i] and then increment i: */
					{
						ihi = (curr_p + pdsum_8[j]);
						pdiff[i] += pdiff_8[j];
						if(pdiff[i] > max_diff)
						{
							max_diff = pdiff[i];
						#if DBG_SIEVE
							printf("pdiff = %d at curr_p = %u\n", 2*max_diff,ihi);
						#endif
						}
						if(++i == nprime)
						{
							break;
						}
					}
				}
				else
				{
					pdiff[i] += pdiff_8[j];
				}
			}
			continue;
		}
		printf("Using first %u odd primes; max gap = %u\n",nprime,2*max_diff);
		printf("max sieving prime = %u\n",ihi);

		ASSERT(HERE, p > thresh, "Mersenne prime exponent must be larger that allowable threshold!");
		ASSERT(HERE, twopmodq32(p-1, p) == 1, "p fails base-2 fprp test!");
		np = 0;	// #primes in the current p-centered cohort
		// find N primes < and > p, compute smoothness norm based on p-1 factorization for each, store each [p,snorm] pair
		f2psp_idx = 0;	// Index to next-expected Fermat base-2 pseudoprime in the precomputed table
		jlo = p-pm_gap; jhi = p+pm_gap;
		// Find right tarting slot in base-2 pseudoprime table:
		while(f2psp[f2psp_idx] < jlo)
		{
			++f2psp_idx;
		}
		for(j = jlo; j <= jhi; j+=2)
		{
			// Do base-2 fprp test of j:
			if(!twopmodq32(j-1,j))
				continue;
			if(j == f2psp[f2psp_idx]) {	// It's a base-2 pseudoprime
				++f2psp_idx;
				continue;
			}
			// j is prime - compute factorization of j-1:
			sdat.p = j;
			pm1 = j - 1;
			printf("%u is prime: factorization of p-1 = ",j);
			ilogn = 1/log(1.0*pm1);	// 1/log(n)
			// We know 2 is a factor; special-case for that:
			nfac = 0;
			u_so_smoove = 0.0;
			curr_p = 2;
			logf = ln2;	// log(factor)
			while((pm1 & 1) == 0) {
				nfac++;
				pm1 >>= 1;
				dtmp = logf*ilogn;
				u_so_smoove += dtmp*dtmp;
			}
			if(nfac > 1) {
				printf("2^%u",nfac);
			} else {
				printf("2");
			}
			curr_p = 3;
			for(m = 0; m < nprime; m++)
			{
				if(pm1 < curr_p*curr_p)	{	// Remaining cofactor must be prime
					sdat.b = pm1;
					printf(".%u",pm1);
					nfac++;
					logf = log(1.0*pm1);	// log(factor)
					dtmp = logf*ilogn;
					u_so_smoove += dtmp*dtmp;
					break;
				}
				k = 0;	// factor multiplicity counter
				while((pm1 % curr_p) == 0) {// curr_p divides (p-1)
					nfac++;	k++;
					pm1 /= curr_p;
					logf = log(1.0*curr_p);	// log(factor)
					dtmp = logf*ilogn;
					u_so_smoove += dtmp*dtmp;
				}
				sdat.b = curr_p;
				if(k > 1) {
					printf(".%u^%u",curr_p,k);
				} else if(k == 1) {
					printf(".%u",curr_p);
				}
				if(pm1 == 1) break;
				curr_p += (pdiff[m] << 1);
			}
			// L2 norm: divide by #factors (multiple-counting repeated factors):
			u_so_smoove = sqrt(u_so_smoove)/nfac;
			sdat.r = u_so_smoove;
			psmooth_vec[np++] = sdat;	// Write completed datum to array or later sorting
			printf("; %u factors, L2 smoothness = %15.13f\n",nfac,u_so_smoove);
		}	// for(j in [p +- pm_gap] loop
		printf("\n");

		// Using array of [p,snorm]-pair structs, sort resulting array-aof-structs by snorm value:
		qsort(psmooth_vec, np, sizeof(struct psmooth), psmooth_cmp_b);
		for(j = 0; j < np; j++) {
			sdat = psmooth_vec[j];
		//	printf("p = %u: B -smoothness = %u\n",sdat.p,sdat.b);
			if(sdat.p == p) {
				printf("B -smoothness: %u is %u of %u, percentile = %5.2f\n",p,j+1,np,100.0*((double)np-j)/np);
				break;
			}
		}
		qsort(psmooth_vec, np, sizeof(struct psmooth), psmooth_cmp_r);
		for(j = 0; j < np; j++) {
			sdat = psmooth_vec[j];
		//	printf("p = %u: L2 smoothness = %15.13f\n",sdat.p,sdat.r);
			if(sdat.p == p) {
				printf("L2-smoothness: %u is %u of %u, percentile = %5.2f\n",p,j+1,np,100.0*((double)np-j)/np);
				break;
			}
		}
		exit(0);	// Only enable this bit of test code when a new M-prime is found
	}	// test_mp_pm1_smooth()

#endif	// ENABLE_MPRIME_PM1_SMOOTH

#undef X64_ASM
#if(defined(CPU_IS_X86_64) && defined(COMPILER_TYPE_GCC) && (OS_BITS == 64))
	#warning util.c: Defining X64_ASM
	#define X64_ASM
	#if FP_MANTISSA_BITS_DOUBLE != 64
		#error x86_64 asm requires FP_MANTISSA_BITS_DOUBLE == 64!
	#endif
#endif

#if defined(USE_GPU) && defined(__CUDACC__)

	#include "gpu_iface.h"

	// Simple vector-add test function:
	__global__ void VecAdd(float* A, float* B, float* C, int N)
	{
		int i = blockDim.x * blockIdx.x + threadIdx.x;
		// Uncomment the if() to Print basic info about threads ... keep I/O reasonable by only
		// doing so for the first 10 of each batch of 2^18:
	//
		if(i%0x3ffff < 10) {
			printf("GPU block %d[dim %d], thread %d ==> seq-thread %d [i%0x3ffff = %d]... \n", blockIdx.x, blockDim.x, threadIdx.x, i,i%0x3ffff);
		}
	//
		if (i < N)
			C[i] = A[i] + B[i];
		else
			printf("GPU block %d[dim %d], thread %d: ERROR: I = %d out of range!\n", blockIdx.x, blockDim.x, threadIdx.x, i);
	}

	// Host code for the VecAdd test:
	void cudaVecAddTest()
	{
		int i, N = 1024*1024;
		size_t size = N * sizeof(float);
		// Allocate input vectors h_A and h_B in host memory
		float* h_A = (float*)malloc(size);
		float* h_B = (float*)malloc(size);
		float* h_C = (float*)malloc(size);
		// Initialize input vectors
		for(i = 0; i < N; ++i) {
			*(h_A+i) = i;
			*(h_B+i) = i*0.1;
		}
		// Allocate vectors in device memory
		float* d_A;
		cudaMalloc(&d_A, size);
		float* d_B;
		cudaMalloc(&d_B, size);
		float* d_C;
		cudaMalloc(&d_C, size);
		// Copy vectors from host memory to device memory
		cudaMemcpy(d_A, h_A, size, cudaMemcpyHostToDevice);
		cudaMemcpy(d_B, h_B, size, cudaMemcpyHostToDevice);
		// Invoke kernel
		int threadsPerBlock = 256;
		int blocksPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock;
		VecAdd<<<blocksPerGrid, threadsPerBlock>>>(d_A, d_B, d_C, N);
		// Copy result from device memory to host memory
		// h_C contains the result in host memory
		cudaMemcpy(h_C, d_C, size, cudaMemcpyDeviceToHost);
		// Free device memory
		cudaFree(d_A);
		cudaFree(d_B);
		cudaFree(d_C);
		// Debug-print results sample:
		for(i = 0; i < 10; ++i) {
			printf("i = %d: Sum = %10.2f\n", i, *(h_C+i));
		}
		printf("...\n");
		for(i = 10; i > 0; --i) {
			printf("i = %d: Sum = %10.2f\n", N-i, *(h_C+N-i));
		}
	}

	/**********************************************************************************************/

	#include "factor.h"
	#include "fac_test_dat96.h"
	#include "twopmodq80.h"

	// Host code for the simpler VecModpow test, same 78-bit [p,q] pair for each thread:
	void cudaVecModpowTest0()
	{
		int i;
		uint32 p, pshift, start_index, zshift, j, lead7;
		uint32 N = 1<<10;
		uint64 k;
		uint96 q96;
		double dbl, rnd;

		// Allocate input vectors (which take the TF p/pshift/zshift/start_index/k data on input) in host memory:
		uint32* h_p     = (uint32*)malloc(N<<2);
		uint32* h_pshft = (uint32*)malloc(N<<2);
		uint32* h_zshft = (uint32*)malloc(N<<2);
		uint32* h_stidx = (uint32*)malloc(N<<2);
		uint64* h_k     = (uint64*)malloc(N<<3);

		p = 16727479;
		k = 7946076362870052;
		// Compute auxiliary TF data:
		pshift = p + 78;
		j = leadz32(pshift);
		/* Extract leftmost 7 bits of pshift (if > 77, use the leftmost 6) and subtract from 96: */
		lead7 = ((pshift<<j) >> 25);
		if(lead7 > 77) {
			lead7 >>= 1;
			start_index =  32-j-6;	/* Use only the leftmost 6 bits */
		} else {
			start_index =  32-j-7;
		}
		zshift = 77 - lead7;
		zshift <<= 1;				/* Doubling the shift count here takes cares of the first SQR_LOHI */
		pshift = ~pshift;

		// Copy to all N vector-input-data:
		for(i = 0; i < N; ++i) {
			*(h_p     + i) = p          ;
			*(h_pshft + i) = pshift     ;
			*(h_zshft + i) = zshift     ;
			*(h_stidx + i) = start_index;
			*(h_k     + i) = k          ;
		}
		printf("Testing %d 78-bit known-factors:\n",N);

		// Initialize output vector (resulting 2^p mod q, in binary "is factor?" form) in host memory:
		uint8*  h_B = (uint8 *)malloc(N);	// Until impl packed-bitmap scheme for device code return values, use byte array for return values
		for(i = 0; i < N; ++i) {
			*(h_B+i) = 0;
		}
//printf("Host code: p = %u, pshift = %u, k = %llu, zshift = %u, start_index = %u\n", p,pshift,h_A[0],zshift,start_index);
		// Allocate vectors in device memory
		uint32 *d_p,*d_pshft,*d_zshft,*d_stidx;
		uint64* d_k;
		cudaMalloc(&d_p    , N<<2);
		cudaMalloc(&d_pshft, N<<2);
		cudaMalloc(&d_zshft, N<<2);
		cudaMalloc(&d_stidx, N<<2);
		cudaMalloc(&d_k    , N<<3);
		uint8 * d_B;
		cudaMalloc(&d_B, N);
		// Copy vectors from host memory to device memory
		cudaMemcpy(d_p    , h_p    , N<<2, cudaMemcpyHostToDevice);
		cudaMemcpy(d_pshft, h_pshft, N<<2, cudaMemcpyHostToDevice);
		cudaMemcpy(d_zshft, h_zshft, N<<2, cudaMemcpyHostToDevice);
		cudaMemcpy(d_stidx, h_stidx, N<<2, cudaMemcpyHostToDevice);
		cudaMemcpy(d_k    , h_k    , N<<3, cudaMemcpyHostToDevice);
		// Do we need to copy the as-yet-uninited output vector to (or just from) the device?
		cudaMemcpy(d_B, h_B, N   , cudaMemcpyHostToDevice);

		// Invoke kernel
		int threadsPerBlock = 256;
		int blocksPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock;

		printf("VecModpow0 with %d 78-bit known-factors [blocksPerGrid = %d, threadsPerBlock = %d]:\n",N, blocksPerGrid, threadsPerBlock);
		VecModpow<<<blocksPerGrid, threadsPerBlock>>>(d_p,d_pshft,d_zshft,d_stidx,d_k, d_B, N);

		// Copy result from device memory to host memory
		// h_B contains the result in host memory
		cudaMemcpy(h_B, d_B, N, cudaMemcpyDeviceToHost);
		// Free device memory
		cudaFree(d_p    );
		cudaFree(d_pshft);
		cudaFree(d_zshft);
		cudaFree(d_stidx);
		cudaFree(d_k    );
		cudaFree(d_B);

		// Reference computation:
		uint64 checksum1, checksum2;
		j = (uint32)twopmodq78_3WORD_DOUBLE(&checksum1, &checksum2, (uint64)p, k);
		ASSERT(HERE, (j == 1), "cudaVecModpowTest0 ref-comp failed!");
		// Test GPU results:
		for(i = 0; i < N; ++i) {
			if(*(h_B + i) != 1) {
				printf("cudaVecModpowTest: Mismatch between Ref and GPU result:\n");
				printf("res[%d] = %d [ref = %d] = 2^p - 1 (mod q) with (p,k) = %u, %ull\n", i,*(h_B + i), j,p,k);
				ASSERT(HERE, *(h_B + i) == 1, "cudaVecModpowTest failed!");
			}
		}
		printf("cudaVecModpowTest with %d test (p,q) pairs succeeded!\n",N);
	}

	// Host code for the VecModpow test:
	void cudaVecModpowTest()
	{
		int i;
		uint32 p, pshift, start_index, zshift, j, lead7;
		uint32 N = 1<<10,nelts;
		uint64 k;
		uint96 q96;
		double dbl, rnd;

		// Allocate input vectors (which take the TF p/pshift/zshift/start_index/k data on input) in host memory:
		uint32* h_p     = (uint32*)malloc(N<<2);
		uint32* h_pshft = (uint32*)malloc(N<<2);
		uint32* h_zshft = (uint32*)malloc(N<<2);
		uint32* h_stidx = (uint32*)malloc(N<<2);
		uint64* h_k     = (uint64*)malloc(N<<3);
		for(i = 0, nelts = 0; i < N; ++i) {
			p = fac96[i].p;
			if(p == 0) {
				break;
			}
			q96.d1 = fac96[i].d1; q96.d0 = fac96[i].d0;
			if((q96.d1 >> 14) != 0) {
				continue;
			}
			// Good to go - compute auxiliary TF data:
			pshift = p + 78;
			j = leadz32(pshift);
			/* Extract leftmost 7 bits of pshift (if > 77, use the leftmost 6) and subtract from 96: */
			lead7 = ((pshift<<j) >> 25);
			if(lead7 > 77) {
				lead7 >>= 1;
				start_index =  32-j-6;	/* Use only the leftmost 6 bits */
			} else {
				start_index =  32-j-7;
			}
			zshift = 77 - lead7;
			zshift <<= 1;				/* Doubling the shift count here takes cares of the first SQR_LOHI */
			pshift = ~pshift;

			// Compute factor k using fast DP math. Integer-truncation-on-store should obviate the need
			// to subtract 1 from q, and (double)q is only accurate to 53 bits to begin with):
			dbl = (double)q96.d0 + (double)q96.d1*TWO64FLOAT;
			dbl /= (2.0*p);
			rnd = DNINT(dbl);
			k = (uint64)rnd;

			*(h_p     + nelts) = p          ;
			*(h_pshft + nelts) = pshift     ;
			*(h_zshft + nelts) = zshift     ;
			*(h_stidx + nelts) = start_index;
			*(h_k     + nelts) = k          ;
//	printf("p[%3d] = %u: pshift = %8u, zshift = %8u, stidx = %2u, k = %llu\n",nelts, p, pshift, zshift, start_index, k);
			++nelts;
		}
		printf("Testing %d 78-bit known-factors:\n",nelts);
		// "Fill in" remaining slots with copy of same datum used in cudaVecModpowTest0:
		p = 16727479;
		k = 7946076362870052;
		// Compute auxiliary TF data:
		pshift = p + 78;
		j = leadz32(pshift);
		/* Extract leftmost 7 bits of pshift (if > 77, use the leftmost 6) and subtract from 96: */
		lead7 = ((pshift<<j) >> 25);
		if(lead7 > 77) {
			lead7 >>= 1;
			start_index =  32-j-6;	/* Use only the leftmost 6 bits */
		} else {
			start_index =  32-j-7;
		}
		zshift = 77 - lead7;
		zshift <<= 1;				/* Doubling the shift count here takes cares of the first SQR_LOHI */
		pshift = ~pshift;
		// Copy to all still-uninited vector-input-data:
		for(i = nelts; i < N; ++i) {
			*(h_p     + i) = p          ;
			*(h_pshft + i) = pshift     ;
			*(h_zshft + i) = zshift     ;
			*(h_stidx + i) = start_index;
			*(h_k     + i) = k          ;
		}

		// Initialize output vector (resulting 2^p mod q, in binary "is factor?" form) in host memory:
		uint8*  h_B = (uint8 *)malloc(N);	// Until impl packed-bitmap scheme for device code return values, use byte array for return values
		for(i = 0; i < N; ++i) {
			*(h_B+i) = 0;
		}
//printf("Host code: p = %u, pshift = %u, k = %llu, zshift = %u, start_index = %u\n", p,pshift,h_A[0],zshift,start_index);
		// Allocate vectors in device memory
		uint32 *d_p,*d_pshft,*d_zshft,*d_stidx;
		uint64* d_k;
		cudaMalloc(&d_p    , N<<2);
		cudaMalloc(&d_pshft, N<<2);
		cudaMalloc(&d_zshft, N<<2);
		cudaMalloc(&d_stidx, N<<2);
		cudaMalloc(&d_k    , N<<3);
		uint8 * d_B;
		cudaMalloc(&d_B, N);
		// Copy vectors from host memory to device memory
		cudaMemcpy(d_p    , h_p    , N<<2, cudaMemcpyHostToDevice);
		cudaMemcpy(d_pshft, h_pshft, N<<2, cudaMemcpyHostToDevice);
		cudaMemcpy(d_zshft, h_zshft, N<<2, cudaMemcpyHostToDevice);
		cudaMemcpy(d_stidx, h_stidx, N<<2, cudaMemcpyHostToDevice);
		cudaMemcpy(d_k    , h_k    , N<<3, cudaMemcpyHostToDevice);
		// Do we need to copy the as-yet-uninited output vector to (or just from) the device?
		cudaMemcpy(d_B, h_B, N   , cudaMemcpyHostToDevice);

		// Invoke kernel
		int threadsPerBlock = 256;
		int blocksPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock;

		printf("VecModpow with %d 78-bit known-factors [blocksPerGrid = %d, threadsPerBlock = %d]:\n",nelts, blocksPerGrid, threadsPerBlock);
		VecModpow<<<blocksPerGrid, threadsPerBlock>>>(d_p,d_pshft,d_zshft,d_stidx,d_k, d_B, nelts);

		// Copy result from device memory to host memory
		// h_B contains the result in host memory
		cudaMemcpy(h_B, d_B, N, cudaMemcpyDeviceToHost);
		// Free device memory
		cudaFree(d_p    );
		cudaFree(d_pshft);
		cudaFree(d_zshft);
		cudaFree(d_stidx);
		cudaFree(d_k    );
		cudaFree(d_B);

		// Reference computation:
		uint64 checksum1, checksum2;
		// Test GPU results:
		for(i = 0; i < nelts; ++i) {
			p = *(h_p + i);
			k = *(h_k + i);
			j = (uint32)twopmodq78_3WORD_DOUBLE(&checksum1, &checksum2, (uint64)p, k);
			if((j != 1) || (*(h_B + i) != 1)) {
				printf("cudaVecModpowTest: Mismatch between Ref and GPU result:\n");
				printf("res[%d] = %d [ref = %d] = 2^p - 1 (mod q) with (p,k) = %u, %ull\n", i,*(h_B + i), j,p,k);
				ASSERT(HERE, (j == 1) && (*(h_B + i) == 1), "cudaVecModpowTest failed!");
			}
		}
		printf("cudaVecModpowTest with %d test (p,q) pairs succeeded!\n",nelts);
	}

#endif

#undef PLATFORM_SKIP_RND_CONST_ENFORCEMENT

/*********************************************************************************/
/* Globals. Unless specified otherwise, these are declared in Mdata.h:           */
/*********************************************************************************/

/* These externs yield 64-mantissa-bit register mode (bits <9:8> = 3),
with IEEE and truncating rounding mode set via bits <11:10> = 0 and 3,
respectively. The other 12 bits are identical to the MSVC defaults: */
unsigned short FPU_64RND = 0x037f, FPU_64CHOP = 0x0f7f;

const int CHAROFFSET = '0';

double RND_A, RND_B;	/* Used for fast NINT emulation; set in util.c. */

double             TWO13FLINV;	/* (double)2^13 inverse */
double TWO25FLOAT, TWO25FLINV;	/* (double)2^25 and inverse */
double TWO26FLOAT, TWO26FLINV;	/* (double)2^26 and inverse */
double TWO32FLOAT, TWO32FLINV;	/* (double)2^32 and inverse */
double TWO48FLOAT, TWO48FLINV;	/* (double)2^48 and inverse */
double TWO50FLOAT, TWO50FLINV;	/* (double)2^50 and inverse */
double TWO51FLOAT, TWO51FLINV;	/* (double)2^51 and inverse */
double TWO52FLOAT, TWO52FLINV;	/* (double)2^52 and inverse */
double TWO53FLOAT, TWO53FLINV;	/* (double)2^53 and inverse */
double TWO54FLOAT;	/* (double)2^54 */
double TWO63FLOAT;	/* (double)2^63 */
double TWO64FLOAT, TWO64FLINV;	/* (double)2^64 and inverse */

int32 DAT_BITS, PAD_BITS;	/* Array padding parameters */

/* Fixed-size (but only necessarily constant during a given FFT-based MUL)
   base for generic FFT-based mul:

	FFT_MUL_BASE = 2^(FFT_MUL_BITS), where FFT_MUL_BITS #def'ed in Mdata.h
*/
double FFT_MUL_BASE, FFT_MUL_BASE_INV;

/***********************/

/*** 11/23/05: MSVC/.NET buggered things up with the second of these tables
     when each table was local to its respective calling function, so moved 'em here: ***/

	/* Table of approximate byte-inverses of 1.{byteval} is here. Since we know
	the input is in [1, 2), we know the multiplicative inverse is in (0.5, 1],
	i.e. we know that the MSB of the inverse (the one immediately right of the
	binary point) is 1. Thus we can use the hidden-bit-is-1 property of inputs
	to also gain a bit of precision in the bytewise approximate inverses, by
	neglecting the leading-order bit - since that one would get stored in the
	hidden-bit slot of the output anyway, this also makes our work easier. */
	/* Unix bc code:
	bc -l
	ibase=2
	obase=2
	d=0.00000001;
	x=1.000000001-d;
	x+=d;1/x
	{256 of these, round 10th bit into MS9, replace MSB by '0x', convert rest to hex}
	*/
	static const uint8 byte_lookup_finvest[256] = {
	0xff,0xfd,0xfb,0xf9,0xf7,0xf5,0xf3,0xf1,0xf0,0xee,0xec,0xea,0xe8,0xe6,0xe5,0xe3,
	0xe1,0xdf,0xdd,0xdc,0xda,0xd8,0xd7,0xd5,0xd3,0xd2,0xd0,0xce,0xcd,0xcb,0xc9,0xc8,
	0xc6,0xc5,0xc3,0xc2,0xc0,0xbf,0xbd,0xbc,0xba,0xb9,0xb7,0xb6,0xb4,0xb3,0xb1,0xb0,
	0xae,0xad,0xac,0xaa,0xa9,0xa7,0xa6,0xa5,0xa3,0xa2,0xa1,0x9f,0x9e,0x9d,0x9c,0x9a,
	0x99,0x98,0x96,0x95,0x94,0x93,0x91,0x90,0x8f,0x8e,0x8d,0x8b,0x8a,0x89,0x88,0x87,
	0x86,0x84,0x83,0x82,0x81,0x80,0x7f,0x7e,0x7c,0x7b,0x7a,0x79,0x78,0x77,0x76,0x75,
	0x74,0x73,0x72,0x71,0x70,0x6f,0x6e,0x6d,0x6c,0x6b,0x6a,0x69,0x68,0x67,0x66,0x65,
	0x64,0x63,0x62,0x61,0x60,0x5f,0x5e,0x5d,0x5c,0x5b,0x5a,0x59,0x58,0x58,0x57,0x56,
	0x55,0x54,0x53,0x52,0x51,0x51,0x50,0x4f,0x4e,0x4d,0x4c,0x4b,0x4b,0x4a,0x49,0x48,
	0x47,0x46,0x46,0x45,0x44,0x43,0x42,0x42,0x41,0x40,0x3f,0x3f,0x3e,0x3d,0x3c,0x3b,
	0x3b,0x3a,0x39,0x38,0x38,0x37,0x36,0x35,0x35,0x34,0x33,0x33,0x32,0x31,0x30,0x30,
	0x2f,0x2e,0x2e,0x2d,0x2c,0x2c,0x2b,0x2a,0x2a,0x29,0x28,0x28,0x27,0x26,0x26,0x25,
	0x24,0x24,0x23,0x22,0x22,0x21,0x20,0x20,0x1f,0x1e,0x1e,0x1d,0x1d,0x1c,0x1b,0x1b,
	0x1a,0x1a,0x19,0x18,0x18,0x17,0x17,0x16,0x15,0x15,0x14,0x14,0x13,0x12,0x12,0x11,
	0x11,0x10,0x10,0x0f,0x0f,0x0e,0x0d,0x0d,0x0c,0x0c,0x0b,0x0b,0x0a,0x0a,0x09,0x09,
	0x08,0x07,0x07,0x06,0x06,0x05,0x05,0x04,0x04,0x03,0x03,0x02,0x02,0x01,0x01,0x00
	};

/***********************/

	/* Table of approximate byte-inverses of 2 * 1.{byteval} is here. Since we know
	the input is in [1, 4), we know the inverse-square-rootis in (0.5, 1],
	i.e. we know that the MSB of the ISQRT (the one immediately right of the
	binary point) is 1.	We cheat a little on the 0 element of the byte table,
	since sqrt(1.000000001) really should give 0x100, not 0xff. But the
	alternative is using uint16s, which doubles the size of the table. */
	/* Unix bc code:
	bc -l
	ibase=2
	obase=2
	d=0.00000001;
	x=1.000000001-d;
	x+=d;1/sqrt(x)
	{768 of these, round 10th bit into MS9, replace MSB by '0x', convert rest to hex}
	*/
	/* Used to store MS 8 non-hidden mantissa bits. We'd need to use a 16-bit int
	to allow for the possibility of a carryout (i.e. result = 256) from rounding
	the 9th-most-significant NHB into the upper 8 (which would involve
	additional logic to handle), we instead deal with the issue of rounding
	by assuming the midpoint - e.g. if truncating to the MS 8 NHBs yields
	a certain integer in [0,255], we assume the resulting roundoff error
	is always 0.5, i.e. our precomputed 1/x values are approximations to
	the resulting midpoints. This also avoids our having to treat an input
	of 1.00000000 as a special case, since we munge that to 1.000000001,
	whose inverse is < 1.0: */
	static const uint8 byte_lookup_fisqrtest[768] = {
	0xff,0xff,0xfe,0xfd,0xfc,0xfb,0xfa,0xf9,0xf8,0xf7,0xf6,0xf5,0xf4,0xf3,0xf2,0xf1,
	0xf0,0xef,0xee,0xee,0xed,0xec,0xeb,0xea,0xe9,0xe8,0xe7,0xe7,0xe6,0xe5,0xe4,0xe3,
	0xe2,0xe1,0xe1,0xe0,0xdf,0xde,0xdd,0xdd,0xdc,0xdb,0xda,0xd9,0xd9,0xd8,0xd7,0xd6,
	0xd5,0xd5,0xd4,0xd3,0xd2,0xd2,0xd1,0xd0,0xcf,0xcf,0xce,0xcd,0xcc,0xcc,0xcb,0xca,
	0xca,0xc9,0xc8,0xc7,0xc7,0xc6,0xc5,0xc5,0xc4,0xc3,0xc3,0xc2,0xc1,0xc1,0xc0,0xbf,
	0xbf,0xbe,0xbd,0xbd,0xbc,0xbb,0xbb,0xba,0xb9,0xb9,0xb8,0xb7,0xb7,0xb6,0xb6,0xb5,
	0xb4,0xb4,0xb3,0xb2,0xb2,0xb1,0xb1,0xb0,0xaf,0xaf,0xae,0xae,0xad,0xac,0xac,0xab,
	0xab,0xaa,0xaa,0xa9,0xa8,0xa8,0xa7,0xa7,0xa6,0xa6,0xa5,0xa5,0xa4,0xa3,0xa3,0xa2,
	0xa2,0xa1,0xa1,0xa0,0xa0,0x9f,0x9f,0x9e,0x9d,0x9d,0x9c,0x9c,0x9b,0x9b,0x9a,0x9a,
	0x99,0x99,0x98,0x98,0x97,0x97,0x96,0x96,0x95,0x95,0x94,0x94,0x93,0x93,0x92,0x92,
	0x91,0x91,0x90,0x90,0x8f,0x8f,0x8f,0x8e,0x8e,0x8d,0x8d,0x8c,0x8c,0x8b,0x8b,0x8a,
	0x8a,0x89,0x89,0x89,0x88,0x88,0x87,0x87,0x86,0x86,0x85,0x85,0x85,0x84,0x84,0x83,
	0x83,0x82,0x82,0x82,0x81,0x81,0x80,0x80,0x7f,0x7f,0x7f,0x7e,0x7e,0x7d,0x7d,0x7d,
	0x7c,0x7c,0x7b,0x7b,0x7a,0x7a,0x7a,0x79,0x79,0x78,0x78,0x78,0x77,0x77,0x76,0x76,
	0x76,0x75,0x75,0x75,0x74,0x74,0x73,0x73,0x73,0x72,0x72,0x72,0x71,0x71,0x70,0x70,
	0x70,0x6f,0x6f,0x6f,0x6e,0x6e,0x6d,0x6d,0x6d,0x6c,0x6c,0x6c,0x6b,0x6b,0x6b,0x6a,
	0x6a,0x6a,0x69,0x69,0x68,0x68,0x68,0x67,0x67,0x67,0x66,0x66,0x66,0x65,0x65,0x65,
	0x64,0x64,0x64,0x63,0x63,0x63,0x62,0x62,0x62,0x61,0x61,0x61,0x60,0x60,0x60,0x5f,
	0x5f,0x5f,0x5e,0x5e,0x5e,0x5d,0x5d,0x5d,0x5d,0x5c,0x5c,0x5c,0x5b,0x5b,0x5b,0x5a,
	0x5a,0x5a,0x59,0x59,0x59,0x58,0x58,0x58,0x58,0x57,0x57,0x57,0x56,0x56,0x56,0x55,
	0x55,0x55,0x55,0x54,0x54,0x54,0x53,0x53,0x53,0x53,0x52,0x52,0x52,0x51,0x51,0x51,
	0x51,0x50,0x50,0x50,0x4f,0x4f,0x4f,0x4f,0x4e,0x4e,0x4e,0x4d,0x4d,0x4d,0x4d,0x4c,
	0x4c,0x4c,0x4c,0x4b,0x4b,0x4b,0x4a,0x4a,0x4a,0x4a,0x49,0x49,0x49,0x49,0x48,0x48,
	0x48,0x48,0x47,0x47,0x47,0x47,0x46,0x46,0x46,0x45,0x45,0x45,0x45,0x44,0x44,0x44,
	0x44,0x43,0x43,0x43,0x43,0x42,0x42,0x42,0x42,0x41,0x41,0x41,0x41,0x40,0x40,0x40,
	0x40,0x3f,0x3f,0x3f,0x3f,0x3f,0x3e,0x3e,0x3e,0x3e,0x3d,0x3d,0x3d,0x3d,0x3c,0x3c,
	0x3c,0x3c,0x3b,0x3b,0x3b,0x3b,0x3a,0x3a,0x3a,0x3a,0x3a,0x39,0x39,0x39,0x39,0x38,
	0x38,0x38,0x38,0x38,0x37,0x37,0x37,0x37,0x36,0x36,0x36,0x36,0x36,0x35,0x35,0x35,
	0x35,0x34,0x34,0x34,0x34,0x34,0x33,0x33,0x33,0x33,0x32,0x32,0x32,0x32,0x32,0x31,
	0x31,0x31,0x31,0x31,0x30,0x30,0x30,0x30,0x30,0x2f,0x2f,0x2f,0x2f,0x2e,0x2e,0x2e,
	0x2e,0x2e,0x2d,0x2d,0x2d,0x2d,0x2d,0x2c,0x2c,0x2c,0x2c,0x2c,0x2b,0x2b,0x2b,0x2b,
	0x2b,0x2a,0x2a,0x2a,0x2a,0x2a,0x29,0x29,0x29,0x29,0x29,0x28,0x28,0x28,0x28,0x28,
	0x28,0x27,0x27,0x27,0x27,0x27,0x26,0x26,0x26,0x26,0x26,0x25,0x25,0x25,0x25,0x25,
	0x24,0x24,0x24,0x24,0x24,0x24,0x23,0x23,0x23,0x23,0x23,0x22,0x22,0x22,0x22,0x22,
	0x22,0x21,0x21,0x21,0x21,0x21,0x20,0x20,0x20,0x20,0x20,0x20,0x1f,0x1f,0x1f,0x1f,
	0x1f,0x1f,0x1e,0x1e,0x1e,0x1e,0x1e,0x1d,0x1d,0x1d,0x1d,0x1d,0x1d,0x1c,0x1c,0x1c,
	0x1c,0x1c,0x1c,0x1b,0x1b,0x1b,0x1b,0x1b,0x1b,0x1a,0x1a,0x1a,0x1a,0x1a,0x1a,0x19,
	0x19,0x19,0x19,0x19,0x19,0x18,0x18,0x18,0x18,0x18,0x18,0x17,0x17,0x17,0x17,0x17,
	0x17,0x16,0x16,0x16,0x16,0x16,0x16,0x15,0x15,0x15,0x15,0x15,0x15,0x15,0x14,0x14,
	0x14,0x14,0x14,0x14,0x13,0x13,0x13,0x13,0x13,0x13,0x13,0x12,0x12,0x12,0x12,0x12,
	0x12,0x11,0x11,0x11,0x11,0x11,0x11,0x11,0x10,0x10,0x10,0x10,0x10,0x10,0x0f,0x0f,
	0x0f,0x0f,0x0f,0x0f,0x0f,0x0e,0x0e,0x0e,0x0e,0x0e,0x0e,0x0e,0x0d,0x0d,0x0d,0x0d,
	0x0d,0x0d,0x0d,0x0c,0x0c,0x0c,0x0c,0x0c,0x0c,0x0c,0x0b,0x0b,0x0b,0x0b,0x0b,0x0b,
	0x0b,0x0a,0x0a,0x0a,0x0a,0x0a,0x0a,0x0a,0x09,0x09,0x09,0x09,0x09,0x09,0x09,0x08,
	0x08,0x08,0x08,0x08,0x08,0x08,0x08,0x07,0x07,0x07,0x07,0x07,0x07,0x07,0x06,0x06,
	0x06,0x06,0x06,0x06,0x06,0x05,0x05,0x05,0x05,0x05,0x05,0x05,0x05,0x04,0x04,0x04,
	0x04,0x04,0x04,0x04,0x04,0x03,0x03,0x03,0x03,0x03,0x03,0x03,0x02,0x02,0x02,0x02,
	0x02,0x02,0x02,0x02,0x01,0x01,0x01,0x01,0x01,0x01,0x01,0x01,0x00,0x00,0x00,0x00
	};


/********FILE-ACCESS-RELATED:**********/

int		file_valid(FILE*fp)
{
	if(fp && !ferror(fp) && !feof(fp))
		return TRUE;
	else
		return FALSE;
}

/* Apparently each implementation uses its own names for these file-internals flags -
TODO: need a portable way of checking file-access mode here.
int		file_valid_for_read (FILE*fp)
{
	blah... (similar to below)
}

int		file_valid_for_write(FILE*fp)
{
	if(fp && (fp->flag & (_EOF | _ERR) == 0) && (fp->flag & _WRITE != 0) )
		return TRUE;
	else
		return FALSE;
}
*/
/******************/

/* Print key platform info, (on x86) set FPU mode, do some basic self-tests: */
void host_init(void)
{
#ifdef MULTITHREAD
	int ncpu, nthr;
#endif
	double dbl;

	/* Various useful precomputed powers of 2 in floating-double form: */
	TWO25FLOAT = (double)0x02000000;
	TWO25FLINV = 1.0/TWO25FLOAT;

	TWO26FLOAT = (double)0x04000000;
	TWO26FLINV = 1.0/TWO26FLOAT;
	dbl = qfdbl(qfmul_pow2(QONE, -26));
	ASSERT(HERE, TWO26FLINV == dbl, "TWO26FLINV!");
	TWO13FLINV = qfdbl(qfmul_pow2(QONE, -13));

	TWO32FLOAT = (double)2.0*0x80000000;
	TWO32FLINV = 1.0/TWO32FLOAT;

	TWO48FLOAT = (double)1.0*0x01000000*0x01000000;
	TWO48FLINV = 1.0/TWO48FLOAT;

	TWO50FLOAT = (double)1.0*0x01000000*0x04000000;
	TWO50FLINV = 1.0/TWO50FLOAT;

	TWO51FLOAT = (double)1.0*0x02000000*0x04000000;
	TWO51FLINV = 1.0/TWO51FLOAT;

	TWO52FLOAT = (double)1.0*0x04000000*0x04000000;
	TWO52FLINV = 1.0/TWO52FLOAT;

	TWO53FLOAT = (double)1.0*0x08000000*0x04000000;
	TWO53FLINV = 1.0/TWO53FLOAT;

	TWO54FLOAT = (double)1.0*0x08000000*0x08000000;	/* (double)2^54 */

	TWO63FLOAT = (double)2.0*0x80000000*0x80000000;

	TWO64FLOAT = (double)4.0*0x80000000*0x80000000;
	TWO64FLINV = 1.0/TWO64FLOAT;

	/* Check qfloat routines (this call is also needed to init various qfloat global constants): */
	printf("INFO: testing qfloat routines...\n");
	qtest();	// 09/23/2012: Move to after above float-consts-inits because of the qfloat/mi64 routines which use those consts.

	/* Use qfloat routines to set the global floating-point constant 1/sqrt(2): */
	ASSERT(HERE, ISRT2 == qfdbl(QISRT2), "1/sqrt2 precision check failed!");		/* 1/sqrt2	*/

#ifdef CPU_IS_X86
	set_x87_fpu_params(FPU_64RND);
#endif
	// ewm [4. Aug 2014] - move below set_x87_fpu_params(), since need rnd-const set for any DNINT-using ref-computations in the GPU self-tests:
	print_host_info();
	check_nbits_in_types();

	/* Test wide-mul routines: */
	printf("INFO: testing IMUL routines...\n");
	ASSERT(HERE, test_mul() == 0, "test_mul() returns nonzero!");

// Quick timings of various mi64 stuff:
#if 0
	printf("INFO: Testing mi64_add speed...\n");

	/*...time-related stuff	*/
	clock_t clock1, clock2;
	double tdiff;

	int i;
	const int n = 1000, iters = 1000000;
	// Allocate the main data arrays, require these to be on 16-byte boundaries to enable SSE2-based addsub:
	uint64 *u = (uint64 *)calloc(n, sizeof(uint64));	ASSERT(HERE, ((uint32)u & 0xf) == 0, "u not 16-byte aligned!");
	uint64 *v = (uint64 *)calloc(n, sizeof(uint64));	ASSERT(HERE, ((uint32)v & 0xf) == 0, "u not 16-byte aligned!");
	uint64 *x = (uint64 *)calloc(n, sizeof(uint64));	ASSERT(HERE, ((uint32)x & 0xf) == 0, "u not 16-byte aligned!");
	uint64 *y = (uint64 *)calloc(n, sizeof(uint64));	ASSERT(HERE, ((uint32)y & 0xf) == 0, "u not 16-byte aligned!");

	/* Init the RNG and the inputs: */
	rng_isaac_init(TRUE);
	for(i = 0; i < n; i++)
	{
		u[i] = rng_isaac_rand();
		v[i] = rng_isaac_rand();
	}

	// First test correctness:
	uint64 cy1 = mi64_add(u,v,x,n);
	uint64 cy2 = mi64_add_ref(u,v,y,n);
	if(cy1 != cy2) {
		printf("Carryout mismatch: cy1 = %llu, cy2 = %llu\n",cy1,cy2);
	//	ASSERT(HERE, 0, "Incorrect mi64_add carryout");	// GCC 4.4.5 builds on my SB give carry-mismatch here ... wtf?
	}
	for(i = 0; i < n; i++)
	{
		if(x[i] != y[i]) {
			printf("Output mismatch: x[%d] = %llu, y[%d] = %llu\n",i,x[i],i,y[i]);
			ASSERT(HERE, 0, "Incorrect mi64_add output element");
		}
	}

	// Now do timing:
	clock1 = clock();
	for(i = 0; i < iters; i++)
	{
		mi64_add(u,v,x,n);
	}
	clock2 = clock();
	tdiff = (double)(clock2 - clock1);
	printf	("mi64_add: Time for %llu limbs =%s\n",(uint64)iters*n, get_time_str(tdiff));
	exit(0);
#endif
	/************************************************************/
	/* Activate these in turn when a new M-prime is discovered: */
	/************************************************************/
//	print_mp_dec();
//	test_mp_pm1_smooth(p);

#ifdef MULTITHREAD

  #ifndef USE_PTHREAD
	#error PTHREAD define barfed - Did you include e.g. '-pthread' in your compile flags?
  #endif

	/* Test Multithreading: */
	ncpu = get_num_cores();
	nthr = 2*ncpu;
	printf("INFO: System has %d available processor cores.\n", ncpu);

  #if 0	// simple pthreading self-test:
	printf("INFO: Testing Multithreading support with %d threads...\n", nthr);
	// Toggle boolean 2nd arg here to enable verbose mode:
	ASSERT(HERE, test_pthreads(nthr,FALSE) == 0, "test_pthreads() returns nonzero!");
  #endif
#endif

// Define TEST_FFT_RADIX at compile time to activate short-length DFT self-test [Must select params in test_fft_radix.c]
#if defined(TEST_FFT_RADIX)	// && defined(USE_SSE2)
	test_fft_radix();
	exit(0);
#endif
}

/***The following 3 routines MUST BE CALLED IN THE SAME ORDER AS IN host_init()!***/

void print_host_info(void)
{
#ifdef USE_GPU
	gpu_config_t gpu_config;
	gpu_info_t ginfo;
	int32 igpu;

	gpu_init(&gpu_config);
	if (gpu_config.num_gpu > 0) {
		printf("Detected %u CUDA-enabled GPU devices.\n", gpu_config.num_gpu);
		for(igpu = 0; igpu < gpu_config.num_gpu; ++igpu) {
			ginfo = gpu_config.gpu_info[igpu];
			printf("GPU #%u: %s v%u.%u\n", igpu, ginfo.name, ginfo.major, ginfo.minor);
			printf("clock_speed = %u MHz\n", ginfo.clockRate/1000);
			printf("num_compute_units = %u\n", ginfo.multiProcessorCount);
			printf("constant_mem_size = %u\n", ginfo.totalConstMem);
			printf("shared_mem_size = %u\n", ginfo.sharedMemPerBlock);
			printf("global_mem_size = %u\n", ginfo.totalGlobalMem);
			printf("registers_per_block = %u\n", ginfo.regsPerBlock);
			printf("max_threads_per_block = %u\n", ginfo.maxThreadsPerBlock);
			printf("can_overlap = %u\n", ginfo.deviceOverlap);
			printf("warp_size = %u\n", ginfo.warpSize);
			printf("max_thread_dim[3] = [%u,%u,%u]\n", ginfo.maxThreadsDim[0], ginfo.maxThreadsDim[1], ginfo.maxThreadsDim[2]);
			printf("max_grid_size[3] = [%u,%u,%u]\n", ginfo.maxGridSize[0], ginfo.maxGridSize[1], ginfo.maxGridSize[2]);
		}
	} else {
		printf("ERROR: No CUDA-enabled GPUs found\n");
		exit(-1);
	}
//	cudaVecAddTest();
	cudaVecModpowTest0();
	cudaVecModpowTest();
//exit(0);

#endif

#if EWM_DEBUG
	printf("INFO: Program compiled with debugging diagnostics ON.\n");
#endif

	printf("CPU Family = %s, OS = %s, %2d-bit Version, compiled with %s, Version %s.\n", CPU_NAME, OS_NAME, OS_BITS, COMPILER_NAME, COMPILER_VERSION);

#if(defined(CPU_IS_X86) || defined(CPU_IS_X86_64) || defined(CPU_IS_X86_64))

//	get_cpu();

/* Enable this call to get gory details: */
	#if(1)
		/* if(1) --> if(0) Enables section below */
	#elif(defined(COMPILER_TYPE_MSVC) || defined(COMPILER_TYPE_ICC))
		cpu_details();
	#endif

  #if(defined(USE_AVX2))

	if(has_avx2()) {
		printf("Info: Build uses AVX2 instruction set.\n");
	} else {
		ASSERT(HERE, 0, "#define USE_AVX2 invoked but no AVX2 support detected on this CPU! Check get_cpuid functionality and CPU type.\n");
	}

  #elif(defined(USE_AVX))

	if(has_avx2()) {
		printf("Info: CPU supports AVX2 instruction set, but using AVX-enabled build.\n");
	} else if(has_avx()) {
		printf("Info: Build uses AVX instruction set.\n");
	} else {
		ASSERT(HERE, 0, "#define USE_AVX invoked but no AVX support detected on this CPU! Check get_cpuid functionality and CPU type.\n");
	}

  #elif(defined(USE_SSE2))
	/* This doesn't work on non-AVX platforms, since XGETBV (needed by has_avx*() functions) does not exist
	if(has_avx2()) {
		printf("Info: CPU supports AVX2 instruction set, but using SSE2-enabled build.\n");
	} else if(has_avx()) {
		printf("Info: CPU supports AVX instruction set, but using SSE2-enabled build.\n");
	} else */
	if(has_sse2()) {
		printf("INFO: Build uses SSE2 instruction set.\n");
	} else {
		ASSERT(HERE, 0, "#define USE_SSE2 invoked but no SSE2 support detected on this CPU! Check get_cpuid functionality and CPU type.\n");
	}

  #else

	if(has_sse2()) {
		printf("Info: CPU supports SSE2 instruction set, but using scalar floating-point build.\n");
	}

  #endif

#endif

#if PFETCH
	printf("INFO: Using prefetch.\n");
#endif

#ifdef MUL_LOHI64_SUBROUTINE
	printf("INFO: Using subroutine form of MUL_LOHI64.\n");
#else
	printf("INFO: Using inline-macro form of MUL_LOHI64.\n");
#endif

#ifdef USE_FMADD
	printf("INFO: Using FMADD-based 100-bit modmul routines for factoring.\n");
#elif(defined(USE_FLOAT))
	printf("INFO: Using floating-double-based modmul routines for factoring.\n");
#endif

}

/*
80x87 FPU Control Register stuff: FPUCTRL is a 16-bit register. Counting from 0-15,
the 4 bits of interest in the present application are bits 8 and 9 which control the
FPU precision, and bits 10 and 11, which control the FPU rounding mode:

	Bits <9:8>:	Value	Precision		<11:10>:	Value	Rounding Mode
				-----	---------					-----	---------------
				00		24 bits						00		==> nearest (i.e. IEEE)
				01		Reserved					01		==> -oo
				10		53 bits						10		==> +oo
				11		64 bits						11		==> 0 (i.e. truncate)

For the purpose of completeness, the other FPU control bits are as follows
( Adapted from http://maven.smith.edu/~thiebaut/ArtOfAssembly/CH14/CH14-3.html ):

	Bits zero through five are the exception masks. These are similar to the
	interrupt enable bit in the 80x86's flags register. If these bits contain
	a one, the corresponding condition is ignored by the 80x87 FPU. However,
	if any bit contains zero, and the corresponding condition occurs, then
	the FPU immediately generates an interrupt so the program can handle the
	degenerate condition.

	Bit zero corresponds to an invalid operation error. This generally occurs
	as the result of a programming error. Problem which raise the invalid
	operation exception include pushing more than eight items onto the stack
	or attempting to pop an item off an empty stack, taking the square root
	of a negative number, or loading a non-empty register.

	Bit one masks the denormalized interrupt which occurs whenever you try to
	manipulate denormalized values. Denormalized values generally occur when
	you load arbitrary extended precision values into the FPU or work with
	very small numbers just beyond the range of the FPU's capabilities.
	Normally, you would probably not enable this exception.

	Bit two masks the zero divide exception. If this bit contains zero, the FPU
	will generate an interrupt if you attempt to divide a nonzero value by zero.
	If you do not enable the zero division exception, the FPU will produce NaN
	(not a number) whenever you perform a zero division.

	Bit three masks the overflow exception. The FPU will raise the overflow
	exception if a calculation overflows or if you attempt to store a value
	which is too large to fit into a destination operand (e.g., storing a large
	extended precision value into a single precision variable).

	Bit four, if set, masks the underflow exception. Underflow occurs when the
	result is too small to fit in the desintation operand. Like overflow, this
	exception can occur whenever you store a small extended precision value into
	a smaller variable (single or double precision) or when the result of
	computation is too small for extended precision.

	Bit five controls whether the precision exception can occur. A precision
	exception occurs whenever the FPU produces an imprecise result, generally
	the result of an internal rounding operation. Although many operations will
	produce an exact result, many more will not. For example, dividing one by
	ten will produce an inexact result. Therefore, this bit is usually one since
	inexact results are very common.

	Bits six and thirteen through fifteen in the control register are currently
	undefined and reserved for future use. Bit seven is the interrupt enable mask,
	but it is only active on the 8087 FPU; a zero in this bit enables 8087
	interrupts and a one disables FPU interrupts.

	The 80x87 provides two instructions, FLDCW (load control word) and FSTCW (store
	control word), that let you load and store the contents of the control register.
	The single operand to these instructions must be a 16 bit memory location. The
	FLDCW instruction loads the control register from the specified memory location,
	FSTCW stores the control register into the specified memory location.
*/

/* On the Itanium, the FPU Control Register is 64-bit, and divided into 5 subfields:

	- A single 6-bit subfield (bits <5:0> containing the 6 trap disable bits;

	- A quartet of 13-bit Floating-Point Status bitfields (numbered 0 through 3),
	all having the same format
	(A summary the precise meaning of the bits seems exorbitantly difficult to find -
	apparently Intel doesn't want users monkeying with these),
	but with only FPSB 0 accessible to the user, via the following intrinsics
	defined in the <xmmintrin.h> header file:

	extern unsigned __int64 _mm_getfpsr(void);
	extern void _mm_setfpsr(unsigned __int64);	(maps to the Maps to the fsetc.sf0 hardware instruction)

	Using ICC v9.0, the getfpsr function returns 0x8A70037F,
	The least-signficant 12 bits of which are identical to the
	default ones we use for the x87, so whatever the internal order
	of the bitfields, it appears that the _mm_getfpsr and _mm_setfpsr
	intrinsics munge things so the user-settable FPUCTRL bits are in
	the same order and locations as on x87.
*/
#ifdef CPU_IS_X86

	void set_x87_fpu_params(unsigned short FPU_MODE)
	{
		/* SSE FPU control word support: */
	#ifdef USE_SSE2
		int oldMXCSR, newMXCSR;
	#endif
	#ifdef CPU_IS_IA64
		int64 FPUCTRL;
	#else
		unsigned short FPUCTRL;
	#endif

		ASSERT(HERE, (FPU_64RND == FPU_MODE) || (FPU_64CHOP == 0x0f7f), "Illegal value of FPU_MODE");

	#ifdef USE_SSE2
		#ifdef COMPILER_TYPE_MSVC

			__asm	stmxcsr oldMXCSR
			newMXCSR = oldMXCSR | 0x8040; // set DAZ and FZ bits
			__asm ldmxcsr newMXCSR

		#elif(defined(COMPILER_TYPE_GCC) || defined(COMPILER_TYPE_SUNC))

			__asm__ volatile ("stmxcsr %0" : "=m" (oldMXCSR) );
			newMXCSR = oldMXCSR | 0x8040; // set DAZ and FZ bits
			__asm__ volatile ("ldmxcsr %0" :: "m" (newMXCSR) );

		#endif
	#endif

		/* Copy the FPU control word set by the compiler to a local variable
		(mainly so can see what the compiler sets), then overwrite with one of the above:
		*/
	#if(defined(COMPILER_TYPE_ICC) && defined(CPU_IS_X86))

		#error Intel C compiler currently unsupported for x86!

	#elif(FP_MANTISSA_BITS_DOUBLE == 64)/* (defined(CPU_IS_X86) || defined(CPU_IS_IA64))) */

	  #ifdef COMPILER_TYPE_GCC
		#warning Setting rnd_const-emulated DNINT for 64-bit x86 register-double significand
	  #endif

		#ifdef CPU_IS_IA64

			#ifndef COMPILER_TYPE_ICC
				#error unsupported compiler type for ia64!
			#endif
			FPUCTRL = _mm_getfpsr();
			info_x87_fpu_ctrl(FPUCTRL);
			/* Just use the same full-16-bit constant on all platforms, to ensure that there
			are no compiler-based differences in the other 12 bits, either: */
			FPUCTRL = (FPUCTRL & 0x0000) + FPU_MODE;
			_mm_setfpsr(FPUCTRL);

		#else

			#if(defined(COMPILER_TYPE_MWERKS) || defined(COMPILER_TYPE_MSVC))

				__asm	fstcw	FPUCTRL
				/*_controlfp(_PC_64, _MCW_PC);*/

			#elif(defined(COMPILER_TYPE_GCC) || defined(COMPILER_TYPE_SUNC))

				__asm__ volatile ("fstcw %0" : "=m" (FPUCTRL) );

			#else
				#error unsupported compiler type for x87!
			#endif
			info_x87_fpu_ctrl((uint64)FPUCTRL);

			/* Irrespective of what values the compiler set for bitfields <9:8> and <11:10>
			of FPUCTRL, set <9:8> = 11 and <11:10> = 00 to get the full 64-mantissa-bit
			precision available on the x87 and to ensure IEEE rounding mode, respectively:
			*/
		  #if 1
			/* Just use the same full-16-bit constant on all platforms, to ensure that there
			are no compiler-based differences in the other 12 bits, either: */
			FPUCTRL = (FPUCTRL & 0x0000) + FPU_MODE;
		  #else
			***obsolete:***
			FPUCTRL &= 0xf0ff;	/* Clear bits 8:11... */
			FPUCTRL |= 0x0300;	/* And set them to the desired value. */
		  #endif

			/* ...and then reload the FPU control word for the changes to take effect. */
			#if(defined(COMPILER_TYPE_MWERKS) || defined(COMPILER_TYPE_MSVC))
				__asm	fldcw	FPUCTRL
			#elif(defined(COMPILER_TYPE_GCC) || defined(COMPILER_TYPE_SUNC))
				__asm__ volatile ("fldcw %0" :: "m" (FPUCTRL) );
			#endif

		#endif	/* endif(CPU_IS_X86...) */

	#endif	/* endif(FP_MANTISSA_BITS_DOUBLE) */
	}

	void info_x87_fpu_ctrl(uint64 FPUCTRL)
	{
	#if EWM_DEBUG
		printf("INFO: x87 FPU Control Word = %16X.\n", (uint64)FPUCTRL);
	#endif

		/* Check bits <9:8>, and warn if the compiler isn't specifying 64-bit precision: */
		switch ((FPUCTRL >> 8) & 0x3) {
		case 0x3:
			break;
		case 0x2:
			printf("INFO: compiler sets x87 FPU to 53-bit mantissa mode. Overriding...Setting to 64-bit mode.\n");
			break;
		case 0x0:
			printf("INFO: compiler sets x87 FPU to 24-bit mantissa mode. Overriding...Setting to 64-bit mode.\n");
			break;
		default:
			printf("INFO: compiler sets x87 FPU to unknown precision. Overriding...Setting to 64-bit mode.\n");
		}

		/* Check bits <11:10>, and warn if the compiler isn't specifying 64-bit precision: */
		switch ((FPUCTRL >> 10) & 0x3) {
		case 0x0:
			break;
		case 0x1:
			printf("INFO: compiler sets x87 FPU to [round ==> -oo] rounding mode. Overriding...Setting to [round ==> nearest].\n");
			break;
		case 0x2:
			printf("INFO: compiler sets x87 FPU to [round ==> +oo] rounding mode. Overriding...Setting to [round ==> nearest].\n");
			break;
		case 0x3:
			printf("INFO: compiler sets x87 FPU to [round ==> 0] (truncate) rounding mode. Overriding...Setting to [round ==> nearest].\n");
			break;
		default:
			ASSERT(HERE, 0,"0");
		}
	}

#endif	// CPU_IS_X86 ?

/******* DEFINE GLOBALS AND TEST RND-CONST FAST-NINT: *******/
void check_nbits_in_types(void)
{
	uint32 i;
	double ftmp, fran, ferr, finv, fsrt;
	double tpi = 3.1415926535897932384;
	double ln2 = LOG2;

	/* Make sure TRUE and FALSE behave as required: */
	ASSERT(HERE, !FALSE && TRUE, "TRUE and FALSE do not behave as required in check_nbits_in_types");

/* Check lengths of basic data types: */
    ASSERT(HERE, sizeof( int8 ) == 1          , "sizeof( int8 ) != 1          ");
    ASSERT(HERE, sizeof(uint8 ) == 1          , "sizeof(uint8 ) != 1          ");
    ASSERT(HERE, sizeof( int16) == 2          , "sizeof( int16) != 2          ");
    ASSERT(HERE, sizeof(uint16) == 2          , "sizeof(uint16) != 2          ");
    ASSERT(HERE, sizeof( int32) == 4          , "sizeof( int32) != 4          ");
    ASSERT(HERE, sizeof(uint32) == 4          , "sizeof(uint32) != 4          ");
    ASSERT(HERE, sizeof( int64) == 8          , "sizeof( int64) != 8          ");
    ASSERT(HERE, sizeof(uint64) == 8          , "sizeof(uint64) != 8          ");
    ASSERT(HERE, sizeof(uint64) >= sizeof(void*), "sizeof(long long) != sizeof(void*)");    /* ALIGN_DOUBLES assumes this. */

/* AltiVec vector types: */
#if(CPU_HAS_ALTIVEC || CPU_IS_CELL)
	ASSERT(HERE, sizeof(vec_uint8X16) == 16 , "sizeof(vec_uint8X16) != 16 ");
	ASSERT(HERE, sizeof(vec_uint16X8) == 16 , "sizeof(vec_uint16x8) != 16 ");
	ASSERT(HERE, sizeof(vec_uint32X4) == 16 , "sizeof(vec_uint32x4) != 16 ");
#endif

/* We move this into a function defined in a separate file in an
attempt to avoid too-clever compilers realizing that RND_A and RND_B
have the same value and optimizing the +RND_A-RND_B sequences below away:
*/
	get_fp_rnd_const(&RND_A, &RND_B);

	/* Attempted workaround for the ICC v9.0 weirdness here: */
	if(DNINT(tpi) != 3.0 || DNINT(ln2) != 1.0)
	{
		if(FP_MANTISSA_BITS_DOUBLE == 64)
		{
			sprintf(cbuf, "WARN: 64-bit rounding constant not behaving as expected - trying 53-bit version.\n");
			DBG_WARN(HERE, cbuf, "", 0);
			/*
			If there are any platforms which fail here not because they need
			the 53-bit form of the rounding constant but rather because they inevitably
			optimize away the rounding add/sub sequence, but which appear to do the
			right thing in the actual ditN_cy_dif1 routines, indicate in the #if here:
			*/
		#if(1)
		  #define PLATFORM_SKIP_RND_CONST_ENFORCEMENT
		#else
		  #define PLATFORM_SKIP_RND_CONST_ENFORCEMENT
			RND_A = 3.0*0x4000000*0x2000000;
			RND_B =12.0*0x2000000*0x1000000;
		#endif
		}
		else
		{
			sprintf(cbuf, "WARN: 53-bit rounding constant not behaving as expected - trying 64-bit version.\n");
			DBG_WARN(HERE, cbuf, "", 0);
			/*
			If there are any platforms which fail here not because they need
			the 53-bit form of the rounding constant but rather because they inevitably
			optimize away the rounding add/sub sequence, but which appear to do the
			right thing in the actual ditN_cy_dif1 routines, indicate in the #if here:
			*/
		#if(1)
		  #define PLATFORM_SKIP_RND_CONST_ENFORCEMENT
		#else
		  #define PLATFORM_SKIP_RND_CONST_ENFORCEMENT
			RND_A = 3.0*0x4000000*0x2000000*0x800;
			RND_B =12.0*0x2000000*0x1000000*0x800;
		#endif
		}
	}

#ifdef PLATFORM_SKIP_RND_CONST_ENFORCEMENT

	sprintf(cbuf,"in check_nbits_in_types: RND_A = %20.3f, pi  = %20.3f,  DNINT(pi ) = %20.3f\n", RND_A, tpi, (double)DNINT(tpi));
	if((double)DNINT(tpi) != 3.0) {
		DBG_WARN(HERE, cbuf, "", TRUE);
	}
	sprintf(cbuf,"in check_nbits_in_types: RND_A = %20.3f, ln2 = %20.3f,  DNINT(ln2) = %20.3f\n", RND_A, ln2, (double)DNINT(ln2));
	if((double)DNINT(ln2) != 1.0) {
		DBG_WARN(HERE, cbuf, "", TRUE);
	}

#else

	sprintf(cbuf,"in check_nbits_in_types: RND_A = %20.3f, pi  = %20.3f,  DNINT(pi ) = %20.3f\n", RND_A, tpi, (double)DNINT(tpi));
	ASSERT(HERE, (double)DNINT(tpi) == 3.0, cbuf);

	sprintf(cbuf,"in check_nbits_in_types: RND_A = %20.3f, ln2 = %20.3f,  DNINT(ln2) = %20.3f\n", RND_A, ln2, (double)DNINT(ln2));
	ASSERT(HERE, (double)DNINT(ln2) == 1.0, cbuf);

#endif

	/* We typically need more information re. the FFT-mul params before being
	able to inteligently set the anti-thrashing array-padding params, so set = -1
	(which is taken to mean uninitialized) here:
	*/
	DAT_BITS = PAD_BITS = (int32)0xffffffff;

	FFT_MUL_BASE = (double)((uint64)1 << FFT_MUL_BITS);

/* Intend to relax this later to allow powers of 2 as large as 2^54: */
ASSERT(HERE, ((uint64)FFT_MUL_BASE >> 16) == 1, "util.c: FFT_MUL_BASE != 2^16");

	ASSERT(HERE, trailz64((uint64)FFT_MUL_BASE) == FFT_MUL_BITS, "mi64_cvt_double_uint64: trailz64((uint64)FFT_MUL_BASE) != FFT_MUL_BITS");
	ASSERT(HERE, DNINT(FFT_MUL_BASE) == FFT_MUL_BASE, "mi64_cvt_double_uint64: FFT_MUL_BASE not pure-integer!");
	ASSERT(HERE, FFT_MUL_BASE < 1.0*0x8000000*0x8000000, "mi64_cvt_double_uint64: FFT_MUL_BASE >= maximum allowed value of 2^54!");

	FFT_MUL_BASE_INV = 1.0/FFT_MUL_BASE;

	/* Test approximate 1/x and 1/sqrt(x) routines: */
	ftmp = finvest(1.5,  8);	/*fprintf(stderr, "finvest(1.5,  8) gives err = %20.10e\n", fabs(ftmp - 0.666666666666667));*/	ASSERT(HERE, fabs(ftmp - 0.666666666666667) < 4e-03, "Unacceptable level of error in finvest() call!");
	ftmp = finvest(1.5, 53);	/*fprintf(stderr, "finvest(1.5, 53) gives err = %20.10e\n", fabs(ftmp - 0.666666666666667));*/	ASSERT(HERE, fabs(ftmp - 0.666666666666667) < 1e-14, "Unacceptable level of error in finvest() call!");
	ftmp = finvest(1.0, 53);	/*fprintf(stderr, "finvest(1.0, 53) gives err = %20.10e\n", fabs(ftmp - 1.000000000000000));*/	ASSERT(HERE, fabs(ftmp - 1.000000000000000) < 1e-14, "Unacceptable level of error in finvest() call!");
	ftmp = finvest(2.0, 53);	/*fprintf(stderr, "finvest(2.0, 53) gives err = %20.10e\n", fabs(ftmp - 0.500000000000000));*/	ASSERT(HERE, fabs(ftmp - 0.500000000000000) < 1e-14, "Unacceptable level of error in finvest() call!");
	ftmp = finvest(0.5, 53);	/*fprintf(stderr, "finvest(0.5, 53) gives err = %20.10e\n", fabs(ftmp - 2.000000000000000));*/	ASSERT(HERE, fabs(ftmp - 2.000000000000000) < 1e-14, "Unacceptable level of error in finvest() call!");
	ftmp = finvest(.75, 53);	/*fprintf(stderr, "finvest(.75, 53) gives err = %20.10e\n", fabs(ftmp - 1.333333333333333));*/	ASSERT(HERE, fabs(ftmp - 1.333333333333333) < 1e-14, "Unacceptable level of error in finvest() call!");
	/* Try some large and small inputs: */
	ftmp = finvest(3.141592653589793e+15, 53);	/*fprintf(stderr, "finvest(3.141592653589793e+15, 53) gives err = %20.10e\n", fabs(ftmp - 3.183098861837907e-16));*/	ASSERT(HERE, fabs(ftmp - 3.183098861837907e-16) < 1e-14, "Unacceptable level of error in finvest() call!");
	ftmp = finvest(3.183098861837907e-16, 53);	/*fprintf(stderr, "finvest(3.183098861837907e-16, 53) gives err = %20.10e\n", fabs(ftmp - 3.141592653589793e+15));*/	ASSERT(HERE, fabs(ftmp - 3.141592653589793e+15) < 1e+00, "Unacceptable level of error in finvest() call!");

	ftmp = fisqrtest(1.5,  8);	/*fprintf(stderr, "fisqrtest(1.5,  8) gives err = %20.10e\n", fabs(ftmp - 0.816496580927726));*/	ASSERT(HERE, fabs(ftmp - 0.816496580927726) < 1e-3 , "Unacceptable level of error in fisqrtest() call!");
	ftmp = fisqrtest(1.5, 53);	/*fprintf(stderr, "fisqrtest(1.5, 53) gives err = %20.10e\n", fabs(ftmp - 0.816496580927726));*/	ASSERT(HERE, fabs(ftmp - 0.816496580927726) < 1e-14, "Unacceptable level of error in fisqrtest() call!");
	ftmp = fisqrtest(1.0, 53);	/*fprintf(stderr, "fisqrtest(1.0, 53) gives err = %20.10e\n", fabs(ftmp - 1.000000000000000));*/	ASSERT(HERE, fabs(ftmp - 1.000000000000000) < 1e-14, "Unacceptable level of error in fisqrtest() call!");
	ftmp = fisqrtest(2.0, 53);	/*fprintf(stderr, "fisqrtest(2.0, 53) gives err = %20.10e\n", fabs(ftmp - 0.707106781186548));*/	ASSERT(HERE, fabs(ftmp - 0.707106781186548) < 1e-14, "Unacceptable level of error in fisqrtest() call!");
	ftmp = fisqrtest(0.5, 53);	/*fprintf(stderr, "fisqrtest(0.5, 53) gives err = %20.10e\n", fabs(ftmp - 1.414213562373095));*/	ASSERT(HERE, fabs(ftmp - 1.414213562373095) < 1e-14, "Unacceptable level of error in fisqrtest() call!");
	ftmp = fisqrtest(0.3, 53);	/*fprintf(stderr, "fisqrtest(0.3, 53) gives err = %20.10e\n", fabs(ftmp - 1.825741858350554));*/	ASSERT(HERE, fabs(ftmp - 1.825741858350554) < 1e-14, "Unacceptable level of error in fisqrtest() call!");
	ftmp = fisqrtest(.25, 53);	/*fprintf(stderr, "fisqrtest(.25, 53) gives err = %20.10e\n", fabs(ftmp - 2.000000000000000));*/	ASSERT(HERE, fabs(ftmp - 2.000000000000000) < 1e-14, "Unacceptable level of error in fisqrtest() call!");
	ftmp = fisqrtest(.75, 53);	/*fprintf(stderr, "fisqrtest(.75, 53) gives err = %20.10e\n", fabs(ftmp - 1.154700538379251));*/	ASSERT(HERE, fabs(ftmp - 1.154700538379251) < 1e-14, "Unacceptable level of error in fisqrtest() call!");
	ftmp = fisqrtest(3.0, 53);	/*fprintf(stderr, "fisqrtest(3.0, 53) gives err = %20.10e\n", fabs(ftmp - 0.577350269189626));*/	ASSERT(HERE, fabs(ftmp - 0.577350269189626) < 1e-14, "Unacceptable level of error in fisqrtest() call!");
	/* Try some large and small inputs: */
	ftmp = fisqrtest(3.141592653589793e+15, 53);	/*fprintf(stderr, "fisqrtest(3.141592653589793e+15, 53); gives err = %20.10e\n", fabs(ftmp - 1.784124116152771e-08));*/	ASSERT(HERE, fabs(ftmp - 1.784124116152771e-08) < 1e-22, "Unacceptable level of error in fisqrtest() call!");
	ftmp = fisqrtest(3.183098861837907e-16, 53);	/*fprintf(stderr, "fisqrtest(3.183098861837907e-16, 53); gives err = %20.10e\n", fabs(ftmp - 5.604991216397928e+07));*/	ASSERT(HERE, fabs(ftmp - 5.604991216397928e+07) < 1e-07, "Unacceptable level of error in fisqrtest() call!");

	/* Now do a whole mess of 'em: */
	rng_isaac_init(TRUE);
	for(i = 0; i < 100000; i++)
	{
		fran = rng_isaac_rand_double();
		fran = fabs(fran);
		if(fran > 0.0)
		{
			ftmp = finvest  (fran, 53);
			finv = 1.0/fran;
			ferr = (ftmp - finv)/(ftmp + finv);
			ASSERT(HERE, fabs(ferr) < 1e-14, "Unacceptable level of error in finvest  () call!");

			ftmp = fisqrtest(fran, 53);
			fsrt = 1.0/sqrt(fran);
			ferr = (ftmp - fsrt)/(ftmp + fsrt);
			ASSERT(HERE, fabs(ferr) < 1e-14, "Unacceptable level of error in fisqrtest() call!");
		}

		fran = rng_isaac_rand_double_norm_pos();
		if(fran < 0.0 || fran >= 1.0)
		{
			sprintf(cbuf, "check_nbits_in_types: rng_isaac_rand_double_norm_pos returns illegal value outside [0, 1): i = %d, %e\n", i,fran);
			ASSERT(HERE, 0, cbuf);
		}

		fran = rng_isaac_rand_double_norm_pm1();
		if(fabs(fran) >= 1.0)
		{
			sprintf(cbuf, "check_nbits_in_types: rng_isaac_rand_double_norm_pm1 returns illegal value outside (-1,+1): i = %d, %e\n", i, fran);
			ASSERT(HERE, 0, cbuf);
		}
	}

#ifdef USE_FMADD
	/* Test any FMADD-based routines, if def'd: */
	printf("Testing MUL50x50 routines...\n");
	test_mul50x50();
#endif

	return;
}

/********************/

/*
Originally wrote this as part of a workaround for the the DEC Unix V4.0 real*16 sincos bug,
which caused incorrect real*16 sincos results when the argument theta = m*pi/4 +- pi/512, m integer.
The more general purpose is to compare a pair of computed sin(theta), cos(theta) values
(either truncated-to-double extended-real ones, or fast-math-library doubles, or whatever)
to the results returned by calls to the standard math library sincos on the host system.
If a discrepancy between the input sincos data and the standard-library ones is detected
which exceeds some threshold (we use 1e-10 at present), print a warning and replace the
inputs with the latter. (If you instead suspect the latter of being the problematic ones,
you'll need to make the appropriate modifications).

Returns the largest absolute difference between real*8 and real*16 sin(theta) and cos(theta)
for the calling values, e.g. if the user wants to save the maximum discrepancy over a series of calls.
*/
#define	USE_DOUBLE_DEFAULT	0

double 	errprint_sincos(double *x, double *y, double theta)
{
	double tmp, adiff, maxdiff = 0.0;

	tmp = cos(theta);
	adiff = fabs(*x - tmp);
	if(adiff > maxdiff) maxdiff = adiff;
	if(adiff > 1e-10)
	{
		fprintf(stderr, "WARNING: real*16 sine error : theta = %20.15f, long double = %20.15f, double = %20.15f", theta,*x,tmp);
	#if USE_DOUBLE_DEFAULT
		fprintf(stderr, " ... using double-precision values instead.\n");
		*x = tmp;
	#else
		fprintf(stderr, "\n");
	#endif
	}

	tmp = sin(theta);
	adiff = fabs(*y - tmp);
	if(adiff > maxdiff) maxdiff = adiff;
	if(adiff > 1e-10)
	{
		fprintf(stderr, "WARNING: real*16 sine error : theta = %20.15f, long double = %20.15f, double = %20.15f", theta,*y,tmp);
	#if USE_DOUBLE_DEFAULT
		fprintf(stderr, " ... using double-precision values instead.\n");
		*y = tmp;
	#else
		fprintf(stderr, "\n");
	#endif
	}
	return adiff;
}

/**********************************/
/******* INFO, WARN ASSERT ********/
/**********************************/

void INFO(long line, char*file, char*info_string, char*info_file, int copy2stderr)
{
	FILE *fp = 0x0;
	if(STRNEQ(info_file, "")) {
		fp = fopen(info_file,"a");
		if(!fp) fprintf(stderr,"WARNING: unable to open file %s in call to DEBUG_INFO.\n", info_file);
	}

	if(fp) {
		fprintf(fp,"INFO: At line %lu of file %s:\n", line, file);
		fprintf(fp,"%s\n", info_string);
		fclose(fp); fp = 0x0;
	}

	if(copy2stderr || !fp) {
		fprintf(stderr,"INFO: At line %lu of file %s:\n", line, file);
		fprintf(stderr,"%s\n", info_string);
		fflush(stderr);
	}
}

void WARN(long line, char*file, char*warn_string, char*warn_file, int copy2stderr)
{
	FILE *fp = 0x0;
	if(STRNEQ(warn_file, "")) {
		fp = fopen(warn_file,"a");
		if(!fp) fprintf(stderr,"WARNING: unable to open file %s in call to DBG_WARN.\n", warn_file);
	}
	if(fp) {
		fprintf(fp,"WARN: At line %lu of file %s:\n", line, file);
		fprintf(fp,"%s\n", warn_string);
		fclose(fp); fp = 0x0;
	}
	if(copy2stderr || !fp) {
		fprintf(stderr,"WARN: At line %lu of file %s:\n", line, file);
		fprintf(stderr,"%s\n", warn_string);
		fflush(stderr);
	}
}

#ifdef	USE_C99
void ASSERT(char*func, long line, char*file, int expr, char*assert_string)
{
	/* Define a convenient spot to set a breakpoint: */
	if(!expr)
	{
		fprintf(stderr,"ERROR: Function %s, at line %lu of file %s\n", func, line, file);
		fprintf(stderr,"Assertion failed: %s\n", assert_string);
		/* Flush all output streams prior to asserting.
		We replace the original assert(0) call with an exit(EXIT_FAILURE),
		since some compilers seem to like to optimize away assertions. */
		fflush(NULL);
		exit(EXIT_FAILURE);
	}
}
#else
void ASSERT(long line, char*file, int expr, char*assert_string)
{
	/* Define a convenient spot to set a breakpoint: */
	if(!expr)
	{
		fprintf(stderr,"ERROR: at line %lu of file %s\n", line, file);
		fprintf(stderr,"Assertion failed: %s\n", assert_string);
		/* Flush all output streams prior to asserting.
		We replace the original assert(0) call with an exit(EXIT_FAILURE),
		since some compilers seem to like to optimize away assertions. */
		fflush(NULL);
		exit(EXIT_FAILURE);
	}
}
#endif

/***************/

/* ewm: Not sure what I intended this for... */
void	VAR_WARN(char *typelist, ...)
{
	char *c;
	 int32 ival;
	uint32 uval;
	double dval;

	va_list varargs;
	va_start(varargs, typelist);
	/* Define a convenient spot to set a breakpoint: */
	for(c = typelist; *c; c++)
	{
		switch(*c)
		{
			case 'i':
				ival = va_arg(varargs, int32);
				break;
			case 'u':
				uval = va_arg(varargs,uint32);
				break;
			case 'd':
				dval = va_arg(varargs,double);
				break;
			default :
				ASSERT(HERE, 0,"0");
				break;
		}
	}
	va_end(varargs);
}

/****************/

/*...take index i of a set of N = 2^k and return the bit-reversed complement integer.
     Since Mlucas isn't restricted to power-of-2 FFT lengths, we don't actually use
     this function much, preferring the generalized-bit-reversal-on-input-index-vector form,
     but include it for purposes of reference/utility-usage.
*/
/*** REMEMBER: reverse() takes its #bits length arg in exponentiated form n = 2^#bits ***/

int reverse(uint32 i, uint32 n)
{
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

/***************/

/* Slow, default-int versions of leading and trailing-zero-counting algorithms: */
/*
int leadz(int x)
{
	uint32 i;
	for(i = 1; i < 32; i++)
	{
		if(x < 0)break;
		x <<= 1;
	}
	return(i-1);
}
*/

uint32 trailz32(uint32 x)
{
	uint32 i;
	if(x == 0) return 32;
	for(i = 1; i < 32; i++)
	{
		if(x & 1)break;
		x >>= 1;
	}
	return(i-1);
}

uint32 trailz64(uint64 x)
{
#ifdef X64_ASM
	int bpos;
	if(x == 0) return 64;
	__asm__ volatile (\
		"bsfq %[__x],%%rax		\n\t"\
		"movl %%eax,%[__bpos]	\n\t"\
		:	/* outputs: none */\
		: [__x] "m" (x)	/* All inputs from memory addresses here */\
		 ,[__bpos] "m" (bpos)	\
		: "cc","memory","rax"	/* Clobbered registers */\
		);
	return bpos;
#else
	uint32 i;
	const uint64 one64 = 1;
	if(x == 0) return 64;
	for(i = 1; i < 64; i++)
	{
		if(x & one64) {
			break;
		}
		x >>= 1;
	}
	return(i-1);
#endif
}


/***********************************************************************************/
/* Fast 32 and 64-bit-specific versions of leading-zero-counting algorithms: */
/*
An illustrated view of the leading-zeros algorithm:
i1 =
00010111001010101001001110101010, has 3 lzeros.
shift := 16
-1>>16 =
00000000000000001111111111111111 < i1, so i1 has  < 16 lzeros, ly stays = 0, shift -= 8
-1>> 8 =
00000000111111111111111111111111 < i1, so i1 has  <  8 lzeros, ly stays = 0, shift -= 4
-1>> 4 =
00001111111111111111111111111111 < i1, so i1 has  <  4 lzeros, ly stays = 0, shift -= 2
-1>> 2 =
00111111111111111111111111111111 >=i1, so i1 has  >= 2 lzeros, ly += 2  = 2, shift += 1
-1>> 3 =
00011111111111111111111111111111 >=i1, so i1 has  >= 3 lzeros, ly += 1  = 3, shift += 0, done.


i2 =
00000000000000000000101110010101, has 20 lzeros.
shift := 16
-1>>16 =
00000000000000001111111111111111 >=i2, so i2 has  >=16 lzeros, ly +=16  =16, shift += 8
-1>>24 =
00000000000000000000000011111111 < i2, so i2 has  < 24 lzeros, ly stays =16, shift -= 4
-1>>20 =
00000000000000000000111111111111 >=i2, so i2 has  >=20 lzeros, ly +=    =20, shift += 2
-1>>22 =
00000000000000000000001111111111 < i2, so i2 has  < 22 lzeros, ly stays =20, shift -= 1
-1>>21 =
00000000000000000000011111111111 < i2, so i2 has  < 21 lzeros, ly stays =20, shift -= 0, done.
*/

/***************/

uint32 leadz32(uint32 i)
{
	uint32 lz, k, shift, ones_mask = 0xFFFFFFFF;
	if(i == 0) return 32;
	if((int32)i < 0) return 0;

	lz    =  0;
	shift = 16;
	k     = 16;
	while(k > 0)
	{
		if( (ones_mask >> shift) < i )
		{
			k >>= 1;
			shift -= k;
		}
		else
		{
			lz += k;
			k >>= 1;
			shift += k;
		}
	}

	DBG_ASSERT(HERE, ( (i << lz) >> lz == i ),"ERROR A in leadz32");
	DBG_ASSERT(HERE, ( i >> (32-lz)    == 0 ),"ERROR B in leadz32");

	return lz;
}

/***************/

uint32 leadz64(uint64 i)
{
	uint32 lz;
#ifdef X64_ASM
	int bpos;
	if(i == 0) return 64;
	if(( int64)i < 0) return 0;
	__asm__ volatile (\
		"bsrq %[__i],%%rax		\n\t"\
		"movl %%eax,%[__bpos]	\n\t"\
		:	/* outputs: none */\
		: [__i] "m" (i)	/* All inputs from memory addresses here */\
		 ,[__bpos] "m" (bpos)	\
		: "cc","memory","rax"	/* Clobbered registers */\
		);
	lz = (63 - bpos);	// BSR returns *index* of leftmost set bit, must subtract from (#bits - 1) to get #lz.
#else
	uint32 k, shift;
	uint64 ones_mask = 0xFFFFFFFFFFFFFFFFull;
	if(i == 0) return 64;
	if(( int64)i < 0) return 0;
	lz    =  0;
	shift = 32;
	k     = 32;
	while(k > 0)
	{
		if( (ones_mask >> shift) < i )
		{
			k >>= 1;
			shift -= k;
		}
		else
		{
			lz += k;
			k >>= 1;
			shift += k;
		}
	}
	DBG_ASSERT(HERE, ( (i << lz) >> lz == i ),"ERROR A in leadz64");
	DBG_ASSERT(HERE, ( i >> (64-lz)    == 0 ),"ERROR B in leadz64");
#endif
	return lz;
}

uint32	leadz128(uint128 i)
{
	if(i.d1)
	{
		return leadz64(i.d1);
	}
	else
	{
		return leadz64(i.d0) + 64;
	}
}

uint32	leadz192(uint192 i)
{
	if(i.d2)
	{
		return leadz64(i.d2);
	}
	else if(i.d1)
	{
		return leadz64(i.d1) + 64;
	}
	else
	{
		return leadz64(i.d0) + 128;
	}
}

uint32	leadz256(uint256 i)
{
	if(i.d3)
	{
		return leadz64(i.d3);
	}
	else if(i.d2)
	{
		return leadz64(i.d2) + 64;
	}
	else if(i.d1)
	{
		return leadz64(i.d1) + 128;
	}
	else
	{
		return leadz64(i.d0) + 192;
	}
}

/***************/

uint32 isPow2(uint32 i32)
{
/*   If the input == 2^p, returns 1; otherwise returns 0. */
	uint32 i;
	ASSERT(HERE, i32, "isPow2: input zero!");	/* Check for zero in debug mode */
	i = trailz32(i32);
	if((i32 >> i) != 1)
		return 0;
	else
		return 1;
}

/***************/

uint32 isPow4(uint32 i32)
{
/*   If the input == 4^p, returns 1; otherwise returns 0. */
	uint32 i;
	ASSERT(HERE, i32, "isPow4: input zero!");	/* Check for zero in debug mode */
	i = trailz32(i32);
	if((i & 1) || (i32 >> i) != 1)
		return 0;
	else
		return 1;
}

/***************/

/* Population count: Returns the number of set bits in a 32-bit int.
   Optimized for the case where the low-order bits are those most-likely to be set.
*/
uint32 popcount32(uint32 num)
{
	uint32 numones = 0;
	while(num)
	{
		if(num & 1)
			numones++;
		num >>= 1;
	}
	return numones;
}

uint64 popcount64(uint64 num)
{
	uint64 numones = 0;
	while(num)
	{
		if(num & 1)
			numones++;
		num >>= 1;
	}
	return numones;
}

/***************/

/* Extract (nbits) bits beginning at position (beg) */

uint32 ibits32(uint32 i, uint32 beg, uint32 nbits)
{
	uint32 ones_mask = 0xFFFFFFFF;
	return ( (i >> beg) & ~(ones_mask << nbits) );
}

uint64 ibits64(uint64 i, uint32 beg, uint32 nbits)
{
	uint64 ib;
	uint64 ones_mask = 0xFFFFFFFFFFFFFFFFull;
	ib = (i >> beg) & ~(ones_mask << nbits);
	return ( ib );
}

/***************/

/* Return (nbits) bits of a 64-bit integer x, starting at bit
(src_bit_start) in a target 64-bit integer y (the return value), starting at bit (tgt_bit_start).
Entire bit-copy range must lie within bits <0:63> of source operand; any bits which
'overhang' the end of the destination operand are discarded.
If bit-index parameters are illegal, asserts.
*/
uint64	getbits64(uint64 x, uint32 src_bit_start, uint32 nbits, uint32 tgt_bit_start)
{
	const uint64 ones_mask = 0xFFFFFFFFFFFFFFFFull;
	uint64 mask;
	ASSERT(HERE, (nbits <= 64) && (src_bit_start+nbits <= 64) && (tgt_bit_start < 64), "Illegal bit-index parameters!");
	if(nbits == 0) return 0;
	mask = (ones_mask >> (64-nbits));
	return ((x >> src_bit_start) & mask) << tgt_bit_start;
}

/* Alternate version of getbits64, here splicing the requested bit into an argument, leaving the surrounding bits unchanged.
The syntax of this version mirrors that of the Fortran-90 MVBITS library function.
*/
void	mvbits64(uint64 x, uint32 src_bit_start, uint32 nbits, uint64*y, uint32 tgt_bit_start)
{
	const uint64 ones_mask = 0xFFFFFFFFFFFFFFFFull;
	uint64 mask;
	ASSERT(HERE, (nbits <= 64) && (src_bit_start+nbits <= 64) && (tgt_bit_start < 64), "Illegal bit-index parameters!");
	if(nbits == 0) return;
	mask = (ones_mask >> (64-nbits));
	/* Zero out the target bits: */
	*y &= ~(mask << tgt_bit_start);
	/* Copy the source bits into the gap: */
	*y += ((x >> src_bit_start) & mask) << tgt_bit_start;
}

/***************/

int pprimeF(uint32 p, uint32 base)
{
/* returns 1 if p is a base-z Fermat pseudoprime, 0 otherwise. */

	uint64 y = 1, n = p-1, flag;
	uint64 z = base;	/* Need a 64-bit to store intermediate products without overflow */

	while(n)
	{
		flag = n & 1;
		n >>= 1;
		if(flag) y = (y*z)%p;
		z = (z*z)%p;
		if(!z) return 0;
	}
	return((int)(y==1));
}

/***************/

int isPRP(uint32 p)
{
	/* TODO: replace/supplement this with a rigorous trial-divide test for p < 2^32 */
	return(pprimeF(p,2) && pprimeF(p,3) && pprimeF(p,5) && pprimeF(p,7) && pprimeF(p,11) && pprimeF(p,13));
}

/*******************/

/* Calculate 2^-p mod q for p, q 32-bit unsigned ints. This can be used (among
other things) to effect a fast Fermat base-2 pseudoprime test, by calling with q = p-1.
*/
uint32 twompmodq32(uint32 p, uint32 q)	// 2^-p % q
{
	 int32 j;
	uint32 lead5, pshift, qhalf, qinv, zshift, start_index, x, lo, hi;

	ASSERT(HERE, (q&1) == 1, "twopmodq32: even modulus!");
	qhalf = q >> 1;	/* = (q-1)/2, since q odd. */

	pshift = p + 32;
	if(pshift < p)	/* Need special-casing for p just below 2^32  - the primes 2^32-(5,17) are good testcases here. */
	{
		j = -1;	/* leadz32(pshift) for 33-bit pshift goes negative */
		/* Extract leftmost 5 bits of pshift: */
		lead5 = 16 + (pshift >> 28);
	}
	else
	{
		/* Find number of leading zeros in p, use it to find the position of the leftmost ones bit: */
		j = leadz32(pshift);
		/* Extract leftmost 5 bits of pshift: */
		lead5 = ((pshift<<j) >> 27);
	}

	start_index = 32-j-5;	/* Leftward bit at which to start the l-r binary powering, assuming
							the leftmost 5 bits have already been processed via a shift (see next). */

	zshift = 31 - lead5;
	zshift <<= 1;		/* Doubling the shift count here takes cares of the first SQR_LOHI */
	pshift = ~pshift;	/* Overflow doesn't matter here, as long as we got the leading 5 bits of pshift right. */

	qinv = (q+q+q) ^ (uint32)2;	/* Overflow doesn't matter here, since we only care about the low 2 bits of 3*q. */

	qinv = qinv*((uint32)2 - q*qinv);
	qinv = qinv*((uint32)2 - q*qinv);
	qinv = qinv*((uint32)2 - q*qinv);

	/* Since zstart is a power of two < 2^32, use a streamlined code sequence for the first iteration: */
	j = start_index-1;

	/* For 64-bit hardware, Make sure we get a 32-bit shift result here by ANDing with 2^32-1: */
	lo = (qinv << zshift) & (uint32)0xffffffff;
	/* Emulate MULH64 here by getting full 64-bit product and right-shifting: */
	lo = (uint32)(((uint64)q * (uint64)lo) >> 32);
	x  = q - lo;

	if((pshift >> j) & (uint32)1)
	{
		DBG_ASSERT(HERE, x < q,"util.c: x < q");
		/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
		if(x > qhalf) {
			x += x;
			x -= q;
		} else {
			x += x;
		}
	}

	for(j = start_index-2; j >= 0; j--)
	{
		/* SQR_LOHI64(x,lo,hi): */
	#ifdef MUL_LOHI32_SUBROUTINE
		MUL_LOHI32(x,x, lo,hi);
		lo *= qinv;
		lo = __MULH32(q, lo);
	#else
		MUL_LOHI32(x,x, lo,hi);
		lo *= qinv;
		MULH32(q,lo, lo);
	#endif

		/* Branchless version is much faster: */
		x = hi - lo + ((-(hi < lo)) & q);

		if((pshift >> j) & (uint32)1)
		{
			x = x + x - ((-(x > qhalf)) & q);
		}
	}
	/*...Double and return.	These are specialized for the case where 2^p == 1 mod q implies divisibility, in which case x = (q+1)/2. */
	return(x + x - ((-(x > qhalf)) & q));
}

int twopmodq32(uint32 p, uint32 q)	// (2^-p % q) == 0
{
	 int32 j;
	uint32 lead5, pshift, qhalf, qinv, zshift, start_index, x, lo, hi;

	ASSERT(HERE, (q&1) == 1, "twopmodq32: even modulus!");
	qhalf = q >> 1;	/* = (q-1)/2, since q odd. */

	pshift = p + 32;
	if(pshift < p)	/* Need special-casing for p just below 2^32  - the primes 2^32-(5,17) are good testcases here. */
	{
		j = -1;	/* leadz32(pshift) for 33-bit pshift goes negative */
		/* Extract leftmost 5 bits of pshift: */
		lead5 = 16 + (pshift >> 28);
	}
	else
	{
		/* Find number of leading zeros in p, use it to find the position of the leftmost ones bit: */
		j = leadz32(pshift);
		/* Extract leftmost 5 bits of pshift: */
		lead5 = ((pshift<<j) >> 27);
	}

	start_index = 32-j-5;	/* Leftward bit at which to start the l-r binary powering, assuming
							the leftmost 5 bits have already been processed via a shift (see next). */

	zshift = 31 - lead5;
	zshift <<= 1;		/* Doubling the shift count here takes cares of the first SQR_LOHI */
	pshift = ~pshift;	/* Overflow doesn't matter here, as long as we got the leading 5 bits of pshift right. */

	/*
	!    Find modular inverse (mod 2^32) of q in preparation for modular multiply.
	!    We use the simple and elegant iterative inversion method of Montgomery,
	!    which amounts to a modular analogue of Newton's method for iterative inversion:
	!
	!    0)   Zinv = Z                   ! Z*Zinv == 1 (mod 2^3)
	!    1)   Zinv = Zinv*(2 - Z*Zinv)   ! Z*Zinv == 1 (mod 2^6), etc.
	!
	!    where the number of correct bits (at the low end) doubles at each step,
	!    all arithmetic is modulo 2^32 and we repeat step (1) until we have the needed 32 bits.
	!
	!    We choose a different starting value of Zinv, XOR(3*Z, 2),
	!    so the first congruence holds modulo 2^4, thus requiring just 3 iterations.
	*/
	qinv = (q+q+q) ^ (uint32)2;	/* Overflow doesn't matter here, since we only care about the low 2 bits of 3*q. */

	qinv = qinv*((uint32)2 - q*qinv);
	qinv = qinv*((uint32)2 - q*qinv);
	qinv = qinv*((uint32)2 - q*qinv);

	/* Since zstart is a power of two < 2^32, use a streamlined code sequence for the first iteration: */
	j = start_index-1;

	/* For 64-bit hardware, Make sure we get a 32-bit shift result here by ANDing with 2^32-1: */
	lo = (qinv << zshift) & (uint32)0xffffffff;
	/* Emulate MULH64 here by getting full 64-bit product and right-shifting: */
	lo = (uint32)(((uint64)q * (uint64)lo) >> 32);
	x  = q - lo;

	if((pshift >> j) & (uint32)1)
	{
		DBG_ASSERT(HERE, x < q,"util.c: x < q");
		/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
		if(x > qhalf)
		{
			x += x;
			x -= q;
		}
		else
		{
			x += x;
		}
	}

	for(j = start_index-2; j >= 0; j--)
	{
		/* SQR_LOHI64(x,lo,hi): */
	#if 0
		uint64 t = (uint64)x * (uint64)x;	hi = (uint32)(t >> 32);	lo = (uint32)(t) & (uint32)0xffffffff;
		lo *= qinv;
		/* MULH32(q,lo,lo): */
		lo = (uint32)(((uint64)q * (uint64)lo) >> 32);
	#else
	  #ifdef MUL_LOHI32_SUBROUTINE
		MUL_LOHI32(x,x, lo,hi);
		lo *= qinv;
		lo = __MULH32(q, lo);
	  #else
		MUL_LOHI32(x,x, lo,hi);
		lo *= qinv;
		MULH32(q,lo, lo);
	  #endif
	#endif

		/* Branchless version is much faster, but less readable, so give the branched one inside a #if 0: */
	#if 0
		x = hi - lo;
		if(x > hi)
		{
			x += q;	/* had a borrow */
		}
	#else
		x = hi - lo + ((-(hi < lo)) & q);
	#endif

		if((pshift >> j) & (uint32)1)
		{
		/* Branchless version is much faster, but less readable, so give the branched one inside a #if 0: */
		#if 0
			if(x > qhalf)	/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
			{
				x = x + x;
				x -= q;
			}
			else
			{
				x = x + x;
			}
		#else
			x = x + x - ((-(x > qhalf)) & q);
		#endif
		}
	}

	/*...Double and return.	These are specialized for the case where 2^p == 1 mod q implies divisibility, in which case x = (q+1)/2. */
	return((int)((x + x - q) == 1));
}

/* Does an 8-fold base-2 PRP test on the prime candidates q0-7. */
int twopmodq32_x8(uint32 q0, uint32 q1, uint32 q2, uint32 q3, uint32 q4, uint32 q5, uint32 q6, uint32 q7)
{
	int retval = 0;
	 int32 j;
	uint32 start_index;
	uint32 lead0, pshift0, qinv0, zshift0, x0, lo0, hi0, qhalf0;
	uint32 lead1, pshift1, qinv1, zshift1, x1, lo1, hi1, qhalf1;
	uint32 lead2, pshift2, qinv2, zshift2, x2, lo2, hi2, qhalf2;
	uint32 lead3, pshift3, qinv3, zshift3, x3, lo3, hi3, qhalf3;
	uint32 lead4, pshift4, qinv4, zshift4, x4, lo4, hi4, qhalf4;
	uint32 lead5, pshift5, qinv5, zshift5, x5, lo5, hi5, qhalf5;
	uint32 lead6, pshift6, qinv6, zshift6, x6, lo6, hi6, qhalf6;
	uint32 lead7, pshift7, qinv7, zshift7, x7, lo7, hi7, qhalf7;

	DBG_ASSERT(HERE, (q0 < q1) && (q1 < q2) && (q2 < q3) && (q3 < q4) && (q4 < q5) && (q5 < q6) && (q6 < q7), "twopmodq32_x8: Inputs nonmonotone!");

	qhalf0 = q0 >> 1;	/* = (q-1)/2, since q odd. */
	qhalf1 = q1 >> 1;
	qhalf2 = q2 >> 1;
	qhalf3 = q3 >> 1;
	qhalf4 = q4 >> 1;
	qhalf5 = q5 >> 1;
	qhalf6 = q6 >> 1;
	qhalf7 = q7 >> 1;

	/* (p[i]-1)+32 = p + [31,33,37,39,41,49,54,60]: */
	pshift0 = q0 + 31;
	pshift1 = q1 + 31;
	pshift2 = q2 + 31;
	pshift3 = q3 + 31;
	pshift4 = q4 + 31;
	pshift5 = q5 + 31;
	pshift6 = q6 + 31;
	pshift7 = q7 + 31;

	/* Find number of leading zeros in p, use it to find the position of the leftmost ones bit: */
	j = leadz32(pshift0);
	if( leadz32(pshift7) != j )	/* Fused 8-fold algo needs all p's to have same bitlength */
	{
		retval  = (uint32)twopmodq32(q0-1, q0);
		retval += (uint32)twopmodq32(q1-1, q1) << 1;
		retval += (uint32)twopmodq32(q2-1, q2) << 2;
		retval += (uint32)twopmodq32(q3-1, q3) << 3;
		retval += (uint32)twopmodq32(q4-1, q4) << 4;
		retval += (uint32)twopmodq32(q5-1, q5) << 5;
		retval += (uint32)twopmodq32(q6-1, q6) << 6;
		retval += (uint32)twopmodq32(q7-1, q7) << 7;
		return retval;
	}

	if(pshift0 < q0)	/* Need special-casing for p just below 2^32  - the primes 2^32-(5,17) are good testcases here. */
	{
		j = -1;	/* leadz32(pshift) for 33-bit pshift goes negative */
		/* Extract leftmost 5 bits of pshift: */
		lead0 = 16 + (pshift0 >> 28);
		lead1 = 16 + (pshift1 >> 28);
		lead2 = 16 + (pshift2 >> 28);
		lead3 = 16 + (pshift3 >> 28);
		lead4 = 16 + (pshift4 >> 28);
		lead5 = 16 + (pshift5 >> 28);
		lead6 = 16 + (pshift6 >> 28);
		lead7 = 16 + (pshift7 >> 28);
	}
	else
	{
		/* Extract leftmost 5 bits of pshift and subtract from 32: */
		lead0 = ((pshift0<<j) >> 27);
		lead1 = ((pshift1<<j) >> 27);
		lead2 = ((pshift2<<j) >> 27);
		lead3 = ((pshift3<<j) >> 27);
		lead4 = ((pshift4<<j) >> 27);
		lead5 = ((pshift5<<j) >> 27);
		lead6 = ((pshift6<<j) >> 27);
		lead7 = ((pshift7<<j) >> 27);
	}

	start_index = 32-j-5;	/* Leftward bit at which to start the l-r binary powering, assuming
				 the leftmost 5 bits have already been processed via a shift (see next). */

	/* Doubling the shift count here takes cares of the first SQR_LOHI */
	zshift0 = 31 - lead0;	zshift0 <<= 1;	pshift0 = ~pshift0;
	zshift1 = 31 - lead1;	zshift1 <<= 1;	pshift1 = ~pshift1;
	zshift2 = 31 - lead2;	zshift2 <<= 1;	pshift2 = ~pshift2;
	zshift3 = 31 - lead3;	zshift3 <<= 1;	pshift3 = ~pshift3;
	zshift4 = 31 - lead4;	zshift4 <<= 1;	pshift4 = ~pshift4;
	zshift5 = 31 - lead5;	zshift5 <<= 1;	pshift5 = ~pshift5;
	zshift6 = 31 - lead6;	zshift6 <<= 1;	pshift6 = ~pshift6;
	zshift7 = 31 - lead7;	zshift7 <<= 1;	pshift7 = ~pshift7;

	/*
	Find modular inverse (mod 2^32) of q in preparation for modular multiply.
	*/
	qinv0 = (q0+q0+q0) ^ (uint32)2;
	qinv1 = (q1+q1+q1) ^ (uint32)2;
	qinv2 = (q2+q2+q2) ^ (uint32)2;
	qinv3 = (q3+q3+q3) ^ (uint32)2;
	qinv4 = (q4+q4+q4) ^ (uint32)2;
	qinv5 = (q5+q5+q5) ^ (uint32)2;
	qinv6 = (q6+q6+q6) ^ (uint32)2;
	qinv7 = (q7+q7+q7) ^ (uint32)2;
	for(j = 0; j < 3; ++j)
	{
		qinv0 = qinv0*((uint32)2 - q0*qinv0);
		qinv1 = qinv1*((uint32)2 - q1*qinv1);
		qinv2 = qinv2*((uint32)2 - q2*qinv2);
		qinv3 = qinv3*((uint32)2 - q3*qinv3);
		qinv4 = qinv4*((uint32)2 - q4*qinv4);
		qinv5 = qinv5*((uint32)2 - q5*qinv5);
		qinv6 = qinv6*((uint32)2 - q6*qinv6);
		qinv7 = qinv7*((uint32)2 - q7*qinv7);
	}

	/* Since zstart is a power of two < 2^32, use a streamlined code sequence for the first iteration: */
	j = start_index-1;

	/* For 64-bit hardware, Make sure we get a 32-bit shift result here by ANDing with 2^32-1: */
	lo0 = (qinv0 << zshift0) & (uint32)0xffffffff;
	lo1 = (qinv1 << zshift1) & (uint32)0xffffffff;
	lo2 = (qinv2 << zshift2) & (uint32)0xffffffff;
	lo3 = (qinv3 << zshift3) & (uint32)0xffffffff;
	lo4 = (qinv4 << zshift4) & (uint32)0xffffffff;
	lo5 = (qinv5 << zshift5) & (uint32)0xffffffff;
	lo6 = (qinv6 << zshift6) & (uint32)0xffffffff;
	lo7 = (qinv7 << zshift7) & (uint32)0xffffffff;

	/* lo = MULH64(q, lo): */
#ifdef MUL_LOHI32_SUBROUTINE
	lo0 = __MULH32(q0, lo0);
	lo1 = __MULH32(q1, lo1);
	lo2 = __MULH32(q2, lo2);
	lo3 = __MULH32(q3, lo3);
	lo4 = __MULH32(q4, lo4);
	lo5 = __MULH32(q5, lo5);
	lo6 = __MULH32(q6, lo6);
	lo7 = __MULH32(q7, lo7);
#else
	MULH32(q0,lo0, lo0);
	MULH32(q1,lo1, lo1);
	MULH32(q2,lo2, lo2);
	MULH32(q3,lo3, lo3);
	MULH32(q4,lo4, lo4);
	MULH32(q5,lo5, lo5);
	MULH32(q6,lo6, lo6);
	MULH32(q7,lo7, lo7);
#endif

	x0  = q0 - lo0;
	x1  = q1 - lo1;
	x2  = q2 - lo2;
	x3  = q3 - lo3;
	x4  = q4 - lo4;
	x5  = q5 - lo5;
	x6  = q6 - lo6;
	x7  = q7 - lo7;

	/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
	if((pshift0 >> j) & (uint32)1){ DBG_ASSERT(HERE, x0 < q0,"util.c: x0 < q0"); x0 = x0 + x0 - ((-(x0 > qhalf0)) & q0); }
	if((pshift1 >> j) & (uint32)1){ DBG_ASSERT(HERE, x1 < q1,"util.c: x1 < q1"); x1 = x1 + x1 - ((-(x1 > qhalf1)) & q1); }
	if((pshift2 >> j) & (uint32)1){ DBG_ASSERT(HERE, x2 < q2,"util.c: x2 < q2"); x2 = x2 + x2 - ((-(x2 > qhalf2)) & q2); }
	if((pshift3 >> j) & (uint32)1){ DBG_ASSERT(HERE, x3 < q3,"util.c: x3 < q3"); x3 = x3 + x3 - ((-(x3 > qhalf3)) & q3); }
	if((pshift4 >> j) & (uint32)1){ DBG_ASSERT(HERE, x4 < q4,"util.c: x4 < q4"); x4 = x4 + x4 - ((-(x4 > qhalf4)) & q4); }
	if((pshift5 >> j) & (uint32)1){ DBG_ASSERT(HERE, x5 < q5,"util.c: x5 < q5"); x5 = x5 + x5 - ((-(x5 > qhalf5)) & q5); }
	if((pshift6 >> j) & (uint32)1){ DBG_ASSERT(HERE, x6 < q6,"util.c: x6 < q6"); x6 = x6 + x6 - ((-(x6 > qhalf6)) & q6); }
	if((pshift7 >> j) & (uint32)1){ DBG_ASSERT(HERE, x7 < q7,"util.c: x7 < q7"); x7 = x7 + x7 - ((-(x7 > qhalf7)) & q7); }

	for(j = start_index-2; j >= 0; j--)
	{
		/* SQR_LOHI64(x,lo,hi): */
	#ifdef MUL_LOHI32_SUBROUTINE
		MUL_LOHI32(x0,x0, lo0,hi0);
		MUL_LOHI32(x1,x1, lo1,hi1);
		MUL_LOHI32(x2,x2, lo2,hi2);
		MUL_LOHI32(x3,x3, lo3,hi3);
		MUL_LOHI32(x4,x4, lo4,hi4);
		MUL_LOHI32(x5,x5, lo5,hi5);
		MUL_LOHI32(x6,x6, lo6,hi6);
		MUL_LOHI32(x7,x7, lo7,hi7);
		lo0 *= qinv0;
		lo1 *= qinv1;
		lo2 *= qinv2;
		lo3 *= qinv3;
		lo4 *= qinv4;
		lo5 *= qinv5;
		lo6 *= qinv6;
		lo7 *= qinv7;
		lo0 = __MULH32(q0, lo0);
		lo1 = __MULH32(q1, lo1);
		lo2 = __MULH32(q2, lo2);
		lo3 = __MULH32(q3, lo3);
		lo4 = __MULH32(q4, lo4);
		lo5 = __MULH32(q5, lo5);
		lo6 = __MULH32(q6, lo6);
		lo7 = __MULH32(q7, lo7);
	#else
		MUL_LOHI32(x0,x0, lo0,hi0);
		MUL_LOHI32(x1,x1, lo1,hi1);
		MUL_LOHI32(x2,x2, lo2,hi2);
		MUL_LOHI32(x3,x3, lo3,hi3);
		MUL_LOHI32(x4,x4, lo4,hi4);
		MUL_LOHI32(x5,x5, lo5,hi5);
		MUL_LOHI32(x6,x6, lo6,hi6);
		MUL_LOHI32(x7,x7, lo7,hi7);
		lo0 *= qinv0;
		lo1 *= qinv1;
		lo2 *= qinv2;
		lo3 *= qinv3;
		lo4 *= qinv4;
		lo5 *= qinv5;
		lo6 *= qinv6;
		lo7 *= qinv7;
		MULH32(q0,lo0, lo0);
		MULH32(q1,lo1, lo1);
		MULH32(q2,lo2, lo2);
		MULH32(q3,lo3, lo3);
		MULH32(q4,lo4, lo4);
		MULH32(q5,lo5, lo5);
		MULH32(q6,lo6, lo6);
		MULH32(q7,lo7, lo7);
	#endif

		/* if(x < 0) x += q; */
		x0 = hi0 - lo0 + ((-(hi0 < lo0)) & q0);
		x1 = hi1 - lo1 + ((-(hi1 < lo1)) & q1);
		x2 = hi2 - lo2 + ((-(hi2 < lo2)) & q2);
		x3 = hi3 - lo3 + ((-(hi3 < lo3)) & q3);
		x4 = hi4 - lo4 + ((-(hi4 < lo4)) & q4);
		x5 = hi5 - lo5 + ((-(hi5 < lo5)) & q5);
		x6 = hi6 - lo6 + ((-(hi6 < lo6)) & q6);
		x7 = hi7 - lo7 + ((-(hi7 < lo7)) & q7);

		/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
		if((pshift0 >> j) & (uint32)1){ DBG_ASSERT(HERE, x0 < q0,"util.c: x0 < q0"); x0 = x0 + x0 - ((-(x0 > qhalf0)) & q0); }
		if((pshift1 >> j) & (uint32)1){ DBG_ASSERT(HERE, x1 < q1,"util.c: x1 < q1"); x1 = x1 + x1 - ((-(x1 > qhalf1)) & q1); }
		if((pshift2 >> j) & (uint32)1){ DBG_ASSERT(HERE, x2 < q2,"util.c: x2 < q2"); x2 = x2 + x2 - ((-(x2 > qhalf2)) & q2); }
		if((pshift3 >> j) & (uint32)1){ DBG_ASSERT(HERE, x3 < q3,"util.c: x3 < q3"); x3 = x3 + x3 - ((-(x3 > qhalf3)) & q3); }
		if((pshift4 >> j) & (uint32)1){ DBG_ASSERT(HERE, x4 < q4,"util.c: x4 < q4"); x4 = x4 + x4 - ((-(x4 > qhalf4)) & q4); }
		if((pshift5 >> j) & (uint32)1){ DBG_ASSERT(HERE, x5 < q5,"util.c: x5 < q5"); x5 = x5 + x5 - ((-(x5 > qhalf5)) & q5); }
		if((pshift6 >> j) & (uint32)1){ DBG_ASSERT(HERE, x6 < q6,"util.c: x6 < q6"); x6 = x6 + x6 - ((-(x6 > qhalf6)) & q6); }
		if((pshift7 >> j) & (uint32)1){ DBG_ASSERT(HERE, x7 < q7,"util.c: x7 < q7"); x7 = x7 + x7 - ((-(x7 > qhalf7)) & q7); }
	}

	/*...Double and return.	These are specialized for the case where 2^p == 1 mod q implies divisibility, in which case x = (q+1)/2. */
	retval += ((x0 + x0 - q0) == 1)     ;
	retval += ((x1 + x1 - q1) == 1) << 1;
	retval += ((x2 + x2 - q2) == 1) << 2;
	retval += ((x3 + x3 - q3) == 1) << 3;
	retval += ((x4 + x4 - q4) == 1) << 4;
	retval += ((x5 + x5 - q5) == 1) << 5;
	retval += ((x6 + x6 - q6) == 1) << 6;
	retval += ((x7 + x7 - q7) == 1) << 7;
	return retval;
}

/*******************/

/* Simple Euclidean GCD for 32-bit unsigned inputs. (Cf. Algorithm 2.1.4 in Crandall/Pomerance.)
For integers x, y with x, y > 0, returns GCD(x,y). If x or y = 0, returns max(x,y).
*/
uint32 gcd32(uint32 x, uint32 y)
{
	uint32 q, f;

	if(y == 0)
	{
		return x;
	}

	while(y)
	{
		/* Find quotient of current x/y and round toward zero: */
		q = x/y;
		/* Find y' and store in temporary: */
		f = x - q*y;
		/* Find x', i.e. move the old value of y into the slots for x: */
		x = y;
		/* New value of y: */
		y = f;
	}

	return(x);
}

uint64 gcd64(uint64 x, uint64 y)
{
	uint64 q, f;

	if(y == 0)
	{
		return x;
	}

	while(y)
	{
		q = x/y;
		f = x - q*y;
		x = y;
		y = f;
	}

	return(x);
}

/*******************/

/* Simple extended Euclidean GCD for 32-bit unsigned inputs. (Cf. Algorithm 2.1.4 in Crandall/Pomerance.)
For integers x, y with x, y > 0, returns integers {a,b,g} such that a*x + b*y = g = GCD(x,y).

When g = 1 and y > 0, the residues a and b are the inverses of x (mod y) and y (mod x), respectively.

The GCD g is the return value of the function; note that the multipliers a and b
overwrite the inputs x and y, so if the original inputs are needed subsequently,
they must be copied prior to calling the function.
*/
uint32 egcd32(uint32 *x, uint32 *y)
{
	uint32 g = *x, w = *y, q;
	int    a = 1, b = 0, u = 0, v = 1;
	/* Sign of these 3 doesn't matter since they're just temporaries: */
	uint32 d, e, f;

	if(*x == *y)
	{
		fprintf(stderr,"ERROR: eGCD of identical arguments x = y = %u is illegal!\n", *x);
		ASSERT(HERE, 0,"0");
	}
	else if(*x == 0)
	{
		fprintf(stderr,"ERROR: eGCD called with zero input: x = %u, y = %u\n", *x, *y);
		ASSERT(HERE, 0,"0");
	}
	else if(*y == 0)
	{
		fprintf(stderr,"ERROR: eGCD called with zero input: x = %u, y = %u\n", *x, *y);
		ASSERT(HERE, 0,"0");
	}

	while(w)
	{
		/* Find quotient of current x/y and round toward zero: */
		q = g/w;
		/* Find (u', v', w') and store in 3 temporaries: */
		d = a - q*u;
		e = b - q*v;
		f = g - q*w;
		/* Find (a', b', g'), i.e. move the old values of (u,v,w) into the slots for (a,b,g): */
		a = u;
		b = v;
		g = w;
		/* Recover new values of (u, v, w) from the temporaries: */
		u = d;
		v = e;
		w = f;
	}

	*x = a;
	*y = b;
	return(g);
}

/*********************/
/*
Finds multiplicative inverse of z (mod n).
*/
int modinv32(uint32 z, uint32 n)
{
	uint32 x = z, y = n;
	uint32 gcd;

	/* 01/26/04: Turns out this sign check isn't necessary, since the eGCD
				routine automatically handles the case x < y:

	if(x < y)
		gcd = egcd32(&y, &x);
	else if(x > y)
	*/

	gcd = egcd32(&x, &y);
	DBG_ASSERT(HERE, gcd == 1,"gcd in modinv21 is non-unity!");

/*printf("modinv(%u, %u) = %u\n", z, n, x);
	if(n > 100)ASSERT(HERE, 0,"0");
*/
	return x;
}

/********************/

/* Complex multiplication */
struct complex cmul(struct complex *a, struct complex *b)
{
	struct complex cout;
	cout.re = (*a).re*(*b).re - (*a).im*(*b).im;
	cout.im = (*a).re*(*b).im + (*a).im*(*b).re;
	return cout;
}

/***********************************************************************************/
/*
Function to reduce x modulo y, where x and y are both 128-bit unsigned integers.
Algorithm is simple-but-slow bitwise shift-and-subtract scheme.
*/
uint128 xmody128(uint128 x, uint128 y)
{
	uint32 lzx, lzy, nshiftl;
	uint128 t;

	/* In preparation for x%y, Find the # of leading zeros in x and y. */
	     if(x.d1)
		lzx = leadz64(x.d1);
	else
		lzx = leadz64(x.d0) + 64;

	     if(y.d1)
		lzy = leadz64(y.d1);
	else
		lzy = leadz64(y.d0) + 64;

	/* X < Y: return unmodified X. */
	if(lzx > lzy)
		return x;

	nshiftl = lzy - lzx;	/* nshiftlr = 64-nshiftl; */

	while(nshiftl)
	{
		/* Use t to store the left-shifted versions of y: */
		LSHIFT128(y, nshiftl, t);

		if(CMPULT128(t, x))
			SUB128(x, t, x);

		/* Right-shift t one place: */
		--nshiftl;
	}
	/* Must ensure that this gets done once even if lzx == lzy: */
	if(CMPULT128(y, x))
		SUB128(x, y, x);

	return x;
}

/***********************************************************************************/
/*
Function to reduce x modulo y, where x and y are both 192-bit unsigned integers.
Algorithm is simple-but-slow bitwise shift-and-subtract scheme.
Returns remainder x mod y; quotient returned in optional pointer argument q.
*/

uint192 xmody192(const uint192 x, const uint192 y, uint192*quot)
{
	uint32 lzx, lzy, nshiftl;
	uint192 r = x, qsh = ONE192, t;

	/* In preparation for x%y, Find the # of leading zeros in x and y. */
	lzx = leadz192(x);
	lzy = leadz192(y);

	/* X < Y: return unmodified X. */
	if(lzx > lzy)
		return r;

	nshiftl = lzy - lzx;
	if(quot) {
		LSHIFT192(qsh, nshiftl, qsh);	// quotient gets built up from sum of left-shifted binary ones.
		quot->d0 = quot->d1 = quot->d2 = 0ull;
	}
/*
printf("x =%20" LLU "*2^128 + %20" LLU "*2^64 + %20" LLU "\n", x.d2, x.d1, x.d0);
printf("y =%20" LLU "*2^128 + %20" LLU "*2^64 + %20" LLU "\n", y.d2, y.d1, y.d0);
printf("nshiftl = %u\n", nshiftl);
*/
	while(nshiftl)
	{
		/* Use t to store the left-shifted versions of y: */
		LSHIFT192(y, nshiftl, t);
/*printf("y<<%u=" LLU "*2^128 + %20" LLU "*2^64 + %20" LLU "\n", nshiftl, t.d2, t.d1, t.d0); */

		if(CMPULT192(t, r))
		{
			SUB192(r, t, r);
			if(quot) {
				ADD192_PTR(quot, (&qsh), quot);
			}
/*printf("r*=%20" LLU "*2^128 + %20" LLU "*2^64 + %20" LLU "\n", r.d2, r.d1, r.d0); */
		}

		/* Right-shift t one place: */
		--nshiftl;
		if(quot) {
			RSHIFT_FAST192(qsh, 1, qsh);
		}
	}
	/* Must ensure that this gets done once even if lzx == lzy: */
	if(CMPULT192(y, r)) {
		SUB192(r, y, r);
		if(quot) {
			ADD192_PTR(quot, (&qsh), quot);
		}
	}
	return r;
}

/***********************************************************************************/
/*
Function to reduce x modulo y, where x and y are both 256-bit unsigned integers.
Algorithm is simple-but-slow bitwise shift-and-subtract scheme.
Returns remainder x mod y; quotient returned in optional pointer argument q.
*/

uint256 xmody256(const uint256 x, const uint256 y, uint256*quot)
{
	uint32 lzx, lzy, nshiftl;
	uint256 r = x, qsh = ONE256, t;

	/* In preparation for x%y, Find the # of leading zeros in x and y. */
	lzx = leadz256(x);
	lzy = leadz256(y);

	/* X < Y: return unmodified X. */
	if(lzx > lzy)
		return r;

	nshiftl = lzy - lzx;
	if(quot) {
		LSHIFT256(qsh, nshiftl, qsh);	// quotient gets built up from sum of left-shifted binary ones.
		quot->d0 = quot->d1 = quot->d2 = quot->d3 = 0ull;
	}

	while(nshiftl)
	{
		/* Use t to store the left-shifted versions of y: */
		LSHIFT256(y, nshiftl, t);

		if(CMPULT256(t, r))
		{
			SUB256(r, t, r);
			if(quot) {
				ADD256_PTR(quot, (&qsh), quot);
			}
		}

		/* Right-shift t one place: */
		--nshiftl;
		if(quot) {
			RSHIFT_FAST256(qsh, 1, qsh);
		}
	}
	/* Must ensure that this gets done once even if lzx == lzy: */
	if(CMPULT256(y, r)) {
		SUB256(r, y, r);
		if(quot) {
			ADD256_PTR(quot, (&qsh), quot);
		}
	}
	return r;
}


/***********************************************************************************/
/*
Divide-with-Remainder of x by y, where x is a 128-bit unsigned (vector) integer and y a 32-bit unsigned scalar.
Returns (x - x%y)/y in x, 32-bit remainder in the function result.

If you only want the remainder, not to perform the divide, call x128_mod_y32 instead.
*/
uint32 x128_div_y32(uint128 *x, uint32 y)
{
	uint64 cy, rem, xlomody, tsum;
	static uint32 ysave = 0;
	static uint64 two64divy, two64mody;

	if(y != ysave)
	{
		ysave = y;
		two64divy = 0x8000000000000000ull/y;
		two64divy = (two64divy + two64divy);

		two64mody = 0x8000000000000000ull%y;
		/* To save a second (expensive) standard library mod call,
		double and subtract y, then re-add y if the result underflows: */
		two64mody = (two64mody + two64mody) - y;
		cy = (two64mody >> 63);
		two64mody += (-cy) & y;
		two64divy += (cy == 0);
/*printf("INIT: two64divy, two64mody = %20llu %20llu\n\n", two64divy, two64mody); */
	}

	/* Divide high digit by y, storing remainder in cy: */
	cy = (x->d1)%y;
	(x->d1) /= y;

	/* Remainder (must calculate this before modifying (x->d0), obviously): */
	xlomody = (x->d0)%y;
	tsum = cy*two64mody + xlomody;
	rem = tsum%y;

	/* Low digit of result: we must separately divide (x->d0) by y
	(making sure to add (x->d0)%y to  cy*two64mody first, so as not to drop a digit)
	because x->d0 may be as large as 2^64-1, and adding cy*two64mody
	prior to dividing risks unsigned integer overflow:
	*/
	(x->d0) = cy*two64divy + tsum/y + (x->d0)/y;
/*printf("%20llu %20llu %2llu %2llu\n", x->d1, x->d0, cy, rem); */
	return (uint32)rem;
}

uint32 x128_mod_y32(uint128 x, uint32 y)
{
	uint64 cy, rem, xlomody, tsum;
	static uint32 ysave = 0;
	static uint64 two64divy, two64mody;

	if(y != ysave)
	{
		ysave = y;
		two64divy = 0x8000000000000000ull/y;
		two64divy = (two64divy + two64divy);

		two64mody = 0x8000000000000000ull%y;
		/* To save a second (expensive) standard library mod call,
		double and subtract y, then re-add y if the result underflows: */
		two64mody = (two64mody + two64mody) - y;
		cy = (two64mody >> 63);
		two64mody += (-cy) & y;
		two64divy += (cy == 0);
	}

	/* Divide high digit by y, storing remainder in cy: */
	cy = (x.d1)%y;
	/* Remainder: */
	xlomody = (x.d0)%y;
	tsum = cy*two64mody + xlomody;
	rem = tsum%y;

	return (uint32)rem;
}

/***********************************************************************************/
/*
Divide-with-Remainder of x by y, where x is a 192-bit unsigned (vector) integer and y a 32-bit unsigned scalar.
Returns (x - x%y)/y in x, 32-bit remainder in the function result.

If you only want the remainder, not to perform the divide, call x192_mod_y32 instead.
*/
uint32 x192_div_y32(uint192 *x, uint32 y)
{
	uint64 cy, rem, xlomody, tsum;
	static uint32 ysave = 0;
	static uint64 two64divy, two64mody;
	uint128 t128;

	if(y != ysave)
	{
		ysave = y;
		two64divy = 0x8000000000000000ull/y;
		two64divy = (two64divy + two64divy);

		two64mody = 0x8000000000000000ull%y;
		/* To save a second (expensive) standard library mod call,
		double and subtract y, then re-add y if the result underflows: */
		two64mody = (two64mody + two64mody) - y;
		cy = (two64mody >> 63);
		two64mody += (-cy) & y;
		two64divy += (cy == 0);
	}

	/* Copy the upper 2 digits into a local uint128, then call 128-bit divide
	on those, return value is the carry into the low digit: */
	t128.d1 = (x->d2);
	t128.d0 = (x->d1);
	cy = x128_div_y32(&t128, y);
	(x->d2) = t128.d1;
	(x->d1) = t128.d0;

	/* Low digit: */
	xlomody = (x->d0)%y;
	tsum = cy*two64mody + xlomody;
	rem = tsum%y;
	(x->d0) = cy*two64divy + tsum/y + (x->d0)/y;

	return (uint32)rem;
}

uint32 x192_mod_y32(uint192 x, uint32 y)
{
	uint64 cy, rem;
	static uint32 ysave = 0;
	static uint64 two64divy, two64mody;
	uint128 t128;

	if(y != ysave)
	{
		ysave = y;
		two64divy = 0x8000000000000000ull/y;
		two64divy = (two64divy + two64divy);

		two64mody = 0x8000000000000000ull%y;
		/* To save a second (expensive) standard library mod call,
		double and subtract y, then re-add y if the result underflows: */
		two64mody = (two64mody + two64mody) - y;
		cy = (two64mody >> 63);
		two64mody += (-cy) & y;
		two64divy += (cy == 0);
	}

	/* Copy the upper 2 digits into a local uint128, then call 128-bit divide
	on those, return value is the carry into the low digit: */
	t128.d1 = (x.d2);
	t128.d0 = (x.d1);
	cy = x128_div_y32(&t128, y);

	/* Low digit: */
	rem = (cy*two64mody + ((x.d0))%y)%y;

	return (uint32)rem;
}

/***********************************************************************************/
/*
Divide-with-Remainder of x by y, where x is a 256-bit unsigned (vector) integer and y a 32-bit unsigned scalar.
Returns (x - x%y)/y in x, 32-bit remainder in the function result.

If you only want the remainder, not to perform the divide, call x256_mod_y32 instead.
*/
uint32 x256_div_y32(uint256 *x, uint32 y)
{
	uint64 cy, rem, xlomody, tsum;
	static uint32 ysave = 0;
	static uint64 two64divy, two64mody;
	uint192 t192;

	if(y != ysave)
	{
		ysave = y;
		two64divy = 0x8000000000000000ull/y;
		two64divy = (two64divy + two64divy);

		two64mody = 0x8000000000000000ull%y;
		/* To save a second (expensive) standard library mod call,
		double and subtract y, then re-add y if the result underflows: */
		two64mody = (two64mody + two64mody) - y;
		cy = (two64mody >> 63);
		two64mody += (-cy) & y;
		two64divy += (cy == 0);
	}

	/* Copy the upper 3 digits into a local uint192, then call 192-bit divide
	on those, return value is the carry into the low digit: */
	t192.d2 = (x->d3);
	t192.d1 = (x->d2);
	t192.d0 = (x->d1);
	cy = x192_div_y32(&t192, y);
	(x->d3) = t192.d2;
	(x->d2) = t192.d1;
	(x->d1) = t192.d0;

	/* Low digit: */
	xlomody = (x->d0)%y;
	tsum = cy*two64mody + xlomody;
	rem = tsum%y;
	(x->d0) = cy*two64divy + tsum/y + (x->d0)/y;

	return (uint32)rem;
}

uint32 x256_mod_y32(uint256 x, uint32 y)
{
	uint64 cy, rem;
	static uint32 ysave = 0;
	static uint64 two64divy, two64mody;
	uint192 t192;

	if(y != ysave)
	{
		ysave = y;
		two64divy = 0x8000000000000000ull/y;
		two64divy = (two64divy + two64divy);

		two64mody = 0x8000000000000000ull%y;
		/* To save a second (expensive) standard library mod call,
		double and subtract y, then re-add y if the result underflows: */
		two64mody = (two64mody + two64mody) - y;
		cy = (two64mody >> 63);
		two64mody += (-cy) & y;
		two64divy += (cy == 0);
	}

	/* Copy the upper 3 digits into a local uint192, then call 192-bit divide
	on those, return value is the carry into the low digit: */
	t192.d2 = (x.d3);
	t192.d1 = (x.d2);
	t192.d0 = (x.d1);
	cy = x192_div_y32(&t192, y);

	/* Low digit: */
	rem = (cy*two64mody + ((x.d0))%y)%y;

	return (uint32)rem;
}

/***********************************************************************************/

/* Need the uint64 ones of these because some compilers (e.g. MSVC, a.k.a .NET)
don't properly print 64-bit ints. */

/*
Returns decimal character representation of a 64-bit unsigned int in char_buf,
and the position of the leftmost nonzero digit (e.g. if the caller wants to print in
left-justified form) in the function result.
*/
int	convert_uint64_base10_char(char char_buf[], uint64 q)
{
	uint32 i, n_dec_digits = 0, curr_digit;
	char c;
	/* 2^64 has 20 decimal digits - assume the user has allocated at least 20+1 for char_buf: */
	uint32 MAX_DIGITS = 20;

	char_buf[MAX_DIGITS-1]='0';
	char_buf[MAX_DIGITS  ]='\0';

	/* Write the decimal digits into the string from right to left.
	This avoids the need to reverse the digits after calculating them.
	*/
	for(i=0; i < MAX_DIGITS; i++)
	{
		/* Needed to cast modulus 10 to uint32 here for result to come out correct: */
		curr_digit = q%(uint32)10;

		/* Only print leading zero if q = 0, in which case we print a right-justifed single zero: */
		if(q != 0 || n_dec_digits == 0)
		{
			c = curr_digit + CHAROFFSET;
			n_dec_digits++;
		}
		else
			c = ' ';

		char_buf[(MAX_DIGITS - 1) - i] = c;

		q /= 10;
	}

	return (int)MAX_DIGITS-n_dec_digits;
}

/********************/

/*
Returns all-caps hexadecimal character representation of a uint64 in char_buf.
*/
int	convert_uint64_base16_char(char char_buf[], uint64 q)
{
	int i;
	const int hex_chars[16] = {'0','1','2','3','4','5','6','7','8','9','A','B','C','D','E','F'};

	for(i=15; i >= 0; i--)
	{
		char_buf[i] = hex_chars[q & 15];
		q >>= 4;
	}
	char_buf[16] = '\0';
	return 0;
}

/*
For really large inputs we'll want to use base-10^19 for our mod, thus processing nearly one 64-bit
chunk at a time and cutting the number of expensive % operations by 19. But this will also require
us to count leading zeros in the leading (leftmost) base-10^19 word, which isn't worth it for small inputs.
*/

/*
Returns decimal character representation of a base-2^64 2-word unsigned int in char_buf,
and the position of the leftmost nonzero digit (e.g. if the caller wants to print in
left-justified form) in the function result.
*/
int	convert_uint128_base10_char(char char_buf[], uint128 q128)
{
	uint32 i, n_dec_digits = 0;
	char c;
	/* 2^128 has 39 decimal digits: */
	uint32 MAX_DIGITS = 39;

	char_buf[MAX_DIGITS-1]='0';
	char_buf[MAX_DIGITS  ]='\0';

	/* Write the decimal digits into the string from right to left.
	This avoids the need to reverse the digits after calculating them.
	*/
	for(i=0; i < MAX_DIGITS; i++)
	{
		/* Only print leading zero if q = 0, in which case we print a right-justifed single zero: */
		/* Since the x***_div_y32 routines return the mod *and* the divided input,
		   don't call the function until *after* performing the if() test:
		*/
		if((q128.d0 || q128.d1) || n_dec_digits == 0)
		{
			c = x128_div_y32(&q128, (uint32)10) + CHAROFFSET;
			n_dec_digits++;
		}
		else
			c = ' ';

		char_buf[(MAX_DIGITS - 1) - i] = c;
	}

	return (int)MAX_DIGITS-n_dec_digits;
}

int	convert_uint96_base10_char(char char_buf[], uint96 q96)
{
	uint128 q128;
	q128.d0 = q96.d0;
	q128.d1 = (uint64)q96.d1;
	return convert_uint128_base10_char(char_buf, q128);
}

int	convert_uint96ptr_base10_char(char char_buf[], uint96*q96)
{
	uint128 q128;
	q128.d0 = q96->d0;
	q128.d1 = (uint64)q96->d1;
	return convert_uint128_base10_char(char_buf, q128);
}

/*
Returns decimal character representation of a base-2^64 3-word unsigned int in char_buf,
and the position of the leftmost nonzero digit (e.g. if the caller wants to print in
left-justified form) in the function result.
*/
int	convert_uint192_base10_char(char char_buf[], uint192 q192)
{
	uint32 i, n_dec_digits = 0;
	char c;
	/* 2^192 has 58 decimal digits: */
	uint32 MAX_DIGITS = 58;

	char_buf[MAX_DIGITS-1]='0';
	char_buf[MAX_DIGITS  ]='\0';

	/* Write the decimal digits into the string from right to left.
	This avoids the need to reverse the digits after calculating them.
	*/
	for(i=0; i < MAX_DIGITS; i++)
	{
		/* Only print leading zero if q = 0, in which case we print a right-justifed single zero: */
		/* Since the x***_div_y32 routines return the mod *and* the divided input,
		   don't call the function until *after* performing the if() test:
		*/
		if((q192.d0 || q192.d1 || q192.d2) || n_dec_digits == 0)
		{
			c = x192_div_y32(&q192, (uint32)10) + CHAROFFSET;
			n_dec_digits++;
		}
		else
			c = ' ';

		char_buf[(MAX_DIGITS - 1) - i] = c;
	}

	return (int)MAX_DIGITS-n_dec_digits;
}

/*
Returns decimal character representation of a base-2^64 4-word unsigned int in char_buf,
and the position of the leftmost nonzero digit (e.g. if the caller wants to print in
left-justified form) in the function result.
*/
int	convert_uint256_base10_char(char char_buf[], uint256 q256)
{
	uint32 i, n_dec_digits = 0;
	char c;
	/* 2^256 has 78 decimal digits: */
	uint32 MAX_DIGITS = 78;

	char_buf[MAX_DIGITS-1]='0';
	char_buf[MAX_DIGITS  ]='\0';

	/* Write the decimal digits into the string from right to left.
	This avoids the need to reverse the digits after calculating them.
	*/
	for(i=0; i < MAX_DIGITS; i++)
	{
		/* Only print leading zero if q = 0, in which case we print a right-justifed single zero: */
		/* Since the x***_div_y32 routines return the mod *and* the divided input,
		   don't call the function until *after* performing the if() test:
		*/
		if((q256.d0 || q256.d1 || q256.d2 || q256.d3) || n_dec_digits == 0)
		{
			c = x256_div_y32(&q256, (uint32)10) + CHAROFFSET;
			n_dec_digits++;
		}
		else
			c = ' ';

		char_buf[(MAX_DIGITS - 1) - i] = c;
	}

	return (int)MAX_DIGITS-n_dec_digits;
}

/********************/
/* Basically a specialized version of the <stdlib.h> strtod function: */
double	convert_base10_char_double (const char*char_buf)
{
	uint64 curr_sum = (uint64)0;
	double curr_mul = 0.0;
	uint32 i;
	int done_with_leading_whitespace = FALSE;
	char c;
	uint64 curr_digit, hi;

	/* Read the decimal digits from the string from left to right,
	skipping any leading whitespace, and stopping if either non-leading
	whitespace or '\0' is encountered. If we encounter a decimal point,
	the curr_mul multiplier is set = 1.0 and multiplied by 0.1 for every
	numeric digit found to the right of the DP.
	*/
	for(i=0; i != 0xffffffff; i++)
	{
		c = char_buf[i];
		if(!isdigit(c))
		{
			if(isspace(c))
			{
				if(done_with_leading_whitespace)
					break;
				else
					continue;
			}

			done_with_leading_whitespace = TRUE;

			if(c == '.')	/* Found a decimal point */
			{
				ASSERT(HERE, curr_mul == 0.0,"curr_mul == 0.0");	/* Make sure this is the first . we've encountered */
				curr_mul = 1.0;
				continue;
			}
			else if(c == '\n' || c == '\0')
			{
				break;
			}
			else
			{
				fprintf(stderr,"convert_base10_char_double: isdigit(c) fails, s = %s, i = %u, c = %c\n", char_buf, i, c);
				ASSERT(HERE, curr_mul == 0.0,"curr_mul == 0.0");
			}
		}
		curr_mul *= 0.1;	/* Only has an effect if we're to the right of the DP */
		curr_digit = (uint64)(c - CHAROFFSET);
		ASSERT(HERE, curr_digit < 10,"convert_base10_char_double: curr_digit < 10");
		/* Store 10*currsum in a 128-bit product, so can check for overflow: */
	#ifdef MUL_LOHI64_SUBROUTINE
		MUL_LOHI64((uint64)10,curr_sum,&curr_sum,&hi);
	#else
		MUL_LOHI64((uint64)10,curr_sum, curr_sum, hi);
	#endif
		if(hi != 0)
		{
			fprintf(stderr, "ERROR: Mul-by-10 overflows in convert_base10_char_double: Offending input string = %s\n", char_buf);
			ASSERT(HERE, 0,"0");
		}
		curr_sum += curr_digit;	/* Since currsum now a multiple of 10, adding a single digit at the low end can't overflow */
	}

	/* If we encountered no DP we simply convert the pure-integer curr_sum to double
	and return that; otherwise we return (double)curr_sum*curr_mul .
	*/
#if 0
	printf("convert_base10_char_double: char_buf = %s, curr_sum = %llu, curr_mul = %lf\n",char_buf, curr_sum, curr_mul);
#endif
	if(curr_mul == 0.0)
	{
		curr_mul = (double)curr_sum;
	}
	else
	{
		curr_mul *= (double)curr_sum;
	}

	return curr_mul;
}

/********************/
/* Basically a 64-bit version of the <stdlib.h> strtoul function: */
uint64 convert_base10_char_uint64 (const char*char_buf)
{
	uint64 curr_sum = (uint64)0;
	uint32 i;
	int done_with_leading_whitespace = FALSE;
	char c;
	uint64 curr_digit, hi;

	/* Read the decimal digits from the string from left to right,
	skipping any leading whitespace, and stopping if either non-leading
	whitespace or '\0' is encountered:
	*/
	for(i=0; i != 0xffffffff; i++)
	{
		c = char_buf[i];
		if(!isdigit(c))
		{
			if(isspace(c))
			{
				if(done_with_leading_whitespace)
					break;
				else
					continue;
			}

			done_with_leading_whitespace = TRUE;

			if(c == '\n' || c == '\0')
			{
				break;
			}
			else
			{
				fprintf(stderr,"convert_base10_char_uint64: isdigit(c) fails, s = %s, i = %u, c = %c\n", char_buf, i, c);
				ASSERT(HERE, 0,"0");
			}
		}
		curr_digit = (uint64)(c - CHAROFFSET);
		ASSERT(HERE, curr_digit < 10,"convert_base10_char_uint64: curr_digit < 10");
		/* Store 10*currsum in a 128-bit product, so can check for overflow: */
	#ifdef MUL_LOHI64_SUBROUTINE
		MUL_LOHI64((uint64)10,curr_sum,&curr_sum,&hi);
	#else
		MUL_LOHI64((uint64)10,curr_sum, curr_sum, hi);
	#endif
		if(hi != 0)
		{
			fprintf(stderr, "ERROR: Mul-by-10 overflows in convert_base10_char_uint64: Offending input string = %s\n", char_buf);
			ASSERT(HERE, 0,"0");
		}
		curr_sum += curr_digit;	/* Since currsum now a multiple of 10, adding a single digit at the low end can't overflow */
	}

	return curr_sum;
}

uint96	convert_base10_char_uint96 (const char*char_buf)
{
	uint96 rslt;
	uint128 t128 = convert_base10_char_uint128(char_buf);
	rslt.d0 = t128.d0;
	rslt.d1 = (uint32)t128.d1;
	return rslt;
}

uint128	convert_base10_char_uint128(const char*char_buf)
{
	const uint32 LEN_MAX = 2;
	uint64 curr_sum[2] = {(uint64)0,(uint64)0};
	uint64 tmp = 0;
	uint128 x128;
	uint32 i, len = 1;
	int done_with_leading_whitespace = FALSE;
	char c;
	uint64 curr_digit;

	/* Read the decimal digits from the string from left to right,
	skipping any leading whitespace, and stopping if either non-leading
	whitespace or '\0' is encountered:
	*/
	for(i=0; i != 0xffffffff; i++)
	{
		c = char_buf[i];
		if(!isdigit(c))
		{
			if(isspace(c))
			{
				if(done_with_leading_whitespace)
					break;
				else
					continue;
			}

			done_with_leading_whitespace = TRUE;

			if(c == '\n' || c == '\0')
			{
				break;
			}
			else
			{
				fprintf(stderr,"convert_base10_char_uint128: isdigit(c) fails, s = %s, i = %u, c = %c\n", char_buf, i, c);
				ASSERT(HERE, 0,"0");
			}
		}
		curr_digit = (uint64)(c - CHAROFFSET);
		ASSERT(HERE, curr_digit < 10,"util.c: curr_digit < 10");
		/* currsum *= 10, and check for overflow: */
		tmp = mi64_mul_scalar(curr_sum, (uint64)10, curr_sum, len);
		if(tmp != 0)
		{
			if(len == LEN_MAX)
			{
				fprintf(stderr, "ERROR: Mul-by-10 overflows in CONVERT_BASE10_CHAR_UINT128: Offending input string = %s\n", char_buf);
				ASSERT(HERE, len <= LEN_MAX,"len <= LEN_MAX");
			}
			curr_sum[len++] = tmp;
		}

		len += mi64_add_scalar(curr_sum, curr_digit, curr_sum, len);
		ASSERT(HERE, len <= LEN_MAX,"len <= LEN_MAX");
	}

	x128.d0 = curr_sum[0];
	x128.d1 = curr_sum[1];
	return x128;
}

uint192	convert_base10_char_uint192(const char*char_buf)
{
	const uint32 LEN_MAX = 3;
	uint64 curr_sum[3] = {(uint64)0,(uint64)0,(uint64)0};
	uint64 tmp = 0;
	uint192 x192;
	uint32 i, len = 1;
	int done_with_leading_whitespace = FALSE;
	char c;
	uint64 curr_digit;

	/* Read the decimal digits from the string from left to right,
	skipping any leading whitespace, and stopping if either non-leading
	whitespace or '\0' is encountered:
	*/
	for(i=0; i != 0xffffffff; i++)
	{
		c = char_buf[i];
		if(!isdigit(c))
		{
			if(isspace(c))
			{
				if(done_with_leading_whitespace)
					break;
				else
					continue;
			}

			done_with_leading_whitespace = TRUE;

			if(c == '\n' || c == '\0')
			{
				break;
			}
			else
			{
				fprintf(stderr,"convert_base10_char_uint192: isdigit(c) fails, s = %s, i = %u, c = %c\n", char_buf, i, c);
				ASSERT(HERE, 0,"0");
			}
		}
		curr_digit = (uint64)(c - CHAROFFSET);
		ASSERT(HERE, curr_digit < 10,"util.c: curr_digit < 10");
		/* currsum *= 10, and check for overflow: */
		tmp = mi64_mul_scalar(curr_sum, (uint64)10, curr_sum, len);
		if(tmp != 0)
		{
			if(len == LEN_MAX)
			{
				fprintf(stderr, "ERROR: Mul-by-10 overflows in CONVERT_BASE10_CHAR_UINT192: Offending input string = %s\n", char_buf);
				ASSERT(HERE, len <= LEN_MAX,"len <= LEN_MAX");
			}
			curr_sum[len++] = tmp;
		}

		len += mi64_add_scalar(curr_sum, curr_digit, curr_sum, len);
		ASSERT(HERE, len <= LEN_MAX,"len <= LEN_MAX");
	}

	x192.d0 = curr_sum[0];
	x192.d1 = curr_sum[1];
	x192.d2 = curr_sum[2];
	return x192;
}

uint256	convert_base10_char_uint256(const char*char_buf)
{
	const uint32 LEN_MAX = 4;
	uint64 curr_sum[4] = {(uint64)0,(uint64)0,(uint64)0,(uint64)0};
	uint64 tmp = 0;
	uint256 x256;
	uint32 i, len = 1;
	int done_with_leading_whitespace = FALSE;
	char c;
	uint64 curr_digit;

	/* Read the decimal digits from the string from left to right,
	skipping any leading whitespace, and stopping if either non-leading
	whitespace or '\0' is encountered:
	*/
	for(i=0; i != 0xffffffff; i++)
	{
		c = char_buf[i];
		if(!isdigit(c))
		{
			if(isspace(c))
			{
				if(done_with_leading_whitespace)
					break;
				else
					continue;
			}

			done_with_leading_whitespace = TRUE;

			if(c == '\n' || c == '\0')
			{
				break;
			}
			else
			{
				fprintf(stderr,"convert_base10_char_uint256: isdigit(c) fails, s = %s, i = %u, c = %c\n", char_buf, i, c);
				ASSERT(HERE, 0,"0");
			}
		}
		curr_digit = (uint64)(c - CHAROFFSET);
		ASSERT(HERE, curr_digit < 10,"util.c: curr_digit < 10");
		/* currsum *= 10, and check for overflow: */
		tmp = mi64_mul_scalar(curr_sum, (uint64)10, curr_sum, len);
		if(tmp != 0)
		{
			if(len == LEN_MAX)
			{
				fprintf(stderr, "ERROR: Mul-by-10 overflows in CONVERT_BASE10_CHAR_UINT256: Offending input string = %s\n", char_buf);
				ASSERT(HERE, len <= LEN_MAX,"len <= LEN_MAX");
			}
			curr_sum[len++] = tmp;
		}

		len += mi64_add_scalar(curr_sum, curr_digit, curr_sum, len);
		ASSERT(HERE, len <= LEN_MAX,"len <= LEN_MAX");
	}

	x256.d0 = curr_sum[0];
	x256.d1 = curr_sum[1];
	x256.d2 = curr_sum[2];
	x256.d3 = curr_sum[3];
	return x256;
}

/***********************/

/* Functions for 96/128/160/192-bit unsigned integer selected-bit testing: */
uint64	TEST_BIT96 (uint96 __x, uint32 __bit)
{
	/* Since call by value, can overwrite __x here: */
	RSHIFT96(__x, __bit, __x);
	return (__x.d0 & 1);
}

uint64	TEST_BIT128(uint128 __x, uint32 __bit)
{
	/* Since call by value, can overwrite __x here: */
	RSHIFT128(__x, __bit, __x);
	return (__x.d0 & 1);
}

uint64	TEST_BIT160(uint160 __x, uint32 __bit)
{
	/* Since call by value, can overwrite __x here: */
	RSHIFT160(__x, __bit, __x);
	return (__x.d0 & 1);
}

uint64	TEST_BIT192(uint192 __x, uint32 __bit)
{
	/* Since call by value, can overwrite __x here: */
	RSHIFT192(__x, __bit, __x);
	return (__x.d0 & 1);
}

uint64	TEST_BIT256(uint256 __x, uint32 __bit)
{
	/* Since call by value, can overwrite __x here: */
	RSHIFT256(__x, __bit, __x);
	return (__x.d0 & 1);
}

/***********************/

/* Given an IEEE-compliant normalized 64-bit float x, generates an approximate
floating-point inverse accurate to at least (numbits) bits of precision. */
double	finvest(double x, uint32 numbits)
{
	/* Used to store MS 8 non-hidden mantissa bits. We'd need to use a 16-bit int
	to allow for the possibility of a carryout (i.e. result = 256) from rounding
	the 9th-most-significant NHB into the upper 8 (which would involve
	additional logic to handle), we instead deal with the issue of rounding
	by assuming the midpoint - e.g. if truncating to the MS 8 NHBs yields
	a certain integer in [0,255], we assume the resulting roundoff error
	is always 0.5, i.e. our precomputed 1/x values are approximations to
	the resulting midpoints. This also avoids our having to treat an input
	of 1.00000000 as a special case, since we munge that to 1.000000001,
	whose inverse is < 1.0: */
	uint32 byteval;
	int ediff;
	uint32 nacc;
	uint64 itmp, mant, exp;
	double ftmp0, ftmp, err_num, err_den;

	/* Max. precision is 53 bits: */
	if(numbits > 53)
	{
		numbits = 53;
	}

	/* Unpack double into a uint64: */
	itmp = *(uint64 *)&x;
	/* Separate upper part of the significand from the sign/exponent fields: */
	exp  = (itmp >> 52) & MASK_EXP;
	mant =  itmp        & MASK_MANT;
	/* Make sure number is normalized: */
	ASSERT(HERE, exp != 0,"finvest: denormalized inputs illegal!");

	/* Store most-significant 8 non-hidden bits: */
	byteval = (mant >> 44) & 0x000000ff;

	/* Munge the exponent to get the inverse's exponent: double-precision
	1.0 has exponent 1023 and is its own inverse, so that is the corner case:
	Numbers in (1.0, 2.0) have exp = 1023 and inverses with exp = 1022, but
	1.0 *exactly* has inverse with exp = 1023. However, our approximate-midpoint
	scheme obviates the need for extra logic to handle this case - 1.0 gets
	approximated as 1.0 + {small}: */
	ediff = (int)exp - 1023;
	exp = (uint64)(1022 - ediff);

	/* Now get the approx-inverse byte and stick it into the most-significant
	8 non-hidden bits of the mantissa field: */
	mant = (uint64)byte_lookup_finvest[byteval] << 44;

	itmp = (itmp & MASK_SIGN) + (exp << 52) + mant;
	ftmp = *(double *)&itmp;

	/* Do as many Newton iterations as required - number of correct
	bits approximately doubles each iteration. The iteration we use
	for y = 1/x is
					y_{n+1} = y_n*[2 - x*y_n] ,
	which is nice, as it involves no divisions.
	*/
	/* Starting # of correct bits from table lookup = 8: */
	nacc = 8;
ftmp0 = ftmp;
	while(nacc < numbits)
	{
		ftmp = ftmp*(2.0 - x*ftmp);
		nacc += nacc;
	}
	err_num = ftmp - ftmp0;
	err_den = ftmp + ftmp0;
	if(fabs(err_num)/fabs(err_den) >= 2e-3)
	{
		sprintf(cbuf, "finvtest: ftmp0 too inaccurate! ftmp = %e, ftmp0 = %e, relerr = %e\n", ftmp, ftmp0,fabs(err_num)/fabs(err_den));
		ASSERT(HERE, 0, cbuf);
	}

	return ftmp;
}

/* Given an IEEE-compliant normalized 64-bit float x, generates an approximate
floating-point inverse square root accurate to at least (numbits) bits of precision.
This routine is very similar to finvest, so see the comments there for details. */
double	fisqrtest(double x, uint32 numbits)
{
	uint32 byteval;
	int ediff;
	uint32 nacc;
	uint64 itmp, mant, exp;
	double ftmp0, ftmp, err_num, err_den;

	/* Max. precision is 53 bits: */
	if(numbits > 53)
	{
		numbits = 53;
	}

	/* Unpack double into a uint64: */
	itmp = *(uint64 *)&x;
	/* Separate upper part of the significand from the sign/exponent fields: */
	exp  = (itmp >> 52) & MASK_EXP;
	mant =  itmp        & MASK_MANT;
	/* Make sure number is normalized: */
	ASSERT(HERE, exp != 0,"finvest: denormalized inputs illegal!");

	/* Store most-significant 9 non-hidden bits - we'll use either all
	or the high 8 of these, depending on the parity of the exponent: */
	byteval = (mant >> 43) & 0x000001ff;

	/* Munge the exponent to get the inverse square root's exponent: double-precision
	1.0 has exponent 1023 and is its own inverse, so that is the corner case:
	Numbers in (1.0, 4.0) have exp in [1023,1024] and ISQRTs with exp = 1022, but
	1.0 *exactly* has ISQRT with exp = 1023. However, our approximate-midpoint
	scheme obviates the need for extra logic to handle this case - 1.0 gets
	approximated as 1.0 + {small}. However, one additional twist in the 1/sqrt
	case is the asymmetry in the handling of ediff: e.g. 2.0 has ediff = +1 but
	maps to 0.707... with exp = 1022 (i.e. we want 1022 - ediff/2 for inputs > 1),
	but e.g. 0.5 has ediff = -1 but maps to 1/sqrt(0.5) = 1.414... with exp = 1023,
	and 0.3 has ediff = -2 and maps to 1/sqrt(0.3) = 1.825... also with exp = 1023,
	and .25 has ediff = -2 and maps to 1/sqrt(.25) = 2.000..., with exp = 1024,
	i.e. we want 1022 - (ediff-1)/2 for inputs < 1.
	*/
	ediff = (int)exp - 1023;	/* 1023 = 0x3ff */
	if(ediff >= 0)
	{
		exp = (uint64)(1022 - ediff/2);

		/* Since we need to handle mantissas in [1, 4), we differentiate via
		inputs in [1,2) and in [2,4) by examining ediff - if it's even it's the
		former interval and we need do nothing; if odd it's the latter and we
		need to add 2 to the floating version of the mantissa, i.e. 0x100 to byteval: */
		if(ediff & 0x1)
		{
			byteval += 0x100;	/* I realize "byteval" is a misnomer in this case... */
		}
		else
			byteval >>= 1;
	}
	else
	{
		exp = (uint64)(1022 - (ediff-1)/2);

		if(ediff & 0x1)
		{
			byteval += 0x100;	/* I realize "byteval" is a misnomer in this case... */
		}
		else
			byteval >>= 1;
	}

	/* Now get the approx-inverse byte and stick it into the most-significant
	8 non-hidden bits of the mantissa field: */

	mant = (uint64)byte_lookup_fisqrtest[byteval] << 44;

	itmp = (itmp & MASK_SIGN) + (exp << 52) + mant;
	ftmp = *(double *)&itmp;

	/* Do as many Newton iterations as required - number of correct
	bits approximately doubles each iteration. The iteration we use
	for y = 1/sqrt(x) is
					y_{n+1} = y_n*[3 - x*(y_n)^2]/2 ,
	which is nice, as it involves no divisions.
	*/
	/* Starting # of correct bits from table lookup = 8: */
	nacc = 8;
ftmp0 = ftmp;
	while(nacc < numbits)
	{
		ftmp = 0.5*ftmp*(3.0 - x*ftmp*ftmp);
		nacc += nacc;
	}
	err_num = ftmp - ftmp0;
	err_den = ftmp + ftmp0;
	if(fabs(err_num)/fabs(err_den) >= 2e-3)
	{
		sprintf(cbuf, "fisqrtest: ftmp0 too inaccurate! ftmp = %e, ftmp0 = %e, relerr = %e\n", ftmp, ftmp0,fabs(err_num)/fabs(err_den));
		ASSERT(HERE, 0, cbuf);
	}

	return ftmp;
}

/********** Testcode and utils for multithreading support *************/

#ifdef MULTITHREAD

  #if 0	// This works on MacOS, but is non-portable:

	int get_num_cores(void)
	{
		/* get the number of CPUs from the system; 'man sysctl' for details */
		int numCPU;	// Under OS X, this needs to be an int (size_t gave garbage results)
		int mib[4];
		size_t len = sizeof(numCPU);
		/* set the mib for hw.ncpu */
		mib[0] = CTL_HW;
		mib[1] = HW_AVAILCPU;  // alternatively, try HW_NCPU;

		sysctl(mib, 2, &numCPU, &len, NULL, 0);

		if( numCPU < 1 )
		{
			mib[1] = HW_NCPU;
			sysctl( mib, 2, &numCPU, &len, NULL, 0 );

			if( numCPU < 1 )
			{
				numCPU = 1;
			}
		}
		return numCPU;
	}

  #else	// This is alleged to be Win/Linux portable: http://stackoverflow.com/questions/4586405/get-number-of-cpus-in-linux-using-c

	#ifdef OS_TYPE_WINDOWS	// NB: Currently only support || builds unde Linux/GCC, but add Win stuff for possible future use

		#include <windows.h>

		#ifndef _SC_NPROCESSORS_ONLN
			SYSTEM_INFO info;
			GetSystemInfo(&info);
			#define sysconf(a) info.dwNumberOfProcessors
			#define _SC_NPROCESSORS_ONLN
		#endif

	#endif

	int get_num_cores(void)
	{
		long nprocs = -1;
		long nprocs_max = -1;

	#ifdef _SC_NPROCESSORS_ONLN

		nprocs = sysconf(_SC_NPROCESSORS_ONLN);
		if(nprocs < 1) {
			fprintf(stderr, "Could not determine number of CPUs online:\n%s\n", strerror (errno));
			exit (EXIT_FAILURE);
		}
		nprocs_max = sysconf(_SC_NPROCESSORS_CONF);
		if (nprocs_max < 1) {
			fprintf(stderr, "Could not determine number of CPUs configured:\n%s\n", strerror (errno));
			exit (EXIT_FAILURE);
		}
	//	printf ("%ld of %ld processors online\n",nprocs, nprocs_max);
	//	exit (EXIT_SUCCESS);

	#else

		fprintf(stderr, "Could not determine number of CPUs");
		exit (EXIT_FAILURE);

	#endif

		return nprocs;
	}

  #endif

	// Simple struct to pass multiple args to the loop/join-test thread function:
	struct do_loop_test_thread_data{
		int tid;
		int ibeg;
		int iend;
		int *retval;
	};

int test_pthreads(int nthreads, int verbose)
{
	// These are collected from a mish-mash of small code samples I used when initially playing with pthreads;
	// collect all the variable decls at top of this function so this will build under strict ansi C style rules.
	int i,ioffset,tid,j,retval[nthreads];
	pthread_t thread[nthreads];
	pthread_attr_t attr;
	int rc;
	void *status;
	int ibig,iinc,isum;	/* ibig = #bigwords in loop-divided-by-threads sequence */
	struct do_loop_test_thread_data tdat[nthreads];
	pthread_t pth = pthread_self();
    int        thr_id;         /* thread ID for a newly created thread */
    pthread_t  p_thread;       /* thread's structure                     */
    int        a = 1;  /* thread 1 identifying number            */
    int        b = 2;  /* thread 2 identifying number            */
	int ncpu = get_num_cores(), nshift, nextra;
	printf("Mlucas running as system-created pthread %u, threading self-test will use %d user-created pthreads.\n", (int)pth, nthreads);
	if(verbose) {
		ASSERT(HERE, nthreads > 0,"Mlucas.c: nthreads > 0");
		if(nthreads > ncpu) {
			printf("WARN: Test using more threads[%d] than there are available CPUs[%d].\n", nthreads, ncpu);
		}
	}
    /* create a pair of threads, each of which will execute a simple timing loop().
	Uncomment the prints in the thread-called function to 'see' the threads executing: */
    thr_id = pthread_create(&p_thread, NULL, ex_loop, (void*)&a);
	/* Thread which prints a hello message - Note the stdout prints resulting from this
	and the surrounding thread-tests may appear in any order, depending on system scheduling of the respective threads: */
	j = pthread_create(&p_thread, NULL, PrintHello, (void *)&b);
	if (j){
		printf("ERROR; return code from pthread_create() is %d\n", j);
		exit(-1);
	}

	/* Initialize and set thread detached attribute */
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	iinc = 10/nthreads;	/* base loop increment; the first [ibig] chunks get one added to this */
	ibig = 10%nthreads;	/* This many of the [j] work chunks will have an extra unit */
	isum = 0;
	/* Populate the thead-specific data structs: */
	for(i = 0; i < nthreads; ++i) {
		tdat[i].tid = i;
		tdat[i].ibeg = isum;
		isum += iinc + (i < ibig);	/* loop increment for current work chunk */
		tdat[i].iend = isum;
		tdat[i].retval = &retval[i];
		if(verbose) printf("INFO: Scheduling thread %d with ibeg = %d, iend = %d\n", i, tdat[i].ibeg, tdat[i].iend);
	}
	/* create nthreads new threads each of which will execute 'do_loop()' over some specified index subrange.
	In order to match the threads executing at any given time to the available CPUs, divide the thread execution
	into [nshift] 'work shifts', each with [ncpu] threads starting and completing their work before the next shift
	comes online:
	*/
	isum = 0;
	nshift = nthreads / ncpu;	// Number of shifts with one thread for each CPU
	for(j = 0; j < nshift; ++j) {
		ioffset = j*ncpu;
		for(i = 0; i < ncpu; ++i) {
			tid = i+ioffset;
			rc = pthread_create(&thread[tid], &attr, do_loop, (void*)(&tdat[tid]));
			if (rc) {
				printf("ERROR; return code from pthread_create() is %d\n", rc);
				exit(-1);
			}
		}
		/* As each thread finishes, add its result into an accumulator in non-blocking fashion (i.e. no mutexes needed): */
		/* Attempting to join returning threads returns error code ESRCH, 'No such process', if there is just one thread in the current team: */
		if(ncpu > 1) {
			for(i = 0; i < ncpu; ++i) {
				tid = i+ioffset;
				rc = pthread_join(thread[tid], &status);
				if (rc) {
					printf("ERROR; return code from pthread_join() is %d\n", rc);
					exit(-1);
				}
				if(verbose) printf("Main: completed join with thread %d having a status of %d\n",tid,(int)status);
				isum += retval[tid];
			}
		}
	}
	// Cleanup pass for cases where ncpu does not divide nthreads
	nextra = (nthreads % ncpu);
	if(nextra != 0) {
		ioffset = j*ncpu;
		for(i = 0; i < nextra; ++i) {
			tid = i+ioffset;
			rc = pthread_create(&thread[tid], &attr, do_loop, (void*)(&tdat[tid]));
			if (rc) {
				printf("ERROR; return code from pthread_create() is %d\n", rc);
				exit(-1);
			}
		}
		/* As each thread finishes, add its result into an accumulator in non-blocking fashion (i.e. no mutexes needed): */
		if(ncpu > 1) {
			for(i = 0; i < ncpu; ++i) {
				tid = i+ioffset;
				rc = pthread_join(thread[tid], &status);
				if (rc) {
					printf("ERROR; return code from pthread_join() is %d\n", rc);
					exit(-1);
				}
				if(verbose) printf("Main: completed join with thread %d having a status of %d\n",tid,(int)status);
				isum += retval[tid];
			}
		}
	}
	/* Free attribute and wait for the other threads */
	pthread_attr_destroy(&attr);

	// 10 sequential iters of test loop yield successive values -1452071552,1390824192,-61247360,-1513318912,1329576832,
	// -122494720,-1574566272,1268329472,-1837420,-1635813632:
	ASSERT(HERE, isum == -1635813632, "retval error!");
    return 0;
}

// Small timing-delay loop test function for pthread stuff:
void* ex_loop(void* data)
{
	int i;                      /* counter, to print numbers */
	int j;                      /* counter, for delay        */
//	int me = *((int*)data);     /* thread identifying number */
	for (i=0; i<10; i++) {
		for (j=0; j<500000; j++) /* delay loop */
			;
	//	printf("'%d' - Got '%d'\n", me, i);
	}
	/* terminate the thread */
	pthread_exit(NULL);
}

// A little hello-world testcode for the pthread stuff:
void *PrintHello(void *threadid)
{
	int tid;
	tid = *((int*)threadid);
//	printf("Hello World! It's me, thread #%ld!\n", tid);
	pthread_exit(NULL);
}

void*
do_loop(void*targ)	// Thread-arg pointer *must* be cast to void and specialized inside the function
{
	struct do_loop_test_thread_data* thread_arg = targ;
	int i;                      /* counter, to print numbers */
	int j;                      /* counter, for delay        */
	int k = 0;	/* accumulator to keep gcc from otimizing away delay-multiply inside test loop */
	ASSERT(HERE, thread_arg != 0x0, "do_loop test function for pthread-test needs live thread_arg pointer!");

#if 0	// BSD thread affinity API barfs in my Mac builds
	cpuset_t *cset;
	pthread_t pth;
	cpuid_t ci;

	cset = cpuset_create();
	if (cset == NULL) {
		ASSERT(HERE, 0, "cpuset_create");
	}
	ci = 0;
	cpuset_set(ci, cset);

	pth = pthread_self();
	error = pthread_setaffinity_np(pth, cpuset_size(cset), cset);
	if (error) {
		ASSERT(HERE, 0, "pthread_setaffinity_np");
	}
	cpuset_destroy(cset);
#endif

//	int me = thread_arg->tid;     /* thread identifying number */
	for (i = thread_arg->ibeg; i < thread_arg->iend; i++)
	{
		for (j=0; j<100000000; j++) /* delay loop */
		{
			k += j*j;
		}
	//	printf("Thread '%d': i = %d, accum = %d\n", me, i, k);
	}
	*(thread_arg->retval) = k;
	pthread_exit(NULL);
}

#endif

/***********************/

char*get_time_str(double tdiff)
{
	static char cbuf[STR_MAX_LEN];
#ifndef MULTITHREAD
	tdiff /= CLOCKS_PER_SEC;	/* NB: CLOCKS_PER_SEC may be a phony value used to scale clock() ranges */
#endif
	sprintf(cbuf, "%2d%1d:%1d%1d:%1d%1d.%1d%1d%1d"
	,(int)tdiff/36000,((int)tdiff%36000)/3600
	,((int)tdiff%3600)/600,((int)tdiff%600)/60
	,((int)tdiff%60)/10,(int)tdiff%10
	,(int)(10*(tdiff-(int)tdiff)),(int)(100*(tdiff-(int)tdiff))%10,(int)(1000*(tdiff-(int)tdiff))%10);
	return cbuf;
}

