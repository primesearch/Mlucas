/*******************************************************************************
*                                                                              *
*   (C) 1997-2018 by Ernst W. Mayer.                                           *
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

#include "align.h"
#include "util.h"
#include "imul_macro.h"
#ifdef TEST_SIMD
	#include "dft_macro.h"
  #ifdef USE_SSE2
	#include "sse2_macro.h"
	#include "radix16_dif_dit_pass_asm.h"
  #endif
#endif
#ifdef USE_GPU
	#include "gpu_iface.h"
#endif

#if 0
	#define USE_FMADD
	#warning USE_FMADD local-defined!
#endif
/**********************************/
/******* INFO, WARN ASSERT ********/
/**********************************/

void INFO(long line, char*file, char*info_string, char*info_file, int copy2stderr) {
	FILE *fp = 0x0;
	if(STRNEQ(info_file, "")) {
		fp = mlucas_fopen(info_file,"a");
		if(!fp) fprintf(stderr,"WARNING: unable to open file %s in call to DEBUG_INFO.\n", info_file);
	}
	if(fp) {
		fprintf(fp,"INFO: At line %lu of file %s:\n", line, file);	fprintf(fp,"%s\n", info_string);	fclose(fp); fp = 0x0;
	}
	if(copy2stderr || !fp) {
		fprintf(stderr,"INFO: At line %lu of file %s:\n", line, file);	fprintf(stderr,"%s\n", info_string);	fflush(stderr);
	}
}

void WARN(long line, char*file, char*warn_string, char*warn_file, int copy2stderr) {
	FILE *fp = 0x0;
	if(STRNEQ(warn_file, "")) {
		fp = mlucas_fopen(warn_file,"a");
		if(!fp) fprintf(stderr,"WARNING: unable to open file %s in call to DBG_WARN.\n", warn_file);
	}
	if(fp) {
		fprintf(fp,"WARN: At line %lu of file %s:\n", line, file);	fprintf(fp,"%s\n", warn_string);	fclose(fp); fp = 0x0;
	}
	if(copy2stderr || !fp) {
		fprintf(stderr,"WARN: At line %lu of file %s:\n", line, file);	fprintf(stderr,"%s\n", warn_string);	fflush(stderr);
	}
}

#ifdef __CUDA_ARCH__
	/* No-op for GPU device-code compiles: */
	__device__ void ASSERT(long line, char*file, int expr, char*assert_string) {}
#else

  #ifdef USE_C99

	void ASSERT(char*func, long line, char*file, int expr, char*assert_string) {
		/* Define a convenient spot to set a breakpoint: */
		if(!expr) {
			fprintf(stderr,"ERROR: Function %s, at line %lu of file %s\n", func, line, file);	fprintf(stderr,"Assertion failed: %s\n", assert_string);
			/* Flush all output streams prior to asserting. We replace the original assert(0) call with
			an exit(EXIT_FAILURE), since some compilers seem to like to optimize away assertions. */
			fflush(NULL);
			exit(EXIT_FAILURE);
		}
	}

  #else

	void ASSERT(long line, char*file, int expr, char*assert_string) {
		/* Define a convenient spot to set a breakpoint: */
		if(!expr) {
			fprintf(stderr,"ERROR: at line %lu of file %s\n", line, file);	fprintf(stderr,"Assertion failed: %s\n", assert_string);
			/* Flush all output streams prior to asserting. We replace the original assert(0) call with
			an exit(EXIT_FAILURE), since some compilers seem to like to optimize away assertions. */
			fflush(NULL);
			exit(EXIT_FAILURE);	// Try to make this line coincide with a line # == 0 (mod 100) to ease breakpointing
		}
	}

  #endif

#endif	// __CUDA_ARCH__ ?

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

/*** Const-arrays used in POPCNT and related bit-counting functions: ***/
// Bytewise POPCNT:
const uint8 pop8[256] = {
	0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
	1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
	1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
	2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8
};
// Bytewise LEADZ:
const uint8 lz8[256] = {
	8,7,6,6,5,5,5,5,4,4,4,4,4,4,4,4,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
	1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
};
// Bytewise TRAILZ:
const uint8 tz8[256] = {
	8,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,5,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
	6,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,5,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
	7,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,5,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
	6,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,5,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0
};

/* Bytewise ith_set_bit - the 8 4-bit subfields of each uint32 are numbered 1-8 and store the
location [0-7] of the 1st-8th set bit in each byte, with hex-F denoting 'no such set bit'.
Here is the simple C-snippet used to auto-generate the array data:
uint32 i,j,bit,x32;
for(i = 0; i < 256; i++) {
	x32 = bit = 0;
	for(j = 0; (j < 8) && (bit < pop8[i]); j++) {
		if((i>>j)&1) x32 += j<<(4*bit++);
	}
	for(j = bit ; j < 8; j++) {
		x32 += 0xf<<(4*j);
	}
	printf("0x%8X,",x32);
}
printf("\n");
*/
const uint32 ith_set_bit8[256] = {
	0xFFFFFFFF,0xFFFFFFF0,0xFFFFFFF1,0xFFFFFF10,0xFFFFFFF2,0xFFFFFF20,0xFFFFFF21,0xFFFFF210,0xFFFFFFF3,0xFFFFFF30,0xFFFFFF31,0xFFFFF310,0xFFFFFF32,0xFFFFF320,0xFFFFF321,0xFFFF3210,
	0xFFFFFFF4,0xFFFFFF40,0xFFFFFF41,0xFFFFF410,0xFFFFFF42,0xFFFFF420,0xFFFFF421,0xFFFF4210,0xFFFFFF43,0xFFFFF430,0xFFFFF431,0xFFFF4310,0xFFFFF432,0xFFFF4320,0xFFFF4321,0xFFF43210,
	0xFFFFFFF5,0xFFFFFF50,0xFFFFFF51,0xFFFFF510,0xFFFFFF52,0xFFFFF520,0xFFFFF521,0xFFFF5210,0xFFFFFF53,0xFFFFF530,0xFFFFF531,0xFFFF5310,0xFFFFF532,0xFFFF5320,0xFFFF5321,0xFFF53210,
	0xFFFFFF54,0xFFFFF540,0xFFFFF541,0xFFFF5410,0xFFFFF542,0xFFFF5420,0xFFFF5421,0xFFF54210,0xFFFFF543,0xFFFF5430,0xFFFF5431,0xFFF54310,0xFFFF5432,0xFFF54320,0xFFF54321,0xFF543210,
	0xFFFFFFF6,0xFFFFFF60,0xFFFFFF61,0xFFFFF610,0xFFFFFF62,0xFFFFF620,0xFFFFF621,0xFFFF6210,0xFFFFFF63,0xFFFFF630,0xFFFFF631,0xFFFF6310,0xFFFFF632,0xFFFF6320,0xFFFF6321,0xFFF63210,
	0xFFFFFF64,0xFFFFF640,0xFFFFF641,0xFFFF6410,0xFFFFF642,0xFFFF6420,0xFFFF6421,0xFFF64210,0xFFFFF643,0xFFFF6430,0xFFFF6431,0xFFF64310,0xFFFF6432,0xFFF64320,0xFFF64321,0xFF643210,
	0xFFFFFF65,0xFFFFF650,0xFFFFF651,0xFFFF6510,0xFFFFF652,0xFFFF6520,0xFFFF6521,0xFFF65210,0xFFFFF653,0xFFFF6530,0xFFFF6531,0xFFF65310,0xFFFF6532,0xFFF65320,0xFFF65321,0xFF653210,
	0xFFFFF654,0xFFFF6540,0xFFFF6541,0xFFF65410,0xFFFF6542,0xFFF65420,0xFFF65421,0xFF654210,0xFFFF6543,0xFFF65430,0xFFF65431,0xFF654310,0xFFF65432,0xFF654320,0xFF654321,0xF6543210,
	0xFFFFFFF7,0xFFFFFF70,0xFFFFFF71,0xFFFFF710,0xFFFFFF72,0xFFFFF720,0xFFFFF721,0xFFFF7210,0xFFFFFF73,0xFFFFF730,0xFFFFF731,0xFFFF7310,0xFFFFF732,0xFFFF7320,0xFFFF7321,0xFFF73210,
	0xFFFFFF74,0xFFFFF740,0xFFFFF741,0xFFFF7410,0xFFFFF742,0xFFFF7420,0xFFFF7421,0xFFF74210,0xFFFFF743,0xFFFF7430,0xFFFF7431,0xFFF74310,0xFFFF7432,0xFFF74320,0xFFF74321,0xFF743210,
	0xFFFFFF75,0xFFFFF750,0xFFFFF751,0xFFFF7510,0xFFFFF752,0xFFFF7520,0xFFFF7521,0xFFF75210,0xFFFFF753,0xFFFF7530,0xFFFF7531,0xFFF75310,0xFFFF7532,0xFFF75320,0xFFF75321,0xFF753210,
	0xFFFFF754,0xFFFF7540,0xFFFF7541,0xFFF75410,0xFFFF7542,0xFFF75420,0xFFF75421,0xFF754210,0xFFFF7543,0xFFF75430,0xFFF75431,0xFF754310,0xFFF75432,0xFF754320,0xFF754321,0xF7543210,
	0xFFFFFF76,0xFFFFF760,0xFFFFF761,0xFFFF7610,0xFFFFF762,0xFFFF7620,0xFFFF7621,0xFFF76210,0xFFFFF763,0xFFFF7630,0xFFFF7631,0xFFF76310,0xFFFF7632,0xFFF76320,0xFFF76321,0xFF763210,
	0xFFFFF764,0xFFFF7640,0xFFFF7641,0xFFF76410,0xFFFF7642,0xFFF76420,0xFFF76421,0xFF764210,0xFFFF7643,0xFFF76430,0xFFF76431,0xFF764310,0xFFF76432,0xFF764320,0xFF764321,0xF7643210,
	0xFFFFF765,0xFFFF7650,0xFFFF7651,0xFFF76510,0xFFFF7652,0xFFF76520,0xFFF76521,0xFF765210,0xFFFF7653,0xFFF76530,0xFFF76531,0xFF765310,0xFFF76532,0xFF765320,0xFF765321,0xF7653210,
	0xFFFF7654,0xFFF76540,0xFFF76541,0xFF765410,0xFFF76542,0xFF765420,0xFF765421,0xF7654210,0xFFF76543,0xFF765430,0xFF765431,0xF7654310,0xFF765432,0xF7654320,0xF7654321,0x76543210
};

// Def ENABLE_MPRIME_PM1_SMOOTH at compile time to enable Mp p-1 smoothness code:
#ifdef ENABLE_MPRIME_PM1_SMOOTH

	/* Here are the 9366 base-2 Fermat pseudoprimes < 2^32 not divisible by 3 or 5: */
	const uint32 fbase2psp[9366] = {
		/* include contents of f2psp[] in f2psp_3_5.h here */	};
	#error Must include contents of f2psp[] in f2psp_3_5.h here!

	#undef psmooth
	struct psmooth {
		uint32 p;
		uint32 b;	// Standard B-smooth measure based on largest prime factor
		double r;	// L2 "roughness" metric in (0,1] defined by L2 norm of log factor sizes
	};

	// Decimal-print m(p) to a file in 100-digit chunks:
	void print_mp_dec(const uint32 p)
	{
		char *str = 0x0, fname[STR_MAX_LEN];	fname[0] = 'p'; fname[1] = '\0';
		FILE*fp = 0x0;
		uint32 i, lenX, lenD, nchars,nc, wrap_every = 100;	// Insert a LF every 100 digits
		uint64 *x,*y,*d,*r;
		uint64 ONES64 = 0xFFFFFFFFFFFFFFFFull;	// In GCC, making this 'const' gives "warning: overflow in implicit constant conversion" wherever it is used.
		i = convert_uint64_base10_char(fname+1, (uint64)p);
		strcpy(fname+1,fname+i+1);
		strcat(fname,"_decimal.txt");
		// Allocate the array containing first M(p) and then subsequent divide-by-10^100 results.
		// Due to the requirement in mi64_div() that dividend and quotient arrays may not point
		// to the same memory, bounce successive-divide results between 2 arrays, x and y:
		lenX = (p>>6);
	//	x = (uint64 *)calloc(lenX + 1, sizeof(uint64));
		x = (uint64 *)calloc(((lenX + 3) & ~3), sizeof(uint64));	// Zero-pad to make multiple of 4, allowing 64-bit DIV algo to use 4-way-folded loops
		memset(x,ONES64,(lenX<<3));	x[lenX++] = (1ull << (p&63)) - 1;
		nchars = ceil(p * log(2.0)/log(10.));
		fprintf(stderr,"Generating decimal printout of M(%u), which has [%u] decimal digits; will write results to file '%s'...\n",p,nchars,fname);

		// Until have generic-FFT-based mi64_divrem algo in place, use mod-10^27, the largest power of 10 whose
		// odd factor (5^27) fits in a uint64, thus allowing the core div-and-mod loops to use 1-word arguments:
		nc = nchars + (nchars/27) + 1;	// Add newlines to count
		str = (char *)calloc(nc, sizeof(char));
		y = (uint64 *)calloc(lenX + 1, sizeof(uint64));
		// 10^100 has 333 bits, thus needs 6 uint64s, as do the mod-10^100 remainders,
		// but we allow the convert_base10_char_mi64() utility to do the allocation of the former for us:
		ASSERT(HERE, 0x0 != (d = convert_base10_char_mi64("1000000000000000000000000000", &lenD)) && (lenD == 2), "0");
		r = (uint64 *)calloc(lenD, sizeof(uint64));
		nc -= 28;		// starting char of first 27-digit chunk
		for(i = 0; ; i+=2) {	// i = #divides counter; do 2 divs per loop exec in attempt to get some modest pipelining
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
			if((i % 1023) == 0)	// 1M digits ~= 37037 loop execs
				fprintf(stderr,"At digit %u of %5.2fM...\n",27*i,(float)nchars/1000000);
		}
		nc = nchars + (nchars/27) + 1;	// Add newlines to count
		str[nc-1] = '\0';
		fp = mlucas_fopen(fname, "w");
		ASSERT(HERE, fp != 0x0, "Null file pointer!");
		fprintf(fp,"%s\n", str);
		fclose(fp);	fp = 0x0;
		fprintf(stderr,"Done writing %s.",fname);
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

	/* Return the heuristic-estimated number of M-primes in the specified interval.
	From Chris Caldwell http://primes.utm.edu/notes/faq/NextMersenne.html page, citing
	the 1980 Lenstra and Pomerance heuristic analyses: "The probability that 2^p-1 is prime
	is about (e^gamma log ap )/(p log 2) where a=2 if p=3 (mod 4) and a=6 if p=1 (mod 4)."
	Need to sum this over odd primes of the specific residue class in the given interval.
	*/
	double est_num_mp_in_interval(const uint32 plo, const uint32 phi)
	{
		const double iln2 = 1.0/LOG2, eGammaIln2 = 1.78107241799019798523*iln2;	// exp(0.57721566490153286060...)/log2
		double expNumPeq1mod4 = 0.0, expNumPeq3mod4 = 0.0;
		// Small-primes-sieving code ripped off from factor.c:
		const uint32 pdsum_8[8] = { 0, 2, 6, 8,12,18,20,26};
		uint32 curr_p,i,ihi,itmp32,maxp,nprime,neq1mod4 = 0,neq3mod4 = 0;
		uint32 fbase2psp_idx = 0;	// Index to next-expected Fermat base-2 pseudoprime in the precomputed table
		if((phi < 3) || (phi < plo)) return 0.0;
		// Pre-procees p < 11, so can start loop with curr_p = 11 == 1 (mod 10), as required by twopmodq32_x8();
		// Note we wait apply the const-multiplier eGammaIln2 to the final 2 summed estimates:
		nprime = 0;	// #odd primes used
		if((plo < 4) && (phi > 2)) { ++nprime;	++neq3mod4;	curr_p = 3;	expNumPeq3mod4 += log(2.0*curr_p)/curr_p; }
		if((plo < 6) && (phi > 4)) { ++nprime;	++neq1mod4;	curr_p = 5;	expNumPeq1mod4 += log(6.0*curr_p)/curr_p; }
		if((plo < 8) && (phi > 6)) { ++nprime;	++neq3mod4;	curr_p = 7;	expNumPeq3mod4 += log(2.0*curr_p)/curr_p; }
		/* Process chunks of length 30, starting with curr_p == 11 (mod 30). Applying the obvious
		divide-by-3,5 mini-sieve, have 8 candidates in each block: curr_p + [ 0, 2, 6, 8,12,18,20,26].
		For example: curr_p = 11 gives the 8 candidates: 11,13,17,19,23,29,31,37.
		*/
		maxp = MIN(phi,0xffffffe3);	// Make sure (curr_p + 29) < 2^32 in our loop
		for(curr_p = 11; curr_p <= maxp; curr_p += 30) {
			/* Do a quick Fermat base-2 compositeness test before invoking the more expensive mod operations: */
			itmp32 = twopmodq32_x8(curr_p, curr_p+ 2, curr_p+ 6, curr_p+ 8, curr_p+12, curr_p+18, curr_p+20, curr_p+26);
			for(i = 0; i < 8; ++i) {
				// It's a PRP: check vs table of known pseudoprimes and (if it's not a PSP) init for the next PSP:
				if((itmp32 >> i)&0x1) {
					ASSERT(HERE, curr_p <= fbase2psp[fbase2psp_idx],"Error in pseudoprime sieve");
					if((curr_p + pdsum_8[i]) == fbase2psp[fbase2psp_idx]) {	// It's a base-2 pseudoprime
						++fbase2psp_idx;
						continue;
					} else {	// It's prime:
						ihi = (curr_p + pdsum_8[i]);
						if(ihi < plo) continue;
						if(ihi > maxp) break;
						++nprime;
					//	printf("At prime = %u, (mod 4) = %u\n",ihi,ihi&3);
						if((ihi&3) == 1) { ++neq1mod4;	expNumPeq1mod4 += log(6.0*ihi)/ihi; }
						if((ihi&3) == 3) { ++neq3mod4;	expNumPeq3mod4 += log(2.0*ihi)/ihi; }
					}
				}
			}
		}
		expNumPeq1mod4 *= eGammaIln2;	expNumPeq3mod4 *= eGammaIln2;
		printf("Using %u odd primes in [%u,%u], of which (%u,%u) == 1,3 (mod 4); Expected #Mp with p == 1,3 (mod 4) = %8.3f, %8.3f\n",nprime,plo,phi,neq1mod4,neq3mod4,expNumPeq1mod4,expNumPeq3mod4);
		printf("Max prime used = %u\n",ihi);
		return expNumPeq1mod4 + expNumPeq3mod4;
	}

	/* Linear least-squares applied to lg(p) for known-M(p) exponents, as described at http://primes.utm.edu/notes/faq/NextMersenne.html
	Jan 2016: See http://www.mersenneforum.org/showthread.php?p=423266#post423266 for results based on latest, M#49
	*/
	void compute_mers_best_fit()
	{
		const double iln2 = 1.0/LOG2;
		double xi, xavg, yavg, num, den, a,b,
			y[] = {2,3,5,7,13,17,19,31,61,89,107,127,521,607,1279,2203,2281,3217,4253,4423,9689,9941
			,11213,19937,21701,23209,44497,86243,110503,132049,216091,756839,859433,1257787,1398269,2976221,3021377,6972593
			,13466917,20996011,24036583,25964951,30402457,32582657,37156667,42643801,43112609,57885161,74207281,0.0};
		int i,j,n,p,eq1mod4 = 0,starts_with_peq2 = (y[0]==2);
		// Convert exponents p (stored as doubles in y-array) into lg(p), and compute averages:
		for(i = 0, yavg = 0.0; y[i] != 0; i++) {
			eq1mod4 += ((uint32)y[i]&3) == 1;
		//	printf("p = %8u, p%%4 = %u\n",(uint32)y[i],(uint32)y[i]&3);
			y[i] = log(y[i])*iln2;	yavg += y[i];
		}	n = i;
		xavg = (1.0 + n)/2;	yavg /= n;	// X-avg uses unit offset
		printf("#M-prime exponents = %u, #==1,3 (mod4) = %u,%u\n",n,eq1mod4,(n-eq1mod4-starts_with_peq2));
		printf("Sample size = %u, xavg = %8.4f, yavg = %8.4f\n",n,xavg,yavg);
		/*
		Linear least squares: Assume best-fit is to line a*x+b, slope a and y-intercept b TBD.
		Each datapoint has x = [index of Mersenne number], hence exact, under assumption of no
		as-yet-undiscivered primes with p less than max_p of of our knowns[] dataset. Ith point
		has 'error' w.r.to best-fit line measured via y-offset (as opposed to, say, normal distance,
		i.e. a total-least-squares approach as to the 'ordinary' one here, which would be more
		appropriate for data with non-exact x-values), di := yi - (a*xi+b). We
		seek a,b such that the sums of the squares of the di-values for our dataset is minimized.

			S = sum_i [yi - (a*xi+b)]^2 = sum_i [yi^2 - 2*yi*(a*xi+b) + (a^2*xi^2 + 2*a*b*xi + b^2)] .

		Take partial derivative of S w.r.to a: dS/da = sum_i [-2*xi*yi + 2*a*xi^2 + 2*b*xi] = 0. [1]

		Take partial derivative of S w.r.to b: dS/db = sum_i [-2*yi + 2*a*xi + 2*b] = 0. [2]

		Since b appears in [2] unmultiplied by xi or yi, can easily solve for it: b = y@ - a*x@, [**]
		where @ denotes the sample mean of the qty in question: x@ = [sum_i xi]/n, y@ = [sum_i yi]/n .
		Substituting the expression for b [**] into [1] we can solve for the slope paramater:

			sum_i [-2*xi*yi + 2*a*xi^2 + 2*b*xi] = 0, div by 2 and sub for b:
		->	sum_i [-xi*yi + a*xi^2 + (y@ - a*x@)*xi] = 0
		->	sum_i [-xi*yi + a*xi^2 + y@*xi - a*x@*xi] = 0
		->	sum_i [-(yi - y@)*xi + a*(xi - x@)*xi] = 0, separate into 2 sums, pull a out of 2nd one and solve for it:

		a = sum_i [(yi - y@)*xi] / sum_i [(xi - x@)*xi]. [*]
		*/
		for(i = 0, num = den = 0.0; i < n; i++) {
			xi = i+1;	num += (y[i] - yavg)*xi;	den += (xi - xavg)*xi;
		}
		a = num/den;	b = yavg - a*xavg;
		printf("Least-squares of full %u-point dataset gives slope = %8.4f, y-intercept = %8.4f\n",n,a,b);

	// Now do another linear regression, this time omitting smallest 10 M-exponents:
		for(i = 10, yavg = 0.0; i < n; i++) {
			yavg += y[i];
		}
		xavg = (1.0 + 10 + n)/2;	yavg /= (n-10);	// X-avg uses unit offset
		printf("Omitting 10 smallest M(p): Sample size = %u, xavg = %8.4f, yavg = %8.4f\n",n-10,xavg,yavg);
		for(i = 10, num = den = 0.0; i < n; i++) {
			xi = i+1;	num += (y[i] - yavg)*xi;	den += (xi - xavg)*xi;
		}
		a = num/den;	b = yavg - a*xavg;
		printf("Least-squares omitting 10 smallest M(p) gives slope = %8.4f, y-intercept = %8.4f\n",a,b);

	// Now do another linear regression, this time omitting smallest 20 M-exponents:
		for(i = 20, yavg = 0.0; i < n; i++) {
			yavg += y[i];
		}
		xavg = (1.0 + 20 + n)/2;	yavg /= (n-20);	// X-avg uses unit offset
		printf("Omitting 20 smallest M(p): Sample size = %u, xavg = %8.4f, yavg = %8.4f\n",n-20,xavg,yavg);
		for(i = 20, num = den = 0.0; i < n; i++) {
			xi = i+1;	num += (y[i] - yavg)*xi;	den += (xi - xavg)*xi;
		}
		a = num/den;	b = yavg - a*xavg;
		printf("Least-squares omitting 20 smallest M(p) gives slope = %8.4f, y-intercept = %8.4f\n",a,b);

	// Now do another linear regression, this time omitting smallest 30 M-exponents:
		for(i = 30, yavg = 0.0; i < n; i++) {
			yavg += y[i];
		}
		xavg = (1.0 + 30 + n)/2;	yavg /= (n-30);	// X-avg uses unit offset
		printf("Omitting 30 smallest M(p): Sample size = %u, xavg = %8.4f, yavg = %8.4f\n",n-30,xavg,yavg);
		for(i = 30, num = den = 0.0; i < n; i++) {
			xi = i+1;	num += (y[i] - yavg)*xi;	den += (xi - xavg)*xi;
		}
		a = num/den;	b = yavg - a*xavg;
		printf("Least-squares omitting 30 smallest M(p) gives slope = %8.4f, y-intercept = %8.4f\n",a,b);

	// Lastly, do another linear regression, this time omitting smallest 40 M-exponents:
		for(i = 40, yavg = 0.0; i < n; i++) {
			yavg += y[i];
		}
		xavg = (1.0 + 40 + n)/2;	yavg /= (n-40);	// X-avg uses unit offset
		printf("Omitting 40 smallest M(p): Sample size = %u, xavg = %8.4f, yavg = %8.4f\n",n-40,xavg,yavg);
		for(i = 40, num = den = 0.0; i < n; i++) {
			xi = i+1;	num += (y[i] - yavg)*xi;	den += (xi - xavg)*xi;
		}
		a = num/den;	b = yavg - a*xavg;
		printf("Least-squares omitting 40 smallest M(p) gives slope = %8.4f, y-intercept = %8.4f\n",a,b);
	}

	void test_mp_pm1_smooth(uint32 p)
	{
		double u_so_smoove, logf, ilogn, dtmp;
		const double ln2 = log(2.0);
		uint32 nprime = 1000, pm_gap = 10000, thresh = 100000;
		uint32 curr_p,fbase2psp_idx,i,ihi,itmp32,j,jlo,jhi,k,max_diff,m,nfac,np,pm1;
		const uint32 pdiff_8[8] = {2,1,2,1,2,3,1,3}, pdsum_8[8] = { 0, 2, 6, 8,12,18,20,26};
		// Compact table storing the (difference/2) between adjacent odd primes.
		unsigned char *pdiff = (unsigned char *)calloc(nprime, sizeof(unsigned char));	// 1000 primes is plenty for this task
		// Struct used for storing smoothness data ... make big enough to store all primes in [p - pm_gap, p + pm_gap] with a safety factor
		struct psmooth sdat;
		// .../10 here is an approximation based on prime density for primes > 100000;
		// note the code uses an interval [p-pm_gap, p+pm_gap], i.e. of length 2*pm_gap, so the calloc needs to be twice pm_gap/10:
		struct psmooth*psmooth_vec = (struct psmooth *)calloc(2*pm_gap/10, sizeof(struct psmooth));

		/* Init first few diffs between 3/5, 5/7, 7/11, so can start loop with curr_p = 11 == 1 (mod 10), as required by twopmodq32_x8(): */
		pdiff[1] = pdiff[2] = 1;
		ihi = curr_p = 11;
		/* Process chunks of length 30, starting with curr_p == 11 (mod 30). Applying the obvious divide-by-3,5 mini-sieve,
		we have 8 candidates in each interval: curr_p + [ 0, 2, 6, 8,12,18,20,26].
		For example: curr_p = 11 gives the 8 candidates: 11,13,17,19,23,29,31,37.
		*/
		fbase2psp_idx = 0;	// Index to next-expected Fermat base-2 pseudoprime in the precomputed table
		for(i = 3; i < nprime; curr_p += 30) {
			/* Make sure (curr_p + 29) < 2^32: */
			if(curr_p > 0xffffffe3) {
				fprintf(stderr,"curr_p overflows 32 bits!");
				nprime = i;
				break;
			}
			/* Do a quick Fermat base-2 compositeness test before invoking the more expensive mod operations: */
			itmp32 = twopmodq32_x8(curr_p, curr_p+ 2, curr_p+ 6, curr_p+ 8, curr_p+12, curr_p+18, curr_p+20, curr_p+26);
			for(j = 0; j < 8; ++j) {
				if((itmp32 >> j)&0x1)	// It's a PRP, so check against the table of known pseudoprimes and
				{						// (if it's not a PSP) init for the next gap
					ASSERT(HERE, curr_p <= fbase2psp[fbase2psp_idx],"Error in pseudoprime sieve");
					if((curr_p + pdsum_8[j]) == fbase2psp[fbase2psp_idx]) {	/* It's a base-2 pseudoprime */
						++fbase2psp_idx;
						pdiff[i] += pdiff_8[j];
						continue;
					} else {	/* It's prime - add final increment to current pdiff[i] and then increment i: */
						ihi = (curr_p + pdsum_8[j]);
						pdiff[i] += pdiff_8[j];
						if(pdiff[i] > max_diff) {
							max_diff = pdiff[i];
						#if DBG_SIEVE
							printf("pdiff = %d at curr_p = %u\n", 2*max_diff,ihi);
						#endif
						}
						if(++i == nprime)
							break;
					}
				} else
					pdiff[i] += pdiff_8[j];
			}
			continue;
		}
		printf("Using first %u odd primes; max gap = %u\n",nprime,2*max_diff);
		printf("max sieving prime = %u\n",ihi);

		ASSERT(HERE, p > thresh, "Mersenne prime exponent must be larger that allowable threshold!");
		ASSERT(HERE, twopmodq32(p-1, p) == 1, "p fails base-2 fprp test!");
		np = 0;	// #primes in the current p-centered cohort
		// find N primes < and > p, compute smoothness norm based on p-1 factorization for each, store each [p,snorm] pair
		fbase2psp_idx = 0;	// Index to next-expected Fermat base-2 pseudoprime in the precomputed table
		jlo = p-pm_gap; jhi = p+pm_gap;
		// Find right tarting slot in base-2 pseudoprime table:
		while(fbase2psp[fbase2psp_idx] < jlo)
			++fbase2psp_idx;
		for(j = jlo; j <= jhi; j+=2) {
			// Do base-2 fprp test of j:
			if(!twopmodq32(j-1,j))
				continue;
			if(j == fbase2psp[fbase2psp_idx]) {	// It's a base-2 pseudoprime
				++fbase2psp_idx;
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
				nfac++;	pm1 >>= 1;	dtmp = logf*ilogn;	u_so_smoove += dtmp*dtmp;
			}
			if(nfac > 1) {
				printf("2^%u",nfac);
			} else {
				printf("2");
			}
			curr_p = 3;
			for(m = 0; m < nprime; m++) {
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
	}	// test_mp_pm1_smooth()

#endif	// ENABLE_MPRIME_PM1_SMOOTH

#if defined(USE_GPU) && defined(__CUDACC__)

	// Simple vector-add test function:
	__global__ void VecAdd(float* A, float* B, float* C, int N)
	{
		int i = blockDim.x * blockIdx.x + threadIdx.x;
		// Uncomment if() to Print basic info about threads ... keep I/O reasonable, only do so for first 10 of each batch of 2^18:
		if(i%0x3ffff < 10)
			printf("GPU block %d[dim %d], thread %d ==> seq-thread %d [i%0x3ffff = %d]... \n", blockIdx.x, blockDim.x, threadIdx.x, i,i%0x3ffff);
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
		float *h_A = (float*)malloc(size), *h_B = (float*)malloc(size), *h_C = (float*)malloc(size);
		// Initialize input vectors
		for(i = 0; i < N; ++i) {
			*(h_A+i) = i;
			*(h_B+i) = i*0.1;
		}
		// Allocate vectors in device memory
		float *d_A, *d_B, *d_C;
		cudaMalloc(&d_A, size);	cudaMalloc(&d_B, size);	cudaMalloc(&d_C, size);
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
		cudaFree(d_A);	cudaFree(d_B);	cudaFree(d_C);
		// Debug-print results sample:
		for(i = 0; i < 10; ++i) {
			printf("i = %d: Sum = %10.2f\n", i, *(h_C+i));
		}
		printf("...\n");
		for(i = 10; i > 0; --i) {
			printf("i = %d: Sum = %10.2f\n", N-i, *(h_C+N-i));
		}
	}

#ifdef USE_FMADD

	__global__ void VecMul50x50_exact(double*a, double*b, double*lo, double*hi)
	{
		int i = blockDim.x * blockIdx.x + threadIdx.x;
		// Exact product a*b = lo + hi:
		hi[i] = a[i] * b[i];
		lo[i] = fma(a[i],b[i], -hi[i]);
	}

	// Host code for the VecAdd test:
	void cudaVecMul50x50Test()
	{
		int i, N = 1024*1024, pow2;
		double pow2_dmult;
		uint64 iax,iay,ialo,iahi;
		size_t size = N * sizeof(double);
		// Allocate input vectors h_A and h_B in host memory
		double *h_A = malloc(size), *h_B = malloc(size), *h_C = malloc(size), *h_D = malloc(size);
		// Allocate vectors in device memory
		double*d _A, *d_B, *d_C, *d_D;
		cudaMalloc(&d_A, size), cudaMalloc(&d_B, size), cudaMalloc(&d_C, size), cudaMalloc(&d_D, size);
		// Assumes rng_isaac_init() has already been called on entry
		pow2_dmult = TWO50FLOAT;	// This must match the loop-starting value of pow2:
		for(pow2 = 50; pow2 < 54; ++pow2)	// Only makes sense to test up the #bits in an IEEE-double mantissa: Any larger and we start losing
		{									// LSBs (I.e. the test may 'succeed' for pow2 > 53, but is only testing the equivalent of pow2 = 53.)
		//	printf("Testing CUDA fma_dmult for %d bits, dmult = %f:\n",pow2,pow2_dmult);
			// Initialize input vectors
			for(i = 0; i < N; ++i) {
				// Input multiplicands in [-2^pow2, +2^pow2]:
				*(h_A+i) = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );
				*(h_B+i) = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );
			}
			// Copy vectors from host memory to device memory
			cudaMemcpy(d_A, h_A, size, cudaMemcpyHostToDevice);
			cudaMemcpy(d_B, h_B, size, cudaMemcpyHostToDevice);
			// Invoke kernel
			int threadsPerBlock = 256;
			int blocksPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock;
			VecMul50x50_exact<<<blocksPerGrid, threadsPerBlock>>>(d_A,d_B, d_C,d_D);
			// Copy result from device memory to host memory
			// h_C contains the result in host memory
			cudaMemcpy(h_C, d_C, size, cudaMemcpyDeviceToHost);
			cudaMemcpy(h_D, d_D, size, cudaMemcpyDeviceToHost);
			// Debug-print results sample:
			for(i = 0; i < N; ++i) {
				iax = ABS(h_A[i]);	iay = ABS(h_B[i]);
			#ifdef MUL_LOHI64_SUBROUTINE
				MUL_LOHI64(iax,iay,&ialo,&iahi);
			#else
				MUL_LOHI64(iax,iay, ialo, iahi);
			#endif
			//	printf("I = %d: x = %f; y = %f; hi,lo = %f,%f\n",i, h_A[i],h_B[i],h_D[i],h_C[i]);
				if(cmp_fma_lohi_vs_exact(h_A[i],h_B[i],h_D[i],h_C[i], iax,iay,iahi,ialo)) {
					printf("ERROR: pow2 = %d, I = %d, outputs differ!\n",pow2,i);
					ASSERT(HERE, 0, "fma_dmult tests failed!");
				}
			}	// i-loop
			pow2_dmult *= 2;
		}	// pow2-loop
		// Free device memory
		cudaFree(d_A);	cudaFree(d_B);	cudaFree(d_C);	cudaFree(d_D);
		printf("CUDA fma_dmult_tests completed successfully!\n");
	}

#endif

	/**********************************************************************************************/

	#include "factor.h"
	#include "fac_test_dat64.h"
	#include "fac_test_dat96.h"
	#include "twopmodq80.h"

	// Host code for the 64-bit VecModpow test:
	void cudaVecModpowTest64()
	{
		int i, nelt64;
		uint64 p, pshift, k,q;
		uint32 start_index, zshift, j,jshift, leadb;
		const uint32 N = 1<<10;
		double dbl, rnd;
		// Allocate input vectors (which take the TF p/pshift/zshift/start_index/k data on input) in host memory:
		uint64 *h_p = malloc(N<<3), *h_pshft = malloc(N<<3), *h_k = malloc(N<<3);
		uint32 *h_zshft = malloc(N<<2), *h_stidx = malloc(N<<2);

		// Do counting pass to set nelt64, the number of 64-bit test data available:
		for(i = 0; i < N; ++i) {
			if(0 == fac64[i].p) break;
		}
		nelt64 = i;
		for(i = 0; i < N; ++i) {
			if(i < nelt64) {
				p = fac64[i].p;
				q = fac64[i].q;
			} else {	// Fill in any remaining slots with 63-bit test data. of which we know we have > (1<<10):
				p = fac63[i-nelt64].p;
				q = fac63[i-nelt64].q;
	//if((i-nelt64) < 10)printf("p[%3d] = %u: q = %llu ... ",i, p, q);
			}
			ASSERT(HERE, p != 0, "p must be nonzero!");
			// Compute auxiliary TF data:
			pshift = p + 64;
			jshift = leadz64(pshift);
			/* Extract leftmost 6 bits of pshift and subtract from 64: */
			leadb = ((pshift<<jshift) >> 58);
			start_index = 64-jshift-6;
			zshift = 63 - leadb;
			zshift <<= 1;				/* Doubling the shift count here takes cares of the first SQR_LOHI */
			pshift = ~pshift;

			// Compute factor k using fast DP math. Integer-truncation-on-store should obviate the need
			// to subtract 1 from q, and (double)q is only accurate to 53 bits to begin with):
			dbl = (double)q;
			dbl /= (2.0*p);
			rnd = DNINT(dbl);
			k = (uint64)rnd;
			ASSERT(HERE, k*(p<<1)+1 == q, "k computed incorrectly!");
			*(h_p     + i) = p          ;	*(h_pshft + i) = pshift     ;	*(h_k + i) = k;
			*(h_zshft + i) = zshift     ;	*(h_stidx + i) = start_index;
//	printf("p[%3d] = %u: pshift = %8u, zshift = %8u, stidx = %2u, k = %llu\n",i, p, pshift, zshift, start_index, k);
		}
		printf("Testing %d = %d 64-bit and %d 63-bit known-factors...",N,nelt64,N-nelt64);

		// Initialize output vector (resulting 2^p mod q, in binary "is factor?" form) in host memory:
		uint8*  h_B = (uint8 *)malloc(N);	// Until impl packed-bitmap scheme for device code return values, use byte array for return values
		for(i = 0; i < N; ++i) {
			*(h_B+i) = 0;
		}
//printf("Host code: p = %u, pshift = %u, k = %llu, zshift = %u, start_index = %u\n", p,pshift,h_A[0],zshift,start_index);
		// Allocate vectors in device memory
		uint64 *d_p,*d_pshft,*d_k;
		uint32 *d_zshft,*d_stidx;
		cudaMalloc(&d_p    , N<<3);	cudaMalloc(&d_pshft, N<<3);	cudaMalloc(&d_k    , N<<3);
		cudaMalloc(&d_zshft, N<<2);	cudaMalloc(&d_stidx, N<<2);
		uint8 * d_B;
		cudaMalloc(&d_B, N);
		// Copy vectors from host memory to device memory
		cudaMemcpy(d_p    , h_p    , N<<3, cudaMemcpyHostToDevice);
		cudaMemcpy(d_pshft, h_pshft, N<<3, cudaMemcpyHostToDevice);
		cudaMemcpy(d_zshft, h_zshft, N<<2, cudaMemcpyHostToDevice);
		cudaMemcpy(d_stidx, h_stidx, N<<2, cudaMemcpyHostToDevice);
		cudaMemcpy(d_k    , h_k    , N<<3, cudaMemcpyHostToDevice);
		// Do we need to copy the as-yet-uninited output vector to (or just from) the device?
	//	cudaMemcpy(d_B, h_B, N   , cudaMemcpyHostToDevice);

		// Invoke kernel
		int threadsPerBlock = 256;
		int blocksPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock;

//	printf("VecModpow64<<< %d, %d >>> with %d 64-bit known-factors:\n", blocksPerGrid, threadsPerBlock, N);
		VecModpow64<<<blocksPerGrid, threadsPerBlock>>>(d_p,d_pshft,d_zshft,d_stidx,d_k, d_B, N);

		// Copy result from device memory to host memory
		// h_B contains the result in host memory
		cudaMemcpy(h_B, d_B, N, cudaMemcpyDeviceToHost);
		// Free device memory
		cudaFree(d_p);	cudaFree(d_pshft);	cudaFree(d_zshft);	cudaFree(d_stidx);	cudaFree(d_k);	cudaFree(d_B);

		// Reference computation to Test GPU results:
		for(i = 0; i < N; ++i) {
			p = *(h_p + i);
			k = *(h_k + i);	q = k*(p<<1)+1;
			j = (uint32)twopmodq64((uint64)p, q);
			if((j != 1) || (*(h_B + i) != 1)) {
				printf("cudaVecModpowTest64: Mismatch between Ref and GPU result:\n");
				printf("res[%d] = %d [ref = %d] = 2^p - 1 (mod q) with (p,q) = %u, %llu\n", i,*(h_B + i), j,p,q);
				ASSERT(HERE, 0, "cudaVecModpowTest64 failed!");
			}
		}
		printf("cudaVecModpowTest64 with %d test (p,q) pairs succeeded!\n",N);
	}

	// Host code for the simpler VecModpow test, same 78-bit [p,q] pair for each thread:
	void cudaVecModpowTest78_0()
	{
		int i;
		uint64 p, pshift, k;
		uint32 start_index, zshift, j, leadb;
		uint32 N = 1<<10;
		// Allocate input vectors (which take the TF p/pshift/zshift/start_index/k data on input) in host memory:
		uint64 *h_p = malloc(N<<3), *h_pshft = malloc(N<<3), *h_k = malloc(N<<3);
		uint32 *h_zshft = malloc(N<<2), *h_stidx = malloc(N<<2);

		p = 16727479;
		k = 7946076362870052ull;
		// Compute auxiliary TF data:
		pshift = p + 78;
		j = leadz64(pshift);
		/* Extract leftmost 7 bits of pshift (if > 77, use the leftmost 6) and subtract from 96: */
		leadb = ((pshift<<j) >> 57);
		if(leadb > 77) {
			leadb >>= 1;
			start_index =  64-j-6;	/* Use only the leftmost 6 bits */
		} else {
			start_index =  64-j-7;
		}
		zshift = 77 - leadb;
		zshift <<= 1;				/* Doubling the shift count here takes cares of the first SQR_LOHI */
		pshift = ~pshift;

		// Copy to all N vector-input-data:
		for(i = 0; i < N; ++i) {
			*(h_p     + i) = p          ;	*(h_pshft + i) = pshift     ;	*(h_k + i) = k;
			*(h_zshft + i) = zshift     ;	*(h_stidx + i) = start_index;
		}
		printf("Testing %d 78-bit known-factors...",N);

		// Initialize output vector (resulting 2^p mod q, in binary "is factor?" form) in host memory:
		uint8*  h_B = (uint8 *)malloc(N);	// Until impl packed-bitmap scheme for device code return values, use byte array for return values
		for(i = 0; i < N; ++i) {
			*(h_B+i) = 0;
		}
//printf("Host code: p = %u, pshift = %u, k = %llu, zshift = %u, start_index = %u\n", p,pshift,h_A[0],zshift,start_index);
		// Allocate vectors in device memory
		uint64 *d_p,*d_pshft,*d_k;
		uint32 *d_zshft,*d_stidx;
		cudaMalloc(&d_p    , N<<3);	cudaMalloc(&d_pshft, N<<3);	cudaMalloc(&d_k    , N<<3);
		cudaMalloc(&d_zshft, N<<2);	cudaMalloc(&d_stidx, N<<2);
		uint8 * d_B;
		cudaMalloc(&d_B, N);
		// Copy vectors from host memory to device memory
		cudaMemcpy(d_p    , h_p    , N<<3, cudaMemcpyHostToDevice);
		cudaMemcpy(d_pshft, h_pshft, N<<3, cudaMemcpyHostToDevice);
		cudaMemcpy(d_zshft, h_zshft, N<<2, cudaMemcpyHostToDevice);
		cudaMemcpy(d_stidx, h_stidx, N<<2, cudaMemcpyHostToDevice);
		cudaMemcpy(d_k    , h_k    , N<<3, cudaMemcpyHostToDevice);
		// Do we need to copy the as-yet-uninited output vector to (or just from) the device?
	//	cudaMemcpy(d_B, h_B, N   , cudaMemcpyHostToDevice);

		// Invoke kernel
		int threadsPerBlock = 256;
		int blocksPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock;

//	printf("VecModpow78<<< %d, %d >>> with %d copies of same 78-bit known-factor [blocksPerGrid = %d, threadsPerBlock = %d]:\n", blocksPerGrid, threadsPerBlock, N);
		VecModpow78<<<blocksPerGrid, threadsPerBlock>>>(d_p,d_pshft,d_zshft,d_stidx,d_k, d_B, N);
//	printf("GPU_TF78<<< %d, %d >>> with %d copies of same 78-bit known-factor:\n", blocksPerGrid, threadsPerBlock, N);
//		GPU_TF78<<<blocksPerGrid, threadsPerBlock>>>(p,pshift,zshift,start_index, d_k, d_B, N);

		// Copy result from device memory to host memory
		// h_B contains the result in host memory
		cudaMemcpy(h_B, d_B, N, cudaMemcpyDeviceToHost);
		// Free device memory
		cudaFree(d_p);	cudaFree(d_pshft);	cudaFree(d_zshft);	cudaFree(d_stidx);	cudaFree(d_k);	cudaFree(d_B);

		// Reference computation:
		j = (uint32)twopmodq78_3WORD_DOUBLE((uint64)p, k);
		ASSERT(HERE, (j == 1), "cudaVecModpowTest78_0 ref-comp failed!");
		// Test GPU results:
		for(i = 0; i < N; ++i) {
			if(*(h_B + i) != 1) {
				printf("cudaVecModpowTest78_0: Mismatch between Ref and GPU result:\n");
				printf("res[%d] = %d [ref = %d] = 2^p - 1 (mod q) with (p,k) = %u, %llu\n", i,*(h_B + i), j,p,k);
				ASSERT(HERE, *(h_B + i) == 1, "cudaVecModpowTest78_0 failed!");
			}
		}
		printf("cudaVecModpowTest78_0 with %d test (p,q) pairs succeeded!\n",N);
	}

	// Host code for the 78-bit VecModpow test:
	void cudaVecModpowTest78()
	{
		int i;
		uint64 p, pshift, k;
		uint32 start_index, zshift, j, leadb;
		uint32 N = 1<<10,nelts;
		uint96 q96;
		double dbl, rnd;
		// Allocate input vectors (which take the TF p/pshift/zshift/start_index/k data on input) in host memory:
		uint64 *h_p = malloc(N<<3), *h_pshft = malloc(N<<3), *h_k = malloc(N<<3);
		uint32 *h_zshft = malloc(N<<2), *h_stidx = malloc(N<<2);

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
			j = leadz64(pshift);
			/* Extract leftmost 7 bits of pshift (if > 77, use the leftmost 6) and subtract from 96: */
			leadb = ((pshift<<j) >> 57);
			if(leadb > 77) {
				leadb >>= 1;
				start_index =  64-j-6;	/* Use only the leftmost 6 bits */
			} else {
				start_index =  64-j-7;
			}
			zshift = 77 - leadb;
			zshift <<= 1;				/* Doubling the shift count here takes cares of the first SQR_LOHI */
			pshift = ~pshift;

			// Compute factor k using fast DP math. Integer-truncation-on-store should obviate the need
			// to subtract 1 from q, and (double)q is only accurate to 53 bits to begin with):
			dbl = (double)q96.d0 + (double)q96.d1*TWO64FLOAT;
			dbl /= (2.0*p);
			rnd = DNINT(dbl);
			k = (uint64)rnd;
			*(h_p     + nelts) = p          ;	*(h_pshft + nelts) = pshift     ;	*(h_k + nelts) = k;
			*(h_zshft + nelts) = zshift     ;	*(h_stidx + nelts) = start_index;
//	printf("p[%3d] = %u: pshift = %8u, zshift = %8u, stidx = %2u, k = %llu\n",nelts, p, pshift, zshift, start_index, k);
			++nelts;
		}
		printf("Testing %d 78-bit known-factors...",nelts);
		// "Fill in" remaining slots with copy of same datum used in cudaVecModpowTest78_0:
		p = 16727479;
		k = 7946076362870052ull;
		// Compute auxiliary TF data:
		pshift = p + 78;
		j = leadz64(pshift);
		/* Extract leftmost 7 bits of pshift (if > 77, use the leftmost 6) and subtract from 96: */
		leadb = ((pshift<<j) >> 57);
		if(leadb > 77) {
			leadb >>= 1;
			start_index =  64-j-6;	/* Use only the leftmost 6 bits */
		} else {
			start_index =  64-j-7;
		}
		zshift = 77 - leadb;
		zshift <<= 1;				/* Doubling the shift count here takes cares of the first SQR_LOHI */
		pshift = ~pshift;
		// Copy to all still-uninited vector-input-data:
		for(i = nelts; i < N; ++i) {
			*(h_p     + i) = p          ;	*(h_pshft + i) = pshift     ;	*(h_k + i) = k;
			*(h_zshft + i) = zshift     ;	*(h_stidx + i) = start_index;
		}

		// Initialize output vector (resulting 2^p mod q, in binary "is factor?" form) in host memory:
		uint8*  h_B = (uint8 *)malloc(N);	// Until impl packed-bitmap scheme for device code return values, use byte array for return values
		for(i = 0; i < N; ++i) {
			*(h_B+i) = 0;
		}
//printf("Host code: p = %u, pshift = %u, k = %llu, zshift = %u, start_index = %u\n", p,pshift,h_A[0],zshift,start_index);
		// Allocate vectors in device memory
		uint64 *d_p,*d_pshft,*d_k;
		uint32 *d_zshft,*d_stidx;
		cudaMalloc(&d_p    , N<<3);	cudaMalloc(&d_pshft, N<<3);	cudaMalloc(&d_k    , N<<3);
		cudaMalloc(&d_zshft, N<<2);	cudaMalloc(&d_stidx, N<<2);
		uint8 * d_B;
		cudaMalloc(&d_B, N);
		// Copy vectors from host memory to device memory
		cudaMemcpy(d_p    , h_p    , N<<3, cudaMemcpyHostToDevice);
		cudaMemcpy(d_pshft, h_pshft, N<<3, cudaMemcpyHostToDevice);
		cudaMemcpy(d_zshft, h_zshft, N<<2, cudaMemcpyHostToDevice);
		cudaMemcpy(d_stidx, h_stidx, N<<2, cudaMemcpyHostToDevice);
		cudaMemcpy(d_k    , h_k    , N<<3, cudaMemcpyHostToDevice);
		// Do we need to copy the as-yet-uninited output vector to (or just from) the device?
	//	cudaMemcpy(d_B, h_B, N   , cudaMemcpyHostToDevice);

		// Invoke kernel
		int threadsPerBlock = 256;
		int blocksPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock;

//		printf("VecModpow78 with %d 78-bit known-factors [blocksPerGrid = %d, threadsPerBlock = %d]:\n",nelts, blocksPerGrid, threadsPerBlock);
		VecModpow78<<<blocksPerGrid, threadsPerBlock>>>(d_p,d_pshft,d_zshft,d_stidx,d_k, d_B, nelts);

		// Copy result from device memory to host memory
		// h_B contains the result in host memory
		cudaMemcpy(h_B, d_B, N, cudaMemcpyDeviceToHost);
		// Free device memory
		cudaFree(d_p);	cudaFree(d_pshft);	cudaFree(d_zshft);	cudaFree(d_stidx);	cudaFree(d_k);	cudaFree(d_B);

		// Reference computation:
		// Test GPU results:
		for(i = 0; i < nelts; ++i) {
			p = *(h_p + i);
			k = *(h_k + i);
			j = (uint32)twopmodq78_3WORD_DOUBLE((uint64)p, k);
			if((j != 1) || (*(h_B + i) != 1)) {
				printf("cudaVecModpowTest78: Mismatch between Ref and GPU result:\n");
				printf("res[%d] = %d [ref = %d] = 2^p - 1 (mod q) with (p,k) = %u, %llu\n", i,*(h_B + i), j,p,k);
				ASSERT(HERE, 0, "cudaVecModpowTest78 failed!");
			}
		}
		printf("cudaVecModpowTest78 with %d test (p,q) pairs succeeded!\n",nelts);
	}

	// Host code for the 96-bit VecModpow test:
	void cudaVecModpowTest96()
	{
		int i;
		uint64 p, pinv, pshift, k, hi64;
		uint32 tmp32, start_index, zshift, j, leadb;
		uint32 N = 1<<10,nelts;
		uint96 q96,x96,pinv96;
		// Allocate input vectors (which take the TF p/pshift/zshift/start_index/k data on input) in host memory:
		uint64 *h_p = malloc(N<<3), *h_pshft = malloc(N<<3), *h_k = malloc(N<<3);
		uint32 *h_zshft = malloc(N<<2), *h_stidx = malloc(N<<2);

		for(i = 0, nelts = 0; i < N; ++i) {
			p = fac96[i].p;
			if(p == 0) {
				break;
			}
			q96.d1 = fac96[i].d1; q96.d0 = fac96[i].d0;
			// Good to go - compute auxiliary TF data:
			pshift = p + 96;
			j = leadz64(pshift);
			// Extract leftmost 7 bits of pshift (if > 95, use the leftmost 6) and subtract from 96:
			leadb = ((pshift<<j) >> 57);
			if(leadb > 95) {
				leadb >>= 1;
				start_index =  64-j-6;	// Use only the leftmost 6 bits
			} else {
				start_index =  64-j-7;
			}
			zshift = 95 - leadb;
			zshift <<= 1;				// Doubling the shift count here takes cares of the first SQR_LOHI
			pshift = ~pshift;

			/* To find the quotient k = (q-1)/(2*p), which may be > 64 bits, use mod-inverse with base 2^96 arithmetic.
			Since the Newtonian mod-inverse algorithm only works for odd inputs, instead of finding (q-1)/(2*p), we find ((q-1)/2)/p.
			First, find inverse (mod 2^96) of p in preparation for modular multiply. See twopmodq96 for an explanation of this:
			*/
			pinv = (p +p +p) ^ 2;
			for(j = 0; j < 3; j++) {
				tmp32 = p * pinv;
				pinv = pinv*(2 - tmp32);
			}
			// One more iteration using uint64 math to get 64-bit inverse:
			pinv96.d0 = (uint64)pinv;	pinv96.d1 = (uint64)0;
			hi64 = (uint64)p * pinv96.d0;
			pinv96.d0 = pinv96.d0*((uint64)2 - hi64);
			// pinv96 has 96 bits, but only the upper 64 get modified here:
		#ifdef MUL_LOHI64_SUBROUTINE
			pinv96.d1 = -pinv96.d0*__MULH64((uint64)p, pinv96.d0);
		#else
			MULH64((uint64)p, pinv96.d0, hi64);
			pinv96.d1 = -pinv96.d0*hi64;
		#endif
			// k is simply the bottom 96 bits of ((q-1)/2)*pinv96:
			x96.d0	= ((q96.d0-1) >> 1) + ((uint64)q96.d1 << 63);	x96.d1	= (q96.d1 >> 1);	// (q-1)/2
			MULL96(x96, pinv96, x96);
			k = x96.d0;
			// Skip any (p,q) pair for which the k > 2^64:
			if(x96.d1 != 0) {	// x128 holds k
			//	printf("Warning: k > 2^64 detected for (p,q) = %u,[%u*2^64 + %llu] ... skipping this datum.\n",p,q96.d1,q96.d0);
				continue;
			}
			*(h_p     + nelts) = p          ;	*(h_pshft + nelts) = pshift     ;	*(h_k + nelts) = k;
			*(h_zshft + nelts) = zshift     ;	*(h_stidx + nelts) = start_index;
//	printf("p[%3d] = %u: pshift = %8u, zshift = %8u, stidx = %2u, k = %llu\n",nelts, p, pshift, zshift, start_index, k);
			++nelts;
		}
		printf("Testing %d 96-bit known-factors...",nelts);
		// "Fill in" remaining slots with copy of same datum used in cudaVecModpowTest96_0:
		p = 16727479;
		k = 7946076362870052ull;
		// Compute auxiliary TF data:
		pshift = p + 96;
		j = leadz32(pshift);
		/* Extract leftmost 7 bits of pshift (if > 85, use the leftmost 6) and subtract from 96: */
		leadb = ((pshift<<j) >> 57);
		if(leadb > 95) {
			leadb >>= 1;
			start_index =  64-j-6;	/* Use only the leftmost 6 bits */
		} else {
			start_index =  64-j-7;
		}
		zshift = 95 - leadb;
		zshift <<= 1;				/* Doubling the shift count here takes cares of the first SQR_LOHI */
		pshift = ~pshift;
		// Copy to all still-uninited vector-input-data:
		for(i = nelts; i < N; ++i) {
			*(h_p     + i) = p          ;	*(h_pshft + i) = pshift     ;	*(h_k + i) = k;
			*(h_zshft + i) = zshift     ;	*(h_stidx + i) = start_index;
		}

		// Initialize output vector (resulting 2^p mod q, in binary "is factor?" form) in host memory:
		uint8*  h_B = (uint8 *)malloc(N);	// Until impl packed-bitmap scheme for device code return values, use byte array for return values
		for(i = 0; i < N; ++i) {
			*(h_B+i) = 0;
		}
//printf("Host code: p = %u, pshift = %u, k = %llu, zshift = %u, start_index = %u\n", p,pshift,h_A[0],zshift,start_index);
		// Allocate vectors in device memory
		uint64 *d_p,*d_pshft,*d_k;
		uint32 *d_zshft,*d_stidx;
		cudaMalloc(&d_p    , N<<3);	cudaMalloc(&d_pshft, N<<3);	cudaMalloc(&d_k    , N<<3);
		cudaMalloc(&d_zshft, N<<2);	cudaMalloc(&d_stidx, N<<2);
		uint8 * d_B;
		cudaMalloc(&d_B, N);
		// Copy vectors from host memory to device memory
		cudaMemcpy(d_p    , h_p    , N<<3, cudaMemcpyHostToDevice);
		cudaMemcpy(d_pshft, h_pshft, N<<3, cudaMemcpyHostToDevice);
		cudaMemcpy(d_zshft, h_zshft, N<<2, cudaMemcpyHostToDevice);
		cudaMemcpy(d_stidx, h_stidx, N<<2, cudaMemcpyHostToDevice);
		cudaMemcpy(d_k    , h_k    , N<<3, cudaMemcpyHostToDevice);
		// Do we need to copy the as-yet-uninited output vector to (or just from) the device?
	//	cudaMemcpy(d_B, h_B, N   , cudaMemcpyHostToDevice);

		// Invoke kernel
		int threadsPerBlock = 256;
		int blocksPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock;

//		printf("VecModpow96 with %d 96-bit known-factors [blocksPerGrid = %d, threadsPerBlock = %d]:\n",nelts, blocksPerGrid, threadsPerBlock);
		VecModpow96<<<blocksPerGrid, threadsPerBlock>>>(d_p,d_pshft,d_zshft,d_stidx,d_k, d_B, nelts);

		// Copy result from device memory to host memory
		// h_B contains the result in host memory
		cudaMemcpy(h_B, d_B, N, cudaMemcpyDeviceToHost);
		// Free device memory
		cudaFree(d_p);	cudaFree(d_pshft);	cudaFree(d_zshft);	cudaFree(d_stidx);	cudaFree(d_k);	cudaFree(d_B);

		// Reference computation:
		// Test GPU results:
		for(i = 0; i < nelts; ++i) {
			p = *(h_p + i);
			k = *(h_k + i);
			q96 = twopmodq96((uint64)p, k);
			j = (q96.d1 == 0) && (q96.d0 == 1);
			if((j != 1) || (*(h_B + i) != 1)) {
				printf("cudaVecModpowTest96: Mismatch between Ref and GPU result:\n");
				printf("res[%d] = %d [ref = %d] = 2^p - 1 (mod q) with (p,k) = %u, %llu\n", i,*(h_B + i), j,p,k);
				ASSERT(HERE, 0, "cudaVecModpowTest96 failed!");
			}
		}
		printf("cudaVecModpowTest96 with %d test (p,q) pairs succeeded!\n",nelts);
	}

#endif

#undef PLATFORM_SKIP_RND_CONST_ENFORCEMENT

/*********************************************************************************/
/* Globals. Unless specified otherwise, these are declared in Mdata.h:           */
/*********************************************************************************/

/* These externs used in x86 builds by util.c:set_x87_fpu_params() yield 64-mantissa-bit
floating-point register mode (bits <9:8> = 3), with IEEE and truncating rounding mode set
via bits <11:10> = 0 and 3, respectively. The other 12 bits are identical to the MSVC defaults:
*/
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

int	file_valid(FILE*fp)
{
	return (fp && !ferror(fp) && !feof(fp));
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
	/*...time-related stuff	*/
	clock_t clock1, clock2;
	double tdiff;

	/* Various useful precomputed powers of 2 in floating-double form: */
	TWO25FLOAT = (double)0x02000000;				TWO25FLINV = 1.0/TWO25FLOAT;
	TWO26FLOAT = (double)0x04000000;				TWO26FLINV = 1.0/TWO26FLOAT;
	dbl = qfdbl(qfmul_pow2(QONE, -26));
	ASSERT(HERE, TWO26FLINV == dbl, "TWO26FLINV!");

	TWO13FLINV = qfdbl(qfmul_pow2(QONE, -13));

	TWO32FLOAT = (double)2.0*0x80000000;			TWO32FLINV = 1.0/TWO32FLOAT;
	TWO48FLOAT = (double)1.0*0x01000000*0x01000000;	TWO48FLINV = 1.0/TWO48FLOAT;
	TWO50FLOAT = (double)1.0*0x01000000*0x04000000;	TWO50FLINV = 1.0/TWO50FLOAT;
	TWO51FLOAT = (double)1.0*0x02000000*0x04000000;	TWO51FLINV = 1.0/TWO51FLOAT;
	TWO52FLOAT = (double)1.0*0x04000000*0x04000000;	TWO52FLINV = 1.0/TWO52FLOAT;
	TWO53FLOAT = (double)1.0*0x08000000*0x04000000;	TWO53FLINV = 1.0/TWO53FLOAT;
	TWO54FLOAT = (double)1.0*0x08000000*0x08000000;
	TWO63FLOAT = (double)2.0*0x80000000*0x80000000;
	TWO64FLOAT = (double)4.0*0x80000000*0x80000000;	TWO64FLINV = 1.0/TWO64FLOAT;

	/* Check qfloat routines (this call is also needed to init various qfloat global constants): */
	printf("INFO: testing qfloat routines...\n");
	qtest();	// 09/23/2012: Move to after above float-consts-inits because of the qfloat/mi64 routines which use those consts.

	/* Use qfloat routines to set the global floating-point constant 1/sqrt(2): */
	ASSERT(HERE, ISRT2 == qfdbl(QISRT2), "1/sqrt2 precision check failed!");
	ASSERT(HERE, SQRT2 == qfdbl(QSQRT2), "  sqrt2 precision check failed!");

#ifdef CPU_IS_X86	// May 2018: It seems I only found need to call this runtime CPU-mode setting in 32-bit x86 mode, not 64-bit. But had occasion
					// to fiddle w/rnd-mode in some x86_64 tests, so changed things so that the function is *defined* in both 32 and 64-bit modes.
	set_x87_fpu_params(FPU_64RND);
#endif
	// ewm [4. Aug 2014] - move below set_x87_fpu_params(), since need rnd-const set for any DNINT-using ref-computations in the GPU self-tests:
	print_host_info();
	check_nbits_in_types();

	/* Test wide-mul routines: */
	printf("INFO: testing IMUL routines...\n");
	ASSERT(HERE, test_mul() == 0, "test_mul() returns nonzero!");

	// Test certain aspects of SIMD functionality (aim is to expand this into a decently comprehensive
	// timing-test-of-key-SIMD-code-constructs suite):
#ifdef TEST_SIMD
	printf("INFO: Timing-testing selected FFT macros...\n");

  #if defined(USE_SSE2) && !defined(USE_AVX)	// 4-DFT is SSE2-only
//	ASSERT(HERE, test_radix4_dft() == 0, "test_radix4_dft() returns nonzero!");
  #endif

//	ASSERT(HERE, test_radix16_dft() == 0, "test_radix16_dft() returns nonzero!");

	#include "radix32_dif_dit_pass_asm.h"	// Commenting this out gives compile error
//	ASSERT(HERE, test_radix32_dft() == 0, "test_radix32_dft() returns nonzero!");

  #ifdef USE_AVX
//	ASSERT(HERE, test_simd_transpose_4x4() == 0, "test_simd_transpose_4x4() returns nonzero!");
  #endif
  #ifdef USE_AVX512
	ASSERT(HERE, test_simd_transpose_8x8() == 0, "test_simd_transpose_8x8() returns nonzero!");
exit(0);
  #endif
#endif

// Quick timings of various mi64 stuff:
#if 0

	rng_isaac_init(TRUE);
	uint32 i32,x32,bit, i,imax;
	uint64 i64;
	imax = 100000000;

	clock1 = clock();
	i64 = 0;
	for(i = 0; i < imax; i++) {
		i64 += rng_isaac_rand();
	}
	clock2 = clock();
	tdiff = (double)(clock2 - clock1);
	printf("Time for %u rng64 calls =%s\n",imax, get_time_str(tdiff));
	ASSERT(HERE, i64 != 0ull,"rng64 sum = 0!");

	clock1 = clock();
	i32 = 0;
	for(i = 0; i < imax; i++) {
		i64 = rng_isaac_rand();
		x32 = (uint32)i64;
		i32 += popcount32(x32);
		i32 += popcount32((uint32)(i64>>16));
		i32 += popcount32((uint32)(i64>>24));
		i32 += popcount32((uint32)(i64>>32));
	}
	clock2 = clock();
	tdiff = (double)(clock2 - clock1);
	printf("Time for %u [rng64 + 4*popcount32()] calls =%s\n",imax, get_time_str(tdiff));
	ASSERT(HERE, i32,"popcount32 sum = 0!");

	clock1 = clock();
	i32 = 0;
	for(i = 0; i < imax; i++) {
		i64 = rng_isaac_rand();
		i32 += popcount64(i64);
		i32 += popcount64(i64>>16);
		i32 += popcount64(i64>>24);
		i32 += popcount64(i64>>32);
	}
	clock2 = clock();
	tdiff = (double)(clock2 - clock1);
	printf("Time for %u [rng64 + 4*popcount64()] calls =%s\n",imax, get_time_str(tdiff));
	ASSERT(HERE, i32,"popcount64 sum = 0!");

	clock1 = clock();
	i32 = 0;
	for(i = 0; i < imax; i++) {
		i64 = rng_isaac_rand();
		x32 = (uint32)i64;
		i32 += leadz32(x32);
		i32 += leadz32((uint32)(i64>>16));
		i32 += leadz32((uint32)(i64>>24));
		i32 += leadz32((uint32)(i64>>32));
	}
	clock2 = clock();
	tdiff = (double)(clock2 - clock1);
	printf("Time for %u [rng64 + 4*leadz32()] calls =%s\n",imax, get_time_str(tdiff));
	ASSERT(HERE, i32,"leadz32 sum = 0!");

	clock1 = clock();
	i32 = 0;
	for(i = 0; i < imax; i++) {
		i64 = rng_isaac_rand();
		i32 += leadz64(i64);
		i32 += leadz64(i64>>16);
		i32 += leadz64(i64>>24);
		i32 += leadz64(i64>>32);
	}
	clock2 = clock();
	tdiff = (double)(clock2 - clock1);
	printf("Time for %u [rng64 + 4*leadz64()] calls =%s\n",imax, get_time_str(tdiff));
	ASSERT(HERE, i32,"leadz64 sum = 0!");

	clock1 = clock();
	i32 = 0;
	for(i = 0; i < imax; i++) {
		i64 = rng_isaac_rand();
		x32 = (uint32)i64;
		i32 += trailz32(x32);
		i32 += trailz32(x32);
		i32 += trailz32((uint32)(i64>>16));
		i32 += trailz32((uint32)(i64>>24));
		i32 += trailz32((uint32)(i64>>32));
	}
	clock2 = clock();
	tdiff = (double)(clock2 - clock1);
	printf("Time for %u [rng64 + 4*trailz32()] calls =%s\n",imax, get_time_str(tdiff));
	ASSERT(HERE, i32,"trailz32 sum = 0!");

	clock1 = clock();
	i32 = 0;
	for(i = 0; i < imax; i++) {
		i64 = rng_isaac_rand();
		i32 += trailz64(i64);
		i32 += trailz64(i64>>16);
		i32 += trailz64(i64>>24);
		i32 += trailz64(i64>>32);
	}
	clock2 = clock();
	tdiff = (double)(clock2 - clock1);
	printf("Time for %u [rng64 + 4*trailz64()] calls =%s\n",imax, get_time_str(tdiff));
	ASSERT(HERE, i32,"trailz64 sum = 0!");
exit(0);
	clock1 = clock();
	for(i = 0; i < imax; i++) {
		uint64 i64 = rng_isaac_rand();
		bit = (i64>>32) & 0x1f;	if(!bit) continue;
		x32 = (uint32)i64;
		int ii = ith_set_bit32(x32,bit);
		if(popcount32(x32) < bit)
			ASSERT(HERE, ii == -1, "[bit]th-bit specifier out of range!");
		else {
			uint32 tmp32 = x32 << (31-ii);
			ASSERT(HERE, tmp32 & 0x80000000,"ith_set_bit64 retval not actually set!");
			ASSERT(HERE, popcount32(tmp32) == bit, "ith_set_bit32 checksum fail!");
		}
	}
	clock2 = clock();
	tdiff = (double)(clock2 - clock1);
	printf("Time for %u ith_set_bit32() calls =%s\n",imax, get_time_str(tdiff));

	clock1 = clock2;
	for(i = 0; i < imax; i++) {
		uint64 i64 = rng_isaac_rand();
		bit = (i64>>32) & 0x3f;	if(!bit) continue;
		int ii = ith_set_bit64(i64,bit);
		if(popcount64(i64) < bit)
			ASSERT(HERE, ii == -1, "[bit]th-bit specifier out of range!");
		else {
			uint64 tmp64 = i64 << (63-ii);
			// Must cast result of AND to 32-bit here (via compare-vs-0) since ASSERT (expr) is 32-bit:
			ASSERT(HERE, (tmp64 & 0x8000000000000000ull) != 0,"ith_set_bit64 retval not actually set!");
			ASSERT(HERE, popcount64(tmp64) == bit, "ith_set_bit64 checksum fail!");
		}
	}
	clock2 = clock();
	tdiff = (double)(clock2 - clock1);
	printf("Time for %u ith_set_bit64() calls =%s\n",imax, get_time_str(tdiff));

	clock1 = clock2;
	uint64 iarr[4];
	for(i = 0; i < imax; i++) {
		iarr[0] = rng_isaac_rand();
		iarr[1] = rng_isaac_rand();
		iarr[2] = rng_isaac_rand();
		iarr[3] = rng_isaac_rand();
		bit = (iarr[0]>>32) & 0xff;	if(!bit) continue;
		int ii = mi64_ith_set_bit(iarr,bit,4);
		if(mi64_popcount(iarr,4) < bit)
			ASSERT(HERE, ii == -1, "[bit]th-bit specifier out of range!");
		else {
			mi64_shl(iarr,iarr,(255-ii),4);
			// Must cast result of AND to 32-bit here (via compare-vs-0) since ASSERT (expr) is 32-bit:
			ASSERT(HERE, (iarr[3] & 0x8000000000000000ull) != 0,"mi64_ith_set_bit64 retval not actually set!");
			ASSERT(HERE, mi64_popcount(iarr,4) == bit, "mi64_ith_set_bit64 checksum fail!");
		}
	}
	clock2 = clock();
	tdiff = (double)(clock2 - clock1);
	printf("Time for %u mi64_ith_set_bit() calls =%s\n",imax, get_time_str(tdiff));
exit(0);

#elif 0
	printf("INFO: Testing mi64_add speed...\n");
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
#ifdef ENABLE_MPRIME_PM1_SMOOTH
	uint32 p = 77232917;	// Jan 2016
	// Uncomment the specific subfunctions you desire:
	est_num_mp_in_interval(0,p);
//	compute_mers_best_fit();
	test_mp_pm1_smooth(p);
//	print_mp_dec(p);
	exit(0);
#endif

#ifdef MULTITHREAD

  #ifndef USE_PTHREAD
	#error PTHREAD define barfed - Did you include e.g. '-pthread' in your compile flags?
  #endif

  #ifdef USE_OMP
	#error Phreads is the only supported multithreading API!
//	ASSERT(HERE, MAX_THREADS = omp_get_num_procs(), "Illegal #Cores value stored in MAX_THREADS");
  #elif(defined(USE_PTHREAD))
	ASSERT(HERE, MAX_THREADS =     get_num_cores(), "Illegal #Cores value stored in MAX_THREADS");
  #else
	#error Unrecognized multithreading model!
  #endif
	// MAX_THREADS based on number of processing cores will most often be a power of 2, but don't assume that.
	ASSERT(HERE, MAX_THREADS > 0,"Mlucas.c: MAX_THREADS must be > 0");

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
#ifdef TEST_FFT_RADIX
  #ifdef USE_SSE2
	#error TEST_FFT_RADIX requires non-SIMD build!
  #endif
	test_fft_radix();
	exit(0);
#endif
}

// Jun 2015: Thanks to Alex Vong for the detective work here - In 32-bit Linux, may need
// to up the stacklimit from the defaults to avoid SIGSEGV faults during running of the
// alloc-heavy self-test suites.
// Since we see such crashes under at least one major distro (Debian),
// invoke the code below for all Linuxes, since at worst it will amount to a no-op.
// OS X shows no issues, but since the implementation here uses only the Posix-compliant
// [get,set]rlimit() utilities (as opposed to the Linux-only prlimit()), could easily extend
// to OS X or, better, "All Posix-compliant OSes", if needed, by changing defined(OS_TYPE_LINUX)
// in the #if function-body wrapper to defined(OS_POSIX_COMPLIANT).

/* If needed, set stack soft limit to stack hard limit and resume execution after registering
the changed setting with the shell, using execvp().

From http://linux.die.net/man/2/setrlimit, here is form of struct returned by getrlimit():
	struct rlimit {
		rlim_t rlim_cur;  // Soft limit
		rlim_t rlim_max;  // Hard limit (ceiling for rlim_cur)
	};
(I use C++ // inside the struct, so as to play nice with my own C-style comment wrapper here).
*/
void set_stacklimit_restart(char *argv[])
{
#if defined(OS_TYPE_LINUX) && defined(CPU_IS_X86)	// CPU_IS_X86 (platform.h) == '32-bit x86'.
	struct rlimit stack_limits;

	if (getrlimit(RLIMIT_STACK, &stack_limits)) {
		fprintf(stderr, "Call to getrlimit() failed.\n");
		ASSERT(HERE, 0, "Exiting.");
	}
	printf("Old stack_limits: cur = %zu, max = %zu, [RLIM_INFINITY = %zu]\n",
	       stack_limits.rlim_cur, stack_limits.rlim_max, RLIM_INFINITY);

	if (stack_limits.rlim_cur == stack_limits.rlim_max)
		return;
	stack_limits.rlim_cur = stack_limits.rlim_max;

	if (setrlimit(RLIMIT_STACK, &stack_limits)) {
		fprintf(stderr, "Call to setrlimit() failed.\n");
		ASSERT(HERE, 0, "Exiting.");
	}
	printf("New stack_limits: cur = %zu, max = %zu\n",
	       stack_limits.rlim_cur, stack_limits.rlim_max);

	if(execvp(argv[0], argv)) {
		fprintf(stderr, "Call to execvp() failed.\n");
		ASSERT(HERE, 0, "Exiting.");
	}
#endif /* CPU_IS_X86  */
}

/***The following 3 routines MUST BE CALLED IN THE SAME ORDER AS IN host_init()!***/

void print_host_info(void)
{
#if defined(USE_GPU) && defined(__CUDACC__)

	gpu_config_t gpu_config;
	gpu_info_t ginfo;
	int32 igpu;

	gpu_init(&gpu_config);
	if (gpu_config.num_gpu > 0) {
		printf("Detected %u CUDA-enabled GPU devices.\n", gpu_config.num_gpu);
		for(igpu = 0; igpu < gpu_config.num_gpu; ++igpu) {
			ginfo = gpu_config.gpu_info[igpu];
			printf("GPU #%u: %s v%u.%u\n", igpu, ginfo.name, ginfo.major, ginfo.minor);
			printf("clockRate = %u MHz\n", ginfo.clockRate/1000);
			printf("multiProcessorCount = %u\n", ginfo.multiProcessorCount);
			printf("totalConstMem = %u\n", ginfo.totalConstMem);
			printf("sharedMemPerBlock = %u\n", ginfo.sharedMemPerBlock);
			printf("totalGlobalMem = %u\n", ginfo.totalGlobalMem);
			printf("reg[ister]sPerBlock = %u\n", ginfo.regsPerBlock);
			printf("maxThreadsPerBlock = %u\n", ginfo.maxThreadsPerBlock);
			printf("deviceOverlap = %u\n", ginfo.deviceOverlap);
			printf("concurrentKernels = %u\n", ginfo.concurrentKernels);
			printf("warpSize = %u\n", ginfo.warpSize);
			printf("max_thread_dim[3] = [%u,%u,%u]\n", ginfo.maxThreadsDim[0], ginfo.maxThreadsDim[1], ginfo.maxThreadsDim[2]);
			printf("max_grid_size[3] = [%u,%u,%u]\n", ginfo.maxGridSize[0], ginfo.maxGridSize[1], ginfo.maxGridSize[2]);
		}
	} else {
		printf("ERROR: No CUDA-enabled GPUs found\n");
		exit(-1);
	}

	// Disable default spin-loop-wait-for-GPU:
	cudaSetDeviceFlags(cudaDeviceBlockingSync);

	cudaError_t cudaError = cudaGetLastError();
	if(cudaError != cudaSuccess)
	{
		printf("ERROR: cudaGetLastError() returned %d: %s\n", cudaError, cudaGetErrorString(cudaError));
		ASSERT(HERE, 0, "gpu_sieve: GPU-side error detected!");
	}

//	cudaVecAddTest();
	cudaVecModpowTest64();
	cudaVecModpowTest78_0();
	cudaVecModpowTest78();
	cudaVecModpowTest96();
  #ifdef USE_FMADD
	cudaVecMul50x50Test();
  #endif
//exit(0);

#endif

#if EWM_DEBUG
	printf("INFO: Program compiled with debugging diagnostics ON.\n");
#endif

	printf("CPU Family = %s, OS = %s, %2d-bit Version, compiled with %s, Version %s.\n", CPU_NAME, OS_NAME, OS_BITS, COMPILER_NAME, COMPILER_VERSION);

#ifdef CPU_IS_ARM_EABI

	// Apr 2018: Due to portability issues, replace the system-headers-based version of the "has advanced SIMD?"
	// check with one based on what amounts to "is the result of 'grep asimd /proc/cpuinfo' empty or not?":
  #if 1

	int has_asimd(void)
	{
		char in_line[STR_MAX_LEN];
		FILE*fp = mlucas_fopen("/proc/cpuinfo", "r");
		ASSERT(HERE, fp != 0x0, "/proc/cpuinfo file not found!");
		while(fgets(in_line, STR_MAX_LEN, fp) != 0x0) {
			if(strstr(in_line, "asimd") != 0)
				return 1;
		}
		fclose(fp);	fp = 0x0;
		return 0;
	}

  #elif __ARM_ARCH >= 8 // Rest of the preprocessor-conditional is the old version:

	#error This system-header-based ARM-has-ASIMD code should not be used!
	// Thanks to Laurent Desnogues for this:
	#include <sys/auxv.h>
	#include <asm/hwcap.h>
	// Check for flag indicating CPU supports advanced SIMD instruction set.
	// NB: For reasons unknown, when I tried putting this function def into get_cpuid.c as I do with the
	// x86 has_sse2() and similar functions, everything compiled fine but got linker errors on the ARM,
	// linker was unable to find the has_asimd() function. After dicking around with that problem for
	// several hours (and hence, several hours too many), tried moving the def here, and it worked:
	int has_asimd(void)
	{
		unsigned long hwcaps = getauxval(AT_HWCAP);
	#ifndef HWCAP_ASIMD	// This is not def'd on pre-ASIMD platforms
		const unsigned long HWCAP_ASIMD = 0;
	#endif
		if (hwcaps & HWCAP_ASIMD) {
			return 1;
		}
		return 0;
	}

  #else	// Pre-v8 ARM does not have above asimd-headers, no point even looking for them

	int has_asimd(void) { return 0; }

  #endif

	if(has_asimd()) {
	#ifdef USE_ARM_V8_SIMD
		printf("INFO: Build uses ARMv8 advanced-SIMD instruction set.\n");
	#else
		printf("INFO: CPU supports ARMv8 advanced-SIMD instruction set, but using scalar floating-point build.\n");
	#endif
	} else {
	#ifdef USE_ARM_V8_SIMD
		ASSERT(HERE, 0, "#define USE_ARM_V8_SIMD invoked but no advanced-SIMD support detected on this CPU!\n");
	#endif
	}

#elif(defined(CPU_IS_X86) || defined(CPU_IS_IA64) || defined(CPU_IS_X86_64))

//	get_cpu();

/* Enable this call to get gory details: */
	#if(1)
		/* if(1) --> if(0) Enables section below */
	#elif(defined(COMPILER_TYPE_MSVC) || defined(COMPILER_TYPE_ICC))
		cpu_details();
	#endif

  #ifdef USE_AVX512

	if(has_avx512()) {
		printf("INFO: Build uses AVX512 instruction set.\n");
	} else {
		#define CPUID(arg1,arg2,ax,bx,cx,dx)\
			__asm__ __volatile__ ("cpuid":\
		"=a" (ax), "=b" (bx), "=c" (cx), "=d" (dx) : "a" (arg1), "c" (arg2));
		uint32 a,b,c,d;
		CPUID(1,0,a,b,c,d);
		printf("has_avx512: CPUID returns [a,b,c,d] = [%8X,%8X,%8X,%8X]\n",a,b,c,d);
		printf("#define USE_AVX512 invoked but no FMA support detected on this CPU! Check get_cpuid functionality and CPU type.\n");
		ASSERT(HERE, 0, "#define USE_AVX512 invoked but no FMA support detected on this CPU! Check get_cpuid functionality and CPU type.\n");
	}

  #elif(defined(USE_AVX2))

	if(has_avx2()) {
		printf("INFO: Build uses AVX2 instruction set.\n");
	} else {
		#define CPUID(arg1,arg2,ax,bx,cx,dx)\
			__asm__ __volatile__ ("cpuid":\
		"=a" (ax), "=b" (bx), "=c" (cx), "=d" (dx) : "a" (arg1), "c" (arg2));
		uint32 a,b,c,d;
		CPUID(1,0,a,b,c,d);
		printf("has_avx2: CPUID returns [a,b,c,d] = [%8X,%8X,%8X,%8X]\n",a,b,c,d);
		printf("#define USE_AVX2 invoked but no FMA support detected on this CPU! Check get_cpuid functionality and CPU type.\n");
		ASSERT(HERE, 0, "#define USE_AVX2 invoked but no FMA support detected on this CPU! Check get_cpuid functionality and CPU type.\n");
	}

  #elif(defined(USE_AVX))

	if(has_avx2()) {
		printf("INFO: CPU supports AVX2 instruction set, but using AVX-enabled build.\n");
	} else if(has_avx()) {
		printf("INFO: Build uses AVX instruction set.\n");
	} else {
		ASSERT(HERE, 0, "#define USE_AVX invoked but no AVX support detected on this CPU! Check get_cpuid functionality and CPU type.\n");
	}

  #elif(defined(USE_SSE2))
	/* This doesn't work on non-AVX platforms, since XGETBV (needed by has_avx*() functions) does not exist
	if(has_avx2()) {
		printf("INFO: CPU supports AVX2 instruction set, but using SSE2-enabled build.\n");
	} else if(has_avx()) {
		printf("INFO: CPU supports AVX instruction set, but using SSE2-enabled build.\n");
	} else */
	if(has_sse2()) {
		printf("INFO: Build uses SSE2 instruction set.\n");
	} else {
		ASSERT(HERE, 0, "#define USE_SSE2 invoked but no SSE2 support detected on this CPU! Check get_cpuid functionality and CPU type.\n");
	}

  #else

	if(has_sse2()) {
		printf("INFO: CPU supports SSE2 instruction set, but using scalar floating-point build.\n");
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

	set_mlucas_path();
	printf("INFO: MLUCAS_PATH is set to \"%s\"\n", MLUCAS_PATH);

	if (MLUCAS_PATH[0] != '\0') {
		if (!mkdir_p(MLUCAS_PATH)) {
			printf("INFO: mkdir -p \"%s\" succeeded\n", MLUCAS_PATH);
		} else {
			fprintf(stderr, "FATAL: mkdir -p \"%s\" failed\n", MLUCAS_PATH);
			ASSERT(HERE, 0, "Exiting.");
		}
	}
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

/* Re.the SIMD control word, from http://softpixel.com/~cwright/programming/simd/sse.php :

	The MXCSR register is a 32-bit register containing flags for control and status information regarding SSE instructions.
	As of SSE3, only bits 0-15 have been defined.
														[EWM: Default value on MSVC/ia32 = 0x1FA0, so the bits marked with [x] are set:]
	Mnemonic	Bit Location	Description				[EWM: Default value on GCC/Core2-ia64 = 0x1FA2 = 1111110100010, bits [y] are set:]
	--------	------------	---------------------	[EWM: Default value on GCC/Haswell-ia64 = 0x9FE2 = 1111110100010, bits [z] are set:]
		FZ		bit 15			Flush To Zero							[z]
		R+		bit<14:13> = 10	Round Positive
		R-		bit<14:13> = 01	Round Negative
		RZ		bit<14:13> = 11	Round To Zero
		RN		bit<14:13> = 00	Round To Nearest		[x]		[y]		[z]
		PM		bit 12			Precision Mask			[x]		[y]		[z]
		UM		bit 11			Underflow Mask			[x]		[y]		[z]
		OM		bit 10			Overflow Mask			[x]		[y]		[z]
		ZM		bit 9			Divide By Zero Mask		[x]		[y]		[z]
		DM		bit 8			Denormal Mask			[x]		[y]		[z]
		IM		bit 7			Invalid Operation Mask	[x]		[y]		[z]
		DAZ		bit 6			Denormals Are Zero						[z]
		PE		bit 5			Precision Flag			[x]		[y]		[z]
		UE		bit 4			Underflow Flag
		OE		bit 3			Overflow Flag
		ZE		bit 2			Divide By Zero Flag
		DE		bit 1			Denormal Flag					[y]		[z]
		IE		bit 0			Invalid Operation Flag
*/
// May 2018: changed things so that the function is *defined* in both 32 and 64-bit modes:
#if defined(CPU_IS_X86) || defined(CPU_IS_X86_64)

	/* Example: To flip bits 13:14 in MXCSR from their default value 00 (round-to-nearest) 11 (round-toward-0):
		uint32 i = 0x00006000; i = x86_simd_mxcsr_toggle(i);
	*/
	// Return current value of the MXCSR
	uint32 x86_simd_mxcsr_getval(void) {
	#ifdef USE_SSE2	// Only supported if SSE is
		uint32 MXCSR_VALUE;
		__asm__ volatile ("stmxcsr %0" : "=m" (MXCSR_VALUE) );
		 return MXCSR_VALUE;
	#else
		return 0;
	#endif
	}
	// Set value of the MXCSR to the specified one. Returns the old value, to support reset-to-default:
	uint32 x86_simd_mxcsr_setval(uint32 MXCSR_NEWVAL) {
	#ifdef USE_SSE2	// Only supported if SSE is
		uint32 MXCSR_OLDVAL;
		__asm__ volatile ("stmxcsr %0" : "=m" (MXCSR_OLDVAL) );
		__asm__ volatile ("ldmxcsr %0" :: "m" (MXCSR_NEWVAL) );
		return MXCSR_OLDVAL;
	#else
		return 0;
	#endif
	}
	// For every set bit in the input XOR-mask, flip the corresponding bit in the MXCSR. Returns the old value.
	uint32 x86_simd_mxcsr_toggle(uint32 MXCSR_MASK) {
	#ifdef USE_SSE2	// Only supported if SSE is
		uint32 MXCSR_OLDVAL,MXCSR_NEWVAL;
		__asm__ volatile ("stmxcsr %0" : "=m" (MXCSR_OLDVAL) );
		MXCSR_NEWVAL = MXCSR_OLDVAL ^ MXCSR_MASK;
		__asm__ volatile ("ldmxcsr %0" :: "m" (MXCSR_NEWVAL) );
		return MXCSR_OLDVAL;
	#else
		return 0;
	#endif
	}

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
		ASSERT(HERE, (FPU_MODE == FPU_64RND) || (FPU_MODE == FPU_64CHOP), "Illegal value of FPU_MODE");

		// Check the SIMD control word:
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
		#warning INFO: Setting rnd_const-emulated DNINT for 64-bit x86 register-double significand
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

/******* DEFINE GLOBALS, CHECK TYPE-LENGTHS AND ENDIANNESS, AND TEST RND-CONST FAST-NINT: *******/
#define FAST_UINT32_MOD	0	// Set = 1 to enable this test

void check_nbits_in_types(void)
{
	uint32 i,j;
	double ftmp, fran, ferr, finv, fsrt;
	double tpi = 3.1415926535897932384;
	double ln2 = LOG2;

	/* Make sure TRUE and FALSE behave as required: */
	ASSERT(HERE, !FALSE && TRUE, "TRUE and FALSE do not behave as required in check_nbits_in_types");

	/* Check lengths of basic data types: */
    ASSERT(HERE, sizeof( int8 ) == 1, "sizeof( int8 ) != 1");
    ASSERT(HERE, sizeof(uint8 ) == 1, "sizeof(uint8 ) != 1");
    ASSERT(HERE, sizeof( int16) == 2, "sizeof( int16) != 2");
    ASSERT(HERE, sizeof(uint16) == 2, "sizeof(uint16) != 2");
    ASSERT(HERE, sizeof( int32) == 4, "sizeof( int32) != 4");
    ASSERT(HERE, sizeof(uint32) == 4, "sizeof(uint32) != 4");
    ASSERT(HERE, sizeof( int64) == 8, "sizeof( int64) != 8");
    ASSERT(HERE, sizeof(uint64) == 8, "sizeof(uint64) != 8");
    ASSERT(HERE, sizeof(uint64) >= sizeof(void*), "sizeof(long long) != sizeof(void*)");    /* ALIGN_DOUBLES assumes this. */

	/* AltiVec vector types: */
#if(CPU_HAS_ALTIVEC || CPU_IS_CELL)
	ASSERT(HERE, sizeof(vec_uint8X16) == 16 , "sizeof(vec_uint8X16) != 16 ");
	ASSERT(HERE, sizeof(vec_uint16X8) == 16 , "sizeof(vec_uint16x8) != 16 ");
	ASSERT(HERE, sizeof(vec_uint32X4) == 16 , "sizeof(vec_uint32x4) != 16 ");
#endif

	uint64 x = 0x0706050403020100ull;
	uint8 *byte_arr = (uint8*)&x;
	// Runtime ordering is little-endian:
	if(byte_arr[0] == 0 && byte_arr[1] == 1 && byte_arr[2] == 2 && byte_arr[3] == 3 && byte_arr[4] == 4 && byte_arr[5] == 5 && byte_arr[6] == 6 && byte_arr[7] == 7) {
	  #ifdef USE_BIG_ENDIAN
		ASSERT(HERE, 0, "USE_BIG_ENDIAN set in platform.h but little-endian detected at runtime!");
	  #endif
	} else if(byte_arr[0] == 7 && byte_arr[1] == 6 && byte_arr[2] == 5 && byte_arr[3] == 4 && byte_arr[4] == 3 && byte_arr[5] == 2 && byte_arr[6] == 1 && byte_arr[7] == 0) {
	  #ifndef USE_BIG_ENDIAN
		ASSERT(HERE, 0, "USE_BIG_ENDIAN not set in platform.h but big-endian detected at runtime!");
	  #endif
	} else {
		ASSERT(HERE, 0, "Endianness detected as neither big nor little-endian at runtime!");
	}

	// Init RNG:
	rng_isaac_init(TRUE);

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

#if 0
#error Code obsolete as of Dec 2015!
	/* We typically need more information re. the FFT-mul params before being
	able to inteligently set the anti-thrashing array-padding params, so set = -1
	(which is taken to mean uninitialized) here:
	*/
	DAT_BITS = PAD_BITS = (int32)0xffffffff;
#else
	// Dec 20155: Subquadratic GCD needs FFT-mul with variable array lengths, thus no longer
	// have the "one FFT length at a time" simplicity of LL-testing where we can fiddle these
	// depending on exponent being tested, so now just set FFT array-padding params right here:
	DAT_BITS = DAT_BITS_DEF;	PAD_BITS = PAD_BITS_DEF;
	printf("Setting DAT_BITS = %d, PAD_BITS = %d\n",DAT_BITS,PAD_BITS);
#endif

	FFT_MUL_BASE = (double)((uint64)1 << FFT_MUL_BITS);
/* Intend to relax this later to allow powers of 2 as large as 2^54: */
ASSERT(HERE, ((uint64)FFT_MUL_BASE >> 16) == 1, "util.c: FFT_MUL_BASE != 2^16");

	ASSERT(HERE, trailz64((uint64)FFT_MUL_BASE) == FFT_MUL_BITS, "mi64_cvt_double_uint64: trailz64((uint64)FFT_MUL_BASE) != FFT_MUL_BITS");
	ASSERT(HERE, DNINT(FFT_MUL_BASE) == FFT_MUL_BASE, "mi64_cvt_double_uint64: FFT_MUL_BASE not pure-integer!");
	ASSERT(HERE, FFT_MUL_BASE < 1.0*0x8000000*0x8000000, "mi64_cvt_double_uint64: FFT_MUL_BASE >= maximum allowed value of 2^54!");
	FFT_MUL_BASE_INV = 1.0/FFT_MUL_BASE;

  #if FAST_UINT32_MOD
	// Test fast 32-bit mod algo:
	#define gen_pinv(p)	(0xffffffff / (p) )

	const uint32 BITS1 = 24, mask1 = 0xffffffff >> (32 - BITS1);	// mask for x-input
	const uint32 BITS2 = 17, mask2 = 0xffffffff >> (32 - BITS2);	// mask for p-input
	uint32 x,p,pinv, r,y, nfail = 0;

	/* Using 10^9 random (x,p) pairs shows no failures up to the full 32 bits, when we fail in cases where pinv = 1.
	Here is a small sampling of the resulting #fails:

		I = 534397: Incorrect a (mod p) result! a = 1698375046, p = 3399405424, pinv = 1: r = 2593936918 [exact = 1698375046]
		I = 534400: Incorrect a (mod p) result! a = 471952975, p = 3563494799, pinv = 1: r = 1203425472 [exact = 471952975]
		I = 534401: Incorrect a (mod p) result! a = 1844700405, p = 3680268453, pinv = 1: r = 2459399248 [exact = 1844700405]
		I = 534409: Incorrect a (mod p) result! a = 190672429, p = 2680969504, pinv = 1: r = 1804670221 [exact = 190672429]
		I = 534410: Incorrect a (mod p) result! a = 2507724634, p = 4081983633, pinv = 1: r = 2720708297 [exact = 2507724634]
		I = 534411: Incorrect a (mod p) result! a = 1772269933, p = 3052535324, pinv = 1: r = 3014701905 [exact = 1772269933]
		I = 534418: Incorrect a (mod p) result! a = 1791901102, p = 4244017845, pinv = 1: r = 1842850553 [exact = 1791901102]

	In each of these cases the correct result is simply the unmodified x-input value, but our algo gives MUL-pair
	result r = 0, thus
		x = x - 0;
		r = x - p;
	and since p > x, the subtract has a borrow, i.e. yields 2^32 + x - p, which is only correct if interpreted as signed.
	To guard against this we simply add logic that checks if the MUL-pair result is 0 and returns x if true. That solves
	part of our problem, but we still must compute r = x - p, now guarding against underflow-on-subtract, as shown below.

	For inputs < 2^31, though, we can just use the simpler 4-step algorithm sans any conditionals except for the final
		r = (r >= p) ? x : r ;
	select step.
	*/
	uint32 ntry = 10000000, nneg = 0;
	printf ("INFO: Performing Fast-uint32-mod test with %u [%u,%u]-bit input pairs ... ",ntry,BITS1,BITS2);
	for(i = 0; i < ntry; i++)	// Cut #test down to 10^7 for production version
	{
		uint64 i64 = rng_isaac_rand();
		x = (uint32)(i64 >> 32);
		p = (uint32)i64;
		x &= mask1;	// Only support inputs of 24-bits or less
		p &= mask2;
		// Must guard against 0 modulus:
		while(p == 0) {
			i64 = rng_isaac_rand();
			p = (uint32)i64 & mask2;
		//	printf ("I = %d: p = 0! Replacing with %u\n",i, p);
		}
		y = x;	// Use x-copy in place of x below, since must leave inputs as-is
		pinv = gen_pinv(p);
		// Compute x % p (we hope) and check:
		r = __MULH32(y, pinv);	//	r = __mulhi (y, pinv);
		r = __MULL32(r, p   );	//	r = r * p;
		if(r != 0) {
			y = y - r;
			r = y - p;
			nneg += (r < p);	// Collect stats on how many of these cases - i.e. where y did need downward adjustment - occur
		//	if(r < p)printf ("I = %d Needed extra sub: a = %u; p = %u; pinv = %u [a/p = %f]: y = %u, r = %u]\n",i, x, p, pinv,(float)x/p, y,r);
			r = (r >= p) ? y : r ;
		} else {
			r = y - p;
			r = (r > y) ? y : r ;	// Must guard vs underflow-on-subtract
		}
		if (r != x%p) {
			++nfail;
			printf ("I = %d: Incorrect a (mod p) result! a = %u, p = %u, pinv = %u: r = %u [exact = %u]\n",i, x, p, pinv, r, x%p);
		}
	}
	printf ("%u cases of %u [%6.2f%%] needed adjustment.\n",nneg,ntry,100.*nneg/(float)ntry);
	ASSERT(HERE, nfail == 0, "Fast-uint32-mod test failed for 1 or more inputs!");
  #endif	// #if FAST_UINT32_MOD ?

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
	for(i = 0; i < 100000; i++)
	{
		fran = rng_isaac_rand_double();
		fran = fabs(fran);
		if(fran > 0.0) {
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
		if(fran < 0.0 || fran >= 1.0) {
			sprintf(cbuf, "check_nbits_in_types: rng_isaac_rand_double_norm_pos returns illegal value outside [0, 1): i = %d, %e\n", i,fran);
			ASSERT(HERE, 0, cbuf);
		}

		fran = rng_isaac_rand_double_norm_pm1();
		if(fabs(fran) >= 1.0) {
			sprintf(cbuf, "check_nbits_in_types: rng_isaac_rand_double_norm_pm1 returns illegal value outside (-1,+1): i = %d, %e\n", i, fran);
			ASSERT(HERE, 0, cbuf);
		}
	}

#ifdef USE_FMADD
#warning USE_FMADD enabled in util.c!
	/* Test any FMADD-based routines, if def'd: */
	printf("INFO: Testing MUL50x50 routines ... ");
	uint32 nerr = test_mul50x50();
	if(!nerr)
		printf("fma_dmult_tests completed successfully!\n");
	else
		ASSERT(HERE, 0, "fma_dmult_tests failed!\n");
#endif

#ifdef USE_FGT61
	printf("INFO: Testing FGT (mod M61) arithmetic ... \n");
	const uint64 q = 0x1FFFFFFFFFFFFFFFull;
	uint64 order, root_re,root_im, re,im;

	// Test various functions defined in fgt_m61.c - modular powering (both that
	// used in root-taking and of the results is a good 'whole code' test.

	// [1] Test out power-of-2 roots in GF(q^2) -
	//	Maximal order (q^2-1) = 2^62 * (2^60-1), allowing power-of-2 roots up to 2^62:
	/* Here are the roots found in the loop below:
		FGT: prim-root of order 2^ 0 = 1 + I*0
		FGT: prim-root of order 2^ 1 = 2305843009213693950 + I*0 == -1 (mod q)
		FGT: prim-root of order 2^ 2 = 0 + I*1
		FGT: prim-root of order 2^ 3 =          1073741824 + I*         1073741824 = 2^30 * (1 + I)
		FGT: prim-root of order 2^ 4 = 1693317751237720973 + I*2283815672160731785
		FGT: prim-root of order 2^ 5 =  697323983679957246 + I*  83304533336094567
		...
		FGT: prim-root of order 2^62 = 1895584235299698857 + I*2150133374943417338
	*/
	for(i = 0; i < 63; i++)
	{
		order = 1ull << i;
		prim_root_q(order, &root_re,&root_im);
	//	printf("FGT: prim-root of order 2^%2u = %llu + I*%llu\n",i, root_re,root_im);
		// Check order-primitivity of roots of order > 1 by powering result up to 2nd order; result must == -1 (mod q):
		if(i > 0) {
			for(j = 1; j < i; j++) {
				csqr_modq(root_re,root_im, &root_re,&root_im);
				root_re = qreduce(root_re);	root_im = qreduce(root_im);	// Only partially reduce intermediates...
			}
			root_re = qreduce_finish(root_re);	root_im = qreduce_finish(root_im);	// ...and then finish reducing here.
			ASSERT(HERE, root_re ==  q-1 && root_im == 0ull, "Bad prim_root_q result!");
		} else {
			ASSERT(HERE, root_re == 1ull && root_im == 0ull, "Bad prim_root_q result!");
		}
	}

#if 0
	// Play with conjugates of both power-of-2 and non-power-of-2 (but still even-order) roots:
	// Power-of-2 roots satisfy simple conjugate rule, modular analog of complex conj(Re,Im) = (Re,-Im):
	order = 16;	prim_root_q(order, &root_re,&root_im);
	pow_modq(order-1, root_re,root_im, &re,&im);
	printf("FGT: prim-root of order %u = %llu + I*%llu, Conjugate = %llu + I*%llu [q-Im = %llu]\n",(uint32)order, root_re,root_im, re,im,q-im);
//	FGT: prim-root of order 16 = 1693317751237720973 + I*2283815672160731785,
//					Conjugate =  1693317751237720973 + I*  22027337052962166 [q-Im = 2283815672160731785]
	ASSERT(HERE, root_re == re && root_im == (q-im), "Bad power-of-2 conjugate!");

	// Non-power-of-2 roots satisfy no simple conjugate rules, so multiply root and its conjugate together as sanity check:
	order = 24;	prim_root_q(order, &root_re,&root_im);
	pow_modq(order-1, root_re,root_im, &re,&im);
	printf("FGT: prim-root of order %u = %llu + I*%llu, Conjugate = %llu + I*%llu [q-Im = %llu]\n",(uint32)order, root_re,root_im, re,im,q-im);
	cmul_modq(root_re,root_im, re,im, &re,&im);
	re = qreduce_full(re);	im = qreduce_full(im);
	ASSERT(HERE, re == 1ull && im == 0ull, "Bad non-power-of-2 conjugate!");
/*
	24th root:
	FGT: prim-root of order 24 = 244692701471512749 + I*2061150307742181202,
					Conjugate = 2061150308815923026 + I*2061150308815923026 [q-Im = 244692700397770925]

	*NOTE* Conjugate has Re == Im! Thus reminiscent of 8th root:
	FGT: prim-root of order 8 =          1073741824 + I*         1073741824 = 2^30 * (1 + I)

	FGT: prim-root of order 3 = 1669582390241348315 + I*0 = R3

	printf("Powers of prim-root:\n");
	re = root_re;	im = root_im;
	for(i = 0; i < order; i++) {
		printf("%2u: %20llu[-= %20llu] + I*%20llu[-= %20llu]\n",i+1, re,q-re,im,q-im);
		cmul_modq(root_re,root_im, re,im, &re,&im);
		re = qreduce_full(re);	im = qreduce_full(im);
	}
Gives
	Powers of prim-root:
	 1:   244692701471512749[-=  2061150307742181202] + I* 2061150307742181202[-=   244692701471512749]	<*** Call thib [a,b]
	 2:                    0[-=  2305843009213693951] + I*  636260618972345636[-=  1669582390241348315]	[0,c]
	 3:           1073741824[-=  2305843008139952127] + I*          1073741824[-=  2305843008139952127]	[d,d]
	 4:  1669582390241348316[-=   636260618972345635] + I*                   0[-=  2305843009213693951]	[e,0]
	 5:   244692700397770925[-=  2061150308815923026] + I* 2061150308815923026[-=   244692700397770925]	[f,-f]
	 6:                    0[-=  2305843009213693951] + I*                   1[-=  2305843009213693950]	[0,1]
	 7:   244692701471512749[-=  2061150307742181202] + I*  244692701471512749[-=  2061150307742181202]	[a,-b]
	 8:  1669582390241348315[-=   636260618972345636] + I*                   0[-=  2305843009213693951]	[-c,0]
	 9:  2305843008139952127[-=           1073741824] + I*          1073741824[-=  2305843008139952127]	[-d,d]
	10:                    0[-=  2305843009213693951] + I* 1669582390241348316[-=   636260618972345635]	[0,e]
	11:   244692700397770925[-=  2061150308815923026] + I*  244692700397770925[-=  2061150308815923026]	[f,f]
	12:  2305843009213693950[-=                    1] + I*                   0[-=  2305843009213693951]	[-1,0]
	13:  2061150307742181202[-=   244692701471512749] + I*  244692701471512749[-=  2061150307742181202]	[-a,b]
	14:                    0[-=  2305843009213693951] + I* 1669582390241348315[-=   636260618972345636]	[0,-c]
	15:  2305843008139952127[-=           1073741824] + I* 2305843008139952127[-=           1073741824]	[-d,-d]
	16:   636260618972345635[-=  1669582390241348316] + I*                   0[-=  2305843009213693951]	[-e,0]
	17:  2061150308815923026[-=   244692700397770925] + I*  244692700397770925[-=  2061150308815923026]	[-f,f]
	18:                    0[-=  2305843009213693951] + I* 2305843009213693950[-=                    1]	[0,-1]
	19:  2061150307742181202[-=   244692701471512749] + I* 2061150307742181202[-=   244692701471512749]	[-a,-b]
	20:   636260618972345636[-=  1669582390241348315] + I*                   0[-=  2305843009213693951]	[c,0]
	21:           1073741824[-=  2305843008139952127] + I* 2305843008139952127[-=           1073741824]	[d,-d]
	22:                    0[-=  2305843009213693951] + I*  636260618972345635[-=  1669582390241348316]	[0,-e]
	23:  2061150308815923026[-=   244692700397770925] + I* 2061150308815923026[-=   244692700397770925]	[-f,-f]
	24:                    1[-=  2305843009213693950] + I*                   0[-=  2305843009213693951]	[1,0]
The four [+-d,+-d] and four powers of I are just the eight 8th roots of unity which are hit on multiple-of-3 index values.
*/
#endif

	// [2] Test out odd-order roots - that means any order dividing (2^60-1) = 3^2.5^2.7.11.13.31.41.61.151.331.1321:
	/* Here are the roots found in the loop below:
		FGT: prim-root of order 3 = 1669582390241348315 + I*0
		FGT: prim-root of order 9 = 1102844585000305877 + I*0
		FGT: prim-root of order 45 = 295230898440480023 + I*0
		...
		FGT: prim-root of order 1152921504606846975 = 1754029865706415802 + I*0
	*/
	const uint32 odd_ord_facs[13] = {3,3,5,5,7,11,13,31,41,61,151,331,1321};
	order = 1ull;
	for(i = 0; i < 13; i++)
	{
		order *= odd_ord_facs[i];
		prim_root_q(order, &root_re,&root_im);
	//	printf("FGT: prim-root of order %llu = %llu + I*%llu\n",order, root_re,root_im);
		ASSERT(HERE, root_im == 0ull, "Odd roots must be strictly real!!");
		// Check order-primitivity of roots by raising result to (order)th power; result must == -1 (mod q):
		pow_modq(order, root_re,root_im, &root_re,&root_im);
		ASSERT(HERE, root_re == 1ull && root_im == 0ull, "Bad prim_root_q result!");
	}
	printf("fgt_m61 tests completed successfully!\n");
#endif

	return;
}

/*
Fast-uint32-mod test:
Any obvious pattern in the need-adjustment cases? Fractional part of a/p shows it's cases where a/p has small frac-part:
I = 19 Needed extra sub: a = 1024576945; p = 12037616; pinv = 356 [a/p = 85.114609]: y = 13417201, r = 1379585]
I = 25 Needed extra sub: a = 230529187; p = 11523083; pinv = 372 [a/p = 20.005859]: y = 11590610, r = 67527]
I = 32 Needed extra sub: a = 211760194; p = 111687; pinv = 38455 [a/p = 1896.014648]: y = 113329, r = 1642]
I = 43 Needed extra sub: a = 1066942432; p = 6128196; pinv = 700 [a/p = 174.103836]: y = 6764524, r = 636328]
I = 50 Needed extra sub: a = 968649863; p = 8139751; pinv = 527 [a/p = 119.002396]: y = 8159245, r = 19494]
I = 87 Needed extra sub: a = 968400190; p = 10515906; pinv = 408 [a/p = 92.089088]: y = 11452744, r = 936838]
I = 94 Needed extra sub: a = 1057486414; p = 15537360; pinv = 276 [a/p = 68.060883]: y = 16483294, r = 945934]
I = 109 Needed extra sub: a = 533959987; p = 14415015; pinv = 297 [a/p = 37.041931]: y = 15019447, r = 604432]
I = 111 Needed extra sub: a = 993195985; p = 3571948; pinv = 1202 [a/p = 278.054443]: y = 3766389, r = 194441]
I = 130 Needed extra sub: a = 992040856; p = 6158255; pinv = 697 [a/p = 161.091217]: y = 6720056, r = 561801]
I = 173 Needed extra sub: a = 770797430; p = 10005071; pinv = 429 [a/p = 77.040680]: y = 10412034, r = 406963]
I = 182 Needed extra sub: a = 811377653; p = 2526617; pinv = 1699 [a/p = 321.132050]: y = 2860213, r = 333596]
I = 188 Needed extra sub: a = 1013613270; p = 8803835; pinv = 487 [a/p = 115.133148]: y = 9976080, r = 1172245]
I = 204 Needed extra sub: a = 862933752; p = 4663527; pinv = 920 [a/p = 185.038864]: y = 4844784, r = 181257]
I = 209 Needed extra sub: a = 381286368; p = 4536865; pinv = 946 [a/p = 84.041817]: y = 4726573, r = 189708]
I = 229 Needed extra sub: a = 311393414; p = 12943483; pinv = 331 [a/p = 24.057930]: y = 13693305, r = 749822]
I = 241 Needed extra sub: a = 875437669; p = 6339203; pinv = 677 [a/p = 138.099014]: y = 6966858, r = 627655]
I = 279 Needed extra sub: a = 860144752; p = 8108626; pinv = 529 [a/p = 106.077744]: y = 8739022, r = 630396]
I = 282 Needed extra sub: a = 874645009; p = 14087847; pinv = 304 [a/p = 62.085072]: y = 15286342, r = 1198495]
I = 295 Needed extra sub: a = 131723624; p = 10972820; pinv = 391 [a/p = 12.004537]: y = 11022604, r = 49784]
I = 312 Needed extra sub: a = 478423934; p = 1696304; pinv = 2531 [a/p = 282.039032]: y = 1762510, r = 66206]
I = 346 Needed extra sub: a = 773372954; p = 7359080; pinv = 583 [a/p = 105.090981]: y = 8028634, r = 669554]
I = 347 Needed extra sub: a = 769693293; p = 3115804; pinv = 1378 [a/p = 247.028793]: y = 3205509, r = 89705]
I = 354 Needed extra sub: a = 753736067; p = 7317291; pinv = 586 [a/p = 103.007530]: y = 7372385, r = 55094]
I = 355 Needed extra sub: a = 402434763; p = 3193892; pinv = 1344 [a/p = 126.001366]: y = 3198263, r = 4371]
I = 373 Needed extra sub: a = 303380021; p = 12635991; pinv = 339 [a/p = 24.009199]: y = 12752228, r = 116237]
I = 375 Needed extra sub: a = 980850028; p = 11514899; pinv = 372 [a/p = 85.180954]: y = 13598512, r = 2083613]
I = 385 Needed extra sub: a = 484277536; p = 4000805; pinv = 1073 [a/p = 121.045021]: y = 4180936, r = 180131]
I = 400 Needed extra sub: a = 453815187; p = 6670148; pinv = 643 [a/p = 68.036751]: y = 6915271, r = 245123]
I = 404 Needed extra sub: a = 867743335; p = 6621183; pinv = 648 [a/p = 131.055634]: y = 6989545, r = 368362]
I = 433 Needed extra sub: a = 866814550; p = 12037278; pinv = 356 [a/p = 72.010841]: y = 12167812, r = 130534]
I = 452 Needed extra sub: a = 914104228; p = 6818552; pinv = 629 [a/p = 134.061340]: y = 7236812, r = 418260]
I = 458 Needed extra sub: a = 392296762; p = 14003966; pinv = 306 [a/p = 28.013262]: y = 14189680, r = 185714]
I = 465 Needed extra sub: a = 504882326; p = 10740929; pinv = 399 [a/p = 47.005463]: y = 10799592, r = 58663]
I = 476 Needed extra sub: a = 255805064; p = 2013610; pinv = 2132 [a/p = 127.038033]: y = 2090204, r = 76594]
I = 540 Needed extra sub: a = 699026720; p = 544411; pinv = 7889 [a/p = 1284.005493]: y = 547407, r = 2996]
I = 549 Needed extra sub: a = 1066225928; p = 8260936; pinv = 519 [a/p = 129.068420]: y = 8826120, r = 565184]
I = 561 Needed extra sub: a = 811854178; p = 4584781; pinv = 936 [a/p = 177.075897]: y = 4932722, r = 347941]
I = 571 Needed extra sub: a = 505541881; p = 11747289; pinv = 365 [a/p = 43.034771]: y = 12155743, r = 408454]
I = 590 Needed extra sub: a = 675714616; p = 6756887; pinv = 635 [a/p = 100.003838]: y = 6782803, r = 25916]
I = 605 Needed extra sub: a = 718292932; p = 10559748; pinv = 406 [a/p = 68.021790]: y = 10789816, r = 230068]
I = 616 Needed extra sub: a = 470903290; p = 7846004; pinv = 547 [a/p = 60.018234]: y = 7989054, r = 143050]
I = 639 Needed extra sub: a = 900990749; p = 13056614; pinv = 328 [a/p = 69.006462]: y = 13140997, r = 84383]
I = 654 Needed extra sub: a = 518909096; p = 7518515; pinv = 571 [a/p = 69.017494]: y = 7650076, r = 131561]
I = 658 Needed extra sub: a = 746682360; p = 9439608; pinv = 454 [a/p = 79.100990]: y = 10392936, r = 953328]
I = 693 Needed extra sub: a = 891783098; p = 14139951; pinv = 303 [a/p = 63.068329]: y = 15106136, r = 966185]
I = 731 Needed extra sub: a = 618218131; p = 13436506; pinv = 319 [a/p = 46.010334]: y = 13575361, r = 138855]
I = 732 Needed extra sub: a = 497047750; p = 4518412; pinv = 950 [a/p = 110.004959]: y = 4540842, r = 22430]
I = 740 Needed extra sub: a = 508301929; p = 10796855; pinv = 397 [a/p = 47.078701]: y = 11646599, r = 849744]
I = 857 Needed extra sub: a = 864091071; p = 8538861; pinv = 502 [a/p = 101.195122]: y = 10204971, r = 1666110]
I = 858 Needed extra sub: a = 131898711; p = 13183606; pinv = 325 [a/p = 10.004752]: y = 13246257, r = 62651]
I = 861 Needed extra sub: a = 887750200; p = 4247290; pinv = 1011 [a/p = 209.015686]: y = 4313880, r = 66590]
I = 864 Needed extra sub: a = 435680256; p = 350783; pinv = 12243 [a/p = 1242.022095]: y = 358553, r = 7770]
I = 891 Needed extra sub: a = 860496863; p = 10227272; pinv = 419 [a/p = 84.137474]: y = 11633287, r = 1406015]
I = 897 Needed extra sub: a = 940041261; p = 10427428; pinv = 411 [a/p = 90.150826]: y = 12000169, r = 1572741]
I = 913 Needed extra sub: a = 1059585665; p = 13749374; pinv = 312 [a/p = 77.064285]: y = 14633241, r = 883867]
I = 929 Needed extra sub: a = 804154824; p = 5824066; pinv = 737 [a/p = 138.074463]: y = 6257782, r = 433716]
I = 964 Needed extra sub: a = 1022970393; p = 5911156; pinv = 726 [a/p = 173.057587]: y = 6251561, r = 340405]
I = 981 Needed extra sub: a = 916753724; p = 11581569; pinv = 370 [a/p = 79.156265]: y = 13391342, r = 1809773]
*/

/********************/

#ifdef USE_FMADD

	// Self-test of fmadd-based 50x50==>100-bit exact integer product algorithms:
	uint32	test_mul50x50()
	{
		int i,j,k, pow2;
		double pow2_dmult,pow2_imult;
		uint32 nerr = 0, itmp32;
	#ifdef USE_AVX2
		const double crnd = 3.0*0x4000000*0x2000000, crnd50 = crnd*TWO50FLOAT;	// Consts used to emulate DNINT(x) and 2^50 * DNINT(x*2^-50)
										// (i.e. round-to-nearest-multiple-of-2^50 ... alas the AVX-512 VRNDSCALEPD instruction only supports
										// round-to-nearest-multiple-of-negative-power-of-2, and said power is further restricted to pow < 16.
		static vec_dbl *sc_arr = 0x0;
		static double *sc_ptr;
		double *tmp, *dptr1,*dptr2,*dptr3,*dptr4, l2lo,l2hi, dblo,dbhi, sqr100lo[4],sqr100hi[4], dtmp,cy_max;
		static double *ax,*bx,*cx,*dx, *ay,*by,*cy,*dy, *alo,*blo,*clo,*dlo, *ahi,*bhi,*chi,*dhi, *acy,*alo_norm,*ahi_norm;
		uint64 itmp64, iax,ibx,icx,idx, iay,iby,icy,idy, ialo,iblo,iclo,idlo, iahi,ibhi,ichi,idhi;
		const double prod1_adj = 3.0;	// Const to multiply by base and add to prod[1] to ensure latter >= 0
		if(!sc_arr) {
			sc_arr = ALLOC_VEC_DBL(sc_arr, 8);
			if(!sc_arr){ sprintf(cbuf, "FATAL: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
			sc_ptr = (double *)ALIGN_VEC_DBL(sc_arr);
			ASSERT(HERE, ((uint32)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");
			/* Remember, rhese are POINTERS-TO-DOUBLES, so need an increment of 4 to span an AVX register: */
			tmp = (double *)sc_ptr;
			ax  = tmp + 0;	bx  = tmp + 1;	cx  = tmp + 2;	dx  = tmp + 3;	tmp += 4;
			ay  = tmp + 0;	by  = tmp + 1;	cy  = tmp + 2;	dy  = tmp + 3;	tmp += 4;
			alo = tmp + 0;	blo = tmp + 1;	clo = tmp + 2;	dlo = tmp + 3;	tmp += 4;
			ahi = tmp + 0;	bhi = tmp + 1;	chi = tmp + 2;	dhi = tmp + 3;	tmp += 4;
			acy = tmp + 0;	tmp += 4;
			alo_norm = tmp + 0;	tmp += 4;
			ahi_norm = tmp + 0;	tmp += 4;
		}
	#endif

		// Assumes rng_isaac_init() has already been called on entry
		pow2_dmult = TWO48FLOAT;	// This must match the loop-starting value of 2^pow2:
		pow2_imult = TWO48FLINV;
		dptr1 = &pow2_dmult;	dptr2 = &pow2_imult;
		dptr3 = &crnd50;		dptr4 = &prod1_adj;
		for(pow2 = 48; pow2 < 54; ++pow2) {
			// Only makes sense to test up the #bits in an IEEE-double mantissa: Any larger and we start losing
			// LSBs (I.e. the test may 'succeed' for pow2 > 53, but is only testing the equivalent of pow2 = 53):
			ASSERT(HERE, pow2 < 54, "No point testing > 53-bit inputs due to loss of LSBs!");
			printf("Testing fma_dmult for %d bits, dmult = %f:\n",pow2,pow2_dmult);
			l2lo = l2hi = cy_max = 0.;	// Init log2-range-bounds-storing vars
			for(j = 0; j < 4; j++) {
				sqr100lo[j] = sqr100hi[j] = 0;
			}

			for(i = 0; i < 100000000; i++)
			{
				// Input multiplicands in [-2^pow2, +2^pow2]:
			#ifdef USE_AVX2
				*ax = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	iax = ABS(*ax);	// Will use positive-signed version of
				*bx = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	ibx = ABS(*bx);	// float-double FMA-pair result to ease
				*cx = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	icx = ABS(*cx);	// comparison with exact integer result
				*dx = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	idx = ABS(*dx);	// computed via unsigned MUL_LOHI64.
				*ay = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	iay = ABS(*ay);
				*by = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	iby = ABS(*by);
				*cy = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	icy = ABS(*cy);
				*dy = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	idy = ABS(*dy);

			  #if 0	/********* This FMA-check stuff is costly; disable if using loop for other stuff *********/

				__asm__ volatile (\
					/* Use ymm[< 16] for broadcast-from-scalar consts to avoid illegal-instruction exceptions (AVX-512F supports b'cast-to-full-width-zmm) */\
					"movq	%[__base],%%rax		\n\t	vbroadcastsd	(%%rax),%%ymm15	\n\t"/* BASE */\
					"movq	%[__binv],%%rbx		\n\t	vbroadcastsd	(%%rbx),%%ymm14	\n\t"/* BINV */\
					"movq	%[__ax] ,%%rax	\n\t"\
					"movq	%[__ay] ,%%rbx	\n\t"\
					"vmovaps	(%%rax),%%ymm2	\n\t"\
					"vmovaps	(%%rbx),%%ymm3	\n\t"\
					"     vmulpd	%%ymm2,%%ymm3,%%ymm1	\n\t"/* hi = fma(a,b, 0  ) - use VMULPD instead of VFMADD231PD, since addend = 0. */
					"vmovaps	%%ymm1,%%ymm0				\n\t"/* cpy hi into lo-reg */\
					"vfmsub231pd	%%ymm2,%%ymm3,%%ymm0	\n\t"/* lo = fma(a,b, -hi) */
					"movq	%[__alo],%%rax	\n\t"\
					"movq	%[__ahi],%%rbx	\n\t"\
					"vmovaps	%%ymm0,(%%rax)	\n\t"/* lo */\
					"vmovaps	%%ymm1,(%%rbx)	\n\t"/* hi */\
					"movq	%[__alo_norm],%%rax	\n\t"\
					"movq	%[__ahi_norm],%%rbx	\n\t"\
					"movq	%[__acy],%%rcx	\n\t"\
					"     vmulpd	%%ymm1,%%ymm14,%%ymm2	\n\t"/* tmp = hi*BINV */\
					"vroundpd	$0,%%ymm2,%%ymm2			\n\t"/* hh = DNINT(hi*BINV) */\
					"vmovaps	%%ymm2,%%ymm3				\n\t"/* cpy hh into cy-reg */\
					"vfnmadd213pd	%%ymm1,%%ymm15,%%ymm3	\n\t"/* cy = FMA(hh,-BASE, hi)	hi - hh*BASE = 'backward carry' from hi into lo, needed for proper base-normalization */\
					"vaddpd			%%ymm3,%%ymm0,%%ymm1	\n\t"/* lo += cy */\
					"vmovaps	%%ymm1,(%%rax)	\n\t"/* lo, base-normalized */\
					"vmovaps	%%ymm2,(%%rbx)	\n\t"/* hh = hi, base-normalized */\
					"vmovaps	%%ymm3,(%%rcx)	\n\t"/* cy */\
					:					/* outputs: none */\
					: [__ax] "m" (ax)	/* All inputs from memory addresses here */\
					 ,[__ay] "m" (ay)\
					 ,[__alo] "m" (alo)\
					 ,[__ahi] "m" (ahi)\
					 ,[__acy] "m" (acy)\
					 ,[__alo_norm] "m" (alo_norm)\
					 ,[__ahi_norm] "m" (ahi_norm)\
					 ,[__base] "m" (dptr1)\
					 ,[__binv] "m" (dptr2)\
					: "cc","memory","rax","rbx","rcx","xmm0","xmm1","xmm2","xmm3","xmm4","xmm14","xmm15"		/* Clobbered registers */\
				);
			//	printf("i = %u: x = %1.0f; %1.0f; lo = %1.0f; hi = %1.0f\n",i,*ax,*ay,*alo,*ahi);
				for(tmp = acy; tmp < acy+4; tmp++) {
					dtmp = *tmp;	// preserve signs
					if(fabs(dtmp) > fabs(cy_max)) cy_max = dtmp;
				}
				// Update log2-range-bounds-storing vars:
				for(tmp = alo; tmp < ahi; tmp++) {
					dblo = fabs(*tmp); dbhi = fabs(*(tmp+4));
					// Odds of a 0 operand are slim, but if() around it anyway:
					if(dblo) { dblo = log(dblo)*ILG2;	if(dblo > l2lo) l2lo = dblo; }
					if(dbhi) { dbhi = log(dbhi)*ILG2;	if(dbhi > l2hi) l2hi = dbhi; }
				}
			  #ifdef MUL_LOHI64_SUBROUTINE
				MUL_LOHI64(iax,iay,&ialo,&iahi);
				MUL_LOHI64(ibx,iby,&iblo,&ibhi);
				MUL_LOHI64(icx,icy,&iclo,&ichi);
				MUL_LOHI64(idx,idy,&idlo,&idhi);
			  #else
				MUL_LOHI64(iax,iay, ialo, iahi);
				MUL_LOHI64(ibx,iby, iblo, ibhi);
				MUL_LOHI64(icx,icy, iclo, ichi);
				MUL_LOHI64(idx,idy, idlo, idhi);
			  #endif
			  /*
				if(pow2 == 53 && i < 100) {
					printf("I = %d: ax = %llu ay = %llu ahi,alo = %f,%f\n",i, *ax,*ay, *ahi,*alo);
					printf("I = %d: bx = %llu by = %llu bhi,blo = %f,%f\n",i, *bx,*by, *bhi,*blo);
					printf("I = %d: cx = %llu cy = %llu chi,clo = %f,%f\n",i, *cx,*cy, *chi,*clo);
					printf("I = %d: dx = %llu dy = %llu dhi,dlo = %f,%f\n",i, *dx,*dy, *dhi,*dlo);
				}
			  */
				if(cmp_fma_lohi_vs_exact(*ax,*ay,*ahi,*alo, iax,iay,iahi,ialo)) { ++nerr; printf("ERROR: pow2 = %d, I = %d, A-outputs differ!\n",pow2,i); ASSERT(HERE, 0, "fma_dmult tests failed!"); }
				if(cmp_fma_lohi_vs_exact(*bx,*by,*bhi,*blo, ibx,iby,ibhi,iblo)) { ++nerr; printf("ERROR: pow2 = %d, I = %d, B-outputs differ!\n",pow2,i); ASSERT(HERE, 0, "fma_dmult tests failed!"); }
				if(cmp_fma_lohi_vs_exact(*cx,*cy,*chi,*clo, icx,icy,ichi,iclo)) { ++nerr; printf("ERROR: pow2 = %d, I = %d, C-outputs differ!\n",pow2,i); ASSERT(HERE, 0, "fma_dmult tests failed!"); }
				if(cmp_fma_lohi_vs_exact(*dx,*dy,*dhi,*dlo, idx,idy,idhi,idlo)) { ++nerr; printf("ERROR: pow2 = %d, I = %d, D-outputs differ!\n",pow2,i); ASSERT(HERE, 0, "fma_dmult tests failed!"); }
			  #elif 0
				#error to-do!
				double r1,r2, lo,hi;
				r1 = rng_isaac_rand_double_norm_pm1() * pow2_dmult;	// in [-2^50, +2^50]
				r2 = rng_isaac_rand_double_norm_pm1() * pow2_dmult;	// in [-2^50, +2^50]
				mul50x50_debug(r1,r2, &lo,&hi);
				printf("mul50x50_: a,b = %llu, %llu\n",*(uint64*)&r1,*(uint64*)&r2);
				printf("mul50x50_: lo = %16llu\n",*(uint64*)alo);
				printf("mul50x50_: hi = %16llu\n",*(uint64*)ahi);
			  #endif

			#endif	// 0?

			/******************** experimental code: Try squaring [lo,hi] (in ymm1,2), sans intermediate base-normalizations: *******************/
				__asm__ volatile (\
					/* Use ymm[< 16] for broadcast-from-scalar consts to avoid illegal-instruction exceptions (AVX-512F supports b'cast-to-full-width-zmm) */\
					"movq	%[__base],%%rax		\n\t	vbroadcastsd	(%%rax),%%ymm15	\n\t"/* BASE */\
					"movq	%[__binv],%%rbx		\n\t	vbroadcastsd	(%%rbx),%%ymm14	\n\t"/* BINV */\
					"movq	%[__crnd50],%%rcx	\n\t	vbroadcastsd	(%%rcx),%%ymm13	\n\t"/* CRND*2^50 */\
					"movq	%[__prod1_adj],%%rdx\n\t	vbroadcastsd	(%%rdx),%%ymm12	\n\t"/* Const to multiply by base and add to prod[1] to ensure latter >= 0 */\
					"movq	%[__ax] ,%%rax	\n\t	vmovaps	(%%rax),%%ymm1			\n\t"/* x.lo */\
					"movq	%[__ay] ,%%rbx	\n\t"/* x.hi */\
					"movq	%[__alo],%%rcx	\n\t"\
					"vaddpd			%%ymm1,%%ymm1,%%ymm2	\n\t"/* 2*lo */\
				/* lo*lo: */
					"     vmulpd	%%ymm1,%%ymm1,%%ymm3	\n\t"/* hi = fma(a,b, 0  ) - use VMULPD instead of VFMADD231PD, since addend = 0. */
					"vmovaps	%%ymm3,%%ymm0				\n\t"/* cpy hi into lo-reg */\
					"vfmsub231pd	%%ymm1,%%ymm1,%%ymm0	\n\t"/* lo = fma(a,b, -hi) */
					"vaddpd			%%ymm13,%%ymm3,%%ymm1	\n\t"\
					"vsubpd			%%ymm13,%%ymm1,%%ymm1	\n\t"/* hh = hi +- crnd50 to round-to-nearest-multiple-of-2^50 */\
					"vsubpd			%%ymm1,%%ymm3,%%ymm3	\n\t"/* hi - hh gives backward-carry... */\
					"vaddpd			%%ymm3,%%ymm0,%%ymm0	\n\t"/* ...which we add to lo (prod[0]). */\
					"vmovaps	%%ymm0,    (%%rcx)	\n\t"/* write prod[0] */\
														/*** prod[0] in ymm0, hi*2^50 in ymm1. ***/\
				/* 2*lo*hi: */
					"     vmulpd	(%%rbx),%%ymm2,%%ymm3	\n\t"/* hi = fma(a,b, 0  ) - use VMULPD instead of VFMADD231PD, since addend = 0. */
					"vmovaps	%%ymm3,%%ymm0				\n\t"/* cpy hi into lo-reg */\
					"vfmsub231pd	(%%rbx),%%ymm2,%%ymm0	\n\t"/* lo = fma(a,b, -hi) */
					"vaddpd			%%ymm13,%%ymm3,%%ymm2	\n\t"\
					"vsubpd			%%ymm13,%%ymm2,%%ymm2	\n\t"/* hh = hi +- crnd50 to round-to-nearest-multiple-of-2^50 */\
					"vfmadd231pd	%%ymm14,%%ymm1,%%ymm0	\n\t"/* Add lo to hi-output of previous lo:hi, pair, which also needs *= binv */
					"vsubpd			%%ymm2,%%ymm3,%%ymm3	\n\t"/* hi - hh gives backward-carry... */\
					"vaddpd			%%ymm3,%%ymm0,%%ymm1	\n\t"/* ...which we add to lo (prod[1]). */\
				"vfmadd231pd	%%ymm12,%%ymm15,%%ymm1	\n\t"/* Add const*base and add to prod[1] to ensure latter >= 0 */
				"vfnmadd231pd	%%ymm12,%%ymm15,%%ymm2	\n\t"/* Must sub same const from prod[2] by way of a carry - but note hh still scaled *= base. */\
					"vmovaps	%%ymm1,0x20(%%rcx)	\n\t"/* write prod[1] */\
														/*** prod[1] in ymm1, hi*2^50 in ymm2. ***/\
				/* hi*hi: */
					"vmovaps	(%%rbx),%%ymm3				\n\t"/* reload hi into ymm3 */\
					"     vmulpd	%%ymm3,%%ymm3,%%ymm0	\n\t"/* hi = fma(a,b, 0  ) - use VMULPD instead of VFMADD231PD, since addend = 0. */
					"vmovaps	%%ymm0,%%ymm1				\n\t"/* cpy hi into lo-reg */\
					"vfmsub231pd	%%ymm3,%%ymm3,%%ymm1	\n\t"/* lo = fma(a,b, -hi) */
					"vaddpd			%%ymm13,%%ymm0,%%ymm3	\n\t"\
					"vsubpd			%%ymm13,%%ymm3,%%ymm3	\n\t"/* hh = hi +- crnd50 to round-to-nearest-multiple-of-2^50 */\
					"vfmadd231pd	%%ymm14,%%ymm2,%%ymm1	\n\t"/* Add lo to hi-output of previous lo:hi, pair, which also needs *= binv */
					"vsubpd			%%ymm3,%%ymm0,%%ymm0	\n\t"/* hi - hh gives backward-carry... */\
					"vaddpd			%%ymm0,%%ymm1,%%ymm2	\n\t"/* ...which we add to lo (prod[2]). */\
					"     vmulpd	%%ymm3,%%ymm14,%%ymm3	\n\t"/* prod[3] = hh*binv */\
					"vmovaps	%%ymm2,0x40(%%rcx)	\n\t"/* write prod[2] */\
					"vmovaps	%%ymm3,0x60(%%rcx)	\n\t"/* write prod[3] */\
														/*** prod[2,3] in ymm2,3. ***/\
					:					/* outputs: none */\
					: [__ax] "m" (ax)	/* All inputs from memory addresses here */\
					 ,[__ay] "m" (ay)\
					 ,[__alo] "m" (alo)\
					 ,[__base] "m" (dptr1)\
					 ,[__binv] "m" (dptr2)\
					 ,[__crnd50] "m" (dptr3)\
					 ,[__prod1_adj] "m" (dptr4)\
					: "cc","memory","rax","rbx","rcx","xmm0","xmm1","xmm2","xmm3","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
				);
			//	printf("i = %u: x0 = %1.0f; x1 = %1.0f; p0-3 = %1.0f,%1.0f,%1.0f,%1.0f\n",i,*ax,*ay,*alo,*ahi);
				// Update log2-range-bounds-storing vars:
				tmp = alo;
				for(j = 0; j < 4; j++) {
					for(k = 0; k < 4; k++,tmp++) {
						dtmp = *tmp;
						if(dtmp < sqr100lo[j]) { sqr100lo[j] = dtmp; }
						if(dtmp > sqr100hi[j]) { sqr100hi[j] = dtmp; }
					}
				}

			}	// i-loop
		#if 0
			printf("\t%u |Outputs|.lo for %d x %d-bit fma_dmult have l2max = %10.7f:\n",i,pow2,pow2,l2lo);
			printf("\t%u |Outputs|.hi for %d x %d-bit fma_dmult have l2max = %10.7f:\n",i,pow2,pow2,l2hi);
			// Use 1.0f as format - .0 means no fractional part, and i/o routines will override the length-1 with actual length:
			if(cy_max > 0) {
				itmp64 = cy_max; itmp32 = trailz64(itmp64); itmp64 >>= itmp32;
				printf("\tcy_max = %1.0f =  %llu * 2^%u\n",cy_max,itmp64,itmp32);
			} else if(cy_max < 0) {
				itmp64 =-cy_max; itmp32 = trailz64(itmp64); itmp64 >>= itmp32;
				printf("\tcy_max = %1.0f = -%llu * 2^%u\n",cy_max,itmp64,itmp32);
			} else {
				printf("\tcy_max =  0\n");
			}
		#endif
			printf("SQR outputs have the following ranges:\n");
			dtmp = pow2_imult;
			for(j = 0; j < 4; j++) {
				printf("\tprod[%u]: [%1.2f, %1.2f] = [%5.2f, %5.2f]*base\n",j,sqr100lo[j],sqr100hi[j],sqr100lo[j]*dtmp,sqr100hi[j]*dtmp);
			}
			printf("\n");

			pow2_dmult *= 2.0;
			pow2_imult *= 0.5;
		}	// pow2-loop
exit(0);
		return nerr;
	}

/* Generic version for non-CUDA-only FMA:

void	mul50x50_debug(double a, double b, double *lo, double *hi)
	{
		// Exact product a*b = lo + hi:
		*hi = fma(a,b, 0   );
		*lo = fma(a,b, -*hi);
	}
*/

	/*
	Example: FMA in-and-outputs are
	x = -163183843911180; y = 1039312675530994; hi,lo = -169599037418760590289388175360, 69115062440,
                          sum to yield exact result:  = -169599037418760590220273112920
                        now write in base-2^64 form:  = -(9193982241*2^64 + 802977670696261464) ,
	We then simply compare result inside () to MUL_LOHI64 result. In practice the assembly requires a bit more work,
	but the principle is simple.
	*/
	int cmp_fma_lohi_vs_exact(double dx, double dy, double dhi, double dlo, uint64 ix, uint64 iy, uint64 ihi, uint64 ilo)
	{
		int retval;
		uint64 i64, e1,e0, m1,m0;
		uint32 s1,s0;
		const uint64 two52 = 0x0010000000000000ull, mmask = 0x000FFFFFFFFFFFFFull;
		uint128 exact;
		const char char_sgn[2] = {' ','-'};
	//printf("I = %d: x = %f; y = %f; hi,lo = %f,%f\n",i, dx,dy, dhi,dlo);
		if(dx == 0. || dy == 0.)	// Comparison algo needs further tweaks to handle 0-result ... not worth coding time.
			return 0;
		s1 = (dhi < 0);	// Sign of product = sign of hi output
		s0 = (dlo < 0);
		if(s1) {	// If product < 0, negate both FMA outputs prior to comparing vs the (unsigned-int) MUL_LOHI64 result
			dhi = -dhi;
			if(dlo != 0.)	// Need to preserve sign of dlo, if = 0
				dlo = -dlo;
		}
		// Extract exp & mantissa fields of the double outputs and restore hidden bits:
		i64 = *(uint64*)&dhi; e1 = (i64>>52)&0x7ff; m1 = (i64&mmask) + two52;
		i64 = *(uint64*)&dlo; e0 = (i64>>52)&0x7ff; m0 = (i64&mmask) +(two52 & (-(dlo != 0.)));
		int nsh1 = e1 - 0x433;	// Shift count of hi-double result = exp - 0x3ff - 52
		int nsh0 = e0 - 0x433;	// Shift count of lo-double result
		exact.d0 = m1; exact.d1 = 0;	LSHIFT128(exact,nsh1, exact);
		if(nsh0 < 0)
			m0 >>= -nsh0;
		else if(nsh0 > 0)
			m0 <<=  nsh0;
		if(s1 ^ s0) {	// Whether to add or sub the lo term depends on the *relative* signs of hi,lo outputs
			i64 = exact.d0; exact.d0 -= m0; exact.d1 -= (exact.d0 > i64);
		} else {
							exact.d0 += m0; exact.d1 += (exact.d0 < m0 );
		}
		retval = ( (ihi != exact.d1) || (ilo != exact.d0) );
		if(retval) {
			printf("In cmp_fma_lohi_vs_exact: FMA-double and pure-int DMUL results differ!\n");
			printf("dx = %f; dy = %f; hi,lo = %f,%f\n",dx,dy, dhi * (1 - 2*(s1 != 0)), dlo * (1 - 2*(s0 != 0)));
			printf("ix = %lld; iy = %lld; ihi,lo = %lld,%llu\n",ix,iy, ihi,ilo);
			printf("Unsigned FMA result: ihi = %llX; ilo = %llX\n",*(uint64*)&dhi,*(uint64*)&dlo);
			printf("nsh1,0 = %d,%d: ehi = %llu; elo = %llu [mlo = %c%llu]\n",nsh1,nsh0,exact.d1,exact.d0, char_sgn[s1 ^ s0],m0);
		}
		return retval;
	}

#endif

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

/******* Supplements to stdio (e.g. binary-formatted output) *********/

void byte_bitstr(const uint8 byte, char*ostr)
{
	const char bytestr[256][10] = {
	"00000000 ","00000001 ","00000010 ","00000011 ","00000100 ","00000101 ","00000110 ","00000111 ","00001000 ","00001001 ","00001010 ","00001011 ","00001100 ","00001101 ","00001110 ","00001111 ",
	"00010000 ","00010001 ","00010010 ","00010011 ","00010100 ","00010101 ","00010110 ","00010111 ","00011000 ","00011001 ","00011010 ","00011011 ","00011100 ","00011101 ","00011110 ","00011111 ",
	"00100000 ","00100001 ","00100010 ","00100011 ","00100100 ","00100101 ","00100110 ","00100111 ","00101000 ","00101001 ","00101010 ","00101011 ","00101100 ","00101101 ","00101110 ","00101111 ",
	"00110000 ","00110001 ","00110010 ","00110011 ","00110100 ","00110101 ","00110110 ","00110111 ","00111000 ","00111001 ","00111010 ","00111011 ","00111100 ","00111101 ","00111110 ","00111111 ",
	"01000000 ","01000001 ","01000010 ","01000011 ","01000100 ","01000101 ","01000110 ","01000111 ","01001000 ","01001001 ","01001010 ","01001011 ","01001100 ","01001101 ","01001110 ","01001111 ",
	"01010000 ","01010001 ","01010010 ","01010011 ","01010100 ","01010101 ","01010110 ","01010111 ","01011000 ","01011001 ","01011010 ","01011011 ","01011100 ","01011101 ","01011110 ","01011111 ",
	"01100000 ","01100001 ","01100010 ","01100011 ","01100100 ","01100101 ","01100110 ","01100111 ","01101000 ","01101001 ","01101010 ","01101011 ","01101100 ","01101101 ","01101110 ","01101111 ",
	"01110000 ","01110001 ","01110010 ","01110011 ","01110100 ","01110101 ","01110110 ","01110111 ","01111000 ","01111001 ","01111010 ","01111011 ","01111100 ","01111101 ","01111110 ","01111111 ",
	"10000000 ","10000001 ","10000010 ","10000011 ","10000100 ","10000101 ","10000110 ","10000111 ","10001000 ","10001001 ","10001010 ","10001011 ","10001100 ","10001101 ","10001110 ","10001111 ",
	"10010000 ","10010001 ","10010010 ","10010011 ","10010100 ","10010101 ","10010110 ","10010111 ","10011000 ","10011001 ","10011010 ","10011011 ","10011100 ","10011101 ","10011110 ","10011111 ",
	"10100000 ","10100001 ","10100010 ","10100011 ","10100100 ","10100101 ","10100110 ","10100111 ","10101000 ","10101001 ","10101010 ","10101011 ","10101100 ","10101101 ","10101110 ","10101111 ",
	"10110000 ","10110001 ","10110010 ","10110011 ","10110100 ","10110101 ","10110110 ","10110111 ","10111000 ","10111001 ","10111010 ","10111011 ","10111100 ","10111101 ","10111110 ","10111111 ",
	"11000000 ","11000001 ","11000010 ","11000011 ","11000100 ","11000101 ","11000110 ","11000111 ","11001000 ","11001001 ","11001010 ","11001011 ","11001100 ","11001101 ","11001110 ","11001111 ",
	"11010000 ","11010001 ","11010010 ","11010011 ","11010100 ","11010101 ","11010110 ","11010111 ","11011000 ","11011001 ","11011010 ","11011011 ","11011100 ","11011101 ","11011110 ","11011111 ",
	"11100000 ","11100001 ","11100010 ","11100011 ","11100100 ","11100101 ","11100110 ","11100111 ","11101000 ","11101001 ","11101010 ","11101011 ","11101100 ","11101101 ","11101110 ","11101111 ",
	"11110000 ","11110001 ","11110010 ","11110011 ","11110100 ","11110101 ","11110110 ","11110111 ","11111000 ","11111001 ","11111010 ","11111011 ","11111100 ","11111101 ","11111110 ","11111111 "
	};
	strcpy(ostr, bytestr[byte]);
}

void	ui32_bitstr(const uint32 ui32, char*ostr)
{
	int i;
	for(i = 0; i < 32; i += 8) {
		// High byte ==> leftmost 8 chars of output string, thus the (24 - i)
		// Also note the padding byte on ostr to account for the blank byte_bitstr output terminators:
		byte_bitstr((uint8)(ui32 >> (24 - i)), ostr + i + (i>>3));
	}
}

void	ui64_bitstr(const uint64 ui64, char*ostr)
{
	int i;
	for(i = 0; i < 64; i += 8) {
		// High byte ==> leftmost 8 chars of output string, thus the (56 - i)
		// Also note the padding byte on ostr to account for the blank byte_bitstr output terminators:
		byte_bitstr((uint8)(ui64 >> (56 - i)), ostr + i + (i>>3));
	}
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

/******* Bit-level utilities: ********/

// 32 and 64-bit leftward circular shift, shift count n assumed unsigned < #bits-in-type:
DEV uint32 cshft32(uint32 x, uint32 n)
{
	if(n)
		return (x << n) + (x >> (32-n));
	else
		return x;
}

DEV uint64 cshft64(uint64 x, uint64 n)
{
	if(n)
		return (x << n) + (x >> (64-n));
	else
		return x;
}

// 32 and 64-bit analogs of the F90 intrinsic ISHFT function:
DEV uint32 ishft32(uint32 x, int shift)
{
	uint32 r;
	if(shift >= 32)
		r  = 0ull;
	else if(shift > 0)
		r  = x << shift;
	else if(shift > -32)
		r  = x >> (-shift);
	else
		r  = 0ull;
	return r;
}

DEV uint64 ishft64(uint64 x, int shift)
{
	uint64 r;
	if(shift > 64)
		r  = 0ull;
	else if(shift > 0)
		r  = x << shift;
	else if(shift > -64)
		r  = x >> (-shift);
	else
		r  = 0ull;
	return r;
}

// Clears [bit]th of array [arr]. No bounds checking is performed.
void bit_clr32(uint32*arr, uint32 bit)
{
#ifdef X32_ASM
	__asm__ volatile (\
		"btrl  %1, %0		\n\t"\
		:	/* outputs: none */\
		:  "m" (*arr), "r" (bit)\
		: "cc","memory"	/* Clobbered registers */\
		);
#else
	uint32 n = bit>>5;
	uint32 mask = ~((uint32)1 << (bit&31));
	arr[n] &= mask;
#endif
}

void bit_clr32_x4(uint32*arr, uint32 bit1, uint32 bit2, uint32 bit3, uint32 bit4)
{
#ifdef X32_ASM
	__asm__ volatile (\
		"btrl  %1, %0		\n\t"\
		"btrl  %2, %0		\n\t"\
		"btrl  %3, %0		\n\t"\
		"btrl  %4, %0		\n\t"\
		:	/* outputs: none */\
		:  "m" (*arr), "r" (bit1), "r" (bit2), "r" (bit3), "r" (bit4)\
		: "cc","memory"	/* Clobbered registers */\
		);
#else
	uint32 n1 = bit1>>5, n2 = bit2>>5, n3 = bit3>>5, n4 = bit4>>5;
	uint32 mask1 = ~((uint32)1 << (bit1&31)),mask2 = ~((uint32)1 << (bit2&31)), mask3 = ~((uint32)1 << (bit3&31)),mask4 = ~((uint32)1 << (bit4&31));
	arr[n1] &= mask1;
	arr[n2] &= mask2;
	arr[n3] &= mask3;
	arr[n4] &= mask4;
#endif
}

// Clears [bit]th (assumed in <0:31>) bit of [n]th word of array [arr]. No bounds checking is performed.
void bit_clr64(uint64*arr, uint64 bit)
{
#ifdef X64_ASM
	__asm__ volatile (\
		"btrq  %0, %1		\n\t"\
		:	/* outputs: none */\
		: "r" (bit), "m" (*arr)\
		: "cc","memory"	/* Clobbered registers */\
		);
#else
	uint64 n = bit>>6;
	uint64 mask = ~((uint64)1 << (bit&63));
	arr[n] &= mask;
#endif
}

void bit_clr64_x4(uint64*arr, uint64 bit1, uint64 bit2, uint64 bit3, uint64 bit4)
{
#ifdef X64_ASM
	__asm__ volatile (\
		"btrq  %1, %0		\n\t"\
		"btrq  %2, %0		\n\t"\
		"btrq  %3, %0		\n\t"\
		"btrq  %4, %0		\n\t"\
		:	/* outputs: none */\
		:  "m" (*arr), "r" (bit1), "r" (bit2), "r" (bit3), "r" (bit4)\
		: "cc","memory"	/* Clobbered registers */\
		);
#else
	uint64 n1 = bit1>>6, n2 = bit2>>6, n3 = bit3>>6, n4 = bit4>>6;
	uint64 mask1 = ~((uint64)1 << (bit1&63)),mask2 = ~((uint64)1 << (bit2&63)), mask3 = ~((uint64)1 << (bit3&63)),mask4 = ~((uint64)1 << (bit4&63));
	arr[n1] &= mask1;
	arr[n2] &= mask2;
	arr[n3] &= mask3;
	arr[n4] &= mask4;
#endif
}

// Counts set bits in 32-bit int x - only tiny speedup on Haswell from using 32-bit POPCNT, so no ASM here:
DEV uint32 popcount32(uint32 x)
{
	uint8 *byte_arr = (uint8*)&x;
	uint32 i,retval = 0;
	for(i = 0; i < 4; i++) {
		retval += pop8[byte_arr[i]];
	}
	return retval;
}

// Counts set bits in 64-bit int x - ~30% speedup on Haswell from using 64-bit POPCNT:
DEV uint32 popcount64(uint64 x)
{
// Intel introduced POPCNT at same time as SSE4.2, but "...not considered part of the SSE4.2
// instruction set; instead, they have their own dedicated CPUID bits to indicate support",
// so avoid dealing with that stupidity by only enabling in AVX-and-beyond builds:
#ifdef USE_AVX
	uint64 i64 = 0;
	__asm__ volatile ("popcntq %1,%0" : "=r" (i64) : "r" (x));
	return (uint32)i64;
#else
	uint8 *byte_arr = (uint8*)&x;
	uint32 i,retval = 0;
	// May 2018: Unrolling the for-loop in favor of an inlined 8-fold sum gave a nice speedup:
  #if 0
	for(i = 0; i < 8; i++) {
		retval += pop8[byte_arr[i]];
	}
  #else
	retval = pop8[byte_arr[0]] + pop8[byte_arr[1]] + pop8[byte_arr[2]] + pop8[byte_arr[3]]
		   + pop8[byte_arr[4]] + pop8[byte_arr[5]] + pop8[byte_arr[6]] + pop8[byte_arr[7]];
  #endif
	return retval;
#endif
}

// Return bit position [0:31] of the [bit = 1:32]th set bits in 32-bit int x, or -1 if there are fewer than
// [bit] set bits in x, or if bit == 0 (i.e. user requests position of nonexistent "0th set bit"):
int ith_set_bit32(uint32 x, uint32 bit)
{
	uint8 curr_byte;
	int curr_pop,i,j,k,retval = 0;
	if(!x || !bit) return -1;
	ASSERT(HERE, bit <= 32, "[bit]th-bit specifier out of range!");
	// Find the byte in which the [bit]th set-bit occurs:
	for(i = 0; i < 32; i += 8) {
		curr_byte = (uint8)(x >> i);
		curr_pop = pop8[curr_byte];
		retval += curr_pop;	// At this point retval stores the popcount-to-date...
		// ... If that >= [bit], replace that with the location of the [bit]th set bit:
		if(retval >= bit) {
			retval -= curr_pop;	// Subtract pop of curr_byte back off
			k = (bit-retval-1)*4;	// Need to sub-1 since e.g. 3rd set-bit is encoded in hex-char 2
			j = (ith_set_bit8[curr_byte] >> k) & 0xf;
			return i+j;
		}
	}
	return -1;
}

// 64-bit version of ith_set_bit32:
// Remember, bit-index in arglist is *unit* offset, i.e. must be in [1:64]!
int ith_set_bit64(uint64 x, uint32 bit)
{
	uint8 curr_byte;
	int curr_pop,i,j,k,retval = 0;
	if(!x || !bit) return -1;
	ASSERT(HERE, bit <= 64, "[bit]th-bit specifier out of range!");
	// Find the byte in which the [bit]th set-bit occurs:
	for(i = 0; i < 64; i += 8) {
		curr_byte = (uint8)(x >> i);
		curr_pop = pop8[curr_byte];
		retval += curr_pop;	// At this point retval stores the popcount-to-date...
		// ... If that >= [bit], replace that with the location of the [bit]th set bit:
		if(retval >= bit) {
			retval -= curr_pop;	// Subtract pop of curr_byte back off
			// On Core2, ith_set_bit8-version cuts 40% off the function runtime vs loop-over-bits-of-curr_byte:
		  #if 0
			for(j = 0; j < 8; j++) {	// ...and step 1-bit-at-a-time until reach desired popcount
				retval += (curr_byte>>j)&1;
				if(retval == bit) return i+j;
			}
		  #else
			k = (bit-retval-1)*4;	// Need to sub-1 since e.g. 3rd set-bit is encoded in hex-char 2
			j = (ith_set_bit8[curr_byte] >> k) & 0xf;
			return i+j;
		  #endif
		}
	}
	return -1;
}

/*** leading and trailing-zero-counting algorithms: ***/
DEV uint32 trailz32(uint32 x)
{
	uint8 *byte_arr = (uint8*)&x, curr_byte;
	uint32 i,retval = 0;
  #ifdef USE_BIG_ENDIAN
	for(i = 0; i < 4; i++) {
		curr_byte = byte_arr[3-i];
		retval += tz8[curr_byte];
		if(curr_byte) break;
	}
  #else
	for(i = 0; i < 4; i++) {
		curr_byte = byte_arr[i];
		retval += tz8[curr_byte];
		if(curr_byte) break;
	}
  #endif
	return retval;
}

DEV uint32 trailz64(uint64 x)
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
	uint8 *byte_arr = (uint8*)&x, curr_byte;
	uint32 i,retval = 0;
  #ifdef USE_BIG_ENDIAN
	for(i = 0; i < 8; i++) {
		curr_byte = byte_arr[7-i];
		retval += tz8[curr_byte];
		if(curr_byte) break;
	}
  #else
	for(i = 0; i < 8; i++) {
		curr_byte = byte_arr[i];
		retval += tz8[curr_byte];
		if(curr_byte) break;
	}
  #endif
	return retval;
#endif
}

/***************/

DEV uint32 leadz32(uint32 x)
{
#ifdef X32_ASM
	uint32 lz;
	int bpos;
	if(x == 0) return 32;
	__asm__ volatile (\
		"bsrl %[__x],%%eax		\n\t"\
		"movl %%eax,%[__bpos]	\n\t"\
		:	/* outputs: none */\
		: [__x] "m" (x)	/* All inputs from memory addresses here */\
		 ,[__bpos] "m" (bpos)	\
		: "cc","memory","eax"	/* Clobbered registers */\
	);
	lz = (31 - bpos);	// BSR returns *index* of leftmost set bit, must subtract from (#bits - 1) to get #lz.
	return lz;
#else
	uint8 *byte_arr = (uint8*)&x, curr_byte;
	uint32 i,retval = 0;
  #ifdef USE_BIG_ENDIAN
	for(i = 0; i < 4; i++) {
		curr_byte = byte_arr[i];
		retval += lz8[curr_byte];
		if(curr_byte) break;
	}
  #else
	for(i = 0; i < 4; i++) {
		curr_byte = byte_arr[3-i];
		retval += lz8[curr_byte];
		if(curr_byte) break;
	}
  #endif
	return retval;
#endif
}

DEV uint32 leadz64(uint64 x)
{
#ifdef X64_ASM
	uint32 lz;
	int bpos;
	if(x == 0) return 64;
	__asm__ volatile (\
		"bsrq %[__x],%%rax		\n\t"\
		"movl %%eax,%[__bpos]	\n\t"\
		:	/* outputs: none */\
		: [__x] "m" (x)	/* All inputs from memory addresses here */\
		 ,[__bpos] "m" (bpos)	\
		: "cc","memory","rax"	/* Clobbered registers */\
	);
	lz = (63 - bpos);	// BSR returns *index* of leftmost set bit, must subtract from (#bits - 1) to get #lz.
	return lz;
#else
	uint8 *byte_arr = (uint8*)&x, curr_byte;
	uint32 i,retval = 0;
  #ifdef USE_BIG_ENDIAN
	for(i = 0; i < 8; i++) {
		curr_byte = byte_arr[i];
		retval += lz8[curr_byte];
		if(curr_byte) break;
	}
  #else
	for(i = 0; i < 8; i++) {
		curr_byte = byte_arr[7-i];
		retval += lz8[curr_byte];
		if(curr_byte) break;
	}
  #endif
	return retval;
#endif
}

DEV uint32	leadz128(uint128 i)
{
	if(i.d1)
		return leadz64(i.d1);
	else
		return leadz64(i.d0) + 64;
}

DEV uint32	leadz192(uint192 i)
{
	if(i.d2)
		return leadz64(i.d2);
	else if(i.d1)
		return leadz64(i.d1) + 64;
	else
		return leadz64(i.d0) + 128;
}

DEV uint32	leadz256(uint256 i)
{
	if(i.d3)
		return leadz64(i.d3);
	else if(i.d2)
		return leadz64(i.d2) + 64;
	else if(i.d1)
		return leadz64(i.d1) + 128;
	else
		return leadz64(i.d0) + 192;
}

/***************/

/* If the input == 2^p, returns 1; otherwise returns 0.
Note: We do not consider 0 as being a power of 2, since there is no number x such that 2^x = 0.
Algo: for unsigned twos-comp int n,
  n&(n-1)  ('is n&(n-1) nonzero?' in boolean terms) is only 0 (false) for n a power-of-2, with exception of n = 0.
Conversely,
!(n&(n-1)) ('is n&(n-1) zero?' in boolean terms) is only !0 (true) for n a power-of-2, with exception of n = 0,
which latter instance is easily special-cased, via n && !(n&(n-1)).

On hardware with a fast pop-count instruction can also consider using (popcount == 1)?, but not worth the tiny cycle savings.
*/
DEV uint32 isPow2(uint32 i32)
{
	return i32 && !(i32&(i32-1));
}

DEV uint64 isPow2_64(uint64 i64)
{
	return i64 && !(i64&(i64-1));
}

/***************/

/* If the input == 4^p, returns 1; otherwise returns 0. */
DEV uint32 isPow4(uint32 i32)
{
	return isPow2(i32) && (i32 & 0x55555555);
}

DEV uint64 isPow4_64(uint64 i64)
{
	return isPow2_64(i64) && (i64 & 0x5555555555555555ull);
}

/***************/

/* Extract (nbits) bits beginning at position (beg) */

DEV uint32 ibits32(uint32 i, uint32 beg, uint32 nbits)
{
	uint32 ones_mask = 0xFFFFFFFF;
	return ( (i >> beg) & ~(ones_mask << nbits) );
}

DEV uint64 ibits64(uint64 i, uint32 beg, uint32 nbits)
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
DEV uint64	getbits64(uint64 x, uint32 src_bit_start, uint32 nbits, uint32 tgt_bit_start)
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
DEV void	mvbits64(uint64 x, uint32 src_bit_start, uint32 nbits, uint64*y, uint32 tgt_bit_start)
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

/* returns 1 if p is a base-z Fermat pseudoprime, 0 otherwise. */
DEV int pprimeF(uint32 p, uint32 base)
{
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

// 64-bit analog of pprimeF:
DEV int pprimeF64(uint64 p, uint64 base)
{
	uint64 y = 1ull, n = p-1ull, flag;
	uint64 z = base;

	while(n)
	{
		flag = n & 1ull;
		n >>= 1;
		if(flag) y = mi64_modmul64(y,z,p);
		z = mi64_modmul64(z,z,p);
		if(!z) return 0;
	}
	return((int)(y==1ull));
}

/***************/

DEV int isPRP(uint32 p)
{
	// Handle even-argument case separately, since the powmod routines may not accept even moduli:
	if((p & 1) == 0)
		return (p == 2);
	/* TODO: replace/supplement this with a rigorous trial-divide test for p < 2^32 */
	return(pprimeF(p,2) && pprimeF(p,3) && pprimeF(p,5) && pprimeF(p,7) && pprimeF(p,11) && pprimeF(p,13));
}

#include "factor.h"	// Needed for twopmodq64() prototype
DEV int isPRP64(uint64 p)
{
	// Handle even-argument case separately, since the powmod routines may not accept even moduli:
	if((p & 1ull) == 0ull)
		return (p == 2ull);
	return twopmodq64(p-1,p) == 1ull;
//	return(pprimeF64(p,2ull) && pprimeF64(p,3ull) && pprimeF64(p,5ull) && pprimeF64(p,7ull) && pprimeF64(p,11ull) && pprimeF64(p,13ull));
}

/*******************/

/* Calculate 2^-p mod q for p, q 32-bit unsigned ints. This can be used (among
other things) to effect a fast Fermat base-2 pseudoprime test, by calling with q = p-1.
*/
// V1 returns the full powmod result:
DEV uint32 twompmodq32(uint32 p, uint32 q)	// 2^-p % q
{
	 int32 j;
	uint32 lead5, pshift, qhalf, qinv, zshift, start_index, x, lo, hi;

	ASSERT(HERE, (q&1) == 1, "twompmodq32: even modulus!");
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
		/* SQR_LOHI32(x,lo,hi): */
		MUL_LOHI32(x,x, lo,hi);
		lo *= qinv;
		MULH32(q,lo, lo);

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

// V2 returns binary == 0? of powmod result:
DEV int twopmodq32(uint32 p, uint32 q)	// (2^-p % q) == 0
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
	} else {
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
		if(x > qhalf) {
			x += x;
			x -= q;
		} else {
			x += x;
		}
	}

	for(j = start_index-2; j >= 0; j--)
	{
		/* SQR_LOHI32(x,lo,hi): */
		MUL_LOHI32(x,x, lo,hi);
		lo *= qinv;
		MULH32(q,lo, lo);

		/* Branchless version is much faster, but less readable, so give the branched one inside a #if 0: */
	#ifdef NOBRANCH
		x = hi - lo + ((-(hi < lo)) & q);
	#else
		x = hi - lo;
		if(x > hi)
			x += q;	/* had a borrow */
	#endif

		if((pshift >> j) & (uint32)1)
		{
		/* Branchless version is much faster, but less readable, so give the branched one inside a #if 0: */
		#ifdef NOBRANCH
			x = x + x - ((-(x > qhalf)) & q);
		#else
			if(x > qhalf) {	/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
				x = x + x;
				x -= q;
			} else {
				x = x + x;
			}
		#endif
		}
	}

	/*...Double and return.	These are specialized for the case where 2^p == 1 mod q implies divisibility, in which case x = (q+1)/2. */
	return((int)((x + x - q) == 1));
}

/* Does an 8-fold base-2 PRP test on the prime candidates q0-7. */
DEV int twopmodq32_x8(uint32 q0, uint32 q1, uint32 q2, uint32 q3, uint32 q4, uint32 q5, uint32 q6, uint32 q7)
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

	/* lo = MULH32(q, lo): */
	MULH32(q0,lo0, lo0);
	MULH32(q1,lo1, lo1);
	MULH32(q2,lo2, lo2);
	MULH32(q3,lo3, lo3);
	MULH32(q4,lo4, lo4);
	MULH32(q5,lo5, lo5);
	MULH32(q6,lo6, lo6);
	MULH32(q7,lo7, lo7);

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
		/* SQR_LOHI32(x,lo,hi): */
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
DEV uint32 gcd32(uint32 x, uint32 y)
{
	uint32 q, f;
	if(!y) return(x);
	while(y) {
		q = x/y;	/* Find quotient of current x/y and round toward zero: */
		f = x - q*y;/* Find y' and store in temporary: */
		x = y;		/* Find x', i.e. move the old value of y into the slots for x: */
		y = f;		/* New value of y: */
	}
	return(x);
}

DEV uint64 gcd64(uint64 x, uint64 y)
{
	uint64 q, f;
	if(!y) return(x);
	while(y) {
		q = x/y;	/* Find quotient of current x/y and round toward zero: */
		f = x - q*y;/* Find y' and store in temporary: */
		x = y;		/* Find x', i.e. move the old value of y into the slots for x: */
		y = f;		/* New value of y: */
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
DEV uint32 egcd32(uint32 *x, uint32 *y)
{
	uint32 g = *x, w = *y, q;
	uint32 a = 1, b = 0, u = 0, v = 1;
	/* Sign of these 3 doesn't matter since they're just temporaries: */
	uint32 d, e, f;

	if(*x == *y) {
		printf("ERROR: eGCD of identical arguments x = y = %u is illegal!\n", *x);	ASSERT(HERE, 0,"0");
	} else if((*x == 0) || (*y == 0)) {
		printf("ERROR: eGCD called with zero input: x = %u, y = %u\n", *x, *y);		ASSERT(HERE, 0,"0");
	}

	while(w)
	{
		// Find quotient of current x/y and round toward zero - makes sense to try to take advantage of the fact
		// that most q's are small (~80% of q's < 4), but in practice I've found that adding even simple logic to
		// special-case for q = 0 (e.g. if(g < w) {d = a; e = b; f = g; } else ...) slows things down:
		q = g/w;
//printf("egcd32: w,q = %d, %d, quotient = %d\n",w,g,q);
//printf("a,b,g = %d,%d,%d\n",a,b,g);
//printf("u,v,w = %d,%d,%d\n",u,v,w);
		/* Find (u', v', w') and store in 3 temporaries: */
		d = a - q*u;
		e = b - q*v;
		f = g - q*w;
//printf("d,e,f = %d,%d,%d\n",d,e,f);
		// Find (a', b', g'), i.e. move the old values of (u,v,w) into the slots for (a,b,g),
		// then recover new values of (u, v, w) from the temporaries:
		a = u;	u = d;
		b = v;	v = e;
		g = w;	w = f;
	}
	if(*y < a)	// E.g. inputs 2,2^31-1 gives a = 3221225473 = (int)-1073741823, need to add to modulus (*y) to get proper mod-inv 1073741824.
		*x = *y + a;
	else
		*x = a;
	*y = b;
	return(g);
}

/*********************/
/*
Finds multiplicative inverse of z (mod n).
*/
DEV int modinv32(uint32 z, uint32 n)
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

// Returns binary-form character representation of a uint64 in char_buf.
// Assumes char_buf has at least 65 bytes allocated.
int	convert_uint64_base2_char(char char_buf[], uint64 q)
{
	int i;
	const int hex_chars[2] = {'0','1'};

	for(i=63; i >= 0; i--)
	{
		char_buf[i] = hex_chars[q & 1];
		q >>= 1;
	}
	char_buf[64] = '\0';
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

/*********************** SIMD functionality/cycle-count tests: **********************************/
#ifdef TEST_SIMD

	// Random (digits of Pi) input data sufficient for 64 AVX1024-sized vec_dbl elements of 16 doubles each:
	const char ran[1024] = {
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
	2,1,7,1,2,2,6,8,0,6,6,1,3,0,0,1,9,2,7,8,7,6,6,1,1,1,9,5,9,0,9,2,1,6,4,2,0,1,9,8,9,3,8,0,9,5,2,5,7,2,0,1,0,6,5,4,8,5,8,6,3,2,7,8
	};

  #ifdef USE_AVX1024
	int	test_simd_transpose_16x16()
	{
		ASSERT(HERE,0,"function not yet supported!");
		return 0;
	}
  #endif

  #ifdef USE_AVX512
	int	test_simd_transpose_8x8()
	{
		/*...time-related stuff	*/
		double clock1, clock2;
		double tdiff, t0,t1,t2,t3;
		int i,imax = 100000001, row,col, nerr;	// Use 10^8 loop execs in effort to yield timing on order of 1 sec on target CPUs
			// Add 1 to make loop count odd, thus result of (imax) successive transposes equivalent to a single one
		const int dim = 64;	// #elements in our matrix, allocate 2x this to allow for real/imag side-by-side variant
		vec_dbl *mem  = ALLOC_VEC_DBL(mem, 2*dim+4);	// Add 4 pads to allow for alignment on up-to-128-byte boundary
		vec_dbl *data = ALIGN_VEC_DBL(mem);
		ASSERT(HERE, ((long)data & 0x1f) == 0, "data not 32-byte aligned!");
		// Init the matrix -  Input matrix has rows containing [0-7][8-15]...[56-63]:
		double *dptr = (double *)data;
		for(i = 0; i < dim; i++) { *(dptr+i) = i; }
	//	printf("Input matrix:\n");
		for(i = 0; i < dim; i += 8) {
			row = i>>3;
		//	printf("Row %u: %3.0f %3.0f %3.0f %3.0f %3.0f %3.0f %3.0f %3.0f\n",row,*(dptr+i),*(dptr+i+1),*(dptr+i+2),*(dptr+i+3),*(dptr+i+4),*(dptr+i+5),*(dptr+i+6),*(dptr+i+7));
		}
	//	printf("\n");

		// Do timing loop using 2 fundamentally different methods of effecting the transpose, the 2nd of
		// which mimics the data movement surrounding the dyadic-square and carry steps of our FFT-mul:

		// [1a] Rowwise-load and in-register data shuffles. On KNL: 45 cycles per loop-exec:
		nerr = 0; clock1 = getRealTime();
		for(i = 0; i < imax; i++) {
			__asm__ volatile (\
				"movq		%[__data],%%rax		\n\t"\
				/* Read in the 8 rows of our input matrix: */\
				"vmovaps		0x000(%%rax),%%zmm0		\n\t"\
				"vmovaps		0x040(%%rax),%%zmm1		\n\t"\
				"vmovaps		0x080(%%rax),%%zmm2		\n\t"\
				"vmovaps		0x0c0(%%rax),%%zmm3		\n\t"\
				"vmovaps		0x100(%%rax),%%zmm4		\n\t"\
				"vmovaps		0x140(%%rax),%%zmm5		\n\t"\
				"vmovaps		0x180(%%rax),%%zmm6		\n\t"\
				"vmovaps		0x1c0(%%rax),%%zmm7		\n\t"\
				/* Transpose uses regs0-7 for data, reg8 for temp: */\
				/* [1] First step is a quartet of [UNPCKLPD,UNPCKHPD] pairs to effect transposed 2x2 submatrices - */\
				/* indices in comments at right are [row,col] pairs, i.e. octal version of linear array indices: */
				"vunpcklpd		%%zmm1,%%zmm0,%%zmm8	\n\t"/* zmm8 = 00 10 02 12 04 14 06 16 */\
				"vunpckhpd		%%zmm1,%%zmm0,%%zmm1	\n\t"/* zmm1 = 01 11 03 13 05 15 07 17 */\
				"vunpcklpd		%%zmm3,%%zmm2,%%zmm0	\n\t"/* zmm0 = 20 30 22 32 24 34 26 36 */\
				"vunpckhpd		%%zmm3,%%zmm2,%%zmm3	\n\t"/* zmm3 = 21 31 23 33 25 35 27 37 */\
				"vunpcklpd		%%zmm5,%%zmm4,%%zmm2	\n\t"/* zmm2 = 40 50 42 52 44 54 46 56 */\
				"vunpckhpd		%%zmm5,%%zmm4,%%zmm5	\n\t"/* zmm5 = 41 51 43 53 45 55 47 57 */\
				"vunpcklpd		%%zmm7,%%zmm6,%%zmm4	\n\t"/* zmm4 = 60 70 62 72 64 74 66 76 */\
				"vunpckhpd		%%zmm7,%%zmm6,%%zmm7	\n\t"/* zmm7 = 61 71 63 73 65 75 67 77 */\
			/**** Getting rid of reg-index-nicifying copies here means Outputs not in 0-7 but in 8,1,0,3,2,5,4,7, with 6 now free ****/\
				/* [2] 1st layer of VSHUFF64x2, 2 outputs each with trailing index pairs [0,4],[1,5],[2,6],[3,7]. */\
				/* Note the imm8 values expressed in terms of 2-bit index subfields again read right-to-left */\
				/* (as for the SHUFPS imm8 values in the AVX 8x8 float code) are 221 = (3,1,3,1) and 136 = (2,0,2,0): */\
				"vshuff64x2	$136,%%zmm0,%%zmm8,%%zmm6	\n\t"/* zmm6 = 00 10 04 14 20 30 24 34 */\
				"vshuff64x2	$221,%%zmm0,%%zmm8,%%zmm0	\n\t"/* zmm0 = 02 12 06 16 22 32 26 36 */\
				"vshuff64x2	$136,%%zmm3,%%zmm1,%%zmm8	\n\t"/* zmm8 = 01 11 05 15 21 31 25 35 */\
				"vshuff64x2	$221,%%zmm3,%%zmm1,%%zmm3	\n\t"/* zmm3 = 03 13 07 17 23 33 27 37 */\
				"vshuff64x2	$136,%%zmm4,%%zmm2,%%zmm1	\n\t"/* zmm1 = 40 50 44 54 60 70 64 74 */\
				"vshuff64x2	$221,%%zmm4,%%zmm2,%%zmm4	\n\t"/* zmm4 = 42 52 46 56 62 72 66 76 */\
				"vshuff64x2	$136,%%zmm7,%%zmm5,%%zmm2	\n\t"/* zmm2 = 41 51 45 55 61 71 65 75 */\
				"vshuff64x2	$221,%%zmm7,%%zmm5,%%zmm7	\n\t"/* zmm7 = 43 53 47 57 63 73 67 77 */\
			/**** Getting rid of reg-index-nicifying copies here means Outputs 8,1,2,5 -> 6,8,1,2, with 5 now free ***/\
				/* [3] Last step in 2nd layer of VSHUFF64x2, now combining reg-pairs sharing same trailing index pairs. */\
				/* Output register indices reflect trailing index of data contained therein: */\
				"vshuff64x2	$136,%%zmm1,%%zmm6,%%zmm5	\n\t"/* zmm5 = 00 10 20 30 40 50 60 70 [row 0 of transpose-matrix] */\
				"vshuff64x2	$221,%%zmm1,%%zmm6,%%zmm1	\n\t"/* zmm1 = 04 14 24 34 44 54 64 74 [row 4 of transpose-matrix] */\
				"vshuff64x2	$136,%%zmm2,%%zmm8,%%zmm6	\n\t"/* zmm6 = 01 11 21 31 41 51 61 71 [row 1 of transpose-matrix] */\
				"vshuff64x2	$221,%%zmm2,%%zmm8,%%zmm2	\n\t"/* zmm2 = 05 15 25 35 45 55 65 75 [row 5 of transpose-matrix] */\
				"vshuff64x2	$136,%%zmm4,%%zmm0,%%zmm8	\n\t"/* zmm8 = 02 12 22 32 42 52 62 72 [row 2 of transpose-matrix] */\
				"vshuff64x2	$221,%%zmm4,%%zmm0,%%zmm4	\n\t"/* zmm4 = 06 16 26 36 46 56 66 76 [row 6 of transpose-matrix] */\
				"vshuff64x2	$136,%%zmm7,%%zmm3,%%zmm0	\n\t"/* zmm0 = 03 13 23 33 43 53 63 73 [row 3 of transpose-matrix] */\
				"vshuff64x2	$221,%%zmm7,%%zmm3,%%zmm7	\n\t"/* zmm7 = 07 17 27 37 47 57 67 77 [row 7 of transpose-matrix] */\
			/**** Getting rid of reg-index-nicifying copies here means Outputs 6,8,0,3 -> 5,6,8,0 with 3 now free ***/\
				/* Write original columns back as rows: */\
				"vmovaps		%%zmm5,0x000(%%rax)		\n\t"\
				"vmovaps		%%zmm6,0x040(%%rax)		\n\t"\
				"vmovaps		%%zmm8,0x080(%%rax)		\n\t"\
				"vmovaps		%%zmm0,0x0c0(%%rax)		\n\t"\
				"vmovaps		%%zmm1,0x100(%%rax)		\n\t"\
				"vmovaps		%%zmm2,0x140(%%rax)		\n\t"\
				"vmovaps		%%zmm4,0x180(%%rax)		\n\t"\
				"vmovaps		%%zmm7,0x1c0(%%rax)		\n\t"\
				:						// outputs: none
				: [__data] "m" (data)	// All inputs from memory addresses here
				: "cc","memory","rax","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8"	// Clobbered registers - use xmm form for compatibility with older versions of clang/gcc
			);
		}
		clock2 = getRealTime();
		tdiff = (double)(clock2 - clock1);
		printf("Method [1a]: Time for %u 8x8 doubles-transposes using in-register shuffles =%s\n",imax, get_time_str(tdiff));
		// Check the result:
	//	printf("Output matrix:\n");
		for(i = 0; i < dim; i += 8) {
			row = i>>3;
		//	printf("Row %u: %3.0f %3.0f %3.0f %3.0f %3.0f %3.0f %3.0f %3.0f\n",row,*(dptr+i),*(dptr+i+1),*(dptr+i+2),*(dptr+i+3),*(dptr+i+4),*(dptr+i+5),*(dptr+i+6),*(dptr+i+7));
			// Expected (transposed-matrix) datum = row + 4*col
			t0 = row; t1 = row+8; t2 = row+16; t3 = row+24;
			nerr += (t0 != *(dptr+i+0)) + (t1 != *(dptr+i+1)) + (t2 != *(dptr+i+2)) + (t3 != *(dptr+i+3));
			t0 += 32; t1 += 32; t2 += 32; t3 += 32;
			nerr += (t0 != *(dptr+i+4)) + (t1 != *(dptr+i+5)) + (t2 != *(dptr+i+6)) + (t3 != *(dptr+i+7));
		}
		if(nerr) printf("Outputs incorrect! #mismatches = %u\n",nerr);

		// [1b] Same as [1a] but with a few reg-copies to make for a nicer indexing pattern. On KNL: 48 cycles per loop-exec:
		for(i = 0; i < dim; i++) { *(dptr+i) = i; }	// Re-init the matrix to be untransposed
		nerr = 0; clock1 = getRealTime();
		for(i = 0; i < imax; i++) {
			__asm__ volatile (\
				"movq		%[__data],%%rax		\n\t"\
				/* Read in the 8 rows of our input matrix: */\
				"vmovaps		0x000(%%rax),%%zmm0		\n\t"\
				"vmovaps		0x040(%%rax),%%zmm1		\n\t"\
				"vmovaps		0x080(%%rax),%%zmm2		\n\t"\
				"vmovaps		0x0c0(%%rax),%%zmm3		\n\t"\
				"vmovaps		0x100(%%rax),%%zmm4		\n\t"\
				"vmovaps		0x140(%%rax),%%zmm5		\n\t"\
				"vmovaps		0x180(%%rax),%%zmm6		\n\t"\
				"vmovaps		0x1c0(%%rax),%%zmm7		\n\t"\
				/* [1] First step is a quartet of [UNPCKLPD,UNPCKHPD] pairs to effect transposed 2x2 submatrices - VUNPCK latency 4-7, rthru = 2: */\
				"vunpckhpd		 %%zmm1,%%zmm0,%%zmm8									\n\t"/* zmm0 = 00 10 02 12 04 14 06 16 [after reg-copy on next line] */\
				"vunpcklpd		 %%zmm1,%%zmm0,%%zmm0 	\n\t	vmovaps	%%zmm8,%%zmm1 	\n\t"/* zmm1 = 01 11 03 13 05 15 07 17 */\
				"vunpckhpd		 %%zmm3,%%zmm2,%%zmm8									\n\t"/* zmm2 = 20 30 22 32 24 34 26 36 */\
				"vunpcklpd		 %%zmm3,%%zmm2,%%zmm2 	\n\t	vmovaps	%%zmm8,%%zmm3 	\n\t"/* zmm3 = 21 31 23 33 25 35 27 37 */\
				"vunpckhpd		 %%zmm5,%%zmm4,%%zmm8									\n\t"/* zmm4 = 40 50 42 52 44 54 46 56 */\
				"vunpcklpd		 %%zmm5,%%zmm4,%%zmm4 	\n\t	vmovaps	%%zmm8,%%zmm5	\n\t"/* zmm5 = 41 51 43 53 45 55 47 57 */\
				"vunpckhpd		 %%zmm7,%%zmm6,%%zmm8									\n\t"/* zmm6 = 60 70 62 72 64 74 66 76 */\
				"vunpcklpd		 %%zmm7,%%zmm6,%%zmm6	\n\t	vmovaps	%%zmm8,%%zmm7	\n\t"/* zmm7 = 61 71 63 73 65 75 67 77 */\
				/* [2] 1st layer of VSHUFF64x2, 2 outputs each with trailing index pairs [0,4],[1,5],[2,6],[3,7] - VSHUFF64x2 latency 4-7, rthru = 2: */\
				"vshuff64x2	$136,%%zmm2,%%zmm0,%%zmm8									\n\t"/* zmm0 = 00 10 04 14 20 30 24 34 */\
				"vshuff64x2	$221,%%zmm2,%%zmm0,%%zmm2 	\n\t	vmovaps	%%zmm8,%%zmm0 	\n\t"/* zmm2 = 02 12 06 16 22 32 26 36 */\
				"vshuff64x2	$136,%%zmm3,%%zmm1,%%zmm8									\n\t"/* zmm1 = 01 11 05 15 21 31 25 35 */\
				"vshuff64x2	$221,%%zmm3,%%zmm1,%%zmm3 	\n\t	vmovaps	%%zmm8,%%zmm1 	\n\t"/* zmm3 = 03 13 07 17 23 33 27 37 */\
				"vshuff64x2	$136,%%zmm6,%%zmm4,%%zmm8									\n\t"/* zmm4 = 40 50 44 54 60 70 64 74 */\
				"vshuff64x2	$221,%%zmm6,%%zmm4,%%zmm6	\n\t	vmovaps	%%zmm8,%%zmm4 	\n\t"/* zmm6 = 42 52 46 56 62 72 66 76 */\
				"vshuff64x2	$136,%%zmm7,%%zmm5,%%zmm8									\n\t"/* zmm5 = 41 51 45 55 61 71 65 75 */\
				"vshuff64x2	$221,%%zmm7,%%zmm5,%%zmm7	\n\t	vmovaps	%%zmm8,%%zmm5	\n\t"/* zmm7 = 43 53 47 57 63 73 67 77 */\
				/* [3] Last step in 2nd layer of VSHUFF64x2, now combining reg-pairs sharing same trailing index pairs: */\
				"vshuff64x2	$136,%%zmm4,%%zmm0,%%zmm8									\n\t"/* zmm0 = 00 10 20 30 40 50 60 70 */\
				"vshuff64x2	$221,%%zmm4,%%zmm0,%%zmm4 	\n\t	vmovaps	%%zmm8,%%zmm0 	\n\t"/* zmm4 = 04 14 24 34 44 54 64 74 */\
				"vshuff64x2	$136,%%zmm5,%%zmm1,%%zmm8									\n\t"/* zmm1 = 01 11 21 31 41 51 61 71 */\
				"vshuff64x2	$221,%%zmm5,%%zmm1,%%zmm5	\n\t	vmovaps	%%zmm8,%%zmm1 	\n\t"/* zmm5 = 05 15 25 35 45 55 65 75 */\
				"vshuff64x2	$136,%%zmm6,%%zmm2,%%zmm8									\n\t"/* zmm2 = 02 12 22 32 42 52 62 72 */\
				"vshuff64x2	$221,%%zmm6,%%zmm2,%%zmm6	\n\t	vmovaps	%%zmm8,%%zmm2 	\n\t"/* zmm6 = 06 16 26 36 46 56 66 76 */\
				"vshuff64x2	$136,%%zmm7,%%zmm3,%%zmm8									\n\t"/* zmm3 = 03 13 23 33 43 53 63 73 */\
				"vshuff64x2	$221,%%zmm7,%%zmm3,%%zmm7	\n\t	vmovaps	%%zmm8,%%zmm3 	\n\t"/* zmm7 = 07 17 27 37 47 57 67 77 */\
				/* Write original columns back as rows: */\
				"vmovaps		%%zmm0,0x000(%%rax)		\n\t"\
				"vmovaps		%%zmm1,0x040(%%rax)		\n\t"\
				"vmovaps		%%zmm2,0x080(%%rax)		\n\t"\
				"vmovaps		%%zmm3,0x0c0(%%rax)		\n\t"\
				"vmovaps		%%zmm4,0x100(%%rax)		\n\t"\
				"vmovaps		%%zmm5,0x140(%%rax)		\n\t"\
				"vmovaps		%%zmm6,0x180(%%rax)		\n\t"\
				"vmovaps		%%zmm7,0x1c0(%%rax)		\n\t"\
				:						// outputs: none
				: [__data] "m" (data)	// All inputs from memory addresses here
				: "cc","memory","rax","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8"	// Clobbered registers - use xmm form for compatibility with older versions of clang/gcc
			);
		}
		clock2 = getRealTime();
		tdiff = (double)(clock2 - clock1);
		printf("Method [1b]: Time for %u 8x8 doubles-transposes using in-register shuffles =%s\n",imax, get_time_str(tdiff));
		// Check the result:
	//	printf("Output matrix:\n");
		for(i = 0; i < dim; i += 8) {
			row = i>>3;
		//	printf("Row %u: %3.0f %3.0f %3.0f %3.0f %3.0f %3.0f %3.0f %3.0f\n",row,*(dptr+i),*(dptr+i+1),*(dptr+i+2),*(dptr+i+3),*(dptr+i+4),*(dptr+i+5),*(dptr+i+6),*(dptr+i+7));
			// Expected (transposed-matrix) datum = row + 4*col
			t0 = row; t1 = row+8; t2 = row+16; t3 = row+24;
			nerr += (t0 != *(dptr+i+0)) + (t1 != *(dptr+i+1)) + (t2 != *(dptr+i+2)) + (t3 != *(dptr+i+3));
			t0 += 32; t1 += 32; t2 += 32; t3 += 32;
			nerr += (t0 != *(dptr+i+4)) + (t1 != *(dptr+i+5)) + (t2 != *(dptr+i+6)) + (t3 != *(dptr+i+7));
		}
		if(nerr) printf("Outputs incorrect! #mismatches = %u\n",nerr);

		/* [1c; Apr 2018] Variant which GeorgeW says saves a few cycles on his Skylake-X:
			"The idea is to use vbroadcastf64x4 to do the 256-bit shuffles. This is more
			load uops, but the masking can be done on either port 0 or port 5which is better
			than the vshuff64x2 it replaces which can only be done on port 5."
		George says this saves a few cycles over the above on Skylake-X, but in my test on KNL it's even slower
		than gather-based variant [2] - seems Agner Fog's 5-cycle latency for vbroadcastf64x4 on KNL is way low.
		KNL: 68 cycles, over 1.5x the 44 cycles for my best 24-shuffle version [1b], and 20% slower than [2]'s 56.
		With shuffle-passes 2 and 3 deleted (i.e. just the load, vbroadcastf64x4 and store steps), get 62 cycles,
		which means vbroadcastf64x4 is roughly as slow as vgatherdpd[!], which differs markedly from the 5-cycle-latency,
		2-per-cycle throughput listed for vbroadcastf64x4-on-KNL in Agner Fog's x86 instruction tables compilation.
		That, coupled with the middle-2-col-pairs-swapped-versus-transpose nature of the output makes it a no-go for me.
		On KNL: 62 cycles per loop-exec:
		*/
		for(i = 0; i < dim; i++) { *(dptr+i) = i; }	// Re-init the matrix to be untransposed
		nerr = 0; clock1 = getRealTime();
		for(i = 0; i < imax; i++) {
			__asm__ volatile (\
				/* Init opmasks k2 and k3 [need only low byte of each]: */\
				"movl	$0xf0,%%eax	\n\t	kmovw	%%eax,%%k1	\n\t"\
				"movl	$0x0f,%%eax	\n\t	kmovw	%%eax,%%k2	\n\t"\
				"movq		%[__data],%%rax		\n\t"\
				/* Read in the 8 rows of our input matrix: */\
				"vmovaps		0x000(%%rax),%%zmm2		\n\t"\
				"vmovaps		0x040(%%rax),%%zmm4		\n\t"\
				"vmovaps		0x080(%%rax),%%zmm5		\n\t"\
				"vmovaps		0x0c0(%%rax),%%zmm3		\n\t"\
				"vmovaps		0x100(%%rax),%%zmm6		\n\t"\
				"vmovaps		0x140(%%rax),%%zmm8		\n\t"\
				"vmovaps		0x180(%%rax),%%zmm0		\n\t"\
				"vmovaps		0x1c0(%%rax),%%zmm7		\n\t"\
				/* [1] Interleave lo/hi halves of rows 0-3 with those of rows 4-7, respectively: */\
				"vbroadcastf64x4	0x100(%%rax),%%zmm2%{%%k1%}	\n\t"/* zmm2 = 00 01 02 03 40 41 42 43 */\
				"vbroadcastf64x4	0x140(%%rax),%%zmm4%{%%k1%}	\n\t"/* zmm4 = 10 11 12 13 50 51 52 53 */\
				"vbroadcastf64x4	0x180(%%rax),%%zmm5%{%%k1%}	\n\t"/* zmm5 = 20 21 22 23 60 61 62 63 */\
				"vbroadcastf64x4	0x1c0(%%rax),%%zmm3%{%%k1%}	\n\t"/* zmm3 = 30 31 32 33 70 71 72 73 */\
				"vbroadcastf64x4	0x020(%%rax),%%zmm6%{%%k2%}	\n\t"/* zmm6 = 04 05 06 07 44 45 46 47 */\
				"vbroadcastf64x4	0x060(%%rax),%%zmm8%{%%k2%}	\n\t"/* zmm8 = 14 15 16 17 54 55 56 57 */\
				"vbroadcastf64x4	0x0a0(%%rax),%%zmm0%{%%k2%}	\n\t"/* zmm0 = 24 25 26 27 64 65 66 67 */\
				"vbroadcastf64x4	0x0e0(%%rax),%%zmm7%{%%k2%}	\n\t"/* zmm7 = 34 35 36 37 74 75 76 77 */\
			/* Now a simple quartet of 4x4 transposes on the resulting four 4x4 submatrices suffices to give the
			desired 8x8 transpose, BUT! - the 4x4 AVX transpose code uses a set of vshufpd (just as our step [2] below)
			followed by a step based on vperm2f128, and there is no 512-bit version of the latter instruction. */\
				/* [2] Use 8 VSHUFPD to effect transposes of the eight 2x2 submatrices: */\
				"vshufpd	$0x00,%%zmm4,%%zmm2,%%zmm1	\n\t"/* zmm1 = 00 10 02 12 40 50 42 52 */\
				"vshufpd	$0xff,%%zmm4,%%zmm2,%%zmm4	\n\t"/* zmm4 = 01 11 03 13 41 51 43 53 */\
				"vshufpd	$0x00,%%zmm3,%%zmm5,%%zmm2	\n\t"/* zmm2 = 20 30 22 32 60 70 62 72 */\
				"vshufpd	$0xff,%%zmm3,%%zmm5,%%zmm3	\n\t"/* zmm3 = 21 31 23 33 61 71 63 73 */\
				"vshufpd	$0x00,%%zmm8,%%zmm6,%%zmm5	\n\t"/* zmm5 = 04 14 06 16 44 54 46 56 */\
				"vshufpd	$0xff,%%zmm8,%%zmm6,%%zmm8	\n\t"/* zmm8 = 05 15 07 17 45 55 47 57 */\
				"vshufpd	$0x00,%%zmm7,%%zmm0,%%zmm6	\n\t"/* zmm6 = 24 34 26 36 64 74 66 76 */\
				"vshufpd	$0xff,%%zmm7,%%zmm0,%%zmm7	\n\t"/* zmm7 = 25 35 27 37 65 75 67 77 */\
				/* [3] Last step is layer of VSHUFF64x2, now combining reg-pairs sharing same trailing index pairs. */\
				/* Note the imm8 values expressed in terms of 2-bit index subfields again read right-to-left (as for the SHUFPS imm7 */\
				/* values in the AVX 8x8 float code) are 0x88 = (2,0,2,0), 0xdd = (3,1,3,1), 0x22 = (0,2,0,2) and 0x77 = (1,3,1,3): */\
				/* Output register indices reflect trailing index of data contained therein: */\
													/****** Output col-pairs [4,5],[2,3] swapped! ******/\
				"vshuff64x2	$0x88,%%zmm2,%%zmm1,%%zmm0	\n\t"/* zmm0 = 00 10 40 50 20 30 60 70 [row 0 of transpose-matrix] */\
				"vshuff64x2	$0xdd,%%zmm2,%%zmm1,%%zmm2	\n\t"/* zmm2 = 02 12 42 52 22 32 62 72 [row 2 of transpose-matrix] */\
				"vshuff64x2	$0x88,%%zmm3,%%zmm4,%%zmm1	\n\t"/* zmm1 = 01 11 41 51 21 31 61 71 [row 1 of transpose-matrix] */\
				"vshuff64x2	$0xdd,%%zmm3,%%zmm4,%%zmm3	\n\t"/* zmm3 = 03 13 43 53 23 33 63 73 [row 3 of transpose-matrix] */\
				"vshuff64x2	$0x88,%%zmm6,%%zmm5,%%zmm4	\n\t"/* zmm4 = 04 14 44 54 24 34 64 74 [row 4 of transpose-matrix] */\
				"vshuff64x2	$0xdd,%%zmm6,%%zmm5,%%zmm6	\n\t"/* zmm6 = 06 16 46 56 26 36 66 76 [row 6 of transpose-matrix] */\
				"vshuff64x2	$0x88,%%zmm7,%%zmm8,%%zmm5	\n\t"/* zmm5 = 05 15 45 55 25 35 65 75 [row 5 of transpose-matrix] */\
				"vshuff64x2	$0xdd,%%zmm7,%%zmm8,%%zmm7	\n\t"/* zmm7 = 07 17 47 57 27 37 67 77 [row 7 of transpose-matrix] */\
				/* Write original columns back as rows: */\
				"vmovaps		%%zmm0,0x000(%%rax)		\n\t"\
				"vmovaps		%%zmm1,0x040(%%rax)		\n\t"\
				"vmovaps		%%zmm2,0x080(%%rax)		\n\t"\
				"vmovaps		%%zmm3,0x0c0(%%rax)		\n\t"\
				"vmovaps		%%zmm4,0x100(%%rax)		\n\t"\
				"vmovaps		%%zmm5,0x140(%%rax)		\n\t"\
				"vmovaps		%%zmm6,0x180(%%rax)		\n\t"\
				"vmovaps		%%zmm7,0x1c0(%%rax)		\n\t"\
				:						// outputs: none
				: [__data] "m" (data)	// All inputs from memory addresses here
				: "cc","memory","rax","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8"	// Clobbered registers - use xmm form for compatibility with older versions of clang/gcc
			);
		}
		clock2 = getRealTime();
		tdiff = (double)(clock2 - clock1);
		// We know this version does not produce a true transpose, so just check timing:
		printf("Method [1c]: Time for %u MIDDLE-COL-PAIR-SWAPPED 8x8 doubles-transposes using in-register shuffles =%s\n",imax, get_time_str(tdiff));

		// [1d] Skylake-X-oriented variant from George. On KNL: 40 cycles per loop-exec, ~15% faster than [1a]:
		for(i = 0; i < dim; i++) { *(dptr+i) = i; }	// Re-init the matrix to be untransposed
		nerr = 0; clock1 = getRealTime();
		for(i = 0; i < imax; i++) {
			__asm__ volatile (\
				"movl	$0b00110011,%%eax	\n\t"/* Constant for vblendmpd instructions goes into mask-reg k1 */\
				"kmovw	%%eax,%%k1			\n\t"\
				/* Init vector index-consts needed by vpermt2pd instructions - if regs were at a premium,
				could also init just prior to [3] and use zmm6,7 to hold index-consts: */\
				"movq	$0x0c040e0608000a02,%%rax	\n\t"/* zmm30 = 8+4 0+4 8+6 0+6 8+0 0+0 8+2 0+2 [msw at left] */\
				"movq	$0x0d050f0709010b03,%%rbx	\n\t"/* zmm31 = 8+5 0+5 8+7 0+7 8+1 0+1 8+3 0+3 */\
					"vmovq		%%rax,%%xmm0 		\n\t"\
					"vmovq		%%rbx,%%xmm1 		\n\t"\
					"vpmovzxbq	%%xmm0,%%zmm30		\n\t"\
					"vpmovzxbq	%%xmm1,%%zmm31		\n\t"\
				"movq		%[__data],%%rax		\n\t"\
				/* Read in the 8 rows of our input matrix: */\
				"vmovaps		0x000(%%rax),%%zmm0		\n\t"\
				"vmovaps		0x040(%%rax),%%zmm1		\n\t"\
				"vmovaps		0x080(%%rax),%%zmm2		\n\t"\
				"vmovaps		0x0c0(%%rax),%%zmm3		\n\t"\
				"vmovaps		0x100(%%rax),%%zmm4		\n\t"\
				"vmovaps		0x140(%%rax),%%zmm5		\n\t"\
				"vmovaps		0x180(%%rax),%%zmm6		\n\t"\
				"vmovaps		0x1c0(%%rax),%%zmm7		\n\t"\
				/* [1] Shuffle the 4-aparts - note the different patterning of the first and second output quartet: */\
				"vshuff64x2	$0b01000100, %%zmm4,	%%zmm0,	%%zmm8 	\n\t"/* 00 01 02 03 40 41 42 43 */\
				"vshuff64x2	$0b11101110, %%zmm4,	%%zmm0,	%%zmm4 	\n\t"/* 04 05 06 07 44 45 46 47 */\
				"vshuff64x2	$0b01000100, %%zmm5,	%%zmm1,	%%zmm9	\n\t"/* 10 11 12 13 50 51 52 53 */\
				"vshuff64x2	$0b11101110, %%zmm5,	%%zmm1,	%%zmm5	\n\t"/* 14 15 16 17 54 55 56 57 */\
				"vshuff64x2	$0b00010001, %%zmm6,	%%zmm2,	%%zmm10	\n\t"/* 22 23 20 21 62 63 60 61 */\
				"vshuff64x2	$0b10111011, %%zmm6,	%%zmm2,	%%zmm6	\n\t"/* 26 27 24 25 66 67 64 65 */\
				"vshuff64x2	$0b00010001, %%zmm7,	%%zmm3,	%%zmm11	\n\t"/* 32 33 30 31 72 73 70 71 */\
				"vshuff64x2	$0b10111011, %%zmm7,	%%zmm3,	%%zmm7	\n\t"/* 36 37 34 35 76 77 74 75 *//* data in 4-11; 0-3 free */\
				/* [2] Blend in the 2-aparts */\
				"vblendmpd	%%zmm8 ,	%%zmm10,	%%zmm0%{%%k1%}	\n\t"/* 00 01 20 21 40 41 60 61 */\
				"vblendmpd	%%zmm10,	%%zmm8 ,	%%zmm8%{%%k1%}	\n\t"/* 22 23 02 03 62 63 42 43 */\
				"vblendmpd	%%zmm4 ,	%%zmm6 ,	%%zmm1%{%%k1%}	\n\t"/* 04 05 24 25 44 45 64 65 */\
				"vblendmpd	%%zmm6 ,	%%zmm4 ,	%%zmm4%{%%k1%}	\n\t"/* 26 27 06 07 66 67 46 47 */\
				"vblendmpd	%%zmm9 ,	%%zmm11,	%%zmm2%{%%k1%}	\n\t"/* 10 11 30 31 50 51 70 71 */\
				"vblendmpd	%%zmm11,	%%zmm9 ,	%%zmm9%{%%k1%}	\n\t"/* 32 33 12 13 72 73 52 53 */\
				"vblendmpd	%%zmm5 ,	%%zmm7 ,	%%zmm3%{%%k1%}	\n\t"/* 14 15 34 35 54 55 74 75 */\
				"vblendmpd	%%zmm7 ,	%%zmm5 ,	%%zmm5%{%%k1%}	\n\t"/* 36 37 16 17 76 77 56 57 *//* data in 0-5,8-9; 6-7,10-11 free */\
				/* [3] Shuffle or permute in the 1-aparts */\
				"vshufpd	$0b00000000,%%zmm2,	%%zmm0,%%zmm10 	\n\t"/* 00 10 20 30 40 50 60 70 */\
				"vshufpd	$0b11111111,%%zmm2,	%%zmm0,%%zmm11 	\n\t"/* 01 11 21 31 41 51 61 71 */\
				"vmovapd	%%zmm8,%%zmm2	\n\t"\
				"vpermt2pd				%%zmm9,	%%zmm30,%%zmm2 	\n\t"/* 02 12 22 32 42 52 62 72 */\
				"vpermt2pd				%%zmm9,	%%zmm31,%%zmm8	\n\t"/* 03 13 23 33 43 53 63 73 */\
				"vshufpd	$0b00000000,%%zmm3,	%%zmm1,%%zmm0 	\n\t"/* 04 14 24 34 44 54 64 74 */\
				"vshufpd	$0b11111111,%%zmm3,	%%zmm1,%%zmm1 	\n\t"/* 05 15 25 35 45 55 65 75 */\
				"vmovapd	%%zmm4,%%zmm3	\n\t"\
				"vpermt2pd				%%zmm5,	%%zmm30,%%zmm3 	\n\t"/* 06 16 26 36 46 56 66 76 */\
				"vpermt2pd				%%zmm5,	%%zmm31,%%zmm4	\n\t"/* 07 17 27 37 47 57 67 77 */\
				/* Write original columns back as rows: */\
				"vmovaps		%%zmm10,0x000(%%rax)		\n\t"\
				"vmovaps		%%zmm11,0x040(%%rax)		\n\t"\
				"vmovaps		%%zmm2 ,0x080(%%rax)		\n\t"\
				"vmovaps		%%zmm8 ,0x0c0(%%rax)		\n\t"\
				"vmovaps		%%zmm0 ,0x100(%%rax)		\n\t"\
				"vmovaps		%%zmm1 ,0x140(%%rax)		\n\t"\
				"vmovaps		%%zmm3 ,0x180(%%rax)		\n\t"\
				"vmovaps		%%zmm4, 0x1c0(%%rax)		\n\t"\
				:						// outputs: none
				: [__data] "m" (data)	// All inputs from memory addresses here
				: "cc","memory","rax","rbx","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11", "xmm30","xmm31"	// Clobbered registers - use xmm form for compatibility with older versions of clang/gcc
			);
		}
		clock2 = getRealTime();
		tdiff = (double)(clock2 - clock1);
		printf("Method [1d]: Time for %u 8x8 doubles-transposes using [vshuff64x2,vblendmpd,vshufpd,vpermt2pd] =%s\n",imax, get_time_str(tdiff));
		// Check the result:
	//	printf("Output matrix:\n");
		for(i = 0; i < dim; i += 8) {
			row = i>>3;
		//	printf("Row %u: %3.0f %3.0f %3.0f %3.0f %3.0f %3.0f %3.0f %3.0f\n",row,*(dptr+i),*(dptr+i+1),*(dptr+i+2),*(dptr+i+3),*(dptr+i+4),*(dptr+i+5),*(dptr+i+6),*(dptr+i+7));
			// Expected (transposed-matrix) datum = row + 4*col
			t0 = row; t1 = row+8; t2 = row+16; t3 = row+24;
			nerr += (t0 != *(dptr+i+0)) + (t1 != *(dptr+i+1)) + (t2 != *(dptr+i+2)) + (t3 != *(dptr+i+3));
			t0 += 32; t1 += 32; t2 += 32; t3 += 32;
			nerr += (t0 != *(dptr+i+4)) + (t1 != *(dptr+i+5)) + (t2 != *(dptr+i+6)) + (t3 != *(dptr+i+7));
		}
		if(nerr) printf("Outputs incorrect! #mismatches = %u\n",nerr);

		// [2a] Columnwise-load-and-rowwise-writeback using AVX512 gather-load functionality. On KNL: 56 cycles per loop-exec.
		/* Compare latency/thruput/ports for shuffle and gather-based versions, using Agner Fog's KNL tables:
		[1] vunpcklpd, vshuff64x2 both have 4-7 cycle latency, one can start every 2 cycles [3,1 on Skylake-X].
								Thus 24 such in sequence with no wait-stalls ==> ~50 cycles, close to what I measure.
		[2] vgatherdpd has 7-cycle latency, no data re. thruput, thus ~60 cycles per loop, again close to that observed.
		Both [1] and [2] use port 5, but on KNL the 'empty' cycles between shuffle-op issues can be used to issue gathers,
		which is what the side-by-ide matrix-pair [2b] variant tests.
		*/
		for(i = 0; i < dim; i++) { *(dptr+i) = i; }	// Re-init the matrix to be untransposed
		nerr = 0; clock1 = getRealTime();
		for(i = 0; i < imax; i++) {	// Nov 2016: 4.3 sec for 10^8 loops @1.3GHz ==> ~7 cycles per gather-load
			__asm__ volatile (\
				"movq		%[__data],%%rax		\n\t"\
				/* Auxiliary register data needed for columnwise loads: */\
			"movq	$0x1c1814100c080400,%%rbx	\n\t"/* 64-bit register w/byte offsets 0x[00,04,08,0c,10,14,18,1c], bytes numbered left-to-right */\
				"vmovq		%%rbx,%%xmm8 		\n\t"/* Copy byte pattern to low qword (64 bits) of ymm8 [NB: avx-512 only supports MOVQ to/from 128-bit vector regs] */\
				"vpmovzxbd	%%xmm8,%%ymm8		\n\t"/* vector-index offsets: ymm8 = 0x[00,04,08,0c,10,14,18,1c] in 32-bit form in low 8 dwords */\
				"vpslld	$4,%%ymm8,%%ymm8		\n\t"/* The above bytewise offsets need scale *16 to get the needed ones - would include but
												e.g. 0x1C<<4 overflows 1 byte), but x86 ISA only permits scale factors 1,2,4,8, so <<= 4 here. */\
			/* Mask-reg zmm9 = 11...11 - this is stupidly zeroed each time we do gather-load, so need to reinit: */\
			"movl	$-1,%%ebx	\n\t"/* Init opmask k1 (Only need the low byte) */\
			/* Gather instruction sets mask-reg = 0, so must re-init opmask prior to each invocation */
			"kmovw	%%ebx,%%k1	\n\t	vgatherdpd 0x00(%%rax,%%ymm8),%%zmm0%{%%k1%}	\n\t"/* Col 0 */\
			"kmovw	%%ebx,%%k1	\n\t	vgatherdpd 0x08(%%rax,%%ymm8),%%zmm1%{%%k1%}	\n\t"/* Col 1 */\
			"kmovw	%%ebx,%%k1	\n\t	vgatherdpd 0x10(%%rax,%%ymm8),%%zmm2%{%%k1%}	\n\t"/* Col 2 */\
			"kmovw	%%ebx,%%k1	\n\t	vgatherdpd 0x18(%%rax,%%ymm8),%%zmm3%{%%k1%}	\n\t"/* Col 3 */\
			"kmovw	%%ebx,%%k1	\n\t	vgatherdpd 0x20(%%rax,%%ymm8),%%zmm4%{%%k1%}	\n\t"/* Col 4 */\
			"kmovw	%%ebx,%%k1	\n\t	vgatherdpd 0x28(%%rax,%%ymm8),%%zmm5%{%%k1%}	\n\t"/* Col 5 */\
			"kmovw	%%ebx,%%k1	\n\t	vgatherdpd 0x30(%%rax,%%ymm8),%%zmm6%{%%k1%}	\n\t"/* Col 6 */\
			"kmovw	%%ebx,%%k1	\n\t	vgatherdpd 0x38(%%rax,%%ymm8),%%zmm7%{%%k1%}	\n\t"/* Col 7 */\
				/* Write original columns back as rows: */\
				"vmovaps	%%zmm0,0x000(%%rax)	\n\t"\
				"vmovaps	%%zmm1,0x040(%%rax)	\n\t"\
				"vmovaps	%%zmm2,0x080(%%rax)	\n\t"\
				"vmovaps	%%zmm3,0x0c0(%%rax)	\n\t"\
				"vmovaps	%%zmm4,0x100(%%rax)	\n\t"\
				"vmovaps	%%zmm5,0x140(%%rax)	\n\t"\
				"vmovaps	%%zmm6,0x180(%%rax)	\n\t"\
				"vmovaps	%%zmm7,0x1c0(%%rax)	\n\t"\
				:						// outputs: none
				: [__data] "m" (data)	// All inputs from memory addresses here
				: "cc","memory","rax","rbx","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8"	// Clobbered registers - use xmm form for compatibility with older versions of clang/gcc
			);
		}
		clock2 = getRealTime();
		tdiff = (double)(clock2 - clock1);
		printf("Method [2a]: Time for %u 8x8 doubles-transposes using gather-loads =%s\n",imax, get_time_str(tdiff));
	//	printf("Output matrix:\n");
		for(i = 0; i < dim; i += 8) {
			row = i>>3;
		//	printf("Row %u: %3.0f %3.0f %3.0f %3.0f %3.0f %3.0f %3.0f %3.0f\n",row,*(dptr+i),*(dptr+i+1),*(dptr+i+2),*(dptr+i+3),*(dptr+i+4),*(dptr+i+5),*(dptr+i+6),*(dptr+i+7));
			// Expected (transposed-matrix) datum = row + 4*col
			t0 = row; t1 = row+8; t2 = row+16; t3 = row+24;
			nerr += (t0 != *(dptr+i+0)) + (t1 != *(dptr+i+1)) + (t2 != *(dptr+i+2)) + (t3 != *(dptr+i+3));
			t0 += 32; t1 += 32; t2 += 32; t3 += 32;
			nerr += (t0 != *(dptr+i+4)) + (t1 != *(dptr+i+5)) + (t2 != *(dptr+i+6)) + (t3 != *(dptr+i+7));
		}
		if(nerr) printf("Outputs incorrect! #mismatches = %u\n",nerr);

		// [2b] Hybrid shuffle/gather side-by-side matrix-pair [2b] variant of [2a]. On KNL this is faster
		// than [1a]+[2a] separately, i.e. we do save some cycles by interleaving the 2 types of transposes,
		// but still needs 1.93x the cycles of the best shuffle-based variant, thus no faster than 2 side-by-side
		// shuffle-transposes, *and* the comparison looks likely to come out even more unfavorably on Skylake-X:
		for(i = 0; i < 2*dim; i++) { *(dptr+i) = i; }	// Re-init the matrix to be untransposed
		nerr = 0; clock1 = getRealTime();
		for(i = 0; i < imax; i++) {
			__asm__ volatile (\
				"movq		%[__data],%%rax		\n\t"\
				/* Lcol: read in 8 rows of matrix use by shuffle-transpose while setting up for gather-loads in rcol: */\
				"vmovaps		0x200(%%rax),%%zmm0		\n\t	movq	$0x1c1814100c080400,%%rbx	\n\t"\
				"vmovaps		0x240(%%rax),%%zmm1		\n\t	vmovq		%%rbx,%%xmm9 		\n\t"\
				"vmovaps		0x280(%%rax),%%zmm2		\n\t	vpmovzxbd	%%xmm9,%%ymm9		\n\t"/* vpmovzxbd CAN ONLY USE XMM0-15! */\
				"vmovaps		0x2c0(%%rax),%%zmm3		\n\t"\
				"vmovaps		0x300(%%rax),%%zmm4		\n\t"\
				"vmovaps		0x340(%%rax),%%zmm5		\n\t"\
				"vmovaps		0x380(%%rax),%%zmm6		\n\t	vpslld	$4,%%ymm9,%%ymm9		\n\t"\
				"vmovaps		0x3c0(%%rax),%%zmm7		\n\t	movl	$-1,%%ebx				\n\t"\
				/* Transpose uses regs0-7 for data, reg8 for temp: */\
									/* Gather instruction sets mask-reg = 0, so must re-init opmask prior to each invocation */
									"kmovw	%%ebx,%%k1	\n\t	vgatherdpd 0x00(%%rax,%%ymm9),%%zmm10%{%%k1%}	\n\t"/* Col 0 */\
				"vunpcklpd		%%zmm1,%%zmm0,%%zmm8	\n\t"\
				"vunpckhpd		%%zmm1,%%zmm0,%%zmm1	\n\t"\
				"vunpcklpd		%%zmm3,%%zmm2,%%zmm0	\n\t"\
									"kmovw	%%ebx,%%k1	\n\t	vgatherdpd 0x08(%%rax,%%ymm9),%%zmm11%{%%k1%}	\n\t"/* Col 1 */\
				"vunpckhpd		%%zmm3,%%zmm2,%%zmm3	\n\t"\
				"vunpcklpd		%%zmm5,%%zmm4,%%zmm2	\n\t"\
				"vunpckhpd		%%zmm5,%%zmm4,%%zmm5	\n\t"\
									"kmovw	%%ebx,%%k1	\n\t	vgatherdpd 0x10(%%rax,%%ymm9),%%zmm12%{%%k1%}	\n\t"/* Col 2 */\
				"vunpcklpd		%%zmm7,%%zmm6,%%zmm4	\n\t"\
				"vunpckhpd		%%zmm7,%%zmm6,%%zmm7	\n\t"\
				"vshuff64x2	$136,%%zmm0,%%zmm8,%%zmm6	\n\t"\
									"kmovw	%%ebx,%%k1	\n\t	vgatherdpd 0x18(%%rax,%%ymm9),%%zmm13%{%%k1%}	\n\t"/* Col 3 */\
				"vshuff64x2	$221,%%zmm0,%%zmm8,%%zmm0	\n\t"\
				"vshuff64x2	$136,%%zmm3,%%zmm1,%%zmm8	\n\t"\
				"vshuff64x2	$221,%%zmm3,%%zmm1,%%zmm3	\n\t"\
									"kmovw	%%ebx,%%k1	\n\t	vgatherdpd 0x20(%%rax,%%ymm9),%%zmm14%{%%k1%}	\n\t"/* Col 4 */\
				"vshuff64x2	$136,%%zmm4,%%zmm2,%%zmm1	\n\t"\
				"vshuff64x2	$221,%%zmm4,%%zmm2,%%zmm4	\n\t"\
				"vshuff64x2	$136,%%zmm7,%%zmm5,%%zmm2	\n\t"\
									"kmovw	%%ebx,%%k1	\n\t	vgatherdpd 0x28(%%rax,%%ymm9),%%zmm15%{%%k1%}	\n\t"/* Col 5 */\
				"vshuff64x2	$221,%%zmm7,%%zmm5,%%zmm7	\n\t"\
				"vshuff64x2	$136,%%zmm1,%%zmm6,%%zmm5	\n\t"/* [output row 0] */\
				"vshuff64x2	$221,%%zmm1,%%zmm6,%%zmm1	\n\t"/* [output row 4] */\
									"kmovw	%%ebx,%%k1	\n\t	vgatherdpd 0x30(%%rax,%%ymm9),%%zmm16%{%%k1%}	\n\t"/* Col 6 */\
				"vshuff64x2	$136,%%zmm2,%%zmm8,%%zmm6	\n\t"/* [output row 1] */\
				"vshuff64x2	$221,%%zmm2,%%zmm8,%%zmm2	\n\t"/* [output row 5] */\
				"vshuff64x2	$136,%%zmm4,%%zmm0,%%zmm8	\n\t"/* [output row 2] */\
									"kmovw	%%ebx,%%k1	\n\t	vgatherdpd 0x38(%%rax,%%ymm9),%%zmm17%{%%k1%}	\n\t"/* Col 7 */\
				"vshuff64x2	$221,%%zmm4,%%zmm0,%%zmm4	\n\t"/* [output row 6] */\
				"vshuff64x2	$136,%%zmm7,%%zmm3,%%zmm0	\n\t"/* [output row 3] */\
				"vshuff64x2	$221,%%zmm7,%%zmm3,%%zmm7	\n\t"/* [output row 7] */\
				/* Write original columns back as rows: */\
				"vmovaps		%%zmm5,0x200(%%rax)		\n\t"\
				"vmovaps		%%zmm6,0x240(%%rax)		\n\t"\
				"vmovaps		%%zmm8,0x280(%%rax)		\n\t"\
				"vmovaps		%%zmm0,0x2c0(%%rax)		\n\t"\
				"vmovaps		%%zmm1,0x300(%%rax)		\n\t"\
				"vmovaps		%%zmm2,0x340(%%rax)		\n\t"\
				"vmovaps		%%zmm4,0x380(%%rax)		\n\t"\
				"vmovaps		%%zmm7,0x3c0(%%rax)		\n\t"\
																"vmovaps	%%zmm10,0x000(%%rax)	\n\t"\
																"vmovaps	%%zmm11,0x040(%%rax)	\n\t"\
																"vmovaps	%%zmm12,0x080(%%rax)	\n\t"\
																"vmovaps	%%zmm13,0x0c0(%%rax)	\n\t"\
																"vmovaps	%%zmm14,0x100(%%rax)	\n\t"\
																"vmovaps	%%zmm15,0x140(%%rax)	\n\t"\
																"vmovaps	%%zmm16,0x180(%%rax)	\n\t"\
																"vmovaps	%%zmm17,0x1c0(%%rax)	\n\t"\
				:						// outputs: none
				: [__data] "m" (data)	// All inputs from memory addresses here
				: "cc","memory","rax","rbx","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15","xmm16","xmm17"	// Clobbered registers - use xmm form for compatibility with older versions of clang/gcc
			);
		}
		clock2 = getRealTime();
		tdiff = (double)(clock2 - clock1);
		printf("Method [2b]: Time for 2 x %u 8x8 doubles-transposes using gather-loads =%s\n",imax, get_time_str(tdiff));
		// Check the result:
	//	printf("Output matrix 1:\n");
		for(i = 0; i < dim; i += 8) {
			row = i>>3;
		//	printf("Row %u: %3.0f %3.0f %3.0f %3.0f %3.0f %3.0f %3.0f %3.0f\n",row,*(dptr+i),*(dptr+i+1),*(dptr+i+2),*(dptr+i+3),*(dptr+i+4),*(dptr+i+5),*(dptr+i+6),*(dptr+i+7));
			// Expected (transposed-matrix) datum = row + 4*col
			t0 = row; t1 = row+8; t2 = row+16; t3 = row+24;
			nerr += (t0 != *(dptr+i+0)) + (t1 != *(dptr+i+1)) + (t2 != *(dptr+i+2)) + (t3 != *(dptr+i+3));
			t0 += 32; t1 += 32; t2 += 32; t3 += 32;
			nerr += (t0 != *(dptr+i+4)) + (t1 != *(dptr+i+5)) + (t2 != *(dptr+i+6)) + (t3 != *(dptr+i+7));
		}
	//	printf("Output matrix 2:\n");
		for(i = dim; i < 2*dim; i += 8) {
			row = i>>3;
		//	printf("Row %2u: %3.0f %3.0f %3.0f %3.0f %3.0f %3.0f %3.0f %3.0f\n",row,*(dptr+i),*(dptr+i+1),*(dptr+i+2),*(dptr+i+3),*(dptr+i+4),*(dptr+i+5),*(dptr+i+6),*(dptr+i+7));
			// Expected (transposed-matrix) datum = row + 4*col
			t0 = row+56; t1 = row+64; t2 = row+72; t3 = row+80;
			nerr += (t0 != *(dptr+i+0)) + (t1 != *(dptr+i+1)) + (t2 != *(dptr+i+2)) + (t3 != *(dptr+i+3));
			t0 += 32; t1 += 32; t2 += 32; t3 += 32;
			nerr += (t0 != *(dptr+i+4)) + (t1 != *(dptr+i+5)) + (t2 != *(dptr+i+6)) + (t3 != *(dptr+i+7));
		}
		if(nerr) printf("Outputs incorrect! #mismatches = %u\n",nerr);

		// [2c] Side-by-side matrix-pair variant of [1d]. On KNL ?
		for(i = 0; i < 2*dim; i++) { *(dptr+i) = i; }	// Re-init the matrix to be untransposed
		nerr = 0; clock1 = getRealTime();
		for(i = 0; i < imax; i++) {
			__asm__ volatile (\
				"movl	$0b00110011,%%eax	\n\t"/* Constant for vblendmpd instructions goes into mask-reg k1 */\
				"kmovw	%%eax,%%k1			\n\t"\
				/* Init vector index-consts needed by vpermt2pd instructions - if regs were at a premium,
				could also init just prior to [3] and use zmm6,7 to hold index-consts: */\
				"movq	$0x0c040e0608000a02,%%rax	\n\t"/* zmm30 = 8+4 0+4 8+6 0+6 8+0 0+0 8+2 0+2 [msw at left] */\
				"movq	$0x0d050f0709010b03,%%rbx	\n\t"/* zmm31 = 8+5 0+5 8+7 0+7 8+1 0+1 8+3 0+3 */\
					"vmovq		%%rax,%%xmm0 		\n\t"\
					"vmovq		%%rbx,%%xmm1 		\n\t"\
					"vpmovzxbq	%%xmm0,%%zmm30		\n\t"\
					"vpmovzxbq	%%xmm1,%%zmm31		\n\t"\
				"movq		%[__data],%%rax		\n\t"\
				/* Read in the 8 rows of our input matrix: */\
				"vmovaps		0x000(%%rax),%%zmm0					\n\t	vmovaps		0x200(%%rax),%%zmm12		\n\t"\
				"vmovaps		0x040(%%rax),%%zmm1					\n\t	vmovaps		0x240(%%rax),%%zmm13		\n\t"\
				"vmovaps		0x080(%%rax),%%zmm2					\n\t	vmovaps		0x280(%%rax),%%zmm14		\n\t"\
				"vmovaps		0x0c0(%%rax),%%zmm3					\n\t	vmovaps		0x2c0(%%rax),%%zmm15		\n\t"\
				"vmovaps		0x100(%%rax),%%zmm4					\n\t	vmovaps		0x300(%%rax),%%zmm16		\n\t"\
				"vmovaps		0x140(%%rax),%%zmm5					\n\t	vmovaps		0x340(%%rax),%%zmm17		\n\t"\
				"vmovaps		0x180(%%rax),%%zmm6					\n\t	vmovaps		0x380(%%rax),%%zmm18		\n\t"\
				"vmovaps		0x1c0(%%rax),%%zmm7					\n\t	vmovaps		0x3c0(%%rax),%%zmm19		\n\t"\
				/* [1] Shuffle the 4-aparts - note the different patterning of the first and second output quartet: */\
				"vshuff64x2	$0b01000100, %%zmm4,	%%zmm0,	%%zmm8 	\n\t	vshuff64x2	$0b01000100, %%zmm16,%%zmm12,	%%zmm20	\n\t"\
				"vshuff64x2	$0b11101110, %%zmm4,	%%zmm0,	%%zmm4 	\n\t	vshuff64x2	$0b11101110, %%zmm16,%%zmm12,	%%zmm16	\n\t"\
				"vshuff64x2	$0b01000100, %%zmm5,	%%zmm1,	%%zmm9	\n\t	vshuff64x2	$0b01000100, %%zmm17,%%zmm13,	%%zmm21	\n\t"\
				"vshuff64x2	$0b11101110, %%zmm5,	%%zmm1,	%%zmm5	\n\t	vshuff64x2	$0b11101110, %%zmm17,%%zmm13,	%%zmm17	\n\t"\
				"vshuff64x2	$0b00010001, %%zmm6,	%%zmm2,	%%zmm10	\n\t	vshuff64x2	$0b00010001, %%zmm18,%%zmm14,	%%zmm22	\n\t"\
				"vshuff64x2	$0b10111011, %%zmm6,	%%zmm2,	%%zmm6	\n\t	vshuff64x2	$0b10111011, %%zmm18,%%zmm14,	%%zmm18	\n\t"\
				"vshuff64x2	$0b00010001, %%zmm7,	%%zmm3,	%%zmm11	\n\t	vshuff64x2	$0b00010001, %%zmm19,%%zmm15,	%%zmm23	\n\t"\
				"vshuff64x2	$0b10111011, %%zmm7,	%%zmm3,	%%zmm7	\n\t	vshuff64x2	$0b10111011, %%zmm19,%%zmm15,	%%zmm19	\n\t"\
				/* [2] Blend in the 2-aparts */\
				"vblendmpd	%%zmm8 ,	%%zmm10,	%%zmm0%{%%k1%}	\n\t	vblendmpd	%%zmm20,	%%zmm22,	%%zmm12%{%%k1%}	\n\t"\
				"vblendmpd	%%zmm10,	%%zmm8 ,	%%zmm8%{%%k1%}	\n\t	vblendmpd	%%zmm22,	%%zmm20,	%%zmm20%{%%k1%}	\n\t"\
				"vblendmpd	%%zmm4 ,	%%zmm6 ,	%%zmm1%{%%k1%}	\n\t	vblendmpd	%%zmm16,	%%zmm18,	%%zmm13%{%%k1%}	\n\t"\
				"vblendmpd	%%zmm6 ,	%%zmm4 ,	%%zmm4%{%%k1%}	\n\t	vblendmpd	%%zmm18,	%%zmm16,	%%zmm16%{%%k1%}	\n\t"\
				"vblendmpd	%%zmm9 ,	%%zmm11,	%%zmm2%{%%k1%}	\n\t	vblendmpd	%%zmm21,	%%zmm23,	%%zmm14%{%%k1%}	\n\t"\
				"vblendmpd	%%zmm11,	%%zmm9 ,	%%zmm9%{%%k1%}	\n\t	vblendmpd	%%zmm23,	%%zmm21,	%%zmm21%{%%k1%}	\n\t"\
				"vblendmpd	%%zmm5 ,	%%zmm7 ,	%%zmm3%{%%k1%}	\n\t	vblendmpd	%%zmm17,	%%zmm19,	%%zmm15%{%%k1%}	\n\t"\
				"vblendmpd	%%zmm7 ,	%%zmm5 ,	%%zmm5%{%%k1%}	\n\t	vblendmpd	%%zmm19,	%%zmm17,	%%zmm17%{%%k1%}	\n\t"\
				/* [3] Shuffle or permute in the 1-aparts */\
				"vshufpd	$0b00000000,%%zmm2,		%%zmm0,%%zmm10 	\n\t	vshufpd	$0b00000000,%%zmm14,	%%zmm12,%%zmm22	\n\t"\
				"vshufpd	$0b11111111,%%zmm2,		%%zmm0,%%zmm11 	\n\t	vshufpd	$0b11111111,%%zmm14,	%%zmm12,%%zmm23	\n\t"\
				"vmovapd	%%zmm8,%%zmm2							\n\t	vmovapd	%%zmm20,%%zmm14	\n\t"\
				"vpermt2pd				%%zmm9,		%%zmm30,%%zmm2 	\n\t	vpermt2pd				%%zmm21,	%%zmm30,%%zmm14	\n\t"\
				"vpermt2pd				%%zmm9,		%%zmm31,%%zmm8	\n\t	vpermt2pd				%%zmm21,	%%zmm31,%%zmm20	\n\t"\
				"vshufpd	$0b00000000,%%zmm3,		%%zmm1,%%zmm0 	\n\t	vshufpd	$0b00000000,%%zmm15,	%%zmm13,%%zmm12	\n\t"\
				"vshufpd	$0b11111111,%%zmm3,		%%zmm1,%%zmm1 	\n\t	vshufpd	$0b11111111,%%zmm15,	%%zmm13,%%zmm13	\n\t"\
				"vmovapd	%%zmm4,%%zmm3							\n\t	vmovapd	%%zmm16,%%zmm15	\n\t"\
				"vpermt2pd				%%zmm5,		%%zmm30,%%zmm3 	\n\t	vpermt2pd				%%zmm17,	%%zmm30,%%zmm15	\n\t"\
				"vpermt2pd				%%zmm5,		%%zmm31,%%zmm4	\n\t	vpermt2pd				%%zmm17,	%%zmm31,%%zmm16	\n\t"\
				/* Write original columns back as rows: */\
				"vmovaps		%%zmm10,0x000(%%rax)				\n\t	vmovaps		%%zmm22,0x200(%%rax)		\n\t"\
				"vmovaps		%%zmm11,0x040(%%rax)				\n\t	vmovaps		%%zmm23,0x240(%%rax)		\n\t"\
				"vmovaps		%%zmm2 ,0x080(%%rax)				\n\t	vmovaps		%%zmm14,0x280(%%rax)		\n\t"\
				"vmovaps		%%zmm8 ,0x0c0(%%rax)				\n\t	vmovaps		%%zmm20,0x2c0(%%rax)		\n\t"\
				"vmovaps		%%zmm0 ,0x100(%%rax)				\n\t	vmovaps		%%zmm12,0x300(%%rax)		\n\t"\
				"vmovaps		%%zmm1 ,0x140(%%rax)				\n\t	vmovaps		%%zmm13,0x340(%%rax)		\n\t"\
				"vmovaps		%%zmm3 ,0x180(%%rax)				\n\t	vmovaps		%%zmm15,0x380(%%rax)		\n\t"\
				"vmovaps		%%zmm4, 0x1c0(%%rax)				\n\t	vmovaps		%%zmm16,0x3c0(%%rax)		\n\t"\
				:						// outputs: none
				: [__data] "m" (data)	// All inputs from memory addresses here
				: "cc","memory","rax","rbx","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15","xmm16","xmm17","xmm18","xmm19","xmm20","xmm21","xmm22","xmm23", "xmm30","xmm31"	// Clobbered registers - use xmm form for compatibility with older versions of clang/gcc
			);
		}
		clock2 = getRealTime();
		tdiff = (double)(clock2 - clock1);
		printf("Method [2c]: Time for 2 x %u 8x8 doubles-transposes using side-by-side impl of algo [1d] =%s\n",imax, get_time_str(tdiff));
		// Check the result:
	//	printf("Output matrix 1:\n");
		for(i = 0; i < dim; i += 8) {
			row = i>>3;
		//	printf("Row %u: %3.0f %3.0f %3.0f %3.0f %3.0f %3.0f %3.0f %3.0f\n",row,*(dptr+i),*(dptr+i+1),*(dptr+i+2),*(dptr+i+3),*(dptr+i+4),*(dptr+i+5),*(dptr+i+6),*(dptr+i+7));
			// Expected (transposed-matrix) datum = row + 4*col
			t0 = row; t1 = row+8; t2 = row+16; t3 = row+24;
			nerr += (t0 != *(dptr+i+0)) + (t1 != *(dptr+i+1)) + (t2 != *(dptr+i+2)) + (t3 != *(dptr+i+3));
			t0 += 32; t1 += 32; t2 += 32; t3 += 32;
			nerr += (t0 != *(dptr+i+4)) + (t1 != *(dptr+i+5)) + (t2 != *(dptr+i+6)) + (t3 != *(dptr+i+7));
		}
	//	printf("Output matrix 2:\n");
		for(i = dim; i < 2*dim; i += 8) {
			row = i>>3;
		//	printf("Row %u: %3.0f %3.0f %3.0f %3.0f %3.0f %3.0f %3.0f %3.0f\n",row,*(dptr+i),*(dptr+i+1),*(dptr+i+2),*(dptr+i+3),*(dptr+i+4),*(dptr+i+5),*(dptr+i+6),*(dptr+i+7));
			// Expected (transposed-matrix) datum = row + 4*col
			t0 = row+56; t1 = row+64; t2 = row+72; t3 = row+80;
			nerr += (t0 != *(dptr+i+0)) + (t1 != *(dptr+i+1)) + (t2 != *(dptr+i+2)) + (t3 != *(dptr+i+3));
			t0 += 32; t1 += 32; t2 += 32; t3 += 32;
			nerr += (t0 != *(dptr+i+4)) + (t1 != *(dptr+i+5)) + (t2 != *(dptr+i+6)) + (t3 != *(dptr+i+7));
		}
		if(nerr) printf("Outputs incorrect! #mismatches = %u\n",nerr);

		return nerr;
	}
  #endif

  #ifdef USE_AVX	//********* Since AVX2 & AVX can run on much of the same hardware, lump them together
					// and differentiate as needed within this outer preprocessor conditional - e.g. the
					// test_simd_transpose_4x4() function uses both, if available, for comparative timings.
	// 4x4 refers to linear memory treated as a 4x4 matrix-of-doubles:
	int	test_simd_transpose_4x4()
	{
		/*...time-related stuff	*/
		double clock1, clock2;
		double tdiff, t0,t1,t2,t3;
		int i,imax = 100000001, row,col, nerr = 0;	// Use 10^8 loop execs in effort to yield timing on order of 1 sec on target CPUs
			// Add 1 to make loop count odd, thus result of (imax) successive transposes equivalent to a single one
		const int dim = 16;		// #elements in our matrix
		vec_dbl *mem  = ALLOC_VEC_DBL(mem, dim+4);	// Add 4 pads to allow for alignment on up-to-128-byte boundary
		vec_dbl *data = ALIGN_VEC_DBL(mem);
		ASSERT(HERE, ((long)data & 0x1f) == 0, "data not 32-byte aligned!");
		// Init the matrix -  Input matrix has rows:
		double *dptr = (double *)data;	//  0, 1, 2, 3
		for(i = 0; i < dim; i++) {		//  4, 5, 6, 7
			*(dptr+i) = i;				//  8, 9,10,11
		}								// 12,13,14,15

		// Do timing loop using 2 fundamentally different methods of effecting the transpose,
		// the 1st of which comes in 2 variants dubbed [1a] and [1b]:

		// [1a] Rowwise-load and in-register data shuffles. On KNL: 23 cycles per loop-exec:
		clock1 = getRealTime();
		for(i = 0; i < imax; i++) {
			__asm__ volatile (\
				"movq		%[__data],%%rax		\n\t"\
				"vmovaps	     (%%rax),%%ymm2			\n\t"\
				"vmovaps	0x020(%%rax),%%ymm1			\n\t"\
				"vshufpd	$15,%%ymm1,%%ymm2,%%ymm3	\n\t"\
				"vshufpd	$0 ,%%ymm1,%%ymm2,%%ymm2	\n\t"\
				"vmovaps	0x040(%%rax),%%ymm4			\n\t"\
				"vmovaps	0x060(%%rax),%%ymm1			\n\t"\
				"vshufpd	$15,%%ymm1,%%ymm4,%%ymm0	\n\t"\
				"vshufpd	$0 ,%%ymm1,%%ymm4,%%ymm4	\n\t"\
				"vperm2f128 $32,%%ymm0,%%ymm3,%%ymm1	\n\t"/* Row 1 */\
				"vperm2f128 $49,%%ymm0,%%ymm3,%%ymm3	\n\t"/* Row 3 */\
				"vperm2f128 $32,%%ymm4,%%ymm2,%%ymm0	\n\t"/* Row 0 */\
				"vperm2f128 $49,%%ymm4,%%ymm2,%%ymm2	\n\t"/* Row 2 */\
				/* Write original columns back as rows: */\
				"vmovaps	%%ymm0,0x00(%%rax)	\n\t"\
				"vmovaps	%%ymm1,0x20(%%rax)	\n\t"\
				"vmovaps	%%ymm2,0x40(%%rax)	\n\t"\
				"vmovaps	%%ymm3,0x60(%%rax)	\n\t"\
				:						// outputs: none
				: [__data] "m" (data)	// All inputs from memory addresses here
				: "cc","memory","rax","xmm0","xmm1","xmm2","xmm3","xmm4"	// Clobbered registers - use xmm form for compatibility with older versions of clang/gcc
			);
		}
		clock2 = getRealTime();
		tdiff = (double)(clock2 - clock1);
		printf("Method [1a]: Time for %u 4x4 matrix-of-doubles transposes =%s\n",imax, get_time_str(tdiff));
		// Check the result:
	//	printf("Output matrix:\n");
		for(i = 0; i < dim; i += 4) {
			row = i>>2;
		//	printf("Row %u: %3.0f %3.0f %3.0f %3.0f\n",row,*(dptr+i),*(dptr+i+1),*(dptr+i+2),*(dptr+i+3));
			t0 = row; t1 = row+4; t2 = row+8; t3 = row+12;	// Expected (transposed-matrix) datum = row + 4*col
			nerr += (t0 != *(dptr+i)) + (t1 != *(dptr+i+1)) + (t2 != *(dptr+i+2)) + (t3 != *(dptr+i+3));
		}
	//	printf("#mismatches = %u\n",nerr);

		// [1b] Rowwise-load and in-register data shuffles, using a different shuffle sequence. On KNL: 24 cycles per loop-exec:
		for(i = 0; i < dim; i++) { *(dptr+i) = i; }	// Re-init the matrix to be untransposed
		clock1 = getRealTime();
		for(i = 0; i < imax; i++) {
			__asm__ volatile (\
				"movq		%[__data],%%rax		\n\t"\
				"vmovaps	    (%%rax),%%xmm0				\n\t"/* r0.lo = 0,1,-,- */\
				"vmovaps	0x20(%%rax),%%xmm4				\n\t"/* r1.lo = 4,5,-,- */\
				"vinsertf128 $1,0x40(%%rax),%%ymm0,%%ymm0	\n\t"/* r0|r2.lo = 0,1,8,9 */\
				"vinsertf128 $1,0x60(%%rax),%%ymm4,%%ymm4	\n\t"/* r1|r3.lo = 4,5,c,d */\
				"vshufpd	$15,%%ymm4,%%ymm0,%%ymm1		\n\t"/* Row 1 = 1,5,9,d */\
				"vshufpd	$0 ,%%ymm4,%%ymm0,%%ymm0		\n\t"/* Row 0 = 0,4,8,c */\
				"vmovaps	0x10(%%rax),%%xmm2				\n\t"/* r0.hi = 2,3,-,- */\
				"vmovaps	0x30(%%rax),%%xmm4				\n\t"/* r1.hi = 6,7,-,- */\
				"vinsertf128 $1,0x50(%%rax),%%ymm2,%%ymm2	\n\t"/* r0|r2.hi = 2,3,a,b */\
				"vinsertf128 $1,0x70(%%rax),%%ymm4,%%ymm4	\n\t"/* r1|r3.hi = 6,7,e,f */\
				"vshufpd	$15,%%ymm4,%%ymm2,%%ymm3		\n\t"/* Row 3 = 3,7,b,f */\
				"vshufpd	$0 ,%%ymm4,%%ymm2,%%ymm2		\n\t"/* Row 2 = 2,6,a,e */\
				/* Write original columns back as rows: */\
				"vmovaps	%%ymm0,0x00(%%rax)	\n\t"\
				"vmovaps	%%ymm1,0x20(%%rax)	\n\t"\
				"vmovaps	%%ymm2,0x40(%%rax)	\n\t"\
				"vmovaps	%%ymm3,0x60(%%rax)	\n\t"\
				:						// outputs: none
				: [__data] "m" (data)	// All inputs from memory addresses here
				: "cc","memory","rax","xmm0","xmm1","xmm2","xmm3","xmm4"	// Clobbered registers - use xmm form for compatibility with older versions of clang/gcc
			);
		}
		clock2 = getRealTime();
		tdiff = (double)(clock2 - clock1);
		printf("Method [1b]: Time for %u 4x4 matrix-of-doubles transposes =%s\n",imax, get_time_str(tdiff));
		// Check the result:
	//	printf("Output matrix:\n");
		for(i = 0; i < dim; i += 4) {
			row = i>>2;
		//	printf("Row %u: %3.0f %3.0f %3.0f %3.0f\n",row,*(dptr+i),*(dptr+i+1),*(dptr+i+2),*(dptr+i+3));
			t0 = row; t1 = row+4; t2 = row+8; t3 = row+12;	// Expected (transposed-matrix) datum = row + 4*col
			nerr += (t0 != *(dptr+i)) + (t1 != *(dptr+i+1)) + (t2 != *(dptr+i+2)) + (t3 != *(dptr+i+3));
		}
	//	printf("#mismatches = %u\n",nerr);

		// [2] Columnwise-load-and-rowwise-writeback using AVX2 gather-load functionality. On KNL: 46 cycles per loop-exec:
	  #ifdef USE_AVX2
	   #ifdef GCC_5PLUS	// gcc 4.x may not support the needed AVX2 instructions (while still being fine for for the FMA
						// instructions used for the FFT), so require an added compile-time define to enable loop [2]
		for(i = 0; i < dim; i++) { *(dptr+i) = i; }	// Re-init the matrix to be untransposed
		clock1 = getRealTime();
		for(i = 0; i < imax; i++) {
			__asm__ volatile (\
				"movq		%[__data],%%rax		\n\t"\
				/* Auxiliary register data needed for columnwise loads: */\
				"movl	$0x60402000,%%ebx		\n\t"/* 32-bit register w/byte offsets [0x00,0x20,0x40,0x60], bytes numbered left-to-right */\
				"vmovd		%%ebx,%%xmm4 		\n\t"/* Copy byte pattern to low dword (32 bits) of ymm4 [NB: avx-512 only supports MOVZX to/from mem-address or 128-bit vector regs] */\
				"vpmovzxbd	%%xmm4,%%ymm4		\n\t"/* vector-index offsets: ymm4 = [0x00,0x20,0x40,0x60] in 32-bit form in low 4 dwords */\
			/* Mask-reg ymm5 = 11...11 - this is stupidly zeroed each time we do gather-load, so need to reinit: */\
			"vpcmpeqd	%%ymm5,%%ymm5,%%ymm5	\n\t	vgatherdpd %%ymm5,0x00(%%rax,%%xmm4),%%ymm0	\n\t"/* Col 0 */\
			"vpcmpeqd	%%ymm5,%%ymm5,%%ymm5	\n\t	vgatherdpd %%ymm5,0x08(%%rax,%%xmm4),%%ymm1	\n\t"/* Col 1 */\
			"vpcmpeqd	%%ymm5,%%ymm5,%%ymm5	\n\t	vgatherdpd %%ymm5,0x10(%%rax,%%xmm4),%%ymm2	\n\t"/* Col 2 */\
			"vpcmpeqd	%%ymm5,%%ymm5,%%ymm5	\n\t	vgatherdpd %%ymm5,0x18(%%rax,%%xmm4),%%ymm3	\n\t"/* Col 3 */\
				/* Write original columns back as rows: */\
				"vmovaps	%%ymm0,0x00(%%rax)	\n\t"\
				"vmovaps	%%ymm1,0x20(%%rax)	\n\t"\
				"vmovaps	%%ymm2,0x40(%%rax)	\n\t"\
				"vmovaps	%%ymm3,0x60(%%rax)	\n\t"\
				:						// outputs: none
				: [__data] "m" (data)	// All inputs from memory addresses here
				: "cc","memory","rax","rbx","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5"	// Clobbered registers - use xmm form for compatibility with older versions of clang/gcc
			);
		}
		clock2 = getRealTime();
		tdiff = (double)(clock2 - clock1);
		printf("Method [2]: Time for %u 4x4 matrix-of-doubles transposes =%s\n",imax, get_time_str(tdiff));
		// Check the result:
	//	printf("Output matrix:\n");
		for(i = 0; i < dim; i += 4) {
			row = i>>2;
		//	printf("Row %u: %3.0f %3.0f %3.0f %3.0f\n",row,*(dptr+i),*(dptr+i+1),*(dptr+i+2),*(dptr+i+3));
			t0 = row; t1 = row+4; t2 = row+8; t3 = row+12;	// Expected (transposed-matrix) datum = row + 4*col
			nerr += (t0 != *(dptr+i)) + (t1 != *(dptr+i+1)) + (t2 != *(dptr+i+2)) + (t3 != *(dptr+i+3));
		}
	//	printf("#mismatches = %u\n",nerr);
	   #endif
	  #endif
		return nerr;
	}
  #endif	// USE_AVX ?

  #ifdef USE_SSE2
	// Timing loop for 128-bit SIMD radix-4 DFT macro:
	int	test_radix4_dft()
	{
		const char func[] = "test_radix4_dft";
	#ifdef USE_AVX	//No AVX support for macros in this function
		return 0;
	#else
		/*...time-related stuff	*/
		double clock1, clock2;
		double tdiff, t0,t1,t2,t3;
		int i,j,dim,imax = 10000000, nerr = 0;	// Expect radix-4 DFT to need ~1/10th the cycles of radix-32, so use 10^7 loop execs
		int p1,p2,p3,p4;
		static vec_dbl *sc_arr = 0x0, *sc_ptr;
		double *add0,*add1,*add2,*add3;	/* Addresses into array sections */
		// Each row of 16 double data corr. to expected outputs from one of up to two vector-complex 4-DFTs being done:
		const double ref1[] = {-5.,40.,163.,139.,106.,45.,38.,-4.,31.,-36.,-115.,-83.,-120.,-45.,-70.,-48.,
								-70.,15.,230.,158.,28.,55.,-14.,-70.,56., 3.,-116.,-102.,-6.,-61.,-68.,30.};
		const double ref2[] = {22.,20.,20.,18.,-18.,-13.,30.,-1.,-26.,-4.,-18.,-12.,-8.,-146.,-28.,72.,
								13.,15.,31.,16.,-17.,27.,45.,28.,6.,-25.,-24.,15.,-6.,-1.,48.,-57.};
		vec_dbl *c_tmp,*s_tmp, *cc0,*two, *r0,*r1,*r2,*r3;
		// Alloc 8 vector-complex elts (16 vec_dbl) per input/output block rather than 4, so can also test two radix-4 DFTs done side-by-side:
		sc_arr = ALLOC_VEC_DBL(sc_arr, 0x42);	if(!sc_arr){ sprintf(cbuf, "FATAL: unable to allocate sc_arr in %s.\n",func); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		sc_ptr = ALIGN_VEC_DBL(sc_arr);
		ASSERT(HERE, ((long)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");
		add0 = sc_ptr;
		add1 = sc_ptr+0x2;
		add2 = sc_ptr+0x4;
		add3 = sc_ptr+0x6;
		r0 = sc_ptr + 0x10;
		r1 = r0 + 0x2;
		r2 = r0 + 0x4;
		r3 = r0 + 0x6;
		cc0 = r0 + 0x10;	// Alloc 8 vector-complex elts rather than 4, so can also test two radix-4 DFTs done side-by-side
		two = r0 + 0x20;	// Similarly alloc 2 sets of 8 vector-complex twiddles
		VEC_DBL_INIT(two,2.0);

		// Do these timing-test DFTs in-place, i.e. p1 = #doubles in a pair of vec_dbl:
		p1 = RE_IM_STRIDE << 1;
		p2 = p1 + p1;
		p3 = p2 + p1;
		p4 = p2 + p2;

		// Twiddles for the purpose of this timing-test can be anything, just set to random digits in [0,9].
		// First twiddle for each of the two 4-DFT datasets = unity, but since macros don't actually do that cmul,
		// only init the non-unity twiddle-triplet needed by each of the up-to-2-independent DFTs we do:
		c_tmp = cc0; s_tmp = c_tmp+1;	/* c0,s0 */
		for(i = 0; i < 6; i+=2, c_tmp+=2, s_tmp+=2) {	// Remaining 3 vector-complex twiddles for each 4-DFT are nontrivial
			VEC_DBL_INIT(c_tmp  , ran[i  ]);	VEC_DBL_INIT(s_tmp  , ran[i+1]);
			VEC_DBL_INIT(c_tmp+6, ran[i+8]);	VEC_DBL_INIT(s_tmp+6, ran[i+9]);
/*
			// Restructure twiddle-muls to use cotangent-scheme:
			ASSERT(HERE, ran[i+1] != 0.0 && ran[i+9] != 0.0,"Need to modify test-twiddles to avoid div-by-0!");
			VEC_DBL_INIT(c_tmp  , ran[i  ]/(double)ran[i+1]);	VEC_DBL_INIT(s_tmp  , ran[i+1]);
			VEC_DBL_INIT(c_tmp+8, ran[i+8]/(double)ran[i+9]);	VEC_DBL_INIT(s_tmp+8, ran[i+9]);
*/
		}
		// Set inputs != 0 to prevent timings being thrown off by any 0-operand arithmetic-shortcuts the CPU may do:
		double *dptr = (double *)r0;
		dim = 8*RE_IM_STRIDE;	// 4 vector-complex data
		// Copy quasirandom digits-of-Pi-data into our vec_dbl inputs:
		for(j = 0; j < dim; j++) { *(add0+j) = ran[j]; }

	// 5 May 2016: 10^7-loop timings, 1-threaded on my 2GHz Core2Duo:
	//																comments
	// no-DFT timing (keeping just 1 address-load of macro)	0.246	49.2 cycles (!! ... that's a lot for loop-control and data-copy)
	// initial timing										0.473	94.6 - 49.2 = 45.4 cycles
	// leaq for add0 + 8*p[1,2,3]							0.462
	// use xmm8,9 to save 6 implied-loads in mulpd			0.454
	// use xmm8,9 to save 2 implied-loads in addpd			0.444
	// replace four add-doublings by 2*x in final butterfly	0.429	36.6 cycles ... Already nearly 20% faster!
	// move loads of xmm8,9 imm before use to ease 2-column	0.452	that is a step back, revert
	// elim 4 redundant loads in p1,3 computation			0.429	no speedup
	// move 1st 2 loads for next combo into end of pvs one	0.440	slower, oddly
	// eliminate 1st twiddle-mul                            0.438	again slower ... now *that* is bizarre
	// Added spill/reload of first output to cut regs 8,9	0.429	nice - now can do 2 such DFTs side-by-side
	// Restructure twiddle-muls to use cotangent-scheme		0.445	slower
	//*** 7 Jul 2017: *** out-of-place DFT allows 1-time init, i.e.
	// no more need to subtract initial-data-copy timing:	0.232	46 cycles - note this is the 'initial timing' 8-register macro above.
	// ARM Neon code timing (1.5 GHz Odroid C2):			0.612	93 cycles, almost exactly 2c the cycle count of SSE2 on Core2
		// Timing loop #1:
		clock1 = getRealTime();
		for(i = 0; i < imax; i++) {
			/* 4-DFT-with-3-complex-twiddles macro. SIMD opcount:
			x86_64 SSE2: 41 MEM (19 load[1 via mem-op in addpd], 14 store, 8 reg-copy), 22 ADDPD, 16 MULPD
			ARM v8 Neon: 11 MEM (7 load-pair, 4 store-pair), 16 FADD, 12 FMUl/FMA
			*/
		#ifdef USE_ARM_V8_SIMD
			__asm__ volatile (\
				"ldr	x0,%[__add0]		\n\t"\
				"ldr	w1,%[__p1]			\n\t"\
				"ldr	w2,%[__p2]			\n\t"\
				"ldr	w3,%[__p3]			\n\t"\
				"ldr	x4,%[__cc0]			\n\t"\
				"ldr	x5,%[__r0]			\n\t"\
				"add	x1, x0,x1,lsl #3	\n\t"\
				"add	x2, x0,x2,lsl #3	\n\t"\
				"add	x3, x0,x3,lsl #3	\n\t"\
				/* SSE2_RADIX_04_DIF_3TWIDDLE(r0,c0): */\
				/* Do	the p0,p2 combo: */\
				"ldp	q4,q5,[x2]			\n\t"\
				"ldp	q8,q9,[x4]			\n\t"/* cc0 */\
				"ldp	q0,q1,[x0]			\n\t"\
				"fmul	v6.2d,v4.2d,v8.2d	\n\t"/* twiddle-mul: */\
				"fmul	v7.2d,v5.2d,v8.2d	\n\t"\
				"fmls	v6.2d,v5.2d,v9.2d	\n\t"\
				"fmla	v7.2d,v4.2d,v9.2d	\n\t"\
				"fsub	v2.2d ,v0.2d,v6.2d	\n\t"/* 2 x 2 complex butterfly: */\
				"fsub	v3.2d ,v1.2d,v7.2d	\n\t"\
				"fadd	v10.2d,v0.2d,v6.2d	\n\t"\
				"fadd	v11.2d,v1.2d,v7.2d	\n\t"\
				/* Do	the p1,3 combo: */\
				"ldp	q8,q9,[x4,#0x40]	\n\t"/* cc0+4 */\
				"ldp	q6,q7,[x3]			\n\t"\
				"fmul	v0.2d,v6.2d,v8.2d	\n\t"/* twiddle-mul: */\
				"fmul	v1.2d,v7.2d,v8.2d	\n\t"\
				"fmls	v0.2d,v7.2d,v9.2d	\n\t"\
				"fmla	v1.2d,v6.2d,v9.2d	\n\t"\
				"ldp	q8,q9,[x4,#0x20]	\n\t"/* cc0+2 */\
				"ldp	q6,q7,[x1]			\n\t"\
				"fmul	v4.2d,v6.2d,v8.2d	\n\t"/* twiddle-mul: */\
				"fmul	v5.2d,v7.2d,v8.2d	\n\t"\
				"fmls	v4.2d,v7.2d,v9.2d	\n\t"\
				"fmla	v5.2d,v6.2d,v9.2d	\n\t"\
				"fadd	v6.2d,v4.2d,v0.2d	\n\t"/* 2 x 2 complex butterfly: */\
				"fadd	v7.2d,v5.2d,v1.2d	\n\t"\
				"fsub	v4.2d,v4.2d,v0.2d	\n\t"\
				"fsub	v5.2d,v5.2d,v1.2d	\n\t"\
				/* Finish radix-4 butterfly and store results: */\
				"fsub	v8.2d,v10.2d,v6.2d	\n\t"\
				"fsub	v9.2d,v11.2d,v7.2d	\n\t"\
				"fsub	v1.2d,v3.2d,v4.2d	\n\t"\
				"fsub	v0.2d,v2.2d,v5.2d	\n\t"\
				"fadd	v6.2d,v6.2d,v10.2d	\n\t"\
				"fadd	v7.2d,v7.2d,v11.2d	\n\t"\
				"fadd	v4.2d,v4.2d,v3.2d	\n\t"\
				"fadd	v5.2d,v5.2d,v2.2d	\n\t"\
				"stp	q6,q7,[x5      ]	\n\t"/* out 0 */\
				"stp	q0,q4,[x5,#0x20]	\n\t"/* out 1 */\
				"stp	q8,q9,[x5,#0x40]	\n\t"/* out 2 */\
				"stp	q5,q1,[x5,#0x60]	\n\t"/* out 3 */\
				:					/* outputs: none */\
				: [__add0] "m" (add0)	/* All inputs from memory addresses here */\
				 ,[__p1] "m" (p1)\
				 ,[__p2] "m" (p2)\
				 ,[__p3] "m" (p3)\
				 ,[__two] "m" (two)\
				 ,[__cc0] "m" (cc0)\
				 ,[__r0] "m" (r0)\
				: "cc","memory","x0","x1","x2","x3","x4","x5","v0","v1","v2","v3","v4","v5","v6","v7","v8","v9","v10","v11"	/* Clobbered registers */\
			);
		#else
			__asm__ volatile (\
				"movq	%[__cc0],%%rsi 		\n\t"\
				"movq	%[__add0],%%rax		\n\t"\
				/* NB: preshifting the p-offsets to save the *= 8 below saves nothing: */\
				"movslq	%[__p1],%%rbx		\n\t"\
				"movslq	%[__p2],%%rcx		\n\t"\
				"movslq	%[__p3],%%rdx		\n\t"\
				"leaq	(%%rax,%%rbx,8),%%rbx	\n\t"\
				"leaq	(%%rax,%%rcx,8),%%rcx	\n\t"\
				"leaq	(%%rax,%%rdx,8),%%rdx	\n\t"\
				/* SSE2_RADIX_04_DIF_3TWIDDLE(r0,c0): */\
				/* Do	the p0,p2 combo: */\
				"movaps	    (%%rcx),%%xmm4	\n\t"\
				"movaps	0x10(%%rcx),%%xmm5	\n\t"\
				"movaps	    (%%rsi),%%xmm2	\n\t"\
				"movaps	0x10(%%rsi),%%xmm3	\n\t"\
				"movaps	%%xmm4,%%xmm6		\n\t"\
				"movaps	%%xmm5,%%xmm7		\n\t"\
				"movaps	    (%%rax),%%xmm0	\n\t"\
				"movaps	0x10(%%rax),%%xmm1	\n\t"\
				"mulpd	%%xmm2,%%xmm4		\n\t"\
				"mulpd	%%xmm2,%%xmm5		\n\t"\
				"mulpd	%%xmm3,%%xmm6		\n\t"\
				"mulpd	%%xmm3,%%xmm7		\n\t"\
				"movaps	%%xmm0,%%xmm2		\n\t"\
				"movaps	%%xmm1,%%xmm3		\n\t"\
				"movq	%[__r0],%%rdi 		\n\t"\
				"addpd	%%xmm6,%%xmm5		\n\t"\
				"subpd	%%xmm7,%%xmm4		\n\t"\
				"addpd	%%xmm4,%%xmm0		\n\t"\
				"addpd	%%xmm5,%%xmm1		\n\t"\
				"subpd	%%xmm4,%%xmm2		\n\t"\
				"subpd	%%xmm5,%%xmm3		\n\t"\
				"movaps	%%xmm0,0x40(%%rdi)	\n\t"/* Spill 1: free up xmm0,1 */\
				"movaps	%%xmm1,0x50(%%rdi)	\n\t"\
				/* Do	the p1,3 combo: */\
				"movaps	0x40(%%rsi),%%xmm0	\n\t"\
				"movaps	0x50(%%rsi),%%xmm1	\n\t"\
				"movaps	    (%%rdx),%%xmm6	\n\t"\
				"movaps	0x10(%%rdx),%%xmm7	\n\t"\
				"movaps	%%xmm6,%%xmm4		\n\t"\
				"movaps	%%xmm7,%%xmm5		\n\t"\
				"mulpd	%%xmm0,%%xmm4		\n\t"\
				"mulpd	%%xmm0,%%xmm5		\n\t"\
				"mulpd	%%xmm1,%%xmm6		\n\t"\
				"mulpd	%%xmm1,%%xmm7		\n\t"\
				"addpd	%%xmm6,%%xmm5		\n\t"\
				"subpd	%%xmm7,%%xmm4		\n\t"\
				"movaps	%%xmm5,0x10(%%rdi)	\n\t"/* Spill 2*/\
				"movaps	%%xmm4,    (%%rdi)	\n\t"\
				"movaps	0x20(%%rsi),%%xmm0	\n\t"\
				"movaps	0x30(%%rsi),%%xmm1	\n\t"\
				"movaps	    (%%rbx),%%xmm6	\n\t"\
				"movaps	0x10(%%rbx),%%xmm7	\n\t"\
				"movaps	%%xmm6,%%xmm4		\n\t"\
				"movaps	%%xmm7,%%xmm5		\n\t"\
				"mulpd	%%xmm0,%%xmm4		\n\t"\
				"mulpd	%%xmm0,%%xmm5		\n\t"\
				"mulpd	%%xmm1,%%xmm6		\n\t"\
				"mulpd	%%xmm1,%%xmm7		\n\t"\
				"movaps	    (%%rdi),%%xmm0	\n\t"/* Restore 2 */\
				"movaps	0x10(%%rdi),%%xmm1	\n\t"\
				"addpd	%%xmm6,%%xmm5		\n\t"\
				"subpd	%%xmm7,%%xmm4		\n\t"\
				"movaps	%%xmm5,%%xmm7		\n\t"\
				"movaps	%%xmm4,%%xmm6		\n\t"\
				"subpd	%%xmm0,%%xmm4		\n\t"\
				"subpd	%%xmm1,%%xmm5		\n\t"\
				"addpd	%%xmm0,%%xmm6		\n\t"\
				"addpd	%%xmm1,%%xmm7		\n\t"\
				/* Finish radix-4 butterfly and store results: */\
				"movq	%[__two],%%rsi		\n\t"\
				"movaps	0x40(%%rdi),%%xmm0	\n\t"/* Restore 1 */\
				"movaps	0x50(%%rdi),%%xmm1	\n\t"\
				"subpd	%%xmm6,%%xmm0		\n\t"\
				"subpd	%%xmm5,%%xmm2		\n\t"\
				"subpd	%%xmm7,%%xmm1		\n\t"\
				"subpd	%%xmm4,%%xmm3		\n\t"\
				"movaps	%%xmm0,0x40(%%rdi)	\n\t	movaps	(%%rsi),%%xmm0	\n\t"/* 2.0 */\
				"movaps	%%xmm2,0x20(%%rdi)	\n\t"\
				"movaps	%%xmm1,0x50(%%rdi)	\n\t"\
				"movaps	%%xmm3,0x70(%%rdi)	\n\t"\
				"mulpd	%%xmm0,%%xmm6		\n\t"\
				"mulpd	%%xmm0,%%xmm5		\n\t"\
				"mulpd	%%xmm0,%%xmm7		\n\t"\
				"mulpd	%%xmm0,%%xmm4		\n\t"\
				"addpd	0x40(%%rdi),%%xmm6	\n\t"\
				"addpd		%%xmm2 ,%%xmm5	\n\t"\
				"addpd		%%xmm1 ,%%xmm7	\n\t"\
				"addpd		%%xmm3 ,%%xmm4	\n\t"\
				"movaps	%%xmm6,    (%%rdi)	\n\t"\
				"movaps	%%xmm5,0x60(%%rdi)	\n\t"\
				"movaps	%%xmm7,0x10(%%rdi)	\n\t"\
				"movaps	%%xmm4,0x30(%%rdi)	\n\t"\
				:					/* outputs: none */\
				: [__add0] "m" (add0)	/* All inputs from memory addresses here */\
				 ,[__p1] "m" (p1)\
				 ,[__p2] "m" (p2)\
				 ,[__p3] "m" (p3)\
				 ,[__two] "m" (two)\
				 ,[__cc0] "m" (cc0)\
				 ,[__r0] "m" (r0)\
				: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
			);
		#endif	// ARM_V8 or X86_64 SIMD?
		}
		clock2 = getRealTime();
		tdiff = (double)(clock2 - clock1);
		printf("%s, loop 1: Time for %u macro calls =%s\n",func, imax, get_time_str(tdiff));
		// Check outputs vs ref-data:
		nerr = 0;
		dptr = (double *)r0;
		for(i = 0; i < 4; i++) {
			j = i<<2;
		//	printf("Out%u Re.[d0,d1] = %16.12f,%16.12f, Im.[d0,d1] = %16.12f,%16.12f\n",i,*dptr,*(dptr+1),*(dptr+2),*(dptr+3));
			nerr += (fabs(*(dptr  ) - ref1[j  ]) > 1e-10);
			nerr += (fabs(*(dptr+1) - ref1[j+1]) > 1e-10);
			nerr += (fabs(*(dptr+2) - ref1[j+2]) > 1e-10);
			nerr += (fabs(*(dptr+3) - ref1[j+3]) > 1e-10);
			dptr += 4;
		}
		ASSERT(HERE, nerr == 0, "Outputs mismatch ref-data!");

	// Timing loop #2 - two radix-4 DFTs (operating on separate data chunks but sharing twiddles) side-by-side:
		/* 6 May 2016, Core2:
		Baseline single-DFT timing of loop #2 = 0.416 sec with-DFT, 0.252 sec sans-DFT ==> 1-DFT = 0.164 sec ==> 32.8 cycles

		*/
		// 7 May 2016: 10^7-loop timings, 1-threaded on my 2GHz Core2Duo:
		//																comments
		// no-DFT timing (keeping just 1 address-load of macro)	0.486	97.2 cycles ... 2x as many data-inits as loop #1
		// initial timing										0.871	126.2 - 49.2 = 77 cycles, vs 36.8 for single-column ASM
		//											*** 2-column actually slows things down in terms of per-cycle throughput! ***
		// rcol down 4 lines to break same-instruction blocks	0.813	66 cycles, big improvement!
		//18 May: 2-col impl of 2x4-DFT-each-with-3-twiddles	0.791	61 cycles, better, need to get under 60

		// 9 May: GCCified George's dual complex 4-DFT macro (not a drop-in replacement candidate for mine
		// (due to its cotangent-twiddles scheme), munged addresses, pasted underneath my opening address-
		// computation block below, just to get a timing:		0.735	50 cycles, yowza!
		//
		// Copy quasirandom digits-of-Pi-data into our vec_dbl inputs:
		for(j = 0; j < dim+dim; j++) { *(add0+j) = ran[j]; }
		int k = 0x60;
		clock1 = getRealTime();
		for(i = 0; i < imax; i++) {
			j = 0x80;	// Set j to the bytewise address offset between the pointers to the first and second DFT's data
						// 0x60 = literal-bytewise address offset between the pointers to the first and second DFT's twiddles:
			SSE2_RADIX_04_DIF_3TWIDDLE_X2(add0,add1,add2,add3,j, two,cc0,k, r0,r1,r2,r3,j)
		}
		clock2 = getRealTime();
		tdiff = (double)(clock2 - clock1);
		printf("%s, loop 2: Time for %u macro calls =%s\n",func, imax, get_time_str(tdiff));
		// Check outputs vs ref-data:
		nerr = 0;
		dptr = (double *)r0;
		for(i = 0; i < 8; i++) {
			j = i<<2;
		//	printf("Out%u Re.[d0,d1] = %16.12f,%16.12f, Im.[d0,d1] = %16.12f,%16.12f\n",i,*dptr,*(dptr+1),*(dptr+2),*(dptr+3));
			nerr += (*dptr != ref1[j]) + (*(dptr+1) != ref1[j+1]) + (*(dptr+2) != ref1[j+2]) + (*(dptr+3) != ref1[j+3]);
			dptr += 4;
		}
		ASSERT(HERE, nerr == 0, "Outputs mismatch ref-data!");

	// Timing loop #3 - single radix-4 DIT DFT:
		dim = 8*RE_IM_STRIDE;	// 4 vector-complex data
		// Copy quasirandom digits-of-Pi-data into our vec_dbl inputs:
		for(j = 0; j < dim; j++) { *(add0+j) = ran[j]; }
		clock1 = getRealTime();
		for(i = 0; i < imax; i++) {
			// Single-4-DFT macro uses same arglist as paired ones, but ignores the args in the j,k,j slots:
			SSE2_RADIX_04_DIT_3TWIDDLE_X1(add0,add1,add2,add3,j, two,cc0,k, r0,r1,r2,r3,j)
		}
		clock2 = getRealTime();
		tdiff = (double)(clock2 - clock1);
		printf("%s, loop 3: Time for %u macro calls =%s\n",func, imax, get_time_str(tdiff));
		// Check outputs vs ref-data:
		nerr = 0;
		dptr = (double *)r0;
		for(i = 0; i < 4; i++) {
			j = i<<2;
		//	printf("Out%u Re.[d0,d1] = %16.12f,%16.12f, Im.[d0,d1] = %16.12f,%16.12f\n",i,*dptr,*(dptr+1),*(dptr+2),*(dptr+3));
			nerr += (*dptr != ref2[j]) + (*(dptr+1) != ref2[j+1]) + (*(dptr+2) != ref2[j+2]) + (*(dptr+3) != ref2[j+3]);
			dptr += 4;
		}
		ASSERT(HERE, nerr == 0, "Outputs mismatch ref-data!");

	// Timing loop #4 - two radix-4 DIT DFTs (operating on separate data chunks but sharing twiddles) side-by-side:
		for(j = 0; j < dim+dim; j++) { *(add0+j) = ran[j]; }
		k = 0x60;
		clock1 = getRealTime();
		for(i = 0; i < imax; i++) {
			j = 0x80;	// Set j to the bytewise address offset between the pointers to the first and second DFT's data
						// 0x60 = literal-bytewise address offset between the pointers to the first and second DFT's twiddles:
			SSE2_RADIX_04_DIT_3TWIDDLE_X2(add0,add1,add2,add3,j, two,cc0,k, r0,r1,r2,r3,j)
		}
		clock2 = getRealTime();
		tdiff = (double)(clock2 - clock1);
		printf("%s, loop 4: Time for %u macro calls =%s\n",func, imax, get_time_str(tdiff));
		// Check outputs vs ref-data:
		nerr = 0;
		dptr = (double *)r0;
		for(i = 0; i < 8; i++) {
			j = i<<2;
		//	printf("Out%u Re.[d0,d1] = %16.12f,%16.12f, Im.[d0,d1] = %16.12f,%16.12f\n",i,*dptr,*(dptr+1),*(dptr+2),*(dptr+3));
			nerr += (*dptr != ref2[j]) + (*(dptr+1) != ref2[j+1]) + (*(dptr+2) != ref2[j+2]) + (*(dptr+3) != ref2[j+3]);
			dptr += 4;
		}
		ASSERT(HERE, nerr == 0, "Outputs mismatch ref-data!");

		free((void *)sc_arr);	sc_arr=0x0;
		return nerr;
	#endif	// USE_AVX?
	}

	// Timing loop for radix-16 DIF macro:
	int	test_radix16_dft()
	{
		const char func[] = "test_radix16_dft";
		/*...time-related stuff	*/
		double clock1, clock2;
		double *dptr, tdiff, rt,it, dtmp, avg_err;
		double t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31;
		int i,j,j1,j2,k,imax = 10000000, nerr = 0;	// Use 10^7 loop execs
		int p1,p2,p3,p4,p5,p6,p7,p8,p9,pA,pB,pC,pD,pE,pF;
		const double c = 0.9238795325112867561281831, s = 0.3826834323650897717284599;	/* exp[i*(twopi/16)]*/
	#ifdef USE_AVX2	// FMA-based DFT needs the tangent
	  #ifdef REFACTOR_4DFT_3TWIDDLE
		#error USE_AVX2 and REFACTOR_4DFT_3TWIDDLE are mutually exclusive preprocessor directives!
	  #endif
		const double tan = 0.41421356237309504879;
	#endif
		// DIF[ref1], DIT[ref2] ref-outputs: cols 0,1 are re,im outputs for scalar-mode,
		// cols [0,1],[2,3] are [re0,im0],[re1,im1] for SSE2 mode,
		// cols [0,1],[2,3],[4,5],[6,7] are [re0,im0],[re1,im1],[re2,im2],[re3,im3] for AVX/AVX2 mode:
		const double ref1[] = {	// DIF ref-outputs:
			 71.836132626544, 51.240113957019,  93.778630517480, 76.253829098865,  84.806628663973, 73.237374006630,  75.813281246317, 71.186791613560,
			 -9.999128163205, -5.041057027004, -11.993749676902, -0.051731207851,  -2.997214610252, -5.029166077852, -11.965466770155, -1.019849970447,
			-14.872954615349,-26.053387304330,  -2.971220181207,  0.996282812198,  -6.000050539329, 22.005043558985,  -5.934761044083, -9.023636204473,
			  5.023325346528, -7.980009625544,  20.909799869969, 19.046570941658,  -0.036810883888,  5.988776862464,  10.025385052093, -5.002254733787,
			  9.056346572622, -1.914547592304,   3.981734242952,  8.184718642903,   5.433348885097,  0.012717649256,  11.329786667132, 10.968444860577,
			-15.077818730218,  7.893027085061,  -5.978667607114,  3.812269869915,  -7.427204810541, 23.968911024385, -11.266969635162, -8.942204365409,
			  0.293315070044,  4.170949230343, -12.723161421476,-21.587788299675,  -9.388763668691, -2.159379128165, -21.281854472092, -0.810792653682,
			  1.740668944464,-14.241457760931,   2.800021476805,-14.482647475474, -16.537232610515,  2.000039861260,   1.280448357857, -9.258323929924,
			 -3.152128563857,-16.771900048992,  -5.521042109482,  7.816522698106,   2.280114913431, 14.541287717387, -10.102184012659, 17.040342071843,
			  2.073046947819,  1.969010269254,  -2.360249253366,-17.814992092284,   8.483022773111, -1.532157855444,  -0.178724865422, -2.336467764655,
			 -9.872285434976,  1.404582139871, -13.035858498590, -8.152619897216,   3.973177562804, -9.433750344665,   8.526197273669, 21.162355221723,
			  3.080286869058,-10.650483193072,   5.076556271954, -1.854679849975,  29.478126055897,-11.446050304373, -30.263630386603, 16.060181777754,
			 14.065076567172,  2.036958078982,  -7.839406695376,-20.597474713295,  10.403805832586, 11.790086387134,  -0.622179305579,  1.390534729868,
			 21.045817946931,  5.610880426309,  -7.416729930697,  2.536423857502,   0.585026163834,  4.233047491042, -11.128509999099, -0.591596670215,
			-12.866425691057,  6.044821934737,  12.200367325354, -1.892360702488, -31.099488921485,  4.039179558341,  11.722605670925, 11.014905859604,
			-14.373275692520, 18.282499430601,  -4.907024330304,-16.212323682887,   8.043515193967, 11.784039593614,  16.046576222863,-15.838429842336
		};
		const double ref2[] = {	// DIT ref-outputs
			 72.000000000000, 51.000000000000,  94.000000000000, 76.000000000000,  85.000000000000, 73.000000000000,  76.000000000000, 71.000000000000,
			 -3.186854935858,-14.814979432853,  16.274310380149,  1.335323711456, -13.073174362886, -2.237484334246,  -1.658016565491, -7.214637436104,
			 -6.722101857669,  7.590944457452,  14.087819349234, -5.303699999882,   2.642267168043,  8.119296138409,   8.194965871283,  3.357676558811,
			 -2.410055786718, 17.860460799853,  -3.115089210758,  8.081417051306,  10.217674711889, 18.135933290759, -23.560603557289, 12.034604903514,
			 -4.013801115483, -8.993853490320,  -3.019938212777,-12.995382764313, -17.975435139150, 16.027592818580,   4.984654315396,-10.007658135448,
			  9.628673950784, -2.618390257617,  -4.977626115382, -7.585510986880,   8.576400244341, 14.894373226740, -32.710474648103,  9.946362646980,
			  1.407075652116, -3.103750917409,  -2.611395268707, 10.713143881330,  19.283439756536,  8.662759113515, -19.031759359980, -4.734395634916,
			  1.513899588864,  8.129850867985,  -9.169123259445,-16.471534197115,   0.296631561689, -6.875402885079,   2.295297934556, -0.702929830780,
			-20.015245660006, -4.938617333789,   7.987690523425, -4.024524829342,   9.033705168679, 10.972336621039,  18.045934640017, 14.944706185410,
			 -9.250056079529, 15.676647091397,   3.540556362652,  4.445231564276, -33.237348371394, -7.299683015832,  16.814644365039, 13.500464006934,
			 18.767722258015, 10.342316476739, -24.117427303634, -6.614666472964,  15.378722566760,  3.819731226687, -10.228223459299, -9.324804950323,
			-12.959770458085,-18.146313510762,  -3.761348735601,-13.576800261814, -11.015270817590,  2.384093214004,   0.269727641998,  2.133525968650,
			  4.013763422689,  2.981560528764,   1.022999041684,  4.995345129289,  11.999872932983, -0.055223113445,   4.981539350928, -4.022967274930,
			 -5.159441616098, -2.202809827852, -14.966371172025,-26.129659478906,   1.670348376017,-13.262428826303, -14.344280414512,  7.844841106696,
			 -1.537138437596,-22.891572112364, -15.313846915632,  9.375246960644, -13.224109940089,  7.363998369751,   5.078609933375, 10.751062960496,
			  5.707698152898,-19.823900762950,   7.918716323359, -6.034843794865,   4.573028564873, 10.363166029094,  -3.060761887966,-13.427985421598
		};

		const int stride = 2*RE_IM_STRIDE, dim = stride<<4;
		double c1,c2,c3,c4,c5,c6,c7,c8,c9,cA,cB,cC,cD,cE,cF, s1,s2,s3,s4,s5,s6,s7,s8,s9,sA,sB,sC,sD,sE,sF;
		static double *a,*a_ptr;	// Dimension = number of scalar-doubles in 16 vector-complex in SIMD build mode
		a_ptr = ALLOC_VEC_DBL(a_ptr, dim/RE_IM_STRIDE);	if(!a_ptr){ sprintf(cbuf, "FATAL: unable to allocate a_ptr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		a     = ALIGN_VEC_DBL(a_ptr);
		ASSERT(HERE, ((long)a & SZ_VDM1) == 0, "a0_ptr not 64-byte aligned!");
	#ifdef USE_SSE2
		const int pfetch_dist = 0;
		int pfetch_addr = 0;	// Don't care about pfetch in this lcal-mem context, so just set these = 0
		static vec_dbl *sc_arr = 0x0, *sc_ptr;
		double *add0,*add1,*add2;	/* Addresses into array sections */
		vec_dbl *c_tmp,*s_tmp, *i0,*i1,*i2,*i3, *o0,*o1,*o2,*o3;
		static vec_dbl *cc0, *ss0, *isrt2, *two, *r00;
		sc_arr = ALLOC_VEC_DBL(sc_arr, 72);	if(!sc_arr){ sprintf(cbuf, "FATAL: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		sc_ptr = ALIGN_VEC_DBL(sc_arr);
		ASSERT(HERE, ((long)sc_ptr & SZ_VDM1) == 0, "sc_ptr not 64-byte aligned!");
		r00 = sc_ptr + 0x00;	  isrt2 = sc_ptr + 0x20;
									cc0 = sc_ptr + 0x21;
									ss0 = sc_ptr + 0x22;
									two = sc_ptr + 0x43;
		/* These remain fixed: */
		VEC_DBL_INIT(isrt2, ISRT2);
	  #if defined(USE_AVX2) && !defined(REFACTOR_4DFT_3TWIDDLE)
		VEC_DBL_INIT(two  , 1.0);
		// cc0,ss0 inited below for AVX2
	  #else
		VEC_DBL_INIT(two  , 2.0);
		VEC_DBL_INIT(cc0  , c);
		VEC_DBL_INIT(ss0  , s);
	  #endif
	#endif	// USE_SSE2 ?

		// Do these timing-test DFTs in-place, i.e. p1 = #doubles in a pair of vec_dbl:
		p1 = dim>>4;	// Set stride equal to the AVX complex-vec_dbl value, 2*4
		p2 = p1 +p1;
		p3 = p2 +p1;
		p4 = p3 +p1;
		p5 = p4 +p1;
		p6 = p5 +p1;
		p7 = p6 +p1;
		p8 = p7 +p1;
		p9 = p8 +p1;
		pA = p9 +p1;
		pB = pA +p1;
		pC = pB +p1;
		pD = pC +p1;
		pE = pD +p1;
		pF = pE +p1;

		// Twiddles for the purpose of this timing-test are w1-15, with w := exp(2*Pi*I/2^14):
		c1 = 0.9999999264657178511447314807; s1 = 0.0003834951875713955890724616812;
		c2 = 0.9999997058628822191602282177; s2 = 0.0007669903187427045269385683580;
		c3 = 0.9999993381915255477888066109; s3 = 0.001150485337113848457071735047;
		c4 = 0.9999988234517019099290257101; s4 = 0.001533980186284765612303697150;
		c5 = 0.9999981616434870076277347923; s5 = 0.001917474809855419109500620455;
		c6 = 0.9999973527669781720689399696; s6 = 0.002300969151425805244235552264;
		c7 = 0.9999963968222943635594898320; s7 = 0.002684463154595961785455992532;
		c8 = 0.9999952938095761715115801256; s8 = 0.003067956762965976270145365491;
		c9 = 0.9999940437289858144220774704; s9 = 0.003451449920135994297977171937;
		cA = 0.9999926465807071398486621178; sA = 0.003834942569706227825960602960;
		cB = 0.9999911023649456243827897550; sB = 0.004218434655276963463076393843;
		cC = 0.9999894110819283736194723572; sC = 0.004601926120448570764901699143;
		cD = 0.9999875727319041221238780943; sD = 0.004985416908821510528222769585;
		cE = 0.9999855873151432333947502950; sE = 0.005368906963996343085634209014;
		cF = 0.9999834548319376998246454755; sF = 0.005752396229573736600123594041;
		// In the refactor-case we restructure these and overwrite w5-7,9-B,D-F as shown:
	#ifdef REFACTOR_4DFT_3TWIDDLE
		// w1.E^2 = [c1,s1].[1,I]/sqrt2 = [c1-s1,c1+s1]*ISRT2;
		// w2.E^4 = [c2,s2].I = [-s2,c2];
		// w3.E^6 = [c3,s3].[-1,I]/sqrt2 = [-c3-s3,c3-s3]*ISRT2
		c5 =  (c1-s1)*ISRT2;	s5 = (c1+s1)*ISRT2;
		c6 = -s2;				s6 = c2;
		c7 = -(c3+s3)*ISRT2;	s7 = (c3-s3)*ISRT2;

		// w1.E   = [c1,s1].[c,s] = [c1.c-s1.s,c1.s+s1.c];
		// w2.E^2 = [c2,s2].[2,I]/sqrt2 = [c2-s2,c2+s2]*ISRT2;
		// w3.E^3 = [c3,s3].[s,c] = [c3.s-s3.c,c3.c+s3.s]
		c9 = c1*c-s1*s;			s9 = c1*s+s1*c;
		cA = (c2-s2)*ISRT2;		sA = (c2+s2)*ISRT2;
		cB = c3*s-s3*c;			sB = c3*c+s3*s;

		// w1.E^3 = [c1,s1].[s,c] = [c1.s-s1.c,c1.c+s1.s];
		// w2.E^6 = [c2,s2].[-1,I]/sqrt2 = [-c2-s2,c2-s2]*ISRT2;
		// w3.-E  = [c3,s3].[-c,-s] = [-c3.c+s3.s,-c3.s-s3.c]
		cD = c1*s-s1*c;			sD = c1*c+s1*s;
		cE = -(c2+s2)*ISRT2;	sE = (c2-s2)*ISRT2;
		cF = -c3*c+s3*s;		sF = -(c3*s+s3*c);
	/*
		// And now do batch-div to setup for 'cotangent'-CMUL scheme:
		c1 /= s1;
		c2 /= s2;
		c3 /= s3;
		c4 /= s4;
		c5 /= s5;
		c6 /= s6;
		c7 /= s7;
		c8 /= s8;
		c9 /= s9;
		cA /= sA;
		cB /= sB;
		cC /= sC;
		cD /= sD;
		cE /= sE;
		cF /= sF;
	*/
	#endif

	#ifdef USE_SSE2

	  #ifdef USE_AVX2	// AVX2/FMA needs tangent-form twiddles:
		// In AVX2/FMA mode, since we need to replace most of the raw sincos data with derived ones,
		// simply place one copy of each scalar-double in a double-sized slot of the local memory.
		// We will be using AVX2/FMA-based Newtonian iterative inversion on the 16 doubles whose
		// multiplicative inverse is needed (the real part of the basic root of unity c and of the
		// 15 complex twiddles, c1-15), so store those in packed form in 4 AVX-register-sized
		// contiguous memory locations, and the others in a separate chunk of memory. After the
		// vector-iterative inversion we'll need to combine the 2 sets of data and place (in suitable
		// vector-register-sized broadcast form) into their final SIMD-suitable memory slots.

	  clock1 = getRealTime();
	  for(i = 0; i < imax; i++) {	// repeater loop
		add0 = (double *)cc0;	// add0 points to 16 cos-data-to-be-inverted; Need a double-ptr on lhs here
		add1 = add0 + 16;	// add1 points to block of memory temporarily used to store the corresponding sine data
		add2 = add0 + 32;	// add2 points to block of memory temporarily used to store the 11 [0-padded to 12]
							//	cosine data which need to be divided by other cosines (i.e. multiplied by inverses)
		/* The add2-addressed cosine ratios are arranged in 3 YMM-register/memory-sized slots like so;
		  once we have filled 4 YYMs with inverses 1/[c3,c1-15] and used those to get the 16 tangents (1st set = 1/c3
		  and discarded) we will do as described in the right column to set up for the cosine-ratios computation:

			double __c31 = __c3/__c1;
			double __c51 = __c5/__c1;
			double __c62 = __c6/__c2;
			[0 pad]						shuffle YMM with 1/[c3,c1,c2,c3] to get 1/[c1,c1,c2,c3], then *= [c3,c5,c6,0]

			double __c73 = __c7/__c3;
			double __c91 = __c9/__c1;
			double __cA2 = __cA/__c2;
			double __cB3 = __cB/__c3;	initialize YMM with 1/[c3,c1,c2,c3], then *= [c7,c9,cA,cB]

			double __cC4 = __cC/__c4;
			double __cD5 = __cD/__c5;
			double __cE6 = __cE/__c6;
			double __cF7 = __cF/__c7;	Multiply YMM with 1/[c4-7] *= [cC-F]
		*/
		// Since tan0 defined as const, use this pair of double slots to hold 1/c3 (via c3,1 on input, then invert c3 and multiply
		// them together), which extra 1/c3 copy saves some really awkward permuting, at least in terms of the idiotic x86 ISA.
		*add0++ = 0.0;	*add1++ = 1.0;
		*add0++ = c1;	// c1, for inversion
		*add1++ = s1;	// s1  slot will hold __r1 = s1 /c1
		*add0++ = c2;	// c2, for inversion
		*add1++ = s2;	// s2  slot will hold __r2 = s2 /c2
		*(add0-3) = c3;	// c3, for inversion ...
		*add0++   = c3;	// place extra copy in 0-slot as described above - put on separate line to avoid ambiguity of *(add0-3) = *add0++ = ...
		*add1++ = s3;	// s3  slot will hold __r3 = s3 /c3
		*add2++ = c3;	// c3, will get multiplied by 1/c1 to yield __c31
		*add0++ = c4;	// c4, for inversion
		*add1++ = s4;	// s4  slot will hold __r4 = s4 /c4
		*add0++ = c5;	// c5, for inversion
		*add1++ = s5;	// s5  slot will hold __r5 = s5 /c5
		*add2++ = c5;	// c5, will get multiplied by 1/c1 to yield __c51
		*add0++ = c6;	// c6, for inversion
		*add1++ = s6;	// s6  slot will hold __r6 = s6 /c6
		*add2++ = c6;	// c6, will get multiplied by 1/c2 to yield __c62
		*add2++ = 0.0;	// 0-pad will get multiplied by 1/c3 term, remains 0-pad.
		*add0++ = c7;	// c7, for inversion
		*add1++ = s7;	// s7  slot will hold __r7 = s7 /c7
		*add2++ = c7;	// c7, will get multiplied by 1/c3 to yield __c73
		*add0++ = c8;	// c8, for inversion
		*add1++ = s8;	// s8  slot will hold __r8 = s8 /c8
		*add0++ = c9;	// c9, for inversion
		*add1++ = s9;	// s9  slot will hold __r9 = s9 /c9
		*add2++ = c9;	// c9, will get multiplied by 1/c1 to yield __c91
		*add0++ = cA;	// c10, for inversion
		*add1++ = sA;	// s10 slot will hold __rA = s10/c10
		*add2++ = cA;	// c10, will get multiplied by 1/c2 to yield __cA2
		*add0++ = cB;	// c11, for inversion
		*add1++ = sB;	// s11 slot will hold __rB = s11/c11
		*add2++ = cB;	// c11, will get multiplied by 1/c3 to yield __cB3
		*add0++ = cC;	// c12, for inversion
		*add1++ = sC;	// s12 slot will hold __rC = s12/c12
		*add2++ = cC;	// c12, will get multiplied by 1/c4 to yield __cC4
		*add0++ = cD;	// c13, for inversion
		*add1++ = sD;	// s13 slot will hold __rD = s13/c13
		*add2++ = cD;	// c13, will get multiplied by 1/c5 to yield __cD5
		*add0++ = cE;	// c14, for inversion
		*add1++ = sE;	// s14 slot will hold __rE = s14/c14
		*add2++ = cE;	// c14, will get multiplied by 1/c6 to yield __cE6
		*add0++ = cF;	// c15, for inversion
		*add1++ = sF;	// s15 slot will hold __rF = s15/c15
		*add2++ = cF;	// c15, will get multiplied by 1/c7 to yield __cF7
		/*
		At this point, the 11 ymm-sized [32-byte] chunks starting at &cc0 contain the following scalar-double data:

		0:	c3,c1-3
		1:	c4-7
		2:	c8-11
		3:	c12-c15
		4:	1.0,s1-3
		5:	s4-7
		6:	s8-11
		7:	s12-s15
		8:	c3,5,6,[0-pad]
		9:	c7,9-B
		A:	cC-F
		*/
		c_tmp = &c; s_tmp = &tan;	// GCC/Clang don't allow cd_address-taking inlined in arglist of macros, so do it here
		RADIX16_COMPUTE_FMA_SINCOS_DIF(cc0,two, c_tmp,s_tmp);
	  }	// repeater loop
		clock2 = getRealTime();
		tdiff = (double)(clock2 - clock1);
		printf("%s: Time for %u FMA DIF tan-twiddles setup calls =%s [tdiff = %20.10e]\n",func, imax, get_time_str(tdiff), tdiff);

	  #elif defined(REFACTOR_4DFT_3TWIDDLE)
		/* Sincos data stored in terms of the following 5 contiguous-data triplets:
			c4,s4, c8,s8, cC,sC
			c1,s1, c2,s2, c3,s3
			c5,s5, c6,s6, c7,s7
			c9,s9, cA,sA, cB,sB
			cD,sD, cE,sE, cF,sF .
		Note that due to my layout of the SSE2_RADIX_04_DIF_3TWIDDLE_X2-macro arglist,
		we need to swap the order of the first 2 sincos-pairs of each triplet:
		*/
		c_tmp = cc0; s_tmp = c_tmp+1;	/* c0,s0 */
		VEC_DBL_INIT(c_tmp, c8);	VEC_DBL_INIT(s_tmp, s8);	c_tmp+=2; s_tmp+=2;
		VEC_DBL_INIT(c_tmp, c4);	VEC_DBL_INIT(s_tmp, s4);	c_tmp+=2; s_tmp+=2;
		VEC_DBL_INIT(c_tmp, cC);	VEC_DBL_INIT(s_tmp, sC);	c_tmp+=2; s_tmp+=2;

		VEC_DBL_INIT(c_tmp, c2);	VEC_DBL_INIT(s_tmp, s2);	c_tmp+=2; s_tmp+=2;
		VEC_DBL_INIT(c_tmp, c1);	VEC_DBL_INIT(s_tmp, s1);	c_tmp+=2; s_tmp+=2;
		VEC_DBL_INIT(c_tmp, c3);	VEC_DBL_INIT(s_tmp, s3);	c_tmp+=2; s_tmp+=2;

		VEC_DBL_INIT(c_tmp, c6);	VEC_DBL_INIT(s_tmp, s6);	c_tmp+=2; s_tmp+=2;
		VEC_DBL_INIT(c_tmp, c5);	VEC_DBL_INIT(s_tmp, s5);	c_tmp+=2; s_tmp+=2;
		VEC_DBL_INIT(c_tmp, c7);	VEC_DBL_INIT(s_tmp, s7);	c_tmp+=2; s_tmp+=2;

		VEC_DBL_INIT(c_tmp, cA);	VEC_DBL_INIT(s_tmp, sA);	c_tmp+=2; s_tmp+=2;
		VEC_DBL_INIT(c_tmp, c9);	VEC_DBL_INIT(s_tmp, s9);	c_tmp+=2; s_tmp+=2;
		VEC_DBL_INIT(c_tmp, cB);	VEC_DBL_INIT(s_tmp, sB);	c_tmp+=2; s_tmp+=2;

		VEC_DBL_INIT(c_tmp, cE);	VEC_DBL_INIT(s_tmp, sE);	c_tmp+=2; s_tmp+=2;
		VEC_DBL_INIT(c_tmp, cD);	VEC_DBL_INIT(s_tmp, sD);	c_tmp+=2; s_tmp+=2;
		VEC_DBL_INIT(c_tmp, cF);	VEC_DBL_INIT(s_tmp, sF);	c_tmp+=2; s_tmp+=2;

	  #else

		// Sincos data stored in BRed form in SSE2 local-data layout:
		c_tmp = cc0 + 0x02; s_tmp = c_tmp+1;	/* c0,s0 */
		rt = 1.0; it = 0.0;
		VEC_DBL_INIT(c_tmp, rt);	VEC_DBL_INIT(s_tmp, it);	c_tmp+=2; s_tmp+=2;
		VEC_DBL_INIT(c_tmp, c8);	VEC_DBL_INIT(s_tmp, s8);	c_tmp+=2; s_tmp+=2;
		VEC_DBL_INIT(c_tmp, c4);	VEC_DBL_INIT(s_tmp, s4);	c_tmp+=2; s_tmp+=2;
		VEC_DBL_INIT(c_tmp, cC);	VEC_DBL_INIT(s_tmp, sC);	c_tmp+=2; s_tmp+=2;
		VEC_DBL_INIT(c_tmp, c2);	VEC_DBL_INIT(s_tmp, s2);	c_tmp+=2; s_tmp+=2;
		VEC_DBL_INIT(c_tmp, cA);	VEC_DBL_INIT(s_tmp, sA);	c_tmp+=2; s_tmp+=2;
		VEC_DBL_INIT(c_tmp, c6);	VEC_DBL_INIT(s_tmp, s6);	c_tmp+=2; s_tmp+=2;
		VEC_DBL_INIT(c_tmp, cE);	VEC_DBL_INIT(s_tmp, sE);	c_tmp+=2; s_tmp+=2;
		VEC_DBL_INIT(c_tmp, c1);	VEC_DBL_INIT(s_tmp, s1);	c_tmp+=2; s_tmp+=2;
		VEC_DBL_INIT(c_tmp, c9);	VEC_DBL_INIT(s_tmp, s9);	c_tmp+=2; s_tmp+=2;
		VEC_DBL_INIT(c_tmp, c5);	VEC_DBL_INIT(s_tmp, s5);	c_tmp+=2; s_tmp+=2;
		VEC_DBL_INIT(c_tmp, cD);	VEC_DBL_INIT(s_tmp, sD);	c_tmp+=2; s_tmp+=2;
		VEC_DBL_INIT(c_tmp, c3);	VEC_DBL_INIT(s_tmp, s3);	c_tmp+=2; s_tmp+=2;
		VEC_DBL_INIT(c_tmp, cB);	VEC_DBL_INIT(s_tmp, sB);	c_tmp+=2; s_tmp+=2;
		VEC_DBL_INIT(c_tmp, c7);	VEC_DBL_INIT(s_tmp, s7);	c_tmp+=2; s_tmp+=2;
		VEC_DBL_INIT(c_tmp, cF);	VEC_DBL_INIT(s_tmp, sF);	c_tmp+=2; s_tmp+=2;

	  #endif

	#endif

		//******************* Timing loop for Radix-16 DIF transform macro: *******************
		clock1 = getRealTime();
		for(i = 0; i < imax; i++)
		{
			// Copy digits of Pi-data into our vec_dbl inputs:
			for(j1 = 0, j2 = 0; j1 < dim; j1 += stride, j2 += 8)	// j2 is base-index into ran[] input array
			{
		/* The normal index-munging takes way too many cycles in this context, so inline it via 8-way loop unroll:
			#ifdef USE_AVX
				j1 = (j & mask02) + br8[j&7];
			#elif defined(USE_SSE2)
				j1 = (j & mask01) + br4[j&3];
			#else
				j1 = j;
			#endif
				a[j1] = ran[j];
		*/
			#ifdef USE_AVX512	// Set this up so that AVX-512 can use the same ref-data as AVX:
				a[j1   ] = ran[j2  ];	/* Re0 */	a[j1+ 4] = 0;	/* Re4 */
				a[j1+ 1] = ran[j2+2];	/* Re1 */	a[j1+ 5] = 0;	/* Re5 */
				a[j1+ 2] = ran[j2+4];	/* Re2 */	a[j1+ 6] = 0;	/* Re6 */
				a[j1+ 3] = ran[j2+6];	/* Re3 */	a[j1+ 7] = 0;	/* Re7 */
				a[j1+ 8] = ran[j2+1];	/* Im0 */	a[j1+12] = 0;	/* Im4 */
				a[j1+ 9] = ran[j2+3];	/* Im1 */	a[j1+13] = 0;	/* Im5 */
				a[j1+10] = ran[j2+5];	/* Im2 */	a[j1+14] = 0;	/* Im6 */
				a[j1+11] = ran[j2+7];	/* Im3 */	a[j1+15] = 0;	/* Im7 */
			#elif defined(USE_AVX)
				a[j1  ] = ran[j2  ];	// Re0
				a[j1+1] = ran[j2+2];	// Re1
				a[j1+2] = ran[j2+4];	// Re2
				a[j1+3] = ran[j2+6];	// Re3
				a[j1+4] = ran[j2+1];	// Im0
				a[j1+5] = ran[j2+3];	// Im1
				a[j1+6] = ran[j2+5];	// Im2
				a[j1+7] = ran[j2+7];	// Im3
			#elif defined(USE_SSE2)
				a[j1  ] = ran[j2  ];	// Re0
				a[j1+1] = ran[j2+2];	// Re1
				a[j1+2] = ran[j2+1];	// Im0
				a[j1+3] = ran[j2+3];	// Im1
				a[j1+4] = ran[j2+4];	// Re2
				a[j1+5] = ran[j2+6];	// Re3
				a[j1+6] = ran[j2+5];	// Im2
				a[j1+7] = ran[j2+7];	// Im3
			#else
				a[j1  ] = ran[j2  ];	// Re0
				a[j1+1] = ran[j2+1];	// Im0
				a[j1+2] = ran[j2+2];	// Re1
				a[j1+3] = ran[j2+3];	// Im1
				a[j1+4] = ran[j2+4];	// Re2
				a[j1+5] = ran[j2+5];	// Im2
				a[j1+6] = ran[j2+6];	// Re3
				a[j1+7] = ran[j2+7];	// Im3
			#endif
			}

			j1 = 0; j2 = RE_IM_STRIDE;

// 19 May 2016: 10^7-loop timings, 1-threaded on my 2GHz Core2Duo:
//												comments
// no-DFT timing						0.969
// initial timing						2.467	(2.467-0.969)*200 = 300 cycles ... SSE2_RADIX_04_DIF_3TWIDDLE_X2-based (4 calls at 60ish cycles each) should beat that handily
// SSE2_RADIX_04_DIF_3TWIDDLE_X2-based	2.730	ack! [And commenting-out the 4 dft-4 macro calls gives 0.980 sec, insign. different from original no-DFT timing, i.e. the C-code pointer math is not at fault.
// fuse four 4-dft macros into radix-16	2.224	(2.224-0.969)*200 = 250 cycles ... that's more like it!
		#ifdef USE_SSE2

		  #ifdef USE_AVX2	// AVX2/FMA needs tangent-form twiddles:

			vec_dbl *add = (vec_dbl *)a;
			SSE2_RADIX16_DIF_TWIDDLE_1(add,p1,p2,p3,p4,p8,pC,r00,isrt2,pfetch_addr,pfetch_dist);
#if 0
10^6-timing:	setup	+=DIF	DIF-only
avx2:			.208	.380	.172 [224 cycles]
avx512:			.296	.472	.176 [229 cycles]	further fiddling with the adressing parts of this macro -> 0.460, 16 cycles faster!
													[avx512 Tan-twiddles precomp = 140 cycles]
#endif
		  #elif defined(REFACTOR_4DFT_3TWIDDLE)

		   #if 1

			vec_dbl *add = (vec_dbl *)a;
			SSE2_RADIX16_DIF_TWIDDLE_V2(add,p1,p2,p3,p4,p8,pC,r00,two,cc0,pfetch_addr,pfetch_dist);
/*
#ifdef USE_AVX
	dptr = (double *)r00;
	printf("Intermediates:\n");
	for(i = 0; i < 16; i++, dptr += 8) {
		printf("%2u Re.[d0-3] = %16.12f,%16.12f,%16.12f,%16.12f, Im.[d0-3] = %16.12f,%16.12f,%16.12f,%16.12f\n",i,*dptr,*(dptr+1),*(dptr+2),*(dptr+3),*(dptr+4),*(dptr+5),*(dptr+6),*(dptr+7));
	}
	exit(0);
#else
	dptr = (double *)r00;
	printf("Intermediates:\n");
	for(i = 0; i < 16; i++, dptr += 4) {
		printf("%2u Re.[d0,d1] = %16.12f,%16.12f, Im.[d0,d1] = %16.12f,%16.12f\n",i,*dptr,*(dptr+1),*(dptr+2),*(dptr+3));
	}
	exit(0);
	#if 0
		SSE2, Intermediates:
		 0 Re.[d0,d1] =  12.996843798630, 24.930865132335, Im.[d0,d1] =   3.041415000035, 24.061237911217
		 1 Re.[d0,d1] =  -1.967770045489, -3.960148397371, Im.[d0,d1] =  -6.012197708235, -5.001442285342
		 2 Re.[d0,d1] =  -0.996872035772, -2.980018327208, Im.[d0,d1] =  -1.023007259458, -6.018361815583
		 3 Re.[d0,d1] =   1.967798282632, -1.990698407756, Im.[d0,d1] =   7.993789967657, -9.041433810292
		 4 Re.[d0,d1] =  17.937037540455, 15.972338992203, Im.[d0,d1] =  20.044361826623, 14.027564530942
		 5 Re.[d0,d1] =   0.026036426331,  0.015342118897, Im.[d0,d1] =  -1.992277237471, -0.012241272931
		 6 Re.[d0,d1] =   4.013790439910, 12.015332706460, Im.[d0,d1] =   2.010786096064, -1.990767874548
		 7 Re.[d0,d1] = -13.976864406696,  3.996986182440, Im.[d0,d1] =  -8.062870685215,  3.975444616537
		 8 Re.[d0,d1] =  15.947773027030, 21.940126671442, Im.[d0,d1] =  19.038248242231, 25.038232997648
		 9 Re.[d0,d1] =   1.022960222909,  2.004558386725, Im.[d0,d1] =  -1.990745621227, 11.007687462509
		10 Re.[d0,d1] =   6.039898671633,  4.022982548356, Im.[d0,d1] =  -9.001451585837,  3.010797836274
		11 Re.[d0,d1] =  -3.010631921572, -7.967667606523, Im.[d0,d1] =   3.953948964833, -7.056718296431
		12 Re.[d0,d1] =  24.987608157941, 30.970725010455, Im.[d0,d1] =   9.067480879677, 13.070512349057
		13 Re.[d0,d1] =  -8.972346109218,  2.036811854087, Im.[d0,d1] =   0.969395724120, -0.007571096299
		14 Re.[d0,d1] =  -2.999955284040, -4.989179812937, Im.[d0,d1] =   0.981587603769, -1.039861018570
		15 Re.[d0,d1] =  -1.015306764683,  3.981642948395, Im.[d0,d1] =   0.981535792434, -0.023080234188
	#endif
#endif
*/
		   #else

			/* Pass 1: */
			j = p2*8;	// Set j to the bytewise address offset between the pointers to the 1st,2nd 4-DFT's inputs
			k = 0x80;	// Set k to the bytewise address offset between the pointers to the 1st,2nd 4-DFT's outputs
			// Pass1-DFTs all use same twiddle-triplet, so lit-byte address offset between ptrs to 1st,2nd DFT's twiddles = 0
			i0 = (vec_dbl *)a; i1 = (vec_dbl *)(a+p4); i2 = (vec_dbl *)(a+p8); i3 = (vec_dbl *)(a+pC);
			o0 = r00; o1 = r00+2; o2 = r00+4; o3 = r00+6;
			c_tmp = cc0;	/* c8,4,C */
			SSE2_RADIX_04_DIF_3TWIDDLE_X2(i0,i1,i2,i3,j, two,c_tmp,   0, o0,o1,o2,o3,k)

			i0 = (vec_dbl *)(a+p1); i1 = (vec_dbl *)(a+p5); i2 = (vec_dbl *)(a+p9); i3 = (vec_dbl *)(a+pD);
			o0 += 16; o1 += 16; o2 += 16; o3 += 16;
			SSE2_RADIX_04_DIF_3TWIDDLE_X2(i0,i1,i2,i3,j, two,c_tmp,   0, o0,o1,o2,o3,k)

			/* Pass 2: */
			j = 0x40;	// Set j to the bytewise address offset between the pointers to the 1st,2nd 4-DFT's inputs
			k = p4*8;	// Set k to the bytewise address offset between the pointers to the 1st,2nd 4-DFT's outputs
			// 0x60 = literal-bytewise address offset between the pointers to the first and second DFT's twiddles = 0:
			i0 = r00; i1 = r00+16; i2 = r00+8; i3 = r00+24;
			o0 = (vec_dbl *)a; o1 = (vec_dbl *)(a+p2); o2 = (vec_dbl *)(a+p1); o3 = (vec_dbl *)(a+p3);
			c_tmp += 6;	/* c2,1,3 */
			SSE2_RADIX_04_DIF_3TWIDDLE_X2(i0,i1,i2,i3,j, two,c_tmp,0x60, o0,o1,o2,o3,k)

			i0 += 2; i1 += 2; i2 += 2; i3 += 2;
			o0 = (vec_dbl *)(a+p8); o1 = (vec_dbl *)(a+pA); o2 = (vec_dbl *)(a+p9); o3 = (vec_dbl *)(a+pB);
			c_tmp += 12;	/* cA,9,B */
			SSE2_RADIX_04_DIF_3TWIDDLE_X2(i0,i1,i2,i3,j, two,c_tmp,0x60, o0,o1,o2,o3,k)
		   #endif

		  #else	// REFACTOR_4DFT_3TWIDDLE = false:

			vec_dbl *add = (vec_dbl *)a;
			SSE2_RADIX16_DIF_TWIDDLE(add,p1,p2,p3,p4,p8,pC,r00,isrt2,pfetch_addr,pfetch_dist);

		  #endif

		#else	// USE_SSE2 = false:

		  #ifdef REFACTOR_4DFT_3TWIDDLE
			/*
			Pass 1:
			y0-3 = radix_4_3twid(x0,4,8,c; w4,8,c)
			y4-7 = radix_4_3twid(x2,6,a,e; w4,8,c)
			y8-b = radix_4_3twid(x1,5,9,d; w4,8,c)
			yc-f = radix_4_3twid(x3,7,b,f; w4,8,c)
			*/
		  	RADIX_04_DIF_3TWIDDLE(a[j1   ],a[j2   ],a[j1+p4],a[j2+p4],a[j1+p8],a[j2+p8],a[j1+pC],a[j2+pC], t0 ,t1 ,t2 ,t3 ,t4 ,t5 ,t6 ,t7 , c4,s4,c8,s8,cC,sC, rt,it);
		  	RADIX_04_DIF_3TWIDDLE(a[j1+p2],a[j2+p2],a[j1+p6],a[j2+p6],a[j1+pA],a[j2+pA],a[j1+pE],a[j2+pE], t8, t9 ,t10,t11,t12,t13,t14,t15, c4,s4,c8,s8,cC,sC, rt,it);
		  	RADIX_04_DIF_3TWIDDLE(a[j1+p1],a[j2+p1],a[j1+p5],a[j2+p5],a[j1+p9],a[j2+p9],a[j1+pD],a[j2+pD], t16,t17,t18,t19,t20,t21,t22,t23, c4,s4,c8,s8,cC,sC, rt,it);
		  	RADIX_04_DIF_3TWIDDLE(a[j1+p3],a[j2+p3],a[j1+p7],a[j2+p7],a[j1+pB],a[j2+pB],a[j1+pF],a[j2+pF], t24,t25,t26,t27,t28,t29,t30,t31, c4,s4,c8,s8,cC,sC, rt,it);

			/*
			Pass 2:
			z0-3 = radix_4_3twid(y0,8,4,c; w1    ,w2      ,w3      )
			z4-7 = radix_4_3twid(y2,a,6,e; w1.E^2,w2.I    ,w3.I.E^2)
			z8-b = radix_4_3twid(y1,9,5,d; w1.E^1,w2.  E^2,w3.  E^3)
			zc-f = radix_4_3twid(y3,b,7,f; w1.E^3,w2.I.E^2,w3.  E  )
			*/
		  	RADIX_04_DIF_3TWIDDLE(t0 ,t1 ,t16,t17,t8 ,t9 ,t24,t25, a[j1   ],a[j2   ],a[j1+p2],a[j2+p2],a[j1+p1],a[j2+p1],a[j1+p3],a[j2+p3], c1,s1, c2,s2, c3,s3, rt,it);
		  	RADIX_04_DIF_3TWIDDLE(t4 ,t5 ,t20,t21,t12,t13,t28,t29, a[j1+p4],a[j2+p4],a[j1+p6],a[j2+p6],a[j1+p5],a[j2+p5],a[j1+p7],a[j2+p7], c5,s5, c6,s6, c7,s7, rt,it);
		  	RADIX_04_DIF_3TWIDDLE(t2 ,t3 ,t18,t19,t10,t11,t26,t27, a[j1+p8],a[j2+p8],a[j1+pA],a[j2+pA],a[j1+p9],a[j2+p9],a[j1+pB],a[j2+pB], c9,s9, cA,sA, cB,sB, rt,it);
		  	RADIX_04_DIF_3TWIDDLE(t6 ,t7 ,t22,t23,t14,t15,t30,t31, a[j1+pC],a[j2+pC],a[j1+pE],a[j2+pE],a[j1+pD],a[j2+pD],a[j1+pF],a[j2+pF], cD,sD, cE,sE, cF,sF, rt,it);

		  #else

			RADIX_16_DIF_TWIDDLE(
				a[j1   ],a[j2   ],a[j1+p1],a[j2+p1],a[j1+p2],a[j2+p2],a[j1+p3],a[j2+p3],a[j1+p4],a[j2+p4],a[j1+p4+p1],a[j2+p4+p1],a[j1+p4+p2],a[j2+p4+p2],a[j1+p4+p3],a[j2+p4+p3],a[j1+p8],a[j2+p8],a[j1+p8+p1],a[j2+p8+p1],a[j1+p8+p2],a[j2+p8+p2],a[j1+p8+p3],a[j2+p8+p3],a[j1+pC],a[j2+pC],a[j1+pC+p1],a[j2+pC+p1],a[j1+pC+p2],a[j2+pC+p2],a[j1+pC+p3],a[j2+pC+p3],
				c1,s1,c2,s2,c3,s3,c4,s4,c5,s5,c6,s6,c7,s7,c8,s8,c9,s9,cA,sA,cB,sB,cC,sC,cD,sD,cE,sE,cF,sF,
				c,s);

		  #endif
/*
printf("DIF Outputs =\n0  [%16.12f,%16.12f]\n1  [%16.12f,%16.12f]\n2  [%16.12f,%16.12f]\n3  [%16.12f,%16.12f]\n4  [%16.12f,%16.12f]\n5  [%16.12f,%16.12f]\n6  [%16.12f,%16.12f]\n7  [%16.12f,%16.12f]\n8  [%16.12f,%16.12f]\n9  [%16.12f,%16.12f]\n10 [%16.12f,%16.12f]\n11 [%16.12f,%16.12f]\n12 [%16.12f,%16.12f]\n13 [%16.12f,%16.12f]\n14 [%16.12f,%16.12f]\n15 [%16.12f,%16.12f]\n",
a[j1   ],a[j2   ],a[j1+p1],a[j2+p1],a[j1+p2],a[j2+p2],a[j1+p3],a[j2+p3],a[j1+p4],a[j2+p4],a[j1+p4+p1],a[j2+p4+p1],a[j1+p4+p2],a[j2+p4+p2],a[j1+p4+p3],a[j2+p4+p3],a[j1+p8],a[j2+p8],a[j1+p8+p1],a[j2+p8+p1],a[j1+p8+p2],a[j2+p8+p2],a[j1+p8+p3],a[j2+p8+p3],a[j1+pC],a[j2+pC],a[j1+pC+p1],a[j2+pC+p1],a[j1+pC+p2],a[j2+pC+p2],a[j1+pC+p3],a[j2+pC+p3]);
exit(0);
*/
		#endif
		}
		clock2 = getRealTime();
		tdiff = (double)(clock2 - clock1);
		printf("%s: Time for %u DIF macro calls =%s [tdiff = %20.10e]\n",func, imax, get_time_str(tdiff), tdiff);

		// Check outputs vs ref-data:
		nerr = 0;	dtmp = avg_err = 0.0;
		for(j1 = 0, j2 = 0; j1 < dim; j1 += stride, j2 += 8)	// j2 is base-index into ref-array
		{
			j = j1+RE_IM_STRIDE;
		#ifdef USE_AVX	// Since we set up AVX-512 mode to only use nonzero data in lower 4 double-slots of each 8-vector, can use same code here:
		//	printf("Out[%2u] Re.[d0-3] = %16.12f,%16.12f,%16.12f,%16.12f, Im.[d0-3] = %16.12f,%16.12f,%16.12f,%16.12f\n",j1/p1,a[j1],a[j1+1],a[j1+2],a[j1+3],a[j],a[j+1],a[j+2],a[j+3]);
		//	printf("Ref[%2u] Re.[d0-3] = %16.12f,%16.12f,%16.12f,%16.12f, Im.[d0-3] = %16.12f,%16.12f,%16.12f,%16.12f\n",j1/p1,ref1[j2],ref1[j2+2],ref1[j2+4],ref1[j2+6],ref1[j2+1],ref1[j2+3],ref1[j2+5],ref1[j2+7]);
			dtmp = fabs(a[j1  ] - ref1[j2  ]); avg_err += dtmp; if(dtmp > 1e-10){ /*printf("error Re.d0\n");*/ nerr++; };
			dtmp = fabs(a[j1+1] - ref1[j2+2]); avg_err += dtmp; if(dtmp > 1e-10){ /*printf("error Re.d1\n");*/ nerr++; };
			dtmp = fabs(a[j1+2] - ref1[j2+4]); avg_err += dtmp; if(dtmp > 1e-10){ /*printf("error Re.d2\n");*/ nerr++; };
			dtmp = fabs(a[j1+3] - ref1[j2+6]); avg_err += dtmp; if(dtmp > 1e-10){ /*printf("error Re.d3\n");*/ nerr++; };
			dtmp = fabs(a[j   ] - ref1[j2+1]); avg_err += dtmp; if(dtmp > 1e-10){ /*printf("error Im.d0\n");*/ nerr++; };
			dtmp = fabs(a[j +1] - ref1[j2+3]); avg_err += dtmp; if(dtmp > 1e-10){ /*printf("error Im.d1\n");*/ nerr++; };
			dtmp = fabs(a[j +2] - ref1[j2+5]); avg_err += dtmp; if(dtmp > 1e-10){ /*printf("error Im.d2\n");*/ nerr++; };
			dtmp = fabs(a[j +3] - ref1[j2+7]); avg_err += dtmp; if(dtmp > 1e-10){ /*printf("error Im.d3\n");*/ nerr++; };
		#elif defined(USE_SSE2)
		//	printf("Out%2u Re.[d0,d1] = %16.12f,%16.12f, Im.[d0,d1] = %16.12f,%16.12f\n",j1/p1,a[j1],a[j1+1],a[j1+2],a[j1+3]);
		//	printf("Ref%2u Re.[d0,d1] = %16.12f,%16.12f, Im.[d0,d1] = %16.12f,%16.12f\n",j1/p1,ref1[j2],ref1[j2+2],ref1[j2+1],ref1[j2+3]);
			dtmp = fabs(a[j1  ] - ref1[j2  ]); avg_err += dtmp; if(dtmp > 1e-10){ /*printf("error Re.d0\n");*/ nerr++; };
			dtmp = fabs(a[j1+1] - ref1[j2+2]); avg_err += dtmp; if(dtmp > 1e-10){ /*printf("error Re.d1\n");*/ nerr++; };
			dtmp = fabs(a[j1+2] - ref1[j2+1]); avg_err += dtmp; if(dtmp > 1e-10){ /*printf("error Im.d0\n");*/ nerr++; };
			dtmp = fabs(a[j1+3] - ref1[j2+3]); avg_err += dtmp; if(dtmp > 1e-10){ /*printf("error Im.d1\n");*/ nerr++; };
		#else
		//	printf("Out%2u Re,Im = %16.12f,%16.12f\n",j1/p1,a[j1],a[j1+1]);
			dtmp = fabs(a[j1  ] - ref1[j2  ]); avg_err += dtmp; if(dtmp > 1e-10){ /*printf("error Re.d0\n");*/ nerr++; };
			dtmp = fabs(a[j1+1] - ref1[j2+1]); avg_err += dtmp; if(dtmp > 1e-10){ /*printf("error Im.d0\n");*/ nerr++; };
		#endif
		}
printf("DIF: nerr = %u, ",nerr);
		ASSERT(HERE, nerr == 0, "DIF Outputs mismatch ref-data!");
		printf("\tSummed roundoff error = %20.10e]\n",avg_err);

		//******************* Timing loop for Radix-16 DIT transform macro: *******************
	  #ifdef USE_AVX2	// AVX2/FMA needs tangent-form twiddles:
		// In AVX2/FMA mode, since we need to replace most of the raw sincos data with derived ones,
		// simply place one copy of each computed double in a double-sized slot of the local memory.
		// We will be using AVX2/FMA-based Newtonian iterative inversion on the 16 doubles whose
		// multiplicative inverse is needed (the real part of the basic root of unity c and of the
		// 15 complex twiddles, c1-15), so store those in packed form in 4 AVX-register-sized
		// contiguous memory locations, and the others in a separate chunk of memory. After the
		// vector-iterative inversion we'll need to combine the 2 sets of data and place (in quadruplicate)
		// into their final SIMD-suitable memory slots.

	  clock1 = getRealTime();
	  for(i = 0; i < imax; i++) {	// repeater loop
		add0 = (double *)cc0;	// add0 points to 16 cos-data-to-be-inverted; Need a double-ptr on lhs here
		add1 = add0 + 16;		// add1 points to block of memory temporarily used to store the corresponding sine data
		*add0++ = c;	// Since tan0 defined as const, we could init these directly, but init with c0,s0 anyway
		*add1++ = s;	// and use result as a check onthe accuracy of the FMA-based Newton iterative inversion.

		*add0++ = c1;	// c1, for inversion
		*add1++ = s1;	// s1  slot will hold __r1 = s1 /c1
		*add0++ = c2;	// c2, for inversion
		*add1++ = s2;	// s2  slot will hold __r2 = s2 /c2
		*add0++ = c3;	// c3, for inversion
		*add1++ = s3;	// s3  slot will hold __r3 = s3 /c3
		*add0++ = c4;	// c4, for inversion
		*add1++ = s4;	// s4  slot will hold __r4 = s4 /c4
		*add0++ = c5;	// c5, for inversion
		*add1++ = s5;	// s5  slot will hold __r5 = s5 /c5
		*add0++ = c6;	// c6, for inversion
		*add1++ = s6;	// s6  slot will hold __r6 = s6 /c6
		*add0++ = c7;	// c7, for inversion
		*add1++ = s7;	// s7  slot will hold __r7 = s7 /c7
		*add0++ = c8;	// c8, for inversion
		*add1++ = s8;	// s8  slot will hold __r8 = s8 /c8
		*add0++ = c9;	// c9, for inversion
		*add1++ = s9;	// s9  slot will hold __r9 = s9 /c9
		*add0++ = cA;	// c10, for inversion
		*add1++ = sA;	// s10 slot will hold __rA = s10/c10
		*add0++ = cB;	// c11, for inversion
		*add1++ = sB;	// s11 slot will hold __rB = s11/c11
		*add0++ = cC;	// c12, for inversion
		*add1++ = sC;	// s12 slot will hold __rC = s12/c12
		*add0++ = cD;	// c13, for inversion
		*add1++ = sD;	// s13 slot will hold __rD = s13/c13
		*add0++ = cE;	// c14, for inversion
		*add1++ = sE;	// s14 slot will hold __rE = s14/c14
		*add0++ = cF;	// c15, for inversion
		*add1++ = sF;	// s15 slot will hold __rF = s15/c15
		/*
		At this point, the 8 ymm-sized [32-byte] chunks starting at &cc0 contain the following scalar-double data:

		0:	c,c1-3		4:	s,s1-3
		1:	c4-7		5:	s4-7
		2:	c8-B		6:	s8-B
		3:	cC-F		7:	sC-F
		*/

		// Now send the cosine terms to the inversion routine, which also does the combine-and-populate-SIMD-slots step.
		RADIX16_COMPUTE_FMA_SINCOS_DIT(cc0,two);
	  }	// repeater loop
		clock2 = getRealTime();
		tdiff = (double)(clock2 - clock1);
		printf("%s: Time for %u FMA DIT tan-twiddles setup calls =%s [tdiff = %20.10e]\n",func, imax, get_time_str(tdiff), tdiff);

	  #endif

		clock1 = getRealTime();
		for(i = 0; i < imax; i++)
		{
			// Copy digits of Pi-data into our vec_dbl inputs:
			for(j1 = 0, j2 = 0; j1 < dim; j1 += stride, j2 += 8)	// j2 is base-index into ran[] input array
			{
		/* The normal index-munging takes way too many cycles in this context, so inline it via 8-way loop unroll:
			#ifdef USE_AVX
				j1 = (j & mask02) + br8[j&7];
			#elif defined(USE_SSE2)
				j1 = (j & mask01) + br4[j&3];
			#else
				j1 = j;
			#endif
				a[j1] = ran[j];
		*/
			#ifdef USE_AVX512	// Set this up so that AVX-512 can use the same ref-data as AVX:
				a[j1   ] = ran[j2  ];	/* Re0 */	a[j1+ 4] = 0;	/* Re4 */
				a[j1+ 1] = ran[j2+2];	/* Re1 */	a[j1+ 5] = 0;	/* Re5 */
				a[j1+ 2] = ran[j2+4];	/* Re2 */	a[j1+ 6] = 0;	/* Re6 */
				a[j1+ 3] = ran[j2+6];	/* Re3 */	a[j1+ 7] = 0;	/* Re7 */
				a[j1+ 8] = ran[j2+1];	/* Im0 */	a[j1+12] = 0;	/* Im4 */
				a[j1+ 9] = ran[j2+3];	/* Im1 */	a[j1+13] = 0;	/* Im5 */
				a[j1+10] = ran[j2+5];	/* Im2 */	a[j1+14] = 0;	/* Im6 */
				a[j1+11] = ran[j2+7];	/* Im3 */	a[j1+15] = 0;	/* Im7 */
			#elif defined(USE_AVX)
				a[j1  ] = ran[j2  ];	// Re0
				a[j1+1] = ran[j2+2];	// Re1
				a[j1+2] = ran[j2+4];	// Re2
				a[j1+3] = ran[j2+6];	// Re3
				a[j1+4] = ran[j2+1];	// Im0
				a[j1+5] = ran[j2+3];	// Im1
				a[j1+6] = ran[j2+5];	// Im2
				a[j1+7] = ran[j2+7];	// Im3
			#elif defined(USE_SSE2)
				a[j1  ] = ran[j2  ];	// Re0
				a[j1+1] = ran[j2+2];	// Re1
				a[j1+2] = ran[j2+1];	// Im0
				a[j1+3] = ran[j2+3];	// Im1
				a[j1+4] = ran[j2+4];	// Re2
				a[j1+5] = ran[j2+6];	// Re3
				a[j1+6] = ran[j2+5];	// Im2
				a[j1+7] = ran[j2+7];	// Im3
			#else
				a[j1  ] = ran[j2  ];	// Re0
				a[j1+1] = ran[j2+1];	// Im0
				a[j1+2] = ran[j2+2];	// Re1
				a[j1+3] = ran[j2+3];	// Im1
				a[j1+4] = ran[j2+4];	// Re2
				a[j1+5] = ran[j2+5];	// Im2
				a[j1+6] = ran[j2+6];	// Re3
				a[j1+7] = ran[j2+7];	// Im3
			#endif
			}

			j1 = 0; j2 = RE_IM_STRIDE;

// 23 May 2016: 10^7-loop timings, 1-threaded on my 2GHz Core2Duo:
//												comments
// no-DFT timing						0.982
// 1x 4-DFT 22 addpd 20 mulpd 31 ld/st	1.202	(1.202-0.982)*200 = 44 cycles [24 cycles for opening 16 addpd, 20 for 3 twiddle-cmul (12 mulpd, 6 addpd)]
// 2x 4-DFT, 2-column side-by-side opt	1.312	66 cycles
// current-code timing					2.422	(2.422-0.978)*200 = 288 cycles
// fuse four 4-dft macros into radix-16	2.233	250 cycles ... still 50 cycles to go, but nice speedup nonetheless.
		#ifdef USE_SSE2

		  #ifdef USE_AVX2	// AVX2/FMA needs tangent-form twiddles:

			vec_dbl *add = (vec_dbl *)a;
			SSE2_RADIX16_DIT_TWIDDLE_1(add,p1,p2,p3,p4,p8,pC,r00,cc0,pfetch_addr,pfetch_dist);
#if 0
10^6-timing:	setup	+=DIT	DIT-only
avx2:			.208	.398	.190 [247 cycles]
avx512:			.296	.489	.193 [251 cycles]	further fiddling with the adressing parts of this macro -> 0.476, 17 cycles faster!
													[avx512 Tan-twiddles precomp = 122 cycles]
#endif
		  #elif defined(REFACTOR_4DFT_3TWIDDLE)

		   #if 1

			vec_dbl *add = (vec_dbl *)a;
			SSE2_RADIX16_DIT_TWIDDLE_V2(add,p1,p2,p3,p4,p8,pC,r00,two,cc0,pfetch_addr,pfetch_dist);

		   #else

			/*
			Pass 1:
			y0-3 = radix_4dit_3twid(x0-3; w1    ,w2    ,w3    )
			y4-7 = radix_4dit_3twid(x4-7; w1.E^2,w2.E^4,w3.E^6)
			y8-b = radix_4dit_3twid(x8-b; w1.E^1,w2.E^2,w3.E^3)
			yc-f = radix_4dit_3twid(xc-f; w1.E^3,w2.E^6,w3.E^9)
			*/
			j = p4*8;	// Set j to the bytewise address offset between the pointers to the 1st,2nd 4-DFT's inputs
			k = 0x80;	// Set k to the bytewise address offset between the pointers to the 1st,2nd 4-DFT's outputs
			// 0x60 = literal-bytewise address offset between the pointers to the first and second DFT's twiddles:
			i0 = (vec_dbl *)a; i1 = (vec_dbl *)(a+p1); i2 = (vec_dbl *)(a+p2); i3 = (vec_dbl *)(a+p3);
			o0 = r00; o1 = r00+2; o2 = r00+4; o3 = r00+6;
			c_tmp = cc0+6;	// c2,1,3
			SSE2_RADIX_04_DIT_3TWIDDLE_X2(i0,i1,i2,i3,j, two,c_tmp,0x60, o0,o1,o2,o3,k)

			i0 = (vec_dbl *)(a+p8); i1 = (vec_dbl *)(a+p9); i2 = (vec_dbl *)(a+pA); i3 = (vec_dbl *)(a+pB);
			o0 += 16; o1 += 16; o2 += 16; o3 += 16;
			c_tmp += 12;		// cA,9,B
			SSE2_RADIX_04_DIT_3TWIDDLE_X2(i0,i1,i2,i3,j, two,c_tmp,0x60, o0,o1,o2,o3,k)

			/*
			Pass 2:
			z0,4,8,c = radix_4dit_3twid(y0,4,8,c; w4,8,c)
			z2,6,a,e = radix_4dit_3twid(y2,6,a,e; w4,8,c)
			z1,5,9,d = radix_4dit_3twid(y1,5,9,d; w4,8,c)
			z3,7,b,f = radix_4dit_3twid(y3,7,b,f; w4,8,c)
			*/
			j = 0x40;	// Set j to the bytewise address offset between the pointers to the 1st,2nd 4-DFT's inputs
			k = p2*8;	// Set k to the bytewise address offset between the pointers to the 1st,2nd 4-DFT's outputs
			c_tmp = cc0;	// c8,4,C
			// Pass1-DFTs all use same twiddle-triplet, so lit-byte address offset between ptrs to 1st,2nd DFT's twiddles = 0
			i0 = r00; i1 = r00+8; i2 = r00+16; i3 = r00+24;
			o0 = (vec_dbl *)a; o1 = (vec_dbl *)(a+p4); o2 = (vec_dbl *)(a+p8); o3 = (vec_dbl *)(a+pC);
			SSE2_RADIX_04_DIT_3TWIDDLE_X2(i0,i1,i2,i3,j, two,c_tmp,   0, o0,o1,o2,o3,k)

			i0 += 2; i1 += 2; i2 += 2; i3 += 2;
			o0 = (vec_dbl *)(a+p1); o1 = (vec_dbl *)(a+p5); o2 = (vec_dbl *)(a+p9); o3 = (vec_dbl *)(a+pD);
			SSE2_RADIX_04_DIT_3TWIDDLE_X2(i0,i1,i2,i3,j, two,c_tmp,   0, o0,o1,o2,o3,k)

		   #endif

		  #else

			vec_dbl *add = (vec_dbl *)a;
			SSE2_RADIX16_DIT_TWIDDLE(add,p1,p2,p3,p4,p8,r00,isrt2,pfetch_addr,pfetch_dist);

		  #endif
/*
	dptr = (double *)r00;
	printf("Intermediates:\n");
	for(i = 0; i < 16; i++, dptr += 4) {
		printf("%2u Re.[d0,d1] = %16.12f,%16.12f, Im.[d0,d1] = %16.12f,%16.12f\n",i,*dptr,*(dptr+1),*(dptr+2),*(dptr+3));
	}
	exit(0);
	#if 0
	SSE2, Intermediates:
	 0 Re.[d0,d1] =  13.000000000000, 25.000000000000, Im.[d0,d1] =  10.000000000000, 16.000000000000
	 1 Re.[d0,d1] =  -2.000383348119, -0.002684466313, Im.[d0,d1] =  -0.999232936091, -6.999999485260
	 2 Re.[d0,d1] =   2.998465136951, -6.998463960403, Im.[d0,d1] =  -2.002300382682,  2.005368343957
	 3 Re.[d0,d1] =  -2.003450132394, -2.008052073743, Im.[d0,d1] =  -2.997697043900, -6.997694396666
	 4 Re.[d0,d1] =  13.000000000000, 26.000000000000, Im.[d0,d1] =  13.000000000000, 20.000000000000
	 5 Re.[d0,d1] =  -4.242098031044,  9.903290617327, Im.[d0,d1] =   1.415840490666,  9.895697799992
	 6 Re.[d0,d1] =   3.008436011095,  1.993863489176, Im.[d0,d1] =  10.997695793535, -8.001531627541
	 7 Re.[d0,d1] =  -5.653596441804, -1.409331530533, Im.[d0,d1] =   2.834933380737,  4.244264911271
	 8 Re.[d0,d1] =  26.000000000000, 26.000000000000, Im.[d0,d1] =  12.000000000000, 19.000000000000
	 9 Re.[d0,d1] =   1.622084999221, -1.464427707060, Im.[d0,d1] =  -3.920311244697,  1.689808122478
	10 Re.[d0,d1] = -11.313705171203,  9.195097171880, Im.[d0,d1] =   0.008677504888,  3.528482393280
	11 Re.[d0,d1] =  -4.358901615977,  2.772957578273, Im.[d0,d1] =   7.937252465573,  1.144860807740
	12 Re.[d0,d1] =  20.000000000000, 17.000000000000, Im.[d0,d1] =  16.000000000000, 21.000000000000
	13 Re.[d0,d1] =   1.433541444084,  7.838131936196, Im.[d0,d1] = -11.311275742731, -3.250182725754
	14 Re.[d0,d1] =  -1.415297834511,  9.897322648581, Im.[d0,d1] =  -1.413128458289, -2.836019109578
	15 Re.[d0,d1] =   9.605892403457, -2.470663184755, Im.[d0,d1] =  10.085971997443,  9.689985728963
	#endif
//
	dptr = (double *)a;
	printf("DIT Outputs:\n");
	for(i = 0; i < 16; i++, dptr += 4) {
		printf("%2u Re.[d0,d1] = %16.12f,%16.12f, Im.[d0,d1] = %16.12f,%16.12f\n",i,*dptr,*(dptr+1),*(dptr+2),*(dptr+3));
	}
	exit(0);
*/
		#else

		  #ifdef REFACTOR_4DFT_3TWIDDLE
			/*
			Pass 1:
			y0-3 = radix_4dit_3twid(x0-3; w1    ,w2    ,w3    )
			y4-7 = radix_4dit_3twid(x4-7; w1.E^2,w2.E^4,w3.E^6)
			y8-b = radix_4dit_3twid(x8-b; w1.E^1,w2.E^2,w3.E^3)
			yc-f = radix_4dit_3twid(xc-f; w1.E^3,w2.E^6,w3.E^9)
			*/
		  	RADIX_04_DIT_3TWIDDLE(a[j1   ],a[j2   ],a[j1+p1],a[j2+p1],a[j1+p2],a[j2+p2],a[j1+p3],a[j2+p3], t0 ,t1 ,t2 ,t3 ,t4 ,t5 ,t6 ,t7 , c1,s1, c2,s2, c3,s3, rt,it);
		  	RADIX_04_DIT_3TWIDDLE(a[j1+p4],a[j2+p4],a[j1+p5],a[j2+p5],a[j1+p6],a[j2+p6],a[j1+p7],a[j2+p7], t8, t9 ,t10,t11,t12,t13,t14,t15, c5,s5, c6,s6, c7,s7, rt,it);
		  	RADIX_04_DIT_3TWIDDLE(a[j1+p8],a[j2+p8],a[j1+p9],a[j2+p9],a[j1+pA],a[j2+pA],a[j1+pB],a[j2+pB], t16,t17,t18,t19,t20,t21,t22,t23, c9,s9, cA,sA, cB,sB, rt,it);
		  	RADIX_04_DIT_3TWIDDLE(a[j1+pC],a[j2+pC],a[j1+pD],a[j2+pD],a[j1+pE],a[j2+pE],a[j1+pF],a[j2+pF], t24,t25,t26,t27,t28,t29,t30,t31, cD,sD, cE,sE, cF,sF, rt,it);
/*
printf("DIT Midputs:\n 0 [%16.12f,%16.12f]\n 1 [%16.12f,%16.12f]\n 2 [%16.12f,%16.12f]\n 3 [%16.12f,%16.12f]\n 4 [%16.12f,%16.12f]\n 5 [%16.12f,%16.12f]\n 6 [%16.12f,%16.12f]\n 7 [%16.12f,%16.12f]\n 8 [%16.12f,%16.12f]\n 9 [%16.12f,%16.12f]\n10 [%16.12f,%16.12f]\n11 [%16.12f,%16.12f]\n12 [%16.12f,%16.12f]\n13 [%16.12f,%16.12f]\n14 [%16.12f,%16.12f]\n15 [%16.12f,%16.12f]\n",
t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31);
exit(0);
*/
			/*
			Pass 2:
			z0,4,8,c = radix_4dit_3twid(y0,4,8,c; w4,8,c)
			z2,6,a,e = radix_4dit_3twid(y2,6,a,e; w4,8,c)
			z1,5,9,d = radix_4dit_3twid(y1,5,9,d; w4,8,c)
			z3,7,b,f = radix_4dit_3twid(y3,7,b,f; w4,8,c)
			*/
		  	RADIX_04_DIT_3TWIDDLE(t0 ,t1 ,t8 ,t9 ,t16,t17,t24,t25, a[j1   ],a[j2   ],a[j1+p4],a[j2+p4],a[j1+p8],a[j2+p8],a[j1+pC],a[j2+pC], c4,s4, c8,s8, cC,sC, rt,it);
		  	RADIX_04_DIT_3TWIDDLE(t4 ,t5 ,t12,t13,t20,t21,t28,t29, a[j1+p2],a[j2+p2],a[j1+p6],a[j2+p6],a[j1+pA],a[j2+pA],a[j1+pE],a[j2+pE], c4,s4, c8,s8, cC,sC, rt,it);
		  	RADIX_04_DIT_3TWIDDLE(t2 ,t3 ,t10,t11,t18,t19,t26,t27, a[j1+p1],a[j2+p1],a[j1+p5],a[j2+p5],a[j1+p9],a[j2+p9],a[j1+pD],a[j2+pD], c4,s4, c8,s8, cC,sC, rt,it);
		  	RADIX_04_DIT_3TWIDDLE(t6 ,t7 ,t14,t15,t22,t23,t30,t31, a[j1+p3],a[j2+p3],a[j1+p7],a[j2+p7],a[j1+pB],a[j2+pB],a[j1+pF],a[j2+pF], c4,s4, c8,s8, cC,sC, rt,it);

		  #else

			RADIX_16_DIT_TWIDDLE(
				a[j1   ],a[j2   ],a[j1+p1],a[j2+p1],a[j1+p2],a[j2+p2],a[j1+p3],a[j2+p3],a[j1+p4],a[j2+p4],a[j1+p4+p1],a[j2+p4+p1],a[j1+p4+p2],a[j2+p4+p2],a[j1+p4+p3],a[j2+p4+p3],a[j1+p8],a[j2+p8],a[j1+p8+p1],a[j2+p8+p1],a[j1+p8+p2],a[j2+p8+p2],a[j1+p8+p3],a[j2+p8+p3],a[j1+pC],a[j2+pC],a[j1+pC+p1],a[j2+pC+p1],a[j1+pC+p2],a[j2+pC+p2],a[j1+pC+p3],a[j2+pC+p3],
				c1,s1,c2,s2,c3,s3,c4,s4,c5,s5,c6,s6,c7,s7,c8,s8,c9,s9,cA,sA,cB,sB,cC,sC,cD,sD,cE,sE,cF,sF,
				c,s);

		  #endif
/*
printf("DIT Outputs:\n 0 [%16.12f,%16.12f]\n 1 [%16.12f,%16.12f]\n 2 [%16.12f,%16.12f]\n 3 [%16.12f,%16.12f]\n 4 [%16.12f,%16.12f]\n 5 [%16.12f,%16.12f]\n 6 [%16.12f,%16.12f]\n 7 [%16.12f,%16.12f]\n 8 [%16.12f,%16.12f]\n 9 [%16.12f,%16.12f]\n10 [%16.12f,%16.12f]\n11 [%16.12f,%16.12f]\n12 [%16.12f,%16.12f]\n13 [%16.12f,%16.12f]\n14 [%16.12f,%16.12f]\n15 [%16.12f,%16.12f]\n",
a[j1   ],a[j2   ],a[j1+p1],a[j2+p1],a[j1+p2],a[j2+p2],a[j1+p3],a[j2+p3],a[j1+p4],a[j2+p4],a[j1+p4+p1],a[j2+p4+p1],a[j1+p4+p2],a[j2+p4+p2],a[j1+p4+p3],a[j2+p4+p3],a[j1+p8],a[j2+p8],a[j1+p8+p1],a[j2+p8+p1],a[j1+p8+p2],a[j2+p8+p2],a[j1+p8+p3],a[j2+p8+p3],a[j1+pC],a[j2+pC],a[j1+pC+p1],a[j2+pC+p1],a[j1+pC+p2],a[j2+pC+p2],a[j1+pC+p3],a[j2+pC+p3]);
exit(0);
*/
		#endif
		}
		clock2 = getRealTime();
		tdiff = (double)(clock2 - clock1);
		printf("%s: Time for %u DIT macro calls =%s [tdiff = %20.10e]\n",func, imax, get_time_str(tdiff), tdiff);
		// Check outputs vs ref-data:
		nerr = 0;
		for(j1 = 0, j2 = 0; j1 < dim; j1 += stride, j2 += 8)	// j2 is base-index into ref-array
		{
			j = j1+RE_IM_STRIDE;
		#ifdef USE_AVX	// Since we set up AVX-512 mode to only use nonzero data in lower 4 double-slots of each 8-vector, can use same code here:
		//	printf("Out[%2u] Re.[d0-3] = %16.12f,%16.12f,%16.12f,%16.12f, Im.[d0-3] = %16.12f,%16.12f,%16.12f,%16.12f\n",j1/p1,a[j1],a[j1+1],a[j1+2],a[j1+3],a[j],a[j+1],a[j+2],a[j+3]);
		//	printf("Ref[%2u] Re.[d0-3] = %16.12f,%16.12f,%16.12f,%16.12f, Im.[d0-3] = %16.12f,%16.12f,%16.12f,%16.12f\n",j1/p1,ref2[j2],ref2[j2+2],ref2[j2+4],ref2[j2+6],ref2[j2+1],ref2[j2+3],ref2[j2+5],ref2[j2+7]);
			dtmp = fabs(a[j1  ] - ref2[j2  ]); avg_err += dtmp; if(dtmp > 1e-10){ /*printf("error Re.d0\n");*/ nerr++; };
			dtmp = fabs(a[j1+1] - ref2[j2+2]); avg_err += dtmp; if(dtmp > 1e-10){ /*printf("error Re.d1\n");*/ nerr++; };
			dtmp = fabs(a[j1+2] - ref2[j2+4]); avg_err += dtmp; if(dtmp > 1e-10){ /*printf("error Re.d2\n");*/ nerr++; };
			dtmp = fabs(a[j1+3] - ref2[j2+6]); avg_err += dtmp; if(dtmp > 1e-10){ /*printf("error Re.d3\n");*/ nerr++; };
			dtmp = fabs(a[j   ] - ref2[j2+1]); avg_err += dtmp; if(dtmp > 1e-10){ /*printf("error Im.d0\n");*/ nerr++; };
			dtmp = fabs(a[j +1] - ref2[j2+3]); avg_err += dtmp; if(dtmp > 1e-10){ /*printf("error Im.d1\n");*/ nerr++; };
			dtmp = fabs(a[j +2] - ref2[j2+5]); avg_err += dtmp; if(dtmp > 1e-10){ /*printf("error Im.d2\n");*/ nerr++; };
			dtmp = fabs(a[j +3] - ref2[j2+7]); avg_err += dtmp; if(dtmp > 1e-10){ /*printf("error Im.d3\n");*/ nerr++; };
		#elif defined(USE_SSE2)
		//	printf("Out%2u Re.[d0,d1] = %16.12f,%16.12f, Im.[d0,d1] = %16.12f,%16.12f\n",j1/p1,a[j1],a[j1+1],a[j1+2],a[j1+3]);
		//	printf("Ref%2u Re.[d0,d1] = %16.12f,%16.12f, Im.[d0,d1] = %16.12f,%16.12f\n",j1/p1,ref1[j2],ref1[j2+2],ref1[j2+1],ref1[j2+3]);
			dtmp = fabs(a[j1  ] - ref2[j2  ]); avg_err += dtmp; if(dtmp > 1e-10){ /*printf("error Re.d0\n");*/ nerr++; };
			dtmp = fabs(a[j1+1] - ref2[j2+2]); avg_err += dtmp; if(dtmp > 1e-10){ /*printf("error Re.d1\n");*/ nerr++; };
			dtmp = fabs(a[j1+2] - ref2[j2+1]); avg_err += dtmp; if(dtmp > 1e-10){ /*printf("error Im.d0\n");*/ nerr++; };
			dtmp = fabs(a[j1+3] - ref2[j2+3]); avg_err += dtmp; if(dtmp > 1e-10){ /*printf("error Im.d1\n");*/ nerr++; };
		#else
		//	printf("Out%2u Re,Im = %16.12f,%16.12f\n",j1/p1,a[j1],a[j1+1]);
			dtmp = fabs(a[j1  ] - ref2[j2  ]); avg_err += dtmp; if(dtmp > 1e-10){ /*printf("error Re.d0\n");*/ nerr++; };
			dtmp = fabs(a[j1+1] - ref2[j2+1]); avg_err += dtmp; if(dtmp > 1e-10){ /*printf("error Im.d0\n");*/ nerr++; };
		#endif
		}
		ASSERT(HERE, nerr == 0, "DIT Outputs mismatch ref-data!");
		printf("\tSummed roundoff error = %20.10e]\n",avg_err);

	#ifdef USE_SSE2
		free((void *)sc_arr);	sc_arr=0x0;
	#endif
		return nerr;
	}
  #endif	// USE_SSE2?

  #ifndef USE_ARM_V8_SIMD
	// Timing loop for radix-32 DIF macro:
	int	test_radix32_dft()
	{
		const char func[] = "test_radix32_dft";
		/*...time-related stuff	*/
		double clock1, clock2;
		double *dptr, tdiff, rt,it, dtmp, avg_err;
		int i,j,j1,j2,k,imax = 1000000, nerr = 0;	// Use 10^6 loop execs
		int p01,p02,p03,p04,p05,p06,p07,p08,p0C,p10,p14,p18,p1C;
		const double c = 0.92387953251128675613, s     = 0.38268343236508977173	/* exp[  i*(twopi/16)]	*/
				,c32_1 = 0.98078528040323044912, s32_1 = 0.19509032201612826784	/* exp(  i*twopi/32), the radix-32 fundamental sincos datum	*/
				,c32_3 = 0.83146961230254523708, s32_3 = 0.55557023301960222473;/* exp(3*i*twopi/32)	*/
		// DIF[ref1], DIT[ref2] ref-outputs: cols 0,1 are re,im outputs for scalar-mode,
		// cols [0,1],[2,3] are [re0,im0],[re1,im1] for SSE2 mode,
		// cols [0,1],[2,3],[4,5],[6,7] are [re0,im0],[re1,im1],[re2,im2],[re3,im3] for AVX/AVX2 mode:
		const double ref1[] = {	// DIF ref-outputs:
			140.024722521050,140.871739019788, 150.122526863264,150.790780470346, 153.197536655879,141.888816564726, 148.200334262373,134.812656830225,
			 -7.022373589968,  3.997069237833,   1.159133789373,-13.938925439525, -13.901657582104,-13.124482281622,  -2.999617826599,  4.019943693420,
			-32.977206808047,-14.181971377571,  -5.941150622203, -3.033534377769, -11.089824908060, 31.985452203788,  -3.669577041080,-34.995027947290,
			  2.947291144157,  5.983799556561,  26.055011959684,  3.038798268224,  -9.000352360218,  3.903863634243,   1.895737733555,  8.877881373357,
			  7.721695233893, 16.796650048533,  -1.649521343162, -0.277518276115,  -2.071817393690, 13.444818786285,  25.317663367727, 37.495197721454,
			-31.977013450316, 25.009310552180, -20.160830489899, -1.920722367161, -11.991767459065, 20.393896379722, -15.580661792276,  4.634000273354,
			 14.178689335366,  2.677210684856,   5.955945647964,-15.836133052721, -11.915297642057, -3.319147296605, -21.843241613928, -8.366802325879,
			-14.001726748089,-28.790061578753,  -3.984933489930, -6.185668330702, -22.209611107602, 29.068526143826, -12.083175812057, -9.916699134160,
			 -6.568428679973,-25.666742139307, -10.520259251446, -4.129936044069, -12.781454011841, -4.869019419839,   3.482712461121, 27.523676200655,
			 -1.057743014962,-14.704293563743,   0.091428302251,-12.571574727013,  -3.287527844215,-14.232470955807,   6.255645869851, -7.339896665050,
			  0.209684129827,  0.547094194837, -10.963984957490, -3.876793343529,  23.218931543837,-18.214148018491,  11.657368908341, 17.991732248441,
			-28.104387550397, -4.400564730738,  25.563818204521, 12.726467336879,  25.580181542195,-14.638111116564, -13.218849203121,  2.040527144615,
			  2.083180593821, 15.885493780831, -23.328762279212,-11.286300882939,   4.860180033335, 12.239324275083,   5.440404097349, -2.575427656550,
			 23.926133714609,  3.108615496842,   1.776500564269,  0.068947971316,  12.647835095707, 29.400844019242,   2.418436230672,-24.083351693267,
			 -0.347455062825,  6.283821567638,  16.843766876447,  4.680856333787, -39.799917123612,  4.754848839740,  21.684071271268, -4.969665935755,
			-21.918629157623, 26.580118490908,  -7.511066523755, -1.759378744471,  -1.946308210422,  5.315482264882,   2.843991214204,-16.363953283434,
			 14.337138897418, -4.962194106347,  17.400826267082, 18.868370862396,  22.650737301208,  5.471644515705,  -3.481341197924, 27.691226024499,
			 21.041509232485,-31.915501765817, -15.142982973770,-20.430421757509,  11.659109059189,-10.116796178117, -11.836298700689, -3.327236349704,
			 37.176413583993, 19.857258124098, -13.987638594313, 17.838096637263,  -5.124495180610, 23.327444083479,   2.804503746718, -0.443969092161,
			  1.406240377551,-26.736449239665, -40.981166734124,  6.201366173317,  -6.432953497570,-15.977458781115,  18.961280110441, 32.499979534635,
			  7.921206529892,-11.357602373281, -34.277534382182,  6.134348984436,  16.337886041884, 12.627179112788,  13.813910679703,  7.030688270899,
			  4.324539783228,-28.533251396752,   4.360438423123,  2.932955337777, -28.430219636280,  5.281214849604,  17.891080818108,  9.968551200424,
			  7.941424352566, -1.696843470649,  15.757691607411,-24.434672430907,   7.914095843634,  9.431160517530,  -1.598367331510, -4.272930898878,
			  2.511150536063, 14.009852378142,  -4.834059506334, -7.550322452784,  37.865137063644,-13.600877868521, -12.482983358963,-12.875845521982,
			 -3.176289269507,-13.256717474827,  -6.547259678267, -1.377382488463, -21.939148433863, 12.042048878725, -14.253417288097, -3.116163694772,
			-13.943561648678,-20.809854886680,  22.715210033873,-17.918211233280,   9.599530649074, 15.040789321374, -22.179669299814,  3.258448604457,
			 17.723097112906, -2.115349487764, -10.249196792606, 10.854902104522,   2.652452128636,  5.342084280612, -15.832686632635, 21.110019793589,
			 12.429420600041,  5.262100716978, -17.218761890202,-23.617070473007,  -0.887637901968, 25.148631396789, -30.012363448624,-40.310433439768,
			 -4.141708826570, 17.834348787749,  11.506602502169, 13.208475215709,   4.843614510496, 26.588636048305,  11.200296759274,  5.815712940036,
			-20.903728146640,-12.233903905367,  11.822574924896,-26.212273886849,  13.968586324442,-16.985128019648, -25.517711247127,-21.063662783561,
			 -4.668605110589,-23.518389239674,  17.853143730060, 10.255539627711,   7.131503684185,-21.262397894420, -32.978064448385, 17.474810119157,
			-31.094680614681, -7.824791900837,  26.314489812508,-29.243065014873,   8.682672815832,  1.643331714301,   9.700588712124, 23.776014448994
		};
		const double ref2[] = {	// DIT ref-outputs
			141.000000000000,140.000000000000, 151.000000000000,150.000000000000, 154.000000000000,141.000000000000, 149.000000000000,134.000000000000,
			 -3.783906043215, -6.736206535300,   4.973271752545, -5.600710881344, -27.516912975494, 10.704249544145,  16.298563940304,-12.532144509101,
			-13.975320751204, 16.069938560493,   3.827475127678,-11.577279846112,  -6.388477515408, 22.564077591341,  -3.986656256409, 36.982179763860,
			  4.156510948913, 41.173817228073, -17.902557040341, 12.221997670440,  17.264655075949, 22.902285126847, -16.713188704199, 17.223442192287,
			 -4.035494871322,-23.135972475146,  -7.957748291754, -5.209624508116, -13.026778182253, 15.312894054796,  -6.331210248048,-11.404515029260,
			 -1.213343589390,  2.270312000675,   2.392595859890,-24.873739661946,  12.007132532990, 23.996630560387, -30.638644917530, -4.507899445454,
			 -4.948713292598,  9.845221290990,  13.661265407876,  5.066750751913,   5.565104645159,-11.248394412931, -12.521767442562,-14.660708821623,
			-11.796906859905, -1.138597253783, -12.919382326987,-22.841118709164, -10.120884281347, -6.529844383861,  -3.726060004042, -9.189821595846,
			-26.975329278755,  8.082797183077,  22.033643988203, 10.932453183120,   9.073588606597, 23.972275440563,  -3.006121794955, -1.990786717330,
			 -9.066596305134, -1.595969993346,  -1.938409248399, -8.058592277841, -26.733887573698, -0.172909814917,   2.941161411616,  3.986398982278,
			 32.576855326640,  4.614723738419, -12.340924255396,-19.211551338795,  35.403183820966, 15.080943005094, -21.406985527751, -9.134686122438,
			 -4.548284998661,-12.297167806507,  -3.624267744685, -8.678742216190, -34.379691370578,  8.355427602184,   3.089654213922, 23.372398473839,
			 -5.905151007980, -1.215478489053,   0.351694332823, 12.776691413054,  21.280023186453, 18.994155125044,   3.489243487555,-20.986842217176,
			  3.602518389431,-14.276719830409, -26.989840512597,-47.671396585122,  -9.847473285215,-22.401510403948, -10.119044829611, 15.796167140639,
			 34.812232991152,-17.729086047369, -25.869666551597, 12.520629440701,  -3.784245396586, -2.029733818180,   2.402928300177, 16.401641652305,
			 18.052147777352,-33.857249786628,   6.801921794274, -6.576455795400,  13.566832313677, 11.095574986667, -14.340101625722,-24.725755028083,
			  2.766779909180,-38.017692314686,  37.011575254755,  1.772934618546,  16.030378227767,  4.901731722027,   3.049030603041,  7.981441748313,
			 -2.730228499662,-22.877430624432,  27.625582010048,  8.102003437602,   1.277400517221,-15.187342092285, -19.625868405112, -1.776741745095,
			  0.525658067551, -0.891291801172,  24.353656294353,  0.820464065546,  11.633979662105, -6.396990493280,  20.190490658467,-30.391285275285,
			 -9.009911878546, -5.397713458813,  11.696339407365,  3.869141878254,   3.252668870016, 13.349874761303, -30.365441147606,  7.032218838475,
			 -3.960443046134,  5.172663690128,   1.790325078625,-20.792517661947, -22.820931786580, 16.882635998245,  16.247377143023, -8.710657249299,
			 20.424243482247, -7.632556998691, -12.288080889385,  9.778300008816,   5.181110844942,  5.760433633129, -34.631930192857, 24.613585612048,
			  7.664220805692,-16.100052979546, -18.783320226262, 16.475099436075,  33.176479850760, 28.370879660324, -25.509413453942,  5.348541054986,
			 14.931180922425, 17.307008785486,  -5.480746581315,-10.068509993107,  10.669638727806, -7.286566226329,   8.364260807000,  7.732785362298,
			-13.165116965253,-17.879588789769,  -6.174617208171,-18.943972717795,   8.981211291404, -2.082749082142,  39.292868290553, 31.639698189163,
			 -9.231165386395, 33.006526880669,   9.123349632636, 16.893393598794, -39.828580130830,-14.182339621860,  30.828764118530, 22.825796980260,
			  5.057098954615, 16.039181372271, -35.856548456771,  6.202346795540,  -4.691311299012, -7.412834751256,   0.892138240784, -9.520577057717,
			-21.518086979503,-23.863875946161,  -4.011715939162,-18.450590208336,  12.326906368091, -3.662946603667,  -2.667379125681,-19.089339157323,
			 13.976462633395,  7.092975107311,   1.677177274040, -2.796344761165,   2.602447851339,-19.120929617852,   6.553117262884, 12.900941352274,
			-13.860571623143,  9.956334468539,  -2.970997394616, -4.569778700177,  13.162621390589, -4.204190717894, -18.569819812753,  0.007457485465,
			-38.057933128419,-27.821062813343,  -4.719711981316,  6.258941911695, -22.560724338994, 16.896478631264,   7.785441577353,  5.052808814922,
			 -6.672156694112, -5.749720392297,   9.001634924866, -5.548569238130,  -4.361598749612,  9.657701141829,   8.205352378866, -2.180603960068
		};

		const int stride = 2*RE_IM_STRIDE, dim = stride<<5, idx[32] = {0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,62};
		double cc[32],ss[32];
		static double *a,*a_ptr;	// Dimension = number of scalar-doubles in 16 vector-complex in SIMD build mode
		a_ptr = ALLOC_VEC_DBL(a_ptr, dim/RE_IM_STRIDE);	if(!a_ptr){ sprintf(cbuf, "FATAL: unable to allocate a_ptr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		a     = ALIGN_VEC_DBL(a_ptr);
		ASSERT(HERE, ((long)a & SZ_VDM1) == 0, "a0_ptr not 64-byte aligned!");
	#ifdef USE_SSE2
		const int pfetch_dist = 0;
		int pfetch_addr = 0;	// Don't care about pfetch in this lcal-mem context, so just set these = 0
		static vec_dbl *sc_arr = 0x0, *sc_ptr;
		double *add0;	/* Addresses into array sections */
		vec_dbl *c_tmp,*s_tmp;
		static vec_dbl *isrt2,*sqrt2, *cc0, *ss0, *cc1, *ss1, *cc3, *ss3, *one,*two, *r00,*r10,*r20,*r30;
		sc_arr = ALLOC_VEC_DBL(sc_arr, 0x90);	if(!sc_arr){ sprintf(cbuf, "FATAL: unable to allocate sc_arr in %s.\n",func); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		sc_ptr = ALIGN_VEC_DBL(sc_arr);
		ASSERT(HERE, ((long)sc_ptr & SZ_VDM1) == 0, "sc_ptr not 64-byte aligned!");
		r00 = sc_ptr;
		r10 = r00 + 0x10;
		r20 = r00 + 0x20;
		r30 = r00 + 0x30;
		isrt2 = r00 + 0x40;
		cc0	  = r00 + 0x41;
		ss0	  = r00 + 0x42;
		cc1	  = r00 + 0x43;
		ss1	  = r00 + 0x44;
		cc3	  = r00 + 0x45;
		ss3	  = r00 + 0x46;
		one   = r00 + 0x87;
		two   = r00 + 0x88;
		sqrt2 = r00 + 0x89;
		/* These remain fixed: */
		VEC_DBL_INIT(isrt2, ISRT2);		VEC_DBL_INIT(sqrt2, SQRT2);
		VEC_DBL_INIT(one  , 1.0  );		VEC_DBL_INIT(two, 2.0  );
		VEC_DBL_INIT(cc0  , c    );		VEC_DBL_INIT(ss0, s    );
		VEC_DBL_INIT(cc1  , c32_1);		VEC_DBL_INIT(ss1, s32_1);
		VEC_DBL_INIT(cc3  , c32_3);		VEC_DBL_INIT(ss3, s32_3);
	#endif	// USE_SSE2 ?

		// Do these timing-test DFTs in-place, i.e. p01 = #doubles in a pair of vec_dbl:
		p01 = RE_IM_STRIDE << 1;
		p02 = p01 +p01;
		p03 = p02 +p01;
		p04 = p03 +p01;
		p05 = p04 +p01;
		p06 = p05 +p01;
		p07 = p06 +p01;
		p08 = p04 +p04;
		p0C = p08 +p04;
		p10 = p0C +p04;
		p14 = p10 +p04;
		p18 = p14 +p04;
		p1C = p18 +p04;

		// Twiddles for the purpose of this timing-test are w1-15, with w := exp(2*Pi*I/2^14):
		cc[0   ] = 1.0                           ; ss[0   ] = 0.0                             ;
		cc[0x01] = 0.9999999264657178511447314807; ss[0x01] = 0.0003834951875713955890724616812;
		cc[0x02] = 0.9999997058628822191602282177; ss[0x02] = 0.0007669903187427045269385683580;
		cc[0x03] = 0.9999993381915255477888066109; ss[0x03] = 0.001150485337113848457071735047;
		cc[0x04] = 0.9999988234517019099290257101; ss[0x04] = 0.001533980186284765612303697150;
		cc[0x05] = 0.9999981616434870076277347923; ss[0x05] = 0.001917474809855419109500620455;
		cc[0x06] = 0.9999973527669781720689399696; ss[0x06] = 0.002300969151425805244235552264;
		cc[0x07] = 0.9999963968222943635594898320; ss[0x07] = 0.002684463154595961785455992532;
		cc[0x08] = 0.9999952938095761715115801256; ss[0x08] = 0.003067956762965976270145365491;
		cc[0x09] = 0.9999940437289858144220774704; ss[0x09] = 0.003451449920135994297977171937;
		cc[0x0A] = 0.9999926465807071398486621178; ss[0x0A] = 0.003834942569706227825960602960;
		cc[0x0B] = 0.9999911023649456243827897550; ss[0x0B] = 0.004218434655276963463076393843;
		cc[0x0C] = 0.9999894110819283736194723572; ss[0x0C] = 0.004601926120448570764901699143;
		cc[0x0D] = 0.9999875727319041221238780943; ss[0x0D] = 0.004985416908821510528222769585;
		cc[0x0E] = 0.9999855873151432333947502950; ss[0x0E] = 0.005368906963996343085634209014;
		cc[0x0F] = 0.9999834548319376998246454755; ss[0x0F] = 0.005752396229573736600123594041;
		cc[0x10] = 0.9999811752826011426569904375; ss[0x10] = 0.006135884649154475359640234589;
		cc[0x11] = 0.9999787486674688119399584425; ss[0x11] = 0.006519372166339468071646855519;
		cc[0x12] = 0.9999761749868975864771644677; ss[0x12] = 0.006902858724729756157652981865;
		cc[0x13] = 0.9999734542412659737751795536; ss[0x13] = 0.007286344267926522047728803984;
		cc[0x14] = 0.9999705864309741099878642477; ss[0x14] = 0.007669828739531097474998305920;
		cc[0x15] = 0.9999675715564437598575211556; ss[0x15] = 0.008053312083144971770110435865;
		cc[0x16] = 0.9999644096181183166528666053; ss[0x16] = 0.008436794242369800155687098562;
		cc[0x17] = 0.9999611006164628021038214358; ss[0x17] = 0.008820275160807412040746750608;
		cc[0x18] = 0.9999576445519638663331209194; ss[0x18] = 0.009203754782059819315102378107;
		cc[0x19] = 0.9999540414251297877847438260; ss[0x19] = 0.009587233049729224643732638124;
		cc[0x1A] = 0.9999502912364904731491606429; ss[0x1A] = 0.009970709907418029761124940991;
		cc[0x1B] = 0.9999463939865974572854009582; ss[0x1B] = 0.01035418529872884376558925796;
		cc[0x1C] = 0.9999423496760239031399400209; ss[0x1C] = 0.01073765916726449141354143107;
		cc[0x1D] = 0.9999381583053646016624044894; ss[0x1D] = 0.01112113145662802141375476658;
		cc[0x1E] = 0.9999338198752359717180973806; ss[0x1E] = 0.01150460211042271472157869220;
		cc[0x1F] = 0.9999293343862760599973422319; ss[0x1F] = 0.01188807107225209283312325799;

	#ifdef USE_SSE2
		/* Sincos data stored in BRed form in SSE2 local-data layout:
		DIF: (cc0,ss0) + 0x[06,26,16,36|0e,2e,1e,3e|0a,2a,1a,3a|12,32,22,42|08,28,18,38|10,30,20,40|0c,2c,1c,3c|14,34,24,44]
			= "  " + 6 + 0x[00,20,10,30|08,28,18,38|04,24,14,34|0c,2c,1c,3c|02,22,12,32|0a,2a,1a,3a|06,26,16,36|0e,2e,1e,3e]
			= "  " + 6 +   [ 0,16, 8,24| 4,20,12,28| 2,18,10,26| 6,22,14,30| 1,17, 9,25| 5,21,13,29| 3,19,11,27| 7,23,15,31], straight bit-reversal.
		These are the addess-offsets of roots 0-31 ... flipping this around and asking what is the roots-order in a linear walk thru memory
		(i.e. the index-offsets or roots 0,1,2,3,... in the last [] sequence above) we get the very same thing, i.e. the index-permutation
		is its own inverse; this is a property of bit-reversal reordering.
		*/
		c_tmp = cc0 + 0x06; s_tmp = c_tmp+1;	/* c0,s0 */
		for(i = 0; i < 32; i++, c_tmp+=2, s_tmp+=2) {
			j = reverse(i,32);
			VEC_DBL_INIT(c_tmp, cc[j]);	VEC_DBL_INIT(s_tmp, ss[j]);
		}
	#endif

		//******************* Timing loop for Radix-32 DIF transform macro: *******************
		clock1 = getRealTime();
		for(i = 0; i < imax; i++) {
			// Copy digits of Pi-data into our vec_dbl inputs:
			for(j1 = 0, j2 = 0; j1 < dim; j1 += stride, j2 += 8)	// j2 is base-index into ran[] input array
			{
			#ifdef USE_AVX512	// Set this up so that AVX-512 can use the same ref-data as AVX:
				a[j1   ] = ran[j2  ];	/* Re0 */	a[j1+ 4] = 0;	/* Re4 */
				a[j1+ 1] = ran[j2+2];	/* Re1 */	a[j1+ 5] = 0;	/* Re5 */
				a[j1+ 2] = ran[j2+4];	/* Re2 */	a[j1+ 6] = 0;	/* Re6 */
				a[j1+ 3] = ran[j2+6];	/* Re3 */	a[j1+ 7] = 0;	/* Re7 */
				a[j1+ 8] = ran[j2+1];	/* Im0 */	a[j1+12] = 0;	/* Im4 */
				a[j1+ 9] = ran[j2+3];	/* Im1 */	a[j1+13] = 0;	/* Im5 */
				a[j1+10] = ran[j2+5];	/* Im2 */	a[j1+14] = 0;	/* Im6 */
				a[j1+11] = ran[j2+7];	/* Im3 */	a[j1+15] = 0;	/* Im7 */
			#elif defined(USE_AVX)
				a[j1  ] = ran[j2  ];	// Re0
				a[j1+1] = ran[j2+2];	// Re1
				a[j1+2] = ran[j2+4];	// Re2
				a[j1+3] = ran[j2+6];	// Re3
				a[j1+4] = ran[j2+1];	// Im0
				a[j1+5] = ran[j2+3];	// Im1
				a[j1+6] = ran[j2+5];	// Im2
				a[j1+7] = ran[j2+7];	// Im3
			#elif defined(USE_SSE2)
				a[j1  ] = ran[j2  ];	// Re0
				a[j1+1] = ran[j2+2];	// Re1
				a[j1+2] = ran[j2+1];	// Im0
				a[j1+3] = ran[j2+3];	// Im1
				a[j1+4] = ran[j2+4];	// Re2
				a[j1+5] = ran[j2+6];	// Re3
				a[j1+6] = ran[j2+5];	// Im2
				a[j1+7] = ran[j2+7];	// Im3
			#else
				a[j1  ] = ran[j2  ];	// Re0
				a[j1+1] = ran[j2+1];	// Im0
				a[j1+2] = ran[j2+2];	// Re1
				a[j1+3] = ran[j2+3];	// Im1
				a[j1+4] = ran[j2+4];	// Re2
				a[j1+5] = ran[j2+5];	// Im2
				a[j1+6] = ran[j2+6];	// Re3
				a[j1+7] = ran[j2+7];	// Im3
			#endif
			}
			j1 = 0; j2 = RE_IM_STRIDE;
		#ifdef USE_SSE2
			SSE2_RADIX32_DIF_TWIDDLE(a,p01,p02,p03,p04,p08,p0C,p10,p18,r00)
		#else
			RADIX_32_DIF_TWIDDLE_OOP( a,idx, a,idx,	// This DFT is in-place
									 cc[0x10],ss[0x10], cc[0x08],ss[0x08], cc[0x18],ss[0x18]
				, cc[0x04],ss[0x04], cc[0x14],ss[0x14], cc[0x0C],ss[0x0C], cc[0x1C],ss[0x1C]
				, cc[0x02],ss[0x02], cc[0x12],ss[0x12], cc[0x0A],ss[0x0A], cc[0x1A],ss[0x1A]
				, cc[0x06],ss[0x06], cc[0x16],ss[0x16], cc[0x0E],ss[0x0E], cc[0x1E],ss[0x1E]
				, cc[0x01],ss[0x01], cc[0x11],ss[0x11], cc[0x09],ss[0x09], cc[0x19],ss[0x19]
				, cc[0x05],ss[0x05], cc[0x15],ss[0x15], cc[0x0D],ss[0x0D], cc[0x1D],ss[0x1D]
				, cc[0x03],ss[0x03], cc[0x13],ss[0x13], cc[0x0B],ss[0x0B], cc[0x1B],ss[0x1B]
				, cc[0x07],ss[0x07], cc[0x17],ss[0x17], cc[0x0F],ss[0x0F], cc[0x1F],ss[0x1F]
			);
		#endif
		}
		clock2 = getRealTime();
		tdiff = (double)(clock2 - clock1);
		printf("%s: Time for %u DIF macro calls =%s [tdiff = %20.10e]\n",func, imax, get_time_str(tdiff), tdiff);

		// Check outputs vs ref-data:
		nerr = 0;	dtmp = avg_err = 0.0;
		for(j1 = 0, j2 = 0; j1 < dim; j1 += stride, j2 += 8)	// j2 is base-index into ref-array
		{
			j = j1+RE_IM_STRIDE;
		#ifdef USE_AVX	// Since we set up AVX-512 mode to only use nonzero data in lower 4 double-slots of each 8-vector, can use same code here:
		//	printf("Out[%2u] Re.[d0-3] = %16.12f,%16.12f,%16.12f,%16.12f, Im.[d0-3] = %16.12f,%16.12f,%16.12f,%16.12f\n",j1/p01,a[j1],a[j1+1],a[j1+2],a[j1+3],a[j],a[j+1],a[j+2],a[j+3]);
		//	printf("Ref[%2u] Re.[d0-3] = %16.12f,%16.12f,%16.12f,%16.12f, Im.[d0-3] = %16.12f,%16.12f,%16.12f,%16.12f\n",j1/p01,ref1[j2],ref1[j2+2],ref1[j2+4],ref1[j2+6],ref1[j2+1],ref1[j2+3],ref1[j2+5],ref1[j2+7]);
			dtmp = fabs(a[j1  ] - ref1[j2  ]); avg_err += dtmp; if(dtmp > 1e-10){ /*printf("error Re.d0\n");*/ nerr++; };
			dtmp = fabs(a[j1+1] - ref1[j2+2]); avg_err += dtmp; if(dtmp > 1e-10){ /*printf("error Re.d1\n");*/ nerr++; };
			dtmp = fabs(a[j1+2] - ref1[j2+4]); avg_err += dtmp; if(dtmp > 1e-10){ /*printf("error Re.d2\n");*/ nerr++; };
			dtmp = fabs(a[j1+3] - ref1[j2+6]); avg_err += dtmp; if(dtmp > 1e-10){ /*printf("error Re.d3\n");*/ nerr++; };
			dtmp = fabs(a[j   ] - ref1[j2+1]); avg_err += dtmp; if(dtmp > 1e-10){ /*printf("error Im.d0\n");*/ nerr++; };
			dtmp = fabs(a[j +1] - ref1[j2+3]); avg_err += dtmp; if(dtmp > 1e-10){ /*printf("error Im.d1\n");*/ nerr++; };
			dtmp = fabs(a[j +2] - ref1[j2+5]); avg_err += dtmp; if(dtmp > 1e-10){ /*printf("error Im.d2\n");*/ nerr++; };
			dtmp = fabs(a[j +3] - ref1[j2+7]); avg_err += dtmp; if(dtmp > 1e-10){ /*printf("error Im.d3\n");*/ nerr++; };
		#elif defined(USE_SSE2)
		//	printf("Out%2u Re.[d0,d1] = %16.12f,%16.12f, Im.[d0,d1] = %16.12f,%16.12f\n",j1/p01,a[j1],a[j1+1],a[j1+2],a[j1+3]);
		//	printf("Ref%2u Re.[d0,d1] = %16.12f,%16.12f, Im.[d0,d1] = %16.12f,%16.12f\n",j1/p01,ref1[j2],ref1[j2+2],ref1[j2+1],ref1[j2+3]);
			nerr += (fabs(a[j1  ] - ref1[j2  ]) > 1e-10);
			nerr += (fabs(a[j1+1] - ref1[j2+2]) > 1e-10);
			nerr += (fabs(a[j1+2] - ref1[j2+1]) > 1e-10);
			nerr += (fabs(a[j1+3] - ref1[j2+3]) > 1e-10);
		#else
		//	printf("Out%2u Re,Im = %16.12f,%16.12f\n",j1/p01,a[j1],a[j1+1]);
			nerr += (fabs(a[j1  ] - ref1[j2  ]) > 1e-10);
			nerr += (fabs(a[j1+1] - ref1[j2+1]) > 1e-10);
		#endif
		}
		ASSERT(HERE, nerr == 0, "DIF Outputs mismatch ref-data!");
		printf("\tSummed roundoff error = %20.10e]\n",avg_err);
	#if 0
		10^6-timing:	setup	+=DIF	DIF-only
		sse2:			.386	.724	.338 [676 cycles]
		avx2:			.058	.530	.472 [614 cycles]
		avx512:			.100	.611	.511 [664 cycles]	Optimized address-comp and using zmm28-31 for consts ==> 585 cycles; still crappy vs
															tangent+FMA - optimized radix-16 @214 cycles, implies target of 535 cycles for radix-32.
	#endif

		//******************* Timing loop for Radix-32 DIT transform macro: *******************
		/* Sincos data in SSE2 local-data layout are a bit funky for DIT:
		DIT: (cc0,ss0) + 0x[06,26,16,36|0e,2e,1e,3e|08,28,18,38|10,30,20,40|0a,2a,1a,3a|12,32,22,42|0c,2c,1c,3c|14,34,24,44].
			= "  " + 6 + 0x[00,20,10,30|08,28,18,38|02,22,12,32|0a,2a,1a,3a|04,24,14,34|0c,2c,1c,3c|06,26,16,36|0e,2e,1e,3e]
			= "  " + 6 +   [ 0,16, 8,24| 4,20,12,28| 1,17, 9,25| 5,21,13,29| 2,18,10,26| 6,22,14,30| 3,19,11,27| 7,23,15,31], swap quartets 2,3 <-> 4,5 relative to BR!
							 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31
			These are the addess-offsets of roots 0-31 ... flipping this around and asking what is the roots-order in a linear walk thru memory
			(i.e. the index-offsets or roots 0,1,2,3,... in the last [] sequence above) we get
							0,16,8,24, 4,12,20,28, 2,10,18,26, 6,14,22,30, 1,9,17,25, 5,13,21,29, 3,11,19,27, 7,15,23,31,
			i.e. the index-permutation *not* is its own inverse, due to our nor longer having a strict bit-reversal reordering. It is this latter perm that apears in dit_funky_br[] below:
		*/
	#ifdef USE_SSE2
		// This 'funky' perm is just the BR-perm with the middle 2 terms of each quartet swapped, i.e. quartet elements appearing in-order:
		const int dit_funky_br[32] = {0,8,16,24, 4,12,20,28, 2,10,18,26, 6,14,22,30, 1,9,17,25, 5,13,21,29, 3,11,19,27, 7,15,23,31};
		c_tmp = cc0 + 0x06; s_tmp = c_tmp+1;	/* c0,s0 */
		for(i = 0; i < 32; i++, c_tmp+=2, s_tmp+=2) {
			j = dit_funky_br[i];
			VEC_DBL_INIT(c_tmp, cc[j]);	VEC_DBL_INIT(s_tmp, ss[j]);
		}
	#endif
		clock1 = getRealTime();
		for(i = 0; i < imax; i++) {
			// Copy digits of Pi-data into our vec_dbl inputs:
			for(j1 = 0, j2 = 0; j1 < dim; j1 += stride, j2 += 8)	// j2 is base-index into ran[] input array
			{
			#ifdef USE_AVX512	// Set this up so that AVX-512 can use the same ref-data as AVX:
				a[j1   ] = ran[j2  ];	/* Re0 */	a[j1+ 4] = 0;	/* Re4 */
				a[j1+ 1] = ran[j2+2];	/* Re1 */	a[j1+ 5] = 0;	/* Re5 */
				a[j1+ 2] = ran[j2+4];	/* Re2 */	a[j1+ 6] = 0;	/* Re6 */
				a[j1+ 3] = ran[j2+6];	/* Re3 */	a[j1+ 7] = 0;	/* Re7 */
				a[j1+ 8] = ran[j2+1];	/* Im0 */	a[j1+12] = 0;	/* Im4 */
				a[j1+ 9] = ran[j2+3];	/* Im1 */	a[j1+13] = 0;	/* Im5 */
				a[j1+10] = ran[j2+5];	/* Im2 */	a[j1+14] = 0;	/* Im6 */
				a[j1+11] = ran[j2+7];	/* Im3 */	a[j1+15] = 0;	/* Im7 */
			#elif defined(USE_AVX)
				a[j1  ] = ran[j2  ];	// Re0
				a[j1+1] = ran[j2+2];	// Re1
				a[j1+2] = ran[j2+4];	// Re2
				a[j1+3] = ran[j2+6];	// Re3
				a[j1+4] = ran[j2+1];	// Im0
				a[j1+5] = ran[j2+3];	// Im1
				a[j1+6] = ran[j2+5];	// Im2
				a[j1+7] = ran[j2+7];	// Im3
			#elif defined(USE_SSE2)
				a[j1  ] = ran[j2  ];	// Re0
				a[j1+1] = ran[j2+2];	// Re1
				a[j1+2] = ran[j2+1];	// Im0
				a[j1+3] = ran[j2+3];	// Im1
				a[j1+4] = ran[j2+4];	// Re2
				a[j1+5] = ran[j2+6];	// Re3
				a[j1+6] = ran[j2+5];	// Im2
				a[j1+7] = ran[j2+7];	// Im3
			#else
				a[j1  ] = ran[j2  ];	// Re0
				a[j1+1] = ran[j2+1];	// Im0
				a[j1+2] = ran[j2+2];	// Re1
				a[j1+3] = ran[j2+3];	// Im1
				a[j1+4] = ran[j2+4];	// Re2
				a[j1+5] = ran[j2+5];	// Im2
				a[j1+6] = ran[j2+6];	// Re3
				a[j1+7] = ran[j2+7];	// Im3
			#endif
			}
			j1 = 0; j2 = RE_IM_STRIDE;
		#ifdef USE_SSE2
			SSE2_RADIX32_DIT_TWIDDLE(a,p01,p02,p03,p04,p05,p06,p07,p08,p10,p18,r00,r10,r20,r30,isrt2)
		#else
			RADIX_32_DIT_TWIDDLE( a,idx, a,idx,	// This DFT is in-place
									 cc[0x10],ss[0x10], cc[0x08],ss[0x08], cc[0x18],ss[0x18]
				, cc[0x04],ss[0x04], cc[0x14],ss[0x14], cc[0x0C],ss[0x0C], cc[0x1C],ss[0x1C]
				, cc[0x02],ss[0x02], cc[0x12],ss[0x12], cc[0x0A],ss[0x0A], cc[0x1A],ss[0x1A]
				, cc[0x06],ss[0x06], cc[0x16],ss[0x16], cc[0x0E],ss[0x0E], cc[0x1E],ss[0x1E]
				, cc[0x01],ss[0x01], cc[0x11],ss[0x11], cc[0x09],ss[0x09], cc[0x19],ss[0x19]
				, cc[0x05],ss[0x05], cc[0x15],ss[0x15], cc[0x0D],ss[0x0D], cc[0x1D],ss[0x1D]
				, cc[0x03],ss[0x03], cc[0x13],ss[0x13], cc[0x0B],ss[0x0B], cc[0x1B],ss[0x1B]
				, cc[0x07],ss[0x07], cc[0x17],ss[0x17], cc[0x0F],ss[0x0F], cc[0x1F],ss[0x1F]
			);
		#endif
		}
		clock2 = getRealTime();
		tdiff = (double)(clock2 - clock1);
		printf("%s: Time for %u DIT macro calls =%s [tdiff = %20.10e]\n",func, imax, get_time_str(tdiff), tdiff);
		// Check outputs vs ref-data:
		nerr = 0;	dtmp = avg_err = 0.0;
		for(j1 = 0, j2 = 0; j1 < dim; j1 += stride, j2 += 8)	// j2 is base-index into ref-array
		{
			j = j1+RE_IM_STRIDE;
		#ifdef USE_AVX	// Since we set up AVX-512 mode to only use nonzero data in lower 4 double-slots of each 8-vector, can use same code here:
		//	printf("Out[%2u] Re.[d0-3] = %16.12f,%16.12f,%16.12f,%16.12f, Im.[d0-3] = %16.12f,%16.12f,%16.12f,%16.12f\n",j1/p01,a[j1],a[j1+1],a[j1+2],a[j1+3],a[j],a[j+1],a[j+2],a[j+3]);
		//	printf("Ref[%2u] Re.[d0-3] = %16.12f,%16.12f,%16.12f,%16.12f, Im.[d0-3] = %16.12f,%16.12f,%16.12f,%16.12f\n",j1/p01,ref2[j2],ref2[j2+2],ref2[j2+4],ref2[j2+6],ref2[j2+1],ref2[j2+3],ref2[j2+5],ref2[j2+7]);
			dtmp = fabs(a[j1  ] - ref2[j2  ]); avg_err += dtmp; if(dtmp > 1e-10){ /*printf("error Re.d0\n");*/ nerr++; };
			dtmp = fabs(a[j1+1] - ref2[j2+2]); avg_err += dtmp; if(dtmp > 1e-10){ /*printf("error Re.d1\n");*/ nerr++; };
			dtmp = fabs(a[j1+2] - ref2[j2+4]); avg_err += dtmp; if(dtmp > 1e-10){ /*printf("error Re.d2\n");*/ nerr++; };
			dtmp = fabs(a[j1+3] - ref2[j2+6]); avg_err += dtmp; if(dtmp > 1e-10){ /*printf("error Re.d3\n");*/ nerr++; };
			dtmp = fabs(a[j   ] - ref2[j2+1]); avg_err += dtmp; if(dtmp > 1e-10){ /*printf("error Im.d0\n");*/ nerr++; };
			dtmp = fabs(a[j +1] - ref2[j2+3]); avg_err += dtmp; if(dtmp > 1e-10){ /*printf("error Im.d1\n");*/ nerr++; };
			dtmp = fabs(a[j +2] - ref2[j2+5]); avg_err += dtmp; if(dtmp > 1e-10){ /*printf("error Im.d2\n");*/ nerr++; };
			dtmp = fabs(a[j +3] - ref2[j2+7]); avg_err += dtmp; if(dtmp > 1e-10){ /*printf("error Im.d3\n");*/ nerr++; };
		#elif defined(USE_SSE2)
		//	printf("Out%2u Re.[d0,d1] = %16.12f,%16.12f, Im.[d0,d1] = %16.12f,%16.12f\n",j1/p01,a[j1],a[j1+1],a[j1+2],a[j1+3]);
		//	printf("Ref%2u Re.[d0,d1] = %16.12f,%16.12f, Im.[d0,d1] = %16.12f,%16.12f\n",j1/p01,ref1[j2],ref1[j2+2],ref1[j2+1],ref1[j2+3]);
			nerr += (fabs(a[j1  ] - ref2[j2  ]) > 1e-10);
			nerr += (fabs(a[j1+1] - ref2[j2+2]) > 1e-10);
			nerr += (fabs(a[j1+2] - ref2[j2+1]) > 1e-10);
			nerr += (fabs(a[j1+3] - ref2[j2+3]) > 1e-10);
		#else
		//	printf("Out%2u Re,Im = %16.12f,%16.12f\n",j1/p01,a[j1],a[j1+1]);
			nerr += (fabs(a[j1  ] - ref2[j2  ]) > 1e-10);
			nerr += (fabs(a[j1+1] - ref2[j2+1]) > 1e-10);
		#endif
		}
		ASSERT(HERE, nerr == 0, "DIT Outputs mismatch ref-data!");
		printf("\tSummed roundoff error = %20.10e]\n",avg_err);
	#if 0
		10^6-timing:	setup	+=DIF	DIF-only
		sse2:			.386	.736	.350 [700 cycles]
		avx2:			.052	.536	.484 [629 cycles]
		avx512:			.106	.630	.524 [681 cycles]	Optimized address-comp and using zmm28-31 for consts ==> 540 cycles; can we get under 500?
	#endif

	#ifdef USE_SSE2
		free((void *)sc_arr);	sc_arr=0x0;
	#endif
		return nerr;
	}
  #endif	// ifndef(USE_ARM_V8_SIMD?)

#endif	// TEST_SIMD ?

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
			for (j=0; j<100000000; j++) {	/* delay loop */
				k += j*j;
			}
		//	printf("Thread '%d': i = %d, accum = %d\n", me, i, k);
		}
		*(thread_arg->retval) = k;
		pthread_exit(NULL);
	}

	/********* Thread-affinity utilities: *********/

	// Parse a single core-affinity-triplet substring and set the corresponding bits in the global CORE_SET bitmap.
	// Returns: #cores specified in the substring.
	uint32 parseAffinityTriplet(char*istr)
	{
		int ncpu = 0, lo = -1,hi = lo,incr = 1, i,bit,word;
		char *char_addr = istr, *endp;
		ASSERT(HERE, char_addr != 0x0, "Null input-string pointer!");
		size_t len = strlen(istr);
		if(len == 0) return 0;	// Allow 0-length input, resulting in no-op
		ASSERT(HERE, len <= STR_MAX_LEN, "Excessive input-substring length!");
		lo = strtoul(char_addr, &endp, 10);	ASSERT(HERE, lo >= 0, "lo-substring not a valid nonnegative number!");
		if(*endp) {
			ASSERT(HERE, *endp == ':', "Non-colon separator in core-affinity-triplet substring!");
			char_addr = endp+1;
			hi = strtoul(char_addr, &endp, 10);
			ASSERT(HERE, hi >= lo, "hi-substring not a valid number >= lo!");
			if(*endp) {
				ASSERT(HERE, *endp == ':', "Non-colon separator in core-affinity-triplet substring!");
				char_addr = endp+1;
				incr = strtoul(char_addr, &endp, 10);
				ASSERT(HERE, incr > 0, "incr-substring not a valid positive number!");
				ASSERT(HERE, *endp == 0x0, "Non-numeric increment substring in core-affinity-triplet substring!");
			} else {
				// If increment (third) argument of triplet omitted, default to incr = 1.
			}
		} else {
			hi = lo;	// If only 'lo' arg of triplet supplied, take hi = lo and incr=1, i.e. add just CPUid = lo to bitmap
		}
		// CPU set encoded by integer-triplet argument corresponds to values of integer loop
		// index i in the C-loop for(i = lo; i < hi; i += incr), excluding loop-exit value of i:
		for(i = lo; i <= hi; i += incr, ncpu++) {
			word = i>>6; bit = i & 63;	ASSERT(HERE, word < MAX_CORES, "Bitmap word exceeds MAX_CORES!");
			if(CORE_SET[word] & (1ull<<bit)) { sprintf(cbuf, "Core %d multiply specified in affinity-setting!",i);	ASSERT(HERE, 0, cbuf); }
			else { CORE_SET[word] |= 1ull<<bit; }
		}
		return ncpu;
	}

	/******************/
	// Parse a single core-affinity-triplet substring and set the corresponding bits in the global CORE_SET bitmap:
	void parseAffinityString(char*istr)
	{
		uint32 ncpu = 0, i,bit,word,nc, core_count_oflow = 0;
		char *char_addr = istr, *cptr;
		ASSERT(HERE, char_addr != 0x0, "Null input-string pointer!");
		size_t len = strlen(istr);	// length, not counting the \0 string terminator
		ASSERT(HERE, len > 0, "Zero input-string length!");
		ASSERT(HERE, len <= STR_MAX_LEN, "Excessive input-string length!");
		// Clear existing core-affinity bitmap:
		for(i = 0; i < MAX_CORES>>6; i++) { CORE_SET[i] = 0ull; }
		// Affinity-triplet substrings are delimited by commas:
		while(0x0 != (cptr = strchr(char_addr,','))) {
			strncpy(cbuf,char_addr,(cptr-char_addr));	cbuf[cptr-char_addr] = '\0';	// Copy substring into cbuf and null-terminate
			ncpu += parseAffinityTriplet(cbuf);
			char_addr = cptr+1;
		}
		ncpu += parseAffinityTriplet(char_addr);	// Final (or only) core-affinity-triplet
		printf("Set affinity for the following %u cores: ",ncpu);
		nc = 0;
		for(i = 0; i < MAX_CORES; i++) {
			word = i>>6; bit = i & 63;
			if(CORE_SET[word] & (1ull<<bit)) {
				++nc;	printf("%u.",i);
				core_count_oflow += (i >= MAX_THREADS);	// Accumulation (rather than simple if() ... = TRUE) allows us to capture #offenders in the numerical value
			}
		}
		printf("\n");
		ASSERT(HERE, nc == ncpu, "Bitmap #set-bits mismatches #cpu!");
		NTHREADS = ncpu;
		if(NTHREADS > MAX_THREADS) {	// Test this first, since if true, it implies truth of the 'else' conditional
		//	fprintf(stderr,"WARN: NTHREADS = %d exceeds number of logical cores = %d ... Affinities for core indices > %d will be set (mod %d).\n",NTHREADS,MAX_THREADS,MAX_THREADS,MAX_THREADS);
			fprintf(stderr,"ERROR: NTHREADS [ = %d] must not exceed those of available logical cores = 0-%d!\n",NTHREADS,MAX_THREADS-1);
			exit(EXIT_FAILURE);
		} else if(core_count_oflow) {	// This can be true even if #threads within bounds, e.g. 2 available cores (with indices 0,1) and user specifies -cpu 0,2
		//	fprintf(stderr,"WARN: %d cores in user-specified core set has index which exceeds number of logical cores = %d ... Affinities for core indices > %d will be set (mod %d).\n",core_count_oflow,MAX_THREADS,MAX_THREADS,MAX_THREADS);
			fprintf(stderr,"ERROR: %d cores in user-specified core set have index exceeding those of available logical cores = 0-%d!\n",core_count_oflow,MAX_THREADS-1);
			exit(EXIT_FAILURE);
		}
	}

#endif	// MULTITHREAD ?

/***********************/

double get_time(double tdiff)
{
#ifndef MULTITHREAD	// In || mode the mod_square routines use getRealTime() to accumulate wall-clock time, thus CLOCKS_PER_SEC not needed
	return tdiff/CLOCKS_PER_SEC;	/* NB: CLOCKS_PER_SEC may be a phony value used to scale clock() ranges */
#else
	return tdiff;
#endif
}

char*get_time_str(double tdiff)
{
	static char cbuf[STR_MAX_LEN];
#ifndef MULTITHREAD	// In || mode the mod_square routines use getRealTime() to accumulate wall-clock time, thus CLOCKS_PER_SEC not needed
	tdiff /= CLOCKS_PER_SEC;	/* NB: CLOCKS_PER_SEC may be a phony value used to scale clock() ranges */
#endif
	sprintf(cbuf, "%2d%1d:%1d%1d:%1d%1d.%1d%1d%1d"
	,(int)tdiff/36000,((int)tdiff%36000)/3600
	,((int)tdiff%3600)/600,((int)tdiff%600)/60
	,((int)tdiff%60)/10,(int)tdiff%10
	,(int)(10*(tdiff-(int)tdiff)),(int)(100*(tdiff-(int)tdiff))%10,(int)(1000*(tdiff-(int)tdiff))%10);
	return cbuf;
}

// EWM: Jun 2015: This code (and related smaller mods elsewhere) due to Alex Vong,
// as part of his Debian-freeware-packaging-of-Mlucas project:

/* MLUCAS_PATH is the prefix of all files opened by mlucas_fopen()
   It must end with a slash, except when it is an empty sting
   For example, "$HOME/.mlucas.d/" is a valid MLUCAS_PATH

   MLUCAS_PATH is the empty string by default,
   user can set its default value by defining cpp macro MLUCAS_DEFAULT_PATH  */
#ifdef MLUCAS_DEFAULT_PATH
char *MLUCAS_PATH = MLUCAS_DEFAULT_PATH;
#else
char *MLUCAS_PATH = "";
#endif

/* Set the global variable MLUCAS_PATH according to 1. the environment variable
   MLUCAS_PATH and 2. the default value of the global variable MLUCAS_PATH
   Notice 1 has precedence over 2 since only 1 can be set at run-time

   Both 1 and 2 will be expanded by the shell so that user can set the default
   value of the global variable MLUCAS_PATH to be something
   like "$HOME/.mlucas.d/" (See cpp macro MLUCAS_DEFAULT_PATH for more details)
   which can be expanded at run-time

   On sucess, set_mlucas_path() returns silently
   On error, set_mlucas_path() prints the cause of error to stderr
   and calls ASSERT(HERE, 0, "Exiting.");

   possible errors:
   unable to allocate buffer
   unable to open pipe
   path is longer than STR_MAX_LEN
   path does not end with a slash  */
void set_mlucas_path(void)
{
	char *mlucas_path;
	char *cmdstr;
	char *expanded_str;
	int  tmp;
	FILE *pipe_ptr;
	size_t bufsize;
	int has_err = FALSE;

	mlucas_path = getenv("MLUCAS_PATH");
	if (mlucas_path != NULL) {
		bufsize = strlen(mlucas_path) + 1;
		MLUCAS_PATH = (char*)malloc(bufsize); /* will not free!  */
		if (MLUCAS_PATH == NULL) {
			fprintf(stderr, "FATAL: unable to allocate buffer MLUCAS_PATH in set_mlucas_path()\n");
			has_err = TRUE;
			goto out_err_check;
		}
		strcpy(MLUCAS_PATH, mlucas_path);
	} else {
		bufsize = strlen(MLUCAS_PATH) + 1;
	}
	bufsize = (bufsize - 1) * 3 + 1;
	mlucas_path = (char*)malloc(bufsize);
	if (mlucas_path == NULL) {
		fprintf(stderr, "FATAL: unable to allocate buffer mlucas_path in set_mlucas_path()\n");
		has_err = TRUE;
		goto out_err_check;
	}

	quote_spaces(mlucas_path, MLUCAS_PATH);
	cmdstr = (char*)malloc(bufsize + strlen("printf \"\""));
	if (cmdstr == NULL) {
		fprintf(stderr, "FATAL: unable to allocate buffer cmdstr in set_mlucas_path()\n");
		has_err = TRUE;
		goto out_mlucas_path;
	}

	strcpy(cmdstr, "printf \"\"");
	strcat(cmdstr, mlucas_path);
	pipe_ptr = popen(cmdstr, "r");
	if (pipe_ptr == NULL) {
		fprintf(stderr, "FATAL: unable to open pipe pipe_ptr in set_mlucas_path()\n");
		has_err = TRUE;
		goto out_cmdstr;
	}

	tmp = getc(pipe_ptr); /* goto out_pipe if shell output nothing  */
	if (tmp == EOF)
		goto out_pipe;
	else
		ungetc(tmp, pipe_ptr);

	expanded_str = (char*)malloc(STR_MAX_LEN + 1); /* do not free!  */
	if (expanded_str == NULL) {
		fprintf(stderr, "FATAL: unable to allocate buffer expanded_str in set_mlucas_path()\n");
		has_err = TRUE;
		goto out_pipe;
	}
	fgets(expanded_str, STR_MAX_LEN + 1, pipe_ptr);
	if (getc(pipe_ptr) != EOF) {
		fprintf(stderr, "FATAL: environment variable MLUCAS_PATH or cpp macro MLUCAS_DEFAULT_PATH is longer than STR_MAX_LEN in set_mlucas_path()\n");
		has_err = TRUE;
		goto out_pipe;
	}
	if (expanded_str[strlen(expanded_str) - 1] != '/') { /* strlen != 0  */
		fprintf(stderr, "FATAL: environment variable MLUCAS_PATH or cpp macro MLUCAS_DEFAULT_PATH does not end with a slash in set_mlucas_path()\n");
		has_err = TRUE;
		goto out_pipe;
	}

	MLUCAS_PATH = expanded_str;
	out_pipe:
	pclose(pipe_ptr);
	out_cmdstr:
	free(cmdstr);
	out_mlucas_path:
	free(mlucas_path);
	out_err_check:
	if (has_err)
		ASSERT(HERE, 0, "Exiting.");
}

/* Double-quote all spaces in the string pointed by src and write it to dest.
   Suppose src needs bufsize b, dest needs at most bufsize (b - 1) * 3 + 1,
   since at worst the whole string would be made of spaces

   example: `I am happy' will be transformed to `I" "am" "Happy'  */
char *quote_spaces(char *dest, char *src)
{
	size_t i;
	size_t j;

	for (i = 0, j = 0; src[i] != '\0'; ++i, ++j) {
		if (src[i] == ' ') {
			dest[j] = '"';
			++j;
			dest[j] = ' ';
			++j;
			dest[j] = '"';
		} else {
			dest[j] = src[i];
		}
	}
	dest[j] = '\0';
	return dest;
}

/* Emulate `mkdir -p path'
   The command either makes directory `path' and all its parent directories
   or does absolutely nothing

   Return 0 if the directory `path' exists and is writable
   Return 1 if the directory does not exist or is not writable  */
int mkdir_p(char *path)
{
	char mlucas_path[STR_MAX_LEN + 1];
	char cmdstr[4 * STR_MAX_LEN + 1];
	char tmp[4 * STR_MAX_LEN + 1] = "";
	char *tok;
	FILE *fp;

	strcpy(mlucas_path, path);
	if (mlucas_path[0] == '\0')
		return 1;
	else if (mlucas_path[0] == '/')
		strcpy(tmp, "/");

	for (tok = strtok(mlucas_path, "/");
	     tok != NULL;
	     tok = strtok(NULL, "/")) {
		shell_quote(cmdstr, tok);
		strcat(tmp, cmdstr);
		strcat(tmp, "/");
		strcpy(cmdstr, "mkdir ");
		strcat(cmdstr, tmp);
		strcat(cmdstr, " 2> /dev/null");
		system(cmdstr);
	}

	strcat(tmp, "_Mlucas_util_c_mkdir_p_tmp");
	strcpy(cmdstr, "printf ");
	strcat(cmdstr, tmp);
	fp = popen(cmdstr, "r");
	if (fp == NULL) {
		fprintf(stderr, "FATAL: unable to open pipe fp in mkdir_p()\n");
		ASSERT(HERE, 0, "Exiting.");
	}
	fgets(tmp, STR_MAX_LEN + 1, fp);
	pclose(fp);

	fp = fopen(tmp, "a");
	if (fp == NULL)
		return 1;
	fclose(fp);

	strcpy(cmdstr, "rm -f ");
	strcat(cmdstr, tmp);
	strcat(cmdstr, " 2> /dev/null");
	system(cmdstr);
	return 0;
}

/* Double-quote all single-quotes in the string
   and single-quote all other characters
   Suppose src needs bufsize b, dest needs bufsize (b - 1) * 3 + 1,
   since all characters needs quoting

   example: 'a' 'b' will be transformed to "'"'a'"'"' '"'"'b'"'"  */
char *shell_quote(char *dest, char *src)
{
	size_t i;
	size_t j;

	for (i = 0, j = 0; src[i] != '\0'; ++i, ++j) {
		if (src[i] == '\'') {
			dest[j] = '"';
			++j;
			dest[j] = '\'';
			++j;
			dest[j] = '"';
		} else {
			dest[j] = '\'';
			++j;
			dest[j] = src[i];
			++j;
			dest[j] = '\'';
		}
	}
	dest[j] = '\0';
	return dest;
}

/* Append path to global variable MLUCAS_PATH to form mlucas_path,
   which is then passed to fopen()

   Since the length of both MLUCAS_PATH and path are at most STR_MAX_LEN,
   we can use strcpy() and strcat() safely  */
FILE *mlucas_fopen(const char *path, const char *mode)
{
	char mlucas_path[2 * STR_MAX_LEN + 1];

	strcpy(mlucas_path, MLUCAS_PATH);
	strcat(mlucas_path, path);
	return fopen(mlucas_path, mode);
}

