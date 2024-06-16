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

#define RADIX 384	// Use #define rather than const int to ensure it's really a compile-time const in the C sense

#define EPS 1e-10

/* When tried to include radix128.h, kept getting this linker error, even after wrapping contents of said header in include-once code:
	ld: duplicate symbol _DFT128_TWIDDLES in radix384_ditN_cy_dif1.o and dft_macro.o for architecture x86_64
	clang: error: linker command failed with exit code 1 (use -v to see invocation)
...so worked around like so: */
#include "radix128.h"
const double RADIX384_DFT128_TWIDDLES[16][14] = {
	{ 1,0, 1,0, 1,0, 1,0, 1,0, 1,0, 1,0 },
	{ 0,1,ISRT2,ISRT2,-ISRT2,ISRT2,c16,s16,-s16,c16,s16,c16,-c16,s16 },
	{ ISRT2,ISRT2,c16,s16,s16,c16,c32_1,s32_1,s32_3,c32_3,c32_3,s32_3,s32_1,c32_1 },
	{ -ISRT2,ISRT2,s16,c16,-c16,-s16,c32_3,s32_3,-c32_1,s32_1,-s32_1,c32_1,-s32_3,-c32_3 },
	{ c16,s16,c32_1,s32_1,c32_3,s32_3,c64_1,s64_1,c64_5,s64_5,c64_3,s64_3,c64_7,s64_7 },
	{ -s16,c16,s32_3,c32_3,-c32_1,s32_1,c64_5,s64_5,-c64_7,s64_7,s64_1,c64_1,-c64_3,-s64_3 },
	{ s16,c16,c32_3,s32_3,-s32_1,c32_1,c64_3,s64_3,s64_1,c64_1,s64_7,c64_7,-s64_5,c64_5 },
	{ -c16,s16,s32_1,c32_1,-s32_3,-c32_3,c64_7,s64_7,-c64_3,-s64_3,-s64_5,c64_5,s64_1,-c64_1 },
	{ c32_1,s32_1, c64_1,s64_1, c64_3,s64_3, c128_1,s128_1, c128_5,s128_5, c128_3,s128_3, c128_7,s128_7 },
	{ -s32_1,c32_1, s64_7,c64_7, -c64_5,s64_5, c128_9,s128_9, -s128_d,c128_d, s128_5,c128_5, -c128_1,s128_1 },
	{ s32_3,c32_3, c64_5,s64_5, s64_1,c64_1, c128_5,s128_5, s128_7,c128_7, c128_f,s128_f, -s128_3,c128_3 },
	{ -c32_3,s32_3, s64_3,c64_3, -c64_7,-s64_7, c128_d,s128_d, -c128_1,-s128_1, -s128_7,c128_7, -s128_5,-c128_5 },
	{ c32_3,s32_3, c64_3,s64_3, s64_7,c64_7, c128_3,s128_3, c128_f,s128_f, c128_9,s128_9, s128_b,c128_b },
	{ -s32_3,c32_3, s64_5,c64_5, -c64_1,-s64_1, c128_b,s128_b, -c128_9,s128_9, -s128_1,c128_1, -c128_d,-s128_d },
	{ s32_1,c32_1, c64_7,s64_7, -s64_5,c64_5, c128_7,s128_7, -s128_3,c128_3, s128_b,c128_b, -c128_f,s128_f },
	{ -c32_1,s32_1, s64_1,c64_1, -s64_3,-c64_3, c128_f,s128_f, -c128_b,-s128_b, -s128_d,c128_d, s128_9,-c128_9 }
};

#ifndef PFETCH_DIST
  #ifdef USE_AVX512
	#define PFETCH_DIST	64	// Feb 2017: Test on KNL point to this as best
  #elif defined(USE_AVX)
	#define PFETCH_DIST	32	// This seems to work best on my Haswell, even though 64 bytes seems more logical in AVX mode
  #else
	#define PFETCH_DIST	32
  #endif
#endif

#ifdef USE_SSE2

	#define EPS 1e-10
  // RADIX = 384 = 0x180:
  // For Mersenne-mod need (16 [SSE2] or 64 [AVX]) + (4 [HIACC] or 40 [LOACC]) added slots for half_arr lookup tables.
  // Max = (40 [SSE2]; 132 [AVX]), add to (half_arr_offset384 + RADIX) to get SIMD value of radix384_creals_in_local_store:
  #ifdef USE_AVX512	// RADIX/8 = 0x30 fewer carry slots than AVX:
	const int half_arr_offset384 = 0x78a;	// + RADIX + 132 (=0x84); Used for thread local-storage-integrity checking
	const int radix384_creals_in_local_store = 0x990;	// ... + 0x180 + 0x084, and round up to nearest multiple of 4
  #elif defined(USE_AVX)
	const int half_arr_offset384 = 0x7ba;	// + RADIX + 132 (=0x84); Used for thread local-storage-integrity checking
	const int radix384_creals_in_local_store = 0x9b0;	// ... + 0x180 + 0x084, and round up to nearest multiple of 4
  #else
	const int half_arr_offset384 = 0x81a;	// + RADIX + 40 = 0x81a + 0x180 + 0x028; Used for thread local-storage-integrity checking
	const int radix384_creals_in_local_store = 0x9c4;	// ...and round up to nearest multiple of 4
  #endif

	#include "sse2_macro_gcc64.h"

#endif	// SSE2

#ifdef USE_PTHREAD

	// Use non-pooled simple spawn/rejoin thread-team model
	#include "threadpool.h"

	struct cy_thread_data_t{
	// int data:
		int iter;
		int tid;
		int ndivr;
		int target_idx, target_set;	// Jun 2018: Add support for residue shift. (Only LL-test needs intervention at carry-loop level).
		double target_cy;

		int khi;
		int i;
		int jstart;
		int jhi;
		int col;
		int co2;
		int co3;
		int sw;
		int nwt;

	// double data:
		double maxerr;
		double scale;
		double prp_mult;

	// pointer data:
		double *arrdat;			/* Main data array */
		double *wt0;
		double *wt1;
		double *wts_mult, *inv_mult;
		int *si;
	#ifdef USE_SSE2
		vec_dbl *r000;
		vec_dbl *half_arr;
	#else
		double *r000;
		double *half_arr;
	#endif
		uint32 bjmodnini;
		int bjmodn0;
	// For large radix0 use thread-local arrays for DWT indices/carries - only caveat is these must be SIMD-aligned:
	#if GCC_EVER_GETS_ITS_ACT_TOGETHER_HERE
	/* Jan 2014: Bloody hell - turns out GCC uses __BIGGEST_ALIGNMENT__ = 16 on x86, which is too small to be useful for avx data!
		int bjmodn[RADIX] __attribute__ ((aligned (32)));
		double cy[RADIX] __attribute__ ((aligned (32)));
	*/
	#else
	// Thus, we are forced to resort to fugly hackage - add pad slots to a garbage-named struct-internal array along with
	// a pointer-to-be-inited-at-runtime, when we set ptr to the lowest-index array element having the desired alginment:
		double *cy;
	  #ifdef USE_AVX512
		double cy_dat[RADIX+8] __attribute__ ((__aligned__(8)));
	  #else
		double cy_dat[RADIX+4] __attribute__ ((__aligned__(8)));	// Enforce min-alignment of 8 bytes in 32-bit builds.
	  #endif
	#endif
	};

#endif

/**************/

int radix384_ditN_cy_dif1(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], double base[], double baseinv[], int iter, double *fracmax, uint64 p)
{
/*
!...Acronym: DWT = Discrete Weighted Transform, DIT = Decimation In Time, DIF = Decimation In Frequency
!
!...Performs a final radix-384 complex DIT pass, an inverse DWT weighting, a carry propagation,
!   a forward DWT weighting, and an initial radix-384 complex DIF pass on the data in the length-N real vector A.
!
!   Data enter and are returned in the A-array.
!
!   See the documentation in mers_mod_square and radix16_dif_pass for further details on the array
!   storage scheme, and radix16_ditN_cy_dif1 for details on the reduced-length weights array scheme.
*/
	const char func[] = "radix384_ditN_cy_dif1";
	const int stride = (int)RE_IM_STRIDE << 1;	// main-array loop stride = 2*RE_IM_STRIDE
#ifndef MULTITHREAD
	const int pfetch_dist = PFETCH_DIST;

	static int dif_i_offsets[128], dif_o_offsets[384];
	static int dit_i_offsets[384], dit_o_offsets[128];
	static int dif_triplets[72];
	static int dit_triplets[24];	// Only need 1/3 as many here as for DIF
#endif

	static double wts_mult[2], inv_mult[2];	// Const wts-multiplier and 2*(its multiplicative inverse)
	double wt_re,wt_im, wi_re,wi_im;	// Fermat-mod/LOACC weights stuff, used in both scalar and SIMD mode
	// Cleanup loop assumes carryins propagate at most 4 words up, but need at least 1 vec_cmplx
	// (2 vec_dbl)'s worth of doubles in wraparound step, hence AVX-512 needs value bumped up:
  #ifdef USE_AVX512
	const int jhi_wrap = 15;
  #else
	const int jhi_wrap =  7;
  #endif
	int NDIVR,i,incr,j,j1,j2,jt,jp,jstart,jhi,full_pass,k,khi,l,l1,l2,outer,nbytes;
	// incr = Carry-chain wts-multipliers recurrence length, which must divide
	// RADIX/[n-wayness of carry macro], e.g. RADIX/[16|8|4] = 48|96|192 for avx512,avx,sse, respectively:
	const int incr_long = 16,incr_med = 8,incr_short = 4;
	// Allows cy-macro error data to be used to fiddle incr on the fly to a smaller, safer value if necessary
	if(USE_SHORT_CY_CHAIN == 0)
		incr = incr_long;
	else if(USE_SHORT_CY_CHAIN == 1)
		incr = incr_med;
	else
		incr = incr_short;
	int k0,k1,k2;
	// Jun 2018: Add support for residue shift. (Only LL-test needs intervention at carry-loop level).
	int target_idx = -1, target_set = 0,tidx_mod_stride;
	double target_cy = 0;
	static double ndivr_inv;
	uint64 itmp64;
	static uint64 psave = 0;
	static uint32 bw,sw,bjmodnini,
		p1,p2,p3,p4,p5,p6,p7,p8,p9,pa,pb,pc,pd,pe,pf,
		p010,p020,p030,p040,p050,p060,p070,p080,p090,p0a0,p0b0,p0c0,p0d0,p0e0,p0f0,p100,p110,p120,p130,p140,p150,p160,p170, nsave = 0;
	static int poff[RADIX>>2];	// Store mults of p4 offset for loop control
	static double radix_inv, n2inv;
#if defined(USE_SSE2) || !defined(MULTITHREAD)
	const double	c3m1= -1.50000000000000000000,	/* cos(twopi/3)-1	*/
					s   =  0.86602540378443864675;
#endif
	double scale, dtmp, maxerr = 0.0;
	// Local storage: We must use an array here because scalars have no guarantees about relative address offsets
	// [and even if those are contiguous-as-hoped-for, they may run in reverse]; Make array type (struct complex)
	// to allow us to use the same offset-indexing as in the original radix-32 in-place DFT macros:
	const double *addr, *addi;
	struct complex t[RADIX], *tptr;
	int err;
	static int first_entry=TRUE;

/*...stuff for the reduced-length DWT weights array is here:	*/
	int n_div_nwt;
  #ifdef USE_AVX512
	double t0,t1,t2,t3;
   #ifdef CARRY_16_WAY
	static struct uint32x16 *n_minus_sil,*n_minus_silp1,*sinwt,*sinwtm1;
   #else
	static struct uint32x8  *n_minus_sil,*n_minus_silp1,*sinwt,*sinwtm1;
   #endif
  #elif defined(USE_AVX)
	static struct uint32x4 *n_minus_sil,*n_minus_silp1,*sinwt,*sinwtm1;
  #endif
#ifndef MULTITHREAD
	int col,co2,co3;
	double wtl,wtlp1,wtn,wtnm1;	/* Mersenne-mod weights stuff */
  #ifndef USE_AVX
	int n_minus_sil,n_minus_silp1,sinwt,sinwtm1;
  #endif
#endif

#ifdef USE_SSE2

  #if !(defined(COMPILER_TYPE_MSVC) || defined(COMPILER_TYPE_GCC))
	#error SSE2 code not supported for this compiler!
  #endif

	static int cslots_in_local_store;
	static vec_dbl *sc_arr = 0x0, *sc_ptr;
	static uint64 *sm_ptr, *sign_mask, *sse_bw, *sse_sw, *sse_n;
	uint64 tmp64;

  #ifdef MULTITHREAD
	static vec_dbl *__r0;	// Base address for discrete per-thread local stores
  #else
	double *add0,*add1,*add2,*add3;
	int *itmp,*itm2;			// Pointer into the bjmodn array
  #endif

	// Uint64 bitmaps for alternate "rounded the other way" copies of sqrt2,isrt2. Default round-to-nearest versions
	// (SQRT2, ISRT2) end in ...3BCD. Since we round these down as ...3BCC90... --> ..3BCC, append _dn to varnames:
	const uint64 sqrt2_dn = 0x3FF6A09E667F3BCCull, isrt2_dn = 0x3FE6A09E667F3BCCull;
	static int *bjmodn;	// Alloc mem for this along with other 	SIMD stuff
	const double crnd = 3.0*0x4000000*0x2000000;
	struct complex *ctmp;	// Hybrid AVX-DFT/SSE2-carry scheme used for Mersenne-mod needs a 2-word-double pointer
	vec_dbl *tmp,*tm1,*tm2;	// Non-static utility ptrs
	static vec_dbl *two,*one,*sqrt2,*isrt2, *cc0, *ss0, *cc1, *ss1, *max_err, *sse2_rnd, *half_arr,
		// ptrs to 16 sets of twiddles shared by the 2nd-half DIF and DIT DFT macros:
		*twid0,*twid1,*twid2,*twid3,*twid4,*twid5,*twid6,*twid7,*twid8,*twid9,*twida,*twidb,*twidc,*twidd,*twide,*twidf,
		*r000,*r080,*r100,	// Head of RADIX*vec_cmplx-sized local store #1
		*s1p000,	// Head of RADIX*vec_cmplx-sized local store #2
		*cy;	// Need RADIX/2 slots for sse2 carries, RADIX/4 for avx
#else
	static int p0123[4];
#endif	// USE_SSE2?

#ifdef MULTITHREAD

	static struct cy_thread_data_t *tdat = 0x0;
	// Threadpool-based dispatch stuff:
	static int main_work_units = 0, pool_work_units = 0;
	static struct threadpool *tpool = 0x0;
	static int task_is_blocking = TRUE;
	static thread_control_t thread_control = {0,0,0};
	// First 3 subfields same for all threads, 4th provides thread-specifc data, will be inited at thread dispatch:
	static task_control_t   task_control = {NULL, (void*)cy384_process_chunk, NULL, 0x0};

#elif !defined(USE_SSE2)

	// Vars needed in scalar mode only:
	const double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
	int m,m2,ntmp;
	double wt,wtinv,wtA,wtB,wtC;	/* Mersenne-mod weights stuff */
  #if PFETCH
	int prefetch_offset;
  #endif
	int bjmodn[RADIX];
	int *itmp,*itm2;	// Pointer into the bjmodn array
	double temp,frac,cy[RADIX],
		t00,t01,t02,t03,t04,t05;

#endif

/*...stuff for the multithreaded implementation is here:	*/
	static uint32 CY_THREADS,pini;
	int ithread,j_jhi;
	uint32 ptr_prod;
	static int *_bjmodnini = 0x0,*_bjmodn[RADIX];
	static int *_i, *_jstart = 0x0, *_jhi = 0x0, *_col = 0x0, *_co2 = 0x0, *_co3 = 0x0;
	static double *_cy[RADIX];
	if(!_jhi) {
		_cy[0] = 0x0;	// First of these used as an "already inited consts?" sentinel, must init = 0x0 at same time do so for non-array static ptrs
	}

	if(MODULUS_TYPE == MODULUS_TYPE_FERMAT)
	{
		ASSERT(0, "Fermat-mod only available for radices 7,8,9,15 and their multiples!");
	}

	// Jan 2018: To support PRP-testing, read the LR-modpow-scalar-multiply-needed bit for the current iteration from the global array:
	double prp_mult = 1.0;
	if((TEST_TYPE & 0xfffffffe) == TEST_TYPE_PRP) {	// Mask off low bit to lump together PRP and PRP-C tests
		i = (iter-1) % ITERS_BETWEEN_CHECKPOINTS;	// Bit we need to read...iter-counter is unit-offset w.r.to iter-interval, hence the -1
		if((BASE_MULTIPLIER_BITS[i>>6] >> (i&63)) & 1)
			prp_mult = PRP_BASE;
	}

/*...change NDIVR and n_div_wt to non-static to work around a gcc compiler bug. */
	NDIVR   = n/RADIX;	ndivr_inv = (double)RADIX/n;
	n_div_nwt = NDIVR >> nwt_bits;

	if((n_div_nwt << nwt_bits) != NDIVR)
	{
		fprintf(stderr,"ERROR: iter = %10d; NWT_BITS does not divide N/RADIX in %s.\n",iter,func);
		err = ERR_SKIP_RADIX_SET;
		return(err);
	}

	if(p != psave || n != nsave
	#ifdef USE_PTHREAD	// Oct 2021: cf. radix176_ditN_cy_dif1.c for why I added this
		|| (tdat != 0x0 && tdat[0].wt1 != wt1)
	#endif
	) {	/* Exponent or array length change triggers re-init */
		first_entry=TRUE;
		/* To-do: Support #thread change here! */
	}

/*...initialize things upon first entry: */

	if(first_entry)
	{
		psave = p;	nsave = n;
		radix_inv = qfdbl(qf_rational_quotient((int64)1, (int64)RADIX));
		n2inv     = qfdbl(qf_rational_quotient((int64)1, (int64)(n/2)));

		bw    = p%n;	/* Number of bigwords in the Crandall/Fagin mixed-radix representation = (Mersenne exponent) mod (vector length).	*/
		sw    = n - bw;	/* Number of smallwords.	*/

	  #ifdef USE_AVX512
	   #ifdef CARRY_16_WAY
		i = 16;
	   #else
		i = 8;
	   #endif
	  #elif defined(USE_AVX)	// AVX LOACC: Make CARRY_8_WAY default here:
		i = 8;
	  #elif defined(USE_SSE2)	// AVX and SSE2 modes use 4-way carry macros
		i = 4;
	  #else	// Scalar-double mode:
		i = 1;
	  #endif

		// For n a power of 2 don't need to worry about 32-bit integer overflow in the sw*NDIVR term,
		// but for non-power-of-2 n we must cast-to-uint64 to avoid such overflows fubaring the result:
		struct qfloat qt,qn;
		qt = i64_to_q(i*(uint64)sw*NDIVR % n);
		qn = i64_to_q((int64) n);
		qt = qfdiv(qt, qn);		// x = (sw*NDIVR (mod n))/n
		qt = qfmul(qt, QLN2);	// x*ln(2)...
		qt = qfexp(qt);			// ...and get 2^x via exp[x*ln(2)].
		wts_mult[0] = qfdbl(qt);		// a = 2^(x/n), with x = sw
		inv_mult[0] = qfdbl(qfinv(qt));	// Double-based inversion (1.0 / wts_mult_a[0]) often gets LSB wrong
		ASSERT(fabs(wts_mult[0]*inv_mult[0] - 1.0) < EPS, "wts_mults fail accuracy check!");
		//curr have w, 2/w, separate-mul-by-1-or-0.5 gives [w,w/2] and [1/w,2/w] for i = 0,1, resp:
		wts_mult[1] = 0.5*wts_mult[0];
		inv_mult[1] = 2.0*inv_mult[0];
		ASSERT(fabs(wts_mult[1]*inv_mult[1] - 1.0) < EPS, "wts_mults fail accuracy check!");

	#ifdef MULTITHREAD

		/* #Chunks ||ized in carry step is ideally a power of 2, so use the largest
		power of 2 that is <= the value of the global NTHREADS (but still <= MAX_THREADS):
		*/
		if(isPow2(NTHREADS))
			CY_THREADS = NTHREADS;
		else {
			i = leadz32(NTHREADS);
			CY_THREADS = (((uint32)NTHREADS << i) & 0x80000000) >> i;
		}

		if(CY_THREADS > MAX_THREADS)
		{
		//	CY_THREADS = MAX_THREADS;
			fprintf(stderr,"WARN: CY_THREADS = %d exceeds number of cores = %d\n", CY_THREADS, MAX_THREADS);
		}
		if(!isPow2(CY_THREADS))		{ WARN(HERE, "CY_THREADS not a power of 2!", "", 1); return(ERR_ASSERT); }
		if(CY_THREADS > 1)
		{
			if(NDIVR    %CY_THREADS != 0) { WARN(HERE, "NDIVR    %CY_THREADS != 0 ... likely more threads than this leading radix can handle.", "", 1); return(ERR_ASSERT); }
			if(n_div_nwt%CY_THREADS != 0) { WARN(HERE, "n_div_nwt%CY_THREADS != 0 ... likely more threads than this leading radix can handle.", "", 1); return(ERR_ASSERT); }
		}

	  #ifdef USE_PTHREAD
		if(tdat == 0x0) {
			j = (uint32)sizeof(struct cy_thread_data_t);
			tdat = (struct cy_thread_data_t *)calloc(CY_THREADS, sizeof(struct cy_thread_data_t));

			// MacOS does weird things with threading (e.g. Idle" main thread burning 100% of 1 CPU)
			// so on that platform try to be clever and interleave main-thread and threadpool-work processing
			#if 0//def OS_TYPE_MACOSX

				if(CY_THREADS > 1) {
					main_work_units = CY_THREADS/2;
					pool_work_units = CY_THREADS - main_work_units;
					ASSERT(0x0 != (tpool = threadpool_init(pool_work_units, MAX_THREADS, pool_work_units, &thread_control)), "threadpool_init failed!");
					printf("radix%d_ditN_cy_dif1: Init threadpool of %d threads\n", RADIX, pool_work_units);
				} else {
					main_work_units = 1;
					printf("radix%d_ditN_cy_dif1: CY_THREADS = 1: Using main execution thread, no threadpool needed.\n", RADIX);
				}

			#else

				pool_work_units = CY_THREADS;
				ASSERT(0x0 != (tpool = threadpool_init(CY_THREADS, MAX_THREADS, CY_THREADS, &thread_control)), "threadpool_init failed!");

			#endif

			fprintf(stderr,"Using %d threads in carry step\n", CY_THREADS);
		}
	  #endif

	#else
		CY_THREADS = 1;
	#endif

	#ifdef USE_PTHREAD
		/* Populate the elements of the thread-specific data structs which don't change after init: */
		for(ithread = 0; ithread < CY_THREADS; ithread++)
		{
			tdat[ithread].iter = iter;
		// int data:
			tdat[ithread].tid = ithread;
			tdat[ithread].ndivr = NDIVR;

			tdat[ithread].sw  = sw;
			tdat[ithread].nwt = nwt;

		// pointer data:
		//	tdat[ithread].arrdat = a;			/* Main data array */
			tdat[ithread].wt0 = wt0;
			tdat[ithread].wt1 = wt1;
			tdat[ithread].wts_mult = wts_mult;
			tdat[ithread].inv_mult = inv_mult;
			tdat[ithread].si  = si;

		// This array pointer must be set based on vec_dbl-sized alignment at runtime for each thread:
			for(l = 0; l < RE_IM_STRIDE; l++) {
				if( ((intptr_t)&tdat[ithread].cy_dat[l] & SZ_VDM1) == 0 ) {
					tdat[ithread].cy = &tdat[ithread].cy_dat[l];
				//	fprintf(stderr,"%d-byte-align cy_dat array at element[%d]\n",SZ_VD,l);
					break;
				}
			}
			ASSERT(l < RE_IM_STRIDE, "Failed to align cy_dat array!");
		}
	#endif

	#ifdef USE_SSE2

		ASSERT(((intptr_t)wt0    & 0x3f) == 0, "wt0[]  not 64-byte aligned!");
		ASSERT(((intptr_t)wt1    & 0x3f) == 0, "wt1[]  not 64-byte aligned!");

		// Use double-complex type size (16 bytes) to alloc a block of local storage
		// consisting of radix384_creals_in_local_store dcomplex and (12+RADIX/2) uint64 element slots per thread
		// (Add as many padding elts to the latter as needed to make it a multiple of 4):
		cslots_in_local_store = radix384_creals_in_local_store + (((12+RADIX/2)/2 + 3) & ~0x3);
		sc_arr = ALLOC_VEC_DBL(sc_arr, cslots_in_local_store*CY_THREADS);	if(!sc_arr){ sprintf(cbuf, "ERROR: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		sc_ptr = ALIGN_VEC_DBL(sc_arr);
		ASSERT(((intptr_t)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");
		sm_ptr = (uint64*)(sc_ptr + radix384_creals_in_local_store);
		ASSERT(((intptr_t)sm_ptr & 0x3f) == 0, "sm_ptr not 64-byte aligned!");

	  #ifdef USE_PTHREAD
		__r0 = sc_ptr;
	  #endif
		tmp = r000 = sc_ptr;
		r080 = tmp + 0x100;
		r100 = tmp + 0x200;
		s1p000 = tmp + 0x300;
		tmp += 0x600;	// sc_ptr += 0x600
		two     = tmp + 0;	// AVX+ versions of radix-8,16,32 twiddleless-DFT macros need consts [2,1,sqrt2,isrt2] quartet laid out thusly
		one     = tmp + 1;
		sqrt2   = tmp + 2;
		isrt2   = tmp + 3;
		cc0		= tmp + 4;
		ss0		= tmp + 5;
		cc1		= tmp + 6;
		ss1		= tmp + 7;
		tmp += 0x08;	// sc_ptr += 0x608
		// Each non-unity root now needs a negated counterpart:
		/* Stupidity: Since a truly general-purpose [in the sense that it can be used for our radix-128 internal-twiddles]
		radix-8 DFT-with-twiddles macro needs 8 in-addresses [corr. to the 8 real parts of the input data], 8 o-addresses,
		and 7 each of cosine and sine data [which cannot be assumed to occur in fixed-stride pairs - cf. our usage of
		SSE2_RADIX8_DIT_TWIDDLE_OOP() below], that hits the GCC hard limit of 30-operands for ASM macros, but we still
		need one more operand for the ISRT2 pointer. Only easy workaround I found for this is to stick a vector-ISRT2 copy
		in between each +-[cc,ss] vector-data pair, thus any time we need a vector-isrt2 for the radix-8 internal twiddles
		we get it at (vec_dbl*)cc-1.	/Stupidity */
		// ptrs to 16 sets (3*7 = 21 vec_dbl data each) of non-unity twiddles shared by the 2nd-half DIF and DIT DFT macros:
		twid0  = tmp + 0x00;
		twid1  = tmp + 0x15;
		twid2  = tmp + 0x2a;
		twid3  = tmp + 0x3f;
		twid4  = tmp + 0x54;
		twid5  = tmp + 0x69;
		twid6  = tmp + 0x7e;
		twid7  = tmp + 0x93;
		twid8  = tmp + 0xa8;
		twid9  = tmp + 0xbd;
		twida  = tmp + 0xd2;
		twidb  = tmp + 0xe7;
		twidc  = tmp + 0xfc;
		twidd  = tmp + 0x111;
		twide  = tmp + 0x126;
		twidf  = tmp + 0x13b;
		tmp += 0x150;	// += 16*21 = 0x150 => sc_ptr + 0x758
	  #ifdef USE_AVX512
		cy = tmp;		tmp += 0x30;	// RADIX/8 vec_dbl slots for carry sub-array
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;
		half_arr= tmp + 0x02;	// += 0x30 + 2 => sc_ptr += 0x78a
	  #elif defined(USE_AVX)
		cy = tmp;		tmp += 0x60;	// RADIX/4 vec_dbl slots for carry sub-array
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// += 0x60 + 2 => sc_ptr += 0x7ba
		// This is where the value of half_arr_offset comes from
		half_arr= tmp + 0x02;	/* This table needs 96 vec_dbl for Mersenne-mod, and 3.5*RADIX[avx] | RADIX[sse2] for Fermat-mod */
	  #else
		cy = tmp;		tmp += 0xc0;	// RADIX/2 vec_dbl slots for carry sub-array
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// += 0xc0 + 2 => sc_ptr += 0x81a
		// This is where the value of half_arr_offset comes from
		half_arr= tmp + 0x02;	/* This table needs 32 x 16 bytes for Mersenne-mod, 2 for Fermat-mod */
	  #endif
//		ASSERT(half_arr_offset == (uint32)(half_arr-sc_ptr), "half_arr_offset mismatches actual!");
		ASSERT((radix384_creals_in_local_store << L2_SZ_VD) >= ((intptr_t)half_arr - (intptr_t)r000) + (20 << L2_SZ_VD), "radix384_creals_in_local_store checksum failed!");

		/* These remain fixed: */
		VEC_DBL_INIT(two  , 2.0  );	VEC_DBL_INIT(one, 1.0  );
		dtmp = *(double *)&sqrt2_dn;	VEC_DBL_INIT(sqrt2, dtmp);
		dtmp = *(double *)&isrt2_dn;	VEC_DBL_INIT(isrt2, dtmp);
		VEC_DBL_INIT(cc0  ,  c16);
		VEC_DBL_INIT(ss0  ,  s16);
		/* cc0 = (cc1+cc2+cc3)/3 - 1; subtract 1 from Nussbaumer's definition in order to ease in-place computation */
		VEC_DBL_INIT(cc1, c3m1);
		VEC_DBL_INIT(ss1, s   );
		// Use loop-based init for code compactness and extensibility to larger pow2 radices:
		for(l = 0; l < 16; l++) {
			j = reverse(l,4);	// twid0-offsets are processed in BR16 order
			tmp = twid0 + (j<<4)+(j<<2)+j;	// Twid-offsets are multiples of 21 vec_dbl
			addr = RADIX384_DFT128_TWIDDLES[l]; addi = addr+1;	// Pointer to required row of 2-D twiddles array
			VEC_DBL_INIT(tmp, ISRT2); tmp++;	VEC_DBL_INIT(tmp,*(addr+0x00)); tmp++;	VEC_DBL_INIT(tmp,*(addi+0x00)); tmp++;
			VEC_DBL_INIT(tmp, ISRT2); tmp++;	VEC_DBL_INIT(tmp,*(addr+0x02)); tmp++;	VEC_DBL_INIT(tmp,*(addi+0x02)); tmp++;
			VEC_DBL_INIT(tmp, ISRT2); tmp++;	VEC_DBL_INIT(tmp,*(addr+0x04)); tmp++;	VEC_DBL_INIT(tmp,*(addi+0x04)); tmp++;
			VEC_DBL_INIT(tmp, ISRT2); tmp++;	VEC_DBL_INIT(tmp,*(addr+0x06)); tmp++;	VEC_DBL_INIT(tmp,*(addi+0x06)); tmp++;
			VEC_DBL_INIT(tmp, ISRT2); tmp++;	VEC_DBL_INIT(tmp,*(addr+0x08)); tmp++;	VEC_DBL_INIT(tmp,*(addi+0x08)); tmp++;
			VEC_DBL_INIT(tmp, ISRT2); tmp++;	VEC_DBL_INIT(tmp,*(addr+0x0a)); tmp++;	VEC_DBL_INIT(tmp,*(addi+0x0a)); tmp++;
			VEC_DBL_INIT(tmp, ISRT2); tmp++;	VEC_DBL_INIT(tmp,*(addr+0x0c)); tmp++;	VEC_DBL_INIT(tmp,*(addi+0x0c)); tmp++;
		}
		/* SSE2 math = 53-mantissa-bit IEEE double-float: */
	  #ifdef USE_AVX512	// In AVX-512 mode, use VRNDSCALEPD for rounding and hijack this vector-data slot for the 4 base/baseinv-consts
		sse2_rnd->d0 = base[0]; sse2_rnd->d1 = baseinv[1]; sse2_rnd->d2 = wts_mult[1]; sse2_rnd->d3 = inv_mult[0];
	  #else
		VEC_DBL_INIT(sse2_rnd, crnd);
	  #endif

		// Propagate the above consts to the remaining threads:
		nbytes = (intptr_t)cy - (intptr_t)two;	// #bytes in 1st of above block of consts
		tmp = two;
		tm2 = tmp + cslots_in_local_store;
		for(ithread = 1; ithread < CY_THREADS; ++ithread) {
			memcpy(tm2, tmp, nbytes);
			tmp = tm2;		tm2 += cslots_in_local_store;
		}
		nbytes = SZ_VD;	// sse2_rnd is a solo (in the SIMD-vector) datum
		tmp = sse2_rnd;
		tm2 = tmp + cslots_in_local_store;
		for(ithread = 1; ithread < CY_THREADS; ++ithread) {
			memcpy(tm2, tmp, nbytes);
			tmp = tm2;		tm2 += cslots_in_local_store;
		}

		/* SSE2 version of the one_half array - we have a 2-bit lookup, low bit is from the low word of the carry pair,
		high bit from the high, i.e. based on this lookup index [listed with LSB at right], we have:

			index	half_lo	half_hi
			00		1.0		1.0
			01		.50		1.0
			10		1.0		.50
			11		.50		.50

		The inverse-weights computation uses a similar table, but with all entries multiplied by .50:

			index2	half_lo	half_hi
			00		.50		.50
			01		.25		.50
			10		.50		.25
			11		.25		.25

		We do similarly for the base[] and baseinv[] table lookups - each of these get 4 further slots in half_arr.
		We also allocate a further 4 16-byte slots [uninitialized] for storage of the wtl,wtn,wtlp1,wtnm1 locals.

		In 4-way SIMD (AVX) mode, we expand this from 2^2 2-vector table entries to 2^4 4-vector entries.
		*/
		tmp = half_arr;

	#ifdef USE_AVX512
		// Each lookup-category in the 'mini-tables' used in AVX mode balloons from 16x32-bytes to 64x64-bytes,
		// so switch to an opmask-based scheme which starts with e.g. a broadcast constant and onditional doubling.
		// Here are the needed consts and opmasks:
		// [1] Fwd-wt multipliers: Init = 0.50 x 8, anytime AVX-style lookup into 1st table below would have bit = 0, double the corr. datum
		// [2] Inv-wt multipliers: Init = 0.25 x 8, anytime AVX-style lookup into 2nd table below would have bit = 0, double the corr. datum
		// [3] Fwd-base mults: Init = base[0] x 8, anytime AVX-style lookup into 3rd table below would have bit = 1, double the corr. datum
		// [4] Inv-base mults: Init = binv[1] x 8, anytime AVX-style lookup into 4th table below would have bit = 0, double the corr. datum
		// [5] [LOACC] Init = wts_mult[1] x 8, anytime AVX-style lookup into 5th table below would have bit = 0, double the corr. datum
		// [6] [LOACC] Init = inv_mult[0] x 8, anytime AVX-style lookup into 6th table below would have bit = 1, double the corr. datum
		nbytes = 0;
	#elif defined(USE_AVX)
		/* Forward-weight multipliers: */
		tmp->d0 = 1.0;	tmp->d1 = 1.0;	tmp->d2 = 1.0;	tmp->d3 = 1.0;	++tmp;
		tmp->d0 = .50;	tmp->d1 = 1.0;	tmp->d2 = 1.0;	tmp->d3 = 1.0;	++tmp;
		tmp->d0 = 1.0;	tmp->d1 = .50;	tmp->d2 = 1.0;	tmp->d3 = 1.0;	++tmp;
		tmp->d0 = .50;	tmp->d1 = .50;	tmp->d2 = 1.0;	tmp->d3 = 1.0;	++tmp;
		tmp->d0 = 1.0;	tmp->d1 = 1.0;	tmp->d2 = .50;	tmp->d3 = 1.0;	++tmp;
		tmp->d0 = .50;	tmp->d1 = 1.0;	tmp->d2 = .50;	tmp->d3 = 1.0;	++tmp;
		tmp->d0 = 1.0;	tmp->d1 = .50;	tmp->d2 = .50;	tmp->d3 = 1.0;	++tmp;
		tmp->d0 = .50;	tmp->d1 = .50;	tmp->d2 = .50;	tmp->d3 = 1.0;	++tmp;
		tmp->d0 = 1.0;	tmp->d1 = 1.0;	tmp->d2 = 1.0;	tmp->d3 = .50;	++tmp;
		tmp->d0 = .50;	tmp->d1 = 1.0;	tmp->d2 = 1.0;	tmp->d3 = .50;	++tmp;
		tmp->d0 = 1.0;	tmp->d1 = .50;	tmp->d2 = 1.0;	tmp->d3 = .50;	++tmp;
		tmp->d0 = .50;	tmp->d1 = .50;	tmp->d2 = 1.0;	tmp->d3 = .50;	++tmp;
		tmp->d0 = 1.0;	tmp->d1 = 1.0;	tmp->d2 = .50;	tmp->d3 = .50;	++tmp;
		tmp->d0 = .50;	tmp->d1 = 1.0;	tmp->d2 = .50;	tmp->d3 = .50;	++tmp;
		tmp->d0 = 1.0;	tmp->d1 = .50;	tmp->d2 = .50;	tmp->d3 = .50;	++tmp;
		tmp->d0 = .50;	tmp->d1 = .50;	tmp->d2 = .50;	tmp->d3 = .50;	++tmp;
		/* Inverse-weight multipliers (only needed for mersenne-mod): */
		tmp->d0 = .50;	tmp->d1 = .50;	tmp->d2 = .50;	tmp->d3 = .50;	++tmp;
		tmp->d0 = .25;	tmp->d1 = .50;	tmp->d2 = .50;	tmp->d3 = .50;	++tmp;
		tmp->d0 = .50;	tmp->d1 = .25;	tmp->d2 = .50;	tmp->d3 = .50;	++tmp;
		tmp->d0 = .25;	tmp->d1 = .25;	tmp->d2 = .50;	tmp->d3 = .50;	++tmp;
		tmp->d0 = .50;	tmp->d1 = .50;	tmp->d2 = .25;	tmp->d3 = .50;	++tmp;
		tmp->d0 = .25;	tmp->d1 = .50;	tmp->d2 = .25;	tmp->d3 = .50;	++tmp;
		tmp->d0 = .50;	tmp->d1 = .25;	tmp->d2 = .25;	tmp->d3 = .50;	++tmp;
		tmp->d0 = .25;	tmp->d1 = .25;	tmp->d2 = .25;	tmp->d3 = .50;	++tmp;
		tmp->d0 = .50;	tmp->d1 = .50;	tmp->d2 = .50;	tmp->d3 = .25;	++tmp;
		tmp->d0 = .25;	tmp->d1 = .50;	tmp->d2 = .50;	tmp->d3 = .25;	++tmp;
		tmp->d0 = .50;	tmp->d1 = .25;	tmp->d2 = .50;	tmp->d3 = .25;	++tmp;
		tmp->d0 = .25;	tmp->d1 = .25;	tmp->d2 = .50;	tmp->d3 = .25;	++tmp;
		tmp->d0 = .50;	tmp->d1 = .50;	tmp->d2 = .25;	tmp->d3 = .25;	++tmp;
		tmp->d0 = .25;	tmp->d1 = .50;	tmp->d2 = .25;	tmp->d3 = .25;	++tmp;
		tmp->d0 = .50;	tmp->d1 = .25;	tmp->d2 = .25;	tmp->d3 = .25;	++tmp;
		tmp->d0 = .25;	tmp->d1 = .25;	tmp->d2 = .25;	tmp->d3 = .25;	++tmp;
		/* Forward-base[] multipliers: */
		tmp->d0 = base   [0];	tmp->d1 = base   [0];	tmp->d2 = base   [0];	tmp->d3 = base   [0];	++tmp;
		tmp->d0 = base   [1];	tmp->d1 = base   [0];	tmp->d2 = base   [0];	tmp->d3 = base   [0];	++tmp;
		tmp->d0 = base   [0];	tmp->d1 = base   [1];	tmp->d2 = base   [0];	tmp->d3 = base   [0];	++tmp;
		tmp->d0 = base   [1];	tmp->d1 = base   [1];	tmp->d2 = base   [0];	tmp->d3 = base   [0];	++tmp;
		tmp->d0 = base   [0];	tmp->d1 = base   [0];	tmp->d2 = base   [1];	tmp->d3 = base   [0];	++tmp;
		tmp->d0 = base   [1];	tmp->d1 = base   [0];	tmp->d2 = base   [1];	tmp->d3 = base   [0];	++tmp;
		tmp->d0 = base   [0];	tmp->d1 = base   [1];	tmp->d2 = base   [1];	tmp->d3 = base   [0];	++tmp;
		tmp->d0 = base   [1];	tmp->d1 = base   [1];	tmp->d2 = base   [1];	tmp->d3 = base   [0];	++tmp;
		tmp->d0 = base   [0];	tmp->d1 = base   [0];	tmp->d2 = base   [0];	tmp->d3 = base   [1];	++tmp;
		tmp->d0 = base   [1];	tmp->d1 = base   [0];	tmp->d2 = base   [0];	tmp->d3 = base   [1];	++tmp;
		tmp->d0 = base   [0];	tmp->d1 = base   [1];	tmp->d2 = base   [0];	tmp->d3 = base   [1];	++tmp;
		tmp->d0 = base   [1];	tmp->d1 = base   [1];	tmp->d2 = base   [0];	tmp->d3 = base   [1];	++tmp;
		tmp->d0 = base   [0];	tmp->d1 = base   [0];	tmp->d2 = base   [1];	tmp->d3 = base   [1];	++tmp;
		tmp->d0 = base   [1];	tmp->d1 = base   [0];	tmp->d2 = base   [1];	tmp->d3 = base   [1];	++tmp;
		tmp->d0 = base   [0];	tmp->d1 = base   [1];	tmp->d2 = base   [1];	tmp->d3 = base   [1];	++tmp;
		tmp->d0 = base   [1];	tmp->d1 = base   [1];	tmp->d2 = base   [1];	tmp->d3 = base   [1];	++tmp;
		/* Inverse-base[] multipliers: */
		tmp->d0 = baseinv[0];	tmp->d1 = baseinv[0];	tmp->d2 = baseinv[0];	tmp->d3 = baseinv[0];	++tmp;
		tmp->d0 = baseinv[1];	tmp->d1 = baseinv[0];	tmp->d2 = baseinv[0];	tmp->d3 = baseinv[0];	++tmp;
		tmp->d0 = baseinv[0];	tmp->d1 = baseinv[1];	tmp->d2 = baseinv[0];	tmp->d3 = baseinv[0];	++tmp;
		tmp->d0 = baseinv[1];	tmp->d1 = baseinv[1];	tmp->d2 = baseinv[0];	tmp->d3 = baseinv[0];	++tmp;
		tmp->d0 = baseinv[0];	tmp->d1 = baseinv[0];	tmp->d2 = baseinv[1];	tmp->d3 = baseinv[0];	++tmp;
		tmp->d0 = baseinv[1];	tmp->d1 = baseinv[0];	tmp->d2 = baseinv[1];	tmp->d3 = baseinv[0];	++tmp;
		tmp->d0 = baseinv[0];	tmp->d1 = baseinv[1];	tmp->d2 = baseinv[1];	tmp->d3 = baseinv[0];	++tmp;
		tmp->d0 = baseinv[1];	tmp->d1 = baseinv[1];	tmp->d2 = baseinv[1];	tmp->d3 = baseinv[0];	++tmp;
		tmp->d0 = baseinv[0];	tmp->d1 = baseinv[0];	tmp->d2 = baseinv[0];	tmp->d3 = baseinv[1];	++tmp;
		tmp->d0 = baseinv[1];	tmp->d1 = baseinv[0];	tmp->d2 = baseinv[0];	tmp->d3 = baseinv[1];	++tmp;
		tmp->d0 = baseinv[0];	tmp->d1 = baseinv[1];	tmp->d2 = baseinv[0];	tmp->d3 = baseinv[1];	++tmp;
		tmp->d0 = baseinv[1];	tmp->d1 = baseinv[1];	tmp->d2 = baseinv[0];	tmp->d3 = baseinv[1];	++tmp;
		tmp->d0 = baseinv[0];	tmp->d1 = baseinv[0];	tmp->d2 = baseinv[1];	tmp->d3 = baseinv[1];	++tmp;
		tmp->d0 = baseinv[1];	tmp->d1 = baseinv[0];	tmp->d2 = baseinv[1];	tmp->d3 = baseinv[1];	++tmp;
		tmp->d0 = baseinv[0];	tmp->d1 = baseinv[1];	tmp->d2 = baseinv[1];	tmp->d3 = baseinv[1];	++tmp;
		tmp->d0 = baseinv[1];	tmp->d1 = baseinv[1];	tmp->d2 = baseinv[1];	tmp->d3 = baseinv[1];	++tmp;
		// In LOACC mode, put wts_mult and their inverses in the first 32 slots below in place of the 1/2-stuff:
		/* wts_mult:*/
		tmp->d0 = wts_mult[0];	tmp->d1 = wts_mult[0];	tmp->d2 = wts_mult[0];	tmp->d3 = wts_mult[0];	++tmp;
		tmp->d0 = wts_mult[1];	tmp->d1 = wts_mult[0];	tmp->d2 = wts_mult[0];	tmp->d3 = wts_mult[0];	++tmp;
		tmp->d0 = wts_mult[0];	tmp->d1 = wts_mult[1];	tmp->d2 = wts_mult[0];	tmp->d3 = wts_mult[0];	++tmp;
		tmp->d0 = wts_mult[1];	tmp->d1 = wts_mult[1];	tmp->d2 = wts_mult[0];	tmp->d3 = wts_mult[0];	++tmp;
		tmp->d0 = wts_mult[0];	tmp->d1 = wts_mult[0];	tmp->d2 = wts_mult[1];	tmp->d3 = wts_mult[0];	++tmp;
		tmp->d0 = wts_mult[1];	tmp->d1 = wts_mult[0];	tmp->d2 = wts_mult[1];	tmp->d3 = wts_mult[0];	++tmp;
		tmp->d0 = wts_mult[0];	tmp->d1 = wts_mult[1];	tmp->d2 = wts_mult[1];	tmp->d3 = wts_mult[0];	++tmp;
		tmp->d0 = wts_mult[1];	tmp->d1 = wts_mult[1];	tmp->d2 = wts_mult[1];	tmp->d3 = wts_mult[0];	++tmp;
		tmp->d0 = wts_mult[0];	tmp->d1 = wts_mult[0];	tmp->d2 = wts_mult[0];	tmp->d3 = wts_mult[1];	++tmp;
		tmp->d0 = wts_mult[1];	tmp->d1 = wts_mult[0];	tmp->d2 = wts_mult[0];	tmp->d3 = wts_mult[1];	++tmp;
		tmp->d0 = wts_mult[0];	tmp->d1 = wts_mult[1];	tmp->d2 = wts_mult[0];	tmp->d3 = wts_mult[1];	++tmp;
		tmp->d0 = wts_mult[1];	tmp->d1 = wts_mult[1];	tmp->d2 = wts_mult[0];	tmp->d3 = wts_mult[1];	++tmp;
		tmp->d0 = wts_mult[0];	tmp->d1 = wts_mult[0];	tmp->d2 = wts_mult[1];	tmp->d3 = wts_mult[1];	++tmp;
		tmp->d0 = wts_mult[1];	tmp->d1 = wts_mult[0];	tmp->d2 = wts_mult[1];	tmp->d3 = wts_mult[1];	++tmp;
		tmp->d0 = wts_mult[0];	tmp->d1 = wts_mult[1];	tmp->d2 = wts_mult[1];	tmp->d3 = wts_mult[1];	++tmp;
		tmp->d0 = wts_mult[1];	tmp->d1 = wts_mult[1];	tmp->d2 = wts_mult[1];	tmp->d3 = wts_mult[1];	++tmp;
		/* inv_mult: */
		tmp->d0 = inv_mult[0];	tmp->d1 = inv_mult[0];	tmp->d2 = inv_mult[0];	tmp->d3 = inv_mult[0];	++tmp;
		tmp->d0 = inv_mult[1];	tmp->d1 = inv_mult[0];	tmp->d2 = inv_mult[0];	tmp->d3 = inv_mult[0];	++tmp;
		tmp->d0 = inv_mult[0];	tmp->d1 = inv_mult[1];	tmp->d2 = inv_mult[0];	tmp->d3 = inv_mult[0];	++tmp;
		tmp->d0 = inv_mult[1];	tmp->d1 = inv_mult[1];	tmp->d2 = inv_mult[0];	tmp->d3 = inv_mult[0];	++tmp;
		tmp->d0 = inv_mult[0];	tmp->d1 = inv_mult[0];	tmp->d2 = inv_mult[1];	tmp->d3 = inv_mult[0];	++tmp;
		tmp->d0 = inv_mult[1];	tmp->d1 = inv_mult[0];	tmp->d2 = inv_mult[1];	tmp->d3 = inv_mult[0];	++tmp;
		tmp->d0 = inv_mult[0];	tmp->d1 = inv_mult[1];	tmp->d2 = inv_mult[1];	tmp->d3 = inv_mult[0];	++tmp;
		tmp->d0 = inv_mult[1];	tmp->d1 = inv_mult[1];	tmp->d2 = inv_mult[1];	tmp->d3 = inv_mult[0];	++tmp;
		tmp->d0 = inv_mult[0];	tmp->d1 = inv_mult[0];	tmp->d2 = inv_mult[0];	tmp->d3 = inv_mult[1];	++tmp;
		tmp->d0 = inv_mult[1];	tmp->d1 = inv_mult[0];	tmp->d2 = inv_mult[0];	tmp->d3 = inv_mult[1];	++tmp;
		tmp->d0 = inv_mult[0];	tmp->d1 = inv_mult[1];	tmp->d2 = inv_mult[0];	tmp->d3 = inv_mult[1];	++tmp;
		tmp->d0 = inv_mult[1];	tmp->d1 = inv_mult[1];	tmp->d2 = inv_mult[0];	tmp->d3 = inv_mult[1];	++tmp;
		tmp->d0 = inv_mult[0];	tmp->d1 = inv_mult[0];	tmp->d2 = inv_mult[1];	tmp->d3 = inv_mult[1];	++tmp;
		tmp->d0 = inv_mult[1];	tmp->d1 = inv_mult[0];	tmp->d2 = inv_mult[1];	tmp->d3 = inv_mult[1];	++tmp;
		tmp->d0 = inv_mult[0];	tmp->d1 = inv_mult[1];	tmp->d2 = inv_mult[1];	tmp->d3 = inv_mult[1];	++tmp;
		tmp->d0 = inv_mult[1];	tmp->d1 = inv_mult[1];	tmp->d2 = inv_mult[1];	tmp->d3 = inv_mult[1];	++tmp;
		nbytes = 96 << L2_SZ_VD;

	  #elif defined(USE_SSE2)

		ctmp = (struct complex *)tmp;
		/* Forward-weight multipliers: */
		ctmp->re = 1.0;	ctmp->im = 1.0;	++ctmp;
		ctmp->re = .50;	ctmp->im = 1.0;	++ctmp;
		ctmp->re = 1.0;	ctmp->im = .50;	++ctmp;
		ctmp->re = .50;	ctmp->im = .50;	++ctmp;
		/* Inverse-weight multipliers (only needed for mersenne-mod): */
		ctmp->re = .50;	ctmp->im = .50;	++ctmp;
		ctmp->re = .25;	ctmp->im = .50;	++ctmp;
		ctmp->re = .50;	ctmp->im = .25;	++ctmp;
		ctmp->re = .25;	ctmp->im = .25;	++ctmp;
		/* Forward-base[] multipliers: */
		ctmp->re = base   [0];	ctmp->im = base   [0];	++ctmp;
		ctmp->re = base   [1];	ctmp->im = base   [0];	++ctmp;
		ctmp->re = base   [0];	ctmp->im = base   [1];	++ctmp;
		ctmp->re = base   [1];	ctmp->im = base   [1];	++ctmp;
		/* Inverse-base[] multipliers: */
		ctmp->re = baseinv[0];	ctmp->im = baseinv[0];	++ctmp;
		ctmp->re = baseinv[1];	ctmp->im = baseinv[0];	++ctmp;
		ctmp->re = baseinv[0];	ctmp->im = baseinv[1];	++ctmp;
		ctmp->re = baseinv[1];	ctmp->im = baseinv[1];	++ctmp;
		// In LOACC mode, put wts_mult and their inverses in the first 8 slots below in place of the 1/2-stuff:
		/* wts_mult:*/
		ctmp->re = wts_mult[0];	ctmp->im = wts_mult[0];	++ctmp;
		ctmp->re = wts_mult[1];	ctmp->im = wts_mult[0];	++ctmp;
		ctmp->re = wts_mult[0];	ctmp->im = wts_mult[1];	++ctmp;
		ctmp->re = wts_mult[1];	ctmp->im = wts_mult[1];	++ctmp;
		/* inv_mult:*/
		ctmp->re = inv_mult[0];	ctmp->im = inv_mult[0];	++ctmp;
		ctmp->re = inv_mult[1];	ctmp->im = inv_mult[0];	++ctmp;
		ctmp->re = inv_mult[0];	ctmp->im = inv_mult[1];	++ctmp;
		ctmp->re = inv_mult[1];	ctmp->im = inv_mult[1];	++ctmp;
		nbytes = 24 << L2_SZ_VD;

	  #endif

		// Propagate the above consts to the remaining threads:
		tmp = half_arr;
		tm2 = tmp + cslots_in_local_store;
		for(ithread = 1; ithread < CY_THREADS; ++ithread) {
			memcpy(tm2, tmp, nbytes);
			tmp = tm2;		tm2 += cslots_in_local_store;
		}

		/* Floating-point sign mask used for FABS on packed doubles: */
		sign_mask = sm_ptr;
		for(i = 0; i < RE_IM_STRIDE; ++i) {
			*(sign_mask+i) = (uint64)0x7FFFFFFFFFFFFFFFull;
		}

		// Set up the SIMD-tupled-32-bit-int SSE constants used by the carry macros:
		sse_bw  = sm_ptr + RE_IM_STRIDE;	// (#doubles in a SIMD complex) x 32-bits = RE_IM_STRIDE x 64-bits
		tmp64 = (uint64)bw;
		tmp64 = tmp64 + (tmp64 << 32);
		for(i = 0; i < RE_IM_STRIDE; ++i) {
			*(sse_bw+i) = tmp64;
		}

		sse_sw  = sse_bw + RE_IM_STRIDE;
		tmp64 = (uint64)sw;
		tmp64 = tmp64 + (tmp64 << 32);
		for(i = 0; i < RE_IM_STRIDE; ++i) {
			*(sse_sw+i) = tmp64;
		}

		sse_n   = sse_sw + RE_IM_STRIDE;
		tmp64 = (uint64)n;
		tmp64 = tmp64 + (tmp64 << 32);
		for(i = 0; i < RE_IM_STRIDE; ++i) {
			*(sse_n +i) = tmp64;
		}

		nbytes = 4 << L2_SZ_VD;

	  #ifdef USE_AVX512
	   #ifdef CARRY_16_WAY
		n_minus_sil   = (struct uint32x16*)sse_n + 1;
		n_minus_silp1 = (struct uint32x16*)sse_n + 2;
		sinwt         = (struct uint32x16*)sse_n + 3;
		sinwtm1       = (struct uint32x16*)sse_n + 4;
		nbytes += 256;
	   #else
		n_minus_sil   = (struct uint32x8 *)sse_n + 1;
		n_minus_silp1 = (struct uint32x8 *)sse_n + 2;
		sinwt         = (struct uint32x8 *)sse_n + 3;
		sinwtm1       = (struct uint32x8 *)sse_n + 4;
		nbytes += 128;
	   #endif
	  #elif defined(USE_AVX)
		n_minus_sil   = (struct uint32x4 *)sse_n + 1;
		n_minus_silp1 = (struct uint32x4 *)sse_n + 2;
		sinwt         = (struct uint32x4 *)sse_n + 3;
		sinwtm1       = (struct uint32x4 *)sse_n + 4;
		nbytes += 64;
	  #endif

		// Propagate the above consts to the remaining threads:
		tmp = (vec_dbl *)sm_ptr;
		tm2 = tmp + cslots_in_local_store;
		for(ithread = 1; ithread < CY_THREADS; ++ithread) {
			memcpy(tm2, tmp, nbytes);
			tmp = tm2;		tm2 += cslots_in_local_store;
		}

	// For large radices, array-access to bjmodn means only init base-ptr here:
	  #ifdef USE_AVX
		bjmodn = (int*)(sinwtm1 + RE_IM_STRIDE);
	  #else
		bjmodn = (int*)(sse_n   + RE_IM_STRIDE);
	  #endif

	#endif	// USE_SSE2

		/*   constant index offsets for array load/stores are here.	*/
		p1 = NDIVR;		p010 = NDIVR<<4;	p100 = NDIVR<<8;
		p2 = p1 + p1;	p020 = p010 + p010;	p110 = p100 + p010;
		p3 = p2 + p1;	p030 = p020 + p010;	p120 = p110 + p010;
		p4 = p3 + p1;	p040 = p030 + p010;	p130 = p120 + p010;
		p5 = p4 + p1;	p050 = p040 + p010;	p140 = p130 + p010;
		p6 = p5 + p1;	p060 = p050 + p010;	p150 = p140 + p010;
		p7 = p6 + p1;	p070 = p060 + p010;	p160 = p150 + p010;
		p8 = p7 + p1;	p080 = p070 + p010;	p170 = p160 + p010;
		p9 = p8 + p1;	p090 = p080 + p010;
		pa = p9 + p1;	p0a0 = p090 + p010;
		pb = pa + p1;	p0b0 = p0a0 + p010;
		pc = pb + p1;	p0c0 = p0b0 + p010;
		pd = pc + p1;	p0d0 = p0c0 + p010;
		pe = pd + p1;	p0e0 = p0d0 + p010;
		pf = pe + p1;	p0f0 = p0e0 + p010;

		p1 += ( (p1 >> DAT_BITS) << PAD_BITS );		p010 += ( (p010 >> DAT_BITS) << PAD_BITS );		p100 += ( (p100 >> DAT_BITS) << PAD_BITS );
		p2 += ( (p2 >> DAT_BITS) << PAD_BITS );		p020 += ( (p020 >> DAT_BITS) << PAD_BITS );		p110 += ( (p110 >> DAT_BITS) << PAD_BITS );
		p3 += ( (p3 >> DAT_BITS) << PAD_BITS );		p030 += ( (p030 >> DAT_BITS) << PAD_BITS );		p120 += ( (p120 >> DAT_BITS) << PAD_BITS );
		p4 += ( (p4 >> DAT_BITS) << PAD_BITS );		p040 += ( (p040 >> DAT_BITS) << PAD_BITS );		p130 += ( (p130 >> DAT_BITS) << PAD_BITS );
		p5 += ( (p5 >> DAT_BITS) << PAD_BITS );		p050 += ( (p050 >> DAT_BITS) << PAD_BITS );		p140 += ( (p140 >> DAT_BITS) << PAD_BITS );
		p6 += ( (p6 >> DAT_BITS) << PAD_BITS );		p060 += ( (p060 >> DAT_BITS) << PAD_BITS );		p150 += ( (p150 >> DAT_BITS) << PAD_BITS );
		p7 += ( (p7 >> DAT_BITS) << PAD_BITS );		p070 += ( (p070 >> DAT_BITS) << PAD_BITS );		p160 += ( (p160 >> DAT_BITS) << PAD_BITS );
		p8 += ( (p8 >> DAT_BITS) << PAD_BITS );		p080 += ( (p080 >> DAT_BITS) << PAD_BITS );		p170 += ( (p170 >> DAT_BITS) << PAD_BITS );
		p9 += ( (p9 >> DAT_BITS) << PAD_BITS );		p090 += ( (p090 >> DAT_BITS) << PAD_BITS );
		pa += ( (pa >> DAT_BITS) << PAD_BITS );		p0a0 += ( (p0a0 >> DAT_BITS) << PAD_BITS );
		pb += ( (pb >> DAT_BITS) << PAD_BITS );		p0b0 += ( (p0b0 >> DAT_BITS) << PAD_BITS );
		pc += ( (pc >> DAT_BITS) << PAD_BITS );		p0c0 += ( (p0c0 >> DAT_BITS) << PAD_BITS );
		pd += ( (pd >> DAT_BITS) << PAD_BITS );		p0d0 += ( (p0d0 >> DAT_BITS) << PAD_BITS );
		pe += ( (pe >> DAT_BITS) << PAD_BITS );		p0e0 += ( (p0e0 >> DAT_BITS) << PAD_BITS );
		pf += ( (pf >> DAT_BITS) << PAD_BITS );		p0f0 += ( (p0f0 >> DAT_BITS) << PAD_BITS );
	#ifndef USE_SSE2
		p0123[0] = 0; p0123[1] = p1; p0123[2] = p2; p0123[3] = p3;
	#endif
		poff[     0] =    0; poff[     1] =      p4; poff[     2] =      p8; poff[     3] =      pc;
		poff[0x04+0] = p010; poff[0x04+1] = p010+p4; poff[0x04+2] = p010+p8; poff[0x04+3] = p010+pc;
		poff[0x08+0] = p020; poff[0x08+1] = p020+p4; poff[0x08+2] = p020+p8; poff[0x08+3] = p020+pc;
		poff[0x0c+0] = p030; poff[0x0c+1] = p030+p4; poff[0x0c+2] = p030+p8; poff[0x0c+3] = p030+pc;
		poff[0x10+0] = p040; poff[0x10+1] = p040+p4; poff[0x10+2] = p040+p8; poff[0x10+3] = p040+pc;
		poff[0x14+0] = p050; poff[0x14+1] = p050+p4; poff[0x14+2] = p050+p8; poff[0x14+3] = p050+pc;
		poff[0x18+0] = p060; poff[0x18+1] = p060+p4; poff[0x18+2] = p060+p8; poff[0x18+3] = p060+pc;
		poff[0x1c+0] = p070; poff[0x1c+1] = p070+p4; poff[0x1c+2] = p070+p8; poff[0x1c+3] = p070+pc;
		poff[0x20+0] = p080; poff[0x20+1] = p080+p4; poff[0x20+2] = p080+p8; poff[0x20+3] = p080+pc;
		poff[0x24+0] = p090; poff[0x24+1] = p090+p4; poff[0x24+2] = p090+p8; poff[0x24+3] = p090+pc;
		poff[0x28+0] = p0a0; poff[0x28+1] = p0a0+p4; poff[0x28+2] = p0a0+p8; poff[0x28+3] = p0a0+pc;
		poff[0x2c+0] = p0b0; poff[0x2c+1] = p0b0+p4; poff[0x2c+2] = p0b0+p8; poff[0x2c+3] = p0b0+pc;
		poff[0x30+0] = p0c0; poff[0x30+1] = p0c0+p4; poff[0x30+2] = p0c0+p8; poff[0x30+3] = p0c0+pc;
		poff[0x34+0] = p0d0; poff[0x34+1] = p0d0+p4; poff[0x34+2] = p0d0+p8; poff[0x34+3] = p0d0+pc;
		poff[0x38+0] = p0e0; poff[0x38+1] = p0e0+p4; poff[0x38+2] = p0e0+p8; poff[0x38+3] = p0e0+pc;
		poff[0x3c+0] = p0f0; poff[0x3c+1] = p0f0+p4; poff[0x3c+2] = p0f0+p8; poff[0x3c+3] = p0f0+pc;
		for(l = 0; l < 32; l++) {
			poff[ 64+l] = poff[l] + p100;
		}

	#ifndef MULTITHREAD
	/*** DIF: ***/
	// Set array offsets for radix-128 inputs [these are same for all 3 such DFTs].
		for(j = 0; j < 0x80; ++j) { dif_i_offsets[j] = j<<1; }
	// For the radix-128 DIF outputs we need the following offsets:
		const int dif_perm0[16] = {0,p1,p2,p3,p5,p4,p7,p6,pa,pb,p9,p8,pf,pe,pc,pd};
		const int dif_perm1[16] = {p5,p4,p7,p6,p2,p3,p1,0,pf,pe,pc,pd,p9,p8,pb,pa};
		const int dif_perm2[16] = {pa,pb,p9,p8,pf,pe,pc,pd,p5,p4,p7,p6,p2,p3,p1,0};
		const int dif_perm3[16] = {pf,pe,pc,pd,p9,p8,pb,pa,p2,p3,p1,0,p7,p6,p4,p5};
		for(j = 0; j < 0x10; ++j) {
			dif_o_offsets[j] = dif_perm0[j];
			dif_o_offsets[j+0x010] = p010 + dif_perm1[j];
			dif_o_offsets[j+0x020] = p020 + dif_perm2[j];
			dif_o_offsets[j+0x030] = p030 + dif_perm3[j];
			dif_o_offsets[j+0x040] = p050 + dif_perm1[j];
			dif_o_offsets[j+0x050] = p040 + dif_perm2[j];
			dif_o_offsets[j+0x060] = p070 + dif_perm3[j];
			dif_o_offsets[j+0x070] = p060 + dif_perm1[j];
		// dif_o_offsets are w.r.to each respective block of 128 complex data, i.e. p-indexing is (mod p80):
			dif_o_offsets[j+0x080] = p020 + dif_perm2[j];
			dif_o_offsets[j+0x090] = p030 + dif_perm3[j];
			dif_o_offsets[j+0x0a0] = p010 + dif_perm1[j];
			dif_o_offsets[j+0x0b0] =    0 + dif_perm2[j];
			dif_o_offsets[j+0x0c0] = p070 + dif_perm3[j];
			dif_o_offsets[j+0x0d0] = p060 + dif_perm1[j];
			dif_o_offsets[j+0x0e0] = p040 + dif_perm2[j];
			dif_o_offsets[j+0x0f0] = p050 + dif_perm3[j];
		// dif_o_offsets are w.r.to each respective block of 128 complex data, i.e. p-indexing is (mod p80):
			dif_o_offsets[j+0x100] = p050 + dif_perm1[j];
			dif_o_offsets[j+0x110] = p040 + dif_perm2[j];
			dif_o_offsets[j+0x120] = p070 + dif_perm3[j];
			dif_o_offsets[j+0x130] = p060 + dif_perm1[j];
			dif_o_offsets[j+0x140] = p020 + dif_perm2[j];
			dif_o_offsets[j+0x150] = p030 + dif_perm3[j];
			dif_o_offsets[j+0x160] = p010 + dif_perm1[j];
			dif_o_offsets[j+0x170] =    0 + dif_perm2[j];
		}
	/*** DIT: ***/
	// Set array offsets for radix-128 outputs [these are same for all 3 such DFTs].
		for(j = 0; j < 0x80; ++j) { dit_o_offsets[j] = j<<1; }
	// For the radix-128 DIT inputs we need the following offsets:
		const int dit_perm0[16] = {0,p1,p3,p2,p7,p6,p5,p4,pf,pe,pd,pc,pb,pa,p9,p8};
		const int dit_perm1[16] = {pf,pe,pd,pc,pb,pa,p9,p8,p7,p6,p5,p4,p3,p2,p1,0};
		const int dit_perm2[16] = {p5,p4,p6,p7,p1,0,p2,p3,p9,p8,pa,pb,pe,pf,pc,pd};
		const int dit_perm3[16] = {p9,p8,pa,pb,pe,pf,pc,pd,p1,0,p2,p3,p6,p7,p4,p5};
		const int dit_perm4[16] = {pa,pb,p8,p9,pc,pd,pf,pe,p2,p3,0,p1,p4,p5,p7,p6};
		const int dit_perm5[16] = {p2,p3,0,p1,p4,p5,p7,p6,pc,pd,pf,pe,p8,p9,pb,pa};
		for(j = 0; j < 0x10; ++j) {
			dit_i_offsets[j] = dit_perm0[j];
			dit_i_offsets[j+0x010] = p010 + dit_perm1[j];
			dit_i_offsets[j+0x020] = p030 + dit_perm1[j];
			dit_i_offsets[j+0x030] = p020 + dit_perm1[j];
			dit_i_offsets[j+0x040] = p070 + dit_perm1[j];
			dit_i_offsets[j+0x050] = p060 + dit_perm1[j];
			dit_i_offsets[j+0x060] = p050 + dit_perm1[j];
			dit_i_offsets[j+0x070] = p040 + dit_perm1[j];
		// dit_i_offsets are w.r.to each respective block of 128 complex data, i.e. p-indexing is (mod p80):
			dit_i_offsets[j+0x080] = p050 + dit_perm2[j];
			dit_i_offsets[j+0x090] = p040 + dit_perm2[j];
			dit_i_offsets[j+0x0a0] = p060 + dit_perm2[j];
			dit_i_offsets[j+0x0b0] = p070 + dit_perm3[j];
			dit_i_offsets[j+0x0c0] = p010 + dit_perm2[j];
			dit_i_offsets[j+0x0d0] =    0 + dit_perm2[j];
			dit_i_offsets[j+0x0e0] = p020 + dit_perm2[j];
			dit_i_offsets[j+0x0f0] = p030 + dit_perm3[j];
		// dit_i_offsets are w.r.to each respective block of 128 complex data, i.e. p-indexing is (mod p80):
			dit_i_offsets[j+0x100] = p020 + dit_perm4[j];
			dit_i_offsets[j+0x110] = p030 + dit_perm5[j];
			dit_i_offsets[j+0x120] =    0 + dit_perm4[j];
			dit_i_offsets[j+0x130] = p010 + dit_perm5[j];
			dit_i_offsets[j+0x140] = p040 + dit_perm4[j];
			dit_i_offsets[j+0x150] = p050 + dit_perm5[j];
			dit_i_offsets[j+0x160] = p070 + dit_perm5[j];
			dit_i_offsets[j+0x170] = p060 + dit_perm5[j];
		}
	// DIF Index-high-bits triplets needed for compact-obj-code scheme:
	  #ifdef USE_SSE2
		k = 0;	// In SIMD mode these are 0xtr-offsets w.r.to a local store:
		dif_triplets[k] = 0x170; dif_triplets[k+1] = 0x0f0; dif_triplets[k+2] = 0x070; k += 3;
		dif_triplets[k] = 0x160; dif_triplets[k+1] = 0x0e0; dif_triplets[k+2] = 0x060; k += 3;
		dif_triplets[k] = 0x150; dif_triplets[k+1] = 0x0d0; dif_triplets[k+2] = 0x050; k += 3;
		dif_triplets[k] = 0x140; dif_triplets[k+1] = 0x0c0; dif_triplets[k+2] = 0x040; k += 3;
		dif_triplets[k] = 0x130; dif_triplets[k+1] = 0x0b0; dif_triplets[k+2] = 0x030; k += 3;
		dif_triplets[k] = 0x120; dif_triplets[k+1] = 0x0a0; dif_triplets[k+2] = 0x020; k += 3;
		dif_triplets[k] = 0x110; dif_triplets[k+1] = 0x090; dif_triplets[k+2] = 0x010; k += 3;
		dif_triplets[k] = 0x100; dif_triplets[k+1] = 0x080; dif_triplets[k+2] =     0; k += 3;
		dif_triplets[k] = 0x0f0; dif_triplets[k+1] = 0x070; dif_triplets[k+2] = 0x170; k += 3;
		dif_triplets[k] = 0x0e0; dif_triplets[k+1] = 0x060; dif_triplets[k+2] = 0x160; k += 3;
		dif_triplets[k] = 0x0d0; dif_triplets[k+1] = 0x050; dif_triplets[k+2] = 0x150; k += 3;
		dif_triplets[k] = 0x0c0; dif_triplets[k+1] = 0x040; dif_triplets[k+2] = 0x140; k += 3;
		dif_triplets[k] = 0x0b0; dif_triplets[k+1] = 0x030; dif_triplets[k+2] = 0x130; k += 3;
		dif_triplets[k] = 0x0a0; dif_triplets[k+1] = 0x020; dif_triplets[k+2] = 0x120; k += 3;
		dif_triplets[k] = 0x090; dif_triplets[k+1] = 0x010; dif_triplets[k+2] = 0x110; k += 3;
		dif_triplets[k] = 0x080; dif_triplets[k+1] =     0; dif_triplets[k+2] = 0x100; k += 3;
		dif_triplets[k] = 0x070; dif_triplets[k+1] = 0x170; dif_triplets[k+2] = 0x0f0; k += 3;
		dif_triplets[k] = 0x060; dif_triplets[k+1] = 0x160; dif_triplets[k+2] = 0x0e0; k += 3;
		dif_triplets[k] = 0x050; dif_triplets[k+1] = 0x150; dif_triplets[k+2] = 0x0d0; k += 3;
		dif_triplets[k] = 0x040; dif_triplets[k+1] = 0x140; dif_triplets[k+2] = 0x0c0; k += 3;
		dif_triplets[k] = 0x030; dif_triplets[k+1] = 0x130; dif_triplets[k+2] = 0x0b0; k += 3;
		dif_triplets[k] = 0x020; dif_triplets[k+1] = 0x120; dif_triplets[k+2] = 0x0a0; k += 3;
		dif_triplets[k] = 0x010; dif_triplets[k+1] = 0x110; dif_triplets[k+2] = 0x090; k += 3;
		dif_triplets[k] =     0; dif_triplets[k+1] = 0x100; dif_triplets[k+2] = 0x080;
		// IN SIMD mode need to double all the above to turn from vec_dbl to vec_cmplx ptr offsets:
		for(l = 0; l < 72; l++) {
			dif_triplets[l] <<= 1;
		}
	  #else
	// Index-high-bits triplets needed for compact-obj-code scheme:
		k = 0;
		dif_triplets[k] = p170; dif_triplets[k+1] = p0f0; dif_triplets[k+2] = p070; k += 3;
		dif_triplets[k] = p160; dif_triplets[k+1] = p0e0; dif_triplets[k+2] = p060; k += 3;
		dif_triplets[k] = p150; dif_triplets[k+1] = p0d0; dif_triplets[k+2] = p050; k += 3;
		dif_triplets[k] = p140; dif_triplets[k+1] = p0c0; dif_triplets[k+2] = p040; k += 3;
		dif_triplets[k] = p130; dif_triplets[k+1] = p0b0; dif_triplets[k+2] = p030; k += 3;
		dif_triplets[k] = p120; dif_triplets[k+1] = p0a0; dif_triplets[k+2] = p020; k += 3;
		dif_triplets[k] = p110; dif_triplets[k+1] = p090; dif_triplets[k+2] = p010; k += 3;
		dif_triplets[k] = p100; dif_triplets[k+1] = p080; dif_triplets[k+2] =    0; k += 3;
		dif_triplets[k] = p0f0; dif_triplets[k+1] = p070; dif_triplets[k+2] = p170; k += 3;
		dif_triplets[k] = p0e0; dif_triplets[k+1] = p060; dif_triplets[k+2] = p160; k += 3;
		dif_triplets[k] = p0d0; dif_triplets[k+1] = p050; dif_triplets[k+2] = p150; k += 3;
		dif_triplets[k] = p0c0; dif_triplets[k+1] = p040; dif_triplets[k+2] = p140; k += 3;
		dif_triplets[k] = p0b0; dif_triplets[k+1] = p030; dif_triplets[k+2] = p130; k += 3;
		dif_triplets[k] = p0a0; dif_triplets[k+1] = p020; dif_triplets[k+2] = p120; k += 3;
		dif_triplets[k] = p090; dif_triplets[k+1] = p010; dif_triplets[k+2] = p110; k += 3;
		dif_triplets[k] = p080; dif_triplets[k+1] =    0; dif_triplets[k+2] = p100; k += 3;
		dif_triplets[k] = p070; dif_triplets[k+1] = p170; dif_triplets[k+2] = p0f0; k += 3;
		dif_triplets[k] = p060; dif_triplets[k+1] = p160; dif_triplets[k+2] = p0e0; k += 3;
		dif_triplets[k] = p050; dif_triplets[k+1] = p150; dif_triplets[k+2] = p0d0; k += 3;
		dif_triplets[k] = p040; dif_triplets[k+1] = p140; dif_triplets[k+2] = p0c0; k += 3;
		dif_triplets[k] = p030; dif_triplets[k+1] = p130; dif_triplets[k+2] = p0b0; k += 3;
		dif_triplets[k] = p020; dif_triplets[k+1] = p120; dif_triplets[k+2] = p0a0; k += 3;
		dif_triplets[k] = p010; dif_triplets[k+1] = p110; dif_triplets[k+2] = p090; k += 3;
		dif_triplets[k] =    0; dif_triplets[k+1] = p100; dif_triplets[k+2] = p080;
	  #endif
	// DIT Index-high-bits triplets needed for compact-obj-code scheme:
	  #ifdef USE_SSE2
		k = 0;	// In SIMD mode these are ptr-offsets w.r.to a local store:
		dit_triplets[k] = 0x070; dit_triplets[k+1] = 0x0f0; dit_triplets[k+2] = 0x170; k += 3;
		dit_triplets[k] = 0x160; dit_triplets[k+1] = 0x060; dit_triplets[k+2] = 0x0e0; k += 3;
		dit_triplets[k] = 0x0d0; dit_triplets[k+1] = 0x150; dit_triplets[k+2] = 0x050; k += 3;
		dit_triplets[k] = 0x040; dit_triplets[k+1] = 0x0c0; dit_triplets[k+2] = 0x140; k += 3;
		dit_triplets[k] = 0x130; dit_triplets[k+1] = 0x030; dit_triplets[k+2] = 0x0b0; k += 3;
		dit_triplets[k] = 0x0a0; dit_triplets[k+1] = 0x120; dit_triplets[k+2] = 0x020; k += 3;
		dit_triplets[k] = 0x010; dit_triplets[k+1] = 0x090; dit_triplets[k+2] = 0x110; k += 3;
		dit_triplets[k] = 0x100; dit_triplets[k+1] =     0; dit_triplets[k+2] = 0x080;
		// IN SIMD mode need to double all the above to turn from vec_dbl to vec_cmplx ptr offsets:
		for(l = 0; l < 24; l++) {
			dit_triplets[l] <<= 1;
		}
	  #else
		// Cf. comments in radix384_dit_pass1 re. this needed reordering:
		k = 0;
		dit_triplets[k] = p070; dit_triplets[k+1] = p0f0; dit_triplets[k+2] = p170; k += 3;
		dit_triplets[k] = p160; dit_triplets[k+1] = p060; dit_triplets[k+2] = p0e0; k += 3;
		dit_triplets[k] = p0d0; dit_triplets[k+1] = p150; dit_triplets[k+2] = p050; k += 3;
		dit_triplets[k] = p040; dit_triplets[k+1] = p0c0; dit_triplets[k+2] = p140; k += 3;
		dit_triplets[k] = p130; dit_triplets[k+1] = p030; dit_triplets[k+2] = p0b0; k += 3;
		dit_triplets[k] = p0a0; dit_triplets[k+1] = p120; dit_triplets[k+2] = p020; k += 3;
		dit_triplets[k] = p010; dit_triplets[k+1] = p090; dit_triplets[k+2] = p110; k += 3;
		dit_triplets[k] = p100; dit_triplets[k+1] =    0; dit_triplets[k+2] = p080;
	  #endif
	#endif	// MULTITHREAD ?

		if(_cy[0])	/* If it's a new exponent of a range test, need to deallocate these. */
		{
			free((void *)_i     ); _i      = 0x0;
			for(i = 0; i < RADIX; i++) {
				free((void *)_bjmodn[i]); _bjmodn[i] = 0x0;
				free((void *)    _cy[i]);     _cy[i] = 0x0;
			}
			free((void *)_jstart ); _jstart  = 0x0;
			free((void *)_jhi    ); _jhi     = 0x0;
			free((void *)_col   ); _col    = 0x0;
			free((void *)_co2   ); _co2    = 0x0;
			free((void *)_co3   ); _co3    = 0x0;
			free((void *)_bjmodnini); _bjmodnini = 0x0;
		}

		ptr_prod = (uint32)0;	/* Store bitmask for allocatable-array ptrs here, check vs 0 after all alloc calls finish */
		j = CY_THREADS*sizeof(int);
		_i       	= (int *)malloc(j);	ptr_prod += (uint32)(_i== 0x0);
		for(i = 0; i < RADIX; i++) {
			_bjmodn[i]	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn[i]== 0x0);
		}
		_jstart  	= (int *)malloc(j);	ptr_prod += (uint32)(_jstart  == 0x0);
		_jhi     	= (int *)malloc(j);	ptr_prod += (uint32)(_jhi     == 0x0);
		_col     	= (int *)malloc(j);	ptr_prod += (uint32)(_col     == 0x0);
		_co2     	= (int *)malloc(j);	ptr_prod += (uint32)(_co2     == 0x0);
		_co3     	= (int *)malloc(j);	ptr_prod += (uint32)(_co3     == 0x0);

		j = CY_THREADS*sizeof(double);
		for(i = 0; i < RADIX; i++) {
			_cy[i]	= (double *)malloc(j);	ptr_prod += (uint32)(_cy[i]== 0x0);
		}

		ASSERT(ptr_prod == 0, "ERROR: unable to allocate one or more auxiliary arrays.");

		/* Create (THREADS + 1) copies of _bjmodnini and use the extra (uppermost) one to store the "master" increment,
		i.e. the one that n2/RADIX-separated FFT outputs need:
		*/
		_bjmodnini = (int *)malloc((CY_THREADS + 1)*sizeof(int));	if(!_bjmodnini){ sprintf(cbuf,"ERROR: unable to allocate array _bjmodnini in %s.\n",func); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		_bjmodnini[0] = 0;
		_bjmodnini[1] = 0;

		jhi = NDIVR/CY_THREADS;

		for(j=0; j < jhi; j++)
		{
			_bjmodnini[1] -= sw; _bjmodnini[1] = _bjmodnini[1] + ( (-(int)((uint32)_bjmodnini[1] >> 31)) & n);
		}

		if(CY_THREADS > 1)
		{
			for(ithread = 2; ithread <= CY_THREADS; ithread++)
			{
				_bjmodnini[ithread] = _bjmodnini[ithread-1] + _bjmodnini[1] - n; _bjmodnini[ithread] = _bjmodnini[ithread] + ( (-(int)((uint32)_bjmodnini[ithread] >> 31)) & n);
			}
		}
		/* Check upper element against scalar value, as precomputed in single-thread mode: */
		bjmodnini=0;
		for(j=0; j < jhi*CY_THREADS; j++)
		{
			bjmodnini -= sw; bjmodnini = bjmodnini + ( (-(int)((uint32)bjmodnini >> 31)) & n);
		}
		ASSERT(_bjmodnini[CY_THREADS] == bjmodnini,"_bjmodnini[CY_THREADS] != bjmodnini");

	#ifdef USE_PTHREAD
		/* Populate the elements of the thread-specific data structs which don't change after init: */
		for(ithread = 0; ithread < CY_THREADS; ithread++)
		{
			tdat[ithread].bjmodnini = _bjmodnini[CY_THREADS];
			tdat[ithread].bjmodn0 = _bjmodnini[ithread];
		#ifdef USE_SSE2
			tdat[ithread].r000 = __r0 + ithread*cslots_in_local_store;
			tdat[ithread].half_arr = (vec_dbl *)((intptr_t)tdat[ithread].r000 + ((intptr_t)half_arr - (intptr_t)r000));
		#else	// In scalar mode use these 2 ptrs to pass the base & baseinv arrays:
			tdat[ithread].r000     = (double *)base;
			tdat[ithread].half_arr = (double *)baseinv;
		#endif	// USE_SSE2
		}
	#endif

		first_entry=FALSE;
	}	/* endif(first_entry) */

	// Jun 2018: If LL test and shift applied, compute target index for data-processing loop.
	// Note that only 1 thread of the carry-processing set will hit the target, but all need the same logic to check for a hit:
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE && TEST_TYPE == TEST_TYPE_PRIMALITY) {
		if(RES_SHIFT) {
			itmp64 = shift_word(a, n, p, RES_SHIFT, 0.0);	// Note return value (specifically high 7 bytes thereof) is an unpadded index
			target_idx = (int)(itmp64 >>  8);	// This still needs to be (mod NDIVR)'ed, but first use unmodded form to compute needed DWT weights
			// Compute wt = 2^(target_idx*sw % n)/n and its reciprocal:
			uint32 sw_idx_modn = ((uint64)target_idx*sw) % n;	// N is 32-bit, so only use 64-bit to hold intermediate product
			double target_wtfwd = pow(2.0, sw_idx_modn*0.5*n2inv);	// 0.5*n2inv = 0.5/(n/2) = 1.0/n
			target_set = target_idx*ndivr_inv;	// Which of the [RADIX] independent sub-carry-chains contains the target index?
			target_idx -= target_set*NDIVR;		// Fast computation of target_idx = (target_idx % NDIVR)
			// Now compute the doubles-pointer offset of the target double w.r.to the SIMD s1p00-... data layout:
			tidx_mod_stride = target_idx & (stride-1);	// Stride a power of 2, so can use AND-minus-1 for mod
			target_idx -= tidx_mod_stride;
		//	printf("Iter %d: cy_shift = %d, target_idx,tidx_mod_stride,target_set = %d,%d,%d\n",iter,(itmp64 & 255),target_idx,tidx_mod_stride,target_set);
		#ifdef USE_AVX512
			tidx_mod_stride = br16[tidx_mod_stride];
		#elif defined(USE_AVX)
			tidx_mod_stride = br8[tidx_mod_stride];
		#elif defined(USE_SSE2)
			tidx_mod_stride = br4[tidx_mod_stride];
		#endif
			target_set = (target_set<<(L2_SZ_VD-2)) + tidx_mod_stride;
			target_cy  = target_wtfwd * (-(int)(2u << (itmp64 & 255)));
		} else {
			target_idx = target_set = 0;
			target_cy = -2.0;
		}
	}

/*...The radix-384 final DIT pass is here.	*/

	/* init carries	*/
	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		for(i = 0; i < RADIX; i++) {
			_cy[i][ithread] = 0;
		}
	}
  #if 0	//ndef USE_SSE2	*** v20: Non-SIMD builds now also support shifted-residue
	/* If an LL test, init the subtract-2: */
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE && TEST_TYPE == TEST_TYPE_PRIMALITY)
	{
		_cy[0][0] = -2;
	}
  #endif
	*fracmax=0;	/* init max. fractional error	*/
	full_pass = 1;	/* set = 1 for normal carry pass, = 0 for wrapper pass	*/
	scale = n2inv;	// init inverse-weight scale factor = 2/n for normal carry pass, 1 for wrapper pass

for(outer=0; outer <= 1; outer++)
{
	_i[0] = 1;		/* Pointer to the BASE and BASEINV arrays. If n does not divide p, lowest-order digit is always a bigword (_i[0] = 1).	*/

	if(CY_THREADS > 1)
	{
		for(ithread = 1; ithread < CY_THREADS; ithread++)
		{
			_i[ithread] = ((uint32)(sw - _bjmodnini[ithread]) >> 31);
		}
	}

	/*
	Moved this inside the outer-loop, so on cleanup pass can use it to reset _col,_co2,_co3 starting values,
	then simply overwrite it with 1 prior to starting the k-loop.
	*/
	khi = n_div_nwt/CY_THREADS;
	j = _bjmodnini[CY_THREADS];
	// Include 0-thread here ... bjmodn terms all 0 for that, but need jhi computed for all threads:
	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		_bjmodn[0][ithread] = _bjmodnini[ithread];
		for(i = 1; i < RADIX; i++) {
			MOD_ADD32(_bjmodn[i-1][ithread], j, n, _bjmodn[i][ithread]);
		}
		_jstart[ithread] = ithread*NDIVR/CY_THREADS;
		if(!full_pass)
			_jhi[ithread] = _jstart[ithread] + jhi_wrap;		/* Cleanup loop assumes carryins propagate at most 4 words up. */
		else
			_jhi[ithread] = _jstart[ithread] + nwt-1;

		_col[ithread] = ithread*(khi*RADIX);			/* col gets incremented by RADIX_VEC[0] on every pass through the k-loop */
		_co2[ithread] = (n>>nwt_bits)-1+RADIX - _col[ithread];	/* co2 gets decremented by RADIX_VEC[0] on every pass through the k-loop */
		_co3[ithread] = _co2[ithread]-RADIX;			/* At the start of each new j-loop, co3=co2-RADIX_VEC[0]	*/
	}

#ifdef USE_PTHREAD
	for(ithread = 0; ithread < CY_THREADS; ++ithread) { tdat[ithread].iter = iter; }
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
		// Carry-injection location for the shifted-residue -2 addend is only needed for full pass:
		if(full_pass) {
			tdat[0].target_idx = target_idx;
			tdat[0].target_set = target_set;
			tdat[0].target_cy  = target_cy;
		} else {
			tdat[0].target_idx = -1;
			tdat[0].target_set = 0;
			tdat[0].target_cy  = 0;
		}
		// Copy to the remaining threads:
		for(ithread = 1; ithread < CY_THREADS; ++ithread) {
			tdat[ithread].target_idx = tdat[0].target_idx;
			tdat[ithread].target_set = tdat[0].target_set;
			tdat[ithread].target_cy  = tdat[0].target_cy;
		}
	}
#endif
#ifdef USE_SSE2

	tmp = max_err;	VEC_DBL_INIT(tmp, 0.0);
	tm2 = tmp + cslots_in_local_store;
	for(ithread = 1; ithread < CY_THREADS; ++ithread) {
		memcpy(tm2, tmp, SZ_VD);
		tmp = tm2;		tm2 += cslots_in_local_store;
	}

#endif	// USE_PTHREAD

	/* Move this cleanup-pass-specific khi setting here, since need regular-pass khi value for above inits: */
	if(!full_pass)
	{
		khi = 1;
	}

#ifdef USE_PTHREAD
	/* Populate the thread-specific data structs - use the invariant terms as memchecks: */
	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		tdat[ithread].iter = iter;
	// int data:
		ASSERT(tdat[ithread].tid == ithread, "thread-local memcheck fail!");
		ASSERT(tdat[ithread].ndivr == NDIVR, "thread-local memcheck fail!");

		tdat[ithread].khi    = khi;
		tdat[ithread].i      = _i[ithread];	/* Pointer to the BASE and BASEINV arrays.	*/
		tdat[ithread].jstart = _jstart[ithread];
		tdat[ithread].jhi    = _jhi[ithread];

		tdat[ithread].col = _col[ithread];
		tdat[ithread].co2 = _co2[ithread];
		tdat[ithread].co3 = _co3[ithread];
		ASSERT(tdat[ithread].sw  == sw, "thread-local memcheck fail!");
		ASSERT(tdat[ithread].nwt == nwt, "thread-local memcheck fail!");

	// double data:
		tdat[ithread].maxerr = 0.0;
		tdat[ithread].scale = scale;
		tdat[ithread].prp_mult = prp_mult;

	// pointer data:
		tdat[ithread].arrdat = a;			/* Main data array */
		ASSERT(tdat[ithread].wt0 == wt0, "thread-local memcheck fail!");
		ASSERT(tdat[ithread].wt1 == wt1, "thread-local memcheck fail!");
		ASSERT(tdat[ithread].si  == si, "thread-local memcheck fail!");
	#ifdef USE_SSE2
		ASSERT(tdat[ithread].r000 == __r0 + ithread*cslots_in_local_store, "thread-local memcheck fail!");
		tmp = tdat[ithread].half_arr;
	  #ifdef USE_AVX512	// In AVX-512 mode, use VRNDSCALEPD for rounding and hijack this vector-data slot for the 4 base/baseinv-consts
		ASSERT(((tmp-1)->d0 == base[0] && (tmp-1)->d1 == baseinv[1] && (tmp-1)->d2 == wts_mult[1] && (tmp-1)->d3 == inv_mult[0]), "thread-local memcheck failed!");
	  #else
		ASSERT(((tmp-1)->d0 == crnd && (tmp-1)->d1 == crnd), "thread-local memcheck failed!");
	  #endif
	  #ifdef USE_AVX512
			/* No-Op */
	  #elif defined(USE_AVX)
		// Grab some elt of base-data [offset by, say, +32] and mpy by its inverse [+16 further]
		dtmp = (tmp+40)->d0 * (tmp+56)->d0;	ASSERT(fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
		dtmp = (tmp+40)->d1 * (tmp+56)->d1;	ASSERT(fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
	  #else	// SSE2:
		dtmp = (tmp+10)->d0 * (tmp+14)->d0;	ASSERT(fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
		dtmp = (tmp+10)->d1 * (tmp+14)->d1;	ASSERT(fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
	  #endif
	#endif
		/* init carries: */
		for(i = 0; i < RADIX; i++) {
			tdat[ithread].cy[i] = _cy[i][ithread];
		}
	}
#endif

/******************* AVX debug stuff: *******************/
#if 0
	int ipad;
	ASSERT(p1 >= 16, "Smallest array-stride must be large enough to hold an AVX-512 vec_cmplx!");
	// Use RNG to populate data array:
	rng_isaac_init(TRUE);
	double dtmp = 1024.0*1024.0*1024.0*1024.0;
	for(i = 0; i < n; i += 16) {
		ipad = i + ( (i >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
		// All the inits are w.r.to an un-SIMD-rearranged ...,re,im,re,im,... pattern:
	#ifdef USE_AVX512
		a[ipad+br16[ 0]] = dtmp*rng_isaac_rand_double_norm_pm1();	// re0
		a[ipad+br16[ 1]] = dtmp*rng_isaac_rand_double_norm_pm1();	// im0
		a[ipad+br16[ 2]] = dtmp*rng_isaac_rand_double_norm_pm1();	// re1
		a[ipad+br16[ 3]] = dtmp*rng_isaac_rand_double_norm_pm1();	// im1
		a[ipad+br16[ 4]] = dtmp*rng_isaac_rand_double_norm_pm1();	// re2
		a[ipad+br16[ 5]] = dtmp*rng_isaac_rand_double_norm_pm1();	// im2
		a[ipad+br16[ 6]] = dtmp*rng_isaac_rand_double_norm_pm1();	// re3
		a[ipad+br16[ 7]] = dtmp*rng_isaac_rand_double_norm_pm1();	// im3
		a[ipad+br16[ 8]] = dtmp*rng_isaac_rand_double_norm_pm1();	// re4
		a[ipad+br16[ 9]] = dtmp*rng_isaac_rand_double_norm_pm1();	// im4
		a[ipad+br16[10]] = dtmp*rng_isaac_rand_double_norm_pm1();	// re5
		a[ipad+br16[11]] = dtmp*rng_isaac_rand_double_norm_pm1();	// im5
		a[ipad+br16[12]] = dtmp*rng_isaac_rand_double_norm_pm1();	// re6
		a[ipad+br16[13]] = dtmp*rng_isaac_rand_double_norm_pm1();	// im6
		a[ipad+br16[14]] = dtmp*rng_isaac_rand_double_norm_pm1();	// re7
		a[ipad+br16[15]] = dtmp*rng_isaac_rand_double_norm_pm1();	// im7
	#elif defined(USE_AVX)
		a[ipad+br8[0]  ] = dtmp*rng_isaac_rand_double_norm_pm1();	// re0
		a[ipad+br8[1]  ] = dtmp*rng_isaac_rand_double_norm_pm1();	// im0
		a[ipad+br8[2]  ] = dtmp*rng_isaac_rand_double_norm_pm1();	// re1
		a[ipad+br8[3]  ] = dtmp*rng_isaac_rand_double_norm_pm1();	// im1
		a[ipad+br8[4]  ] = dtmp*rng_isaac_rand_double_norm_pm1();	// re2
		a[ipad+br8[5]  ] = dtmp*rng_isaac_rand_double_norm_pm1();	// im2
		a[ipad+br8[6]  ] = dtmp*rng_isaac_rand_double_norm_pm1();	// re3
		a[ipad+br8[7]  ] = dtmp*rng_isaac_rand_double_norm_pm1();	// im3
		a[ipad+br8[0]+8] = dtmp*rng_isaac_rand_double_norm_pm1();	// re4
		a[ipad+br8[1]+8] = dtmp*rng_isaac_rand_double_norm_pm1();	// im4
		a[ipad+br8[2]+8] = dtmp*rng_isaac_rand_double_norm_pm1();	// re5
		a[ipad+br8[3]+8] = dtmp*rng_isaac_rand_double_norm_pm1();	// im5
		a[ipad+br8[4]+8] = dtmp*rng_isaac_rand_double_norm_pm1();	// re6
		a[ipad+br8[5]+8] = dtmp*rng_isaac_rand_double_norm_pm1();	// im6
		a[ipad+br8[6]+8] = dtmp*rng_isaac_rand_double_norm_pm1();	// re7
		a[ipad+br8[7]+8] = dtmp*rng_isaac_rand_double_norm_pm1();	// im7
	#elif defined(USE_SSE2)
		a[ipad+br4[0]   ] = dtmp*rng_isaac_rand_double_norm_pm1();	// re0
		a[ipad+br4[1]   ] = dtmp*rng_isaac_rand_double_norm_pm1();	// im0
		a[ipad+br4[2]   ] = dtmp*rng_isaac_rand_double_norm_pm1();	// re1
		a[ipad+br4[3]   ] = dtmp*rng_isaac_rand_double_norm_pm1();	// im1
		a[ipad+br4[0]+ 4] = dtmp*rng_isaac_rand_double_norm_pm1();	// re2
		a[ipad+br4[1]+ 4] = dtmp*rng_isaac_rand_double_norm_pm1();	// im2
		a[ipad+br4[2]+ 4] = dtmp*rng_isaac_rand_double_norm_pm1();	// re3
		a[ipad+br4[3]+ 4] = dtmp*rng_isaac_rand_double_norm_pm1();	// im3
		a[ipad+br4[0]+ 8] = dtmp*rng_isaac_rand_double_norm_pm1();	// re4
		a[ipad+br4[1]+ 8] = dtmp*rng_isaac_rand_double_norm_pm1();	// im4
		a[ipad+br4[2]+ 8] = dtmp*rng_isaac_rand_double_norm_pm1();	// re5
		a[ipad+br4[3]+ 8] = dtmp*rng_isaac_rand_double_norm_pm1();	// im5
		a[ipad+br4[0]+12] = dtmp*rng_isaac_rand_double_norm_pm1();	// re6
		a[ipad+br4[1]+12] = dtmp*rng_isaac_rand_double_norm_pm1();	// im6
		a[ipad+br4[2]+12] = dtmp*rng_isaac_rand_double_norm_pm1();	// re7
		a[ipad+br4[3]+12] = dtmp*rng_isaac_rand_double_norm_pm1();	// im7
	#else
		a[ipad+0   ] = dtmp*rng_isaac_rand_double_norm_pm1();	// re0
		a[ipad+1   ] = dtmp*rng_isaac_rand_double_norm_pm1();	// im0
		a[ipad+2   ] = dtmp*rng_isaac_rand_double_norm_pm1();	// re1
		a[ipad+3   ] = dtmp*rng_isaac_rand_double_norm_pm1();	// im1
		a[ipad+0+ 4] = dtmp*rng_isaac_rand_double_norm_pm1();	// re2
		a[ipad+1+ 4] = dtmp*rng_isaac_rand_double_norm_pm1();	// im2
		a[ipad+2+ 4] = dtmp*rng_isaac_rand_double_norm_pm1();	// re3
		a[ipad+3+ 4] = dtmp*rng_isaac_rand_double_norm_pm1();	// im3
		a[ipad+0+ 8] = dtmp*rng_isaac_rand_double_norm_pm1();	// re4
		a[ipad+1+ 8] = dtmp*rng_isaac_rand_double_norm_pm1();	// im4
		a[ipad+2+ 8] = dtmp*rng_isaac_rand_double_norm_pm1();	// re5
		a[ipad+3+ 8] = dtmp*rng_isaac_rand_double_norm_pm1();	// im5
		a[ipad+0+12] = dtmp*rng_isaac_rand_double_norm_pm1();	// re6
		a[ipad+1+12] = dtmp*rng_isaac_rand_double_norm_pm1();	// im6
		a[ipad+2+12] = dtmp*rng_isaac_rand_double_norm_pm1();	// re7
		a[ipad+3+12] = dtmp*rng_isaac_rand_double_norm_pm1();	// im7
	#endif
  #if 0	// print DFT inputs in linear array fashion, 1st w.r.to non-SIMD {re,im,re,im,...} layout, then w.r.to actual SIMD layout:
	#ifdef USE_AVX512
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i, 0, a[ipad+br16[ 0]],ipad+ 0, a[ipad+ 0]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i, 1, a[ipad+br16[ 1]],ipad+ 1, a[ipad+ 1]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i, 2, a[ipad+br16[ 2]],ipad+ 2, a[ipad+ 2]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i, 3, a[ipad+br16[ 3]],ipad+ 3, a[ipad+ 3]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i, 4, a[ipad+br16[ 4]],ipad+ 4, a[ipad+ 4]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i, 5, a[ipad+br16[ 5]],ipad+ 5, a[ipad+ 5]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i, 6, a[ipad+br16[ 6]],ipad+ 6, a[ipad+ 6]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i, 7, a[ipad+br16[ 7]],ipad+ 7, a[ipad+ 7]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i, 8, a[ipad+br16[ 8]],ipad+ 8, a[ipad+ 8]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i, 9, a[ipad+br16[ 9]],ipad+ 9, a[ipad+ 9]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i,10, a[ipad+br16[10]],ipad+10, a[ipad+10]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i,11, a[ipad+br16[11]],ipad+11, a[ipad+11]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i,12, a[ipad+br16[12]],ipad+12, a[ipad+12]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i,13, a[ipad+br16[13]],ipad+13, a[ipad+13]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i,14, a[ipad+br16[14]],ipad+14, a[ipad+14]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i,15, a[ipad+br16[15]],ipad+15, a[ipad+15]);
	#elif defined(USE_AVX)
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i,0  ,a[ipad+br8[0]  ],ipad+0  ,a[ipad+0  ]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i,1  ,a[ipad+br8[1]  ],ipad+1  ,a[ipad+1  ]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i,2  ,a[ipad+br8[2]  ],ipad+2  ,a[ipad+2  ]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i,3  ,a[ipad+br8[3]  ],ipad+3  ,a[ipad+3  ]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i,4  ,a[ipad+br8[4]  ],ipad+4  ,a[ipad+4  ]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i,5  ,a[ipad+br8[5]  ],ipad+5  ,a[ipad+5  ]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i,6  ,a[ipad+br8[6]  ],ipad+6  ,a[ipad+6  ]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i,7  ,a[ipad+br8[7]  ],ipad+7  ,a[ipad+7  ]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i,0+8,a[ipad+br8[0]+8],ipad+0+8,a[ipad+0+8]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i,1+8,a[ipad+br8[1]+8],ipad+1+8,a[ipad+1+8]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i,2+8,a[ipad+br8[2]+8],ipad+2+8,a[ipad+2+8]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i,3+8,a[ipad+br8[3]+8],ipad+3+8,a[ipad+3+8]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i,4+8,a[ipad+br8[4]+8],ipad+4+8,a[ipad+4+8]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i,5+8,a[ipad+br8[5]+8],ipad+5+8,a[ipad+5+8]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i,6+8,a[ipad+br8[6]+8],ipad+6+8,a[ipad+6+8]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i,7+8,a[ipad+br8[7]+8],ipad+7+8,a[ipad+7+8]);
	#elif defined(USE_SSE2)
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i,0   ,a[ipad+br4[0]   ],ipad+0   ,a[ipad+0   ]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i,1   ,a[ipad+br4[1]   ],ipad+1   ,a[ipad+1   ]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i,2   ,a[ipad+br4[2]   ],ipad+2   ,a[ipad+2   ]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i,3   ,a[ipad+br4[3]   ],ipad+3   ,a[ipad+3   ]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i,0+ 4,a[ipad+br4[0]+ 4],ipad+0+ 4,a[ipad+0+ 4]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i,1+ 4,a[ipad+br4[1]+ 4],ipad+1+ 4,a[ipad+1+ 4]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i,2+ 4,a[ipad+br4[2]+ 4],ipad+2+ 4,a[ipad+2+ 4]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i,3+ 4,a[ipad+br4[3]+ 4],ipad+3+ 4,a[ipad+3+ 4]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i,0+ 8,a[ipad+br4[0]+ 8],ipad+0+ 8,a[ipad+0+ 8]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i,1+ 8,a[ipad+br4[1]+ 8],ipad+1+ 8,a[ipad+1+ 8]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i,2+ 8,a[ipad+br4[2]+ 8],ipad+2+ 8,a[ipad+2+ 8]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i,3+ 8,a[ipad+br4[3]+ 8],ipad+3+ 8,a[ipad+3+ 8]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i,0+12,a[ipad+br4[0]+12],ipad+0+12,a[ipad+0+12]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i,1+12,a[ipad+br4[1]+12],ipad+1+12,a[ipad+1+12]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i,2+12,a[ipad+br4[2]+12],ipad+2+12,a[ipad+2+12]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i,3+12,a[ipad+br4[3]+12],ipad+3+12,a[ipad+3+12]);
	#else
		printf("A[%3d][%2d] = %20.10e; ipad: %2d\n",i,0   ,a[ipad+0   ],ipad+0   );
		printf("A[%3d][%2d] = %20.10e; ipad: %2d\n",i,1   ,a[ipad+1   ],ipad+1   );
		printf("A[%3d][%2d] = %20.10e; ipad: %2d\n",i,2   ,a[ipad+2   ],ipad+2   );
		printf("A[%3d][%2d] = %20.10e; ipad: %2d\n",i,3   ,a[ipad+3   ],ipad+3   );
		printf("A[%3d][%2d] = %20.10e; ipad: %2d\n",i,0+ 4,a[ipad+0+ 4],ipad+0+ 4);
		printf("A[%3d][%2d] = %20.10e; ipad: %2d\n",i,1+ 4,a[ipad+1+ 4],ipad+1+ 4);
		printf("A[%3d][%2d] = %20.10e; ipad: %2d\n",i,2+ 4,a[ipad+2+ 4],ipad+2+ 4);
		printf("A[%3d][%2d] = %20.10e; ipad: %2d\n",i,3+ 4,a[ipad+3+ 4],ipad+3+ 4);
		printf("A[%3d][%2d] = %20.10e; ipad: %2d\n",i,0+ 8,a[ipad+0+ 8],ipad+0+ 8);
		printf("A[%3d][%2d] = %20.10e; ipad: %2d\n",i,1+ 8,a[ipad+1+ 8],ipad+1+ 8);
		printf("A[%3d][%2d] = %20.10e; ipad: %2d\n",i,2+ 8,a[ipad+2+ 8],ipad+2+ 8);
		printf("A[%3d][%2d] = %20.10e; ipad: %2d\n",i,3+ 8,a[ipad+3+ 8],ipad+3+ 8);
		printf("A[%3d][%2d] = %20.10e; ipad: %2d\n",i,0+12,a[ipad+0+12],ipad+0+12);
		printf("A[%3d][%2d] = %20.10e; ipad: %2d\n",i,1+12,a[ipad+1+12],ipad+1+12);
		printf("A[%3d][%2d] = %20.10e; ipad: %2d\n",i,2+12,a[ipad+2+12],ipad+2+12);
		printf("A[%3d][%2d] = %20.10e; ipad: %2d\n",i,3+12,a[ipad+3+12],ipad+3+12);
	#endif
	if(i+16 >= n) exit(0);	// If printing the above inputs, exit immediately.
  #endif
	}
#endif
/********************************************************/

#ifdef USE_PTHREAD

	// If also using main thread to do work units, that task-dispatch occurs after all the threadpool-task launches:
	for(ithread = 0; ithread < pool_work_units; ithread++)
	{
		task_control.data = (void*)(&tdat[ithread]);
		threadpool_add_task(tpool, &task_control, task_is_blocking);

#else

	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		if(full_pass) maxerr = 0.0;
	#ifdef USE_SSE2
	//	VEC_DBL_INIT(max_err, 0.0);	*** must do this in conjunction with thread-local-data-copy
	#endif

		i      = _i[ithread];	/* Pointer to the BASE and BASEINV arrays.	*/
		jstart = _jstart[ithread];
		jhi    = _jhi[ithread];

		col = _col[ithread];
		co2 = _co2[ithread];
		co3 = _co3[ithread];

		for(l = 0; l < RADIX; l++) {
			bjmodn[l] = _bjmodn[l][ithread];
		}
		/* init carries	*/
	#ifdef USE_AVX512
		tmp = cy;
		for(l = 0; l < RADIX; l += 8, ++tmp) {
			tmp->d0 = _cy[l  ][ithread];
			tmp->d1 = _cy[l+1][ithread];
			tmp->d2 = _cy[l+2][ithread];
			tmp->d3 = _cy[l+3][ithread];
			tmp->d4 = _cy[l+4][ithread];
			tmp->d5 = _cy[l+5][ithread];
			tmp->d6 = _cy[l+6][ithread];
			tmp->d7 = _cy[l+7][ithread];
		}
	#elif defined(USE_AVX)	// AVX and AVX2 both use 256-bit registers
		tmp = cy;
		for(l = 0; l < RADIX; l += 4, ++tmp) {
			tmp->d0 = _cy[l  ][ithread];
			tmp->d1 = _cy[l+1][ithread];
			tmp->d2 = _cy[l+2][ithread];
			tmp->d3 = _cy[l+3][ithread];
		}
	#elif defined(USE_SSE2)
		tmp = cy;
		for(l = 0; l < RADIX; l += 2, ++tmp) {
			tmp->d0 = _cy[l  ][ithread];
			tmp->d1 = _cy[l+1][ithread];
		}
	#else
		for(l = 0; l < RADIX; l++) {
			cy[l] = _cy[l][ithread];
		}
	#endif

		/********************************************************************************/
		/* This main loop is same for un-and-multithreaded, so stick into a header file */
		/* (can't use a macro because of the #if-enclosed stuff).                       */
		/********************************************************************************/
		#include "radix384_main_carry_loop.h"

	#ifdef USE_AVX512
		tmp = cy;
		for(l = 0; l < RADIX; l += 8, ++tmp) {
			_cy[l  ][ithread] = tmp->d0;
			_cy[l+1][ithread] = tmp->d1;
			_cy[l+2][ithread] = tmp->d2;
			_cy[l+3][ithread] = tmp->d3;
			_cy[l+4][ithread] = tmp->d4;
			_cy[l+5][ithread] = tmp->d5;
			_cy[l+6][ithread] = tmp->d6;
			_cy[l+7][ithread] = tmp->d7;
		}
		if(full_pass) {
			t0 = MAX(max_err->d0,max_err->d1);
			t1 = MAX(max_err->d2,max_err->d3);
			t2 = MAX(max_err->d4,max_err->d5);
			t3 = MAX(max_err->d6,max_err->d7);
			maxerr = MAX( MAX(t0,t1), MAX(t2,t3) );
		}
	#elif defined(USE_AVX)	// AVX and AVX2 both use 256-bit registers
		tmp = cy;
		for(l = 0; l < RADIX; l += 4, ++tmp) {
			_cy[l  ][ithread] = tmp->d0;
			_cy[l+1][ithread] = tmp->d1;
			_cy[l+2][ithread] = tmp->d2;
			_cy[l+3][ithread] = tmp->d3;
		}
		if(full_pass) maxerr = MAX( MAX(max_err->d0,max_err->d1) , MAX(max_err->d2,max_err->d3) );
	#elif defined(USE_SSE2)
		tmp = cy;
		for(l = 0; l < RADIX; l += 2, ++tmp) {
			_cy[l  ][ithread] = tmp->d0;
			_cy[l+1][ithread] = tmp->d1;
		}
		if(full_pass) maxerr = MAX(max_err->d0,max_err->d1);
	#else
		for(l = 0; l < RADIX; l++) {
			_cy[l][ithread] = cy[l];
		}
	#endif

  #endif	// #ifdef USE_PTHREAD

	}	/******* END OF PARALLEL FOR-LOOP ********/

#ifdef USE_PTHREAD	// End of threadpool-based dispatch: Add a small wait-loop to ensure all threads complete

  #if 0//def OS_TYPE_MACOSX

	/*** Main execution thread executes remaining chunks in serial fashion (but in || with the pool threads): ***/
	for(j = 0; j < main_work_units; ++j)
	{
	//	printf("adding main task %d\n",j + pool_work_units);
		ASSERT(0x0 == cy384_process_chunk( (void*)(&tdat[j + pool_work_units]) ), "Main-thread task failure!");
	}

  #endif

	struct timespec ns_time;	// We want a sleep interval of 0.1 mSec here...
	ns_time.tv_sec  =      0;	// (time_t)seconds - Don't use this because under OS X it's of type __darwin_time_t, which is long rather than double as under most linux distros
	ns_time.tv_nsec = 100000;	// (long)nanoseconds - Get our desired 0.1 mSec as 10^5 nSec here

	while(tpool && tpool->free_tasks_queue.num_tasks != pool_work_units) {
		ASSERT(0 == mlucas_nanosleep(&ns_time), "nanosleep fail!");
	}

	/* Copy the thread-specific output carry data back to shared memory: */
	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		if(maxerr < tdat[ithread].maxerr) {
			maxerr = tdat[ithread].maxerr;
		}
		for(l = 0; l < RADIX; l++) {
			_cy[l][ithread] = tdat[ithread].cy[l];
		}
	}
#endif

	if(full_pass) {
	//	printf("Iter = %d, maxerr = %20.15f\n",iter,maxerr);
	} else {
		break;
	}

	/*   Wraparound carry cleanup loop is here:

	The cleanup carries from the end of each length-N/RADIX set of contiguous data into the begining of the next
	can all be neatly processed as follows:

	(1) Invert the forward DIF FFT of the first block of RADIX complex elements in A and unweight;
	(2) Propagate cleanup carries among the real and imaginary parts of the RADIX outputs of (1);
	(3) Reweight and perform a forward DIF FFT on the result of (2);
	(4) If any of the exit carries from (2) are nonzero, advance to the next RADIX elements and repeat (1-4).
	*/
	for(l = 0; l < RADIX; l++) {
		t[l].re = _cy[l][CY_THREADS - 1];
	}
	for(ithread = CY_THREADS - 1; ithread > 0; ithread--)
	{
		for(l = 0; l < RADIX; l++) {
			_cy[l][ithread] = _cy[l][ithread-1];
		}
	}
	_cy[0][0] =+t[RADIX-1].re;	/* ...The wraparound carry is here: */
	for(l = 1; l < RADIX; l++) {
		_cy[l][0] = t[l-1].re;
	}

	full_pass = 0;
	scale = prp_mult = 1;
	j_jhi = jhi_wrap;

	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		for(j = ithread*pini; j <= ithread*pini + j_jhi; j++)
		{
			// Generate padded version of j, since prepadding pini is thread-count unsafe:
			j1 = j + ( (j >> DAT_BITS) << PAD_BITS );
			for(l = 0; l < RADIX>>2; l++) {
				jt = j1 + poff[l];
				a[jt   ] *= radix_inv;
				a[jt+p1] *= radix_inv;
				a[jt+p2] *= radix_inv;
				a[jt+p3] *= radix_inv;
			}
		}
	}
}	/* endfor(outer) */

	dtmp = 0;
	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		for(l = 0; l < RADIX; l++) {
			dtmp += fabs(_cy[l][ithread]);
		}
		*fracmax = maxerr;
	}
	if(dtmp != 0.0)
	{
		sprintf(cbuf,"ERROR: iter = %10d; nonzero exit carry in %s - input wordsize may be too small.\n",iter,func);
		mlucas_fprint(cbuf,INTERACT);
		err = ERR_CARRY;
		return(err);
	}

	return(0);
}

/***************/

void radix384_dif_pass1(double a[], int n)
{
/*
!...Acronym: DIF = Decimation In Frequency
!
!...Subroutine to perform an initial radix-384 complex DIF FFT pass on the data in the length-N real vector A.
!
!   See the documentation in radix16_dif_pass for further details on storage and indexing.
!   See the documentation in radix3_dif_pass for details on the radix-3 subtransforms.
!   See the documentation in radix384_dif_pass for details on the twiddleless coprime-radix-DFTS setup.
*/
	int j,j1,j2,jp,jt;
	int k,l,l1,l2,k0,k1,k2;
	static int dif_triplets[72];
	// In order to preserve length-2 numeric-index-offset property here, use hex for digit:
	static int NDIVR, first_entry=TRUE,
		p1,p2,p3,p4,p5,p6,p7,p8,p9,pa,pb,pc,pd,pe,pf,
	     p010,p020,p030,p040,p050,p060,p070,p080,p090,p0a0,p0b0,p0c0,p0d0,p0e0,p0f0,p100,p110,p120,p130,p140,p150,p160,p170;
	const double c3m1 = -1.50000000000000000000, s = 0.86602540378443864675;	// cos(twopi/3)-1, sin(twopi/3)
	static int i_offsets[128], o_offsets[384];
	// Local storage: We must use an array here because scalars have no guarantees about relative address offsets
	// [and even if those are contiguous-as-hoped-for, they may run in reverse]; Make array type (struct complex)
	// to allow us to use the same offset-indexing as in the original radix-32 in-place DFT macros:
	struct complex t[RADIX];
	double t00,t01,t02,t03,t04,t05;

	if(!first_entry && (n/RADIX) != NDIVR)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		first_entry=FALSE;
		NDIVR = n/RADIX;

		/*   constant index offsets for array load/stores are here.	*/
		p1 = NDIVR;		p010 = NDIVR<<4;	p100 = NDIVR<<8;
		p2 = p1 + p1;	p020 = p010 + p010;	p110 = p100 + p010;
		p3 = p2 + p1;	p030 = p020 + p010;	p120 = p110 + p010;
		p4 = p3 + p1;	p040 = p030 + p010;	p130 = p120 + p010;
		p5 = p4 + p1;	p050 = p040 + p010;	p140 = p130 + p010;
		p6 = p5 + p1;	p060 = p050 + p010;	p150 = p140 + p010;
		p7 = p6 + p1;	p070 = p060 + p010;	p160 = p150 + p010;
		p8 = p7 + p1;	p080 = p070 + p010;	p170 = p160 + p010;
		p9 = p8 + p1;	p090 = p080 + p010;
		pa = p9 + p1;	p0a0 = p090 + p010;
		pb = pa + p1;	p0b0 = p0a0 + p010;
		pc = pb + p1;	p0c0 = p0b0 + p010;
		pd = pc + p1;	p0d0 = p0c0 + p010;
		pe = pd + p1;	p0e0 = p0d0 + p010;
		pf = pe + p1;	p0f0 = p0e0 + p010;

		p1 += ( (p1 >> DAT_BITS) << PAD_BITS );		p010 += ( (p010 >> DAT_BITS) << PAD_BITS );		p100 += ( (p100 >> DAT_BITS) << PAD_BITS );
		p2 += ( (p2 >> DAT_BITS) << PAD_BITS );		p020 += ( (p020 >> DAT_BITS) << PAD_BITS );		p110 += ( (p110 >> DAT_BITS) << PAD_BITS );
		p3 += ( (p3 >> DAT_BITS) << PAD_BITS );		p030 += ( (p030 >> DAT_BITS) << PAD_BITS );		p120 += ( (p120 >> DAT_BITS) << PAD_BITS );
		p4 += ( (p4 >> DAT_BITS) << PAD_BITS );		p040 += ( (p040 >> DAT_BITS) << PAD_BITS );		p130 += ( (p130 >> DAT_BITS) << PAD_BITS );
		p5 += ( (p5 >> DAT_BITS) << PAD_BITS );		p050 += ( (p050 >> DAT_BITS) << PAD_BITS );		p140 += ( (p140 >> DAT_BITS) << PAD_BITS );
		p6 += ( (p6 >> DAT_BITS) << PAD_BITS );		p060 += ( (p060 >> DAT_BITS) << PAD_BITS );		p150 += ( (p150 >> DAT_BITS) << PAD_BITS );
		p7 += ( (p7 >> DAT_BITS) << PAD_BITS );		p070 += ( (p070 >> DAT_BITS) << PAD_BITS );		p160 += ( (p160 >> DAT_BITS) << PAD_BITS );
		p8 += ( (p8 >> DAT_BITS) << PAD_BITS );		p080 += ( (p080 >> DAT_BITS) << PAD_BITS );		p170 += ( (p170 >> DAT_BITS) << PAD_BITS );
		p9 += ( (p9 >> DAT_BITS) << PAD_BITS );		p090 += ( (p090 >> DAT_BITS) << PAD_BITS );
		pa += ( (pa >> DAT_BITS) << PAD_BITS );		p0a0 += ( (p0a0 >> DAT_BITS) << PAD_BITS );
		pb += ( (pb >> DAT_BITS) << PAD_BITS );		p0b0 += ( (p0b0 >> DAT_BITS) << PAD_BITS );
		pc += ( (pc >> DAT_BITS) << PAD_BITS );		p0c0 += ( (p0c0 >> DAT_BITS) << PAD_BITS );
		pd += ( (pd >> DAT_BITS) << PAD_BITS );		p0d0 += ( (p0d0 >> DAT_BITS) << PAD_BITS );
		pe += ( (pe >> DAT_BITS) << PAD_BITS );		p0e0 += ( (p0e0 >> DAT_BITS) << PAD_BITS );
		pf += ( (pf >> DAT_BITS) << PAD_BITS );		p0f0 += ( (p0f0 >> DAT_BITS) << PAD_BITS );

	// Set array offsets for radix-128 inputs [these are same for all 3 such DFTs].
		for(j = 0; j < 0x80; ++j) { i_offsets[j] = j<<1; }

	// For the radix-128 outputs we need the following offsets:
	#if 0	// For initial run, used to determine the needed output index scramble; same ascending 128-set reused 3 times:
		o_offsets[0x00] =  0;
		o_offsets[0x01] = p1;
		o_offsets[0x02] = p2;
		o_offsets[0x03] = p3;
		o_offsets[0x04] = p4;
		o_offsets[0x05] = p5;
		o_offsets[0x06] = p6;
		o_offsets[0x07] = p7;
		o_offsets[0x08] = p8;
		o_offsets[0x09] = p9;
		o_offsets[0x0a] = pa;
		o_offsets[0x0b] = pb;
		o_offsets[0x0c] = pc;
		o_offsets[0x0d] = pd;
		o_offsets[0x0e] = pe;
		o_offsets[0x0f] = pf;
		for(j = 0; j < 0x10; ++j) { o_offsets[j+0x010] = o_offsets[j] + p010; }
		for(j = 0; j < 0x10; ++j) { o_offsets[j+0x020] = o_offsets[j] + p020; }
		for(j = 0; j < 0x10; ++j) { o_offsets[j+0x030] = o_offsets[j] + p030; }
		for(j = 0; j < 0x10; ++j) { o_offsets[j+0x040] = o_offsets[j] + p040; }
		for(j = 0; j < 0x10; ++j) { o_offsets[j+0x050] = o_offsets[j] + p050; }
		for(j = 0; j < 0x10; ++j) { o_offsets[j+0x060] = o_offsets[j] + p060; }
		for(j = 0; j < 0x10; ++j) { o_offsets[j+0x070] = o_offsets[j] + p070; }

		for(j = 0; j < 0x80; ++j) { o_offsets[j+0x080] = o_offsets[j]; }
		for(j = 0; j < 0x80; ++j) { o_offsets[j+0x100] = o_offsets[j]; }
	#else
	/*	Required output index permutation = [
		000,001,002,003,005,004,007,006,00a,00b,009,008,00f,00e,00c,00d = p000 + p01235476ab98fecd[perm0]
		015,014,017,016,012,013,011,010,01f,01e,01c,01d,019,018,01b,01a = p010 + p54762310fecd98ba[perm1]
		02a,02b,029,028,02f,02e,02c,02d,025,024,027,026,022,023,021,020 = p020 + pab98fecd54762310[perm2]
		03f,03e,03c,03d,039,038,03b,03a,032,033,031,030,037,036,034,035 = p030 + pfecd98ba23107645[perm3]
		055,054,057,056,052,053,051,050,05f,05e,05c,05d,059,058,05b,05a = p050 + [perm1]
		04a,04b,049,048,04f,04e,04c,04d,045,044,047,046,042,043,041,040 = p040 + [perm2]
		07f,07e,07c,07d,079,078,07b,07a,072,073,071,070,077,076,074,075 = p070 + [perm3]
		065,064,067,066,062,063,061,060,06f,06e,06c,06d,069,068,06b,06a = p060 + [perm1]
		12a,12b,129,128,12f,12e,12c,12d,125,124,127,126,122,123,121,120 = p120 + [perm2]	<<< In order to effect
		13f,13e,13c,13d,139,138,13b,13a,132,133,131,130,137,136,134,135 = p130 + [perm3]	<<< the swap of p100-0f0
		115,114,117,116,112,113,111,110,11f,11e,11c,11d,119,118,11b,11a = p110 + [perm1]	<<< with p080-0f0 in the
		10a,10b,109,108,10f,10e,10c,10d,105,104,107,106,102,103,101,100 = p100 + [perm2]	<<< final pair of 128-DIF
		17f,17e,17c,17d,179,178,17b,17a,172,173,171,170,177,176,174,175 = p170 + [perm3]	<<< outputs, swap a+j1+p80
		165,164,167,166,162,163,161,160,16f,16e,16c,16d,169,168,16b,16a = p160 + [perm1]	<<< a+j1+p100 in that pair
		14a,14b,149,148,14f,14e,14c,14d,145,144,147,146,142,143,141,140 = p140 + [perm2]	<<< of function calls, w.r.to
		15f,15e,15c,15d,159,158,15b,15a,152,153,151,150,157,156,154,155 = p150 + [perm3]	<<< the initial perm-finding run.
		0d5,0d4,0d7,0d6,0d2,0d3,0d1,0d0,0df,0de,0dc,0dd,0d9,0d8,0db,0da = p0d0 + [perm1]
		0ca,0cb,0c9,0c8,0cf,0ce,0cc,0cd,0c5,0c4,0c7,0c6,0c2,0c3,0c1,0c0 = p0c0 + [perm2]
		0ff,0fe,0fc,0fd,0f9,0f8,0fb,0fa,0f2,0f3,0f1,0f0,0f7,0f6,0f4,0f5 = p0f0 + [perm3]
		0e5,0e4,0e7,0e6,0e2,0e3,0e1,0e0,0ef,0ee,0ec,0ed,0e9,0e8,0eb,0ea = p0e0 + [perm1]
		0aa,0ab,0a9,0a8,0af,0ae,0ac,0ad,0a5,0a4,0a7,0a6,0a2,0a3,0a1,0a0 = p0a0 + [perm2]
		0bf,0be,0bc,0bd,0b9,0b8,0bb,0ba,0b2,0b3,0b1,0b0,0b7,0b6,0b4,0b5 = p0b0 + [perm3]
		095,094,097,096,092,093,091,090,09f,09e,09c,09d,099,098,09b,09a = p090 + [perm1]
		08a,08b,089,088,08f,08e,08c,08d,085,084,087,086,082,083,081,080 = p080 + [perm2]
		] */
		const int perm0[16] = {0,p1,p2,p3,p5,p4,p7,p6,pa,pb,p9,p8,pf,pe,pc,pd};
		const int perm1[16] = {p5,p4,p7,p6,p2,p3,p1,0,pf,pe,pc,pd,p9,p8,pb,pa};
		const int perm2[16] = {pa,pb,p9,p8,pf,pe,pc,pd,p5,p4,p7,p6,p2,p3,p1,0};
		const int perm3[16] = {pf,pe,pc,pd,p9,p8,pb,pa,p2,p3,p1,0,p7,p6,p4,p5};
		for(j = 0; j < 0x10; ++j) { o_offsets[j] = perm0[j]; }
		for(j = 0; j < 0x10; ++j) { o_offsets[j+0x010] = p010 + perm1[j]; }
		for(j = 0; j < 0x10; ++j) { o_offsets[j+0x020] = p020 + perm2[j]; }
		for(j = 0; j < 0x10; ++j) { o_offsets[j+0x030] = p030 + perm3[j]; }
		for(j = 0; j < 0x10; ++j) { o_offsets[j+0x040] = p050 + perm1[j]; }
		for(j = 0; j < 0x10; ++j) { o_offsets[j+0x050] = p040 + perm2[j]; }
		for(j = 0; j < 0x10; ++j) { o_offsets[j+0x060] = p070 + perm3[j]; }
		for(j = 0; j < 0x10; ++j) { o_offsets[j+0x070] = p060 + perm1[j]; }
	// o_offsets are w.r.to each respective block of 128 complex data, i.e. p-indexing is (mod p80):
		for(j = 0; j < 0x10; ++j) { o_offsets[j+0x080] = p020 + perm2[j]; }
		for(j = 0; j < 0x10; ++j) { o_offsets[j+0x090] = p030 + perm3[j]; }
		for(j = 0; j < 0x10; ++j) { o_offsets[j+0x0a0] = p010 + perm1[j]; }
		for(j = 0; j < 0x10; ++j) { o_offsets[j+0x0b0] =    0 + perm2[j]; }
		for(j = 0; j < 0x10; ++j) { o_offsets[j+0x0c0] = p070 + perm3[j]; }
		for(j = 0; j < 0x10; ++j) { o_offsets[j+0x0d0] = p060 + perm1[j]; }
		for(j = 0; j < 0x10; ++j) { o_offsets[j+0x0e0] = p040 + perm2[j]; }
		for(j = 0; j < 0x10; ++j) { o_offsets[j+0x0f0] = p050 + perm3[j]; }
	// o_offsets are w.r.to each respective block of 128 complex data, i.e. p-indexing is (mod p80):
		for(j = 0; j < 0x10; ++j) { o_offsets[j+0x100] = p050 + perm1[j]; }
		for(j = 0; j < 0x10; ++j) { o_offsets[j+0x110] = p040 + perm2[j]; }
		for(j = 0; j < 0x10; ++j) { o_offsets[j+0x120] = p070 + perm3[j]; }
		for(j = 0; j < 0x10; ++j) { o_offsets[j+0x130] = p060 + perm1[j]; }
		for(j = 0; j < 0x10; ++j) { o_offsets[j+0x140] = p020 + perm2[j]; }
		for(j = 0; j < 0x10; ++j) { o_offsets[j+0x150] = p030 + perm3[j]; }
		for(j = 0; j < 0x10; ++j) { o_offsets[j+0x160] = p010 + perm1[j]; }
		for(j = 0; j < 0x10; ++j) { o_offsets[j+0x170] =    0 + perm2[j]; }
	#endif
	// Index-high-bits triplets needed for compact-obj-code scheme:
		k = 0;
		dif_triplets[k] = p170; dif_triplets[k+1] = p0f0; dif_triplets[k+2] = p070; k += 3;
		dif_triplets[k] = p160; dif_triplets[k+1] = p0e0; dif_triplets[k+2] = p060; k += 3;
		dif_triplets[k] = p150; dif_triplets[k+1] = p0d0; dif_triplets[k+2] = p050; k += 3;
		dif_triplets[k] = p140; dif_triplets[k+1] = p0c0; dif_triplets[k+2] = p040; k += 3;
		dif_triplets[k] = p130; dif_triplets[k+1] = p0b0; dif_triplets[k+2] = p030; k += 3;
		dif_triplets[k] = p120; dif_triplets[k+1] = p0a0; dif_triplets[k+2] = p020; k += 3;
		dif_triplets[k] = p110; dif_triplets[k+1] = p090; dif_triplets[k+2] = p010; k += 3;
		dif_triplets[k] = p100; dif_triplets[k+1] = p080; dif_triplets[k+2] =    0; k += 3;
		dif_triplets[k] = p0f0; dif_triplets[k+1] = p070; dif_triplets[k+2] = p170; k += 3;
		dif_triplets[k] = p0e0; dif_triplets[k+1] = p060; dif_triplets[k+2] = p160; k += 3;
		dif_triplets[k] = p0d0; dif_triplets[k+1] = p050; dif_triplets[k+2] = p150; k += 3;
		dif_triplets[k] = p0c0; dif_triplets[k+1] = p040; dif_triplets[k+2] = p140; k += 3;
		dif_triplets[k] = p0b0; dif_triplets[k+1] = p030; dif_triplets[k+2] = p130; k += 3;
		dif_triplets[k] = p0a0; dif_triplets[k+1] = p020; dif_triplets[k+2] = p120; k += 3;
		dif_triplets[k] = p090; dif_triplets[k+1] = p010; dif_triplets[k+2] = p110; k += 3;
		dif_triplets[k] = p080; dif_triplets[k+1] =    0; dif_triplets[k+2] = p100; k += 3;
		dif_triplets[k] = p070; dif_triplets[k+1] = p170; dif_triplets[k+2] = p0f0; k += 3;
		dif_triplets[k] = p060; dif_triplets[k+1] = p160; dif_triplets[k+2] = p0e0; k += 3;
		dif_triplets[k] = p050; dif_triplets[k+1] = p150; dif_triplets[k+2] = p0d0; k += 3;
		dif_triplets[k] = p040; dif_triplets[k+1] = p140; dif_triplets[k+2] = p0c0; k += 3;
		dif_triplets[k] = p030; dif_triplets[k+1] = p130; dif_triplets[k+2] = p0b0; k += 3;
		dif_triplets[k] = p020; dif_triplets[k+1] = p120; dif_triplets[k+2] = p0a0; k += 3;
		dif_triplets[k] = p010; dif_triplets[k+1] = p110; dif_triplets[k+2] = p090; k += 3;
		dif_triplets[k] =    0; dif_triplets[k+1] = p100; dif_triplets[k+2] = p080;
	}

/*...The radix-384 pass is here.	*/

	for(j = 0; j < NDIVR; j += 2)
	{
	#ifdef USE_AVX512
		j1 = (j & mask03) + br16[j&15];
	#elif defined(USE_AVX)
		j1 = (j & mask02) + br8[j&7];
	#elif defined(USE_SSE2)
		j1 = (j & mask01) + br4[j&3];
	#else
		j1 = j;
	#endif
		j1 += ( (j1>> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
		j2 = j1 + RE_IM_STRIDE;

	/*
	Twiddleless version arranges 128 sets of radix-3 DFT inputs as follows: 0 in upper left corner,
	decrement 128 horizontally and 3 vertically, all (mod 384). Break into 4 cols for readability,
	and insert blank rows to reveal 16-macro-calls subgroupings exploited by compact-obj-code scheme below:

		DIF/DIT input-scramble array = [
		000,100,080	<*** Wrap this around to end...

		170,0f0,070 + pd		110,090,010 + pd		0b0,030,130 + pd		050,150,0d0 + pd
		170,0f0,070 + pa		110,090,010 + pa		0b0,030,130 + pa		050,150,0d0 + pa
		170,0f0,070 + p7		110,090,010 + p7		0b0,030,130 + p7		050,150,0d0 + p7
		170,0f0,070 + p4		110,090,010 + p4		0b0,030,130 + p4		050,150,0d0 + p4
		170,0f0,070 + p1		110,090,010 + p1		0b0,030,130 + p1		050,150,0d0 + p1
		160,0e0,060 + pe		100,080,000 + pe		0a0,020,120 + pe		040,140,0c0 + pe
		160,0e0,060 + pb		100,080,000 + pb		0a0,020,120 + pb		040,140,0c0 + pb
		160,0e0,060 + p8		100,080,000 + p8		0a0,020,120 + p8		040,140,0c0 + p8
		160,0e0,060 + p5		100,080,000 + p5		0a0,020,120 + p5		040,140,0c0 + p5
		160,0e0,060 + p2		100,080,000 + p2		0a0,020,120 + p2		040,140,0c0 + p2
		150,0d0,050 + pf		0f0,070,170 + pf		090,010,110 + pf		030,130,0b0 + pf
		150,0d0,050 + pc		0f0,070,170 + pc		090,010,110 + pc		030,130,0b0 + pc
		150,0d0,050 + p9		0f0,070,170 + p9		090,010,110 + p9		030,130,0b0 + p9
		150,0d0,050 + p6		0f0,070,170 + p6		090,010,110 + p6		030,130,0b0 + p6
		150,0d0,050 + p3		0f0,070,170 + p3		090,010,110 + p3		030,130,0b0 + p3
		150,0d0,050 + p0		0f0,070,170 + p0		090,010,110 + p0		030,130,0b0 + p0

		140,0c0,040 + pd		0e0,060,160 + pd		080,000,100 + pd		020,120,0a0 + pd
		140,0c0,040 + pa		0e0,060,160 + pa		080,000,100 + pa		020,120,0a0 + pa
		140,0c0,040 + p7		0e0,060,160 + p7		080,000,100 + p7		020,120,0a0 + p7
		140,0c0,040 + p4		0e0,060,160 + p4		080,000,100 + p4		020,120,0a0 + p4
		140,0c0,040 + p1		0e0,060,160 + p1		080,000,100 + p1		020,120,0a0 + p1
		130,0b0,030 + pe		0d0,050,150 + pe		070,170,0f0 + pe		010,110,090 + pe
		130,0b0,030 + pb		0d0,050,150 + pb		070,170,0f0 + pb		010,110,090 + pb
		130,0b0,030 + p8		0d0,050,150 + p8		070,170,0f0 + p8		010,110,090 + p8
		130,0b0,030 + p5		0d0,050,150 + p5		070,170,0f0 + p5		010,110,090 + p5
		130,0b0,030 + p2		0d0,050,150 + p2		070,170,0f0 + p2		010,110,090 + p2
		120,0a0,020 + pf		0c0,040,140 + pf		060,160,0e0 + pf		000,100,080 + pf
		120,0a0,020 + pc		0c0,040,140 + pc		060,160,0e0 + pc		000,100,080 + pc
		120,0a0,020 + p9		0c0,040,140 + p9		060,160,0e0 + p9		000,100,080 + p9
		120,0a0,020 + p6		0c0,040,140 + p6		060,160,0e0 + p6		000,100,080 + p6
		120,0a0,020 + p3		0c0,040,140 + p3		060,160,0e0 + p3		000,100,080 + p3
		120,0a0,020 + p0		0c0,040,140 + p0		060,160,0e0 + p0		000,100,080 + p0	<*** 0-term wrapped ***
		[cont. in col2]			[cont. in col3]			[cont. in col4]

	To handle the wraparound of the 0-term we just need a bit of mod-128 indexing magic.
	*/
	/*...gather the needed data (384 64-bit complex) and do 128 radix-3 transforms - We want unit-strides in the radix384-DFT macro, so use large output strides here: */
		// Loop-based compact-obj-code impl exploits above index pattern to group interior [sandwiched between single leading
		// and 15 trailing, which also get fused into a final group] macro calls into sets of 16 with neatly cutoff index groupings:
		l = 1; l1 = l+128; l2 = l+256;	// Skip 0-term, which gets saved for wraparound
		for(k = 0; k < 72; k += 3) {
			k0 = dif_triplets[k]; k1 = dif_triplets[k+1]; k2 = dif_triplets[k+2]; k += 3;
			jt = j1 + pd; jp = j2 + pd;	RADIX_03_DFT(s,c3m1, a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2], t00,t01,t02,t03,t04,t05, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im); ++l; ++l1; ++l2;
			jt = j1 + pa; jp = j2 + pa;	RADIX_03_DFT(s,c3m1, a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2], t00,t01,t02,t03,t04,t05, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im); ++l; ++l1; ++l2;
			jt = j1 + p7; jp = j2 + p7;	RADIX_03_DFT(s,c3m1, a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2], t00,t01,t02,t03,t04,t05, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im); ++l; ++l1; ++l2;
			jt = j1 + p4; jp = j2 + p4;	RADIX_03_DFT(s,c3m1, a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2], t00,t01,t02,t03,t04,t05, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im); ++l; ++l1; ++l2;
			jt = j1 + p1; jp = j2 + p1;	RADIX_03_DFT(s,c3m1, a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2], t00,t01,t02,t03,t04,t05, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im); ++l; ++l1; ++l2;
			k0 = dif_triplets[k]; k1 = dif_triplets[k+1]; k2 = dif_triplets[k+2]; k += 3;
			jt = j1 + pe; jp = j2 + pe;	RADIX_03_DFT(s,c3m1, a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2], t00,t01,t02,t03,t04,t05, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im); ++l; ++l1; ++l2;
			jt = j1 + pb; jp = j2 + pb;	RADIX_03_DFT(s,c3m1, a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2], t00,t01,t02,t03,t04,t05, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im); ++l; ++l1; ++l2;
			jt = j1 + p8; jp = j2 + p8;	RADIX_03_DFT(s,c3m1, a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2], t00,t01,t02,t03,t04,t05, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im); ++l; ++l1; ++l2;
			jt = j1 + p5; jp = j2 + p5;	RADIX_03_DFT(s,c3m1, a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2], t00,t01,t02,t03,t04,t05, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im); ++l; ++l1; ++l2;
			jt = j1 + p2; jp = j2 + p2;	RADIX_03_DFT(s,c3m1, a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2], t00,t01,t02,t03,t04,t05, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im); ++l; ++l1; ++l2;
			k0 = dif_triplets[k]; k1 = dif_triplets[k+1]; k2 = dif_triplets[k+2];	// Skip k-incr here since loop control handles this one
			jt = j1 + pf; jp = j2 + pf;	RADIX_03_DFT(s,c3m1, a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2], t00,t01,t02,t03,t04,t05, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im); ++l; ++l1; ++l2;
			jt = j1 + pc; jp = j2 + pc;	RADIX_03_DFT(s,c3m1, a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2], t00,t01,t02,t03,t04,t05, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im); ++l; ++l1; ++l2;
			jt = j1 + p9; jp = j2 + p9;	RADIX_03_DFT(s,c3m1, a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2], t00,t01,t02,t03,t04,t05, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im); ++l; ++l1; ++l2;
			jt = j1 + p6; jp = j2 + p6;	RADIX_03_DFT(s,c3m1, a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2], t00,t01,t02,t03,t04,t05, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im); ++l; ++l1; ++l2;
			jt = j1 + p3; jp = j2 + p3;	RADIX_03_DFT(s,c3m1, a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2], t00,t01,t02,t03,t04,t05, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im); ++l; l &= 0x7f; l1 = l+128; l2 = l+256;	// <*** needed for final-loop-pass wraparound of 0-term
										RADIX_03_DFT(s,c3m1, a[j1+k0],a[j2+k0],a[j1+k1],a[j2+k1],a[j1+k2],a[j2+k2], t00,t01,t02,t03,t04,t05, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im); ++l; ++l1; ++l2;
		}

	/*...and now do 3 radix-128 transforms, first assuming in-order output index offsets, then using the
	resulting data-mismatch tables produced by test_fft_radix() to derive the needed output permutation:
	*/
		// NOTE: Due to casting-to-double of inputs temp-array pointer Radix-32 macro doesn't play nice with r-array
		// offsets of form "r+32", so make address-taking explicit and wrapping in () prior to macro code execution:
		//
		// Since i_data are in a local-complex scratch array, need to override any non-unity SIMD-re/im data stride with 1 for inputs:
		//                                         vvv
		RADIX_128_DIF((double *)(t+0x000),i_offsets,1, (a+j1     ),o_offsets      ,RE_IM_STRIDE);	// Inputs in t[ 00- 7f]
		RADIX_128_DIF((double *)(t+0x080),i_offsets,1, (a+j1+p100),o_offsets+0x080,RE_IM_STRIDE);	// Inputs in t[ 80- ff]
		RADIX_128_DIF((double *)(t+0x100),i_offsets,1, (a+j1+p080),o_offsets+0x100,RE_IM_STRIDE);	// Inputs in t[100-17f]
	}
}

/***************/

void radix384_dit_pass1(double a[], int n)
{
/*
!...Acronym: DIT = Decimation In Time
!
!...Subroutine to perform an initial radix-384 complex DIT FFT pass on the data in the length-N real vector A.
!
!   See the documentation in radix16_dif_pass for further details on storage and indexing.
*/
	int j,j1,j2,jp,jt;
	int k,l,l1,l2,k0,k1,k2;
	static int dit_triplets[24];	// Only need 1/3 as many here as for DIF
	// In order to preserve length-2 numeric-index-offset property here, use hex for digit:
	static int NDIVR, first_entry=TRUE,
		p1,p2,p3,p4,p5,p6,p7,p8,p9,pa,pb,pc,pd,pe,pf,
	     p010,p020,p030,p040,p050,p060,p070,p080,p090,p0a0,p0b0,p0c0,p0d0,p0e0,p0f0,p100,p110,p120,p130,p140,p150,p160,p170;
	const double c3m1 = -1.50000000000000000000, s = 0.86602540378443864675;	// cos(twopi/3)-1, sin(twopi/3)
	static int i_offsets[384], o_offsets[128];
	// Local storage: We must use an array here because scalars have no guarantees about relative address offsets
	// [and even if those are contiguous-as-hoped-for, they may run in reverse]; Make array type (struct complex)
	// to allow us to use the same offset-indexing as in the original radix-32 in-place DFT macros:
	struct complex t[RADIX];
	double t00,t01,t02,t03,t04,t05;

	if(!first_entry && (n/RADIX) != NDIVR)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		first_entry=FALSE;
		NDIVR = n/RADIX;

		/*   constant index offsets for array load/stores are here.	*/
		p1 = NDIVR;		p010 = NDIVR<<4;	p100 = NDIVR<<8;
		p2 = p1 + p1;	p020 = p010 + p010;	p110 = p100 + p010;
		p3 = p2 + p1;	p030 = p020 + p010;	p120 = p110 + p010;
		p4 = p3 + p1;	p040 = p030 + p010;	p130 = p120 + p010;
		p5 = p4 + p1;	p050 = p040 + p010;	p140 = p130 + p010;
		p6 = p5 + p1;	p060 = p050 + p010;	p150 = p140 + p010;
		p7 = p6 + p1;	p070 = p060 + p010;	p160 = p150 + p010;
		p8 = p7 + p1;	p080 = p070 + p010;	p170 = p160 + p010;
		p9 = p8 + p1;	p090 = p080 + p010;
		pa = p9 + p1;	p0a0 = p090 + p010;
		pb = pa + p1;	p0b0 = p0a0 + p010;
		pc = pb + p1;	p0c0 = p0b0 + p010;
		pd = pc + p1;	p0d0 = p0c0 + p010;
		pe = pd + p1;	p0e0 = p0d0 + p010;
		pf = pe + p1;	p0f0 = p0e0 + p010;

		p1 += ( (p1 >> DAT_BITS) << PAD_BITS );		p010 += ( (p010 >> DAT_BITS) << PAD_BITS );		p100 += ( (p100 >> DAT_BITS) << PAD_BITS );
		p2 += ( (p2 >> DAT_BITS) << PAD_BITS );		p020 += ( (p020 >> DAT_BITS) << PAD_BITS );		p110 += ( (p110 >> DAT_BITS) << PAD_BITS );
		p3 += ( (p3 >> DAT_BITS) << PAD_BITS );		p030 += ( (p030 >> DAT_BITS) << PAD_BITS );		p120 += ( (p120 >> DAT_BITS) << PAD_BITS );
		p4 += ( (p4 >> DAT_BITS) << PAD_BITS );		p040 += ( (p040 >> DAT_BITS) << PAD_BITS );		p130 += ( (p130 >> DAT_BITS) << PAD_BITS );
		p5 += ( (p5 >> DAT_BITS) << PAD_BITS );		p050 += ( (p050 >> DAT_BITS) << PAD_BITS );		p140 += ( (p140 >> DAT_BITS) << PAD_BITS );
		p6 += ( (p6 >> DAT_BITS) << PAD_BITS );		p060 += ( (p060 >> DAT_BITS) << PAD_BITS );		p150 += ( (p150 >> DAT_BITS) << PAD_BITS );
		p7 += ( (p7 >> DAT_BITS) << PAD_BITS );		p070 += ( (p070 >> DAT_BITS) << PAD_BITS );		p160 += ( (p160 >> DAT_BITS) << PAD_BITS );
		p8 += ( (p8 >> DAT_BITS) << PAD_BITS );		p080 += ( (p080 >> DAT_BITS) << PAD_BITS );		p170 += ( (p170 >> DAT_BITS) << PAD_BITS );
		p9 += ( (p9 >> DAT_BITS) << PAD_BITS );		p090 += ( (p090 >> DAT_BITS) << PAD_BITS );
		pa += ( (pa >> DAT_BITS) << PAD_BITS );		p0a0 += ( (p0a0 >> DAT_BITS) << PAD_BITS );
		pb += ( (pb >> DAT_BITS) << PAD_BITS );		p0b0 += ( (p0b0 >> DAT_BITS) << PAD_BITS );
		pc += ( (pc >> DAT_BITS) << PAD_BITS );		p0c0 += ( (p0c0 >> DAT_BITS) << PAD_BITS );
		pd += ( (pd >> DAT_BITS) << PAD_BITS );		p0d0 += ( (p0d0 >> DAT_BITS) << PAD_BITS );
		pe += ( (pe >> DAT_BITS) << PAD_BITS );		p0e0 += ( (p0e0 >> DAT_BITS) << PAD_BITS );
		pf += ( (pf >> DAT_BITS) << PAD_BITS );		p0f0 += ( (p0f0 >> DAT_BITS) << PAD_BITS );

	// Set array offsets for radix-128 outputs [these are same for all 3 such DFTs].
		for(j = 0; j < 0x80; ++j) { o_offsets[j] = j<<1; }

		const int perm0[16] = {0,p1,p3,p2,p7,p6,p5,p4,pf,pe,pd,pc,pb,pa,p9,p8};
		const int perm1[16] = {pf,pe,pd,pc,pb,pa,p9,p8,p7,p6,p5,p4,p3,p2,p1,0};
		const int perm2[16] = {p5,p4,p6,p7,p1,0,p2,p3,p9,p8,pa,pb,pe,pf,pc,pd};
		const int perm3[16] = {p9,p8,pa,pb,pe,pf,pc,pd,p1,0,p2,p3,p6,p7,p4,p5};
		const int perm4[16] = {pa,pb,p8,p9,pc,pd,pf,pe,p2,p3,0,p1,p4,p5,p7,p6};
		const int perm5[16] = {p2,p3,0,p1,p4,p5,p7,p6,pc,pd,pf,pe,p8,p9,pb,pa};
		for(j = 0; j < 0x10; ++j) { i_offsets[j] = perm0[j]; }
		for(j = 0; j < 0x10; ++j) { i_offsets[j+0x010] = p010 + perm1[j]; }
		for(j = 0; j < 0x10; ++j) { i_offsets[j+0x020] = p030 + perm1[j]; }
		for(j = 0; j < 0x10; ++j) { i_offsets[j+0x030] = p020 + perm1[j]; }
		for(j = 0; j < 0x10; ++j) { i_offsets[j+0x040] = p070 + perm1[j]; }
		for(j = 0; j < 0x10; ++j) { i_offsets[j+0x050] = p060 + perm1[j]; }
		for(j = 0; j < 0x10; ++j) { i_offsets[j+0x060] = p050 + perm1[j]; }
		for(j = 0; j < 0x10; ++j) { i_offsets[j+0x070] = p040 + perm1[j]; }
	// i_offsets are w.r.to each respective block of 128 complex data, i.e. p-indexing is (mod p80):
		for(j = 0; j < 0x10; ++j) { i_offsets[j+0x080] = p050 + perm2[j]; }
		for(j = 0; j < 0x10; ++j) { i_offsets[j+0x090] = p040 + perm2[j]; }
		for(j = 0; j < 0x10; ++j) { i_offsets[j+0x0a0] = p060 + perm2[j]; }
		for(j = 0; j < 0x10; ++j) { i_offsets[j+0x0b0] = p070 + perm3[j]; }
		for(j = 0; j < 0x10; ++j) { i_offsets[j+0x0c0] = p010 + perm2[j]; }
		for(j = 0; j < 0x10; ++j) { i_offsets[j+0x0d0] =    0 + perm2[j]; }
		for(j = 0; j < 0x10; ++j) { i_offsets[j+0x0e0] = p020 + perm2[j]; }
		for(j = 0; j < 0x10; ++j) { i_offsets[j+0x0f0] = p030 + perm3[j]; }
	// i_offsets are w.r.to each respective block of 128 complex data, i.e. p-indexing is (mod p80):
		for(j = 0; j < 0x10; ++j) { i_offsets[j+0x100] = p020 + perm4[j]; }
		for(j = 0; j < 0x10; ++j) { i_offsets[j+0x110] = p030 + perm5[j]; }
		for(j = 0; j < 0x10; ++j) { i_offsets[j+0x120] =    0 + perm4[j]; }
		for(j = 0; j < 0x10; ++j) { i_offsets[j+0x130] = p010 + perm5[j]; }
		for(j = 0; j < 0x10; ++j) { i_offsets[j+0x140] = p040 + perm4[j]; }
		for(j = 0; j < 0x10; ++j) { i_offsets[j+0x150] = p050 + perm5[j]; }
		for(j = 0; j < 0x10; ++j) { i_offsets[j+0x160] = p070 + perm5[j]; }
		for(j = 0; j < 0x10; ++j) { i_offsets[j+0x170] = p060 + perm5[j]; }

	// Index-high-bits triplets needed for compact-obj-code scheme:
		k = 0;
	#if 0
		dit_triplets[k] =    0; dit_triplets[k+1] = p080; dit_triplets[k+2] = p100; k += 3;
		dit_triplets[k] = p010; dit_triplets[k+1] = p090; dit_triplets[k+2] = p110; k += 3;
		dit_triplets[k] = p020; dit_triplets[k+1] = p0a0; dit_triplets[k+2] = p120; k += 3;
		dit_triplets[k] = p030; dit_triplets[k+1] = p0b0; dit_triplets[k+2] = p130; k += 3;
		dit_triplets[k] = p040; dit_triplets[k+1] = p0c0; dit_triplets[k+2] = p140; k += 3;
		dit_triplets[k] = p050; dit_triplets[k+1] = p0d0; dit_triplets[k+2] = p150; k += 3;
		dit_triplets[k] = p060; dit_triplets[k+1] = p0e0; dit_triplets[k+2] = p160; k += 3;
		dit_triplets[k] = p070; dit_triplets[k+1] = p0f0; dit_triplets[k+2] = p170;
	#else
		/* Cf. comments below re. this needed reordering w.r.to the orginal ordering in the middle column:
		Loop	p*0-offsets	need reorder
		----	-----------	------------
		0		070,0f0,170	-
		1		060,0e0,160	160,060,0e0
		2		050,0d0,150	0d0,150,050
		3		040,0c0,140	-
		4		030,0b0,130	130,030,0b0
		5		020,0a0,120	0a0,120,020
		6		010,090,110	-
		7		000,080,100	100,000,080
		*/
		dit_triplets[k] = p070; dit_triplets[k+1] = p0f0; dit_triplets[k+2] = p170; k += 3;
		dit_triplets[k] = p160; dit_triplets[k+1] = p060; dit_triplets[k+2] = p0e0; k += 3;
		dit_triplets[k] = p0d0; dit_triplets[k+1] = p150; dit_triplets[k+2] = p050; k += 3;
		dit_triplets[k] = p040; dit_triplets[k+1] = p0c0; dit_triplets[k+2] = p140; k += 3;
		dit_triplets[k] = p130; dit_triplets[k+1] = p030; dit_triplets[k+2] = p0b0; k += 3;
		dit_triplets[k] = p0a0; dit_triplets[k+1] = p120; dit_triplets[k+2] = p020; k += 3;
		dit_triplets[k] = p010; dit_triplets[k+1] = p090; dit_triplets[k+2] = p110; k += 3;
		dit_triplets[k] = p100; dit_triplets[k+1] =    0; dit_triplets[k+2] = p080;
	#endif
	}

/*...The radix-384 pass is here.	*/

	for(j = 0; j < NDIVR; j += 2)
	{
	#ifdef USE_AVX512
		j1 = (j & mask03) + br16[j&15];
	#elif defined(USE_AVX)
		j1 = (j & mask02) + br8[j&7];
	#elif defined(USE_SSE2)
		j1 = (j & mask01) + br4[j&3];
	#else
		j1 = j;
	#endif
		j1 += ( (j1>> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
		j2 = j1 + RE_IM_STRIDE;
	/*
	Twiddleless version uses same linear-index-vector-form permutation as in DIF -
	Remember, inputs to DIT are bit-reversed, so using output of test_fft_radix(),
	store the 3 index-offset 128-tets going into the radix-128 DFTs in the i_offset array.
	Combined DIT input-scramble array =
	000,001,003,002,007,006,005,004,00f,00e,00d,00c,00b,00a,009,008 = p000 + p01327654fedcba98[perm0]
	01f,01e,01d,01c,01b,01a,019,018,017,016,015,014,013,012,011,010 = p010 + pfedcba9876543210[perm1]
	03f,03e,03d,03c,03b,03a,039,038,037,036,035,034,033,032,031,030 = p030 + perm1
	02f,02e,02d,02c,02b,02a,029,028,027,026,025,024,023,022,021,020 = p020 + perm1
	07f,07e,07d,07c,07b,07a,079,078,077,076,075,074,073,072,071,070 = p070 + perm1
	06f,06e,06d,06c,06b,06a,069,068,067,066,065,064,063,062,061,060 = p060 + perm1
	05f,05e,05d,05c,05b,05a,059,058,057,056,055,054,053,052,051,050 = p050 + perm1
	04f,04e,04d,04c,04b,04a,049,048,047,046,045,044,043,042,041,040 = p040 + perm1

	0d5,0d4,0d6,0d7,0d1,0d0,0d2,0d3,0d9,0d8,0da,0db,0de,0df,0dc,0dd = p0d0 + p5467102398abefcd[perm2]
	0c5,0c4,0c6,0c7,0c1,0c0,0c2,0c3,0c9,0c8,0ca,0cb,0ce,0cf,0cc,0cd = p0c0 + perm2
	0e5,0e4,0e6,0e7,0e1,0e0,0e2,0e3,0e9,0e8,0ea,0eb,0ee,0ef,0ec,0ed = p0e0 + perm2
	0f9,0f8,0fa,0fb,0fe,0ff,0fc,0fd,0f1,0f0,0f2,0f3,0f6,0f7,0f4,0f5 = p0f0 + p98abefcd10236745[perm3]
	095,094,096,097,091,090,092,093,099,098,09a,09b,09e,09f,09c,09d = p090 + perm2
	085,084,086,087,081,080,082,083,089,088,08a,08b,08e,08f,08c,08d = p080 + perm2
	0a5,0a4,0a6,0a7,0a1,0a0,0a2,0a3,0a9,0a8,0aa,0ab,0ae,0af,0ac,0ad = p0a0 + perm2
	0b9,0b8,0ba,0bb,0be,0bf,0bc,0bd,0b1,0b0,0b2,0b3,0b6,0b7,0b4,0b5 = p0b0 + perm3

	12a,12b,128,129,12c,12d,12f,12e,122,123,120,121,124,125,127,126 = p120 + pab89cdfe23014576[perm4]
	132,133,130,131,134,135,137,136,13c,13d,13f,13e,138,139,13b,13a = p130 + p23014576cdfe89ba[perm5]
	10a,10b,108,109,10c,10d,10f,10e,102,103,100,101,104,105,107,106 = p100 + perm4
	112,113,110,111,114,115,117,116,11c,11d,11f,11e,118,119,11b,11a = p110 + perm5
	14a,14b,148,149,14c,14d,14f,14e,142,143,140,141,144,145,147,146 = p140 + perm4
	152,153,150,151,154,155,157,156,15c,15d,15f,15e,158,159,15b,15a = p150 + perm5
	172,173,170,171,174,175,177,176,17c,17d,17f,17e,178,179,17b,17a = p170 + perm5
	162,163,160,161,164,165,167,166,16c,16d,16f,16e,168,169,16b,16a = p160 + perm5, thus perm0-5 get used 1,7,6,2,3,5 times, respectively.
	*/
	/*...gather the needed data (384 64-bit complex) and do 3 radix-128 transforms,	*/
		// Since o_data are in a local-complex scratch array, need to override any non-unity SIMD-re/im data stride with 1 for outputs:
		//                                                                                   vvv
		RADIX_128_DIT((a+j1     ),i_offsets      ,RE_IM_STRIDE, (double *)(t+0x000),o_offsets,1);	// Outputs in t[ 00- 7f]
		RADIX_128_DIT((a+j1+p080),i_offsets+0x080,RE_IM_STRIDE, (double *)(t+0x080),o_offsets,1);	// Outputs in t[ 80- ff]
		RADIX_128_DIT((a+j1+p100),i_offsets+0x100,RE_IM_STRIDE, (double *)(t+0x100),o_offsets,1);	// Outputs in t[100-17f]
	/*
	We start by putting output triplets into [0,p080,p100]-strided array locations, with the offset of the first elt. of each such triplet incrementing by p1. That gives
	Required output index permutation = [
	000,0ff,07e,17d,0fc,07b,17a,0f9,078,177,0f6,075,174,0f3,072,171,0f0,06f,16e,0ed,06c,16b,0ea,069,168,0e7,066,165,0e4,063,162,0e1,060,15f,0de,05d,15c,0db,05a,159,0d8,057,156,0d5,054,153,0d2,051,150,0cf,04e,14d,0cc,04b,14a,0c9,048,147,0c6,045,144,0c3,042,141,0c0,03f,13e,0bd,03c,13b,0ba,039,138,0b7,036,135,0b4,033,132,0b1,030,12f,0ae,02d,12c,0ab,02a,129,0a8,027,126,0a5,024,123,0a2,021,120,09f,01e,11d,09c,01b,11a,099,018,117,096,015,114,093,012,111,090,00f,10e,08d,00c,10b,08a,009,108,087,006,105,084,003,102,081,
	100,07f,17e,0fd,07c,17b,0fa,079,178,0f7,076,175,0f4,073,172,0f1,070,16f,0ee,06d,16c,0eb,06a,169,0e8,067,166,0e5,064,163,0e2,061,160,0df,05e,15d,0dc,05b,15a,0d9,058,157,0d6,055,154,0d3,052,151,0d0,04f,14e,0cd,04c,14b,0ca,049,148,0c7,046,145,0c4,043,142,0c1,040,13f,0be,03d,13c,0bb,03a,139,0b8,037,136,0b5,034,133,0b2,031,130,0af,02e,12d,0ac,02b,12a,0a9,028,127,0a6,025,124,0a3,022,121,0a0,01f,11e,09d,01c,11b,09a,019,118,097,016       ,115,094,013,112,091,010,10f,08e,00d,10c,08b,00a,109,088,007,106,085,004,103,082,001,
	080,17f,0fe,07d,17c,0fb,07a,179,0f8,077,176,0f5,074,173,0f2,071,170,0ef,06e,16d,0ec,06b,16a,0e9,068,167,0e6,065,164,0e3,062,161,0e0,05f,15e,0dd,05c,15b,0da,059,158,0d7,056,155,0d4,053,152,0d1,050,14f,0ce,04d,14c,0cb,04a,149,0c8,047,146,0c5,044,143,0c2,041,140,0bf,03e,13d,0bc,03b,13a,0b9,038,137,0b6,035,134,0b3,032,131,0b0,02f,12e,0ad,02c,12b,0aa,029,128,0a7,026,125,0a4,023,122,0a1,020,11f,09e,01d,11c,09b,01a,119,098,017,116,095,014,113,092,011,110,08f,00e,10d,08c,00b,10a,089,008,107,086,005,104,083,002,101,
	]
	That means the first-triplet ouptuts (leftmost 3-elt col) are in the right spot, just need to swap 100 <--> 080, i.e. last 2 elts of
	that triplet. For the remaining triplets we notice that the high-indexed elt. of each follows a descending pattern, 17f,17e,...,002,001.
	So need to rejigger the output-triplet indexing in the loop to mirror that; to handle the wraparound of the 0-term we just need a bit of mod-128 indexing magic.
	That gives, with | marking boundaries of data blocks processed by each loop pass:
	Required output index permutation = [
	000,081,102,003,084,105,006,087,108,009,08a,10b,00c,08d,10e,00f|090,111,012,093,114,015,096,117|018,099,11a,01b,09c,11d,01e,09f|120,021,0a2,123,024,0a5,126,027,0a8,129,02a,0ab,12c,02d,0ae,12f|030,0b1,132,033,0b4,135,036,0b7,138,039,0ba,13b,03c,0bd,13e,03f|0c0,141,042,0c3,144,045,0c6,147,048,0c9,14a,04b,0cc,14d,04e,0cf|150,051,0d2,153,054,0d5,156,057,0d8,159,05a,0db,15c,05d,0de,15f|060,0e1,162,063,0e4,165,066,0e7,168,069,0ea,16b,06c,0ed,16e,06f|0f0,171,072,0f3,174,075,0f6,177,078,0f9,17a,07b,0fc,17d,07e,0ff,
	100,001,082,103,004,085,106,007,088,109,00a,08b,10c,00d,08e,10f|010,091,112,013,094,115,016,097|118,019,09a,11b,01c,09d,11e,01f|0a0,121,022,0a3,124,025,0a6,127,028,0a9,12a,02b,0ac,12d,02e,0af|130,031,0b2,133,034,0b5,136,037,0b8,139,03a,0bb,13c,03d,0be,13f|040,0c1,142,043,0c4,145,046,0c7,148,049,0ca,14b,04c,0cd,14e,04f|0d0,151,052,0d3,154,055,0d6,157,058,0d9,15a,05b,0dc,15d,05e,0df|160,061,0e2,163,064,0e5,166,067,0e8,169,06a,0eb,16c,06d,0ee,16f|070,0f1,172,073,0f4,175,076,0f7,178,079,0fa,17b,07c,0fd,17e,07f,
	080,101,002,083,104,005,086,107,008,089,10a,00b,08c,10d,00e,08f|110,011,092,113,014,095,116,017|098,119,01a,09b,11c,01d,09e,11f|020,0a1,122,023,0a4,125,026,0a7,128,029,0aa,12b,02c,0ad,12e,02f|0b0,131,032,0b3,134,035,0b6,137,038,0b9,13a,03b,0bc,13d,03e,0bf|140,041,0c2,143,044,0c5,146,047,0c8,149,04a,0cb,14c,04d,0ce,14f|050,0d1,152,053,0d4,155,056,0d7,158,059,0da,15b,05c,0dd,15e,05f|0e0,161,062,0e3,164,065,0e6,167,068,0e9,16a,06b,0ec,16d,06e,0ef|170,071,0f2,173,074,0f5,176,077,0f8,179,07a,0fb,17c,07d,0fe,17f,
	]
	Now consider the k-index triplets of the outputs of the 16 3-DFTs on the first loop execution:
		Loop processes triplets:	Need triplets:	Needed output k-indices order:
		0: l = 1, k0-2 =  7f, ff,17f	f,7,17		1,0,2
		1: l = 2, k0-2 =  7e, fe,17e	7,17,f		0,2,1
		2: l = 3, k0-2 =  7d, fd,17d	17,f,7		2,1,0
		4: l = 5, k0-2 =  7b, fb,17b	7,17,f		0,2,1
		5: l = 6, k0-2 =  7a, fa,17a	17,f,7		2,1,0
		6: l = 7, k0-2 =  79, f9,179	f,7,17		1,0,2
		7: l = 8, k0-2 =  78, f8,178	7,17,f		0,2,1
		3: l = 4, k0-2 =  7c, fc,17c	f,7,17		1,0,2
		8: l = 9, k0-2 =  77, f7,177	17,f,7		2,1,0
		9: l = a, k0-2 =  76, f6,176	f,7,17		1,0,2
		a: l = b, k0-2 =  75, f5,175	7,17,f		0,2,1
		b: l = c, k0-2 =  74, f4,174	17,f,7		2,1,0
		c: l = d, k0-2 =  73, f3,173	f,7,17		1,0,2
		d: l = e, k0-2 =  72, f2,172	7,17,f		0,2,1
		e: l = f, k0-2 =  71, f1,171	17,f,7		2,1,0
		f: l =10, k0-2 =  70, f0,170	f,7,17		1,0,2
	That gets us closer, but it turns out that each loop pass has a different needed pattern of such [0,1,2]-perms:
	Required output index permutation = [
								Loop 7:															Loop 6:															Loop 5:															Loop 4:															Loop 3:															Loop 2:															Loop 1:															Loop 0:
	100,101,102,103,104,105,106,107,108,109,10a,10b,10c,10d,10e,10f|010,011,012,013,014,015,016,017,018,019,01a,01b,01c,01d,01e,01f|0a0,0a1,0a2,0a3,0a4,0a5,0a6,0a7,0a8,0a9,0aa,0ab,0ac,0ad,0ae,0af|130,131,132,133,134,135,136,137,138,139,13a,13b,13c,13d,13e,13f|040,041,042,043,044,045,046,047,048,049,04a,04b,04c,04d,04e,04f|0d0,0d1,0d2,0d3,0d4,0d5,0d6,0d7,0d8,0d9,0da,0db,0dc,0dd,0de,0df|160,161,162,163,164,165,166,167,168,169,16a,16b,16c,16d,16e,16f|070,071,072,073,074,075,076,077,078,079,07a,07b,07c,07d,07e,07f,
	000,001,002,003,004,005,006,007,008,009,00a,00b,00c,00d,00e,00f|090,091,092,093,094,095,096,097,098,099,09a,09b,09c,09d,09e,09f|120,121,122,123,124,125,126,127,128,129,12a,12b,12c,12d,12e,12f|030,031,032,033,034,035,036,037,038,039,03a,03b,03c,03d,03e,03f|0c0,0c1,0c2,0c3,0c4,0c5,0c6,0c7,0c8,0c9,0ca,0cb,0cc,0cd,0ce,0cf|150,151,152,153,154,155,156,157,158,159,15a,15b,15c,15d,15e,15f|060,061,062,063,064,065,066,067,068,069,06a,06b,06c,06d,06e,06f|0f0,0f1,0f2,0f3,0f4,0f5,0f6,0f7,0f8,0f9,0fa,0fb,0fc,0fd,0fe,0ff,
	080,081,082,083,084,085,086,087,088,089,08a,08b,08c,08d,08e,08f|110,111,112,113,114,115,116,117,118,119,11a,11b,11c,11d,11e,11f|020,021,022,023,024,025,026,027,028,029,02a,02b,02c,02d,02e,02f|0b0,0b1,0b2,0b3,0b4,0b5,0b6,0b7,0b8,0b9,0ba,0bb,0bc,0bd,0be,0bf|140,141,142,143,144,145,146,147,148,149,14a,14b,14c,14d,14e,14f|050,051,052,053,054,055,056,057,058,059,05a,05b,05c,05d,05e,05f|0e0,0e1,0e2,0e3,0e4,0e5,0e6,0e7,0e8,0e9,0ea,0eb,0ec,0ed,0ee,0ef|170,171,172,173,174,175,176,177,178,179,17a,17b,17c,17d,17e,17f,
	]
	Thus we need the following block-level swaps of the p*0 indices - remember, w.r.to the above ascending-index data ordering
	the loop pass runs from 0-7 right-to-left:
		Loop	p*0-offsets		need reorder (simply the leftmost 3-elt col in each of the above loop data-blocks)
		----	-----------		------------
		0		070,0f0,170		-
		1		060,0e0,160		160,060,0e0
		2		050,0d0,150		0d0,150,050
		3		040,0c0,140		-
		4		030,0b0,130		130,030,0b0
		5		020,0a0,120		0a0,120,020
		6		010,090,110		-
		7		000,080,100		100,000,080
	*/
	/*...and now do 128 radix-3 transforms: */
		// Loop-based compact-obj-code impl exploits above index pattern to group interior [sandwiched
		// between single leading and 15 trailing] macro calls into sets of 16 with neatly cutoff index groupings:
		#if 0
		l = 0; l1 = l+128; l2 = l+256;
		for(k = 0; k < 24; k += 3) {
			k0 = dit_triplets[k]; k1 = dit_triplets[k+1]; k2 = dit_triplets[k+2];
										RADIX_03_DFT(s,c3m1, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im, t00,t01,t02,t03,t04,t05, a[j1+k0],a[j2+k0],a[j1+k1],a[j2+k1],a[j1+k2],a[j2+k2]); ++l; ++l1; ++l2;
			jt = j1 + p1; jp = j2 + p1;	RADIX_03_DFT(s,c3m1, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im, t00,t01,t02,t03,t04,t05, a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2]); ++l; ++l1; ++l2;
			jt = j1 + p2; jp = j2 + p2;	RADIX_03_DFT(s,c3m1, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im, t00,t01,t02,t03,t04,t05, a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2]); ++l; ++l1; ++l2;
			jt = j1 + p3; jp = j2 + p3;	RADIX_03_DFT(s,c3m1, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im, t00,t01,t02,t03,t04,t05, a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2]); ++l; ++l1; ++l2;
			jt = j1 + p4; jp = j2 + p4;	RADIX_03_DFT(s,c3m1, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im, t00,t01,t02,t03,t04,t05, a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2]); ++l; ++l1; ++l2;
			jt = j1 + p5; jp = j2 + p5;	RADIX_03_DFT(s,c3m1, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im, t00,t01,t02,t03,t04,t05, a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2]); ++l; ++l1; ++l2;
			jt = j1 + p6; jp = j2 + p6;	RADIX_03_DFT(s,c3m1, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im, t00,t01,t02,t03,t04,t05, a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2]); ++l; ++l1; ++l2;
			jt = j1 + p7; jp = j2 + p7;	RADIX_03_DFT(s,c3m1, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im, t00,t01,t02,t03,t04,t05, a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2]); ++l; ++l1; ++l2;
			jt = j1 + p8; jp = j2 + p8;	RADIX_03_DFT(s,c3m1, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im, t00,t01,t02,t03,t04,t05, a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2]); ++l; ++l1; ++l2;
			jt = j1 + p9; jp = j2 + p9;	RADIX_03_DFT(s,c3m1, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im, t00,t01,t02,t03,t04,t05, a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2]); ++l; ++l1; ++l2;
			jt = j1 + pa; jp = j2 + pa;	RADIX_03_DFT(s,c3m1, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im, t00,t01,t02,t03,t04,t05, a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2]); ++l; ++l1; ++l2;
			jt = j1 + pb; jp = j2 + pb;	RADIX_03_DFT(s,c3m1, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im, t00,t01,t02,t03,t04,t05, a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2]); ++l; ++l1; ++l2;
			jt = j1 + pc; jp = j2 + pc;	RADIX_03_DFT(s,c3m1, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im, t00,t01,t02,t03,t04,t05, a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2]); ++l; ++l1; ++l2;
			jt = j1 + pd; jp = j2 + pd;	RADIX_03_DFT(s,c3m1, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im, t00,t01,t02,t03,t04,t05, a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2]); ++l; ++l1; ++l2;
			jt = j1 + pe; jp = j2 + pe;	RADIX_03_DFT(s,c3m1, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im, t00,t01,t02,t03,t04,t05, a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2]); ++l; ++l1; ++l2;
			jt = j1 + pf; jp = j2 + pf;	RADIX_03_DFT(s,c3m1, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im, t00,t01,t02,t03,t04,t05, a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2]); ++l; ++l1; ++l2;
		}
		#else
		l = 1; l1 = l+128; l2 = l+256;	// Skip 0-term, which gets saved for wraparound
		for(k = 0; k < 24; k += 3) {
			k0 = dit_triplets[k]; k1 = dit_triplets[k+1]; k2 = dit_triplets[k+2];
			jt = j1 + pf; jp = j2 + pf;	RADIX_03_DFT(s,c3m1, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im, t00,t01,t02,t03,t04,t05, a[jt+k1],a[jp+k1],a[jt+k0],a[jp+k0],a[jt+k2],a[jp+k2]);	/*printf("%2u.%x: %2x, k0-2 = %3x,%3x,%3x, a0-5 = %15.10f,%15.10f, %15.10f,%15.10f, %15.10f,%15.10f\n",k/3,0x0,l,(jt+k0)/p1,(jt+k1)/p1,(jt+k2)/p1,a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2]);*/ ++l; ++l1; ++l2;
			jt = j1 + pe; jp = j2 + pe;	RADIX_03_DFT(s,c3m1, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im, t00,t01,t02,t03,t04,t05, a[jt+k0],a[jp+k0],a[jt+k2],a[jp+k2],a[jt+k1],a[jp+k1]);	/*printf("%2u.%x: %2x, k0-2 = %3x,%3x,%3x, a0-5 = %15.10f,%15.10f, %15.10f,%15.10f, %15.10f,%15.10f\n",k/3,0x1,l,(jt+k0)/p1,(jt+k1)/p1,(jt+k2)/p1,a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2]);*/ ++l; ++l1; ++l2;
			jt = j1 + pd; jp = j2 + pd;	RADIX_03_DFT(s,c3m1, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im, t00,t01,t02,t03,t04,t05, a[jt+k2],a[jp+k2],a[jt+k1],a[jp+k1],a[jt+k0],a[jp+k0]);	/*printf("%2u.%x: %2x, k0-2 = %3x,%3x,%3x, a0-5 = %15.10f,%15.10f, %15.10f,%15.10f, %15.10f,%15.10f\n",k/3,0x2,l,(jt+k0)/p1,(jt+k1)/p1,(jt+k2)/p1,a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2]);*/ ++l; ++l1; ++l2;
			jt = j1 + pc; jp = j2 + pc;	RADIX_03_DFT(s,c3m1, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im, t00,t01,t02,t03,t04,t05, a[jt+k1],a[jp+k1],a[jt+k0],a[jp+k0],a[jt+k2],a[jp+k2]);	/*printf("%2u.%x: %2x, k0-2 = %3x,%3x,%3x, a0-5 = %15.10f,%15.10f, %15.10f,%15.10f, %15.10f,%15.10f\n",k/3,0x3,l,(jt+k0)/p1,(jt+k1)/p1,(jt+k2)/p1,a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2]);*/ ++l; ++l1; ++l2;
			jt = j1 + pb; jp = j2 + pb;	RADIX_03_DFT(s,c3m1, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im, t00,t01,t02,t03,t04,t05, a[jt+k0],a[jp+k0],a[jt+k2],a[jp+k2],a[jt+k1],a[jp+k1]);	/*printf("%2u.%x: %2x, k0-2 = %3x,%3x,%3x, a0-5 = %15.10f,%15.10f, %15.10f,%15.10f, %15.10f,%15.10f\n",k/3,0x4,l,(jt+k0)/p1,(jt+k1)/p1,(jt+k2)/p1,a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2]);*/ ++l; ++l1; ++l2;
			jt = j1 + pa; jp = j2 + pa;	RADIX_03_DFT(s,c3m1, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im, t00,t01,t02,t03,t04,t05, a[jt+k2],a[jp+k2],a[jt+k1],a[jp+k1],a[jt+k0],a[jp+k0]);	/*printf("%2u.%x: %2x, k0-2 = %3x,%3x,%3x, a0-5 = %15.10f,%15.10f, %15.10f,%15.10f, %15.10f,%15.10f\n",k/3,0x5,l,(jt+k0)/p1,(jt+k1)/p1,(jt+k2)/p1,a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2]);*/ ++l; ++l1; ++l2;
			jt = j1 + p9; jp = j2 + p9;	RADIX_03_DFT(s,c3m1, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im, t00,t01,t02,t03,t04,t05, a[jt+k1],a[jp+k1],a[jt+k0],a[jp+k0],a[jt+k2],a[jp+k2]);	/*printf("%2u.%x: %2x, k0-2 = %3x,%3x,%3x, a0-5 = %15.10f,%15.10f, %15.10f,%15.10f, %15.10f,%15.10f\n",k/3,0x6,l,(jt+k0)/p1,(jt+k1)/p1,(jt+k2)/p1,a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2]);*/ ++l; ++l1; ++l2;
			jt = j1 + p8; jp = j2 + p8;	RADIX_03_DFT(s,c3m1, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im, t00,t01,t02,t03,t04,t05, a[jt+k0],a[jp+k0],a[jt+k2],a[jp+k2],a[jt+k1],a[jp+k1]);	/*printf("%2u.%x: %2x, k0-2 = %3x,%3x,%3x, a0-5 = %15.10f,%15.10f, %15.10f,%15.10f, %15.10f,%15.10f\n",k/3,0x7,l,(jt+k0)/p1,(jt+k1)/p1,(jt+k2)/p1,a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2]);*/ ++l; ++l1; ++l2;
			jt = j1 + p7; jp = j2 + p7;	RADIX_03_DFT(s,c3m1, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im, t00,t01,t02,t03,t04,t05, a[jt+k2],a[jp+k2],a[jt+k1],a[jp+k1],a[jt+k0],a[jp+k0]);	/*printf("%2u.%x: %2x, k0-2 = %3x,%3x,%3x, a0-5 = %15.10f,%15.10f, %15.10f,%15.10f, %15.10f,%15.10f\n",k/3,0x8,l,(jt+k0)/p1,(jt+k1)/p1,(jt+k2)/p1,a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2]);*/ ++l; ++l1; ++l2;
			jt = j1 + p6; jp = j2 + p6;	RADIX_03_DFT(s,c3m1, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im, t00,t01,t02,t03,t04,t05, a[jt+k1],a[jp+k1],a[jt+k0],a[jp+k0],a[jt+k2],a[jp+k2]);	/*printf("%2u.%x: %2x, k0-2 = %3x,%3x,%3x, a0-5 = %15.10f,%15.10f, %15.10f,%15.10f, %15.10f,%15.10f\n",k/3,0x9,l,(jt+k0)/p1,(jt+k1)/p1,(jt+k2)/p1,a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2]);*/ ++l; ++l1; ++l2;
			jt = j1 + p5; jp = j2 + p5;	RADIX_03_DFT(s,c3m1, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im, t00,t01,t02,t03,t04,t05, a[jt+k0],a[jp+k0],a[jt+k2],a[jp+k2],a[jt+k1],a[jp+k1]);	/*printf("%2u.%x: %2x, k0-2 = %3x,%3x,%3x, a0-5 = %15.10f,%15.10f, %15.10f,%15.10f, %15.10f,%15.10f\n",k/3,0xa,l,(jt+k0)/p1,(jt+k1)/p1,(jt+k2)/p1,a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2]);*/ ++l; ++l1; ++l2;
			jt = j1 + p4; jp = j2 + p4;	RADIX_03_DFT(s,c3m1, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im, t00,t01,t02,t03,t04,t05, a[jt+k2],a[jp+k2],a[jt+k1],a[jp+k1],a[jt+k0],a[jp+k0]);	/*printf("%2u.%x: %2x, k0-2 = %3x,%3x,%3x, a0-5 = %15.10f,%15.10f, %15.10f,%15.10f, %15.10f,%15.10f\n",k/3,0xb,l,(jt+k0)/p1,(jt+k1)/p1,(jt+k2)/p1,a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2]);*/ ++l; ++l1; ++l2;
			jt = j1 + p3; jp = j2 + p3;	RADIX_03_DFT(s,c3m1, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im, t00,t01,t02,t03,t04,t05, a[jt+k1],a[jp+k1],a[jt+k0],a[jp+k0],a[jt+k2],a[jp+k2]);	/*printf("%2u.%x: %2x, k0-2 = %3x,%3x,%3x, a0-5 = %15.10f,%15.10f, %15.10f,%15.10f, %15.10f,%15.10f\n",k/3,0xc,l,(jt+k0)/p1,(jt+k1)/p1,(jt+k2)/p1,a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2]);*/ ++l; ++l1; ++l2;
			jt = j1 + p2; jp = j2 + p2;	RADIX_03_DFT(s,c3m1, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im, t00,t01,t02,t03,t04,t05, a[jt+k0],a[jp+k0],a[jt+k2],a[jp+k2],a[jt+k1],a[jp+k1]);	/*printf("%2u.%x: %2x, k0-2 = %3x,%3x,%3x, a0-5 = %15.10f,%15.10f, %15.10f,%15.10f, %15.10f,%15.10f\n",k/3,0xd,l,(jt+k0)/p1,(jt+k1)/p1,(jt+k2)/p1,a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2]);*/ ++l; ++l1; ++l2;
			jt = j1 + p1; jp = j2 + p1;	RADIX_03_DFT(s,c3m1, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im, t00,t01,t02,t03,t04,t05, a[jt+k2],a[jp+k2],a[jt+k1],a[jp+k1],a[jt+k0],a[jp+k0]);	/*printf("%2u.%x: %2x, k0-2 = %3x,%3x,%3x, a0-5 = %15.10f,%15.10f, %15.10f,%15.10f, %15.10f,%15.10f\n",k/3,0xe,l,(jt+k0)/p1,(jt+k1)/p1,(jt+k2)/p1,a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2]);*/ ++l; l &= 0x07f; l1 = l+128; l2 = l+256;	// <*** needed for final-loop-pass wraparound of 0-term
			jt = j1     ; jp = j2     ;	RADIX_03_DFT(s,c3m1, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im, t00,t01,t02,t03,t04,t05, a[jt+k1],a[jp+k1],a[jt+k0],a[jp+k0],a[jt+k2],a[jp+k2]);	/*printf("%2u.%x: %2x, k0-2 = %3x,%3x,%3x, a0-5 = %15.10f,%15.10f, %15.10f,%15.10f, %15.10f,%15.10f\n",k/3,0xf,l,(jt+k0)/p1,(jt+k1)/p1,(jt+k2)/p1,a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2]);*/ ++l; ++l1; ++l2;
		}
		#endif
	}
}

/******************** Multithreaded function body - NO STATIC VARS BELOW THIS POINT!: ***************************/

#ifdef USE_PTHREAD

	#ifndef COMPILER_TYPE_GCC
		#error pthreaded carry code requires GCC build!
	#endif

	void*
	cy384_process_chunk(void*targ)	// Thread-arg pointer *must* be cast to void and specialized inside the function
	{
	const char func[] = "radix384_ditN_cy_dif1";
		struct cy_thread_data_t* thread_arg = targ;	// Move to top because scalar-mode carry pointers taken directly from it
		double *addr;
		const int pfetch_dist = PFETCH_DIST;
		const int stride = (int)RE_IM_STRIDE << 1;	// main-array loop stride = 2*RE_IM_STRIDE
		uint32 p1,p2,p3,p4,p5,p6,p7,p8,p9,pa,pb,pc,pd,pe,pf,
		     p010,p020,p030,p040,p050,p060,p070,p080,p090,p0a0,p0b0,p0c0,p0d0,p0e0,p0f0,p100,p110,p120,p130,p140,p150,p160,p170;
		int poff[RADIX>>2];
		double wt_re,wt_im, wi_re,wi_im;	// Fermat-mod/LOACC weights stuff, used in both scalar and SIMD mode

		int dif_i_offsets[128], dif_o_offsets[384];
		int dit_i_offsets[384], dit_o_offsets[128];
		int dif_triplets[72];
		int dit_triplets[24];	// Only need 1/3 as many here as for DIF

		int incr,j,j1,j2,k,l,l1,l2,k0,k1,k2;
		// incr = Carry-chain wts-multipliers recurrence length, which must divide
		// RADIX/[n-wayness of carry macro], e.g. RADIX/[16|8|4] = 48|96|192 for avx512,avx,sse, respectively:
		const int incr_long = 16,incr_med = 8,incr_short = 4;
		// Allows cy-macro error data to be used to fiddle incr on the fly to a smaller, safer value if necessary
		if(USE_SHORT_CY_CHAIN == 0)
			incr = incr_long;
		else if(USE_SHORT_CY_CHAIN == 1)
			incr = incr_med;
		else
			incr = incr_short;
	#ifdef USE_AVX512
		double t0,t1,t2,t3;
	  #ifdef CARRY_16_WAY
		struct uint32x16 *n_minus_sil,*n_minus_silp1,*sinwt,*sinwtm1;
	  #else
		struct uint32x8  *n_minus_sil,*n_minus_silp1,*sinwt,*sinwtm1;
	  #endif
	#elif defined(USE_AVX)
		struct uint32x4 *n_minus_sil,*n_minus_silp1,*sinwt,*sinwtm1;
	#else
		int n_minus_sil,n_minus_silp1,sinwt,sinwtm1;
		double wtl,wtlp1,wtn,wtnm1;	/* Mersenne-mod weights stuff */
	#endif

	#ifdef USE_SSE2

		const double crnd = 3.0*0x4000000*0x2000000;
		double *add0,*add1,*add2,*add3;
		int *bjmodn;	// Alloc mem for this along with other 	SIMD stuff
		vec_dbl *tmp,*tm1,*tm2;	// utility ptrs
		int *itmp,*itm2;			// Pointer into the bjmodn array
		struct complex *ctmp;	// Hybrid AVX-DFT/SSE2-carry scheme used for Mersenne-mod needs a 2-word-double pointer
		vec_dbl *two,*one,*sqrt2,*isrt2, *cc0, *ss0, *cc1, *ss1, *max_err, *sse2_rnd, *half_arr,
			// ptrs to 16 sets of twiddles shared by the 2nd-half DIF and DIT DFT macros:
			*twid0,*twid1,*twid2,*twid3,*twid4,*twid5,*twid6,*twid7,*twid8,*twid9,*twida,*twidb,*twidc,*twidd,*twide,*twidf,
			*r000,*r080,*r100,	// Head of RADIX*vec_cmplx-sized local store #1
			*s1p000,	// Head of RADIX*vec_cmplx-sized local store #2
			*cy;	// Need RADIX/2 slots for sse2 carries, RADIX/4 for avx
		double dtmp;
		uint64 *sign_mask, *sse_bw, *sse_sw, *sse_n;

	#else

		const double	c3m1= -1.50000000000000000000,	/* cos(twopi/3)-1	*/
						s   =  0.86602540378443864675;	/* sin(twopi/3)		*/
		double *base, *baseinv;
		int p0123[4];
		const double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
		int jt,jp,m,m2;
		double wt,wtinv,wtA,wtB,wtC;	/* Mersenne-mod weights stuff */
		int bjmodn[RADIX];	// Thread only carries a base datum here, must alloc a local array for remaining values
		double *cy = thread_arg->cy, temp,frac,
			t00,t01,t02,t03,t04,t05;

		// Local storage: We must use an array here because scalars have no guarantees about relative address offsets
		// [and even if those are contiguous-as-hoped-for, they may run in reverse]; Make array type (struct complex)
		// to allow us to use the same offset-indexing as in the original radix-32 in-place DFT macros:
		struct complex t[RADIX], *tptr;
		int *itmp;	// Pointer into the bjmodn array

	#endif

	// int data:
		int tid = thread_arg->tid;
		int iter = thread_arg->iter;
		int NDIVR = thread_arg->ndivr;
		int n = NDIVR*RADIX;
		int target_idx = thread_arg->target_idx;
		int target_set = thread_arg->target_set;
		double target_cy  = thread_arg->target_cy;
		int khi    = thread_arg->khi;
		int i      = thread_arg->i;	/* Pointer to the BASE and BASEINV arrays.	*/
		int jstart = thread_arg->jstart;
		int jhi    = thread_arg->jhi;
		int col = thread_arg->col;
		int co2 = thread_arg->co2;
		int co3 = thread_arg->co3;
		int sw  = thread_arg->sw;
		int nwt = thread_arg->nwt;

	// double data:
		double maxerr = thread_arg->maxerr;
		double scale = thread_arg->scale;	int full_pass = scale < 0.5;
		double prp_mult = thread_arg->prp_mult;

	// pointer data:
		double *a = thread_arg->arrdat;
		double *wt0 = thread_arg->wt0;
		double *wt1 = thread_arg->wt1;
		double *wts_mult = thread_arg->wts_mult;	// Const Intra-block wts-multiplier...
		double *inv_mult = thread_arg->inv_mult;	// ...and 2*(its multiplicative inverse).
		ASSERT(fabs(wts_mult[0]*inv_mult[0] - 1.0) < EPS, "wts_mults fail accuracy check!");
		ASSERT(fabs(wts_mult[1]*inv_mult[1] - 1.0) < EPS, "wts_mults fail accuracy check!");
		int *si = thread_arg->si;

		/*   constant index offsets for array load/stores are here.	*/
		p1 = NDIVR;		p010 = NDIVR<<4;	p100 = NDIVR<<8;
		p2 = p1 + p1;	p020 = p010 + p010;	p110 = p100 + p010;
		p3 = p2 + p1;	p030 = p020 + p010;	p120 = p110 + p010;
		p4 = p3 + p1;	p040 = p030 + p010;	p130 = p120 + p010;
		p5 = p4 + p1;	p050 = p040 + p010;	p140 = p130 + p010;
		p6 = p5 + p1;	p060 = p050 + p010;	p150 = p140 + p010;
		p7 = p6 + p1;	p070 = p060 + p010;	p160 = p150 + p010;
		p8 = p7 + p1;	p080 = p070 + p010;	p170 = p160 + p010;
		p9 = p8 + p1;	p090 = p080 + p010;
		pa = p9 + p1;	p0a0 = p090 + p010;
		pb = pa + p1;	p0b0 = p0a0 + p010;
		pc = pb + p1;	p0c0 = p0b0 + p010;
		pd = pc + p1;	p0d0 = p0c0 + p010;
		pe = pd + p1;	p0e0 = p0d0 + p010;
		pf = pe + p1;	p0f0 = p0e0 + p010;

		p1 += ( (p1 >> DAT_BITS) << PAD_BITS );		p010 += ( (p010 >> DAT_BITS) << PAD_BITS );		p100 += ( (p100 >> DAT_BITS) << PAD_BITS );
		p2 += ( (p2 >> DAT_BITS) << PAD_BITS );		p020 += ( (p020 >> DAT_BITS) << PAD_BITS );		p110 += ( (p110 >> DAT_BITS) << PAD_BITS );
		p3 += ( (p3 >> DAT_BITS) << PAD_BITS );		p030 += ( (p030 >> DAT_BITS) << PAD_BITS );		p120 += ( (p120 >> DAT_BITS) << PAD_BITS );
		p4 += ( (p4 >> DAT_BITS) << PAD_BITS );		p040 += ( (p040 >> DAT_BITS) << PAD_BITS );		p130 += ( (p130 >> DAT_BITS) << PAD_BITS );
		p5 += ( (p5 >> DAT_BITS) << PAD_BITS );		p050 += ( (p050 >> DAT_BITS) << PAD_BITS );		p140 += ( (p140 >> DAT_BITS) << PAD_BITS );
		p6 += ( (p6 >> DAT_BITS) << PAD_BITS );		p060 += ( (p060 >> DAT_BITS) << PAD_BITS );		p150 += ( (p150 >> DAT_BITS) << PAD_BITS );
		p7 += ( (p7 >> DAT_BITS) << PAD_BITS );		p070 += ( (p070 >> DAT_BITS) << PAD_BITS );		p160 += ( (p160 >> DAT_BITS) << PAD_BITS );
		p8 += ( (p8 >> DAT_BITS) << PAD_BITS );		p080 += ( (p080 >> DAT_BITS) << PAD_BITS );		p170 += ( (p170 >> DAT_BITS) << PAD_BITS );
		p9 += ( (p9 >> DAT_BITS) << PAD_BITS );		p090 += ( (p090 >> DAT_BITS) << PAD_BITS );
		pa += ( (pa >> DAT_BITS) << PAD_BITS );		p0a0 += ( (p0a0 >> DAT_BITS) << PAD_BITS );
		pb += ( (pb >> DAT_BITS) << PAD_BITS );		p0b0 += ( (p0b0 >> DAT_BITS) << PAD_BITS );
		pc += ( (pc >> DAT_BITS) << PAD_BITS );		p0c0 += ( (p0c0 >> DAT_BITS) << PAD_BITS );
		pd += ( (pd >> DAT_BITS) << PAD_BITS );		p0d0 += ( (p0d0 >> DAT_BITS) << PAD_BITS );
		pe += ( (pe >> DAT_BITS) << PAD_BITS );		p0e0 += ( (p0e0 >> DAT_BITS) << PAD_BITS );
		pf += ( (pf >> DAT_BITS) << PAD_BITS );		p0f0 += ( (p0f0 >> DAT_BITS) << PAD_BITS );
	#ifndef USE_SSE2
		p0123[0] = 0; p0123[1] = p1; p0123[2] = p2; p0123[3] = p3;
	#endif
		poff[     0] =    0; poff[     1] =      p4; poff[     2] =      p8; poff[     3] =      pc;
		poff[0x04+0] = p010; poff[0x04+1] = p010+p4; poff[0x04+2] = p010+p8; poff[0x04+3] = p010+pc;
		poff[0x08+0] = p020; poff[0x08+1] = p020+p4; poff[0x08+2] = p020+p8; poff[0x08+3] = p020+pc;
		poff[0x0c+0] = p030; poff[0x0c+1] = p030+p4; poff[0x0c+2] = p030+p8; poff[0x0c+3] = p030+pc;
		poff[0x10+0] = p040; poff[0x10+1] = p040+p4; poff[0x10+2] = p040+p8; poff[0x10+3] = p040+pc;
		poff[0x14+0] = p050; poff[0x14+1] = p050+p4; poff[0x14+2] = p050+p8; poff[0x14+3] = p050+pc;
		poff[0x18+0] = p060; poff[0x18+1] = p060+p4; poff[0x18+2] = p060+p8; poff[0x18+3] = p060+pc;
		poff[0x1c+0] = p070; poff[0x1c+1] = p070+p4; poff[0x1c+2] = p070+p8; poff[0x1c+3] = p070+pc;
		poff[0x20+0] = p080; poff[0x20+1] = p080+p4; poff[0x20+2] = p080+p8; poff[0x20+3] = p080+pc;
		poff[0x24+0] = p090; poff[0x24+1] = p090+p4; poff[0x24+2] = p090+p8; poff[0x24+3] = p090+pc;
		poff[0x28+0] = p0a0; poff[0x28+1] = p0a0+p4; poff[0x28+2] = p0a0+p8; poff[0x28+3] = p0a0+pc;
		poff[0x2c+0] = p0b0; poff[0x2c+1] = p0b0+p4; poff[0x2c+2] = p0b0+p8; poff[0x2c+3] = p0b0+pc;
		poff[0x30+0] = p0c0; poff[0x30+1] = p0c0+p4; poff[0x30+2] = p0c0+p8; poff[0x30+3] = p0c0+pc;
		poff[0x34+0] = p0d0; poff[0x34+1] = p0d0+p4; poff[0x34+2] = p0d0+p8; poff[0x34+3] = p0d0+pc;
		poff[0x38+0] = p0e0; poff[0x38+1] = p0e0+p4; poff[0x38+2] = p0e0+p8; poff[0x38+3] = p0e0+pc;
		poff[0x3c+0] = p0f0; poff[0x3c+1] = p0f0+p4; poff[0x3c+2] = p0f0+p8; poff[0x3c+3] = p0f0+pc;
		for(l = 0; l < 32; l++) {
			poff[ 64+l] = poff[l] + p100;
		}

	/*** DIF: ***/
	// Set array offsets for radix-128 inputs [these are same for all 3 such DFTs].
		for(j = 0; j < 0x80; ++j) { dif_i_offsets[j] = j<<1; }
	// For the radix-128 DIF outputs we need the following offsets:
		const int dif_perm0[16] = {0,p1,p2,p3,p5,p4,p7,p6,pa,pb,p9,p8,pf,pe,pc,pd};
		const int dif_perm1[16] = {p5,p4,p7,p6,p2,p3,p1,0,pf,pe,pc,pd,p9,p8,pb,pa};
		const int dif_perm2[16] = {pa,pb,p9,p8,pf,pe,pc,pd,p5,p4,p7,p6,p2,p3,p1,0};
		const int dif_perm3[16] = {pf,pe,pc,pd,p9,p8,pb,pa,p2,p3,p1,0,p7,p6,p4,p5};
		for(j = 0; j < 0x10; ++j) {
			dif_o_offsets[j] = dif_perm0[j];
			dif_o_offsets[j+0x010] = p010 + dif_perm1[j];
			dif_o_offsets[j+0x020] = p020 + dif_perm2[j];
			dif_o_offsets[j+0x030] = p030 + dif_perm3[j];
			dif_o_offsets[j+0x040] = p050 + dif_perm1[j];
			dif_o_offsets[j+0x050] = p040 + dif_perm2[j];
			dif_o_offsets[j+0x060] = p070 + dif_perm3[j];
			dif_o_offsets[j+0x070] = p060 + dif_perm1[j];
		// dif_o_offsets are w.r.to each respective block of 128 complex data, i.e. p-indexing is (mod p80):
			dif_o_offsets[j+0x080] = p020 + dif_perm2[j];
			dif_o_offsets[j+0x090] = p030 + dif_perm3[j];
			dif_o_offsets[j+0x0a0] = p010 + dif_perm1[j];
			dif_o_offsets[j+0x0b0] =    0 + dif_perm2[j];
			dif_o_offsets[j+0x0c0] = p070 + dif_perm3[j];
			dif_o_offsets[j+0x0d0] = p060 + dif_perm1[j];
			dif_o_offsets[j+0x0e0] = p040 + dif_perm2[j];
			dif_o_offsets[j+0x0f0] = p050 + dif_perm3[j];
		// dif_o_offsets are w.r.to each respective block of 128 complex data, i.e. p-indexing is (mod p80):
			dif_o_offsets[j+0x100] = p050 + dif_perm1[j];
			dif_o_offsets[j+0x110] = p040 + dif_perm2[j];
			dif_o_offsets[j+0x120] = p070 + dif_perm3[j];
			dif_o_offsets[j+0x130] = p060 + dif_perm1[j];
			dif_o_offsets[j+0x140] = p020 + dif_perm2[j];
			dif_o_offsets[j+0x150] = p030 + dif_perm3[j];
			dif_o_offsets[j+0x160] = p010 + dif_perm1[j];
			dif_o_offsets[j+0x170] =    0 + dif_perm2[j];
		}
	/*** DIT: ***/
	// Set array offsets for radix-128 outputs [these are same for all 3 such DFTs].
		for(j = 0; j < 0x80; ++j) { dit_o_offsets[j] = j<<1; }
	// For the radix-128 DIT inputs we need the following offsets:
		const int dit_perm0[16] = {0,p1,p3,p2,p7,p6,p5,p4,pf,pe,pd,pc,pb,pa,p9,p8};
		const int dit_perm1[16] = {pf,pe,pd,pc,pb,pa,p9,p8,p7,p6,p5,p4,p3,p2,p1,0};
		const int dit_perm2[16] = {p5,p4,p6,p7,p1,0,p2,p3,p9,p8,pa,pb,pe,pf,pc,pd};
		const int dit_perm3[16] = {p9,p8,pa,pb,pe,pf,pc,pd,p1,0,p2,p3,p6,p7,p4,p5};
		const int dit_perm4[16] = {pa,pb,p8,p9,pc,pd,pf,pe,p2,p3,0,p1,p4,p5,p7,p6};
		const int dit_perm5[16] = {p2,p3,0,p1,p4,p5,p7,p6,pc,pd,pf,pe,p8,p9,pb,pa};
		for(j = 0; j < 0x10; ++j) {
			dit_i_offsets[j] = dit_perm0[j];
			dit_i_offsets[j+0x010] = p010 + dit_perm1[j];
			dit_i_offsets[j+0x020] = p030 + dit_perm1[j];
			dit_i_offsets[j+0x030] = p020 + dit_perm1[j];
			dit_i_offsets[j+0x040] = p070 + dit_perm1[j];
			dit_i_offsets[j+0x050] = p060 + dit_perm1[j];
			dit_i_offsets[j+0x060] = p050 + dit_perm1[j];
			dit_i_offsets[j+0x070] = p040 + dit_perm1[j];
		// dit_i_offsets are w.r.to each respective block of 128 complex data, i.e. p-indexing is (mod p80):
			dit_i_offsets[j+0x080] = p050 + dit_perm2[j];
			dit_i_offsets[j+0x090] = p040 + dit_perm2[j];
			dit_i_offsets[j+0x0a0] = p060 + dit_perm2[j];
			dit_i_offsets[j+0x0b0] = p070 + dit_perm3[j];
			dit_i_offsets[j+0x0c0] = p010 + dit_perm2[j];
			dit_i_offsets[j+0x0d0] =    0 + dit_perm2[j];
			dit_i_offsets[j+0x0e0] = p020 + dit_perm2[j];
			dit_i_offsets[j+0x0f0] = p030 + dit_perm3[j];
		// dit_i_offsets are w.r.to each respective block of 128 complex data, i.e. p-indexing is (mod p80):
			dit_i_offsets[j+0x100] = p020 + dit_perm4[j];
			dit_i_offsets[j+0x110] = p030 + dit_perm5[j];
			dit_i_offsets[j+0x120] =    0 + dit_perm4[j];
			dit_i_offsets[j+0x130] = p010 + dit_perm5[j];
			dit_i_offsets[j+0x140] = p040 + dit_perm4[j];
			dit_i_offsets[j+0x150] = p050 + dit_perm5[j];
			dit_i_offsets[j+0x160] = p070 + dit_perm5[j];
			dit_i_offsets[j+0x170] = p060 + dit_perm5[j];
		}
	// DIF Index-high-bits triplets needed for compact-obj-code scheme:
	  #ifdef USE_SSE2
		k = 0;	// In SIMD mode these are 0xtr-offsets w.r.to a local store:
		dif_triplets[k] = 0x170; dif_triplets[k+1] = 0x0f0; dif_triplets[k+2] = 0x070; k += 3;
		dif_triplets[k] = 0x160; dif_triplets[k+1] = 0x0e0; dif_triplets[k+2] = 0x060; k += 3;
		dif_triplets[k] = 0x150; dif_triplets[k+1] = 0x0d0; dif_triplets[k+2] = 0x050; k += 3;
		dif_triplets[k] = 0x140; dif_triplets[k+1] = 0x0c0; dif_triplets[k+2] = 0x040; k += 3;
		dif_triplets[k] = 0x130; dif_triplets[k+1] = 0x0b0; dif_triplets[k+2] = 0x030; k += 3;
		dif_triplets[k] = 0x120; dif_triplets[k+1] = 0x0a0; dif_triplets[k+2] = 0x020; k += 3;
		dif_triplets[k] = 0x110; dif_triplets[k+1] = 0x090; dif_triplets[k+2] = 0x010; k += 3;
		dif_triplets[k] = 0x100; dif_triplets[k+1] = 0x080; dif_triplets[k+2] =     0; k += 3;
		dif_triplets[k] = 0x0f0; dif_triplets[k+1] = 0x070; dif_triplets[k+2] = 0x170; k += 3;
		dif_triplets[k] = 0x0e0; dif_triplets[k+1] = 0x060; dif_triplets[k+2] = 0x160; k += 3;
		dif_triplets[k] = 0x0d0; dif_triplets[k+1] = 0x050; dif_triplets[k+2] = 0x150; k += 3;
		dif_triplets[k] = 0x0c0; dif_triplets[k+1] = 0x040; dif_triplets[k+2] = 0x140; k += 3;
		dif_triplets[k] = 0x0b0; dif_triplets[k+1] = 0x030; dif_triplets[k+2] = 0x130; k += 3;
		dif_triplets[k] = 0x0a0; dif_triplets[k+1] = 0x020; dif_triplets[k+2] = 0x120; k += 3;
		dif_triplets[k] = 0x090; dif_triplets[k+1] = 0x010; dif_triplets[k+2] = 0x110; k += 3;
		dif_triplets[k] = 0x080; dif_triplets[k+1] =     0; dif_triplets[k+2] = 0x100; k += 3;
		dif_triplets[k] = 0x070; dif_triplets[k+1] = 0x170; dif_triplets[k+2] = 0x0f0; k += 3;
		dif_triplets[k] = 0x060; dif_triplets[k+1] = 0x160; dif_triplets[k+2] = 0x0e0; k += 3;
		dif_triplets[k] = 0x050; dif_triplets[k+1] = 0x150; dif_triplets[k+2] = 0x0d0; k += 3;
		dif_triplets[k] = 0x040; dif_triplets[k+1] = 0x140; dif_triplets[k+2] = 0x0c0; k += 3;
		dif_triplets[k] = 0x030; dif_triplets[k+1] = 0x130; dif_triplets[k+2] = 0x0b0; k += 3;
		dif_triplets[k] = 0x020; dif_triplets[k+1] = 0x120; dif_triplets[k+2] = 0x0a0; k += 3;
		dif_triplets[k] = 0x010; dif_triplets[k+1] = 0x110; dif_triplets[k+2] = 0x090; k += 3;
		dif_triplets[k] =     0; dif_triplets[k+1] = 0x100; dif_triplets[k+2] = 0x080;
		// IN SIMD mode need to double all the above to turn from vec_dbl to vec_cmplx ptr offsets:
		for(l = 0; l < 72; l++) {
			dif_triplets[l] <<= 1;
		}
	  #else
	// Index-high-bits triplets needed for compact-obj-code scheme:
		k = 0;
		dif_triplets[k] = p170; dif_triplets[k+1] = p0f0; dif_triplets[k+2] = p070; k += 3;
		dif_triplets[k] = p160; dif_triplets[k+1] = p0e0; dif_triplets[k+2] = p060; k += 3;
		dif_triplets[k] = p150; dif_triplets[k+1] = p0d0; dif_triplets[k+2] = p050; k += 3;
		dif_triplets[k] = p140; dif_triplets[k+1] = p0c0; dif_triplets[k+2] = p040; k += 3;
		dif_triplets[k] = p130; dif_triplets[k+1] = p0b0; dif_triplets[k+2] = p030; k += 3;
		dif_triplets[k] = p120; dif_triplets[k+1] = p0a0; dif_triplets[k+2] = p020; k += 3;
		dif_triplets[k] = p110; dif_triplets[k+1] = p090; dif_triplets[k+2] = p010; k += 3;
		dif_triplets[k] = p100; dif_triplets[k+1] = p080; dif_triplets[k+2] =    0; k += 3;
		dif_triplets[k] = p0f0; dif_triplets[k+1] = p070; dif_triplets[k+2] = p170; k += 3;
		dif_triplets[k] = p0e0; dif_triplets[k+1] = p060; dif_triplets[k+2] = p160; k += 3;
		dif_triplets[k] = p0d0; dif_triplets[k+1] = p050; dif_triplets[k+2] = p150; k += 3;
		dif_triplets[k] = p0c0; dif_triplets[k+1] = p040; dif_triplets[k+2] = p140; k += 3;
		dif_triplets[k] = p0b0; dif_triplets[k+1] = p030; dif_triplets[k+2] = p130; k += 3;
		dif_triplets[k] = p0a0; dif_triplets[k+1] = p020; dif_triplets[k+2] = p120; k += 3;
		dif_triplets[k] = p090; dif_triplets[k+1] = p010; dif_triplets[k+2] = p110; k += 3;
		dif_triplets[k] = p080; dif_triplets[k+1] =    0; dif_triplets[k+2] = p100; k += 3;
		dif_triplets[k] = p070; dif_triplets[k+1] = p170; dif_triplets[k+2] = p0f0; k += 3;
		dif_triplets[k] = p060; dif_triplets[k+1] = p160; dif_triplets[k+2] = p0e0; k += 3;
		dif_triplets[k] = p050; dif_triplets[k+1] = p150; dif_triplets[k+2] = p0d0; k += 3;
		dif_triplets[k] = p040; dif_triplets[k+1] = p140; dif_triplets[k+2] = p0c0; k += 3;
		dif_triplets[k] = p030; dif_triplets[k+1] = p130; dif_triplets[k+2] = p0b0; k += 3;
		dif_triplets[k] = p020; dif_triplets[k+1] = p120; dif_triplets[k+2] = p0a0; k += 3;
		dif_triplets[k] = p010; dif_triplets[k+1] = p110; dif_triplets[k+2] = p090; k += 3;
		dif_triplets[k] =    0; dif_triplets[k+1] = p100; dif_triplets[k+2] = p080;
	  #endif
	// DIT Index-high-bits triplets needed for compact-obj-code scheme:
	  #ifdef USE_SSE2
		k = 0;	// In SIMD mode these are ptr-offsets w.r.to a local store:
		dit_triplets[k] = 0x070; dit_triplets[k+1] = 0x0f0; dit_triplets[k+2] = 0x170; k += 3;
		dit_triplets[k] = 0x160; dit_triplets[k+1] = 0x060; dit_triplets[k+2] = 0x0e0; k += 3;
		dit_triplets[k] = 0x0d0; dit_triplets[k+1] = 0x150; dit_triplets[k+2] = 0x050; k += 3;
		dit_triplets[k] = 0x040; dit_triplets[k+1] = 0x0c0; dit_triplets[k+2] = 0x140; k += 3;
		dit_triplets[k] = 0x130; dit_triplets[k+1] = 0x030; dit_triplets[k+2] = 0x0b0; k += 3;
		dit_triplets[k] = 0x0a0; dit_triplets[k+1] = 0x120; dit_triplets[k+2] = 0x020; k += 3;
		dit_triplets[k] = 0x010; dit_triplets[k+1] = 0x090; dit_triplets[k+2] = 0x110; k += 3;
		dit_triplets[k] = 0x100; dit_triplets[k+1] =     0; dit_triplets[k+2] = 0x080;
		// IN SIMD mode need to double all the above to turn from vec_dbl to vec_cmplx ptr offsets:
		for(l = 0; l < 24; l++) {
			dit_triplets[l] <<= 1;
		}
	  #else
		// Cf. comments in radix384_dit_pass1 re. this needed reordering:
		k = 0;
		dit_triplets[k] = p070; dit_triplets[k+1] = p0f0; dit_triplets[k+2] = p170; k += 3;
		dit_triplets[k] = p160; dit_triplets[k+1] = p060; dit_triplets[k+2] = p0e0; k += 3;
		dit_triplets[k] = p0d0; dit_triplets[k+1] = p150; dit_triplets[k+2] = p050; k += 3;
		dit_triplets[k] = p040; dit_triplets[k+1] = p0c0; dit_triplets[k+2] = p140; k += 3;
		dit_triplets[k] = p130; dit_triplets[k+1] = p030; dit_triplets[k+2] = p0b0; k += 3;
		dit_triplets[k] = p0a0; dit_triplets[k+1] = p120; dit_triplets[k+2] = p020; k += 3;
		dit_triplets[k] = p010; dit_triplets[k+1] = p090; dit_triplets[k+2] = p110; k += 3;
		dit_triplets[k] = p100; dit_triplets[k+1] =    0; dit_triplets[k+2] = p080;
	  #endif

	#ifdef USE_SSE2
		tmp = r000 = thread_arg->r000;
		r080 = tmp + 0x100;
		r100 = tmp + 0x200;
		s1p000 = tmp + 0x300;
		tmp += 0x600;	// sc_ptr += 0x600
		two     = tmp + 0;	// AVX+ versions of radix-8,16,32 twiddleless-DFT macros need consts [2,1,sqrt2,isrt2] quartet laid out thusly
		one     = tmp + 1;
		sqrt2   = tmp + 2;
		isrt2   = tmp + 3;
		cc0		= tmp + 4;
		ss0		= tmp + 5;
		cc1		= tmp + 6;
		ss1		= tmp + 7;
		tmp += 0x08;	// sc_ptr += 0x608
		// Each non-unity root now needs a negated counterpart:
		/* Stupidity: Since a truly general-purpose [in the sense that it can be used for our radix-128 internal-twiddles]
		radix-8 DFT-with-twiddles macro needs 8 in-addresses [corr. to the 8 real parts of the input data], 8 o-addresses,
		and 7 each of cosine and sine data [which cannot be assumed to occur in fixed-stride pairs - cf. our usage of
		SSE2_RADIX8_DIT_TWIDDLE_OOP() below], that hits the GCC hard limit of 30-operands for ASM macros, but we still
		need one more operand for the ISRT2 pointer. Only easy workaround I found for this is to stick a vector-ISRT2 copy
		in between each +-[cc,ss] vector-data pair, thus any time we need a vector-isrt2 for the radix-8 internal twiddles
		we get it at (vec_dbl*)cc-1.	/Stupidity */
		// ptrs to 16 sets (3*7 = 21 vec_dbl data each) of non-unity twiddles shared by the 2nd-half DIF and DIT DFT macros:
		twid0  = tmp + 0x00;
		twid1  = tmp + 0x15;
		twid2  = tmp + 0x2a;
		twid3  = tmp + 0x3f;
		twid4  = tmp + 0x54;
		twid5  = tmp + 0x69;
		twid6  = tmp + 0x7e;
		twid7  = tmp + 0x93;
		twid8  = tmp + 0xa8;
		twid9  = tmp + 0xbd;
		twida  = tmp + 0xd2;
		twidb  = tmp + 0xe7;
		twidc  = tmp + 0xfc;
		twidd  = tmp + 0x111;
		twide  = tmp + 0x126;
		twidf  = tmp + 0x13b;
		tmp += 0x150;	// += 16*21 = 0x150 => sc_ptr + 0x758
	  #ifdef USE_AVX512
		cy = tmp;		tmp += 0x30;	// RADIX/8 vec_dbl slots for carry sub-array
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;
		half_arr= tmp + 0x02;	// += 0x30 + 2 => sc_ptr += 0x78a
	  #elif defined(USE_AVX)
		cy = tmp;		tmp += 0x60;	// RADIX/4 vec_dbl slots for carry sub-array
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// += 0x60 + 2 => sc_ptr += 0x7ba
		// This is where the value of half_arr_offset comes from
		half_arr= tmp + 0x02;	/* This table needs 96 vec_dbl for Mersenne-mod, and 3.5*RADIX[avx] | RADIX[sse2] for Fermat-mod */
	  #else
		cy = tmp;		tmp += 0xc0;	// RADIX/2 vec_dbl slots for carry sub-array
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// += 0xc0 + 2 => sc_ptr += 0x81a
		// This is where the value of half_arr_offset comes from
		half_arr= tmp + 0x02;	/* This table needs 32 x 16 bytes for Mersenne-mod, 2 for Fermat-mod */
	  #endif

		ASSERT((r000 == thread_arg->r000), "thread-local memcheck failed!");
		ASSERT((half_arr == thread_arg->half_arr), "thread-local memcheck failed!");
	  #ifndef USE_AVX512	// In AVX-512 mode, use VRNDSCALEPD for rounding and hijack this vector-data slot for the 4 base/baseinv-consts:
		ASSERT((sse2_rnd->d0 == crnd && sse2_rnd->d1 == crnd), "thread-local memcheck failed!");
	  #endif
		tmp = half_arr;
	  #ifdef USE_AVX512
		/* No-Op */
	  #elif defined(USE_AVX)
		// Grab some elt of base-data [offset by, say, +32] and mpy by its inverse [+16 further]
		dtmp = (tmp+40)->d0 * (tmp+56)->d0;	ASSERT(fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
		dtmp = (tmp+40)->d1 * (tmp+56)->d1;	ASSERT(fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
	  #else	// SSE2:
		dtmp = (tmp+10)->d0 * (tmp+14)->d0;	ASSERT(fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
		dtmp = (tmp+10)->d1 * (tmp+14)->d1;	ASSERT(fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
	  #endif

		VEC_DBL_INIT(max_err, 0.0);

		sign_mask = (uint64*)(r000 + radix384_creals_in_local_store);
		sse_bw  = sign_mask + RE_IM_STRIDE;	// (  #doubles in a SIMD complex) x 32-bits = RE_IM_STRIDE x 64-bits
		sse_sw  = sse_bw    + RE_IM_STRIDE;
		sse_n   = sse_sw    + RE_IM_STRIDE;

	  #ifdef USE_AVX512
	   #ifdef CARRY_16_WAY
		n_minus_sil   = (struct uint32x16*)sse_n + 1;
		n_minus_silp1 = (struct uint32x16*)sse_n + 2;
		sinwt         = (struct uint32x16*)sse_n + 3;
		sinwtm1       = (struct uint32x16*)sse_n + 4;
	   #else
		n_minus_sil   = (struct uint32x8 *)sse_n + 1;
		n_minus_silp1 = (struct uint32x8 *)sse_n + 2;
		sinwt         = (struct uint32x8 *)sse_n + 3;
		sinwtm1       = (struct uint32x8 *)sse_n + 4;
	   #endif
	  #elif defined(USE_AVX)
		n_minus_sil   = (struct uint32x4 *)sse_n + 1;
		n_minus_silp1 = (struct uint32x4 *)sse_n + 2;
		sinwt         = (struct uint32x4 *)sse_n + 3;
		sinwtm1       = (struct uint32x4 *)sse_n + 4;
	  #endif
	  #ifdef USE_AVX
		bjmodn = (int*)(sinwtm1 + RE_IM_STRIDE);
	  #else
		bjmodn = (int*)(sse_n + RE_IM_STRIDE);
	  #endif

	#else

		// In scalar mode use these 2 ptrs to pass the base & baseinv arrays:
		base    = (double *)thread_arg->r000    ;
		baseinv = (double *)thread_arg->half_arr;

	#endif	// USE_SSE2 ?

		/* Init DWT-indices: */
		uint32 bjmodnini = thread_arg->bjmodnini;
		bjmodn[0] = thread_arg->bjmodn0;
		for(l = 1; l < RADIX; l++) {	// must use e.g. l for loop idx here as i is used for dwt indexing
			MOD_ADD32(bjmodn[l-1], bjmodnini, n, bjmodn[l]);
		}

		/* init carries	*/
		addr = thread_arg->cy;
	#ifdef USE_AVX512
		tmp = cy;
		for(l = 0; l < RADIX; l += 8, ++tmp) {
			tmp->d0 = *(addr+l  );
			tmp->d1 = *(addr+l+1);
			tmp->d2 = *(addr+l+2);
			tmp->d3 = *(addr+l+3);
			tmp->d4 = *(addr+l+4);
			tmp->d5 = *(addr+l+5);
			tmp->d6 = *(addr+l+6);
			tmp->d7 = *(addr+l+7);
		}
	#elif defined(USE_AVX)
		tmp = cy;
		for(l = 0; l < RADIX; l += 4, ++tmp) {
			tmp->d0 = *(addr+l  );
			tmp->d1 = *(addr+l+1);
			tmp->d2 = *(addr+l+2);
			tmp->d3 = *(addr+l+3);
		}
	#elif defined(USE_SSE2)
		tmp = cy;
		for(l = 0; l < RADIX; l += 2, ++tmp) {
			tmp->d0 = *(addr+l  );
			tmp->d1 = *(addr+l+1);
		}
	#elif 0	// No_op in scalar case, since carry pattern matches that of thread data
		for(l = 0; l < RADIX; l++) {
			cy[l] = *(addr+l);
		}
	#endif

		/********************************************************************************/
		/* This main loop is same for un-and-multithreaded, so stick into a header file */
		/* (can't use a macro because of the #if-enclosed stuff).                       */
		/********************************************************************************/
		#include "radix384_main_carry_loop.h"

		/* At end of each thread-processed work chunk, dump the
		carryouts into their non-thread-private array slots:
		*/
		addr = thread_arg->cy;
	#ifdef USE_AVX512
		tmp = cy;
		for(l = 0; l < RADIX; l += 8, ++tmp) {
			*(addr+l  ) = tmp->d0;
			*(addr+l+1) = tmp->d1;
			*(addr+l+2) = tmp->d2;
			*(addr+l+3) = tmp->d3;
			*(addr+l+4) = tmp->d4;
			*(addr+l+5) = tmp->d5;
			*(addr+l+6) = tmp->d6;
			*(addr+l+7) = tmp->d7;
		}
		t0 = MAX(max_err->d0,max_err->d1);
		t1 = MAX(max_err->d2,max_err->d3);
		t2 = MAX(max_err->d4,max_err->d5);
		t3 = MAX(max_err->d6,max_err->d7);
		maxerr = MAX( MAX(t0,t1), MAX(t2,t3) );
	#elif defined(USE_AVX)
		tmp = cy;
		for(l = 0; l < RADIX; l += 4, ++tmp) {
			*(addr+l  ) = tmp->d0;
			*(addr+l+1) = tmp->d1;
			*(addr+l+2) = tmp->d2;
			*(addr+l+3) = tmp->d3;
		}
		maxerr = MAX( MAX(max_err->d0,max_err->d1) , MAX(max_err->d2,max_err->d3) );
	#elif defined(USE_SSE2)
		tmp = cy;
		for(l = 0; l < RADIX; l += 2, ++tmp) {
			*(addr+l  ) = tmp->d0;
			*(addr+l+1) = tmp->d1;
		}
		maxerr = MAX(max_err->d0,max_err->d1);
	#elif 0	// No_op in scalar case, since carry pattern matches that of thread data
		for(l = 0; l < RADIX; l++) {
			*(addr+l) = cy[l];
		}
	#endif

		/* Since will lose separate maxerr values when threads are merged, save them after each pass. */
		if(thread_arg->maxerr < maxerr)
		{
			thread_arg->maxerr = maxerr;
		}

		return 0x0;
	}
#endif

#undef RADIX
#undef PFETCH_DIST
