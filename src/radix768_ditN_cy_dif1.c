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
#include "radix256.h"

#define RADIX 768	// Use #define rather than const int to ensure it's really a compile-time const in the C sense

#define EPS 1e-10

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

  // For Mersenne-mod need (16 [SSE2] or 64 [AVX]) + (4 [HIACC] or 40 [LOACC]) added slots for half_arr lookup tables.
  // Max = (40 [SSE2]; 132 [AVX]), add to (half_arr_offset768 + RADIX) to get SIMD value of radix768_creals_in_local_store:
  #ifdef USE_AVX512	// RADIX/8 = 0x60 fewer carry slots than AVX:
	const int half_arr_offset768 = 0xe47;	// + RADIX + 132 (=0x84); Used for thread local-storage-integrity checking
	const int radix768_creals_in_local_store = 0x11cc;	// ... + 0x300 + 0x084, and round up to nearest multiple of 4
  #elif defined(USE_AVX)
	const int half_arr_offset768 = 0xea7;	// + RADIX + 132 (=0x84); Used for thread local-storage-integrity checking
	const int radix768_creals_in_local_store = 0x122c;	// ... + 0x300 + 0x084, and round up to nearest multiple of 4
  #else
	const int half_arr_offset768 = 0xf67;	// + RADIX + 40 = 0xf67 + 0x300 + 0x028; Used for thread local-storage-integrity checking
	const int radix768_creals_in_local_store = 0x1290;	// ...and round up to nearest multiple of 4
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

int radix768_ditN_cy_dif1(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], double base[], double baseinv[], int iter, double *fracmax, uint64 p)
{
/*
!...Acronym: DWT = Discrete Weighted Transform, DIT = Decimation In Time, DIF = Decimation In Frequency
!
!...Performs a final radix-768 complex DIT pass, an inverse DWT weighting, a carry propagation,
!   a forward DWT weighting, and an initial radix-768 complex DIF pass on the data in the length-N real vector A.
!
!   Data enter and are returned in the A-array.
!
!   See the documentation in mers_mod_square and radix16_dif_pass for further details on the array
!   storage scheme, and radix16_ditN_cy_dif1 for details on the reduced-length weights array scheme.
*/
	const char func[] = "radix768_ditN_cy_dif1";
	const int stride = (int)RE_IM_STRIDE << 1;	// main-array loop stride = 2*RE_IM_STRIDE
#ifndef MULTITHREAD
  #ifdef USE_SSE2
	const int pfetch_dist = PFETCH_DIST;
  #endif

	static int dif_offsets_lo[64];	// 4 subsets of 16
	// Bitfields encoding the sequence of the dif_offsets_lo subset0-3 vectors to use for each radix-256 DIF's outputs:
	const uint32 dif_idx1 = 0xe79e79e4,dif_idx2 = 0x79e79e79,dif_idx3 = 0x9e79e79e;
	static int dif_offsets_hi1[16],dif_offsets_hi2[16],dif_offsets_hi3[16];

	static int dit_offsets_lo[96];	// 6 subsets of 16
	// Bitfields encoding the sequence of the dit_offsets_lo subset0-5 vectors to use for each radix-256 DIT's outputs:
	const uint32 dit_idx1 = 0x55555554,dit_idx2 = 0x44404040,dit_idx3 = 0x54445444;
	static int dit_offsets_hi1[16],dit_offsets_hi2[16],dit_offsets_hi3[16];

	static int dif_triplets[144], dit_triplets[48];
#endif

	static double wts_mult[2], inv_mult[2];	// Const wts-multiplier and 2*(its multiplicative inverse)
  #if !defined(MULTITHREAD) && !defined(USE_SSE2)
	double wt_re,wt_im, wi_re,wi_im;	// Fermat-mod/LOACC weights stuff, used in both scalar and SIMD mode
  #endif
	// Cleanup loop assumes carryins propagate at most 4 words up, but need at least 1 vec_cmplx
	// (2 vec_dbl)'s worth of doubles in wraparound step, hence AVX-512 needs value bumped up:
  #ifdef USE_AVX512
	const int jhi_wrap = 15;
  #else
	const int jhi_wrap =  7;
  #endif
	int NDIVR,i,j,j1,jt,jhi,full_pass,khi,l,outer;
  #if !defined(MULTITHREAD) && !defined(USE_SSE2)
	int j2,jp;
  #endif
  #ifndef MULTITHREAD
	int jstart,k,l1,l2;
  #endif
  #ifdef USE_SSE2
	int nbytes;
  #endif
#if !defined(MULTITHREAD) && defined(USE_SSE2)
	// incr = Carry-chain wts-multipliers recurrence length, which must divide
	// RADIX/[n-wayness of carry macro], e.g. RADIX/[16|8|4] = 48|96|192 for avx512,avx,sse, respectively:
	int incr;
	const int incr_long = 16,incr_med = 8,incr_short = 4;
  // Have no specialized HIACC carry macro in USE_AVX512 and ARMv8 SIMD, so use "goes to 11" in LOACC mode via an incr_hiacc = 2:
  #if defined(USE_AVX512) || defined(USE_ARM_V8_SIMD)
	const int incr_hiacc = 2;
  #else
	const int incr_hiacc = 0;
  #endif
	// Allows cy-macro error data to be used to fiddle incr on the fly to a smaller, safer value if necessary
	if(USE_SHORT_CY_CHAIN == 0)
		incr = incr_long;
	else if(USE_SHORT_CY_CHAIN == 1)
		incr = incr_med;
	else if(USE_SHORT_CY_CHAIN == 2)
		incr = incr_short;
	else
		incr = incr_hiacc;
#endif // !MULTITHREAD && USE_SSE2

#ifndef MULTITHREAD
	int k0,k1,k2;
#endif
	// Jun 2018: Add support for residue shift. (Only LL-test needs intervention at carry-loop level).
	int target_idx = -1, target_set = 0,tidx_mod_stride;
	double target_cy = 0;
	static double ndivr_inv;
	uint64 itmp64;
	static uint64 psave = 0;
	static uint32 bw,sw,bjmodnini,
		p1,p2,p3,p4,p5,p6,p7,p8,p9,pa,pb,pc,pd,pe,pf,
		p10,p20,p30,p40,p50,p60,p70,p80,p90,pa0,pb0,pc0,pd0,pe0,pf0,
		p100,p110,p120,p130,p140,p150,p160,p170,p180,p190,p1a0,p1b0,p1c0,p1d0,p1e0,p1f0,
		p200,p210,p220,p230,p240,p250,p260,p270,p280,p290,p2a0,p2b0,p2c0,p2d0,p2e0,p2f0, nsave = 0;
	static int poff[RADIX>>2];	// Store mults of p4 offset for loop control
	static double radix_inv, n2inv;
#if defined(USE_SSE2) || !defined(MULTITHREAD)
	const double	c3m1= -1.50000000000000000000,	/* cos(twopi/3)-1	*/
					s   =  0.86602540378443864675;
#endif
	// FMA-based DFT needs the tangent:
#if defined(USE_AVX2) && !defined(USE_IMCI512)
	static double tan = 0.41421356237309504879;
#endif
	double scale, dtmp, maxerr = 0.0;
	// Local storage: We must use an array here because scalars have no guarantees about relative address offsets
	// [and even if those are contiguous-as-hoped-for, they may run in reverse]; Make array type (struct complex)
	// to allow us to use the same offset-indexing as in the original radix-32 in-place DFT macros:
	struct complex t[RADIX];
	int err;
	static int first_entry=TRUE;

/*...stuff for the reduced-length DWT weights array is here:	*/
	int n_div_nwt;
#ifndef MULTITHREAD
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
	int col,co2,co3;
  #ifndef USE_AVX
	double wtl,wtlp1,wtn,wtnm1;	/* Mersenne-mod weights stuff */
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
	double *addr, *add0,*add1,*add2,*add3;
	int *itmp;			// Pointer into the bjmodn array
   #if defined(USE_AVX) && !defined(USE_AVX512)
	int *itm2;			// Pointer into the bjmodn array
   #endif
  #endif

	// Uint64 bitmaps for alternate "rounded the other way" copies of sqrt2,isrt2. Default round-to-nearest versions
	// (SQRT2, ISRT2) end in ...3BCD. Since we round these down as ...3BCC90... --> ..3BCC, append _dn to varnames:
	const uint64 sqrt2_dn = 0x3FF6A09E667F3BCCull, isrt2_dn = 0x3FE6A09E667F3BCCull;
  #ifndef MULTITHREAD
	static int *bjmodn;	// Alloc mem for this along with other 	SIMD stuff
  #endif
  #ifndef USE_AVX512
	const double crnd = 3.0*0x4000000*0x2000000;
  #endif
  #ifndef USE_AVX
	struct complex *ctmp;	// Hybrid AVX-DFT/SSE2-carry scheme used for Mersenne-mod needs a 2-word-double pointer
  #endif
	vec_dbl *tmp,*tm2;	// Non-static utility ptrs
  #ifndef MULTITHREAD
	vec_dbl *tm1;	// Non-static utility ptrs
  #endif
	static vec_dbl *two,*one,*sqrt2,*isrt2, *cc0, *ss0, *cc1, *ss1, *max_err, *sse2_rnd, *half_arr,
		// ptrs to 16 sets of twiddles shared by the 2nd-half DIF and DIT DFT macros:
		*twid0,*twid1,*twid2,*twid3,*twid4,*twid5,*twid6,*twid7,*twid8,*twid9,*twida,*twidb,*twidc,*twidd,*twide,*twidf,
		*r000,			// Head of RADIX*vec_cmplx-sized local store #1
	  #ifndef MULTITHREAD
		*r100,*r200,	// Head of RADIX*vec_cmplx-sized local store #1
		*s1p000,		// Head of RADIX*vec_cmplx-sized local store #2
	  #endif
		*cy;	// Need RADIX/2 slots for sse2 carries, RADIX/4 for avx
#else
  #ifndef MULTITHREAD
	static int p0123[4];
  #endif
#endif	// USE_SSE2?

#ifdef MULTITHREAD

	static struct cy_thread_data_t *tdat = 0x0;
	// Threadpool-based dispatch stuff:
  #if 0//def OS_TYPE_MACOSX
	static int main_work_units = 0;
  #endif
	static int pool_work_units = 0;
	static struct threadpool *tpool = 0x0;
	static int task_is_blocking = TRUE;
	static thread_control_t thread_control = {0,0,0};
	// First 3 subfields same for all threads, 4th provides thread-specifc data, will be inited at thread dispatch:
	static task_control_t   task_control = {NULL, (void*)cy768_process_chunk, NULL, 0x0};

#elif !defined(USE_SSE2)

	// Vars needed in scalar mode only:
	const double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
	int m,m2;
	double wt,wtinv,wtA,wtB,wtC;	/* Mersenne-mod weights stuff */
	double *addr;
	int bjmodn[RADIX];
	int *itmp;	// Pointer into the bjmodn array
	double temp,frac,cy[RADIX],
		t00,t01,t02,t03,t04,t05;

	static int t_offsets_lo[16], t_offsets_hi[16];

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
		// consisting of radix768_creals_in_local_store dcomplex and (12+RADIX/2) uint64 element slots per thread
		// (Add as many padding elts to the latter as needed to make it a multiple of 4):
		cslots_in_local_store = radix768_creals_in_local_store + (((12+RADIX/2)/2 + 3) & ~0x3);
		sc_arr = ALLOC_VEC_DBL(sc_arr, cslots_in_local_store*CY_THREADS);	if(!sc_arr){ sprintf(cbuf, "ERROR: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		sc_ptr = ALIGN_VEC_DBL(sc_arr);
		ASSERT(((intptr_t)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");
		sm_ptr = (uint64*)(sc_ptr + radix768_creals_in_local_store);
		ASSERT(((intptr_t)sm_ptr & 0x3f) == 0, "sm_ptr not 64-byte aligned!");

	  #ifdef USE_PTHREAD
		__r0 = sc_ptr;
	  #endif
		tmp = r000 = sc_ptr;
	  #ifndef MULTITHREAD
		r100 = tmp + 0x200;
		r200 = tmp + 0x400;
	  #endif
		tmp += 0x600;
	  #ifndef MULTITHREAD
						s1p000 = tmp;
	  #endif
		tmp += 0x600;	// sc_ptr += 0xc00
		two     = tmp + 0;	// AVX+ versions of radix-8,16,32 twiddleless-DFT macros need consts [2,1,sqrt2,isrt2] quartet laid out thusly
		one     = tmp + 1;
		sqrt2   = tmp + 2;
		isrt2   = tmp + 3;
		cc0		= tmp + 4;
		ss0		= tmp + 5;
		cc1		= tmp + 6;
		ss1		= tmp + 7;
		tmp += 0x08;	// sc_ptr += 0xc08
		// ptrs to 15 sets (30 vec_dbl data each) of non-unity twiddles shared by the 2nd-half DIF and DIT DFT macros:
		twid0  = tmp + 0x00;
		twid1  = tmp + 0x1e;
		twid2  = tmp + 0x3c;
		twid3  = tmp + 0x5a;
		twid4  = tmp + 0x78;
		twid5  = tmp + 0x96;
		twid6  = tmp + 0xb4;
		twid7  = tmp + 0xd2;
		twid8  = tmp + 0xf0;
		twid9  = tmp + 0x10e;
		twida  = tmp + 0x12c;
		twidb  = tmp + 0x14a;
		twidc  = tmp + 0x168;
		twidd  = tmp + 0x186;
		twide  = tmp + 0x1a4;
		twidf  = tmp + 0x1c2;
		tmp += 0x1e0;	// += 15*30 => sc_ptr += 0xde8
	  #ifdef USE_AVX512
		cy = tmp;		tmp += 0x60;	// RADIX/8 vec_dbl slots for carry sub-array
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;
		half_arr= tmp + 0x02;
	  #elif defined(USE_AVX)
		cy = tmp;		tmp += 0xc0;	// RADIX/4 vec_dbl slots for carry sub-array
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// += 0xc0 + 2 => sc_ptr += 0xeaa
		// This is where the value of half_arr_offset comes from
		half_arr= tmp + 0x02;	/* This table needs 96 vec_dbl for Mersenne-mod, and 3.5*RADIX[avx] | RADIX[sse2] for Fermat-mod */
	  #else
		cy = tmp;		tmp += 0x180;	// RADIX/2 vec_dbl slots for carry sub-array
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// += 0x180 + 2 => sc_ptr += 0xf6a
		// This is where the value of half_arr_offset comes from
		half_arr= tmp + 0x02;	/* This table needs 32 x 16 bytes for Mersenne-mod, 2 for Fermat-mod */
	  #endif
//		ASSERT(half_arr_offset == (uint32)(half_arr-sc_ptr), "half_arr_offset mismatches actual!");
		ASSERT((radix768_creals_in_local_store << L2_SZ_VD) >= ((intptr_t)half_arr - (intptr_t)r000) + (20 << L2_SZ_VD), "radix768_creals_in_local_store checksum failed!");

		/* These remain fixed: */
		VEC_DBL_INIT(two  , 2.0  );	VEC_DBL_INIT(one, 1.0  );
		dtmp = *(double *)&sqrt2_dn;	VEC_DBL_INIT(sqrt2, dtmp);
		dtmp = *(double *)&isrt2_dn;	VEC_DBL_INIT(isrt2, dtmp);
		VEC_DBL_INIT(cc0  ,  c16);
		VEC_DBL_INIT(ss0  ,  s16);
		/* cc0 = (cc1+cc2+cc3)/3 - 1; subtract 1 from Nussbaumer's definition in order to ease in-place computation */
		VEC_DBL_INIT(cc1, c3m1);
		VEC_DBL_INIT(ss1, s   );
		/* SSE2 math = 53-mantissa-bit IEEE double-float: */
	  #ifdef USE_AVX512	// In AVX-512 mode, use VRNDSCALEPD for rounding and hijack this vector-data slot for the 4 base/baseinv-consts
		sse2_rnd->d0 = base[0]; sse2_rnd->d1 = baseinv[1]; sse2_rnd->d2 = wts_mult[1]; sse2_rnd->d3 = inv_mult[0];
	  #else
		VEC_DBL_INIT(sse2_rnd, crnd);
	  #endif

		// ptrs to 16 sets (30 vec_dbl data each) of non-unity twiddles shared by the 2nd-half DIF and DIT DFT macros.
		// Since we copied the init-blocks here from the code below in which the twiddle-sets appear in BR order, init same way:
		// The first 2 sets (twid0/8) are processed in non-FMA fashion by both DIF/DIT macros, so init in non-FMA fashion:
		tmp = twid0; VEC_DBL_INIT(tmp,1); ++tmp; VEC_DBL_INIT(tmp,0); ++tmp; VEC_DBL_INIT(tmp,1); ++tmp; VEC_DBL_INIT(tmp,0); ++tmp; VEC_DBL_INIT(tmp,1); ++tmp; VEC_DBL_INIT(tmp,0); ++tmp; VEC_DBL_INIT(tmp,1); ++tmp; VEC_DBL_INIT(tmp,0); ++tmp; VEC_DBL_INIT(tmp,1); ++tmp; VEC_DBL_INIT(tmp,0); ++tmp; VEC_DBL_INIT(tmp,1); ++tmp; VEC_DBL_INIT(tmp,0); ++tmp; VEC_DBL_INIT(tmp,1); ++tmp; VEC_DBL_INIT(tmp,0); ++tmp; VEC_DBL_INIT(tmp,1); ++tmp; VEC_DBL_INIT(tmp,0); ++tmp; VEC_DBL_INIT(tmp,1); ++tmp; VEC_DBL_INIT(tmp,0); ++tmp; VEC_DBL_INIT(tmp,1); ++tmp; VEC_DBL_INIT(tmp,0); ++tmp; VEC_DBL_INIT(tmp,1); ++tmp; VEC_DBL_INIT(tmp,0); ++tmp; VEC_DBL_INIT(tmp,1); ++tmp; VEC_DBL_INIT(tmp,0); ++tmp; VEC_DBL_INIT(tmp,1); ++tmp; VEC_DBL_INIT(tmp,0); ++tmp; VEC_DBL_INIT(tmp,1); ++tmp; VEC_DBL_INIT(tmp,0); ++tmp; VEC_DBL_INIT(tmp,1); ++tmp; VEC_DBL_INIT(tmp,0); ++tmp;
		tmp = twid8; VEC_DBL_INIT(tmp,0); ++tmp; VEC_DBL_INIT(tmp,1); ++tmp; VEC_DBL_INIT(tmp, ISRT2); ++tmp; VEC_DBL_INIT(tmp,ISRT2); ++tmp; VEC_DBL_INIT(tmp, -ISRT2); ++tmp; VEC_DBL_INIT(tmp,ISRT2); ++tmp; VEC_DBL_INIT(tmp, c16); ++tmp; VEC_DBL_INIT(tmp,s16); ++tmp; VEC_DBL_INIT(tmp, -s16); ++tmp; VEC_DBL_INIT(tmp,c16); ++tmp; VEC_DBL_INIT(tmp, s16); ++tmp; VEC_DBL_INIT(tmp,c16); ++tmp; VEC_DBL_INIT(tmp, -c16); ++tmp; VEC_DBL_INIT(tmp,s16); ++tmp; VEC_DBL_INIT(tmp, c32_1); ++tmp; VEC_DBL_INIT(tmp,s32_1); ++tmp; VEC_DBL_INIT(tmp, -s32_1); ++tmp; VEC_DBL_INIT(tmp,c32_1); ++tmp; VEC_DBL_INIT(tmp, s32_3); ++tmp; VEC_DBL_INIT(tmp,c32_3); ++tmp; VEC_DBL_INIT(tmp, -c32_3); ++tmp; VEC_DBL_INIT(tmp,s32_3); ++tmp; VEC_DBL_INIT(tmp, c32_3); ++tmp; VEC_DBL_INIT(tmp,s32_3); ++tmp; VEC_DBL_INIT(tmp, -s32_3); ++tmp; VEC_DBL_INIT(tmp,c32_3); ++tmp; VEC_DBL_INIT(tmp, s32_1); ++tmp; VEC_DBL_INIT(tmp,c32_1); ++tmp; VEC_DBL_INIT(tmp, -c32_1); ++tmp; VEC_DBL_INIT(tmp,s32_1);

		// The remaining 14 sets are inited differently depending on whether SIMD+FMA is used:
	  #if !defined(USE_AVX2) || defined(USE_IMCI512)

		tmp = twid4; VEC_DBL_INIT(tmp,ISRT2); ++tmp; VEC_DBL_INIT(tmp,ISRT2); ++tmp; VEC_DBL_INIT(tmp, c16); ++tmp; VEC_DBL_INIT(tmp,s16); ++tmp; VEC_DBL_INIT(tmp, s16); ++tmp; VEC_DBL_INIT(tmp,c16); ++tmp; VEC_DBL_INIT(tmp, c32_1); ++tmp; VEC_DBL_INIT(tmp,s32_1); ++tmp; VEC_DBL_INIT(tmp, s32_3); ++tmp; VEC_DBL_INIT(tmp,c32_3); ++tmp; VEC_DBL_INIT(tmp, c32_3); ++tmp; VEC_DBL_INIT(tmp,s32_3); ++tmp; VEC_DBL_INIT(tmp, s32_1); ++tmp; VEC_DBL_INIT(tmp,c32_1); ++tmp; VEC_DBL_INIT(tmp, c64_1); ++tmp; VEC_DBL_INIT(tmp,s64_1); ++tmp; VEC_DBL_INIT(tmp, s64_7); ++tmp; VEC_DBL_INIT(tmp,c64_7); ++tmp; VEC_DBL_INIT(tmp, c64_5); ++tmp; VEC_DBL_INIT(tmp,s64_5); ++tmp; VEC_DBL_INIT(tmp, s64_3); ++tmp; VEC_DBL_INIT(tmp,c64_3); ++tmp; VEC_DBL_INIT(tmp, c64_3); ++tmp; VEC_DBL_INIT(tmp,s64_3); ++tmp; VEC_DBL_INIT(tmp, s64_5); ++tmp; VEC_DBL_INIT(tmp,c64_5); ++tmp; VEC_DBL_INIT(tmp, c64_7); ++tmp; VEC_DBL_INIT(tmp,s64_7); ++tmp; VEC_DBL_INIT(tmp, s64_1); ++tmp; VEC_DBL_INIT(tmp,c64_1);
		tmp = twidc; VEC_DBL_INIT(tmp,-ISRT2); ++tmp; VEC_DBL_INIT(tmp,ISRT2); ++tmp; VEC_DBL_INIT(tmp, s16); ++tmp; VEC_DBL_INIT(tmp,c16); ++tmp; VEC_DBL_INIT(tmp, -c16); ++tmp; VEC_DBL_INIT(tmp,-s16); ++tmp; VEC_DBL_INIT(tmp, c32_3); ++tmp; VEC_DBL_INIT(tmp,s32_3); ++tmp; VEC_DBL_INIT(tmp, -c32_1); ++tmp; VEC_DBL_INIT(tmp,s32_1); ++tmp; VEC_DBL_INIT(tmp, -s32_1); ++tmp; VEC_DBL_INIT(tmp,c32_1); ++tmp; VEC_DBL_INIT(tmp, -s32_3); ++tmp; VEC_DBL_INIT(tmp,-c32_3); ++tmp; VEC_DBL_INIT(tmp, c64_3); ++tmp; VEC_DBL_INIT(tmp,s64_3); ++tmp; VEC_DBL_INIT(tmp, -c64_5); ++tmp; VEC_DBL_INIT(tmp,s64_5); ++tmp; VEC_DBL_INIT(tmp, s64_1); ++tmp; VEC_DBL_INIT(tmp,c64_1); ++tmp; VEC_DBL_INIT(tmp, -c64_7); ++tmp; VEC_DBL_INIT(tmp,-s64_7); ++tmp; VEC_DBL_INIT(tmp, s64_7); ++tmp; VEC_DBL_INIT(tmp,c64_7); ++tmp; VEC_DBL_INIT(tmp, -c64_1); ++tmp; VEC_DBL_INIT(tmp,-s64_1); ++tmp; VEC_DBL_INIT(tmp, -s64_5); ++tmp; VEC_DBL_INIT(tmp,c64_5); ++tmp; VEC_DBL_INIT(tmp, -s64_3); ++tmp; VEC_DBL_INIT(tmp,-c64_3);
		tmp = twid2; VEC_DBL_INIT(tmp,c16); ++tmp; VEC_DBL_INIT(tmp,s16); ++tmp; VEC_DBL_INIT(tmp, c32_1); ++tmp; VEC_DBL_INIT(tmp,s32_1); ++tmp; VEC_DBL_INIT(tmp, c32_3); ++tmp; VEC_DBL_INIT(tmp,s32_3); ++tmp; VEC_DBL_INIT(tmp, c64_1); ++tmp; VEC_DBL_INIT(tmp,s64_1); ++tmp; VEC_DBL_INIT(tmp, c64_5); ++tmp; VEC_DBL_INIT(tmp,s64_5); ++tmp; VEC_DBL_INIT(tmp, c64_3); ++tmp; VEC_DBL_INIT(tmp,s64_3); ++tmp; VEC_DBL_INIT(tmp, c64_7); ++tmp; VEC_DBL_INIT(tmp,s64_7); ++tmp; VEC_DBL_INIT(tmp, c128_1); ++tmp; VEC_DBL_INIT(tmp,s128_1); ++tmp; VEC_DBL_INIT(tmp, c128_9); ++tmp; VEC_DBL_INIT(tmp,s128_9); ++tmp; VEC_DBL_INIT(tmp, c128_5); ++tmp; VEC_DBL_INIT(tmp,s128_5); ++tmp; VEC_DBL_INIT(tmp, c128_d); ++tmp; VEC_DBL_INIT(tmp,s128_d); ++tmp; VEC_DBL_INIT(tmp, c128_3); ++tmp; VEC_DBL_INIT(tmp,s128_3); ++tmp; VEC_DBL_INIT(tmp, c128_b); ++tmp; VEC_DBL_INIT(tmp,s128_b); ++tmp; VEC_DBL_INIT(tmp, c128_7); ++tmp; VEC_DBL_INIT(tmp,s128_7); ++tmp; VEC_DBL_INIT(tmp, c128_f); ++tmp; VEC_DBL_INIT(tmp,s128_f);
		tmp = twida; VEC_DBL_INIT(tmp,-s16); ++tmp; VEC_DBL_INIT(tmp,c16); ++tmp; VEC_DBL_INIT(tmp, s32_3); ++tmp; VEC_DBL_INIT(tmp,c32_3); ++tmp; VEC_DBL_INIT(tmp, -c32_1); ++tmp; VEC_DBL_INIT(tmp,s32_1); ++tmp; VEC_DBL_INIT(tmp, c64_5); ++tmp; VEC_DBL_INIT(tmp,s64_5); ++tmp; VEC_DBL_INIT(tmp, -c64_7); ++tmp; VEC_DBL_INIT(tmp,s64_7); ++tmp; VEC_DBL_INIT(tmp, s64_1); ++tmp; VEC_DBL_INIT(tmp,c64_1); ++tmp; VEC_DBL_INIT(tmp, -c64_3); ++tmp; VEC_DBL_INIT(tmp,-s64_3); ++tmp; VEC_DBL_INIT(tmp, c128_5); ++tmp; VEC_DBL_INIT(tmp,s128_5); ++tmp; VEC_DBL_INIT(tmp, -s128_d); ++tmp; VEC_DBL_INIT(tmp,c128_d); ++tmp; VEC_DBL_INIT(tmp, s128_7); ++tmp; VEC_DBL_INIT(tmp,c128_7); ++tmp; VEC_DBL_INIT(tmp, -c128_1); ++tmp; VEC_DBL_INIT(tmp,-s128_1); ++tmp; VEC_DBL_INIT(tmp, c128_f); ++tmp; VEC_DBL_INIT(tmp,s128_f); ++tmp; VEC_DBL_INIT(tmp, -c128_9); ++tmp; VEC_DBL_INIT(tmp,s128_9); ++tmp; VEC_DBL_INIT(tmp, -s128_3); ++tmp; VEC_DBL_INIT(tmp,c128_3); ++tmp; VEC_DBL_INIT(tmp, -c128_b); ++tmp; VEC_DBL_INIT(tmp,-s128_b);
		tmp = twid6; VEC_DBL_INIT(tmp,s16); ++tmp; VEC_DBL_INIT(tmp,c16); ++tmp; VEC_DBL_INIT(tmp, c32_3); ++tmp; VEC_DBL_INIT(tmp,s32_3); ++tmp; VEC_DBL_INIT(tmp, -s32_1); ++tmp; VEC_DBL_INIT(tmp,c32_1); ++tmp; VEC_DBL_INIT(tmp, c64_3); ++tmp; VEC_DBL_INIT(tmp,s64_3); ++tmp; VEC_DBL_INIT(tmp, s64_1); ++tmp; VEC_DBL_INIT(tmp,c64_1); ++tmp; VEC_DBL_INIT(tmp, s64_7); ++tmp; VEC_DBL_INIT(tmp,c64_7); ++tmp; VEC_DBL_INIT(tmp, -s64_5); ++tmp; VEC_DBL_INIT(tmp,c64_5); ++tmp; VEC_DBL_INIT(tmp, c128_3); ++tmp; VEC_DBL_INIT(tmp,s128_3); ++tmp; VEC_DBL_INIT(tmp, s128_5); ++tmp; VEC_DBL_INIT(tmp,c128_5); ++tmp; VEC_DBL_INIT(tmp, c128_f); ++tmp; VEC_DBL_INIT(tmp,s128_f); ++tmp; VEC_DBL_INIT(tmp, -s128_7); ++tmp; VEC_DBL_INIT(tmp,c128_7); ++tmp; VEC_DBL_INIT(tmp, c128_9); ++tmp; VEC_DBL_INIT(tmp,s128_9); ++tmp; VEC_DBL_INIT(tmp, -s128_1); ++tmp; VEC_DBL_INIT(tmp,c128_1); ++tmp; VEC_DBL_INIT(tmp, s128_b); ++tmp; VEC_DBL_INIT(tmp,c128_b); ++tmp; VEC_DBL_INIT(tmp, -s128_d); ++tmp; VEC_DBL_INIT(tmp,c128_d);
		tmp = twide; VEC_DBL_INIT(tmp,-c16); ++tmp; VEC_DBL_INIT(tmp,s16); ++tmp; VEC_DBL_INIT(tmp, s32_1); ++tmp; VEC_DBL_INIT(tmp,c32_1); ++tmp; VEC_DBL_INIT(tmp, -s32_3); ++tmp; VEC_DBL_INIT(tmp,-c32_3); ++tmp; VEC_DBL_INIT(tmp, c64_7); ++tmp; VEC_DBL_INIT(tmp,s64_7); ++tmp; VEC_DBL_INIT(tmp, -c64_3); ++tmp; VEC_DBL_INIT(tmp,-s64_3); ++tmp; VEC_DBL_INIT(tmp, -s64_5); ++tmp; VEC_DBL_INIT(tmp,c64_5); ++tmp; VEC_DBL_INIT(tmp, s64_1); ++tmp; VEC_DBL_INIT(tmp,-c64_1); ++tmp; VEC_DBL_INIT(tmp, c128_7); ++tmp; VEC_DBL_INIT(tmp,s128_7); ++tmp; VEC_DBL_INIT(tmp, -c128_1); ++tmp; VEC_DBL_INIT(tmp,s128_1); ++tmp; VEC_DBL_INIT(tmp, -s128_3); ++tmp; VEC_DBL_INIT(tmp,c128_3); ++tmp; VEC_DBL_INIT(tmp, -s128_5); ++tmp; VEC_DBL_INIT(tmp,-c128_5); ++tmp; VEC_DBL_INIT(tmp, s128_b); ++tmp; VEC_DBL_INIT(tmp,c128_b); ++tmp; VEC_DBL_INIT(tmp, -c128_d); ++tmp; VEC_DBL_INIT(tmp,-s128_d); ++tmp; VEC_DBL_INIT(tmp, -c128_f); ++tmp; VEC_DBL_INIT(tmp,s128_f); ++tmp; VEC_DBL_INIT(tmp, s128_9); ++tmp; VEC_DBL_INIT(tmp,-c128_9);
		tmp = twid1; VEC_DBL_INIT(tmp,c32_1); ++tmp; VEC_DBL_INIT(tmp,s32_1); ++tmp; VEC_DBL_INIT(tmp, c64_1); ++tmp; VEC_DBL_INIT(tmp,s64_1); ++tmp; VEC_DBL_INIT(tmp, c64_3); ++tmp; VEC_DBL_INIT(tmp,s64_3); ++tmp; VEC_DBL_INIT(tmp, c128_1); ++tmp; VEC_DBL_INIT(tmp,s128_1); ++tmp; VEC_DBL_INIT(tmp, c128_5); ++tmp; VEC_DBL_INIT(tmp,s128_5); ++tmp; VEC_DBL_INIT(tmp, c128_3); ++tmp; VEC_DBL_INIT(tmp,s128_3); ++tmp; VEC_DBL_INIT(tmp, c128_7); ++tmp; VEC_DBL_INIT(tmp,s128_7); ++tmp; VEC_DBL_INIT(tmp, c256_01); ++tmp; VEC_DBL_INIT(tmp,s256_01); ++tmp; VEC_DBL_INIT(tmp, c256_09); ++tmp; VEC_DBL_INIT(tmp,s256_09); ++tmp; VEC_DBL_INIT(tmp, c256_05); ++tmp; VEC_DBL_INIT(tmp,s256_05); ++tmp; VEC_DBL_INIT(tmp, c256_0d); ++tmp; VEC_DBL_INIT(tmp,s256_0d); ++tmp; VEC_DBL_INIT(tmp, c256_03); ++tmp; VEC_DBL_INIT(tmp,s256_03); ++tmp; VEC_DBL_INIT(tmp, c256_0b); ++tmp; VEC_DBL_INIT(tmp,s256_0b); ++tmp; VEC_DBL_INIT(tmp, c256_07); ++tmp; VEC_DBL_INIT(tmp,s256_07); ++tmp; VEC_DBL_INIT(tmp, c256_0f); ++tmp; VEC_DBL_INIT(tmp,s256_0f);
		tmp = twid9; VEC_DBL_INIT(tmp,-s32_1); ++tmp; VEC_DBL_INIT(tmp,c32_1); ++tmp; VEC_DBL_INIT(tmp, s64_7); ++tmp; VEC_DBL_INIT(tmp,c64_7); ++tmp; VEC_DBL_INIT(tmp, -c64_5); ++tmp; VEC_DBL_INIT(tmp,s64_5); ++tmp; VEC_DBL_INIT(tmp, c128_9); ++tmp; VEC_DBL_INIT(tmp,s128_9); ++tmp; VEC_DBL_INIT(tmp, -s128_d); ++tmp; VEC_DBL_INIT(tmp,c128_d); ++tmp; VEC_DBL_INIT(tmp, s128_5); ++tmp; VEC_DBL_INIT(tmp,c128_5); ++tmp; VEC_DBL_INIT(tmp, -c128_1); ++tmp; VEC_DBL_INIT(tmp,s128_1); ++tmp; VEC_DBL_INIT(tmp, c256_09); ++tmp; VEC_DBL_INIT(tmp,s256_09); ++tmp; VEC_DBL_INIT(tmp, -s256_11); ++tmp; VEC_DBL_INIT(tmp,c256_11); ++tmp; VEC_DBL_INIT(tmp, s256_13); ++tmp; VEC_DBL_INIT(tmp,c256_13); ++tmp; VEC_DBL_INIT(tmp, -c256_0b); ++tmp; VEC_DBL_INIT(tmp,s256_0b); ++tmp; VEC_DBL_INIT(tmp, c256_1b); ++tmp; VEC_DBL_INIT(tmp,s256_1b); ++tmp; VEC_DBL_INIT(tmp, -c256_1d); ++tmp; VEC_DBL_INIT(tmp,s256_1d); ++tmp; VEC_DBL_INIT(tmp, s256_01); ++tmp; VEC_DBL_INIT(tmp,c256_01); ++tmp; VEC_DBL_INIT(tmp, -c256_07); ++tmp; VEC_DBL_INIT(tmp,-s256_07);
		tmp = twid5; VEC_DBL_INIT(tmp,s32_3); ++tmp; VEC_DBL_INIT(tmp,c32_3); ++tmp; VEC_DBL_INIT(tmp, c64_5); ++tmp; VEC_DBL_INIT(tmp,s64_5); ++tmp; VEC_DBL_INIT(tmp, s64_1); ++tmp; VEC_DBL_INIT(tmp,c64_1); ++tmp; VEC_DBL_INIT(tmp, c128_5); ++tmp; VEC_DBL_INIT(tmp,s128_5); ++tmp; VEC_DBL_INIT(tmp, s128_7); ++tmp; VEC_DBL_INIT(tmp,c128_7); ++tmp; VEC_DBL_INIT(tmp, c128_f); ++tmp; VEC_DBL_INIT(tmp,s128_f); ++tmp; VEC_DBL_INIT(tmp, -s128_3); ++tmp; VEC_DBL_INIT(tmp,c128_3); ++tmp; VEC_DBL_INIT(tmp, c256_05); ++tmp; VEC_DBL_INIT(tmp,s256_05); ++tmp; VEC_DBL_INIT(tmp, s256_13); ++tmp; VEC_DBL_INIT(tmp,c256_13); ++tmp; VEC_DBL_INIT(tmp, c256_19); ++tmp; VEC_DBL_INIT(tmp,s256_19); ++tmp; VEC_DBL_INIT(tmp, -s256_01); ++tmp; VEC_DBL_INIT(tmp,c256_01); ++tmp; VEC_DBL_INIT(tmp, c256_0f); ++tmp; VEC_DBL_INIT(tmp,s256_0f); ++tmp; VEC_DBL_INIT(tmp, s256_09); ++tmp; VEC_DBL_INIT(tmp,c256_09); ++tmp; VEC_DBL_INIT(tmp, s256_1d); ++tmp; VEC_DBL_INIT(tmp,c256_1d); ++tmp; VEC_DBL_INIT(tmp, -s256_0b); ++tmp; VEC_DBL_INIT(tmp,c256_0b);
		tmp = twidd; VEC_DBL_INIT(tmp,-c32_3); ++tmp; VEC_DBL_INIT(tmp,s32_3); ++tmp; VEC_DBL_INIT(tmp, s64_3); ++tmp; VEC_DBL_INIT(tmp,c64_3); ++tmp; VEC_DBL_INIT(tmp, -c64_7); ++tmp; VEC_DBL_INIT(tmp,-s64_7); ++tmp; VEC_DBL_INIT(tmp, c128_d); ++tmp; VEC_DBL_INIT(tmp,s128_d); ++tmp; VEC_DBL_INIT(tmp, -c128_1); ++tmp; VEC_DBL_INIT(tmp,-s128_1); ++tmp; VEC_DBL_INIT(tmp, -s128_7); ++tmp; VEC_DBL_INIT(tmp,c128_7); ++tmp; VEC_DBL_INIT(tmp, -s128_5); ++tmp; VEC_DBL_INIT(tmp,-c128_5); ++tmp; VEC_DBL_INIT(tmp, c256_0d); ++tmp; VEC_DBL_INIT(tmp,s256_0d); ++tmp; VEC_DBL_INIT(tmp, -c256_0b); ++tmp; VEC_DBL_INIT(tmp,s256_0b); ++tmp; VEC_DBL_INIT(tmp, -s256_01); ++tmp; VEC_DBL_INIT(tmp,c256_01); ++tmp; VEC_DBL_INIT(tmp, -s256_17); ++tmp; VEC_DBL_INIT(tmp,-c256_17); ++tmp; VEC_DBL_INIT(tmp, s256_19); ++tmp; VEC_DBL_INIT(tmp,c256_19); ++tmp; VEC_DBL_INIT(tmp, -c256_0f); ++tmp; VEC_DBL_INIT(tmp,-s256_0f); ++tmp; VEC_DBL_INIT(tmp, -s256_1b); ++tmp; VEC_DBL_INIT(tmp,c256_1b); ++tmp; VEC_DBL_INIT(tmp, s256_03); ++tmp; VEC_DBL_INIT(tmp,-c256_03);
		tmp = twid3; VEC_DBL_INIT(tmp,c32_3); ++tmp; VEC_DBL_INIT(tmp,s32_3); ++tmp; VEC_DBL_INIT(tmp, c64_3); ++tmp; VEC_DBL_INIT(tmp,s64_3); ++tmp; VEC_DBL_INIT(tmp, s64_7); ++tmp; VEC_DBL_INIT(tmp,c64_7); ++tmp; VEC_DBL_INIT(tmp, c128_3); ++tmp; VEC_DBL_INIT(tmp,s128_3); ++tmp; VEC_DBL_INIT(tmp, c128_f); ++tmp; VEC_DBL_INIT(tmp,s128_f); ++tmp; VEC_DBL_INIT(tmp, c128_9); ++tmp; VEC_DBL_INIT(tmp,s128_9); ++tmp; VEC_DBL_INIT(tmp, s128_b); ++tmp; VEC_DBL_INIT(tmp,c128_b); ++tmp; VEC_DBL_INIT(tmp, c256_03); ++tmp; VEC_DBL_INIT(tmp,s256_03); ++tmp; VEC_DBL_INIT(tmp, c256_1b); ++tmp; VEC_DBL_INIT(tmp,s256_1b); ++tmp; VEC_DBL_INIT(tmp, c256_0f); ++tmp; VEC_DBL_INIT(tmp,s256_0f); ++tmp; VEC_DBL_INIT(tmp, s256_19); ++tmp; VEC_DBL_INIT(tmp,c256_19); ++tmp; VEC_DBL_INIT(tmp, c256_09); ++tmp; VEC_DBL_INIT(tmp,s256_09); ++tmp; VEC_DBL_INIT(tmp, s256_1f); ++tmp; VEC_DBL_INIT(tmp,c256_1f); ++tmp; VEC_DBL_INIT(tmp, c256_15); ++tmp; VEC_DBL_INIT(tmp,s256_15); ++tmp; VEC_DBL_INIT(tmp, s256_13); ++tmp; VEC_DBL_INIT(tmp,c256_13);
		tmp = twidb; VEC_DBL_INIT(tmp,-s32_3); ++tmp; VEC_DBL_INIT(tmp,c32_3); ++tmp; VEC_DBL_INIT(tmp, s64_5); ++tmp; VEC_DBL_INIT(tmp,c64_5); ++tmp; VEC_DBL_INIT(tmp, -c64_1); ++tmp; VEC_DBL_INIT(tmp,-s64_1); ++tmp; VEC_DBL_INIT(tmp, c128_b); ++tmp; VEC_DBL_INIT(tmp,s128_b); ++tmp; VEC_DBL_INIT(tmp, -c128_9); ++tmp; VEC_DBL_INIT(tmp,s128_9); ++tmp; VEC_DBL_INIT(tmp, -s128_1); ++tmp; VEC_DBL_INIT(tmp,c128_1); ++tmp; VEC_DBL_INIT(tmp, -c128_d); ++tmp; VEC_DBL_INIT(tmp,-s128_d); ++tmp; VEC_DBL_INIT(tmp, c256_0b); ++tmp; VEC_DBL_INIT(tmp,s256_0b); ++tmp; VEC_DBL_INIT(tmp, -c256_1d); ++tmp; VEC_DBL_INIT(tmp,s256_1d); ++tmp; VEC_DBL_INIT(tmp, s256_09); ++tmp; VEC_DBL_INIT(tmp,c256_09); ++tmp; VEC_DBL_INIT(tmp, -c256_0f); ++tmp; VEC_DBL_INIT(tmp,-s256_0f); ++tmp; VEC_DBL_INIT(tmp, s256_1f); ++tmp; VEC_DBL_INIT(tmp,c256_1f); ++tmp; VEC_DBL_INIT(tmp, -c256_07); ++tmp; VEC_DBL_INIT(tmp,s256_07); ++tmp; VEC_DBL_INIT(tmp, -s256_0d); ++tmp; VEC_DBL_INIT(tmp,c256_0d); ++tmp; VEC_DBL_INIT(tmp, -s256_1b); ++tmp; VEC_DBL_INIT(tmp,-c256_1b);
		tmp = twid7; VEC_DBL_INIT(tmp,s32_1); ++tmp; VEC_DBL_INIT(tmp,c32_1); ++tmp; VEC_DBL_INIT(tmp, c64_7); ++tmp; VEC_DBL_INIT(tmp,s64_7); ++tmp; VEC_DBL_INIT(tmp, -s64_5); ++tmp; VEC_DBL_INIT(tmp,c64_5); ++tmp; VEC_DBL_INIT(tmp, c128_7); ++tmp; VEC_DBL_INIT(tmp,s128_7); ++tmp; VEC_DBL_INIT(tmp, -s128_3); ++tmp; VEC_DBL_INIT(tmp,c128_3); ++tmp; VEC_DBL_INIT(tmp, s128_b); ++tmp; VEC_DBL_INIT(tmp,c128_b); ++tmp; VEC_DBL_INIT(tmp, -c128_f); ++tmp; VEC_DBL_INIT(tmp,s128_f); ++tmp; VEC_DBL_INIT(tmp, c256_07); ++tmp; VEC_DBL_INIT(tmp,s256_07); ++tmp; VEC_DBL_INIT(tmp, s256_01); ++tmp; VEC_DBL_INIT(tmp,c256_01); ++tmp; VEC_DBL_INIT(tmp, s256_1d); ++tmp; VEC_DBL_INIT(tmp,c256_1d); ++tmp; VEC_DBL_INIT(tmp, -s256_1b); ++tmp; VEC_DBL_INIT(tmp,c256_1b); ++tmp; VEC_DBL_INIT(tmp, c256_15); ++tmp; VEC_DBL_INIT(tmp,s256_15); ++tmp; VEC_DBL_INIT(tmp, -s256_0d); ++tmp; VEC_DBL_INIT(tmp,c256_0d); ++tmp; VEC_DBL_INIT(tmp, s256_0f); ++tmp; VEC_DBL_INIT(tmp,c256_0f); ++tmp; VEC_DBL_INIT(tmp, -c256_17); ++tmp; VEC_DBL_INIT(tmp,s256_17);
		tmp = twidf; VEC_DBL_INIT(tmp,-c32_1); ++tmp; VEC_DBL_INIT(tmp,s32_1); ++tmp; VEC_DBL_INIT(tmp, s64_1); ++tmp; VEC_DBL_INIT(tmp,c64_1); ++tmp; VEC_DBL_INIT(tmp, -s64_3); ++tmp; VEC_DBL_INIT(tmp,-c64_3); ++tmp; VEC_DBL_INIT(tmp, c128_f); ++tmp; VEC_DBL_INIT(tmp,s128_f); ++tmp; VEC_DBL_INIT(tmp, -c128_b); ++tmp; VEC_DBL_INIT(tmp,-s128_b); ++tmp; VEC_DBL_INIT(tmp, -s128_d); ++tmp; VEC_DBL_INIT(tmp,c128_d); ++tmp; VEC_DBL_INIT(tmp, s128_9); ++tmp; VEC_DBL_INIT(tmp,-c128_9); ++tmp; VEC_DBL_INIT(tmp, c256_0f); ++tmp; VEC_DBL_INIT(tmp,s256_0f); ++tmp; VEC_DBL_INIT(tmp, -c256_07); ++tmp; VEC_DBL_INIT(tmp,-s256_07); ++tmp; VEC_DBL_INIT(tmp, -s256_0b); ++tmp; VEC_DBL_INIT(tmp,c256_0b); ++tmp; VEC_DBL_INIT(tmp, s256_03); ++tmp; VEC_DBL_INIT(tmp,-c256_03); ++tmp; VEC_DBL_INIT(tmp, s256_13); ++tmp; VEC_DBL_INIT(tmp,c256_13); ++tmp; VEC_DBL_INIT(tmp, -s256_1b); ++tmp; VEC_DBL_INIT(tmp,-c256_1b); ++tmp; VEC_DBL_INIT(tmp, -c256_17); ++tmp; VEC_DBL_INIT(tmp,s256_17); ++tmp; VEC_DBL_INIT(tmp, c256_1f); ++tmp; VEC_DBL_INIT(tmp,-s256_1f);

	  #else // USE_AVX2 = true:

		// Precompute the FMA-modified twiddles for the 2nd-pass radix-16 DFTs:
		#ifdef USE_FMA
			#error USE_FMA flag not supported in SIMD mode - to use FMA under AVX2/FMA3, define *only* USE_AVX2!
		#endif

		#include "radix16_dif_dit_pass_asm.h"	// Need this for FMA_TWIDDLE_FIDDLE macro

		FMA_TWIDDLE_FIDDLE(
			ISRT2,ISRT2, c16,s16, s16,c16, c32_1,s32_1, s32_3,c32_3, c32_3,s32_3, s32_1,c32_1, c64_1,s64_1, s64_7,c64_7, c64_5,s64_5, s64_3,c64_3, c64_3,s64_3, s64_5,c64_5, c64_7,s64_7, s64_1,c64_1	,c16,tan,
			twid4)
		FMA_TWIDDLE_FIDDLE(
			-ISRT2,ISRT2, s16,c16, -c16,-s16, c32_3,s32_3, -c32_1,s32_1, -s32_1,c32_1, -s32_3,-c32_3, c64_3,s64_3, -c64_5,s64_5, s64_1,c64_1, -c64_7,-s64_7, s64_7,c64_7, -c64_1,-s64_1, -s64_5,c64_5, -s64_3,-c64_3	,c16,tan,
			twidc)
		FMA_TWIDDLE_FIDDLE(
			c16,s16, c32_1,s32_1, c32_3,s32_3, c64_1,s64_1, c64_5,s64_5, c64_3,s64_3, c64_7,s64_7, c128_1,s128_1, c128_9,s128_9, c128_5,s128_5, c128_d,s128_d, c128_3,s128_3, c128_b,s128_b, c128_7,s128_7, c128_f,s128_f	,c16,tan,
			twid2)
		FMA_TWIDDLE_FIDDLE(
			-s16,c16, s32_3,c32_3, -c32_1,s32_1, c64_5,s64_5, -c64_7,s64_7, s64_1,c64_1, -c64_3,-s64_3, c128_5,s128_5, -s128_d,c128_d, s128_7,c128_7, -c128_1,-s128_1, c128_f,s128_f, -c128_9,s128_9, -s128_3,c128_3, -c128_b,-s128_b	,c16,tan,
			twida)
		FMA_TWIDDLE_FIDDLE(
			s16,c16, c32_3,s32_3, -s32_1,c32_1, c64_3,s64_3, s64_1,c64_1, s64_7,c64_7, -s64_5,c64_5, c128_3,s128_3, s128_5,c128_5, c128_f,s128_f, -s128_7,c128_7, c128_9,s128_9, -s128_1,c128_1, s128_b,c128_b, -s128_d,c128_d	,c16,tan,
			twid6)
		FMA_TWIDDLE_FIDDLE(
			-c16,s16, s32_1,c32_1, -s32_3,-c32_3, c64_7,s64_7, -c64_3,-s64_3, -s64_5,c64_5, s64_1,-c64_1, c128_7,s128_7, -c128_1,s128_1, -s128_3,c128_3, -s128_5,-c128_5, s128_b,c128_b, -c128_d,-s128_d, -c128_f,s128_f, s128_9,-c128_9	,c16,tan,
			twide)
		FMA_TWIDDLE_FIDDLE(
			c32_1,s32_1, c64_1,s64_1, c64_3,s64_3, c128_1,s128_1, c128_5,s128_5, c128_3,s128_3, c128_7,s128_7, c256_01,s256_01, c256_09,s256_09, c256_05,s256_05, c256_0d,s256_0d, c256_03,s256_03, c256_0b,s256_0b, c256_07,s256_07, c256_0f,s256_0f	,c16,tan,
			twid1)
		FMA_TWIDDLE_FIDDLE(
			-s32_1,c32_1, s64_7,c64_7, -c64_5,s64_5, c128_9,s128_9, -s128_d,c128_d, s128_5,c128_5, -c128_1,s128_1, c256_09,s256_09, -s256_11,c256_11, s256_13,c256_13, -c256_0b,s256_0b, c256_1b,s256_1b, -c256_1d,s256_1d, s256_01,c256_01, -c256_07,-s256_07	,c16,tan,
			twid9)
		FMA_TWIDDLE_FIDDLE(
			s32_3,c32_3, c64_5,s64_5, s64_1,c64_1, c128_5,s128_5, s128_7,c128_7, c128_f,s128_f, -s128_3,c128_3, c256_05,s256_05, s256_13,c256_13, c256_19,s256_19, -s256_01,c256_01, c256_0f,s256_0f, s256_09,c256_09, s256_1d,c256_1d, -s256_0b,c256_0b	,c16,tan,
			twid5)
		FMA_TWIDDLE_FIDDLE(
			-c32_3,s32_3, s64_3,c64_3, -c64_7,-s64_7, c128_d,s128_d, -c128_1,-s128_1, -s128_7,c128_7, -s128_5,-c128_5, c256_0d,s256_0d, -c256_0b,s256_0b, -s256_01,c256_01, -s256_17,-c256_17, s256_19,c256_19, -c256_0f,-s256_0f, -s256_1b,c256_1b, s256_03,-c256_03	,c16,tan,
			twidd)
		FMA_TWIDDLE_FIDDLE(
			c32_3,s32_3, c64_3,s64_3, s64_7,c64_7, c128_3,s128_3, c128_f,s128_f, c128_9,s128_9, s128_b,c128_b, c256_03,s256_03, c256_1b,s256_1b, c256_0f,s256_0f, s256_19,c256_19, c256_09,s256_09, s256_1f,c256_1f, c256_15,s256_15, s256_13,c256_13	,c16,tan,
			twid3)
		FMA_TWIDDLE_FIDDLE(
			-s32_3,c32_3, s64_5,c64_5, -c64_1,-s64_1, c128_b,s128_b, -c128_9,s128_9, -s128_1,c128_1, -c128_d,-s128_d, c256_0b,s256_0b, -c256_1d,s256_1d, s256_09,c256_09, -c256_0f,-s256_0f, s256_1f,c256_1f, -c256_07,s256_07, -s256_0d,c256_0d, -s256_1b,-c256_1b	,c16,tan,
			twidb)
		FMA_TWIDDLE_FIDDLE(
			s32_1,c32_1, c64_7,s64_7, -s64_5,c64_5, c128_7,s128_7, -s128_3,c128_3, s128_b,c128_b, -c128_f,s128_f, c256_07,s256_07, s256_01,c256_01, s256_1d,c256_1d, -s256_1b,c256_1b, c256_15,s256_15, -s256_0d,c256_0d, s256_0f,c256_0f, -c256_17,s256_17	,c16,tan,
			twid7)
		FMA_TWIDDLE_FIDDLE(
			-c32_1,s32_1, s64_1,c64_1, -s64_3,-c64_3, c128_f,s128_f, -c128_b,-s128_b, -s128_d,c128_d, s128_9,-c128_9, c256_0f,s256_0f, -c256_07,-s256_07, -s256_0b,c256_0b, s256_03,-c256_03, s256_13,c256_13, -s256_1b,-c256_1b, -c256_17,s256_17, c256_1f,-s256_1f	,c16,tan,
			twidf)

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
		#ifndef MULTITHREAD
		n_minus_sil   = (struct uint32x16*)sse_n + 1;
		n_minus_silp1 = (struct uint32x16*)sse_n + 2;
		sinwt         = (struct uint32x16*)sse_n + 3;
		sinwtm1       = (struct uint32x16*)sse_n + 4;
		#endif
		nbytes += 256;
	   #else
		#ifndef MULTITHREAD
		n_minus_sil   = (struct uint32x8 *)sse_n + 1;
		n_minus_silp1 = (struct uint32x8 *)sse_n + 2;
		sinwt         = (struct uint32x8 *)sse_n + 3;
		sinwtm1       = (struct uint32x8 *)sse_n + 4;
		#endif
		nbytes += 128;
	   #endif
	  #elif defined(USE_AVX)
	   #ifndef MULTITHREAD
		n_minus_sil   = (struct uint32x4 *)sse_n + 1;
		n_minus_silp1 = (struct uint32x4 *)sse_n + 2;
		sinwt         = (struct uint32x4 *)sse_n + 3;
		sinwtm1       = (struct uint32x4 *)sse_n + 4;
	   #endif
		nbytes += 64;
	  #endif

		// Propagate the above consts to the remaining threads:
		tmp = (vec_dbl *)sm_ptr;
		tm2 = tmp + cslots_in_local_store;
		for(ithread = 1; ithread < CY_THREADS; ++ithread) {
			memcpy(tm2, tmp, nbytes);
			tmp = tm2;		tm2 += cslots_in_local_store;
		}

	 #ifndef MULTITHREAD
	// For large radices, array-access to bjmodn means only init base-ptr here:
	  #ifdef USE_AVX
		bjmodn = (int*)(sinwtm1 + RE_IM_STRIDE);
	  #else
		bjmodn = (int*)(sse_n   + RE_IM_STRIDE);
	  #endif
	 #endif

	#endif	// USE_SSE2

		/*   constant index offsets for array load/stores are here.	*/
		pini = NDIVR/CY_THREADS;
		p1 = NDIVR;		p10 = NDIVR<<4;		p100 = NDIVR<<8;		p200 = NDIVR<<9;
		p2 = p1 + p1;	p20 = p10 + p10;	p110 = p100 + p10;		p210 = p200 + p10;
		p3 = p2 + p1;	p30 = p20 + p10;	p120 = p110 + p10;		p220 = p210 + p10;
		p4 = p3 + p1;	p40 = p30 + p10;	p130 = p120 + p10;		p230 = p220 + p10;
		p5 = p4 + p1;	p50 = p40 + p10;	p140 = p130 + p10;		p240 = p230 + p10;
		p6 = p5 + p1;	p60 = p50 + p10;	p150 = p140 + p10;		p250 = p240 + p10;
		p7 = p6 + p1;	p70 = p60 + p10;	p160 = p150 + p10;		p260 = p250 + p10;
		p8 = p7 + p1;	p80 = p70 + p10;	p170 = p160 + p10;		p270 = p260 + p10;
		p9 = p8 + p1;	p90 = p80 + p10;	p180 = p170 + p10;		p280 = p270 + p10;
		pa = p9 + p1;	pa0 = p90 + p10;	p190 = p180 + p10;		p290 = p280 + p10;
		pb = pa + p1;	pb0 = pa0 + p10;	p1a0 = p190 + p10;		p2a0 = p290 + p10;
		pc = pb + p1;	pc0 = pb0 + p10;	p1b0 = p1a0 + p10;		p2b0 = p2a0 + p10;
		pd = pc + p1;	pd0 = pc0 + p10;	p1c0 = p1b0 + p10;		p2c0 = p2b0 + p10;
		pe = pd + p1;	pe0 = pd0 + p10;	p1d0 = p1c0 + p10;		p2d0 = p2c0 + p10;
		pf = pe + p1;	pf0 = pe0 + p10;	p1e0 = p1d0 + p10;		p2e0 = p2d0 + p10;
											p1f0 = p1e0 + p10;		p2f0 = p2e0 + p10;
		p1 += ( (p1 >> DAT_BITS) << PAD_BITS );
		p2 += ( (p2 >> DAT_BITS) << PAD_BITS );
		p3 += ( (p3 >> DAT_BITS) << PAD_BITS );
		p4 += ( (p4 >> DAT_BITS) << PAD_BITS );
		p5 += ( (p5 >> DAT_BITS) << PAD_BITS );
		p6 += ( (p6 >> DAT_BITS) << PAD_BITS );
		p7 += ( (p7 >> DAT_BITS) << PAD_BITS );
		p8 += ( (p8 >> DAT_BITS) << PAD_BITS );
		p9 += ( (p9 >> DAT_BITS) << PAD_BITS );
		pa += ( (pa >> DAT_BITS) << PAD_BITS );
		pb += ( (pb >> DAT_BITS) << PAD_BITS );
		pc += ( (pc >> DAT_BITS) << PAD_BITS );
		pd += ( (pd >> DAT_BITS) << PAD_BITS );
		pe += ( (pe >> DAT_BITS) << PAD_BITS );
		pf += ( (pf >> DAT_BITS) << PAD_BITS );
		p10 += ( (p10 >> DAT_BITS) << PAD_BITS );
		p20 += ( (p20 >> DAT_BITS) << PAD_BITS );
		p30 += ( (p30 >> DAT_BITS) << PAD_BITS );
		p40 += ( (p40 >> DAT_BITS) << PAD_BITS );
		p50 += ( (p50 >> DAT_BITS) << PAD_BITS );
		p60 += ( (p60 >> DAT_BITS) << PAD_BITS );
		p70 += ( (p70 >> DAT_BITS) << PAD_BITS );
		p80 += ( (p80 >> DAT_BITS) << PAD_BITS );
		p90 += ( (p90 >> DAT_BITS) << PAD_BITS );
		pa0 += ( (pa0 >> DAT_BITS) << PAD_BITS );
		pb0 += ( (pb0 >> DAT_BITS) << PAD_BITS );
		pc0 += ( (pc0 >> DAT_BITS) << PAD_BITS );
		pd0 += ( (pd0 >> DAT_BITS) << PAD_BITS );
		pe0 += ( (pe0 >> DAT_BITS) << PAD_BITS );
		pf0 += ( (pf0 >> DAT_BITS) << PAD_BITS );
		p100 += ( (p100 >> DAT_BITS) << PAD_BITS );
		p110 += ( (p110 >> DAT_BITS) << PAD_BITS );
		p120 += ( (p120 >> DAT_BITS) << PAD_BITS );
		p130 += ( (p130 >> DAT_BITS) << PAD_BITS );
		p140 += ( (p140 >> DAT_BITS) << PAD_BITS );
		p150 += ( (p150 >> DAT_BITS) << PAD_BITS );
		p160 += ( (p160 >> DAT_BITS) << PAD_BITS );
		p170 += ( (p170 >> DAT_BITS) << PAD_BITS );
		p180 += ( (p180 >> DAT_BITS) << PAD_BITS );
		p190 += ( (p190 >> DAT_BITS) << PAD_BITS );
		p1a0 += ( (p1a0 >> DAT_BITS) << PAD_BITS );
		p1b0 += ( (p1b0 >> DAT_BITS) << PAD_BITS );
		p1c0 += ( (p1c0 >> DAT_BITS) << PAD_BITS );
		p1d0 += ( (p1d0 >> DAT_BITS) << PAD_BITS );
		p1e0 += ( (p1e0 >> DAT_BITS) << PAD_BITS );
		p1f0 += ( (p1f0 >> DAT_BITS) << PAD_BITS );
		p200 += ( (p200 >> DAT_BITS) << PAD_BITS );
		p210 += ( (p210 >> DAT_BITS) << PAD_BITS );
		p220 += ( (p220 >> DAT_BITS) << PAD_BITS );
		p230 += ( (p230 >> DAT_BITS) << PAD_BITS );
		p240 += ( (p240 >> DAT_BITS) << PAD_BITS );
		p250 += ( (p250 >> DAT_BITS) << PAD_BITS );
		p260 += ( (p260 >> DAT_BITS) << PAD_BITS );
		p270 += ( (p270 >> DAT_BITS) << PAD_BITS );
		p280 += ( (p280 >> DAT_BITS) << PAD_BITS );
		p290 += ( (p290 >> DAT_BITS) << PAD_BITS );
		p2a0 += ( (p2a0 >> DAT_BITS) << PAD_BITS );
		p2b0 += ( (p2b0 >> DAT_BITS) << PAD_BITS );
		p2c0 += ( (p2c0 >> DAT_BITS) << PAD_BITS );
		p2d0 += ( (p2d0 >> DAT_BITS) << PAD_BITS );
		p2e0 += ( (p2e0 >> DAT_BITS) << PAD_BITS );
		p2f0 += ( (p2f0 >> DAT_BITS) << PAD_BITS );
	#if !defined(USE_SSE2) && !defined(MULTITHREAD)
		p0123[0] = 0; p0123[1] = p1; p0123[2] = p2; p0123[3] = p3;
	#endif
		poff[     0] =   0; poff[     1] =     p4; poff[     2] =     p8; poff[     3] =     pc;
		poff[0x04+0] = p10; poff[0x04+1] = p10+p4; poff[0x04+2] = p10+p8; poff[0x04+3] = p10+pc;
		poff[0x08+0] = p20; poff[0x08+1] = p20+p4; poff[0x08+2] = p20+p8; poff[0x08+3] = p20+pc;
		poff[0x0c+0] = p30; poff[0x0c+1] = p30+p4; poff[0x0c+2] = p30+p8; poff[0x0c+3] = p30+pc;
		poff[0x10+0] = p40; poff[0x10+1] = p40+p4; poff[0x10+2] = p40+p8; poff[0x10+3] = p40+pc;
		poff[0x14+0] = p50; poff[0x14+1] = p50+p4; poff[0x14+2] = p50+p8; poff[0x14+3] = p50+pc;
		poff[0x18+0] = p60; poff[0x18+1] = p60+p4; poff[0x18+2] = p60+p8; poff[0x18+3] = p60+pc;
		poff[0x1c+0] = p70; poff[0x1c+1] = p70+p4; poff[0x1c+2] = p70+p8; poff[0x1c+3] = p70+pc;
		poff[0x20+0] = p80; poff[0x20+1] = p80+p4; poff[0x20+2] = p80+p8; poff[0x20+3] = p80+pc;
		poff[0x24+0] = p90; poff[0x24+1] = p90+p4; poff[0x24+2] = p90+p8; poff[0x24+3] = p90+pc;
		poff[0x28+0] = pa0; poff[0x28+1] = pa0+p4; poff[0x28+2] = pa0+p8; poff[0x28+3] = pa0+pc;
		poff[0x2c+0] = pb0; poff[0x2c+1] = pb0+p4; poff[0x2c+2] = pb0+p8; poff[0x2c+3] = pb0+pc;
		poff[0x30+0] = pc0; poff[0x30+1] = pc0+p4; poff[0x30+2] = pc0+p8; poff[0x30+3] = pc0+pc;
		poff[0x34+0] = pd0; poff[0x34+1] = pd0+p4; poff[0x34+2] = pd0+p8; poff[0x34+3] = pd0+pc;
		poff[0x38+0] = pe0; poff[0x38+1] = pe0+p4; poff[0x38+2] = pe0+p8; poff[0x38+3] = pe0+pc;
		poff[0x3c+0] = pf0; poff[0x3c+1] = pf0+p4; poff[0x3c+2] = pf0+p8; poff[0x3c+3] = pf0+pc;
		for(l = 0; l < 64; l++) {
			poff[ 64+l] = poff[l] + p100;
			poff[128+l] = poff[l] + p200;
		}

	#ifndef MULTITHREAD
	// Cf. radix768_dit_pass1() for details on the indexing scheme used here:
	// Set dit_offsets_lo for 1st set of radix-256 DIT inputs:
		dit_offsets_lo[0x00] =  0;	dit_offsets_lo[0x10] = pf;
		dit_offsets_lo[0x01] = p1;	dit_offsets_lo[0x11] = pe;
		dit_offsets_lo[0x02] = p3;	dit_offsets_lo[0x12] = pd;
		dit_offsets_lo[0x03] = p2;	dit_offsets_lo[0x13] = pc;
		dit_offsets_lo[0x04] = p7;	dit_offsets_lo[0x14] = pb;
		dit_offsets_lo[0x05] = p6;	dit_offsets_lo[0x15] = pa;
		dit_offsets_lo[0x06] = p5;	dit_offsets_lo[0x16] = p9;
		dit_offsets_lo[0x07] = p4;	dit_offsets_lo[0x17] = p8;
		dit_offsets_lo[0x08] = pf;	dit_offsets_lo[0x18] = p7;
		dit_offsets_lo[0x09] = pe;	dit_offsets_lo[0x19] = p6;
		dit_offsets_lo[0x0a] = pd;	dit_offsets_lo[0x1a] = p5;
		dit_offsets_lo[0x0b] = pc;	dit_offsets_lo[0x1b] = p4;
		dit_offsets_lo[0x0c] = pb;	dit_offsets_lo[0x1c] = p3;
		dit_offsets_lo[0x0d] = pa;	dit_offsets_lo[0x1d] = p2;
		dit_offsets_lo[0x0e] = p9;	dit_offsets_lo[0x1e] = p1;
		dit_offsets_lo[0x0f] = p8;	dit_offsets_lo[0x1f] =  0;

	// Set dit_offsets for 2nd set of radix-256 DIT inputs:
		dit_offsets_lo[0x20] = p5;	dit_offsets_lo[0x30] = p9;
		dit_offsets_lo[0x21] = p4;	dit_offsets_lo[0x31] = p8;
		dit_offsets_lo[0x22] = p6;	dit_offsets_lo[0x32] = pa;
		dit_offsets_lo[0x23] = p7;	dit_offsets_lo[0x33] = pb;
		dit_offsets_lo[0x24] = p1;	dit_offsets_lo[0x34] = pe;
		dit_offsets_lo[0x25] =  0;	dit_offsets_lo[0x35] = pf;
		dit_offsets_lo[0x26] = p2;	dit_offsets_lo[0x36] = pc;
		dit_offsets_lo[0x27] = p3;	dit_offsets_lo[0x37] = pd;
		dit_offsets_lo[0x28] = p9;	dit_offsets_lo[0x38] = p1;
		dit_offsets_lo[0x29] = p8;	dit_offsets_lo[0x39] =  0;
		dit_offsets_lo[0x2a] = pa;	dit_offsets_lo[0x3a] = p2;
		dit_offsets_lo[0x2b] = pb;	dit_offsets_lo[0x3b] = p3;
		dit_offsets_lo[0x2c] = pe;	dit_offsets_lo[0x3c] = p6;
		dit_offsets_lo[0x2d] = pf;	dit_offsets_lo[0x3d] = p7;
		dit_offsets_lo[0x2e] = pc;	dit_offsets_lo[0x3e] = p4;
		dit_offsets_lo[0x2f] = pd;	dit_offsets_lo[0x3f] = p5;

	// Set dit_offsets for 3rd set of radix-256 DIT inputs:
		dit_offsets_lo[0x40] = pa;	dit_offsets_lo[0x50] = p2;
		dit_offsets_lo[0x41] = pb;	dit_offsets_lo[0x51] = p3;
		dit_offsets_lo[0x42] = p8;	dit_offsets_lo[0x52] =  0;
		dit_offsets_lo[0x43] = p9;	dit_offsets_lo[0x53] = p1;
		dit_offsets_lo[0x44] = pc;	dit_offsets_lo[0x54] = p4;
		dit_offsets_lo[0x45] = pd;	dit_offsets_lo[0x55] = p5;
		dit_offsets_lo[0x46] = pf;	dit_offsets_lo[0x56] = p7;
		dit_offsets_lo[0x47] = pe;	dit_offsets_lo[0x57] = p6;
		dit_offsets_lo[0x48] = p2;	dit_offsets_lo[0x58] = pc;
		dit_offsets_lo[0x49] = p3;	dit_offsets_lo[0x59] = pd;
		dit_offsets_lo[0x4a] =  0;	dit_offsets_lo[0x5a] = pf;
		dit_offsets_lo[0x4b] = p1;	dit_offsets_lo[0x5b] = pe;
		dit_offsets_lo[0x4c] = p4;	dit_offsets_lo[0x5c] = p8;
		dit_offsets_lo[0x4d] = p5;	dit_offsets_lo[0x5d] = p9;
		dit_offsets_lo[0x4e] = p7;	dit_offsets_lo[0x5e] = pb;
		dit_offsets_lo[0x4f] = p6;	dit_offsets_lo[0x5f] = pa;

	// ...and one distinct high-part vector for each radix-256 DIT:
		dit_offsets_hi1[0x0] =   0;	dit_offsets_hi2[0x0] = p50;	dit_offsets_hi3[0x0] = pa0;
		dit_offsets_hi1[0x1] = p10;	dit_offsets_hi2[0x1] = p40;	dit_offsets_hi3[0x1] = pb0;
		dit_offsets_hi1[0x2] = p30;	dit_offsets_hi2[0x2] = p60;	dit_offsets_hi3[0x2] = p80;
		dit_offsets_hi1[0x3] = p20;	dit_offsets_hi2[0x3] = p70;	dit_offsets_hi3[0x3] = p90;
		dit_offsets_hi1[0x4] = p70;	dit_offsets_hi2[0x4] = p10;	dit_offsets_hi3[0x4] = pc0;
		dit_offsets_hi1[0x5] = p60;	dit_offsets_hi2[0x5] =   0;	dit_offsets_hi3[0x5] = pd0;
		dit_offsets_hi1[0x6] = p50;	dit_offsets_hi2[0x6] = p20;	dit_offsets_hi3[0x6] = pf0;
		dit_offsets_hi1[0x7] = p40;	dit_offsets_hi2[0x7] = p30;	dit_offsets_hi3[0x7] = pe0;
		dit_offsets_hi1[0x8] = pf0;	dit_offsets_hi2[0x8] = p90;	dit_offsets_hi3[0x8] = p20;
		dit_offsets_hi1[0x9] = pe0;	dit_offsets_hi2[0x9] = p80;	dit_offsets_hi3[0x9] = p30;
		dit_offsets_hi1[0xa] = pd0;	dit_offsets_hi2[0xa] = pa0;	dit_offsets_hi3[0xa] =   0;
		dit_offsets_hi1[0xb] = pc0;	dit_offsets_hi2[0xb] = pb0;	dit_offsets_hi3[0xb] = p10;
		dit_offsets_hi1[0xc] = pb0;	dit_offsets_hi2[0xc] = pe0;	dit_offsets_hi3[0xc] = p40;
		dit_offsets_hi1[0xd] = pa0;	dit_offsets_hi2[0xd] = pf0;	dit_offsets_hi3[0xd] = p50;
		dit_offsets_hi1[0xe] = p90;	dit_offsets_hi2[0xe] = pc0;	dit_offsets_hi3[0xe] = p70;
		dit_offsets_hi1[0xf] = p80;	dit_offsets_hi2[0xf] = pd0;	dit_offsets_hi3[0xf] = p60;

	  #ifndef USE_SSE2
		// Idx offsets w.r.to r-array are const and shared by both sets of radix-256 transforms:
		// low parts in first 16 slots, high parts in next 16:
		t_offsets_lo[0x0] = 0x00<<1;	t_offsets_hi[0x0] = 0x00<<1;
		t_offsets_lo[0x1] = 0x01<<1;	t_offsets_hi[0x1] = 0x10<<1;
		t_offsets_lo[0x2] = 0x02<<1;	t_offsets_hi[0x2] = 0x20<<1;
		t_offsets_lo[0x3] = 0x03<<1;	t_offsets_hi[0x3] = 0x30<<1;
		t_offsets_lo[0x4] = 0x04<<1;	t_offsets_hi[0x4] = 0x40<<1;
		t_offsets_lo[0x5] = 0x05<<1;	t_offsets_hi[0x5] = 0x50<<1;
		t_offsets_lo[0x6] = 0x06<<1;	t_offsets_hi[0x6] = 0x60<<1;
		t_offsets_lo[0x7] = 0x07<<1;	t_offsets_hi[0x7] = 0x70<<1;
		t_offsets_lo[0x8] = 0x08<<1;	t_offsets_hi[0x8] = 0x80<<1;
		t_offsets_lo[0x9] = 0x09<<1;	t_offsets_hi[0x9] = 0x90<<1;
		t_offsets_lo[0xa] = 0x0a<<1;	t_offsets_hi[0xa] = 0xa0<<1;
		t_offsets_lo[0xb] = 0x0b<<1;	t_offsets_hi[0xb] = 0xb0<<1;
		t_offsets_lo[0xc] = 0x0c<<1;	t_offsets_hi[0xc] = 0xc0<<1;
		t_offsets_lo[0xd] = 0x0d<<1;	t_offsets_hi[0xd] = 0xd0<<1;
		t_offsets_lo[0xe] = 0x0e<<1;	t_offsets_hi[0xe] = 0xe0<<1;
		t_offsets_lo[0xf] = 0x0f<<1;	t_offsets_hi[0xf] = 0xf0<<1;
	  #endif

	// Cf. radix768_dif_pass1() for details on the indexing scheme used here:
	// Set dit_offsets_lo for the 4 subvectors shared by the 3 sets of radix-256 DIT inputs:
		dif_offsets_lo[0x00] =  0;	dif_offsets_lo[0x10] = p5;	dif_offsets_lo[0x20] = pa;	dif_offsets_lo[0x30] = pf;
		dif_offsets_lo[0x01] = p1;	dif_offsets_lo[0x11] = p4;	dif_offsets_lo[0x21] = pb;	dif_offsets_lo[0x31] = pe;
		dif_offsets_lo[0x02] = p2;	dif_offsets_lo[0x12] = p7;	dif_offsets_lo[0x22] = p9;	dif_offsets_lo[0x32] = pc;
		dif_offsets_lo[0x03] = p3;	dif_offsets_lo[0x13] = p6;	dif_offsets_lo[0x23] = p8;	dif_offsets_lo[0x33] = pd;
		dif_offsets_lo[0x04] = p5;	dif_offsets_lo[0x14] = p2;	dif_offsets_lo[0x24] = pf;	dif_offsets_lo[0x34] = p9;
		dif_offsets_lo[0x05] = p4;	dif_offsets_lo[0x15] = p3;	dif_offsets_lo[0x25] = pe;	dif_offsets_lo[0x35] = p8;
		dif_offsets_lo[0x06] = p7;	dif_offsets_lo[0x16] = p1;	dif_offsets_lo[0x26] = pc;	dif_offsets_lo[0x36] = pb;
		dif_offsets_lo[0x07] = p6;	dif_offsets_lo[0x17] =  0;	dif_offsets_lo[0x27] = pd;	dif_offsets_lo[0x37] = pa;
		dif_offsets_lo[0x08] = pa;	dif_offsets_lo[0x18] = pf;	dif_offsets_lo[0x28] = p5;	dif_offsets_lo[0x38] = p2;
		dif_offsets_lo[0x09] = pb;	dif_offsets_lo[0x19] = pe;	dif_offsets_lo[0x29] = p4;	dif_offsets_lo[0x39] = p3;
		dif_offsets_lo[0x0a] = p9;	dif_offsets_lo[0x1a] = pc;	dif_offsets_lo[0x2a] = p7;	dif_offsets_lo[0x3a] = p1;
		dif_offsets_lo[0x0b] = p8;	dif_offsets_lo[0x1b] = pd;	dif_offsets_lo[0x2b] = p6;	dif_offsets_lo[0x3b] =  0;
		dif_offsets_lo[0x0c] = pf;	dif_offsets_lo[0x1c] = p9;	dif_offsets_lo[0x2c] = p2;	dif_offsets_lo[0x3c] = p7;
		dif_offsets_lo[0x0d] = pe;	dif_offsets_lo[0x1d] = p8;	dif_offsets_lo[0x2d] = p3;	dif_offsets_lo[0x3d] = p6;
		dif_offsets_lo[0x0e] = pc;	dif_offsets_lo[0x1e] = pb;	dif_offsets_lo[0x2e] = p1;	dif_offsets_lo[0x3e] = p4;
		dif_offsets_lo[0x0f] = pd;	dif_offsets_lo[0x1f] = pa;	dif_offsets_lo[0x2f] =  0;	dif_offsets_lo[0x3f] = p5;

	// ...and one distinct high-part vector for each radix-256 DFT:
		dif_offsets_hi1[0x0] =   0;	dif_offsets_hi2[0x0] = p50;	dif_offsets_hi3[0x0] = pa0;
		dif_offsets_hi1[0x1] = p10;	dif_offsets_hi2[0x1] = p40;	dif_offsets_hi3[0x1] = pb0;
		dif_offsets_hi1[0x2] = p20;	dif_offsets_hi2[0x2] = p70;	dif_offsets_hi3[0x2] = p90;
		dif_offsets_hi1[0x3] = p30;	dif_offsets_hi2[0x3] = p60;	dif_offsets_hi3[0x3] = p80;
		dif_offsets_hi1[0x4] = p50;	dif_offsets_hi2[0x4] = p20;	dif_offsets_hi3[0x4] = pf0;
		dif_offsets_hi1[0x5] = p40;	dif_offsets_hi2[0x5] = p30;	dif_offsets_hi3[0x5] = pe0;
		dif_offsets_hi1[0x6] = p70;	dif_offsets_hi2[0x6] = p10;	dif_offsets_hi3[0x6] = pc0;
		dif_offsets_hi1[0x7] = p60;	dif_offsets_hi2[0x7] =   0;	dif_offsets_hi3[0x7] = pd0;
		dif_offsets_hi1[0x8] = pa0;	dif_offsets_hi2[0x8] = pf0;	dif_offsets_hi3[0x8] = p50;
		dif_offsets_hi1[0x9] = pb0;	dif_offsets_hi2[0x9] = pe0;	dif_offsets_hi3[0x9] = p40;
		dif_offsets_hi1[0xa] = p90;	dif_offsets_hi2[0xa] = pc0;	dif_offsets_hi3[0xa] = p70;
		dif_offsets_hi1[0xb] = p80;	dif_offsets_hi2[0xb] = pd0;	dif_offsets_hi3[0xb] = p60;
		dif_offsets_hi1[0xc] = pf0;	dif_offsets_hi2[0xc] = p90;	dif_offsets_hi3[0xc] = p20;
		dif_offsets_hi1[0xd] = pe0;	dif_offsets_hi2[0xd] = p80;	dif_offsets_hi3[0xd] = p30;
		dif_offsets_hi1[0xe] = pc0;	dif_offsets_hi2[0xe] = pb0;	dif_offsets_hi3[0xe] = p10;
		dif_offsets_hi1[0xf] = pd0;	dif_offsets_hi2[0xf] = pa0;	dif_offsets_hi3[0xf] =   0;

	// DIF Index-high-bits triplets needed for compact-obj-code scheme:
	  #ifdef USE_SSE2
		k = 0;	// In SIMD mode these are 0xtr-offsets w.r.to a local store:
		dif_triplets[k] = 0x2f0; dif_triplets[k+1] = 0x1f0; dif_triplets[k+2] = 0x0f0; k += 3;
		dif_triplets[k] = 0x2e0; dif_triplets[k+1] = 0x1e0; dif_triplets[k+2] = 0x0e0; k += 3;
		dif_triplets[k] = 0x2d0; dif_triplets[k+1] = 0x1d0; dif_triplets[k+2] = 0x0d0; k += 3;
		dif_triplets[k] = 0x2c0; dif_triplets[k+1] = 0x1c0; dif_triplets[k+2] = 0x0c0; k += 3;
		dif_triplets[k] = 0x2b0; dif_triplets[k+1] = 0x1b0; dif_triplets[k+2] = 0x0b0; k += 3;
		dif_triplets[k] = 0x2a0; dif_triplets[k+1] = 0x1a0; dif_triplets[k+2] = 0x0a0; k += 3;
		dif_triplets[k] = 0x290; dif_triplets[k+1] = 0x190; dif_triplets[k+2] = 0x090; k += 3;
		dif_triplets[k] = 0x280; dif_triplets[k+1] = 0x180; dif_triplets[k+2] = 0x080; k += 3;
		dif_triplets[k] = 0x270; dif_triplets[k+1] = 0x170; dif_triplets[k+2] = 0x070; k += 3;
		dif_triplets[k] = 0x260; dif_triplets[k+1] = 0x160; dif_triplets[k+2] = 0x060; k += 3;
		dif_triplets[k] = 0x250; dif_triplets[k+1] = 0x150; dif_triplets[k+2] = 0x050; k += 3;
		dif_triplets[k] = 0x240; dif_triplets[k+1] = 0x140; dif_triplets[k+2] = 0x040; k += 3;
		dif_triplets[k] = 0x230; dif_triplets[k+1] = 0x130; dif_triplets[k+2] = 0x030; k += 3;
		dif_triplets[k] = 0x220; dif_triplets[k+1] = 0x120; dif_triplets[k+2] = 0x020; k += 3;
		dif_triplets[k] = 0x210; dif_triplets[k+1] = 0x110; dif_triplets[k+2] = 0x010; k += 3;
		dif_triplets[k] = 0x200; dif_triplets[k+1] = 0x100; dif_triplets[k+2] = 0x000; k += 3;
		dif_triplets[k] = 0x1f0; dif_triplets[k+1] = 0x0f0; dif_triplets[k+2] = 0x2f0; k += 3;
		dif_triplets[k] = 0x1e0; dif_triplets[k+1] = 0x0e0; dif_triplets[k+2] = 0x2e0; k += 3;
		dif_triplets[k] = 0x1d0; dif_triplets[k+1] = 0x0d0; dif_triplets[k+2] = 0x2d0; k += 3;
		dif_triplets[k] = 0x1c0; dif_triplets[k+1] = 0x0c0; dif_triplets[k+2] = 0x2c0; k += 3;
		dif_triplets[k] = 0x1b0; dif_triplets[k+1] = 0x0b0; dif_triplets[k+2] = 0x2b0; k += 3;
		dif_triplets[k] = 0x1a0; dif_triplets[k+1] = 0x0a0; dif_triplets[k+2] = 0x2a0; k += 3;
		dif_triplets[k] = 0x190; dif_triplets[k+1] = 0x090; dif_triplets[k+2] = 0x290; k += 3;
		dif_triplets[k] = 0x180; dif_triplets[k+1] = 0x080; dif_triplets[k+2] = 0x280; k += 3;
		dif_triplets[k] = 0x170; dif_triplets[k+1] = 0x070; dif_triplets[k+2] = 0x270; k += 3;
		dif_triplets[k] = 0x160; dif_triplets[k+1] = 0x060; dif_triplets[k+2] = 0x260; k += 3;
		dif_triplets[k] = 0x150; dif_triplets[k+1] = 0x050; dif_triplets[k+2] = 0x250; k += 3;
		dif_triplets[k] = 0x140; dif_triplets[k+1] = 0x040; dif_triplets[k+2] = 0x240; k += 3;
		dif_triplets[k] = 0x130; dif_triplets[k+1] = 0x030; dif_triplets[k+2] = 0x230; k += 3;
		dif_triplets[k] = 0x120; dif_triplets[k+1] = 0x020; dif_triplets[k+2] = 0x220; k += 3;
		dif_triplets[k] = 0x110; dif_triplets[k+1] = 0x010; dif_triplets[k+2] = 0x210; k += 3;
		dif_triplets[k] = 0x100; dif_triplets[k+1] = 0x000; dif_triplets[k+2] = 0x200; k += 3;
		dif_triplets[k] = 0x0f0; dif_triplets[k+1] = 0x2f0; dif_triplets[k+2] = 0x1f0; k += 3;
		dif_triplets[k] = 0x0e0; dif_triplets[k+1] = 0x2e0; dif_triplets[k+2] = 0x1e0; k += 3;
		dif_triplets[k] = 0x0d0; dif_triplets[k+1] = 0x2d0; dif_triplets[k+2] = 0x1d0; k += 3;
		dif_triplets[k] = 0x0c0; dif_triplets[k+1] = 0x2c0; dif_triplets[k+2] = 0x1c0; k += 3;
		dif_triplets[k] = 0x0b0; dif_triplets[k+1] = 0x2b0; dif_triplets[k+2] = 0x1b0; k += 3;
		dif_triplets[k] = 0x0a0; dif_triplets[k+1] = 0x2a0; dif_triplets[k+2] = 0x1a0; k += 3;
		dif_triplets[k] = 0x090; dif_triplets[k+1] = 0x290; dif_triplets[k+2] = 0x190; k += 3;
		dif_triplets[k] = 0x080; dif_triplets[k+1] = 0x280; dif_triplets[k+2] = 0x180; k += 3;
		dif_triplets[k] = 0x070; dif_triplets[k+1] = 0x270; dif_triplets[k+2] = 0x170; k += 3;
		dif_triplets[k] = 0x060; dif_triplets[k+1] = 0x260; dif_triplets[k+2] = 0x160; k += 3;
		dif_triplets[k] = 0x050; dif_triplets[k+1] = 0x250; dif_triplets[k+2] = 0x150; k += 3;
		dif_triplets[k] = 0x040; dif_triplets[k+1] = 0x240; dif_triplets[k+2] = 0x140; k += 3;
		dif_triplets[k] = 0x030; dif_triplets[k+1] = 0x230; dif_triplets[k+2] = 0x130; k += 3;
		dif_triplets[k] = 0x020; dif_triplets[k+1] = 0x220; dif_triplets[k+2] = 0x120; k += 3;
		dif_triplets[k] = 0x010; dif_triplets[k+1] = 0x210; dif_triplets[k+2] = 0x110; k += 3;
		dif_triplets[k] = 0x000; dif_triplets[k+1] = 0x200; dif_triplets[k+2] = 0x100;
		// IN SIMD mode need to double all the above to turn from vec_dbl to vec_cmplx ptr offsets:
		for(l = 0; l < 144; l++) {
			dif_triplets[l] <<= 1;
		}
	  #else
		k = 0;
		dif_triplets[k] = p2f0; dif_triplets[k+1] = p1f0; dif_triplets[k+2] = pf0; k += 3;
		dif_triplets[k] = p2e0; dif_triplets[k+1] = p1e0; dif_triplets[k+2] = pe0; k += 3;
		dif_triplets[k] = p2d0; dif_triplets[k+1] = p1d0; dif_triplets[k+2] = pd0; k += 3;
		dif_triplets[k] = p2c0; dif_triplets[k+1] = p1c0; dif_triplets[k+2] = pc0; k += 3;
		dif_triplets[k] = p2b0; dif_triplets[k+1] = p1b0; dif_triplets[k+2] = pb0; k += 3;
		dif_triplets[k] = p2a0; dif_triplets[k+1] = p1a0; dif_triplets[k+2] = pa0; k += 3;
		dif_triplets[k] = p290; dif_triplets[k+1] = p190; dif_triplets[k+2] = p90; k += 3;
		dif_triplets[k] = p280; dif_triplets[k+1] = p180; dif_triplets[k+2] = p80; k += 3;
		dif_triplets[k] = p270; dif_triplets[k+1] = p170; dif_triplets[k+2] = p70; k += 3;
		dif_triplets[k] = p260; dif_triplets[k+1] = p160; dif_triplets[k+2] = p60; k += 3;
		dif_triplets[k] = p250; dif_triplets[k+1] = p150; dif_triplets[k+2] = p50; k += 3;
		dif_triplets[k] = p240; dif_triplets[k+1] = p140; dif_triplets[k+2] = p40; k += 3;
		dif_triplets[k] = p230; dif_triplets[k+1] = p130; dif_triplets[k+2] = p30; k += 3;
		dif_triplets[k] = p220; dif_triplets[k+1] = p120; dif_triplets[k+2] = p20; k += 3;
		dif_triplets[k] = p210; dif_triplets[k+1] = p110; dif_triplets[k+2] = p10; k += 3;
		dif_triplets[k] = p200; dif_triplets[k+1] = p100; dif_triplets[k+2] =   0; k += 3;
		dif_triplets[k] = p1f0; dif_triplets[k+1] = pf0; dif_triplets[k+2] = p2f0; k += 3;
		dif_triplets[k] = p1e0; dif_triplets[k+1] = pe0; dif_triplets[k+2] = p2e0; k += 3;
		dif_triplets[k] = p1d0; dif_triplets[k+1] = pd0; dif_triplets[k+2] = p2d0; k += 3;
		dif_triplets[k] = p1c0; dif_triplets[k+1] = pc0; dif_triplets[k+2] = p2c0; k += 3;
		dif_triplets[k] = p1b0; dif_triplets[k+1] = pb0; dif_triplets[k+2] = p2b0; k += 3;
		dif_triplets[k] = p1a0; dif_triplets[k+1] = pa0; dif_triplets[k+2] = p2a0; k += 3;
		dif_triplets[k] = p190; dif_triplets[k+1] = p90; dif_triplets[k+2] = p290; k += 3;
		dif_triplets[k] = p180; dif_triplets[k+1] = p80; dif_triplets[k+2] = p280; k += 3;
		dif_triplets[k] = p170; dif_triplets[k+1] = p70; dif_triplets[k+2] = p270; k += 3;
		dif_triplets[k] = p160; dif_triplets[k+1] = p60; dif_triplets[k+2] = p260; k += 3;
		dif_triplets[k] = p150; dif_triplets[k+1] = p50; dif_triplets[k+2] = p250; k += 3;
		dif_triplets[k] = p140; dif_triplets[k+1] = p40; dif_triplets[k+2] = p240; k += 3;
		dif_triplets[k] = p130; dif_triplets[k+1] = p30; dif_triplets[k+2] = p230; k += 3;
		dif_triplets[k] = p120; dif_triplets[k+1] = p20; dif_triplets[k+2] = p220; k += 3;
		dif_triplets[k] = p110; dif_triplets[k+1] = p10; dif_triplets[k+2] = p210; k += 3;
		dif_triplets[k] = p100; dif_triplets[k+1] =   0; dif_triplets[k+2] = p200; k += 3;
		dif_triplets[k] = pf0; dif_triplets[k+1] = p2f0; dif_triplets[k+2] = p1f0; k += 3;
		dif_triplets[k] = pe0; dif_triplets[k+1] = p2e0; dif_triplets[k+2] = p1e0; k += 3;
		dif_triplets[k] = pd0; dif_triplets[k+1] = p2d0; dif_triplets[k+2] = p1d0; k += 3;
		dif_triplets[k] = pc0; dif_triplets[k+1] = p2c0; dif_triplets[k+2] = p1c0; k += 3;
		dif_triplets[k] = pb0; dif_triplets[k+1] = p2b0; dif_triplets[k+2] = p1b0; k += 3;
		dif_triplets[k] = pa0; dif_triplets[k+1] = p2a0; dif_triplets[k+2] = p1a0; k += 3;
		dif_triplets[k] = p90; dif_triplets[k+1] = p290; dif_triplets[k+2] = p190; k += 3;
		dif_triplets[k] = p80; dif_triplets[k+1] = p280; dif_triplets[k+2] = p180; k += 3;
		dif_triplets[k] = p70; dif_triplets[k+1] = p270; dif_triplets[k+2] = p170; k += 3;
		dif_triplets[k] = p60; dif_triplets[k+1] = p260; dif_triplets[k+2] = p160; k += 3;
		dif_triplets[k] = p50; dif_triplets[k+1] = p250; dif_triplets[k+2] = p150; k += 3;
		dif_triplets[k] = p40; dif_triplets[k+1] = p240; dif_triplets[k+2] = p140; k += 3;
		dif_triplets[k] = p30; dif_triplets[k+1] = p230; dif_triplets[k+2] = p130; k += 3;
		dif_triplets[k] = p20; dif_triplets[k+1] = p220; dif_triplets[k+2] = p120; k += 3;
		dif_triplets[k] = p10; dif_triplets[k+1] = p210; dif_triplets[k+2] = p110; k += 3;
		dif_triplets[k] =   0; dif_triplets[k+1] = p200; dif_triplets[k+2] = p100;
	  #endif
	// DIT Index-high-bits triplets needed for compact-obj-code scheme:
	  #ifdef USE_SSE2
		k = 0;	// In SIMD mode these are ptr-offsets w.r.to a local store:
		dit_triplets[k] = 0x0f0; dit_triplets[k+1] = 0x1f0; dit_triplets[k+2] = 0x2f0; k += 3;
		dit_triplets[k] = 0x1e0; dit_triplets[k+1] = 0x2e0; dit_triplets[k+2] = 0x0e0; k += 3;
		dit_triplets[k] = 0x2d0; dit_triplets[k+1] = 0x0d0; dit_triplets[k+2] = 0x1d0; k += 3;
		dit_triplets[k] = 0x0c0; dit_triplets[k+1] = 0x1c0; dit_triplets[k+2] = 0x2c0; k += 3;
		dit_triplets[k] = 0x1b0; dit_triplets[k+1] = 0x2b0; dit_triplets[k+2] = 0x0b0; k += 3;
		dit_triplets[k] = 0x2a0; dit_triplets[k+1] = 0x0a0; dit_triplets[k+2] = 0x1a0; k += 3;
		dit_triplets[k] = 0x090; dit_triplets[k+1] = 0x190; dit_triplets[k+2] = 0x290; k += 3;
		dit_triplets[k] = 0x180; dit_triplets[k+1] = 0x280; dit_triplets[k+2] = 0x080; k += 3;
		dit_triplets[k] = 0x270; dit_triplets[k+1] = 0x070; dit_triplets[k+2] = 0x170; k += 3;
		dit_triplets[k] = 0x060; dit_triplets[k+1] = 0x160; dit_triplets[k+2] = 0x260; k += 3;
		dit_triplets[k] = 0x150; dit_triplets[k+1] = 0x250; dit_triplets[k+2] = 0x050; k += 3;
		dit_triplets[k] = 0x240; dit_triplets[k+1] = 0x040; dit_triplets[k+2] = 0x140; k += 3;
		dit_triplets[k] = 0x030; dit_triplets[k+1] = 0x130; dit_triplets[k+2] = 0x230; k += 3;
		dit_triplets[k] = 0x120; dit_triplets[k+1] = 0x220; dit_triplets[k+2] = 0x020; k += 3;
		dit_triplets[k] = 0x210; dit_triplets[k+1] = 0x010; dit_triplets[k+2] = 0x110; k += 3;
		dit_triplets[k] = 0x000; dit_triplets[k+1] = 0x100; dit_triplets[k+2] = 0x200;
		// IN SIMD mode need to double all the above to turn from vec_dbl to vec_cmplx ptr offsets:
		for(l = 0; l < 48; l++) {
			dit_triplets[l] <<= 1;
		}
	  #else
		k = 0;
		dit_triplets[k] =  pf0; dit_triplets[k+1] = p1f0; dit_triplets[k+2] = p2f0; k += 3;
		dit_triplets[k] = p1e0; dit_triplets[k+1] = p2e0; dit_triplets[k+2] =  pe0; k += 3;
		dit_triplets[k] = p2d0; dit_triplets[k+1] =  pd0; dit_triplets[k+2] = p1d0; k += 3;
		dit_triplets[k] =  pc0; dit_triplets[k+1] = p1c0; dit_triplets[k+2] = p2c0; k += 3;
		dit_triplets[k] = p1b0; dit_triplets[k+1] = p2b0; dit_triplets[k+2] =  pb0; k += 3;
		dit_triplets[k] = p2a0; dit_triplets[k+1] =  pa0; dit_triplets[k+2] = p1a0; k += 3;
		dit_triplets[k] =  p90; dit_triplets[k+1] = p190; dit_triplets[k+2] = p290; k += 3;
		dit_triplets[k] = p180; dit_triplets[k+1] = p280; dit_triplets[k+2] =  p80; k += 3;
		dit_triplets[k] = p270; dit_triplets[k+1] =  p70; dit_triplets[k+2] = p170; k += 3;
		dit_triplets[k] =  p60; dit_triplets[k+1] = p160; dit_triplets[k+2] = p260; k += 3;
		dit_triplets[k] = p150; dit_triplets[k+1] = p250; dit_triplets[k+2] =  p50; k += 3;
		dit_triplets[k] = p240; dit_triplets[k+1] =  p40; dit_triplets[k+2] = p140; k += 3;
		dit_triplets[k] =  p30; dit_triplets[k+1] = p130; dit_triplets[k+2] = p230; k += 3;
		dit_triplets[k] = p120; dit_triplets[k+1] = p220; dit_triplets[k+2] =  p20; k += 3;
		dit_triplets[k] = p210; dit_triplets[k+1] =  p10; dit_triplets[k+2] = p110; k += 3;
		dit_triplets[k] =    0; dit_triplets[k+1] = p100; dit_triplets[k+2] = p200;
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

/*...The radix-768 final DIT pass is here.	*/

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
		#include "radix768_main_carry_loop.h"

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
		ASSERT(0x0 == cy768_process_chunk( (void*)(&tdat[j + pool_work_units]) ), "Main-thread task failure!");
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

void radix768_dif_pass1(double a[], int n)
{
/*
!...Acronym: DIF = Decimation In Frequency
!
!...Subroutine to perform an initial radix-768 complex DIF FFT pass on the data in the length-N real vector A.
!
!   See the documentation in radix16_dif_pass for further details on storage and indexing.
!   See the documentation in radix3_dif_pass for details on the radix-3 subtransforms.
!   See the documentation in radix768_dif_pass for details on the twiddleless coprime-radix-DFTS setup.
*/
	int j,j1,j2,jp,jt;
	int k,l,l1,l2,k0,k1,k2;
	static int dif_triplets[144];
	// In order to preserve length-2 numeric-index-offset property here, use hex for digit:
	static int NDIVR, first_entry=TRUE,
		p1,p2,p3,p4,p5,p6,p7,p8,p9,pa,pb,pc,pd,pe,pf,
		     p10,p20,p30,p40,p50,p60,p70,p80,p90,pa0,pb0,pc0,pd0,pe0,pf0,
		p100,p110,p120,p130,p140,p150,p160,p170,p180,p190,p1a0,p1b0,p1c0,p1d0,p1e0,p1f0,
		p200,p210,p220,p230,p240,p250,p260,p270,p280,p290,p2a0,p2b0,p2c0,p2d0,p2e0,p2f0;
	const double c3m1 = -1.50000000000000000000, s = 0.86602540378443864675;	// cos(twopi/3)-1, sin(twopi/3)
	static int i_offsets_lo[16], i_offsets_hi[16];
	static int o_offsets_lo[64];	// 4 subsets of 16
	// Bitfields encoding the sequence of the o_offsets_lo subset0-3 vectors to use for each radix-256 DFT's outputs:
	const uint32 o_idx1 = 0xe79e79e4,o_idx2 = 0x79e79e79,o_idx3 = 0x9e79e79e;
	static int o_offsets_hi1[16],o_offsets_hi2[16],o_offsets_hi3[16];
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
		p1 = NDIVR;		p10 = NDIVR<<4;		p100 = NDIVR<<8;		p200 = NDIVR<<9;
		p2 = p1 + p1;	p20 = p10 + p10;	p110 = p100 + p10;		p210 = p200 + p10;
		p3 = p2 + p1;	p30 = p20 + p10;	p120 = p110 + p10;		p220 = p210 + p10;
		p4 = p3 + p1;	p40 = p30 + p10;	p130 = p120 + p10;		p230 = p220 + p10;
		p5 = p4 + p1;	p50 = p40 + p10;	p140 = p130 + p10;		p240 = p230 + p10;
		p6 = p5 + p1;	p60 = p50 + p10;	p150 = p140 + p10;		p250 = p240 + p10;
		p7 = p6 + p1;	p70 = p60 + p10;	p160 = p150 + p10;		p260 = p250 + p10;
		p8 = p7 + p1;	p80 = p70 + p10;	p170 = p160 + p10;		p270 = p260 + p10;
		p9 = p8 + p1;	p90 = p80 + p10;	p180 = p170 + p10;		p280 = p270 + p10;
		pa = p9 + p1;	pa0 = p90 + p10;	p190 = p180 + p10;		p290 = p280 + p10;
		pb = pa + p1;	pb0 = pa0 + p10;	p1a0 = p190 + p10;		p2a0 = p290 + p10;
		pc = pb + p1;	pc0 = pb0 + p10;	p1b0 = p1a0 + p10;		p2b0 = p2a0 + p10;
		pd = pc + p1;	pd0 = pc0 + p10;	p1c0 = p1b0 + p10;		p2c0 = p2b0 + p10;
		pe = pd + p1;	pe0 = pd0 + p10;	p1d0 = p1c0 + p10;		p2d0 = p2c0 + p10;
		pf = pe + p1;	pf0 = pe0 + p10;	p1e0 = p1d0 + p10;		p2e0 = p2d0 + p10;
											p1f0 = p1e0 + p10;		p2f0 = p2e0 + p10;
		p1 += ( (p1 >> DAT_BITS) << PAD_BITS );
		p2 += ( (p2 >> DAT_BITS) << PAD_BITS );
		p3 += ( (p3 >> DAT_BITS) << PAD_BITS );
		p4 += ( (p4 >> DAT_BITS) << PAD_BITS );
		p5 += ( (p5 >> DAT_BITS) << PAD_BITS );
		p6 += ( (p6 >> DAT_BITS) << PAD_BITS );
		p7 += ( (p7 >> DAT_BITS) << PAD_BITS );
		p8 += ( (p8 >> DAT_BITS) << PAD_BITS );
		p9 += ( (p9 >> DAT_BITS) << PAD_BITS );
		pa += ( (pa >> DAT_BITS) << PAD_BITS );
		pb += ( (pb >> DAT_BITS) << PAD_BITS );
		pc += ( (pc >> DAT_BITS) << PAD_BITS );
		pd += ( (pd >> DAT_BITS) << PAD_BITS );
		pe += ( (pe >> DAT_BITS) << PAD_BITS );
		pf += ( (pf >> DAT_BITS) << PAD_BITS );
		p10 += ( (p10 >> DAT_BITS) << PAD_BITS );
		p20 += ( (p20 >> DAT_BITS) << PAD_BITS );
		p30 += ( (p30 >> DAT_BITS) << PAD_BITS );
		p40 += ( (p40 >> DAT_BITS) << PAD_BITS );
		p50 += ( (p50 >> DAT_BITS) << PAD_BITS );
		p60 += ( (p60 >> DAT_BITS) << PAD_BITS );
		p70 += ( (p70 >> DAT_BITS) << PAD_BITS );
		p80 += ( (p80 >> DAT_BITS) << PAD_BITS );
		p90 += ( (p90 >> DAT_BITS) << PAD_BITS );
		pa0 += ( (pa0 >> DAT_BITS) << PAD_BITS );
		pb0 += ( (pb0 >> DAT_BITS) << PAD_BITS );
		pc0 += ( (pc0 >> DAT_BITS) << PAD_BITS );
		pd0 += ( (pd0 >> DAT_BITS) << PAD_BITS );
		pe0 += ( (pe0 >> DAT_BITS) << PAD_BITS );
		pf0 += ( (pf0 >> DAT_BITS) << PAD_BITS );
		p100 += ( (p100 >> DAT_BITS) << PAD_BITS );
		p110 += ( (p110 >> DAT_BITS) << PAD_BITS );
		p120 += ( (p120 >> DAT_BITS) << PAD_BITS );
		p130 += ( (p130 >> DAT_BITS) << PAD_BITS );
		p140 += ( (p140 >> DAT_BITS) << PAD_BITS );
		p150 += ( (p150 >> DAT_BITS) << PAD_BITS );
		p160 += ( (p160 >> DAT_BITS) << PAD_BITS );
		p170 += ( (p170 >> DAT_BITS) << PAD_BITS );
		p180 += ( (p180 >> DAT_BITS) << PAD_BITS );
		p190 += ( (p190 >> DAT_BITS) << PAD_BITS );
		p1a0 += ( (p1a0 >> DAT_BITS) << PAD_BITS );
		p1b0 += ( (p1b0 >> DAT_BITS) << PAD_BITS );
		p1c0 += ( (p1c0 >> DAT_BITS) << PAD_BITS );
		p1d0 += ( (p1d0 >> DAT_BITS) << PAD_BITS );
		p1e0 += ( (p1e0 >> DAT_BITS) << PAD_BITS );
		p1f0 += ( (p1f0 >> DAT_BITS) << PAD_BITS );
		p200 += ( (p200 >> DAT_BITS) << PAD_BITS );
		p210 += ( (p210 >> DAT_BITS) << PAD_BITS );
		p220 += ( (p220 >> DAT_BITS) << PAD_BITS );
		p230 += ( (p230 >> DAT_BITS) << PAD_BITS );
		p240 += ( (p240 >> DAT_BITS) << PAD_BITS );
		p250 += ( (p250 >> DAT_BITS) << PAD_BITS );
		p260 += ( (p260 >> DAT_BITS) << PAD_BITS );
		p270 += ( (p270 >> DAT_BITS) << PAD_BITS );
		p280 += ( (p280 >> DAT_BITS) << PAD_BITS );
		p290 += ( (p290 >> DAT_BITS) << PAD_BITS );
		p2a0 += ( (p2a0 >> DAT_BITS) << PAD_BITS );
		p2b0 += ( (p2b0 >> DAT_BITS) << PAD_BITS );
		p2c0 += ( (p2c0 >> DAT_BITS) << PAD_BITS );
		p2d0 += ( (p2d0 >> DAT_BITS) << PAD_BITS );
		p2e0 += ( (p2e0 >> DAT_BITS) << PAD_BITS );
		p2f0 += ( (p2f0 >> DAT_BITS) << PAD_BITS );

	// Set array offsets for radix-256 inputs [these are same for all 3 such DFTs].
		// low parts in first 16 slots, high parts in next 16:
		i_offsets_lo[0x0] = 0x00<<1;	i_offsets_hi[0x0] = 0x00<<1;
		i_offsets_lo[0x1] = 0x01<<1;	i_offsets_hi[0x1] = 0x10<<1;
		i_offsets_lo[0x2] = 0x02<<1;	i_offsets_hi[0x2] = 0x20<<1;
		i_offsets_lo[0x3] = 0x03<<1;	i_offsets_hi[0x3] = 0x30<<1;
		i_offsets_lo[0x4] = 0x04<<1;	i_offsets_hi[0x4] = 0x40<<1;
		i_offsets_lo[0x5] = 0x05<<1;	i_offsets_hi[0x5] = 0x50<<1;
		i_offsets_lo[0x6] = 0x06<<1;	i_offsets_hi[0x6] = 0x60<<1;
		i_offsets_lo[0x7] = 0x07<<1;	i_offsets_hi[0x7] = 0x70<<1;
		i_offsets_lo[0x8] = 0x08<<1;	i_offsets_hi[0x8] = 0x80<<1;
		i_offsets_lo[0x9] = 0x09<<1;	i_offsets_hi[0x9] = 0x90<<1;
		i_offsets_lo[0xa] = 0x0a<<1;	i_offsets_hi[0xa] = 0xa0<<1;
		i_offsets_lo[0xb] = 0x0b<<1;	i_offsets_hi[0xb] = 0xb0<<1;
		i_offsets_lo[0xc] = 0x0c<<1;	i_offsets_hi[0xc] = 0xc0<<1;
		i_offsets_lo[0xd] = 0x0d<<1;	i_offsets_hi[0xd] = 0xd0<<1;
		i_offsets_lo[0xe] = 0x0e<<1;	i_offsets_hi[0xe] = 0xe0<<1;
		i_offsets_lo[0xf] = 0x0f<<1;	i_offsets_hi[0xf] = 0xf0<<1;

	// For the radix-256 outputs we need 4 distinct low-part vectors, which go into 4 16-element subsets of a [64]-vec.
	// These get used in the following sequence (encoded as a 32-bit int, 16 subfields of 2 bits each) for each DFT:
	// DFT #1: [0],[1],[2],[3],[1],[2],[3],[1],[2],[3],[1],[2],[3],[1],[2],[3] = 0xe79e79e4
	// DFT #2: [1],[2],[3],[1],[2],[3],[1],[2],[3],[1],[2],[3],[1],[2],[3],[1] = 0x79e79e79
	// DFT #3: [2],[3],[1],[2],[3],[1],[2],[3],[1],[2],[3],[1],[2],[3],[1],[2] = 0x9e79e79e
	// I.e. we start with vec0, use that just once, then cycle through 1,2,3 the rest if the way.
		o_offsets_lo[0x00] =  0;	o_offsets_lo[0x10] = p5;	o_offsets_lo[0x20] = pa;	o_offsets_lo[0x30] = pf;
		o_offsets_lo[0x01] = p1;	o_offsets_lo[0x11] = p4;	o_offsets_lo[0x21] = pb;	o_offsets_lo[0x31] = pe;
		o_offsets_lo[0x02] = p2;	o_offsets_lo[0x12] = p7;	o_offsets_lo[0x22] = p9;	o_offsets_lo[0x32] = pc;
		o_offsets_lo[0x03] = p3;	o_offsets_lo[0x13] = p6;	o_offsets_lo[0x23] = p8;	o_offsets_lo[0x33] = pd;
		o_offsets_lo[0x04] = p5;	o_offsets_lo[0x14] = p2;	o_offsets_lo[0x24] = pf;	o_offsets_lo[0x34] = p9;
		o_offsets_lo[0x05] = p4;	o_offsets_lo[0x15] = p3;	o_offsets_lo[0x25] = pe;	o_offsets_lo[0x35] = p8;
		o_offsets_lo[0x06] = p7;	o_offsets_lo[0x16] = p1;	o_offsets_lo[0x26] = pc;	o_offsets_lo[0x36] = pb;
		o_offsets_lo[0x07] = p6;	o_offsets_lo[0x17] =  0;	o_offsets_lo[0x27] = pd;	o_offsets_lo[0x37] = pa;
		o_offsets_lo[0x08] = pa;	o_offsets_lo[0x18] = pf;	o_offsets_lo[0x28] = p5;	o_offsets_lo[0x38] = p2;
		o_offsets_lo[0x09] = pb;	o_offsets_lo[0x19] = pe;	o_offsets_lo[0x29] = p4;	o_offsets_lo[0x39] = p3;
		o_offsets_lo[0x0a] = p9;	o_offsets_lo[0x1a] = pc;	o_offsets_lo[0x2a] = p7;	o_offsets_lo[0x3a] = p1;
		o_offsets_lo[0x0b] = p8;	o_offsets_lo[0x1b] = pd;	o_offsets_lo[0x2b] = p6;	o_offsets_lo[0x3b] =  0;
		o_offsets_lo[0x0c] = pf;	o_offsets_lo[0x1c] = p9;	o_offsets_lo[0x2c] = p2;	o_offsets_lo[0x3c] = p7;
		o_offsets_lo[0x0d] = pe;	o_offsets_lo[0x1d] = p8;	o_offsets_lo[0x2d] = p3;	o_offsets_lo[0x3d] = p6;
		o_offsets_lo[0x0e] = pc;	o_offsets_lo[0x1e] = pb;	o_offsets_lo[0x2e] = p1;	o_offsets_lo[0x3e] = p4;
		o_offsets_lo[0x0f] = pd;	o_offsets_lo[0x1f] = pa;	o_offsets_lo[0x2f] =  0;	o_offsets_lo[0x3f] = p5;

	// ...and one distinct high-part vector for each radix-256 DFT:
		o_offsets_hi1[0x0] =   0;	o_offsets_hi2[0x0] = p50;	o_offsets_hi3[0x0] = pa0;
		o_offsets_hi1[0x1] = p10;	o_offsets_hi2[0x1] = p40;	o_offsets_hi3[0x1] = pb0;
		o_offsets_hi1[0x2] = p20;	o_offsets_hi2[0x2] = p70;	o_offsets_hi3[0x2] = p90;
		o_offsets_hi1[0x3] = p30;	o_offsets_hi2[0x3] = p60;	o_offsets_hi3[0x3] = p80;
		o_offsets_hi1[0x4] = p50;	o_offsets_hi2[0x4] = p20;	o_offsets_hi3[0x4] = pf0;
		o_offsets_hi1[0x5] = p40;	o_offsets_hi2[0x5] = p30;	o_offsets_hi3[0x5] = pe0;
		o_offsets_hi1[0x6] = p70;	o_offsets_hi2[0x6] = p10;	o_offsets_hi3[0x6] = pc0;
		o_offsets_hi1[0x7] = p60;	o_offsets_hi2[0x7] =   0;	o_offsets_hi3[0x7] = pd0;
		o_offsets_hi1[0x8] = pa0;	o_offsets_hi2[0x8] = pf0;	o_offsets_hi3[0x8] = p50;
		o_offsets_hi1[0x9] = pb0;	o_offsets_hi2[0x9] = pe0;	o_offsets_hi3[0x9] = p40;
		o_offsets_hi1[0xa] = p90;	o_offsets_hi2[0xa] = pc0;	o_offsets_hi3[0xa] = p70;
		o_offsets_hi1[0xb] = p80;	o_offsets_hi2[0xb] = pd0;	o_offsets_hi3[0xb] = p60;
		o_offsets_hi1[0xc] = pf0;	o_offsets_hi2[0xc] = p90;	o_offsets_hi3[0xc] = p20;
		o_offsets_hi1[0xd] = pe0;	o_offsets_hi2[0xd] = p80;	o_offsets_hi3[0xd] = p30;
		o_offsets_hi1[0xe] = pc0;	o_offsets_hi2[0xe] = pb0;	o_offsets_hi3[0xe] = p10;
		o_offsets_hi1[0xf] = pd0;	o_offsets_hi2[0xf] = pa0;	o_offsets_hi3[0xf] =   0;

	// Index-high-bits triplets needed for compact-obj-code scheme:
		k = 0;
		dif_triplets[k] = p2f0; dif_triplets[k+1] = p1f0; dif_triplets[k+2] =  pf0; k += 3;
		dif_triplets[k] = p2e0; dif_triplets[k+1] = p1e0; dif_triplets[k+2] =  pe0; k += 3;
		dif_triplets[k] = p2d0; dif_triplets[k+1] = p1d0; dif_triplets[k+2] =  pd0; k += 3;
		dif_triplets[k] = p2c0; dif_triplets[k+1] = p1c0; dif_triplets[k+2] =  pc0; k += 3;
		dif_triplets[k] = p2b0; dif_triplets[k+1] = p1b0; dif_triplets[k+2] =  pb0; k += 3;
		dif_triplets[k] = p2a0; dif_triplets[k+1] = p1a0; dif_triplets[k+2] =  pa0; k += 3;
		dif_triplets[k] = p290; dif_triplets[k+1] = p190; dif_triplets[k+2] =  p90; k += 3;
		dif_triplets[k] = p280; dif_triplets[k+1] = p180; dif_triplets[k+2] =  p80; k += 3;
		dif_triplets[k] = p270; dif_triplets[k+1] = p170; dif_triplets[k+2] =  p70; k += 3;
		dif_triplets[k] = p260; dif_triplets[k+1] = p160; dif_triplets[k+2] =  p60; k += 3;
		dif_triplets[k] = p250; dif_triplets[k+1] = p150; dif_triplets[k+2] =  p50; k += 3;
		dif_triplets[k] = p240; dif_triplets[k+1] = p140; dif_triplets[k+2] =  p40; k += 3;
		dif_triplets[k] = p230; dif_triplets[k+1] = p130; dif_triplets[k+2] =  p30; k += 3;
		dif_triplets[k] = p220; dif_triplets[k+1] = p120; dif_triplets[k+2] =  p20; k += 3;
		dif_triplets[k] = p210; dif_triplets[k+1] = p110; dif_triplets[k+2] =  p10; k += 3;
		dif_triplets[k] = p200; dif_triplets[k+1] = p100; dif_triplets[k+2] =    0; k += 3;
		dif_triplets[k] = p1f0; dif_triplets[k+1] =  pf0; dif_triplets[k+2] = p2f0; k += 3;
		dif_triplets[k] = p1e0; dif_triplets[k+1] =  pe0; dif_triplets[k+2] = p2e0; k += 3;
		dif_triplets[k] = p1d0; dif_triplets[k+1] =  pd0; dif_triplets[k+2] = p2d0; k += 3;
		dif_triplets[k] = p1c0; dif_triplets[k+1] =  pc0; dif_triplets[k+2] = p2c0; k += 3;
		dif_triplets[k] = p1b0; dif_triplets[k+1] =  pb0; dif_triplets[k+2] = p2b0; k += 3;
		dif_triplets[k] = p1a0; dif_triplets[k+1] =  pa0; dif_triplets[k+2] = p2a0; k += 3;
		dif_triplets[k] = p190; dif_triplets[k+1] =  p90; dif_triplets[k+2] = p290; k += 3;
		dif_triplets[k] = p180; dif_triplets[k+1] =  p80; dif_triplets[k+2] = p280; k += 3;
		dif_triplets[k] = p170; dif_triplets[k+1] =  p70; dif_triplets[k+2] = p270; k += 3;
		dif_triplets[k] = p160; dif_triplets[k+1] =  p60; dif_triplets[k+2] = p260; k += 3;
		dif_triplets[k] = p150; dif_triplets[k+1] =  p50; dif_triplets[k+2] = p250; k += 3;
		dif_triplets[k] = p140; dif_triplets[k+1] =  p40; dif_triplets[k+2] = p240; k += 3;
		dif_triplets[k] = p130; dif_triplets[k+1] =  p30; dif_triplets[k+2] = p230; k += 3;
		dif_triplets[k] = p120; dif_triplets[k+1] =  p20; dif_triplets[k+2] = p220; k += 3;
		dif_triplets[k] = p110; dif_triplets[k+1] =  p10; dif_triplets[k+2] = p210; k += 3;
		dif_triplets[k] = p100; dif_triplets[k+1] =    0; dif_triplets[k+2] = p200; k += 3;
		dif_triplets[k] =  pf0; dif_triplets[k+1] = p2f0; dif_triplets[k+2] = p1f0; k += 3;
		dif_triplets[k] =  pe0; dif_triplets[k+1] = p2e0; dif_triplets[k+2] = p1e0; k += 3;
		dif_triplets[k] =  pd0; dif_triplets[k+1] = p2d0; dif_triplets[k+2] = p1d0; k += 3;
		dif_triplets[k] =  pc0; dif_triplets[k+1] = p2c0; dif_triplets[k+2] = p1c0; k += 3;
		dif_triplets[k] =  pb0; dif_triplets[k+1] = p2b0; dif_triplets[k+2] = p1b0; k += 3;
		dif_triplets[k] =  pa0; dif_triplets[k+1] = p2a0; dif_triplets[k+2] = p1a0; k += 3;
		dif_triplets[k] =  p90; dif_triplets[k+1] = p290; dif_triplets[k+2] = p190; k += 3;
		dif_triplets[k] =  p80; dif_triplets[k+1] = p280; dif_triplets[k+2] = p180; k += 3;
		dif_triplets[k] =  p70; dif_triplets[k+1] = p270; dif_triplets[k+2] = p170; k += 3;
		dif_triplets[k] =  p60; dif_triplets[k+1] = p260; dif_triplets[k+2] = p160; k += 3;
		dif_triplets[k] =  p50; dif_triplets[k+1] = p250; dif_triplets[k+2] = p150; k += 3;
		dif_triplets[k] =  p40; dif_triplets[k+1] = p240; dif_triplets[k+2] = p140; k += 3;
		dif_triplets[k] =  p30; dif_triplets[k+1] = p230; dif_triplets[k+2] = p130; k += 3;
		dif_triplets[k] =  p20; dif_triplets[k+1] = p220; dif_triplets[k+2] = p120; k += 3;
		dif_triplets[k] =  p10; dif_triplets[k+1] = p210; dif_triplets[k+2] = p110; k += 3;
		dif_triplets[k] =    0; dif_triplets[k+1] = p200; dif_triplets[k+2] = p100;

	}

/*...The radix-768 pass is here.	*/

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
	Twiddleless version arranges 256 sets of radix-3 DFT inputs as follows: 0 in upper left corner,
	decrement 256 horizontally and 3 vertically, all (mod 768). Break into 4 cols for readability,
	and insert blank rows to reveal 16-macro-calls subgroupings exploited by compact-obj-code scheme below:
		000,200,100 + p0	<*** Wrap this around to end...

		2f0,1f0,0f0 + pd		230,130,030 + pd		170,070,270 + pd		0b0,2b0,1b0 + pd
		2f0,1f0,0f0 + pa		230,130,030 + pa		170,070,270 + pa		0b0,2b0,1b0 + pa
		2f0,1f0,0f0 + p7		230,130,030 + p7		170,070,270 + p7		0b0,2b0,1b0 + p7
		2f0,1f0,0f0 + p4		230,130,030 + p4		170,070,270 + p4		0b0,2b0,1b0 + p4
		2f0,1f0,0f0 + p1		230,130,030 + p1		170,070,270 + p1		0b0,2b0,1b0 + p1
		2e0,1e0,0e0 + pe		220,120,020 + pe		160,060,260 + pe		0a0,2a0,1a0 + pe
		2e0,1e0,0e0 + pb		220,120,020 + pb		160,060,260 + pb		0a0,2a0,1a0 + pb
		2e0,1e0,0e0 + p8		220,120,020 + p8		160,060,260 + p8		0a0,2a0,1a0 + p8
		2e0,1e0,0e0 + p5		220,120,020 + p5		160,060,260 + p5		0a0,2a0,1a0 + p5
		2e0,1e0,0e0 + p2		220,120,020 + p2		160,060,260 + p2		0a0,2a0,1a0 + p2
		2d0,1d0,0d0 + pf		210,110,010 + pf		150,050,250 + pf		090,290,190 + pf
		2d0,1d0,0d0 + pc		210,110,010 + pc		150,050,250 + pc		090,290,190 + pc
		2d0,1d0,0d0 + p9		210,110,010 + p9		150,050,250 + p9		090,290,190 + p9
		2d0,1d0,0d0 + p6		210,110,010 + p6		150,050,250 + p6		090,290,190 + p6
		2d0,1d0,0d0 + p3		210,110,010 + p3		150,050,250 + p3		090,290,190 + p3
		2d0,1d0,0d0 + p0		210,110,010 + p0		150,050,250 + p0		090,290,190 + p0

		2c0,1c0,0c0 + pd		200,100,000 + pd		140,040,240 + pd		080,280,180 + pd
		2c0,1c0,0c0 + pa		200,100,000 + pa		140,040,240 + pa		080,280,180 + pa
		2c0,1c0,0c0 + p7		200,100,000 + p7		140,040,240 + p7		080,280,180 + p7
		2c0,1c0,0c0 + p4		200,100,000 + p4		140,040,240 + p4		080,280,180 + p4
		2c0,1c0,0c0 + p1		200,100,000 + p1		140,040,240 + p1		080,280,180 + p1
		2b0,1b0,0b0 + pe		1f0,0f0,2f0 + pe		130,030,230 + pe		070,270,170 + pe
		2b0,1b0,0b0 + pb		1f0,0f0,2f0 + pb		130,030,230 + pb		070,270,170 + pb
		2b0,1b0,0b0 + p8		1f0,0f0,2f0 + p8		130,030,230 + p8		070,270,170 + p8
		2b0,1b0,0b0 + p5		1f0,0f0,2f0 + p5		130,030,230 + p5		070,270,170 + p5
		2b0,1b0,0b0 + p2		1f0,0f0,2f0 + p2		130,030,230 + p2		070,270,170 + p2
		2a0,1a0,0a0 + pf		1e0,0e0,2e0 + pf		120,020,220 + pf		060,260,160 + pf
		2a0,1a0,0a0 + pc		1e0,0e0,2e0 + pc		120,020,220 + pc		060,260,160 + pc
		2a0,1a0,0a0 + p9		1e0,0e0,2e0 + p9		120,020,220 + p9		060,260,160 + p9
		2a0,1a0,0a0 + p6		1e0,0e0,2e0 + p6		120,020,220 + p6		060,260,160 + p6
		2a0,1a0,0a0 + p3		1e0,0e0,2e0 + p3		120,020,220 + p3		060,260,160 + p3
		2a0,1a0,0a0 + p0		1e0,0e0,2e0 + p0		120,020,220 + p0		060,260,160 + p0

		290,190,090 + pd		1d0,0d0,2d0 + pd		110,010,210 + pd		050,250,150 + pd
		290,190,090 + pa		1d0,0d0,2d0 + pa		110,010,210 + pa		050,250,150 + pa
		290,190,090 + p7		1d0,0d0,2d0 + p7		110,010,210 + p7		050,250,150 + p7
		290,190,090 + p4		1d0,0d0,2d0 + p4		110,010,210 + p4		050,250,150 + p4
		290,190,090 + p1		1d0,0d0,2d0 + p1		110,010,210 + p1		050,250,150 + p1
		280,180,080 + pe		1c0,0c0,2c0 + pe		100,000,200 + pe		040,240,140 + pe
		280,180,080 + pb		1c0,0c0,2c0 + pb		100,000,200 + pb		040,240,140 + pb
		280,180,080 + p8		1c0,0c0,2c0 + p8		100,000,200 + p8		040,240,140 + p8
		280,180,080 + p5		1c0,0c0,2c0 + p5		100,000,200 + p5		040,240,140 + p5
		280,180,080 + p2		1c0,0c0,2c0 + p2		100,000,200 + p2		040,240,140 + p2
		270,170,070 + pf		1b0,0b0,2b0 + pf		0f0,2f0,1f0 + pf		030,230,130 + pf
		270,170,070 + pc		1b0,0b0,2b0 + pc		0f0,2f0,1f0 + pc		030,230,130 + pc
		270,170,070 + p9		1b0,0b0,2b0 + p9		0f0,2f0,1f0 + p9		030,230,130 + p9
		270,170,070 + p6		1b0,0b0,2b0 + p6		0f0,2f0,1f0 + p6		030,230,130 + p6
		270,170,070 + p3		1b0,0b0,2b0 + p3		0f0,2f0,1f0 + p3		030,230,130 + p3
		270,170,070 + p0		1b0,0b0,2b0 + p0		0f0,2f0,1f0 + p0		030,230,130 + p0

		260,160,060 + pd		1a0,0a0,2a0 + pd		0e0,2e0,1e0 + pd		020,220,120 + pd
		260,160,060 + pa		1a0,0a0,2a0 + pa		0e0,2e0,1e0 + pa		020,220,120 + pa
		260,160,060 + p7		1a0,0a0,2a0 + p7		0e0,2e0,1e0 + p7		020,220,120 + p7
		260,160,060 + p4		1a0,0a0,2a0 + p4		0e0,2e0,1e0 + p4		020,220,120 + p4
		260,160,060 + p1		1a0,0a0,2a0 + p1		0e0,2e0,1e0 + p1		020,220,120 + p1
		250,150,050 + pe		190,090,290 + pe		0d0,2d0,1d0 + pe		010,210,110 + pe
		250,150,050 + pb		190,090,290 + pb		0d0,2d0,1d0 + pb		010,210,110 + pb
		250,150,050 + p8		190,090,290 + p8		0d0,2d0,1d0 + p8		010,210,110 + p8
		250,150,050 + p5		190,090,290 + p5		0d0,2d0,1d0 + p5		010,210,110 + p5
		250,150,050 + p2		190,090,290 + p2		0d0,2d0,1d0 + p2		010,210,110 + p2
		240,140,040 + pf		180,080,280 + pf		0c0,2c0,1c0 + pf		000,200,100 + pf
		240,140,040 + pc		180,080,280 + pc		0c0,2c0,1c0 + pc		000,200,100 + pc
		240,140,040 + p9		180,080,280 + p9		0c0,2c0,1c0 + p9		000,200,100 + p9
		240,140,040 + p6		180,080,280 + p6		0c0,2c0,1c0 + p6		000,200,100 + p6
		240,140,040 + p3		180,080,280 + p3		0c0,2c0,1c0 + p3		000,200,100 + p3
		240,140,040 + p0		180,080,280 + p0		0c0,2c0,1c0 + p0		000,200,100 + p0	<*** 0-term wrapped ***
		[cont. in col2]			[cont. in col3]			[cont. in col4]

	To handle the wraparound of the 0-term we just need a bit of mod-256 indexing magic.
	*/
	/*...gather the needed data (768 64-bit complex) and do 256 radix-3 transforms - We want unit-strides in the radix768-DFT macro, so use large output strides here: */
		// Loop-based compact-obj-code impl exploits above index pattern to group interior [sandwiched between single leading
		// and 15 trailing, which also get fused into a final group] macro calls into sets of 16 with neatly cutoff index groupings:
		l = 1; l1 = l+256; l2 = l+512;	// Skip 0-term, which gets saved for wraparound
		for(k = 0; k < 144; k += 3) {
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
			jt = j1 + p3; jp = j2 + p3;	RADIX_03_DFT(s,c3m1, a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2], t00,t01,t02,t03,t04,t05, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im); ++l; l &= 0xff; l1 = l+256; l2 = l+512;	// <*** needed for final-loop-pass wraparound of 0-term
										RADIX_03_DFT(s,c3m1, a[j1+k0],a[j2+k0],a[j1+k1],a[j2+k1],a[j1+k2],a[j2+k2], t00,t01,t02,t03,t04,t05, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im); ++l; ++l1; ++l2;
		}

	/*...and now do 3 radix-256 transforms, first assuming in-order output index offsets, then using the
	resulting data-mismatch tables produced by test_fft_radix() to derive the needed output permutation:
	*/
		// NOTE: Due to casting-to-double of inputs temp-array pointer Radix-32 macro doesn't play nice with r-array
		// offsets of form "r+32", so make address-taking explicit and wrapping in () prior to macro code execution:
		//
		// Since i_data are in a local-complex scratch array, need to override any non-unity SIMD-re/im data stride with 1:
		//                               vvv
		RADIX_256_DIF((double *)(t+0x000),1,i_offsets_lo,i_offsets_hi, (a+j1     ),RE_IM_STRIDE,o_offsets_lo,o_idx1,o_offsets_hi1);	/* Inputs in t[ 00- ff] */
		RADIX_256_DIF((double *)(t+0x100),1,i_offsets_lo,i_offsets_hi, (a+j1+p200),RE_IM_STRIDE,o_offsets_lo,o_idx2,o_offsets_hi2);	/* Inputs in t[100-1ff] */
		RADIX_256_DIF((double *)(t+0x200),1,i_offsets_lo,i_offsets_hi, (a+j1+p100),RE_IM_STRIDE,o_offsets_lo,o_idx3,o_offsets_hi3);	/* Inputs in t[200-2ff] */
	}
}

/***************/

void radix768_dit_pass1(double a[], int n)
{
/*
!...Acronym: DIT = Decimation In Time
!
!...Subroutine to perform an initial radix-768 complex DIT FFT pass on the data in the length-N real vector A.
!
!   See the documentation in radix16_dif_pass for further details on storage and indexing.
*/
	int j,j1,j2,jp,jt;
	int k,l,l1,l2,k0,k1,k2;
	static int dit_triplets[48];	// Only need 1/3 as many here as for DIF
	// In order to preserve length-2 numeric-index-offset property here, use hex for digit:
	static int NDIVR, first_entry=TRUE,
		p1,p2,p3,p4,p5,p6,p7,p8,p9,pa,pb,pc,pd,pe,pf,
		     p10,p20,p30,p40,p50,p60,p70,p80,p90,pa0,pb0,pc0,pd0,pe0,pf0,
		p100,p110,p120,p130,p140,p150,p160,p170,p180,p190,p1a0,p1b0,p1c0,p1d0,p1e0,p1f0,
		p200,p210,p220,p230,p240,p250,p260,p270,p280,p290,p2a0,p2b0,p2c0,p2d0,p2e0,p2f0;
	const double c3m1 = -1.50000000000000000000, s = 0.86602540378443864675;	// cos(twopi/3)-1, sin(twopi/3)
	static int o_offsets_lo[16], o_offsets_hi[16];

	static int i_offsets_lo[96];	// 6 subsets of 16
	// Bitfields encoding the sequence of the i_offsets_lo subset0-3 vectors to use for each radix-256 DFT's outputs:
	const uint32 i_idx1 = 0x55555554,i_idx2 = 0x44404040,i_idx3 = 0x54445444;
	static int i_offsets_hi1[16],i_offsets_hi2[16],i_offsets_hi3[16];
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
		p1 = NDIVR;		p10 = NDIVR<<4;		p100 = NDIVR<<8;		p200 = NDIVR<<9;
		p2 = p1 + p1;	p20 = p10 + p10;	p110 = p100 + p10;		p210 = p200 + p10;
		p3 = p2 + p1;	p30 = p20 + p10;	p120 = p110 + p10;		p220 = p210 + p10;
		p4 = p3 + p1;	p40 = p30 + p10;	p130 = p120 + p10;		p230 = p220 + p10;
		p5 = p4 + p1;	p50 = p40 + p10;	p140 = p130 + p10;		p240 = p230 + p10;
		p6 = p5 + p1;	p60 = p50 + p10;	p150 = p140 + p10;		p250 = p240 + p10;
		p7 = p6 + p1;	p70 = p60 + p10;	p160 = p150 + p10;		p260 = p250 + p10;
		p8 = p7 + p1;	p80 = p70 + p10;	p170 = p160 + p10;		p270 = p260 + p10;
		p9 = p8 + p1;	p90 = p80 + p10;	p180 = p170 + p10;		p280 = p270 + p10;
		pa = p9 + p1;	pa0 = p90 + p10;	p190 = p180 + p10;		p290 = p280 + p10;
		pb = pa + p1;	pb0 = pa0 + p10;	p1a0 = p190 + p10;		p2a0 = p290 + p10;
		pc = pb + p1;	pc0 = pb0 + p10;	p1b0 = p1a0 + p10;		p2b0 = p2a0 + p10;
		pd = pc + p1;	pd0 = pc0 + p10;	p1c0 = p1b0 + p10;		p2c0 = p2b0 + p10;
		pe = pd + p1;	pe0 = pd0 + p10;	p1d0 = p1c0 + p10;		p2d0 = p2c0 + p10;
		pf = pe + p1;	pf0 = pe0 + p10;	p1e0 = p1d0 + p10;		p2e0 = p2d0 + p10;
											p1f0 = p1e0 + p10;		p2f0 = p2e0 + p10;
		p1 += ( (p1 >> DAT_BITS) << PAD_BITS );
		p2 += ( (p2 >> DAT_BITS) << PAD_BITS );
		p3 += ( (p3 >> DAT_BITS) << PAD_BITS );
		p4 += ( (p4 >> DAT_BITS) << PAD_BITS );
		p5 += ( (p5 >> DAT_BITS) << PAD_BITS );
		p6 += ( (p6 >> DAT_BITS) << PAD_BITS );
		p7 += ( (p7 >> DAT_BITS) << PAD_BITS );
		p8 += ( (p8 >> DAT_BITS) << PAD_BITS );
		p9 += ( (p9 >> DAT_BITS) << PAD_BITS );
		pa += ( (pa >> DAT_BITS) << PAD_BITS );
		pb += ( (pb >> DAT_BITS) << PAD_BITS );
		pc += ( (pc >> DAT_BITS) << PAD_BITS );
		pd += ( (pd >> DAT_BITS) << PAD_BITS );
		pe += ( (pe >> DAT_BITS) << PAD_BITS );
		pf += ( (pf >> DAT_BITS) << PAD_BITS );
		p10 += ( (p10 >> DAT_BITS) << PAD_BITS );
		p20 += ( (p20 >> DAT_BITS) << PAD_BITS );
		p30 += ( (p30 >> DAT_BITS) << PAD_BITS );
		p40 += ( (p40 >> DAT_BITS) << PAD_BITS );
		p50 += ( (p50 >> DAT_BITS) << PAD_BITS );
		p60 += ( (p60 >> DAT_BITS) << PAD_BITS );
		p70 += ( (p70 >> DAT_BITS) << PAD_BITS );
		p80 += ( (p80 >> DAT_BITS) << PAD_BITS );
		p90 += ( (p90 >> DAT_BITS) << PAD_BITS );
		pa0 += ( (pa0 >> DAT_BITS) << PAD_BITS );
		pb0 += ( (pb0 >> DAT_BITS) << PAD_BITS );
		pc0 += ( (pc0 >> DAT_BITS) << PAD_BITS );
		pd0 += ( (pd0 >> DAT_BITS) << PAD_BITS );
		pe0 += ( (pe0 >> DAT_BITS) << PAD_BITS );
		pf0 += ( (pf0 >> DAT_BITS) << PAD_BITS );
		p100 += ( (p100 >> DAT_BITS) << PAD_BITS );
		p110 += ( (p110 >> DAT_BITS) << PAD_BITS );
		p120 += ( (p120 >> DAT_BITS) << PAD_BITS );
		p130 += ( (p130 >> DAT_BITS) << PAD_BITS );
		p140 += ( (p140 >> DAT_BITS) << PAD_BITS );
		p150 += ( (p150 >> DAT_BITS) << PAD_BITS );
		p160 += ( (p160 >> DAT_BITS) << PAD_BITS );
		p170 += ( (p170 >> DAT_BITS) << PAD_BITS );
		p180 += ( (p180 >> DAT_BITS) << PAD_BITS );
		p190 += ( (p190 >> DAT_BITS) << PAD_BITS );
		p1a0 += ( (p1a0 >> DAT_BITS) << PAD_BITS );
		p1b0 += ( (p1b0 >> DAT_BITS) << PAD_BITS );
		p1c0 += ( (p1c0 >> DAT_BITS) << PAD_BITS );
		p1d0 += ( (p1d0 >> DAT_BITS) << PAD_BITS );
		p1e0 += ( (p1e0 >> DAT_BITS) << PAD_BITS );
		p1f0 += ( (p1f0 >> DAT_BITS) << PAD_BITS );
		p200 += ( (p200 >> DAT_BITS) << PAD_BITS );
		p210 += ( (p210 >> DAT_BITS) << PAD_BITS );
		p220 += ( (p220 >> DAT_BITS) << PAD_BITS );
		p230 += ( (p230 >> DAT_BITS) << PAD_BITS );
		p240 += ( (p240 >> DAT_BITS) << PAD_BITS );
		p250 += ( (p250 >> DAT_BITS) << PAD_BITS );
		p260 += ( (p260 >> DAT_BITS) << PAD_BITS );
		p270 += ( (p270 >> DAT_BITS) << PAD_BITS );
		p280 += ( (p280 >> DAT_BITS) << PAD_BITS );
		p290 += ( (p290 >> DAT_BITS) << PAD_BITS );
		p2a0 += ( (p2a0 >> DAT_BITS) << PAD_BITS );
		p2b0 += ( (p2b0 >> DAT_BITS) << PAD_BITS );
		p2c0 += ( (p2c0 >> DAT_BITS) << PAD_BITS );
		p2d0 += ( (p2d0 >> DAT_BITS) << PAD_BITS );
		p2e0 += ( (p2e0 >> DAT_BITS) << PAD_BITS );
		p2f0 += ( (p2f0 >> DAT_BITS) << PAD_BITS );

	// Set array offsets for radix-256 outputs [these are same for all 3 such DFTs].
		// low parts in first 16 slots, high parts in next 16:
		o_offsets_lo[0x0] = 0x00<<1;	o_offsets_hi[0x0] = 0x00<<1;
		o_offsets_lo[0x1] = 0x01<<1;	o_offsets_hi[0x1] = 0x10<<1;
		o_offsets_lo[0x2] = 0x02<<1;	o_offsets_hi[0x2] = 0x20<<1;
		o_offsets_lo[0x3] = 0x03<<1;	o_offsets_hi[0x3] = 0x30<<1;
		o_offsets_lo[0x4] = 0x04<<1;	o_offsets_hi[0x4] = 0x40<<1;
		o_offsets_lo[0x5] = 0x05<<1;	o_offsets_hi[0x5] = 0x50<<1;
		o_offsets_lo[0x6] = 0x06<<1;	o_offsets_hi[0x6] = 0x60<<1;
		o_offsets_lo[0x7] = 0x07<<1;	o_offsets_hi[0x7] = 0x70<<1;
		o_offsets_lo[0x8] = 0x08<<1;	o_offsets_hi[0x8] = 0x80<<1;
		o_offsets_lo[0x9] = 0x09<<1;	o_offsets_hi[0x9] = 0x90<<1;
		o_offsets_lo[0xa] = 0x0a<<1;	o_offsets_hi[0xa] = 0xa0<<1;
		o_offsets_lo[0xb] = 0x0b<<1;	o_offsets_hi[0xb] = 0xb0<<1;
		o_offsets_lo[0xc] = 0x0c<<1;	o_offsets_hi[0xc] = 0xc0<<1;
		o_offsets_lo[0xd] = 0x0d<<1;	o_offsets_hi[0xd] = 0xd0<<1;
		o_offsets_lo[0xe] = 0x0e<<1;	o_offsets_hi[0xe] = 0xe0<<1;
		o_offsets_lo[0xf] = 0x0f<<1;	o_offsets_hi[0xf] = 0xf0<<1;

	// For the radix-256 inputs we need 6 distinct low-part vectors, which go into 6 16-element subsets of a [96]-vec.
	// These get used in the following sequence (encoded as a 32-bit int, 16 subfields of 2 bits each) for each DFT.
	// NOTE that the multi-16-vector scheme we devised for radix-256 DIF formally permits just 4 distinct low-part 16-vectors
	// but that assumed that any of the 4 could appear in a given radix-256 DIF. Here we have 6 distinct 16-vectors
	// but each radix-256 DIT uses just 2 of them, thus we can deploy the same scheme, but sending the pointer to the
	// low one of the 16-element vector pair used by the current radix-256 DIF, thus each 2-bit index-offset subfield
	// of our uint32 bitfields take just values 0 and 1, odd-order bits unused:
	// DFT #1:       [0],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1] <<< [left->right  = 0x55555554
	// DFT #2: [2] + [0],[0],[0],[1],[0],[0],[0],[1],[0],[0],[0],[1],[0],[1],[0],[1] <<< in terms of   = 0x44404040
	// DFT #3: [4] + [0],[1],[0],[1],[0],[1],[1],[1],[0],[1],[0],[1],[0],[1],[1],[1] <<< significance] = 0x54445444 .

	// Set i_offsets_lo for 1st set of radix-256 DIT inputs:
		i_offsets_lo[0x00] =  0;	i_offsets_lo[0x10] = pf;
		i_offsets_lo[0x01] = p1;	i_offsets_lo[0x11] = pe;
		i_offsets_lo[0x02] = p3;	i_offsets_lo[0x12] = pd;
		i_offsets_lo[0x03] = p2;	i_offsets_lo[0x13] = pc;
		i_offsets_lo[0x04] = p7;	i_offsets_lo[0x14] = pb;
		i_offsets_lo[0x05] = p6;	i_offsets_lo[0x15] = pa;
		i_offsets_lo[0x06] = p5;	i_offsets_lo[0x16] = p9;
		i_offsets_lo[0x07] = p4;	i_offsets_lo[0x17] = p8;
		i_offsets_lo[0x08] = pf;	i_offsets_lo[0x18] = p7;
		i_offsets_lo[0x09] = pe;	i_offsets_lo[0x19] = p6;
		i_offsets_lo[0x0a] = pd;	i_offsets_lo[0x1a] = p5;
		i_offsets_lo[0x0b] = pc;	i_offsets_lo[0x1b] = p4;
		i_offsets_lo[0x0c] = pb;	i_offsets_lo[0x1c] = p3;
		i_offsets_lo[0x0d] = pa;	i_offsets_lo[0x1d] = p2;
		i_offsets_lo[0x0e] = p9;	i_offsets_lo[0x1e] = p1;
		i_offsets_lo[0x0f] = p8;	i_offsets_lo[0x1f] =  0;

	// Set i_offsets for 2nd set of radix-256 DIT inputs:
		i_offsets_lo[0x20] = p5;	i_offsets_lo[0x30] = p9;
		i_offsets_lo[0x21] = p4;	i_offsets_lo[0x31] = p8;
		i_offsets_lo[0x22] = p6;	i_offsets_lo[0x32] = pa;
		i_offsets_lo[0x23] = p7;	i_offsets_lo[0x33] = pb;
		i_offsets_lo[0x24] = p1;	i_offsets_lo[0x34] = pe;
		i_offsets_lo[0x25] =  0;	i_offsets_lo[0x35] = pf;
		i_offsets_lo[0x26] = p2;	i_offsets_lo[0x36] = pc;
		i_offsets_lo[0x27] = p3;	i_offsets_lo[0x37] = pd;
		i_offsets_lo[0x28] = p9;	i_offsets_lo[0x38] = p1;
		i_offsets_lo[0x29] = p8;	i_offsets_lo[0x39] =  0;
		i_offsets_lo[0x2a] = pa;	i_offsets_lo[0x3a] = p2;
		i_offsets_lo[0x2b] = pb;	i_offsets_lo[0x3b] = p3;
		i_offsets_lo[0x2c] = pe;	i_offsets_lo[0x3c] = p6;
		i_offsets_lo[0x2d] = pf;	i_offsets_lo[0x3d] = p7;
		i_offsets_lo[0x2e] = pc;	i_offsets_lo[0x3e] = p4;
		i_offsets_lo[0x2f] = pd;	i_offsets_lo[0x3f] = p5;

	// Set i_offsets for 3rd set of radix-256 DIT inputs:
		i_offsets_lo[0x40] = pa;	i_offsets_lo[0x50] = p2;
		i_offsets_lo[0x41] = pb;	i_offsets_lo[0x51] = p3;
		i_offsets_lo[0x42] = p8;	i_offsets_lo[0x52] =  0;
		i_offsets_lo[0x43] = p9;	i_offsets_lo[0x53] = p1;
		i_offsets_lo[0x44] = pc;	i_offsets_lo[0x54] = p4;
		i_offsets_lo[0x45] = pd;	i_offsets_lo[0x55] = p5;
		i_offsets_lo[0x46] = pf;	i_offsets_lo[0x56] = p7;
		i_offsets_lo[0x47] = pe;	i_offsets_lo[0x57] = p6;
		i_offsets_lo[0x48] = p2;	i_offsets_lo[0x58] = pc;
		i_offsets_lo[0x49] = p3;	i_offsets_lo[0x59] = pd;
		i_offsets_lo[0x4a] =  0;	i_offsets_lo[0x5a] = pf;
		i_offsets_lo[0x4b] = p1;	i_offsets_lo[0x5b] = pe;
		i_offsets_lo[0x4c] = p4;	i_offsets_lo[0x5c] = p8;
		i_offsets_lo[0x4d] = p5;	i_offsets_lo[0x5d] = p9;
		i_offsets_lo[0x4e] = p7;	i_offsets_lo[0x5e] = pb;
		i_offsets_lo[0x4f] = p6;	i_offsets_lo[0x5f] = pa;

	// ...and one distinct high-part vector for each radix-256 DIT:
		i_offsets_hi1[0x0] =   0;	i_offsets_hi2[0x0] = p50;	i_offsets_hi3[0x0] = pa0;
		i_offsets_hi1[0x1] = p10;	i_offsets_hi2[0x1] = p40;	i_offsets_hi3[0x1] = pb0;
		i_offsets_hi1[0x2] = p30;	i_offsets_hi2[0x2] = p60;	i_offsets_hi3[0x2] = p80;
		i_offsets_hi1[0x3] = p20;	i_offsets_hi2[0x3] = p70;	i_offsets_hi3[0x3] = p90;
		i_offsets_hi1[0x4] = p70;	i_offsets_hi2[0x4] = p10;	i_offsets_hi3[0x4] = pc0;
		i_offsets_hi1[0x5] = p60;	i_offsets_hi2[0x5] =   0;	i_offsets_hi3[0x5] = pd0;
		i_offsets_hi1[0x6] = p50;	i_offsets_hi2[0x6] = p20;	i_offsets_hi3[0x6] = pf0;
		i_offsets_hi1[0x7] = p40;	i_offsets_hi2[0x7] = p30;	i_offsets_hi3[0x7] = pe0;
		i_offsets_hi1[0x8] = pf0;	i_offsets_hi2[0x8] = p90;	i_offsets_hi3[0x8] = p20;
		i_offsets_hi1[0x9] = pe0;	i_offsets_hi2[0x9] = p80;	i_offsets_hi3[0x9] = p30;
		i_offsets_hi1[0xa] = pd0;	i_offsets_hi2[0xa] = pa0;	i_offsets_hi3[0xa] =   0;
		i_offsets_hi1[0xb] = pc0;	i_offsets_hi2[0xb] = pb0;	i_offsets_hi3[0xb] = p10;
		i_offsets_hi1[0xc] = pb0;	i_offsets_hi2[0xc] = pe0;	i_offsets_hi3[0xc] = p40;
		i_offsets_hi1[0xd] = pa0;	i_offsets_hi2[0xd] = pf0;	i_offsets_hi3[0xd] = p50;
		i_offsets_hi1[0xe] = p90;	i_offsets_hi2[0xe] = pc0;	i_offsets_hi3[0xe] = p70;
		i_offsets_hi1[0xf] = p80;	i_offsets_hi2[0xf] = pd0;	i_offsets_hi3[0xf] = p60;

	// Index-high-bits triplets needed for compact-obj-code scheme:
		k = 0;
		dit_triplets[k] =  pf0; dit_triplets[k+1] = p1f0; dit_triplets[k+2] = p2f0; k += 3;
		dit_triplets[k] = p1e0; dit_triplets[k+1] = p2e0; dit_triplets[k+2] =  pe0; k += 3;
		dit_triplets[k] = p2d0; dit_triplets[k+1] =  pd0; dit_triplets[k+2] = p1d0; k += 3;
		dit_triplets[k] =  pc0; dit_triplets[k+1] = p1c0; dit_triplets[k+2] = p2c0; k += 3;
		dit_triplets[k] = p1b0; dit_triplets[k+1] = p2b0; dit_triplets[k+2] =  pb0; k += 3;
		dit_triplets[k] = p2a0; dit_triplets[k+1] =  pa0; dit_triplets[k+2] = p1a0; k += 3;
		dit_triplets[k] =  p90; dit_triplets[k+1] = p190; dit_triplets[k+2] = p290; k += 3;
		dit_triplets[k] = p180; dit_triplets[k+1] = p280; dit_triplets[k+2] =  p80; k += 3;
		dit_triplets[k] = p270; dit_triplets[k+1] =  p70; dit_triplets[k+2] = p170; k += 3;
		dit_triplets[k] =  p60; dit_triplets[k+1] = p160; dit_triplets[k+2] = p260; k += 3;
		dit_triplets[k] = p150; dit_triplets[k+1] = p250; dit_triplets[k+2] =  p50; k += 3;
		dit_triplets[k] = p240; dit_triplets[k+1] =  p40; dit_triplets[k+2] = p140; k += 3;
		dit_triplets[k] =  p30; dit_triplets[k+1] = p130; dit_triplets[k+2] = p230; k += 3;
		dit_triplets[k] = p120; dit_triplets[k+1] = p220; dit_triplets[k+2] =  p20; k += 3;
		dit_triplets[k] = p210; dit_triplets[k+1] =  p10; dit_triplets[k+2] = p110; k += 3;
		dit_triplets[k] =    0; dit_triplets[k+1] = p100; dit_triplets[k+2] = p200;
	}

/*...The radix-768 pass is here.	*/

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
	store the 3 index-offset 256-tets going into the radix-256 DFTs in the i_offset array.

	Radix-3 high-part p-index offsets tabulated below - Break into 4 cols for readability,
	and insert blank rows to reveal 16-macro-calls subgroupings exploited by the compact-obj-code scheme.
	Notice that each 16-triplet set has 3 distinct triplet patterns repeating in "bookend fashion",
	ABCABCABCACBABCA, with the first and last triplet the same [A]. 15 such loop passes means 45 triplets
	which are precomputed, same number as for the DIF, even though the triplet patterns in each loop pass
	differ for the 2 flavors of DFT. But we can cut that in 1/3 because we note that the 3 triplet patterns
	A,B,c within each loop pass are simply the 3 distinct circular-perms of the A-triplet:

		000,100,200 + p0	<*** Wrap this around to end...

		0f0,1f0,2f0 + pf		1b0,2b0,0b0 + pf		270,070,170 + pf		030,130,230 + pf
		1f0,2f0,0f0 + pe		2b0,0b0,1b0 + pe		070,170,270 + pe		130,230,030 + pe
		2f0,0f0,1f0 + pd		0b0,1b0,2b0 + pd		170,270,070 + pd		230,030,130 + pd
		0f0,1f0,2f0 + pc		1b0,2b0,0b0 + pc		270,070,170 + pc		030,130,230 + pc
		1f0,2f0,0f0 + pb		2b0,0b0,1b0 + pb		070,170,270 + pb		130,230,030 + pb
		2f0,0f0,1f0 + pa		0b0,1b0,2b0 + pa		170,270,070 + pa		230,030,130 + pa
		0f0,1f0,2f0 + p9		1b0,2b0,0b0 + p9		270,070,170 + p9		030,130,230 + p9
		1f0,2f0,0f0 + p8		2b0,0b0,1b0 + p8		070,170,270 + p8		130,230,030 + p8
		2f0,0f0,1f0 + p7		0b0,1b0,2b0 + p7		170,270,070 + p7		230,030,130 + p7
		0f0,1f0,2f0 + p6		1b0,2b0,0b0 + p6		270,070,170 + p6		030,130,230 + p6
		1f0,2f0,0f0 + p5		2b0,0b0,1b0 + p5		070,170,270 + p5		130,230,030 + p5
		2f0,0f0,1f0 + p4		0b0,1b0,2b0 + p4		170,270,070 + p4		230,030,130 + p4
		0f0,1f0,2f0 + p3		1b0,2b0,0b0 + p3		270,070,170 + p3		030,130,230 + p3
		1f0,2f0,0f0 + p2		2b0,0b0,1b0 + p2		070,170,270 + p2		130,230,030 + p2
		2f0,0f0,1f0 + p1		0b0,1b0,2b0 + p1		170,270,070 + p1		230,030,130 + p1
		0f0,1f0,2f0 + p0		1b0,2b0,0b0 + p0		270,070,170 + p0		030,130,230 + p0

		1e0,2e0,0e0 + pf		2a0,0a0,1a0 + pf		060,160,260 + pf		120,220,020 + pf
		2e0,0e0,1e0 + pe		0a0,1a0,2a0 + pe		160,260,060 + pe		220,020,120 + pe
		0e0,1e0,2e0 + pd		1a0,2a0,0a0 + pd		260,060,160 + pd		020,120,220 + pd
		1e0,2e0,0e0 + pc		2a0,0a0,1a0 + pc		060,160,260 + pc		120,220,020 + pc
		2e0,0e0,1e0 + pb		0a0,1a0,2a0 + pb		160,260,060 + pb		220,020,120 + pb
		0e0,1e0,2e0 + pa		1a0,2a0,0a0 + pa		260,060,160 + pa		020,120,220 + pa
		1e0,2e0,0e0 + p9		2a0,0a0,1a0 + p9		060,160,260 + p9		120,220,020 + p9
		2e0,0e0,1e0 + p8		0a0,1a0,2a0 + p8		160,260,060 + p8		220,020,120 + p8
		0e0,1e0,2e0 + p7		1a0,2a0,0a0 + p7		260,060,160 + p7		020,120,220 + p7
		1e0,2e0,0e0 + p6		2a0,0a0,1a0 + p6		060,160,260 + p6		120,220,020 + p6
		2e0,0e0,1e0 + p5		0a0,1a0,2a0 + p5		160,260,060 + p5		220,020,120 + p5
		0e0,1e0,2e0 + p4		1a0,2a0,0a0 + p4		260,060,160 + p4		020,120,220 + p4
		1e0,2e0,0e0 + p3		2a0,0a0,1a0 + p3		060,160,260 + p3		120,220,020 + p3
		2e0,0e0,1e0 + p2		0a0,1a0,2a0 + p2		160,260,060 + p2		220,020,120 + p2
		0e0,1e0,2e0 + p1		1a0,2a0,0a0 + p1		260,060,160 + p1		020,120,220 + p1
		1e0,2e0,0e0 + p0		2a0,0a0,1a0 + p0		060,160,260 + p0		120,220,020 + p0

		2d0,0d0,1d0 + pf		090,190,290 + pf		150,250,050 + pf		210,010,110 + pf
		0d0,1d0,2d0 + pe		190,290,090 + pe		250,050,150 + pe		010,110,210 + pe
		1d0,2d0,0d0 + pd		290,090,190 + pd		050,150,250 + pd		110,210,010 + pd
		2d0,0d0,1d0 + pc		090,190,290 + pc		150,250,050 + pc		210,010,110 + pc
		0d0,1d0,2d0 + pb		190,290,090 + pb		250,050,150 + pb		010,110,210 + pb
		1d0,2d0,0d0 + pa		290,090,190 + pa		050,150,250 + pa		110,210,010 + pa
		2d0,0d0,1d0 + p9		090,190,290 + p9		150,250,050 + p9		210,010,110 + p9
		0d0,1d0,2d0 + p8		190,290,090 + p8		250,050,150 + p8		010,110,210 + p8
		1d0,2d0,0d0 + p7		290,090,190 + p7		050,150,250 + p7		110,210,010 + p7
		2d0,0d0,1d0 + p6		090,190,290 + p6		150,250,050 + p6		210,010,110 + p6
		0d0,1d0,2d0 + p5		190,290,090 + p5		250,050,150 + p5		010,110,210 + p5
		1d0,2d0,0d0 + p4		290,090,190 + p4		050,150,250 + p4		110,210,010 + p4
		2d0,0d0,1d0 + p3		090,190,290 + p3		150,250,050 + p3		210,010,110 + p3
		0d0,1d0,2d0 + p2		190,290,090 + p2		250,050,150 + p2		010,110,210 + p2
		1d0,2d0,0d0 + p1		290,090,190 + p1		050,150,250 + p1		110,210,010 + p1
		2d0,0d0,1d0 + p0		090,190,290 + p0		150,250,050 + p0		210,010,110 + p0

		0c0,1c0,2c0 + pf		180,280,080 + pf		240,040,140 + pf		000,100,200 + pf
		1c0,2c0,0c0 + pe		280,080,180 + pe		040,140,240 + pe		100,200,000 + pe
		2c0,0c0,1c0 + pd		080,180,280 + pd		140,240,040 + pd		200,000,100 + pd
		0c0,1c0,2c0 + pc		180,280,080 + pc		240,040,140 + pc		000,100,200 + pc
		1c0,2c0,0c0 + pb		280,080,180 + pb		040,140,240 + pb		100,200,000 + pb
		2c0,0c0,1c0 + pa		080,180,280 + pa		140,240,040 + pa		200,000,100 + pa
		0c0,1c0,2c0 + p9		180,280,080 + p9		240,040,140 + p9		000,100,200 + p9
		1c0,2c0,0c0 + p8		280,080,180 + p8		040,140,240 + p8		100,200,000 + p8
		2c0,0c0,1c0 + p7		080,180,280 + p7		140,240,040 + p7		200,000,100 + p7
		0c0,1c0,2c0 + p6		180,280,080 + p6		240,040,140 + p6		000,100,200 + p6
		1c0,2c0,0c0 + p5		280,080,180 + p5		040,140,240 + p5		100,200,000 + p5
		2c0,0c0,1c0 + p4		080,180,280 + p4		140,240,040 + p4		200,000,100 + p4
		0c0,1c0,2c0 + p3		180,280,080 + p3		240,040,140 + p3		000,100,200 + p3
		1c0,2c0,0c0 + p2		280,080,180 + p2		040,140,240 + p2		100,200,000 + p2
		2c0,0c0,1c0 + p1		080,180,280 + p1		140,240,040 + p1		200,000,100 + p1
		0c0,1c0,2c0 + p0		180,280,080 + p0		240,040,140 + p0		000,100,200 + p0	<*** 0-term wrapped ***
		[cont. in col2]			[cont. in col3]			[cont. in col4]

	To handle the wraparound of the 0-term we just need a bit of mod-256 indexing magic.
	*/
	/*...gather the needed data (768 64-bit complex) and do 3 radix-256 transforms,	*/

		// Since o_data are in a local-complex scratch array, need to override any non-unity SIMD-re/im data stride with 1:
		//                                                                                              vvv
		RADIX_256_DIT((a+j1     ),RE_IM_STRIDE,(int *)(i_offsets_lo+0x00),i_idx1,i_offsets_hi1, (double *)(t+0x000),1,o_offsets_lo,o_offsets_hi);	/* Outputs in t[ 00- ff] */
		RADIX_256_DIT((a+j1+p200),RE_IM_STRIDE,(int *)(i_offsets_lo+0x20),i_idx2,i_offsets_hi2, (double *)(t+0x100),1,o_offsets_lo,o_offsets_hi);	/* Outputs in t[100-1ff] */
		RADIX_256_DIT((a+j1+p100),RE_IM_STRIDE,(int *)(i_offsets_lo+0x40),i_idx3,i_offsets_hi3, (double *)(t+0x200),1,o_offsets_lo,o_offsets_hi);	/* Outputs in t[200-2ff] */

	/*...and now do 256 radix-3 transforms: */
		// Loop-based compact-obj-code impl exploits above index pattern to group interior [sandwiched
		// between single leading and 15 trailing] macro calls into sets of 16 with neatly cutoff index groupings:
		l = 1; l1 = l+256; l2 = l+512;	// Skip 0-term, which gets saved for wraparound
		for(k = 0; k < 48; k += 3) {
			k0 = dit_triplets[k]; k1 = dit_triplets[k+1]; k2 = dit_triplets[k+2];
			jt = j1 + pf; jp = j2 + pf;	RADIX_03_DFT(s,c3m1, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im, t00,t01,t02,t03,t04,t05, a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2]); ++l; ++l1; ++l2;
			jt = j1 + pe; jp = j2 + pe;	RADIX_03_DFT(s,c3m1, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im, t00,t01,t02,t03,t04,t05, a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2],a[jt+k0],a[jp+k0]); ++l; ++l1; ++l2;
			jt = j1 + pd; jp = j2 + pd;	RADIX_03_DFT(s,c3m1, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im, t00,t01,t02,t03,t04,t05, a[jt+k2],a[jp+k2],a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1]); ++l; ++l1; ++l2;
			jt = j1 + pc; jp = j2 + pc;	RADIX_03_DFT(s,c3m1, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im, t00,t01,t02,t03,t04,t05, a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2]); ++l; ++l1; ++l2;
			jt = j1 + pb; jp = j2 + pb;	RADIX_03_DFT(s,c3m1, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im, t00,t01,t02,t03,t04,t05, a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2],a[jt+k0],a[jp+k0]); ++l; ++l1; ++l2;
			jt = j1 + pa; jp = j2 + pa;	RADIX_03_DFT(s,c3m1, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im, t00,t01,t02,t03,t04,t05, a[jt+k2],a[jp+k2],a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1]); ++l; ++l1; ++l2;
			jt = j1 + p9; jp = j2 + p9;	RADIX_03_DFT(s,c3m1, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im, t00,t01,t02,t03,t04,t05, a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2]); ++l; ++l1; ++l2;
			jt = j1 + p8; jp = j2 + p8;	RADIX_03_DFT(s,c3m1, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im, t00,t01,t02,t03,t04,t05, a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2],a[jt+k0],a[jp+k0]); ++l; ++l1; ++l2;
			jt = j1 + p7; jp = j2 + p7;	RADIX_03_DFT(s,c3m1, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im, t00,t01,t02,t03,t04,t05, a[jt+k2],a[jp+k2],a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1]); ++l; ++l1; ++l2;
			jt = j1 + p6; jp = j2 + p6;	RADIX_03_DFT(s,c3m1, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im, t00,t01,t02,t03,t04,t05, a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2]); ++l; ++l1; ++l2;
			jt = j1 + p5; jp = j2 + p5;	RADIX_03_DFT(s,c3m1, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im, t00,t01,t02,t03,t04,t05, a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2],a[jt+k0],a[jp+k0]); ++l; ++l1; ++l2;
			jt = j1 + p4; jp = j2 + p4;	RADIX_03_DFT(s,c3m1, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im, t00,t01,t02,t03,t04,t05, a[jt+k2],a[jp+k2],a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1]); ++l; ++l1; ++l2;
			jt = j1 + p3; jp = j2 + p3;	RADIX_03_DFT(s,c3m1, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im, t00,t01,t02,t03,t04,t05, a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2]); ++l; ++l1; ++l2;
			jt = j1 + p2; jp = j2 + p2;	RADIX_03_DFT(s,c3m1, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im, t00,t01,t02,t03,t04,t05, a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2],a[jt+k0],a[jp+k0]); ++l; ++l1; ++l2;
			jt = j1 + p1; jp = j2 + p1;	RADIX_03_DFT(s,c3m1, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im, t00,t01,t02,t03,t04,t05, a[jt+k2],a[jp+k2],a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1]); ++l; l &= 0xff; l1 = l+256; l2 = l+512;	// <*** needed for final-loop-pass wraparound of 0-term
										RADIX_03_DFT(s,c3m1, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im, t00,t01,t02,t03,t04,t05, a[j1+k0],a[j2+k0],a[j1+k1],a[j2+k1],a[j1+k2],a[j2+k2]); ++l; ++l1; ++l2;
		}
	}
}

/******************** Multithreaded function body - NO STATIC VARS BELOW THIS POINT!: ***************************/

#ifdef USE_PTHREAD

	#ifndef COMPILER_TYPE_GCC
		#error pthreaded carry code requires GCC build!
	#endif

	void*
	cy768_process_chunk(void*targ)	// Thread-arg pointer *must* be cast to void and specialized inside the function
	{
		struct cy_thread_data_t* thread_arg = targ;	// Move to top because scalar-mode carry pointers taken directly from it
	#if defined(USE_FMA)
		const double tan = 0.41421356237309504879;
	#endif
		double *addr;
	#ifdef USE_SSE2
		const int pfetch_dist = PFETCH_DIST;
	#endif
		const int stride = (int)RE_IM_STRIDE << 1;	// main-array loop stride = 2*RE_IM_STRIDE
		uint32 p1,p2,p3,p4,p5,p6,p7,p8,p9,pa,pb,pc,pd,pe,pf,
		     p10,p20,p30,p40,p50,p60,p70,p80,p90,pa0,pb0,pc0,pd0,pe0,pf0,
		p100,p110,p120,p130,p140,p150,p160,p170,p180,p190,p1a0,p1b0,p1c0,p1d0,p1e0,p1f0,
		p200,p210,p220,p230,p240,p250,p260,p270,p280,p290,p2a0,p2b0,p2c0,p2d0,p2e0,p2f0;
		int poff[RADIX>>2];
	#ifndef USE_SSE2
		double wt_re,wt_im, wi_re,wi_im;	// Fermat-mod/LOACC weights stuff, used in both scalar and SIMD mode
	#endif

		int dif_offsets_lo[64];	// 4 subsets of 16
		// Bitfields encoding the sequence of the dif_offsets_lo subset0-3 vectors to use for each radix-256 DIF's outputs:
		const uint32 dif_idx1 = 0xe79e79e4,dif_idx2 = 0x79e79e79,dif_idx3 = 0x9e79e79e;
		int dif_offsets_hi1[16],dif_offsets_hi2[16],dif_offsets_hi3[16];

		int dit_offsets_lo[96];	// 6 subsets of 16
		// Bitfields encoding the sequence of the dit_offsets_lo subset0-5 vectors to use for each radix-256 DIT's outputs:
		const uint32 dit_idx1 = 0x55555554,dit_idx2 = 0x44404040,dit_idx3 = 0x54445444;
		int dit_offsets_hi1[16],dit_offsets_hi2[16],dit_offsets_hi3[16];

		int dif_triplets[144], dit_triplets[48];

		int j,j1,k,l,l1,l2,k0,k1,k2;
	#ifndef USE_SSE2
		int j2;
	#endif
	#ifdef USE_SSE2
		// incr = Carry-chain wts-multipliers recurrence length, which must divide
		// RADIX/[n-wayness of carry macro], e.g. RADIX/[16|8|4] = 48|96|192 for avx512,avx,sse, respectively:
		int incr;
		const int incr_long = 16,incr_med = 8,incr_short = 4;
	  // Have no specialized HIACC carry macro in USE_AVX512 and ARMv8 SIMD, so use "goes to 11" in LOACC mode via an incr_hiacc = 2:
	  #if defined(USE_AVX512) || defined(USE_ARM_V8_SIMD)
		const int incr_hiacc = 2;
	  #else
		const int incr_hiacc = 0;
	  #endif
		// Allows cy-macro error data to be used to fiddle incr on the fly to a smaller, safer value if necessary
		if(USE_SHORT_CY_CHAIN == 0)
			incr = incr_long;
		else if(USE_SHORT_CY_CHAIN == 1)
			incr = incr_med;
		else if(USE_SHORT_CY_CHAIN == 2)
			incr = incr_short;
		else
			incr = incr_hiacc;
	#endif

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

	  #ifndef USE_AVX512
		const double crnd = 3.0*0x4000000*0x2000000;
	  #endif
		double *add0,*add1,*add2,*add3;
		int *bjmodn;	// Alloc mem for this along with other 	SIMD stuff
		vec_dbl *tmp,*tm1,*tm2;	// utility ptrs
		int *itmp;			// Pointer into the bjmodn array
	  #if defined(USE_AVX) && !defined(USE_AVX512)
		int *itm2;			// Pointer into the bjmodn array
	  #endif
	  #ifndef USE_AVX
		struct complex *ctmp;	// Hybrid AVX-DFT/SSE2-carry scheme used for Mersenne-mod needs a 2-word-double pointer
	  #endif
		vec_dbl *two,/* *one,*sqrt2, */*isrt2, /* *cc0, *ss0, */ *cc1, /* *ss1, */ *max_err,
		  #ifndef USE_AVX512
			*sse2_rnd,
		  #endif
			*half_arr,
			// ptrs to 16 sets of twiddles shared by the 2nd-half DIF and DIT DFT macros:
			*twid0,//*twid1,*twid2,*twid3,*twid4,*twid5,*twid6,*twid7,*twid8,*twid9,*twida,*twidb,*twidc,*twidd,*twide,*twidf,
			*r000,*r100,*r200,	// Head of RADIX*vec_cmplx-sized local store #1
			*s1p000,	// Head of RADIX*vec_cmplx-sized local store #2
			*cy;	// Need RADIX/2 slots for sse2 carries, RADIX/4 for avx
	  #ifndef USE_AVX512
		double dtmp;
	  #endif
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

		int t_offsets_lo[16], t_offsets_hi[16];
		// Local storage: We must use an array here because scalars have no guarantees about relative address offsets
		// [and even if those are contiguous-as-hoped-for, they may run in reverse]; Make array type (struct complex)
		// to allow us to use the same offset-indexing as in the original radix-32 in-place DFT macros:
		struct complex t[RADIX];
		int *itmp;	// Pointer into the bjmodn array

	#endif

	// int data:
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
		double scale = thread_arg->scale;
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
		p1 = NDIVR;		p10 = NDIVR<<4;		p100 = NDIVR<<8;		p200 = NDIVR<<9;
		p2 = p1 + p1;	p20 = p10 + p10;	p110 = p100 + p10;		p210 = p200 + p10;
		p3 = p2 + p1;	p30 = p20 + p10;	p120 = p110 + p10;		p220 = p210 + p10;
		p4 = p3 + p1;	p40 = p30 + p10;	p130 = p120 + p10;		p230 = p220 + p10;
		p5 = p4 + p1;	p50 = p40 + p10;	p140 = p130 + p10;		p240 = p230 + p10;
		p6 = p5 + p1;	p60 = p50 + p10;	p150 = p140 + p10;		p250 = p240 + p10;
		p7 = p6 + p1;	p70 = p60 + p10;	p160 = p150 + p10;		p260 = p250 + p10;
		p8 = p7 + p1;	p80 = p70 + p10;	p170 = p160 + p10;		p270 = p260 + p10;
		p9 = p8 + p1;	p90 = p80 + p10;	p180 = p170 + p10;		p280 = p270 + p10;
		pa = p9 + p1;	pa0 = p90 + p10;	p190 = p180 + p10;		p290 = p280 + p10;
		pb = pa + p1;	pb0 = pa0 + p10;	p1a0 = p190 + p10;		p2a0 = p290 + p10;
		pc = pb + p1;	pc0 = pb0 + p10;	p1b0 = p1a0 + p10;		p2b0 = p2a0 + p10;
		pd = pc + p1;	pd0 = pc0 + p10;	p1c0 = p1b0 + p10;		p2c0 = p2b0 + p10;
		pe = pd + p1;	pe0 = pd0 + p10;	p1d0 = p1c0 + p10;		p2d0 = p2c0 + p10;
		pf = pe + p1;	pf0 = pe0 + p10;	p1e0 = p1d0 + p10;		p2e0 = p2d0 + p10;
											p1f0 = p1e0 + p10;		p2f0 = p2e0 + p10;
		p1 += ( (p1 >> DAT_BITS) << PAD_BITS );
		p2 += ( (p2 >> DAT_BITS) << PAD_BITS );
		p3 += ( (p3 >> DAT_BITS) << PAD_BITS );
		p4 += ( (p4 >> DAT_BITS) << PAD_BITS );
		p5 += ( (p5 >> DAT_BITS) << PAD_BITS );
		p6 += ( (p6 >> DAT_BITS) << PAD_BITS );
		p7 += ( (p7 >> DAT_BITS) << PAD_BITS );
		p8 += ( (p8 >> DAT_BITS) << PAD_BITS );
		p9 += ( (p9 >> DAT_BITS) << PAD_BITS );
		pa += ( (pa >> DAT_BITS) << PAD_BITS );
		pb += ( (pb >> DAT_BITS) << PAD_BITS );
		pc += ( (pc >> DAT_BITS) << PAD_BITS );
		pd += ( (pd >> DAT_BITS) << PAD_BITS );
		pe += ( (pe >> DAT_BITS) << PAD_BITS );
		pf += ( (pf >> DAT_BITS) << PAD_BITS );
		p10 += ( (p10 >> DAT_BITS) << PAD_BITS );
		p20 += ( (p20 >> DAT_BITS) << PAD_BITS );
		p30 += ( (p30 >> DAT_BITS) << PAD_BITS );
		p40 += ( (p40 >> DAT_BITS) << PAD_BITS );
		p50 += ( (p50 >> DAT_BITS) << PAD_BITS );
		p60 += ( (p60 >> DAT_BITS) << PAD_BITS );
		p70 += ( (p70 >> DAT_BITS) << PAD_BITS );
		p80 += ( (p80 >> DAT_BITS) << PAD_BITS );
		p90 += ( (p90 >> DAT_BITS) << PAD_BITS );
		pa0 += ( (pa0 >> DAT_BITS) << PAD_BITS );
		pb0 += ( (pb0 >> DAT_BITS) << PAD_BITS );
		pc0 += ( (pc0 >> DAT_BITS) << PAD_BITS );
		pd0 += ( (pd0 >> DAT_BITS) << PAD_BITS );
		pe0 += ( (pe0 >> DAT_BITS) << PAD_BITS );
		pf0 += ( (pf0 >> DAT_BITS) << PAD_BITS );
		p100 += ( (p100 >> DAT_BITS) << PAD_BITS );
		p110 += ( (p110 >> DAT_BITS) << PAD_BITS );
		p120 += ( (p120 >> DAT_BITS) << PAD_BITS );
		p130 += ( (p130 >> DAT_BITS) << PAD_BITS );
		p140 += ( (p140 >> DAT_BITS) << PAD_BITS );
		p150 += ( (p150 >> DAT_BITS) << PAD_BITS );
		p160 += ( (p160 >> DAT_BITS) << PAD_BITS );
		p170 += ( (p170 >> DAT_BITS) << PAD_BITS );
		p180 += ( (p180 >> DAT_BITS) << PAD_BITS );
		p190 += ( (p190 >> DAT_BITS) << PAD_BITS );
		p1a0 += ( (p1a0 >> DAT_BITS) << PAD_BITS );
		p1b0 += ( (p1b0 >> DAT_BITS) << PAD_BITS );
		p1c0 += ( (p1c0 >> DAT_BITS) << PAD_BITS );
		p1d0 += ( (p1d0 >> DAT_BITS) << PAD_BITS );
		p1e0 += ( (p1e0 >> DAT_BITS) << PAD_BITS );
		p1f0 += ( (p1f0 >> DAT_BITS) << PAD_BITS );
		p200 += ( (p200 >> DAT_BITS) << PAD_BITS );
		p210 += ( (p210 >> DAT_BITS) << PAD_BITS );
		p220 += ( (p220 >> DAT_BITS) << PAD_BITS );
		p230 += ( (p230 >> DAT_BITS) << PAD_BITS );
		p240 += ( (p240 >> DAT_BITS) << PAD_BITS );
		p250 += ( (p250 >> DAT_BITS) << PAD_BITS );
		p260 += ( (p260 >> DAT_BITS) << PAD_BITS );
		p270 += ( (p270 >> DAT_BITS) << PAD_BITS );
		p280 += ( (p280 >> DAT_BITS) << PAD_BITS );
		p290 += ( (p290 >> DAT_BITS) << PAD_BITS );
		p2a0 += ( (p2a0 >> DAT_BITS) << PAD_BITS );
		p2b0 += ( (p2b0 >> DAT_BITS) << PAD_BITS );
		p2c0 += ( (p2c0 >> DAT_BITS) << PAD_BITS );
		p2d0 += ( (p2d0 >> DAT_BITS) << PAD_BITS );
		p2e0 += ( (p2e0 >> DAT_BITS) << PAD_BITS );
		p2f0 += ( (p2f0 >> DAT_BITS) << PAD_BITS );
	#ifndef USE_SSE2
		p0123[0] = 0; p0123[1] = p1; p0123[2] = p2; p0123[3] = p3;
	#endif
		poff[     0] =   0; poff[     1] =     p4; poff[     2] =     p8; poff[     3] =     pc;
		poff[0x04+0] = p10; poff[0x04+1] = p10+p4; poff[0x04+2] = p10+p8; poff[0x04+3] = p10+pc;
		poff[0x08+0] = p20; poff[0x08+1] = p20+p4; poff[0x08+2] = p20+p8; poff[0x08+3] = p20+pc;
		poff[0x0c+0] = p30; poff[0x0c+1] = p30+p4; poff[0x0c+2] = p30+p8; poff[0x0c+3] = p30+pc;
		poff[0x10+0] = p40; poff[0x10+1] = p40+p4; poff[0x10+2] = p40+p8; poff[0x10+3] = p40+pc;
		poff[0x14+0] = p50; poff[0x14+1] = p50+p4; poff[0x14+2] = p50+p8; poff[0x14+3] = p50+pc;
		poff[0x18+0] = p60; poff[0x18+1] = p60+p4; poff[0x18+2] = p60+p8; poff[0x18+3] = p60+pc;
		poff[0x1c+0] = p70; poff[0x1c+1] = p70+p4; poff[0x1c+2] = p70+p8; poff[0x1c+3] = p70+pc;
		poff[0x20+0] = p80; poff[0x20+1] = p80+p4; poff[0x20+2] = p80+p8; poff[0x20+3] = p80+pc;
		poff[0x24+0] = p90; poff[0x24+1] = p90+p4; poff[0x24+2] = p90+p8; poff[0x24+3] = p90+pc;
		poff[0x28+0] = pa0; poff[0x28+1] = pa0+p4; poff[0x28+2] = pa0+p8; poff[0x28+3] = pa0+pc;
		poff[0x2c+0] = pb0; poff[0x2c+1] = pb0+p4; poff[0x2c+2] = pb0+p8; poff[0x2c+3] = pb0+pc;
		poff[0x30+0] = pc0; poff[0x30+1] = pc0+p4; poff[0x30+2] = pc0+p8; poff[0x30+3] = pc0+pc;
		poff[0x34+0] = pd0; poff[0x34+1] = pd0+p4; poff[0x34+2] = pd0+p8; poff[0x34+3] = pd0+pc;
		poff[0x38+0] = pe0; poff[0x38+1] = pe0+p4; poff[0x38+2] = pe0+p8; poff[0x38+3] = pe0+pc;
		poff[0x3c+0] = pf0; poff[0x3c+1] = pf0+p4; poff[0x3c+2] = pf0+p8; poff[0x3c+3] = pf0+pc;
		for(l = 0; l < 64; l++) {
			poff[ 64+l] = poff[l] + p100;
			poff[128+l] = poff[l] + p200;
		}

	// Cf. radix768_dit_pass1() for details on the indexing scheme used here:
	// Set dit_offsets_lo for 1st set of radix-256 DIT inputs:
		dit_offsets_lo[0x00] =  0;	dit_offsets_lo[0x10] = pf;
		dit_offsets_lo[0x01] = p1;	dit_offsets_lo[0x11] = pe;
		dit_offsets_lo[0x02] = p3;	dit_offsets_lo[0x12] = pd;
		dit_offsets_lo[0x03] = p2;	dit_offsets_lo[0x13] = pc;
		dit_offsets_lo[0x04] = p7;	dit_offsets_lo[0x14] = pb;
		dit_offsets_lo[0x05] = p6;	dit_offsets_lo[0x15] = pa;
		dit_offsets_lo[0x06] = p5;	dit_offsets_lo[0x16] = p9;
		dit_offsets_lo[0x07] = p4;	dit_offsets_lo[0x17] = p8;
		dit_offsets_lo[0x08] = pf;	dit_offsets_lo[0x18] = p7;
		dit_offsets_lo[0x09] = pe;	dit_offsets_lo[0x19] = p6;
		dit_offsets_lo[0x0a] = pd;	dit_offsets_lo[0x1a] = p5;
		dit_offsets_lo[0x0b] = pc;	dit_offsets_lo[0x1b] = p4;
		dit_offsets_lo[0x0c] = pb;	dit_offsets_lo[0x1c] = p3;
		dit_offsets_lo[0x0d] = pa;	dit_offsets_lo[0x1d] = p2;
		dit_offsets_lo[0x0e] = p9;	dit_offsets_lo[0x1e] = p1;
		dit_offsets_lo[0x0f] = p8;	dit_offsets_lo[0x1f] =  0;

	// Set dit_offsets for 2nd set of radix-256 DIT inputs:
		dit_offsets_lo[0x20] = p5;	dit_offsets_lo[0x30] = p9;
		dit_offsets_lo[0x21] = p4;	dit_offsets_lo[0x31] = p8;
		dit_offsets_lo[0x22] = p6;	dit_offsets_lo[0x32] = pa;
		dit_offsets_lo[0x23] = p7;	dit_offsets_lo[0x33] = pb;
		dit_offsets_lo[0x24] = p1;	dit_offsets_lo[0x34] = pe;
		dit_offsets_lo[0x25] =  0;	dit_offsets_lo[0x35] = pf;
		dit_offsets_lo[0x26] = p2;	dit_offsets_lo[0x36] = pc;
		dit_offsets_lo[0x27] = p3;	dit_offsets_lo[0x37] = pd;
		dit_offsets_lo[0x28] = p9;	dit_offsets_lo[0x38] = p1;
		dit_offsets_lo[0x29] = p8;	dit_offsets_lo[0x39] =  0;
		dit_offsets_lo[0x2a] = pa;	dit_offsets_lo[0x3a] = p2;
		dit_offsets_lo[0x2b] = pb;	dit_offsets_lo[0x3b] = p3;
		dit_offsets_lo[0x2c] = pe;	dit_offsets_lo[0x3c] = p6;
		dit_offsets_lo[0x2d] = pf;	dit_offsets_lo[0x3d] = p7;
		dit_offsets_lo[0x2e] = pc;	dit_offsets_lo[0x3e] = p4;
		dit_offsets_lo[0x2f] = pd;	dit_offsets_lo[0x3f] = p5;

	// Set dit_offsets for 3rd set of radix-256 DIT inputs:
		dit_offsets_lo[0x40] = pa;	dit_offsets_lo[0x50] = p2;
		dit_offsets_lo[0x41] = pb;	dit_offsets_lo[0x51] = p3;
		dit_offsets_lo[0x42] = p8;	dit_offsets_lo[0x52] =  0;
		dit_offsets_lo[0x43] = p9;	dit_offsets_lo[0x53] = p1;
		dit_offsets_lo[0x44] = pc;	dit_offsets_lo[0x54] = p4;
		dit_offsets_lo[0x45] = pd;	dit_offsets_lo[0x55] = p5;
		dit_offsets_lo[0x46] = pf;	dit_offsets_lo[0x56] = p7;
		dit_offsets_lo[0x47] = pe;	dit_offsets_lo[0x57] = p6;
		dit_offsets_lo[0x48] = p2;	dit_offsets_lo[0x58] = pc;
		dit_offsets_lo[0x49] = p3;	dit_offsets_lo[0x59] = pd;
		dit_offsets_lo[0x4a] =  0;	dit_offsets_lo[0x5a] = pf;
		dit_offsets_lo[0x4b] = p1;	dit_offsets_lo[0x5b] = pe;
		dit_offsets_lo[0x4c] = p4;	dit_offsets_lo[0x5c] = p8;
		dit_offsets_lo[0x4d] = p5;	dit_offsets_lo[0x5d] = p9;
		dit_offsets_lo[0x4e] = p7;	dit_offsets_lo[0x5e] = pb;
		dit_offsets_lo[0x4f] = p6;	dit_offsets_lo[0x5f] = pa;

	// ...and one distinct high-part vector for each radix-256 DIT:
		dit_offsets_hi1[0x0] =   0;	dit_offsets_hi2[0x0] = p50;	dit_offsets_hi3[0x0] = pa0;
		dit_offsets_hi1[0x1] = p10;	dit_offsets_hi2[0x1] = p40;	dit_offsets_hi3[0x1] = pb0;
		dit_offsets_hi1[0x2] = p30;	dit_offsets_hi2[0x2] = p60;	dit_offsets_hi3[0x2] = p80;
		dit_offsets_hi1[0x3] = p20;	dit_offsets_hi2[0x3] = p70;	dit_offsets_hi3[0x3] = p90;
		dit_offsets_hi1[0x4] = p70;	dit_offsets_hi2[0x4] = p10;	dit_offsets_hi3[0x4] = pc0;
		dit_offsets_hi1[0x5] = p60;	dit_offsets_hi2[0x5] =   0;	dit_offsets_hi3[0x5] = pd0;
		dit_offsets_hi1[0x6] = p50;	dit_offsets_hi2[0x6] = p20;	dit_offsets_hi3[0x6] = pf0;
		dit_offsets_hi1[0x7] = p40;	dit_offsets_hi2[0x7] = p30;	dit_offsets_hi3[0x7] = pe0;
		dit_offsets_hi1[0x8] = pf0;	dit_offsets_hi2[0x8] = p90;	dit_offsets_hi3[0x8] = p20;
		dit_offsets_hi1[0x9] = pe0;	dit_offsets_hi2[0x9] = p80;	dit_offsets_hi3[0x9] = p30;
		dit_offsets_hi1[0xa] = pd0;	dit_offsets_hi2[0xa] = pa0;	dit_offsets_hi3[0xa] =   0;
		dit_offsets_hi1[0xb] = pc0;	dit_offsets_hi2[0xb] = pb0;	dit_offsets_hi3[0xb] = p10;
		dit_offsets_hi1[0xc] = pb0;	dit_offsets_hi2[0xc] = pe0;	dit_offsets_hi3[0xc] = p40;
		dit_offsets_hi1[0xd] = pa0;	dit_offsets_hi2[0xd] = pf0;	dit_offsets_hi3[0xd] = p50;
		dit_offsets_hi1[0xe] = p90;	dit_offsets_hi2[0xe] = pc0;	dit_offsets_hi3[0xe] = p70;
		dit_offsets_hi1[0xf] = p80;	dit_offsets_hi2[0xf] = pd0;	dit_offsets_hi3[0xf] = p60;

	  #ifndef USE_SSE2
		// Idx offsets w.r.to r-array are const and shared by both sets of radix-256 transforms:
		// low parts in first 16 slots, high parts in next 16:
		t_offsets_lo[0x0] = 0x00<<1;	t_offsets_hi[0x0] = 0x00<<1;
		t_offsets_lo[0x1] = 0x01<<1;	t_offsets_hi[0x1] = 0x10<<1;
		t_offsets_lo[0x2] = 0x02<<1;	t_offsets_hi[0x2] = 0x20<<1;
		t_offsets_lo[0x3] = 0x03<<1;	t_offsets_hi[0x3] = 0x30<<1;
		t_offsets_lo[0x4] = 0x04<<1;	t_offsets_hi[0x4] = 0x40<<1;
		t_offsets_lo[0x5] = 0x05<<1;	t_offsets_hi[0x5] = 0x50<<1;
		t_offsets_lo[0x6] = 0x06<<1;	t_offsets_hi[0x6] = 0x60<<1;
		t_offsets_lo[0x7] = 0x07<<1;	t_offsets_hi[0x7] = 0x70<<1;
		t_offsets_lo[0x8] = 0x08<<1;	t_offsets_hi[0x8] = 0x80<<1;
		t_offsets_lo[0x9] = 0x09<<1;	t_offsets_hi[0x9] = 0x90<<1;
		t_offsets_lo[0xa] = 0x0a<<1;	t_offsets_hi[0xa] = 0xa0<<1;
		t_offsets_lo[0xb] = 0x0b<<1;	t_offsets_hi[0xb] = 0xb0<<1;
		t_offsets_lo[0xc] = 0x0c<<1;	t_offsets_hi[0xc] = 0xc0<<1;
		t_offsets_lo[0xd] = 0x0d<<1;	t_offsets_hi[0xd] = 0xd0<<1;
		t_offsets_lo[0xe] = 0x0e<<1;	t_offsets_hi[0xe] = 0xe0<<1;
		t_offsets_lo[0xf] = 0x0f<<1;	t_offsets_hi[0xf] = 0xf0<<1;
	  #endif

	// Cf. radix768_dif_pass1() for details on the indexing scheme used here:
	// Set dit_offsets_lo for the 4 subvectors shared by the 3 sets of radix-256 DIT inputs:
		dif_offsets_lo[0x00] =  0;	dif_offsets_lo[0x10] = p5;	dif_offsets_lo[0x20] = pa;	dif_offsets_lo[0x30] = pf;
		dif_offsets_lo[0x01] = p1;	dif_offsets_lo[0x11] = p4;	dif_offsets_lo[0x21] = pb;	dif_offsets_lo[0x31] = pe;
		dif_offsets_lo[0x02] = p2;	dif_offsets_lo[0x12] = p7;	dif_offsets_lo[0x22] = p9;	dif_offsets_lo[0x32] = pc;
		dif_offsets_lo[0x03] = p3;	dif_offsets_lo[0x13] = p6;	dif_offsets_lo[0x23] = p8;	dif_offsets_lo[0x33] = pd;
		dif_offsets_lo[0x04] = p5;	dif_offsets_lo[0x14] = p2;	dif_offsets_lo[0x24] = pf;	dif_offsets_lo[0x34] = p9;
		dif_offsets_lo[0x05] = p4;	dif_offsets_lo[0x15] = p3;	dif_offsets_lo[0x25] = pe;	dif_offsets_lo[0x35] = p8;
		dif_offsets_lo[0x06] = p7;	dif_offsets_lo[0x16] = p1;	dif_offsets_lo[0x26] = pc;	dif_offsets_lo[0x36] = pb;
		dif_offsets_lo[0x07] = p6;	dif_offsets_lo[0x17] =  0;	dif_offsets_lo[0x27] = pd;	dif_offsets_lo[0x37] = pa;
		dif_offsets_lo[0x08] = pa;	dif_offsets_lo[0x18] = pf;	dif_offsets_lo[0x28] = p5;	dif_offsets_lo[0x38] = p2;
		dif_offsets_lo[0x09] = pb;	dif_offsets_lo[0x19] = pe;	dif_offsets_lo[0x29] = p4;	dif_offsets_lo[0x39] = p3;
		dif_offsets_lo[0x0a] = p9;	dif_offsets_lo[0x1a] = pc;	dif_offsets_lo[0x2a] = p7;	dif_offsets_lo[0x3a] = p1;
		dif_offsets_lo[0x0b] = p8;	dif_offsets_lo[0x1b] = pd;	dif_offsets_lo[0x2b] = p6;	dif_offsets_lo[0x3b] =  0;
		dif_offsets_lo[0x0c] = pf;	dif_offsets_lo[0x1c] = p9;	dif_offsets_lo[0x2c] = p2;	dif_offsets_lo[0x3c] = p7;
		dif_offsets_lo[0x0d] = pe;	dif_offsets_lo[0x1d] = p8;	dif_offsets_lo[0x2d] = p3;	dif_offsets_lo[0x3d] = p6;
		dif_offsets_lo[0x0e] = pc;	dif_offsets_lo[0x1e] = pb;	dif_offsets_lo[0x2e] = p1;	dif_offsets_lo[0x3e] = p4;
		dif_offsets_lo[0x0f] = pd;	dif_offsets_lo[0x1f] = pa;	dif_offsets_lo[0x2f] =  0;	dif_offsets_lo[0x3f] = p5;

	// ...and one distinct high-part vector for each radix-256 DFT:
		dif_offsets_hi1[0x0] =   0;	dif_offsets_hi2[0x0] = p50;	dif_offsets_hi3[0x0] = pa0;
		dif_offsets_hi1[0x1] = p10;	dif_offsets_hi2[0x1] = p40;	dif_offsets_hi3[0x1] = pb0;
		dif_offsets_hi1[0x2] = p20;	dif_offsets_hi2[0x2] = p70;	dif_offsets_hi3[0x2] = p90;
		dif_offsets_hi1[0x3] = p30;	dif_offsets_hi2[0x3] = p60;	dif_offsets_hi3[0x3] = p80;
		dif_offsets_hi1[0x4] = p50;	dif_offsets_hi2[0x4] = p20;	dif_offsets_hi3[0x4] = pf0;
		dif_offsets_hi1[0x5] = p40;	dif_offsets_hi2[0x5] = p30;	dif_offsets_hi3[0x5] = pe0;
		dif_offsets_hi1[0x6] = p70;	dif_offsets_hi2[0x6] = p10;	dif_offsets_hi3[0x6] = pc0;
		dif_offsets_hi1[0x7] = p60;	dif_offsets_hi2[0x7] =   0;	dif_offsets_hi3[0x7] = pd0;
		dif_offsets_hi1[0x8] = pa0;	dif_offsets_hi2[0x8] = pf0;	dif_offsets_hi3[0x8] = p50;
		dif_offsets_hi1[0x9] = pb0;	dif_offsets_hi2[0x9] = pe0;	dif_offsets_hi3[0x9] = p40;
		dif_offsets_hi1[0xa] = p90;	dif_offsets_hi2[0xa] = pc0;	dif_offsets_hi3[0xa] = p70;
		dif_offsets_hi1[0xb] = p80;	dif_offsets_hi2[0xb] = pd0;	dif_offsets_hi3[0xb] = p60;
		dif_offsets_hi1[0xc] = pf0;	dif_offsets_hi2[0xc] = p90;	dif_offsets_hi3[0xc] = p20;
		dif_offsets_hi1[0xd] = pe0;	dif_offsets_hi2[0xd] = p80;	dif_offsets_hi3[0xd] = p30;
		dif_offsets_hi1[0xe] = pc0;	dif_offsets_hi2[0xe] = pb0;	dif_offsets_hi3[0xe] = p10;
		dif_offsets_hi1[0xf] = pd0;	dif_offsets_hi2[0xf] = pa0;	dif_offsets_hi3[0xf] =   0;

	// DIF Index-high-bits triplets needed for compact-obj-code scheme:
	#ifdef USE_SSE2
		k = 0;	// In SIMD mode these are 0xtr-offsets w.r.to a local store:
		dif_triplets[k] = 0x2f0; dif_triplets[k+1] = 0x1f0; dif_triplets[k+2] = 0x0f0; k += 3;
		dif_triplets[k] = 0x2e0; dif_triplets[k+1] = 0x1e0; dif_triplets[k+2] = 0x0e0; k += 3;
		dif_triplets[k] = 0x2d0; dif_triplets[k+1] = 0x1d0; dif_triplets[k+2] = 0x0d0; k += 3;
		dif_triplets[k] = 0x2c0; dif_triplets[k+1] = 0x1c0; dif_triplets[k+2] = 0x0c0; k += 3;
		dif_triplets[k] = 0x2b0; dif_triplets[k+1] = 0x1b0; dif_triplets[k+2] = 0x0b0; k += 3;
		dif_triplets[k] = 0x2a0; dif_triplets[k+1] = 0x1a0; dif_triplets[k+2] = 0x0a0; k += 3;
		dif_triplets[k] = 0x290; dif_triplets[k+1] = 0x190; dif_triplets[k+2] = 0x090; k += 3;
		dif_triplets[k] = 0x280; dif_triplets[k+1] = 0x180; dif_triplets[k+2] = 0x080; k += 3;
		dif_triplets[k] = 0x270; dif_triplets[k+1] = 0x170; dif_triplets[k+2] = 0x070; k += 3;
		dif_triplets[k] = 0x260; dif_triplets[k+1] = 0x160; dif_triplets[k+2] = 0x060; k += 3;
		dif_triplets[k] = 0x250; dif_triplets[k+1] = 0x150; dif_triplets[k+2] = 0x050; k += 3;
		dif_triplets[k] = 0x240; dif_triplets[k+1] = 0x140; dif_triplets[k+2] = 0x040; k += 3;
		dif_triplets[k] = 0x230; dif_triplets[k+1] = 0x130; dif_triplets[k+2] = 0x030; k += 3;
		dif_triplets[k] = 0x220; dif_triplets[k+1] = 0x120; dif_triplets[k+2] = 0x020; k += 3;
		dif_triplets[k] = 0x210; dif_triplets[k+1] = 0x110; dif_triplets[k+2] = 0x010; k += 3;
		dif_triplets[k] = 0x200; dif_triplets[k+1] = 0x100; dif_triplets[k+2] = 0x000; k += 3;
		dif_triplets[k] = 0x1f0; dif_triplets[k+1] = 0x0f0; dif_triplets[k+2] = 0x2f0; k += 3;
		dif_triplets[k] = 0x1e0; dif_triplets[k+1] = 0x0e0; dif_triplets[k+2] = 0x2e0; k += 3;
		dif_triplets[k] = 0x1d0; dif_triplets[k+1] = 0x0d0; dif_triplets[k+2] = 0x2d0; k += 3;
		dif_triplets[k] = 0x1c0; dif_triplets[k+1] = 0x0c0; dif_triplets[k+2] = 0x2c0; k += 3;
		dif_triplets[k] = 0x1b0; dif_triplets[k+1] = 0x0b0; dif_triplets[k+2] = 0x2b0; k += 3;
		dif_triplets[k] = 0x1a0; dif_triplets[k+1] = 0x0a0; dif_triplets[k+2] = 0x2a0; k += 3;
		dif_triplets[k] = 0x190; dif_triplets[k+1] = 0x090; dif_triplets[k+2] = 0x290; k += 3;
		dif_triplets[k] = 0x180; dif_triplets[k+1] = 0x080; dif_triplets[k+2] = 0x280; k += 3;
		dif_triplets[k] = 0x170; dif_triplets[k+1] = 0x070; dif_triplets[k+2] = 0x270; k += 3;
		dif_triplets[k] = 0x160; dif_triplets[k+1] = 0x060; dif_triplets[k+2] = 0x260; k += 3;
		dif_triplets[k] = 0x150; dif_triplets[k+1] = 0x050; dif_triplets[k+2] = 0x250; k += 3;
		dif_triplets[k] = 0x140; dif_triplets[k+1] = 0x040; dif_triplets[k+2] = 0x240; k += 3;
		dif_triplets[k] = 0x130; dif_triplets[k+1] = 0x030; dif_triplets[k+2] = 0x230; k += 3;
		dif_triplets[k] = 0x120; dif_triplets[k+1] = 0x020; dif_triplets[k+2] = 0x220; k += 3;
		dif_triplets[k] = 0x110; dif_triplets[k+1] = 0x010; dif_triplets[k+2] = 0x210; k += 3;
		dif_triplets[k] = 0x100; dif_triplets[k+1] = 0x000; dif_triplets[k+2] = 0x200; k += 3;
		dif_triplets[k] = 0x0f0; dif_triplets[k+1] = 0x2f0; dif_triplets[k+2] = 0x1f0; k += 3;
		dif_triplets[k] = 0x0e0; dif_triplets[k+1] = 0x2e0; dif_triplets[k+2] = 0x1e0; k += 3;
		dif_triplets[k] = 0x0d0; dif_triplets[k+1] = 0x2d0; dif_triplets[k+2] = 0x1d0; k += 3;
		dif_triplets[k] = 0x0c0; dif_triplets[k+1] = 0x2c0; dif_triplets[k+2] = 0x1c0; k += 3;
		dif_triplets[k] = 0x0b0; dif_triplets[k+1] = 0x2b0; dif_triplets[k+2] = 0x1b0; k += 3;
		dif_triplets[k] = 0x0a0; dif_triplets[k+1] = 0x2a0; dif_triplets[k+2] = 0x1a0; k += 3;
		dif_triplets[k] = 0x090; dif_triplets[k+1] = 0x290; dif_triplets[k+2] = 0x190; k += 3;
		dif_triplets[k] = 0x080; dif_triplets[k+1] = 0x280; dif_triplets[k+2] = 0x180; k += 3;
		dif_triplets[k] = 0x070; dif_triplets[k+1] = 0x270; dif_triplets[k+2] = 0x170; k += 3;
		dif_triplets[k] = 0x060; dif_triplets[k+1] = 0x260; dif_triplets[k+2] = 0x160; k += 3;
		dif_triplets[k] = 0x050; dif_triplets[k+1] = 0x250; dif_triplets[k+2] = 0x150; k += 3;
		dif_triplets[k] = 0x040; dif_triplets[k+1] = 0x240; dif_triplets[k+2] = 0x140; k += 3;
		dif_triplets[k] = 0x030; dif_triplets[k+1] = 0x230; dif_triplets[k+2] = 0x130; k += 3;
		dif_triplets[k] = 0x020; dif_triplets[k+1] = 0x220; dif_triplets[k+2] = 0x120; k += 3;
		dif_triplets[k] = 0x010; dif_triplets[k+1] = 0x210; dif_triplets[k+2] = 0x110; k += 3;
		dif_triplets[k] = 0x000; dif_triplets[k+1] = 0x200; dif_triplets[k+2] = 0x100;
		// IN SIMD mode need to double all the above to turn from vec_dbl to vec_cmplx ptr offsets:
		for(l = 0; l < 144; l++) {
			dif_triplets[l] <<= 1;
		}
	#else
		k = 0;
		dif_triplets[k] = p2f0; dif_triplets[k+1] = p1f0; dif_triplets[k+2] = pf0; k += 3;
		dif_triplets[k] = p2e0; dif_triplets[k+1] = p1e0; dif_triplets[k+2] = pe0; k += 3;
		dif_triplets[k] = p2d0; dif_triplets[k+1] = p1d0; dif_triplets[k+2] = pd0; k += 3;
		dif_triplets[k] = p2c0; dif_triplets[k+1] = p1c0; dif_triplets[k+2] = pc0; k += 3;
		dif_triplets[k] = p2b0; dif_triplets[k+1] = p1b0; dif_triplets[k+2] = pb0; k += 3;
		dif_triplets[k] = p2a0; dif_triplets[k+1] = p1a0; dif_triplets[k+2] = pa0; k += 3;
		dif_triplets[k] = p290; dif_triplets[k+1] = p190; dif_triplets[k+2] = p90; k += 3;
		dif_triplets[k] = p280; dif_triplets[k+1] = p180; dif_triplets[k+2] = p80; k += 3;
		dif_triplets[k] = p270; dif_triplets[k+1] = p170; dif_triplets[k+2] = p70; k += 3;
		dif_triplets[k] = p260; dif_triplets[k+1] = p160; dif_triplets[k+2] = p60; k += 3;
		dif_triplets[k] = p250; dif_triplets[k+1] = p150; dif_triplets[k+2] = p50; k += 3;
		dif_triplets[k] = p240; dif_triplets[k+1] = p140; dif_triplets[k+2] = p40; k += 3;
		dif_triplets[k] = p230; dif_triplets[k+1] = p130; dif_triplets[k+2] = p30; k += 3;
		dif_triplets[k] = p220; dif_triplets[k+1] = p120; dif_triplets[k+2] = p20; k += 3;
		dif_triplets[k] = p210; dif_triplets[k+1] = p110; dif_triplets[k+2] = p10; k += 3;
		dif_triplets[k] = p200; dif_triplets[k+1] = p100; dif_triplets[k+2] =   0; k += 3;
		dif_triplets[k] = p1f0; dif_triplets[k+1] = pf0; dif_triplets[k+2] = p2f0; k += 3;
		dif_triplets[k] = p1e0; dif_triplets[k+1] = pe0; dif_triplets[k+2] = p2e0; k += 3;
		dif_triplets[k] = p1d0; dif_triplets[k+1] = pd0; dif_triplets[k+2] = p2d0; k += 3;
		dif_triplets[k] = p1c0; dif_triplets[k+1] = pc0; dif_triplets[k+2] = p2c0; k += 3;
		dif_triplets[k] = p1b0; dif_triplets[k+1] = pb0; dif_triplets[k+2] = p2b0; k += 3;
		dif_triplets[k] = p1a0; dif_triplets[k+1] = pa0; dif_triplets[k+2] = p2a0; k += 3;
		dif_triplets[k] = p190; dif_triplets[k+1] = p90; dif_triplets[k+2] = p290; k += 3;
		dif_triplets[k] = p180; dif_triplets[k+1] = p80; dif_triplets[k+2] = p280; k += 3;
		dif_triplets[k] = p170; dif_triplets[k+1] = p70; dif_triplets[k+2] = p270; k += 3;
		dif_triplets[k] = p160; dif_triplets[k+1] = p60; dif_triplets[k+2] = p260; k += 3;
		dif_triplets[k] = p150; dif_triplets[k+1] = p50; dif_triplets[k+2] = p250; k += 3;
		dif_triplets[k] = p140; dif_triplets[k+1] = p40; dif_triplets[k+2] = p240; k += 3;
		dif_triplets[k] = p130; dif_triplets[k+1] = p30; dif_triplets[k+2] = p230; k += 3;
		dif_triplets[k] = p120; dif_triplets[k+1] = p20; dif_triplets[k+2] = p220; k += 3;
		dif_triplets[k] = p110; dif_triplets[k+1] = p10; dif_triplets[k+2] = p210; k += 3;
		dif_triplets[k] = p100; dif_triplets[k+1] =   0; dif_triplets[k+2] = p200; k += 3;
		dif_triplets[k] = pf0; dif_triplets[k+1] = p2f0; dif_triplets[k+2] = p1f0; k += 3;
		dif_triplets[k] = pe0; dif_triplets[k+1] = p2e0; dif_triplets[k+2] = p1e0; k += 3;
		dif_triplets[k] = pd0; dif_triplets[k+1] = p2d0; dif_triplets[k+2] = p1d0; k += 3;
		dif_triplets[k] = pc0; dif_triplets[k+1] = p2c0; dif_triplets[k+2] = p1c0; k += 3;
		dif_triplets[k] = pb0; dif_triplets[k+1] = p2b0; dif_triplets[k+2] = p1b0; k += 3;
		dif_triplets[k] = pa0; dif_triplets[k+1] = p2a0; dif_triplets[k+2] = p1a0; k += 3;
		dif_triplets[k] = p90; dif_triplets[k+1] = p290; dif_triplets[k+2] = p190; k += 3;
		dif_triplets[k] = p80; dif_triplets[k+1] = p280; dif_triplets[k+2] = p180; k += 3;
		dif_triplets[k] = p70; dif_triplets[k+1] = p270; dif_triplets[k+2] = p170; k += 3;
		dif_triplets[k] = p60; dif_triplets[k+1] = p260; dif_triplets[k+2] = p160; k += 3;
		dif_triplets[k] = p50; dif_triplets[k+1] = p250; dif_triplets[k+2] = p150; k += 3;
		dif_triplets[k] = p40; dif_triplets[k+1] = p240; dif_triplets[k+2] = p140; k += 3;
		dif_triplets[k] = p30; dif_triplets[k+1] = p230; dif_triplets[k+2] = p130; k += 3;
		dif_triplets[k] = p20; dif_triplets[k+1] = p220; dif_triplets[k+2] = p120; k += 3;
		dif_triplets[k] = p10; dif_triplets[k+1] = p210; dif_triplets[k+2] = p110; k += 3;
		dif_triplets[k] =   0; dif_triplets[k+1] = p200; dif_triplets[k+2] = p100;
	#endif
	// DIT Index-high-bits triplets needed for compact-obj-code scheme:
	#ifdef USE_SSE2
		k = 0;	// In SIMD mode these are ptr-offsets w.r.to a local store:
		dit_triplets[k] = 0x0f0; dit_triplets[k+1] = 0x1f0; dit_triplets[k+2] = 0x2f0; k += 3;
		dit_triplets[k] = 0x1e0; dit_triplets[k+1] = 0x2e0; dit_triplets[k+2] = 0x0e0; k += 3;
		dit_triplets[k] = 0x2d0; dit_triplets[k+1] = 0x0d0; dit_triplets[k+2] = 0x1d0; k += 3;
		dit_triplets[k] = 0x0c0; dit_triplets[k+1] = 0x1c0; dit_triplets[k+2] = 0x2c0; k += 3;
		dit_triplets[k] = 0x1b0; dit_triplets[k+1] = 0x2b0; dit_triplets[k+2] = 0x0b0; k += 3;
		dit_triplets[k] = 0x2a0; dit_triplets[k+1] = 0x0a0; dit_triplets[k+2] = 0x1a0; k += 3;
		dit_triplets[k] = 0x090; dit_triplets[k+1] = 0x190; dit_triplets[k+2] = 0x290; k += 3;
		dit_triplets[k] = 0x180; dit_triplets[k+1] = 0x280; dit_triplets[k+2] = 0x080; k += 3;
		dit_triplets[k] = 0x270; dit_triplets[k+1] = 0x070; dit_triplets[k+2] = 0x170; k += 3;
		dit_triplets[k] = 0x060; dit_triplets[k+1] = 0x160; dit_triplets[k+2] = 0x260; k += 3;
		dit_triplets[k] = 0x150; dit_triplets[k+1] = 0x250; dit_triplets[k+2] = 0x050; k += 3;
		dit_triplets[k] = 0x240; dit_triplets[k+1] = 0x040; dit_triplets[k+2] = 0x140; k += 3;
		dit_triplets[k] = 0x030; dit_triplets[k+1] = 0x130; dit_triplets[k+2] = 0x230; k += 3;
		dit_triplets[k] = 0x120; dit_triplets[k+1] = 0x220; dit_triplets[k+2] = 0x020; k += 3;
		dit_triplets[k] = 0x210; dit_triplets[k+1] = 0x010; dit_triplets[k+2] = 0x110; k += 3;
		dit_triplets[k] = 0x000; dit_triplets[k+1] = 0x100; dit_triplets[k+2] = 0x200;
		// IN SIMD mode need to double all the above to turn from vec_dbl to vec_cmplx ptr offsets:
		for(l = 0; l < 48; l++) {
			dit_triplets[l] <<= 1;
		}
	#else
		k = 0;
		dit_triplets[k] =  pf0; dit_triplets[k+1] = p1f0; dit_triplets[k+2] = p2f0; k += 3;
		dit_triplets[k] = p1e0; dit_triplets[k+1] = p2e0; dit_triplets[k+2] =  pe0; k += 3;
		dit_triplets[k] = p2d0; dit_triplets[k+1] =  pd0; dit_triplets[k+2] = p1d0; k += 3;
		dit_triplets[k] =  pc0; dit_triplets[k+1] = p1c0; dit_triplets[k+2] = p2c0; k += 3;
		dit_triplets[k] = p1b0; dit_triplets[k+1] = p2b0; dit_triplets[k+2] =  pb0; k += 3;
		dit_triplets[k] = p2a0; dit_triplets[k+1] =  pa0; dit_triplets[k+2] = p1a0; k += 3;
		dit_triplets[k] =  p90; dit_triplets[k+1] = p190; dit_triplets[k+2] = p290; k += 3;
		dit_triplets[k] = p180; dit_triplets[k+1] = p280; dit_triplets[k+2] =  p80; k += 3;
		dit_triplets[k] = p270; dit_triplets[k+1] =  p70; dit_triplets[k+2] = p170; k += 3;
		dit_triplets[k] =  p60; dit_triplets[k+1] = p160; dit_triplets[k+2] = p260; k += 3;
		dit_triplets[k] = p150; dit_triplets[k+1] = p250; dit_triplets[k+2] =  p50; k += 3;
		dit_triplets[k] = p240; dit_triplets[k+1] =  p40; dit_triplets[k+2] = p140; k += 3;
		dit_triplets[k] =  p30; dit_triplets[k+1] = p130; dit_triplets[k+2] = p230; k += 3;
		dit_triplets[k] = p120; dit_triplets[k+1] = p220; dit_triplets[k+2] =  p20; k += 3;
		dit_triplets[k] = p210; dit_triplets[k+1] =  p10; dit_triplets[k+2] = p110; k += 3;
		dit_triplets[k] =    0; dit_triplets[k+1] = p100; dit_triplets[k+2] = p200;
	#endif

	#ifdef USE_SSE2
		tmp = r000 = thread_arg->r000;
		r100 = tmp + 0x200;
		r200 = tmp + 0x400;
		tmp += 0x600;	s1p000 = tmp;
		tmp += 0x600;	// r000 += 0xc00
		two     = tmp + 0;	// AVX+ versions of radix-8,16,32 twiddleless-DFT macros need consts [2,1,sqrt2,isrt2] quartet laid out thusly
		//one     = tmp + 1;
		//sqrt2   = tmp + 2;
		isrt2   = tmp + 3;
		//cc0		= tmp + 4;
		//ss0		= tmp + 5;
		cc1		= tmp + 6;
		//ss1		= tmp + 7;
		tmp += 0x08;	// sc_ptr += 0xc08
		// ptrs to 15 sets (30 vec_dbl data each) of non-unity twiddles shared by the 2nd-half DIF and DIT DFT macros:
		twid0  = tmp + 0x00;
		//twid1  = tmp + 0x1e;
		//twid2  = tmp + 0x3c;
		//twid3  = tmp + 0x5a;
		//twid4  = tmp + 0x78;
		//twid5  = tmp + 0x96;
		//twid6  = tmp + 0xb4;
		//twid7  = tmp + 0xd2;
		//twid8  = tmp + 0xf0;
		//twid9  = tmp + 0x10e;
		//twida  = tmp + 0x12c;
		//twidb  = tmp + 0x14a;
		//twidc  = tmp + 0x168;
		//twidd  = tmp + 0x186;
		//twide  = tmp + 0x1a4;
		//twidf  = tmp + 0x1c2;
		tmp += 0x1e0;	// += 15*30 => sc_ptr += 0xde8
	  #ifdef USE_AVX512
		cy = tmp;		tmp += 0x60;	// RADIX/8 vec_dbl slots for carry sub-array
		max_err = tmp + 0x00;
		//sse2_rnd= tmp + 0x01;
		half_arr= tmp + 0x02;
	  #elif defined(USE_AVX)
		cy = tmp;		tmp += 0xc0;	// RADIX/4 vec_dbl slots for carry sub-array
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// += 0xc0 + 2 => sc_ptr += 0xeaa
		// This is where the value of half_arr_offset comes from
		half_arr= tmp + 0x02;	/* This table needs 68 vec_dbl for Mersenne-mod, and 3.5*RADIX[avx] | RADIX[sse2] for Fermat-mod */
	  #else
		cy = tmp;		tmp += 0x180;	// RADIX/2 vec_dbl slots for carry sub-array
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// += 0x180 + 2 => sc_ptr += 0xf6a
		// This is where the value of half_arr_offset comes from
		half_arr= tmp + 0x02;	/* This table needs 20 x 16 bytes for Mersenne-mod, 2 for Fermat-mod */
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

		sign_mask = (uint64*)(r000 + radix768_creals_in_local_store);
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
		#include "radix768_main_carry_loop.h"

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
